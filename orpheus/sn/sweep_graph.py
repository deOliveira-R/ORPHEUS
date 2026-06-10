r"""Per-octant causal cell DAG for the 2-D Cartesian wavefront sweep.

This module ships the **§15A.2 "upwind trace complex / causal transport
DAG / direction sweep ordering"** primitive (Grand Report v3 lines
2137-2171) for 2-D Cartesian SN, lifting the per-call ``_diag_cache``
precompute that previously lived inside
:func:`orpheus.sn.sweep._sweep_2d_wavefront` (lines 766-785) to
mesh-time work, and replacing the per-ordinate Python loop over
ordinates with a per-octant batched ``apply``.

The graph is a **derived object** — it depends only on cell topology
(``Mesh2D``) and the octant sign convention. It does NOT depend on
fluxes, sources, or iteration state, so it is built ONCE per
``(SNMesh, octant)`` pair at mesh construction (Wave 2 / C2.4) and
reused across every source iteration / Krylov matvec / outer iteration.

Architectural framing (Cardinal Rule 2)
=======================================

Per the project memory note ``project_moc_structure.md`` and
:doc:`../sn/spatial/cell_update`, this DAG abstraction is **SN-specific
by design**. MoC will define its own analog (per-ray traversal — fiber
bundles + solution sheaves, NOT a topological sort over a cell graph).
There is no shared ``SweepGraph`` Protocol because there is no shared
mathematical structure.

Tensorial framing (Wave 2 plan, §15A.2)
=======================================

Per octant, the streaming inverse :math:`L^{-1}_{\text{oct}}` acts on
the ``(octant_ordinates × cells × groups)`` tensor::

    L⁻¹_oct : (N_oct, nx, ny, ng) -> (N_oct, nx, ny, ng)

The anti-diagonal schedule vectorises ALL THREE axes inside one
``update_batch`` call — no Python loop over ordinates, no Python loop
over cells along a diagonal, no Python loop over groups. The outer
sweep iterates octants (4 in 2-D — structural) and topological levels
per octant (structural — sweep-DAG data dependency). The ordinate
axis is internal to every ``apply`` / ``einsum`` invocation.

References
==========

* Grand Report v3 §15A.2 lines 2137-2171 — "upwind trace complex /
  causal transport DAG / direction sweep ordering" primitive
  description with the ``assert_*`` invariant set this module's
  tests pin (``assert_upwind_orientation``,
  ``assert_topologically_sorted``, ``assert_face_pairing_consistent``,
  ``assert_boundary_trace_classification``,
  ``assert_cycles_are_declared``).
* Wave 2 plan ``.claude/plans/transient-giggling-cake.md`` C2.3 —
  this module's design + dispatch boundary with ``CellUpdate``.
* :class:`~orpheus.sn.spatial.cell_update.SweepCellSlice` — the
  per-level packet ``apply`` builds.
* :meth:`~orpheus.sn.spatial.cell_update.CellUpdateBase.update_batch`
  — the closure-pluggable consumer (Wave 2 / C2.2).

See also
========

* Future direction (out of scope): assemble per-octant lower-triangular
  ``L_oct`` as ``scipy.sparse.csr_matrix`` and replace the
  forward-substitution ``apply`` with
  ``scipy.sparse.linalg.spsolve_triangular``. Trades memory for a
  compiled inner kernel; worth it for very large 3-D problems. File a
  follow-up when the 3-D consumer arrives.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np

from orpheus.sn.spatial.cell_update import CellUpdateBase, SweepCellSlice


# ═══════════════════════════════════════════════════════════════════════
# OctantLabel — 2-D in-plane octant signature
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class OctantLabel:
    r"""Octant signature for the dimension-generic wavefront sweep.

    Carries one direction sign per spatial axis,
    ``signs[axis] ∈ {-1, 0, +1}``, so a single type labels a 1-D
    (``(±1,)``), 2-D (``(±1, ±1)``), or 3-D (``(±1, ±1, ±1)``) octant.
    The out-of-plane ``sign_z`` of an ``S²`` ordinate label may be
    dropped by the 2-D Cartesian orchestration (the in-plane sweep is
    invariant under it); multiple ordinates that project to the same
    in-plane ``signs`` share a single :class:`SweepDependencyGraph`.

    A label whose signs are *all* zero denotes the degenerate
    no-streaming set of ordinates (e.g. the pure-:math:`z` ordinates in
    2-D); the wavefront sweep handles those via a ``Q / Σ_t``
    short-circuit and builds no graph for it (:attr:`streams` is
    ``False``).

    Frozen + slotted: hashable for use as a ``dict`` key (the SNMesh
    holds ``Mapping[OctantLabel, SweepDependencyGraph]``).
    """

    signs: tuple[int, ...]

    def __post_init__(self) -> None:
        # Coerce any sequence to a tuple of ints (hashability + the
        # ``signs[axis]`` contract) BEFORE validating.
        object.__setattr__(self, "signs", tuple(int(s) for s in self.signs))
        for axis, s in enumerate(self.signs):
            if s not in (-1, 0, +1):
                raise ValueError(
                    f"signs[{axis}] must be in {{-1, 0, +1}}; got {s}"
                )

    @property
    def ndim(self) -> int:
        """Number of spatial axes this octant labels."""
        return len(self.signs)

    @property
    def streams(self) -> bool:
        """``False`` for the all-zero degenerate label; ``True`` otherwise."""
        return any(s != 0 for s in self.signs)

    # ── 2-D in-plane convenience: the orchestration that is still
    #    explicitly 2-D (the moving-frontier window walk, the matvec
    #    twin, the schedule) reads ``sign_x`` / ``sign_y`` and the
    #    ``streams_in_2d`` alias. These retire as those call sites go
    #    d-generic in a later C3 stage.
    @property
    def sign_x(self) -> int:
        """In-plane x-sign (``signs[0]``); valid for ``ndim ≥ 1``."""
        return self.signs[0]

    @property
    def sign_y(self) -> int:
        """In-plane y-sign (``signs[1]``); valid for ``ndim ≥ 2``."""
        return self.signs[1]

    @property
    def streams_in_2d(self) -> bool:
        """Deprecated alias for :attr:`streams` (2-D call-site compat)."""
        return self.streams


# ═══════════════════════════════════════════════════════════════════════
# _MovingFrontier — the rolling 2-diagonal interior face cochain (storage-B)
# ═══════════════════════════════════════════════════════════════════════


class _MovingFrontier:
    r"""The rolling 2-diagonal interior face cochain — storage-B's moving frontier.

    The interior per-ordinate face flux (the cochain :math:`C^1_{\rm int}`,
    :class:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux` in its
    full-field realization) need only be held on the ACTIVE anti-diagonal
    frontier during a sweep: a face produced at level ``k`` is consumed at
    ``k+1``, so two diagonals (ping-ponged by parity ``k % 2``) suffice. This
    shrinks the interior backing from ``O(N·ng·nx·ny)`` to ``O(N·ng·nx)``.

    The frontier is FAST *and* expressive because a diagonal's slots are
    CONTIGUOUS (``local_i ∈ [i0, i1]``): the gather is a basic-slice zero-copy
    VIEW and the advance a slice-assign — never the fancy 2-array index the
    full-field grid-diagonal walk is forced into (which copies every level).
    Measured ~0.77× the full-field walk time (a SPEEDUP, not a memory-vs-time
    trade), with the abstraction itself free (~1.00× vs the inline form).

    Layout (per octant, ``win_x`` carries a GHOST column 0 for the x-inflow):

    * ``win_x`` ``(N_oct, ng, 2, nx+1)`` — cell ``local_i`` READS x-column
      ``local_i`` (its upstream neighbour's x-outflow, OR the ghost for the
      ``local_i==0`` inflow edge) and WRITES column ``local_i+1``.
    * ``win_y`` ``(N_oct, ng, 2, nx)`` — cell ``local_i`` reads + writes column
      ``local_i`` across parities; the y-inflow edge cell's column is seeded so
      the gather slice includes it.

    The API is the cochain trace algebra: :meth:`seed_x` / :meth:`seed_y` are
    the :math:`\iota_*` inflow injection (per edge level), :meth:`incoming` the
    gather, :meth:`emit` the advance. The :math:`\iota^*` outflow shed reads the
    just-emitted ``out_x`` / ``out_y`` and is done by the caller (it targets the
    boundary trace, not the interior frontier).
    """

    __slots__ = ("win_x", "win_y")

    def __init__(self, N_oct: int, ng: int, nx: int) -> None:
        # +1 x-column = the ghost holding the per-level x-inflow.
        self.win_x = np.zeros((N_oct, ng, 2, nx + 1))
        self.win_y = np.zeros((N_oct, ng, 2, nx))

    def seed_x(self, prev: int, x_inflow_col: np.ndarray) -> None:
        r""":math:`\iota_*` — inject the x-inflow into the ghost column (read by
        the ``local_i==0`` edge cell as its upstream-x slot)."""
        self.win_x[:, :, prev, 0] = x_inflow_col

    def seed_y(self, prev: int, i1: int, y_inflow_col: np.ndarray) -> None:
        r""":math:`\iota_*` — inject the y-inflow into the edge cell's column
        ``i1`` (the level's last slot; the cell-below it does not exist)."""
        self.win_y[:, :, prev, i1] = y_inflow_col

    def incoming(
        self, prev: int, i0: int, i1: int,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Gather the level's incoming x/y face flux as basic-slice VIEWS
        (zero-copy): x-column ``local_i`` = ``[i0, i1]``, y-column ``[i0, i1]``."""
        return (
            self.win_x[:, :, prev, i0:i1 + 1],
            self.win_y[:, :, prev, i0:i1 + 1],
        )

    def emit(
        self, cur: int, i0: int, i1: int,
        out_x: np.ndarray, out_y: np.ndarray,
    ) -> None:
        """Advance the frontier: scatter the level's outgoing faces by
        slice-assign (x-column ``local_i+1`` = ``[i0+1, i1+1]``; y ``[i0, i1]``)."""
        self.win_x[:, :, cur, i0 + 1:i1 + 2] = out_x
        self.win_y[:, :, cur, i0:i1 + 1] = out_y


# ═══════════════════════════════════════════════════════════════════════
# SweepDependencyGraph — per-octant causal cell DAG (2-D Cartesian)
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class SweepDependencyGraph:
    r"""Per-octant causal cell DAG for the 2-D Cartesian wavefront sweep.

    The DAG topology is a regular pattern: cells :math:`(i, j)` and
    :math:`(i', j')` are mutually independent at the same topological
    level iff :math:`i + j = i' + j'` (under the canonical
    :math:`(+1, +1)` orientation; reversed indices for negative-sign
    octants). This collapses to ``nx + ny - 1`` topological levels,
    each carrying an anti-diagonal of cells whose updates are
    mutually independent.

    The closed-form precompute on a Cartesian grid lives entirely
    inside :meth:`from_cartesian_2d` and never appears in the sweep
    loop. This is structural, not hand-rolled — the "library
    version" (a generic topological-sort over an explicit DAG) would
    be over-engineering for a regular pattern that collapses to ~5
    lines of arange + diagonal extraction.

    The :meth:`apply` method does NOT inline the per-cell update
    math — it dispatches to ``cell_update.update_batch(slice_args)``,
    so the WDD / LD / EC / Step closure is pluggable. This aligns
    with the strategy contract in
    :mod:`orpheus.sn.spatial.cell_update`.

    Attributes
    ----------
    label :
        :class:`OctantLabel` carried for diagnostics + key matching.
    levels :
        Tuple of ``(ii, jj)`` pairs, one per topological level. Each
        is a pair of ``int`` arrays of equal length ``n_diag``
        carrying the cell-i / cell-j indices on that level. Walked
        in topological order.
    face_in_x, face_out_x, face_in_y, face_out_y :
        Per-axis face-index offsets. ``face_in`` is the offset added
        to ``ii`` (or ``jj``) to get the incoming-face index;
        ``face_out`` for the outgoing face. For ``sx ≥ 0``:
        ``face_in_x = 0`` (the "lower" x-face), ``face_out_x = 1``
        (the "upper"). For ``sx < 0``: reversed.

    Notes
    -----
    Frozen + (no slots — frozen suffices for hashability and
    immutability; ``slots=True`` would conflict with ``levels`` being
    a tuple-of-array, which numpy can store in slot but is awkward).

    Iteration over ``levels`` is the topological order of the upwind
    DAG — structural, not a smoking gun. Within each level the cell
    update is delegated to :meth:`CellUpdateBase.update_batch`, so
    the ordinate axis (``N_oct``) is internal to every ``apply``
    call: there is no ``for n in range(N_oct)`` anywhere in this
    module.
    """

    label: OctantLabel
    levels: tuple[tuple[np.ndarray, ...], ...]
    face_in: tuple[int, ...]
    face_out: tuple[int, ...]
    spatial_shape: tuple[int, ...]
    # ── Phase 5b storage-B: mesh-time moving-frontier window metadata ──
    # The d=2 OPTIMIZATION ONLY (``None`` for ``d != 2``): the d=1 cumprod
    # scan and the d≥3 full-field walk do not use the moving-frontier
    # window. The rolling 2-diagonal window walk keys the interior face
    # cochain by the cell's LOCAL sweep-order slot (``local_i ∈ [0, nx)``)
    # instead of its global face position — so the backing is O(N·ng·nx)
    # not O(N·ng·nx·ny). ``window_slots[k]`` is the per-level ``local_i``
    # array (= the slot to write ``win_*[k%2, slot]`` and read
    # ``win_y[(k-1)%2, slot]`` / ``win_x[(k-1)%2, slot-1]``); the companion
    # ``local_j`` is ``k - window_slots[k]`` (construction invariant), so
    # the domain-edge masks derive in the walk. ``window_edges[k]`` is the
    # mesh-time 4-tuple ``(has_x_in, has_x_out, has_y_in, has_y_out)`` of
    # Python bools (the touched cell is at a FIXED position — x-inflow at 0,
    # x-outflow at -1, y-inflow at -1, y-outflow at 0 — so the walk guards
    # each edge op with a precomputed bool). Derived mesh-time (no flux
    # dependence) — L16: zero per-sweep recompute.
    window_slots: tuple[np.ndarray, ...] | None
    window_edges: tuple[tuple[bool, ...], ...] | None

    # ── 2-D in-plane convenience accessors (compat; retire with the
    #    d-generic orchestration). The window walk, the matvec twin, and
    #    the legacy-equivalence test read these. ──
    @property
    def ndim(self) -> int:
        """Number of spatial axes (``len(spatial_shape)``)."""
        return len(self.spatial_shape)

    @property
    def nx(self) -> int:
        return self.spatial_shape[0]

    @property
    def ny(self) -> int:
        return self.spatial_shape[1]

    @property
    def face_in_x(self) -> int:
        return self.face_in[0]

    @property
    def face_out_x(self) -> int:
        return self.face_out[0]

    @property
    def face_in_y(self) -> int:
        return self.face_in[1]

    @property
    def face_out_y(self) -> int:
        return self.face_out[1]

    # ── Construction ────────────────────────────────────────────────

    @classmethod
    def from_cartesian(
        cls,
        shape: tuple[int, ...],
        *,
        label: OctantLabel,
    ) -> "SweepDependencyGraph":
        r"""Build the per-octant upwind DAG for a ``d``-dim Cartesian grid.

        Dimension-generic (``d = len(shape) ∈ {1, 2, 3}``). The cells on
        topological level ``k`` are those whose per-octant local indices
        satisfy :math:`\sum_a \mathrm{local}_a = k` — the anti-hyperplane
        of the index lattice. There are :math:`\sum_a (n_a - 1) + 1`
        levels. Per axis ``a``, ``local_a`` walks ``0 → n_a-1`` forward
        when ``signs[a] ≥ 0`` and reversed when ``signs[a] < 0``; the
        local index translates to the global cell index via
        ``axis_map[a][local_a]``.

        Within a level the cells are ordered **C-major** over the index
        lattice (axis 0 slowest). At ``d = 2`` this reproduces the legacy
        ascending-``local_i`` anti-diagonal order **bit-for-bit**, so a
        2-D graph built here is identical to the retired
        ``from_cartesian_2d`` (pinned by
        ``test_d2_from_cartesian_matches_legacy``). At ``d = 1`` the DAG
        is a pure chain (level ``ℓ`` holds only cell ``ℓ``); at ``d = 3``
        each level is a 2-D lattice slice of the ``i+j+k = ℓ`` hyperplane.

        Parameters
        ----------
        shape :
            Per-axis cell counts ``(n_0, …, n_{d-1})``.
        label :
            Octant signature with ``len(signs) == d`` and at least one
            non-zero sign (the DAG is undefined for the all-zero
            degenerate label — building one raises ``ValueError``).

        Raises
        ------
        ValueError
            If ``label.streams`` is ``False`` (all signs zero), or if
            ``len(label.signs) != len(shape)``.
        """
        if not label.streams:
            raise ValueError(
                f"Cannot build a SweepDependencyGraph for the all-zero "
                f"degenerate octant {label!r}; it has no streaming and the "
                "wavefront sweep handles it via a Q/Σ_t short-circuit."
            )
        shape = tuple(int(n) for n in shape)
        signs = label.signs
        if len(signs) != len(shape):
            raise ValueError(
                f"label ndim {len(signs)} != shape ndim {len(shape)} "
                f"(signs={signs}, shape={shape})"
            )
        d = len(shape)

        # Per-axis face-index offsets (a cell's incoming / outgoing face
        # along each axis). Forward orientation puts the incoming face at
        # the lower index (0); a negative sign flips it.
        face_in = tuple(0 if s >= 0 else 1 for s in signs)
        face_out = tuple(1 if s >= 0 else 0 for s in signs)

        # Per-axis local→global index map: forward arange for s ≥ 0,
        # reversed for s < 0 (the per-octant sweep direction).
        axis_map = [
            np.arange(n) if s >= 0 else np.arange(n)[::-1]
            for n, s in zip(shape, signs)
        ]

        # Anti-hyperplane levels: enumerate the local index lattice in
        # C-major order, group by index-sum. Boolean masking preserves the
        # C-major within-level order ⟹ d=2 bit-identical to legacy. Each
        # level is the ndim-tuple of equal-length global-index arrays.
        local = np.indices(shape).reshape(d, -1)        # (d, prod), C-major
        level_of = local.sum(axis=0)                    # (prod,)
        n_levels = sum(n - 1 for n in shape) + 1
        levels: tuple[tuple[np.ndarray, ...], ...] = tuple(
            tuple(axis_map[a][local[a, level_of == k]] for a in range(d))
            for k in range(n_levels)
        )

        window_slots, window_edges = cls._window_metadata(shape, n_levels)

        return cls(
            label=label,
            levels=levels,
            face_in=face_in,
            face_out=face_out,
            spatial_shape=shape,
            window_slots=window_slots,
            window_edges=window_edges,
        )

    @staticmethod
    def _window_metadata(
        shape: tuple[int, ...], n_levels: int,
    ) -> tuple[
        tuple[np.ndarray, ...] | None, tuple[tuple[bool, ...], ...] | None
    ]:
        """Storage-B moving-frontier window metadata — the d=2 optimization.

        Returns ``(None, None)`` for ``d != 2`` (the d=1 cumprod scan and
        the d≥3 full-field walk do not use the moving-frontier window).
        For ``d == 2`` this is the legacy per-level ``local_i`` slot array
        + the ``(has_x_in, has_x_out, has_y_in, has_y_out)`` edge bools,
        built exactly as the retired ``from_cartesian_2d`` did (so the
        window walk stays bit-identical). The slot ``local_i`` is
        contiguous ⟹ basic-slice addressable; the edge bools mark the
        FIXED domain-edge positions (x-inflow at 0, x-outflow at -1,
        y-inflow at -1, y-outflow at 0) so the walk guards each edge op
        without a per-level mask.
        """
        if len(shape) != 2:
            return None, None
        nx, ny = shape
        window_slots: list[np.ndarray] = []
        window_edges: list[tuple[bool, ...]] = []
        for k in range(n_levels):
            i_start = max(0, k - ny + 1)
            i_end = min(nx - 1, k)
            window_slots.append(np.arange(i_start, i_end + 1))
            window_edges.append((
                i_start == 0,            # has_x_in  (cell at position 0)
                i_end == nx - 1,         # has_x_out (cell at position -1)
                i_end == k,              # has_y_in  (cell at position -1)
                k >= ny - 1,             # has_y_out (cell at position 0)
            ))
        return tuple(window_slots), tuple(window_edges)

    @classmethod
    def from_cartesian_2d(
        cls,
        *,
        nx: int,
        ny: int,
        label: OctantLabel,
    ) -> "SweepDependencyGraph":
        """Deprecated 2-D alias for :meth:`from_cartesian`.

        Retained one cycle for the 2-D call sites + the legacy-equivalence
        pin; retires when the geometry octant build switches to
        ``from_cartesian`` + ``itertools.product`` (C3.3).
        """
        return cls.from_cartesian((nx, ny), label=label)

    # ── Internal: shared per-level slice builder ────────────────────

    def _make_slice(
        self,
        ii: np.ndarray,
        jj: np.ndarray,
        *,
        psi_x: np.ndarray,
        psi_y: np.ndarray,
        Q: np.ndarray,
        sig_t: np.ndarray,
        str_x: np.ndarray,
        str_y: np.ndarray,
        psi_avg_probe: np.ndarray | None = None,
    ) -> SweepCellSlice:
        """Build the per-level :class:`SweepCellSlice` for this graph's
        face-offset convention.

        Shared by :meth:`apply` (solve → ``update_batch``) and
        :meth:`residual` (apply → ``residual_batch``) so the level +
        face-index wiring (``face_in_x`` / ``face_out_x`` / …) is a
        single source of truth — the two walks differ only in their
        accumulation, not in the slice they build. ``psi_avg_probe`` is
        ``None`` for the solve direction (``update_batch`` ignores it).
        """
        return SweepCellSlice(
            ii=ii, jj=jj,
            face_in_x_idx=ii + self.face_in_x,
            face_out_x_idx=ii + self.face_out_x,
            face_in_y_idx=jj + self.face_in_y,
            face_out_y_idx=jj + self.face_out_y,
            psi_x=psi_x, psi_y=psi_y,
            Q=Q, sig_t=sig_t, str_x=str_x, str_y=str_y,
            psi_avg_probe=psi_avg_probe,
        )

    # ── Public API ─────────────────────────────────────────────────

    def apply(
        self,
        *,
        cell_update: CellUpdateBase,
        psi_x_octant: np.ndarray,        # (N_oct, ng, nx+1, ny) — mutated in place
        psi_y_octant: np.ndarray,        # (N_oct, ng, nx, ny+1) — mutated in place
        Q_octant: np.ndarray,            # (N_oct or 1, ng, nx, ny)
        sig_t: np.ndarray,               # (ng, nx, ny)
        str_x_octant: np.ndarray,        # (N_oct, nx)
        str_y_octant: np.ndarray,        # (N_oct, ny)
        weights_octant: np.ndarray,      # (N_oct,)
        angular_flux_octant: np.ndarray, # (N_oct, ng, nx, ny) — written
        scalar_flux_buf: np.ndarray,     # (ng, nx, ny) — accumulated into
    ) -> None:
        r"""Walk the topological levels and accumulate scalar + angular flux.

        Forward-substitute along the DAG, vectorised over
        ``(N_oct, n_diag, ng)`` per level. The ordinate axis
        (``N_oct``) is INTERNAL to every ``apply`` call: there is no
        ``for n in range(N_oct)`` anywhere in this method.

        Mutation contract
        -----------------

        * ``psi_x_octant`` / ``psi_y_octant`` — outgoing face fluxes
          are scattered into these buffers in place by
          :meth:`CellUpdateBase.update_batch`. Caller is responsible
          for seeding the incoming-face entries (BC apply happens
          one level up, in the wavefront sweep) and for scattering
          the post-sweep buffers back into the persistent BC state.
        * ``angular_flux_octant`` — written at every level's cells.
          Caller is responsible for scattering back into the global
          ``angular_flux`` buffer keyed by octant indices.
        * ``scalar_flux_buf`` — accumulated into via
          ``scalar_flux_buf[ii, jj, :] += weighted_sum_psi_avg``.
          Caller seeds the buffer (typically zero on the first
          octant of an outer iteration) and reads the result at the
          end of all octants.

        Parameters
        ----------
        cell_update :
            The closure strategy. Must override
            :meth:`update_batch` (DD does; Step / LD do not yet).
        psi_x_octant, psi_y_octant :
            Per-octant face buffers; mutated in place.
        Q_octant :
            Per-octant per-cell volumetric source, already
            weight-normalised. Leading axis is ``1`` for
            isotropic-only sweeps or ``N_oct`` when an aniso
            component is folded in.
        sig_t :
            Per-cell per-group total cross section.
        str_x_octant, str_y_octant :
            Per-octant streaming coefficients ``2|μ_axis|/Δaxis``.
        weights_octant :
            Per-octant ordinate weights for scalar-flux reduction.
        angular_flux_octant :
            Per-octant angular-flux output buffer.
        scalar_flux_buf :
            Global scalar-flux accumulator.
        """
        for ii, jj in self.levels:
            slice_args = self._make_slice(
                ii, jj,
                psi_x=psi_x_octant, psi_y=psi_y_octant,
                Q=Q_octant, sig_t=sig_t,
                str_x=str_x_octant, str_y=str_y_octant,
            )
            # Issue #196 PR-INDEX-5: psi_avg returned principled
            # ``(N_oct, ng, n_diag)`` (ordinate, group, anti-diagonal).
            psi_avg = cell_update.update_batch(slice_args)  # (N_oct, ng, n_diag)
            # Scalar-flux accumulation — sum over ordinates with
            # quadrature weights. ``weights_octant`` is already
            # weight-normalised by the caller (the wavefront sweep
            # divides by Σw before invoking apply), so this is a
            # plain weighted sum.
            scalar_flux_buf[:, ii, jj] += np.einsum(
                "ngd,n->gd", psi_avg, weights_octant,
            )
            # Angular-flux scatter — write all ``(N_oct, ng, n_diag)``
            # values into the per-octant principled angular-flux buffer.
            # Indices ``ii`` and ``jj`` broadcast naturally with the
            # leading ``:`` (octant) and the second ``:`` (group).
            angular_flux_octant[:, :, ii, jj] = psi_avg

    def residual(
        self,
        *,
        cell_update: CellUpdateBase,
        psi_x_octant: np.ndarray,            # (N_oct, ng, nx+1, ny) — mutated in place
        psi_y_octant: np.ndarray,            # (N_oct, ng, nx, ny+1) — mutated in place
        psi_avg_probe_octant: np.ndarray,    # (N_oct, ng, nx, ny) — the probe (read)
        Q_octant: np.ndarray,                # (N_oct or 1, ng, nx, ny)
        sig_t: np.ndarray,                   # (ng, nx, ny)
        str_x_octant: np.ndarray,            # (N_oct, nx)
        str_y_octant: np.ndarray,            # (N_oct, ny)
        residual_octant: np.ndarray,         # (N_oct, ng, nx, ny) — written
    ) -> None:
        r"""Walk the topological levels accumulating the operator residual.

        The **apply-direction** companion of :meth:`apply`: where ``apply``
        forward-substitutes the *solve* (``update_batch``) and reduces to
        scalar + angular flux, ``residual`` forward-substitutes the *apply*
        (:meth:`CellUpdateBase.residual_batch`) and writes the per-cell
        operator residual :math:`(L+C)\,\overline\psi - q` into
        ``residual_octant``. The 2-D Cartesian matvec
        ``StreamingOperator.apply`` (Wave O #208 O.4b) drives this so it
        shares the SAME wavefront DAG + closure as the sweep — matvec and
        sweep are one discretization (L21), no FD twin path.

        The edge-flux reconstruction is load-bearing: each cell's incoming
        face flux is the upstream cell's outgoing flux via the diamond
        closure :math:`\psi^{\rm out} = 2\overline\psi^{\rm probe} -
        \psi^{\rm in}` (scattered into ``psi_x_octant`` / ``psi_y_octant``
        by ``residual_batch`` in place, propagated from the seeded boundary
        inflow along the wavefront). This is exactly why the matvec needs
        the topological walk rather than a per-cell formula.

        Mutation contract
        -----------------

        * ``psi_x_octant`` / ``psi_y_octant`` — outgoing face fluxes
          scattered in place (from the probe). Caller seeds the
          incoming-face entries (BC apply happens one level up, in the
          matvec) exactly as for :meth:`apply`.
        * ``residual_octant`` — written at every level's cells; caller
          scatters it back into the global bulk-residual buffer keyed by
          octant indices.

        Parameters mirror :meth:`apply`, with ``psi_avg_probe_octant``
        (the apply target) replacing the flux-reduction outputs
        (``weights_octant`` / ``scalar_flux_buf`` / ``angular_flux_octant``).
        """
        for ii, jj in self.levels:
            slice_args = self._make_slice(
                ii, jj,
                psi_x=psi_x_octant, psi_y=psi_y_octant,
                Q=Q_octant, sig_t=sig_t,
                str_x=str_x_octant, str_y=str_y_octant,
                psi_avg_probe=psi_avg_probe_octant,
            )
            residual_octant[:, :, ii, jj] = cell_update.residual_batch(slice_args)

    # ── Phase 5b storage-B: rolling moving-frontier window walks ────────
    #
    # The PRODUCTION walks. They carry the interior face cochain on a
    # :class:`_MovingFrontier` (rolling 2-diagonal window, O(N·ng·nx) not
    # O(N·ng·nx·ny)) and ADVANCE it level by level — basic-slice views/assigns
    # over the contiguous slots, FASTER than the full-field grid-diagonal walk
    # (~0.77×). The cell math is the SAME shared kernel
    # (DiamondDifference.cell_kernel_batch / residual_kernel_batch) the
    # full-field walks (apply / residual) use — so window and full-field cannot
    # drift mathematically (proven bit-identical by the window≡full oracle test).
    # Domain-edge inflow is ι_*-seeded onto the frontier per level (from the
    # passed inflow arrays); domain-edge outflow is ι*-shed into the capture
    # arrays at the level that produces it, BEFORE its frontier slot recycles.

    def _bounds(self, slot: np.ndarray) -> tuple[int, int]:
        """Contiguous local-i bounds ``(i0, i1)`` of a level's slot range."""
        return int(slot[0]), int(slot[-1])

    def apply_windowed(
        self,
        *,
        cell_update: CellUpdateBase,
        inflow_x: np.ndarray,            # (N_oct, ng, ny) — domain x-inflow by global y-cell
        inflow_y: np.ndarray,            # (N_oct, ng, nx) — domain y-inflow by global x-cell
        Q_octant: np.ndarray,            # (N_oct or 1, ng, nx, ny)
        sig_t: np.ndarray,               # (ng, nx, ny)
        str_x_octant: np.ndarray,        # (N_oct, nx)
        str_y_octant: np.ndarray,        # (N_oct, ny)
        weights_octant: np.ndarray,      # (N_oct,)
        capture_x: np.ndarray,           # (N_oct, ng, ny) — domain x-outflow, written
        capture_y: np.ndarray,           # (N_oct, ng, nx) — domain y-outflow, written
        angular_flux_octant: np.ndarray | None = None,  # (N_oct, ng, nx, ny) — written (angular mode)
        scalar_flux_buf: np.ndarray | None = None,       # (ng, nx, ny) — accumulated (angular mode)
        moment_buf: np.ndarray | None = None,            # (L+1, 2L+1, ng, nx, ny) — accumulated (moment mode)
        Y_octant: np.ndarray | None = None,              # (N_oct, L+1, 2L+1) — octant harmonics (moment mode)
    ) -> None:
        r"""Storage-B solve walk: advance a :class:`_MovingFrontier` along the
        anti-diagonals, sheding domain-edge outflow as it goes.

        The window-backed twin of :meth:`apply`. The interior face cochain + the
        cell math are bit-identical to it (same shared cell kernel, same
        anti-diagonal order); proven by the ``window ≡ full-field`` oracle test.
        The boundary inflow is ι_*-seeded onto the frontier from ``inflow_x`` /
        ``inflow_y`` (the octant's incoming domain-edge trace); the outflow is
        ι*-shed into ``capture_x`` / ``capture_y``.

        Output representation (Phase 5c) — exactly ONE of two modes, selected by
        which output buffer is given (dependency injection, mirroring the
        ``reflect=None`` idiom in :func:`_sweep_2d_scheduled`):

        * **angular** (``moment_buf is None``) — write the full per-ordinate
          ``angular_flux_octant[:, :, ii, jj] = psi_avg`` per level and
          accumulate the scalar ``scalar_flux_buf``. The full-angular OUTPUT a
          reconstruction / Krylov / full-angular SI consumer needs.
        * **moment** (``moment_buf`` given) — project per anti-diagonal directly
          into the harmonic moment tensor
          :math:`\phi_\ell^m \mathrel{+}= \sum_n w_n Y_\ell^m(\hat\Omega_n)
          \psi_n` (``Y_octant`` the octant's harmonics), NEVER materializing the
          full angular field (the Phase 5c ~3× linear peak-memory win — the
          persistent windowed SI iterate is already moments, 5a). The scalar is
          subsumed (``moment_buf[0, 0]`` IS :math:`\phi_0^0`, since
          :math:`Y_0^0=1`), so no separate ``scalar_flux_buf`` accumulation. The
          cross-octant ``+=`` reorders the ordinate sum vs the post-sweep flat
          :class:`~orpheus.numerics.projection.MomentProjection` reduce ⇒
          principled-equivalence, NOT bit-identity (``vv-principles``
          §"Bit-identity vs principled-equivalence"; de-risk ≤ 4 ULP). The
          per-cell ``w·Y·psi`` fold matches ``MomentProjection.apply``
          term-for-term — only the accumulation order differs.
        """
        if self.window_slots is None or self.window_edges is None:
            raise RuntimeError(
                "apply_windowed is the d=2 storage-B optimization; this graph "
                f"has ndim={self.ndim} (no window metadata) — use the "
                "full-field apply() for d != 2."
            )
        N_oct, ng = inflow_x.shape[0], inflow_x.shape[1]
        frontier = _MovingFrontier(N_oct, ng, self.nx)
        for k, ((ii, jj), slot, edges) in enumerate(zip(
            self.levels, self.window_slots, self.window_edges,
        )):
            has_x_in, has_x_out, has_y_in, has_y_out = edges
            i0, i1 = self._bounds(slot)
            cur, prev = k % 2, (k - 1) % 2
            if has_x_in:
                frontier.seed_x(prev, inflow_x[:, :, jj[0]])
            if has_y_in:
                frontier.seed_y(prev, i1, inflow_y[:, :, ii[-1]])
            in_x, in_y = frontier.incoming(prev, i0, i1)
            psi_avg, out_x, out_y = cell_update.cell_kernel_batch(
                psi_in_x=in_x, psi_in_y=in_y,
                sx=str_x_octant[:, ii][:, None, :],
                sy=str_y_octant[:, jj][:, None, :],
                sigt_cells=sig_t[:, ii, jj], Q_cells=Q_octant[:, :, ii, jj],
            )
            frontier.emit(cur, i0, i1, out_x, out_y)
            if has_x_out:
                capture_x[:, :, jj[-1]] = out_x[:, :, -1]   # ι* x-outflow @ pos -1
            if has_y_out:
                capture_y[:, :, ii[0]] = out_y[:, :, 0]     # ι* y-outflow @ pos 0
            if moment_buf is None:
                scalar_flux_buf[:, ii, jj] += np.einsum(
                    "ngd,n->gd", psi_avg, weights_octant,
                )
                angular_flux_octant[:, :, ii, jj] = psi_avg
            else:
                moment_buf[:, :, :, ii, jj] += np.einsum(
                    "nlm,ngd,n->lmgd", Y_octant, psi_avg, weights_octant,
                )

    def residual_windowed(
        self,
        *,
        cell_update: CellUpdateBase,
        inflow_x: np.ndarray,             # (N_oct, ng, ny)
        inflow_y: np.ndarray,             # (N_oct, ng, nx)
        psi_avg_probe_octant: np.ndarray, # (N_oct, ng, nx, ny) — the probe
        Q_octant: np.ndarray,             # (N_oct or 1, ng, nx, ny)
        sig_t: np.ndarray,                # (ng, nx, ny)
        str_x_octant: np.ndarray,         # (N_oct, nx)
        str_y_octant: np.ndarray,         # (N_oct, ny)
        residual_octant: np.ndarray,      # (N_oct, ng, nx, ny) — written
        capture_x: np.ndarray,            # (N_oct, ng, ny) — written
        capture_y: np.ndarray,            # (N_oct, ng, nx) — written
    ) -> None:
        r"""Storage-B apply (matvec) walk: the :class:`_MovingFrontier`-backed
        twin of :meth:`residual`. Reconstructs edge fluxes from the PROBE along
        the rolling frontier (``psi_out = 2·psi_bar − psi_in``) and writes the
        per-cell operator residual; sheds the domain-edge outflow. The matvec
        residual output stays full ``(N_oct, ng, nx, ny)`` (Krylov consumes
        the full bulk residual) — only the interior FACE cochain is windowed.
        """
        if self.window_slots is None or self.window_edges is None:
            raise RuntimeError(
                "residual_windowed is the d=2 storage-B optimization; this "
                f"graph has ndim={self.ndim} (no window metadata) — use the "
                "full-field residual() for d != 2."
            )
        N_oct, ng = inflow_x.shape[0], inflow_x.shape[1]
        frontier = _MovingFrontier(N_oct, ng, self.nx)
        for k, ((ii, jj), slot, edges) in enumerate(zip(
            self.levels, self.window_slots, self.window_edges,
        )):
            has_x_in, has_x_out, has_y_in, has_y_out = edges
            i0, i1 = self._bounds(slot)
            cur, prev = k % 2, (k - 1) % 2
            if has_x_in:
                frontier.seed_x(prev, inflow_x[:, :, jj[0]])
            if has_y_in:
                frontier.seed_y(prev, i1, inflow_y[:, :, ii[-1]])
            in_x, in_y = frontier.incoming(prev, i0, i1)
            res, out_x, out_y = cell_update.residual_kernel_batch(
                psi_bar=psi_avg_probe_octant[:, :, ii, jj],
                psi_in_x=in_x, psi_in_y=in_y,
                sx=str_x_octant[:, ii][:, None, :],
                sy=str_y_octant[:, jj][:, None, :],
                sigt_cells=sig_t[:, ii, jj], Q_cells=Q_octant[:, :, ii, jj],
            )
            frontier.emit(cur, i0, i1, out_x, out_y)
            if has_x_out:
                capture_x[:, :, jj[-1]] = out_x[:, :, -1]
            if has_y_out:
                capture_y[:, :, ii[0]] = out_y[:, :, 0]
            residual_octant[:, :, ii, jj] = res


__all__ = [
    "OctantLabel",
    "SweepDependencyGraph",
]
