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

from dataclasses import dataclass

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
# _FrontierPlan / _MovingFrontier — the rolling (d−1)-frontier face cochain
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class _LevelFrontier:
    r"""Mesh-time slab addressing for ONE anti-hyperplane level (storage-B).

    All :math:`d`-dependent index arithmetic for one level of the rolling
    :math:`(d{-}1)`-frontier window, precomputed so the walk
    (:meth:`SweepDependencyGraph.apply_windowed` /
    :meth:`~SweepDependencyGraph.residual_windowed`) is dimension-agnostic
    (L16: zero per-sweep recompute).

    The frontier slab is indexed by the :math:`(d{-}1)` FREE *local*
    coordinates (the first :math:`d{-}1` axes); the determined axis (the last)
    is the parity-ping-ponged one.

    Attributes
    ----------
    read :
        The slab gather selector — the level's cells' free-coordinate
        positions, a :class:`slice` per coordinate when the free region is
        box-contiguous (``d ≤ 2`` — a 1-D anti-hyperplane is an interval) and
        fancy index arrays when it is a simplex (``d ≥ 3``).  Applied as
        ``win[a][:, :, prev, *read] → (N_oct, ng, n_diag)``.
    write :
        Per-axis slab scatter selector.  For a FREE axis it is ``read`` shifted
        ``+1`` on that axis's own coordinate (cell ``f`` writes its high-face
        for downstream cell ``f + e_a``); for the DETERMINED axis it equals
        ``read`` (the progression is the parity roll).
    seed :
        Per-axis domain-edge inflow op: ``None`` (no edge cells this level) or
        ``(slab_target, inflow_source)`` — scatter ``inflow[a][*inflow_source]``
        into ``win[a][prev][*slab_target]`` before the gather.
    shed :
        Per-axis domain-edge outflow op: ``None`` or ``(out_mask, capture_target)``
        — gather the kernel's outgoing ``out_faces[a][..., out_mask]`` into
        ``capture[a][*capture_target]`` after the kernel.
    """

    read: tuple
    write: tuple
    seed: tuple
    shed: tuple


@dataclass(frozen=True)
class _FrontierPlan:
    r"""The whole-sweep rolling :math:`(d{-}1)`-frontier window plan (storage-B).

    A face produced at level ``k`` is consumed at ``k+1``, so the interior face
    cochain (the cochain :math:`C^1_{\rm int}`,
    :class:`~orpheus.transport.fields.wavefront_flux.WavefrontFlux` in its
    full-field realization) need only be held on the two active levels — a
    rolling slab ping-ponged by parity ``k % 2``.  The slab is the
    :math:`(d{-}1)`-dim free bounding box, shrinking the interior backing from
    ``O(N·ng·∏ n_a)`` to ``O(N·ng·∏_{a<d−1} n_a)`` (the ``~3×`` peak-memory win
    Phase 5b measured at d=2; it grows with the angular order and generalises
    to any ``d``).

    The frontier is the :math:`(d{-}1)`-dim rolling slab: a *point* at d=1
    (``frontier_dim == 0``, the degenerate base), a *line* at d=2, a *surface*
    at d=3.  At d=2 the per-level :attr:`~_LevelFrontier.read` selector is a
    contiguous :class:`slice` ⟹ the gather is a basic-slice zero-copy VIEW and
    the advance a slice-assign — the ``~0.77×`` contiguity SPEEDUP (a free win,
    not a memory-vs-time trade).  At d≥3 the anti-hyperplane is a simplex (not a
    box), so the selector is a fancy index (copies) — the memory win still
    holds; whether the speed win survives the d=3 *surface* is the one
    measured-cost question deferred to profiling (no 3-D quadrature yet; the
    correctness gate is the synthetic ``window ≡ full`` admission).

    The slab per axis ``a``: a FREE axis (``a < det``) carries a ``+1`` ghost
    along its own coordinate (cell ``f`` reads its low-``a`` face at ``f``, the
    ghost holding the domain inflow for the ``f_a == 0`` edge, and writes its
    high-``a`` face at ``f + e_a``); the DETERMINED axis (``a == det``) carries
    no ghost (read coordinate == write coordinate, the progression is parity).
    """

    spatial_shape: tuple[int, ...]
    free_bbox: tuple[int, ...]
    det: int
    levels: tuple[_LevelFrontier, ...]

    @property
    def is_point(self) -> bool:
        """``True`` for the d=1 degenerate frontier (a point — no free axes)."""
        return len(self.free_bbox) == 0


class _MovingFrontier:
    r"""The rolling :math:`(d{-}1)`-frontier interior face cochain (storage-B).

    Holds one slab per spatial axis (a FREE axis ghosted ``+1`` on its own
    coordinate; the DETERMINED axis plain, parity-rolled), driven by a mesh-time
    :class:`_FrontierPlan`.  The ping-pong is valid at every ``d``: a face
    produced at level ``k`` (parity ``k % 2``) is read at level ``k+1`` (parity
    ``(k+1−1) % 2 == k % 2``), and every read hits a slab position the
    immediately-prior level wrote (cell ``c`` reads its low-``a`` face from the
    upstream cell ``c − e_a``, which is on level ``k−1``) or a seeded domain
    inflow.

    The API is the cochain trace algebra: :meth:`seed` is the :math:`\iota_*`
    inflow injection (per edge level), :meth:`incoming` the gather, :meth:`emit`
    the advance, :meth:`shed` the :math:`\iota^*` domain-outflow capture.
    """

    __slots__ = ("_win", "_plan")

    def __init__(self, N_oct: int, ng: int, plan: _FrontierPlan) -> None:
        free = plan.free_bbox
        det = plan.det
        win = []
        for a in range(len(plan.spatial_shape)):
            shp = list(free)
            if a < det:                          # FREE axis: +1 ghost on coord a
                shp[a] = plan.spatial_shape[a] + 1
            win.append(np.zeros((N_oct, ng, 2, *shp)))
        self._win = tuple(win)
        self._plan = plan

    def seed(self, prev: int, lvl: int, inflow: tuple[np.ndarray, ...]) -> None:
        r""":math:`\iota_*` — inject each axis's domain inflow into the slab at
        this level's edge cells (so the subsequent :meth:`incoming` gather reads
        them).  A no-op for an axis with no edge cell on this level."""
        for a, op in enumerate(self._plan.levels[lvl].seed):
            if op is None:
                continue
            slab_target, inflow_source = op
            self._win[a][(slice(None), slice(None), prev) + slab_target] = (
                inflow[a][(slice(None), slice(None)) + inflow_source]
            )

    def incoming(self, prev: int, lvl: int) -> tuple[np.ndarray, ...]:
        """Gather the level's incoming face flux per axis (a zero-copy VIEW at
        d≤2 via the contiguous slice; a fancy-index copy at d≥3)."""
        read = self._plan.levels[lvl].read
        faces = []
        for a in range(len(self._win)):
            g = self._win[a][(slice(None), slice(None), prev) + read]
            if self._plan.is_point:              # d=1: restore the n_diag axis
                g = g[:, :, None]
            faces.append(g)
        return tuple(faces)

    def emit(
        self, cur: int, lvl: int, out_faces: tuple[np.ndarray, ...],
    ) -> None:
        """Advance the frontier: scatter each axis's outgoing face flux into the
        downstream slot (free axis ``f + e_a``; determined axis ``f``)."""
        write = self._plan.levels[lvl].write
        for a in range(len(self._win)):
            out = out_faces[a][:, :, 0] if self._plan.is_point else out_faces[a]
            self._win[a][(slice(None), slice(None), cur) + write[a]] = out

    def shed(
        self, lvl: int, out_faces: tuple[np.ndarray, ...],
        capture: tuple[np.ndarray, ...],
    ) -> None:
        r""":math:`\iota^*` — capture each axis's domain outflow (the high-edge
        cells' outgoing faces) into the per-axis ``capture`` arrays.  A no-op for
        an axis with no outflow-edge cell on this level."""
        for a, op in enumerate(self._plan.levels[lvl].shed):
            if op is None:
                continue
            out_mask, capture_target = op
            # d=1 degenerate point: the single edge cell, n_diag axis squeezed
            # (the domain outflow has no perpendicular coordinate).
            shed_faces = (
                out_faces[a][:, :, 0] if self._plan.is_point
                else out_faces[a][:, :, out_mask]
            )
            capture[a][(slice(None), slice(None)) + capture_target] = shed_faces


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
    inside :meth:`from_cartesian` and never appears in the sweep
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
    # ── Phase 5b storage-B: mesh-time rolling (d−1)-frontier window plan ──
    # The production-optimization addressing for the ``apply_windowed`` /
    # ``residual_windowed`` walks (the ``MovingFrontierWindow`` strategy): the
    # interior face cochain is held on a rolling (d−1)-frontier slab, keyed by
    # the cell's LOCAL sweep-order free coordinates rather than its global face
    # position — so the backing is the (d−1)-slab ``O(N·ng·∏_{a<d−1} n_a)`` not
    # the full ``O(N·ng·∏ n_a)``.  All d-dependent index arithmetic is
    # precomputed here (the per-level read/write selectors + the domain-edge
    # seed/shed index maps) so the walk is dimension-agnostic — L16: zero
    # per-sweep recompute.  Built for every streaming Cartesian graph (d ≥ 1);
    # the d=1 instance is the degenerate ``frontier_dim == 0`` point (the window
    # is not the production default at d=1 — the cumprod scan is — but the plan
    # is built so the window's d=1 capability is verifiable, the governing
    # principle "construct general, select narrow").
    window_plan: _FrontierPlan

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
        2-D graph built here matches the legacy 2-D anti-diagonal
        construction (pinned by ``test_d2_from_cartesian_matches_legacy``,
        a hand-derived frozen golden). At ``d = 1`` the DAG
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
        # Per level: the LOCAL sweep-order coords (for the window slab) AND the
        # GLOBAL cell coords (axis_map applied; for the source / XS / flux /
        # boundary indexing).  Both share the C-major within-level order.
        local_levels = tuple(
            tuple(local[a, level_of == k] for a in range(d))
            for k in range(n_levels)
        )
        levels: tuple[tuple[np.ndarray, ...], ...] = tuple(
            tuple(axis_map[a][loc] for a, loc in enumerate(local_lvl))
            for local_lvl in local_levels
        )

        window_plan = cls._build_frontier_plan(shape, levels, local_levels)

        return cls(
            label=label,
            levels=levels,
            face_in=face_in,
            face_out=face_out,
            spatial_shape=shape,
            window_plan=window_plan,
        )

    @staticmethod
    def _build_frontier_plan(
        shape: tuple[int, ...],
        levels: tuple[tuple[np.ndarray, ...], ...],
        local_levels: tuple[tuple[np.ndarray, ...], ...],
    ) -> _FrontierPlan:
        r"""Build the rolling :math:`(d{-}1)`-frontier window plan (storage-B).

        Precomputes, for every anti-hyperplane level, the slab read/write
        selectors + the domain-edge seed/shed index maps — all the
        :math:`d`-dependent index arithmetic, so the
        :meth:`apply_windowed` / :meth:`residual_windowed` walks are
        dimension-agnostic (L16: zero per-sweep recompute).  ``levels`` carries
        the GLOBAL cell coords (per-octant sweep direction), ``local_levels``
        the LOCAL sweep-order coords (the slab is local-indexed) — both in the
        same C-major within-level order.

        The determined axis is the LAST (``det = d − 1``); the free coordinates
        are the first :math:`d{-}1`.  At ``d ≤ 2`` the level's free region is an
        interval, so :attr:`~_LevelFrontier.read` is a contiguous :class:`slice`
        (the d=2 zero-copy contiguity, reproduced byte-for-byte — the bit-id
        anchor); at ``d ≥ 3`` it is a fancy index of the simplex.  At ``d == 1``
        the free box is empty (the ``frontier_dim == 0`` degenerate point).
        """
        d = len(shape)
        # ``det`` is the SINGLE SOURCE OF TRUTH for the free/determined axis
        # partition: the determined (parity-rolled) axis is the LAST, the free
        # axes the first ``det``.  If a later stage makes the determined axis a
        # policy (non-last), this prefix slice becomes the ``b != det``
        # comprehension and every ``a < det`` test follows — change ``det`` only.
        det = d - 1
        free_bbox = tuple(shape[:det])                  # () at d=1
        box_contiguous = det <= 1                        # interval ⟺ d ≤ 2

        plan_levels: list[_LevelFrontier] = []
        for gcell, lcell in zip(levels, local_levels):
            lfree = lcell[:det]                          # (d−1)-tuple LOCAL free coords
            # READ — the level's free-coord positions (slab gather).
            if box_contiguous:
                # d ≤ 2: each free coord is a contiguous ascending range ⟹ a
                # zero-copy slice (the d=2 contiguity speedup, preserved).
                read = tuple(
                    slice(int(c[0]), int(c[-1]) + 1) for c in lfree
                )
            else:
                read = tuple(lfree)                      # d ≥ 3: fancy simplex index
            # WRITE — per axis.  Free axis: read shifted +1 on coord a.
            # Determined axis: read (the progression is the parity roll).
            write: list[tuple] = []
            for a in range(d):
                if a < det:
                    if box_contiguous:
                        w = list(read)
                        w[a] = slice(read[a].start + 1, read[a].stop + 1)
                        write.append(tuple(w))
                    else:
                        w = list(lfree)
                        w[a] = lfree[a] + 1
                        write.append(tuple(w))
                else:
                    write.append(read)
            # SEED — per axis: the inflow-edge cells (LOCAL coord_a == 0).
            seed: list[tuple | None] = []
            for a in range(d):
                mask = lcell[a] == 0
                if not mask.any():
                    seed.append(None)
                    continue
                slab_target = tuple(c[mask] for c in lfree)
                inflow_source = tuple(
                    gcell[b][mask] for b in range(d) if b != a
                )
                seed.append((slab_target, inflow_source))
            # SHED — per axis: the outflow-edge cells (LOCAL coord_a == n_a−1).
            shed: list[tuple | None] = []
            for a in range(d):
                mask = lcell[a] == shape[a] - 1
                if not mask.any():
                    shed.append(None)
                    continue
                capture_target = tuple(
                    gcell[b][mask] for b in range(d) if b != a
                )
                shed.append((mask, capture_target))
            plan_levels.append(_LevelFrontier(
                read=read, write=tuple(write),
                seed=tuple(seed), shed=tuple(shed),
            ))

        return _FrontierPlan(
            spatial_shape=tuple(shape),
            free_bbox=free_bbox,
            det=det,
            levels=tuple(plan_levels),
        )

    # ── Internal: shared per-level slice builder ────────────────────

    def _make_slice(
        self,
        cell_idx: tuple[np.ndarray, ...],
        *,
        psi_faces: tuple[np.ndarray, ...],
        Q: np.ndarray,
        sig_t: np.ndarray,
        str_axes: tuple[np.ndarray, ...],
        psi_avg_probe: np.ndarray | None = None,
    ) -> SweepCellSlice:
        """Build the per-level :class:`SweepCellSlice` for this graph's
        face-offset convention (dimension-generic).

        Shared by :meth:`apply` (solve → ``update_batch``) and
        :meth:`residual` (apply → ``residual_batch``) so the level +
        per-axis face-index wiring (``cell_idx[a] + face_in[a]`` /
        ``cell_idx[a] + face_out[a]``) is a single source of truth — the
        two walks differ only in their accumulation, not in the slice they
        build. The per-axis tuples (``cell_idx``, ``psi_faces``,
        ``str_axes``) are positional-by-axis, matching the graph's
        :attr:`face_in` / :attr:`face_out` offset order. ``psi_avg_probe``
        is ``None`` for the solve direction (``update_batch`` ignores it).
        """
        return SweepCellSlice(
            cell_idx=cell_idx,
            face_in_idx=tuple(
                cell_idx[a] + self.face_in[a] for a in range(self.ndim)
            ),
            face_out_idx=tuple(
                cell_idx[a] + self.face_out[a] for a in range(self.ndim)
            ),
            psi_faces=psi_faces,
            Q=Q, sig_t=sig_t, str_axes=str_axes,
            psi_avg_probe=psi_avg_probe,
        )

    # ── Public API ─────────────────────────────────────────────────

    def apply(
        self,
        *,
        cell_update: CellUpdateBase,
        psi_faces_octant: tuple[np.ndarray, ...],  # d face buffers — mutated in place
        Q_octant: np.ndarray,            # (N_oct or 1, ng, *spatial)
        sig_t: np.ndarray,               # (ng, *spatial)
        str_axes_octant: tuple[np.ndarray, ...],   # d streaming arrays, each (N_oct, n_a)
        weights_octant: np.ndarray,      # (N_oct,)
        angular_flux_octant: np.ndarray, # (N_oct, ng, *spatial) — written
        scalar_flux_buf: np.ndarray,     # (ng, *spatial) — accumulated into
    ) -> None:
        r"""Walk the topological levels and accumulate scalar + angular flux.

        Forward-substitute along the DAG, vectorised over
        ``(N_oct, n_diag, ng)`` per level. The ordinate axis
        (``N_oct``) is INTERNAL to every ``apply`` call: there is no
        ``for n in range(N_oct)`` anywhere in this method. Dimension-generic
        (``d = ndim``): each level is the per-axis cell-index tuple
        ``cell_idx`` (length ``d``), and the scalar-/angular-flux scatters use
        ``(…, *cell_idx)`` advanced indexing — at ``d = 2`` exactly the legacy
        ``[:, ii, jj]`` / ``[:, :, ii, jj]`` slices.

        Mutation contract
        -----------------

        * ``psi_faces_octant`` — the per-axis face buffers; outgoing face
          fluxes are scattered into them in place by
          :meth:`CellUpdateBase.update_batch`. Caller is responsible for
          seeding the incoming-face entries (BC apply happens one level up,
          in the wavefront sweep) and for scattering the post-sweep buffers
          back into the persistent BC state. The tuple is positional-by-axis
          (born from ``WavefrontFlux.face(a)`` over ``WavefrontFlux.axes`` at
          the orchestrator).
        * ``angular_flux_octant`` — written at every level's cells. Caller is
          responsible for scattering back into the global ``angular_flux``
          buffer keyed by octant indices.
        * ``scalar_flux_buf`` — accumulated into via
          ``scalar_flux_buf[(:, *cell_idx)] += weighted_sum_psi_avg``. Caller
          seeds the buffer (typically zero on the first octant of an outer
          iteration) and reads the result at the end of all octants.

        Parameters
        ----------
        cell_update :
            The closure strategy. Must override
            :meth:`update_batch` (DD does; Step / LD do not yet).
        psi_faces_octant :
            Per-octant per-axis face buffers; mutated in place.
        Q_octant :
            Per-octant per-cell volumetric source, already
            weight-normalised. Leading axis is ``1`` for
            isotropic-only sweeps or ``N_oct`` when an aniso
            component is folded in.
        sig_t :
            Per-cell per-group total cross section.
        str_axes_octant :
            Per-octant per-axis streaming coefficients ``2|μ_a|/Δa``.
        weights_octant :
            Per-octant ordinate weights for scalar-flux reduction.
        angular_flux_octant :
            Per-octant angular-flux output buffer.
        scalar_flux_buf :
            Global scalar-flux accumulator.
        """
        for cell_idx in self.levels:
            slice_args = self._make_slice(
                cell_idx,
                psi_faces=psi_faces_octant,
                Q=Q_octant, sig_t=sig_t, str_axes=str_axes_octant,
            )
            # Issue #196 PR-INDEX-5: psi_avg returned principled
            # ``(N_oct, ng, n_diag)`` (ordinate, group, anti-hyperplane).
            psi_avg = cell_update.update_batch(slice_args)  # (N_oct, ng, n_diag)
            # Scalar-flux accumulation — sum over ordinates with
            # quadrature weights. ``weights_octant`` is already
            # weight-normalised by the caller (the wavefront sweep
            # divides by Σw before invoking apply), so this is a
            # plain weighted sum. ``(:, *cell_idx)`` selects (ng, n_diag).
            scalar_flux_buf[(slice(None), *cell_idx)] += np.einsum(
                "ngd,n->gd", psi_avg, weights_octant,
            )
            # Angular-flux scatter — write all ``(N_oct, ng, n_diag)``
            # values into the per-octant principled angular-flux buffer.
            # The cell-index tuple broadcasts naturally with the leading
            # ``:`` (octant) and the second ``:`` (group).
            angular_flux_octant[(slice(None), slice(None), *cell_idx)] = psi_avg

    def residual(
        self,
        *,
        cell_update: CellUpdateBase,
        psi_faces_octant: tuple[np.ndarray, ...],  # d face buffers — mutated in place
        psi_avg_probe_octant: np.ndarray,    # (N_oct, ng, *spatial) — the probe (read)
        Q_octant: np.ndarray,                # (N_oct or 1, ng, *spatial)
        sig_t: np.ndarray,                   # (ng, *spatial)
        str_axes_octant: tuple[np.ndarray, ...],   # d streaming arrays, each (N_oct, n_a)
        residual_octant: np.ndarray,         # (N_oct, ng, *spatial) — written
    ) -> None:
        r"""Walk the topological levels accumulating the operator residual.

        The **apply-direction** companion of :meth:`apply`: where ``apply``
        forward-substitutes the *solve* (``update_batch``) and reduces to
        scalar + angular flux, ``residual`` forward-substitutes the *apply*
        (:meth:`CellUpdateBase.residual_batch`) and writes the per-cell
        operator residual :math:`(L+C)\,\overline\psi - q` into
        ``residual_octant``. The Cartesian matvec ``StreamingOperator.apply``
        (Wave O #208 O.4b) drives this so it shares the SAME wavefront DAG +
        closure as the sweep — matvec and sweep are one discretization (L21),
        no FD twin path. Dimension-generic exactly as :meth:`apply` (each
        level is the per-axis ``cell_idx`` tuple; the residual scatter uses
        ``(:, :, *cell_idx)``).

        The edge-flux reconstruction is load-bearing: each cell's incoming
        face flux is the upstream cell's outgoing flux via the diamond
        closure :math:`\psi^{\rm out}_a = 2\overline\psi^{\rm probe} -
        \psi^{\rm in}_a` (scattered into the per-axis ``psi_faces_octant``
        buffers by ``residual_batch`` in place, propagated from the seeded
        boundary inflow along the wavefront). This is exactly why the matvec
        needs the topological walk rather than a per-cell formula.

        Mutation contract
        -----------------

        * ``psi_faces_octant`` — outgoing face fluxes scattered in place
          (from the probe). Caller seeds the incoming-face entries (BC apply
          happens one level up, in the matvec) exactly as for :meth:`apply`.
        * ``residual_octant`` — written at every level's cells; caller
          scatters it back into the global bulk-residual buffer keyed by
          octant indices.

        Parameters mirror :meth:`apply`, with ``psi_avg_probe_octant``
        (the apply target) replacing the flux-reduction outputs
        (``weights_octant`` / ``scalar_flux_buf`` / ``angular_flux_octant``).
        """
        for cell_idx in self.levels:
            slice_args = self._make_slice(
                cell_idx,
                psi_faces=psi_faces_octant,
                Q=Q_octant, sig_t=sig_t, str_axes=str_axes_octant,
                psi_avg_probe=psi_avg_probe_octant,
            )
            residual_octant[(slice(None), slice(None), *cell_idx)] = (
                cell_update.residual_batch(slice_args)
            )

    # ── Phase 5b storage-B: rolling (d−1)-frontier window walks ─────────
    #
    # The PRODUCTION walks (the ``MovingFrontierWindow`` strategy). They carry
    # the interior face cochain on a :class:`_MovingFrontier` (rolling
    # (d−1)-frontier, ``O(N·ng·∏_{a<d−1} n_a)`` not the full ``O(N·ng·∏ n_a)``)
    # and ADVANCE it level by level via the mesh-time :class:`_FrontierPlan`.
    # At d=2 the frontier is a contiguous line ⟹ basic-slice views/assigns,
    # FASTER than the full-field grid-diagonal walk (~0.77×); at d≥3 it is a
    # simplex (fancy index — the memory win holds, the speed win is the deferred
    # profiling question).  The cell math is the SAME shared kernel
    # (DiamondDifference.cell_kernel_batch / residual_kernel_batch) the
    # full-field walks (apply / residual) use — so window and full-field cannot
    # drift mathematically (proven by the d=1/d=2/d=3 ``window ≡ full`` tests).
    # Domain-edge inflow is ι_*-seeded onto the frontier per level (from the
    # passed per-axis ``inflow`` tuple); domain-edge outflow is ι*-shed into the
    # per-axis ``capture`` tuple at the level that produces it, BEFORE its
    # frontier slot recycles.  The walk is dimension-agnostic: the per-axis
    # ``inflow`` / ``str_axes_octant`` / ``capture`` tuples and the d-generic
    # ``cell_idx`` scatter reproduce the legacy 2-D ``ii`` / ``jj`` access at
    # d=2 byte-for-byte (the bit-identity anchor).

    def apply_windowed(
        self,
        *,
        cell_update: CellUpdateBase,
        inflow: tuple[np.ndarray, ...],          # d-tuple — per-axis domain inflow
        Q_octant: np.ndarray,                    # (N_oct or 1, ng, *spatial)
        sig_t: np.ndarray,                       # (ng, *spatial)
        str_axes_octant: tuple[np.ndarray, ...], # d-tuple, each (N_oct, n_a)
        weights_octant: np.ndarray,              # (N_oct,)
        capture: tuple[np.ndarray, ...],         # d-tuple — per-axis domain outflow, written
        angular_flux_octant: np.ndarray | None = None,  # (N_oct, ng, *spatial) — written (angular mode)
        scalar_flux_buf: np.ndarray | None = None,       # (ng, *spatial) — accumulated (angular mode)
        moment_buf: np.ndarray | None = None,            # (L+1, 2L+1, ng, *spatial) — accumulated (moment mode)
        Y_octant: np.ndarray | None = None,              # (N_oct, L+1, 2L+1) — octant harmonics (moment mode)
    ) -> None:
        r"""Storage-B solve walk: advance a :class:`_MovingFrontier` along the
        anti-hyperplanes, sheding domain-edge outflow as it goes.

        The window-backed twin of :meth:`apply`. The interior face cochain + the
        cell math are equal to it (same shared cell kernel, same anti-hyperplane
        order); proven by the ``window ≡ full-field`` oracle test (bit-identical
        at d=2; the d=1/d=3 synthetic admission). The per-axis domain inflow is
        ι_*-seeded onto the frontier from the ``inflow`` tuple (each octant's
        incoming domain-edge trace); the outflow is ι*-shed into ``capture``.

        Output representation (Phase 5c) — exactly ONE of two modes, selected by
        which output buffer is given (dependency injection, mirroring the
        ``reflect=None`` idiom in :func:`_sweep_2d_scheduled`):

        * **angular** (``moment_buf is None``) — write the full per-ordinate
          ``angular_flux_octant[:, :, *cell_idx] = psi_avg`` per level and
          accumulate the scalar ``scalar_flux_buf``. The full-angular OUTPUT a
          reconstruction / Krylov / full-angular SI consumer needs.
        * **moment** (``moment_buf`` given) — project per anti-hyperplane
          directly into the harmonic moment tensor
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
        N_oct, ng = inflow[0].shape[0], inflow[0].shape[1]
        d = self.ndim
        frontier = _MovingFrontier(N_oct, ng, self.window_plan)
        for k, cell_idx in enumerate(self.levels):
            cur, prev = k % 2, (k - 1) % 2
            frontier.seed(prev, k, inflow)
            in_faces = frontier.incoming(prev, k)
            psi_avg, out_faces = cell_update.cell_kernel_batch(
                psi_in=in_faces,
                s_axes=tuple(
                    str_axes_octant[a][:, cell_idx[a]][:, None, :]
                    for a in range(d)
                ),
                sigt_cells=sig_t[(slice(None), *cell_idx)],
                Q_cells=Q_octant[(slice(None), slice(None), *cell_idx)],
            )
            frontier.emit(cur, k, out_faces)
            frontier.shed(k, out_faces, capture)
            if moment_buf is None:
                scalar_flux_buf[(slice(None), *cell_idx)] += np.einsum(
                    "ngd,n->gd", psi_avg, weights_octant,
                )
                angular_flux_octant[
                    (slice(None), slice(None), *cell_idx)
                ] = psi_avg
            else:
                moment_buf[
                    (slice(None), slice(None), slice(None), *cell_idx)
                ] += np.einsum("nlm,ngd,n->lmgd", Y_octant, psi_avg, weights_octant)

    def residual_windowed(
        self,
        *,
        cell_update: CellUpdateBase,
        inflow: tuple[np.ndarray, ...],          # d-tuple — per-axis domain inflow
        psi_avg_probe_octant: np.ndarray,        # (N_oct, ng, *spatial) — the probe
        Q_octant: np.ndarray,                    # (N_oct or 1, ng, *spatial)
        sig_t: np.ndarray,                       # (ng, *spatial)
        str_axes_octant: tuple[np.ndarray, ...], # d-tuple, each (N_oct, n_a)
        residual_octant: np.ndarray,             # (N_oct, ng, *spatial) — written
        capture: tuple[np.ndarray, ...],         # d-tuple — per-axis domain outflow, written
    ) -> None:
        r"""Storage-B apply (matvec) walk: the :class:`_MovingFrontier`-backed
        twin of :meth:`residual`. Reconstructs edge fluxes from the PROBE along
        the rolling frontier (``psi_out = 2·psi_bar − psi_in``) and writes the
        per-cell operator residual; sheds the domain-edge outflow. The matvec
        residual output stays full ``(N_oct, ng, *spatial)`` (Krylov consumes
        the full bulk residual) — only the interior FACE cochain is windowed.
        """
        N_oct, ng = inflow[0].shape[0], inflow[0].shape[1]
        d = self.ndim
        frontier = _MovingFrontier(N_oct, ng, self.window_plan)
        for k, cell_idx in enumerate(self.levels):
            cur, prev = k % 2, (k - 1) % 2
            frontier.seed(prev, k, inflow)
            in_faces = frontier.incoming(prev, k)
            res, out_faces = cell_update.residual_kernel_batch(
                psi_bar=psi_avg_probe_octant[(slice(None), slice(None), *cell_idx)],
                psi_in=in_faces,
                s_axes=tuple(
                    str_axes_octant[a][:, cell_idx[a]][:, None, :]
                    for a in range(d)
                ),
                sigt_cells=sig_t[(slice(None), *cell_idx)],
                Q_cells=Q_octant[(slice(None), slice(None), *cell_idx)],
            )
            frontier.emit(cur, k, out_faces)
            frontier.shed(k, out_faces, capture)
            residual_octant[(slice(None), slice(None), *cell_idx)] = res


__all__ = [
    "OctantLabel",
    "SweepDependencyGraph",
]
