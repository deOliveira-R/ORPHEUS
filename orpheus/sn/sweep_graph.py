r"""Per-octant causal cell DAG for the Cartesian wavefront sweep (d-generic).

This module ships the **§15A.2 "upwind trace complex / causal transport
DAG / direction sweep ordering"** primitive (Grand Report v3 lines
2137-2171) for Cartesian SN, lifting the per-call ``_diag_cache``
precompute that historically lived inside the 2-D sweep body (now the
Jacobi entry :func:`orpheus.sn.loss_representation._sweep_jacobi`) to
mesh-time work, and replacing the per-ordinate Python loop over
ordinates with a per-octant batched ``apply``.

The graph is a **derived object** — it depends only on cell topology
(``Mesh2D``) and the octant sign convention. It does NOT depend on
fluxes, sources, or iteration state, so it is built ONCE per
``(shape, octant)`` pair in the :meth:`SweepDependencyGraph.for_shape`
cache (S6.4(c); historically at mesh construction, Wave 2 / C2.4) and
reused across every source iteration / Krylov matvec / outer iteration.

Architectural framing (Cardinal Rule 2)
=======================================

Per the project memory note ``project_moc_structure.md`` and
:doc:`../sn/spatial/scheme`, this DAG abstraction is **SN-specific
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
cell-kernel call (since S6.4(e): the level operation dispatching
``cell_kernel_batch`` / ``residual_kernel_batch``) — no Python loop
over ordinates, no Python loop over cells along a diagonal, no Python
loop over groups. The outer sweep iterates octants (4 in 2-D —
structural) and topological levels per octant (structural — sweep-DAG
data dependency). The ordinate axis is internal to every kernel /
``einsum`` invocation.

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
  this module's design + dispatch boundary with ``DiscretizationScheme``.
* :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.cell_kernel_batch`
  / :meth:`~orpheus.sn.spatial.scheme.DiscretizationSchemeBase.residual_kernel_batch`
  — the closure-pluggable kernel pair the level operations dispatch
  (since S6.4(e); historically the ``SweepCellSlice``-packeted
  ``update_batch``, Wave 2 / C2.2).

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

import itertools
from dataclasses import dataclass
from functools import lru_cache
from types import MappingProxyType

import numpy as np

from orpheus.numerics.moment_layout import is_moment_valued_by_rank
from orpheus.sn.spatial._ubld import face_moment_tail, octant_moment_frame_signs
from orpheus.sn.spatial.scheme import DiscretizationSchemeBase


def _reframe(
    arr: np.ndarray,
    frame_signs: "np.ndarray | None",
    *,
    is_moment_valued: bool,
) -> np.ndarray:
    r"""Apply the sweep⇄global moment-frame involution to a moment array.

    Reframe iff this is a moment buffer at a non-trivial closure.

    ``frame_signs`` is the ``2^d``-length involution from
    :func:`~orpheus.sn.spatial._ubld.octant_moment_frame_signs` (#240 D5b-S3),
    and is ``None`` for a single-moment closure (DD/Step) — the
    already-typed ``frame_signs is None ⟺ DD/Step`` gate.  ``is_moment_valued``
    states, from the array's TYPED ORIGIN, whether ``arr`` actually carries the
    trailing ``2^d`` moment axis (a genuine iterate / moment source / residual)
    or is a flat scalar buffer (a matvec-zero / flat external source, whose
    only populated moment is the sign-invariant average).  Only a moment array
    at a multi-moment closure gets re-signed; everything else passes through
    untouched, so DD/Step stay byte-identical and a flat source is never
    broadcast into a spurious moment axis.  The map is its own inverse
    (global→sweep on input == sweep→global on output).

    #246: the old guard keyed the moment-axis question on a COINCIDENTAL
    trailing-axis length (``arr.shape[-1] != frame_signs.shape[0]``) — at d=2
    a non-moment array whose trailing axis is coincidentally ``2^d == 4`` would
    mis-fire the sign-flip.  ``is_moment_valued`` keys the question on intent
    (the caller knows the array's typed origin), not on a coincidental size.
    """
    if frame_signs is None or not is_moment_valued:
        return arr
    return arr * frame_signs


# ═══════════════════════════════════════════════════════════════════════
# OctantLabel — in-plane octant signature (one sign per mesh axis)
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True, slots=True)
class OctantLabel:
    r"""Octant signature for the dimension-generic wavefront sweep.

    Carries one direction sign per spatial axis,
    ``signs[axis] ∈ {-1, 0, +1}``, so a single type labels a 1-D
    (``(±1,)``), 2-D (``(±1, ±1)``), or 3-D (``(±1, ±1, ±1)``) octant.
    Ordinate-label signs beyond the mesh's spatial ``ndim`` (e.g. the
    ``sign_z`` of an ``S²`` ordinate over a 2-D mesh) are projected out
    by the schedule (:func:`orpheus.sn.sweep_schedule._octant_sweep` —
    the SOLE in-plane projection site; the in-plane sweep is invariant
    under the out-of-plane signs); multiple ordinates that project to
    the same in-plane ``signs`` share a single
    :class:`SweepDependencyGraph`.

    A label whose signs are *all* zero denotes the degenerate
    no-streaming set of ordinates (e.g. the pure-:math:`z` ordinates in
    2-D); the wavefront sweep handles those via a ``Q / Σ_t``
    short-circuit and builds no graph for it (:attr:`streams` is
    ``False``).

    Frozen + slotted: hashable for use as a mapping key (the per-shape
    family :meth:`SweepDependencyGraph.for_shape` is keyed by it).
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


# ═══════════════════════════════════════════════════════════════════════
# _FrontierPlan / _MovingFrontier — the rolling (d−1)-frontier face cochain
# ═══════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class _LevelFrontier:
    r"""Mesh-time slab addressing for ONE anti-hyperplane level (storage-B).

    All :math:`d`-dependent index arithmetic for one level of the rolling
    :math:`(d{-}1)`-frontier window, precomputed so the walk
    (:meth:`SweepDependencyGraph.walk_windowed` with either level
    operation) is dimension-agnostic
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
        Domain-edge inflow ops, ONE axis-tagged record
        ``(axis, slab_target, inflow_source)`` per axis whose inflow edge is
        present on this level (NOT a d-tuple padded with ``None`` — most
        interior levels carry an empty tuple): scatter
        ``inflow[axis][*inflow_source]`` into ``win[axis][prev][*slab_target]``
        before the gather.
    shed :
        Domain-edge outflow ops, same axis-tagged-record form
        ``(axis, out_mask, capture_target)`` per present outflow edge: gather
        the kernel's outgoing ``out_faces[axis][..., out_mask]`` into
        ``capture[axis][*capture_target]`` after the kernel.
    """

    read: tuple
    write: tuple
    seed: tuple
    shed: tuple


@dataclass(frozen=True)
class _FrontierPlan:
    r"""The whole-sweep rolling :math:`(d{-}1)`-frontier window plan (storage-B).

    A face produced at level ``k`` is consumed at ``k+1``, so the interior face
    cochain (the cochain :math:`C^1_{\rm int}`; its full-field realization is
    the per-octant ``_octant_face_cochain`` buffers — historically the typed
    ``WavefrontFlux``, retired S6.4(f)) need only be held on the two active
    levels — a
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
    measured-cost question deferred to profiling (no 3-axis MESH yet — the
    angular quadrature is already 3-cosine with all 8 sign-octants; the
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

    __slots__ = ("_win", "_plan", "_face_moments")

    def __init__(
        self, N_oct: int, ng: int, plan: _FrontierPlan, n_face_moments: int = 1,
    ) -> None:
        free = plan.free_bbox
        det = plan.det
        # The interior face cochain carries a trailing ``2^{d-1}``-moment axis
        # for a multi-moment closure (LD's bilinear face — #240 D5b); the
        # slopeless cell-average closures (DD/Step) leave ``n_face_moments == 1``
        # and the trailing axis is ABSENT (no length-1 axis appended), keeping
        # the rank-r face buffers byte-identical.  The seed/incoming/emit/shed
        # selectors index the leading (N_oct, ng, parity, *slab) axes only, so a
        # trailing moment axis rides along untouched.
        tail = face_moment_tail(n_face_moments)
        win = []
        for a, n_a in enumerate(plan.spatial_shape):
            shp = list(free)
            if a < det:                          # FREE axis: +1 ghost on coord a
                shp[a] = n_a + 1
            win.append(np.zeros((N_oct, ng, 2, *shp, *tail)))
        self._win = tuple(win)
        self._plan = plan
        self._face_moments = n_face_moments

    def seed(self, prev: int, lvl: int, inflow: tuple[np.ndarray, ...]) -> None:
        r""":math:`\iota_*` — inject the domain inflow into the slab at this
        level's edge cells (so the subsequent :meth:`incoming` gather reads
        them).  Iterates only the axes that HAVE an inflow edge on this level."""
        for axis, slab_target, inflow_source in self._plan.levels[lvl].seed:
            self._win[axis][(slice(None), slice(None), prev) + slab_target] = (
                inflow[axis][(slice(None), slice(None)) + inflow_source]
            )

    def incoming(self, prev: int, lvl: int) -> tuple[np.ndarray, ...]:
        """Gather the level's incoming face flux per axis (a zero-copy VIEW at
        d≤2 via the contiguous slice; a fancy-index copy at d≥3)."""
        gate = (slice(None), slice(None), prev) + self._plan.levels[lvl].read
        if self._plan.is_point:                  # d=1: restore the n_diag axis
            return tuple(win_a[gate][:, :, None] for win_a in self._win)
        return tuple(win_a[gate] for win_a in self._win)

    def emit(
        self, cur: int, lvl: int, out_faces: tuple[np.ndarray, ...],
    ) -> None:
        """Advance the frontier: scatter each axis's outgoing face flux into the
        downstream slot (free axis ``f + e_a``; determined axis ``f``).

        The three per-axis collections — the slabs, the per-axis write
        selectors, and the kernel's outgoing faces — are advanced together
        (``zip``), the determined-vs-free distinction already baked into each
        selector by the plan.
        """
        point = self._plan.is_point
        for win_a, write_a, out_a in zip(
            self._win, self._plan.levels[lvl].write, out_faces,
        ):
            win_a[(slice(None), slice(None), cur) + write_a] = (
                out_a[:, :, 0] if point else out_a
            )

    def shed(
        self, lvl: int, out_faces: tuple[np.ndarray, ...],
        capture: tuple[np.ndarray, ...],
    ) -> None:
        r""":math:`\iota^*` — capture the domain outflow (the high-edge cells'
        outgoing faces) into the per-axis ``capture`` arrays.  Iterates only the
        axes that HAVE an outflow edge on this level.  At the d=1 degenerate
        point the single edge cell's n_diag axis is squeezed (the domain
        outflow has no perpendicular coordinate)."""
        point = self._plan.is_point
        for axis, out_mask, capture_target in self._plan.levels[lvl].shed:
            shed_faces = (
                out_faces[axis][:, :, 0] if point
                else out_faces[axis][:, :, out_mask]
            )
            capture[axis][(slice(None), slice(None)) + capture_target] = shed_faces


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

    The walks do NOT inline the per-cell update math — the level
    operation dispatches the strategy's kernel pair
    (``cell_kernel_batch`` / ``residual_kernel_batch``), so the
    WDD / LD / EC / Step closure is pluggable. This aligns with the
    strategy contract in :mod:`orpheus.sn.spatial.scheme`.

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
    update is delegated to the strategy's batched kernel pair via the
    level operation, so the ordinate axis (``N_oct``) is internal to
    every kernel call: there is no ``for n in range(N_oct)`` anywhere
    in this module.
    """

    label: OctantLabel
    levels: tuple[tuple[np.ndarray, ...], ...]
    face_in: tuple[int, ...]
    face_out: tuple[int, ...]
    spatial_shape: tuple[int, ...]
    # ── Phase 5b storage-B: mesh-time rolling (d−1)-frontier window plan ──
    # The production-optimization addressing for the ``walk_windowed``
    # walk (the ``MovingFrontierWindow`` strategy): the
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
    def for_shape(
        cls, spatial_shape: tuple[int, ...],
    ) -> "MappingProxyType[OctantLabel, SweepDependencyGraph]":
        r"""The per-octant DAG family for a Cartesian grid shape — CACHED.

        One :class:`SweepDependencyGraph` per streaming octant: the
        :math:`2^d` sign signatures ``itertools.product((-1, +1), repeat=d)``
        over the ``d = len(spatial_shape)`` spatial axes.  A 1-D chain gets
        2 graphs (``±x``); a 2-D mesh the 4 in-plane octants.  (Pure-z
        ordinates with an all-zero in-plane sign are handled by the walk's
        ``Q/Σ_t`` short-circuit and have no entry here.)

        The graphs depend ONLY on cell topology + octant sign convention —
        independent of fluxes / sources / iteration state / cross sections —
        so they are cached **per shape** and shared by every consumer on a
        same-shape mesh (S6.4(c): the DAG is OWNED by the
        ``_DAGWavefront`` representation family through this accessor, not
        by the mesh — a mesh is pure geometry; DAG-free representations
        never mention this substrate).

        At ``d = 2`` ``itertools.product`` yields the octants in the same
        order as the historical nested ``sx, sy`` loop — ``(-1,-1), (-1,+1),
        (+1,-1), (+1,+1)`` — so the dict is built bit-for-bit as the
        retired mesh-construction build did.

        The returned mapping is a read-only :class:`~types.MappingProxyType`
        — it is shared across callers, so mutation is unrepresentable (the
        small LRU bound only re-pays the construction cost when many
        distinct shapes interleave).
        """
        return _graphs_for_shape(tuple(int(n) for n in spatial_shape))

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
        :meth:`walk_windowed` walks are
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
            # SEED — only the axes whose inflow edge (LOCAL coord_a == 0) is
            # PRESENT on this level (an axis-tagged record per real edge, NOT a
            # d-tuple padded with ``None``): the walk iterates the edges, never
            # the axes — most interior levels touch no edge at all.
            seed: list[tuple] = []
            for a in range(d):
                mask = lcell[a] == 0
                if not mask.any():
                    continue
                slab_target = tuple(c[mask] for c in lfree)
                inflow_source = tuple(
                    gcell[b][mask] for b in range(d) if b != a
                )
                seed.append((a, slab_target, inflow_source))
            # SHED — only the axes whose outflow edge (LOCAL coord_a == n_a−1) is
            # present on this level (same axis-tagged-record form as SEED).
            shed: list[tuple] = []
            for a in range(d):
                mask = lcell[a] == shape[a] - 1
                if not mask.any():
                    continue
                capture_target = tuple(
                    gcell[b][mask] for b in range(d) if b != a
                )
                shed.append((a, mask, capture_target))
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

    # ── The level walks — TWO storage policies × ONE direction object ──
    #
    # S6.4(e): the four direction×storage product methods (``apply`` /
    # ``residual`` full-field + ``apply_windowed`` / ``residual_windowed``)
    # COLLAPSED into two storage walks parameterized by a LEVEL OPERATION
    # object (``_CellSolve`` | ``_CellResidual`` — the solve/apply direction
    # fork, bottoming out in the diamond kernel pair).  The walk owns the
    # level loop + storage (frontier or full cochain) + the per-level
    # operand extraction; the level op owns the cell algebra dispatch + the
    # per-level emit.  Direction is an OBJECT, never a boolean flag.

    def walk_full(
        self,
        *,
        level_op: "_CellSolve | _CellResidual",
        psi_faces_octant: tuple[np.ndarray, ...],  # d face buffers — mutated in place
        Q_octant: np.ndarray,                      # (N_oct or 1, ng, *spatial)
        sig_t: np.ndarray,                         # (ng, *spatial)
        str_axes_octant: tuple[np.ndarray, ...],   # d arrays, each (N_oct, n_a)
    ) -> None:
        r"""Full-cochain level walk — the verification-oracle storage policy.

        Walks the topological levels carrying the COMPLETE per-axis interior
        face cochain (every face retained — the fuller view): per level,
        gather the incoming faces at ``cell_idx[a] + face_in[a]``, run the
        level operation, scatter the outgoing faces at
        ``cell_idx[a] + face_out[a]``.  The per-axis gather/scatter is the
        d-generic advanced-index loop (at d = 2 byte-identical to the legacy
        ``psi_x[:, :, face, jj]`` / ``psi_y[:, :, ii, face]`` slices).

        Mutation contract: ``psi_faces_octant`` carries the seeded in-edges on
        entry and the walked cochain (out-edges included) on exit; the level
        op writes its own output buffers.
        """
        d = self.ndim
        for cell_idx in self.levels:
            psi_in = tuple(
                psi_faces_octant[a][
                    _cell_face_selector(cell_idx, a, cell_idx[a] + self.face_in[a])
                ]
                for a in range(d)
            )
            psi_out = level_op.cell(
                cell_idx,
                psi_in=psi_in,
                s_axes=tuple(
                    str_axes_octant[a][:, cell_idx[a]][:, None, :] for a in range(d)
                ),
                reaction_xs=sig_t[(slice(None), *cell_idx)],
                Q_cells=Q_octant[(slice(None), slice(None), *cell_idx)],
            )
            for a in range(d):
                psi_faces_octant[a][
                    _cell_face_selector(cell_idx, a, cell_idx[a] + self.face_out[a])
                ] = psi_out[a]

    def walk_windowed(
        self,
        *,
        level_op: "_CellSolve | _CellResidual",
        inflow: tuple[np.ndarray, ...],            # d-tuple — per-axis domain inflow
        Q_octant: np.ndarray,                      # (N_oct or 1, ng, *spatial)
        sig_t: np.ndarray,                         # (ng, *spatial)
        str_axes_octant: tuple[np.ndarray, ...],   # d-tuple, each (N_oct, n_a)
        capture: tuple[np.ndarray, ...],           # d-tuple — per-axis domain outflow, written
        n_face_moments: int = 1,                   # trailing 2^{d-1} face-moment axis (LD multi-D)
    ) -> None:
        r"""Rolling-frontier level walk — the storage-B PRODUCTION policy.

        Advances a :class:`_MovingFrontier` (the rolling :math:`(d{-}1)`-
        frontier, ``O(N·ng·∏_{a<d−1} n_a)`` instead of the full cochain)
        along the anti-hyperplanes: per level, ι_*-seed the frontier's
        domain-edge slots from ``inflow``, gather the incoming faces, run the
        level operation, advance the frontier with the outgoing faces, and
        ι*-shed the domain-edge outflow into ``capture`` before its slot
        recycles.  Pinned bit-identical to :meth:`walk_full` by the
        ``window ≡ full`` oracles (d=1/d=2/d=3) — same cell math, same level
        order, different storage.

        ``n_face_moments`` is the per-face transverse moment count
        :math:`(\text{per\_axis})^{d-1}` (DD/Step: 1; LD-2D: 2) — the frontier
        slabs carry a trailing moment axis when ``> 1`` (#240 D5b), absent
        otherwise (DD byte-identical).
        """
        N_oct, ng = inflow[0].shape[0], inflow[0].shape[1]
        frontier = _MovingFrontier(N_oct, ng, self.window_plan, n_face_moments)
        for k, cell_idx in enumerate(self.levels):
            cur, prev = k % 2, (k - 1) % 2
            frontier.seed(prev, k, inflow)
            psi_in = frontier.incoming(prev, k)
            psi_out = level_op.cell(
                cell_idx,
                psi_in=psi_in,
                # s_axes — the level's streaming coefficient per axis: the
                # per-axis streaming array and the level's per-axis cell index
                # advanced TOGETHER (the latent zip, not range(d) + indexing).
                s_axes=tuple(
                    s[:, c][:, None, :]
                    for s, c in zip(str_axes_octant, cell_idx)
                ),
                reaction_xs=sig_t[(slice(None), *cell_idx)],
                Q_cells=Q_octant[(slice(None), slice(None), *cell_idx)],
            )
            frontier.emit(cur, k, psi_out)
            frontier.shed(k, psi_out, capture)


def _cell_face_selector(
    cell_idx: tuple[np.ndarray, ...], axis: int, face_idx: np.ndarray,
) -> tuple:
    r"""Advanced-index tuple selecting axis ``axis``'s face for one level.

    Starts from the cell-centred selector ``(:, :, *cell_idx)`` (the leading
    ``(:, :)`` skips the ``(N_oct, ng)`` prefix) and replaces axis ``axis``'s
    in-cell index with the face index ``face_idx`` at lattice position
    ``2 + axis``.  At ``d = 2`` this is ``(:, :, face, jj)`` for axis 0 and
    ``(:, :, ii, face)`` for axis 1 — byte-identical to the legacy
    hand-written face slices.  (S6.4(e): relocated from the retired
    ``DiamondDifference`` storage adapters — face addressing is the WALK's
    storage concern, not the discretization's.)
    """
    sel = [slice(None), slice(None), *cell_idx]
    sel[2 + axis] = face_idx
    return tuple(sel)


# ═══════════════════════════════════════════════════════════════════════
# The level operations — the direction fork, as OBJECTS (S6.4(e))
# ═══════════════════════════════════════════════════════════════════════
#
# ONE of these is constructed per octant walk; the storage walks call
# ``level_op.cell(...)`` per level.  The solve/apply fork bottoms out in the
# discretization's pure kernel pair
# (:meth:`DiamondDifference.cell_kernel_batch` /
# :meth:`~DiamondDifference.residual_kernel_batch`); the level op adds only
# the direction's per-level EMIT.  The emit expressions + their order are
# bit-identity-load-bearing — relocated VERBATIM from the four retired walk
# methods.


@dataclass(frozen=True)
class _CellSolve:
    r"""SOLVE-direction level operation: solve the WDD balance for
    :math:`\bar\psi`, emit angular or harmonic-moment output.

    Output mode is the Phase-5c DI (exactly one of angular / moment is
    wired — the caller passes the mode's buffers, the established
    output-DI idiom):

    * **angular** — write ``angular_flux_octant[:, :, *cell_idx]`` and
      accumulate ``scalar_flux_buf[(:, *cell_idx)] += Σ_n w_n ψ_n``;
    * **moment** — accumulate the harmonic tensor
      :math:`\phi_\ell^m \mathrel{+}= \sum_n w_n Y_\ell^m \psi_n`
      per level (the full angular field never materializes).
    """

    scheme: DiscretizationSchemeBase
    weights_octant: np.ndarray                    # (N_oct,)
    angular_flux_octant: "np.ndarray | None" = None
    scalar_flux_buf: "np.ndarray | None" = None
    moment_buf: "np.ndarray | None" = None
    Y_octant: "np.ndarray | None" = None
    moment_frame_signs: "np.ndarray | None" = None  # (2^d,) sweep⇄global slope frame

    def __post_init__(self) -> None:
        # Exactly ONE output mode must be wired (mirrors _SweepEmit's guard
        # one level up) — a half-wired mode would otherwise crash deep inside
        # the walk, far from the construction site.
        angular = (
            self.angular_flux_octant is not None
            and self.scalar_flux_buf is not None
        )
        moment = (self.moment_buf is not None) and (self.Y_octant is not None)
        if angular == moment:
            raise ValueError(
                "_CellSolve: exactly ONE output mode must be wired — either "
                "(angular_flux_octant AND scalar_flux_buf) or "
                "(moment_buf AND Y_octant)."
            )

    def cell(
        self,
        cell_idx: tuple[np.ndarray, ...],
        *,
        psi_in: tuple[np.ndarray, ...],
        s_axes: tuple[np.ndarray, ...],
        reaction_xs: np.ndarray,
        Q_cells: np.ndarray,
    ) -> tuple[np.ndarray, ...]:
        # The cell kernel consumes/produces the moment vector in the per-ordinate
        # SWEEP frame.  The iterate (φ̂) + its scattering source ``Σ_s·φ̂`` live in
        # the GLOBAL frame, so a GENUINE moment source ``Q_cells`` is mapped
        # global→sweep on input and the emitted ``psi_avg`` sweep→global on output
        # (#240 D5b-S3 — the diffusion-limit root cause: backward-ordinate slopes
        # are sign-flipped between the two frames; summing them un-corrected
        # cancels φ̂).  ``_reframe`` applies the ``2^d`` involution
        # (``octant_moment_frame_signs``) ONLY to a genuine moment array (trailing
        # axis == 2^d); a FLAT scalar source (matvec zero / flat external)
        # populates only the average moment (sign +1) and is frame-invariant, so
        # it passes through untouched.  DD/Step (``None``) → byte-identical.  The
        # OUTGOING FACE (``psi_out``) stays in the sweep frame — it propagates
        # along the wavefront, never crosses into the global-frame iterate.
        moment_valued = self.scheme.is_multi_moment
        # The cell SOLVE source can arrive EITHER moment-lifted (production
        # ``solve_sn_fixed_source`` runs ``_lift_external_source_to_moments``,
        # so ``Q_cells`` is ``(N_oct, ng, n_diag, 2^d)`` — its trailing axis IS
        # the moment axis, slopes carry the scattering source ``Σ_s·φ̂`` that
        # needs global→sweep re-signing) OR flat (the lower-level
        # ``FullFieldWavefront``/``MovingFrontierWindow`` ``sweep`` API passes a
        # flat ``(N_oct, ng, n_diag)`` source that ``_ubld_system`` lifts onto
        # slot 0 itself).  A flat source has only the sign-invariant average
        # moment, so it must NOT be reframed.  Discriminate by RANK against
        # ``reaction_xs`` via the shared S4-safe
        # :func:`~orpheus.numerics.moment_layout.is_moment_valued_by_rank`
        # (single-sourced with ``_moment_broadcast_sigma``; a coincidental
        # ``n_diag == 2^d`` would mis-fire a trailing-length probe, the rank
        # cannot).  At DD/Step ``frame_signs`` is None → no-op regardless.
        source_is_moment = is_moment_valued_by_rank(Q_cells, reaction_xs)
        psi_avg, psi_out = self.scheme.cell_kernel_batch(
            psi_in=psi_in, s_axes=s_axes, reaction_xs=reaction_xs,
            Q_cells=_reframe(
                Q_cells, self.moment_frame_signs, is_moment_valued=source_is_moment,
            ),
        )
        # ``psi_avg`` is the cell-solve moment output — a genuine ``2^d`` moment
        # vector at a multi-moment closure (LD), a scalar at DD/Step.
        psi_avg = _reframe(
            psi_avg, self.moment_frame_signs, is_moment_valued=moment_valued,
        )
        # Multi-moment closures (LD's bilinear UBLD in EVERY d ≥ 1) return a
        # trailing ``2^d``-moment ``psi_avg`` and the spatial-moment iterate φ̂
        # accumulates ALL of it (#240 D5b-S3 — Increment C carries φ̂ between
        # sweeps so the scattering-slope source ``Σ_s·φ̂`` couples globally).
        # The emit einsums are spatial-moment-axis-AGNOSTIC (``ng...`` /
        # ``nlm,ng...``) — the same convention :meth:`_SweepEmit.pure_z` uses —
        # so the trailing ``2^d`` moment axis (the cell-level ``d`` plus the
        # spatial-moment ``p``) rides through with no per-axis branch.
        # DiamondDifference at ``per_axis == 1`` returns a rank-3 (scalar)
        # ``psi_avg`` → the trailing pack is just the cell axis ``d``, the
        # original byte-identical reduction (the negative control).
        cell = (slice(None), *cell_idx)
        cell_g = (slice(None), *cell)
        if self.moment_buf is None:
            self.scalar_flux_buf[cell] += np.einsum(
                "ng...,n->g...", psi_avg, self.weights_octant,
            )
            self.angular_flux_octant[cell_g] = psi_avg
        else:
            self.moment_buf[(slice(None), *cell_g)] += np.einsum(
                "nlm,ng...,n->lmg...", self.Y_octant, psi_avg, self.weights_octant,
            )
        return psi_out


@dataclass(frozen=True)
class _CellResidual:
    r"""APPLY-direction level operation: evaluate the operator residual
    :math:`r = (\Sigma_t + \sum_a s_a)\,\bar\psi - (Q + \sum_a s_a\,
    \psi^{\rm in}_a)` at the GIVEN probe, reconstructing the outgoing faces
    from it (``out = 2ψ̄ − in``) so the edge fluxes propagate along the
    wavefront exactly as the sweep's do — matvec ≡ sweep, ONE discretization
    (L21).
    """

    scheme: DiscretizationSchemeBase
    psi_avg_probe_octant: np.ndarray              # (N_oct, ng, *spatial[, 2^d]) — read
    residual_octant: np.ndarray                   # (N_oct, ng, *spatial[, 2^d]) — written
    moment_frame_signs: "np.ndarray | None" = None  # (2^d,) sweep⇄global slope frame

    def cell(
        self,
        cell_idx: tuple[np.ndarray, ...],
        *,
        psi_in: tuple[np.ndarray, ...],
        s_axes: tuple[np.ndarray, ...],
        reaction_xs: np.ndarray,
        Q_cells: np.ndarray,
    ) -> tuple[np.ndarray, ...]:
        # The unified moment matvec (#240 D5b-S3): the probe carries the full
        # ``2^d``-moment spatial iterate (a trailing axis on the bulk field at a
        # multi-moment closure), so the per-cell slice feeds the bilinear UBLD
        # residual its full moment vector — the d≥2 raise is RETIRED.  At a
        # single-moment closure (DD/Step) there is no trailing axis and the slice
        # yields the scalar probe byte-identically (the negative control).
        #
        # The probe (the iterate) is GLOBAL-frame; the residual kernel works in
        # the per-ordinate SWEEP frame (the same frame the cell SOLVE uses), so
        # the probe + source map global→sweep on input and the residual maps
        # sweep→global on output — the matvec twin of the SOLVE's frame map
        # (#240 D5b-S3 root cause; ``moment_frame_signs`` is the same ``2^d``
        # involution).  DD/Step (``None``) → byte-identical.  The OUTGOING FACE
        # stays sweep-frame (it propagates along the wavefront).
        cell_g = (slice(None), slice(None), *cell_idx)
        signs = self.moment_frame_signs
        moment_valued = self.scheme.is_multi_moment
        residual, psi_out = self.scheme.residual_kernel_batch(
            # The probe (the iterate) carries the full ``2^d`` moment axis at a
            # multi-moment closure (LD); a scalar probe at DD/Step.
            psi_bar=_reframe(
                self.psi_avg_probe_octant[cell_g], signs,
                is_moment_valued=moment_valued,
            ),
            psi_in=psi_in, s_axes=s_axes,
            reaction_xs=reaction_xs,
            # The matvec source is an all-zero buffer (the operator action
            # ``(L+C)ψ̄`` carries no volumetric source) — identically zero, so
            # the involution is a no-op regardless; pass ``False`` honestly.
            Q_cells=_reframe(Q_cells, signs, is_moment_valued=False),
        )
        # The residual is a genuine ``2^d`` moment vector at a multi-moment
        # closure (LD), a scalar at DD/Step.
        self.residual_octant[cell_g] = _reframe(
            residual, signs, is_moment_valued=moment_valued,
        )
        return psi_out



@lru_cache(maxsize=8)
def _graphs_for_shape(
    spatial_shape: tuple[int, ...],
) -> "MappingProxyType[OctantLabel, SweepDependencyGraph]":
    """The cached body of :meth:`SweepDependencyGraph.for_shape` (module-level
    so ``lru_cache`` keys on the shape tuple alone, not a class argument).

    Returns a :class:`~types.MappingProxyType` — the family is SHARED across
    every same-shape consumer, so mutation is made unrepresentable rather
    than merely discouraged (the proxy itself is the cached object, keeping
    same-shape identity stable)."""
    return MappingProxyType({
        OctantLabel(signs): SweepDependencyGraph.from_cartesian(
            spatial_shape, label=OctantLabel(signs),
        )
        for signs in itertools.product((-1, +1), repeat=len(spatial_shape))
    })


__all__ = [
    "OctantLabel",
    "SweepDependencyGraph",
]
