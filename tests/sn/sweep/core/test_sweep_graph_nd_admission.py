r"""C3 synthetic-N-D shape-admission pins for ``SweepDependencyGraph``.

THIS FILE IS A TEST-ARCHITECT STUB — UNCOMMITTED, FOR MAIN-AGENT REVIEW.

It pins the STRUCTURAL contract that the dimension-generic
``SweepDependencyGraph.from_cartesian(shape)`` admits ``d = 1`` (a REAL
production compute path under the spine+optimization layering — see
below) AND ``d = 3`` shapes WITHOUT a 3-D quadrature or a converged
3-D solve (C3 invariant #3: "synthetic-3-D SHAPE ADMISSION only").
The assertions are shape/count/lattice-level facts on the upwind DAG +
the diamond cell kernel's streaming-sum reduction — NOT a flux value.

Spine + pluggable-optimizations (the EXPANDED C3 architecture)
--------------------------------------------------------------
The wavefront DAG full-field walk (``from_cartesian(shape)`` for
d = 1, 2, 3 + the d-generic ``DiamondDifference.cell_kernel_batch``
summing ``denom = σ_t + Σ_axis s_axis``) is the dimension-agnostic
SPINE and the verification ORACLE — ONE code path for every d. Two
optimizations plug in WITHOUT changing the spine's value:

* d = 1: the cumprod parallel-prefix scan (``ordinate_scan`` /
  ``_sweep_1d_unified``) — default-selected for d=1 production. The
  spine instantiated at ``from_cartesian((nx,))`` MUST compute the
  same 1-D recurrence (proven by
  ``tests/sn/sweep/core/test_wavefront_cumprod_equivalence.py`` — the
  separate oracle file, since it lives near the 1-D sweep math).
* d ≥ 2: the ``_MovingFrontier`` storage-B window — pinned ≡ full-field
  by ``tests/sn/sweep/cartesian_2d/test_2d_full_field_oracle.py``.

So d = 1 here is NOT a synthetic-shape stub: it is a real compute path
whose DAG structure (NEW-3 below: ``2^1 = 2`` octants, ``nx`` levels,
one cell per level, single-axis ``denom = σ_t + s_x``) the B6 pins
nail down. The d = 3 pins (B2–B5) remain pure shape admission.

These are ``@pytest.mark.foundation`` software invariants (no theory
``:label:``) per ``vv-principles`` L&V level taxonomy: the d=3 DAG is
a data-structure contract, not an equation claim. They carry NO
``verifies(...)``.

``-O`` discipline (vv Mode 8): EVERY assertion below is a
``np.testing.*`` / ``pytest.fail`` FUNCTION CALL, never a bare
``assert`` — so it fires under the canonical ``python -O`` invocation.
Do NOT rewrite any of these as bare ``assert`` (it would compile to a
NO-OP and silently pass).

API ASSUMPTION (the C3 carve target — confirm the final spelling at
review): ``SweepDependencyGraph.from_cartesian(shape, *, label)`` where
``shape`` is a ``tuple[int, ...]`` (``len == d``) and ``label`` is a
``d``-tuple sign signature (the dimension-generic successor of the
2-element :class:`OctantLabel`). At ``d == 2`` it MUST reproduce
``from_cartesian_2d(nx=…, ny=…, label=OctantLabel(sx, sy))`` — that
d=2 equivalence is pinned by ``test_d2_from_cartesian_matches_legacy``
below.

If the C3 carve names the generic builder or the generic label
differently, ONLY the three module-level helpers
(``_build_nd``, ``_levels_of``, ``_octant_labels``) need updating; the
assertion bodies are spelling-independent.
"""
from __future__ import annotations

import itertools

import numpy as np
import pytest

from orpheus.sn.spatial.diamond import DiamondDifference
from orpheus.sn.sweep_graph import OctantLabel, SweepDependencyGraph


pytestmark = pytest.mark.foundation


# ``from_cartesian`` (the d-generic spine builder) is the C3 carve
# target — NOT yet landed (only ``from_cartesian_2d`` exists). Every pin
# that BUILDS a generic DAG is gated ``xfail(strict=False)`` on its
# absence so the C3 test set stays GREEN pre-carve and FLIPS to xpass
# automatically the moment the builder lands (``feedback_vv_tagging``
# idiom: strict=False + reason pointing at the unlocking API). The pure
# hand-oracle pins (octant counts, kernel closures) need no builder and
# run TODAY.
_HAS_GENERIC_BUILDER = hasattr(SweepDependencyGraph, "from_cartesian")
_needs_spine = pytest.mark.xfail(
    not _HAS_GENERIC_BUILDER,
    reason=(
        "SweepDependencyGraph.from_cartesian (the C3 d-generic spine "
        "builder) not yet landed — unlocks when the carve generalises "
        "from_cartesian_2d. strict=False ⟹ auto-xpass on landing."
    ),
    strict=False,
)


# ─── C3-carve adapters (UPDATE ONLY THESE if the final API differs) ───


def _octant_labels(ndim: int):
    """The ``2^ndim`` sign signatures, dimension-generic."""
    return list(itertools.product((-1, +1), repeat=ndim))


def _build_nd(shape: tuple[int, ...], signs: tuple[int, ...]):
    """Build the per-octant DAG for an ``ndim``-tuple ``shape``.

    C3 API: ``from_cartesian(shape, *, label=OctantLabel(signs))`` where
    ``OctantLabel`` carries the d-tuple ``signs`` signature.
    """
    return SweepDependencyGraph.from_cartesian(shape, label=OctantLabel(signs))


def _levels_of(graph) -> tuple:
    """The per-level cell-index tuple of the built graph (``graph.levels``)."""
    return graph.levels


# ─── B1 — octant enumeration: count == 2^d ──────────────────────────


@pytest.mark.parametrize("ndim,expected", [(1, 2), (2, 4), (3, 8)])
def test_octant_count_is_2_pow_ndim(ndim, expected):
    """B1: ``itertools.product((-1,+1), repeat=d)`` yields ``2^d`` octants."""
    labels = _octant_labels(ndim)
    np.testing.assert_array_equal(
        len(labels), expected,
        err_msg=f"d={ndim}: octant count {len(labels)} != 2^{ndim}={expected}",
    )
    # Distinctness — no duplicate sign tuple.
    if len(set(labels)) != expected:
        pytest.fail(f"d={ndim}: octant labels not distinct: {labels}")


# ─── B2 — the DAG builds at d=3 (synthetic shape admission) ──────────


@_needs_spine
@pytest.mark.parametrize("shape", [(2, 3, 2), (3, 3, 3), (4, 2, 3)])
def test_d3_dag_builds_for_every_octant(shape):
    """B2: ``from_cartesian`` builds a DAG for all 8 d=3 octants, no raise."""
    for signs in _octant_labels(3):
        graph = _build_nd(shape, signs)  # MUST NOT raise
        if graph is None:
            pytest.fail(f"shape={shape} signs={signs}: build returned None")


# ─── B3 — number of DAG levels == Σ(n_axis − 1) + 1 ─────────────────


@_needs_spine
@pytest.mark.parametrize("shape,n_levels", [
    ((2, 3, 2), 5),     # (2-1)+(3-1)+(2-1)+1 = 1+2+1+1 = 5
    ((3, 3, 3), 7),     # (2)+(2)+(2)+1 = 7
    ((4, 2, 3), 7),     # (3)+(1)+(2)+1 = 7
])
def test_d3_level_count(shape, n_levels):
    """B3: the anti-hyperplane DAG has ``Σ(n_axis−1)+1`` levels at d=3.

    Generalises the d=2 ``nx+ny−1``. The canonical ``(+1,+1,+1)``
    orientation; reversed-sign octants have the SAME level count
    (index reversal is a relabel, not a re-partition).
    """
    for signs in _octant_labels(3):
        graph = _build_nd(shape, signs)
        got = len(_levels_of(graph))
        np.testing.assert_array_equal(
            got, n_levels,
            err_msg=(
                f"shape={shape} signs={signs}: {got} levels, "
                f"expected Σ(n−1)+1 = {n_levels}"
            ),
        )


# ─── B4 — diagonal-hyperplane membership i+j+k == ℓ + total coverage ─


@_needs_spine
def test_d3_level_membership_is_index_sum_hyperplane():
    """B4: under ``(+1,+1,+1)`` each level ``ℓ`` holds exactly the cells
    with ``i+j+k == ℓ``; every cell appears exactly once; counts match
    the closed-form anti-hyperplane cardinality.

    For ``(2,3,2)`` the per-level counts are ``[1, 3, 4, 3, 1]``
    (sum = 12 = nx·ny·nz). This is the d=3 successor of the d=2
    anti-diagonal ``i+j == k`` coverage (``test_every_cell_visited_once``).
    """
    shape = (2, 3, 2)
    expected_counts = [1, 3, 4, 3, 1]
    graph = _build_nd(shape, (+1, +1, +1))
    levels = _levels_of(graph)

    # Level membership: every cell on level ℓ satisfies i+j+k == ℓ.
    seen: set[tuple[int, ...]] = set()
    counts: list[int] = []
    for ell, idx_arrays in enumerate(levels):
        # idx_arrays is the ndim-tuple of equal-length index arrays.
        cols = [np.asarray(a) for a in idx_arrays]
        n_cells = cols[0].size
        counts.append(int(n_cells))
        for c in range(n_cells):
            cell = tuple(int(col[c]) for col in cols)
            if sum(cell) != ell:
                pytest.fail(
                    f"cell {cell} on level {ell} has index-sum "
                    f"{sum(cell)} != {ell}"
                )
            if cell in seen:
                pytest.fail(f"cell {cell} visited twice (levels not disjoint)")
            seen.add(cell)

    np.testing.assert_array_equal(
        counts, expected_counts,
        err_msg=f"per-level cell counts {counts} != {expected_counts}",
    )
    np.testing.assert_array_equal(
        len(seen), int(np.prod(shape)),
        err_msg=f"total cells {len(seen)} != nx·ny·nz = {int(np.prod(shape))}",
    )


# ─── B5 — diamond kernel admits a 3-tuple of streaming coefficients ──


def test_diamond_kernel_admits_three_axis_streaming_sum():
    r"""B5: the WDD cell kernel's ``denom = Σ_t + Σ_axis s_axis`` and
    ``psi_avg = (Q + Σ_axis s_axis·in_axis)/denom`` reduce correctly when
    given a 3-tuple of per-axis streaming coefficients.

    This is the d-generic successor of the explicit ``sx``/``sy`` 2-axis
    form. The structural claim: the kernel sums over ``d`` axes and the
    per-axis closure ``psi_out_axis = 2·psi_avg − psi_in_axis`` holds for
    each of the 3 axes independently. NO solve — a single hand-checkable
    cell evaluation.

    ASSUMED C3 API: ``DiamondDifference`` exposes a d-generic batch kernel
    taking ``psi_in: tuple/stack over axes`` and ``s: tuple/stack over
    axes``. If the carve keeps named ``cell_kernel_batch`` but generalises
    its signature, update the CALL below only; the hand-computed oracle
    (``_hand``) is API-independent and the load-bearing reference.
    """
    sigt = 0.5
    s = (0.3, 0.7, 1.1)                  # 3-axis streaming coeffs
    psi_in = (1.0, 2.0, 0.5)            # incoming face flux per axis
    Q = 0.8

    # Structurally-independent hand oracle (scalar, d-generic by sum).
    denom_ref = sigt + sum(s)
    psi_avg_ref = (Q + sum(si * ai for si, ai in zip(s, psi_in))) / denom_ref
    psi_out_ref = tuple(2.0 * psi_avg_ref - ai for ai in psi_in)

    np.testing.assert_allclose(
        denom_ref, 2.6, rtol=0, atol=1e-15,
        err_msg="hand-oracle denom drifted (sanity)",
    )

    # The REAL d-generic kernel (C3.2): promote each scalar to the
    # (N_oct=1, ng=1, n_diag=1) batch shape — ``s`` is (N_oct, 1, n_diag).
    dd = DiamondDifference()
    psi_in_b = tuple(np.full((1, 1, 1), a) for a in psi_in)
    s_b = tuple(np.full((1, 1, 1), si) for si in s)
    psi_avg, psi_out = dd.cell_kernel_batch(
        psi_in=psi_in_b, s_axes=s_b,
        sigt_cells=np.full((1, 1), sigt),      # (ng, n_diag)
        Q_cells=np.full((1, 1, 1), Q),         # (N_oct, ng, n_diag)
    )
    np.testing.assert_allclose(
        psi_avg.ravel()[0], psi_avg_ref, rtol=0, atol=1e-14,
        err_msg="d=3 kernel psi_avg != hand-oracle (denom = sigt + Σ_a s_a)",
    )
    if len(psi_out) != 3:
        pytest.fail(f"d=3 kernel returned {len(psi_out)} outgoing faces != 3")
    for a, (out_a, ref_a) in enumerate(zip(psi_out, psi_out_ref)):
        np.testing.assert_allclose(
            out_a.ravel()[0], ref_a, rtol=0, atol=1e-14,
            err_msg=f"d=3 kernel psi_out[{a}] != 2·psi_avg − psi_in[{a}]",
        )


# ─── B6 — d=1 is a REAL spine compute path (NEW-3) ──────────────────
#
# d = 1 is NOT synthetic shape-admission: it is the production wavefront
# spine instantiated at one axis, validated against the cumprod
# optimization in the companion file
# ``test_wavefront_cumprod_equivalence.py``. These pins nail the d = 1
# DAG STRUCTURE (the structural precondition the equivalence oracle
# assumes): 2 octants, nx levels, one cell per level, single-axis
# closure.


def test_d1_octant_count_is_two():
    """B6a: d=1 has exactly ``2^1 = 2`` octants (the ±μ sweep directions).

    This is the degenerate-but-real base case of B1. A 1-D sweep has
    one streaming axis ⟹ two sign signatures ``(-1,)`` and ``(+1,)``.
    """
    labels = _octant_labels(1)
    np.testing.assert_array_equal(
        labels, [(-1,), (+1,)],
        err_msg=f"d=1 octant labels {labels} != [(-1,), (+1,)]",
    )


@_needs_spine
@pytest.mark.parametrize("nx", [5, 8, 12])
def test_d1_dag_has_nx_levels_one_cell_each(nx):
    """B6b: ``from_cartesian((nx,))`` builds an ``nx``-level DAG with
    exactly one cell per level, for BOTH octants.

    The d=1 specialisation of B3 (``Σ(n_axis−1)+1 = (nx−1)+1 = nx``)
    AND B4 (per-level membership ``i == ℓ`` ⟹ singleton levels). A 1-D
    wavefront is a pure chain: level ``ℓ`` holds only cell ``i = ℓ``.
    Each level being a singleton is the structural reason the spine's
    per-level ``cell_kernel_batch`` reduces to the same scalar
    recurrence the cumprod scan computes (the equivalence-oracle
    precondition).
    """
    for signs in _octant_labels(1):
        graph = _build_nd((nx,), signs)
        levels = _levels_of(graph)
        np.testing.assert_array_equal(
            len(levels), nx,
            err_msg=f"nx={nx} signs={signs}: {len(levels)} levels != nx={nx}",
        )
        seen: set[int] = set()
        for ell, idx_arrays in enumerate(levels):
            cols = [np.asarray(a) for a in idx_arrays]
            np.testing.assert_array_equal(
                cols[0].size, 1,
                err_msg=(
                    f"nx={nx} signs={signs} level {ell}: "
                    f"{cols[0].size} cells, expected 1 (singleton chain)"
                ),
            )
            # Canonical (+1,) orientation: level ℓ holds cell i == ℓ.
            cell_i = int(cols[0][0])
            if signs == (+1,) and cell_i != ell:
                pytest.fail(
                    f"nx={nx} (+1,): level {ell} holds cell {cell_i} != {ell}"
                )
            if cell_i in seen:
                pytest.fail(f"nx={nx} signs={signs}: cell {cell_i} visited twice")
            seen.add(cell_i)
        np.testing.assert_array_equal(
            len(seen), nx,
            err_msg=f"nx={nx} signs={signs}: {len(seen)} cells != nx={nx}",
        )


def test_d1_diamond_kernel_single_axis_streaming_sum():
    r"""B6c: the d-generic cell kernel at d=1 reduces to the known 1-D DD
    closure: ``denom = σ_t + s_x``, ``psi_avg = (Q + s_x·in_x)/denom``,
    ``psi_out = 2·psi_avg − in_x``.

    The d=1 specialisation of B5: ``Σ_axis`` over ONE axis. This is the
    SAME scalar recurrence (per cell) that ``ordinate_scan`` computes in
    the cumprod path — but here we pin only that the spine kernel
    EVALUATES it correctly, NOT the FP-association order (that is the
    nulp-bounded equivalence in the oracle file). Hand oracle, no solve.

    ASSUMED C3 API: same d-generic kernel as B5; update the CALL only.
    """
    sigt = 0.5
    s = (0.7,)                          # single-axis streaming coeff
    psi_in = (1.3,)
    Q = 0.8

    denom_ref = sigt + sum(s)           # 1.2
    psi_avg_ref = (Q + s[0] * psi_in[0]) / denom_ref
    psi_out_ref = 2.0 * psi_avg_ref - psi_in[0]

    np.testing.assert_allclose(
        denom_ref, 1.2, rtol=0, atol=1e-15,
        err_msg="d=1 hand-oracle denom drifted (sanity)",
    )
    # The REAL d-generic kernel at d=1 (single streaming axis).
    dd = DiamondDifference()
    psi_avg, psi_out = dd.cell_kernel_batch(
        psi_in=(np.full((1, 1, 1), psi_in[0]),),
        s_axes=(np.full((1, 1, 1), s[0]),),
        sigt_cells=np.full((1, 1), sigt),
        Q_cells=np.full((1, 1, 1), Q),
    )
    np.testing.assert_allclose(
        psi_avg.ravel()[0], psi_avg_ref, rtol=0, atol=1e-14,
        err_msg="d=1 kernel psi_avg != hand-oracle (denom = sigt + s_x)",
    )
    if len(psi_out) != 1:
        pytest.fail(f"d=1 kernel returned {len(psi_out)} faces != 1")
    np.testing.assert_allclose(
        psi_out[0].ravel()[0], psi_out_ref, rtol=0, atol=1e-14,
        err_msg="d=1 kernel psi_out != 2·psi_avg − psi_in",
    )


# ─── B7 — the d-generic WALK (apply / residual) runs at d=1 and d=3 ──
#
# B5/B6 pin the d-generic cell KERNEL; B7 pins the d-generic graph WALK
# end-to-end (``SweepDependencyGraph.apply`` / ``.residual``). The check is
# the geometry-agnostic matvec==sweep round-trip: the operator residual of
# the swept solution vanishes (``test_sweep_graph.py``'s d=2 contract lifted
# to d=1 and d=3). This is what makes "d=1/d=3 CAN walk" (the C3.2b
# deliverable) a TESTED fact, not just a kernel admission — it exercises the
# per-axis advanced-index gather/scatter (``_cell_face_selector``) and the
# ``(…, *cell_idx)`` flux scatter at d ≠ 2.


def _walk_roundtrip_residual(shape, signs, *, ng, N_oct, seed):
    """Solve via ``graph.apply``, then re-apply the operator at the swept
    solution via ``graph.residual`` from the SAME boundary seeds; return the
    per-cell residual (must vanish — the apply↔solve round-trip).

    Geometry-agnostic: ``shape`` is any ``d``-tuple, the face buffers are
    sized ``n_a + 1`` along axis ``a``, and the streaming coefficients /
    source / ordinate count are arbitrary (NO real quadrature — this is the
    synthetic-shape admission of the WALK, the d≠2 sibling of B5/B6).
    """
    rng = np.random.default_rng(seed)
    d = len(shape)

    def _face_shape(a):
        s = list(shape)
        s[a] += 1
        return (N_oct, ng, *s)

    psi_faces = tuple(rng.standard_normal(_face_shape(a)) for a in range(d))
    Q = rng.standard_normal((N_oct, ng, *shape))
    sig_t = rng.uniform(0.1, 0.5, size=(ng, *shape))
    str_axes = tuple(
        rng.uniform(0.1, 1.0, size=(N_oct, shape[a])) for a in range(d)
    )
    weights = rng.uniform(0.5, 1.5, size=N_oct)
    graph = _build_nd(shape, signs)

    # SOLVE — forward-substitute update_batch along the DAG.
    ang = np.zeros((N_oct, ng, *shape))
    scal = np.zeros((ng, *shape))
    graph.apply(
        cell_update=DiamondDifference(),
        psi_faces_octant=tuple(pf.copy() for pf in psi_faces),
        Q_octant=Q, sig_t=sig_t,
        str_axes_octant=str_axes,
        weights_octant=weights,
        angular_flux_octant=ang, scalar_flux_buf=scal,
    )
    # APPLY at the swept solution — FRESH copies with the SAME inflow seeds
    # (apply only scatters OUTGOING faces; the seeds carry the inflow).
    residual = np.zeros((N_oct, ng, *shape))
    graph.residual(
        cell_update=DiamondDifference(),
        psi_faces_octant=tuple(pf.copy() for pf in psi_faces),
        psi_avg_probe_octant=ang,
        Q_octant=Q, sig_t=sig_t,
        str_axes_octant=str_axes,
        residual_octant=residual,
    )
    return residual


@_needs_spine
@pytest.mark.parametrize("signs", _octant_labels(1))
def test_d1_walk_residual_vanishes_at_apply_solution(signs):
    """B7(d=1): the d=1 ``graph.apply`` / ``graph.residual`` walks run (the
    1-D chain) and the operator residual of the swept solution vanishes — the
    matvec==sweep contract on the single-axis spine. Proves the d=1 WALK (not
    just the kernel) is wired; the cumprod oracle (separate file) then pins
    the d=1 spine ≡ the cumprod optimization (C3.4)."""
    residual = _walk_roundtrip_residual(
        (6,), signs, ng=2, N_oct=4, seed=101 + signs[0],
    )
    np.testing.assert_allclose(
        residual, 0.0, atol=1e-11,
        err_msg=f"d=1 signs={signs}: round-trip residual did not vanish",
    )


@_needs_spine
@pytest.mark.parametrize("signs", _octant_labels(3))
def test_d3_walk_residual_vanishes_at_apply_solution(signs):
    """B7(d=3): the d=3 walks run on a synthetic shape (NO 3-D quadrature —
    arbitrary N_oct / streaming / source) and the round-trip residual
    vanishes. Proves the anti-hyperplane forward substitution + the 3-axis
    diamond closure compose into a self-consistent WALK at d=3 (C3 invariant
    #3: synthetic-3-D shape admission, here lifted from the kernel to the
    walk)."""
    residual = _walk_roundtrip_residual(
        (3, 2, 3), signs, ng=2, N_oct=3, seed=303 + sum(signs),
    )
    np.testing.assert_allclose(
        residual, 0.0, atol=1e-11,
        err_msg=f"d=3 signs={signs}: round-trip residual did not vanish",
    )


# ─── C — d=2 equivalence: from_cartesian(2) == from_cartesian_2d ─────


@_needs_spine
@pytest.mark.parametrize("nx,ny", [(3, 4), (4, 3), (5, 7), (2, 3)])
@pytest.mark.parametrize("sx,sy", [(+1, +1), (+1, -1), (-1, +1), (-1, -1)])
def test_d2_from_cartesian_matches_legacy(nx, ny, sx, sy):
    """C: the generic ``from_cartesian((nx,ny))`` reproduces the OLD
    ``from_cartesian_2d`` level structure + face conventions bit-for-bit.

    Retire-the-predecessor-with-a-pin (``feedback_retirement_means_test_migration``):
    captures the legacy structure BEFORE it is retired and asserts the new
    builder matches. Non-square ``(5,7)`` / ``(4,3)`` cases catch an x↔y
    axis swap in the generic level partition (AI failure mode 2) that a
    square grid cannot see.

    NOTE: if C3 retires ``from_cartesian_2d`` entirely, replace ``legacy``
    with a frozen golden captured at the PRE-C3 commit (the structure is
    small — levels + 4 face ints — and serialises cleanly to an ``.npz``).
    """
    legacy = SweepDependencyGraph.from_cartesian_2d(
        nx=nx, ny=ny, label=OctantLabel((sx, sy)),
    )
    generic = _build_nd((nx, ny), (sx, sy))

    # Same number of levels.
    np.testing.assert_array_equal(
        len(generic.levels), len(legacy.levels),
        err_msg=f"({nx},{ny}) signs=({sx},{sy}): level count drift",
    )
    # Same cells per level, in the same per-octant traversal order.
    for k, (leg_lvl, gen_lvl) in enumerate(zip(legacy.levels, generic.levels)):
        leg_ii, leg_jj = leg_lvl
        gen_axes = gen_lvl
        np.testing.assert_array_equal(
            np.asarray(gen_axes[0]), np.asarray(leg_ii),
            err_msg=f"level {k}: x-index drift (axis swap?)",
        )
        np.testing.assert_array_equal(
            np.asarray(gen_axes[1]), np.asarray(leg_jj),
            err_msg=f"level {k}: y-index drift (axis swap?)",
        )
    # Same face-offset conventions (the upwind orientation crosswalk).
    np.testing.assert_array_equal(
        [generic.face_in_x, generic.face_out_x,
         generic.face_in_y, generic.face_out_y],
        [legacy.face_in_x, legacy.face_out_x,
         legacy.face_in_y, legacy.face_out_y],
        err_msg="face-offset convention drift between generic and legacy",
    )
