r"""P0 performance gates — the composition-cost floor for the
operator · splitting · realization campaign (plan §5, constraint O-6).

Normative spec: ``.claude/plans/campaign_verification_plan.md`` §5 (P-1
call-count scaling, P-2 normalised throughput, P-3 allocation ceiling,
P-4 baseline capture) and §6 O-6 (**stage 3 must not merge without these
green and the P-4 baseline captured on the immediately preceding
commit**).

Why this file exists (lesson L16, the precedent)
================================================

A refactor moved hoisted vectorised work into a per-cell Python fold and
cost **10-20x on slab / ~6x on the full suite with every correctness gate
green** — the L0 streaming-equilibrium anchor was 26/26 PASS while the SN
suite ran three hours to 28 %.  Stage 3's ``A_ij = R_i A J_j`` block
composition is the same shape of risk: every ``R_i`` / ``J_j`` is an
opportunity to materialise a slice per block per cell.  No correctness
gate in this repo measures wall clock, so the regression is structurally
invisible to them.

The three legs, and which one is the catcher
============================================

**P-1 (call-count scaling) is the catcher.**  It wraps the batched leaf
kernel in-process and counts entries while the mesh is refined.  It is a
Mode-11-style wrap: it simultaneously proves the measured path *is
executed* and that its call arity is *structural, not per-cell*.  It is
deterministic — no timing noise, no machine dependence, no CI flake.
P-2 (wall clock) and P-3 (allocation) are corroborating legs; only P-3
shares P-1's exact reproducibility.

Which leaf kernel — the spec's illustrative name is STALE
=========================================================

Plan §5.2 names ``DiamondDifference.update_batch``.  **That method does
not exist**: it was renamed at S6.4(e) to the packet-free pair
``cell_kernel_batch`` (solve direction) / ``residual_kernel_batch``
(apply direction); ``git grep update_batch orpheus`` returns one
docstring mention and no code.  The genuine batched leaf on the path
``A.apply(x)`` takes is

    :meth:`orpheus.transport.spatial.diamond.DiamondDifference.residual_kernel_batch`

reached from ``StreamingCollisionOperator.apply`` ->
``loss_representation.loss_action`` -> (1-D) ``_OneDimScanWalk._apply_walk``
per-cell ``visit`` / (2-D Cartesian) ``ScanMarch`` row-march.
:func:`test_p1_wrap_fires_and_siblings_stay_cold` is the non-vacuity
proof: the wrap fires (count > 0) on a real ``A.apply(x)`` on BOTH
geometries, *and* the three plausible sibling kernels (``cell_kernel_batch``,
the per-cell ``residual``, the per-cell ``update``) stay at exactly 0 —
so the counted method is THE apply-direction leaf, not merely *a* method
that happens to run.

The pre-carve finding (do NOT read this file as "the count is invariant")
=========================================================================

Measured on this tree, the two production geometries behave DIFFERENTLY,
and only one of them is mesh-invariant:

* **2-D Cartesian** (``ScanMarch``, the production default since the S6.9
  Fork-B2 flip): ``calls = 8 * ny`` — **exactly invariant in nx** (the
  x-axis is inside the Blelloch ``ordinate_scan``), invariant in the
  quadrature order and in the group count, linear in ny (the marched
  axis).  ``calls / n_cells`` therefore HALVES under isotropic
  refinement.  This is the hoisted regime P-1a gates with a strict ``==``.
* **1-D slab** (``CumprodScan``): ``calls = 2 * nx`` — **linear in the
  cell count**.  The 1-D *apply* direction is ALREADY a per-cell Python
  fold on this tree (by design: ``_OneDimScanWalk``'s docstring calls it
  "the apply-loop frame ... ONE orientation-parametrized per-cell loop"),
  while the 1-D *solve* direction on the same mesh is the vectorised
  scan at a mesh-invariant **2** ``ordinate_scan`` calls.  That asymmetry
  is the L16 shape already present, pre-carve.  P-1d pins it as a
  regression floor (not as an invariance claim), and P-1e ships the
  solve-direction invariance as the **control leg** proving the counting
  instrument can see invariance when it is there.

Every gate's teeth, mutation-verified in process
===============================================

Each row was applied by rebinding a class attribute in process (NEVER a
``git checkout`` — this tree carries uncommitted-by-policy files) and the
named gate confirmed RED, with the unmutated control confirmed GREEN in
the same run (the ERR-067 closure-trap discipline: a still-broken
baseline mimics "caught").

===============================  ==========================================
mutation                         gate -> observed RED
===============================  ==========================================
per-cell (column) fold           P-1a  80,80,80,80 -> 720,1360,2640,5200
per-cell (column) fold           P-1c  calls/cell 1.00,0.50,0.25 -> 9.0,8.5
per-ordinate fold                P-1b  80,80,80 -> 640,1920,3840
per-ordinate fold                P-1d  40,80,160,320 -> 200,400,800,1600
per-cell ``ordinate_scan`` fold  P-1e  2,2,2,2 -> 42,82,162,322
per-cell fold                    P-2   ratio 23.4 -> 362.1
+1 held full-field temporary     P-3   3.031 -> 4.035
counted method never entered     P-1 non-vacuity -> "NEVER FIRED ... VACUOUS"
``_COST_NX + 1``                 P-4   fingerprint drift message
===============================  ==========================================

**The headline measurement**: the per-cell fold is EXACTLY value-identical
(``np.array_equal`` True, ``max|diff|`` 0.0e+00) — the DD residual kernel
is cell-local, so splitting its batch column-wise changes the arity and
nothing else.  No value gate, at any tolerance, on any reference, can
ever see this regression class.  That is why P-1 is the catcher and P-2
is only corroboration.

Fixtures are FROZEN and LOCAL — deliberately
============================================

The package's shared ``_config`` supplies the campaign's *conventions*
(direct ``Mixture`` construction with P1 + asymmetric ``SigS``; the one
posing site; the carrier-generic probe state) and this module imports
exactly those.  It does NOT use ``_config``'s meshes: a performance
baseline's **fixture is part of the baseline**, and ``_config``'s meshes
are shared with the correctness gates, which are free to retune them (a
Mode-9 fix bumping S4 -> S6, a third region) — that would silently
invalidate every committed constant below.  So the sizes live here,
frozen, and :func:`test_p4_fixture_fingerprint_matches_the_baseline`
reds with a "re-capture the baseline" message if any of them drifts.

Mandatory configuration (plan §8), and what a cost gate actually needs
=====================================================================

Every fixture here is >=2G, heterogeneous (2 regions), non-uniform ``h``,
non-square ``nx != ny`` in 2-D, mixed reflective/vacuum BC, P1
anisotropic scattering built through the direct ``Mixture`` route, and
probed with a fixed-seed random non-flat state (bulk AND trace).  For the
*cost* legs (P-2/P-3) the binding reason is term ACTIVATION: a P0-only or
zero-``SigS`` mixture makes ``S.apply`` a near-no-op and the measured
cost stops covering the operator whose block-slicing stage 3 introduces.
For the *count* leg (P-1) the arity is provably independent of the
materials, the group count and the quadrature order — measured, and
pinned by :func:`test_p1_call_count_is_invariant_in_ordinates_and_groups`.

Run (HOST -> ``.venv/bin/python``; canonical ``-O``; SERIAL, xdist is
unstable on this venv)::

    .venv/bin/python -O -m pytest \
      tests/sn/architecture/test_composition_cost.py -q -p no:cacheprovider

vv-principles compliance
========================

* **Mode 8**: every assertion is ``np.testing.assert_*`` / ``pytest.fail``
  — never a bare ``assert``, in test bodies or helpers alike, so every
  gate fires under the canonical ``-O``.
* **Mode 11**: the counted method is proved to be on the production path
  before any count is asserted (:func:`test_p1_wrap_fires_and_siblings_stay_cold`).
  A counter that stays 0 because production routes around the wrap makes
  the whole gate vacuous; that is the single most likely way this gate
  fails silently.
* **foundation, no ``verifies(...)``**: these are software/architecture
  invariants with no theory ``:label:`` (the verifies-perp-level doctrine).
* **No ``catches(ERR-NNN)``**: L16 is a lessons entry, not a catalogued
  ERR; attaching a coverage claim without a mutation-verified red would
  be a phantom edge.
"""
from __future__ import annotations

import functools
import time
import tracemalloc
from contextlib import contextmanager
from typing import TYPE_CHECKING, Any, Iterator

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn import loss_representation as _loss_representation
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.sn.operators.streaming import StreamingCollisionOperator
from orpheus.transport.spatial.diamond import DiamondDifference
from tests.sn.architecture._config import (
    anisotropic_mixture, random_state, record_for, system_a,
)

if TYPE_CHECKING:
    from orpheus.data.macro_xs.mixture import Mixture
    from orpheus.numerics.coupled_system import CoupledField
    from orpheus.numerics.operator import LinearOperator

pytestmark = pytest.mark.foundation


# ═══════════════════════════════════════════════════════════════════════
# P-4 — the pre-carve baseline.  Every constant below was MEASURED on the
# commit named in the provenance banner; a perf gate whose baseline was
# captured AFTER the regression is worthless (plan §5.2 P-4).
# ═══════════════════════════════════════════════════════════════════════

# ---------------------------------------------------------------------
# PROVENANCE — measured on branch refactor/operator-strategy-layers
# @ 361471ec (``git diff 360a8087..361471ec -- orpheus`` is EMPTY, so
# every number below is equally the 360a8087 reading), host .venv
# Py 3.14.3 (darwin 25.4.0), SERIAL (xdist is unstable here), python -O.
#
#   CONTENTION-IMMUNE (final): every P-1 call count; the P-3 allocation
#   ratio.  Both are exact integers / exactly reproducible byte counts —
#   they do not move with machine load.
#
#   PROVISIONAL — recapture on a quiescent tree: _COMPOSITION_OVERHEAD_MAX
#   only.  The host was CONTENDED during capture (concurrent edits +
#   builds), so the wall-clock ratio it derives from is not trustworthy
#   as a tight bound.  The measured in-process band was 23.35 - 24.41
#   (6 trials, min-of-7 apply / min-of-25 calibration); the constant is
#   set at 3x the top of that band so a contended run cannot false-red.
#   TIGHTEN IT once a quiescent measurement exists.
# ---------------------------------------------------------------------

#: 1-D walk-DAG legs (one inward, one outward) — the slab apply calls the
#: leaf kernel once per (leg, cell).  MEASURED: 40/80/160/320 at
#: nx = 20/40/80/160; invariant in the quadrature order (S2..S16) and in ng.
_SLAB_CALLS_PER_CELL = 2

#: Octant-sweep units of the quadrature (the 2^3 sign combinations of
#: (mu_x, mu_y, mu_z)) — the 2-D row-march calls the leaf kernel once per
#: (octant, y-row).  MEASURED: 80 at ny=10 for nx in {8,16,32}; 88 at
#: ny=11; 160 at ny=20; 320 at ny=40.
_CART2D_CALLS_PER_ROW = 8

#: ``ordinate_scan`` calls in ONE 1-D sweep — "Exactly 2 scan calls per
#: sweep" (``loss_representation`` line ~3894), i.e. one per chain.
#: MEASURED mesh-invariant: 2 at nx = 20/40/80/160.
_SLAB_SOLVE_SCAN_CALLS = 2

#: P-2 CLAIM: one ``A.apply(x)`` on the frozen cost fixture costs at most
#: this multiple of a dense ``(m x m)`` BLAS contraction carrying the SAME
#: nominal FLOP count (``m = round(sqrt(_CALIBRATION_FLOPS_PER_DOF *
#: n_dof))``).  It is a CLAIM about relative cost, not a magic number:
#: dividing by an in-process calibration normalises out the host's raw
#: numpy throughput, which an absolute-ms gate (the
#: ``tests/sn/sweep/core/test_cache.py:553`` precedent) cannot do.
#: MEASURED PRE-CARVE: 23.4 - 24.4.  PROVISIONAL — see the banner.
_COMPOSITION_OVERHEAD_MAX = 72.0

#: Nominal FLOPs per degree of freedom for a DD cell balance (a handful of
#: fused multiply-adds per DOF).  Sets the calibration size so the dense
#: reference does comparable work; it is a fixed convention, not a tuned
#: parameter — changing it changes what the ratio MEANS.
_CALIBRATION_FLOPS_PER_DOF = 16

#: P-3 CLAIM: peak bytes allocated by ONE ``A.apply(x)``, divided by
#: ``n_dof * 8``.  A value of 1.0 is the structural floor (the output
#: field itself); the measured ~3.03 is roughly three full-field
#: temporaries.  Catches "materialise a slice per block per cell"
#: directly, and unlike wall clock it is exactly reproducible.
#: MEASURED PRE-CARVE: 3.0299 - 3.0328 (spread 0.03 %) — FINAL, not
#: provisional: this leg is contention-immune.
#:
#: 4.0, NOT the plan's suggested 6.0.  The suggestion predates any
#: measurement; with the baseline at 3.03 and a 0.03 % spread, the
#: sensitivity is exactly one held full-field temporary.  MEASURED
#: (extra full-field temporaries held alive across one apply):
#:
#:      +0 -> 3.041   +1 -> 4.035   +2 -> 5.033   +3 -> 6.033
#:
#: so 4.0 reds on the FIRST extra field while 6.0 would need THREE.  The
#: block composition materialising one field per block is the exact
#: mechanism this leg exists to catch, so the tighter constant is the
#: honest one; a principled change that legitimately costs a field must
#: re-baseline consciously (lesson L7).
_ALLOC_FACTOR_MAX = 4.0

#: Timing reps.  ``min`` is the estimator, NEVER the mean: the minimum is
#: the only order statistic that is not contaminated by a descheduling
#: event, and this host is not quiescent.
_N_REPS_APPLY = 7
_N_REPS_CALIBRATION = 25


# ═══════════════════════════════════════════════════════════════════════
# The FROZEN fixtures.  Sizes are part of the baseline (see the module
# docstring); the fingerprint gate reds if any of them drifts.
# ═══════════════════════════════════════════════════════════════════════

#: Two heterogeneous 2G materials, P1 anisotropic, ASYMMETRIC ``SigS`` in
#: both Legendre orders — so ``S.apply`` does real work in the cost legs
#: and no term the block slicing will touch is nulled (plan §8).
_REGION_FUEL: tuple[Any, ...] = (
    [1.4, 2.2], [[0.35, 0.28], [0.03, 0.95]], [[0.09, 0.06], [0.01, 0.22]],
)
_REGION_MODERATOR: tuple[Any, ...] = (
    [0.8, 2.7], [[0.24, 0.44], [0.01, 1.20]], [[0.06, 0.10], [0.00, 0.33]],
)

#: The P-1 2-D ladder.  ``ny`` is FIXED while ``nx`` refines: nx is the
#: SCANNED axis, so this is the rung set on which the count must be
#: exactly invariant (P-1a, the catcher).
_LADDER_NY = 10
_LADDER_NX = (8, 16, 32, 64)

#: Isotropic refinement — n_cells x4 per rung while the count only x2, so
#: calls-per-cell halves (P-1c).
_LADDER_ISOTROPIC = ((8, 10), (16, 20), (32, 40))

#: The P-1 slab ladder (the finding + its control leg).
_LADDER_SLAB_NX = (20, 40, 80, 160)

#: The FROZEN cost fixture (P-2/P-3).  Large enough that fixed overhead
#: does not dominate the ratio, small enough to stay well under a second.
_COST_NX, _COST_NY = 32, 40

#: Fingerprint of the cost fixture — the scalars every committed constant
#: above was measured against.  MEASURED at capture.
_COST_FINGERPRINT = {
    "n_cells": 1280,
    "n_ordinates": 24,
    "n_groups": 2,
    "n_regions": 2,
    "n_dof": 68352,
}


def _stretched_edges(length: float, n_cells: int, *, ratio: float) -> np.ndarray:
    r"""Geometrically graded cell edges on ``[0, length]`` — NON-uniform ``h``.

    Plan §8: a uniform mesh lets a constant metric cancel and hides an
    ``x <-> y`` swap.  Grading is mild (``ratio`` per cell) so the fixture
    stays well-conditioned; the widths span ``ratio ** (n_cells - 1)``.
    """
    widths = ratio ** np.arange(n_cells, dtype=float)
    return np.concatenate(([0.0], np.cumsum(widths))) * (length / widths.sum())


def _two_region_materials() -> "dict[int, Mixture]":
    return {
        0: anisotropic_mixture(*_REGION_FUEL),
        1: anisotropic_mixture(*_REGION_MODERATOR),
    }


def _ng_generic_materials(ng: int) -> "dict[int, Mixture]":
    r"""Two heterogeneous regions at an arbitrary group count.

    Used ONLY by the group-invariance rung of P-1 (ng in {1, 4}), where
    the frozen 2G pair is not spellable.  Still P1 (an ``sig_s1`` block
    is mandatory — every fixture here is posed at ``scattering_order=1``,
    and a P0-only mixture would raise on the missing Legendre block) and
    still asymmetric off-diagonal, so nothing the arity depends on is
    special-cased relative to the frozen pair.
    """
    def block(diag: float, off: float) -> "list[list[float]]":
        matrix = np.full((ng, ng), off) + np.eye(ng) * diag
        matrix[np.tril_indices(ng, -1)] *= 0.1      # asymmetric up/down scatter
        return matrix.tolist()

    return {
        0: anisotropic_mixture([1.4] * ng, block(0.30, 0.10), block(0.08, 0.02)),
        1: anisotropic_mixture([0.8] * ng, block(0.20, 0.06), block(0.05, 0.01)),
    }


def _cart2d(nx: int, ny: int, *, ng: int = 2, sn_order: int = 4) -> SNMesh:
    r"""2-D Cartesian: non-square, heterogeneous along x, non-uniform ``h``,
    mixed reflective/vacuum BC, ``level_symmetric`` (genuine ``mu_y``;
    avoids the #214 ``mu_y == 0`` GL rank mismatch AND the ERR-056
    axis-aligned degeneracy).

    Shaped after ``tests/sn/operators/test_removal_form_matvec_sweep.py``'s
    ``_cart2d``, widened with the mesh ladder + grading this file needs.
    """
    mat_map = np.zeros((nx, ny), dtype=int)
    mat_map[nx // 2:, :] = 1
    mesh = Mesh2D(
        edges_x=_stretched_edges(2.0, nx, ratio=1.02),
        edges_y=_stretched_edges(3.0, ny, ratio=1.03),
        mat_map=mat_map,
        coord=CoordSystem.CARTESIAN,
        bc_xmin=BC("reflective"), bc_xmax=BC("vacuum"),
        bc_ymin=BC("reflective"), bc_ymax=BC("vacuum"),
    )
    materials = _two_region_materials() if ng == 2 else _ng_generic_materials(ng)
    return SNMesh(mesh, Quadrature.level_symmetric(sn_order=sn_order), materials)


def _slab(nx: int, *, n_ord: int = 8) -> SNMesh:
    r"""1-D Cartesian slab: heterogeneous (2 regions), non-uniform ``h``,
    mixed reflective/vacuum BC, 2G, P1.

    Shaped after the same file's ``_slab``, widened with grading and the
    two-region split.
    """
    mat_ids = np.zeros(nx, dtype=int)
    mat_ids[nx // 2:] = 1
    mesh = Mesh1D(
        edges=_stretched_edges(4.0, nx, ratio=1.01),
        mat_ids=mat_ids,
        coord=CoordSystem.CARTESIAN,
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    return SNMesh(
        mesh, Quadrature.gauss_legendre(n_ordinates=n_ord),
        _two_region_materials(),
    )


def _posed_apply(
    sn_mesh: SNMesh, *, seed: int = 20260728,
) -> "tuple[LinearOperator, CoupledField]":
    r"""``(A, x)`` for the production posed loss on ``sn_mesh``.

    ``A`` is the record's own ``loss`` grid — the object stage 3 replaces
    with ``A_ij = R_i A J_j``, so it is the object whose cost must be
    gated.  ``x`` is a fixed-seed random NON-flat state on ``A``'s own
    carrier (bulk AND trace populated: a flat probe nulls the streaming
    coupling and a zero trace nulls ``B``).
    """
    record = record_for(sn_mesh, scattering_order=1)
    return record.loss, random_state(record, seed=seed)


def _n_dof(state: "CoupledField") -> int:
    return int(state.to_flat().size)


# ═══════════════════════════════════════════════════════════════════════
# The in-process wrap (vv Mode 11).  Counts entries into a production
# method without touching production code.
# ═══════════════════════════════════════════════════════════════════════


class _CallCounter:
    """Mutable entry counter handed out by :func:`_count_calls`."""

    __slots__ = ("n",)

    def __init__(self) -> None:
        self.n = 0


@contextmanager
def _count_calls(owner: Any, method_name: str) -> "Iterator[_CallCounter]":
    r"""Wrap ``owner.method_name`` in process and count entries.

    Strictly stronger than a file-write probe (vv Mode 11 sharpening): it
    runs IN the test process, on the SAME object the production path
    constructs, so a path that routes around the wrapped method leaves
    the counter at 0 and cannot fake it.  ``owner`` may be a class or a
    module — the call sites of interest are a bound method
    (``scheme.residual_kernel_batch``) and a module-level function
    (``loss_representation.ordinate_scan``), and attribute rebinding
    covers both.

    The original is restored in ``finally``, so a failing assertion
    inside the ``with`` cannot leak a wrapped kernel into a later test.
    """
    counter = _CallCounter()
    original = getattr(owner, method_name)

    @functools.wraps(original)
    def _wrapped(*args: Any, **kwargs: Any) -> Any:
        counter.n += 1
        return original(*args, **kwargs)

    setattr(owner, method_name, _wrapped)
    try:
        yield counter
    finally:
        setattr(owner, method_name, original)


def _leaf_calls_for_apply(sn_mesh: SNMesh, *, method: str = "residual_kernel_batch") -> int:
    """Entries into ``DiamondDifference.<method>`` during ONE ``A.apply(x)``.

    One un-counted warm-up apply first, so a lazily-built cache (the
    two-stratum sweep cache, the per-octant DAG) is not counted as sweep
    work.  The counted apply is therefore the steady-state arity.
    """
    A, x = _posed_apply(sn_mesh)
    A.apply(x)                                   # warm the derived caches
    with _count_calls(DiamondDifference, method) as counter:
        A.apply(x)
    return counter.n


# ═══════════════════════════════════════════════════════════════════════
# P-1 — call-count scaling.  THE CATCHER.
# ═══════════════════════════════════════════════════════════════════════


def test_p1_wrap_fires_and_siblings_stay_cold() -> None:
    r"""[foundation] The wrap is NON-VACUOUS: the counted method is THE
    apply-direction leaf kernel (vv Mode 11).

    Two claims, and the second is the one that makes the first mean
    something:

    1. ``residual_kernel_batch`` fires (> 0) on a real ``A.apply(x)``, on
       BOTH production geometries.  A counter that stays 0 because the
       production path routes around the wrap makes every count assertion
       below vacuous — a green gate that measures nothing.
    2. The three plausible SIBLING kernels stay at exactly 0:
       ``cell_kernel_batch`` (the SOLVE-direction twin), and the per-cell
       ``residual`` / ``update`` reference contracts.  Without this, "the
       wrap fired" would be consistent with the apply routing through a
       different kernel that merely also touches this one.

    This test is what would have caught the spec's stale
    ``update_batch`` name: that method does not exist on the tree, so a
    gate written against it would have raised ``AttributeError`` — or,
    worse, if a same-named method had survived elsewhere, counted 0
    forever.
    """
    for label, sn_mesh in (
        ("slab", _slab(_LADDER_SLAB_NX[0])),
        ("cart2d", _cart2d(*_LADDER_ISOTROPIC[0])),
    ):
        fired = _leaf_calls_for_apply(sn_mesh)
        if fired <= 0:
            pytest.fail(
                f"{label}: the P-1 wrap on "
                f"DiamondDifference.residual_kernel_batch NEVER FIRED during "
                f"A.apply(x) (count={fired}).  The gate is VACUOUS — the "
                f"production apply routes around the wrapped method.  Find "
                f"the method the path really enters and wrap THAT "
                f"(vv Mode 11)."
            )
        for sibling in ("cell_kernel_batch", "residual", "update"):
            cold = _leaf_calls_for_apply(sn_mesh, method=sibling)
            if cold != 0:
                pytest.fail(
                    f"{label}: sibling kernel {sibling!r} fired {cold} times "
                    f"during A.apply(x).  The apply direction is supposed to "
                    f"enter residual_kernel_batch and nothing else; if this "
                    f"changed, P-1's counted method is no longer THE leaf "
                    f"and every rung below must be re-derived."
                )


def test_p1_scanned_axis_call_count_is_mesh_invariant() -> None:
    r"""[foundation] **THE CATCHER.** 2-D Cartesian: refining the SCANNED
    axis must not change the leaf-kernel call count, at all.

    ``ScanMarch`` (the 2-D production default) is ``scan(x) o march(y)``:
    the whole x-axis is inside one Blelloch ``ordinate_scan`` per row, so
    the leaf kernel is entered once per (octant, y-row) and ``nx`` is
    invisible to the arity.  Refining ``nx`` 8 -> 64 multiplies the cell
    count by 8 and MUST leave the count bit-identical.

    This is the L16 regression in its exact shape: the moment a block
    composition (``A_ij = R_i A J_j``) pulls the x-axis out of the scan
    and into a per-cell Python fold, this count scales 8x with the mesh
    and the gate reds — deterministically, with no timing noise.

    TEETH (plan §7 M-20), mutation-verified in process: re-dispatching
    the batched call once per column of the level slice takes the counts
    from ``{8: 80, 16: 80, 32: 80, 64: 80}`` to
    ``{8: 720, 16: 1360, 32: 2640, 64: 5200}`` -> RED.  **That mutation
    is EXACTLY value-identical** (``np.array_equal`` True, ``max|diff|``
    0.0e+00 — the DD residual kernel is cell-local, so a column-wise
    split changes nothing but the arity).  That measurement is the whole
    argument for this file: the L16 regression class cannot be seen by
    any value gate, at any tolerance, ever.  Only the count moves.
    """
    counts = {
        nx: _leaf_calls_for_apply(_cart2d(nx, _LADDER_NY))
        for nx in _LADDER_NX
    }
    expected = _CART2D_CALLS_PER_ROW * _LADDER_NY
    baseline = counts[_LADDER_NX[0]]
    np.testing.assert_array_equal(
        np.asarray([counts[nx] for nx in _LADDER_NX]),
        np.full(len(_LADDER_NX), expected),
        err_msg=(
            f"leaf-kernel call count is NOT invariant under scanned-axis "
            f"(nx) refinement: {counts} (expected {expected} on every rung, "
            f"= _CART2D_CALLS_PER_ROW({_CART2D_CALLS_PER_ROW}) * "
            f"ny({_LADDER_NY})).  Growth with nx means a hoisted vectorised "
            f"call has been replaced by a per-cell Python fold — the L16 "
            f"regression.  Ratio nx=64/nx=8: "
            f"{counts[_LADDER_NX[-1]] / max(baseline, 1):.2f}x."
        ),
    )


def test_p1_call_count_is_invariant_in_ordinates_and_groups() -> None:
    r"""[foundation] The ORDINATE and GROUP axes are internal to every
    kernel invocation — refining either must not change the arity.

    ``sweep_graph``'s design statement: "the anti-diagonal schedule
    vectorises ALL THREE axes inside one cell-kernel call ... NO Python
    loop over ordinates, cells-along-a-diagonal, or groups."  This gate
    pins the ordinate and group halves of that claim (the cell half is
    P-1a).  S2 -> S6 triples the ordinate count; ng 1 -> 4 quadruples the
    group count; the count must not move.

    A block composition that slices per group (a plausible ``R_i`` /
    ``J_j`` indexing choice) reds the ng leg by 4x while every correctness
    gate stays green.
    """
    per_order = {
        order: _leaf_calls_for_apply(_cart2d(8, _LADDER_NY, sn_order=order))
        for order in (2, 4, 6)
    }
    np.testing.assert_array_equal(
        np.asarray(list(per_order.values())),
        np.full(len(per_order), _CART2D_CALLS_PER_ROW * _LADDER_NY),
        err_msg=(
            f"call count moved with the quadrature order: {per_order} — the "
            f"ordinate axis has left the batched kernel and become a Python "
            f"loop."
        ),
    )
    per_groups = {
        ng: _leaf_calls_for_apply(_cart2d(8, _LADDER_NY, ng=ng))
        for ng in (1, 2, 4)
    }
    np.testing.assert_array_equal(
        np.asarray(list(per_groups.values())),
        np.full(len(per_groups), _CART2D_CALLS_PER_ROW * _LADDER_NY),
        err_msg=(
            f"call count moved with the group count: {per_groups} — the "
            f"group axis has left the batched kernel and become a Python "
            f"loop."
        ),
    )


def test_p1_calls_per_cell_decays_under_isotropic_refinement() -> None:
    r"""[foundation] Under isotropic refinement the cell count grows 4x per
    rung while the call count grows at most 2x — so calls-per-cell DECAYS.

    The complement of P-1a: it admits the honest ``ny`` linearity (the
    marched axis genuinely costs one kernel entry per row) while still
    forbidding the per-cell regime, in which calls-per-cell would be
    CONSTANT.  Stated as a strict inequality on the ratio, so it survives
    a future schedule change that alters the constant without reverting
    to a fold.
    """
    rows = []
    for nx, ny in _LADDER_ISOTROPIC:
        calls = _leaf_calls_for_apply(_cart2d(nx, ny))
        rows.append((nx, ny, nx * ny, calls, calls / (nx * ny)))
    table = "  ".join(
        f"({nx}x{ny}: {calls} calls, {per:.4f}/cell)"
        for nx, ny, _, calls, per in rows
    )
    for (nx0, ny0, cells0, calls0, per0), (nx1, ny1, cells1, calls1, per1) in zip(
        rows, rows[1:],
    ):
        if not cells1 == 4 * cells0:
            pytest.fail(
                f"fixture bug: isotropic refinement must quadruple the cell "
                f"count; got {cells0} -> {cells1}"
            )
        if not per1 < 0.75 * per0:
            pytest.fail(
                f"calls-per-cell did not decay under isotropic refinement "
                f"({nx0}x{ny0} -> {nx1}x{ny1}): {per0:.4f} -> {per1:.4f} "
                f"/cell (call count {calls0} -> {calls1} while cells "
                f"{cells0} -> {cells1}).  A constant calls-per-cell IS the "
                f"per-cell Python fold (L16).  Full table: {table}"
            )


def test_p1_slab_apply_is_a_per_cell_fold_regression_floor() -> None:
    r"""[foundation] **FINDING, pinned as a floor — NOT an invariance claim.**

    Tracked as `#323 <https://github.com/deOliveira-R/ORPHEUS/issues/323>`_
    (measurement-gated, like #227 — 1-D meshes are small, so the carve may
    not be worth its cost until a profile says so).  **Whoever closes #323
    will get a RED from this gate.  That is intended** — see the two-sided
    reddening below.

    The 1-D *apply* direction is ALREADY a per-cell Python fold on this
    tree: ``_OneDimScanWalk``'s apply-loop frame visits every cell of
    every walk-DAG leg and calls the batched leaf once per visit, so
    ``calls = 2 * nx`` exactly (2 = the inward + outward legs; measured
    independent of the quadrature order and the group count).  The
    *solve* direction on the SAME mesh is fully hoisted — see the control
    leg :func:`test_p1_slab_solve_call_count_is_mesh_invariant`.

    This gate therefore pins the per-cell ARITY at exactly 2, which
    reddens on BOTH sides deliberately:

    * arity > 2 — a further degradation (a per-ordinate or per-group
      inner loop inside the visit) that would be invisible to every
      correctness gate;
    * arity < 2 — an IMPROVEMENT (hoisting the 1-D apply onto the scan,
      as the solve already is).  That is a welcome change and it must
      still red here, so the baseline is re-captured consciously rather
      than drifting (lesson L7: the tolerance IS the claim).

    TEETH, mutation-verified: a per-ORDINATE fold inside the per-cell
    visit takes the counts from ``{20: 40, 40: 80, 80: 160, 160: 320}``
    to ``{20: 200, 40: 400, 80: 800, 160: 1600}`` -> RED.

    HONEST SCOPE: the *cell*-axis fold that reddens P-1a is INERT here,
    and that is the finding rather than a gap.  MEASURED: the slab hands
    the leaf ``psi_bar`` of shape ``(4, 2, 1)`` — the cell axis is
    already width 1 — where the 2-D row-march hands it ``(3, 2, nx)``.
    There is no batched cell axis left in 1-D to fold; the axes that
    remain batched (ordinate, group) are the ones this leg can still
    lose, and it catches those.
    """
    counts = {nx: _leaf_calls_for_apply(_slab(nx)) for nx in _LADDER_SLAB_NX}
    np.testing.assert_array_equal(
        np.asarray([counts[nx] for nx in _LADDER_SLAB_NX]),
        np.asarray([_SLAB_CALLS_PER_CELL * nx for nx in _LADDER_SLAB_NX]),
        err_msg=(
            f"slab apply arity moved off the pinned {_SLAB_CALLS_PER_CELL} "
            f"leaf-kernel calls per cell: {counts}.  If the count DROPPED, "
            f"the 1-D apply was hoisted onto the scan (good — re-capture "
            f"_SLAB_CALLS_PER_CELL and say so in the closeout).  If it "
            f"ROSE, a per-ordinate or per-group inner loop has appeared "
            f"inside the per-cell visit (L16)."
        ),
    )


def test_p1_slab_solve_call_count_is_mesh_invariant() -> None:
    r"""[foundation] **CONTROL LEG** — the counting instrument CAN see
    mesh-invariance, so the finding above is a property of the code.

    On the very same slab meshes, the SOLVE direction's Blelloch scan is
    entered exactly ``_SLAB_SOLVE_SCAN_CALLS`` times regardless of ``nx``
    ("Exactly 2 scan calls per sweep" — one per chain).  Without this
    leg, a reader of the per-cell-fold finding could not distinguish "the
    1-D apply is folded" from "this harness cannot observe invariance in
    1-D" (lesson L6: a finding needs its control).

    It is also a live tripwire on the hoisted half: if the 1-D SOLVE ever
    grows a per-cell loop, this reds while every correctness gate stays
    green — which is precisely the L16 event.

    TEETH, mutation-verified: re-dispatching the chain scan once per cell
    takes the counts from ``{20: 2, 40: 2, 80: 2, 160: 2}`` to
    ``{20: 42, 40: 82, 80: 162, 160: 322}`` -> RED.
    """
    counts: dict[int, int] = {}
    for nx in _LADDER_SLAB_NX:
        record = record_for(_slab(nx), scattering_order=1)
        implicit = record.implicit_operator
        if not isinstance(implicit, StreamingCollisionOperator):
            pytest.fail(
                f"fixture bug: a seedless slab must carry the bare (L+C) as "
                f"its implicit operator; got {type(implicit).__name__}.  The "
                f"scan-count claim below is about the 1-D sweep, not a "
                f"coupled back-substitution."
            )
        # ``system_a`` is the package's ONE narrowing of a coupled member to
        # its real composite type — the bare ``(L+C)``'s domain is the
        # System-A ``FullField``, not the coupled carrier.
        rhs = system_a(random_state(record, seed=5))
        implicit.solve(rhs)                          # warm the caches
        with _count_calls(_loss_representation, "ordinate_scan") as counter:
            implicit.solve(rhs)
        counts[nx] = counter.n
    np.testing.assert_array_equal(
        np.asarray(list(counts.values())),
        np.full(len(counts), _SLAB_SOLVE_SCAN_CALLS),
        err_msg=(
            f"the 1-D SOLVE direction's scan count is no longer "
            f"mesh-invariant: {counts} (expected "
            f"{_SLAB_SOLVE_SCAN_CALLS} on every rung).  The hoisted "
            f"vectorised axis has been folded back into a per-cell loop — "
            f"the exact L16 regression, on the one path that was still "
            f"hoisted."
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# P-2 — normalised throughput.  NEVER an absolute ms threshold.
# ═══════════════════════════════════════════════════════════════════════


def _best_of(fn, n_reps: int) -> float:
    """Minimum wall time over ``n_reps`` calls, after one warm-up.

    ``min``, never the mean: the minimum is the only order statistic not
    contaminated by a descheduling event, and this host is not quiescent.
    """
    fn()
    best = float("inf")
    for _ in range(n_reps):
        t0 = time.perf_counter()
        fn()
        best = min(best, time.perf_counter() - t0)
    return best


def test_p2_composition_overhead_vs_in_process_calibration() -> None:
    r"""[foundation] One ``A.apply(x)`` costs at most
    ``_COMPOSITION_OVERHEAD_MAX`` x a dense contraction of the same
    nominal FLOP count, measured IN THE SAME PROCESS.

    The claim, stated plainly: *the block composition costs at most N x
    the raw contraction of the same FLOP count on this host.*  The
    in-process calibration is what makes it machine-independent — the
    existing ``tests/sn/sweep/core/test_cache.py:553`` gate
    (``elapsed_per_sweep_ms < 2.0``) is both the precedent for having a
    throughput gate at all and the cautionary tale for how to write one:
    an absolute millisecond threshold encodes the machine it was captured
    on, and either false-reds on a slow host or goes permanently green on
    a fast one.

    PROVISIONAL: the pre-carve band (23.4 - 24.4) was captured on a
    CONTENDED host; ``_COMPOSITION_OVERHEAD_MAX`` is set at ~3x that so a
    contended run cannot false-red.  Re-measure on a quiescent tree and
    tighten.  If this is the ONLY red leg, re-run it alone before
    believing it — P-1 and P-3 are the contention-immune legs.
    """
    A, x = _posed_apply(_cart2d(_COST_NX, _COST_NY))
    n_dof = _n_dof(x)
    A.apply(x)

    side = int(round(np.sqrt(_CALIBRATION_FLOPS_PER_DOF * n_dof)))
    rng = np.random.default_rng(20260728)
    dense = rng.standard_normal((side, side))
    vector = rng.standard_normal(side)

    t_apply = _best_of(lambda: A.apply(x), _N_REPS_APPLY)
    t_calibration = _best_of(
        lambda: np.einsum("ij,j->i", dense, vector), _N_REPS_CALIBRATION,
    )
    if not t_calibration > 0.0:
        pytest.fail("calibration timed at zero — the clock is unusable here")
    ratio = t_apply / t_calibration
    print(
        f"\nP-2: n_dof={n_dof} calib={side}x{side} "
        f"apply={t_apply * 1e3:.3f} ms calib={t_calibration * 1e6:.1f} us "
        f"ratio={ratio:.2f} (max {_COMPOSITION_OVERHEAD_MAX})"
    )
    if not ratio < _COMPOSITION_OVERHEAD_MAX:
        pytest.fail(
            f"composition overhead {ratio:.2f}x exceeds the "
            f"{_COMPOSITION_OVERHEAD_MAX}x claim (pre-carve band 23.4-24.4). "
            f"apply={t_apply * 1e3:.3f} ms against a {side}x{side} dense "
            f"contraction at {t_calibration * 1e6:.1f} us.  Check P-1 first: "
            f"if a call count also scaled, this is the L16 fold and the "
            f"timing is corroboration, not the finding.  If P-1 is green, "
            f"re-run this leg ALONE on a quiescent host before believing it."
        )


# ═══════════════════════════════════════════════════════════════════════
# P-3 — allocation ceiling.  Exactly reproducible; contention-immune.
# ═══════════════════════════════════════════════════════════════════════


def test_p3_peak_allocation_per_dof() -> None:
    r"""[foundation] Peak bytes allocated by ONE ``A.apply(x)``, per DOF-word.

    ``tracemalloc`` peak / ``(n_dof * 8)``.  1.0 is the structural floor
    (the output field itself); the pre-carve ~3.03 is about three
    full-field temporaries.  This is the leg that catches "materialise a
    slice per block per cell" DIRECTLY — the mechanism stage 3's
    ``A_ij = R_i A J_j`` composition makes newly expressible — and unlike
    wall clock it is exactly reproducible (measured spread across six
    trials: 0.03 %), so it is contention-immune and its baseline is final.

    TEETH, mutation-verified: holding ONE extra full-field temporary
    alive across the apply reads 4.035 -> RED (two: 5.033; three: 6.033).
    Sensitivity is therefore exactly one field.  Note the metric is a
    PEAK, so a short-lived transient smaller than the existing peak is
    invisible by construction (a 0.6 MB scratch per call moved the
    reading by 0.0006); the mechanism it does see — temporaries that stay
    alive across the composition — is the one stage 3 introduces.
    """
    A, x = _posed_apply(_cart2d(_COST_NX, _COST_NY))
    n_dof = _n_dof(x)
    A.apply(x)                                   # warm every derived cache

    tracemalloc.start()
    try:
        tracemalloc.reset_peak()
        A.apply(x)
        peak = tracemalloc.get_traced_memory()[1]
    finally:
        tracemalloc.stop()

    factor = peak / (n_dof * 8)
    print(
        f"\nP-3: peak={peak} B  n_dof={n_dof}  "
        f"peak/(n_dof*8)={factor:.4f} (max {_ALLOC_FACTOR_MAX})"
    )
    if not factor < _ALLOC_FACTOR_MAX:
        pytest.fail(
            f"peak allocation {factor:.4f} x (n_dof * 8) exceeds the "
            f"{_ALLOC_FACTOR_MAX} x ceiling (pre-carve 3.03).  peak={peak} B "
            f"for n_dof={n_dof}.  The composition is materialising "
            f"per-block temporaries; this leg is exactly reproducible, so "
            f"it is a real regression, not noise."
        )


# ═══════════════════════════════════════════════════════════════════════
# P-4 — the baseline's integrity.  A constant measured against a fixture
# that has since moved is worse than no constant.
# ═══════════════════════════════════════════════════════════════════════


def test_p4_fixture_fingerprint_matches_the_baseline() -> None:
    r"""[foundation] The cost fixture still has the shape every committed
    constant was measured against.

    A performance baseline is a pair (number, fixture).  If the fixture
    drifts — a rung resized, the quadrature order bumped, a region
    dropped — the committed constants silently start describing a
    different problem and both P-2 and P-3 become unfalsifiable in the
    direction that matters.  This gate makes that drift LOUD, with the
    remedy in the message.
    """
    sn_mesh = _cart2d(_COST_NX, _COST_NY)
    _, x = _posed_apply(sn_mesh)
    measured = {
        "n_cells": int(np.prod(sn_mesh.spatial_shape)),
        "n_ordinates": int(sn_mesh.quad.N),
        "n_groups": int(sn_mesh.ng),
        # Materials as the OPERATORS resolve them (the typed per-material
        # index map), not as the raw ``mat_map`` declares them — a region
        # that no cell carries is not heterogeneity, and this accessor is
        # mesh-dimension-agnostic where ``mat_map`` is 2-D-only.
        "n_regions": len(sn_mesh.material_xs_field().cells_by_material),
        "n_dof": _n_dof(x),
    }
    if measured != _COST_FINGERPRINT:
        pytest.fail(
            f"the P-2/P-3 cost fixture has drifted from the baseline it was "
            f"measured against.\n  expected {_COST_FINGERPRINT}\n  measured "
            f"{measured}\nRE-CAPTURE the baseline (P-4): re-measure "
            f"_COMPOSITION_OVERHEAD_MAX and _ALLOC_FACTOR_MAX against the "
            f"new fixture on a QUIESCENT host and update the provenance "
            f"banner, or restore the fixture."
        )
