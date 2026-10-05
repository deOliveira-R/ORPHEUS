r"""L14 four-leg standoff — slab + cylinder via the typed operator algebra.

Per ``.claude/lessons.md`` L14, solver correctness is four-legged:

1.  **Algorithm 1 (Krylov via ``(L + C)``) ≡ structurally-independent reference.**
2.  **Algorithm 2 (sweep via ``sweep_once`` — the production strategy path) ≡ structurally-independent reference.**
3.  **Algorithm 1 ≡ Algorithm 2** (twin-path agreement).
4.  **All three under mesh refinement** (right rate to right limit).

Two algorithms agreeing is necessary but NOT sufficient — both can be
equally wrong.  Post-D-K (commit ``dadf4e8``), ``solve_sn`` routes
through ``StreamingOperator + CollisionOperator`` =
:class:`StreamingCollisionOperator`; the Krylov path uses GMRES on
``StreamingCollisionOperator.apply``, left-preconditioned by the sweep since
#200 (2026-10-04; an identity before, when its inner iteration count grew with
the mesh). The two paths' FIXED POINTS share no cell kernel: the Krylov
matvec's per-cell algebra is ``DiamondDifference.residual_kernel_batch``, the
sweep's is ``affine_scan_coefficients``; the preconditioner reads the sweep's
but moves only GMRES's trajectory. So a defect in one kernel moves one path's
k only (a matvec defect also trips the sweep's convergence-claim guard, which
re-measures the sweep's residual through the matvec).

References (semi-analytical pillar per ``vv-principles``):

- **Slab** — the Case singular-eigenfunction reference
  ``sn_slab_1eg_2rg_S8`` (Case & Zweifel 1967; Gülderen & Türeci 2023;
  Garis & Sjöstrand 1990).  1G two-region reflective, S8.  The
  reference is mesh-independent and mathematically self-contained;
  ``solve_sn`` at the matching quadrature order must converge to it
  as h → 0.
- **Cylinder** — the trajectory-resolvent Variant α Green's-function
  reference ``solve_greens_function_cylinder_mr``.  3-region 2G
  reflective ABA layout.  Uses ray-traced chord integration in 3-D
  space; structurally independent of every discrete primitive in the
  SN code (no shared FP path, no shared redist closure, no shared
  boundary recurrence).

Sweep route (``inner_solver="source_iteration"``) routes through the
``(L+C)`` strategy sweep / ``DiscretizationScheme.update``, which uses WDD
(Cartesian) / the same per-cell algebra as
:func:`transport_operator_matvec_unified` (curvilinear).  (The operator-free
``transport_sweep`` entry retired at the coupled-block campaign step 6.)
"""
from __future__ import annotations

import contextlib
import functools
import math

import pytest

from orpheus.derivations.reference_values import continuous_get
from orpheus.geometry import BC, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh1D, Mesher
from orpheus.sn import solve_sn
from orpheus.numerics.quadrature import Quadrature
from tests.gates.derivations._trajectory_resolvent_ladders import (
    CYLINDER_3REG_SN_4X8_K_STEP,
    tolerance_for,
)
from orpheus.geometry import CoordSystem
from orpheus.numerics.observable import Eigenvalue
from tests.gates.sn.verification.analytical._aba_reference import (
    aba_materials,
    aba_reference,
    aba_uniform_width_mesh,
    assert_cylinder_record,
    awaits_cylinder_bound,
    verify_cylinder_k,
)


# Post-D-K (commit ``dadf4e8``), the within-group loss composite
# ``L + C`` (:func:`build_streaming_collision` →
# :class:`StreamingCollisionOperator` = ``StreamingOperator +
# MultiplicationOperator``) routes through
# :func:`transport_operator_matvec_unified` natively for 1-D slab /
# sphere / cylinder and through
# the representation's ``loss_action`` (which since S6.3 lives on
# the loss representation, off the operator; ``ScanMarch`` default
# since S6.9) for 2-D Cartesian.
# No monkey-patch is required.


# ═══════════════════════════════════════════════════════════════════════
# Cylinder fixtures
# ═══════════════════════════════════════════════════════════════════════
def _cylinder_k_ref() -> float:
    """The shared cylinder reference's eigenvalue, uncertified (solved once per session; the problem is ``_aba_reference``'s)."""
    return aba_reference(CoordSystem.CYLINDRICAL).read(Eigenvalue()).value


def _k(solution) -> float:
    """Production's eigenvalue, through the reading verb."""
    return solution.read(Eigenvalue()).value


@functools.cache
def _solve_cyl_via_krylov_unified(nx: int):
    mesh, materials = aba_uniform_width_mesh(CoordSystem.CYLINDRICAL, nx), dict(aba_materials())
    quad = Quadrature.folded_product(n_mu=4, n_phi=8)
    sol = solve_sn(
            materials=materials, mesh=mesh, quadrature=quad,
            inner_solver="krylov",
            max_outer=200, keff_tol=1e-7, flux_tol=1e-7,
            max_inner=200, inner_tol=1e-9,
        )
    return sol


@functools.cache
def _solve_cyl_via_sweep(nx: int):
    mesh, materials = aba_uniform_width_mesh(CoordSystem.CYLINDRICAL, nx), dict(aba_materials())
    quad = Quadrature.folded_product(n_mu=4, n_phi=8)
    sol = solve_sn(
        materials=materials, mesh=mesh, quadrature=quad,
        inner_solver="source_iteration",
        max_outer=500, keff_tol=1e-7, flux_tol=1e-7,
        max_inner=500, inner_tol=1e-9,
    )
    return sol


# ═══════════════════════════════════════════════════════════════════════
# Cylinder L1 standoff
# ═══════════════════════════════════════════════════════════════════════
# Reference tolerance, once the reference is certified: ``tolerance_for(e, None)``
# of ``tests/gates/derivations/_trajectory_resolvent_ladders.py``, the
# reference assumed at the floor, for the folded-4x8 solve's own k error
# e = 5.9e-4 ([M] 2026-09-26), giving 2e-3. It was 3e-2, justified as
# "trajectory_resolvent's quadrature error budget at n_r=24"; that budget was
# the reference's one-spline emission density (ERR-090), and the reference
# still carries no certificate (#566, #516), so the reference legs are strict
# xfails on the verbs' refusal, with a RECORD companion.
# Twin-path tolerance: 1e-5 rel (both algorithms converge to the same
# discrete fixed point at the matched solver tolerances).
_CYL_REF_RTOL = tolerance_for(CYLINDER_3REG_SN_4X8_K_STEP, None)
_CYL_TWIN_RTOL = 1.0e-5
_CYL_SUPPORTS = (
    "tests/gates/derivations/test_trajectory_resolvent_regionwise_source.py"
    "::test_mr_oracle_first_leg_matches_the_line_integral[cylinder]",
)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.rests_on(*_CYL_SUPPORTS)
@awaits_cylinder_bound
def test_cylinder_l1_sweep_vs_trajectory_resolvent() -> None:
    r"""**Cylinder Leg 2** — sweep ≡ trajectory_resolvent reference.

    Production source-iteration path (the ``(L+C)`` strategy sweep /
    ``DiscretizationScheme.update``) on the 3-region 2G ABA cylinder. No shim:
    the sweep already uses the WDD-correct per-cell algebra that the
    unified matvec also wraps. Strict ``xfail`` on the verbs' refusal: the
    reference's family derives no bound (#566, #516), so it has no
    certificate; ``test_cylinder_l1_reference_record`` keeps the reading live.
    """
    verify_cylinder_k(_solve_cyl_via_sweep(nx=40), _CYL_REF_RTOL)


@pytest.mark.l1
@pytest.mark.slow
def test_cylinder_l1_sweep_vs_krylov_twin_path() -> None:
    r"""**Cylinder Leg 3** — sweep ≡ Krylov-via-unified twin-path agreement.

    Both algorithms drive the SAME continuous-equation discrete fixed
    point; at matched solver tolerances they MUST agree at sub-ULP-of-
    iteration drift. Disagreement here is the L14 manifestation-#6
    signature: same equation, two algorithmic paths, two answers.

    HISTORY — un-xfailed at S6.4 (2026-06-11): the pre-D-K divergence
    (rel ≈ 4e-3 at nx=40, the cell-centre-proxy Carlson seed) was
    HEALED when D-K retargeted ``solve_sn`` onto the B1''-aware (L+C)
    algebra, exactly as the original xfail reason predicted ("should be
    re-validated and likely flipped green").  The strict xfail did its
    job: the first full slow-suite run after the heal reported
    XPASS(strict), and two independent executions (assertions live
    under ``-O`` — pytest rewrites test-module asserts) confirmed
    sweep ≡ Krylov ≡ trajectory reference at every leg.
    """
    k_sweep = _k(_solve_cyl_via_sweep(nx=40))
    k_krylov = _k(_solve_cyl_via_krylov_unified(nx=40))
    rel = abs(k_sweep - k_krylov) / k_sweep
    assert rel < _CYL_TWIN_RTOL, (
        f"cylinder twin-path disagreement: "
        f"k_sweep={k_sweep:.10f}, k_krylov={k_krylov:.10f}, rel={rel:.3e}"
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.parametrize("nx", [20, 40])
def test_cylinder_l1_refinement_both_paths(nx: int) -> None:
    r"""**Cylinder Leg 4, twin half** — sweep ≡ Krylov-via-unified at each refinement.

    For each nx ∈ {20, 40} the sweep and the Krylov-via-unified paths
    must reach the same discrete fixed point (:data:`_CYL_TWIN_RTOL`). The
    reference half of this leg is
    ``test_cylinder_l1_refinement_against_reference``, split off on
    2026-09-26 when the reference legs became strict xfails on #516, so that
    this assertion stays live.

    WHY NOT nx = 80 (dropped 2026-10-04): the twin claim does not depend on
    the mesh, and nx = 80 was its least discriminating member at its highest
    cost. ``[M]`` 2026-10-04 (``scratch/reference_architecture/p3/perf_l1/
    p2_cyl_timing.log``): the twin gap reads 4.2e-11, 4.7e-11 and 5.1e-11
    relative at nx = 20, 40, 80, six orders below the band at every level, so
    the third level adds no reading the first two lack; a defect local to one
    cell (the axis cell, an interface cell) moves k in proportion to that
    cell's weight, which is largest on the coarsest mesh; and the nx = 80
    Krylov solve cost 184 s of the row's 192 s, unpreconditioned (before #200). The
    refinement ladder that needs three levels is the reference half's "right
    rate to right limit", which keeps nx = 80.

    HISTORY — un-xfailed at S6.4 (2026-06-11) alongside
    ``test_cylinder_l1_sweep_vs_krylov_twin_path`` (same heal: the
    D-K retargeting of ``solve_sn`` onto the B1''-aware (L+C) algebra
    removed the pre-D-K Carlson cell-centre-proxy divergence; validated
    by two independent XPASS(strict) executions at all three nx).
    """
    k_sweep = _k(_solve_cyl_via_sweep(nx=nx))
    k_krylov = _k(_solve_cyl_via_krylov_unified(nx=nx))
    rel_twin = abs(k_sweep - k_krylov) / k_sweep
    assert rel_twin < _CYL_TWIN_RTOL, (
        f"cylinder nx={nx}: twin-path rel={rel_twin:.3e} ≥ {_CYL_TWIN_RTOL:.0e}"
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.rests_on(*_CYL_SUPPORTS)
@awaits_cylinder_bound
@pytest.mark.parametrize("nx", [20, 40, 80])
def test_cylinder_l1_refinement_against_reference(nx: int) -> None:
    r"""**Cylinder Leg 4, reference half** — both paths against the reference at each refinement.

    "Right rate to right limit": at nx ∈ {20, 40, 80} both algorithms must
    agree with the trajectory_resolvent reference to :data:`_CYL_REF_RTOL`.
    Strict ``xfail`` on the verbs' refusal: the reference has no certificate
    (#566, #516).

    Each path is solved when its comparison is reached, not before: while the
    reference is uncertified the first comparison refuses, so the Krylov solve
    (184 s at nx = 80 before #200's preconditioner) is never paid for a refusal;
    once the family is certified every solve runs, as before.
    """
    for solve in (_solve_cyl_via_sweep, _solve_cyl_via_krylov_unified):
        verify_cylinder_k(solve(nx=nx), _CYL_REF_RTOL)


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.rests_on(*_CYL_SUPPORTS)
def test_cylinder_l1_reference_record() -> None:
    r"""RECORD: the reference's k and the sweep's k (nx = 40), as they read today.

    Not verification: it keeps the cylinder reference legs live while they
    are strict xfails, and reddens when either side moves (an SN change, or
    the reference's #516/#566 repair, after which the family is certified and
    the xfails lifted).
    """
    assert_cylinder_record({"k_ref": _cylinder_k_ref(), "standoff_sweep_k_nx40": _k(_solve_cyl_via_sweep(nx=40))})


# ═══════════════════════════════════════════════════════════════════════
# Slab fixtures
# ═══════════════════════════════════════════════════════════════════════


def _build_slab_2region_mesh(n_per: int) -> tuple[Mesh1D, dict, int]:
    r"""Build the ``sn_slab_1eg_2rg_S8`` Case singular-eigenfunction mesh.

    Returns ``(mesh, materials, N_ord)`` for direct use by ``solve_sn``.
    """
    ref = continuous_get("sn_slab_1eg_2rg_S8")
    geom = ref.problem.geometry_params
    materials = ref.problem.materials
    H_A = float(geom["fuel_height"])
    H_B = float(geom["refl_height"])
    N_ord = int(geom["n_ordinates"])
    slab = StructuredGeometry.slab(
        (0.0, H_A, H_A + H_B), (0, 1), left=BC.reflective, right=BC.reflective,
    )
    mesh = Mesher(slab).partition(CellsByCount.uniform_width(n_per)).mesh
    return mesh, materials, N_ord


def _slab_k_ref() -> float:
    k_eff = continuous_get("sn_slab_1eg_2rg_S8").k_eff
    if k_eff is None:
        raise ValueError("sn_slab_1eg_2rg_S8 carries no eigenvalue: the slab reference answers the k question")
    return float(k_eff)


#: The inner tolerance of each inner solver on the slab: the Krylov solve (GMRES
#: on the unified matvec) and the source iteration (the ``(L+C)`` sweep), the
#: ones the rows always used; the outer tolerance is ``keff_tol = 1e-12`` for both.
_SLAB_INNER_TOL = {"krylov": 1e-9, "source_iteration": 1e-12}
_SLAB_PATHS = tuple(_SLAB_INNER_TOL)


@functools.cache
def _slab_k(inner_solver: str, n_per: int) -> float:
    """Production's k on the Case problem with ``n_per`` cells per region, by ``inner_solver``.

    Cached as the cylinder solves are: each (path, mesh) is solved once per
    session, and every row reads the same solves, so no problem is solved twice.
    """
    mesh, materials, N_ord = _build_slab_2region_mesh(n_per=n_per)
    sol = solve_sn(
        materials, mesh, Quadrature.gauss_legendre(N_ord),
        inner_solver=inner_solver, max_outer=500, max_inner=500,
        keff_tol=1e-12, inner_tol=_SLAB_INNER_TOL[inner_solver],
    )
    return float(sol.outcome.keff)


# ═══════════════════════════════════════════════════════════════════════
# Slab L1 standoff
# ═══════════════════════════════════════════════════════════════════════
# The ladder: 10, 20 and 40 cells per region. Every row reads these six solves
# (two paths, three meshes) and no others.
#
# ``[M]`` 2026-10-04 (``scratch/reference_architecture/p3/perf_l1/`` p4, p5, p7,
# and ``scratch/reference_architecture/p3/gates_repair/battery/l1_*.log``),
# |k − k_Case| on the ladder; production reads the same on both paths:
#
#   ===================================  ========  ========  ========
#   mutation (path it moves)             n_per=10  n_per=20  n_per=40
#   ===================================  ========  ========  ========
#   none: diamond difference             3.45e-05  8.63e-06  2.16e-06
#   step scheme (Krylov matvec)          2.47e-04  1.04e-04  4.66e-05
#   Σ_t × (1 + 1e-3) (Krylov matvec)     4.15e-04  4.41e-04  4.48e-04
#   ERR-025's a (sweep scan)             6.87e-02  7.11e-02  7.23e-02
#   ===================================  ========  ========  ========
#
# Production's observed order is 2.00 on both halvings; the step scheme's is
# 1.25 and 1.16, and the wrong-limit mutations have none. The rows these
# replace sat at n_per = 160 with a 2e-5 band, which the step-scheme matvec
# PASSED (1.06e-5 at 160): an O(h) scheme reaches any band on a fine enough
# mesh, so a band at one fine mesh cannot tell first order from second.
#
# A matvec mutation also reddens the sweep rows, through production's own
# guard rather than their assertions: the source iteration re-measures its
# claimed convergence with the honest residual ‖Aψ − q‖, which reads the
# matvec, and raises ``ConvergenceClaimError``. ERR-025's coefficient in the
# scan does NOT trip that guard; it reddens the sweep rows on their assertions.
_SLAB_LADDER = (10, 20, 40)
#: The order every path must show on each halving of the ladder: diamond
#: difference is second order (2.00 measured), the step scheme first (1.25).
_SLAB_MIN_ORDER = 1.8
#: |k − k_Case| at n_per = 40: about nine times production's 2.16e-6, and below
#: the step-scheme matvec's 4.66e-5.
_SLAB_REF_ABSTOL_AT_40 = 2.0e-5
#: Twin-path tolerance, relative: both paths converge to one discrete fixed point.
_SLAB_TWIN_RTOL = 1.0e-5

_HERE = "tests/gates/sn/verification/analytical/test_l1_standoff_slab_cylinder.py"


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.parametrize("inner_solver", _SLAB_PATHS)
def test_slab_l1_order_against_case(inner_solver: str) -> None:
    r"""**Slab Leg 4** — each path converges to the Case reference at second order.

    The observed order :math:`p = \log_2(e_h / e_{h/2})` of the error against
    the Case reference, on each halving of the ladder 10 → 20 → 40 cells per
    region, is at least :data:`_SLAB_MIN_ORDER` (measured 2.00 on both paths:
    diamond difference keeps its second order across the material interface,
    which falls on a cell edge). A consistent first-order scheme (the step
    scheme, 1.25 and 1.16) and a wrong limit (the error stops falling) both
    redden here at every mesh size; neither can hide behind a finer mesh.
    """
    k_ref = _slab_k_ref()
    errors = [abs(_slab_k(inner_solver, n) - k_ref) for n in _SLAB_LADDER]
    orders = [math.log2(coarse / fine) for coarse, fine in zip(errors, errors[1:])]
    assert min(orders) >= _SLAB_MIN_ORDER, (
        f"slab {inner_solver}: observed orders {[f'{p:.3f}' for p in orders]} on n_per {_SLAB_LADDER} "
        f"(|Δk| {[f'{e:.3e}' for e in errors]}), below {_SLAB_MIN_ORDER}"
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.rests_on(f"{_HERE}::test_slab_l1_order_against_case[krylov]")
def test_slab_l1_krylov_via_unified_vs_case() -> None:
    r"""**Slab Leg 1** — Krylov via the unified matvec ≡ Case reference.

    Case singular-eigenfunction reference ``sn_slab_1eg_2rg_S8`` at the
    matching S8 quadrature. The Krylov path is production's
    ``inner_solver="krylov"``: GMRES on the unified diamond-difference
    matvec, with no patch. Diamond difference converges to the Case reference
    at O(h²) (Leg 4); at n_per = 40 it reaches |Δk| = 2.16e-6 (``[M]``
    2026-10-04), inside :data:`_SLAB_REF_ABSTOL_AT_40`, which the step-scheme
    matvec (4.66e-5) and a matvec Σ_t off by 1e-3 (4.48e-4) both exceed.
    """
    k_ref = _slab_k_ref()
    k_krylov = _slab_k("krylov", 40)
    abs_err = abs(k_krylov - k_ref)
    assert abs_err < _SLAB_REF_ABSTOL_AT_40, (
        f"slab Krylov-via-unified vs Case reference at n_per=40: "
        f"k_krylov={k_krylov:.10f}, k_ref={k_ref:.10f}, |Δ|={abs_err:.3e}"
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.catches("ERR-025")
@pytest.mark.rests_on(f"{_HERE}::test_slab_l1_order_against_case[source_iteration]")
def test_slab_l1_sweep_vs_case() -> None:
    r"""**Slab Leg 2** — the sweep ≡ Case reference.

    Production's ``inner_solver="source_iteration"`` on the ``(L+C)`` sweep,
    whose per-cell algebra is ``affine_scan_coefficients``: at n_per = 40 it
    reaches |Δk| = 2.16e-6 (``[M]`` 2026-10-04), inside
    :data:`_SLAB_REF_ABSTOL_AT_40`.

    ERR-025's catcher in this module: re-dropping that defect's attenuation
    coefficient a = 2μ/(2μ + Δx·Σ_t) into the scan moves this k by 7.2e-2
    (``[M]``), while the Krylov rows do not move at all (their matvec does not
    read the scan coefficients), so the marker that sat on Leg 1 had decayed.
    """
    k_ref = _slab_k_ref()
    k_sweep = _slab_k("source_iteration", 40)
    abs_err = abs(k_sweep - k_ref)
    assert abs_err < _SLAB_REF_ABSTOL_AT_40, (
        f"slab sweep vs Case reference at n_per=40: "
        f"k_sweep={k_sweep:.10f}, k_ref={k_ref:.10f}, |Δ|={abs_err:.3e}"
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.parametrize("n_per", _SLAB_LADDER)
def test_slab_l1_sweep_vs_krylov_twin_path(n_per: int) -> None:
    r"""**Slab Leg 3** — sweep ≡ Krylov-via-unified at each mesh of the ladder.

    Both algorithms drive the SAME diamond-difference discrete fixed point on
    the same mesh, through cell kernels that share no code (the sweep's scan
    coefficients, the matvec's residual kernel). At matched solver tolerances
    they agree far inside :data:`_SLAB_TWIN_RTOL`. Disagreement signals
    algorithmic divergence (the L14 #6 manifestation pattern). The claim does
    not depend on the mesh, so it reads the ladder's solves; it formerly sat at
    n_per = 80, a mesh no other row solves.
    """
    k_sweep = _slab_k("source_iteration", n_per)
    k_krylov = _slab_k("krylov", n_per)
    rel = abs(k_sweep - k_krylov) / k_sweep
    assert rel < _SLAB_TWIN_RTOL, (
        f"slab n_per={n_per} twin-path disagreement: "
        f"k_sweep={k_sweep:.10f}, k_krylov={k_krylov:.10f}, rel={rel:.3e}"
    )
