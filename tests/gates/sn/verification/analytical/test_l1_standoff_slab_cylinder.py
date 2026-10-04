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
``StreamingCollisionOperator.apply`` with the sweep as preconditioner.

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
@pytest.mark.parametrize("nx", [20, 40, 80])
def test_cylinder_l1_refinement_both_paths(nx: int) -> None:
    r"""**Cylinder Leg 4, twin half** — sweep ≡ Krylov-via-unified at each refinement.

    For each nx ∈ {20, 40, 80} the sweep and the Krylov-via-unified paths
    must reach the same discrete fixed point (:data:`_CYL_TWIN_RTOL`). The
    reference half of this leg is
    ``test_cylinder_l1_refinement_against_reference``, split off on
    2026-09-26 when the reference legs became strict xfails on #516, so that
    this assertion stays live.

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
    """
    for solution in (_solve_cyl_via_sweep(nx=nx), _solve_cyl_via_krylov_unified(nx=nx)):
        verify_cylinder_k(solution, _CYL_REF_RTOL)


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
    return float(continuous_get("sn_slab_1eg_2rg_S8").k_eff)


def _solve_slab_via_krylov_unified(n_per: int) -> float:
    mesh, materials, N_ord = _build_slab_2region_mesh(n_per=n_per)
    quad = Quadrature.gauss_legendre(N_ord)
    sol = solve_sn(
            materials, mesh, quad,
            inner_solver="krylov",
            max_outer=500, max_inner=500,
            keff_tol=1e-12, inner_tol=1e-9,
        )
    return float(sol.outcome.keff)


def _solve_slab_via_sweep(n_per: int) -> float:
    mesh, materials, N_ord = _build_slab_2region_mesh(n_per=n_per)
    quad = Quadrature.gauss_legendre(N_ord)
    sol = solve_sn(
        materials, mesh, quad,
        inner_solver="source_iteration",
        max_outer=500, max_inner=500,
        keff_tol=1e-12, inner_tol=1e-12,
    )
    return float(sol.outcome.keff)


# ═══════════════════════════════════════════════════════════════════════
# Slab L1 standoff
# ═══════════════════════════════════════════════════════════════════════
# Reference tolerance budget (per existing ``test_sn_2region_reflective_case_eigenvalue``):
#   |Δk| < 1e-5 at n_per=320 (sweep finest-mesh budget).
#   At n_per=160 the production sweep converges to |Δk| ~ 1e-7 typically
#   (still O(h) toward the Case reference); allow |Δk| < 2e-5 for safety
#   margin at the coarser meshes used here.
# Twin-path tolerance: 1e-6 rel.
_SLAB_REF_ABSTOL_AT_160 = 2.0e-5
_SLAB_TWIN_RTOL = 1.0e-5


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.catches("ERR-025")
def test_slab_l1_krylov_via_unified_vs_case() -> None:
    r"""**Slab Leg 1** — Krylov via unified matvec ≡ Case reference.

    Case singular-eigenfunction reference ``sn_slab_1eg_2rg_S8`` at the
    matching S8 quadrature.  Krylov is monkey-patched to route through
    the unified WDD matvec (vs the legacy 1st-order FD).  WDD converges
    to the Case reference at O(h) (piecewise-Σ interface degrades from
    nominal O(h²)).  At n_per=160 the production sweep reaches |Δk| ~
    1e-7; we allow 2e-5 for tolerance budget headroom.
    """
    k_ref = _slab_k_ref()
    k_krylov = _solve_slab_via_krylov_unified(n_per=160)
    abs_err = abs(k_krylov - k_ref)
    assert abs_err < _SLAB_REF_ABSTOL_AT_160, (
        f"slab Krylov-via-unified vs Case reference: "
        f"k_krylov={k_krylov:.10f}, k_ref={k_ref:.10f}, |Δ|={abs_err:.3e}"
    )


@pytest.mark.l1
@pytest.mark.slow
def test_slab_l1_sweep_vs_krylov_twin_path() -> None:
    r"""**Slab Leg 3** — sweep ≡ Krylov-via-unified twin-path agreement.

    Both algorithms drive the SAME WDD discrete fixed point on the same
    mesh.  At matched solver tolerances they agree at sub-ULP-of-
    iteration drift.  Disagreement signals algorithmic divergence (the
    L14 #6 manifestation pattern).
    """
    k_sweep = _solve_slab_via_sweep(n_per=80)
    k_krylov = _solve_slab_via_krylov_unified(n_per=80)
    rel = abs(k_sweep - k_krylov) / k_sweep
    assert rel < _SLAB_TWIN_RTOL, (
        f"slab twin-path disagreement: "
        f"k_sweep={k_sweep:.10f}, k_krylov={k_krylov:.10f}, rel={rel:.3e}"
    )


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.parametrize("n_per", [40, 80, 160])
def test_slab_l1_refinement_both_paths(n_per: int) -> None:
    r"""**Slab Leg 4** — both paths converge to Case ref + agree at each refinement.

    For each n_per ∈ {40, 80, 160}, both sweep and Krylov-via-unified
    must agree with the Case reference (decreasing error under
    refinement) and with each other (twin-path agreement).  The existing
    sweep refinement test in ``test_heterogeneous_transport.py`` covers
    n_per ∈ {20, 40, 80, 160, 320} for the sweep alone; this test adds
    the Krylov-via-unified leg.
    """
    k_ref = _slab_k_ref()
    k_sweep = _solve_slab_via_sweep(n_per=n_per)
    k_krylov = _solve_slab_via_krylov_unified(n_per=n_per)
    abs_err_sweep = abs(k_sweep - k_ref)
    abs_err_krylov = abs(k_krylov - k_ref)
    rel_twin = abs(k_sweep - k_krylov) / k_sweep
    # Each path must converge toward the Case reference.  Looser
    # tolerance at coarser meshes (O(h) interface error scales linearly):
    # 5e-4 at n=40, 2.5e-4 at n=80, 2e-5 at n=160.
    tolerance = {40: 5.0e-4, 80: 2.5e-4, 160: 2.0e-5}[n_per]
    assert abs_err_sweep < tolerance, (
        f"slab n_per={n_per}: sweep |Δk|={abs_err_sweep:.3e} ≥ {tolerance:.0e}"
    )
    assert abs_err_krylov < tolerance, (
        f"slab n_per={n_per}: krylov |Δk|={abs_err_krylov:.3e} ≥ {tolerance:.0e}"
    )
    assert rel_twin < _SLAB_TWIN_RTOL, (
        f"slab n_per={n_per}: twin-path rel={rel_twin:.3e} ≥ {_SLAB_TWIN_RTOL:.0e}"
    )
