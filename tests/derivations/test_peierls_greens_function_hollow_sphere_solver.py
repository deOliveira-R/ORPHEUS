r"""Phase-3C-1 numerical-implementation gates for the **hollow sphere**
Variant α Green's function reference (1G + MG, homogeneous shell,
INDEPENDENT specular reflectivities :math:`\alpha_{\rm in} \in [0, 1]`
and :math:`\alpha_{\rm out} \in [0, 1]`, rank-2 boundary-to-boundary
scattering resolvent on through-rays + rank-1 outer-only closure).

Mirrors :mod:`.test_peierls_greens_function_slab_asymmetric_solver`
(Phase-3B asymmetric slab) with the curvilinear 2-surface generalisation
(impact-parameter phase-space partition).

Test strategy
-------------

1. **V_α1 closed-shell exactness (load-bearing)**: closed shell with
   :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1` gives
   :math:`k_{\rm eff} = k_\infty` to machine precision. This is the
   structural composability check — the impact-parameter partition
   must compose cleanly with the rank-2 closure on the through-ray
   subset. Failure here would HALT Phase 3C-1.

2. **R_in → 0 limit**: hollow with :math:`R_{\rm in} = 10^{-3} R_{\rm
   out}` at :math:`\alpha_{\rm in} = 0, \alpha_{\rm out} = 0` MUST
   converge to solid-sphere Variant α at the same outer radius and
   :math:`\alpha = 0`. The structural identity that the inner cavity
   contributes zero phase-space measure as :math:`R_{\rm in} \to 0`.

3. **Asymmetric reflective-inner / vacuum-outer**:
   :math:`\alpha_{\rm in} = 1, \alpha_{\rm out} = 0` — physically
   sensible flux profile, k_eff < k_inf.

4. **Asymmetric vacuum-inner / reflective-outer**:
   :math:`\alpha_{\rm in} = 0, \alpha_{\rm out} = 1` — the cavity-
   absorber case. Through-rays are lost to the cavity; outer-only
   rays bounce normally. k_eff < k_inf.

5. **Symmetric reflective intermediate**: :math:`\alpha_{\rm in} =
   \alpha_{\rm out} = 0.5` — well-behaved rank-2 closure under
   symmetric BCs.

6. **Multi-group asymmetric scattering at closed shell**: 2G with
   asymmetric scattering matrix at :math:`\alpha_{\rm in} = \alpha_{
   \rm out} = 1` reduces to :math:`k_\infty` from
   :func:`kinf_homogeneous` to :math:`\le 10^{-9}` (ERR-002 anti-
   pattern detector).

7. **MG G=1 reduction**: at G=1, MG solver bit-equal to 1G solver at
   vacuum-vacuum BC.

8. **Convergence floor**: pin the achieved accuracy at successive
   quadrature orders for the load-bearing reflective-inner / vacuum-
   outer case.

9. **Sphere/cylinder/slab regression**: not in this file directly,
   but the full acceptance gate ensures these tests still pass at
   bit-equal accuracy floors.

References
----------

- :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 3C-1
  hollow sphere plan.
- :mod:`.test_peierls_greens_function_hollow_sphere_symbolic` —
  V_α1_hollow_sph/V_α2_hollow_sph/V_α3_hollow_sph SymPy gates.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import kinf_homogeneous
from orpheus.derivations.continuous.peierls_greens_function.greens_function import (
    solve_greens_function_sphere,
)
from orpheus.derivations.continuous.peierls_greens_function.greens_function_hollow_sphere import (
    solve_greens_function_hollow_sphere,
    solve_greens_function_hollow_sphere_mg,
)


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — fuel-A-like XS analogous to sphere/slab tests
# ═══════════════════════════════════════════════════════════════════════


@pytest.fixture(scope="module")
def fuelA_thin_shell():
    """fuel-A-like XS, R_in=0.5, R_out=2.5 (τ shell ≈ 1)."""
    return {
        "R_in": 0.5, "R_out": 2.5,
        "sigma_t": 0.5, "sigma_s": 0.38, "nu_sigma_f": 0.025,
    }


@pytest.fixture(scope="module")
def fuelA_moderate_shell():
    """fuel-A-like XS, R_in=1.0, R_out=3.0 (moderate shell)."""
    return {
        "R_in": 1.0, "R_out": 3.0,
        "sigma_t": 0.5, "sigma_s": 0.38, "nu_sigma_f": 0.025,
    }


# ═══════════════════════════════════════════════════════════════════════
# THE LOAD-BEARING V_α1 CLOSED-SHELL TEST (composability check)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-hollow-sph-architecture")
def test_v_alpha1_hollow_sph_numerical_closed_BC_gives_kinf(
    fuelA_moderate_shell,
):
    r"""L1 — closed shell :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1`
    gives :math:`k_{\rm eff} = k_\infty` exactly (load-bearing).

    The structural composability check: BOTH the outer-only branch
    (rank-1, :math:`b > R_{\rm in}`) AND the through-ray branch
    (rank-2, :math:`b \le R_{\rm in}`) must independently produce
    :math:`\psi = q/\Sigma_t` for the closed-shell eigenmode. Failure
    here would mean the impact-parameter partition does not compose
    cleanly with the rank-2 closure — a structural bug that would
    HALT Phase 3C-1 per plan.
    """
    fix = fuelA_moderate_shell
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_hollow_sphere(
        **fix, alpha_in=1.0, alpha_out=1.0,
        n_r=12, n_mu=16, n_traj_quad=32,
        max_iter=50, tol=1e-12,
    )
    assert res.converged, (
        f"closed shell α_in=α_out=1: did not converge; iter={res.iterations}"
    )
    assert res.iterations <= 2, (
        f"closed shell: should converge in ≤2 iter from const initial "
        f"guess; got {res.iterations}"
    )
    np.testing.assert_allclose(
        res.k_eff, k_inf, rtol=1e-10,
        err_msg=(
            f"closed shell α_in=α_out=1: k_eff = {res.k_eff} "
            f"≠ k_inf = {k_inf}"
        ),
    )
    rel_phi_spread = res.phi.std() / res.phi.mean()
    assert rel_phi_spread < 1e-10, (
        f"closed shell: φ should be uniform; got "
        f"std/mean = {rel_phi_spread:.2e}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-hollow-sph-architecture")
def test_closed_shell_two_thicknesses(fuelA_thin_shell, fuelA_moderate_shell):
    r"""L1 — closed shell gives k_inf at thin and moderate optical
    depth.

    The rank-2 + rank-1 composite closure must hold across different
    shell thicknesses when both surfaces are fully reflective.
    """
    for fix in (fuelA_thin_shell, fuelA_moderate_shell):
        k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])
        res = solve_greens_function_hollow_sphere(
            **fix, alpha_in=1.0, alpha_out=1.0,
            n_r=10, n_mu=12, n_traj_quad=24,
            max_iter=20, tol=1e-12,
        )
        np.testing.assert_allclose(
            res.k_eff, k_inf, rtol=1e-10,
            err_msg=(
                f"R_in={fix['R_in']}, R_out={fix['R_out']}: "
                f"k_eff = {res.k_eff} ≠ k_inf = {k_inf}"
            ),
        )


# ═══════════════════════════════════════════════════════════════════════
# R_in → 0 SOLID-SPHERE LIMIT TEST
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-hollow-sph-impact-parameter-partition")
def test_R_in_to_zero_limit_matches_solid_sphere_vacuum():
    r"""L1 — :math:`R_{\rm in} \to 0` limit converges to solid-sphere
    vacuum.

    With :math:`R_{\rm in} = 10^{-3} R_{\rm out}` and
    :math:`\alpha_{\rm in} = 0, \alpha_{\rm out} = 0`, the through-ray
    subset (:math:`b \le R_{\rm in} = 10^{-3} R_{\rm out}`) has
    near-zero phase-space measure. The outer-only subset dominates,
    and the closure on it is structurally identical to the solid-
    sphere vacuum closure.

    This is the canonical structural-independence check for the
    impact-parameter partition: the partition's behaviour at the
    R_in → 0 limit must reduce to the solid-sphere result (which is
    itself verified independently against PS-1982).

    Plan target: 1e-3 rel agreement. Achieved (this prototype): ≤ 1e-7
    even at coarse quadrature, well below target.
    """
    sigma_t, sigma_s, nu_sigma_f = 0.5, 0.38, 0.025
    R_out = 5.0

    res_hollow = solve_greens_function_hollow_sphere(
        R_in=1e-3, R_out=R_out,
        sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha_in=0.0, alpha_out=0.0,
        n_r=24, n_mu=24, n_traj_quad=48,
        max_iter=300, tol=1e-9,
    )
    res_solid = solve_greens_function_sphere(
        R=R_out,
        sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha=0.0,
        n_r=24, n_mu=24, n_traj_quad=48,
        max_iter=300, tol=1e-9,
    )
    assert res_hollow.converged and res_solid.converged
    rel_diff = abs(res_hollow.k_eff - res_solid.k_eff) / res_solid.k_eff
    assert rel_diff < 1e-7, (
        f"R_in → 0 limit failed: hollow k_eff = {res_hollow.k_eff:.10e}, "
        f"solid k_eff = {res_solid.k_eff:.10e}, "
        f"rel diff = {rel_diff:.3e} > 1e-7"
    )


# ═══════════════════════════════════════════════════════════════════════
# Asymmetric BC: reflective inner, vacuum outer
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-hollow-sph-architecture")
def test_reflective_inner_vacuum_outer_physical_sanity(fuelA_moderate_shell):
    r"""L1 — :math:`\alpha_{\rm in} = 1, \alpha_{\rm out} = 0` —
    flux peaks near the reflective inner wall.

    Reflective inner + vacuum outer: through-rays bounce off the inner
    wall and escape at the outer wall after one shell traversal. The
    flux should peak NEAR the reflective inner surface (hottest
    region) and decay toward the vacuum outer surface. k_eff < k_inf
    due to outer-surface leakage.

    The flux profile may exhibit a small interior maximum slightly
    away from R_in due to spherical-shell volumetric weighting (more
    available shell volume at intermediate r). The strict "monotone
    decreasing from R_in" expectation does not strictly hold; the
    physical sanity check is that phi(R_in) > phi(R_out).
    """
    fix = fuelA_moderate_shell
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_hollow_sphere(
        **fix, alpha_in=1.0, alpha_out=0.0,
        n_r=24, n_mu=32, n_traj_quad=48,
        max_iter=300, tol=1e-10,
    )
    assert res.converged
    assert 0 < res.k_eff < k_inf, (
        f"refl-vac: k_eff={res.k_eff} should satisfy 0 < k_eff < "
        f"k_inf={k_inf}"
    )
    # Inner-surface flux should exceed outer-surface flux.
    assert res.phi[0] > res.phi[-1], (
        f"refl-vac: phi(R_in≈)={res.phi[0]:.3e} should exceed "
        f"phi(R_out≈)={res.phi[-1]:.3e} (reflective wall has more flux)"
    )
    # Ratio should be substantial (factor ≥ 2) for τ_shell = 1.
    assert res.phi[0] / res.phi[-1] > 2.0, (
        f"refl-vac: phi ratio inner/outer = "
        f"{res.phi[0]/res.phi[-1]:.3f} should be > 2 for "
        f"reflective-inner vacuum-outer"
    )


# ═══════════════════════════════════════════════════════════════════════
# Asymmetric BC: vacuum inner (cavity absorber), reflective outer
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-hollow-sph-architecture")
def test_vacuum_inner_reflective_outer_cavity_absorber(fuelA_moderate_shell):
    r"""L1 — :math:`\alpha_{\rm in} = 0, \alpha_{\rm out} = 1` —
    cavity acts as a perfect absorber.

    Through-rays (:math:`b \le R_{\rm in}`) are absorbed at the inner
    surface (lost to the cavity). Outer-only rays
    (:math:`b > R_{\rm in}`) bounce normally between two outer-surface
    points (rank-1 closure with :math:`\alpha = 1`).

    k_eff is below k_inf because some flux is lost to the cavity, but
    the rank-2 closure remains well-defined (det = 1 when α_in = 0).
    The flux profile should peak somewhere in the shell interior
    (away from the absorbing inner surface, biased toward the
    reflective outer surface).
    """
    fix = fuelA_moderate_shell
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_hollow_sphere(
        **fix, alpha_in=0.0, alpha_out=1.0,
        n_r=24, n_mu=32, n_traj_quad=48,
        max_iter=300, tol=1e-10,
    )
    assert res.converged
    assert 0 < res.k_eff < k_inf, (
        f"vac-refl (cavity absorber): k_eff={res.k_eff} should "
        f"satisfy 0 < k_eff < k_inf={k_inf}"
    )
    # Cavity absorber: flux at outer (reflective) wall exceeds flux
    # at inner (absorbing) wall.
    assert res.phi[-1] > res.phi[0], (
        f"vac-refl: phi(R_out≈)={res.phi[-1]:.3e} should exceed "
        f"phi(R_in≈)={res.phi[0]:.3e} (cavity absorbs flux)"
    )


# ═══════════════════════════════════════════════════════════════════════
# Symmetric reflective intermediate
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-hollow-sph-architecture")
def test_symmetric_reflective_intermediate_alpha(fuelA_moderate_shell):
    r"""L1 — :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0.5` converges
    to a sensible eigenvalue.

    Symmetric reflective at intermediate α. Validates the rank-2
    closure under symmetric BCs at non-trivial α (away from the
    rank-2 → rank-1 reductions at α = 0 and α = 1). Sanity bounds:
    k_eff > 0 (physical) and k_eff < k_inf (some leakage, since
    reflectivity < 1 means not all flux returns).
    """
    fix = fuelA_moderate_shell
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])
    alpha = 0.5

    res = solve_greens_function_hollow_sphere(
        **fix, alpha_in=alpha, alpha_out=alpha,
        n_r=20, n_mu=24, n_traj_quad=48,
        max_iter=300, tol=1e-10,
    )
    assert res.converged
    assert 0 < res.k_eff < k_inf, (
        f"sym refl α={alpha}: k_eff={res.k_eff} should satisfy "
        f"0 < k_eff < k_inf={k_inf}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Vacuum-vacuum: phase-space partition fully tested via R_in → 0 above
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-hollow-sph-architecture")
def test_vacuum_vacuum_BC_leaks(fuelA_moderate_shell):
    r"""L1 — :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0` (vacuum-
    vacuum) leaks; rank-2 closure correctly bypasses to bare first-leg.

    Both surfaces absorbing: k_eff < k_inf, and the rank-2 prototype
    handles this with no special-case branch (the leading-α factors
    in the closure zero out psi_surf).
    """
    fix = fuelA_moderate_shell
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_hollow_sphere(
        **fix, alpha_in=0.0, alpha_out=0.0,
        n_r=20, n_mu=24, n_traj_quad=48,
        max_iter=300, tol=1e-10,
    )
    assert res.converged
    assert 0 < res.k_eff < k_inf, (
        f"vac-vac: k_eff={res.k_eff} should satisfy 0 < k_eff < "
        f"k_inf={k_inf}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Multi-group asymmetric scattering at α_in=α_out=1
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-hollow-sph-architecture")
def test_mg_2g_asymmetric_scattering_at_closed_shell_gives_kinf():
    r"""L1 — 2G asymmetric scattering at :math:`\alpha_{\rm in} =
    \alpha_{\rm out} = 1` MUST give k_inf to ≤ 1e-9.

    Anti-pattern detector for the SigS-transpose convention drift
    (ERR-002). Asymmetric scattering matrix is the discriminator
    (1G is degenerate, symmetric matrices are blind).
    """
    sigma_t = np.array([1.0, 0.6])
    sigma_s = np.array([
        [0.2, 0.5],
        [0.05, 0.3],
    ])
    nu_sigma_f = np.array([0.05, 0.4])
    chi = np.array([1.0, 0.0])
    k_inf = kinf_homogeneous(sigma_t, sigma_s, nu_sigma_f, chi)

    res = solve_greens_function_hollow_sphere_mg(
        R_in=1.0, R_out=3.0,
        sigma_t=sigma_t, sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f, chi=chi,
        alpha_in=1.0, alpha_out=1.0,
        n_r=12, n_mu=16, n_traj_quad=32, max_iter=100, tol=1e-11,
    )
    assert res.converged
    rel_err = abs(res.k_eff - k_inf) / k_inf
    assert rel_err < 1e-9, (
        f"MG 2G closed shell: k_eff = {res.k_eff} ≠ k_inf = {k_inf} "
        f"by {rel_err:.3e}"
    )
    for g in range(2):
        spread = res.phi_g[g].std() / res.phi_g[g].mean()
        assert spread < 1e-6, (
            f"MG g={g}: φ_g spread {spread:.2e} not uniform"
        )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-hollow-sph-architecture")
def test_mg_2g_asymmetric_BC_with_asymmetric_scattering():
    r"""L1 — 2G asymmetric scattering with :math:`\alpha_{\rm in} = 1,
    \alpha_{\rm out} = 0.3` (mixed BC) converges sensibly.

    Stress-tests the rank-2 MG path at a non-trivial corner of the
    BC parameter square where neither closed-shell nor vacuum-vacuum
    simplifications apply. Asserts k_eff < k_inf (leakage) and
    k_eff > 0 (physical).
    """
    sigma_t = np.array([1.0, 0.6])
    sigma_s = np.array([
        [0.2, 0.5],
        [0.05, 0.3],
    ])
    nu_sigma_f = np.array([0.05, 0.4])
    chi = np.array([1.0, 0.0])
    k_inf = kinf_homogeneous(sigma_t, sigma_s, nu_sigma_f, chi)

    res = solve_greens_function_hollow_sphere_mg(
        R_in=1.0, R_out=3.0,
        sigma_t=sigma_t, sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f, chi=chi,
        alpha_in=1.0, alpha_out=0.3,
        n_r=20, n_mu=24, n_traj_quad=48, max_iter=300, tol=1e-9,
    )
    assert res.converged
    assert 0 < res.k_eff < k_inf, (
        f"MG 2G mixed BC: k_eff={res.k_eff} should satisfy 0 < k_eff "
        f"< k_inf={k_inf}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-hollow-sph-architecture")
def test_mg_g1_reduces_to_1g_solver(fuelA_moderate_shell):
    r"""L1 — MG solver at G=1 reproduces the 1G solver bit-equal at
    :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0` vacuum.

    Architectural consistency check.
    """
    fix = fuelA_moderate_shell
    res_1g = solve_greens_function_hollow_sphere(
        **fix, alpha_in=0.0, alpha_out=0.0,
        n_r=12, n_mu=16, n_traj_quad=32,
        max_iter=200, tol=1e-10,
    )
    res_mg = solve_greens_function_hollow_sphere_mg(
        R_in=fix["R_in"], R_out=fix["R_out"],
        sigma_t=np.array([fix["sigma_t"]]),
        sigma_s=np.array([[fix["sigma_s"]]]),
        nu_sigma_f=np.array([fix["nu_sigma_f"]]),
        chi=np.array([1.0]),
        alpha_in=0.0, alpha_out=0.0,
        n_r=12, n_mu=16, n_traj_quad=32,
        max_iter=200, tol=1e-10,
    )
    assert res_1g.converged and res_mg.converged
    np.testing.assert_allclose(
        res_mg.k_eff, res_1g.k_eff, rtol=1e-9,
        err_msg=(
            f"MG G=1 reduction: MG k_eff = {res_mg.k_eff} != 1G "
            f"k_eff = {res_1g.k_eff}"
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# Convergence floor — research-grade gate
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-hollow-sph-architecture")
def test_reflective_inner_vacuum_outer_convergence_floor():
    r"""L1 — :math:`\alpha_{\rm in} = 1, \alpha_{\rm out} = 0` k_eff
    convergence floor (slow, **research-grade gate**).

    Pins the achieved accuracy floor for the asymmetric reflective-
    vacuum case across a quadrature ladder. Convergence between the
    finest pair of grids serves as the proxy for the true k_eff.

    Quadrature ladder (n_r, n_mu, n_traj_quad):

    - (16, 16, 32) — coarse
    - (24, 24, 48) — medium
    - (32, 32, 64) — fine
    - (48, 48, 96) — research-grade
    """
    sigma_t, sigma_s, nu_sigma_f = 0.5, 0.38, 0.025
    R_in, R_out = 1.0, 3.0

    orders = [(16, 16, 32), (24, 24, 48), (32, 32, 64), (48, 48, 96)]
    k_vals = []

    for n_r, n_mu, n_t in orders:
        res = solve_greens_function_hollow_sphere(
            R_in=R_in, R_out=R_out,
            sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
            alpha_in=1.0, alpha_out=0.0,
            n_r=n_r, n_mu=n_mu, n_traj_quad=n_t,
            max_iter=500, tol=1e-11,
        )
        assert res.converged
        k_vals.append(res.k_eff)

    # Convergence at finest pair: tighter than 1e-3.
    finest_pair_consistency = (
        abs(k_vals[-1] - k_vals[-2]) / k_vals[-1]
    )
    assert finest_pair_consistency < 1e-3, (
        f"refl-vac convergence at finest pair: "
        f"{orders[-2]} → {orders[-1]} differs by "
        f"{finest_pair_consistency:.3e}"
    )
