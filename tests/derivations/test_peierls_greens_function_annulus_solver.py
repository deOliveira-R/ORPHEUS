r"""Phase-3C-2 numerical-implementation gates for the **annulus**
(hollow cylinder) Variant α Green's function reference (1G + MG,
homogeneous shell, INDEPENDENT specular reflectivities :math:`\alpha_{
\rm in}, \alpha_{\rm out} \in [0, 1]`, rank-2 boundary-to-boundary
scattering resolvent on through-rays + rank-1 outer-only closure with
cylinder 3D angular phase-space).

Mirrors :mod:`.test_peierls_greens_function_hollow_sphere_solver`
(Phase-3C-1 hollow sphere) with the cylinder 3D angular-correction
lift. Per Issue #129, this prototype stays **angle-resolved** —
DO NOT cross-check via the planar slab limit (the cylinder Bickley-
Naylor :math:`\mathrm{Ki}_n` form pre-integrates axial and gives ~22%
mismatch). Verification uses self-consistency, R_in → 0 reduction to
the verified solid cylinder, and exact analytical invariants
(:math:`k_\infty`).

Test strategy
-------------

1. **V_α1 closed-shell exactness (load-bearing)**: closed annulus
   with :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1` gives
   :math:`k_{\rm eff} = k_\infty` to machine precision. Structural
   composability check on the cylinder's 3D angular phase-space.
   Failure here HALTs Phase 3C-2.

2. **R_in → 0 limit**: annulus with :math:`R_{\rm in} = 10^{-3}
   R_{\rm out}` at vacuum-vacuum MUST converge to solid-cylinder
   vacuum. The phase-space partition reduces correctly as the
   through-ray subset's measure shrinks.

3. **Reflective-inner / vacuum-outer**: :math:`\alpha_{\rm in} = 1,
   \alpha_{\rm out} = 0`. Physical sanity — flux peaks near the
   reflective inner wall.

4. **Vacuum-inner / reflective-outer (cavity absorber)**:
   :math:`\alpha_{\rm in} = 0, \alpha_{\rm out} = 1`. Through-rays
   absorbed at inner; outer-only rays bounce normally.

5. **Symmetric reflective at intermediate α**: :math:`\alpha_{\rm in}
   = \alpha_{\rm out} = 0.5`. Rank-2 closure under symmetric BCs at
   non-trivial α.

6. **Multi-group asymmetric scattering at closed shell**: 2G with
   asymmetric Σ_s at :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1`
   reduces to :math:`k_\infty` to ≤ 1e-9 (ERR-002 anti-pattern
   detector).

7. **MG G=1 reduction**: at G=1, MG solver bit-equal to 1G solver.

8. **Convergence floor**: pin achieved accuracy at successive
   quadrature orders for the load-bearing reflective-inner /
   vacuum-outer case.

References
----------

- :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 3C-2
  annulus plan.
- :mod:`.test_peierls_greens_function_annulus_symbolic` — V_α1_annulus
  / V_α2_annulus / V_α2_annulus.aux / V_α3_annulus SymPy gates.
- :mod:`.test_peierls_greens_function_hollow_sphere_solver` — Phase-3C-1
  hollow sphere numerical gates (the analog template).
- :mod:`.test_peierls_greens_function_cylinder_solver` — Phase-1 solid
  cylinder numerical gates (the R_in → 0 reduction reference).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.eigenvalue import kinf_homogeneous
from orpheus.derivations.continuous.peierls_greens_function.greens_function_annulus import (
    solve_greens_function_annulus,
    solve_greens_function_annulus_mg,
)
from orpheus.derivations.continuous.peierls_greens_function.greens_function_cylinder import (
    solve_greens_function_cylinder,
)


# ═══════════════════════════════════════════════════════════════════════
# Fixtures — fuel-A-like XS
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
@pytest.mark.verifies("peierls-greens-annulus-architecture")
def test_v_alpha1_annulus_numerical_closed_BC_gives_kinf(
    fuelA_moderate_shell,
):
    r"""L1 — closed annulus :math:`\alpha_{\rm in} = \alpha_{\rm out}
    = 1` gives :math:`k_{\rm eff} = k_\infty` exactly (load-bearing).

    The structural composability check on the cylinder's 3D angular
    phase-space: BOTH the outer-only branch (rank-1, :math:`b > R_{
    \rm in}`) AND the through-ray branch (rank-2, :math:`b \le R_{\rm
    in}`) must independently produce :math:`\psi = q/\Sigma_t` for the
    closed-shell eigenmode, with the cylinder axial correction
    :math:`1/\sqrt{1 - \mu_{\rm axial}^2}` folded into the chord. A
    bug in the axial-correction lift would surface here.
    """
    fix = fuelA_moderate_shell
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_annulus(
        **fix, alpha_in=1.0, alpha_out=1.0,
        n_r=10, n_mu_axial=12, n_phi_az=24, n_traj_quad=24,
        max_iter=20, tol=1e-12,
    )
    assert res.converged, (
        f"closed annulus α_in=α_out=1: did not converge; "
        f"iter={res.iterations}"
    )
    assert res.iterations <= 2, (
        f"closed annulus: should converge in ≤2 iter from const initial "
        f"guess; got {res.iterations}"
    )
    np.testing.assert_allclose(
        res.k_eff, k_inf, rtol=1e-10,
        err_msg=(
            f"closed annulus α_in=α_out=1: k_eff = {res.k_eff} "
            f"≠ k_inf = {k_inf}"
        ),
    )
    rel_phi_spread = res.phi.std() / res.phi.mean()
    assert rel_phi_spread < 1e-10, (
        f"closed annulus: φ should be uniform; got "
        f"std/mean = {rel_phi_spread:.2e}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-annulus-architecture")
def test_closed_annulus_two_thicknesses(
    fuelA_thin_shell, fuelA_moderate_shell,
):
    r"""L1 — closed annulus gives k_inf at thin and moderate optical
    depth.

    The rank-2 + rank-1 composite closure must hold across different
    shell thicknesses on the cylinder geometry.
    """
    for fix in (fuelA_thin_shell, fuelA_moderate_shell):
        k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])
        res = solve_greens_function_annulus(
            **fix, alpha_in=1.0, alpha_out=1.0,
            n_r=8, n_mu_axial=10, n_phi_az=20, n_traj_quad=20,
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
# R_in → 0 SOLID-CYLINDER LIMIT TEST
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-annulus-impact-parameter-partition")
def test_R_in_to_zero_limit_matches_solid_cylinder_vacuum():
    r"""L1 — :math:`R_{\rm in} \to 0` limit converges to solid-cylinder
    vacuum.

    With :math:`R_{\rm in} = 10^{-3} R_{\rm out}` and :math:`\alpha_{
    \rm in} = 0, \alpha_{\rm out} = 0`, the through-ray subset
    (:math:`b \le R_{\rm in} = 10^{-3} R_{\rm out}`) has near-zero
    phase-space measure. The outer-only subset dominates, and the
    closure on it is structurally identical to the solid-cylinder
    vacuum closure.

    Plan target: 1e-3. Achieved (this prototype): ~3-5e-6.
    """
    sigma_t, sigma_s, nu_sigma_f = 0.5, 0.38, 0.025
    R_out = 5.0

    res_ann = solve_greens_function_annulus(
        R_in=1e-3, R_out=R_out,
        sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha_in=0.0, alpha_out=0.0,
        n_r=20, n_mu_axial=20, n_phi_az=24, n_traj_quad=32,
        max_iter=300, tol=1e-9,
    )
    res_cyl = solve_greens_function_cylinder(
        R=R_out,
        sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
        alpha=0.0,
        n_r=20, n_mu_axial=20, n_phi_az=24, n_traj_quad=32,
        max_iter=300, tol=1e-9,
    )
    assert res_ann.converged and res_cyl.converged
    rel_diff = abs(res_ann.k_eff - res_cyl.k_eff) / res_cyl.k_eff
    assert rel_diff < 1e-4, (
        f"R_in → 0 limit failed: annulus k_eff = {res_ann.k_eff:.10e}, "
        f"solid k_eff = {res_cyl.k_eff:.10e}, "
        f"rel diff = {rel_diff:.3e} > 1e-4 (plan target was 1e-3)"
    )


# ═══════════════════════════════════════════════════════════════════════
# Asymmetric BC: reflective inner, vacuum outer
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-annulus-architecture")
def test_reflective_inner_vacuum_outer_physical_sanity(
    fuelA_moderate_shell,
):
    r"""L1 — :math:`\alpha_{\rm in} = 1, \alpha_{\rm out} = 0` —
    flux peaks near the reflective inner wall.

    Reflective inner + vacuum outer: through-rays bounce off the inner
    wall and escape at the outer wall after one shell traversal. Flux
    peaks NEAR the reflective inner surface; k_eff < k_inf due to
    outer-surface leakage.
    """
    fix = fuelA_moderate_shell
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_annulus(
        **fix, alpha_in=1.0, alpha_out=0.0,
        n_r=20, n_mu_axial=20, n_phi_az=32, n_traj_quad=32,
        max_iter=300, tol=1e-10,
    )
    assert res.converged
    assert 0 < res.k_eff < k_inf, (
        f"refl-vac: k_eff={res.k_eff} should satisfy 0 < k_eff < "
        f"k_inf={k_inf}"
    )
    assert res.phi[0] > res.phi[-1], (
        f"refl-vac: phi(R_in≈)={res.phi[0]:.3e} should exceed "
        f"phi(R_out≈)={res.phi[-1]:.3e} (reflective inner has more flux)"
    )
    assert res.phi[0] / res.phi[-1] > 2.0, (
        f"refl-vac: phi ratio inner/outer = "
        f"{res.phi[0]/res.phi[-1]:.3f} should be > 2 for "
        f"reflective-inner vacuum-outer"
    )


# ═══════════════════════════════════════════════════════════════════════
# Asymmetric BC: vacuum inner (cavity absorber), reflective outer
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-annulus-architecture")
def test_vacuum_inner_reflective_outer_cavity_absorber(
    fuelA_moderate_shell,
):
    r"""L1 — :math:`\alpha_{\rm in} = 0, \alpha_{\rm out} = 1` —
    cavity acts as a perfect absorber.

    Through-rays (:math:`b \le R_{\rm in}`) are absorbed at the inner
    surface (lost to the cavity). Outer-only rays (:math:`b > R_{\rm
    in}`) bounce normally between two outer-cylinder points (rank-1
    closure with :math:`\alpha = 1`).

    k_eff is below k_inf because some flux is lost to the cavity.
    Flux profile peaks somewhere outward of the absorbing inner
    surface.
    """
    fix = fuelA_moderate_shell
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_annulus(
        **fix, alpha_in=0.0, alpha_out=1.0,
        n_r=20, n_mu_axial=20, n_phi_az=32, n_traj_quad=32,
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
@pytest.mark.verifies("peierls-greens-annulus-architecture")
def test_symmetric_reflective_intermediate_alpha(fuelA_moderate_shell):
    r"""L1 — :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0.5` converges
    to a sensible eigenvalue.

    Symmetric reflective at intermediate α. Validates the rank-2
    closure under symmetric BCs at non-trivial α (away from the
    rank-2 → rank-1 reductions at α = 0 and α = 1).
    """
    fix = fuelA_moderate_shell
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])
    alpha = 0.5

    res = solve_greens_function_annulus(
        **fix, alpha_in=alpha, alpha_out=alpha,
        n_r=16, n_mu_axial=16, n_phi_az=24, n_traj_quad=32,
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
@pytest.mark.verifies("peierls-greens-annulus-architecture")
def test_vacuum_vacuum_BC_leaks(fuelA_moderate_shell):
    r"""L1 — :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0` (vacuum-
    vacuum) leaks; rank-2 closure correctly bypasses to bare first-leg.

    Both surfaces absorbing: k_eff < k_inf, and the rank-2 prototype
    handles this with no special-case branch (the leading-α factors
    in the closure zero out psi_surf).
    """
    fix = fuelA_moderate_shell
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_annulus(
        **fix, alpha_in=0.0, alpha_out=0.0,
        n_r=16, n_mu_axial=16, n_phi_az=24, n_traj_quad=32,
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
@pytest.mark.verifies("peierls-greens-annulus-architecture")
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

    res = solve_greens_function_annulus_mg(
        R_in=1.0, R_out=3.0,
        sigma_t=sigma_t, sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f, chi=chi,
        alpha_in=1.0, alpha_out=1.0,
        n_r=10, n_mu_axial=12, n_phi_az=24, n_traj_quad=24,
        max_iter=100, tol=1e-11,
    )
    assert res.converged
    rel_err = abs(res.k_eff - k_inf) / k_inf
    assert rel_err < 1e-9, (
        f"MG 2G closed annulus: k_eff = {res.k_eff} ≠ k_inf = {k_inf} "
        f"by {rel_err:.3e}"
    )
    for g in range(2):
        spread = res.phi_g[g].std() / res.phi_g[g].mean()
        assert spread < 1e-6, (
            f"MG g={g}: φ_g spread {spread:.2e} not uniform"
        )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-annulus-architecture")
def test_mg_2g_asymmetric_BC_with_asymmetric_scattering():
    r"""L1 — 2G asymmetric scattering with :math:`\alpha_{\rm in} = 1,
    \alpha_{\rm out} = 0.3` (mixed BC) converges sensibly.

    Stress-tests the rank-2 MG path on the cylinder geometry at a
    non-trivial corner where neither closed-shell nor vacuum-vacuum
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

    res = solve_greens_function_annulus_mg(
        R_in=1.0, R_out=3.0,
        sigma_t=sigma_t, sigma_s=sigma_s,
        nu_sigma_f=nu_sigma_f, chi=chi,
        alpha_in=1.0, alpha_out=0.3,
        n_r=16, n_mu_axial=16, n_phi_az=24, n_traj_quad=32,
        max_iter=300, tol=1e-9,
    )
    assert res.converged
    assert 0 < res.k_eff < k_inf, (
        f"MG 2G mixed BC: k_eff={res.k_eff} should satisfy 0 < k_eff "
        f"< k_inf={k_inf}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-greens-annulus-architecture")
def test_mg_g1_reduces_to_1g_solver(fuelA_moderate_shell):
    r"""L1 — MG solver at G=1 reproduces the 1G solver bit-equal at
    :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0` vacuum.

    Architectural consistency check.
    """
    fix = fuelA_moderate_shell
    res_1g = solve_greens_function_annulus(
        **fix, alpha_in=0.0, alpha_out=0.0,
        n_r=10, n_mu_axial=12, n_phi_az=24, n_traj_quad=24,
        max_iter=200, tol=1e-10,
    )
    res_mg = solve_greens_function_annulus_mg(
        R_in=fix["R_in"], R_out=fix["R_out"],
        sigma_t=np.array([fix["sigma_t"]]),
        sigma_s=np.array([[fix["sigma_s"]]]),
        nu_sigma_f=np.array([fix["nu_sigma_f"]]),
        chi=np.array([1.0]),
        alpha_in=0.0, alpha_out=0.0,
        n_r=10, n_mu_axial=12, n_phi_az=24, n_traj_quad=24,
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
# Convergence floor — research-grade gate (slow)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-annulus-architecture")
def test_reflective_inner_vacuum_outer_convergence_floor():
    r"""L1 — :math:`\alpha_{\rm in} = 1, \alpha_{\rm out} = 0` k_eff
    convergence floor (slow, **research-grade gate**).

    Pins the achieved accuracy floor for the asymmetric reflective-
    vacuum case across a quadrature ladder. Convergence between the
    finest pair of grids serves as the proxy for the true k_eff.

    Quadrature ladder (n_r, n_mu, n_phi, n_traj_quad):

    - (10, 10, 16, 16) — coarse
    - (14, 14, 20, 24) — medium
    - (18, 18, 24, 32) — fine
    - (24, 24, 32, 48) — research-grade
    """
    sigma_t, sigma_s, nu_sigma_f = 0.5, 0.38, 0.025
    R_in, R_out = 1.0, 3.0

    orders = [
        (10, 10, 16, 16),
        (14, 14, 20, 24),
        (18, 18, 24, 32),
        (24, 24, 32, 48),
    ]
    k_vals = []

    for n_r, n_mu, n_phi, n_t in orders:
        res = solve_greens_function_annulus(
            R_in=R_in, R_out=R_out,
            sigma_t=sigma_t, sigma_s=sigma_s, nu_sigma_f=nu_sigma_f,
            alpha_in=1.0, alpha_out=0.0,
            n_r=n_r, n_mu_axial=n_mu, n_phi_az=n_phi, n_traj_quad=n_t,
            max_iter=500, tol=1e-11,
        )
        assert res.converged
        k_vals.append(res.k_eff)

    finest_pair_consistency = (
        abs(k_vals[-1] - k_vals[-2]) / k_vals[-1]
    )
    assert finest_pair_consistency < 1e-3, (
        f"refl-vac convergence at finest pair: "
        f"{orders[-2]} → {orders[-1]} differs by "
        f"{finest_pair_consistency:.3e}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Task-3 (V&V hardening) — off-diagonal intermediate α convergence
# ═══════════════════════════════════════════════════════════════════════
#
# Annulus has 3D phase space (r, μ_axial, φ_az) so off-diagonal
# convergence studies are expensive. One test case: α_in=0.4,
# α_out=0.7 at fuel-A τ_R≈1.5.
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-annulus-architecture")
def test_off_diagonal_intermediate_alpha_annulus():
    r"""L1 (slow) — annulus at off-diagonal intermediate α: convergence
    to-self gate.

    Self-consistency-only test (no external reference). Single
    parametrization: :math:`\alpha_{\rm in} = 0.4, \alpha_{\rm out}
    = 0.7`, fuel-A τ_R≈1.5. The 3D phase space (r, μ_axial,
    φ_az) makes a finer ladder expensive; the existing
    :func:`test_reflective_inner_vacuum_outer_convergence_floor`
    uses a (24, 24, 32, 48) cap and accepts a 1e-3 finest-pair gate.

    This test mirrors that pattern at the off-diagonal interior point.
    """
    fix = dict(
        R_in=1.0, R_out=3.0,
        sigma_t=0.5, sigma_s=0.38, nu_sigma_f=0.025,
    )
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    # Convergence is geometric-class noisy on annulus due to:
    #   - 3D phase space (r, μ_axial, φ_az);
    #   - impact-parameter partition surface b = R_in;
    #   - cylinder-specific 3D chord lift τ_annulus / sqrt(1 - μ²).
    # Measured 2026-05-02 at α=(0.4, 0.7), τ_R≈1.5:
    #   (14,14,20,24) → (18,18,24,32): 6e-4 rel
    #   (18,18,24,32) → (24,24,32,48): 2.9e-3 rel
    #   (24,24,32,48) → (32,28,40,64): 1.9e-3 rel (k oscillating
    #                                              about a mean)
    orders = [
        (18, 18, 24, 32),
        (24, 24, 32, 48),
        (32, 28, 40, 64),
    ]
    k_vals = []
    for n_r, n_mu, n_phi, n_t in orders:
        res = solve_greens_function_annulus(
            **fix, alpha_in=0.4, alpha_out=0.7,
            n_r=n_r, n_mu_axial=n_mu, n_phi_az=n_phi, n_traj_quad=n_t,
            max_iter=500, tol=1e-10,
        )
        assert res.converged, (
            f"annulus off-diagonal α=(0.4, 0.7): "
            f"order {(n_r, n_mu, n_phi, n_t)} did not converge"
        )
        k_vals.append(res.k_eff)

    # (c) Sensible band.
    k_finest = k_vals[-1]
    assert 0.0 < k_finest < k_inf, (
        f"annulus off-diagonal: k_eff = {k_finest} outside (0, k_inf)"
    )

    # (b) Finest-pair self-consistency ≤ 5e-3. Looser than slab-asym
    # (2e-5 achievable) and hollow sphere (1e-3 achievable) because
    # the 3D phase-space + impact-parameter partition + axial chord
    # lift compound the quadrature noise. Convergence appears
    # oscillatory about a mean rather than monotone-asymptotic.
    finest_pair = abs(k_vals[-1] - k_vals[-2]) / k_vals[-1]
    assert finest_pair < 5e-3, (
        f"annulus off-diagonal α=(0.4, 0.7): finest-pair "
        f"self-consistency = {finest_pair:.3e} exceeds 5e-3 gate. "
        f"k_vals = {k_vals}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Task-4 (V&V hardening) — grazing-ray stability
# ═══════════════════════════════════════════════════════════════════════
#
# Annulus grazing-ray loci:
#   - b → R_in (impact parameter tangent to inner cylinder).
#   - μ_axial → ±1 (axial-only ray, finite chord but with the
#     τ_annulus = τ_hollow_sph / sqrt(1 - μ_axial²) chord lift).
#
# We refine n_phi_az (controls b-resolution near R_in) and verify
# stability. The μ_axial ladder is fixed at a moderate value because
# pre-existing tests already exercise μ_axial → ±1 indirectly via
# the chord-lift formula in the closed-shell V_α1_annulus test.
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.slow
@pytest.mark.verifies("peierls-greens-annulus-architecture")
def test_grazing_ray_stability_annulus():
    r"""L1 (slow) — annulus grazing-ray (b → R_in and μ_axial → ±1)
    stability under :math:`n_\phi` refinement.

    Refines :math:`n_\phi \in \{16, 24, 32\}` at fixed
    :math:`(n_r, n_{\mu\,{\rm axial}}, n_{\rm traj})` and α=0.5 on
    both surfaces. Verify stability (no NaN, no oscillation,
    no blow-up).

    Configuration: fuel-A τ_R≈1.5 (R_in=1, R_out=3, σ_t=0.5),
    α_in = α_out = 0.5.
    """
    fix = dict(
        R_in=1.0, R_out=3.0,
        sigma_t=0.5, sigma_s=0.38, nu_sigma_f=0.025,
    )

    n_phi_ladder = [16, 24, 32]
    k_vals = []
    for n_phi in n_phi_ladder:
        res = solve_greens_function_annulus(
            **fix, alpha_in=0.5, alpha_out=0.5,
            n_r=14, n_mu_axial=14, n_phi_az=n_phi, n_traj_quad=32,
            max_iter=500, tol=1e-10,
        )
        assert res.converged, (
            f"annulus grazing-ray: n_phi={n_phi} did not converge"
        )
        assert np.isfinite(res.k_eff), (
            f"annulus grazing-ray: n_phi={n_phi} non-finite k_eff"
        )
        assert np.all(np.isfinite(res.psi)), (
            f"annulus grazing-ray: n_phi={n_phi} non-finite ψ"
        )
        k_vals.append(res.k_eff)

    # 3D phase-space + impact-parameter partition + chord-lift =
    # geometric-class quadrature error dominates. Achieved at
    # α=(0.5, 0.5), τ_R≈1.5:
    #   n_phi 16 → 24: ~4.9e-3 rel
    #   n_phi 24 → 32: ~1.1e-3 rel
    # k is monotone-converging across the ladder. The 6e-3 gate
    # captures the 16→24 step plus margin and gates against actual
    # blow-up (NaN, sign flips, large oscillations).
    for i in range(len(n_phi_ladder) - 1):
        rel = abs(k_vals[i + 1] - k_vals[i]) / k_vals[-1]
        assert rel < 6e-3, (
            f"annulus grazing-ray: n_phi={n_phi_ladder[i]} → "
            f"{n_phi_ladder[i+1]} differs by {rel:.3e}, exceeds "
            f"6e-3 stability gate. k_vals = {k_vals}"
        )

    # Monotonicity check: differences should shrink (convergence)
    # OR sign-flip with shrinking magnitude (overshoot converging).
    diffs = [k_vals[i + 1] - k_vals[i] for i in range(len(k_vals) - 1)]
    diff_mags = [abs(d) for d in diffs]
    # Either no sign flip (monotone) OR magnitude must shrink across
    # the ladder.
    sign_change = any(diffs[i] * diffs[i + 1] < 0 for i in range(len(diffs) - 1))
    if sign_change:
        assert diff_mags[-1] < diff_mags[0], (
            f"annulus grazing-ray: k_eff oscillates without "
            f"convergence. k_vals = {k_vals}"
        )
