r"""B4 numerical-implementation gate for the Plan 2 Variant α prototype.

V_α1 algebraically proves that for homogeneous sphere with perfect
specular BC and isotropic scattering, the rank-1 isotropic mode is the
unique eigenmode with :math:`k_{\rm eff} = k_\infty = \nu\Sigma_f /
\Sigma_a`. These tests verify the **numerical implementation** in
:mod:`orpheus.derivations.continuous.peierls.greens_function` matches
that algebraic identity within quadrature error.

Test strategy:

1. **Constant initial guess** (V_α1.numerical) — :math:`\psi_0 = 1`
   should land on the eigenmode in one iteration; k_eff = k_inf to
   machine precision.
2. **Non-uniform initial guess** (sinusoidal perturbation) — power
   iteration should converge to the same isotropic eigenmode and
   recover k_inf within 0.05 %.
3. **Two thicknesses** — thin (τ_R = 2.5) and moderate (τ_R = 5)
   sphere should both yield k_eff = k_inf since closed sphere has no
   leakage. Validates that the bounce-sum geometric series is wired
   correctly across the τ_R range.

Predecessor:
:mod:`tests.derivations.test_peierls_greens_function_symbolic`
(V_α1, V_α2, V_α3 SymPy gates).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.peierls.greens_function import (
    solve_greens_function_specular_sphere,
)


@pytest.fixture(scope="module")
def fuelA_like_thin_sphere():
    """Same fuel-A-like XS used by the existing peierls_specular tests:
    σ_t = 0.5, σ_s = 0.38, νσ_f = 0.025, R = 5; τ_R = 2.5."""
    return {
        "R": 5.0, "sigma_t": 0.5, "sigma_s": 0.38, "nu_sigma_f": 0.025,
    }


@pytest.fixture(scope="module")
def fuelA_like_moderate_sphere():
    """Moderate-thickness sphere: same XS, R = 10; τ_R = 5."""
    return {
        "R": 10.0, "sigma_t": 0.5, "sigma_s": 0.38, "nu_sigma_f": 0.025,
    }


@pytest.mark.foundation
def test_v_alpha1_numerical_constant_initial_guess(fuelA_like_thin_sphere):
    r"""V_α1.numerical — constant ψ_0 lands on eigenmode in 1 iteration.

    Validates the bounce-sum + first-leg trajectory machinery against
    the V_α1 algebraic identity :math:`(K \cdot 1)(r,\mu) = \omega_0`:

    - Convergence in 1 iteration.
    - k_eff = k_inf to machine precision (no quadrature error for
      constant trial — the trajectory integrand is itself just a
      constant times :math:`e^{-\Sigma_t s}`, integrable in closed
      form).
    - φ(r) uniform to machine precision.
    """
    fix = fuelA_like_thin_sphere
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    res = solve_greens_function_specular_sphere(
        **fix, n_r=12, n_mu=12, n_traj_quad=24, max_iter=50, tol=1e-12,
    )

    assert res.converged, (
        f"V_α1 numerical: iteration did not converge in 50 iters; "
        f"k_eff = {res.k_eff:.8f}"
    )
    assert res.iterations <= 2, (
        f"V_α1 numerical: should converge in ≤2 iter from const initial "
        f"guess (rank-1 eigenmode); got {res.iterations}"
    )
    np.testing.assert_allclose(
        res.k_eff, k_inf, rtol=1e-12,
        err_msg=(
            f"V_α1 numerical: k_eff = {res.k_eff} ≠ k_inf = {k_inf} "
            f"to machine precision"
        ),
    )
    # φ should be uniform.
    rel_phi_spread = res.phi.std() / res.phi.mean()
    assert rel_phi_spread < 1e-10, (
        f"V_α1 numerical: φ should be uniform (rank-1 eigenmode), "
        f"got std/mean = {rel_phi_spread:.2e}"
    )


@pytest.mark.foundation
def test_v_alpha1_numerical_nonuniform_initial_guess(
    fuelA_like_thin_sphere,
):
    r"""V_α1.numerical — non-uniform ψ_0 converges to rank-1 eigenmode.

    Power iteration starting from a sinusoidal perturbation
    :math:`\psi_0(r, \mu) = 1 + 0.5\,\sin(\pi r/R) \cos(\pi \mu)` should
    converge to the constant rank-1 isotropic eigenmode within ~5–10
    iterations and recover k_eff = k_inf within 0.05 % (the Plan 2 B4
    accuracy target).
    """
    fix = fuelA_like_thin_sphere
    k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])

    n_r, n_mu = 16, 16
    R = fix["R"]

    # Build the same grid the solver uses, so we can craft a matching
    # non-uniform initial guess.
    r_pts, _ = np.polynomial.legendre.leggauss(n_r)
    r_nodes = R * 0.5 * (r_pts + 1.0)
    mu_pts, _ = np.polynomial.legendre.leggauss(n_mu)
    mu_nodes = 0.5 * (mu_pts + 1.0)
    R_grid, MU_grid = np.meshgrid(r_nodes, mu_nodes, indexing="ij")
    psi0 = 1.0 + 0.5 * np.sin(np.pi * R_grid / R) * np.cos(np.pi * MU_grid)

    res = solve_greens_function_specular_sphere(
        **fix, n_r=n_r, n_mu=n_mu, n_traj_quad=32,
        max_iter=100, tol=1e-10, initial_psi=psi0,
    )

    assert res.converged, (
        f"V_α1 nonuniform: iteration did not converge in 100 iters; "
        f"k_eff = {res.k_eff:.8f}, iter = {res.iterations}"
    )

    rel_err = abs(res.k_eff - k_inf) / k_inf
    assert rel_err < 5e-4, (
        f"V_α1 nonuniform: k_eff = {res.k_eff:.8f} differs from "
        f"k_inf = {k_inf:.8f} by {rel_err*100:.4f} %, exceeds 0.05 % target"
    )

    # φ should converge to uniform (rank-1 eigenmode).
    rel_phi_spread = res.phi.std() / res.phi.mean()
    assert rel_phi_spread < 1e-3, (
        f"V_α1 nonuniform: φ should converge to uniform; got "
        f"std/mean = {rel_phi_spread:.2e}"
    )


@pytest.mark.foundation
def test_v_alpha1_numerical_two_thicknesses(
    fuelA_like_thin_sphere, fuelA_like_moderate_sphere,
):
    r"""V_α1.numerical — k_eff = k_inf at two τ_R values.

    Closed sphere has no leakage regardless of thickness. The bounce-
    sum machinery (and the geometric :math:`T(\mu) = 1/(1 - e^{-\Sigma_t
    L_p})` factor) must produce the same k_inf at both thin
    (τ_R = 2.5) and moderate (τ_R = 5) optical depths.
    """
    for fix in (fuelA_like_thin_sphere, fuelA_like_moderate_sphere):
        k_inf = fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])
        res = solve_greens_function_specular_sphere(
            **fix, n_r=12, n_mu=12, n_traj_quad=24,
            max_iter=20, tol=1e-12,
        )
        np.testing.assert_allclose(
            res.k_eff, k_inf, rtol=1e-10,
            err_msg=(
                f"τ_R = {fix['sigma_t']*fix['R']}: k_eff = {res.k_eff} "
                f"≠ k_inf = {k_inf}"
            ),
        )
