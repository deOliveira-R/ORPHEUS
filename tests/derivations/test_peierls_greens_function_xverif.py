r"""Plan-2 B5 cross-verification — Variant α vs Phase 4 specular_multibounce
vs white_hebert on the closed homogeneous fuel-A sphere.

The load-bearing claim of Plan 2 is that **Variant α gives the exact
:math:`k_{\rm eff} = k_\infty` for the closed homogeneous sphere with
specular BC** (V_α1 algebraic identity), while Phase 4
``specular_multibounce`` and ``white_hebert`` carry a small rank-N
truncation error that converges as N grows.

This module pins:

1. **Variant α produces k_inf to machine precision** (V_α1.numerical
   recap, also in :mod:`.test_peierls_greens_function_solver`).
2. **Phase 4 specular_multibounce at rank-1 ≡ white_hebert** (already
   pinned by :func:`test_specular_multibounce_rank1_equals_hebert`,
   recapped here as part of the cross-check matrix).
3. **Phase 4 specular_multibounce converges to Variant α as rank grows**:
   :math:`|k_{4,N} - k_\alpha| / k_\alpha` decreases monotonically (or
   approximately so) as :math:`N \in \{1, 2, 3\}`.

Reference values (fuel-A-like 1G, σ_t = 0.5, σ_s = 0.38, νσ_f = 0.025,
R = 5):

- :math:`k_\infty = 0.025 / 0.12 = 0.2083\overline{3}`
- Variant α k_eff = k_inf to ≤ 1e-10 (V_α1 numerical)
- Phase 4 N = 1 ≡ white_hebert ≈ k_inf · (1 ± 0.5 %) (existing pin)
- Phase 4 N = 3 ≈ k_inf · (1 ± 0.2 %) (existing pin)

Predecessor:

- :mod:`.test_peierls_greens_function_solver` (V_α1.numerical gates)
- :mod:`.test_peierls_specular_bc` (Phase 4 closure tests)
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.peierls.geometry import (
    SPHERE_1D,
    solve_peierls_1g,
)
from orpheus.derivations.continuous.peierls.greens_function import (
    solve_greens_function_sphere,
)


@pytest.fixture(scope="module")
def fuelA_thin_sphere_1G():
    sig_t = np.array([0.5])
    sig_s = np.array([[0.38]])
    nu_sig_f = np.array([0.025])
    radii = np.array([5.0])
    k_inf = float(nu_sig_f[0]) / (float(sig_t[0]) - float(sig_s[0, 0]))
    return {
        "sig_t": sig_t, "sig_s": sig_s, "nu_sig_f": nu_sig_f,
        "radii": radii, "k_inf": k_inf,
    }


@pytest.mark.foundation
def test_b5_variant_alpha_gives_k_inf_exactly(fuelA_thin_sphere_1G):
    r"""B5.A — Variant α gives :math:`k_{\rm eff} = k_\infty` to
    machine precision (V_α1.numerical recap).

    This is the **truth source** for the cross-verification matrix.
    All other closures (white_hebert, specular_multibounce at any
    rank) are approximations of this exact result and should approach
    it as their rank or quadrature parameters tighten.
    """
    fix = fuelA_thin_sphere_1G
    res = solve_greens_function_sphere(
        R=float(fix["radii"][0]),
        sigma_t=float(fix["sig_t"][0]),
        sigma_s=float(fix["sig_s"][0, 0]),
        nu_sigma_f=float(fix["nu_sig_f"][0]),
        n_r=16, n_mu=16, n_traj_quad=32,
        max_iter=50, tol=1e-12,
    )
    np.testing.assert_allclose(
        res.k_eff, fix["k_inf"], rtol=1e-10,
        err_msg=(
            f"B5.A: Variant α k_eff = {res.k_eff} ≠ k_inf = "
            f"{fix['k_inf']} to 1e-10 tolerance"
        ),
    )


@pytest.mark.foundation
def test_b5_phase4_rank1_equals_white_hebert(fuelA_thin_sphere_1G):
    r"""B5.B — Phase 4 specular_multibounce at rank-1 equals white_hebert
    bit-for-bit (recap of existing pin).

    Algebraic identity: at :math:`N = 1` the construction
    :math:`K_{\rm bc}^{\rm spec,mb} = G\,R\,(I - T R)^{-1}\,P` collapses
    to :math:`G\,(1/(1-P_{ss}))\,P` since :math:`R = [[1]]` and
    :math:`T_{00} = P_{ss}` (V_α2). Recapped here so the B5
    cross-verification matrix is self-contained.
    """
    fix = fuelA_thin_sphere_1G
    sol_mb = solve_peierls_1g(
        SPHERE_1D, fix["radii"], fix["sig_t"], fix["sig_s"],
        fix["nu_sig_f"], boundary="specular_multibounce", n_bc_modes=1,
        p_order=4, n_panels_per_region=2,
        n_angular=24, n_rho=24, n_surf_quad=24, dps=20,
        tol=1e-10,
    )
    sol_heb = solve_peierls_1g(
        SPHERE_1D, fix["radii"], fix["sig_t"], fix["sig_s"],
        fix["nu_sig_f"], boundary="white_hebert", n_bc_modes=1,
        p_order=4, n_panels_per_region=2,
        n_angular=24, n_rho=24, n_surf_quad=24, dps=20,
        tol=1e-10,
    )
    np.testing.assert_allclose(
        sol_mb.k_eff, sol_heb.k_eff, rtol=1e-8, atol=1e-10,
        err_msg=(
            f"B5.B: Phase 4 N=1 = {sol_mb.k_eff}, white_hebert = "
            f"{sol_heb.k_eff} — should be bit-equal"
        ),
    )


@pytest.mark.foundation
def test_b5_phase4_converges_toward_variant_alpha(fuelA_thin_sphere_1G):
    r"""B5.C — Phase 4 specular_multibounce error wrt Variant α
    decreases as rank N grows.

    Variant α gives :math:`k_\alpha = k_\infty` exactly. Phase 4 at
    rank N has error :math:`\epsilon_N = |k_{4,N} - k_\alpha|/k_\alpha`.
    For thin sphere (τ_R = 2.5):

    - :math:`\epsilon_1 \lesssim 0.5 \%` (rank-1 ≡ white_hebert)
    - :math:`\epsilon_3 \lesssim 0.2 \%` (rank-3 closer to convergence)
    - :math:`\epsilon_3 < \epsilon_1` (rank convergence)

    The pinned ratios are the documented Phase 4 accuracy figures
    against the V_α1-exact reference. Variant α therefore provides a
    cleaner reference than Phase 4 specular_multibounce for closed-
    sphere homogeneous specular work.
    """
    fix = fuelA_thin_sphere_1G

    # Variant α reference (exact = k_inf).
    res_alpha = solve_greens_function_sphere(
        R=float(fix["radii"][0]),
        sigma_t=float(fix["sig_t"][0]),
        sigma_s=float(fix["sig_s"][0, 0]),
        nu_sigma_f=float(fix["nu_sig_f"][0]),
        n_r=16, n_mu=16, n_traj_quad=32,
        max_iter=50, tol=1e-12,
    )
    k_alpha = res_alpha.k_eff

    # Phase 4 at rank N=1, 2, 3.
    sol = {}
    for N in (1, 2, 3):
        sol[N] = solve_peierls_1g(
            SPHERE_1D, fix["radii"], fix["sig_t"], fix["sig_s"],
            fix["nu_sig_f"], boundary="specular_multibounce", n_bc_modes=N,
            p_order=4, n_panels_per_region=2,
            n_angular=24, n_rho=24, n_surf_quad=24, dps=20,
            tol=1e-10,
        )

    eps_1 = abs(sol[1].k_eff - k_alpha) / k_alpha
    eps_2 = abs(sol[2].k_eff - k_alpha) / k_alpha
    eps_3 = abs(sol[3].k_eff - k_alpha) / k_alpha

    # Pin the absolute error figures against the V_α1-exact reference.
    assert eps_1 < 5e-3, (
        f"B5.C: Phase 4 N=1 error wrt Variant α = {eps_1*100:.4f} %, "
        f"exceeds 0.5 % expectation"
    )
    assert eps_3 < 2e-3, (
        f"B5.C: Phase 4 N=3 error wrt Variant α = {eps_3*100:.4f} %, "
        f"exceeds 0.2 % expectation"
    )

    # Pin the rank-convergence direction: N=3 should be closer to
    # Variant α than N=1.
    assert eps_3 <= eps_1, (
        f"B5.C: Phase 4 should converge toward Variant α as rank grows; "
        f"got eps_1 = {eps_1*100:.4f} % vs eps_3 = {eps_3*100:.4f} % "
        f"(rank-3 not closer than rank-1)"
    )
