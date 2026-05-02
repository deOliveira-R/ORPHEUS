r"""Plan-2 follow-on A1 — Variant α prototype with vacuum BC (α=0).

Closed sphere (α=1) is degenerate (rank-1 isotropic eigenmode trivially
gives k_eff=k_inf, V_α1 algebraic identity). Vacuum BC breaks the
degeneracy: spatial leakage produces a non-trivial fundamental mode
(:math:`\phi` peaked near center, vanishing at surface) and
:math:`k_{\rm eff} < k_\infty`. **This case stress-tests the
trajectory machinery** that the closed-sphere tests do not exercise.

Test strategy (L0 sanity checks; the L1 PS-1982 cross-check lives in
``test_peierls_greens_function_xverif_ps1982.py``):

1. **Vacuum k_eff is strictly less than k_inf** (leakage exists).
2. **Increasing R reduces leakage**: k_eff increases monotonically
   with R, asymptoting to k_inf as R → ∞.
3. **Continuity in α**: k_eff(α=0) < k_eff(α=0.5) < k_eff(α=1) for
   any subcritical configuration.
4. **Closed-sphere limit**: at α=1, vacuum-branch reduces to the
   existing rank-1 eigenmode (k_eff = k_inf, machine precision).
5. **Spatial mode is non-trivial for α=0**: φ varies across the
   sphere (std/mean > 0.1), i.e., not the rank-1 isotropic eigenmode.
6. **L0 self-consistency**: vacuum k_eff iteration converges and
   produces positive φ everywhere.

Predecessor:

- :mod:`.test_peierls_greens_function_solver` (closed-sphere V_α1.numerical)
- :mod:`.test_peierls_greens_function_symbolic` (V_α1, V_α2, V_α3 SymPy)

References (see V_α3 SymPy for the algebraic foundation of α=0
reduction; the bounce sum :math:`\alpha B/(1 - \alpha e^{-\Sigma_t L_p})`
collapses to 0 at :math:`\alpha = 0`, so the operator reduces to the
first-leg trajectory integral alone).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.peierls.greens_function import (
    solve_greens_function_sphere,
)


@pytest.fixture(scope="module")
def fuelA_thin_vacuum():
    """fuel-A-like XS, R=5, τ_R=2.5 — substantial leakage expected."""
    return {
        "R": 5.0, "sigma_t": 0.5, "sigma_s": 0.38, "nu_sigma_f": 0.025,
    }


@pytest.fixture(scope="module")
def fuelA_thick_vacuum():
    """fuel-A-like XS, R=20, τ_R=10 — minimal leakage expected."""
    return {
        "R": 20.0, "sigma_t": 0.5, "sigma_s": 0.38, "nu_sigma_f": 0.025,
    }


def _k_inf(fix: dict) -> float:
    return fix["nu_sigma_f"] / (fix["sigma_t"] - fix["sigma_s"])


@pytest.mark.foundation
def test_a1_vacuum_k_eff_less_than_k_inf(fuelA_thin_vacuum):
    """A1.1 — vacuum k_eff < k_inf (leakage exists)."""
    fix = fuelA_thin_vacuum
    res = solve_greens_function_sphere(
        **fix, alpha=0.0, n_r=24, n_mu=24, n_traj_quad=64,
        max_iter=300, tol=1e-8,
    )
    assert res.converged, (
        f"A1.1: vacuum iteration did not converge in 300 iter; "
        f"k_eff = {res.k_eff:.6f}, iter = {res.iterations}"
    )
    k_inf = _k_inf(fix)
    assert res.k_eff < k_inf, (
        f"A1.1: vacuum k_eff = {res.k_eff:.6f} should be < k_inf = "
        f"{k_inf:.6f}"
    )
    # Leakage should be non-trivial for τ_R = 2.5.
    assert res.k_eff < 0.9 * k_inf, (
        f"A1.1: τ_R=2.5 vacuum should have substantial leakage; got "
        f"k_eff/k_inf = {res.k_eff/k_inf:.4f}"
    )


@pytest.mark.foundation
def test_a1_vacuum_thick_sphere_approaches_k_inf(fuelA_thick_vacuum):
    """A1.2 — thicker sphere with vacuum BC has less leakage.

    At τ_R = 10 (10 mean free paths radially) k_eff/k_inf should be
    significantly larger than the τ_R = 2.5 result. For fuel-A-like
    XS at R=20, diffusion estimate is k_eff/k_inf ≈ 1/(1+L²B²) with
    L² = D/Σ_a = 5.56, B² = (π/(R+d))² ≈ 0.0215, giving ratio ≈ 0.89.
    Transport gives a similar value; pin in [0.85, 1.0].
    """
    fix = fuelA_thick_vacuum
    res = solve_greens_function_sphere(
        **fix, alpha=0.0, n_r=24, n_mu=24, n_traj_quad=64,
        max_iter=300, tol=1e-8,
    )
    assert res.converged
    k_inf = _k_inf(fix)
    ratio = res.k_eff / k_inf
    assert 0.85 < ratio < 1.0, (
        f"A1.2: τ_R=10 vacuum should give k_eff/k_inf in [0.85, 1.0); "
        f"got {ratio:.4f}"
    )


@pytest.mark.foundation
def test_a1_alpha_continuity(fuelA_thin_vacuum):
    """A1.3 — k_eff(α) monotonic increasing in α ∈ [0, 1].

    Fewer reflections (smaller α) means more leakage, lower k_eff.
    """
    fix = fuelA_thin_vacuum
    k_eff_by_alpha = {}
    for alpha in (0.0, 0.25, 0.5, 0.75, 1.0):
        res = solve_greens_function_sphere(
            **fix, alpha=alpha, n_r=20, n_mu=20, n_traj_quad=48,
            max_iter=300, tol=1e-8,
        )
        assert res.converged, (
            f"α={alpha}: did not converge ({res.iterations} iter)"
        )
        k_eff_by_alpha[alpha] = res.k_eff

    alphas = sorted(k_eff_by_alpha.keys())
    k_values = [k_eff_by_alpha[a] for a in alphas]

    # Monotonic non-decreasing.
    for i in range(len(k_values) - 1):
        assert k_values[i] <= k_values[i + 1] + 1e-8, (
            f"A1.3: k_eff(α={alphas[i]}) = {k_values[i]:.6f} > "
            f"k_eff(α={alphas[i+1]}) = {k_values[i+1]:.6f} — "
            "should be monotonic non-decreasing in α"
        )

    # Endpoints: α=0 < k_inf, α=1 = k_inf.
    k_inf = _k_inf(fix)
    np.testing.assert_allclose(k_eff_by_alpha[1.0], k_inf, rtol=1e-8)
    assert k_eff_by_alpha[0.0] < k_inf


@pytest.mark.foundation
def test_a1_vacuum_spatial_mode_nontrivial(fuelA_thin_vacuum):
    """A1.4 — vacuum BC produces spatially non-trivial eigenmode.

    For closed sphere the fundamental mode is rank-1 isotropic
    (constant φ); for vacuum the leakage at the surface forces φ → 0
    at r → R, giving a non-trivial spatial profile. The trajectory +
    no-bounce machinery is only stress-tested when this profile is
    non-trivial.
    """
    fix = fuelA_thin_vacuum
    res = solve_greens_function_sphere(
        **fix, alpha=0.0, n_r=24, n_mu=24, n_traj_quad=64,
        max_iter=300, tol=1e-8,
    )
    assert res.converged

    # φ should be positive everywhere (no negative flux artefacts).
    assert (res.phi > 0).all(), (
        f"A1.4: φ should be positive everywhere; got "
        f"min = {res.phi.min():.4e}"
    )

    # Spatial variation: std/mean significantly nonzero. For fuel-A-like
    # at τ_R=2.5 the corrected (backward-leg) trajectory gives a moderate
    # spatial mode with std/mean ≈ 0.45 (peak at r=0 ≈ 5.7 × surface φ).
    rel_spread = res.phi.std() / res.phi.mean()
    assert rel_spread > 0.3, (
        f"A1.4: vacuum φ should be spatially non-trivial; got "
        f"std/mean = {rel_spread:.3f}"
    )

    # Centre-to-surface ratio: should be substantial (vacuum drains
    # surface flux). At τ_R=2.5 expect ratio > 3.
    centre_to_surface = res.phi[0] / res.phi[-1]
    assert centre_to_surface > 3.0, (
        f"A1.4: vacuum φ should be peaked at centre, depleted at surface; "
        f"got φ(centre)/φ(surface) = {centre_to_surface:.3f}"
    )


@pytest.mark.foundation
def test_a1_alpha_one_unchanged_from_specular_branch(fuelA_thin_vacuum):
    """A1.5 — α=1 (default) reproduces the V_α1.numerical result.

    Sanity check that adding the α-parametrisation didn't break the
    closed-sphere case. At α=1 (default) the bounce-sum closure
    :math:`\\alpha B / (1 - \\alpha e^{-\\Sigma_t L_p})` reduces to
    :math:`B / (1 - e^{-\\Sigma_t L_p}) = T(\\mu_{\\rm surf}) \\cdot B`
    (the original specular form), so V_α1.numerical k_eff = k_inf
    must still hold.
    """
    fix = fuelA_thin_vacuum
    k_inf = _k_inf(fix)
    res = solve_greens_function_sphere(
        **fix, alpha=1.0, n_r=12, n_mu=12, n_traj_quad=24,
        max_iter=50, tol=1e-12,
    )
    np.testing.assert_allclose(
        res.k_eff, k_inf, rtol=1e-10,
        err_msg=f"A1.5: α=1 closed sphere should give k_eff=k_inf",
    )
    assert res.iterations <= 2
