"""Verify the Monte Carlo solver against homogeneous eigenvalues (statistical)."""

import pytest

from derivations import get
from monte_carlo import MCParams, solve_monte_carlo


@pytest.mark.parametrize("case_name,n_active", [
    ("mc_cyl1D_1eg_1rg", 500),
    ("mc_cyl1D_2eg_1rg", 200),
])
def test_mc_zscore(case_name, n_active):
    """MC eigenvalue must be within 5 sigma of the analytical reference."""
    case = get(case_name)
    mix = next(iter(case.materials.values()))
    materials = {0: mix, 1: mix, 2: mix}

    params = MCParams(
        n_neutrons=200, n_inactive=50, n_active=n_active, seed=42,
    )
    result = solve_monte_carlo(materials, params)
    z_score = abs(result.keff - case.k_inf) / max(result.sigma, 1e-10)

    assert z_score < 5.0, (
        f"{case_name}: k_mc={result.keff:.6f} +/- {result.sigma:.5f}, "
        f"k_ref={case.k_inf:.6f}, z={z_score:.2f}"
    )


@pytest.mark.slow
@pytest.mark.parametrize("case_name", ["mc_cyl1D_1eg_1rg"])
def test_mc_high_stats(case_name):
    """With more histories, MC should be within 3 sigma."""
    case = get(case_name)
    mix = next(iter(case.materials.values()))
    materials = {0: mix, 1: mix, 2: mix}

    params = MCParams(
        n_neutrons=500, n_inactive=100, n_active=2000, seed=42,
    )
    result = solve_monte_carlo(materials, params)
    z_score = abs(result.keff - case.k_inf) / max(result.sigma, 1e-10)

    assert z_score < 3.0, (
        f"High-stats {case_name}: k_mc={result.keff:.6f} +/- {result.sigma:.5f}, "
        f"k_ref={case.k_inf:.6f}, z={z_score:.2f}"
    )
