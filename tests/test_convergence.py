"""Cross-solver consistency: SN vs CP reference values."""

import pytest

from derivations import get
from derivations._xs_library import get_mixture
from sn_1d import GaussLegendreQuadrature, Slab1DGeometry, solve_sn_1d


@pytest.mark.slow
def test_sn_approaches_cp_reference():
    """Fine-mesh SN 1D should approach the CP eigenvalue for slab geometry.

    The gap is due to the white-BC approximation in CP (~1% for moderate
    optical thicknesses). The SN value with reflective BCs is more accurate.
    """
    cp_ref = get("cp_slab_1eg_2rg")

    fuel = get_mixture("A", "1g")
    mod = get_mixture("B", "1g")
    materials = {2: fuel, 0: mod}

    geom = Slab1DGeometry.from_benchmark(
        n_fuel=40, n_mod=40, t_fuel=0.5, t_mod=0.5,
    )
    quad = GaussLegendreQuadrature.gauss_legendre(32)
    result = solve_sn_1d(
        materials, geom, quad,
        max_outer=300, max_inner=500, inner_tol=1e-10,
    )

    gap = abs(result.keff - cp_ref.k_inf)
    assert gap < 0.02, (
        f"SN-CP gap too large: {gap:.4f} "
        f"(SN={result.keff:.8f}, CP={cp_ref.k_inf:.8f})"
    )
