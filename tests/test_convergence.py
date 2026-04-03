"""Cross-solver consistency: SN vs CP reference values."""

import pytest

from derivations import get
from sn_1d import GaussLegendreQuadrature, Slab1DGeometry, solve_sn_1d
from derivations.cp_slab import _XS_A, _XS_B, _make_mixture


@pytest.mark.slow
def test_sn_approaches_cp_reference():
    """Fine-mesh SN 1D should approach the CP eigenvalue for slab geometry.

    The gap is due to the white-BC approximation in CP (~1% for moderate
    optical thicknesses). The SN value with reflective BCs is more accurate.
    """
    cp_ref = get("cp_slab_1eg_2rg")

    fuel = _make_mixture(
        _XS_A["sig_t_1g"], _XS_A["sig_c_1g"],
        _XS_A["sig_f_1g"], _XS_A["nu_1g"],
        _XS_A["chi_1g"], _XS_A["sig_s_1g"],
    )
    mod = _make_mixture(
        _XS_B["sig_t_1g"], _XS_B["sig_c_1g"],
        _XS_B["sig_f_1g"], _XS_B["nu_1g"],
        _XS_B["chi_1g"], _XS_B["sig_s_1g"],
    )
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
