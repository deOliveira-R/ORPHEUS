"""Verify the slab collision probability solver against analytical CP eigenvalues."""

import pytest

from derivations import get
from collision_probability import SlabGeometry, solve_cp_slab


@pytest.mark.parametrize("case_name", [
    "cp_slab_1eg_2rg",
    "cp_slab_2eg_2rg",
])
def test_slab_cp_eigenvalue(case_name):
    """Slab CP solver must match the analytical CP eigenvalue."""
    case = get(case_name)
    geom = SlabGeometry.default_pwr(**case.geom_params)
    result = solve_cp_slab(case.materials, geom, keff_tol=1e-7, flux_tol=1e-6)

    err = abs(result.keff - case.k_inf)
    assert err < 1e-6, (
        f"{case_name}: solver={result.keff:.10f} "
        f"analytical={case.k_inf:.10f} err={err:.2e}"
    )
