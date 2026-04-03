"""Verify the cylindrical collision probability solver against analytical CP eigenvalues."""

import pytest

from derivations import get
from collision_probability import CPGeometry, solve_cp_concentric


@pytest.mark.parametrize("case_name", [
    "cp_cyl1D_1eg_2rg",
    "cp_cyl1D_2eg_2rg",
])
def test_cylinder_cp_eigenvalue(case_name):
    """Cylindrical CP solver must match the analytical CP eigenvalue."""
    case = get(case_name)
    geom = CPGeometry.default_pwr(**case.geom_params)
    result = solve_cp_concentric(case.materials, geom)

    err = abs(result.keff - case.k_inf)
    assert err < 1e-5, (
        f"{case_name}: solver={result.keff:.10f} "
        f"analytical={case.k_inf:.10f} err={err:.2e}"
    )
