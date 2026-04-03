"""Verify the Method of Characteristics solver against homogeneous eigenvalues."""

import pytest

from derivations import get
from method_of_characteristics import MoCGeometry, solve_moc


@pytest.mark.parametrize("case_name", [
    "moc_cyl1D_1eg_1rg",
    "moc_cyl1D_2eg_1rg",
])
def test_moc_homogeneous(case_name):
    """MOC with homogeneous fill must match infinite-medium eigenvalue."""
    case = get(case_name)
    mix = next(iter(case.materials.values()))
    materials = {0: mix, 1: mix, 2: mix}
    geom = MoCGeometry.default_pwr()
    result = solve_moc(materials, geom, max_outer=200)

    err = abs(result.keff - case.k_inf)
    assert err < 1e-4, (
        f"{case_name}: solver={result.keff:.6f} "
        f"analytical={case.k_inf:.6f} err={err:.2e}"
    )
