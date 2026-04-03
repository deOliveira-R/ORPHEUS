"""Verify the cylindrical collision probability solver against analytical CP eigenvalues."""

import numpy as np
import pytest

from derivations import get
from collision_probability import CPGeometry, solve_cp_concentric


@pytest.mark.parametrize("case_name", [
    "cp_cyl1D_1eg_1rg",
    "cp_cyl1D_1eg_2rg",
    "cp_cyl1D_1eg_4rg",
    "cp_cyl1D_2eg_1rg",
    "cp_cyl1D_2eg_2rg",
    "cp_cyl1D_2eg_4rg",
    "cp_cyl1D_4eg_1rg",
    "cp_cyl1D_4eg_2rg",
    "cp_cyl1D_4eg_4rg",
])
def test_cylinder_cp_eigenvalue(case_name):
    """Cylindrical CP solver must match the analytical CP eigenvalue."""
    case = get(case_name)
    gp = case.geom_params
    radii = np.array(gp["radii"])
    mat_ids_list = gp["mat_ids"]

    # Build CPGeometry manually from radii and mat_ids
    r_cell = radii[-1]
    r_inner = np.zeros(len(radii))
    r_inner[1:] = radii[:-1]
    volumes = np.pi * (radii**2 - r_inner**2)

    geom = CPGeometry(
        r_fuel=radii[0],
        r_clad=radii[min(1, len(radii) - 1)],
        r_cell=r_cell,
        n_fuel=0, n_clad=0, n_cool=0,
        radii=radii,
        volumes=volumes,
        mat_ids=np.array(mat_ids_list),
        N=len(radii),
    )
    result = solve_cp_concentric(case.materials, geom)

    err = abs(result.keff - case.k_inf)
    assert err < 1e-5, (
        f"{case_name}: solver={result.keff:.10f} "
        f"analytical={case.k_inf:.10f} err={err:.2e}"
    )
