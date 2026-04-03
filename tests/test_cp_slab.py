"""Verify the slab collision probability solver against analytical CP eigenvalues."""

import numpy as np
import pytest

from derivations import get
from collision_probability import SlabGeometry, solve_cp_slab


@pytest.mark.parametrize("case_name", [
    "cp_slab_1eg_1rg",
    "cp_slab_1eg_2rg",
    "cp_slab_1eg_4rg",
    "cp_slab_2eg_1rg",
    "cp_slab_2eg_2rg",
    "cp_slab_2eg_4rg",
    "cp_slab_4eg_1rg",
    "cp_slab_4eg_2rg",
    "cp_slab_4eg_4rg",
])
def test_slab_cp_eigenvalue(case_name):
    """Slab CP solver must match the analytical CP eigenvalue."""
    case = get(case_name)
    gp = case.geom_params
    geom = SlabGeometry(
        n_fuel=0, n_clad=0, n_cool=0,
        thicknesses=np.array(gp["thicknesses"]),
        mat_ids=np.array(gp["mat_ids"]),
        N=len(gp["thicknesses"]),
    )
    result = solve_cp_slab(case.materials, geom, keff_tol=1e-7, flux_tol=1e-6)

    err = abs(result.keff - case.k_inf)
    assert err < 1e-6, (
        f"{case_name}: solver={result.keff:.10f} "
        f"analytical={case.k_inf:.10f} err={err:.2e}"
    )
