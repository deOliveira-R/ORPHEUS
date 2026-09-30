"""Verify the spherical collision probability solver against analytical CP eigenvalues."""

import pytest

from orpheus.derivations import get
from orpheus.geometry import BC, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesher
from orpheus.cp.solver import solve_cp

pytestmark = [pytest.mark.l1, pytest.mark.verifies(
    "collision-rate",
    "chord-length",
    "self-sph",
    "second-diff-sph",
    "one-group-kinf",
    "matrix-eigenvalue",
    "mg-balance",
)]


@pytest.mark.parametrize("case_name", [
    "cp_sph1D_1eg_1rg",
    "cp_sph1D_1eg_2rg",
    "cp_sph1D_1eg_4rg",
    "cp_sph1D_2eg_1rg",
    "cp_sph1D_2eg_2rg",
    "cp_sph1D_2eg_4rg",
    "cp_sph1D_4eg_1rg",
    "cp_sph1D_4eg_2rg",
    "cp_sph1D_4eg_4rg",
])
def test_sphere_cp_eigenvalue(case_name):
    """Spherical CP solver must match the analytical CP eigenvalue."""
    case = get(case_name)
    gp = case.geom_params
    geometry = StructuredGeometry.sphere(
        (0.0, *gp["radii"]), gp["mat_ids"], outer=BC.white,
    )
    mesh = Mesher(geometry).partition(CellsByCount.uniform_width(1)).mesh
    result = solve_cp(case.materials, mesh)

    err = abs(result.keff - case.k_inf)
    assert err < 1e-5, (
        f"{case_name}: solver={result.keff:.10f} "
        f"analytical={case.k_inf:.10f} err={err:.2e}"
    )
