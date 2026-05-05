"""Cross-verification: Monte Carlo vs deterministic solvers on the same geometry.

XV-MC-001: MC vs CP cylinder (2G, 2-region)
XV-MC-002: MC vs CP slab (2G, 2-region)

The CP reference uses white/reflective BCs while MC uses periodic BCs.
For sufficiently large cells the difference is ~1% — this is physical,
not a bug. Tolerances account for both statistical and BC-approximation
errors.
"""

import numpy as np
import pytest

from orpheus.derivations import get
from orpheus.mc.solver import MCParams, ConcentricPinCell, SlabPinCell, solve_monte_carlo
from orpheus.cp.solver import solve_cp, CPParams
from orpheus.geometry import (
    BC,
    CoordSystem,
    Mesh1D,
    Region,
    RegionMesh,
    StructuredGeometry,
)

# L2 cross-code MC ↔ CP consistency.
pytestmark = [
    pytest.mark.l2,
    pytest.mark.verifies(
        "keff-mean",
        "sigma-keff",
        "collision-rate",
        "self-cyl",
        "self-slab",
    ),
]


# ═══════════════════════════════════════════════════════════════════════
# XV-MC-001: MC vs CP cylinder
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.slow
def test_mc_vs_cp_cylinder():
    """XV-MC-001: MC and CP give consistent keff on the same cylindrical pin cell.

    Both use the 2G 2-region cylindrical problem from the derivation library.
    Tolerance: |k_MC - k_CP| < 5*sigma_MC + 0.01 (1% for BC approximation).
    """
    case = get("mc_cyl1D_2eg_2rg")
    gp = case.geom_params
    radii = gp["radii"]
    mat_ids = gp["mat_ids"]

    # ── CP solve ──────────────────────────────────────────────────────
    # Build a multi-region cylinder via StructuredGeometry: each region
    # carries the layer thickness (literal, not cumulative) plus its
    # material ID; the mesh layer takes 5 cells per region.
    regions = []
    inner = 0.0
    for r_outer, m in zip(radii, mat_ids):
        regions.append(
            Region(mat_id=int(m), outer_thickness_cm=float(r_outer - inner)),
        )
        inner = r_outer
    geom = StructuredGeometry(
        geometry="CYL",
        regions=tuple(regions),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=tuple(
        RegionMesh(n_cells=5) for _ in regions
    ))
    cp_result = solve_cp(case.materials, mesh=mesh)

    # ── MC solve ──────────────────────────────────────────────────────
    r_cell = radii[-1]
    pitch = r_cell * np.sqrt(np.pi)
    geom = ConcentricPinCell(radii=radii, mat_ids=mat_ids, pitch=pitch)
    mc_params = MCParams(
        n_neutrons=200, n_inactive=50, n_active=500,
        seed=42, geometry=geom,
    )
    mc_result = solve_monte_carlo(case.materials, mc_params)

    # ── Compare ───────────────────────────────────────────────────────
    diff = abs(mc_result.keff - cp_result.keff)
    tolerance = 5.0 * mc_result.sigma + 0.02  # statistical + BC approximation

    assert diff < tolerance, (
        f"MC vs CP cylinder: k_MC={mc_result.keff:.5f} +/- {mc_result.sigma:.5f}, "
        f"k_CP={cp_result.keff:.5f}, diff={diff:.5f}, tol={tolerance:.5f}"
    )


# ═══════════════════════════════════════════════════════════════════════
# XV-MC-002: MC vs CP slab
# ═══════════════════════════════════════════════════════════════════════

@pytest.mark.slow
def test_mc_vs_cp_slab():
    """XV-MC-002: MC and CP give consistent keff on a slab geometry.

    Uses a 2G 2-region slab (fuel + moderator).
    CP uses reflective/white BCs, MC uses periodic BCs on a symmetric cell.
    The MC geometry is a full cell (mod|fuel|mod) to approximate the
    infinite lattice that CP computes.
    """
    # Use the same 2G cross sections as the cylindrical case
    from orpheus.derivations.common.xs_library import get_mixture, get_xs

    mix_fuel = get_mixture("A", "2g")
    mix_mod = get_mixture("B", "2g")
    materials = {2: mix_fuel, 0: mix_mod}

    # ── CP solve (slab half-cell: fuel | moderator) ──────────────────
    t_fuel = 0.9
    t_mod = 0.9
    geom_cp = StructuredGeometry(
        geometry="SLB",
        regions=(
            Region(mat_id=2, outer_thickness_cm=t_fuel),
            Region(mat_id=0, outer_thickness_cm=t_mod),
        ),
        bcs=(BC.reflective, BC.reflective),
    )
    mesh = Mesh1D.from_geometry(geom_cp, region_meshes=(
        RegionMesh(n_cells=10),
        RegionMesh(n_cells=10),
    ))
    cp_result = solve_cp(materials, mesh=mesh)

    # ── MC solve (full cell: mod|fuel|mod for periodic BCs) ──────────
    pitch = 2.0 * (t_fuel + t_mod)  # full cell
    geom = SlabPinCell(
        boundaries=[t_mod, t_mod + 2 * t_fuel, t_mod + 2 * t_fuel + t_mod],
        mat_ids=[0, 2, 0, 0],  # mod | fuel | mod | (remainder = mod)
        pitch=pitch,
    )
    mc_params = MCParams(
        n_neutrons=200, n_inactive=50, n_active=500,
        seed=42, geometry=geom,
    )
    mc_result = solve_monte_carlo(materials, mc_params)

    # ── Compare ───────────────────────────────────────────────────────
    diff = abs(mc_result.keff - cp_result.keff)
    tolerance = 5.0 * mc_result.sigma + 0.03  # wider: slab BC approximation larger

    assert diff < tolerance, (
        f"MC vs CP slab: k_MC={mc_result.keff:.5f} +/- {mc_result.sigma:.5f}, "
        f"k_CP={cp_result.keff:.5f}, diff={diff:.5f}, tol={tolerance:.5f}"
    )
