#!/usr/bin/env python3
"""Monte Carlo simulation using the unified geometry scheme.

Demonstrates the MCMesh augmented geometry, which wraps a base Mesh1D
(from :meth:`Mesh1D.from_geometry`) and provides point-wise material
lookup for delta-tracking — the same pattern used by CPMesh and SNMesh.

Two cases:
  1. Cylindrical pin cell via StructuredGeometry.wigner_seitz_pin_cell
  2. Cartesian slab via StructuredGeometry.pwr_slab_half_cell
"""

from pathlib import Path

from orpheus.geometry import Mesh1D, RegionMesh, StructuredGeometry
from orpheus.data.macro_xs.recipes import borated_water, uo2_fuel, zircaloy_clad
from orpheus.mc.solver import MCMesh, MCParams, solve_monte_carlo
from plotting import plot_mc_keff, plot_mc_spectrum

OUTPUT = Path("results_geometry")


def run_case(name: str, mc_mesh, materials, n_neutrons=100,
             n_inactive=50, n_active=500):
    """Run a single MC case and report results."""
    print(f"\n{'=' * 70}")
    print(f"  {name}")
    print(f"{'=' * 70}")
    print(f"  Mesh: {mc_mesh.mesh.N} cells, {mc_mesh.mesh.coord.value}")
    print(f"  Pitch: {mc_mesh.pitch} cm")
    print(f"  Neutrons/cycle: {n_neutrons}")
    print(f"  Cycles: {n_inactive} inactive + {n_active} active\n")

    params = MCParams(
        n_neutrons=n_neutrons,
        n_inactive=n_inactive,
        n_active=n_active,
        geometry=mc_mesh,
    )
    result = solve_monte_carlo(materials, params)

    print(f"\n  keff = {result.keff:.5f} +/- {result.sigma:.5f}")
    return result


def main():
    # Material cross sections
    fuel = uo2_fuel(temp_K=900)
    clad = zircaloy_clad(temp_K=600)
    cool = borated_water(temp_K=600, pressure_MPa=16.0, boron_ppm=4000)
    materials = {2: fuel, 1: clad, 0: cool}

    pitch = 3.6

    # ── Case 1: Cylindrical pin cell (Wigner-Seitz) ─────────────────
    geom_cyl = StructuredGeometry.wigner_seitz_pin_cell(
        r_fuel=0.9, r_clad=1.1, pitch=pitch,
    )
    mesh_cyl = Mesh1D.from_geometry(geom_cyl, region_meshes=(
        RegionMesh(n_cells=10),  # fuel
        RegionMesh(n_cells=3),   # clad
        RegionMesh(n_cells=7),   # cool
    ))
    mc_cyl = MCMesh(mesh_cyl, pitch=pitch)
    result_cyl = run_case("Cylindrical pin cell (MCMesh)", mc_cyl, materials)

    # ── Case 2: Cartesian slab ──────────────────────────────────────
    geom_slab = StructuredGeometry.pwr_slab_half_cell(
        fuel_half=0.9, clad_thick=0.2, cool_thick=0.7,
    )
    mesh_slab = Mesh1D.from_geometry(geom_slab, region_meshes=(
        RegionMesh(n_cells=10),  # fuel
        RegionMesh(n_cells=3),   # clad
        RegionMesh(n_cells=7),   # cool
    ))
    mc_slab = MCMesh(mesh_slab, pitch=pitch)
    result_slab = run_case("Cartesian slab (MCMesh)", mc_slab, materials)

    # ── Save plots (last case wins — both ran successfully) ────────
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plot_mc_keff(result_cyl, OUTPUT)
    plot_mc_spectrum(result_cyl, OUTPUT)
    print(f"\n  Plots saved to {OUTPUT.resolve()}/")


if __name__ == "__main__":
    main()
