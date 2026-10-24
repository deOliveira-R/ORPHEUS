#!/usr/bin/env python3
"""Monte Carlo simulation using the unified geometry scheme.

Demonstrates the MCMesh augmented geometry, which wraps a base Mesh1D
(built by :class:`~orpheus.mesh.mesher.Mesher`) and provides point-wise material
lookup for delta-tracking — the same pattern used by CPMesh and SNProblem.

Two cases:
  1. Cylindrical pin cell via StructuredGeometry.wigner_seitz_pin_cell: MC REFUSES the specular wall declared on a cylinder
     (``Refused(TOPOLOGY, ..., #547)``): what MC computes there is a periodic equal-area square lattice cell (a torus wrap), which
     no one-dimensional declaration describes. The case shows the refusal.
  2. Cartesian slab, a translation deck on each face (the periodic wrap MC computes), built with StructuredGeometry.from_thicknesses.
"""

from pathlib import Path

from orpheus.geometry import CoordSystem, StructuredGeometry, TranslationDeck
from orpheus.mesh import CellsByCount, Mesher
from orpheus.data.macro_xs.recipes import borated_water, uo2_fuel, zircaloy_clad
from orpheus.mc.solver import MCMesh, MCParams, solve_monte_carlo
from orpheus.transport.boundary_conditions import BoundaryRefused, WallCondition
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
    mesh_cyl = Mesher(geom_cyl).partition((
        CellsByCount.uniform_volume(10),  # fuel
        CellsByCount.uniform_volume(3),   # clad
        CellsByCount.uniform_volume(7),   # cool
    )).mesh
    try:
        # the specular wall is the declaration MoC realizes on this cell; MC refuses it (#547)
        MCMesh(mesh_cyl, pitch=pitch, boundary_conditions={0: WallCondition.specular(1.0)})
    except BoundaryRefused as refusal:
        print(f"\n  Cylindrical pin cell (MCMesh): refused, {refusal.admission.reason} (#{refusal.admission.issue})")

    # ── Case 2: Cartesian slab ──────────────────────────────────────
    geom_slab = StructuredGeometry.from_thicknesses(
        coord=CoordSystem.CARTESIAN, thicknesses=(0.9, 0.2, 0.7), mat_ids=(2, 1, 0),
        face_roles=(TranslationDeck(), TranslationDeck()),
    )
    mesh_slab = Mesher(geom_slab).partition((
        CellsByCount.uniform_volume(10),  # fuel
        CellsByCount.uniform_volume(3),   # clad
        CellsByCount.uniform_volume(7),   # cool
    )).mesh
    mc_slab = MCMesh(mesh_slab, pitch=pitch)
    result_slab = run_case("Cartesian slab (MCMesh)", mc_slab, materials)

    # ── Save plots (the slab case: the cylinder case is refused) ───
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plot_mc_keff(result_slab, OUTPUT)
    plot_mc_spectrum(result_slab, OUTPUT)
    print(f"\n  Plots saved to {OUTPUT.resolve()}/")


if __name__ == "__main__":
    main()
