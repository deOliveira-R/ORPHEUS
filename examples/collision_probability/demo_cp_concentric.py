#!/usr/bin/env python3
"""Run the collision probability transport calculation for a PWR pin cell.

Uses the same materials as the Discrete Ordinates (SN) exercise to allow
direct comparison of eigenvalues.

Supports both Jacobi (default) and Gauss-Seidel solver modes.
Pass --gauss-seidel to use the GS group sweep with inner iterations.

Reference results:
    SN (slab geometry):   keff = 1.04188
    MC (slab geometry):   keff = 1.03484 +/- 0.00192
"""

import argparse
from pathlib import Path

from orpheus.data.macro_xs.recipes import borated_water, uo2_fuel, zircaloy_clad
from orpheus.geometry import Mesh1D, RegionMesh, StructuredGeometry
from orpheus.cp.solver import CPParams, solve_cp
from plotting import (
    plot_cp_convergence,
    plot_cp_geometry,
    plot_cp_inner_iterations,
    plot_cp_radial_flux,
    plot_cp_spectra,
)

OUTPUT = Path("results")


def main():
    parser = argparse.ArgumentParser(description="CP solver for PWR pin cell")
    parser.add_argument("--gauss-seidel", action="store_true",
                        help="Use Gauss-Seidel group sweep with inner iterations")
    args = parser.parse_args()

    solver_mode = "gauss_seidel" if args.gauss_seidel else "jacobi"

    print("=" * 70)
    print(f"COLLISION PROBABILITY — PWR PIN CELL (mode: {solver_mode})")
    print("=" * 70)

    # 1. Build per-material macroscopic cross sections (same as SN)
    fuel = uo2_fuel(temp_K=900)
    clad = zircaloy_clad(temp_K=600)
    cool = borated_water(temp_K=600, pressure_MPa=16.0, boron_ppm=4000)
    materials = {2: fuel, 1: clad, 0: cool}

    # 2. Set up Wigner-Seitz cylindrical geometry, then discretize.
    #    Geometry layer: shape + BCs only. Mesh layer: per-region cell
    #    counts and discretization scheme.
    geom = StructuredGeometry.wigner_seitz_pin_cell(
        r_fuel=0.9, r_clad=1.1, pitch=3.6,
    )
    mesh = Mesh1D.from_geometry(geom, region_meshes=(
        RegionMesh(n_cells=10),  # fuel — equal-volume default
        RegionMesh(n_cells=3),   # clad
        RegionMesh(n_cells=7),   # coolant
    ))
    params = CPParams(solver_mode=solver_mode)

    n_fuel = (mesh.mat_ids == 2).sum()
    n_clad = (mesh.mat_ids == 1).sum()
    n_cool = (mesh.mat_ids == 0).sum()

    print(f"\n  Geometry: r_cell = {mesh.edges[-1]:.3f} cm")
    print(f"  Sub-regions: {n_fuel} fuel + {n_clad} clad "
          f"+ {n_cool} cool = {mesh.N} total")
    print()

    # 3. Solve
    result = solve_cp(materials, mesh, params)

    # 4. Report
    print(f"\n  keff = {result.keff:.5f}  (SN slab reference: 1.04188)")
    print(f"  Outer iterations: {len(result.keff_history)}")
    if result.residual_history:
        print(f"  Final balance residual: {result.residual_history[-1]:.2e}")
    # #340 N4: the iteration tree replaced the retired ``n_inner`` array.
    # Each inner child names its own group, so the worst group needs no axis
    # convention -- and the record also says whether it CONVERGED.
    inners = list(result.record.children)
    if inners:
        ng = len({child.label for child in inners})
        last_outer = inners[-ng:]
        counts = [child.n_iterations for child in last_outer]
        print(f"  Inner iterations (last outer): "
              f"max={max(counts)}, mean={sum(counts) / len(counts):.1f}")
        worst = max(last_outer, key=lambda child: child.n_iterations)
        print(f"  Most inner iterations in {worst.label} "
              f"({worst.n_iterations} iters, {worst.status})")
    failure = result.record.first_failure
    if failure is not None:
        # ``covering`` converts the projection (in the level's own iteration
        # units) into the units the knob takes.  Identity for CP, whose
        # ``max_inner`` is one within-group pass; it is the Krylov arm where
        # the two differ (#349).
        needed = failure.projected_iterations()
        print(f"  [!] NOT fully converged -- {failure.label} "
              f"{failure.status}; set {failure.budget.name}="
              f"{failure.budget.covering(needed) if needed else '?'}")
    print(f"  Wall time: {result.elapsed_seconds:.1f}s")

    # 5. Plots
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plot_cp_geometry(mesh, OUTPUT)
    plot_cp_convergence(result, OUTPUT)
    plot_cp_spectra(result, OUTPUT)
    plot_cp_radial_flux(result, OUTPUT)
    plot_cp_inner_iterations(result, OUTPUT)
    print(f"\n  Plots saved to {OUTPUT.resolve()}/")


if __name__ == "__main__":
    main()