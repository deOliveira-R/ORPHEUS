#!/usr/bin/env python3
"""Run the 1D SN transport calculation for a PWR pin cell slab.

True 1D version: same materials and geometry as the 2D demo, but using
Gauss-Legendre quadrature on the 1D slab transport equation.
No pseudo-2D artifacts (no y-direction, no Lebedev quadrature).

Reference MATLAB result (2D Lebedev, 110 ordinates, P0):
    keff = 1.04188
"""

import sys; sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))
import numpy as np

from data.macro_xs.recipes import borated_water, uo2_fuel, zircaloy_clad
from sn_geometry import CartesianMesh
from sn_quadrature import GaussLegendre1D
from sn_solver import solve_sn


def main():
    print("=" * 70)
    print("DISCRETE ORDINATES — PWR PIN CELL (1D SLAB)")
    print("=" * 70)

    # 1. Build per-material macroscopic cross sections (same as 2D demo)
    fuel = uo2_fuel(temp_K=900)
    clad = zircaloy_clad(temp_K=600)
    cool = borated_water(temp_K=600, pressure_MPa=16.0, boron_ppm=4000)
    materials = {2: fuel, 1: clad, 0: cool}

    # 2. Set up 1D slab: [fuel×5, clad×1, cool×4], delta=0.2 cm
    widths = np.full(10, 0.2)
    mat_ids = np.array([2, 2, 2, 2, 2, 1, 0, 0, 0, 0], dtype=int)
    mesh = CartesianMesh.from_slab_1d(widths, mat_ids)

    # 3. S16 Gauss-Legendre quadrature
    n_ord = 16
    gl = GaussLegendre1D.create(n_ord)

    print(f"\n  Slab: {mesh.nx} cells, delta = {widths[0]:.2f} cm")
    print(f"  Ordinates: S{n_ord} ({gl.N} GL points)")
    print(f"  Inner solver: BiCGSTAB (direct transport operator)")
    print(f"  Scattering anisotropy: P0")
    print()

    result = solve_sn(
        materials, mesh, gl,
        inner_solver="bicgstab",
        max_outer=500,
        max_inner=2000,
        inner_tol=1e-4,
    )

    # 4. Report
    print(f"\n  keff = {result.keff:.5f}")
    print(f"  Outer iterations: {len(result.keff_history)}")
    print(f"  Wall time: {result.elapsed_seconds:.1f}s")


if __name__ == "__main__":
    main()
