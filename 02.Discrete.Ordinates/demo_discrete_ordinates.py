#!/usr/bin/env python3
"""Discrete Ordinates (SN) transport calculation for a PWR pin cell slab.

Solves the 1D multi-group neutron transport equation with reflective
boundary conditions using S16 Gauss-Legendre quadrature, P0 scattering,
and BiCGSTAB inner iterations (direct transport operator).

Geometry: 10-cell slab (fuel × 5 + clad × 1 + coolant × 4), δ = 0.2 cm.
"""

import sys; sys.path.insert(0, str(__import__('pathlib').Path(__file__).resolve().parent.parent))
from pathlib import Path

import numpy as np

from data.macro_xs.recipes import borated_water, uo2_fuel, zircaloy_clad
from sn_geometry import CartesianMesh
from sn_quadrature import GaussLegendre1D
from sn_solver import solve_sn
from plotting import plot_do_convergence, plot_do_spectra

OUTPUT = Path("results")


def plot_spatial_flux_1d(result, output_dir):
    """Plot thermal / resonance / fast flux along the slab."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mesh = result.geometry
    x = np.cumsum(mesh.dx) - mesh.dx / 2  # cell centres

    sf = result.scalar_flux[:, 0, :]  # (nx, ng), squeeze ny=1
    ng = sf.shape[1]

    FI_T = sf[:, :50].sum(axis=1) if ng > 50 else sf.sum(axis=1)
    FI_R = sf[:, 50:287].sum(axis=1) if ng > 287 else np.zeros(mesh.nx)
    FI_F = sf[:, 287:].sum(axis=1) if ng > 287 else np.zeros(mesh.nx)

    fig, ax = plt.subplots()
    ax.plot(x, FI_F, "-or", label="Fast", markersize=4)
    ax.plot(x, FI_R, "-og", label="Resonance", markersize=4)
    ax.plot(x, FI_T, "-ob", label="Thermal", markersize=4)
    ax.set_xlabel("Distance from cell centre (cm)")
    ax.set_ylabel("Neutron flux (a.u.)")
    ax.legend(loc="upper left")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(output_dir / "DO_05_flux_cell.pdf")
    plt.close(fig)


def main():
    print("=" * 70)
    print("DISCRETE ORDINATES — PWR PIN CELL (1D SLAB)")
    print("=" * 70)

    # 1. Build per-material macroscopic cross sections
    fuel = uo2_fuel(temp_K=900)
    clad = zircaloy_clad(temp_K=600)
    cool = borated_water(temp_K=600, pressure_MPa=16.0, boron_ppm=4000)
    materials = {2: fuel, 1: clad, 0: cool}

    # 2. Set up 1D slab geometry: fuel(5) + clad(1) + coolant(4), δ=0.2 cm
    widths = np.full(10, 0.2)
    mat_ids = np.array([2, 2, 2, 2, 2, 1, 0, 0, 0, 0], dtype=int)
    mesh = CartesianMesh.from_slab_1d(widths, mat_ids)

    # 3. Angular quadrature
    n_ord = 16
    quad = GaussLegendre1D.create(n_ord)

    print(f"\n  Slab: {mesh.nx} cells, delta = {widths[0]:.2f} cm")
    print(f"  Layout: fuel(5) + clad(1) + coolant(4)")
    print(f"  Ordinates: S{n_ord} ({quad.N} Gauss-Legendre points)")
    print(f"  Inner solver: BiCGSTAB (direct transport operator)")
    print(f"  Scattering anisotropy: P0")
    print()

    # 4. Solve
    result = solve_sn(
        materials, mesh, quad,
        inner_solver="bicgstab",
        max_outer=500,
        max_inner=2000,
        inner_tol=1e-4,
    )

    # 5. Report
    print(f"\n  keff = {result.keff:.5f}")
    print(f"  Outer iterations: {len(result.keff_history)}")
    print(f"  Wall time: {result.elapsed_seconds:.1f}s")

    # 6. Plots
    OUTPUT.mkdir(parents=True, exist_ok=True)
    plot_do_convergence(result, OUTPUT)
    plot_do_spectra(result, materials, OUTPUT)
    plot_spatial_flux_1d(result, OUTPUT)
    print(f"\n  Plots saved to {OUTPUT.resolve()}/")


if __name__ == "__main__":
    main()
