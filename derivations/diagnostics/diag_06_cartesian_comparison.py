"""Diagnostic: Compare spherical vs cartesian sweep for constant source.

Created by numerics-investigator on 2026-04-16.
If cartesian gives flat flux but spherical doesn't, the bug is in the
curvilinear-specific parts (redistribution or geometry factors).
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep import transport_sweep


def run_sweep_test(geometry, nx=10, N_ord=4, R=5.0, max_iter=100, tol=1e-14):
    """Run converging sweeps for a given geometry."""
    if geometry == "spherical":
        coord = CoordSystem.SPHERICAL
    elif geometry == "cylindrical":
        coord = CoordSystem.CYLINDRICAL
    else:
        coord = CoordSystem.CARTESIAN

    mesh = Mesh1D(
        edges=np.linspace(0.0 if geometry != "cartesian" else 0.0, R, nx + 1),
        mat_ids=np.ones(nx, dtype=int),
        coord=coord,
        bc_left=BC.reflective,
        bc_right=BC.reflective,
    )
    quad = GaussLegendre1D.create(N_ord)
    sn_mesh = SNMesh(mesh, quad)

    ng = 1
    sig_t = np.full((nx, 1, ng), 1.0)
    Q_iso = np.full((nx, 1, ng), 1.0)
    psi_bc = {}

    for it in range(max_iter):
        _, phi = transport_sweep(Q_iso, sig_t, sn_mesh, psi_bc)
        if it > 0:
            res = np.linalg.norm(phi - phi_old) / max(np.linalg.norm(phi), 1e-30)
            if res < tol:
                break
        phi_old = phi.copy()

    return phi[:, 0, 0], it + 1


def test_geometry_comparison():
    """Compare all geometries."""
    nx = 10
    N_ord = 4
    R = 5.0
    expected = 1.0  # Q / Sig_t

    print(f"\n{'='*70}")
    print(f"Geometry comparison: constant source, reflective BCs")
    print(f"Expected phi = {expected:.6f}")

    for geom in ["cartesian", "cylindrical", "spherical"]:
        try:
            phi, n_iter = run_sweep_test(geom, nx, N_ord, R)
            print(f"\n  {geom:12s}: converged in {n_iter} iters")
            print(f"    phi range: [{phi.min():.8f}, {phi.max():.8f}]")
            print(f"    phi mean:  {phi.mean():.8f}")
            print(f"    max |err|: {np.max(np.abs(phi/expected - 1)):.6e}")
        except Exception as e:
            print(f"\n  {geom:12s}: FAILED — {e}")

    print(f"{'='*70}")


if __name__ == "__main__":
    test_geometry_comparison()
