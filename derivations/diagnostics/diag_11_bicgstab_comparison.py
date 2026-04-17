"""Diagnostic: Compare sweep vs BiCGSTAB for fixed-source spherical problem.

Created by numerics-investigator on 2026-04-16.
If BiCGSTAB gives flat flux but sweep doesn't, the sweep has a bug
in how it iterates the angular redistribution.
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep import transport_sweep
from orpheus.sn.operator import (
    build_equation_map_spherical,
    build_transport_linear_operator_spherical,
    angular_flux_to_scalar,
    solution_to_angular_flux_spherical,
)
from scipy.sparse.linalg import bicgstab


def test_bicgstab_vs_sweep():
    """Solve the same constant-source problem via BiCGSTAB and compare."""
    nx = 10
    R = 10.0
    N_ord = 4

    mesh = Mesh1D(
        edges=np.linspace(0.0, R, nx + 1),
        mat_ids=np.ones(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC.reflective,
        bc_right=BC.reflective,
    )
    quad = GaussLegendre1D.create(N_ord)
    sn_mesh = SNMesh(mesh, quad)

    N = quad.N
    ng = 1
    sig_t_val = 1.0
    Q_val = 1.0
    W = quad.weights.sum()

    sig_t = np.full((nx, 1, ng), sig_t_val)
    Q_iso_3d = np.full((nx, 1, ng), Q_val)

    # Method 1: Sweep-based source iteration
    psi_bc = {}
    for it in range(200):
        _, phi_sweep = transport_sweep(Q_iso_3d, sig_t, sn_mesh, psi_bc)
        if it > 0:
            res = np.linalg.norm(phi_sweep - phi_old) / max(np.linalg.norm(phi_sweep), 1e-30)
            if res < 1e-14:
                break
        phi_old = phi_sweep.copy()

    phi_sweep_1d = phi_sweep[:, 0, 0]

    # Method 2: BiCGSTAB
    eq_map = build_equation_map_spherical(nx, quad, ng)
    T_op = build_transport_linear_operator_spherical(
        eq_map, quad, sig_t,
        nx, ng,
        sn_mesh.face_areas,
        sn_mesh.volumes,
        sn_mesh.alpha_half,
        sn_mesh.redist_dAw,
        sn_mesh.tau_mm,
    )

    # RHS: isotropic source / sum(w), per ordinate
    # Use the equation map to build the RHS correctly
    rhs = np.zeros((ng, eq_map.n_eq))
    eq_idx = 0
    for ix in range(nx):
        for n in range(N):
            # Skip incoming at outer boundary (reflective BC)
            if ix == nx - 1 and quad.mu_x[n] < -1e-15:
                continue
            rhs[:, eq_idx] = Q_val / W
            eq_idx += 1
    rhs = rhs.ravel(order='F')

    solution, info = bicgstab(T_op, rhs, rtol=1e-14, maxiter=1000)
    if info != 0:
        print(f"  BiCGSTAB did not converge: info={info}")

    fi = solution_to_angular_flux_spherical(solution, eq_map, quad, nx, ng)
    phi_bicgstab = angular_flux_to_scalar(fi, quad, nx, 1, ng)
    phi_bicgstab_1d = phi_bicgstab[:, 0, 0]

    expected = Q_val / sig_t_val

    print(f"\n{'='*70}")
    print(f"BiCGSTAB vs Sweep comparison (nx={nx}, S{N_ord})")
    print(f"Expected phi = {expected:.6f}")
    print(f"\n  {'Cell':>4s}  {'Sweep':>12s}  {'BiCGSTAB':>12s}  {'Expected':>12s}")
    for i in range(nx):
        print(f"  {i:4d}  {phi_sweep_1d[i]:12.8f}  {phi_bicgstab_1d[i]:12.8f}  {expected:12.8f}")

    print(f"\n  Sweep max |err|:    {np.max(np.abs(phi_sweep_1d/expected - 1)):.6e}")
    print(f"  BiCGSTAB max |err|: {np.max(np.abs(phi_bicgstab_1d/expected - 1)):.6e}")
    print(f"{'='*70}")


if __name__ == "__main__":
    test_bicgstab_vs_sweep()
