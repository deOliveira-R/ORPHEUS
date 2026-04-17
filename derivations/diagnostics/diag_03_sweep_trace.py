"""Diagnostic: Trace a single sweep with constant source, check cell 0.

Created by numerics-investigator on 2026-04-16.
The sweep iterates to a converged solution that is NOT flat, despite
flat flux being the analytical solution. Trace one sweep to understand
why the iteration diverges from flat.
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep import transport_sweep


def trace_single_sweep():
    """Run ONE sweep with constant source, trace cell 0 for each ordinate."""
    nx = 5
    R = 5.0
    mesh = Mesh1D(
        edges=np.linspace(0.0, R, nx + 1),
        mat_ids=np.ones(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC.reflective,
        bc_right=BC.reflective,
    )
    quad = GaussLegendre1D.create(4)  # S4
    sn_mesh = SNMesh(mesh, quad)

    N = quad.N
    ng = 1
    sig_t_val = 1.0
    Q_val = 1.0
    W = quad.weights.sum()

    sig_t = np.full((nx, 1, ng), sig_t_val)
    Q_iso = np.full((nx, 1, ng), Q_val)  # just the external source, no scattering

    # Initialize with flat angular flux at the boundary
    psi_bc = {}

    # Run ONE sweep
    angular_flux, scalar_flux = transport_sweep(Q_iso, sig_t, sn_mesh, psi_bc)

    expected_psi = Q_val / (W * sig_t_val)
    expected_phi = Q_val / sig_t_val

    print(f"\n{'='*70}")
    print(f"Single sweep trace (nx={nx}, S4, c=0, no BC history)")
    print(f"Expected psi = {expected_psi:.6f}, expected phi = {expected_phi:.6f}")
    print(f"\nPer-ordinate angular flux at each cell:")

    A = sn_mesh.face_areas
    V = sn_mesh.volumes[:, 0]
    mu = quad.mu_x
    weights = quad.weights

    for n in range(N):
        psi_cells = angular_flux[n, :, 0, 0]
        print(f"  n={n} (mu={mu[n]:+.6f}, w={weights[n]:.6f}): "
              f"psi = {psi_cells}")

    print(f"\nScalar flux per cell:")
    for i in range(nx):
        phi = scalar_flux[i, 0, 0]
        print(f"  cell {i}: phi = {phi:.8f}  (expected {expected_phi:.6f}, "
              f"err = {phi/expected_phi - 1:+.6f})")

    # Now run multiple sweeps to converge (no scattering, just repeated sweeps)
    print(f"\nIterative convergence (Q fixed, no scattering):")
    psi_bc2 = {}
    for it in range(20):
        angular_flux2, scalar_flux2 = transport_sweep(Q_iso, sig_t, sn_mesh, psi_bc2)
        phi_0 = scalar_flux2[0, 0, 0]
        phi_mid = scalar_flux2[nx//2, 0, 0]
        phi_last = scalar_flux2[-1, 0, 0]
        print(f"  iter {it+1:2d}: phi[0]={phi_0:.8f}  phi[mid]={phi_mid:.8f}  "
              f"phi[last]={phi_last:.8f}")

    # Check: after convergence, what is the boundary flux?
    if "bc_sph" in psi_bc2:
        bc = psi_bc2["bc_sph"]
        print(f"\nBoundary flux at r=R (per ordinate, stored for reflection):")
        for n in range(N):
            print(f"  n={n} (mu={mu[n]:+.6f}): bc = {bc[n, 0]:.8f}")

    print(f"{'='*70}")


def trace_with_source_iteration():
    """Run full source iteration (with scattering) and watch convergence."""
    nx = 5
    R = 5.0
    mesh = Mesh1D(
        edges=np.linspace(0.0, R, nx + 1),
        mat_ids=np.ones(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC.reflective,
        bc_right=BC.reflective,
    )
    quad = GaussLegendre1D.create(4)
    sn_mesh = SNMesh(mesh, quad)

    N = quad.N
    ng = 1
    sig_t_val = 1.0
    c = 0.5
    sig_s = c * sig_t_val
    Q_val = 1.0
    W = quad.weights.sum()

    sig_t = np.full((nx, 1, ng), sig_t_val)
    expected = Q_val / (sig_t_val * (1 - c))

    phi = np.ones((nx, 1, ng))
    psi_bc = {}

    print(f"\n{'='*70}")
    print(f"Source iteration with scattering (c={c})")
    print(f"Expected phi = {expected:.6f}")

    for it in range(100):
        phi_old = phi.copy()
        Q_iso = np.full((nx, 1, ng), Q_val) + sig_s * phi
        _, phi = transport_sweep(Q_iso, sig_t, sn_mesh, psi_bc)

        res = np.linalg.norm(phi - phi_old) / max(np.linalg.norm(phi), 1e-30)
        if it < 5 or it % 10 == 0 or res < 1e-12:
            print(f"  iter {it+1:3d}: phi[0]={phi[0,0,0]:.8f}  "
                  f"phi[mid]={phi[nx//2,0,0]:.8f}  "
                  f"phi[-1]={phi[-1,0,0]:.8f}  res={res:.2e}")
        if res < 1e-12:
            break

    print(f"{'='*70}")


if __name__ == "__main__":
    trace_single_sweep()
    trace_with_source_iteration()
