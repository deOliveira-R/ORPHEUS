"""Diagnostic: Start the sweep from flat flux and check if it preserves it.

Created by numerics-investigator on 2026-04-16.
If starting from flat preserves flat, the sweep's issue is in its
ITERATIVE convergence path, not in the balance equation itself.
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep import transport_sweep, _sweep_1d_spherical


def test_flat_start():
    """Initialize boundary fluxes to the flat-flux value and sweep once."""
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
    psi_flat = Q_val / (W * sig_t_val)

    sig_t = np.full((nx, 1, ng), sig_t_val)
    Q_iso = np.full((nx, 1, ng), Q_val)

    # Initialize boundary flux to the flat value
    psi_bc = {
        "bc_sph": np.full((N, ng), psi_flat),
    }

    # Run ONE sweep
    angular, scalar = transport_sweep(Q_iso, sig_t, sn_mesh, psi_bc)
    phi = scalar[:, 0, 0]

    expected = Q_val / sig_t_val

    print(f"\n{'='*70}")
    print(f"Sweep from flat start (psi_flat={psi_flat:.6f})")
    print(f"Expected phi = {expected:.6f}")
    print(f"\nAfter 1 sweep:")
    for i in range(nx):
        print(f"  cell {i:2d}: phi={phi[i]:.10f}  err={phi[i]/expected - 1:+.6e}")

    print(f"\n  Max |err|: {np.max(np.abs(phi/expected - 1)):.6e}")

    # Check: does a second sweep preserve it?
    angular2, scalar2 = transport_sweep(Q_iso, sig_t, sn_mesh, psi_bc)
    phi2 = scalar2[:, 0, 0]
    print(f"\nAfter 2 sweeps:")
    for i in range(nx):
        print(f"  cell {i:2d}: phi={phi2[i]:.10f}  err={phi2[i]/expected - 1:+.6e}")

    print(f"\n  Change from sweep 1 to 2: {np.max(np.abs(phi2 - phi)):.6e}")
    print(f"{'='*70}")


if __name__ == "__main__":
    test_flat_start()
