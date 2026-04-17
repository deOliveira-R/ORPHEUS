"""Diagnostic: Check global conservation — integrate the balance over all
ordinates and all cells. For a constant-source, reflective-BC problem,
the total source should equal the total absorption.

Created by numerics-investigator on 2026-04-16.
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep import transport_sweep


def test_conservation():
    """Check volume-weighted flux integral = total source / Sig_t."""
    nx = 20
    N_ord = 8
    R = 10.0

    mesh = Mesh1D(
        edges=np.linspace(0.0, R, nx + 1),
        mat_ids=np.ones(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC.reflective,
        bc_right=BC.reflective,
    )
    quad = GaussLegendre1D.create(N_ord)
    sn_mesh = SNMesh(mesh, quad)

    sig_t_val = 1.0
    Q_val = 1.0
    ng = 1
    sig_t = np.full((nx, 1, ng), sig_t_val)
    Q_iso = np.full((nx, 1, ng), Q_val)
    psi_bc = {}

    # Converge
    for it in range(200):
        _, phi = transport_sweep(Q_iso, sig_t, sn_mesh, psi_bc)
        if it > 0:
            res = np.linalg.norm(phi - phi_old) / max(np.linalg.norm(phi), 1e-30)
            if res < 1e-14:
                break
        phi_old = phi.copy()

    V = sn_mesh.volumes[:, 0]
    phi_1d = phi[:, 0, 0]

    # Conservation check:
    # Total source = sum(Q * V) = Q * V_total
    # Total absorption = sum(Sig_t * phi * V)
    # These should be equal (reflective BCs = no leakage)
    total_source = Q_val * V.sum()
    total_absorption = sig_t_val * np.sum(phi_1d * V)
    total_volume = V.sum()

    # Volume-weighted average flux
    phi_vol_avg = np.sum(phi_1d * V) / V.sum()
    expected_phi = Q_val / sig_t_val

    print(f"\n{'='*70}")
    print(f"Conservation check (nx={nx}, S{N_ord})")
    print(f"  Total volume:     {total_volume:.6f}")
    print(f"  Total source:     {total_source:.6f}")
    print(f"  Total absorption: {total_absorption:.6f}")
    print(f"  Source/Absorption: {total_source/total_absorption:.10f}")
    print(f"  Volume-avg phi:   {phi_vol_avg:.10f} (expected {expected_phi})")
    print(f"  Relative error in vol-avg: {phi_vol_avg/expected_phi - 1:.6e}")
    print(f"\n  phi profile (unweighted):")
    for i in range(nx):
        r_mid = 0.5 * (mesh.edges[i] + mesh.edges[i+1])
        print(f"    cell {i:2d} (r={r_mid:5.2f}): phi={phi_1d[i]:.8f}  "
              f"V={V[i]:.4f}  phi*V={phi_1d[i]*V[i]:.4f}")
    print(f"{'='*70}")


if __name__ == "__main__":
    test_conservation()
