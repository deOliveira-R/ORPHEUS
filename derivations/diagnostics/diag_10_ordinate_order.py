"""Diagnostic: Test if reversing ordinate order changes the solution.

Created by numerics-investigator on 2026-04-16.
If the solution depends on ordinate order, the angular redistribution
coupling is the source of the error.
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh


def sweep_with_order(sn_mesh, quad, Q_iso_1d, sig_t_1d, psi_bc, reverse=False):
    """Sweep with configurable ordinate processing order."""
    nx = sn_mesh.nx
    N = quad.N
    ng = 1
    mu = quad.mu_x
    weights = quad.weights
    W = weights.sum()
    ref = quad.reflection_index("x")

    A = sn_mesh.face_areas
    V = sn_mesh.volumes[:, 0]
    alpha = sn_mesh.alpha_half
    dAw = sn_mesh.redist_dAw
    tau = sn_mesh.tau_mm

    is_vacuum_outer = sn_mesh.bc_right == "vacuum"

    if "bc_sph" not in psi_bc:
        psi_bc["bc_sph"] = np.zeros((N, ng))
    bc_outer = psi_bc["bc_sph"]

    psi_angle = np.zeros((nx, ng))
    scalar_flux = np.zeros((nx, ng))

    weight_norm = 1.0 / W
    QV_iso = Q_iso_1d * V[:, None] * weight_norm

    ordinate_order = range(N - 1, -1, -1) if reverse else range(N)

    for n in ordinate_order:
        mu_n = mu[n]
        abs_mu = abs(mu_n)
        w_n = weights[n]

        if reverse:
            # Reversed: process from n=N-1 to n=0
            # α_in is α_{n+1/2} (the "higher" edge)
            # α_out is α_{n-1/2} (the "lower" edge)
            alpha_in = alpha[n + 1]
            alpha_out = alpha[n]
        else:
            # Normal order
            alpha_in = alpha[n]
            alpha_out = alpha[n + 1]

        tau_n = tau[n]
        c_out = alpha_out / tau_n if tau_n > 0 else 0.0
        c_in = (1.0 - tau_n) / tau_n * alpha_out + alpha_in if tau_n > 0 else alpha_in

        if mu_n < 0:
            psi_spatial_in = bc_outer[ref[n]].copy() if not is_vacuum_outer else np.zeros(ng)
            for i in range(nx - 1, -1, -1):
                A_in = A[i + 1]
                A_out = A[i]
                dA_w = dAw[i, n]

                denom = 2.0 * abs_mu * A_out + dA_w * c_out + sig_t_1d[i] * V[i]
                numer = QV_iso[i] + abs_mu * (A_in + A_out) * psi_spatial_in + dA_w * c_in * psi_angle[i]

                psi = numer / denom
                psi_spatial_out = 2.0 * psi - psi_spatial_in
                psi_angle[i] = (psi - (1.0 - tau_n) * psi_angle[i]) / tau_n

                scalar_flux[i] += w_n * psi
                psi_spatial_in = psi_spatial_out
        else:
            psi_spatial_in = np.zeros(ng)
            for i in range(nx):
                A_in = A[i]
                A_out = A[i + 1]
                dA_w = dAw[i, n]

                denom = 2.0 * abs_mu * A_out + dA_w * c_out + sig_t_1d[i] * V[i]
                numer = QV_iso[i] + abs_mu * (A_in + A_out) * psi_spatial_in + dA_w * c_in * psi_angle[i]

                psi = numer / denom
                psi_spatial_out = 2.0 * psi - psi_spatial_in
                psi_angle[i] = (psi - (1.0 - tau_n) * psi_angle[i]) / tau_n

                scalar_flux[i] += w_n * psi
                psi_spatial_in = psi_spatial_out

            bc_outer[n] = psi_spatial_out

    return scalar_flux


def test_ordinate_order():
    """Compare standard vs reversed ordinate order."""
    nx = 5
    R = 5.0
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

    Q_1d = np.ones((nx, 1))
    sig_t_1d = np.ones((nx, 1))

    # Standard order
    psi_bc1 = {}
    for it in range(100):
        phi1 = sweep_with_order(sn_mesh, quad, Q_1d, sig_t_1d, psi_bc1, reverse=False)
        if it > 0:
            res = np.linalg.norm(phi1 - phi1_old) / max(np.linalg.norm(phi1), 1e-30)
            if res < 1e-14:
                break
        phi1_old = phi1.copy()

    # Reversed order
    psi_bc2 = {}
    for it in range(100):
        phi2 = sweep_with_order(sn_mesh, quad, Q_1d, sig_t_1d, psi_bc2, reverse=True)
        if it > 0:
            res = np.linalg.norm(phi2 - phi2_old) / max(np.linalg.norm(phi2), 1e-30)
            if res < 1e-14:
                break
        phi2_old = phi2.copy()

    print(f"\n{'='*70}")
    print(f"Ordinate order comparison (nx={nx}, S{N_ord})")
    print(f"Expected phi = 1.0 everywhere")
    print(f"  {'Cell':>4s}  {'Standard':>12s}  {'Reversed':>12s}  {'Diff':>12s}")
    for i in range(nx):
        print(f"  {i:4d}  {phi1[i,0]:12.8f}  {phi2[i,0]:12.8f}  {phi1[i,0]-phi2[i,0]:12.6e}")
    print(f"{'='*70}")


if __name__ == "__main__":
    test_ordinate_order()
