"""Diagnostic: Component isolation — zero out redistribution vs streaming.

Created by numerics-investigator on 2026-04-16.
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh


def manual_spherical_sweep(Q_iso_1d, sig_t_1d, sn_mesh, quad, psi_bc,
                           zero_redistribution=False, zero_streaming=False):
    """Manual sweep with option to zero out components."""
    nx = sn_mesh.nx
    N = quad.N
    ng = Q_iso_1d.shape[1]
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
    angular_flux = np.zeros((N, nx, ng))

    weight_norm = 1.0 / W
    QV_iso = Q_iso_1d * V[:, None] * weight_norm

    for n in range(N):
        mu_n = mu[n]
        abs_mu = abs(mu_n)
        w_n = weights[n]
        alpha_in = alpha[n]
        alpha_out = alpha[n + 1]
        tau_n = tau[n]

        if zero_redistribution:
            c_out = 0.0
            c_in = 0.0
        else:
            c_out = alpha_out / tau_n
            c_in = (1.0 - tau_n) / tau_n * alpha_out + alpha_in

        if mu_n < 0:
            if is_vacuum_outer:
                psi_spatial_in = np.zeros(ng)
            else:
                psi_spatial_in = bc_outer[ref[n]].copy()

            for i in range(nx - 1, -1, -1):
                A_in = A[i + 1]
                A_out = A[i]
                dA_w = dAw[i, n] if not zero_redistribution else 0.0

                if zero_streaming:
                    stream_denom = 0.0
                    stream_numer = 0.0
                else:
                    stream_denom = 2.0 * abs_mu * A_out
                    stream_numer = abs_mu * (A_in + A_out) * psi_spatial_in

                denom = stream_denom + dA_w * c_out + sig_t_1d[i] * V[i]
                numer = QV_iso[i] + stream_numer + dA_w * c_in * psi_angle[i]

                psi = numer / denom

                if not zero_streaming:
                    psi_spatial_out = 2.0 * psi - psi_spatial_in
                else:
                    psi_spatial_out = psi_spatial_in

                if not zero_redistribution:
                    psi_angle[i] = (psi - (1.0 - tau_n) * psi_angle[i]) / tau_n

                angular_flux[n, i] = psi
                scalar_flux[i] += w_n * psi

                psi_spatial_in = psi_spatial_out

        else:
            psi_spatial_in = np.zeros(ng)

            for i in range(nx):
                A_in = A[i]
                A_out = A[i + 1]
                dA_w = dAw[i, n] if not zero_redistribution else 0.0

                if zero_streaming:
                    stream_denom = 0.0
                    stream_numer = 0.0
                else:
                    stream_denom = 2.0 * abs_mu * A_out
                    stream_numer = abs_mu * (A_in + A_out) * psi_spatial_in

                denom = stream_denom + dA_w * c_out + sig_t_1d[i] * V[i]
                numer = QV_iso[i] + stream_numer + dA_w * c_in * psi_angle[i]

                psi = numer / denom

                if not zero_streaming:
                    psi_spatial_out = 2.0 * psi - psi_spatial_in
                else:
                    psi_spatial_out = psi_spatial_in

                if not zero_redistribution:
                    psi_angle[i] = (psi - (1.0 - tau_n) * psi_angle[i]) / tau_n

                angular_flux[n, i] = psi
                scalar_flux[i] += w_n * psi

                psi_spatial_in = psi_spatial_out

            bc_outer[n] = psi_spatial_out

    return angular_flux, scalar_flux


def run_test(label, zero_redistribution=False, zero_streaming=False,
             nx=10, N_ord=4, R=5.0, max_iter=100, tol=1e-14):
    """Run converging manual sweeps."""
    mesh = Mesh1D(
        edges=np.linspace(0.0, R, nx + 1),
        mat_ids=np.ones(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC.reflective,
        bc_right=BC.reflective,
    )
    quad = GaussLegendre1D.create(N_ord)
    sn_mesh = SNMesh(mesh, quad)

    ng = 1
    Q_1d = np.ones((nx, ng))
    sig_t_1d = np.ones((nx, ng))
    psi_bc = {}

    phi_old = np.zeros((nx, ng))
    for it in range(max_iter):
        _, phi = manual_spherical_sweep(
            Q_1d, sig_t_1d, sn_mesh, quad, psi_bc,
            zero_redistribution=zero_redistribution,
            zero_streaming=zero_streaming,
        )
        res = np.linalg.norm(phi - phi_old) / max(np.linalg.norm(phi), 1e-30)
        phi_old = phi.copy()
        if res < tol:
            break

    print(f"  {label:30s}: iters={it+1:3d}  "
          f"phi[0]={phi[0,0]:.6f}  phi[mid]={phi[nx//2,0]:.6f}  "
          f"phi[-1]={phi[-1,0]:.6f}  range=[{phi[:,0].min():.6f}, {phi[:,0].max():.6f}]")


def test_component_isolation():
    """Zero out components one at a time."""
    print(f"\n{'='*70}")
    print(f"Component isolation: expected phi = 1.0 everywhere")
    run_test("Full sweep (baseline)")
    run_test("No redistribution (alpha=0)", zero_redistribution=True)
    # Note: zeroing streaming doesn't make physical sense for this test
    # as it removes the spatial coupling entirely
    print(f"{'='*70}")


if __name__ == "__main__":
    test_component_isolation()
