"""Diagnostic: Check for negative face fluxes and trace psi_angle evolution.

Created by numerics-investigator on 2026-04-16.
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh


def instrumented_sweep(Q_iso_1d, sig_t_1d, sn_mesh, quad, psi_bc):
    """Sweep with full diagnostics on face fluxes."""
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

    neg_spatial = 0
    neg_angular = 0

    for n in range(N):
        mu_n = mu[n]
        abs_mu = abs(mu_n)
        w_n = weights[n]
        alpha_in = alpha[n]
        alpha_out = alpha[n + 1]
        tau_n = tau[n]
        c_out = alpha_out / tau_n
        c_in = (1.0 - tau_n) / tau_n * alpha_out + alpha_in

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
                psi_angle_out = (psi - (1.0 - tau_n) * psi_angle[i]) / tau_n

                if np.any(psi_spatial_out < -1e-15):
                    neg_spatial += 1
                if np.any(psi_angle_out < -1e-15):
                    neg_angular += 1

                psi_angle[i] = psi_angle_out
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
                psi_angle_out = (psi - (1.0 - tau_n) * psi_angle[i]) / tau_n

                if np.any(psi_spatial_out < -1e-15):
                    neg_spatial += 1
                if np.any(psi_angle_out < -1e-15):
                    neg_angular += 1

                psi_angle[i] = psi_angle_out
                scalar_flux[i] += w_n * psi
                psi_spatial_in = psi_spatial_out

            bc_outer[n] = psi_spatial_out

    return scalar_flux, neg_spatial, neg_angular


def test_negatives():
    """Check for negative face fluxes in converged solution."""
    nx = 10
    R = 10.0
    N_ord = 8

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
    psi_bc = {}

    print(f"\n{'='*70}")
    print(f"Negative face flux check (nx={nx}, S{N_ord})")

    for it in range(50):
        phi, neg_s, neg_a = instrumented_sweep(Q_1d, sig_t_1d, sn_mesh, quad, psi_bc)
        if neg_s > 0 or neg_a > 0:
            print(f"  iter {it+1}: neg_spatial={neg_s}, neg_angular={neg_a}, phi[0]={phi[0,0]:.6f}")
        if it > 0:
            res = np.linalg.norm(phi - phi_old) / max(np.linalg.norm(phi), 1e-30)
            if res < 1e-14:
                print(f"  Converged at iter {it+1}")
                break
        phi_old = phi.copy()

    print(f"  Final phi[0] = {phi[0,0]:.8f}, phi[-1] = {phi[-1,0]:.8f}")
    print(f"{'='*70}")


def test_per_ordinate_consistency_converged():
    """For the converged solution, check per-ordinate streaming + redistribution.

    For flat flux, streaming_n + redist_n = 0 for all n.
    For the (wrong) converged solution, what is the per-ordinate residual?
    """
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

    N = quad.N
    mu = quad.mu_x
    weights = quad.weights
    W = weights.sum()
    A = sn_mesh.face_areas
    V = sn_mesh.volumes[:, 0]
    alpha = sn_mesh.alpha_half
    dAw = sn_mesh.redist_dAw
    tau = sn_mesh.tau_mm
    ref = quad.reflection_index("x")

    sig_t_val = 1.0
    Q_val = 1.0
    S = Q_val / W

    # Converge with the actual sweep
    from orpheus.sn.sweep import transport_sweep
    sig_t = np.full((nx, 1, 1), sig_t_val)
    Q_iso = np.full((nx, 1, 1), Q_val)
    psi_bc = {}
    for it in range(100):
        angular, scalar = transport_sweep(Q_iso, sig_t, sn_mesh, psi_bc)
        if it > 0:
            res = np.linalg.norm(scalar - scalar_old) / max(np.linalg.norm(scalar), 1e-30)
            if res < 1e-14:
                break
        scalar_old = scalar.copy()

    psi = angular[:, :, 0, 0]  # (N, nx)
    bc_outer = psi_bc["bc_sph"][:, 0]

    print(f"\n{'='*70}")
    print(f"Per-ordinate balance check at cell 0")
    print(f"  Converged phi[0] = {scalar[0,0,0]:.8f} (expected {Q_val/sig_t_val})")

    # Reconstruct face fluxes by re-running the sweep logic
    psi_angle = np.zeros(nx)

    for n in range(N):
        mu_n = mu[n]
        abs_mu = abs(mu_n)
        alpha_in = alpha[n]
        alpha_out = alpha[n+1]
        tau_n = tau[n]

        if mu_n < 0:
            psi_spatial_in = bc_outer[ref[n]]
            for i in range(nx-1, -1, -1):
                psi_n_i = psi[n, i]
                A_outer = A[i+1]
                A_inner = A[i]

                psi_spatial_out = 2.0 * psi_n_i - psi_spatial_in
                psi_angle_old = psi_angle[i]
                psi_angle_out = (psi_n_i - (1.0 - tau_n) * psi_angle_old) / tau_n

                if i == 0:
                    streaming = mu_n * (A_outer * psi_spatial_in - A_inner * psi_spatial_out)
                    redist = dAw[i, n] * (alpha_out * psi_angle_out - alpha_in * psi_angle_old)
                    collision = sig_t_val * V[i] * psi_n_i
                    source = S * V[i]
                    balance = streaming + redist + collision - source
                    print(f"  n={n} (mu={mu_n:+.6f}): "
                          f"psi={psi_n_i:.6f}, psi_in={psi_spatial_in:.6f}, "
                          f"psi_out={psi_spatial_out:.6f}")
                    print(f"    stream={streaming:.6e}, redist={redist:.6e}, "
                          f"coll={collision:.6e}, src={source:.6e}, "
                          f"balance={balance:.6e}")
                    print(f"    psi_angle_in={psi_angle_old:.6f}, "
                          f"psi_angle_out={psi_angle_out:.6f}")

                psi_angle[i] = psi_angle_out
                psi_spatial_in = psi_spatial_out

        else:
            psi_spatial_in = 0.0
            for i in range(nx):
                psi_n_i = psi[n, i]
                A_inner = A[i]
                A_outer = A[i+1]

                psi_spatial_out = 2.0 * psi_n_i - psi_spatial_in
                psi_angle_old = psi_angle[i]
                psi_angle_out = (psi_n_i - (1.0 - tau_n) * psi_angle_old) / tau_n

                if i == 0:
                    streaming = mu_n * (A_outer * psi_spatial_out - A_inner * psi_spatial_in)
                    redist = dAw[i, n] * (alpha_out * psi_angle_out - alpha_in * psi_angle_old)
                    collision = sig_t_val * V[i] * psi_n_i
                    source = S * V[i]
                    balance = streaming + redist + collision - source
                    print(f"  n={n} (mu={mu_n:+.6f}): "
                          f"psi={psi_n_i:.6f}, psi_in={psi_spatial_in:.6f}, "
                          f"psi_out={psi_spatial_out:.6f}")
                    print(f"    stream={streaming:.6e}, redist={redist:.6e}, "
                          f"coll={collision:.6e}, src={source:.6e}, "
                          f"balance={balance:.6e}")
                    print(f"    psi_angle_in={psi_angle_old:.6f}, "
                          f"psi_angle_out={psi_angle_out:.6f}")

                psi_angle[i] = psi_angle_out
                psi_spatial_in = psi_spatial_out

    print(f"{'='*70}")


if __name__ == "__main__":
    test_negatives()
    test_per_ordinate_consistency_converged()
