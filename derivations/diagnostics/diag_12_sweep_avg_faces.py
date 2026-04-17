"""Diagnostic: Sweep using arithmetic average face fluxes instead of DD/WDD.

Created by numerics-investigator on 2026-04-16.
Tests whether the DD/WDD closure is the source of the error.
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh


def sweep_jacobi(sn_mesh, quad, Q_iso_1d, sig_t_1d, psi_old):
    """One Jacobi-like sweep: compute new psi from old psi using the
    BiCGSTAB-style face flux approximation (arithmetic averages).

    This is NOT a traditional sweep — it's a simultaneous update.
    """
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

    weight_norm = 1.0 / W
    S = Q_iso_1d * weight_norm  # source density per ordinate (nx, ng)

    psi_new = np.zeros((N, nx, ng))

    for n in range(N):
        mu_n = mu[n]
        tau_n = tau[n]

        for i in range(nx):
            # Spatial face fluxes: arithmetic average of cell centers
            # Right face (i+1/2)
            if i < nx - 1:
                psi_right = 0.5 * (psi_old[n, i] + psi_old[n, i + 1])
            else:
                # Outer boundary: outgoing = cell value, incoming = reflected
                if mu_n > 1e-15:
                    psi_right = psi_old[n, i]
                else:
                    psi_right = psi_old[ref[n], i]

            # Left face (i-1/2)
            if i > 0:
                psi_left = 0.5 * (psi_old[n, i - 1] + psi_old[n, i])
            else:
                # r = 0: A[0] = 0, so this doesn't matter
                psi_left = 0.0

            streaming = mu_n * (A[i + 1] * psi_right - A[i] * psi_left) / V[i]

            # Angular face fluxes: M-M weighted interpolation
            if n < N - 1:
                psi_angle_right = tau_n * psi_old[n + 1, i] + (1.0 - tau_n) * psi_old[n, i]
            else:
                psi_angle_right = psi_old[n, i]

            if n > 0:
                psi_angle_left = tau[n - 1] * psi_old[n, i] + (1.0 - tau[n - 1]) * psi_old[n - 1, i]
            else:
                psi_angle_left = psi_old[n, i]

            redistribution = dAw[i, n] * (alpha[n + 1] * psi_angle_right
                                          - alpha[n] * psi_angle_left) / V[i]

            collision = sig_t_1d[i] * psi_old[n, i]

            # T*psi = S => psi_new = (S - streaming - redistribution) / Sig_t
            # Actually for Jacobi: psi_new = (S*V - (streaming+redistribution)*V) / (Sig_t*V)
            # But we can't just do this because the operator couples psi[n,i] to neighbors
            # Let me think...
            # The operator is: T*psi = S, where T = streaming + redistribution + collision
            # For a point relaxation: psi_new[n,i] = (S - off_diag) / diag
            # The diagonal of T is Sig_t (collision part) + self-coupling from streaming and redistribution

            # Actually, let's just do Richardson iteration: psi_new = psi_old + omega*(S - T*psi_old)/Sig_t
            rhs = S[i]
            residual = rhs - (streaming + redistribution + collision)
            psi_new[n, i] = psi_old[n, i] + residual / sig_t_1d[i]

    return psi_new


def test_jacobi_iteration():
    """Run Jacobi-style iteration with arithmetic face fluxes."""
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
    Q_1d = np.ones((nx, ng))
    sig_t_1d = np.ones((nx, ng))

    # Start from flat guess
    W = quad.weights.sum()
    psi = np.full((N, nx, ng), 1.0 / W)

    print(f"\n{'='*70}")
    print(f"Jacobi iteration with arithmetic face fluxes (BiCGSTAB-style)")

    for it in range(500):
        psi_new = sweep_jacobi(sn_mesh, quad, Q_1d, sig_t_1d, psi)
        res = np.linalg.norm(psi_new - psi) / max(np.linalg.norm(psi), 1e-30)

        # Compute scalar flux
        phi = np.zeros(nx)
        for n in range(N):
            phi += quad.weights[n] * psi_new[n, :, 0]

        if it < 5 or it % 50 == 0 or res < 1e-12:
            print(f"  iter {it+1:3d}: phi[0]={phi[0]:.8f}  phi[mid]={phi[nx//2]:.8f}  "
                  f"phi[-1]={phi[-1]:.8f}  res={res:.2e}")

        if res < 1e-12:
            break

        psi = psi_new

    print(f"\nExpected phi = 1.0 everywhere")
    print(f"{'='*70}")


if __name__ == "__main__":
    test_jacobi_iteration()
