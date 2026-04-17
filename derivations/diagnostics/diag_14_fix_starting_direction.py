"""Diagnostic: Fix the starting-direction angular face flux initialization.

Created by numerics-investigator on 2026-04-16.
Hypothesis: the WDD closure propagates psi_angle_in = 0 through ordinate
sequence, corrupting all subsequent ordinates. Fix: after computing the
first ordinate, set psi_angle = psi (cell average) instead of WDD output.
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh


def sweep_fixed(sn_mesh, quad, Q_1d, sig_t_1d, psi_bc, fix_start=True):
    """Spherical sweep with optional starting-direction fix."""
    nx = sn_mesh.nx
    N = quad.N
    ng = Q_1d.shape[1]
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
    angular_flux = np.zeros((N, nx, 1, ng))

    weight_norm = 1.0 / W
    QV_iso = Q_1d * V[:, None] * weight_norm

    for n in range(N):
        mu_n = mu[n]
        abs_mu = abs(mu_n)
        w_n = weights[n]
        alpha_in = alpha[n]
        alpha_out = alpha[n + 1]
        tau_n = tau[n]
        c_out = alpha_out / tau_n
        c_in = (1.0 - tau_n) / tau_n * alpha_out + alpha_in

        # Starting direction fix: for the first ordinate, the balance
        # equation is independent of psi_angle_in (alpha_in=0 nullifies it).
        # But the WDD closure psi_angle_out = (psi - (1-tau)*psi_angle_in)/tau
        # IS sensitive to psi_angle_in. When psi_angle_in=0 and tau=0.5,
        # psi_angle_out = 2*psi, which corrupts subsequent ordinates.
        #
        # Fix: use psi_angle_in = psi (cell average) for the starting
        # direction, so WDD gives psi_angle_out = psi (preserving flat flux).
        is_start = (n == 0)

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

                if fix_start and is_start:
                    # For starting direction: set psi_angle = psi (cell avg)
                    # This makes WDD output = psi for flat flux
                    psi_angle[i] = psi
                else:
                    psi_angle[i] = (psi - (1.0 - tau_n) * psi_angle[i]) / tau_n

                angular_flux[n, i, 0] = psi
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

                if fix_start and is_start:
                    psi_angle[i] = psi
                else:
                    psi_angle[i] = (psi - (1.0 - tau_n) * psi_angle[i]) / tau_n

                angular_flux[n, i, 0] = psi
                scalar_flux[i] += w_n * psi
                psi_spatial_in = psi_spatial_out

            bc_outer[n] = psi_spatial_out

    return angular_flux, scalar_flux[:, None, :]


def run_test(label, fix_start, nx=10, N_ord=8, R=10.0, max_iter=200, tol=1e-14):
    """Run converging sweeps."""
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

    phi_old = np.zeros((nx, 1, ng))
    for it in range(max_iter):
        _, phi = sweep_fixed(sn_mesh, quad, Q_1d, sig_t_1d, psi_bc, fix_start=fix_start)
        res = np.linalg.norm(phi - phi_old) / max(np.linalg.norm(phi), 1e-30)
        phi_old = phi.copy()
        if res < tol:
            break

    phi_1d = phi[:, 0, 0]
    expected = 1.0
    max_err = np.max(np.abs(phi_1d / expected - 1))

    print(f"  {label:40s}: iters={it+1:3d}  "
          f"phi[0]={phi_1d[0]:.8f}  phi[mid]={phi_1d[nx//2]:.8f}  "
          f"phi[-1]={phi_1d[-1]:.8f}  max|err|={max_err:.6e}")

    return phi_1d


def test_starting_direction_fix():
    """Compare original vs fixed sweep."""
    print(f"\n{'='*70}")
    print(f"Starting-direction fix comparison (expected phi = 1.0)")

    run_test("Original (psi_angle_init = 0)", fix_start=False)
    run_test("Fixed (psi_angle = psi at n=0)", fix_start=True)

    # Test scaling with fix
    print(f"\nScaling analysis with fix:")
    print(f"  {'Cells':>6s}  {'Max|err|':>12s}  {'Ratio':>8s}")
    prev_err = None
    for nx in [5, 10, 20, 40]:
        phi = run_test(f"  fix, nx={nx}", fix_start=True, nx=nx)
        err = np.max(np.abs(phi / 1.0 - 1))
        ratio = prev_err / err if prev_err is not None and err > 0 else float('nan')
        print(f"  {nx:6d}  {err:12.6e}  {ratio:8.2f}")
        prev_err = err

    print(f"{'='*70}")


if __name__ == "__main__":
    test_starting_direction_fix()
