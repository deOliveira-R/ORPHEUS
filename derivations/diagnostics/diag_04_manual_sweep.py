"""Diagnostic: Manual step-by-step sweep at cell 0 to find the bug.

Created by numerics-investigator on 2026-04-16.
Reproduce the sweep calculation by hand for cell 0 and compare with
the actual sweep output.
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep import transport_sweep


def manual_sweep():
    """Reproduce the sweep manually for a 3-cell spherical mesh."""
    nx = 3
    R = 3.0
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
    Q_val = 1.0
    W = quad.weights.sum()
    weight_norm = 1.0 / W

    A = sn_mesh.face_areas
    V = sn_mesh.volumes[:, 0]
    alpha = sn_mesh.alpha_half
    dAw = sn_mesh.redist_dAw
    tau = sn_mesh.tau_mm
    mu = quad.mu_x
    weights = quad.weights
    ref = quad.reflection_index("x")

    sig_t = np.full((nx, 1, ng), sig_t_val)
    Q_iso = np.full((nx, 1, ng), Q_val)

    print(f"\n{'='*70}")
    print(f"Manual sweep: nx={nx}, R={R}, S4")
    print(f"  W = {W:.6f}")
    print(f"  A = {A}")
    print(f"  V = {V}")
    print(f"  alpha = {alpha}")
    print(f"  tau = {tau}")
    print(f"  mu = {mu}")
    print(f"  weights = {weights}")
    print(f"  ref = {ref}")
    print(f"  dAw[0,:] = {dAw[0,:]}")
    print()

    # Now reproduce the sweep manually
    QV_iso = Q_val * V * weight_norm  # (nx,)
    psi_angle = np.zeros(nx)
    bc_outer = np.zeros((N,))
    is_vacuum_outer = False

    scalar_flux_manual = np.zeros(nx)
    angular_flux_manual = np.zeros((N, nx))

    for n in range(N):
        mu_n = mu[n]
        abs_mu = abs(mu_n)
        w_n = weights[n]
        alpha_in = alpha[n]
        alpha_out = alpha[n + 1]
        tau_n = tau[n]
        c_out = alpha_out / tau_n
        c_in = (1.0 - tau_n) / tau_n * alpha_out + alpha_in

        print(f"Ordinate n={n}: mu={mu_n:+.6f}, w={w_n:.6f}")
        print(f"  alpha_in={alpha_in:.6f}, alpha_out={alpha_out:.6f}, tau={tau_n:.6f}")
        print(f"  c_out={c_out:.6f}, c_in={c_in:.6f}")

        if mu_n < 0:
            # Inward sweep
            if is_vacuum_outer:
                psi_spatial_in = 0.0
            else:
                psi_spatial_in = bc_outer[ref[n]]

            print(f"  INWARD: psi_spatial_in (from BC) = {psi_spatial_in:.8f}")

            for i in range(nx - 1, -1, -1):
                A_in = A[i + 1]
                A_out = A[i]
                dA_w = dAw[i, n]

                denom = 2.0 * abs_mu * A_out + dA_w * c_out + sig_t_val * V[i]
                numer = (QV_iso[i]
                         + abs_mu * (A_in + A_out) * psi_spatial_in
                         + dA_w * c_in * psi_angle[i])

                psi = numer / denom
                psi_spatial_out = 2.0 * psi - psi_spatial_in
                psi_angle[i] = (psi - (1.0 - tau_n) * psi_angle[i]) / tau_n

                angular_flux_manual[n, i] = psi
                scalar_flux_manual[i] += w_n * psi

                print(f"    cell {i}: A_in={A_in:.4f}, A_out={A_out:.4f}, "
                      f"dAw={dA_w:.4f}")
                print(f"      denom = {denom:.6f}")
                print(f"        2|mu|A_out={2*abs_mu*A_out:.6f}, "
                      f"dAw*c_out={dA_w*c_out:.6f}, "
                      f"sig_t*V={sig_t_val*V[i]:.6f}")
                print(f"      numer = {numer:.6f}")
                print(f"        QV={QV_iso[i]:.6f}, "
                      f"|mu|(A_in+A_out)*psi_in={abs_mu*(A_in+A_out)*psi_spatial_in:.6f}, "
                      f"dAw*c_in*psi_angle={dA_w*c_in*psi_angle[i]:.6f}")
                print(f"      psi = {psi:.8f}")
                print(f"      psi_spatial_out = {psi_spatial_out:.8f}")
                print(f"      psi_angle[{i}] (updated) = {psi_angle[i]:.8f}")

                psi_spatial_in = psi_spatial_out

        else:
            # Outward sweep
            psi_spatial_in = 0.0
            print(f"  OUTWARD: psi_spatial_in = 0 (r=0)")

            for i in range(nx):
                A_in = A[i]
                A_out = A[i + 1]
                dA_w = dAw[i, n]

                denom = 2.0 * abs_mu * A_out + dA_w * c_out + sig_t_val * V[i]
                numer = (QV_iso[i]
                         + abs_mu * (A_in + A_out) * psi_spatial_in
                         + dA_w * c_in * psi_angle[i])

                psi = numer / denom
                psi_spatial_out = 2.0 * psi - psi_spatial_in
                psi_angle[i] = (psi - (1.0 - tau_n) * psi_angle[i]) / tau_n

                angular_flux_manual[n, i] = psi
                scalar_flux_manual[i] += w_n * psi

                print(f"    cell {i}: A_in={A_in:.4f}, A_out={A_out:.4f}, "
                      f"dAw={dA_w:.4f}")
                print(f"      denom = {denom:.6f}")
                print(f"      numer = {numer:.6f}")
                print(f"      psi = {psi:.8f}")
                print(f"      psi_spatial_out = {psi_spatial_out:.8f}")
                print(f"      psi_angle[{i}] (updated) = {psi_angle[i]:.8f}")

                psi_spatial_in = psi_spatial_out

            bc_outer[n] = psi_spatial_out
            print(f"  bc_outer[{n}] = {psi_spatial_out:.8f}")

        print()

    # Compare with actual sweep
    psi_bc_actual = {}
    angular_actual, scalar_actual = transport_sweep(Q_iso, sig_t, sn_mesh, psi_bc_actual)

    print(f"\n{'='*70}")
    print(f"Comparison: manual vs actual sweep")
    print(f"  {'Cell':>4s}  {'Manual phi':>12s}  {'Actual phi':>12s}  {'Diff':>12s}")
    for i in range(nx):
        m = scalar_flux_manual[i]
        a = scalar_actual[i, 0, 0]
        print(f"  {i:4d}  {m:12.8f}  {a:12.8f}  {m-a:12.2e}")

    print(f"\nExpected phi = {Q_val / sig_t_val:.6f} everywhere")
    print(f"{'='*70}")


if __name__ == "__main__":
    manual_sweep()
