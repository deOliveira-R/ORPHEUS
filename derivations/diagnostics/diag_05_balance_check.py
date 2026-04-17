"""Diagnostic: Check if the converged sweep solution satisfies the balance equation.

Created by numerics-investigator on 2026-04-16.
If the converged solution satisfies the balance equation, the issue is
that the balance equation itself doesn't have a unique flat solution.
If it doesn't satisfy the balance equation, the sweep has a bug.
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep import transport_sweep


def converge_sweep(nx=5, N_ord=4, R=5.0, max_iter=200, tol=1e-14):
    """Converge the sweep with constant source, no scattering."""
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

    for it in range(max_iter):
        ang_old = None
        angular, scalar = transport_sweep(Q_iso, sig_t, sn_mesh, psi_bc)
        if it > 0:
            res = np.linalg.norm(scalar - scalar_old) / max(np.linalg.norm(scalar), 1e-30)
            if res < tol:
                break
        scalar_old = scalar.copy()

    return angular, scalar, sn_mesh, psi_bc, quad


def test_balance_residual():
    """Check if converged solution satisfies the discrete balance equation."""
    nx = 5
    N_ord = 4
    R = 5.0

    angular, scalar, sn_mesh, psi_bc, quad = converge_sweep(nx, N_ord, R)

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
    S = Q_val / W  # isotropic source density per ordinate

    # Extract angular flux (N, nx)
    psi = angular[:, :, 0, 0]

    # Get boundary flux
    bc_outer = psi_bc["bc_sph"][:, 0]

    print(f"\n{'='*70}")
    print(f"Balance residual check for converged solution")
    print(f"Expected phi = {Q_val/sig_t_val:.6f}")
    print(f"Converged phi: {scalar[:, 0, 0]}")
    print(f"\nBoundary flux at r=R:")
    for n in range(N):
        print(f"  n={n} (mu={mu[n]:+.6f}): bc={bc_outer[n]:.8f}")

    # For each ordinate and cell, check the balance equation:
    # mu_n * [A_{i+1/2}*psi_{i+1/2} - A_{i-1/2}*psi_{i-1/2}]
    # + (dAw) * [alpha_{n+1/2}*psi_{n+1/2} - alpha_{n-1/2}*psi_{n-1/2}]
    # + Sig_t * V * psi = S * V
    #
    # We need to reconstruct the face fluxes from the cell-center values.
    # The sweep uses DD: psi = (psi_in + psi_out)/2
    # So psi_out = 2*psi - psi_in (spatial)
    # And WDD: psi_angle_out = (psi - (1-tau)*psi_angle_in) / tau

    # Reconstruct by re-running the sweep logic to get face fluxes
    psi_angle = np.zeros(nx)

    print(f"\nPer-ordinate balance residuals:")
    max_residual = 0.0

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

                # DD face flux
                psi_spatial_out = 2.0 * psi_n_i - psi_spatial_in

                # WDD angular face flux
                psi_angle_in = psi_angle[i]  # alpha_{n-1/2} side
                psi_angle_out = (psi_n_i - (1.0 - tau_n) * psi_angle_in) / tau_n

                # Balance: mu_n * [A_out*psi_out - A_in*psi_in]  (for inward)
                # Wait: the equation uses A_{i+1/2} and A_{i-1/2} with mu_n (signed)
                # For standard indexing: i+1/2 is outer face, i-1/2 is inner face
                A_outer = A[i+1]  # A_{i+1/2}
                A_inner = A[i]    # A_{i-1/2}

                # For inward: psi at A_outer = psi_spatial_in, psi at A_inner = psi_spatial_out
                streaming = mu_n * (A_outer * psi_spatial_in - A_inner * psi_spatial_out)

                # Redistribution
                redist = dAw[i, n] * (alpha_out * psi_angle_out - alpha_in * psi_angle_in)

                collision = sig_t_val * V[i] * psi_n_i
                source = S * V[i]

                residual = streaming + redist + collision - source

                if abs(residual) > 1e-8:
                    print(f"  n={n} cell={i}: stream={streaming:.6e} "
                          f"redist={redist:.6e} coll={collision:.6e} "
                          f"src={source:.6e} RESIDUAL={residual:.6e}")

                max_residual = max(max_residual, abs(residual))

                psi_angle[i] = psi_angle_out
                psi_spatial_in = psi_spatial_out

        else:
            psi_spatial_in = 0.0  # r=0

            for i in range(nx):
                psi_n_i = psi[n, i]

                psi_spatial_out = 2.0 * psi_n_i - psi_spatial_in

                psi_angle_in = psi_angle[i]
                psi_angle_out = (psi_n_i - (1.0 - tau_n) * psi_angle_in) / tau_n

                A_outer = A[i+1]
                A_inner = A[i]

                # For outward: psi at A_inner = psi_spatial_in, psi at A_outer = psi_spatial_out
                streaming = mu_n * (A_outer * psi_spatial_out - A_inner * psi_spatial_in)

                redist = dAw[i, n] * (alpha_out * psi_angle_out - alpha_in * psi_angle_in)

                collision = sig_t_val * V[i] * psi_n_i
                source = S * V[i]

                residual = streaming + redist + collision - source

                if abs(residual) > 1e-8:
                    print(f"  n={n} cell={i}: stream={streaming:.6e} "
                          f"redist={redist:.6e} coll={collision:.6e} "
                          f"src={source:.6e} RESIDUAL={residual:.6e}")

                max_residual = max(max_residual, abs(residual))

                psi_angle[i] = psi_angle_out
                psi_spatial_in = psi_spatial_out

    print(f"\nMax balance residual: {max_residual:.6e}")
    if max_residual < 1e-8:
        print(f"Converged solution DOES satisfy the balance equation.")
        print(f"=> The balance equation itself does NOT have flat flux as unique solution!")
    else:
        print(f"Converged solution does NOT satisfy the balance equation.")
        print(f"=> The sweep has a convergence bug.")
    print(f"{'='*70}")


def test_balance_for_flat():
    """What if we plug FLAT flux into the balance? Does it satisfy it?

    This tests the DISCRETE balance equation (with DD and WDD closures),
    not the continuous one.
    """
    nx = 5
    N_ord = 4
    R = 5.0

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

    sig_t_val = 1.0
    Q_val = 1.0
    S = Q_val / W
    psi_flat = S / sig_t_val

    print(f"\n{'='*70}")
    print(f"Flat flux balance check (with DD+WDD closures)")
    print(f"psi_flat = {psi_flat:.8f}")

    # For flat flux: all face fluxes (spatial and angular) = psi_flat
    # DD: psi_out = 2*psi - psi_in = 2*psi_flat - psi_flat = psi_flat. OK.
    # WDD: psi_angle_out = (psi - (1-tau)*psi_angle_in) / tau = psi_flat. OK.

    psi_angle_val = psi_flat  # all angular face fluxes are flat

    max_residual = 0.0
    for n in range(N):
        mu_n = mu[n]
        alpha_in = alpha[n]
        alpha_out = alpha[n+1]

        for i in range(nx):
            A_outer = A[i+1]
            A_inner = A[i]

            # All face fluxes = psi_flat
            streaming = mu_n * (A_outer - A_inner) * psi_flat
            redist = dAw[i, n] * (alpha_out - alpha_in) * psi_flat
            collision = sig_t_val * V[i] * psi_flat
            source = S * V[i]

            residual = streaming + redist + collision - source

            max_residual = max(max_residual, abs(residual))
            if abs(residual) > 1e-10:
                print(f"  n={n} (mu={mu_n:+.4f}) cell={i}: "
                      f"stream={streaming:.6e} redist={redist:.6e} "
                      f"coll={collision:.6e} src={source:.6e} "
                      f"RESIDUAL={residual:.6e}")

    print(f"\nMax residual for flat flux: {max_residual:.6e}")
    if max_residual < 1e-10:
        print(f"Flat flux satisfies discrete balance (DD+WDD).")
    else:
        print(f"Flat flux does NOT satisfy discrete balance!")
    print(f"{'='*70}")


if __name__ == "__main__":
    test_balance_for_flat()
    test_balance_residual()
