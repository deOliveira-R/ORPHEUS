"""Diagnostic: Investigate the r=0 cell balance equation.

Created by numerics-investigator on 2026-04-16.
The diverging error at cell 0 with refinement suggests the balance
equation is wrong at/near the origin.
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh


def test_origin_cell_geometry():
    """Print geometric quantities near r=0 to check for degeneracies."""
    nx = 10
    R = 10.0
    mesh = Mesh1D(
        edges=np.linspace(0.0, R, nx + 1),
        mat_ids=np.ones(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC.reflective,
        bc_right=BC.reflective,
    )
    quad = GaussLegendre1D.create(8)
    sn_mesh = SNMesh(mesh, quad)

    A = sn_mesh.face_areas
    V = sn_mesh.volumes[:, 0]
    alpha = sn_mesh.alpha_half
    dAw = sn_mesh.redist_dAw
    tau = sn_mesh.tau_mm

    print(f"\n{'='*70}")
    print(f"Geometric quantities for nx={nx}, R={R}")
    print(f"\nFace areas A[i] (should be 4*pi*r^2):")
    for i in range(min(5, len(A))):
        r = mesh.edges[i]
        print(f"  A[{i}] = {A[i]:.6f}  (r={r:.2f}, 4*pi*r^2={4*np.pi*r**2:.6f})")

    print(f"\nVolumes V[i] (should be (4/3)*pi*(r_out^3 - r_in^3)):")
    for i in range(min(5, nx)):
        r_in, r_out = mesh.edges[i], mesh.edges[i+1]
        V_exact = (4.0/3.0) * np.pi * (r_out**3 - r_in**3)
        print(f"  V[{i}] = {V[i]:.6f}  (exact={V_exact:.6f})")

    print(f"\nAlpha coefficients (dome shape):")
    for n in range(len(alpha)):
        print(f"  alpha[{n}+1/2] = {alpha[n]:.6f}")

    print(f"\nTau (Morel-Montry) weights:")
    for n in range(len(tau)):
        print(f"  tau[{n}] = {tau[n]:.6f}")

    print(f"\ndAw[cell=0, ordinate n] = (A_out - A_in) / w_n:")
    for n in range(quad.N):
        print(f"  dAw[0, {n}] = {dAw[0, n]:.6f}  (mu={quad.mu_x[n]:+.6f}, w={quad.weights[n]:.6f})")

    # Key check: for cell 0, A_in = A[0] = 0 (r=0)
    # So dA = A[1] - A[0] = A[1]
    # For the OUTWARD sweep (mu>0), the balance for cell 0:
    # denom = 2*|mu|*A[1] + dAw*c_out + sig_t*V[0]
    # numer = QV[0] + |mu|*(0 + A[1])*0  + dAw*c_in*psi_angle[0]
    #       = QV[0] + 0 + dAw*c_in*psi_angle[0]
    # On the first ordinate (n=0, mu<0, inward sweep):
    # A_in = A[1], A_out = A[0] = 0
    # denom = 2*|mu|*0 + dAw*c_out + sig_t*V[0] = dAw*c_out + sig_t*V[0]
    # numer = QV[0] + |mu|*(A[1]+0)*psi_in + dAw*c_in*psi_angle[0]
    # Note: 2*|mu|*A_out = 0 for cell 0 inward! The spatial streaming OUT
    # at r=0 vanishes because A[0]=0. This is CORRECT physically.

    print(f"\nCell 0 balance check:")
    print(f"  A[0] = {A[0]:.10f} (should be 0)")
    print(f"  A[1] = {A[1]:.6f}")
    print(f"  V[0] = {V[0]:.6f}")
    print(f"  For inward sweep (mu<0) at cell 0: denom includes 2*|mu|*A[0] = 0")
    print(f"  For outward sweep (mu>0) at cell 0: psi_in = 0 (no incoming)")
    print(f"{'='*70}")


def test_single_sweep_balance():
    """Do ONE sweep and check the balance equation residual per cell per ordinate."""
    nx = 10
    R = 10.0
    mesh = Mesh1D(
        edges=np.linspace(0.0, R, nx + 1),
        mat_ids=np.ones(nx, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_left=BC.reflective,
        bc_right=BC.reflective,
    )
    quad = GaussLegendre1D.create(4)  # S4 for clarity
    sn_mesh = SNMesh(mesh, quad)

    N = quad.N
    ng = 1
    sig_t_val = 1.0
    Q_val = 1.0

    A = sn_mesh.face_areas
    V = sn_mesh.volumes[:, 0]
    alpha = sn_mesh.alpha_half
    dAw = sn_mesh.redist_dAw
    tau = sn_mesh.tau_mm
    mu = quad.mu_x
    weights = quad.weights
    W = weights.sum()

    # For a FLAT flux test: psi = const = psi_flat for all ordinates and cells.
    # The balance equation should be satisfied exactly.
    # With psi_flat, the streaming term should be:
    # mu_n * [A_{out}*psi - A_{in}*psi] = mu_n * psi * (A_out - A_in) = mu_n * psi * dA
    # The redistribution term should be:
    # (dA/w_n) * [alpha_{n+1/2}*psi_angle_out - alpha_{n-1/2}*psi_angle_in]
    # For flat flux, if psi_angle is also flat, then:
    # (dA/w_n) * psi * (alpha_{n+1/2} - alpha_{n-1/2}) = (dA/w_n) * psi * (-w_n * mu_n) = -dA * mu_n * psi

    # So: streaming + redistribution = mu_n*psi*dA - dA*mu_n*psi = 0
    # And: Sig_t * V * psi = S * V, giving psi = S/Sig_t = Q/(W*Sig_t)

    psi_flat = Q_val / (W * sig_t_val)
    print(f"\n{'='*70}")
    print(f"Flat flux test: psi_flat = Q/(W*Sig_t) = {psi_flat:.8f}")
    print(f"Checking balance residual per cell per ordinate:")

    max_residual = 0.0
    for n in range(N):
        mu_n = mu[n]
        w_n = weights[n]
        alpha_in = alpha[n]
        alpha_out = alpha[n+1]
        tau_n = tau[n]

        for i in range(nx):
            # Face fluxes for flat case
            A_inner = A[i]
            A_outer = A[i+1]
            dA = A_outer - A_inner

            # Spatial streaming: mu_n * (A_outer * psi - A_inner * psi)
            streaming = mu_n * dA * psi_flat

            # Angular redistribution: (dA/w_n) * (alpha_out * psi_angle_out - alpha_in * psi_angle_in)
            # For flat flux, psi_angle should also be flat.
            # alpha_{n+1/2} - alpha_{n-1/2} = -w_n * mu_n
            # So: (dA/w_n) * psi * (alpha_out - alpha_in) = (dA/w_n) * psi * (-w_n * mu_n)
            #   = -dA * mu_n * psi
            redistribution = -dA * mu_n * psi_flat

            collision = sig_t_val * V[i] * psi_flat
            source = (Q_val / W) * V[i]

            residual = streaming + redistribution + collision - source
            if abs(residual) > max_residual:
                max_residual = abs(residual)

            if abs(residual) > 1e-10:
                print(f"  cell={i}, n={n} (mu={mu_n:+.4f}): "
                      f"stream={streaming:.6e}, redist={redistribution:.6e}, "
                      f"coll={collision:.6e}, src={source:.6e}, "
                      f"RESIDUAL={residual:.6e}")

    print(f"  Max residual: {max_residual:.6e}")
    if max_residual < 1e-10:
        print(f"  PASS: flat flux satisfies the balance equation analytically.")
    else:
        print(f"  FAIL: flat flux does NOT satisfy the balance equation!")

    # Now check: does the WDD balance reproduce this for flat angular flux?
    print(f"\n  Checking WDD closure for flat angular face flux:")
    psi_angle_flat = psi_flat  # flat angular face flux
    for n in range(N):
        alpha_in = alpha[n]
        alpha_out = alpha[n+1]
        tau_n = tau[n]

        # WDD closure: psi_angle_out = (psi - (1-tau)*psi_angle_in) / tau
        psi_angle_out = (psi_flat - (1.0 - tau_n) * psi_angle_flat) / tau_n
        # For flat: should give psi_angle_out = psi_flat
        err = psi_angle_out - psi_flat
        if abs(err) > 1e-15:
            print(f"    n={n}: psi_angle_out = {psi_angle_out:.10f}, expected {psi_flat:.10f}, err={err:.2e}")

    # Check DD closure
    print(f"\n  Checking DD closure for flat spatial face flux:")
    psi_spatial_in = psi_flat
    psi_dd_out = 2.0 * psi_flat - psi_spatial_in
    print(f"    DD: psi_out = 2*psi - psi_in = {psi_dd_out:.10f}, expected {psi_flat:.10f}")

    print(f"{'='*70}")


if __name__ == "__main__":
    test_origin_cell_geometry()
    test_single_sweep_balance()
