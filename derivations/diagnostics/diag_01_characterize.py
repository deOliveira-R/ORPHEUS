"""Diagnostic: Characterize spherical sweep failure for constant source.

Created by numerics-investigator on 2026-04-16.
If this test catches a real bug, promote to tests/sn/test_spherical.py.
"""
import numpy as np
import sys
sys.path.insert(0, "/workspaces/ORPHEUS")

from orpheus.geometry.mesh import Mesh1D, BC
from orpheus.geometry.coord import CoordSystem
from orpheus.sn.quadrature import GaussLegendre1D
from orpheus.sn.geometry import SNMesh
from orpheus.sn.sweep import transport_sweep


def make_spherical_mesh(nx, R=10.0):
    """Uniform spherical mesh from r=0 to r=R."""
    edges = np.linspace(0.0, R, nx + 1)
    mat_ids = np.ones(nx, dtype=int)
    return Mesh1D(
        edges=edges,
        mat_ids=mat_ids,
        coord=CoordSystem.SPHERICAL,
        bc_left=BC.reflective,
        bc_right=BC.reflective,
    )


def run_fixed_source_sweep(nx, N_ord, sig_t_val=1.0, Q_val=1.0, c=0.5,
                           max_iter=200, tol=1e-12, R=10.0):
    """Run source iteration with the spherical sweep.

    Constant isotropic source Q, constant Sig_t, scattering ratio c.
    Reflective BCs at outer, r=0 at inner.
    Expected: phi = Q / (Sig_t * (1-c)) everywhere.
    """
    mesh = make_spherical_mesh(nx, R)
    quad = GaussLegendre1D.create(N_ord)
    sn_mesh = SNMesh(mesh, quad)

    ng = 1
    N = quad.N

    sig_t = np.full((nx, 1, ng), sig_t_val)
    sig_s_val = c * sig_t_val

    phi = np.ones((nx, 1, ng))
    psi_bc = {}

    expected_phi = Q_val / (sig_t_val * (1.0 - c))

    residuals = []
    for it in range(max_iter):
        phi_old = phi.copy()

        # Isotropic source = external + scattering
        Q_iso = np.full((nx, 1, ng), Q_val) + sig_s_val * phi

        _, phi = transport_sweep(Q_iso, sig_t, sn_mesh, psi_bc)

        norm = np.linalg.norm(phi)
        if norm > 0:
            res = np.linalg.norm(phi - phi_old) / norm
            residuals.append(res)
            if res < tol:
                break

    return phi[:, 0, 0], expected_phi, it + 1, residuals


def test_characterize_spherical_sweep():
    """Characterize the failure: constant source, reflective BCs, spherical."""
    nx = 20
    N_ord = 8

    phi, expected, n_iter, residuals = run_fixed_source_sweep(nx, N_ord)

    print(f"\n{'='*60}")
    print(f"Spherical sweep: constant source diagnostic")
    print(f"  nx={nx}, S{N_ord}, reflective BCs")
    print(f"  Expected phi = {expected:.6f}")
    print(f"  Converged in {n_iter} iterations")
    if residuals:
        print(f"  Final residual: {residuals[-1]:.2e}")
    print(f"  phi range: [{phi.min():.6f}, {phi.max():.6f}]")
    print(f"  phi mean:  {phi.mean():.6f}")
    print(f"  Rel err range: [{(phi.min()/expected - 1)*100:.2f}%, {(phi.max()/expected - 1)*100:.2f}%]")

    for i in range(min(5, len(phi))):
        print(f"    cell {i:3d}: phi={phi[i]:.8f}  err={phi[i]/expected - 1:+.6f}")
    print(f"    ...")
    for i in range(max(len(phi)-5, 5), len(phi)):
        print(f"    cell {i:3d}: phi={phi[i]:.8f}  err={phi[i]/expected - 1:+.6f}")
    print(f"{'='*60}")


def test_no_scattering():
    """c=0 — single sweep should give phi=Q/Sig_t (but BCs need iteration)."""
    nx = 20
    N_ord = 8

    phi, expected, n_iter, residuals = run_fixed_source_sweep(
        nx, N_ord, c=0.0, max_iter=100
    )

    print(f"\n{'='*60}")
    print(f"Spherical sweep: c=0 (no scattering)")
    print(f"  Expected phi = {expected:.6f}")
    print(f"  n_iter = {n_iter}")
    print(f"  phi range: [{phi.min():.6f}, {phi.max():.6f}]")
    print(f"  phi mean:  {phi.mean():.6f}")
    for i in range(min(5, len(phi))):
        print(f"    cell {i:3d}: phi={phi[i]:.8f}  err={phi[i]/expected - 1:+.6f}")
    print(f"{'='*60}")


def test_scaling():
    """Error vs refinement."""
    N_ord = 8
    c = 0.5
    expected = 1.0 / (1.0 * (1.0 - c))

    print(f"\n{'='*60}")
    print(f"Scaling analysis: spherical sweep (expected={expected})")
    print(f"  {'Cells':>6s}  {'Mean':>10s}  {'Min':>10s}  {'Max':>10s}  {'RelErr':>10s}  {'Ratio':>8s}")

    prev_err = None
    for nx in [5, 10, 20, 40]:
        phi, _, n_iter, _ = run_fixed_source_sweep(nx, N_ord)
        err = abs(phi.mean() / expected - 1)
        ratio = prev_err / err if prev_err is not None and err > 0 else float('nan')
        print(f"  {nx:6d}  {phi.mean():10.6f}  {phi.min():10.6f}  {phi.max():10.6f}  {err:10.6e}  {ratio:8.2f}")
        prev_err = err

    print(f"{'='*60}")


if __name__ == "__main__":
    test_characterize_spherical_sweep()
    test_no_scattering()
    test_scaling()
