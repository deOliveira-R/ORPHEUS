"""Diagnostic: cylinder Phase G Step 2 minimal reproducer.

Created by numerics-investigator on 2026-05-13.

Minimal cylinder configuration that demonstrates the 3-way disagreement:

- SI ``source_iteration`` path
- Krylov ``apply-matvec`` path
- Analytical closed-form (Pomraning isotropy on homogeneous reflective cylinder)

Mixture B 1G (Σ_t=2, Σ_s=1.9, c=0.95). Reflective cylinder R=2 cm. Isotropic
source Q=1 → φ_analytic = 10. With ProductQuadrature(n_mu=2, n_phi=2),
Σw = 4π, so ψ_n_analytic = 10/(4π) ≈ 0.7957747 for every ordinate.

If this diagnostic finds disagreement between SI vs analytical OR Krylov
vs analytical OR SI vs Krylov, the path containing the larger error is
the defective production code.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import ProductQuadrature
from orpheus.sn.solver import solve_sn_fixed_source


def _build(n_cells: int, n_mu: int, n_phi: int):
    fuel = get_mixture("B", "1g")
    R = 2.0
    geom = StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=R),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=n_cells),),
    )
    quad = ProductQuadrature.create(n_mu=n_mu, n_phi=n_phi)
    N = quad.N
    nx = mesh.N
    ng = 1
    Q = np.ones((N, nx, 1, ng))
    return fuel, mesh, quad, Q, N, nx, ng


def _run(inner_solver: str, n_cells=2, n_mu=2, n_phi=2):
    fuel, mesh, quad, Q, N, nx, ng = _build(n_cells, n_mu, n_phi)
    result = solve_sn_fixed_source(
        materials={0: fuel}, mesh=mesh, quadrature=quad, external_source=Q,
        boundary_condition="reflective", inner_solver=inner_solver,
    )
    psi = result.angular_flux[:, :, 0, :]  # (N, nx, ng)
    sf = result.scalar_flux[:, 0, :]       # (nx, ng)
    return psi, sf, quad


def test_minimal_2x4_three_way():
    """2x4 (n_cells=2, n_mu=2, n_phi=2 → N=4) cylinder probe."""
    psi_si, sf_si, quad = _run("source_iteration", 2, 2, 2)
    psi_k, sf_k, _ = _run("krylov", 2, 2, 2)

    psi_ref = 10.0 / quad.weights.sum()  # 10 / 4π
    phi_ref = 10.0

    print("\n=== Cylinder Phase G Step 2 minimal 2x4 (n_cells=2, N=4) ===")
    print(f"weights = {quad.weights}")
    print(f"weights.sum() = {quad.weights.sum():.6f} (= 4π = {4*np.pi:.6f})")
    print(f"mu_x (radial cosines) = {quad.mu_x}")
    print(f"mu_z (axial cosines)  = {quad.mu_z}")
    print(f"n_levels = {quad.n_levels}")
    print(f"level_indices = {quad.level_indices}")
    print(f"level_mu (z) = {quad.level_mu}")

    print(f"\nAnalytical ψ_n per ordinate per cell = {psi_ref:.6f}")
    print(f"Analytical φ per cell = {phi_ref:.6f}")

    print(f"\nSI ψ shape = {psi_si.shape}")
    print(f"SI ψ[:, :, 0] (N x nx):")
    print(psi_si[:, :, 0])
    print(f"SI φ[:, 0] (nx,): {sf_si[:, 0]}")

    print(f"\nKrylov ψ[:, :, 0]:")
    print(psi_k[:, :, 0])
    print(f"Krylov φ[:, 0]: {sf_k[:, 0]}")

    err_si = np.max(np.abs(psi_si - psi_ref))
    err_k = np.max(np.abs(psi_k - psi_ref))
    err_si_vs_k = np.max(np.abs(psi_si - psi_k))

    print(f"\nmax|ψ_SI - ψ_ref|     = {err_si:.6e}   (rel: {err_si/psi_ref:.3%})")
    print(f"max|ψ_K - ψ_ref|      = {err_k:.6e}   (rel: {err_k/psi_ref:.3%})")
    print(f"max|ψ_SI - ψ_K|       = {err_si_vs_k:.6e}")

    # Both should converge to the analytical answer at machine precision.
    # If either does not, that path is defective.
    print(f"\nSI converges to analytical? {err_si < 1e-9}")
    print(f"Krylov converges to analytical? {err_k < 1e-9}")
    print(f"SI agrees with Krylov?         {err_si_vs_k < 1e-9}")


if __name__ == "__main__":
    test_minimal_2x4_three_way()
