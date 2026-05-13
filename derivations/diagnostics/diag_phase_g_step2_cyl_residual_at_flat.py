"""Diagnostic: cylinder apply-matvec residual at the candidate flat fixed point.

Created by numerics-investigator on 2026-05-13.

Question: if we feed ψ = analytical-flat to the cylindrical apply-matvec
operator, what residual does it produce? Compare to sphere.

If sphere residual ≈ 0 but cylinder residual >> 0, the cylinder apply-matvec
has a structural defect even before fixed-source convergence.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import GaussLegendre1D, ProductQuadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    SNStreamingOperator,
    build_equation_map_spherical,
    build_equation_map_cylindrical,
)


def _build_curvilinear(geom_type, n_cells, n_ord_or_mu, n_phi=4):
    fuel = get_mixture("B", "1g")
    R = 2.0
    geom = StructuredGeometry(
        geometry=geom_type,
        regions=(Region(mat_id=0, outer_thickness_cm=R),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=n_cells),),
    )
    if geom_type == "SPH":
        quad = GaussLegendre1D.create(n_ordinates=n_ord_or_mu)
    else:
        quad = ProductQuadrature.create(n_mu=n_ord_or_mu, n_phi=n_phi)
    sn_mesh = SNMesh(mesh, quad)
    return fuel, mesh, quad, sn_mesh


def _evaluate_residual(geom_type, n_cells, n_ord_or_mu, n_phi=4):
    """Compute L·ψ_flat for both sphere and cylinder; compare to Σ_t·ψ."""
    fuel, mesh, quad, sn_mesh = _build_curvilinear(
        geom_type, n_cells, n_ord_or_mu, n_phi,
    )
    nx = mesh.N
    ng = 1
    N = quad.N

    # Cross-section field
    sig_t = np.full((nx, 1, ng), 2.0)
    sigma_a = 0.1
    Q_ext = 1.0
    Sw = quad.weights.sum()
    psi_flat_val = (Q_ext / sigma_a) / Sw

    # Build the packed solution vector with ψ = ψ_flat everywhere
    if geom_type == "SPH":
        eq_map = build_equation_map_spherical(nx, quad, ng)
    else:
        eq_map = build_equation_map_cylindrical(nx, quad, ng)

    sol_flat = np.full(eq_map.n_unknowns, psi_flat_val)

    # Apply L via SNStreamingOperator
    L = SNStreamingOperator(sn_mesh=sn_mesh, sig_t=sig_t)
    lhs = L.apply(sol_flat)

    # Expected: at flat ψ, L = streaming + redistribution + collision.
    # For flat ψ: ∇·(Ωψ) = ψ ∇·Ω = 0 (free streaming)
    #             angular redistribution sum should also be 0
    #             collision is Σ_t·ψ
    # So L·ψ_flat should equal Σ_t·ψ_flat = 2 · 0.7957747 = 1.5915
    # or for sphere: 2 · 5 = 10.0
    expected_lhs = 2.0 * psi_flat_val

    residual = lhs - expected_lhs

    return lhs, expected_lhs, residual, eq_map, sn_mesh


def test_residuals_at_flat():
    print("\n=== Residual L·ψ_flat − Σ_t·ψ_flat ===")
    print("If geometry+redistribution operators are consistent, this should be ~0.")
    print()
    print("--- SPHERE ---")
    for n_ord in [4, 8]:
        for n_cells in [10, 20, 40]:
            lhs, exp, res, eq_map, sn_mesh = _evaluate_residual(
                "SPH", n_cells=n_cells, n_ord_or_mu=n_ord,
            )
            psi_flat = exp / 2.0
            print(f"  n_ord={n_ord}, n_cells={n_cells}: "
                  f"||res||_inf = {np.max(np.abs(res)):.4e}  "
                  f"(rel: {np.max(np.abs(res))/psi_flat/2.0:.2e})")
    print()
    print("--- CYLINDER ---")
    for n_mu in [2, 4]:
        for n_phi in [2, 4]:
            for n_cells in [10, 20, 40]:
                lhs, exp, res, eq_map, sn_mesh = _evaluate_residual(
                    "CYL", n_cells=n_cells, n_ord_or_mu=n_mu, n_phi=n_phi,
                )
                psi_flat = exp / 2.0
                print(f"  n_mu={n_mu}, n_phi={n_phi}, n_cells={n_cells}: "
                      f"||res||_inf = {np.max(np.abs(res)):.4e}  "
                      f"(rel: {np.max(np.abs(res))/psi_flat/2.0:.2e})")
    print()
    print("--- CYLINDER detailed (n_mu=2 n_phi=2 n_cells=5) ---")
    lhs, exp, res, eq_map, sn_mesh = _evaluate_residual(
        "CYL", n_cells=5, n_ord_or_mu=2, n_phi=2,
    )
    from orpheus.sn.operator import solution_to_angular_flux_cylindrical
    fi_lhs = solution_to_angular_flux_cylindrical(lhs, eq_map, sn_mesh.quad, 5, 1)
    psi_flat = exp / 2.0
    print(f"  expected: every entry = {2.0*psi_flat:.6f} = Σ_t·ψ_flat")
    print(f"  L·ψ at flat (N=4 rows × nx=5 cols):")
    print(fi_lhs[0, :, :, 0])
    print(f"\n  residual map (L·ψ − Σ_t·ψ):")
    print(fi_lhs[0, :, :, 0] - 2.0 * psi_flat)


if __name__ == "__main__":
    test_residuals_at_flat()
