"""Diagnostic: full Krylov solve on cylinder with patched Carlson seed.

Created by numerics-investigator on 2026-05-13.

Plug the fixed-norm MM closure into SN solver and run end-to-end.
Both Krylov and SI need separate fixes (the SI sweep has its own
Carlson seed branch).
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import ProductQuadrature
from orpheus.sn.geometry import SNMesh

# Import the patched closure
import sys
sys.path.insert(0, "derivations/diagnostics")
from diag_phase_g_step2_cyl_carlson_seed_fix import FixedCylinderMM


def _run_krylov_with_fixed_closure(n_cells, n_mu, n_phi):
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
    sn_mesh = SNMesh(mesh, quad, pole_angular_closure=FixedCylinderMM())

    # Build the SN solver pieces manually (use SNStreamingOperator + scipy Krylov)
    from orpheus.sn.operator import (
        SNStreamingOperator,
        build_equation_map_cylindrical,
        solution_to_angular_flux_cylindrical,
    )

    ng = 1
    sig_t = np.full((n_cells, 1, ng), 2.0)
    sigma_s = 1.9
    eq_map = build_equation_map_cylindrical(n_cells, quad, ng)

    L = SNStreamingOperator(sn_mesh=sn_mesh, sig_t=sig_t)
    n_unk = L.n_unknowns

    from scipy.sparse.linalg import LinearOperator, gmres

    def matvec(psi_vec):
        # (L - S) ψ; S is isotropic scattering: S·ψ = (Σ_s/Σw) · Σ_n w_n ψ_n
        # applied to each ordinate
        lpsi = L.apply(psi_vec)
        # Build scalar flux from solution vector
        fi = solution_to_angular_flux_cylindrical(
            psi_vec, eq_map, quad, n_cells, ng,
        )
        phi_0 = np.einsum("n,gnxy->gxy", quad.weights, fi)  # (ng, nx, 1)
        # Scattering source per-ordinate at cell i: (Σ_s/Σw) · φ_0
        Spsi = np.zeros_like(fi)
        for n in range(quad.N):
            Spsi[:, n, :, 0] = sigma_s * phi_0[:, :, 0] / quad.weights.sum()
        # Pack back into solution vector
        Spsi_vec = np.zeros_like(psi_vec)
        flux_view = Spsi_vec.reshape(ng, eq_map.n_eq, order="F")
        for k in range(eq_map.n_eq):
            flux_view[:, k] = Spsi[:, eq_map.ordinate[k], eq_map.ix[k], 0]
        return lpsi - Spsi_vec

    A = LinearOperator((n_unk, n_unk), matvec=matvec, dtype=float)
    # RHS: per-ordinate external source Q_n = Q_ext / Σw (W normalisation)
    Q_per_ord = 1.0 / quad.weights.sum()
    rhs = np.full(n_unk, Q_per_ord)

    psi_vec, info = gmres(A, rhs, rtol=1e-12, atol=1e-12, maxiter=200)
    fi = solution_to_angular_flux_cylindrical(psi_vec, eq_map, quad, n_cells, ng)
    psi = fi[0, :, :, 0]                 # (N, nx)
    phi = np.einsum("n,nx->x", quad.weights, psi)  # (nx,)

    psi_ref = 10.0 / quad.weights.sum()
    err_psi = np.max(np.abs(psi - psi_ref))
    err_phi = np.max(np.abs(phi - 10.0))
    return psi, phi, err_psi, err_phi, info


def test_full_solve_with_fix():
    print("\n=== Cylinder fixed-source Krylov solve with patched Carlson seed ===")
    print(f"Expected ψ_n = 0.7957747 per ordinate per cell, φ = 10 per cell\n")
    for n_mu in [2, 4]:
        for n_phi in [2, 4]:
            print(f"\n  ProductQuadrature(n_mu={n_mu}, n_phi={n_phi})")
            for n_cells in [5, 10, 20, 40, 80]:
                try:
                    psi, phi, err_psi, err_phi, info = _run_krylov_with_fixed_closure(
                        n_cells, n_mu, n_phi,
                    )
                    print(f"    n_cells={n_cells:3d}: "
                          f"err_ψ = {err_psi:.3e}, err_φ = {err_phi:.3e}, "
                          f"GMRES info={info}")
                except Exception as e:
                    print(f"    n_cells={n_cells:3d}: EXCEPTION {e!s}")


if __name__ == "__main__":
    test_full_solve_with_fix()
