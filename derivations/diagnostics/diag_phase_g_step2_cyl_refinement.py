"""Diagnostic: cylinder Phase G Step 2 convergence-under-refinement.

Created by numerics-investigator on 2026-05-13.

Probe whether SI and Krylov converge to the analytical ψ_n = φ/Σw = 0.7957747
under mesh refinement on the homogeneous reflective cylinder L0 test.
A mesh-independent residual is the structural-defect fingerprint.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import ProductQuadrature
from orpheus.sn.solver import solve_sn_fixed_source


def _run(inner_solver, n_cells, n_mu, n_phi):
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
    Q = np.ones((quad.N, mesh.N, 1, 1))
    result = solve_sn_fixed_source(
        materials={0: fuel}, mesh=mesh, quadrature=quad, external_source=Q,
        boundary_condition="reflective", inner_solver=inner_solver,
    )
    psi = result.angular_flux[:, :, 0, :]
    sf = result.scalar_flux[:, 0, :]
    return psi, sf, quad


def test_refinement():
    n_mu, n_phi = 2, 2
    print(f"\nProductQuadrature(n_mu={n_mu}, n_phi={n_phi}); N=4")
    print(f"Analytical ψ_n = 10/(4π) = {10/(4*np.pi):.10f}")
    print(f"Analytical φ   = 10")
    print()
    print(f"{'n_cells':<8} {'err_SI':<14} {'err_K':<14} "
          f"{'err_SI_vs_K':<14}  "
          f"phi(i=0)_SI  phi(i=0)_K")
    print("-" * 90)
    for n in [2, 5, 10, 20, 40, 80]:
        try:
            psi_si, sf_si, quad = _run("source_iteration", n, n_mu, n_phi)
            psi_k, sf_k, _ = _run("krylov", n, n_mu, n_phi)
            psi_ref = 10.0 / quad.weights.sum()
            err_si = np.max(np.abs(psi_si - psi_ref))
            err_k = np.max(np.abs(psi_k - psi_ref))
            err_si_vs_k = np.max(np.abs(psi_si - psi_k))
            print(f"{n:<8} {err_si:<14.6e} {err_k:<14.6e} "
                  f"{err_si_vs_k:<14.6e}  "
                  f"{sf_si[0,0]:.4f}     {sf_k[0,0]:.4f}")
        except Exception as e:
            print(f"{n:<8} EXCEPTION: {e!s}")


if __name__ == "__main__":
    test_refinement()
