"""Diagnostic: cylinder refinement convergence post-fix.

Created by numerics-investigator on 2026-05-13.

Compare pre-fix vs post-fix residuals on flat ψ AND end-to-end fixed-source
convergence under mesh refinement n_cells ∈ {20, 40, 80, 160}.
"""
from __future__ import annotations

import numpy as np

from orpheus.derivations.common.xs_library import get_mixture
from orpheus.geometry import BC, Mesh1D, Region, RegionMesh, StructuredGeometry
from orpheus.sn.quadrature import ProductQuadrature
from orpheus.sn.geometry import SNMesh
from orpheus.sn.operator import (
    build_equation_map_cylindrical,
    transport_operator_matvec_cylindrical,
    solution_to_angular_flux_cylindrical,
    SNStreamingOperator,
)
from orpheus.sn.spatial.pole_angular_closure import MorelMontryAngularSweep

import sys
sys.path.insert(0, "derivations/diagnostics")
from diag_phase_g_step2_cyl_carlson_seed_fix import FixedCylinderMM


def test_residual_at_flat_grid():
    fuel = get_mixture("B", "1g")
    R = 2.0
    geom = StructuredGeometry(
        geometry="CYL",
        regions=(Region(mat_id=0, outer_thickness_cm=R),),
        bcs=(BC.reflective,),
    )

    print("\n=== Residual L·ψ − Σ_t·ψ on flat ψ — mesh refinement ===")
    print(f"{'closure':<28} {'n_mu':<5} {'n_phi':<5} "
          + " ".join(f"{n:>12}" for n in [20, 40, 80, 160]))
    print("-" * 100)
    for closure_name, closure_cls in [
        ("MorelMontry (PRE-FIX)", MorelMontryAngularSweep),
        ("FixedCylinderMM (FIX)", FixedCylinderMM),
    ]:
        for n_mu in [2, 4]:
            for n_phi in [2, 4]:
                row = [closure_name, str(n_mu), str(n_phi)]
                for n_cells in [20, 40, 80, 160]:
                    mesh = Mesh1D.from_geometry(
                        geom, region_meshes=(RegionMesh(n_cells=n_cells),),
                    )
                    quad = ProductQuadrature.create(n_mu=n_mu, n_phi=n_phi)
                    sn_mesh = SNMesh(
                        mesh, quad, pole_angular_closure=closure_cls(),
                    )
                    nx = n_cells
                    ng = 1
                    sig_t = np.full((nx, 1, ng), 2.0)
                    psi_flat = 10.0 / quad.weights.sum()
                    eq_map = build_equation_map_cylindrical(nx, quad, ng)
                    sol_flat = np.full(eq_map.n_unknowns, psi_flat)
                    lhs = transport_operator_matvec_cylindrical(
                        sol_flat, eq_map, quad, sig_t, nx, ng,
                        sn_mesh.reduced.face_areas, sn_mesh.volumes,
                        sn_mesh.reduced.alpha_per_level,
                        sn_mesh.reduced.redist_dAw_per_level,
                        sn_mesh.reduced.tau_mm_per_level,
                        sn_mesh=sn_mesh,
                        pole_angular_closure=closure_cls(),
                    )
                    res = np.max(np.abs(lhs - 2.0 * psi_flat))
                    row.append(f"{res:.3e}")
                print(f"{row[0]:<28} {row[1]:<5} {row[2]:<5} "
                      + " ".join(f"{r:>12}" for r in row[3:]))


def test_endtoend_refinement():
    """End-to-end fixed-source solve under refinement (full Krylov with fix)."""
    sys.path.insert(0, "derivations/diagnostics")
    from diag_phase_g_step2_cyl_full_solve_with_fix import _run_krylov_with_fixed_closure
    print("\n\n=== End-to-end Krylov fixed-source with FixedCylinderMM ===")
    print(f"{'n_cells':<10} " + " ".join(
        f"{f'{nmu}x{nphi} err_ψ':<14}" for nmu in [2, 4] for nphi in [2, 4]
    ))
    print("-" * 80)
    for n_cells in [20, 40, 80, 160]:
        row = [str(n_cells)]
        for n_mu in [2, 4]:
            for n_phi in [2, 4]:
                psi, phi, err_psi, err_phi, info = _run_krylov_with_fixed_closure(
                    n_cells, n_mu, n_phi,
                )
                row.append(f"{err_psi:.3e}")
        print(f"{row[0]:<10} " + " ".join(f"{r:<14}" for r in row[1:]))


if __name__ == "__main__":
    test_residual_at_flat_grid()
    test_endtoend_refinement()
