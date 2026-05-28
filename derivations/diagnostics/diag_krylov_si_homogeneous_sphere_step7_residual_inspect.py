"""Diagnostic Step 7: residual inspection — is GMRES actually converging?

Created by numerics-investigator on 2026-05-28.

Step 6 found tight tolerances rescue Krylov AT n_cells=10 but NOT at
n_cells=20.  Hypothesis: at n=20 GMRES hits the max_inner=1000 cap
without converging, returns an unconverged result that bypasses the
``rtol`` decision.

This probe runs the SAME flow as `_solve_krylov` but exposes the
residual history per outer iteration AND the scipy.gmres ``info`` flag.
"""
from __future__ import annotations

import warnings

import numpy as np


def _run_krylov_inspect(*, n_cells, inner_tol, max_inner, max_outer=80):
    """Run manual Krylov power iteration with residual reporting."""
    from orpheus.derivations.common.xs_library import get_mixture
    from orpheus.geometry import (
        BC, Mesh1D, Region, RegionMesh, StructuredGeometry,
    )
    from orpheus.numerics.iteration import KrylovAcceleration
    from orpheus.numerics.operator import ZeroOperator
    from orpheus.numerics.quadrature import Quadrature
    from orpheus.sn.angular_flux import AngularFlux as L1AngularFlux
    from orpheus.sn.geometry import SNMesh
    from orpheus.sn.operator import CollisionOperator, StreamingOperator
    from orpheus.sn.solver import SNSolver
    from orpheus.transport.sources import PerOrdinateSource

    fuel = get_mixture("A", "2g")
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(Region(mat_id=0, outer_thickness_cm=2.0),),
        bcs=(BC.reflective,),
    )
    mesh = Mesh1D.from_geometry(
        geom, region_meshes=(RegionMesh(n_cells=n_cells),),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    sn_mesh = SNMesh(mesh, quad, {0: fuel})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        solver = SNSolver(
            sn_mesh=sn_mesh, scattering_order=0, inner_solver="krylov",
            max_inner=max_inner, inner_tol=inner_tol,
        )

    L_leaf = StreamingOperator(sn_mesh, solver.mat_xs.total_cross_section)
    C_t = CollisionOperator(sn_mesh, solver.mat_xs.total_cross_section)
    LC = L_leaf + C_t
    N = sn_mesh.quad.N
    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, solver.ng
    krylov = KrylovAcceleration(
        LC, solver.scattering_op, ZeroOperator(),
        preconditioner=lambda q: q,
        tol=inner_tol, max_iter=max_inner,
        restart=min(50, N * ng * nx * ny),
    )

    phi = solver.initial_flux_distribution()
    keff = 1.0
    psi_typed_warm = None
    for n_outer in range(max_outer):
        fis = solver.compute_fission_source(phi, keff)
        q_ext_per_ord = PerOrdinateSource.from_isotropic(fis, sn_mesh)
        q_ext_typed = L1AngularFlux(q_ext_per_ord.values, sn_mesh)
        psi_typed, residuals = krylov.solve(
            q_ext_typed, initial_guess=psi_typed_warm,
        )
        psi_typed_warm = psi_typed
        phi = psi_typed.integrate_angular().values
        keff_new = solver.compute_keff(phi)
        n_inner = len(residuals)
        last_res = residuals[-1] if residuals else float("nan")
        print(
            f"  outer={n_outer:3d}  keff={keff_new:.10f}  "
            f"|dk|={abs(keff_new-keff):.3e}  "
            f"n_inner={n_inner:5d}  last_res={last_res:.3e}"
        )
        if abs(keff_new - keff) < 1e-10:
            return keff_new, n_outer + 1
        keff = keff_new
    return keff, max_outer


def test_step7_n10_tight():
    print("\n--- n=10, tol=1e-12, max_inner=1000 ---")
    keff, n = _run_krylov_inspect(n_cells=10, inner_tol=1e-12, max_inner=1000, max_outer=30)
    print(f"FINAL: keff={keff:.10f}  err={abs(keff - 1.875):.3e}  n_outer={n}")


def test_step7_n20_tight():
    print("\n--- n=20, tol=1e-12, max_inner=1000 ---")
    keff, n = _run_krylov_inspect(n_cells=20, inner_tol=1e-12, max_inner=1000, max_outer=30)
    print(f"FINAL: keff={keff:.10f}  err={abs(keff - 1.875):.3e}  n_outer={n}")


def test_step7_n20_loose():
    print("\n--- n=20, tol=1e-8, max_inner=200 (production default) ---")
    keff, n = _run_krylov_inspect(n_cells=20, inner_tol=1e-8, max_inner=200, max_outer=30)
    print(f"FINAL: keff={keff:.10f}  err={abs(keff - 1.875):.3e}  n_outer={n}")


if __name__ == "__main__":
    test_step7_n10_tight()
    test_step7_n20_tight()
    test_step7_n20_loose()
