"""Diagnostic Step 8: GMRES restart size and info-flag check.

Created by numerics-investigator on 2026-05-28.

Step 7 showed at n_cells=20 GMRES residuals saturate around 6.5 — i.e.,
**GMRES is not converging**.  Production
:class:`~orpheus.numerics.iteration.KrylovAcceleration` (a) hardcodes
``restart=min(50, N*ng*nx*ny)`` (line 716 of solver.py + iteration.py),
and (b) DISCARDS the scipy ``info`` flag (``solution, _info = spla.gmres(...)``,
line 735 of iteration.py).

This step toggles restart size to see if the GMRES result actually
converges given enough Krylov subspace, AND inspects the scipy info
flag from a wrapped call.
"""
from __future__ import annotations

import warnings

import numpy as np


def _run_krylov_with_restart(*, n_cells, restart, inner_tol=1e-12, max_inner=2000, max_outer=30):
    """Manual Krylov with explicit restart size."""
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
    krylov = KrylovAcceleration(
        LC, solver.scattering_op, ZeroOperator(),
        preconditioner=lambda q: q,
        tol=inner_tol, max_iter=max_inner,
        restart=restart,
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
        if abs(keff_new - keff) < 1e-10:
            return keff_new, n_outer + 1, len(residuals), residuals[-1] if residuals else float("nan")
        keff = keff_new
    return keff, max_outer, len(residuals), residuals[-1] if residuals else float("nan")


def test_step8_restart_sweep_n20():
    """Sweep GMRES restart at n_cells=20 to see if convergence is reachable."""
    kinf = 1.875
    print()
    print(f"  {'restart':>8s} {'keff':>15s} {'err':>10s} {'last_res':>10s}")
    # The full angular-flux vector size at n_cells=20, N=8, ng=2:
    # bulk = 8*2*20 = 320 + boundary face slots = ~16 → ~336.
    # Try restart values that should encompass full Krylov subspace.
    for restart in [50, 100, 200, 400, 800]:
        keff, n_out, n_inner, last_res = _run_krylov_with_restart(
            n_cells=20, restart=restart,
        )
        print(f"  {restart:>8d} {keff:>15.10f} {abs(keff - kinf):>10.3e} {last_res:>10.3e}")


def test_step8_info_flag_inspect():
    """Inspect scipy.gmres info flag directly for n=20 production-default case."""
    import scipy.sparse.linalg as spla
    from orpheus.derivations.common.xs_library import get_mixture
    from orpheus.geometry import (
        BC, Mesh1D, Region, RegionMesh, StructuredGeometry,
    )
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
        geom, region_meshes=(RegionMesh(n_cells=20),),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    sn_mesh = SNMesh(mesh, quad, {0: fuel})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        solver = SNSolver(
            sn_mesh=sn_mesh, scattering_order=0, inner_solver="krylov",
            max_inner=200, inner_tol=1e-8,
        )
    L_leaf = StreamingOperator(sn_mesh, solver.mat_xs.total_cross_section)
    C_t = CollisionOperator(sn_mesh, solver.mat_xs.total_cross_section)
    LC = L_leaf + C_t

    phi = solver.initial_flux_distribution()
    fis = solver.compute_fission_source(phi, 1.0)
    q_ext_per_ord = PerOrdinateSource.from_isotropic(fis, sn_mesh)
    q_ext_typed = L1AngularFlux(q_ext_per_ord.values, sn_mesh)

    # Build raw scipy GMRES on (L+C-S)
    def matvec(psi_flat):
        psi_t = L1AngularFlux.from_flat_with_traces(psi_flat, sn_mesh)
        # (L+C-S) psi
        Lpsi = LC.apply(psi_t)
        Spsi = solver.scattering_op.apply(psi_t)
        out = Lpsi - Spsi
        return out.to_flat_with_traces()

    b_flat = q_ext_typed.to_flat_with_traces()
    n_flat = b_flat.size
    A = spla.LinearOperator((n_flat, n_flat), matvec=matvec, dtype=float)

    print()
    print(f"n_flat (full ravel size) = {n_flat}")
    print(f"||b||_2 = {np.linalg.norm(b_flat):.3e}")
    print()
    for restart in [50, 100, 200, 320]:
        residuals = []
        sol, info = spla.gmres(
            A, b_flat, rtol=1e-12, atol=0.0,
            maxiter=200, restart=min(restart, n_flat),
            callback=lambda r: residuals.append(float(np.asarray(r))),
            callback_type='pr_norm',
        )
        true_res = np.linalg.norm(b_flat - matvec(sol))
        print(
            f"  restart={restart:4d}  info={info:4d}  n_callbacks={len(residuals):5d}  "
            f"true_||r||={true_res:.3e}  last_pr_norm={residuals[-1] if residuals else float('nan'):.3e}"
        )


if __name__ == "__main__":
    test_step8_info_flag_inspect()
    test_step8_restart_sweep_n20()
