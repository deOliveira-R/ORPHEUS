"""Diagnostic Step 3: replicate L1 anchor exactly, then diff vs production.

Created by numerics-investigator on 2026-05-28.

The L1 anchor ``test_identity_preconditioner_recovers_kinf[sphere]`` PASSES
with the same operators/quadrature/geometry/material but
``SNSolver._solve_krylov`` FAILS (1.402 vs 1.875 = k_inf).  Differences:

  D1. ``inner_tol``     : L1=1e-12, prod=1e-8
  D2. ``max_inner``     : L1=300,   prod=200
  D3. ``keff_tol`` outer: L1=1e-10, prod=1e-12
  D4. Outer rescale     : L1=NONE,  prod=power_iteration divides by P
  D5. RHS wrapper       : L1=L1 AngularFlux,  prod=TimedFullField composite
  D6. q_ext.boundary    : L1=N/A,    prod=zeros_for_sn_mesh
  D7. initial_guess     : L1=psi_typed_warm (LOCAL),
                          prod=self._psi_typed (instance cache)
  D8. keff convergence  : L1 returns once |dk|<tol; prod requires |dk| and
                          flux-residual <tol

Step 1 confirmed: tightening inner_tol alone does NOT fix Krylov
(1.4019239124 at both 1e-8 and 1e-12).  So D1 is NOT causal.

This step toggles D4 (the rescale) inside an otherwise-L1-style loop
and reports.
"""
from __future__ import annotations

import warnings

import numpy as np


def _run_L1_style(*, do_rescale: bool, max_outer: int = 80, keff_tol: float = 1e-10):
    """L1 anchor flow exactly, with optional production-rate rescale."""
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
        geom, region_meshes=(RegionMesh(n_cells=20),),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    sn_mesh = SNMesh(mesh, quad, {0: fuel})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        solver = SNSolver(
            sn_mesh=sn_mesh, scattering_order=0,
            inner_solver="krylov",
            max_inner=300, inner_tol=1e-12,
        )

    L_leaf = StreamingOperator(sn_mesh, solver.mat_xs.total_cross_section)
    C_t = CollisionOperator(sn_mesh, solver.mat_xs.total_cross_section)
    LC = L_leaf + C_t

    N = sn_mesh.quad.N
    nx, ny, ng = sn_mesh.nx, sn_mesh.ny, solver.ng
    krylov = KrylovAcceleration(
        LC, solver.scattering_op, ZeroOperator(),
        preconditioner=lambda q: q,
        tol=1e-12, max_iter=300,
        restart=min(50, N * ng * nx * ny),
    )

    phi = solver.initial_flux_distribution()
    keff = 1.0
    psi_typed_warm: L1AngularFlux | None = None
    keff_history = []
    for n_outer in range(max_outer):
        keff_old = keff
        fis = solver.compute_fission_source(phi, keff)
        q_ext_per_ord = PerOrdinateSource.from_isotropic(fis, sn_mesh)
        q_ext_typed = L1AngularFlux(q_ext_per_ord.values, sn_mesh)
        psi_typed, _residuals = krylov.solve(
            q_ext_typed, initial_guess=psi_typed_warm,
        )
        psi_typed_warm = psi_typed
        phi = psi_typed.integrate_angular().values
        # OPTIONAL production-rate rescale (like power_iteration)
        if do_rescale:
            p = float(solver.compute_production_rate(phi))
            if p > 0.0:
                phi = phi / p
        keff_new = solver.compute_keff(phi)
        keff_history.append(keff_new)
        if abs(keff_new - keff) < keff_tol:
            keff = keff_new
            return keff, n_outer + 1, keff_history
        keff = keff_new
    return keff, max_outer, keff_history


def test_step3_L1_no_rescale_recovers_kinf():
    """Baseline: L1 anchor flow (no rescale) recovers k_inf=1.875."""
    keff, n, _ = _run_L1_style(do_rescale=False)
    print(f"\n[step3-L1] no-rescale: keff={keff:.10f}  err={abs(keff - 1.875):.3e}  n={n}")


def test_step3_L1_with_rescale_keff():
    """Toggle: same L1 flow + production-rate rescale of phi."""
    keff, n, _ = _run_L1_style(do_rescale=True)
    print(f"\n[step3-L1] WITH rescale: keff={keff:.10f}  err={abs(keff - 1.875):.3e}  n={n}")


if __name__ == "__main__":
    test_step3_L1_no_rescale_recovers_kinf()
    test_step3_L1_with_rescale_keff()
