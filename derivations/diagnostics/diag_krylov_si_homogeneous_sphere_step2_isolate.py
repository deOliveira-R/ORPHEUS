"""Diagnostic Step 2: factor-isolate the SI-vs-Krylov mismatch.

Created by numerics-investigator on 2026-05-28.

Replicates the L1 anchor's manual outer power loop (no production_rate
rescale, manual warm start) using BOTH inner solvers, varying:

  axis-A: production-rate rescale ON / OFF
  axis-B: warm-start cache ON / OFF
  axis-C: inner_tol = 1e-8 / 1e-12

This is the cross of the three suspect factors from the brief.

The L1 anchor `test_identity_preconditioner_recovers_kinf[sphere]` passes
at rescale=OFF, warm-start=local-no-clear, tol=1e-12.  So we already
know `(OFF, local, 1e-12)` works.

If this test catches a real bug, promote to ``tests/sn/test_sn_solver.py``.
"""
from __future__ import annotations

import warnings

import numpy as np


def _build_problem():
    """Sphere, mixture A 2g, n_cells=20, GL n_ord=8."""
    from orpheus.derivations.common.xs_library import get_mixture
    from orpheus.geometry import (
        BC, Mesh1D, Region, RegionMesh, StructuredGeometry,
    )
    from orpheus.numerics.quadrature import Quadrature

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
    return {0: fuel}, mesh, quad


def _manual_outer(
    *,
    inner_solver: str,
    inner_tol: float,
    do_rescale: bool,
    warm_start: bool,
    max_outer: int = 80,
    keff_tol: float = 1e-10,
):
    """Manual outer power loop with explicit rescale / warm-start switches.

    Replicates ``power_iteration`` minus optional production-rate rescale,
    optional warm-start cache.  Returns (keff, n_outer, keff_history,
    inner_iter_counts).
    """
    from orpheus.sn.mesh.augmented_mesh import SNMesh
    from orpheus.sn.solver import SNSolver

    materials, mesh, quad = _build_problem()
    sn_mesh = SNMesh(mesh, quad, materials)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        solver = SNSolver(
            sn_mesh=sn_mesh, scattering_order=0,
            inner_solver=inner_solver,
            keff_tol=keff_tol, flux_tol=1e-10,
            max_inner=300, inner_tol=inner_tol,
        )

    phi = solver.initial_flux_distribution()
    keff = 1.0
    keff_history = []
    for n_outer in range(max_outer):
        keff_old = keff
        fis = solver.compute_fission_source(phi, keff)
        # solve_fixed_source uses solver._psi_typed as warm-start cache.
        # If warm_start=False, clear it before each call.
        if not warm_start:
            solver._psi_typed = None
        phi = solver.solve_fixed_source(fis, phi)
        if do_rescale:
            p = float(solver.compute_production_rate(phi))
            if p > 0:
                phi = phi / p
        keff = solver.compute_keff(phi)
        keff_history.append(keff)
        if abs(keff - keff_old) < keff_tol:
            return keff, n_outer + 1, keff_history

    return keff, max_outer, keff_history


def test_step2_si_passes_at_all_combinations():
    """SI is robust to rescale/warm-start variations."""
    kinf = 1.875

    combos = [
        ("rescale=ON  warm=ON  tol=1e-8 ", True, True, 1e-8),
        ("rescale=ON  warm=ON  tol=1e-12", True, True, 1e-12),
        ("rescale=OFF warm=ON  tol=1e-8 ", False, True, 1e-8),
        ("rescale=OFF warm=OFF tol=1e-8 ", False, False, 1e-8),
    ]

    print()
    for label, rescale, warm, tol in combos:
        keff, n, _ = _manual_outer(
            inner_solver="source_iteration",
            inner_tol=tol, do_rescale=rescale, warm_start=warm,
        )
        err = abs(keff - kinf)
        print(f"[step2-SI ] {label}  keff={keff:.10f}  err={err:.3e}  n={n}")


def test_step2_krylov_factor_isolation():
    """Cross product of (rescale, warm-start, inner_tol) for Krylov."""
    kinf = 1.875

    combos = [
        # Each row: (label, do_rescale, warm_start, inner_tol)
        ("rescale=ON  warm=ON  tol=1e-8 ", True,  True,  1e-8),   # production default
        ("rescale=ON  warm=ON  tol=1e-12", True,  True,  1e-12),
        ("rescale=OFF warm=ON  tol=1e-8 ", False, True,  1e-8),
        ("rescale=OFF warm=ON  tol=1e-12", False, True,  1e-12),  # ≈ L1 anchor
        ("rescale=ON  warm=OFF tol=1e-8 ", True,  False, 1e-8),
        ("rescale=ON  warm=OFF tol=1e-12", True,  False, 1e-12),
        ("rescale=OFF warm=OFF tol=1e-8 ", False, False, 1e-8),
        ("rescale=OFF warm=OFF tol=1e-12", False, False, 1e-12),
    ]

    print()
    for label, rescale, warm, tol in combos:
        keff, n, _ = _manual_outer(
            inner_solver="krylov",
            inner_tol=tol, do_rescale=rescale, warm_start=warm,
        )
        err = abs(keff - kinf)
        print(f"[step2-KR ] {label}  keff={keff:.10f}  err={err:.3e}  n={n}")


if __name__ == "__main__":
    test_step2_si_passes_at_all_combinations()
    test_step2_krylov_factor_isolation()
