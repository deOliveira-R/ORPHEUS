"""Diagnostic: characterize the curvilinear MMS convergence-rate gap.

Reproduces the ~1.26 order figure from Issue #168. Runs the spherical
isotropic MMS at nc in {10, 20, 40, 80} with three solver paths and
tabulates the orders.

Created by numerics-investigator on 2026-05-09. NOT a test promotion
candidate; this is a reproduction probe.

Run with: pytest scratch/derivations/diagnostics/diag_issue168_01_characterize.py -v -s
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    build_spherical_mms_case,
    build_cylindrical_mms_case,
)
from orpheus.sn import solve_sn_fixed_source


def _l2_1d(phi_num, phi_ref, volumes):
    diff = phi_num - phi_ref
    return float(np.sqrt(np.sum(volumes * diff * diff)))


def _convergence_table(case, n_cells, inner_solver):
    errs = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        kwargs = dict(
            max_inner=200, inner_tol=1e-9,
        )
        if inner_solver is not None:
            kwargs["inner_solver"] = inner_solver
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q, **kwargs
        )
        phi_num = result.scalar_flux[:, 0, 0]
        phi_ref = case.phi_exact(mesh.centers)
        errs.append(_l2_1d(phi_num, phi_ref, mesh.volumes))
    errs = np.asarray(errs)
    orders = np.log2(errs[:-1] / errs[1:])
    return errs, orders


def test_spherical_mms_orders_three_solver_paths():
    """Reproduce the three-solver-path convergence comparison.

    Expected (per Issue #168 narrative):
      - source_iteration (sweep WDD): orders diverge (negative or near-zero)
      - krylov on apply (Round 3 BC fix): ~1.26 — first-order-with-correction
      - true O(h^2) is NOT reached by either path
    """
    case = build_spherical_mms_case(n_ordinates=8)
    n_cells = [10, 20]
    print()
    print(f"=== Spherical MMS (R={case.radius}, sigma_t={case.sigma_t}, "
          f"sigma_s={case.sigma_s}, n_ord={len(case.quadrature.mu_x)}) ===")
    print(f"  ansatz: phi(r) = sin(pi*r/R), vacuum BC at r=R")
    print()

    for inner_solver in ("krylov",):
        label = "default (source_iteration WDD)" if inner_solver is None else "krylov on apply"
        try:
            errs, orders = _convergence_table(case, n_cells, inner_solver)
            print(f"  Solver = {label}")
            for nc, e in zip(n_cells, errs):
                print(f"    nc={nc:3d}  L2 err = {e:.4e}")
            print(f"    Orders (n -> 2n): {orders}")
            print()
        except Exception as exc:
            print(f"  Solver = {label}: FAILED with {type(exc).__name__}: {exc}")
            print()


def test_cylindrical_mms_orders_two_solver_paths():
    """Same but cylindrical Product quadrature."""
    case = build_cylindrical_mms_case()
    n_cells = [10, 20, 40]   # cyl is slower; smaller sweep
    print()
    print(f"=== Cylindrical MMS (R={case.radius}) ===")
    print()

    for inner_solver in (None, "krylov"):
        label = "default (source_iteration WDD)" if inner_solver is None else "krylov on apply"
        try:
            errs, orders = _convergence_table(case, n_cells, inner_solver)
            print(f"  Solver = {label}")
            for nc, e in zip(n_cells, errs):
                print(f"    nc={nc:3d}  L2 err = {e:.4e}")
            print(f"    Orders: {orders}")
            print()
        except Exception as exc:
            print(f"  Solver = {label}: FAILED with {type(exc).__name__}: {exc}")
            print()
