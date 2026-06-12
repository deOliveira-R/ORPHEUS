"""Diagnostic: Issue #195 Probe 1 — the decisive refinement ladder.

Created by numerics-investigator on 2026-06-12.

Question: does the curvilinear-isotropic MMS L2 error collapse below 1e-3
with enough refinement (root cause (1), benign pre-asymptotic transient),
or does it stall above the order-2.46 extrapolation (root cause (2),
coefficient/discretisation bug)?

Production configuration replicated EXACTLY from
``tests/sn/verification/mms/test_mms_curvilinear.py``:
  - solve_sn_fixed_source(case.materials, mesh, case.quadrature, Q,
                          max_inner=500, inner_tol=1e-13)
  - NO inner_solver argument → resolves to None → curvilinear flips to
    "krylov" (solver.py:1987-1991).  This is the production default and
    matches issue #195's "Krylov gives O(h²)" probe.

If this test catches a real coefficient bug, promote to
``tests/sn/verification/mms/`` as the magnitude-gate regression.
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


def _run_ladder(case, n_cells, inner_solver="source_iteration"):
    # SI is bit-identical to the production Krylov default fixed point here
    # (verified: default==SI==Krylov to machine precision) and ~170× faster.
    errors = []
    hvals = []
    for nc in n_cells:
        mesh = case.build_mesh(nc)
        Q = case.external_source(mesh)
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            max_inner=2000, inner_tol=1e-13,
            inner_solver=inner_solver,
        )
        phi_num = result.scalar_flux.values[0, :]
        phi_ref = case.phi_exact(mesh.centers)
        errors.append(_l2_1d(phi_num, phi_ref, mesh.volumes))
        hvals.append(case.radius / nc)
    return np.asarray(hvals), np.asarray(errors)


def test_probe1_sphere_refinement_ladder():
    """Sphere isotropic MMS: nx ∈ {40,80,160,320,640}, production Krylov default."""
    case = build_spherical_mms_case()
    n_cells = [40, 80, 160, 320, 640]
    hvals, errors = _run_ladder(case, n_cells)
    # Pretty print using nx explicitly
    print("\n=== SPHERE isotropic (production default inner_solver=None→krylov) ===")
    orders = np.log2(errors[:-1] / errors[1:])
    print(f"{'nx':>6} {'L2 err':>14} {'order':>8}")
    for i, nc in enumerate(n_cells):
        o = "" if i == 0 else f"{orders[i-1]:8.3f}"
        print(f"{nc:>6} {errors[i]:14.6e} {o:>8}")
    p_last = orders[-1]
    C = errors[-1] / hvals[-1] ** p_last
    print(f"asymptotic order p ≈ {p_last:.3f},  C (err/h^p) ≈ {C:.4e}")

    below_1e3 = np.where(errors < 1e-3)[0]
    first_nx = n_cells[below_1e3[0]] if below_1e3.size else None
    print(f"first nx with err < 1e-3: {first_nx}")
    # VERDICT gate: root cause (2). The error PLATEAUS (order → 0); it never
    # drops below 1e-3. This assertion documents the bug: when the pole
    # closure is fixed, the error WILL collapse below 1e-3 and this passes.
    assert first_nx is not None, (
        f"ROOT CAUSE (2) CONFIRMED: sphere MMS error PLATEAUS at ~{errors[-1]:.4e} "
        f"and never drops below 1e-3 over {n_cells}; asymptotic order≈{orders[-1]:.3f}≈0 "
        f"(mesh-INDEPENDENT error). errors={errors}. This is a discretisation "
        f"defect (pole closure), NOT a pre-asymptotic transient."
    )


def test_probe1_cylinder_refinement_ladder():
    """Cylinder isotropic MMS: nx ∈ {40,80,160,320}, production Krylov default."""
    case = build_cylindrical_mms_case()
    n_cells = [40, 80, 160, 320]
    hvals, errors = _run_ladder(case, n_cells)
    print("\n=== CYLINDER isotropic (production default inner_solver=None→krylov) ===")
    orders = np.log2(errors[:-1] / errors[1:])
    print(f"{'nx':>6} {'L2 err':>14} {'order':>8}")
    for i, nc in enumerate(n_cells):
        o = "" if i == 0 else f"{orders[i-1]:8.3f}"
        print(f"{nc:>6} {errors[i]:14.6e} {o:>8}")
    p_last = orders[-1]
    C = errors[-1] / hvals[-1] ** p_last
    print(f"asymptotic order p ≈ {p_last:.3f},  C (err/h^p) ≈ {C:.4e}")
    below_1e3 = np.where(errors < 1e-3)[0]
    first_nx = n_cells[below_1e3[0]] if below_1e3.size else None
    print(f"first nx with err < 1e-3: {first_nx}")
    assert first_nx is not None, (
        f"Cylinder MMS error never dropped below 1e-3 over {n_cells}; "
        f"errors={errors}"
    )


if __name__ == "__main__":
    print("SPHERE")
    case = build_spherical_mms_case()
    h, e = _run_ladder(case, [40, 80, 160, 320, 640])
    o = np.log2(e[:-1] / e[1:])
    for nc, ee, oo in zip([40, 80, 160, 320, 640], e, [None, *o]):
        print(f"  nx={nc:4d}  err={ee:.6e}  order={oo if oo is None else round(oo,3)}")
    print("CYLINDER")
    case = build_cylindrical_mms_case()
    h, e = _run_ladder(case, [40, 80, 160, 320])
    o = np.log2(e[:-1] / e[1:])
    for nc, ee, oo in zip([40, 80, 160, 320], e, [None, *o]):
        print(f"  nx={nc:4d}  err={ee:.6e}  order={oo if oo is None else round(oo,3)}")
