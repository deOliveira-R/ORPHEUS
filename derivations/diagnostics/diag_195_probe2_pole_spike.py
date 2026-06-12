"""Diagnostic: Issue #195 Probe 2 — error localisation + pole-spike persistence.

Created by numerics-investigator on 2026-06-12.

VERDICT EVIDENCE for root cause (2): the curvilinear-isotropic MMS L2 error
PLATEAUS (does not converge) because a pole-cell flux spike GROWS under
refinement.  This script pins the spike:

  Sphere err[0] (pole cell):  nx 40→640:  0.160 → 0.264 → 0.346 → 0.402 → 0.437
  phi_num[0] grows 0.20→0.44 while ref[0] (correctly) shrinks 0.039→0.0024.
  L2 plateaus ~0.0413 because V_pole ~ r²·dr ~ h³ shrinks as the spike grows.

A bulk-distributed C·h² profile would support root cause (1); this
boundary(pole)-concentrated O(1) spike that GROWS with refinement is the
fingerprint of a discretisation defect at the pole cell (root cause 2).

This is a LIVE ERR-026 manifestation on a NON-FLAT radial profile — the
Phase D Carlson seed closed the per-ordinate FLAT-flux identity, but the
sin(πr/R) profile exposes the pole-cell assembly residual that does not
cancel in the scalar moment.

If promoted: a regression gate that the pole-cell error must DECREASE under
refinement (it currently grows) — catches the wrong-fixed-point class.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.continuous.mms.sn import (
    build_spherical_mms_case,
    build_cylindrical_mms_case,
)
from orpheus.sn import solve_sn_fixed_source


def _solve_profile(case, nc):
    mesh = case.build_mesh(nc)
    Q = case.external_source(mesh)
    r = solve_sn_fixed_source(
        case.materials, mesh, case.quadrature, Q,
        max_inner=2000, inner_tol=1e-13,
        inner_solver="source_iteration",  # == Krylov fixed point, 170× faster
    )
    phi = r.scalar_flux.values[0, :]
    ref = case.phi_exact(mesh.centers)
    return mesh, phi, ref


def _pole_error_ladder(case, ncs):
    pole_errs = []
    L2s = []
    for nc in ncs:
        mesh, phi, ref = _solve_profile(case, nc)
        err = phi - ref
        pole_errs.append(abs(err[0]))
        L2s.append(float(np.sqrt(np.sum(mesh.volumes * err * err))))
    return np.asarray(pole_errs), np.asarray(L2s)


@pytest.mark.parametrize("builder,name", [
    (build_spherical_mms_case, "sphere"),
    (build_cylindrical_mms_case, "cylinder"),
])
def test_probe2_pole_spike_grows_under_refinement(builder, name):
    """Demonstrate the pole-cell error GROWS (the root-cause-2 fingerprint).

    A correct discretisation would have pole_err DECREASE with refinement.
    This test FAILS on current main (pole error grows) — the desired
    behaviour after the pole-closure fix is pole_err[-1] < pole_err[0].
    """
    case = builder()
    ncs = [40, 80, 160, 320]
    pole_errs, L2s = _pole_error_ladder(case, ncs)
    print(f"\n{name}: pole_err ladder = {pole_errs}")
    print(f"{name}: L2 ladder       = {L2s}")
    # The DIAGNOSTIC assertion (currently FAILS on main, documents the bug):
    assert pole_errs[-1] < pole_errs[0], (
        f"{name}: pole-cell error GROWS under refinement "
        f"({pole_errs[0]:.4e} → {pole_errs[-1]:.4e}) — ERR-026 pole-closure "
        f"defect on non-flat profile (root cause 2). The L2 plateau "
        f"({L2s[0]:.4e} → {L2s[-1]:.4e}) is a consequence, not a transient."
    )


if __name__ == "__main__":
    for builder, name in [
        (build_spherical_mms_case, "SPHERE"),
        (build_cylindrical_mms_case, "CYLINDER"),
    ]:
        case = builder()
        print(f"\n========== {name} ==========")
        ncs = [40, 80, 160, 320, 640]
        print(f"{'nx':>5} {'r_0':>8} {'phi[0]':>10} {'ref[0]':>10} "
              f"{'pole_err':>11} {'L2':>11} {'L2 order':>9}")
        prevL2 = None
        for nc in ncs:
            mesh, phi, ref = _solve_profile(case, nc)
            err = phi - ref
            L2 = float(np.sqrt(np.sum(mesh.volumes * err * err)))
            o = "" if prevL2 is None else f"{np.log2(prevL2/L2):9.3f}"
            print(f"{nc:>5} {mesh.centers[0]:8.4f} {phi[0]:10.5f} {ref[0]:10.5f} "
                  f"{abs(err[0]):11.4e} {L2:11.4e} {o:>9}")
            prevL2 = L2
