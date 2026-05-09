"""Diagnostic: confirm the monkey-patch is being consumed.

The Option A patch in diag_issue168_02 sets opmod.transport_operator_matvec_spherical
to a different function. SNStreamingOperator.apply() calls
``transport_operator_matvec_spherical(...)`` from within the operator module
itself, which resolves it via the module's global namespace. The patch
should therefore be picked up. Verify directly.
"""
from __future__ import annotations

import numpy as np
import orpheus.sn.operator as opmod
from orpheus.derivations.continuous.mms.sn import build_spherical_mms_case
from orpheus.sn import solve_sn_fixed_source


def test_patch_is_consumed():
    """Replace the matvec with a marker; ensure SNStreamingOperator.apply
    routes through the patched function."""
    case = build_spherical_mms_case(n_ordinates=4)
    mesh = case.build_mesh(8)
    Q = case.external_source(mesh)

    call_count = [0]
    original = opmod.transport_operator_matvec_spherical

    def tracking_wrapper(*args, **kwargs):
        call_count[0] += 1
        return original(*args, **kwargs)

    opmod.transport_operator_matvec_spherical = tracking_wrapper
    try:
        result = solve_sn_fixed_source(
            case.materials, mesh, case.quadrature, Q,
            inner_solver="krylov",
            max_inner=20, inner_tol=1e-6,
        )
    finally:
        opmod.transport_operator_matvec_spherical = original

    print(f"\nCall count after solve: {call_count[0]}")
    assert call_count[0] > 0, "patched matvec was never called"
