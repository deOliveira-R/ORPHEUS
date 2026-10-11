r"""Foundation regression net for direct math-heart-class construction.

Phase D
-------

The pre-Phase-D regression net verified that polymorphic dispatch on
the ``TransportSolver`` Protocol agreed with the per-method adapter
dispatch. The Protocol was retired in Phase D as part of the
architectural reset (it conflated continuous reference generators
with discrete production solvers, which have functionally different
roles); the agreement contract still matters and lives here.

The 2 tests in this file:

1. ``test_dispatch_fn_slab_matches_adapter`` — running
   ``MomentSpace(geometry=..., materials=...).solve_critical()``
   produces the same critical half-thickness as
   ``FNSlabAdapter().solve(case)``.
2. ``test_dispatch_fn_sphere_matches_adapter`` — sphere counterpart.

Three rows compared ``Billiard`` with the trajectory-resolvent adapters until
step (e2) of the characteristic-reference campaign deleted both (they were
tautological: the adapter called the solver ``Billiard`` dispatched to). The
characteristic adapter has no second construction route to compare with: it
poses the one door, ``characteristic_reference``.
"""
from __future__ import annotations

import warnings

import pytest

from orpheus.derivations.continuous.fn_method.moment_space import MomentSpace
from orpheus.geometry.structured_geometry import StructuredGeometry

from .adapters import (
    FNSlabAdapter,
    FNSphereAdapter,
)
from .cases import (
    BARE_CRITICAL_SLAB_CASES,
    BARE_CRITICAL_SPHERE_CASES,
)


pytestmark = [pytest.mark.foundation]


# ----------------------------------------------------------------------
# Helpers — resolve the StructuredGeometry for a CrossMethodCase
# ----------------------------------------------------------------------


def _structured_geom_for(case) -> StructuredGeometry:
    """Resolve the :class:`StructuredGeometry` for a CrossMethodCase.

    Prefers the inline ``structured_geometry`` override (when set);
    otherwise builds it from ``case.registry_case.to_geometry()``.
    """
    if case.structured_geometry is not None:
        return case.structured_geometry
    return case.registry_case.to_geometry()


# ----------------------------------------------------------------------
# Foundation gate 1 — fn_slab adapter ↔ MomentSpace direct
# ----------------------------------------------------------------------


def test_dispatch_fn_slab_matches_adapter():
    r"""``FNSlabAdapter`` agrees with ``MomentSpace.solve_critical``."""
    case = next(
        c for c in BARE_CRITICAL_SLAB_CASES if "Ua-1-0-SL" in c.case_id
    )

    adapter = FNSlabAdapter(n_modes=8)
    res_adapter = adapter.solve(case)

    geom = _structured_geom_for(case)
    moment = MomentSpace(
        geometry=geom,
        materials=case.registry_case.materials,
        fn_order=adapter.n_modes,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol_protocol = moment.solve_critical()

    assert sol_protocol.parameter_kind == "half_thickness_mfp"
    diff = abs(sol_protocol.parameter_value - res_adapter.value)
    assert diff < 1e-12, (
        f"fn_slab adapter ({res_adapter.value}) vs MomentSpace direct "
        f"({sol_protocol.parameter_value}) drift {diff:.3e} > 1e-12"
    )


# ----------------------------------------------------------------------
# Foundation gate 2 — fn_sphere adapter ↔ MomentSpace direct
# ----------------------------------------------------------------------


def test_dispatch_fn_sphere_matches_adapter():
    r"""``FNSphereAdapter`` agrees with ``MomentSpace.solve_critical``."""
    case = next(
        c for c in BARE_CRITICAL_SPHERE_CASES if "Ua-1-0-SP" in c.case_id
    )

    adapter = FNSphereAdapter(n_modes=8)
    res_adapter = adapter.solve(case)

    geom = _structured_geom_for(case)
    moment = MomentSpace(
        geometry=geom,
        materials=case.registry_case.materials,
        fn_order=adapter.n_modes,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol_protocol = moment.solve_critical()

    assert sol_protocol.parameter_kind == "radius_mfp"
    diff = abs(sol_protocol.parameter_value - res_adapter.value)
    assert diff < 1e-12, (
        f"fn_sphere adapter ({res_adapter.value}) vs MomentSpace "
        f"direct ({sol_protocol.parameter_value}) drift {diff:.3e}"
    )
