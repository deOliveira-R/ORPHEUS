r"""Foundation regression net for the TransportSolver Protocol.

This file is the **safety net** required by the test-architect spec
§ "File 4". It pins that polymorphic dispatch on the
:class:`~orpheus.derivations.common.solver_protocol.TransportSolver`
Protocol agrees with the per-method adapter dispatch in
:mod:`tests.cross_method.adapters` for the canonical case set.

Why a regression net
--------------------

The cross-method gates in :mod:`tests.cross_method.test_eigenvalue`
have exercised the per-method adapters (``FNSlabAdapter`` /
``TrajectoryResolventSlabAdapter`` / ``FNSphereAdapter`` /
``TrajectoryResolventSphereAdapter`` / ``TrajectoryResolventSphereClosedAdapter``)
for years. The adapters do unit conversion (mfp ↔ cm, half-thickness
↔ full slab) and parameter selection (n_modes for fn_method,
n_r/n_mu for trajectory_resolvent). They are the load-bearing
extraction layer and stay; what's new is the Protocol that lets a
test body dispatch on either adapter or directly on ``Billiard`` /
``MomentSpace``.

The 5 tests in this file:

1. ``test_polymorphic_dispatch_fn_slab_matches_adapter`` — running
   ``MomentSpace.from_problem(...).solve_critical()`` produces the
   same critical half-thickness as ``FNSlabAdapter().solve(case)``.
2. ``test_polymorphic_dispatch_fn_sphere_matches_adapter`` — sphere
   counterpart.
3. ``test_polymorphic_dispatch_trajectory_resolvent_slab_matches_adapter``
   — Billiard slab vs ``TrajectoryResolventSlabAdapter``.
4. ``test_polymorphic_dispatch_trajectory_resolvent_sphere_matches_adapter``
   — Billiard sphere vs ``TrajectoryResolventSphereAdapter``.
5. ``test_polymorphic_dispatch_closed_sphere_matches_adapter`` —
   Billiard closed-sphere vs ``TrajectoryResolventSphereClosedAdapter``.

Foundation-tier
---------------

Module-level ``pytestmark = [pytest.mark.foundation]`` (no
``verifies(...)``); these are software invariants pinning that the
adapter layer and the Protocol layer agree.

Reduced quadrature
------------------

To keep the foundation suite under 5 s aggregate per
:doc:`/skills/algebra-of-record` § "Cost discipline", the tests
configure each solver at SMALL grid sizes (n_r=12, n_mu=12,
n_traj_quad=24 for trajectory_resolvent; fn_order=8 for fn_method).
The agreement check uses bit-equality where the math is determinate
and a loose 1e-3 numerical floor where iteration / quadrature
precision is at play.

References
----------

* :doc:`/theory/transport_solver_protocol` § "Protocol vs adapter
  layer" — why both layers coexist.
"""
from __future__ import annotations

import warnings

import pytest

from orpheus.derivations.common.solver_protocol import TransportSolver
from orpheus.derivations.continuous.fn_method.moment_space import MomentSpace
from orpheus.derivations.continuous.trajectory_resolvent import Billiard

from .adapters import (
    FNSlabAdapter,
    FNSphereAdapter,
    TrajectoryResolventSlabAdapter,
    TrajectoryResolventSphereAdapter,
    TrajectoryResolventSphereClosedAdapter,
    _geometry_spec_for,
)
from .cases import (
    BARE_CRITICAL_SLAB_CASES,
    BARE_CRITICAL_SPHERE_CASES,
    CLOSED_SPHERE_KINF_CASES,
)


pytestmark = [pytest.mark.foundation]


# ----------------------------------------------------------------------
# Foundation gate 1 — Protocol conformance of the math-heart classes
# ----------------------------------------------------------------------


def test_polymorphic_dispatch_protocol_conformance():
    r"""Every Protocol-conforming class actually conforms.

    Sanity gate: catches a regression where someone removes
    ``method_name`` or accidentally types ``materials`` / ``geometry_spec``
    incorrectly. ``Billiard`` and ``MomentSpace`` MUST satisfy the
    Protocol; cross-method tests rely on it.
    """
    # Pick the canonical sphere case (Sood Ua-1-0-SP, c=1.30).
    case = next(
        c for c in BARE_CRITICAL_SPHERE_CASES if "Ua-1-0-SP" in c.case_id
    )
    spec = _geometry_spec_for(case)

    # Construct a Billiard via the production-protocol signature.
    billiard = Billiard.from_problem(
        materials=case.registry_case.materials,
        geometry_spec=spec,
        alpha=spec.bc_right.to_alpha(),
    )
    assert isinstance(billiard, TransportSolver), (
        "Billiard does not satisfy TransportSolver Protocol"
    )
    assert billiard.method_name == "trajectory_resolvent"

    # Construct a MomentSpace via the production-protocol signature.
    moment = MomentSpace.from_problem(
        materials=case.registry_case.materials,
        geometry=spec,
        fn_order=8,
    )
    assert isinstance(moment, TransportSolver), (
        "MomentSpace does not satisfy TransportSolver Protocol"
    )
    assert moment.method_name == "fn_method"


# ----------------------------------------------------------------------
# Foundation gate 2 — fn_slab adapter ↔ MomentSpace direct
# ----------------------------------------------------------------------


def test_polymorphic_dispatch_fn_slab_matches_adapter():
    r"""``FNSlabAdapter`` agrees with ``MomentSpace.solve_critical``.

    Run the Sood Ua-1-0-SL c=1.30 case through both the adapter and
    the direct ``MomentSpace`` Protocol path. Both routes converge on
    the SAME function-level call (``solve_fn_slab_bare_critical``),
    so the agreement is bit-equal modulo the route's own floating-
    point operations on c (inputs are identical).
    """
    case = next(
        c for c in BARE_CRITICAL_SLAB_CASES if "Ua-1-0-SL" in c.case_id
    )

    adapter = FNSlabAdapter(n_modes=8)
    res_adapter = adapter.solve(case)

    # Direct path: build MomentSpace via production-protocol inputs.
    moment = MomentSpace.from_problem(
        materials=case.registry_case.materials,
        geometry=_geometry_spec_for(case),
        fn_order=adapter.n_modes,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sol_protocol = moment.solve_critical()

    # The adapter returns ``a_critical_mfp``; the Protocol surface
    # returns the same value via ``parameter_value`` with
    # ``parameter_kind == "half_thickness_mfp"``.
    assert sol_protocol.parameter_kind == "half_thickness_mfp"
    diff = abs(sol_protocol.parameter_value - res_adapter.value)
    # Both routes go through the same function-level API, so agreement
    # should be bit-equal — loose 1e-12 floor leaves room for any
    # incidental ULP drift.
    assert diff < 1e-12, (
        f"fn_slab adapter ({res_adapter.value}) vs MomentSpace direct "
        f"({sol_protocol.parameter_value}) drift {diff:.3e} > 1e-12"
    )


# ----------------------------------------------------------------------
# Foundation gate 3 — fn_sphere adapter ↔ MomentSpace direct
# ----------------------------------------------------------------------


def test_polymorphic_dispatch_fn_sphere_matches_adapter():
    r"""``FNSphereAdapter`` agrees with ``MomentSpace.solve_critical``."""
    case = next(
        c for c in BARE_CRITICAL_SPHERE_CASES if "Ua-1-0-SP" in c.case_id
    )

    adapter = FNSphereAdapter(n_modes=8)
    res_adapter = adapter.solve(case)

    moment = MomentSpace.from_problem(
        materials=case.registry_case.materials,
        geometry=_geometry_spec_for(case),
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


# ----------------------------------------------------------------------
# Foundation gate 4 — trajectory_resolvent_slab ↔ Billiard direct
# ----------------------------------------------------------------------


def test_polymorphic_dispatch_trajectory_resolvent_slab_matches_adapter():
    r"""``TrajectoryResolventSlabAdapter`` agrees with ``Billiard.solve_critical``.

    Both paths land at the same ``solve_greens_function_slab`` call
    with the same arguments; ``k_eff`` agreement should be bit-equal.

    Quadrature is reduced to (n_x=8, n_mu=8, n_traj_quad=16) to keep
    the foundation suite fast (default would be n_x=48, n_mu=128,
    n_traj_quad=96 ≈ 30 s). The agreement check is structural; we
    don't assert against truth — that's the L1 truth gate's job in
    ``test_eigenvalue.py``.
    """
    case = next(
        c for c in BARE_CRITICAL_SLAB_CASES if "Ua-1-0-SL" in c.case_id
    )

    # Reduced-quadrature adapter for foundation cost discipline.
    adapter = TrajectoryResolventSlabAdapter(
        n_x=8, n_mu=8, n_traj_quad=16, max_iter=100, tol=1e-7,
    )
    res_adapter = adapter.solve(case)

    # Direct Billiard via the production-protocol signature; thread
    # through the SAME quadrature so the two routes agree bit-equal.
    spec = _geometry_spec_for(case)
    alpha = spec.bc_right.to_alpha()

    billiard = Billiard.from_problem(
        materials=case.registry_case.materials,
        geometry_spec=spec,
        alpha=alpha,
        quadrature={"n_x": 8, "n_mu": 8, "n_traj_quad": 16},
    )
    sol_protocol = billiard.solve_critical(max_iter=100, tol=1e-7)

    assert sol_protocol.eigenvalue_kind == "k_eff"
    diff = abs(sol_protocol.eigenvalue - res_adapter.value)
    assert diff < 1e-14, (
        f"trajectory_resolvent_slab adapter ({res_adapter.value}) vs "
        f"Billiard direct ({sol_protocol.eigenvalue}) drift {diff:.3e}"
    )


# ----------------------------------------------------------------------
# Foundation gate 5 — trajectory_resolvent_sphere ↔ Billiard direct
# ----------------------------------------------------------------------


def test_polymorphic_dispatch_trajectory_resolvent_sphere_matches_adapter():
    r"""``TrajectoryResolventSphereAdapter`` agrees with ``Billiard``."""
    case = next(
        c for c in BARE_CRITICAL_SPHERE_CASES if "Ua-1-0-SP" in c.case_id
    )

    adapter = TrajectoryResolventSphereAdapter(
        n_r=12, n_mu=12, n_traj_quad=24, max_iter=100, tol=1e-8,
    )
    res_adapter = adapter.solve(case)

    spec = _geometry_spec_for(case)
    alpha = spec.bc_right.to_alpha()

    billiard = Billiard.from_problem(
        materials=case.registry_case.materials,
        geometry_spec=spec,
        alpha=alpha,
        quadrature={"n_r": 12, "n_mu": 12, "n_traj_quad": 24},
    )
    sol_protocol = billiard.solve_critical(max_iter=100, tol=1e-8)

    assert sol_protocol.eigenvalue_kind == "k_eff"
    diff = abs(sol_protocol.eigenvalue - res_adapter.value)
    assert diff < 1e-14, (
        f"trajectory_resolvent_sphere adapter ({res_adapter.value}) "
        f"vs Billiard direct ({sol_protocol.eigenvalue}) drift "
        f"{diff:.3e}"
    )


# ----------------------------------------------------------------------
# Foundation gate 6 — closed_sphere_kinf ↔ Billiard direct
# ----------------------------------------------------------------------


def test_polymorphic_dispatch_closed_sphere_matches_adapter():
    r"""``TrajectoryResolventSphereClosedAdapter`` agrees with ``Billiard``.

    The closed-sphere case (α=1) carries inline materials and
    geometry_spec (no registry case). Both the adapter and the
    direct Billiard path consume the same ``case.materials``
    (``dict[int, Mixture]``) and ``case.geometry_spec`` pair via
    the production-protocol signature.
    """
    case = next(iter(CLOSED_SPHERE_KINF_CASES))

    adapter = TrajectoryResolventSphereClosedAdapter(
        n_r=12, n_mu=12, n_traj_quad=24, max_iter=50, tol=1e-12,
    )
    res_adapter = adapter.solve(case)

    # Closed-sphere uses inline materials (no registry case).
    spec = _geometry_spec_for(case)
    alpha = spec.bc_right.to_alpha()

    billiard = Billiard.from_problem(
        materials=case.materials,
        geometry_spec=spec,
        alpha=alpha,
        quadrature={"n_r": 12, "n_mu": 12, "n_traj_quad": 24},
    )
    sol_protocol = billiard.solve_critical(max_iter=50, tol=1e-12)

    diff = abs(sol_protocol.eigenvalue - res_adapter.value)
    assert diff < 1e-14, (
        f"closed_sphere adapter ({res_adapter.value}) vs Billiard "
        f"direct ({sol_protocol.eigenvalue}) drift {diff:.3e}"
    )
