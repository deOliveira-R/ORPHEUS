r"""Foundation regression net for direct math-heart-class construction.

Phase D
-------

The pre-Phase-D regression net verified that polymorphic dispatch on
the ``TransportSolver`` Protocol agreed with the per-method adapter
dispatch. The Protocol was retired in Phase D as part of the
architectural reset (it conflated continuous reference generators
with discrete production solvers, which have functionally different
roles); the agreement contract still matters and lives here.

The 5 tests in this file:

1. ``test_dispatch_fn_slab_matches_adapter`` — running
   ``MomentSpace(geometry=..., materials=...).solve_critical()``
   produces the same critical half-thickness as
   ``FNSlabAdapter().solve(case)``.
2. ``test_dispatch_fn_sphere_matches_adapter`` — sphere counterpart.
3. ``test_dispatch_trajectory_resolvent_slab_matches_adapter`` —
   ``Billiard`` slab vs ``TrajectoryResolventSlabAdapter``.
4. ``test_dispatch_trajectory_resolvent_sphere_matches_adapter`` —
   ``Billiard`` sphere vs ``TrajectoryResolventSphereAdapter``.
5. ``test_dispatch_closed_sphere_matches_adapter`` — ``Billiard``
   closed-sphere vs ``TrajectoryResolventSphereClosedAdapter``.
"""
from __future__ import annotations

import warnings

import pytest

from orpheus.derivations.continuous.fn_method.moment_space import MomentSpace
from orpheus.derivations.continuous.trajectory_resolvent import Billiard
from orpheus.geometry.mesh import BC
from orpheus.geometry.structured_geometry import (
    Region,
    StructuredGeometry,
)

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
# Helpers — bridge legacy GeometrySpec cases to StructuredGeometry
# ----------------------------------------------------------------------


def _structured_geom_from_case(spec) -> StructuredGeometry:
    """Build a StructuredGeometry equivalent to a registry case's GeometrySpec.

    The registry / cross_method cases still carry GeometrySpec
    (Phase F territory). Translate to the new geometry layer so the
    direct math-heart-class construction path can consume them.
    """
    tag_map = {"slab": "SLB", "sphere": "SPH", "cylinder": "CYL"}
    tag = tag_map[spec.geometry]
    if tag == "SLB":
        # GeometrySpec.domain_extent_cm is the FULL slab width
        # (2 * critical_dimension_cm), which is what
        # StructuredGeometry stores end-to-end as well.
        return StructuredGeometry(
            geometry="SLB",
            regions=(
                Region(mat_id=0, outer_thickness_cm=spec.domain_extent_cm),
            ),
            bcs=(spec.bc_left, spec.bc_right),
        )
    # CYL / SPH: single-region radius, single outer BC.
    return StructuredGeometry(
        geometry=tag,
        regions=(
            Region(mat_id=0, outer_thickness_cm=spec.domain_extent_cm),
        ),
        bcs=(spec.bc_right,),
    )


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

    geom = _structured_geom_from_case(_geometry_spec_for(case))
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

    geom = _structured_geom_from_case(_geometry_spec_for(case))
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


# ----------------------------------------------------------------------
# Foundation gate 3 — trajectory_resolvent_slab ↔ Billiard direct
# ----------------------------------------------------------------------


def test_dispatch_trajectory_resolvent_slab_matches_adapter():
    r"""``TrajectoryResolventSlabAdapter`` agrees with ``Billiard.solve_critical``."""
    case = next(
        c for c in BARE_CRITICAL_SLAB_CASES if "Ua-1-0-SL" in c.case_id
    )

    adapter = TrajectoryResolventSlabAdapter(
        n_x=8, n_mu=8, n_traj_quad=16, max_iter=100, tol=1e-7,
    )
    res_adapter = adapter.solve(case)

    spec = _geometry_spec_for(case)
    alpha = spec.bc_right.to_alpha()
    geom = _structured_geom_from_case(spec)

    billiard = Billiard(
        geometry=geom,
        materials=case.registry_case.materials,
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
# Foundation gate 4 — trajectory_resolvent_sphere ↔ Billiard direct
# ----------------------------------------------------------------------


def test_dispatch_trajectory_resolvent_sphere_matches_adapter():
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
    geom = _structured_geom_from_case(spec)

    billiard = Billiard(
        geometry=geom,
        materials=case.registry_case.materials,
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
# Foundation gate 5 — closed_sphere_kinf ↔ Billiard direct
# ----------------------------------------------------------------------


def test_dispatch_closed_sphere_matches_adapter():
    r"""``TrajectoryResolventSphereClosedAdapter`` agrees with ``Billiard``."""
    case = next(iter(CLOSED_SPHERE_KINF_CASES))

    adapter = TrajectoryResolventSphereClosedAdapter(
        n_r=12, n_mu=12, n_traj_quad=24, max_iter=50, tol=1e-12,
    )
    res_adapter = adapter.solve(case)

    spec = _geometry_spec_for(case)
    alpha = spec.bc_right.to_alpha()
    geom = StructuredGeometry(
        geometry="SPH",
        regions=(
            Region(mat_id=0, outer_thickness_cm=spec.domain_extent_cm),
        ),
        # Closed sphere → outer is reflective (the inline cases set
        # bc_right=BC.reflective; the adapter routes alpha=1 from
        # bc_right.to_alpha()). StructuredGeometry single-endpoint
        # carries only the outer BC.
        bcs=(BC.reflective,),
    )

    billiard = Billiard(
        geometry=geom,
        materials=case.materials,
        alpha=alpha,
        quadrature={"n_r": 12, "n_mu": 12, "n_traj_quad": 24},
    )
    sol_protocol = billiard.solve_critical(max_iter=50, tol=1e-12)

    diff = abs(sol_protocol.eigenvalue - res_adapter.value)
    assert diff < 1e-14, (
        f"closed_sphere adapter ({res_adapter.value}) vs Billiard "
        f"direct ({sol_protocol.eigenvalue}) drift {diff:.3e}"
    )
