r"""Tests for ``SNMesh._resolve_bcs`` Wave-8 realizer wiring (C8.3).

The Wave-8 SNMesh routes BC resolution through
:class:`SNBoundaryRealizer` for Cartesian meshes. Each face attribute
(``bc_xmin`` / ``bc_xmax`` / ``bc_ymin`` / ``bc_ymax`` / ``bc_left`` /
``bc_right``) becomes a :class:`_BoundBoundaryOperator` shim wrapping
the 1-arg realized :class:`LinearOperator`. The shim's
``apply(psi, *_extra, **_kw)`` forwards to the realized op's 1-arg
``apply(psi)`` for backward-compat with the 13 production call sites
that still pass ``(psi, quad)``. Wave 9 migrates those.

Curvilinear meshes (spherical / cylindrical) skip the realizer
because :class:`InflowTraceSpace` is Cartesian-only today (Wave 2
deferred curvilinear). The legacy 2-arg ``BoundaryTraceLaw.apply``
path runs unchanged, preserving bit-identity for the curvilinear
regression snapshots.

V&V tags
--------
``@pytest.mark.l1`` — the wiring assertions are cross-implementation
checks (Wave-5 realizer dispatch + Wave-2 trace-mask construction +
Wave-8 shim composition produce the same observable per-face apply).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.geometry.boundary import (
    ReflectiveBoundary,
    VacuumInflow,
)
from orpheus.geometry.boundary._bound_compat import _BoundBoundaryOperator
from orpheus.sn.geometry import SNMesh
from orpheus.sn.quadrature import (
    GaussLegendre1D,
    LebedevSphere,
    LevelSymmetricSN,
)


pytestmark = pytest.mark.l1


# ─────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────


@pytest.fixture
def quad_2d():
    """LebedevSphere(17) — the canonical 2-D Cartesian quadrature."""
    return LebedevSphere.create(17)


@pytest.fixture
def quad_1d():
    """GaussLegendre1D(8) — 1-D Cartesian / curvilinear quadrature."""
    return GaussLegendre1D.create(8)


# ─────────────────────────────────────────────────────────────────────
# 2-D Cartesian: realizer wiring + §16A.5 inflow-only mask
# ─────────────────────────────────────────────────────────────────────


def test_2d_cartesian_vacuum_xmin_masks_only_inflow(quad_2d):
    """Vacuum on xmin: the realized shim's apply zeros ONLY the inflow
    ordinates per §16A.5 (mu_x > 0 ordinates entering the domain).
    Outflow + tangential ordinates pass through unchanged. This is the
    Wave-8 semantic correction relative to the legacy zeros-all body.
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, 5), edges_y=np.linspace(0, 1, 4),
        mat_map=np.zeros((4, 3), dtype=int),
        bc_xmin=BC("vacuum"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    sn = SNMesh(mesh, quad_2d)
    assert isinstance(sn.bc_xmin, _BoundBoundaryOperator)
    assert sn.bc_xmin == "vacuum"

    rng = np.random.default_rng(42)
    psi = rng.uniform(0.5, 2.0, size=(quad_2d.N, 3, 2))
    out = sn.bc_xmin.apply(psi, quad_2d)
    inflow = np.flatnonzero(quad_2d.mu_x > 1e-12)
    non_inflow = np.setdiff1d(np.arange(quad_2d.N), inflow)
    np.testing.assert_array_equal(out[inflow], 0.0)
    np.testing.assert_array_equal(out[non_inflow], psi[non_inflow])


def test_2d_cartesian_reflective_ymax_returns_permutation(quad_2d):
    """Reflective on ymax: the realized shim's apply returns
    ``psi[ref]`` where ref = quad.reflection_index("y"). Equivalent
    to the legacy SpecularBoundaryOperator(axis="y") output.
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, 5), edges_y=np.linspace(0, 1, 4),
        mat_map=np.zeros((4, 3), dtype=int),
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    sn = SNMesh(mesh, quad_2d)
    assert isinstance(sn.bc_ymax, _BoundBoundaryOperator)
    assert sn.bc_ymax == "reflective"

    rng = np.random.default_rng(1)
    psi = rng.standard_normal(size=(quad_2d.N, 4, 2))
    ref = quad_2d.reflection_index("y")
    expected = psi[ref]
    np.testing.assert_array_equal(sn.bc_ymax.apply(psi, quad_2d), expected)


def test_2d_cartesian_construction_populates_traces(quad_2d):
    """SNMesh on a Cartesian mesh populates :attr:`_inflow_trace` and
    :attr:`_outflow_trace`. Carries the per-face inflow indices used
    by the realizer.
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, 5), edges_y=np.linspace(0, 1, 4),
        mat_map=np.zeros((4, 3), dtype=int),
    )
    sn = SNMesh(mesh, quad_2d)
    assert sn._inflow_trace is not None
    assert sn._outflow_trace is not None
    assert sn._inflow_trace.face_names == ("xmin", "xmax", "ymin", "ymax")
    # xmin: inflow ordinates have mu_x > 0
    np.testing.assert_array_equal(
        sn._inflow_trace.inflow_indices_for_face("xmin"),
        np.flatnonzero(quad_2d.mu_x > 1e-12),
    )


# ─────────────────────────────────────────────────────────────────────
# 1-D Cartesian: realizer wiring through bc_left / bc_right
# ─────────────────────────────────────────────────────────────────────


def test_1d_cartesian_vacuum_right_masks_only_inflow(quad_1d):
    """1-D slab with vacuum on right: the realized shim zeros only the
    inflow rows (mu_x < 0 at the right face).
    """
    mesh = Mesh1D(
        edges=np.linspace(0, 2, 9),
        mat_ids=np.zeros(8, dtype=int),
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    sn = SNMesh(mesh, quad_1d)
    assert isinstance(sn.bc_right, _BoundBoundaryOperator)
    assert sn.bc_right == "vacuum"

    psi = np.arange(quad_1d.N * 2, dtype=float).reshape(quad_1d.N, 2)
    out = sn.bc_right.apply(psi, quad_1d)
    inflow = np.flatnonzero(quad_1d.mu_x < -1e-12)
    non_inflow = np.setdiff1d(np.arange(quad_1d.N), inflow)
    np.testing.assert_array_equal(out[inflow], 0.0)
    np.testing.assert_array_equal(out[non_inflow], psi[non_inflow])


def test_1d_cartesian_y_face_placeholders_unwrapped(quad_1d):
    """1-D meshes carry bare :class:`ReflectiveBoundary` placeholders
    on the degenerate y faces — no realizer routing because the 1-D
    trace space cannot service ``inflow_indices_for_face('ymin')``.
    These placeholders are never consumed by 1-D sweeps; the
    invariant pins they remain bare BoundaryTraceLaw instances.
    """
    mesh = Mesh1D(
        edges=np.linspace(0, 1, 5),
        mat_ids=np.zeros(4, dtype=int),
    )
    sn = SNMesh(mesh, quad_1d)
    assert isinstance(sn.bc_ymin, ReflectiveBoundary)
    assert isinstance(sn.bc_ymax, ReflectiveBoundary)
    assert sn.bc_ymin.axis == "y"
    assert sn.bc_ymax.axis == "y"
    # And they are NOT shim-wrapped.
    assert not isinstance(sn.bc_ymin, _BoundBoundaryOperator)


# ─────────────────────────────────────────────────────────────────────
# 1-D curvilinear: legacy path (no realizer)
# ─────────────────────────────────────────────────────────────────────


def test_1d_spherical_vacuum_uses_legacy_path(quad_1d):
    """Spherical vacuum: ``InflowTraceSpace`` raises NotImplementedError
    on curvilinear, so ``_resolve_one`` falls back to the bare
    :class:`VacuumInflow` legacy law. ``apply(psi, quad)`` returns the
    legacy ``np.zeros_like`` body.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.1, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_right=BC("vacuum"),
    )
    sn = SNMesh(mesh, quad_1d)
    assert isinstance(sn.bc_right, VacuumInflow)
    # NOT shim-wrapped.
    assert not isinstance(sn.bc_right, _BoundBoundaryOperator)
    # No trace spaces are built for curvilinear.
    assert sn._inflow_trace is None
    assert sn._outflow_trace is None

    psi = np.random.default_rng(2).standard_normal(size=(quad_1d.N, 2))
    out = sn.bc_right.apply(psi, quad_1d)
    np.testing.assert_array_equal(out, np.zeros_like(psi))


def test_1d_cylindrical_reflective_uses_legacy_path():
    """Cylindrical reflective: same as spherical — bare
    :class:`ReflectiveBoundary` instance, no shim, no realizer.
    Cylindrical needs a level-structured quadrature (LevelSymmetricSN
    or ProductQuadrature); the legacy path is used regardless of
    quadrature family.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.1, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"), bc_right=BC("reflective"),
    )
    quad = LevelSymmetricSN.create(sn_order=4)
    sn = SNMesh(mesh, quad)
    assert isinstance(sn.bc_left, ReflectiveBoundary)
    assert sn.bc_left.axis == "x"
    assert not isinstance(sn.bc_left, _BoundBoundaryOperator)
    assert sn._inflow_trace is None


# ─────────────────────────────────────────────────────────────────────
# Registry contract
# ─────────────────────────────────────────────────────────────────────


def test_registry_contains_only_vacuum_and_reflective():
    """The SN registry pins exactly the kinds the legacy SNMesh
    accepted today. Adding ``white`` / ``periodic`` / ``albedo`` etc.
    requires sweep-side wiring out of Wave 8 scope.
    """
    assert set(SNMesh.BOUNDARY_OPERATOR_REGISTRY) == {"vacuum", "reflective"}


def test_unknown_bc_kind_raises_valueerror(quad_1d):
    """Unsupported BC kind raises ``ValueError`` listing the supported
    set. Pinned for the BC-resolution diagnostic contract.
    """
    mesh = Mesh1D(
        edges=np.linspace(0, 1, 5), mat_ids=np.zeros(4, dtype=int),
        bc_left=BC("periodic"),
    )
    with pytest.raises(ValueError, match="'reflective'.*'vacuum'"):
        SNMesh(mesh, quad_1d)
