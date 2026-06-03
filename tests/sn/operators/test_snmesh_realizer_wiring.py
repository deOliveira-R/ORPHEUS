r"""Tests for ``SNMesh._resolve_bcs`` Wave-8 + C188.3 realizer wiring.

The Wave-8 SNMesh routes BC resolution through
:class:`SNBoundaryRealizer` for every supported mesh. Each face
attribute (``bc_xmin`` / ``bc_xmax`` / ``bc_ymin`` / ``bc_ymax`` /
``bc_left`` / ``bc_right``) becomes a :class:`_BoundBoundaryOperator`
shim wrapping the 1-arg realized :class:`LinearOperator`. The shim's
``apply(psi, *_extra, **_kw)`` forwards to the realized op's 1-arg
``apply(psi)`` for backward-compat with the 13 production call sites
that still pass ``(psi, quad)``. Wave 9 migrates those.

Issue #188 / C188.3 (this revision): the curvilinear bypass branch
in ``_resolve_one`` is gone. With the unified
:class:`~orpheus.numerics.spaces.trace_space.TraceSpace`'s curvilinear
support, 1-D spherical and 1-D cylindrical meshes route through the
SAME realizer-then-shim path as Cartesian meshes — but a solid
sphere / cylinder has only the outer (``xmax``) boundary face; the
pole r=0 is the angular closure's regularity condition, so
``bc_left`` / ``bc_xmin`` are ``None``. The 1-D y-face placeholders
use a minimal method space (no trace; the realizer's
:class:`ReflectiveBoundary` branch does not consume inflow_indices).
For GL1D the y-reflection permutation is the identity, so the
realized op is a no-op.

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
from orpheus.geometry.boundary._bound_compat import _BoundBoundaryOperator
from orpheus.numerics.operator import (
    IdentityOperator,
    IncomingOrdinateMaskTensor,
    LinearOperator,
    PermutationOperator,
    TensorProductOperator,
)
from orpheus.sn.geometry import SNMesh
from orpheus.numerics.quadrature import Quadrature
from tests.sn._test_helpers import placeholder_materials


pytestmark = pytest.mark.l1


def _angular_factor(op: LinearOperator) -> LinearOperator:
    """Return the angular factor of a realized boundary operator.

    Wave T step T.1 lifted every realized boundary operator into a
    ``TensorProductOperator`` of the form ``<angular-op> ⊗ Identity``
    (the angular law on the ordinate axis, identity on the spatial
    axis). Before T.1 the realized ``inner`` WAS the bare angular op.
    This accessor returns the angular factor (``ops[0]``) when the
    inner is a tensor product, else the bare op — so the structural
    assertions below survive the T.1 lift without being blind to it.
    """
    if isinstance(op, TensorProductOperator):
        # ``<angular-op> ⊗ Identity``: ops[1] MUST be the spatial
        # identity, ops[0] the angular law the test inspects.
        assert len(op.ops) == 2, f"expected a 2-factor TP, got {len(op.ops)}"
        assert isinstance(op.ops[1], IdentityOperator), (
            f"TP spatial factor is {type(op.ops[1]).__name__}, not Identity"
        )
        return op.ops[0]
    return op


# ─────────────────────────────────────────────────────────────────────
# Fixtures
# ─────────────────────────────────────────────────────────────────────


@pytest.fixture
def quad_2d():
    """LebedevSphere(17) — the canonical 2-D Cartesian quadrature."""
    return Quadrature.lebedev(17)


@pytest.fixture
def quad_1d():
    """GaussLegendre1D(8) — 1-D Cartesian / curvilinear quadrature."""
    return Quadrature.gauss_legendre(8)


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
    sn = SNMesh(mesh, quad_2d, placeholder_materials())
    assert isinstance(sn.bc_xmin, _BoundBoundaryOperator)
    assert sn.bc_xmin == "vacuum"

    rng = np.random.default_rng(42)
    psi = rng.uniform(0.5, 2.0, size=(quad_2d.N, 3, 2))
    out = sn.bc_xmin.apply(psi)
    inflow = np.flatnonzero(quad_2d.mu_x > 1e-12)
    non_inflow = np.setdiff1d(np.arange(quad_2d.N), inflow)
    np.testing.assert_array_equal(out[inflow], 0.0)
    np.testing.assert_array_equal(out[non_inflow], psi[non_inflow])


def test_2d_cartesian_reflective_ymax_returns_permutation(quad_2d):
    """Reflective on ymax: the realized shim's apply returns
    ``psi[ref]`` where ref = quad.reflection_index("y"). Equivalent
    to the legacy ReflectiveBoundary(axis="y") output.
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, 5), edges_y=np.linspace(0, 1, 4),
        mat_map=np.zeros((4, 3), dtype=int),
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    sn = SNMesh(mesh, quad_2d, placeholder_materials())
    assert isinstance(sn.bc_ymax, _BoundBoundaryOperator)
    assert sn.bc_ymax == "reflective"

    rng = np.random.default_rng(1)
    psi = rng.standard_normal(size=(quad_2d.N, 4, 2))
    ref = quad_2d.reflection_index("y")
    expected = psi[ref]
    np.testing.assert_array_equal(sn.bc_ymax.apply(psi), expected)


def test_2d_cartesian_construction_populates_trace(quad_2d):
    """SNMesh on a Cartesian mesh populates the unified :attr:`_trace`.
    Carries the per-face inflow indices used by the realizer.
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, 5), edges_y=np.linspace(0, 1, 4),
        mat_map=np.zeros((4, 3), dtype=int),
    )
    sn = SNMesh(mesh, quad_2d, placeholder_materials())
    assert sn._trace is not None
    assert sn._trace.face_names == ("xmin", "xmax", "ymin", "ymax")
    # xmin: inflow ordinates have mu_x > 0
    np.testing.assert_array_equal(
        sn._trace.inflow_indices_for_face("xmin"),
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
    sn = SNMesh(mesh, quad_1d, placeholder_materials())
    assert isinstance(sn.bc_right, _BoundBoundaryOperator)
    assert sn.bc_right == "vacuum"

    psi = np.arange(quad_1d.N * 2, dtype=float).reshape(quad_1d.N, 2)
    out = sn.bc_right.apply(psi)
    inflow = np.flatnonzero(quad_1d.mu_x < -1e-12)
    non_inflow = np.setdiff1d(np.arange(quad_1d.N), inflow)
    np.testing.assert_array_equal(out[inflow], 0.0)
    np.testing.assert_array_equal(out[non_inflow], psi[non_inflow])


def test_1d_cartesian_y_face_placeholders_realized_through_minimal_space(quad_1d):
    """Issue #188 / C188.3: 1-D meshes route the degenerate y-face
    :class:`ReflectiveBoundary` placeholders through
    :class:`SNBoundaryRealizer` with :meth:`SNMethodSpace.minimal` —
    the 1-D trace space cannot service
    ``inflow_indices_for_face('ymin')`` (its face_names are
    ``("xmin", "xmax")`` only) but the realizer's
    :class:`ReflectiveBoundary` branch does NOT read inflow_indices,
    only ``law.axis`` and ``quad.reflection_index``. For
    :class:`GaussLegendre1D`, ``reflection_index("y")`` returns the
    identity permutation (every ordinate is its own partner because
    ``mu_y == 0``), so the realized op is an identity
    :class:`PermutationOperator` — a no-op 1-arg
    :class:`LinearOperator`. The shim carries no bound quadrature
    (``_quadrature is None``) since the realized op is already 1-arg.
    """
    mesh = Mesh1D(
        edges=np.linspace(0, 1, 5),
        mat_ids=np.zeros(4, dtype=int),
    )
    sn = SNMesh(mesh, quad_1d, placeholder_materials())
    assert isinstance(sn.bc_ymin, _BoundBoundaryOperator)
    assert isinstance(sn.bc_ymax, _BoundBoundaryOperator)
    # Inner is the REALIZED PermutationOperator (NOT the raw law),
    # post-T.1 lifted into ``Permutation ⊗ Identity``.
    ymin_perm = _angular_factor(sn.bc_ymin.inner)
    ymax_perm = _angular_factor(sn.bc_ymax.inner)
    assert isinstance(ymin_perm, PermutationOperator)
    assert isinstance(ymax_perm, PermutationOperator)
    # The permutation is the identity for GL1D y-reflection.
    np.testing.assert_array_equal(ymin_perm.perm, np.arange(quad_1d.N))
    np.testing.assert_array_equal(ymax_perm.perm, np.arange(quad_1d.N))
    # No bound quadrature — the realized op is already 1-arg.
    # Issue #176 / C176.1 dropped the _quadrature attribute entirely.
    assert not hasattr(sn.bc_ymin, "_quadrature")
    assert not hasattr(sn.bc_ymax, "_quadrature")
    # Kind preserved for the legacy string-compare surface.
    assert sn.bc_ymin == "reflective"
    assert sn.bc_ymax == "reflective"
    # No-op identity apply: any psi goes through unchanged.
    psi = np.arange(quad_1d.N * 3, dtype=float).reshape(quad_1d.N, 3)
    np.testing.assert_array_equal(sn.bc_ymin.apply(psi), psi)
    np.testing.assert_array_equal(sn.bc_ymax.apply(psi), psi)


# ─────────────────────────────────────────────────────────────────────
# 1-D curvilinear: realizer path (Issue #188 / C188.3)
# ─────────────────────────────────────────────────────────────────────


def test_1d_spherical_vacuum_routes_through_realizer(quad_1d):
    """Spherical vacuum routes through :class:`SNBoundaryRealizer`. A
    solid sphere has exactly ONE boundary — the outer radius
    (``xmax``); the pole r=0 is the angular closure's regularity
    condition, not a BC face. The unified :class:`TraceSpace` therefore
    carries only the ``xmax`` face, and ``bc_left`` / ``bc_xmin`` are
    ``None`` (no inner boundary). The realizer's vacuum branch returns
    an :class:`IncomingOrdinateMaskTensor` over the per-face inflow
    indices (mu_x < 0 at the outer face). Per §16A.5 the mask zeros
    ONLY the inflow rows, leaving outflow rows untouched.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.1, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_right=BC("vacuum"),
    )
    sn = SNMesh(mesh, quad_1d, placeholder_materials())
    # Realizer path: shim wraps a realized 1-arg op, post-T.1 lifted
    # into ``IncomingOrdinateMaskTensor ⊗ Identity``.
    assert isinstance(sn.bc_right, _BoundBoundaryOperator)
    assert isinstance(
        _angular_factor(sn.bc_right.inner), IncomingOrdinateMaskTensor
    )
    # Issue #176 / C176.1 dropped the _quadrature attribute entirely.
    assert not hasattr(sn.bc_right, "_quadrature")
    assert sn.bc_right == "vacuum"
    # ONE boundary: no inner-face operator at the pole.
    assert sn.bc_left is None
    assert sn.bc_xmin is None
    # The unified trace carries only the outer ``xmax`` face.
    assert sn._trace is not None
    assert sn._trace.face_names == ("xmax",)
    # The outer-face inflow indices are the mu_x < 0 ordinates
    # (Ω · n_xmax = +mu_x; inflow iff Ω · n < 0).
    expected_inflow = np.flatnonzero(quad_1d.mu_x < -1e-12)
    np.testing.assert_array_equal(
        sn._trace.inflow_indices_for_face("xmax"),
        expected_inflow,
    )

    # §16A.5 inflow-only mask: zeros inflow rows, preserves outflow.
    psi = np.arange(quad_1d.N * 2, dtype=float).reshape(quad_1d.N, 2) + 1.0
    out = sn.bc_right.apply(psi)
    non_inflow = np.setdiff1d(np.arange(quad_1d.N), expected_inflow)
    np.testing.assert_array_equal(out[expected_inflow], 0.0)
    np.testing.assert_array_equal(out[non_inflow], psi[non_inflow])


def test_1d_cylindrical_one_boundary_outer_reflective():
    """A solid cylinder has ONE boundary — the outer radius (``xmax``).
    Any ``bc_left`` declaration at the pole r=0 is moot: the centreline
    is a geometry-forced symmetry handled by the angular pole closure,
    not an externally-imposed BC. So ``bc_left`` / ``bc_xmin`` are
    ``None``, and only the outer reflective BC is realized. The
    :class:`ReflectiveBoundary` branch produces a
    :class:`PermutationOperator` over ``quad.reflection_index("x")``;
    the shim wraps it with no bound quadrature.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.1, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"), bc_right=BC("reflective"),
    )
    quad = Quadrature.level_symmetric(sn_order=4)
    sn = SNMesh(mesh, quad, placeholder_materials())
    # ONE boundary: no inner-face operator at the pole (the bc_left
    # declaration on the mesh is ignored — the axis is the pole
    # closure's regularity condition, always symmetric by geometry).
    assert sn.bc_left is None
    assert sn.bc_xmin is None
    assert sn._trace is not None
    assert sn._trace.face_names == ("xmax",)
    # Outer face: realizer path; shim wraps a realized 1-arg op,
    # post-T.1 lifted into ``Permutation ⊗ Identity``.
    assert isinstance(sn.bc_right, _BoundBoundaryOperator)
    outer_perm = _angular_factor(sn.bc_right.inner)
    assert isinstance(outer_perm, PermutationOperator)
    assert not hasattr(sn.bc_right, "_quadrature")
    assert sn.bc_right == "reflective"
    np.testing.assert_array_equal(
        outer_perm.perm, quad.reflection_index("x"),
    )

    # Bit-equivalence: the shim's 1-arg apply matches the
    # ReflectiveBoundary semantics — psi[reflection_index].
    rng = np.random.default_rng(7)
    psi = rng.standard_normal(size=(quad.N, 2))
    np.testing.assert_array_equal(
        sn.bc_right.apply(psi),
        psi[quad.reflection_index("x")],
    )


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
        SNMesh(mesh, quad_1d, placeholder_materials())
