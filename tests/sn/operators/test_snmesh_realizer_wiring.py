r"""Tests for ``SNMesh._resolve_bcs`` Wave-8 + C188.3 + C4 realizer wiring.

The Wave-8 SNMesh routes BC resolution through
:class:`SNBoundaryRealizer` for every supported mesh. C4 (#220) made
the resolution surface dimension-generic: ONE loop over
:attr:`SNMesh.face_labels` populates the face-name-keyed
:attr:`SNMesh.bc` dict (``sn.bc["xmin"]`` …), whose keys equal
``boundary_face_layout.faces`` by construction (both derived from
``face_labels`` through the single-sourced
:attr:`FaceLabel.face_name` crosswalk). Each entry is a
:class:`_BoundBoundaryOperator` shim wrapping the 1-arg realized
:class:`LinearOperator`. The pre-C4 named attributes (``bc_xmin`` …
``bc_ymax``, ``bc_left`` / ``bc_right`` aliases, degenerate 1-D
y-placeholders) are retired — negatives pinned below.

Issue #188 / C188.3: the curvilinear bypass branch in
``_resolve_one`` is gone. With the unified
:class:`~orpheus.numerics.spaces.trace_space.TraceSpace`'s curvilinear
support, 1-D spherical and 1-D cylindrical meshes route through the
SAME realizer-then-shim path as Cartesian meshes — but a solid
sphere / cylinder has only the outer (``xmax``) boundary face; the
pole r=0 is the angular closure's regularity condition, so the
``bc`` dict has NO pole entry (structurally absent, not ``None``).

V&V tags
--------
``@pytest.mark.l1`` — the wiring assertions are cross-implementation
checks (Wave-5 realizer dispatch + Wave-2 trace-mask construction +
Wave-8 shim composition produce the same observable per-face apply).

Verification design (C4):
``.claude/agent-memory/test-architect/c4_snmesh_bc_dict_verification.md``.
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
from orpheus.sn.mesh.augmented_mesh import SNMesh
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
    assert isinstance(sn.bc["xmin"], _BoundBoundaryOperator)
    assert sn.bc["xmin"] == "vacuum"

    rng = np.random.default_rng(42)
    psi = rng.uniform(0.5, 2.0, size=(quad_2d.N, 3, 2))
    out = sn.bc["xmin"].apply(psi)
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
    assert isinstance(sn.bc["ymax"], _BoundBoundaryOperator)
    assert sn.bc["ymax"] == "reflective"

    rng = np.random.default_rng(1)
    psi = rng.standard_normal(size=(quad_2d.N, 4, 2))
    ref = quad_2d.reflection_index("y")
    expected = psi[ref]
    np.testing.assert_array_equal(sn.bc["ymax"].apply(psi), expected)


def test_2d_reflective_y_face_builds_y_axis_permutation(quad_2d):
    """A y-face's reflective law reflects across the Y axis — structurally.

    The pre-C4 ``_resolve_one`` mapped every non-y face to axis ``"x"``
    via a hand-listed membership test (``"y" if face in ("ymin",
    "ymax") else "x"``) — correct at d≤2 by string coincidence, but a
    z-face would have silently built the X-axis permutation (the wrong
    reflection partner — vv Mode-9 class). Post-C4 the axis IS the
    label's own ``AXIS_NAMES[axis_index]``. Pin the d=2 observable:
    the realized ymin permutation equals ``reflection_index("y")``,
    and the x/y permutations differ under Lebedev (else the pin would
    be vacuous).
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, 5), edges_y=np.linspace(0, 1, 4),
        mat_map=np.zeros((4, 3), dtype=int),
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    sn = SNMesh(mesh, quad_2d, placeholder_materials())
    if np.array_equal(
        quad_2d.reflection_index("x"), quad_2d.reflection_index("y")
    ):
        pytest.fail("vacuous pin: x and y reflection maps coincide")
    for face, axis in (("ymin", "y"), ("ymax", "y"), ("xmin", "x")):
        perm = _angular_factor(sn.bc[face].inner)
        assert isinstance(perm, PermutationOperator)
        np.testing.assert_array_equal(
            perm.perm, quad_2d.reflection_index(axis),
        )


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
# 1-D Cartesian: realizer wiring through bc["xmin"] / bc["xmax"]
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
    assert isinstance(sn.bc["xmax"], _BoundBoundaryOperator)
    assert sn.bc["xmax"] == "vacuum"

    psi = np.arange(quad_1d.N * 2, dtype=float).reshape(quad_1d.N, 2)
    out = sn.bc["xmax"].apply(psi)
    inflow = np.flatnonzero(quad_1d.mu_x < -1e-12)
    non_inflow = np.setdiff1d(np.arange(quad_1d.N), inflow)
    np.testing.assert_array_equal(out[inflow], 0.0)
    np.testing.assert_array_equal(out[non_inflow], psi[non_inflow])


# ─────────────────────────────────────────────────────────────────────
# C4 (#220): the bc-dict inventory IS the face inventory
# ─────────────────────────────────────────────────────────────────────


def test_bc_inventory_equals_face_layout_across_geometries(quad_1d, quad_2d):
    """``set(sn.bc) == set(boundary_face_layout.faces)`` — by construction.

    Both sides derive from ``face_labels`` through the SAME
    ``FaceLabel.face_name`` crosswalk, so the BC inventory and the
    face layout cannot drift. Entry counts: slab 2, 2-D Cartesian 4,
    sphere 1, cylinder 1 (issue #220 acceptance). The 1-D inventories
    carry NO y-entries — the pre-C4 degenerate y-placeholders (a
    realized no-op ``ReflectiveBoundary(axis="y")`` pair no production
    code ever read) are retired, dict misses fail loud below.
    """
    slab = SNMesh(
        Mesh1D(edges=np.linspace(0, 1, 5), mat_ids=np.zeros(4, dtype=int)),
        quad_1d, placeholder_materials(),
    )
    two_d = SNMesh(
        Mesh2D(edges_x=np.linspace(0, 1, 5), edges_y=np.linspace(0, 1, 4),
               mat_map=np.zeros((4, 3), dtype=int)),
        quad_2d, placeholder_materials(),
    )
    sphere = SNMesh(
        Mesh1D(edges=np.linspace(0.1, 1.0, 6), mat_ids=np.zeros(5, dtype=int),
               coord=CoordSystem.SPHERICAL),
        quad_1d, placeholder_materials(),
    )
    cylinder = SNMesh(
        Mesh1D(edges=np.linspace(0.1, 1.0, 6), mat_ids=np.zeros(5, dtype=int),
               coord=CoordSystem.CYLINDRICAL),
        Quadrature.level_symmetric(sn_order=4), placeholder_materials(),
    )
    expected = {
        "slab": (slab, {"xmin", "xmax"}),
        "2d": (two_d, {"xmin", "xmax", "ymin", "ymax"}),
        "sphere": (sphere, {"xmax"}),
        "cylinder": (cylinder, {"xmax"}),
    }
    for name, (sn, faces) in expected.items():
        if set(sn.bc) != faces:
            pytest.fail(f"{name}: bc keys {set(sn.bc)} != {faces}")
        if set(sn.bc) != set(sn.boundary_face_layout.faces):
            pytest.fail(f"{name}: bc keys drift from boundary_face_layout")


def test_bc_dict_misses_and_retired_attributes_fail_loud(quad_1d):
    """Negatives: a face that doesn't exist is a ``KeyError`` (plain
    dict — no masking default), and the retired named attributes are
    ``AttributeError`` (no silent ``None``-shim survives; a read-through
    ``@property`` reappearing would be a deprecation outliving its
    cycle).
    """
    slab = SNMesh(
        Mesh1D(edges=np.linspace(0, 1, 5), mat_ids=np.zeros(4, dtype=int)),
        quad_1d, placeholder_materials(),
    )
    with pytest.raises(KeyError):
        slab.bc["ymin"]
    for retired in (
        "bc_left", "bc_right", "bc_xmin", "bc_xmax", "bc_ymin", "bc_ymax",
    ):
        with pytest.raises(AttributeError):
            getattr(slab, retired)


# ─────────────────────────────────────────────────────────────────────
# 1-D curvilinear: realizer path (Issue #188 / C188.3)
# ─────────────────────────────────────────────────────────────────────


def test_1d_spherical_vacuum_routes_through_realizer(quad_1d):
    """Spherical vacuum routes through :class:`SNBoundaryRealizer`. A
    solid sphere has exactly ONE boundary — the outer radius
    (``xmax``); the pole r=0 is the angular closure's regularity
    condition, not a BC face. The unified :class:`TraceSpace` therefore
    carries only the ``xmax`` face, and the ``bc`` dict has NO pole
    entry (structurally absent). The realizer's vacuum branch returns
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
    assert isinstance(sn.bc["xmax"], _BoundBoundaryOperator)
    assert isinstance(
        _angular_factor(sn.bc["xmax"].inner), IncomingOrdinateMaskTensor
    )
    # Issue #176 / C176.1 dropped the _quadrature attribute entirely.
    assert not hasattr(sn.bc["xmax"], "_quadrature")
    assert sn.bc["xmax"] == "vacuum"
    # ONE boundary: no inner-face entry at the pole (C4 — structurally
    # absent, not None).
    assert set(sn.bc) == {"xmax"}
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
    out = sn.bc["xmax"].apply(psi)
    non_inflow = np.setdiff1d(np.arange(quad_1d.N), expected_inflow)
    np.testing.assert_array_equal(out[expected_inflow], 0.0)
    np.testing.assert_array_equal(out[non_inflow], psi[non_inflow])


def test_1d_cylindrical_one_boundary_outer_reflective():
    """A solid cylinder has ONE boundary — the outer radius (``xmax``).
    Any ``bc_left`` declaration at the pole r=0 is moot: the centreline
    is a geometry-forced symmetry handled by the angular pole closure,
    not an externally-imposed BC. So the ``bc`` dict has no pole
    entry, and only the outer reflective BC is realized. The
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
    # ONE boundary: no inner-face entry at the pole (the bc_left
    # declaration on the mesh is ignored — the axis is the pole
    # closure's regularity condition, always symmetric by geometry).
    assert set(sn.bc) == {"xmax"}
    assert sn._trace is not None
    assert sn._trace.face_names == ("xmax",)
    # Outer face: realizer path; shim wraps a realized 1-arg op,
    # post-T.1 lifted into ``Permutation ⊗ Identity``.
    assert isinstance(sn.bc["xmax"], _BoundBoundaryOperator)
    outer_perm = _angular_factor(sn.bc["xmax"].inner)
    assert isinstance(outer_perm, PermutationOperator)
    assert not hasattr(sn.bc["xmax"], "_quadrature")
    assert sn.bc["xmax"] == "reflective"
    np.testing.assert_array_equal(
        outer_perm.perm, quad.reflection_index("x"),
    )

    # Bit-equivalence: the shim's 1-arg apply matches the
    # ReflectiveBoundary semantics — psi[reflection_index].
    rng = np.random.default_rng(7)
    psi = rng.standard_normal(size=(quad.N, 2))
    np.testing.assert_array_equal(
        sn.bc["xmax"].apply(psi),
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
