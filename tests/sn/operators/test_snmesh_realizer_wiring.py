r"""Tests for the SNMesh BC-resolution wiring (Wave 8 + C188.3 + C4 + #290 P7b).

The Wave-8 SNMesh routes BC resolution through
:class:`SNBoundaryRealizer` for every supported mesh. C4 (#220) made
the resolution surface dimension-generic: ONE loop over the face
labels populates the face-name-keyed :attr:`SNMesh.bc` dict
(``sn.bc["xmin"]`` …), whose keys equal
``boundary_face_layout.faces`` by construction (both derived from
``face_labels`` through the single-sourced
:attr:`FaceLabel.face_name` crosswalk). Since #290 P7b that loop is
the ONE shared ``TransportMethod`` body
(:func:`orpheus.transport.method.resolve_boundary_conditions`), and
``SNMesh.realize_boundary_law`` is the SN arm it dispatches. Each
entry is a :class:`_BoundBoundaryOperator` shim wrapping the 1-arg
realized :class:`LinearOperator`. The pre-C4 named attributes
(``bc_xmin`` … ``bc_ymax``, ``bc_left`` / ``bc_right`` aliases,
degenerate 1-D y-placeholders) are retired — negatives pinned below.

Issue #188 / C188.3: the curvilinear bypass branch in
``_resolve_one`` is gone. With the unified
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`'s curvilinear
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

from orpheus.geometry.boundary import ReflectiveBoundary, VacuumInflow
from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.geometry.boundary._bound_compat import _BoundBoundaryOperator
from orpheus.numerics.operator import (
    IdentityOperator,
    LinearOperator,
    PermutationOperator,
    TensorProductOperator,
    ZeroMorphism,
)
from orpheus.sn.mesh.augmented_mesh import SNMesh
from orpheus.numerics.quadrature import Quadrature
from tests._harness.references import mirror_partner_indices
from tests.sn._test_helpers import local_positions, placeholder_materials


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
# 2-D Cartesian: realizer wiring + the narrowed zero map Γ₊ → Γ₋
# ─────────────────────────────────────────────────────────────────────


def test_2d_cartesian_vacuum_xmin_is_the_zero_map(quad_2d):
    r"""Vacuum on xmin: the realized shim maps :math:`\Gamma_+ \to \Gamma_-`
    and its whole image is zero.

    RE-POSED at campaign phase **B3.2** (C-1). The pre-B3.2 claim was
    "zeros ONLY the inflow ordinates; outflow + tangential pass through
    unchanged" — the ``IncomingOrdinateMaskTensor`` contract. With the law's
    domain narrowed there is nothing to pass through: vacuum's entire content
    is :math:`R = 0`, so it realizes as the honest zero map and the preserved
    rows leave the picture rather than being justified.
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, 5), edges_y=np.linspace(0, 1, 4),
        mat_map=np.zeros((4, 3), dtype=int),
        bc_xmin=BC("vacuum"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    sn = SNMesh(mesh, quad_2d, placeholder_materials())
    assert isinstance(sn.bc["xmin"], _BoundBoundaryOperator)
    assert isinstance(sn.bc["xmin"].law, VacuumInflow)

    rng = np.random.default_rng(42)
    # S4-amendment: the probe is an HONEST Γ₊ element — (ng=1, ny=3)
    # trailing, the declared trace shape (the bound ZeroMorphism emits its
    # DECLARED codomain zero, not an echo of arbitrary probe trailing).
    psi = rng.uniform(0.5, 2.0, size=(quad_2d.N, 1, 3))
    inflow = np.flatnonzero(quad_2d.mu_x > 1e-12)
    outflow = np.flatnonzero(quad_2d.mu_x < -1e-12)
    out = sn.bc["xmin"].apply(psi[outflow])
    assert out.shape == (inflow.size, 1, 3), (
        f"vacuum on xmin emitted {out.shape}; the narrowed codomain is Γ₋."
    )
    np.testing.assert_array_equal(out, 0.0)


def test_2d_cartesian_reflective_ymax_returns_narrowed_permutation(quad_2d):
    r"""Reflective on ymax: the realized shim returns ``psi[ref][inflow]``,
    the pre-B3.2 full-face gather RESTRICTED to :math:`\Gamma_-`.

    RE-POSED at B3.2. The reference is the pre-B3.2 full-face gather with the
    inflow rows selected — which is exactly the bit-identity claim the phase
    makes at the mesh-wired shim — and its partner map is the independent
    geometric reference (§7d.3), built without ``to_local`` and without any
    production pairing derivation, so neither a remap error nor a pairing
    error can cancel against it.
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, 5), edges_y=np.linspace(0, 1, 4),
        mat_map=np.zeros((4, 3), dtype=int),
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    sn = SNMesh(mesh, quad_2d, placeholder_materials())
    assert isinstance(sn.bc["ymax"], _BoundBoundaryOperator)
    assert isinstance(sn.bc["ymax"].law, ReflectiveBoundary)

    rng = np.random.default_rng(1)
    psi = rng.standard_normal(size=(quad_2d.N, 4, 2))
    inflow = sn.angular_trace.inflow_indices_for_face("ymax")
    outflow = sn.angular_trace.outflow_indices_for_face("ymax")
    expected = psi[mirror_partner_indices(quad_2d, "y")][inflow]
    np.testing.assert_array_equal(sn.bc["ymax"].apply(psi[outflow]), expected)


def test_2d_reflective_y_face_builds_y_axis_permutation(quad_2d):
    r"""A y-face's reflective law reflects across the Y axis — structurally.

    The pre-C4 ``_resolve_one`` mapped every non-y face to axis ``"x"``
    via a hand-listed membership test (``"y" if face in ("ymin",
    "ymax") else "x"``) — correct at d≤2 by string coincidence, but a
    z-face would have silently built the X-axis permutation (the wrong
    reflection partner — vv Mode-9 class). Post-C4 the axis IS the
    label's own ``AXIS_NAMES[axis_index]``.

    RE-POSED at **B3.2**: the realized permutation now lives on the REDUCED
    ordinate axis, so the pinned value is
    ``local_positions(mirror_partner[inflow], outflow)`` rather than the
    full-face partner map itself. The non-vacuity guard is re-posed with
    it, and STRENGTHENED — the pre-B3.2 guard ("x and y reflection maps
    differ") is no longer sufficient, because the narrowing could in principle
    collapse two distinct full-face maps onto the same reduced one. The guard
    asked here is the one that actually matters: on THIS face, does the WRONG
    axis produce a different narrowed permutation (or refuse outright)?
    """
    mesh = Mesh2D(
        edges_x=np.linspace(0, 1, 5), edges_y=np.linspace(0, 1, 4),
        mat_map=np.zeros((4, 3), dtype=int),
        bc_xmin=BC("reflective"), bc_xmax=BC("reflective"),
        bc_ymin=BC("reflective"), bc_ymax=BC("reflective"),
    )
    sn = SNMesh(mesh, quad_2d, placeholder_materials())
    for face, axis, wrong_axis in (
        ("ymin", "y", "x"), ("ymax", "y", "x"), ("xmin", "x", "y"),
    ):
        inflow = sn.angular_trace.inflow_indices_for_face(face)
        outflow = sn.angular_trace.outflow_indices_for_face(face)
        expected = local_positions(
            mirror_partner_indices(quad_2d, axis)[inflow], outflow,
        )
        # Non-vacuity, per face and IN THE NARROWED COORDINATES: the wrong
        # axis must be distinguishable here, else the pin proves nothing.
        try:
            wrong = local_positions(
                mirror_partner_indices(quad_2d, wrong_axis)[inflow], outflow,
            )
        except KeyError:
            pass  # the wrong axis maps Γ₋ off Γ₊ entirely — production refuses
        else:
            if np.array_equal(wrong, expected):
                pytest.fail(
                    f"vacuous pin at {face!r}: the {axis!r} and {wrong_axis!r} "
                    f"reflection maps agree once narrowed to this face's "
                    f"half-traces."
                )
        perm = _angular_factor(sn.bc[face].inner)
        assert isinstance(perm, PermutationOperator)
        np.testing.assert_array_equal(perm.perm, expected)


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


def test_1d_cartesian_vacuum_right_is_the_zero_map(quad_1d):
    r"""1-D slab with vacuum on right: the realized shim is the zero map
    :math:`\Gamma_+ \to \Gamma_-`.

    RE-POSED at **B3.2** from "zeros only the inflow rows, passes the rest
    through" — see the 2-D sibling for why the pass-through claim retired.
    """
    mesh = Mesh1D(
        edges=np.linspace(0, 2, 9),
        mat_ids=np.zeros(8, dtype=int),
        bc_left=BC("reflective"), bc_right=BC("vacuum"),
    )
    sn = SNMesh(mesh, quad_1d, placeholder_materials())
    assert isinstance(sn.bc["xmax"], _BoundBoundaryOperator)
    assert isinstance(sn.bc["xmax"].law, VacuumInflow)

    # S4-amendment: an honest Γ₊ element carries the declared ng=1 tail.
    psi = np.arange(quad_1d.N, dtype=float).reshape(quad_1d.N, 1)
    inflow = np.flatnonzero(quad_1d.mu_x < -1e-12)
    outflow = np.flatnonzero(quad_1d.mu_x > +1e-12)
    out = sn.bc["xmax"].apply(psi[outflow])
    assert out.shape == (inflow.size, 1)
    np.testing.assert_array_equal(out, 0.0)


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
        Quadrature.folded_product(n_mu=4, n_phi=8), placeholder_materials(),
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
    r"""Spherical vacuum routes through :class:`SNBoundaryRealizer`. A
    solid sphere has exactly ONE boundary — the outer radius
    (``xmax``); the pole r=0 is the angular closure's regularity
    condition, not a BC face. The unified :class:`AngularTraceSpace` therefore
    carries only the ``xmax`` face, and the ``bc`` dict has NO pole
    entry (structurally absent).

    RE-POSED at **B3.2**: the realizer's vacuum branch returns the zero map
    :math:`\Gamma_+ \to \Gamma_-` (a :class:`ZeroMorphism`), not an
    ``IncomingOrdinateMaskTensor`` lifted into a tensor product. The
    routing claim — curvilinear meshes go through the SAME realizer path as
    Cartesian ones — is what this test is for and it is unchanged; only the
    object at the end of that path moved.
    """
    mesh = Mesh1D(
        edges=np.linspace(0.1, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.SPHERICAL,
        bc_right=BC("vacuum"),
    )
    sn = SNMesh(mesh, quad_1d, placeholder_materials())
    # Realizer path: shim wraps a realized 1-arg op — since B3.2 the honest
    # ``Γ₊ → Γ₋`` zero map (no tensor-product lift: there is no full-face
    # projector left to decompose).
    assert isinstance(sn.bc["xmax"], _BoundBoundaryOperator)
    assert isinstance(sn.bc["xmax"].inner, ZeroMorphism)
    # Issue #176 / C176.1 dropped the _quadrature attribute entirely.
    assert not hasattr(sn.bc["xmax"], "_quadrature")
    assert isinstance(sn.bc["xmax"].law, VacuumInflow)
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

    # B3.2: the zero map onto Γ₋ — the whole image, not a masked full face.
    # S4-amendment: an honest Γ₊ element carries the declared ng=1 tail.
    psi = np.arange(quad_1d.N, dtype=float).reshape(quad_1d.N, 1) + 1.0
    outflow = sn._trace.outflow_indices_for_face("xmax")
    out = sn.bc["xmax"].apply(psi[outflow])
    assert out.shape == (expected_inflow.size, 1)
    np.testing.assert_array_equal(out, 0.0)


def test_1d_cylindrical_one_boundary_outer_reflective():
    """A solid cylinder has ONE boundary — the outer radius (``xmax``).
    Any ``bc_left`` declaration at the pole r=0 is moot: the centreline
    is a geometry-forced symmetry handled by the angular pole closure,
    not an externally-imposed BC. So the ``bc`` dict has no pole
    entry, and only the outer reflective BC is realized. The
    :class:`ReflectiveBoundary` branch produces a
    :class:`PermutationOperator`; the shim wraps it with no bound quadrature.

    RE-POSED at **B3.2**: the permutation is on the REDUCED ordinate axis, so
    its table is ``local_positions(mirror_partner[inflow], outflow)`` and its
    length is :math:`|\\Gamma_+|`. Both the structural pin and the value pin
    move with it; the reference stays the pre-B3.2 full-face gather,
    restricted, with the partner map from the independent geometric
    reference (§7d.3).
    """
    mesh = Mesh1D(
        edges=np.linspace(0.1, 1.0, 6),
        mat_ids=np.zeros(5, dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        bc_left=BC("reflective"), bc_right=BC("reflective"),
    )
    quad = Quadrature.folded_product(n_mu=4, n_phi=8)
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
    assert isinstance(sn.bc["xmax"].law, ReflectiveBoundary)
    inflow = sn._trace.inflow_indices_for_face("xmax")
    outflow = sn._trace.outflow_indices_for_face("xmax")
    np.testing.assert_array_equal(
        outer_perm.perm,
        local_positions(mirror_partner_indices(quad, "x")[inflow], outflow),
    )
    assert outer_perm.perm.size == outflow.size < quad.N

    # Bit-equivalence: the shim's 1-arg apply matches the
    # ReflectiveBoundary semantics — ``psi[mirror_partner]`` RESTRICTED to
    # Γ₋, which is the pre-B3.2 value this face's consumer actually read.
    rng = np.random.default_rng(7)
    psi = rng.standard_normal(size=(quad.N, 2))
    np.testing.assert_array_equal(
        sn.bc["xmax"].apply(psi[outflow]),
        psi[mirror_partner_indices(quad, "x")][inflow],
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
