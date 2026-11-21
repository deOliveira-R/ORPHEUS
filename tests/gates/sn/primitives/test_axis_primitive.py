r"""Foundation pins for the :class:`Axis1D` primitive (R-1 Phase A C1).

F0.1 — :class:`AxisMesh` endpoints are ``("min", "max")``.
F0.2 — :class:`RadialAxisMesh` endpoints are ``("outer",)``.
F0.3 — :meth:`SNProblem.from_axes` 2-D shape + 4 face labels.
F0.4 — Spherical/cylindrical 1-D mesh has 1 face label (the pole is NOT a face).
F0.5 — :func:`face_shape` derivation matches per-axis cell counts.
F0.6 — :func:`face_outflow_ordinates` mask matches the inline expression at
       the 5+ legacy sites (``np.where(sign * mu_axis > 1e-15)[0]``).
F0.7 — Synthetic 3-D admission: ``(AxisMesh(5), AxisMesh(7), AxisMesh(9))``
       produces 6 face labels via the pure shape functions (no SNProblem).
F0.8 — The production full field ``SNProblem.full_field_space`` of a problem
       posed by :meth:`SNProblem.from_axes` has the closed-form size
       ``N·ng·∏cells + Σ_faces N·ng·∏(face cells)`` (every ordinate on every
       traced face) on a slab, a solid sphere, a 2-D and a 3-D Cartesian mesh,
       two groups each. The legs read the interior block, the trace block and
       their sum, so they pin both the trace's per-face ordinate and group
       counts and that ``from_axes`` poses on the caller's axes, quadrature and
       materials. ``[M]`` 2026-11-21: a ``from_axes`` that drops the caller's
       quadrature for a fixed one on a Cartesian mesh reddens the slab, 2-D and
       3-D legs; at d = 1 it reddens no test of the 13 other files that call
       ``from_axes`` (``-m "not slow"``). (Until the W0 re-pose these legs compared a
       test-side outflow-only counter, the retired ``n_unknowns_flat``, with a
       count of the same outflow sets; no production object has that size.)
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.data.macro_xs.synthetic import make_mixture
from orpheus.mesh import (
    Axis1D,
    AxisCoord,
    AxisMesh,
    FaceLabel,
    RadialAxisMesh,
    face_labels,
    face_outflow_ordinates,
    face_shape,
    spatial_shape,
)
from orpheus.sn.problem import SNProblem
from orpheus.numerics.quadrature import Quadrature
from orpheus.geometry.deck import MirrorDeck
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.vacuum()}



# Every F0.x pin here is a foundation invariant on the dim-agnostic
# axis primitive — it pins a software contract (axis-tuple shape
# algebra, face-label inventory, pack-vector sizing) rather than
# verifying a documented equation, so it carries @foundation and
# does NOT carry @verifies(...). See vv-principles §"V&V level
# taxonomy" — foundation is orthogonal to L0/L1/L2 and pins the
# software invariants that the verification ladder depends on.
pytestmark = [pytest.mark.foundation]

# The law every axis here declares where its subject is not the boundary: every
# consumer resolved an undeclared face to it before step 3c (#405).
_REFLECTIVE = MirrorDeck()


# ─── Fixtures ────────────────────────────────────────────────────────────


def _one_group_mixture():
    """A trivial 1-group, 1-region material used for shape-only tests."""
    return make_mixture(
        sig_t=np.array([1.0]),
        sig_c=np.array([0.5]),
        sig_f=np.array([0.0]),
        nu=np.array([0.0]),
        chi=np.zeros(1),  # non-fissile ⇒ null spectrum (S10a __post_init__ guard)
        sig_s=np.array([[0.5]]),
    )


def _level_symmetric_quad_2d(order: int = 4):
    """Minimal 2-D quadrature (LevelSymmetric S4)."""
    from orpheus.numerics.quadrature import Quadrature
    return Quadrature.level_symmetric(sn_order=order)


# ─── F0.1 ────────────────────────────────────────────────────────────────


def test_f0_1_axismesh_endpoints_min_max() -> None:
    """AxisMesh's two endpoints are ``min`` and ``max`` (Cartesian convention)."""
    ax = AxisMesh(edges=np.linspace(0.0, 1.0, 16), role_low=_REFLECTIVE, role_high=_REFLECTIVE)
    assert ax.endpoints == ("min", "max")
    assert ax.coord == AxisCoord.CARTESIAN
    assert ax.n == 15  # n derived from len(edges) - 1
    assert ax.roles == {"min": _REFLECTIVE, "max": _REFLECTIVE}  # the declaration, keyed by endpoint


def test_f0_1_axismesh_endpoints_relabelable() -> None:
    """Endpoint labels can be overridden for slab-conventional naming."""
    ax = AxisMesh(
        edges=np.linspace(0.0, 1.0, 6),
        label_low="left", label_high="right",
        role_low=0, role_high=MirrorDeck(),
    )
    assert ax.endpoints == ("left", "right")
    assert ax.roles == {"left": 0, "right": MirrorDeck()}


def test_f0_1_axismesh_monotonicity_validated() -> None:
    """Non-monotonic edges raise at construction."""
    with pytest.raises(ValueError, match="monotonically increasing"):
        AxisMesh(edges=np.array([0.0, 1.0, 0.5]), role_low=_REFLECTIVE, role_high=_REFLECTIVE)


# ─── F0.2 ────────────────────────────────────────────────────────────────


def test_f0_2_radial_axismesh_one_endpoint_outer() -> None:
    """RadialAxisMesh has exactly one endpoint ``outer`` (pole excluded)."""
    ax = RadialAxisMesh(
        edges=np.linspace(0.0, 2.0, 11),
        coord=AxisCoord.RADIAL_SPHERICAL,
        role_outer=0,
    )
    assert ax.endpoints == ("outer",)
    assert ax.n == 10
    assert ax.roles == {"outer": 0}


def test_f0_2_radial_axismesh_label_overridable() -> None:
    """The outer-endpoint label can be customised."""
    ax = RadialAxisMesh(
        edges=np.linspace(0.0, 1.0, 6),
        coord=AxisCoord.RADIAL_CYLINDRICAL,
        role_outer=0,
        label_outer="surface",
    )
    assert ax.endpoints == ("surface",)
    assert ax.roles == {"surface": 0}


def test_f0_2_radial_axismesh_rejects_cartesian_coord() -> None:
    """RadialAxisMesh refuses non-radial coords (Pattern 4)."""
    with pytest.raises(ValueError, match="RADIAL_SPHERICAL or RADIAL_CYLINDRICAL"):
        RadialAxisMesh(
            edges=np.linspace(0.0, 1.0, 6),
            coord=AxisCoord.CARTESIAN, role_outer=_REFLECTIVE,
        )


# ─── F0.3 ────────────────────────────────────────────────────────────────


def test_f0_3_from_axes_2d_cartesian_shape_and_face_labels() -> None:
    """SNProblem.from_axes 2-D Cart: shape (nx, ny) + 4 face labels."""
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 11), role_low=_REFLECTIVE, role_high=_REFLECTIVE),  # n=10
        AxisMesh(edges=np.linspace(0.0, 1.0, 8), role_low=_REFLECTIVE, role_high=_REFLECTIVE),   # n=7
    )
    mesh = SNProblem.from_axes(
        axes,
        quadrature=_level_symmetric_quad_2d(order=4),
        materials={0: _one_group_mixture()}, boundary_conditions=conditions_for(axes, _CONDITIONS),
    )
    assert mesh.ndim == 2
    assert mesh.spatial_shape == (10, 7)
    assert mesh.face_labels == (
        FaceLabel(0, "min"),
        FaceLabel(0, "max"),
        FaceLabel(1, "min"),
        FaceLabel(1, "max"),
    )


# ─── F0.4 ────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("coord,make_quad", [
    (AxisCoord.RADIAL_SPHERICAL, lambda: Quadrature.gauss_legendre(n_ordinates=8)),
    (AxisCoord.RADIAL_CYLINDRICAL, lambda: Quadrature.folded_product(n_mu=2, n_phi=4)),
])
def test_f0_4_solid_radial_mesh_has_one_face_label(coord, make_quad) -> None:
    """Sphere / cylinder SNProblem has exactly 1 face label (``outer``).

    Sphere uses the polar GL quadrature; cylinder requires a level-
    structured quadrature (the folded product family) because the
    cylindrical streaming reduction depends on level partitioning.
    """
    axes = (
        RadialAxisMesh(
            edges=np.linspace(0.0, 1.0, 21),
            coord=coord,
            role_outer=0,
        ),
    )
    mesh = SNProblem.from_axes(
        axes,
        quadrature=make_quad(),
        materials={0: _one_group_mixture()}, boundary_conditions=conditions_for(axes, _CONDITIONS),
    )
    assert mesh.ndim == 1
    assert mesh.spatial_shape == (20,)
    assert mesh.face_labels == (FaceLabel(0, "outer"),)
    # The pole is structurally absent — no FaceLabel(0, "pole") or
    # similar. The Carlson coupled-pole sweep in
    # MorelMontryAngularSweep carries the pole's regularity.


def test_f0_4_slab_1d_has_two_face_labels() -> None:
    """1-D Cartesian slab has 2 face labels (``min``, ``max``)."""
    axes = (AxisMesh(edges=np.linspace(0.0, 1.0, 16),
                     role_low=0, role_high=0),)
    mesh = SNProblem.from_axes(
        axes,
        quadrature=Quadrature.gauss_legendre(n_ordinates=8),
        materials={0: _one_group_mixture()}, boundary_conditions=conditions_for(axes, _CONDITIONS),
    )
    assert mesh.face_labels == (FaceLabel(0, "min"), FaceLabel(0, "max"))


# ─── F0.5 ────────────────────────────────────────────────────────────────


def test_f0_5_face_shape_2d_each_face() -> None:
    """face_shape: the face on axis-i has shape over the OTHER axis."""
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),   # n=5
        AxisMesh(edges=np.linspace(0.0, 1.0, 8), role_low=_REFLECTIVE, role_high=_REFLECTIVE),   # n=7
    )
    assert face_shape(axes, FaceLabel(0, "min")) == (7,)
    assert face_shape(axes, FaceLabel(0, "max")) == (7,)
    assert face_shape(axes, FaceLabel(1, "min")) == (5,)
    assert face_shape(axes, FaceLabel(1, "max")) == (5,)


def test_f0_5_face_shape_3d_synthetic() -> None:
    """3-D synthetic: face on axis-1 has shape over (axis-0, axis-2)."""
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),   # n=5
        AxisMesh(edges=np.linspace(0.0, 1.0, 8), role_low=_REFLECTIVE, role_high=_REFLECTIVE),   # n=7
        AxisMesh(edges=np.linspace(0.0, 1.0, 10), role_low=_REFLECTIVE, role_high=_REFLECTIVE),  # n=9
    )
    assert face_shape(axes, FaceLabel(0, "min")) == (7, 9)
    assert face_shape(axes, FaceLabel(1, "min")) == (5, 9)
    assert face_shape(axes, FaceLabel(2, "min")) == (5, 7)


def test_f0_5_face_shape_radial_solid_no_inner_face() -> None:
    """Solid radial mesh has only the outer face; shape is empty tuple."""
    axes = (RadialAxisMesh(
        edges=np.linspace(0.0, 1.0, 11),
        coord=AxisCoord.RADIAL_SPHERICAL, role_outer=0,
    ),)
    labels = face_labels(axes)
    assert labels == (FaceLabel(0, "outer"),)
    assert face_shape(axes, labels[0]) == ()


# ─── F0.6 ────────────────────────────────────────────────────────────────


def test_f0_6_face_outflow_ordinates_matches_inline_expression_1d() -> None:
    """Outflow ordinate mask matches np.where(sign * mu > 1e-15)."""
    axes = (AxisMesh(edges=np.linspace(0.0, 1.0, 11), role_low=_REFLECTIVE, role_high=_REFLECTIVE),)
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    out_max = face_outflow_ordinates(axes, FaceLabel(0, "max"), quad)
    out_min = face_outflow_ordinates(axes, FaceLabel(0, "min"), quad)
    np.testing.assert_array_equal(out_max, np.where(quad.mu_x > 1e-15)[0])
    np.testing.assert_array_equal(out_min, np.where(-quad.mu_x > 1e-15)[0])
    # In 1-D GL even-N the two are complementary halves of (0..N-1).
    assert out_max.size == out_min.size == quad.n_ordinates // 2


def test_f0_6_face_outflow_ordinates_2d_LS4() -> None:
    """2-D LS4: outflow at each face is the half-sphere of that axis."""
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),
        AxisMesh(edges=np.linspace(0.0, 1.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),
    )
    quad = _level_symmetric_quad_2d(order=4)
    out_xmax = face_outflow_ordinates(axes, FaceLabel(0, "max"), quad)
    out_xmin = face_outflow_ordinates(axes, FaceLabel(0, "min"), quad)
    out_ymax = face_outflow_ordinates(axes, FaceLabel(1, "max"), quad)
    out_ymin = face_outflow_ordinates(axes, FaceLabel(1, "min"), quad)
    np.testing.assert_array_equal(out_xmax, np.where(quad.mu_x > 1e-15)[0])
    np.testing.assert_array_equal(out_xmin, np.where(-quad.mu_x > 1e-15)[0])
    np.testing.assert_array_equal(out_ymax, np.where(quad.mu_y > 1e-15)[0])
    np.testing.assert_array_equal(out_ymin, np.where(-quad.mu_y > 1e-15)[0])


def test_f0_6_face_outflow_ordinates_radial() -> None:
    """Solid sphere outer-face outflow uses the polar GL cosines."""
    axes = (RadialAxisMesh(
        edges=np.linspace(0.0, 1.0, 11),
        coord=AxisCoord.RADIAL_SPHERICAL, role_outer=0,
    ),)
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    out_outer = face_outflow_ordinates(axes, FaceLabel(0, "outer"), quad)
    np.testing.assert_array_equal(out_outer, np.where(quad.mu_x > 1e-15)[0])


# ─── F0.7 ────────────────────────────────────────────────────────────────


def test_f0_7_synthetic_3d_admission_six_face_labels() -> None:
    """Synthetic 3-D axis tuple admits the dim-agnostic shape primitives.

    C1 ships the pure shape functions in ``orpheus.mesh.axis``; no
    ``Mesh3D`` dataclass exists yet (D9 of the ultraplan). The 3-D
    admission gate exercises these functions on a synthetic axis tuple
    without instantiating an :class:`SNProblem` — that lands in a
    followup once :class:`Mesh3D` is in tree.
    """
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),   # n=5
        AxisMesh(edges=np.linspace(0.0, 1.0, 8), role_low=_REFLECTIVE, role_high=_REFLECTIVE),   # n=7
        AxisMesh(edges=np.linspace(0.0, 1.0, 10), role_low=_REFLECTIVE, role_high=_REFLECTIVE),  # n=9
    )
    assert spatial_shape(axes) == (5, 7, 9)
    labels = face_labels(axes)
    assert labels == (
        FaceLabel(0, "min"), FaceLabel(0, "max"),
        FaceLabel(1, "min"), FaceLabel(1, "max"),
        FaceLabel(2, "min"), FaceLabel(2, "max"),
    )
    assert len(labels) == 6


def test_f0_7_synthetic_3d_face_outflow_axis_beyond_quad_dim_is_empty() -> None:
    r"""Axes beyond the quadrature's intrinsic dim have empty outflow.

    GL1D is the slab quadrature; its :class:`DiscreteMeasure` has
    ``nodes`` of shape ``(N,)`` (1-D scalar polar :math:`\mu`).
    Pairing GL1D with a 3-D axis tuple is unphysical, but the
    dim-agnostic shape primitive must handle it gracefully — for
    any axis index :math:`i > 0` on a 1-D measure (or :math:`i \ge
    d` on a d-dim measure), ``Quadrature.axis_cosines`` returns
    zeros, and ``np.where(\pm 1 \cdot 0 > 10^{-15})[0]`` is empty.

    This is the Pattern 7 verification: the shape primitive reads
    ``quad.measure.nodes`` directly (the single source of truth);
    there is no per-axis name dispatch and no getattr fallback. The
    SHAPE of the measure (``nodes.ndim`` and ``nodes.shape[1]`` if
    multi-dim) decides which axis indices return real cosines vs
    zeros.
    """
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),
        AxisMesh(edges=np.linspace(0.0, 1.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),
        AxisMesh(edges=np.linspace(0.0, 1.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    # GL1D's measure is 1-D scalar — only axis 0 carries real cosines.
    assert quad.measure.nodes.ndim == 1
    out_xmax = face_outflow_ordinates(axes, FaceLabel(0, "max"), quad)
    assert out_xmax.size == quad.n_ordinates // 2   # axis 0 is the polar μ — real cosines
    # Axes 1 and 2 are beyond the quadrature's intrinsic dim — outflow is empty.
    for axis_idx in (1, 2):
        out_max = face_outflow_ordinates(axes, FaceLabel(axis_idx, "max"), quad)
        out_min = face_outflow_ordinates(axes, FaceLabel(axis_idx, "min"), quad)
        assert out_max.size == 0
        assert out_min.size == 0


def test_f0_7_synthetic_3d_face_outflow_axis_2_ls4_native_mu_z() -> None:
    """LS4 (and other 3-D quadratures) already expose ``mu_z`` natively.

    For a 3-D axis tuple paired with LS4, axis-2 outflow counts the
    ordinates with ``mu_z`` strictly outward at each endpoint — the
    same shape contract as axis-0 / axis-1.
    """
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),
        AxisMesh(edges=np.linspace(0.0, 1.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),
        AxisMesh(edges=np.linspace(0.0, 1.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),
    )
    quad = _level_symmetric_quad_2d(order=4)
    out_zmax = face_outflow_ordinates(axes, FaceLabel(2, "max"), quad)
    out_zmin = face_outflow_ordinates(axes, FaceLabel(2, "min"), quad)
    np.testing.assert_array_equal(out_zmax, np.where(quad.mu_z > 1e-15)[0])
    np.testing.assert_array_equal(out_zmin, np.where(-quad.mu_z > 1e-15)[0])
    # LS4 is symmetric over the unit sphere: half the ordinates have
    # mu_z > 0, half have mu_z < 0.
    assert out_zmax.size == out_zmin.size == quad.n_ordinates // 2


# ─── F0.8 ────────────────────────────────────────────────────────────────
#
# The production full field of the angular flux is the direct sum
# ``V_bulk ⊕ V_trace`` (``SNProblem.full_field_space``): the interior carries
# every ordinate of every group in every cell, and EVERY traced face carries
# EVERY ordinate of every group on each of its face cells (inflow, outflow and
# tangential alike; the inflow / outflow split is a selector over the trace's
# signed Ω·n, not a size). Each leg poses the problem through
# ``SNProblem.from_axes`` and asserts the three sizes against a closed-form
# count written here from the cell counts, the face inventory and the
# quadrature's ordinate count; nothing in the count calls the code that builds
# the spaces. Two groups, so a dropped energy axis on either block is visible
# (one group nulls the group factor).


def _two_group_mixture():
    """A 2-group, 1-region, non-fissile mixture with downscatter; sizes only, gated on its own balance."""
    sig_t = np.array([1.0, 2.0])
    sig_c = np.array([0.5, 1.0])
    sig_s = np.array([[0.4, 0.1], [0.0, 1.0]])  # [from, to]
    # the mixture's consistency identity (``vv-testing``): sigma_t = sigma_c + sigma_f + sum_to SigS[from, to]
    np.testing.assert_allclose(sig_t, sig_c + sig_s.sum(axis=1), rtol=0, atol=0)
    return make_mixture(sig_t=sig_t, sig_c=sig_c, sig_f=np.zeros(2), nu=np.zeros(2), chi=np.zeros(2), sig_s=sig_s)


def _pose(axes: tuple[Axis1D, ...], quad: Quadrature) -> SNProblem:
    return SNProblem.from_axes(
        axes, quadrature=quad, materials={0: _two_group_mixture()},
        boundary_conditions=conditions_for(axes, _CONDITIONS),
    )


def _require_sizes(problem: SNProblem, *, interior: int, trace: int) -> None:
    """The interior block, the trace block and their direct sum, each against its closed-form count."""
    space = problem.full_field_space
    got = (space.interior_space.n_coordinates, space.trace_space.n_coordinates, space.n_coordinates)
    if got != (interior, trace, interior + trace):
        pytest.fail(f"full field (interior, trace, sum) reads {got}; the closed form is {(interior, trace, interior + trace)}")


@pytest.mark.rests_on(
    "tests/gates/transport/test_phase_spaces.py::test_the_face_layout_is_the_hand_layout",
    "tests/gates/mesh/test_trace_face_shapes.py::test_the_traced_faces_are_the_hand_inventory",
)
def test_f0_8_full_field_size_slab() -> None:
    """Slab, 10 cells, GL8 (8 ordinates), 2 groups: the two point faces each carry all 8 ordinates."""
    axes = (AxisMesh(edges=np.linspace(0.0, 1.0, 11), role_low=_REFLECTIVE, role_high=_REFLECTIVE),)
    N, ng = 8, 2  # a Gauss-Legendre rule of n points has n ordinates
    _require_sizes(
        _pose(axes, Quadrature.gauss_legendre(n_ordinates=8)),
        interior=N * ng * 10,
        trace=2 * N * ng,  # xmin + xmax, each a point (face shape ())
    )


@pytest.mark.rests_on(
    "tests/gates/transport/test_phase_spaces.py::test_the_face_layout_is_the_hand_layout",
    "tests/gates/mesh/test_trace_face_shapes.py::test_the_traced_faces_are_the_hand_inventory",
)
def test_f0_8_full_field_size_sphere() -> None:
    """Solid sphere, 20 shells, GL8, 2 groups: one traced face (the outer surface; the pole is not a face)."""
    axes = (RadialAxisMesh(edges=np.linspace(0.0, 1.0, 21), coord=AxisCoord.RADIAL_SPHERICAL, role_outer=0),)
    N, ng = 8, 2
    _require_sizes(
        _pose(axes, Quadrature.gauss_legendre(n_ordinates=8)),
        interior=N * ng * 20,
        trace=N * ng,  # xmax only, a point
    )


@pytest.mark.rests_on(
    "tests/gates/transport/test_phase_spaces.py::test_the_face_layout_is_the_hand_layout",
    "tests/gates/mesh/test_trace_face_shapes.py::test_the_traced_faces_are_the_hand_inventory",
)
def test_f0_8_full_field_size_2d_cartesian() -> None:
    """2-D Cartesian 5 x 7, level-symmetric S4, 2 groups: x faces span the 7 y-cells, y faces the 5 x-cells."""
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),  # n=5
        AxisMesh(edges=np.linspace(0.0, 1.0, 8), role_low=_REFLECTIVE, role_high=_REFLECTIVE),  # n=7
    )
    N, ng = 4 * (4 + 2), 2  # a level-symmetric S_n rule has n(n+2) ordinates over the sphere
    _require_sizes(
        _pose(axes, _level_symmetric_quad_2d(order=4)),
        interior=N * ng * 5 * 7,
        trace=2 * N * ng * 7 + 2 * N * ng * 5,  # xmin + xmax over y, ymin + ymax over x
    )


@pytest.mark.rests_on(
    "tests/gates/transport/test_phase_spaces.py::test_the_face_layout_is_the_hand_layout",
    "tests/gates/mesh/test_trace_face_shapes.py::test_the_traced_faces_are_the_hand_inventory",
)
def test_f0_8_full_field_size_3d_cartesian() -> None:
    """3-D Cartesian 5 x 7 x 9, level-symmetric S4, 2 groups: six faces, each spanning the other two axes."""
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),   # n=5
        AxisMesh(edges=np.linspace(0.0, 1.0, 8), role_low=_REFLECTIVE, role_high=_REFLECTIVE),   # n=7
        AxisMesh(edges=np.linspace(0.0, 1.0, 10), role_low=_REFLECTIVE, role_high=_REFLECTIVE),  # n=9
    )
    N, ng = 4 * (4 + 2), 2
    _require_sizes(
        _pose(axes, _level_symmetric_quad_2d(order=4)),
        interior=N * ng * 5 * 7 * 9,
        trace=2 * N * ng * (7 * 9 + 5 * 9 + 5 * 7),  # x faces over (y, z), y faces over (x, z), z faces over (x, y)
    )
