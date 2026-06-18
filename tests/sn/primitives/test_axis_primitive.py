r"""Foundation pins for the :class:`Axis1D` primitive (R-1 Phase A C1).

F0.1 — :class:`AxisMesh` endpoints are ``("min", "max")``.
F0.2 — :class:`RadialAxisMesh` endpoints are ``("outer",)``.
F0.3 — :meth:`SNMesh.from_axes` 2-D shape + 4 face labels.
F0.4 — Spherical/cylindrical 1-D mesh has 1 face label (the pole is NOT a face).
F0.5 — :func:`face_shape` derivation matches per-axis cell counts.
F0.6 — :func:`face_outflow_ordinates` mask matches the inline expression at
       the 5+ legacy sites (``np.where(sign * mu_axis > 1e-15)[0]``).
F0.7 — Synthetic 3-D admission: ``(AxisMesh(5), AxisMesh(7), AxisMesh(9))``
       produces 6 face labels via the pure shape functions (no SNMesh).
F0.8 — :attr:`SNMesh.n_unknowns_flat` matches manual computation across
       slab, sphere, 2-D Cartesian, and synthetic 3-D.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC
from orpheus.sn.axis import (
    AxisCoord,
    AxisMesh,
    FaceLabel,
    RadialAxisMesh,
    face_labels,
    face_outflow_ordinates,
    face_shape,
    n_unknowns_flat,
    spatial_shape,
)
from orpheus.sn.geometry import SNMesh
from orpheus.numerics.quadrature import Quadrature


# Every F0.x pin here is a foundation invariant on the dim-agnostic
# axis primitive — it pins a software contract (axis-tuple shape
# algebra, face-label inventory, pack-vector sizing) rather than
# verifying a documented equation, so it carries @foundation and
# does NOT carry @verifies(...). See vv-principles §"V&V level
# taxonomy" — foundation is orthogonal to L0/L1/L2 and pins the
# software invariants that the verification ladder depends on.
pytestmark = [pytest.mark.foundation]


# ─── Fixtures ────────────────────────────────────────────────────────────


def _one_group_mixture():
    """A trivial 1-group, 1-region material used for shape-only tests."""
    return make_mixture(
        sig_t=np.array([1.0]),
        sig_c=np.array([0.5]),
        sig_f=np.array([0.0]),
        nu=np.array([0.0]),
        chi=np.array([1.0]),
        sig_s=np.array([[0.5]]),
    )


def _level_symmetric_quad_2d(order: int = 4):
    """Minimal 2-D quadrature (LevelSymmetric S4)."""
    from orpheus.numerics.quadrature import Quadrature
    return Quadrature.level_symmetric(sn_order=order)


# ─── F0.1 ────────────────────────────────────────────────────────────────


def test_f0_1_axismesh_endpoints_min_max() -> None:
    """AxisMesh's two endpoints are ``min`` and ``max`` (Cartesian convention)."""
    ax = AxisMesh(edges=np.linspace(0.0, 1.0, 16))
    assert ax.endpoints == ("min", "max")
    assert ax.coord == AxisCoord.CARTESIAN
    assert ax.n == 15  # n derived from len(edges) - 1
    assert ax.bc == {"min": None, "max": None}


def test_f0_1_axismesh_endpoints_relabelable() -> None:
    """Endpoint labels can be overridden for slab-conventional naming."""
    ax = AxisMesh(
        edges=np.linspace(0.0, 1.0, 6),
        label_low="left", label_high="right",
        bc_low=BC("vacuum"), bc_high=BC("reflective"),
    )
    assert ax.endpoints == ("left", "right")
    assert ax.bc == {"left": BC("vacuum"), "right": BC("reflective")}


def test_f0_1_axismesh_monotonicity_validated() -> None:
    """Non-monotonic edges raise at construction."""
    with pytest.raises(ValueError, match="monotonically increasing"):
        AxisMesh(edges=np.array([0.0, 1.0, 0.5]))


# ─── F0.2 ────────────────────────────────────────────────────────────────


def test_f0_2_radial_axismesh_one_endpoint_outer() -> None:
    """RadialAxisMesh has exactly one endpoint ``outer`` (pole excluded)."""
    ax = RadialAxisMesh(
        edges=np.linspace(0.0, 2.0, 11),
        coord=AxisCoord.RADIAL_SPHERICAL,
        bc_outer=BC("vacuum"),
    )
    assert ax.endpoints == ("outer",)
    assert ax.n == 10
    assert ax.bc == {"outer": BC("vacuum")}


def test_f0_2_radial_axismesh_label_overridable() -> None:
    """The outer-endpoint label can be customised."""
    ax = RadialAxisMesh(
        edges=np.linspace(0.0, 1.0, 6),
        coord=AxisCoord.RADIAL_CYLINDRICAL,
        bc_outer=BC("vacuum"),
        label_outer="surface",
    )
    assert ax.endpoints == ("surface",)
    assert ax.bc == {"surface": BC("vacuum")}


def test_f0_2_radial_axismesh_rejects_cartesian_coord() -> None:
    """RadialAxisMesh refuses non-radial coords (Pattern 4)."""
    with pytest.raises(ValueError, match="RADIAL_SPHERICAL or RADIAL_CYLINDRICAL"):
        RadialAxisMesh(
            edges=np.linspace(0.0, 1.0, 6),
            coord=AxisCoord.CARTESIAN,
        )


# ─── F0.3 ────────────────────────────────────────────────────────────────


def test_f0_3_from_axes_2d_cartesian_shape_and_face_labels() -> None:
    """SNMesh.from_axes 2-D Cart: shape (nx, ny) + 4 face labels."""
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 11)),  # n=10
        AxisMesh(edges=np.linspace(0.0, 1.0, 8)),   # n=7
    )
    mesh = SNMesh.from_axes(
        axes,
        quadrature=_level_symmetric_quad_2d(order=4),
        materials={0: _one_group_mixture()},
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
    (AxisCoord.RADIAL_CYLINDRICAL, lambda: Quadrature.product(n_mu=2, n_phi=4)),
])
def test_f0_4_solid_radial_mesh_has_one_face_label(coord, make_quad) -> None:
    """Sphere / cylinder SNMesh has exactly 1 face label (``outer``).

    Sphere uses the polar GL quadrature; cylinder requires a level-
    structured product quadrature (LS or Product) because the
    cylindrical streaming reduction depends on level partitioning.
    """
    axes = (
        RadialAxisMesh(
            edges=np.linspace(0.0, 1.0, 21),
            coord=coord,
            bc_outer=BC("vacuum"),
        ),
    )
    mesh = SNMesh.from_axes(
        axes,
        quadrature=make_quad(),
        materials={0: _one_group_mixture()},
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
                     bc_low=BC("vacuum"), bc_high=BC("vacuum")),)
    mesh = SNMesh.from_axes(
        axes,
        quadrature=Quadrature.gauss_legendre(n_ordinates=8),
        materials={0: _one_group_mixture()},
    )
    assert mesh.face_labels == (FaceLabel(0, "min"), FaceLabel(0, "max"))


# ─── F0.5 ────────────────────────────────────────────────────────────────


def test_f0_5_face_shape_2d_each_face() -> None:
    """face_shape: the face on axis-i has shape over the OTHER axis."""
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6)),   # n=5
        AxisMesh(edges=np.linspace(0.0, 1.0, 8)),   # n=7
    )
    assert face_shape(axes, FaceLabel(0, "min")) == (7,)
    assert face_shape(axes, FaceLabel(0, "max")) == (7,)
    assert face_shape(axes, FaceLabel(1, "min")) == (5,)
    assert face_shape(axes, FaceLabel(1, "max")) == (5,)


def test_f0_5_face_shape_3d_synthetic() -> None:
    """3-D synthetic: face on axis-1 has shape over (axis-0, axis-2)."""
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6)),   # n=5
        AxisMesh(edges=np.linspace(0.0, 1.0, 8)),   # n=7
        AxisMesh(edges=np.linspace(0.0, 1.0, 10)),  # n=9
    )
    assert face_shape(axes, FaceLabel(0, "min")) == (7, 9)
    assert face_shape(axes, FaceLabel(1, "min")) == (5, 9)
    assert face_shape(axes, FaceLabel(2, "min")) == (5, 7)


def test_f0_5_face_shape_radial_solid_no_inner_face() -> None:
    """Solid radial mesh has only the outer face; shape is empty tuple."""
    axes = (RadialAxisMesh(
        edges=np.linspace(0.0, 1.0, 11),
        coord=AxisCoord.RADIAL_SPHERICAL,
    ),)
    labels = face_labels(axes)
    assert labels == (FaceLabel(0, "outer"),)
    assert face_shape(axes, labels[0]) == ()


# ─── F0.6 ────────────────────────────────────────────────────────────────


def test_f0_6_face_outflow_ordinates_matches_inline_expression_1d() -> None:
    """Outflow ordinate mask matches np.where(sign * mu > 1e-15)."""
    axes = (AxisMesh(edges=np.linspace(0.0, 1.0, 11)),)
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    out_max = face_outflow_ordinates(axes, FaceLabel(0, "max"), quad)
    out_min = face_outflow_ordinates(axes, FaceLabel(0, "min"), quad)
    np.testing.assert_array_equal(out_max, np.where(quad.mu_x > 1e-15)[0])
    np.testing.assert_array_equal(out_min, np.where(-quad.mu_x > 1e-15)[0])
    # In 1-D GL even-N the two are complementary halves of (0..N-1).
    assert out_max.size == out_min.size == quad.N // 2


def test_f0_6_face_outflow_ordinates_2d_LS4() -> None:
    """2-D LS4: outflow at each face is the half-sphere of that axis."""
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6)),
        AxisMesh(edges=np.linspace(0.0, 1.0, 6)),
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
        coord=AxisCoord.RADIAL_SPHERICAL,
    ),)
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    out_outer = face_outflow_ordinates(axes, FaceLabel(0, "outer"), quad)
    np.testing.assert_array_equal(out_outer, np.where(quad.mu_x > 1e-15)[0])


# ─── F0.7 ────────────────────────────────────────────────────────────────


def test_f0_7_synthetic_3d_admission_six_face_labels() -> None:
    """Synthetic 3-D axis tuple admits the dim-agnostic shape primitives.

    C1 ships the pure shape functions in ``orpheus.sn.axis``; no
    ``Mesh3D`` dataclass exists yet (D9 of the ultraplan). The 3-D
    admission gate exercises these functions on a synthetic axis tuple
    without instantiating an :class:`SNMesh` — that lands in a
    followup once :class:`Mesh3D` is in tree.
    """
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6)),   # n=5
        AxisMesh(edges=np.linspace(0.0, 1.0, 8)),   # n=7
        AxisMesh(edges=np.linspace(0.0, 1.0, 10)),  # n=9
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
        AxisMesh(edges=np.linspace(0.0, 1.0, 6)),
        AxisMesh(edges=np.linspace(0.0, 1.0, 6)),
        AxisMesh(edges=np.linspace(0.0, 1.0, 6)),
    )
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    # GL1D's measure is 1-D scalar — only axis 0 carries real cosines.
    assert quad.measure.nodes.ndim == 1
    out_xmax = face_outflow_ordinates(axes, FaceLabel(0, "max"), quad)
    assert out_xmax.size == quad.N // 2   # axis 0 is the polar μ — real cosines
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
        AxisMesh(edges=np.linspace(0.0, 1.0, 6)),
        AxisMesh(edges=np.linspace(0.0, 1.0, 6)),
        AxisMesh(edges=np.linspace(0.0, 1.0, 6)),
    )
    quad = _level_symmetric_quad_2d(order=4)
    out_zmax = face_outflow_ordinates(axes, FaceLabel(2, "max"), quad)
    out_zmin = face_outflow_ordinates(axes, FaceLabel(2, "min"), quad)
    np.testing.assert_array_equal(out_zmax, np.where(quad.mu_z > 1e-15)[0])
    np.testing.assert_array_equal(out_zmin, np.where(-quad.mu_z > 1e-15)[0])
    # LS4 is symmetric over the unit sphere: half the ordinates have
    # mu_z > 0, half have mu_z < 0.
    assert out_zmax.size == out_zmin.size == quad.N // 2


# ─── F0.8 ────────────────────────────────────────────────────────────────


def test_f0_8_n_unknowns_flat_slab() -> None:
    """Slab 1-D: n_cells + outflow at min + outflow at max."""
    axes = (AxisMesh(edges=np.linspace(0.0, 1.0, 11)),)  # n=10
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    mesh = SNMesh.from_axes(
        axes, quadrature=quad, materials={0: _one_group_mixture()},
    )
    N, ng = quad.N, 1
    n_max = np.where(quad.mu_x > 1e-15)[0].size
    n_min = np.where(-quad.mu_x > 1e-15)[0].size
    expected = N * ng * 10 + (n_max + n_min) * ng * 1  # face_shape=()
    assert mesh.n_unknowns_flat == expected
    assert n_unknowns_flat(axes, quad, ng) == expected


def test_f0_8_n_unknowns_flat_sphere() -> None:
    """Sphere 1-D: n_cells + outflow at outer (single face)."""
    axes = (RadialAxisMesh(
        edges=np.linspace(0.0, 1.0, 21),
        coord=AxisCoord.RADIAL_SPHERICAL,
        bc_outer=BC("vacuum"),
    ),)
    quad = Quadrature.gauss_legendre(n_ordinates=8)
    mesh = SNMesh.from_axes(
        axes, quadrature=quad, materials={0: _one_group_mixture()},
    )
    N, ng = quad.N, 1
    n_outer = np.where(quad.mu_x > 1e-15)[0].size
    expected = N * ng * 20 + n_outer * ng * 1
    assert mesh.n_unknowns_flat == expected


def test_f0_8_n_unknowns_flat_2d_cartesian() -> None:
    """2-D Cart: n_cells + (xmin + xmax) faces over y + (ymin + ymax) over x."""
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6)),  # n=5
        AxisMesh(edges=np.linspace(0.0, 1.0, 8)),  # n=7
    )
    quad = _level_symmetric_quad_2d(order=4)
    mesh = SNMesh.from_axes(
        axes, quadrature=quad, materials={0: _one_group_mixture()},
    )
    N, ng = quad.N, 1
    n_xpos = np.where(quad.mu_x > 1e-15)[0].size
    n_xneg = np.where(-quad.mu_x > 1e-15)[0].size
    n_ypos = np.where(quad.mu_y > 1e-15)[0].size
    n_yneg = np.where(-quad.mu_y > 1e-15)[0].size
    expected = (
        N * ng * 5 * 7                            # cells
        + (n_xneg + n_xpos) * ng * 7              # xmin + xmax over y
        + (n_yneg + n_ypos) * ng * 5              # ymin + ymax over x
    )
    assert mesh.n_unknowns_flat == expected


def test_f0_8_n_unknowns_flat_synthetic_3d() -> None:
    """Synthetic 3-D (no SNMesh): n_unknowns_flat from pure function.

    LS4 already exposes ``mu_z`` as a concrete-class attribute (the
    Protocol declaration catches up in C2), so axis-2 outflow IS
    non-empty for this 3-D admission gate.
    """
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 6)),   # n=5
        AxisMesh(edges=np.linspace(0.0, 1.0, 8)),   # n=7
        AxisMesh(edges=np.linspace(0.0, 1.0, 10)),  # n=9
    )
    quad = _level_symmetric_quad_2d(order=4)
    ng = 1
    N = quad.N
    n_xpos = np.where(quad.mu_x > 1e-15)[0].size
    n_xneg = np.where(-quad.mu_x > 1e-15)[0].size
    n_ypos = np.where(quad.mu_y > 1e-15)[0].size
    n_yneg = np.where(-quad.mu_y > 1e-15)[0].size
    n_zpos = np.where(quad.mu_z > 1e-15)[0].size
    n_zneg = np.where(-quad.mu_z > 1e-15)[0].size
    expected = (
        N * ng * 5 * 7 * 9
        + (n_xneg + n_xpos) * ng * 7 * 9
        + (n_yneg + n_ypos) * ng * 5 * 9
        + (n_zneg + n_zpos) * ng * 5 * 7
    )
    assert n_unknowns_flat(axes, quad, ng) == expected
