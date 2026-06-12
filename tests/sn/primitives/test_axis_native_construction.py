r"""C5.1 axis-primary inversion gates (#225).

``SNMesh.from_axes`` stores the caller's axes VERBATIM — the pre-C5.1
implementation synthesized a legacy ``Mesh1D`` / ``Mesh2D`` from the
axes and the constructor re-derived the axes from it (an
axes → mesh → axes round-trip that silently reset custom endpoint
labels). Both construction surfaces now funnel into ONE body
(``SNMesh._init_core``) whose shape metadata derives from the axes.

The equivalence gates pin the d≤2 invariant the inversion must
preserve: an SNMesh built from an axis tuple carries metadata
byte-identical (``np.testing.assert_array_equal`` — exact, not
approximate) to one built from the equivalent legacy mesh. The
``np.diff(ax.edges)`` widths spelling is bitwise identical to the
legacy ``Mesh1D.widths`` / ``Mesh2D.dx`` / ``Mesh2D.dy`` because those
are exactly ``np.diff`` of the same edge arrays (mesh.py:287, :567,
:572).

Assertions use ``np.testing.assert_*`` / ``pytest.raises`` exclusively
— the canonical invocation is ``python -O -m pytest``, which strips
bare ``assert`` statements (vv-principles Mode 8).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.derivations.common.xs_library import make_mixture
from orpheus.geometry import BC, CoordSystem, Mesh1D, Mesh2D
from orpheus.numerics.quadrature import Quadrature
from orpheus.sn.axis import AxisCoord, AxisMesh, RadialAxisMesh, coord_system
from orpheus.sn.geometry import SNMesh

pytestmark = [pytest.mark.foundation]


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


_MATERIALS = {0: _one_group_mixture()}


# ─── Equivalence: axis-native ≡ legacy at d≤2 (byte-identical) ───────────


def test_d2_metadata_byte_identical_axis_vs_legacy() -> None:
    """The SAME 2-D mesh via from_axes and via Mesh2D carries identical metadata.

    nx≠ny and asymmetric BCs so an axis swap (Mode 2) is a detectable
    difference, not a symmetric coincidence.
    """
    edges_x = np.linspace(0.0, 2.0, 5)   # nx=4
    edges_y = np.linspace(0.0, 3.0, 8)   # ny=7
    mat_map = np.arange(4 * 7).reshape(4, 7) % 1  # zeros, but shaped
    quad = Quadrature.lebedev(17)

    legacy = SNMesh(
        Mesh2D(
            edges_x=edges_x, edges_y=edges_y, mat_map=mat_map,
            coord=CoordSystem.CARTESIAN,
            bc_xmin=BC("vacuum"), bc_xmax=BC("vacuum"),
            bc_ymin=BC("reflective"), bc_ymax=BC("vacuum"),
        ),
        quad, _MATERIALS,
    )
    native = SNMesh.from_axes(
        (
            AxisMesh(edges=edges_x, bc_low=BC("vacuum"), bc_high=BC("vacuum")),
            AxisMesh(edges=edges_y, bc_low=BC("reflective"), bc_high=BC("vacuum")),
        ),
        quad, _MATERIALS, mat_map=mat_map,
    )

    np.testing.assert_equal(native.spatial_shape, legacy.spatial_shape)
    np.testing.assert_equal(native.coord, legacy.coord)
    np.testing.assert_equal((native.nx, native.ny), (legacy.nx, legacy.ny))
    for axis in range(2):
        np.testing.assert_array_equal(
            native._axis_widths[axis], legacy._axis_widths[axis],
        )
        np.testing.assert_array_equal(
            native.streaming(axis), legacy.streaming(axis),
        )
    np.testing.assert_array_equal(native.volumes, legacy.volumes)
    np.testing.assert_array_equal(native.mat_map, legacy.mat_map)
    np.testing.assert_equal(sorted(native.bc), sorted(legacy.bc))
    for face in legacy.bc:
        np.testing.assert_equal(native.bc[face].kind, legacy.bc[face].kind)
    np.testing.assert_equal(
        native.boundary_face_layout.faces.keys(),
        legacy.boundary_face_layout.faces.keys(),
    )


def test_1d_slab_metadata_byte_identical_axis_vs_legacy() -> None:
    """The SAME slab via from_axes and via Mesh1D carries identical metadata."""
    edges = np.linspace(0.0, 4.0, 9)
    quad = Quadrature.gauss_legendre(n_ordinates=8)

    legacy = SNMesh(
        Mesh1D(
            edges=edges, mat_ids=np.zeros(8, dtype=int),
            coord=CoordSystem.CARTESIAN,
            bc_left=BC("vacuum"), bc_right=BC("reflective"),
        ),
        quad, _MATERIALS,
    )
    native = SNMesh.from_axes(
        (AxisMesh(edges=edges, bc_low=BC("vacuum"), bc_high=BC("reflective")),),
        quad, _MATERIALS,
    )

    np.testing.assert_equal(native.spatial_shape, legacy.spatial_shape)
    np.testing.assert_array_equal(native._axis_widths[0], legacy._axis_widths[0])
    np.testing.assert_array_equal(native.dx, legacy.dx)
    np.testing.assert_array_equal(native.volumes, legacy.volumes)
    np.testing.assert_array_equal(native.areas, legacy.areas)
    np.testing.assert_array_equal(native.mat_map, legacy.mat_map)
    np.testing.assert_equal(sorted(native.bc), ["xmax", "xmin"])
    np.testing.assert_equal(native.bc["xmin"].kind, "vacuum")
    np.testing.assert_equal(native.bc["xmax"].kind, "reflective")


# ─── Verbatim axes (the round-trip is gone) ──────────────────────────────


def test_from_axes_stores_axes_verbatim() -> None:
    """The supplied axis OBJECTS are the mesh's axes — no re-derivation."""
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 5)),
        AxisMesh(edges=np.linspace(0.0, 1.0, 4)),
    )
    mesh = SNMesh.from_axes(
        axes, Quadrature.lebedev(17), _MATERIALS,
    )
    for supplied, stored in zip(axes, mesh.axes):
        if stored is not supplied:
            pytest.fail(
                "from_axes re-derived its axes instead of storing the "
                "caller's objects verbatim (the pre-C5.1 round-trip)."
            )


def test_from_axes_custom_labels_fail_loud() -> None:
    """Custom endpoint labels RAISE at construction (pre-C5.1: silently reset).

    The bc dict is keyed by ``FaceLabel.face_name``, which has no
    rendering for a non-canonical endpoint (C4 doctrine). The
    pre-C5.1 round-trip masked this by silently normalizing the
    labels back to ``min``/``max`` — a lossy rewrite of the caller's
    declaration.
    """
    axes = (
        AxisMesh(
            edges=np.linspace(0.0, 1.0, 5),
            label_low="left", label_high="right",
        ),
    )
    with pytest.raises(ValueError, match="face name"):
        SNMesh.from_axes(
            axes, Quadrature.gauss_legendre(n_ordinates=4), _MATERIALS,
        )


# ─── Curvilinear axis-native path ────────────────────────────────────────


@pytest.mark.parametrize("coord,sys,make_quad", [
    (
        AxisCoord.RADIAL_SPHERICAL,
        CoordSystem.SPHERICAL,
        lambda: Quadrature.gauss_legendre(n_ordinates=8),
    ),
    (
        AxisCoord.RADIAL_CYLINDRICAL,
        CoordSystem.CYLINDRICAL,
        lambda: Quadrature.product(n_mu=2, n_phi=4),
    ),
])
def test_from_axes_curvilinear_keeps_mesh1d_reduced_path(
    coord, sys, make_quad,
) -> None:
    """Axis-native curvilinear still feeds the Mesh1D-bound reduced operators.

    The 1-D reduced streaming constructors
    (:func:`orpheus.geometry.reduced_operator.spherical_streaming` /
    ``cylindrical_streaming``) are the genuine remaining ``Mesh1D``
    consumers — ``from_axes`` synthesizes the adapter for them and the
    result must match the legacy surface.
    """
    edges = np.linspace(0.0, 1.0, 6)
    quad = make_quad()
    native = SNMesh.from_axes(
        (RadialAxisMesh(edges=edges, coord=coord, bc_outer=BC("vacuum")),),
        quad, _MATERIALS,
    )
    legacy = SNMesh(
        Mesh1D(
            edges=edges, mat_ids=np.zeros(5, dtype=int), coord=sys,
            bc_left=None, bc_right=BC("vacuum"),
        ),
        quad, _MATERIALS,
    )
    np.testing.assert_equal(native.coord, sys)
    np.testing.assert_equal(native.curvature, legacy.curvature)
    if native.reduced is None:
        pytest.fail("axis-native curvilinear mesh did not build .reduced")
    if not isinstance(native.mesh, Mesh1D):
        pytest.fail(
            "axis-native curvilinear mesh must synthesize the Mesh1D "
            "adapter the reduced operators consume"
        )
    np.testing.assert_array_equal(native.volumes, legacy.volumes)
    np.testing.assert_array_equal(native.areas, legacy.areas)
    np.testing.assert_equal(sorted(native.bc), sorted(legacy.bc))


# ─── d≥3 admission gate (flips green in C5.5) ────────────────────────────


def test_d3_admission_gate_names_c55() -> None:
    """3-axis construction refuses loudly, naming the admission substep.

    This is the inverse of red-first gate C5-G1
    (``test_from_axes_d3_constructs``): until C5.2–C5.4 dissolve the
    mesh-adapter-bound trace / volume-measure / windowing seams, a
    3-axis SNMesh would be a half-constructed object. This pin is
    DELETED when C5.5 lands and the constructive gate replaces it.
    """
    axes = (AxisMesh(edges=np.linspace(0.0, 1.0, 4)),) * 3
    with pytest.raises(NotImplementedError, match="C5.5"):
        SNMesh.from_axes(
            axes, Quadrature.level_symmetric(sn_order=4), _MATERIALS,
        )


# ─── coord_system pure primitive ─────────────────────────────────────────


def test_coord_system_primitive() -> None:
    """Single-axis maps; multi-axis all-Cartesian; mixed refuses."""
    cart = AxisMesh(edges=np.linspace(0.0, 1.0, 3))
    sph = RadialAxisMesh(
        edges=np.linspace(0.0, 1.0, 3), coord=AxisCoord.RADIAL_SPHERICAL,
    )
    cyl = RadialAxisMesh(
        edges=np.linspace(0.0, 1.0, 3), coord=AxisCoord.RADIAL_CYLINDRICAL,
    )
    np.testing.assert_equal(coord_system((cart,)), CoordSystem.CARTESIAN)
    np.testing.assert_equal(coord_system((sph,)), CoordSystem.SPHERICAL)
    np.testing.assert_equal(coord_system((cyl,)), CoordSystem.CYLINDRICAL)
    np.testing.assert_equal(
        coord_system((cart, cart, cart)), CoordSystem.CARTESIAN,
    )
    with pytest.raises(NotImplementedError, match="all-Cartesian"):
        coord_system((cart, sph))
