r"""C5.1 axis-primary inversion gates (#225).

``SNProblem.from_axes`` stores the caller's axes VERBATIM — the pre-C5.1
implementation synthesized a legacy ``Mesh1D`` / ``Mesh2D`` from the
axes and the constructor re-derived the axes from it (an
axes → mesh → axes round-trip that silently reset custom endpoint
labels). Both construction surfaces now funnel into ONE body
(``SNProblem._init_core``) whose shape metadata derives from the axes.

The equivalence gates pin the d≤2 invariant the inversion must
preserve: an SNProblem built from an axis tuple carries metadata
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

from orpheus.data.macro_xs.synthetic import make_mixture
from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.mesh import CellsByCount, Mesh1D, Mesh2D, Mesher
from orpheus.numerics.quadrature import Quadrature
from orpheus.mesh import AxisCoord, AxisMesh, RadialAxisMesh, coord_system
from orpheus.sn.problem import SNProblem
from orpheus.geometry.deck import MirrorDeck
from orpheus.transport.boundary_conditions import WallCondition
from tests.gates._exterior import conditions_for

#: Exterior conditions of the tags this file's meshes carry (wave 2: a face's role is a deck or an int tag; the condition a tag
#: stands for is declared at the solve, keyword ``boundary_conditions``).
_CONDITIONS = {0: WallCondition.vacuum()}


pytestmark = [pytest.mark.foundation]

# The law every axis here declares where its subject is not the boundary: every
# consumer resolved an undeclared face to it before step 3c (#405).
_REFLECTIVE = MirrorDeck()

#: The reflective role of a CURVED body's outer face: a tag (rule R2), declared as the unit specular wall by the table ``_CURVED``.
_REFLECTIVE_CURVED = 1
_CURVED = {1: WallCondition.specular(1.0)}


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

    legacy = SNProblem(
        (_carrier1 := Mesh2D(
            edges_x=edges_x, edges_y=edges_y, mat_map=mat_map,
            coord=CoordSystem.CARTESIAN,
            face_roles={
                "xmin": 0,
                "xmax": 0,
                "ymin": MirrorDeck(),
                "ymax": 0,
            },
        )),
        quad, _MATERIALS, boundary_conditions=conditions_for(_carrier1, _CONDITIONS),
    )
    native = SNProblem.from_axes(
        (_carrier2 := (
            AxisMesh(edges=edges_x, role_low=0, role_high=0),
            AxisMesh(edges=edges_y, role_low=MirrorDeck(), role_high=0),
        )),
        quad, _MATERIALS, mat_map=mat_map, boundary_conditions=conditions_for(_carrier2, _CONDITIONS),
    )

    np.testing.assert_equal(native.spatial_shape, legacy.spatial_shape)
    np.testing.assert_equal(native.coord, legacy.coord)
    np.testing.assert_equal(native.nx, legacy.nx)
    for axis in range(2):
        np.testing.assert_array_equal(
            native.axis_widths[axis], legacy.axis_widths[axis],
        )
        np.testing.assert_array_equal(
            native.streaming(axis), legacy.streaming(axis),
        )
    np.testing.assert_array_equal(native.volumes, legacy.volumes)
    np.testing.assert_array_equal(native.mat_map, legacy.mat_map)
    np.testing.assert_equal(sorted(native.bc), sorted(legacy.bc))
    for face in legacy.bc:
        np.testing.assert_equal(type(native.bc[face]), type(legacy.bc[face]))
    np.testing.assert_equal(
        native.boundary_face_layout.faces.keys(),
        legacy.boundary_face_layout.faces.keys(),
    )


def test_1d_slab_metadata_byte_identical_axis_vs_legacy() -> None:
    """The SAME slab via from_axes and via Mesh1D carries identical metadata."""
    edges = np.linspace(0.0, 4.0, 9)
    quad = Quadrature.gauss_legendre(n_ordinates=8)

    legacy_mesh = Mesher(StructuredGeometry.slab(
        (0.0, 4.0), (0,), left=0, right=MirrorDeck(),
    )).partition(CellsByCount.uniform_width(8)).mesh
    legacy = SNProblem(
        legacy_mesh, quad, _MATERIALS, boundary_conditions=conditions_for(legacy_mesh, _CONDITIONS),
    )
    native = SNProblem.from_axes(
        (_carrier3 := (AxisMesh(edges=edges, role_low=0, role_high=MirrorDeck()),)),
        quad, _MATERIALS, boundary_conditions=conditions_for(_carrier3, _CONDITIONS),
    )

    np.testing.assert_equal(native.spatial_shape, legacy.spatial_shape)
    np.testing.assert_array_equal(native.axis_widths[0], legacy.axis_widths[0])
    np.testing.assert_array_equal(native.volumes, legacy.volumes)
    np.testing.assert_array_equal(native.areas, legacy.areas)
    np.testing.assert_array_equal(native.mat_map, legacy.mat_map)
    np.testing.assert_equal(sorted(native.bc), ["xmax", "xmin"])
    np.testing.assert_equal(native.bc["xmin"].is_vacuum, True)
    np.testing.assert_equal(
        (native.bc["xmax"].permutes_ordinates, native.bc["xmax"].return_fraction), (True, 1.0),
    )  # the mirror deck: a permutation returning everything


# ─── Verbatim axes (the round-trip is gone) ──────────────────────────────


def test_from_axes_stores_axes_verbatim() -> None:
    """The supplied axis OBJECTS are the mesh's axes — no re-derivation."""
    axes = (
        AxisMesh(edges=np.linspace(0.0, 1.0, 5), role_low=_REFLECTIVE, role_high=_REFLECTIVE),
        AxisMesh(edges=np.linspace(0.0, 1.0, 4), role_low=_REFLECTIVE, role_high=_REFLECTIVE),
    )
    mesh = SNProblem.from_axes(
        axes, Quadrature.lebedev(17), _MATERIALS, boundary_conditions=conditions_for(axes, _CONDITIONS),
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
            edges=np.linspace(0.0, 1.0, 5), role_low=_REFLECTIVE, role_high=_REFLECTIVE,
            label_low="left", label_high="right",
        ),
    )
    with pytest.raises(ValueError, match="face name"):
        SNProblem.from_axes(
            axes, Quadrature.gauss_legendre(n_ordinates=4), _MATERIALS, boundary_conditions=conditions_for(axes, _CONDITIONS),
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
        lambda: Quadrature.folded_product(n_mu=2, n_phi=4),
    ),
])
def test_from_axes_curvilinear_keeps_mesh1d_reduced_path(
    coord, sys, make_quad,
) -> None:
    """Axis-native curvilinear still feeds the Mesh1D-bound reduced operators.

    The 1-D reduced streaming constructors
    (:func:`orpheus.sn.mesh.reduced_operator.spherical_streaming` /
    ``cylindrical_streaming``) are the genuine remaining ``Mesh1D``
    consumers — ``from_axes`` synthesizes the adapter for them and the
    result must match the legacy surface.
    """
    edges = np.linspace(0.0, 1.0, 6)
    quad = make_quad()
    native = SNProblem.from_axes(
        (_carrier4 := (RadialAxisMesh(edges=edges, coord=coord, role_outer=0),)),
        quad, _MATERIALS, boundary_conditions=conditions_for(_carrier4, _CONDITIONS),
    )
    legacy_mesh = Mesher(StructuredGeometry(
        coord=sys, breakpoints=(0.0, 1.0), mat_ids=(0,), face_roles=(0,),
    )).partition(CellsByCount.uniform_width(5)).mesh
    legacy = SNProblem(
        legacy_mesh, quad, _MATERIALS, boundary_conditions=conditions_for(legacy_mesh, _CONDITIONS),
    )
    np.testing.assert_equal(native.coord, sys)
    np.testing.assert_equal(native.coord, legacy.coord)
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


# ─── C5.5 — d=3 axis-native admission (the constructive gates) ──────────


def _d3_axes() -> tuple:
    """nx=3, ny=4, nz=5 with distinct extents — Mode-2 asymmetry."""
    return (
        AxisMesh(edges=np.linspace(0.0, 1.0, 4), role_low=_REFLECTIVE, role_high=_REFLECTIVE),
        AxisMesh(edges=np.linspace(0.0, 2.0, 5), role_low=_REFLECTIVE, role_high=_REFLECTIVE),
        AxisMesh(edges=np.linspace(0.0, 3.0, 6), role_low=_REFLECTIVE, role_high=_REFLECTIVE),
    )


def test_from_axes_d3_constructs() -> None:
    """C5-G1 (the red-first sentinel, green since C5.5): 3-axis admission.

    A 3-axis Cartesian SNProblem constructs through the SAME generic body
    as d≤2 — mesh-adapter-free from birth (``mesh is None``).
    """
    m = SNProblem.from_axes(
        (_carrier5 := _d3_axes()), Quadrature.level_symmetric(sn_order=4), _MATERIALS, boundary_conditions=conditions_for(_carrier5, _CONDITIONS),
    )
    np.testing.assert_equal(m.ndim, 3)
    np.testing.assert_equal(m.spatial_shape, (3, 4, 5))
    np.testing.assert_equal(m.mesh, None)
    np.testing.assert_equal(m.is_cartesian, True)
    np.testing.assert_equal(m.is_1d, False)
    np.testing.assert_equal(m.nx, 3)


def test_d3_cartesian_builds_trace_and_bc_inventory() -> None:
    """C5-G7: six faces on trace, layout, and the bc dict — one inventory."""
    m = SNProblem.from_axes(
        (_carrier6 := _d3_axes()), Quadrature.level_symmetric(sn_order=4), _MATERIALS, boundary_conditions=conditions_for(_carrier6, _CONDITIONS),
    )
    six = ("xmin", "xmax", "ymin", "ymax", "zmin", "zmax")
    np.testing.assert_equal(tuple(m.angular_trace.layout.faces), six)
    np.testing.assert_equal(tuple(m.boundary_face_layout.faces), six)
    np.testing.assert_equal(sorted(m.bc), sorted(six))
    np.testing.assert_equal(
        tuple(label.face_name for label in m.face_labels), six,
    )


def test_d3_volumes_cartesian_outer_product() -> None:
    """C5-G5: axis-native volumes equal the hand outer product of widths.

    The oracle recomputes from the EDGES (structurally independent of
    the implementation's reduce-outer call): V[i,j,k] = Δx_i·Δy_j·Δz_k.
    """
    axes = _d3_axes()
    m = SNProblem.from_axes(
        axes, Quadrature.level_symmetric(sn_order=4), _MATERIALS, boundary_conditions=conditions_for(axes, _CONDITIONS),
    )
    dx, dy, dz = (np.diff(ax.edges) for ax in axes)
    oracle = np.einsum("i,j,k->ijk", dx, dy, dz)
    np.testing.assert_equal(m.volumes.shape, (3, 4, 5))
    np.testing.assert_allclose(m.volumes, oracle, rtol=1e-15)
    for axis in range(3):
        np.testing.assert_array_equal(
            m.axis_widths[axis], np.diff(axes[axis].edges),
        )


def test_volume_measure_d3_integrates_to_total_volume() -> None:
    """C5-G12: ∫ 1 dV per group equals the box volume Lx·Ly·Lz (closed form)."""
    m = SNProblem.from_axes(
        (_carrier7 := _d3_axes()), Quadrature.level_symmetric(sn_order=4), _MATERIALS, boundary_conditions=conditions_for(_carrier7, _CONDITIONS),
    )
    ones = np.ones((int(np.prod(m.spatial_shape)), m.ng))
    total = m.volume_measure(ones)
    np.testing.assert_allclose(total, [1.0 * 2.0 * 3.0] * m.ng, rtol=1e-13)


# ─── coord_system pure primitive ─────────────────────────────────────────


def test_coord_system_primitive() -> None:
    """Single-axis maps; multi-axis all-Cartesian; mixed refuses."""
    cart = AxisMesh(edges=np.linspace(0.0, 1.0, 3), role_low=_REFLECTIVE, role_high=_REFLECTIVE)
    sph = RadialAxisMesh(
        edges=np.linspace(0.0, 1.0, 3), coord=AxisCoord.RADIAL_SPHERICAL, role_outer=_REFLECTIVE_CURVED,
    )
    cyl = RadialAxisMesh(
        edges=np.linspace(0.0, 1.0, 3), coord=AxisCoord.RADIAL_CYLINDRICAL, role_outer=_REFLECTIVE_CURVED,
    )
    np.testing.assert_equal(coord_system((cart,)), CoordSystem.CARTESIAN)
    np.testing.assert_equal(coord_system((sph,)), CoordSystem.SPHERICAL)
    np.testing.assert_equal(coord_system((cyl,)), CoordSystem.CYLINDRICAL)
    np.testing.assert_equal(
        coord_system((cart, cart, cart)), CoordSystem.CARTESIAN,
    )
    with pytest.raises(NotImplementedError, match="all-Cartesian"):
        coord_system((cart, sph))


# ─── C5.2 — phantom-shim retirement + native volume_measure ─────────────


def _slab_sn() -> SNProblem:
    return SNProblem.from_axes(
        (_carrier8 := (AxisMesh(edges=np.linspace(0.0, 4.0, 9), role_low=_REFLECTIVE, role_high=_REFLECTIVE),)),
        Quadrature.gauss_legendre(n_ordinates=8), _MATERIALS, boundary_conditions=conditions_for(_carrier8, _CONDITIONS),
    )


def test_c52_retired_shims_fail_loud() -> None:
    """``ny``/``dx``/``dy`` raise AttributeError (C5.2 retirement, #225).

    ``ny``/``dy`` LIED at d=1 (phantom ``1``/``[1.0]`` — the #214 bug
    class) and underspecify at d≥3; ``dx`` was a duplicate spelling of
    ``axis_widths[0]``. ``nx`` survives as ``spatial_shape[0]`` sugar.
    """
    sn = _slab_sn()
    np.testing.assert_equal(sn.nx, sn.spatial_shape[0])
    for retired in ("ny", "dx", "dy"):
        with pytest.raises(AttributeError):
            getattr(sn, retired)
    mat_xs = sn.mat_xs
    np.testing.assert_equal(mat_xs.spatial_shape, sn.spatial_shape)
    for retired in ("nx", "ny"):
        with pytest.raises(AttributeError):
            getattr(mat_xs, retired)


def test_volume_measure_d2_delegates_byte_identical() -> None:
    """C5-G13: ``SNProblem.volume_measure`` ≡ the legacy dataclass's measure.

    The SN-side consumers (keff production/absorption rates) now read
    ``problem.volume_measure``; while the mesh adapter is present the
    property must integrate byte-identically to
    ``problem.mesh.volume_measure`` (same atoms, same construction).
    """
    edges_x = np.linspace(0.0, 2.0, 5)
    edges_y = np.linspace(0.0, 3.0, 8)
    sn = SNProblem.from_axes(
        (_carrier9 := (
            AxisMesh(edges=edges_x, role_low=_REFLECTIVE, role_high=_REFLECTIVE),
            AxisMesh(edges=edges_y, role_low=_REFLECTIVE, role_high=_REFLECTIVE),
        )),
        Quadrature.lebedev(17), _MATERIALS, boundary_conditions=conditions_for(_carrier9, _CONDITIONS),
    )
    rng = np.random.default_rng(42)
    values = rng.random((int(np.prod(sn.spatial_shape)), sn.ng))
    np.testing.assert_array_equal(
        sn.volume_measure(values), sn.mesh.volume_measure(values),
    )


# ─── C5.3 — geometry-blind trace; 2-D cylindrical refused at construction ─


def test_2d_cylindrical_mesh_refused_at_construction() -> None:
    """A 2-D cylindrical Mesh2D cannot become an SNProblem (C5.3 pin, #225).

    Migrated from tests/gates/numerics/test_trace_space.py
    ``test_2d_cylindrical_raises``: the trace space is geometry-blind
    now (it never sees a mesh), so the refusal lives where the geometry
    enters — the axis conversion at SNProblem construction.
    """
    mesh2d_cyl = Mesh2D(
        edges_x=np.linspace(0.0, 1.0, 4),
        edges_y=np.linspace(0.0, 1.0, 4),
        mat_map=np.zeros((3, 3), dtype=int),
        coord=CoordSystem.CYLINDRICAL,
        face_roles={
            "xmax": _REFLECTIVE_CURVED,
            "ymin": MirrorDeck(),
            "ymax": MirrorDeck(),
        },
    )
    with pytest.raises(NotImplementedError, match="non-Cartesian"):
        SNProblem(mesh2d_cyl, Quadrature.lebedev(17), _MATERIALS, boundary_conditions=conditions_for(mesh2d_cyl, _CURVED))


def test_d2_trace_builds_for_every_constructible_geometry() -> None:
    """C5-G8: trace is UNCONDITIONAL — every constructible SNProblem has one."""
    slab = _slab_sn()
    np.testing.assert_equal(sorted(slab.angular_trace.layout.faces), ["xmax", "xmin"])
    sphere = SNProblem.from_axes(
        (_carrier10 := (RadialAxisMesh(
            edges=np.linspace(0.0, 1.0, 4),
            coord=AxisCoord.RADIAL_SPHERICAL, role_outer=_REFLECTIVE_CURVED,
        ),)),
        Quadrature.gauss_legendre(n_ordinates=8), _MATERIALS, boundary_conditions=conditions_for(_carrier10, _CURVED),
    )
    np.testing.assert_equal(tuple(sphere.angular_trace.layout.faces), ("xmax",))
