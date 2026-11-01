r"""The transport phase spaces built from their ingredients (``orpheus.transport.phase_spaces``, layer DAG Phase 4).

Each constructor is checked against a REFERENCE assembled here by a different route from the same verified
primitives (``instrument-doctrine`` X4: the numerics space classes are trusted below the line and gated in
``tests/gates/numerics/``; what each constructor OWNS is the composition, and that is what the reference spells by
hand):

* ``angular_bulk_space(mesh, quadrature)``: axis 0 is the quadrature's ordinate axis (its weights are the rule's
  ``w_n``), axes 1.. are the material mesh's own ``bulk_space`` axes, the SAME objects; its metric diagonal is
  ``w_n x 1 x V_cell`` with ``V_cell`` computed here from the cell edges (slab ``h``, cylinder ``pi (r_o^2 - r_i^2)``,
  sphere ``4 pi / 3 (r_o^3 - r_i^3)``, rectangle ``dx dy``), not read from the mesh;
* ``angular_face_layout(mesh, quadrature, scheme)``: the hand layout ``(face, (N, ng, *face shape, *tail))`` with the
  tail ``(2,)`` exactly on a 2-D LD face (``p^{d-1} = 2``) and absent elsewhere (a 1-D face is a point);
* ``angular_trace_space(...)``: equal to ``AngularTraceSpace.from_quadrature_and_layout(quadrature, <hand layout>)``;
* ``scalar_trace_space(mesh)``: equal to ``ScalarTraceSpace.for_faces(<hand faces>, ng, <hand areas>)``, with the
  areas from the geometry (slab 1, cylinder ``2 pi R``, sphere ``4 pi R^2``; a hollow body's inner surface carries no
  trace, #511); a 2-D mesh is refused by the constructor itself (SCOPE-BOUNDARY, #685).

The constructors are raw: two calls give EQUAL spaces with equal hashes; ``is`` is not asserted (a problem's cached
property is where one instance per problem lives, gated in ``tests/gates/sn/mesh/test_angular_bulk_space.py`` and
``tests/gates/diffusion/test_augmented_mesh.py``).

Population: slab, solid and hollow sphere and cylinder, a non-square rectangle (``nx != ny``), DD and LD where the
scheme admits the chart, and the four two-region monomorphic meshes (non-uniform cells). Named first reds, one
in-process arm each (``scratch/layer_dag/phase4/battery/p4_arms.py``): ``bulk_uniform_ordinate_weights``,
``bulk_copies_the_mesh_axes``, ``layout_drops_the_moment_tail``, ``layout_group_before_ordinate``,
``scalar_area_keyed_on_max``, ``scalar_area_constant``.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from orpheus.numerics.face_layout import FaceLayout
from orpheus.numerics.spaces.angular_trace_space import AngularTraceSpace
from orpheus.numerics.spaces.scalar_trace_space import ScalarTraceSpace
from orpheus.transport.phase_spaces import angular_bulk_space, angular_face_layout, angular_trace_space, scalar_trace_space
from orpheus.transport.spatial import DiamondDifference, LinearDiscontinuous
from tests.gates.transport import _phase_space_cases as cases

pytestmark = pytest.mark.foundation

_CURVED = ("sphere", "hollow_sphere", "cylinder", "hollow_cylinder")


def _hand_volumes(name: str, mesh) -> np.ndarray:
    """Cell volumes from the cell edges and the chart, never from ``mesh.volumes``."""
    if name == "cart2d":
        dx = np.diff(mesh.mesh.edges_x)
        dy = np.diff(mesh.mesh.edges_y)
        return np.outer(dx, dy)
    e = np.asarray(mesh.mesh.edges, float)
    if "sphere" in name:
        return 4.0 * math.pi / 3.0 * np.diff(e ** 3)
    if "cylinder" in name:
        return math.pi * np.diff(e ** 2)
    return np.diff(e)


#: The hand traced faces: (face name, its spatial shape), in the canonical order.
_HAND_FACES = {
    "slab": (("xmin", ()), ("xmax", ())),
    "sphere": (("xmax", ()),), "hollow_sphere": (("xmax", ()),),
    "cylinder": (("xmax", ()),), "hollow_cylinder": (("xmax", ()),),
    "cart2d": (("xmin", (2,)), ("xmax", (2,)), ("ymin", (3,)), ("ymax", (3,))),
}
#: The hand outer-surface areas (R = 1): slab 1 per face, cylinder 2 pi R, sphere 4 pi R^2.
_HAND_AREAS = {"slab": 1.0, "sphere": 4.0 * math.pi, "hollow_sphere": 4.0 * math.pi,
               "cylinder": 2.0 * math.pi, "hollow_cylinder": 2.0 * math.pi}


def _schemes(name: str):
    """DD on every chart; LD where it has a moment mass (a slab and a rectangle; a radial chart refuses it)."""
    return [DiamondDifference()] + ([] if name in _CURVED else [LinearDiscontinuous()])


def _hand_layout(name: str, mesh, quadrature, scheme) -> FaceLayout:
    tail = (2,) if (isinstance(scheme, LinearDiscontinuous) and mesh.ndim == 2) else ()
    return FaceLayout.from_named_shapes(
        [(face, (quadrature.n_ordinates, mesh.ng, *shape, *tail)) for face, shape in _HAND_FACES[name]]
    )


_CASES = [(g, type(s).__name__) for g in cases.GEOMETRIES for s in _schemes(g)]


def _case(geometry: str, scheme_name: str):
    mesh, quadrature = cases.GEOMETRIES[geometry](), cases.QUADRATURES[geometry]()
    scheme = {"DiamondDifference": DiamondDifference, "LinearDiscontinuous": LinearDiscontinuous}[scheme_name]()
    return mesh, quadrature, scheme


# ── angular_bulk_space ──


@pytest.mark.parametrize("geometry", list(cases.GEOMETRIES))
def test_the_angular_bulk_is_the_ordinate_axis_on_the_meshs_own_axes(geometry: str) -> None:
    mesh, quadrature = cases.GEOMETRIES[geometry](), cases.QUADRATURES[geometry]()
    space = angular_bulk_space(mesh, quadrature)
    axes, bulk_axes = space.axes, mesh.bulk_space.axes
    if axes is None or bulk_axes is None or len(axes) != 1 + len(bulk_axes):
        pytest.fail(f"{geometry}: the angular bulk is not one axis plus the material mesh's bulk axes")
    for k, (a, b) in enumerate(zip(axes[1:], bulk_axes)):
        if a is not b:
            pytest.fail(f"{geometry}: axis {k + 1} is not the material mesh's own bulk axis object")
    head = axes[0]
    if head.label != "angular" or head.shape != (quadrature.n_ordinates,):
        pytest.fail(f"{geometry}: the leading axis is {head.label!r} {head.shape}")
    np.testing.assert_array_equal(np.asarray(head.weights, float), np.asarray(quadrature.weights, float))
    if space.shape != (quadrature.n_ordinates, mesh.ng, *mesh.spatial_shape):
        pytest.fail(f"{geometry}: shape {space.shape}")


@pytest.mark.parametrize("geometry", list(cases.GEOMETRIES))
def test_the_angular_bulk_metric_is_the_ordinate_weight_times_the_cell_volume(geometry: str) -> None:
    mesh, quadrature = cases.GEOMETRIES[geometry](), cases.QUADRATURES[geometry]()
    diag = angular_bulk_space(mesh, quadrature).metric_diagonal()
    w = np.asarray(quadrature.weights, float)
    V = _hand_volumes(geometry, mesh)
    expected = w.reshape((-1,) + (1,) * (1 + V.ndim)) * np.ones((1, mesh.ng) + (1,) * V.ndim) * V[None, None]
    np.testing.assert_allclose(np.asarray(diag).reshape(expected.shape), expected, rtol=8 * np.finfo(float).eps, atol=0)


@pytest.mark.parametrize("name", ["slab", "sphere", "cylinder", "cart2d"])
def test_the_angular_bulk_on_the_monomorphic_meshes(name: str) -> None:
    """The two-region, non-uniform meshes (the sphere's volumes span three decades): the same two laws."""
    mesh, quadrature = cases.monomorphic()[name]
    space = angular_bulk_space(mesh, quadrature)
    axes, bulk_axes = space.axes, mesh.bulk_space.axes
    if axes is None or bulk_axes is None:
        pytest.fail(f"{name}: a space is not axis-built")
    for a, b in zip(axes[1:], bulk_axes):
        if a is not b:
            pytest.fail(f"{name}: a bulk axis is not the mesh's own object")
    w = np.asarray(quadrature.weights, float)
    V = np.asarray(mesh.volumes, float)  # the edge-derived reference is the GEOMETRY rows' job; here the cells vary
    expected = w.reshape((-1,) + (1,) * (1 + V.ndim)) * np.ones((1, mesh.ng) + (1,) * V.ndim) * V[None, None]
    np.testing.assert_allclose(np.asarray(space.metric_diagonal()).reshape(expected.shape), expected, rtol=8 * np.finfo(float).eps, atol=0)


@pytest.mark.parametrize("geometry", list(cases.GEOMETRIES))
def test_two_calls_build_equal_bulks(geometry: str) -> None:
    mesh, quadrature = cases.GEOMETRIES[geometry](), cases.QUADRATURES[geometry]()
    a, b = angular_bulk_space(mesh, quadrature), angular_bulk_space(mesh, quadrature)
    if a != b or hash(a) != hash(b):
        pytest.fail(f"{geometry}: two calls on one mesh build unequal bulks")


def test_a_different_mesh_builds_a_different_bulk() -> None:
    """The content discriminator: equal shapes, unequal cell volumes, unequal spaces."""
    q = cases.QUADRATURES["slab"]()
    if angular_bulk_space(cases.slab(), q) == angular_bulk_space(cases.stretched_slab(), q):
        pytest.fail("a slab and a slab of twice the width build equal angular bulks")


# ── angular_face_layout and angular_trace_space ──


@pytest.mark.parametrize("geometry,scheme_name", _CASES, ids=[f"{g}-{s}" for g, s in _CASES])
@pytest.mark.rests_on("tests/gates/mesh/test_trace_face_shapes.py::test_the_traced_faces_are_the_hand_inventory")
def test_the_face_layout_is_the_hand_layout(geometry: str, scheme_name: str) -> None:
    mesh, quadrature, scheme = _case(geometry, scheme_name)
    got = angular_face_layout(mesh, quadrature, scheme)
    expected = _hand_layout(geometry, mesh, quadrature, scheme)
    got_slots = [(key, slot.shape) for key, slot in got.faces.items()]
    expected_slots = [(key, slot.shape) for key, slot in expected.faces.items()]
    # per-slot SHAPES in order: ``structure`` holds (key, offset, size) only, blind to a transposed slot of one size
    if got_slots != expected_slots:
        pytest.fail(f"{geometry}/{scheme_name}: layout {got_slots}, the hand layout {expected_slots}")


@pytest.mark.parametrize("geometry,scheme_name", _CASES, ids=[f"{g}-{s}" for g, s in _CASES])
@pytest.mark.rests_on("tests/gates/mesh/test_trace_face_shapes.py::test_the_traced_faces_are_the_hand_inventory")
def test_the_angular_trace_is_the_numerics_trace_over_the_hand_layout(geometry: str, scheme_name: str) -> None:
    mesh, quadrature, scheme = _case(geometry, scheme_name)
    got = angular_trace_space(mesh, quadrature, scheme)
    expected = AngularTraceSpace.from_quadrature_and_layout(quadrature, _hand_layout(geometry, mesh, quadrature, scheme))
    if got != expected or hash(got) != hash(expected):
        pytest.fail(f"{geometry}/{scheme_name}: the angular trace differs from the numerics trace over the hand layout")
    np.testing.assert_array_equal(np.asarray(got.metric_diagonal()), np.asarray(expected.metric_diagonal()))


# ── scalar_trace_space ──


@pytest.mark.parametrize("geometry", [g for g in cases.GEOMETRIES if g != "cart2d"])
@pytest.mark.rests_on("tests/gates/mesh/test_trace_face_shapes.py::test_the_traced_faces_are_the_hand_inventory")
def test_the_scalar_trace_is_the_numerics_trace_over_the_hand_faces_and_areas(geometry: str) -> None:
    mesh = cases.GEOMETRIES[geometry]()
    got = scalar_trace_space(mesh)
    faces = list(_HAND_FACES[geometry])
    expected = ScalarTraceSpace.for_faces(faces, mesh.ng, {face: _HAND_AREAS[geometry] for face, _ in faces})
    if got != expected or hash(got) != hash(expected):
        pytest.fail(f"{geometry}: the scalar trace differs from the numerics trace over the hand faces and areas")
    np.testing.assert_allclose(np.asarray(got.metric_diagonal()), np.asarray(expected.metric_diagonal()), rtol=4 * np.finfo(float).eps, atol=0)


def test_the_scalar_trace_refuses_a_rectangle_naming_its_scope_boundary() -> None:
    """A 2-D mesh's faces carry a per-face-cell surface measure the scalar trace does not hold yet: the constructor
    refuses it by name (SCOPE-BOUNDARY, #685), before any ``areas`` read."""
    with pytest.raises(ValueError, match="scalar_trace_space: a 2-D mesh has a per-face-cell surface measure"):
        scalar_trace_space(cases.cart2d())


# ── Moved from tests/gates/numerics/test_axis_marginal.py by layer DAG Phase 4 (2026-11-01) ──


def test_the_marginal_of_the_angular_bulk_is_the_scalar_bulk() -> None:
    """The angular marginal of the angular bulk IS the material mesh's scalar bulk, not a lookalike: the two share
    the energy and spatial axes verbatim (this module's docstring: "the scalar bulk is this space minus axis 0").
    Its subject is a transport space, so it left the numerics nodal battery, which keeps the retraction's laws."""
    from orpheus.numerics.axis_marginal import AxisRetractionOperator

    mesh = cases.slab(nx=5)
    R = AxisRetractionOperator.of(angular_bulk_space(mesh, cases.QUADRATURES["slab"]()), "angular")
    if R.codomain != mesh.bulk_space:
        pytest.fail(f"the angular marginal must BE the scalar bulk space; got {R.codomain!r} vs {mesh.bulk_space!r}")


# ── The numerics twins match the production shapes (tests/gates/numerics/_phase_space_twins.py's claim) ──


@pytest.mark.parametrize("name", ["slab", "sphere", "cylinder", "cart2d"])
def test_the_numerics_twins_have_the_production_shapes(name: str) -> None:
    """The numerics law files pose their laws on twins of these spaces built from numerics vocabulary only. This row
    holds the twins to the production products on the four monomorphic meshes: the same axis labels, shapes and
    origin kinds, the same ordinate axis, spatial weights equal to 8 eps (a hand volume formula and the mesh's own
    arithmetic differ by about 1e-15 on a slab, so the content is not asserted bit-equal), and the LD trial widening
    the same moment axis. First red: a twin built on another chart's volume formula."""
    from tests.gates.numerics import _phase_space_twins as twins

    mesh, quadrature = cases.monomorphic()[name]
    spatial = (
        twins.spatial_axis_2d(*twins.MONOMORPHIC_RECTANGLE) if name == "cart2d"
        else twins.spatial_axis_1d(name, twins.MONOMORPHIC_EDGES)
    )
    twin = twins.twin(quadrature, spatial)
    production = angular_bulk_space(mesh, quadrature)
    pairs = [(twin.angular_bulk, production)]
    if name in ("slab", "cart2d"):
        pairs.append((twin.ld_trial(), LinearDiscontinuous().trial_space(production, mesh.axes)))
    for got, expected in pairs:
        got_axes, expected_axes = got.axes, expected.axes
        if got_axes is None or expected_axes is None or len(got_axes) != len(expected_axes):
            pytest.fail(f"{name}: the twin's axes do not match production's")
        for a, b in zip(got_axes, expected_axes):
            if (a.label, a.shape, type(a.origin)) != (b.label, b.shape, type(b.origin)):
                pytest.fail(f"{name}: axis {a.label} {a.shape} {type(a.origin).__name__} vs {b.label} {b.shape} {type(b.origin).__name__}")
            if (a.weights is None) != (b.weights is None):
                pytest.fail(f"{name}: axis {a.label}: one side carries weights, the other none")
            if a.weights is not None and b.weights is not None:
                np.testing.assert_allclose(np.asarray(a.weights, float), np.asarray(b.weights, float), rtol=8 * np.finfo(float).eps, atol=0)
    twin_axes, production_axes = twin.angular_bulk.axes, production.axes
    if twin_axes is None or production_axes is None or twin_axes[0] != production_axes[0]:
        pytest.fail(f"{name}: the ordinate axes differ")
