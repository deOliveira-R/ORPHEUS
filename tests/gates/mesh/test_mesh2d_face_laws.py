r"""The face laws of both meshes, and the one element parser (P1 step 3c, #405).

``Mesh2D(edges_x, edges_y, mat_map, *, face_laws, coord=CARTESIAN)`` declares one
law per boundary face (the user's ruling of 2026-09-29; the four ``bc_*`` fields
retire), and since the ruling of 2026-09-30 ``Mesh1D`` and ``Mesh2D`` store one
value, :class:`~orpheus.mesh.face_laws.FaceLaws`: a frozen, ordered, picklable
mapping from face name to law over the mesh's
:func:`~orpheus.mesh.face_laws.face_inventory`. The inventory rule: along the
first axis the faces are ``CoordSystem.boundary_points`` of its edges (both ends
on a slab, on (x, y) and on a hollow body; only the outer surface on a solid
radial body, whose centre is interior), along every further axis both ends,
each named by ``FaceLabel(axis_index, endpoint).face_name``. The positional
tuple, ``Mesh1D.boundary_faces`` (now ``boundary_points``) and
``Mesh2D.boundary_faces`` are retired. The named cells carry the model's law
(``TestTheNamedCells``).

Claim kind: THEOREM for every row. The expected inventories are written by hand
from the topology law, independently of ``face_inventory``; the one-parser and
one-inventory rows are ROUTE gates (a counting spy on every binding of
``parse_boundary_law`` / ``face_inventory``), because a value or a message
alone cannot tell the one function from a twin that copied it.

First reds, by class (each named at its row, and each measured in the step-3c
batteries, ``scratch/reference_architecture/p1step3c/battery.md`` and
``battery2.md``):

* inventory: a solid (r, z) mesh given ``("xmin", "xmax")`` along x reds
  ``test_the_inventory``'s solid row and ``TestTheOneInventoryRule``;
* refusals: the pre-3c ``Mesh2D`` (four ``bc_*`` fields defaulting to ``None``)
  reds every row of ``TestRefusals``; ``face_laws`` given a default reds
  ``test_face_laws_is_required_and_keyword_only``;
* storage and value: a ``FaceLaws`` that stores (``__setitem__``), iterates in
  the caller's order, compares by identity or loses its ``__reduce__`` reds the
  matching ``TestFaceLaws`` row and the storage rows;
* the one parser: an inline ``isinstance`` check in any of the five classes in
  place of the ``parse_boundary_law`` call reds that class's row of
  ``TestTheOneElementParser::test_every_declared_law_passes_the_one_parser``.
"""
from __future__ import annotations

import dataclasses
import inspect
import sys
from collections.abc import Callable

import numpy as np
import pytest

import orpheus.geometry.structured_geometry as structured_geometry
from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.geometry.boundary import AlbedoBoundary, ReflectiveBoundary, VacuumInflow
from orpheus.mesh import AxisCoord, AxisMesh, Mesh1D, Mesh2D, RadialAxisMesh
from orpheus.mesh.face_laws import FaceLaws, face_inventory

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/mesh/test_mesh2d_face_laws.py"
_NONE_IN_MESH1D = "tests/gates/mesh/test_mesh1d.py::TestConstructionLaws::test_refusal[none-law]"

_CART = CoordSystem.CARTESIAN
_CYL = CoordSystem.CYLINDRICAL

#: (edges_x, the inventory written by hand from the topology law).
_INVENTORIES = {
    "xy": (_CART, (0.0, 1.0, 2.0), ("xmin", "xmax", "ymin", "ymax")),
    "rz-solid": (_CYL, (0.0, 0.5, 1.0), ("xmax", "ymin", "ymax")),
    "rz-hollow": (_CYL, (0.25, 0.5, 1.0), ("xmin", "xmax", "ymin", "ymax")),
}

_EDGES_Y = (0.0, 1.0, 2.0, 3.0)

#: Four distinct laws, so a permutation of faces cannot read equal.
_DISTINCT = {
    "xmin": BC.vacuum,
    "xmax": BC.reflective,
    "ymin": BC("albedo", {"albedo": 0.3}),
    "ymax": ReflectiveBoundary(axis="y", albedo=1.0),
}


def _laws(faces, source=_DISTINCT) -> dict:
    return {face: source[face] for face in faces}


def _mesh2d(case: str = "xy", **overrides) -> Mesh2D:
    coord, edges_x, faces = _INVENTORIES[case]
    kwargs = {"face_laws": _laws(faces), "coord": coord} | overrides
    mat_map = np.zeros((len(edges_x) - 1, len(_EDGES_Y) - 1), dtype=int)
    return Mesh2D(np.asarray(edges_x), np.asarray(_EDGES_Y), mat_map, **kwargs)


# ─────────────────────────────────────────────────────────────────────
# The face inventory
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.parametrize("case", list(_INVENTORIES))
def test_the_inventory(case):
    """The faces are ``boundary_points`` along x and both ends along y, named by
    ``FaceLabel.face_name``. First red: the solid row, if the x faces stop
    reading ``boundary_points`` (a solid (r, z) given an ``xmin``)."""
    coord, edges_x, expected = _INVENTORIES[case]
    mesh = _mesh2d(case)
    assert isinstance(mesh.face_laws, FaceLaws)
    assert tuple(mesh.face_laws) == expected
    # the premise the hand-written inventory rests on, read from the one law
    n_x = len(coord.boundary_points(edges_x[0], edges_x[-1]))
    assert len(expected) == n_x + 2


def test_the_retired_fields_are_gone():
    """The four ``bc_*`` fields retire, and ``boundary_faces`` with them
    (``tuple(mesh.face_laws)`` are the names; ``retirement-audit`` D.18's
    witness). First red: any of them re-added."""
    names = {f.name for f in dataclasses.fields(Mesh2D)}
    assert names == {"edges_x", "edges_y", "mat_map", "face_laws", "coord"}, names
    for name in ("bc_xmin", "bc_xmax", "bc_ymin", "bc_ymax", "boundary_faces"):
        assert not hasattr(Mesh2D, name)
        assert not hasattr(_mesh2d(), name)


def test_face_laws_is_required_and_keyword_only():
    """First red: a default on ``face_laws`` (``required``), or the field moved
    before ``KW_ONLY`` (``keyword-only``)."""
    param = inspect.signature(Mesh2D).parameters["face_laws"]
    assert param.kind is inspect.Parameter.KEYWORD_ONLY
    assert param.default is inspect.Parameter.empty


# ─────────────────────────────────────────────────────────────────────
# Refusals, one negative each, each keyed to the input that triggers it
# ─────────────────────────────────────────────────────────────────────


class TestRefusals:
    """Each row is red on the pre-3c ``Mesh2D`` (which admitted ``None`` and had
    no ``face_laws``). The inventory is named in every count refusal."""

    def test_omitted(self):
        mat_map = np.zeros((2, 3), dtype=int)
        with pytest.raises(TypeError, match="face_laws"):
            Mesh2D(np.array([0.0, 1.0, 2.0]), np.asarray(_EDGES_Y), mat_map)  # type: ignore[call-arg]  # the refusal input

    def test_positional(self):
        mat_map = np.zeros((2, 3), dtype=int)
        laws = _laws(_INVENTORIES["xy"][2])
        with pytest.raises(TypeError, match="positional argument"):
            Mesh2D(np.array([0.0, 1.0, 2.0]), np.asarray(_EDGES_Y), mat_map, laws)  # type: ignore[misc]  # the refusal input

    @pytest.mark.parametrize("case", list(_INVENTORIES))
    def test_a_missing_face(self, case):
        _, _, faces = _INVENTORIES[case]
        laws = _laws(faces)
        del laws["ymax"]
        with pytest.raises(ValueError, match="Mesh2D.face_laws") as info:
            _mesh2d(case, face_laws=laws)
        for face in faces:  # the message names the inventory
            assert repr(face) in str(info.value), (face, str(info.value))

    @pytest.mark.parametrize("case", list(_INVENTORIES))
    def test_an_extra_face(self, case):
        _, _, faces = _INVENTORIES[case]
        laws = _laws(faces) | {"zmin": BC.vacuum}
        with pytest.raises(ValueError, match="Mesh2D.face_laws") as info:
            _mesh2d(case, face_laws=laws)
        for face in faces:
            assert repr(face) in str(info.value), (face, str(info.value))

    def test_no_law_on_the_inner_surface_of_a_hollow_rz_mesh(self):
        """A hollow (r, z) without ``xmin``: the reason names the inner surface."""
        laws = _laws(("xmax", "ymin", "ymax"))
        with pytest.raises(ValueError, match="has an inner surface, which needs its own law"):
            _mesh2d("rz-hollow", face_laws=laws)

    def test_a_law_on_the_axis_of_a_solid_rz_mesh(self):
        """``xmin`` on a solid (r, z): the axis r = 0 is interior and carries no
        law. The reason is keyed: a hollow mesh missing ``ymax`` does not say it."""
        laws = _laws(("xmin", "xmax", "ymin", "ymax"))
        with pytest.raises(ValueError, match="carries no law") as info:
            _mesh2d("rz-solid", face_laws=laws)
        assert "'xmax', 'ymin', 'ymax'" in str(info.value)
        missing = _laws(_INVENTORIES["rz-hollow"][2])
        del missing["ymax"]
        with pytest.raises(ValueError) as other:
            _mesh2d("rz-hollow", face_laws=missing)
        assert "carries no law" not in str(other.value)

    @pytest.mark.parametrize("face", ["xmin", "ymax"])
    def test_none(self, face):
        laws = _laws(_INVENTORIES["xy"][2]) | {face: None}
        with pytest.raises(TypeError, match="None is not a boundary law") as info:
            _mesh2d(face_laws=laws)
        assert repr(face) in str(info.value)

    @pytest.mark.parametrize("bad", ["vacuum", 0, object()], ids=["string", "int", "object"])
    def test_a_non_law_object(self, bad):
        laws = _laws(_INVENTORIES["xy"][2]) | {"ymin": bad}
        with pytest.raises(TypeError, match="must be a BC tag or a BoundaryTraceLaw") as info:
            _mesh2d(face_laws=laws)
        assert "'ymin'" in str(info.value)

    def test_the_positional_tuple(self):
        """The positional tuple is retired on both meshes (``Mesh1D``'s row is
        ``tests/gates/mesh/test_mesh1d.py`` ``tuple-laws``)."""
        with pytest.raises(TypeError, match="Mesh2D.face_laws is a mapping from face name to law"):
            _mesh2d(face_laws=tuple(_DISTINCT.values()))


# ─────────────────────────────────────────────────────────────────────
# Storage
# ─────────────────────────────────────────────────────────────────────


def test_the_storage_is_read_only():
    """First red: storing the caller's dict (or a copy that is a ``dict``)."""
    caller = _laws(_INVENTORIES["xy"][2])
    mesh = _mesh2d(face_laws=caller)
    with pytest.raises(TypeError):
        mesh.face_laws["xmin"] = BC.reflective  # type: ignore[index]  # the refusal input
    caller["xmin"] = BC.reflective  # the caller's mapping is not the mesh's
    assert mesh.face_laws["xmin"] is _DISTINCT["xmin"]
    with pytest.raises(dataclasses.FrozenInstanceError):
        mesh.face_laws = {}  # type: ignore[misc]  # the refusal input


@pytest.mark.parametrize("case", list(_INVENTORIES))
def test_the_storage_is_in_inventory_order(case):
    """Given in reverse, stored in inventory order, each law the object declared.
    First red: storing in the caller's key order."""
    _, _, faces = _INVENTORIES[case]
    reversed_laws = {face: _DISTINCT[face] for face in reversed(faces)}
    mesh = _mesh2d(case, face_laws=reversed_laws)
    assert tuple(mesh.face_laws) == faces
    for face in faces:
        assert mesh.face_laws[face] is _DISTINCT[face], face


def test_a_typed_law_is_admitted_beside_a_tag():
    """A ``BoundaryTraceLaw`` instance and a ``BC`` tag on one mesh, each stored
    as the object declared. First red: a parser admitting only ``BC``."""
    laws = {"xmin": VacuumInflow(), "xmax": BC.vacuum, "ymin": AlbedoBoundary(albedo=0.5), "ymax": BC.reflective}
    mesh = _mesh2d(face_laws=laws)
    for face, law in laws.items():
        assert mesh.face_laws[face] is law


# ─────────────────────────────────────────────────────────────────────
# The one element parser: Mesh1D, Mesh2D, AxisMesh, RadialAxisMesh and
# StructuredGeometry all declare through ``parse_boundary_law``
# ─────────────────────────────────────────────────────────────────────

_SLAB_EDGES = np.array([0.0, 0.5, 2.0])
_SPH_EDGES = np.array([0.0, 0.5, 2.0])


def _build_mesh1d(low, high):
    return Mesh1D(coord=_CART, edges=_SLAB_EDGES, volumes=_CART.measure(_SLAB_EDGES),
                  mat_ids=np.zeros(2, dtype=int), face_laws={"xmin": low, "xmax": high})


def _build_mesh2d(low, high):
    return _mesh2d(face_laws={"xmin": low, "xmax": high, "ymin": BC.vacuum, "ymax": BC.vacuum})


def _build_axis(low, high):
    return AxisMesh(edges=_SLAB_EDGES, bc_low=low, bc_high=high)


def _build_radial(low, high):
    del high  # one endpoint
    return RadialAxisMesh(edges=_SPH_EDGES, coord=AxisCoord.RADIAL_SPHERICAL, bc_outer=low)


def _build_geometry(low, high):
    return StructuredGeometry(coord=_CART, breakpoints=(0.0, 2.0), mat_ids=(0,), boundaries=(low, high))


#: (owner, the builder, the slot name the refusal names).
_OWNERS: list[tuple[str, Callable[[object, object], object], str]] = [
    ("Mesh1D", _build_mesh1d, "Mesh1D.face_laws['xmin']"),
    ("Mesh2D", _build_mesh2d, "Mesh2D.face_laws['xmin']"),
    ("AxisMesh", _build_axis, "AxisMesh.bc_low"),
    ("RadialAxisMesh", _build_radial, "RadialAxisMesh.bc_outer"),
    ("StructuredGeometry", _build_geometry, "StructuredGeometry.boundaries[0]"),
]


class TestTheOneElementParser:
    """One check of a single boundary declaration (the design of step 3c).

    Two legs: the ROUTE (every declared law passes through
    ``parse_boundary_law``, counted by a spy on every module binding of it) and
    the MESSAGE (``None`` is refused with the parser's fragment, raised in the
    parser's own frame, naming the owner's slot)."""

    @pytest.mark.parametrize("owner, build, slot", _OWNERS, ids=[o[0] for o in _OWNERS])
    def test_every_declared_law_passes_the_one_parser(self, owner, build, slot, monkeypatch):
        del slot
        original = structured_geometry.parse_boundary_law
        seen: list[int] = []

        def spy(law, where):
            seen.append(id(law))
            return original(law, where)

        rebound = 0
        for module in list(sys.modules.values()):
            if getattr(module, "parse_boundary_law", None) is original:
                monkeypatch.setattr(module, "parse_boundary_law", spy)
                rebound += 1
        # the spy is installed where the constructors read it: the defining
        # module, the mesh module and the axis module at least
        assert rebound >= 1 and structured_geometry.parse_boundary_law is spy, rebound
        low, high = BC("albedo", {"albedo": 0.25}), ReflectiveBoundary(axis="x", albedo=1.0)
        build(low, high)
        assert id(low) in seen, f"{owner}: its law never reached parse_boundary_law (rebound {rebound})"
        if owner != "RadialAxisMesh":
            assert id(high) in seen, f"{owner}: its second law never reached parse_boundary_law"

    @pytest.mark.rests_on(_NONE_IN_MESH1D)
    @pytest.mark.parametrize("owner, build, slot", _OWNERS, ids=[o[0] for o in _OWNERS])
    def test_none_is_refused_by_the_one_parser(self, owner, build, slot):
        with pytest.raises(TypeError, match="is not a boundary law") as info:
            build(None, BC.vacuum)
        assert slot in str(info.value), (owner, str(info.value))
        assert info.traceback[-1].name == "parse_boundary_law", info.traceback[-1].name

    @pytest.mark.parametrize("owner, build, slot", _OWNERS, ids=[o[0] for o in _OWNERS])
    def test_a_non_law_is_refused_by_the_one_parser(self, owner, build, slot):
        with pytest.raises(TypeError, match="must be a BC tag or a BoundaryTraceLaw") as info:
            build("vacuum", BC.vacuum)
        assert slot in str(info.value), (owner, str(info.value))
        assert info.traceback[-1].name == "parse_boundary_law", info.traceback[-1].name


# ─────────────────────────────────────────────────────────────────────
# One inventory rule for both meshes (the user's ruling of 2026-09-30:
# one FaceLaws value for Mesh1D and Mesh2D)
# ─────────────────────────────────────────────────────────────────────

#: (coord, first-axis edges, the first-axis faces written by hand).
_FIRST_AXIS = {
    "slab": (_CART, (0.0, 1.0, 2.0), ("xmin", "xmax")),
    "solid-cylinder": (_CYL, (0.0, 0.5, 1.0), ("xmax",)),
    "hollow-cylinder": (_CYL, (0.25, 0.5, 1.0), ("xmin", "xmax")),
}


def _mesh1d(coord, edges, laws):
    e = np.asarray(edges, dtype=float)
    return Mesh1D(coord=coord, edges=e, volumes=coord.measure(e), mat_ids=np.zeros(len(e) - 1, dtype=int),
                  face_laws=laws)


class TestTheOneInventoryRule:
    """``Mesh1D`` and ``Mesh2D`` on the same coordinate system and first-axis
    edges name the same first-axis faces, and a 2-D mesh adds exactly
    ``ymin``, ``ymax``. Two legs: the VALUE (against the hand-written
    inventory) and the ROUTE (both constructors ask ``face_inventory``).
    First red: either mesh spelling its own inventory (a ``Mesh2D`` whose
    solid radial axis is given ``xmin``/``xmax``, the pre-FaceLaws
    ``_faces_2d`` shape, or a constructor bypassing ``face_inventory``)."""

    @pytest.mark.parametrize("case", list(_FIRST_AXIS))
    def test_both_meshes_name_the_same_first_axis_faces(self, case):
        coord, edges, faces = _FIRST_AXIS[case]
        one = _mesh1d(coord, edges, _laws(faces))
        two = Mesh2D(np.asarray(edges), np.asarray(_EDGES_Y), np.zeros((2, 3), dtype=int),
                     face_laws=_laws(faces + ("ymin", "ymax")), coord=coord)
        assert tuple(one.face_laws) == faces
        assert tuple(two.face_laws) == faces + ("ymin", "ymax")
        assert isinstance(one.face_laws, FaceLaws) and isinstance(two.face_laws, FaceLaws)

    def test_both_constructors_ask_the_one_rule(self, monkeypatch):
        import orpheus.mesh.face_laws as face_laws_module

        original = face_laws_module.face_inventory
        dimensions: list[int] = []

        def spy(coord, first_axis_edges, dimension):
            dimensions.append(dimension)
            return original(coord, first_axis_edges, dimension)

        rebound = 0
        for module in list(sys.modules.values()):
            if getattr(module, "face_inventory", None) is original:
                monkeypatch.setattr(module, "face_inventory", spy)
                rebound += 1
        assert rebound >= 2, rebound  # the defining module and the mesh module at least
        _mesh1d(_CART, (0.0, 1.0), _laws(("xmin", "xmax")))
        assert dimensions == [1], dimensions
        _mesh2d()
        assert dimensions == [1, 2], dimensions

    @pytest.mark.parametrize("case", list(_FIRST_AXIS))
    def test_the_rule_itself(self, case):
        coord, edges, faces = _FIRST_AXIS[case]
        assert face_inventory(coord, np.asarray(edges), 1) == faces
        assert face_inventory(coord, np.asarray(edges), 2) == faces + ("ymin", "ymax")
        assert face_inventory(coord, np.asarray(edges), 3) == faces + ("ymin", "ymax", "zmin", "zmax")


# ─────────────────────────────────────────────────────────────────────
# FaceLaws: the value both meshes store
# ─────────────────────────────────────────────────────────────────────

_XY = ("xmin", "xmax", "ymin", "ymax")


class TestFaceLaws:
    """The laws of the type (``vv-testing``: a math-bearing type ships the test
    of its defining laws): inventory order, immutability, mapping equality,
    pickle round trip. First reds: ``__iter__`` over the caller's order (the
    order row); a ``__setattr__`` or ``__setitem__`` that stores (the
    immutability row); an identity ``__eq__`` (the equality row); a
    ``__reduce__`` dropped on a slotted class, or one that rebuilds without the
    laws (the pickle rows)."""

    def test_order_is_the_inventory(self):
        given = {face: _DISTINCT[face] for face in reversed(_XY)}
        laws = FaceLaws.over(_XY, given, "where")
        assert tuple(laws) == _XY
        assert list(laws.items()) == [(face, _DISTINCT[face]) for face in _XY]
        assert len(laws) == 4

    def test_immutability(self):
        laws = FaceLaws.over(_XY, _laws(_XY), "where")
        with pytest.raises(TypeError):
            laws["xmin"] = BC.reflective  # type: ignore[index]  # the refusal input
        with pytest.raises(AttributeError):
            laws._items = ()  # type: ignore[misc]  # the refusal input
        with pytest.raises(AttributeError):
            laws.extra = 1  # type: ignore[attr-defined]  # the refusal input
        assert laws["xmin"] is _DISTINCT["xmin"]
        with pytest.raises(KeyError):
            laws["zmin"]

    def test_mapping_equality(self):
        a = FaceLaws.over(_XY, _laws(_XY), "where")
        b = FaceLaws.over(_XY, {face: _DISTINCT[face] for face in reversed(_XY)}, "where")
        assert a == b and a is not b
        assert a == dict(_laws(_XY))  # a mapping's equality
        other = _laws(_XY) | {"ymin": BC("albedo", {"albedo": 0.4})}
        assert a != FaceLaws.over(_XY, other, "where")
        assert a != FaceLaws.over(("xmin", "xmax"), _laws(("xmin", "xmax")), "where")

    def test_the_refusals_name_the_where(self):
        with pytest.raises(TypeError, match="where-tag is a mapping from face name to law"):
            FaceLaws.over(_XY, tuple(_DISTINCT.values()), "where-tag")
        with pytest.raises(ValueError, match="where-tag: this cartesian mesh has the boundary faces"):
            FaceLaws.over(_XY, _laws(("xmin", "xmax", "ymin")), "where-tag")

    def test_pickle_round_trip_of_the_value(self):
        import pickle

        laws = FaceLaws.over(_XY, {"xmin": BC.vacuum, "xmax": AlbedoBoundary(albedo=0.5),
                                   "ymin": BC("albedo", {"albedo": 0.3}), "ymax": VacuumInflow()}, "where")
        back = pickle.loads(pickle.dumps(laws))
        assert type(back) is FaceLaws
        assert back == laws and tuple(back) == tuple(laws)

    def test_pickle_round_trip_of_both_meshes(self):
        import pickle

        one = _mesh1d(_CYL, (0.25, 0.5, 1.0), {"xmin": BC.reflective, "xmax": AlbedoBoundary(albedo=0.5)})
        back1 = pickle.loads(pickle.dumps(one))
        assert back1 == one  # Mesh1D's bitwise equality, laws included
        two = _mesh2d("rz-hollow", face_laws={"xmin": BC.reflective, "xmax": BC.vacuum,
                                               "ymin": AlbedoBoundary(albedo=0.5), "ymax": VacuumInflow()})
        back2 = pickle.loads(pickle.dumps(two))
        assert back2.face_laws == two.face_laws and tuple(back2.face_laws) == tuple(two.face_laws)
        assert back2.coord is two.coord
        for name in ("edges_x", "edges_y", "mat_map"):
            np.testing.assert_array_equal(getattr(back2, name), getattr(two, name))


# ─────────────────────────────────────────────────────────────────────
# The named cells: the law is the model's (the user's ruling of 2026-09-30)
# ─────────────────────────────────────────────────────────────────────


class TestTheNamedCells:
    """``wigner_seitz_pin_cell`` is the cylindricalised lattice cell, white on
    its outer surface; ``pwr_slab_half_cell`` is bounded by two symmetry
    planes, reflective on both. Another law is another body
    (``StructuredGeometry.cylinder`` / ``.slab``). First reds: a
    ``boundaries=`` keyword re-added (the signature row); the model's law
    changed (the law row)."""

    @pytest.mark.parametrize("name", ["wigner_seitz_pin_cell", "pwr_slab_half_cell"])
    def test_no_law_keyword(self, name):
        parameters = inspect.signature(getattr(StructuredGeometry, name)).parameters
        assert "boundaries" not in parameters, name
        assert inspect.Parameter.VAR_KEYWORD not in {p.kind for p in parameters.values()}, name

    def test_the_models_laws(self):
        assert StructuredGeometry.wigner_seitz_pin_cell().boundaries == (BC("white"),)
        assert StructuredGeometry.pwr_slab_half_cell().boundaries == (BC("reflective"), BC("reflective"))
