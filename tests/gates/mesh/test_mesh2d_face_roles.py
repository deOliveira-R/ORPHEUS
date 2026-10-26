r"""The face roles of both meshes, and the one element parser (P1 step 3c, #405; re-spelled onto roles at wave 2, slice S1).

``Mesh2D(edges_x, edges_y, mat_map, *, face_roles, coord=CARTESIAN)`` declares one
role per boundary face (the user's ruling of 2026-09-29; the four ``bc_*`` fields
retire), and since the ruling of 2026-09-30 ``Mesh1D`` and ``Mesh2D`` store one
value, :class:`~orpheus.mesh.face_roles.FaceRoles`: a frozen, ordered, picklable
mapping from face name to role (a deck or an ``int`` tag) over the mesh's
:func:`~orpheus.mesh.face_roles.face_inventory`. The inventory rule: along the
first axis the faces are ``CoordSystem.boundary_points`` of its edges (both ends
on a slab, on (x, y) and on a hollow body; only the outer surface on a solid
radial body, whose centre is interior), along every further axis both ends,
each named by ``FaceLabel(axis_index, endpoint).face_name``. The positional
tuple, ``Mesh1D.boundary_faces`` (now ``boundary_points``) and
``Mesh2D.boundary_faces`` are retired. Two geometry claims that lived here, the inventory rule itself and the
named cells' roles, moved to ``tests/gates/geometry/test_geometry_faces.py`` in the layer-DAG re-filing (2026-10-26);
the two-mesh rows of ``TestTheOneInventoryRule`` stay. Wave 2 renamed ``face_laws`` to ``face_roles``, the law arm to ``Deck | int``,
and the one parser ``parse_boundary_law`` to :func:`~orpheus.geometry.deck.parse_face_role`; the deck rows and the
retirement rows live in ``test_face_roles.py``.

Claim kind: THEOREM for every row. The expected inventories are written by hand
from the topology law, independently of ``face_inventory``; the one-parser and
one-inventory rows are ROUTE gates (a counting spy on every binding of
``parse_face_role`` / ``face_inventory``), because a value or a message
alone cannot tell the one function from a twin that copied it.

First reds, by class (each named at its row, and each measured in the step-3c
batteries, ``scratch/reference_architecture/p1step3c/battery.md`` and
``battery2.md``):

* inventory: a solid (r, z) mesh given ``("xmin", "xmax")`` along x reds
  ``test_the_inventory``'s solid row and ``TestTheOneInventoryRule``;
* refusals: the pre-3c ``Mesh2D`` (four ``bc_*`` fields defaulting to ``None``)
  reds every row of ``TestRefusals``; ``face_roles`` given a default reds
  ``test_face_roles_is_required_and_keyword_only``;
* storage and value: a ``FaceRoles`` that stores (``__setitem__``), iterates in
  the caller's order, compares by identity or loses its ``__reduce__`` reds the
  matching ``TestFaceRoles`` row and the storage rows;
* the one parser: an inline ``isinstance`` check in any of the five classes in
  place of the ``parse_face_role`` call reds that class's row of
  ``TestTheOneElementParser::test_every_declared_role_passes_the_one_parser``.
"""
from __future__ import annotations

import dataclasses
import inspect
import sys
from collections.abc import Callable

import numpy as np
import pytest

import orpheus.geometry.deck as deck_module
from orpheus.geometry import CoordSystem, MirrorDeck, StructuredGeometry, TranslationDeck
from orpheus.mesh import AxisCoord, AxisMesh, Mesh1D, Mesh2D, RadialAxisMesh
from orpheus.mesh.face_roles import FaceRoles
from tests.gates.geometry._face_fixtures import FIRST_AXIS as _FIRST_AXIS

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/mesh/test_mesh2d_face_roles.py"
_NONE_IN_MESH1D = "tests/gates/mesh/test_mesh1d.py::TestConstructionLaws::test_refusal[none-role]"

_CART = CoordSystem.CARTESIAN
_CYL = CoordSystem.CYLINDRICAL

#: (edges_x, the inventory written by hand from the topology law).
_INVENTORIES = {
    "xy": (_CART, (0.0, 1.0, 2.0), ("xmin", "xmax", "ymin", "ymax")),
    "rz-solid": (_CYL, (0.0, 0.5, 1.0), ("xmax", "ymin", "ymax")),
    "rz-hollow": (_CYL, (0.25, 0.5, 1.0), ("xmin", "xmax", "ymin", "ymax")),
}

_EDGES_Y = (0.0, 1.0, 2.0, 3.0)

#: Four distinct tags, so a permutation of faces cannot read equal (a deck on a radial face is refused, so the
#: shared fixture is tags; the decks have their own rows).
_DISTINCT = {"xmin": 0, "xmax": 1, "ymin": 2, "ymax": 3}


def _roles(faces, source=_DISTINCT) -> dict:
    return {face: source[face] for face in faces}


def _mesh2d(case: str = "xy", **overrides) -> Mesh2D:
    coord, edges_x, faces = _INVENTORIES[case]
    kwargs = {"face_roles": _roles(faces), "coord": coord} | overrides
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
    assert isinstance(mesh.face_roles, FaceRoles)
    assert tuple(mesh.face_roles) == expected
    # the premise the hand-written inventory rests on, read from the one law
    n_x = len(coord.boundary_points(edges_x[0], edges_x[-1]))
    assert len(expected) == n_x + 2


def test_the_retired_fields_are_gone():
    """The four ``bc_*`` fields retire, and ``boundary_faces`` with them
    (``tuple(mesh.face_roles)`` are the names; ``retirement-audit`` D.18's
    witness). First red: any of them re-added."""
    names = {f.name for f in dataclasses.fields(Mesh2D)}
    assert names == {"edges_x", "edges_y", "mat_map", "face_roles", "coord"}, names
    for name in ("bc_xmin", "bc_xmax", "bc_ymin", "bc_ymax", "boundary_faces"):
        assert not hasattr(Mesh2D, name)
        assert not hasattr(_mesh2d(), name)


def test_face_roles_is_required_and_keyword_only():
    """First red: a default on ``face_roles`` (``required``), or the field moved
    before ``KW_ONLY`` (``keyword-only``)."""
    param = inspect.signature(Mesh2D).parameters["face_roles"]
    assert param.kind is inspect.Parameter.KEYWORD_ONLY
    assert param.default is inspect.Parameter.empty


# ─────────────────────────────────────────────────────────────────────
# Refusals, one negative each, each keyed to the input that triggers it
# ─────────────────────────────────────────────────────────────────────


class TestRefusals:
    """Each row is red on the pre-3c ``Mesh2D`` (which admitted ``None`` and had
    no ``face_roles``). The inventory is named in every count refusal."""

    def test_omitted(self):
        mat_map = np.zeros((2, 3), dtype=int)
        with pytest.raises(TypeError, match="face_roles"):
            Mesh2D(np.array([0.0, 1.0, 2.0]), np.asarray(_EDGES_Y), mat_map)  # type: ignore[call-arg]  # the refusal input

    def test_positional(self):
        mat_map = np.zeros((2, 3), dtype=int)
        roles = _roles(_INVENTORIES["xy"][2])
        with pytest.raises(TypeError, match="positional argument"):
            Mesh2D(np.array([0.0, 1.0, 2.0]), np.asarray(_EDGES_Y), mat_map, roles)  # type: ignore[misc]  # the refusal input

    @pytest.mark.parametrize("case", list(_INVENTORIES))
    def test_a_missing_face(self, case):
        _, _, faces = _INVENTORIES[case]
        roles = _roles(faces)
        del roles["ymax"]
        with pytest.raises(ValueError, match="Mesh2D.face_roles") as info:
            _mesh2d(case, face_roles=roles)
        for face in faces:  # the message names the inventory
            assert repr(face) in str(info.value), (face, str(info.value))

    @pytest.mark.parametrize("case", list(_INVENTORIES))
    def test_an_extra_face(self, case):
        _, _, faces = _INVENTORIES[case]
        roles = _roles(faces) | {"zmin": 0}
        with pytest.raises(ValueError, match="Mesh2D.face_roles") as info:
            _mesh2d(case, face_roles=roles)
        for face in faces:
            assert repr(face) in str(info.value), (face, str(info.value))

    def test_no_role_on_the_inner_surface_of_a_hollow_rz_mesh(self):
        """A hollow (r, z) without ``xmin``: the reason names the inner surface."""
        roles = _roles(("xmax", "ymin", "ymax"))
        with pytest.raises(ValueError, match="has an inner surface, which needs its own role"):
            _mesh2d("rz-hollow", face_roles=roles)

    def test_a_role_on_the_axis_of_a_solid_rz_mesh(self):
        """``xmin`` on a solid (r, z): the axis r = 0 is interior and carries no
        role. The reason is keyed: a hollow mesh missing ``ymax`` does not say it."""
        roles = _roles(("xmin", "xmax", "ymin", "ymax"))
        with pytest.raises(ValueError, match="carries no role") as info:
            _mesh2d("rz-solid", face_roles=roles)
        assert "'xmax', 'ymin', 'ymax'" in str(info.value)
        missing = _roles(_INVENTORIES["rz-hollow"][2])
        del missing["ymax"]
        with pytest.raises(ValueError) as other:
            _mesh2d("rz-hollow", face_roles=missing)
        assert "carries no role" not in str(other.value)

    @pytest.mark.parametrize("face", ["xmin", "ymax"])
    def test_none(self, face):
        roles = _roles(_INVENTORIES["xy"][2]) | {face: None}
        with pytest.raises(TypeError, match="a face role is a MirrorDeck, a TranslationDeck or an int tag, got None") as info:
            _mesh2d(face_roles=roles)
        assert repr(face) in str(info.value)

    @pytest.mark.parametrize("bad", ["vacuum", 1.5, True, object()], ids=["string", "float", "bool", "object"])
    def test_a_non_role_object(self, bad):
        roles = _roles(_INVENTORIES["xy"][2]) | {"ymin": bad}
        with pytest.raises(TypeError, match="a face role is a MirrorDeck, a TranslationDeck or an int tag") as info:
            _mesh2d(face_roles=roles)
        assert "'ymin'" in str(info.value)

    def test_the_positional_tuple(self):
        """The positional tuple is retired on both meshes (``Mesh1D``'s row is
        ``tests/gates/mesh/test_mesh1d.py`` ``tuple-laws``)."""
        with pytest.raises(TypeError, match="Mesh2D.face_roles is a mapping from face name to role"):
            _mesh2d(face_roles=tuple(_DISTINCT.values()))


# ─────────────────────────────────────────────────────────────────────
# Storage
# ─────────────────────────────────────────────────────────────────────


def test_the_storage_is_read_only():
    """First red: storing the caller's dict (or a copy that is a ``dict``)."""
    caller = _roles(_INVENTORIES["xy"][2])
    mesh = _mesh2d(face_roles=caller)
    with pytest.raises(TypeError):
        mesh.face_roles["xmin"] = 9  # type: ignore[index]  # the refusal input
    caller["xmin"] = 9  # the caller's mapping is not the mesh's
    assert mesh.face_roles["xmin"] == _DISTINCT["xmin"]
    with pytest.raises(dataclasses.FrozenInstanceError):
        mesh.face_roles = {}  # type: ignore[misc]  # the refusal input


@pytest.mark.parametrize("case", list(_INVENTORIES))
def test_the_storage_is_in_inventory_order(case):
    """Given in reverse, stored in inventory order, each role as declared.
    First red: storing in the caller's key order."""
    _, _, faces = _INVENTORIES[case]
    reversed_roles = {face: _DISTINCT[face] for face in reversed(faces)}
    mesh = _mesh2d(case, face_roles=reversed_roles)
    assert tuple(mesh.face_roles) == faces
    for face in faces:
        assert mesh.face_roles[face] == _DISTINCT[face], face


def test_a_deck_is_admitted_beside_a_tag():
    """A deck and a tag on one mesh, each stored as declared (the translation pair closes on ``x``).
    First red: a parser admitting only ints."""
    roles = {"xmin": TranslationDeck(), "xmax": TranslationDeck(), "ymin": MirrorDeck(), "ymax": 3}
    mesh = _mesh2d(face_roles=roles)
    for face, role in roles.items():
        assert mesh.face_roles[face] == role


# ─────────────────────────────────────────────────────────────────────
# The one element parser: Mesh1D, Mesh2D, AxisMesh, RadialAxisMesh and
# StructuredGeometry all declare through ``parse_face_role``
# ─────────────────────────────────────────────────────────────────────

_SLAB_EDGES = np.array([0.0, 0.5, 2.0])
_SPH_EDGES = np.array([0.0, 0.5, 2.0])


def _build_mesh1d(low, high):
    return Mesh1D(coord=_CART, edges=_SLAB_EDGES, volumes=_CART.measure(_SLAB_EDGES),
                  region_ids=np.zeros(2, dtype=int), region_materials=(0,), face_roles={"xmin": low, "xmax": high})


def _build_mesh2d(low, high):
    return _mesh2d(face_roles={"xmin": low, "xmax": high, "ymin": 2, "ymax": 2})


def _build_axis(low, high):
    return AxisMesh(edges=_SLAB_EDGES, role_low=low, role_high=high)


def _build_radial(low, high):
    del high  # one endpoint
    return RadialAxisMesh(edges=_SPH_EDGES, coord=AxisCoord.RADIAL_SPHERICAL, role_outer=low)


def _build_geometry(low, high):
    return StructuredGeometry(coord=_CART, breakpoints=(0.0, 2.0), mat_ids=(0,), face_roles=(low, high))


#: (owner, the builder, the slot name the refusal names).
_OWNERS: list[tuple[str, Callable[[object, object], object], str]] = [
    ("Mesh1D", _build_mesh1d, "Mesh1D.face_roles['xmin']"),
    ("Mesh2D", _build_mesh2d, "Mesh2D.face_roles['xmin']"),
    ("AxisMesh", _build_axis, "AxisMesh.role_low"),
    ("RadialAxisMesh", _build_radial, "RadialAxisMesh.role_outer"),
    ("StructuredGeometry", _build_geometry, "StructuredGeometry.face_roles['xmin']"),
]


class TestTheOneElementParser:
    """One check of a single boundary declaration (the design of step 3c).

    Two legs: the ROUTE (every declared role passes through
    ``parse_face_role``, counted by a spy on every module binding of it) and
    the MESSAGE (``None`` is refused with the parser's fragment, raised in the
    parser's own frame, naming the owner's slot)."""

    @pytest.mark.parametrize("owner, build, slot", _OWNERS, ids=[o[0] for o in _OWNERS])
    def test_every_declared_role_passes_the_one_parser(self, owner, build, slot, monkeypatch):
        del slot
        original = deck_module.parse_face_role
        seen: list[int] = []

        def spy(role, where):
            seen.append(id(role))
            return original(role, where)

        rebound = 0
        for module in list(sys.modules.values()):
            if getattr(module, "parse_face_role", None) is original:
                monkeypatch.setattr(module, "parse_face_role", spy)
                rebound += 1
        # the spy is installed where the constructors read it: the defining
        # module and the axis module at least
        assert rebound >= 2 and deck_module.parse_face_role is spy, rebound
        low, high = 777, MirrorDeck()
        build(low, high)
        assert id(low) in seen, f"{owner}: its role never reached parse_face_role (rebound {rebound})"
        if owner != "RadialAxisMesh":
            assert id(high) in seen, f"{owner}: its second role never reached parse_face_role"

    @pytest.mark.rests_on(_NONE_IN_MESH1D)
    @pytest.mark.parametrize("owner, build, slot", _OWNERS, ids=[o[0] for o in _OWNERS])
    def test_none_is_refused_by_the_one_parser(self, owner, build, slot):
        with pytest.raises(TypeError, match="a face role is a MirrorDeck, a TranslationDeck or an int tag, got None") as info:
            build(None, 0)
        assert slot in str(info.value), (owner, str(info.value))
        assert info.traceback[-1].name == "parse_face_role", info.traceback[-1].name

    @pytest.mark.parametrize("owner, build, slot", _OWNERS, ids=[o[0] for o in _OWNERS])
    def test_a_non_role_is_refused_by_the_one_parser(self, owner, build, slot):
        with pytest.raises(TypeError, match="a face role is a MirrorDeck, a TranslationDeck or an int tag") as info:
            build("vacuum", 0)
        assert slot in str(info.value), (owner, str(info.value))
        assert info.traceback[-1].name == "parse_face_role", info.traceback[-1].name


# ─────────────────────────────────────────────────────────────────────
# One inventory rule for both meshes (the user's ruling of 2026-09-30:
# one FaceRoles value for Mesh1D and Mesh2D)
# ─────────────────────────────────────────────────────────────────────

#: (coord, first-axis edges, the first-axis faces written by hand): ``_FIRST_AXIS``, imported from
#: ``tests/gates/geometry/_face_fixtures.py``, which the rule's own row
#: (``tests/gates/geometry/test_geometry_faces.py::TestTheOneInventoryRule``) reads too.


def _mesh1d(coord, edges, roles):
    e = np.asarray(edges, dtype=float)
    return Mesh1D(coord=coord, edges=e, volumes=coord.measure(e), region_ids=np.zeros(len(e) - 1, dtype=int), region_materials=(0,),
                  face_roles=roles)


class TestTheOneInventoryRule:
    """``Mesh1D`` and ``Mesh2D`` on the same coordinate system and first-axis
    edges name the same first-axis faces, and a 2-D mesh adds exactly
    ``ymin``, ``ymax``. Two legs: the VALUE (against the hand-written
    inventory) and the ROUTE (both constructors ask ``face_inventory``).
    First red: either mesh spelling its own inventory (a ``Mesh2D`` whose
    solid radial axis is given ``xmin``/``xmax``, the pre-FaceRoles
    ``_faces_2d`` shape, or a constructor bypassing ``face_inventory``)."""

    @pytest.mark.parametrize("case", list(_FIRST_AXIS))
    def test_both_meshes_name_the_same_first_axis_faces(self, case):
        coord, edges, faces = _FIRST_AXIS[case]
        one = _mesh1d(coord, edges, _roles(faces))
        two = Mesh2D(np.asarray(edges), np.asarray(_EDGES_Y), np.zeros((2, 3), dtype=int),
                     face_roles=_roles(faces + ("ymin", "ymax")), coord=coord)
        assert tuple(one.face_roles) == faces
        assert tuple(two.face_roles) == faces + ("ymin", "ymax")
        assert isinstance(one.face_roles, FaceRoles) and isinstance(two.face_roles, FaceRoles)

    def test_both_constructors_ask_the_one_rule(self, monkeypatch):
        import orpheus.mesh.face_roles as face_roles_module

        original = face_roles_module.face_inventory
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
        _mesh1d(_CART, (0.0, 1.0), _roles(("xmin", "xmax")))
        assert dimensions == [1], dimensions
        _mesh2d()
        assert dimensions == [1, 2], dimensions



# ─────────────────────────────────────────────────────────────────────
# FaceRoles: the value both meshes store
# ─────────────────────────────────────────────────────────────────────

_XY = ("xmin", "xmax", "ymin", "ymax")


class TestFaceRoles:
    """The laws of the type (``vv-testing``: a math-bearing type ships the test
    of its defining laws): inventory order, immutability, mapping equality,
    pickle round trip. First reds: ``__iter__`` over the caller's order (the
    order row); a ``__setattr__`` or ``__setitem__`` that stores (the
    immutability row); an identity ``__eq__`` (the equality row); a
    ``__reduce__`` dropped on a slotted class, or one that rebuilds without the
    roles (the pickle rows)."""

    def test_order_is_the_inventory(self):
        given = {face: _DISTINCT[face] for face in reversed(_XY)}
        roles = FaceRoles.over(_XY, given, "where", coord=_CART)
        assert tuple(roles) == _XY
        assert list(roles.items()) == [(face, _DISTINCT[face]) for face in _XY]
        assert len(roles) == 4

    def test_immutability(self):
        roles = FaceRoles.over(_XY, _roles(_XY), "where", coord=_CART)
        with pytest.raises(TypeError):
            roles["xmin"] = 9  # type: ignore[index]  # the refusal input
        with pytest.raises(AttributeError):
            roles._items = ()  # type: ignore[misc]  # the refusal input
        with pytest.raises(AttributeError):
            roles.extra = 1  # type: ignore[attr-defined]  # the refusal input
        assert roles["xmin"] == _DISTINCT["xmin"]
        with pytest.raises(KeyError):
            roles["zmin"]

    def test_mapping_equality(self):
        a = FaceRoles.over(_XY, _roles(_XY), "where", coord=_CART)
        b = FaceRoles.over(_XY, {face: _DISTINCT[face] for face in reversed(_XY)}, "where", coord=_CART)
        assert a == b and a is not b
        # Content identity (#405 P1 step 5): equal to another FaceRoles of the
        # same items in any order, and no longer to a plain dict (until then
        # a Mapping's equality made it so); its items are the dict.
        assert a != dict(_roles(_XY)) and dict(a) == dict(_roles(_XY))
        assert hash(a) == hash(b)
        other = _roles(_XY) | {"ymin": 4}
        assert a != FaceRoles.over(_XY, other, "where", coord=_CART)
        assert a != FaceRoles.over(("xmin", "xmax"), _roles(("xmin", "xmax")), "where", coord=_CART)

    def test_the_refusals_name_the_where(self):
        with pytest.raises(TypeError, match="where-tag is a mapping from face name to role"):
            FaceRoles.over(_XY, tuple(_DISTINCT.values()), "where-tag", coord=_CART)
        with pytest.raises(ValueError, match="where-tag: this cartesian mesh has the boundary faces"):
            FaceRoles.over(_XY, _roles(("xmin", "xmax", "ymin")), "where-tag", coord=_CART)

    def test_pickle_round_trip_of_the_value(self):
        import pickle

        roles = FaceRoles.over(_XY, {"xmin": MirrorDeck(), "xmax": MirrorDeck(), "ymin": TranslationDeck(), "ymax": TranslationDeck()},
                               "where", coord=_CART)
        back = pickle.loads(pickle.dumps(roles))
        assert type(back) is FaceRoles
        assert back == roles and tuple(back) == tuple(roles)

    def test_pickle_round_trip_of_both_meshes(self):
        import pickle

        one = _mesh1d(_CYL, (0.25, 0.5, 1.0), {"xmin": 1, "xmax": 0})
        back1 = pickle.loads(pickle.dumps(one))
        assert back1 == one  # Mesh1D's bitwise equality, roles included
        two = _mesh2d("rz-hollow", face_roles={"xmin": 1, "xmax": 0, "ymin": MirrorDeck(), "ymax": MirrorDeck()})
        back2 = pickle.loads(pickle.dumps(two))
        assert back2.face_roles == two.face_roles and tuple(back2.face_roles) == tuple(two.face_roles)
        assert back2.coord is two.coord
        for name in ("edges_x", "edges_y", "mat_map"):
            np.testing.assert_array_equal(getattr(back2, name), getattr(two, name))
