r"""The face roles of both meshes, the geometry's symmetry on a mesh, and the retirement of ``FaceLaws``
(wave 2 of the layer DAG, WF2.3 slice S1; ``.claude/plans/wave2_boundary_carve.md`` sections 2.2, 2.3 and 6.1).

A mesh stores one :class:`~orpheus.mesh.face_roles.FaceRoles` (a frozen, ordered mapping from face name to a
``Deck | int`` role over the mesh's face inventory) where it stored ``FaceLaws``; ``Mesh1D`` also answers
``symmetry``, the colouring group of its REGIONS, so that a mesh and its geometry give one group whatever the cell
partition. Claims, with layer and kind (``vv-principles``):

* ``TestFaceRolesOver``: THEOREM rows of the parse (the inventory rule, the deck admission per axis and chart, the
  pairing closure, the bool refusal) with the old refusal messages kept; the inventory is written by hand from the
  topology law, as ``test_mesh2d_face_roles.py`` does;
* ``TestMeshesStoreFaceRoles`` (gates 3 and 4, mesh half): content identity of a role, the Mesher's pass-through,
  the factory's single ``role``; REFERENCE (the geometry's own roles and tags are the expected values);
* ``TestTheSymmetryIsTheGeometrys`` (gate 6, law L-C input): THEOREM (refinement independence: the group is a
  function of the colouring and the end roles) against a REFERENCE that shares nothing with the mesh route, a direct
  call of ``colouring_symmetry`` on the geometry's raw breakpoints, ids and roles; the activation is a palindromic
  triple ``(0, 1, 0)`` whose interior breakpoints are the regions' edges, with an UNSYMMETRIC cell partition, so a
  reading of the cell edges instead of the region edges is a different group;
* ``TestFaceLawsIsRetired``: STRUCTURE (no dual spelling);
* ``TestFaceRoles``: THEOREM rows of a mesh's stored roles, through the mesher (1-D) and the bare constructor (2-D):
  one role per face, a deck beside a tag, every 2-D face admitting a tag, the frozen field, and the refusal of a
  string role. It lived in ``tests/gates/geometry/test_geometry.py`` until the layer-DAG re-filing
  (``.claude/plans/layer_dag.md``, ruling 21).

Every new attribute is read at run time, so each row's red is its own. First reds at the base tip ``dec1f9ff``:
``orpheus.mesh.face_roles`` does not exist (every ``TestFaceRolesOver`` row); ``Mesh1D(... face_roles=...)`` is an
unexpected keyword (``TestMeshesStoreFaceRoles``); ``Mesh1D`` has no ``symmetry`` (gate 6); ``orpheus.mesh.face_laws``
and ``FaceLaws`` exist (``TestFaceLawsIsRetired``).

Mutations, for the reviewers (each on a COPY of the module named):

* ``FaceRoles.over`` skipping ``admitted_decks`` (the radial deck rows), skipping ``check_pairings_close`` (the
  translation-pair rows), a default of ``coord=CARTESIAN`` (the signature row), losing the ``_why`` suffixes;
* ``Mesh1D.symmetry`` that passes the CELL edges and the cell material ids instead of the region edges and region
  materials is bit-neutral (the run merge removes the difference, by design), so the mutants that bite are: always
  the identity (reds the symmetric rows), always the mirror (reds ``(0, 1, 2)`` and the unequal-roles row), the roles
  dropped from the call (reds the ``1`` versus ``2`` row), the group of ``mesh.edges`` ends swapped with the roles'
  order;
* ``Mesher.partition`` passing the geometry's roles in a different order from the inventory (reds the hollow rows).
"""

from __future__ import annotations

import ast
import importlib
import importlib.util
import inspect
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import orpheus.mesh as mesh_package
from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.geometry.colouring import UncomputedSymmetries, colouring_symmetry
from orpheus.geometry.deck import DeckError, MirrorDeck, TranslationDeck, tags_of
from orpheus.mesh import CellsByCount, Mesh1D, Mesh2D, Mesher
from orpheus.values.content import content_digest
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/mesh/test_face_roles.py"
_GEOMETRY_ROLES = "tests/gates/geometry/test_face_roles.py"
_ROOT = Path(__file__).resolve().parents[3]

_CART, _CYL, _SPH = CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL
_REFUSAL = r"N\(G_c\)/G_c"


class _Lazy:
    def __init__(self, name: str) -> None:
        self._name = name

    def __getattr__(self, attr: str) -> Any:
        return getattr(importlib.import_module(self._name), attr)


face_roles: Any = _Lazy("orpheus.mesh.face_roles")


# ── FaceRoles.over ───────────────────────────────────────────────────────────

#: (coord, the inventory written by hand from the topology law).
_INVENTORIES = {
    "slab": (_CART, ("xmin", "xmax")),
    "solid cylinder (r)": (_CYL, ("xmax",)),
    "hollow sphere (r)": (_SPH, ("xmin", "xmax")),
}


@pytest.mark.rests_on(_GEOMETRY_ROLES + "::TestACurvedFaceCarriesNoDeck::test_the_slab_admits_both_decks_and_closes_the_pairing")
class TestFaceRolesOver:
    """The parse of a mesh's roles over its inventory."""

    def test_coord_is_required_and_keyword_only(self) -> None:
        """C9: a default of ``CARTESIAN`` would admit a mirror on a radial face for a caller that omits it."""
        parameter = inspect.signature(face_roles.FaceRoles.over).parameters["coord"]
        require(parameter.kind is inspect.Parameter.KEYWORD_ONLY, f"coord is {parameter.kind}")
        require(parameter.default is inspect.Parameter.empty, f"coord has the default {parameter.default!r}")
        with pytest.raises(TypeError, match="coord"):
            face_roles.FaceRoles.over(("xmin", "xmax"), {"xmin": 0, "xmax": 0}, "m.face_roles")

    @pytest.mark.parametrize("label", list(_INVENTORIES))
    def test_tags_build_in_inventory_order(self, label: str) -> None:
        coord, inventory = _INVENTORIES[label]
        given = {face: k for k, face in enumerate(reversed(inventory))}  # the caller's order is the reverse
        built = face_roles.FaceRoles.over(inventory, given, "m.face_roles", coord=coord)
        require(tuple(built) == inventory, f"iteration is the inventory order: {tuple(built)}")
        require(dict(built) == given, "and the roles are the caller's")
        require(built.tags == frozenset(given.values()), f"tags {built.tags}")

    def test_a_deck_is_admitted_where_the_chart_admits_it(self) -> None:
        slab = face_roles.FaceRoles.over(("xmin", "xmax"), {"xmin": MirrorDeck(), "xmax": 4}, "m.face_roles", coord=_CART)
        require(slab["xmin"] == MirrorDeck() and slab.tags == frozenset({4}), f"{dict(slab)}")
        rz = ("xmax", "ymin", "ymax")
        axial = face_roles.FaceRoles.over(rz, {"xmax": 0, "ymin": MirrorDeck(), "ymax": MirrorDeck()}, "m.face_roles", coord=_CYL)
        require(axial["ymin"] == MirrorDeck(), "an axial face of an (r, z) mesh admits the mirror")
        face_roles.FaceRoles.over(rz, {"xmax": 0, "ymin": TranslationDeck(), "ymax": TranslationDeck()}, "m.face_roles", coord=_CYL)

    @pytest.mark.parametrize("deck", [MirrorDeck(), TranslationDeck()], ids=["mirror", "translation"])
    @pytest.mark.parametrize("label", ["solid cylinder (r)", "hollow sphere (r)"])
    def test_a_deck_on_a_radial_face_is_refused(self, label: str, deck: Any) -> None:
        coord, inventory = _INVENTORIES[label]
        control = {face: 0 for face in inventory}
        face_roles.FaceRoles.over(inventory, control, "m.face_roles", coord=coord)
        for face in inventory:
            with pytest.raises(DeckError, match=_REFUSAL):
                face_roles.FaceRoles.over(inventory, {**control, face: deck}, "m.face_roles", coord=coord)

    def test_a_translation_needs_its_opposite(self) -> None:
        pair = {"xmin": TranslationDeck(), "xmax": TranslationDeck()}
        face_roles.FaceRoles.over(("xmin", "xmax"), pair, "m.face_roles", coord=_CART)
        for other in (0, MirrorDeck()):
            with pytest.raises(DeckError, match="TranslationDeck"):
                face_roles.FaceRoles.over(("xmin", "xmax"), {"xmin": TranslationDeck(), "xmax": other}, "m.face_roles", coord=_CART)
        four = ("xmin", "xmax", "ymin", "ymax")
        face_roles.FaceRoles.over(four, {"xmin": TranslationDeck(), "xmax": TranslationDeck(), "ymin": 1, "ymax": 2}, "m.face_roles", coord=_CART)
        with pytest.raises(DeckError, match="ymin"):
            face_roles.FaceRoles.over(four, {"xmin": 0, "xmax": 0, "ymin": TranslationDeck(), "ymax": 2}, "m.face_roles", coord=_CART)

    def test_the_refusals_keep_their_vocabulary(self) -> None:
        over = face_roles.FaceRoles.over
        with pytest.raises(ValueError, match="boundary faces"):
            over(("xmin", "xmax"), {"xmin": 0}, "m.face_roles", coord=_CART)
        with pytest.raises(ValueError, match="interior point"):
            over(("xmax",), {"xmin": 0, "xmax": 0}, "m.face_roles", coord=_CYL)
        with pytest.raises(ValueError, match="inner surface"):
            over(("xmin", "xmax"), {"xmax": 0}, "m.face_roles", coord=_SPH)
        with pytest.raises(TypeError, match="mapping"):
            over(("xmin", "xmax"), (0, 0), "m.face_roles", coord=_CART)
        for bad in (True, np.False_, None, 1.0, "vacuum"):
            with pytest.raises(TypeError, match=r"m\.face_roles\['xmax'\]"):
                over(("xmin", "xmax"), {"xmin": 0, "xmax": bad}, "m.face_roles", coord=_CART)

    def test_it_is_a_content_value(self) -> None:
        over = face_roles.FaceRoles.over
        a = over(("xmin", "xmax"), {"xmin": MirrorDeck(), "xmax": 1}, "w", coord=_CART)
        b = over(("xmin", "xmax"), {"xmax": 1, "xmin": MirrorDeck()}, "w", coord=_CART)
        require(a == b and hash(a) == hash(b) and content_digest(a) == content_digest(b), "the order the caller gives is not content")
        c = over(("xmin", "xmax"), {"xmin": 1, "xmax": 1}, "w", coord=_CART)
        require(a != c and content_digest(a) != content_digest(c), "a deck is not the tag")
        require(a != {"xmin": MirrorDeck(), "xmax": 1}, "and it is equal only to another FaceRoles")


# ── the meshes ───────────────────────────────────────────────────────────────

_EDGES = np.array([0.0, 1.0, 2.0])


def _mesh1d(left: object = 0, right: object = 1) -> Mesh1D:
    geometry = StructuredGeometry.slab((0.0, 1.0, 2.0), (0, 1), left=left, right=right)  # type: ignore[arg-type]
    return Mesher(geometry).partition(CellsByCount.uniform_width(2)).mesh


def _mesh2d(roles: dict[str, object] | None = None, coord: CoordSystem = _CART) -> Mesh2D:
    given = {"xmin": 0, "xmax": 1, "ymin": 2, "ymax": 3} if roles is None else roles
    return Mesh2D(_EDGES, np.array([0.0, 1.0]), np.array([[0], [1]]), face_roles=given, coord=coord)  # type: ignore[call-arg]


@pytest.mark.rests_on(_HERE + "::TestFaceRolesOver::test_it_is_a_content_value")
class TestMeshesStoreFaceRoles:
    """``Mesh1D`` and ``Mesh2D`` store ``face_roles``; the Mesher lifts the geometry's roles; the factory takes ``role``."""

    def test_the_mesh_stores_a_face_roles_value(self) -> None:
        for mesh in (_mesh1d(MirrorDeck(), 5), _mesh2d()):
            require(isinstance(mesh.face_roles, face_roles.FaceRoles), f"{type(mesh).__name__}.face_roles is {type(mesh.face_roles).__name__}")
        one = _mesh1d(MirrorDeck(), 5)
        require(dict(one.face_roles) == {"xmin": MirrorDeck(), "xmax": 5}, f"{dict(one.face_roles)}")
        require(one.outer_role == 5, "outer_role is the xmax role")
        require(dict(_mesh2d().face_roles) == {"xmin": 0, "xmax": 1, "ymin": 2, "ymax": 3}, "the four faces")

    def test_a_role_is_content_of_the_mesh(self) -> None:
        require(_mesh1d() == _mesh1d() and content_digest(_mesh1d()) == content_digest(_mesh1d()), "control: two builds are one mesh")
        for label, other in (("a tag", _mesh1d(left=7)), ("a deck", _mesh1d(left=MirrorDeck())), ("the outer tag", _mesh1d(right=2))):
            require(other != _mesh1d() and content_digest(other) != content_digest(_mesh1d()), f"{label} moves the mesh")
        require(_mesh2d() == _mesh2d() and _mesh2d() != _mesh2d({"xmin": 0, "xmax": 1, "ymin": 2, "ymax": 4}), "one role of four, Mesh2D")
        require(_mesh2d({"xmin": 0, "xmax": 1, "ymin": 2, "ymax": 3}) == _mesh2d({"ymax": 3, "ymin": 2, "xmax": 1, "xmin": 0}), "the order is not content")

    def test_the_old_spellings_do_not_construct(self) -> None:
        control = _mesh2d()  # the new spelling builds
        require(control.face_roles["xmin"] == 0, "control")
        with pytest.raises(TypeError, match="face_laws"):
            Mesh2D(_EDGES, np.array([0.0, 1.0]), np.array([[0], [1]]), face_laws={"xmin": 0, "xmax": 1, "ymin": 2, "ymax": 3})  # type: ignore[call-arg]
        require(not hasattr(control, "face_laws"), "no face_laws attribute on Mesh2D")
        mesh = _mesh1d()
        require(not hasattr(mesh, "face_laws") and not hasattr(mesh, "outer_law"), "no face_laws or outer_law on Mesh1D")
        parameter = inspect.signature(Mesh2D).parameters["face_roles"]
        require(parameter.kind is inspect.Parameter.KEYWORD_ONLY and parameter.default is inspect.Parameter.empty, "face_roles is required and keyword-only")

    def test_a_radial_mesh_refuses_a_deck(self) -> None:
        with pytest.raises(DeckError, match=_REFUSAL):
            _mesh2d({"xmax": MirrorDeck(), "ymin": 0, "ymax": 0}, coord=_CYL)
        built = _mesh2d({"xmax": 3, "ymin": MirrorDeck(), "ymax": MirrorDeck()}, coord=_CYL)
        require(tuple(built.face_roles) == ("xmax", "ymin", "ymax"), "an (r, z) mesh admits a mirror on its axial faces")

    @pytest.mark.parametrize(
        ("coord", "breakpoints", "kwargs", "expected"),
        [
            (_CART, (0.0, 1.0, 2.0), {"left": MirrorDeck(), "right": 3}, {"xmin": MirrorDeck(), "xmax": 3}),
            (_CYL, (0.0, 1.0, 2.0), {"outer": 4}, {"xmax": 4}),
            (_SPH, (0.5, 1.0, 2.0), {"outer": 0, "inner": 1}, {"xmin": 1, "xmax": 0}),
        ],
        ids=["slab", "solid cylinder", "hollow sphere"],
    )
    def test_the_mesher_lifts_the_geometrys_roles(self, coord: CoordSystem, breakpoints: tuple[float, ...], kwargs: dict[str, Any], expected: dict[str, Any]) -> None:
        build = {_CART: StructuredGeometry.slab, _CYL: StructuredGeometry.cylinder, _SPH: StructuredGeometry.sphere}[coord]
        geometry = build(breakpoints, (0, 1), **kwargs)
        mesh = Mesher(geometry).partition(CellsByCount.uniform_width(2)).mesh
        require(dict(mesh.face_roles) == expected, f"{dict(mesh.face_roles)} vs {expected}")
        require(tags_of(mesh.face_roles.values()) == geometry.tags, f"mesh tags {tags_of(mesh.face_roles.values())} vs geometry tags {geometry.tags}")

    def test_the_factory_takes_one_role(self) -> None:
        from orpheus.mesh import pwr_pin_2d

        for role in (0, MirrorDeck(), TranslationDeck()):
            mesh = pwr_pin_2d(role=role, n_cells=4)
            require(dict(mesh.face_roles) == {face: role for face in ("xmin", "xmax", "ymin", "ymax")}, f"{role!r}: {dict(mesh.face_roles)}")
        with pytest.raises(TypeError, match="law"):
            pwr_pin_2d(law=0, n_cells=4)  # type: ignore[call-arg]
        with pytest.raises(TypeError, match="role"):
            pwr_pin_2d(n_cells=4)  # type: ignore[call-arg]


# ── gate 6 ───────────────────────────────────────────────────────────────────

def _geometry(breakpoints: tuple[float, ...], ids: tuple[int, ...], low: object, high: object) -> StructuredGeometry:
    return StructuredGeometry.slab(breakpoints, ids, left=low, right=high)  # type: ignore[arg-type]


def _meshed(geometry: StructuredGeometry, how: str) -> Mesh1D:
    n = len(geometry.intervals)
    if how == "refined":
        return Mesher(geometry).partition(CellsByCount.uniform_width(2)).refine(2).mesh
    rule = CellsByCount.uniform_width(3) if how == "uniform" else tuple(
        CellsByCount.uniform_width(k) for k in (2, 3, 7)[:n]  # cells per interval: the cell edges are not symmetric
    )
    return Mesher(geometry).partition(rule).mesh


def _reference(geometry: StructuredGeometry) -> Any:
    """The expected group, by a direct call on the geometry's raw data (no mesh anywhere on this route)."""
    return colouring_symmetry(
        _CART, (geometry.breakpoints[0], geometry.breakpoints[-1]), geometry.breakpoints, geometry.mat_ids,
        geometry.face_roles[0], geometry.face_roles[1],
    )


#: (label, breakpoints, ids, low role, high role, number of group elements). The counts are written by hand:
#: a palindromic colouring with equal end roles has the identity and the mirror, anything else the identity only.
_SLABS = [
    ("one material, one tag", (0.0, 1.0, 2.0), (0, 0), 1, 1, 2),
    ("palindrome, equal tags", (0.0, 1.0, 2.0, 3.0), (0, 1, 0), 1, 1, 2),
    ("palindrome, two mirrors", (0.0, 1.0, 2.0, 3.0), (0, 1, 0), MirrorDeck(), MirrorDeck(), 2),
    ("the float trap (0, .1, .2, .3)", (0.0, 0.1, 0.2, 0.3), (0, 1, 0), 1, 1, 2),
    ("asymmetric materials", (0.0, 1.0, 2.0, 3.0), (0, 1, 2), 1, 1, 1),
    ("asymmetric breakpoints", (0.0, 1.0, 2.5, 3.0), (0, 1, 0), 1, 1, 1),
    ("unequal tags", (0.0, 1.0, 2.0, 3.0), (0, 1, 0), 1, 2, 1),
    ("a mirror and a tag", (0.0, 1.0, 2.0, 3.0), (0, 1, 0), MirrorDeck(), 1, 1),
]


@pytest.mark.rests_on(_HERE + "::TestMeshesStoreFaceRoles::test_the_mesher_lifts_the_geometrys_roles", "tests/gates/geometry/test_colouring_symmetry.py::test_colouring_symmetry_reads_runs_and_roles")
class TestTheSymmetryIsTheGeometrys:
    """Gate 6: ``Mesh1D.symmetry == StructuredGeometry.symmetry`` for every partition of the same geometry."""

    @pytest.mark.parametrize("label", [row[0] for row in _SLABS])
    def test_the_geometry_gives_the_reference_group(self, label: str) -> None:
        (_, breakpoints, ids, low, high, count), = [row for row in _SLABS if row[0] == label]
        geometry = _geometry(breakpoints, ids, low, high)
        require(len(geometry.symmetry.elements) == count, f"{label}: {len(geometry.symmetry.elements)} elements, expected {count}")
        require(geometry.symmetry == _reference(geometry), f"{label}: the geometry's group is the direct call's")
        require(geometry.symmetry.elements[0].is_identity, "the identity is first")

    @pytest.mark.parametrize("how", ["uniform", "unsymmetric", "refined"])
    @pytest.mark.parametrize("label", [row[0] for row in _SLABS])
    def test_the_mesh_gives_the_same_group(self, label: str, how: str) -> None:
        (_, breakpoints, ids, low, high, count), = [row for row in _SLABS if row[0] == label]
        geometry = _geometry(breakpoints, ids, low, high)
        mesh = _meshed(geometry, how)
        require(mesh.symmetry == geometry.symmetry, f"{label} / {how}: mesh and geometry disagree")
        require(len(mesh.symmetry.elements) == count, f"{label} / {how}: {len(mesh.symmetry.elements)} elements, expected {count}")

    def test_the_mirror_swaps_the_faces(self) -> None:
        mirror = _meshed(_geometry((0.0, 1.0, 2.0), (0, 0), 1, 1), "unsymmetric").symmetry.elements[1]
        require(dict(mirror.faces) == {"xmin": "xmax", "xmax": "xmin"}, f"{dict(mirror.faces)}")

    def test_a_radial_mesh_has_the_identity_only(self) -> None:
        geometry = StructuredGeometry.cylinder((0.0, 1.0, 2.0), (0, 0), outer=0)
        mesh = Mesher(geometry).partition(CellsByCount.uniform_width(3)).mesh
        require(len(mesh.symmetry.elements) == 1 and mesh.symmetry == geometry.symmetry, "W is trivial on a radial chart")

    def test_the_2d_mesh_says_it_did_not_compute(self) -> None:
        require(_mesh2d().symmetry == UncomputedSymmetries("Mesh2D"), f"{_mesh2d().symmetry!r}")
        require(not isinstance(_mesh2d().symmetry, type(_mesh1d().symmetry)), "and it is not a group of one element")


# ── the retirement ───────────────────────────────────────────────────────────


class TestFaceLawsIsRetired:
    """No dual spelling. Layer: architecture; kind STRUCTURE."""

    def test_the_module_and_the_class_are_gone(self) -> None:
        require(importlib.util.find_spec("orpheus.mesh.face_laws") is None, "orpheus.mesh.face_laws is renamed, not copied")
        require(not hasattr(mesh_package, "FaceLaws") and "FaceLaws" not in mesh_package.__all__, "FaceLaws is not exported")
        require(hasattr(mesh_package, "FaceRoles") and "FaceRoles" in mesh_package.__all__, "FaceRoles is")

    def test_the_inventory_has_one_home(self) -> None:
        from orpheus.mesh import face_inventory

        require(face_roles.face_inventory is face_inventory, "orpheus.mesh.face_inventory is the module's own")
        require(face_inventory(_CYL, np.array([0.0, 1.0]), 2) == ("xmax", "ymin", "ymax"), "unchanged: a solid (r, z) mesh")

    def test_no_source_names_the_retired_spellings(self) -> None:
        """No ``FaceLaws``, ``face_laws`` or ``outer_law`` identifier in ``orpheus/mesh`` (an ``ast`` walk; strings and
        docstrings excluded, because prose may say what was retired)."""
        found: list[str] = []
        files = sorted((_ROOT / "orpheus" / "mesh").rglob("*.py"))
        require(len(files) >= 7, f"{len(files)} files scanned")
        for path in files:
            for node in ast.walk(ast.parse(path.read_text())):
                names = []
                if isinstance(node, ast.Name):
                    names = [node.id]
                elif isinstance(node, ast.Attribute):
                    names = [node.attr]
                elif isinstance(node, ast.arg):
                    names = [node.arg]
                elif isinstance(node, ast.keyword) and node.arg:
                    names = [node.arg]
                elif isinstance(node, ast.alias):
                    names = [node.name.split(".")[-1]]
                found += [f"{path.name}:{getattr(node, 'lineno', '?')} {n}" for n in names if n in {"FaceLaws", "face_laws", "outer_law", "bc_low", "bc_high", "bc_outer"}]
        require(not found, "retired identifiers remain:\n  " + "\n  ".join(found))


# ═══════════════════════════════════════════════════════════════════════
# The mesh's face roles
# ═══════════════════════════════════════════════════════════════════════

class TestFaceRoles:
    """The mesh's face roles (the ``BC`` declaration dataclass and its ``to_alpha`` bridge retired in WF2.4)."""

    # ── Mesh1D face roles ────────────────────────────────────────────
    # RE-POSED at P1 step 3b: ``bc_left``/``bc_right`` and the ``None``
    # default retired; a mesh carries one role per boundary face (a deck or a tag, wave 2 S1). The
    # ``None`` default row and the "backward compatible" row are deleted
    # (their subject was the retired default); ``None`` is refused on a face
    # by ``tests/gates/mesh/test_mesh1d.py`` (S3.10).

    def test_mesh1d_bc_explicit(self):
        mesh = Mesher(StructuredGeometry.slab(
            (0.0, 1.0, 2.0), (0, 1), left=MirrorDeck(), right=0,
        )).partition(CellsByCount.uniform_width(1)).mesh
        assert dict(mesh.face_roles) == {"xmin": MirrorDeck(), "xmax": 0}
        assert mesh.outer_role == 0
        assert mesh.boundary_points == (0.0, 2.0)

    def test_mesh1d_bc_frozen(self):
        mesh = Mesher(StructuredGeometry.slab(
            (0.0, 1.0, 2.0), (0, 1), left=MirrorDeck(), right=MirrorDeck(),
        )).partition(CellsByCount.uniform_width(1)).mesh
        with pytest.raises(AttributeError):
            mesh.face_roles = {"xmin": 0, "xmax": 0}  # type: ignore[misc]  # the frozen field is the subject

    def test_mesh1d_bc_invalid_type_raises(self):
        with pytest.raises(TypeError, match="a face role is a MirrorDeck, a TranslationDeck or an int tag"):
            Mesh1D(
                coord=CoordSystem.CARTESIAN, edges=np.array([0.0, 1.0]),
                volumes=np.array([1.0]), region_ids=np.array([0]), region_materials=(0,),
                face_roles={"xmin": "vacuum", "xmax": 0},  # type: ignore[dict-item]  # a refusal input
            )

    def test_mesh1d_bc_accepts_a_deck_beside_a_tag(self):
        """The declaration channel: a deck and a tag are both legal face roles.

        The positive leg of the row above (``vv`` anti-pattern #11). A prescribed
        inflow, whose content is a function, has no role spelling: the problem
        binds it to a tag (a response is not the geometry's).
        """
        mesh = Mesher(StructuredGeometry.slab(
            (0.0, 1.0), (0,), left=7, right=MirrorDeck(),
        )).partition(CellsByCount.uniform_width(1)).mesh
        assert mesh.face_roles["xmin"] == 7 and mesh.face_roles["xmax"] == MirrorDeck()

    # ── Mesh2D face roles ────────────────────────────────────────────

    def test_mesh2d_bc_explicit(self):
        mesh = Mesh2D(
            edges_x=[0, 1], edges_y=[0, 1], mat_map=np.array([[0]]),
            face_roles={
                "xmin": MirrorDeck(),
                "xmax": 0,
                "ymin": MirrorDeck(),
                "ymax": 0,
            },
        )
        assert mesh.face_roles["xmin"] == MirrorDeck()
        assert mesh.face_roles["xmax"] == 0

    def test_mesh2d_bc_invalid_type_raises(self):
        with pytest.raises(TypeError, match=r"face_roles\['xmin'\]: a face role is a MirrorDeck"):
            Mesh2D(
                edges_x=[0, 1], edges_y=[0, 1], mat_map=np.array([[0]]),
                face_roles={
                    "xmin": "reflective",  # type: ignore[dict-item]  # the refusal input
                    "xmax": MirrorDeck(),
                    "ymin": MirrorDeck(),
                    "ymax": MirrorDeck(),
                },
            )

    def test_mesh2d_bc_accepts_a_tag_on_every_face(self):
        """The 2-D arm of the declaration channel — all four faces.

        Parameterised over every face because the four roles were validated
        by a hand-written loop until the channel landed; a widening that
        reached only ``xmin`` would pass a single-face row.
        """
        for face in ("xmin", "xmax", "ymin", "ymax"):
            face_roles: dict[str, object] = dict.fromkeys(("xmin", "xmax", "ymin", "ymax"), MirrorDeck())
            face_roles[face] = 5
            mesh = Mesh2D(
                edges_x=[0, 1], edges_y=[0, 1], mat_map=np.array([[0]]),
                face_roles=face_roles,
            )
            assert mesh.face_roles[face] == 5, face
