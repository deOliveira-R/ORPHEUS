r"""The inner face of a hollow cylinder or sphere is a TAG (D5), on the geometry, the mesh, the axis and the adapter
(wave 2 of the layer DAG, WF2.3 slice S1; ``.claude/plans/wave2_boundary_carve.md`` sections 2.3, 6.0 R6 and 6.1
gate 5).

Before the carve the radial axis had one law slot, ``bc_outer``; the inner law of a hollow body was dropped, and the
adapter either refused any inner law that was not reflective (``NotImplementedError``, 15 of 15 constructible
non-reflective typed laws) or synthesized a reflective cavity for the one it admitted. A reflecting curved wall is a
response (a tag bound to a specular wall), never a deck, so the axis carries the inner tag as DATA
(``RadialAxisMesh.role_inner``), the adapter round-trips it, and nothing is synthesized. The radial axis still has no
inner-surface TRACE (``SCOPE-BOUNDARY[guard]``, #511): the tag is carried, not realized, and is not in
``face_labels``.

Claims (layer: model; kind THEOREM for the round trip and the refusals, STRUCTURE for the synthesis): each refusal
row sits beside the construction that must succeed, so a red cannot be a constructor that refuses everything and
cannot be the unexpected keyword of the base tip hiding behind a ``pytest.raises(TypeError)``.

First reds at the base tip ``dec1f9ff``: ``Mesher`` on a hollow cylinder with ``inner=1`` raises ``TypeError`` from
``parse_boundary_law`` (an int is not a ``BC``); ``RadialAxisMesh`` has no ``role_outer`` or ``role_inner``;
``legacy_adapter.py:171`` returns ``BC("reflective")`` for the hollow inner face and ``:44`` refuses any other inner
law; ``_refuse_an_inner_law_the_radial_axis_drops`` exists.

Mutations, for the reviewers: the adapter returning ``role_outer`` for both faces (the round trip of ``inner=7``);
the adapter returning the string ``"reflective"`` or ``1`` for the inner face (the synthesis row and the tag-7 row);
``RadialAxisMesh`` accepting ``role_inner=None`` on a hollow edge array (the refusal row); ``admitted_decks`` not
consulted for ``role_inner`` (the deck row); ``face_labels`` listing the inner role (the inventory row).
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.geometry.deck import DeckError, MirrorDeck, TranslationDeck
from orpheus.mesh import (
    AxisCoord,
    AxisMesh,
    CellsByCount,
    Mesher,
    RadialAxisMesh,
    axes_from_legacy_mesh,
    face_labels,
    legacy_adapter,
    legacy_mesh_from_axes,
)
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/mesh/test_radial_inner_role.py"
_GEOMETRY_ROLES = "tests/gates/geometry/test_face_roles.py"
_ROOT = Path(__file__).resolve().parents[3]

_HOLLOW = np.array([0.5, 1.0, 1.5, 2.0])
_SOLID = np.array([0.0, 1.0, 1.5, 2.0])
_RADIAL = {
    "cylinder": (StructuredGeometry.cylinder, AxisCoord.RADIAL_CYLINDRICAL),
    "sphere": (StructuredGeometry.sphere, AxisCoord.RADIAL_SPHERICAL),
}
_REFUSED = (ValueError, TypeError)


def _axis(edges: np.ndarray, coord: AxisCoord, **roles: object) -> RadialAxisMesh:
    return RadialAxisMesh(edges=edges, coord=coord, **roles)  # type: ignore[arg-type]


def _mesh(kind: str, hollow: bool, **roles: int) -> Any:
    build, _ = _RADIAL[kind]
    breakpoints = (0.5, 1.0, 2.0) if hollow else (0.0, 1.0, 2.0)
    geometry = build(breakpoints, (0, 1), **roles)
    return Mesher(geometry).partition(CellsByCount.uniform_width(2)).mesh


@pytest.mark.rests_on(_GEOMETRY_ROLES + "::TestACurvedFaceCarriesNoDeck::test_a_hollow_body_either_face")
@pytest.mark.parametrize("kind", list(_RADIAL))
class TestTheHollowInnerFaceIsATag:
    def test_the_geometry_and_the_mesher_carry_it(self, kind: str) -> None:
        build, _ = _RADIAL[kind]
        geometry = build((0.5, 1.0, 2.0), (0, 1), outer=0, inner=1)
        require(geometry.face_roles == (1, 0), f"inner first: {geometry.face_roles!r}")
        mesh = Mesher(geometry).partition(CellsByCount.uniform_width(2)).mesh
        require(dict(mesh.face_roles) == {"xmin": 1, "xmax": 0}, f"{dict(mesh.face_roles)}")
        require(tuple(mesh.face_roles) == ("xmin", "xmax"), "in inventory order")

    def test_the_axis_carries_it_and_the_trace_inventory_does_not(self, kind: str) -> None:
        _, coord = _RADIAL[kind]
        hollow = _axis(_HOLLOW, coord, role_outer=0, role_inner=1)
        require(hollow.role_outer == 0 and hollow.role_inner == 1, f"{hollow.role_outer!r}, {hollow.role_inner!r}")
        solid = _axis(_SOLID, coord, role_outer=0)
        require(solid.role_inner is None, "a solid axis has no inner face")
        require(hollow.roles["outer"] == 0 and solid.roles == {"outer": 0}, f"{hollow.roles}, {solid.roles}")
        labels = face_labels((hollow,))
        require(len(labels) == 1, f"the radial axis has no inner-surface trace (#511): {labels}")
        require(len(face_labels((solid,))) == 1, "and the solid one has the same single face")

    def test_the_round_trip_keeps_the_inner_tag(self, kind: str) -> None:
        for inner in (1, 7, -2):
            mesh = _mesh(kind, hollow=True, outer=0, inner=inner)
            (axis,) = axes_from_legacy_mesh(mesh)
            require(isinstance(axis, RadialAxisMesh) and axis.role_inner == inner and axis.role_outer == 0, f"inner={inner}: {axis!r}")
            back = legacy_mesh_from_axes((axis,))
            require(dict(back.face_roles) == {"xmin": inner, "xmax": 0}, f"inner={inner}: {dict(back.face_roles)}")
            require(np.array_equal(back.edges, mesh.edges), "and the same cells")

    def test_the_adapter_synthesizes_no_cavity(self, kind: str) -> None:
        _, coord = _RADIAL[kind]
        built = legacy_mesh_from_axes((_axis(_HOLLOW, coord, role_outer=3, role_inner=5),))
        require(built.face_roles["xmin"] == 5, f"the inner role is the axis's, got {built.face_roles['xmin']!r}")
        solid = legacy_mesh_from_axes((_axis(_SOLID, coord, role_outer=3),))
        require(tuple(solid.face_roles) == ("xmax",) and solid.face_roles["xmax"] == 3, f"{dict(solid.face_roles)}")

    def test_a_solid_body_round_trips_without_an_inner_role(self, kind: str) -> None:
        mesh = _mesh(kind, hollow=False, outer=4)
        (axis,) = axes_from_legacy_mesh(mesh)
        require(axis.role_inner is None and axis.role_outer == 4, f"{axis!r}")
        require(dict(legacy_mesh_from_axes((axis,)).face_roles) == {"xmax": 4}, "and back")


@pytest.mark.rests_on(_HERE + "::TestTheHollowInnerFaceIsATag::test_the_axis_carries_it_and_the_trace_inventory_does_not")
@pytest.mark.parametrize("kind", list(_RADIAL))
class TestTheInnerRoleIsRequiredExactlyWhenHollow:
    def test_a_hollow_axis_without_an_inner_role_is_refused(self, kind: str) -> None:
        _, coord = _RADIAL[kind]
        _axis(_HOLLOW, coord, role_outer=0, role_inner=1)  # the control builds
        with pytest.raises(_REFUSED, match="role_inner"):
            _axis(_HOLLOW, coord, role_outer=0)
        with pytest.raises(_REFUSED, match="role_inner"):
            _axis(_HOLLOW, coord, role_outer=0, role_inner=None)

    def test_a_solid_axis_with_an_inner_role_is_refused(self, kind: str) -> None:
        _, coord = _RADIAL[kind]
        _axis(_SOLID, coord, role_outer=0)  # the control builds
        with pytest.raises(_REFUSED, match="role_inner"):
            _axis(_SOLID, coord, role_outer=0, role_inner=1)

    @pytest.mark.parametrize("deck", [MirrorDeck(), TranslationDeck()], ids=["mirror", "translation"])
    def test_no_radial_face_is_a_deck(self, kind: str, deck: Any) -> None:
        _, coord = _RADIAL[kind]
        _axis(_HOLLOW, coord, role_outer=0, role_inner=1)
        with pytest.raises(DeckError, match=r"N\(G_c\)/G_c"):
            _axis(_HOLLOW, coord, role_outer=0, role_inner=deck)
        with pytest.raises(DeckError, match=r"N\(G_c\)/G_c"):
            _axis(_HOLLOW, coord, role_outer=deck, role_inner=1)
        with pytest.raises(DeckError, match=r"N\(G_c\)/G_c"):
            _axis(_SOLID, coord, role_outer=deck)

    def test_a_role_is_a_tag_or_nothing(self, kind: str) -> None:
        _, coord = _RADIAL[kind]
        _axis(_SOLID, coord, role_outer=0)  # the control builds
        for bad in (True, 1.5, "vacuum", None):
            with pytest.raises(TypeError, match="role_outer"):
                _axis(_SOLID, coord, role_outer=bad)


@pytest.mark.rests_on(_HERE + "::TestTheHollowInnerFaceIsATag::test_the_axis_carries_it_and_the_trace_inventory_does_not")
class TestTheCartesianAxis:
    def test_it_carries_roles_under_the_new_names(self) -> None:
        axis = AxisMesh(edges=np.array([0.0, 1.0, 2.0]), role_low=MirrorDeck(), role_high=3)
        require(axis.role_low == MirrorDeck() and axis.role_high == 3, f"{axis!r}")
        require(axis.roles == {"min": MirrorDeck(), "max": 3}, f"{axis.roles}")
        renamed = AxisMesh(edges=np.array([0.0, 1.0]), role_low=1, role_high=2, label_low="left", label_high="right")
        require(renamed.roles == {"left": 1, "right": 2}, "the labels key the roles")
        with pytest.raises(TypeError, match="role_low"):
            AxisMesh(edges=np.array([0.0, 1.0]), role_low=True, role_high=2)
        with pytest.raises(TypeError, match="role_high"):
            AxisMesh(edges=np.array([0.0, 1.0]), role_low=1, role_high=None)  # type: ignore[call-arg]

    def test_the_pairing_closes_on_the_axis(self) -> None:
        """A translation on one face needs a translation on the other; the control pairs two."""
        edges = np.array([0.0, 1.0])
        built = AxisMesh(edges=edges, role_low=TranslationDeck(), role_high=TranslationDeck())
        require(built.role_low == TranslationDeck(), "control: a translation pair builds")
        with pytest.raises(DeckError):
            AxisMesh(edges=edges, role_low=TranslationDeck(), role_high=0)
        with pytest.raises(DeckError, match="xmin.*TranslationDeck"):
            AxisMesh(edges=edges, role_low=TranslationDeck(), role_high=0, label_low="left", label_high="right")

    def test_the_old_names_are_gone(self) -> None:
        axis = AxisMesh(edges=np.array([0.0, 1.0]), role_low=1, role_high=2)
        gone = [n for n in ("bc_low", "bc_high", "bc") if hasattr(axis, n)]
        require(not gone, f"AxisMesh still has {gone}")
        radial = _axis(_SOLID, AxisCoord.RADIAL_CYLINDRICAL, role_outer=0)
        require(not any(hasattr(radial, n) for n in ("bc_outer", "bc")), "RadialAxisMesh still has bc_outer or bc")
        with pytest.raises(TypeError, match="bc_low"):
            AxisMesh(edges=np.array([0.0, 1.0]), role_low=1, role_high=2, bc_low=1)  # type: ignore[call-arg]

    def test_a_slab_mesh_round_trips_its_roles(self) -> None:
        geometry = StructuredGeometry.slab((0.0, 1.0, 2.0), (0, 1), left=MirrorDeck(), right=6)
        mesh = Mesher(geometry).partition(CellsByCount.uniform_width(2)).mesh
        (axis,) = axes_from_legacy_mesh(mesh)
        require(axis.roles == {"min": MirrorDeck(), "max": 6}, f"{axis.roles}")
        require(dict(legacy_mesh_from_axes((axis,)).face_roles) == {"xmin": MirrorDeck(), "xmax": 6}, "and back")


class TestTheAdapterHoldsNoCavityAndNoGuard:
    """STRUCTURE: the synthetic reflective cavity and its guard are gone from ``legacy_adapter.py``."""

    def test_the_guard_and_the_literal_are_gone(self) -> None:
        require(not hasattr(legacy_adapter, "_refuse_an_inner_law_the_radial_axis_drops"), "the guard is retired")
        source = (_ROOT / "orpheus" / "mesh" / "legacy_adapter.py").read_text()
        tree = ast.parse(source)
        strings = [n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str)]
        # Docstrings may say what was retired; code strings may not name a reflective cavity.
        code_strings = [s for s in strings if len(s) < 40]
        require("reflective" not in code_strings, "no code literal 'reflective' (the synthesized cavity)")
        names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)} | {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
        require("BC" not in names, "and no BC")
        require(not hasattr(legacy_adapter, "_face_laws_of_axis"), "the helper is renamed with its content")
