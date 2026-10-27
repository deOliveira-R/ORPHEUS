r"""The face roles of a structured geometry, and the retirement of the response vocabulary from geometry and mesh
(wave 2 of the layer DAG, WF2.3 slice S1; ``.claude/plans/wave2_boundary_carve.md`` sections 2.2, 6.1 and 8).

Each face of a ``StructuredGeometry`` carries a ROLE, ``Deck | int``: a deck says the domain continues across the
face (a mirror, a translation), an integer is a boundary tag whose condition a problem binds later. Geometry names no
response: no law, no ``BC``. The claims, each with its kind and layer (``vv-principles``):

* gate 1, ``test_a_curved_face_carries_no_deck``: THEOREM (``W = N(G_c)/G_c`` is trivial on a 1-D radial chart, so a
  reflecting curved wall is a response, a tagged face with a specular wall condition) plus the translation pairing;
  each refusal row sits beside the construction that must succeed (a positive control per row, so a red cannot be
  a constructor that refuses everything);
* gate 2, ``test_geometry_and_mesh_name_no_response``: a STRUCTURE gate (an ``ast`` scan over the geometry and mesh
  packages, validated by a planted import of each shape) and law L-F (a response-shaped value as a role is a
  ``TypeError``);
* gate 3, ``test_face_roles_are_content`` (law L-J): THEOREM of the content-identity convention;
* gate 4, ``test_tags_are_ints_like_material_ids``: REFERENCE (the material-id twin rule);
* gate 7, ``test_the_presets``: RECORD of two declared models (the Wigner-Seitz cell binds tag 0 to a white wall
  elsewhere; the PWR half cell is a lattice's symmetry);
* the retirement rows (``test_the_retired_spellings_are_gone``) and the two delegation rows
  (``test_the_mirror_twins_delegate``): STRUCTURE gates for the deletions and the single-sourced axis-letter to
  mirror spelling (X4).

Gates 5, 6 and 8 live in ``tests/gates/mesh/test_radial_inner_role.py``, ``tests/gates/mesh/test_face_roles.py``
and the content-identity rosters. The module reads every NEW attribute at run time, so each row's red is its own.

First reds at the base tip ``dec1f9ff`` (each row's own reason): ``StructuredGeometry.slab(... left=MirrorDeck())``
raises ``TypeError("... must be a BC tag or a BoundaryTraceLaw")``, never ``DeckError`` (gate 1); six files under
``orpheus/geometry`` and ``orpheus/mesh`` import ``orpheus.geometry.boundary`` (gate 2, the population); the keyword
``face_roles`` does not exist (gates 2 L-F, 3, 4, 7); ``orpheus.geometry.BC`` exists and ``parse_boundary_law(s)``
are defined (retirement); ``_reflections('q')`` says ``axis must be x/y/z`` and both twins keep their own bodies
(delegation).

Mutations, for the reviewers (each on a COPY, never by ``git checkout``):

* gate 1: ``admitted_decks`` returning the Cartesian set for a radial first axis (the radial rows); a geometry that
  skips ``check_pairings_close`` (the ``slab(left=TranslationDeck(), right=0)`` row); a geometry that checks the
  decks of the OUTER face only (the hollow ``inner=`` rows);
* gate 2: re-adding ``from orpheus.geometry.boundary import BC`` under ``if TYPE_CHECKING:`` in any mesh module (the
  scan sees annotation-only edges); a ``parse_face_role`` that falls back to ``BC``-like objects having a ``kind``;
* gate 3: a geometry whose digest ignores ``face_roles`` (the digest legs); a role compared by ``type`` alone (the
  deck-versus-tag leg);
* gate 4: a ``parse_face_role`` that accepts ``bool`` (``True`` becomes tag 1); a ``tags`` property returning the
  decks too;
* gate 7: ``wigner_seitz_pin_cell`` left at ``outer=BC("white")`` is a ``TypeError`` before the row reads it.
"""

from __future__ import annotations

import ast
import importlib.util
import inspect
import re
import types
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import orpheus.geometry as geometry_package
from orpheus.geometry import CoordSystem, StructuredGeometry, structured_geometry, symmetry
from orpheus.geometry.deck import DeckError, MirrorDeck, TranslationDeck
from orpheus.transport.boundary_conditions import WallCondition
from orpheus.values.content import content_digest
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/geometry/test_face_roles.py"
_DECK = "tests/gates/geometry/test_deck.py"
_ROOT = Path(__file__).resolve().parents[3]

_SPH, _CYL, _CART = CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL, CoordSystem.CARTESIAN
_SOLID = (0.0, 1.0, 2.0)
_HOLLOW = (0.5, 1.0, 2.0)
_IDS = (0, 1)
_REFUSAL = r"N\(G_c\)/G_c"


class _LawShaped:
    """A response-shaped stand-in (a ``kind`` and ``params``, the old ``BC`` surface): not a deck and not an int."""

    kind = "vacuum"
    params: dict[str, float] = {}


# ── gate 1 ───────────────────────────────────────────────────────────────────


@pytest.mark.rests_on(_DECK + "::test_a_curved_face_admits_no_deck", _DECK + "::test_pairings_close")
class TestACurvedFaceCarriesNoDeck:
    """Gate 1 (G2): a deck on a radial face is a ``DeckError`` naming ``W = N(G_c)/G_c``; the tag beside it builds.

    Layer: model; kind THEOREM. The radial first axis admits no deck (``admitted_decks(coord, 0)`` is empty), a
    translation requires a translation on the opposite face. Red at base: every deck is a ``TypeError`` from
    ``parse_boundary_law`` ("must be a BC tag or a BoundaryTraceLaw"), never a ``DeckError``.
    """

    @pytest.mark.parametrize("deck", [MirrorDeck(), TranslationDeck()], ids=["mirror", "translation"])
    @pytest.mark.parametrize("coord", [_CYL, _SPH], ids=["cylinder", "sphere"])
    def test_a_solid_outer_face(self, coord: CoordSystem, deck: Any) -> None:
        build = StructuredGeometry.cylinder if coord is _CYL else StructuredGeometry.sphere
        built = build(_SOLID, _IDS, outer=0)
        require(built.face_roles == (0,), f"the control: a solid {coord.name.lower()} with outer tag 0 builds, got {built.face_roles!r}")
        with pytest.raises(DeckError, match=_REFUSAL):
            build(_SOLID, _IDS, outer=deck)

    @pytest.mark.parametrize("deck", [MirrorDeck(), TranslationDeck()], ids=["mirror", "translation"])
    @pytest.mark.parametrize("coord", [_CYL, _SPH], ids=["cylinder", "sphere"])
    @pytest.mark.parametrize("where", ["inner", "outer"])
    def test_a_hollow_body_either_face(self, coord: CoordSystem, deck: Any, where: str) -> None:
        """The inner face of a hollow body is a tag (D5) and neither of its faces is a deck."""
        build = StructuredGeometry.cylinder if coord is _CYL else StructuredGeometry.sphere
        built = build(_HOLLOW, _IDS, outer=0, inner=1)
        require(built.face_roles == (1, 0), f"the control: tags on both faces build, inner first, got {built.face_roles!r}")
        roles = {"inner": 1, "outer": 0, where: deck}
        with pytest.raises(DeckError, match=_REFUSAL):
            build(_HOLLOW, _IDS, outer=roles["outer"], inner=roles["inner"])

    def test_the_slab_admits_both_decks_and_closes_the_pairing(self) -> None:
        """A Cartesian face admits both decks; a translation needs a translation on the opposite face."""
        mirror = StructuredGeometry.slab((0.0, 1.0, 2.0), _IDS, left=MirrorDeck(), right=0)
        require(mirror.face_roles == (MirrorDeck(), 0), f"mirror on the left, tag on the right: {mirror.face_roles!r}")
        torus = StructuredGeometry.slab((0.0, 1.0, 2.0), _IDS, left=TranslationDeck(), right=TranslationDeck())
        require(torus.face_roles == (TranslationDeck(), TranslationDeck()), f"the periodic slab: {torus.face_roles!r}")
        for left, right in ((TranslationDeck(), 0), (0, TranslationDeck()), (TranslationDeck(), MirrorDeck())):
            with pytest.raises(DeckError, match="TranslationDeck"):
                StructuredGeometry.slab((0.0, 1.0, 2.0), _IDS, left=left, right=right)


# ── gate 2 ───────────────────────────────────────────────────────────────────


def _boundary_imports(root: Path, packages: tuple[str, ...], law_names: frozenset[str]) -> tuple[list[str], int]:
    """The imports of the response vocabulary (the ``orpheus.geometry.boundary`` package, or a law name taken from the
    ``orpheus.geometry`` facade) in the modules under ``root/orpheus/<package>``, outside ``geometry/boundary/``.

    Returns ``(findings, number of files scanned)``. An ``ast`` walk, so a function-local import, an import under
    ``TYPE_CHECKING``, an aliased import and a relative import are all seen; a string passed to ``importlib`` is not
    (the stated exclusion).
    """
    findings: list[str] = []
    scanned = 0
    for package in packages:
        for path in sorted((root / "orpheus" / package).rglob("*.py")):
            relative = path.relative_to(root)
            if relative.parts[:3] == ("orpheus", "geometry", "boundary"):
                continue
            scanned += 1
            parts = list(relative.with_suffix("").parts)
            module_package = parts[:-1] if parts[-1] != "__init__" else parts[:-1]
            for node in ast.walk(ast.parse(path.read_text())):
                if isinstance(node, ast.Import):
                    findings += [
                        f"{relative}:{node.lineno} import {a.name}" for a in node.names
                        if a.name == "orpheus.geometry.boundary" or a.name.startswith("orpheus.geometry.boundary.")
                    ]
                elif isinstance(node, ast.ImportFrom):
                    base = node.module or ""
                    if node.level:
                        anchor = module_package[: len(module_package) - (node.level - 1)]
                        base = ".".join(anchor + ([node.module] if node.module else []))
                    imported = {a.name for a in node.names}
                    if base == "orpheus.geometry.boundary" or base.startswith("orpheus.geometry.boundary."):
                        findings.append(f"{relative}:{node.lineno} from {base} import {sorted(imported)}")
                    elif base == "orpheus.geometry" and (imported & (law_names | {"boundary"})):
                        findings.append(f"{relative}:{node.lineno} from {base} import {sorted(imported & (law_names | {'boundary'}))}")
    return findings, scanned


def _the_response_vocabulary() -> frozenset[str]:
    """The names the facade must not export: the package is gone (WF2.4), so the module path is the whole rule and ``BC`` the one name."""
    return frozenset({"BC"})


_PLANTED = {
    "absolute from-import": ("geometry", "from orpheus.geometry.boundary import BC\n"),
    "plain import": ("geometry", "import orpheus.geometry.boundary.reflective\n"),
    "a submodule": ("mesh", "from orpheus.geometry.boundary._factors import x\n"),
    "the facade, a law name": ("mesh", "from orpheus.geometry import BC, CoordSystem\n"),
    "the facade, the package": ("mesh", "from orpheus.geometry import boundary\n"),
    "relative, geometry": ("geometry", "from .boundary import BC\n"),
    "relative, up from mesh": ("mesh", "from ..geometry.boundary import BoundaryTraceLaw\n"),
    "under TYPE_CHECKING": ("mesh", "from typing import TYPE_CHECKING\nif TYPE_CHECKING:\n    from orpheus.geometry.boundary import BC\n"),
    "function-local": ("geometry", "def f():\n    from orpheus.geometry.boundary import BoundaryTraceLaw\n"),
}


class TestGeometryAndMeshNameNoResponse:
    """Gate 2 (law L-F): neither package imports the response vocabulary, and a response is not a role.

    Layer: architecture (``item2_carve.md`` forbidden edges 1 and 2); kind STRUCTURE for the scan, THEOREM for L-F.
    The scan is validated first, by a planted import of each of the nine shapes it must see and a tree of clean imports (and a
    file inside ``boundary/``) it must not flag.
    Population at base: the findings below name 10 import statements in 7 modules (``geometry/__init__``,
    ``geometry/structured_geometry`` x3, ``mesh/axis`` x2, ``mesh/face_laws``, ``mesh/factories``,
    ``mesh/legacy_adapter``, ``mesh/structured``). Exclusions: ``orpheus/geometry/boundary/`` itself (deleted at
    WF2.4), string-spelled imports.
    """

    @pytest.mark.parametrize("label", list(_PLANTED))
    def test_the_scan_sees_every_planted_shape(self, tmp_path: Path, label: str) -> None:
        package, source = _PLANTED[label]
        (tmp_path / "orpheus" / package).mkdir(parents=True)
        (tmp_path / "orpheus" / package / "planted.py").write_text(source)
        findings, scanned = _boundary_imports(tmp_path, ("geometry", "mesh"), frozenset({"BC", "BoundaryTraceLaw"}))
        require(scanned == 1 and len(findings) == 1, f"{label}: scanned {scanned}, found {findings}")

    def test_the_scan_ignores_what_it_must(self, tmp_path: Path) -> None:
        (tmp_path / "orpheus" / "geometry" / "boundary").mkdir(parents=True)
        (tmp_path / "orpheus" / "mesh").mkdir(parents=True)
        (tmp_path / "orpheus" / "geometry" / "boundary" / "inside.py").write_text("from .reflective import ReflectiveBoundary\nfrom orpheus.geometry.boundary import BC\n")
        (tmp_path / "orpheus" / "geometry" / "clean.py").write_text("from orpheus.geometry.deck import MirrorDeck\nfrom orpheus.geometry import CoordSystem\n")
        (tmp_path / "orpheus" / "mesh" / "clean.py").write_text("from orpheus.geometry.colouring import colouring_symmetry\nfrom .axis import AxisMesh\n")
        findings, scanned = _boundary_imports(tmp_path, ("geometry", "mesh"), frozenset({"BC"}))
        require(scanned == 2 and findings == [], f"scanned {scanned}, found {findings}")

    def test_no_import_of_the_response_vocabulary(self) -> None:
        findings, scanned = _boundary_imports(_ROOT, ("geometry", "mesh"), _the_response_vocabulary())
        require(scanned >= 20, f"the scan reads the two packages: {scanned} files")
        require(findings == [], "geometry and mesh name no response; imports found:\n  " + "\n  ".join(findings))

    def test_the_facade_exports_no_law(self) -> None:
        # Module objects are left out of the vars() leg: Python binds the submodule ``boundary`` on the package whenever
        # any module imports it (this file does), which says nothing about what the facade exports.
        bound = {n for n, v in vars(geometry_package).items() if not n.startswith("_") and not isinstance(v, types.ModuleType)}
        exported = set(geometry_package.__all__) | bound
        leaked = sorted((exported & _the_response_vocabulary()) | ({"boundary"} & set(geometry_package.__all__)))
        require(not leaked, f"orpheus.geometry exports the response vocabulary: {leaked}")

    def test_the_boundary_package_is_gone(self) -> None:
        """WF2.4 row: the retired package has no module spec (it was a strict xfail until the deletion)."""
        require(importlib.util.find_spec("orpheus.geometry.boundary") is None, "orpheus.geometry.boundary still exists")

    @pytest.mark.parametrize(
        "stand_in",
        [_LawShaped(), _LawShaped, "vacuum", None, 1.0, WallCondition.vacuum(), WallCondition.specular(1.0)],
        ids=["law-shaped", "law-class", "string", "none", "float", "wall-vacuum", "wall-specular"],
    )
    def test_a_response_is_not_a_role(self, stand_in: object) -> None:
        """L-F. The control sits first: the same call with a tag builds. Red at base: the keyword does not exist.

        The rows are duck-typed stand-ins of the retired law surface (a ``kind`` and ``params``): the threat is a role parse
        that falls back to any object with a ``kind``. (The two rows that built the real ``VacuumInflow`` and
        ``ReflectiveBoundary`` retired with their classes in WF2.4; their successors are the two ``WallCondition`` rows, the response
        type that replaced the law family: a response is not a role.)
        """
        built = StructuredGeometry(coord=_CART, breakpoints=(0.0, 1.0), mat_ids=(0,), face_roles=(0, 1))
        require(built.face_roles == (0, 1), f"the control builds: {built.face_roles!r}")
        with pytest.raises(TypeError, match=re.escape(repr(stand_in))):
            StructuredGeometry(coord=_CART, breakpoints=(0.0, 1.0), mat_ids=(0,), face_roles=(stand_in, 1))  # type: ignore[call-arg]


# ── gate 3 ───────────────────────────────────────────────────────────────────


def _slab(left: object = 0, right: object = 1) -> StructuredGeometry:
    return StructuredGeometry.slab((0.0, 1.0, 2.0), _IDS, left=left, right=right)  # type: ignore[arg-type]


@pytest.mark.rests_on(_HERE + "::TestACurvedFaceCarriesNoDeck::test_the_slab_admits_both_decks_and_closes_the_pairing")
def test_face_roles_are_content() -> None:
    """Gate 3 (law L-J). Layer: model; kind THEOREM (content identity).

    A deck and a tag are different declarations; the faces' order is content for a geometry (inner first, a tuple);
    a role moves the digest and the equality. The controls: equal builds are equal and share a digest.
    """
    require(_slab() == _slab() and content_digest(_slab()) == content_digest(_slab()), "control: two builds are one value")
    require(hash(_slab()) == hash(_slab()), "control: and hash alike")
    variants = {
        "a mirror, not tag 0": _slab(left=MirrorDeck()),
        "tag 2, not tag 0": _slab(left=2),
        "swapped tags": _slab(left=1, right=0),
        "a translation pair": _slab(left=TranslationDeck(), right=TranslationDeck()),
        "a mirror pair": _slab(left=MirrorDeck(), right=MirrorDeck()),
        "mirror on the right": _slab(left=0, right=MirrorDeck()),
    }
    base = _slab()
    for label, other in variants.items():
        require(other != base, f"{label}: not equal to the base")
        require(content_digest(other) != content_digest(base), f"{label}: the digest moves")
    require(len({content_digest(v) for v in variants.values()}) == len(variants), "the six variants have six digests")
    require(_slab(left=MirrorDeck()) != _slab(left=1), "a deck is not a tag")
    require(content_digest(_slab(left=MirrorDeck())) != content_digest(_slab(left=1)), "nor in digest")
    require(_slab(left=np.int64(3)) == _slab(left=3), "a numpy integer tag is the int tag")
    sphere = StructuredGeometry.sphere(_SOLID, _IDS, outer=0)
    require(sphere != StructuredGeometry.sphere(_SOLID, _IDS, outer=1), "a radial outer tag is content")


# ── gate 4 ───────────────────────────────────────────────────────────────────


@pytest.mark.rests_on(_DECK + "::test_parse_face_role", _HERE + "::test_face_roles_are_content")
def test_tags_are_ints_like_material_ids() -> None:
    """Gate 4. Layer: model; kind REFERENCE (the material-id twin rule, ``orpheus/data/materials.py``).

    A ``bool`` role is refused; a numpy integer is the int; negative tags are tags; ``tags`` is the set of tags the
    faces carry (decks excluded). Red at base: no ``face_roles`` and no ``tags``.
    """
    require(_slab(left=3, right=3).tags == frozenset({3}), "one tag on both faces")
    require(_slab(left=MirrorDeck(), right=4).tags == frozenset({4}), "a deck carries no tag")
    require(_slab(left=MirrorDeck(), right=MirrorDeck()).tags == frozenset(), "two decks, no tags")
    require(_slab(left=-2, right=np.int32(5)).tags == frozenset({-2, 5}), "negative and numpy tags")
    require(all(type(t) is int for t in _slab(left=np.int64(1), right=2).tags), "tags are plain ints")
    for bad in (True, False, np.True_):
        with pytest.raises(TypeError, match="face_roles"):
            _slab(left=bad)
    with pytest.raises(TypeError, match="face_roles"):
        _slab(right=1.5)
    require(StructuredGeometry.sphere(_HOLLOW, _IDS, outer=2, inner=7).tags == frozenset({2, 7}), "a hollow body: both tags")


# ── gate 7 ───────────────────────────────────────────────────────────────────


@pytest.mark.rests_on(_HERE + "::test_face_roles_are_content")
def test_the_presets() -> None:
    """Gate 7. Layer: model; kind RECORD of two declarations.

    The Wigner-Seitz cell declares tag 0 (the model binds it to a white wall, which a geometry cannot hold); the PWR
    half cell is a lattice's symmetry planes on both faces.
    """
    wigner_seitz = StructuredGeometry.wigner_seitz_pin_cell()
    require(wigner_seitz.face_roles == (0,), f"{wigner_seitz.face_roles!r}")
    half_cell = StructuredGeometry.pwr_slab_half_cell()
    require(half_cell.face_roles == (MirrorDeck(), MirrorDeck()), f"{half_cell.face_roles!r}")
    require(half_cell.tags == frozenset() and wigner_seitz.tags == frozenset({0}), "the tags each preset carries")
    require("WallCondition.white" in (StructuredGeometry.wigner_seitz_pin_cell.__doc__ or ""), "the docstring says the model binds the wall")


# ── the retirements and the delegations ─────────────────────────────────────


@pytest.mark.rests_on(_HERE + "::TestGeometryAndMeshNameNoResponse::test_no_import_of_the_response_vocabulary")
def test_the_retired_spellings_are_gone() -> None:
    """No dual spelling: the old geometry API stops working in S1. Layer: architecture; kind STRUCTURE.

    ``boundaries`` is no field, the keyword and the property are gone, ``parse_boundary_law(s)`` are undefined and
    the facade's ``BC`` is not exported. Each row names the old spelling it rejects; the controls are the new ones.
    """
    import dataclasses

    fields = [f.name for f in dataclasses.fields(StructuredGeometry)]
    require("face_roles" in fields and "boundaries" not in fields, f"the fields are {fields}")
    structured = structured_geometry
    for name in ("parse_boundary_law", "parse_boundary_laws"):
        require(not hasattr(structured, name), f"{name} is retired")
        require(name not in getattr(structured, "__all__", ()), f"{name} is not exported")
    require(not hasattr(geometry_package, "BC") and "BC" not in geometry_package.__all__, "the facade drops BC")
    with pytest.raises(TypeError, match="boundaries"):
        StructuredGeometry(coord=_CART, breakpoints=(0.0, 1.0), mat_ids=(0,), boundaries=(0, 1))  # type: ignore[call-arg]


def _delegation_body(module: types.ModuleType, name: str) -> list[ast.stmt]:
    tree = ast.parse(Path(str(module.__file__)).read_text())
    (fn,) = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == name]
    body = fn.body
    if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) and isinstance(body[0].value.value, str):
        body = body[1:]
    return body


@pytest.mark.parametrize(
    ("module", "name"),
    [
        (symmetry, "_reflections"),
    ],
    ids=["symmetry"],
)
def test_the_mirror_twins_delegate(module: types.ModuleType, name: str) -> None:
    """X4: one axis-letter to mirror spelling. Layer: architecture; kind STRUCTURE plus a message fragment.

    The body of each twin, past its docstring, is ONE statement and its only call into the motion algebra is
    ``mirror_motion``: no ``AXIS_INDEX`` lookup, no ``np.zeros``, no ``RigidMotion.reflection`` of its own. A bad
    axis keeps the fragment "axis must be one of" (the pinned vocabulary). Red at base: both bodies are the
    original twins, and ``_reflections`` says "axis must be x/y/z".
    """
    body = _delegation_body(module, name)
    require(len(body) == 1 and isinstance(body[0], ast.Return), f"{name}: one return statement, got {[type(s).__name__ for s in body]}")
    names = {n.id for n in ast.walk(body[0]) if isinstance(n, ast.Name)}
    attrs = {n.attr for n in ast.walk(body[0]) if isinstance(n, ast.Attribute)}
    require("mirror_motion" in names | attrs, f"{name} calls mirror_motion: names {sorted(names)}, attributes {sorted(attrs)}")
    require(not ({"AXIS_INDEX", "zeros", "reflection", "eye"} & (names | attrs)), f"{name} re-spells the mirror: {sorted(names | attrs)}")
    twin = getattr(module, name)
    arguments = inspect.signature(twin).parameters
    with pytest.raises(ValueError, match="axis must be one of"):
        twin("q", 2) if len(arguments) == 2 else twin("q")
