r"""The reference reads its walls from geometry's decks and the specification's amplitudes (wave 2, S5; plan 2.9 and gates 3 of 6.5; decision 54).

Specified by the test-architect (2026-10-23, ``.claude/plans/wave2_boundary_carve.md`` 6.5). The characteristic reference's wall
table ``Walls`` stops reading a typed boundary law through its two factors and reads two plain things: a face role that is a
DECK (``MirrorDeck`` is a mirror returning everything along the reflected line; ``TranslationDeck`` a wrap re-entering at the
opposite wall) and, for a face role that is an integer TAG, the ``WallAmplitudes(specular, diffuse)`` the specification's
``exterior`` states for it. ``Wall`` COMPOSES ``WallAmplitudes`` (``Wall(breakpoint, amplitudes, partner)``, with ``specular`` and
``diffuse`` read-only properties of ``amplitudes``), so the bounds have ONE owner. ``specular_albedos`` and ``require_vacuum``
read the same pair ``(geometry, exterior)``. The three reference consumers (``Spectrum``, ``MomentSpace``, ``BasisSpace``)
gain an ``exterior`` field beside ``geometry``.

The gate is an equivalence between two spellings of one population, grounded twice, independently (``instrument-doctrine`` X4):

* **A hand rule per declaration** (``test_walls_of_reads_each_body_as_the_wall_its_physics_names``): a vacuum returns nothing, a
  mirror returns along the reflected line, a partial specular returns the stated fraction specularly, a white wall
  isotropically, a wrap re-enters at the opposite wall; the expected ``(breakpoint, specular, diffuse, partner)`` of each body is
  written from that rule by :func:`_expected` (the rule, not a table of the answers), never through ``Walls``.
* **The HEAD reading** (``test_walls_of_reads_what_head_read``): the walls ``Walls.of(geometry)`` returned at ``main`` ``5d04c296``
  for the same 19 bodies spelled as the retired string tags, captured once into ``reference_carry_head.json``
  (:mod:`tests.gates.derivations._capture_reference_carry_head`). The population (19 bodies: a 3-region slab under 8 wall pairs, a
  3-region sphere under 4, a 2-region cylinder under 3, a hollow sphere under 3, a hollow cylinder under 1) is
  :mod:`tests.gates.derivations._reference_carry_cases`; it spans every wall kind in every chart a reference serves, the
  mirror on a curved face (a tag bound to ``(1, 0)``: the carve's meaning change) and both orderings of a diffuse and a specular
  wall on one body.

Then the new types' own laws: ``Wall``'s composition and its two read-only properties; the R-8b re-baseline (the bounds
refusal is a ``ValueError`` from ``WallAmplitudes``, fragment "not a physical wall" kept, no longer a ``NotImplementedError``);
the combined wall, unchanged, ``NotImplementedError`` naming #650; ``specular_albedos(geometry, exterior, *, owner)`` and
``require_vacuum`` on the pair; the ``exterior`` field of the three consumers (read from the source by AST, so a red is the
field and not an import chain); and the Sood registry's finite cases carrying the one tag ``0`` that the consumers' ``exterior``
binds to vacuum (rule R1: equal declarations share a tag).

Declared red at the base tip (``de354496``), each on its own line: ``orpheus.derivations.continuous.characteristic`` imports at the
base (its ``walls.py`` still reads ``orpheus.geometry.boundary``), so the walls rows record ``TypeError: Walls.of() takes 2
positional arguments but 3 were given`` and the ``Wall`` rows the old four-float signature; ``common.reference_body`` cannot be
imported (``ImportError: cannot import name 'BC'``), so the albedo rows record that line; the consumer rows record the missing
field. The mutations the reviewers apply to the carved tree: (1) ``MirrorDeck`` read as vacuum (``Wall(k, WallAmplitudes(0, 0), k)``) reddens
the mirror bodies of both walls rows; (2) the wrap's partner left at the wall itself reddens the periodic body; (3) the diffuse
amplitude read into ``specular`` reddens every white body; (4) the amplitudes read from the WRONG tag (``exterior[0]`` always)
reddens every two-tag body (``slab-partial-vacuum``, ``hollow-sphere-partial-mirror``), which a one-tag-only fixture cannot see
(Mode 12: the one-tag bodies annihilate that degree of freedom, the two-tag ones carry it); (5) ``Wall`` re-given float fields
reddens the composition rows. The activation of each row is stated: no row nulls the heterogeneity (every body is 2 or 3 regions)
and the walls are read independent of the group count (a wall carries no energy index), which is the reference's declared
blindness here and is why the SOLVES live in ``test_reference_carry.py`` on a two-group body.
"""

from __future__ import annotations

import ast
import dataclasses
import json
from pathlib import Path

import pytest

from tests.gates._content_identity_helpers import require
from tests.gates.derivations import _reference_carry_cases as cases
from tests.gates.derivations import _reference_carry_spelling as spelling

pytestmark = pytest.mark.foundation

_HERE_DIR = Path(__file__).resolve().parent
_ROOT = _HERE_DIR.parents[2]
_FIXTURE = _HERE_DIR / "reference_carry_head.json"
_HERE = "tests/gates/derivations/test_walls_of_reads_decks_and_amplitudes.py::"
_LEAF = "tests/gates/specification/test_exterior.py::test_wall_amplitudes_a_stated_wall_constructs"
_BODY_IDS = [c[0] for c in cases.CASES]


def _head() -> dict:
    return json.loads(_FIXTURE.read_text())


def _expected(chart: str, breakpoints: tuple[float, ...], faces: tuple[cases.Decl, ...]) -> list[tuple[int, float, float, int]]:
    """The walls a body's declarations name, from the physics rule, inner first: ``(breakpoint, specular, diffuse, partner)``.

    A solid radial body has the one wall at its outer breakpoint ``n`` (its centre is an interior point); a slab and a hollow
    body have walls at ``0`` and ``n``. A wrap re-enters at the opposite wall; every other wall returns to itself.
    """
    n = len(breakpoints) - 1
    solid_radial = chart != "slab" and breakpoints[0] == 0.0
    ends = (n,) if solid_radial else (0, n)
    rows = []
    for end, decl in zip(ends, faces, strict=True):
        if decl[0] == "periodic":
            rows.append((end, 1.0, 0.0, n - end))
        else:
            specular, diffuse = cases.amplitudes(decl)
            rows.append((end, specular, diffuse, end))
    return rows


def _read(walls) -> list[tuple[int, float, float, int]]:
    return [(w.breakpoint, w.specular, w.diffuse, w.partner) for w in walls.walls]


def _walls_of(chart: str, breakpoints, faces):
    from orpheus.derivations.continuous.characteristic import Walls

    geometry, exterior = spelling.geometry_and_exterior(chart, breakpoints, faces)
    return Walls.of(geometry, exterior)


def test_the_head_fixture_is_the_population_it_claims() -> None:
    """The fixture holds exactly the 19 bodies of the table, each with its walls; the hand rule is run beside it.

    Positive control and population in one: the fixture and the table name the same bodies (a body added to one and not the other
    reds here), and the hand rule reproduces the fixture's walls for every body (so the fixture is read as what it is, the HEAD
    reading of the SAME declarations, before any carved code runs). Not red at base: it reads no production."""
    head = _head()
    require(sorted(k for k in head if not k.startswith("_")) == sorted(_BODY_IDS), "the fixture and the case table name different bodies")
    require(len(_BODY_IDS) == 19, f"the population is {len(_BODY_IDS)} bodies, not 19")
    for name, chart, bp, faces in cases.CASES:
        require(head[name]["walls"] == [list(r) for r in _expected(chart, bp, faces)], f"{name}: HEAD read {head[name]['walls']}")


@pytest.mark.rests_on(_LEAF, _HERE + "test_the_head_fixture_is_the_population_it_claims")
@pytest.mark.parametrize(("chart", "breakpoints", "faces"), [c[1:] for c in cases.CASES], ids=_BODY_IDS)
def test_walls_of_reads_each_body_as_the_wall_its_physics_names(chart, breakpoints, faces) -> None:
    """``Walls.of(geometry, exterior)`` of each body is the hand-ruled wall table: decks and amplitudes, nothing else.

    Equality of ``(breakpoint, specular, diffuse, partner)`` per wall, exact (an amplitude is the stated float, a partner an
    index). The two-tag bodies are the discriminators of Mode 12 (see the module docstring, mutation (4))."""
    got = _read(_walls_of(chart, breakpoints, faces))
    require(got == _expected(chart, breakpoints, faces), f"{got}")


@pytest.mark.rests_on(_HERE + "test_walls_of_reads_each_body_as_the_wall_its_physics_names")
@pytest.mark.parametrize(("name", "chart", "breakpoints", "faces"), cases.CASES, ids=_BODY_IDS)
def test_walls_of_reads_what_head_read(name, chart, breakpoints, faces) -> None:
    """The carved reading equals ``main`` ``5d04c296``'s reading of the same declarations, wall by wall, exactly.

    The independent second ground (the first is the hand rule above): a regression that moved one amplitude, one partner or one
    wall in any body is a red here and in the rule row, and a carve that moved BOTH the rule's reading and the HEAD reading the
    same way is impossible, the HEAD reading having been captured before the carve existed."""
    got = [list(r) for r in _read(_walls_of(chart, breakpoints, faces))]
    require(got == _head()[name]["walls"], f"{name}: {got} against HEAD {_head()[name]['walls']}")


def test_walls_of_takes_the_geometry_and_the_exterior_and_a_deck_only_body_needs_none() -> None:
    """The signature is ``Walls.of(geometry, exterior)``; a deck-only slab is read with the empty exterior.

    Reads ``inspect.signature`` (the specified pair), then the two decks' walls: a mirror slab ``(1, 0)`` at each end with itself as
    partner, a wrap slab ``(1, 0)`` with the opposite end as partner. Red at base: the import."""
    import inspect

    from orpheus.derivations.continuous.characteristic import Walls
    from orpheus.geometry import MirrorDeck, StructuredGeometry, TranslationDeck

    names = [p for p in inspect.signature(Walls.of).parameters]
    require(names == ["geometry", "exterior"], f"Walls.of takes {names}")
    mirror = StructuredGeometry.slab((0.0, 1.0, 2.0), (0, 1), left=MirrorDeck(), right=MirrorDeck())
    require(_read(Walls.of(mirror, {})) == [(0, 1.0, 0.0, 0), (2, 1.0, 0.0, 2)], "a mirror slab")
    wrap = StructuredGeometry.slab((0.0, 1.0, 2.0), (0, 1), left=TranslationDeck(), right=TranslationDeck())
    require(_read(Walls.of(wrap, {})) == [(0, 1.0, 0.0, 2), (2, 1.0, 0.0, 0)], "a wrap slab")


# ── Wall composes WallAmplitudes (decision 54) ───────────────────────────


def _amplitudes(specular: float, diffuse: float):
    from orpheus.specification.exterior import WallAmplitudes

    return WallAmplitudes(specular, diffuse)


@pytest.mark.rests_on(_LEAF)
def test_wall_composes_wall_amplitudes_and_exposes_read_only_coordinates() -> None:
    """``Wall(breakpoint, amplitudes, partner)``: the field set, the properties, equality by content, no setter.

    The bounds have ONE owner because the field holds the value that owns them. The properties return the held value's
    coordinates; assigning to either raises. Equality is by content (two walls of equal amplitudes are equal, a different amplitude
    differs). Mutation (5): ``Wall`` given back float fields reddens the first two assertions. Red at base: the import."""
    from orpheus.derivations.continuous.characteristic import Wall

    require([f.name for f in dataclasses.fields(Wall)] == ["breakpoint", "amplitudes", "partner"], "the fields of Wall")
    amplitudes = _amplitudes(0.37, 0.0)
    wall = Wall(3, amplitudes, 3)
    require(wall.amplitudes is amplitudes, "Wall does not hold the amplitudes it was given")
    require((wall.specular, wall.diffuse) == (0.37, 0.0), f"the properties read {(wall.specular, wall.diffuse)}")
    with pytest.raises(AttributeError):
        wall.specular = 0.5  # pyright: ignore[reportAttributeAccessIssue]
    with pytest.raises(AttributeError):
        wall.diffuse = 0.5  # pyright: ignore[reportAttributeAccessIssue]
    require(wall == Wall(3, _amplitudes(0.37, 0.0), 3), "two walls of equal content differ")
    require(wall != Wall(3, _amplitudes(0.38, 0.0), 3) and wall != Wall(3, _amplitudes(0.0, 0.37), 3), "a different amplitude is the same wall")


@pytest.mark.rests_on(_HERE + "test_wall_composes_wall_amplitudes_and_exposes_read_only_coordinates")
def test_the_old_four_float_wall_spelling_is_refused_as_a_type() -> None:
    """``Wall(3, 0.6, 1.0, 3)`` (the pre-S5 spelling, ``breakpoint, specular, diffuse, partner``) is a ``TypeError``.

    32 call sites were re-spelled by script (5 in ``walls.py``, 27 in 5 test files); a site the script missed fails loudly
    here instead of building a wall whose ``amplitudes`` is a float. The row also fixes that ``Wall`` takes a ``WallAmplitudes`` and
    not a bare pair. Red at base: the import."""
    from orpheus.derivations.continuous.characteristic import Wall

    with pytest.raises(TypeError):
        Wall(3, 0.6, 1.0, 3)  # pyright: ignore[reportCallIssue]
    with pytest.raises(TypeError):
        Wall(3, (0.6, 0.0), 3)  # pyright: ignore[reportArgumentType]


@pytest.mark.rests_on(_HERE + "test_wall_composes_wall_amplitudes_and_exposes_read_only_coordinates", "tests/gates/specification/test_exterior.py::test_wall_amplitudes_a_value_outside_the_bounds_is_a_parse_refusal")
@pytest.mark.parametrize("specular,diffuse", [(1.0, 1.0), (0.0, 1.0 + 1e-12), (0.6, 1.0 + 1e-9), (1.0 + 1e-12, 0.0), (-1e-12, 0.0)])
def test_an_out_of_bounds_wall_is_a_value_error_not_a_scope_refusal(specular: float, diffuse: float) -> None:
    """R-8b: the bounds refusal is ``ValueError`` (fragment "not a physical wall"), no longer ``NotImplementedError``.

    The rows are the three re-pinned sites (``test_characteristic_assembly.py`` ``returns_twice (1, 1)`` and
    ``diffuse_above_one (0, 1 + 1e-12)``; ``test_characteristic_walls.py`` ``(0.6, 1 + 1e-9)``) plus the two edges of the
    bound at 1 ulp scale. The TYPE is asserted exactly (``type(error) is ValueError``): a ``NotImplementedError`` would be a
    scope refusal that a later widening lifts, and an out-of-bounds wall is not liftable (no widening makes it physical).
    Mutation: ``Wall`` re-introduces its own ``_refuse``. Red at base: the import."""
    from orpheus.derivations.continuous.characteristic import Wall

    # The refusal is raised by the argument's own construction, so the row first requires ``Wall`` to take a ``WallAmplitudes``
    # (at the base it takes four floats and this line is the red; without it the row is green for the wrong object).
    require(Wall(3, _amplitudes(0.37, 0.0), 3).specular == 0.37, "Wall does not compose a WallAmplitudes")
    with pytest.raises(ValueError, match="not a physical wall") as caught:
        Wall(3, _amplitudes(specular, diffuse), 3)
    require(type(caught.value) is ValueError, f"the refusal is a {type(caught.value).__name__}, not a ValueError")


@pytest.mark.rests_on(_HERE + "test_wall_composes_wall_amplitudes_and_exposes_read_only_coordinates")
def test_a_combined_wall_is_still_the_scope_refusal_naming_650() -> None:
    """``Wall(.., WallAmplitudes(0.5, 0.5), ..)`` raises ``NotImplementedError`` naming #650, the one refusal that IS a scope boundary.

    The counterpart of the row above: the two refusals keep their two types, so a ``ValueError`` here would mean the combined wall
    had been demoted to a bad input and #650's widening could not retire it. Red at base: the import."""
    from orpheus.derivations.continuous.characteristic import Wall

    require(Wall(0, _amplitudes(0.5, 0.0), 0).specular == 0.5, "Wall does not compose a WallAmplitudes")  # red at base: the old signature
    with pytest.raises(NotImplementedError, match="#650"):
        Wall(0, _amplitudes(0.5, 0.5), 0)


# ── specular_albedos and require_vacuum on (geometry, exterior) ──────────


@pytest.mark.rests_on(_HERE + "test_walls_of_reads_each_body_as_the_wall_its_physics_names")
@pytest.mark.parametrize(("name", "chart", "breakpoints", "faces"), [c for c in cases.CASES if c[0] in cases.SPECULAR_ONLY], ids=cases.SPECULAR_ONLY)
def test_specular_albedos_reads_the_pair_geometry_and_exterior(name, chart, breakpoints, faces) -> None:
    """``specular_albedos(geometry, exterior, *, owner)`` is the by-hand albedos of the body, and what HEAD returned for it.

    A mirror deck reads ``1.0``; a tag reads its ``exterior[tag].specular`` when its ``diffuse`` is 0. Population: the
    specular-only bodies of the table (a vacuum, a mirror or a partial specular wall on each boundary point). Red at base: the
    import."""
    from orpheus.derivations.common.reference_body import specular_albedos

    geometry, exterior = spelling.geometry_and_exterior(chart, breakpoints, faces)
    got = specular_albedos(geometry, exterior, owner="carry")
    by_hand = tuple(cases.albedo(d) for d in faces)
    require(got == by_hand, f"{name}: read {got}, by hand {by_hand}")
    require(list(got) == _head()[name]["albedos"], f"{name}: read {got}, HEAD {_head()[name]['albedos']}")


@pytest.mark.rests_on(_HERE + "test_specular_albedos_reads_the_pair_geometry_and_exterior")
@pytest.mark.parametrize(
    ("name", "chart", "breakpoints", "faces"),
    [c for c in cases.CASES if c[0] not in cases.SPECULAR_ONLY],
    ids=[c[0] for c in cases.CASES if c[0] not in cases.SPECULAR_ONLY],
)
def test_specular_albedos_refuses_a_diffuse_wall_and_a_wrap_through_the_door(name, chart, breakpoints, faces) -> None:
    """A body with a diffuse wall or a wrap is refused through ``refuse_unserved`` (``NotImplementedError`` naming #536).

    The refusals are the ones HEAD made for the same declarations (a white law, a periodic law); the NEGATIVE leg of the row
    above, over the complement of its population. Red at base: the import."""
    from orpheus.derivations.common.reference_body import specular_albedos

    geometry, exterior = spelling.geometry_and_exterior(chart, breakpoints, faces)
    with pytest.raises(NotImplementedError, match="#536"):
        specular_albedos(geometry, exterior, owner="carry")


@pytest.mark.rests_on(_HERE + "test_specular_albedos_reads_the_pair_geometry_and_exterior")
def test_require_vacuum_reads_the_same_pair() -> None:
    """``require_vacuum(geometry, exterior, *, owner, missing)`` passes a vacuum body and refuses a mirror one, by the same albedos.

    Positive leg: ``slab-vacuum-vacuum`` returns ``None``. Negative legs: a mirror deck and a partial tag are refused with the
    caller's ``missing`` text in the message. Red at base: the import."""
    from orpheus.derivations.common.reference_body import require_vacuum

    by_name = {c[0]: c for c in cases.CASES}
    for name, refused in (("slab-vacuum-vacuum", False), ("slab-mirror-vacuum", True), ("slab-partial-vacuum", True), ("sphere-vacuum", False), ("sphere-mirror", True)):
        _n, chart, bp, faces = by_name[name]
        geometry, exterior = spelling.geometry_and_exterior(chart, bp, faces)
        if refused:
            with pytest.raises(NotImplementedError, match="the missing solver"):
                require_vacuum(geometry, exterior, owner="carry", missing="the missing solver")
        else:
            require_vacuum(geometry, exterior, owner="carry", missing="the missing solver")


# ── the three consumers gain ``exterior``; the Sood registry's finite cases carry the vacuum tag ─────


_CONSUMERS = (
    ("orpheus/derivations/continuous/singular_eigenfunction/spectrum.py", "Spectrum"),
    ("orpheus/derivations/continuous/fn_method/moment_space.py", "MomentSpace"),
    ("orpheus/derivations/continuous/galerkin_spectral/basis_space.py", "BasisSpace"),
)


def _class_field_names(source: str, name: str) -> list[str]:
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, ast.ClassDef) and node.name == name:
            return [n.target.id for n in node.body if isinstance(n, ast.AnnAssign) and isinstance(n.target, ast.Name)]
    raise AssertionError(f"no class {name}")


@pytest.mark.parametrize(("path", "cls"), _CONSUMERS, ids=[c[1] for c in _CONSUMERS])
def test_each_reference_consumer_has_an_exterior_field_beside_the_geometry(path: str, cls: str) -> None:
    """``Spectrum``, ``MomentSpace`` and ``BasisSpace`` declare ``exterior`` immediately after ``geometry`` (an AST read, import-free).

    Positive control: the reader finds ``geometry`` in each class (so a class the reader mis-parses cannot pass by being empty),
    and finds nothing in an absent class. The field is a plain-data field of the reference's own value, never a production
    one. Red at base: no ``exterior`` field."""
    fields = _class_field_names((_ROOT / path).read_text(), cls)
    require("geometry" in fields, f"control: {cls} fields read as {fields}")
    with pytest.raises(AssertionError, match="no class"):
        _class_field_names("class A: pass", "B")
    require("exterior" in fields, f"{cls} has no field exterior: {fields}")
    require(fields.index("exterior") == fields.index("geometry") + 1, f"exterior does not follow geometry in {cls}: {fields}")


@pytest.mark.rests_on(_HERE + "test_each_reference_consumer_has_an_exterior_field_beside_the_geometry")
def test_the_finite_sood_cases_carry_the_one_vacuum_tag() -> None:
    """Every finite case of the Sood registry builds a geometry whose faces all carry the one tag ``0`` (rule R1: equal declarations share a tag).

    The registry states its vacuum boundary as an integer tag, and the case's consumers state the binding
    ``{0: WallAmplitudes(0, 0)}``; the tag NUMBER is what a consumer's ``exterior`` keys on, so it is the contract this row
    pins. Population: every ``SOOD2003_CASES`` member whose ``geometry_kind`` is not ``infinite`` (counted in the message: 25 of the 47 cases at the base) and
    that carries a critical dimension. Red at base: the import (the registry builds a vacuum string-tag geometry)."""
    from orpheus.derivations.continuous.sood_registry import SOOD2003_CASES

    finite = [c for c in SOOD2003_CASES.values() if c.geometry_kind != "infinite" and c.truth.critical_dimension_mfp is not None]
    require(len(finite) > 5, f"control: only {len(finite)} finite cases found")
    bad = [c.case_id for c in finite if c.to_geometry().tags != frozenset({0}) or any(isinstance(r, int) is False for r in c.to_geometry().face_roles)]
    require(not bad, f"{len(bad)} of {len(finite)} finite cases do not carry exactly the tag 0 on every face: {bad[:5]}")
