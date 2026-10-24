r"""The specification states its exterior as plain data: ``GeometrySpecification.exterior`` (wave 2, S5; plan 2.9 and gate 1 of 6.5; law L-J).

Specified by the test-architect (2026-10-23, ``.claude/plans/wave2_boundary_carve.md`` 6.5, decisions D2, D3 and 54). A geometry's
faces carry a role: a deck (a symmetry of the domain, ``MirrorDeck`` or ``TranslationDeck``, which states its own response and
needs no entry) or an integer TAG of the exterior. The reference specification binds each tag to a ``WallAmplitudes(specular,
diffuse)`` in a trailing field ``exterior: Mapping[int, WallAmplitudes]``, defaulting to the empty mapping, so a deck-only
geometry (and every construction that carries no tag) keeps its positional three-argument form. The claim is that the
specification is plain data and COMPLETE: the exact cover of the tags, no vocabulary for what the reference cannot answer.

The laws gated, each a THEOREM about construction (the specification, not a solve):

* **L-J1 (exact cover).** ``exterior.keys() == geometry.tags``, as sets of integers: a tag with no amplitudes is refused, an entry
  with no tag is refused, each by NUMBER in the message; a tag shared by two faces needs ONE entry (the cover is by tag, never by
  face); a deck-only geometry takes the empty exterior and refuses any entry. ``test_a_tag_without_amplitudes_is_refused``,
  ``test_an_entry_without_a_tag_is_refused`` and ``test_the_cover_is_by_tag_not_by_face`` are the three legs, with the accepting
  leg (``test_an_exact_cover_constructs_and_reads_back``) as the positive control the refusals are read against.
* **L-J2 (the value is the closed one).** An entry is a ``WallAmplitudes``: a ``WallCondition`` (production), a bare pair, a float
  or ``None`` is a ``TypeError``, so the production vocabulary (a prescribed inflow, a combined chain, a zero-flux
  condition) cannot be written here at all, and the specification imports nothing from transport. The unwritable wall that D3
  names is unwritable by TYPE: ``WallAmplitudes(0.5, 0.5)`` raises ``NotImplementedError`` naming #650 before any
  specification sees it (``test_a_combined_wall_cannot_reach_the_specification``), and there is no field for an inflow
  (``test_the_specification_has_no_field_for_what_the_reference_cannot_answer``).
* **L-J3 (content).** ``exterior`` is part of the specification's identity: two specifications differing only in it have
  different digests, one differing only in insertion order of the same entries has the same digest and is equal, and the
  amplitudes' two coordinates are not interchangeable (``(0.3, 0)`` and ``(0, 0.3)`` are two keys). A copy of the passed
  mapping is held (a caller's later mutation does not change the specification) and the held mapping refuses assignment.

Declared red at the base tip (``de354496``), the input that reddens each row (``instrument-doctrine`` X1): every row builds a
specification with an ``exterior`` keyword, which ``GeometrySpecification`` does not accept at the base (``TypeError``), so the
FIRST red is that kwarg and the base cannot tell the rows apart; their named reasons are proven by the mutations the reviewers
apply to the carved tree, per row: (1) the exact-cover check removed from ``__post_init__`` reddens the two refusal rows (a
missing tag then passes, and the build fails later in ``Walls.of``), (2) the cover compared by FACE count instead of tag set
reddens ``test_the_cover_is_by_tag_not_by_face``, (3) ``exterior`` left out of the specification's content parts reddens the digest
row (the planted-extra-field control of the plan: a field the digest ignores is exactly this), (4) the entry type unchecked
reddens the ``TypeError`` row, (5) the mapping stored by reference reddens the copy row. The ``WallAmplitudes`` rows
themselves (bounds, combined wall, one owner) are ``tests/gates/specification/test_exterior.py`` (L5, landed), on which every
row here rests.

The specification digest pin (the plan's R-2 re-baseline, gate 7 of 6.5) is NOT here: its bytes are the build's (the schema is
final only once the field exists), recorded in ``test_content_identity_specification.py::test_s8_5_the_digest_bytes_are_recorded``
at the commit that lands the field.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Mapping

import pytest

from orpheus.data.cells import CellCoefficient, Channel
from orpheus.data.materials import Materials
from orpheus.geometry import MirrorDeck, StructuredGeometry
from orpheus.numerics.question import Eigen
from orpheus.specification import GeometrySpecification
from orpheus.specification.exterior import WallAmplitudes
from orpheus.values.content import content_digest
from tests.gates._content_identity_helpers import require
from tests.gates.specification._fixtures import fuel, moderator, slab2

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/specification/test_exterior_on_the_specification.py::"
_EXTERIOR = "tests/gates/specification/test_exterior.py::"
_LEAF = _EXTERIOR + "test_wall_amplitudes_a_stated_wall_constructs"
_VACUUM = WallAmplitudes(0.0, 0.0)


def _question() -> Eigen:
    return Eigen(CellCoefficient.every(Channel.FISSION_EMISSION))


def _materials() -> Materials:
    return Materials({0: fuel(), 1: moderator()})


def _spec(geometry: StructuredGeometry, **exterior_kw) -> GeometrySpecification:
    """The specification of ``geometry``; the exterior is passed by keyword only when the row states one."""
    return GeometrySpecification(_materials(), geometry, _question(), **exterior_kw)


def _deck_only() -> StructuredGeometry:
    return StructuredGeometry.slab((0.0, 1.0, 3.0), (0, 1), left=MirrorDeck(), right=MirrorDeck())


def _two_tag_hollow_sphere() -> StructuredGeometry:
    return StructuredGeometry.sphere((0.4, 1.0, 3.0), (0, 1), inner=1, outer=2)


def _shared_tag_slab() -> StructuredGeometry:
    return StructuredGeometry.slab((0.0, 1.0, 3.0), (0, 1), left=0, right=0)


@pytest.mark.rests_on(_LEAF)
def test_the_exterior_is_the_trailing_field_and_defaults_to_the_empty_mapping() -> None:
    """``exterior`` is the LAST field (the three-argument positional form keeps working) and a deck-only geometry takes none.

    The default is empty, not ``None`` (the plan: "a trailing field defaulting to an empty ``FrozenMapping``"), so every
    reader reads a mapping. Red at base: no such field (``fields`` ends at ``question``)."""
    names = [f.name for f in dataclasses.fields(GeometrySpecification)]
    require(names[:3] == ["materials", "geometry", "question"], f"the leading fields moved: {names}")
    require(names[-1] == "exterior", f"exterior is not the trailing field: {names}")
    spec = _spec(_deck_only())
    require(isinstance(spec.exterior, Mapping), f"the default exterior is a {type(spec.exterior).__name__}, not a Mapping")
    require(len(spec.exterior) == 0, f"the default exterior is not empty: {dict(spec.exterior)}")


@pytest.mark.rests_on(_LEAF, _HERE + "test_the_exterior_is_the_trailing_field_and_defaults_to_the_empty_mapping")
@pytest.mark.parametrize(
    ("make_geometry", "exterior"),
    [
        (slab2, {0: _VACUUM}),
        (slab2, {0: WallAmplitudes(0.37, 0.0)}),
        (_two_tag_hollow_sphere, {1: WallAmplitudes(0.0, 0.6), 2: WallAmplitudes(1.0, 0.0)}),
        (_deck_only, {}),
    ],
    ids=["slab-vacuum", "slab-partial", "hollow-sphere-two-tags", "deck-only"],
)
def test_an_exact_cover_constructs_and_reads_back(make_geometry, exterior) -> None:
    """The POSITIVE leg: ``exterior.keys() == geometry.tags`` constructs, and the stated amplitudes read back unchanged.

    The control the refusals below are read against: without an accepting leg a refusal row passes a specification that
    refuses everything. Red at base: the ``exterior`` keyword."""
    geometry = make_geometry()
    spec = _spec(geometry, exterior=exterior)
    require(set(spec.exterior) == set(geometry.tags), f"the tags {set(geometry.tags)} read back as {set(spec.exterior)}")
    for tag, wall in exterior.items():
        require(spec.exterior[tag].specular == wall.specular and spec.exterior[tag].diffuse == wall.diffuse, f"tag {tag} reads back {spec.exterior[tag]!r}")


@pytest.mark.rests_on(_HERE + "test_an_exact_cover_constructs_and_reads_back")
@pytest.mark.parametrize(
    ("make_geometry", "exterior", "missing"),
    [
        (slab2, {}, 0),
        (_two_tag_hollow_sphere, {1: _VACUUM}, 2),
        (_two_tag_hollow_sphere, {2: _VACUUM}, 1),
    ],
    ids=["slab-no-entry", "hollow-sphere-outer-only-missing", "hollow-sphere-inner-missing"],
)
def test_a_tag_without_amplitudes_is_refused(make_geometry, exterior, missing: int) -> None:
    """A tag a face carries and the exterior does not state is a ``ValueError`` naming that tag's number.

    The cover is checked at construction (D3: today's late refusal moves to the specification's door). The message names the
    tag NUMBER, so the diagnosis is a lookup. Mutation: drop the cover check, and the first row builds a specification that
    only fails inside ``Walls.of``. Red at base: the ``exterior`` keyword."""
    geometry = make_geometry()
    with pytest.raises(ValueError, match=rf"\b{missing}\b"):
        _spec(geometry, exterior=exterior)


@pytest.mark.rests_on(_HERE + "test_an_exact_cover_constructs_and_reads_back")
@pytest.mark.parametrize(
    ("make_geometry", "exterior", "extra"),
    [
        (slab2, {0: _VACUUM, 5: _VACUUM}, 5),
        (_deck_only, {7: _VACUUM}, 7),
        (_two_tag_hollow_sphere, {1: _VACUUM, 2: _VACUUM, 3: _VACUUM}, 3),
    ],
    ids=["slab-extra", "deck-only-takes-none", "hollow-sphere-extra"],
)
def test_an_entry_without_a_tag_is_refused(make_geometry, exterior, extra: int) -> None:
    """An entry for a tag no face carries is a ``ValueError`` naming it; a deck-only geometry takes the empty exterior only.

    The mirror of the previous refusal (an exact cover has two sides; the production side of the same law is
    ``UnusedConditionError``). Mutation: check only ``geometry.tags <= exterior.keys()`` (the missing side), and every row
    here builds. Red at base: the ``exterior`` keyword."""
    geometry = make_geometry()
    with pytest.raises(ValueError, match=rf"\b{extra}\b"):
        _spec(geometry, exterior=exterior)


@pytest.mark.rests_on(_HERE + "test_an_exact_cover_constructs_and_reads_back")
def test_the_cover_is_by_tag_not_by_face() -> None:
    """A slab whose two faces carry the SAME tag takes ONE entry; two entries (one a tag no face carries) are refused.

    The discriminator between a cover by tag (the specified law) and a cover by face count (a plausible wrong one, ``len(exterior)
    == number of tagged faces``): the shared-tag slab has two tagged faces and one tag. Red at base: the ``exterior`` keyword."""
    geometry = _shared_tag_slab()
    require(geometry.tags == frozenset({0}), f"control: the fixture carries tags {set(geometry.tags)}, not {{0}}")
    _spec(geometry, exterior={0: WallAmplitudes(0.5, 0.0)})
    with pytest.raises(ValueError):
        _spec(geometry, exterior={0: _VACUUM, 1: _VACUUM})


@pytest.mark.rests_on(_HERE + "test_an_exact_cover_constructs_and_reads_back")
@pytest.mark.parametrize(
    "value",
    [(0.3, 0.0), 0.3, None, "vacuum", [0.3, 0.0]],
    ids=["a-pair", "a-float", "none", "a-string", "a-list"],
)
def test_an_entry_that_is_not_wall_amplitudes_is_a_type_error(value) -> None:
    """The value of an entry is a ``WallAmplitudes``: nothing else is admitted, which is what keeps the vocabulary closed.

    Mutation: leave the entry unchecked and a bare pair builds a specification whose digest hashes a tuple. Red at base: the
    ``exterior`` keyword, through the inline accepting control: at the base an unexpected-keyword ``TypeError`` would satisfy
    ``pytest.raises(TypeError)`` for the wrong reason (``vv-principles`` mode 8(4)), so the row first requires the SAME
    keyword to be accepted with a legal value."""
    _spec(slab2(), exterior={0: _VACUUM})
    with pytest.raises(TypeError):
        _spec(slab2(), exterior={0: value})


@pytest.mark.rests_on(_HERE + "test_an_entry_that_is_not_wall_amplitudes_is_a_type_error")
def test_a_production_condition_cannot_be_stated_in_the_specification() -> None:
    """A ``WallCondition`` (the production vocabulary: a chain sum, an inflow, a displacement) is a ``TypeError`` as an entry.

    The wall a reference cannot answer (a prescribed inflow, a combined chain) is unspellable here by TYPE, with no ``isinstance``
    on the production class (the specification imports nothing from transport: the value is simply not a ``WallAmplitudes``).
    The production type is built through its public facades (``vacuum``, ``uniform``). Red at base: the ``exterior`` keyword."""
    from orpheus.transport.boundary_conditions import WallCondition

    _spec(slab2(), exterior={0: _VACUUM})  # the accepting control: an unexpected-keyword TypeError must not pass the rows below
    for condition in (WallCondition.vacuum(), WallCondition.specular(0.5), WallCondition.uniform(1.0)):
        with pytest.raises(TypeError):
            _spec(slab2(), exterior={0: condition})


@pytest.mark.rests_on(_LEAF, _EXTERIOR + "test_wall_amplitudes_a_combined_wall_is_a_scope_boundary_naming_650")
def test_a_combined_wall_cannot_reach_the_specification() -> None:
    """``WallAmplitudes(0.5, 0.5)`` raises ``NotImplementedError`` naming #650, and no specification is ever built from it.

    The value is refused where it is built, not where it is placed (the specification's cover cannot be asked to refuse a wall
    it never sees); an exterior assembled in one expression propagates that refusal. This row is GREEN at the base (the leaf
    landed at L5 and the refusal is raised while the argument is evaluated, before any ``exterior`` keyword is reached): it is
    the live control that the closed vocabulary's refusal precedes the specification's, and reddens only if the leaf's
    combined-wall refusal is weakened (mutation: ``WallAmplitudes`` accepts ``(0.5, 0.5)``)."""
    with pytest.raises(NotImplementedError, match="#650"):
        _spec(slab2(), exterior={0: WallAmplitudes(0.5, 0.5)})


@pytest.mark.rests_on(_HERE + "test_the_exterior_is_the_trailing_field_and_defaults_to_the_empty_mapping")
def test_the_specification_has_no_field_for_what_the_reference_cannot_answer() -> None:
    """The fields are exactly ``(materials, geometry, question, exterior)``: no inflow, no zero flux, no boundary law field.

    D3 by construction: a specification posing a prescribed inflow, a combined wall or a zero flux cannot be written because
    there is no field. The exact field set is asserted (a leftover ``boundary_conditions`` or a new ``inflow`` field reddens)
    and ``exterior``'s annotation is the closed value (``WallAmplitudes``), read from the dataclass's field type string. Red at
    base: the field set lacks ``exterior``."""
    fields = {f.name: f for f in dataclasses.fields(GeometrySpecification)}
    require(list(fields) == ["materials", "geometry", "question", "exterior"], f"the field set is {list(fields)}")
    annotation = str(fields["exterior"].type)
    require("WallAmplitudes" in annotation, f"exterior is annotated {annotation!r}, not a mapping to WallAmplitudes")


@pytest.mark.rests_on(_HERE + "test_an_exact_cover_constructs_and_reads_back")
def test_two_specifications_differing_only_in_the_exterior_have_different_digests() -> None:
    """``exterior`` is content: a wall's coordinates move the digest, the two coordinates are not interchangeable, insertion order does not.

    Legs (each a digest comparison on slab2 / the two-tag hollow sphere, everything else held equal): vacuum against 0.3
    specular (differ); 0.3 specular against 0.3 diffuse (differ: the coordinates are ordered); the same stated exterior built
    twice (equal digest, equal specifications: the positive control, a digest that is a per-instance counter fails it); the
    two-tag exterior in either insertion order (equal digest: a mapping's order is not content). Mutation (3) of the module
    docstring: ``exterior`` left out of the content parts, which reddens the first two legs and nothing else. Red at base: the
    ``exterior`` keyword."""
    vacuum = _spec(slab2(), exterior={0: _VACUUM})
    specular = _spec(slab2(), exterior={0: WallAmplitudes(0.3, 0.0)})
    diffuse = _spec(slab2(), exterior={0: WallAmplitudes(0.0, 0.3)})
    require(content_digest(vacuum) != content_digest(specular), "a vacuum wall and a 0.3 specular wall share a digest")
    require(content_digest(specular) != content_digest(diffuse), "a 0.3 specular wall and a 0.3 diffuse wall share a digest")
    again = _spec(slab2(), exterior={0: WallAmplitudes(0.3, 0.0)})
    require(content_digest(specular) == content_digest(again) and specular == again, "one stated exterior built twice is two keys")
    forward = _spec(_two_tag_hollow_sphere(), exterior={1: _VACUUM, 2: WallAmplitudes(1.0, 0.0)})
    backward = _spec(_two_tag_hollow_sphere(), exterior={2: WallAmplitudes(1.0, 0.0), 1: _VACUUM})
    require(content_digest(forward) == content_digest(backward) and forward == backward, "the insertion order of the exterior is content")


@pytest.mark.rests_on(_HERE + "test_an_exact_cover_constructs_and_reads_back")
def test_the_specification_holds_its_own_copy_of_a_frozen_exterior() -> None:
    """A caller mutating the mapping it passed does not change the specification, and the held mapping refuses assignment.

    A content-keyed object that aliased a mutable dict would change its own digest after construction. Red at base: the
    ``exterior`` keyword. Mutation (5): store ``exterior`` as given."""
    passed = {0: WallAmplitudes(0.3, 0.0)}
    spec = _spec(slab2(), exterior=passed)
    before = content_digest(spec)
    passed[0] = WallAmplitudes(0.9, 0.0)
    passed[4] = _VACUUM
    require(spec.exterior[0].specular == 0.3 and 4 not in spec.exterior, "the specification shares the caller's mapping")
    require(content_digest(spec) == before, "the digest moved with the caller's mapping")
    with pytest.raises(TypeError):
        spec.exterior[0] = _VACUUM  # pyright: ignore[reportIndexIssue]
