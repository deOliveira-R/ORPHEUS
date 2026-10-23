r"""The algebra of the shared condition: the cone as a VALUE predicate (wave 2, slice S2a, gates 4 to 8).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 2.4 ("the algebra, and ``displacement.py``", decision 50),
section 3 (laws L-B, L-F, L-J), 6.2 S2a. The toy that measured these rules is
``scratch/layer_dag/wave2_spec/_d3_cone_arith.py`` (G1 to G9); the expected values below are transcribed from it,
not from the prose.

The claims (``vv-principles``), all THEOREM on the type's own coordinates ``(s, d, uniform, sampled)``:

* gate 4, ``test_the_algebra_and_the_cone``: ``+``, ``-``, unary ``-`` and scalar ``*`` act on the pair
  ``(network, inflow)`` and exit through ONE function: a ``WallCondition`` iff the result lies in the cone ``K``
  (every network weight ``>= 0``, ``math.fsum(w * unit_return) <= 1.0``, every inflow weight ``>= 0``), else a
  ``WallDisplacement``. The sign of a scalar and the presence of ``-`` decide nothing. Each identity is stated with
  its quantifier: exact on dyadic coordinates, ONE ulp for the ``0.1, 0.2, 0.3`` triple (a wall written two ways can
  differ in digest: a cache miss, never a wrong answer).
* gate 5, ``test_the_source_is_summed_with_the_response``: the affine part ``q`` travels with the linear part ``R``
  through every combination (the defect of DRAFT 1: an algebra on responses that dropped the source).
* gate 6, ``test_displacements_are_values_not_conditions``: a displacement is the affine hull of ``K x (positive
  cone)``, closed under the same operators, never bindable, never equal to the condition of the same coordinates.
* gate 7, ``test_a_condition_is_not_a_deck_slot`` and ``test_a_deck_is_not_an_operand`` (law L-F): no value of the
  condition types puts a response in the deck slot, no deck is an operand of the algebra, and the refusal comes from
  Python's operator protocol (no production frame in the traceback).
* gate 8, ``test_inflow_frames_and_terms``: the frame is the type (``FaceFrameInflow`` is face-local,
  ``LabFrameInflow`` is not), the content is open, and ``invariant_under`` is conservative and deterministic.

Activation (which term each fixture activates and which it nulls, ``vv-principles`` mode 7): the cone grid activates
both chains and the source; the ``(0.6, 0.6)`` fixtures activate the cone sum (a single-chain fixture nulls it); the
``0.1, 0.2, 0.3`` triple activates the rounding that the dyadic fixture nulls; the framed profile activates the
sampled terms that a uniform-only fixture nulls.

First reds, at the base tip ``4ec95bc6``: every row, ``ModuleNotFoundError: No module named
'orpheus.transport.boundary_conditions'``. The input that reddens gate 4 once the module exists, ``[M]`` at the base
tip, read from the retired vocabulary: ``LawSum(LawScaled(0.6, spec), LawScaled(0.6, white))`` constructs (a
super-Markov wall) and ``(a + b) + c == a + (b + c)`` is ``False`` for three laws (``_composition.py:206-222``).

Mutations (each on a COPY of the module; the reviewers' battery):

* gate 4: drop the cone parse in ``WallCondition.__post_init__``; a ``__add__`` that sums the networks and drops the
  sources (the V1 defect: gate 5 reddens too); a copy that decides the sort by the scalar's sign (DRAFT 2's rule),
  which raises on ``2.0 * specular(0.6)`` where the predicate returns a displacement; a copy that keeps ``(0.0, f)``
  in ``sampled`` (``0.0 * framed(f) == vacuum()`` reddens); a copy that compares the cone sum with ``<`` (the
  ``0.9 + 0.1`` row).
* gate 6: let ``BoundaryConditions`` accept a ``WallDisplacement``; make ``of`` return a ``WallCondition`` when the
  combination lies in ``K``.
* gate 7: give ``MirrorDeck`` an ``__rmul__`` (``0.5 * MirrorDeck()`` stops raising); let ``WallCondition.__add__``
  accept any object with a ``network`` attribute.
* gate 8: make ``LabFrameInflow`` report ``face_local`` True; make ``invariant_under`` ignore ``is_identity``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from orpheus.geometry.colouring import colouring_symmetry
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.deck import MirrorDeck, TranslationDeck, parse_face_role
from orpheus.values.content import ContentlessError, content_digest
from tests.gates._content_identity_helpers import require
from tests.gates.transport._boundary_conditions_helpers import (
    OpaqueSampler,
    ProfileSampler,
    bc,
    cone_grid,
    function_space,
    local,
    no_orpheus_frame,
    node,
    wall,
    weights,
)

pytestmark = pytest.mark.foundation

_HERE = "test_wall_condition_algebra.py"
_PARSE = node("test_wall_condition.py", "test_wall_condition_parse")
_ALGEBRA = node(_HERE, "test_the_algebra_and_the_cone")
_DECK = "tests/gates/geometry/test_deck.py"
_COLOURING = "tests/gates/geometry/test_colouring_symmetry.py::test_colouring_symmetry_reads_runs_and_roles"


def _spec(a: float = 1.0):
    return bc.WallCondition.specular(a)


def _white(a: float = 1.0):
    return bc.WallCondition.white(a)


def _is_condition(x) -> bool:
    return type(x) is bc.WallCondition


def _is_displacement(x) -> bool:
    return type(x) is bc.WallDisplacement


def _one_ulp(a: float, b: float) -> bool:
    return abs(a - b) <= math.ulp(max(abs(a), abs(b)))


# ── gate 4: the algebra and the cone (laws L-B, L-J) ─────────────────────────


@pytest.mark.rests_on(_PARSE)
def test_the_algebra_and_the_cone() -> None:
    """Every identity of the section, with its stated quantifier; the value predicate with its measured cases
    (G1 to G9 of the toy); direct construction outside ``K`` raises; a sampled term is dropped at weight 0.

    Exact (on every admitted value): ``a + b == b + a`` over the cone grid; ``vacuum() + x == x``; ``1.0 * x == x``;
    ``0.0 * x == vacuum()``; ``0.5 * specular(1.0) == specular(0.5)``; ``specular(0.25) + specular(0.5) ==
    specular(0.75)``; the dyadic association ``(a + b) + c == a + (b + c)`` on ``0.25, 0.5, 0.125``. One ulp per
    coordinate: the ``0.1, 0.2, 0.3`` triple, where the digests differ (``[M]`` ``(0.1 + 0.2) + 0.3 =
    0.6000000000000001``, ``0.1 + (0.2 + 0.3) = 0.6``: the limit is stated, not hidden).

    Value predicate (decision 50): ``specular(0.5) + (white(0.3) - white(0.3)) == specular(0.5)`` and is a
    condition (G1); ``2.0 * specular(0.6)`` is a displacement (G2, which discriminates the value rule from the sign
    rule); ``specular(0.5) - specular(0.3)`` is a condition (G3); ``white(0.3) - specular(0.5)`` and ``0.6 specular
    + 0.6 white`` are displacements (G4, G6); ``0.0 * framed(f) == vacuum()`` (G5); ``-a`` and ``-0.5 * a`` are
    displacements for EVERY ``a`` with a nonzero coordinate, and ``-vacuum()`` is ``vacuum()`` (a condition); a
    displacement plus its inverse step re-enters ``K`` as a condition (G7); ``specular(1.2)``, ``white(-0.1)`` and a
    direct two-chain ``(0.6, 0.6)`` raise (G8); ``specular(0.9) + white(0.1)`` is a condition by the float cone sum
    (G9). A test-local condition with ``return_fraction < 0`` defines no operator, so it cannot enter a sum.

    Reddens on: the mutations of the module docstring.
    """
    grid = cone_grid()
    walls = [wall(s, d) for s, d in grid]
    require(len(walls) == 48, f"the grid has {len(walls)} walls, not 48")

    # exact identities, over the whole grid
    vac = bc.WallCondition.vacuum()
    for x in walls:
        require(vac + x == x and x + vac == x, f"vacuum() + x != x at {x.cone_point}")
        require(1.0 * x == x and x * 1.0 == x, f"1.0 * x != x at {x.cone_point}")
        require(0.0 * x == vac and _is_condition(0.0 * x), f"0.0 * x != vacuum() at {x.cone_point}")
    for x, y in zip(walls, reversed(walls)):
        require((x + y) == (y + x), f"a + b != b + a at {x.cone_point}, {y.cone_point}")
        require(type(x + y) is type(y + x), "the sort of a + b depends on the order")
    require(0.5 * _spec(1.0) == _spec(0.5), "0.5 * specular(1.0) != specular(0.5)")
    require(_spec(0.25) + _spec(0.5) == _spec(0.75), "the dyadic merge is not exact")
    a, b, c = _spec(0.25), _white(0.5), _spec(0.125)
    require((a + b) + c == a + (b + c), "the dyadic association is not exact")
    require(((a + b) + c).cone_point == (0.375, 0.5), f"the dyadic sum reads {((a + b) + c).cone_point}")
    require(math.fsum([0.25, 0.5, 0.125]) == 0.875, "control: the dyadic triple sums exactly")
    require((a + b).return_fraction == a.return_fraction + b.return_fraction, "return_fraction is not additive on dyadic coordinates")

    # the one-ulp leg: the associativity of general doubles is NOT exact, and the limit is stated
    x1, x2, x3 = _spec(0.1), _spec(0.2), _spec(0.3)
    left, right = (x1 + x2) + x3, x1 + (x2 + x3)
    require(_is_condition(left) and _is_condition(right), "the 0.1, 0.2, 0.3 sums lie in K")
    require(0.1 + 0.2 + 0.3 != 0.1 + (0.2 + 0.3), "control: the two associations of the triple differ in the double")
    for lc, rc in zip(left.cone_point, right.cone_point):
        require(_one_ulp(lc, rc), f"the associations differ by more than one ulp: {lc!r} vs {rc!r}")
    require(left != right and content_digest(left) != content_digest(right),
            "the two spellings of one physical wall coincide in digest: associativity became exact (rational coordinates?)")

    # the value predicate
    g1 = _spec(0.5) + (_white(0.3) - _white(0.3))
    require(g1 == _spec(0.5) and _is_condition(g1), "G1: specular(0.5) + (white(0.3) - white(0.3)) is not specular(0.5)")
    g2 = 2.0 * _spec(0.6)
    require(_is_displacement(g2), f"G2: 2.0 * specular(0.6) is a {type(g2).__name__}")
    g3 = _spec(0.5) - _spec(0.3)
    require(_is_condition(g3), f"G3: specular(0.5) - specular(0.3) is a {type(g3).__name__}")
    require(_is_displacement(_white(0.3) - _spec(0.5)), "G4: white(0.3) - specular(0.5) is not a displacement")
    require(_is_displacement(0.6 * _spec(1.0) + 0.6 * _white(1.0)), "G6: 0.6 specular + 0.6 white is not a displacement")
    require(_is_condition(_spec(0.9) + _white(0.1)), "G9: specular(0.9) + white(0.1) is not a condition (the float cone sum)")
    require(_is_displacement(_spec(0.9) + _white(0.1) + bc.WallCondition.from_amplitudes(0.0, 2.0**-52)),
            "one double step past 1 is a displacement")
    step_back = -0.6 * _white(1.0)
    g7 = (0.6 * _spec(1.0) + 0.6 * _white(1.0)) + step_back
    require(_is_condition(g7) and g7 == _spec(0.6), f"G7: a displacement plus its inverse step is {type(g7).__name__} {g7!r}")

    negated = [x for x in walls if x.cone_point != (0.0, 0.0)] + [bc.WallCondition.inflow(2.0), bc.WallCondition.framed(bc.FaceFrameInflow(ProfileSampler()))]
    require(len(negated) == 49, f"the negation population has {len(negated)} members, not 49 (47 nonzero grid points + 2 source-only walls)")
    for x in negated:
        require(_is_displacement(-x), f"-a is a {type(-x).__name__} at {x.cone_point}")
        require(_is_displacement(-0.5 * x), f"-0.5 * a is a {type(-0.5 * x).__name__} at {x.cone_point}")
    require(-vac == vac and _is_condition(-vac), "-vacuum() is vacuum(), a condition")

    # direct construction outside K raises, and so do the facades
    empty = weights()
    with pytest.raises(bc.SubmarkovViolationError):
        bc.WallCondition(weights(specular=0.6, diffuse=0.6))
    with pytest.raises(bc.BoundaryResponseNotPositiveError):
        bc.WallCondition(weights(specular=-0.1))
    with pytest.raises(bc.BoundaryError):
        bc.WallCondition(empty, bc.InflowTerms(-1.0, ()))
    with pytest.raises(bc.BoundaryError):
        bc.WallCondition(empty, bc.InflowTerms(0.0, ((-1.0, bc.FaceFrameInflow(ProfileSampler())),)))
    require(_is_condition(bc.WallCondition(empty)), "the inflow defaults to empty")

    # the zero-weight drop reaches the sampled terms
    framed = bc.WallCondition.framed(bc.FaceFrameInflow(OpaqueSampler()))
    dropped = 0.0 * framed
    require(dropped == vac and dropped.source.sampled == () and _is_condition(dropped), "G5: 0.0 * framed(f) keeps a zero-weight term")

    # a method-only condition defines no operator
    probe = local().NegativeProbe()
    for thunk in (lambda: probe + vac, lambda: vac + probe, lambda: probe - vac, lambda: -probe, lambda: 2.0 * probe, lambda: probe * 2.0):
        with pytest.raises(TypeError):
            thunk()


# ── gate 5: the source travels with the response ─────────────────────────────


@pytest.mark.rests_on(_ALGEBRA)
def test_the_source_is_summed_with_the_response() -> None:
    """Decision 32 (elegance V1, structure frame 1): ``0.5 * inflow(2.0) + specular(0.5)`` has
    ``source.uniform == 1.0`` and ``cone_point == (0.5, 0.0)``; ``inflow(1.0) + inflow(2.0) == inflow(3.0)``;
    ``0.5 * framed(FaceFrameInflow(f))`` neither raises nor loses the weight (the open-sampler case that raised in
    DRAFT 1's digest sort); a wall holding a sampled term has IDENTITY equality (it is contentless) and its digest
    raises ``ContentlessError``; ``a + b == b + a`` is asserted only for uniform-only sources.

    Reddens on: a ``__add__`` that sums the networks and drops the sources (the first assertion); a scalar ``*``
    that scales the network and not the source; a copy sorting sampled terms by digest (the open-sampler row).
    """
    mixed = 0.5 * bc.WallCondition.inflow(2.0) + _spec(0.5)
    require(_is_condition(mixed), f"the mixed wall is a {type(mixed).__name__}")
    require(mixed.source.uniform == 1.0 and mixed.cone_point == (0.5, 0.0), f"the mixed wall reads {mixed.cone_point}, {mixed.source}")
    require(bc.WallCondition.inflow(1.0) + bc.WallCondition.inflow(2.0) == bc.WallCondition.inflow(3.0), "inflow(1) + inflow(2) != inflow(3)")
    require((_spec(0.5) + bc.WallCondition.inflow(1.0)) == (bc.WallCondition.inflow(1.0) + _spec(0.5)), "a + b != b + a with uniform sources")
    require(2.0 * mixed == bc.WallCondition(weights(specular=1.0), bc.InflowTerms(2.0, ())), "scaling a mixed wall lost its source")
    difference = _spec(0.5) - bc.WallCondition.inflow(1.0)
    require(_is_displacement(difference) and difference.inflow.uniform == -1.0, "a - b does not subtract the sources")

    sampler = OpaqueSampler()
    profile = bc.FaceFrameInflow(sampler)
    half = 0.5 * bc.WallCondition.framed(profile)
    require(_is_condition(half), f"0.5 * framed(f) is a {type(half).__name__}")
    require(len(half.source.sampled) == 1, f"the sampled terms are {half.source.sampled}")
    weight, held = half.source.sampled[0]
    require(weight == 0.5 and held.sampler is sampler, f"the weight or the profile was lost: {weight}, {held}")
    require(bc.WallCondition.framed(profile, 2.0).source.sampled[0][0] == 2.0, "framed(p, w) does not store w")
    require(half.cone_point == (0.0, 0.0), "a framed wall has the empty network")
    require(half == half, "a wall equals itself")
    other = 0.5 * bc.WallCondition.framed(bc.FaceFrameInflow(OpaqueSampler()))
    require(not (half == other), "two walls holding two opaque samplers compare equal: identity equality was lost")
    with pytest.raises(ContentlessError):
        content_digest(half)
    summed = half + _spec(0.25)
    require(_is_condition(summed) and summed.source.sampled[0][0] == 0.5 and summed.cone_point == (0.25, 0.0), "adding a response to a framed wall lost a part")


# ── gate 6: displacements are values, not conditions ─────────────────────────


@pytest.mark.rests_on(_ALGEBRA)
def test_displacements_are_values_not_conditions() -> None:
    """``BoundaryConditions({1: 2.0 * specular(0.6)})`` raises ``TypeError`` naming ``WallDisplacement``;
    ``WallDisplacement.of((1.0, a), (-1.0, b))`` is a displacement for the operands of the two negative-coefficient
    tests of ``test_law_composition.py`` (``specular(1.0)``, ``white(1.0)``) AND for an in-``K`` pair (``of`` is
    ALWAYS a displacement); it is closed under ``+``, ``-``, scalar ``*``, is not a ``BoundaryCondition``, and does
    not equal (nor hash-equal) the condition of the same coordinates; a displacement plus a step that re-enters
    ``K`` is a ``WallCondition`` (G7).

    Reddens on: ``BoundaryConditions`` accepting a displacement; ``of`` returning a ``WallCondition`` when the
    combination lies in ``K``; a displacement whose digest equals its condition's.
    """
    over = 2.0 * _spec(0.6)
    with pytest.raises(TypeError, match="WallDisplacement"):
        bc.BoundaryConditions({1: over})
    spec1, white1 = _spec(1.0), _white(1.0)
    pinned = bc.WallDisplacement.of((1.0, spec1), (-1.0, white1))
    require(_is_displacement(pinned), "of((1, spec), (-1, white)) is not a displacement")
    require(pinned.network.cone_point == (1.0, -1.0), f"the signed coordinates read {pinned.network.cone_point}")
    require(_is_displacement(spec1 - white1), "the arithmetic spelling spec - white leaves K for these operands")
    require(_is_displacement(-(spec1 + white1)), "the arithmetic spelling -(spec + white) leaves K for these operands")
    require(_is_displacement(-spec1), "the arithmetic spelling -spec leaves K")

    inside = bc.WallDisplacement.of((0.5, _spec(0.5)), (0.25, _white(0.5)))
    require(_is_displacement(inside), "of(...) over an in-K combination is not a displacement")
    require(inside.network.cone_point == (0.25, 0.125), f"the in-K displacement reads {inside.network.cone_point}")
    twin = bc.WallCondition.from_amplitudes(0.25, 0.125)
    require(inside != twin and not (inside == twin), "a displacement equals the condition of its coordinates")
    require(content_digest(inside) != content_digest(twin), "a displacement and its condition share a digest")
    require(hash(inside) != hash(twin), "a displacement and its condition share a hash")
    require(not isinstance(inside, bc.BoundaryCondition), "a displacement is a BoundaryCondition")
    with pytest.raises(TypeError, match="WallDisplacement"):
        bc.BoundaryConditions({1: inside})

    closed = [pinned + pinned, pinned - pinned, 2.0 * pinned, -pinned, pinned + inside, pinned * 0.5]
    # `pinned - pinned` has all-zero coordinates: still a displacement? it lies in K, so by the value predicate it is a vacuum condition
    for result in closed:
        require(type(result) in (bc.WallCondition, bc.WallDisplacement), f"a displacement operator returned a {type(result).__name__}")
    require(_is_displacement(pinned + pinned) and _is_displacement(2.0 * pinned) and _is_displacement(-pinned), "a displacement operator left the displacements outright")
    require(_is_condition(pinned - pinned) and (pinned - pinned) == bc.WallCondition.vacuum(), "d - d re-enters K as vacuum()")

    step = -0.6 * _white(1.0)
    base = 0.6 * _spec(1.0) + 0.6 * _white(1.0)
    require(_is_displacement(base) and _is_displacement(step), "G7's operands are displacements")
    back = base + step
    require(_is_condition(back) and back == _spec(0.6), f"G7: the re-entry is {type(back).__name__}")
    entered = _spec(0.5) + bc.WallDisplacement.of((-0.25, _spec(1.0)))
    require(_is_condition(entered) and entered == _spec(0.25), "a WallCondition plus a displacement that re-enters K is a condition")


# ── gate 7: the two slots are disjoint (law L-F) ─────────────────────────────


@pytest.mark.rests_on(_PARSE, f"{_DECK}::test_parse_face_role")
def test_a_condition_is_not_a_deck_slot() -> None:
    """``BoundaryConditions({1: <any non-condition>})`` raises ``TypeError`` for a deck, a ``str`` and a number
    (the declared tag-to-condition map is not a role map); no value of the condition types is a face role
    (``parse_face_role`` refuses a ``WallCondition``); ``Binding`` has no mixed mapping, only ``decks`` and ``tags``.

    Reddens on: a copy of ``declared.py`` that accepts any value; a copy of the geometry parse that accepts a
    ``BoundaryCondition``.
    """
    import dataclasses

    population = [MirrorDeck(), TranslationDeck(), "vacuum", 3, 0.5, None, ("vacuum",)]
    for value in population:
        with pytest.raises(TypeError):
            bc.BoundaryConditions({1: value})
    require(len(population) >= 7, "the non-condition population shrank")
    with pytest.raises(TypeError):
        parse_face_role(bc.WallCondition.vacuum(), "a face role")
    names = [f.name for f in dataclasses.fields(bc.Binding)]
    require(names == ["decks", "tags", "declared", "symmetry"], f"Binding holds {names}: no mixed role mapping may exist")


@pytest.mark.rests_on(_PARSE, f"{_DECK}::test_a_deck_is_its_type_and_has_no_field")
def test_a_deck_is_not_an_operand() -> None:
    """``WallCondition.vacuum() + 3`` and ``0.5 * MirrorDeck()``, ``MirrorDeck() + WallCondition.vacuum()``,
    ``WallCondition.vacuum() + TranslationDeck()`` raise ``TypeError`` from Python's operator protocol: the
    traceback holds no frame of ``orpheus/`` (the deck defines no operator, the condition returns
    ``NotImplemented``). The positive control: ``0.5 * WallCondition.vacuum()`` works.

    Reddens on: giving ``MirrorDeck`` an ``__rmul__`` (``0.5 * MirrorDeck()`` stops raising), or
    ``WallCondition.__add__`` accepting any object.
    """
    require(0.5 * bc.WallCondition.vacuum() == bc.WallCondition.vacuum(), "control: the operator path itself works")
    cases = {
        "vacuum + 3": lambda: bc.WallCondition.vacuum() + 3,
        "0.5 * MirrorDeck": lambda: 0.5 * MirrorDeck(),
        "MirrorDeck + vacuum": lambda: MirrorDeck() + bc.WallCondition.vacuum(),
        "vacuum + TranslationDeck": lambda: bc.WallCondition.vacuum() + TranslationDeck(),
        "vacuum * MirrorDeck": lambda: bc.WallCondition.vacuum() * MirrorDeck(),  # type: ignore[operator]
        "vacuum - str": lambda: bc.WallCondition.vacuum() - "x",  # type: ignore[operator]
    }
    for label, thunk in cases.items():
        with pytest.raises(TypeError) as caught:
            thunk()
        require(no_orpheus_frame(caught.value), f"{label}: the TypeError came from a production frame, not from the operator protocol")


# ── gate 8: frames and terms ─────────────────────────────────────────────────


@pytest.mark.rests_on(_PARSE, _COLOURING)
def test_inflow_frames_and_terms() -> None:
    """The frame is the type: ``FaceFrameInflow(f)`` and the uniform number are ``face_local``, ``LabFrameInflow(f)``
    is not; a condition's ``invariant_under(identity)`` is True and under a non-identity ``FaceIsometry`` is False
    iff it holds a lab-frame term; ``InflowTerms.empty()`` is ``is_zero`` and face-local; the default body of the
    sampling verb is ``q = uniform * ones + sum_i w_i sampler_i.sample(space)`` (not W-normalised, one call per
    profile); a method-only test condition defines no operator.

    Reddens on: ``LabFrameInflow`` reporting ``face_local`` True; ``invariant_under`` ignoring ``is_identity``
    (it would then say True for a lab-frame term under the mirror); ``sample_terms`` dropping the uniform part or
    weighting by a quadrature.
    """
    sampler = OpaqueSampler()
    face = bc.WallCondition.framed(bc.FaceFrameInflow(sampler))
    lab = bc.WallCondition.framed(bc.LabFrameInflow(sampler))
    uniform = bc.WallCondition.inflow(1.0)
    require(face.face_local is True and uniform.face_local is True and bc.WallCondition.vacuum().face_local is True, "a face-frame, uniform or empty source is face-local")
    require(lab.face_local is False, "a lab-frame source is not face-local")
    require((lab + uniform).face_local is False and (face + uniform).face_local is True, "face_local does not follow the terms through +")
    require(bc.InflowTerms.empty().is_zero and bc.InflowTerms.empty().face_local, "the empty inflow is zero and face-local")
    require(not uniform.source.is_zero and not face.source.is_zero, "a nonzero source reads as zero")

    sym = colouring_symmetry(CoordSystem.CARTESIAN, (0.0, 2.0), (0.0, 1.0, 2.0), (0, 0), 0, 0)
    require(len(sym.elements) == 2, "the symmetric slab has the mirror")
    identity, mirror = sym.elements
    require(identity.is_identity and not mirror.is_identity, "the group's elements are the identity and the mirror")
    for cond in (face, uniform, lab, bc.WallCondition.vacuum(), _spec(0.5)):
        require(cond.invariant_under(identity) is True, "every condition is invariant under the identity")
    require(lab.invariant_under(mirror) is False, "a lab-frame term is not invariant under a non-identity isometry")
    for cond in (face, uniform, bc.WallCondition.vacuum(), _spec(0.5), _white(0.5)):
        require(cond.invariant_under(mirror) is True, "a face-local condition is invariant under every isometry")

    space = function_space(3)
    require(space.shape == (3,), f"control: the space has shape {space.shape}")
    only_uniform = bc.sample_terms(bc.InflowTerms(2.0, ()), space)
    require(np.array_equal(only_uniform, 2.0 * np.ones((3,))), f"a uniform inflow samples to {only_uniform}")
    a, b = OpaqueSampler(0.5), OpaqueSampler(4.0)
    terms = bc.InflowTerms(1.0, ((0.25, bc.FaceFrameInflow(a)), (0.5, bc.LabFrameInflow(b))))
    q = bc.sample_terms(terms, space)
    require(np.array_equal(q, (1.0 + 0.25 * 0.5 + 0.5 * 4.0) * np.ones((3,))), f"q = uniform + sum w_i sample_i read {q}")
    require(len(a.calls) == 1 and len(b.calls) == 1 and a.calls[0] is space, "each profile's sampler is called once, with the space")

    inert = local().Probe()
    with pytest.raises(TypeError):
        inert + inert  # type: ignore[operator]  # a method-only condition defines no operator
