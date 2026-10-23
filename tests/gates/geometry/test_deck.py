r"""The deck family: a closed sum type of the transformations that glue a face to its partner
(wave 2 of the layer DAG, WF2.2 unit L1; ``.claude/plans/wave2_boundary_carve.md`` sections 1.2, 2.1, 5.1).

The claims, each with its layer and its kind (``vv-principles``):

* gate 1, ``test_isometry_classification_is_closed``: THEOREM (Chasles, the six conjugacy classes of
  non-identity isometries of ``E^n`` for ``n <= 3``), stated over the type's own enumeration;
* gate 2, ``test_every_deck_match_ends_in_assert_never``: a STRUCTURE gate (an AST scan), the proof that a
  forgotten member is a type error and not a runtime refusal;
* gate 3, ``test_a_curved_face_admits_no_deck``: THEOREM (``W = N(G_c)/G_c`` is trivial on a 1-D radial
  chart), REFERENCE for the table;
* gates 4 and 5: the parse of a face role and the closure of a translation pairing (REFERENCE: the
  material-id twin rule, ``orpheus/data/materials.py:105``);
* gate 6, ``test_pairing_of_and_its_adjoint``: THEOREM (law L-A ii: the adjoint of a deck is the deck of its
  partner), over EVERY (deck, face) row;
* gate 10: identity of the content types (``MirrorDeck() == MirrorDeck()``, the roster rows live in
  ``test_content_identity_geometry.py``).

The module imports nothing from ``orpheus.geometry.deck`` at collection: the production module does not exist at
the base tip, and a collection error would hide which gate fails for which reason. Every row reads it through
:data:`deck` at run time, so each red names its own line (``ModuleNotFoundError`` at the row's first use, or the
row's own assertion once the module exists).

Mutations (each applied to a COPY of ``deck.py``, never by ``git checkout``; the reviewers' battery):

* gate 1: delete the ``GLIDE_REFLECTION`` row of ``DEFERRED_DECKS``; add a seventh ``IsometryClass`` member;
* gate 2: delete the ``case _: assert_never(deck)`` arm of ``pairing_of``;
* gate 3: return the Cartesian set for a radial axis;
* gate 4: accept ``bool`` in ``parse_face_role`` (the material-id twin rule);
* gate 5: delete the closure check of ``check_pairings_close``;
* gate 6: make ``translation_motion`` ignore its face; make ``FacePairing.inverse`` keep the motion; delete the
  involution invariant of ``FacePairing``.
"""

from __future__ import annotations

import ast
import dataclasses
import importlib
import re
from pathlib import Path
from types import SimpleNamespace
from typing import Any, get_args

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem
from orpheus.geometry.faces import FACE_NAMES, face_normal, face_opposite
from orpheus.geometry.transformation import RigidMotion
from orpheus.values.content import content_digest
from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/geometry/test_deck.py"
_REFLECTION_IS_AN_INVOLUTION = "tests/gates/geometry/test_transformation.py::test_c1_a_reflection_is_an_involution"
_ROOT = Path(__file__).resolve().parents[3]


class _Lazy:
    """The production module, imported at first attribute read (see the module docstring)."""

    def __init__(self, name: str) -> None:
        self._name = name

    def __getattr__(self, attr: str) -> Any:
        return getattr(importlib.import_module(self._name), attr)


deck: Any = _Lazy("orpheus.geometry.deck")

#: The four classes the deferral table names and does not mint, with the stem each future type's name carries.
_DEFERRED_STEMS = {
    "ROTATION": "RotationDeck",
    "GLIDE_REFLECTION": "GlideDeck",
    "SCREW_MOTION": "ScrewDeck",
    "ROTOREFLECTION": "RotoreflectionDeck",
}


# ── gate 1: the classification is closed ─────────────────────────────────────


@pytest.mark.rests_on(_REFLECTION_IS_AN_INVOLUTION)
def test_isometry_classification_is_closed() -> None:
    """Six classes (Chasles); 2, 4 and 6 of them in dimension 1, 2 and 3; every class is minted XOR deferred.

    Reddens on: a class missing from both the minted types and ``DEFERRED_DECKS`` (delete the
    ``GLIDE_REFLECTION`` row), a seventh member (the count and the dimension table), a class both minted and
    deferred (the disjointness), a deferral with no trigger (the machinery that mints it is unnamed).
    """
    classes = set(deck.IsometryClass)
    require(len(deck.IsometryClass) == 6, f"IsometryClass has {len(deck.IsometryClass)} members, not 6")
    require(
        {c.name for c in classes} == {"TRANSLATION", "REFLECTION", "ROTATION", "GLIDE_REFLECTION", "SCREW_MOTION", "ROTOREFLECTION"},
        f"the six classes are named {sorted(c.name for c in classes)}",
    )
    by_dim = deck.ISOMETRY_CLASSES_BY_DIMENSION
    require(set(by_dim) == {1, 2, 3}, f"the dimension table covers {sorted(by_dim)}")
    require([len(by_dim[d]) for d in (1, 2, 3)] == [2, 4, 6], f"sizes {[len(by_dim[d]) for d in (1, 2, 3)]}, not 2, 4, 6")
    require(by_dim[1] < by_dim[2] < by_dim[3], "the classes of dimension n are a proper subset of those of n + 1")
    require(by_dim[3] == frozenset(classes), "dimension 3 carries all six")
    require(
        by_dim[1] == frozenset({deck.IsometryClass.TRANSLATION, deck.IsometryClass.REFLECTION}),
        "dimension 1 is the translation and the reflection",
    )

    minted = {member.isometry_class for member in get_args(deck.DeckTransformation)}
    deferred = set(deck.DEFERRED_DECKS)
    require(len(get_args(deck.DeckTransformation)) == 2, "exactly two minted deck types")
    require(minted == {deck.IsometryClass.REFLECTION, deck.IsometryClass.TRANSLATION}, f"minted classes {minted}")
    require(minted.isdisjoint(deferred), f"a class is both minted and deferred: {minted & deferred}")
    require(minted | deferred == classes, f"classes neither minted nor deferred: {classes - (minted | deferred)}")
    require(deck.MirrorDeck.isometry_class is deck.IsometryClass.REFLECTION, "MirrorDeck is the reflection")
    require(deck.TranslationDeck.isometry_class is deck.IsometryClass.TRANSLATION, "TranslationDeck is the translation")

    for cls, row in deck.DEFERRED_DECKS.items():
        require(row.isometry_class is cls, f"the deferral row of {cls} is keyed by another class")
        require(isinstance(row.trigger, str) and row.trigger.strip(), f"{cls}: the deferral names no trigger")
        require(row.name.startswith(_DEFERRED_STEMS[cls.name]), f"{cls}: the future type is named {row.name!r}")
        require(dataclasses.is_dataclass(row) and row.__dataclass_params__.frozen, "a deferral row is a frozen value")


# ── gate 2: the match is exhaustive by type, not by runtime refusal ──────────

_DECK_ALIAS = re.compile(r"\bDeckTransformation\b")


def _matches_on_a_deck(source: str) -> tuple[int, list[str]]:
    """``(population, violations)`` of one source text: the ``match`` statements whose SUBJECT is a name annotated
    with the closed alias ``DeckTransformation`` (a parameter or a local annotation), and those among them whose
    last ``case`` is not a guard-free catch-all calling ``assert_never``.

    Outside the population by this predicate, stated: a subject typed ``object`` (``parse_face_role``, which ends
    in a ``TypeError``) and a subject typed ``FaceRole`` (the alias with ``int``).
    """
    population, violations = 0, []
    for fn in ast.walk(ast.parse(source)):
        if not isinstance(fn, ast.FunctionDef | ast.AsyncFunctionDef):
            continue
        typed = {
            a.arg
            for a in (*fn.args.posonlyargs, *fn.args.args, *fn.args.kwonlyargs)
            if a.annotation is not None and _DECK_ALIAS.search(ast.unparse(a.annotation))
        }
        for node in ast.walk(fn):
            if isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name) and _DECK_ALIAS.search(ast.unparse(node.annotation)):
                typed.add(node.target.id)
        for node in ast.walk(fn):
            if not (isinstance(node, ast.Match) and isinstance(node.subject, ast.Name) and node.subject.id in typed):
                continue
            population += 1
            last = node.cases[-1]
            catch_all = isinstance(last.pattern, ast.MatchAs) and last.pattern.pattern is None and last.guard is None
            call = last.body[-1].value if isinstance(last.body[-1], ast.Expr) else None
            names = {call.func.id} if isinstance(call, ast.Call) and isinstance(call.func, ast.Name) else set()
            names |= {call.func.attr} if isinstance(call, ast.Call) and isinstance(call.func, ast.Attribute) else set()
            if not (catch_all and "assert_never" in names):
                violations.append(f"{fn.name}:{node.lineno}")
    return population, violations


_PLANTED_BAD = '''
def f(deck: DeckTransformation):
    match deck:
        case MirrorDeck():
            return 1
        case TranslationDeck():
            return 2
'''
_PLANTED_GOOD = '''
def f(deck: DeckTransformation):
    match deck:
        case MirrorDeck():
            return 1
        case TranslationDeck():
            return 2
        case _:
            assert_never(deck)
'''
_PLANTED_OUTSIDE = '''
def parse(role: object, where: str):
    match role:
        case MirrorDeck():
            return role
        case _:
            raise TypeError(where)
'''
_PLANTED_RAISE_ARM = '''
def f(deck: DeckTransformation):
    match deck:
        case MirrorDeck():
            return 1
        case _:
            raise ValueError("unregistered")
'''


def test_every_deck_match_ends_in_assert_never() -> None:
    """Every ``match`` over a ``DeckTransformation`` in ``orpheus/`` ends in ``case _: assert_never(...)``, so a new
    deck type (the day ``RotationDeck`` is minted) fails pyright at every method that forgot it, never a runtime
    ``Refused(UNREGISTERED)``.

    Population at L1: ``pairing_of``'s one match (``[M]`` counted below, over every ``.py`` under ``orpheus/``).
    Controls (the scanner must read each shape): a planted match without the arm is a violation; a planted arm that
    raises instead of ``assert_never`` is a violation; the planted correct match is clean; a match on an
    ``object``-typed subject is outside the population (the exclusion of the claim).

    Reddens on: a copy of ``deck.py`` whose ``pairing_of`` lacks the arm (a violation), and at the base tip on
    the population (0 matches: nothing to be exhaustive over).
    """
    require(_matches_on_a_deck(_PLANTED_BAD) == (1, ["f:3"]), "control: the planted match without the arm is a violation")
    require(_matches_on_a_deck(_PLANTED_RAISE_ARM) == (1, ["f:3"]), "control: an arm that raises is not assert_never")
    require(_matches_on_a_deck(_PLANTED_GOOD) == (1, []), "control: the planted correct match is clean")
    require(_matches_on_a_deck(_PLANTED_OUTSIDE) == (0, []), "control: an object-typed subject is outside the population")

    population, violations, sites = 0, [], []
    files = sorted((_ROOT / "orpheus").rglob("*.py"))
    require(len(files) > 100 and (_ROOT / "orpheus/geometry/faces.py") in files, f"the scan read {len(files)} files")
    for path in files:
        n, bad = _matches_on_a_deck(path.read_text(encoding="utf-8"))
        population += n
        violations += [f"{path.relative_to(_ROOT)}:{b}" for b in bad]
        if n:
            sites.append(str(path.relative_to(_ROOT)))
    require(population >= 1, "no match over a DeckTransformation exists in orpheus/: nothing is exhaustive")
    require("orpheus/geometry/deck.py" in sites, f"pairing_of's match was not found; sites: {sites}")
    require(not violations, f"a match over a DeckTransformation lacks `case _: assert_never(...)`: {violations}")


# ── gate 3: a curved face admits no deck ─────────────────────────────────────


@pytest.mark.rests_on(_HERE + "::test_isometry_classification_is_closed")
def test_a_curved_face_admits_no_deck() -> None:
    """``W = N(G_c)/G_c`` is trivial for the 1-D cylinder and sphere, so the radial face admits no deck (a reflecting
    curved wall is a response, never a deck); every other axis admits the mirror and the translation.

    Reddens on: a copy returning the Cartesian set for a radial axis.
    """
    both = frozenset({deck.MirrorDeck, deck.TranslationDeck})
    require(deck.admitted_decks(CoordSystem.CYLINDRICAL, 0) == frozenset(), "a cylinder's radial face admits a deck")
    require(deck.admitted_decks(CoordSystem.SPHERICAL, 0) == frozenset(), "a sphere's radial face admits a deck")
    require(deck.admitted_decks(CoordSystem.CARTESIAN, 0) == both, "a slab's face admits the mirror and the translation")
    require(deck.admitted_decks(CoordSystem.CYLINDRICAL, 1) == both, "a cylinder's axial face admits both")
    require(deck.admitted_decks(CoordSystem.SPHERICAL, 1) == both, "a sphere's non-radial axis admits both")
    require(deck.admitted_decks(CoordSystem.CARTESIAN, 2) == both, "a slab's third axis admits both")
    require(isinstance(deck.admitted_decks(CoordSystem.CARTESIAN, 0), frozenset), "the admitted set is a frozenset")


# ── gate 4: the parse of a face role ─────────────────────────────────────────


@pytest.mark.rests_on(_HERE + "::test_isometry_classification_is_closed")
def test_parse_face_role() -> None:
    """A face role is a deck or an ``int`` tag; a ``bool`` is not a tag (the material-id twin rule).

    Reddens on: a copy accepting ``bool`` (``True`` would parse as tag 1); a copy that returns a numpy integer
    uncoerced (``type(...) is int`` fails); a copy that refuses a negative tag.
    """
    mirror, translation = deck.MirrorDeck(), deck.TranslationDeck()
    require(deck.parse_face_role(mirror, "w") is mirror, "a mirror deck passes through")
    require(deck.parse_face_role(translation, "w") is translation, "a translation deck passes through")
    for tag in (0, 1, 7, -3, np.int64(2), np.int32(5)):
        parsed = deck.parse_face_role(tag, "w")
        require(parsed == int(tag) and type(parsed) is int, f"{tag!r} parsed as {parsed!r} of {type(parsed).__name__}")

    refused: list[tuple[str, object]] = [
        ("True", True),
        ("False", False),
        ("np.bool_", np.True_),
        ("a BC", BC("reflective")),
        ("a BC-like object", SimpleNamespace(kind="reflective")),
        ("the string reflective", "reflective"),
        ("None", None),
        ("1.0", 1.0),
        ("np.float64", np.float64(1.0)),
        ("the class, not an instance", deck.MirrorDeck),
    ]
    for label, role in refused:
        with pytest.raises(TypeError, match="face_roles\\[3\\]") as info:
            deck.parse_face_role(role, "face_roles[3]")
        require(info.type is TypeError, f"{label}: refused with {info.type.__name__}")


# ── gate 5: a translation pairing closes ─────────────────────────────────────


@pytest.mark.rests_on(_HERE + "::test_parse_face_role")
def test_pairings_close() -> None:
    """A translation on a face requires a translation on the opposite face (an involution), else ``DeckError``.

    Reddens on: a copy whose ``check_pairings_close`` is the empty function (the first input passes).
    """
    t, m = deck.TranslationDeck(), deck.MirrorDeck()
    for roles in (
        {"xmin": t, "xmax": 0},
        {"xmin": t, "xmax": m},
        {"xmin": m, "xmax": t},
        {"xmin": 3, "xmax": t},
        {"ymin": t, "ymax": 1, "xmin": 0, "xmax": 0},
    ):
        with pytest.raises(deck.DeckError, match="where-token"):
            deck.check_pairings_close(roles, "where-token")
    for roles in (
        {"xmin": t, "xmax": t},
        {"xmin": m, "xmax": 0},
        {"xmin": m, "xmax": m},
        {"xmin": 1, "xmax": 1},
        {"ymin": t, "ymax": t, "xmin": 0, "xmax": 2},
    ):
        require(deck.check_pairings_close(roles, "where-token") is None, f"a closed pairing {roles} was refused")
    require(issubclass(deck.DeckError, ValueError), "DeckError is a ValueError")


# ── gate 6: pairing_of, and the adjoint of a deck is its partner's deck ──────


@pytest.mark.rests_on(_REFLECTION_IS_AN_INVOLUTION)
def test_pairing_of_and_its_adjoint() -> None:
    """Law L-A (ii): ``pairing_of(d, f).inverse == pairing_of(d, partner)`` for EVERY (deck, face, dimension) row.

    The translation is SIGNED BY THE FACE: it carries the partner onto the face, so ``g_xmin = -e_x`` and
    ``g_xmax = +e_x`` and the groupoid law closes. A face-blind unit translation (the ``PairedDeck.wrap`` spelling)
    makes the inverse of the ``xmin`` pairing the ``xmin`` pairing itself, and the row below false as written.

    Reddens on: a copy whose ``translation_motion`` ignores the face; a copy whose ``inverse`` keeps the motion;
    a copy whose mirror pairs a face with its opposite; a copy that drops the involution invariant.
    """
    mirror, translation = deck.MirrorDeck(), deck.TranslationDeck()
    rows = 0
    for dimension in (1, 2, 3):
        for face in FACE_NAMES:
            axis, sign = face_normal(face)
            if axis >= dimension:
                continue
            m = deck.pairing_of(mirror, face, dimension=dimension)
            require(m.face == face and m.partner == face and m.is_self_paired, f"{face}: a mirror pairs a face with itself")
            require(m.permutes_ordinates, f"{face}: a mirror permutes ordinates")
            expected = np.eye(dimension)
            expected[axis, axis] = -1.0
            require(np.array_equal(m.motion.linear, expected), f"{face}/{dimension}: the mirror's linear part is {m.motion.linear}")
            require(np.array_equal(m.motion.translation, np.zeros(dimension)), f"{face}: the mirror's plane position is stored")
            require(np.array_equal(m.motion.linear @ m.motion.linear, np.eye(dimension)), f"{face}: not an involution")
            require(m.inverse == m, f"{face}: the mirror is its own adjoint")

            t = deck.pairing_of(translation, face, dimension=dimension)
            partner = face_opposite(face)
            require(t.face == face and t.partner == partner and not t.is_self_paired, f"{face}: partner {t.partner}")
            require(not t.permutes_ordinates and t.motion.is_translation, f"{face}: a translation moves no direction")
            vector = t.motion.translation
            expected_sign = np.zeros(dimension)
            expected_sign[axis] = sign
            require(
                np.array_equal(np.sign(vector), expected_sign),
                f"{face}: the translation {vector} does not point along face_normal ({sign:+d} e_{axis})",
            )
            require(t.inverse == deck.pairing_of(translation, partner, dimension=dimension), f"{face}: the adjoint is not the partner's deck")
            require(t.inverse.inverse == t, f"{face}: the inverse is not an involution")
            rows += 2
    require(rows == 2 * (2 + 4 + 6), f"the loop read {rows} rows, expected {2 * (2 + 4 + 6)}")

    # the two spellings of the pole: the named rows of the spec
    xmax = deck.pairing_of(mirror, "xmax", dimension=3)
    require(xmax.is_self_paired and xmax.permutes_ordinates, "pairing_of(MirrorDeck, xmax) is self-paired")
    xmin_t = deck.pairing_of(translation, "xmin", dimension=3)
    xmax_t = deck.pairing_of(translation, "xmax", dimension=3)
    require(xmin_t.partner == "xmax" and np.array_equal(np.sign(xmin_t.motion.translation), [-1.0, 0.0, 0.0]), "xmin: -e_x")
    require(np.array_equal(np.sign(xmax_t.motion.translation), [1.0, 0.0, 0.0]), "xmax: +e_x")

    # the motion builders are the ONE spelling
    require(np.array_equal(deck.mirror_motion("y", 3).linear, np.diag([1.0, -1.0, 1.0])), "mirror_motion('y', 3)")
    require(np.array_equal(deck.mirror_motion("x", 1).linear, [[-1.0]]), "mirror_motion('x', 1)")
    require(
        np.array_equal(deck.translation_motion("xmin", 2).translation, -deck.translation_motion("xmax", 2).translation),
        "translation_motion is odd under the face's opposite",
    )
    for build, args in ((deck.mirror_motion, ("w", 3)), (deck.mirror_motion, ("z", 2)), (deck.translation_motion, ("zmax", 2))):
        with pytest.raises(ValueError):
            build(*args)

    # the parsed invariant: a face paired with itself carries an involution, never a quarter turn
    turn = RigidMotion(np.array([[0.0, -1.0], [1.0, 0.0]]))  # a rotation of order 4
    with pytest.raises(ValueError):
        deck.FacePairing("xmin", "xmin", turn)

    # the whole law, not the linear part: each of these constructed before the amendment
    identity2 = RigidMotion.identity(2)
    bad = {
        "self-paired by the identity": ("xmin", "xmin", identity2),
        "self-paired by a translation": ("xmin", "xmin", deck.translation_motion("xmin", 2)),
        "a mirror on the wrong axis": ("xmin", "xmin", deck.mirror_motion("y", 2)),
        "a mirror carrying a position": ("xmin", "xmin", RigidMotion.reflection(normal=np.array([1.0, 0.0]), offset=1.0)),
        "a name that is no face": ("banana", "kiwi", identity2),
        "distinct faces glued by the identity": ("xmin", "xmax", identity2),
        "a translation to a non-opposite face": ("xmin", "ymax", deck.translation_motion("xmin", 2)),
        "a translation with the wrong sign": ("xmin", "xmax", deck.translation_motion("xmax", 2)),
        "a translation along the wrong axis": ("xmin", "xmax", deck.translation_motion("ymin", 2)),
        "a quarter turn between distinct faces (a deferred deck)": ("xmin", "xmax", turn),
        "a quarter turn between adjacent faces (a deferred deck)": ("xmin", "ymin", turn),
        "the mirror of the face's own axis between opposite faces": ("xmin", "xmax", deck.mirror_motion("x", 2)),
    }
    for name, args in bad.items():
        with pytest.raises(deck.DeckError):
            deck.FacePairing(*args)
            raise AssertionError(f"{name}: constructed")  # noqa: B904 (a red)
    require(deck.FacePairing("xmin", "xmax", deck.translation_motion("xmin", 2)).partner == "xmax", "positive control: the signed translation")
    require(deck.FacePairing("ymax", "ymax", deck.mirror_motion("y", 2)).is_self_paired, "positive control: the mirror of its own axis")


# ── gate 10 (the part that needs no roster): a deck is its type ──────────────


@pytest.mark.rests_on(_HERE + "::test_isometry_classification_is_closed")
def test_a_deck_is_its_type_and_has_no_field() -> None:
    """``MirrorDeck() == MirrorDeck()``, one digest, one set member; unequal to ``TranslationDeck()``; no field.

    The face fixes the plane, the partner and the direction, so a deck stores nothing and the wrong-axis
    declaration is unspellable.

    Reddens on: a copy giving ``MirrorDeck`` an ``axis`` field (the field tuple), a copy with ``eq=True`` (the
    aggregator ``test_s5_9_every_content_class_takes_the_mixins_dunders`` reddens beside this row), a copy whose
    two types share a qualname-free digest.
    """
    for cls in (deck.MirrorDeck, deck.TranslationDeck):
        require(dataclasses.fields(cls) == (), f"{cls.__name__} carries fields {dataclasses.fields(cls)}")
        require(cls() == cls(), f"{cls.__name__}() != {cls.__name__}()")
        require(content_digest(cls()) == content_digest(cls()), f"{cls.__name__}: the digest is not one value")
        require(hash(cls()) == hash(cls()) and len({cls(), cls()}) == 1, f"{cls.__name__}: not one set member")
    require(deck.MirrorDeck() != deck.TranslationDeck(), "a mirror equals a translation")
    require(content_digest(deck.MirrorDeck()) != content_digest(deck.TranslationDeck()), "one digest for two decks")
    require(len({deck.MirrorDeck(), deck.TranslationDeck()}) == 2, "a set merged the two decks")
