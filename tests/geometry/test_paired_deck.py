r"""The genuinely-paired deck element's construction invariant (step G6.3-7).

:class:`~orpheus.geometry.boundary.PairedDeck` is :math:`G` for a face paired
with a **DISTINCT** face — the second of Poincaré's two cases and the exact
complement of :class:`~orpheus.geometry.boundary.SelfPairedDeck`. It retires
``SpatialWrap``, whose ``axis: str`` could name neither a lattice vector nor a
rotation angle.

The one claim that needs BOTH types to be visible at once
---------------------------------------------------------

The two guards are complements of a single criterion — the **fixed subspace**:

.. math::

   \text{self-paired} \iff Q \text{ linear} \wedge
   \dim \mathrm{Fix}(Q) \ge d - 1

so a motion is admitted by *exactly one* of the two types, in every dimension.
No test living in one type's module can see that: a guard that drifts into
**overlap** (both accept) or into a **gap** (neither accepts) leaves every
single-type positive and negative leg green. :class:`TestTheGuardsPartitionE3`
is the row that closes it, and it is the reason this module imports both types.

⭐ Involution is NECESSARY but NOT SUFFICIENT for self-pairing
--------------------------------------------------------------

This is the load-bearing row of the taxonomy, and it is the one an
``element_order()``-shaped guard gets backwards. :math:`E(3)` has FOUR
involution families and the last two map a face to its **opposite**:

===============  =====  ====  ===  ==============  ============
element          order   det   fix  pairing         type
===============  =====  ====  ===  ==============  ============
identity             1    +1    3  self-paired     SelfPaired
reflection           2    −1    2  self-paired     SelfPaired
**half-turn**    **2**    +1    1  face → opposite **Paired**
**inversion**    **2**    −1    0  face → opposite **Paired**
===============  =====  ====  ===  ==============  ============

`[M]` 2026-08-07 against the live guards, reproduced by
:class:`TestTheGuardTable` below. The half-turn row carries its own
``element_order() == 2`` activation guard: if a future edit made the fixture
a non-involution, the row would still pass while proving nothing, so the
witness asserts it IS an involution before asserting it is admitted here.

What is NOT gated here, deliberately
------------------------------------

* **The realized arrow.** Whether the deck element's *action on the trace* is
  derived from the motion — and in which direction the gather runs — is the
  SN realization tier: ``tests/sn/operators/test_deck_kernel.py``.
* **``catches("ERR-NNN")``.** ERR-073 (a non-injective node match certified as
  a permutation) and ERR-074 (an uncertified reflection table) are defects of
  the node-matching layer, not of this type. A marker here would be a
  same-area tag, which ``vv-principles`` forbids — a ``catches`` marker is a
  coverage CLAIM. The teeth for both live in the SN module.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.boundary import (
    BoundaryGeometryMap,
    PairedDeck,
    PeriodicBoundary,
    SelfPairedDeck,
)
from orpheus.geometry.boundary._errors import BoundaryError
from orpheus.geometry.transformation import RigidMotion
from orpheus.numerics.face_layout import AXIS_NAMES, face_opposite

pytestmark = pytest.mark.foundation


# ---------------------------------------------------------------------------
# The motion inventory.
#
# Built inside lambdas, never at import: a parametrize list that CONSTRUCTS
# the subject runs at COLLECTION time, so a regression turns every gate in the
# file into one collection error whose diagnosis is "1 error" rather than the
# named invariant that broke. (The sibling ``test_self_paired_deck.py``
# measured exactly that and made its table lazy for the same reason.)
# ---------------------------------------------------------------------------

#: ``name -> (factory, admitted_by_paired_deck)``. Every entry is an element of
#: :math:`E(3)`; the flag is the *expected* verdict of ``PairedDeck``'s guard,
#: and by the partition law it is the negation of ``SelfPairedDeck``'s.
_MOTIONS: dict[str, tuple] = {
    # --- self-paired: linear AND fixes a hyperplane pointwise -------------
    "identity": (lambda: RigidMotion.identity(3), False),
    "mirror-x": (lambda: RigidMotion.reflection(normal=[1.0, 0.0, 0.0]), False),
    "mirror-diagonal": (
        lambda: RigidMotion.reflection(normal=[1.0, 1.0, 0.0]), False,
    ),
    # --- genuinely paired --------------------------------------------------
    "rotation-90": (
        lambda: RigidMotion.rotation_about_axis(
            axis=[0.0, 0.0, 1.0], angle=np.pi / 2
        ),
        True,
    ),
    # ⭐ an INVOLUTION that is still face-swapping — the row that refutes an
    #    ``element_order`` guard.
    "half-turn": (
        lambda: RigidMotion.rotation_about_axis(
            axis=[0.0, 0.0, 1.0], angle=np.pi
        ),
        True,
    ),
    # ⭐ the second involution family: fixes only the origin.
    "inversion": (lambda: RigidMotion.inversion(3), True),
    "translation-x": (
        lambda: RigidMotion.translation_by([1.0, 0.0, 0.0]), True,
    ),
    "translation-diagonal": (
        lambda: RigidMotion.translation_by([1.0, 1.0, 0.0]), True,
    ),
    # A glide FIXES a 2-plane setwise (so it clears the fixed-set clause) and
    # is refused by SelfPairedDeck on LINEARITY alone — the partner witness to
    # the inversion, which is refused on the fixed-set clause alone.
    "glide": (
        lambda: RigidMotion.translation_by([0.0, 1.0, 0.0])
        @ RigidMotion.reflection(normal=[1.0, 0.0, 0.0]),
        True,
    ),
}

#: Motions in dimensions other than 3, for the partition law only. The
#: criterion is written ``dim Fix >= d - 1``, so it is dimension-generic and a
#: hard-coded ``2`` would pass every 3-D row in this file.
_LOW_DIM_MOTIONS: dict[str, tuple] = {
    "1d-identity": (lambda: RigidMotion.identity(1), False),
    "1d-reflection": (lambda: RigidMotion.reflection(normal=[1.0]), False),
    "1d-translation": (lambda: RigidMotion.translation_by([1.0]), True),
    "2d-identity": (lambda: RigidMotion.identity(2), False),
    "2d-reflection": (
        lambda: RigidMotion.reflection(normal=[1.0, 0.0]), False,
    ),
    "2d-half-turn": (lambda: RigidMotion(-np.eye(2)), True),
    "2d-translation": (
        lambda: RigidMotion.translation_by([1.0, 0.0]), True,
    ),
}


def _admits(cls, motion: RigidMotion) -> bool:
    """``True`` iff ``cls`` constructs from ``motion`` without raising."""
    try:
        cls(motion)
    except (TypeError, ValueError):
        return False
    return True


# ============================================================================
# P1 — the guard table, both directions
# ============================================================================


class TestTheGuardTable:
    r"""Which motions :class:`PairedDeck` admits, and which it refuses.

    Positive AND negative legs (``vv-principles`` anti-pattern #11): a guard
    tested only against broken instances validates the *raising*, never the
    *claim*. Here the claim is a partition of :math:`E(3)`, so both halves
    are load-bearing and the refusals must name a REASON, not merely raise.
    """

    @pytest.mark.parametrize(
        "name", [n for n, (_, ok) in _MOTIONS.items() if ok]
    )
    def test_a_genuinely_paired_element_is_ADMITTED(self, name: str) -> None:
        """Rotation, half-turn, inversion, translation and glide construct."""
        motion = _MOTIONS[name][0]()
        deck = PairedDeck(motion)
        assert deck.motion is motion, (
            "the deck stores the element it was given — a copy would make "
            "identity comparisons against the caller's motion fail silently"
        )

    @pytest.mark.parametrize(
        "name", [n for n, (_, ok) in _MOTIONS.items() if not ok]
    )
    def test_a_SELF_paired_element_is_REFUSED(self, name: str) -> None:
        r"""The identity and every mirror belong to the other type.

        ``match="fix"`` keys the refusal to the fixed-subspace criterion. A
        bare ``pytest.raises(ValueError)`` would also be satisfied by an
        unrelated failure inside ``RigidMotion``, which is the ERR-class this
        campaign has hit twice (lessons: a refusal gate without a MESSAGE leg
        is teeth-less).
        """
        motion = _MOTIONS[name][0]()
        with pytest.raises(ValueError, match="fixes a subspace"):
            PairedDeck(motion)

    def test_the_HALF_TURN_is_an_involution_and_still_admitted(self) -> None:
        r"""⭐ The row that refutes an ``element_order`` guard.

        A guard spelled ``element_order() not in (1, 2)`` would refuse the
        half-turn and the inversion — the two elements this type exists to
        carry — while admitting every mirror, which belongs to the other
        type. It would be wrong in BOTH directions at once.

        The ``element_order() == 2`` assertion is this row's activation
        guard: without it, a fixture edit that made the witness a
        non-involution would leave the gate green and vacuous.
        """
        for name in ("half-turn", "inversion"):
            motion = _MOTIONS[name][0]()
            assert motion.element_order() == 2, (
                f"{name} must be an involution for this row to be the "
                f"discriminating one it claims to be; got order "
                f"{motion.element_order()}"
            )
            deck = PairedDeck(motion)          # MUST NOT raise
            assert deck.motion.fixed_subspace_dimension < 2, (
                f"{name} must fix LESS than a hyperplane, or it is not a "
                f"face-swapping element and this row proves nothing"
            )

    def test_the_two_clauses_are_separately_witnessed(self) -> None:
        r"""The guard is a CONJUNCTION, so each clause needs its own witness.

        ``PairedDeck`` refuses ``is_linear ∧ fix ≥ d−1``. The glide fails the
        first conjunct while passing the second (``fix == 2``); the inversion
        passes the first while failing the second (``fix == 0``). Both are
        admitted here, and each is admitted *for a different reason* — which
        is what proves the guard is not accidentally a one-clause test.
        """
        glide = _MOTIONS["glide"][0]()
        assert not glide.is_linear
        assert glide.fixed_subspace_dimension == 2, (
            "the glide must CLEAR the fixed-set clause, or it is not the "
            "witness this row needs"
        )
        inversion = _MOTIONS["inversion"][0]()
        assert inversion.is_linear
        assert inversion.fixed_subspace_dimension == 0
        PairedDeck(glide)
        PairedDeck(inversion)


# ============================================================================
# P2 / P11 — the partition law, and the CONTROL through the other type
# ============================================================================


class TestTheGuardsPartitionE3:
    r"""⭐ Every rigid motion is admitted by **exactly one** deck type.

    The law that neither type's own module can state, because each type's
    tests are consistent with its own guard whatever that guard says. Swept
    over dimensions 1, 2 and 3: the criterion is ``dim Fix >= d − 1``, so a
    hard-coded ``2`` would satisfy every 3-D row in this file and fail only
    here.

    `[M]` 2026-08-07, run over this module **and**
    ``test_self_paired_deck.py`` together — the exclusivity claim, bounded
    honestly rather than asserted:

    ==========================================  ==================  =========================
    mutation                                    reds HERE           reds in test_self_paired_deck
    ==========================================  ==================  =========================
    ``PairedDeck`` threshold ``fix >= d``       P1 (2) + **this (4)**  **0** — invisible there
    (OVERLAP: both types accept the mirror)
    ``SelfPairedDeck`` guard fully relaxed      **this only (9)**   5
    ``SelfPairedDeck`` guard → "is it an        **this only (3)**   2
    involution?" (the converse)
    ==========================================  ==================  =========================

    So the precise claim is two-part: this class is the **only** catcher for
    the overlap/gap class anywhere in the tree, and within THIS module it is
    the only row that can see the sibling type's guard move at all. A weaker
    reading — "no other test in the tree sees a sibling drift" — is false and
    is not what this class is credited with.
    """

    @pytest.mark.parametrize(
        "name", list(_MOTIONS) + list(_LOW_DIM_MOTIONS)
    )
    def test_exactly_one_type_admits_each_motion(self, name: str) -> None:
        """No overlap, no gap — the two guards are complements."""
        factory, paired_expected = (_MOTIONS | _LOW_DIM_MOTIONS)[name]
        motion = factory()
        paired = _admits(PairedDeck, motion)
        self_paired = _admits(SelfPairedDeck, motion)
        assert paired != self_paired, (
            f"{name}: PairedDeck admits={paired}, SelfPairedDeck "
            f"admits={self_paired}. The two guards are complements of ONE "
            f"criterion (dim Fix >= d-1 for a LINEAR element), so a motion "
            f"admitted by both is an overlap and one admitted by neither is "
            f"a gap — either way a deck element has no home or two homes."
        )
        assert paired is paired_expected, (
            f"{name}: PairedDeck admits={paired}, expected "
            f"{paired_expected}. The taxonomy moved, not just the guard."
        )

    def test_the_partition_is_not_vacuous_in_either_direction(self) -> None:
        """Both sides of the partition are populated by the fixture set.

        A partition law is trivially satisfied by a table that is entirely
        one-sided; this row pins that the inventory exercises both.
        """
        table = _MOTIONS | _LOW_DIM_MOTIONS
        paired = {n for n, (_, ok) in table.items() if ok}
        self_paired = set(table) - paired
        assert len(paired) >= 3 and len(self_paired) >= 3, (
            f"the motion inventory must populate both sides of the "
            f"partition; got {len(paired)} paired / {len(self_paired)} "
            f"self-paired"
        )


# ============================================================================
# P3 / P4 / P5 — the ``wrap`` constructor
# ============================================================================


class TestTheWrapConstructor:
    r""":meth:`PairedDeck.wrap` — the periodic translation, named by an axis.

    The axis LETTER survives as a *constructor argument* (ergonomics) while
    the stored datum is the MOTION. That is the whole retirement: a letter is
    fine as a way to say which unit vector, and useless as the field a
    realizer reads.
    """

    @pytest.mark.parametrize("axis", list(AXIS_NAMES))
    def test_wrap_stores_the_UNIT_translation_along_the_axis(
        self, axis: str
    ) -> None:
        r"""``t == ê_axis``, exactly — and the length is deliberately 1.

        The wrap's *direction* is intrinsic to the law; its *length* is the
        domain extent, i.e. configuration. And the length is bit-identically
        invisible to realization, because ``on_directions`` drops the
        translation — so a stored physical length would be a Mode-12
        designed-green field. ``array_equal``: these are constructed
        constants, not measured quantities (lessons: bit-exactness is earned
        per law; a tolerance here would admit the bug).
        """
        deck = PairedDeck.wrap(axis=axis)
        expected = np.eye(3)[AXIS_NAMES.index(axis)]
        np.testing.assert_array_equal(deck.motion.translation, expected)
        np.testing.assert_array_equal(deck.motion.linear, np.eye(3))

    def test_the_offset_is_invisible_which_is_WHY_the_unit_is_enough(
        self,
    ) -> None:
        r"""`[M]` The Mode-12 argument, asserted rather than asserted-about.

        Two wraps of different length act identically on DIRECTIONS, bit for
        bit. If this ever fails, the length is observable after all and the
        unit-vector rationale needs rewriting, not the assertion.
        """
        unit = PairedDeck.wrap(axis="x").motion
        long = RigidMotion.translation_by([17.0, 0.0, 0.0])
        probe = np.array([[0.3, -0.7, 0.6], [-1.0, 0.0, 0.0]])
        np.testing.assert_array_equal(
            unit.on_directions(probe), long.on_directions(probe)
        )
        assert unit.is_translation and long.is_translation

    def test_a_misnamed_axis_is_refused_at_CONSTRUCTION(self) -> None:
        """Not at realization, holding a face layout it should never reach."""
        with pytest.raises(ValueError, match="axis must be one of"):
            PairedDeck.wrap(axis="q")

    def test_an_axis_outside_the_dimension_is_refused(self) -> None:
        """A z-wrap has no meaning on a 1-D slab."""
        with pytest.raises(ValueError, match="out of range"):
            PairedDeck.wrap(axis="z", dimension=1)

    @pytest.mark.parametrize(
        "kwargs, pattern",
        [
            ({"axis": "q"}, "axis must be one of"),
            ({"axis": "z", "dimension": 1}, "out of range"),
        ],
    )
    def test_the_refusal_SHAPE_matches_the_self_paired_sibling(
        self, kwargs: dict, pattern: str
    ) -> None:
        r"""⭐ Both constructors refuse the same inputs, with the same messages.

        The two axis-resolving classmethods (:meth:`PairedDeck.wrap` and
        :meth:`SelfPairedDeck.mirror`) are the tier's only two entry points
        that take a letter. If one validates and the other does not, a caller
        learns the convention from whichever it happened to use — which is
        the two-sources-of-truth this campaign closes everywhere else.

        Compared as a pair rather than duplicated per class, so the row reds
        when they DIVERGE rather than needing to be edited in two files.
        """
        for factory in (PairedDeck.wrap, SelfPairedDeck.mirror):
            with pytest.raises(ValueError, match=pattern):
                factory(**kwargs)

    def test_a_non_motion_is_refused_with_a_type_error(self) -> None:
        """The deck element is the transformation, never a name for one.

        Passing the retired ``axis: str`` value directly is the exact
        regression this type exists to prevent, so the witness IS ``"x"``.
        """
        with pytest.raises(TypeError, match="RigidMotion"):
            PairedDeck("x")  # type: ignore[arg-type]

    def test_two_wraps_of_the_same_axis_compare_EQUAL(self) -> None:
        r"""Value semantics, which ``_LAW_SPECS`` in
        ``test_boundary_factors.py`` compares against with ``==``.

        Not free: :class:`RigidMotion` is declared ``eq=False`` and supplies
        a bit-exact ``__eq__`` of its own. Two separately-constructed unit
        translations are bit-identical, so the equality holds — and a wrap
        along a DIFFERENT axis must not.
        """
        assert PairedDeck.wrap(axis="x") == PairedDeck.wrap(axis="x")
        assert PairedDeck.wrap(axis="x") != PairedDeck.wrap(axis="y")
        assert PairedDeck.wrap(axis="x") != SelfPairedDeck.mirror(axis="x")


# ============================================================================
# P6 / P7 / P8 — ``domain_face``: the partner, derived from the motion
# ============================================================================


class TestDomainFaceNamesThePartner:
    r""":meth:`PairedDeck.domain_face` — the only
    :class:`BoundaryGeometryMap` whose answer is not ``face``.

    Derived from the motion, never stored: which face is the partner depends
    on where the law is installed (configuration), while "wrap along x" is
    intrinsic. That is B0's rule, and it is why the answer is a METHOD taking
    the installation face.
    """

    @pytest.mark.parametrize("axis", list(AXIS_NAMES))
    def test_the_partner_is_the_OPPOSITE_face_on_axis(self, axis: str) -> None:
        r"""``domain_face(f) == face_opposite(f)`` for both on-axis faces.

        The reference is :func:`~orpheus.numerics.face_layout.face_opposite`,
        the tree's single home for that convention — so this row pins the
        DERIVATION (translation ⟹ opposite face) rather than re-spelling a
        name table that could drift from it.
        """
        deck = PairedDeck.wrap(axis=axis)
        for suffix in ("min", "max"):
            face = f"{axis}{suffix}"
            assert deck.domain_face(face) == face_opposite(face)
            assert deck.domain_face(face) != face, (
                "a genuine pair's domain is a DIFFERENT face — returning "
                "`face` is the self-paired answer and would realize the wrap "
                "as a bare identity on the wrong half-trace (ERR-041 class)"
            )

    @pytest.mark.parametrize(
        "wrap_axis, face",
        [("x", "ymin"), ("x", "zmax"), ("y", "xmin"), ("z", "ymax")],
    )
    def test_an_OFF_AXIS_face_is_refused_naming_both_axes(
        self, wrap_axis: str, face: str
    ) -> None:
        r"""A wrap along ``y`` installed on an ``x`` face identifies nothing.

        The refusal message must name the wrap's axis AND the face's axis:
        the caller's error is that the two disagree, and a message naming
        only one of them cannot say which to change. ``exc.value.law`` is
        asserted because an untyped raise loses the attribution the whole
        :class:`BoundaryError` channel exists to carry.
        """
        deck = PairedDeck.wrap(axis=wrap_axis)
        with pytest.raises(BoundaryError) as excinfo:
            deck.domain_face(face)
        message = str(excinfo.value)
        assert f"axis={wrap_axis!r}" in message and face in message, (
            f"the refusal must name the wrap's axis and the installed face; "
            f"got: {message}"
        )
        assert excinfo.value.law == "periodic"

    def test_an_unknown_face_name_is_refused_as_a_BoundaryError(self) -> None:
        """Not as the raw ``ValueError`` from the face-name parser.

        The caller asked a boundary question; the diagnosis must arrive on
        the boundary channel with its ``law`` attribution, not as an
        unattributed complaint from a name table two layers down.
        """
        with pytest.raises(BoundaryError, match="cannot name a partner"):
            PairedDeck.wrap(axis="x").domain_face("bogus")

    def test_a_ROTATION_has_no_face_name_derivation_and_says_so(self) -> None:
        r"""⭐ #178, not a silent wrap.

        A sector rotation is a legitimate deck ELEMENT — this type admits it
        — but an axis-aligned box has no face-NAME derivation for it, so
        ``domain_face`` must refuse rather than guess. The refusal is keyed
        to *"not a translation"*: a blanket "the message mentions #178" would
        pin the same reason on the diagonal-translation row below, whose
        defect is different.
        """
        rotation = PairedDeck(_MOTIONS["rotation-90"][0]())
        with pytest.raises(BoundaryError, match="not a translation"):
            rotation.domain_face("xmin")

    def test_a_DIAGONAL_translation_has_no_face_name_derivation_either(
        self,
    ) -> None:
        r"""The second #178 arm, with its own distinct reason.

        A translation along ``(1, 1, 0)`` IS a translation — it clears the
        first precondition — and still carries no face of an axis-aligned box
        onto another. If the ``is_translation`` precondition were dropped, a
        rotation would reach the single-axis test, whose ``flatnonzero``
        finds a zero translation and would blame *this* reason for *that*
        motion. Two distinct ``match=`` strings are what keep the two
        diagnoses apart.
        """
        diagonal = PairedDeck(_MOTIONS["translation-diagonal"][0]())
        with pytest.raises(BoundaryError, match="not along a single axis"):
            diagonal.domain_face("xmin")

    def test_domain_face_is_an_INVOLUTION_on_the_pair(self) -> None:
        """``domain_face(domain_face(f)) == f`` — the pair is a pair.

        Cheap, and it is the structural statement that the two faces are
        each other's partner rather than both pointing at a third.
        """
        deck = PairedDeck.wrap(axis="x")
        assert deck.domain_face(deck.domain_face("xmin")) == "xmin"
        assert deck.domain_face(deck.domain_face("xmax")) == "xmax"


# ============================================================================
# P9 / P10 — the properties the motion makes DERIVABLE
# ============================================================================


class TestThePropertiesAreDerivedFromTheMotion:
    r"""``permutes_ordinates`` and ``is_adjointable`` — read off, not declared.

    ⭐ ``permutes_ordinates`` is the property a hard-coded ``return False``
    gets right for every law the tree SHIPS and wrong for the latent
    rotational sector BC (#178). The wrap alone cannot discriminate the two
    implementations, which is why every row here carries a rotation witness.
    """

    @pytest.mark.parametrize(
        "name, expected",
        [
            ("translation-x", False),
            ("translation-diagonal", False),
            ("glide", False),          # a glide's linear part is a REFLECTION…
            ("rotation-90", True),
            ("half-turn", True),
            ("inversion", True),
        ],
    )
    def test_permutes_ordinates_answers_the_motion(
        self, name: str, expected: bool
    ) -> None:
        r"""``permutes_ordinates ⟺ the linear part is not the identity``.

        ⚠ The ``glide`` row is deliberately ``False``-expected and is the
        one that would catch a naive ``not is_linear`` spelling: a glide has
        a reflection as its linear part, so it DOES move angle, and
        ``is_translation`` (linear part ``== I``) is the predicate that gets
        it right where ``is_linear`` (translation part ``== 0``) does not.
        The two predicates are complements of each other on the diagonal and
        disagree exactly here.
        """
        motion = _MOTIONS[name][0]()
        deck = PairedDeck(motion)
        # ⚠ the glide's linear part is a reflection => it DOES permute.
        truth = not np.array_equal(motion.linear, np.eye(motion.dimension))
        assert deck.permutes_ordinates is truth, (
            f"{name}: permutes_ordinates={deck.permutes_ordinates} but the "
            f"linear part is "
            f"{'the identity' if not truth else 'non-trivial'} — the "
            f"property must be READ OFF the motion, not declared beside it"
        )
        if name != "glide":
            assert deck.permutes_ordinates is expected

    def test_the_wrap_does_not_permute_and_the_rotation_does(self) -> None:
        r"""The discriminating pair, stated as one assertion.

        A ``return False`` implementation passes the first and fails the
        second; a ``return True`` implementation does the reverse. Only the
        pair pins the derivation.
        """
        assert PairedDeck.wrap(axis="x").permutes_ordinates is False
        assert PairedDeck(
            _MOTIONS["rotation-90"][0]()
        ).permutes_ordinates is True

    @pytest.mark.parametrize("name", ["translation-x", "rotation-90"])
    def test_it_conforms_to_the_geometry_map_protocol(self, name: str) -> None:
        """Runtime structural conformance, and the adjointability theorem.

        ``is_adjointable`` is a THEOREM for a genuine deck transformation —
        the composition operator of a bijection is invertible and
        measure-preservation makes that inverse the transpose — so ``True``
        here is not a declaration that could drift. (Its predecessor
        ``SpatialWrap`` answered ``False`` until B3.4c, reporting an
        implementation gap in a slot whose contract is a property of the MAP.)
        """
        deck = PairedDeck(_MOTIONS[name][0]())
        assert isinstance(deck, BoundaryGeometryMap)
        assert deck.is_adjointable is True


# ============================================================================
# P-wire — the production law reaches this type
# ============================================================================


@pytest.mark.parametrize("axis", list(AXIS_NAMES))
def test_periodic_boundary_hands_out_a_paired_deck(axis: str) -> None:
    r"""``PeriodicBoundary(axis).geometry_map`` IS the wrap, valued.

    Compared by VALUE, not by ``isinstance``: since the tier has two deck
    types and one of them is now parameterised by a motion, an isinstance
    check would be blind to the AXIS entirely — the exact blindness
    ``test_boundary_factors.py`` records for the collapsed self-paired half.
    """
    law = PeriodicBoundary(axis=axis)
    assert law.geometry_map == PairedDeck.wrap(axis=axis)
    assert law.geometry_map.domain_face(f"{axis}min") == f"{axis}max"


# ============================================================================
# P12 — the retirement
# ============================================================================


@pytest.mark.parametrize(
    "module_name",
    ["orpheus.geometry.boundary", "orpheus.geometry.boundary._factors"],
)
def test_spatial_wrap_is_RETIRED_with_no_alias(module_name: str) -> None:
    r"""``SpatialWrap`` is gone — from the package facade AND its defining
    module.

    The only gate that can see a compatibility shim being re-added, which the
    campaign's retirement rule forbids ("deprecation aliases live for one
    merge cycle only", and this one is authorised for zero). Both modules are
    swept because a shim can be installed at either level and the facade's
    ``__all__`` would not necessarily advertise it.

    ``hasattr`` rather than a ``from … import`` inside ``pytest.raises``: the
    ``from``-form resolves through the same ``getattr`` (module-level
    ``__getattr__`` included), so it tests nothing extra — and it is a
    statically-unresolvable import that a type checker must flag.
    """
    import importlib

    module = importlib.import_module(module_name)
    assert not hasattr(module, "SpatialWrap"), (
        f"{module_name} still exposes SpatialWrap; it was retired onto "
        f"PairedDeck at G6.3 step 7, and a surviving attribute means an alias "
        f"was re-added — the two-spellings state the retirement removed"
    )
    assert "SpatialWrap" not in getattr(module, "__all__", ()), (
        f"{module_name}.__all__ still advertises the retired name"
    )
