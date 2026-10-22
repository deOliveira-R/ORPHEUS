r"""CS5 — an axis can name the object that minted it (the provenance doctrine).

One doctrine, one file (the ``test_space_of_axes.py`` file-level-doctrine
precedent): since CS5 (2026-08-29) the mint routes THROUGH the generating
object (``Axis.of_measure(measure, label)`` / ``Axis.of_quadrature(quad,
label)``), which records itself on the axis: provenance, **never identity**.
Since S7 of the numerics untangling (ruling 64) the record is
``Axis.origin``, a closed sum (``FromMeasure``, ``FromQuadrature``,
``FromBasis``, ``FromFrame``, ``Declared``) read by ``match``, and the axis
keeps its support and atoms as content (ruling 71), so an axis is no longer
a forgetful image of its measure: what the origin adds is only what the
measure cannot answer (the level fibration).

The gates here are the landed machinery's ONLY witnesses — ``[M]`` a
three-arm in-process mutation battery over the 184-test anchor set found
**0 genuine catchers** before this module landed (the one red was a
cross-process harness artefact; ``scratch/cs5_verification_plan.md`` §1).

G5 (the generator-less refusal) and G7 (the route keystone) deliberately do
NOT live here yet: no solve-time consumer reads a generator until the
streaming plan's P4-remainder phase re-points the producer, and a gate that
lands before the case it catches ships green and unfalsifiable
(plan-authoring §6c). Their specs are in the verification plan.
"""

from __future__ import annotations

import numpy as np
import numpy.testing as npt
import pytest

from orpheus.geometry.manifold import RealSpace
from orpheus.numerics.axis import Axis, BasisKind, Declared, FromMeasure, FromQuadrature, quadrature_refusal
from orpheus.numerics.quadrature.directional import Quadrature
from orpheus.numerics.space import FunctionSpace

pytestmark = pytest.mark.foundation


def _require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


def _point() -> Axis:
    """The quotient spatial point (counting weight — the density convention)."""
    return Axis("spatial", (1,), origin=Declared(BasisKind.NODAL))


# The shipped angular family — EXHAUSTIVE by construction, not a ladder
# (vv-principles #31's finite-roster corollary): these are the FIVE
# ``Quadrature`` classmethod factories — [M] 2026-08-29 (archivist census)
# ``[n for n, v in vars(Quadrature).items() if isinstance(v, classmethod)]``
# = gauss_legendre, lebedev, level_symmetric, product, folded_product; the
# roster's first cut omitted folded_product, the sigma_y-folded CARRYING
# rule the curvilinear MMS builders default to — i.e. the member with the
# richest ``level_indices``, on exactly the axis this roster gates. Only
# the cylindrical-compatible rules carry a LevelStructure; on the others
# ``level_indices`` is ``[arange(N)]`` — which is the point: the accessor
# answers the fourth name on every rule, degenerately where the rule has
# one level.
_RULES = [
    ("gauss_legendre(4)", lambda: Quadrature.gauss_legendre(4)),
    ("level_symmetric(4)", lambda: Quadrature.level_symmetric(4)),
    ("product(4,8)", lambda: Quadrature.product(4, 8)),
    ("folded_product(4,8)", lambda: Quadrature.folded_product(4, 8)),
    ("lebedev(5)", lambda: Quadrature.lebedev(5)),
]


class TestG1GeneratorIsProvenanceNotIdentity:
    """G1 — the exclusion from ``Axis.content_parts``: a RULING (provenance is not identity).

    (The parts were a hand-written ``_identity_key`` tuple until #405 P1 step 5,
    2026-10-02, when ``Axis`` moved onto the one content encoder.)

    The negative legs (a DIFFERENT weights/kind is still unequal) already
    ship as ``test_axis.py::test_weights_are_part_of_identity`` and
    ``::test_basis_kind_participates_in_identity`` — cited, not
    duplicated; this class is the positive-exclusion half of that pair.
    """

    @pytest.mark.parametrize("name,factory", _RULES)
    def test_distinct_instance_content_equal_generators_are_ONE_axis(
        self, name, factory
    ):
        """G1a — provenance is not identity, and #403 never reaches the axis.

        Mutation (measured before step 5, on the then ``_identity_key``): put
        ``generator`` in the key → ``a1 == a2`` RAISES ``ValueError`` (ndarray
        truth ambiguity through ``DiscreteMeasure.__eq__``). Since step 5 the
        same mutation of ``content_parts`` makes the generator a contentless
        part, so ``a1 == a2`` falls back to identity (False) and ``hash``
        raises: still red here.
        """
        q1, q2 = factory(), factory()
        _require(q1 is not q2, "the two rules must be distinct instances")
        a1, a2 = Axis.of_quadrature(q1, "angular"), Axis.of_quadrature(q2, "angular")
        _require(a1 == a2, f"{name}: content-equal generators must not split identity")
        _require(hash(a1) == hash(a2), f"{name}: hash must agree with eq")
        _require(len({a1, a2}) == 1, f"{name}: one axis, one set member")

    def test_the_origin_arm_is_provenance_and_the_atoms_are_content(self):
        """G1b — re-posed at S7 (ruling 71; until then "a generator-ful axis
        equals the generator-LESS one", which INVERTED). The arm is
        provenance: a quadrature and its measure mint one axis. The atoms are
        content: the same weights declared with no atoms are another axis."""
        q = Quadrature.gauss_legendre(4)
        minted = Axis.of_quadrature(q, "angular")
        through_measure = Axis.of_measure(q.measure, "angular")
        _require(
            minted == through_measure and hash(minted) == hash(through_measure),
            "FromQuadrature(q) and FromMeasure(q.measure) must mint one axis",
        )
        _require(len({minted, through_measure}) == 1, "one axis, one set member")
        bare = Axis(
            "angular", (q.n_ordinates,), weights=np.asarray(q.weights, float),
            origin=Declared(BasisKind.NODAL),
        )
        _require(not (minted == bare), "an atomless axis must not pass for the sampled one")

    def test_the_exclusion_is_a_RULING_and_the_measure_leg_still_bites(self):
        """G1c ⭐ — WHY the generator stays out of the key, restated 2026-09-12.

        Until the consumers campaign's S1b a ``Quadrature`` was UNHASHABLE
        (``frozen=False`` dataclass with ``eq=True``), so an identity key
        containing the generator would have made ``hash(Axis)`` RAISE — the
        exclusion was a necessity. S1b gave ``Quadrature`` CONTENT identity
        (R-cc3: it is a generating datum of every SN problem), so that leg
        INVERTED: two rules of equal content now hash equal. The exclusion
        survives as a RULING — provenance is not identity (G1a/G1b above
        are its gates: an axis re-weighted through one generator equals
        the same weights minted through another) — and the second leg
        still bites: a ``DiscreteMeasure`` (frozen ``eq=True`` over ndarrays)
        raises the truth-value ambiguity on ``==``, so a key holding the
        MEASURE rather than its bytes would still make ``Axis.__eq__`` raise.
        """
        q1, q2 = Quadrature.gauss_legendre(4), Quadrature.gauss_legendre(4)
        _require(q1 is not q2 and q1 == q2 and hash(q1) == hash(q2),
                 "S1b: two rules of equal content are one quadrature (content identity)")
        with pytest.raises(ValueError, match="ambiguous"):
            _ = q1.measure == q2.measure


class TestG3MintFidelity:
    """G3 — ``Axis.of_quadrature(quad, label)`` ≡ the literal mint, plus provenance.

    Independence, stated honestly: the two sides are NOT independent on
    the weight VALUES (``q.weights`` is ``measure.weights``, which is
    what ``Axis.of_measure`` reads). What this gate genuinely pins
    is the THREADING — label, shape spelling (``(n_points,)`` vs
    ``(quad.n_ordinates,)``: same value, different producer chain), ``kind``, the
    ``weights=`` wiring, and that ``replace`` preserves the canonicalized
    bytes. The independent anchor for the weight VALUES is the
    quadrature weight-sum / moment-exactness suite
    (``tests/gates/numerics/test_quadrature_directional.py``).
    """

    @pytest.mark.parametrize("name,factory", _RULES)
    def test_quad_axis_is_the_literal_construction_plus_provenance(
        self, name, factory
    ):
        q = factory()
        minted = Axis.of_quadrature(q, "angular")
        literal = Axis(
            "angular", (q.n_ordinates,), weights=np.asarray(q.weights, float),
            origin=Declared(BasisKind.NODAL),
        )
        _require(minted.label == literal.label == "angular", f"{name}: label")
        _require(minted.shape == literal.shape == (q.n_ordinates,), f"{name}: shape")
        # the MODAL mutation's catcher:
        _require(minted.kind is BasisKind.NODAL, f"{name}: kind must be NODAL")
        _require(
            (minted.weights is None) == (literal.weights is None),
            f"{name}: canonicalization diverged between mint and literal",
        )
        if minted.weights is not None:
            assert literal.weights is not None
            _require(
                minted.weights.tobytes() == literal.weights.tobytes(),
                f"{name}: the weight BYTES feed the digest — array_equal "
                f"is not enough",
            )
            _require(
                not minted.weights.flags.writeable,
                f"{name}: the read-only canonicalization must survive replace()",
            )
        # the drop-the-upgrade mutation's catcher:
        _require(
            isinstance(minted.origin, FromQuadrature) and minted.origin.quadrature is q,
            f"{name}: the origin is the QUADRATURE, not the measure "
            f"(got {type(minted.origin).__name__})",
        )
        # the atoms and the support are the rule's (ruling 71)
        npt.assert_array_equal(minted.nodes, q.measure.nodes, err_msg=f"{name}: atoms")
        _require(minted.support == q.measure.support, f"{name}: support")


def _quadrature_through(space: FunctionSpace) -> Quadrature:
    """The quadrature recovered THROUGH a space: its angular axis's origin."""
    match space.axis("angular").origin:
        case FromQuadrature(quadrature):
            return quadrature
        case origin:
            pytest.fail(f"the angular axis was minted through {type(origin).__name__}")


class TestG4TheFourNamesAnswerThroughTheSpace:
    """G4 ⭐ — the done-when's core, verbatim: ``mu_x`` / ``eta`` / ``mu_z``
    / ``level_indices`` are reachable from a SPACE, with no quadrature in
    the caller's hand. A numerics-layer gate because the CONTRACT is
    numerics-layer — it builds a ``FunctionSpace``, not an ``SNProblem``.

    ⚠ ``eta`` rows are CONTRACT rows, not re-point coverage: [M] ``eta``
    has no solve-time consumer outside the mesh (verification plan R6).
    """

    @pytest.mark.parametrize("name,factory", _RULES)
    @pytest.mark.parametrize("attr", ["mu_x", "eta", "mu_z"])
    def test_direction_cosines_answer_through_the_space(self, name, factory, attr):
        q = factory()
        space = FunctionSpace.of_axes(Axis.of_quadrature(q, "angular"), _point())
        g = _quadrature_through(space)
        npt.assert_array_equal(
            getattr(g, attr), getattr(q, attr),
            err_msg=f"{name}.{attr} lost through the space",
        )

    @pytest.mark.parametrize("name,factory", _RULES)
    def test_level_indices_answers_through_the_space(self, name, factory):
        """⛔ The name a bare DiscreteMeasure CANNOT answer.

        ``level_indices`` lives on the LevelStructure side-channel, not
        in ``(nodes, weights, support, invariance_group, exactness)`` —
        which is exactly why the angular generator is the Quadrature and
        not the measure (``scratch/cs5_ground_measure.md`` §L). This row
        is that ruling's witness.
        """
        q = factory()
        g = _quadrature_through(FunctionSpace.of_axes(Axis.of_quadrature(q, "angular"), _point()))
        got, want = g.level_indices, q.level_indices
        _require(len(got) == len(want), f"{name}: level count")
        for k, (a, b) in enumerate(zip(got, want)):
            npt.assert_array_equal(a, b, err_msg=f"{name} level {k}")

    @pytest.mark.parametrize("name,factory", _RULES)
    def test_the_full_generator_surface_answers_through_the_space(
        self, name, factory
    ):
        """The done-when's four names are the HARD ones, not the whole
        reach-past set — [M] ``N`` (12 reads) and ``weights`` (12) reach
        past too (ground memo §I). The generator answers the FULL
        ``AngularMeasure`` surface, ``level_structure`` included, so the
        P4-remainder re-point never finds a name this contract lacks.
        """
        q = factory()
        g = _quadrature_through(FunctionSpace.of_axes(Axis.of_quadrature(q, "angular"), _point()))
        _require(g.n_ordinates == q.n_ordinates, f"{name}: N")
        npt.assert_array_equal(g.weights, q.weights, err_msg=f"{name}: weights")
        _require(
            g.level_structure is q.level_structure,
            f"{name}: level_structure must be the rule's own side-channel "
            f"(None on slab/2-D rules is the honest reading)",
        )

    def test_the_bare_measure_cannot_answer_the_fourth_name(self):
        """The refutation leg — do NOT weaken the done-when to three.

        A ``DiscreteMeasure``-generated axis answers three of four; the
        fourth is the whole reason ``Axis.of_quadrature`` exists. If a
        measure starts answering ``level_indices``, this design must be
        re-ruled — that is a signal, not a convenience.
        """
        q = Quadrature.level_symmetric(4)
        measure_axis = Axis.of_measure(q.measure, "angular")
        _require(
            isinstance(measure_axis.origin, FromMeasure)
            and not hasattr(measure_axis.origin.measure, "level_indices"),
            "a DiscreteMeasure must NOT answer level_indices — if it "
            "starts to, Axis.of_quadrature's raison d'etre changed and this "
            "design must be re-ruled",
        )


class TestG8TheMintIsASectionOfTheForgetfulMap:
    """G8 ⭐ — the intrinsic law: for every axis ``a`` minted from a measure or
    a quadrature, re-minting from its origin's object gives ``a`` back (the mint
    is a SECTION, so what the axis keeps is recoverable).

    ⚠ Scoped honestly: the law ranges over GENERATOR-FUL axes only. A
    generator-less axis (the homogeneous counting point, every
    EnergyAxis, the MODAL moment axis) is outside the domain, not a
    counterexample. The spatial round-trip rows live beside G6 in
    ``test_space_of_axes.py`` — the rank-d arm is generator-less by
    contract (the CS2 seam), so only the rank-1 arm is in the law's
    domain.

    Mutation: mint at a shape other than ``(generator.n_points,)``, or
    hand-pass a ``generator=`` that did not produce the axis — this is
    the gate that catches the rank-d mismatch (verification plan R1).
    """

    @pytest.mark.parametrize("name,factory", _RULES)
    def test_angular_round_trip(self, name, factory):
        a = Axis.of_quadrature(factory(), "angular")
        g = _quadrature_through(FunctionSpace.of_axes(a))
        _require(Axis.of_quadrature(g, a.label) == a, f"{name}: not a section")
        again = Axis.of_quadrature(g, a.label).origin
        _require(
            isinstance(again, FromQuadrature) and again.quadrature is g,
            f"{name}: the round-trip must preserve the quadrature itself",
        )

    def test_measure_round_trip(self):
        """The NODAL-measure leg of the same law (spatial-shaped)."""
        from orpheus.numerics.measure import DiscreteMeasure

        m = DiscreteMeasure(
            nodes=np.array([0.5, 1.5, 2.5]),
            weights=np.array([0.2, 0.3, 0.5]),
            support=RealSpace(1),
        )
        a = Axis.of_measure(m, "spatial")
        _require(isinstance(a.origin, FromMeasure) and a.origin.measure is m, "the mint must record its measure")
        _require(Axis.of_measure(m, a.label) == a, "measure mint: not a section")


class TestG5TheQuadratureRefusalIsOneSpelling:
    """P4-remainder G5, re-posed at S7: each consumer that needs the generating
    quadrature matches ``case FromQuadrature(quadrature)`` itself, and its other
    arm raises :func:`~orpheus.numerics.axis.quadrature_refusal`, the one
    spelling of the message, naming the axis, the asker and the origin it got.
    Per-consumer refusal rows live beside their consumers
    (``tests/gates/sn/mesh/test_reduced_operator.py``,
    ``tests/gates/sn/sweep/curvilinear/test_angular_closure.py``); the census of
    the ``case FromQuadrature`` sites is ``test_axis_identity.py`` O3 (it
    replaced G5_6b's ``generator_as`` call-site census).
    """

    def test_a_declared_axis_is_refused_with_both_names(self):
        q = Quadrature.gauss_legendre(4)
        bare = Axis(
            "angular", (q.n_ordinates,), weights=np.asarray(q.weights, float),
            origin=Declared(BasisKind.NODAL),
        )
        message = str(quadrature_refusal(bare, "somebody"))
        for fragment in ("angular", "somebody", "minted through", "Declared"):
            _require(fragment in message, f"the refusal must carry {fragment!r}: {message}")

    def test_a_measure_minted_axis_is_refused_naming_its_origin(self):
        """A measure-minted axis is not a quadrature-minted one — the
        consumers' match is a type claim, not a None check (parse, don't
        validate). The message names what it got."""
        from orpheus.numerics.measure import DiscreteMeasure

        m = DiscreteMeasure(
            nodes=np.array([0.5, 1.5]), weights=np.array([0.4, 0.6]),
            support=RealSpace(1),
        )
        _require("FromMeasure" in str(quadrature_refusal(Axis.of_measure(m, "spatial"), "test")),
                 "the refusal must name the origin it got")

    def test_fragment_disjointness_with_the_space_refusal(self):
        """G5.4 (L43c) — the neighbouring refusal on the same consumer
        path spells 'axis lookup:'; this message must not, so a pin on
        either can never match the other's raise."""
        q = Quadrature.gauss_legendre(4)
        bare = Axis(
            "angular", (q.n_ordinates,), weights=np.asarray(q.weights, float),
            origin=Declared(BasisKind.NODAL),
        )
        _require("axis lookup:" not in str(quadrature_refusal(bare, "x")), "fragments must be disjoint")

    def test_G5_6a_the_shipped_atomless_axes_are_Declared(self):
        """The inventory row: the counting point and every EnergyAxis are
        ``Declared`` (ruling 82), so a quadrature consumer handed one refuses."""
        from orpheus.numerics.axis import EnergyAxis

        for ax in (Axis("spatial", (1,), origin=Declared(BasisKind.NODAL)), EnergyAxis.synthetic(3)):
            _require(isinstance(ax.origin, Declared), f"{ax.label}: not Declared")
            _require(ax.support is None and ax.nodes is None, f"{ax.label}: an atomless axis has no atoms")
