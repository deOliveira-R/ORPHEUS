r"""The admission mechanism: the lattice, the registration, the one door, the functor (wave 2, slice S2a, gate 13).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 2.4 ("``admission.py``", "``realization.py``", decisions 48,
52, 53, 56), law L-D (admission equals realization) and law L-L (the realizer is a homomorphism of the algebra), 6.2
S2a gate 13 (a) to (i).

The mechanism is gated on a TEST-LOCAL table (``_admission_helpers.make_table``), never on a method's: the real
tables are slices S2b, S3 and S4, each appending its rows to ``test_admission_matrix.py``. The claims, kind and
layer per row (``vv-principles``):

* (a) ``test_plumbing_needs_an_issue``, ``test_the_lattice_of_admissions``: THEOREM on the lattice
  ``Realized(True) > Realized(False) > Approximated > Refused``, defined ONCE: ``Admission.meet`` is the lesser, the
  ``Approximated`` names and issues union, the first ``Refused`` wins, the empty meet is ``Realized(True)``; and a
  refusal on the plumbing axis owes the issue that tracks the missing machinery.
* (b) ``test_registration_by_class_keyword``: STRUCTURE. A realizer registers at class creation (a class keyword), a
  duplicate, a key outside the closed set, a sealed table and a second ``DeckRealizer`` are refused, an unsealed
  table refuses to answer, the state cone is a REQUIRED keyword.
* (c) ``test_one_door`` (law L-D): THEOREM over a grid. ``admit`` is ``Refused`` iff ``realize`` raises
  ``BoundaryRefused``, which carries the ``Refused``; ``require`` is ``admit`` plus that raise; no realizer body
  refuses on its own. The checker has a planted positive control (a realizer that raises on its own).
* (d) ``test_exact_type_lookup``: a subclass is a new declaration and owes its own entry.
* (e) ``test_the_unregistered_door_names_the_axis``: a condition with ``return_fraction < 0`` on a non-negative-flux
  table is ``Refused(STATE_CONE)`` without the table importing the condition's type; otherwise ``Refused(PLUMBING)``.
* (f) ``test_catalogue_is_an_upper_bound``: no sampled ``admit`` exceeds its realizer's ``status`` in the lattice.
* (g) ``test_the_deck_arm``: the deck realizer receives the pairing the geometry derives, and its exhaustive match
  ends in ``assert_never``.
* (h) ``test_the_functor_law`` (law L-L): the image of the empty network is the ``ZeroNetwork`` realizer's, the image
  of one chain is that chain's realizer's with the weight beside it, the image of a sum is the sum of the images
  (in ``CHAIN_ORDER``), the zero network is consulted only for an empty network, the source only for a nonzero one.
* (i) ``test_face_site``: the three geometric site classes are derived once, from ``(face, coord, dimension,
  hollow_inner)``.
* (j) ``test_the_source_verb``: ``sample_inflow`` is ``require`` plus the registered ``SourceRealizer.sample``,
  whose default is ``sample_terms``.

Activation: the SN-like table activates every refusal axis the section lists and the three site classes, so a
mutation of any arm of the meet, the door or the lookup moves some row; a one-unit table would null the meet.

First reds, at the base tip ``4ec95bc6``: ``ModuleNotFoundError`` on every row (the first read of the package). The
inputs that redden each row once the package exists are the mutations below, each applied to a COPY of the module:

* (a) accept ``Refused(PLUMBING, r)`` without an issue; make ``meet`` return the greater; make the empty meet
  ``Realized(False)``; let the meet of two ``Approximated`` keep one name.
* (b) a no-op ``seal``; skip the duplicate check; accept any class as a key; let ``RealizerTable("m")`` default the
  state cone; allow two deck realizers.
* (c) a ``realize`` that skips ``admit``; a ``require`` that does not raise; a ``BoundaryRefused`` without the
  admission.
* (d) look the realizer up through the MRO (``isinstance`` scan).
* (e) swap the two refusals; consult the table's cone for a wall; omit the issue of the plumbing refusal.
* (f) a realizer whose ``admit`` returns ``Realized(True)`` under an ``Approximated`` ``status``; ``catalogue`` that
  returns the realizers' sampled admissions rather than their statuses.
* (g) a ``DeckRealizer.admit`` without the ``assert_never`` arm; a table passing the face instead of the pairing.
* (h) consult the ``ZeroNetwork`` realizer for every wall; sum the images in insertion order; drop the weight when
  calling a chain's realizer; consult the source realizer for a zero source.
* (i) classify the hollow inner face as the outer; treat a Cartesian axis of a radial mesh as curved.
* (j) a ``sample_inflow`` that calls the sampler before ``require``.
"""

from __future__ import annotations

import dataclasses
import itertools

import numpy as np
import pytest

from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.deck import MirrorDeck, TranslationDeck, pairing_of
from tests.gates._content_identity_helpers import require
from tests.gates.transport._admission_helpers import SITES, SPACE, Img, make_table, site, values
from tests.gates.transport._boundary_conditions_helpers import bc, function_space, local, node

pytestmark = pytest.mark.foundation

_HERE = "test_admission_mechanism.py"
_ALGEBRA = node("test_wall_condition_algebra.py", "test_the_algebra_and_the_cone")
_LATTICE = node(_HERE, "test_the_lattice_of_admissions")
_REGISTRATION = node(_HERE, "test_registration_by_class_keyword")
_DOOR = node(_HERE, "test_one_door")
_PAIRING = "tests/gates/geometry/test_deck.py::test_pairing_of_and_its_adjoint"

#: What a refused registration may raise (the spec says "raises"; the sealed table's exception is named).
_REFUSED_REGISTRATION = (TypeError, ValueError, KeyError)


def _rank(admission) -> int:
    if isinstance(admission, bc.Realized):
        return 3 if admission.faithful else 2
    if isinstance(admission, bc.Approximated):
        return 1
    require(isinstance(admission, bc.Refused), f"{admission!r} is not an admission")
    return 0


# ── (a) the lattice ──────────────────────────────────────────────────────────


def test_plumbing_needs_an_issue() -> None:
    """``Refused(PLUMBING, "r")`` raises; with an issue it constructs; the three other axes construct without one;
    ``Approximated`` holds frozensets and ``Approximated.of(name, issue)`` is the one-member spelling;
    ``Realized(True)`` and ``Realized(False)`` differ; the members are frozen; the axes are exactly the four of the
    section.

    Reddens on: a ``Refused`` accepting PLUMBING with no issue (a refusal with no theory behind it and no tracker).
    """
    Axis = bc.RefusalAxis
    require({a.name for a in Axis} == {"ANGULAR_RESOLUTION", "TOPOLOGY", "STATE_CONE", "PLUMBING"}, f"the axes are {[a.name for a in Axis]}")
    with pytest.raises((ValueError, TypeError)):
        bc.Refused(Axis.PLUMBING, "r")
    require(bc.Refused(Axis.PLUMBING, "r", 511).issue == 511, "a plumbing refusal keeps its issue")
    for axis in (Axis.ANGULAR_RESOLUTION, Axis.TOPOLOGY, Axis.STATE_CONE):
        require(bc.Refused(axis, "r").issue is None, f"{axis.name} owes no issue")
    approx = bc.Approximated(frozenset({"a square cell"}), frozenset({547}))
    require(isinstance(approx.names, frozenset) and isinstance(approx.issues, frozenset), "Approximated holds frozensets")
    require(bc.Approximated.of("a square cell", 547) == approx, "Approximated.of is the one-member spelling")
    require(bc.Realized(True) != bc.Realized(False) and bc.Realized(True) == bc.Realized(True), "tier 1 and tier 2 differ")
    for member in (approx, bc.Realized(True), bc.Refused(Axis.TOPOLOGY, "r")):
        with pytest.raises(dataclasses.FrozenInstanceError):
            member.reason = "x"  # type: ignore[misc, attr-defined]  # the mutation input


def test_the_lattice_of_admissions() -> None:
    """``Realized(True) > Realized(False) > Approximated > Refused``; ``Admission.meet(*members)`` is the lesser; the
    meet of two ``Approximated`` is ``Approximated(x.names | y.names, x.issues | y.issues)``; the meet of two
    ``Refused`` is the FIRST argument's; the meet of no members is ``Realized(True)``. Over the seven-member grid:
    idempotent, ``Realized(True)`` is the identity, associative on all 343 triples, commutative except between two
    distinct refusals (the documented deterministic tie-break).

    Reddens on: a ``meet`` returning the greater; the empty meet ``Realized(False)``; an ``Approximated`` meet
    keeping one name; a ``Refused`` tie-break by the second argument.
    """
    meet = bc.Admission.meet
    Axis = bc.RefusalAxis
    r1, r0 = bc.Realized(True), bc.Realized(False)
    ax, ay = bc.Approximated.of("x", 547), bc.Approximated.of("y", 548)
    f1, f2, f3 = bc.Refused(Axis.STATE_CONE, "a"), bc.Refused(Axis.TOPOLOGY, "b"), bc.Refused(Axis.PLUMBING, "c", 511)
    grid = [r1, r0, ax, ay, f1, f2, f3]

    require(meet() == r1, "the meet of no members is Realized(True)")
    union = bc.Approximated(frozenset({"x", "y"}), frozenset({547, 548}))
    cases = [
        ((r1, r0), r0), ((r0, r1), r0), ((r1, ax), ax), ((r0, ax), ax), ((ax, ay), union), ((ay, ax), union),
        ((ax, f1), f1), ((f1, ax), f1), ((r1, f2), f2), ((f1, f2), f1), ((f2, f1), f2), ((f3, f1), f3),
        ((r1, ax, ay, r0), union), ((r1, ax, f2, ay), f2), ((ax, f1, f2), f1),
    ]
    for members, expected in cases:
        require(meet(*members) == expected, f"meet{tuple(type(m).__name__ for m in members)} = {meet(*members)!r}, expected {expected!r}")
    for member in grid:
        require(meet(member) == member and meet(member, member) == member, f"{member!r} is not idempotent")
        require(meet(r1, member) == member and meet(member, r1) == member, f"Realized(True) is not the identity at {member!r}")
    for a, b, c in itertools.product(grid, repeat=3):
        require(meet(meet(a, b), c) == meet(a, meet(b, c)), f"meet is not associative at {a!r}, {b!r}, {c!r}")
    for a, b in itertools.product(grid, repeat=2):
        both_refused = isinstance(a, bc.Refused) and isinstance(b, bc.Refused) and a != b
        if not both_refused:
            require(meet(a, b) == meet(b, a), f"meet is not commutative at {a!r}, {b!r}")
    for a, b in itertools.product(grid, repeat=2):
        require(_rank(meet(a, b)) <= min(_rank(a), _rank(b)), f"the meet of {a!r}, {b!r} is above one of them")


# ── (b) registration ─────────────────────────────────────────────────────────


@pytest.mark.rests_on(_LATTICE)
def test_registration_by_class_keyword() -> None:
    """A realizer registers at class creation. A duplicate key, a key outside ``get_args(Chain)``, ``ZeroNetwork``
    and the ``BoundaryCondition`` subclasses (an ``int``, a factor, ``InflowTerms`` for a response realizer) and a
    second ``SourceRealizer`` or ``DeckRealizer`` are refused; a method-only condition type is a key; registration
    after ``seal()`` raises ``SealedTableError``; an unsealed table refuses ``admit``, ``realize`` and ``catalogue``
    with ``TableNotSealedError``; ``RealizerTable("m")`` without ``state_cone`` is a ``TypeError``; two tables do not
    share a registry; the source realizer registers under ``InflowTerms``.

    Reddens on: the mutations of the module docstring, (b).
    """
    R = bc.Realized
    built = make_table(seal=False)
    t = built.table

    with pytest.raises(bc.TableNotSealedError):
        t.admit(bc.WallCondition.vacuum(), site("CART"))
    with pytest.raises(bc.TableNotSealedError):
        t.realize(bc.WallCondition.vacuum(), site("CART"), SPACE)
    with pytest.raises(bc.TableNotSealedError):
        t.catalogue()

    with pytest.raises(_REFUSED_REGISTRATION):
        class DuplicateSpec(bc.ResponseRealizer, table=t, realizes=bc.SpecularChain):
            status = R(True)

            def realize(self, unit, weight, site_, space):
                return Img(())

    for bad_key in (int, bc.MirrorPermutation, bc.InflowTerms, bc.IsotropicSection):
        with pytest.raises(_REFUSED_REGISTRATION):
            class BadKey(bc.ResponseRealizer, table=t, realizes=bad_key):  # type: ignore[arg-type, misc]
                status = R(True)

                def realize(self, unit, weight, site_, space):
                    return Img(())

    with pytest.raises(_REFUSED_REGISTRATION):
        class SecondSource(bc.SourceRealizer, table=t):
            status = R(True)

    with pytest.raises(_REFUSED_REGISTRATION):
        class SecondDeck(bc.DeckRealizer, table=t):
            def admit(self, deck, pairing, site_, space=None):
                return R(True)

            def realize(self, deck, pairing, site_, space):
                return Img(())

    method_only = local().Probe
    class ProbeRealizer(bc.ResponseRealizer, table=t, realizes=method_only):  # a method-only condition type is a key
        status = R(True)

        def realize(self, unit, weight, site_, space):
            return Img((("probe", weight),))

    t.seal()
    cat = t.catalogue()
    require(set(cat) == {bc.SpecularChain, bc.DiffuseChain, bc.ZeroNetwork, bc.InflowTerms, method_only}, f"the keys are {sorted(k.__name__ for k in cat)}")

    with pytest.raises(bc.SealedTableError):
        class Late(bc.ResponseRealizer, table=t, realizes=local().NegativeProbe):
            status = R(True)

            def realize(self, unit, weight, site_, space):
                return Img(())
    with pytest.raises(bc.SealedTableError):
        class LateSource(bc.SourceRealizer, table=t):
            status = R(True)
    require(local().NegativeProbe not in t.catalogue(), "a refused late registration still entered the table")

    with pytest.raises(TypeError):
        bc.RealizerTable("m")  # type: ignore[call-arg]  # the cone is a required keyword
    other = make_table(method="other")
    require(set(other.table.catalogue()) == {bc.SpecularChain, bc.DiffuseChain, bc.ZeroNetwork, bc.InflowTerms}, "a second table shares the first's registry")
    require(method_only not in other.table.catalogue(), "the first table's late key leaked into the second")


# ── (c) the one door ─────────────────────────────────────────────────────────


def _door_disagreements(table, grid) -> list[str]:
    """Every (value, site) where ``admit`` is Refused iff ``realize`` and ``require`` raise ``BoundaryRefused`` fails,
    or the raised error does not carry the very admission ``admit`` returned."""
    out: list[str] = []
    for label, value, where in grid:
        admission = table.admit(value, where)
        refused = isinstance(admission, bc.Refused)
        for verb, call in (("realize", lambda: table.realize(value, where, SPACE)), ("require", lambda: table.require(value, where))):
            try:
                call()
                raised, error = False, None
            except bc.BoundaryRefused as err:
                raised, error = True, err
            if raised != refused:
                out.append(f"{label} @ {where.site_class.name}: admit refused={refused} but {verb} raised={raised}")
            if error is not None and error.admission != admission:
                out.append(f"{label} @ {where.site_class.name}: {verb}'s error carries {error.admission!r}, admit returned {admission!r}")
    return out


@pytest.mark.rests_on(_REGISTRATION, _ALGEBRA)
def test_one_door() -> None:
    """Law L-D over the grid (seven value classes and two method-only conditions, times the three sites): ``admit``
    is ``Refused`` iff ``realize`` and ``require`` raise ``BoundaryRefused``, whose ``admission`` is the very value
    ``admit`` returned; the message is ``f"{method} realizes {sorted(admitted names)}; not {type(declared).__name__}
    ({axis}: {reason})"`` (its fragments); ``method`` and ``face`` are carried; ``require`` returns ``None`` when it
    does not raise; ``space`` reaches the realizer's ``admit`` when given and is ``None`` otherwise.

    The positive control (X1): a planted table whose realizer REFUSES ON ITS OWN in ``realize`` while ``admit`` says
    ``Realized`` is reported by the same checker; the honest table is clean. Population: the grid holds refused and
    admitted cells of both verbs (a gate over a grid with one outcome cannot fail).

    Reddens on: a ``realize`` that skips ``admit`` (the refused cells then build an image); a ``require`` that never
    raises; a ``BoundaryRefused`` built without the admission.
    """
    built = make_table()
    t = built.table
    probe, negative = local().Probe(), local().NegativeProbe()
    grid = [(label, value, site(kind)) for label, value in values().items() for kind in SITES]
    grid += [(f"{name} (no entry)", cond, site(kind)) for name, cond in (("probe", probe), ("negative", negative)) for kind in SITES]
    require(len(grid) == 27, f"the grid has {len(grid)} cells")
    outcomes = {isinstance(t.admit(v, s), bc.Refused) for _, v, s in grid}
    require(outcomes == {True, False}, "the grid has a single outcome: the iff is untested")
    require(_door_disagreements(t, grid) == [], f"admit and realize disagree: {_door_disagreements(t, grid)}")

    # the planted positive control: a realizer that refuses on its own
    planted = bc.RealizerTable("planted", state_cone=bc.StateCone.NONNEGATIVE_FLUX)

    class OwnRefusal(bc.BoundaryRefused):
        """Built without the base constructor, whose signature the section does not fix."""

        def __init__(self, admission):
            Exception.__init__(self, "refused in the body")
            self.declaration, self.method, self.face, self.admission = None, "planted", "xmax", admission

    class SelfRefusing(bc.ResponseRealizer, table=planted, realizes=bc.SpecularChain):
        status = bc.Realized(True)

        def realize(self, unit, weight, site_, space):
            raise OwnRefusal(bc.Refused(bc.RefusalAxis.TOPOLOGY, "r"))

    planted.seal()
    caught = _door_disagreements(planted, [("specular", bc.WallCondition.specular(0.5), site("CART"))])
    require(caught != [], "control: a realizer that refuses in its own body is not reported")

    where = site("RAD_OUTER")
    white = bc.WallCondition.white(0.5)
    with pytest.raises(bc.BoundaryRefused) as caught_error:
        t.realize(white, where, SPACE)
    err = caught_error.value
    require(isinstance(err, bc.BoundaryError) and err.method == "mech" and err.face == "xmax", f"the error carries method {err.method!r}, face {err.face!r}")
    require(isinstance(err.admission, bc.Refused) and err.admission.axis is bc.RefusalAxis.PLUMBING and err.admission.issue == 653, f"the admission is {err.admission!r}")
    text = str(err)
    for fragment in ("mech realizes", "not WallCondition", "PLUMBING", err.admission.reason, "SpecularChain", "DiffuseChain", "ZeroNetwork", "InflowTerms"):
        require(fragment in text, f"the message lacks {fragment!r}: {text}")
    require(t.require(bc.WallCondition.vacuum(), site("CART")) is None, "require returns None for an admitted value")

    built.spies.admit.clear()
    t.admit(bc.WallCondition.specular(0.5), site("CART"))
    t.admit(bc.WallCondition.specular(0.5), site("CART"), SPACE)
    require([entry[4] for entry in built.spies.admit] == [None, SPACE], f"the space reached the realizer's admit as {[entry[4] for entry in built.spies.admit]}")


# ── (d) exact-type lookup ────────────────────────────────────────────────────


@pytest.mark.rests_on(_REGISTRATION)
def test_exact_type_lookup() -> None:
    """A subclass of a method-only condition is a NEW declaration: with an entry for the parent only, the parent is
    realized and the child is refused (a method-only condition of another method must not be realized silently as
    its parent).

    Reddens on: a lookup through the MRO (``isinstance`` scan): the child is then admitted as its parent.
    """
    parent, child = local().Probe, local().ProbeChild
    table = bc.RealizerTable("exact", state_cone=bc.StateCone.NONNEGATIVE_FLUX)

    class ProbeR(bc.ResponseRealizer, table=table, realizes=parent):
        status = bc.Realized(True)

        def realize(self, unit, weight, site_, space):
            return Img((("probe", weight),))

    table.seal()
    where = site("CART")
    require(table.admit(parent(), where) == bc.Realized(True), "the registered type is not admitted")
    adm = table.admit(child(), where)
    require(isinstance(adm, bc.Refused), f"the subclass was admitted as its parent: {adm!r}")
    require(table.realize(parent(), where, SPACE) == Img((("probe", 1.0),)), "the parent's image is the realizer's, with weight 1.0")
    with pytest.raises(bc.BoundaryRefused):
        table.realize(child(), where, SPACE)


# ── (e) the unregistered door ────────────────────────────────────────────────


@pytest.mark.rests_on(_REGISTRATION)
def test_the_unregistered_door_names_the_axis() -> None:
    """A condition with ``return_fraction < 0`` and no entry is ``Refused(STATE_CONE)`` on a ``NONNEGATIVE_FLUX``
    table and ``Refused(PLUMBING)`` on a ``SIGNED_CURRENT`` one; a condition with ``return_fraction >= 0`` and no
    entry is ``Refused(PLUMBING)`` (carrying an issue, since the axis owes one), and the message names the admitted
    keys; a registered ``NegativeProbe`` is admitted on a signed table; a wall whose chain has no entry
    (``white`` on a specular-only table; ``vacuum`` on a table without the zero network) is ``Refused(PLUMBING)``.

    Reddens on: swapping the two refusals; reading the cone for a wall; a plumbing refusal without an issue (which
    raises at construction, so the door itself would raise).
    """
    where = site("CART")
    negative, probe = local().NegativeProbe(), local().Probe()
    nonneg = make_table().table
    signed = bc.RealizerTable("signed", state_cone=bc.StateCone.SIGNED_CURRENT)
    signed.seal()

    adm = nonneg.admit(negative, where)
    require(isinstance(adm, bc.Refused) and adm.axis is bc.RefusalAxis.STATE_CONE, f"a negative return on a non-negative table is {adm!r}")
    adm = signed.admit(negative, where)
    require(isinstance(adm, bc.Refused) and adm.axis is bc.RefusalAxis.PLUMBING and adm.issue is not None, f"a negative return on a signed table with no entry is {adm!r}")
    for table in (nonneg, signed):
        adm = table.admit(probe, where)
        require(isinstance(adm, bc.Refused) and adm.axis is bc.RefusalAxis.PLUMBING and adm.issue is not None, f"a non-negative return with no entry is {adm!r}")
    with pytest.raises(bc.BoundaryRefused) as caught:
        nonneg.realize(probe, where, SPACE)
    require("Probe" in str(caught.value) and "SpecularChain" in str(caught.value), f"the message does not name the type and the admitted keys: {caught.value}")

    registered = bc.RealizerTable("registered", state_cone=bc.StateCone.SIGNED_CURRENT)

    class NegativeR(bc.ResponseRealizer, table=registered, realizes=type(negative)):
        status = bc.Realized(True)

        def realize(self, unit, weight, site_, space):
            return Img((("negative", weight),))

    registered.seal()
    require(registered.admit(negative, where) == bc.Realized(True), "a registered signed condition is admitted on a signed table")

    spec_only = bc.RealizerTable("spec only", state_cone=bc.StateCone.NONNEGATIVE_FLUX)

    class OnlySpec(bc.ResponseRealizer, table=spec_only, realizes=bc.SpecularChain):
        status = bc.Realized(True)

        def realize(self, unit, weight, site_, space):
            return Img((("spec", weight),))

    spec_only.seal()
    require(spec_only.admit(bc.WallCondition.specular(0.5), where) == bc.Realized(True), "the registered chain is admitted")
    for label, wall in (("white", bc.WallCondition.white(0.5)), ("vacuum", bc.WallCondition.vacuum()), ("inflow", bc.WallCondition.inflow(1.0))):
        adm = spec_only.admit(wall, where)
        require(isinstance(adm, bc.Refused) and adm.axis is bc.RefusalAxis.PLUMBING, f"{label} on a specular-only table is {adm!r}")


# ── (f) the catalogue ────────────────────────────────────────────────────────


def _catalogue_violations(table, grid) -> list[str]:
    """Every sampled ``admit`` that exceeds its realizer's catalogue entry in the lattice (``faithful`` included)."""
    catalogue = table.catalogue()
    out = []
    unit_of = {bc.WallCondition.specular(0.5): bc.SpecularChain, bc.WallCondition.white(0.5): bc.DiffuseChain}
    for _, value, where in grid:
        key = unit_of.get(value)
        if key is None or key not in catalogue:
            continue
        admitted, bound = table.admit(value, where), catalogue[key]
        if _rank(admitted) > _rank(bound):
            out.append(f"{key.__name__} @ {where.site_class.name}: admit {admitted!r} exceeds status {bound!r}")
    return out


@pytest.mark.rests_on(_DOOR)
def test_catalogue_is_an_upper_bound() -> None:
    """``catalogue()`` maps each registered key to its realizer's ``status``; every sampled ``admit`` is at most that
    status in the lattice (and ``faithful`` is never above the catalogue's). The honest table is clean; a planted
    realizer whose ``admit`` returns ``Realized(True)`` under an ``Approximated`` status is reported (the control).

    Reddens on: a realizer whose admit exceeds its status; a ``catalogue`` returning something other than the status.
    """
    built = make_table()
    cat = built.table.catalogue()
    require(set(cat) == {bc.SpecularChain, bc.DiffuseChain, bc.ZeroNetwork, bc.InflowTerms}, f"the keys are {sorted(k.__name__ for k in cat)}")
    require(all(cat[k] == bc.Realized(True) for k in cat), f"the statuses are {cat}")
    grid = [(label, v, site(kind)) for label, v in values().items() for kind in SITES]
    require(_catalogue_violations(built.table, grid) == [], f"the catalogue is not an upper bound: {_catalogue_violations(built.table, grid)}")

    lying = bc.RealizerTable("lying", state_cone=bc.StateCone.NONNEGATIVE_FLUX)

    class Lying(bc.ResponseRealizer, table=lying, realizes=bc.SpecularChain):
        status = bc.Approximated.of("a named approximation", 547)

        def admit(self, unit, weight, site_, space=None):
            return bc.Realized(True)

        def realize(self, unit, weight, site_, space):
            return Img(())

    lying.seal()
    require(lying.catalogue()[bc.SpecularChain] == bc.Approximated.of("a named approximation", 547), "the catalogue is not the status")
    require(_catalogue_violations(lying, grid) != [], "control: an admit above its status is not reported")


# ── (g) the deck arm ─────────────────────────────────────────────────────────


@pytest.mark.rests_on(_DOOR, _PAIRING)
def test_the_deck_arm() -> None:
    """The deck realizer receives ``pairing_of(deck, site.face, dimension=site.dimension)`` (the geometry's, not a
    re-derivation); ``admit`` is its match (the mirror admitted, the translation refused on the topology axis);
    ``realize`` goes through the same door; a foreign object reaches ``assert_never`` (``AssertionError``).

    Reddens on: a deck realizer without the ``assert_never`` arm (the foreign row returns instead of raising); a
    table passing the face where the pairing belongs.
    """
    built = make_table()
    t, where = built.table, site("CART")
    mirror, translation = MirrorDeck(), TranslationDeck()
    require(t.admit(mirror, where) == bc.Realized(True), "the mirror is not admitted")
    refused = t.admit(translation, where)
    require(isinstance(refused, bc.Refused) and refused.axis is bc.RefusalAxis.TOPOLOGY, f"the translation is {refused!r}")
    admit_calls = [entry for entry in built.spies.deck if entry[0] == "admit"]
    require(len(admit_calls) == 2, f"the deck realizer was called {len(admit_calls)} times")
    for _, deck, pairing, seen_site, space in admit_calls:
        require(pairing == pairing_of(deck, where.face, dimension=where.dimension), f"the pairing passed for {type(deck).__name__} is {pairing!r}")
        require(seen_site == where and space is None, "the site or the space was not passed through")
    require(t.realize(mirror, where, SPACE) == Img((("deck:MirrorDeck", 1.0),)), "the mirror's image is the deck realizer's")
    with pytest.raises(bc.BoundaryRefused):
        t.realize(translation, where, SPACE)
    realized = [entry for entry in built.spies.deck if entry[0] == "realize"]
    require(len(realized) == 1 and realized[0][4] is SPACE, "realize reached the deck realizer a number of times other than once (the refused deck must not)")

    foreign = built.DeckR()
    with pytest.raises(AssertionError):
        foreign.admit(object(), pairing_of(mirror, where.face, dimension=1), where)  # type: ignore[arg-type]  # the exhaustiveness input


# ── (h) the functor ──────────────────────────────────────────────────────────


@pytest.mark.rests_on(_DOOR, _ALGEBRA)
def test_the_functor_law() -> None:
    """Law L-L on the test table. The image of the empty network is the ``ZeroNetwork`` realizer's (a pure inflow's
    linear image too); of one chain, that chain's realizer's called with the CHAIN INSTANCE and its weight; of a
    sum, the sum of the images, in ``CHAIN_ORDER``; ``real(c * x)`` has the weights of ``c`` times ``x``'s;
    ``real(a + b) == real(a) + real(b)`` for dyadic coordinates; the source is NOT part of the image. In ``admit``:
    the zero network is consulted iff the network is empty, whatever the source; the source realizer iff the source
    is nonzero; the answer is the meet over the consulted units, and the weight reaches each chain's ``admit``.

    Reddens on: consulting the zero network for every wall (the variant table whose zero row refuses); summing the
    images in insertion order; dropping the weight; consulting the source for a zero source (the variant whose
    source row refuses).
    """
    where, W = site("CART"), bc.WallCondition
    built = make_table()
    t, spies = built.table, built.spies

    require(t.realize(W.vacuum(), where, SPACE) == Img((("zero", 1.0),)), "the empty network's image is not the zero realizer's")
    require(t.realize(W.inflow(2.0), where, SPACE) == Img((("zero", 1.0),)), "a pure inflow's linear image is not the zero map")
    require(t.realize(W.specular(0.5) + W.inflow(2.0), where, SPACE) == Img((("spec", 0.5),)), "the source entered the image")
    require(t.realize(W.specular(0.25), where, SPACE) == Img((("spec", 0.25),)), "one chain's image")
    mixed = W.from_amplitudes(0.25, 0.5)
    both = t.realize(mixed, where, SPACE)
    require(both == Img((("spec", 0.25), ("diff", 0.5))), f"the two-chain image is {both!r}, not the spec-then-diffuse sum")
    require(both.terms[0][0] == "spec", "the images are summed in CHAIN_ORDER")
    last = [entry for entry in spies.realize if entry[0] == "spec"][-1]
    require(last[1] == bc.SpecularChain() and last[2] == 0.25 and last[3] == where and last[4] is SPACE, f"the chain realizer was called with {last[1:]}")

    for c in (0.5, 0.25, 2.0):
        x = W.from_amplitudes(0.125, 0.25)
        scaled = t.realize(c * x, where, SPACE)
        require(scaled == t.realize(x, where, SPACE).scaled(c), f"real({c} x) != {c} real(x): {scaled!r}")
    a, b = W.specular(0.25), W.from_amplitudes(0.125, 0.5)
    require(t.realize(a + b, where, SPACE) == t.realize(a, where, SPACE) + t.realize(b, where, SPACE), "real(a + b) != real(a) + real(b)")

    # admit: the meet over the consulted units
    spies.admit.clear()
    t.admit(W.specular(0.5), where)
    require(spies.consulted("zero") == 0 and spies.consulted("source") == 0 and spies.consulted("spec") == 1, f"a specular wall consulted {[e[0] for e in spies.admit]}")
    require(spies.admit[0][2] == 0.5, "the weight did not reach the chain's admit")
    spies.admit.clear()
    t.admit(W.inflow(1.0), where)
    require([e[0] for e in spies.admit] == ["zero", "source"], f"a pure inflow consulted {[e[0] for e in spies.admit]}")
    spies.admit.clear()
    t.admit(W.vacuum(), where)
    require([e[0] for e in spies.admit] == ["zero"], f"vacuum consulted {[e[0] for e in spies.admit]}")
    spies.admit.clear()
    t.admit(mixed + W.inflow(1.0), where)
    require([e[0] for e in spies.admit] == ["spec", "diff", "source"], f"a mixed wall over a source consulted {[e[0] for e in spies.admit]}")

    # the discriminating variants: a refusing zero row and a refusing source row
    refuse = lambda site_: bc.Refused(bc.RefusalAxis.PLUMBING, "x", 900)  # noqa: E731
    zero_refuses = make_table(zero_admit=refuse).table
    require(zero_refuses.admit(W.specular(0.5), where) == bc.Realized(True), "the zero network was consulted for a wall with a network")
    require(isinstance(zero_refuses.admit(W.vacuum(), where), bc.Refused), "vacuum does not consult the zero network")
    require(isinstance(zero_refuses.admit(W.specular(0.5) * 0.0 + W.inflow(1.0), where), bc.Refused), "a pure inflow does not consult the zero network")
    source_refuses = make_table(source_admit=refuse).table
    require(source_refuses.admit(W.white(0.5), where) == bc.Realized(True), "the source was consulted for a wall with no source")
    require(isinstance(source_refuses.admit(W.white(0.5) + W.inflow(1.0), where), bc.Refused), "a nonzero source is not consulted")

    # the meet over the units: one tier-3 chain among faithful ones
    approx = bc.RealizerTable("approx", state_cone=bc.StateCone.NONNEGATIVE_FLUX)

    class ApproxSpec(bc.ResponseRealizer, table=approx, realizes=bc.SpecularChain):
        status = bc.Approximated.of("a square cell", 547)

        def realize(self, unit, weight, site_, space):
            return Img((("spec", weight),))

    class FaithfulDiff(bc.ResponseRealizer, table=approx, realizes=bc.DiffuseChain):
        status = bc.Realized(True)

        def realize(self, unit, weight, site_, space):
            return Img((("diff", weight),))

    approx.seal()
    require(approx.admit(mixed, where) == bc.Approximated.of("a square cell", 547), "the meet over a tier-3 chain and a faithful one is not the approximation")
    require(approx.admit(W.white(0.5), where) == bc.Realized(True), "the faithful chain alone is faithful")


# ── (i) the face site ────────────────────────────────────────────────────────


def test_face_site() -> None:
    """``FaceSite`` is constructible from ``(face, coord, dimension, hollow_inner)`` before any method space exists;
    ``is_curved`` and ``site_class`` are derived once: a Cartesian face is ``CART`` and not curved; the outer face of a
    cylinder or sphere is ``RAD_OUTER`` and curved; the inner face of a hollow body is ``RAD_INNER`` and curved; a
    Cartesian-axis face (z) of a radial (r, z) mesh is ``CART`` and not curved. The three classes are exactly the
    members of ``SiteClass``.

    Reddens on: classifying the hollow inner face as the outer; treating the z faces of a radial mesh as curved.
    """
    Cart, Cyl, Sph = CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL
    rows = [
        (("xmin", Cart, 1, False), False, "CART"), (("xmax", Cart, 1, False), False, "CART"), (("ymax", Cart, 2, False), False, "CART"),
        (("xmax", Cyl, 1, False), True, "RAD_OUTER"), (("xmax", Sph, 1, False), True, "RAD_OUTER"),
        (("xmin", Cyl, 1, True), True, "RAD_INNER"), (("xmin", Sph, 1, True), True, "RAD_INNER"),
        (("xmax", Cyl, 2, False), True, "RAD_OUTER"), (("ymin", Cyl, 2, False), False, "CART"), (("ymax", Cyl, 2, False), False, "CART"),
    ]
    seen = set()
    for args, curved, klass in rows:
        s = bc.FaceSite(*args)
        require(s.is_curved is curved, f"{args}: is_curved {s.is_curved}, expected {curved}")
        require(s.site_class is bc.SiteClass[klass], f"{args}: site_class {s.site_class}, expected {klass}")
        seen.add(s.site_class.name)
    require(seen == {m.name for m in bc.SiteClass} == {"CART", "RAD_OUTER", "RAD_INNER"}, f"the classes are {sorted(m.name for m in bc.SiteClass)}, the rows saw {sorted(seen)}")
    require(bc.FaceSite("xmax", Cart, 1, False) == bc.FaceSite("xmax", Cart, 1, False), "FaceSite is a value")
    with pytest.raises(dataclasses.FrozenInstanceError):
        bc.FaceSite("xmax", Cart, 1, False).face = "xmin"  # type: ignore[misc]


# ── (j) the source verb ──────────────────────────────────────────────────────


@pytest.mark.rests_on(_DOOR)
def test_the_source_verb() -> None:
    """``sample_inflow(terms, site, space)`` is ``require`` then the registered ``SourceRealizer.sample``; the default
    ``sample`` is ``sample_terms(terms, space)``; a test-local override is the one called (it receives the terms, the
    site and the space); a refused source raises ``BoundaryRefused`` BEFORE the sampler is read.

    Reddens on: a ``sample_inflow`` that samples before ``require``; one that bypasses the registered realizer; a
    default body other than ``sample_terms``.
    """
    space = function_space(3)
    terms = bc.InflowTerms(2.0, ())
    built = make_table()
    q = built.table.sample_inflow(terms, site("CART"), space)
    require(np.array_equal(q, bc.sample_terms(terms, space)) and np.array_equal(q, 2.0 * np.ones(3)), f"the default sample is {q}")
    require(len(built.spies.sample) == 1 and built.spies.sample[0][2] is space, "the registered realizer's sample was not the one called")

    marker = np.array([7.0, 8.0, 9.0])
    override = make_table(sample_override=lambda terms_, site_, space_: marker)
    out = override.table.sample_inflow(terms, site("CART"), space)
    require(np.array_equal(out, marker) and override.spies.sample[0][0] is terms and override.spies.sample[0][1] == site("CART"), "the override was not called with (terms, site, space)")

    sampler_calls: list[object] = []

    class Counting:
        def sample(self, space_):
            sampler_calls.append(space_)
            return np.ones(space_.shape)

    sampled = bc.InflowTerms(0.0, ((1.0, bc.FaceFrameInflow(Counting())),))
    with pytest.raises(bc.BoundaryRefused) as caught:
        built.table.sample_inflow(sampled, site("RAD_OUTER"), space)
    require(caught.value.admission.axis is bc.RefusalAxis.PLUMBING and sampler_calls == [], f"the sampler was read before the refusal: {sampler_calls}")
    require(len(built.spies.sample) == 1, "the realizer's sample was called for a refused source")
    require(np.array_equal(built.table.sample_inflow(sampled, site("CART"), space), np.ones(3)) and len(sampler_calls) == 1, "an admitted sampled source is not sampled once")
