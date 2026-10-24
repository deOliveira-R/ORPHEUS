r"""Content identity of the transport boundary-conditions package: the condition, the network, the terms, the
displacement, the declared map and the binding (wave 2, slice S2a, gate 15).

Spec: ``.claude/plans/wave2_boundary_carve.md`` section 5.0 ("the content-identity rosters") and 6.2 S2a gate 15, with
``.claude/plans/reference_p1_spec.md`` section 1.5 for the legs. The encoder and the mixin are gated in
``tests/gates/values/test_content_identity.py``, whose population gate S5.9 ("every concrete ``ContentIdentity`` class
has a roster entry") imports :data:`ROSTER` from here (the aggregator edit of the same commit).

The six TYPES whose identity is their content get a full entry (equal content is one value, a change to any part
moves the digest, the pickle round trip): ``WallCondition``, ``ResponseNetwork``, ``WallDisplacement``,
``BoundaryConditions``, ``Binding``, ``InflowTerms``. The frozen dataclasses that are only PARTS (the factor and
chain markers, the zero marker, the two profiles, the admission members, ``FaceSite``) take
``content_identity=False``: they keep their dataclass equality and are content through the encoder's frozen-dataclass
rule, so a part carrying one moves its holder's digest.

Equal-content pairs worth their rows (each names the law it witnesses): the albedo twin gone (``specular(0.0)``,
``white(0.0)`` and ``vacuum()``: law L-B), a facade against the constructor (``white(a)`` against
``from_amplitudes(0.0, a)``: law L-G), the order of a network's weights, a zero weight, ``-0.0``.

First red, at the base tip ``4ec95bc6``: the module fails to import (``ModuleNotFoundError`` on the package). Once the
package exists the rows red on the legacy behaviour each rejects: a network that keeps a zero weight (the zero pair),
a map that keeps the caller's order (the reordering pair), a binding that does not hold its symmetry as a part (the
symmetry leg), a displacement whose digest equals its condition's (gate 6's row, not here), a profile with no
content moving its holder (the sampler leg of ``InflowTerms``).

Mutations (on a COPY of the module): drop ``symmetry`` from ``Binding`` (the part-population row reds); make a
``ResponseNetwork`` store the caller's order; declare ``InflowTerms`` with ``eq=True`` (the aggregator's
``test_s5_9`` and ``check_population`` red).
"""

from __future__ import annotations

import inspect
from typing import get_args

import numpy as np
import pytest

from orpheus.geometry.colouring import Symmetries, UncomputedSymmetries
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.deck import DeckTransformation, MirrorDeck
from orpheus.values.content import ContentIdentity, FrozenMapping
from tests.gates._content_identity_helpers import (
    Entry,
    check_equal_pair,
    check_perturbation,
    check_pickle,
    check_population,
    leg,
    pair_ids,
    param_id,
    perturbation_ids,
    require,
)
from tests.gates.transport._boundary_conditions_helpers import PACKAGE, ProfileSampler, bc, slab_mesh, walk_subclasses, weights

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/transport/test_content_identity_transport.py"
_ENCODER = "tests/gates/values/test_content_identity.py"
_S51 = f"{_ENCODER}::test_s5_1_digests_and_hashes_are_seed_stable"

W = bc.WallCondition


# ── the full entries ─────────────────────────────────────────────────────────


def _wall(network=None, inflow=None):
    return W(network if network is not None else weights(specular=0.5), inflow if inflow is not None else bc.InflowTerms.empty())


_WALL_CONDITION = Entry(
    cls=bc.WallCondition, base=_wall, parts=("network", "inflow"),
    perturb={
        "network": (
            leg("another weight", lambda: _wall(weights(specular=0.25))),
            leg("the other chain", lambda: _wall(weights(diffuse=0.5))),
            leg("both chains", lambda: _wall(weights(specular=0.5, diffuse=0.25))),
        ),
        "inflow": (
            leg("a uniform source", lambda: _wall(inflow=bc.InflowTerms(1.0, ()))),
            leg("a sampled source", lambda: _wall(inflow=bc.InflowTerms(0.0, ((0.5, bc.FaceFrameInflow(ProfileSampler(2.0))),)))),
        ),
    },
    pairs=(
        ("two builds", _wall, _wall),
        ("the facade against the constructor", lambda: W.specular(0.5), _wall),
        ("white against its amplitudes (law L-G)", lambda: W.white(0.5), lambda: W.from_amplitudes(0.0, 0.5)),
        ("the albedo twin is gone (law L-B)", lambda: W.specular(0.0), lambda: W.white(0.0)),
        ("a zero albedo is vacuum", lambda: W.white(0.0), lambda: W.vacuum()),
        ("a zero inflow is vacuum", lambda: W.uniform(0.0), lambda: W.vacuum()),
        ("the dyadic merge is exact", lambda: W.specular(0.25) + W.specular(0.5), lambda: W.specular(0.75)),
    ),
)

_RESPONSE_NETWORK = Entry(
    cls=bc.ResponseNetwork, base=lambda: weights(specular=0.5, diffuse=0.25), parts=("weights",),
    perturb={
        "weights": (
            leg("a weight moved", lambda: weights(specular=0.5, diffuse=0.125)),
            leg("a chain absent", lambda: weights(specular=0.5)),
            leg("the weights swapped between the chains", lambda: weights(specular=0.25, diffuse=0.5)),
            leg("a signed weight", lambda: weights(specular=0.5, diffuse=-0.25)),
        ),
    },
    pairs=(
        ("two builds", lambda: weights(specular=0.5, diffuse=0.25), lambda: weights(specular=0.5, diffuse=0.25)),
        ("the order of the weights", lambda: weights(specular=0.5, diffuse=0.25),
         lambda: bc.ResponseNetwork(FrozenMapping({bc.DiffuseChain(): 0.25, bc.SpecularChain(): 0.5}))),
        ("a zero weight is not stored", lambda: weights(specular=0.5, diffuse=0.0), lambda: weights(specular=0.5)),
        ("-0.0 is 0.0", lambda: weights(specular=-0.0), lambda: weights()),
    ),
)


def _displacement(*terms):
    return bc.WallDisplacement.of(*(terms or ((2.0, W.specular(0.6)),)))


_WALL_DISPLACEMENT = Entry(
    cls=bc.WallDisplacement, base=_displacement, parts=("network", "inflow"),
    perturb={
        "network": (
            leg("another weight", lambda: _displacement((2.0, W.specular(0.5)))),
            leg("the other chain", lambda: _displacement((2.0, W.white(0.6)))),
        ),
        "inflow": (
            leg("a uniform source", lambda: _displacement((2.0, W.specular(0.6)), (1.0, W.uniform(1.0)))),
        ),
    },
    pairs=(
        ("two builds", _displacement, _displacement),
        ("the order of the terms", lambda: _displacement((1.0, W.specular(0.25)), (2.0, W.uniform(0.5))),
         lambda: _displacement((2.0, W.uniform(0.5)), (1.0, W.specular(0.25)))),
    ),
)


def _declared(mapping=None):
    return bc.BoundaryConditions(mapping if mapping is not None else {0: W.vacuum(), 1: W.white(0.5)})


_BOUNDARY_CONDITIONS = Entry(
    cls=bc.BoundaryConditions, base=_declared, parts=("conditions",),
    perturb={
        "conditions": (
            leg("a condition moved", lambda: _declared({0: W.vacuum(), 1: W.white(0.25)})),
            leg("a tag moved", lambda: _declared({0: W.vacuum(), 2: W.white(0.5)})),
            leg("a tag absent", lambda: _declared({0: W.vacuum()})),
        ),
    },
    pairs=(
        ("two builds", _declared, _declared),
        ("the order of declaration", _declared, lambda: bc.BoundaryConditions({1: W.white(0.5), 0: W.vacuum()})),
        ("an np.integer tag", _declared, lambda: bc.BoundaryConditions({np.int64(0): W.vacuum(), 1: W.white(0.5)})),  # type: ignore[dict-item]
    ),
)


def _binding(
    *,
    decks: FrozenMapping[str, DeckTransformation] | None = None,
    tags: FrozenMapping[str, int] | None = None,
    declared: bc.BoundaryConditions | None = None,
    symmetry: Symmetries | None = None,
):
    return bc.Binding(
        decks=decks if decks is not None else FrozenMapping({"xmin": MirrorDeck()}),
        tags=tags if tags is not None else FrozenMapping({"xmax": 0}),
        declared=declared if declared is not None else bc.BoundaryConditions({0: W.white(0.5)}),
        symmetry=symmetry if symmetry is not None else slab_mesh({"xmin": MirrorDeck(), "xmax": 0}).symmetry,
    )


_BINDING = Entry(
    cls=bc.Binding, base=_binding, parts=("decks", "tags", "declared", "symmetry"),
    perturb={
        "decks": (
            leg("no deck", lambda: _binding(decks=FrozenMapping({}), tags=FrozenMapping({"xmin": 0, "xmax": 0}), symmetry=slab_mesh({"xmin": 0, "xmax": 0}).symmetry), "tags", "symmetry"),
            leg("a second deck", lambda: _binding(decks=FrozenMapping({"xmin": MirrorDeck(), "ymin": MirrorDeck()}), symmetry=UncomputedSymmetries("Mesh2D")), "symmetry"),
        ),
        "tags": (
            leg("a second tagged face", lambda: _binding(tags=FrozenMapping({"xmax": 0, "ymax": 0}), symmetry=UncomputedSymmetries("Mesh2D")), "symmetry"),
        ),
        "declared": (
            leg("another condition", lambda: _binding(declared=bc.BoundaryConditions({0: W.white(0.25)}))),
        ),
        "symmetry": (
            leg("the mirror group", lambda: _binding(symmetry=slab_mesh({"xmin": 0, "xmax": 0}).symmetry)),
            leg("an uncomputed group", lambda: _binding(symmetry=UncomputedSymmetries("Mesh2D"))),
        ),
    },
    pairs=(
        ("two builds", _binding, _binding),
        ("bind against the constructor", _binding, lambda: bc.bind(slab_mesh({"xmin": MirrorDeck(), "xmax": 0}), {0: W.white(0.5)})),
    ),
)


def _terms(*, uniform: float = 1.0, sampled: tuple[tuple[float, bc.InflowProfile], ...] | None = None):
    return bc.InflowTerms(uniform=uniform, sampled=sampled if sampled is not None else ((0.5, bc.FaceFrameInflow(ProfileSampler(2.0))),))


_INFLOW_TERMS = Entry(
    cls=bc.InflowTerms, base=_terms, parts=("uniform", "sampled"),
    perturb={
        "uniform": (leg("another uniform weight", lambda: _terms(uniform=2.0)),),
        "sampled": (
            leg("another weight", lambda: _terms(sampled=((0.25, bc.FaceFrameInflow(ProfileSampler(2.0))),))),
            leg("another frame", lambda: _terms(sampled=((0.5, bc.LabFrameInflow(ProfileSampler(2.0))),))),
            leg("another sampler", lambda: _terms(sampled=((0.5, bc.FaceFrameInflow(ProfileSampler(3.0))),))),
            leg("no sampled term", lambda: _terms(sampled=())),
        ),
    },
    pairs=(
        ("two builds", _terms, _terms),
        ("a zero weight is not stored", lambda: _terms(sampled=((0.0, bc.FaceFrameInflow(ProfileSampler(2.0))),)), lambda: _terms(sampled=())),
        ("-0.0 is 0.0", lambda: _terms(uniform=-0.0, sampled=()), lambda: _terms(uniform=0.0, sampled=())),
    ),
)

# ── the parts: frozen dataclasses with dataclass equality ────────────────────


def _part(cls):
    """A fieldless part: ``base`` builds it, there is nothing to perturb, two builds are equal."""
    return Entry(cls=cls, base=cls, parts=(), perturb={}, pairs=(("two builds", cls, cls),), content_identity=False)


_PARTS = (
    _part(bc.MirrorPermutation), _part(bc.CurrentRetraction), _part(bc.IsotropicSection),
    _part(bc.SpecularChain), _part(bc.DiffuseChain), _part(bc.ZeroNetwork),
    Entry(
        cls=bc.FaceFrameInflow, base=lambda: bc.FaceFrameInflow(ProfileSampler(1.0)), parts=("sampler",),
        perturb={"sampler": (leg("another sampler", lambda: bc.FaceFrameInflow(ProfileSampler(2.0))),)},
        pairs=(("two builds", lambda: bc.FaceFrameInflow(ProfileSampler(1.0)), lambda: bc.FaceFrameInflow(ProfileSampler(1.0))),),
        content_identity=False,
    ),
    Entry(
        cls=bc.LabFrameInflow, base=lambda: bc.LabFrameInflow(ProfileSampler(1.0)), parts=("sampler",),
        perturb={"sampler": (leg("another sampler", lambda: bc.LabFrameInflow(ProfileSampler(2.0))),)},
        pairs=(("two builds", lambda: bc.LabFrameInflow(ProfileSampler(1.0)), lambda: bc.LabFrameInflow(ProfileSampler(1.0))),),
        content_identity=False,
    ),
    Entry(
        cls=bc.Realized, base=lambda: bc.Realized(True), parts=("faithful",),
        perturb={"faithful": (leg("tier 2", lambda: bc.Realized(False)),)},
        pairs=(("two builds", lambda: bc.Realized(True), lambda: bc.Realized(True)),), content_identity=False,
    ),
    Entry(
        cls=bc.Approximated, base=lambda: bc.Approximated.of("a square cell", 547), parts=("approximations",),
        perturb={
            "approximations": (
                leg("another name", lambda: bc.Approximated.of("a torus", 547)),
                leg("another issue", lambda: bc.Approximated.of("a square cell", 548)),
                leg("a second pair", lambda: bc.Approximated(frozenset({("a square cell", 547), ("a torus", 548)}))),
            ),
        },
        pairs=(("two builds", lambda: bc.Approximated.of("a square cell", 547), lambda: bc.Approximated.of("a square cell", 547)),),
        content_identity=False,
    ),
    Entry(
        cls=bc.Refused, base=lambda: bc.Refused(bc.RefusalAxis.TOPOLOGY, "no coupling", 547), parts=("axis", "reason", "issue"),
        perturb={
            "axis": (leg("another axis", lambda: bc.Refused(bc.RefusalAxis.STATE_CONE, "no coupling", 547)),),
            "reason": (leg("another reason", lambda: bc.Refused(bc.RefusalAxis.TOPOLOGY, "other", 547)),),
            "issue": (leg("another issue", lambda: bc.Refused(bc.RefusalAxis.TOPOLOGY, "no coupling", 548)),
                      leg("no issue", lambda: bc.Refused(bc.RefusalAxis.TOPOLOGY, "no coupling", None))),
        },
        pairs=(("two builds", lambda: bc.Refused(bc.RefusalAxis.TOPOLOGY, "no coupling", 547), lambda: bc.Refused(bc.RefusalAxis.TOPOLOGY, "no coupling", 547)),),
        content_identity=False,
    ),
    Entry(
        cls=bc.FaceSite, base=lambda: bc.FaceSite("xmax", CoordSystem.CYLINDRICAL, 1), parts=("face", "coord", "dimension"),
        perturb={
            "face": (leg("another face", lambda: bc.FaceSite("ymax", CoordSystem.CYLINDRICAL, 1)),),
            "coord": (leg("another chart", lambda: bc.FaceSite("xmax", CoordSystem.SPHERICAL, 1)),),
            "dimension": (leg("another dimension", lambda: bc.FaceSite("xmax", CoordSystem.CYLINDRICAL, 2)),),
        },
        pairs=(("two builds", lambda: bc.FaceSite("xmax", CoordSystem.CYLINDRICAL, 1), lambda: bc.FaceSite("xmax", CoordSystem.CYLINDRICAL, 1)),),
        content_identity=False,
    ),
)

ROSTER: tuple[Entry, ...] = (
    _WALL_CONDITION, _RESPONSE_NETWORK, _WALL_DISPLACEMENT, _BOUNDARY_CONDITIONS, _BINDING, _INFLOW_TERMS, *_PARTS,
)


# ── the population is the TYPE, not a hand list ──────────────────────────────


def test_the_roster_covers_the_package() -> None:
    """Every concrete ``ContentIdentity`` class DEFINED in ``orpheus.transport.boundary_conditions`` has a
    content-identity entry, and every chain and factor class has an entry (as a part). The denominator is the set of
    classes the interpreter holds (``__subclasses__`` walked from ``ContentIdentity``), filtered to the package's
    modules and to the concrete ones, so a type added later and not rostered reds here and in the aggregator.

    Population: at least the six types of the section (so a walk that finds nothing, or a filter that drops them,
    cannot pass).
    """
    found = {c for c in walk_subclasses(ContentIdentity) if c.__module__.startswith(PACKAGE) and not inspect.isabstract(c)}
    names = {c.__name__ for c in found}
    wanted = {"WallCondition", "ResponseNetwork", "WallDisplacement", "BoundaryConditions", "Binding", "InflowTerms"}
    require(wanted <= names, f"the walk lost {sorted(wanted - names)}: it found {sorted(names)}")
    rostered = {e.cls for e in ROSTER if e.content_identity}
    require(found <= rostered, f"content types with no roster entry: {sorted(c.__qualname__ for c in found - rostered)}")
    require(rostered <= found, f"roster entries that are not content types of the package: {sorted(c.__qualname__ for c in rostered - found)}")
    parts = {e.cls for e in ROSTER if not e.content_identity}
    for family in (get_args(bc.Chain), get_args(bc.Factor)):
        require(set(family) <= parts, f"chain or factor classes with no part entry: {sorted(c.__name__ for c in set(family) - parts)}")
    require(len(ROSTER) == 18, f"the roster holds {len(ROSTER)} entries, not the 18 it was written with (6 types, 12 parts)")


def test_the_content_types_take_the_mixins_dunders() -> None:
    """The conformance rule of the aggregator's S5.9 (``tests/gates/values/test_content_identity.py``), restated for
    the package so that slice S2a's own targeted set reads it: that file cannot be collected until slice S5 (its
    rosters import the retired ``BC``). Every concrete ``ContentIdentity`` class of the package takes ``__eq__`` and
    ``__hash__`` from the mixin, and is a frozen dataclass declared with ``eq=False`` (a dataclass with ``eq=True``
    generates its own ``__eq__``/``__hash__`` and shadows the digest).

    Reddens on: ``InflowTerms`` (or any of the six) declared ``@dataclass(frozen=True)``.
    """
    found = sorted((c for c in walk_subclasses(ContentIdentity) if c.__module__.startswith(PACKAGE) and not inspect.isabstract(c)), key=lambda c: c.__name__)
    require(len(found) == 6, f"the package defines {len(found)} concrete content types, not the 6 of the section: {[c.__name__ for c in found]}")
    problems = []
    for cls in found:
        owner = lambda name: next(k for k in cls.__mro__ if name in vars(k))  # noqa: E731
        if owner("__eq__") is not ContentIdentity:
            problems.append(f"{cls.__name__}.__eq__ comes from {owner('__eq__').__qualname__}")
        if cls.__hash__ is not ContentIdentity.__hash__:
            problems.append(f"{cls.__name__}.__hash__ comes from {owner('__hash__').__qualname__}")
        params = getattr(cls, "__dataclass_params__", None)
        if params is None or params.eq or not params.frozen:
            problems.append(f"{cls.__name__} is not a frozen eq=False dataclass: {params}")
    require(not problems, "\n".join(problems))


@pytest.mark.rests_on(_S51)
@pytest.mark.parametrize("entry", ROSTER, ids=lambda e: e.id)
def test_s5_3_population(entry: Entry) -> None:
    check_population(entry)


@pytest.mark.rests_on(f"{_HERE}::test_s5_3_population")
@pytest.mark.parametrize(
    "entry,part,the_leg",
    [pytest.param(e, p, l, id=param_id(e.id, p, l[0])) for e, p, l in perturbation_ids(ROSTER)],
)
def test_s5_3_a_moved_part_moves_the_digest(entry: Entry, part: str, the_leg) -> None:
    check_perturbation(entry, part, the_leg)


@pytest.mark.rests_on(_S51)
@pytest.mark.parametrize(
    "entry,pair", [pytest.param(e, p, id=param_id(e.id, p[0])) for e, p in pair_ids(ROSTER)]
)
def test_s5_2_equal_content_is_one_value(entry: Entry, pair) -> None:
    check_equal_pair(entry, pair)


@pytest.mark.rests_on(_S51)
@pytest.mark.parametrize("entry", [e for e in ROSTER if e.pickles], ids=lambda e: e.id)
def test_s5_2_pickle_round_trip(entry: Entry) -> None:
    check_pickle(entry)
