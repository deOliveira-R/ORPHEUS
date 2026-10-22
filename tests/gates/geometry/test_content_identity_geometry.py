r"""Content identity of the geometry layer: ``BC``, the boundary laws and their
parts, the law algebra, ``StructuredGeometry`` (#405 P1 step 5).

Spec: ``.claude/plans/reference_p1_spec.md`` §1.5, the geometry legs of S5.2,
S5.3 and S5.4, and S5.5 (a law with no content refuses the digest and is equal
only to itself). The encoder and the mixin are gated in
``tests/gates/values/test_content_identity.py``.

First red, measured on ``main`` ``1dc31163``: the module fails to import
(``orpheus.values.content`` does not exist). With the encoder alone added,
the rows red on the legacy behaviour each rejects: ``hash(BC.vacuum)`` raises
``TypeError: unhashable type: 'dict'``; ``BC.params`` is a mutable dict, the
shared constants' included; ``hash(StructuredGeometry(...))`` raises when a
boundary is a ``BC``; ``PrescribedInflow`` over a plain source object is equal
to a second law over the same source and hashes by the source's id.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from orpheus.geometry import BC, CoordSystem, StructuredGeometry
from orpheus.geometry.boundary import (
    AlbedoBoundary,
    ConstantInflowSource,
    IsotropicReturn,
    LawScaled,
    LawSum,
    NoSource,
    PeriodicBoundary,
    PrescribedInflow,
    ReflectiveBoundary,
    SpecularReturn,
    VacuumInflow,
    WhiteBoundary,
    ZeroFluxBoundary,
)
from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.transformation import RigidMotion
from orpheus.values.content import ContentlessError, content_digest, encode_exact
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

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/geometry/test_content_identity_geometry.py"
_ENCODER = "tests/gates/values/test_content_identity.py"
_S51 = f"{_ENCODER}::test_s5_1_digests_and_hashes_are_seed_stable"
_S54 = f"{_ENCODER}::TestS54EncoderCanonicalForms"

_SPH, _CYL = CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL


def _bc() -> BC:
    return BC("albedo", {"albedo": 0.3, "beta": 2.0})


def _geometry(**over) -> StructuredGeometry:
    """A hollow sphere, two intervals, a typed inner law and a ``BC`` outer law."""
    kw = dict(
        coord=_SPH, breakpoints=(0.5, 1.0, 2.0), mat_ids=(0, 1),
        boundaries=(AlbedoBoundary(0.3), BC.vacuum),
    )
    kw.update(over)
    return StructuredGeometry(**kw)  # type: ignore[arg-type]  # kwargs assembled for the perturbation table


# ── The roster ──────────────────────────────────────────────────────────────

_BC = Entry(
    cls=BC,
    base=_bc,
    parts=("kind", "params"),
    perturb={
        "kind": (leg("another kind", lambda: BC("partial", {"albedo": 0.3, "beta": 2.0})),),
        "params": (
            leg("a value moved one ulp", lambda: BC("albedo", {"albedo": float(np.nextafter(0.3, 1.0)), "beta": 2.0})),
            leg("a key added", lambda: BC("albedo", {"albedo": 0.3, "beta": 2.0, "gamma": 0.0})),
            leg("a key renamed", lambda: BC("albedo", {"alpha": 0.3, "beta": 2.0})),
            leg("the values swapped between keys", lambda: BC("albedo", {"albedo": 2.0, "beta": 0.3})),
        ),
    },
    pairs=(
        ("two builds", _bc, _bc),
        ("params insertion order", _bc, lambda: BC("albedo", {"beta": 2.0, "albedo": 0.3})),
        ("int vs float", lambda: BC("albedo", {"albedo": 1}), lambda: BC("albedo", {"albedo": 1.0})),
        ("numpy scalar", lambda: BC("albedo", {"albedo": 0.3}), lambda: BC("albedo", {"albedo": np.float64(0.3)})),
        ("the constant vs a build", lambda: BC.vacuum, lambda: BC("vacuum")),
        ("default params vs empty", lambda: BC("vacuum"), lambda: BC("vacuum", {})),
        ("signed zero", lambda: BC("albedo", {"albedo": -0.0}), lambda: BC("albedo", {"albedo": 0.0})),
    ),
)

_VACUUM = Entry(
    cls=VacuumInflow, base=VacuumInflow, parts=(), perturb={},
    pairs=(("two builds", VacuumInflow, VacuumInflow),),
)
_REFLECTIVE = Entry(
    cls=ReflectiveBoundary, base=lambda: ReflectiveBoundary("x"), parts=("axis",),
    perturb={"axis": (leg("y", lambda: ReflectiveBoundary("y")),)},
    pairs=(("two builds", lambda: ReflectiveBoundary("x"), lambda: ReflectiveBoundary(axis="x")),),
)
_WHITE = Entry(
    cls=WhiteBoundary, base=lambda: WhiteBoundary("x", +1, 1.0),
    parts=("axis", "outward_sign", "albedo"),
    perturb={
        "axis": (leg("y", lambda: WhiteBoundary("y", +1, 1.0)),),
        "outward_sign": (leg("-1", lambda: WhiteBoundary("x", -1, 1.0)),),
        "albedo": (leg("0.5", lambda: WhiteBoundary("x", +1, 0.5)),),
    },
    pairs=(("int vs float albedo", lambda: WhiteBoundary("x", 1, 1), lambda: WhiteBoundary("x", 1, 1.0)),),
)
_PERIODIC = Entry(
    cls=PeriodicBoundary, base=lambda: PeriodicBoundary("x"), parts=("axis",),
    perturb={"axis": (leg("y", lambda: PeriodicBoundary("y")),)},
    pairs=(("two builds", lambda: PeriodicBoundary("x"), lambda: PeriodicBoundary("x")),),
)
_ALBEDO = Entry(
    cls=AlbedoBoundary, base=lambda: AlbedoBoundary(0.3), parts=("albedo", "reemission"),
    perturb={
        "albedo": (leg("0.4", lambda: AlbedoBoundary(0.4)),),
        "reemission": (
            leg("specular", lambda: AlbedoBoundary(0.3, SpecularReturn("x"))),
            leg("isotropic", lambda: AlbedoBoundary(0.3, IsotropicReturn("x", +1))),
        ),
    },
    pairs=(
        ("int vs float", lambda: AlbedoBoundary(1), lambda: AlbedoBoundary(1.0)),
        ("numpy scalar", lambda: AlbedoBoundary(0.3), lambda: AlbedoBoundary(np.float64(0.3))),
        ("signed zero", lambda: AlbedoBoundary(-0.0), lambda: AlbedoBoundary(0.0)),
        ("a closure built twice", lambda: AlbedoBoundary(0.3, IsotropicReturn("x", 1)),
         lambda: AlbedoBoundary(0.3, IsotropicReturn("x", +1))),
    ),
)
_PRESCRIBED = Entry(
    cls=PrescribedInflow, base=lambda: PrescribedInflow(NoSource()), parts=("_source",),
    perturb={
        "_source": (
            leg("a constant source", lambda: PrescribedInflow(ConstantInflowSource(1.0))),
        ),
    },
    pairs=(
        ("two builds", lambda: PrescribedInflow(NoSource()), lambda: PrescribedInflow()),
        ("int vs float constant", lambda: PrescribedInflow(ConstantInflowSource(2)),
         lambda: PrescribedInflow(ConstantInflowSource(2.0))),
    ),
)
_ZERO_FLUX = Entry(
    cls=ZeroFluxBoundary, base=ZeroFluxBoundary, parts=(), perturb={},
    pairs=(("two builds", ZeroFluxBoundary, ZeroFluxBoundary),),
)
_LAW_SUM = Entry(
    cls=LawSum, base=lambda: LawSum(AlbedoBoundary(0.3), VacuumInflow()), parts=("a", "b"),
    perturb={
        "a": (leg("another albedo", lambda: LawSum(AlbedoBoundary(0.4), VacuumInflow())),),
        "b": (leg("zero flux", lambda: LawSum(AlbedoBoundary(0.3), ZeroFluxBoundary())),),
    },
    pairs=(("two builds", lambda: LawSum(AlbedoBoundary(0.3), VacuumInflow()),
            lambda: LawSum(AlbedoBoundary(0.3), VacuumInflow())),),
)
_LAW_SCALED = Entry(
    cls=LawScaled, base=lambda: LawScaled(0.5, AlbedoBoundary(0.3)), parts=("scalar", "inner"),
    perturb={
        "scalar": (leg("0.25", lambda: LawScaled(0.25, AlbedoBoundary(0.3))),),
        "inner": (leg("another albedo", lambda: LawScaled(0.5, AlbedoBoundary(0.4))),),
    },
    pairs=(("int vs float scalar", lambda: LawScaled(1, AlbedoBoundary(0.3)),
            lambda: LawScaled(1.0, AlbedoBoundary(0.3))),),
)

# The laws' PARTS: frozen dataclasses that are content through the encoder's
# frozen-dataclass rule and keep their dataclass equality (main's ruling,
# 2026-10-02). Their own fields are the population.
_SPECULAR_RETURN = Entry(
    cls=SpecularReturn, base=lambda: SpecularReturn("x"), parts=("axis",), content_identity=False,
    perturb={"axis": (leg("y", lambda: SpecularReturn("y")),)},
)
_ISOTROPIC_RETURN = Entry(
    cls=IsotropicReturn, base=lambda: IsotropicReturn("x", +1), parts=("axis", "outward_sign"),
    content_identity=False,
    perturb={
        "axis": (leg("y", lambda: IsotropicReturn("y", +1)),),
        "outward_sign": (leg("-1", lambda: IsotropicReturn("x", -1)),),
    },
)
_NO_SOURCE = Entry(cls=NoSource, base=NoSource, parts=(), perturb={}, content_identity=False)
_CONSTANT_SOURCE = Entry(
    cls=ConstantInflowSource, base=lambda: ConstantInflowSource(1.0), parts=("value",),
    content_identity=False,
    perturb={"value": (leg("2.0", lambda: ConstantInflowSource(2.0)),)},
)

_GEOMETRY = Entry(
    cls=StructuredGeometry, base=_geometry,
    parts=("coord", "breakpoints", "mat_ids", "boundaries"),
    perturb={
        "coord": (leg("cylinder", lambda: _geometry(coord=_CYL)),),
        "breakpoints": (
            leg("an interior breakpoint", lambda: _geometry(breakpoints=(0.5, 1.5, 2.0))),
            leg("the outer breakpoint one ulp", lambda: _geometry(breakpoints=(0.5, 1.0, float(np.nextafter(2.0, 3.0))))),
        ),
        "mat_ids": (leg("a material id", lambda: _geometry(mat_ids=(0, 2))),),
        "boundaries": (
            leg("the inner law", lambda: _geometry(boundaries=(AlbedoBoundary(0.4), BC.vacuum))),
            leg("the outer tag vs its typed law", lambda: _geometry(boundaries=(AlbedoBoundary(0.3), VacuumInflow()))),
        ),
    },
    pairs=(
        ("two builds", _geometry, _geometry),
        ("int breakpoints", _geometry, lambda: _geometry(breakpoints=(0.5, 1, 2))),
        ("a BC built, not the constant", _geometry,
         lambda: _geometry(boundaries=(AlbedoBoundary(0.3), BC("vacuum")))),
    ),
)

# ── C1 (#405 P3 step 2c): paste into tests/gates/geometry/test_content_identity_geometry.py ──────────────────────
# Imports to add: from orpheus.geometry.chart import Chart; from orpheus.geometry.chord import ConcentricPartition;
#                 from orpheus.geometry.transformation import RigidMotion; from orpheus.values.content import encode_exact
# and the two entries appended to ROSTER: (..., _GEOMETRY, _RIGID_MOTION, _PARTITION).

_Q90 = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])            # a quarter turn about z, exact


def _motion(**over) -> RigidMotion:
    kw = dict(linear=_Q90, translation=np.array([0.5, -0.25, 1.0]))
    kw.update(over)
    return RigidMotion(**kw)  # type: ignore[arg-type]  # kwargs assembled for the perturbation table


_RIGID_MOTION = Entry(
    cls=RigidMotion, base=_motion, parts=("linear", "translation"),
    perturb={
        "linear": (
            leg("a reflection for the turn", lambda: _motion(linear=np.diag([1.0, 1.0, -1.0]))),
            leg("the other quarter turn", lambda: _motion(linear=_Q90.T)),
        ),
        "translation": (
            leg("one component one ulp", lambda: _motion(translation=np.array([0.5, -0.25, float(np.nextafter(1.0, 2.0))]))),
            leg("components swapped", lambda: _motion(translation=np.array([-0.25, 0.5, 1.0]))),
        ),
    },
    pairs=(
        ("two builds", _motion, _motion),
        ("int arrays", _motion, lambda: _motion(linear=_Q90.astype(int), translation=np.array([0.5, -0.25, 1]))),
        ("signed zero in the translation", lambda: _motion(translation=np.array([-0.0, 0.0, 0.0])),
         lambda: _motion(translation=np.zeros(3))),
        ("the identity, named and built", lambda: RigidMotion.identity(3), lambda: RigidMotion(np.eye(3))),
        ("an omitted translation is the zero one", lambda: RigidMotion(_Q90), lambda: RigidMotion(_Q90, np.zeros(3))),
    ),
)


def _partition(**over) -> ConcentricPartition:
    kw = dict(chart=Chart(_SPH), breakpoints=(0.0, 0.5, 1.5, 2.0), pose=RigidMotion.translation_by(np.array([0.0, 0.0, 1.0])))
    kw.update(over)
    return ConcentricPartition(**kw)  # type: ignore[arg-type]  # kwargs assembled for the perturbation table


_PARTITION = Entry(
    cls=ConcentricPartition, base=_partition, parts=("chart", "breakpoints", "pose"),
    perturb={
        "chart": (leg("cylinder", lambda: _partition(chart=Chart(_CYL))),),
        "breakpoints": (
            leg("an interior breakpoint", lambda: _partition(breakpoints=(0.0, 0.6, 1.5, 2.0))),
            leg("the outer breakpoint one ulp", lambda: _partition(breakpoints=(0.0, 0.5, 1.5, float(np.nextafter(2.0, 3.0))))),
        ),
        "pose": (leg("the identity pose", lambda: _partition(pose=RigidMotion.identity(3))),),
    },
    pairs=(
        ("two builds", _partition, _partition),
        ("int breakpoints", _partition, lambda: _partition(breakpoints=(0, 0.5, 1.5, 2))),
        ("a list of breakpoints", _partition, lambda: _partition(breakpoints=[0.0, 0.5, 1.5, 2.0])),
        ("the default pose is the identity", lambda: ConcentricPartition(Chart(_SPH), (0.0, 1.0)),
         lambda: ConcentricPartition(Chart(_SPH), (0.0, 1.0), RigidMotion.identity(3))),
        ("of a geometry vs built", lambda: ConcentricPartition.of(_geometry()),
         lambda: ConcentricPartition(Chart(_SPH), (0.5, 1.0, 2.0))),
    ),
)




ROSTER: tuple[Entry, ...] = (
    _BC, _VACUUM, _REFLECTIVE, _WHITE, _PERIODIC, _ALBEDO, _PRESCRIBED, _ZERO_FLUX,
    _LAW_SUM, _LAW_SCALED, _SPECULAR_RETURN, _ISOTROPIC_RETURN, _NO_SOURCE,
    _CONSTANT_SOURCE, _GEOMETRY, _RIGID_MOTION, _PARTITION,
)


def test_the_roster_covers_the_law_registry() -> None:
    """The law population is the REGISTRY, not a hand list: every registered
    law has a roster entry (``[M]`` 7 keys at 1dc31163)."""
    from orpheus.geometry.boundary import BoundaryTraceLaw

    registered = set(BoundaryTraceLaw.registry.values())
    covered = {entry.cls for entry in ROSTER}
    require(len(registered) == 7, f"the registry holds {len(registered)} laws, not 7")
    require(registered <= covered, f"laws with no roster entry: {registered - covered}")


# ── S5.3 ─────────────────────────────────────────────────────────────────────


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


# ── S5.2 ─────────────────────────────────────────────────────────────────────


@pytest.mark.rests_on(_S51)
@pytest.mark.parametrize(
    "entry,pair", [pytest.param(e, p, id=param_id(e.id, p[0])) for e, p in pair_ids(ROSTER)]
)
def test_s5_2_equal_content_is_one_value(entry: Entry, pair) -> None:
    """First red with the encoder present: every ``BC`` pair, and the geometry
    pairs, raise ``TypeError: unhashable type: 'dict'`` at ``hash``."""
    check_equal_pair(entry, pair)


@pytest.mark.rests_on(_S51)
@pytest.mark.parametrize("entry", [e for e in ROSTER if e.pickles], ids=lambda e: e.id)
def test_s5_2_pickle_round_trip(entry: Entry) -> None:
    check_pickle(entry)


@pytest.mark.rests_on(_S51)
def test_s5_2_a_bc_cannot_be_mutated() -> None:
    """``params`` is frozen: a write raises and the shared constants cannot
    drift. ``[M]`` today ``BC("albedo", {...}).params["x"] = 1.0`` succeeds,
    and so would a write into ``BC.vacuum.params`` (one dict shared by every
    declaration of vacuum)."""
    bc = _bc()
    with pytest.raises(TypeError):
        bc.params["gamma"] = 1.0  # type: ignore[index]  # the refused write is the subject
    with pytest.raises(TypeError):
        BC.vacuum.params["albedo"] = 0.5  # type: ignore[index]  # idem
    require(dict(BC.vacuum.params) == {}, "the shared constant drifted")
    require(dict(bc.params) == {"albedo": 0.3, "beta": 2.0}, "the params read back")


@pytest.mark.rests_on(_S51)
def test_s5_2_the_tag_and_its_typed_law_stay_distinct() -> None:
    """A ``BC`` tag and the typed law it resolves to are different types, so
    different values (``[M]`` today ``geometry(BC.vacuum) !=
    geometry(VacuumInflow())``; the carve must keep it)."""
    require(BC.vacuum != VacuumInflow(), "a tag is not its typed law")
    require(content_digest(BC.vacuum) != content_digest(VacuumInflow()), "nor its digest")
    require(len({BC.vacuum, VacuumInflow()}) == 2, "a set keeps both")


def test_the_string_arm_is_kept() -> None:
    """The user's ruling (kept by the plan): ``VacuumInflow() == "vacuum"`` and
    ``ReflectiveBoundary() == "reflective"`` stay True. The arm breaks the
    eq/hash contract (``hash(VacuumInflow()) != hash("vacuum")``), which is the
    known, filed defect; this row pins the ruling, not the defect."""
    require(VacuumInflow() == "vacuum", "the vacuum string arm")
    require(ReflectiveBoundary("x") == "reflective", "the reflective string arm")
    require(not (VacuumInflow() == "reflective"), "the arm is keyed by kind")
    require(ReflectiveBoundary("x") != ReflectiveBoundary("y"), "law-vs-law goes by content")


# ── S5.4, the geometry legs ──────────────────────────────────────────────────


@pytest.mark.rests_on(_S54)
@pytest.mark.parametrize(
    "build,fragment",
    [
        pytest.param(lambda: BC("albedo", {"albedo": float("nan")}), r"albedo", id="BC-param"),
        pytest.param(lambda: AlbedoBoundary(float("nan")), r"albedo", id="AlbedoBoundary"),
        pytest.param(lambda: WhiteBoundary("x", 1, float("nan")), r"albedo", id="WhiteBoundary"),
        pytest.param(lambda: ConstantInflowSource(float("nan")), r"value", id="ConstantInflowSource"),
    ],
)
def test_s5_4_a_nan_value_cannot_be_constructed(build, fragment: str) -> None:
    """Parse at the boundary (main's ruling, 2026-10-02): NaN is refused at
    CONSTRUCTION with a ``ValueError`` naming the field or parameter, so ``==``
    and ``hash`` never raise on a constructible value. ``[M]`` at 1dc31163
    ``AlbedoBoundary(nan)`` and ``BC(..., {"albedo": nan})`` construct. The
    builder is called inside ``raises`` and no digest is taken: a refusal
    only at the digest (the encoder's backstop) leaves this row red."""
    with pytest.raises(ValueError, match=fragment):
        build()


@pytest.mark.rests_on(f"{_HERE}::test_s5_4_a_nan_value_cannot_be_constructed")
def test_s5_4_a_nan_value_never_reaches_equality() -> None:
    """The consequence the placement buys: every constructible law compares
    and hashes without raising. Witness: the NaN refusal moved to the digest
    would make the row above red, not this one; this row pins that ``==`` on
    two constructible values is total."""
    a, b = AlbedoBoundary(0.3), BC("albedo", {"albedo": 0.3})
    require((a == a) and not (a == b) and isinstance(hash(a) ^ hash(b), int), "== and hash are total")


@pytest.mark.rests_on(_S51)
@pytest.mark.parametrize(
    "value,param_type",
    [pytest.param("2.0", "str", id="a-str"), pytest.param(True, "bool", id="a-bool")],
)
def test_s5_4_a_bc_param_is_a_real_number(value: object, param_type: str) -> None:
    """``BC`` params are real numbers only (main's ruling, 2026-10-02): a
    ``str`` or a ``bool`` value is refused at construction with a
    ``TypeError`` naming the parameter. ``[M]`` at 1dc31163
    ``BC("a", {"name": "x"})`` constructs; 0 of the 51 literal ``BC(kind,
    {...})`` calls in ``orpheus/``, ``tests/``, ``examples/`` pass a str."""
    with pytest.raises(TypeError, match="beta"):
        BC("albedo", {"albedo": 0.3, "beta": value})  # pyright: ignore[reportArgumentType]  # the refused input is the subject


# ── S5.5: a law with no content ──────────────────────────────────────────────


def _manufactured_source():
    """The in-tree contentless source: the test-side ``_ManufacturedFaceInflow``
    (``tests/gates/sn/verification/analytical/test_mms_declared_inflow.py``), a
    plain class holding a mutable list. ``[M]`` it hashes by id today, so
    ``_law_key`` keeps the law itself and never reaches its ``id()`` arm."""
    from tests.gates.sn.verification.analytical.test_mms_declared_inflow import (
        _ManufacturedFaceInflow,
    )

    case = SimpleNamespace(quadrature=SimpleNamespace(weights=np.ones(2)), n_groups=1)
    return _ManufacturedFaceInflow(case=case, x_face=0.0)


_CONTENTLESS_SOURCES = [
    pytest.param(_manufactured_source, r"_ManufacturedFaceInflow", id="manufactured-inflow"),
    pytest.param(lambda: (lambda space: np.zeros(space.shape)), r"function", id="a-function"),
]


@pytest.mark.rests_on(f"{_HERE}::test_s5_3_population")
@pytest.mark.parametrize("make_source,fragment", _CONTENTLESS_SOURCES)
class TestS55ContentlessLaw:
    """A ``PrescribedInflow`` over a source with no content has no digest, no
    hash, and is equal only to itself (identity is the honest equality of a
    value with no content). ``[M]`` first red with the encoder present: two
    laws over ONE source object compare equal and hash (by the source's id)."""

    def test_the_digest_is_refused_naming_the_part(self, make_source, fragment: str) -> None:
        law = PrescribedInflow(make_source())
        with pytest.raises(ContentlessError, match=rf"_source.*{fragment}|{fragment}.*_source"):
            content_digest(law)
        require(issubclass(ContentlessError, TypeError), "a contentless value is a TypeError")

    def test_it_is_unhashable(self, make_source, fragment: str) -> None:
        law = PrescribedInflow(make_source())
        with pytest.raises(TypeError, match="unhashable"):
            hash(law)

    def test_it_is_equal_only_to_itself(self, make_source, fragment: str) -> None:
        source = make_source()
        law, other = PrescribedInflow(source), PrescribedInflow(source)
        require(law == law, "reflexive")
        require(not (law == other), "two laws over one contentless source are two values")
        require(law != other, "!= agrees")
        require(not (law == PrescribedInflow(NoSource())), "nor equal to a content law")

    def test_it_propagates_to_a_geometry(self, make_source, fragment: str) -> None:
        """A geometry holding the law inherits the refusal, with the path."""
        geometry = StructuredGeometry(
            coord=CoordSystem.CARTESIAN, breakpoints=(0.0, 1.0), mat_ids=(0,),
            boundaries=(PrescribedInflow(make_source()), BC.vacuum),
        )
        with pytest.raises(ContentlessError, match=r"boundaries"):
            content_digest(geometry)
        with pytest.raises(TypeError, match="unhashable"):
            hash(geometry)
        require(geometry == geometry, "reflexive")


@pytest.mark.rests_on(f"{_HERE}::test_s5_3_population")
def test_s5_5_content_sources_digest() -> None:
    """The positive legs: the shipped sources are content."""
    for law in (PrescribedInflow(NoSource()), PrescribedInflow(ConstantInflowSource(1.5))):
        require(len(content_digest(law)) == 32, f"{law!r} digests")
        require(hash(law) == hash(PrescribedInflow(law._source)), f"{law!r} hashes by content")


# ── C1 rows beyond the roster ──────────────────────────────────────────────────────────────────────────────────

_C1 = "C1 (#405 P3 step 2c, scratch/characteristic_architecture/p3/c1/gate_spec.md)"


@pytest.mark.rests_on(f"{_HERE}::test_s5_2_equal_content_is_one_value")
@pytest.mark.parametrize("pair", [
    pytest.param((lambda: _motion(translation=np.array([-0.0, 0.0, 0.0])), lambda: _motion(translation=np.zeros(3))), id="motion-signed-zero"),
    pytest.param((lambda: _motion(linear=_Q90.astype(int)), _motion), id="motion-int-linear"),
    pytest.param((lambda: RigidMotion(_Q90), lambda: RigidMotion(_Q90, np.zeros(3))), id="motion-omitted-translation"),
    pytest.param((lambda: _partition(breakpoints=(0, 0.5, 1.5, 2)), _partition), id="partition-int-breakpoints"),
    pytest.param((lambda: ConcentricPartition(Chart(_SPH), (0.0, 1.0)),
                  lambda: ConcentricPartition(Chart(_SPH), (0.0, 1.0), RigidMotion.identity(3))), id="partition-default-pose"),
    pytest.param((lambda: RigidMotion.inversion(3), lambda: RigidMotion(np.diag([-1.0, -1.0, -1.0]))), id="motion-inversion-signed-zero-linear"),
])
def test_c1_one_value_is_one_memo_key(pair) -> None:
    """[C1-K] Two spellings of one value give one EXACT encoding (``encode_exact``, the traced memo's key), not only one
    content digest: the exact key separates ``-0.0`` from ``0.0``, an int from a float and an int array from a float one,
    so each constructor must canonicalise before it stores. ``RigidMotion`` folds both parts through
    ``canonical_reals`` (float64, ``-0.0`` onto ``0.0``); ``ConcentricPartition`` stores its breakpoints as a tuple of
    floats. The inversion row is the LINEAR part's fold (qa, 2026-10-15): ``RigidMotion.inversion(3)`` produces
    ``-0.0`` off the diagonal, and dropping that fold reddened 0 of 492 rows before this row; it asserts the digest too.

    First reds: the signed-zero row on a ``RigidMotion`` that stores the given translation (today's
    ``__post_init__`` keeps ``-0.0``; its ``_exact_key`` canonicalised at comparison only, ``[M]`` 2026-10-14
    ``probe_q1_simulated.log``); the int row on a partition that stores its breakpoints as given (arm
    ``breakpoints-as-given``)."""
    a, b = pair[0](), pair[1]()
    require(a is not b, "activation: one object")
    require(content_digest(a) == content_digest(b), f"{type(a).__qualname__}: one value, two digests")
    require(encode_exact(a) == encode_exact(b), f"{type(a).__qualname__}: one value, two memo keys")


@pytest.mark.rests_on(f"{_HERE}::test_s5_2_pickle_round_trip")
def test_c1_a_rigid_motion_owns_read_only_arrays() -> None:
    """[C1-RO] A motion's arrays are read-only (a write raises), and they are its OWN: building one from a caller's
    writeable float64 arrays leaves the caller's arrays writeable, and a later write into them does not move the motion.
    The content encoder refuses a writeable array as a part, so without this a partition posed by a motion has no
    digest. First reds: today's tree (``linear`` is writeable, the caller's array is aliased:
    ``np.asarray(linear, dtype=float)`` returns it); arm ``freeze-in-place`` (``setflags`` on the aliased input: the
    caller's array becomes read-only)."""
    q, t = _Q90.copy(), np.array([0.5, -0.25, 1.0])
    m = RigidMotion(q, t)
    for name in ("linear", "translation"):
        part = getattr(m, name)
        require(not part.flags.writeable, f"RigidMotion.{name} is writeable")
        with pytest.raises(ValueError, match="read-only"):
            part[0] = 7.0
    require(q.flags.writeable and t.flags.writeable, "the caller's arrays were frozen in place")
    require(not np.shares_memory(m.linear, q) and not np.shares_memory(m.translation, t), "the motion aliases the caller's arrays")
    q[0, 0], t[0] = 3.0, 3.0
    require(m == RigidMotion(_Q90, np.array([0.5, -0.25, 1.0])), "a write into the caller's array moved the motion")


@pytest.mark.rests_on(f"{_HERE}::test_s5_2_pickle_round_trip")
@pytest.mark.parametrize("linear,translation,fragment", [
    pytest.param(np.eye(3), np.array([np.nan, 0.0, 0.0]), r"the translation, the entry \(0,\) is NaN", id="a-NaN-translation"),
    pytest.param(np.eye(3), np.array([np.inf, 0.0, 0.0]), "translation must be finite", id="an-infinite-translation"),
    pytest.param(np.array([[np.nan, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]), np.zeros(3),
                 r"the linear part, the entry \(0, 0\) is NaN", id="a-NaN-linear-part"),
    pytest.param(np.array([[np.inf, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]), np.zeros(3),
                 "linear part must be finite", id="an-infinite-linear-part"),
])
def test_c1_a_rigid_motion_refuses_a_non_finite_part(linear, translation, fragment: str) -> None:
    """[C1-NF] A non-finite linear part or translation is refused at construction, naming it. Content identity refuses
    NaN (it is not equal to itself), so without this law a NaN motion would construct and then fail at ``hash``, far
    from its cause. The NaN linear part also passes today's orthogonality guard (``nan > atol`` is False).
    Each input is keyed to its own guard's fragment: a NaN is refused by ``canonical_reals`` (BEFORE the orthogonality
    guard, which ``nan > atol`` would pass), an infinity by the finiteness guard. First reds: today's tree (all
    construct; ``[M]`` 2026-10-14 ``probe_q1_simulated.log``, the translation); arm ``finite-guard-removed`` (the two
    infinity rows)."""
    with pytest.raises(ValueError, match=fragment):
        RigidMotion(linear, translation)


@pytest.mark.rests_on(f"{_HERE}::test_s5_2_pickle_round_trip")
def test_c1_one_equality_for_a_rigid_motion() -> None:
    """[C1-EQ] ``RigidMotion`` has one equality, the encoder's: its ``_exact_key`` is retired (X4: the bytes key and the
    digest were two definitions of one equality), and ``-0.0``/``0.0`` and two dimensions behave as the retired key did
    (the dimension is the shape of the encoded arrays). The mixin's dunders are S5.9's leg; this row adds the retirement
    and the two behaviours the retired key carried by hand. First red: today's tree (``_exact_key`` defined)."""
    require(not hasattr(RigidMotion, "_exact_key"), "RigidMotion._exact_key survives beside the encoder")
    require(RigidMotion.identity(2) != RigidMotion.identity(3), "two dimensions compare equal")
    require(len({RigidMotion.identity(2), RigidMotion.identity(3)}) == 2, "a set merged two dimensions")
