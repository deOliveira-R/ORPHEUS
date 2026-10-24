r"""Content identity of the geometry layer: ``StructuredGeometry``, its decks, rigid motions, partitions and colouring
symmetry (#405 P1 step 5).

Spec: ``.claude/plans/reference_p1_spec.md`` §1.5, the geometry legs of S5.2 and S5.3. The encoder and the mixin are
gated in ``tests/gates/values/test_content_identity.py``.

The law rows of the original file (``BC``, the boundary laws, their parts, the law algebra, the contentless source of S5.5, and the S5.4 NaN and non-real-parameter refusals of those types)
retired with ``orpheus/geometry/boundary/`` in WF2.4 of the wave-2 carve: their subject types no longer exist. Their successors
are values of the transport layer (``WallCondition`` and its network), gated in ``tests/gates/transport/``.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem, StructuredGeometry
from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.colouring import ColouringSymmetry, FaceIsometry, UncomputedSymmetries
from orpheus.geometry.deck import MirrorDeck, TranslationDeck
from orpheus.geometry.transformation import RigidMotion
from orpheus.values.content import FrozenMapping, content_digest, encode_exact
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

_SPH, _CYL = CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL


def _geometry(**over) -> StructuredGeometry:
    """A hollow sphere, two intervals, the tags 1 (inner) and 0 (outer)."""
    kw = dict(
        coord=_SPH, breakpoints=(0.5, 1.0, 2.0), mat_ids=(0, 1),
        face_roles=(1, 0),
    )
    kw.update(over)
    return StructuredGeometry(**kw)  # type: ignore[arg-type]  # kwargs assembled for the perturbation table


# ── The roster ──────────────────────────────────────────────────────────────

_GEOMETRY = Entry(
    cls=StructuredGeometry, base=_geometry,
    parts=("coord", "breakpoints", "mat_ids", "face_roles"),
    perturb={
        "coord": (leg("cylinder", lambda: _geometry(coord=_CYL)),),
        "breakpoints": (
            leg("an interior breakpoint", lambda: _geometry(breakpoints=(0.5, 1.5, 2.0))),
            leg("the outer breakpoint one ulp", lambda: _geometry(breakpoints=(0.5, 1.0, float(np.nextafter(2.0, 3.0))))),
        ),
        "mat_ids": (leg("a material id", lambda: _geometry(mat_ids=(0, 2))),),
        "face_roles": (
            leg("the inner tag", lambda: _geometry(face_roles=(2, 0))),
            leg("the outer tag", lambda: _geometry(face_roles=(1, 3))),
            leg("the tags swapped between the faces", lambda: _geometry(face_roles=(0, 1))),
        ),
    },
    pairs=(
        ("two builds", _geometry, _geometry),
        ("int breakpoints", _geometry, lambda: _geometry(breakpoints=(0.5, 1, 2))),
        ("a numpy-integer tag vs an int tag", _geometry, lambda: _geometry(face_roles=(np.int64(1), np.int64(0)))),
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


# ── Wave 2, unit L1: the decks and the colouring symmetry (``wave2_boundary_carve.md`` §5.0, §5.1 gate 10) ─────────
# A deck has no field (the face supplies the plane and the partner), so the two decks are parts-free entries;
# they are distinct values by TYPE (``test_deck.py::test_a_deck_is_its_type_and_has_no_field``).

_MIRROR_DECK = Entry(
    cls=MirrorDeck, base=MirrorDeck, parts=(), perturb={},
    pairs=(("two builds", MirrorDeck, MirrorDeck),),
)
_TRANSLATION_DECK = Entry(
    cls=TranslationDeck, base=TranslationDeck, parts=(), perturb={},
    pairs=(("two builds", TranslationDeck, TranslationDeck),),
)


def _swap(left: str = "xmin", right: str = "xmax") -> FrozenMapping:
    return FrozenMapping({left: right, right: left})


def _face_isometry(**over) -> FaceIsometry:
    """The mirror ``x -> 3 - x`` of a slab ``[0, 3]`` and the permutation it induces on the end faces."""
    kw = dict(motion=RigidMotion.reflection(normal=np.array([1.0]), offset=1.5), faces=_swap())
    kw.update(over)
    return FaceIsometry(**kw)  # type: ignore[arg-type]  # kwargs assembled for the perturbation table


_FACE_ISOMETRY = Entry(
    cls=FaceIsometry, base=_face_isometry, parts=("motion", "faces"),
    perturb={
        "motion": (
            leg("the plane moved one ulp", lambda: _face_isometry(
                motion=RigidMotion.reflection(normal=np.array([1.0]), offset=float(np.nextafter(1.5, 2.0))))),
            leg("a translation for the reflection", lambda: _face_isometry(motion=RigidMotion.translation_by(np.array([3.0])))),
        ),
        "faces": (
            leg("the identity permutation", lambda: _face_isometry(faces=FrozenMapping({"xmin": "xmin", "xmax": "xmax"}))),
            leg("another axis", lambda: _face_isometry(faces=_swap("ymin", "ymax"))),
        ),
    },
    pairs=(
        ("two builds", _face_isometry, _face_isometry),
        ("the permutation's insertion order", _face_isometry,
         lambda: _face_isometry(faces=FrozenMapping({"xmax": "xmin", "xmin": "xmax"}))),
    ),
)


def _identity_element() -> FaceIsometry:
    return FaceIsometry(RigidMotion.identity(1), FrozenMapping({"xmin": "xmin", "xmax": "xmax"}))


def _colouring(**over) -> ColouringSymmetry:
    kw = dict(elements=(_identity_element(), _face_isometry()))
    kw.update(over)
    return ColouringSymmetry(**kw)  # type: ignore[arg-type]  # kwargs assembled for the perturbation table


_COLOURING_SYMMETRY = Entry(
    cls=ColouringSymmetry, base=_colouring, parts=("elements",),
    perturb={
        "elements": (
            leg("the identity alone", lambda: _colouring(elements=(_identity_element(),))),
            leg("another mirror", lambda: _colouring(elements=(_identity_element(), _face_isometry(
                motion=RigidMotion.reflection(normal=np.array([1.0]), offset=2.0))))),
        ),
    },
    pairs=(("two builds", _colouring, _colouring),),
)
_UNCOMPUTED_SYMMETRIES = Entry(
    cls=UncomputedSymmetries, base=lambda: UncomputedSymmetries("Mesh2D"), parts=("carrier",),
    perturb={"carrier": (leg("another carrier", lambda: UncomputedSymmetries("Mesh3D")),)},
    pairs=(("two builds", lambda: UncomputedSymmetries("Mesh2D"), lambda: UncomputedSymmetries("Mesh2D")),),
)


ROSTER: tuple[Entry, ...] = (
    _GEOMETRY, _RIGID_MOTION, _PARTITION,
    _MIRROR_DECK, _TRANSLATION_DECK, _FACE_ISOMETRY, _COLOURING_SYMMETRY, _UNCOMPUTED_SYMMETRIES,
)


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
