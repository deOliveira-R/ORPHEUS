r"""A curved periodic declaration is a geometry ``DeckError``, and the reference's own guards stay live for a directly built ``Walls`` (wave 2, S5; plan 6.5 gate 4).

Specified by the test-architect (2026-10-23, ``.claude/plans/wave2_boundary_carve.md`` 6.5). Before the carve the reference refused
a wrap on a radial body (and an unpaired wrap) inside ``Walls.__post_init__`` (``NotImplementedError``), after the geometry had
accepted a ``PeriodicBoundary`` law on any face. After slice S1 the geometry itself cannot say it: a ``TranslationDeck`` on a
radial face is a ``DeckError`` (a translation maps a radial chart's level sets onto none of them), an unpaired wrap is a
``DeckError`` (a translation identifies two faces, so both carry it), and a ``MirrorDeck`` on a radial face is a ``DeckError`` too
(the radial mirror is not a symmetry of the radial chart's orbit space: the specular wall of a curved body is the TAG bound to the
amplitudes ``(1, 0)``, the carve's meaning change of Appendix B). So the three refusal rows of
``test_characteristic_walls.py`` (``wrap_on_hollow_sphere``, ``wrap_on_hollow_cylinder``, ``wrap_on_solid_sphere``, lines 212-216)
and the two unpaired-wrap rows MOVE: their subject is no longer constructible at the reference, so they are re-posed here at the
geometry, and the old rows are DEAD (deleted, not repaired by passing a new argument; ``test-architect`` section 4).

**Land status.** The geometry rows are GREEN at the base tip (``de354496``): S1 landed the ``DeckError``. They are a REGRESSION fence
for the refusals the reference stopped making, not a first red: ``[M]`` ``python -O -m pytest`` on this file at the base reads 9
passed, 2 failed (the two rows that call the reference with its new signatures). The named mutation that gives each geometry row its
tooth, applied to the carved tree in place of a first red (``instrument-doctrine`` X1): drop the refusal of a ``TranslationDeck`` on a
radial axis from ``orpheus/geometry/deck.py::parse_face_role_on`` (the ``DeckError`` raised at ``deck.py:370``, which reads
``admitted_decks(coord, axis)``): the radial rows redden; drop ``check_pairings_close`` (``deck.py:342-353``): the unpaired rows
redden; make ``admitted_decks`` of a radial chart admit ``MirrorDeck``: the curved-mirror rows redden. The positive legs (a tag on
the same faces builds; a wrap or a mirror on the Cartesian faces builds) are what keep a refusal that refuses everything from
passing.

**The battery, measured** (``[M]`` 2026-10-23, ``python -O -m pytest -p`` plugin rebinding in every ``sys.modules`` binding, rebind count
printed; scratch ``test-architect`` scratchpad ``ta_s5/mut_deck.py``; the two rows that call the reference are red at the base by
signature and are excluded from the reading). A1 (``admitted_decks`` of a radial chart admits the decks): the two hollow-body wrap rows and
the two mirror rows redden, 4 of 9 geometry rows. A2 (``check_pairings_close`` a no-op): the unpaired-slab row reddens, 1 of 9. A3
(positive control: every tag refused): all 9 redden. A4 (A1 and A2 together): all 9 redden. The reading that matters: a wrap on a
SOLID radial body, and a single wrap face on a hollow body, are refused by TWO independent guards (the chart's admission AND the pairing
closure, because a lone face has no partner), so no single-guard arm reddens those four rows and only A4 does: they pin the refusal,
not either guard alone (a twin predicate guards the path, ``vv-principles`` #17(d)), and the hollow-body both-faces rows plus the
mirror rows are what isolate the admission guard.

**The surviving guard.** ``Walls.__post_init__`` keeps its two structural checks (an unpaired wrap, a wrap on a radial chart),
because a ``Walls`` is a public value that can be built directly, not only through ``Walls.of``, and the guard is live there: it is
NOT DEAD, and the row below pins it on the new ``Wall`` signature (red at the base, on the old one). The refusal is
``NotImplementedError`` with the pre-carve fragments ("periodic wrap on a cylinder or sphere", "whose partner does not wrap back"), kept
because the message strings are an API once a test pins them (``retirement-audit`` B.11).
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.gates._content_identity_helpers import require

pytestmark = pytest.mark.foundation

_HERE = "tests/gates/derivations/test_curved_wraps_are_deck_errors.py::"
_HOLLOW = (0.4, 1.0, 2.0)
_SOLID = (0.0, 1.0, 2.0)


def _deck_error():
    from orpheus.geometry.deck import DeckError

    return DeckError


def _radial(coord: str, breakpoints, mats, **faces):
    from orpheus.geometry import StructuredGeometry

    build = StructuredGeometry.sphere if coord == "sphere" else StructuredGeometry.cylinder
    return build(breakpoints, mats, **faces)


@pytest.mark.parametrize(
    ("coord", "breakpoints", "mats", "faces"),
    [
        ("sphere", _HOLLOW, (0, 1), ("inner", "outer")),
        ("cylinder", _HOLLOW, (0, 1), ("inner", "outer")),
        ("sphere", _SOLID, (0, 1), ("outer",)),
        ("cylinder", _SOLID, (0, 1), ("outer",)),
    ],
    ids=["hollow-sphere", "hollow-cylinder", "solid-sphere", "solid-cylinder"],
)
def test_a_wrap_on_a_radial_face_is_a_deck_error_at_the_geometry(coord, breakpoints, mats, faces) -> None:
    """``TranslationDeck`` on every radial face of a body: ``DeckError`` at ``StructuredGeometry`` construction, never at ``Walls``.

    The first three ids are the three re-posed rows of ``test_characteristic_walls.py:212-216``; the fourth completes the
    chart set (a solid cylinder, which the old file had no row for). Positive leg inside the row: the same body with its
    faces as TAGS builds, so the refusal is the deck's and not the body's."""
    from orpheus.geometry import TranslationDeck

    _radial(coord, breakpoints, mats, **{f: 0 for f in faces})
    with pytest.raises(_deck_error()):
        _radial(coord, breakpoints, mats, **{f: TranslationDeck() for f in faces})


@pytest.mark.parametrize("face", ["inner", "outer"])
def test_a_single_wrap_face_on_a_hollow_body_is_a_deck_error(face: str) -> None:
    """One ``TranslationDeck`` on a hollow body (the other face a tag) is refused: a radial face admits no translation at all.

    The old ``wrap_on_hollow_*`` rows wrapped BOTH faces, which cannot tell a per-face refusal from a pairing refusal; this row
    wraps one, so the refusal is the face's admissibility."""
    from orpheus.geometry import TranslationDeck

    other = "outer" if face == "inner" else "inner"
    with pytest.raises(_deck_error()):
        _radial("sphere", _HOLLOW, (0, 1), **{face: TranslationDeck(), other: 0})


@pytest.mark.parametrize(("coord", "breakpoints", "faces"), [("sphere", _HOLLOW, ("inner", "outer")), ("cylinder", _SOLID, ("outer",))], ids=["hollow-sphere", "solid-cylinder"])
def test_a_mirror_deck_on_a_radial_face_is_a_deck_error_and_its_tag_builds(coord, breakpoints, faces) -> None:
    """The radial mirror is a TAG bound to ``WallAmplitudes(1, 0)``, never a ``MirrorDeck``: the meaning change of Appendix B.

    A reflective declaration on a curved face used to be a law; it is now unspellable as a deck and spelled as a tag whose
    wall the specification states. Both legs: ``MirrorDeck`` refused, the tag built."""
    from orpheus.geometry import MirrorDeck

    _radial(coord, breakpoints, (0, 1), **{f: 0 for f in faces})
    with pytest.raises(_deck_error()):
        _radial(coord, breakpoints, (0, 1), **{f: MirrorDeck() for f in faces})


def test_an_unpaired_wrap_on_a_slab_is_a_deck_error_at_the_geometry() -> None:
    """A ``TranslationDeck`` on one slab face whose partner carries anything else is a ``DeckError``: a wrap identifies two faces.

    The two re-posed rows ``wrap_not_paired`` (wrap left, vacuum right) and ``wrap_not_paired_right`` (reflective left, wrap
    right) of ``test_characteristic_walls.py``, plus the tag partner. The positive legs: both faces wrapped builds, and a
    mirror next to a tag builds."""
    from orpheus.geometry import MirrorDeck, StructuredGeometry, TranslationDeck

    bp, mats = (0.0, 1.0, 2.0), (0, 1)
    StructuredGeometry.slab(bp, mats, left=TranslationDeck(), right=TranslationDeck())
    StructuredGeometry.slab(bp, mats, left=MirrorDeck(), right=0)
    for left, right in ((TranslationDeck(), 0), (MirrorDeck(), TranslationDeck()), (0, TranslationDeck())):
        with pytest.raises(_deck_error()):
            StructuredGeometry.slab(bp, mats, left=left, right=right)


def test_walls_of_never_meets_a_curved_wrap() -> None:
    """The refusal precedes the reference: no ``Walls.of`` call can be made on a curved wrap, because there is no geometry to pass.

    Stated as the row's own sequence: constructing the geometry raises, so the reference's door is never reached with it. The
    live control: the reference reads the same body with a tag. Red at base: ``Walls.of`` takes one argument (the old
    signature raises ``TypeError`` for the two-argument call this row makes on the legal body)."""
    from orpheus.derivations.continuous.characteristic import Walls
    from orpheus.geometry import TranslationDeck
    from orpheus.specification.exterior import WallAmplitudes

    legal = _radial("sphere", _HOLLOW, (0, 1), inner=0, outer=0)
    walls = Walls.of(legal, {0: WallAmplitudes(0.0, 0.0)})
    require(len(walls.walls) == 2, f"the legal hollow body reads {len(walls.walls)} walls")
    with pytest.raises(_deck_error()):
        Walls.of(_radial("sphere", _HOLLOW, (0, 1), inner=TranslationDeck(), outer=TranslationDeck()), {})


@pytest.mark.rests_on(_HERE + "test_a_wrap_on_a_radial_face_is_a_deck_error_at_the_geometry")
def test_a_directly_built_walls_keeps_its_structural_guards_on_the_new_wall() -> None:
    """``Walls`` built by hand still refuses an unpaired wrap and a wrap on a radial chart (``NotImplementedError``, fragments kept).

    The re-posed ``test_a_directly_built_walls_refuses_what_walls_of_refuses``: the same five constructions on the new
    ``Wall(breakpoint, WallAmplitudes, partner)`` signature, with its last line (the out-of-bounds amplitude) moved to
    ``test_walls_of_reads_decks_and_amplitudes.py`` (R-8b). The guards are live (a ``Walls`` is public and can be built
    without the geometry), so NOT DEAD. Mutation: drop the ``acts_on_kept_space`` refusal or the partner check from
    ``Walls.__post_init__``: the radial and the unpaired legs redden. Red at base: ``Wall`` takes four floats."""
    from orpheus.derivations.continuous.characteristic import Wall, Walls
    from orpheus.geometry.chart import Chart
    from orpheus.geometry.coord import CoordSystem
    from orpheus.specification.exterior import WallAmplitudes

    mirror = WallAmplitudes(1.0, 0.0)
    vacuum = WallAmplitudes(0.0, 0.0)
    slab, sphere = Chart(CoordSystem.CARTESIAN), Chart(CoordSystem.SPHERICAL)
    with pytest.raises(ValueError, match="walls sit at distinct breakpoints"):
        Walls((Wall(1, vacuum, 1),), 3, slab)
    with pytest.raises(ValueError, match="walls sit at distinct breakpoints"):
        Walls((Wall(0, vacuum, 0), Wall(0, mirror, 0)), 3, slab)
    with pytest.raises(NotImplementedError, match="whose partner does not wrap back"):
        Walls((Wall(0, mirror, 3), Wall(3, mirror, 3)), 3, slab)
    with pytest.raises(NotImplementedError, match="periodic wrap on a cylinder or sphere"):
        Walls((Wall(0, mirror, 2), Wall(2, mirror, 0)), 2, sphere)
    wrapped = Walls((Wall(0, mirror, 3), Wall(3, mirror, 0)), 3, slab)
    require(wrapped.partner_at(np.array([0, 3])).tolist() == [3, 0], "a paired wrap on a slab is refused")
