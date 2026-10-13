"""Gates for the image's exact level (#590): a line's half-chords formed from a radius and its exact half-chord there.

P1 step (b), rung 5b widened (``.claude/plans/characteristic_reference_architecture.md``, "rung 5b: the reviews, two
rulings, and the sketch of the widened rung", ruled 2026-10-08, K): a line stores its impact parameter b through its
moment, to an ulp, so near a tangency, where the half-chord y is below sqrt(2 r eps(r)), b rounds the half-chord away
(qa's F3: the wall reading drifts as 3e-15/(1 - a) at every resolution). ``RadialImage`` carries a level
(r*, y*), by default (b, 0), and every half-chord is ``half_chord_at(r) = sqrt((r - r*)(r + r*) + y*^2)``. Spec:
``scratch/characteristic_architecture/p1_step_b5b/spec.md`` §8, rows K1-K6.

Every reference here is written in the test: HEAD's (``a5492017``) half-chord expression, copied, and the
closed-form crossing parameters of a canonical line. First reds: battery arms in
``scratch/characteristic_architecture/p1_step_b5b/ta/battery/``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem
from orpheus.geometry.chart import AxialImage, Chart, RadialImage
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.line import Line

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")

_RADII = (0.0, 0.5, 1.5, 2.0)
_SLAB = (0.0, 0.4, 1.5, 2.3)


def _partition(coord: CoordSystem, breakpoints=_RADII) -> ConcentricPartition:
    return ConcentricPartition(Chart(coord), tuple(breakpoints))


def _lines(coord: CoordSystem, seed: int = 7, n: int = 400) -> Line:
    """Seeded lines through the body, plus the hard members: b exactly on each breakpoint, b = 0, and (on the
    cylinder) lines along the axis."""
    rng = np.random.default_rng(seed)
    direction = rng.normal(size=(n, 3))
    direction /= np.linalg.norm(direction, axis=1)[:, None]
    point = rng.uniform(-2.2, 2.2, size=(n, 3))
    lines = [Line.through(point, direction)]
    for b in (*_RADII, np.nextafter(0.5, 0.0), np.nextafter(2.0, 0.0)):
        lines.append(Line.through(np.array([[0.0, b, 0.0]]), np.array([[1.0, 0.0, 0.0]])))
    if coord is CoordSystem.CYLINDRICAL:
        lines.append(Line.through(np.array([[0.7, 0.2, 0.0]]), np.array([[0.0, 0.0, 1.0]])))
    return Line(direction=np.concatenate([l.direction for l in lines]), moment=np.concatenate([l.moment for l in lines]))


def _head_half_chord(image: RadialImage, r: np.ndarray) -> np.ndarray:
    """HEAD ``a5492017``'s ``_radial_chord`` half-chord, copied: sqrt((r - b)(r + b)) where b < r and not parallel."""
    b = image.impact_parameter
    crossed = (b[..., None] < r) & ~image.parallel[..., None]
    return np.sqrt(np.where(crossed, (r - b[..., None]) * (r + b[..., None]), 0.0))


@pytest.mark.foundation
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL], ids=["sphere", "cylinder"])
def test_the_default_level_is_heads_arithmetic_bit_for_bit(coord) -> None:
    """[K2; foundation] ``Chart.image``'s default level (r* = b, y* = 0) forms every half-chord bit for bit as HEAD
    did, and the chord's crossing parameters are ``parameter_at(+-h)`` of those half-chords, ``array_equal``, on
    seeded lines and the hard members (b on every breakpoint, an ulp below the interface and the wall, b = 0, lines
    along the cylinder's axis).

    Every other caller of ``Chart.image`` is thereby unchanged. The hard members carry b EXACTLY on 0, 0.5, 1.5 and 2.0
    (asserted), where a tangency is decided. First reds: the default level y* = 1e-160 in place of 0; a tangent line
    b = r_k given a crossing (``_radial_chord``'s ``h > 0``; the same predicate inside ``half_chord_at`` is
    designed-green, since an uncrossed and a tangent radius both have h = 0).
    """
    partition = _partition(coord)
    r = np.asarray(_RADII)
    lines = _lines(coord)
    # HEAD's chord, replicated: the canonical image, its crossings, then the caller's shift
    canonical = lines.moved_by(partition.pose.inverse())
    image = partition.chart.image(canonical)
    assert isinstance(image, RadialImage)
    assert np.isin(np.asarray(_RADII), image.impact_parameter).all(), "the premise: b exactly on every breakpoint"
    head = _head_half_chord(image, r)
    assert np.array_equal(image.half_chord_at(r), head)
    shift = lines.parameter_of(partition.pose.on_points(canonical.foot))
    expected = image.parameter_at(np.concatenate([-head[..., ::-1], head], axis=-1)) + shift[..., None]
    crossed = head > 0.0
    present = np.concatenate([crossed[..., ::-1], crossed], axis=-1)
    chord = partition.chord(lines)
    assert np.array_equal(chord.crossings.present, present)
    assert np.array_equal(chord.crossings.parameter, expected)


@pytest.mark.foundation
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL], ids=["sphere", "cylinder"])
@pytest.mark.catches("ERR-106")
def test_the_half_chord_at_the_level_is_the_level_half_chord_exactly(coord) -> None:
    """[K1; foundation] ``half_chord_at(level) == level_half_chord`` exactly, for levels at every breakpoint and
    half-chords from 1e-300 to 1: sqrt(fl(y^2)) is y for every y whose square neither underflows nor overflows (a
    correctly rounded square and root), and (r* - r*) = 0. Also through ``chord(line, level=...)``: the image the
    chord carries holds the level. First red: the level dropped by ``chord`` (it forms from b)."""
    partition = _partition(coord)
    rng = np.random.default_rng(3)
    radius = rng.choice(np.asarray(_RADII[1:]), size=64)
    half_chord = 10.0 ** rng.uniform(-150, 0, size=64)
    b = np.sqrt(np.maximum((radius - half_chord) * (radius + half_chord), 0.0))
    lines = Line.through(np.stack([np.zeros(64), b, np.zeros(64)], axis=-1), np.tile([1.0, 0.0, 0.0], (64, 1)))
    chord = partition.chord(lines, level=(radius, half_chord))
    image = chord.image
    assert isinstance(image, RadialImage)
    assert image.level is not None and image.level_half_chord is not None
    assert np.array_equal(image.level, radius) and np.array_equal(image.level_half_chord, half_chord)
    assert np.array_equal(image.half_chord_at(radius[:, None])[:, 0], half_chord)


@pytest.mark.foundation
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL], ids=["sphere", "cylinder"])
@pytest.mark.catches("ERR-106")
def test_a_line_whose_rounded_b_is_the_level_still_crosses_it(coord) -> None:
    """[K3; foundation] A line grazing the wall with half-chord y in {1e-9, 1e-12, 1e-15}: its stored b loses y (the
    default image's half-chord at R is off by orders of magnitude, or absent where b rounds onto R: the premise,
    asserted); passed with its level (R, y), ``chord`` crosses R at closest +- y / |P Omega| (relative 1e-15 of y).
    One line has b = R exactly (no crossing without the level: HEAD's tangent line). First red: the level not
    threaded into ``_radial_chord``."""
    partition = _partition(coord)
    R = _RADII[-1]
    y = np.array([1e-9, 1e-12, 1e-15, 1e-9])
    b = np.sqrt((R - y) * (R + y))
    b[-1] = R
    lines = Line.through(np.stack([np.zeros(4), b, np.zeros(4)], axis=-1), np.tile([1.0, 0.0, 0.0], (4, 1)))
    blind = partition.chord(lines).image
    assert isinstance(blind, RadialImage)
    lost = blind.half_chord_at(np.array([R]))[:, 0]
    assert lost[-1] == 0.0 and np.all(np.abs(lost[:-1] / y[:-1] - 1.0) > 1e-3), ("the premise", lost)
    chord = partition.chord(lines, level=(np.full(4, R), y))
    image = chord.image
    assert isinstance(image, RadialImage)
    outer = chord.crossings.present & (chord.crossings.breakpoint == len(_RADII) - 1)
    assert np.all(outer.sum(axis=-1) == 2)
    t = np.where(outer, chord.crossings.parameter, np.nan)
    closest = image.closest_approach
    assert np.all(np.abs((np.nanmax(t, axis=-1) - closest) / y - 1.0) < 1e-6)   # the difference of two O(1) parameters
    assert np.all(np.abs((closest - np.nanmin(t, axis=-1)) / y - 1.0) < 1e-6)
    assert np.array_equal(image.half_chord_at(np.array([R]))[:, 0], y)


@pytest.mark.foundation
def test_the_axial_images_parameters_at_its_levels_are_its_chords_crossings() -> None:
    """[K4; foundation] On the slab, ``AxialImage.parameters_at(r)`` at the breakpoints equals the chord's crossing
    parameters bitwise (both directions, seeded lines): one definition of the parameter at a level on both image
    classes. First red: the chord's crossing formed by a second expression ((r - c)/rate re-spelt as r/rate - c/rate)."""
    partition = _partition(CoordSystem.CARTESIAN, _SLAB)
    lines = _lines(CoordSystem.CARTESIAN)
    keep = lines.direction[:, 0] != 0.0
    lines = Line(direction=lines.direction[keep], moment=lines.moment[keep])
    chord = partition.chord(lines)
    canonical = lines.moved_by(partition.pose.inverse())
    image = partition.chart.image(canonical)
    assert isinstance(image, AxialImage)
    shift = lines.parameter_of(partition.pose.on_points(canonical.foot))
    expected = image.parameters_at(np.asarray(_SLAB)[chord.crossings.breakpoint]) + shift[..., None]
    assert np.array_equal(chord.crossings.parameter, expected)


@pytest.mark.foundation
def test_a_level_passed_for_an_axial_image_is_refused() -> None:
    """[K5; foundation] ``chord(line, level=...)`` on the slab, whose image has no impact parameter and no
    half-chord, is refused naming the level. First red: the level silently ignored."""
    partition = _partition(CoordSystem.CARTESIAN, _SLAB)
    lines = Line.through(np.array([[0.7, 0.0, 0.0]]), np.array([[0.6, 0.8, 0.0]]))
    with pytest.raises(ValueError, match="a slab's lines have neither"):
        partition.chord(lines, level=(np.array([2.3]), np.array([0.1])))


@pytest.mark.foundation
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL], ids=["sphere", "cylinder"])
def test_the_parameters_at_a_level_are_the_chords_crossings(coord) -> None:
    """[K6; foundation] ``RadialImage.parameters_at(r_k, side)`` equals the chord's crossing at r_k on that side,
    bitwise, with and without a level: the reading's parameters and the chord's crossings are one expression (the
    elegance review's V2). First red: the reading's parameter re-spelt (``point_parameters``' old sqrt((c - b)(c + b)))."""
    partition = _partition(coord)
    R = _RADII[-1]
    y = np.array([0.3, 1e-6, 1e-12])
    b = np.sqrt((R - y) * (R + y))
    lines = Line.through(np.stack([np.zeros(3), b, np.zeros(3)], axis=-1), np.tile([1.0, 0.0, 0.0], (3, 1)))
    for level in (None, (np.full(3, R), y)):
        chord = partition.chord(lines, level=level)
        present = chord.crossings.present & (chord.crossings.breakpoint == len(_RADII) - 1)
        if not present.any():
            continue                                            # the default image of a grazing line: K3's premise
        image = chord.image
        assert isinstance(image, RadialImage)
        both = image.parameters_at(np.array([R]), np.array([-1.0]))[..., 0], image.parameters_at(np.array([R]), np.array([1.0]))[..., 0]
        rows = present.any(axis=-1)
        t = np.where(present, chord.crossings.parameter, np.nan)[rows]
        assert np.array_equal(np.nanmin(t, axis=-1), both[0][rows])
        assert np.array_equal(np.nanmax(t, axis=-1), both[1][rows])


@pytest.mark.foundation
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL], ids=["sphere", "cylinder"])
def test_a_level_that_disagrees_with_the_impact_parameter_is_refused(coord) -> None:
    """[K7; foundation] ``at_level`` holds the pair and the line's b to one line: |r*^2 - y*^2 - b^2| within 64 ulp of
    r*^2. Positive leg: the consistent pairs of K3 (y from 0.3 to 1e-15, b formed from them) are accepted. Negative
    leg: the same lines with y* moved by 1e-12 relative at y = 0.3 (r*^2 - y*^2 moves by 1.8e-13, above the 5.7e-14
    band at r* = 2), and with r* moved to the next radius, are refused naming the disagreement. First red: the check
    deleted (the wrong level is accepted and its half-chords read)."""
    chart = Chart(coord)
    R = _RADII[-1]
    y = np.array([0.3, 1e-6, 1e-12, 1e-15])
    b = np.sqrt((R - y) * (R + y))
    lines = Line.through(np.stack([np.zeros(4), b, np.zeros(4)], axis=-1), np.tile([1.0, 0.0, 0.0], (4, 1)))
    image = chart.image(lines)
    assert isinstance(image, RadialImage)
    image.at_level(np.full(4, R), y)                                       # accepted
    for radius, half_chord in ((np.full(4, R), y * np.array([1.0 + 1e-12, 1.0, 1.0, 1.0])),
                               (np.array([1.5, R, R, R]), y)):
        with pytest.raises(ValueError, match="level disagrees with its impact parameter"):
            image.at_level(radius, half_chord)


# ── KG1, KG2: the lines read from a crossing, and advanced along themselves (P3 step 3f, route B) ──────────────────
#
# Route B integrates each piece by the Gauss rule of its attenuation weight, its nodes held as OFFSETS from the piece's
# attenuated end; the source is read at those offsets through ``Chord.rays_at`` (a crossing, anchored exactly) and
# ``Rays.advanced`` (a slot end, advanced from its crossing) (the user's rulings of 2026-10-13). Spec:
# ``scratch/characteristic_architecture/p3/product/gate_spec.md`` (R3, KG1; KG2). The reference is written here in
# mpmath (80 digits): the line's point P and unit direction Omega promoted exactly, the crossing's parameter from the
# closed form (t* + side sqrt(r^2 - b^2)/|P Omega| on the cylinder and the sphere, (r - P_x)/Omega_x on the slab) in
# the reference's own parameterisation, and every coordinate |P_kept + t Omega_kept|; only DISTANCES along the line
# are compared, so the line's own parameter origin never enters (X4, input axis). First reds: in-process arms of
# ``scratch/characteristic_architecture/p3/product/ta/kg_arms.py`` (``KG_ARM=<arm> pytest -p kg_arms``), named per row.

import mpmath as mp  # noqa: E402

from orpheus.geometry.chord import Rays  # noqa: E402

_KG_EPS = float(np.finfo(float).eps)
#: KG1, KG2's band: 4 eps relative. `[M]` 2026-10-13 (the coordinator's probe): 1.00 eps worst on KG1's rows, 0.83 eps
#: on an advanced point; a half-chord is two products, a sum and a root of correctly rounded operations.
_KG_BAND = 4


def _kg_chart(chart: str) -> CoordSystem:
    return {"sphere": CoordSystem.SPHERICAL, "cylinder": CoordSystem.CYLINDRICAL, "slab": CoordSystem.CARTESIAN}[chart]


def _kg_chord(chart: str, point, direction, level=None):
    partition = _partition(_kg_chart(chart), _SLAB if chart == "slab" else _RADII)
    line = Line.through(np.array([point], dtype=float), np.array([direction], dtype=float))
    lv = None if level is None else (np.array([level[0]]), np.array([level[1]]))
    return partition.chord(line, lv)


def _kg_column(chord, radius: float, sense: int) -> np.ndarray:
    """The crossing column of ``radius`` with ``sense`` on the (one) line, found from the breakpoint VALUES."""
    r = np.asarray(_SLAB if isinstance(chord.image, AxialImage) else _RADII)[chord.crossings.breakpoint[0]]
    hit = np.flatnonzero(chord.crossings.present[0] & (r == radius) & (chord.crossings.sense[0] == sense))
    assert hit.size == 1, (radius, sense, hit)
    return hit


def _kg_reference(chart: str, point, direction, radius: float, sense: int, distance: float) -> float:
    """The orbit coordinate at ``distance`` along the line from its crossing of ``radius`` with ``sense``, 80 digits."""
    with mp.workdps(80):
        P = [mp.mpf(repr(float(v))) for v in point]
        O = [mp.mpf(repr(float(v))) for v in direction]
        d = mp.mpf(repr(float(distance)))
        if chart == "slab":
            return float(mp.mpf(repr(radius)) + O[0] * d)
        keep = 3 if chart == "sphere" else 2
        v2 = sum(O[i] ** 2 for i in range(keep))
        t_star = -sum(P[i] * O[i] for i in range(keep)) / v2
        b2 = sum(P[i] ** 2 for i in range(keep)) - sum(P[i] * O[i] for i in range(keep)) ** 2 / v2
        t_c = t_star + sense * mp.sqrt(mp.mpf(repr(radius)) ** 2 - b2) / mp.sqrt(v2)
        return float(mp.sqrt(sum((P[i] + (t_c + d) * O[i]) ** 2 for i in range(keep))))


def _kg_close(got: float, want: float, what) -> None:
    assert abs(got - want) <= _KG_BAND * _KG_EPS * abs(want), (what, got, want, abs(got - want) / abs(want))


def _unit(x: float, z: float = 0.0) -> tuple[float, float, float]:
    return x, float(np.sqrt((1.0 - x) * (1.0 + x) - z * z)), z


#: (id, chart, the line's point (its parameter origin elsewhere than the crossing), direction, crossing radius, sense,
#: the level passed to the chord or None)
_KG1_ROWS = [
    ("slab-face-exit", "slab", (0.7, 0.0, 0.0), (-0.6, 0.8, 0.0), 0.0, -1, None),
    ("slab-face-grazing", "slab", (1.9, 0.0, 0.0), _unit(-0.05), 0.0, -1, None),
    ("slab-interface", "slab", (0.1, 0.0, 0.0), _unit(0.3), 0.4, 1, None),
    ("sphere-wall-exit-b0.7", "sphere", (0.7, -0.4, 0.0), (0.0, 1.0, 0.0), 2.0, 1, None),
    ("sphere-interface-inbound", "sphere", (0.3, -1.9, 0.0), (0.0, 1.0, 0.0), 0.5, -1, None),
    ("cylinder-interface-3d", "cylinder", (0.3, -0.9, 0.4), (0.0, 0.6, 0.8), 1.5, 1, None),
    # grazing the wall at half-chord 2e-12, the chord minted with its exact level (as ``rays_from`` mints it)
    ("sphere-wall-grazing-level", "sphere", (2.0, 0.0, 0.0), _unit(-1e-12), 2.0, -1, (2.0, 2.0 * 1e-12)),
]
_KG1_DISTANCES = (0.0, 1e-15, -1e-15, 1e-12, -1e-12, 1e-9, -1e-6, 1e-3, -0.1)


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "point", "direction", "radius", "sense", "level"), [r[1:] for r in _KG1_ROWS],
                         ids=[r[0] for r in _KG1_ROWS])
def test_a_line_read_from_its_crossing_forms_every_distance_from_the_surface(chart, point, direction, radius, sense,
                                                                           level) -> None:
    """[KG1; foundation, REFERENCE] ``chord.rays_at(k)`` reads its line from crossing k: its coordinate is exactly
    r_k, ``orbit_coordinate_at(0) == r_k``, and ``orbit_coordinate_at(l)`` meets the mpmath coordinate at distance l
    from the crossing within 4 eps relative for l from -0.1 to 1e-3 (both signs: the piece behind an exit, the piece
    ahead of an entry). On chords whose parameter origin is not the crossing. First reds `[M]` 2026-10-13
    (``kg_arms.py``, ``battery/runs/kg_<arm>.log``): ``param_route`` (the coordinate through the line parameter) the 2
    slab-face rows (8e13 to 4.5e15 eps: the parameter's rounding read as a coordinate); ``no_relevel`` (the image
    left at the line's own level) the 3 radial rows whose chord carries no level; ``b_relevel`` (re-levelled from b,
    sqrt((r - b)(r + b))) the grazing row. DECLARED blind: the radial rows to ``param_route`` (c ~ r_k, the round trip
    costs a relative ulp); the grazing row to ``no_relevel`` (its chord is minted with the exact level already); the
    slab rows to both re-level arms (an axial image has no level)."""
    chord = _kg_chord(chart, point, direction, level)
    k = _kg_column(chord, radius, sense)
    rays = chord.rays_at(k)
    assert isinstance(rays, Rays)
    assert float(rays.orbit_coordinate[0]) == radius
    ells = np.array([_KG1_DISTANCES])
    got = np.asarray(rays.orbit_coordinate_at(ells), dtype=float)[0]
    assert got[0] == radius, (got[0], radius)
    for ell, g in zip(_KG1_DISTANCES[1:], got[1:]):
        _kg_close(float(g), _kg_reference(chart, point, direction, radius, sense, ell), ell)


#: (id, chart, point, direction, crossing radius, sense, advance d, reads e after it)
_KG2_ROWS = [
    ("sphere-exit-behind", "sphere", (0.7, -0.4, 0.0), (0.0, 1.0, 0.0), 2.0, 1, -0.3, (0.0, -0.2, 1e-9, 0.15)),
    ("sphere-entry-ahead", "sphere", (0.7, -0.4, 0.0), (0.0, 1.0, 0.0), 2.0, -1, 0.25, (0.0, 1e-12, -0.2)),
    # past the closest approach: entry at 2 on the b = 0.7 line, half-chord 1.874; d = 2.5 is 0.626 beyond it
    ("sphere-entry-past-closest", "sphere", (0.7, -0.4, 0.0), (0.0, 1.0, 0.0), 2.0, -1, 2.5, (0.0, -0.7, 0.3)),
    ("sphere-exit-back-past-closest", "sphere", (0.7, -0.4, 0.0), (0.0, 1.0, 0.0), 2.0, 1, -3.0, (0.0, 0.2)),
    ("cylinder-3d-tiny-advance", "cylinder", (0.3, -0.9, 0.4), (0.0, 0.6, 0.8), 1.5, 1, -1e-9, (0.0, -1e-9, 0.1)),
    ("cylinder-3d-past-closest", "cylinder", (0.3, -0.9, 0.4), (0.0, 0.6, 0.8), 1.5, -1, 2.6, (0.0, -0.5)),
    ("slab-face-behind", "slab", (0.7, 0.0, 0.0), (-0.6, 0.8, 0.0), 0.0, -1, -0.5, (0.0, 1e-9, -0.2)),
]


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "point", "direction", "radius", "sense", "d", "reads"), [r[1:] for r in _KG2_ROWS],
                         ids=[r[0] for r in _KG2_ROWS])
@pytest.mark.rests_on("tests/gates/geometry/test_chord_level.py::"
                      "test_a_line_read_from_its_crossing_forms_every_distance_from_the_surface")
def test_a_line_advanced_from_its_crossing_reads_the_coordinates_beyond_it(chart, point, direction, radius, sense, d,
                                                                           reads) -> None:
    """[KG2; foundation, REFERENCE] ``rays_at(k).advanced(d)``: its coordinate is the mpmath coordinate at distance d
    from the crossing, and ``orbit_coordinate_at(e)`` from it the coordinate at d + e, within 4 eps of the position's
    scale r_k + |d| (+ |e|), absolute: the advance spends one rounding of the position s' = +-y + |P Omega| d, so near
    the closest approach, where |s'| is small against the parts it was formed from, the RELATIVE error of c is not
    bounded by eps (`[M]` 2026-10-13: 6.4 eps relative at c = 0.313 on ``cylinder-3d-past-closest``; the docstring of
    ``advanced`` states the rounding); d of both
    signs, tiny (1e-9) and past the closest approach (the side flips: asserted). First reds `[M]` 2026-10-13
    (``kg_arms.py``): ``adv_no_relevel`` (the new point read with the crossing's level) 6 of 6 radial rows;
    ``adv_no_flip`` (the side kept) the 3 rows past the closest approach; ``no_relevel`` 6 of 6 radial rows (the
    crossing's level is where the advance starts). DECLARED blind: the slab row to every arm but none (its advanced
    point sits at 0.3, away from c = 0, so ``param_route`` costs an ulp there; KG1's slab rows catch it)."""
    chord = _kg_chord(chart, point, direction)
    k = _kg_column(chord, radius, sense)
    moved = chord.rays_at(k).advanced(np.array([d]))
    scale = radius + abs(d)                                    # the advanced position's scale, |s'| <= r_k + |d|
    want = _kg_reference(chart, point, direction, radius, sense, d)
    assert abs(float(moved.orbit_coordinate[0]) - want) <= _KG_BAND * _KG_EPS * scale, (
        "the point", (float(moved.orbit_coordinate[0]) - want) / (_KG_EPS * scale))
    got = np.asarray(moved.orbit_coordinate_at(np.array([reads])), dtype=float)[0]
    for e, g in zip(reads, got):
        want = _kg_reference(chart, point, direction, radius, sense, d + e)
        assert abs(float(g) - want) <= _KG_BAND * _KG_EPS * (scale + abs(e)), (e, (float(g) - want) / (_KG_EPS * scale))
    if chart != "slab":
        flipped = abs(d) > 1.0 and np.sign(d) != sense
        assert (float(moved.side[0]) == -sense) if flipped else (float(moved.side[0]) == sense), (moved.side, sense, d)


@pytest.mark.foundation
@pytest.mark.parametrize("row", [r[0] for r in _KG2_ROWS])
@pytest.mark.rests_on("tests/gates/geometry/test_chord_level.py::"
                      "test_a_line_advanced_from_its_crossing_reads_the_coordinates_beyond_it")
def test_advancing_and_returning_lands_on_the_crossing(row: str) -> None:
    """[KG2b; foundation, THEOREM] ``rays_at(k).advanced(d).advanced(-d)`` reads the crossing's coordinate r_k within
    _KG_RETURN eps of max(r_k, |P Omega| |d|) (each advance spends one rounding of the position, of size eps times
    |s'| <= r_k + |P Omega| |d|). First reds `[M]`: ``adv_no_relevel`` 6 of 7, ``adv_no_flip`` the 3 rows past the
    closest approach, ``no_relevel`` 6 of 7 (the slab row DECLARED blind to all three)."""
    _, chart, point, direction, radius, sense, d, _ = next(r for r in _KG2_ROWS if r[0] == row)
    chord = _kg_chord(chart, point, direction)
    k = _kg_column(chord, radius, sense)
    back = chord.rays_at(k).advanced(np.array([d])).advanced(np.array([-d]))
    scale = max(radius, abs(d))
    assert abs(float(back.orbit_coordinate[0]) - radius) <= _KG_RETURN * _KG_EPS * scale, (
        float(back.orbit_coordinate[0]) - radius) / (_KG_EPS * scale)
    assert float(np.asarray(back.orbit_coordinate_at(np.array([[0.0]])))[0, 0]) == float(back.orbit_coordinate[0])


#: KG2b's band, in eps of max(r_k, |d|). `[M]` 2026-10-13: <= 0.67 on the seven rows (KG2's on the advanced point,
#: <= 0.49 of r_k + |d|).
_KG_RETURN = 4


@pytest.mark.foundation
def test_a_line_is_not_read_from_a_crossing_it_does_not_make() -> None:
    """[KG3; foundation] A sphere line at b = 0.7 never meets r = 0.5: ``rays_at`` on that column is refused, naming
    the reason; a column it does make is not (the positive leg). First red `[M]`: the presence check removed
    (``no_refusal``), 1 of 1."""
    chord = _kg_chord("sphere", (0.7, -0.4, 0.0), (0.0, 1.0, 0.0))
    r = np.asarray(_RADII)[chord.crossings.breakpoint[0]]
    absent = np.flatnonzero(~chord.crossings.present[0] & (r == 0.5))
    assert absent.size >= 1
    with pytest.raises(ValueError, match="a line is read from a crossing it makes"):
        chord.rays_at(absent[:1])
    chord.rays_at(_kg_column(chord, 2.0, 1))
