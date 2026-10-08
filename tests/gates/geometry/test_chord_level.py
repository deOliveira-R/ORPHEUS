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
