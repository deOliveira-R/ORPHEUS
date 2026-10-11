"""Gates for :meth:`orpheus.geometry.chord.Chord.lengths_between`, the slots of a line cut to ``[start, stop]`` (#578, P2).

Verification spec: ``scratch/characteristic_architecture/p2/gate_spec.md`` §P.1b
(rows LB1-LB7). The verb is the kernel's one definition of a run along a line:
:meth:`~orpheus.geometry.chord.Chord.lengths_beyond` is ``lengths_between(start,
inf)``, so the five half-line rows of ``test_chord.py`` (C5, the posed
half-line, the hand-written multi-region table, the half-line from the cavity,
the mixed batch) pin this verb's ``stop = inf`` face with hand-written values
and stay as they are.

Every reference is built HERE from the fixture's construction data: the
crossings of a line with the circles are the roots of the line's kept-space
quadratic ``|P(x + s d)|^2 = r_k^2`` in mpmath at 40 dps (not the kernel's
level half-chords, X4), and a piece's region is read at its midpoint.

First reds (``scratch/characteristic_architecture/p2/ta/lb/battery``, each
arm a source transform of the verb installed on the class): see each row.
"""
from __future__ import annotations

import warnings

import mpmath as mp
import numpy as np
import pytest

from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.line import Line

_EPS = np.finfo(float).eps
_CHORD = "tests/gates/geometry/test_chord.py::"
_HERE = "tests/gates/geometry/test_chord_lengths_between.py::"

_SOLID = (0.0, 0.3, 1.1, 2.0)
_HOLLOW = (0.4, 1.1, 2.0)


def _partition(coord: CoordSystem, breakpoints) -> ConcentricPartition:
    return ConcentricPartition(Chart(coord), breakpoints)


def _kept(coord: CoordSystem) -> int:
    """The number of leading coordinates the orbit coordinate reads: 3 on the sphere, 2 on the cylinder."""
    return 3 if coord is CoordSystem.SPHERICAL else 2


def _region_lengths_mp(coord: CoordSystem, breakpoints, p, d, s0: float, s1: float) -> list:
    """The length of ``[s0, s1]`` (distances from ``p`` along the unit ``d``) in each region ``0 .. n-1`` and in the cavity, mpmath."""
    mp.mp.dps = 40
    k = _kept(coord)
    P = [mp.mpf(float(v)) for v in p]
    D = [mp.mpf(float(v)) for v in d]
    a = sum(D[j] ** 2 for j in range(k))
    b = 2 * sum(P[j] * D[j] for j in range(k))
    c = sum(P[j] ** 2 for j in range(k))
    lo, hi = mp.mpf(s0), mp.mpf(s1)
    cuts = [lo, hi]
    for r in breakpoints:
        disc = b * b - 4 * a * (c - mp.mpf(r) ** 2)
        if disc > 0:
            for s in ((-b - mp.sqrt(disc)) / (2 * a), (-b + mp.sqrt(disc)) / (2 * a)):
                if lo < s < hi:
                    cuts.append(s)
    cuts.sort()
    n = len(breakpoints) - 1
    out = [mp.mpf(0)] * (n + 1)                       # regions 0 .. n-1, then the cavity (code n)
    for u, v in zip(cuts[:-1], cuts[1:]):
        m = (u + v) / 2
        rho = mp.sqrt(c + b * m + a * m * m)
        if rho < breakpoints[0]:
            out[n] += v - u
        elif rho <= breakpoints[-1]:
            region = next(j for j in range(n) if rho <= breakpoints[j + 1])
            out[region] += v - u
    return out


def _random_lines(rng, n: int):
    p = np.stack([rng.uniform(-2.2, 2.2, n), rng.uniform(-2.2, 2.2, n), rng.uniform(-1.0, 1.0, n)], axis=-1)
    d = rng.normal(size=(n, 3))
    d /= np.linalg.norm(d, axis=1)[:, None]
    return p, d


# ── LB1: per region, against the roots ────────────────────────────────────


@pytest.mark.l0
@pytest.mark.verifies("geometry-chord-segment-lengths")
@pytest.mark.rests_on(_CHORD + "test_the_multi_region_segments_are_the_hand_written_table")
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL], ids=["sphere", "cylinder"])
@pytest.mark.parametrize("breakpoints", [_SOLID, _HOLLOW], ids=["solid", "hollow"])
def test_the_lengths_between_two_parameters_are_the_roots_pieces(coord: CoordSystem, breakpoints) -> None:
    """LB1: 200 random lines (seed 578) cut at random ``start < stop`` within 3 of the base point, summed per region and
    per cavity, against the mpmath pieces between the roots.

    Band per region ``16 eps (R + |start| + |stop|)`` on the kernel's parameters (`[M]` the prototype over 3000 lines:
    worst 53 eps absolute, parameters up to about 5); and no cut slot exceeds its traversed length. First reds: the cut
    formed from ``stop`` alone (``max(begin, start)`` replaced by ``begin``); ``min(end, stop)`` spelled ``max``.
    """
    rng = np.random.default_rng(578)
    part = _partition(coord, breakpoints)
    n = 200
    p, d = _random_lines(rng, n)
    line = Line.through(p, d)
    ch = part.chord(line)
    t = line.parameter_of(p)
    u, v = rng.uniform(-3.0, 3.0, (2, n))
    s0, s1 = np.minimum(u, v), np.maximum(u, v)
    got = ch.lengths_between(t + s0, t + s1)
    over = got > ch.traversed_length
    assert not np.any(over), f"{int(over.sum())} cut slots exceed their traversed length"
    R = breakpoints[-1]
    for i in range(n):
        want = _region_lengths_mp(coord, breakpoints, p[i], d[i], s0[i], s1[i])
        band = 16 * _EPS * (R + abs(t[i] + s0[i]) + abs(t[i] + s1[i]))
        for region, w in enumerate(want):
            g = float(np.sum(got[i][ch.slot_region[i] == region]))
            assert abs(g - float(w)) <= band, (
                f"line {i}, region {region}: {g!r} against {mp.nstr(w, 17)} (band {band:.1e})"
            )


# ── LB2: the half-line is the stop = +inf face ──────────────────────────────


@pytest.mark.foundation
@pytest.mark.rests_on(_CHORD + "test_a_half_line_keeps_the_slots_beyond_its_start",
                      _CHORD + "test_a_half_line_starting_in_the_cavity",
                      _CHORD + "test_lengths_beyond_on_a_batch_mixing_parallel_and_crossing_lines")
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL, CoordSystem.CARTESIAN],
                         ids=["sphere", "cylinder", "slab"])
def test_a_stop_at_infinity_is_the_half_line_bit_for_bit(coord: CoordSystem) -> None:
    """LB2: ``lengths_between(start, inf)`` keeps every slot wholly beyond ``start`` at its cancellation-free traversed
    length, bit for bit, and cuts only the slot holding ``start``, as ``end - start``.

    Asserted from the traversed lengths directly (not through ``lengths_beyond``, which is the same call): a slot whose
    ``begin >= start`` is ``array_equal`` to ``traversed_length``; the cut slot is ``min(end - start, traversed)``.
    First red: the wholly-inside branch dropped, every slot formed as ``min(end, stop) - max(begin, start)``
    (`[M]` the prototype: 0 of 4 cases bit-identical without the branch).
    """
    rng = np.random.default_rng(579)
    bp = (-0.7, 0.3, 1.1, 2.0) if coord is CoordSystem.CARTESIAN else _SOLID
    part = _partition(coord, bp)
    p, d = _random_lines(rng, 500)
    line = Line.through(p, d)
    ch = part.chord(line)
    t = line.parameter_of(p)
    got = ch.lengths_between(t, np.full(t.shape, np.inf))
    begin = ch.crossings.parameter[..., :-1]
    end = ch.crossings.parameter[..., 1:]
    whole = begin >= t[:, None]
    np.testing.assert_array_equal(got[whole], ch.traversed_length[whole])
    cut = np.minimum(np.clip(end - t[:, None], 0.0, None), ch.traversed_length)
    np.testing.assert_array_equal(got[~whole], cut[~whole])
    assert np.count_nonzero(whole & (ch.traversed_length > 0.0)) > 100, "the draw must hold traversed whole slots"


# ── LB3: interior cuts, one per slot kind ──────────────────────────────────


@pytest.mark.l0
@pytest.mark.verifies("geometry-chord-segment-lengths")
@pytest.mark.rests_on(_HERE + "test_the_lengths_between_two_parameters_are_the_roots_pieces")
def test_an_interior_cut_in_each_slot_kind() -> None:
    """LB3: hand algebra on the hollow sphere (0.4, 1.1, 2.0), the line ``x = 0.2`` along ``e_y`` (its parameter is y):
    slots inbound region 1 ``[-h_2, -h_1]``, region 0 ``[-h_1, -h_0]``, the cavity ``[-h_0, h_0]``, then outbound
    region 0 and region 1, with ``h_k = sqrt(r_k^2 - 0.04)``.

    (a) ``[-1.5, -1.2]`` inside the inbound region-1 slot; (b) ``[-0.2, 0.1]`` inside the cavity; (c) ``[1.2, 1.5]``
    inside the outbound region-1 slot; each 0.3 in its slot and 0 elsewhere. (d) the solid sphere's line ``x = 0.7``:
    region 1 holds its closest approach, its inbound and outbound slots split at ``t* = 0``; ``[-0.25, 0.4]`` gives
    0.25 inbound and 0.4 outbound. (e) start and stop on the kernel's own crossings ``-h_1`` and ``+h_0``: the slots
    between keep their traversed lengths bit for bit, the others are 0.
    First red: ``min(end, stop)`` spelled ``max(end, stop)``.
    """
    hollow = _partition(CoordSystem.SPHERICAL, _HOLLOW)
    line = Line.through(np.array([0.2, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]))
    ch = hollow.chord(line)
    np.testing.assert_array_equal(ch.slot_region, [1, 0, 2, 0, 1])
    for (s0, s1), slot in (((-1.5, -1.2), 0), ((-0.2, 0.1), 2), ((1.2, 1.5), 4)):
        got = ch.lengths_between(np.array(s0), np.array(s1))
        want = np.zeros(5)
        want[slot] = 0.3
        np.testing.assert_allclose(got, want, rtol=0, atol=4 * _EPS)
    solid = _partition(CoordSystem.SPHERICAL, _SOLID)
    turning = solid.chord(Line.through(np.array([0.7, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])))
    np.testing.assert_array_equal(turning.slot_region, [2, 1, 0, 3, 0, 1, 2])
    np.testing.assert_allclose(turning.lengths_between(np.array(-0.25), np.array(0.4)),
                               [0.0, 0.25, 0.0, 0.0, 0.0, 0.4, 0.0], rtol=0, atol=4 * _EPS)
    q = ch.crossings.parameter
    on = ch.lengths_between(q[1], q[3])                 # from the inbound crossing of r_1 to the outbound one of r_0
    np.testing.assert_array_equal(on, np.where(np.isin(np.arange(5), [1, 2]), ch.traversed_length, 0.0))


# ── LB4: parallel lines ─────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.rests_on(_CHORD + "test_a_line_parallel_to_the_orbit_space_is_inside_one_region_or_on_an_interface")
def test_a_parallel_line_holds_the_interval_in_its_home_slot() -> None:
    """LB4: an axial line of the cylinder at rho = 0.7 (region 1) and a slab line with ``Omega_x = 0`` at x = 0.5
    (region 1): the home slot holds ``stop - start`` = 3.5 exactly, every other slot 0; a cylinder line ON r = 1.1 has
    every slot 0 (assigned to neither side, as ``slot_length`` rules). A subnormal projected speed is no parallel line:
    its finite interval stays finite.
    First red: the home slot returned as ``slot_length`` (infinite) whatever the interval.
    """
    cyl = _partition(CoordSystem.CYLINDRICAL, _SOLID)
    slab = _partition(CoordSystem.CARTESIAN, (-0.7, 0.3, 1.1, 2.0))
    for part, point, direction in ((cyl, (0.7, 0.0, 2.0), (0.0, 0.0, 1.0)), (slab, (0.5, 0.0, 0.0), (0.0, 1.0, 0.0))):
        line = Line.through(np.array(point), np.array(direction))
        ch = part.chord(line)
        got = ch.lengths_between(np.array(-1.0), np.array(2.5))
        assert int(np.count_nonzero(got)) == 1 and float(np.sum(got)) == 3.5, got
        assert int(ch.slot_region[got > 0.0][0]) == 1
    on = cyl.chord(Line.through(np.array([1.1, 0.0, 0.0]), np.array([0.0, 0.0, 1.0])))
    np.testing.assert_array_equal(on.lengths_between(np.array(-1.0), np.array(2.5)), np.zeros(7))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        slow = cyl.chord(Line.through(np.array([0.7, 0.0, 0.0]), np.array([5e-320, 0.0, 1.0])))
        got = slow.lengths_between(np.array(-1.0), np.array(2.5))
    assert np.all(np.isfinite(got)) and float(np.sum(got)) == 3.5, got


# ── LB5: an empty interval ──────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_lengths_between_two_parameters_are_the_roots_pieces")
def test_an_empty_interval_has_no_length() -> None:
    """LB5: ``start == stop`` and ``start > stop`` give every slot exactly 0, on crossing and parallel lines alike.
    First red: the clip at 0 removed (negative lengths where ``stop < start``).
    """
    rng = np.random.default_rng(580)
    part = _partition(CoordSystem.SPHERICAL, _HOLLOW)
    p, d = _random_lines(rng, 100)
    line = Line.through(p, d)
    ch = part.chord(line)
    t = line.parameter_of(p)
    np.testing.assert_array_equal(ch.lengths_between(t, t), np.zeros(ch.slot_region.shape))
    np.testing.assert_array_equal(ch.lengths_between(t + 0.5, t - 0.5), np.zeros(ch.slot_region.shape))
    axial = _partition(CoordSystem.CYLINDRICAL, _SOLID).chord(Line.through(np.array([0.7, 0, 0.0]), np.array([0, 0, 1.0])))
    np.testing.assert_array_equal(axial.lengths_between(np.array(2.0), np.array(1.0)), np.zeros(7))


# ── LB6: additivity ─────────────────────────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_lengths_between_two_parameters_are_the_roots_pieces")
def test_two_adjacent_intervals_add_to_their_union() -> None:
    """LB6: ``between(a, b) + between(b, c) == between(a, c)`` for ``a < b < c``: bit for bit on every slot that does
    not hold ``b``, within ``8 eps (|a| + |c| + R)`` on the one that does.
    First red: an off-by-one in the cut mask (the slot after the one holding an end treated as cut).
    """
    rng = np.random.default_rng(581)
    part = _partition(CoordSystem.CYLINDRICAL, _HOLLOW)
    p, d = _random_lines(rng, 300)
    line = Line.through(p, d)
    ch = part.chord(line)
    t = line.parameter_of(p)
    a, m, c = np.sort(rng.uniform(-3.0, 3.0, (3, 300)), axis=0) + t
    left, right, whole = ch.lengths_between(a, m), ch.lengths_between(m, c), ch.lengths_between(a, c)
    begin, end = ch.crossings.parameter[..., :-1], ch.crossings.parameter[..., 1:]
    holds_m = (begin < m[:, None]) & (m[:, None] < end) & (ch.traversed_length > 0.0)
    np.testing.assert_array_equal((left + right)[~holds_m], whole[~holds_m])
    band = 8 * _EPS * (np.abs(a) + np.abs(c) + 2.0)[:, None]
    assert np.all(np.abs(left + right - whole) <= band)
    assert np.count_nonzero(holds_m) > 100, "the draw must cut slots in their interior"


# ── LB7: several distances per line, as the caller's batch shape ───────────


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_lengths_between_two_parameters_are_the_roots_pieces")
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL], ids=["sphere", "cylinder"])
def test_several_stops_per_line_broadcast_as_the_callers_batch(coord: CoordSystem) -> None:
    """LB7 (the coordinator's ruling, 2026-10-10): the verb adds no axis of its own; lines minted ``(n, 1)`` with
    ``stop`` ``(n, q)`` give ``(n, q, S)``, ``array_equal`` to the q calls on ``(n,)`` lines. A ``stop`` that does not
    broadcast against the lines' batch (``(n + 1,)`` against ``(n,)`` lines) raises.
    The q distances span every slot kind of the hollow body (inside, through the cavity, beyond the exit).
    First reds: the stop's trailing axis expanded twice (a shape error); the stop flattened before use (the refusal
    leg passes silently).
    """
    rng = np.random.default_rng(582)
    part = _partition(coord, _HOLLOW)
    n, q = 40, 6
    p, d = _random_lines(rng, n)
    batch = Line.through(p[:, None, :], d[:, None, :])
    ch = part.chord(batch)
    t = batch.parameter_of(p[:, None, :])
    dist = np.array([0.1, 0.6, 1.3, 2.2, 3.5, 6.0])
    got = ch.lengths_between(t, t + dist)
    assert got.shape == (n, q, ch.slot_region.shape[-1])
    flat = Line.through(p, d)
    one = part.chord(flat)
    t1 = flat.parameter_of(p)
    for j in range(q):
        np.testing.assert_array_equal(got[:, j, :], one.lengths_between(t1, t1 + dist[j]))
    with pytest.raises(ValueError, match="broadcast"):
        one.lengths_between(t1, np.zeros(n + 1))


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_lengths_between_two_parameters_are_the_roots_pieces")
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL], ids=["sphere", "cylinder"])
def test_a_cut_slot_one_rounding_short_of_whole_never_exceeds_its_length(coord: CoordSystem) -> None:
    """LB1b, the clip's witness (the elegance review's m2): with ``start`` one rounding inside a slot's opening crossing
    and ``stop`` on its closing one (and the mirror), the cut ``end - start`` formed by subtraction can exceed the
    slot's cancellation-free traversed length by a rounding; the verb clips it there. Asserted: no slot exceeds its
    traversed length, on a fixture where the unclipped cut does exceed it (counted, the activation premise).
    `[M]` the elegance review: the clip binds on 54 000 of 1 179 480 slots with ends one ulp apart, by <= 1.6e-15;
    uniform draws (LB1) miss it. First red: the upper clip removed.
    """
    rng = np.random.default_rng(583)
    part = _partition(coord, _HOLLOW)
    p, d = _random_lines(rng, 400)
    ch = part.chord(Line.through(p, d))
    begin, end = ch.crossings.parameter[..., :-1], ch.crossings.parameter[..., 1:]
    traversed = ch.traversed_length
    binding = 0
    for s in range(begin.shape[-1]):
        live = traversed[:, s] > 0.0
        for t0, t1 in ((np.nextafter(begin[:, s], np.inf), end[:, s]), (begin[:, s], np.nextafter(end[:, s], -np.inf))):
            got = ch.lengths_between(t0, t1)[:, s]
            assert np.all(got[live] <= traversed[live, s]), "a cut slot exceeds its traversed length"
            binding += int(np.count_nonzero((np.minimum(end[:, s], t1) - np.maximum(begin[:, s], t0) > traversed[:, s])
                                            & live))
    assert binding > 0, "the fixture never makes the unclipped cut exceed the traversed length: the row is vacuous"
