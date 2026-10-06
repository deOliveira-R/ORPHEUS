"""Gates for :attr:`orpheus.geometry.chord.Chord.transits`: the maximal runs of a line inside the domain.

Migration step (a) of P1 of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "P1 API sketch"
item 1 and the ruling "Transits: in the kernel"; verification spec
``scratch/characteristic_architecture/p1_verification_spec.md`` row A1, whose
hand-counted wall table this file realises on the kernel's slots).

The definition gated: a transit is a maximal run of slots in which every
TRAVERSED slot (``slot_length > 0``) has an interior region code (``< n``);
an untraversed slot never breaks a run; a transit holds at least one
traversed slot. Transits are ordered along the line; a wall is the
breakpoint index (0 or n) of the crossing that opens the first traversed
slot or closes the last one, never a position in the tuple. An absent
transit is the empty slot range at the end of the chord
(``first_slot == stop_slot ==`` the slot count ``S``) with the wall code
``n + 1`` (ruled 2026-10-06 by the main agent, answering the gate's
question: ``n + 1`` would be a real slot on a radial chart).

Every expected value is counted by hand from the slot layout in
:class:`~orpheus.geometry.chord.Chord`'s docstring: on the sphere and the
cylinder the slots are the regions inbound ``n-1 .. 0``, the cavity, the
regions outbound ``0 .. n-1``, and the crossing that opens slot ``i`` is
crossing ``i``, on breakpoints ``n, .., 0, 0, .., n``; on the slab the slots
are the regions in the order the line meets them.

Fixtures (as ``test_chord.py``): ``SPH`` (0, 0.3, 1.1, 2.0) solid,
``SPH_H`` (0.4, 1.1, 2.0) hollow, the cylinders alike, ``SLB``
(-0.7, 0.3, 1.1, 2.0).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.chart import Chart, RadialImage
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.line import Line

_HERE = "tests/gates/geometry/test_chord_transits.py::"
_CHORD = "tests/gates/geometry/test_chord.py::"

_SPH_BP = (0.0, 0.3, 1.1, 2.0)
_SPH_H_BP = (0.4, 1.1, 2.0)
_SLB_BP = (-0.7, 0.3, 1.1, 2.0)

SPH = ConcentricPartition(Chart(CoordSystem.SPHERICAL), _SPH_BP)
SPH_H = ConcentricPartition(Chart(CoordSystem.SPHERICAL), _SPH_H_BP)
CYL = ConcentricPartition(Chart(CoordSystem.CYLINDRICAL), _SPH_BP)
CYL_H = ConcentricPartition(Chart(CoordSystem.CYLINDRICAL), _SPH_H_BP)
SLB = ConcentricPartition(Chart(CoordSystem.CARTESIAN), _SLB_BP)

_BELOW_R0 = float(np.nextafter(0.4, 0.0))       # one ulp inside the cavity's radius
_BELOW_R = float(np.nextafter(2.0, 0.0))        # one ulp inside the outer surface

# Hand-counted transits, (first_slot, stop_slot, entry_wall, exit_wall) per transit,
# keyed by impact parameter b. Solid (n = 3, S = 7): one transit (0, 7, 3, 3) for
# every b < R, whatever the regions it misses (their slots are untraversed and
# never break the run; the cavity slot 3 is untraversed on a solid body). Hollow
# (n = 2, S = 5): b < r_0 traverses the cavity slot 2, so two transits, the first
# closing on the inner wall 0 and the second opening on it.
_ONE_SOLID = [(0, 7, 3, 3)]
_ONE_HOLLOW = [(0, 5, 2, 2)]
_TWO_HOLLOW = [(0, 2, 2, 0), (3, 5, 0, 2)]
_SOLID_ROWS = [  # (id, b, transits)
    ("b0_centre", 0.0, _ONE_SOLID),
    ("b0.2_region0", 0.2, _ONE_SOLID),
    ("b0.3_interior_tangency", 0.3, _ONE_SOLID),
    ("b0.7_region1", 0.7, _ONE_SOLID),
    ("b1.1_interior_tangency", 1.1, _ONE_SOLID),
    ("b1.5_outer_region_only", 1.5, _ONE_SOLID),
    ("b_one_ulp_below_R", _BELOW_R, _ONE_SOLID),
    ("bR_outer_tangency", 2.0, []),
    ("b2.5_misses", 2.5, []),
]
_HOLLOW_ROWS = [
    ("b0_centre", 0.0, _TWO_HOLLOW),
    ("b0.2_through_cavity", 0.2, _TWO_HOLLOW),
    ("b_one_ulp_below_r0", _BELOW_R0, _TWO_HOLLOW),
    ("b_r0_cavity_tangency", 0.4, _ONE_HOLLOW),
    ("b0.7_region0", 0.7, _ONE_HOLLOW),
    ("b1.1_interior_tangency", 1.1, _ONE_HOLLOW),
    ("b1.5_outer_region_only", 1.5, _ONE_HOLLOW),
    ("bR_outer_tangency", 2.0, []),
]
_AXIAL = [0.0, 0.8, -0.6]                         # the cylinder's Omega_z: in-plane, up, down


def _line(p, omega) -> Line:
    return Line.through(np.array(p, dtype=float), np.array(omega, dtype=float))


def _cases():
    """(id, partition, point, direction, expected) for every hand-counted row."""
    rows = []
    for name, part_s, part_c, table in (("solid", SPH, CYL, _SOLID_ROWS), ("hollow", SPH_H, CYL_H, _HOLLOW_ROWS)):
        for rid, b, expected in table:
            # The point (b, 0, 0) is the closest approach and the in-plane direction is along y, so the kernel's
            # impact parameter is b bit for bit (as test_chord.py's tangency rows), except one case measured:
            # a tilted cylinder line one ulp below R, whose Plucker foot rounds b up to R (`_expected_at`).
            rows.append((f"sphere_{name}_{rid}", part_s, (b, 0.0, 0.0), (0.0, 1.0, 0.0), expected))
            for wz in _AXIAL:
                s = float(np.sqrt((1.0 - wz) * (1.0 + wz)))
                rows.append((f"cylinder_{name}_{rid}_wz{wz}", part_c, (b, 0.0, 0.0), (0.0, s, wz), expected))
    # The slab: one transit from the wall the line meets first to the other.
    rows.append(("slab_rising", SLB, (0.5, 0.0, 0.0), (0.6, 0.8, 0.0), [(0, 3, 0, 3)]))
    rows.append(("slab_falling", SLB, (0.5, 0.0, 0.0), (-0.6, 0.8, 0.0), [(0, 3, 3, 0)]))
    rows.append(("slab_rising_from_outside", SLB, (-5.0, 2.0, 1.0), (0.6, 0.0, -0.8), [(0, 3, 0, 3)]))
    # A near-axial cylinder line: |P Omega| = 1e-200, every traversed length ~1e200 and finite, still one transit.
    rows.append(("cylinder_near_axial", CYL, (0.7, -0.4, 0.0), (0.0, 1e-200, 1.0), _ONE_SOLID))
    # Parallel lines (|P Omega| = 0): inside a region (an infinite slot), on an interface, outside. None has a transit.
    rows.append(("cylinder_parallel_in_region1", CYL, (0.7, 0.0, 0.0), (0.0, 0.0, 1.0), []))
    rows.append(("cylinder_parallel_on_axis", CYL, (0.0, 0.0, 0.0), (0.0, 0.0, -1.0), []))
    rows.append(("cylinder_parallel_in_interface", CYL, (1.1, 0.0, 0.0), (0.0, 0.0, 1.0), []))
    rows.append(("cylinder_parallel_in_cavity", CYL_H, (0.2, 0.0, 0.0), (0.0, 0.0, 1.0), []))
    rows.append(("cylinder_parallel_outside", CYL, (2.5, 0.0, 0.0), (0.0, 0.0, 1.0), []))
    rows.append(("slab_parallel_in_region1", SLB, (0.5, 0.0, 0.0), (0.0, 0.6, 0.8), []))
    rows.append(("slab_parallel_in_interface", SLB, (1.1, 0.0, 0.0), (0.0, 1.0, 0.0), []))
    rows.append(("slab_parallel_in_outer_wall", SLB, (2.0, 0.0, 0.0), (0.0, 0.0, 1.0), []))
    return rows


_CASES = _cases()


def _radial_b(image, index=()) -> float:
    """The impact parameter a radial chord's image carries (the slab's axial image has none)."""
    if not isinstance(image, RadialImage):
        pytest.fail("a curvilinear chord's image is radial")
    return float(image.impact_parameter[index])


def _hand_rule(part: ConcentricPartition, b: float) -> list:
    """The table's rule, stated once by hand: through the cavity two, through the body one, else none."""
    r0, rn = part.breakpoints[0], part.breakpoints[-1]
    if b >= rn:
        return []
    if r0 > 0.0 and b < r0:
        return _TWO_HOLLOW
    return _ONE_HOLLOW if r0 > 0.0 else _ONE_SOLID


def _expected_at(part: ConcentricPartition, kernel_b: float, b: float, expected: list) -> list:
    """The hand count at the chord's own impact parameter.

    The transits read the chord, so the stratum side is the chord's ``b``. Where
    the line's representation rounds ``b`` (``[M]`` 2026-10-06: a cylinder line
    tilted by ``Omega_z`` in {0.8, -0.6} through ``(R - ulp, 0, 0)`` reads
    ``b = R``), the count follows the rounded value by the same hand rule, and
    the rounding itself is bounded to one ulp (the chord's gates own it).
    """
    assert _hand_rule(part, b) == expected, "the table and its rule disagree"
    if kernel_b == b:
        return expected
    assert abs(kernel_b - b) <= np.spacing(b), f"the chord's b {kernel_b!r} is more than an ulp from {b!r}"
    return _hand_rule(part, kernel_b)


def _assert_transits(tr, expected, n: int, slots: int, where: str) -> None:
    """The transits at one line equal the hand table; absent columns carry the absent codes."""
    got = [
        (int(tr.first_slot[j]), int(tr.stop_slot[j]), int(tr.entry_wall[j]), int(tr.exit_wall[j]))
        for j in range(2) if bool(tr.present[j])
    ]
    assert got == expected, f"{where}: transits {got} != hand count {expected}"
    assert [bool(p) for p in tr.present] == [True] * len(expected) + [False] * (2 - len(expected)), (
        f"{where}: present transits are not the leading columns"
    )
    for j in range(len(expected), 2):
        assert (int(tr.first_slot[j]), int(tr.stop_slot[j])) == (slots, slots), f"{where}: absent slot range"
        assert (int(tr.entry_wall[j]), int(tr.exit_wall[j])) == (n + 1, n + 1), f"{where}: absent wall code"


@pytest.mark.foundation
@pytest.mark.rests_on(_CHORD + "test_the_crossings_carry_the_region_they_enter",
                      _CHORD + "test_a_tangency_is_not_a_crossing",
                      _CHORD + "test_the_centre_of_a_solid_body_is_never_crossed",
                      _CHORD + "test_the_slab_chord_in_both_orientations",
                      _CHORD + "test_a_line_parallel_to_the_orbit_space_is_inside_one_region_or_on_an_interface")
@pytest.mark.parametrize(("part", "p", "omega", "expected"), [c[1:] for c in _CASES], ids=[c[0] for c in _CASES])
def test_the_transits_match_the_hand_counted_walls(part: ConcentricPartition, p, omega, expected) -> None:
    """A1: every chart, every rank, both sides of every tangency stratum; walls by hand.

    Solid body: one transit from the outer wall to the outer wall for every
    ``b < R``, including ``b`` at an interior tangency (0.3, 1.1: the region
    inside is untraversed, which never breaks the run) and one ulp below ``R``;
    none at ``b = R`` (a tangency is not a crossing) and beyond. Hollow body:
    two transits for ``b < r_0`` and one ulp below it, the first closing on the
    inner wall ``0`` and the second opening on it; ONE at ``b = r_0`` exactly (the
    cavity slot is untraversed). The cylinder repeats every row at
    ``Omega_z`` in {0, 0.8, -0.6} (the in-plane impact parameter decides; the
    axial tilt only scales the lengths). The slab: ``(0, n)`` rising, ``(n, 0)``
    falling. Parallel lines, inside a region (an infinite slot), on the axis,
    in an interface, in the cavity or outside: none.
    First reds (run, ``[M]`` 2026-10-06, battery
    ``scratch/characteristic_architecture/p1_step_a/battery``, arms T1-T5):
    (a) an untraversed exterior slot counted as breaking the run
    (the solid body's untraversed cavity splits every ``b < r_1`` line into two
    transits; the hollow body at ``b = r_0`` reads two); (b) walls read by tuple
    position (``entry_wall = first_slot``: the sphere's 0 for 3); (c) a parallel
    line inside a region taken as a transit through its infinite slot; (d) the
    absent slot code ``n + 1`` (slot 4 of 7 on ``SPH``, a real slot).
    """
    ch = part.chord(_line(p, omega))
    if part.chart.acts_on_kept_space and not bool(ch.parallel):
        expected = _expected_at(part, _radial_b(ch.image), p[0], expected)
    _assert_transits(ch.transits, expected, part.n_regions, ch.slot_region.shape[-1], "single line")


def _batch(part: ConcentricPartition, rows) -> tuple[Line, list]:
    pts = np.array([r[1] for r in rows], dtype=float).reshape(2, -1, 3)
    oms = np.array([r[2] for r in rows], dtype=float).reshape(2, -1, 3)
    return Line.through(pts, oms), [r[3] for r in rows]


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_transits_match_the_hand_counted_walls")
@pytest.mark.parametrize("part", [CYL_H, CYL, SLB], ids=["cylinder_hollow", "cylinder_solid", "slab"])
def test_a_batch_of_mixed_lines_reads_each_line_as_it_reads_alone(part: ConcentricPartition) -> None:
    """A1, batch shape ``(2, k)``: every rank, parallel and interface lines mixed in one call.

    The batch is the hand table's rows on this partition, reshaped ``(2, k)``
    (so a broadcast of one line's answer over the batch, or a transposed
    leading axis, puts some row's transits on another's line); each element
    must equal its own hand count.
    First reds (run): every transit arm T1-T5 reds 2 or 3 of the 3 batches.
    Not run: a reduction over the wrong axis of the ``(..., 2, S)`` member mask.
    """
    rows = [(c[0], c[2], c[3], c[4]) for c in _CASES if c[1] is part]
    if len(rows) % 2:
        rows = rows[:-1]
    assert len(rows) >= 6, f"the batch must mix ranks: {len(rows)} rows"
    line, expected = _batch(part, rows)
    ch = part.chord(line)
    tr = ch.transits
    for field in ("first_slot", "stop_slot", "entry_wall", "exit_wall", "present"):
        assert getattr(tr, field).shape == (*line.shape, 2), f"{field} is (..., 2)"
    flat = np.ndindex(*line.shape)
    points = [r[1] for r in rows]
    for index, exp, p in zip(flat, expected, points, strict=True):
        if part.chart.acts_on_kept_space and not bool(ch.parallel[index]):
            exp = _expected_at(part, _radial_b(ch.image, index), p[0], exp)
        sub = type(tr)(*(getattr(tr, f)[index] for f in ("first_slot", "stop_slot", "entry_wall", "exit_wall", "present")))
        _assert_transits(sub, exp, part.n_regions, ch.slot_region.shape[-1], f"batch element {index}")
    ranks = {len(e) for e in expected}
    assert len(ranks) >= 2, f"the batch must mix ranks, got {ranks}"


def _runs_by_hand(slot_length: np.ndarray, slot_region: np.ndarray, breakpoint: np.ndarray, n: int, parallel: bool):
    """The definition, one slot at a time: no vectorisation, no cumulative sums."""
    if parallel:
        return []
    runs, current = [], None
    for i, (length, region) in enumerate(zip(slot_length, slot_region)):
        if not length > 0.0:
            continue                                    # untraversed: never breaks a run
        if region < n:
            current = [i, i] if current is None else [current[0], i]
        elif current is not None:
            runs.append(current)
            current = None
    if current is not None:
        runs.append(current)
    return [(a, b + 1, int(breakpoint[a]), int(breakpoint[b + 1])) for a, b in runs]


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_transits_match_the_hand_counted_walls")
@pytest.mark.parametrize("part", [SPH, SPH_H, CYL, CYL_H, SLB], ids=["sphere", "sphere_hollow", "cylinder", "cylinder_hollow", "slab"])
def test_the_transits_are_the_maximal_runs_of_traversed_interior_slots(part: ConcentricPartition) -> None:
    """A1, the law over 1500 seeded lines: the definition evaluated slot by slot.

    The population draws feet in ``[-2.5, 2.5]^3`` and isotropic directions,
    plus 1 in 5 lines with the axial or slab-parallel component zeroed (parallel
    lines) on the cylinder and the slab. Per line, the runs are computed here
    from ``slot_length > 0`` and ``slot_region < n`` by a plain loop; they must
    equal the present transits, every wall must be 0 or ``n`` (a wall is a
    boundary point, never an interface), the first and last slot of a transit
    are traversed, and two transits are ordered and disjoint.
    First reds: as the hand table, on whichever lines of the population carry
    them; the population reports how many lines had 0, 1 and 2 transits and
    requires each count non-zero where the partition admits it (X1: a
    population that never draws a hollow body's through-cavity line cannot see
    the run split).
    """
    rng = np.random.default_rng(20261006)
    feet = rng.uniform(-2.5, 2.5, (1500, 3))
    om = rng.normal(size=(1500, 3))
    if part.chart.kept_columns < 3:                     # a sphere has no parallel direction
        om[::5, : part.chart.kept_columns] = 0.0
    om /= np.linalg.norm(om, axis=1, keepdims=True)
    ch = part.chord(Line.through(feet, om))
    tr = ch.transits
    n = part.n_regions
    slots = ch.slot_region.shape[-1]
    bp = np.broadcast_to(ch.crossings.breakpoint, ch.crossings.parameter.shape)
    counts = {0: 0, 1: 0, 2: 0}
    for i in range(1500):
        expected = _runs_by_hand(ch.slot_length[i], ch.slot_region[i], bp[i], n, bool(ch.parallel[i]))
        sub = type(tr)(*(getattr(tr, f)[i] for f in ("first_slot", "stop_slot", "entry_wall", "exit_wall", "present")))
        _assert_transits(sub, expected, n, slots, f"line {i}")
        for _, _, entry, exit_ in expected:
            assert entry in (0, n) and exit_ in (0, n), f"line {i}: a wall at an interior breakpoint"
        counts[len(expected)] += 1
    assert counts[0] > 0 and counts[1] > 0, f"population {counts}"
    if part.chart.acts_on_kept_space and part.breakpoints[0] > 0.0:
        assert counts[2] > 0, f"no through-cavity line drawn on a hollow body: {counts}"
    else:
        assert counts[2] == 0, f"two transits without a cavity: {counts}"


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_transits_match_the_hand_counted_walls")
def test_an_absent_transit_cannot_be_read_as_a_slot_or_a_wall() -> None:
    """The absent codes are out of range: ``S`` for the slots, ``n + 1`` for the walls.

    On ``SPH`` (``n = 3``, ``S = 7``) a line at ``b = 0.7`` has one transit; its
    second column is absent. Indexing ``slot_length`` with the absent
    ``first_slot`` raises; the absent slot slice is empty; indexing a
    per-breakpoint wall table (length ``n + 1``) with either absent wall
    raises. A line missing the body has both columns absent.
    The positive leg: the present transit's codes index both tables.
    First red (run, T4): the slot fields filled with ``n + 1 = 4`` (slot 4 is
    the outbound region-0 slot: a silent read of a real length).

    Designed-green, measured: a tangency counted as a crossing in the kernel
    (``b <= r_k``, arm T6) reds NO row of this file, because the crossing
    pair it adds bounds a zero-length slot, which the transits ignore by
    definition; ``test_chord.py::test_a_tangency_is_not_a_crossing`` and three
    sibling rows catch it.
    """
    ch = SPH.chord(_line((0.7, -0.4, 0.0), (0.0, 1.0, 0.0)))
    tr = ch.transits
    wall_table = np.zeros(SPH.n_regions + 1)
    present, absent = 0, 1
    assert bool(tr.present[present]) and not bool(tr.present[absent])
    _ = ch.slot_length[int(tr.first_slot[present])], wall_table[int(tr.entry_wall[present])]
    with pytest.raises(IndexError):
        _ = ch.slot_length[int(tr.first_slot[absent])]
    assert ch.slot_length[int(tr.first_slot[absent]): int(tr.stop_slot[absent])].size == 0
    for wall in (tr.entry_wall[absent], tr.exit_wall[absent]):
        with pytest.raises(IndexError):
            _ = wall_table[int(wall)]
    miss = SPH.chord(_line((2.5, 0.0, 0.0), (0.0, 1.0, 0.0))).transits
    assert not miss.present.any()
    with pytest.raises(IndexError):
        _ = ch.slot_length[int(miss.first_slot[0])]
