"""Gates for the partition's refinement verbs: ``ConcentricPartition.refined_at`` and ``ConcentricPartition.coarsening``.

P3 step 3e of the characteristic-reference campaign (``.claude/plans/characteristic_reference_architecture.md``; the
elegance review's S2): the target-graded reading chords its lines through the basis's panels refined at the point's
level, and every consumer of such a chord maps the refined regions back onto the coarse ones (``PanelBasis.region_of_panel``,
``TraversalRule._panel_of_region``, the reading's split cross sections). The two verbs are the kernel's; their laws:

* ``refined_at(c)``: the same chart and pose, the breakpoints with ``c`` inserted in order, the same ends; the
  partition itself when ``c`` is already a breakpoint; refused outside ``[r_0, r_n]``;
* ``coarsening(coarse)``: the coarse region holding each fine region, read here from each fine region's midpoint by
  a hand loop; the identity on itself; refused unless the partition refines ``coarse`` (the same chart, the same
  ends, every coarse breakpoint among its own).

First reds: ``scratch/characteristic_architecture/p3/step3/ta/e/`` (``battery.log``).
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry import CoordSystem
from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition

_BODIES = {
    CoordSystem.SPHERICAL: (0.0, 0.5, 1.5, 2.0),
    CoordSystem.CYLINDRICAL: (0.4, 0.5, 1.5, 2.0),
    CoordSystem.CARTESIAN: (0.0, 0.4, 1.5, 2.3),
}
_COORDS = list(_BODIES)


def _partition(coord: CoordSystem, breakpoints=None) -> ConcentricPartition:
    return ConcentricPartition(Chart(coord), tuple(_BODIES[coord] if breakpoints is None else breakpoints))


def _by_hand(fine: tuple[float, ...], coarse: tuple[float, ...]) -> list[int]:
    """The coarse region holding each fine region's midpoint."""
    out = []
    for a, b in zip(fine[:-1], fine[1:]):
        mid = (a + b) / 2.0
        out.append(next(k for k in range(len(coarse) - 1) if coarse[k] <= mid <= coarse[k + 1]))
    return out


@pytest.mark.foundation
@pytest.mark.parametrize("coord", _COORDS, ids=[c.name for c in _COORDS])
def test_a_refinement_inserts_the_level_and_keeps_the_chart_the_pose_and_the_ends(coord) -> None:
    """[RF1] ``refined_at(c)`` at an interior level and at a level beside a breakpoint (one ulp above it) inserts ``c``
    in order, keeps the chart, the pose and both ends; at a breakpoint (an end or an interface) it is the partition
    itself (``is``). First reds: the level appended unsorted; a new partition minted at a breakpoint (a duplicate
    breakpoint, an empty region)."""
    partition = _partition(coord)
    ends = partition.breakpoints
    for c in (0.77, float(np.nextafter(ends[1], np.inf))):
        refined = partition.refined_at(c)
        assert refined.breakpoints == tuple(sorted((*ends, c)))
        assert refined.chart == partition.chart and refined.pose == partition.pose
        assert (refined.breakpoints[0], refined.breakpoints[-1]) == (ends[0], ends[-1])
    for c in ends:
        assert partition.refined_at(c) is partition


@pytest.mark.foundation
@pytest.mark.parametrize("coord", _COORDS, ids=[c.name for c in _COORDS])
def test_a_refinement_outside_the_ends_is_refused(coord) -> None:
    """[RF2] A level outside ``[r_0, r_n]`` (an ulp beyond either end) is refused, naming the interval. First red: the
    guard removed (the refinement gains a new end, and every wall re-keyed by its index moves)."""
    partition = _partition(coord)
    lo, hi = partition.breakpoints[0], partition.breakpoints[-1]
    for c in (float(np.nextafter(hi, np.inf)), float(np.nextafter(lo, -np.inf)) if lo > 0.0 else -1e-12):
        with pytest.raises(ValueError, match="a refinement keeps the ends"):
            partition.refined_at(c)


@pytest.mark.foundation
@pytest.mark.parametrize("coord", _COORDS, ids=[c.name for c in _COORDS])
@pytest.mark.rests_on("tests/gates/geometry/test_chord_refinement.py::test_a_refinement_inserts_the_level_and_keeps_the_chart_the_pose_and_the_ends")
def test_a_refinement_coarsens_onto_what_it_refines(coord) -> None:
    """[RF3] ``coarsening`` is the coarse region of each fine region (the hand loop on midpoints), ``array_equal``, for
    one and two inserted levels; it is the identity on the partition itself. First reds: ``side="left"`` (a fine
    region starting on a coarse breakpoint reads the region below, and the first wraps to -1: qa's ``panel_left``
    arm, 331 reds in the suite)."""
    coarse = _partition(coord)
    np.testing.assert_array_equal(coarse.coarsening(coarse), np.arange(coarse.n_regions))
    for levels in ((0.77,), (0.45, 1.9)):
        fine = coarse
        for c in levels:
            fine = fine.refined_at(c)
        np.testing.assert_array_equal(fine.coarsening(coarse), _by_hand(fine.breakpoints, coarse.breakpoints))


@pytest.mark.foundation
@pytest.mark.parametrize("coord", _COORDS, ids=[c.name for c in _COORDS])
def test_a_partition_coarsens_only_onto_one_it_refines(coord) -> None:
    """[RF4] Refused, each by the one fragment: another chart; other ends (the coarse body ending further out); a coarse
    breakpoint missing from the fine partition. Positive leg: RF3. First red: the guard removed (a map onto a
    partition the chord does not refine, read as panel indices)."""
    coarse = _partition(coord)
    fine = coarse.refined_at(0.77)
    other = next(c for c in _COORDS if c is not coord)
    wrong = [
        ConcentricPartition(Chart(other), coarse.breakpoints),
        _partition(coord, (*coarse.breakpoints[:-1], coarse.breakpoints[-1] + 0.5)),
        _partition(coord, (coarse.breakpoints[0], 1.0, coarse.breakpoints[-1])),
    ]
    for target in wrong:
        with pytest.raises(ValueError, match="coarsens only onto one it refines"):
            fine.coarsening(target)
