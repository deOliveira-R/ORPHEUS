"""Gates for :class:`orpheus.geometry.chord.ConcentricPartition` and its :class:`Chord`.

The geometric kernel seed (``.claude/plans/characteristic_reference_architecture.md``,
P0; verification spec ``scratch/characteristic_architecture/seed_verification_spec.md``
§5-§7a and §8, re-posed onto the built API: the directional locator is retired,
its rows live on ``Crossings``; the region of closest approach is two slots).

Every expected value is a closed form evaluated here, by hand or in mpmath at
50 dps from the fixture's own construction data (the impact parameter ``b`` is
the value the fixture was built at, never the kernel's). The chord of a
concentric body is symmetric about its closest approach, so a full-line
region/length assertion cannot see a reversed orientation (spec §2 regime 10):
orientation is gated by the crossings' surfaces and senses, by half-lines, and
by the cylinder's crossing parameters (whose closest approach is not the foot).
A total-length assertion cannot see a misplaced interior crossing (it
telescopes), so every slot carries its own reference.

Fixtures (spec §3): ``SPH`` (0, 0.3, 1.1, 2.0) solid, ``SPH_H`` (0.4, 1.1, 2.0)
hollow, the cylinders alike, ``SLB`` (-0.7, 0.3, 1.1, 2.0) asymmetric about 0.
"""
from __future__ import annotations

import math
import warnings

import mpmath as mp
import numpy as np
import pytest

from orpheus.geometry.chart import Chart
from orpheus.geometry.chart import RadialImage
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.line import Line
from orpheus.geometry.transformation import RigidMotion

_EPS = np.finfo(float).eps
_HERE = "tests/gates/geometry/test_chord.py::"
_LINE = "tests/gates/geometry/test_line.py::"
_CHART = "tests/gates/geometry/test_chart.py::"

_SPH_BP = (0.0, 0.3, 1.1, 2.0)
_SPH_H_BP = (0.4, 1.1, 2.0)
_SLB_BP = (-0.7, 0.3, 1.1, 2.0)


def _partition(coord: CoordSystem, breakpoints, pose: RigidMotion | None = None) -> ConcentricPartition:
    if pose is None:
        return ConcentricPartition(Chart(coord), breakpoints)
    return ConcentricPartition(Chart(coord), breakpoints, pose)


SPH = _partition(CoordSystem.SPHERICAL, _SPH_BP)
SPH_H = _partition(CoordSystem.SPHERICAL, _SPH_H_BP)
CYL = _partition(CoordSystem.CYLINDRICAL, _SPH_BP)
CYL_H = _partition(CoordSystem.CYLINDRICAL, _SPH_H_BP)
SLB = _partition(CoordSystem.CARTESIAN, _SLB_BP)
# The exterior codes are per partition and out of range for a per-region table
# (ruled 2026-10-05): inner n, outer n + 1; for n = 3 (SPH, CYL, SLB) 3 and 4, for n = 2 (the hollow ones) 2 and 3.
I3, O3 = 3, 4
I2, O2 = 2, 3


def _chord(partition: ConcentricPartition, p, omega):
    return partition.chord(Line.through(np.array(p, dtype=float), np.array(omega, dtype=float)))


def _impact_parameter(ch) -> np.ndarray:
    """The chord's impact parameter, which a radial image always carries (the slab's axial image has none)."""
    if not isinstance(ch.image, RadialImage):
        pytest.fail("a curvilinear chord's image is radial")
    return ch.image.impact_parameter


def _half_chords(ch) -> np.ndarray | None:
    """``h_k = sqrt((r_k - b)(r_k + b))`` where ``b < r_k``, else 0, derived HERE from ``b`` (the chord stores none); ``None`` on the slab."""
    if not isinstance(ch.image, RadialImage):
        return None
    r = np.asarray(ch.partition.breakpoints)
    b = ch.image.impact_parameter[..., None]
    return np.sqrt(np.where(b < r, (r - b) * (r + b), 0.0))


def _reference_slots(breakpoints, b: float, speed) -> list:
    """The slot lengths of a line at impact parameter ``b``, in mpmath, from the half-chords.

    Slots: regions inbound ``n-1 .. 0``, the cavity, regions outbound ``0 .. n-1``.
    A region crossed on both of its circles has ``h_{j+1} - h_j`` per side (an
    exact difference at 50 dps); the region of closest approach ``h_{j+1}`` per
    side; the cavity ``2 h_0``. Every length is divided by the projected speed.
    """
    mp.mp.dps = 50
    r = [mp.mpf(v) for v in breakpoints]
    bb = mp.mpf(b)
    h = [mp.sqrt(rk * rk - bb * bb) if bb < rk else mp.mpf(0) for rk in r]
    n = len(r) - 1
    side = []
    for j in range(n):
        if bb < r[j]:
            side.append(h[j + 1] - h[j])
        elif bb < r[j + 1]:
            side.append(h[j + 1])
        else:
            side.append(mp.mpf(0))
    cavity = 2 * h[0] if bb < r[0] else mp.mpf(0)
    s = mp.mpf(speed)
    return [x / s for x in side[::-1] + [cavity] + side]


def _assert_slots(got: np.ndarray, expected: list, nulp: float) -> None:
    for g, e in zip(got, expected, strict=True):
        if e == 0:
            assert g == 0.0, f"a slot the line does not traverse has length {g!r}"
        else:
            assert abs(mp.mpf(float(g)) - e) <= nulp * mp.mpf(np.spacing(float(e))), (
                f"slot length {g!r} vs {mp.nstr(e, 20)}"
            )


# ── R1/R2/C1: slot lengths against the closed form ─────────────────────────

_IN_PLANE = [  # (id, breakpoints, b)
    ("SPH_LC_b0", _SPH_BP, 0.0),
    ("SPH_LB_b0.7", _SPH_BP, 0.7),
    ("SPH_b0.2", _SPH_BP, 0.2),
    ("SPH_H_LB_b0.7", _SPH_H_BP, 0.7),
    ("SPH_H_LCAV_b0.2", _SPH_H_BP, 0.2),
    ("thin_shell_b0.5", (0.0, 1.1 * (1.0 - 1e-8), 1.1, 2.0), 0.5),
]
_SLOT_NULP = 4.0


@pytest.mark.l0
@pytest.mark.verifies("geometry-chord-segment-lengths")
@pytest.mark.rests_on(_CHART + "test_the_orbit_coordinate_of_each_chart",
                      _LINE + "test_the_line_does_not_depend_on_the_base_point_named")
@pytest.mark.parametrize(("bp", "b"), [(bp, b) for _, bp, b in _IN_PLANE], ids=[i for i, _, _ in _IN_PLANE])
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL], ids=["sphere", "cylinder"])
def test_every_slot_length_matches_the_closed_form(coord: CoordSystem, bp, b: float) -> None:
    """C1: each slot's length is its closed form, to 4 ulp of that slot.

    The kernel spells a shell's segment in the conditioned form
    ``(r_{k+1}^2 - r_k^2)/(h_{k+1} + h_k)``; the reference is the difference of
    half-chords at 50 dps, a different spelling of the same closed form. The
    thin shell (``r_1 = 1.1 (1 - 1e-8)``) is the conditioning row: by a
    difference of float half-chords it loses ``[M]`` 2e-8 relative (spec §2
    regime 4), far outside 4 ulp.
    First reds: a shell segment as a difference of float half-chords (thin
    shell); the region of closest approach given ``2 h`` in each slot; a slot
    filled for a region the line misses.
    """
    z = 1.5 if coord is CoordSystem.CYLINDRICAL else 0.0     # the sphere's b would include z
    ch = _chord(_partition(coord, bp), (b, -0.4, z), (0.0, 1.0, 0.0))
    _assert_slots(ch.slot_length, _reference_slots(bp, b, 1.0), _SLOT_NULP)


_CYL_DIRECTIONS = [
    (0.0, 0.6, 0.8), (0.0, 0.8, -0.6), (0.0, 1.0, 0.0),
    (0.0, 1e-6, float(np.sqrt(1.0 - 1e-12))),
]


@pytest.mark.l0
@pytest.mark.verifies("geometry-cylinder-axial-factor", "geometry-chord-segment-lengths")
@pytest.mark.rests_on(_HERE + "test_every_slot_length_matches_the_closed_form",
                      _CHART + "test_the_projected_speed_is_summed_from_components")
@pytest.mark.parametrize("z", [0.0, -3.7], ids=["z0", "z-3.7"])
@pytest.mark.parametrize("omega", _CYL_DIRECTIONS, ids=["s0.6", "s0.8_down", "s1", "near_axial"])
def test_a_cylinder_slot_is_the_in_plane_length_times_the_obliquity(omega, z: float) -> None:
    """C6: the 3-D length of a cylinder slot is its in-plane length over ``|P Omega|``.

    The directions are not closed under ``Omega_z -> -Omega_z`` and their
    ``|P Omega|`` are distinct (0.6, 0.8, 1, 1e-6), so a reversal or permutation
    of the obliquity across directions changes some row (the hoist's blind spot:
    a symmetric axial rule hid a reversed lift). Reference: the in-plane closed
    form divided by ``sqrt(Omega_x^2 + Omega_y^2)`` in mpmath.
    First reds: the obliquity multiplied instead of divided; omitted;
    ``|P Omega|^2`` as ``1 - Omega_z^2`` (the near-axial row, ``[M]`` 8.9e-5).
    """
    mp.mp.dps = 50
    om = np.array(omega)
    speed = mp.sqrt(mp.mpf(om[0]) ** 2 + mp.mpf(om[1]) ** 2)
    ch = _chord(CYL, (0.7, 0.0, z), om)
    _assert_slots(ch.slot_length, _reference_slots(_SPH_BP, 0.7, speed), 8.0)


# ── C3 / crossings: regions from the crossing order ────────────────────────

_SLOT_REGIONS = np.array([2, 1, 0, I3, 0, 1, 2])


@pytest.mark.l0
@pytest.mark.verifies("geometry-crossing-order", "geometry-line-crossing-law")
@pytest.mark.rests_on(_HERE + "test_every_slot_length_matches_the_closed_form")
def test_the_crossings_carry_the_region_they_enter() -> None:
    """C3: each present crossing's breakpoint, sense and region entered, by hand.

    ``L_C`` on ``SPH`` (b = 0): crossings at r_3, r_2, r_1 inward entering
    2, 1, 0, then r_1, r_2, r_3 outward entering 1, 2, outside; r_0 = 0 is never
    crossed. ``L_CAV`` on ``SPH_H`` (b = 0.2): inward through r_2, r_1, r_0
    entering 1, 0, the cavity; outward entering 0, 1, outside.
    The slot layout is fixed: regions inbound 2, 1, 0, the cavity, outbound 0, 1, 2.
    First reds: entering a circle inward moves one region OUT (every row);
    the cavity crossing labelled region 0.
    """
    ch = _chord(SPH, (0.0, -0.5, 0.0), (0.0, 1.0, 0.0))
    np.testing.assert_array_equal(ch.slot_region, _SLOT_REGIONS)
    c = ch.crossings
    on = np.broadcast_to(c.present, c.parameter.shape)
    np.testing.assert_array_equal(np.broadcast_to(c.breakpoint, on.shape)[on], [3, 2, 1, 1, 2, 3])
    np.testing.assert_array_equal(np.broadcast_to(c.sense, on.shape)[on], [-1, -1, -1, 1, 1, 1])
    np.testing.assert_array_equal(np.broadcast_to(c.region_entered, on.shape)[on], [2, 1, 0, 1, 2, O3])
    np.testing.assert_allclose(c.parameter[on], [-2.0, -1.1, -0.3, 0.3, 1.1, 2.0], rtol=0, atol=2 * _EPS)

    ch = _chord(SPH_H, (0.2, 0.0, 0.0), (0.0, 1.0, 0.0))
    np.testing.assert_array_equal(ch.slot_region, [1, 0, I2, 0, 1])
    c = ch.crossings
    on = np.broadcast_to(c.present, c.parameter.shape)
    np.testing.assert_array_equal(np.broadcast_to(c.breakpoint, on.shape)[on], [2, 1, 0, 0, 1, 2])
    np.testing.assert_array_equal(np.broadcast_to(c.region_entered, on.shape)[on], [1, 0, I2, 0, 1, O2])


_SURFACE_C = 16.0  # [M] 2026-10-05, seed 20261005: largest ratio 1.08 (sphere), 1.00, 0.77; restructured kernel 1.25, 0.96, 0.78; 10x rounded up


@pytest.mark.l0
@pytest.mark.verifies("geometry-crossing-order", "geometry-line-crossing-law")
@pytest.mark.rests_on(_HERE + "test_the_crossings_carry_the_region_they_enter")
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL, CoordSystem.CARTESIAN],
                         ids=["sphere", "cylinder", "slab"])
def test_every_crossing_lies_on_its_surface_and_moves_the_way_its_sense_says(coord: CoordSystem) -> None:
    """THEOREM over 400 seeded lines in a posed partition: crossings are where ``c = r_k``, in the direction ``sense``.

    Three routes that share no arithmetic with the chord's closed form: the
    point ``line.at(t)`` moved back by the pose and read by the chart's orbit
    coordinate must equal ``r_k``; ``c(t + d) - c(t - d)`` must have the sign
    ``sense``; the bare locator at ``t + d`` must return ``region_entered``.
    This is the orientation gate: a crossing parameter of the wrong sign still
    lies on the right surface on the sphere (a palindrome), and is caught by the
    sense and the region entered. Lines with ``h_min < 1e-3 R`` or
    ``|P Omega| < 1e-2`` are excluded (the problem is ill-conditioned there:
    spec §8); the count kept is asserted.
    First reds: the closest approach ``t*`` with the wrong sign (cylinder: the
    points leave their surfaces); the crossing parameters negated (sphere: the
    senses read backwards); ``region_entered`` off by one.
    """
    rng = np.random.default_rng(20261005)
    bp = _SLB_BP if coord is CoordSystem.CARTESIAN else _SPH_BP
    pose = RigidMotion(
        RigidMotion.rotation_about_axis(axis=(1.0, -2.0, 0.5), angle=0.7).linear, np.array([0.4, -1.2, 2.5]),
    )
    part = _partition(coord, bp, pose)
    om = rng.normal(size=(400, 3))
    om /= np.linalg.norm(om, axis=1, keepdims=True)
    line = Line.through(pose.on_points(rng.uniform(-2.5, 2.5, (400, 3))), pose.on_directions(om))
    ch = part.chord(line)
    keep = ch.projected_speed >= 1e-2
    hc = _half_chords(ch)
    if hc is not None:
        h = np.where(hc > 0, hc, np.inf)
        keep &= (h.min(axis=1) >= 2e-3) | ~np.isfinite(h.min(axis=1))
    kept = int(keep.sum())
    assert kept >= 300, f"only {kept} of 400 lines kept"
    c = ch.crossings
    shape = c.parameter.shape
    present = np.broadcast_to(c.present, shape) & keep[:, None]
    k = np.broadcast_to(c.breakpoint, shape)[present]
    sense = np.broadcast_to(c.sense, shape)[present]
    entered = np.broadcast_to(c.region_entered, shape)[present]
    t = c.parameter[present]
    rows = np.nonzero(present)[0]
    sub = Line(direction=line.direction[rows], moment=line.moment[rows])    # the CALLER's lines
    back = pose.inverse()
    rk = np.asarray(bp)[k]
    scale = np.linalg.norm(sub.foot, axis=1) + np.abs(t) + 4.0
    radial = part.chart.orbit_coordinate(back.on_points(sub.at(t)))
    assert np.all(np.abs(radial - rk) <= _SURFACE_C * _EPS * scale / ch.projected_speed[rows])
    d = 1e-6
    ahead = part.chart.orbit_coordinate(back.on_points(sub.at(t + d)))
    behind = part.chart.orbit_coordinate(back.on_points(sub.at(t - d)))
    np.testing.assert_array_equal(np.sign(ahead - behind), sense)
    np.testing.assert_array_equal(part.region_containing(ahead), entered)


@pytest.mark.l0
@pytest.mark.verifies("geometry-line-crossing-law")
@pytest.mark.rests_on(_HERE + "test_every_crossing_lies_on_its_surface_and_moves_the_way_its_sense_says")
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL, CoordSystem.CARTESIAN],
                         ids=["sphere", "cylinder", "slab"])
def test_the_orbit_coordinate_along_a_chord_is_the_charts(coord: CoordSystem) -> None:
    """``Chord.orbit_coordinate_at(t)`` equals the chart's coordinate of ``pose^-1(line.at(t))``, posed, at 400 x 5 seeded parameters.

    Two routes: the chord's closed form ``sqrt(b^2 + (|P Omega|(t - t*))^2)``
    (curvilinear) and the chart's map on the moved point. The band is
    ``16 eps (|foot| + |t| + 4)``. At the crossings it reads ``r_k``.
    First reds: ``t*`` not shifted onto the caller's line; ``|P Omega|``
    omitted from the travelled distance.
    """
    rng = np.random.default_rng(20261005)
    bp = _SLB_BP if coord is CoordSystem.CARTESIAN else _SPH_BP
    pose = _rot((1.0, -2.0, 0.5), 0.7, (0.4, -1.2, 2.5))
    part = _partition(coord, bp, pose)
    om = rng.normal(size=(400, 3))
    om /= np.linalg.norm(om, axis=1, keepdims=True)
    line = Line.through(pose.on_points(rng.uniform(-2.5, 2.5, (400, 3))), pose.on_directions(om))
    ch = part.chord(line)
    keep = ~np.asarray(ch.parallel)
    t = rng.uniform(-4, 4, (400, 5))
    got = ch.orbit_coordinate_at(t)
    pts = line.foot[:, None, :] + t[..., None] * line.direction[:, None, :]
    expected = part.chart.orbit_coordinate(pose.inverse().on_points(pts))
    band = 16 * _EPS * (np.linalg.norm(line.foot, axis=1)[:, None] + np.abs(t) + 4.0)
    assert np.all(np.abs(got - expected)[keep] <= band[keep])
    present = np.broadcast_to(ch.crossings.present, ch.crossings.parameter.shape)
    at_crossings = ch.orbit_coordinate_at(ch.crossings.parameter)
    rk = np.broadcast_to(np.asarray(bp)[ch.crossings.breakpoint], present.shape)
    scale = np.linalg.norm(line.foot, axis=1)[:, None] + np.abs(ch.crossings.parameter) + 4.0
    assert np.all((np.abs(at_crossings - rk) <= 16 * _EPS * scale)[present])


# ── R3 / R4 / R5: tangency, the centre, conditioning ───────────────────────


@pytest.mark.l0
@pytest.mark.verifies("geometry-line-crossing-law")
@pytest.mark.rests_on(_HERE + "test_every_slot_length_matches_the_closed_form")
def test_a_tangency_is_not_a_crossing() -> None:
    """R3: a line at ``b = r_k`` exactly does not cross ``c = r_k`` and stays in the outer region.

    ``L_T`` on ``SPH`` (b = r_2 = 1.1): only region 2 has length, in its two
    slots of ``sqrt(2.79)`` each, no crossing of r_2 is present. ``L_TR``
    (b = R): every slot 0, no crossing. ``L_TI`` on ``SPH_H`` (b = r_0 = 0.4): the
    cavity slot is 0, no crossing of r_0, region 0 has ``sqrt(1.21 - 0.16)`` per
    side. On the cylinder the same with an oblique direction (in-plane tangency).
    ``b`` equals ``r_k`` bit for bit (``|foot| = b`` exactly for these data), so
    the discriminant is exactly 0.
    First red: tangency tested with ``b <= r_k`` (a crossing pair at one point:
    two present crossings of r_2 on ``L_T``).
    """
    mp.mp.dps = 50
    for part, b, k, ref in (
        (SPH, 1.1, 2, _reference_slots(_SPH_BP, 1.1, 1)),
        (SPH, 2.0, 3, _reference_slots(_SPH_BP, 2.0, 1)),
        (SPH_H, 0.4, 0, _reference_slots(_SPH_H_BP, 0.4, 1)),
    ):
        ch = _chord(part, (b, 0.0, 0.0), (0.0, 1.0, 0.0))
        assert _impact_parameter(ch) == b
        _assert_slots(ch.slot_length, ref, _SLOT_NULP)
        c = ch.crossings
        on = np.broadcast_to(c.present, c.parameter.shape)
        assert not np.any(on & (np.broadcast_to(c.breakpoint, on.shape) == k))
    ch = _chord(SPH, (1.1, 0.0, 0.0), (0.0, 1.0, 0.0))
    np.testing.assert_array_equal(ch.slot_length > 0, ch.slot_region == 2)
    ch = _chord(CYL, (1.1, 0.0, 5.0), (0.0, 0.6, 0.8))
    np.testing.assert_array_equal(ch.slot_length > 0, ch.slot_region == 2)
    _assert_slots(ch.slot_length, _reference_slots(_SPH_BP, 1.1, mp.mpf("0.6")), 8.0)


@pytest.mark.l0
@pytest.mark.verifies("geometry-line-crossing-law")
@pytest.mark.rests_on(_HERE + "test_every_slot_length_matches_the_closed_form")
@pytest.mark.parametrize("part", [SPH, CYL], ids=["sphere", "cylinder"])
def test_the_centre_of_a_solid_body_is_never_crossed(part: ConcentricPartition) -> None:
    """R4: ``r_0 = 0`` is the singular stratum, not a surface.

    A line through the centre (b = 0) crosses r_3, r_2, r_1 twice each and r_0
    never; region 0 has 0.3 per side, split at the centre; the cavity slot is 0.
    First red: ``r_0`` included in the crossing set (a crossing pair at ``t*``,
    region entered ``INNER_EXTERIOR`` on a solid body).
    """
    ch = _chord(part, (0.0, -0.5, 0.0), (0.0, 1.0, 0.0))
    np.testing.assert_allclose(ch.slot_length, [0.9, 0.8, 0.3, 0.0, 0.3, 0.8, 0.9], rtol=4 * _EPS, atol=0)
    c = ch.crossings
    on = np.broadcast_to(c.present, c.parameter.shape)
    assert int(on.sum()) == 6
    assert not np.any(on & (np.broadcast_to(c.breakpoint, on.shape) == 0))
    assert not np.any(np.broadcast_to(c.region_entered, on.shape)[on] == part.inner_exterior)


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_every_slot_length_matches_the_closed_form")
def test_the_half_chord_is_well_conditioned_at_near_tangency() -> None:
    """R5: at ``b = R (1 - 1e-12)`` the half-chord is within 2 ulp of its exact value.

    ``sqrt(R*R - b*b)`` errs there by up to ``[M]`` 2.8e14 ulp (spec §2 regime 3);
    ``sqrt((R - b)(R + b))`` by at most 1.16 ulp. The promised invariant
    (lengths within 4 ulp, C1) and the realised quality are separate rows
    (anti-pattern #16).
    First red: the naive half-chord.
    """
    mp.mp.dps = 50
    R = 2.0
    b = R * (1.0 - 1e-12)
    ch = _chord(_partition(CoordSystem.SPHERICAL, (0.0, R)), (b, 0.0, 0.0), (0.0, 1.0, 0.0))
    exact = mp.sqrt(mp.mpf(R) ** 2 - mp.mpf(b) ** 2)
    for got in (ch.slot_length[0], ch.slot_length[2]):
        assert abs(mp.mpf(float(got)) - exact) <= 2 * mp.mpf(np.spacing(float(exact)))


# ── C5: half-lines ─────────────────────────────────────────────────────────


@pytest.mark.l0
@pytest.mark.verifies("geometry-chord-segment-lengths")
@pytest.mark.rests_on(_HERE + "test_the_crossings_carry_the_region_they_enter")
def test_a_half_line_keeps_the_slots_beyond_its_start() -> None:
    """C5: ``lengths_beyond(start)`` from an interior point, both orientations.

    From ``(0, 0.7, 0)`` in region 1 of ``SPH``: towards the centre, region 1
    for 0.4, region 0 for 0.3 + 0.3, region 1 for 0.8, region 2 for 0.9; away
    from it, region 1 for 0.4 and region 2 for 0.9. The two orientations differ,
    so the chord's palindrome cannot hide a reversal here.
    On ``SLB`` from x = 0 along (0.6, 0.8, 0): region 0 for 0.3/0.6, then
    0.8/0.6 and 0.9/0.6; reversed: region 0 for 0.7/0.6 only.
    First reds: ``lengths_beyond`` keeping the slots BEFORE the start; the cut
    slot given its full length.
    """
    p = np.array([0.0, 0.7, 0.0])
    for omega, expected in (
        ((0.0, -1.0, 0.0), [0.0, 0.4, 0.3, 0.0, 0.3, 0.8, 0.9]),
        ((0.0, 1.0, 0.0), [0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.9]),
    ):
        line = Line.through(p, np.array(omega))
        got = SPH.chord(line).lengths_beyond(line.parameter_of(p))
        np.testing.assert_allclose(got, expected, rtol=0, atol=8 * _EPS)
    p = np.zeros(3)
    for omega, regions, expected in (
        ((0.6, 0.8, 0.0), [0, 1, 2], [0.3 / 0.6, 0.8 / 0.6, 0.9 / 0.6]),
        ((-0.6, 0.8, 0.0), [2, 1, 0], [0.0, 0.0, 0.7 / 0.6]),
    ):
        line = Line.through(p, np.array(omega))
        ch = SLB.chord(line)
        np.testing.assert_array_equal(ch.slot_region, regions)
        np.testing.assert_allclose(ch.lengths_beyond(line.parameter_of(p)), expected, rtol=0, atol=8 * _EPS)


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_a_half_line_keeps_the_slots_beyond_its_start")
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL, CoordSystem.CARTESIAN],
                         ids=["sphere", "cylinder", "slab"])
def test_a_posed_half_line_reads_its_start_on_the_callers_line(coord: CoordSystem) -> None:
    """C5 under a pose: ``chord.lengths_beyond(line.parameter_of(p))`` on the CALLER's line is the canonical answer.

    The pose is a rotation by sqrt(2) rad about a skew axis with a translation,
    so the parameter shift onto the caller's line is rounded: the band is
    ``16 eps (|p| + 4)``. With the identity pose the
    caller's line and the chord's canonical line share a foot, which is why only
    a posed partition can see the parameter frame.
    First red (measured 2026-10-05, before repair (a)): the chord's parameters were measured on the
    canonical line, so the caller's start was off by the foot's shift (sphere:
    region 2 inbound read 0.1 and region 1 read 0.8 where 0 and 0.4 are right).
    """
    bp = _SLB_BP if coord is CoordSystem.CARTESIAN else _SPH_BP
    M = _rot((1.0, -2.0, 0.5), np.sqrt(2.0), (1.0, 0.5, 0.25))
    canonical, posed = _partition(coord, bp), _partition(coord, bp, M)
    p = np.array([0.0, 0.7, 0.0]) if coord is not CoordSystem.CARTESIAN else np.zeros(3)
    # No direction parallel to the slab's planes: a rounded rotation makes an
    # exactly parallel line 3e-17 off parallel, with slots of length 1e16 (the
    # obliquity's own conditioning, spec §8), not the canonical infinite slot.
    for omega in ((0.28, -0.96, 0.0), (0.6, 0.0, -0.8), (-0.6, 0.8, 0.0)):
        om = np.array(omega)
        here, there = Line.through(p, om), Line.through(M.on_points(p), M.on_directions(om))
        expected = canonical.chord(here).lengths_beyond(here.parameter_of(p))
        got = posed.chord(there).lengths_beyond(there.parameter_of(M.on_points(p)))
        finite = np.isfinite(expected)
        np.testing.assert_array_equal(np.isfinite(got), finite)
        np.testing.assert_allclose(got[finite], expected[finite], rtol=0, atol=16 * _EPS * (np.linalg.norm(p) + 4.0))


# ── the multi-region segment table, by hand ──────────────────────────────

#: The body of the table: a solid sphere or cylinder with breakpoints (0, 0.5, 1.5, 2.0), regions 0, 1, 2.
_MR_BP = (0.0, 0.5, 1.5, 2.0)


def _half_chord(r: float, b: float) -> float:
    """The half-chord of the circle of radius ``r`` at impact parameter ``b``, ``sqrt((r - b)(r + b))``."""
    return math.sqrt((r - b) * (r + b))


def _outbound_by_hand(b: float) -> tuple[float, float, float]:
    """The OUTBOUND orbit-space length of each region (0, 1, 2) of a line at impact parameter ``b`` through ``_MR_BP``.

    Written from the half-chords, case by case: the region of closest approach holds its half-chord, every region
    outside it the difference of its two half-chords, and a region inside it 0. At ``b = 0.5`` and ``b = 1.5`` the
    line is tangent to an interface: a tangency is not a crossing, so the region inside is 0.
    """
    if b < 0.5:
        return (_half_chord(0.5, b), _half_chord(1.5, b) - _half_chord(0.5, b), _half_chord(2.0, b) - _half_chord(1.5, b))
    if b < 1.5:
        return (0.0, _half_chord(1.5, b), _half_chord(2.0, b) - _half_chord(1.5, b))
    return (0.0, 0.0, _half_chord(2.0, b))


#: The impact parameters: the centre, one per region, and the two interfaces exactly (tangent lines).
_MR_B = [("b0", 0.0), ("b0.3", 0.3), ("b0.5_tangent_r1", 0.5), ("b1.1", 1.1), ("b1.5_tangent_r2", 1.5), ("b1.9", 1.9)]
#: (id, coordinate system, direction, in-plane speed |P Omega|): the cylinder's 3-D lengths are the table's over it.
_MR_LINES = [
    ("sphere", CoordSystem.SPHERICAL, (0.0, 1.0, 0.0), 1.0),
    ("cylinder_wz0", CoordSystem.CYLINDRICAL, (0.0, 1.0, 0.0), 1.0),
    ("cylinder_wz0.8", CoordSystem.CYLINDRICAL, (0.0, 0.6, 0.8), 0.6),
]
#: The band of each slot: ``_MR_NULP`` ulp of the slot, plus the slot's sensitivity to the impact parameter times
#: ``_MR_B_ULP`` ulp of ``b``. The line normalises its direction, so on the oblique cylinder the chord's impact
#: parameter is 2 ulp from the ``b`` the table was written at (``[M]`` 2026-10-10: 1.9000000000000004), and near a
#: tangency a slot amplifies that by ``b^2 / h^2`` (9.3 at ``b = 1.9``: 10 ulp in the slot). ``[M]`` 2026-10-10
#: (``scratch/characteristic_architecture/p1_step_e/ta_e1b/consumers/segment_table_ulps.py``): with the
#: sensitivity term the largest full-chord gap is 0.41 of the band over the 18 lines (10 ulp raw).
_MR_NULP, _MR_B_ULP = 4.0, 4.0


def _slot_band(b: float, speed: float, slot: int) -> float:
    """The band of outbound slot ``slot`` at ``b``: ulp of the slot plus its central-difference slope in ``b`` times ulp of ``b``."""
    w = _outbound_by_hand(b)[slot] / speed
    if b == 0.0:
        return _MR_NULP * float(np.spacing(w))
    step = 1e-7 * b
    slope = (_outbound_by_hand(b + step)[slot] - _outbound_by_hand(b - step)[slot]) / (2.0 * step * speed)
    return _MR_NULP * float(np.spacing(w)) + abs(slope) * _MR_B_ULP * float(np.spacing(b))


@pytest.mark.l0
@pytest.mark.verifies("geometry-chord-segment-lengths", "peierls-greens-cylinder-mr-trajectory-segments",
                      "peierls-greens-cylinder-trajectory")
@pytest.mark.rests_on(_HERE + "test_every_slot_length_matches_the_closed_form",
                      _HERE + "test_a_cylinder_slot_is_the_in_plane_length_times_the_obliquity",
                      _HERE + "test_the_crossings_carry_the_region_they_enter",
                      _HERE + "test_a_half_line_keeps_the_slots_beyond_its_start")
@pytest.mark.parametrize("b", [b for _, b in _MR_B], ids=[i for i, _ in _MR_B])
@pytest.mark.parametrize(("coord", "omega", "speed"), [r[1:] for r in _MR_LINES], ids=[r[0] for r in _MR_LINES])
def test_the_multi_region_segments_are_the_hand_written_table(coord: CoordSystem, omega, speed: float, b: float) -> None:
    """The kernel's segments through a three-region body against a hand-written table of closed-form chord lengths.

    The successor of ``test_kernel_corroboration.py::test_the_backward_segments_agree_with_variant_alpha``
    (P1 step (e), ruling 4 of 2026-10-10): the old trajectory-resolvent oracle was the only other multi-region
    segment spelling, so the corroboration is replaced by a REFERENCE written in this test
    (:func:`_outbound_by_hand`). Three legs per line:

    1. the full chord: the slot regions ``(2, 1, 0, cavity, 0, 1, 2)`` and each slot's length, inbound the mirror of
       outbound, divided by the in-plane speed on the cylinder;
    2. the backward first leg from the foot of the perpendicular, ``lengths_beyond(0)``: the inbound slots 0 and the
       outbound slots whole (the half-line's orientation, which a full chord's palindrome hides);
    3. on the lines through region 1 (b < 1.5), the half-line from the point 0.8 back in the plane (radius
       ``sqrt(b^2 + 0.64)``, in region 1): region 1's inbound slot cut at the point, every later slot whole.

    The reference is the difference of float half-chords, a spelling the kernel does not use (it forms
    ``(r_{k+1}^2 - r_k^2) / (h_{k+1} + h_k)``); the band is :func:`_slot_band`. An untraversed slot is 0 exactly.
    First reds (``consumers/segment_battery.py``): the region of closest approach given ``2 h`` per side; the
    cylinder's obliquity dropped; the outbound slots written in the inbound order.

    ``peierls-greens-cylinder-mr-trajectory-segments`` (the in-plane conic ``r(s)^2 = r_0^2 - 2 r_0 cos(phi) s + s^2``)
    is verified at its roots: each slot end is where the conic meets a breakpoint, on legs 2 and 3 from two starting
    radii. ``peierls-greens-cylinder-trajectory`` (``L_2D = r cos(phi) + sqrt(R^2 - r^2 sin^2(phi))``, ``L_3D = L_2D /
    |P Omega|``) is leg 3's slots summed: from the start ``r cos(phi) = 0.8`` and ``r sin(phi) = b``, they telescope to
    ``(0.8 + sqrt(R^2 - b^2)) / |P Omega|``, each slot asserted on its own, the obliquity on the Omega_z = 0.8 rows.
    """
    part = _partition(coord, _MR_BP)
    lengths = np.array(_outbound_by_hand(b)) / speed
    bands = np.array([_slot_band(b, speed, slot) for slot in range(3)])
    p = np.array([b, 0.0, 0.0])
    line = Line.through(p, np.array(omega))
    ch = part.chord(line)
    np.testing.assert_array_equal(ch.slot_region, [2, 1, 0, I3, 0, 1, 2])

    def assert_lengths(got: np.ndarray, want: np.ndarray, band: np.ndarray, what: str) -> None:
        for slot, (g, w, tol) in enumerate(zip(got, want, band, strict=True)):
            if w == 0.0:
                assert g == 0.0, f"{what}, slot {slot}: the line does not traverse it, got {g!r}"
            else:
                assert abs(g - w) <= tol, f"{what}, slot {slot}: kernel {g!r}, hand {w!r}, gap {abs(g - w) / tol:.2f} of the band"

    zero = np.zeros(1)
    assert_lengths(ch.slot_length, np.concatenate([lengths[::-1], zero, lengths]),
                   np.concatenate([bands[::-1], zero, bands]), "full chord")
    assert_lengths(ch.lengths_beyond(line.parameter_of(p)), np.concatenate([np.zeros(4), lengths]),
                   np.concatenate([np.zeros(4), bands]), "from the foot")
    if b < 1.5:
        start = p - 0.8 / speed * np.array(omega)
        h1 = lengths[0]                                          # 0 where the line misses region 0 (b >= 0.5)
        want = np.array([0.0, 0.8 / speed - h1, h1, 0.0, *lengths])
        cut_band = _MR_NULP * float(np.spacing(want[1])) + bands[0]
        assert_lengths(ch.lengths_beyond(line.parameter_of(start)), want,
                       np.array([0.0, cut_band, bands[0], 0.0, *bands]), "from 0.8 back")


# ── C7: the slab ───────────────────────────────────────────────────────────


@pytest.mark.l0
@pytest.mark.verifies("geometry-chord-segment-lengths", "geometry-cylinder-axial-factor", "peierls-greens-slab-trajectory")
@pytest.mark.rests_on(_CHART + "test_the_projected_speed_is_summed_from_components")
def test_the_slab_chord_in_both_orientations() -> None:
    """C7: the slab's crossings ``t_k = (r_k - x_foot)/Omega_x`` and lengths ``(r_{k+1} - r_k)/|Omega_x|``.

    ``SLB`` is asymmetric about 0 (r_0 = -0.7), so the two orientations differ.
    Along (0.6, 0.8, 0) through the origin: regions 0, 1, 2, lengths 1/0.6,
    0.8/0.6, 0.9/0.6, crossings at r_k/0.6 entering 0, 1, 2, outside. Reversed:
    regions 2, 1, 0, crossings at -r_k/0.6 in the order r_3 .. r_0, entering
    2, 1, 0, the inner exterior.
    First reds: ``|Omega_x|`` replaced by ``Omega_x`` (negative lengths on the
    reversed line); breakpoints assumed non-negative.
    ``peierls-greens-slab-trajectory`` (``L_first = (x - x_0)/mu`` for mu > 0, ``(x_n - x)/|mu|`` for mu < 0) is the
    crossing parameter of the wall behind the point x = 0 on each orientation (``-r_0/0.6`` and ``r_3/0.6``).
    """
    r = np.array(_SLB_BP)
    for omega, order, entered in (
        ((0.6, 0.8, 0.0), [0, 1, 2], [0, 1, 2, O3]),
        ((-0.6, 0.8, 0.0), [2, 1, 0], [2, 1, 0, I3]),
    ):
        ch = _chord(SLB, (0.0, 0.0, 0.0), omega)
        np.testing.assert_array_equal(ch.slot_region, order)
        np.testing.assert_allclose(ch.slot_length, np.diff(r)[order] / 0.6, rtol=2 * _EPS, atol=0)
        c = ch.crossings
        ks = [0, 1, 2, 3] if omega[0] > 0 else [3, 2, 1, 0]
        np.testing.assert_array_equal(c.breakpoint, ks)
        np.testing.assert_array_equal(c.region_entered, entered)
        np.testing.assert_allclose(c.parameter, r[ks] / omega[0], rtol=2 * _EPS, atol=0)
        assert np.all(np.diff(c.parameter) > 0)


# ── C8: parallel lines ─────────────────────────────────────────────────────


@pytest.mark.l0
@pytest.mark.verifies("geometry-crossing-order", "geometry-line-crossing-law")
def test_a_line_parallel_to_the_orbit_space_is_inside_one_region_or_on_an_interface() -> None:
    """C8: ``|P Omega| = 0``: no crossing; one infinite slot, or ``interface = k`` with every slot 0.

    Cylinder, ``Omega = e_z``: at rho = 0.7 the outbound slot of region 1 is
    infinite; at rho = 1.1 (exactly r_2) ``interface == 2`` and every slot is 0;
    at rho = 2.5 every slot is 0 and there is no interface; on the axis of the
    solid body region 0 is infinite and the axis is NOT an interface (r_0 = 0 is
    the singular stratum); on the axis of the hollow body the cavity slot is
    infinite. Slab, ``Omega = e_y``: x = 0.5 gives region 1 infinite; x = 0.3,
    -0.7 and 2.0 give interfaces 1, 0 and 3.
    First reds: a line on a breakpoint assigned to its inner side (inner-owns
    folded in, against the ruling of 2026-10-05); the solid axis typed as
    interface 0.
    """
    def infinite_regions(ch):
        return sorted(int(v) for v in np.asarray(ch.slot_region)[np.isinf(ch.slot_length)])

    none = None                                   # stands for the partition's no_interface code (out of range)
    for part, rho, infinite, iface in (
        (CYL, 0.7, [1], none), (CYL, 1.1, [], 2), (CYL, 2.5, [], none),
        (CYL, 0.0, [0], none), (CYL_H, 0.0, [I2], none), (CYL_H, 0.4, [], 0),
    ):
        ch = _chord(part, (rho, 0.0, 3.0), (0.0, 0.0, 1.0))
        assert bool(ch.parallel) and int(ch.interface) == (part.no_interface if iface is None else iface)
        assert infinite_regions(ch) == infinite
        assert np.all(np.isinf(ch.slot_length) | (ch.slot_length == 0.0))
        assert not np.any(np.broadcast_to(ch.crossings.present, ch.crossings.parameter.shape))
    for x, infinite, iface in ((0.5, [1], None), (0.3, [], 1), (-0.7, [], 0), (2.0, [], 3), (2.5, [], None)):
        ch = _chord(SLB, (x, 0.0, 0.0), (0.0, 1.0, 0.0))
        assert bool(ch.parallel) and int(ch.interface) == (SLB.no_interface if iface is None else iface)
        assert infinite_regions(ch) == infinite
    assert SLB.no_interface == len(_SLB_BP) and CYL_H.no_interface == len(_SPH_H_BP)    # out of range for the breakpoints


# ── C4: the route: no slot is located ──────────────────────────────────────


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_crossings_carry_the_region_they_enter",
                      _HERE + "test_a_line_parallel_to_the_orbit_space_is_inside_one_region_or_on_an_interface")
def test_no_crossing_line_reads_the_bare_locator(monkeypatch: pytest.MonkeyPatch) -> None:
    """C4 (route gate, anti-pattern #26): a crossing line's slots and crossings do not depend on point location.

    The bare locator is replaced by a decoy that answers the outer exterior
    everywhere. Every non-parallel chord is unchanged, array for array (its
    regions come from the crossing order: ruled 2026-10-05). Activation leg: a
    parallel line, whose home region IS located, moves (its infinite slot
    disappears), so the decoy reaches the kernel.
    First red: slots labelled by locating their midpoints.
    """
    lines = [
        (SPH, (0.0, -0.5, 0.0), (0.0, 1.0, 0.0)), (SPH_H, (0.2, 0.0, 0.0), (0.0, 1.0, 0.0)),
        (CYL, (0.7, 0.0, 1.0), (0.0, 0.6, 0.8)), (SLB, (0.0, 0.0, 0.0), (-0.6, 0.8, 0.0)),
    ]
    # Materialised before the patch: slot_region and slot_length are properties, read lazily.
    honest = [
        (c.slot_region.copy(), c.slot_length.copy(), c.crossings.region_entered.copy())
        for c in (_chord(*row) for row in lines)
    ]
    parallel_honest_infinite = bool(np.isinf(_chord(CYL, (0.7, 0.0, 3.0), (0.0, 0.0, 1.0)).slot_length).any())

    def decoy(self, orbit_coordinate):
        return np.full(np.shape(orbit_coordinate), self.outer_exterior)

    monkeypatch.setattr(ConcentricPartition, "region_containing", decoy)
    for row, (region, length, entered) in zip(lines, honest):
        after = _chord(*row)
        np.testing.assert_array_equal(after.slot_region, region)
        np.testing.assert_array_equal(after.slot_length, length)
        np.testing.assert_array_equal(after.crossings.region_entered, entered)
    parallel_decoyed = _chord(CYL, (0.7, 0.0, 3.0), (0.0, 0.0, 1.0))
    assert parallel_honest_infinite and not np.isinf(parallel_decoyed.slot_length).any()


# ── C9 / C10: invariance ───────────────────────────────────────────────────

# Measured on the kernel (seed 20261005, 2000 lines per chart, 2026-10-05): the
# largest ratio of a slot's change to eps * scale was 4.99 (sphere, group) and
# 4.10 (cylinder, pose) on the first kernel; 4.30 (cylinder flip) and 3.15
# (cylinder, pose) on the restructured one. The constant is 10x the largest, rounded up.
_INVARIANCE_C = 50.0


def _scale(ch, p: np.ndarray) -> np.ndarray:
    """``R (|p| + R) / (h_min |P Omega|)`` per line (curvilinear); ``(|p| + R) / |Omega_x|^2`` (slab). Spec §8."""
    R = 2.0
    pn = np.linalg.norm(p, axis=1)
    hc = _half_chords(ch)
    if hc is None:
        return (pn + R) / ch.projected_speed ** 2
    h = np.where(hc > 0, hc, np.inf).min(axis=1)
    return R * (pn + R) / h / ch.projected_speed


def _invariance_lines(part: ConcentricPartition, n: int = 2000):
    rng = np.random.default_rng(20261005)
    om = rng.normal(size=(n, 3))
    om /= np.linalg.norm(om, axis=1, keepdims=True)
    p = rng.uniform(-3, 3, (n, 3))
    line = Line.through(p, om)
    ch = part.chord(line)
    keep = ch.projected_speed >= 1e-3
    hc = _half_chords(ch)
    if hc is not None:
        r = np.asarray(part.breakpoints)[1:]
        b = _impact_parameter(ch)
        h = np.where(hc > 0, hc, np.inf).min(axis=1)
        keep &= (np.abs(b[:, None] - r).min(axis=1) >= 2e-6) & (h >= 2e-3) & np.isfinite(h)
    idx = np.nonzero(keep)[0]
    return Line(direction=line.direction[idx], moment=line.moment[idx]), p[idx]


def _rot(axis, angle, shift=(0.0, 0.0, 0.0)):
    return RigidMotion(RigidMotion.rotation_about_axis(axis=axis, angle=angle).linear, np.array(shift, dtype=float))


_GROUP = [
    ("sphere_rotation", CoordSystem.SPHERICAL, _rot((1.0, 2.0, -0.5), np.sqrt(2.0)), True),
    ("sphere_mirror", CoordSystem.SPHERICAL, RigidMotion(RigidMotion.reflection(normal=(0.3, -1.0, 0.2)).linear), True),
    ("sphere_translated", CoordSystem.SPHERICAL, RigidMotion.translation_by([0.3, 0.0, 0.0]), False),
    ("cylinder_rotation_shift", CoordSystem.CYLINDRICAL, _rot((0.0, 0.0, 1.0), np.sqrt(2.0), (0.0, 0.0, 3.7)), True),
    ("cylinder_flip", CoordSystem.CYLINDRICAL, _rot((1.0, 1.0, 0.0), np.pi), True),
    ("cylinder_tilted", CoordSystem.CYLINDRICAL, _rot((1.0, 0.0, 0.0), 0.4), False),
    ("slab_rotation_shift", CoordSystem.CARTESIAN, _rot((1.0, 0.0, 0.0), np.sqrt(2.0), (0.0, 1.3, -0.4)), True),
    ("slab_tilted", CoordSystem.CARTESIAN, _rot((0.0, 0.0, 1.0), 0.4), False),
]


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_every_slot_length_matches_the_closed_form",
                      _CHART + "test_membership_in_the_symmetry_group_is_invariance_of_the_orbit_coordinate")
@pytest.mark.parametrize(("coord", "motion", "member"), [g[1:] for g in _GROUP], ids=[g[0] for g in _GROUP])
def test_the_chord_is_invariant_under_the_symmetry_group(coord: CoordSystem, motion: RigidMotion, member: bool) -> None:
    """C9 (THEOREM): a line moved by an element of ``G_c`` has the same chord.

    Labels exact; lengths within ``50 eps R (|p| + R)/(h_min |P Omega|)`` (the
    problem's conditioning at tangency, spec §8; constant re-measured on the
    kernel). Lines near a tangency (``|b - r_k| < 2e-6``) or with ``h_min < 2e-3``
    are excluded, so no slot can appear or vanish by rounding. Negative legs
    (a motion outside ``G_c``) must move some slot by more than 1e-6: the gate
    is loaded (anti-pattern #19).
    First reds: the sphere's coordinate read about a displaced centre; the
    cylinder's including z.
    """
    part = _partition(coord, _SLB_BP if coord is CoordSystem.CARTESIAN else _SPH_BP)
    assert part.chart.contains(motion) is member
    line, p = _invariance_lines(part)
    assert line.shape[0] >= 500, line.shape
    before = part.chord(line)
    after = part.chord(line.moved_by(motion))
    if member:
        np.testing.assert_array_equal(after.slot_region, before.slot_region)
        drift = np.abs(after.slot_length - before.slot_length).max(axis=1)
        assert np.all(drift <= _INVARIANCE_C * _EPS * _scale(before, p))
    else:
        assert np.max(np.abs(after.slot_length - before.slot_length)) > 1e-6


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_chord_is_invariant_under_the_symmetry_group")
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL, CoordSystem.CARTESIAN],
                         ids=["sphere", "cylinder", "slab"])
def test_the_chord_is_unchanged_when_the_line_and_the_partition_move_together(coord: CoordSystem) -> None:
    """C10 (THEOREM): ``chord(pose; M line) == chord(identity; line)`` for a general rigid motion ``M``.

    The partition is posed by ``M`` (a rotation by sqrt(2) rad about a skew
    axis, with a translation) and the line is moved by the same ``M``. Labels
    exact; lengths within the C9 band. The posed partition's bare locator
    answers the canonical one's at the moved points. Negative leg: the line
    moved and the partition not (or the converse) changes the chord.
    First reds: the pose applied instead of its inverse; the pose ignored.
    """
    bp = _SLB_BP if coord is CoordSystem.CARTESIAN else _SPH_BP
    M = _rot((1.0, -2.0, 0.5), np.sqrt(2.0), (0.4, -1.2, 2.5))
    canonical, posed = _partition(coord, bp), _partition(coord, bp, M)
    line, p = _invariance_lines(canonical)
    before = canonical.chord(line)
    after = posed.chord(line.moved_by(M))
    np.testing.assert_array_equal(after.slot_region, before.slot_region)
    drift = np.abs(after.slot_length - before.slot_length).max(axis=1)
    assert np.all(drift <= _INVARIANCE_C * _EPS * _scale(before, p))
    assert np.max(np.abs(posed.chord(line).slot_length - before.slot_length)) > 1e-6
    x = np.random.default_rng(20261005).uniform(-2.5, 2.5, (200, 3))
    c = canonical.chart.orbit_coordinate(x)
    away = np.abs(c[:, None] - np.asarray(bp)).min(axis=1) > 1e-9        # no point within rounding of a surface
    assert away.sum() >= 190
    np.testing.assert_array_equal(posed.region_at(M.on_points(x[away])), canonical.region_at(x[away]))


# ── L1: the bare locator ───────────────────────────────────────────────────

_UP, _DOWN = np.inf, -np.inf


def _bracket(v: float):
    return [np.nextafter(v, _DOWN), v, np.nextafter(v, _UP)]


@pytest.mark.foundation
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL], ids=["sphere", "cylinder"])
def test_the_bare_locator_is_inner_owns_on_the_closed_domain(coord: CoordSystem) -> None:
    """L1 (ruled 2026-10-05): region j is ``(r_j, r_{j+1}]``, region 0 is ``[r_0, r_1]``; below and above are typed exteriors.

    Exact on-surface points on every axis (``sqrt(r*r) == r`` in binary64) and
    the neighbouring float on each side of every breakpoint, so the comparison
    operator itself is pinned, not only its side. Solid: the centre is region 0.
    Hollow: ``r_0`` is region 0, the float below it and the centre are the
    inner exterior. The cylinder's points carry z in {-3.7, 1e3}.
    First reds: outer-owns (``searchsorted`` side flipped); the domain open at
    R or at r_0; the exteriors merged.
    """
    solid, hollow = _partition(coord, _SPH_BP), _partition(coord, _SPH_H_BP)
    expected_solid = {0.0: 0, 0.3: [0, 0, 1], 1.1: [1, 1, 2], 2.0: [2, 2, O3]}
    expected_hollow = {0.0: I2, 0.4: [I2, 0, 0], 1.1: [0, 0, 1], 2.0: [1, 1, O2]}
    for part, table in ((solid, expected_solid), (hollow, expected_hollow)):
        for v, exp in table.items():
            values = [v] if not isinstance(exp, list) else _bracket(v)
            exps = [exp] if not isinstance(exp, list) else exp
            for value, e in zip(values, exps):
                np.testing.assert_array_equal(part.region_containing(np.array(value)), e)
                axes = ((value, 0.0, -3.7), (0.0, value, 1e3)) if coord is CoordSystem.CYLINDRICAL else (
                    (value, 0.0, 0.0), (0.0, value, 0.0), (0.0, 0.0, -value))
                for x in axes:
                    np.testing.assert_array_equal(part.region_at(np.array(x)), e)


@pytest.mark.foundation
def test_the_slab_bare_locator_is_inner_owns_on_the_closed_domain() -> None:
    """L1, slab: ``x = -0.7`` is region 0 and the float below it the inner exterior; ``x = 2.0`` is region 2."""
    for v, exp in ((-0.7, [I3, 0, 0]), (0.3, [0, 0, 1]), (1.1, [1, 1, 2]), (2.0, [2, 2, O3])):
        for value, e in zip(_bracket(v), exp):
            np.testing.assert_array_equal(SLB.region_containing(np.array(value)), e)
            np.testing.assert_array_equal(SLB.region_at(np.array([value, 9.0, -9.0])), e)
    np.testing.assert_array_equal(SLB.region_containing(np.array([-5.0, -0.3])), [I3, 0])    # the slab admits negatives


@pytest.mark.foundation
@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf], ids=["nan", "inf", "neg_inf"])
def test_a_non_finite_orbit_coordinate_is_never_a_region_index(bad: float) -> None:
    """A NaN has no region and is no exterior: refused at the boundary.

    First red (measured 2026-10-05, before the fix): ``region_containing(nan)``
    on ``SPH`` returned 3 = ``n_regions``, an index that is neither a region
    nor an exterior code.
    """
    with pytest.raises(ValueError, match="a non-finite orbit coordinate is in no region"):
        SPH.region_containing(np.array([0.5, bad]))


@pytest.mark.foundation
@pytest.mark.parametrize("part", [SPH, CYL], ids=["sphere", "cylinder"])
def test_a_negative_orbit_coordinate_is_refused_where_the_coordinate_is_a_distance(part: ConcentricPartition) -> None:
    """A distance is never negative: ``region_containing(-0.3)`` is refused on the cylinder and the sphere.

    Before the restructure it returned the inner exterior, a cavity a solid body
    does not have (qa F6). The slab admits negatives
    (``test_the_slab_bare_locator_is_inner_owns_on_the_closed_domain``); ``-0.0``
    is not negative and is region 0.
    """
    with pytest.raises(ValueError, match="an orbit coordinate that is a distance is never negative"):
        part.region_containing(np.array([0.5, -0.3]))
    np.testing.assert_array_equal(part.region_containing(np.array([-0.0])), [0])


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_crossings_carry_the_region_they_enter")
def test_an_exterior_code_cannot_index_a_per_region_table() -> None:
    """The exterior codes are out of range for a table with one entry per region (ruled 2026-10-05; qa F1).

    qa's input: a hollow sphere (0.5, 1, 2), the line through the centre along z,
    ``sigma_t = [1, 2]``. ``sigma_t[slot_region]`` raises ``IndexError`` (the
    cavity slot's code is ``n = 2``). Extended explicitly with zeros for the two
    exteriors, the optical depth is ``2*1 + 1*0.5 + 0*1 + 1*0.5 + 2*1 = 5.0``;
    the pre-restructure codes (-1, -2) read the cavity as region 1 and gave 7.0.
    First red: exterior codes inside ``[-n, n)`` (any negative code indexes silently).
    """
    part = _partition(CoordSystem.SPHERICAL, (0.5, 1.0, 2.0))
    ch = _chord(part, (0.0, 0.0, -0.3), (0.0, 0.0, 1.0))
    sigma_t = np.array([1.0, 2.0])
    assert part.inner_exterior == part.n_regions == 2 and part.outer_exterior == 3
    with pytest.raises(IndexError):
        sigma_t[ch.slot_region]
    extended = np.append(sigma_t, [0.0, 0.0])
    tau = float(np.sum(extended[ch.slot_region] * ch.slot_length))
    assert tau == pytest.approx(5.0, rel=8 * _EPS)


@pytest.mark.l0
@pytest.mark.verifies("geometry-chord-segment-lengths")
@pytest.mark.rests_on(_HERE + "test_a_half_line_keeps_the_slots_beyond_its_start")
@pytest.mark.parametrize("part", [SPH_H, CYL_H], ids=["sphere", "cylinder"])
def test_a_half_line_starting_in_the_cavity(part: ConcentricPartition) -> None:
    """C5, qa gap Q2: a half-line from inside a hollow body's cavity, both orientations.

    From ``(0, 0.1, 0)`` on ``(0.4, 1.1, 2.0)`` (b = 0): along ``+e_y`` the cavity
    for 0.3, region 0 for 0.7, region 1 for 0.9; along ``-e_y`` the cavity for 0.5,
    then 0.7 and 0.9. The start is off the closest approach, so the two
    orientations cut the cavity slot differently.
    First red: the cavity slot's start (the inbound crossing of r_0) with its
    sign flipped (qa's arm Q2: 0 rows red before this gate).
    """
    p = np.array([0.0, 0.1, 0.0])
    for omega, expected in (((0.0, 1.0, 0.0), [0.0, 0.0, 0.3, 0.7, 0.9]), ((0.0, -1.0, 0.0), [0.0, 0.0, 0.5, 0.7, 0.9])):
        line = Line.through(p, np.array(omega))
        ch = part.chord(line)
        np.testing.assert_array_equal(ch.slot_region, [1, 0, I2, 0, 1])
        np.testing.assert_allclose(ch.lengths_beyond(line.parameter_of(p)), expected, rtol=0, atol=8 * _EPS)


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_a_half_line_keeps_the_slots_beyond_its_start",
                      _HERE + "test_a_line_parallel_to_the_orbit_space_is_inside_one_region_or_on_an_interface")
def test_lengths_beyond_on_a_batch_mixing_parallel_and_crossing_lines() -> None:
    """qa gap Q9 and F2: ``lengths_beyond`` on one batch holding parallel and crossing lines, warnings as errors.

    Cylinder ``CYL``, starts at each line's base point. Row 0: along ``e_z`` at
    rho = 0.7, parallel inside region 1: unbounded beyond its start, one infinite
    slot. Row 1: along ``-e_z`` on r_2 = 1.1: on an interface, every slot 0.
    Row 2: from ``(0, 0.7, 0)`` along ``-e_y``: 0.4, 0.3 + 0.3, 0.8, 0.9 as in the
    solid half-line row. Row 3: along ``e_y`` from the same point: 0.4, 0.9.
    First reds: a parallel line's half-line returned as 0 (qa's arm Q9); a
    ``-inf + inf`` evaluated before masking (a RuntimeWarning, raised here).
    """
    p = np.array([[0.7, 0.0, 2.0], [1.1, 0.0, 2.0], [0.0, 0.7, 0.0], [0.0, 0.7, 0.0]])
    om = np.array([[0.0, 0.0, 1.0], [0.0, 0.0, -1.0], [0.0, -1.0, 0.0], [0.0, 1.0, 0.0]])
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        line = Line.through(p, om)
        ch = CYL.chord(line)
        got = ch.lengths_beyond(line.parameter_of(p))
    # Row 0: one infinite slot, in region 1 (which of region 1's two slots holds it is not the contract).
    assert int(np.isinf(got[0]).sum()) == 1 and int(ch.slot_region[0][np.isinf(got[0])][0]) == 1
    assert np.all(got[0][np.isfinite(got[0])] == 0.0)
    np.testing.assert_array_equal(got[1], np.zeros(7))
    expected = np.array([[0.0, 0.4, 0.3, 0.0, 0.3, 0.8, 0.9], [0.0, 0.0, 0.0, 0.0, 0.0, 0.4, 0.9]])
    np.testing.assert_allclose(got[2:], expected, rtol=0, atol=8 * _EPS)


@pytest.mark.foundation
@pytest.mark.parametrize("part", [SPH, SLB], ids=["radial", "axial"])
def test_every_crossings_field_has_the_batch_shape(part: ConcentricPartition) -> None:
    """qa gap (d): every ``Crossings`` field is ``(..., m)``, for batch shapes (), (4,) and (2, 3).

    ``m = 2(n + 1)`` on a radial image, ``n + 1`` on an axial one; ``n_regions``
    is the partition's. First red: a field left ``(m,)`` and broadcast only
    at the consumer.
    """
    rng = np.random.default_rng(20261005)
    m = 2 * part.n_regions + 2 if part.chart.acts_on_kept_space else part.n_regions + 1
    for shape in ((), (4,), (2, 3)):
        om = rng.normal(size=(*shape, 3))
        om /= np.linalg.norm(om, axis=-1, keepdims=True)
        c = part.chord(Line.through(rng.uniform(-1, 1, (*shape, 3)), om)).crossings
        assert c.n_regions == part.n_regions
        for name in ("parameter", "breakpoint", "sense", "present", "region_entered"):
            assert np.shape(getattr(c, name)) == (*shape, m), name


@pytest.mark.foundation
@pytest.mark.parametrize("omega_x", [5e-324, 1e-310], ids=["smallest_subnormal", "subnormal_1e-310"])
def test_a_subnormal_projected_speed_lifts_to_infinite_lengths_never_nan(omega_x: float) -> None:
    """At a subnormal ``|P Omega|`` the obliquity overflows to infinity: a traversed slot
    is infinitely long, an untraversed one stays 0 (never ``0 * inf = nan``), the impact
    parameter is finite, and no crossing parameter is NaN (never ``-inf + inf``). First
    reds, measured 2026-10-05: scaling the orbit lengths by a bare product (5e-324: an
    untraversed slot NaN); the closest approach as ``t* = -(foot . P Omega)/|P Omega|^2``
    then ``foot + t* P Omega`` (1e-310: ``-inf * 0``, b NaN; the archivist's input).
    """
    part = ConcentricPartition(Chart(CoordSystem.CYLINDRICAL), (0.0, 1.0, 2.0))
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        chord = part.chord(Line.through(np.array([0.5, 0.0, 0.0]), np.array([omega_x, 0.0, 1.0])))
        lengths, starts = chord.slot_length, chord.slot_start
        parameters = chord.crossings.parameter
    assert chord.image.impact_parameter == 0.0
    assert not np.isnan(lengths).any() and not np.isnan(starts).any() and not np.isnan(parameters).any()
    np.testing.assert_array_equal(lengths, [np.inf, np.inf, 0.0, np.inf, np.inf])
