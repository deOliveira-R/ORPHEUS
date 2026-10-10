"""Gates for the line part of the characteristic closure (:class:`~orpheus.derivations.continuous.characteristic.LinePeriod`).

P1 step (b), first rung, of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "P1 step (b),
first rung: API sketch", item 2). Verification spec
``scratch/characteristic_architecture/p1_verification_spec.md``, rows A1, A3,
A6, A2 (its optical-depth half), B1, B3, B4a (closure leg), B5b (closure
legs), re-specified onto the built names, plus the near-lossless row.

The ladder, bottom up (``rests_on`` on each row):

1. the kernel's transits (``tests/gates/geometry/test_chord_transits.py``) and
   the walls (``test_characteristic_walls.py``);
2. the period: a hand-counted table of (rank, traversals, exit walls,
   amplitudes) [A1, A6], and the period against the PHYSICAL unfolded path,
   marched by reflecting (or translating) the 3-D line at each wall in the
   test [A3, the structural reference for the period's combinatorial rule];
3. the optical depth of each traversal against closed-form segment lengths
   computed in mpmath [A2];
4. the cycle's least solution against the unfolded wall-by-wall sum on
   seeded abstract data [B1], its exact edges [B3, B4a, B5b] and the
   near-lossless line [1 - P without cancellation].

The albedo pairing (spec §0): the inflow to traversal k + 1 carries the
amplitude of the wall traversal k EXITS at, which is the wall at which the
backward path from k + 1 reflects.

Rows at ``l0`` with ``verifies`` (labels on
``docs/theory/references/characteristic.rst``): ``characteristic-transit-rank``
for ``test_the_period_matches_the_hand_counted_table`` and
``test_the_period_is_the_physically_unfolded_path``;
``characteristic-closure`` for ``test_the_inflow_is_the_unfolded_wall_by_wall_sum``,
``test_an_absorbing_wall_zeroes_exactly_the_inflow_it_feeds``,
``test_every_amplitude_zero_adds_nothing``,
``test_the_least_solution_on_a_lossless_trapped_line`` and
``test_a_nearly_lossless_line_keeps_its_digits``. The three optical-depth
functions (12 rows) are ``l0`` with ``verifies("characteristic-traversal-integrals")``
since rung 2 (`[M]` 2026-10-06, their witnesses: the depth scaled by 1 + 1e-12,
the crossing mask dropped, and a void given a floor of Sigma; together they
redden 12 of 12, ``scratch/characteristic_architecture/p1_step_b2/battery``
arms O1-O3). Every other row is ``foundation``.
"""
from __future__ import annotations

from typing import Any

import mpmath as mp
import numpy as np
import pytest

from orpheus.derivations.continuous.characteristic import LinePeriod, TrappedSource, Walls
from orpheus.geometry.boundary import (
    BC,
    AlbedoBoundary,
    PeriodicBoundary,
    ReflectiveBoundary,
    SpecularReturn,
    VacuumInflow,
)
from orpheus.geometry.chart import RadialImage
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.line import Line
from orpheus.geometry.structured_geometry import StructuredGeometry

_HERE = "tests/gates/derivations/test_characteristic_closure.py::"
_WALLS = "tests/gates/derivations/test_characteristic_walls.py::"
_TRANSITS = "tests/gates/geometry/test_chord_transits.py::"
_CHORD = "tests/gates/geometry/test_chord.py::"

_EPS = float(np.finfo(float).eps)


def _spec(a: float):
    return AlbedoBoundary(a, SpecularReturn(axis="x"))


_SOLID_BP = (0.0, 0.3, 1.1, 2.0)
_HOLLOW_BP = (0.4, 1.1, 2.0)
_SLAB_BP = (-0.7, 0.3, 1.1, 2.0)
_A_IN, _A_OUT = 0.3, 0.6           # distinct amplitudes, so a pairing swap moves every reading
_BELOW_R0 = float(np.nextafter(0.4, 0.0))
_BELOW_R = float(np.nextafter(2.0, 0.0))


def _body(coord: str, *, inner: Any = None, outer: Any = None, left: Any = None, right: Any = None) -> StructuredGeometry:
    match coord:
        case "sphere_solid":
            return StructuredGeometry.sphere(_SOLID_BP, (0, 1, 2), outer=outer)
        case "cylinder_solid":
            return StructuredGeometry.cylinder(_SOLID_BP, (0, 1, 2), outer=outer)
        case "sphere_hollow":
            return StructuredGeometry.sphere(_HOLLOW_BP, (0, 1), inner=inner, outer=outer)
        case "cylinder_hollow":
            return StructuredGeometry.cylinder(_HOLLOW_BP, (0, 1), inner=inner, outer=outer)
        case "slab":
            return StructuredGeometry.slab(_SLAB_BP, (0, 1, 2), left=left, right=right)
    raise AssertionError(coord)


def _period(geometry: StructuredGeometry, points, directions) -> LinePeriod:
    line = Line.through(np.asarray(points, dtype=float), np.asarray(directions, dtype=float))
    return LinePeriod.of(ConcentricPartition.of(geometry).chord(line), Walls.of(geometry))


def _in_plane(b: float, wz: float):
    """The point (b, 0, 0) and a direction with in-plane part along y: the impact parameter is b bit for bit."""
    return (b, 0.0, 0.0), (0.0, float(np.sqrt((1.0 - wz) * (1.0 + wz))), wz)


# ── 2. the period: the hand-counted table [A1, A6] ───────────────────────
#
# Per row: the traversals (transit, reversed, exit wall, amplitude), in order.
_SOLID = [(0, False, 3, _A_OUT)]
_HOLLOW_ONE = [(0, False, 2, _A_OUT)]
_HOLLOW_TWO = [(0, False, 0, _A_IN), (1, False, 2, _A_OUT)]   # first leg ends on the inner wall


def _table():
    rows = []
    for wz in (0.0, 0.8, -0.6):
        charts = [("sphere", 0.0)] if wz == 0.0 else []
        charts.append(("cylinder", wz))
        for chart, z in charts:
            solid = f"{chart}_solid"
            for b, expected in ((0.0, _SOLID), (0.2, _SOLID), (0.3, _SOLID), (0.7, _SOLID), (1.1, _SOLID),
                                (1.5, _SOLID), (_BELOW_R, _SOLID), (2.0, []), (2.5, [])):
                if b == _BELOW_R and z != 0.0:
                    continue      # a tilted line's Plucker foot rounds b up to R (test_chord_transits' measured case)
                rows.append((f"{solid}_b{b}_wz{z}", _body(solid, outer=_spec(_A_OUT)), *_in_plane(b, z), expected))
            hollow = f"{chart}_hollow"
            for b, expected in ((0.0, _HOLLOW_TWO), (0.2, _HOLLOW_TWO), (_BELOW_R0, _HOLLOW_TWO),
                                (0.4, _HOLLOW_ONE), (0.7, _HOLLOW_ONE), (1.1, _HOLLOW_ONE), (1.5, _HOLLOW_ONE),
                                (2.0, [])):
                geometry = _body(hollow, inner=_spec(_A_IN), outer=_spec(_A_OUT))
                rows.append((f"{hollow}_b{b}_wz{z}", geometry, *_in_plane(b, z), expected))
    mirrors = _body("slab", left=_spec(_A_IN), right=_spec(_A_OUT))
    rows.append(("slab_mirrors_rising", mirrors, (0.5, 0.0, 0.0), (0.6, 0.8, 0.0),
                 [(0, False, 3, _A_OUT), (0, True, 0, _A_IN)]))
    rows.append(("slab_mirrors_falling", mirrors, (0.5, 0.0, 0.0), (-0.3, 0.0, float(np.sqrt(0.91))),
                 [(0, False, 0, _A_IN), (0, True, 3, _A_OUT)]))
    periodic = _body("slab", left=PeriodicBoundary(axis="x"), right=BC("periodic"))
    rows.append(("slab_periodic_rising", periodic, (0.5, 0.0, 0.0), (0.6, 0.8, 0.0), [(0, False, 3, 1.0)]))
    rows.append(("slab_periodic_falling", periodic, (0.5, 0.0, 0.0), (-0.3, 0.0, float(np.sqrt(0.91))),
                 [(0, False, 0, 1.0)]))
    rows.append(("slab_vacuum_left_mirror_right", _body("slab", left=VacuumInflow(), right=BC.reflective),
                 (0.5, 0.0, 0.0), (0.6, 0.8, 0.0), [(0, False, 3, 1.0), (0, True, 0, 0.0)]))
    # Parallel lines: no transit, rank 0, amplitude 0.
    rows.append(("cylinder_axial", _body("cylinder_solid", outer=_spec(_A_OUT)), (0.7, 0.0, 0.0), (0.0, 0.0, 1.0), []))
    rows.append(("slab_parallel", mirrors, (0.5, 0.0, 0.0), (0.0, 0.6, 0.8), []))
    return rows


_TABLE = _table()


@pytest.mark.l0
@pytest.mark.verifies("characteristic-transit-rank", "peierls-greens-hollow-sph-impact-parameter-partition",
                      "peierls-greens-annulus-impact-parameter-partition")
@pytest.mark.parametrize(("geometry", "point", "direction", "expected"), [r[1:] for r in _TABLE],
                         ids=[r[0] for r in _TABLE])
@pytest.mark.rests_on(_TRANSITS + "test_the_transits_match_the_hand_counted_walls",
                      _WALLS + "test_every_law_reads_as_the_wall_its_physics_names")
def test_the_period_matches_the_hand_counted_table(geometry, point, direction, expected) -> None:
    """[A1, A6] Rank, traversals, exit walls and amplitudes of each line's period, counted by hand.

    The rank is never tagged: 1 on a solid body, on a shell's ray missing the
    cavity (b = r_0 included: a tangency is not a crossing) and on the periodic
    slab; 2 through the cavity (one ulp inside r_0 included) and between slab
    mirrors; 0 on a parallel line and on a tangent or missing line. The
    periodic slab continues WITHOUT reversal. First reds: (a) rank 2 on every
    ray of a hollow body; (b) the periodic face read as a mirror (rank 2, the
    transit reversed); (c) the amplitude of the ENTRY wall (the pairing
    swapped); (d) the reversed candidate preferred over the forward one (the
    convention the sketch ruled).
    """
    period = _period(geometry, point, direction)
    m = len(expected)
    assert int(period.rank) == m
    np.testing.assert_array_equal(period.present, [k < m for k in range(2)])
    got = [(int(period.transit[k]), bool(period.reversed[k]), int(period.exit_wall[k]), float(period.amplitude[k]))
           for k in range(m)]
    assert got == expected, f"{got} != {expected}"
    np.testing.assert_array_equal(period.amplitude[m:], 0.0)
    n_plus_1 = len(geometry.breakpoints)
    np.testing.assert_array_equal(period.exit_wall[m:], n_plus_1)
    np.testing.assert_array_equal(period.entry_wall[m:], n_plus_1)
    if m:
        # The cycle closes: each traversal enters at the partner of the wall the previous one exits.
        walls = Walls.of(geometry)
        exits = period.exit_wall[:m]
        np.testing.assert_array_equal(period.entry_wall[:m], walls.partner_at(np.roll(exits, 1)))


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_period_matches_the_hand_counted_table")
def test_a_batch_of_lines_reads_each_line_as_it_reads_alone() -> None:
    """The period of a batch equals the per-line periods, for a batch mixing ranks 0, 1 and 2."""
    geometry = _body("sphere_hollow", inner=_spec(_A_IN), outer=_spec(_A_OUT))
    bs = [0.2, 0.7, 2.5, 0.0, _BELOW_R0, 0.4]
    batch = _period(geometry, [(b, 0.0, 0.0) for b in bs], [(0.0, 1.0, 0.0)] * len(bs))
    for i, b in enumerate(bs):
        alone = _period(geometry, (b, 0.0, 0.0), (0.0, 1.0, 0.0))
        for field in ("transit", "reversed", "amplitude", "present", "exit_wall", "entry_wall"):
            np.testing.assert_array_equal(getattr(batch, field)[i], getattr(alone, field), err_msg=f"{field} b={b}")


# ── 2. the period against the physically unfolded path [A3] ──────────────


def _radial_normal(coord: str, x: np.ndarray) -> np.ndarray:
    if coord.startswith("sphere"):
        return x / np.linalg.norm(x)
    if coord.startswith("cylinder"):
        p = np.array([x[0], x[1], 0.0])
        return p / np.linalg.norm(p)
    return np.array([1.0, 0.0, 0.0])


def _traversed(chord, transit: int):
    """(regions, lengths) of the traversed slots of one transit of a single-line chord, in line order."""
    t = chord.transits
    first, stop = int(t.first_slot[transit]), int(t.stop_slot[transit])
    length = chord.slot_length[first:stop]
    keep = length > 0.0
    return chord.slot_region[first:stop][keep], length[keep]


def _march(coord: str, geometry: StructuredGeometry, point, direction, steps: int):
    """The physically unfolded path: at each exit, reflect the 3-D line by Householder (or translate it across a wrap).

    Returns, per step, the (regions, lengths, exit wall, line) of the transit
    the path runs along. No ``LinePeriod`` and no wall table on this side:
    the partner of a wrap is the opposite face, written here.
    """
    partition = ConcentricPartition.of(geometry)
    periodic = isinstance(geometry.boundaries[0], PeriodicBoundary)
    width = geometry.breakpoints[-1] - geometry.breakpoints[0]
    line = Line.through(np.asarray(point, dtype=float), np.asarray(direction, dtype=float))
    start = None                                   # the parameter where the current transit begins
    out = []
    for _ in range(steps):
        chord = partition.chord(line)
        t = chord.transits
        present = [j for j in range(2) if bool(t.present[j])]
        if start is None:
            j = present[0]
        else:
            starts = [float(chord.crossings.parameter[int(t.first_slot[i])]) for i in present]
            j = present[int(np.argmin([abs(s - start) for s in starts]))]
        regions, lengths = _traversed(chord, j)
        exit_wall = int(t.exit_wall[j])
        out.append((regions, lengths, exit_wall, line))
        t_exit = float(chord.crossings.parameter[int(t.stop_slot[j])])
        x = line.at(np.array(t_exit))
        omega = np.asarray(line.direction, dtype=float)
        if periodic:
            x = x - np.sign(omega[0]) * np.array([width, 0.0, 0.0])
        else:
            n = _radial_normal(coord, x)
            omega = omega - 2.0 * float(omega @ n) * n
        line = Line.through(x, omega)
        start = float(line.parameter_of(x))
    return out


_MARCH = [  # (id, chart, geometry, point, direction)
    ("sphere_hollow_through_cavity", "sphere_hollow", _body("sphere_hollow", inner=_spec(_A_IN), outer=_spec(_A_OUT)),
     (0.23, 0.0, 0.0), (0.0, 1.0, 0.0)),
    ("sphere_hollow_shell_ray", "sphere_hollow", _body("sphere_hollow", inner=_spec(_A_IN), outer=_spec(_A_OUT)),
     (0.0, 0.0, 0.0), (0.0, 0.0, 1.0)),       # through the centre: b = 0
    ("sphere_solid_region1", "sphere_solid", _body("sphere_solid", outer=_spec(_A_OUT)),
     (0.13, 0.62, -0.2), (0.48, 0.0, float(np.sqrt((1.0 - 0.48) * (1.0 + 0.48))))),
    ("cylinder_hollow_oblique", "cylinder_hollow", _body("cylinder_hollow", inner=_spec(_A_IN), outer=_spec(_A_OUT)),
     (0.17, 0.05, 0.0), (0.36, 0.48, 0.8)),
    ("cylinder_hollow_oblique_shell", "cylinder_hollow",
     _body("cylinder_hollow", inner=_spec(_A_IN), outer=_spec(_A_OUT)), (0.83, 0.1, 0.0), (0.0, 0.8, -0.6)),
    ("slab_mirrors_rising", "slab", _body("slab", left=_spec(_A_IN), right=_spec(_A_OUT)),
     (0.5, 0.0, 0.0), (0.6, 0.8, 0.0)),
    ("slab_mirrors_falling", "slab", _body("slab", left=_spec(_A_IN), right=_spec(_A_OUT)),
     (0.5, 0.0, 0.0), (-0.3, 0.0, float(np.sqrt(0.91)))),
    ("slab_periodic_rising", "slab", _body("slab", left=PeriodicBoundary(axis="x"), right=PeriodicBoundary(axis="x")),
     (0.5, 0.0, 0.0), (0.6, 0.8, 0.0)),
    ("slab_periodic_falling", "slab", _body("slab", left=PeriodicBoundary(axis="x"), right=PeriodicBoundary(axis="x")),
     (0.5, 0.0, 0.0), (-0.3, 0.0, float(np.sqrt(0.91)))),
]


def _unit(v):
    v = np.asarray(v, dtype=float)
    return v / np.linalg.norm(v)


@pytest.mark.l0
@pytest.mark.verifies("characteristic-transit-rank")
@pytest.mark.parametrize(("chart", "geometry", "point", "direction"), [r[1:] for r in _MARCH],
                         ids=[r[0] for r in _MARCH])
@pytest.mark.rests_on(_HERE + "test_the_period_matches_the_hand_counted_table",
                      _CHORD + "test_the_chord_is_invariant_under_the_symmetry_group")
def test_the_period_is_the_physically_unfolded_path(chart, geometry, point, direction) -> None:
    """[A3, and the reference for the period] Each traversal is the transit the reflected 3-D line runs along.

    The test marches the line: at each exit it reflects the direction by
    ``Omega - 2 (Omega . n) n`` with the wall normal it builds (radial on the
    sphere, its in-plane part on the cylinder, ``e_x`` on the slab), or
    translates the exit point across a periodic slab. Per step it asserts the
    reflection invariants (A3: the impact parameter and ``|P Omega|`` kept,
    ``Omega_z`` kept on the cylinder, ``Omega_x`` negated by a slab mirror and
    kept by a wrap), and that the traversed (region, length) sequence and the
    exit wall equal the period's traversal read off the ORIGINAL chord
    (reversed where the period says so). After the rank's steps the march is
    back on traversal 0. Lengths to 64 ulp relative (the reflected line's
    impact parameter re-derived through a rounded exit point).

    First reds: the slab mirror's traversal 1 read forward (the region order
    0, 1, 2 against 2, 1, 0); the periodic slab read as a mirror (reversed, and
    ``Omega_x`` kept by the march); the hollow rank collapsed to 1. Preferring the
    reversed candidate over the forward one also reds here (`[M]` battery arm
    C4: 5 of 9 rows, the radial ones), through the rank it changes, not through
    the per-slot order.
    """
    period = _period(geometry, point, direction)
    m = int(period.rank)
    assert m >= 1
    original = ConcentricPartition.of(geometry).chord(
        Line.through(np.asarray(point, dtype=float), np.asarray(direction, dtype=float)))
    steps = _march(chart, geometry, point, _unit(direction), m + 1)
    image0 = original.image
    for k, (regions, lengths, exit_wall, line) in enumerate(steps):
        slot = k % m
        want_regions, want_lengths = _traversed(original, int(period.transit[slot]))
        if bool(period.reversed[slot]):
            want_regions, want_lengths = want_regions[::-1], want_lengths[::-1]
        np.testing.assert_array_equal(regions, want_regions, err_msg=f"step {k}: regions")
        np.testing.assert_allclose(lengths, want_lengths, rtol=64 * _EPS, atol=0.0, err_msg=f"step {k}: lengths")
        assert exit_wall == int(period.exit_wall[slot]), (k, exit_wall, period.exit_wall)
        # A3: the invariants the reflection keeps.
        chord = ConcentricPartition.of(geometry).chord(line)
        np.testing.assert_allclose(chord.projected_speed, original.projected_speed, rtol=8 * _EPS)
        omega, omega0 = np.asarray(line.direction), _unit(direction)
        if chart == "slab":
            sign = 1.0 if isinstance(geometry.boundaries[0], PeriodicBoundary) or k % 2 == 0 else -1.0
            np.testing.assert_allclose(omega[0], sign * omega0[0], rtol=4 * _EPS)
        else:
            assert isinstance(chord.image, RadialImage) and isinstance(image0, RadialImage)
            np.testing.assert_allclose(chord.image.impact_parameter, image0.impact_parameter,
                                       rtol=64 * _EPS, atol=64 * _EPS * geometry.breakpoints[-1])
            if chart.startswith("cylinder"):
                np.testing.assert_allclose(omega[2], omega0[2], rtol=4 * _EPS, atol=4 * _EPS)


# ── 3. the optical depth of each traversal [A2, the tau half] ────────────

mp.mp.dps = 40


def _h(r: float, b: float):
    return mp.sqrt(max(mp.mpf(r) ** 2 - mp.mpf(b) ** 2, 0))


def _radial_depths(bp, sigma, b: float, speed: float, rank: int):
    """Closed-form optical depths of the transits of a radial line at impact parameter b (mpmath)."""
    n = len(bp) - 1
    half = mp.mpf(0)                       # one side of the closest approach, from the outer wall inward
    for j in range(n):
        lo, hi = bp[j], bp[j + 1]
        if b >= hi:
            continue
        inner = _h(lo, b) if b < lo else mp.mpf(0)
        half += mp.mpf(sigma[j]) * (_h(hi, b) - inner)
    half /= mp.mpf(speed)
    return [half, half] if rank == 2 else [2 * half]


_SIGMA_SOLID = (0.7, 1.9, 0.4)
_SIGMA_HOLLOW = (1.3, 0.0)                 # an outer void shell: its slots contribute exactly 0
_SIGMA_SLAB = (0.5, 1.7, 0.9)


_DEPTH_ROWS = [  # (id, chart, b, wz, rank)
    ("sphere_solid_b0", "sphere_solid", 0.0, 0.0, 1),
    ("sphere_solid_b0.2_region0", "sphere_solid", 0.2, 0.0, 1),
    ("sphere_solid_b0.7_region1", "sphere_solid", 0.7, 0.0, 1),
    ("sphere_solid_b1.5_outer_only", "sphere_solid", 1.5, 0.0, 1),
    ("cylinder_solid_b0.7_oblique", "cylinder_solid", 0.7, 0.8, 1),
    ("sphere_hollow_b0.2_through_cavity", "sphere_hollow", 0.2, 0.0, 2),
    ("sphere_hollow_b0.7_shell", "sphere_hollow", 0.7, 0.0, 1),
    ("cylinder_hollow_b0.2_through_cavity_oblique", "cylinder_hollow", 0.2, -0.6, 2),
    ("cylinder_hollow_b1.5_void_shell_only", "cylinder_hollow", 1.5, 0.8, 1),
]


@pytest.mark.l0
@pytest.mark.verifies("characteristic-traversal-integrals")
@pytest.mark.parametrize(("chart", "b", "wz", "rank"), [r[1:] for r in _DEPTH_ROWS], ids=[r[0] for r in _DEPTH_ROWS])
@pytest.mark.rests_on(_HERE + "test_the_period_matches_the_hand_counted_table",
                      _CHORD + "test_every_slot_length_matches_the_closed_form")
def test_the_optical_depth_of_each_traversal_is_the_closed_form(chart: str, b: float, wz: float, rank: int) -> None:
    """[A2, tau] tau_k = sum over the traversal's slots of Sigma_t(region) x length, against mpmath segment lengths.

    Multi-region with distinct Sigma_t, a void shell (Sigma_t = 0), oblique
    cylinder lines (lengths over |P Omega|), and rank 2 through the cavity,
    where each transit's tau is ONE side only. 8 ulp relative. First reds: (a)
    the sum taken over every slot of the chord instead of the transit's (a
    rank-2 tau doubled; green on rank 1, declared); (b) Sigma_t read in
    reversed region order (red on the multi-region rows); (c) the exteriors
    given a non-zero Sigma_t (green: the cavity slot lies outside every
    transit, declared).
    """
    sigma = _SIGMA_HOLLOW if "hollow" in chart else _SIGMA_SOLID
    bp = _HOLLOW_BP if "hollow" in chart else _SOLID_BP
    outer = _spec(_A_OUT)
    geometry = _body(chart, inner=_spec(_A_IN) if "hollow" in chart else None, outer=outer)
    period = _period(geometry, *_in_plane(b, wz))
    assert int(period.rank) == rank
    speed = float(np.sqrt((1.0 - wz) * (1.0 + wz)))
    want = [float(v) for v in _radial_depths(bp, sigma, b, speed, rank)]
    tau = period.optical_depth(np.array(sigma))
    assert tau.shape == (2,)
    np.testing.assert_allclose(tau[:rank], want, rtol=8 * _EPS, atol=0.0)
    np.testing.assert_array_equal(tau[rank:], 0.0)


@pytest.mark.l0
@pytest.mark.verifies("characteristic-traversal-integrals")
@pytest.mark.parametrize("direction", [(0.6, 0.8, 0.0), (-0.3, 0.0, float(np.sqrt(0.91)))], ids=["rising", "falling"])
@pytest.mark.rests_on(_HERE + "test_the_period_matches_the_hand_counted_table")
def test_the_slab_optical_depth_is_the_widths_over_the_cosine(direction) -> None:
    """[A2, tau] Each slab traversal (forward and reversed) crosses every region once: sum Sigma_j w_j / |Omega_x|."""
    geometry = _body("slab", left=_spec(_A_IN), right=_spec(_A_OUT))
    period = _period(geometry, (0.5, 0.0, 0.0), direction)
    widths = np.diff(_SLAB_BP)
    want = float(sum(mp.mpf(s) * mp.mpf(w) for s, w in zip(_SIGMA_SLAB, widths)) / abs(mp.mpf(direction[0])))
    np.testing.assert_allclose(period.optical_depth(np.array(_SIGMA_SLAB)), [want, want], rtol=8 * _EPS, atol=0.0)


@pytest.mark.l0
@pytest.mark.verifies("characteristic-traversal-integrals")
def test_a_parallel_line_has_no_optical_depth() -> None:
    """A rank-0 line's tau is 0 in both columns, never the infinite slot's inf x Sigma or a NaN."""
    period = _period(_body("cylinder_solid", outer=_spec(_A_OUT)), (0.7, 0.0, 0.0), (0.0, 0.0, 1.0))
    np.testing.assert_array_equal(period.optical_depth(np.array(_SIGMA_SOLID)), [0.0, 0.0])


# ── 4. the cycle's least solution [B1] ───────────────────────────────────


def _unfolded(amplitude, tau, outflow, m: int):
    """The wall-by-wall sum, mpmath: in_k = sum_{j>=1} a_{k-j} B_{k-j} prod_{i=1}^{j-1} a_{k-i} e^{-tau_{k-i}}.

    Marches the backward path one wall at a time (indices mod m) and stops
    when a term falls below 1e-30 of the partial sum; no geometric-series
    division on this side.
    """
    a = [mp.mpf(float(v)) for v in amplitude[:m]]
    g = [a[i] * mp.exp(-mp.mpf(float(tau[i]))) for i in range(m)]
    out = []
    for k in range(m):
        total, weight, j = mp.mpf(0), mp.mpf(1), 1
        while True:
            i = (k - j) % m
            term = weight * a[i] * mp.mpf(float(outflow[i]))
            total += term
            weight *= g[i]
            j += 1
            if weight == 0 or (j > 2 and abs(weight) < mp.mpf("1e-30") * max(abs(total), mp.mpf("1e-300"))):
                break
        out.append(total)
    return out


def _period_with(amplitudes, m: int, batch: int) -> LinePeriod:
    """A period of rank m on a batch of identical lines, its amplitudes set by the walls' laws.

    m = 1: a solid sphere's line (amplitude a_out); m = 2: a hollow sphere's
    through-cavity line (amplitudes a_in, a_out in traversal order).
    """
    if m == 1:
        geometry = _body("sphere_solid", outer=_spec(amplitudes[0]))
        b = 0.7
    else:
        geometry = _body("sphere_hollow", inner=_spec(amplitudes[0]), outer=_spec(amplitudes[1]))
        b = 0.2
    period = _period(geometry, [(b, 0.0, 0.0)] * batch, [(0.0, 1.0, 0.0)] * batch)
    np.testing.assert_array_equal(period.amplitude[0, :m], amplitudes[:m])
    return period


_AMPLITUDES = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (0.3, 0.6), (0.6, 0.3), (1.0, 0.85), (0.47, 1.0), (0.99, 0.12),
               (0.5, 0.5)]


@pytest.mark.l0
@pytest.mark.verifies("characteristic-closure")
@pytest.mark.catches("ERR-035")
@pytest.mark.parametrize("m", [1, 2])
@pytest.mark.parametrize("amplitudes", _AMPLITUDES, ids=[f"a{a0}_{a1}" for a0, a1 in _AMPLITUDES])
@pytest.mark.rests_on(_HERE + "test_the_period_matches_the_hand_counted_table")
def test_the_inflow_is_the_unfolded_wall_by_wall_sum(amplitudes, m: int) -> None:
    """[B1] The closed cycle equals the unfolded backward path, on seeded abstract (tau, B), 16 ulp relative.

    40 seeded draws per amplitude pair: tau in [0, 3] with exact zeros mixed
    in, B > 0 with a trailing basis axis of 3. The cycle product is kept
    below 0.9 (the unfolded sum converges; the lossless and near-lossless
    edges are their own rows). Distinct amplitudes, so the pairing is
    visible. First reds: (a) the ERR-035 denominator (a_0^2 e^{-2 tau_0} for
    the cycle product); (b) the once-around cross term dropped; (c) the
    pairing swapped (the inflow to k scaled by a_k).

    ERR-035 itself (the symmetric slab's closure by analogy with the sphere:
    the out-and-back integral with no amplitude at the inner reflection) reds
    8 of the 9 rank-2 rows, every pair but the vacuum (0, 0); the equal pair
    (0.5, 0.5) is its own regime (``[M]`` 2026-10-10, arm ``err035`` of
    ``scratch/characteristic_architecture/p1_step_e/ta_e1b/battery``). The
    rank-1 rows stay green under it: the old rank-1 closure was honest.
    """
    rng = np.random.default_rng(20261006 + 7 * m + int(100 * amplitudes[0]) + int(1000 * amplitudes[1]))
    draws = 40
    period = _period_with(np.array(amplitudes), m, draws)
    tau = rng.uniform(0.0, 3.0, size=(draws, 2))
    tau[rng.uniform(size=(draws, 2)) < 0.2] = 0.0
    a = np.array(amplitudes[:m])
    product = np.prod(a) * np.exp(-tau[:, :m].sum(axis=1))
    tau[product > 0.9, :m] += 0.2                 # keep the cycle product below 0.9 (product <= e^{-0.2 m} there)
    outflow = rng.uniform(0.1, 2.0, size=(draws, 2, 3))
    got = period.inflow(tau, outflow)
    assert got.shape == (draws, 2, 3)
    np.testing.assert_array_equal(got[:, m:], 0.0)
    for d in range(draws):
        for c in range(3):
            want = [float(v) for v in _unfolded(period.amplitude[d], tau[d], outflow[d, :, c], m)]
            scale = max(max(abs(w) for w in want), 1e-300)
            np.testing.assert_allclose(got[d, :m, c], want, rtol=16 * _EPS, atol=16 * _EPS * scale,
                                       err_msg=f"draw {d}, column {c}")
    # The trailing axis is a broadcast: one column alone reads the same.
    np.testing.assert_array_equal(period.inflow(tau, outflow[..., 1]), got[..., 1])


# ── 4. the exact edges [B3, B4a, B5b] ────────────────────────────────────


@pytest.mark.l0
@pytest.mark.verifies("characteristic-closure")
@pytest.mark.parametrize("case", ["slab_left_absorbing", "slab_right_absorbing", "hollow_inner_absorbing"])
@pytest.mark.rests_on(_HERE + "test_the_inflow_is_the_unfolded_wall_by_wall_sum")
def test_an_absorbing_wall_zeroes_exactly_the_inflow_it_feeds(case: str) -> None:
    """[B3] With the wall the backward path reaches first absorbing (amplitude 0), that traversal's inflow is 0 bitwise.

    The other traversal's inflow is then exactly its predecessor's outflow
    (the absorbing wall cuts the cycle: 1 - P = 1). Slab with vacuum left,
    mirror right (and the mirror image); hollow with an absorbing inner wall
    and a mirror outside. First red: the pairing swapped (the amplitude of the
    wall the traversal ENTERS at), which puts the mirror's 1 on the cut.
    """
    match case:
        case "slab_left_absorbing":
            period = _period(_body("slab", left=VacuumInflow(), right=ReflectiveBoundary(axis="x")),
                             (0.5, 0.0, 0.0), (0.6, 0.8, 0.0))
            absorbing = 0
        case "slab_right_absorbing":
            period = _period(_body("slab", left=ReflectiveBoundary(axis="x"), right=BC.vacuum),
                             (0.5, 0.0, 0.0), (0.6, 0.8, 0.0))
            absorbing = 3
        case "hollow_inner_absorbing":
            period = _period(_body("sphere_hollow", inner=BC.vacuum, outer=BC.reflective),
                             (0.2, 0.0, 0.0), (0.0, 1.0, 0.0))
            absorbing = 0
        case _:
            raise AssertionError(case)
    assert int(period.rank) == 2
    cut = int(np.flatnonzero(period.exit_wall == absorbing)[0])    # the traversal exiting at the absorbing wall
    fed = (cut + 1) % 2                                            # the traversal whose inflow it feeds
    tau = np.array([0.37, 1.21])
    outflow = np.array([[0.83, 2.5], [1.9, 0.04]])                 # (2 traversals, 2 basis columns)
    got = period.inflow(tau, outflow)
    np.testing.assert_array_equal(got[fed], [0.0, 0.0])
    np.testing.assert_array_equal(got[cut], outflow[fed])          # amplitude 1, nothing returns around


@pytest.mark.l0
@pytest.mark.verifies("characteristic-closure")
@pytest.mark.parametrize("geometry", [
    _body("slab", left=VacuumInflow(), right=AlbedoBoundary(0.0)),
    _body("sphere_hollow", inner=BC.vacuum, outer=BC.vacuum),
    _body("sphere_solid", outer=VacuumInflow()),
], ids=["slab", "sphere_hollow", "sphere_solid"])
@pytest.mark.rests_on(_HERE + "test_the_inflow_is_the_unfolded_wall_by_wall_sum",
                      _WALLS + "test_every_law_reads_as_the_wall_its_physics_names")
def test_every_amplitude_zero_adds_nothing(geometry) -> None:
    """[B4a, closure leg] Vacuum everywhere: the period's amplitudes are 0 and the inflow is 0 bitwise, finite.

    A batch of lines of every rank on the body (the absorbing walls still
    define the period). First reds: a vacuum wall read with a default
    amplitude 1; the amplitude left out of the returned outflow.
    """
    bs = [0.0, 0.2, 0.7, 1.5, 2.5]
    direction = (0.6, 0.8, 0.0) if geometry.coord.name == "CARTESIAN" else (0.0, 1.0, 0.0)
    period = _period(geometry, [(b, 0.0, 0.0) for b in bs], [direction] * len(bs))
    np.testing.assert_array_equal(period.amplitude, 0.0)
    rng = np.random.default_rng(4)
    got = period.inflow(rng.uniform(0.0, 3.0, size=(len(bs), 2)), rng.uniform(0.1, 2.0, size=(len(bs), 2, 4)))
    np.testing.assert_array_equal(got, 0.0)


@pytest.mark.l0
@pytest.mark.verifies("characteristic-closure")
@pytest.mark.rests_on(_HERE + "test_the_inflow_is_the_unfolded_wall_by_wall_sum",
                      _HERE + "test_the_optical_depth_of_each_traversal_is_the_closed_form")
def test_the_least_solution_on_a_lossless_trapped_line() -> None:
    """[B5b, closure legs] A trapped line with no source carries 0 exactly; with a source it is refused; else finite.

    The body: a solid sphere (0, 1, 2) with Sigma_t = (1, 0), an outer void
    shell, and a mirror outside. A line with b in (1, 2) runs in the void
    only: tau = 0 exactly from ``optical_depth``, amplitude 1, the cycle
    product 1. Batched with a line through the material (b = 0.5) and a
    two-traversal lossless line beside a lossy one (the refusal is the
    trapped line's, never the void's). First reds: (a) the closed form at
    1 - P = 0 (NaN on the trapped line); (b) the refusal removed (a source
    on the trapped line read as 0); (c) the refusal widened to any lossless
    traversal (the mixed m = 2 line refused).
    """
    geometry = StructuredGeometry.sphere((0.0, 1.0, 2.0), (0, 1), outer=BC.reflective)
    period = _period(geometry, [(0.5, 0.0, 0.0), (1.5, 0.0, 0.0)], [(0.0, 1.0, 0.0)] * 2)
    tau = period.optical_depth(np.array([1.0, 0.0]))
    assert tau[1, 0] == 0.0 and tau[0, 0] > 0.0
    source_free = np.array([[[0.7], [0.0]], [[0.0], [0.0]]])
    got = period.inflow(tau, source_free)
    np.testing.assert_array_equal(got[1], 0.0)
    assert np.all(np.isfinite(got)) and got[0, 0, 0] > 0.0
    with pytest.raises(TrappedSource, match="lossless trapped line"):
        period.inflow(tau, np.array([[[0.7], [0.0]], [[1e-3], [0.0]]]))
    # A two-traversal cycle lossless on one side only is not trapped: finite and positive.
    mirrors = _period(_body("slab", left=BC.reflective, right=BC.reflective), (0.5, 0.0, 0.0), (0.6, 0.8, 0.0))
    mixed = mirrors.inflow(np.array([0.0, 0.5]), np.array([1.0, 1.0]))
    assert np.all(np.isfinite(mixed)) and np.all(mixed > 0.0)
    # Both sides lossless: trapped again, 0 without a source, refused with one.
    np.testing.assert_array_equal(mirrors.inflow(np.array([0.0, 0.0]), np.array([0.0, 0.0])), 0.0)
    with pytest.raises(TrappedSource, match="lossless trapped line"):
        mirrors.inflow(np.array([0.0, 0.0]), np.array([0.0, 2.0]))


@pytest.mark.l0
@pytest.mark.verifies("characteristic-closure")
@pytest.mark.parametrize(("amplitudes", "tau"), [
    ((1.0, 1.0), (1e-12, 2.5e-12)),
    ((1.0 - 1e-13, 1.0), (3e-13, 4e-13)),
    ((1.0, 1.0), (1e-12, None)),
], ids=["m2_lossless_walls", "m2_wall_loss_dominant", "m1"])
@pytest.mark.rests_on(_HERE + "test_the_inflow_is_the_unfolded_wall_by_wall_sum")
def test_a_nearly_lossless_line_keeps_its_digits(amplitudes, tau) -> None:
    """1 - P is formed without cancellation: tau ~ 1e-12 against the exact closed form (mpmath), 4 ulp.

    First red: ``1 - exp(sum log a - sum tau)`` in place of ``-expm1(...)``,
    which loses about 5 digits on 1 - P (``[M]`` relative error 2.2e-5 at tau = 1e-12,
    6.3e-6 on this row's period; both far above the 4-ulp band).
    """
    m = 1 if tau[1] is None else 2
    period = _period_with(np.array(amplitudes), m, 1)
    tau_arr = np.array([tau[0], tau[1] if m == 2 else 0.0])[None, :]
    outflow = np.array([[1.3, 0.7]])
    got = period.inflow(tau_arr, outflow)[0]
    a = [mp.mpf(float(v)) for v in period.amplitude[0, :m]]
    g = [a[i] * mp.exp(-mp.mpf(float(tau_arr[0, i]))) for i in range(m)]
    one_minus = 1 - (g[0] * g[1] if m == 2 else g[0])
    b = [mp.mpf(float(v)) for v in outflow[0, :m]]
    if m == 1:
        want = [a[0] * b[0] / one_minus]
    else:
        want = [(a[1] * b[1] + g[1] * a[0] * b[0]) / one_minus, (a[0] * b[0] + g[0] * a[1] * b[1]) / one_minus]
    np.testing.assert_allclose(got[:m], [float(w) for w in want], rtol=4 * _EPS, atol=0.0)


@pytest.mark.foundation
@pytest.mark.parametrize("length", [2, 4], ids=["n_minus_1", "n_plus_1"])
def test_a_sigma_t_of_the_wrong_length_is_refused(length: int) -> None:
    """``optical_depth`` takes one total cross section per region; a wrong length never broadcasts or truncates."""
    period = _period(_body("slab", left=BC.reflective, right=BC.reflective), (0.5, 0.0, 0.0), (0.6, 0.8, 0.0))
    with pytest.raises(ValueError, match="one total cross section per region"):
        period.optical_depth(np.linspace(0.3, 1.2, length))


@pytest.mark.foundation
@pytest.mark.rests_on(_HERE + "test_the_slab_optical_depth_is_the_widths_over_the_cosine")
def test_a_subnormal_cosine_through_a_void_reads_infinite_depth_and_finite_inflow() -> None:
    """Omega_x = 5e-324 through the slab (0.5, 0, 0.9): every slot is infinitely long.

    tau = +inf exactly on both traversals (the true depth, never a NaN from
    the void's inf x 0), so the inflow is the one-bounce return, a_{k-1} B_{k-1},
    finite. With every region void tau = 0 exactly and the line is a trapped
    one. The leg ``no RuntimeWarning from closure.py``: the masked product must
    not be evaluated at inf x 0 (`[M]` 2026-10-06 it was: ``invalid value in
    multiply``). The kernel's own overflow warning (``chord.py``) is outside
    this rung and filtered out.
    """
    import warnings

    geometry = _body("slab", left=BC.reflective, right=BC.reflective)
    period = _period(geometry, (0.5, 0.0, 0.0), (5e-324, 1.0, 0.0))
    with warnings.catch_warnings():
        warnings.filterwarnings("error", category=RuntimeWarning, module=r".*characteristic\.closure")
        tau = period.optical_depth(np.array([0.5, 0.0, 0.9]))
        got = period.inflow(tau, np.array([1.0, 2.0]))
        void = period.optical_depth(np.zeros(3))
        trapped = period.inflow(void, np.zeros(2))
    np.testing.assert_array_equal(tau, [np.inf, np.inf])
    np.testing.assert_array_equal(got, [2.0, 1.0])
    np.testing.assert_array_equal(void, [0.0, 0.0])
    np.testing.assert_array_equal(trapped, [0.0, 0.0])


@pytest.mark.foundation
def test_a_chord_whose_walls_disagree_with_its_transits_is_refused(monkeypatch: pytest.MonkeyPatch) -> None:
    """The two ``RuntimeError`` guards of ``LinePeriod.of``, each by its own fragment.

    "enters at the partner": reachable from valid inputs, a radial chord (one
    transit, entering and leaving at breakpoint 3) read through a periodic
    slab's walls (partner of 3 is 0, where nothing enters). "did not close":
    unreachable from a valid ``Walls`` (its invariants make the partner map an
    involution), so ``Walls.partner_at`` is monkeypatched (on the frozen class) to send both slab walls to
    breakpoint 3: the period runs 0 forward, 0 reversed, 0 reversed, and never
    returns to its first traversal.
    """
    periodic = Walls.of(_body("slab", left=BC("periodic"), right=BC("periodic")))
    radial = ConcentricPartition.of(_body("sphere_solid", outer=BC.reflective)).chord(
        Line.through(np.array([0.7, 0.0, 0.0]), np.array([0.0, 1.0, 0.0])))
    with pytest.raises(RuntimeError, match="enters at the partner"):
        LinePeriod.of(radial, periodic)
    mirrors = Walls.of(_body("slab", left=BC.reflective, right=BC.reflective))
    chord = ConcentricPartition.of(_body("slab", left=BC.reflective, right=BC.reflective)).chord(
        Line.through(np.array([0.5, 0.0, 0.0]), np.array([0.6, 0.8, 0.0])))
    assert int(LinePeriod.of(chord, mirrors).rank) == 2                 # the honest walls close
    monkeypatch.setattr(Walls, "partner_at", lambda self, bp, where=None: np.full(np.shape(bp), 3))
    with pytest.raises(RuntimeError, match="did not close"):
        LinePeriod.of(chord, mirrors)


# ── 7. an arriving flux (rung 3) [IN1-IN4] ───────────────────────────────
#
# P1 step (b), third rung (the plan's API sketch item 4; spec
# ``scratch/characteristic_architecture/p1_step_b3/spec.md`` rows IN1-IN4).
# ``arriving`` is s_k, a flux injected at traversal k's entry from outside the
# line part (a diffuse wall's re-entry), not multiplied by the wall's amplitude.
# The cycle in_k = a_{k-1} (e^{-tau_{k-1}} in_{k-1} + B_{k-1}) + s_k is solved here
# by hand, in mpmath: rank 1, in_0 = (a_0 B_0 + s_0) / (1 - g_0); rank 2,
# in_0 = (a_1 B_1 + s_0 + g_1 (a_0 B_0 + s_1)) / (1 - g_0 g_1), g = a e^{-tau}.


def _arriving_closed_form(amplitude, tau, outflow, arriving, m: int):
    a = [mp.mpf(float(v)) for v in amplitude[:m]]
    g = [a[i] * mp.exp(-mp.mpf(float(tau[i]))) for i in range(m)]
    B = [mp.mpf(float(v)) for v in outflow[:m]]
    s = [mp.mpf(float(v)) for v in arriving[:m]]
    if m == 1:
        return [(a[0] * B[0] + s[0]) / (1 - g[0])]
    return [(a[(k + 1) % 2] * B[(k + 1) % 2] + s[k] + g[(k + 1) % 2] * (a[k] * B[k] + s[(k + 1) % 2])) / (1 - g[0] * g[1])
            for k in range(2)]


def _arriving_periods():
    periodic = _body("slab", left=PeriodicBoundary(axis="x"), right=BC("periodic"))
    partial = _body("slab", left=_spec(0.3), right=_spec(0.8))
    return {
        "sphere_solid_a0.6": (_period_with((0.6,), 1, 20), 1),
        "slab_periodic_a1": (_period(periodic, [(0.5, 0.0, 0.0)] * 20, [(0.6, 0.8, 0.0)] * 20), 1),
        "sphere_hollow_cavity_a0.3_0.8": (_period_with((0.3, 0.8), 2, 20), 2),
        "slab_partial_mirrors_0.3_0.8": (_period(partial, [(0.5, 0.0, 0.0)] * 20, [(0.6, 0.8, 0.0)] * 20), 2),
    }


_ARRIVING = _arriving_periods()


@pytest.mark.l0
@pytest.mark.verifies("characteristic-closure")
@pytest.mark.parametrize("name", list(_ARRIVING))
@pytest.mark.rests_on(_HERE + "test_the_inflow_is_the_unfolded_wall_by_wall_sum")
def test_an_arriving_flux_enters_the_cycle_at_its_traversal(name: str) -> None:
    """[IN1, IN2] ``inflow(tau, B, arriving=s)`` is the cycle solved by hand, rank 1 and rank 2, 20 seeded draws, 8 ulp.

    A trailing basis axis of 5; s_0 != s_1 in every rank-2 draw. `[M]` 2026-10-06
    on the spec's prototype: 1.0, 1.6, 0.8, 0.9 ulp. First reds: (a) s_k injected at
    traversal k + 1's entry (the shift dropped: `[M]` 0.25 hollow, 0.46 slab, and
    BLIND at rank 1, where the shift is the identity: the rank-1 rows cannot see
    it, declared); (b) s multiplied by the wall's amplitude (`[M]` 0.27 sphere, 0.31
    hollow, 0.14 slab; BLIND on the periodic slab, amplitude 1, declared).
    """
    period, m = _ARRIVING[name]
    rng = np.random.default_rng(20261006 + m + len(name))
    draws = period.present.shape[0]
    tau = rng.uniform(0.05, 3.0, size=(draws, 2)) * period.present
    outflow = rng.uniform(0.1, 2.0, size=(draws, 2, 5)) * period.present[..., None]
    arriving = rng.uniform(0.1, 2.0, size=(draws, 2, 5)) * period.present[..., None]
    if m == 2:
        assert np.all(arriving[:, 0] != arriving[:, 1])
    got = period.inflow(tau, outflow, arriving)
    np.testing.assert_array_equal(got[:, m:], 0.0)
    for d in range(draws):
        for c in range(5):
            want = [float(v) for v in _arriving_closed_form(period.amplitude[d], tau[d], outflow[d, :, c], arriving[d, :, c], m)]
            np.testing.assert_allclose(got[d, :m, c], want, rtol=8 * _EPS, atol=0.0, err_msg=f"draw {d}, column {c}")


@pytest.mark.foundation
@pytest.mark.parametrize("name", list(_ARRIVING))
@pytest.mark.rests_on(_HERE + "test_an_arriving_flux_enters_the_cycle_at_its_traversal")
def test_a_zero_arriving_flux_is_no_arriving_flux_bitwise(name: str) -> None:
    """[IN3] ``arriving=0`` returns ``arriving=None``'s result bit for bit (a zero added in floating point is exact).

    ``arriving=None`` is gated against the rung-2 closed forms by every other
    inflow row of this file, unchanged. First red: ``arriving=None`` routed
    through an expression regrouped from the ``arriving`` path (the bits move).
    """
    period, _m = _ARRIVING[name]
    rng = np.random.default_rng(7)
    draws = period.present.shape[0]
    tau = rng.uniform(0.05, 3.0, size=(draws, 2)) * period.present
    outflow = rng.uniform(0.1, 2.0, size=(draws, 2, 3)) * period.present[..., None]
    np.testing.assert_array_equal(period.inflow(tau, outflow, np.zeros_like(outflow)), period.inflow(tau, outflow))


@pytest.mark.l0
@pytest.mark.verifies("characteristic-closure")
@pytest.mark.rests_on(_HERE + "test_the_least_solution_on_a_lossless_trapped_line")
def test_an_arriving_flux_on_a_lossless_trapped_line_is_refused() -> None:
    """[IN4] The ``TrappedSource`` refusal covers an arriving flux: a lossless trapped line with s != 0 and no source raises.

    With s = 0 there the inflow is 0, no raise; the material line beside it is
    finite. First red: the trapped test reading the outflow only (the arriving
    flux added after the test): the first leg returns inf or a silent 0.
    """
    geometry = StructuredGeometry.sphere((0.0, 1.0, 2.0), (0, 1), outer=BC.reflective)
    period = _period(geometry, [(0.5, 0.0, 0.0), (1.5, 0.0, 0.0)], [(0.0, 1.0, 0.0)] * 2)
    tau = period.optical_depth(np.array([1.0, 0.0]))
    no_source = np.zeros((2, 2, 1))
    with pytest.raises(TrappedSource, match="lossless trapped line"):
        period.inflow(tau, no_source, np.array([[[0.0], [0.0]], [[0.4], [0.0]]]))
    got = period.inflow(tau, no_source, np.array([[[0.4], [0.0]], [[0.0], [0.0]]]))
    np.testing.assert_array_equal(got[1], 0.0)
    assert np.isfinite(got[0, 0, 0]) and got[0, 0, 0] > 0.0
