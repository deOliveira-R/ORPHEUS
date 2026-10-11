r"""Gates for the kernel's rays read from points, :meth:`ConcentricPartition.rays_from`, and the chord verbs the ray
layers read: ``line_integral``, ``cosine_at`` and :func:`region_half_lengths` (#578, P2 fix round).

The references are built here in mpmath (40 digits) from the ray's own algebra: the point ``c e_x`` and the unit
direction ``Omega`` give ``|P(c e_x + s Omega)|^2 = r^2``, a quadratic in the distance ``s`` (the slab's linear
``(r - c)/Omega_x``), whose roots are the crossings; a piece between two crossings is attributed to the region holding
its midpoint. Nothing here reads the kernel's level half-chords or slot bookkeeping (X4, derivation axis).

First reds: ``scratch/characteristic_architecture/p2/ta/battery`` (``run.sh "<arms>" <files>``), named per row.
"""
from __future__ import annotations

from typing import Any

import mpmath as mp
import numpy as np
import pytest

from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition, region_half_lengths
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.line import Line

_EPS = np.finfo(float).eps
_HERE = "tests/gates/geometry/test_rays.py::"
_CHORD = "tests/gates/geometry/test_chord.py::"
_BETWEEN = "tests/gates/geometry/test_chord_lengths_between.py::"

_SOLID = (0.0, 0.6, 1.3, 2.0)
_HOLLOW = (0.25, 0.6, 1.3, 2.0)
_SLAB = (-0.7, 0.3, 1.1, 2.0)
_CASES = [(CoordSystem.SPHERICAL, _SOLID), (CoordSystem.SPHERICAL, _HOLLOW), (CoordSystem.CYLINDRICAL, _SOLID),
          (CoordSystem.CYLINDRICAL, _HOLLOW), (CoordSystem.CARTESIAN, _SLAB)]
_IDS = ["sphere", "sphere-hollow", "cylinder", "cylinder-hollow", "slab"]


def _kept(coord: CoordSystem) -> int:
    return {CoordSystem.SPHERICAL: 3, CoordSystem.CYLINDRICAL: 2, CoordSystem.CARTESIAN: 1}[coord]


def _walk_mp(coord: CoordSystem, bps, c: float, omega) -> tuple[list, Any, int | None]:
    """The oracle: the pieces ``(start, end, region)`` of the ray from ``c e_x`` along ``omega`` until it leaves the
    domain ``[r_0, r_n]`` (region an index ``0 .. n-1``), the distance to that exit, and the wall it leaves through
    (the breakpoint index: 0 the cavity or the slab's lower face, n the outer surface; ``None`` for a ray that never
    enters, as from a point outside heading away)."""
    mp.mp.dps = 40
    k = _kept(coord)
    P = [mp.mpf(float(c)), mp.mpf(0), mp.mpf(0)]
    D = [mp.mpf(float(v)) for v in omega]
    r = [mp.mpf(float(x)) for x in bps]
    if coord is CoordSystem.CARTESIAN:
        roots = [(rk - P[0]) / D[0] for rk in r] if D[0] != 0 else []
        pos = (lambda s: P[0] + s * D[0])
    else:
        a = sum(D[j] ** 2 for j in range(k))
        b = 2 * sum(P[j] * D[j] for j in range(k))
        c0 = sum(P[j] ** 2 for j in range(k))
        roots = []
        for rk in r:
            disc = b * b - 4 * a * (c0 - rk * rk)
            if a > 0 and disc > 0:
                roots += [(-b - mp.sqrt(disc)) / (2 * a), (-b + mp.sqrt(disc)) / (2 * a)]
        pos = (lambda s: mp.sqrt(c0 + b * s + a * s * s))
    cuts = sorted({s for s in roots if s > 0})
    pieces: list = []
    start = mp.mpf(0)
    for end in [*cuts, (cuts[-1] if cuts else mp.mpf(0)) + 10]:
        m = (start + end) / 2
        p = pos(m)
        if p < r[0] or p > r[-1]:
            wall = 0 if p < r[0] else len(r) - 1
            return pieces, start, (wall if pieces or start > 0 or True else None)
        region = next(j for j in range(len(r) - 1) if p <= r[j + 1])
        pieces.append((start, end, region))
        start = end
    return pieces, start, None


def _directions(coord: CoordSystem, rng, n: int) -> np.ndarray:
    d = rng.normal(size=(n, 3))
    if coord is CoordSystem.CARTESIAN:
        d[:, 0] = np.where(np.abs(d[:, 0]) < 1e-3, 1e-3, d[:, 0])
    return d / np.linalg.norm(d, axis=1)[:, None]


def _points(bps) -> list[float]:
    """Points on every breakpoint and inside every region (never in the cavity)."""
    inside = [(lo + hi) / 2 for lo, hi in zip(bps[:-1], bps[1:])]
    return sorted({*map(float, bps), *inside} - ({0.0} if bps[0] == 0.0 else set()) | ({0.0} if bps[0] == 0.0 else set()))


# ── rays_from, ahead, to_wall, wall ─────────────────────────────────────────


@pytest.mark.l0
@pytest.mark.verifies("geometry-chord-segment-lengths")
@pytest.mark.rests_on(_BETWEEN + "test_the_lengths_between_two_parameters_are_the_roots_pieces",
                      _CHORD + "test_the_multi_region_segments_are_the_hand_written_table")
@pytest.mark.parametrize(("coord", "bps"), _CASES, ids=_IDS)
def test_the_run_to_the_first_wall_and_the_wall_are_the_oracles(coord: CoordSystem, bps) -> None:
    """``to_wall`` per region and its total, and ``wall``, against the oracle's walk to the first wall, from points on
    every breakpoint and inside every region, along 40 random directions each (seed 578).

    Band per region ``64 eps (r_n + 1) / |P Omega|`` (the parameters' scale lifted by the obliquity); the wall
    exactly. First reds: ``to_wall`` taken as ``ahead`` (the cavity and the far side counted); the first transit
    chosen by its opening crossing instead of its closing one (a point on a wall facing outward put in the next
    transit); ``wall`` read from the transit's entry wall.
    """
    rng = np.random.default_rng(578)
    part = ConcentricPartition(Chart(coord), bps)
    n = part.n_regions
    for c in _points(bps):
        omega = _directions(coord, rng, 40)
        rays = part.rays_from(np.full(40, c), omega)
        speed = np.maximum(Chart(coord).projected_speed(omega), 1e-300)
        for i in range(40):
            pieces, exit_, wall = _walk_mp(coord, bps, c, omega[i])
            band = 64 * _EPS * (abs(bps[-1]) + 1.0) / speed[i]
            for region in range(n):
                want = float(sum(e - s for s, e, reg in pieces if reg == region))
                got = float(np.sum(rays.to_wall[i][rays.chord.slot_region[i] == region]))
                assert abs(got - want) <= band, (coord, bps, c, omega[i], region, got, want)
            want_wall = part.no_interface if (wall is None or not pieces) else wall
            if not pieces and wall is not None:
                want_wall = wall                                  # on a wall facing out: in the transit ending there
            assert int(rays.wall[i]) == want_wall, (coord, bps, c, omega[i], int(rays.wall[i]), want_wall)


@pytest.mark.l0
@pytest.mark.rests_on(_HERE + "test_the_run_to_the_first_wall_and_the_wall_are_the_oracles")
@pytest.mark.parametrize(("coord", "bps"), _CASES, ids=_IDS)
def test_a_point_on_a_wall_sits_on_its_crossing(coord: CoordSystem, bps) -> None:
    """On every breakpoint, inward and outward and (curvilinear) tangent: the point's parameter is EXACTLY the
    parameter of its crossing of that breakpoint (formed from the same level, #590), so a point on a wall facing out
    of the domain leaves at once (``to_wall`` 0, ``wall`` that breakpoint) and one facing in has nothing behind it.
    `[M]` 2026-10-10, before the parameter was formed from the level: a 1-ulp cavity slot ahead of an observer on
    the cavity wall looking outward. First red: the point's parameter from ``line.parameter_of``.
    """
    part = ConcentricPartition(Chart(coord), bps)
    n = part.n_regions
    for k, r in enumerate(bps):
        if coord is not CoordSystem.CARTESIAN and r == 0.0:
            continue                                                # the centre is a stratum, never a wall
        mus = (0.6, 0.40267903915837233, 0.9793148738097697, 0.9606546472234416)   # the last three: `[M]` 2026-10-10 reds
        dirs = [np.array([m, np.sqrt((1 - m) * (1 + m)), 0.0]) for m in mus]
        rays = part.rays_from(np.full(2 * len(mus), r), np.stack([*dirs, *(d * [-1, 1, 1] for d in dirs)]))
        q, bp, present = rays.chord.crossings.parameter, rays.chord.crossings.breakpoint, rays.chord.crossings.present
        for i in range(2 * len(mus)):
            on = present[i] & (bp[i] == k)
            assert np.any(rays.start[i] == q[i][on]), (coord, r, i, rays.start[i], q[i][on])
        leaves = (k == n and True) or (k == 0)                      # outward at the outer wall, inward at the inner one
        out, inward = slice(0, len(mus)), slice(len(mus), 2 * len(mus))
        if k == n:
            assert np.all(rays.to_wall[out].sum(axis=-1) == 0.0) and np.all(rays.wall[out] == n)
            assert np.all(rays.ahead[out].sum(axis=-1) == 0.0)
        if k == 0 and bps[0] > 0.0 or (coord is CoordSystem.CARTESIAN and k == 0):
            assert np.all(rays.to_wall[inward].sum(axis=-1) == 0.0) and np.all(rays.wall[inward] == 0)
        if k == 0 and bps[0] > 0.0 and coord is not CoordSystem.CARTESIAN:
            assert np.all(rays.wall[out] == n), "a point on the cavity wall looking outward sees the cavity"
        del leaves
    if coord is not CoordSystem.CARTESIAN:
        r = bps[-1]
        for b in [x for x in bps[1:-1] if x > 0.0]:
            t = -float(np.sqrt(1.0 - (b / r) ** 2))
            tangent = part.rays_from(np.array(r), np.array([t, float(np.sqrt(1 - t * t)), 0.0]))
            assert np.isfinite(tangent.ahead).all() and float(tangent.ahead.sum()) > 0.0


@pytest.mark.l0
@pytest.mark.rests_on(_HERE + "test_a_point_on_a_wall_sits_on_its_crossing")
def test_a_slab_point_on_a_face_sits_on_its_crossing() -> None:
    """The slab's twin of :func:`test_a_point_on_a_wall_sits_on_its_crossing`, over 2000 random 3-D directions per
    face (seed 1): the point's parameter equals its crossing's bit for bit, and a point facing out of the slab leaves
    at once through that face. Found 2026-10-11 writing these gates: 1306 of 8000 differed by a rounding (the axial
    image added its shift to the coordinate, the crossings to the parameter), and an outward point on face 0 read
    ``wall = no_interface``; fixed the same day (``AxialImage.parameter_origin``). First red: the point's parameter
    formed through ``line.parameter_of`` (the curvilinear twin of this arm reds its row above)."""
    rng = np.random.default_rng(1)
    part = ConcentricPartition(Chart(CoordSystem.CARTESIAN), _SLAB)
    for k, c in enumerate(_SLAB):
        omega = _directions(CoordSystem.CARTESIAN, rng, 2000)
        rays = part.rays_from(np.full(2000, c), omega)
        on = rays.chord.crossings.present & (rays.chord.crossings.breakpoint == k)
        hits = np.any((rays.chord.crossings.parameter == rays.start[:, None]) & on, axis=-1)
        assert np.all(hits), f"face {k}: {int(np.sum(~hits))} of 2000 points off their crossing"
        outward = (omega[:, 0] < 0) if k == 0 else (omega[:, 0] > 0) if k == len(_SLAB) - 1 else np.zeros(2000, bool)
        assert np.all(rays.wall[outward] == k)


@pytest.mark.foundation
@pytest.mark.parametrize(("coord", "bps"), _CASES, ids=_IDS)
def test_the_orbit_coordinate_at_the_point_is_the_points(coord: CoordSystem, bps) -> None:
    """On the cylinder and the sphere ``Rays.orbit_coordinate_at(0)`` is the point's coordinate bit for bit (formed
    from the point's own position ``+-y* + |P Omega| l``), on every breakpoint and inside every region, along 200
    directions (`[M]` the elegance review: the hypot form exceeded R at the point in 2123 of 20001); along the run
    ahead it exceeds the outer radius by at most one rounding (`[M]` 2026-10-11: 4.4e-16 at the exit end). On the
    slab it is ``c + Omega_x l``, exact at 0 too (until 2026-10-11 it read through the parameter: 4 of 200 off by 1 ulp).
    First red: the coordinate read through the parameter, ``chord.orbit_coordinate_at(start + l)``.
    """
    rng = np.random.default_rng(579)
    part = ConcentricPartition(Chart(coord), bps)
    for c in _points(bps):
        omega = _directions(coord, rng, 200)
        rays = part.rays_from(np.full(200, c), omega)
        at0 = rays.orbit_coordinate_at(np.zeros((200, 1)))[..., 0]
        np.testing.assert_array_equal(at0, np.full(200, c))
        if coord is CoordSystem.CARTESIAN:
            continue
        length = rays.ahead.sum(axis=-1)
        along = rays.orbit_coordinate_at(length[:, None] * np.linspace(0.0, 1.0, 9))
        assert np.all(along <= bps[-1] + 2 * np.spacing(bps[-1])), float(np.max(along) - bps[-1])


# ── line_integral ───────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_a_per_region_field_of_the_wrong_shape_is_refused() -> None:
    """``line_integral`` takes one value per region, shape ``(n,)``: one short (it would empty the outermost region),
    one long (it would fill the cavity) and a column ``(n, 1)`` are refused with "one value per region".
    First red: the shape check removed (`[M]` the elegance review, before the verb: one short read 2.856 for 5.914).
    """
    part = ConcentricPartition(Chart(CoordSystem.SPHERICAL), _HOLLOW)
    rays = part.rays_from(np.array([1.0]), np.array([[-0.6, 0.8, 0.0]]))
    for bad in (np.ones(2), np.ones(4), np.ones((3, 1))):
        with pytest.raises(ValueError, match="one value per region"):
            rays.chord.line_integral(bad, rays.ahead)
    assert float(rays.chord.line_integral(np.ones(3), rays.ahead)[0]) > 0.0


@pytest.mark.l0
@pytest.mark.rests_on(_HERE + "test_the_run_to_the_first_wall_and_the_wall_are_the_oracles")
def test_the_line_integral_is_the_oracles_and_broadcasts_over_runs() -> None:
    """``sum f_region l`` over ``ahead`` against the oracle's pieces (exteriors void), and over a ``(k, S)`` stack of runs
    per line (``lengths_between`` to k distances) equal to the k one-run calls bit for bit.
    First reds: the cavity given the innermost region's value; the run axis broadcast against the lines' axis.
    """
    rng = np.random.default_rng(580)
    f = np.array([0.7, 3.1, 0.2])
    for coord in (CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL):
        part = ConcentricPartition(Chart(coord), _HOLLOW)
        omega = _directions(coord, rng, 30)
        rays = part.rays_from(np.full(30, 1.0), omega)
        got = rays.chord.line_integral(f, rays.ahead)
        for i in range(30):
            mp.mp.dps = 40
            k = _kept(coord)
            a = sum(mp.mpf(float(v)) ** 2 for v in omega[i][:k])
            b = 2 * mp.mpf(float(omega[i][0]))
            want = mp.mpf(0)
            roots = []
            for rk in _HOLLOW:
                disc = b * b - 4 * a * (1 - mp.mpf(rk) ** 2)
                if disc > 0:
                    roots += [(-b - mp.sqrt(disc)) / (2 * a), (-b + mp.sqrt(disc)) / (2 * a)]
            cuts = sorted({s for s in roots if s > 0})
            start = mp.mpf(0)
            for end in cuts:
                p = mp.sqrt(1 + b * (start + end) / 2 + a * ((start + end) / 2) ** 2)
                if _HOLLOW[0] < p < _HOLLOW[-1]:
                    want += mp.mpf(float(f[next(j for j in range(3) if p <= _HOLLOW[j + 1])])) * (end - start)
                start = end
            speed = float(np.linalg.norm(omega[i][:k]))
            assert abs(float(got[i]) - float(want)) <= 64 * _EPS * 4.0 * f.max() / max(speed, 1e-300), (coord, i)
        dist = np.array([0.2, 0.9, 1.7, 3.5])
        stack = np.stack([rays.lengths_to(np.full(30, x)) for x in dist], axis=-2)        # (30, 4, S)
        many = rays.chord.line_integral(f, stack)
        assert many.shape == (30, 4)
        for j, x in enumerate(dist):
            np.testing.assert_array_equal(many[:, j], rays.chord.line_integral(f, rays.lengths_to(np.full(30, x))))


@pytest.mark.foundation
def test_a_void_slot_of_infinite_length_adds_nothing() -> None:
    """A line parallel to the cylinder's axis inside a void region has one infinite slot: the integral of a field that
    is 0 there is exactly 0 (never ``0 * inf = nan``), and of one that is not, +inf.
    First red: the product formed unmasked (``(length * value).sum()``).
    """
    part = ConcentricPartition(Chart(CoordSystem.CYLINDRICAL), _SOLID)
    line = Line.through(np.array([1.0, 0.0, 0.0]), np.array([0.0, 0.0, 1.0]))        # region 1, parallel
    chord = part.chord(line)
    assert np.isinf(chord.slot_length).any()
    with np.errstate(all="raise"):
        assert float(chord.line_integral(np.array([0.7, 0.0, 0.2]), chord.slot_length)) == 0.0
        assert float(chord.line_integral(np.array([0.7, 3.1, 0.2]), chord.slot_length)) == np.inf


# ── cosine_at, region_half_lengths ──────────────────────────────────────────


@pytest.mark.l0
@pytest.mark.parametrize("coord", [CoordSystem.SPHERICAL, CoordSystem.CYLINDRICAL, CoordSystem.CARTESIAN],
                         ids=["sphere", "cylinder", "slab"])
def test_the_cosine_at_a_surface_is_the_angle_to_its_normal(coord: CoordSystem) -> None:
    """``cosine_at(r)`` is ``|Omega . n|`` at the line's crossing of ``c = r``, ``n`` the surface's normal there
    (radial in the kept coordinates on the sphere and the cylinder, ``e_x`` on the slab), computed here from the
    crossing point in mpmath; 0 where the line misses the surface. 3-D directions, so ``|P Omega| < 1`` on the cylinder.
    First reds: the cosine without the projected speed (``h / r``, right only at ``|P Omega| = 1``); the sine ``b / r``.
    """
    rng = np.random.default_rng(581)
    part = ConcentricPartition(Chart(coord), _SLAB if coord is CoordSystem.CARTESIAN else _SOLID)
    omega = _directions(coord, rng, 200)
    rays = part.rays_from(np.full(200, 1.0), omega)
    radius = 1.3
    got = rays.chord.cosine_at(radius)
    k = _kept(coord)
    for i in range(200):
        if coord is CoordSystem.CARTESIAN:
            want = abs(float(omega[i][0]))
        else:
            mp.mp.dps = 40
            D = [mp.mpf(float(v)) for v in omega[i]]
            a = sum(D[j] ** 2 for j in range(k))
            disc = (2 * D[0]) ** 2 - 4 * a * (1 - mp.mpf(radius) ** 2)
            if disc <= 0:
                assert got[i] == 0.0
                continue
            s = (-2 * D[0] + mp.sqrt(disc)) / (2 * a)
            x = [(1 if j == 0 else 0) + s * D[j] for j in range(k)]
            want = float(abs(sum(x[j] * D[j] for j in range(k))) / mp.mpf(radius))
        assert abs(float(got[i]) - want) <= 64 * _EPS, (coord, i, float(got[i]), want)


@pytest.mark.l0
@pytest.mark.verifies("geometry-chord-segment-lengths")
def test_the_region_half_lengths_are_the_differences_of_half_chords() -> None:
    """``region_half_lengths(r, h)`` against ``sqrt(r_{k+1}^2 - b^2) - sqrt(r_k^2 - b^2)`` in mpmath (the region of
    closest approach its outer half-chord; 0 where unreached), within 8 eps relative on every nonzero entry, over
    2000 impact parameters plus near-tangent ones (r_k (1 - 10^-u), u in [1, 15]) on a thin-shell partition.
    First red: the shells formed as the difference ``h_{k+1} - h_k``.
    """
    mp.mp.dps = 50
    rng = np.random.default_rng(582)
    r = np.array([0.0, 1.0, 1.0 + 1e-6, 2.0])
    b = np.concatenate([rng.uniform(0, 2, 2000), (r[1:] * (1 - 10.0 ** -rng.uniform(1, 15, (60, 3)))).ravel()])
    h = np.sqrt(np.where(b[:, None] < r, (r - b[:, None]) * (r + b[:, None]), 0.0))
    got = region_half_lengths(r, h)
    worst = 0.0
    for i, bb in enumerate(b):
        B = mp.mpf(float(bb))
        H = [mp.sqrt(mp.mpf(float(x)) ** 2 - B * B) if B < x else mp.mpf(0) for x in r]
        for k in range(3):
            want = H[k + 1] - H[k] if B < r[k] else (H[k + 1] if B < r[k + 1] else mp.mpf(0))
            if want == 0:
                assert got[i, k] == 0.0
            else:
                worst = max(worst, float(abs(mp.mpf(float(got[i, k])) - want) / want))
    assert worst <= 8 * _EPS, worst
