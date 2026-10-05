"""L1 gates for the measure on lines: chord lengths integrated over the beam density.

The geometric kernel seed (``.claude/plans/characteristic_reference_architecture.md``,
P0; verification spec ``scratch/characteristic_architecture/seed_verification_spec.md``
§7b, rows M1-M4). The invariant measure on lines is ``dA_perp dOmega``; for one
direction it is a density over the impact parameter ``b >= 0``
(``Chart.beam_density``, ``S_{d-2} b^{d-2} |P Omega|``: ``2 pi b`` on the sphere,
``2 |P Omega|`` per unit height on the cylinder, the two lines at ``+-b``;
``|Omega_x|`` per unit area on the slab). Every integral runs over ``b`` in ``[0, R]``. Integrating a region's chord length
against it gives the region's measure; averaging over directions gives Cauchy's
mean chord ``4V/S``.

**The quadrature is the test's, not the kernel's** (ruled 2026-10-05: geometry
does not import ``derivations/``). A chord length has a square-root endpoint
at every radius, so plain Gauss-Legendre in ``b`` converges as ``n^-3``
(``[M]`` 4.1e-5 at 16 nodes per piece); the rule here is Gauss-Legendre in
``theta`` on each piece ``[r_{j-1}, r_j]`` with ``b = r_j sin(theta)``, which
absorbs that endpoint (``[M]`` 2.7e-16 at 16 nodes). Its tolerance is the
reduction depth: ``n * P * eps`` for ``n`` nodes on ``P`` pieces.

The references are typed closed forms (``pi (r_{j+1}^2 - r_j^2)``,
``(4/3) pi (r_{j+1}^3 - r_j^3)``, ``4V/S``); a second leg compares the areas and
volumes with ``Chart.measure``, which derives them from ``T(r) = r^d``, a
different derivation sharing no identity with a chord integral.
"""
from __future__ import annotations

import numpy as np
import pytest
from numpy.polynomial.legendre import leggauss

from orpheus.geometry.chart import Chart
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.coord import CoordSystem
from orpheus.geometry.line import Line

pytestmark = pytest.mark.l1

_EPS = np.finfo(float).eps
_N = 16
_CHORD = "tests/gates/geometry/test_chord.py::"
_CHART = "tests/gates/geometry/test_chart.py::"
_SOLID = (0.0, 0.3, 1.1, 2.0)
_HOLLOW = (0.4, 1.1, 2.0)


def _impact_rule(breakpoints) -> tuple[np.ndarray, np.ndarray]:
    """Nodes and weights on ``b`` in ``[0, R]``: per piece ``[lo, hi]``, Gauss-Legendre in ``theta`` with ``b = hi sin(theta)``.

    The pieces are ``[0, r_1], [r_1, r_2], ...`` (a hollow body adds ``[0, r_0]``).
    """
    x, w = leggauss(_N)
    edges = [0.0, *[r for r in breakpoints if r > 0.0]]
    nodes, weights = [], []
    for lo, hi in zip(edges, edges[1:]):
        a, b = np.arcsin(lo / hi), np.pi / 2
        theta = 0.5 * (b - a) * x + 0.5 * (b + a)
        nodes.append(hi * np.sin(theta))
        weights.append(0.5 * (b - a) * w * hi * np.cos(theta))
    return np.concatenate(nodes), np.concatenate(weights)


def _region_lengths(partition: ConcentricPartition, b: np.ndarray, omega: np.ndarray) -> np.ndarray:
    """The 3-D length of each region (and the cavity, last column) on the lines at impact parameters ``b``.

    The lines pass through ``(b, 0, 0)`` with direction ``omega``, whose in-plane
    part is along ``e_y``, so the impact parameter is ``b`` by construction.
    """
    line = Line.through(np.stack([b, 0 * b, 0 * b], axis=-1), np.broadcast_to(omega, (b.size, 3)))
    ch = partition.chord(line)
    n = partition.n_regions
    out = np.zeros((b.size, n + 1))
    for j in range(n):
        out[:, j] = np.sum(np.where(ch.slot_region == j, ch.slot_length, 0.0), axis=-1)
    out[:, n] = np.sum(np.where(ch.slot_region == partition.inner_exterior, ch.slot_length, 0.0), axis=-1)
    return out


def _tolerance(breakpoints) -> float:
    pieces = len([r for r in breakpoints if r > 0.0])
    return _N * pieces * _EPS


@pytest.mark.verifies("geometry-measure-on-lines")
@pytest.mark.rests_on(_CHORD + "test_every_slot_length_matches_the_closed_form",
                      _CHART + "test_the_beam_density_of_each_chart")
@pytest.mark.parametrize("bp", [_SOLID, _HOLLOW], ids=["solid", "hollow"])
@pytest.mark.parametrize(
    "omega", [(0.0, 1.0, 0.0), (0.0, 0.6, 0.8), (0.0, 0.8, -0.6), (0.0, 0.28, 0.96)],
    ids=["in_plane", "s0.6", "s0.8_down", "s0.28"],
)
def test_the_beam_integral_of_a_cylinder_region_is_its_area(bp, omega) -> None:
    """M1/M3: ``int_0^R l_j(b) 2 |P Omega| db`` is the area of annulus ``j``, for every direction.

    The 3-D chord carries ``1/|P Omega|`` and the beam density ``|P Omega|``; they
    cancel only if both are right, so the oblique directions (``|P Omega|`` 0.6,
    0.8, 0.28, with ``Omega_z`` of both signs) test the obliquity against the
    measure. The cavity of the hollow body integrates to ``pi r_0^2``.
    First reds: the obliquity multiplied; the density without ``|P Omega|`` (every
    oblique row); the density without ``S_0 = 2`` (every row: the beam over ``b >= 0``
    counts the two lines at ``+-b``).
    """
    part = ConcentricPartition(Chart(CoordSystem.CYLINDRICAL), bp)
    om = np.array(omega)
    b, w = _impact_rule(bp)
    density = part.chart.beam_density(b, om)
    got = np.sum((w * density)[:, None] * _region_lengths(part, b, om), axis=0)
    r = np.asarray(bp)
    areas = np.pi * (r[1:] ** 2 - r[:-1] ** 2)
    cavity = np.pi * r[0] ** 2
    np.testing.assert_allclose(got[:-1], areas, rtol=_tolerance(bp), atol=0)
    np.testing.assert_allclose(got[-1], cavity, rtol=_tolerance(bp), atol=_tolerance(bp))
    np.testing.assert_allclose(areas, part.chart.measure(r), rtol=4 * _EPS, atol=0)


@pytest.mark.verifies("geometry-measure-on-lines")
@pytest.mark.rests_on(_CHORD + "test_every_slot_length_matches_the_closed_form",
                      _CHART + "test_the_beam_density_of_each_chart")
@pytest.mark.parametrize("bp", [_SOLID, _HOLLOW], ids=["solid", "hollow"])
def test_the_beam_integral_of_a_spherical_shell_is_its_volume(bp) -> None:
    """M2: ``int_0^R l_j(b) 2 pi b db`` is the volume of shell ``j``; the cavity gives ``(4/3) pi r_0^3``.

    First red: the sphere's density without its factor ``b`` (the planar
    measure; refuted 2026-10-05 for the sphere, Cauchy's 4/3 against pi/2).
    """
    part = ConcentricPartition(Chart(CoordSystem.SPHERICAL), bp)
    om = np.array([0.0, 1.0, 0.0])
    b, w = _impact_rule(bp)
    got = np.sum((w * part.chart.beam_density(b, om))[:, None] * _region_lengths(part, b, om), axis=0)
    r = np.asarray(bp)
    volumes = (4.0 / 3.0) * np.pi * (r[1:] ** 3 - r[:-1] ** 3)
    np.testing.assert_allclose(got[:-1], volumes, rtol=_tolerance(bp), atol=0)
    np.testing.assert_allclose(got[-1], (4.0 / 3.0) * np.pi * r[0] ** 3, rtol=_tolerance(bp), atol=_tolerance(bp))
    np.testing.assert_allclose(volumes, part.chart.measure(r), rtol=4 * _EPS, atol=0)


@pytest.mark.verifies("geometry-measure-on-lines")
@pytest.mark.rests_on(_CHORD + "test_the_slab_chord_in_both_orientations")
def test_the_beam_integral_of_a_slab_region_is_its_width() -> None:
    """M1, slab: per unit area, ``l_j |Omega_x|`` is the width ``r_{j+1} - r_j`` for every non-parallel direction.

    First red: the slab's density taken as 1 (each row reads the width over ``|Omega_x|``).
    """
    part = ConcentricPartition(Chart(CoordSystem.CARTESIAN), (-0.7, 0.3, 1.1, 2.0))
    r = np.asarray(part.breakpoints)
    for omega in ((0.6, 0.8, 0.0), (-0.28, 0.0, 0.96), (1.0, 0.0, 0.0)):
        om = np.array(omega)
        ch = part.chord(Line.through(np.zeros(3), om))
        lengths = np.array([np.sum(ch.slot_length[ch.slot_region == j]) for j in range(3)])
        np.testing.assert_allclose(lengths * part.chart.beam_density(np.array(0.0), om), np.diff(r), rtol=4 * _EPS, atol=0)


# [M] 2026-10-05: the direction-averaged rows err by 3.3e-16 (cylinder) and 1.6e-15
# (slab) relative with the 32-point rules; 10x the larger, rounded up.
_CAUCHY_RTOL = 2e-14


def _direction_rule(upper: float = np.pi, n: int = 32) -> tuple[np.ndarray, np.ndarray]:
    """Polar angles and weights for ``int dOmega = 2 pi int_0^upper sin(theta) dtheta`` (the azimuth is a symmetry of the body).

    The slab's density ``|cos(theta)|`` has a kink at ``pi/2``, so its rule is
    the hemisphere ``[0, pi/2]`` (the other is its mirror image); the cylinder's
    ``sin(theta)`` is smooth on ``[0, pi]``.
    """
    x, w = leggauss(n)
    theta = 0.5 * upper * (x + 1.0)
    return theta, 2.0 * np.pi * 0.5 * upper * w * np.sin(theta)


@pytest.mark.verifies("geometry-cauchy-mean-chord")
@pytest.mark.rests_on("tests/gates/geometry/test_line_measure.py::test_the_beam_integral_of_a_cylinder_region_is_its_area",
                      "tests/gates/geometry/test_line_measure.py::test_the_beam_integral_of_a_spherical_shell_is_its_volume")
def test_cauchys_mean_chord_is_four_volume_over_surface() -> None:
    """M4: the mean chord over isotropic uniform lines meeting a convex body is ``4V/S``.

    Sphere of radius 2: ``4R/3``. Infinite cylinder of radius 2, isotropic 3-D
    lines: ``2R`` per unit height. Slab of width 2.7 (per unit area, both faces):
    ``2 L``.
    Scope (X4, said here): for a body of revolution the sphere row is M2 divided
    by the projected area, so it is not independent of M2; it adds only the
    measure of the lines that hit. The cylinder and slab rows add the direction
    average (Cauchy's mean projected area ``S/4``), which is independent content:
    with the cylinder's density taken as 1 the cylinder row reads ``[M]`` 2.467 R
    (the cross-domain review's measurement), not 2 R.
    The direction rule (32-point Gauss-Legendre in ``theta``, smooth integrands)
    is measured (``_CAUCHY_RTOL``).
    """
    R = 2.0
    b, w = _impact_rule((0.0, R))
    sphere = ConcentricPartition(Chart(CoordSystem.SPHERICAL), (0.0, R))
    om = np.array([0.0, 1.0, 0.0])
    dens = sphere.chart.beam_density(b, om)
    total = np.sum(w * dens * _region_lengths(sphere, b, om)[:, 0])
    hit = np.sum(w * dens)
    np.testing.assert_allclose(total / hit, 4.0 * R / 3.0, rtol=_tolerance((0.0, R)), atol=0)

    cylinder = ConcentricPartition(Chart(CoordSystem.CYLINDRICAL), (0.0, R))
    theta, wt = _direction_rule()
    num = den = 0.0
    for th, wth in zip(theta, wt):
        o = np.array([0.0, np.sin(th), np.cos(th)])
        d = cylinder.chart.beam_density(b, o)
        num += wth * np.sum(w * d * _region_lengths(cylinder, b, o)[:, 0])
        den += wth * np.sum(w * d)
    np.testing.assert_allclose(num / den, 2.0 * R, rtol=_CAUCHY_RTOL, atol=0)

    slab = ConcentricPartition(Chart(CoordSystem.CARTESIAN), (-0.7, 2.0))
    num = den = 0.0
    for th, wth in zip(*_direction_rule(np.pi / 2)):
        o = np.array([np.cos(th), np.sin(th), 0.0])
        ch = slab.chord(Line.through(np.zeros(3), o))
        d = float(slab.chart.beam_density(np.array(0.0), o))
        num += wth * d * float(np.sum(ch.slot_length))
        den += wth * d
    np.testing.assert_allclose(num / den, 2.0 * 2.7, rtol=_CAUCHY_RTOL, atol=0)
