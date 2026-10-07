"""Gates for :meth:`orpheus.geometry.chart.Chart.line_domain` and its :class:`LineDomain`.

P1 step (b), third rung, of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "P1 step (b), third
rung: API sketch", item 2, and the user's ruling of 2026-10-06 that the
cylinder's second coordinate is the polar angle theta = arccos mu_z).
Verification spec ``scratch/characteristic_architecture/p1_step_b3/spec.md``,
rows LD1-LD5.

The object: the oriented lines of space modulo the chart's group, a box in at
most two coordinates, with the invariant measure dA_perp dOmega counted per
unit of the discarded measure. Its density is the beam density
(``Chart.beam_density``, gated in ``test_chart.py``) times the direction
measure the quotient folds into one point: 4 pi on the sphere, 2 pi x 2 sin
theta on the cylinder (the azimuth, the fold theta <-> pi - theta, and
dOmega = sin theta dtheta dphi), 2 pi on the slab. In closed form 8 pi^2 b,
8 pi sin^2 theta and 2 pi |mu|, typed here.

Independence (X4). LD3 (Cauchy) integrates the domain's density times the 3-D
chord length the P0 kernel returns (``ConcentricPartition.chord``, gated against
closed forms in ``test_chord.py``) with a rule written here, and compares with
4 pi times ``Chart.measure``, a derivation sharing no identity with a chord
integral. Stabiliser, declared: on the cylinder the chord carries 1/|P Omega|
(the kernel's obliquity) and the density |P Omega|, so their product is blind to
a common error in both; LD2 pins the density pointwise.

Levels: LD2 ``l0`` with ``verifies("geometry-measure-on-lines")``; LD3 ``l1``
with ``verifies("geometry-line-domain")`` (minted 2026-10-07); the others are
``foundation``.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from orpheus.geometry.chart import Chart, LineShape
from orpheus.geometry.chord import ConcentricPartition
from orpheus.geometry.coord import CoordSystem

_EPS = float(np.finfo(float).eps)
_HERE = "tests/gates/geometry/test_line_domain.py::"
_CHART = "tests/gates/geometry/test_chart.py::"
_CHORD = "tests/gates/geometry/test_chord.py::"
_MEASURE = "tests/gates/geometry/test_line_measure.py::"

_SLAB, _CYL, _SPH = (Chart(c) for c in (CoordSystem.CARTESIAN, CoordSystem.CYLINDRICAL, CoordSystem.SPHERICAL))

#: (id, chart, shape, axes, bounds): the table by hand, from the groups. O(3) folds a line's orientation; D_inf_h
#: (with the axial translations) folds the azimuth and theta <-> pi - theta; O(2)_x fixes the kept space, so a
#: slab line and its reverse are two orbits (mu and -mu).
_TABLE = [
    ("sphere", _SPH, LineShape.IMPACT, ("impact",), ((0.0, math.inf),)),
    ("cylinder", _CYL, LineShape.IMPACT_POLAR, ("impact", "polar_angle"), ((0.0, math.inf), (0.0, math.pi / 2))),
    ("slab", _SLAB, LineShape.COSINE, ("cosine",), ((-1.0, 1.0),)),
]


@pytest.mark.foundation
@pytest.mark.parametrize(("chart", "shape", "axes", "bounds"), [t[1:] for t in _TABLE], ids=[t[0] for t in _TABLE])
def test_the_line_domain_of_each_chart(chart, shape, axes, bounds) -> None:
    """[LD1] Shape, axes and bounds per chart, by hand; every ``LineShape`` member is produced by exactly one chart.

    First reds: the cylinder's polar bound pi (no fold of theta and pi - theta;
    LD3 reds by a factor 2 too); the second axis left as the axial cosine.
    """
    domain = chart.line_domain()
    assert domain.shape is shape
    assert domain.axes == axes
    assert domain.bounds == bounds
    produced = [c.line_domain().shape for _, c, *_ in _TABLE]
    assert sorted(s.name for s in produced) == sorted(s.name for s in LineShape)


def _grid(chart: Chart) -> np.ndarray:
    match chart.line_domain().shape:
        case LineShape.IMPACT:
            return np.array([[0.0], [0.3], [1.9], [37.0]])
        case LineShape.IMPACT_POLAR:
            b, t = np.meshgrid([0.0, 0.3, 1.9, 37.0], [0.0, 0.2, 0.7, 1.3, math.pi / 2], indexing="ij")
            return np.stack([b.ravel(), t.ravel()], axis=-1)
        case LineShape.COSINE:
            return np.array([[-1.0], [-0.3], [0.0], [0.7], [1.0]])
    raise AssertionError(chart)


def _by_hand(chart: Chart, q: np.ndarray) -> np.ndarray:
    match chart.line_domain().shape:
        case LineShape.IMPACT:
            return 8.0 * math.pi ** 2 * q[:, 0]
        case LineShape.IMPACT_POLAR:
            return 8.0 * math.pi * np.sin(q[:, 1]) ** 2
        case LineShape.COSINE:
            return 2.0 * math.pi * np.abs(q[:, 0])
    raise AssertionError(chart)


@pytest.mark.l0
@pytest.mark.verifies("geometry-measure-on-lines")
@pytest.mark.parametrize("chart", [t[1] for t in _TABLE], ids=[t[0] for t in _TABLE])
@pytest.mark.rests_on(_HERE + "test_the_line_domain_of_each_chart", _CHART + "test_the_beam_density_of_each_chart")
def test_the_density_is_the_beam_density_times_the_folded_directions(chart) -> None:
    """[LD2] density = 8 pi^2 b (sphere), 8 pi sin^2 theta (cylinder), 2 pi |mu| (slab); and = beam density x fold.

    The closed forms typed here, 8 ulp; the factorisation leg reads
    ``Chart.beam_density`` of the domain's own lines times the fold (4 pi, 4 pi
    sin theta, 2 pi), 4 ulp. Not divided by 4 pi (that is the consumer's). First
    reds: the sphere's fold 2 pi; the cylinder without one sin theta (the
    Jacobian of mu_z -> theta, or the beam's |P Omega|); the slab's density
    without |mu|; a division by 4 pi here.
    """
    domain = chart.line_domain()
    q = _grid(chart)
    got = domain.density(q)
    want = _by_hand(chart, q)
    np.testing.assert_allclose(got, want, rtol=8 * _EPS, atol=8 * _EPS * float(np.max(want)))
    impact = q[:, 0] if domain.shape is not LineShape.COSINE else np.zeros(len(q))
    beam = chart.beam_density(impact, domain.lines(q).direction)
    fold = {LineShape.IMPACT: 4.0 * math.pi + 0.0 * q[:, 0], LineShape.COSINE: 2.0 * math.pi + 0.0 * q[:, 0],
            LineShape.IMPACT_POLAR: 4.0 * math.pi * np.sin(q[:, -1])}[domain.shape]
    np.testing.assert_allclose(got, beam * fold, rtol=4 * _EPS, atol=4 * _EPS * float(np.max(want)))


_N = 32


def _pieces_in_b(breakpoints) -> tuple[np.ndarray, np.ndarray]:
    """Nodes and weights on b in [0, R]: per piece [lo, hi], Gauss-Legendre in phi with b = hi sin(phi) (``test_line_measure``)."""
    x, w = np.polynomial.legendre.leggauss(_N)
    edges = [0.0, *[r for r in breakpoints if r > 0.0]]
    nodes, weights = [], []
    for lo, hi in zip(edges, edges[1:]):
        a, b = math.asin(lo / hi), math.pi / 2
        phi = 0.5 * (b - a) * x + 0.5 * (b + a)
        nodes.append(hi * np.sin(phi))
        weights.append(0.5 * (b - a) * w * hi * np.cos(phi))
    return np.concatenate(nodes), np.concatenate(weights)


def _test_rule(chart: Chart, breakpoints) -> tuple[np.ndarray, np.ndarray]:
    """The test's own rule on the domain (never the reference's ``LineRule``)."""
    x, w = np.polynomial.legendre.leggauss(_N)
    match chart.line_domain().shape:
        case LineShape.IMPACT:
            b, wb = _pieces_in_b(breakpoints)
            return b[:, None], wb
        case LineShape.IMPACT_POLAR:
            b, wb = _pieces_in_b(breakpoints)
            theta, wt = math.pi / 4 * (x + 1.0), math.pi / 4 * w
            bb, tt = np.meshgrid(b, theta, indexing="ij")
            return np.stack([bb.ravel(), tt.ravel()], axis=-1), np.outer(wb, wt).ravel()
        case LineShape.COSINE:
            half = 0.5 * (x + 1.0)
            return np.concatenate([-half, half])[:, None], np.concatenate([0.5 * w, 0.5 * w])
    raise AssertionError(chart)


_BODIES = [
    ("sphere_solid", _SPH, (0.0, 0.5, 1.5, 2.0)),
    ("sphere_hollow", _SPH, (0.4, 0.5, 1.5, 2.0)),
    ("cylinder_solid", _CYL, (0.0, 0.5, 1.5, 2.0)),
    ("cylinder_hollow", _CYL, (0.4, 0.5, 1.5, 2.0)),
    ("slab", _SLAB, (0.0, 0.4, 1.5, 2.3)),
]


@pytest.mark.l1
@pytest.mark.verifies("geometry-line-domain")
@pytest.mark.parametrize(("chart", "breakpoints"), [b[1:] for b in _BODIES], ids=[b[0] for b in _BODIES])
@pytest.mark.rests_on(_HERE + "test_the_density_is_the_beam_density_times_the_folded_directions",
                      _CHORD + "test_every_slot_length_matches_the_closed_form",
                      _MEASURE + "test_cauchys_mean_chord_is_four_volume_over_surface")
def test_cauchys_formula_on_the_line_domain(chart, breakpoints) -> None:
    """[LD3] The integral over the domain of density x the 3-D chord length through the body is 4 pi x its measure.

    The chord lengths are the kernel's (every slot inside the body: the cavity of a
    hollow body is outside it); the measure is ``Chart.measure`` (the one
    definition); the rule is this file's: per piece b = r_j sin(phi) (absorbing
    each chord's square-root end), Gauss-Legendre in theta, and in mu on each
    sign. 64 n P eps (the reduction depth; `[M]` on the spec's prototype with the
    reference's own rule: 0 sphere, 7.1e-13 cylinder at 16, 2.2e-16 slab). First
    reds (`[M]` the spec's arms): the sphere's fold 2 pi: -0.50; the cylinder's fold 2 pi: -0.50;
    the cylinder density without |P Omega|: +0.53; the slab without |mu|: +13.
    """
    domain = chart.line_domain()
    q, w = _test_rule(chart, breakpoints)
    partition = ConcentricPartition(chart, breakpoints)
    chord = partition.chord(domain.lines(q))
    inside = (chord.slot_region < partition.n_regions) & np.isfinite(chord.slot_length)
    length = np.sum(np.where(inside, chord.slot_length, 0.0), axis=-1)
    total = float(np.sum(w * domain.density(q) * length))
    want = 4.0 * math.pi * float(np.sum(chart.measure(np.asarray(breakpoints))))
    pieces = len([r for r in breakpoints if r > 0.0])
    assert abs(total / want - 1.0) <= 64 * _N * pieces * _EPS, f"{total / want - 1.0:.2e}"


@pytest.mark.foundation
@pytest.mark.parametrize("chart", [t[1] for t in _TABLE], ids=[t[0] for t in _TABLE])
@pytest.mark.rests_on(_HERE + "test_the_line_domain_of_each_chart",
                      "tests/gates/geometry/test_chart_directions.py::test_the_impact_parameter_is_the_kernels_and_the_closed_form")
def test_the_representative_lines_carry_their_coordinates(chart) -> None:
    """[LD4] b of each representative line under ``Chart.image`` is the coordinate b to 2 ulp; cos(theta) its axial
    cosine, bitwise; the slab line's Omega_x is mu, bitwise; every direction a unit vector.

    Not bitwise in b (`[M]` 2026-10-06: 1 ulp at b = 0.3, theta = 0.2): a ``Line``
    stores its moment p x Omega and returns its foot as Omega x m, so the foot's
    b is 0.3 (Omega_x^2 + Omega_z^2) rounded, and the docstring of
    ``LineDomain.lines`` ("bit for bit") claims more than the kernel's Plücker
    representation delivers. First red: the representative foot at (b, 0, 0) with an in-plane direction
    having an x component (the image's b is then not the coordinate).
    """
    domain = chart.line_domain()
    q = _grid(chart)
    lines = domain.lines(q)
    np.testing.assert_allclose(np.linalg.norm(lines.direction, axis=-1), 1.0, rtol=0.0, atol=2 * _EPS)
    match domain.shape:
        case LineShape.IMPACT | LineShape.IMPACT_POLAR:
            np.testing.assert_allclose(chart.image(lines).impact_parameter, q[:, 0], rtol=2 * _EPS, atol=0.0)
            if domain.shape is LineShape.IMPACT_POLAR:
                np.testing.assert_array_equal(lines.direction[:, 2], np.cos(q[:, 1]))
        case LineShape.COSINE:
            np.testing.assert_array_equal(lines.direction[:, 0], q[:, 0])


_BAD = [
    ("sphere_b_negative", _SPH, [[-1e-12]], "lies in"),
    ("sphere_b_nan", _SPH, [[math.nan]], "lies in"),
    ("sphere_b_inf", _SPH, [[math.inf]], "lies in"),
    ("cylinder_theta_above", _CYL, [[0.3, math.pi / 2 + 1e-12]], "lies in"),
    ("cylinder_theta_below", _CYL, [[0.3, -1e-12]], "lies in"),
    ("slab_mu_above", _SLAB, [[1.0 + 1e-12]], "lies in"),
    ("cylinder_one_coordinate", _CYL, [[0.3]], "have shape"),
    ("sphere_two_coordinates", _SPH, [[0.3, 0.1]], "have shape"),
]


@pytest.mark.foundation
@pytest.mark.parametrize("verb", ["lines", "density"])
@pytest.mark.parametrize(("chart", "bad", "fragment"), [b[1:] for b in _BAD], ids=[b[0] for b in _BAD])
def test_a_coordinate_off_the_domain_is_refused(chart, bad, fragment, verb) -> None:
    """[LD5] A coordinate outside its bound or not finite, or of the wrong width, is refused by its own fragment.

    The two fragments are disjoint (asserted once). First red: either guard
    deleted (its legs return a line or a density silently).
    """
    assert "lies in" not in "have shape" and "have shape" not in "lies in"
    with pytest.raises(ValueError, match=fragment):
        getattr(chart.line_domain(), verb)(np.array(bad))
