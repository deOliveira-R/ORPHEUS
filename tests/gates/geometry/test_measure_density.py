"""Gates for the density of the one measure: :meth:`MeasureCoordinate.derivative`, :meth:`CoordSystem.measure_density`, :meth:`Chart.measure_density`.

P1 step (b), third rung, of the characteristic-reference campaign
(``.claude/plans/characteristic_reference_architecture.md``, "P1 step (b), third
rung: API sketch", item 1). Verification spec
``scratch/characteristic_architecture/p1_step_b3/spec.md``, rows MD1 and MD2
(MD3, the route through the basis and the wall coupling, is in
``tests/gates/derivations/test_characteristic_assembly.py``: it reads
``derivations``, which this layer does not import).

The measure has one definition, ``c (T(b) - T(a))`` with ``T(r) = r^p``
(``tests/gates/mesh/test_partition.py::TestTheOneMeasure``). Its density
``c T'(r)`` is that definition's derivative, and by the coarea formula the area
of the level set at ``r``: ``4 pi r^2`` (sphere), ``2 pi r`` per unit height
(cylinder), 1 per unit transverse area (slab).

Independence (X4). MD1's expectations are typed per chart by hand, never from
``measure_constant``. MD2 is a fundamental-theorem identity: two-point
Gauss-Legendre of the density (exact for its degree <= 2) against the measure
in mpmath at 40 digits; it shares the definition, not the derivative. Its
tolerance carries the conditioning ``kappa = max(|a|, |b|) / (b - a)``: the
measure ``c (b^p - a^p)`` cancels on a narrow interval far from 0 (`[M]`
2026-10-06, the spec's P-5: 2568 ulp raw on the cylinder over 2000 random
intervals, 0.9 ulp over ``1 + kappa``).

Levels: ``foundation`` until the archivist mints ``geometry-measure-density``
(the spec's MD1 and MD2 are then ``l0``).
"""
from __future__ import annotations

import math

import mpmath as mp
import numpy as np
import pytest

from orpheus.geometry.chart import Chart
from orpheus.geometry.coord import CoordSystem, MeasureCoordinate

_EPS = float(np.finfo(float).eps)
_ONE_MEASURE = "tests/gates/mesh/test_partition.py::TestTheOneMeasure::test_the_measure_is_c_times_the_difference_of_T"
_HERE = "tests/gates/geometry/test_measure_density.py::"

#: (coordinate system, the density by hand, the measure's constant by hand, the exponent p)
_CHARTS = [
    (CoordSystem.CARTESIAN, lambda r: 1.0 + 0.0 * r, 1, 1),
    (CoordSystem.CYLINDRICAL, lambda r: 2.0 * math.pi * r, mp.pi, 2),
    (CoordSystem.SPHERICAL, lambda r: 4.0 * math.pi * r * r, 4 * mp.pi / 3, 3),
]
_IDS = ["slab", "cylinder", "sphere"]
_R = np.array([0.0, 0.37, 1.1, 1.9, 2.6e3])


@pytest.mark.foundation
@pytest.mark.parametrize(("coord", "by_hand", "_c", "p"), _CHARTS, ids=_IDS)
@pytest.mark.rests_on(_ONE_MEASURE)
def test_the_measure_density_is_each_charts_level_set_area(coord, by_hand, _c, p) -> None:
    """[MD1] T'(r) = p r^(p-1); c T'(r) = 1, 2 pi r, 4 pi r^2; the chart's verb is the system's, bitwise.

    The densities are typed here per chart (not from ``measure_constant``). 4 ulp
    relative (one power, two products). First reds: ``derivative`` returning
    p r^p (one power too many); the factor p dropped; ``Chart.measure_density``
    re-derived instead of delegating (its leg is ``array_equal``).
    """
    np.testing.assert_allclose(MeasureCoordinate(p).derivative(_R), p * _R ** (p - 1), rtol=2 * _EPS, atol=0.0)
    got = coord.measure_density(_R)
    np.testing.assert_allclose(got, by_hand(_R), rtol=4 * _EPS, atol=0.0)
    np.testing.assert_array_equal(Chart(coord).measure_density(_R), got)
    assert got.shape == _R.shape


@pytest.mark.foundation
@pytest.mark.parametrize(("coord", "_by_hand", "c", "p"), _CHARTS, ids=_IDS)
@pytest.mark.rests_on(_HERE + "test_the_measure_density_is_each_charts_level_set_area", _ONE_MEASURE)
def test_the_density_integrates_to_the_one_measure(coord, _by_hand, c, p) -> None:
    """[MD2] Two-point Gauss-Legendre of ``measure_density`` on [a, b] equals c (b^p - a^p), 2000 seeded intervals.

    Against mpmath at 40 digits: 4 (1 + kappa) ulp (`[M]` 2026-10-06 on the
    spec's prototype, 1.1 (1 + kappa)); against ``Chart.measure([a, b])``, the
    production definition: 8 (1 + kappa) ulp. The two-point rule is exact for a
    polynomial of degree <= 3, so the gap is rounding only. First reds: the
    density's exponent p instead of p - 1 (O(1)); the constant c dropped.
    """
    rng = np.random.default_rng(20261006 + p)
    ends = np.sort(rng.uniform(0.0, 3.0, size=(2000, 2)), axis=-1)
    x, w = np.polynomial.legendre.leggauss(2)
    a, b = ends[:, :1], ends[:, 1:]
    nodes = 0.5 * (b - a) * x + 0.5 * (b + a)
    got = np.sum(0.5 * (b - a) * w * coord.measure_density(nodes), axis=-1)
    kappa = np.max(np.abs(ends), axis=-1) / (ends[:, 1] - ends[:, 0])
    with mp.workdps(40):
        exact = np.array([float(c * (mp.mpf(hi) ** p - mp.mpf(lo) ** p)) for lo, hi in ends])
    worst = np.max(np.abs(got - exact) / np.abs(exact) / (1.0 + kappa))
    assert worst <= 4 * _EPS, f"{worst / _EPS:.1f} ulp over (1 + kappa)"
    measure = np.array([Chart(coord).measure(np.array([lo, hi]))[0] for lo, hi in ends])
    assert np.max(np.abs(got - measure) / np.abs(measure) / (1.0 + kappa)) <= 8 * _EPS
