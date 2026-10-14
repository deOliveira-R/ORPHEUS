"""Gates for the Gauss rule of the attenuation weight (:mod:`~orpheus.derivations.continuous.characteristic.attenuation_rule`).

P3 step 3f of the characteristic-reference campaign, route B (the user's ruling of 2026-10-13): along a piece of a line
every integral is a smooth factor against the attenuation to one of the piece's ends, integrated by the m-point Gauss
rule of that weight, its nodes held as offsets from the attenuated end. Spec:
``scratch/characteristic_architecture/p3/product/gate_spec.md`` (REVISION 2: RB1-RB3, PT1).

The ladder: the reference's own anchor (RB0: the Bessel moment equals the by-parts closed form); the rule's
construction laws (RB1 positivity, RB2 constants); its exactness (RB3, l0); the triangle built on it (PT1, l0, the
orientation gate at the weight level). References: ``_characteristic_product_mp.py`` (mpmath, no import of the code
under test). First reds: in-process arms of ``scratch/characteristic_architecture/p3/product/ta/carve/rule_arms.py``
(``RULE_ARM=<arm> pytest -p rule_arms``), named per row.
"""
from __future__ import annotations

import mpmath as mp
import numpy as np
import pytest

from orpheus.derivations.continuous.characteristic.attenuation_rule import (
    BAND,
    MOST_NODES,
    PANEL_WIDTH,
    TRUNCATION,
    AttenuationRule,
    zeroth_moment,
)

from . import _characteristic_product_mp as R

pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning")

_HERE = "tests/gates/derivations/test_characteristic_attenuation_rule.py::"
_EPS = float(np.finfo(float).eps)
_MS = (1, 2, 3, 4, 6, 8, 12, 16, 17, 20, 24, 28, 32, 36, 40)
_XS = (0.0, 5e-324, 1e-300, 1e-8, 1e-3, 0.1, 1.0, 3.0, 10.0, 16.0, 30.0, 100.0, 1e3, 1e5)
#: RB3's band, in eps of the zeroth moment: ``attenuation_rule.BAND`` (64 eps = 1.4e-14), spelled once in production
#: and also the law's target (``grading.ALONG_LINE_ERROR = BAND``, the user's ruling of 2026-10-13: the node-count law
#: spends no more than the rule's own band). `[M]` 2026-10-13 on the round-2 build over this file's grid
#: (``scratch/.../p3/product/ta/carve/probe_rb_worst.log``): worst 60.1 eps (m = 28, x = 16, k = 34); RB2 worst 1.53 ulp.
_RB3_BAND = BAND / float(np.finfo(float).eps)


def _grid(m: int) -> np.ndarray:
    """``_XS``, m/2, m, 2m, and the construction's two boundaries approached to 1e-12 from both sides: 2x = PANEL_WIDTH
    (one panel in s below) and 2x = TRUNCATION (the measure truncated above)."""
    edges = []
    for two_x in (PANEL_WIDTH, TRUNCATION):
        x = two_x / 2.0
        edges += [x * (1 - 1e-12), x, x * (1 + 1e-12)]
    return np.array(sorted(set(_XS + tuple(edges) + (m / 2.0, float(m), 2.0 * m))))


def _zeroth(x: float) -> float:
    with mp.workdps(40):
        return 2.0 if x == 0.0 else float(-mp.expm1(-2 * mp.mpf(repr(float(x)))) / mp.mpf(repr(float(x))))


_RB0 = [(k, x) for k in (0, 1, 5, 11, 23, 31) for x in (0.0, 1e-6, 0.3, 4.0, 60.0, 900.0)]


@pytest.mark.foundation
@pytest.mark.parametrize(("k", "x"), _RB0, ids=[f"k{k}-x{x:g}" for k, x in _RB0])
def test_the_bessel_moment_is_the_by_parts_closed_form(k: int, x: float) -> None:
    """[RB0; foundation, the reference's anchor] 2 ĩ_k(x) from mpmath's Bessel function equals the by-parts closed form
    of ∫P_k e^{-x(1-s)} within (8 + k) eps: RB3 reads the first, which is fast at k ~ 80. A test of the reference."""
    want, got = R.legendre_toward_end(k, x), R.legendre_moment(k, x)
    assert abs(got - want) <= (8 + k) * _EPS * max(abs(float(want)), 1e-290), (got, want)


@pytest.mark.foundation
@pytest.mark.parametrize("m", _MS, ids=[f"m{m}" for m in _MS])
def test_the_rule_has_positive_weights_and_distinct_offsets_inside(m: int) -> None:
    """[RB1; foundation, THEOREM] A Gauss rule of a positive weight: m distinct offsets strictly inside (0, 2),
    ASCENDING from the attenuated end at every x (since round 2 every rule is built in one scaled variable; the
    small-x branch used to return them descending), and m positive, finite weights, count m, at every x of the grid.
    First reds: ``extra0`` (the measure discretised by m points per panel); ``descending`` (the nodes reversed)."""
    xs = _grid(m)
    rule = AttenuationRule.of(xs, m)
    assert rule.offset.shape == rule.weight.shape == (xs.size, m)
    np.testing.assert_array_equal(rule.count, m)
    for x, d, w in zip(xs, rule.offset, rule.weight):
        assert np.all(np.isfinite(d)) and np.all(np.isfinite(w)), (x, d, w)
        assert np.all(w > 0.0), (x, w)
        gaps = np.diff(np.asarray(d, dtype=float))
        assert np.all(d > 0.0) and np.all(d < 2.0) and np.all(gaps > 0.0), (x, d)          # ascending, distinct


@pytest.mark.l0
@pytest.mark.verifies("characteristic-traversal-integrals")
@pytest.mark.parametrize("m", _MS, ids=[f"m{m}" for m in _MS])
@pytest.mark.rests_on(_HERE + "test_the_rule_has_positive_weights_and_distinct_offsets_inside")
def test_the_rule_integrates_a_constant_to_two_ulp(m: int) -> None:
    """[RB2; l0, THEOREM] sum_q w_q = (1 - e^{-2x})/x (2 at x = 0) within 2 ulp, and ``zeroth_moment`` is that closed
    form within 1 ulp: the weights are rescaled to it, so the one-level rule integrates the attenuation itself to the
    ulp (the invariant G1 asserted of the retired exponential grading; its marker moved here, retirement-audit C.13).
    T7's 4-ulp leg rests on it. First red: ``no_rescale`` (the rescale line removed)."""
    xs = _grid(m)
    rule = AttenuationRule.of(xs, m)
    for x, w, z in zip(xs, rule.weight, zeroth_moment(xs)):
        want = _zeroth(float(x))
        assert abs(z - want) <= _EPS * want, (x, z, want)
        assert abs(w.sum() - want) <= 2 * _EPS * want, (x, (w.sum() - want) / (_EPS * want))


@pytest.mark.l0
@pytest.mark.parametrize("m", _MS, ids=[f"m{m}" for m in _MS])
@pytest.mark.rests_on(_HERE + "test_the_rule_has_positive_weights_and_distinct_offsets_inside",
                      _HERE + "test_the_bessel_moment_is_the_by_parts_closed_form")
def test_the_rule_is_exact_on_every_legendre_polynomial_below_degree_two_m(m: int) -> None:
    """[RB3; l0, REFERENCE] sum_q w_q P_k(1 - delta_q) = ∫P_k(s) e^{-x(1-s)} ds for k <= 2m - 1, within 64 eps of the
    zeroth moment (|P_k| <= 1, so it bounds ∫omega|P_k|), the sum formed in mpmath from the double offsets and weights
    (``legendre_sums``: the row measures the rule, not its evaluation; `[M]` a double P_71 alone cost 74.8 eps at
    x = 1e5). Reference: 2 ĩ_k through mpmath's Bessel function (RB0). The small-x branch (2x <= PANEL_WIDTH) is
    built in the offset too since the carve's review: the grid's two sides of 2x = 4 decide it. First reds:
    ``in_s`` (the offsets formed through s, delta = 1 - (1 - u/x)); ``panel_wide`` (one panel up to 2x = 40);
    ``extra0``. (The Stieltjes re-orthogonalisation had no witness at m <= 80 and was removed, the user's ruling of
    2026-10-13.)"""
    xs = _grid(m)
    rule = AttenuationRule.of(xs, m)
    for x, d, w in zip(xs, rule.offset, rule.weight):
        nu0 = _zeroth(float(x))
        sums = R.legendre_sums(d, w, 2 * m)
        for k in range(2 * m):
            want = R.legendre_moment(k, float(x))
            assert abs(sums[k] - want) <= _RB3_BAND * _EPS * nu0, (x, k, (sums[k] - want) / (_EPS * nu0))


@pytest.mark.foundation
def test_rules_of_several_counts_in_one_batch_are_each_rule_alone() -> None:
    """[RB4; foundation, THEOREM] ``AttenuationRule.of(x, counts)`` with mixed counts holds, per entry, the rule of its
    count built alone, ``array_equal`` on its live nodes, zero weights on the padding and the last offset repeated;
    ``count`` is the requested one. First red: the padding carrying the next rule's weights."""
    xs = np.array([0.0, 0.7, 3.0, 25.0, 900.0, 0.7])
    counts = np.array([3, 7, 1, 12, 5, 12])
    batch = AttenuationRule.of(xs, counts)
    np.testing.assert_array_equal(batch.count, counts)
    for i, (x, m) in enumerate(zip(xs, counts)):
        alone = AttenuationRule.of(np.array([x]), int(m))
        np.testing.assert_array_equal(batch.offset[i, :m], alone.offset[0])
        np.testing.assert_array_equal(batch.weight[i, :m], alone.weight[0])
        np.testing.assert_array_equal(batch.weight[i, m:], 0.0)
        np.testing.assert_array_equal(batch.offset[i, m:], alone.offset[0, -1])


@pytest.mark.foundation
def test_a_rule_at_a_negative_width_or_no_nodes_is_refused() -> None:
    """[RB5; foundation] A negative or non-finite optical half-width, and a count below 1, are refused, each naming its
    reason; a valid call is not (the positive leg). First red: either check removed."""
    with pytest.raises(ValueError, match="finite and non-negative"):
        AttenuationRule.of(np.array([-1e-3]), 3)
    with pytest.raises(ValueError, match="finite and non-negative"):
        AttenuationRule.of(np.array([np.inf]), 3)
    with pytest.raises(ValueError, match="at least one node"):
        AttenuationRule.of(np.array([1.0]), 0)
    AttenuationRule.of(np.array([1.0]), 1)


@pytest.mark.foundation
def test_a_rule_past_the_measured_node_count_is_refused() -> None:
    """[RB6; foundation, SCOPE-BOUNDARY] ``AttenuationRule.of`` refuses a count above ``MOST_NODES`` (40), naming the
    boundary: RB3 measured the rule within its band up to 40 nodes only (`[M]` 2026-10-13: 69.8 eps at m = 56, past
    the 64-eps band, ``carve/probe_reorth.log``). At exactly 40 it builds (the positive leg; RB1-RB3 gate m = 40).
    First red: the boundary removed."""
    assert MOST_NODES == 40
    with pytest.raises(ValueError, match="a scope boundary: MOST_NODES"):
        AttenuationRule.of(np.array([1.0, 3.0]), np.array([12, MOST_NODES + 1]))
    AttenuationRule.of(np.array([1.0]), MOST_NODES)


#: Two distinct cubics: f . V g != g . V f, so the triangle's orientation is visible.
_F, _G = [0.4, -1.1, 0.7, 0.9], [1.3, 0.5, -0.8, 0.2]
_TRI_XS = (0.0, 1e-6, 0.3, 4.0, 30.0, 100.0, 1e3)


@pytest.mark.l0
@pytest.mark.parametrize("m", (2, 4, 7, 12), ids=["m2", "m4", "m7", "m12"])
@pytest.mark.rests_on(_HERE + "test_the_rule_is_exact_on_every_legendre_polynomial_below_degree_two_m")
def test_the_triangle_rule_is_the_closed_form_double_integral_in_its_orientation(m: int) -> None:
    """[PT1; l0, REFERENCE: the orientation gate at the weight level] ``AttenuationRule.of(x, m).triangle()``:
    sum_n w_n f(1 - lead_n) g(1 - lag_n) equals ∫f(s)∫_{-1}^{s} g(s') e^{-x(s-s')} ds' ds (by parts in mpmath) within
    32 eps of max|f| max|g| times the triangle's mass; cubics from m = 4 (G of degree 7 outer, degree 6 inner), lines
    at m = 2; the swapped pair differs by more than 1e-6. The lag point is upstream: lag >= lead. First red:
    ``transpose`` (lead and lag exchanged)."""
    deg = 3 if m >= 4 else 1
    f, g = np.array(_F[: deg + 1]), np.array(_G[: deg + 1])
    rule = AttenuationRule.of(np.asarray(_TRI_XS), m).triangle()
    assert np.all(rule.lag >= rule.lead) and np.all(rule.weight >= 0.0)
    for i, x in enumerate(_TRI_XS):
        want = R.triangle(f, g, x)
        swapped = R.triangle(g, f, x)
        assert abs(want - swapped) >= 1e-6 * abs(want), "the pair is symmetric: the row is vacuous"
        got = float(np.sum(rule.weight[i] * np.polynomial.polynomial.polyval(1.0 - rule.lead[i], f)
                           * np.polynomial.polynomial.polyval(1.0 - rule.lag[i], g)))
        mass = 2.0 * min(2.0, 1.0 / max(x, 1e-300))
        bound = float(np.abs(_F).sum() * np.abs(_G).sum() * mass)
        assert abs(got - want) <= 32 * _EPS * bound, (x, got, want, abs(got - want) / (_EPS * bound))



# ── the law's two counting primitives (grading) ───────────────────────────────────────────────────────────────────

from orpheus.derivations.continuous.characteristic.grading import ALONG_LINE_ERROR, exact_nodes, gauss_nodes  # noqa: E402


@pytest.mark.foundation
def test_the_exact_node_count_is_the_fewest_gauss_nodes_exact_on_a_degree() -> None:
    """[LW1; foundation, THEOREM] ``exact_nodes(k)`` = ceil((k + 1)/2), the fewest nodes of a Gauss rule exact on
    degree k (2m - 1 >= k), for k = 0..80, checked against that definition and, at k <= 15, against the rule itself:
    the rule of ``exact_nodes(k)`` nodes integrates P_k exactly (RB3's band) and the rule of one fewer does not
    (P_{2m-1} with m - 1 nodes misses by O(1) of its moment at x = 0.7). First red: (k + 1) // 2."""
    k = np.arange(81)
    np.testing.assert_array_equal(exact_nodes(k), np.ceil((k + 1) / 2.0).astype(int))
    x = 0.7
    for degree in range(16):
        m = int(exact_nodes(degree))
        rule = AttenuationRule.of(np.array([x]), m)
        got = R.legendre_sums(rule.offset[0], rule.weight[0], degree + 1)[degree]
        assert abs(got - R.legendre_moment(degree, x)) <= _RB3_BAND * _EPS * _zeroth(x), degree
        if m > 1 and degree == 2 * m - 1:
            fewer = AttenuationRule.of(np.array([x]), m - 1)
            miss = R.legendre_sums(fewer.offset[0], fewer.weight[0], degree + 1)[degree] - R.legendre_moment(degree, x)
            assert abs(miss) > 1e-6 * _zeroth(x), (degree, miss)


@pytest.mark.foundation
def test_the_node_count_law_refuses_a_singularity_on_the_interval_or_a_bad_majorant() -> None:
    """[LW2; foundation] ``gauss_nodes`` refuses rho_sing <= 1 (the singularity on the interval: no ellipse inside
    it) and a majorant that is not finite and positive, each naming its reason; a valid call returns at least
    ``exact_nodes(1)`` and grows as rho_sing shrinks toward 1 (the positive legs). Its target defaults to
    ``ALONG_LINE_ERROR`` = 64 eps (the user's ruling of 2026-10-13). First reds: either refusal removed."""
    assert ALONG_LINE_ERROR == 64.0 * _EPS
    with pytest.raises(ValueError, match="it lies on the interval"):
        gauss_nodes(np.array([1.0]), lambda rho: np.ones_like(rho))
    with pytest.raises(ValueError, match="finite and positive"):
        gauss_nodes(np.array([3.0]), lambda rho: np.full_like(rho, np.inf))
    with pytest.raises(ValueError, match="finite and positive"):
        gauss_nodes(np.array([3.0]), lambda rho: np.zeros_like(rho))
    counts = gauss_nodes(np.array([10.0, 3.0, 1.2]), lambda rho: np.ones_like(rho))
    assert np.all(counts >= 1) and counts[0] < counts[1] < counts[2], counts
