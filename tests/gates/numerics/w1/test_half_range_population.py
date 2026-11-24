r"""G02b: no shipped (rule, axis) pair is refused or misread by the exact half-range (N02, spec 7.3).

The population is the 17 rules of the census, each at axes 0, 1, 2 (51 pairs), plus a rotated-rule stub whose cosines are round-off
sized. The property holds on today's tree (the exact classifier is not an input of any shipped rule), so this gate is green by design
and its first red is its mutation: a classifier that reads ``|cos| < 1e-12`` as tangential reddens the stub row, and no shipped
``|cos|`` lies below its measured floor, so the stub is the only input that can red a tolerance.

Every row reads integer index sets and a ``min`` over cosines compared with fixed bounds ``1e-10`` and ``1e-3``: no tolerance of this
gate is near a value that a different LAPACK or libm could move (the rules are built from exact arithmetic on roots of unity and
Gauss nodes, and the floor sits two decades above round-off); each row is host-independent.

First red today: import-only (``half_range`` does not exist). Behavioural witness: the reflection-fixed-point owner
(``ordinate_pairing``) refuses 3 of the 29 genuine pairs of the ``q1`` population, which is why it is not the half-range owner.
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.gates import _half_range_helpers as H

pytestmark = pytest.mark.foundation

_ROTATED = [-0.5, 1e-13, -1e-13, 0.0, 0.5]
_ROTATED_WANT = ((1, 4), (0, 2), (3,))


def _owner():
    from orpheus.numerics.tangential import half_range

    return half_range


def _pairs():
    r"""``(name, axis, quadrature)`` for the 51 pairs; each rule is built once."""
    return [(name, axis, quad) for name, make in H.census_rules() for quad in (make(),) for axis in range(3)]


def kernel_minus_tangential(cosines: np.ndarray, weights: np.ndarray) -> list[int]:
    r"""Indices where the zero set of ``|c| w`` differs from the zero set of ``c``: the kernel of ``Diag(|cos| w)`` against the half-range's tangential set."""
    kernel = np.abs(cosines) * weights == 0.0
    tangential = cosines == 0.0
    return np.flatnonzero(kernel != tangential).tolist()


def test_the_population_is_the_census() -> None:
    pairs = _pairs()
    assert len(pairs) == 51, f"{len(pairs)} pairs, wanted 51 (17 rules, 3 axes)"
    assert len({name for name, _, _ in pairs}) == 17
    assert set(H.Q1_RULES) <= {name for name, _, _ in pairs}


def test_no_shipped_rule_axis_pair_is_refused() -> None:
    r"""``half_range`` partitions the ordinates for both signs on 51 of 51 pairs; ``tangential`` is exactly ``mu == 0.0``; the 8 identically
    zero pairs (``gauss_legendre`` at axes 1 and 2) return an all-tangential ``HalfRange``."""
    half_range = _owner()
    all_zero: list[tuple[str, int]] = []
    genuine: list[tuple[str, int]] = []
    bad: list[str] = []
    for name, axis, quad in _pairs():
        mu = np.asarray(quad.mean_axis_cosine(axis), dtype=np.float64)
        n = int(quad.n_ordinates)
        for sign in (+1, -1):
            half = half_range(quad, axis, sign=sign)
            bad += H.partition_violations(half, n, f"{name} axis {axis} sign {sign:+d}")
            if tuple(half.tangential.indices) != tuple(np.flatnonzero(mu == 0.0).tolist()):
                bad.append(f"{name} axis {axis}: tangential is not mu == 0.0")
        (genuine if np.any(mu) else all_zero).append((name, axis))
        if not np.any(mu) and len(half_range(quad, axis).tangential.indices) != n:
            bad.append(f"{name} axis {axis}: an identically zero pair is not all tangential")
    assert not bad, bad
    assert sorted(all_zero) == sorted((f"gl{n}", a) for n in (2, 4, 8, 16) for a in (1, 2)), all_zero
    assert len(genuine) == 43, len(genuine)


def test_the_q1_population_has_29_genuine_pairs() -> None:
    r"""The 11 rules of the qa probe ``q1`` give 33 pairs, 29 not identically zero: the population of the claim the pairing owner refuses 3 of."""
    genuine = 0
    total = 0
    for name, axis, quad in _pairs():
        if name in H.Q1_RULES:
            total += 1
            genuine += int(np.any(np.asarray(quad.mean_axis_cosine(axis))))
    assert (total, genuine) == (33, 29)


def test_the_round_off_gap_is_empty_and_the_floor_is_genuine() -> None:
    r"""The durable form of the retired band gate: no nonzero ``|cos|`` is ``<= 1e-10`` on the 51 pairs, and the smallest is ``> 1e-3``."""
    nonzero = np.concatenate(
        [np.abs(mu[mu != 0.0]) for _, axis, quad in _pairs() for mu in (np.asarray(quad.mean_axis_cosine(axis), dtype=np.float64),)]
    )
    assert nonzero.size > 0
    assert not np.any(nonzero <= 1e-10), f"{int(np.sum(nonzero <= 1e-10))} nonzero cosines lie in (0, 1e-10]"
    assert float(nonzero.min()) > 1e-3, f"smallest nonzero |cos| is {float(nonzero.min())!r}"


def test_the_kernel_of_the_trace_form_is_the_tangential_set() -> None:
    r"""0 of the 17 rules carries an exactly zero weight, so the zero set of ``|cos| w`` equals the tangential set on 51 of 51 pairs."""
    bad = []
    for name, axis, quad in _pairs():
        w = np.asarray(quad.weights, dtype=np.float64)
        assert not np.any(w == 0.0), f"{name} carries a zero weight"
        diff = kernel_minus_tangential(np.asarray(quad.mean_axis_cosine(axis), dtype=np.float64), w)
        if diff:
            bad.append((name, axis, diff))
    assert not bad, bad


def test_the_kernel_predicate_reads_a_zero_weight_on_a_genuine_ordinate() -> None:
    r"""Control: a synthetic rule with one zero weight on a non-tangential ordinate reads one disagreement; a clean rule reads none."""
    c = np.array([-0.5, 0.0, 0.5])
    assert kernel_minus_tangential(c, np.array([1.0, 1.0, 0.0])) == [2]
    assert kernel_minus_tangential(c, np.array([1.0, 1.0, 1.0])) == []
    assert kernel_minus_tangential(np.array([5e-324, 0.5]), np.array([0.25, 1.0])) == [0], "5e-324 * 0.25 underflows to 0"


def stub_violations(owner) -> list[str]:
    h = owner(H.CosineStub(_ROTATED), 0)
    got = (tuple(h.outflow.indices), tuple(h.inflow.indices), tuple(h.tangential.indices))
    return [] if got == _ROTATED_WANT else [f"rotated stub: {got}, wanted {_ROTATED_WANT}"]


def test_a_rotated_rule_is_read_exactly() -> None:
    r"""A round-off sized cosine of either sign is genuine: outflow ``(1, 4)``, inflow ``(0, 2)``, tangential ``(3,)``."""
    assert stub_violations(_owner()) == []


def test_the_rotated_stub_reads_a_tolerance_classifier_as_wrong() -> None:
    r"""Mutation: ``abs(mu) < 1e-12`` as tangential reddens the stub row (the cosines ``1e-13`` are below it). The ``1e-15`` cutoff of
    the retired selectors does not: it is below the stub, and no shipped cosine is near it either, which is why only G02a's round-off
    stub reads it."""
    assert stub_violations(H.banded_owner(1e-12)) != []
    assert stub_violations(H.banded_owner(1e-15)) == []
    assert stub_violations(H.exact_owner) == []
