r"""G03l (moved from N03): slot-wise constructibility of a tensor of subquotients, exhaustively (N04, spec 7.5; law L7).

All 2187 pairs of coordinate subquotients of a 3-coordinate axis ``Diag[1, 2, 3]`` (27 states) and a 4-coordinate axis ``Diag[1, 2, 3, 0]`` (81 states, a kernel on the 4-axis only).
``of_axes(f1, f2).riesz(on=both)`` is accepted iff both slots construct, and the accepted set equals the independent HAND rule ``S meet ker G inside K`` per slot (written here from the
definitions, with its own set arithmetic and none of production's predicates); every accepted space has a nondegenerate product form ``G1 (x) G2`` (the hand rule on the quotient: no slot
has a nonzero kernel left); the pairs whose product form is nondegenerate but that the slot rule refuses (1971 of 2187 agree, 216 do not) all have a zero-dimensional factor.

First red at the base tip: import-only (``TensorSpace`` is absent). Mutation: ``Factor._require_constructible`` a no-op reddens (513 accepted spaces whose product form is degenerate, on the
spec's measurement); the fixture MUST carry a kernel: on ``Diag[1, 2, 3]`` and ``Diag[1, 2, 3, 4]`` the same mutation is green. Exact on any host: set arithmetic and refusals.
"""

from __future__ import annotations

import itertools

import pytest

from tests.gates.numerics.w1 import _factor_probe as FP
from tests.gates.numerics.w1 import _tensor_probe as P

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _slot_rule(S: tuple[int, ...], K: tuple[int, ...], kernel: set[int]) -> bool:
    """The hand rule: the induced form on S/K is nondegenerate iff S meet ker G lies inside K."""
    return (set(S) & kernel) <= set(K)


def _states(n: int) -> list[tuple[tuple[int, ...], tuple[int, ...]]]:
    return FP.pairs(n)


@pytest.mark.parametrize("w3, w4", [((1.0, 2.0, 3.0), (1.0, 2.0, 3.0, 0.0)), ((1.0, -2.0, 3.0), (1.0, 2.0, -3.0, 0.0))], ids=["positive", "signed"])
def test_slotwise_constructibility_is_the_set_rule(w3: tuple[float, ...], w4: tuple[float, ...]) -> None:
    from orpheus.numerics import factor as F

    ax3, ax4 = P.diag_axis(w3, "x"), P.diag_axis(w4, "y")
    ker3, ker4 = set(), {3}
    TS = P.tensor_space()
    s3, s4 = _states(3), _states(4)
    _require((len(s3), len(s4), len(s3) * len(s4)) == (27, 81, 2187), "the enumeration is not 27 x 81 = 2187")
    accepted = nondegenerate_product = agree = exceptions = 0
    zero_dim_exceptions = 0
    wrong: list[str] = []
    for (S1, K1), (S2, K2) in itertools.product(s3, s4):
        rule = _slot_rule(S1, K1, ker3) and _slot_rule(S2, K2, ker4)
        product_nondegenerate = (set(S1) == set(K1)) or (set(S2) == set(K2)) or rule
        built = True
        try:
            f1 = F.Factor(ax3, F.Variance(F.Kind.MEASURE, F.VarianceChart.DENSITY), FP.sk_of(3, S1, K1))
            f2 = F.Factor(ax4, F.Variance(F.Kind.MEASURE, F.VarianceChart.DENSITY), FP.sk_of(4, S2, K2))
            TS.of_axes(f1, f2).riesz(on={"x", "y"})
        except ValueError:
            built = False
        if built != rule:
            wrong.append(f"S1={S1} K1={K1} S2={S2} K2={K2}: built {built}, hand rule {rule}")
        accepted += int(built)
        if built:
            _require(product_nondegenerate, f"an accepted space has a degenerate product form: {S1}/{K1} x {S2}/{K2}")
            nondegenerate_product += 1
        agree += int(built == product_nondegenerate)
        if built != product_nondegenerate:
            exceptions += 1
            zero_dim_exceptions += int(set(S1) == set(K1) or set(S2) == set(K2))
    _require(not wrong, f"{len(wrong)} of 2187 pairs disagree with the hand rule: {wrong[:3]}")
    _require(accepted == 27 * 54, f"{accepted} accepted; the hand rule gives 27 x 54 = {27 * 54}")
    _require((agree, exceptions) == (1971, 216), f"the slot rule and the product form agree on {agree} and differ on {exceptions}; the spec states 1971 and 216")
    _require(zero_dim_exceptions == exceptions, f"{exceptions - zero_dim_exceptions} of the {exceptions} exceptions have no zero-dimensional factor")
