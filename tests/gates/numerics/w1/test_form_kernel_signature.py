r"""G03e: a form reads its kernel and its signature exactly, and a factor refuses a negative weight where it matters (N03, spec 7.4).

The kernel of ``Diag([1, 0, 2])`` is ``Coord(3, (1,))`` exactly and ``signature == (2, 1, 0)``; a signed ``Diag([1, -1, 2])`` is ADMITTED and reads
``signature == (2, 0, 1)``; a ``Factor`` on it is refused with a message naming the Krein arm and #706 when the negative weight lies on ``S minus K``
(``S = (0, 1, 2)``, ``K = ()``) and when it lies inside ``K`` (``S = (0, 1, 2)``, ``K = (1,)``: the coset minimum over ``K`` of a negative weight is minus
infinity, so no induced form exists), and is ACCEPTED when it lies outside ``S`` (``S = (0, 2)``: the restriction does not read it; the orchestrator's P1 of
2026-11-25, amending R2).
``DiagonalMetric.counting((2, 3))`` has weights ``ones((2, 3))``; a ``DenseMetric`` with a ``Spanned`` kernel of dimension 1 on 3 coordinates has
``signature == (2, 1, 0)``.

Exact on any host (signs, counts and index sets). First red at ``ff305e0f``: a negative weight is accepted and ``rank((3,))`` reads 3; ``kernel``,
``signature`` and ``counting`` are ``AttributeError``s. Mutations: ``signature`` counts a negative weight as positive (reddens); ``Factor`` admits an
indefinite form (reddens the refusal leg); the refusal read on the whole axis form instead of ``S`` (reddens the accepting case); the refusal read on ``S minus K`` only (reddens the inside-``K`` case); the refusal read on ``self.sk.S`` instead of the density partner's ``S``
(reddens the signed-diagonal population below on the integrated variances only).
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.gates.numerics.w1 import _factor_probe as P

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def test_form_kernel_and_signature() -> None:
    from orpheus.numerics import factor as F
    from orpheus.numerics.declaration import Coord, Spanned
    from orpheus.numerics.metric import DenseMetric, DiagonalMetric, Face

    k = DiagonalMetric(np.array([1.0, 0.0, 2.0]))
    _require(k.kernel == Coord(3, (1,)), f"kernel {k.kernel}")
    _require(k.signature == (2, 1, 0), f"signature {k.signature}")
    signed = DiagonalMetric(np.array([1.0, -1.0, 2.0]))
    _require(signed.signature == (2, 0, 1), f"a signed diagonal form reads signature {signed.signature}, expected (2, 0, 1)")
    _require(np.array_equal(DiagonalMetric.counting((2, 3)).weights, np.ones((2, 3))), "counting((2, 3)) is not ones((2, 3))")
    g = np.diag([1.0, 1.0, 0.0])
    dense = DenseMetric(g, kernel=Spanned(3, ("k", 1), 1), face=Face.METRIC)
    _require(dense.signature == (2, 1, 0), f"a dense form with a Spanned kernel of dimension 1 reads signature {dense.signature}")
    ax = P.diag_axis([1.0, -1.0, 2.0])
    v = F.Variance(F.Kind.MEASURE, F.VarianceChart.DENSITY)
    try:
        F.Factor(ax, v, P.sk_of(3, (0, 1, 2), ()))
    except ValueError as e:
        _require("Krein" in str(e) and "706" in str(e), f"the refusal does not name the Krein arm and #706: {e}")
    else:
        pytest.fail("a factor with a negative weight on S minus K was constructed")
    for S, K in (((0, 1, 2), (1,)), ((0, 1), (1,))):
        try:
            F.Factor(ax, v, P.sk_of(3, S, K))
        except ValueError as e:
            _require("Krein" in str(e) and "706" in str(e), f"the refusal (S={S}, K={K}) does not name the Krein arm and #706: {e}")
        else:
            pytest.fail(f"a factor with the negative weight inside K was constructed (S={S}, K={K})")
    try:
        F.Factor(ax, v, P.sk_of(3, (0, 2), ()))
    except ValueError as e:
        pytest.fail(f"a negative weight outside S was refused: {e}")


def test_a_factored_metric_reads_the_signature_of_its_blocks() -> None:
    """A tensor product of forms reads ``signature`` and ``negative`` from its blocks (a product of two directions is negative when their signs
    differ): ``Diag([1, -1]) (x) Diag([1, 1])`` reads ``(2, 0, 2)`` and ``negative == 2``; ``Diag([1, -1]) (x) Diag([1, -1])`` reads ``(2, 0, 2)``;
    a zero weight multiplies out as a zero direction (``Diag([1, 0]) (x)`` a Euclidean 3-block reads ``(3, 3, 0)``), and a positive product reads
    ``negative == 0``. Exact on any host (integer counts). First red at the opening of this gate: ``negative`` read 0 and ``signature`` raised
    ``AttributeError``. Mutation: ``negative`` left at the base default 0 reddens the first leg."""
    from orpheus.numerics.metric import DiagonalMetric, FactoredMetric

    def product(*blocks):
        return FactoredMetric(tuple(((len(b.weights),), b) if b is not None else ((3,), None) for b in blocks))

    mixed = product(DiagonalMetric(np.array([1.0, -1.0])), DiagonalMetric(np.array([1.0, 1.0])))
    _require(mixed.signature == (2, 0, 2) and mixed.negative == 2, f"signed (x) positive reads {mixed.signature}, negative {mixed.negative}")
    both = product(DiagonalMetric(np.array([1.0, -1.0])), DiagonalMetric(np.array([1.0, -1.0])))
    _require(both.signature == (2, 0, 2), f"signed (x) signed reads {both.signature}")
    zero = product(DiagonalMetric(np.array([1.0, 0.0])), None)
    _require(zero.signature == (3, 3, 0), f"a zero direction (x) a Euclidean 3-block reads {zero.signature}")
    _require(product(DiagonalMetric(np.array([1.0, 2.0])), DiagonalMetric(np.array([3.0, 4.0]))).negative == 0, "a positive product reads negative directions")


def test_the_signed_diagonal_population_constructs_exactly_where_the_hand_rule_says() -> None:
    """The signed-diagonal population, on both ends: the axis with weights ``[1, -1, 2, 0]``, every ``Coord`` pair ``K ⊆ S`` on its 4 indices (81
    pairs) times the 4 variances (324 states), each compared with a rule written here from the definitions. The density partner's subquotient is
    ``(S, K)`` for a density chart and ``(Ann K, Ann S)`` (complements on the index set) for an integrated chart. The state constructs iff the
    partner is constructible (the zero-weight index 3 lies in the partner's ``S`` only inside its ``K``: ``S ∩ ker G ⊆ K``) and no weight on the
    partner's ``S`` is negative (index 1 is not in the partner's ``S``). A refusal for the negative weight names the Krein arm and #706. Exact:
    index sets and signs. Mutation: the refusal read on ``self.sk.S`` instead of the partner's ``S`` disagrees on the integrated states."""
    from orpheus.numerics import factor as F

    negative, zero = 1, 3
    constructed = refused = agree = 0
    for sk, axis in P.signed_diagonal_states():
        S, K = set(sk.S.indices), set(sk.K.indices)
        everything = set(range(4))
        for v in P.variances():
            pS, pK = (S, K) if v.chart is F.VarianceChart.DENSITY else (everything - K, everything - S)
            neg_in_s = negative in pS
            want = (zero not in pS or zero in pK) and not neg_in_s
            try:
                F.Factor(axis, v, sk)
                got, why = True, ""
            except ValueError as e:
                got, why = False, str(e)
            constructed += got
            refused += not got
            ok = got == want
            if not got and neg_in_s:
                ok = ok and "Krein" in why and "706" in why
            agree += ok
    _require(constructed + refused == 324, f"{constructed + refused} states, expected 324")
    _require(
        agree == 324,
        f"{agree} of 324 states agree with the hand rule ({constructed} constructed, {refused} refused)",
    )
    _require(constructed > 0 and refused > 0, f"the population is one-sided ({constructed} constructed, {refused} refused)")
