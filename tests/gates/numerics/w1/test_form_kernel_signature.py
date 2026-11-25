r"""G03e: a form reads its kernel and its signature exactly, and a factor admits a signed diagonal form wherever its induced form is invertible (N03b, spec 7.4b).

The kernel of ``Diag([1, 0, 2])`` is ``Coord(3, (1,))`` exactly and ``signature == (2, 1, 0)``; a signed ``Diag([1, -1, 2])`` reads
``signature == (2, 0, 1)``. A ``Factor`` on the signed form is constructed for every ``Coord`` pair ``(S, K)`` whose kernel condition holds: the
form's kernel is its exact zero set, and ``Diag([1, -1, 2])`` has none. ``Factor.signature`` is ``(n_+, n_0, n_-)`` read on the induced form, the
weights on ``S minus K`` and zero elsewhere, and ``Factor.positive`` is ``signature[2] == 0``: ``S = (0, 1, 2)`` with ``K = ()`` reads
``(2, 0, 1)`` and not positive; with ``K = (1,)`` it reads ``(2, 1, 0)`` and positive; ``S = (0, 1)`` with ``K = (1,)`` reads ``(1, 2, 0)``; and
``S = (0, 2)`` with ``K = ()`` reads ``(2, 1, 0)`` and positive, its induced weights ``[1, 0, 2]``. A formless factor has no ``signature``
(``AttributeError``). ``DiagonalMetric.counting((2, 3))`` has weights ``ones((2, 3))``; a ``DenseMetric`` with a ``Spanned`` kernel of
dimension 1 on 3 coordinates has ``signature == (2, 1, 0)``.

Exact on any host (signs, counts and index sets). Every expected value is written here by hand from the definitions. First red at the base of
N03b: the first three states are refused (a negative weight on ``S`` is refused by ``Factor``) and ``Factor.signature`` does not exist. Mutations:
reinstalling the negative-weight refusal in ``Factor._require_constructible`` reddens the first state; ``signature`` counting a negative weight
as positive reddens the ``(2, 0, 1)`` row; ``positive`` hard-wired to True reddens the first state; the zero weights read as nonzero in the
count reddens the ``(2, 1, 0)`` rows; the density partner's ``S`` read from ``self.sk.S`` instead of the partner's reddens the signed diagonal
population below on the integrated variances only.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest

from tests.gates.numerics.w1 import _factor_probe as P

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def test_form_kernel_and_signature() -> None:
    from orpheus.numerics import factor as F
    from orpheus.numerics.axis import Axis, BasisKind, Declared, EnergyAxis
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
    # (S, K, the signature and the positive flag by hand: n_+ and n_- count the signs of the weights on S minus K, n_0 the other indices)
    rows = (
        ((0, 1, 2), (), (2, 0, 1), False),
        ((0, 1, 2), (1,), (2, 1, 0), True),
        ((0, 1), (1,), (1, 2, 0), True),
        ((0, 2), (), (2, 1, 0), True),
    )
    for S, K, want_signature, want_positive in rows:
        try:
            f = F.Factor(ax, v, P.sk_of(3, S, K))
        except ValueError as e:
            pytest.fail(f"a factor on the signed diagonal form with S={S}, K={K} was refused: {e}")
        got = getattr(f, "signature", None)
        _require(got == want_signature, f"S={S}, K={K}: signature {got}, expected {want_signature}")
        _require(getattr(f, "positive", None) is want_positive, f"S={S}, K={K}: positive {getattr(f, 'positive', None)}, expected {want_positive}")
    induced = F.Factor(ax, v, P.sk_of(3, (0, 2), ())).form.apply(np.ones(3))
    _require(np.array_equal(np.asarray(induced), np.array([1.0, 0.0, 2.0])), f"the induced weights of S=(0, 2), K=() read {induced}, expected [1, 0, 2]")
    _require(float(np.sum(induced)) == 3.0, f"the three induced entries sum to {np.sum(induced)}, expected 3")
    formless = F.Factor(EnergyAxis.synthetic(3), F.FormlessVariance.MEASURE, F.Subquotient.whole(3))
    try:
        formless.signature
    except AttributeError:
        pass
    else:
        pytest.fail("a formless factor has a signature")
    formed = Axis("s", (3,), origin=Declared(BasisKind.NODAL))
    _require(F.Factor(formed, v, F.Subquotient.whole(3)).signature == (3, 0, 0), "a counting factor does not read (3, 0, 0)")


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


def _hand_rule(weights: tuple[float, ...], sk: Any, variance: Any, F: Any) -> tuple[bool, np.ndarray | None, tuple[int, int, int] | None, bool | None]:
    r"""``(constructs, induced density weights, signature, positive)`` of one state, from the definitions and with no production predicate.

    The density partner of the state is ``(S, K)`` on a density chart and ``(Ann K, Ann S)`` on an integrated chart (complements on the index
    set). The state constructs iff no zero-weight index lies in the partner's ``S`` outside its ``K`` (``S ∩ ker G ⊆ K``, the kernel of a
    diagonal form being its exact zero set); the induced weights are the weights on ``S minus K`` and zero elsewhere.
    """
    n = len(weights)
    everything = set(range(n))
    S, K = set(sk.S.indices), set(sk.K.indices)
    pS, pK = (S, K) if variance.chart is F.VarianceChart.DENSITY else (everything - K, everything - S)
    if any(weights[i] == 0.0 and i not in pK for i in pS):
        return False, None, None, None
    live = sorted(pS - pK)
    induced = np.array([weights[i] if i in live else 0.0 for i in range(n)])
    plus = sum(1 for i in live if weights[i] > 0.0)
    minus = sum(1 for i in live if weights[i] < 0.0)
    return True, induced, (plus, n - plus - minus, minus), minus == 0


def _population(weights: tuple[float, ...]) -> dict[str, int]:
    r"""Compare production with the hand rule on the 81 ``Coord`` pairs times the 4 variances of an axis with these weights.

    Counts: states, constructed (production), refused, hand-constructed, and the states on which production agrees with the hand rule on the
    verdict, on the induced form (the bytes of ``form.apply(ones)``, read on the density end and as reciprocals on the integrated end), on
    the signature, on the positive flag, and on the refusal being a ``ValueError``; and the constructed states on which ``dual`` and ``riesz``
    are involutions (L3).
    """
    from orpheus.numerics import factor as F

    axis = P.diag_axis(list(weights))
    out = dict.fromkeys(("states", "constructed", "refused", "hand", "verdict", "form", "signature", "positive", "dual_dual", "riesz_riesz", "flips"), 0)
    for S, K in P.pairs(len(weights)):
        sk = P.sk_of(len(weights), S, K)
        for v in P.variances():
            want, induced, want_sig, want_pos = _hand_rule(weights, sk, v, F)
            out["states"] += 1
            out["hand"] += want
            try:
                f = F.Factor(axis, v, sk)
            except ValueError:
                out["refused"] += 1
                out["verdict"] += not want
                continue
            out["constructed"] += 1
            out["verdict"] += want
            if not want or induced is None:
                continue
            hand_form = induced if v.chart is F.VarianceChart.DENSITY else np.array([1.0 / x if x != 0.0 else 0.0 for x in induced])
            out["form"] += bool(np.array_equal(np.asarray(f.form.apply(np.ones(len(weights)))), hand_form))
            out["signature"] += getattr(f, "signature", None) == want_sig
            out["positive"] += getattr(f, "positive", None) is want_pos
            out["dual_dual"] += f.dual().dual() == f
            out["riesz_riesz"] += f.riesz().riesz() == f
            # the flips by hand: dual reverses kind and chart, riesz reverses the chart only, and both carry the subquotient (Ann K, Ann S)
            d, r = f.dual(), f.riesz()
            everything = set(range(len(weights)))
            want_sk = (everything - set(sk.K.indices), everything - set(sk.S.indices))
            out["flips"] += (
                d.kind is not v.kind
                and d.chart is not v.chart
                and r.kind is v.kind
                and r.chart is not v.chart
                and all((set(x.sk.S.indices), set(x.sk.K.indices)) == want_sk for x in (d, r))
            )
    return out


def test_the_signed_diagonal_population_constructs_exactly_where_the_hand_rule_says() -> None:
    """The signed-diagonal population on both ends: the axis with weights ``[1, -1, 2, 0]``, every ``Coord`` pair ``K ⊆ S`` on its 4 indices (81
    pairs) times the 4 variances (324 states), each compared with the hand rule above. The negative weight at index 1 refuses nothing: 216 states
    construct and 108 are refused, exactly the states a zero weight (index 3) lies in the partner's ``S`` outside its ``K``. On every constructed
    state the induced form, the signature and the ``positive`` flag equal the hand rule, and ``dual`` and ``riesz`` are involutions that flip what the hand rule says (L3: ``dual`` reverses the kind and the chart, ``riesz`` the chart, both take the subquotient to ``(Ann K, Ann S)``). The
    positive control ``[1, 2, 3, 0]`` has the same zero index and reads the same 216 and 324 of 324. Exact: index sets, signs and the bytes of
    reciprocals of small weights. Mutation: a negative-weight refusal on ``S`` of the density partner reads 180 of 324 verdicts (the tip) and a
    refusal on ``self.sk.S`` disagrees on the integrated states."""
    for weights in ((1.0, -1.0, 2.0, 0.0), (1.0, 2.0, 3.0, 0.0)):
        c = _population(weights)
        label = f"weights {list(weights)}"
        _require(c["states"] == 324, f"{label}: {c['states']} states, expected 324")
        _require(c["hand"] == 216, f"{label}: the hand rule constructs {c['hand']} states, expected 216")
        _require(
            c["verdict"] == 324 and c["constructed"] == 216 and c["refused"] == 108,
            f"{label}: {c['verdict']} of 324 verdicts agree with the hand rule ({c['constructed']} constructed, {c['refused']} refused)",
        )
        _require(c["form"] == 216, f"{label}: the induced form agrees on {c['form']} of 216 constructed states")
        _require(c["signature"] == 216, f"{label}: the signature agrees on {c['signature']} of 216 constructed states")
        _require(c["positive"] == 216, f"{label}: the positive flag agrees on {c['positive']} of 216 constructed states")
        _require(
            c["dual_dual"] == 216 and c["riesz_riesz"] == 216 and c["flips"] == 216,
            f"{label}: dual o dual = id on {c['dual_dual']}, riesz o riesz = id on {c['riesz_riesz']} and the flips equal the hand rule on {c['flips']} of 216",
        )
