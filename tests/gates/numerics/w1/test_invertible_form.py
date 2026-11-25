r"""G03s and G03t: positivity is asked by the verbs that need it, and a dense form on a factor has no singular state (N03b, spec 7.4b).

G03s. The factor requires an INVERTIBLE induced form. On a dense form it asks for no invertibility test of its own (a floating-point pivot test is not
exact and no tolerance decides structure): the form is positive semi-definite where the formal path built it, its kernel is its declared kernel, and on
such a form the kernel condition ``S ∩ ker G ⊆ K`` implies an invertible induced form on ``S/K``. The gate states that implication as a law over
two positive semi-definite dense forms with the declared kernel ``Coord(4, (3,))``, 324 states each (the 81 ``Coord`` pairs times the 4 variances):
152 construct, on every constructed state the live block of the induced density form (the indices ``S minus K`` of the density partner) is an
invertible symmetric matrix whose smallest eigenvalue magnitude is at least a quarter of its largest, and the factor reads ``signature[2] == 0`` and
``positive``. The fixtures are integer matrices whose live-block eigenvalues are of integer scale (the measured margin is about 0.33), so the count
and the margin are exact on any host and no row reads an eigenvalue that can be zero to rounding.

The ``negative`` and ``not_positive`` legs read the DECLARED signature of a dense form (its ``n_-`` is 0 by the declaration), so they can redden only by
breaking ``Factor.positive`` or ``Factor.signature`` themselves, not the dense law; the law's teeth are the singular-block count, the margin leg, and the
``n_+`` leg, which compares ``signature[0]`` with an eigenvalue count of the live block (exact on these forms).

G03t. ``FunctionSpace.norm`` is the verb that needs a positive form. On the measure space of ``lebedev(13)`` at ``L = 2`` (the signature
``(66, 0, 8)``: 8 negative weights) it raises ``ValueError`` naming the Krein arm and ``#706``; on ``lebedev(11)`` and on a counting space it returns
as before. The raise is a real ``raise`` that survives ``python -O`` (an ``ast`` leg reads the helper's body: a ``Raise`` node and no ``assert``). The refusal is reached through three more routes, each a leg: the dual of
the signed space (``DualMetric.negative`` reads its primal's), and the products of the signed space with a positive one in both block orders and with
itself (``FactoredMetric.negative`` reads the Kronecker signature, which a dropped term silences in the ``signed * positive`` order). The signed count
8 is also read from the rule's own negative weights, not from the metric.

Exact on any host (counts, index sets, a raise). At the base of N03b (``c18f156f``) ``Factor.signature`` did not exist (the signature leg of G03s
read ``AttributeError``) and ``norm`` returned ``3.5449077018110318`` on the signed space for the vector of ones (what it reads depends on the
vector: the inner product of other vectors can be negative and ``norm`` then reads ``nan``). Mutations: G03s, ``Factor._require_constructible``
made a no-op, which reddens the count (260 of 324 states construct, 108 of them with a singular live block); G03t, the ``raise`` replaced by an
``assert`` (the behavioural leg reddens under ``-O`` and the ``ast`` leg reddens on the source), or the call deleted from ``norm`` (the signed space
then returns).
"""

from __future__ import annotations

import ast
import inspect
from typing import Any

import numpy as np
import pytest

from tests.gates.numerics.w1 import _factor_probe as P

pytestmark = pytest.mark.foundation

#: Two positive semi-definite dense forms with the kernel ``Coord(4, (3,))`` and integer entries.
_PSD_FORMS = (
    ("PSD A", [[2.0, 1.0, 0.0, 0.0], [1.0, 2.0, 0.0, 0.0], [0.0, 0.0, 3.0, 0.0], [0.0, 0.0, 0.0, 0.0]]),
    ("PSD B", [[2.0, 1.0, 0.0, 0.0], [1.0, 5.0, 0.0, 0.0], [0.0, 0.0, 3.0, 0.0], [0.0, 0.0, 0.0, 0.0]]),
)
_KERNEL_INDEX = 3


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def _dense_axis(matrix: list[list[float]]) -> Any:
    from orpheus.numerics.axis import Axis, BasisKind, Declared
    from orpheus.numerics.declaration import Coord
    from orpheus.numerics.metric import DenseMetric, Face

    return Axis("h", (4,), origin=Declared(BasisKind.NODAL, DenseMetric(np.array(matrix), kernel=Coord(4, (_KERNEL_INDEX,)), face=Face.METRIC)))


def _readings(matrix: list[list[float]]) -> dict[str, Any]:
    """Construct the 324 states of one dense form and read each constructed one. RETURNS counts and the worst margin; asserts nothing."""
    from orpheus.numerics import factor as F

    axis = _dense_axis(matrix)
    everything = set(range(4))
    out: dict[str, Any] = {"states": 0, "constructed": 0, "singular": 0, "worst_margin": np.inf, "negative": 0, "not_positive": 0, "hand_agree": 0, "plus_disagree": 0}
    for S, K in P.pairs(4):
        sk = P.sk_of(4, S, K)
        for v in P.variances():
            pS, pK = (set(S), set(K)) if v.chart is F.VarianceChart.DENSITY else (everything - set(K), everything - set(S))
            live = sorted(pS - pK)
            # the hand verdict: the kernel index may lie in the partner's S only inside its K, and a dense form has no zero-dimensional state
            hand = (_KERNEL_INDEX not in pS or _KERNEL_INDEX in pK) and bool(live)
            out["states"] += 1
            try:
                f = F.Factor(axis, v, sk)
            except ValueError:
                out["hand_agree"] += not hand
                continue
            out["hand_agree"] += hand
            out["constructed"] += 1
            density = f.riesz() if v.chart is F.VarianceChart.INTEGRATED else f
            block = np.column_stack([np.asarray(density.form.apply(e)) for e in np.eye(4)])[np.ix_(live, live)]
            eig = np.linalg.eigvalsh((block + block.T) / 2.0)
            out["singular"] += bool(np.abs(eig).min() == 0.0)
            out["worst_margin"] = min(out["worst_margin"], float(np.abs(eig).min() / np.abs(eig).max()))
            signature = getattr(f, "signature", None)
            # the independent count: positive eigenvalues of the live block (exact on these integer-scale forms, no eigenvalue is near zero)
            out["plus_disagree"] += signature is None or signature[0] != int((eig > 0.0).sum())
            out["negative"] += signature is None or signature[2] != 0
            out["not_positive"] += getattr(f, "positive", None) is not True
    return out


def test_psd_dense_forms_have_no_singular_state() -> None:
    for label, matrix in _PSD_FORMS:
        r = _readings(matrix)
        _require(r["states"] == 324, f"{label}: {r['states']} states, expected 324")
        _require(r["constructed"] == 152, f"{label}: {r['constructed']} of 324 states construct, expected 152")
        _require(r["hand_agree"] == 324, f"{label}: {r['hand_agree']} of 324 verdicts agree with the hand rule")
        _require(r["singular"] == 0, f"{label}: {r['singular']} of 152 constructed states have a singular live block")
        _require(r["worst_margin"] >= 0.25, f"{label}: the smallest eigenvalue magnitude is {r['worst_margin']:.3f} of the largest, expected at least 0.25")
        _require(r["plus_disagree"] == 0, f"{label}: {r['plus_disagree']} constructed states read an n_+ that is not the live block's count of positive eigenvalues")
        _require(r["negative"] == 0, f"{label}: {r['negative']} constructed states read no signature or a negative direction")
        _require(r["not_positive"] == 0, f"{label}: {r['not_positive']} constructed states do not read positive")


def _signed_measure_space(order: int) -> Any:
    from orpheus.numerics import frame_faces as ff
    from orpheus.numerics.quadrature import Quadrature

    return ff.measure_space(ff.angular_frame(Quadrature.lebedev(order), 2))


def test_norm_raises_on_a_signed_diagonal_space() -> None:
    from orpheus.numerics import metric as M
    from orpheus.numerics.space import FunctionSpace

    from orpheus.numerics.quadrature import Quadrature

    signed = _signed_measure_space(13)
    resolved = signed._resolved_metric
    rule_negative = int((np.asarray(Quadrature.lebedev(13).weights) < 0.0).sum())
    _require(rule_negative == 8, f"the rule lebedev(13) has {rule_negative} negative weights, expected 8")
    _require(resolved.signature == (66, 0, rule_negative), f"the measure space of lebedev(13) reads signature {resolved.signature}, expected (66, 0, {rule_negative})")
    try:
        signed.norm(np.ones(signed.shape))
    except ValueError as e:
        _require("Krein" in str(e) and "#706" in str(e), f"the norm refusal does not name the Krein arm and #706: {e}")
    else:
        pytest.fail("norm returned on a space whose metric has 8 negative weights")
    positive = _signed_measure_space(11)
    _require(positive._resolved_metric.signature == (50, 0, 0), f"lebedev(11) reads signature {positive._resolved_metric.signature}")
    _require(np.isfinite(positive.norm(np.ones(positive.shape))), "norm on the positive lebedev(11) space is not finite")
    counting = FunctionSpace("counting", (5,))
    _require(counting.norm(np.ones(5)) == float(np.sqrt(5.0)), "norm on a counting space is not sqrt(5)")
    # the other routes to the same refusal: the dual space and the products in both block orders
    for label, space in (("dual", signed.dual()), ("signed*positive", signed * positive), ("positive*signed", positive * signed), ("signed*signed", signed * signed)):
        try:
            value = space.norm(np.ones(space.shape))
        except ValueError as e:
            _require("Krein" in str(e) and "#706" in str(e), f"{label}: the norm refusal does not name the Krein arm and #706: {e}")
        else:
            pytest.fail(f"{label}: norm returned {value} on a space with a negative direction")
    _require(np.isfinite((positive * positive).norm(np.ones((positive * positive).shape))), "norm on the product of two positive spaces is not finite")
    # a real raise, not an assert: the helper's body holds a Raise node and no Assert
    tree = ast.parse(inspect.getsource(M))
    helpers = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "_require_positive"]
    _require(len(helpers) == 1, f"{len(helpers)} definitions of _require_positive in metric.py, expected 1")
    body = list(ast.walk(helpers[0]))
    _require(any(isinstance(n, ast.Raise) for n in body), "_require_positive holds no raise")
    _require(not any(isinstance(n, ast.Assert) for n in body), "_require_positive holds an assert, which python -O strips")
