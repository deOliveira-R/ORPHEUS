r"""G03c: two variance types, one for a formless axis and one for a formed axis (N03, spec 7.4; ruled by the user, 2026-11-24).

``FormlessVariance`` has exactly the members ``FUNCTION`` and ``MEASURE``, ``dual`` swaps them, and it has no ``riesz``. ``Variance`` has the four
values. ``Factor(EnergyAxis, Variance(k, c), whole)`` raises ``TypeError`` for each of the four values; ``Factor(formed axis, FormlessVariance.<m>,
whole)`` raises ``TypeError``; all four ``Variance`` values construct on a formed axis. ``Factor(EnergyAxis, FormlessVariance.MEASURE, whole).riesz()``
raises ``AttributeError`` from the variance (the method exists on ``Factor``; the state is unspellable at the variance and refused at the factor
until the factor splits by form, issue #757). ``EnergyAxis.synthetic(3).form`` is ``NoForm`` and a formless ``Declared`` axis's form is
``Diag(ones)``: the counting default is one line of ``Axis.__post_init__``, which the energy arm's origin (it states ``NoForm``) never reaches.

Exact on any host. First red at ``ff305e0f``: import-only (``orpheus.numerics.factor``; ``EnergyAxis`` has no ``form``). Mutations: give
``FormlessVariance`` a ``riesz`` (reddens), drop the agreement of the variance type with the axis's form (the energy axis accepts a ``Variance``),
make the energy arm's ``Declared`` state no form (the counting default is reached and the form is ``Diag(ones)``).
"""

from __future__ import annotations

import numpy as np
import pytest

from tests.gates.numerics.w1 import _factor_probe as P

pytestmark = pytest.mark.foundation


def _require(cond: object, msg: str) -> None:
    if not cond:
        pytest.fail(msg)


def test_two_variance_types() -> None:
    from orpheus.numerics import factor as F
    from orpheus.numerics.axis import Axis, BasisKind, Declared, EnergyAxis
    from orpheus.numerics.metric import DiagonalMetric, NoForm

    e = EnergyAxis.synthetic(3)
    whole = F.Subquotient.whole(3)
    _require([m.name for m in F.FormlessVariance] == ["FUNCTION", "MEASURE"], "FormlessVariance does not have exactly FUNCTION and MEASURE")
    _require(not hasattr(F.FormlessVariance.FUNCTION, "riesz"), "FormlessVariance has a riesz")
    _require(F.FormlessVariance.FUNCTION.dual() is F.FormlessVariance.MEASURE and F.FormlessVariance.MEASURE.dual() is F.FormlessVariance.FUNCTION, "dual does not swap")
    _require(len(P.variances()) == 4 and len(set(P.variances())) == 4, "Variance does not have four distinct values")
    for v in P.variances():
        try:
            F.Factor(e, v, whole)
        except TypeError:
            pass
        else:
            pytest.fail(f"the energy axis accepted {v}")
    formed = Axis("s", (3,), origin=Declared(BasisKind.NODAL))
    for m in F.FormlessVariance:
        try:
            F.Factor(formed, m, whole)
        except TypeError:
            pass
        else:
            pytest.fail(f"a formed axis accepted FormlessVariance.{m.name}")
    for v in P.variances():
        F.Factor(formed, v, whole)  # all four construct
    f = F.Factor(e, F.FormlessVariance.MEASURE, whole)
    try:
        f.riesz()
    except AttributeError:
        pass
    else:
        pytest.fail("riesz on a formless factor returned")
    _require(f.dual().variance is F.FormlessVariance.FUNCTION, "dual of a formless MEASURE factor is not FUNCTION")
    _require(isinstance(e.form, NoForm), f"EnergyAxis.synthetic(3).form is {type(e.form).__name__}, not NoForm")
    counting = Axis("t", (3,), origin=Declared(BasisKind.NODAL)).form
    _require(isinstance(counting, DiagonalMetric) and np.array_equal(counting.weights, np.ones(3)), "a formless Declared axis does not read Diag(ones)")
