r"""``outer(reconstruction, functional)`` — the rank-1 dyad's defining laws (L0).

Intrinsic-property gate for the universal rank-1 constructor
:func:`~orpheus.numerics.operator.outer` (and its product type
:class:`~orpheus.numerics.operator.RankOneOperator`). Tests the dyad's
DEFINING laws, not merely its usage:

* **Action law** — ``outer(v, w).apply(x) == v * w.evaluate(x)`` (the matvec IS
  ``reconstruction * functional.evaluate``, routed through the functional).
* **Rank-1 structure** — the per-fiber output lies in ``span{v}``: for any two
  carriers, ``apply(x1)/⟨w,x1⟩ == apply(x2)/⟨w,x2⟩ == v``.
* **Bridge** — ``RankOneOperator(v, w)`` and ``outer(v, w)`` are the same dyad.
* **Capability** — ``{CAP_APPLY}`` only (rank-1 is non-invertible by
  construction).
* **Linearity** — the dyad is linear in its carrier (the row-factor is a genuine
  linear functional).

The value references are hand outer products / Python-loop dots (structurally
independent of the SUT where 0-ULP is not claimed). vv Mode-8: structural gates
via :func:`_require`, value gates via ``np.testing.*``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.functional import InnerProductFunctional
from orpheus.numerics.operator import CAP_APPLY, CAP_SOLVE, RankOneOperator, outer

pytestmark = pytest.mark.l0


def _require(condition: bool, message: str) -> None:
    if not condition:
        pytest.fail(message)


def _hand_dot_axis0(w: np.ndarray, x: np.ndarray) -> np.ndarray:
    """``(k,m) -> (m,)`` contraction over axis 0 by an explicit Python loop."""
    k, m = w.shape
    return np.array([sum(w[i, j] * x[i, j] for i in range(k)) for j in range(m)])


class TestDyadActionLaw:
    r"""``outer(v, w).apply(x) == v · ⟨w, x⟩``."""

    def test_action_equals_reconstruction_times_evaluate(self):
        v = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])  # (2, 3) column
        w_arr = np.array([[0.5, 1.0, 2.0], [1.5, 0.0, 0.25]])  # (2, 3) row weight
        x = np.array([[2.0, 1.0, 0.5], [3.0, 4.0, 8.0]])
        w = InnerProductFunctional(w_arr, axis=0)
        op = outer(v, w)

        got = op.apply(x)
        # 0-ULP: the matvec IS reconstruction * functional.evaluate.
        np.testing.assert_array_equal(got, v * w.evaluate(x))
        # structurally-independent: v * (hand Python-loop dot).
        ref = v * _hand_dot_axis0(w_arr, x)[None, :]
        np.testing.assert_array_almost_equal_nulp(got, ref, nulp=8)


class TestRankOneStructure:
    r"""The per-fiber output lies in ``span{v}`` — the defining rank-1 property."""

    def test_output_is_proportional_to_reconstruction(self):
        v = np.array([[1.0], [3.0], [5.0]])  # (3, 1) column, single fiber
        w = InnerProductFunctional(np.array([[2.0], [0.0], [1.0]]), axis=0)
        x1 = np.array([[1.0], [9.0], [2.0]])
        x2 = np.array([[4.0], [1.0], [7.0]])
        op = outer(v, w)

        inner1 = float(w.evaluate(x1).squeeze())
        inner2 = float(w.evaluate(x2).squeeze())
        _require(inner1 != 0.0 and inner2 != 0.0, "test inputs must give non-zero inner products.")
        # apply(x)/⟨w,x⟩ == v for every x → the column space is 1-D (= span{v}).
        np.testing.assert_allclose(op.apply(x1) / inner1, v, rtol=1e-13)
        np.testing.assert_allclose(op.apply(x2) / inner2, v, rtol=1e-13)


class TestOuterEqualsRankOneOperator:
    r"""``outer(v, w)`` is the :class:`RankOneOperator` of its factors."""

    def test_outer_returns_rank_one_operator(self):
        v = np.array([[1.0], [2.0]])
        w = InnerProductFunctional(np.array([[3.0], [4.0]]), axis=0)
        op = outer(v, w)
        _require(isinstance(op, RankOneOperator), "outer() must return a RankOneOperator.")
        _require(op.reconstruction is v, "the reconstruction must be the column passed to outer().")
        _require(op.functional is w, "the functional must be the row passed to outer().")
        # The explicit constructor and the verb agree.
        x = np.array([[5.0], [6.0]])
        np.testing.assert_array_equal(op.apply(x), RankOneOperator(v, w).apply(x))


class TestCapability:
    r"""Rank-1 advertises ``apply`` only — non-invertible by construction."""

    def test_capabilities_apply_only(self):
        op = outer(np.array([[1.0], [2.0]]), InnerProductFunctional(np.array([[3.0], [4.0]]), axis=0))
        _require(op.capabilities == frozenset({CAP_APPLY}), f"expected {{CAP_APPLY}}; got {op.capabilities}.")
        _require(CAP_SOLVE not in op.capabilities, "a rank-1 operator must NOT advertise solve.")


class TestLinearity:
    r"""The dyad is linear in its carrier."""

    def test_linear_in_carrier(self):
        v = np.array([[1.0], [2.0], [3.0]])
        w = InnerProductFunctional(np.array([[0.7], [0.2], [1.1]]), axis=0)
        op = outer(v, w)
        x1 = np.array([[1.0], [2.0], [3.0]])
        x2 = np.array([[0.5], [-1.0], [4.0]])
        a, b = 2.3, -0.7
        np.testing.assert_allclose(op.apply(a * x1 + b * x2), a * op.apply(x1) + b * op.apply(x2), rtol=1e-12)
