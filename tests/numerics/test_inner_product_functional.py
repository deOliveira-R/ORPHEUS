r"""``InnerProductFunctional`` — the generic co-vector :math:`\langle w,\cdot\rangle` (L0).

Intrinsic-property gate (``feedback_test_intrinsic_properties``: every
math-bearing type ships a test of its defining laws). The laws:

* **Evaluate law** — ``evaluate(x) == (w*x).sum(axis, keepdims=True)``, pinned
  against an explicit Python-loop dot (a structurally-independent reference,
  no numpy reduction shared with the SUT).
* **Axis parametrisation** — contracting a different axis changes both the
  value and the surviving shape.
* **Linearity in the carrier** — the co-vector is a genuine *linear* functional.
* **Category membership** — it satisfies the :class:`Functional` Protocol and
  is NOT a :class:`LinearOperator` (no ``apply``/``capabilities``).

vv Mode-8: structural gates route through :func:`_require` (a function call,
fires under ``python -O``); value gates through ``np.testing.*``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.functional import Functional, InnerProductFunctional
from orpheus.numerics.operator import LinearOperator

pytestmark = pytest.mark.l0


def _require(condition: bool, message: str) -> None:
    """Fail with ``message`` if false. Fires under ``-O`` (not a bare assert)."""
    if not condition:
        pytest.fail(message)


def _hand_dot_axis0(w: np.ndarray, x: np.ndarray) -> np.ndarray:
    """``(n,m) -> (m,)`` contraction over axis 0 by an explicit Python loop."""
    n, m = w.shape
    return np.array([sum(w[k, j] * x[k, j] for k in range(n)) for j in range(m)])


def _hand_dot_axis1(w: np.ndarray, x: np.ndarray) -> np.ndarray:
    """``(n,m) -> (n,)`` contraction over axis 1 by an explicit Python loop."""
    n, m = w.shape
    return np.array([sum(w[i, k] * x[i, k] for k in range(m)) for i in range(n)])


class TestEvaluateLaw:
    """``evaluate(x) == ⟨w, x⟩`` over the contracted axis, keepdims."""

    def test_evaluate_matches_hand_loop_axis0(self):
        w = np.array([[2.0, 3.0, 5.0], [0.2, 0.3, 0.5]])  # (2, 3)
        x = np.array([[1.1, 1.3, 1.7], [2.2, 2.4, 2.6]])
        got = InnerProductFunctional(w, axis=0).evaluate(x)
        _require(got.shape == (1, 3), f"keepdims must leave (1,*comp); got {got.shape}.")
        np.testing.assert_allclose(got[0], _hand_dot_axis0(w, x), rtol=1e-13)

    def test_axis_parametrisation_changes_value_and_shape(self):
        w = np.arange(1, 13, dtype=float).reshape(3, 4)
        x = np.arange(1, 13, dtype=float).reshape(3, 4) * 0.1
        got0 = InnerProductFunctional(w, axis=0).evaluate(x)  # (1, 4)
        got1 = InnerProductFunctional(w, axis=1).evaluate(x)  # (3, 1)
        _require(
            got0.shape == (1, 4) and got1.shape == (3, 1),
            f"axis 0 → (1,4), axis 1 → (3,1); got {got0.shape}, {got1.shape}.",
        )
        np.testing.assert_allclose(got0[0], _hand_dot_axis0(w, x), rtol=1e-13)
        np.testing.assert_allclose(got1[:, 0], _hand_dot_axis1(w, x), rtol=1e-13)


class TestLinearity:
    """The functional is linear in its carrier argument."""

    def test_linear_in_carrier(self):
        w = np.array([[0.7], [0.2], [1.1]])
        x1 = np.array([[1.0], [2.0], [3.0]])
        x2 = np.array([[0.5], [-1.0], [4.0]])
        f = InnerProductFunctional(w, axis=0)
        a, b = 2.3, -0.7
        lhs = f.evaluate(a * x1 + b * x2)
        rhs = a * f.evaluate(x1) + b * f.evaluate(x2)
        np.testing.assert_allclose(lhs, rhs, rtol=1e-12)


class TestCategoryMembership:
    """A Functional, not a LinearOperator — the §5.6 category partition."""

    def test_satisfies_functional_not_linear_operator(self):
        f = InnerProductFunctional(np.array([[1.0], [2.0]]), axis=0)
        _require(
            isinstance(f, Functional),
            "InnerProductFunctional must satisfy the Functional Protocol (evaluate).",
        )
        _require(
            not isinstance(f, LinearOperator),
            "A Functional must NOT satisfy LinearOperator (it has no apply/"
            "capabilities surface) — the category partition is structural.",
        )
