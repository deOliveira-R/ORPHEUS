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
* **Predicates** — rank-1 is STRUCTURALLY non-invertible (no ``solve``, no
  ``inverse()``; ``is_invertible`` False); adjointable exactly when the row
  is an :class:`InnerProductFunctional` (the dual dyad ``|w⟩⟨v|``).
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
from orpheus.numerics.operator import (
    MissingAdjoint,
    RankOneOperator,
    outer,
)

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


class TestPredicates:
    r"""Rank-1 has ``apply`` + ``apply_transpose`` (the dual dyad), never ``solve``/``inverse``.

    A dyad is **structurally non-invertible** (its kernel is the row's
    orthogonal complement along the contracted axis, so no ``solve`` and no
    ``inverse()`` are declared — misuse is a static error, Design C) but it
    HAS a **transpose**: ``(|v⟩⟨w|)ᵀ = |w⟩⟨v|``. The transpose exists iff the
    row is a genuine co-vector (:class:`InnerProductFunctional`) whose weight
    is the dual column; an opaque / nonlinear functional has no dual column,
    so its dyad refuses ``apply_transpose`` with :class:`MissingAdjoint` (the
    honest refusal). Campaign #276 (adjoint transport) added the transpose on
    the primitive — F† falls out of it.
    """

    def test_dyad_is_structurally_non_invertible(self):
        op = outer(np.array([[1.0], [2.0]]), InnerProductFunctional(np.array([[3.0], [4.0]]), axis=0))
        _require(
            callable(getattr(op, "apply", None)),
            "a rank-1 dyad must expose apply.",
        )
        _require(not op.is_invertible, "a rank-1 operator must NOT be invertible.")
        _require(
            not hasattr(op, "solve"),
            "a rank-1 operator must not declare solve (structural refusal is static).",
        )
        _require(
            not hasattr(op, "inverse"),
            "a rank-1 operator must not declare inverse() (structural refusal is static).",
        )

    def test_inner_product_dyad_is_adjointable(self):
        """The IPF-rowed dyad must ADVERTISE its working transpose: the
        ``is_adjointable`` predicate is the runtime successor of the
        pre-carve ``CAP_APPLY_TRANSPOSE`` membership this file pinned."""
        op = outer(np.array([[1.0], [2.0]]), InnerProductFunctional(np.array([[3.0], [4.0]]), axis=0))
        _require(
            op.is_adjointable,
            "a rank-1 dyad with an InnerProductFunctional row must be "
            "adjointable (the dual dyad |w⟩⟨v| exists and apply_transpose works).",
        )

    def test_opaque_functional_dyad_is_apply_only(self):
        # A row with no dual column (NOT an InnerProductFunctional) yields an
        # apply-only dyad — apply_transpose is honestly unavailable (it raises
        # MissingAdjoint, the ADJOINT-axis refusal).
        class _OpaqueFunctional:
            def evaluate(self, x, /):
                return np.asarray(x).sum(axis=0, keepdims=True)

        op = outer(np.array([[1.0], [2.0]]), _OpaqueFunctional())
        _require(
            not op.is_adjointable,
            "an opaque-functional dyad has no dual column → must NOT be adjointable.",
        )
        _require(not op.is_invertible, "an opaque-functional dyad is still non-invertible.")
        with pytest.raises(MissingAdjoint):
            op.apply_transpose(np.array([[1.0], [1.0]]))


class TestDyadTranspose:
    r"""``(|v⟩⟨w|)ᵀ = |w⟩⟨v|`` — the dual-dyad transpose law (campaign #276, A1)."""

    def test_apply_transpose_is_dual_dyad(self):
        # Forward |v⟩⟨w| x = v · ⟨w, x⟩; transpose |w⟩⟨v| y = w · ⟨v, y⟩, the
        # column and the row swapped, contracting over the SAME axis.
        v = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])  # (2, 3) column
        w_arr = np.array([[0.5, 1.0, 2.0], [1.5, 0.0, 0.25]])  # (2, 3) row weight
        op = outer(v, InnerProductFunctional(w_arr, axis=0))
        y = np.array([[2.0, 1.0, 0.5], [3.0, 4.0, 8.0]])

        got = op.apply_transpose(y)
        # structurally-independent: w · (hand Python-loop dot of v with y).
        ref = w_arr * _hand_dot_axis0(v, y)[None, :]
        np.testing.assert_array_almost_equal_nulp(got, ref, nulp=8)
        # the transpose IS the dual dyad outer(w, ⟨v|).apply — single source of truth.
        dual = outer(w_arr, InnerProductFunctional(v, axis=0))
        np.testing.assert_array_equal(got, dual.apply(y))

    def test_transpose_defining_inner_product_identity(self):
        # ⟨A x, y⟩ == ⟨x, Aᵀ y⟩ (full Euclidean) — the transpose-DEFINING identity.
        rng = np.random.default_rng(20260628)
        v = rng.uniform(size=(3, 4))
        w_arr = rng.uniform(size=(3, 4))
        op = outer(v, InnerProductFunctional(w_arr, axis=0))
        x = rng.uniform(size=(3, 4))
        y = rng.uniform(size=(3, 4))
        lhs = float((op.apply(x) * y).sum())            # ⟨A x, y⟩
        rhs = float((x * op.apply_transpose(y)).sum())  # ⟨x, Aᵀ y⟩
        np.testing.assert_allclose(lhs, rhs, rtol=1e-12)


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
