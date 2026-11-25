r"""Foundation tests for :mod:`orpheus.numerics.space`'s algebraic
constructions — :class:`TensorSpace`, :class:`DualSpace`, and
the :meth:`TensorSpace.__mul__` / :meth:`TensorSpace.dual` verbs.

Depth B step D-B (load-bearing). Pins the grand-report §6.1 / §15
tensor-product algebra at the L1 layer; consumed by Wave T (see
``.claude/plans/wave_t_tensor_network.md``) to type the codomains of
the boundary-realizer, fission, scattering, and streaming operators
that are rewired to ``TensorProductOperator`` and sums of them (``SumOfTensorProductsOperator``
until S2, now an ``OperatorSum`` of tensor products).

The invariants tested below are the type-system gates Wave T will
exercise — associativity of ``*``, shape composition, inner-product
factorisation, dual idempotency. Any failure here breaks the Wave T
operator-algebra rewires by construction.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.metric import FactoredMetric
from orpheus.numerics.space import (
    DualSpace,
    FunctionSpace,
    TensorSpace,
)
from orpheus.numerics.factor import Kind, Variance, VarianceChart
from orpheus.numerics.metric import DiagonalMetric


# ─────────────────────────────────────────────────────────────────────
# TensorSpace construction + identity
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_mul_returns_tensor_product_space():
    a = TensorSpace.declared("X", (4,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    b = TensorSpace.declared("G", (2,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    tp = a * b
    assert isinstance(tp, TensorSpace)
    assert tp.shape == (4, 2)
    assert tp.factors == a.factors + b.factors


@pytest.mark.foundation
def test_mul_three_factors_associativity_left():
    """``(A * B) * C`` flattens to a 3-factor product (no nesting)."""
    a = TensorSpace.declared("X", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    b = TensorSpace.declared("Omega", (8,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    c = TensorSpace.declared("G", (2,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    tp = (a * b) * c
    assert isinstance(tp, TensorSpace)
    assert tp.factors == a.factors + b.factors + c.factors
    assert tp.shape == (3, 8, 2)


@pytest.mark.foundation
def test_mul_three_factors_associativity_right():
    """``A * (B * C)`` also flattens to a 3-factor product."""
    a = TensorSpace.declared("X", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    b = TensorSpace.declared("Omega", (8,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    c = TensorSpace.declared("G", (2,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    tp = a * (b * c)
    assert isinstance(tp, TensorSpace)
    assert tp.factors == a.factors + b.factors + c.factors
    assert tp.shape == (3, 8, 2)


@pytest.mark.foundation
def test_mul_associativity_produces_identical_results():
    """``(A*B)*C`` and ``A*(B*C)`` compare equal at the space level."""
    a = TensorSpace.declared("X", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    b = TensorSpace.declared("Omega", (8,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    c = TensorSpace.declared("G", (2,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    left = (a * b) * c
    right = a * (b * c)
    assert left == right


@pytest.mark.foundation
def test_mul_name_format_joins_the_labels():
    """The name reads the slots in order — ``X(4,)*G(2,)`` for ``X * G`` (a tensor space's derived display name)."""
    a = TensorSpace.declared("X", (4,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    b = TensorSpace.declared("G", (2,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    tp = a * b
    assert tp.name.startswith("X(4,)*G(2,)#")


@pytest.mark.foundation
def test_mul_rejects_non_function_space():
    a = TensorSpace.declared("X", (4,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    with pytest.raises(TypeError):
        _ = a * 5  # int is not FunctionSpace


@pytest.mark.foundation
def test_the_empty_product_is_refused():
    a = TensorSpace.declared("X", (4,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    with pytest.raises(ValueError, match="at least one factor"):
        TensorSpace.of_axes()


# ─────────────────────────────────────────────────────────────────────
# Inner product factorisation
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_inner_product_factorises_euclidean():
    r"""For Euclidean factors,
    :math:`\langle x \otimes y, a \otimes b\rangle = \langle x, a\rangle \cdot \langle y, b\rangle`.

    With no inner-product weights, the TP space is Euclidean too.
    """
    a = TensorSpace.declared("X", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    b = TensorSpace.declared("G", (2,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    tp = a * b
    # Euclidean: TP weights stay None.
    assert tp.inner_product_weights is None
    # Element-level identity: <x⊗y, a⊗b> = (Σ x·a)(Σ y·b).
    x = np.array([1.0, 2.0, 3.0])
    y = np.array([4.0, 5.0])
    a_vec = np.array([2.0, 1.0, -1.0])
    b_vec = np.array([1.0, 3.0])
    # Tensor product values: shape (3, 2).
    lhs = tp.inner_product(np.outer(x, y), np.outer(a_vec, b_vec))
    rhs = a.inner_product(x, a_vec) * b.inner_product(y, b_vec)
    assert lhs == pytest.approx(rhs)


@pytest.mark.foundation
def test_inner_product_factorises_weighted():
    r"""For weighted factors, the tensor-product inner product is the
    factor inner products multiplied together.

    Pins the §15 identity that
    :math:`\langle \cdot, \cdot \rangle_{V_1 \otimes V_2}` factorises
    as
    :math:`\langle \cdot, \cdot \rangle_{V_1} \cdot \langle \cdot, \cdot \rangle_{V_2}`.
    """
    w_a = np.array([1.0, 2.0, 3.0])
    w_b = np.array([0.5, 0.5])
    a = TensorSpace.declared("X", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY), form=DiagonalMetric(w_a))
    b = TensorSpace.declared("G", (2,), Variance(Kind.FUNCTION, VarianceChart.DENSITY), form=DiagonalMetric(w_b))
    tp = a * b
    # The product carries its metric FACTORED (CS4c step 6 item 6.2a) —
    # never the outer product as a stored tensor; the VALUES it applies
    # are that outer product.
    assert tp.inner_product_weights is None
    assert isinstance(tp._resolved_metric, FactoredMetric)
    np.testing.assert_array_almost_equal(
        tp.apply_metric(np.ones((3, 2))), np.outer(w_a, w_b),
    )
    # Factorisation identity.
    x = np.array([1.0, 1.0, 1.0])
    y = np.array([1.0, 1.0])
    a_vec = np.array([2.0, 2.0, 2.0])
    b_vec = np.array([1.0, 1.0])
    lhs = tp.inner_product(np.outer(x, y), np.outer(a_vec, b_vec))
    rhs = a.inner_product(x, a_vec) * b.inner_product(y, b_vec)
    assert lhs == pytest.approx(rhs)


@pytest.mark.foundation
def test_inner_product_mixed_euclidean_and_weighted():
    """Mixing one Euclidean factor with one weighted factor produces a
    weighted TP (the Euclidean side contributes ones-of-shape)."""
    a = TensorSpace.declared("X", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))  # Euclidean
    w_b = np.array([1.0, 2.0])
    b = TensorSpace.declared("G", (2,), Variance(Kind.FUNCTION, VarianceChart.DENSITY), form=DiagonalMetric(w_b))
    tp = a * b
    assert tp.inner_product_weights is None
    assert isinstance(tp._resolved_metric, FactoredMetric)
    applied = tp.apply_metric(np.ones((3, 2)))
    # Each row equals w_b (since the Euclidean factor contributes 1).
    for i in range(3):
        np.testing.assert_array_almost_equal(applied[i], w_b)


# ─────────────────────────────────────────────────────────────────────
# Equality / hashing
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_two_tps_with_same_factors_compare_equal():
    a = TensorSpace.declared("X", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    b = TensorSpace.declared("G", (2,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    tp1 = a * b
    tp2 = TensorSpace.of_axes(*a.factors, *b.factors)
    assert tp1 == tp2


@pytest.mark.foundation
def test_tp_usable_as_dict_key():
    a = TensorSpace.declared("X", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    b = TensorSpace.declared("G", (2,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    tp = a * b
    d = {tp: "value"}
    assert d[a * b] == "value"


# ─────────────────────────────────────────────────────────────────────
# Dual space
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_dual_returns_dual_space():
    a = TensorSpace.declared("X", (4,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    a_dual = a.dual()
    assert isinstance(a_dual, TensorSpace)
    assert a_dual.shape == a.shape
    assert a_dual != a
    assert a_dual.factors == tuple(f.dual() for f in a.factors)


@pytest.mark.foundation
def test_dual_name_keeps_the_label_and_moves_the_digest():
    a = TensorSpace.declared("X", (4,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    assert a.dual().name.startswith("X(4,)#")
    assert a.dual().name != a.name


@pytest.mark.foundation
def test_dual_idempotent():
    r""":math:`V^{**} = V` (Riesz identification, double-dual is the primal)."""
    a = TensorSpace.declared("X", (4,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    a_doubledual = a.dual().dual()
    assert a_doubledual is a  # exact identity, not just equal


@pytest.mark.foundation
def test_the_dual_carries_the_inverse_weights():
    """The dual carries the inverse Riesz form G⁻¹ (S6a, ruling 26, #617),
    read through its primal: it copies no weights of its own, and its metric
    divides by w, bit for bit the primal's masked reciprocal. Until S6a this
    row pinned the opposite (``dual.inner_product_weights == w``, the dual
    carrying G).

    Units are NOT a space property (View-G, issues #205 / #207) — they
    live on the field role-leaf, so there is nothing unit-like for the
    dual to carry here.
    """
    w = np.array([1.0, 2.0, 3.0])
    a = TensorSpace.declared("X", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY), form=DiagonalMetric(w))
    a_dual = a.dual()
    assert a_dual.inner_product_weights is None
    x = np.array([1.0, 1.0, 1.0])
    np.testing.assert_array_equal(a_dual.apply_metric(x), x / w)


# ─────────────────────────────────────────────────────────────────────
# Repr
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_tp_repr_joins_the_factor_labels():
    a = TensorSpace.declared("X", (3,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    b = TensorSpace.declared("G", (2,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    r = repr(a * b)
    assert "TensorSpace" in r
    assert "X(3,)*G(2,)" in r


@pytest.mark.foundation
def test_dual_repr_names_its_class_and_label():
    a = TensorSpace.declared("X", (4,), Variance(Kind.FUNCTION, VarianceChart.DENSITY))
    r = repr(a.dual())
    assert "TensorSpace" in r
    assert "X(4,)" in r
