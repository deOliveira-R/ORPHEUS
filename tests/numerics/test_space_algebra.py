r"""Foundation tests for :mod:`orpheus.numerics.space`'s algebraic
constructions — :class:`TensorProductSpace`, :class:`DualSpace`, and
the :meth:`FunctionSpace.__mul__` / :meth:`FunctionSpace.dual` dunders.

Depth B step D-B (load-bearing). Pins the grand-report §6.1 / §15
tensor-product algebra at the L1 layer; consumed by Wave T (see
``.claude/plans/wave_t_tensor_network.md``) to type the codomains of
the boundary-realizer, fission, scattering, and streaming operators
that are rewired to ``TensorProductOperator`` / ``SumOfTensorProductsOperator``.

The invariants tested below are the type-system gates Wave T will
exercise — associativity of ``*``, shape composition, inner-product
factorisation, dual idempotency. Any failure here breaks the Wave T
operator-algebra rewires by construction.
"""
from __future__ import annotations

import numpy as np
import pint
import pytest

from orpheus.numerics.space import (
    DualSpace,
    FunctionSpace,
    TensorProductSpace,
)


# ─────────────────────────────────────────────────────────────────────
# Test fixtures
# ─────────────────────────────────────────────────────────────────────


@pytest.fixture(scope="module")
def ureg() -> pint.UnitRegistry:
    return pint.UnitRegistry()


# ─────────────────────────────────────────────────────────────────────
# TensorProductSpace construction + identity
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_mul_returns_tensor_product_space():
    a = FunctionSpace(name="X", shape=(4,))
    b = FunctionSpace(name="G", shape=(2,))
    tp = a * b
    assert isinstance(tp, TensorProductSpace)
    assert tp.shape == (4, 2)
    assert tp.factors == (a, b)


@pytest.mark.foundation
def test_mul_three_factors_associativity_left():
    """``(A * B) * C`` flattens to a 3-factor product (no nesting)."""
    a = FunctionSpace(name="X", shape=(3,))
    b = FunctionSpace(name="Omega", shape=(8,))
    c = FunctionSpace(name="G", shape=(2,))
    tp = (a * b) * c
    assert isinstance(tp, TensorProductSpace)
    assert tp.factors == (a, b, c)
    assert tp.shape == (3, 8, 2)


@pytest.mark.foundation
def test_mul_three_factors_associativity_right():
    """``A * (B * C)`` also flattens to a 3-factor product."""
    a = FunctionSpace(name="X", shape=(3,))
    b = FunctionSpace(name="Omega", shape=(8,))
    c = FunctionSpace(name="G", shape=(2,))
    tp = a * (b * c)
    assert isinstance(tp, TensorProductSpace)
    assert tp.factors == (a, b, c)
    assert tp.shape == (3, 8, 2)


@pytest.mark.foundation
def test_mul_associativity_produces_identical_results():
    """``(A*B)*C`` and ``A*(B*C)`` compare equal at the space level."""
    a = FunctionSpace(name="X", shape=(3,))
    b = FunctionSpace(name="Omega", shape=(8,))
    c = FunctionSpace(name="G", shape=(2,))
    left = (a * b) * c
    right = a * (b * c)
    assert left == right


@pytest.mark.foundation
def test_mul_name_format_uses_otimes():
    """The name reads as the math — ``X ⊗ G`` for ``X * G``."""
    a = FunctionSpace(name="X", shape=(4,))
    b = FunctionSpace(name="G", shape=(2,))
    tp = a * b
    assert "X" in tp.name and "G" in tp.name
    assert "⊗" in tp.name  # the actual U+2297 character


@pytest.mark.foundation
def test_mul_rejects_non_function_space():
    a = FunctionSpace(name="X", shape=(4,))
    with pytest.raises(TypeError):
        _ = a * 5  # int is not FunctionSpace


@pytest.mark.foundation
def test_from_factors_requires_at_least_two_factors():
    a = FunctionSpace(name="X", shape=(4,))
    with pytest.raises(ValueError, match="at least 2 factors"):
        TensorProductSpace.from_factors((a,))


# ─────────────────────────────────────────────────────────────────────
# Inner product factorisation
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_inner_product_factorises_euclidean():
    r"""For Euclidean factors,
    :math:`\langle x \otimes y, a \otimes b\rangle = \langle x, a\rangle \cdot \langle y, b\rangle`.

    With no inner-product weights, the TP space is Euclidean too.
    """
    a = FunctionSpace(name="X", shape=(3,))
    b = FunctionSpace(name="G", shape=(2,))
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
    a = FunctionSpace(name="X", shape=(3,), inner_product_weights=w_a)
    b = FunctionSpace(name="G", shape=(2,), inner_product_weights=w_b)
    tp = a * b
    # TP weights = outer product.
    assert tp.inner_product_weights is not None
    assert tp.inner_product_weights.shape == (3, 2)
    np.testing.assert_array_almost_equal(
        tp.inner_product_weights, np.outer(w_a, w_b),
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
    a = FunctionSpace(name="X", shape=(3,))  # Euclidean
    w_b = np.array([1.0, 2.0])
    b = FunctionSpace(name="G", shape=(2,), inner_product_weights=w_b)
    tp = a * b
    assert tp.inner_product_weights is not None
    assert tp.inner_product_weights.shape == (3, 2)
    # Each row equals w_b (since Euclidean factor contributes 1).
    for i in range(3):
        np.testing.assert_array_almost_equal(
            tp.inner_product_weights[i], w_b,
        )


# ─────────────────────────────────────────────────────────────────────
# Units composition
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_tp_units_multiply(ureg):
    r"""Units of :math:`V_1 \otimes V_2` are the product of factor units."""
    a = FunctionSpace(name="X", shape=(3,), units=ureg.Unit("cm**-2"))
    b = FunctionSpace(name="T", shape=(2,), units=ureg.Unit("s**-1"))
    tp = a * b
    assert tp.units is not None
    # Dimensionality is [length]^-2 · [time]^-1.
    expected_dim = (ureg.Unit("cm**-2") * ureg.Unit("s**-1")).dimensionality
    assert tp.units.dimensionality == expected_dim


@pytest.mark.foundation
def test_tp_units_none_if_any_factor_unitless(ureg):
    """If any factor lacks units, the TP product is None (cannot infer)."""
    a = FunctionSpace(name="X", shape=(3,), units=ureg.Unit("cm**-2"))
    b = FunctionSpace(name="G", shape=(2,))  # no units
    tp = a * b
    assert tp.units is None


# ─────────────────────────────────────────────────────────────────────
# Equality / hashing
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_two_tps_with_same_factors_compare_equal():
    a = FunctionSpace(name="X", shape=(3,))
    b = FunctionSpace(name="G", shape=(2,))
    tp1 = a * b
    tp2 = TensorProductSpace.from_factors((a, b))
    assert tp1 == tp2


@pytest.mark.foundation
def test_tp_usable_as_dict_key():
    a = FunctionSpace(name="X", shape=(3,))
    b = FunctionSpace(name="G", shape=(2,))
    tp = a * b
    d = {tp: "value"}
    assert d[a * b] == "value"


# ─────────────────────────────────────────────────────────────────────
# Dual space
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_dual_returns_dual_space():
    a = FunctionSpace(name="X", shape=(4,))
    a_dual = a.dual()
    assert isinstance(a_dual, DualSpace)
    assert a_dual.shape == a.shape
    assert a_dual.primal is a


@pytest.mark.foundation
def test_dual_name_appends_star():
    a = FunctionSpace(name="X", shape=(4,))
    assert "X" in a.dual().name
    assert "*" in a.dual().name


@pytest.mark.foundation
def test_dual_idempotent():
    r""":math:`V^{**} = V` (Riesz identification, double-dual is the primal)."""
    a = FunctionSpace(name="X", shape=(4,))
    a_doubledual = a.dual().dual()
    assert a_doubledual is a  # exact identity, not just equal


@pytest.mark.foundation
def test_dual_preserves_weights_and_units(ureg):
    """The L²-Riesz dual carries the same weights and units as the primal."""
    w = np.array([1.0, 2.0, 3.0])
    a = FunctionSpace(
        name="X", shape=(3,),
        inner_product_weights=w,
        units=ureg.Unit("cm**-2"),
    )
    a_dual = a.dual()
    np.testing.assert_array_equal(a_dual.inner_product_weights, w)
    assert a_dual.units == a.units


# ─────────────────────────────────────────────────────────────────────
# Repr
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.foundation
def test_tp_repr_includes_otimes_name():
    a = FunctionSpace(name="X", shape=(3,))
    b = FunctionSpace(name="G", shape=(2,))
    r = repr(a * b)
    assert "TensorProductSpace" in r
    assert "⊗" in r


@pytest.mark.foundation
def test_dual_repr_includes_star():
    a = FunctionSpace(name="X", shape=(4,))
    r = repr(a.dual())
    assert "DualSpace" in r
    assert "X*" in r
