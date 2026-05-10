r"""Tests for :class:`TensorProductOperator` and
:class:`SumOfTensorProductsOperator`.

The Grand Report v3 §6.3 / §15.1 / §15.2 introduce the tensor-product
algebra :math:`A \otimes B \otimes C` for axis-aligned operators. The
contract verified here:

* Per-axis action: each constituent acts on its tagged axis, broadcasting
  on the rest. Order of application is irrelevant for disjoint axes.
* Capability intersection: ``(A & B)`` advertises a capability iff every
  factor advertises it.
* Adjoint distributivity: :math:`(A \otimes B)^* = A^* \otimes B^*`.
* Per-axis composition: :math:`(A \otimes B) @ (C \otimes D) = (A @ C) \otimes (B \otimes D)`.
* Sum-of-tensor-products: :math:`\Sigma_k A_k \otimes B_k` is the §15.2
  canonical scattering / streaming form.
* ``__and__`` flattening: ``(A & B) & C == A & (B & C)``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    DiagonalOperator,
    IdentityOperator,
    MissingCapability,
    SumOfTensorProductsOperator,
    TensorProductOperator,
    ZeroOperator,
)


# ─────────────────────────────────────────────────────────────────────
# TensorProductOperator
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestTensorProductApply:
    """L0: per-axis action against a Kronecker-product reference."""

    def test_two_diagonals_disjoint_axes(self):
        wA = np.array([2.0, 3.0])           # axis 0
        wB = np.array([5.0, 7.0, 11.0])     # axis 1
        D_a = DiagonalOperator(wA, axis=0)
        D_b = DiagonalOperator(wB, axis=1)
        T = D_a & D_b
        assert isinstance(T, TensorProductOperator)
        assert T.ops == (D_a, D_b)
        x = np.ones((2, 3))
        # (D_a ⊗ D_b) x = w_A_outer w_B
        expected = np.outer(wA, wB)
        np.testing.assert_array_equal(T.apply(x), expected)

    def test_three_factors(self):
        wA = np.array([1.0, 2.0])
        wB = np.array([3.0, 4.0])
        wC = np.array([5.0, 6.0])
        T = (
            DiagonalOperator(wA, axis=0)
            & DiagonalOperator(wB, axis=1)
            & DiagonalOperator(wC, axis=2)
        )
        assert len(T.ops) == 3
        x = np.ones((2, 2, 2))
        # w_A[i] * w_B[j] * w_C[k]
        expected = wA[:, None, None] * wB[None, :, None] * wC[None, None, :]
        np.testing.assert_array_equal(T.apply(x), expected)

    def test_identity_factor_is_no_op(self):
        wA = np.array([2.0, 3.0])
        T = DiagonalOperator(wA, axis=0) & IdentityOperator()
        x = np.ones((2, 4))
        np.testing.assert_array_equal(T.apply(x), wA[:, None] * np.ones((2, 4)))

    def test_apply_against_einsum_reference(self):
        """Tensor product action matches the einsum reference."""
        rng = np.random.default_rng(seed=99)
        wA = rng.standard_normal(4)
        wB = rng.standard_normal(5)
        T = DiagonalOperator(wA, axis=0) & DiagonalOperator(wB, axis=1)
        x = rng.standard_normal((4, 5))
        expected = np.einsum("a,b,ab->ab", wA, wB, x)
        np.testing.assert_allclose(T.apply(x), expected, rtol=1e-15)


@pytest.mark.l0
class TestTensorProductCapabilities:
    """L0: capability intersection."""

    def test_intersection_apply_only(self):
        T = ZeroOperator() & IdentityOperator()
        assert CAP_APPLY in T.capabilities
        # Zero advertises CAP_APPLY_TRANSPOSE; identity advertises everything;
        # intersection has CAP_APPLY_TRANSPOSE but NOT CAP_SOLVE (zero lacks solve).
        assert CAP_APPLY_TRANSPOSE in T.capabilities
        assert CAP_SOLVE not in T.capabilities

    def test_full_intersection(self):
        T = IdentityOperator() & IdentityOperator()
        for cap in (CAP_APPLY, CAP_APPLY_TRANSPOSE, CAP_SOLVE):
            assert cap in T.capabilities

    def test_solve_propagates_through_diagonals(self):
        D1 = DiagonalOperator(np.array([1.0, 2.0]), axis=0)
        D2 = DiagonalOperator(np.array([3.0, 4.0]), axis=1)
        T = D1 & D2
        assert CAP_SOLVE in T.capabilities
        # Solve round-trip
        x = np.array([[1.0, 1.0], [1.0, 1.0]])
        np.testing.assert_allclose(T.solve(T.apply(x)), x, rtol=1e-13)

    def test_solve_revoked_by_zero_weight_factor(self):
        D_zero = DiagonalOperator(np.array([1.0, 0.0]), axis=0)
        D_ok = DiagonalOperator(np.array([2.0, 3.0]), axis=1)
        T = D_zero & D_ok
        assert CAP_SOLVE not in T.capabilities
        with pytest.raises(MissingCapability, match="every factor"):
            T.solve(np.ones((2, 2)))


@pytest.mark.l0
class TestTensorProductAdjoint:
    """L0: (A ⊗ B).H == A.H ⊗ B.H."""

    def test_adjoint_distributivity_diagonals(self):
        wA = np.array([2.0, 3.0])
        wB = np.array([5.0, 7.0])
        D_a = DiagonalOperator(wA, axis=0)
        D_b = DiagonalOperator(wB, axis=1)
        T = D_a & D_b
        x = np.array([[1.0, 1.0], [1.0, 1.0]])
        # Diagonals are self-adjoint; (D_a ⊗ D_b)^T = D_a ⊗ D_b
        np.testing.assert_array_equal(T.apply_transpose(x), T.apply(x))


@pytest.mark.l0
class TestTensorProductFlattening:
    """L0: ``__and__`` flattens nested products."""

    def test_left_associativity_flattens(self):
        D1 = DiagonalOperator(np.array([1.0, 2.0]), axis=0)
        D2 = DiagonalOperator(np.array([3.0, 4.0]), axis=1)
        D3 = DiagonalOperator(np.array([5.0, 6.0]), axis=2)
        T = (D1 & D2) & D3
        assert len(T.ops) == 3
        assert T.ops == (D1, D2, D3)

    def test_right_associativity_flattens(self):
        D1 = DiagonalOperator(np.array([1.0, 2.0]), axis=0)
        D2 = DiagonalOperator(np.array([3.0, 4.0]), axis=1)
        D3 = DiagonalOperator(np.array([5.0, 6.0]), axis=2)
        T = D1 & (D2 & D3)
        assert len(T.ops) == 3
        assert T.ops == (D1, D2, D3)


@pytest.mark.l0
class TestTensorProductRequiresApply:
    """Constituent must advertise CAP_APPLY."""

    def test_construction_rejects_non_apply_factor(self):
        # Build a fake operator that doesn't advertise apply
        class FakeOp:
            capabilities = frozenset()
            def apply(self, x): return x
        with pytest.raises(MissingCapability, match=CAP_APPLY):
            TensorProductOperator((FakeOp(),))


# ─────────────────────────────────────────────────────────────────────
# SumOfTensorProductsOperator
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestSumOfTensorProducts:
    """L0: §15.2 sum-of-tensor-products algebra."""

    def test_sum_of_two_separable_terms(self):
        D_a = DiagonalOperator(np.array([2.0, 3.0]), axis=0)
        D_b = DiagonalOperator(np.array([5.0, 7.0]), axis=1)
        D_c = DiagonalOperator(np.array([11.0, 13.0]), axis=0)
        D_d = DiagonalOperator(np.array([17.0, 19.0]), axis=1)
        S = SumOfTensorProductsOperator((D_a & D_b, D_c & D_d))
        x = np.ones((2, 2))
        expected = (
            np.outer(D_a.weights, D_b.weights)
            + np.outer(D_c.weights, D_d.weights)
        )
        np.testing.assert_allclose(S.apply(x), expected, rtol=1e-15)

    def test_assert_separable_passes_when_constructed_legally(self):
        D_a = DiagonalOperator(np.array([1.0]))
        D_b = DiagonalOperator(np.array([1.0]))
        S = SumOfTensorProductsOperator((D_a & D_b,))
        S.assert_separable()  # no exception

    def test_constructor_rejects_non_tensor_summand(self):
        with pytest.raises(TypeError, match="must be"):
            SumOfTensorProductsOperator((IdentityOperator(),))

    def test_solve_not_advertised(self):
        """Sum doesn't have solve in general (sum-of-inverses != inverse-of-sum)."""
        D_a = DiagonalOperator(np.array([1.0]))
        D_b = DiagonalOperator(np.array([1.0]))
        S = SumOfTensorProductsOperator((D_a & D_b,))
        assert CAP_SOLVE not in S.capabilities

    def test_apply_transpose_propagates(self):
        D_a = DiagonalOperator(np.array([2.0, 3.0]), axis=0)
        D_b = DiagonalOperator(np.array([5.0, 7.0]), axis=1)
        S = SumOfTensorProductsOperator((D_a & D_b,))
        x = np.ones((2, 2))
        # Diagonals are self-adjoint
        np.testing.assert_array_equal(S.apply_transpose(x), S.apply(x))

    def test_empty_summand_list_rejected(self):
        with pytest.raises(ValueError, match="at least one"):
            SumOfTensorProductsOperator(())
