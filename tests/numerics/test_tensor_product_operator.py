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
    RankOneOperator,
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


# ─────────────────────────────────────────────────────────────────────
# RankOneOperator (Wave T step T.2 — fission kernel primitive)
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestRankOneOperator:
    """L0: rank-1 outer product :math:`|\\ell\\rangle\\langle r|` on a
    tagged axis with spatial-parametrised vectors.  Introduced for the
    multigroup fission emission :math:`F = \\chi \\otimes \\nu\\Sigma_f`
    (Grand Report v3 §15 line 2008).
    """

    def test_apply_matches_hand_computed_outer_product(self):
        """1-D case: ``out_i = ell_i · sum_j r_j · x_j`` (the pure
        :math:`|\\ell\\rangle\\langle r|` action)."""
        ell = np.array([1.0, 2.0, 3.0])
        r = np.array([4.0, 5.0, 6.0])
        x = np.array([0.5, -0.25, 0.75])
        R = RankOneOperator(ell, r, axis=0)
        inner = float(r @ x)  # = 4·0.5 + 5·(-0.25) + 6·0.75 = 5.25
        expected = ell * inner
        np.testing.assert_array_equal(R.apply(x), expected)

    def test_apply_multi_dim_spatial_parametrisation(self):
        """3-D case (group × spatial × spatial) — fission's actual
        ``(ng, nx, ny)`` shape.  ``out[g, x, y] = chi[g, x, y] ·
        sum_g' (nu_sigma_f[g', x, y] · phi[g', x, y])``.

        ``np.einsum`` and ``.sum(axis=0)`` may pick different
        pairwise-reduction trees (1 ULP drift over a depth-``ng``
        sum); compared at ``nulp=4`` per the project's
        algebra-of-record FP discipline.
        """
        ng, nx, ny = 4, 3, 2
        rng = np.random.default_rng(11)
        chi = rng.uniform(0.1, 1.0, size=(ng, nx, ny))
        nu_sigma_f = rng.uniform(0.01, 0.5, size=(ng, nx, ny))
        phi = rng.uniform(0.1, 2.0, size=(ng, nx, ny))

        R = RankOneOperator(chi, nu_sigma_f, axis=0)
        out = R.apply(phi)

        # Hand-computed per-cell rank-1 action.
        inner = np.einsum("gxy,gxy->xy", nu_sigma_f, phi)
        expected = chi * inner[None, :, :]
        np.testing.assert_array_almost_equal_nulp(out, expected, nulp=4)

    def test_capabilities_apply_only(self):
        """Rank-1 ops advertise CAP_APPLY only (non-invertible by
        construction; ``apply_transpose`` not yet a consumer)."""
        R = RankOneOperator(np.ones(3), np.ones(3))
        assert R.capabilities == frozenset({CAP_APPLY})

    def test_constructor_rejects_shape_mismatch(self):
        with pytest.raises(ValueError, match="matching shapes"):
            RankOneOperator(np.ones((3,)), np.ones((4,)))

    def test_constructor_rejects_scalar(self):
        with pytest.raises(ValueError, match="ScaledOperator"):
            RankOneOperator(np.array(1.0), np.array(2.0))

    def test_constructor_rejects_axis_out_of_range(self):
        with pytest.raises(ValueError, match="axis=5 out of range"):
            RankOneOperator(np.ones(3), np.ones(3), axis=5)

    def test_apply_rejects_shape_mismatch(self):
        R = RankOneOperator(np.ones(3), np.ones(3))
        with pytest.raises(ValueError, match="expected"):
            R.apply(np.ones(5))

    def test_negative_axis_normalised(self):
        """``axis=-1`` on a 3-D op resolves to ``axis=2``."""
        ell = np.ones((2, 3, 4))
        r = np.ones((2, 3, 4))
        R = RankOneOperator(ell, r, axis=-1)
        assert R.axis == 2

    def test_compose_with_identity_via_tensor_product(self):
        """``RankOneOperator & IdentityOperator`` is a valid 2-factor
        TP whose apply is bit-identical to the bare RankOneOperator
        (the §16A.10 separable form for fission's Wave T T.2 lift)."""
        ng, nx, ny = 3, 2, 2
        rng = np.random.default_rng(23)
        ell = rng.uniform(0.1, 1.0, size=(ng, nx, ny))
        r = rng.uniform(0.1, 1.0, size=(ng, nx, ny))
        phi = rng.uniform(0.1, 2.0, size=(ng, nx, ny))

        bare = RankOneOperator(ell, r, axis=0)
        wrapped = RankOneOperator(ell, r, axis=0) & IdentityOperator()

        assert isinstance(wrapped, TensorProductOperator)
        assert wrapped.capabilities == frozenset({CAP_APPLY})
        np.testing.assert_array_equal(wrapped.apply(phi), bare.apply(phi))
