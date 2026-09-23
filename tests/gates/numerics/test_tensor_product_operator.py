r"""Tests for :class:`TensorProductOperator` and
:class:`SumOfTensorProductsOperator`.

The Grand Report v3 §6.3 / §15.1 / §15.2 introduce the tensor-product
algebra :math:`A \otimes B \otimes C` for axis-aligned operators. The
contract verified here:

* Per-axis action: each constituent acts on its tagged axis, broadcasting
  on the rest. Order of application is irrelevant for disjoint axes.
* Predicate intersection: ``(A & B)`` is invertible / adjointable iff every
  factor is (the factor-wise recursion; carve P4 — the frozenset retired),
  and the inverse is ALGEBRA-CLOSED: ``(A ⊗ B)⁻¹ = A⁻¹ ⊗ B⁻¹``, so solving
  is spelled ``.inverse().apply(b)`` (``solve`` retired).
* Adjoint distributivity: :math:`(A \otimes B)^* = A^* \otimes B^*`.
* Per-axis composition: :math:`(A \otimes B) @ (C \otimes D) = (A @ C) \otimes (B \otimes D)`.
* Sum-of-tensor-products: :math:`\Sigma_k A_k \otimes B_k` is the §15.2
  canonical scattering / streaming form.
* ``__and__`` flattening: ``(A & B) & C == A & (B & C)``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.functional import InnerProductFunctional
from orpheus.numerics.operator import (
    DiagonalOperator,
    IdentityOperator,
    IncompatibleOperatorComposition,
    NotInvertible,
    PermutationOperator,
    RankOneOperator,
    SumOfTensorProductsOperator,
    TensorProductOperator,
    ZeroOperator,
    outer,
)
from orpheus.numerics.space import FunctionSpace


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
class TestTensorProductPredicates:
    """L0: predicate intersection (invertible/adjointable iff every factor)."""

    def test_intersection_adjointable_not_invertible(self):
        T = ZeroOperator() & IdentityOperator()
        # Zero is adjointable (self-adjoint) but not invertible; Identity is
        # both — the factor-wise intersection keeps adjointable, drops
        # invertible.
        assert T.is_adjointable
        assert not T.is_invertible

    def test_full_intersection(self):
        T = IdentityOperator() & IdentityOperator()
        assert T.is_invertible
        assert T.is_adjointable

    def test_invertibility_propagates_through_diagonals(self):
        """Rewired at carve P4: ``TensorProductOperator.solve`` retired
        (algebra-closed inverse ``(A⊗B)⁻¹ = A⁻¹⊗B⁻¹``); solving is
        ``.inverse().apply`` — same round-trip value, same tolerance."""
        D1 = DiagonalOperator(np.array([1.0, 2.0]), axis=0)
        D2 = DiagonalOperator(np.array([3.0, 4.0]), axis=1)
        T = D1 & D2
        assert T.is_invertible
        assert not hasattr(T, "solve")  # retired at carve P4
        # Inverse round-trip (the factor-wise inverse applies each factor's
        # division in stored order — bit-identical to the retired solve).
        x = np.array([[1.0, 1.0], [1.0, 1.0]])
        np.testing.assert_allclose(T.inverse().apply(T.apply(x)), x, rtol=1e-13)

    def test_inverse_revoked_by_zero_weight_factor(self):
        D_zero = DiagonalOperator(np.array([1.0, 0.0]), axis=0)
        D_ok = DiagonalOperator(np.array([2.0, 3.0]), axis=1)
        T = D_zero & D_ok
        assert not T.is_invertible
        with pytest.raises(NotInvertible, match="every factor"):
            T.inverse()


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
    """Constituent must expose ``apply`` (the eager composition guard)."""

    def test_construction_rejects_non_apply_factor(self):
        # A fake operator that genuinely LACKS apply (post carve P4 the
        # ctor guard probes the verb itself — ``callable(op.apply)`` —
        # not a capability advertisement).
        class FakeOp:
            pass

        with pytest.raises(TypeError, match="apply"):
            TensorProductOperator((FakeOp(),))


# ─────────────────────────────────────────────────────────────────────
# The spaces (G6.3 step 8.0) — resolved by AGREEMENT, never by position
# ─────────────────────────────────────────────────────────────────────


_V = FunctionSpace(name="V", shape=(3,))
_W = FunctionSpace(name="W", shape=(3,))


def _bound_swap(domain, codomain):
    """A leaf that genuinely maps ``domain -> codomain`` (a 3-cycle on axis 0)."""
    return PermutationOperator(
        np.array([1, 2, 0]), axis=0, domain=domain, codomain=codomain,
    )


@pytest.mark.l0
@pytest.mark.catches("ERR-076")
class TestTensorProductSpaces:
    r"""L0: :math:`A \otimes I` is bound exactly where :math:`A` is.

    Before step 8.0 the product derived nothing, so a binding real at the
    inner factor was invisible at the composite — and since
    ``AdjointOperator`` reads the spaces to apply the metrics, ``.H`` also
    degraded silently to the Euclidean transpose.
    """

    def test_the_lone_bound_factor_supplies_both_ends(self):
        r"""``is``-identity, and it names WHICH end is which.

        A ``domain``/``codomain`` swap changes no arithmetic and no shape —
        both ends are ``(3,)`` here, mirroring the real case where
        :math:`|\Gamma_+| = |\Gamma_-|` on every reachable face — so an
        identity assertion naming the ends is the only instrument that can
        see it.
        """
        A = _bound_swap(_V, _W)
        T = A & IdentityOperator()
        assert T.domain is _V
        assert T.codomain is _W

    def test_the_derived_spaces_do_not_depend_on_factor_ORDER(self):
        r"""⭐ The law that forced agreement over position.

        The factors commute by contract — the class docstring's "the order
        does not matter" — so a position-based rule (``ops[0].domain`` /
        ``ops[-1].codomain``, which is what :class:`OperatorProduct` correctly
        uses) would give an order-DEPENDENT answer for an order-INDEPENDENT
        operator. Here it would have handed ``I & A`` a ``None`` domain.
        """
        A = _bound_swap(_V, _W)
        forward, reversed_ = A & IdentityOperator(), IdentityOperator() & A
        assert forward.domain is reversed_.domain is _V
        assert forward.codomain is reversed_.codomain is _W

    def test_silent_factors_leave_the_product_unbound(self):
        """All-silent stays ``None`` — the transitional contract (#330).

        Periodic realizes to ``I & I`` and the fission kernel to an unbound
        rank-one dyad; neither declares a space, and neither may be given one
        by inference.
        """
        assert (IdentityOperator() & IdentityOperator()).domain is None
        assert (IdentityOperator() & IdentityOperator()).codomain is None
        unbound = outer(
            np.arange(3.0), InnerProductFunctional(np.arange(3.0), axis=0)
        ) & IdentityOperator()
        assert unbound.domain is None and unbound.codomain is None

    def test_factors_that_AGREE_compose(self):
        """The positive leg: agreement is not merely "no disagreement"."""
        T = _bound_swap(_V, _W) & _bound_swap(_V, _W)
        assert T.domain is _V and T.codomain is _W

    @pytest.mark.parametrize(
        "second, role",
        [
            (lambda: _bound_swap(_W, _W), "domain"),
            (lambda: _bound_swap(_V, _V), "codomain"),
        ],
    )
    def test_factors_that_DISAGREE_are_refused(self, second, role):
        """Whole-space bindings that differ cannot both be the product's.

        The refusal is the honest answer, not a fallback: picking one would
        make the composite advertise a space it does not have. The message
        names the machinery a legitimate per-axis case would need.

        ⚠ The match carries the OWNER and the plural — ``codomain`` contains
        ``domain``, so a bare role would let the domain row pass on a codomain
        refusal and the parametrization would gate one branch twice; and the
        law is shared with :class:`OperatorSum`, so the owner is what says
        which composite refused.
        """
        with pytest.raises(
            IncompatibleOperatorComposition,
            match=f"TensorProductOperator requires equal {role}s",
        ):
            _ = _bound_swap(_V, _W) & second()

    def test_the_inverse_INVERTS_the_binding(self):
        r""":math:`(A \otimes I)^{-1} : W \to V` — ends swapped, not dropped."""
        inv = (_bound_swap(_V, _W) & IdentityOperator()).inverse()
        assert inv.domain is _W
        assert inv.codomain is _V

    def test_the_adjoint_SWAPS_the_binding(self):
        r""":math:`(A \otimes I)^{*} : W \to V`, via the ``AdjointOperator``."""
        adj = (_bound_swap(_V, _W) & IdentityOperator()).adjoint()
        assert adj.domain is _W
        assert adj.codomain is _V


@pytest.mark.l0
class TestSumOfTensorProductsSpaces:
    """L0: a sum of tensor products is a SUM — same agreement law."""

    def test_the_summands_supply_the_spaces(self):
        S = SumOfTensorProductsOperator(
            (_bound_swap(_V, _W) & IdentityOperator(),
             _bound_swap(_V, _W) & IdentityOperator())
        )
        assert S.domain is _V and S.codomain is _W

    def test_disagreeing_summands_are_refused(self):
        with pytest.raises(
            IncompatibleOperatorComposition,
            match="SumOfTensorProductsOperator requires equal codomains",
        ):
            SumOfTensorProductsOperator(
                (_bound_swap(_V, _W) & IdentityOperator(),
                 _bound_swap(_V, _V) & IdentityOperator())
            )


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

    def test_not_invertible_no_solve(self):
        """A sum of tensor products is not invertible (sum-of-inverses !=
        inverse-of-sum; ``is_invertible`` inherits the base ``False``), and
        the ``solve`` verb does not exist anywhere on the type (carve P4)."""
        D_a = DiagonalOperator(np.array([1.0]))
        D_b = DiagonalOperator(np.array([1.0]))
        S = SumOfTensorProductsOperator((D_a & D_b,))
        assert not S.is_invertible
        assert not hasattr(S, "solve")  # never existed post carve P4
        assert not hasattr(S, "inverse")  # structural refusal is static

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
        R = RankOneOperator(ell, InnerProductFunctional(r, axis=0))
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

        R = RankOneOperator(chi, InnerProductFunctional(nu_sigma_f, axis=0))
        out = R.apply(phi)

        # Hand-computed per-cell rank-1 action.
        inner = np.einsum("gxy,gxy->xy", nu_sigma_f, phi)
        expected = chi * inner[None, :, :]
        np.testing.assert_array_almost_equal_nulp(out, expected, nulp=4)

    def test_structurally_non_invertible_with_working_transpose(self):
        """Rank-1 ops are non-invertible by construction (no ``solve``, no
        ``inverse()`` — structural refusal is static, Design C), but the dyad
        HAS a transpose: the dual dyad ``|w⟩⟨v|`` (campaign #276), available
        when the row is an ``InnerProductFunctional`` (the predicate pin for
        that arm lives in ``test_outer_dyad.py::TestPredicates``)."""
        R = RankOneOperator(np.ones(3), InnerProductFunctional(np.ones(3)))
        assert not R.is_invertible
        assert not hasattr(R, "solve")
        assert not hasattr(R, "inverse")
        # The transpose VERB works for the IPF row (the dual dyad).
        np.testing.assert_array_equal(
            R.apply_transpose(np.full(3, 2.0)), np.full(3, 6.0)
        )

    def test_apply_rejects_shape_mismatch(self):
        """A carrier that doesn't align with the row co-vector raises.

        Post-refactor the contraction is owned by the functional, so a
        mismatched carrier surfaces as numpy's broadcast error inside
        ``InnerProductFunctional.evaluate`` (``(weight * x).sum(...)``) —
        the behaviour (a mismatched apply raises) is preserved; the error
        is now numpy's, since the functional owns the contraction.
        """
        R = RankOneOperator(np.ones(3), InnerProductFunctional(np.ones(3)))
        with pytest.raises(ValueError, match="broadcast"):
            R.apply(np.ones(5))

    def test_apply_negative_axis_contracts_last_axis(self):
        """``axis=-1`` on a 3-D row co-vector contracts the trailing axis.

        Post-refactor there is no ``RankOneOperator.axis`` attribute — the
        contraction axis lives on the functional, and numpy resolves
        ``axis=-1`` natively in ``InnerProductFunctional.evaluate``. This
        behavioural row pins that ``apply`` contracts the LAST axis (the
        old test asserted the normalised ``R.axis == 2``; the equivalent
        post-refactor claim is the trailing-axis action).
        """
        ell = np.ones((2, 3, 4))
        r = np.arange(2 * 3 * 4, dtype=float).reshape(2, 3, 4)
        x = np.ones((2, 3, 4))
        R = RankOneOperator(ell, InnerProductFunctional(r, axis=-1))
        # ⟨r, x⟩ over the last axis, keepdims → (2, 3, 1); broadcast by ell.
        inner = (r * x).sum(axis=-1, keepdims=True)
        np.testing.assert_array_equal(R.apply(x), ell * inner)

    def test_compose_with_identity_via_tensor_product(self):
        """``RankOneOperator & IdentityOperator`` is a valid 2-factor
        TP whose apply is bit-identical to the bare RankOneOperator
        (the §16A.10 separable form for fission's Wave T T.2 lift)."""
        ng, nx, ny = 3, 2, 2
        rng = np.random.default_rng(23)
        ell = rng.uniform(0.1, 1.0, size=(ng, nx, ny))
        r = rng.uniform(0.1, 1.0, size=(ng, nx, ny))
        phi = rng.uniform(0.1, 2.0, size=(ng, nx, ny))

        bare = outer(ell, InnerProductFunctional(r, axis=0))
        wrapped = outer(ell, InnerProductFunctional(r, axis=0)) & IdentityOperator()

        assert isinstance(wrapped, TensorProductOperator)
        # The dyad factor is non-invertible → the factor-wise TP recursion
        # drops invertibility (its adjointable arm is pinned below).
        assert not wrapped.is_invertible
        np.testing.assert_array_equal(wrapped.apply(phi), bare.apply(phi))

    def test_dyad_identity_tensor_product_is_adjointable(self):
        """``(|v⟩⟨w| & I).is_adjointable`` — the TP recursion over an
        IPF-rowed dyad (working dual-dyad transpose) and Identity must
        advertise the transpose, as the pre-carve caps intersection did."""
        wrapped = outer(
            np.ones((2, 2)), InnerProductFunctional(np.ones((2, 2)), axis=0)
        ) & IdentityOperator()
        assert wrapped.is_adjointable
