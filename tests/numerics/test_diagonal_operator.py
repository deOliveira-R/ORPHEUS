r"""Tests for :class:`orpheus.numerics.operator.DiagonalOperator`.

The operator is the N-D broadcast engine: it multiplies a carrier
tensor by a coefficient occupying a sub-product of the carrier's axes,
broadcasting over the complementary ``broadcast_axes`` (i.e.
``np.expand_dims(coeff, broadcast_axes) * carrier``). The 1-D
``(weights, axis)`` form is the rank-1 special case (one axis,
rank-agnostic broadcast over the rest) and stays bit-identical.

Verifies the operator's invariants:

* Self-adjointness: ``apply == apply_transpose`` for real coefficients.
* 1-D special case: apply along axis 0 and non-zero axes (broadcasting).
* N-D broadcast oracle: the leading-axis ``σ_t`` form
  ``DiagonalOperator(sigma, broadcast_axes=(0,))`` reduces EXACTLY to
  ``sigma[None] * carrier`` (the bit-identity hinge for the transport
  ``MultiplicationOperator(σ_t)`` promotion), and the 1-D form is
  bit-identical to the old ``_reshape(w, axis) * x``.
* ``solve`` round-trip: ``solve(apply(x)) == x`` when no zero entries
  (``DiagonalOperator`` KEEPS its native ``solve`` — carve P4 Design B).
* Invertibility: ``is_invertible`` iff the coefficient is all non-zero
  (the value-dependent arm; the guard raises :class:`NotInvertible`).
* Composition under the operator algebra (``+``, ``-``, ``*``, ``@``).
* Construction from a ``DiscreteMeasure`` via :meth:`from_measure`.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.manifold import UNIT_INTERVAL
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.operator import (
    DiagonalOperator,
    IdentityOperator,
    NotInvertible,
    OperatorProduct,
    OperatorSum,
)


@pytest.mark.l0
class TestDiagonalApplyShape:
    """L0: shape and basic numerical correctness of `apply`."""

    def test_axis0_1d_input(self):
        w = np.array([2.0, 3.0, 5.0])
        D = DiagonalOperator(w)
        x = np.array([1.0, 1.0, 1.0])
        np.testing.assert_array_equal(D.apply(x), w)

    def test_axis0_2d_input_broadcasts_other_axes(self):
        w = np.array([2.0, 3.0])
        D = DiagonalOperator(w, axis=0)
        x = np.ones((2, 4))
        # Each row scaled by w[i]
        expected = np.array([[2.0] * 4, [3.0] * 4])
        np.testing.assert_array_equal(D.apply(x), expected)

    def test_axis_nonzero(self):
        w = np.array([2.0, 3.0, 5.0])
        D = DiagonalOperator(w, axis=2)
        x = np.ones((4, 5, 3))
        # Last axis scaled
        expected = np.broadcast_to(w, (4, 5, 3))
        np.testing.assert_array_equal(D.apply(x), expected)

    def test_axis_size_mismatch_raises(self):
        w = np.array([1.0, 2.0])
        D = DiagonalOperator(w, axis=0)
        with pytest.raises(ValueError, match="axis size"):
            D.apply(np.ones(5))


@pytest.mark.l0
class TestDiagonalSelfAdjoint:
    """L0: real-valued diagonal is self-adjoint — apply == apply_transpose."""

    def test_self_adjoint_axis0(self):
        rng = np.random.default_rng(seed=7)
        w = rng.standard_normal(6)
        D = DiagonalOperator(w)
        x = rng.standard_normal(6)
        np.testing.assert_array_equal(D.apply(x), D.apply_transpose(x))

    def test_self_adjoint_multiaxis(self):
        rng = np.random.default_rng(seed=11)
        w = rng.standard_normal(4)
        D = DiagonalOperator(w, axis=1)
        x = rng.standard_normal((3, 4, 5))
        np.testing.assert_array_equal(D.apply(x), D.apply_transpose(x))


@pytest.mark.l0
class TestDiagonalInvertibilityAndSolve:
    def test_nonzero_weights_are_invertible_and_adjointable(self):
        D = DiagonalOperator(np.array([1.0, 2.0]))
        assert D.is_adjointable
        assert D.is_invertible

    def test_zero_weight_revokes_invertibility(self):
        D = DiagonalOperator(np.array([1.0, 0.0, 2.0]))
        assert not D.is_invertible

    def test_solve_round_trip(self):
        rng = np.random.default_rng(seed=42)
        w = rng.uniform(0.1, 5.0, size=8)
        D = DiagonalOperator(w)
        x = rng.standard_normal(8)
        np.testing.assert_allclose(D.solve(D.apply(x)), x, rtol=1e-13)

    def test_solve_with_zero_weight_raises(self):
        D = DiagonalOperator(np.array([1.0, 0.0, 2.0]))
        with pytest.raises(NotInvertible, match="non-zero coefficient"):
            D.solve(np.ones(3))


@pytest.mark.l0
class TestDiagonalComposition:
    """L0: composes correctly under the operator algebra."""

    def test_sum_of_two_diagonals(self):
        D1 = DiagonalOperator(np.array([1.0, 2.0]))
        D2 = DiagonalOperator(np.array([10.0, 20.0]))
        S = D1 + D2
        assert isinstance(S, OperatorSum)
        x = np.array([1.0, 1.0])
        np.testing.assert_array_equal(S.apply(x), np.array([11.0, 22.0]))

    def test_product_of_two_diagonals(self):
        D1 = DiagonalOperator(np.array([1.0, 2.0]))
        D2 = DiagonalOperator(np.array([3.0, 4.0]))
        P = D1 @ D2
        assert isinstance(P, OperatorProduct)
        x = np.array([1.0, 1.0])
        # D1 ∘ D2 = pointwise (1*3, 2*4) = (3, 8)
        np.testing.assert_array_equal(P.apply(x), np.array([3.0, 8.0]))

    def test_scaled_diagonal(self):
        D = DiagonalOperator(np.array([1.0, 2.0]))
        S = 2.5 * D
        x = np.array([1.0, 1.0])
        np.testing.assert_allclose(S.apply(x), np.array([2.5, 5.0]))

    def test_diagonal_minus_identity_acts_as_w_minus_one(self):
        D = DiagonalOperator(np.array([3.0, 5.0]))
        I = IdentityOperator()
        diff = D - I
        x = np.array([1.0, 1.0])
        np.testing.assert_array_equal(diff.apply(x), np.array([2.0, 4.0]))


@pytest.mark.l0
class TestDiagonalFromMeasure:
    def test_from_discrete_measure(self):
        nodes = np.array([0.1, 0.5, 0.9])
        weights = np.array([0.2, 0.4, 0.4])
        measure = DiscreteMeasure(nodes=nodes, weights=weights, support=UNIT_INTERVAL)
        D = DiagonalOperator.from_measure(measure)
        np.testing.assert_array_equal(D.weights, weights)
        x = np.array([1.0, 1.0, 1.0])
        np.testing.assert_array_equal(D.apply(x), weights)

    def test_from_measure_axis_argument(self):
        nodes = np.array([0.1, 0.5])
        weights = np.array([0.5, 0.5])
        measure = DiscreteMeasure(nodes=nodes, weights=weights, support=UNIT_INTERVAL)
        D = DiagonalOperator.from_measure(measure, axis=1)
        x = np.ones((3, 2, 4))
        out = D.apply(x)
        # axis=1 of size 2 scaled by (0.5, 0.5)
        np.testing.assert_allclose(out, 0.5)


@pytest.mark.l0
class TestDiagonalRejectsBadInput:
    def test_2d_coefficient_without_broadcast_axes_rejected(self):
        """The 1-D special case (no ``broadcast_axes``) requires a 1-D
        coefficient. An N-D coefficient must declare the carrier axes it
        broadcasts over."""
        with pytest.raises(ValueError, match="1-D special case"):
            DiagonalOperator(np.ones((3, 3)))

    def test_coefficient_rank_must_match_complement_rank(self):
        """NEW contract: with ``broadcast_axes`` the coefficient rank
        MUST equal ``carrier.ndim - len(broadcast_axes)`` (NOT "1-D").
        A (ng, nx) coefficient broadcasting over the leading ordinate
        axis demands a 3-D ``(N, ng, nx)`` carrier; a rank-2 carrier is
        rejected (expected coefficient rank 1 ≠ actual rank 2)."""
        D = DiagonalOperator(np.ones((2, 5)), broadcast_axes=(0,))
        with pytest.raises(ValueError, match="rank-1 coefficient"):
            D.apply(np.ones((2, 5)))  # carrier rank 2 ⇒ expects rank-1 coeff
        # And it works for the matching rank-3 carrier:
        np.testing.assert_array_equal(
            D.apply(np.ones((3, 2, 5))),
            np.ones((3, 2, 5)),
        )

    def test_broadcast_axes_must_be_distinct(self):
        with pytest.raises(ValueError, match="distinct"):
            DiagonalOperator(np.ones((2, 3)), broadcast_axes=(0, 0))


@pytest.mark.foundation
class TestDiagonalBroadcastOracle:
    r"""Foundation: the N-D broadcast engine reduces exactly to the two
    target forms.

    Bit-identity here is **inherited verification**, not coincidence:
    the engine performs the SAME single broadcast-multiply
    (``np.expand_dims(coeff, broadcast_axes) * carrier``) as the
    reference — reduction depth 1, no summation — so
    :func:`np.testing.assert_array_equal` (0 ULP) is the correct gate
    and IS expected. A non-bit-identical result would only be
    acceptable if it arose from a MORE principled construction (the
    user's principled-equivalence directive); a pure broadcast can
    produce no such re-association, so a reviewer MUST NOT silently
    loosen this to ``allclose``/``nulp``.
    """

    def test_leading_axis_equals_sigma_None_times_carrier(self):
        r"""The σ_t case S3b needs: a ``(ng, *spatial)`` coefficient
        broadcast over the LEADING ordinate axis of a ``(N, ng, *spatial)``
        carrier reduces EXACTLY to ``sigma[None] * carrier``.

        The carrier is built with **nx ≠ ny** (``(N, 2, 5, 3)``): a
        transposed / wrong-axis broadcast (vv mode #2 variable-swap)
        would either raise (rank/shape mismatch) or mis-scale, so this
        regime discriminates the axis-ordering bug the equal-nx-ny case
        is blind to.
        """
        rng = np.random.default_rng(seed=2025)
        n_ord, ng, nx, ny = 6, 2, 5, 3
        sigma = rng.uniform(0.1, 3.0, size=(ng, nx, ny))
        psi = rng.standard_normal((n_ord, ng, nx, ny))

        D = DiagonalOperator(sigma, broadcast_axes=(0,))
        # The exact form the promoted MultiplicationOperator
        # must reproduce.
        np.testing.assert_array_equal(D.apply(psi), sigma[None] * psi)

    def test_one_d_axis_mode_equals_legacy_reshape(self):
        """The 1-D special case is bit-identical to the old
        ``_reshape(w, axis) * x`` broadcast — verified against an
        independent hand-built reshape for several axes/ranks."""
        rng = np.random.default_rng(seed=7)

        def legacy_reshape_apply(w, axis, x):
            shape = [1] * x.ndim
            shape[axis] = -1
            return w.reshape(shape) * x

        for axis, carrier_shape in [
            (0, (4,)),
            (0, (4, 7)),
            (1, (3, 4, 5)),
            (2, (3, 4, 5)),
        ]:
            n = carrier_shape[axis]
            w = rng.standard_normal(n)
            x = rng.standard_normal(carrier_shape)
            D = DiagonalOperator(w, axis=axis)
            np.testing.assert_array_equal(
                D.apply(x), legacy_reshape_apply(w, axis, x)
            )

    def test_non_leading_broadcast_axis(self):
        """The engine is not hardcoded to ``axis=0``: a ``(N, nx)``
        coefficient broadcasting over a MIDDLE axis ``(N, ng, nx)``
        reduces to ``coeff[:, None, :] * carrier``."""
        rng = np.random.default_rng(seed=13)
        n_ord, ng, nx = 4, 3, 5
        coeff = rng.uniform(0.1, 2.0, size=(n_ord, nx))
        carrier = rng.standard_normal((n_ord, ng, nx))
        D = DiagonalOperator(coeff, broadcast_axes=(1,))
        np.testing.assert_array_equal(
            D.apply(carrier), coeff[:, None, :] * carrier
        )
