r"""Tests for :class:`orpheus.numerics.operator.DiagonalOperator`.

Verifies the operator's invariants:

* Self-adjointness: ``apply == apply_transpose`` for real weights.
* Apply along axis 0 and along non-zero axes (broadcasting).
* ``solve`` round-trip: ``solve(apply(x)) == x`` when no zero weights.
* Capability set: ``CAP_SOLVE`` advertised iff weights are all non-zero.
* Composition under the operator algebra (``+``, ``-``, ``*``, ``@``).
* Construction from a ``DiscreteMeasure`` via :meth:`from_measure`.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    CAP_SOLVE,
    DiagonalOperator,
    IdentityOperator,
    MissingCapability,
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
class TestDiagonalCapabilitiesAndSolve:
    def test_nonzero_weights_advertise_solve(self):
        D = DiagonalOperator(np.array([1.0, 2.0]))
        assert CAP_APPLY in D.capabilities
        assert CAP_APPLY_TRANSPOSE in D.capabilities
        assert CAP_SOLVE in D.capabilities

    def test_zero_weight_revokes_solve(self):
        D = DiagonalOperator(np.array([1.0, 0.0, 2.0]))
        assert CAP_APPLY in D.capabilities
        assert CAP_SOLVE not in D.capabilities

    def test_solve_round_trip(self):
        rng = np.random.default_rng(seed=42)
        w = rng.uniform(0.1, 5.0, size=8)
        D = DiagonalOperator(w)
        x = rng.standard_normal(8)
        np.testing.assert_allclose(D.solve(D.apply(x)), x, rtol=1e-13)

    def test_solve_with_zero_weight_raises(self):
        D = DiagonalOperator(np.array([1.0, 0.0, 2.0]))
        with pytest.raises(MissingCapability, match="non-zero weights"):
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
        measure = DiscreteMeasure(nodes=nodes, weights=weights, space="[0,1]")
        D = DiagonalOperator.from_measure(measure)
        np.testing.assert_array_equal(D.weights, weights)
        x = np.array([1.0, 1.0, 1.0])
        np.testing.assert_array_equal(D.apply(x), weights)

    def test_from_measure_axis_argument(self):
        nodes = np.array([0.1, 0.5])
        weights = np.array([0.5, 0.5])
        measure = DiscreteMeasure(nodes=nodes, weights=weights, space="[0,1]")
        D = DiagonalOperator.from_measure(measure, axis=1)
        x = np.ones((3, 2, 4))
        out = D.apply(x)
        # axis=1 of size 2 scaled by (0.5, 0.5)
        np.testing.assert_allclose(out, 0.5)


@pytest.mark.l0
class TestDiagonalRejectsBadInput:
    def test_2d_weights_rejected(self):
        with pytest.raises(ValueError, match="1-D"):
            DiagonalOperator(np.ones((3, 3)))
