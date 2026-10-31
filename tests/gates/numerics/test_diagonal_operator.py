r"""Tests for :class:`orpheus.numerics.elementary_operators.DiagonalOperator`.

The operator multiplies an operand of one fixed domain ``shape`` by a
coefficient occupying the positions ``spans`` of that shape, constant over
the complementary positions (its ``broadcast_axes``):
``np.expand_dims(coeff, broadcast_axes) * x``. The domain shape is fixed
when the operator is built, and every operand is parsed at that exact shape
(:func:`orpheus.numerics.coordinate_space.coordinates`); a one-position
weight is the placement ``spans=(k,)``.

Verifies the operator's invariants:

* Self-adjointness: ``apply == apply_transpose`` for real coefficients.
* One-position placement: apply along position 0 and non-zero positions.
* Placement oracle: the leading-position ``σ_t`` placement
  ``DiagonalOperator(sigma, shape=(N, ng, nx, ny), spans=(1, 2, 3))``
  reduces EXACTLY to ``sigma[None] * carrier`` (the bit-identity hinge for
  the transport ``MultiplicationOperator(σ_t)``), and the one-position
  placement is bit-identical to the hand-built ``_reshape(w, axis) * x``.
* ``solve`` round-trip: ``solve(apply(x)) == x`` when no zero entries
  (``DiagonalOperator`` KEEPS its native ``solve`` — carve P4 Design B).
* Invertibility: ``is_invertible`` iff the coefficient is all non-zero
  (the value-dependent arm; the guard raises :class:`NotInvertible`).
* Composition under the operator algebra (``+``, ``-``, ``*``, ``@``).
* Construction from a ``DiscreteMeasure`` via :meth:`from_measure`.
* Refusals: an operand of another shape (:class:`NotInCoordinateSpace`,
  at run time) and a malformed placement (``ValueError``, at construction).
  The construction refusal over every mis-filled placement is
  ``test_coordinate_space.py``, ``TestAnEngineIsMintedWithItsDomainShape``.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.geometry.manifold import UNIT_INTERVAL
from orpheus.numerics.coordinate_space import NotInCoordinateSpace
from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.elementary_operators import DiagonalOperator
from orpheus.numerics.operator import (
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
        D = DiagonalOperator(w, shape=(3,), spans=(0,))
        x = np.array([1.0, 1.0, 1.0])
        np.testing.assert_array_equal(D.apply(x), w)

    def test_axis0_2d_input_broadcasts_other_axes(self):
        w = np.array([2.0, 3.0])
        D = DiagonalOperator(w, shape=(2, 4), spans=(0,))
        x = np.ones((2, 4))
        # Each row scaled by w[i]
        expected = np.array([[2.0] * 4, [3.0] * 4])
        np.testing.assert_array_equal(D.apply(x), expected)

    def test_axis_nonzero(self):
        w = np.array([2.0, 3.0, 5.0])
        D = DiagonalOperator(w, shape=(4, 5, 3), spans=(2,))
        x = np.ones((4, 5, 3))
        # Last axis scaled
        expected = np.broadcast_to(w, (4, 5, 3))
        np.testing.assert_array_equal(D.apply(x), expected)

    def test_an_operand_of_another_shape_is_refused(self):
        w = np.array([1.0, 2.0])
        D = DiagonalOperator(w, shape=(2,), spans=(0,))
        with pytest.raises(NotInCoordinateSpace, match="exact shape"):
            D.apply(np.ones(5))


@pytest.mark.l0
class TestDiagonalSelfAdjoint:
    """L0: real-valued diagonal is self-adjoint — apply == apply_transpose."""

    def test_self_adjoint_axis0(self):
        rng = np.random.default_rng(seed=7)
        w = rng.standard_normal(6)
        D = DiagonalOperator(w, shape=(6,), spans=(0,))
        x = rng.standard_normal(6)
        np.testing.assert_array_equal(D.apply(x), D.apply_transpose(x))

    def test_self_adjoint_multiaxis(self):
        rng = np.random.default_rng(seed=11)
        w = rng.standard_normal(4)
        D = DiagonalOperator(w, shape=(3, 4, 5), spans=(1,))
        x = rng.standard_normal((3, 4, 5))
        np.testing.assert_array_equal(D.apply(x), D.apply_transpose(x))


@pytest.mark.l0
class TestDiagonalInvertibilityAndSolve:
    def test_nonzero_weights_are_invertible_and_adjointable(self):
        D = DiagonalOperator(np.array([1.0, 2.0]), shape=(2,), spans=(0,))
        assert D.is_adjointable
        assert D.is_invertible

    def test_zero_weight_revokes_invertibility(self):
        D = DiagonalOperator(np.array([1.0, 0.0, 2.0]), shape=(3,), spans=(0,))
        assert not D.is_invertible

    def test_solve_round_trip(self):
        rng = np.random.default_rng(seed=42)
        w = rng.uniform(0.1, 5.0, size=8)
        D = DiagonalOperator(w, shape=(8,), spans=(0,))
        x = rng.standard_normal(8)
        np.testing.assert_allclose(D.solve(D.apply(x)), x, rtol=1e-13)

    def test_solve_with_zero_weight_raises(self):
        D = DiagonalOperator(np.array([1.0, 0.0, 2.0]), shape=(3,), spans=(0,))
        with pytest.raises(NotInvertible, match="non-zero coefficient"):
            D.solve(np.ones(3))


@pytest.mark.l0
class TestDiagonalComposition:
    """L0: composes correctly under the operator algebra."""

    def test_sum_of_two_diagonals(self):
        D1 = DiagonalOperator(np.array([1.0, 2.0]), shape=(2,), spans=(0,))
        D2 = DiagonalOperator(np.array([10.0, 20.0]), shape=(2,), spans=(0,))
        S = D1 + D2
        assert isinstance(S, OperatorSum)
        x = np.array([1.0, 1.0])
        np.testing.assert_array_equal(S.apply(x), np.array([11.0, 22.0]))

    def test_product_of_two_diagonals(self):
        D1 = DiagonalOperator(np.array([1.0, 2.0]), shape=(2,), spans=(0,))
        D2 = DiagonalOperator(np.array([3.0, 4.0]), shape=(2,), spans=(0,))
        P = D1 @ D2
        assert isinstance(P, OperatorProduct)
        x = np.array([1.0, 1.0])
        # D1 ∘ D2 = pointwise (1*3, 2*4) = (3, 8)
        np.testing.assert_array_equal(P.apply(x), np.array([3.0, 8.0]))

    def test_scaled_diagonal(self):
        D = DiagonalOperator(np.array([1.0, 2.0]), shape=(2,), spans=(0,))
        S = 2.5 * D
        x = np.array([1.0, 1.0])
        np.testing.assert_allclose(S.apply(x), np.array([2.5, 5.0]))

    def test_diagonal_minus_identity_acts_as_w_minus_one(self):
        D = DiagonalOperator(np.array([3.0, 5.0]), shape=(2,), spans=(0,))
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
        D = DiagonalOperator.from_measure(measure, shape=(3,))
        np.testing.assert_array_equal(D.weights, weights)
        x = np.array([1.0, 1.0, 1.0])
        np.testing.assert_array_equal(D.apply(x), weights)

    def test_from_measure_axis_argument(self):
        nodes = np.array([0.1, 0.5])
        weights = np.array([0.5, 0.5])
        measure = DiscreteMeasure(nodes=nodes, weights=weights, support=UNIT_INTERVAL)
        D = DiagonalOperator.from_measure(measure, shape=(3, 2, 4), axis=1)
        x = np.ones((3, 2, 4))
        out = D.apply(x)
        # axis=1 of size 2 scaled by (0.5, 0.5)
        np.testing.assert_allclose(out, 0.5)


@pytest.mark.l0
class TestDiagonalRejectsBadInput:
    def test_a_coefficient_needs_one_position_per_axis(self):
        """A rank-2 coefficient placed on one position is refused when
        built: ``spans`` names one domain position per coefficient axis."""
        with pytest.raises(ValueError, match="DiagonalOperator"):
            DiagonalOperator(np.ones((3, 3)), shape=(3,), spans=(0,))

    def test_an_operand_shaped_like_the_coefficient_is_refused(self):
        """A ``(ng, nx)`` coefficient placed on the trailing positions of an
        ``(N, ng, nx)`` domain admits only an ``(N, ng, nx)`` operand; an
        operand shaped like the coefficient itself, which numpy would
        broadcast, is refused at the exact-shape parse."""
        D = DiagonalOperator(np.ones((2, 5)), shape=(3, 2, 5), spans=(1, 2))
        with pytest.raises(NotInCoordinateSpace, match="exact shape"):
            D.apply(np.ones((2, 5)))
        # And it works for the matching rank-3 carrier:
        np.testing.assert_array_equal(
            D.apply(np.ones((3, 2, 5))),
            np.ones((3, 2, 5)),
        )

    def test_spans_must_be_distinct(self):
        with pytest.raises(ValueError, match="distinct"):
            DiagonalOperator(np.ones((2, 2)), shape=(2, 3), spans=(0, 0))


@pytest.mark.foundation
class TestDiagonalBroadcastOracle:
    r"""Foundation: the placed engine reduces exactly to the two target
    forms.

    Bit-identity here is **inherited verification**, not coincidence:
    the engine performs the SAME single multiply of the placed coefficient
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
        placed on the trailing positions of a ``(N, ng, *spatial)`` domain
        (constant over the LEADING ordinate position) reduces EXACTLY to
        ``sigma[None] * carrier``.

        The carrier is built with **nx ≠ ny** (``(N, 2, 5, 3)``): a
        transposed / wrong-position placement (vv mode #2 variable-swap)
        would either raise (shape mismatch) or mis-scale, so this
        regime discriminates the axis-ordering bug the equal-nx-ny case
        is blind to.
        """
        rng = np.random.default_rng(seed=2025)
        n_ord, ng, nx, ny = 6, 2, 5, 3
        sigma = rng.uniform(0.1, 3.0, size=(ng, nx, ny))
        psi = rng.standard_normal((n_ord, ng, nx, ny))

        D = DiagonalOperator(sigma, shape=(n_ord, ng, nx, ny), spans=(1, 2, 3))
        # The exact form the promoted MultiplicationOperator
        # must reproduce.
        np.testing.assert_array_equal(D.apply(psi), sigma[None] * psi)

    def test_one_position_placement_equals_hand_built_reshape(self):
        """The one-position placement ``spans=(axis,)`` is bit-identical to
        the hand-built ``_reshape(w, axis) * x`` — one operator minted per
        domain shape, for several positions and ranks."""
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
            D = DiagonalOperator(w, shape=carrier_shape, spans=(axis,))
            np.testing.assert_array_equal(
                D.apply(x), legacy_reshape_apply(w, axis, x)
            )

    def test_non_leading_broadcast_axis(self):
        """The engine is not hardcoded to a leading constant position: a
        ``(N, nx)`` coefficient placed on positions ``(0, 2)`` of an
        ``(N, ng, nx)`` domain, constant over the MIDDLE position, reduces
        to ``coeff[:, None, :] * carrier``."""
        rng = np.random.default_rng(seed=13)
        n_ord, ng, nx = 4, 3, 5
        coeff = rng.uniform(0.1, 2.0, size=(n_ord, nx))
        carrier = rng.standard_normal((n_ord, ng, nx))
        D = DiagonalOperator(coeff, shape=(n_ord, ng, nx), spans=(0, 2))
        np.testing.assert_array_equal(
            D.apply(carrier), coeff[:, None, :] * carrier
        )
