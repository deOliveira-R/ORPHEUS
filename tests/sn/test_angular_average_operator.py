"""Tests for AngularAverageOperator (Wave 1 / C1.1).

The operator is the §15.2 geometric projection G_diff in the white
BC tensor decomposition R_white = G_diff ⊗ α. It implements the
cosine-weighted Lambertian average over an outgoing hemisphere; the
incoming flux at a white surface equals α × (this average).
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.operator import CAP_APPLY
from orpheus.sn.angular_operator import AngularAverageOperator
from orpheus.sn.quadrature import (
    GaussLegendre1D,
    LebedevSphere,
    LevelSymmetricSN,
)


# ─────────────────────────────────────────────────────────────────────
# Outgoing-hemisphere mask correctness (L0, parametrized over 3 quadratures)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestOutgoingHemisphereMask:
    """The cos_w array zeroes contributions from non-outgoing ordinates."""

    def test_lebedev_xmax_outgoing_only(self):
        quad = LebedevSphere.create(17)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        # Apply to a flux that is +1 on outgoing (mu_x > 0) and -1 on incoming.
        # Cosine-weighted average must equal +1 (only outgoing contributes).
        psi = np.where(quad.mu_x > 0, +1.0, -1.0)
        result = op.apply(psi)
        # All ordinates carry the same broadcast value.
        np.testing.assert_allclose(result, +1.0, rtol=1e-14)

    def test_level_symmetric_ymin_outgoing_only(self):
        quad = LevelSymmetricSN.create(4)
        op = AngularAverageOperator.from_quadrature(quad, axis="y", outward_sign=-1)
        # ymin face: outgoing = mu_y < 0.
        psi = np.where(quad.mu_y < 0, +1.0, -1.0)
        result = op.apply(psi)
        np.testing.assert_allclose(result, +1.0, rtol=1e-14)

    def test_gauss_legendre_1d_xmax(self):
        quad = GaussLegendre1D.create(8)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        psi = np.where(quad.mu_x > 0, +1.0, -1.0)
        result = op.apply(psi)
        np.testing.assert_allclose(result, +1.0, rtol=1e-14)


# ─────────────────────────────────────────────────────────────────────
# Conservation (L1)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l1
class TestCosineWeightedCurrentConservation:
    r"""The outgoing current equals the (cosine-weighted) average × normalization.

    By construction, the operator returns ``ψ_avg = (Σ w·|μ|·ψ) / (Σ w·|μ|)``
    broadcast over all ordinates. The outgoing cosine-weighted
    current of the OUTPUT (with the SAME outgoing-cos weights) is
    therefore ``ψ_avg × norm = Σ w·|μ|·ψ`` — exactly the outgoing
    cosine-weighted current of the INPUT. This is the white BC
    conservation property: outgoing current → incoming current with
    no loss (under α=1).
    """

    def test_outgoing_current_equals_incoming_lebedev(self):
        quad = LebedevSphere.create(17)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        rng = np.random.default_rng(42)
        psi = rng.uniform(0.5, 1.5, size=quad.N)  # positive flux
        result = op.apply(psi)
        outgoing_mask = quad.mu_x > 0
        cos_w_out = quad.weights * quad.mu_x  # outward_sign=+1
        cos_w_out = np.where(outgoing_mask, cos_w_out, 0.0)
        j_in = (cos_w_out * psi).sum()
        # The result is constant in ordinate (Lambertian); compute
        # the outgoing cosine-weighted current of the broadcast result.
        j_out = (cos_w_out * result).sum()
        np.testing.assert_allclose(j_in, j_out, rtol=1e-13)

    def test_outgoing_current_with_anisotropic_input(self):
        quad = LevelSymmetricSN.create(6)
        op = AngularAverageOperator.from_quadrature(quad, axis="z", outward_sign=+1)
        # Strongly anisotropic input — exercises the cosine weighting.
        psi = np.exp(quad.mu_z) * (quad.mu_z + 2.0)
        result = op.apply(psi)
        outgoing_mask = quad.mu_z > 0
        cos_w = np.where(outgoing_mask, quad.weights * quad.mu_z, 0.0)
        j_in = (cos_w * psi).sum()
        j_out = (cos_w * result).sum()
        np.testing.assert_allclose(j_in, j_out, rtol=1e-13)


# ─────────────────────────────────────────────────────────────────────
# Self-adjointness under cosine-weighted inner product (L1)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l1
class TestSelfAdjointnessCosineWeighted:
    r"""Under the cosine-weighted inner product over the outgoing
    hemisphere ⟨x, y⟩_w = Σ w·|μ|·x_n·y_n, the operator is
    self-adjoint: ⟨A x, y⟩_w == ⟨x, A y⟩_w.

    Mechanically: (A x)_n = (Σ_m w_m μ_m x_m) / norm = c (scalar).
    Then ⟨A x, y⟩_w = c · (Σ w μ y) and ⟨x, A y⟩_w = (Σ w μ x) · c' where
    c' = (Σ w μ y) / norm. Both equal (Σ w μ x)(Σ w μ y) / norm.

    Tolerance ULP-bound × (reduction depth N).
    """

    def test_self_adjoint_lebedev(self):
        quad = LebedevSphere.create(17)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        rng = np.random.default_rng(0)
        x = rng.uniform(0, 1, size=quad.N)
        y = rng.uniform(0, 1, size=quad.N)
        outgoing = quad.mu_x > 0
        cos_w = np.where(outgoing, quad.weights * quad.mu_x, 0.0)
        Ax = op.apply(x)
        Ay = op.apply(y)
        lhs = (cos_w * Ax * y).sum()
        rhs = (cos_w * x * Ay).sum()
        np.testing.assert_allclose(lhs, rhs, rtol=1e-13)

    def test_self_adjoint_level_symmetric(self):
        quad = LevelSymmetricSN.create(4)
        op = AngularAverageOperator.from_quadrature(quad, axis="y", outward_sign=-1)
        rng = np.random.default_rng(1)
        x = rng.uniform(0, 1, size=quad.N)
        y = rng.uniform(0, 1, size=quad.N)
        outgoing = (-1.0 * quad.mu_y) > 0
        cos_w = np.where(outgoing, quad.weights * (-1.0 * quad.mu_y), 0.0)
        Ax = op.apply(x)
        Ay = op.apply(y)
        lhs = (cos_w * Ax * y).sum()
        rhs = (cos_w * x * Ay).sum()
        np.testing.assert_allclose(lhs, rhs, rtol=1e-13)


# ─────────────────────────────────────────────────────────────────────
# Capability set (L0)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestCapabilitySet:
    def test_caps(self):
        quad = LebedevSphere.create(5)
        op = AngularAverageOperator.from_quadrature(quad, axis="z", outward_sign=+1)
        assert op.capabilities == frozenset({CAP_APPLY})


# ─────────────────────────────────────────────────────────────────────
# Axis selection — z on Lebedev OK, z on GL raises (L0)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestAxisDispatch:
    def test_z_on_lebedev_succeeds(self):
        quad = LebedevSphere.create(11)
        op = AngularAverageOperator.from_quadrature(quad, axis="z", outward_sign=+1)
        assert op.n_ordinates == quad.N

    def test_z_on_gauss_legendre_raises(self):
        quad = GaussLegendre1D.create(8)
        with pytest.raises(ValueError, match="mu_z"):
            AngularAverageOperator.from_quadrature(quad, axis="z", outward_sign=+1)

    def test_invalid_axis_raises(self):
        quad = LebedevSphere.create(5)
        with pytest.raises(ValueError, match="Unknown axis"):
            AngularAverageOperator.from_quadrature(quad, axis="w", outward_sign=+1)

    def test_invalid_outward_sign_raises(self):
        quad = LebedevSphere.create(5)
        with pytest.raises(ValueError, match="outward_sign"):
            AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=0)


# ─────────────────────────────────────────────────────────────────────
# Defensive copy — operator is decoupled from the source quadrature (L0)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestDefensiveCopy:
    def test_quadrature_reference_not_held(self):
        quad = LebedevSphere.create(5)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        # No quadrature attribute on the operator (no held reference).
        assert not hasattr(op, "quadrature")
        assert not hasattr(op, "_quadrature")
        # The internal cos_w is its own copy.
        assert op._cos_w is not quad.weights
        assert op._cos_w is not quad.mu_x

    def test_output_is_fresh_array(self):
        """Calling code may mutate the output without affecting the
        operator's internal state."""
        quad = LebedevSphere.create(5)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        psi = np.ones(quad.N)
        result1 = op.apply(psi)
        result1[0] = 999.0
        result2 = op.apply(psi)
        # result2 must NOT see the mutation (operator state unaffected).
        assert result2[0] != 999.0


# ─────────────────────────────────────────────────────────────────────
# Apply input-shape validation (L0)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestInputShape:
    def test_wrong_n_ordinates_raises(self):
        quad = LebedevSphere.create(5)  # N=14
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        with pytest.raises(ValueError, match="psi.shape"):
            op.apply(np.ones(13))

    def test_multi_dim_input_broadcasts(self):
        """Shape (N_ord, 5, 3) input — average is taken over leading axis only,
        broadcast back."""
        quad = LebedevSphere.create(11)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        rng = np.random.default_rng(7)
        psi = rng.uniform(0, 1, size=(quad.N, 5, 3))
        result = op.apply(psi)
        assert result.shape == psi.shape
        # All leading-axis slices identical.
        for n in range(quad.N):
            np.testing.assert_allclose(result[n], result[0], rtol=1e-14)


# ─────────────────────────────────────────────────────────────────────
# Bit-equivalence vs legacy WhiteBoundaryOperator (L0)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestLegacyBitEquivalence:
    """The body is lifted verbatim from WhiteBoundaryOperator.apply;
    verify bit-identity on a non-trivial case so Wave 7 can confidently
    delegate ``WhiteBoundary`` to ``AngularAverageOperator`` via the
    SN realizer with no semantic drift."""

    def test_matches_white_bc_lebedev17(self):
        from orpheus.geometry.boundary import WhiteBoundaryOperator
        quad = LebedevSphere.create(17)
        bc = WhiteBoundaryOperator(axis="x", outward_sign=+1, albedo=1.0)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        rng = np.random.default_rng(123)
        psi = rng.uniform(0, 2, size=(quad.N, 5, 3))
        np.testing.assert_array_equal(op.apply(psi), bc.apply(psi, quad))
