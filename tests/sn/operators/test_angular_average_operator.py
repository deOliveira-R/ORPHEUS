"""Tests for AngularAverageOperator (Wave 1 / C1.1).

The operator implements the cosine-weighted Lambertian average over an
outgoing hemisphere; the incoming flux at a white surface equals
α × (this average).

**It is the RESPONSE, not a geometry.** Campaign phase B3.0 corrected
the older reading — "the §15.2 geometric projection G_diff in the
decomposition R_white = G_diff ⊗ α" — on a decidable criterion: a
geometry map is the composition operator of a measure-preserving
bijection, hence **multiplicative** (``G(ψφ) = (Gψ)(Gφ)``), and an
average is never that. So white's factors are ``G = IdentityMap``
(a white face fixes no geometry) and
``R = LambertianReemission(α)``, which is what this operator realizes.

The misreading was invisible for a structural reason worth keeping in
mind here: a rank-one response **annihilates** ``G`` entirely
(``R∘G = u ⊗ (Gᵀv)`` and the Lambertian's ``v = |Ω·n|`` is invariant
under both the mirror and the periodic translation), so white's
geometry slot had no observable consequence and the physics drifted
into it. Nothing this file measures changed when the factors moved.
"""

from __future__ import annotations

import numpy as np
import pytest

from orpheus.sn.boundary.angular import AngularAverageOperator
from orpheus.numerics.quadrature import Quadrature


# ─────────────────────────────────────────────────────────────────────
# Outgoing-hemisphere mask correctness (L0, parametrized over 3 quadratures)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestOutgoingHemisphereMask:
    """The cos_w array zeroes contributions from non-outgoing ordinates."""

    def test_lebedev_xmax_outgoing_only(self):
        quad = Quadrature.lebedev(17)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        # Apply to a flux that is +1 on outgoing (mu_x > 0) and -1 on incoming.
        # Cosine-weighted average must equal +1 (only outgoing contributes).
        psi = np.where(quad.mu_x > 0, +1.0, -1.0)
        result = op.apply(psi)
        # All ordinates carry the same broadcast value.
        np.testing.assert_allclose(result, +1.0, rtol=1e-14)

    def test_level_symmetric_ymin_outgoing_only(self):
        quad = Quadrature.level_symmetric(4)
        op = AngularAverageOperator.from_quadrature(quad, axis="y", outward_sign=-1)
        # ymin face: outgoing = mu_y < 0.
        psi = np.where(quad.mu_y < 0, +1.0, -1.0)
        result = op.apply(psi)
        np.testing.assert_allclose(result, +1.0, rtol=1e-14)

    def test_gauss_legendre_1d_xmax(self):
        quad = Quadrature.gauss_legendre(8)
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
        quad = Quadrature.lebedev(17)
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
        quad = Quadrature.level_symmetric(6)
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
        quad = Quadrature.lebedev(17)
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
        quad = Quadrature.level_symmetric(4)
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
# Structural predicates (L0)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestPredicates:
    def test_apply_only(self):
        # Rank-deficient projection: neither invertible nor adjointable —
        # apply is the only verb (apply itself is universal, guarded at
        # composition time).
        quad = Quadrature.lebedev(5)
        op = AngularAverageOperator.from_quadrature(quad, axis="z", outward_sign=+1)
        assert not op.is_invertible and not op.is_adjointable


# ─────────────────────────────────────────────────────────────────────
# Axis selection — z on Lebedev OK, z on GL raises (L0)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestAxisDispatch:
    def test_z_on_lebedev_succeeds(self):
        quad = Quadrature.lebedev(11)
        op = AngularAverageOperator.from_quadrature(quad, axis="z", outward_sign=+1)
        assert op.n_ordinates == quad.N

    def test_z_on_gauss_legendre_raises(self):
        quad = Quadrature.gauss_legendre(8)
        # 1-D Gauss-Legendre carries mu_z=zeros(N) (the axis is degenerate),
        # so AngularAverageOperator.from_quadrature reaches the "no outgoing
        # ordinates" guard rather than the "requires mu_z" early-return.
        # Both messages encode the same defect; pin the actual one.
        with pytest.raises(ValueError, match="no outgoing"):
            AngularAverageOperator.from_quadrature(quad, axis="z", outward_sign=+1)

    def test_invalid_axis_raises(self):
        quad = Quadrature.lebedev(5)
        with pytest.raises(ValueError, match="Unknown axis"):
            AngularAverageOperator.from_quadrature(quad, axis="w", outward_sign=+1)

    def test_invalid_outward_sign_raises(self):
        quad = Quadrature.lebedev(5)
        with pytest.raises(ValueError, match="outward_sign"):
            AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=0)


# ─────────────────────────────────────────────────────────────────────
# Defensive copy — operator is decoupled from the source quadrature (L0)
# ─────────────────────────────────────────────────────────────────────

@pytest.mark.l0
class TestDefensiveCopy:
    def test_quadrature_reference_not_held(self):
        quad = Quadrature.lebedev(5)
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
        quad = Quadrature.lebedev(5)
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
        quad = Quadrature.lebedev(5)  # N=14
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        with pytest.raises(ValueError, match="psi.shape"):
            op.apply(np.ones(13))

    def test_multi_dim_input_broadcasts(self):
        """Shape (N_ord, 5, 3) input — average is taken over leading axis only,
        broadcast back."""
        quad = Quadrature.lebedev(11)
        op = AngularAverageOperator.from_quadrature(quad, axis="x", outward_sign=+1)
        rng = np.random.default_rng(7)
        psi = rng.uniform(0, 1, size=(quad.N, 5, 3))
        result = op.apply(psi)
        assert result.shape == psi.shape
        # All leading-axis slices identical.
        for n in range(quad.N):
            np.testing.assert_allclose(result[n], result[0], rtol=1e-14)


# ─────────────────────────────────────────────────────────────────────
# Legacy-vs-realiser bit-equivalence test removed in Issue #186 (B3 + β2)
# ─────────────────────────────────────────────────────────────────────
#
# The pre-#186 TestLegacyBitEquivalence class compared
# ``AngularAverageOperator.apply(psi)`` to ``WhiteBoundary.apply(psi, quad)``.
# That legacy method no longer exists; the comparison would now be
# circular (the realiser-path test pins the same equivalence in
# tests/sn/test_sn_boundary_realizer.py). The Wave-1 cosine-weighted
# current conservation + self-adjointness tests above are the
# structurally-independent references for AngularAverageOperator
# correctness.
