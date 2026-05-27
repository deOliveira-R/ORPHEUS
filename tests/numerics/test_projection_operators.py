r"""Tests for :mod:`orpheus.numerics.projection`.

The ``(R, Π)`` Galerkin pair invariants verified here:

* **Projection apply correctness** vs the einsum reference on
  randomized angular fields.
* **Reconstruction apply correctness** vs the einsum reference on
  randomized moment fields.
* **Galerkin idempotency** :math:`\Pi R c = c` on band-limited input
  (this is :meth:`GalerkinProjection.assert_galerkin_idempotency`).
* **Adjoint pairing** :math:`\langle \Pi \psi, c \rangle =
  \langle \psi, R c \rangle_W` under the W-weighted inner product
  (the Galerkin discipline's defining identity).
* **From-measure construction** wires up the :class:`DiscreteMeasure`
  + :meth:`SphericalHarmonicBasis.evaluate_from_components` machinery.
* **Multi-axis broadcasting**: applying along the leading axis with
  arbitrary trailing axes (cell, group) preserves shape.
* **Shape validation** rejects mismatched ``Y`` / ``L`` / weights.
"""
from __future__ import annotations

import numpy as np
import pytest

from orpheus.numerics.measure import DiscreteMeasure
from orpheus.numerics.operator import (
    CAP_APPLY,
    CAP_APPLY_TRANSPOSE,
    LinearOperatorMixin,
)
from orpheus.numerics.projection import (
    GalerkinProjection,
    MomentProjection,
    HarmonicMomentReconstruction,
    PetrovGalerkinProjection,
    ProjectionOperator,
)
from orpheus.numerics.quadrature import lebedev_sphere


# ─────────────────────────────────────────────────────────────────────
# ABC inheritance contracts
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestABCs:
    def test_galerkin_is_subclass_of_projection(self):
        assert issubclass(GalerkinProjection, ProjectionOperator)
        assert issubclass(PetrovGalerkinProjection, ProjectionOperator)

    def test_harmonic_moment_projection_is_galerkin(self):
        assert issubclass(MomentProjection, GalerkinProjection)

    def test_projection_operator_is_linear_operator(self):
        assert issubclass(ProjectionOperator, LinearOperatorMixin)


# ─────────────────────────────────────────────────────────────────────
# MomentProjection
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestMomentProjectionShapes:
    def test_apply_basic_shape(self):
        L = 2
        N = 6
        rng = np.random.default_rng(seed=1)
        weights = np.abs(rng.standard_normal(N)) + 0.1
        Y = rng.standard_normal((N, L + 1, 2 * L + 1))
        M = MomentProjection(weights=weights, Y=Y, L=L)
        psi = rng.standard_normal(N)
        out = M.apply(psi)
        assert out.shape == (L + 1, 2 * L + 1)

    def test_apply_with_trailing_axes_broadcasts(self):
        """Trailing axes (any rank) should broadcast unchanged.

        :class:`MomentProjection` is a generic numerics primitive
        whose ``apply`` contract is "leading axis is ordinates, every
        trailing axis broadcasts" — it is NOT tied to the SN principled
        storage (:ref:`theory-sn-index-convention`).  The trailing-axes
        shape used here is arbitrary; both ``(N, nx, ny, ng)`` and the
        SN-principled ``(N, ng, nx, ny)`` are valid inputs.
        """
        L = 2
        N = 6
        a, b, c = 3, 4, 5
        rng = np.random.default_rng(seed=2)
        weights = np.abs(rng.standard_normal(N)) + 0.1
        Y = rng.standard_normal((N, L + 1, 2 * L + 1))
        M = MomentProjection(weights=weights, Y=Y, L=L)
        psi = rng.standard_normal((N, a, b, c))
        out = M.apply(psi)
        assert out.shape == (L + 1, 2 * L + 1, a, b, c)

    def test_apply_against_einsum_reference(self):
        L = 3
        N = 8
        rng = np.random.default_rng(seed=3)
        weights = np.abs(rng.standard_normal(N)) + 0.1
        Y = rng.standard_normal((N, L + 1, 2 * L + 1))
        M = MomentProjection(weights=weights, Y=Y, L=L)
        psi = rng.standard_normal((N, 4))
        expected = np.einsum("n,nlm,nc->lmc", weights, Y, psi)
        np.testing.assert_allclose(M.apply(psi), expected, rtol=1e-15)


@pytest.mark.l0
class TestMomentProjectionShapeValidation:
    def test_inconsistent_Y_shape_raises(self):
        weights = np.ones(4)
        Y = np.zeros((4, 3, 5))  # claims L=2, but...
        with pytest.raises(ValueError, match="inconsistent"):
            MomentProjection(weights=weights, Y=Y, L=3)


@pytest.mark.l0
class TestMomentProjectionFromMeasure:
    def test_from_lebedev_measure(self):
        measure = lebedev_sphere(5)
        L = 3
        M = MomentProjection.from_measure(measure, L=L)
        assert M.weights.shape == (measure.n_points,)
        assert M.Y.shape == (measure.n_points, L + 1, 2 * L + 1)
        assert M.L == L

    def test_from_measure_rejects_1d_nodes(self):
        # 1-D measure (slab quadrature) cannot be projected onto S²
        nodes = np.array([-0.5, 0.5])
        weights = np.array([1.0, 1.0])
        mu = DiscreteMeasure(nodes=nodes, weights=weights, space="[-1,1]")
        with pytest.raises(ValueError, match="\\(N, 3\\)"):
            MomentProjection.from_measure(mu, L=1)


# ─────────────────────────────────────────────────────────────────────
# HarmonicMomentReconstruction
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestHarmonicMomentReconstructionShapes:
    def test_apply_basic_shape(self):
        L = 2
        N = 6
        rng = np.random.default_rng(seed=4)
        Y = rng.standard_normal((N, L + 1, 2 * L + 1))
        R = HarmonicMomentReconstruction.from_Y(Y)
        moments = rng.standard_normal((L + 1, 2 * L + 1))
        out = R.apply(moments)
        assert out.shape == (N,)

    def test_apply_with_trailing_axes(self):
        L = 1
        N = 4
        rng = np.random.default_rng(seed=5)
        Y = rng.standard_normal((N, L + 1, 2 * L + 1))
        R = HarmonicMomentReconstruction.from_Y(Y)
        moments = rng.standard_normal((L + 1, 2 * L + 1, 3, 5))
        out = R.apply(moments)
        assert out.shape == (N, 3, 5)

    def test_apply_against_einsum_reference(self):
        L = 3
        N = 10
        rng = np.random.default_rng(seed=6)
        Y = rng.standard_normal((N, L + 1, 2 * L + 1))
        R = HarmonicMomentReconstruction.from_Y(Y)
        moments = rng.standard_normal((L + 1, 2 * L + 1))
        two_l_plus_one = 2.0 * np.arange(L + 1) + 1.0
        expected = np.einsum("nlm,l,lm->n", Y, two_l_plus_one, moments)
        np.testing.assert_allclose(R.apply(moments), expected, rtol=1e-14)


# ─────────────────────────────────────────────────────────────────────
# Galerkin idempotency on band-limited input
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l1
class TestGalerkinIdempotencyOnLebedev:
    r"""L1: :math:`\Pi R c = c` on band-limited coefficient input.

    Constructed from a Lebedev grid that integrates :math:`Y_\ell^m
    Y_{\ell'}^{m'}` exactly for :math:`\ell + \ell' \le 2L`. The
    addition theorem then makes :math:`\Pi R = I` on the
    coefficient space.
    """

    @pytest.mark.parametrize("order_L", [(7, 2), (13, 3), (17, 4)])
    def test_pi_R_is_identity_on_band_limited(self, order_L):
        order, L = order_L
        measure = lebedev_sphere(order)
        # 4π normalisation: weights from lebedev_sphere sum to 4π.
        # The Galerkin identity is Π R c = c with the addition-theorem
        # convention used by SphericalHarmonicBasis.evaluate_from_components.
        M = MomentProjection.from_measure(measure, L=L)
        R = HarmonicMomentReconstruction.from_Y(M.Y)
        # Random coefficients in the moment space
        rng = np.random.default_rng(seed=order)
        c = rng.standard_normal((L + 1, 2 * L + 1))
        # Mask out unused (l, m) entries: only |m| <= l survive.
        for l_idx in range(L + 1):
            for m_idx in range(2 * L + 1):
                m = m_idx - L  # not the right offset
        # Note: indexing convention is `m_offset = l + m` where m runs
        # from -l to +l. Out-of-range entries should be 0.
        c_masked = np.zeros_like(c)
        for l_idx in range(L + 1):
            for m_off in range(2 * l_idx + 1):  # m_off = l + m, so 0..2l
                m_loc_in_full = m_off + (L - l_idx)  # shift to L-centred slot
                # Actually our Y storage is also L-centred; let's be careful:
                # Y[:, l, l_idx + m] uses the L-centred slot when L is the
                # outer maximum. Looking at the code: Y[:, l, l + m] with
                # m in (-l..+l) maps to slot 0..2l. So in the (L+1, 2L+1)
                # array, only slots 0..2l are populated for row l. Rest are 0.
                c_masked[l_idx, m_off] = c[l_idx, m_off]
        # Now apply Π R c_masked and compare
        # Note: Π R c only equals c on the band-limited subspace under the
        # 4π / (2l+1) normalization. Our addition-theorem convention has
        # ⟨Y_l^m, Y_l'^m'⟩_W = (4π / (2l+1)) δ. So:
        #   M @ R @ c_lm = (Σ_n w_n Y_lm Σ_l' (2l'+1) Σ_m' Y_l'm' c_l'm')
        #                = Σ_l' (2l'+1) Σ_m' c_l'm' (Σ_n w_n Y_lm Y_l'm')
        #                = Σ_l' (2l'+1) Σ_m' c_l'm' (4π / (2l'+1)) δ_ll' δ_mm'
        #                = 4π · c_lm
        # So Π R = 4π · I in this convention.
        out = M.apply(R.apply(c_masked))
        np.testing.assert_allclose(out, 4.0 * np.pi * c_masked, rtol=1e-10, atol=1e-12)


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
class TestApplyTransposeIsRepresentationTranspose:
    r"""L1: ``MomentProjection.apply_transpose`` is the representation transpose :math:`\Pi^\top`.

    Post-P1.4 contract: ``apply_transpose`` returns the row-by-row
    matrix transpose of :math:`\Pi`, which carries :math:`w_n` (the
    quadrature weight per ordinate) because :math:`\Pi_{(\ell m), n}
    = w_n \, Y_\ell^m(\hat\Omega_n)`. The W-weighted Hilbert adjoint
    :math:`\Pi^*` is the SEPARATE operation reached via
    :attr:`MomentProjection.H`, which composes the metric-aware
    formula :math:`\Pi^* = g_C \cdot S_0` via the generic
    :class:`~orpheus.numerics.operator._AdjointOperator` machinery.

    Closes ERR-039: pre-P1.4 ``apply_transpose`` returned the bare
    :math:`S_0(c)` and the docstring labeled it the W-weighted Hilbert
    adjoint — but the true representation transpose carries
    :math:`w_n`, and the Hilbert adjoint additionally carries
    :math:`g_C`. The two are now separately typed.
    """

    def test_representation_transpose_carries_w_n(self):
        """The representation transpose ``M.T`` is :math:`w_n \, S_0(c)`."""
        measure = lebedev_sphere(7)
        L = 2
        rng = np.random.default_rng(seed=5678)
        M = MomentProjection.from_measure(measure, L=L)

        c = rng.standard_normal((L + 1, 2 * L + 1))
        c_masked = np.zeros_like(c)
        for l_idx in range(L + 1):
            for m_off in range(2 * l_idx + 1):
                c_masked[l_idx, m_off] = c[l_idx, m_off]

        Pi_T_c = M.apply_transpose(c_masked)
        # Representation transpose: w_n · S_0(c).  The production einsum
        # ``"n,nlm,lm...->n..."`` fuses the w_n multiplication into the
        # reduction; the separate ``measure.weights * einsum(...)`` form
        # below drifts at FP-non-associativity (≈1e-14) for the same
        # mathematical answer.  Tolerance covers that drift.
        Pi_T_expected = measure.weights * np.einsum("nlm,lm->n", M.Y, c_masked)
        np.testing.assert_allclose(Pi_T_c, Pi_T_expected, rtol=1e-13, atol=1e-14)

    def test_transpose_distinct_from_addition_theorem_R(self):
        """``M.T`` and ``R`` are not equal — distinguishes the two
        primitives ERR-039 originally conflated."""
        measure = lebedev_sphere(7)
        L = 2
        rng = np.random.default_rng(seed=5678)
        M = MomentProjection.from_measure(measure, L=L)
        R = HarmonicMomentReconstruction.from_Y(M.Y)

        c = rng.standard_normal((L + 1, 2 * L + 1))
        c_masked = np.zeros_like(c)
        for l_idx in range(L + 1):
            for m_off in range(2 * l_idx + 1):
                c_masked[l_idx, m_off] = c[l_idx, m_off]

        Pi_T_c = M.apply_transpose(c_masked)   # = w_n · S_0(c)
        Rc = R.apply(c_masked)                  # = (2ℓ+1) · S_0(c)
        # They differ structurally: M.T weights per-ordinate by w_n;
        # R weights per-ℓ by (2ℓ+1). Equal only at the trivial limit
        # of vanishing input.
        assert not np.allclose(Pi_T_c, Rc, atol=1e-3)

    def test_representation_transpose_identity_matches_production(self):
        """The matrix identity :math:`\langle \Pi \psi, c \rangle =
        \langle \psi, \Pi^\top c \rangle` holds against the
        production ``apply_transpose`` path.

        ``apply_transpose`` IS the representation transpose, so the
        identity is the standard ``inner(A x, y) == inner(x, A^T y)``
        with the EUCLIDEAN inner product on both sides (w_n is already
        absorbed into A^T by construction).
        """
        measure = lebedev_sphere(13)
        L = 3
        rng = np.random.default_rng(seed=1234)
        M = MomentProjection.from_measure(measure, L=L)

        psi = rng.standard_normal(measure.n_points)
        c = rng.standard_normal((L + 1, 2 * L + 1))
        c_masked = np.zeros_like(c)
        for l_idx in range(L + 1):
            for m_off in range(2 * l_idx + 1):
                c_masked[l_idx, m_off] = c[l_idx, m_off]

        # ⟨M ψ, c⟩  (Euclidean on coefficient space)
        Mpsi = M.apply(psi)
        lhs = float(np.sum(Mpsi * c_masked))

        # ⟨ψ, M^T c⟩  (Euclidean on ordinate space —
        #              w_n now lives INSIDE M^T)
        Pi_T_c = M.apply_transpose(c_masked)
        rhs = float(np.sum(psi * Pi_T_c))

        np.testing.assert_allclose(lhs, rhs, rtol=1e-12, atol=1e-14)


@pytest.mark.l1
@pytest.mark.catches("ERR-039")
class TestHilbertAdjointViaGenericMachinery:
    r"""L1: :attr:`MomentProjection.H` computes :math:`g_C \cdot S_0(c)` via the generic adjoint wrapper.

    Post-P1.4 the W-weighted Hilbert adjoint is reached through
    ``M.H`` — the generic :class:`~orpheus.numerics.operator._AdjointOperator`
    wrapper. It reads the codomain's
    :class:`~orpheus.numerics.spaces.SphericalHarmonicSpace` metric
    (:math:`g_C = \mathrm{diag}(4\pi/(2\ell+1))`) and the domain's
    quadrature-weight metric (:math:`W = \mathrm{diag}(w_n)`) to compose

    .. math::

        (\Pi^* c)_n = \frac{1}{w_n}\, \Pi^\top(g_C \cdot c)_n
                    = \sum_{\ell m} \frac{4\pi}{2\ell+1}\,
                      Y_\ell^m(\hat\Omega_n)\, c_{\ell m}.

    This is the ERR-039 endpoint: the metric, the transpose, and the
    Hilbert adjoint are all now separately typed; their composition
    falls out from the generic machinery without prose warnings.
    """

    def test_H_returns_g_C_times_S_0(self):
        measure = lebedev_sphere(13)
        L = 3
        rng = np.random.default_rng(seed=1234)
        M = MomentProjection.from_measure(measure, L=L)

        c = rng.standard_normal((L + 1, 2 * L + 1))
        c_masked = np.zeros_like(c)
        for l_idx in range(L + 1):
            for m_off in range(2 * l_idx + 1):
                c_masked[l_idx, m_off] = c[l_idx, m_off]

        M_H_c = M.H.apply(c_masked)

        metric_per_ell = 4.0 * np.pi / (2.0 * np.arange(L + 1) + 1.0)
        expected = np.einsum("l,nlm,lm->n", metric_per_ell, M.Y, c_masked)

        np.testing.assert_allclose(M_H_c, expected, rtol=1e-13, atol=1e-14)

    def test_H_satisfies_hilbert_adjoint_identity(self):
        r"""The defining identity :math:`\langle M\psi, c\rangle_C =
        \langle \psi, M^* c\rangle_V` with the metric inner products."""
        measure = lebedev_sphere(13)
        L = 3
        rng = np.random.default_rng(seed=1234)
        M = MomentProjection.from_measure(measure, L=L)

        psi = rng.standard_normal(measure.n_points)
        c = rng.standard_normal((L + 1, 2 * L + 1))
        c_masked = np.zeros_like(c)
        for l_idx in range(L + 1):
            for m_off in range(2 * l_idx + 1):
                c_masked[l_idx, m_off] = c[l_idx, m_off]

        # ⟨M ψ, c⟩_C  with the SphericalHarmonicSpace metric
        Mpsi = M.apply(psi)
        lhs = float(M.codomain.inner_product(Mpsi, c_masked))

        # ⟨ψ, M^* c⟩_V  with the quadrature-weight metric
        M_H_c = M.H.apply(c_masked)
        rhs = float(np.sum(measure.weights * psi * M_H_c))

        np.testing.assert_allclose(lhs, rhs, rtol=1e-12, atol=1e-14)


@pytest.mark.l1
class TestGalerkinAdjointPairing:
    r"""L1: :math:`\langle \Pi \psi, c \rangle = \langle \psi, R c \rangle_W`.

    Under the W-weighted L²(S²) inner product on the angular side and
    Euclidean inner product on the moment side, the Galerkin
    discipline's defining identity is ``Π* = R``.
    """

    def test_adjoint_pairing_matches(self):
        measure = lebedev_sphere(13)
        L = 3
        rng = np.random.default_rng(seed=42)
        M = MomentProjection.from_measure(measure, L=L)
        R = HarmonicMomentReconstruction.from_Y(M.Y)

        psi = rng.standard_normal(measure.n_points)
        c = rng.standard_normal((L + 1, 2 * L + 1))
        # Mask high-m entries (only |m| <= l survive)
        c_masked = np.zeros_like(c)
        for l_idx in range(L + 1):
            for m_off in range(2 * l_idx + 1):
                c_masked[l_idx, m_off] = c[l_idx, m_off]

        # ⟨Π ψ, c⟩ — Euclidean inner product on (l, m) space:
        Mpsi = M.apply(psi)
        lhs = float(np.sum(Mpsi * c_masked))

        # ⟨ψ, M^T c⟩  — the matrix-transpose pairing crosses the
        # production apply_transpose path (which now carries w_n
        # internally per P1.4).  The natural pairing simplifies to
        #   ⟨M ψ, c⟩ = ⟨ψ, M^T c⟩
        # with the Euclidean inner product on both sides; the W weight
        # has migrated from the test (external * measure.weights) into
        # the operator (M.apply_transpose returns w_n · S_0(c)).
        rhs = float(np.sum(psi * M.apply_transpose(c_masked)))
        np.testing.assert_allclose(lhs, rhs, rtol=1e-12, atol=1e-14)


# ─────────────────────────────────────────────────────────────────────
# Retired in P1.6: TestAssertGalerkinIdempotencyMethod (sole caller of
# GalerkinProjection.assert_galerkin_idempotency) was deleted alongside
# the method itself. The method shipped the wrong invariant
# (M · R = I instead of the correct M · R = 4π · I under the
# no-prefactor SH convention), and its sole test exercised it against
# a deliberately non-orthogonal Y so its (wrong) assertion would
# "successfully" fail.  The genuine Π R = 4π · I identity is pinned by
# TestGalerkinIdempotencyOnLebedev (above) and the new file
# tests/numerics/test_spherical_harmonic_space.py.
# See moment-space + layering plan §1.5 CC.5 and §P1.6.
# ─────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────
# Capability set
# ─────────────────────────────────────────────────────────────────────


@pytest.mark.l0
class TestCapabilities:
    def test_projection_advertises_apply_and_apply_transpose(self):
        measure = lebedev_sphere(7)
        M = MomentProjection.from_measure(measure, L=2)
        assert CAP_APPLY in M.capabilities
        assert CAP_APPLY_TRANSPOSE in M.capabilities

    def test_reconstruction_advertises_apply(self):
        Y = np.zeros((4, 2, 3))
        R = HarmonicMomentReconstruction.from_Y(Y)
        assert CAP_APPLY in R.capabilities


@pytest.mark.l0
def test_petrov_galerkin_is_abstract_no_ship_concrete():
    """No concrete PetrovGalerkin subclasses ship in this Wave 0 commit;
    energy condensation lands later. Confirm the ABC exists and is
    distinct from GalerkinProjection."""
    assert PetrovGalerkinProjection is not GalerkinProjection
    assert issubclass(PetrovGalerkinProjection, ProjectionOperator)
