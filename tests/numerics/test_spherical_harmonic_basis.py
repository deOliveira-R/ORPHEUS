r"""Tests for :class:`orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`.

The class's :meth:`evaluate` / :meth:`evaluate_from_components` methods
return the :math:`(N, L+1, 2L+1)` table of real spherical harmonics
under the no-:math:`4\pi/(2\ell+1)`-prefactor convention used by the
SN Galerkin reconstruction. The invariants below are the test hooks
the architecture catalog (§B-table) names for the harmonic basis:

* ``Y_0^0 = 1`` (P0 normalisation).
* Hard-coded P1 values: :math:`Y_1^{-1} = \mu_z`, :math:`Y_1^0 = \mu_x`,
  :math:`Y_1^{+1} = \mu_y`.
* Addition theorem :math:`\sum_m Y_\ell^m(\Omega) Y_\ell^m(\Omega') =
  P_\ell(\Omega \cdot \Omega')` — the identity consumed by the Pℓ
  source build.
* Discrete orthogonality on a Lebedev grid:
  :math:`\sum_n w_n Y_\ell^m(\Omega_n) Y_{\ell'}^{m'}(\Omega_n) =
  (4\pi / (2\ell+1)) \delta_{\ell\ell'} \delta_{mm'}` for
  :math:`\ell + \ell' \le \deg`.
"""
from __future__ import annotations

import numpy as np
import pytest
from scipy.special import eval_legendre

from orpheus.numerics.basis.spherical_harmonic_basis import (
    SphericalHarmonicBasis,
)
from orpheus.numerics.quadrature import lebedev_sphere


@pytest.mark.l0
class TestEvaluateRealSHBasic:
    """L0: hard-coded P0/P1 invariants."""

    def test_l_negative_rejected_at_construction(self):
        """Negative L is semantically meaningless for spherical harmonics;
        :class:`SphericalHarmonicBasis` rejects it at ``__post_init__``.

        Tightens the contract relative to the pre-P3.2 ``evaluate_real_sh``
        free function, which silently returned an empty array on
        ``L < 0`` — a permissive behaviour with no consumer."""
        with pytest.raises(ValueError, match="L must be non-negative"):
            SphericalHarmonicBasis(L=-1)

    def test_l_zero_p0_unity(self):
        mu_x = np.array([0.5, -0.5, 1.0])
        Y = SphericalHarmonicBasis(L=0).evaluate_from_components(
            mu_x, np.zeros(3), np.zeros(3),
        )
        assert Y.shape == (3, 1, 1)
        np.testing.assert_array_equal(Y[:, 0, 0], 1.0)

    def test_l_one_p1_hardcoded(self):
        mu_x = np.array([0.6, -0.3, 0.8])
        mu_y = np.array([0.0, 0.4, -0.5])
        mu_z = np.array([0.5, 0.7, 0.1])
        Y = SphericalHarmonicBasis(L=1).evaluate_from_components(mu_x, mu_y, mu_z)
        assert Y.shape == (3, 2, 3)
        np.testing.assert_array_equal(Y[:, 0, 0], 1.0)
        np.testing.assert_array_equal(Y[:, 1, 0], mu_z)   # m = -1
        np.testing.assert_array_equal(Y[:, 1, 1], mu_x)   # m =  0
        np.testing.assert_array_equal(Y[:, 1, 2], mu_y)   # m = +1


@pytest.mark.l1
class TestAdditionTheorem:
    r"""L1: :math:`\sum_m Y_\ell^m(\Omega) Y_\ell^m(\Omega') =
    P_\ell(\Omega \cdot \Omega')`.

    This is the structural identity the no-:math:`4\pi/(2\ell+1)`
    convention is *defined by*. Verifying it on a few random ordinate
    pairs gates the entire scaling factor.
    """

    @pytest.mark.parametrize("L", [2, 3, 4, 5])
    def test_addition_theorem_pairs(self, L):
        rng = np.random.default_rng(seed=2026 + L)
        basis = SphericalHarmonicBasis(L=L)
        # Sample two random unit vectors per trial, 8 trials.
        for trial in range(8):
            v = rng.standard_normal((2, 3))
            v /= np.linalg.norm(v, axis=1, keepdims=True)
            mu_x = v[:, 0]
            mu_y = v[:, 1]
            mu_z = v[:, 2]
            Y = basis.evaluate_from_components(mu_x, mu_y, mu_z)
            # cos(angle) = Ω · Ω'
            cos_alpha = float(np.dot(v[0], v[1]))
            for l in range(L + 1):
                lhs = float(np.sum(Y[0, l, : 2 * l + 1] * Y[1, l, : 2 * l + 1]))
                rhs = float(eval_legendre(l, cos_alpha))
                assert abs(lhs - rhs) < 1e-12, (
                    f"addition theorem mismatch at L={L}, l={l}, "
                    f"trial={trial}: lhs={lhs}, rhs={rhs}"
                )


@pytest.mark.l1
class TestDiscreteOrthogonalityOnLebedev:
    r"""L1: discrete orthogonality on a Lebedev quadrature.

    For the no-:math:`4\pi/(2\ell+1)`-prefactor convention,

    .. math::

        \sum_n w_n Y_\ell^m(\Omega_n) Y_{\ell'}^{m'}(\Omega_n)
        = \frac{4\pi}{2\ell+1} \delta_{\ell\ell'} \delta_{mm'}

    holds for :math:`\ell + \ell' \le \mathrm{deg}` where ``deg`` is
    the spherical-harmonic exactness degree of the Lebedev rule.
    """

    @pytest.mark.parametrize("order", [3, 5, 7])
    def test_lebedev_orthogonality(self, order):
        measure = lebedev_sphere(order)
        nodes = measure.nodes  # (N, 3)
        mu_x, mu_y, mu_z = nodes[:, 0], nodes[:, 1], nodes[:, 2]
        w = measure.weights
        deg = measure.degree_of_exactness
        L = deg // 2  # safe order so ℓ+ℓ' ≤ deg
        Y = SphericalHarmonicBasis(L=L).evaluate_from_components(mu_x, mu_y, mu_z)
        # Build flat (Y_lm) basis indexed by (l, m_offset).
        flat = []
        labels = []
        for l in range(L + 1):
            for m in range(-l, l + 1):
                flat.append(Y[:, l, l + m])
                labels.append((l, m))
        flat = np.array(flat)  # (n_basis, N)
        gram = np.einsum("an,bn,n->ab", flat, flat, w)
        for a, (la, ma) in enumerate(labels):
            for b, (lb, mb) in enumerate(labels):
                expected = (4.0 * np.pi / (2 * la + 1)) if (la == lb and ma == mb) else 0.0
                assert abs(gram[a, b] - expected) < 1e-10, (
                    f"orthogonality fails at "
                    f"(l,m)=({la},{ma}) vs ({lb},{mb}): "
                    f"got {gram[a, b]}, expected {expected}"
                )


@pytest.mark.l0
def test_on_axis_safe_no_division_by_zero():
    """When sin(theta) ≈ 0 (mu_x = ±1), the cos_phi/sin_phi expressions
    must NOT divide by zero. The implementation uses ``np.where`` to
    sentinel; verify L=3 evaluates without warnings on a pure-x input."""
    mu_x = np.array([1.0, -1.0, 0.0])
    mu_y = np.array([0.0, 0.0, 1.0])
    mu_z = np.array([0.0, 0.0, 0.0])
    with np.errstate(divide="raise", invalid="raise"):
        Y = SphericalHarmonicBasis(L=3).evaluate_from_components(mu_x, mu_y, mu_z)
    assert np.all(np.isfinite(Y))
