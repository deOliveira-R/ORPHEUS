r"""Paired symbolic-vs-production contract test for the cylinder
mode-N Knyazev :math:`\mathrm{Ki}_{2+k}` expansion.

Math-origin pattern (mirroring
:func:`orpheus.derivations.sn_balance.derive_cumprod_recurrence`): the
SymPy derivation in :mod:`orpheus.derivations.peierls_cylinder_knyazev`
is the **source of truth** for the rank-:math:`N` cylinder
:math:`P_{\rm esc}^{(n,3d)}` and :math:`G_{\rm bc}^{(n,3d)}`
primitives. This test consumes the SymPy origin as a contract and pins:

1. The polar integral identity :eq:`peierls-cyl-knyazev-polar-id`
   (numerical match to :func:`ki_n_mp`).
2. The :math:`1/\pi` and :math:`4/\pi` prefactors for P and G.
3. The mode-0 reduction of the Knyazev primitive equals the legacy
   :func:`compute_P_esc` cylinder branch bit-exactly.
4. The mode-1, 2, 3 lambdified symbolic expansion equals the
   production :func:`compute_P_esc_cylinder_3d_mode` to within
   quadrature precision.

The verified equation is :eq:`peierls-cyl-3d-mode-formula` (cylinder
mode-N :math:`P_{\rm esc}^{(n,3d)}` formula, Knyazev expansion).

The canonical derivation sources are
``orpheus/derivations/peierls_cylinder_knyazev.py:1`` (SymPy module
docstring) and ``orpheus/derivations/peierls_geometry.py:1927`` (the
shipped :func:`compute_P_esc_cylinder_3d_mode`).
"""
from __future__ import annotations

import numpy as np
import pytest
import sympy as sp

from orpheus.derivations._kernels import ki_n_mp
from orpheus.derivations._shifted_legendre import (
    shifted_legendre_monomial_coefs,
)
from orpheus.derivations.peierls_cylinder_knyazev import (
    derive_g_prefactor,
    derive_p_prefactor,
    derive_polar_integral_identity,
)
from orpheus.derivations.peierls_geometry import (
    CYLINDER_1D,
    compute_P_esc,
    compute_P_esc_cylinder_3d_mode,
)


# ═══════════════════════════════════════════════════════════════════════
# Derivation-level identities
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-cyl-3d-mode-formula")
@pytest.mark.parametrize("k", [0, 1, 2, 3])
@pytest.mark.parametrize("x_val", [0.1, 0.5, 1.0, 5.0])
def test_polar_integral_identity(k, x_val):
    r"""Verify :eq:`peierls-cyl-knyazev-polar-id`:
    :math:`\int_0^{\pi/2} \sin^{k+1}\theta\,e^{-x/\sin\theta}\,
    \mathrm d\theta = \mathrm{Ki}_{k+2}(x)` numerically.

    Compares the symbolic integrand from
    :func:`derive_polar_integral_identity` (lambdified and integrated
    via :func:`scipy.integrate.quad`) against the canonical
    :func:`ki_n_mp` Bickley-Naylor implementation.
    """
    from scipy.integrate import quad

    sym = derive_polar_integral_identity(k)
    integrand = sp.lambdify(
        (sym["theta_p"], sym["x"]),
        sym["theta_form"],
        "math",
    )
    val_quad, _ = quad(integrand, 0.0, np.pi / 2, args=(x_val,),
                       epsabs=1e-13, epsrel=1e-12)
    val_ki = float(ki_n_mp(k + 2, x_val, precision_digits=30))
    assert val_quad == pytest.approx(val_ki, rel=1e-10), (
        f"Polar integral identity fails at k={k}, x={x_val}: "
        f"quad={val_quad:.16e}, Ki_{k+2}={val_ki:.16e}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-cyl-3d-mode-formula")
def test_p_prefactor_is_1_over_pi():
    r"""The cylinder P-prefactor from
    :func:`derive_p_prefactor` is :math:`1/\pi` after simplification.
    """
    val = derive_p_prefactor()
    diff = sp.simplify(val - sp.Rational(1) / sp.pi)
    assert diff == 0, (
        f"P prefactor != 1/π: got {val}, diff = {diff}"
    )


@pytest.mark.l1
@pytest.mark.verifies("peierls-cyl-3d-mode-formula")
def test_g_prefactor_is_4_over_pi():
    r"""The cylinder G-prefactor from
    :func:`derive_g_prefactor` is :math:`4/\pi` after simplification.
    """
    val = derive_g_prefactor()
    diff = sp.simplify(val - sp.Rational(4) / sp.pi)
    assert diff == 0, (
        f"G prefactor != 4/π: got {val}, diff = {diff}"
    )


# ═══════════════════════════════════════════════════════════════════════
# Production parity — n = 0 reduces to compute_P_esc cylinder branch
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.l1
@pytest.mark.verifies("peierls-cyl-3d-mode-formula")
def test_compute_P_esc_cylinder_3d_mode_n0_matches_compute_P_esc():
    r"""For :math:`n = 0`, :math:`c_0^0 = 1` and only :math:`k = 0`
    survives in the Knyazev sum, so
    :math:`P_{\rm esc}^{(0,3d)} \equiv (1/\pi)\!\int_0^\pi
    \mathrm{Ki}_2(\tau_{\rm 2D})\,\mathrm d\omega`. The legacy
    :func:`compute_P_esc` cylinder branch uses the same integrand
    (its ``escape_kernel_mp`` for cylinder is
    :math:`(1/\pi)\,\mathrm{Ki}_2`). Therefore mode-0 Knyazev must
    bit-equal :func:`compute_P_esc` to within quadrature precision.
    """
    radii = np.array([5.0])
    sig_t = np.array([1.0])
    r_nodes = np.linspace(0.1, 4.9, 5)

    P_legacy = compute_P_esc(
        CYLINDER_1D, r_nodes, radii, sig_t,
        n_angular=32, dps=25,
    )
    P_mode0 = compute_P_esc_cylinder_3d_mode(
        CYLINDER_1D, r_nodes, radii, sig_t, n_mode=0,
        n_angular=32, dps=25,
    )
    np.testing.assert_allclose(
        P_mode0, P_legacy, rtol=1e-12, atol=1e-13,
        err_msg=(
            "compute_P_esc_cylinder_3d_mode(n=0) drifted from legacy "
            "compute_P_esc cylinder branch"
        ),
    )


# ═══════════════════════════════════════════════════════════════════════
# Production parity — n >= 1 matches lambdified symbolic expansion
# ═══════════════════════════════════════════════════════════════════════


def _knyazev_p_symbolic_value(
    r_i: float, R: float, sig_t: float, n_mode: int,
    n_quad: int = 256,
) -> float:
    r"""Reference SymPy → numerical evaluation of
    :math:`P_{\rm esc}^{(n,3d)}(r_i)` via direct quadrature on the
    symbolic expansion at the homogeneous limit.

    Uses :func:`shifted_legendre_monomial_coefs` (the canonical
    coefficient source) and :func:`ki_n_mp` (the canonical Bickley
    primitive). Computed at high :math:`n_{\rm quad}` to give a
    reference for the production primitive; quadrature precision
    will dominate the difference at the production-default
    :math:`n_{\rm angular} = 32`.
    """
    from scipy.integrate import quad

    coefs = shifted_legendre_monomial_coefs(n_mode)

    def integrand(omega: float) -> float:
        cos_om = float(np.cos(omega))
        sin_om = float(np.sin(omega))
        disc = max(R * R - r_i * r_i * sin_om * sin_om, 0.0)
        d_2d = -r_i * cos_om + np.sqrt(disc)
        if d_2d <= 0.0:
            return 0.0
        tau = sig_t * d_2d
        mu_2d = (r_i * cos_om + d_2d) / R
        return sum(
            c_k * (mu_2d ** k) * float(ki_n_mp(k + 2, tau, 30))
            for k, c_k in enumerate(coefs) if c_k != 0.0
        )

    val, _ = quad(integrand, 0.0, np.pi, limit=n_quad, epsabs=1e-13,
                  epsrel=1e-12)
    return float(val / np.pi)


@pytest.mark.l1
@pytest.mark.verifies("peierls-cyl-3d-mode-formula")
@pytest.mark.parametrize("n_mode", [1, 2, 3])
def test_compute_P_esc_cylinder_3d_mode_n123_matches_symbolic_expansion(n_mode):
    r"""For :math:`n \in \{1, 2, 3\}`, the production
    :func:`compute_P_esc_cylinder_3d_mode` matches the lambdified
    symbolic Knyazev expansion at a homogeneous cylinder fixture
    (:math:`R = 5\,{\rm cm}`, :math:`\Sigma_t = 1`).

    The reference uses :func:`scipy.integrate.quad` (adaptive) on the
    full symbolic integrand
    :math:`(1/\pi)\!\int_0^\pi \sum_k c_n^k\,\mu_{\rm 2D}^k\,
    \mathrm{Ki}_{k+2}(\tau_{\rm 2D})\,\mathrm d\omega`; the
    production uses :func:`per_observer_angular_assembly` with
    :math:`n_{\rm angular} = 64` panels per region. The mismatch
    floor is the production quadrature precision (~1e-9 to 1e-10
    for smooth cylinder integrands at this :math:`n_{\rm angular}`).
    """
    radii = np.array([5.0])
    sig_t = np.array([1.0])
    r_nodes = np.array([0.5, 2.5, 4.5])

    P_prod = compute_P_esc_cylinder_3d_mode(
        CYLINDER_1D, r_nodes, radii, sig_t, n_mode=n_mode,
        n_angular=64, dps=25,
    )
    P_ref = np.array([
        _knyazev_p_symbolic_value(float(r_i), float(radii[0]),
                                  float(sig_t[0]), n_mode)
        for r_i in r_nodes
    ])
    np.testing.assert_allclose(
        P_prod, P_ref, rtol=2e-8, atol=1e-10,
        err_msg=(
            f"compute_P_esc_cylinder_3d_mode(n={n_mode}) drifted "
            f"from lambdified symbolic Knyazev expansion."
        ),
    )
