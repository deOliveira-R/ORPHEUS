r"""Atalay 1997 Eq 12 dispersion relation — discrete Case eigenvalue
:math:`\nu_0 = i u_0` for one-pair-of-modes regime.

For the linearly anisotropic transport operator (Atalay Eq 1), the
dispersion relation Eq 11

.. math::

   \Lambda(\nu) = d(\nu^2)\big[1 - c\nu\,\mathrm{tanh}^{-1}(1/\nu)\big]
                  - 3 f_1 (1-c)^2 \nu^2 = 0

becomes, on substituting :math:`\nu_0 = i u_0` (purely imaginary,
:math:`u_0 > 0`) and using
:math:`\nu_0\,\mathrm{tanh}^{-1}(1/\nu_0) = u_0\,\mathrm{atan}(1/u_0)`
(real), the **real** transcendental equation

.. math::

   \Lambda_{\rm real}(u_0; c, f_1) =
       -3 f_1 u_0^2 [u_0\,\mathrm{atan}(1/u_0) - 1]\,c^2
       + \{3 f_1 u_0^2 [u_0\,\mathrm{atan}(1/u_0) - 1] -
            u_0\,\mathrm{atan}(1/u_0)\}\,c
       + 1 = 0 .

This is the **Eq 12** quadratic-in-c form, but here we use it as a
transcendental in :math:`u_0` for fixed :math:`(c, f_1)`. Root-find
on :math:`u_0 > 0` via :class:`scipy.optimize.brentq`.

Alternative
-----------

Atalay's Table 1 gives :math:`\bar\nu(c)` for :math:`f_1 = 0`. The
isotropic limit reduces to
:math:`c u_0 \mathrm{atan}(1/u_0) = 1`, the classical Case-Zweifel
relation, whose root we use as the initial bracket for the
anisotropic root-find.

Limit checks
------------

* :math:`c \to 1^+`: :math:`u_0 \to \infty`. Root behaves as
  :math:`u_0 \sim \sqrt{1/(3(c-1))}` for isotropic — the diffusion
  limit of Case eigenvalue.
* :math:`c \to c_{\max} = 1 + 1/(3 f_1)`: :math:`u_0 \to 0^+`. At
  this boundary the dispersion relation has a double root, and beyond
  it complex eigenvalues appear (Atalay Eq 5).
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq


@dataclass(frozen=True)
class CaseDispersionResult:
    """Result of :func:`case_atalay_u0`.

    Attributes
    ----------
    u0 : float
        Magnitude of the purely-imaginary discrete Case eigenvalue,
        :math:`\\nu_0 = i u_0`. Real positive.
    c : float
        Mean number of secondaries per collision.
    f1 : float
        Linearly-anisotropic mean cosine.
    Lambda_residual : float
        :math:`\\Lambda_{\\rm real}(u_0)` at the solved root. Should
        be small (≲ 1e-12).
    bracket_lo : float
        Lower bracket value used by Brent's method.
    bracket_hi : float
        Upper bracket value used by Brent's method.
    """

    u0: float
    c: float
    f1: float
    Lambda_residual: float
    bracket_lo: float
    bracket_hi: float


def _u_atan_inv_u(u: float) -> float:
    r"""Compute :math:`u \cdot \mathrm{atan}(1/u)` for :math:`u > 0`.

    For :math:`u \to \infty`, :math:`u\,\mathrm{atan}(1/u) \to 1` (the
    leading asymptote :math:`1/u \cdot u = 1`).

    For :math:`u \to 0^+`, :math:`u\,\mathrm{atan}(1/u) \to 0` (the
    :math:`\mathrm{atan}` saturates at :math:`\pi/2`, multiplied by
    :math:`u \to 0`).
    """
    if u <= 0:
        raise ValueError(f"u must be positive, got {u}")
    return u * math.atan(1.0 / u)


def _Lambda_real(u0: float, c: float, f1: float) -> float:
    r"""Atalay Eq 12 evaluated at :math:`\nu_0 = i u_0` (purely imaginary).

    Returns the real-valued residual. A root :math:`u_0` of this
    function is the magnitude of the discrete Case eigenvalue.
    """
    L = _u_atan_inv_u(u0)
    Q = L - 1.0  # u_0 · atan(1/u_0) - 1; negative for u_0 > 0 (since atan < π/2 for u_0 > 0... wait, see below)
    # Note: for u_0 ∈ (0, ∞), atan(1/u_0) ∈ (0, π/2). So u_0·atan(1/u_0) is
    # positive; for u_0 = 1, value is π/4 ≈ 0.785; for u_0 = 0.5, value is
    # 0.5 · atan(2) ≈ 0.554; for u_0 = 0.1, value is 0.1 · atan(10) ≈ 0.147;
    # for u_0 = 10, value is 10 · atan(0.1) ≈ 0.997 → asymptote 1.
    # So Q < 0 for u_0 ∈ (0, ∞) and Q → -1 as u_0 → 0, Q → 0⁻ as u_0 → ∞.
    return -3 * f1 * u0**2 * Q * c**2 + (3 * f1 * u0**2 * Q - L) * c + 1.0


def case_atalay_dispersion_quadratic_coeffs(
    u0: float, c: float, f1: float,
) -> tuple[float, float, float]:
    r"""Return the (a, b, c_const) coefficients of Atalay Eq 12 evaluated
    at :math:`\nu_0 = i u_0`.

    Eq 12 reads :math:`a c^2 + b c + c_{const} = 0` where
    :math:`a = 3 f_1 \nu_0^2 [\nu_0 \mathrm{tanh}^{-1}(1/\nu_0) - 1]`,
    :math:`b = -[3 f_1 \nu_0^2 (\dots) + \nu_0\mathrm{tanh}^{-1}(1/\nu_0)]`,
    :math:`c_{const} = 1`. With :math:`\nu_0 = i u_0`, :math:`\nu_0^2 = -u_0^2`
    and :math:`\nu_0\mathrm{tanh}^{-1}(1/\nu_0) = u_0 \mathrm{atan}(1/u_0)` (real).
    Hence:

    .. math::

       a &= -3 f_1 u_0^2 [u_0\mathrm{atan}(1/u_0) - 1] \\
       b &= +3 f_1 u_0^2 [u_0\mathrm{atan}(1/u_0) - 1] - u_0\mathrm{atan}(1/u_0) \\
       c_{const} &= 1 .

    Useful for cross-checks (the c value satisfies the quadratic).
    """
    L = _u_atan_inv_u(u0)
    Q = L - 1.0
    a = -3 * f1 * u0**2 * Q
    b = 3 * f1 * u0**2 * Q - L
    return (a, b, 1.0)


def case_atalay_u0(
    c: float,
    f1: float = 0.0,
    *,
    u0_min: float = 1e-9,
    u0_max: float = 1.0e6,
    bisect_tol: float = 1e-13,
) -> CaseDispersionResult:
    r"""Find the discrete Case eigenvalue :math:`u_0 > 0` for given
    :math:`(c, f_1)`.

    Solves :math:`\Lambda_{\rm real}(u_0; c, f_1) = 0` by Brent's method.

    Parameters
    ----------
    c : float
        Mean number of secondaries per collision. Must satisfy
        :math:`c > 1` (multiplying medium) and
        :math:`c \le 1 + 1/(3 f_1)` (validity bound, Atalay Eq 5).
    f1 : float, default 0.0
        Linearly-anisotropic mean cosine. ``f1 = 0`` is isotropic.
    u0_min, u0_max : float, default 1e-9, 1e6
        Bracket for the Brent root-find. Must contain a sign change.
    bisect_tol : float, default 1e-13
        Absolute tolerance on :math:`u_0`.

    Returns
    -------
    :class:`CaseDispersionResult`

    Raises
    ------
    ValueError
        If :math:`c \le 1`, if :math:`c` violates the validity bound
        (Atalay Eq 5), or if no sign change of :math:`\Lambda_{\rm real}`
        is found in :math:`[u_{0\,\min}, u_{0\,\max}]`.
    """
    if c <= 1.0:
        raise ValueError(f"case_atalay_u0 requires c > 1, got c={c}")
    if f1 < 0.0:
        raise ValueError(f"f1 must be ≥ 0, got {f1}")
    if f1 > 0:
        c_max = 1.0 + 1.0 / (3.0 * f1)
        if c > c_max:
            raise ValueError(
                f"c={c} violates Atalay Eq 5 validity bound "
                f"c <= 1 + 1/(3 f_1) = {c_max} for f_1={f1}. "
                f"Outside this band complex eigenvalues appear "
                f"(Dahl-Sjöstrand 1979 / Kohut 1993)."
            )

    f_lo = _Lambda_real(u0_min, c, f1)
    f_hi = _Lambda_real(u0_max, c, f1)
    if f_lo * f_hi > 0:
        raise ValueError(
            f"No sign change of Λ_real(u_0) in [{u0_min}, {u0_max}] "
            f"for (c={c}, f_1={f1}). "
            f"Λ({u0_min})={f_lo}, Λ({u0_max})={f_hi}. "
            f"Try widening the bracket."
        )

    u0 = brentq(_Lambda_real, u0_min, u0_max, args=(c, f1), xtol=bisect_tol)
    residual = _Lambda_real(u0, c, f1)

    return CaseDispersionResult(
        u0=float(u0),
        c=float(c),
        f1=float(f1),
        Lambda_residual=float(residual),
        bracket_lo=float(u0_min),
        bracket_hi=float(u0_max),
    )


def nu_bar_atalay(c: float, f1: float, u0: float) -> float:
    r"""Compute :math:`\bar\nu = \gamma^{(1)}/\gamma^{(0)}` (Atalay Eq 24).

    For purely imaginary :math:`\nu_0 = i u_0` the moment ratio is
    real-valued. Atalay's Table 1 tabulates :math:`\bar\nu(c)` for
    :math:`f_1 = 0`; this function reproduces those values via the
    integral definition (Eq 25):

    .. math::

       \gamma^{(n)} = \int_0^1 d\mu\,\mu^n\,\gamma(\mu)
                    = \int_0^1 d\mu\,\mu^n \cdot
                       \frac{c \mu / 2}{(1-c)(1 - c f_1)(\nu_0^2 - \mu^2) X(-\mu)} .

    The X-function carries the load-bearing :math:`\mu`-dependence;
    this routine delegates to :mod:`.x_function`.

    For the **isotropic limit** :math:`f_1 = 0`, the bracket
    :math:`(1 - c f_1)` reduces to 1 and :math:`X(-\mu)` is the
    classical isotropic X-function. For Atalay Table 1 row
    :math:`c = 1.30`: :math:`\bar\nu = 0.666526`.

    Parameters
    ----------
    c : float
        Mean number of secondaries per collision.
    f1 : float
        Linearly-anisotropic mean cosine.
    u0 : float
        :math:`u_0 = |\nu_0|`. Must be the root of Atalay Eq 12 for
        consistency.

    Returns
    -------
    nu_bar : float
        :math:`\bar\nu(c, f_1, u_0)`.

    Notes
    -----
    Implemented via the gamma-integral form (Eq 23):

    .. math::

       \gamma(\mu) = \frac{c \mu}{2 (1-c)(1 - c f_1)(\nu_0^2 - \mu^2) X(-\mu)} .

    The :math:`\nu_0^2 - \mu^2 = -u_0^2 - \mu^2` is negative for all
    real :math:`\mu \in (0, 1)`, giving :math:`\gamma(\mu) < 0`.
    Atalay's Table 1 reports positive :math:`\bar\nu` — the ratio
    :math:`\gamma^{(1)}/\gamma^{(0)}` is positive because the negative
    factor cancels.
    """
    from scipy.integrate import quad

    from .half_range import _build_X_minus_nu_cache

    # γ(μ) = c μ / [2 (1-c)(1 - c f_1)(ν_0² - μ²) X(-μ)]
    #      = c μ / [2 (1-c)(1 - c f_1)(-u_0² - μ²) X(-μ)]
    # The factor (1 - c f_1)(1 - c) c/2 is μ-independent; cancels in the ratio.
    # γ^(n) = ∫_0^1 dμ μ^n γ(μ).
    # γ^(0) = ∫ μ / [(-u₀² - μ²) X(-μ)] dμ  (the μ from γ itself).
    # γ^(1) = ∫ μ² / [(-u₀² - μ²) X(-μ)] dμ.
    # ν̄ = γ^(1)/γ^(0).

    # Use the spline-cached X(-μ) to avoid re-evaluating mpmath at each
    # quadrature node (otherwise this routine takes ~30 s per call).
    X_spline = _build_X_minus_nu_cache(c=c, f1=f1, u0=u0, n_grid=64, dps=15)

    def integrand(mu, n):
        X_mmu_real = float(X_spline(mu))
        return mu**n / ((-u0**2 - mu**2) * X_mmu_real)

    g0_value, _ = quad(integrand, 0, 1, args=(1,))
    g1_value, _ = quad(integrand, 0, 1, args=(2,))
    return float(g1_value / g0_value)
