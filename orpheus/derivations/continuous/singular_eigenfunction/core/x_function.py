r"""Atalay 1997 Eq 40 X-function — Wiener-Hopf X-function for
linearly-anisotropic scattering.

.. math::

   X(\mu) = \exp\!\Bigg\{ -\frac{c}{2} \int_0^1 d\nu\, g_1(c,\nu)\,
       \Big[d^2(\nu^2)\Big(1 + \frac{c\nu^2}{1-\nu^2}\Big)
            + 3 f_1 (1-c)^2 \nu^2 d(-\nu^2)\Big]\,\ln(\nu - \mu) \Bigg\}

with

* :math:`d(ab) = 1 + 3 f_1 (1-c) ab` (Atalay Eq 10),
* :math:`g_1(c, \nu) = \nu/N(\nu)` (Atalay Eq 35),
* :math:`N(\nu) = \nu \{\lambda^2(\nu) + [\pi c \nu / 2 \cdot d(\nu^2)]^2\}` (Eq 36),
* :math:`\lambda(\nu) = d(\nu^2)(1 - c\nu\,\mathrm{tanh}^{-1}\nu)
       - 3 f_1 (1-c)^2 \nu^2` (Atalay Eq 9; **Note:** the inverse hyperbolic
  tangent ``tanh⁻¹ ν`` for :math:`\nu \in (0, 1)` is real-valued — this
  is the form that appears in Eq 9, distinct from
  :math:`\nu \mathrm{tanh}^{-1}(1/\nu)` in Eq 11).

For the criticality computation we need :math:`X(\mu)` at:

* :math:`\mu = \pm i u_0` (purely imaginary; non-singular).
* :math:`\mu = -\nu` for :math:`\nu \in (0, 1)` — argument is in
  :math:`(-1, 0)`; integrand is non-singular.

Direct evaluation at :math:`\mu \in (0, 1)` would trigger the
logarithmic singularity at :math:`\nu = \mu` and require split-quadrature.

Numerical implementation
------------------------

Two backends:

* **scipy** (default): :func:`scipy.integrate.quad` with split-at-0.5
  for ``mu_arg`` real and not in (0, 1). ~50× faster than mpmath at
  comparable accuracy (1e-6 absolute on X).
* **mpmath** (``backend = "mpmath"``): kept for arbitrary-precision
  diagnostics. Slow (~0.15 s per call at dps=15).

For the criticality computation we use scipy by default; mpmath is
the fallback for cross-verification.
"""
from __future__ import annotations

import math
import warnings

import mpmath as mp
import numpy as np
from scipy.integrate import quad


def _d_aniso(a, b, *, c, f1):
    r"""Atalay Eq 10: :math:`d(ab) = 1 + 3 f_1 (1-c) ab`.

    Works for both float and mpmath inputs.
    """
    return 1.0 + 3.0 * f1 * (1.0 - c) * a * b


def _lambda_aniso(nu, *, c, f1):
    r"""Atalay Eq 9: :math:`\lambda(\nu) = d(\nu^2)(1 - c\nu\,\mathrm{tanh}^{-1}\nu)
    - 3 f_1 (1-c)^2 \nu^2`.

    For :math:`\nu \in (0, 1)`, :math:`\mathrm{tanh}^{-1}\nu` is real;
    this routine handles both float (using ``math.atanh``) and mpmath
    inputs.
    """
    nu_sq = nu * nu
    d2 = _d_aniso(nu_sq, 1.0, c=c, f1=f1)
    if isinstance(nu, (mp.mpf, mp.mpc)):
        atanh_nu = mp.atanh(nu)
    else:
        atanh_nu = math.atanh(nu)
    return d2 * (1 - c * nu * atanh_nu) - 3.0 * f1 * (1.0 - c)**2 * nu_sq


def _N_aniso(nu, *, c, f1):
    r"""Atalay Eq 36: :math:`N(\nu) = \nu \{ \lambda^2(\nu) +
    [\pi c\nu/2 \cdot d(\nu^2)]^2 \}`."""
    lam = _lambda_aniso(nu, c=c, f1=f1)
    d2 = _d_aniso(nu * nu, 1.0, c=c, f1=f1)
    if isinstance(nu, (mp.mpf, mp.mpc)):
        pi_term = mp.pi * c * nu / 2 * d2
    else:
        pi_term = math.pi * c * nu / 2 * d2
    return nu * (lam * lam + pi_term * pi_term)


def _g1_aniso(nu, *, c, f1):
    r"""Atalay Eq 35: :math:`g_1(c, \nu) = \nu / N(\nu) = 1 / [\lambda^2 + (\pi c\nu/2 \cdot d(\nu^2))^2]`."""
    return nu / _N_aniso(nu, c=c, f1=f1)


def _X_integrand_real(nu: float, mu_arg, c: float, f1: float):
    r"""X-function integrand for real or complex :math:`\mu_{\rm arg}`,
    evaluated in float arithmetic (used by scipy backend).

    For complex ``mu_arg`` (e.g. :math:`\pm i u_0`), this returns a
    complex value; scipy.integrate.quad cannot handle complex integrands
    directly, so we split into real-imaginary parts at the call site.
    """
    # Guard against scipy.quad evaluating at the endpoints ν = 0 or ν = 1
    # (the integrand is integrable but the bracket has a 1/(1-ν²) factor).
    if nu >= 1.0 - 1e-300:
        nu = 1.0 - 1e-12
    if nu <= 1e-300:
        nu = 1e-12
    d_pos = _d_aniso(nu * nu, 1.0, c=c, f1=f1)         # d(ν²)
    d_neg = _d_aniso(-nu * nu, 1.0, c=c, f1=f1)        # d(-ν²)
    bracket = (
        d_pos * d_pos * (1 + c * nu * nu / (1 - nu * nu))
        + 3.0 * f1 * (1.0 - c)**2 * nu * nu * d_neg
    )
    g1 = _g1_aniso(nu, c=c, f1=f1)

    # Compute ln(ν - μ_arg).
    z = nu - mu_arg
    if isinstance(z, complex):
        # log of complex: real = (1/2) log|z|², imag = arg(z).
        return g1 * bracket * (
            0.5 * math.log(z.real * z.real + z.imag * z.imag)
            + 1j * math.atan2(z.imag, z.real)
        )
    else:
        return g1 * bracket * math.log(z)


def atalay_X_function(
    mu_arg,
    *,
    c: float,
    f1: float = 0.0,
    u0: float | None = None,
    dps: int = 15,
    maxdegree: int = 8,
    backend: str = "scipy",
) -> complex:
    r"""Compute Atalay Eq 40 X-function at argument :math:`\mu_{\rm arg}`.

    Parameters
    ----------
    mu_arg : complex or float
        Argument of the X-function.
    c : float
        Mean number of secondaries per collision.
    f1 : float, default 0.0
        Linearly-anisotropic mean cosine.
    u0 : float, optional
        :math:`u_0 = |\nu_0|`; not used directly (the integrand does
        not depend on it).
    dps : int, default 15
        mpmath precision (``backend = "mpmath"``).
    maxdegree : int, default 8
        mpmath adaptive quadrature max order (``backend = "mpmath"``).
    backend : str, default "scipy"
        Either ``"scipy"`` (fast; ~1e-6 abs accuracy) or ``"mpmath"``
        (slow but arbitrary precision).

    Returns
    -------
    X(mu_arg) : complex
    """
    if c <= 1.0:
        raise ValueError(f"atalay_X_function requires c > 1, got c={c}")

    if backend == "scipy":
        return _atalay_X_function_scipy(mu_arg, c=c, f1=f1)
    elif backend == "mpmath":
        return _atalay_X_function_mpmath(
            mu_arg, c=c, f1=f1, dps=dps, maxdegree=maxdegree,
        )
    else:
        raise ValueError(f"Unknown backend {backend!r}")


def _atalay_X_function_scipy(mu_arg, *, c: float, f1: float) -> complex:
    """scipy backend: split-at-0.5 quadrature in float arithmetic."""
    is_complex = isinstance(mu_arg, complex) and mu_arg.imag != 0.0
    is_real_in_open_unit = (
        (not isinstance(mu_arg, complex) or mu_arg.imag == 0.0)
        and 0.0 < (mu_arg.real if isinstance(mu_arg, complex) else mu_arg) < 1.0
    )

    if is_complex:
        # Real and imag parts of integrand.
        def integrand_re(nu):
            return _X_integrand_real(nu, mu_arg, c, f1).real

        def integrand_im(nu):
            return _X_integrand_real(nu, mu_arg, c, f1).imag

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=Warning)
            re_part, _ = quad(
                integrand_re, 0.0, 1.0, points=[0.5],
                limit=100, epsabs=1e-12, epsrel=1e-12,
            )
            im_part, _ = quad(
                integrand_im, 0.0, 1.0, points=[0.5],
                limit=100, epsabs=1e-12, epsrel=1e-12,
            )
        log_X = -c / 2.0 * (re_part + 1j * im_part)
        return complex(math.exp(log_X.real) * math.cos(log_X.imag),
                       math.exp(log_X.real) * math.sin(log_X.imag))

    if is_real_in_open_unit:
        # Singular at ν = mu_arg. Split-quadrature.
        mu_real = mu_arg.real if isinstance(mu_arg, complex) else float(mu_arg)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=Warning)
            integral_value, _ = quad(
                lambda nu: _X_integrand_real(nu, mu_real, c, f1).real,
                0.0, 1.0, points=[mu_real],
                limit=100, epsabs=1e-12, epsrel=1e-12,
            )
        return complex(math.exp(-c / 2.0 * integral_value), 0.0)

    # Real argument outside (0, 1) — non-singular, single integral.
    mu_real = float(mu_arg.real) if isinstance(mu_arg, complex) else float(mu_arg)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=Warning)
        integral_value, _ = quad(
            lambda nu: _X_integrand_real(nu, mu_real, c, f1).real,
            0.0, 1.0, points=[0.5],
            limit=100, epsabs=1e-12, epsrel=1e-12,
        )
    return complex(math.exp(-c / 2.0 * integral_value), 0.0)


def _atalay_X_function_mpmath(
    mu_arg, *, c: float, f1: float, dps: int, maxdegree: int,
) -> complex:
    """mpmath backend: arbitrary-precision integration."""
    with mp.workdps(dps):
        if isinstance(mu_arg, complex):
            mu_mp = mp.mpc(mu_arg.real, mu_arg.imag)
        else:
            mu_mp = mp.mpf(float(mu_arg))

        c_mp = mp.mpf(c)
        f1_mp = mp.mpf(f1)

        def integrand_mp(nu):
            d_pos = _d_aniso(nu * nu, mp.mpf(1), c=c_mp, f1=f1_mp)
            d_neg = _d_aniso(-nu * nu, mp.mpf(1), c=c_mp, f1=f1_mp)
            bracket = (
                d_pos * d_pos * (1 + c_mp * nu * nu / (1 - nu * nu))
                + 3 * f1_mp * (1 - c_mp)**2 * nu * nu * d_neg
            )
            g1 = _g1_aniso(nu, c=c_mp, f1=f1_mp)
            return g1 * bracket * mp.log(nu - mu_mp)

        is_real_in_unit = (
            (not isinstance(mu_arg, complex) or mu_arg.imag == 0.0)
            and 0.0 < (mu_arg.real if isinstance(mu_arg, complex) else mu_arg) < 1.0
        )
        if is_real_in_unit:
            integral = mp.quad(integrand_mp, [0, mu_mp, 1], maxdegree=maxdegree)
        else:
            integral = mp.quad(integrand_mp, [0, 1], maxdegree=maxdegree)

        log_X = -c_mp / 2 * integral
        X = mp.exp(log_X)

    return complex(X)
