r"""Atalay 1997 half-range integral moments :math:`K_j(c, R)` (slab,
Eq 38) and :math:`L_j(c, R)` (sphere, Eq 51), and the T-functions
:math:`T(R, \mu)` (Eq 33) and :math:`T_1(R, \mu)` (Eq 50).

These are the load-bearing kernels of the criticality conditions
Eq 39 / 46 (slab) and Eq 52 / 54 (sphere).

T-functions (overflow-safe form)
--------------------------------

* **Slab** (Atalay Eq 33): :math:`T(R, \mu) = (R - e^{-2d/\mu})/(R e^{-2d/\mu} - 1)`.
  Limits: :math:`T(R=0, \mu) = e^{-2d/\mu}` (vacuum), :math:`T(R=1, \mu) = -1`.
* **Sphere** (Atalay Eq 50): :math:`T_1(R, \mu) = (R + e^{-2d/\mu})/(R e^{-2d/\mu} + 1)`.
  Limits: :math:`T_1(R=0, \mu) = e^{-2d/\mu}`, :math:`T_1(R=1, \mu) = +1`.

K_j and L_j moments
-------------------

* **Slab** (Atalay Eq 38):

  .. math::

     K_j(c, R, d) = \int_0^1 d\nu\,\nu^j\,
                 \frac{g_1(c, \nu)}{\bar\nu}\,
                 \big(c\nu/2\big) X^2(-\nu)\,\Lambda(\infty)\,d(-\nu\bar\nu)\,T(R, \nu)

  where :math:`\Lambda(\infty) = (1-c)(1 - c f_1)` (Atalay Eq 37).

* **Sphere** (Atalay Eq 51): same as slab with :math:`T \to T_1`.

Performance
-----------

X(-ν) is precomputed on a Gauss-Legendre grid at moderate ``dps``
(15) and shared across all K_j evaluations for fixed (c, f_1, u_0).
This avoids the nested-quadrature blow-up that would result from
evaluating X(-ν) inside each K_j quadrature node at ``dps = 30``.
Atalay's published tabulation (6 digits) requires ~1e-6 absolute on
K_j; ``dps = 15`` on X is more than sufficient for that.
"""
from __future__ import annotations

import math
from typing import Callable

import mpmath as mp
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline


def T_function(R: float, mu: float, d_thick: float) -> float:
    r"""Atalay Eq 33 slab T-function (overflow-safe form):
    :math:`T(R, \mu) = (R - e^{-2d/\mu})/(R e^{-2d/\mu} - 1)`.
    """
    expmf = math.exp(-2.0 * d_thick / mu)
    num = R - expmf
    den = R * expmf - 1.0
    return num / den


def T1_function(R: float, mu: float, d_thick: float) -> float:
    r"""Atalay Eq 50 sphere :math:`T_1`-function (overflow-safe form):
    :math:`T_1(R, \mu) = (R + e^{-2d/\mu})/(R e^{-2d/\mu} + 1)`.
    """
    expmf = math.exp(-2.0 * d_thick / mu)
    num = R + expmf
    den = R * expmf + 1.0
    return num / den


# ─── Cached X-function evaluations on Gauss-Legendre grid ───


_X_CACHE_KEY: tuple | None = None
_X_CACHE: tuple[np.ndarray, np.ndarray, CubicSpline] | None = None


def _build_X_minus_nu_cache(
    c: float, f1: float, u0: float, *, n_grid: int = 64, dps: int = 15,
) -> CubicSpline:
    r"""Precompute :math:`X(-\nu)` on a fine grid in :math:`\nu \in (0, 1)`
    and return a CubicSpline interpolant.

    The interpolant is reused for all K_j / L_j evaluations at fixed
    ``(c, f_1, u_0)``. ``n_grid = 64`` and ``dps = 15`` give ~1e-6
    accuracy on K_j, more than sufficient for matching Atalay's
    6-digit tabulation.
    """
    global _X_CACHE_KEY, _X_CACHE
    key = (c, f1, u0, n_grid, dps)
    if _X_CACHE_KEY == key and _X_CACHE is not None:
        return _X_CACHE[2]

    from .x_function import atalay_X_function

    # Use a fine Chebyshev-like grid clustered toward 1 (where X varies most).
    # Map uniform t ∈ [0, 1] → ν via ν = (1 - cos(π t)) / 2 (clusters near 0 and 1).
    t = np.linspace(0.0, 1.0, n_grid)
    nus = 0.5 * (1.0 - np.cos(math.pi * t))
    # Avoid endpoints (X has integrable but slow-decaying behaviour at 0, 1).
    nus = nus[1:-1]
    Xs = np.empty(len(nus))
    for i, nu in enumerate(nus):
        x_val = atalay_X_function(-float(nu), c=c, f1=f1, u0=u0, dps=dps, maxdegree=10)
        Xs[i] = x_val.real
    spline = CubicSpline(nus, Xs, extrapolate=True)
    _X_CACHE_KEY = key
    _X_CACHE = (nus, Xs, spline)
    return spline


def _atalay_K_or_L_moment_value(
    *,
    c: float,
    f1: float,
    u0: float,
    nu_bar: float,
    R: float,
    d_thick: float,
    j: int,
    geometry_sign: int,
    n_grid: int = 64,
    x_dps: int = 15,
    quad_limit: int = 80,
    quad_epsabs: float = 1e-10,
    quad_epsrel: float = 1e-10,
) -> float:
    r"""Compute :math:`K_j` (slab) or :math:`L_j` (sphere) by precomputing
    :math:`X(-\nu)` on a grid and using scipy.integrate.quad over
    :math:`\nu \in (0, 1)`.
    """
    from .x_function import _d_aniso, _g1_aniso, _N_aniso, _lambda_aniso

    X_spline = _build_X_minus_nu_cache(c=c, f1=f1, u0=u0, n_grid=n_grid, dps=x_dps)

    Lambda_inf = (1.0 - c) * (1.0 - c * f1)
    T_fn: Callable[[float], float] = (
        (lambda nu: T_function(R, nu, d_thick))
        if geometry_sign == 1
        else (lambda nu: T1_function(R, nu, d_thick))
    )

    def integrand(nu: float) -> float:
        # X(-ν) from spline.
        X_val = float(X_spline(nu))
        X2 = X_val * X_val

        # d(-ν ν̄) = 1 + 3 f_1 (1-c) (-ν ν̄).
        d_neg = 1.0 + 3.0 * f1 * (1.0 - c) * (-nu * nu_bar)

        # g_1(c, ν) = 1 / [λ²(ν) + (πcν/2 d(ν²))²]
        # Compute λ and d(ν²) directly in numpy.
        if 0.0 < nu < 1.0:
            atanh_nu = 0.5 * math.log((1.0 + nu) / (1.0 - nu))
        else:
            atanh_nu = math.atanh(nu)
        d_nu_sq = 1.0 + 3.0 * f1 * (1.0 - c) * nu * nu  # d(ν²)
        lam = d_nu_sq * (1.0 - c * nu * atanh_nu) - 3.0 * f1 * (1.0 - c)**2 * nu * nu
        pi_term = math.pi * c * nu / 2.0 * d_nu_sq
        g1 = 1.0 / (lam * lam + pi_term * pi_term)

        T_val = T_fn(nu)
        return (
            nu**j * (g1 / nu_bar) * (c * nu / 2.0) * X2 * Lambda_inf * d_neg * T_val
        )

    result, _err = quad(
        integrand, 0.0, 1.0,
        limit=quad_limit, epsabs=quad_epsabs, epsrel=quad_epsrel,
    )
    return float(result)


def atalay_K_moments(
    *,
    c: float,
    f1: float,
    u0: float,
    nu_bar: float,
    R: float,
    d_thick: float,
    j_max: int = 2,
    n_grid: int = 64,
    x_dps: int = 15,
    dps: int | None = None,        # accepted for API compat; ignored
    maxdegree: int | None = None,  # accepted for API compat; ignored
) -> dict[int, float]:
    r"""Compute Atalay Eq 38 :math:`K_j(c, R, d)` for slab criticality.

    Parameters
    ----------
    c, f1 : float
        Material parameters.
    u0 : float
        :math:`u_0 = |\nu_0|`.
    nu_bar : float
        :math:`\bar\nu = \gamma^{(1)}/\gamma^{(0)}`.
    R : float
        Reflection coefficient (R = 0 vacuum, R = 1 perfect reflector).
    d_thick : float
        Slab half-thickness :math:`d` in mfp.
    j_max : int, default 2
        Maximum moment order. Atalay Eq 39 needs :math:`K_0, K_1, K_2`.
    n_grid : int, default 64
        Grid points for X(-ν) precomputation.
    x_dps : int, default 15
        mpmath precision for X(-ν) precomputation.
    dps, maxdegree : optional
        Accepted for API compatibility; ignored by the new spline-cached
        implementation.

    Returns
    -------
    Mapping from ``j`` (int) to :math:`K_j` (float).
    """
    return {
        j: _atalay_K_or_L_moment_value(
            c=c, f1=f1, u0=u0, nu_bar=nu_bar, R=R, d_thick=d_thick,
            j=j, geometry_sign=+1, n_grid=n_grid, x_dps=x_dps,
        )
        for j in range(j_max + 1)
    }


def atalay_L_moments(
    *,
    c: float,
    f1: float,
    u0: float,
    nu_bar: float,
    R: float,
    d_thick: float,
    j_max: int = 2,
    n_grid: int = 64,
    x_dps: int = 15,
    dps: int | None = None,        # accepted for API compat; ignored
    maxdegree: int | None = None,  # accepted for API compat; ignored
) -> dict[int, float]:
    r"""Compute Atalay Eq 51 :math:`L_j(c, R, d)` for sphere criticality.

    Same signature and structure as :func:`atalay_K_moments` but uses
    :math:`T_1` (Atalay Eq 50) in place of :math:`T` (Eq 33).
    """
    return {
        j: _atalay_K_or_L_moment_value(
            c=c, f1=f1, u0=u0, nu_bar=nu_bar, R=R, d_thick=d_thick,
            j=j, geometry_sign=-1, n_grid=n_grid, x_dps=x_dps,
        )
        for j in range(j_max + 1)
    }


def clear_X_cache() -> None:
    """Clear the X(-ν) cache. Useful for tests that change (c, f_1, u_0)."""
    global _X_CACHE_KEY, _X_CACHE
    _X_CACHE_KEY = None
    _X_CACHE = None
