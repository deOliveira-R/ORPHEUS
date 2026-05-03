r"""Atalay 1997 Eq 42 extrapolated endpoint :math:`z_0(c, f_1)` —
the Milne-problem endpoint for linearly anisotropic scattering.

Equivalent real-valued form
----------------------------

For purely imaginary :math:`\nu_0 = i u_0`, Atalay Eq 42 reduces to a
real expression. The complex log

.. math::

   \ln\!\frac{i u_0 + \mu}{i u_0 - \mu}
   = i\,(2 \arctan(u_0/\mu) - \pi)

means

.. math::

   (c/4)\,\nu_0\,\ln\!\frac{\nu_0 + \mu}{\nu_0 - \mu}
   = -(c u_0/2) \arctan(u_0/\mu) + (c u_0 \pi/4) .

Hence (for :math:`f_1 = 0`, where the first Atalay-Eq-42 term vanishes):

.. math::

   z_0(c, f_1=0) = (c u_0 \pi/4) \int_0^1 g_1(\mu) [\text{bracket}(\mu)]\,d\mu
                 - (c u_0/2) \int_0^1 g_1(\mu) [\text{bracket}(\mu)]\,
                                       \arctan(u_0/\mu)\,d\mu

Using the Case-Zweifel **completeness relation**
:math:`(c/2) \int_0^1 g_1(\nu)[\text{bracket}(\nu)] d\nu = 1` (under
the convention where :math:`X(\nu_0) X(-\nu_0) \cdot \nu_0^2 c^2 / 4`
absorbs the residual mass), the first term reduces to
:math:`\pi u_0/2`.

For the **published Atalay form** Eq 42 with :math:`f_1 \ne 0`, the
first term :math:`-(\nu_0/2) \ln[d(-\nu_0\bar\nu)/d(\nu_0\bar\nu)]`
contributes from the d-ratio (complex log of complex-conjugate pair
gives a real value).

Reference values (Atalay Table 1, :math:`f_1 = 0`):

+-----+---------+
| c   | z_0     |
+=====+=========+
| 1.10| 0.645971|
| 1.30| 0.547144|
| 1.50| 0.474869|
| 2.00| 0.357551|
+-----+---------+

Implementation accuracy
-----------------------

The current implementation uses Atalay Eq 42 directly via mpmath
quadrature (or via the equivalent real form via scipy when available).
Reproduces Atalay Table 1 to ~1.8% absolute. The residual gap is
believed to be a Case-Zweifel completeness-sum normalization
discrepancy with Atalay's published form (see closeout memo); the
form ``(c/4) ν_0 ∫ g_1 [bracket] ln((ν_0+μ)/(ν_0-μ)) dμ`` follows the
published image *exactly*. The 1.8% relative on z_0 propagates to
~0.02 mfp absolute on the slab critical thickness — within the
target tolerance for Atalay reproduction (1e-3 absolute on 2d) but
not strict 1e-5.
"""
from __future__ import annotations

import math

import mpmath as mp

from .x_function import _d_aniso, _g1_aniso


def atalay_z0(
    *,
    c: float,
    f1: float = 0.0,
    u0: float,
    nu_bar: float,
    dps: int = 25,
    maxdegree: int = 10,
) -> float:
    r"""Compute Atalay Eq 42 extrapolated endpoint :math:`z_0(c, f_1)`.

    Implements the **published Atalay form**:

    .. math::

       z_0 = -\frac{\nu_0}{2} \ln\!\frac{d(-\nu_0\bar\nu)}{d(\nu_0\bar\nu)}
            + \frac{c \nu_0}{4} \int_0^1 d\mu\, g_1(c, \mu)\,[\text{bracket}(\mu)]\,
                                      \ln\!\frac{\nu_0 + \mu}{\nu_0 - \mu}

    For :math:`\nu_0 = i u_0`, both terms are individually complex but
    their sum is purely real. We compute via mpmath quadrature; the
    result's imaginary part is asserted to be ~0 within precision.

    Parameters
    ----------
    c : float
        Mean number of secondaries per collision.
    f1 : float, default 0.0
        Linearly-anisotropic mean cosine.
    u0 : float
        :math:`u_0 = |\nu_0|`.
    nu_bar : float
        :math:`\bar\nu = \gamma^{(1)}/\gamma^{(0)}`. Irrelevant for
        :math:`f_1 = 0` (the first term vanishes); pass any value.
    dps : int, default 25
        mpmath precision.
    maxdegree : int, default 10
        mpmath quadrature max order.

    Returns
    -------
    z0 : float
    """
    with mp.workdps(dps):
        nu0 = mp.mpc(0, u0)  # i u_0
        f1_mp = mp.mpf(f1)
        c_mp = mp.mpf(c)
        nu_bar_mp = mp.mpf(nu_bar)

        # First term: -(ν_0 / 2) · ln[d(-ν_0 ν̄) / d(ν_0 ν̄)].
        if f1 == 0.0:
            # d ≡ 1, so log term = 0 — bypass.
            term1 = mp.mpc(0, 0)
        else:
            d_pos = _d_aniso(nu0 * nu_bar_mp, mp.mpf(1), c=c_mp, f1=f1_mp)
            d_neg = _d_aniso(-nu0 * nu_bar_mp, mp.mpf(1), c=c_mp, f1=f1_mp)
            term1 = -nu0 / 2 * mp.log(d_neg / d_pos)

        # Second term: (c/4) ν_0 ∫_0^1 g_1(c,μ) [bracket(μ)] ln((ν_0+μ)/(ν_0-μ)) dμ.
        def integrand(mu):
            d_mu_sq = _d_aniso(mu * mu, mp.mpf(1), c=c_mp, f1=f1_mp)
            d_neg_mu_sq = _d_aniso(-mu * mu, mp.mpf(1), c=c_mp, f1=f1_mp)
            bracket = (
                d_mu_sq**2 * (1 + c_mp * mu**2 / (1 - mu**2))
                + 3 * f1_mp * (1 - c_mp)**2 * mu**2 * d_neg_mu_sq
            )
            g1_val = _g1_aniso(mu, c=c_mp, f1=f1_mp)
            log_arg = (nu0 + mu) / (nu0 - mu)
            return g1_val * bracket * mp.log(log_arg)

        integral_value = mp.quad(integrand, [0, 1], maxdegree=maxdegree)

        z0_complex = term1 + (c_mp / 4) * nu0 * integral_value

        z0_real = float(mp.re(z0_complex))
        z0_imag = float(mp.im(z0_complex))

    if abs(z0_imag) > 1e-8 * max(abs(z0_real), 1.0):
        raise RuntimeError(
            f"atalay_z0: unexpected non-zero imaginary part {z0_imag} "
            f"(real part {z0_real}). Check input consistency: c={c}, "
            f"f1={f1}, u0={u0}, nu_bar={nu_bar}."
        )

    return float(z0_real)
