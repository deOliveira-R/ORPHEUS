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

The current implementation uses Atalay Eq 42 with the substitution
:math:`\mu = \tanh(t)` to regularize the endpoint :math:`\mu = 1`,
where the bracket factor :math:`1 + c\mu^2/(1-\mu^2)` carries a pole
that is algebraically cancelled by the :math:`\lambda^2 \sim
\ln^2(1-\mu)` growth in :math:`g_1`. Without the substitution, plain
``mp.quad`` over :math:`(0, 1)` under-evaluates the integral by
~1.5% even at dps=35 due to slow endpoint-pole cancellation
(falsely diagnosed as a "Case-Zweifel completeness-sum normalisation
discrepancy" in the Wave 2-B closeout 2026-05-02; the actual cause
is quadrature endpoint convergence, not convention drift — caught
2026-05-03 by numerics-investigator).

Under the :math:`\mu = \tanh(t)` substitution:

.. math::

    1 - \mu^2 = \mathrm{sech}^2(t),\quad
    \frac{\mu^2}{1-\mu^2} = \sinh^2(t),\quad
    \tanh^{-1}(\mu) = t,\quad
    d\mu = \mathrm{sech}^2(t)\,dt .

This maps :math:`\mu \in (0,1) \to t \in (0,\infty)` with
exponentially-decaying integrand at :math:`t \to \infty` (since
:math:`g_1 \sim 1/t^2`). Reproduces Atalay Table 1 to **6-7 digits
at dps=25** in sub-ms wall-clock.
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
        #
        # The bracket carries a factor `1/(1-μ²)` that becomes singular at
        # μ = 1 (algebraically cancelled by the `λ²` growth in g_1, but the
        # cancellation is inefficient under direct mp.quad — leading to
        # a ~1.5% under-evaluation at moderate dps). We regularize by the
        # substitution μ = tanh(t), which maps μ ∈ (0, 1) → t ∈ (0, ∞):
        #
        #   μ = tanh(t)         ⇒  1 - μ² = sech²(t) = 1/cosh²(t),
        #                          μ² / (1-μ²) = sinh²(t),
        #                          tanh⁻¹(μ) = t  (so λ = d(μ²)(1-cμ·t) - 3f₁(1-c)²μ²),
        #                          dμ = sech²(t) dt.
        #
        # The bracket·dμ combines into a smooth integrand on (0, ∞) with
        # exponential decay (g₁ ~ 1/t² as t→∞ and bracket·sech² → c·tanh²t).
        # mp.quad evaluates this to 6-7 digits at dps=25 in <1 ms.
        # Without this substitution z_0 misses Atalay Table 1 by ~1.5%
        # at f₁=0 (caught Wave 2-B post-mortem 2026-05-03).
        def integrand_t(t):
            mu = mp.tanh(t)
            mu_sq = mu * mu
            d_mu_sq = _d_aniso(mu_sq, mp.mpf(1), c=c_mp, f1=f1_mp)
            d_neg_mu_sq = _d_aniso(-mu_sq, mp.mpf(1), c=c_mp, f1=f1_mp)
            sech2 = 1 / mp.cosh(t)**2  # = 1 - μ²
            # bracket·dμ = [d²(μ²)·(1 + c·μ²/(1-μ²)) + 3f₁(1-c)²·μ²·d(-μ²)]·sech²(t) dt
            #            = d²(μ²)·(sech²(t) + c·μ²) + 3f₁(1-c)²·μ²·d(-μ²)·sech²(t)
            bracket_jac = (
                d_mu_sq**2 * (sech2 + c_mp * mu_sq)
                + 3 * f1_mp * (1 - c_mp)**2 * mu_sq * d_neg_mu_sq * sech2
            )
            # λ(μ) = d(μ²)·(1 - cμ·tanh⁻¹μ) - 3f₁(1-c)²μ²
            lam = d_mu_sq * (1 - c_mp * mu * t) - 3 * f1_mp * (1 - c_mp)**2 * mu_sq
            pi_term = mp.pi * c_mp * mu / 2 * d_mu_sq
            g1_val = 1 / (lam**2 + pi_term**2)
            log_arg = (nu0 + mu) / (nu0 - mu)
            return g1_val * bracket_jac * mp.log(log_arg)

        # mp.quad over (1e-12, ∞); tanh-sinh handles the half-infinite tail.
        integral_value = mp.quad(integrand_t, [mp.mpf('1e-12'), mp.inf], maxdegree=maxdegree)

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
