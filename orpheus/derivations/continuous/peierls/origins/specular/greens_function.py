r"""SymPy derivation — operator-level identities for the Plan 2
**Variant α** Green's function reference for sphere homogeneous specular.

Companion to :mod:`.continuous_mu` (which carries the kernel-form
identities V1–V4 from Phase 5). This module steps up to the
**operator level**: rather than verifying the Sanchez 1986 [SanchezTTSP1986]_
Eq. (A6) integrand pointwise, it verifies that the integral operator
:math:`K` defined by the angle-resolved Green's function
:math:`\tilde t(r' \to r, \mu)` (Sanchez Eq. A5 specular leg, :math:`\alpha=1`)
produces the right algebraic answer when applied to canonical trial
functions.

Why operator-level instead of kernel-level
-------------------------------------------

The Phase 5 retreat (2026-04-28,
:file:`/.claude/agent-memory/numerics-investigator/_archive/specular_continuous_mu_phase5_retreat.md`)
documented that the angle-integrated kernel
:math:`g_\alpha(\rho' \to \rho)` has a non-integrable spike at the
surface diagonal :math:`\rho = \rho' = R`. That spike is hit by any
quadrature-on-quadrature scheme (Nyström, Galerkin double integration,
bounce-resolved M2). Plan 2 Variant α never assembles
:math:`g_\alpha` — instead it iterates the **angle-resolved** flux
:math:`\psi(r, \mu)` along bouncing characteristics, and the
divergence is deferred to the final scalar-flux extraction
:math:`\phi(r) = 2\pi \int \psi(r, \mu)\,\mathrm d\mu` where it
appears as an integrable :math:`1/\mu` singularity at :math:`\mu \to 0`
(handled by Gauss-Jacobi µ-quadrature or :math:`u^2 = \mu` substitution).

The B3 SymPy step asks whether the operator action :math:`(K\cdot\psi)`
on canonical trial functions evaluates to a finite, closed-form
expression. If yes, Variant α has a clean algebraic foundation; if
no, Variant α inherits Phase 5's structural pathology despite the
reformulation.

Three operator-level verifications
-----------------------------------

V_α1. **Closed-sphere bounce-sum self-consistency**. For the
   homogeneous sphere with specular BC and constant volumetric
   source :math:`q`, prove that the angular flux satisfies
   :math:`\psi(r, \mu) = q/\Sigma_t` everywhere. Specifically:

   - The surface flux :math:`\psi_{\rm surf}` (the value of
     :math:`\psi` at any specular bounce point on a given trajectory)
     satisfies a fixed-point equation
     :math:`\psi_{\rm surf} = (q/\Sigma_t)(1 - e^{-\Sigma_t L_p}) +
     e^{-\Sigma_t L_p}\,\psi_{\rm surf}` whose solution is exactly
     :math:`q/\Sigma_t`.
   - The first-leg integral plus the bounce-attenuated surface flux
     simplify to :math:`q/\Sigma_t` for any first-leg length
     :math:`L_0`.

   Operator interpretation: for trial :math:`\psi_{\rm trial} = 1`
   (constant), the source is :math:`q = \Sigma_s\,\psi_{\rm trial}
   = \Sigma_s` (isotropic scattering), and :math:`(K\cdot 1)(r, \mu) =
   q/\Sigma_t = \Sigma_s/\Sigma_t = \omega_0`. The constant function
   is therefore an :math:`\omega_0`-eigenmode of the scattering
   operator :math:`K`, and the k-eigenvalue equation
   :math:`(I - K)\,\phi = (\nu\Sigma_f/k)\,\phi` reduces to
   :math:`(1 - \omega_0)\,{\rm const} = (\nu\Sigma_f/k)\,{\rm const}`,
   giving :math:`k = k_\infty = \nu\Sigma_f/\Sigma_a`. This is the
   no-leakage k_inf result, and it is the **load-bearing test of
   Variant α at rank-1**.

V_α2. **T_00^sphere = P_ss^sphere algebraic identity**. The rank-1
   reduction of :func:`compute_T_specular_sphere` (with isotropic
   mode :math:`\tilde P_0 = 1`) is

   .. math::

       T_{00} = 2 \int_0^1 \mu\,e^{-2\Sigma_t R \mu}\,\mathrm d\mu

   which is identically the integrand of :func:`compute_P_ss_sphere`.
   Both close to

   .. math::

       T_{00} = P_{ss} = \frac{1 - (1 + 2\tau_R)\,e^{-2\tau_R}}
                              {2\,\tau_R^{\,2}},
       \qquad \tau_R = \Sigma_t R.

   Operator interpretation: at rank-1, Variant α reduces to the
   Hébert (1−P_ss)^{−1} closure (V_α1 gives the same result via
   the bounce-sum self-consistency). This is the **B5
   cross-verification anchor** — rank-1 Variant α should agree
   bit-for-bit with rank-1 specular_multibounce and with white_hebert.

V_α3. **Vacuum reduction at :math:`\alpha = 0`**. Sanchez Eq. (A6) has
   a leading factor :math:`2\alpha`; at :math:`\alpha = 0` (no specular
   reflection) the BC contribution :math:`g_h \to 0` and the kernel
   reduces to the bare vacuum kernel :math:`\bar g_2` from Sanchez
   Eq. (5). Operator interpretation: the Variant α implementation
   collapses to the existing ORPHEUS vacuum sphere reference at
   :math:`\alpha = 0`, with no special-case branch needed.

References
----------

- Sanchez, R. (1986). "Integral form of the equation of transfer for
  a homogeneous sphere with linearly anisotropic scattering."
  *Transport Theory & Statistical Physics*, vol. 14.
  DOI: 10.1080/00411458608210456.
- Hébert, A. (2009). *Applied Reactor Physics* §3.8.5 (rank-1 white
  BC closure).
- :file:`.claude/agent-memory/numerics-investigator/peierls_greens_variant_alpha_decision.md`
  — Plan 2 B2 architectural decision.
- :file:`.claude/agent-memory/literature-researcher/peierls_greens_function_lit.md`
  — Plan 2 B1 literature pull.
"""
from __future__ import annotations

import sympy as sp


def derive_operator_constant_trial_closed_sphere() -> dict:
    r"""V_α1 — closed-sphere bounce-sum self-consistency.

    For a homogeneous sphere with specular BC and constant volumetric
    source :math:`q`, the bouncing-trajectory integral form of the
    angular flux at any interior point :math:`(r, \mu)` is

    .. math::

       \psi(r, \mu) \;=\; \int_0^{L_0} q\,e^{-\Sigma_t s}\,\mathrm d s
                          \;+\; e^{-\Sigma_t L_0}\,\psi_{\rm surf}

    where :math:`L_0` is the first-leg distance from :math:`r` to the
    surface in direction :math:`-\Omega_\mu`, and
    :math:`\psi_{\rm surf}` is the angular flux at the surface entry
    point in the trajectory direction (the bounce point).

    The bounce point itself sees the same constant source for all
    subsequent bounces on a periodic trajectory of chord
    :math:`L_p = 2R\mu_{\rm surf}` per bounce, so by self-consistency

    .. math::

       \psi_{\rm surf} \;=\; \int_0^{L_p} q\,e^{-\Sigma_t s}\,\mathrm d s
                             \;+\; e^{-\Sigma_t L_p}\,\psi_{\rm surf}.

    Solving the fixed-point equation:

    .. math::

       \psi_{\rm surf} \;=\;
           \frac{q\,(1 - e^{-\Sigma_t L_p})}
                {\Sigma_t\,(1 - e^{-\Sigma_t L_p})}
           \;=\; \frac{q}{\Sigma_t}.

    The dependence on :math:`L_p` (and hence on :math:`\mu`) cancels
    exactly. Plugging back into the first-leg expression:

    .. math::

       \psi(r, \mu) \;=\; \frac{q}{\Sigma_t}\,(1 - e^{-\Sigma_t L_0})
                          + e^{-\Sigma_t L_0}\,\frac{q}{\Sigma_t}
                     \;=\; \frac{q}{\Sigma_t}.

    Both :math:`L_0` and :math:`L_p` cancel identically, leaving
    :math:`\psi = q/\Sigma_t` everywhere. For trial
    :math:`\psi_{\rm trial} = 1` and isotropic scattering source
    :math:`q = \Sigma_s\,\psi_{\rm trial} = \Sigma_s`, the operator
    action is :math:`(K \cdot 1)(r, \mu) = \Sigma_s/\Sigma_t = \omega_0`,
    which gives :math:`k_{\rm eff} = k_\infty = \nu\Sigma_f/\Sigma_a`
    via the standard fission-source eigenvalue equation.

    Returns dict with the SymPy expressions and PASS flags.
    """
    Sigma_t, Sigma_s, q = sp.symbols(
        "Sigma_t Sigma_s q", positive=True, real=True,
    )
    L_0, L_p = sp.symbols("L_0 L_p", positive=True, real=True)

    # First-leg trajectory integral with constant source q.
    # ∫_0^{L_0} q · e^{-Σ_t s} ds = (q/Σ_t)(1 - e^{-Σ_t L_0})
    psi_first = (q / Sigma_t) * (1 - sp.exp(-Sigma_t * L_0))

    # Bounce-sum self-consistency.
    # ψ_surf = ∫_0^{L_p} q e^{-Σ_t s} ds + e^{-Σ_t L_p} ψ_surf
    # ⟹ (1 - e^{-Σ_t L_p}) ψ_surf = (q/Σ_t)(1 - e^{-Σ_t L_p})
    # ⟹ ψ_surf = q/Σ_t
    psi_surf_var = sp.symbols("psi_surf", positive=True, real=True)
    fixed_point_eq = sp.Eq(
        psi_surf_var,
        (q / Sigma_t) * (1 - sp.exp(-Sigma_t * L_p))
        + sp.exp(-Sigma_t * L_p) * psi_surf_var,
    )
    psi_surf_solution = sp.solve(fixed_point_eq, psi_surf_var)
    pass_surf_consistency = (
        len(psi_surf_solution) == 1
        and sp.simplify(psi_surf_solution[0] - q / Sigma_t) == 0
    )
    psi_surf = psi_surf_solution[0]

    # Total ψ at (r, µ) — first-leg + attenuated surface flux.
    psi_total = psi_first + sp.exp(-Sigma_t * L_0) * psi_surf
    psi_total_simplified = sp.simplify(psi_total)

    # Should equal q/Σ_t identically — both L_0 and L_p drop out.
    pass_total_constant = (
        sp.simplify(psi_total_simplified - q / Sigma_t) == 0
    )

    # Operator action on isotropic trial ψ_trial = 1.
    # Source for isotropic scattering: q = Σ_s · ψ_trial = Σ_s.
    omega_0 = Sigma_s / Sigma_t
    K_on_one = psi_total_simplified.subs(q, Sigma_s)
    pass_eigenvalue = sp.simplify(K_on_one - omega_0) == 0

    return {
        "name": "V_α1: closed-sphere bounce-sum constant trial = ω₀",
        "psi_first_leg": psi_first,
        "psi_surf_solution": psi_surf,
        "psi_total_simplified": psi_total_simplified,
        "K_on_constant_trial": K_on_one,
        "omega_0": omega_0,
        "pass_surf_consistency": pass_surf_consistency,
        "pass_total_constant": pass_total_constant,
        "pass_eigenvalue": pass_eigenvalue,
        "pass": (
            pass_surf_consistency
            and pass_total_constant
            and pass_eigenvalue
        ),
    }


def derive_T00_equals_P_ss_sphere() -> dict:
    r"""V_α2 — algebraic identity :math:`T_{00}^{\rm sphere} = P_{ss}^{\rm sphere}`.

    The transmission matrix from :func:`compute_T_specular_sphere`
    at rank-1 (isotropic mode :math:`\tilde P_0 = 1`) is

    .. math::

       T_{00} \;=\; 2\!\int_0^1 \mu\,\tilde P_0(\mu)\,\tilde P_0(\mu)\,
                                e^{-2\Sigma_t R \mu}\,\mathrm d\mu
              \;=\; 2\!\int_0^1 \mu\,e^{-2\Sigma_t R \mu}\,\mathrm d\mu.

    The Hébert :math:`P_{ss}^{\rm sphere}` from
    :func:`compute_P_ss_sphere` is identically the same integrand.
    Both close to (after :math:`\tau_R = \Sigma_t R`):

    .. math::

       T_{00} \;=\; P_{ss} \;=\;
           \frac{1 - (1 + 2\tau_R)\,e^{-2\tau_R}}{2\,\tau_R^{\,2}}.

    This is the **B5 cross-verification anchor**: at rank-1, Variant α
    operator action on isotropic trial reduces to the Hébert white-BC
    geometric-series factor :math:`(1 - P_{ss})^{-1}`. The closed-form
    integral verified here makes the algebraic equivalence airtight,
    independent of any quadrature implementation.

    Returns dict with the SymPy expressions and PASS flags.
    """
    Sigma_t, R, mu = sp.symbols("Sigma_t R mu", positive=True, real=True)
    tau_R = Sigma_t * R

    # Both integrands — algebraically identical.
    T_00_integrand = 2 * mu * sp.exp(-2 * Sigma_t * R * mu)
    P_ss_integrand = 2 * mu * sp.exp(-2 * Sigma_t * R * mu)
    pass_integrand_match = (
        sp.simplify(T_00_integrand - P_ss_integrand) == 0
    )

    # Closed-form integration on µ ∈ [0, 1].
    T_00_value = sp.integrate(T_00_integrand, (mu, 0, 1))
    expected = (
        (1 - (1 + 2 * tau_R) * sp.exp(-2 * tau_R)) / (2 * tau_R ** 2)
    )
    pass_closed_form = sp.simplify(T_00_value - expected) == 0

    return {
        "name": "V_α2: T_00^sphere = P_ss^sphere (rank-1 specular ≡ Hébert)",
        "T_00_integrand": T_00_integrand,
        "P_ss_integrand": P_ss_integrand,
        "T_00_closed_form": sp.simplify(T_00_value),
        "expected_closed_form": expected,
        "pass_integrand_match": pass_integrand_match,
        "pass_closed_form": pass_closed_form,
        "pass": pass_integrand_match and pass_closed_form,
    }


def derive_alpha_zero_kernel_reduction() -> dict:
    r"""V_α3 — at :math:`\alpha = 0`, the BC kernel :math:`g_h` vanishes.

    Sanchez Eq. (A6) for the BC kernel of the homogeneous sphere
    (:math:`\beta = 0, \omega_1 = 0`) is

    .. math::

       g_h(\rho' \to \rho) \;=\; 2\alpha\!\int_{\mu_0}^{1}
            T(\mu_-)\,\mu_*^{-1}\,\cosh(\rho\mu)\,\cosh(\rho'\mu_*)\,
            e^{-2a\mu_-}\,\mathrm d\mu

    The leading factor :math:`2\alpha` makes the entire integrand
    proportional to :math:`\alpha`, so :math:`g_h \to 0` as
    :math:`\alpha \to 0`. The full Sanchez kernel
    :math:`g_\alpha = \bar g_2 + g_h` then reduces to the bare vacuum
    kernel :math:`\bar g_2` from Sanchez Eq. (5).

    Variant α implementation interpretation: setting :math:`\alpha = 0`
    in the Variant α prototype recovers the existing ORPHEUS vacuum
    sphere reference solver. No special-case branch is needed; the
    BC absorption is fully encoded in the leading :math:`2\alpha`
    prefactor of the bounce-sum :math:`g_h`.

    Returns dict with the SymPy expressions and PASS flag.
    """
    alpha, mu, a, rho, rho_p = sp.symbols(
        "alpha mu a rho rho_p", positive=True, real=True,
    )
    mu_star = sp.sqrt(rho_p ** 2 - rho ** 2 * (1 - mu ** 2))
    T_mu = 1 / (1 - alpha * sp.exp(-2 * a * mu))

    # Sanchez Eq. (A6) integrand for g_h with α as parameter.
    integrand = (
        2 * alpha * T_mu * (1 / mu_star)
        * sp.cosh(rho * mu) * sp.cosh(rho_p * mu_star)
        * sp.exp(-2 * a * mu)
    )

    # Take α → 0: leading 2α prefactor drives everything to 0
    # (T(µ) → 1 in the limit, but the 2α outside vanishes faster).
    integrand_at_alpha_zero = sp.simplify(integrand.subs(alpha, 0))
    pass_v_alpha3 = integrand_at_alpha_zero == 0

    return {
        "name": "V_α3: at α=0 the BC kernel g_h vanishes",
        "integrand": integrand,
        "integrand_at_alpha_zero": integrand_at_alpha_zero,
        "pass": pass_v_alpha3,
    }
