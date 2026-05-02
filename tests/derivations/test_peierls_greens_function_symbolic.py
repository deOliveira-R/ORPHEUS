r"""Paired symbolic-vs-textbook contract test for the Plan 2
**Variant α** Green's function reference (sphere homogeneous specular).

Math-origin pattern: the SymPy derivation in
:mod:`orpheus.derivations.continuous.peierls.origins.specular.greens_function`
is the source of truth for the **operator-level** identities of
Variant α. These verifications gate any Variant α prototype
implementation (Plan 2 step B4).

Three operator-level verifications
-----------------------------------

V_α1. **Closed-sphere bounce-sum self-consistency**. For homogeneous
    sphere with specular BC and constant volumetric source :math:`q`,
    the angular flux satisfies :math:`\psi(r, \mu) = q/\Sigma_t`
    everywhere — the surface fixed-point equation has a unique
    solution :math:`\psi_{\rm surf} = q/\Sigma_t` independent of the
    bounce period :math:`L_p`, and the first-leg + attenuated-surface
    sum cancels :math:`L_0` exactly. For trial :math:`\psi_{\rm trial}
    = 1`, the operator action is :math:`(K\cdot 1)(r,\mu) = \omega_0`,
    giving :math:`k_{\rm eff} = k_\infty = \nu\Sigma_f/\Sigma_a`.
    **Load-bearing test of Variant α at rank-1.**

V_α2. **T_00^sphere = P_ss^sphere**. The rank-1 reduction of
    :func:`compute_T_specular_sphere` (isotropic mode :math:`\tilde
    P_0 = 1`) is identically :func:`compute_P_ss_sphere`. Both close
    to :math:`(1 - (1+2\tau_R)e^{-2\tau_R})/(2\tau_R^2)`. **B5
    cross-verification anchor**: rank-1 Variant α ≡ rank-1
    specular_multibounce ≡ Hébert white_hebert.

V_α3. **Vacuum reduction at :math:`\alpha = 0`**. Sanchez Eq. (A6)
    leading factor :math:`2\alpha` makes :math:`g_h \to 0` at
    :math:`\alpha = 0`, so the kernel reduces to the bare vacuum
    kernel :math:`\bar g_2`. No special-case branch needed in the
    Variant α implementation.

Predecessors:

- :mod:`.test_peierls_specular_continuous_mu_symbolic` — kernel-form
  identities V1–V4 for Sanchez (A6) integrand (Phase 5).
- :file:`.claude/agent-memory/numerics-investigator/peierls_greens_variant_alpha_decision.md`
  — Plan 2 B2 architectural decision.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
  DOI: 10.1080/00411458608210456.
- Hébert, A. (2009). *Applied Reactor Physics* §3.8.5.
"""
from __future__ import annotations

import pytest
import sympy as sp

from orpheus.derivations.continuous.peierls.origins.specular import (
    derive_T00_equals_P_ss_sphere,
    derive_alpha_zero_kernel_reduction,
    derive_operator_constant_trial_closed_sphere,
)


# ═══════════════════════════════════════════════════════════════════════
# V_α1 — closed-sphere bounce-sum self-consistency on constant trial
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_v_alpha1_surface_fixed_point_solves_to_q_over_sigma_t():
    r"""V_α1.a — surface fixed-point gives :math:`\psi_{\rm surf} = q/\Sigma_t`.

    The bounce self-consistency equation
    :math:`\psi_{\rm surf} = (q/\Sigma_t)(1 - e^{-\Sigma_t L_p}) +
    e^{-\Sigma_t L_p}\,\psi_{\rm surf}` has a unique solution
    independent of the bounce period :math:`L_p`.
    """
    result = derive_operator_constant_trial_closed_sphere()
    assert result["pass_surf_consistency"], (
        f"V_α1 surface fixed-point failed: solution = "
        f"{result['psi_surf_solution']}"
    )


@pytest.mark.foundation
def test_v_alpha1_total_psi_is_independent_of_first_leg():
    r"""V_α1.b — total :math:`\psi(r,\mu) = q/\Sigma_t` everywhere.

    The first-leg contribution :math:`(q/\Sigma_t)(1 - e^{-\Sigma_t
    L_0})` plus the attenuated-surface contribution :math:`e^{-\Sigma_t
    L_0}\,\psi_{\rm surf}` cancels the :math:`L_0`-dependence
    identically.
    """
    result = derive_operator_constant_trial_closed_sphere()
    assert result["pass_total_constant"], (
        f"V_α1 total-ψ-constant failed: psi_total = "
        f"{result['psi_total_simplified']}"
    )


@pytest.mark.foundation
def test_v_alpha1_operator_on_constant_gives_omega_0():
    r"""V_α1.c — :math:`(K \cdot 1)(r,\mu) = \omega_0 = \Sigma_s/\Sigma_t`.

    For ψ_trial = 1 (isotropic constant) with isotropic-scattering
    source :math:`q = \Sigma_s\,\psi_{\rm trial} = \Sigma_s`, the
    operator action is :math:`q/\Sigma_t = \Sigma_s/\Sigma_t = \omega_0`.
    The k_inf identity follows: :math:`(1 - \omega_0)\,\phi =
    (\nu\Sigma_f/k)\,\phi` ⟹ :math:`k = \nu\Sigma_f/\Sigma_a`.
    """
    result = derive_operator_constant_trial_closed_sphere()
    assert result["pass_eigenvalue"], (
        f"V_α1 operator-eigenvalue failed: K·1 = "
        f"{result['K_on_constant_trial']}, expected = "
        f"{result['omega_0']}"
    )


@pytest.mark.foundation
def test_v_alpha1_overall_pass():
    """V_α1 — composite gate."""
    result = derive_operator_constant_trial_closed_sphere()
    assert result["pass"], f"V_α1 composite failed: {result}"


# ═══════════════════════════════════════════════════════════════════════
# V_α2 — T_00^sphere = P_ss^sphere algebraic identity
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_v_alpha2_T00_matches_hebert_via_matrix_path():
    r"""V_α2.a — Path A: :math:`T_{00}` from the transfer-matrix
    definition reduces to the Hébert closed form.

    SymPy integrates :math:`T_{00} = 2\!\int_0^1 \mu\,P_0(\mu)^2\,
    e^{-2\Sigma_t R\mu}\,\mathrm d\mu` over :math:`\mu \in [0, 1]`
    starting from the matrix-element construction (with
    :math:`P_0(\mu)` kept symbolic). The result must equal
    :math:`(1 - (1 + 2\tau_R)e^{-2\tau_R})/(2\tau_R^2)`.
    """
    result = derive_T00_equals_P_ss_sphere()
    assert result["pass_T00_matches_hebert"], (
        f"V_α2 Path A failed: T_00 closed form = "
        f"{result['T_00_closed_form']}, Hébert = "
        f"{result['hebert_closed_form']}"
    )


@pytest.mark.foundation
def test_v_alpha2_Pss_matches_hebert_via_polar_path():
    r"""V_α2.b — Path B: :math:`P_{ss}` from the escape-probability
    polar integral reduces to the same Hébert closed form.

    SymPy integrates :math:`P_{ss} = 2\!\int_0^{\pi/2}
    \cos\theta'\sin\theta'\,e^{-2\Sigma_t R \cos\theta'}\,\mathrm d\theta'`
    over :math:`\theta' \in [0, \pi/2]` — a different domain with a
    different Jacobian (:math:`\sin\theta'` solid-angle measure) from
    Path A. The independent integration must converge to the same
    Hébert form. This is the structural-independence pillar of V_α2:
    two different mathematical objects (transfer-matrix element vs
    escape probability) reduce to the same closed form via SymPy
    integrating two distinct integrals.
    """
    result = derive_T00_equals_P_ss_sphere()
    assert result["pass_Pss_matches_hebert"], (
        f"V_α2 Path B failed: P_ss closed form = "
        f"{result['P_ss_closed_form']}, Hébert = "
        f"{result['hebert_closed_form']}"
    )


@pytest.mark.foundation
def test_v_alpha2_T00_equals_Pss_closed_forms():
    r"""V_α2.c — closed-form match :math:`T_{00} = P_{ss}` after
    independent SymPy integration.

    Cross-check that Paths A and B produce numerically-equal SymPy
    closed forms (transitively from V_α2.a and V_α2.b, but verified
    here as a direct algebraic identity to catch any SymPy
    simplification asymmetry).
    """
    result = derive_T00_equals_P_ss_sphere()
    assert result["pass_T00_equals_Pss"], (
        f"V_α2 cross-path identity failed: "
        f"T_00 = {result['T_00_closed_form']}, "
        f"P_ss = {result['P_ss_closed_form']}"
    )


@pytest.mark.foundation
def test_v_alpha2_overall_pass():
    """V_α2 — composite gate."""
    result = derive_T00_equals_P_ss_sphere()
    assert result["pass"], f"V_α2 composite failed: {result}"


# ═══════════════════════════════════════════════════════════════════════
# V_α3 — vacuum reduction at α=0
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_v_alpha3_g_h_vanishes_at_alpha_zero():
    r"""V_α3 — Sanchez (A6) BC kernel :math:`g_h \to 0` at :math:`\alpha = 0`.

    The leading :math:`2\alpha` prefactor of Sanchez Eq. (A6) ensures
    that the BC kernel vanishes identically at zero specular
    reflection, recovering the vacuum sphere kernel
    :math:`g_\alpha = \bar g_2`. The Variant α prototype therefore
    inherits the existing ORPHEUS vacuum reference at :math:`\alpha
    = 0` with no special-case logic.
    """
    result = derive_alpha_zero_kernel_reduction()
    assert result["pass"], (
        f"V_α3 failed: integrand at α=0 = "
        f"{result['integrand_at_alpha_zero']}"
    )
