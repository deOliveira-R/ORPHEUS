r"""Paired symbolic-vs-textbook contract test for the **cylinder**
Variant α Green's function reference (Phase 1 standalone).

Math-origin pattern: the SymPy derivation in
:mod:`orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function_cylinder`
is the source of truth for the operator-level identities of cylinder
Variant α. Mirrors the sphere symbolic test gates with cylinder-specific
geometry and the SAME bounce-sum closure algebra.

Cylinder Variant α phase-space and conserved invariants
--------------------------------------------------------

Phase space :math:`(r, \mu_{\rm axial}, \varphi_{\rm az})`. Specular
reflection on an infinite cylinder preserves both :math:`\mu_{\rm axial}`
(cosine to axis) and :math:`b = r\,|\sin\varphi_{\rm az}|` (impact
parameter). The 3D bounce-period chord is

.. math::

   L_{\rm period}(b, \mu_{\rm axial}) =
       \frac{2\sqrt{R^2 - b^2}}{\sqrt{1 - \mu_{\rm axial}^2}}.

This is the load-bearing geometry — the bounce-sum closure is
parametrised by :math:`L_{\rm period}` rather than the sphere's
:math:`2 R \mu_{\rm surf}`.

Three operator-level verifications
-----------------------------------

V_α1_cyl. **Closed-cylinder bounce-sum self-consistency**. Constant
    source :math:`q` produces :math:`\psi(r, \mu_{\rm axial},
    \varphi_{\rm az}) = q/\Sigma_t` everywhere. Algebraically identical
    to V_α1 sphere — both reduce via the surface fixed-point
    :math:`\psi_{\rm surf} = q/\Sigma_t` independent of period chord
    length.

V_α2_cyl. **T_00^cyl = P_ss^cyl**. Rank-1 cylinder Knyazev T-matrix
    integrand identically equals cylinder Hébert :math:`P_{ss}` integrand
    (same :math:`\mathrm{Ki}_3` kernel, same cosα weight, same prefactor).

V_α3_cyl. **Vacuum reduction at :math:`\alpha = 0`**. Surface fixed-point
    closure carries leading factor :math:`\alpha` so :math:`\psi_{\rm
    surf} \to 0` at :math:`\alpha = 0`. No special-case branch needed.

Predecessor / sibling tests:

- :mod:`.test_peierls_greens_function_symbolic` — sphere V_α1/V_α2/V_α3.

References
----------

- :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 1
  cylinder Variant α plan.
"""
from __future__ import annotations

import pytest
import sympy as sp

from orpheus.derivations.continuous.peierls_greens_function.origins.specular import (
    derive_T00_equals_P_ss_cylinder,
    derive_alpha_zero_kernel_reduction_cylinder,
    derive_bounce_period_chord_cylinder,
    derive_operator_constant_trial_closed_cylinder,
)


# ═══════════════════════════════════════════════════════════════════════
# V_α1_cyl — closed-cylinder bounce-sum self-consistency on constant trial
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-cylinder-trajectory")
def test_v_alpha1_cyl_surface_fixed_point_solves_to_q_over_sigma_t():
    r"""V_α1_cyl.a — surface fixed-point gives :math:`\psi_{\rm surf}
    = q/\Sigma_t`.

    The bounce self-consistency equation for cylinder
    :math:`\psi_{\rm surf} = (q/\Sigma_t)(1 - e^{-\Sigma_t L_{\rm period}})
    + e^{-\Sigma_t L_{\rm period}}\,\psi_{\rm surf}` has a unique
    solution independent of the bounce period :math:`L_{\rm period}`.
    The algebra is structurally identical to V_α1 sphere — only the
    chord formula :math:`L_{\rm period}(b, \mu_{\rm axial})` differs.

    Verifies the first-leg trajectory eq.
    :eq:`peierls-greens-cylinder-trajectory`: the surface fixed-point
    closure cancels the first-leg :math:`L_0`-dependence by
    construction.
    """
    result = derive_operator_constant_trial_closed_cylinder()
    assert result["pass_surf_consistency"], (
        f"V_α1_cyl surface fixed-point failed: solution = "
        f"{result['psi_surf_solution']}"
    )


@pytest.mark.foundation
def test_v_alpha1_cyl_total_psi_is_independent_of_first_leg():
    r"""V_α1_cyl.b — total :math:`\psi(r, \mu_{\rm axial},
    \varphi_{\rm az}) = q/\Sigma_t` everywhere.

    First-leg :math:`(q/\Sigma_t)(1 - e^{-\Sigma_t L_0})` plus
    attenuated-surface :math:`e^{-\Sigma_t L_0}\,\psi_{\rm surf}` cancels
    the :math:`L_0`-dependence identically.
    """
    result = derive_operator_constant_trial_closed_cylinder()
    assert result["pass_total_constant"], (
        f"V_α1_cyl total-ψ-constant failed: psi_total = "
        f"{result['psi_total_simplified']}"
    )


@pytest.mark.foundation
def test_v_alpha1_cyl_operator_on_constant_gives_omega_0():
    r"""V_α1_cyl.c — :math:`(K \cdot 1) = \omega_0 = \Sigma_s/\Sigma_t`.

    For ψ_trial = 1 (isotropic constant) with isotropic-scattering
    source :math:`q = \Sigma_s\,\psi_{\rm trial} = \Sigma_s`, the
    operator action is :math:`q/\Sigma_t = \Sigma_s/\Sigma_t = \omega_0`.
    The k_inf identity follows: :math:`(1 - \omega_0)\,\phi =
    (\nu\Sigma_f/k)\,\phi` ⟹ :math:`k = \nu\Sigma_f/\Sigma_a`.
    """
    result = derive_operator_constant_trial_closed_cylinder()
    assert result["pass_eigenvalue"], (
        f"V_α1_cyl operator-eigenvalue failed: K·1 = "
        f"{result['K_on_constant_trial']}, expected = "
        f"{result['omega_0']}"
    )


@pytest.mark.foundation
def test_v_alpha1_cyl_overall_pass():
    """V_α1_cyl — composite gate."""
    result = derive_operator_constant_trial_closed_cylinder()
    assert result["pass"], f"V_α1_cyl composite failed: {result}"


# ═══════════════════════════════════════════════════════════════════════
# V_α1_cyl.geometry — bounce-period 3D chord formula
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies(
    "peierls-greens-cylinder-bounce-period",
    "peierls-greens-cylinder-in-plane-speed",
    "peierls-greens-cylinder-impact-parameter",
)
def test_v_alpha1_cyl_bounce_period_chord_two_derivations_agree():
    r"""V_α1_cyl.geometry — :math:`L_{\rm period} = 2\sqrt{R^2-b^2} /
    \sqrt{1-\mu_{\rm axial}^2}`.

    Two independent derivations (impact-parameter form and
    surface-tangent form) of the cylinder bounce-period 3D chord
    must agree. Without this geometric foundation, the bounce-sum
    closure has the wrong period and V_α1_cyl's algebraic cancellation
    breaks.

    Verifies three coupled cylinder geometry equations:

    - :eq:`peierls-greens-cylinder-impact-parameter` —
      :math:`b = r|\sin\varphi_{\rm az}|` is the conserved impact
      parameter, used in both derivations of :math:`L_{\rm period}`.
    - :eq:`peierls-greens-cylinder-in-plane-speed` —
      :math:`s_{\rm in\!-\!plane} = \sqrt{1-\mu_{\rm axial}^2}` is
      the in-plane velocity fraction; the :math:`1/s_{\rm in\!-\!plane}`
      factor in the 3D chord is what couples axial and in-plane geometry.
    - :eq:`peierls-greens-cylinder-bounce-period` —
      :math:`L_{\rm period}(b, \mu_{\rm axial}) = 2\sqrt{R^2-b^2}/
      \sqrt{1-\mu_{\rm axial}^2}` itself, the load-bearing geometry
      identity.
    """
    result = derive_bounce_period_chord_cylinder()
    assert result["pass"], (
        f"V_α1_cyl bounce-period chord disagrees between derivations: "
        f"v1 = {result['L_period_via_b']}, v2 = "
        f"{result['L_period_via_alpha']}"
    )


# ═══════════════════════════════════════════════════════════════════════
# V_α2_cyl — T_00^cyl = P_ss^cyl algebraic identity
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_v_alpha2_cyl_integrands_are_identical():
    r"""V_α2_cyl.a — :math:`T_{00}^{\rm cyl}` and :math:`P_{ss}^{\rm cyl}`
    have identical integrand.

    Both reduce to :math:`\cos\alpha\,\mathrm{Ki}_3(2\Sigma_t R \cos\alpha)`
    on :math:`\alpha \in [0, \pi/2]` for homogeneous cylinder. Same
    Bickley-Naylor kernel from polar integration with :math:`\sin^2\beta`
    weight; same in-plane chord :math:`2R\cos\alpha`; same prefactor
    :math:`4/\pi`.
    """
    result = derive_T00_equals_P_ss_cylinder()
    assert result["pass_integrand_match"], (
        f"V_α2_cyl integrand-match failed: T_00 - P_ss = "
        f"{sp.simplify(result['T_00_integrand'] - result['P_ss_integrand'])}"
    )


@pytest.mark.foundation
def test_v_alpha2_cyl_full_form_match():
    r"""V_α2_cyl.b — full integrals (with prefactor :math:`4/\pi` and
    bounds :math:`[0, \pi/2]`) are symbolically equal.

    The cylinder analogue of the sphere :math:`T_{00} = P_{ss}` closed
    form. Unlike sphere, the cylinder integral involves
    :math:`\mathrm{Ki}_3` which has no elementary closed form — but
    the symbolic equality of the two integrals is provable.
    """
    result = derive_T00_equals_P_ss_cylinder()
    assert result["pass_full_match"], (
        f"V_α2_cyl full-form match failed: T_00 - P_ss = "
        f"{sp.simplify(result['T_00_full'] - result['P_ss_full'])}"
    )


@pytest.mark.foundation
def test_v_alpha2_cyl_overall_pass():
    """V_α2_cyl — composite gate."""
    result = derive_T00_equals_P_ss_cylinder()
    assert result["pass"], f"V_α2_cyl composite failed: {result}"


# ═══════════════════════════════════════════════════════════════════════
# V_α3_cyl — vacuum reduction at α=0
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
def test_v_alpha3_cyl_psi_surf_vanishes_at_alpha_zero():
    r"""V_α3_cyl — surface fixed-point closure :math:`\psi_{\rm surf}
    \to 0` at :math:`\alpha = 0`.

    Leading factor :math:`\alpha` in :math:`\psi_{\rm surf} = \alpha
    B / (1 - \alpha e^{-\Sigma_t L_{\rm period}})` ensures the
    surface-flux contribution vanishes identically at zero specular
    reflection, recovering the vacuum cylinder kernel via the
    first-leg integral alone. The cylinder Variant α prototype
    therefore handles vacuum BC with no special-case branch.
    """
    result = derive_alpha_zero_kernel_reduction_cylinder()
    assert result["pass"], (
        f"V_α3_cyl failed: psi_surf at α=0 = "
        f"{result['psi_surf_at_alpha_zero']}, limit = "
        f"{result['psi_surf_limit']}"
    )
