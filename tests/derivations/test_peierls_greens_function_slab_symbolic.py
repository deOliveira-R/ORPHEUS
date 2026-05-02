r"""Paired symbolic-vs-textbook contract test for the **slab** Variant α
Green's function reference (Phase 3A standalone, symmetric reflective
specular BC).

Math-origin pattern: the SymPy derivation in
:mod:`orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function_slab`
is the source of truth for the operator-level identities of slab
Variant α. Mirrors the sphere/cylinder symbolic test gates with slab-
specific geometry and the **2-bounces-per-period** trajectory structure.

Slab Variant α phase-space and conserved invariants
----------------------------------------------------

Phase space :math:`(x, \mu)` with :math:`x \in [0, L]` (Cartesian
position) and :math:`\mu \in [-1, 1]` (signed cosine wrt outward wall
normal at :math:`x = L`). Specular reflection on a symmetric slab
preserves :math:`|\mu|`. Trajectories alternate between the two parallel
walls; one full period contains TWO surface reflections and traverses

.. math::

   L_{\rm period}^{\rm slab}(\mu) = \frac{2L}{|\mu|}.

This is the load-bearing structural difference from sphere/cylinder
(both 1-bounce-per-period). The per-period reflection product is
:math:`\alpha^2`, while the leading factor in the surface fixed-point
closure remains :math:`\alpha^1` (single reflection at the FIRST
surface arrival).

Three operator-level verifications
-----------------------------------

V_α1_slab. **Closed-slab bounce-sum self-consistency**. Constant
   source :math:`q` produces :math:`\psi(x, \mu) = q/\Sigma_t`
   everywhere. Algebraically identical to V_α1 sphere/cylinder — both
   reduce via the surface fixed-point :math:`\psi_{\rm surf} =
   q/\Sigma_t` independent of period chord length.

V_α2_slab. **T_00^slab = P_ss^slab = 2 E_3(Σ_t·L)**. Two structurally-
   independent derivational paths: per-face transfer-matrix definition
   (µ-domain) vs polar escape-probability integral (θ-domain). SymPy
   chokes on the direct integration (`exp(-τ/µ)` Add-as-integer
   nseries error) but reduces via substitution `u = 1/µ` to the
   canonical 2·E_3 form. Numerical structural independence verified
   via arbitrary-precision mpmath quadrature on both ORIGINAL
   integrands.

V_α3_slab. **Vacuum reduction at :math:`\alpha = 0`**. Surface fixed-
   point closure carries leading :math:`\alpha`, so :math:`\psi_{\rm
   surf} \to 0` at :math:`\alpha = 0`. The :math:`\alpha^2` inside the
   geometric resolvent does not affect the limit (both vanish smoothly).

Predecessor / sibling tests:

- :mod:`.test_peierls_greens_function_symbolic` — sphere V_α1/V_α2/V_α3.
- :mod:`.test_peierls_greens_function_cylinder_symbolic` — cylinder
  V_α1_cyl/V_α2_cyl/V_α3_cyl.

References
----------

- :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 3A
  slab Variant α plan.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.peierls_greens_function.origins.specular import (
    derive_T00_equals_P_ss_slab,
    derive_alpha_zero_kernel_reduction_slab,
    derive_operator_constant_trial_closed_slab,
)


# ═══════════════════════════════════════════════════════════════════════
# V_α1_slab — closed-slab bounce-sum self-consistency on constant trial
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-slab-trajectory")
def test_v_alpha1_slab_surface_fixed_point_solves_to_q_over_sigma_t():
    r"""V_α1_slab.a — surface fixed-point gives :math:`\psi_{\rm surf}
    = q/\Sigma_t`.

    The bounce self-consistency equation for slab
    :math:`\psi_{\rm surf} = (q/\Sigma_t)(1 - e^{-\Sigma_t L_{\rm
    period}}) + e^{-\Sigma_t L_{\rm period}}\,\psi_{\rm surf}` has a
    unique solution independent of the period chord :math:`L_{\rm
    period} = 2L/|\mu|`. The algebra is structurally identical to
    V_α1 sphere/cylinder — only the chord formula differs.

    Verifies the slab trajectory eq.
    :eq:`peierls-greens-slab-trajectory`: the surface fixed-point
    closure cancels the first-leg :math:`L_{\rm first}`-dependence by
    construction.
    """
    result = derive_operator_constant_trial_closed_slab()
    assert result["pass_surf_consistency"], (
        f"V_α1_slab surface fixed-point failed: solution = "
        f"{result['psi_surf_solution']}"
    )


@pytest.mark.foundation
def test_v_alpha1_slab_total_psi_is_independent_of_first_leg():
    r"""V_α1_slab.b — total :math:`\psi(x, \mu) = q/\Sigma_t` everywhere.

    First-leg :math:`(q/\Sigma_t)(1 - e^{-\Sigma_t L_{\rm first}})`
    plus attenuated-surface :math:`e^{-\Sigma_t L_{\rm first}}\,
    \psi_{\rm surf}` cancels the :math:`L_{\rm first}`-dependence
    identically.
    """
    result = derive_operator_constant_trial_closed_slab()
    assert result["pass_total_constant"], (
        f"V_α1_slab total-ψ-constant failed: psi_total = "
        f"{result['psi_total_simplified']}"
    )


@pytest.mark.foundation
def test_v_alpha1_slab_operator_on_constant_gives_omega_0():
    r"""V_α1_slab.c — :math:`(K \cdot 1) = \omega_0 = \Sigma_s/\Sigma_t`.

    For ψ_trial = 1 (isotropic constant) with isotropic-scattering
    source :math:`q = \Sigma_s\,\psi_{\rm trial} = \Sigma_s`, the
    operator action is :math:`q/\Sigma_t = \Sigma_s/\Sigma_t = \omega_0`.
    The k_inf identity follows: :math:`(1 - \omega_0)\,\phi =
    (\nu\Sigma_f/k)\,\phi` ⟹ :math:`k = \nu\Sigma_f/\Sigma_a`.
    """
    result = derive_operator_constant_trial_closed_slab()
    assert result["pass_eigenvalue"], (
        f"V_α1_slab operator-eigenvalue failed: K·1 = "
        f"{result['K_on_constant_trial']}, expected = "
        f"{result['omega_0']}"
    )


@pytest.mark.foundation
def test_v_alpha1_slab_overall_pass():
    """V_α1_slab — composite gate."""
    result = derive_operator_constant_trial_closed_slab()
    assert result["pass"], f"V_α1_slab composite failed: {result}"


# ═══════════════════════════════════════════════════════════════════════
# V_α2_slab — T_00^slab = P_ss^slab = 2 E_3(Σ_t·L)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-slab-T")
def test_v_alpha2_slab_substitution_algebra_holds():
    r"""V_α2_slab.a — the SymPy substitution :math:`u = 1/\mu` on the
    Path A integrand reduces to :math:`2 u^{-3} e^{-\tau u}`, which is
    the integrand of :math:`E_3(\tau)` definitionally.

    This is the symbolic component of V_α2_slab — SymPy verifies the
    algebra of the change-of-variable (Path A's transfer-matrix
    integrand mapped into the E_3 definition domain). Combined with
    the Path B substitution check, the closed form
    :math:`T_{00}^{\rm slab} = 2 E_3(\Sigma_t L)` is established.
    """
    result = derive_T00_equals_P_ss_slab()
    assert result["pass_substitution_algebra"], (
        f"V_α2_slab Path A substitution algebra failed: u-integrand "
        f"differs from expected 2 u^{{-3}} e^{{-τu}}: "
        f"got {result['T_00_u_integrand_after_sub']}, expected "
        f"{result['T_00_u_integrand_expected']}"
    )


@pytest.mark.foundation
def test_v_alpha2_slab_path_b_polar_substitution_to_path_a():
    r"""V_α2_slab.b — Path B's polar :math:`\theta`-integrand
    reduces to Path A's µ-integrand under :math:`\mu = \cos\theta`.

    Symbolically verifies that the polar escape-probability integral
    construction (with :math:`\sin\theta` Jacobian + :math:`\cos\theta`
    cosine weight) maps onto the same µ-domain integrand as the
    transfer-matrix definition. Structural-independence claim: the
    two paths START in different mathematical settings (matrix
    element vs surface escape) and CONVERGE to the same integrand
    after change of variable.
    """
    result = derive_T00_equals_P_ss_slab()
    assert result["pass_substitution_to_T00"], (
        f"V_α2_slab Path B polar substitution failed: substituted "
        f"integrand differs from Path A T_00 integrand"
    )


@pytest.mark.foundation
def test_v_alpha2_slab_closed_form_equals_canonical():
    r"""V_α2_slab.c — Path A closed form equals canonical 2·E_3(Σ_t·L).

    SymPy verifies :math:`T_{00}^{\rm slab} = 2 E_3(\Sigma_t L)`
    symbolically (both sides constructed as ``2 * sp.expint(3, ·)``).
    The canonical literature form is the slab analogue of the sphere
    Hébert closed form :math:`(1 - (1+2\tau_R)e^{-2\tau_R})/(2\tau_R^2)`.
    """
    result = derive_T00_equals_P_ss_slab()
    assert result["pass_T00_equals_canonical"], (
        f"V_α2_slab closed form ≠ canonical 2·E_3: got "
        f"{result['T_00_closed_form']} vs {result['canonical_closed_form']}"
    )


@pytest.mark.foundation
def test_v_alpha2_slab_numerical_path_independence():
    r"""V_α2_slab.d — arbitrary-precision mpmath cross-check at multiple
    :math:`\tau_L` values.

    The structurally-independent NUMERICAL cross-check: mpmath quad
    on the ORIGINAL Path A µ-integrand (:math:`2\mu e^{-\tau/\mu}`)
    and the ORIGINAL Path B θ-integrand (:math:`2\cos\theta\sin\theta
    e^{-\tau/\cos\theta}`) — both integrated DIRECTLY without going
    through the SymPy substitution algebra — must agree with the
    canonical :math:`2 E_3(\tau)` to absolute tolerance 1e-12 across
    six :math:`\tau_L` values spanning thin (0.1) to thick (10).

    This is the load-bearing structural-independence evidence: SymPy
    substitution is one path; mpmath direct quadrature is the
    independent verification.
    """
    result = derive_T00_equals_P_ss_slab()
    assert result["pass_numerical"], (
        f"V_α2_slab numerical cross-check failed at some τ_L: "
        f"max_abs_diff per τ_L = "
        f"{[(d['tau_L'], d['max_abs_diff']) for d in result['numerical_check']]}"
    )


@pytest.mark.foundation
def test_v_alpha2_slab_overall_pass():
    """V_α2_slab — composite gate."""
    result = derive_T00_equals_P_ss_slab()
    assert result["pass"], f"V_α2_slab composite failed: {result}"


# ═══════════════════════════════════════════════════════════════════════
# V_α3_slab — vacuum reduction at α=0
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-slab-architecture")
def test_v_alpha3_slab_psi_surf_vanishes_at_alpha_zero():
    r"""V_α3_slab — surface fixed-point closure :math:`\psi_{\rm surf}
    \to 0` at :math:`\alpha = 0`.

    Leading factor :math:`\alpha` in :math:`\psi_{\rm surf} = \alpha
    B / (1 - \alpha^2 e^{-\Sigma_t L_{\rm period}})` ensures the
    surface-flux contribution vanishes identically at zero specular
    reflection, recovering the vacuum slab kernel via the first-leg
    integral alone. The 2-bounce-per-period :math:`\alpha^2` in the
    denominator does NOT obstruct the limit — both numerator and
    denominator vanish smoothly. The slab Variant α prototype handles
    vacuum BC with no special-case branch.
    """
    result = derive_alpha_zero_kernel_reduction_slab()
    assert result["pass"], (
        f"V_α3_slab failed: psi_surf at α=0 = "
        f"{result['psi_surf_at_alpha_zero']}, limit = "
        f"{result['psi_surf_limit']}"
    )
