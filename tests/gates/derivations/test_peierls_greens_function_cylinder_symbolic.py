r"""Paired symbolic-vs-textbook contract test for the **cylinder**
Variant α Green's function reference (Phase 1 standalone).

Math-origin pattern: the SymPy derivation in
:mod:`orpheus.derivations.continuous.characteristic.origins.specular.greens_function_cylinder`
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

**Retired in step (e2) of the characteristic-reference campaign** (2026-10-10):
the 4 V_α1_cyl rows, which re-stated the sphere's surface fixed point
(``test_peierls_greens_function_symbolic.py``), with their ``derive_*`` function.
Their label ``peierls-greens-cylinder-trajectory`` moved to the kernel's chord
row ``tests/gates/geometry/test_chord.py::test_the_multi_region_segments_are_the_hand_written_table``.
"""
from __future__ import annotations

import pytest
import sympy as sp

from orpheus.derivations.continuous.characteristic.origins.specular import (
    derive_T00_equals_P_ss_cylinder,
    derive_alpha_zero_kernel_reduction_cylinder,
    derive_bounce_period_chord_cylinder,
    derive_homogeneous_limit_reducibility_cylinder_mr,
    derive_piecewise_3d_optical_depth_cylinder_mr,
    derive_two_region_constant_source_consistency_cylinder_mr,
)


# ═══════════════════════════════════════════════════════════════════════
# V_α1_cyl — closed-cylinder bounce-sum self-consistency on constant trial
# ═══════════════════════════════════════════════════════════════════════


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


# ═══════════════════════════════════════════════════════════════════════
# Phase 1b — multi-region cylinder Variant α SymPy identities
# (Branch-1 algebra-of-record backing the Phase 1b numpy solver)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-cylinder-mr-piecewise-tau")
def test_v_alpha1_cyl_mr_homogeneous_reducibility():
    r"""V_α1_cyl_mr — piecewise τ + B accumulator reduces to MG at the
    uniform-material limit.

    SymPy-symbolic ground for Gate 1 of the cylinder MR verification
    plan: with all per-region :math:`\Sigma_{t,k} = \Sigma_t` equal,
    the piecewise sum :math:`\sum_k \Sigma_{t,k} \Delta\ell_k`
    collapses exactly to :math:`\Sigma_t L_{\rm period}`. Likewise the
    piecewise source-line integral telescopes to a single integral
    over :math:`[0, L_{\rm period}]`. Both identities are algebraic
    necessities of the MR machinery; failing either would prevent
    Gate 1's numerical bit-identity from holding.
    """
    result = derive_homogeneous_limit_reducibility_cylinder_mr()
    assert result["pass_tau_reducibility"], (
        f"V_α1_cyl_mr τ reducibility failed: diff = {result['diff_tau']}"
    )
    assert result["pass_B_reducibility"], (
        f"V_α1_cyl_mr B reducibility failed: diff = {result['diff_B']}"
    )
    assert result["pass"], f"V_α1_cyl_mr composite failed: {result}"


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-cylinder-mr-piecewise-tau")
def test_v_alpha1_cyl_mr_piecewise_3d_optical_depth():
    r"""V_α1_cyl_mr.b — piecewise 3D optical depth factors the axial
    Jacobian out of the segment sum.

    :math:`\sum_k \Sigma_{t,k}\,\Delta\ell_{2D,k}/\sqrt{1-\mu_{\rm
    axial}^2}` = :math:`\bigl(\sum_k \Sigma_{t,k}\,\Delta\ell_{2D,k}
    \bigr)/\sqrt{1-\mu_{\rm axial}^2}` is symbolically exact: the
    axial-Jacobian factor :math:`1/\sqrt{1-\mu_{\rm axial}^2}` is
    independent of segment index :math:`k` so it factors out of the
    sum. This identity prevents a class of bugs where the lift is
    applied per-segment differently — they would all give the same
    answer when properly factored, but bug-detectable when only one
    segment receives the lift.
    """
    result = derive_piecewise_3d_optical_depth_cylinder_mr()
    assert result["pass"], (
        f"V_α1_cyl_mr.b 3D-Jacobian factoring failed: diff = "
        f"{result['diff']}"
    )


@pytest.mark.foundation
@pytest.mark.verifies(
    "peierls-greens-cylinder-mr-homogeneous-reduction",
    "peierls-greens-cylinder-mr-bounce-sum-piecewise",
)
def test_v_alpha1_cyl_mr_two_region_constant_source_homogeneous_limit():
    r"""V_α1_cyl_mr.q — 2-region piecewise τ surface fixed-point reduces
    to V_α1_cyl at uniform-material limit.

    The strongest Branch-1 algebraic identity backing Gate 1 + Gate 7
    of the MR verification plan: a 2-region cylinder with
    constant-source :math:`q` has surface-flux closure
    :math:`\psi_{\rm surf} = \alpha B / (1 - \alpha e^{-\tau_{\rm
    period}})` where :math:`\tau_{\rm period} = \Sigma_{t,1}\ell_1
    + \Sigma_{t,2}\ell_2` is the piecewise period optical depth. At
    :math:`\Sigma_{t,1} = \Sigma_{t,2} = \Sigma_t` and :math:`\alpha
    = 1`, SymPy proves :math:`\psi_{\rm surf} = q/\Sigma_t` exactly —
    the V_α1_cyl identity lifts verbatim into MR.
    """
    result = derive_two_region_constant_source_consistency_cylinder_mr()
    assert result["pass"], (
        f"V_α1_cyl_mr.q failed: ψ_surf at homog limit = "
        f"{result['psi_surf_homog']}, expected = {result['expected']}, "
        f"diff = {result['diff']}"
    )
