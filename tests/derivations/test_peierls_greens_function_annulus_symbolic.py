r"""Phase-3C-2 symbolic-vs-textbook contract tests for the **annulus**
(hollow cylinder) Variant α Green's function reference (rank-2
boundary-to-boundary scattering resolvent on through-rays + rank-1
outer-only closure, independent :math:`\alpha_{\rm in}, \alpha_{\rm out}
\in [0, 1]`, cylinder 3D angular phase-space).

Math-origin pattern: the SymPy derivation in
:mod:`orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function_annulus`
is the source of truth for the operator-level identities of annulus
Variant α. Mirrors the Phase-3C-1 hollow sphere symbolic test gates with
the cylinder 3D angular-correction lift (the
:math:`1/\sqrt{1 - \mu_{\rm axial}^2}` axial factor on the shell chord).

Three operator-level verifications + 1 chord-algebra auxiliary
---------------------------------------------------------------

V_α1_annulus. **Closed-shell bounce-sum self-consistency at corner**
   :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1`. With both surfaces
   fully reflective and constant volumetric source :math:`q`, BOTH
   the outer-only branch (rank-1, :math:`b > R_{\rm in}`) AND the
   through-ray branch (rank-2, :math:`b \le R_{\rm in}`) independently
   produce :math:`\psi(r, \mu_{\rm axial}, \varphi_{\rm az}) =
   q/\Sigma_t`. Composability of the impact-parameter partition with
   the rank-2 closure on the cylinder's 3D angular phase-space.

V_α2_annulus. **Rank-2 resolvent algebra**. Direct SymPy matrix
   inversion verifies the canonical 2x2 form, the determinant
   :math:`1 - \alpha_{\rm in}\,\alpha_{\rm out}\,e^{-2\tau_{\rm step}}`,
   and four reductions: symmetric, vacuum inner (cavity absorber),
   vacuum outer (reflective cavity), vacuum-vacuum :math:`T = I`.

V_α3_annulus. **Vacuum reduction** at :math:`\alpha_{\rm in} =
   \alpha_{\rm out} = 0`. Closure vanishes; bare first-leg kernel.

V_α2_annulus.aux. **3D chord scaling identity**:
   :math:`\tau_{\rm step}^{\rm annulus}(b, \mu_{\rm axial}) =
   \tau_{\rm step}^{\rm hollow\,sph}(b) / \sqrt{1 - \mu_{\rm axial}^2}`.
   The load-bearing cylinder-3D-axial-correction algebra that lifts
   the hollow-sphere shell chord meaning to the annulus.

References
----------

- :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 3C-2
  annulus plan.
- :mod:`.test_peierls_greens_function_hollow_sphere_symbolic` —
  Phase-3C-1 hollow sphere V_α (the rank-2 + impact-parameter template
  the annulus closure inherits at the operator-symbol level).
- :mod:`.test_peierls_greens_function_cylinder_symbolic` — Phase-1
  cylinder V_α (the cylinder 3D angular framework + outer-only rank-1
  template).
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.peierls_greens_function.origins.specular import (
    derive_3d_chord_scaling_annulus,
    derive_alpha_zero_kernel_reduction_annulus,
    derive_operator_constant_trial_closed_annulus,
    derive_rank2_resolvent_annulus,
)


# ═══════════════════════════════════════════════════════════════════════
# V_α1_annulus — closed-shell bounce-sum self-consistency
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-annulus-architecture")
def test_v_alpha1_annulus_outer_only_closed_BC_equals_q_over_sigma_t():
    r"""V_α1_annulus.a — outer-only branch (:math:`b > R_{\rm in}`) at
    :math:`\alpha_{\rm out} = 1` gives :math:`\psi_{\rm surf} =
    q/\Sigma_t` on a constant source.

    Identical to V_α1 cylinder — outer-only rays in the annulus are
    structurally identical to solid-cylinder rays at the same outer
    radius and impact parameter (with the axial-correction lift built
    into the period chord). The inner cavity is "invisible" to these
    trajectories.
    """
    result = derive_operator_constant_trial_closed_annulus()
    assert result["pass_outer_closed_equals_q_over_sigma_t"], (
        f"V_α1_annulus outer-only at α_out=1 failed: "
        f"got {result['psi_surf_outer_closed']}"
    )


@pytest.mark.foundation
def test_v_alpha1_annulus_psi_in_out_at_closed_BC_equals_q_over_sigma_t():
    r"""V_α1_annulus.b — through-ray branch (:math:`b \le R_{\rm in}`)
    at :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1` gives
    :math:`\psi_{\rm in}^{\rm out} = q/\Sigma_t`.

    The rank-2 closure on a constant source with both surfaces fully
    reflective collapses (via :math:`(1 - e^{-\tau})(1 + e^{-\tau}) =
    1 - e^{-2\tau}`) to the chord-length-independent
    :math:`q/\Sigma_t`. The 3D-lifted chord meaning :math:`\tau_{\rm
    step}^{\rm annulus} = \tau_{\rm step}^{\rm hollow\,sph} /
    \sqrt{1 - \mu_{\rm axial}^2}` does not change the algebraic
    identity — the closure is invariant to the chord meaning.
    """
    result = derive_operator_constant_trial_closed_annulus()
    assert result["pass_psi_in_out_closed_equals_q_over_sigma_t"], (
        f"V_α1_annulus ψ_in^out at closed BC failed: "
        f"got {result['psi_in_out_closed']}"
    )


@pytest.mark.foundation
def test_v_alpha1_annulus_psi_out_in_at_closed_BC_equals_q_over_sigma_t():
    r"""V_α1_annulus.c — through-ray :math:`\psi_{\rm out}^{\rm in} =
    q/\Sigma_t` at closed BC.

    Symmetric counterpart of the :math:`\psi_{\rm in}^{\rm out}`
    identity — both surface flux components must independently equal
    :math:`q/\Sigma_t`.
    """
    result = derive_operator_constant_trial_closed_annulus()
    assert result["pass_psi_out_in_closed_equals_q_over_sigma_t"], (
        f"V_α1_annulus ψ_out^in at closed BC failed: "
        f"got {result['psi_out_in_closed']}"
    )


@pytest.mark.foundation
def test_v_alpha1_annulus_operator_outer_branch_gives_omega_0():
    r"""V_α1_annulus.d — operator action on isotropic trial gives
    :math:`(K \cdot 1)_{\rm outer} = \omega_0 = \Sigma_s/\Sigma_t` on
    the outer-only subset.
    """
    result = derive_operator_constant_trial_closed_annulus()
    assert result["pass_eigenvalue_outer"], (
        f"V_α1_annulus K·1 (outer) failed: got "
        f"{result['K_outer_on_constant_trial']}, expected "
        f"{result['omega_0']}"
    )


@pytest.mark.foundation
def test_v_alpha1_annulus_operator_through_branch_gives_omega_0():
    r"""V_α1_annulus.e — operator action on isotropic trial gives
    :math:`(K \cdot 1)_{\rm through} = \omega_0 = \Sigma_s/\Sigma_t`
    on the through-ray subset.

    Combined with the outer-branch identity, this proves the closed-
    annulus eigenvalue is :math:`k_{\rm eff} = k_\infty =
    \nu\Sigma_f/\Sigma_a` IDENTICALLY across BOTH phase-space subsets
    on the cylinder's 3D angular phase-space. Failure here means the
    impact-parameter partition does not compose cleanly with the
    rank-2 closure when the cylinder axial correction is folded in —
    a structural bug that would HALT Phase 3C-2.
    """
    result = derive_operator_constant_trial_closed_annulus()
    assert result["pass_eigenvalue_inner"], (
        f"V_α1_annulus K·1 (through-ray) failed: got "
        f"{result['K_inner_on_constant_trial']}, expected "
        f"{result['omega_0']}"
    )


@pytest.mark.foundation
def test_v_alpha1_annulus_leaky_corner_strictly_below_q_over_sigma_t():
    r"""V_α1_annulus.f — at :math:`(\alpha_{\rm in}, \alpha_{\rm out})
    = (1/2, 1)` the through-ray surface flux is strictly less than
    :math:`q/\Sigma_t`.

    Leakage at the partial-reflection inner surface keeps the through-
    ray surface flux below :math:`q/\Sigma_t`.
    """
    result = derive_operator_constant_trial_closed_annulus()
    assert result["pass_leaky_below"], (
        f"V_α1_annulus leaky-mode below q/Σ_t failed: "
        f"sample = {result['psi_sample_leaky']}"
    )


@pytest.mark.foundation
def test_v_alpha1_annulus_overall_pass():
    """V_α1_annulus — composite gate."""
    result = derive_operator_constant_trial_closed_annulus()
    assert result["pass"], f"V_α1_annulus composite failed: {result}"


# ═══════════════════════════════════════════════════════════════════════
# V_α2_annulus — rank-2 resolvent T = (I - S)^{-1}
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-annulus-through-rank2")
def test_v_alpha2_annulus_determinant_canonical_form():
    r"""V_α2_annulus.a — :math:`\det(I - S) = 1 - \alpha_{\rm in}\,
    \alpha_{\rm out}\,e^{-2\tau_{\rm step}}` via SymPy direct matrix
    determinant.

    Same algebraic structure as the hollow-sphere / slab-asymmetric
    determinant; only the meaning of :math:`\tau_{\rm step}` changes
    (cylinder 3D-lifted shell chord vs spherical shell chord vs slab
    chord). The singular locus :math:`\alpha_{\rm in}\,\alpha_{\rm out}\,
    e^{-2\tau_{\rm step}} = 1` is reachable only at :math:`\alpha_{\rm
    in} = \alpha_{\rm out} = 1` and :math:`\tau_{\rm step} = 0`
    (tangent rays at :math:`b = R_{\rm in}` ∩ axial-grazing
    :math:`\mu_{\rm axial} \to \pm 1` would ALSO drive :math:`\tau_{\rm
    step} \to \infty` via the divergent :math:`1/s_{\rm ip}` factor —
    the OPPOSITE corner from the singular locus, so the singular
    locus on the annulus is the same isolated point as the hollow
    sphere).
    """
    result = derive_rank2_resolvent_annulus()
    assert result["pass_det"], (
        f"V_α2_annulus determinant failed: got {result['det_M']}, "
        f"expected {result['det_canonical']}"
    )


@pytest.mark.foundation
def test_v_alpha2_annulus_canonical_T_form():
    r"""V_α2_annulus.b — the rank-2 resolvent matches the canonical
    closed form
    :math:`T = (1/\det)\,\bigl[[1, \alpha_{\rm in} e^{-\tau_{\rm step}}],
    [\alpha_{\rm out} e^{-\tau_{\rm step}}, 1]\bigr]`.

    SymPy's ``Matrix.inv()`` produces the inverse, and each entry is
    matched against the canonical form via ``sp.simplify``.
    """
    result = derive_rank2_resolvent_annulus()
    assert result["pass_T_form"], (
        f"V_α2_annulus canonical T form failed: "
        f"T = {result['T_resolvent']}"
    )


@pytest.mark.foundation
def test_v_alpha2_annulus_symmetric_reduction():
    r"""V_α2_annulus.c — at :math:`\alpha_{\rm in} = \alpha_{\rm out}
    = \alpha`, :math:`T_{11} + T_{12} = 1/(1 - \alpha\,e^{-\tau_{
    \rm step}})`.

    The rank-2 → rank-1 collapse on a constant source with symmetric
    BC. Same algebraic identity as V_α2_hollow_sph.c.
    """
    result = derive_rank2_resolvent_annulus()
    assert result["pass_symmetric_simplification"], (
        f"V_α2_annulus symmetric reduction failed."
    )


@pytest.mark.foundation
def test_v_alpha2_annulus_alpha_in_zero_reduction():
    r"""V_α2_annulus.d — at :math:`\alpha_{\rm in} = 0` (cavity
    absorber), the determinant is 1 and only :math:`T_{21} = \alpha_{
    \rm out}\,e^{-\tau_{\rm step}}` is nonzero off-diagonal.

    Physical meaning: with the inner surface absorbing perfectly,
    through-rays that reach the inner cavity are lost. This is the
    annulus analog of the hollow-sphere cavity-absorber case.
    """
    result = derive_rank2_resolvent_annulus()
    assert result["pass_alpha_in_zero"], (
        f"V_α2_annulus α_in=0 reduction failed: "
        f"T = {result['T_inner_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha2_annulus_alpha_out_zero_reduction():
    r"""V_α2_annulus.e — at :math:`\alpha_{\rm out} = 0` (vacuum
    outer), the determinant is 1 and only :math:`T_{12} = \alpha_{
    \rm in}\,e^{-\tau_{\rm step}}` is nonzero off-diagonal.

    Physical meaning: with the outer surface in vacuum, through-rays
    that reach the outer surface escape. The reflective inner can
    re-launch them once but they leak at the next outer arrival.
    Outer-only rays escape at outer too.
    """
    result = derive_rank2_resolvent_annulus()
    assert result["pass_alpha_out_zero"], (
        f"V_α2_annulus α_out=0 reduction failed: "
        f"T = {result['T_outer_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha2_annulus_vacuum_vacuum_identity_reduction():
    r"""V_α2_annulus.f — at :math:`\alpha_{\rm in} = \alpha_{\rm out}
    = 0`, :math:`T = I`.

    Vacuum-vacuum BC removes the bouncing geometric series entirely;
    the closure reduces to :math:`\psi_{\rm surf} = T \cdot \alpha\,B
    = I \cdot 0 = 0`. The interior reconstruction is the bare
    first-leg integral.
    """
    result = derive_rank2_resolvent_annulus()
    assert result["pass_vacuum_identity"], (
        f"V_α2_annulus α_in=α_out=0 → T=I reduction failed: "
        f"T = {result['T_vac_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha2_annulus_overall_pass():
    """V_α2_annulus — composite gate."""
    result = derive_rank2_resolvent_annulus()
    assert result["pass"], f"V_α2_annulus composite failed: {result}"


# ═══════════════════════════════════════════════════════════════════════
# V_α2_annulus.aux — 3D chord scaling identity
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-annulus-3d-chord-scaling")
def test_v_alpha2_annulus_aux_3d_chord_scaling_step():
    r"""V_α2_annulus.aux.a — :math:`\tau_{\rm step}^{\rm annulus} =
    \tau_{\rm step}^{\rm hollow\,sph} / \sqrt{1 - \mu_{\rm axial}^2}`.

    The load-bearing chord-algebra identity that takes the hollow-
    sphere shell chord meaning to the cylinder 3D angular phase-space.
    Issue #129 angle-resolved discipline encoded as a SymPy identity:
    the cylinder Bickley-Naylor pre-integrated form would NOT factor
    cleanly like this; only the explicit 3D-lifted form does.
    """
    result = derive_3d_chord_scaling_annulus()
    assert result["pass_step_ratio"], (
        f"V_α2_annulus.aux 3D chord scaling ratio failed: "
        f"τ_annulus / τ_hollow_sph = {result['ratio']}"
    )


@pytest.mark.foundation
def test_v_alpha2_annulus_aux_3d_chord_scaling_period():
    r"""V_α2_annulus.aux.b — outer-only bounce period :math:`L_{\rm
    period}^{\rm 3D, annulus} = L_{\rm period}^{\rm 2D, cyl} /
    \sqrt{1 - \mu_{\rm axial}^2}`.

    Same axial-correction lift applies to the outer-only branch's
    bounce-period chord. This is the cylinder Phase-1
    :func:`derive_bounce_period_chord_cylinder` identity at the
    annulus outer surface.
    """
    result = derive_3d_chord_scaling_annulus()
    assert result["pass_period_ratio"], (
        f"V_α2_annulus.aux 3D bounce-period chord scaling failed."
    )


@pytest.mark.foundation
def test_v_alpha2_annulus_aux_in_plane_match():
    r"""V_α2_annulus.aux.c — at :math:`\mu_{\rm axial} = 0` (in-plane
    motion), the annulus 3D shell chord equals the hollow-sphere
    shell chord.

    Sanity check: when the trajectory is purely in-plane (no axial
    component), :math:`\sqrt{1 - \mu_{\rm axial}^2} = 1` and the
    cylinder 3D chord coincides with the in-plane 2D chord, which
    matches the hollow-sphere shell chord at the same impact
    parameter.
    """
    result = derive_3d_chord_scaling_annulus()
    assert result["pass_in_plane_match"], (
        f"V_α2_annulus.aux in-plane match failed."
    )


@pytest.mark.foundation
def test_v_alpha2_annulus_aux_overall_pass():
    """V_α2_annulus.aux — composite gate."""
    result = derive_3d_chord_scaling_annulus()
    assert result["pass"], (
        f"V_α2_annulus.aux composite failed: {result}"
    )


# ═══════════════════════════════════════════════════════════════════════
# V_α3_annulus — vacuum reduction at α_in = α_out = 0
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-annulus-architecture")
def test_v_alpha3_annulus_psi_surf_vanishes_at_vacuum_vacuum():
    r"""V_α3_annulus — surface fluxes
    :math:`\psi_{\rm in}^{\rm out} = \psi_{\rm out}^{\rm in} = 0`
    at :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0`.

    Both leading factors zero the entire closure components. The
    annulus prototype handles vacuum-vacuum BC with no special-case
    branch.
    """
    result = derive_alpha_zero_kernel_reduction_annulus()
    assert result["pass_substitution"] and result["pass_limit"], (
        f"V_α3_annulus vacuum-vacuum reduction failed: "
        f"ψ_in^out at zero = {result['psi_in_out_at_zero']}, "
        f"ψ_out^in at zero = {result['psi_out_in_at_zero']}"
    )


@pytest.mark.foundation
def test_v_alpha3_annulus_cavity_absorber_reduction():
    r"""V_α3_annulus — at :math:`\alpha_{\rm in} = 0` with
    :math:`\alpha_{\rm out} \in (0, 1]` (cavity absorber),
    :math:`\psi_{\rm in}^{\rm out} = 0` and
    :math:`\psi_{\rm out}^{\rm in} = \alpha_{\rm out}\,B_{\rm out}`.

    Through-rays that reach the inner cavity are lost; the bouncing
    chain breaks. The outer-surface outgoing flux receives only the
    direct one-transit contribution, without geometric-series
    amplification. det → 1.
    """
    result = derive_alpha_zero_kernel_reduction_annulus()
    assert result["pass_cavity_absorber"], (
        f"V_α3_annulus cavity-absorber reduction failed: "
        f"ψ_in^out = {result['psi_in_at_inner_vac']}, "
        f"ψ_out^in = {result['psi_out_at_inner_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha3_annulus_reflective_cavity_reduction():
    r"""V_α3_annulus — at :math:`\alpha_{\rm in} = 1, \alpha_{\rm out}
    = 0` (reflective cavity / vacuum outer),
    :math:`\psi_{\rm in}^{\rm out} = \alpha_{\rm in}\,B_{\rm in}`
    and :math:`\psi_{\rm out}^{\rm in} = 0`.

    Through-rays escape at the outer surface; the inner reflective
    surface provides only one launch per particle. det → 1.
    """
    result = derive_alpha_zero_kernel_reduction_annulus()
    assert result["pass_reflective_cavity"], (
        f"V_α3_annulus reflective-cavity reduction failed: "
        f"ψ_in^out = {result['psi_in_at_outer_vac']}, "
        f"ψ_out^in = {result['psi_out_at_outer_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha3_annulus_overall_pass():
    """V_α3_annulus — composite gate."""
    result = derive_alpha_zero_kernel_reduction_annulus()
    assert result["pass"], f"V_α3_annulus composite failed: {result}"
