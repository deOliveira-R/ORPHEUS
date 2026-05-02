r"""Phase-3C-1 symbolic-vs-textbook contract tests for the **hollow
sphere** Variant α Green's function reference (rank-2 boundary-to-
boundary scattering resolvent on through-rays + rank-1 outer-only
closure, independent :math:`\alpha_{\rm in}, \alpha_{\rm out} \in
[0, 1]`).

Math-origin pattern: the SymPy derivation in
:mod:`orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function_hollow_sphere`
is the source of truth for the operator-level identities of hollow-
sphere Variant α. Mirrors the Phase-3B asymmetric-slab symbolic test
gates with the curvilinear 2-surface generalisation (impact-parameter
phase-space partition).

Three operator-level verifications
-----------------------------------

V_α1_hollow_sph. **Closed-shell bounce-sum self-consistency at
   corner** :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1`. With both
   surfaces fully reflective and constant volumetric source :math:`q`,
   BOTH the outer-only branch (rank-1, :math:`b > R_{\rm in}`) AND
   the through-ray branch (rank-2, :math:`b \le R_{\rm in}`)
   independently produce :math:`\psi(r, \mu) = q/\Sigma_t`. The
   composability of the impact-parameter partition with the rank-2
   closure is the load-bearing structural claim of Phase 3C-1.

V_α2_hollow_sph. **Rank-2 resolvent algebra** specialised to the
   shell. Direct SymPy matrix inversion of :math:`(I - S)` with
   :math:`S = \mathrm{antidiag}(\alpha_{\rm in}\,e^{-\tau_{\rm step}},
   \alpha_{\rm out}\,e^{-\tau_{\rm step}})`. Verifies the canonical
   :math:`2 \times 2` form, the determinant
   :math:`1 - \alpha_{\rm in}\,\alpha_{\rm out}\,e^{-2\tau_{\rm step}}`,
   and four special-case reductions: symmetric :math:`\alpha_{\rm in}
   = \alpha_{\rm out} = \alpha`, vacuum inner (cavity absorber),
   vacuum outer (reflective cavity), vacuum-vacuum :math:`T = I`.

V_α3_hollow_sph. **Vacuum reduction at**
   :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0`. Both leading-α
   factors zero the closure on each surface component; the kernel
   reduces to bare first-leg integrals on both phase-space subsets
   (rank-1 outer-only AND rank-2 through-ray) with no special-case
   branch needed.

Predecessor / sibling tests
---------------------------

- :mod:`.test_peierls_greens_function_slab_asymmetric_symbolic` —
  Phase-3B asymmetric slab V_α (the rank-2 template the hollow-sphere
  rank-2 closure is lifted from).
- :mod:`.test_peierls_greens_function_symbolic` — sphere V_α1/V_α2/V_α3
  (the rank-1 template the outer-only branch reuses).
- :mod:`.test_peierls_greens_function_cylinder_symbolic` — cylinder
  V_α1/V_α2/V_α3.

References
----------

- :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 3C-1
  hollow sphere plan.
- :file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
  — frame match: rank-2 BIE block resolvent at 2-surface topologies.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.peierls_greens_function.origins.specular import (
    derive_alpha_zero_kernel_reduction_hollow_sphere,
    derive_operator_constant_trial_closed_hollow_sphere,
    derive_rank2_resolvent_hollow_sphere,
)


# ═══════════════════════════════════════════════════════════════════════
# V_α1_hollow_sph — closed-shell bounce-sum self-consistency
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-hollow-sph-architecture")
def test_v_alpha1_hollow_sph_outer_only_closed_BC_equals_q_over_sigma_t():
    r"""V_α1_hollow_sph.a — outer-only branch (:math:`b > R_{\rm in}`)
    at :math:`\alpha_{\rm out} = 1` gives :math:`\psi_{\rm surf} =
    q/\Sigma_t` on a constant source.

    Identical to V_α1 sphere — outer-only rays in the hollow sphere
    are structurally identical to solid-sphere rays at the same outer
    radius and impact parameter. The inner cavity is "invisible" to
    these trajectories. The rank-1 closure
    :math:`\psi_{\rm surf} = B/(1 - e^{-\tau_{\rm period}})` collapses
    to :math:`q/\Sigma_t` by the same identity that drives V_α1 sphere.
    """
    result = derive_operator_constant_trial_closed_hollow_sphere()
    assert result["pass_outer_closed_equals_q_over_sigma_t"], (
        f"V_α1_hollow_sph outer-only at α_out=1 failed: "
        f"got {result['psi_surf_outer_closed']}"
    )


@pytest.mark.foundation
def test_v_alpha1_hollow_sph_psi_in_out_at_closed_BC_equals_q_over_sigma_t():
    r"""V_α1_hollow_sph.b — through-ray branch (:math:`b \le R_{\rm
    in}`) at :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1` gives
    :math:`\psi_{\rm in}^{\rm out} = q/\Sigma_t`.

    The rank-2 closure on a constant source with both surfaces fully
    reflective collapses (via :math:`(1 - e^{-\tau})(1 + e^{-\tau}) =
    1 - e^{-2\tau}`) to the chord-length-independent
    :math:`q/\Sigma_t`. Algebraically identical to V_α1_slab_asym
    closed-BC corner — the curvilinear chord :math:`\tau_{\rm step}`
    plays the role of the slab :math:`\tau`.
    """
    result = derive_operator_constant_trial_closed_hollow_sphere()
    assert result["pass_psi_in_out_closed_equals_q_over_sigma_t"], (
        f"V_α1_hollow_sph ψ_in^out at closed BC failed: "
        f"got {result['psi_in_out_closed']}"
    )


@pytest.mark.foundation
def test_v_alpha1_hollow_sph_psi_out_in_at_closed_BC_equals_q_over_sigma_t():
    r"""V_α1_hollow_sph.c — through-ray :math:`\psi_{\rm out}^{\rm in}
    = q/\Sigma_t` at closed BC.

    Symmetric counterpart of the :math:`\psi_{\rm in}^{\rm out}`
    identity — both surface flux components must independently equal
    :math:`q/\Sigma_t` at closed BC. The structural composability with
    the outer-only branch is what makes the impact-parameter partition
    sound.
    """
    result = derive_operator_constant_trial_closed_hollow_sphere()
    assert result["pass_psi_out_in_closed_equals_q_over_sigma_t"], (
        f"V_α1_hollow_sph ψ_out^in at closed BC failed: "
        f"got {result['psi_out_in_closed']}"
    )


@pytest.mark.foundation
def test_v_alpha1_hollow_sph_operator_outer_branch_gives_omega_0():
    r"""V_α1_hollow_sph.d — operator action on isotropic trial gives
    :math:`(K \cdot 1)_{\rm outer} = \omega_0 = \Sigma_s/\Sigma_t` on
    the outer-only subset.
    """
    result = derive_operator_constant_trial_closed_hollow_sphere()
    assert result["pass_eigenvalue_outer"], (
        f"V_α1_hollow_sph K·1 (outer) failed: got "
        f"{result['K_outer_on_constant_trial']}, expected "
        f"{result['omega_0']}"
    )


@pytest.mark.foundation
def test_v_alpha1_hollow_sph_operator_through_branch_gives_omega_0():
    r"""V_α1_hollow_sph.e — operator action on isotropic trial gives
    :math:`(K \cdot 1)_{\rm through} = \omega_0 = \Sigma_s/\Sigma_t`
    on the through-ray subset.

    Combined with the outer-branch identity, this proves the closed-
    shell eigenvalue is :math:`k_{\rm eff} = k_\infty =
    \nu\Sigma_f/\Sigma_a` IDENTICALLY across BOTH phase-space
    subsets. Failure here means the impact-parameter partition does
    not compose cleanly with the rank-2 closure — a structural bug
    that would HALT Phase 3C-1.
    """
    result = derive_operator_constant_trial_closed_hollow_sphere()
    assert result["pass_eigenvalue_inner"], (
        f"V_α1_hollow_sph K·1 (through-ray) failed: got "
        f"{result['K_inner_on_constant_trial']}, expected "
        f"{result['omega_0']}"
    )


@pytest.mark.foundation
def test_v_alpha1_hollow_sph_leaky_corner_strictly_below_q_over_sigma_t():
    r"""V_α1_hollow_sph.f — at :math:`(\alpha_{\rm in}, \alpha_{\rm
    out}) = (1/2, 1)` the through-ray surface flux is strictly less
    than :math:`q/\Sigma_t`.

    Leakage at the partial-reflection inner surface keeps the through-
    ray surface flux below :math:`q/\Sigma_t`. The rank-2 form's
    distinction from the closed-shell eigenmode.
    """
    result = derive_operator_constant_trial_closed_hollow_sphere()
    assert result["pass_leaky_below"], (
        f"V_α1_hollow_sph leaky-mode below q/Σ_t failed: "
        f"sample = {result['psi_sample_leaky']}"
    )


@pytest.mark.foundation
def test_v_alpha1_hollow_sph_overall_pass():
    """V_α1_hollow_sph — composite gate."""
    result = derive_operator_constant_trial_closed_hollow_sphere()
    assert result["pass"], f"V_α1_hollow_sph composite failed: {result}"


# ═══════════════════════════════════════════════════════════════════════
# V_α2_hollow_sph — rank-2 resolvent T = (I - S)^{-1}
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-hollow-sph-through-rank2")
def test_v_alpha2_hollow_sph_determinant_canonical_form():
    r"""V_α2_hollow_sph.a — :math:`\det(I - S) = 1 - \alpha_{\rm in}\,
    \alpha_{\rm out}\,e^{-2\tau_{\rm step}}` via SymPy direct matrix
    determinant.

    Same algebraic structure as the slab-asymmetric determinant; only
    the meaning of :math:`\tau_{\rm step}` changes (curvilinear shell
    chord vs slab chord). The singular locus :math:`\alpha_{\rm in}\,
    \alpha_{\rm out}\,e^{-2\tau_{\rm step}} = 1` is reachable only at
    :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1` and
    :math:`\tau_{\rm step} = 0` (tangent rays at :math:`b = R_{\rm
    in}` in a closed shell).
    """
    result = derive_rank2_resolvent_hollow_sphere()
    assert result["pass_det"], (
        f"V_α2_hollow_sph determinant failed: got {result['det_M']}, "
        f"expected {result['det_canonical']}"
    )


@pytest.mark.foundation
def test_v_alpha2_hollow_sph_canonical_T_form():
    r"""V_α2_hollow_sph.b — the rank-2 resolvent matches the canonical
    closed form
    :math:`T = (1/\det)\,\bigl[[1, \alpha_{\rm in} e^{-\tau_{\rm step}}],
    [\alpha_{\rm out} e^{-\tau_{\rm step}}, 1]\bigr]`.

    SymPy's ``Matrix.inv()`` produces the inverse, and each entry is
    matched against the canonical form via ``sp.simplify``.
    """
    result = derive_rank2_resolvent_hollow_sphere()
    assert result["pass_T_form"], (
        f"V_α2_hollow_sph canonical T form failed: "
        f"T = {result['T_resolvent']}"
    )


@pytest.mark.foundation
def test_v_alpha2_hollow_sph_symmetric_reduction():
    r"""V_α2_hollow_sph.c — at :math:`\alpha_{\rm in} = \alpha_{\rm
    out} = \alpha`, :math:`T_{11} + T_{12} = 1/(1 - \alpha\,
    e^{-\tau_{\rm step}})`.

    The rank-2 → rank-1 collapse on a constant source with symmetric
    BC. Same algebraic identity as V_α2_slab_asym.c.
    """
    result = derive_rank2_resolvent_hollow_sphere()
    assert result["pass_symmetric_simplification"], (
        f"V_α2_hollow_sph symmetric reduction failed."
    )


@pytest.mark.foundation
def test_v_alpha2_hollow_sph_alpha_in_zero_reduction():
    r"""V_α2_hollow_sph.d — at :math:`\alpha_{\rm in} = 0` (cavity
    absorber), the determinant is 1 and only :math:`T_{21} = \alpha_{
    \rm out}\,e^{-\tau_{\rm step}}` is nonzero off-diagonal.

    Physical meaning: with the inner surface absorbing perfectly, the
    cavity acts as a "particle sink." Through-rays that reach the
    inner surface are lost; the bouncing geometric series is broken
    after the first traversal.
    """
    result = derive_rank2_resolvent_hollow_sphere()
    assert result["pass_alpha_in_zero"], (
        f"V_α2_hollow_sph α_in=0 reduction failed: "
        f"T = {result['T_inner_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha2_hollow_sph_alpha_out_zero_reduction():
    r"""V_α2_hollow_sph.e — at :math:`\alpha_{\rm out} = 0` (vacuum
    outer), the determinant is 1 and only :math:`T_{12} = \alpha_{
    \rm in}\,e^{-\tau_{\rm step}}` is nonzero off-diagonal.

    Physical meaning: with the outer surface in vacuum, through-rays
    that reach the outer surface escape. The reflective inner surface
    can re-launch them once but they leak out at the next outer
    arrival. (The outer-only subset has :math:`\alpha_{\rm out} = 0`
    rank-1 closure with :math:`T = 1` — outer-only rays escape on the
    first arrival.)
    """
    result = derive_rank2_resolvent_hollow_sphere()
    assert result["pass_alpha_out_zero"], (
        f"V_α2_hollow_sph α_out=0 reduction failed: "
        f"T = {result['T_outer_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha2_hollow_sph_vacuum_vacuum_identity_reduction():
    r"""V_α2_hollow_sph.f — at :math:`\alpha_{\rm in} = \alpha_{\rm
    out} = 0`, :math:`T = I`.

    Vacuum-vacuum BC removes the bouncing geometric series entirely;
    the closure reduces to :math:`\psi_{\rm surf} = T \cdot \alpha\,
    B = I \cdot 0 = 0`. The interior reconstruction is the bare
    first-leg integral.
    """
    result = derive_rank2_resolvent_hollow_sphere()
    assert result["pass_vacuum_identity"], (
        f"V_α2_hollow_sph α_in=α_out=0 → T=I reduction failed: "
        f"T = {result['T_vac_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha2_hollow_sph_overall_pass():
    """V_α2_hollow_sph — composite gate."""
    result = derive_rank2_resolvent_hollow_sphere()
    assert result["pass"], f"V_α2_hollow_sph composite failed: {result}"


# ═══════════════════════════════════════════════════════════════════════
# V_α3_hollow_sph — vacuum reduction at α_in = α_out = 0
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-hollow-sph-architecture")
def test_v_alpha3_hollow_sph_psi_surf_vanishes_at_vacuum_vacuum():
    r"""V_α3_hollow_sph — surface fluxes
    :math:`\psi_{\rm in}^{\rm out} = \psi_{\rm out}^{\rm in} = 0`
    at :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0`.

    Both leading factors :math:`\alpha_{\rm in}` (in
    :math:`\psi_{\rm in}^{\rm out}`) and :math:`\alpha_{\rm out}` (in
    :math:`\psi_{\rm out}^{\rm in}`) zero the entire closure
    components. The hollow-sphere prototype handles vacuum-vacuum BC
    with no special-case branch on the through-ray subset.
    """
    result = derive_alpha_zero_kernel_reduction_hollow_sphere()
    assert result["pass_substitution"] and result["pass_limit"], (
        f"V_α3_hollow_sph vacuum-vacuum reduction failed: "
        f"ψ_in^out at zero = {result['psi_in_out_at_zero']}, "
        f"ψ_out^in at zero = {result['psi_out_in_at_zero']}"
    )


@pytest.mark.foundation
def test_v_alpha3_hollow_sph_cavity_absorber_reduction():
    r"""V_α3_hollow_sph — at :math:`\alpha_{\rm in} = 0` with
    :math:`\alpha_{\rm out} \in (0, 1]` (cavity absorber),
    :math:`\psi_{\rm in}^{\rm out} = 0` and
    :math:`\psi_{\rm out}^{\rm in} = \alpha_{\rm out}\,B_{\rm out}`.

    With the inner surface absorbing perfectly, no through-ray ever
    returns to the inner surface (the bouncing chain breaks after
    one absorbed inner arrival). The outer-surface outgoing flux
    receives only the direct contribution :math:`\alpha_{\rm out}\,
    B_{\rm out}` (one transit + one reflection), without the
    geometric-series amplification. det → 1 since
    :math:`\alpha_{\rm in}\,\alpha_{\rm out}\,e^{-2\tau_{\rm step}} =
    0`.
    """
    result = derive_alpha_zero_kernel_reduction_hollow_sphere()
    assert result["pass_cavity_absorber"], (
        f"V_α3_hollow_sph cavity-absorber reduction failed: "
        f"ψ_in^out = {result['psi_in_at_inner_vac']}, "
        f"ψ_out^in = {result['psi_out_at_inner_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha3_hollow_sph_reflective_cavity_reduction():
    r"""V_α3_hollow_sph — at :math:`\alpha_{\rm in} = 1, \alpha_{\rm
    out} = 0` (reflective cavity / vacuum outer),
    :math:`\psi_{\rm in}^{\rm out} = \alpha_{\rm in}\,B_{\rm in}`
    and :math:`\psi_{\rm out}^{\rm in} = 0`.

    Through-rays escape at the outer surface (no return); the inner
    reflective surface provides only one launch per particle. det → 1.
    """
    result = derive_alpha_zero_kernel_reduction_hollow_sphere()
    assert result["pass_reflective_cavity"], (
        f"V_α3_hollow_sph reflective-cavity reduction failed: "
        f"ψ_in^out = {result['psi_in_at_outer_vac']}, "
        f"ψ_out^in = {result['psi_out_at_outer_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha3_hollow_sph_overall_pass():
    """V_α3_hollow_sph — composite gate."""
    result = derive_alpha_zero_kernel_reduction_hollow_sphere()
    assert result["pass"], f"V_α3_hollow_sph composite failed: {result}"
