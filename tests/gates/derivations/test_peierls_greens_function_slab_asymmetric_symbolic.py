r"""Phase-3B symbolic-vs-textbook contract tests for the **asymmetric**
slab Variant α Green's function reference (rank-2 boundary-to-boundary
scattering resolvent, independent :math:`\alpha_L, \alpha_R \in [0, 1]`).

Math-origin pattern: the SymPy derivation in
:mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_slab_asymmetric`
is the source of truth for the operator-level identities of asymmetric
slab Variant α. Mirrors the Phase-3A symmetric-slab symbolic test gates
with the rank-2 generalisation.

Three operator-level verifications
-----------------------------------

V_α1_slab_asym. **Closed-slab bounce-sum self-consistency at corner
   α_L=α_R=1** (rank-2 form). Constant source :math:`q` produces
   :math:`\psi(x, \mu) = q/\Sigma_t` everywhere when both walls are
   fully reflective. At :math:`(\alpha_L, \alpha_R) \in (0, 1)^2` with
   either factor :math:`< 1`, the surface flux is **strictly less**
   than :math:`q/\Sigma_t` (leakage) — a leaky-mode sample point is
   verified.

V_α2_slab_asym. **Rank-2 resolvent algebra**. Direct SymPy matrix
   inversion of :math:`(I - S)` with :math:`S = \mathrm{antidiag}(
   \alpha_L\,e^{-\tau}, \alpha_R\,e^{-\tau})`. Verifies the canonical
   2x2 form, the determinant :math:`1 - \alpha_L\,\alpha_R\,e^{-2\tau}`,
   and three special-case reductions: symmetric :math:`\alpha_L =
   \alpha_R`, one-vacuum-wall :math:`\alpha_L = 0` or :math:`\alpha_R
   = 0`, and vacuum-vacuum :math:`T = I`.

V_α3_slab_asym. **Vacuum reduction at α_L=α_R=0**. The rank-2 closure
   collapses to the bare first-leg integral with no special-case
   branch; one-vacuum-wall reductions also verified.

Predecessor / sibling tests
---------------------------

- :mod:`.test_trajectory_resolvent_slab_symbolic` — Phase-3A
  symmetric slab V_α1_slab/V_α2_slab/V_α3_slab.
- :mod:`.test_trajectory_resolvent_symbolic` — sphere V_α1/V_α2/V_α3.

References
----------

- :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 3B
  asymmetric slab plan.
- :file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
  — frame match: the rank-2 resolvent is the predicted generalisation
  of Variant α to two-surface topologies with independent BCs.
"""
from __future__ import annotations

import pytest

from orpheus.derivations.continuous.trajectory_resolvent.origins.specular import (
    derive_alpha_zero_kernel_reduction_slab_asymmetric,
    derive_operator_constant_trial_closed_slab_asymmetric,
    derive_rank2_resolvent_slab_asymmetric,
)


# ═══════════════════════════════════════════════════════════════════════
# V_α1_slab_asym — closed-slab bounce-sum self-consistency (α_L=α_R=1)
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-slab-asym-closure")
def test_v_alpha1_slab_asym_psi_L_plus_at_closed_BC_equals_q_over_sigma_t():
    r"""V_α1_slab_asym.a — :math:`\psi_L^+(\mu) = q/\Sigma_t` at closed
    BC :math:`\alpha_L = \alpha_R = 1`.

    The rank-2 closure on a constant source :math:`q` with both walls
    fully reflective (:math:`\alpha_L = \alpha_R = 1`) reduces (via
    the algebraic identity :math:`(1-e^{-\tau})(1+e^{-\tau}) =
    1 - e^{-2\tau}`) to the chord-length-independent
    :math:`\psi_L^+ = q/\Sigma_t`. This is the closed-asymmetric-slab
    analog of V_α1_slab's surface fixed-point and the load-bearing
    rank-2-form algebraic identity.
    """
    result = derive_operator_constant_trial_closed_slab_asymmetric()
    assert result["pass_psi_L_closed_equals_q_over_sigma_t"], (
        f"V_α1_slab_asym ψ_L^+ at closed BC failed: "
        f"got {result['psi_L_plus_closed']}"
    )


@pytest.mark.foundation
def test_v_alpha1_slab_asym_psi_R_minus_at_closed_BC_equals_q_over_sigma_t():
    r"""V_α1_slab_asym.b — :math:`\psi_R^-(\mu) = q/\Sigma_t` at closed
    BC :math:`\alpha_L = \alpha_R = 1`.

    Symmetric counterpart of the :math:`\psi_L^+` identity — both
    surface flux components must independently equal :math:`q/\Sigma_t`
    at closed BC.
    """
    result = derive_operator_constant_trial_closed_slab_asymmetric()
    assert result["pass_psi_R_closed_equals_q_over_sigma_t"], (
        f"V_α1_slab_asym ψ_R^- at closed BC failed: "
        f"got {result['psi_R_minus_closed']}"
    )


@pytest.mark.foundation
def test_v_alpha1_slab_asym_total_psi_constant_at_closed_BC():
    r"""V_α1_slab_asym.c — :math:`\psi(x, \mu) = q/\Sigma_t` everywhere
    at closed BC.

    First-leg attenuation cancels against attenuated surface flux,
    making the interior reconstruction independent of position
    :math:`x` and direction :math:`\mu`. Algebraically identical to
    V_α1_slab.
    """
    result = derive_operator_constant_trial_closed_slab_asymmetric()
    assert result["pass_psi_total_closed_constant"], (
        f"V_α1_slab_asym total ψ constant failed: "
        f"got {result['psi_total_closed']}"
    )


@pytest.mark.foundation
def test_v_alpha1_slab_asym_operator_on_constant_gives_omega_0():
    r"""V_α1_slab_asym.d — operator action on isotropic trial gives
    :math:`(K \cdot 1) = \omega_0 = \Sigma_s/\Sigma_t` at closed BC.

    Yields :math:`k_{\rm eff} = k_\infty = \nu\Sigma_f/\Sigma_a` for
    the closed asymmetric slab — the rank-2 framework reproduces the
    closed-slab eigenvalue invariant when both reflectivities are 1.
    """
    result = derive_operator_constant_trial_closed_slab_asymmetric()
    assert result["pass_eigenvalue"], (
        f"V_α1_slab_asym K·1 failed: got "
        f"{result['K_on_constant_trial_closed']}, expected "
        f"{result['omega_0']}"
    )


@pytest.mark.foundation
def test_v_alpha1_slab_asym_leaky_corner_strictly_below_q_over_sigma_t():
    r"""V_α1_slab_asym.e — at :math:`(\alpha_L, \alpha_R) = (1/2, 1)`
    the surface flux is strictly less than :math:`q/\Sigma_t`.

    Leakage out the partial-reflection wall keeps the surface flux
    below :math:`q/\Sigma_t` whenever **either** wall is sub-unity.
    This is the rank-2 form's signature distinction from the closed-
    slab special case — leaky modes ARE NOT eigenmodes of the
    closed-slab operator.
    """
    result = derive_operator_constant_trial_closed_slab_asymmetric()
    assert result["pass_leaky_below"], (
        f"V_α1_slab_asym leaky-mode below q/Σ_t failed: "
        f"sample at (α_L=1/2, α_R=1, τ=1, q/Σ_t=1) = "
        f"{result['psi_sample_leaky']}"
    )


@pytest.mark.foundation
def test_v_alpha1_slab_asym_overall_pass():
    """V_α1_slab_asym — composite gate."""
    result = derive_operator_constant_trial_closed_slab_asymmetric()
    assert result["pass"], f"V_α1_slab_asym composite failed: {result}"


# ═══════════════════════════════════════════════════════════════════════
# V_α2_slab_asym — rank-2 resolvent T = (I - S)^{-1}
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-slab-asym-resolvent")
def test_v_alpha2_slab_asym_determinant_canonical_form():
    r"""V_α2_slab_asym.a — :math:`\det(I - S) = 1 - \alpha_L\,
    \alpha_R\,e^{-2\tau}` via SymPy direct matrix determinant.

    The rank-2 monodromy :math:`S = \mathrm{antidiag}(\alpha_L\,
    e^{-\tau}, \alpha_R\,e^{-\tau})` yields :math:`\det(I - S) =
    1 - (\alpha_L\,e^{-\tau})(\alpha_R\,e^{-\tau}) = 1 - \alpha_L\,
    \alpha_R\,e^{-2\tau}`. This determinant carries the entire
    closed-form algebra of the rank-2 resolvent — the singular locus
    :math:`\alpha_L\alpha_R\,e^{-2\tau} = 1` is reachable only at
    :math:`\alpha_L = \alpha_R = 1, \tau = 0` (grazing rays in a
    closed slab), which is the same divergent locus as the rank-1
    symmetric form.
    """
    result = derive_rank2_resolvent_slab_asymmetric()
    assert result["pass_det"], (
        f"V_α2_slab_asym determinant failed: got {result['det_M']}, "
        f"expected {result['det_canonical']}"
    )


@pytest.mark.foundation
def test_v_alpha2_slab_asym_canonical_T_form():
    r"""V_α2_slab_asym.b — the rank-2 resolvent matches the canonical
    closed form

    .. math::

       T = \frac{1}{\det}\,\begin{pmatrix} 1 & \alpha_L\,e^{-\tau}
       \\ \alpha_R\,e^{-\tau} & 1 \end{pmatrix}.

    SymPy's ``Matrix.inv()`` produces the inverse, and each entry is
    matched against the canonical form via ``sp.simplify``.
    """
    result = derive_rank2_resolvent_slab_asymmetric()
    assert result["pass_T_form"], (
        f"V_α2_slab_asym canonical T form failed: T = {result['T_resolvent']}"
    )


@pytest.mark.foundation
def test_v_alpha2_slab_asym_symmetric_reduction():
    r"""V_α2_slab_asym.c — at :math:`\alpha_L = \alpha_R = \alpha`,
    :math:`T_{11} + T_{12} = 1/(1 - \alpha\,e^{-\tau})`.

    The rank-2 closure on a constant source — at symmetric BC and
    factoring out the leading :math:`\alpha\,B` — produces the scalar
    factor :math:`(T_{11} + T_{12})`, which simplifies to
    :math:`1/(1 - \alpha\,e^{-\tau})`. This is **NOT** the same
    arithmetic form as the rank-1 symmetric :math:`1/(1 - \alpha^2\,
    e^{-2\tau})`; the bridge is the structural identity
    :math:`B_{\rm period} = (1 + e^{-\tau})\,B_{\rm single\,transit}`
    on a constant source. SymPy verifies the rank-2 simplification
    here; the numerical reduce-to-symmetric consistency check is
    deferred to the solver-level test in
    ``test_trajectory_resolvent_slab_asymmetric_solver.py``.
    """
    result = derive_rank2_resolvent_slab_asymmetric()
    assert result["pass_symmetric_simplification"], (
        f"V_α2_slab_asym symmetric reduction failed."
    )


@pytest.mark.foundation
def test_v_alpha2_slab_asym_alpha_L_zero_reduction():
    r"""V_α2_slab_asym.d — at :math:`\alpha_L = 0`, the determinant is
    1 and only :math:`T_{21} = \alpha_R\,e^{-\tau}` is nonzero
    off-diagonal.

    With no left-wall reflection, the right wall outgoing flux receives
    NO contribution from the right-wall-injected B integral via the
    bouncing geometric series — only the direct :math:`B_{LR}` term.
    """
    result = derive_rank2_resolvent_slab_asymmetric()
    assert result["pass_alpha_L_zero"], (
        f"V_α2_slab_asym α_L=0 reduction failed: "
        f"T = {result['T_left_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha2_slab_asym_alpha_R_zero_reduction():
    r"""V_α2_slab_asym.e — at :math:`\alpha_R = 0`, the determinant is
    1 and only :math:`T_{12} = \alpha_L\,e^{-\tau}` is nonzero
    off-diagonal.

    Mirror of the :math:`\alpha_L = 0` case. Combined with V_α2_slab_asym.d,
    these establish the rank-2 form's correct reduction on the
    one-vacuum-wall edges of the BC parameter square.
    """
    result = derive_rank2_resolvent_slab_asymmetric()
    assert result["pass_alpha_R_zero"], (
        f"V_α2_slab_asym α_R=0 reduction failed: "
        f"T = {result['T_right_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha2_slab_asym_vacuum_vacuum_identity_reduction():
    r"""V_α2_slab_asym.f — at :math:`\alpha_L = \alpha_R = 0`,
    :math:`T = I`.

    Vacuum-vacuum BC removes the bouncing geometric series entirely:
    the rank-2 resolvent collapses to the identity, and the closure
    reduces to :math:`\psi_{\rm surf} = T \cdot \alpha\,B = I \cdot 0
    = 0`. The interior reconstruction is the bare first-leg integral.
    """
    result = derive_rank2_resolvent_slab_asymmetric()
    assert result["pass_vacuum_identity"], (
        f"V_α2_slab_asym α_L=α_R=0 → T=I reduction failed: "
        f"T = {result['T_vac_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha2_slab_asym_overall_pass():
    """V_α2_slab_asym — composite gate."""
    result = derive_rank2_resolvent_slab_asymmetric()
    assert result["pass"], f"V_α2_slab_asym composite failed: {result}"


# ═══════════════════════════════════════════════════════════════════════
# V_α3_slab_asym — vacuum reduction at α_L = α_R = 0
# ═══════════════════════════════════════════════════════════════════════


@pytest.mark.foundation
@pytest.mark.verifies("peierls-greens-slab-asym-architecture")
def test_v_alpha3_slab_asym_psi_surf_vanishes_at_vacuum_vacuum():
    r"""V_α3_slab_asym — surface fluxes :math:`\psi_L^+ = \psi_R^- = 0`
    at :math:`\alpha_L = \alpha_R = 0`.

    Both leading factors :math:`\alpha_L` (in :math:`\psi_L^+`) and
    :math:`\alpha_R` (in :math:`\psi_R^-`) zero the entire closure
    components. The rank-2 prototype handles vacuum-vacuum BC with no
    special-case branch. The interior reconstruction reduces to the
    bare first-leg integral, equivalent to the Phase-3A symmetric
    α=0 case.
    """
    result = derive_alpha_zero_kernel_reduction_slab_asymmetric()
    assert result["pass_substitution"] and result["pass_limit"], (
        f"V_α3_slab_asym vacuum-vacuum reduction failed: "
        f"ψ_L^+ at zero = {result['psi_L_plus_at_zero']}, "
        f"ψ_R^- at zero = {result['psi_R_minus_at_zero']}"
    )


@pytest.mark.foundation
def test_v_alpha3_slab_asym_one_vacuum_wall_left_reduction():
    r"""V_α3_slab_asym — at :math:`\alpha_L = 0` with :math:`\alpha_R
    \in (0, 1]`, :math:`\psi_L^+ = 0` and :math:`\psi_R^- = \alpha_R\,
    B_{RL}`.

    With no left-wall reflection, no flux ever returns from the left
    wall, so :math:`\psi_L^+ = 0` identically. The right-wall outgoing
    flux receives only the direct :math:`B_{RL}` contribution
    (:math:`\alpha_R\,B_{RL}` after the right-wall reflection),
    without the bouncing geometric series amplification (det → 1 when
    :math:`\alpha_L = 0`). This is the canonical method-of-images
    setup before the symmetry-extension test.
    """
    result = derive_alpha_zero_kernel_reduction_slab_asymmetric()
    assert result["pass_one_vacuum_wall_left"], (
        f"V_α3_slab_asym α_L=0 reduction failed: "
        f"ψ_L^+ = {result['psi_L_at_left_vac']}, "
        f"ψ_R^- = {result['psi_R_at_left_vac']}"
    )


@pytest.mark.foundation
def test_v_alpha3_slab_asym_overall_pass():
    """V_α3_slab_asym — composite gate."""
    result = derive_alpha_zero_kernel_reduction_slab_asymmetric()
    assert result["pass"], f"V_α3_slab_asym composite failed: {result}"
