r"""SymPy derivation — operator-level identities for the **hollow sphere**
Variant α Green's function reference (3D radial, independent specular
reflectivities :math:`\alpha_{\rm in} \in [0, 1]` on the inner cavity
surface :math:`r = R_{\rm in}` and :math:`\alpha_{\rm out} \in [0, 1]`
on the outer surface :math:`r = R_{\rm out}`).

Phase-3C-1 extension of the slab-asymmetric Variant α (per
:file:`/.claude/plans/peierls-greens-cylinder-and-2bc.md`). Mirrors the
Phase-3B asymmetric-slab proofs in
:mod:`.greens_function_slab_asymmetric` with the **rank-2 boundary-to-
boundary scattering resolvent** specialised to a CURVILINEAR 2-surface
geometry.

Phase-space partition by impact parameter
------------------------------------------

The hollow sphere has two surfaces (:math:`r = R_{\rm in}` and
:math:`r = R_{\rm out}`) but a particle's trajectory line interacts
with at most ONE OR TWO of them depending on its **impact parameter**
:math:`b = r\sqrt{1 - \mu^2}`. Two cases:

- :math:`b > R_{\rm in}` — outer-only rays. The straight-line
  trajectory does not intersect the inner sphere; the particle
  bounces between two points on the OUTER surface only. Closure
  rank-1 (orbit-space class: one-surface-compact), structurally
  identical to a solid sphere ray at the same outer radius and
  impact parameter.
- :math:`b \le R_{\rm in}` — through-rays. The trajectory line
  crosses the inner cavity. Under interpretation (A) — see plan §
  "Hollow sphere geometry" — the inner surface is a specular reflector
  with reflectivity :math:`\alpha_{\rm in}`, so the particle bounces
  alternately between the inner and outer surfaces. Closure rank-2
  (orbit-space class: two-surface), structurally identical to the
  asymmetric slab with curved chord algebra: both have an M/G
  interval with two distinct BC endpoints, only the chord algebra
  and the G-equivariant lift back to the higher-dim M differ.

Interpretation (A) treats :math:`\alpha_{\rm in} = 0` as
"particle is lost to the cavity" (perfect inner absorber) and
:math:`\alpha_{\rm in} = 1` as "perfect inner reflector." This makes
the hollow sphere a clean **orbit-space M/G analog** of the
asymmetric slab. See Sphinx §\ ``orbit-space-m-g-classification``
for the structural signature that drives this analog.

Through-ray rank-2 monodromy
----------------------------

For :math:`b \le R_{\rm in}`, the surface state is

.. math::

   \psi_{\rm surf}(b, \mu) = \begin{pmatrix}
       \psi_{\rm in}^{\rm out}(b, \mu) \\
       \psi_{\rm out}^{\rm in}(b, \mu)
   \end{pmatrix}

(outgoing flux at the inner surface in the outward direction, and
outgoing flux at the outer surface in the inward direction). One step
of the bouncing trajectory transports :math:`\psi_{\rm in}^{\rm out}`
to :math:`\psi_{\rm out}^{\rm in}` (transit + reflection at the outer
wall) and :math:`\psi_{\rm out}^{\rm in}` to
:math:`\psi_{\rm in}^{\rm out}` (transit + reflection at the inner
wall) with **single-transit** optical depth

.. math::

   \tau_{\rm step}(b) = \Sigma_t \cdot \bigl(\sqrt{R_{\rm out}^2 - b^2}
                                            - \sqrt{R_{\rm in}^2 - b^2}\bigr).

The monodromy operator is

.. math::

   S(\alpha_{\rm in}, \alpha_{\rm out}, \tau_{\rm step}) =
       \begin{pmatrix}
           0                                             & \alpha_{\rm in}\,e^{-\tau_{\rm step}} \\
           \alpha_{\rm out}\,e^{-\tau_{\rm step}}        & 0
       \end{pmatrix},

and the rank-2 resolvent is :math:`T = (I - S)^{-1}`.

The diagonal of :math:`S` is identically zero — one step never returns
to the same surface. The off-diagonal coupling is encoded by the
SHELL-MATERIAL chord (cavity transit is in vacuum, optical-depth zero,
and is "behind the mirror" under interpretation (A)). The key
structural identity carrying this geometry into the slab framework:

.. math::

   \boxed{\tau_{\rm step}^{\rm hollow\,sph}(b)
       = \Sigma_t \bigl(\sqrt{R_{\rm out}^2 - b^2}
                       - \sqrt{R_{\rm in}^2 - b^2}\bigr).}

This is the slab :math:`\tau = \Sigma_t L /|\mu|` with the curvilinear
shell-traversal length playing the role of :math:`L /|\mu|`.

Outer-only chord identity
-------------------------

For :math:`b > R_{\rm in}`, the bounce period is a single full chord
on the outer sphere with length

.. math::

   L_{\rm period}(b) = 2\sqrt{R_{\rm out}^2 - b^2}, \qquad
   \tau_{\rm period}(b) = \Sigma_t \cdot 2\sqrt{R_{\rm out}^2 - b^2}.

This is identical to a solid-sphere chord at the same outer radius and
impact parameter — the inner cavity is "invisible" to outer-only rays.
The closure is the rank-1 :math:`T = 1/(1 - \alpha_{\rm out}\,e^{
-\tau_{\rm period}})` form.

Three operator-level verifications
-----------------------------------

V_α1_hollow_sph. **Closed-shell bounce-sum self-consistency at
   corner** :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1`. With both
   surfaces fully reflective and constant volumetric source :math:`q`,
   BOTH the outer-only branch (rank-1) AND the through-ray branch
   (rank-2) MUST produce :math:`\psi(r, \mu) = q/\Sigma_t` everywhere.
   This is the structural composability check — the impact-parameter
   partition must compose cleanly with the rank-2 closure.

V_α2_hollow_sph. **Rank-2 resolvent algebra** specialised to the
   shell. Direct SymPy inversion of :math:`(I - S)` produces the
   canonical 2x2 form (same as slab asymmetric, with the curvilinear
   :math:`\tau_{\rm step}`). Reductions:

   - :math:`\alpha_{\rm in} = \alpha_{\rm out} = \alpha` (symmetric
     reflective): the rank-2 form. Special: at :math:`R_{\rm in} \to 0`
     and constant α_in = α_out = α, the through-ray subset has measure
     zero and the outer-only branch is the entire phase space —
     structural reduction to solid-sphere Variant α.
   - :math:`\alpha_{\rm in} = 1, \alpha_{\rm out} = 0`: reflective
     inner / vacuum outer.
   - :math:`\alpha_{\rm in} = 0, \alpha_{\rm out} = 1`: vacuum inner
     (cavity absorber) / reflective outer.
   - :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0`: vacuum-vacuum.

V_α3_hollow_sph. **Vacuum reduction at**
   :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0`. Both leading-:math:`\alpha`
   factors are zero; the outer-only :math:`T = 1` and the through-ray
   :math:`T = I`; the kernel reduces to bare first-leg integrals on
   both subsets. NO special-case branch needed — the BC absorption is
   fully encoded by the leading-:math:`\alpha` factors.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
- :file:`/.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 3C-1
  hollow sphere plan (interpretation (A) for the inner-surface BC).
- :mod:`.greens_function_slab_asymmetric` — Phase-3B asymmetric slab
  V_α derivations (the slab template lifted here to the curvilinear
  2-surface case).
- :mod:`.greens_function` — sphere V_α derivations (the rank-1
  outer-only template — re-used at :math:`b > R_{\rm in}`).
- :file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
  — frame match: rank-2 BIE block resolvent at 2-surface topologies.
"""
from __future__ import annotations

import sympy as sp


def derive_operator_constant_trial_closed_hollow_sphere() -> dict:
    r"""V_α1_hollow_sph — closed-shell bounce-sum self-consistency at
    corner :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1`.

    Verifies BOTH rank-1 (outer-only) AND rank-2 (through-ray) closure
    branches independently produce :math:`\psi(r, \mu) = q/\Sigma_t`
    on a constant volumetric source :math:`q` when both surfaces are
    fully reflective. The composability of the impact-parameter
    partition with the rank-2 closure is the load-bearing structural
    claim.

    **Outer-only branch** (:math:`b > R_{\rm in}`). Identical to V_α1
    sphere — rank-1 closure on a constant source. The full
    bounce-period chord :math:`L_p = 2\sqrt{R_{\rm out}^2 - b^2}` gives

    .. math::

       B(b; q) = \int_0^{L_p} q\,e^{-\Sigma_t s}\,\mathrm d s
               = \frac{q}{\Sigma_t}(1 - e^{-\Sigma_t L_p}),

    and the rank-1 surface flux at :math:`\alpha_{\rm out} = 1` is

    .. math::

       \psi_{\rm surf}(b) = \frac{B(b; q)}{1 - e^{-\Sigma_t L_p}}
                           = \frac{q}{\Sigma_t}.

    **Through-ray branch** (:math:`b \le R_{\rm in}`). Identical
    algebra to V_α1_slab_asym. The single-transit B integrals on a
    constant source :math:`q` are

    .. math::

       B_{\rm in}(b; q) = B_{\rm out}(b; q)
            = \int_0^{\tau_{\rm step}/\Sigma_t} q\,e^{-\Sigma_t s}\,\mathrm d s
            = \frac{q}{\Sigma_t}(1 - e^{-\tau_{\rm step}}),

    where :math:`\tau_{\rm step} = \Sigma_t \bigl(\sqrt{R_{\rm out}^2
    - b^2} - \sqrt{R_{\rm in}^2 - b^2}\bigr)`. The rank-2 closure at
    :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1` reduces (via the
    same algebraic identity :math:`(1 - e^{-\tau})(1 + e^{-\tau}) =
    1 - e^{-2\tau}`) to

    .. math::

       \psi_{\rm in}^{\rm out}(b) = \psi_{\rm out}^{\rm in}(b)
            = \frac{q}{\Sigma_t}.

    The interior reconstruction at any :math:`(r, \mu)` cancels first-
    leg attenuation against the surface flux — structurally identical
    to V_α1 on sphere/cylinder/slab.

    **Operator action on isotropic trial** :math:`\psi_{\rm trial} = 1`
    with :math:`q = \Sigma_s` (closed BC):

    .. math::

       (K \cdot 1)(r, \mu) = \frac{\Sigma_s}{\Sigma_t} = \omega_0,

    yielding :math:`k_{\rm eff} = k_\infty = \nu\Sigma_f/\Sigma_a`.

    Returns dict with the SymPy expressions and PASS flags.
    """
    Sigma_t, Sigma_s, q = sp.symbols(
        "Sigma_t Sigma_s q", positive=True, real=True,
    )
    R_in, R_out, b = sp.symbols("R_in R_out b", positive=True, real=True)
    alpha_in, alpha_out = sp.symbols(
        "alpha_in alpha_out", nonnegative=True, real=True,
    )

    # Outer-only branch: bounce-period chord = 2·sqrt(R_out^2 - b^2).
    L_period_outer_only = 2 * sp.sqrt(R_out**2 - b**2)
    tau_period = Sigma_t * L_period_outer_only

    # B(b; q) on a constant source q.
    B_outer_only = (q / Sigma_t) * (1 - sp.exp(-tau_period))

    # Rank-1 surface flux at α_out = 1 (closed outer).
    # T = 1/(1 - e^{-τ_period}); ψ_surf = α B T = B / (1 - e^{-τ_period})
    # = (q/Σ_t) · (1 - e^{-τ_period}) / (1 - e^{-τ_period}) = q/Σ_t.
    psi_surf_outer_closed = sp.simplify(
        B_outer_only / (1 - sp.exp(-tau_period))
    )
    pass_outer_closed_equals_q_over_sigma_t = (
        sp.simplify(psi_surf_outer_closed - q / Sigma_t) == 0
    )

    # Through-ray branch: τ_step = Σ_t · (sqrt(R_out^2 - b^2) - sqrt(R_in^2 - b^2)).
    tau_step = Sigma_t * (sp.sqrt(R_out**2 - b**2) - sp.sqrt(R_in**2 - b**2))

    # Single-transit B integrals on constant source q.
    B_in = (q / Sigma_t) * (1 - sp.exp(-tau_step))
    B_out = B_in  # by symmetry on a constant source

    # Rank-2 closure: same as slab asymmetric.
    e_step = sp.exp(-tau_step)
    det = 1 - alpha_in * alpha_out * e_step**2
    psi_in_out = (alpha_in * B_in
                  + alpha_in * e_step * alpha_out * B_out) / det
    psi_out_in = (alpha_out * e_step * alpha_in * B_in
                  + alpha_out * B_out) / det

    # Closed-BC corner: α_in = α_out = 1.
    psi_in_out_closed = sp.simplify(
        psi_in_out.subs([(alpha_in, 1), (alpha_out, 1)])
    )
    pass_psi_in_out_closed_equals_q_over_sigma_t = (
        sp.simplify(psi_in_out_closed - q / Sigma_t) == 0
    )

    psi_out_in_closed = sp.simplify(
        psi_out_in.subs([(alpha_in, 1), (alpha_out, 1)])
    )
    pass_psi_out_in_closed_equals_q_over_sigma_t = (
        sp.simplify(psi_out_in_closed - q / Sigma_t) == 0
    )

    # Operator action on isotropic trial ψ_trial = 1 → q = Σ_s, both
    # branches give Σ_s/Σ_t = ω_0.
    omega_0 = Sigma_s / Sigma_t
    K_outer = psi_surf_outer_closed.subs(q, Sigma_s)
    K_inner = psi_in_out_closed.subs(q, Sigma_s)
    pass_eigenvalue_outer = sp.simplify(K_outer - omega_0) == 0
    pass_eigenvalue_inner = sp.simplify(K_inner - omega_0) == 0

    # Leaky-mode sample: at α_in = 1/2, α_out = 1, b = R_in/2,
    # R_in = 1, R_out = 2, Σ_t = 1, q = 1: should give surface flux
    # strictly less than 1 = q/Σ_t.
    psi_sample = psi_in_out.subs([
        (alpha_in, sp.Rational(1, 2)),
        (alpha_out, 1),
        (R_in, 1), (R_out, 2),
        (b, sp.Rational(1, 2)),
        (Sigma_t, 1),
        (q, 1),
    ])
    psi_sample_simplified = sp.simplify(psi_sample)
    pass_leaky_below = bool(
        float(psi_sample_simplified.evalf()) < 1.0
    )

    return {
        "name": (
            "V_α1_hollow_sph: closed-shell bounce-sum self-consistency "
            "at α_in=α_out=1 (BOTH outer-only AND through-ray branches "
            "give q/Σ_t); leaky case strictly below"
        ),
        "tau_period_outer_only": tau_period,
        "tau_step": tau_step,
        "B_outer_only": B_outer_only,
        "B_in": B_in,
        "B_out": B_out,
        "psi_surf_outer_closed": psi_surf_outer_closed,
        "psi_in_out_closed": psi_in_out_closed,
        "psi_out_in_closed": psi_out_in_closed,
        "K_outer_on_constant_trial": K_outer,
        "K_inner_on_constant_trial": K_inner,
        "omega_0": omega_0,
        "psi_sample_leaky": psi_sample_simplified,
        "pass_outer_closed_equals_q_over_sigma_t":
            pass_outer_closed_equals_q_over_sigma_t,
        "pass_psi_in_out_closed_equals_q_over_sigma_t":
            pass_psi_in_out_closed_equals_q_over_sigma_t,
        "pass_psi_out_in_closed_equals_q_over_sigma_t":
            pass_psi_out_in_closed_equals_q_over_sigma_t,
        "pass_eigenvalue_outer": pass_eigenvalue_outer,
        "pass_eigenvalue_inner": pass_eigenvalue_inner,
        "pass_leaky_below": pass_leaky_below,
        "pass": (
            pass_outer_closed_equals_q_over_sigma_t
            and pass_psi_in_out_closed_equals_q_over_sigma_t
            and pass_psi_out_in_closed_equals_q_over_sigma_t
            and pass_eigenvalue_outer
            and pass_eigenvalue_inner
            and pass_leaky_below
        ),
    }


def derive_rank2_resolvent_hollow_sphere() -> dict:
    r"""V_α2_hollow_sph — rank-2 resolvent :math:`T = (I - S)^{-1}` for
    the hollow-sphere through-ray monodromy.

    Identical algebraic structure to V_α2_slab_asym — the rank-2 form
    is geometry-independent at the operator-symbol level. SymPy
    inverts :math:`(I - S)` directly with
    :math:`S = \mathrm{antidiag}(\alpha_{\rm in}\,e^{-\tau_{\rm step}},
    \alpha_{\rm out}\,e^{-\tau_{\rm step}})` and verifies:

    - :math:`\det(I - S) = 1 - \alpha_{\rm in}\,\alpha_{\rm out}\,
      e^{-2\tau_{\rm step}}`.
    - Canonical 2x2 form.
    - Reduction to symmetric :math:`\alpha_{\rm in} = \alpha_{\rm out}
      = \alpha` (the rank-2 → rank-1 collapse).
    - Reduction at :math:`\alpha_{\rm in} = 0` (vacuum inner): the
      cavity absorber case — through-rays are absorbed at the inner
      surface, no return.
    - Reduction at :math:`\alpha_{\rm out} = 0` (vacuum outer): the
      reflective-cavity case — through-rays escape at the outer
      surface, no return.
    - Reduction at :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0`:
      :math:`T = I` (vacuum-vacuum).

    The hollow-sphere-specific algebraic feature is that
    :math:`\tau_{\rm step}` is the SHELL-traversal length (curvilinear
    chord through the shell only — the cavity contributes no optical
    depth under interpretation (A)). At the algebraic level :math:`T`
    has the same form; only the chord meaning differs.

    Returns dict with the SymPy expressions and PASS flags.
    """
    alpha_in, alpha_out, tau = sp.symbols(
        "alpha_in alpha_out tau", nonnegative=True, real=True,
    )
    e_tau = sp.exp(-tau)

    S = sp.Matrix([[0, alpha_in * e_tau],
                   [alpha_out * e_tau, 0]])
    I2 = sp.eye(2)
    M = I2 - S
    T = M.inv()
    det_M = M.det()
    det_canonical = 1 - alpha_in * alpha_out * e_tau ** 2
    pass_det = sp.simplify(det_M - det_canonical) == 0

    T_canonical = sp.Matrix([
        [1, alpha_in * e_tau],
        [alpha_out * e_tau, 1],
    ]) / det_canonical
    pass_T_form = all(
        sp.simplify(T[i, j] - T_canonical[i, j]) == 0
        for i in range(2) for j in range(2)
    )

    # --- Reduction to symmetric α_in = α_out = α --------------------
    alpha_sym = sp.symbols("alpha", nonnegative=True, real=True)
    T_sym = T.subs([(alpha_in, alpha_sym), (alpha_out, alpha_sym)])

    # T_11 + T_12 = (1 + α·e^{-τ})/(1 - α²·e^{-2τ}) = 1/(1 - α·e^{-τ})
    sum_T11_T12 = sp.simplify(T_sym[0, 0] + T_sym[0, 1])
    expected_simplified = 1 / (1 - alpha_sym * e_tau)
    pass_symmetric_simplification = sp.simplify(
        sum_T11_T12 - expected_simplified
    ) == 0

    # --- Reduction at α_in = 0 (vacuum inner / cavity absorber) ------
    T_inner_vac = T.subs(alpha_in, 0)
    T_inner_vac_simp = sp.simplify(T_inner_vac)
    pass_alpha_in_zero = (
        sp.simplify(T_inner_vac_simp[0, 0] - 1) == 0
        and sp.simplify(T_inner_vac_simp[1, 1] - 1) == 0
        and sp.simplify(T_inner_vac_simp[0, 1] - 0) == 0
        and sp.simplify(T_inner_vac_simp[1, 0] - alpha_out * e_tau) == 0
    )

    # --- Reduction at α_out = 0 (vacuum outer / reflective cavity) ---
    T_outer_vac = T.subs(alpha_out, 0)
    T_outer_vac_simp = sp.simplify(T_outer_vac)
    pass_alpha_out_zero = (
        sp.simplify(T_outer_vac_simp[0, 0] - 1) == 0
        and sp.simplify(T_outer_vac_simp[1, 1] - 1) == 0
        and sp.simplify(T_outer_vac_simp[0, 1] - alpha_in * e_tau) == 0
        and sp.simplify(T_outer_vac_simp[1, 0] - 0) == 0
    )

    # --- Reduction at α_in = α_out = 0: T = I ------------------------
    T_vac_vac = T.subs([(alpha_in, 0), (alpha_out, 0)])
    pass_vacuum_identity = sp.simplify(T_vac_vac - I2) == sp.zeros(2, 2)

    return {
        "name": (
            "V_α2_hollow_sph: rank-2 resolvent T = (I-S)^{-1} via direct "
            "matrix inversion + reductions to symmetric / one-vacuum-"
            "surface / vacuum-vacuum corners"
        ),
        "S_monodromy": S,
        "T_resolvent": T,
        "T_canonical": T_canonical,
        "det_M": det_M,
        "det_canonical": det_canonical,
        "T_symmetric": T_sym,
        "T_inner_vac": T_inner_vac_simp,
        "T_outer_vac": T_outer_vac_simp,
        "T_vac_vac": T_vac_vac,
        "pass_det": pass_det,
        "pass_T_form": pass_T_form,
        "pass_symmetric_simplification": pass_symmetric_simplification,
        "pass_alpha_in_zero": pass_alpha_in_zero,
        "pass_alpha_out_zero": pass_alpha_out_zero,
        "pass_vacuum_identity": pass_vacuum_identity,
        "pass": (
            pass_det
            and pass_T_form
            and pass_symmetric_simplification
            and pass_alpha_in_zero
            and pass_alpha_out_zero
            and pass_vacuum_identity
        ),
    }


def derive_alpha_zero_kernel_reduction_hollow_sphere() -> dict:
    r"""V_α3_hollow_sph — at :math:`\alpha_{\rm in} = \alpha_{\rm out}
    = 0`, the rank-2 resolvent collapses to the identity AND the
    rank-1 outer-only resolvent collapses to 1; the kernel reduces to
    the bare first-leg integral on BOTH phase-space subsets.

    For the through-ray subset (:math:`b \le R_{\rm in}`): both
    leading factors :math:`\alpha_{\rm in}` (in :math:`\psi_{\rm in}^
    {\rm out}`) and :math:`\alpha_{\rm out}` (in :math:`\psi_{\rm out}^
    {\rm in}`) are zero, so :math:`\psi_{\rm surf} = (0, 0)^T` and the
    interior reconstruction is purely the first-leg integral.

    For the outer-only subset (:math:`b > R_{\rm in}`): the rank-1
    closure :math:`\psi_{\rm surf} = \alpha_{\rm out}\,B\,T = 0` at
    :math:`\alpha_{\rm out} = 0`, so the interior reconstruction is
    purely the first-leg integral.

    The vacuum-inner / reflective-outer (:math:`\alpha_{\rm in} = 0,
    \alpha_{\rm out} = 1`) corner is the **cavity absorber** case:
    through-rays are lost to the cavity (det → 1, ψ_in^out = 0,
    ψ_out^in = α_out · B_RL = B_RL — the right-wall outgoing flux
    receives only the direct one-transit-source contribution, no
    bouncing). Outer-only rays bounce normally at α_out = 1 between
    outer-and-outer.

    The reflective-inner / vacuum-outer (:math:`\alpha_{\rm in} = 1,
    \alpha_{\rm out} = 0`) corner is the **reflective cavity**:
    through-rays escape at the outer surface (det → 1, ψ_out^in = 0,
    ψ_in^out = α_in · B_LR = B_LR). Outer-only rays escape at outer
    too (α_out = 0).

    Returns dict with the SymPy expressions and PASS flag.
    """
    alpha_in, alpha_out, tau = sp.symbols(
        "alpha_in alpha_out tau", nonnegative=True, real=True,
    )
    B_in_sym, B_out_sym = sp.symbols("B_in B_out", real=True)

    e_tau = sp.exp(-tau)
    det = 1 - alpha_in * alpha_out * e_tau ** 2
    # Closure: B_in feeds ψ_in^out, B_out feeds ψ_out^in.
    psi_in_out = (alpha_in * B_in_sym
                  + alpha_in * e_tau * alpha_out * B_out_sym) / det
    psi_out_in = (alpha_out * e_tau * alpha_in * B_in_sym
                  + alpha_out * B_out_sym) / det

    # --- Vacuum-vacuum: α_in = α_out = 0 ---
    psi_in_out_at_zero = sp.simplify(
        psi_in_out.subs([(alpha_in, 0), (alpha_out, 0)])
    )
    psi_out_in_at_zero = sp.simplify(
        psi_out_in.subs([(alpha_in, 0), (alpha_out, 0)])
    )
    pass_substitution = (
        psi_in_out_at_zero == 0 and psi_out_in_at_zero == 0
    )

    psi_in_lim = sp.limit(
        sp.limit(psi_in_out, alpha_in, 0), alpha_out, 0,
    )
    psi_out_lim = sp.limit(
        sp.limit(psi_out_in, alpha_in, 0), alpha_out, 0,
    )
    pass_limit = (
        sp.simplify(psi_in_lim) == 0 and sp.simplify(psi_out_lim) == 0
    )

    # --- Vacuum-inner / reflective-outer: α_in = 0, α_out = 1 ---
    # Through-rays absorbed at inner surface (cavity absorber).
    # det → 1; ψ_in^out = 0; ψ_out^in = α_out · B_out = 1 · B_out = B_out.
    psi_in_at_inner_vac = sp.simplify(psi_in_out.subs(alpha_in, 0))
    psi_out_at_inner_vac = sp.simplify(psi_out_in.subs(alpha_in, 0))
    pass_cavity_absorber = (
        psi_in_at_inner_vac == 0
        and sp.simplify(psi_out_at_inner_vac - alpha_out * B_out_sym) == 0
    )

    # --- Reflective-inner / vacuum-outer: α_in = 1, α_out = 0 ---
    psi_in_at_outer_vac = sp.simplify(psi_in_out.subs(alpha_out, 0))
    psi_out_at_outer_vac = sp.simplify(psi_out_in.subs(alpha_out, 0))
    pass_reflective_cavity = (
        sp.simplify(psi_in_at_outer_vac - alpha_in * B_in_sym) == 0
        and psi_out_at_outer_vac == 0
    )

    return {
        "name": (
            "V_α3_hollow_sph: rank-2 closure vanishes at α_in=α_out=0; "
            "cavity-absorber and reflective-cavity reductions correct"
        ),
        "psi_in_out": psi_in_out,
        "psi_out_in": psi_out_in,
        "psi_in_out_at_zero": psi_in_out_at_zero,
        "psi_out_in_at_zero": psi_out_in_at_zero,
        "psi_in_lim": psi_in_lim,
        "psi_out_lim": psi_out_lim,
        "psi_in_at_inner_vac": psi_in_at_inner_vac,
        "psi_out_at_inner_vac": psi_out_at_inner_vac,
        "psi_in_at_outer_vac": psi_in_at_outer_vac,
        "psi_out_at_outer_vac": psi_out_at_outer_vac,
        "pass_substitution": pass_substitution,
        "pass_limit": pass_limit,
        "pass_cavity_absorber": pass_cavity_absorber,
        "pass_reflective_cavity": pass_reflective_cavity,
        "pass": (
            pass_substitution
            and pass_limit
            and pass_cavity_absorber
            and pass_reflective_cavity
        ),
    }
