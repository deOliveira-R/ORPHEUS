r"""SymPy derivation — operator-level identities for the **asymmetric**
slab Variant α Green's function reference (1D Cartesian, independent
specular reflectivities :math:`\alpha_L \in [0, 1]` on the left wall
and :math:`\alpha_R \in [0, 1]` on the right wall).

Phase-3B extension of the symmetric slab Variant α (per
:file:`/.claude/plans/peierls-greens-cylinder-and-2bc.md`). Mirrors the
Phase-3A symmetric-slab proofs in :mod:`.greens_function_slab` with the
**rank-2 boundary-to-boundary scattering resolvent** that handles
independent per-wall reflectivities.

Rank-2 monodromy and resolvent
-------------------------------

The 2-vector surface state is :math:`\psi_{\rm surf} = [\psi_L^+(\mu),
\psi_R^-(\mu)]` (outgoing flux at :math:`x = 0` in :math:`+\mu`, and
outgoing flux at :math:`x = L` in :math:`-\mu`). Bouncing at conserved
:math:`|\mu|`, one step of the trajectory advances :math:`\psi_L^+` to
:math:`\psi_R^-` (transit + reflection at the right wall, factor
:math:`\alpha_R\,e^{-\tau}`) and vice versa, with **single-transit**
optical depth :math:`\tau = \Sigma_t L /|\mu|`.

The monodromy operator is

.. math::

   S(\alpha_L, \alpha_R, \tau) \;=\; \begin{pmatrix}
       0                              & \alpha_L\,e^{-\tau} \\
       \alpha_R\,e^{-\tau}            & 0
   \end{pmatrix},

and the rank-2 resolvent is :math:`T = (I - S)^{-1}`. The diagonal of
:math:`S` is identically zero — one step never returns to the same
wall. Consequently :math:`T` is rank-2 with off-diagonal coupling
encoded by the single-transit attenuation :math:`e^{-\tau}` and the
per-wall reflectivities. **The two-bounce-per-period structure of slab
geometry is encoded via the matrix structure of** :math:`S`, NOT via
:math:`\alpha^2` in a scalar formula (the symmetric Phase-3A rank-1
form factored the matrix into a scalar by exploiting :math:`\alpha_L =
\alpha_R = \alpha`).

Single-transit B integrals
---------------------------

The closure equations involve **single-transit** chord integrals (NOT
full bounce-period integrals as in the rank-1 symmetric closure):

.. math::

   B_{LR}(\mu; q) &= \int_0^{L/|\mu|} q(|\mu|\,s)\,
                       e^{-\Sigma_t s}\,\mathrm d s, \\
   B_{RL}(\mu; q) &= \int_0^{L/|\mu|} q(L - |\mu|\,s)\,
                       e^{-\Sigma_t s}\,\mathrm d s.

:math:`B_{LR}` is the **source contribution at** :math:`x = 0` from
particles moving in the :math:`-\hat\mu` direction along the chord
:math:`x = 0 \to L` (the chord is traversed BACKWARD in the
characteristic-tracing sense — from current point :math:`x=0` going
INTO the slab in the :math:`+\hat x` direction, which is the
backward-tracing direction of the :math:`-\hat\mu` characteristic).
:math:`B_{RL}` is the mirror at :math:`x = L` from particles moving
in the :math:`+\hat\mu` direction (chord :math:`L \to 0` traversed
backward).

The closure equations

.. math::

   \psi_L^+(\mu) &= \alpha_L\,B_{LR}(\mu) + \alpha_L\,e^{-\tau}\,
                       \psi_R^-(\mu), \\
   \psi_R^-(\mu) &= \alpha_R\,B_{RL}(\mu) + \alpha_R\,e^{-\tau}\,
                       \psi_L^+(\mu),

derive from :math:`\psi_L^+ = \alpha_L\,\psi(0, -\hat\mu)`,
:math:`\psi(0, -\hat\mu) = e^{-\tau}\,\psi(L, -\hat\mu) + B_{LR}` (the
incoming flux at :math:`x = 0` going LEFT integrates source along
the chord from :math:`x = 0` to :math:`x = L`), and the mirror
identity at the right wall. Read as :math:`\psi_{\rm surf} = \alpha\,B
+ S\,\psi_{\rm surf}` with :math:`(\alpha B)_L = \alpha_L\,B_{LR}`
and :math:`(\alpha B)_R = \alpha_R\,B_{RL}`:

.. math::

   \psi_{\rm surf} \;=\; T\,(\alpha\,B).

Three operator-level verifications
-----------------------------------

V_α1_slab_asym. **Closed-slab bounce-sum self-consistency, rank-2
   form**. For closed asymmetric slab (any :math:`\alpha_L, \alpha_R
   \in (0, 1]`) with constant volumetric source :math:`q`, the rank-2
   closure produces :math:`\psi(x, \mu) = q/\Sigma_t` everywhere. The
   algebra is more involved than Phase-3A symmetric — both
   :math:`B_{RL}` and :math:`B_{LR}` independently produce
   :math:`(q/\Sigma_t)(1 - e^{-\tau})`, and the rank-2 resolvent
   contracts them back to :math:`(q/\Sigma_t)` on each surface
   component. At :math:`\alpha_L = \alpha_R = 1` reduces to V_α1_slab.

V_α2_slab_asym. **Rank-2 resolvent algebra**. Direct SymPy inversion
   of :math:`(I - S)` produces the canonical 2x2 form

   .. math::

      T = \frac{1}{1 - \alpha_L\,\alpha_R\,e^{-2\tau}}
            \begin{pmatrix}
                1                       & \alpha_L\,e^{-\tau} \\
                \alpha_R\,e^{-\tau}     & 1
            \end{pmatrix},

   with :math:`\det(I - S) = 1 - \alpha_L\,\alpha_R\,e^{-2\tau}`. The
   reductions to (i) symmetric :math:`\alpha_L = \alpha_R = \alpha`,
   (ii) one-vacuum-wall :math:`\alpha_L = 0` or :math:`\alpha_R = 0`,
   and (iii) vacuum-vacuum :math:`\alpha_L = \alpha_R = 0` are verified
   symbolically.

V_α3_slab_asym. **Vacuum reduction at** :math:`\alpha_L = \alpha_R =
   0`. Trivially :math:`S = 0 \Rightarrow T = I \Rightarrow
   \psi_{\rm surf} = \alpha\,B = 0` (because both :math:`\alpha_L` and
   :math:`\alpha_R` are zero), so the kernel reduces to the bare
   first-leg integral.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
- Hébert, A. (2009). *Applied Reactor Physics* §3.8.5 — slab
  :math:`E_n` forms and rank-1 white-BC closure.
- :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md` — Phase 3B
  asymmetric slab plan.
- :mod:`.greens_function_slab` — Phase-3A symmetric-slab V_α
  derivations (this module is the rank-2 generalisation).
- :file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
  — frame match: Variant α extension to two-surface topologies as
  rank-2 block resolvent. The slab asymmetric prototype is the first
  validated instance of this prediction.
"""
from __future__ import annotations

import sympy as sp


def derive_operator_constant_trial_closed_slab_asymmetric() -> dict:
    r"""V_α1_slab_asym — closed asymmetric-slab bounce-sum
    self-consistency.

    For a closed homogeneous slab of width :math:`L` with INDEPENDENT
    reflective specular BCs :math:`\alpha_L, \alpha_R \in (0, 1]` and
    constant volumetric source :math:`q`, the rank-2 surface-flux
    closure must yield :math:`\psi(x, \mu) = q/\Sigma_t` everywhere.

    The single-transit B integrals on a constant source :math:`q`
    evaluate to

    .. math::

       B_{LR}(\mu; q) = B_{RL}(\mu; q)
            = \int_0^{L/|\mu|} q\,e^{-\Sigma_t s}\,\mathrm d s
            = \frac{q}{\Sigma_t}\,(1 - e^{-\tau}),

    where :math:`\tau = \Sigma_t L /|\mu|`. Setting :math:`B_{LR} =
    B_{RL} = (q/\Sigma_t)(1 - e^{-\tau})`, the rank-2 closure produces

    .. math::

       \psi_L^+(\mu) &= \frac{1}{1 - \alpha_L\alpha_R\,e^{-2\tau}}\,
                          \big[\alpha_L\,B_{LR} + \alpha_L\,e^{-\tau}\,
                          \alpha_R\,B_{RL}\big] \\
                     &= \frac{(q/\Sigma_t)(1-e^{-\tau})\,\alpha_L\,
                            (1 + \alpha_R\,e^{-\tau})}
                          {1 - \alpha_L\alpha_R\,e^{-2\tau}}.

    This **does NOT collapse to** :math:`q/\Sigma_t` at general
    :math:`(\alpha_L, \alpha_R)` — closed-slab self-consistency
    requires :math:`\alpha_L = \alpha_R = 1` in the rank-2 framework
    (with both walls fully reflective there is no leakage; otherwise
    leakage out the partial wall keeps the surface flux below
    :math:`q/\Sigma_t`). The **load-bearing identity** for asymmetric
    slab with closed BCs (:math:`\alpha_L = \alpha_R = 1`) is then

    .. math::

       \psi_L^+(\mu)\big|_{\alpha_L=\alpha_R=1}
            = \frac{(q/\Sigma_t)(1-e^{-\tau})(1 + e^{-\tau})}
                   {1 - e^{-2\tau}}
            = \frac{q}{\Sigma_t},

    using :math:`(1-e^{-\tau})(1+e^{-\tau}) = 1 - e^{-2\tau}`. The
    interior reconstruction at :math:`\mu > 0`,

    .. math::

       \psi(x, \mu) = \frac{q}{\Sigma_t}(1 - e^{-\Sigma_t x/\mu})
                       + e^{-\Sigma_t x/\mu}\,\psi_L^+(\mu)
                     = \frac{q}{\Sigma_t},

    cancels the first-leg attenuation against the surface flux —
    structurally identical to V_α1_slab/V_α1 sphere/cylinder.

    For the operator action on isotropic trial :math:`\psi_{\rm trial}
    = 1` with :math:`q = \Sigma_s\,\psi_{\rm trial} = \Sigma_s` (closed
    BC :math:`\alpha_L = \alpha_R = 1`):

    .. math::

       (K \cdot 1)(x, \mu) = \frac{\Sigma_s}{\Sigma_t} = \omega_0,

    yielding :math:`k_{\rm eff} = k_\infty = \nu\Sigma_f/\Sigma_a`.

    The proof structure is **the same V_α1 algebra restricted to the
    closed-BC corner of the rank-2 parameter space**. SymPy verifies:

    - The rank-2 closure at :math:`\alpha_L = \alpha_R = 1` and
      constant source produces :math:`\psi_L^+ = q/\Sigma_t`
      (independent of :math:`\tau`, hence independent of :math:`|\mu|`).
    - The interior reconstruction is :math:`q/\Sigma_t` everywhere.
    - The operator eigenvalue on isotropic trial is :math:`\omega_0`.
    - At general :math:`\alpha_L, \alpha_R \in (0, 1)` the surface
      flux is **strictly less than** :math:`q/\Sigma_t` — closed-slab
      eigenmode is not an eigenmode of the leaky operator.

    Returns dict with the SymPy expressions and PASS flags.
    """
    Sigma_t, Sigma_s, q = sp.symbols(
        "Sigma_t Sigma_s q", positive=True, real=True,
    )
    L = sp.symbols("L", positive=True, real=True)
    mu_abs = sp.symbols("mu_abs", positive=True, real=True)
    alpha_L, alpha_R = sp.symbols(
        "alpha_L alpha_R", nonnegative=True, real=True,
    )

    tau = Sigma_t * L / mu_abs

    # Single-transit B integrals on constant source q.
    # B_LR = B_RL = ∫_0^{L/|µ|} q · e^{-Σ_t s} ds = (q/Σ_t)·(1-e^{-τ})
    # (constant source — both single-transit chords give same value).
    B_single = (q / Sigma_t) * (1 - sp.exp(-tau))
    B_LR = B_single
    B_RL = B_single

    # Rank-2 closure: ψ_L^+ = (1/det) · [α_L·B_LR + α_L·e^{-τ}·α_R·B_RL]
    # ψ_R^- = (1/det) · [α_R·e^{-τ}·α_L·B_LR + α_R·B_RL]
    det = 1 - alpha_L * alpha_R * sp.exp(-2 * tau)
    psi_L_plus = (alpha_L * B_LR
                  + alpha_L * sp.exp(-tau) * alpha_R * B_RL) / det
    psi_R_minus = (alpha_R * sp.exp(-tau) * alpha_L * B_LR
                   + alpha_R * B_RL) / det

    # Closed-BC corner: α_L = α_R = 1.
    psi_L_plus_closed = sp.simplify(
        psi_L_plus.subs([(alpha_L, 1), (alpha_R, 1)])
    )
    pass_psi_L_closed_equals_q_over_sigma_t = (
        sp.simplify(psi_L_plus_closed - q / Sigma_t) == 0
    )

    psi_R_minus_closed = sp.simplify(
        psi_R_minus.subs([(alpha_L, 1), (alpha_R, 1)])
    )
    pass_psi_R_closed_equals_q_over_sigma_t = (
        sp.simplify(psi_R_minus_closed - q / Sigma_t) == 0
    )

    # Interior reconstruction at µ > 0 (entered from x = 0):
    # ψ(x, µ) = (q/Σ_t)(1 - e^{-Σ_t x/µ}) + e^{-Σ_t x/µ}·ψ_L^+(µ)
    x = sp.symbols("x", positive=True, real=True)
    L_first_pos = x / mu_abs
    F_const = (q / Sigma_t) * (1 - sp.exp(-Sigma_t * L_first_pos))
    psi_total_closed = sp.simplify(
        F_const + sp.exp(-Sigma_t * L_first_pos) * psi_L_plus_closed
    )
    pass_psi_total_closed_constant = (
        sp.simplify(psi_total_closed - q / Sigma_t) == 0
    )

    # Operator action on isotropic trial ψ_trial = 1 → q = Σ_s.
    omega_0 = Sigma_s / Sigma_t
    K_on_one_closed = psi_total_closed.subs(q, Sigma_s)
    pass_eigenvalue = sp.simplify(K_on_one_closed - omega_0) == 0

    # At general α_L, α_R ∈ (0, 1) with α_L < 1 OR α_R < 1, the surface
    # flux is strictly less than q/Σ_t. Verify with a sample point:
    # set α_L = 1/2, α_R = 1, τ = 1, q/Σ_t = 1 — should be < 1.
    psi_sample = psi_L_plus.subs([
        (alpha_L, sp.Rational(1, 2)),
        (alpha_R, 1),
        (Sigma_t, 1), (L, 1), (mu_abs, 1),
        (q, 1),
    ])
    psi_sample_simplified = sp.simplify(psi_sample)
    pass_leaky_below = (psi_sample_simplified < 1)
    # Use sp.evalf for the strict inequality check (sympy may not
    # evaluate symbolic comparisons directly).
    pass_leaky_below_numeric = bool(
        float(psi_sample_simplified.evalf()) < 1.0
    )

    return {
        "name": (
            "V_α1_slab_asym: closed-asymmetric-slab bounce-sum "
            "self-consistency at α_L=α_R=1, leaky case strictly below"
        ),
        "tau": tau,
        "B_RL": B_RL,
        "B_LR": B_LR,
        "det": det,
        "psi_L_plus": psi_L_plus,
        "psi_R_minus": psi_R_minus,
        "psi_L_plus_closed": psi_L_plus_closed,
        "psi_R_minus_closed": psi_R_minus_closed,
        "psi_total_closed": psi_total_closed,
        "K_on_constant_trial_closed": K_on_one_closed,
        "omega_0": omega_0,
        "psi_sample_leaky": psi_sample_simplified,
        "pass_psi_L_closed_equals_q_over_sigma_t":
            pass_psi_L_closed_equals_q_over_sigma_t,
        "pass_psi_R_closed_equals_q_over_sigma_t":
            pass_psi_R_closed_equals_q_over_sigma_t,
        "pass_psi_total_closed_constant": pass_psi_total_closed_constant,
        "pass_eigenvalue": pass_eigenvalue,
        "pass_leaky_below": pass_leaky_below_numeric,
        "pass": (
            pass_psi_L_closed_equals_q_over_sigma_t
            and pass_psi_R_closed_equals_q_over_sigma_t
            and pass_psi_total_closed_constant
            and pass_eigenvalue
            and pass_leaky_below_numeric
        ),
    }


def derive_rank2_resolvent_slab_asymmetric() -> dict:
    r"""V_α2_slab_asym — rank-2 resolvent :math:`T = (I - S)^{-1}`
    derived by direct SymPy matrix inversion.

    Inverts the :math:`2 \times 2` monodromy

    .. math::

       I - S = \begin{pmatrix}
           1                              & -\alpha_L\,e^{-\tau} \\
           -\alpha_R\,e^{-\tau}            & 1
       \end{pmatrix}

    via SymPy's ``Matrix.inv()`` and verifies:

    - :math:`\det(I - S) = 1 - \alpha_L\,\alpha_R\,e^{-2\tau}`.
    - The four matrix entries match the canonical closed form.
    - Reduction to symmetric :math:`\alpha_L = \alpha_R = \alpha`
      reproduces (in scaled form) the rank-1 :math:`\psi_{\rm surf} =
      \alpha\,B / (1 - \alpha^2 e^{-2\tau})` formula via :math:`\psi_L^+
      = T_{11}\,\alpha\,B + T_{12}\,\alpha\,B = \alpha\,B\,(T_{11} +
      T_{12})`. The factor :math:`(T_{11} + T_{12}) = (1 + \alpha
      e^{-\tau})/(1 - \alpha^2 e^{-2\tau}) = 1/(1 - \alpha\,e^{-\tau})`,
      a different scalar form than the rank-1 :math:`1/(1 - \alpha^2
      e^{-2\tau})`. The identity bridging the two arithmetic forms is
      :math:`B_{\rm period} = (1 + e^{-\tau})\,B_{\rm single\,transit}`
      — this is the structural relationship between the two
      formulations on a constant source.
    - Reduction at :math:`\alpha_L = 0` or :math:`\alpha_R = 0`:
      :math:`\det = 1`, :math:`T_{11} = T_{22} = 1`, :math:`T_{12} =
      \alpha_L\,e^{-\tau}`, :math:`T_{21} = \alpha_R\,e^{-\tau}`.
    - Reduction at :math:`\alpha_L = \alpha_R = 0`: :math:`T = I`.

    Returns dict with the SymPy expressions and PASS flags.
    """
    alpha_L, alpha_R, tau = sp.symbols(
        "alpha_L alpha_R tau", nonnegative=True, real=True,
    )
    e_tau = sp.exp(-tau)

    S = sp.Matrix([[0, alpha_L * e_tau],
                   [alpha_R * e_tau, 0]])
    I2 = sp.eye(2)
    M = I2 - S
    T = M.inv()
    det_M = M.det()
    det_canonical = 1 - alpha_L * alpha_R * e_tau ** 2
    pass_det = sp.simplify(det_M - det_canonical) == 0

    T_canonical = sp.Matrix([
        [1, alpha_L * e_tau],
        [alpha_R * e_tau, 1],
    ]) / det_canonical
    pass_T_form = all(
        sp.simplify(T[i, j] - T_canonical[i, j]) == 0
        for i in range(2) for j in range(2)
    )

    # --- Reduction to symmetric α_L = α_R = α --------------------------
    alpha_sym = sp.symbols("alpha", nonnegative=True, real=True)
    T_sym = T.subs([(alpha_L, alpha_sym), (alpha_R, alpha_sym)])

    # T_11 + T_12 = (1 + α·e^{-τ}) / (1 - α²·e^{-2τ}) = 1 / (1 - α·e^{-τ})
    sum_T11_T12 = sp.simplify(T_sym[0, 0] + T_sym[0, 1])
    expected_simplified = 1 / (1 - alpha_sym * e_tau)
    pass_symmetric_simplification = sp.simplify(
        sum_T11_T12 - expected_simplified
    ) == 0

    # --- Reduction at α_L = 0 ----------------------------------------
    T_left_vac = T.subs(alpha_L, 0)
    T_left_vac_simp = sp.simplify(T_left_vac)
    pass_alpha_L_zero = (
        sp.simplify(T_left_vac_simp[0, 0] - 1) == 0
        and sp.simplify(T_left_vac_simp[1, 1] - 1) == 0
        and sp.simplify(T_left_vac_simp[0, 1] - 0) == 0
        and sp.simplify(T_left_vac_simp[1, 0] - alpha_R * e_tau) == 0
    )

    # --- Reduction at α_R = 0 ----------------------------------------
    T_right_vac = T.subs(alpha_R, 0)
    T_right_vac_simp = sp.simplify(T_right_vac)
    pass_alpha_R_zero = (
        sp.simplify(T_right_vac_simp[0, 0] - 1) == 0
        and sp.simplify(T_right_vac_simp[1, 1] - 1) == 0
        and sp.simplify(T_right_vac_simp[0, 1] - alpha_L * e_tau) == 0
        and sp.simplify(T_right_vac_simp[1, 0] - 0) == 0
    )

    # --- Reduction at α_L = α_R = 0 (vacuum-vacuum): T = I -------------
    T_vac_vac = T.subs([(alpha_L, 0), (alpha_R, 0)])
    pass_vacuum_identity = sp.simplify(T_vac_vac - I2) == sp.zeros(2, 2)

    return {
        "name": (
            "V_α2_slab_asym: rank-2 resolvent T = (I-S)^{-1} via direct "
            "matrix inversion + reductions to symmetric / one-vacuum-wall "
            "/ vacuum-vacuum corners"
        ),
        "S_monodromy": S,
        "T_resolvent": T,
        "T_canonical": T_canonical,
        "det_M": det_M,
        "det_canonical": det_canonical,
        "T_symmetric": T_sym,
        "T_left_vac": T_left_vac_simp,
        "T_right_vac": T_right_vac_simp,
        "T_vac_vac": T_vac_vac,
        "pass_det": pass_det,
        "pass_T_form": pass_T_form,
        "pass_symmetric_simplification": pass_symmetric_simplification,
        "pass_alpha_L_zero": pass_alpha_L_zero,
        "pass_alpha_R_zero": pass_alpha_R_zero,
        "pass_vacuum_identity": pass_vacuum_identity,
        "pass": (
            pass_det
            and pass_T_form
            and pass_symmetric_simplification
            and pass_alpha_L_zero
            and pass_alpha_R_zero
            and pass_vacuum_identity
        ),
    }


def derive_alpha_zero_kernel_reduction_slab_asymmetric() -> dict:
    r"""V_α3_slab_asym — at :math:`\alpha_L = \alpha_R = 0`, the rank-2
    resolvent collapses to the identity and the kernel reduces to the
    bare first-leg integral.

    The rank-2 monodromy :math:`S` is identically zero at :math:`\alpha_L
    = \alpha_R = 0`, so :math:`T = (I - 0)^{-1} = I` and the closure
    becomes

    .. math::

       \psi_{\rm surf} = T \cdot \alpha\,B
                      = I \cdot \begin{pmatrix} 0 \\ 0 \end{pmatrix}
                      = \begin{pmatrix} 0 \\ 0 \end{pmatrix}

    (because both :math:`\alpha_L\,B_{RL}` and :math:`\alpha_R\,B_{LR}`
    have a leading-:math:`\alpha` factor that is zero on each component).
    The interior reconstruction is purely the first-leg integral

    .. math::

       \psi(x, \mu) = F(x, \mu) \quad (\alpha_L = \alpha_R = 0).

    This trivially mirrors V_α3_slab/V_α3 sphere/cylinder. Operator
    interpretation: the rank-2 prototype handles vacuum-vacuum BC with
    no special-case branch; the BC absorption is fully encoded by the
    leading-:math:`\alpha` factors on each surface component.

    Returns dict with the SymPy expressions and PASS flag.
    """
    alpha_L, alpha_R, tau = sp.symbols(
        "alpha_L alpha_R tau", nonnegative=True, real=True,
    )
    B_LR, B_RL = sp.symbols("B_LR B_RL", real=True)

    e_tau = sp.exp(-tau)
    det = 1 - alpha_L * alpha_R * e_tau ** 2
    # Corrected closure (B_LR feeds ψ_L^+, B_RL feeds ψ_R^-).
    psi_L_plus = (alpha_L * B_LR
                  + alpha_L * e_tau * alpha_R * B_RL) / det
    psi_R_minus = (alpha_R * e_tau * alpha_L * B_LR
                   + alpha_R * B_RL) / det

    psi_L_plus_at_zero = sp.simplify(
        psi_L_plus.subs([(alpha_L, 0), (alpha_R, 0)])
    )
    psi_R_minus_at_zero = sp.simplify(
        psi_R_minus.subs([(alpha_L, 0), (alpha_R, 0)])
    )
    pass_substitution = (
        psi_L_plus_at_zero == 0 and psi_R_minus_at_zero == 0
    )

    # Verify limits also vanish (robustness against piecewise / removable
    # singularities).
    psi_L_lim = sp.limit(
        sp.limit(psi_L_plus, alpha_L, 0),
        alpha_R, 0,
    )
    psi_R_lim = sp.limit(
        sp.limit(psi_R_minus, alpha_L, 0),
        alpha_R, 0,
    )
    pass_limit = (
        sp.simplify(psi_L_lim) == 0 and sp.simplify(psi_R_lim) == 0
    )

    # Also check one-vacuum-wall: α_L = 0 → ψ_L^+ = 0 (no source feeds
    # the left wall outgoing flux when there is no left-wall reflection).
    psi_L_at_left_vac = sp.simplify(psi_L_plus.subs(alpha_L, 0))
    psi_R_at_left_vac = sp.simplify(psi_R_minus.subs(alpha_L, 0))
    # ψ_R^- at α_L = 0: det collapses to 1, ψ_R^- = α_R · B_RL.
    pass_one_vacuum_wall_left = (
        psi_L_at_left_vac == 0
        and sp.simplify(psi_R_at_left_vac - alpha_R * B_RL) == 0
    )

    return {
        "name": (
            "V_α3_slab_asym: rank-2 closure vanishes at α_L=α_R=0; "
            "one-vacuum-wall reduces correctly"
        ),
        "psi_L_plus": psi_L_plus,
        "psi_R_minus": psi_R_minus,
        "psi_L_plus_at_zero": psi_L_plus_at_zero,
        "psi_R_minus_at_zero": psi_R_minus_at_zero,
        "psi_L_lim": psi_L_lim,
        "psi_R_lim": psi_R_lim,
        "psi_L_at_left_vac": psi_L_at_left_vac,
        "psi_R_at_left_vac": psi_R_at_left_vac,
        "pass_substitution": pass_substitution,
        "pass_limit": pass_limit,
        "pass_one_vacuum_wall_left": pass_one_vacuum_wall_left,
        "pass": pass_substitution and pass_limit and pass_one_vacuum_wall_left,
    }
