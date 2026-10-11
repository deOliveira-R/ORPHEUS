r"""SymPy derivation — operator-level identities for the **annulus**
(hollow cylinder) Variant α Green's function reference (3D angle-resolved
:math:`(r, \mu_{\rm axial}, \varphi_{\rm az})`, with **independent**
specular reflectivities :math:`\alpha_{\rm in} \in [0, 1]` on the inner
cavity surface :math:`r = R_{\rm in}` and :math:`\alpha_{\rm out} \in
[0, 1]` on the outer surface :math:`r = R_{\rm out}`).

Phase-3C-2 extension of the Phase-3C-1 hollow sphere (see
:mod:`.greens_function_hollow_sphere`) and Phase-1 cylinder (see
:mod:`.greens_function_cylinder`). This is the **last two-surface
orbit-space class instance** in the Variant α plan; once it ships,
the 6-geometry × 2-orbit-space-class family (sphere, cylinder, slab,
slab-asym, hollow sphere, annulus) is complete on the unified
rank-1 / rank-2 framework (the retired ``variant_alpha_core``; the
characteristic reference's :mod:`~orpheus.derivations.continuous.characteristic.closure` is its successor). See Sphinx
§\ ``orbit-space-m-g-classification`` for the structural signature.

What is new vs hollow sphere
----------------------------

The annulus is the cylindrical analog of the hollow sphere. The
**closure algebra is byte-equal-shared** with hollow sphere — both
sit on the rank-2 boundary-to-boundary scattering resolvent
:math:`T = (I - S)^{-1}` with antidiagonal :math:`S(\alpha_{\rm in},
\alpha_{\rm out}, \tau_{\rm step})`. Only the **chord algebra**
differs:

- **Phase-space**: cylinder uses :math:`(r, \mu_{\rm axial},
  \varphi_{\rm az})` instead of sphere's :math:`(r, \mu)`. The two
  conserved quantities under specular reflection at radial-normal
  surfaces are :math:`\mu_{\rm axial}` (axial cosine — the bounce
  point's mirror plane is perpendicular to the cylinder axis) and
  :math:`b = r\,|\sin\varphi_{\rm az}|` (in-plane impact parameter —
  preserved across all bounces).

- **3D chord scaling**: the cylinder's in-plane chord
  :math:`L_{\rm 2D}(b) = \sqrt{R^2 - b^2} - \sqrt{R_{\rm in}^2 - b^2}`
  must be lifted to 3D arclength
  :math:`L_{\rm 3D}(b, \mu_{\rm axial}) = L_{\rm 2D}(b) / \sqrt{1 -
  \mu_{\rm axial}^2}` because the 3D ray inclines off the in-plane
  cross-section by the axial cosine. The optical depth is
  :math:`\tau_{\rm step}^{\rm annulus} = \Sigma_t \cdot L_{\rm 3D}`,
  i.e. the hollow-sphere shell-traversal length **plus** the
  :math:`1/\sqrt{1 - \mu_{\rm axial}^2}` axial correction.

Phase-space partition by impact parameter
------------------------------------------

For each interior phase-space point :math:`(r, \mu_{\rm axial},
\varphi_{\rm az})` with :math:`r \in [R_{\rm in}, R_{\rm out}]`,
the **conserved in-plane impact parameter** is
:math:`b = r\,|\sin\varphi_{\rm az}|`. Two cases:

- :math:`b > R_{\rm in}` — outer-only rays. The 2D in-plane
  projection ray does NOT intersect the inner-radius circle; the
  particle bounces between two points on the OUTER cylinder only.
  Closure rank-1 (orbit-space class: one-surface-compact),
  structurally identical to a solid-cylinder ray at the same outer
  radius and impact parameter.
- :math:`b \le R_{\rm in}` — through-rays. The 2D in-plane ray
  crosses the inner cavity circle. Under interpretation (A) (inner
  surface as specular reflector with reflectivity :math:`\alpha_{\rm
  in}`), the particle bounces alternately inner ↔ outer.
  Closure rank-2 (orbit-space class: two-surface — M/G interval
  with two distinct BC endpoints).

Through-ray rank-2 monodromy
----------------------------

For :math:`b \le R_{\rm in}`, the surface state is

.. math::

   \psi_{\rm surf}(b, \mu_{\rm axial}) = \begin{pmatrix}
       \psi_{\rm in}^{\rm out}(b, \mu_{\rm axial}) \\
       \psi_{\rm out}^{\rm in}(b, \mu_{\rm axial})
   \end{pmatrix}

and one step of the bouncing trajectory transports
:math:`\psi_{\rm in}^{\rm out}` to :math:`\psi_{\rm out}^{\rm in}`
(and vice versa) with **single-transit** optical depth

.. math::
   :label: peierls-greens-annulus-tau-step

   \tau_{\rm step}^{\rm annulus}(b, \mu_{\rm axial})
        = \Sigma_t \cdot
          \frac{\sqrt{R_{\rm out}^2 - b^2} - \sqrt{R_{\rm in}^2 - b^2}}
               {\sqrt{1 - \mu_{\rm axial}^2}}.

Compare hollow sphere :math:`\tau_{\rm step}^{\rm hollow\,sph}(b) =
\Sigma_t \cdot (\sqrt{R_{\rm out}^2 - b^2} - \sqrt{R_{\rm in}^2 -
b^2})` — the annulus inherits the hollow-sphere shell-traversal
chord and adds the cylinder axial correction
:math:`1/\sqrt{1 - \mu_{\rm axial}^2}`. The monodromy operator and
resolvent algebra are otherwise identical.

Outer-only chord identity
-------------------------

For :math:`b > R_{\rm in}`, the bounce period is a single full chord
on the outer cylinder cross-section in 2D, lifted to 3D:

.. math::

   L_{\rm period}^{\rm annulus}(b, \mu_{\rm axial})
        = \frac{2\sqrt{R_{\rm out}^2 - b^2}}{\sqrt{1 - \mu_{\rm axial}^2}},

which is the Phase-1 cylinder bounce-period chord identity (recall
:func:`derive_bounce_period_chord_cylinder` for the closed solid
cylinder). The closure is the rank-1
:math:`T = 1/(1 - \alpha_{\rm out}\,e^{-\Sigma_t L_{\rm period}})`
form.

Three operator-level verifications
-----------------------------------

V_α1_annulus. **Closed-shell bounce-sum self-consistency at corner**
   :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1`. With both surfaces
   fully reflective and constant volumetric source :math:`q`, BOTH
   the outer-only branch (rank-1) AND the through-ray branch (rank-2)
   MUST produce :math:`\psi(r, \mu_{\rm axial}, \varphi_{\rm az}) =
   q/\Sigma_t` everywhere. This is the structural composability check
   on the cylinder's 3D angular phase-space — the impact-parameter
   partition must compose cleanly with the rank-2 closure AND survive
   the axial-correction lift.

V_α2_annulus. **Rank-2 resolvent algebra** specialised to the annulus
   shell. SymPy verifies that the rank-2 inversion
   :math:`T = (I - S)^{-1}` produces the canonical 2x2 form when
   :math:`\tau_{\rm step}` is the annulus 3D optical depth (Eq.
   :eq:`peierls-greens-annulus-tau-step`). The algebraic structure is
   the same as hollow sphere (the rank-2 resolvent is geometry-
   independent at the operator-symbol level — only the chord meaning
   differs), but we re-verify in the cylinder's 3D phase-space to lock
   the correspondence at the symbol level.

V_α3_annulus. **Vacuum reduction at**
   :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0`. Both leading-α
   factors zero the closure on each surface component; the kernel
   reduces to bare first-leg integrals on both phase-space subsets
   (rank-1 outer-only AND rank-2 through-ray) with no special-case
   branch needed.

References
----------

- Sanchez, R. (1986). *Transp. Theor. Stat. Phys.* 14.
- :file:`/.claude/plans/archive/peierls-greens-cylinder-and-2bc.md` — Phase 3C-2
  annulus plan.
- :mod:`.greens_function_hollow_sphere` — Phase-3C-1 hollow sphere V_α
  derivations (the rank-2 + impact-parameter template lifted here to
  the cylinder cross-section).
- :mod:`.greens_function_cylinder` — Phase-1 cylinder V_α derivations
  (the cylinder 3D angular framework + rank-1 outer-only template the
  outer-only branch reuses).
- :mod:`.greens_function_slab_asymmetric` — Phase-3B asymmetric slab
  V_α derivations (the rank-2 closure operator the through-ray branch
  reuses).
- :file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
  — frame match: rank-2 BIE block resolvent at 2-surface topologies.
- Issue #129 — cylinder Bickley-Naylor planar-limit physics, motivating
  the angle-resolved (NOT axial-pre-integrated) discipline maintained
  here.
"""
from __future__ import annotations

import sympy as sp


def derive_operator_constant_trial_closed_annulus() -> dict:
    r"""V_α1_annulus — closed-shell bounce-sum self-consistency at
    corner :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1`.

    Verifies BOTH rank-1 (outer-only) AND rank-2 (through-ray) closure
    branches independently produce :math:`\psi(r, \mu_{\rm axial},
    \varphi_{\rm az}) = q/\Sigma_t` on a constant volumetric source
    :math:`q` when both surfaces are fully reflective. The composability
    of the impact-parameter partition with the rank-2 closure on the
    cylinder's 3D angular phase-space is the load-bearing structural
    claim.

    The algebraic identity is the same as hollow sphere — only the
    chord-meaning differs (we replace the hollow-sphere shell chord
    :math:`\Delta = \sqrt{R_{\rm out}^2 - b^2} - \sqrt{R_{\rm in}^2 -
    b^2}` with its 3D-lifted form
    :math:`\Delta_{\rm 3D} = \Delta / \sqrt{1 - \mu_{\rm axial}^2}`).

    **Outer-only branch** (:math:`b > R_{\rm in}`). Rank-1 closure on
    a constant source. The full bounce-period chord :math:`L_p^{\rm 3D}
    = 2\sqrt{R_{\rm out}^2 - b^2} / \sqrt{1 - \mu_{\rm axial}^2}` gives

    .. math::

       B(b, \mu_{\rm axial}; q) = \int_0^{L_p^{\rm 3D}}
            q\,e^{-\Sigma_t s}\,\mathrm d s
            = \frac{q}{\Sigma_t}(1 - e^{-\Sigma_t L_p^{\rm 3D}}),

    and the rank-1 surface flux at :math:`\alpha_{\rm out} = 1` is

    .. math::

       \psi_{\rm surf}(b, \mu_{\rm axial})
            = \frac{B}{1 - e^{-\Sigma_t L_p^{\rm 3D}}}
            = \frac{q}{\Sigma_t}.

    **Through-ray branch** (:math:`b \le R_{\rm in}`). Rank-2 closure
    with single-transit 3D-lifted optical depth
    :math:`\tau_{\rm step} = \Sigma_t (\sqrt{R_{\rm out}^2 - b^2}
    - \sqrt{R_{\rm in}^2 - b^2})/\sqrt{1 - \mu_{\rm axial}^2}`. The
    same algebraic identity :math:`(1 - e^{-\tau})(1 + e^{-\tau}) =
    1 - e^{-2\tau}` collapses both surface fluxes to :math:`q/\Sigma_t`
    at :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1`.

    **Operator action on isotropic trial** :math:`\psi_{\rm trial} = 1`
    with :math:`q = \Sigma_s` (closed BC):

    .. math::

       (K \cdot 1)(r, \mu_{\rm axial}, \varphi_{\rm az})
            = \frac{\Sigma_s}{\Sigma_t} = \omega_0,

    yielding :math:`k_{\rm eff} = k_\infty = \nu\Sigma_f/\Sigma_a`.

    Returns dict with the SymPy expressions and PASS flags.
    """
    Sigma_t, Sigma_s, q = sp.symbols(
        "Sigma_t Sigma_s q", positive=True, real=True,
    )
    R_in, R_out, b = sp.symbols("R_in R_out b", positive=True, real=True)
    mu_ax = sp.symbols("mu_axial", real=True)
    alpha_in, alpha_out = sp.symbols(
        "alpha_in alpha_out", nonnegative=True, real=True,
    )

    # In-plane factor s_ip = sqrt(1 - μ_axial²); 3D arclength = (2D arc) / s_ip.
    s_ip = sp.sqrt(1 - mu_ax ** 2)

    # ===================== Outer-only branch (rank-1) =====================
    # Bounce-period 3D chord (in-plane antipodal lifted by 1/s_ip).
    L_period_3D_outer_only = 2 * sp.sqrt(R_out ** 2 - b ** 2) / s_ip
    tau_period = Sigma_t * L_period_3D_outer_only

    # B(b, μ_axial; q) on a constant source q.
    B_outer_only = (q / Sigma_t) * (1 - sp.exp(-tau_period))

    # Rank-1 surface flux at α_out = 1 (closed outer).
    psi_surf_outer_closed = sp.simplify(
        B_outer_only / (1 - sp.exp(-tau_period))
    )
    pass_outer_closed_equals_q_over_sigma_t = (
        sp.simplify(psi_surf_outer_closed - q / Sigma_t) == 0
    )

    # ===================== Through-ray branch (rank-2) =====================
    # Single-transit 3D shell-traversal length.
    L_step_3D = (
        sp.sqrt(R_out ** 2 - b ** 2) - sp.sqrt(R_in ** 2 - b ** 2)
    ) / s_ip
    tau_step = Sigma_t * L_step_3D

    # Single-transit B integrals on a constant source q (symmetric).
    B_in = (q / Sigma_t) * (1 - sp.exp(-tau_step))
    B_out = B_in

    # Rank-2 closure: same algebraic form as hollow sphere /
    # slab-asymmetric — only τ_step differs.
    e_step = sp.exp(-tau_step)
    det = 1 - alpha_in * alpha_out * e_step ** 2
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

    # Operator action on isotropic trial ψ_trial = 1 → q = Σ_s.
    omega_0 = Sigma_s / Sigma_t
    K_outer = psi_surf_outer_closed.subs(q, Sigma_s)
    K_inner = psi_in_out_closed.subs(q, Sigma_s)
    pass_eigenvalue_outer = sp.simplify(K_outer - omega_0) == 0
    pass_eigenvalue_inner = sp.simplify(K_inner - omega_0) == 0

    # Leaky-mode sample: at α_in = 1/2, α_out = 1, b = R_in/2,
    # R_in = 1, R_out = 2, μ_axial = 1/3, Σ_t = 1, q = 1: should be
    # strictly less than 1 = q/Σ_t.
    psi_sample = psi_in_out.subs([
        (alpha_in, sp.Rational(1, 2)),
        (alpha_out, 1),
        (R_in, 1), (R_out, 2),
        (b, sp.Rational(1, 2)),
        (mu_ax, sp.Rational(1, 3)),
        (Sigma_t, 1),
        (q, 1),
    ])
    psi_sample_simplified = sp.simplify(psi_sample)
    pass_leaky_below = bool(
        float(psi_sample_simplified.evalf()) < 1.0
    )

    return {
        "name": (
            "V_α1_annulus: closed-shell bounce-sum self-consistency at "
            "α_in=α_out=1 (BOTH outer-only AND through-ray branches give "
            "q/Σ_t on the cylinder 3D angular phase-space); leaky case "
            "strictly below"
        ),
        "tau_period_outer_only": tau_period,
        "tau_step": tau_step,
        "L_period_3D_outer_only": L_period_3D_outer_only,
        "L_step_3D": L_step_3D,
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


def derive_3d_chord_scaling_annulus() -> dict:
    r"""V_α2_annulus.aux — algebraic verification that the annulus
    3D chord scales as the hollow-sphere shell chord divided by
    :math:`\sqrt{1 - \mu_{\rm axial}^2}`.

    Verifies:

    .. math::

       \tau_{\rm step}^{\rm annulus}(b, \mu_{\rm axial})
            = \frac{\tau_{\rm step}^{\rm hollow\,sph}(b)}
                   {\sqrt{1 - \mu_{\rm axial}^2}}

    where :math:`\tau_{\rm step}^{\rm hollow\,sph}(b) = \Sigma_t
    (\sqrt{R_{\rm out}^2 - b^2} - \sqrt{R_{\rm in}^2 - b^2})`.

    This is the load-bearing chord-algebra identity that takes the
    Phase-3C-1 hollow sphere chord meaning to the cylinder 3D angular
    phase-space (Issue #129 angle-resolved discipline). It also
    verifies the in-plane → 3D arclength conversion for the outer-only
    branch:

    .. math::

       L_{\rm period}^{\rm 3D, annulus}(b, \mu_{\rm axial})
            = \frac{L_{\rm period}^{\rm 2D, cyl}(b)}
                   {\sqrt{1 - \mu_{\rm axial}^2}}.

    Returns dict with the SymPy expressions and PASS flag.
    """
    Sigma_t = sp.symbols("Sigma_t", positive=True, real=True)
    R_in, R_out, b = sp.symbols("R_in R_out b", positive=True, real=True)
    mu_ax = sp.symbols("mu_axial", real=True)
    s_ip = sp.sqrt(1 - mu_ax ** 2)

    # Hollow-sphere shell chord (no axial dependence — sphere has
    # μ ≡ μ_axial and chord = √(R²-b²) − √(R_in²-b²) period of half).
    L_step_hollow_sph = sp.sqrt(R_out ** 2 - b ** 2) - sp.sqrt(R_in ** 2 - b ** 2)
    tau_step_hollow_sph = Sigma_t * L_step_hollow_sph

    # Annulus 3D shell-traversal chord (this same quantity in 2D, lifted
    # to 3D arclength by the axial correction).
    L_step_3D_annulus = (
        sp.sqrt(R_out ** 2 - b ** 2) - sp.sqrt(R_in ** 2 - b ** 2)
    ) / s_ip
    tau_step_annulus = Sigma_t * L_step_3D_annulus

    # Identity: τ_annulus / τ_hollow_sph = 1 / s_ip.
    ratio = sp.simplify(tau_step_annulus / tau_step_hollow_sph)
    pass_step_ratio = sp.simplify(ratio - 1 / s_ip) == 0

    # Outer-only bounce period: 2D cyl is 2·sqrt(R²-b²); 3D is divided
    # by s_ip.
    L_period_2D_cyl = 2 * sp.sqrt(R_out ** 2 - b ** 2)
    L_period_3D_annulus = L_period_2D_cyl / s_ip
    pass_period_ratio = sp.simplify(
        L_period_3D_annulus * s_ip - L_period_2D_cyl
    ) == 0

    # Sanity: at μ_axial = 0 (in-plane motion), s_ip = 1, so 3D chord
    # = 2D chord.
    L_step_in_plane = sp.simplify(L_step_3D_annulus.subs(mu_ax, 0))
    pass_in_plane_match = sp.simplify(
        L_step_in_plane - L_step_hollow_sph
    ) == 0

    return {
        "name": (
            "V_α2_annulus.aux: 3D chord scaling identity τ_annulus = "
            "τ_hollow_sph / sqrt(1 - μ_axial²)"
        ),
        "tau_step_hollow_sph": tau_step_hollow_sph,
        "tau_step_annulus": tau_step_annulus,
        "L_period_2D_cyl": L_period_2D_cyl,
        "L_period_3D_annulus": L_period_3D_annulus,
        "ratio": ratio,
        "pass_step_ratio": pass_step_ratio,
        "pass_period_ratio": pass_period_ratio,
        "pass_in_plane_match": pass_in_plane_match,
        "pass": (
            pass_step_ratio
            and pass_period_ratio
            and pass_in_plane_match
        ),
    }
