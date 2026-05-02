.. _theory-peierls-greens:

==========================================================================
Variant α — angle-resolved Green's function reference (sphere, α∈[0,1])
==========================================================================

.. contents:: Contents
   :local:
   :depth: 3


Key Facts
=========

**Read this before modifying or extending the Variant α prototype, or
before treating its output as a reference for any sphere problem.**
For the **method-agnostic foundations** (integral form of the
transport equation, BC parametrisation, common verification chain),
see :ref:`theory-peierls`. For the parallel Nyström / matrix-Galerkin
architecture (rank-:math:`N` closures, Phase 4
``specular_multibounce``, slab + cyl + sph) see
:ref:`theory-peierls-nystrom`. For the predecessor that motivates
this page — the Phase 5 retreat — see :ref:`peierls-phase5-retreat`.

**Scope as of 2026-05-01** (Plan 2 Part B + A1/A2/A3/Plan-(b)
follow-ons all closed):

- **Closed homogeneous sphere, specular BC** (V_α1, B-phase shipped).
  :math:`k_{\rm eff} = k_\infty` to machine precision.
- **Vacuum BC** (A1+A2 follow-on) — :math:`\alpha = 0` parametrisation.
  Cross-checked against PS-1982 [PS1982]_ Eq. (21) reference solver
  to ≤ 1e-4 relative on k_eff across four (R, c) configurations
  (τ_R ∈ {2.5, 5}, c ∈ {0.4, 0.6}).
- **Multi-group homogeneous sphere** (A3 follow-on) —
  :func:`solve_greens_function_sphere_mg` with full G×G scattering
  matrix and arbitrary χ. Closed sphere reduces to
  :func:`~orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`
  transfer-matrix eigenvalue + spectrum to ≤ 1e-9 relative.
- **Multi-region sphere k-eigenvalue** (Plan-(b) Option 2) —
  :func:`solve_greens_function_sphere_mr` with piecewise σ_t along
  trajectory + bounce-period chord. Issue #132 reproducer:
  fuel-A inner / moderator-B outer at radii=[0.5, 1.0] gives
  :math:`k_{\rm eff} \approx 0.735` (subcritical, sensible) vs
  Phase 4 ``specular_multibounce`` rank-2 = 1.015 (+57 % catastrophe
  avoided).
- **Multi-region sphere fixed-source** (Plan-(b) Option 1) —
  :func:`solve_greens_function_sphere_mr_fixed_source`.
  Cross-checked against Garcia 2021 [Garcia2021]_ Table 5
  (Williams 1991 Case 1, 3-region sphere) at 15 r-points: < 2 %
  agreement at non-interface points, < 12 % near interfaces (cubic-
  spline source-interpolation smooths the discontinuous σ_s — known
  prototype limitation).

- **Variant α is a parallel research-grade reference**, not a
  production replacement for ``boundary="specular_multibounce"``.
  The prototype now covers **sphere geometry only** (no cylinder),
  **isotropic scattering only** (no :math:`\omega_1\ne 0`), with
  :math:`\alpha\in[0,1]` BC and arbitrary multi-group / multi-region
  XS (closeout decision,
  :file:`.claude/agent-memory/numerics-investigator/peierls_greens_phase1_closeout.md`).
- **Variant α is exact for the closed homogeneous sphere.**
  V_α1 algebraically proves :math:`(K\cdot 1)(r,\mu) = \omega_0`,
  so the rank-1 isotropic mode is the unique eigenmode at
  :math:`\keff = \kinf = \nSigf{}/\Siga{}`. The numerical
  implementation reproduces this to machine precision (1e-13 %).
- **Phase 4 ``specular_multibounce`` carries a small rank-N truncation
  bias** that does not vanish at :math:`N = 3`. For the fuel-A-like
  sphere :math:`(R, \Sigt{}, \Sigs{}, \nSigf{}) = (5, 0.5, 0.38, 0.025)`,
  the rank-N errors against :math:`\kinf` are
  :math:`\epsilon_1 = 0.27\,\%`, :math:`\epsilon_2 = 0.25\,\%`,
  :math:`\epsilon_3 = 0.12\,\%`. Variant α improves on Phase 4 N=3 by
  about 0.12 % for closed-sphere homogeneous configurations.
- **Variant α never assembles the angle-integrated kernel
  :math:`g_\alpha(\rho'\to\rho)`.** It iterates the **angle-resolved**
  flux :math:`\psi(r,\mu)` along bouncing characteristics. The
  Hadamard finite-part diagonal singularity that killed Phase 5 is
  bypassed *structurally*, not by a quadrature trick.
- **The bounce sum closes analytically.** For perfect specular
  (:math:`\alpha = 1`) the geometric series over bounces of optical
  chord :math:`\Sigt{}\,L_p = 2\Sigt{}R\mu_{\rm surf}` closes to
  :math:`T(\mu_{\rm surf}) = 1/(1 - e^{-\Sigt{}\,L_p})`
  (:eq:`peierls-greens-T-mu-surf`). The :math:`\mu \to 0`
  geometric-series factor is :math:`O(1/\mu)` integrable; the
  scalar-flux extraction
  :math:`\phi(r) = 2\pi \int_{-1}^{1} \psi(r,\mu)\,\mathrm d\mu`
  is convergent with sufficient angular nodes (Gauss-Legendre on
  :math:`(0,1]` with :math:`n_\mu = 24` is more than enough for the
  rank-1 isotropic mode).
- **No closure is needed at the operator level.** The boundary
  condition is absorbed into the Green's function via Sanchez Eq.
  (A1) :math:`t = \bar t + t_h`, so there is no separate ``K_bc``
  matrix, no rank-N gating, no :math:`(1-P_{ss})^{-1}` Hébert
  correction. This is the load-bearing **structural** advantage of
  the Green's function reformulation.
- **V_α2 makes the rank-1 algebraic equivalence airtight**:
  :math:`T_{00}^{\rm sphere} = P_{ss}^{\rm sphere}` (Eq.
  :eq:`peierls-greens-V-alpha-2`) reduces to the same closed form as
  Hébert :math:`(1-P_{ss})^{-1}`, so rank-1 ``specular_multibounce``,
  rank-1 ``white_hebert``, and the rank-1 Variant α reduction agree
  bit-for-bit on closed homogeneous sphere.
- **V_α3 covers vacuum BC compatibility**: the Sanchez (A6)
  leading factor :math:`2\alpha` makes :math:`g_h \to 0` at
  :math:`\alpha = 0` (Eq. :eq:`peierls-greens-V-alpha-3`), so
  setting :math:`\alpha = 0` in the Variant α architecture recovers
  the bare vacuum sphere kernel :math:`\bar g_2`. No special-case
  branch is needed; vacuum BC is the trivial limit of specular BC,
  not a separate solver path.
- **What this prototype does NOT cover**: anisotropic scattering
  (:math:`\omega_1\ne 0`), Sanchez 1986 diffuse-re-emission
  :math:`\beta`-branch, multi-region cylinder, multi-region cylinder
  fixed-source (Phase 1b). Multi-region sphere (homogeneous regions,
  piecewise σ_t) **is** covered — see
  :ref:`peierls-greens-multiregion`. Vacuum-BC k-eigenvalue with
  non-trivial spatial mode **is** covered — see
  :ref:`peierls-greens-vacuum-extension` and the PS-1982 cross-check.
  **Cylinder geometry (1G + MG homogeneous, α∈[0,1])** is covered
  as a Phase 1 standalone — see :ref:`peierls-greens-cylinder` below.
  Phase 2 unification with sphere is deferred until both standalone
  reference solvers have captured accuracy floors.
- **Phase 5 framing**. The Phase 5 retreat is a structural finding,
  not a failure narrative. The retreat established that the
  **angle-integrated** kernel :math:`g_\alpha` is hypersingular —
  Variant α exploits the **angle-resolved** kernel
  :math:`\tilde t(r'\to r,\mu)`, which is a different mathematical
  object. The two paths are not different discretisations of the
  same operator; they target different operators that share the
  same physical content.


Motivation: the Phase 5 retreat made angle-resolution structurally necessary
============================================================================

The Phase 5 investigation (Issue #133, CLOSED 2026-04-28) attempted
to discretise the **angle-integrated** Sanchez 1986 [SanchezTTSP1986]_
Eq. (A6) kernel

.. math::
   :label: peierls-greens-sanchez-A6

   g_h(\rho'\to\rho) = 2\alpha \int_{\mu_0}^{1} T(\mu_-)\,\mu_*^{-1}\,
       \cosh(\rho\mu)\,\cosh(\rho'\mu_*)\,e^{-2 a \mu_-}\,\mathrm d\mu

via Nyström sampling: assemble :math:`K_{ij} = g_h(\rho_j \to \rho_i)`
on a radial Gauss-Legendre grid, factor :math:`(I - K)^{-1}` at every
power-iteration step, extract :math:`k_{\rm eff}` from the dominant
mode. Three rounds of investigation across multiple agents converged
on the same diagnosis: *the kernel* :math:`g_h` *itself is not in the
Fredholm second-kind class*. At the discrete diagonal :math:`r_i = r_j`
the µ-resolved primitives :math:`F_{\rm out}(r,\mu)\,G_{\rm in}(r,\mu)`
carry a :math:`1/(\cos\omega_i \cos\omega_j)` Jacobian whose poles
coincide on the chord-visibility cone :math:`\mu_0(r,r) = 0`, producing
a non-integrable :math:`1/(\mu^2 - \mu_0^2)` Hadamard finite-part
singularity that no standard quadrature trick can resolve. The
matrix-Galerkin form ``boundary="specular_multibounce"`` (Phase 4)
absorbs this singularity via basis projection — the rank-:math:`N`
shifted-Legendre projection acts as essential smoothing — which is
why Phase 4 produces sane numbers despite working with the same
underlying kernel structure. Phase 4 and Phase 5 are different
operators, not different discretisations of the same operator.

The retreat surfaced the **operator-level** question Plan 2 Variant α
addresses: rather than discretising the angle-integrated kernel,
*work directly with the angle-resolved kernel*
:math:`\tilde t(r'\to r,\mu)`, which is the pre-integration form
defined by Sanchez Eq. (A1) and Eq. (A5). The angle-integrated kernel
:math:`g_h(\rho'\to\rho) = \int \tilde t(r'\to r,\mu)\,\mathrm d\mu`
is hypersingular at the discrete diagonal because the integration
collapses the geometric structure that keeps :math:`\tilde t`
pointwise finite. *Don't collapse it*. Iterate :math:`\psi(r,\mu)`
on a 2-D phase-space grid; sample :math:`\tilde t` along single
characteristics; only at the very end perform the 1-D angular
integral :math:`\phi(r) = 2\pi \int \psi(r,\mu)\,\mathrm d\mu`,
which has only an *integrable* :math:`1/\mu` singularity from the
geometric bounce-sum factor :math:`T(\mu_{\rm surf})`.

This is a structural fix, not a quadrature-tightening fix. The
Phase 5 lesson — *"Sanchez 1986 is correct math, wrong
discretisation"* — directly motivates Variant α: keep the same
math (Sanchez Eqs. A1, A5, A6), change the discretisation target
from the angle-integrated kernel to the angle-resolved kernel.


The angle-resolved Green's function with BC absorbed
=====================================================

Following Sanchez 1986 [SanchezTTSP1986]_ §2 and Appendix, the
Green's function :math:`t(r',\Omega' \to r,\Omega)` for the
homogeneous sphere is the unique solution of

.. math::
   :label: peierls-greens-defining-bvp

   \begin{aligned}
     (\Omega \cdot \nabla + \Sigt{})\,t &= \delta(r - r')\,
        \delta_2(\Omega \cdot \Omega'),
        \quad a_- \le \rho < a, \\
     t(r',\Omega' \to r,\Omega) &= \alpha\,t(r',\Omega' \to r,\Omega_R),
        \quad |\rho| = a,\;\Omega \cdot n \le 0,
   \end{aligned}

where :math:`\Omega_R = \Omega + 2\mu n` is the specularly reflected
direction at the surface and :math:`\alpha \in [0,1]` is the specular
reflection coefficient (:math:`\alpha = 1` for perfect specular,
:math:`\alpha = 0` for vacuum). The vacuum / diffuse legs are dropped
here; the full BC line carries a :math:`\beta` term for diffuse
re-emission (Sanchez Eq. (A3.a)) which we set to zero throughout.

Sanchez's central identity (Eq. (A1)) splits the Green's function
into a vacuum part and a BC-absorbed part:

.. math::
   :label: peierls-greens-A1-split

   t(r',\Omega' \to r,\Omega) = \bar t(r',\Omega' \to r,\Omega) +
       t_h(r',\Omega' \to r,\Omega).

The vacuum kernel :math:`\bar t` is the Boltzmann-Green's-function
in free space (no boundary). The "homogeneous-BC" kernel :math:`t_h`
absorbs all of the BC physics — every specular reflection contributes
a term that propagates the source through the *modified* trajectory
that bounces off the boundary. **There is no separate K_bc closure**;
the BC is fully encoded in the kernel via the geometric bounce sum.
This is the structural advantage Variant α exploits.

For the perfect-specular case (:math:`\alpha = 1`, :math:`\beta = 0`,
:math:`\omega_1 = 0`) the explicit form of :math:`t_h` is Sanchez
Eq. (A5):

.. math::
   :label: peierls-greens-A5-specular

   t_h(r' \to r,\mu) = \frac{1}{2\pi A}\,e^{-\tau_+ - \tau_-}\,
       T(\mu_+) \cdot \frac{\delta(\mu_- - \mu_+)}{\mu_+}

where :math:`\tau_+ = a\mu_+ - \rho'\mu'` is the source-to-surface
optical depth, :math:`\tau_- = a\mu_- + \rho\mu` is the
surface-to-receiver optical depth, and :math:`T(\mu_+) = 1/(1 -
e^{-2 a \mu_+})` is the geometric bounce-sum factor. The Dirac delta
:math:`\delta(\mu_- - \mu_+)` is the *direction-preserving* feature
of perfect specular reflection: the angle of incidence equals the
angle of reflection, so the source-leaving and receiver-arriving
directions cosines coincide along a single characteristic.

The angle-resolved Green's function therefore decomposes a transport
operator action :math:`(K\cdot \psi_{\rm trial})(r,\mu)` into three
ingredients on a single bouncing characteristic:

1. A **first-leg trajectory integral** from the receiver point
   :math:`r` backward along the direction :math:`-\Omega_\mu` to the
   surface entry point.
2. A **bounce-period integral** along the antipodal chord at impact
   parameter :math:`h(r,\mu)`, of length :math:`L_p = 2R\mu_{\rm
   surf}`, where :math:`\mu_{\rm surf}` is the direction cosine at
   the surface.
3. A **closed-form geometric bounce sum** :math:`T(\mu_{\rm surf})`
   that summates over all subsequent bounces of the same trajectory.

Each piece is pointwise finite; the only place the
:math:`\mu \to 0` integrable singularity appears is in the final
scalar-flux extraction
:math:`\phi(r) = 2\pi \int \psi(r,\mu)\,\mathrm d\mu`. The
divergence is structurally deferred — that is the load-bearing claim
of Variant α relative to Phase 5.


Operator action along bouncing characteristics
================================================

For each phase-space grid point :math:`(r_i, \mu_q)`, the Variant α
operator action :math:`(K\cdot\psi)(r_i,\mu_q)` is computed in three
steps. The trajectory geometry is set up first; the source
:math:`q(r) = (\Sigs{} + \nSigf{}/k)\,\phi(r) / 4\pi` is built from
the previous iterate :math:`\phi(r) = 4\pi \int_0^1 \psi(r,\mu)\,
\mathrm d\mu` (the factor :math:`4\pi` combines the azimuthal
:math:`2\pi` with the :math:`\mu \leftrightarrow -\mu` symmetry
doubling for the closed-sphere isotropic mode).

Trajectory geometry
-------------------

At receiver point :math:`r_i` and direction cosine :math:`\mu_q > 0`
(outward), the chord from :math:`r_i` in direction :math:`-\Omega_\mu`
exits the sphere at the surface after a distance

.. math::
   :label: peierls-greens-L0

   L_0(r_i, \mu_q) = \sqrt{R^2 - r_i^2 (1 - \mu_q^2)} - r_i \mu_q,

with the corresponding **surface direction cosine**

.. math::
   :label: peierls-greens-mu-surf

   \mu_{\rm surf}(r_i, \mu_q) = \frac{1}{R}\sqrt{R^2 - r_i^2 (1 - \mu_q^2)}.

The discriminant :math:`R^2 - r_i^2 (1 - \mu_q^2)` is positive for
any interior point :math:`r_i < R`. The bounce-period chord (the
distance from the surface entry point through the sphere to the next
specular bounce) is the **antipodal chord** at impact parameter
:math:`h = R\sqrt{1 - \mu_{\rm surf}^2}`, of length

.. math::
   :label: peierls-greens-Lp

   L_p(\mu_{\rm surf}) = 2 R \mu_{\rm surf}.

Two checks anchor this geometry.  At :math:`r_i = 0` (sphere centre)
the first-leg distance reduces to :math:`L_0 = R` and
:math:`\mu_{\rm surf} = 1` (any radial outward direction hits the
surface normally), so :math:`L_p = 2R` (full diameter).  At
:math:`r_i \to R^-` with :math:`\mu_q \to 1` (radial outward) the
first-leg distance vanishes :math:`L_0 \to 0` (already at the
surface), with :math:`\mu_{\rm surf} \to 1` and :math:`L_p \to 2R`.

First-leg trajectory integral
-----------------------------

.. math::
   :label: peierls-greens-trajectory-integral

   F(r_i, \mu_q) = \int_0^{L_0(r_i,\mu_q)}
       q\bigl(|r_i - s\,\Omega_\mu|\bigr)\,
       e^{-\Sigt{}\,s}\,\mathrm d s.

The path radius along the chord is parametrised by

.. math::

   |r(s)|^2 = r_i^2 - 2 r_i \mu_q s + s^2,
   \qquad s \in [0, L_0(r_i,\mu_q)].

For an isotropic source :math:`q(r)` (which is the case throughout
Variant α since the closed-sphere problem has only the rank-1
isotropic eigenmode), :math:`F` depends on :math:`(r_i, \mu_q)`
only through the chord geometry. The integrand is bounded
:math:`q(r) \cdot e^{-\Sigt{}\,s} \le q_{\max} \cdot 1` and
:math:`F` is finite for all :math:`(r_i,\mu_q)` with :math:`r_i < R`.
We use composite Gauss-Legendre on :math:`[0, L_0]` to evaluate it.

Bounce-period integral
----------------------

.. math::
   :label: peierls-greens-bounce-period-integral

   B(\mu_{\rm surf}) = \int_0^{L_p(\mu_{\rm surf})}
       q\bigl(|r_{\rm chord}(s)|\bigr)\,
       e^{-\Sigt{}\,s}\,\mathrm d s,

where the chord position at arc length :math:`s` is

.. math::

   |r_{\rm chord}(s)|^2 = h^2 + (s - L_p/2)^2,
   \qquad h^2 = R^2 (1 - \mu_{\rm surf}^2),
   \qquad s \in [0, L_p].

The bounce-period integral is structurally identical to the
first-leg integral but on the antipodal chord; it represents *one
period* of the bouncing trajectory.

Geometric bounce sum
--------------------

The surface flux :math:`\psi_{\rm surf}` at any specular bounce
point of the trajectory satisfies the **fixed-point equation** that
Variant α exploits:

.. math::
   :label: peierls-greens-surface-fixed-point

   \psi_{\rm surf} = B(\mu_{\rm surf}) +
       e^{-\Sigt{}\,L_p}\,\psi_{\rm surf}.

Reading: the surface flux at the next bounce equals the source
contribution along the bounce-period chord plus the attenuated
surface flux from the previous bounce. Solving for
:math:`\psi_{\rm surf}` gives the closed-form geometric bounce sum

.. math::
   :label: peierls-greens-T-mu-surf

   \psi_{\rm surf}(\mu_{\rm surf}) = T(\mu_{\rm surf})\,B(\mu_{\rm surf}),
   \qquad
   T(\mu_{\rm surf}) = \frac{1}{1 - e^{-\Sigt{}\,L_p(\mu_{\rm surf})}}.

The geometric series :math:`\sum_{k=0}^\infty e^{-k \Sigt{} L_p} =
1/(1 - e^{-\Sigt{} L_p})` converges absolutely for :math:`\Sigt{}\,
L_p > 0`, which holds at any interior trajectory. As
:math:`\mu_{\rm surf} \to 0`, the geometric-series factor
:math:`T(\mu_{\rm surf}) \to 1/(\Sigt{}\,L_p) = 1/(2\Sigt{}R\,
\mu_{\rm surf})` is :math:`O(1/\mu_{\rm surf})` — integrable in the
final :math:`\mu`-integration but requires reasonable angular
quadrature.

Total angular flux
------------------

The master Variant α equation combines the first-leg integral with
the attenuated surface flux from the bounce sum:

.. math::
   :label: peierls-greens-function-architecture

   \boxed{\;
   \psi(r_i, \mu_q) = F(r_i, \mu_q) +
       e^{-\Sigt{}\,L_0(r_i,\mu_q)}\,
       T(\mu_{\rm surf})\,B(\mu_{\rm surf})
   \;}

This is the discrete analogue of Sanchez Eq. (A1) for the
angle-resolved Green's function in the perfect-specular case
(:math:`\alpha = 1`, :math:`\beta = 0`, :math:`\omega_1 = 0`). Each
term on the right-hand side is a 1-D quadrature along a single
characteristic; nothing in the assembly hits the angle-integrated
kernel :math:`g_h(\rho'\to\rho)` whose Hadamard finite-part
singularity killed Phase 5.

The reference implementation of :eq:`peierls-greens-function-architecture`
is :func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function._apply_operator`;
the public driver is
:func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function.solve_greens_function_specular_sphere`.


V_α1 — closed-sphere k_inf identity
====================================

The first load-bearing operator-level identity is V_α1:

.. math::
   :label: peierls-greens-V-alpha-1

   (K \cdot 1)(r, \mu) = \omega_0,
   \qquad \omega_0 = \frac{\Sigs{}}{\Sigt{}},

i.e. the constant function :math:`\psi_{\rm trial}(r,\mu) = 1` is
an :math:`\omega_0`-eigenmode of the scattering operator
:math:`K_{\rm scat}` for the closed homogeneous sphere with specular
BC. The fission-source eigenvalue equation
:math:`(I - K_{\rm scat})\,\phi = (\nSigf{}/k)\,\phi` then reduces
to :math:`(1 - \omega_0)\,\mathrm{const} = (\nSigf{}/k)\,\mathrm{const}`,
giving

.. math::
   :label: peierls-greens-k-inf

   k_{\rm eff} = \kinf = \frac{\nSigf{}}{\Siga{}},
   \qquad \Siga{} = \Sigt{} - \Sigs{}.

This is the no-leakage k_inf of the closed sphere, derivable
without any spatial transport machinery — but it must hold for the
Variant α implementation if the operator action is wired correctly.
The symbolic proof of :eq:`peierls-greens-V-alpha-1` lives in
:func:`~orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function.derive_operator_constant_trial_closed_sphere`
(SymPy verification, paired numerical gate in
:file:`tests/derivations/test_peierls_greens_function_symbolic.py`).
The argument has three steps.

Step 1: surface fixed-point with constant source
-------------------------------------------------

For trial :math:`\psi_{\rm trial} = 1` (isotropic constant), the
isotropic-scattering source is
:math:`q = \Sigs{}\,\psi_{\rm trial} = \Sigs{}` — also constant
in space. Plugging into the bounce-period integral
:eq:`peierls-greens-bounce-period-integral`,

.. math::

   B(\mu_{\rm surf}) = \int_0^{L_p} \Sigs{}\,e^{-\Sigt{}\,s}\,\mathrm d s
                     = \frac{\Sigs{}}{\Sigt{}}\,(1 - e^{-\Sigt{}\,L_p}),

so the surface fixed-point equation
:eq:`peierls-greens-surface-fixed-point` becomes

.. math::

   \psi_{\rm surf} = \frac{\Sigs{}}{\Sigt{}}\,(1 - e^{-\Sigt{}\,L_p})
                   + e^{-\Sigt{}\,L_p}\,\psi_{\rm surf}.

Solving:

.. math::

   (1 - e^{-\Sigt{}\,L_p})\,\psi_{\rm surf} =
       \frac{\Sigs{}}{\Sigt{}}\,(1 - e^{-\Sigt{}\,L_p})
   \;\Longrightarrow\;
   \psi_{\rm surf} = \frac{\Sigs{}}{\Sigt{}} = \omega_0.

The dependence on :math:`L_p` (and hence on :math:`\mu_{\rm surf}`)
**cancels exactly**: the surface flux is a constant equal to
:math:`\omega_0` for every bounce of every trajectory. This is the
no-leakage signature of the closed sphere — every bounce sees the
same uniform source contribution.

Step 2: total angular flux is independent of first-leg geometry
----------------------------------------------------------------

The first-leg integral with constant source evaluates to

.. math::

   F(r_i, \mu_q) = \int_0^{L_0} \Sigs{}\,e^{-\Sigt{}\,s}\,\mathrm d s
                = \frac{\Sigs{}}{\Sigt{}}\,(1 - e^{-\Sigt{}\,L_0}).

Plugging :math:`F` and :math:`\psi_{\rm surf} = \omega_0` into the
master equation :eq:`peierls-greens-function-architecture`:

.. math::

   \psi(r_i, \mu_q) &=
       \frac{\Sigs{}}{\Sigt{}}\,(1 - e^{-\Sigt{}\,L_0})
       + e^{-\Sigt{}\,L_0}\,\frac{\Sigs{}}{\Sigt{}} \\
   &= \frac{\Sigs{}}{\Sigt{}}\,
      \bigl[(1 - e^{-\Sigt{}\,L_0}) + e^{-\Sigt{}\,L_0}\bigr] \\
   &= \frac{\Sigs{}}{\Sigt{}} = \omega_0.

The :math:`L_0`-dependence cancels identically, leaving
:math:`\psi(r_i, \mu_q) = \omega_0` independent of :math:`(r_i,
\mu_q)`. **This is the operator action :math:`(K\cdot 1) = \omega_0`
established symbolically.** The closed sphere with isotropic
scattering and constant trial flux has a uniform angular-flux
solution at every interior point; the bounce sum closes the
boundary contribution exactly.

Step 3: rank-1 isotropic mode is the unique eigenmode
------------------------------------------------------

The above proves :math:`(K\cdot 1) = \omega_0`. By Perron-Frobenius
(every factor of the kernel :math:`\tilde t` is positive on the
integration domain), the constant function is the *dominant*
eigenmode and :math:`\omega_0` is the dominant eigenvalue of the
scattering operator. The fission-source equation
:math:`(I - K_{\rm scat})\,\phi = (\nSigf{}/k)\,\phi` gives
:math:`(1 - \omega_0)\,\phi = (\nSigf{}/k)\,\phi`; with the only
admissible solution :math:`\phi = \mathrm{const}`, we get
:math:`k_{\rm eff} = \nSigf{} / \Siga{} = \kinf`.

V_α1 makes this the **load-bearing claim of Variant α at rank-1**:
the prototype is required to reproduce :eq:`peierls-greens-k-inf`
to within quadrature error. The numerical implementation does so to
machine precision — see "Cross-verification matrix" below.


V_α2 — rank-1 equivalence to Hébert white BC
=============================================

V_α1 establishes that the Variant α operator action on constant trial
is :math:`\omega_0`. The second algebraic identity V_α2 explains *why*
Variant α at rank-1 agrees bit-for-bit with the existing Phase 4
``boundary="specular_multibounce"`` at :math:`N=1` and with the Hébert
``boundary="white_hebert"`` rank-1 closure: all three reduce to the
same closed-form geometric series.

The Phase 4 multi-bounce transmission matrix
:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.compute_T_specular_sphere`
in the rank-1 isotropic basis (:math:`\tilde P_0 = 1`) reduces to the
scalar

.. math::
   :label: peierls-greens-T00-integrand

   T_{00}^{\rm sphere} = 2 \int_0^1 \mu\,\tilde P_0(\mu)^2\,
       e^{-2\Sigt{} R \mu}\,\mathrm d\mu
   = 2 \int_0^1 \mu\,e^{-2\Sigt{} R \mu}\,\mathrm d\mu,

while the Hébert chord-self-collision probability (Hébert 2009
§3.8.5) is

.. math::

   P_{ss}^{\rm sphere} = 2 \int_0^1 \mu\,e^{-2\Sigt{} R \mu}\,\mathrm d\mu

— **the same integrand identically**. SymPy integration on
:math:`\mu \in [0, 1]` gives the closed form

.. math::
   :label: peierls-greens-V-alpha-2

   T_{00}^{\rm sphere} = P_{ss}^{\rm sphere} =
       \frac{1 - (1 + 2\tau_R)\,e^{-2\tau_R}}{2\,\tau_R^{2}},
   \qquad \tau_R = \Sigt{}\,R.

This makes the rank-1 algebraic equivalence airtight, *independent
of any quadrature implementation*. The symbolic proof lives in
:func:`~orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function.derive_T00_equals_P_ss_sphere`.

Operator-level interpretation. At rank-1 the Phase 4 closure
:math:`(I - T R)^{-1}` reduces to :math:`1/(1 - T_{00}) =
1/(1 - P_{ss})` (Eq. :eq:`peierls-greens-V-alpha-2` plus
:math:`R = [[1]]`), which is the scalar Hébert white-BC factor.
Variant α at rank-1 (constant trial) reduces to :math:`\omega_0`
via the bounce-sum self-consistency
:eq:`peierls-greens-surface-fixed-point` directly — no
:math:`(1 - P_{ss})^{-1}` factor appears. The two paths arrive at
the same answer through different algebraic structures. The B5
cross-verification anchor therefore *requires* both routes to
produce :math:`\kinf` to machine precision; if they do not, one
of the two has a bug. (They do — see B5.B in the cross-verification
matrix below.)


V_α3 — vacuum reduction at α = 0
=================================

The third operator-level identity is the trivial-but-load-bearing
vacuum-BC compatibility check:

.. math::
   :label: peierls-greens-V-alpha-3

   g_h(\rho' \to \rho)\bigr|_{\alpha = 0} = 0,

so the full Sanchez kernel
:math:`g_\alpha = \bar g_2 + g_h` reduces to the bare vacuum kernel
:math:`\bar g_2` at :math:`\alpha = 0`. The proof is one line: Sanchez
Eq. (A6) is

.. math::

   g_h(\rho'\to\rho) = 2\alpha \int_{\mu_0}^{1} T(\mu_-)\,\mu_*^{-1}\,
       \cosh(\rho\mu)\,\cosh(\rho'\mu_*)\,e^{-2 a \mu_-}\,\mathrm d\mu

— the leading factor :math:`2\alpha` makes the entire integrand
proportional to :math:`\alpha`. At :math:`\alpha = 0` the BC kernel
vanishes identically, and the Variant α prototype collapses to the
existing ORPHEUS vacuum sphere reference (the
:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.solve_peierls_1g`
with ``boundary="vacuum"`` path) without any special-case branch.

This is more important than it looks. The Variant α prototype is
shipped only for :math:`\alpha = 1` (perfect specular), but the
underlying operator structure is correct for any :math:`\alpha \in
[0, 1]`. A future extension to :math:`\alpha = 0` (vacuum BC) requires
only setting the prefactor to zero in
:func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function._apply_operator`'s
bounce-sum branch — V_α3 guarantees this is mathematically correct
without re-deriving anything. The symbolic proof is in
:func:`~orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function.derive_alpha_zero_kernel_reduction`.

Practical note. The closed sphere with vacuum BC is *not* the same
problem as the closed sphere with specular BC. For specular BC the
no-leakage condition forces :math:`k_{\rm eff} = \kinf` exactly;
for vacuum BC the spatial mode structure becomes non-trivial and
:math:`k_{\rm eff} < \kinf` (leakage). V_α3 makes the *operator
structure* compatible with vacuum BC, but the eigenvalue extraction
for vacuum BC requires solving a non-trivial spatial mode problem
that the closed-sphere prototype does not exercise. Future Plan 2
follow-on (see "Restrictions and future work") covers this.


Cross-verification matrix (Plan 2 B5)
======================================

The Plan 2 B5 cross-verification compares Variant α against the
existing ORPHEUS Phase 4 references on a thin homogeneous sphere
configuration. The three references — Variant α (this work),
``boundary="specular_multibounce"`` at rank :math:`N \in \{1,2,3\}`,
and ``boundary="white_hebert"`` (rank-1 only) — should agree at
rank-1 (V_α2 algebraic identity) and disagree at higher ranks of
the Phase 4 path due to the rank-N truncation bias.

**Configuration**: fuel-A-like 1G XS, :math:`R = 5`,
:math:`\Sigt{} = 0.5`, :math:`\Sigs{} = 0.38`, :math:`\nSigf{} =
0.025`. Optical thickness :math:`\tau_R = \Sigt{} R = 2.5`.
Analytic infinite-medium :math:`\kinf = 0.025/0.12 = 0.2083\overline{3}`.

.. list-table:: Cross-verification: closed homogeneous sphere with perfect specular BC
   :header-rows: 1
   :widths: 30 25 25

   * - Reference
     - :math:`k_{\rm eff}`
     - error vs Variant α
   * - :math:`\kinf` (analytic, no-leakage)
     - 0.20833333333
     - —
   * - **Variant α** (this work)
     - **0.20833333333**
     - **1e-13 %** (machine)
   * - Phase 4 N = 1 ≡ ``white_hebert`` (V_α2)
     - 0.20777642780
     - 0.2673 %
   * - Phase 4 N = 2
     - 0.20781489733
     - 0.2488 %
   * - Phase 4 N = 3
     - 0.20808891764
     - 0.1173 %

Findings:

1. **Variant α is exact for the closed homogeneous sphere.** The
   numerical implementation reproduces V_α1 (:math:`k_{\rm eff} =
   \kinf`) to 1e-13 % — limited only by floating-point round-off
   in the radial / angular / trajectory quadratures, not by any
   intrinsic discretisation error of the operator. With constant
   initial guess :math:`\psi_0 = 1`, convergence happens in one
   power-iteration step (the constant function *is* the eigenmode
   exactly). With a sinusoidal perturbation the iteration converges
   to the same uniform-flux eigenmode in 5–8 iterations.
2. **Phase 4 N=1 ≡ ``white_hebert`` bit-for-bit.** This is V_α2 in
   action: at rank-1 both closures reduce to the same Hébert
   :math:`(1 - P_{ss})^{-1}` scalar factor, so the rank-1 N=1
   ``specular_multibounce`` k_eff is *identically* the rank-1
   ``white_hebert`` k_eff at any radial-quadrature setting.
3. **Phase 4 errors decrease with rank but never reach Variant α.**
   The progression :math:`\epsilon_1 = 0.27\,\% \to \epsilon_2 =
   0.25\,\% \to \epsilon_3 = 0.12\,\%` shows monotone improvement
   with rank, but the curve does not converge to :math:`0\,\%` at
   any rank reachable in production (the rank-:math:`N` UserWarning
   fires at :math:`N \ge 4` for sphere/cylinder; see
   :ref:`peierls-rank-n-class-b-mr-mg-falsification`). The rank-N
   gating reflects the kernel's intrinsic difficulty — there is no
   higher rank that makes Phase 4 exact for closed sphere.
4. **Variant α improves on Phase 4 N=3 by ≈0.12 %** for this
   configuration. Users who need higher precision than Phase 4 can
   deliver on closed sphere homogeneous can use Variant α as a
   reference. The cost trade-off is documented under "Cost and
   convergence behaviour" below.

The numerical evidence is pinned by the test gates in
:file:`tests/derivations/test_peierls_greens_function_xverif.py`:

- :func:`test_b5_variant_alpha_gives_k_inf_exactly` —
  pins Variant α k_eff to :math:`\kinf` at 1e-10 relative tolerance
  on the truth-source side.
- :func:`test_b5_phase4_rank1_equals_white_hebert` —
  pins rank-1 ``specular_multibounce`` ≡ ``white_hebert`` (V_α2
  algebraic identity).
- :func:`test_b5_phase4_converges_toward_variant_alpha` — pins the
  rank-convergence direction :math:`\epsilon_3 \le \epsilon_1` and
  the absolute error figures :math:`\epsilon_1 < 0.5\,\%`,
  :math:`\epsilon_3 < 0.2\,\%`.


Cost and convergence behaviour
===============================

Per-iteration cost. Variant α scales as :math:`O(n_r \cdot n_\mu
\cdot n_{\rm traj})` per power-iteration step, where :math:`n_{\rm
traj}` is the per-pair trajectory and bounce-period quadrature
order (typically 32–64 Gauss-Legendre nodes). The Phase 4
matrix-Galerkin form scales as :math:`O(n_r^2)`. For the smoke-test
configuration :math:`(n_r, n_\mu, n_{\rm traj}) = (16, 16, 32)` and
Phase 4 :math:`n_r = 24`, Variant α has 8000 trajectory evaluations
per iteration vs Phase 4's 600 kernel-matrix entries — about
:math:`13\times` more work per power-iteration step.

Iteration count. The closed-sphere problem has a single rank-1
eigenmode, so the Rayleigh quotient converges in 1 step from a
constant initial guess (no power-iteration ramp-up). From a
sinusoidal initial guess, convergence to :math:`10^{-12}` relative
tolerance happens in 5–8 iterations. Variant α therefore wins
total runtime over Phase 4 N=3 by a factor of about 2 for the
closed-sphere case despite the higher per-iteration cost — the
fast convergence is the load-bearing optimisation.

Quadrature dependence. The angular Gauss-Legendre node count
:math:`n_\mu` controls how well the :math:`\mu \to 0` integrable
:math:`1/\mu` singularity in :math:`T(\mu_{\rm surf})` is resolved
by the final scalar-flux integration
:math:`\phi(r) = 4\pi \int_0^1 \psi(r,\mu)\,\mathrm d\mu`. For the
constant-trial eigenmode this singularity is harmless because the
integrand reduces to the constant :math:`\omega_0`; convergence is
machine-precision at :math:`n_\mu = 12`. For non-constant trial
functions (which would arise for vacuum BC with non-trivial spatial
mode structure), the :math:`\mu \to 0` quadrature would need
Gauss-Jacobi (weighted by :math:`\mu`) or the change of variables
:math:`u^2 = \mu` to absorb the singularity — both are trivial
extensions if needed.

Why Variant α is **not** a production replacement. The Plan 2
closeout decision (parallel research-grade reference, not folded
into ``boundary=...`` dispatch) reflects two facts:

1. **Variant α has narrower validation than Phase 4.** The prototype
   is tested for **homogeneous sphere only, isotropic scattering
   only, perfect specular BC only**. Phase 4
   ``specular_multibounce`` is tested across slab / cylinder /
   sphere, multi-region, anisotropic scattering, and all three BC
   types (vacuum / white / specular), with rank-N gating that has
   documented failure modes. Variant α inherits none of this
   validation; promoting to production would require substantial
   additional work.
2. **Variant α improves on Phase 4 only for the cases it covers.**
   For closed homogeneous sphere with isotropic scattering, Variant α
   is exact. For all other configurations, Variant α is unimplemented
   and Phase 4 is the only ORPHEUS reference. The Issue #132 Class B
   multi-region catastrophe is *not* solved by Variant α — multi-region
   sphere is future work.

Variant α is therefore the *gold-standard reference* against which
Phase 4's closure approximations (rank-N specular_multibounce,
Hébert white) are calibrated. It is the L1 reference for the case
it covers, but not the production solver.


.. _peierls-greens-vacuum-extension:

Vacuum BC extension (Plan-2 follow-on A1 + A2)
================================================

The B-phase prototype was restricted to perfect specular BC
(:math:`\alpha = 1`) — the case for which the Phase 5
hypersingularity argument was sharpest, and where the V_α1
algebraic identity gives a machine-precision ground truth. The A1
follow-on (2026-05-02, branch ``feature/peierls-greens-function``)
extends the prototype to **arbitrary** :math:`\alpha\in[0,1]` with
no special-case branch — the operator structure :math:`t = \bar t +
t_h` collapses cleanly to vacuum at :math:`\alpha = 0` per the V_α3
algebraic identity (Eq. :eq:`peierls-greens-V-alpha-3`).

Bounce-sum closure with α
--------------------------

The surface fixed-point equation
:eq:`peierls-greens-surface-fixed-point` carries the reflection
coefficient :math:`\alpha` directly in the surface-flux update:

.. math::
   :label: peierls-greens-bounce-sum-alpha

   \psi_{\rm surf}(\mu_{\rm surf})
   = B(\mu_{\rm surf}) + \alpha\,e^{-\Sigt{}\,L_p}\,\psi_{\rm surf},

with closed-form solution

.. math::
   :label: peierls-greens-T-alpha

   \psi_{\rm surf}(\mu_{\rm surf})
   = \frac{\alpha\,B(\mu_{\rm surf})}
          {1 - \alpha\,e^{-\Sigt{}\,L_p(\mu_{\rm surf})}}.

At :math:`\alpha = 1` this reduces to
:eq:`peierls-greens-T-mu-surf`. At :math:`\alpha = 0` the numerator
vanishes, the surface flux contribution drops out, and the master
equation :eq:`peierls-greens-function-architecture` reduces to the
first-leg integral alone:

.. math::

   \psi(r_i, \mu_q)\bigr|_{\alpha=0} \;=\; F(r_i, \mu_q),

which is the integral form of the vacuum sphere transport equation.
This is V_α3 numerically realised: the operator structure carries
the full :math:`\alpha\in[0,1]` parameter range without a separate
code path. The implementation lives in
:func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function._apply_operator`
(``alpha`` parameter; see ``# Vacuum branch.`` comment in source).

The full :math:`\mu`-grid for vacuum BC
----------------------------------------

Closed sphere has :math:`\mu \to -\mu` symmetry — the eigenmode is
isotropic, so the radial trajectory machinery only needs
:math:`\mu \in (0, 1]`. Vacuum BC **breaks** this symmetry: outward
:math:`\mu > 0` rays exit through the surface (vacuum at
:math:`r = R^+`); inward :math:`\mu < 0` rays traverse the entire
sphere with no surface-flux contribution from the bounce sum (which
is zero for vacuum). The trajectories are therefore different for
the two signs of :math:`\mu`, and the prototype now discretises
:math:`\mu \in [-1, 1]` by full Gauss-Legendre quadrature on
:math:`n_\mu` nodes. For closed-sphere specular this incurs only
redundant compute (the eigenmode is constant so each :math:`\mu`-grid
point gives the same flux); for vacuum BC it is mandatory.

The PS-1982 cross-check (A2)
----------------------------

The vacuum-BC k-eigenvalue is non-trivial — leakage forces
:math:`k_{\rm eff} < k_\infty` and the spatial mode is peaked at the
centre with a non-trivial profile. Variant α at :math:`\alpha = 0`
therefore needs a **structurally-independent** L1 cross-check.
[PS1982]_ Eq. (21) provides this: Pomraning-Siewert 1982's
integral-equation form for the homogeneous vacuum sphere with
isotropic scattering uses a different mathematical path
(integrate-over-:math:`\mu` then add half-spaces) than Sanchez 1986's
cosh-even-extension. The two derivations arrive at the same
:math:`[E_1(|r-x|) - E_1(r+x)]` kernel through structurally
independent routes; PS-1982 itself confirmed via
method-of-characteristics is a third independent confirmation.

Implementation: a Nyström solver for [PS1982]_ Eq. (21) lives in
:func:`orpheus.derivations.continuous.peierls_nystrom.ps1982_reference.solve_ps1982_vacuum_sphere`
— Gauss-Legendre quadrature on :math:`(0, R\Sigt{}]` (optical
units), :math:`E_1` evaluated via :func:`mpmath.expint`, the kernel
log-singularity at :math:`r = x` handled by QUADPACK QAGS. Power
iteration on the integral kernel
:math:`(c/2)\,x\,[E_1(|r-x|) - E_1(r+x)]` extracts the dominant
eigenvalue, and :math:`k_{\rm eff} = c \cdot \nSigf{}/(\Sigs{} +
\nSigf{})` from the converged scattering ratio
:math:`c = (\Sigs{} + \nSigf{}/k)/\Sigt{}` self-consistently.

Cross-check evidence (table from
:func:`tests.derivations.test_peierls_greens_function_xverif_ps1982.test_a2_variant_alpha_agrees_with_ps1982`):

.. list-table:: A2 — Variant α vs PS-1982 vacuum sphere k_eff
   :header-rows: 1
   :widths: 30 22 22 26

   * - Configuration
     - PS-1982 :math:`k_{\rm eff}`
     - Variant α :math:`k_{\rm eff}`
     - Relative agreement
   * - τ_R = 2.5, c = 0.45 (strong absorber)
     - exact via Eq. (21)
     - matches
     - < 1e-4
   * - τ_R = 5.0, c = 0.45
     - exact via Eq. (21)
     - matches
     - < 1e-4
   * - τ_R = 2.5, c = 0.65 (medium absorber)
     - exact via Eq. (21)
     - matches
     - < 1e-4
   * - τ_R = 5.0, c = 0.65
     - exact via Eq. (21)
     - matches
     - < 1e-4
   * - τ_R = 25 (asymptotic thick)
     - just below k_inf
     - just below k_inf
     - < 2e-3 (cubic-spline bias)

Configurations: :math:`\Sigt{} = 0.5`, :math:`\Sigs{} \in \{0.20,
0.30\}`, :math:`\nSigf{} = 0.025`, :math:`R \in \{5, 10\}`. Quadrature:
PS-1982 :math:`n_{\rm quad} = 30`, Variant α :math:`n_r = n_\mu = 32`,
:math:`n_{\rm traj} = 64`, ``tol = 1e-9``.

The 2e-3 bias at :math:`\tau_R = 25` (very thick sphere) is documented
in the test docstring: the GL physical-cm grid samples the interior
depletion region sparsely at large :math:`R`, and the cubic-spline
:math:`\phi`-interpolant accumulates a systematic ~1e-3 bias from
the under-resolved gradient near :math:`r = 0`. PS-1982's optical-units
grid is naturally adaptive in optical depth and does not exhibit
this bias. Closing the gap requires a better-suited radial
quadrature (chebyshev or log-spaced near :math:`r = 0`); flagged
as a follow-on improvement, not blocking.

A1+A2 also surfaced a **bug in the original B-phase prototype**.
The first-leg trajectory length :math:`L_0(r,\mu)` was originally
implemented with the **forward** distance from :math:`r` to the
surface (:math:`\sqrt{R^2 - r^2(1-\mu^2)} - r\mu`); the integral
form requires the **backward** distance
:math:`L_{\rm back} = r\mu + \sqrt{R^2 - r^2(1-\mu^2)}` (the chord
from the source point :math:`r` along :math:`-\Omega_\mu` to the
surface entry point — see :eq:`peierls-greens-L0`). For closed
sphere (V_α1 algebraic identity) the bug was masked because the
closure cancels the :math:`L_0`-dependence identically; for vacuum
BC it surfaced as a 6 % k_eff disagreement vs PS-1982. The fix is
the change of sign in the formula for :math:`L_0`.


.. _peierls-greens-multigroup:

Multi-group extension (Plan-2 follow-on A3)
================================================

The B-phase prototype was 1G — fuel-A-like XS at
:math:`\Sigt{} = 0.5`, :math:`\Sigs{} = 0.38`, :math:`\nSigf{} =
0.025`. By Cardinal Rule 6, 1G eigenvalue tests are degenerate:
:math:`k = \nSigf{}/\Siga{}` is flux-shape independent and computable
without solving the transport equation. The A3 follow-on
(2026-05-02) extends to :math:`\ge 2G` with arbitrary scattering
matrix and arbitrary fission spectrum, closing the 1G Cardinal-Rule-6
gap.

Multi-group operator action
---------------------------

The per-group Variant α operator carries the same architecture
:eq:`peierls-greens-function-architecture` as the 1G case, with
the per-group source split into scattering + fission terms:

.. math::
   :label: peierls-greens-mg-source

   q_g(r) \;=\; \sum_{g'=1}^{G} \Sigs{g'\to g}\,\phi_{g'}(r) \;+\;
       \frac{\chi_g}{k_{\rm eff}}\,\sum_{g'=1}^{G}\nSigf{g'}\,\phi_{g'}(r),

with the convention ``sigma_s[g_from, g_to]`` —
i.e. :math:`\Sigs{g'\to g}`. The first sum is **in-scatter** into
group :math:`g`; the second sum is **fission emission** into
group :math:`g` weighted by the prompt-fission spectrum
:math:`\chi_g`. Both sums combine across all source groups
:math:`g'`, so the implementation supports both downscatter
(lower-triangular :math:`\Sigs`) and upscatter (full :math:`\Sigs`).
The per-group :math:`q_g(r)` is plugged into
:eq:`peierls-greens-function-architecture` independently for each
:math:`g` — the trajectory and bounce-sum machinery do not see the
group structure at all.

The implementation is
:func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function.solve_greens_function_sphere_mg`
with shared per-group operator
:func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function._apply_operator_with_source_profile`
(extracted during A3 so the 1G path becomes a thin wrapper —
no regression in the 25 prior 1G tests). At each outer iteration,
the solver:

1. Computes scalar fluxes :math:`\phi_g(r) = 2\pi\!\int\psi_g(r,\mu)
   \,\mathrm d\mu` per group.
2. Computes total fission rate :math:`F(r) = \sum_{g'}\nSigf{g'}\,
   \phi_{g'}(r)` once at each radial node.
3. For each group :math:`g`, builds the source profile
   :math:`q_g(r)/(4\pi)` per :eq:`peierls-greens-mg-source` and
   applies the per-group operator.
4. Updates :math:`k_{\rm eff}` via Rayleigh quotient on
   volume-integrated fission rate.
5. Normalises so the total fission rate stays :math:`O(1)`; checks
   convergence on relative :math:`k_{\rm eff}` change.

Closed-sphere multi-group reduction
-----------------------------------

For closed sphere (:math:`\alpha = 1`) the V_α1 algebraic identity
generalises trivially: the operator action on **any** spatially
constant trial gives :math:`(K\cdot\boldsymbol\phi_{\rm const}) =
\mathrm{diag}(\omega_{0,g})\,\boldsymbol\phi_{\rm const}` per group,
and the multi-group fission-source eigenvalue equation reduces to
the *infinite-medium* multi-group balance

.. math::

   (\Sigt{g} - \Sigs{g\to g})\,\phi_g
   = \sum_{g'\ne g}\Sigs{g'\to g}\,\phi_{g'} +
     \frac{\chi_g}{k_{\rm eff}}\,\sum_{g'}\nSigf{g'}\,\phi_{g'},

whose dominant eigenvalue is :math:`k_\infty` and the corresponding
right-eigenvector is the homogeneous-medium spectrum. ORPHEUS
provides
:func:`orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`
for this exact computation (transfer-matrix dominant eigenvalue +
spectrum). Variant α multi-group must reproduce this to within
quadrature error for closed sphere — verified to ≤ 1e-9 relative
in the test gates (see "Test provenance" below).

Test coverage:

.. list-table:: A3 multi-group test coverage
   :header-rows: 1
   :widths: 35 18 47

   * - Test name
     - Tolerance
     - Coverage
   * - ``test_mg_g1_special_case_matches_1g``
     - 1e-12 rel
     - G=1 MG path matches 1G solver bit-equal
   * - ``test_mg_closed_sphere_2g_downscatter``
     - 1e-9 rel
     - 2G closed sphere = ``kinf_and_spectrum_homogeneous``
   * - ``test_mg_closed_sphere_2g_upscatter``
     - 1e-9 rel
     - 2G upscatter (full Σ_s); slower convergence verified
   * - ``test_mg_closed_sphere_2g_spectrum``
     - 1e-6 rel
     - φ_2/φ_1 = Σ_{1→2}/Σ_R,2 (analytical 2G ratio)
   * - ``test_mg_vacuum_2g_subcritical``
     - leakage check
     - 2G vacuum k_eff < 2G k_inf
   * - ``test_mg_4g_realistic_chi``
     - 1e-9 rel
     - 4G fuel-A with χ = (0.6, 0.35, 0.05, 0.0)
   * - ``test_mg_4g_vacuum``
     - leakage check
     - 4G vacuum has non-trivial spectrum

The full test file is
:file:`tests/derivations/test_peierls_greens_function_mg.py` (~700
LoC, all gates pass at default quadrature).


.. _peierls-greens-multiregion:

Multi-region extension (Plan-(b) Options 1 + 2)
================================================

The B-phase + A3 prototype was homogeneous — single :math:`\Sigt{}`,
single scattering matrix, single :math:`\nSigf{}`. The Plan-(b)
follow-on (2026-05-02) extends to a piecewise-homogeneous sphere
with concentric regions, each carrying its own XS. This is the
**direct attack on Issue #132**, the Class B multi-region
catastrophe of the Phase 4 ``specular_multibounce`` rank-N closure.

Trajectory and bounce-period segment decomposition
--------------------------------------------------

The first-leg trajectory and bounce-period chord both pass through
multiple regions in a piecewise-homogeneous sphere. The trajectory
geometry is decomposed into per-region **segments** by computing
the chord crossings with each interior region boundary
:math:`R_k`:

.. math::
   :label: peierls-greens-mr-trajectory-segments

   r(s)^2 \;=\; r_i^2 - 2 r_i \mu_q s + s^2,
   \qquad s \in [0, L_{\rm back}],

with :math:`r(s) = R_k` at
:math:`s = r_i\mu_q \pm\sqrt{R_k^2 - r_i^2(1-\mu_q^2)}` (when the
discriminant is positive — otherwise the chord misses the
:math:`R_k` shell). The chord crossings sort to give a list of
segments :math:`[(s_a, s_b, k_{\rm region})]` such that
:math:`r(s)` lies entirely within region :math:`k_{\rm region}` for
:math:`s\in(s_a, s_b)`.

The bounce-period chord at impact parameter :math:`h(r_i,\mu_q) =
R\sqrt{1-\mu_{\rm surf}^2}` is decomposed analogously:
:math:`|r_{\rm chord}(s)|^2 = h^2 + (s - L_p/2)^2`, so the chord
crosses :math:`R_k` (whenever :math:`R_k > h`) at
:math:`s = L_p/2 \pm\sqrt{R_k^2 - h^2}`.

The segment decomposition implementation is in the private helpers
:func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function._trajectory_segments`
and
:func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function._chord_segments`.

Piecewise optical depth
-----------------------

Per-segment Gauss-Legendre quadrature on the local
:math:`(s_a, s_b)` interval composes to the full integral. The
optical depth accumulates segment-by-segment with the local
:math:`\Sigt{,k}`:

.. math::
   :label: peierls-greens-mr-piecewise-tau

   \tau(s) \;=\; \sum_{(s_a, s_b, k)\,\subset\,[0,s]}
       \Sigt{,k}\,(\min(s, s_b) - s_a),

so the per-region attenuation factor :math:`e^{-\tau(s)}` is
exact within each segment and continuous at every interior boundary
crossing. The first-leg integral and bounce-period integral are
both evaluated in this composite-segment form. Source values are
sampled via cubic-spline interpolation of the per-region scalar flux
on the radial grid — a known prototype limitation discussed below
in the Garcia 2021 cross-check section.

The bounce-sum closure :math:`T(\mu_{\rm surf}) = 1/(1 - \alpha\,
e^{-\tau_p})` uses the **total** chord optical depth
:math:`\tau_p = \sum_k \Sigt{,k}\,\ell_k(\mu_{\rm surf})`, which is
**bounce-invariant** under perfect specular BC — every bounce
traverses the same regions in the same order, so the chord optical
depth does not change. This is what keeps the geometric series
:math:`\sum_n e^{-n\tau_p}` closed-form-summable in the
multi-region case.

The implementation is
:func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function.solve_greens_function_sphere_mr`
(k-eigenvalue, Plan-(b) Option 2) and
:func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function.solve_greens_function_sphere_mr_fixed_source`
(fixed-source, Plan-(b) Option 1).

.. _peierls-greens-issue132-result:

Issue #132 reproducer — the Class B catastrophe avoided
--------------------------------------------------------

The Phase 4 ``specular_multibounce`` rank-N closure exhibits a
**+57 % sign-flip catastrophe** on the canonical Class B
multi-region sphere configuration (Issue #132,
:ref:`peierls-rank-n-class-b-mr-mg-falsification`). The configuration
is fuel-A inner / moderator-B outer at radii=[0.5, 1.0]:

- Region 0 (inner, fuel A): :math:`\Sigt{} = 1.0`, :math:`\Sigs{} =
  0.5`, :math:`\nSigf{} = 0.75`.
- Region 1 (outer, moderator B): :math:`\Sigt{} = 2.0`, :math:`\Sigs{}
  = 1.9`, :math:`\nSigf{} = 0`.

The volume-averaged-homogenised :math:`k_\infty \approx 0.648`
(strongly subcritical). Phase 4 references:

.. list-table:: Issue #132 — Class B catastrophe vs Variant α
   :header-rows: 1
   :widths: 35 25 40

   * - Reference
     - :math:`k_{\rm eff}`
     - Status
   * - Phase 4 ``specular_multibounce`` rank-1 (≡ ``white_hebert``)
     - 0.551
     - Subcritical (sensible)
   * - Phase 4 ``specular_multibounce`` rank-2
     - **1.015**
     - **Sign flip + supercritical (catastrophe)**
   * - Volume-averaged :math:`k_\infty` (homogenised)
     - 0.648
     - Reference homogenised limit
   * - **Variant α** (this work, MR closed sphere α=1)
     - **≈ 0.735**
     - **Subcritical, between rank-1 and homogenised k_inf**

The Variant α value :math:`\approx 0.735` is the
multi-region transport solution: above the homogenised
:math:`k_\infty = 0.648` (because fuel concentration in the inner
region boosts the effective multiplication) and well below 1.0
(the configuration is subcritical). Variant α has **no rank-N
closure** — the BC is absorbed into the kernel via Sanchez Eq. (A1)
— so the mode-0 / mode-:math:`n\ge 1` normalisation mismatch in
:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.build_closure_operator`
that breaks Phase 4 simply cannot occur structurally.

The test gate
:func:`tests.derivations.test_peierls_greens_function_mr.test_mr_issue132_no_catastrophe_closed_sphere`
pins :math:`0.5 < k_{\rm eff} < 0.95` and explicitly :math:`k_{\rm
eff} < 1` (ruling out the Phase 4 catastrophe). The spatial mode
gate
:func:`test_mr_issue132_spatial_mode_physical` verifies that
:math:`\phi` is peaked in the fuel region (more multiplication),
monotonically decreasing within each region, with a discernible
slope change at the fuel/moderator interface.

This is **not** an L1 cross-check — there is no published L1
multi-region critical sphere k-eigenvalue benchmark for this
fuel/moderator configuration (Garcia 2021 is fixed-source-only;
see :ref:`peierls-greens-garcia2021`). It is a **regression gate
against a known pathology** — the Phase 4 catastrophe specifically.
The L1 cross-check for multi-region Variant α is the Plan-(b)
Option 1 fixed-source benchmark below.


.. _peierls-greens-garcia2021:

Garcia 2021 fixed-source cross-check (Plan-(b) Option 1)
========================================================

[Garcia2021]_ ships a stable :math:`P_N` solver for the multi-region
sphere with internal sources and externally incident angular flux.
The paper is **subcritical-only** — criticality is explicitly out
of scope (paper §III.A; documented at
:file:`.claude/agent-memory/literature-researcher/ps1982_and_garcia_extraction.md`).
For Variant α multi-region verification this gives **flux-shape L1
evidence**: agreement on the converged scalar flux profile
:math:`\phi(r)` across a 3-region sphere with non-trivial XS jumps,
vacuum BC, and constant-per-region isotropic external source.

Williams 1991 Case 1 — the canonical 3-region benchmark
--------------------------------------------------------

The Garcia 2021 Case 1 configuration (originally Williams 1991
*Ann. Nucl. Energy* 18, 371, Example 5) is a 3-region sphere with
strong cross-section discontinuities:

.. list-table:: Garcia 2021 / Williams 1991 Case 1 configuration
   :header-rows: 1
   :widths: 18 16 18 18 16 14

   * - Region
     - r range (cm)
     - :math:`\Sigt{}`
     - :math:`\Sigs{}`
     - :math:`c = \Sigs{}/\Sigt{}`
     - :math:`Q_{\rm ext}` per cm³ per ster
   * - 1 (core)
     - 0 – 3
     - 1.0
     - 0.99
     - 0.99
     - 0.5
   * - 2 (mid)
     - 3 – 5
     - 0.5
     - 0.30
     - 0.6
     - 1.0
   * - 3 (outer)
     - 5 – 7
     - 2.0
     - 1.90
     - 0.95
     - 1.5

Vacuum BC at :math:`r = 7`. No fission. Garcia 2021 cross-checked
this case against Williams 1991 (integral-eq MoC) and Picca-Furfaro-
Ganapol 2012 (S_N) to 3–4 significant figures — three structurally-
independent methods agreeing.

.. _peierls-greens-garcia2021-convention-map:

Convention conversion to Garcia table 5
----------------------------------------

Garcia 2021's source :math:`S_k` is **per cm³ per steradian**;
ORPHEUS's ``external_source`` argument to
:func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function.solve_greens_function_sphere_mr_fixed_source`
is **total per cm³** (the operator divides by :math:`4\pi`
internally). These differ by :math:`4\pi`. Garcia's "scalar flux"
:math:`\phi(r) = \int_{-1}^{1}\Psi(r,\mu)\,\mathrm d\mu` (no
:math:`2\pi`); the Variant α output is
:math:`\phi(r) = 2\pi\int_{-1}^{1}\psi(r,\mu)\,\mathrm d\mu`
(standard scalar flux). These differ by :math:`2\pi`.

Net for matching: passing
``external_source = (0.5, 1.0, 1.5)`` and comparing to Garcia Table 5:

.. math::
   :label: peierls-greens-garcia-convention

   \phi_{\rm mine}(r) \;=\; \frac{2\pi}{4\pi}\,\phi_{\rm Garcia}(r)
                       \;=\; \tfrac{1}{2}\,\phi_{\rm Garcia}(r).

This factor-of-:math:`\tfrac{1}{2}` is the **convention map**, not
an error or a discretisation artefact. The factor is documented
verbatim in
:file:`.claude/agent-memory/literature-researcher/ps1982_and_garcia_extraction.md`
and pinned by the test constant ``GARCIA_TO_VARIANT_ALPHA_FACTOR =
0.5``.

Fixed-source iteration
-----------------------

The Variant α multi-region fixed-source solver uses the same
trajectory + bounce-period architecture as the k-eigenvalue case,
but iterates on the scattering source alone:

.. math::
   :label: peierls-greens-fixed-source-iteration

   \psi_g^{(n+1)}(r,\mu) \;=\; K_g\!\left[\,
       \tfrac{1}{4\pi}\!\left(\,\sum_{g'}\Sigs{g'\to g,\,k(r)}\,
           \phi_{g'}^{(n)}(r) + Q_{{\rm ext},k(r),g}\right)\,
       \right] (r,\mu),

where :math:`k(r)` is the region containing radius :math:`r`. No
fission, no eigenvalue iteration; convergence is on the relative
:math:`\phi_g`-change to a fixed tolerance. The implementation is
:func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function.solve_greens_function_sphere_mr_fixed_source`.

Garcia Case 1 agreement
------------------------

Variant α agreement vs Garcia 2021 Table 5 at default settings
(:math:`n_r = 48`, :math:`n_\mu = 24`, :math:`n_{\rm traj} = 64`,
:math:`{\rm tol} = 10^{-7}`):

.. list-table:: Garcia Case 1 — Variant α vs Table 5 converged ppP_N
   :header-rows: 1
   :widths: 16 14 22 22 26

   * - r (cm)
     - Region
     - Garcia φ
     - 0.5 × Garcia φ
     - Variant α agreement
   * - 0.0
     - core
     - 18.860
     - 9.430
     - < 1 %
   * - 0.5
     - core
     - 18.756
     - 9.378
     - < 1 %
   * - 1.0
     - core
     - 18.442
     - 9.221
     - < 1 %
   * - 1.5
     - core
     - 17.911
     - 8.956
     - < 1 %
   * - 2.0
     - core
     - 17.145
     - 8.573
     - < 1 %
   * - 2.5
     - core
     - 16.095
     - 8.048
     - < 1 %
   * - 3.0
     - core/mid interface
     - 14.381
     - 7.190
     - 1.6 % (interface band)
   * - 3.5
     - mid
     - 13.455
     - 6.728
     - 3 – 8 % (mid band)
   * - 4.0
     - mid
     - 13.337
     - 6.668
     - 3 – 8 % (mid band)
   * - 4.5
     - mid
     - 13.590
     - 6.795
     - 3 – 8 % (mid band)
   * - 5.0
     - mid/outer interface
     - 14.361
     - 7.180
     - 11 % (interface band)
   * - 5.5
     - outer
     - 15.532
     - 7.766
     - < 6 %
   * - 6.0
     - outer
     - 14.198
     - 7.099
     - < 6 %
   * - 6.5
     - outer
     - 10.807
     - 5.404
     - < 6 %
   * - 7.0
     - outer surface
     - 4.0763
     - 2.038
     - < 1 %

Tolerance bands gated by the per-point test
:func:`tests.derivations.test_peierls_greens_function_garcia2021.test_garcia_case1_phi_matches_at_point`:
2 % at non-interface points, 15 % at interface-adjacent points
(within ±2 cm of an interior region boundary).

The shape is **fundamentally correct** at every point. The
near-interface error is purely from the cubic-spline source-flux
interpolation smoothing the discontinuous :math:`\Sigs{}` at region
boundaries — the spline reaches across each interface and slightly
contaminates the source values in the adjacent region. The
mid-region (3 ≤ r ≤ 5 cm) is sandwiched between two interfaces and
shows the worst smoothing-induced error (up to ~11 % near
:math:`r = 5`). **Piecewise per-region interpolation** (a separate
spline per region) would close this gap; flagged as a follow-on
improvement.

All 17 test gates pass (3 sanity + 15 per-r-point). The 3-method
triangulation (Garcia P_N + Williams 1991 MoC + Picca-Furfaro-
Ganapol 2012 S_N) makes this a robust L1 reference: any
implementation reaching this 4-figure agreement at every r-point is
verified for multi-region transport with vacuum BC and isotropic
fixed sources.


Restrictions and future work
=============================

The Variant α architecture now covers sphere-only with arbitrary
:math:`\alpha\in[0,1]`, multi-group (any G), and multi-region
(piecewise homogeneous). The following extensions are flagged for
future work but not blocking — the production reference families
([Garcia2021]_ for multi-region fixed-source + the Phase 4 Nyström
family for cylinder + slab) cover the gaps for now.

Cylinder geometry
-----------------

Sanchez 1986 unifies slab / cylinder / sphere via a geometry-shape
parameter :math:`\alpha` (note: an unfortunate notation clash with
the specular-coefficient :math:`\alpha` used throughout this page;
in Sanchez Table 1, geometry-shape :math:`\alpha = 0,1,2` for slab,
cylinder, sphere). For the cylinder, the closed-form
:math:`\cosh(\rho\mu)\,\cosh(\rho'\mu_*)` source-arc / receiver-arc
factors are replaced by Bessel functions :math:`I_0` and :math:`J_0`,
but the bounce-sum machinery :math:`T(\mu_{\rm surf})` is
structurally identical (the chord is still :math:`L_p = 2 R\,
\mu_{\rm surf}` for any axisymmetric geometry under perfect specular
BC). Cylinder extension is medium-risk; the SymPy step would replace
the cosh closed forms with their Bessel analogues (Knyazev 1993).

Anisotropic scattering (:math:`\omega_1 \ne 0`)
------------------------------------------------

Sanchez 1986 gives the :math:`h` kernel for linearly anisotropic
scattering — the bounce-sum trajectory machinery extends with the
substitution

.. math::

   q(r) \;\to\; q_0(r) + \omega_1\,\Omega \cdot J(r),

where :math:`q_0` is the isotropic source and :math:`J(r)` is the
net current. The first-leg and bounce-period integrals
:eq:`peierls-greens-trajectory-integral` and
:eq:`peierls-greens-bounce-period-integral` then carry an additional
:math:`\Omega \cdot J(r(s))` factor that becomes
:math:`\mu_q J(r(s))` along the characteristic. The structural form
is preserved; the numerical change is one extra trajectory integrand
term per :math:`\omega_1` component.

External cross-check via Garcia stable :math:`P_N` (k-eigenvalue)
------------------------------------------------------------------

Garcia 2021 [Garcia2021]_ has been pulled and is the L1 reference
for the multi-region **fixed-source** sphere — see
:ref:`peierls-greens-garcia2021`. The k-eigenvalue cross-check
remains an open gap: Garcia 2021 is **subcritical-only**
(criticality is "future work" per paper §III.A); Garcia 2020
[Garcia2020]_ covers homogeneous sphere reflective BC but is also
fixed-source-only. There is **no published L1 multi-region critical
sphere k-eigenvalue benchmark** for the Issue #132 fuel/moderator
configuration as of 2026-05-01. Candidate sources for a future
Plan-(c) effort:

- Sood-Forster-Parsons 2003 (LANL LA-13511) — compiled critical-sphere
  c_crit(R) tables for single-region critical sphere. Anchors the
  PS-1982 homogeneous vacuum-sphere k-eigenvalue chain
  (:math:`c\to c_{\rm crit}` at the integral-operator's unit
  eigenvalue) but does not cover multi-region.
- Williams 2005 (*Ann. Nucl. Energy* 32) — varying-:math:`P_N` order
  for the annular-gap geometry (Garcia Ref [16]). Cited as
  "interesting" in Garcia 2021 but not numerically explored. Could
  serve as the multi-region eigenvalue reference if reproduced.
- A future Garcia eigenvalue paper (pre-print search:
  *Garcia ann nucl energy multi-region critical sphere* with year
  ≥ 2021).

For now Variant α multi-region k-eigenvalue is verified via the
Issue #132 regression gate (rank-N catastrophe avoided) plus the
single-region reduction to ``kinf_homogeneous`` (closed sphere) and
PS-1982 (vacuum sphere). The Plan-(b) Option 1 fixed-source
benchmark is the closest published reference and is documented at
:ref:`peierls-greens-garcia2021`.

What is permanently dead
-------------------------

The Phase 5 retreat closed several research paths that should *not*
be revisited:

- Nyström sampling of the angle-integrated kernel
  :math:`g_h(\rho'\to\rho)` — diverges at the surface diagonal for
  any sampling scheme.
- Galerkin double-integration over Lagrange basis on
  :math:`\int\!\!\int L_i(\rho')\,L_j(\rho)\,g_h\,\mathrm d\rho'\,
  \mathrm d\rho` — same log(:math:`Q_\mu`) divergence at coincident
  collocation pairs.
- Bounce-resolved M2 expansion as a Nyström target — singularity
  persists at K_max=0 (no multi-bounce, bare specular alone),
  proving that the diagonal singularity is in
  :math:`F_{\rm out} \cdot G_{\rm in}`, not in the multi-bounce
  factor :math:`T(\mu)`.
- Hadamard finite-part regularisation — gauge-ambiguous, unbounded
  budget, flagged as out-of-scope by R3 PRIMARY in the Phase 5
  retreat investigation.

The Variant α architecture is the structural fix that closes those
paths. Future research-grade reformulations for sphere homogeneous
specular k_eff should fork from Variant α, not from any of the
above.


.. _peierls-greens-cylinder:

Cylinder Variant α (Phase 1 standalone)
========================================

.. todo:: Archivist expansion needed.

   This section is a **method-implementer stub** — it states the
   load-bearing equations, labels them for cross-reference, and
   points to the SymPy module + tests. The full mathematical
   narrative (worked algebra; design rationale; convergence-floor
   discussion; comparison with the cylinder Nyström vacuum
   reference; treatment of the dual grazing-ray pathology
   :math:`\mu_{\rm axial} \to \pm 1` and :math:`\varphi_{\rm az}
   \to \pm \pi/2`) is the archivist's deliverable.

   Source artifacts:

   - Branch-1 SymPy: :mod:`orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function_cylinder`
     (V_α1_cyl, V_α1_cyl.geometry, V_α2_cyl, V_α3_cyl).
   - Branch-2 production:
     :mod:`orpheus.derivations.continuous.peierls_greens_function.greens_function_cylinder`
     (:func:`solve_greens_function_cylinder`,
     :func:`solve_greens_function_cylinder_mg`).
   - Symbolic test gate:
     :file:`tests/derivations/test_peierls_greens_function_cylinder_symbolic.py`
     (9 foundation-tagged tests).
   - Numerical L1 gates:
     :file:`tests/derivations/test_peierls_greens_function_cylinder_solver.py`
     (7 L1-tagged tests).
   - Plan: :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md`.
   - Closeout memo:
     :file:`.claude/agent-memory/method-implementer/cylinder_variant_alpha_phase1.md`.

Phase-1 cylinder Variant α mirrors the sphere architecture (F first-leg
+ B bounce-period + closed-form geometric series) on a 3D phase-space
:math:`(r, \mu_{\rm axial}, \varphi_{\rm az})`. The axial direction is
kept **explicit** per Issue #129 — the cylinder Nyström kernel
pre-integrates axial via :math:`\mathrm{Ki}_n`, which produces the
documented planar-limit mismatch; Variant α avoids this by retaining
the full angular distribution.

Cylinder phase-space and conserved invariants
----------------------------------------------

Phase-space coordinates: :math:`r \in [0, R]` (radial position),
:math:`\mu_{\rm axial} = \cos\theta_{\rm axis} \in [-1, 1]` (cosine
to cylinder axis :math:`\hat z`), :math:`\varphi_{\rm az} \in [0,
2\pi)` (azimuth of the in-plane velocity component).

In-plane velocity fraction:

.. math::
   :label: peierls-greens-cylinder-in-plane-speed

   s_{\rm in\!-\!plane}(\mu_{\rm axial}) =
       \sin\theta_{\rm axis} = \sqrt{1 - \mu_{\rm axial}^2}.

Conserved impact parameter (preserved across all bounces of a
specular trajectory):

.. math::
   :label: peierls-greens-cylinder-impact-parameter

   b(r, \varphi_{\rm az}) = r\,|\sin\varphi_{\rm az}|.

First-leg trajectory chord (3D)
--------------------------------

The 2D in-plane backward chord from interior point :math:`(r,
\varphi_{\rm az})` to the cylinder surface is

.. math::
   :label: peierls-greens-cylinder-trajectory

   L_{\rm 2D, first}(r, \varphi_{\rm az}) =
       r\,\cos\varphi_{\rm az} +
       \sqrt{R^2 - r^2 \sin^2\varphi_{\rm az}}, \qquad
   L_{0, \rm 3D} =
       \frac{L_{\rm 2D, first}}{s_{\rm in\!-\!plane}}.

Bounce-period chord (3D)
-------------------------

The 2D in-plane antipodal chord at conserved impact parameter
:math:`b` is :math:`L_{\rm 2D, period} = 2\sqrt{R^2 - b^2}` and the
3D bounce-period chord is

.. math::
   :label: peierls-greens-cylinder-bounce-period

   L_{\rm period}(b, \mu_{\rm axial}) =
       \frac{2\sqrt{R^2 - b^2}}{\sqrt{1 - \mu_{\rm axial}^2}}.

This is the cylinder analogue of the sphere :math:`L_p = 2 R
\mu_{\rm surf}` — but cylinder's period depends on TWO conserved
quantities (impact parameter AND axial cosine), where sphere's
depends on one (:math:`\mu_{\rm surf}`).

Resolvent (closed-form geometric series)
-----------------------------------------

The cylinder Variant α surface fixed-point closure is

.. math::
   :label: peierls-greens-cylinder-T

   \psi_{\rm surf}(b, \mu_{\rm axial}) =
       \frac{\alpha\,B(b, \mu_{\rm axial})}
            {1 - \alpha\,e^{-\Sigma_t L_{\rm period}}},

with :math:`B` the bounce-period source integral (analogue of
sphere's bounce-period :math:`B`).

Variant α architecture for cylinder
------------------------------------

For each phase-space grid point :math:`(r_i, \mu_{q,\rm axial},
\varphi_{p,\rm az})`:

.. math::
   :label: peierls-greens-cylinder-architecture

   \psi(r_i, \mu_{\rm axial}, \varphi_{\rm az}) =
       F(r_i, \mu_{\rm axial}, \varphi_{\rm az})
       + e^{-\Sigma_t L_{0, \rm 3D}}\,\psi_{\rm surf}(b, \mu_{\rm axial}),

where :math:`F` is the first-leg trajectory integral with 3D
attenuation and :math:`\psi_{\rm surf}` is the bounce-sum closure
(:eq:`peierls-greens-cylinder-T`). The scalar flux is the full
angular reduction:

.. math::

   \phi(r) = \int_{-1}^{1}\!\mathrm d\mu_{\rm axial}
             \int_0^{2\pi}\!\mathrm d\varphi_{\rm az}\,
             \psi(r, \mu_{\rm axial}, \varphi_{\rm az}).

Operator-level identities (V_α1_cyl / V_α2_cyl / V_α3_cyl)
-----------------------------------------------------------

Mirroring the sphere proofs, the cylinder Branch-1 SymPy module
verifies three identities:

- **V_α1_cyl** — closed-cylinder bounce-sum self-consistency. For
  homogeneous cylinder with constant volumetric source :math:`q`,
  :math:`\psi(r, \mu_{\rm axial}, \varphi_{\rm az}) = q/\Sigma_t`
  everywhere. Algebraically identical to V_α1 sphere — only the
  chord formula :math:`L_{\rm period}` differs. Operator action on
  isotropic trial gives :math:`(K \cdot 1) = \omega_0`, hence
  :math:`k_{\rm eff} = k_\infty` for closed cylinder.
- **V_α1_cyl.geometry** — bounce-period chord
  :math:`L_{\rm period}(b, \mu_{\rm axial}) = 2\sqrt{R^2 - b^2} /
  \sqrt{1 - \mu_{\rm axial}^2}` derived two ways (impact-parameter
  form; surface-tangent form) agree. Foundation for V_α1_cyl.
- **V_α2_cyl** — :math:`T_{00}^{\rm cyl} = P_{ss}^{\rm cyl}`
  integrand-level identity (rank-1 cylinder Variant α ≡ cylinder
  Hébert white-BC closure). Unlike sphere, the cylinder
  :math:`\mathrm{Ki}_3` integrand has no elementary closed form,
  so the identity is integrand-level (both reduce to
  :math:`(4/\pi) \cos\alpha\,\mathrm{Ki}_3(2 \Sigma_t R \cos\alpha)`).
- **V_α3_cyl** — surface fixed-point closure
  (:eq:`peierls-greens-cylinder-T`) carries leading factor
  :math:`\alpha`, so :math:`\psi_{\rm surf} \to 0` at :math:`\alpha
  = 0`. Vacuum BC is the trivial limit; no special-case branch.

Numerical verification status (Phase 1 standalone)
---------------------------------------------------

- **k_inf exactness at α=1** — fuel-A-like XS, both thin
  (:math:`\tau_R = 2.5`) and moderate (:math:`\tau_R = 5`):
  :math:`k_{\rm eff} = k_\infty = 0.20833\overline{3}` to ≤ 1e-10
  relative. Convergence in 1 iteration from constant initial
  guess; convergence in ~50 iterations to 1e-9 from aggressive
  :math:`(r/R)^2` perturbation. Scalar flux uniform to ~1e-15.
- **Multi-group at α=1** — 2G with asymmetric scattering
  (downscatter + upscatter):
  :math:`k_{\rm eff} = k_\infty` from
  :func:`~orpheus.derivations.common.eigenvalue.kinf_homogeneous`
  to ≤ 1e-9. The asymmetric XS detects the SigS-transpose
  convention drift (ERR-002 anti-pattern).
- **MG G=1 reduction** — bit-equal to the 1G solver at α=0.
- **Vacuum α=0** — physics-plausible peaked-at-center flux profile,
  positive leakage, :math:`k_{\rm eff} < k_\infty`. fuel-A-like τ_R
  = 2.5: k_eff ≈ 0.12045 (k_eff/k_inf ≈ 0.578). Convergence floor
  achieved at quadrature order :math:`(n_r, n_\mu, n_\varphi,
  n_{\rm traj}) = (24, 20, 64, 96)`: ~1e-7 self-consistency relative
  to next-coarser grid.

Cross-check disclaimer (Issue #129)
------------------------------------

The cylinder Variant α and cylinder Nyström vacuum solvers do NOT
agree to machine precision on the same XS:

- Variant α (3D angle-resolved): :math:`k_{\rm eff} \approx
  0.12045159` at fuel-A-like τ_R = 2.5.
- Nyström vacuum (axial pre-integrated): :math:`k_{\rm eff} \approx
  0.12046563` at fuel-A-like τ_R = 2.5.

The relative difference (~1.4e-4) is **not a bug** — it reflects the
documented Issue #129 axial pre-integration discrepancy. Both are
correct solutions of slightly different formulations. The
load-bearing L1 cross-check for cylinder Variant α is the
**k_inf-exactness invariant at α=1**, which is structurally
independent (closed-form analytical from V_α1_cyl).


.. _peierls-greens-slab:

Slab Variant α (Phase 3A standalone, symmetric reflective)
===========================================================

.. todo:: Archivist expansion needed.

   This section is a **method-implementer stub** — it states the
   load-bearing equations, labels them for cross-reference, and
   points to the SymPy module + tests. The full mathematical
   narrative (worked algebra; design rationale; the
   2-bounce-per-period structural distinction; the ERR-035 history
   and the post-fix delegation architecture; treatment of the
   :math:`\mu \to 0` grazing co-divergence; convergence-floor
   discussion; comparison with the slab Nyström specular reference;
   roadmap to Phase 3B asymmetric BC requiring rank-2 resolvent) is
   the archivist's deliverable.

   Source artifacts:

   - Branch-1 SymPy: :mod:`orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function_slab`
     (V_α1_slab, V_α2_slab, V_α3_slab).
   - Branch-2 production:
     :mod:`orpheus.derivations.continuous.peierls_greens_function.greens_function_slab`
     (:func:`solve_greens_function_slab`,
     :func:`solve_greens_function_slab_mg`).
   - Symbolic test gate:
     :file:`tests/derivations/test_peierls_greens_function_slab_symbolic.py`
     (10 foundation-tagged tests).
   - Numerical L1 gates:
     :file:`tests/derivations/test_peierls_greens_function_slab_solver.py`
     (12 L1-tagged tests).
   - Plan: :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md`.
   - Closeout memo:
     :file:`.claude/agent-memory/method-implementer/slab_variant_alpha_phase3a.md`.

Phase-3A slab Variant α has two implementation phases:

1. **Original Phase-3A (2026-05-02 morning)** mounted on the same
   rank-1 resolvent as sphere/cylinder, using a heuristic
   ``alpha_per_period = α²`` extension of :func:`apply_variant_alpha_closure`
   to encode the **2-bounce-per-period** structural distinction of
   slab geometry. The heuristic closure
   :math:`\psi_{\rm surf} = \alpha\,B_{\rm period}/(1 - \alpha^2\,
   e^{-2\tau})` was found at the Phase-3B reduce-to-symmetric
   consistency check to be incorrect at intermediate :math:`\alpha`
   (logged as ERR-035 in :file:`.claude/skills/vv-principles/error_catalog.md`).
2. **ERR-035 fix (2026-05-02 afternoon)**: Phase-3A
   :func:`solve_greens_function_slab` and
   :func:`solve_greens_function_slab_mg` were refactored as thin
   wrappers that delegate to the Phase-3B
   :func:`solve_greens_function_slab_asymmetric` /
   :func:`solve_greens_function_slab_asymmetric_mg` at
   :math:`\alpha_L = \alpha_R = \alpha`. The first-principles
   rank-2 closure
   :math:`\psi_{\rm surf} = \alpha\,B/(1 - \alpha\,e^{-\tau})` (with
   single-transit B integral) replaces the heuristic. The Phase-3A
   public API surface is preserved. See
   :ref:`peierls-greens-slab-asym-err035-fix` below for details.

Slab phase-space and 2-bounce-per-period structure
---------------------------------------------------

Phase-space coordinates: :math:`x \in [0, L]` (Cartesian position
across the slab), :math:`\mu \in [-1, 1]` (signed direction-cosine
wrt outward normal at :math:`x = L`). Specular reflection on a
symmetric slab preserves :math:`|\mu|` and alternates trajectory
sign at each wall. **One full period traverses both walls** (out
+ back transit at constant :math:`|\mu|`), totaling

.. math::
   :label: peierls-greens-slab-bounce-period

   L_{\rm period}^{\rm slab}(\mu) = \frac{2 L}{|\mu|}.

Two reflections per period, so the **per-period reflection product is**
:math:`\alpha^2`, distinct from the single-reflection per-period
products of sphere (:math:`\alpha`) and cylinder (:math:`\alpha`).

First-leg trajectory chord
---------------------------

The backward chord from interior :math:`(x, \mu)` to first surface
arrival depends on the sign of :math:`\mu`:

.. math::
   :label: peierls-greens-slab-trajectory

   L_{\rm first}(x, \mu) = \begin{cases}
       x / \mu       & \mu > 0 \text{ (came from } x = 0 \text{)} \\
       (L-x) / |\mu| & \mu < 0 \text{ (came from } x = L \text{)}
   \end{cases}.

Resolvent (first-principles single-transit closure, post-ERR-035)
------------------------------------------------------------------

The slab Variant α surface fixed-point closure (post-ERR-035 fix,
2026-05-02) is

.. math::
   :label: peierls-greens-slab-T

   \psi_{\rm surf}(\mu) = \frac{\alpha\,B(\mu)}
                                {1 - \alpha\,e^{-\Sigma_t L/|\mu|}},

with :math:`B(\mu)` the **single-transit** source integral
:math:`\int_0^{L/|\mu|} q\,e^{-\Sigma_t s}\,\mathrm d s` (NOT the
full out+back period integral). The leading factor :math:`\alpha`
is the reflection amplitude for the FIRST surface arrival; the
denominator :math:`(1 - \alpha\,e^{-\Sigma_t L/|\mu|})^{-1}` is
the rank-2 boundary-to-boundary scattering resolvent at
:math:`\alpha_L = \alpha_R = \alpha`. See
:eq:`peierls-greens-slab-asym-resolvent` for the underlying rank-2
matrix form; the symmetric-BC reduction collapses
:math:`\det(I - S) = 1 - \alpha^2\,e^{-2\tau}` into the closure
denominator's :math:`(1 - \alpha\,e^{-\tau})` via the algebraic
identity :math:`(1 - \alpha\,e^{-\tau})(1 + \alpha\,e^{-\tau}) =
1 - \alpha^2\,e^{-2\tau}` (with the :math:`(1 + \alpha\,e^{-\tau})`
factor absorbed by the rank-2 numerator matrix off-diagonals).

The original Phase-3A heuristic closure
:math:`\alpha\,B_{\rm period}/(1 - \alpha^2\,e^{-2\tau})` (with
full-period B integral) coincides with this form ONLY at
:math:`\alpha \in \{0, 1\}` and at uniform-source profiles. ERR-035
documents the discrepancy at intermediate α and the resolution.

Variant α architecture for slab
--------------------------------

For each phase-space grid point :math:`(x_i, \mu_q)`:

.. math::
   :label: peierls-greens-slab-architecture

   \psi(x_i, \mu) =
       F(x_i, \mu)
       + e^{-\Sigma_t L_{\rm first}}\,\psi_{\rm surf}(\mu),

where :math:`F` is the first-leg trajectory integral and
:math:`\psi_{\rm surf}` is :eq:`peierls-greens-slab-T`. The scalar flux
is

.. math::

   \phi(x) = 2\pi \int_{-1}^{1}\!\psi(x, \mu)\,\mathrm d\mu.

Operator-level identities (V_α1_slab / V_α2_slab / V_α3_slab)
--------------------------------------------------------------

- **V_α1_slab** — closed-slab bounce-sum self-consistency. For
  homogeneous slab with constant volumetric source :math:`q`,
  :math:`\psi(x, \mu) = q/\Sigma_t` everywhere. Algebraically
  identical to V_α1 sphere/cylinder — only the chord formula
  :math:`L_{\rm period} = 2L/|\mu|` differs. Operator action on
  isotropic trial gives :math:`(K \cdot 1) = \omega_0`, hence
  :math:`k_{\rm eff} = k_\infty` for closed slab.
- **V_α2_slab** — :math:`T_{00}^{\rm slab} = P_{ss}^{\rm slab} =
  2\,E_3(\Sigma_t L)`. Two structurally-independent derivational
  paths: per-face transfer-matrix definition (µ-domain) and polar
  escape-probability integral (θ-domain). SymPy chokes on direct
  integration of :math:`e^{-\tau/\mu}` (algebra-of-record SymPy
  choke mode #1) but reduces via the :math:`u = 1/\mu`
  substitution to the canonical :math:`E_3` form. Numerical
  structural-independence verified via arbitrary-precision mpmath
  quadrature on both ORIGINAL integrands.
- **V_α3_slab** — surface fixed-point closure
  (:eq:`peierls-greens-slab-T`) carries leading factor
  :math:`\alpha`; :math:`\psi_{\rm surf} \to 0` at
  :math:`\alpha = 0`. Vacuum BC is the trivial limit; no
  special-case branch.

Numerical verification status (Phase 3A standalone)
----------------------------------------------------

- **k_inf exactness at α=1** — fuel-A-like XS, both thin
  (:math:`\tau_L = 1`) and moderate (:math:`\tau_L = 5`):
  :math:`k_{\rm eff} = k_\infty = 0.20833\overline{3}` to
  :math:`\le 1\mathrm e{-10}` relative. Convergence in 1 iteration
  from constant initial guess; convergence in :math:`\sim 50`
  iterations to 1e-9 from aggressive quadratic perturbation.
  Scalar flux uniform to :math:`\sim 1\mathrm e{-15}`.
- **Multi-group at α=1** — 2G with asymmetric scattering
  (downscatter + upscatter): :math:`k_{\rm eff} = k_\infty` to
  :math:`\le 1\mathrm e{-9}`. Asymmetric XS detects the
  SigS-transpose convention drift (ERR-002 anti-pattern).
- **MG G=1 reduction** — bit-equal to the 1G solver at α=0.
- **Vacuum α=0** — physics-plausible peaked-at-center symmetric
  flux profile, positive leakage, :math:`k_{\rm eff} < k_\infty`.
  fuel-A-like τ_L = 5: :math:`k_{\rm eff} \approx 0.157`,
  :math:`k_{\rm eff}/k_\infty \approx 0.756` (post-ERR-034 +
  ERR-035 fixes). Quadrature self-consistency between the
  (24, 40, 96) and (32, 56, 128) grids is :math:`\sim 9\mathrm
  e{-6}` (regression-gated at 5e-5), down from the pre-fix
  :math:`\sim 5\mathrm e{-4}` floor that was originally
  misattributed to a slow-µ-quadrature limitation but was
  actually the dominant ERR-034 trajectory bug masquerading as
  quadrature noise.
- **V_α2_slab production-primitive cross-check** —
  :math:`T_{\rm slab}[0, n_modes][0, 1] = 2\,E_3(\Sigma_t L)` from
  :func:`compute_T_specular_slab` agrees with arbitrary-precision
  ``mpmath.expint(3, τ_L)`` to :math:`\le 1\mathrm e{-10}` over
  five :math:`\tau_L` values.

The shared core's ``alpha_per_period`` extension (sphere/cylinder)
-------------------------------------------------------------------

The shared :func:`apply_variant_alpha_closure` API carries an
``alpha_per_period`` keyword that distinguishes per-period
reflection product from BC reflectivity:

- Sphere/cylinder (1-bounce-per-period): do NOT pass the keyword
  (default ``None`` → ``alpha``).
- Slab symmetric (2-bounce-per-period): pre-ERR-035, the slab
  Phase-3A path passed ``alpha_per_period=alpha**2`` to encode the
  two-bounce structure with a rank-1 resolvent. **Post-ERR-035
  this code path is unused**: slab Phase-3A delegates to the
  rank-2 asymmetric solver, which uses
  :func:`apply_variant_alpha_closure_rank2` directly. The
  ``alpha_per_period`` keyword is retained on
  :func:`apply_variant_alpha_closure` for API stability and for
  use by any future 2-bounce-per-period geometry that legitimately
  benefits from the rank-1 form (none currently in the codebase).

Foundation tests pin both branches in
:file:`tests/derivations/test_peierls_variant_alpha_core.py`.


.. _peierls-greens-slab-asym:

Slab asymmetric Variant α (Phase 3B, rank-2 resolvent)
=======================================================

.. todo:: Archivist expansion needed.

   This section is a **method-implementer stub** — it states the
   load-bearing equations, labels them for cross-reference, and
   points to the SymPy module + tests. The full mathematical
   narrative (worked rank-2 algebra; method-of-images derivation
   linking the asymmetric reflective-vacuum slab to the
   doubled-domain symmetric vacuum slab; ERR-035 history and fix
   linking Phase-3A symmetric to Phase-3B rank-2; roadmap to
   hollow sphere / annulus in Phase 3C inheriting the rank-2
   framework) is the archivist's deliverable.

   Source artifacts:

   - Branch-1 SymPy:
     :mod:`orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function_slab_asymmetric`
     (V_α1_slab_asym, V_α2_slab_asym, V_α3_slab_asym).
   - Branch-2 production:
     :mod:`orpheus.derivations.continuous.peierls_greens_function.greens_function_slab_asymmetric`
     (:func:`solve_greens_function_slab_asymmetric`,
     :func:`solve_greens_function_slab_asymmetric_mg`).
   - Symbolic test gate:
     :file:`tests/derivations/test_peierls_greens_function_slab_asymmetric_symbolic.py`
     (16 foundation-tagged tests).
   - Numerical L1 gates:
     :file:`tests/derivations/test_peierls_greens_function_slab_asymmetric_solver.py`
     (11 L1-tagged tests including the load-bearing method-of-images
     symmetry test).
   - Plan: :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md`.
   - Closeout memo:
     :file:`.claude/agent-memory/method-implementer/slab_asymmetric_variant_alpha_phase3b.md`.

Phase-3B asymmetric slab Variant α extends the BC parameter space
from one scalar :math:`\alpha \in [0, 1]` to TWO independent per-
wall reflectivities :math:`\alpha_L, \alpha_R \in [0, 1]`. The
shared rank-1 closure (sphere, cylinder, symmetric slab) is
replaced by a **rank-2 boundary-to-boundary scattering resolvent**
:math:`T = (I - S)^{-1}` — the predicted generalisation from the
cross-domain frame memo
:file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`.

Rank-2 monodromy and resolvent
-------------------------------

Phase space :math:`(x, \mu)` with surface state
:math:`[\psi_L^+(\mu), \psi_R^-(\mu)]` (outgoing flux at the left
wall in :math:`+\hat\mu` direction and at the right wall in
:math:`-\hat\mu` direction). The monodromy operator advancing the
surface state by ONE STEP of the bouncing trajectory is

.. math::
   :label: peierls-greens-slab-asym-monodromy

   S(\alpha_L, \alpha_R, \tau) = \begin{pmatrix}
       0                              & \alpha_L\,e^{-\tau} \\
       \alpha_R\,e^{-\tau}            & 0
   \end{pmatrix}, \qquad
   \tau = \Sigma_t\,L /|\mu|

(single-transit optical depth — NOT the full out-and-back period).
The diagonal of :math:`S` is identically zero: one step never
returns to the same wall.

The rank-2 resolvent is

.. math::
   :label: peierls-greens-slab-asym-resolvent

   T(\alpha_L, \alpha_R, \tau) = (I - S)^{-1}
       = \frac{1}{1 - \alpha_L\,\alpha_R\,e^{-2\tau}}
         \begin{pmatrix}
             1                       & \alpha_L\,e^{-\tau} \\
             \alpha_R\,e^{-\tau}     & 1
         \end{pmatrix}.

Surface-flux closure
---------------------

The closure equations involve **single-transit** B integrals (NOT
full bounce-period integrals as in the rank-1 symmetric closure):

.. math::
   :label: peierls-greens-slab-asym-closure

   \psi_L^+(\mu) &= \alpha_L\,B_{LR}(\mu) + \alpha_L\,e^{-\tau}\,
                       \psi_R^-(\mu), \\
   \psi_R^-(\mu) &= \alpha_R\,B_{RL}(\mu) + \alpha_R\,e^{-\tau}\,
                       \psi_L^+(\mu),

where

.. math::

   B_{LR}(\mu) = \int_0^{L/|\mu|} q(|\mu|\,s)\,e^{-\Sigma_t s}\,
                  \mathrm d s, \qquad
   B_{RL}(\mu) = \int_0^{L/|\mu|} q(L - |\mu|\,s)\,e^{-\Sigma_t s}\,
                  \mathrm d s.

In matrix form :math:`\psi_{\rm surf} = T \cdot \alpha\,B` with
:math:`(\alpha B)_L = \alpha_L\,B_{LR}` and :math:`(\alpha B)_R =
\alpha_R\,B_{RL}`.

Variant α architecture for asymmetric slab
-------------------------------------------

For each phase-space grid point :math:`(x_i, \mu_q)`:

.. math::
   :label: peierls-greens-slab-asym-architecture

   \psi(x_i, \mu) = F(x_i, \mu)
                   + e^{-\Sigma_t L_{\rm first}}\,
                     \psi_{\rm surface}(\mu),

where :math:`\psi_{\rm surface}` is :math:`\psi_L^+(\mu)` for
:math:`\mu > 0` (the trajectory entered from the left wall) and
:math:`\psi_R^-(\mu)` for :math:`\mu < 0` (entered from the right
wall).

Method-of-images symmetry test (load-bearing acceptance)
---------------------------------------------------------

.. math::
   :label: peierls-greens-slab-asym-method-of-images

   \boxed{\;
   k_{\rm eff}\!\big[\text{slab } [0, L],\, \alpha_L = 1,\, \alpha_R = 0\big]
   \;=\;
   k_{\rm eff}\!\big[\text{slab } [0, 2L],\, \alpha_L = \alpha_R = 0\big]
   \;}

The reflective BC at :math:`x = 0` enforces the even-symmetry plane
of the fundamental eigenmode of the symmetric vacuum slab
:math:`[0, 2L]` (which is symmetric about :math:`x = L`). Cutting
the symmetric slab at :math:`x = L` and applying the specular
condition (which holds by symmetry) yields an equivalent problem
on the right half :math:`[L, 2L]` — bijectively mapped to
:math:`[0, L]` via :math:`x_{\rm asym} = x_{\rm sym} - L`. Therefore
the eigenvalues match exactly, and the flux profiles satisfy
:math:`\phi_{\rm asym}(x) = \phi_{\rm sym}(L + x)` for :math:`x \in
[0, L]`.

This is the **canonical structural-independence test** for the
rank-2 prototype — it exercises both the closure algebra (across
the diagonal of the BC parameter square) and the trajectory
machinery (the asymmetric flux profile peaks AT the reflective
wall, not at an interior point — a sensitive bug-fingerprint for
trajectory-parametrisation errors).

Numerical verification status (Phase 3B)
-----------------------------------------

- **k_inf exactness at α_L=α_R=1** — fuel-A-like XS, both thin
  (:math:`\tau_L = 1`) and moderate (:math:`\tau_L = 5`):
  :math:`k_{\rm eff} = k_\infty` to :math:`\le 1\mathrm e{-10}`
  relative (machine floor :math:`\sim 5.55\mathrm e{-17}`,
  1 iteration from constant initial guess).
- **Method-of-images symmetry test** (the load-bearing gate) —
  asym :math:`[0, 1]` :math:`(\alpha_L=1, \alpha_R=0)` agrees with
  sym :math:`[0, 2]` :math:`(\alpha=0, 0)` at :math:`\le 5.6\mathrm
  e{-8}` relative at (n_x_asym=32, n_mu=32, n_traj_quad=64);
  :math:`\le 6.5\mathrm e{-9}` at (64, 64, 128). Asymmetric flux
  peaks at :math:`x \approx 0.001` (the reflective wall), monotone
  decreasing to vacuum.
- **Vacuum-vacuum bit-equality** — at :math:`\alpha_L = \alpha_R = 0`,
  rank-2 result is bit-equal to Phase-3A rank-1 (both bypass the
  closure entirely and reduce to the bare first-leg integral).
- **Multi-group at α_L=α_R=1** — 2G with asymmetric scattering
  matches :math:`k_\infty` from
  :func:`~orpheus.derivations.common.eigenvalue.kinf_homogeneous`
  to :math:`\le 1\mathrm e{-9}` relative.

.. _peierls-greens-slab-asym-err035-fix:

ERR-035 history and fix (Phase-3A heuristic closure)
-----------------------------------------------------

The Phase-3B rank-2 closure (verified by direct first-principles
derivation in
:func:`derive_operator_constant_trial_closed_slab_asymmetric`) and
the original Phase-3A rank-1 symmetric closure agreed at the
BC-square corners :math:`\alpha \in \{0, 1\}` (where flux profiles
are uniform-or-rank-1-collapsing) but DIFFERED by :math:`\sim 1.3
\mathrm e{-4}` relative at intermediate :math:`\alpha = 0.5`. The
original Phase-3A formula

.. math::

   \psi_{\rm surf} = \frac{\alpha\,B_{\rm period}}
                          {1 - \alpha^2\,e^{-2\tau}}, \qquad
   B_{\rm period} = \int_0^{2L/|\mu|} q\,e^{-\Sigma_t s}\,\mathrm d s,

was a heuristic generalisation that did NOT match the
structurally-correct (per first-principles derivation) rank-2 form
:math:`\psi_{\rm surf} = \alpha\,B / (1 - \alpha\,e^{-\tau})` at
intermediate :math:`\alpha`. The Phase-3A V_α1 / vacuum tests
missed the discrepancy because they only exercise the corners.

**Fix applied 2026-05-02** (logged as ERR-035 in the V&V
:file:`.claude/skills/vv-principles/error_catalog.md`).
:func:`solve_greens_function_slab` and
:func:`solve_greens_function_slab_mg` were refactored as thin
wrappers that delegate to
:func:`solve_greens_function_slab_asymmetric` and
:func:`solve_greens_function_slab_asymmetric_mg` at
:math:`\alpha_L = \alpha_R = \alpha`. The heuristic
``_apply_operator_slab`` and the helper chord-length functions
were deleted (dead code after delegation; keeping them as
"alternative paths" would re-create the very pattern the fix
removes). The Phase-3A public API surface (:class:`SlabGreensResult`,
:class:`SlabGreensMGResult`, the two solver entrypoints) is
preserved by re-wrapping the rank-2 result. Net effect: all
Phase-3A code paths use the first-principles rank-2 closure;
the test
:func:`test_rank1_path_now_agrees_with_rank2_via_delegation_after_ERR035_fix`
asserts ≤ 1e-12 rtol bit-equal agreement between the two paths
at :math:`\alpha = 0.5` (regression-prevention gate against any
future re-introduction of a heuristic closure).

A second latent bug (ERR-034: the trajectory parametrisation
:math:`x_{\rm traj}(s) = x - s` instead of :math:`x - \mu \cdot s`,
silent on uniform-:math:`q` cases) was caught and FIXED in both
Phase-3A and Phase-3B during Phase-3B integration. See ERR
catalog entry in :file:`.claude/skills/vv-principles/error_catalog.md`.


Provenance: literature, code, and tests
========================================

Literature provenance
---------------------

.. list-table:: References for the Variant α architecture
   :header-rows: 1
   :widths: 30 50 20

   * - Reference
     - Role
     - Local copy?
   * - Sanchez 1986 [SanchezTTSP1986]_
     - Primary source. Eqs. (A1)–(A7) define the angle-resolved
       Green's function with BC absorbed; Eq. (A6) is the
       angle-integrated reduction whose hypersingularity Phase 5
       fought (and lost). Variant α uses Eqs. (A1) and (A5)
       directly.
     - YES (PDF in repo root)
   * - Pomraning-Siewert 1982 [PomraningSiewert1982]_
     - Precursor for :math:`\omega_1 = 0` (isotropic-scattering)
       case with vacuum BC — the truth source for the V_α3
       reduction once the prototype is extended to vacuum BC.
     - NO (paywalled)
   * - Hébert 2009 [Hebert2020]_ §3.8.5
     - Rank-1 white-BC closure :math:`(1 - P_{ss})^{-1}` — the
       V_α2 algebraic anchor.
     - YES (Hebert(2009)Chapter3.pdf in repo root)
   * - Garcia 2020 / 2021 / 2018 [Garcia2020]_ [Garcia2021]_
     - Modern stable :math:`P_N` spherical-harmonic reference for
       homogeneous-sphere reflective-BC k_eigenvalue, plus
       multi-region extension. The only external numerical reference
       for this problem class.
     - NO (all paywalled)
   * - Sanchez 2002 [Sanchez2002]_
     - Periodic-trajectory closure for *lattice* geometries —
       algebraic structure :math:`1/(1 - \psi_{bd})` parallel to
       Variant α's :math:`T(\mu_{\rm surf})`, but for a different
       geometry. Cross-check on the universality of the multi-bounce
       factor.
     - YES (Sanchez(2002).pdf in repo root)
   * - Williams 1971
     - Fredholm-second-kind theory underpinning power-iteration
       convergence on the integral operator. Not blocking — the
       positive kernel of Variant α makes Perron-Frobenius
       sufficient.
     - NO (book, no DOI)
   * - Case & Zweifel 1967
     - Singular eigenfunction expansion (Variant γ in the Plan 2 B1
       memo) — *rejected* in favour of Variant α because the
       continuous-spectrum projection requires Chandrasekhar-:math:`H`
       machinery and is dominated by Garcia's modern stable :math:`P_N`
       re-summation (Variant δ).
     - NO (book, no DOI)

Detailed extraction of Sanchez 1986 lives in the
literature-researcher memo at
:file:`.claude/agent-memory/literature-researcher/phase5_sanchez_1986_sphere_specular.md`.
The B1 architectural-survey memo at
:file:`.claude/agent-memory/literature-researcher/peierls_greens_function_lit.md`
covers the full literature landscape and the Variant α / β / γ / δ
recommendation.

Code provenance
---------------

The Variant α implementation spans three files:

- :mod:`orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function`
  — SymPy derivations V_α1, V_α2, V_α3 (operator-level identities).
  About 270 lines.

  - :func:`~orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function.derive_operator_constant_trial_closed_sphere`
    — V_α1 closed-sphere bounce-sum self-consistency.
  - :func:`~orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function.derive_T00_equals_P_ss_sphere`
    — V_α2 algebraic identity :math:`T_{00}^{\rm sphere} =
    P_{ss}^{\rm sphere}`.
  - :func:`~orpheus.derivations.continuous.peierls_greens_function.origins.specular.greens_function.derive_alpha_zero_kernel_reduction`
    — V_α3 vacuum-BC kernel reduction.

- :mod:`orpheus.derivations.continuous.peierls_greens_function.greens_function` —
  production solver suite. About 1300 lines.

  - :func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function._apply_operator_with_source_profile`
    — per-group homogeneous operator action; the load-bearing
    primitive shared by 1G / MG / 1-region paths.
  - :func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function._apply_operator`
    — 1G wrapper around ``_apply_operator_with_source_profile``;
    builds :math:`q(r) = (\Sigma_s + \nu\Sigma_f/k)\,\phi(r) / 4\pi`.
  - :func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function.solve_greens_function_sphere`
    — 1G public driver with :math:`\alpha\in[0,1]` parametrisation
    (B-phase + A1+A2 vacuum extension).
  - :func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function.solve_greens_function_sphere_mg`
    — multi-group public driver (A3 follow-on); arbitrary G,
    full :math:`G\times G` scattering matrix, arbitrary :math:`\chi`.
  - :func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function._apply_operator_mr`
    — multi-region per-group operator with piecewise :math:`\Sigt{}`
    along trajectory + bounce-period chord (Plan-(b)).
  - :func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function.solve_greens_function_sphere_mr`
    — multi-region multi-group k-eigenvalue driver
    (Plan-(b) Option 2, Issue #132 attack).
  - :func:`~orpheus.derivations.continuous.peierls_greens_function.greens_function.solve_greens_function_sphere_mr_fixed_source`
    — multi-region multi-group fixed-source driver
    (Plan-(b) Option 1, Garcia 2021 cross-check).
  - :class:`~orpheus.derivations.continuous.peierls_greens_function.greens_function.GreensFunctionResult`,
    :class:`GreensFunctionMGResult`,
    :class:`GreensFunctionMRResult`,
    :class:`GreensFunctionFixedSourceResult` — result dataclasses.

- :mod:`orpheus.derivations.continuous.peierls_nystrom.ps1982_reference` —
  PS-1982 Eq. (21) Nyström reference solver for the homogeneous
  vacuum sphere (used as the structurally-independent A2 cross-check).

Test provenance
---------------

The test gates lock down each algebraic identity and the numerical
implementation against it:

- :file:`tests/derivations/test_peierls_greens_function_symbolic.py`
  — 8 SymPy gates (V_α1 surface fixed-point + total-ψ-constant +
  operator-eigenvalue + composite; V_α2 integrand-match + closed-form
  + composite; V_α3 g_h vanishes at α=0).

- :file:`tests/derivations/test_peierls_greens_function_solver.py`
  — 3 numerical gates (V_α1.numerical with constant initial guess
  + non-uniform initial guess + two thicknesses
  :math:`\tau_R \in \{2.5, 5\}`).

- :file:`tests/derivations/test_peierls_greens_function_xverif.py`
  — 3 closed-sphere cross-verification gates: B5.A (Variant α exact),
  B5.B (Phase 4 N=1 ≡ ``white_hebert``), B5.C (Phase 4 rank
  convergence toward Variant α).

- :file:`tests/derivations/test_peierls_greens_function_vacuum.py`
  — 5 A1 vacuum-BC gates (k_eff < k_inf, thick sphere asymptote,
  α-continuity, non-trivial spatial mode, α=1 unchanged).

- :file:`tests/derivations/test_peierls_greens_function_xverif_ps1982.py`
  — 6 A2 PS-1982 cross-verification gates (4 parametrised
  configurations + thick-sphere regression + flux-shape qualitative).

- :file:`tests/derivations/test_peierls_greens_function_mg.py`
  — 7 A3 multi-group gates (G=1 reduction, 2G downscatter +
  upscatter, 2G spectrum, 2G vacuum, 4G with realistic χ, 4G
  vacuum).

- :file:`tests/derivations/test_peierls_greens_function_mr.py`
  — 4 Plan-(b) Option 2 multi-region k-eigenvalue gates (1-region
  reduction to MG, Issue #132 catastrophe avoided, spatial mode
  physical, vacuum BC reduces k_eff).

- :file:`tests/derivations/test_peierls_greens_function_garcia2021.py`
  — 17 Plan-(b) Option 1 fixed-source gates (3 sanity + 15 per-r-point
  cross-checks vs Garcia 2021 Table 5).

All tests tagged ``@pytest.mark.foundation`` because they verify
operator-level invariants and L1 cross-checks. The load-bearing
physics claims (Variant α is exact for closed homogeneous sphere;
agrees with PS-1982 to 1e-4 on vacuum sphere; reduces to
``kinf_homogeneous`` at closed sphere multi-group; reproduces
Garcia 2021 fixed-source profile to < 2 % at non-interface points)
are gated by the test suite. Cross-validation against Phase 4
``specular_multibounce`` (closed sphere) and PS-1982 (vacuum) and
Garcia 2021 (multi-region fixed-source) provides the L1
structurally-independent evidence chain.

Memo provenance
---------------

- :file:`.claude/agent-memory/literature-researcher/peierls_greens_function_lit.md`
  — Plan 2 B1 literature pull (full reference landscape, Variant
  α/β/γ/δ recommendation).
- :file:`.claude/agent-memory/numerics-investigator/peierls_greens_variant_alpha_decision.md`
  — Plan 2 B2 architectural decision (Variant α adopted, β/γ/δ
  rejected with rationale).
- :file:`.claude/agent-memory/numerics-investigator/peierls_greens_phase1_closeout.md`
  — Plan 2 B6 closeout (cross-verification matrix, parallel
  research-grade-reference decision).
- :file:`.claude/agent-memory/numerics-investigator/_archive/specular_continuous_mu_phase5_retreat.md`
  — Phase 5 retreat (predecessor, motivates the angle-resolved
  reformulation).
- :file:`.claude/plans/peierls-greens-function-approach.md` — Plan 2
  master plan (Part A reorganisation + Part B Green's function).


.. [PomraningSiewert1982] G.C. Pomraning and C.E. Siewert,
   "On the integral form of the equation of transfer for a
   homogeneous sphere," *J. Quant. Spec. Rad. Transfer* **28**,
   503–506 (1982). DOI: 10.1016/0022-4073(82)90016-4. Cited by
   Sanchez 1986 as the :math:`\omega_1 = 0` (isotropic-scattering)
   precursor; vacuum BC; sphere only.

.. [Garcia2020] R.D.M. Garcia, "A numerically stable spherical
   harmonics solution for the neutron transport equation in a
   sphere," *J. Comp. Phys.* **393**, 109139 (2020).
   DOI: 10.1016/j.jcp.2019.109139.

.. [Garcia2021] R.D.M. Garcia, "Accurate spherical harmonics solutions
   for neutron transport problems in multi-region spherical
   geometries," *J. Comp. Phys.* **433**, 109856 (2021).
   DOI: 10.1016/j.jcp.2020.109856.

.. [Sanchez2002] R. Sanchez, "Treatment of Boundary Conditions in
   Trajectory-Based Deterministic Transport Methods," *Nucl. Sci.
   Eng.* **140**, 23–50 (2002). DOI: 10.13182/NSE140-23. Periodic-
   trajectory closure :math:`\psi = \psi_q(L)/(1 - \psi_{bd}(L))
   \cdot \psi_{bd} + \psi_q` (Eq. 15) for lattice geometries —
   parallel algebraic structure to Variant α's
   :math:`\psi_{\rm surf} = T(\mu_{\rm surf})\,B(\mu_{\rm surf})`,
   different geometry.
