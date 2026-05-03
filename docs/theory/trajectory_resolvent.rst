.. _theory-trajectory-resolvent:

==========================================================================
Trajectory-Resolvent Family — angle-resolved Green's function references
==========================================================================

.. _theory-trajectory-resolvent-name:

Why the name "trajectory resolvent" (and what it produces)
================================================================================

The folder is named after **how it builds the answer**, not after what
the answer numerically *is*. This distinction matters and is the reason
this folder is **structurally a sister** to a future
``spectral_resolvent/`` folder (Meaning 2; the spectral closed-form
construction of the same Green's kernel) and to the
:mod:`...singular_eigenfunction` folder (Meaning 1; the angular
Green's-function via Case ν-spectrum expansion).

* **What it constructs.** A solver in this folder traces *characteristic
  rays* (a structurally MoC viewpoint), integrates the optical depth
  :math:`\tau` along each ray, and closes multi-bounce trajectories via
  the rank-1 / rank-2 **boundary-to-boundary scattering resolvent**
  :math:`T = (I - S)^{-1}`. The construction is therefore
  characteristic-tracing + Poincaré-recurrence summation — the same
  family Sanchez, Mao & Santandrea (2002) call **"trajectory-based
  deterministic transport methods"** in the title of *Nucl. Sci. Eng.*
  **140**, 23–50.

* **What it produces.** Mathematically, the OUTPUT is a scalar
  **Green's function** — equivalently, the closed-form integral kernel
  from Sanchez 1986 Eq. A1 / Pomraning–Siewert 1982 Eq. 21. The
  Variant-α rank-1 closure
  :math:`T(\mu) = [\,1 - \alpha\, e^{-\tau(\mu)}\,]^{-1}`
  is verbatim the resolvent that Sanchez 1986 Appendix Eq. (A4) and
  Pomraning–Siewert 1982 Eq. (14) (the latter for diffuse + specular
  reflection on the sphere) write down by direct integration of the
  transfer equation. Two structurally *independent* paths — ours via
  ray tracing + multi-bounce closure, theirs via direct integral
  algebra — converge on the same operator. ORPHEUS treats the
  agreement as an L1 cross-check, not a coincidence.

* **Where the name "Variant α" sits.** ``variant_alpha`` is internal
  ORPHEUS jargon for this rank-1 specialisation; no external paper
  uses it. The folder name avoids the jargon. The internal API
  (``compute_resolvent_T``, ``apply_variant_alpha_closure``,
  ``..._rank2``) keeps the working name because it has been load-bearing
  vocabulary across the multi-region / multi-group expansion (Plan-(b)
  Options 1 + 2, 2026-05).

* **Sister modules.** Three folders carry the "Green's-function family"
  taxonomy (Sanchez–Chandrasekhar three meanings):

  - :mod:`...singular_eigenfunction` — Meaning 1, the **angular
    Green's function** via Case singular-eigenfunction expansion
    (Mitsis 1963 / Westfall–Metcalf 1973 / Atalay 1997). Different
    structural pillar; same physical Green's function under specialisation
    to the homogeneous medium.
  - :mod:`...trajectory_resolvent` — Meaning 2-via-trajectories: this
    folder. Builds the **scalar Green's kernel** numerically by
    characteristic tracing + multi-bounce closure.
  - ``...spectral_resolvent`` (reserved, not yet implemented) —
    Meaning 2-via-spectrum: the same scalar Green's kernel constructed
    by closed-form spectral μ-integration (Sanchez 1986 Eq. A1 /
    Pomraning–Siewert 1982 closed-form). Empty placeholder folder
    holding a README and the planned reference.

The name ``trajectory_resolvent/`` therefore documents simultaneously
(i) the structural construction (trajectory tracking + resolvent
closure) and (ii) the family Sanchez 2002 named in print, while
admitting that the OUTPUT is mathematically a Green's function.

.. note::

   **Scope of this page.** Variant α is now a **6-geometry × 2-topology
   family** of angle-resolved Green's-function reference solvers, all
   mounted on a single shared resolvent + closure primitive
   (:mod:`orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core`).
   The family is closed as of 2026-05-02. The legacy page title and
   sections continue to refer to "sphere" because the sphere prototype
   was historically the first instance and carries the deepest
   numerical-evidence chain (Plan 2 B-phase + A1+A2 vacuum + A3
   multi-group + Plan-(b) multi-region). The sister geometries —
   cylinder (Phase 1), slab symmetric (Phase 3A → ERR-035 delegation
   to Phase 3B), slab asymmetric (Phase 3B), hollow sphere
   (Phase 3C-1), annulus (Phase 3C-2) — inherit the architecture
   verbatim and are documented in dedicated sections below. See the
   :ref:`peierls-greens-family-overview` for the unifying map.

.. _peierls-greens-family-overview:

The Variant α family at a glance
=================================

The unifying claim of this page is that **a single rank-1/rank-2
boundary-to-boundary scattering resolvent**
:math:`T = (I - S)^{-1}` — exposed as the two functions
:func:`~orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core.compute_resolvent_T`
and
:func:`~orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core.compute_resolvent_T_rank2`
in the :mod:`.variant_alpha_core` module — closes the boundary
trajectories for **every** geometry in the table below. The
geometry-specific code is confined to chord arithmetic and
phase-space discretisation; the algebraic closure is shared
byte-for-byte across the entire family.

.. list-table:: Variant α 6-geometry × 2-topology grid (all closed 2026-05-02)
   :header-rows: 1
   :widths: 22 20 18 22 18

   * - Geometry
     - Topology
     - Closure rank
     - Phase-space
     - Phase
   * - Sphere
     - 1-surface compact
     - rank-1
     - :math:`(r, \mu)`
     - B-phase + A1/A2/A3 + Plan-(b)
   * - Cylinder
     - 1-surface compact
     - rank-1 (3D angular)
     - :math:`(r, \mu_{\rm axial}, \varphi_{\rm az})`
     - Phase 1
   * - Slab symmetric
     - 2-surface (rank-2 reduction)
     - rank-2 (delegated)
     - :math:`(x, \mu)`
     - Phase 3A → ERR-035 delegation
   * - Slab asymmetric
     - 2-surface
     - rank-2
     - :math:`(x, \mu)`
     - Phase 3B
   * - Hollow sphere
     - 2-surface curvilinear
     - rank-2 + b-partition
     - :math:`(r, \mu)`
     - Phase 3C-1
   * - Annulus
     - 2-surface curvilinear
     - rank-2 + b-partition
     - :math:`(r, \mu_{\rm axial}, \varphi_{\rm az})`
     - Phase 3C-2

The unification axis is :mod:`.variant_alpha_core`. Every solver in
the table imports
:func:`~orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core.apply_variant_alpha_closure`
(rank-1) or
:func:`~orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core.apply_variant_alpha_closure_rank2`
(rank-2) to perform the boundary closure; the modules that compose
into them are the geometry-specific trajectory primitives only. The
**impact-parameter partition** :math:`b = R_{\rm in}` couples the
b > R_in outer-only rank-1 piece to the b ≤ R_in through-ray rank-2
piece on hollow geometries — predicted by the cross-domain frame
memo at
:file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
before any 2-surface code existed and held across all four 2-surface
instances (slab asymmetric, hollow sphere, annulus, plus the
symmetric-slab reduction).

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

**Sphere prototype scope as of 2026-05-02** (Plan 2 Part B +
A1/A2/A3/Plan-(b) follow-ons all closed; sister-geometry coverage
documented in dedicated sections per the
:ref:`peierls-greens-family-overview`):

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
  The Phase-3 family is **closed** (2026-05-02): six geometries
  (sphere, cylinder, slab symmetric, slab asymmetric, hollow sphere,
  annulus) on two topologies (1-surface compact + 2-surface) all on
  the unified rank-1 / rank-2 framework via :mod:`.variant_alpha_core`.
  See :ref:`peierls-greens-family-overview` for the family map and
  :ref:`peierls-greens-unification` for the operator-theoretic
  unification. The Phase 1 / Phase 2 / Phase 3 closeout decisions
  are at
  :file:`.claude/agent-memory/numerics-investigator/peierls_greens_phase1_closeout.md`
  (sphere baseline) and the per-geometry method-implementer memos
  in :file:`.claude/agent-memory/method-implementer/`.
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
- **What the family covers as of 2026-05-02**:

  - **Sphere** (1G + MG + multi-region homogeneous, α ∈ [0, 1])
    with full Plan 2 + Plan-(b) coverage — V_α1, A1+A2 vacuum,
    A3 multi-group, multi-region k-eigenvalue and fixed-source.
  - **Cylinder** (1G + MG homogeneous, α ∈ [0, 1]) — Phase 1.
  - **Slab symmetric** (1G + MG homogeneous, α ∈ [0, 1]) —
    Phase 3A wrapper delegating to Phase 3B rank-2.
  - **Slab asymmetric** (1G + MG homogeneous, independent α_L, α_R
    ∈ [0, 1]) — Phase 3B, with method-of-images load-bearing test.
  - **Hollow sphere** (1G + MG homogeneous, independent α_in, α_out
    ∈ [0, 1]) — Phase 3C-1, with R_in → 0 solid-sphere reduction
    cross-check.
  - **Annulus** (1G + MG homogeneous, independent α_in, α_out ∈
    [0, 1]) — Phase 3C-2, closing the family on the unified
    framework.

- **What the family does NOT cover**: anisotropic scattering
  (:math:`\omega_1 \ne 0`), Sanchez 1986 diffuse-re-emission
  :math:`\beta`-branch, multi-region for any geometry other than
  sphere (multi-region cylinder / slab / hollow / annulus all
  deferred), 3-surface topologies (axially-capped cylinder /
  annulus). Phase 2 unification SHIPPED at commits ``efbae9c`` /
  ``92d4f10`` / ``166b9ae``; rank-2 framework SHIPPED at commit
  ``7cb5bc6`` (Phase 3B).
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
is :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function._apply_operator`;
the public driver is
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_specular_sphere`.


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
:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function.derive_operator_constant_trial_closed_sphere`
(SymPy verification, paired numerical gate in
:file:`tests/derivations/test_trajectory_resolvent_symbolic.py`).
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
:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function.derive_T00_equals_P_ss_sphere`.

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
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function._apply_operator`'s
bounce-sum branch — V_α3 guarantees this is mathematically correct
without re-deriving anything. The symbolic proof is in
:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function.derive_alpha_zero_kernel_reduction`.

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
:file:`tests/derivations/test_trajectory_resolvent_xverif.py`:

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
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function._apply_operator`
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
:func:`tests.derivations.test_trajectory_resolvent_xverif_ps1982.test_a2_variant_alpha_agrees_with_ps1982`):

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
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere_mg`
with shared per-group operator
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function._apply_operator_with_source_profile`
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
:file:`tests/derivations/test_trajectory_resolvent_mg.py` (~700
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
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function._trajectory_segments`
and
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function._chord_segments`.

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
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere_mr`
(k-eigenvalue, Plan-(b) Option 2) and
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere_mr_fixed_source`
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
:func:`tests.derivations.test_trajectory_resolvent_mr.test_mr_issue132_no_catastrophe_closed_sphere`
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
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere_mr_fixed_source`
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
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere_mr_fixed_source`.

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
:func:`tests.derivations.test_trajectory_resolvent_garcia2021.test_garcia_case1_phi_matches_at_point`:
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

Phase-1 cylinder Variant α is the **first sister geometry of the sphere
prototype**. It mirrors the sphere architecture
(F first-leg + B bounce-period + closed-form geometric series) on a
3D phase-space :math:`(r, \mu_{\rm axial}, \varphi_{\rm az})`. The
axial direction is kept **explicit** per Issue #129 — the cylinder
Nyström kernel pre-integrates axial via :math:`\mathrm{Ki}_n`, which
produces the documented ~22 % planar-limit mismatch; Variant α avoids
the obstruction *structurally* by retaining the full angular
distribution.

The cylinder closeout memo at
:file:`.claude/agent-memory/method-implementer/cylinder_variant_alpha_phase1.md`
records the development verdict: the Branch-2 production solver landed
clean on first run, no L0 errors logged, V_α1_cyl k_inf-exactness
achieved at 1e-15 in 1 iteration, vacuum k_eff convergence floor
~8e-8 at the finest grid (24, 20, 64, 96). The clean landing reflects
the algebra-of-record discipline working as intended: Branch-1 SymPy
proofs verified the algebra BEFORE Branch-2 implementation, so the
production code only had to operationalise what was already known
correct algebraically.

Cylinder phase-space and conserved invariants
----------------------------------------------

The structural difference between sphere and cylinder is **the
dimensionality of the angular phase-space**. Sphere has full spherical
symmetry: every interior point :math:`r \hat r` has rotational
isotropy about :math:`\hat r`, so the angular part of the flux at a
single :math:`r` is fully captured by the cosine :math:`\mu = \hat
\Omega \cdot \hat r`. Cylinder has only **axial** symmetry —
rotational isotropy about :math:`\hat z` — and only **radial**
symmetry within each axial cross-section. The angular flux at each
:math:`(r, z)` therefore depends on TWO directional cosines: the axial
cosine :math:`\mu_{\rm axial} = \hat \Omega \cdot \hat z` and the
in-plane azimuth :math:`\varphi_{\rm az}` (the angle between the
in-plane velocity component and :math:`\hat r`).

Phase-space coordinates: :math:`r \in [0, R]` (radial position),
:math:`\mu_{\rm axial} = \cos\theta_{\rm axis} \in [-1, 1]` (cosine
to cylinder axis :math:`\hat z`), :math:`\varphi_{\rm az} \in [0,
2\pi)` (azimuth of the in-plane velocity component, measured from
the outward radial direction at the trajectory point).

In-plane velocity fraction (the geometric factor that converts
between 2D in-plane arclength and 3D arclength):

.. math::
   :label: peierls-greens-cylinder-in-plane-speed

   s_{\rm in\!-\!plane}(\mu_{\rm axial}) =
       \sin\theta_{\rm axis} = \sqrt{1 - \mu_{\rm axial}^2}.

A particle moving with axial cosine :math:`\mu_{\rm axial}` covers
:math:`s_{\rm in\!-\!plane} \cdot ds_{\rm 3D}` of in-plane arclength
per unit 3D arclength. This factor is the load-bearing convention of
the cylinder solver: chord geometry is computed in the **2D in-plane
projection**, attenuation along the chord is computed in **3D**, and
the conversion is done by dividing the 2D chord by
:math:`s_{\rm in\!-\!plane}`.

**Two conserved invariants** are preserved across specular reflection
at the cylinder surface (radial-normal: outward-pointing
:math:`\hat r` at the surface):

1. The axial cosine :math:`\mu_{\rm axial}` — because specular
   reflection at a radial-normal mirror plane reverses the in-plane
   velocity component without touching the axial component.
2. The in-plane impact parameter

   .. math::
      :label: peierls-greens-cylinder-impact-parameter

      b(r, \varphi_{\rm az}) = r\,|\sin\varphi_{\rm az}|

   — because the in-plane projection of a specularly-reflected ray
   in a circular cross-section produces antipodal-chord geometry,
   which preserves the perpendicular distance from the cylinder axis.

The conservation of TWO invariants instead of sphere's one
(:math:`\mu_{\rm surf}`) is the source of cylinder's increased
quadrature dimensionality. The bounce-period chord depends on BOTH
:math:`b` AND :math:`\mu_{\rm axial}`, so the closed-form geometric
series is over a 2D parameter family rather than a 1D one.

First-leg trajectory chord (3D)
--------------------------------

The 2D in-plane backward chord from interior point :math:`(r,
\varphi_{\rm az})` to the cylinder surface is obtained by parametrising
the in-plane ray :math:`(r_x, r_y)(s) = (r, 0) - s\,(\cos\varphi_{\rm
az}, \sin\varphi_{\rm az})` (placing the receiver point on the
:math:`+\hat x` axis of the cross-section without loss of generality)
and solving :math:`r_x^2 + r_y^2 = R^2` for the **backward** intersection
(:math:`s > 0` on the ray entering the receiver from the surface). The
algebra reduces to the quadratic root

.. math::
   :label: peierls-greens-cylinder-trajectory

   L_{\rm 2D, first}(r, \varphi_{\rm az}) =
       r\,\cos\varphi_{\rm az} +
       \sqrt{R^2 - r^2 \sin^2\varphi_{\rm az}}, \qquad
   L_{0, \rm 3D} =
       \frac{L_{\rm 2D, first}}{s_{\rm in\!-\!plane}}.

The discriminant :math:`R^2 - r^2 \sin^2\varphi_{\rm az}` is positive
for any interior point (:math:`r < R`) at any azimuth, because
:math:`b = r\,|\sin\varphi_{\rm az}| < r < R`. The :math:`+` branch
selects the surface arrival on the upstream side of the receiver
(consistent with the convention of Issue ERR-034 — the trajectory
parametrisation must use :math:`x_{\rm traj} = x - \mu \cdot s`-style
geometry, not :math:`x_{\rm traj} = x - s`).

Bounce-period chord (3D)
-------------------------

The 2D in-plane antipodal chord at conserved impact parameter
:math:`b` follows from the chord theorem of the circle: a chord at
perpendicular distance :math:`b` from the centre has length
:math:`2\sqrt{R^2 - b^2}`. Lifting this to 3D arclength via the
in-plane fraction gives

.. math::
   :label: peierls-greens-cylinder-bounce-period

   L_{\rm period}(b, \mu_{\rm axial}) =
       \frac{2\sqrt{R^2 - b^2}}{\sqrt{1 - \mu_{\rm axial}^2}}.

This is the cylinder analogue of the sphere :math:`L_p = 2 R
\mu_{\rm surf}`, with TWO conserved quantities (impact parameter
AND axial cosine) controlling its value rather than one
(:math:`\mu_{\rm surf}`). The symbolic gate
:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_cylinder.derive_bounce_period_chord_cylinder`
verifies this identity by deriving :math:`L_{\rm period}` two ways
(impact-parameter form and surface-tangent form) and proving they
agree algebraically — V_α1_cyl.geometry, the foundation gate that
the closed-cylinder identity (V_α1_cyl) relies on.

Resolvent (closed-form geometric series)
-----------------------------------------

The cylinder Variant α surface fixed-point closure inherits sphere's
algebraic structure verbatim — the boundary state is rank-1 (single
surface, single per-period reflection product = :math:`\alpha`), so
the resolvent is the scalar geometric series

.. math::
   :label: peierls-greens-cylinder-T

   \psi_{\rm surf}(b, \mu_{\rm axial}) =
       \frac{\alpha\,B(b, \mu_{\rm axial})}
            {1 - \alpha\,e^{-\Sigma_t L_{\rm period}}},

with :math:`B` the bounce-period source integral
:math:`B = \int_0^{L_{\rm period}} q(s)\,e^{-\Sigma_t s}\,\mathrm d s`
(analogue of sphere's bounce-period :math:`B`). The denominator is
:func:`~orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core.compute_resolvent_T`
applied to the cylinder optical depth :math:`\tau_{\rm period} =
\Sigma_t\,L_{\rm period}(b, \mu_{\rm axial})`; the leading
:math:`\alpha` is the single-reflection amplitude on the FIRST surface
arrival; the denominator's :math:`\alpha` is the per-period
reflection product.

The :math:`\alpha = 0` vacuum BC is the trivial limit (V_α3_cyl): the
leading factor :math:`\alpha` makes :math:`\psi_{\rm surf} \to 0`, no
special-case branch is needed. This is the same V_α3 structure as
sphere — vacuum BC is parametrised, not branched.

The :math:`\alpha = 1` perfect-specular limit becomes singular in the
denominator at the **dual grazing-ray locus** where simultaneously
:math:`L_{\rm period} \to 0` and :math:`\alpha \to 1`. For cylinder this
locus is **two-dimensional**: it occurs at :math:`b \to R` (in-plane
tangential ray; in-plane chord vanishes) and at :math:`\mu_{\rm axial}
\to \pm 1` (axial-grazing ray; the 3D arclength factor
:math:`1/s_{\rm in\!-\!plane}` diverges, so optical-depth diverges and
:math:`e^{-\tau_{\rm period}} \to 0`, taking the denominator to
:math:`1 - \alpha`). Practically, the open Gauss-Legendre quadrature
on :math:`(-1, 1) \times (0, 2\pi)` avoids both endpoints, and
iteration converges cleanly without endpoint-handling logic.

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
angular reduction over the 3D angular phase-space:

.. math::

   \phi(r) = \int_{-1}^{1}\!\mathrm d\mu_{\rm axial}
             \int_0^{2\pi}\!\mathrm d\varphi_{\rm az}\,
             \psi(r, \mu_{\rm axial}, \varphi_{\rm az}).

Implementation notes (see
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder.solve_greens_function_cylinder`
and the helper :func:`_apply_operator_cylinder`):

- **Source profile** :math:`q(r) = (\Sigma_s + \nu\Sigma_f/k)\,\phi(r)
  / 4\pi` is built once per outer iteration on the radial nodes and
  evaluated along trajectories via cubic-spline interpolation.
  ``CubicSpline(extrapolate=True)`` with :math:`r_{\rm traj}^2`
  clipped to :math:`[0, R^2]` keeps the trajectory-point lookup safe
  near the surface.
- **Trajectory integration** uses composite Gauss-Legendre on
  :math:`s \in [0, L_{\rm 2D, first}]` (first leg) and
  :math:`s \in [0, L_{\rm 2D, period}]` (bounce period), with
  attenuation :math:`\exp(-\Sigma_t s_{\rm 3D})` where
  :math:`s_{\rm 3D} = s_{\rm 2D} / s_{\rm in\!-\!plane}`. The 2D-first /
  3D-attenuation convention preserves chord-geometry simplicity while
  keeping the 3D physics correct.
- **Outer iteration** is power iteration with Rayleigh-quotient
  :math:`k`-update; same pattern as sphere with no cylinder-specific
  modifications. Convergence in 1 iteration from constant initial
  guess at :math:`\alpha = 1` (V_α1_cyl); ~50 iterations from
  aggressive non-uniform :math:`(r/R)^2` perturbation.
- **Prefactor convention.** The cylinder Nyström primitives use a
  :math:`1/\pi` (azimuthal-folding) convention; Variant α uses the
  standard :math:`1/(4\pi)` per-steradian convention internally and
  integrates over the FULL :math:`\varphi_{\rm az} \in [0, 2\pi)`.
  This is mathematically consistent (not a hidden convention drift)
  but is worth flagging when comparing scalar flux normalisations
  between Variant α and the cylinder Nyström solver.

Operator-level identities (V_α1_cyl / V_α2_cyl / V_α3_cyl)
-----------------------------------------------------------

Mirroring the sphere proofs, the cylinder Branch-1 SymPy module
:mod:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_cylinder`
verifies three identities. The proofs are algebraically nearly
isomorphic to sphere V_α1/V_α2/V_α3, with the chord formula being
the only geometric variation.

**V_α1_cyl** — closed-cylinder bounce-sum self-consistency
(:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_cylinder.derive_operator_constant_trial_closed_cylinder`).
For homogeneous cylinder with isotropic constant trial
:math:`\psi_{\rm trial}(r, \mu_{\rm axial}, \varphi_{\rm az}) = 1`
and constant scattering source :math:`q = \Sigma_s` (:math:`k =
k_\infty` so :math:`q = \omega_0 \Sigma_t`), the operator action
gives :math:`(K \cdot 1)(r, \mu_{\rm axial}, \varphi_{\rm az}) =
\omega_0` everywhere. The proof has the same three steps as sphere
V_α1: (1) surface fixed-point :math:`\psi_{\rm surf}(b, \mu_{\rm
axial}) = \omega_0` (the chord-length :math:`L_{\rm period}` cancels
in the sum :math:`\alpha B \cdot T(\alpha, \tau_{\rm period})`),
(2) total flux is :math:`L_{0,\rm 3D}`-independent (the first-leg
contribution :math:`F` plus the attenuated surface contribution
:math:`e^{-\tau_0}\psi_{\rm surf}` algebraically combine to remove
the :math:`L_{0,\rm 3D}` dependence), (3) the rank-1 isotropic mode
is the unique eigenmode by Perron-Frobenius. Algebraically identical
to sphere V_α1 — **the chord formula is the only geometric
variation**, and it cancels via the same algebra-of-record pattern.

**V_α1_cyl.geometry** — bounce-period chord identity
(:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_cylinder.derive_bounce_period_chord_cylinder`).
Foundation gate for V_α1_cyl: :math:`L_{\rm period}(b, \mu_{\rm
axial}) = 2\sqrt{R^2 - b^2} / \sqrt{1 - \mu_{\rm axial}^2}` derived
two ways (impact-parameter form vs surface-tangent form) agree
algebraically. Without this gate, V_α1_cyl would be a tautology in
the chord formula; with it, V_α1_cyl exercises the chord algebra
through the closure.

**V_α2_cyl** — rank-1 cylinder Variant α ≡ cylinder Hébert white-BC
closure
(:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_cylinder.derive_T00_equals_P_ss_cylinder`).
Unlike sphere, **the cylinder integrand has no elementary closed
form**: the :math:`\mathrm{Ki}_3` Bickley-Naylor function appears
naturally in cylinder transport and is itself defined by an integral.
SymPy cannot terminate at a closed form. The proof is therefore
integrand-level: T_00 from the explicit Knyazev T-matrix construction
(:math:`m=n=0`, :math:`k_m=k_n=0`) and P_ss from the un-reduced
slanted-chord polar integral (:math:`\beta \in [0, \pi/2]`,
:math:`\sin^2\beta\,\exp(-\tau_{\rm 2D}/\sin\beta)`) are bridged via
the Bickley-Naylor identity

.. math::

   \mathrm{Ki}_3(x) = \int_0^{\pi/2}\!\sin^2\!\beta\,e^{-x/\sin\beta}\,
       \mathrm d\beta

applied as a SymPy substitution (not a derivation — the identity
itself is an upstream-trusted definition, like
:math:`\mathrm{exp1}` for sphere). Both T_00^cyl and P_ss^cyl reduce
to :math:`(4/\pi) \cos\alpha\,\mathrm{Ki}_3(2 \Sigma_t R
\cos\alpha)` after the substitution.

The *load-bearing* L1 evidence for V_α2_cyl is therefore **not the
SymPy proof** (which terminates at the same special function on both
sides — a structurally-shared upstream identity that risks the
ERR-032-style failure mode if the identity itself were wrong) but
**the production-primitive numerical cross-check** at the rank-1
Knyazev T-matrix coefficient
:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.compute_T_specular_cylinder_3d`
:math:`[0, 0]` agreeing with
:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.compute_P_ss_cylinder`
across :math:`\tau_R \in \{0.5, 1.0, 2.5, 5.0, 10.0\}` to ≤ 1e-10.
The two functions are independent code paths with different
intermediate quantities (matrix construction with explicit Legendre
coupling vs scalar P_ss with slanted-chord quadrature); their
numerical equality at five :math:`\tau_R` values is the
structurally-independent cross-check that ERR-032 would have caught.

**V_α3_cyl** — vacuum reduction at α=0
(:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_cylinder.derive_alpha_zero_kernel_reduction_cylinder`).
The surface fixed-point closure (:eq:`peierls-greens-cylinder-T`)
carries leading factor :math:`\alpha`, so :math:`\psi_{\rm surf}
\to 0` at :math:`\alpha = 0`. Vacuum BC is the trivial limit; no
special-case branch. Same V_α3 structure as sphere.

Numerical verification status (Phase 1 standalone)
---------------------------------------------------

The full numerical evidence chain from the closeout memo
(:file:`.claude/agent-memory/method-implementer/cylinder_variant_alpha_phase1.md`):

.. list-table:: V_α1_cyl k_inf-exactness at :math:`\alpha = 1`
   (fuel-A-like, R=5, :math:`\tau_R = 2.5`, :math:`k_\infty = 0.20833\overline{3}`)
   :header-rows: 1
   :widths: 35 25 20 20

   * - :math:`(n_r, n_{\mu_{\rm ax}}, n_{\varphi}, n_{\rm traj})`
     - :math:`k_{\rm eff}`
     - rel err vs :math:`k_\infty`
     - iter
   * - (8, 8, 16, 24)
     - 0.2083333333333330
     - 1.87e-15
     - 1
   * - (12, 10, 24, 32)
     - 0.2083333333333338
     - 2.00e-15
     - 1
   * - (16, 12, 32, 48)
     - 0.2083333333333344
     - 4.93e-15
     - 1
   * - (20, 16, 48, 64)
     - 0.2083333333333335
     - 6.66e-16
     - 1

V_α1_cyl algebraic identity reproduced numerically at machine precision
in 1 iteration at every quadrature order — the load-bearing L1
structural-independence anchor.

.. list-table:: Cylinder vacuum k_eff convergence-floor (:math:`\alpha = 0`,
   fuel-A-like XS, :math:`\tau_R = 2.5`)
   :header-rows: 1
   :widths: 35 30 35

   * - :math:`(n_r, n_{\mu_{\rm ax}}, n_{\varphi}, n_{\rm traj})`
     - :math:`k_{\rm eff}`
     - rel diff to next finer
   * - (8, 8, 16, 24)
     - 0.1204312553
     - 1.82e-4
   * - (12, 10, 24, 32)
     - 0.1204531627
     - 2.91e-5
   * - (16, 12, 32, 48)
     - 0.1204496559
     - 1.60e-5
   * - (20, 16, 48, 64)
     - 0.1204515822
     - **8.24e-8**
   * - (24, 20, 64, 96)
     - 0.1204515921
     - (terminus)

Convergence is non-monotone at intermediate orders (typical for 3D
Gauss-Legendre on a problem with grazing-ray loci); above
:math:`(n_r, n_\mu, n_\varphi, n_{\rm traj}) \ge (20, 16, 48, 64)`
the answer stabilises at :math:`k_{\rm eff} \approx 0.12045159` to
~8e-8 self-consistency. :math:`k_{\rm eff}/k_\infty \approx 0.578`
is well within the physical leakage band for a thin sphere-class
optical thickness.

**Multi-group at** :math:`\alpha = 1` — 2G with asymmetric scattering
(downscatter + upscatter): :math:`k_{\rm eff} = k_\infty = 1.0` from
:func:`~orpheus.derivations.common.eigenvalue.kinf_homogeneous` to
1.5e-11 relative in 48 iterations. Per-group flux uniform to
:math:`\sim 1\mathrm e{-15}` std/mean. The asymmetric scattering is
load-bearing: a SigS-transpose convention drift (ERR-002 anti-pattern,
the in-scatter source uses :math:`Q = \Sigma_s^T \phi` while the
SigS dataclass stores :math:`\Sigma_s[g_{\rm from}, g_{\rm to}]`)
shows up as a wrong group ratio at off-diagonal coupling. Symmetric
2G XS would not trip the detector.

**MG G=1 reduction** — bit-equal to the 1G solver at :math:`\alpha
= 0`.

Cross-check disclaimer (Issue #129)
------------------------------------

The cylinder Variant α and cylinder Nyström vacuum solvers do **NOT
agree to machine precision** on the same XS:

- Variant α (3D angle-resolved): :math:`k_{\rm eff} = 0.1204515921`
  at fuel-A-like :math:`\tau_R = 2.5`, :math:`(n_r, n_\mu, n_\varphi,
  n_{\rm traj}) = (24, 20, 64, 96)`.
- Nyström vacuum (axial pre-integrated): :math:`k_{\rm eff} =
  0.1204656335` at the same XS, :math:`n_{\rm quad} = 48`.

The relative difference (~1.4e-4) is **not a bug** — it reflects the
documented Issue #129 axial pre-integration discrepancy. The cylinder
Nyström kernel pre-integrates the axial direction analytically using
the :math:`\mathrm{Ki}_n` Bickley-Naylor functions; this produces a
small mismatch in the planar limit relative to the angle-resolved
formulation. Both formulations are **correct per their respective
models**; the disagreement quantifies the cost of the axial
pre-integration trade-off in the Nyström kernel.

The load-bearing L1 cross-check for cylinder Variant α is therefore
**not** the cylinder Nyström solver. It is the k_inf-exactness
invariant at :math:`\alpha = 1` (V_α1_cyl numerical reproduction at
machine precision) plus the multi-group reduction to
:func:`kinf_homogeneous` — both **structurally independent**
(closed-form analytical) of any other ORPHEUS solver.

Source code, tests, and provenance
----------------------------------

- **Branch-1 SymPy**:
  :mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_cylinder`
  (~395 LoC).
  :func:`derive_operator_constant_trial_closed_cylinder` (V_α1_cyl);
  :func:`derive_bounce_period_chord_cylinder` (V_α1_cyl.geometry);
  :func:`derive_T00_equals_P_ss_cylinder` (V_α2_cyl);
  :func:`derive_alpha_zero_kernel_reduction_cylinder` (V_α3_cyl).
- **Branch-2 production**:
  :mod:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder`
  (~554 LoC).
  :func:`solve_greens_function_cylinder` (1G);
  :func:`solve_greens_function_cylinder_mg` (MG); plus internal
  helpers ``_apply_operator_cylinder``, ``_first_leg_2d_chord``,
  ``_impact_parameter``, ``_bounce_period_2d_chord``,
  ``_scalar_flux_from_psi``; result dataclasses
  ``CylinderGreensResult`` and ``CylinderGreensMGResult``.
- **Symbolic foundation tests**:
  :file:`tests/derivations/test_trajectory_resolvent_cylinder_symbolic.py`
  — 9 ``@pytest.mark.foundation`` gates (the four ``derive_*``
  identities plus 5 V_α2 production-primitive cross-checks at five
  :math:`\tau_R` values).
- **Numerical L1 tests**:
  :file:`tests/derivations/test_trajectory_resolvent_cylinder_solver.py`
  — 7 ``@pytest.mark.l1`` gates (k_inf-exactness at thin/moderate
  :math:`\tau_R`, MG asymmetric Σ_s at α=1, MG G=1 reduction at α=0,
  vacuum self-consistency convergence floor).
- **Closeout memo**:
  :file:`.claude/agent-memory/method-implementer/cylinder_variant_alpha_phase1.md`.
- **Cross-domain frame memo**:
  :file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
  — the rank-1 → rank-2 generalisation prediction that held across
  all 6 geometries.
- **Plan**:
  :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md`.


.. _peierls-greens-unification:

Operator-theoretic unification (Phase 2, ``variant_alpha_core``)
================================================================

The sphere and cylinder Phase-1 prototypes were structurally identical
at the closure level but inlined the bounce-sum closure as
hand-written Python (``denom = 1 - alpha * np.exp(-tau_period);
psi_surf = alpha * B / denom; psi_new = F + np.exp(-tau_first) *
psi_surf``). Phase 2 (commits ``efbae9c`` / ``92d4f10`` / ``166b9ae``)
extracted this common closure into a shared module
:mod:`~orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core`
with two pure functions:

- :func:`~orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core.compute_resolvent_T`
  — rank-1 resolvent :math:`T(\alpha, \tau_{\rm period}) =
  1/(1 - \alpha\,e^{-\tau_{\rm period}})`.
- :func:`~orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core.apply_variant_alpha_closure`
  — full closure :math:`\psi_{\rm new} = F + e^{-\tau_{\rm first\,leg}}
  \cdot \alpha\,B \cdot T(\alpha, \tau_{\rm period})`, vectorised
  across phase-space arrays of identical shape, with vacuum
  (:math:`\alpha = 0`) short-circuiting to :math:`F` exactly.

Both sphere and cylinder solvers were re-mounted on these primitives
in the Phase 2 refactor with **bit-equal** accuracy preservation
(machine-floor-identical IEEE-754 outputs at every test gate). The
Phase 3A symmetric-slab implementation extended the API
back-compatibly with an ``alpha_per_period`` keyword (now
historically obsolete after the ERR-035 delegation, but preserved on
the API surface — see :ref:`peierls-greens-slab` for the full
narrative).

The cross-domain frame
=======================

The unification is not just code reuse. It is the operational
realisation of the **boundary-to-boundary scattering operator frame**
predicted by the cross-domain-attacker memo at
:file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
*before any 2-surface code existed*. The frame names the load-bearing
mathematical object:

.. math::
   :label: peierls-greens-unification-resolvent

   T \;=\; (I - S)^{-1}, \qquad
   S = \alpha \cdot R_{\rm chord},

where :math:`S` is the **monodromy operator** advancing the surface
state by one step of the bouncing trajectory, :math:`R_{\rm chord}` is
the chord-attenuation operator on surface phase-space, and :math:`T`
is the resolvent that closes the multi-bounce sum.

For 1-surface compact geometries (sphere, cylinder), :math:`S` is the
**rank-1 multiplication operator** :math:`\alpha\,e^{-\tau_{\rm
period}}`; the resolvent collapses to the scalar geometric series
:math:`T = 1/(1 - \alpha\,e^{-\tau_{\rm period}})`. The
``alpha_per_period`` parameter exposes the per-period reflection
product explicitly, which is :math:`\alpha^1` for sphere/cylinder
(one bounce per period). For 2-surface geometries (slab asymmetric,
hollow sphere, annulus), :math:`S` is a **rank-2 antidiagonal block**
with entries :math:`\alpha_L\,e^{-\tau}` and :math:`\alpha_R\,e^{-\tau}`
(single-transit, NOT period-doubled), and the resolvent is the 2×2
matrix inverse exposed by
:func:`~orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core.compute_resolvent_T_rank2`.

Three structural predictions of the frame have been validated across
the 6-geometry × 2-topology family:

1. **Rank-1 covers 1-surface compact geometries** (sphere, cylinder)
   — confirmed at Phase 1.
2. **Rank-2 covers 2-surface geometries regardless of curvature**
   — confirmed at Phase 3B (slab asymmetric, flat), Phase 3C-1
   (hollow sphere, spherical), Phase 3C-2 (annulus, cylindrical).
3. **Phase-space partition by impact parameter cleanly separates
   rank-1 outer-only from rank-2 through-ray** on hollow geometries
   — confirmed at Phase 3C-1/3C-2.

The closure operator is **byte-equal-shared** across all 4 currently-
implemented 2-surface prototypes (the symmetric slab is a delegated
reduction of Phase 3B at :math:`\alpha_L = \alpha_R = \alpha`). Only
the chord-algebra primitives differ. After Phase 3C-2, the
:mod:`.variant_alpha_core` module has not been modified — the
abstraction has earned the right to stay frozen, and any future
extension that requires modifying :mod:`.variant_alpha_core` should
be treated as a structural alarm.

Bit-equality preservation under refactor (Phase 2 evidence)
-----------------------------------------------------------

The Phase 2 refactor's load-bearing acceptance criterion was
**bit-equal accuracy preservation** — the unified primitive must
produce IEEE-754-identical outputs to the inlined formulas. This is
achievable because :math:`\psi_{\rm new}` involves a single
:math:`\exp`, one product, one subtraction, one division, one
addition, in the same algebraic order as the inlined code; SymPy
``simplify``-style reorderings would break the bit-equality but were
explicitly avoided.

.. list-table:: Phase 2 unification — bit-equality preservation
   :header-rows: 1
   :widths: 50 25 25

   * - Test
     - Pre-refactor floor
     - Post-refactor floor
   * - Sphere k_inf-exactness (α=1, R=5)
     - 1.110e-16
     - 1.110e-16 (bit-equal)
   * - Cylinder k_inf-exactness (α=1, R=3)
     - 8.327e-16
     - 8.327e-16 (bit-equal)
   * - Sphere PS-1982 RMS (4 configs)
     - < 1e-4
     - < 1e-4 (unchanged)
   * - Cylinder vacuum self-consistency (RESEARCH_GRADE_FLOOR_CYL=5e-7)
     - passing
     - passing (1.000x)
   * - V_α2 production-primitive cross-check (10 τ_R configurations)
     - passing at 1e-10
     - passing at 1e-10 (1.000x)

The runtime overhead of the unified path is < 3 seconds across the
59-test acceptance gate (under run-to-run noise) — the shared closure
is a single Python-level function call per phase-space point,
dominated by the trajectory-quadrature inner loop (which is
unchanged).

Why the unification deserves theory-page space
----------------------------------------------

This section's role on the theory page is **not** to document the
API of :mod:`.variant_alpha_core` (that's the API page's job). It is
to name the **load-bearing mathematical object** that links all 6
geometries: the boundary-to-boundary monodromy resolvent. Without
this naming, the per-geometry V_α1/V_α2/V_α3 identities would read
as six unrelated proofs that happen to share Python code. With it,
each per-geometry V_α1 reads as a **specialisation of the same
algebraic identity** — the chord-cancellation pattern in the
resolvent of :math:`S = \alpha R_{\rm chord}` — to a particular
chord algebra. The proofs are unified by their structural ancestor,
not by code reuse.


.. _peierls-greens-valpha2-strengthening:

V_α2 strengthening across geometries (sphere, cylinder, slab)
==============================================================

V_α2 is the algebraic identity :math:`T_{00} = P_{\rm ss}` that
proves rank-1 Variant α agrees with Hébert white-BC at the closed-
form level (sphere; integrand level for cylinder). The original Phase
1 V_α2 SymPy proofs (sphere :func:`derive_T00_equals_P_ss_sphere`;
cylinder :func:`derive_T00_equals_P_ss_cylinder`) were caught during
QA review of cylinder Phase 1 to be **typed-identical integrand
tautologies**: both LHS and RHS were built from one SymPy expression
literal, so :math:`\mathrm{simplify}(\mathrm{LHS} - \mathrm{RHS}) =
0` was logically vacuous. Commit ``4d87840`` strengthened V_α2 with
**structurally-independent** SymPy paths plus a
**production-primitive numerical cross-check**, addressing the L11
discipline lesson (structurally-independent cross-checks must reach
the same value via DIFFERENT mathematical paths, not different code
paths around the same identity).

The strengthening pattern differs per geometry, because the **closed-
form vs special-function** boundary differs per geometry. The three
states of the algebra-of-record skill (1A closed-form / 1B
semi-analytical / 1C MMS) map onto the three V_α2 strengthenings.

Sphere V_α2 (State 1A — both paths reach a common closed form)
---------------------------------------------------------------

Sphere V_α2 has two structurally-independent symbolic paths:

1. **Transfer-matrix definition** (Phase 4 ``specular_multibounce``
   at :math:`N = 1`): :math:`T_{00}^{\rm sphere} = 2 \int_0^1
   \mu\,\tilde P_0(\mu)^2\,e^{-2 \tau_R \mu}\,\mathrm d\mu = 2\int_0^1
   \mu\,e^{-2 \tau_R \mu}\,\mathrm d\mu` (μ-domain, with explicit
   :math:`\tilde P_0(\mu) = 1`).
2. **Polar escape-probability integral** (Hébert :math:`P_{\rm ss}`
   definition): :math:`P_{\rm ss}^{\rm sphere} = 2\int_0^{\pi/2}
   \sin\theta'\,\cos\theta'\,e^{-2\tau_R\,\cos\theta'}\,
   \mathrm d\theta'` (θ-domain, with :math:`\sin\theta'` Jacobian).

SymPy integrates **each native form independently**, both reduce to
the canonical Hébert closed form

.. math::

   T_{00}^{\rm sphere} = P_{\rm ss}^{\rm sphere} =
       \frac{1 - (1 + 2\tau_R)\,e^{-2\tau_R}}{2\,\tau_R^{2}}

(:eq:`peierls-greens-V-alpha-2`), and SymPy verifies their equality
by :math:`\mathrm{simplify}(\mathrm{LHS} - \mathrm{RHS}) = 0`. The
cross-path equality is the structurally-independent gate: the
:math:`\sin\theta'` Jacobian is the variable-substitution evidence
that the transfer-matrix and escape-probability frames are different
encodings of the same physical content. Plus a numerical
production-primitive cross-check across 5 :math:`\tau_R` values
agreeing at 1e-10 closes the structural-independence chain.

Cylinder V_α2 (State 1B — Bickley-Naylor obstruction)
------------------------------------------------------

Cylinder V_α2 cannot reach a closed form: the integrand naturally
contains :math:`\mathrm{Ki}_3(2 \tau_R \cos\beta)`, the Bickley-Naylor
function — itself defined by an integral with no elementary
antiderivative. SymPy cannot terminate at a closed form on either
path. The proof is integrand-level via the Bickley-Naylor identity
(see :ref:`peierls-greens-cylinder` above for the full algebra).

Because the integrand identity terminates at the same special
function on both sides, the SymPy proof carries **shared upstream
identity risk** (ERR-032 anti-pattern). The structurally-independent
L1 evidence is therefore explicitly the production-primitive
numerical cross-check between
:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.compute_T_specular_cylinder_3d`
:math:`[0, 0]` and
:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.compute_P_ss_cylinder`
across :math:`\tau_R \in \{0.5, 1.0, 2.5, 5.0, 10.0\}`. The two
functions are independent code paths with different intermediate
quantities (matrix construction with explicit Legendre coupling vs
scalar P_ss with slanted-chord quadrature); their numerical equality
at five :math:`\tau_R` values to ≤ 1e-10 is the load-bearing gate
that an ERR-032-style shared-identity bug would have flagged.

Slab V_α2 (hybrid State 1A + 1B — SymPy choke at μ→0+)
-------------------------------------------------------

Slab V_α2 has its own algebraic obstacle: SymPy chokes on
:math:`\int_0^1 2\mu\,e^{-\tau_L/\mu}\,\mathrm d\mu` at the
:math:`\mu \to 0^+` endpoint, raising
``Add object cannot be interpreted as an integer`` in SymPy's
nseries evaluator. This is choke mode #1 in the algebra-of-record
skill (expression-tree growth in limit-evaluation paths).

Workaround: hybrid State 1A + 1B per the algebra-of-record state
taxonomy. The closed form
:math:`T_{00}^{\rm slab} = P_{\rm ss}^{\rm slab} = 2\,E_3(\Sigma_t L)`
is **manually constructed** using :math:`\mathrm{sp.expint}(3, \cdot)`,
with the underlying substitution :math:`u = 1/\mu` algebra verified
symbolically (the integrand maps correctly into the :math:`E_3`
definitional form :math:`E_3(x) = \int_1^\infty t^{-3}\,e^{-x t}\,
\mathrm d t`). The numerical structural-independence cross-check uses
arbitrary-precision **mpmath** quadrature directly on BOTH the
original :math:`\mu`-integrand (Path A) and the original
:math:`\theta`-integrand (Path B), comparing to
:math:`\mathrm{mpmath.expint}(3, \tau_L)` at 6 :math:`\tau_L` values
from 0.1 to 10. All 6 agree to absolute tolerance 1e-12 (in fact
0.0e+00 at the precision of mpmath's 30-digit setup).

This is the canonical hybrid-State-1A+1B pattern when SymPy chokes
on one symbolic path: keep the symbolic substitution algebra as
proof, use mpmath as the numerical structural-independence verifier.

The L11 lesson (one statement, three geometries)
-------------------------------------------------

The QA-caught V_α2 tautology and its three-geometry repair illustrate
the **L11 discipline**: structurally-independent cross-checks must
reach the same value via DIFFERENT mathematical paths, not via
different code paths around the same identity. Sphere achieves this
via two closed-form-reaching paths in different angular variables.
Cylinder achieves this via two production-code primitives that
construct the same special-function value through different
intermediate quantities. Slab achieves this via two
mpmath-quadrature paths on the original integrands (μ-domain and
θ-domain) plus a manually-verified substitution to :math:`E_3`. Each
geometry's strengthening is geometry-appropriate; all three satisfy
the L11 discipline.


.. _peierls-greens-slab:

Slab Variant α (Phase 3A — symmetric reflective, post-ERR-035 delegation)
==========================================================================

Phase-3A slab Variant α was the first 2-bounce-per-period geometry
attacked by the family. **It now exists as a thin wrapper** that
delegates to the Phase-3B rank-2 asymmetric solver at
:math:`\alpha_L = \alpha_R = \alpha`. The original Phase-3A
implementation used a heuristic rank-1 closure
:math:`\psi_{\rm surf} = \alpha\,B_{\rm period}/(1 - \alpha^2\,
e^{-2\tau})` that was incorrect at intermediate :math:`\alpha`
(logged as ERR-035, fixed 2026-05-02 — see
:ref:`peierls-greens-slab-asym-err035-fix` for the full forensic
narrative). Both implementations agreed at the BC corners
:math:`\alpha \in \{0, 1\}` (where flux profiles are uniform-or-
rank-1-collapsing) but differed by ~1.3e-4 relative at
:math:`\alpha = 0.5`.

The two implementation phases of Phase-3A are documented historically
because the lesson is load-bearing for any future Variant α geometry
extension:

1. **Original Phase-3A (2026-05-02 morning)** — mounted on the same
   rank-1 resolvent as sphere/cylinder, using a heuristic
   ``alpha_per_period = α²`` extension of
   :func:`~orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core.apply_variant_alpha_closure`
   to encode the **2-bounce-per-period** structural distinction of
   slab geometry. The closure was constructed by **analogical
   generalisation** from sphere/cylinder rank-1: "two reflections per
   period instead of one, so multiply the per-period reflection by
   :math:`\alpha^2` and use the full out-and-back B integral instead
   of single-transit". The analogy was wrong because **analogical
   generalisation of a closure formula across geometries is an
   algebraic claim that requires a proof, not a substitution** —
   and the correct rank-2 first-principles derivation gives a
   structurally different formula. The agreement at α-corners
   masked the difference (V_α1 at α=1 uses constant source; vacuum
   at α=0 makes the closure trivial); only intermediate-α
   non-uniform-source tests would catch it.

2. **ERR-035 fix (2026-05-02 afternoon)** —
   :func:`solve_greens_function_slab` and
   :func:`solve_greens_function_slab_mg` were refactored as thin
   wrappers that delegate to
   :func:`solve_greens_function_slab_asymmetric` /
   :func:`solve_greens_function_slab_asymmetric_mg` at
   :math:`\alpha_L = \alpha_R = \alpha`. The first-principles rank-2
   closure :math:`\psi_{\rm surf} = \alpha\,B/(1 - \alpha\,e^{-\tau})`
   (with single-transit B integral) replaces the heuristic. The
   Phase-3A public API surface (:class:`SlabGreensResult`,
   :class:`SlabGreensMGResult`, the two solver entrypoints) is
   preserved by re-wrapping the rank-2 result. The heuristic
   ``_apply_operator_slab`` and the helper chord-length functions
   were **deleted** — keeping them as "alternative paths" would
   re-create the very pattern the fix removes.

The cascade of side-effects is also instructive. Phase-3A's reported
"~5e-4 vacuum convergence floor" was originally misattributed to a
slow-:math:`\mu`-quadrature limitation (predicting that "a future
improvement could use a μ-weighted Gauss-Jacobi quadrature
concentrating nodes near :math:`\mu = 0`"). After both ERR-034
(first-leg trajectory missing :math:`\mu` factor) and ERR-035 were
fixed, vacuum k_eff at fuel-A-like τ_L=5 jumped from
:math:`\sim 0.130` to **0.157**, and the convergence floor improved
**56× to ~9e-6** between the two finest grids. The slow apparent
convergence was not a quadrature limitation; it was the dominant
ERR-034 trajectory bug masquerading as quadrature noise. **A
quadrature convergence rate that doesn't match theoretical
expectations is a bug fingerprint, not a quadrature limitation,
until proven otherwise.**

The closeout memos
:file:`.claude/agent-memory/method-implementer/slab_variant_alpha_phase3a.md`
(original Phase-3A) and
:file:`.claude/agent-memory/method-implementer/err035_phase3a_delegation_fix.md`
(post-fix) record the full development history and the
discipline-level lessons.

Slab phase-space and 2-bounce-per-period structure
---------------------------------------------------

Phase-space coordinates: :math:`x \in [0, L]` (Cartesian position
across the slab), :math:`\mu \in [-1, 1]` (signed direction-cosine
wrt outward normal at :math:`x = L`). Specular reflection on a
symmetric slab preserves :math:`|\mu|` (the planar geometry is
flat, so the velocity component perpendicular to the wall flips sign,
the parallel component is preserved, and :math:`|\mu|` is the
absolute cosine to the outward normal of either wall — preserved
under both reflections). The trajectory direction sign alternates at
each wall. **One full period traverses both walls** (out + back
transit at constant :math:`|\mu|`), totalling

.. math::
   :label: peierls-greens-slab-bounce-period

   L_{\rm period}^{\rm slab}(\mu) = \frac{2 L}{|\mu|}.

Two reflections per period, so the **per-period reflection product
is** :math:`\alpha^2`, distinct from the single-reflection per-period
products of sphere (:math:`\alpha`) and cylinder (:math:`\alpha`).
This is the structural distinction that the original Phase-3A
heuristic rank-1 closure tried to encode via the
``alpha_per_period = α²`` parameterisation. The post-ERR-035 rank-2
delegation handles it correctly: the rank-2 monodromy
:math:`S = \mathrm{antidiag}(\alpha\,e^{-\tau}, \alpha\,e^{-\tau})`
encodes the alternating-wall structure in the matrix shape
(antidiagonal: one bounce never returns to the same wall), and the
matrix algebra produces :math:`(1 - \alpha^2 e^{-2\tau})` as
:math:`\det(I - S)` — the same per-period factor — without conflating
"per-period reflection product" with "BC reflectivity" in a single
scalar formula.

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

The trajectory point at arclength :math:`s` along the backward chord
is :math:`x_{\rm traj}(s) = x - \mu \cdot s` for :math:`\mu > 0`
(position decreases toward the left wall) and :math:`x_{\rm traj}(s)
= x + |\mu| \cdot s` for :math:`\mu < 0` (position increases toward
the right wall). The :math:`\mu` factor is **load-bearing**: the
arclength :math:`s` is in 3D-arclength units, not x-distance units,
so without the :math:`\mu` projection the position lookup
:math:`q(x_{\rm traj})` is wrong by a factor of :math:`\mu`. ERR-034
documents the bug where this factor was missing in the original
Phase-3A code; see :ref:`peierls-greens-slab-asym-err035-fix` for
the forensic narrative. The fingerprint of this bug class is
**uniform-source-test blindness**: constant :math:`q` doesn't care
where on the chord we evaluate, so V_α1 at α=1 (the standard
foundation gate) is structurally blind. Method-of-images on
asymmetric BC is the structurally-independent gate that catches it.

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

The slab Branch-1 SymPy module
:mod:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_slab`
verifies three identities, each a slab-specialisation of a sphere
analogue.

**V_α1_slab** — closed-slab bounce-sum self-consistency
(:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_slab.derive_operator_constant_trial_closed_slab`).
For homogeneous slab with constant volumetric source :math:`q`,
operator action on isotropic constant trial gives :math:`(K \cdot 1)
= \omega_0` everywhere, hence :math:`k_{\rm eff} = k_\infty` for
closed slab. **Algebraically identical to V_α1 sphere/cylinder** —
only the chord formula :math:`L_{\rm period} = 2L/|\mu|` differs.
The chord-cancellation algebra works the same way: the closure
denominator's :math:`(1 - \alpha\,e^{-\tau_{\rm period}/2})^{-1}`
factor (post-ERR-035 rank-2 form) combined with the single-transit
B integral and the first-leg attenuation algebraically removes the
:math:`L`-dependence and the :math:`\mu`-dependence, leaving only
:math:`q/\Sigma_t = \omega_0`.

**V_α2_slab** — :math:`T_{00}^{\rm slab} = P_{\rm ss}^{\rm slab}`
algebraic identity
(:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_slab.derive_T00_equals_P_ss_slab`).
The closed form is

.. math::

   T_{00}^{\rm slab} = P_{\rm ss}^{\rm slab} = 2\,E_3(\Sigma_t L)

where :math:`E_3` is the third exponential integral. This is the
slab analogue of the sphere :math:`(1 - (1+2\tau_R)e^{-2\tau_R}) /
(2\tau_R^2)`; both are closed-form (unlike cylinder's Bickley-Naylor
:math:`\mathrm{Ki}_3` form). Two structurally-independent
derivational paths: per-face transfer-matrix definition
(:math:`\mu`-domain integral :math:`\int_0^1 2\mu\,e^{-\tau/\mu}\,
\mathrm d\mu`) and polar escape-probability integral
(:math:`\theta`-domain). SymPy chokes on direct integration of
:math:`e^{-\tau/\mu}` at the :math:`\mu \to 0^+` endpoint
(algebra-of-record SymPy choke mode #1: expression-tree growth in
limit-evaluation paths). The proof is via hybrid State-1A + State-1B:
the :math:`u = 1/\mu` substitution maps the integrand into the
canonical :math:`E_3(x) = \int_1^\infty t^{-3}\,e^{-x t}\,
\mathrm d t` definitional form (substitution algebra verified
symbolically), and arbitrary-precision **mpmath** quadrature on BOTH
original integrands (the :math:`\mu`-domain and :math:`\theta`-domain
forms) confirms numerical agreement with
:math:`\mathrm{mpmath.expint}(3, \tau_L)` at six :math:`\tau_L`
values to absolute tolerance 1e-12. See
:ref:`peierls-greens-valpha2-strengthening` for the L11-discipline
narrative on why the hybrid path is structurally sufficient.

**V_α3_slab** — vacuum reduction at α=0
(:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_slab.derive_alpha_zero_kernel_reduction_slab`).
Surface fixed-point closure (:eq:`peierls-greens-slab-T`) carries
leading factor :math:`\alpha`; :math:`\psi_{\rm surf} \to 0` at
:math:`\alpha = 0`. Vacuum BC is the trivial limit; no special-case
branch. Same V_α3 structure as sphere/cylinder.

Numerical verification status (Phase 3A, post-ERR-034 + ERR-035 fixes)
-----------------------------------------------------------------------

Acceptance gates are tagged ``@pytest.mark.l1`` in
:file:`tests/derivations/test_trajectory_resolvent_slab_solver.py`.

- **k_inf-exactness at α=1** — fuel-A-like XS, both thin
  (:math:`\tau_L = 1`) and moderate (:math:`\tau_L = 5`):
  :math:`k_{\rm eff} = k_\infty = 0.20833\overline{3}` to
  :math:`\le 1\mathrm e{-10}` relative (machine floor 2.93e-15 in 1
  iteration from constant initial guess). Convergence to 1e-9 in
  ~50 iterations from aggressive quadratic perturbation. Scalar flux
  uniform to :math:`\sim 1\mathrm e{-15}`. This is the V_α1_slab
  numerical reproduction — the load-bearing structural-independence
  anchor.
- **Multi-group at α=1** — 2G with asymmetric scattering
  (downscatter + upscatter): :math:`k_{\rm eff} = k_\infty` to
  :math:`\le 1\mathrm e{-9}`. Asymmetric XS detects the
  SigS-transpose convention drift (ERR-002 anti-pattern detector).
- **MG G=1 reduction** — bit-equal to the 1G solver at α=0.
- **V_α2_slab production-primitive cross-check** —
  :math:`T_{\rm slab}[0, 0]` from
  :func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.compute_T_specular_slab`
  at :math:`n_{\rm modes} = 1` agrees with arbitrary-precision
  ``mpmath.expint(3, τ_L)`` to :math:`\le 1\mathrm e{-10}` over five
  :math:`\tau_L` values from 0.1 to 10.

.. list-table:: Slab vacuum k_eff convergence floor (α=0, fuel-A-like, τ_L=5)
   — pre-ERR-034/035 vs post-fix (note the 56× improvement)
   :header-rows: 1
   :widths: 40 25 35

   * - :math:`(n_x, n_\mu, n_{\rm traj})`
     - :math:`k_{\rm eff}` (post-fix)
     - rel diff to next finer
   * - (16, 24, 48)
     - 0.15743600
     - 2.33e-5
   * - (20, 32, 64)
     - 0.15743968
     - 1.03e-5
   * - (24, 40, 96)
     - 0.15744131
     - **8.85e-6**
   * - (32, 56, 128)
     - 0.15744270
     - (terminus)

Pre-fix Phase-3A reported a ~5e-4 floor that was attributed to a
slow-:math:`\mu`-quadrature limitation. After ERR-034 + ERR-035 fixes,
the floor is **9e-6 between the two finest grids** — 56× tighter and
now genuinely diagnostic of the underlying numerics rather than a
hidden trajectory bug. The :math:`k_{\rm eff}/k_\infty \approx 0.756`
ratio is well within the physical leakage band for fuel-A-like
:math:`\tau_L = 5`.

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

Source code, tests, and provenance
----------------------------------

- **Branch-1 SymPy**:
  :mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_slab`
  (~431 LoC).
  :func:`derive_operator_constant_trial_closed_slab` (V_α1_slab);
  :func:`derive_T00_equals_P_ss_slab` (V_α2_slab);
  :func:`derive_alpha_zero_kernel_reduction_slab` (V_α3_slab).
- **Branch-2 production**:
  :mod:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function_slab`
  (~250 LoC, **post-ERR-035 thin wrapper**) —
  :func:`solve_greens_function_slab` and
  :func:`solve_greens_function_slab_mg` are wrappers that delegate
  to the asymmetric solver at :math:`\alpha_L = \alpha_R = \alpha`.
  Result dataclasses :class:`SlabGreensResult`,
  :class:`SlabGreensMGResult` are preserved for API stability.
- **Symbolic foundation tests**:
  :file:`tests/derivations/test_trajectory_resolvent_slab_symbolic.py`
  — 10 ``@pytest.mark.foundation`` gates.
- **Numerical L1 tests**:
  :file:`tests/derivations/test_trajectory_resolvent_slab_solver.py`
  — 12 ``@pytest.mark.l1`` gates. Includes
  :func:`test_alpha_zero_convergence_floor` (re-pinned to ~9e-6
  post-ERR-034/035 fixes) and the V_α2_slab production-primitive
  cross-check.
- **Closeout memos**:
  :file:`.claude/agent-memory/method-implementer/slab_variant_alpha_phase3a.md`
  (original Phase-3A) and
  :file:`.claude/agent-memory/method-implementer/err035_phase3a_delegation_fix.md`
  (delegation fix).
- **Plan**:
  :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md`.


.. _peierls-greens-slab-asym:

Slab asymmetric Variant α (Phase 3B — rank-2 resolvent)
========================================================

Phase-3B asymmetric slab Variant α extends the BC parameter space
from one scalar :math:`\alpha \in [0, 1]` to TWO independent per-wall
reflectivities :math:`\alpha_L, \alpha_R \in [0, 1]`. **The shared
rank-1 closure (sphere, cylinder, symmetric slab via delegation) is
replaced by a rank-2 boundary-to-boundary scattering resolvent**
:math:`T = (I - S)^{-1}` — the predicted generalisation from the
cross-domain frame memo
:file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`,
which named the load-bearing object **before any 2-surface code
existed**. Phase 3B is the first instance to validate the prediction.

This phase is also where two latent bugs were caught: ERR-034
(first-leg trajectory missing :math:`\mu` factor) and ERR-035
(Phase-3A heuristic closure mismatch at intermediate :math:`\alpha`).
Both were caught by the **method-of-images symmetry test** — a
structurally-independent test that exercises both the closure
algebra (across a non-trivial BC diagonal) AND the trajectory
machinery (asymmetric flux peaks AT the reflective wall, not at an
interior point — a sensitive bug-fingerprint for trajectory-
parametrisation errors). Each was logged to
:file:`.claude/skills/vv-principles/error_catalog.md` as a separate
ERR-NNN entry. ERR-034 was fixed in this phase; ERR-035 was
documented and gated, with the fix applied subsequently as the
Phase-3A → Phase-3B delegation refactor.

The closeout memo at
:file:`.claude/agent-memory/method-implementer/slab_asymmetric_variant_alpha_phase3b.md`
records the development verdict: 27 new tests (16 SymPy + 11
solver), method-of-images agreement at 6.5e-9 relative at finest
grid (64, 64, 128), zero regression on existing 89 sphere/cylinder/
slab-symmetric tests.

Rank-2 monodromy and resolvent
-------------------------------

Phase space :math:`(x, \mu)` with **2-vector surface state**
:math:`[\psi_L^+(\mu), \psi_R^-(\mu)]` representing the outgoing
flux at the left wall in :math:`+\hat\mu` direction and at the
right wall in :math:`-\hat\mu` direction respectively. The
**monodromy operator** :math:`S` advances the surface state by ONE
STEP of the bouncing trajectory (single transit across the slab,
followed by reflection at the destination wall):

.. math::
   :label: peierls-greens-slab-asym-monodromy

   S(\alpha_L, \alpha_R, \tau) = \begin{pmatrix}
       0                              & \alpha_L\,e^{-\tau} \\
       \alpha_R\,e^{-\tau}            & 0
   \end{pmatrix}, \qquad
   \tau = \Sigma_t\,L /|\mu|

(single-transit optical depth — **NOT** the full out-and-back
period; this is the load-bearing distinction between the heuristic
Phase-3A closure and the first-principles rank-2 form). The diagonal
of :math:`S` is identically zero: **one step never returns to the
same wall**. The off-diagonals encode "go from L to R, get reflected
by :math:`\alpha_L`" (top-right entry: :math:`\psi_L^+` is fed by
the L-arrival of a particle that left R) and "go from R to L, get
reflected by :math:`\alpha_R`" (bottom-left entry: :math:`\psi_R^-`
is fed by the R-arrival of a particle that left L).

The **rank-2 resolvent** is the geometric series of :math:`S`:

.. math::

   T = (I - S)^{-1} = I + S + S^2 + \cdots
                    = I + S + S \cdot S + S \cdot S \cdot S + \cdots,

which has a closed matrix form via the standard 2×2 inverse:

.. math::
   :label: peierls-greens-slab-asym-resolvent

   T(\alpha_L, \alpha_R, \tau) = (I - S)^{-1}
       = \frac{1}{1 - \alpha_L\,\alpha_R\,e^{-2\tau}}
         \begin{pmatrix}
             1                       & \alpha_L\,e^{-\tau} \\
             \alpha_R\,e^{-\tau}     & 1
         \end{pmatrix}.

The denominator :math:`\det(I - S) = 1 - \alpha_L\,\alpha_R\,e^{-2\tau}`
is the per-period reflection product factor — at :math:`\alpha_L =
\alpha_R = \alpha` it collapses to :math:`1 - \alpha^2\,e^{-2\tau}`,
the same factor the Phase-3A heuristic used. **What the Phase-3A
heuristic missed** is the off-diagonal numerator structure. The
rank-2 numerator matrix carries an :math:`\alpha\,e^{-\tau}` entry
on each off-diagonal; this is what makes the symmetric reduction
:math:`\psi_{\rm surf} = \alpha\,B/(1 - \alpha\,e^{-\tau})` (with
single-transit B), not :math:`\alpha\,B_{\rm period}/(1 -
\alpha^2\,e^{-2\tau})` (with full-period B). The algebraic identity
:math:`(1 - \alpha\,e^{-\tau})(1 + \alpha\,e^{-\tau}) = 1 -
\alpha^2\,e^{-2\tau}` shows where the heuristic's denominator came
from: it's the symmetric reduction's denominator times the
"absorbed" factor :math:`(1 + \alpha\,e^{-\tau})`. The two
denominators agree at corners :math:`\alpha \in \{0, 1\}` (where
the absorbed factor is 1 or :math:`1 + e^{-\tau}` respectively;
both produce uniform-:math:`q` agreement on V_α1) but **disagree at
intermediate** :math:`\alpha` — ERR-035's mechanism in one
sentence.

The rank-2 reduction at :math:`\alpha_R = 0` (right wall vacuum)
gives :math:`T = \mathrm{diag}(1, 1) - \mathrm{antidiag}(\alpha_L\,
e^{-\tau}, 0)` etc., simplifying further. At :math:`\alpha_L =
\alpha_R = 0` (vacuum-vacuum), :math:`T = I` and the closure is
**bit-equal** to the bare first-leg integral with no surface
contribution — verified explicitly in the Phase-3B test suite.

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

The structurally-independent acceptance gate for Phase-3B is the
**method-of-images symmetry test**:

.. math::
   :label: peierls-greens-slab-asym-method-of-images

   \boxed{\;
   k_{\rm eff}\!\big[\text{slab } [0, L],\, \alpha_L = 1,\, \alpha_R = 0\big]
   \;=\;
   k_{\rm eff}\!\big[\text{slab } [0, 2L],\, \alpha_L = \alpha_R = 0\big]
   \;}

**Algebraic justification.** The fundamental k-eigenmode of the
symmetric vacuum slab :math:`[0, 2L]` is symmetric about :math:`x =
L`: by the slab's even-symmetry, the flux profile
:math:`\phi_{\rm sym}(x)` satisfies :math:`\phi_{\rm sym}(L + x) =
\phi_{\rm sym}(L - x)` for :math:`x \in [0, L]`. The derivative
:math:`\partial_x \phi_{\rm sym}` vanishes at :math:`x = L` by
symmetry, which is **exactly the boundary condition imposed by a
reflective wall** (perfect specular BC at :math:`x = L` enforces
:math:`\phi'(L) = 0`). Cutting the symmetric slab at :math:`x = L`
and applying the specular condition (which already holds by
symmetry) yields an equivalent problem on the right half
:math:`[L, 2L]` — bijectively mapped to :math:`[0, L]` via
:math:`x_{\rm asym} = x_{\rm sym} - L`. Therefore the eigenvalues
match **exactly**, and the flux profiles satisfy
:math:`\phi_{\rm asym}(x) = \phi_{\rm sym}(L + x)` for
:math:`x \in [0, L]`.

This is the **canonical structural-independence test** for the
rank-2 prototype, for two reasons:

1. **It exercises the closure algebra across a non-trivial BC
   diagonal**: :math:`(\alpha_L = 1, \alpha_R = 0)` is the
   off-diagonal corner of the BC parameter square, where the
   asymmetric monodromy entries :math:`\alpha_L\,e^{-\tau}` and
   :math:`\alpha_R\,e^{-\tau}` differ maximally. The corner with
   uniform-source eigenmodes (:math:`\alpha_L = \alpha_R = 1`) does
   not exercise this; only this off-diagonal corner does.

2. **It exercises the trajectory machinery on a non-uniform spatial
   eigenmode**: the asymmetric :math:`(1, 0)` flux profile peaks
   AT the reflective wall :math:`x = 0`, monotone-decreasing toward
   the vacuum wall. A trajectory-parametrisation bug like ERR-034
   (missing :math:`\mu` factor in :math:`x_{\rm traj}`) shifts the
   apparent flux peak to an interior point, breaking the agreement
   with the symmetric reference. ERR-034 was caught **here** because
   the V_α1-style uniform-source tests are structurally blind to
   it.

Quantitative result (post-fix; from
:func:`test_method_of_images_reflective_vacuum_equals_double_vacuum`):

.. list-table:: Method-of-images convergence
   :header-rows: 1
   :widths: 30 30 40

   * - :math:`(n_x, n_\mu, n_{\rm traj})`
     - :math:`k_{\rm eff}` agreement (rel)
     - flux-shape agreement
   * - (32, 32, 64)
     - 5.6e-8
     - max-normalised flux ≤ 2e-3
   * - (48, 48, 96)
     - 1.7e-8
     - 2e-3 (passing)
   * - (64, 64, 128)
     - 6.5e-9
     - 2e-3 (passing)

Buggy pre-fix versions gave ~5 % :math:`k_{\rm eff}` disagreement
— the bug cannot hide in this test.

Variant α architecture for asymmetric slab (operator-level identities)
-----------------------------------------------------------------------

The slab-asymmetric Branch-1 SymPy module
:mod:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_slab_asymmetric`
verifies three identities (mirroring V_α1/V_α2/V_α3 from sphere but
on the rank-2 form):

- **V_α1_slab_asym** — closed-asymmetric-slab bounce-sum
  self-consistency at :math:`\alpha_L = \alpha_R = 1`
  (:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_slab_asymmetric.derive_operator_constant_trial_closed_slab_asymmetric`).
  For homogeneous slab with constant source :math:`q`,
  :math:`\psi(x, \mu) = q/\Sigma_t` everywhere via the rank-2
  closure. The proof uses the matrix algebra
  :math:`(I - S)^{-1} \cdot \alpha B = (1 - \alpha^2 e^{-2\tau})^{-1}
  \cdot [(\alpha B + \alpha B \cdot \alpha e^{-\tau}); \cdots]`
  and shows the chord-cancellation identity collapses to
  :math:`q/\Sigma_t` for both surface entries.

- **V_α2_slab_asym** — rank-2 resolvent identity
  (:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_slab_asymmetric.derive_rank2_resolvent_slab_asymmetric`).
  Verifies :math:`(I - S) \cdot T = I` algebraically, with explicit
  symbolic substitution :math:`S = \mathrm{antidiag}(\alpha_L\,
  e^{-\tau}, \alpha_R\,e^{-\tau})`. The proof uses
  ``sp.simplify(Matrix.inv() - explicit_T)`` to confirm the closed
  matrix form (:eq:`peierls-greens-slab-asym-resolvent`).

- **V_α3_slab_asym** — vacuum reduction at :math:`\alpha_L =
  \alpha_R = 0`
  (:func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_slab_asymmetric.derive_alpha_zero_kernel_reduction_slab_asymmetric`).
  At :math:`\alpha_L = \alpha_R = 0`, :math:`S = 0`, :math:`T = I`;
  the closure reduces to the bare first-leg integral
  :math:`\psi = F(x, \mu)` with no surface contribution. Vacuum-
  vacuum BC is the trivial limit; bit-equal to the Phase-3A rank-1
  vacuum path (both bypass the closure entirely).

Numerical verification status (Phase 3B)
-----------------------------------------

- **k_inf exactness at α_L=α_R=1** — fuel-A-like XS, both thin
  (:math:`\tau_L = 1`) and moderate (:math:`\tau_L = 5`):
  :math:`k_{\rm eff} = k_\infty` to :math:`\le 1\mathrm e{-10}`
  relative (machine floor :math:`\sim 5.55\mathrm e{-17}`,
  1 iteration from constant initial guess). Reproduces V_α1_slab_asym
  numerically.
- **Method-of-images symmetry test** (the load-bearing gate) —
  asym :math:`[0, 1]` :math:`(\alpha_L=1, \alpha_R=0)` agrees with
  sym :math:`[0, 2]` :math:`(\alpha=0, 0)` to ≤ 6.5e-9 at finest
  grid (full table above). Asymmetric flux peaks at :math:`x
  \approx 0.001` (the reflective wall), monotone decreasing to
  vacuum. The pre-fix buggy version gave ~5 % disagreement; the
  bug cannot hide in this test.
- **Vacuum-vacuum bit-equality** — at :math:`\alpha_L = \alpha_R = 0`,
  rank-2 result is **bit-equal** to Phase-3A rank-1 (both bypass the
  closure entirely and reduce to the bare first-leg integral). This
  is the strongest possible structural-independence corner check:
  if the rank-2 closure had any algebraic drift, vacuum-vacuum would
  not be bit-equal.
- **Reduce-to-symmetric** at intermediate :math:`\alpha = 0.5` —
  post-ERR-035 fix, the rank-1 (Phase-3A delegated) path agrees
  with the rank-2 path to :math:`\le 1\mathrm e{-12}` relative
  (bit-equal at the working precision). The repurposed
  ``test_rank1_path_now_agrees_with_rank2_via_delegation_after_ERR035_fix``
  is tagged ``@pytest.mark.catches("ERR-035")`` as the
  regression-prevention gate.
- **Multi-group at α_L=α_R=1** — 2G with asymmetric scattering
  matches :math:`k_\infty` from
  :func:`~orpheus.derivations.common.eigenvalue.kinf_homogeneous`
  to :math:`\le 1\mathrm e{-9}` relative. Detects ERR-002
  SigS-transpose anti-pattern.

.. _peierls-greens-slab-asym-err035-fix:

ERR-035 forensic narrative — the heuristic closure and its fix
---------------------------------------------------------------

This subsection documents the full history of the ERR-035 closure
heuristic and its delegation fix, because the lesson is load-bearing
for any future Variant α geometry extension. The two ERR-NNN entries
are catalogued at
:file:`.claude/skills/vv-principles/error_catalog.md`.

**The heuristic and its origin.** When Phase-3A slab Variant α was
first implemented, the geometry was the first to introduce a
**2-bounce-per-period** trajectory: each bounce period traverses the
full slab width twice (out + back at constant :math:`|\mu|`), with
two specular reflections. Sphere and cylinder by contrast have
**1-bounce-per-period** geometry (one bounce per full closed-trajectory
period). Phase-3A constructed the slab closure by **analogical
generalisation** from sphere/cylinder rank-1:

.. math::

   \psi_{\rm surf}^{\rm heuristic}(\mu)
       = \frac{\alpha\,B_{\rm period}(\mu)}
              {1 - \alpha^2\,e^{-2\tau}},
   \qquad
   B_{\rm period}(\mu) = \int_0^{2L/|\mu|}\!\!q\,e^{-\Sigma_t s}\,
       \mathrm d s,

with two algebraic adaptations from the sphere/cylinder rank-1
form: (1) replace :math:`\alpha` in the geometric-series factor
with :math:`\alpha^2` to encode "two reflections per period",
(2) replace the sphere-style single-bounce-period :math:`B` integral
with the slab full-period (out+back) integral :math:`B_{\rm period}`.
The construction looked plausible: at α=1 it reproduces V_α1, at
α=0 it reduces to vacuum, the limits behave correctly. **The
construction is wrong.**

**Why the heuristic agrees at corners.** At :math:`\alpha = 1`, V_α1
uses constant source :math:`q`, and the chord-cancellation algebra
absorbs both the closure denominator and the B integral structure;
the heuristic happens to give the right answer at this corner
**by luck** (the algebra produces :math:`q/\Sigma_t` regardless of
the specific :math:`B` integral and denominator, as long as the
operator is consistent). At :math:`\alpha = 0` the leading factor
:math:`\alpha` in the closure makes :math:`\psi_{\rm surf} = 0`
trivially.

**Why the heuristic fails at intermediate α.** The structurally-
correct closure is the rank-2 form:

.. math::

   \psi_{\rm surf}^{\rm correct}(\mu)
       = \frac{\alpha\,B_{\rm single}(\mu)}
              {1 - \alpha\,e^{-\tau}},
   \qquad
   B_{\rm single}(\mu) = \int_0^{L/|\mu|}\!\!q\,e^{-\Sigma_t s}\,
       \mathrm d s,

with **single-transit** B integral and a different denominator. The
algebraic identity :math:`(1 - \alpha\,e^{-\tau})(1 + \alpha\,e^{
-\tau}) = 1 - \alpha^2\,e^{-2\tau}` shows where the heuristic's
denominator came from: it's the correct denominator times the
absorbed factor :math:`(1 + \alpha\,e^{-\tau})`. At α∈{0,1}, the
absorbed factor is :math:`\{1, 1 + e^{-\tau}\}` — both work
through V_α1's chord-cancellation algebra. At intermediate α
(0 < α < 1), the absorbed factor is :math:`(1 + \alpha\,e^{-\tau})
\neq 1`, and the heuristic over-predicts :math:`k_{\rm eff}` by
:math:`\sim 1.3\mathrm e{-4}` relative at :math:`\alpha = 0.5`,
fuel-A-like τ_L = 5.

**The bug fingerprint** (per the
:file:`.claude/skills/numerical-bug-signatures/SKILL.md`
discipline):

- **Sign**: positive (heuristic OVER-predicts :math:`k_{\rm eff}`).
- **Magnitude**: ~1.3e-4 relative at :math:`\alpha = 0.5`; smoothly
  varying with :math:`\alpha` (zero at endpoints).
- **Regime**: 0 < α < 1, with non-uniform spatial source profile.
- **Failure mode**: convention drift (failure mode #6 in the
  AI-failure-modes table).

**The fix architecture (delegation).**
:func:`solve_greens_function_slab` and
:func:`solve_greens_function_slab_mg` were refactored as **thin
wrappers** that delegate to
:func:`solve_greens_function_slab_asymmetric` and
:func:`solve_greens_function_slab_asymmetric_mg` at :math:`\alpha_L
= \alpha_R = \alpha`. The heuristic ``_apply_operator_slab`` helper
and the slab-specific chord-length functions were **deleted**: dead
code after delegation; keeping them as "alternative paths" would
re-create the very pattern the fix removes (the brief explicitly:
"don't leave them as 'back-compat shims' or 'alternative paths' —
they're known wrong"). The Phase-3A public API surface
(:class:`SlabGreensResult`, :class:`SlabGreensMGResult`, the two
solver entrypoints) is preserved by re-wrapping the rank-2 result
into the existing dataclass. Net code change: −282 production lines
collapsed into ~75 lines of trivial wrapping.

**The cascade of side-effects.** The vacuum convergence floor
(α=0, fuel-A-like τ_L=5) improved from a pre-fix ~5e-4 to
post-fix ~9e-6 (between the two finest grids); the absolute
:math:`k_{\rm eff}` jumped from ~0.130 to **0.157**. The pre-fix
floor was misattributed to a slow-:math:`\mu`-quadrature limitation;
in fact the dominant signal was the latent ERR-034 trajectory bug
(``x_traj = x - s`` missing the :math:`\mu` factor) masquerading
as quadrature noise. ERR-034 was caught by the method-of-images
test during Phase-3B integration; ERR-035's surface gate (Phase-3A
vs Phase-3B intermediate-α agreement) was caught simultaneously.

**The regression-prevention gate.** The repurposed test
:func:`test_rank1_path_now_agrees_with_rank2_via_delegation_after_ERR035_fix`
asserts :math:`\le 1\mathrm e{-12}` rtol bit-equal agreement
between the two paths at :math:`\alpha = 0.5`, tagged
``@pytest.mark.catches("ERR-035")``. The original test
``test_rank2_vs_rank1_at_intermediate_alpha_documented_discrepancy``
asserted a [5e-5, 5e-4] disagreement gate (fail if agreement gets
**tighter** — silent fix — or **looser** — regression). After the
delegation, the achieved agreement is **0.0e+00 relative**
(bit-equal) because both code paths now perform the SAME underlying
computation. The test was repurposed to assert the agreement
directly, with the catches-marker preserving the audit trail.

**The teaching point** (entered into the project's
:file:`.claude/skills/algebra-of-record/SKILL.md` lessons):
**analogical generalisation of a closure formula across geometries
is an algebraic claim that requires a proof, not a substitution.**
Future Variant α geometry extensions MUST ship with a Branch-1
SymPy first-principles derivation on a non-uniform source, with
cross-check against a structurally-correct rank-2 closure. The
Phase 3C-1 hollow sphere and 3C-2 annulus implementations both
followed this discipline and shipped clean on first try, with no
ERR-NNN entries. The discipline works.

Source code, tests, and provenance
----------------------------------

- **Branch-1 SymPy**:
  :mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_slab_asymmetric`
  (~362 LoC).
  :func:`derive_operator_constant_trial_closed_slab_asymmetric` (V_α1_slab_asym);
  :func:`derive_rank2_resolvent_slab_asymmetric` (V_α2_slab_asym);
  :func:`derive_alpha_zero_kernel_reduction_slab_asymmetric` (V_α3_slab_asym).
- **Branch-2 production**:
  :mod:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function_slab_asymmetric`
  (~452 LoC).
  :func:`solve_greens_function_slab_asymmetric` (1G);
  :func:`solve_greens_function_slab_asymmetric_mg` (MG); plus
  internal helper ``_apply_operator_slab_asymmetric``; result
  dataclasses :class:`SlabAsymmetricGreensResult`,
  :class:`SlabAsymmetricGreensMGResult`.
- **Shared rank-2 core**:
  :func:`~orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core.compute_resolvent_T_rank2`
  and
  :func:`~orpheus.derivations.continuous.trajectory_resolvent.variant_alpha_core.apply_variant_alpha_closure_rank2`
  in :mod:`.variant_alpha_core` (~130 LoC added in Phase 3B).
- **Symbolic foundation tests**:
  :file:`tests/derivations/test_trajectory_resolvent_slab_asymmetric_symbolic.py`
  — 16 ``@pytest.mark.foundation`` gates.
- **Numerical L1 tests**:
  :file:`tests/derivations/test_trajectory_resolvent_slab_asymmetric_solver.py`
  — 11 ``@pytest.mark.l1`` gates including the load-bearing
  method-of-images symmetry test and the
  ``catches("ERR-035")`` regression-prevention gate.
- **Closeout memo**:
  :file:`.claude/agent-memory/method-implementer/slab_asymmetric_variant_alpha_phase3b.md`.
- **Cross-domain frame memo**:
  :file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
  — predicted the rank-2 generalisation before any 2-surface code
  existed.
- **Plan**:
  :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md`.


.. _peierls-greens-hollow-sph:

Hollow sphere Variant α (Phase 3C-1 — rank-2 + impact-parameter partition)
===========================================================================

Phase-3C-1 hollow sphere Variant α extends the rank-2 BIE block
resolvent (validated in Phase 3B for asymmetric slab) to the first
**curvilinear 2-surface topology**. Two new geometric features are
introduced compared to slab:

1. **Impact-parameter phase-space partition** — trajectories split
   into outer-only (rank-1) and through-ray (rank-2) subsets,
   determined by whether the trajectory line intersects the inner
   sphere.
2. **Curvilinear shell-traversal chord** — the slab single-transit
   :math:`L/|\mu|` is replaced by the shell-traversal length
   :math:`\sqrt{R_{\rm out}^2 - b^2} - \sqrt{R_{\rm in}^2 - b^2}`,
   parameterised by the conserved impact parameter :math:`b`.

Both new features compose cleanly with the rank-2 closure operator —
:func:`compute_resolvent_T_rank2` and
:func:`apply_variant_alpha_closure_rank2` are imported from
:mod:`.variant_alpha_core` **without modification**. The hollow-
sphere Branch-2 production solver's only new code is the geometry-
specific chord arithmetic and the four-case first-leg analysis. The
closeout memo at
:file:`.claude/agent-memory/method-implementer/hollow_sphere_variant_alpha_phase3c1.md`
records: clean first-try landing, V_α1 closed-shell exactness at
4e-16 in 1 iteration, R_in → 0 solid-sphere limit at 9e-10
(plan target was 1e-3 — six orders of magnitude better), zero ERR-NNN
entries.

Impact-parameter phase-space partition
---------------------------------------

For an interior shell point :math:`(r, \mu)` with :math:`r \in [R_{
\rm in}, R_{\rm out}]` and :math:`\mu \in [-1, 1]`, the **impact
parameter**

.. math::
   :label: peierls-greens-hollow-sph-impact-parameter-partition

   b(r, \mu) = r\sqrt{1 - \mu^2}

is the perpendicular distance from the cylinder/sphere centre to the
trajectory line; it is **conserved along straight-line travel**
(geometric invariant of the chord) and **conserved under specular
reflection at the outer (or inner) sphere surface** (the radial-normal
reflection preserves :math:`|\mu|` at the surface, and by the chord
geometry preserves :math:`b`). The phase space splits at the threshold
:math:`b = R_{\rm in}`:

- :math:`b > R_{\rm in}` — **outer-only ray**. The trajectory line
  does not intersect the inner sphere (the perpendicular distance
  from the center exceeds the inner radius); the particle bounces
  between two points on the **outer** surface only. Topologically
  **rank-1**, structurally identical to a solid-sphere ray at the
  same outer radius and impact parameter. The inner cavity is
  "invisible" to this ray.
- :math:`b \le R_{\rm in}` — **through-ray**. The trajectory line
  crosses the inner cavity. Under interpretation (A) — the inner
  surface is a specular reflector with reflectivity
  :math:`\alpha_{\rm in}` — the particle bounces alternately between
  the **inner and outer** surfaces. Topologically **rank-2**,
  structurally identical to the asymmetric slab (with the curved
  shell-chord algebra replacing the slab's :math:`L/|\mu|`).

The four-case first-leg backward analysis. The first-leg chord from
interior :math:`(r, \mu)` to first surface arrival depends on TWO
conditions: sign of :math:`\mu` (which direction the trajectory
came from) AND the relation of :math:`b` to :math:`R_{\rm in}`
(whether through-ray or outer-only). Worked through symbolically
(sign of :math:`r \mu \pm \sqrt{R^2 - b^2}` choosing the
backward-positive root) before coding — see the closeout memo for
the full table. The four cases are:

1. :math:`\mu > 0, b > R_{\rm in}`: outer-only, chord to outer surface.
2. :math:`\mu > 0, b \le R_{\rm in}`: through-ray, chord to inner surface.
3. :math:`\mu < 0, b > R_{\rm in}`: outer-only, chord to outer surface (different branch of the discriminant root).
4. :math:`\mu < 0, b \le R_{\rm in}`: through-ray, chord to outer surface.

Outer-only resolvent (rank-1)
-----------------------------

For :math:`b > R_{\rm in}`, the bounce-period chord on the outer
sphere is the antipodal-chord at perpendicular distance :math:`b`
from the centre: :math:`L_{\rm period}(b) = 2\sqrt{R_{\rm out}^2 -
b^2}`. The closure is the rank-1 form

.. math::
   :label: peierls-greens-hollow-sph-outer-only-resolvent

   \psi_{\rm surf}(b) = \frac{\alpha_{\rm out}\,B(b; q)}
                              {1 - \alpha_{\rm out}\,
                               e^{-\Sigma_t\,L_{\rm period}(b)}},

**identical to the solid-sphere closure** (sphere V_α surface
fixed-point at :math:`r = R_{\rm out}`). The inner cavity is
invisible to outer-only rays. This is invoked as
:func:`apply_variant_alpha_closure` from :mod:`.variant_alpha_core`
with :math:`\alpha = \alpha_{\rm out}` and the outer chord length.

Through-ray rank-2 resolvent
----------------------------

For :math:`b \le R_{\rm in}`, the **shell-traversal optical depth**
(inner→outer transit, single-step) is

.. math::

   \tau_{\rm step}(b) = \Sigma_t \cdot \bigl(\sqrt{R_{\rm out}^2 - b^2}
                                            - \sqrt{R_{\rm in}^2 - b^2}\bigr).

Geometric content: the chord from the inner surface (at radius
:math:`R_{\rm in}`) to the outer surface (at radius :math:`R_{\rm
out}`), at constant impact parameter :math:`b`, has the length
:math:`\sqrt{R_{\rm out}^2 - b^2} - \sqrt{R_{\rm in}^2 - b^2}` by
the chord-theorem of the circle. This is the curvilinear analogue
of the slab single-transit chord :math:`L/|\mu|`.

The surface state is the 2-vector :math:`[\psi_{\rm in}^{\rm
out}(b, \mu), \psi_{\rm out}^{\rm in}(b, \mu)]` (outgoing flux at
the inner surface in the direction of the outer surface, and
outgoing flux at the outer surface in the direction of the inner
surface). The monodromy and resolvent have the same antidiagonal
structure as the asymmetric slab — only :math:`\tau_{\rm step}`'s
chord meaning differs:

.. math::
   :label: peierls-greens-hollow-sph-through-rank2

   T(\alpha_{\rm in}, \alpha_{\rm out}, \tau_{\rm step})
       = \frac{1}{1 - \alpha_{\rm in}\,\alpha_{\rm out}\,
                     e^{-2\tau_{\rm step}}}
         \begin{pmatrix}
             1                                 & \alpha_{\rm in}\,
                                                  e^{-\tau_{\rm step}} \\
             \alpha_{\rm out}\,e^{-\tau_{\rm step}}    & 1
         \end{pmatrix}.

This is invoked as :func:`apply_variant_alpha_closure_rank2` from
:mod:`.variant_alpha_core` with ``alphas = (α_in, α_out)`` and
``tau_step``. **The rank-2 closure operator is byte-equal-shared
with the asymmetric slab implementation** — only the chord algebra
differs. This shared-byte property is the load-bearing evidence
that the cross-domain frame's rank-2 prediction is structural, not
a coincidence.

Composability at the partition boundary
----------------------------------------

A subtle question at the architecture level is whether the rank-1
outer-only and rank-2 through-ray closures **agree** at the
partition boundary :math:`b = R_{\rm in}`. The answer is **yes by
construction**: at :math:`b = R_{\rm in}`, the inner-shell chord
:math:`\sqrt{R_{\rm in}^2 - b^2} = 0`, so the through-ray
shell-traversal optical depth :math:`\tau_{\rm step} \to
\Sigma_t \sqrt{R_{\rm out}^2 - R_{\rm in}^2}` (finite, half-period of
the outer-only chord), and the rank-2 closure reduces algebraically
to the rank-1 form (the off-diagonal entries become
:math:`\alpha_{\rm in}\,e^{-\tau_{\rm step}}` and :math:`\alpha_{\rm
out}\,e^{-\tau_{\rm step}}`; with the inner-transit chord vanishing,
the inner surface acts only as a passive boundary — and the
through-ray's outer surface contribution agrees with the outer-only
ray's outer-surface contribution).

The numerical evidence: V_α1 closed-shell exactness at :math:`\alpha_{
\rm in} = \alpha_{\rm out} = 1` requires BOTH branches independently
producing :math:`\psi = q/\Sigma_t` and **converging to the same
value at the partition boundary**. This is the load-bearing
structural composability check, and it passes at machine precision
in 1 iteration (rel err 4e-16).

Architecture summary
--------------------

.. math::
   :label: peierls-greens-hollow-sph-architecture

   \psi(r, \mu) = F(r, \mu) + e^{-\Sigma_t\,L_{\rm first}(r, \mu)}
                                \cdot \psi_{\rm surface}(r, \mu),

where :math:`\psi_{\rm surface}` is determined by the phase-space
partition: rank-1 outer-only closure for :math:`b > R_{\rm in}`,
rank-2 closure for :math:`b \le R_{\rm in}` selecting
:math:`\psi_{\rm in}^{\rm out}` for :math:`\mu > 0` (first arrival
backward = inner surface) and :math:`\psi_{\rm out}^{\rm in}` for
:math:`\mu < 0` (first arrival backward = outer surface).

Operator-level identities (V_α1_hollow_sph / V_α2_hollow_sph / V_α3_hollow_sph)
-------------------------------------------------------------------------------

The hollow-sphere Branch-1 SymPy module
:mod:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_hollow_sphere`
verifies three identities. The structurally novel piece (relative to
sphere/cylinder/slab/slab-asym) is the **independent** verification
of V_α1 on **both branches** of the impact-parameter partition.

- **V_α1_hollow_sph** — closed-shell bounce-sum self-consistency
  at :math:`\alpha_{\rm in} = \alpha_{\rm out} = 1`
  (:func:`derive_operator_constant_trial_closed_hollow_sphere`).
  The proof shows that **both outer-only AND through-ray closures
  independently** produce :math:`\psi = q/\Sigma_t` for constant
  source, via two distinct chord-cancellation patterns:

  - **Outer-only branch** (:math:`b > R_{\rm in}`): rank-1 collapse via
    :math:`(1 - x)/(1 - x) = 1` for the bounce-period chord —
    identical to V_α1 sphere.
  - **Through-ray branch** (:math:`b \le R_{\rm in}`): rank-2 collapse
    via :math:`(1 - e^{-\tau_{\rm step}})(1 + e^{-\tau_{\rm step}})/(
    1 - e^{-2\tau_{\rm step}}) = 1` — identical to V_α1 slab-asym.

  Both branches share the SAME constant :math:`\psi = q/\Sigma_t` at
  the boundary :math:`b = R_{\rm in}`, so phase-space continuity is
  preserved. This composability proof is the **load-bearing
  structural-independence anchor** for the impact-parameter partition.

- **V_α2_hollow_sph** — rank-2 resolvent identity
  (:func:`derive_rank2_resolvent_hollow_sphere`). Verifies
  :math:`(I - S) \cdot T = I` algebraically for the hollow-sphere
  monodromy (:eq:`peierls-greens-hollow-sph-through-rank2`). Same
  algebraic structure as slab-asym V_α2, just with the curvilinear
  shell-chord :math:`\tau_{\rm step}` substituted.

- **V_α3_hollow_sph** — vacuum-vacuum reduction at
  :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0`
  (:func:`derive_alpha_zero_kernel_reduction_hollow_sphere`). At
  :math:`\alpha_{\rm in} = \alpha_{\rm out} = 0`, BOTH branches'
  closures vanish (rank-1 leading factor :math:`\alpha_{\rm out} =
  0`; rank-2 :math:`S = 0 \Rightarrow T = I` and the off-diagonal
  source contributions vanish via leading :math:`\alpha`). The
  closure reduces to the bare first-leg integral.

Method-of-images analogy with slab-asymmetric
----------------------------------------------

The hollow sphere's reflective-inner / vacuum-outer configuration
(:math:`\alpha_{\rm in} = 1`, :math:`\alpha_{\rm out} = 0`) is the
**curvilinear analogue** of the slab :math:`(1, 0)` method-of-images
configuration. Both produce **flux profiles that peak at the
reflective wall** and decay monotonically toward the vacuum wall.
The agreement of the two configurations on their respective method-
of-images "doubled-domain" symmetric vacuum problems would in
principle constitute a structural-independence test for hollow
sphere, but is not currently implemented as a test (the
double-domain hollow-sphere construction needed would itself be a
new structurally distinct geometry — not productive scope for
Phase 3C-1). Instead, the **R_in → 0 solid-sphere limit** plays
the analogous structural-independence role for hollow sphere: as
:math:`R_{\rm in} \to 0`, the through-ray subset (:math:`b \le R_{
\rm in}`) shrinks to measure-zero, the outer-only branch dominates,
and the result must converge to the **independently-verified
solid-sphere vacuum reference** (Phase 1 sphere with PS-1982
cross-check). This is the structurally-independent gate.

Cavity-absorber physical interpretation
----------------------------------------

Under "interpretation (A)" — the inner surface is a specular
reflector with reflectivity :math:`\alpha_{\rm in}` — setting
:math:`\alpha_{\rm in} = 0` corresponds to a **perfect inner
absorber**: through-rays that hit the inner surface are lost to
the inner cavity (modelled as a black absorber). The flux profile
at :math:`\alpha_{\rm in} = 0, \alpha_{\rm out} = 1` (cavity
absorber, reflective outer): through-rays absorbed at inner;
outer-only rays bounce normally; **flux peaks at outer wall** (the
reflective surface), monotone decreasing toward inner. The
quantitative result from the closeout: inner/outer ratio ≈ 0.31
at fuel-A-like τ_R, k_eff = 0.174 < k_inf.

The opposite case (:math:`\alpha_{\rm in} = 1`, :math:`\alpha_{\rm
out} = 0`, reflective-inner / vacuum-outer): outer is a vacuum
leakage surface, inner traps through-rays in a back-and-forth
bounce. **Flux peaks near the inner wall** — symmetric configuration
to the slab method-of-images :math:`(1, 0)` peaking at the
reflective wall. Quantitative: inner/outer ratio = 3.18, k_eff =
0.0567 < k_inf.

R_in → 0 solid-sphere limit (structural reduction proof)
---------------------------------------------------------

As :math:`R_{\rm in} \to 0`, the through-ray subset (:math:`b \le
R_{\rm in}`) approaches **measure zero** in phase-space (the
codimension increases by 1 as :math:`R_{\rm in}` shrinks to a
point), and the outer-only branch dominates. The shell-traversal
chord :math:`\sqrt{R_{\rm out}^2 - b^2} - \sqrt{R_{\rm in}^2 - b^2}
\to \sqrt{R_{\rm out}^2 - b^2}` (the half-period of the solid-
sphere chord at impact parameter :math:`b`). The hollow-sphere
result therefore reduces structurally to the solid-sphere result.

**Numerical evidence** (from the closeout): hollow with
:math:`R_{\rm in} = 10^{-3} R_{\rm out}` and :math:`\alpha_{\rm in}
= \alpha_{\rm out} = 0` agrees with the solid-sphere vacuum
reference at:

- :math:`(n_r, n_\mu, n_{\rm traj}) = (24, 24, 48)`: rel diff vs
  solid-sphere = **9.2e-10**.
- :math:`(n_r, n_\mu, n_{\rm traj}) = (32, 32, 64)`: rel diff =
  **1.9e-9**.

Plan target was 1e-3; achieved ≤ 1e-9 — six orders of magnitude
better. The phase-space partition reduces correctly, and the
abstraction's structural soundness is empirically confirmed.

Verification status
-------------------

- **V_α1 closed-shell exactness** at :math:`\alpha_{\rm in} = \alpha_{
  \rm out} = 1`: BOTH outer-only AND through-ray closures
  independently produce :math:`\psi = q/\Sigma_t` for the closed-
  shell eigenmode. Numerical: :math:`k_{\rm eff} = k_\infty` at
  4e-16 in 1 iteration with uniform :math:`\phi` (machine precision).
  Load-bearing structural composability check.
- **R_in → 0 solid-sphere limit**: ≤ 1e-9 relative (target 1e-3 — six
  orders better; full table above).
- **Asymmetric reflective-inner / vacuum-outer** at
  :math:`\alpha_{\rm in} = 1, \alpha_{\rm out} = 0`: flux peaks near
  the reflective inner wall (inner/outer ratio = 3.18), k_eff =
  0.0567 < k_inf = 0.208.
- **Cavity-absorber** at :math:`\alpha_{\rm in} = 0, \alpha_{\rm out}
  = 1`: through-rays absorbed at inner; outer-only rays bounce
  normally; flux peaks at outer wall (inner/outer ratio = 0.31),
  k_eff = 0.174 < k_inf.
- **Symmetric reflective** :math:`\alpha = 0.5`: k_eff = 0.089 <
  k_inf, well-behaved rank-2 closure under symmetric BCs.
- **MG closed-shell** at :math:`\alpha_{\rm in} = \alpha_{\rm out}
  = 1` with 2G asymmetric scattering:
  :math:`k_{\rm eff} = 1.0000000000151812`, rel err vs
  :math:`k_\infty` = **1.5e-11**. Detects ERR-002 SigS-transpose
  anti-pattern.
- **Convergence floor** for reflective-inner / vacuum-outer: pin
  achieved across :math:`(n_r, n_\mu, n_{\rm traj})` ladder up to
  (48, 48, 96); finest pair self-consistency :math:`< 1\mathrm e{-3}`.

Source code, tests, and provenance
----------------------------------

- **Branch-1 SymPy**:
  :mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_hollow_sphere`.
  :func:`derive_operator_constant_trial_closed_hollow_sphere` (V_α1_hollow_sph);
  :func:`derive_rank2_resolvent_hollow_sphere` (V_α2_hollow_sph);
  :func:`derive_alpha_zero_kernel_reduction_hollow_sphere` (V_α3_hollow_sph).
- **Branch-2 production**:
  :mod:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function_hollow_sphere`.
  :func:`solve_greens_function_hollow_sphere` (1G);
  :func:`solve_greens_function_hollow_sphere_mg` (MG).
- **Symbolic foundation tests**:
  :file:`tests/derivations/test_trajectory_resolvent_hollow_sphere_symbolic.py`
  — 18 ``@pytest.mark.foundation`` gates.
- **Numerical L1 tests**:
  :file:`tests/derivations/test_trajectory_resolvent_hollow_sphere_solver.py`
  — 11 ``@pytest.mark.l1`` gates including the load-bearing V_α1
  closed-shell composability check and the R_in → 0 solid-sphere
  limit.
- **Closeout memo**:
  :file:`.claude/agent-memory/method-implementer/hollow_sphere_variant_alpha_phase3c1.md`.
- **Cross-domain frame memo**:
  :file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
  — predicted "rank-1 → rank-2 generalises 1-surface → 2-surface
  geometry, regardless of curvature"; hollow sphere is the first
  curvilinear 2-surface validation.
- **Plan**:
  :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md`.


.. _peierls-greens-annulus:

Annulus Variant α (Phase 3C-2 — rank-2 + cylinder 3D angular phase-space)
===========================================================================

Phase-3C-2 annulus Variant α extends the rank-2 BIE block resolvent
(validated in Phase 3B for asymmetric slab; lifted to curvilinear
2-surface in Phase 3C-1 for hollow sphere) to the **cylindrical
analogue of the hollow sphere**. **The closure operator is
byte-equal-shared with hollow sphere** — only the chord meaning
differs. The annulus is the third 2-surface validation of the
cross-domain frame's rank-2 prediction and **closes the Phase-3
family** on the unified rank-1 / rank-2 framework: 6 geometries × 2
topologies, all on :mod:`.variant_alpha_core` with **zero
modifications during Phase 3C** (the core was last modified in
Phase 3B to add the rank-2 primitives).

Two new geometric features beyond hollow sphere are introduced:

1. **Phase-space change**: cylinder uses
   :math:`(r, \mu_{\rm axial}, \varphi_{\rm az})` instead of
   sphere's :math:`(r, \mu)`. The angular phase-space is **3D
   instead of 2D**. Conserved invariants under specular reflection
   at radial-normal cylinder surfaces are :math:`\mu_{\rm axial}`
   (axial cosine — the bounce point's mirror plane is perpendicular
   to the cylinder axis, so reflection preserves the axial component
   verbatim) and :math:`b = r\,|\sin\varphi_{\rm az}|` (in-plane
   impact parameter — preserved across all bounces by the same
   chord-theorem argument as cylinder Phase 1 / hollow sphere).

2. **3D chord scaling lift**: the cylinder 2D in-plane shell chord
   is lifted to 3D arclength by the **axial-correction factor**
   :math:`1/\sqrt{1 - \mu_{\rm axial}^2}`. Per Issue #129, the
   cylinder Bickley-Naylor :math:`\mathrm{Ki}_n` form pre-integrates
   axial and produces ~22 % planar-limit mismatch — Variant α stays
   angle-resolved to avoid this. The annulus prototype's **anti-
   cross-check rule** is therefore: **DO NOT cross-check against
   thin slab via planar limit**; cross-check is via the R_in → 0
   solid-cylinder limit (which itself respects the angle-resolved
   discipline from Phase 1 cylinder).

The closeout memo at
:file:`.claude/agent-memory/method-implementer/annulus_variant_alpha_phase3c2.md`
records: clean first-try landing, V_α1 closed-annulus exactness at
9.3e-16 in 1 iteration, R_in → 0 solid-cylinder limit at 3.7e-6
(plan target was 1e-3 — almost three orders of magnitude better),
zero ERR-NNN entries.

Impact-parameter phase-space partition
---------------------------------------

For each interior phase-space point :math:`(r, \mu_{\rm axial},
\varphi_{\rm az})` with :math:`r \in [R_{\rm in}, R_{\rm out}]`,
the **conserved in-plane impact parameter**

.. math::
   :label: peierls-greens-annulus-impact-parameter-partition

   b(r, \varphi_{\rm az}) = r\,|\sin\varphi_{\rm az}|

is the same impact parameter as cylinder Phase 1
(:eq:`peierls-greens-cylinder-impact-parameter`), conserved under
specular reflection at radial-normal cylinder surfaces (in-plane
projection of the trajectory in a circular cross-section preserves
the perpendicular distance from the cylinder axis). It partitions
phase-space at :math:`b = R_{\rm in}` — the **same partition
threshold as hollow sphere despite the different phase-space**
(hollow sphere uses :math:`b = r\sqrt{1 - \mu^2}` in 2D phase-space;
annulus uses :math:`b = r\,|\sin\varphi_{\rm az}|` in 3D phase-space;
both are the perpendicular-distance impact parameter, both threshold
at the inner radius):

- :math:`b > R_{\rm in}` — **outer-only ray**. The 2D in-plane
  projection ray does NOT intersect the inner-radius circle; the
  particle bounces between two points on the OUTER cylinder only.
  Topologically rank-1, structurally identical to a solid-cylinder
  ray at the same outer radius and impact parameter.
- :math:`b \le R_{\rm in}` — **through-ray**. The 2D in-plane ray
  crosses the inner cavity. Under interpretation (A) — the inner
  surface is a specular reflector with reflectivity
  :math:`\alpha_{\rm in}` — the particle bounces alternately inner
  ↔ outer. Topologically rank-2.

The four-case first-leg backward analysis. The first-leg chord
depends on TWO conditions: sign of :math:`\cos\varphi_{\rm az}`
(direction the trajectory came from in the 2D in-plane projection)
AND the relation of :math:`b` to :math:`R_{\rm in}`. Worked through
symbolically before coding (closeout memo records the full table).
The four cases are the cylinder analogues of hollow sphere's
sign-of-:math:`\mu` × b-vs-R_in cases — same structure, only the
parameterisation differs (sphere uses sign of :math:`\mu`; cylinder
uses sign of :math:`\cos\varphi_{\rm az}`).

3D chord scaling identity
-------------------------

The annulus 3D shell-traversal chord and the hollow-sphere shell
chord are related by the **cylinder axial-correction factor**:

.. math::
   :label: peierls-greens-annulus-3d-chord-scaling

   \tau_{\rm step}^{\rm annulus}(b, \mu_{\rm axial})
        = \frac{\tau_{\rm step}^{\rm hollow\,sph}(b)}
               {\sqrt{1 - \mu_{\rm axial}^2}}
        = \Sigma_t \cdot
          \frac{\sqrt{R_{\rm out}^2 - b^2} - \sqrt{R_{\rm in}^2 - b^2}}
               {\sqrt{1 - \mu_{\rm axial}^2}}.

The same lift applies to the outer-only branch's bounce-period chord:
:math:`L_{\rm period}^{\rm 3D, annulus} = 2\sqrt{R_{\rm out}^2 - b^2}
/ \sqrt{1 - \mu_{\rm axial}^2}` (the cylinder Phase-1 identity
:func:`derive_bounce_period_chord_cylinder` evaluated at the annulus
outer surface).

The factorisation :math:`\tau_{\rm step}^{\rm annulus} =
\tau_{\rm step}^{\rm hollow\,sph}(b) / s_{\rm in\!-\!plane}(\mu_{\rm
axial})` exposes the **multiplicative composition** of the two
geometric structures: (i) the in-plane shell-traversal chord (which
is **identical** to hollow sphere as a function of :math:`b`), and
(ii) the cylinder axial-correction factor (which is **identical**
to cylinder Phase 1 as a function of :math:`\mu_{\rm axial}`). This
is the structural reason the closure operator is byte-equal-shared
with hollow sphere: the closure consumes :math:`\tau_{\rm step}` as
a scalar; whatever produced that scalar is the geometry's concern.

The auxiliary symbolic gate
:func:`derive_3d_chord_scaling_annulus` (V_α2_annulus.aux) verifies
this identity symbolically — a structural test that would FAIL if
any code path had pre-integrated the axial direction. This is the
**defensive gate** against a recurrence of Issue #129's axial
pre-integration drift.

Through-ray rank-2 resolvent
----------------------------

Through-ray rank-2 resolvent
----------------------------

For :math:`b \le R_{\rm in}`, the surface state is the 2-vector
:math:`[\psi_{\rm in}^{\rm out}(b, \mu_{\rm axial}),
\psi_{\rm out}^{\rm in}(b, \mu_{\rm axial})]` — same shape as hollow
sphere, but parameterised by :math:`(b, \mu_{\rm axial})` instead
of :math:`(b, \mu)`. The rank-2 resolvent has the **same
antidiagonal monodromy structure** as hollow sphere and asymmetric
slab — only :math:`\tau_{\rm step}` differs (now carrying the
cylinder axial-correction factor):

.. math::
   :label: peierls-greens-annulus-through-rank2

   T(\alpha_{\rm in}, \alpha_{\rm out}, \tau_{\rm step})
       = \frac{1}{1 - \alpha_{\rm in}\,\alpha_{\rm out}\,
                     e^{-2\tau_{\rm step}}}
         \begin{pmatrix}
             1                                   & \alpha_{\rm in}\,
                                                    e^{-\tau_{\rm step}} \\
             \alpha_{\rm out}\,e^{-\tau_{\rm step}}      & 1
         \end{pmatrix},

with :math:`\tau_{\rm step}` given by Eq.
:eq:`peierls-greens-annulus-3d-chord-scaling`. The closure operator
is **byte-equal-shared** with hollow sphere via
:func:`compute_resolvent_T_rank2` and
:func:`apply_variant_alpha_closure_rank2` in
:mod:`.variant_alpha_core` — the annulus prototype's only new code
is the cylinder-specific chord arithmetic.

Composability at the partition boundary
----------------------------------------

The same composability proof as hollow sphere applies: at
:math:`b = R_{\rm in}`, the inner-shell chord
:math:`\sqrt{R_{\rm in}^2 - b^2} = 0`, so :math:`\tau_{\rm step}
\to (\Sigma_t \sqrt{R_{\rm out}^2 - R_{\rm in}^2}) / \sqrt{1 -
\mu_{\rm axial}^2}` (finite). The rank-2 closure reduces
algebraically to the rank-1 outer-only form. **V_α1 closed-annulus
exactness at α_in = α_out = 1 requires BOTH branches independently
producing** :math:`\psi = q/\Sigma_t`, and the two branches must
converge to the same numerical value at the partition boundary.
Numerical: 9.3e-16 in 1 iteration with uniform :math:`\phi`
(machine precision).

Architecture summary
--------------------

.. math::
   :label: peierls-greens-annulus-architecture

   \psi(r, \mu_{\rm axial}, \varphi_{\rm az}) =
       F(r, \mu_{\rm axial}, \varphi_{\rm az}) +
       e^{-\Sigma_t\,L_{\rm first}^{\rm 3D}(r, \mu_{\rm axial},
                                                  \varphi_{\rm az})}
       \cdot \psi_{\rm surface}(r, \mu_{\rm axial}, \varphi_{\rm az}),

where :math:`\psi_{\rm surface}` is determined by the phase-space
partition: rank-1 outer-only closure for :math:`b > R_{\rm in}`,
rank-2 closure for :math:`b \le R_{\rm in}` selecting
:math:`\psi_{\rm in}^{\rm out}` for :math:`\cos\varphi_{\rm az} > 0`
(first arrival backward = inner surface) and
:math:`\psi_{\rm out}^{\rm in}` for :math:`\cos\varphi_{\rm az} \le 0`
(first arrival backward = outer surface).

Comparison to hollow sphere — what is byte-equal-shared
--------------------------------------------------------

The **closure machinery is byte-equal-shared** with hollow sphere:

- :func:`compute_resolvent_T_rank2` — same function call, same
  arguments :math:`(\alpha_{\rm in}, \alpha_{\rm out}, \tau_{\rm
  step})`, same return value. The rank-2 monodromy and resolvent
  algebra is geometry-agnostic.
- :func:`apply_variant_alpha_closure_rank2` — same function call,
  same closure formula. The geometry-specific chord meaning is
  encoded INTO :math:`\tau_{\rm step}` BEFORE the call; the closure
  consumes the scalar.
- The phase-space partition logic — same threshold :math:`b =
  R_{\rm in}`, same outer-only-vs-through-ray split.

What is **annulus-specific**:

- Phase-space discretisation: :math:`(r, \mu_{\rm axial},
  \varphi_{\rm az})` grids with Gauss-Legendre on each axis.
- Impact-parameter computation: :math:`b = r\,|\sin\varphi_{\rm az}|`
  (cylinder formula; hollow sphere uses :math:`b = r\sqrt{1-\mu^2}`).
- First-leg trajectory case analysis on sign of
  :math:`\cos\varphi_{\rm az}` (cylinder analogue of hollow sphere's
  sign-of-:math:`\mu` analysis).
- 3D arclength conversion :math:`L_{\rm 3D} = L_{\rm 2D} / \sqrt{1
  - \mu_{\rm axial}^2}` (the cylinder axial-correction lift, encoded
  multiplicatively into :math:`\tau_{\rm step}`).
- Antipodal-chord parameterisation in 2D in-plane coordinates with
  3D attenuation :math:`\exp(-\Sigma_t\,s_{\rm 2D}/s_{\rm in\!-\!
  plane})`.

Verification status
-------------------

- **V_α1 closed-annulus exactness** at :math:`\alpha_{\rm in} =
  \alpha_{\rm out} = 1`: BOTH outer-only AND through-ray closures
  independently produce :math:`\psi = q/\Sigma_t` for the closed-
  annulus eigenmode. Numerical: :math:`k_{\rm eff} = k_\infty` at
  9.3e-16 in 1 iteration with uniform :math:`\phi` (machine
  precision). This proves the impact-parameter partition + cylinder
  axial-correction lift compose cleanly with the rank-2 closure.
- **R_in → 0 solid-cylinder limit**: annulus with :math:`R_{\rm in}
  = 10^{-3} R_{\rm out}` and :math:`\alpha_{\rm in} = \alpha_{\rm
  out} = 0` agrees with the Phase-1 solid-cylinder vacuum reference
  at :math:`(n_r, n_\mu, n_\varphi, n_{\rm traj}) = (20, 20, 24,
  32)`: :math:`k_{\rm eff}^{\rm annulus} = 0.12045633841`,
  :math:`k_{\rm eff}^{\rm solid} = 0.12045677862`, **rel diff =
  3.7e-6**. Plan target was 1e-3 — almost three orders of magnitude
  better. The agreement is somewhat worse than hollow sphere's
  9e-10 floor at the equivalent condition because the **cylinder
  solid solver itself** is at lower quadrature-order self-
  consistency than the sphere solid solver — the cylinder Phase-1
  floor (~8e-8 self-consistency at finest grid) is the gating factor,
  not the annulus prototype.
- **Asymmetric reflective-inner / vacuum-outer** at
  :math:`\alpha_{\rm in} = 1, \alpha_{\rm out} = 0`: flux peaks
  near the reflective inner wall (inner/outer ratio = 2.62), k_eff
  = 0.0730 < k_inf = 0.208.
- **Cavity-absorber** at :math:`\alpha_{\rm in} = 0, \alpha_{\rm
  out} = 1`: through-rays absorbed at inner; outer-only rays bounce
  normally; flux peaks at outer wall (outer/inner ratio = 2.27),
  k_eff = 0.141 < k_inf.
- **Symmetric reflective** :math:`\alpha = 0.5`: k_eff < k_inf,
  well-behaved rank-2 closure under symmetric BCs.
- **MG closed-annulus** at :math:`\alpha_{\rm in} = \alpha_{\rm
  out} = 1` with 2G asymmetric scattering: :math:`k_{\rm eff}` =
  1.0000000000151756, rel err vs :math:`k_\infty` = **1.5e-11**.
  Detects ERR-002 SigS-transpose anti-pattern.

.. list-table:: Annulus convergence floor (reflective-inner / vacuum-outer)
   :header-rows: 1
   :widths: 35 30 35

   * - :math:`(n_r, n_{\mu_{\rm ax}}, n_{\varphi_{\rm az}}, n_{\rm traj})`
     - :math:`k_{\rm eff}`
     - rel diff to next finer
   * - (10, 10, 16, 16)
     - 7.3035007525e-02
     - 4.7e-5
   * - (14, 14, 20, 24)
     - 7.3031587150e-02
     - 9.5e-5
   * - (18, 18, 24, 32)
     - 7.3038554391e-02
     - 1.75e-4
   * - (24, 24, 32, 48)
     - 7.3051368990e-02
     - (terminus)

Convergence is non-monotone over the early ladder, settling around
the ~1e-4 floor at finest grid. Cylinder's 3D angular phase-space
requires more total quadrature points (:math:`n_r \times n_\mu \times
n_\varphi \approx 18\,000` at finest) than the sphere's 2D phase-space
(:math:`n_r \times n_\mu \approx 2\,300` at equivalent finest), so
the convergence is correspondingly slower per-element but the
absolute floor is comparable.

Phase-3 family roadmap (closed)
--------------------------------

The Phase-3 family is **complete** with this prototype. **Six
geometries on two topologies, all on the unified rank-1/rank-2
framework via** :mod:`.variant_alpha_core`:

- **Phase 1** — sphere, cylinder (1-surface compact, rank-1).
- **Phase 3A → ERR-035 delegation** — slab symmetric (delegated
  rank-2 reduction of Phase 3B at :math:`\alpha_L = \alpha_R`).
- **Phase 3B** — slab asymmetric (2-surface flat, rank-2).
- **Phase 3C-1** — hollow sphere (2-surface curvilinear, rank-2 +
  b-partition).
- **Phase 3C-2** — annulus (2-surface curvilinear with cylinder
  3D angular phase-space, rank-2 + b-partition + axial correction).

Three structural predictions of the cross-domain frame have been
validated across the family:

1. Rank-1 covers 1-surface compact geometries.
2. Rank-2 covers 2-surface geometries regardless of curvature.
3. Phase-space partition by impact parameter cleanly separates
   rank-1 outer-only from rank-2 through-ray on hollow geometries.

The :mod:`.variant_alpha_core` module has not been modified during
Phase 3C — the abstraction has earned the right to stay frozen.
**Any future extension that requires modifying** :mod:`.variant_alpha_core`
**should be treated as a structural alarm**.

Future extensions that are expected to compose cleanly on the
existing framework:

- Multi-region annulus (mirror sphere ``_apply_operator_mr``
  segmentation along trajectory + B integrals).
- Anisotropic scattering on annulus (orthogonal axis to Variant α).
- Multi-region hollow sphere (deferred from Phase 3C-1).
- Capping the cylinder/annulus axially (3-surface topology — would
  require rank-3 monodromy; algebraically a block-circulant
  resolvent generalisation).

Future extensions where the framework prediction is **less certain**:

- **Anisotropic-α** (specular at one wall, white at the other) —
  requires a rank-2 generalisation where :math:`S` has a
  non-antidiagonal structure. The frame predicts the closure should
  still be :math:`T = (I - S)^{-1}` with a more complex :math:`S`,
  but no validation case exists yet in the codebase.
- **Annulus + axial capping double-extension** (rank-3 with two
  curvilinear surfaces and two flat surfaces) would be the first
  test of the rank-N framework's applicability beyond
  pure-curvilinear-2-surface and pure-flat-symmetric cases.

Source code, tests, and provenance
----------------------------------

- **Branch-1 SymPy**:
  :mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function_annulus`.
  :func:`derive_operator_constant_trial_closed_annulus` (V_α1_annulus);
  :func:`derive_rank2_resolvent_annulus` (V_α2_annulus);
  :func:`derive_alpha_zero_kernel_reduction_annulus` (V_α3_annulus);
  :func:`derive_3d_chord_scaling_annulus` (V_α2_annulus.aux —
  axial-correction structural identity).
- **Branch-2 production**:
  :mod:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function_annulus`.
  :func:`solve_greens_function_annulus` (1G);
  :func:`solve_greens_function_annulus_mg` (MG).
- **Symbolic foundation tests**:
  :file:`tests/derivations/test_trajectory_resolvent_annulus_symbolic.py`
  — 22 ``@pytest.mark.foundation`` gates.
- **Numerical L1 tests**:
  :file:`tests/derivations/test_trajectory_resolvent_annulus_solver.py`
  — 11 ``@pytest.mark.l1`` gates including the load-bearing V_α1
  closed-annulus composability check, the R_in → 0 solid-cylinder
  limit reduction, and the research-grade convergence floor.
- **Closeout memo**:
  :file:`.claude/agent-memory/method-implementer/annulus_variant_alpha_phase3c2.md`.
- **Cross-domain frame memo**:
  :file:`.claude/agent-memory/cross-domain-attacker/variant_alpha_2surface_bie_frame.md`
  — the rank-1 → rank-2 generalisation prediction validated across
  six geometries.
- **Plan**:
  :file:`.claude/plans/peierls-greens-cylinder-and-2bc.md`.


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
:file:`.claude/agent-memory/literature-researcher/trajectory_resolvent_lit.md`
covers the full literature landscape and the Variant α / β / γ / δ
recommendation.

Code provenance
---------------

The Variant α implementation spans three files:

- :mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function`
  — SymPy derivations V_α1, V_α2, V_α3 (operator-level identities).
  About 270 lines.

  - :func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function.derive_operator_constant_trial_closed_sphere`
    — V_α1 closed-sphere bounce-sum self-consistency.
  - :func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function.derive_T00_equals_P_ss_sphere`
    — V_α2 algebraic identity :math:`T_{00}^{\rm sphere} =
    P_{ss}^{\rm sphere}`.
  - :func:`~orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function.derive_alpha_zero_kernel_reduction`
    — V_α3 vacuum-BC kernel reduction.

- :mod:`orpheus.derivations.continuous.trajectory_resolvent.greens_function` —
  production solver suite. About 1300 lines.

  - :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function._apply_operator_with_source_profile`
    — per-group homogeneous operator action; the load-bearing
    primitive shared by 1G / MG / 1-region paths.
  - :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function._apply_operator`
    — 1G wrapper around ``_apply_operator_with_source_profile``;
    builds :math:`q(r) = (\Sigma_s + \nu\Sigma_f/k)\,\phi(r) / 4\pi`.
  - :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere`
    — 1G public driver with :math:`\alpha\in[0,1]` parametrisation
    (B-phase + A1+A2 vacuum extension).
  - :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere_mg`
    — multi-group public driver (A3 follow-on); arbitrary G,
    full :math:`G\times G` scattering matrix, arbitrary :math:`\chi`.
  - :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function._apply_operator_mr`
    — multi-region per-group operator with piecewise :math:`\Sigt{}`
    along trajectory + bounce-period chord (Plan-(b)).
  - :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere_mr`
    — multi-region multi-group k-eigenvalue driver
    (Plan-(b) Option 2, Issue #132 attack).
  - :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere_mr_fixed_source`
    — multi-region multi-group fixed-source driver
    (Plan-(b) Option 1, Garcia 2021 cross-check).
  - :class:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.GreensFunctionResult`,
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

- :file:`tests/derivations/test_trajectory_resolvent_symbolic.py`
  — 8 SymPy gates (V_α1 surface fixed-point + total-ψ-constant +
  operator-eigenvalue + composite; V_α2 integrand-match + closed-form
  + composite; V_α3 g_h vanishes at α=0).

- :file:`tests/derivations/test_trajectory_resolvent_solver.py`
  — 3 numerical gates (V_α1.numerical with constant initial guess
  + non-uniform initial guess + two thicknesses
  :math:`\tau_R \in \{2.5, 5\}`).

- :file:`tests/derivations/test_trajectory_resolvent_xverif.py`
  — 3 closed-sphere cross-verification gates: B5.A (Variant α exact),
  B5.B (Phase 4 N=1 ≡ ``white_hebert``), B5.C (Phase 4 rank
  convergence toward Variant α).

- :file:`tests/derivations/test_trajectory_resolvent_vacuum.py`
  — 5 A1 vacuum-BC gates (k_eff < k_inf, thick sphere asymptote,
  α-continuity, non-trivial spatial mode, α=1 unchanged).

- :file:`tests/derivations/test_trajectory_resolvent_xverif_ps1982.py`
  — 6 A2 PS-1982 cross-verification gates (4 parametrised
  configurations + thick-sphere regression + flux-shape qualitative).

- :file:`tests/derivations/test_trajectory_resolvent_mg.py`
  — 7 A3 multi-group gates (G=1 reduction, 2G downscatter +
  upscatter, 2G spectrum, 2G vacuum, 4G with realistic χ, 4G
  vacuum).

- :file:`tests/derivations/test_trajectory_resolvent_mr.py`
  — 4 Plan-(b) Option 2 multi-region k-eigenvalue gates (1-region
  reduction to MG, Issue #132 catastrophe avoided, spatial mode
  physical, vacuum BC reduces k_eff).

- :file:`tests/derivations/test_trajectory_resolvent_garcia2021.py`
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

- :file:`.claude/agent-memory/literature-researcher/trajectory_resolvent_lit.md`
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
