.. _theory-peierls:

==================================================================
Peierls integral transport — method-agnostic foundations
==================================================================

.. contents:: Contents
   :local:
   :depth: 3


Key Facts
=========

**This page is the index for the two ORPHEUS Peierls implementation
families.** Both implementations solve the same continuous physics
(integral form of the steady-state monoenergetic transport equation
with isotropic scattering) but discretise it through *fundamentally
different operators*. Read this page first, then follow one of the
two forward links:

- :ref:`theory-peierls-nystrom` — the production
  **Nyström / matrix-Galerkin** architecture (slab + cyl + sphere).
  Discretises the **angle-integrated** kernel
  :math:`g_d(\rho'\to\rho)` on a radial Nyström grid; assembles the
  matrix :math:`K_{ij} = w_j\,g_d(r_j\to r_i)`; closes the BC via
  the tensor network :math:`K_{\rm bc} = G\cdot R\cdot P` with rank
  :math:`N` Marshak / F.4 / specular closures.
- :ref:`theory-peierls-greens` — the research-grade **Green's
  function (Variant α)** architecture (sphere homogeneous +
  multi-region, parametrised by reflectivity :math:`\alpha\in[0,1]`).
  Iterates the **angle-resolved** Green's function
  :math:`\tilde t(r'\to r,\mu)` along bouncing characteristics; the
  angle-integrated kernel is **never assembled**; the BC is absorbed
  into the kernel via Sanchez 1986 Eq. (A1).

The two architectures are *not* different discretisations of the
same operator. They target different operators that share the same
physical content:

- The Nyström operator collapses the angular structure into a 1-D
  scalar-flux integral equation; the angular details enter only
  through the pre-integrated kernel :math:`g_d`. This is the natural
  form for matrix-Galerkin closure (Hébert 2009 Eq. 3.323; rank-:math:`N`
  per-face F.4) and the production path for vacuum / white BC.
- The Green's function operator keeps the full angular variable
  :math:`\mu` and integrates it only at the very end of the iteration
  via :math:`\phi(r) = 2\pi\int\psi(r,\mu)\,\mathrm d\mu`. This is the
  natural form for problems where the BC absorption is non-trivial
  (specular multi-bounce; multi-region with strong scatterer outside
  fuel) — the cases that strain the Nyström rank-:math:`N` closures.

The Phase 5 retreat (Issue #133, CLOSED 2026-04-28; documented at
:ref:`peierls-phase5-retreat`) established that the
**angle-integrated** kernel :math:`g_\alpha` for specular BC is
*hypersingular* (Hadamard finite-part) — Nyström sampling diverges
at the diagonal and no quadrature trick rescues it. The
matrix-Galerkin form's mode-mixing absorbs the singularity via basis
projection, which is why ``boundary="specular_multibounce"`` works
(rank-N gating reflects the kernel's intrinsic difficulty, not a
basis-truncation artefact). The Green's function reformulation is
the structural fix: by working with the angle-resolved kernel along
characteristics, the Hadamard singularity is bypassed *structurally*
rather than by a quadrature trick.

For the canonical question "what reference do we ship for problem X?"
see the capability matrix at :ref:`theory-peierls-capabilities`
(inside :ref:`theory-peierls-nystrom`). The Nyström family is the
production reference for nearly all configurations; Variant α is the
research-grade reference for closed-sphere specular and multi-region
sphere where the Phase 4 rank-:math:`N` Marshak closure has
documented failure modes (Issue #132 — see
:ref:`peierls-rank-n-class-b-mr-mg-falsification`).


Motivation: why integral-transport / Peierls instead of S_N or MoC
==================================================================

The Peierls integral form is the spatial-integral analogue of the
discrete-ordinates / method-of-characteristics (MoC) sweeps used
elsewhere in ORPHEUS. The angular and free-flight degrees of freedom
have already been integrated out into the kernel
:math:`g_d(\rho'\to\rho)`; what remains is a single integral equation
in the spatial coordinate alone.

Why this matters for verification:

1. **Reference quality.** Adaptive ``mpmath.quad`` evaluation of
   :math:`g_d` gives machine-precision kernel values. Nyström assembly
   then produces a Peierls reference whose discretisation error is
   bounded by the radial-grid + closure-quadrature error, *not* by
   any angular quadrature error (which has already been folded into
   the kernel analytically). This is the architectural reason
   ``orpheus.derivations.continuous.peierls`` is the L1 reference
   class for the discrete CP / S_N / MoC solvers in
   ``orpheus.cp`` / ``orpheus.sn`` / ``orpheus.moc``.

2. **Geometry uniformity.** Slab, cylinder, and sphere reduce to the
   *same* polar-form integral equation
   :math:numref:`peierls-unified` once the kernel pre-integration is
   carried out. The remaining geometry-specific work is four
   primitives — angular measure, kernel function, ray-boundary
   distance, and source position — plugged into a shared scaffolding.
   This unification is the load-bearing observation of
   :ref:`theory-peierls-nystrom`.

3. **No ray effects.** The angular integration has already been
   performed analytically. There is no discrete-ordinates ray-effect
   contamination — the artefact pattern that plagues S_N at coarse
   :math:`N` and dominates the lattice-boundary signature at low
   :math:`S_N` orders. Peierls is the structurally correct reference
   for tests that probe ray-effect-related solver bugs.

4. **No flat-source assumption.** Standard CP (the
   :mod:`orpheus.cp` module) integrates the same equation under the
   *flat-source* assumption — the source is region-averaged before
   the kernel pre-integration. This collapses the second integration
   to a closed-form :math:`\mathrm{Ki}_n` / :math:`E_n` second
   difference at the cost of region-averaging error. The Peierls
   Nyström reference keeps the source pointwise (Lagrange-interpolated
   on a radial grid), so its error is one order-of-magnitude tighter
   than CP for the same regional discretisation.

What integral transport is *not* good for: realistic 2-D / 3-D
geometries, anisotropic scattering of high order, time-dependent
problems. Those go through MoC / S_N. The Peierls reference is the
1-D analytical-quadrature backbone behind ORPHEUS's verification
chain, not a production solver for survey calculations.


The transport equation in integral form
=========================================

The starting point is the steady-state, monoenergetic, isotropic-
emission transport equation for the angular flux
:math:`\psi(\mathbf r,\Omega)`,

.. math::
   :label: peierls-boltzmann

   \Omega\cdot\nabla\psi(\mathbf r,\Omega) + \Sigt{}\,
       \psi(\mathbf r,\Omega) = \frac{1}{4\pi}\,Q(\mathbf r),

where :math:`Q(\mathbf r) = \Sigs{}\,\phi(\mathbf r) +
(\nSigf{}/k)\,\phi(\mathbf r)` is the isotropic emission density
(scattering + fission, k-eigenvalue normalisation), and
:math:`\phi(\mathbf r) = \int\psi(\mathbf r,\Omega)\,\mathrm d\Omega`
is the scalar flux. Anisotropic scattering and the
:math:`\omega_1\,\Omega\cdot J(\mathbf r)` linear-anisotropic-source
term (Sanchez 1986 Eq. (1)) are out of scope on this page; both
implementations share the isotropic restriction, but the structural
arguments below extend without modification to anisotropic scattering
once the sources include the appropriate angular harmonics.

Integration along the characteristic
:math:`\mathbf r'(s) = \mathbf r - s\Omega` from the upstream
boundary to the receiver point, with the assumption of vacuum
boundary on the upstream end (the BC is restored in
:ref:`peierls-bc-foundations` below), gives the **integral form**:

.. math::
   :label: peierls-integral-form

   \psi(\mathbf r,\Omega) = \int_0^\infty
       \frac{Q(\mathbf r - s\Omega)}{4\pi}\,
       e^{-\Sigt{}\,s}\,\mathrm d s.

The scalar flux :math:`\phi(\mathbf r) = \int\psi(\mathbf r,\Omega)\,
\mathrm d\Omega` then satisfies

.. math::
   :label: peierls-3d

   \phi(\mathbf r) = \int Q(\mathbf r')\,
       \frac{e^{-\Sigt{}\,|\mathbf r-\mathbf r'|}}
            {4\pi |\mathbf r-\mathbf r'|^2}\,\mathrm d^3\mathbf r',

after recognising that :math:`\mathrm d^3\mathbf r' =
s^2\,\mathrm d s\,\mathrm d\Omega` and the kernel
:math:`\kappa_3 = e^{-\Sigt{} R}/(4\pi R^2)` is the *3-D point
kernel* (the Green's function of the streaming-plus-collision
operator with isotropic source). The full derivation lives in
:ref:`theory-peierls-nystrom` Section 1; what matters for this page
is that all 1-D integral equations below are *projections* of
:math:numref:`peierls-3d` along the symmetry directions of each
geometry.

Slab reduction (:math:`E_n` machinery)
--------------------------------------

For slab geometry, integrating :math:numref:`peierls-3d` over the
two transverse coordinates :math:`(y, z)` collapses the kernel to
the **first exponential integral** :math:`E_1`,

.. math::
   :label: peierls-slab-foundations

   \Sigt{}\,\phi(x) = \frac{1}{2}\int_0^L
       Q(x')\,E_1(\Sigt{}\,|x-x'|)\,\mathrm d x',

with :math:`E_n(\tau) = \int_1^\infty t^{-n}\,e^{-\tau t}\,\mathrm d t`.
The slab kernel has a logarithmic singularity at :math:`x = x'`
(:math:`E_1(\tau) \sim -\ln\tau` as :math:`\tau\to 0^+`), which is
integrable but requires either composite quadrature with diagonal
subtraction or — the path ORPHEUS takes — adaptive ``mpmath.quad``
with the singularity handled by the integrator's tanh-sinh rule.

Cylinder reduction (Bickley-Naylor :math:`\mathrm{Ki}_n` machinery)
-------------------------------------------------------------------

For an infinitely long cylinder with axisymmetric source, integrating
:math:numref:`peierls-3d` over the axial coordinate :math:`z` collapses
the kernel to the **Bickley-Naylor function**
:math:`\mathrm{Ki}_n`,

.. math::
   :label: peierls-cyl-foundations

   \Sigt{}\,\phi(r) = \int_0^R Q(r')\,\frac{r'}{|r-r'|}\,
       \mathrm{Ki}_1(\Sigt{}\,|r-r'|)\,\mathrm d r'\,\cdots
       \;+\; \text{angular details},

with :math:`\mathrm{Ki}_n(\tau) = \int_0^{\pi/2} \cos^{n-1}\theta\,
e^{-\tau/\cos\theta}\,\mathrm d\theta` (the equivalent of the slab
:math:`E_n` in cylinder geometry; see Sanchez 1982 §IV). The
cylinder kernel has an *integrable* singularity at :math:`r = r'`
(:math:`\mathrm{Ki}_1(\tau) \sim \pi/2 - 2\tau` as :math:`\tau\to
0^+` — finite limit). The full polar-form expression, including the
in-plane azimuthal angle :math:`\beta` and the chord parametrisation,
lives in :ref:`theory-peierls-nystrom` Section 2 and the derivation
notes at :ref:`theory-peierls-capabilities`.

Sphere reduction (radial-integration / cosh-extension)
-------------------------------------------------------

For sphere geometry, the radial-symmetry integration of
:math:numref:`peierls-3d` produces an integral equation in
:math:`r\phi(r)` whose kernel is the **bare exponential**
:math:`e^{-\Sigt{}\,|r-r'|}`. Two equivalent reductions appear in
the literature:

- The **Pomraning-Siewert 1982** form (:math:`r\phi`, integration
  over :math:`\mu` then half-space addition):

  .. math::
     :label: peierls-sph-ps1982-foundations

     r\phi(r) = \int_0^R x\,Q(x)\,
         \bigl[\,E_1(\Sigt{}\,|r-x|) - E_1(\Sigt{}\,(r+x))\,\bigr]\,
         \mathrm d x,

  derived in [PS1982]_ Eq. (21) (vacuum BC, isotropic source,
  homogeneous medium). The :math:`r\phi(r)` device makes the radial
  problem behave like a slab on :math:`r \in [-R, R]` with
  reflective symmetry at :math:`r = 0`.

- The **Sanchez 1986** form (cosh-even-extension):

  .. math::

     \bar g_2(\rho'\to\rho) = (\rho'/\rho)\,
         \bigl[\,E_1(\Sigt{}\,|r-r'|) - E_1(\Sigt{}\,(r+r'))\,\bigr]
         /2,

  obtained by extending the sphere problem to a slab on :math:`r \in
  [-R, R]` via the cosh-symmetric trial functions. The two forms are
  algebraically equivalent up to a factor of 2 (PS multiplies by
  :math:`r\phi(r)`, Sanchez gives the angular kernel directly).

Both forms agree at the closed-form level. The two derivations are
**structurally independent** — PS-1982's path is integrate-over-:math:`\mu`-
then-add-half-spaces; Sanchez 1986's path is cosh-even-extension. PS-1982
explicitly cites a method-of-characteristics derivation as a third
independent confirmation. This *triangulation* is what makes the
sphere Peierls kernel a robust L1 reference: any two of the three
derivations agreeing rules out reference contamination.

The polar form for sphere with bare :math:`e^{-\Sigt{}\,R}` 3-D
point kernel, used by both ORPHEUS implementations, is documented in
detail at :ref:`theory-peierls-nystrom` Section 3 (the unified polar
form) and at :ref:`theory-peierls-greens` (Variant α trajectory
geometry).


.. _peierls-bc-foundations:

Boundary conditions parametrised by α + β (Sanchez 1986 Eq. A3.a)
=================================================================

The **most general** linear boundary condition for the integral
transport equation, following Sanchez 1986 [SanchezTTSP1986]_ Eq.
(A3.a), parametrises the BC by two coefficients :math:`\alpha` and
:math:`\beta`:

.. math::
   :label: peierls-bc-general

   \psi(\mathbf r_b,\Omega) = K(\Omega) +
       \alpha\,\psi(\mathbf r_b,\Omega_R) +
       \beta\,\chi(\Omega)\!\int_{\Omega'\cdot n>0}
            \psi(\mathbf r_b,\Omega')\,(\Omega'\cdot n)\,\mathrm d\Omega',
       \qquad \Omega\cdot n \le 0,

where :math:`\Omega_R = \Omega + 2\mu\,n` is the specularly-reflected
direction, :math:`n` is the outward normal at the boundary point
:math:`\mathbf r_b`, :math:`K(\Omega)` is the externally-imposed
incident flux, and :math:`\chi(\Omega)` is a normalised diffuse
re-emission angular distribution
(:math:`\int\chi(\Omega)(\Omega\cdot n)\,\mathrm d\Omega = 1` per
[SanchezTTSP1986]_ Eq. (A3.b)).

Five canonical BC choices fall out of :math:numref:`peierls-bc-general`:

.. list-table:: Boundary conditions in the (α, β, K, χ) parametrisation
   :header-rows: 1
   :widths: 18 8 8 22 44

   * - Name
     - α
     - β
     - K, χ
     - Physical content
   * - Vacuum
     - 0
     - 0
     - K = 0
     - No reflection, no diffuse re-emission, no incident flux.
       Outgoing rays escape.
   * - White (isotropic)
     - 0
     - 1
     - χ ∝ 1
     - Diffuse isotropic re-emission of all incoming flux. The
       textbook "albedo = 1, isotropic" closure used by flat-source
       CP. Rank-1 Mark in :ref:`theory-peierls-nystrom` §8;
       Hébert 2009 §3.8.5 gives the rank-1 :math:`(1-P_{ss})^{-1}`
       scalar form.
   * - Specular
     - 1
     - 0
     - K = 0
     - Mirror reflection: angle of incidence = angle of reflection.
       Closed sphere with no leakage; rank-1 isotropic mode →
       :math:`k_{\rm eff} = k_\infty` exactly. The case Variant α
       is shipped for.
   * - Partial-albedo
     - α∈[0,1]
     - β∈[0,1−α]
     - any
     - Convex combination of vacuum / specular / diffuse. The
       :math:`\alpha`-parameter in
       :func:`~orpheus.derivations.continuous.peierls.greens_function.solve_greens_function_sphere`
       interpolates the full :math:`\alpha\in[0,1]` range.
   * - Periodic
     - n/a
     - n/a
     - n/a
     - Lattice geometries; Sanchez 2002
       :math:`\psi = \psi_q(L)/(1-\psi_{bd}(L))\cdot\psi_{bd} + \psi_q`
       (Eq. 15) closure. Not shipped in either ORPHEUS Peierls
       family.

Both ORPHEUS Peierls families implement subsets of this
parametrisation:

- **Nyström / matrix-Galerkin** (:ref:`theory-peierls-nystrom`):
  vacuum, white (rank-1 Mark, F.4 rank-2 per-face), specular
  (rank-:math:`N`, multi-bounce-corrected). All implemented as
  separate ``boundary=`` strings. The :math:`(\alpha,\beta)`
  parametrisation is *not* exposed at the public API — instead each
  closure is hard-coded as a discrete kernel-builder.
- **Green's function (Variant α)** (:ref:`theory-peierls-greens`):
  vacuum and specular as the **two endpoints of a single
  :math:`\alpha`-parametrised solver**; partial-albedo
  :math:`\alpha\in(0,1)` is reachable without a separate code path.
  This is the load-bearing structural advantage of the Green's
  function reformulation: the BC is encoded in the kernel via Sanchez
  Eq. (A1) :math:`t = \bar t + t_h`, so the closure question
  *dissolves* — there is no separate ``K_bc`` matrix, no rank-:math:`N`
  gating, no :math:`(1-P_{ss})^{-1}` scalar factor at the operator
  level.

The :math:`\beta`-branch (diffuse re-emission) is **not** shipped
in either family. Sanchez 1986 Eq. (A6) carries the full
:math:`(\alpha,\beta)` kernel, but the prototype Green's function
solver in
:func:`~orpheus.derivations.continuous.peierls.greens_function.solve_greens_function_sphere`
is restricted to :math:`\beta = 0`. Adding :math:`\beta`-support is
flagged as future work in :ref:`theory-peierls-greens`.


Two architectural choices for the discretisation
=================================================

The two implementations partition cleanly along the operator they
discretise. Use this comparison table to decide which forward-link
to follow:

.. list-table:: Nyström / matrix-Galerkin vs Green's function (Variant α)
   :header-rows: 1
   :widths: 30 35 35

   * - Property
     - :ref:`theory-peierls-nystrom` (Nyström)
     - :ref:`theory-peierls-greens` (Green's function)
   * - Operator discretised
     - Angle-integrated kernel
       :math:`g_d(\rho'\to\rho)` — assembled as matrix
       :math:`K_{ij} = w_j\,g_d(r_j\to r_i)`
     - Angle-resolved Green's function
       :math:`\tilde t(r'\to r,\mu)` — sampled along single
       characteristics; iterated as
       :math:`\psi^{(n+1)} = K[\psi^{(n)}]`
   * - Spatial state
     - 1-D :math:`\phi(r)` on a radial Nyström grid
     - 2-D :math:`\psi(r,\mu)` on a phase-space grid
   * - BC handling
     - Separate closure tensor
       :math:`K_{\rm bc} = G\cdot R\cdot P` with rank-:math:`N`
       Marshak / F.4 / specular variants
     - BC absorbed into kernel via Sanchez Eq. (A1); no separate
       :math:`K_{\rm bc}`; bounce sum closed analytically as
       :math:`T(\mu_{\rm surf}) = 1/(1 - \alpha\,e^{-\Sigt{}\,L_p})`
   * - Geometries
     - Slab + cylinder + sphere; hollow + solid topology classes
     - Sphere only (homogeneous + multi-region)
   * - Closure types
     - vacuum, white rank-1 Mark, white_f4 (rank-2 per-face),
       white_hebert, specular rank-:math:`N`, specular_multibounce
     - vacuum (:math:`\alpha = 0`), specular (:math:`\alpha = 1`),
       partial-albedo (:math:`\alpha\in(0,1)`); β-branch not shipped
   * - Multi-region
     - Slab only (rank-2 per-face); cyl/sph rank-:math:`N`-per-face
       falsified (Issue #133, see
       :ref:`peierls-rank-n-per-face-closeout`)
     - Sphere yes (Plan-(b) Option 1 fixed-source + Option 2
       k-eigenvalue); piecewise :math:`\Sigma_t` along trajectory
       and bounce-period chord
   * - Multi-group
     - Production via
       :func:`~orpheus.derivations.continuous.peierls.geometry.solve_peierls_mg`
       (Issue #104); shipped registry rows for slab + hollow
       cyl/sph 2G
     - Production via
       :func:`~orpheus.derivations.continuous.peierls.greens_function.solve_greens_function_sphere_mg`
       (closed sphere reduces to
       :func:`~orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`
       transfer-matrix dominant eigenvalue)
   * - Anisotropic scattering
     - Issue #112 / #100 open; rank-1 Mark scalar only
     - Future work; Sanchez 1986 :math:`h`-kernel
       :math:`q\to q_0 + \omega_1\,\Omega\cdot J` extends along
       the trajectory machinery
   * - Production status
     - **Production reference for all configurations** (slab + cyl
       + sphere; vacuum + white + specular)
     - **Research-grade reference** for closed sphere homogeneous
       (V_α1 exact = :math:`k_\infty`); attacks Issue #132 Class B
       MR catastrophe via piecewise-:math:`\Sigma_t` extension
   * - Where it shines
     - Vacuum / white BC at any geometry; rank-1 Mark for solid
       Class B; F.4 rank-2 for hollow Class A
     - Closed sphere homogeneous (exact), multi-region sphere
       (no rank-:math:`N` closure pathology), sphere fixed-source
       cross-checked against Garcia 2021 stable-:math:`P_N`
   * - Where it fails
     - Specular / hypersingular kernel at the diagonal — Phase 5
       retreat (:ref:`peierls-phase5-retreat`) closed Nyström
       sampling of :math:`g_\alpha`. Multi-region sphere has the
       Issue #132 mode-0 / mode-:math:`n\ge 1` normalisation
       mismatch.
     - Cylinder geometry (not implemented). Anisotropic scattering
       (not implemented). β-branch diffuse re-emission (not
       implemented). Cubic-spline source-interpolation smooths
       discontinuous :math:`\sigma_s` at multi-region interfaces
       — accounts for ~12 % near-interface error vs Garcia 2021.

The two families are **complementary**, not competing. The
Nyström family is the production reference for nearly all
verification chains in ORPHEUS; Variant α is the parallel
research-grade reference that closes the cases Phase 4 cannot
handle correctly (closed-sphere specular exact; multi-region sphere
without the mode-mixing pathology).

For the verification matrix that maps each shipped configuration to
its production reference family, see
:ref:`theory-peierls-capabilities`. Cardinal Rule 1 (correctness is
critical): never use Variant α as a reference for a case it does not
cover (cylinder, anisotropic scattering); never trust the rank-N
Marshak closure for Class B multi-region (Issue #132 documents the
+57 % catastrophe).


Common verification chain
==========================

Both families share the same V&V framework, the same literature
references, and the same per-case verification chain. The references
below are *method-agnostic*: they apply equally to the Nyström and
Green's function paths.

The k-inf transfer-matrix anchor
---------------------------------

For closed sphere with isotropic scattering and any :math:`(\alpha,\beta)`
BC that admits the constant function as an eigenmode (specular at
:math:`\alpha = 1`; white at :math:`\alpha = 0`, :math:`\beta = 1`),
the eigenvalue collapses to :math:`k_{\rm eff} = k_\infty =
\nSigf{}/\Siga{}` *independent of the discretisation*. Both ORPHEUS
families must reproduce this identity to within their respective
quadrature errors — Variant α at machine precision (V_α1 algebraic
identity), Phase 4 ``specular_multibounce`` at rank-:math:`N`
truncation error :math:`\sim 0.27\,\%/0.25\,\%/0.12\,\%` for
:math:`N\in\{1,2,3\}`. The :math:`k_\infty` reference is
:func:`~orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`
(transfer-matrix dominant eigenvalue + spectrum, multi-group).

The Pomraning-Siewert 1982 vacuum sphere
-----------------------------------------

[PS1982]_ Eq. (21) is the **structurally-independent** vacuum-sphere
reference for both families. Its derivation path
(integrate-over-:math:`\mu` then add half-spaces) is genuinely
different from the Sanchez 1986 cosh-even-extension path; PS-1982
itself confirmed via method-of-characteristics gives a third
independent confirmation of the same kernel. The local PDF copy is
``1982NSE80-481.pdf`` in the repo root (commercial title:
*JQSRT* 28(6), 503-506 (1982)).

ORPHEUS implementation:
:func:`orpheus.derivations.continuous.peierls.ps1982_reference.solve_ps1982_vacuum_sphere`
— Nyström quadrature on the
:math:`[E_1(|r-x|) - E_1(r+x)]` kernel via
:func:`mpmath.expint`. Used as the L1 cross-check truth source for
Variant α vacuum BC; also serves as a candidate L1 reference for the
Nyström family's vacuum sphere closure when high precision is needed.

The Sanchez 1986 / 2002 architecture references
------------------------------------------------

[SanchezTTSP1986]_ is the load-bearing literature for both families.
Its Appendix gives the Green's function :math:`t(r',\Omega'\to r,
\Omega)` for sphere with the full
:math:`(\alpha,\beta)`-parametrised BC; the angle-integrated reduction
:math:`g_\alpha(\rho'\to\rho) = \int t\,\mathrm d\Omega` is what the
Nyström family discretises (Phase 4 matrix-Galerkin) and what the
Green's function family *avoids* discretising (Variant α works with
:math:`t` directly). Both perspectives are correct math; the
discretisation choice is what makes them different.

[Sanchez2002]_ extends the trajectory architecture to *lattice*
geometries via the periodic-trajectory closure
:math:`\psi = \psi_q(L)/(1-\psi_{bd}(L))\cdot\psi_{bd} + \psi_q`
(Eq. 15) — algebraically parallel to Variant α's
:math:`\psi_{\rm surf} = T(\mu_{\rm surf})\,B(\mu_{\rm surf})`
bounce sum, but for a different geometry. Cross-check on the
universality of the multi-bounce factor.

Hébert 2009 and Stamm'ler 1983
-------------------------------

[Hebert2020]_ Chapter 3 is the textbook reference for the
**rank-1** white-BC closure :math:`(1-P_{ss})^{-1}` (§3.8.5;
Eq. 3.323 = :math:numref:`hebert-3-323`) and the rank-2 per-face
F.4 closure (§3.8.4). [Stamm1983]_ Chapter 4 gives the same closure
in a different notation (Stamm'ler Eq. 34 = Hébert 3.323). Variant α's
V_α2 algebraic identity :math:`T_{00}^{\rm sphere} =
P_{ss}^{\rm sphere}` (Eq. :math:numref:`peierls-greens-V-alpha-2` in
:ref:`theory-peierls-greens`) explains *why* Variant α at rank-1
agrees bit-for-bit with the existing Phase 4
``boundary="specular_multibounce"`` at :math:`N=1` and with
``boundary="white_hebert"`` rank-1: all three reduce to the same
closed-form geometric series.

[Stamm1983]_ Chapter 6 is the architectural reference for the
**multi-group** Peierls form — group-wise iteration with
cross-group + fission source. The same form is implemented in
both ORPHEUS families
(:func:`~orpheus.derivations.continuous.peierls.geometry.solve_peierls_mg`
for the Nyström family;
:func:`~orpheus.derivations.continuous.peierls.greens_function.solve_greens_function_sphere_mg`
for the Green's function family).

The Garcia 2018 / 2020 / 2021 stable :math:`P_N` family
--------------------------------------------------------

Garcia's stable-:math:`P_N` family of papers ([Garcia2020]_,
[Garcia2021]_; the 2018 JCTT predecessor that establishes the
exterior-of-sphere stability) is the only modern external numerical
reference in the published literature for sphere homogeneous +
multi-region transport with reflective BC. Garcia 2021 specifically
covers the multi-region sphere with internal sources and provides
4-significant-figure converged scalar-flux profiles — the L1
flux-shape reference for Variant α multi-region (Plan-(b) Option 1;
documented at :ref:`theory-peierls-greens`). Garcia 2021 is
*subcritical-only*: criticality is explicitly out of scope of the
2021 paper (§III.A); for k-eigenvalue cross-checks the chain extends
through :func:`~orpheus.derivations.common.eigenvalue.kinf_and_spectrum_homogeneous`
(homogeneous limit) plus PS-1982 (vacuum BC).

Per-case routing
-----------------

The verification chain (which reference for which case) is:

.. list-table:: Per-case verification chain
   :header-rows: 1
   :widths: 30 35 35

   * - Configuration
     - Reference solver
     - Truth source
   * - Closed homogeneous sphere, specular
     - Variant α (machine-exact)
     - :math:`k_\infty = \nSigf{}/\Siga{}` analytic
   * - Vacuum sphere, homogeneous
     - PS-1982 Eq. (21) Nyström
     - PS-1982 / Sanchez 1986 / MoC triangulation
   * - Multi-region sphere fixed-source
     - Variant α + cubic-spline source
     - Garcia 2021 Table 5 (Williams 1991 Case 1)
   * - Multi-region sphere k-eigenvalue
     - Variant α MR
     - Issue #132 reproducer; no published L1
       eigenvalue reference yet (Garcia 2021 subcritical-only;
       Sood 2003 covers single-region critical sphere)
   * - Slab / cyl / sph hollow rank-2 white
     - Phase 4 F.4 (Stamm'ler Eq. 34)
     - :math:`k_\infty` Wigner-Seitz identity (slab exact;
       cyl/sph carry residual scalar-mode error from
       :ref:`peierls-rank-n-per-face-closeout`)
   * - Solid Class B sphere/cylinder, specular
     - Phase 4 ``specular_multibounce`` rank-:math:`N\le 3`
       (UserWarning at :math:`N\ge 4`)
     - Phase 4 N=3 vs Variant α (sphere homogeneous;
       cylinder no Variant α reference)


Provenance: local PDFs and codebase pointers
=============================================

The literature PDFs cited above are stored in the repo root
``/scratch/literature/`` (and a few at the top level):

- ``1982NSE80-481.pdf`` — Pomraning-Siewert 1982 (vacuum sphere
  integral form). Local copy.
- ``Hebert(2009)Chapter3.pdf`` — Hébert 2009 Ch. 3 (CP, integral
  transport, white-BC closures). Local copy.
- ``Sanchez(2002).pdf`` — Sanchez 2002 (periodic-trajectory closure,
  lattice). Local copy.
- ``Stammler(1983)Chapter4.pdf`` + ``Stammler(1983)Chapter6.pdf`` —
  Stamm'ler 1983 (rank-2 per-face closure; multi-group). Local
  copies.
- ``Stacey(2007)Chapter9.pdf`` — Stacey 2007 (textbook intro to
  integral transport). Local copy.
- ``Ligou(1982)Chapter8.pdf`` — Ligou 1982 (collision probability,
  multi-region). Local copy.
- ``Hebert(2009)Chapter3.pdf`` — Hébert 2009 Ch. 3. Local copy.
- ``Approximate Solutions of the Two-Dimensional Integral Transport
  Equation by Collision Probability Methods.pdf`` — Sanchez 1982
  (anisotropic 2D CP). Local copy.
- ``Collision Probabilities for Finite Cylinders and Cuboids.pdf`` —
  finite-axial-extent cylinder + cuboid CP. Local copy.
- ``Integral form of the equation of transfer for a homogeneous
  sphere with linearly anisotropic scattering.pdf`` — Sanchez 1986
  precursor (linearly anisotropic scattering). Local copy.
- Sanchez 1986 [SanchezTTSP1986]_ — primary; PDF in repo root as
  ``Sanchez(1986)TTSP14.pdf``.
- Garcia 2018 / 2020 / 2021 [Garcia2020]_ [Garcia2021]_ — paywalled;
  no local copy.
- Mitsis 1963 (ANL-6787) — pseudo-slab equivalence for vacuum sphere;
  numerical reference. No local copy.
- Sood-Forster-Parsons 2003 (LA-13511) — analytical benchmark
  test set (critical sphere c_crit(R) tables). No local copy.

Codebase pointers:

- :mod:`orpheus.derivations.continuous.peierls.geometry` — Nyström
  family unified solver
  (:func:`~orpheus.derivations.continuous.peierls.geometry.solve_peierls_1g`,
  :func:`~orpheus.derivations.continuous.peierls.geometry.solve_peierls_mg`).
- :mod:`orpheus.derivations.continuous.peierls.greens_function` —
  Green's function family
  (:func:`~orpheus.derivations.continuous.peierls.greens_function.solve_greens_function_sphere`,
  ``_mg``, ``_mr``, ``_mr_fixed_source``).
- :mod:`orpheus.derivations.continuous.peierls.ps1982_reference` —
  PS-1982 reference solver (vacuum sphere only).
- :mod:`orpheus.derivations.continuous.peierls.origins.specular.greens_function`
  — SymPy derivations V_α1, V_α2, V_α3 for the Green's function
  family.
- :file:`docs/theory/_peierls_capability_matrix.inc.rst` —
  auto-generated from
  :func:`orpheus.derivations.continuous.peierls.cases.capability_rows`;
  rebuilt every Sphinx build by
  :mod:`tools.verification.generate_peierls_matrix`.


Provenance: literature references
==================================

.. [PS1982] G.C. Pomraning and C.E. Siewert, "On the integral form of
   the equation of transfer for a homogeneous sphere," *J. Quant.
   Spec. Rad. Transfer* **28**, 503–506 (1982).
   DOI: 10.1016/0022-4073(82)90016-4. The vacuum-sphere reduction
   used by both ORPHEUS Peierls families as a structurally-independent
   L1 cross-check.
