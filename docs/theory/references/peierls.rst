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
- :ref:`theory-trajectory-resolvent` — the research-grade **Green's
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
at the diagonal and no :term:`quadrature` trick rescues it. The
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
discrete-ordinates / method-of-characteristics (MoC) :term:`sweeps <sweep>` used
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
   ``orpheus.derivations.continuous.peierls_nystrom`` is the L1 reference
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
   performed analytically. There is no :term:`discrete-ordinates <ordinate>` ray-effect
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
emission transport equation for the :term:`angular flux`
:math:`\psi(\mathbf r,\Omega)`,

.. math::
   :label: peierls-boltzmann

   \Omega\cdot\nabla\psi(\mathbf r,\Omega) + \Sigt{}\,
       \psi(\mathbf r,\Omega) = \frac{1}{4\pi}\,Q(\mathbf r),

.. (vv-status rationale) governing: Governing equation defining the problem class — the entire solver is the verification. Governing Boltzmann transport equation for the steady-state monoenergetic angular flux — the entire transport solver IS the verification (mirror of the documented `boltzmann` label).
.. vv-status: peierls-boltzmann documented


where :math:`Q(\mathbf r) = \Sigs{}\,\phi(\mathbf r) +
(\nSigf{}/k)\,\phi(\mathbf r)` is the isotropic emission density
(scattering + fission, k-eigenvalue normalisation), and
:math:`\phi(\mathbf r) = \int\psi(\mathbf r,\Omega)\,\mathrm d\Omega`
is the :term:`scalar flux`. Anisotropic scattering and the
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

.. (vv-status rationale) derivation: Intermediate step in a derivation chain whose capstone identity is tested elsewhere. Integral form of the transport equation along characteristics — derivation step toward the volume-integrated peierls forms; all downstream Peierls reductions are verified at L1 by peierls-unified.
.. vv-status: peierls-integral-form documented


The scalar flux :math:`\phi(\mathbf r) = \int\psi(\mathbf r,\Omega)\,
\mathrm d\Omega` then satisfies

.. math::
   :label: peierls-3d

   \phi(\mathbf r) = \int Q(\mathbf r')\,
       \frac{e^{-\Sigt{}\,|\mathbf r-\mathbf r'|}}
            {4\pi |\mathbf r-\mathbf r'|^2}\,\mathrm d^3\mathbf r',

.. (vv-status rationale) derivation: Intermediate step in a derivation chain whose capstone identity is tested elsewhere. 3D point-kernel form of the volume-integrated transport equation — definitional reduction of the integral form.
.. vv-status: peierls-3d documented


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

.. (vv-status rationale) derivation: Intermediate step in a derivation chain whose capstone identity is tested elsewhere. Slab Peierls reduction Σ_t φ(x) = (1/2) ∫ Q E_1 dx' — the specialised kernel is verified by TestSlabKernelRowSum at L1.
.. vv-status: peierls-slab-foundations documented


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

.. (vv-status rationale) derivation: Intermediate step in a derivation chain whose capstone identity is tested elsewhere. Cylinder Bickley-Naylor reduction with Ki_n machinery — the specialised kernel is verified by TestCylinderKernelRowSum at L1.
.. vv-status: peierls-cyl-foundations documented


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

  .. (vv-status rationale) derivation: Intermediate step in a derivation chain whose capstone identity is tested elsewhere. PS-1982 sphere reduction r·φ(r) = ∫ x Q (E_1(|r-x|) - E_1(r+x)) — the specialised kernel is verified by TestSphereKernelRowSum at L1 and the PS-1982 cross-check.
  .. vv-status: peierls-sph-ps1982-foundations documented


  derived in :cite:`PomraningSiewert1982` Eq. (21) (vacuum BC, isotropic source,
  homogeneous medium). The :math:`r\phi(r)` device makes the radial
  problem behave like a slab on :math:`r \in [-R, R]` with
  reflective symmetry at :math:`r = 0`.

- The **Sanchez 1986** form (cosh-even-extension):

  .. math::
     :label: peierls-sph-sanchez-foundations

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
form) and at :ref:`theory-trajectory-resolvent` (Variant α trajectory
geometry).


.. _peierls-bc-foundations:

Boundary conditions parametrised by α + β (Sanchez 1986 Eq. A3.a)
=================================================================

The **most general** linear boundary condition for the integral
transport equation, following Sanchez 1986 :cite:`SanchezTTSP1986` Eq.
(A3.a), parametrises the BC by two coefficients :math:`\alpha` and
:math:`\beta`:

.. math::
   :label: peierls-bc-general

   \psi(\mathbf r_b,\Omega) = K(\Omega) +
       \alpha\,\psi(\mathbf r_b,\Omega_R) +
       \beta\,\chi(\Omega)\!\int_{\Omega'\cdot n>0}
            \psi(\mathbf r_b,\Omega')\,(\Omega'\cdot n)\,\mathrm d\Omega',
       \qquad \Omega\cdot n \le 0,

.. (vv-status rationale) transcription: Verbatim transcription of a literature reference. Sanchez 1986 (A3.a) general BC parametrisation (α + β + K) — transcription of the literature definition; no separate identity to verify (specific BCs are tested individually).
.. vv-status: peierls-bc-general documented


where :math:`\Omega_R = \Omega + 2\mu\,n` is the specularly-reflected
direction, :math:`n` is the outward normal at the boundary point
:math:`\mathbf r_b`, :math:`K(\Omega)` is the externally-imposed
incident flux, and :math:`\chi(\Omega)` is a normalised diffuse
re-emission angular distribution
(:math:`\int\chi(\Omega)(\Omega\cdot n)\,\mathrm d\Omega = 1` per
:cite:`SanchezTTSP1986` Eq. (A3.b)).

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
       :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere`
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
- **Green's function (Variant α)** (:ref:`theory-trajectory-resolvent`):
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
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere`
is restricted to :math:`\beta = 0`. Adding :math:`\beta`-support is
flagged as future work in :ref:`theory-trajectory-resolvent`.


Peierls integral equation reference
====================================

The collision probability method discretises the integral transport
equation by assuming flat sources within each region.  The **Peierls
integral equation** is the continuous (pre-discretisation) form of that
equation, and solving it at high quadrature resolution provides an
independent reference for verifying the CP solver's flux profiles.

Derivation from the 1-D transport equation
-------------------------------------------

The starting point is the steady-state, isotropic-source transport
equation for a 1-D slab in energy group :math:`g`:

.. math::

   \mu\,\frac{\partial\psi_g}{\partial x}(x,\mu)
   + \Sigt{g}(x)\,\psi_g(x,\mu)
   = \frac{q_g(x)}{2}

where :math:`\psi_g(x,\mu)` is the angular flux,
:math:`\mu = \cos\theta \in [-1, 1]` is the direction cosine, and
:math:`q_g(x)` is the isotropic source (scattering + fission).

For :math:`\mu > 0`, the formal solution (integrating factor) for a
slab :math:`[0, L]` with vacuum boundary at :math:`x = 0` gives:

.. math::

   \psi_g(x,\mu) = \frac{1}{2\mu}\int_0^x
   \exp\!\Bigl(-\frac{\tau_g(x',x)}{\mu}\Bigr)\,q_g(x')\,\mathrm{d}x'

where
:math:`\tau_g(x',x) = \int_{x'}^{x}\Sigt{g}(s)\,\mathrm{d}s` is the
optical path.  A symmetric expression holds for :math:`\mu < 0`.

The scalar flux is obtained by integrating over angle:

.. math::

   \varphi_g(x) = \int_{-1}^{1}\psi_g(x,\mu)\,\mathrm{d}\mu
   = \frac{1}{2}\int_0^L q_g(x')\!\int_0^1
   \frac{1}{\mu}\exp\!\Bigl(-\frac{\tau_g(x,x')}{\mu}\Bigr)\,\mathrm{d}\mu
   \;\mathrm{d}x'

The inner angular integral is exactly the **first exponential integral**:

.. math::

   E_1(z) = \int_0^1 \frac{1}{\mu}\,e^{-z/\mu}\,\mathrm{d}\mu

(Case & Zweifel 1967, Ch. 2, Eq. 2.67; equivalently,
:math:`E_1(z) = \int_1^\infty t^{-1} e^{-zt}\,\mathrm{d}t`.)

Substituting yields the **Peierls integral equation** for the scalar flux:

.. math::
   :label: peierls-equation

   \varphi(x)
   = \tfrac12 \int_0^L E_1\!\bigl(\tau(x,x')\bigr)\,q(x')\,\mathrm{d}x'
     + \varphi_{\rm bc}(x)

where :math:`q(x') = \sum_{g'}\Sigma_{s,g'\to g}\varphi_{g'}(x')
+ \chi_g\sum_{g'}\nu\Sigma_{f,g'}\varphi_{g'}(x')/k` is the total
isotropic source and the optical path is
:math:`\tau(x,x') = \int_{\min(x,x')}^{\max(x,x')}\Sigma_t(s)\,\mathrm{d}s`.

.. note::

   The equation above is written for :math:`\varphi` (scalar flux), NOT
   :math:`\Sigt{}\varphi` (collision rate).  This matters for the
   operator form: the LHS is :math:`\varphi`, not
   :math:`\Sigt{}\varphi`, so the identity operator appears in the
   eigenvalue system :math:`(I - K_s)\varphi = (1/k)\,K_f\varphi`.
   An earlier version incorrectly placed :math:`\Sigt{}` on the LHS,
   producing a ``diag(Sigma_t) - K_scatter`` operator that broke the
   conservation property (see :ref:`peierls-conservation` below).

Why :math:`E_1` and not :math:`E_3`
------------------------------------

The CP method uses the :math:`E_3` function (double antiderivative of
:math:`E_1`) to integrate the kernel analytically over flat-source
regions (:eq:`dd-slab`, :eq:`dc-slab`, :eq:`self-slab`).  A Peierls
reference that uses :math:`E_1` directly provides **genuinely
independent** verification:

- **Different special function**: :math:`E_3` is computed from :math:`E_1`
  via recursion (:math:`E_{n+1}(z) = [e^{-z} - z\,E_n(z)]/n`).  A
  systematic error in the :math:`E_3` recursion (wrong sign, off-by-one
  in the index) would affect both the CP matrices and any reference
  built from the same recursion.

- **No flat-source assumption**: the Peierls Nystrom solver resolves the
  flux within each region as a high-order polynomial (degree
  :math:`p-1` per GL panel), while CP assumes piecewise-constant flux.
  This makes the Peierls reference sensitive to flux shape errors that
  are invisible to a CP-vs-CP comparison.

- **Common-mode failure prevention**: the existing CP eigenvalue tests
  compare the CP solver against analytically-constructed CP matrices
  (same :math:`E_3` kernel, same flat-source discretisation).  A sign
  error in the :math:`E_3` second-difference formula would cancel
  between the test and the code.  The Peierls reference would catch it.

Nystrom method details
----------------------

The Nystrom method converts the integral equation into a matrix
eigenvalue problem by approximating the integral with a quadrature rule
:cite:`Kress2014`.  For the Peierls equation, the key challenge is the
logarithmic singularity in :math:`E_1` at :math:`x = x'`.

**Composite Gauss--Legendre panels.**
The slab :math:`[0, L]` is divided into :math:`n_{\rm panels}` panels
per material region.  Within each panel, :math:`p` Gauss--Legendre
nodes and weights are computed.  The total number of quadrature nodes is
:math:`N = n_{\rm regions} \times n_{\rm panels} \times p`.

**Singularity treatment.**
The :math:`E_1` kernel has a logarithmic singularity at :math:`x = x'`:

.. math::
   :label: e1-decomposition

   E_1(z) = \bigl[-\ln z - \gamma\bigr] + R(z),
   \qquad R(z) \equiv E_1(z) + \ln z + \gamma,\quad R(0) = 0.

The remainder :math:`R(z)` is a smooth (analytic) function that
vanishes at the origin. This decomposition motivates the classical
**singularity-subtraction** approach (used in the original
implementation), in which diagonal panels split the kernel into a
smooth :math:`R` part (handled by GL weights) and a
:math:`-\ln|x_i - x'|` part (handled by product-integration weights).
Off-diagonal panels would then use standard GL weights on the smooth
:math:`E_1` integrand.

**Why the classical singularity-subtraction scheme failed.** Issue
#113 closed the slab K-matrix cross-panel logarithmic-singularity
bug — the naive GL-collocation scheme made two assumptions that
turned out to be insufficient (off-diagonal GL exact for smooth
integrands; :math:`R(\tau)` smooth across the diagonal panel) and
hid the resulting ~1 % error inside row-sum conservation because
:math:`\sum_j L_j(x') = 1` cancels basis-individual :math:`K[i,j]`
errors. The forensic — including the basis-aware ERR-027 /
ERR-028 catalog entries and the rank-N "passes at mode 0, fails
at mode n > 0" diagnostic signature — is preserved in
`Issue #113 <https://github.com/deOliveira-R/ORPHEUS/issues/113>`_.
The unified ``_pg.solve_peierls_*`` adaptive-quadrature K-matrix
assembly described next is the production path.

**Unified basis-aware assembly** (current implementation,
:func:`~orpheus.derivations.continuous.peierls_nystrom.slab._basis_kernel_weights`).
Every :math:`K[i, j]` is computed directly as

.. math::
   :label: peierls-slab-nystrom

   K[i, j] \;=\; \tfrac{1}{2} \int_{p_a}^{p_b}
                  E_1\!\bigl(\tau(x_i, x')\bigr)\,L_j(x')\,\mathrm{d}x'

via adaptive :func:`mpmath.quad`, with the subdivision hint
:math:`[p_a, x_i, p_b]` when :math:`x_i` lies inside the source
panel (same-panel case). This single code path:

- handles the integrable log singularity at :math:`x'=x_i` via the
  subdivision hint (mpmath resolves it to machine precision);
- handles the derivative kink of :math:`R(\tau)` in :math:`x'` — GL
  on each smooth half of the subdivided panel converges spectrally;
- eliminates the off-diagonal-panel quadrature error — adaptive
  refinement resolves :math:`E_1`'s non-polynomial structure on
  arbitrary panel pairs without relying on near-log assumptions.

The implementation exactly mirrors the adaptive reference
:func:`~orpheus.derivations.continuous.peierls_nystrom.reference.slab_K_vol_element`, so
the production code and the reference agree to machine
:math:`\mathrm{dps}` by construction.

**Why adaptive quadrature over the alternatives.**
Four strategies exist for the Peierls log singularity :cite:`Atkinson1997`:

1. *Graded meshes* (cluster GL nodes near the diagonal): algebraic
   convergence only; many nodes needed for high accuracy.
2. *IMT transformation* (Iri–Moriguti–Takahashi double-exponential):
   effective for endpoint singularities but hard to control for the
   travelling singularity at :math:`x'=x_i`.
3. *Singularity subtraction + product integration*: the classical
   approach. Exact for the isolated log, but still requires the
   surrounding integrand to be polynomial-representable — which is
   false for the full :math:`E_1`-times-:math:`L_j` product (see
   bugs above).
4. *Adaptive :func:`mpmath.quad` with explicit subdivision hints*
   (used here): handles all three non-smoothness sources — log
   singularity at :math:`x'=x_i`, derivative kink of :math:`R` at
   :math:`x'=x_i`, and non-polynomial :math:`E_1` decay across the
   source panel — in one uniform code path. The cost is higher per
   weight than GL product integration, but this is a reference-quality
   solver; correctness beats speed.

White boundary conditions
-------------------------

For an infinite lattice (white BC), the re-entrant isotropic flux is
derived by requiring that the angular flux entering each face equals
the isotropic average of the flux exiting the opposite face.  The
result is a separable rank-2 perturbation of the interior kernel:

.. math::
   :label: peierls-white-bc

   \varphi_{\rm bc}(x) = \frac{1}{1 - T^2}\sum_j w_j\,q_j\,
   \bigl[E_2(\tau_{0,i})(E_2(\tau_{0,j}) + T\,E_2(\tau_{L,j}))
   + E_2(\tau_{L,i})(E_2(\tau_{L,j}) + T\,E_2(\tau_{0,j}))\bigr]

where :math:`T = 2\,E_3(\tau_{\rm total})` is the slab transmission
probability,
:math:`\tau_{0,i} = \int_0^{x_i}\Sigt{}(s)\,\mathrm{d}s`, and
:math:`\tau_{L,i} = \int_{x_i}^{L}\Sigt{}(s)\,\mathrm{d}s`.

The :math:`E_2` factors represent the uncollided angular flux
arriving at point :math:`x_i` from a face source:
:math:`E_2(\tau) = \int_0^1 e^{-\tau/\mu}\,\mathrm{d}\mu`.  The
:math:`1/(1-T^2)` denominator sums the geometric series of
multiple slab traversals.

.. _peierls-conservation:

Conservation property
---------------------

For the :math:`\varphi` equation (identity on the LHS, not
:math:`\Sigt{}`), the Nystrom kernel satisfies the following
conservation identity for constant unit flux with white BC:

.. math::
   :label: peierls-slab-row-sum-identity

   \sum_{j=1}^{N} K_{\rm total}[i, j]\,w_j = \frac{1}{\Sigt{i}}
   \qquad\text{(row-sum identity)}

where :math:`K_{\rm total}` includes both the interior :math:`E_1`
kernel and the white-BC re-entry kernel.

**Physical interpretation**: a spatially uniform, isotropically
scattering medium with :math:`\Sigs{} = \Sigt{}` (no absorption)
must produce :math:`\varphi(x) = \text{const}` everywhere.
Substituting :math:`q(x') = \Sigt{}\,\varphi_0` into
:eq:`peierls-equation` and requiring
:math:`\varphi(x_i) = \varphi_0` gives the row-sum identity above.

.. warning::

   An earlier version of the Peierls solver used
   :math:`\Sigt{}\,\varphi(x)` on the LHS instead of :math:`\varphi(x)`.
   This corresponds to a different normalisation where the eigenvalue
   system is :math:`(\mathrm{diag}(\Sigt{}) - K_s)\varphi
   = (1/k)\,K_f\varphi`.  The row-sum identity then becomes
   :math:`\sum_j K[i,j]\,w_j = 1`, which ALSO holds --- but the
   operator itself is wrong because the :math:`\Sigt{}` factor changes
   the eigenvectors.  The identity-LHS form is the physically correct
   normalisation for computing the scalar flux.  This was a debugging
   insight during Phase 4.1 development.

Relationship to CP
------------------

The CP flat-source approximation integrates the :math:`E_1` kernel
analytically over each region to obtain the :math:`E_3` second-difference
formulae (:eq:`dd-slab`, :eq:`dc-slab`, :eq:`self-slab`).  The
Peierls reference bypasses this integration and solves the full
integral equation numerically, making it an independent check on
the CP method's spatial discretisation.

Multi-group scattering convention
---------------------------------

See :ref:`peierls-scattering-convention` (in the Peierls unified
theory page) for the canonical, project-wide statement of the
``sig_s[g_src, g_dst]`` convention, which the CP / Peierls
drivers and the XS library
(:mod:`orpheus.derivations.common.xs_library`) all follow. In the Peierls
slab assembly loop
(:func:`orpheus.derivations.continuous.peierls_nystrom.slab._build_system_matrices`)
the scatter kernel for the equation in group ``ge`` sums over
source groups ``gs`` via ``sig_s_at_node[j][gs][ge]`` =
:math:`\Sigma_{s,\,gs \to ge}` — first index source, second
destination.

.. warning::

   A naive reading of ``sig_s[gs][ge]`` might suggest "from ge to gs"
   (confusing row/column with from/to).  The correct convention is
   ``sig_s[from, to]``, consistent with the ``Mixture.SigS`` storage
   where ``SigS[g_from, g_to]`` and the in-scatter source is
   ``Q = SigS^T @ phi``.  See the scattering convention note in the
   :ref:`key-facts-cp` section and the authoritative
   :ref:`peierls-scattering-convention` note.

Numerical evidence
------------------

The Nystrom eigenvalue converges to the CP eigenvalue as the quadrature
is refined.  The following table shows the 2-group, 2-region slab case
(materials A + B, white BC, :math:`p`-point GL per panel):

.. list-table:: Nystrom k-eigenvalue convergence (2G, 2-region slab)
   :header-rows: 1
   :widths: 15 15 25 25

   * - Panels/region
     - :math:`p`
     - Nystrom :math:`\keff`
     - Relative error vs CP
   * - 4
     - 4
     - 1.2149
     - :math:`1.4 \times 10^{-2}`
   * - 8
     - 6
     - 1.2281
     - :math:`3.3 \times 10^{-3}`
   * - 8
     - 8
     - 1.2307
     - :math:`1.2 \times 10^{-3}`

The convergence is slow because the :math:`E_1` kernel's logarithmic
singularity limits the global convergence rate even with
product-integration handling on the diagonal panel.  For verification
purposes, the 2% agreement at moderate resolution is sufficient to
confirm that both methods solve the same integral transport equation.

The registered reference in the test suite
(``continuous_cases()`` in :mod:`orpheus.derivations.continuous.peierls_nystrom.slab`) uses a
lightweight configuration (4 panels |times| 4 GL points per region,
20-digit ``mpmath`` precision) to keep import time fast.  Tests that
need higher accuracy should call
``_build_peierls_slab_case()``
directly with larger parameters.

Performance note
----------------

The ``mpmath`` solver is CPU-bound on the :math:`O(N^2)` kernel
assembly and :math:`O(N^3)` LU factorisation at arbitrary precision:

- **4 panels |times| 4 points** per region (2G 2R): ~6 s
- **8 panels |times| 6 points** per region (2G 2R): ~120 s
- **16 panels |times| 6 points** per region (2G 2R): ~10+ min

For this reason, the slow 2G 2-region convergence test is marked
``@pytest.mark.slow`` and excluded from routine CI.

.. seealso::

   :mod:`orpheus.derivations.continuous.peierls_nystrom.slab` — Nystrom solver implementation.

   :class:`orpheus.derivations.continuous.peierls_nystrom.slab.PeierlsSlabSolution` — result container
   with barycentric interpolation for flux evaluation at arbitrary points.

   ``tests/derivations/test_peierls_convergence.py`` — L0 self-convergence
   and eigenvalue agreement tests.

   ``tests/cp/test_peierls_flux.py`` — L1 CP flux convergence against
   the Peierls reference.


Peierls integral equation reference — cylinder
===============================================

The slab Peierls reference above verifies the CP flat-source
discretisation in Cartesian 1-D. The **cylindrical** Peierls reference
serves the same role for ``cyl1D`` meshes: it solves the integral
transport equation on a bare or concentric-annulus cylinder at
arbitrary quadrature order, providing an independent numerical
reference against :func:`~orpheus.cp.solver.solve_cp` and the
analytical CP eigenvalue in :mod:`orpheus.derivations.continuous.flat_source_cp.cylinder`.

Unlike the slab, the cylinder's kernel is not an exponential integral
:math:`E_n` but the **Bickley--Naylor function** :math:`\mathrm{Ki}_1`
that arises from integrating the 3-D point-kernel over the infinite
axial direction. Unlike the slab, the boundary is a continuous lateral
surface, not a pair of discrete faces, so white-BC closure is not
a rank-2 outer product. And unlike the slab, a direct Nyström
discretisation of the canonical Sanchez--McCormick chord form picks
up a **non-integrable coincident singularity** that the slab's
product-integration trick does not cure — this motivates the
reformulation described below.

The implementation lives in :mod:`orpheus.derivations.continuous.peierls_nystrom.cylinder`.
This section documents the mathematics, the formulation choice
(including the dead-end that was tried first), and the verification
evidence.

Canonical chord form and why it is not what the code solves
------------------------------------------------------------

Integrating the 3-D point kernel :math:`e^{-\tau R}/(4\pi R^{2})`
over the infinite axial coordinate :math:`z` yields the 2-D
transverse Green's function

.. math::
   :label: peierls-cylinder-green-2d

   G_{\rm 2D}(|\mathbf{r} - \mathbf{r}'|)
     \;=\; \frac{\mathrm{Ki}_1(\tau)}{2\pi\,|\mathbf{r}-\mathbf{r}'|},
   \qquad
   \tau \;=\; \int_{\mathbf{r}'}^{\mathbf{r}} \Sigma_t(\mathbf{s})\,\mathrm{d}\ell.

.. vv-status: peierls-cylinder-green-2d documented

The **pointwise** scalar-flux form of the Peierls integral equation
on a bare cylinder of radius :math:`R` is therefore

.. math::

   \varphi(\mathbf{r})
     \;=\; \frac{1}{2\pi}\!\iint_{\rm disc}
       \frac{\mathrm{Ki}_1\!\bigl(\tau(\mathbf{r},\mathbf{r}')\bigr)}
            {|\mathbf{r}-\mathbf{r}'|}\,q(\mathbf{r}')\,\mathrm{d}^{2}r'
     \;+\; \varphi_{\rm bc}(\mathbf{r}).

The classical textbook presentation (:cite:`Sanchez1982` §IV.A,
Eqs. 47--49; :cite:`Stamm1983` §6.2--6.3; :cite:`Hebert2020` Eqs. 3.95--3.110)
rotates the 2-D integral to the **chord** coordinate
system :math:`(y, r')`, where :math:`y` is the perpendicular distance
from the cylinder axis to the straight-line trajectory through
:math:`\mathbf{r}` and :math:`\mathbf{r}'`. Expressing
:math:`\mathrm{d}^{2}r'` as :math:`\bigl(r'/\sqrt{r'^{2}-y^{2}}\bigr)\,
\mathrm{d}r'\,\mathrm{d}y` on each branch and pairing the two
branches gives

.. math::
   :label: peierls-cylinder-chord-form

   \Sigma_t(r)\,\varphi(r)
     \;=\; \frac{1}{\pi}
       \int_{0}^{\min(r,R)}\!\mathrm{d}y
       \int_{y}^{R}
         \bigl[\mathrm{Ki}_1(\tau^{+}) + \mathrm{Ki}_1(\tau^{-})\bigr]\,
         \frac{q(r')\,r'}{\sqrt{r'^{2}-y^{2}}}\,\mathrm{d}r'
     \;+\; S_{\rm bc}(r).

.. warning::

   A derivation shortcut that keeps only the :math:`r'` Jacobian
   :math:`r'/\sqrt{r'^{2}-y^{2}}` (as the Phase-4.2 literature
   sweep initially reported) **is missing a factor**. Computing the
   branch sum :math:`|\mathrm{d}\alpha_{+}/\mathrm{d}y| +
   |\mathrm{d}\alpha_{-}/\mathrm{d}y|` for the
   :math:`(r', \alpha) \to (y, r')` transformation — where
   :math:`\alpha` is the chord-angle coordinate at the source point
   — gives :math:`2/\sqrt{\min(r,r')^{2} - y^{2}}`, **not**
   :math:`2` as the one-sided Jacobian would suggest. The correct
   combined Jacobian is therefore

   .. math::

      \frac{1}{\sqrt{(r^{2}-y^{2})(r'^{2}-y^{2})}},

   with a **second** integrable singularity at :math:`y = r`
   (co-located with the :math:`y = r'` root of the :math:`r'`-side
   factor when :math:`r = r'`). The unchecked "simplified" form with
   only :math:`r'/\sqrt{r'^{2}-y^{2}}` would amount to a mass-loss
   bug of the same flavour as missing the :math:`\Delta A/w_m`
   redistribution factor in a cylindrical :math:`S_N` sweep: the
   integrand does not reproduce the infinite-medium identity
   :math:`\sum_j K_{ij}\,\Sigma_t(r_j) = \Sigma_t(r_i)`.

The chord form above therefore has **two coincident endpoint
singularities** (at :math:`y = r` *and* :math:`y = r'`). The slab's
product-integration recipe absorbs one singularity (log at
:math:`x = x'`) analytically against a Lagrange basis; it does not
generalise to two coincident inverse-square-root singularities
sitting at nested quadrature endpoints. Attempting it produced
numerical divergence of the row-sum identity under refinement —
the kernel matrix simply does not converge for a
moderate-precision radial grid.

Formulation pivot: polar coordinates centred at the observer
-------------------------------------------------------------

Rather than patching the chord form, the implementation uses the
**equivalent polar form** centred on the observer. Let
:math:`\beta \in [0, 2\pi]` be the azimuth from the outward radial
direction at :math:`\mathbf{r}`, and :math:`\rho \ge 0` the distance
along the ray at angle :math:`\beta`. The source position is

.. math::
   :label: peierls-cylinder-r-prime

   r'(r, \rho, \beta) \;=\; \sqrt{r^{2} + 2r\rho\cos\beta + \rho^{2}}.

.. vv-status: peierls-cylinder-r-prime documented

Because the 2-D area element is :math:`\rho\,\mathrm{d}\rho\,
\mathrm{d}\beta` and the Green's function carries
:math:`1/|\mathbf{r} - \mathbf{r}'| = 1/\rho`, the :math:`\rho`
factor **cancels** and the integrand becomes smooth:

.. math::
   :label: peierls-cylinder-polar

   \varphi(r)
     \;=\; \frac{1}{\pi}\!
       \int_{0}^{\pi}\!\mathrm{d}\beta\!
       \int_{0}^{\rho_{\max}(r,\beta)}\!\!
         \mathrm{Ki}_1\!\bigl(\tau(r, \rho, \beta)\bigr)\,
         q\bigl(r'(r, \rho, \beta)\bigr)\,\mathrm{d}\rho
     \;+\; \varphi_{\rm bc}(r).

.. vv-status: peierls-cylinder-polar documented

The prefactor :math:`1/\pi` absorbs the :math:`1/(2\pi)` of the 2-D
Green's function plus a factor of 2 from :math:`\beta \to -\beta`
symmetry that folds :math:`[0, 2\pi] \to [0, \pi]`. The upper
radial limit is the intersection of the ray with the cylinder
boundary,

.. math::
   :label: peierls-cylinder-rho-max

   \rho_{\max}(r, \beta)
     \;=\; -r\cos\beta
         + \sqrt{r^{2}\cos^{2}\beta + R^{2} - r^{2}}.

.. vv-status: peierls-cylinder-rho-max documented

Writing the identity-LHS form used by the eigenvalue driver, the
canonical cylindrical Peierls equation solved by this module is

.. math::
   :label: peierls-cylinder-equation

   \Sigma_t(r_i)\,\varphi(r_i)
     \;=\; \frac{\Sigma_t(r_i)}{\pi}\!
       \int_{0}^{\pi}\!\mathrm{d}\beta\!
       \int_{0}^{\rho_{\max}(r_i,\beta)}\!\!
         \mathrm{Ki}_1\!\bigl(\tau(r_i, \rho, \beta)\bigr)\,
         q\!\bigl(r'(r_i, \rho, \beta)\bigr)\,\mathrm{d}\rho
     \;+\; S_{\rm bc}(r_i).

.. vv-status: peierls-cylinder-equation documented

The Sanchez tie-point and row-sum-identity tests currently carry
``@pytest.mark.verifies("peierls-equation", ...)`` (the slab label).
Retrofitting those decorators to point at the cylinder/sphere-specific
labels is tracked under `Issue #142
<https://github.com/deOliveira-R/ORPHEUS/issues/142>`_; until then
this equation is marked ``documented`` rather than ``tested`` to keep
the orphan gate honest.

.. note::

   The polar formulation is **mathematically equivalent** to the
   Sanchez chord form — the underlying integral equation is the
   same — but the :math:`\rho\,\mathrm{d}\rho\,\mathrm{d}\beta`
   area element absorbs the :math:`1/\rho` of the Green's function
   and thus the integrand
   :math:`\mathrm{Ki}_1(\tau(\rho))\,q(r'(\rho))` is **regular on
   the whole integration domain**. Ordinary tensor-product
   Gauss--Legendre quadrature converges spectrally. This is the
   dominant motivation for the pivot: the chord form trades one
   integrable singularity for two, the polar form eliminates them.

Why :math:`\mathrm{Ki}_1` and not :math:`\mathrm{Ki}_3`
-------------------------------------------------------

The flat-source CP method (:mod:`orpheus.derivations.continuous.flat_source_cp.cylinder`)
uses the :math:`\mathrm{Ki}_3` kernel because it averages the
pointwise :math:`\mathrm{Ki}_1` kernel twice — once over the
source region :math:`j` and once over the target region :math:`i`
— producing the second-difference formula

.. math::

   P_{ij} \;\propto\; \mathrm{Ki}_3(\text{gap})
                    - \mathrm{Ki}_3(\text{gap} + \tau_i)
                    - \mathrm{Ki}_3(\text{gap} + \tau_j)
                    + \mathrm{Ki}_3(\text{gap} + \tau_i + \tau_j).

See :eq:`ki3-def` for the :math:`\mathrm{Ki}_n` definition and
:eq:`chord-length` for the chord geometry underlying
:math:`P_{ij}`.

The Peierls reference solves for the **pointwise** flux
:math:`\varphi(r_i)`, not the region-average collision rate, so it
uses :math:`\mathrm{Ki}_1` directly. This is the independent-kernel
property that makes the reference useful:

- The :math:`\mathrm{Ki}_n` kernels are computed by recursion
  :math:`\mathrm{Ki}_{n+1}(x) = \int_x^\infty \mathrm{Ki}_n(t)\,
  \mathrm{d}t`, so a sign error or off-by-one index in the recursion
  would affect every :math:`\mathrm{Ki}_n` built from
  :math:`\mathrm{Ki}_1`. The CP test suite verifies the CP solver
  against analytically-constructed CP matrices (same
  :math:`\mathrm{Ki}_3`, same flat-source discretisation), so a
  systematic :math:`\mathrm{Ki}` bug could cancel between test and
  code. The Peierls reference uses :math:`\mathrm{Ki}_1` and a
  polynomial source representation, breaking the common-mode path.
- Because the Peierls Nyström operator resolves the flux as a
  piecewise polynomial of degree :math:`p-1` on each radial panel,
  it is sensitive to flux-shape errors that are invisible to a
  CP-vs-CP comparison with flat sources.

The canonical :math:`\mathrm{Ki}_n` recurrence evaluator lives in
:func:`~orpheus.derivations.common.kernels.ki_n_mp` (arbitrary precision
via :func:`mpmath.quad` on the A&S integral form). The
double-precision fast path for :math:`\mathrm{Ki}_3` goes through
:func:`~orpheus.derivations.continuous.flat_source_cp.geometry._ki3_mp` — a Chebyshev
interpolant built from ``ki_n_mp`` at module load (Phase B.4,
commit ``6badbe5``, Issue #94). The legacy ``BickleyTables``
20 000-point tabulation this replaced is documented historically
in :ref:`ki-table-construction`.

Nyström assembly in polar coordinates
-------------------------------------

The discretisation has three nested quadrature layers, each chosen
to match a distinct piece of the integrand's structure.

**Radial grid (composite Gauss--Legendre on** :math:`[0, R]`\ **).**
The r-axis is partitioned into the :math:`N_{\rm reg}` annular
regions :math:`[r_{k-1}, r_k]`, :math:`r_0 = 0`, :math:`r_{N} = R`.
Each region carries :math:`n_{\rm panels}` panels, each carrying
:math:`p` Gauss--Legendre nodes. Panel breakpoints coincide with
annular radii so the emission density :math:`q(r')`, which is
piecewise-smooth but has **slope discontinuities** at material
boundaries, is represented by a piecewise polynomial of degree
:math:`p-1`. This mirrors the slab composite-panel strategy.
The total number of radial Nyström unknowns is
:math:`N = N_{\rm reg} \times n_{\rm panels} \times p`. The builder
is :func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.composite_gl_r`
(shared with the sphere; the ``composite_gl_y`` alias was retired in
favour of the single ``composite_gl_r`` name).

**Azimuthal quadrature (Gauss--Legendre on** :math:`[0, \pi]`\ **).**
With :math:`n_\beta` nodes and weights :math:`w_{\beta,k}`; the
physical interval :math:`[0, 2\pi]` is folded to :math:`[0, \pi]`
by the :math:`\beta \to -\beta` symmetry already absorbed into the
:math:`1/\pi` prefactor.

**Ray-distance quadrature (Gauss--Legendre on**
:math:`[0, \rho_{\max}(r_i, \beta_k)]`\ **).**
With :math:`n_\rho` nodes per ray. The upper limit depends on both
the observer radius :math:`r_i` and the direction :math:`\beta_k`,
so the ρ-quadrature is **remapped per (i, k)** from the reference
interval :math:`[-1, 1]`. For a homogeneous bare cylinder, a fixed
:math:`\rho`-scale would under-resolve rays near the tangent
direction and over-resolve radial rays; the per-ray remapping
gives uniform relative accuracy.

**Source interpolation by Lagrange basis.**
Because the source point :math:`r'(r_i, \rho_m, \beta_k)` is
**generally not a radial quadrature node**, the emission density
at the source is expressed via the panel-local Lagrange basis:

.. math::

   q\bigl(r'_{ikm}\bigr)
     \;=\; \sum_{j=1}^{N} L_j(r'_{ikm})\,q_j,

where :math:`L_j` is the degree-:math:`(p-1)` Lagrange polynomial
supported only on the panel containing :math:`r'_{ikm}`
(piecewise-polynomial representation matching the composite GL
radial mesh). The basis is built by
:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.lagrange_basis_on_panels`.
Two properties are
enforced by L0 foundation tests:

- **Partition of unity**: :math:`\sum_j L_j(r) = 1` for any
  :math:`r \in [0, R]`. Tested in
  ``TestLagrangeBasisOnPanels.test_partition_of_unity``.
- **Polynomial reproduction**: for any polynomial of degree
  :math:`< p`, :math:`\sum_j p(r_j)\,L_j(r) = p(r)` exactly.
  Tested in
  ``TestLagrangeBasisOnPanels.test_reproduces_polynomial``.

**Assembled matrix.**
Substituting :math:`q(r'_{ikm}) = \sum_j L_j(r'_{ikm})\,q_j`
into :eq:`peierls-cylinder-equation` gives the identity-LHS form

.. math::
   :label: peierls-cylinder-nystrom

   \Sigma_t(r_i)\,\varphi_i
     \;=\; \sum_{j=1}^{N} K_{ij}\,q_j + S_{\rm bc}(r_i),

.. vv-status: peierls-cylinder-nystrom documented

with

.. math::

   K_{ij}
     \;=\; \frac{\Sigma_t(r_i)}{\pi}
       \sum_{k=1}^{n_\beta}\!\sum_{m=1}^{n_\rho}
         w_{\beta,k}\,w_{\rho,m}(r_i,\beta_k)\,
         \mathrm{Ki}_1(\tau_{ikm})\,L_j(r'_{ikm}).

The kernel matrix is assembled by ``build_volume_kernel`` in
:mod:`orpheus.derivations.continuous.peierls_nystrom.cylinder`. The per-sample optical
depth :math:`\tau_{ikm}` is computed by ``_optical_depth_along_ray``,
which walks annular boundary crossings as described next.

Ray optical-depth walker
------------------------

The optical depth along the ray from :math:`r_i` in direction
:math:`\beta` over distance :math:`\rho`,

.. math::
   :label: peierls-cylinder-ray-optical-depth

   \tau(r_i, \rho, \beta)
     \;=\; \int_{0}^{\rho}
       \Sigma_t\!\bigl(r'(r_i, s, \beta)\bigr)\,\mathrm{d}s,

.. vv-status: peierls-cylinder-ray-optical-depth documented

is piecewise-constant in the integrand. The boundary crossings
:math:`|\mathbf{r}(s)|^{2} = r_k^{2}` give the quadratic

.. math::

   s^{2} + 2 r_i \cos\beta\,s + (r_i^{2} - r_k^{2}) \;=\; 0,

whose roots :math:`s = -r_i \cos\beta \pm
\sqrt{r_i^{2}\cos^{2}\beta + r_k^{2} - r_i^{2}}` are the entry and
exit points for the ray crossing annulus :math:`k`. The walker
sorts all such roots in :math:`(0, \rho)`, evaluates
:math:`r_{\rm mid}` on each segment, and accumulates
:math:`\Sigma_{t,k}\cdot\Delta s`. The homogeneous case
(:math:`N_{\rm reg} = 1`) short-circuits to
:math:`\tau = \Sigma_t\,\rho` for speed.

The walker is L0-verified against closed-form traversals in
``tests/derivations/test_peierls_cylinder_multi_region.py``:
``TestOpticalDepthAlongRay`` covers the homogeneous short-circuit,
a ray staying in the outer annulus, a ray crossing one inner
boundary, a ray through the axis traversing three annular
segments, and a tangent ray that grazes the inner annulus.

Relationship to the :math:`\tau^{\pm}` chord walker
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The module historically exposed a separate walker
``optical_depths_pm`` that computed the same-side (:math:`\tau^{+}`)
and through-centre (:math:`\tau^{-}`) branches of the chord form for
reference / verification purposes. It was retired in favour of the
unified :meth:`~orpheus.derivations.continuous.peierls_nystrom.geometry.CurvilinearGeometry.optical_depth_along_ray`,
which now carries the piecewise chord-integration for both branches.
The chord form is not used by ``build_volume_kernel`` (which operates
in polar coordinates) but is the primitive that would be needed for
any future Schur-complemented white-BC boundary closure (see
:ref:`peierls-cylinder-white-bc` below), where the relevant
variable is the chord impact parameter :math:`y`, not the
observer-centred :math:`\rho`. Its L0 tests live in
``tests/derivations/test_peierls_cylinder_geometry.py``.

.. _peierls-cylinder-row-sum:

Row-sum identity — homogeneous vs multi-region
----------------------------------------------

The row-sum identity is the single most diagnostic
self-consistency check for the Peierls operator. It isolates the
prefactor, the kernel normalisation, and the quadrature against
the multiplicative-factor class of bugs that otherwise only show
up as a biased eigenvalue. The cylindrical identity is subtler
than the slab's because of how :math:`\Sigma_t` interacts with
the ray integral.

**Homogeneous cylinder.** For a bare homogeneous cylinder of
radius :math:`R`, the infinite-medium identity for the
identity-LHS kernel :eq:`peierls-cylinder-nystrom` is

.. math::

   \sum_{j=1}^{N} K_{ij} \;=\; \Sigma_t(r_i) \qquad (R \to \infty).

The finite-cylinder deficit
:math:`\Sigma_t - \sum_j K_{ij}` equals the uncollided escape
probability :math:`\Sigma_t\,P_{\rm esc}(r_i)` times :math:`\Sigma_t`
(a standard result: :cite:`BellGlasstone1970` §2.7; :cite:`Hebert2020`
Eq. 3.101), and for :math:`R = 10` MFP this deficit is
:math:`< 10^{-3}` at :math:`r_i \le R/2`. Tested in
``TestRowSumIdentity.test_interior_row_sum_equals_sigma_t`` in
``tests/derivations/test_peierls_cylinder_prefactor.py``.

**Multi-region cylinder.** The naive "apply :math:`K` to
:math:`q \equiv 1`" identity **fails** when :math:`\Sigma_t` is
piecewise-constant across annuli. The reason is visible in the
change-of-variables :math:`u = \tau(\rho)`:

.. math::

   \int_{0}^{\rho_{\max}}\!\mathrm{Ki}_1\!\bigl(\tau(\rho)\bigr)\,
     \mathrm{d}\rho
     \;=\; \int_{0}^{\tau_{\max}}\!
       \frac{\mathrm{Ki}_1(u)}{\Sigma_t\!\bigl(r'(u)\bigr)}\,\mathrm{d}u.

The :math:`1/\Sigma_t` in the integrand depends on **where along
the ray** the source point sits, so the identity collapses only
if :math:`\Sigma_t` is constant.

The correct multi-region identity is obtained by applying
:math:`K` to the source :math:`q = \Sigma_t` — physically, the
pure-scatter emission density that sustains a spatially uniform
flux :math:`\varphi \equiv 1`:

.. math::
   :label: peierls-cylinder-row-sum-identity

   \sum_{j=1}^{N} K_{ij}\,\Sigma_t(r_j) \;=\; \Sigma_t(r_i)
   \qquad\text{(multi-region, } R \to \infty\text{)}.

.. vv-status: peierls-cylinder-row-sum-identity documented

The :math:`\Sigma_t(r_j)` factor absorbs the :math:`1/\Sigma_t`
left behind by the change of variables, restoring
:math:`\int \mathrm{Ki}_1(u)\,\mathrm{d}u = 1` independently of
:math:`\Sigma_t` variation along the ray. This is the identity
actually tested in
``TestMultiRegionKernel.test_K_applied_to_sig_t_gives_local_sig_t``
— it applies :math:`K` to :math:`q_j = \Sigma_t(r_j)` and
verifies recovery of :math:`\Sigma_t(r_i)` at every interior
observer, to 0.5 % accuracy on an
:math:`(r_1, R) = (3, 10)`-MFP two-annulus configuration with
:math:`\Sigma_{t,{\rm inner}} = 0.8, \Sigma_{t,{\rm outer}} = 1.4`.

.. warning::

   A test that applied :math:`K` to :math:`\mathbf{1}` and
   compared to :math:`\Sigma_t(r_i)` would silently fail for the
   multi-region case even when the implementation is correct.
   The row sum
   :math:`\sum_j K_{ij}` instead equals a ray-path-weighted
   average of :math:`1/\Sigma_t`, which is not a local quantity.
   This is the reason the multi-region test file applies
   :math:`K` to :math:`\Sigma_t`, not to :math:`\mathbf{1}`.

.. _peierls-cylinder-white-bc:

Boundary conditions
-------------------

**Vacuum.** :math:`S_{\rm bc} \equiv 0`. The kernel :math:`K` is
the full operator; the eigenvalue problem is

.. math::

   \bigl[\mathrm{diag}(\Sigma_t) - K\,\mathrm{diag}(\Sigma_s)\bigr]
     \varphi \;=\; \frac{1}{k}\,K\,\mathrm{diag}(\nu\Sigma_f)\,\varphi,

solved by fission-source power iteration in
:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.solve_peierls_1g`
(with ``geometry=_pg.CYLINDER_1D`` and ``boundary="vacuum"``); the
sphere-/cylinder-specific façades in
:mod:`orpheus.derivations.continuous.peierls_nystrom.cylinder` /
:mod:`orpheus.derivations.continuous.peierls_nystrom.sphere` are registry-only
shims after the Issue #138 collapse. This is the closure that is
currently implemented; it is used for the Sanchez tie-point
verification below.

**White / reflective.** The cylinder's lateral surface is a
continuous 1-parameter family of boundary points; every
(impact-parameter :math:`y`, travel-direction) pair produces a
distinct re-entering ray. The slab's rank-2 :math:`E_2` outer
product — which comes from exactly two boundary faces — does **not**
generalise. Instead, the white-BC closure is a rank-:math:`N_y`
**dense** Schur-complement block, where :math:`N_y` is the
number of :math:`y`-quadrature nodes used to represent the
outgoing/incoming currents on the cylinder surface.

The two implementable options are:

(a) **Reduced form.** Eliminate the boundary-current unknowns
    by Schur complement; the result is an effective boundary
    source :math:`S_{\rm bc}(r_i)` that is a smooth integral
    operator of the volume unknowns :math:`\varphi_j`. :cite:`Sanchez1982`
    §IV.B.3 uses this form.
(b) **Full form.** Keep the boundary currents :math:`J^\pm(y_\ell)`
    as explicit unknowns and solve the coupled
    :math:`(\varphi, J^{+}, J^{-})` block system. :cite:`Hebert2020`
    uses this form because the coupling block is trivially
    populated from the :math:`\tau^{\pm}` walker.

For :math:`N_y = O(50\text{--}100)` the dense block is a
:math:`50 \times 50` to :math:`100 \times 100` matrix, negligible
relative to the :math:`O(N^{3})` radial LU factorisation.

.. note::

   The cylinder white-BC closure is **not yet implemented**;
   :func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.solve_peierls_1g`
   with ``geometry=_pg.CYLINDER_1D`` handles vacuum BC only. The
   :math:`\tau^{\pm}` walker (``optical_depths_pm``) that lives
   alongside ``build_volume_kernel`` is the primitive needed for
   either option (a) or option (b); it is retained for the planned
   extension. Tracked under `Issue #103
   <https://github.com/deOliveira-R/ORPHEUS/issues/103>`_
   (higher-rank N1 white-BC closure for cylinder + sphere).

Verification evidence
---------------------

Three independent checks gate the Peierls cylinder implementation.

**Sanchez--McCormick 1982 tie-point.** For a bare 1-group
homogeneous cylinder with :math:`\Sigma_t = 1` cm⁻¹,
:math:`(\Sigma_s, \nu\Sigma_f) = (0.5, 0.75)` — giving
:math:`k_\infty = \nu\Sigma_f/\Sigma_a = 1.5` — :cite:`Sanchez1982`
Table IV reports a critical radius :math:`R = 1.9798` cm.
At that :math:`R` the present solver gives

.. math::

   k_{\rm eff}(R = 1.9798) \;=\; 1.00421 \pm 10^{-5}

under polar-quadrature refinement
:math:`(n_\beta, n_\rho) = (20, 20) \to (32, 32)`. The 0.42 %
offset from unity reflects the ambiguous scatter / fission split
in the Sanchez ``c = 1.5`` problem definition: the 1-group
:math:`k_{\rm eff}` is **not invariant** under the split at fixed
:math:`k_\infty`, because :math:`\Sigma_s` enters the resolvent
:math:`(\Sigma_t\mathbf{I} - K\Sigma_s)^{-1}` separately from the
fission source. The reference-split lit cross-check is tracked
under `Issue #144
<https://github.com/deOliveira-R/ORPHEUS/issues/144>`_; until
then the test gate is set to a 1 % tolerance in
``TestSanchezTiePoint.test_k_eff_at_R_equals_1_dot_9798``, which
is robust to this ambiguity but tight enough to catch any
multiplicative-factor regression.

**Vacuum-BC thick-cylinder limit.** As :math:`R \to \infty`,
leakage vanishes and :math:`k_{\rm eff} \to k_\infty`. At
:math:`R = 30` MFP the solver reaches
:math:`k_{\rm eff} \approx 1.49` against :math:`k_\infty = 1.5`
(0.5 % gap), and :math:`k_{\rm eff}(R)` is **monotone increasing**
for :math:`R \in \{1.5, 3, 6, 12, 24\}` MFP. Tested in
``TestVacuumBCThickLimit``.

**Row-sum identity.** For the homogeneous :math:`R = 10` MFP
configuration, :math:`\max_i |\Sigma_t - \sum_j K_{ij}| < 10^{-3}`
at :math:`r_i \le R/2`, and the deficit decays monotonically
toward :math:`1` as :math:`r_i \to R` (escape probability rises
at the surface). The multi-region identity
:math:`\sum_j K_{ij}\,\Sigma_t(r_j) = \Sigma_t(r_i)` holds to
0.5 % on a :math:`(\Sigma_{t,{\rm inner}},
\Sigma_{t,{\rm outer}}) = (0.8, 1.4)` two-annulus problem.
Tested in
``TestRowSumIdentity`` and ``TestMultiRegionKernel`` in
``tests/derivations/test_peierls_cylinder_prefactor.py`` and
``tests/derivations/test_peierls_cylinder_multi_region.py``.

.. list-table:: Cylindrical Peierls verification summary
   :header-rows: 1
   :widths: 35 25 20 20

   * - Check
     - Tolerance
     - Status
     - Identity
   * - Sanchez tie-point :math:`R = 1.9798`
     - :math:`|k_{\rm eff} - 1| < 2\times10^{-2}`
     - passes at :math:`10^{-2}`
     - :eq:`peierls-cylinder-equation`
   * - Thick limit :math:`R = 30` MFP
     - :math:`|k_{\rm eff} - k_\infty|/k_\infty < 10^{-2}`
     - passes at :math:`5\times10^{-3}`
     - vacuum-BC fixed point
   * - Row sum (homogeneous)
     - :math:`<10^{-3}` at :math:`r_i \le R/2`
     - passes
     - :eq:`peierls-cylinder-nystrom`
   * - Row sum (multi-region)
     - :math:`<5\times10^{-3}` bulk interior
     - passes
     - :eq:`peierls-cylinder-row-sum-identity`

Relationship to the CP flat-source cylinder solver
--------------------------------------------------

The CP flat-source method for the cylinder
(:func:`~orpheus.cp.solver.solve_cp` on ``cyl1D`` meshes;
:mod:`orpheus.derivations.continuous.flat_source_cp.cylinder`) integrates the
:math:`\mathrm{Ki}_1` kernel analytically over each annulus to
produce the :math:`\mathrm{Ki}_3` second-difference formula
quoted above. The Peierls reference **bypasses that integration**
entirely:

- the kernel is :math:`\mathrm{Ki}_1`, not :math:`\mathrm{Ki}_3`,
- the spatial representation is a piecewise polynomial of degree
  :math:`p - 1` per panel, not a piecewise constant,
- the ray integration is performed numerically in polar
  coordinates, not analytically over rectangular annular regions.

So the two methods share almost nothing except the underlying
integral equation. A sign error, off-by-one index, or factor-of-2
in the :math:`\mathrm{Ki}_3` second-difference formula — which
would cancel between the CP solver and a CP-self-verification
test — would be caught by the Peierls reference. Conversely, a
systematic error in :math:`\mathrm{Ki}_1` evaluation would be
caught by the CP eigenvalue tests (which use pre-tabulated
:math:`\mathrm{Ki}_3` values). The two together triangulate the
cylindrical integral-transport stack.

Numerical cost
--------------

The :math:`(\beta, \rho)` tensor-product quadrature is the dominant
cost. For each observer :math:`r_i` and each :math:`\beta_k`, the
kernel assembly evaluates :math:`\mathrm{Ki}_1` at :math:`n_\rho`
points via :func:`~orpheus.derivations.common.kernels.ki_n_mp` (mpmath
at ``dps`` precision), which is :math:`O(N \cdot n_\beta \cdot
n_\rho)` kernel evaluations. For :math:`N = 10` radial nodes,
:math:`(n_\beta, n_\rho) = (24, 24)`, ``dps = 20``, kernel
assembly takes :math:`\approx 3` s on current hardware; eigenvalue
power iteration is a further :math:`O(N^{3})` LU per iteration,
typically converging in 20--30 iterations to
:math:`10^{-10}` eigenvalue tolerance.

Short-circuit: the homogeneous single-region branch of
``_optical_depth_along_ray`` returns :math:`\Sigma_t \rho` without
sorting crossings, making the bare-cylinder case
:math:`\sim\!2\times` faster than the multi-region path.

.. seealso::

   :mod:`orpheus.derivations.continuous.peierls_nystrom.cylinder` — registry-only
   module; binds the cylinder ``GEOMETRY`` singleton and ships the
   ``_build_peierls_cylinder_*_case`` continuous-reference
   constructors. The Nyström solver implementation lives in
   :mod:`~orpheus.derivations.continuous.peierls_nystrom.geometry`.

   :class:`~orpheus.derivations.continuous.peierls_nystrom.geometry.PeierlsSolution`
   — canonical result container with radial node positions, flux
   values, :math:`k_{\rm eff}`, and ``geometry_kind`` discriminator.
   Same dataclass for slab / cylinder / sphere.

   :func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.solve_peierls_1g`
   — 1-group eigenvalue driver (call with
   ``geometry=_pg.CYLINDER_1D`` and ``boundary="vacuum"`` for the
   vacuum-BC closure described above).

   ``tests/derivations/test_peierls_cylinder_geometry.py`` — L0
   tests for ``composite_gl_y`` and ``optical_depths_pm``.

   ``tests/derivations/test_peierls_cylinder_prefactor.py`` — L0
   row-sum-identity tests (homogeneous).

   ``tests/derivations/test_peierls_cylinder_multi_region.py`` —
   L0 multi-region optical-depth walker, Lagrange-basis
   foundation tests, and the multi-region
   :math:`\sum_j K_{ij}\,\Sigma_t(r_j) = \Sigma_t(r_i)` identity.

   ``tests/derivations/test_peierls_cylinder_eigenvalue.py`` — L1
   Sanchez tie-point and thick-cylinder limit eigenvalue tests.


Peierls integral equation reference — sphere
=============================================

The cylindrical Peierls reference above verifies the CP flat-source
discretisation for 1-D radial geometries with 2-D translational
symmetry. The **spherical** Peierls reference closes the trio for
``sph1D`` meshes: it solves the 3-D integral transport equation on a
bare or concentric-shell sphere at arbitrary quadrature order,
providing an independent numerical reference against
:func:`~orpheus.cp.solver.solve_cp` on spherical meshes and the
flat-source spherical CP of :mod:`orpheus.derivations.continuous.flat_source_cp.sphere`.

Unlike the slab (:math:`E_1` from :math:`y, z`-integration) and the
cylinder (:math:`\mathrm{Ki}_1` from :math:`z`-integration), the
sphere **does not reduce dimensions** — the native 3-D point kernel
:math:`e^{-\tau}/(4\pi\,d^{2})` is already 3-D. Rotational symmetry
about the centre is not a translation and cannot be used to collapse
the kernel; it only constrains the *source field* to depend on
:math:`|\mathbf r'| = r'` rather than on all three coordinates. The
polar-form pivot still buys the Jacobian cancellation (the
:math:`\rho^{2}` polar volume element cancels the :math:`1/d^{2}` of
the Green's function exactly), but the resulting reduced kernel is
the bare exponential :math:`e^{-\tau}`, not any :math:`E_n` or
:math:`\mathrm{Ki}_n` function. This makes the sphere both the
simplest kernel to evaluate (plain ``np.exp``) and the most
informative cross-check for a common-mode :math:`\mathrm{Ki}_n`
recursion bug: any factor-of-two, off-by-one, or sign error in the
Bickley recurrence that would affect the cylinder cancels out of the
sphere, and vice-versa.

The implementation lives in :mod:`orpheus.derivations.continuous.peierls_nystrom.sphere`
— a thin facade over the unified
:class:`~orpheus.derivations.continuous.peierls_nystrom.geometry.CurvilinearGeometry`
(``kind = "sphere-1d"``). The sphere shares the
Lagrange-basis, composite-GL radial quadrature, optical-depth walker,
and power-iteration primitives with the cylinder verbatim; the only
geometry-specific ingredients are the kernel :math:`\kappa_d(\tau) =
e^{-\tau}`, the angular weight :math:`W_\Omega(\theta) =
\sin\theta`, the prefactor :math:`C_d = 1/2`, and the
surface-divisor :math:`A_d = R^{2}` in the rank-1 white-BC
closure. The architectural separation is documented at full length
in :doc:`/theory/references/peierls_nystrom`.

Derivation from the 3-D point kernel
-------------------------------------

Starting from the 3-D Green's function for the one-speed, isotropic
integral transport equation,

.. math::
   :label: peierls-sphere-green-3d

   G_{\rm 3D}(\mathbf r, \mathbf r')
     \;=\; \frac{e^{-\tau(\mathbf r, \mathbf r')}}
                {4\pi\,|\mathbf r - \mathbf r'|^{2}},
   \qquad
   \tau \;=\; \int_{\mathbf r'}^{\mathbf r}\!\Sigma_t(\mathbf s)\,
     \mathrm d\ell,

.. vv-status: peierls-sphere-green-3d documented

the **pointwise** scalar-flux form of the Peierls integral equation
on a bare sphere of radius :math:`R` is

.. math::

   \varphi(\mathbf r)
     \;=\; \iiint_{\rm ball}\!
       \frac{e^{-\tau(\mathbf r, \mathbf r')}}
            {4\pi\,|\mathbf r - \mathbf r'|^{2}}\,
         q(\mathbf r')\,\mathrm d^{3}r'
     \;+\; \varphi_{\rm bc}(\mathbf r).

This is the **native** 3-D kernel. Compare with the slab
:eq:`peierls-equation` (whose kernel is :math:`E_1`, obtained by
integrating the point kernel over two transverse dimensions) and
the cylinder :eq:`peierls-cylinder-green-2d` (whose kernel is
:math:`\mathrm{Ki}_1`, obtained by integrating over one axial
direction). The sphere's
point kernel **cannot be pre-integrated** because there is no
translational symmetry to exploit — a radial 1-D problem inherits
only rotational symmetry from the embedding space, and rotations
move every point on a shell to every other point on the same shell
without ever crossing the shell boundary.

.. note::

   The monotone progression of dimensional reductions — two for the
   slab, one for the cylinder, zero for the sphere — is the defining
   feature of the trio. :doc:`/theory/references/peierls_nystrom` §2 tabulates the
   reduced kernels side-by-side.

Observer-centred polar form and Jacobian cancellation
------------------------------------------------------

The native 3-D integral above is not directly tractable: the
:math:`1/|\mathbf r - \mathbf r'|^{2}` singularity at
:math:`\mathbf r' = \mathbf r` is a **volume singularity**
(non-integrable without the Jacobian of an appropriate coordinate
change). Rather than attempt a chord/impact-parameter formulation
analogous to the cylinder's (which, on the sphere, would introduce
*three* coincident endpoint singularities at a point), the
implementation uses the **equivalent polar form** centred on the
observer.

Let :math:`\theta \in [0, \pi]` be the polar angle from the outward
radial direction at :math:`\mathbf r`, :math:`\phi \in [0, 2\pi]`
the azimuth, and :math:`\rho \ge 0` the distance along the ray from
:math:`\mathbf r` in direction :math:`(\theta, \phi)`. The source
position is

.. math::
   :label: peierls-sphere-r-prime

   r'(r, \rho, \theta) \;=\;
     \sqrt{r^{2} + 2 r \rho \cos\theta + \rho^{2}},

.. vv-status: peierls-sphere-r-prime documented

**identical** to the cylinder case
:eq:`peierls-cylinder-r-prime` — the 1-D radial chord algebra does
not care whether the surrounding source field is
:math:`2`-D-symmetric (cylinder) or :math:`3`-D-symmetric
(sphere). That is the architectural insight that let
:class:`~orpheus.derivations.continuous.peierls_nystrom.geometry.CurvilinearGeometry`
share ray primitives between the two geometries.

The 3-D volume element in observer-centred coordinates is
:math:`\rho^{2}\,\sin\theta\,\mathrm d\rho\,\mathrm d\theta\,
\mathrm d\phi`. Combined with the 3-D Green's function the
integrand becomes

.. math::

   \frac{e^{-\tau}}{4\pi\,\rho^{2}}\,\cdot\,
   \rho^{2}\,\sin\theta\,\mathrm d\rho\,\mathrm d\theta\,
   \mathrm d\phi
   \;=\;
   \frac{\sin\theta}{4\pi}\,e^{-\tau}\,\mathrm d\rho\,
   \mathrm d\theta\,\mathrm d\phi.

The :math:`\rho^{2}` polar volume factor **cancels the**
:math:`1/\rho^{2}` **of the Green's function exactly** — the
Jacobian-cancellation trick derived in full generality in
:doc:`/theory/references/peierls_nystrom` §3. The integrand is now a bounded,
polynomial-smooth function on the whole integration domain (modulo
:math:`\Sigma_t` jumps, which are handled by the composite radial
grid described below). Because the source field is radially
symmetric, nothing in the integrand depends on :math:`\phi`, so the
azimuthal integral collapses to a factor of :math:`2\pi`:

.. math::
   :label: peierls-sphere-polar

   \varphi(r)
     \;=\; \frac{1}{2}\!
       \int_{0}^{\pi}\!\sin\theta\,\mathrm d\theta\!
       \int_{0}^{\rho_{\max}(r, \theta)}\!\!
         e^{-\tau(r, \rho, \theta)}\,
         q\!\bigl(r'(r, \rho, \theta)\bigr)\,\mathrm d\rho
     \;+\; \varphi_{\rm bc}(r).

.. vv-status: peierls-sphere-polar documented

The prefactor :math:`1/2` is :math:`2\pi / (4\pi) = 1/2`, and the
:math:`\sin\theta` weight remains from the solid-angle element
:math:`\mathrm d\Omega = \sin\theta\,\mathrm d\theta\,\mathrm d\phi`.

.. note::

   **No** :math:`\beta \to -\beta` **folding is needed** on the
   sphere. The cylindrical polar form
   :eq:`peierls-cylinder-polar` folds :math:`\beta \in [0, 2\pi]`
   to :math:`[0, \pi]` (buying a factor of 2 into the
   :math:`1/\pi` prefactor) because the integrand there is
   :math:`\beta`-symmetric. On the sphere the polar angle
   :math:`\theta \in [0, \pi]` already covers the full hemisphere
   of ray directions at the observer, and :math:`\sin\theta \ge 0`
   on that interval — no folding is available or needed.

The upper radial limit :math:`\rho_{\max}` is the intersection of
the ray with the outer sphere,

.. math::
   :label: peierls-sphere-rho-max

   \rho_{\max}(r, \theta)
     \;=\; -r\cos\theta
         + \sqrt{r^{2}\cos^{2}\theta + R^{2} - r^{2}},

.. vv-status: peierls-sphere-rho-max documented

**identical** to the cylinder :eq:`peierls-cylinder-rho-max`.
Verified by ``TestSphereRhoMax`` in
``tests/derivations/test_peierls_sphere_geometry.py``, which covers
the radial-outward ray (:math:`\rho_{\max} = R - r`), the
radial-inward through-diameter ray (:math:`\rho_{\max} = R + r`),
the tangential ray from the centre (:math:`\rho_{\max} = R`), and
the observer-on-surface outward ray (:math:`\rho_{\max} = 0`).

Writing the identity-LHS form used by the eigenvalue driver, the
canonical spherical Peierls equation solved by this module is

.. math::
   :label: peierls-sphere-equation

   \Sigma_t(r_i)\,\varphi(r_i)
     \;=\; \frac{\Sigma_t(r_i)}{2}\!
       \int_{0}^{\pi}\!\sin\theta\,\mathrm d\theta\!
       \int_{0}^{\rho_{\max}(r_i, \theta)}\!\!
         e^{-\tau(r_i, \rho, \theta)}\,
         q\!\bigl(r'(r_i, \rho, \theta)\bigr)\,\mathrm d\rho
     \;+\; S_{\rm bc}(r_i).

.. (vv-status rationale) derivation: Intermediate step in a derivation chain whose capstone identity is tested elsewhere. Observer-centred polar form of the Peierls sphere equation — intermediate derivation step before the Nyström assembly. The assembled system is verified by peierls-unified row-sum and elementwise tests in test_peierls_reference.py.
.. vv-status: peierls-sphere-equation documented


.. vv-status: peierls-sphere-equation tested

The sphere test files carry
``@pytest.mark.verifies("peierls-unified")``, which is the coarse
label shared across the whole unified polar-form implementation
while finer-grained per-equation labels are retrofitted in a
follow-up V&V harness pass. The ``vv-status: ... tested`` annotation
above reflects that coverage.

Why :math:`e^{-\tau}` and not :math:`E_3` / :math:`\mathrm{Ki}_3`
------------------------------------------------------------------

The flat-source CP method for the sphere
(:mod:`orpheus.derivations.continuous.flat_source_cp.sphere`) uses a second-difference
formula in the :math:`E_3` function (see :eq:`second-diff-sph` above)
because it averages the pointwise :math:`e^{-\tau}` kernel twice —
once over the source region :math:`j` and once over the target
region :math:`i` — and the **double antiderivative** of
:math:`e^{-\tau}/d^{2}` along a chord produces exactly
:math:`E_3(\tau)`. This is the sphere-specific instance of the
general second-difference identity:

.. list-table:: Pointwise kernels vs flat-source second-differences
   :header-rows: 1
   :widths: 18 30 30

   * - Geometry
     - Pointwise kernel
     - Flat-source second-difference
   * - Slab
     - :math:`E_1`
     - :math:`E_3`
   * - Cylinder
     - :math:`\mathrm{Ki}_1`
     - :math:`\mathrm{Ki}_3`
   * - Sphere
     - :math:`e^{-\tau}`
     - :math:`E_3`

The sphere is alone in the right column: the second antiderivative
of :math:`e^{-\tau}` involves :math:`E_1` and :math:`E_2` depending
on which combination of chord endpoints is averaged, and the
particular combination that appears for a concentric-shell spherical
geometry happens to collapse to :math:`E_3`. Full derivation in
``orpheus/derivations/continuous/flat_source_cp/sphere.py``; see :eq:`second-diff-sph` and
:eq:`self-sph` for the resulting CP matrix elements, and
:eq:`rcp-from-double-antideriv` for the general second-difference
identity that specialises to the sphere via the same
double-antiderivative argument used for the slab and cylinder.

The Peierls reference solves for the **pointwise** flux
:math:`\varphi(r_i)` — not a region-averaged collision rate — so it
uses the raw :math:`e^{-\tau}` kernel directly and **bypasses** the
analytic chord averaging that produces :math:`E_3`. This makes the
sphere Peierls a clean cross-check on the CP flat-source
construction:

- A sign error, off-by-one, or factor-of-two in the
  :math:`E_n(\tau)` recurrence for :math:`n \ge 2` would affect the
  CP :math:`E_3` second-difference formula *but not* the Peierls
  :math:`e^{-\tau}` evaluation (which is just ``np.exp``). A
  systematic :math:`E_n` bug would therefore cancel between CP
  solver and CP-self-verification test but be caught cleanly by
  the Peierls reference.
- Because the Peierls Nyström operator resolves the flux as a
  piecewise polynomial of degree :math:`p - 1` on each radial
  panel, it is sensitive to flux-shape errors that are invisible to
  a flat-source CP-vs-CP comparison. The Phase-A Peierls-vs-CP
  flux-shape test
  (``TestCPvsPeierlsSphereAtThickR.test_flux_shape_agrees_at_thick_R``
  in ``tests/cp/test_peierls_sphere_flux.py``) is the first-order
  check that the CP flat-source approximation recovers the correct
  pointwise flux in the thick-sphere limit where the approximation
  is asymptotically exact.

The :math:`e^{-\tau}` kernel also avoids the common-mode
:math:`\mathrm{Ki}_n` recursion path: the cylinder Peierls depends
on :func:`~orpheus.derivations.common.kernels.ki_n_mp` (via the
mpmath-backed :math:`\mathrm{Ki}_1` evaluator, as the Phase B.4
retirement of ``BickleyTables`` routes all cylindrical kernel
evaluations through a single canonical primitive), but a sphere
Peierls run makes **zero** calls into any :math:`\mathrm{Ki}_n`
code path — just :func:`numpy.exp`. The two references triangulate
the integral-transport stack from orthogonal angles.

Nyström assembly in polar coordinates
-------------------------------------

The sphere discretisation mirrors the cylinder's three-layer polar
quadrature; each layer is dispatched through the unified
:class:`~orpheus.derivations.continuous.peierls_nystrom.geometry.CurvilinearGeometry`
with ``kind = "sphere-1d"``.

**Radial grid (composite Gauss–Legendre on** :math:`[0, R]`\ **).**
The r-axis is partitioned into :math:`N_{\rm reg}` concentric
shells :math:`[r_{k-1}, r_k]`, :math:`r_0 = 0`, :math:`r_N = R`.
Each shell carries :math:`n_{\rm panels}` panels, each carrying
:math:`p` Gauss–Legendre nodes. Panel breakpoints coincide with
shell radii so the emission density :math:`q(r')` — which is
piecewise-smooth but has slope discontinuities at material
boundaries — is represented by a piecewise polynomial of degree
:math:`p - 1`. This is the same strategy as the cylinder. The
total number of radial Nyström unknowns is
:math:`N = N_{\rm reg} \cdot n_{\rm panels} \cdot p`. The builder is
:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.composite_gl_r`
(shared with the cylinder; the sphere module in
:mod:`orpheus.derivations.continuous.peierls_nystrom.sphere` calls it
directly).

Verified by ``TestSphereCompositeRadialGL`` in
``tests/derivations/test_peierls_sphere_geometry.py``: the weighted
integrals :math:`\int_0^R 1\,\mathrm dr = R` and
:math:`\int_0^R 4\pi r^{2}\,\mathrm dr = \tfrac{4}{3}\pi R^{3}`
recover the analytic values to machine precision under the
composite grid.

**Polar quadrature (Gauss–Legendre on** :math:`[0, \pi]`\ **).**
With :math:`n_\theta` nodes and weights :math:`w_{\theta, k}`;
the physical interval :math:`[0, \pi]` is **not folded** (see the
note under :eq:`peierls-sphere-polar`). The :math:`\sin\theta`
factor from :eq:`peierls-sphere-equation` is applied inside the
assembly by
:meth:`CurvilinearGeometry.angular_weight`; for
``kind = "sphere-1d"`` this returns :math:`\sin\theta` evaluated at
each node, distinguishing the sphere from the cylinder's constant
:math:`W_\Omega = 1`.

**Ray-distance quadrature (Gauss–Legendre on**
:math:`[0, \rho_{\max}(r_i, \theta_k)]`\ **).**
With :math:`n_\rho` nodes per ray. The upper limit depends on both
the observer radius and the polar angle, so the
:math:`\rho`-quadrature is **remapped per** :math:`(i, k)` from the
reference interval :math:`[-1, 1]` to :math:`[0, \rho_{\max}]`.
For a homogeneous bare sphere, a fixed :math:`\rho`-scale would
under-resolve the long rays that pass near the diameter and
over-resolve short outward rays near the surface; the per-ray
remapping gives uniform relative accuracy.

**Source interpolation by Lagrange basis.**
Because :math:`r'(r_i, \rho_m, \theta_k)` is generally not a radial
quadrature node, the emission density at the source point is
expressed via the panel-local Lagrange basis:

.. math::

   q\!\bigl(r'_{ikm}\bigr)
     \;=\; \sum_{j=1}^{N} L_j\!\bigl(r'_{ikm}\bigr)\,q_j,

where :math:`L_j` is the degree-:math:`(p-1)` Lagrange polynomial
supported only on the panel containing :math:`r'_{ikm}`. The basis
is shared with the cylinder
(:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.lagrange_basis_on_panels`);
partition of unity
and polynomial reproduction are L0-verified in the cylinder's
``TestLagrangeBasisOnPanels`` and carry over to the sphere case
without modification.

**Assembled matrix.**
Substituting into :eq:`peierls-sphere-equation` gives the
identity-LHS form

.. math::
   :label: peierls-sphere-nystrom

   \Sigma_t(r_i)\,\varphi_i
     \;=\; \sum_{j=1}^{N} K_{ij}\,q_j + S_{\rm bc}(r_i),

.. (vv-status rationale) derivation: Intermediate step in a derivation chain whose capstone identity is tested elsewhere. Discretised Nyström assembly of the sphere equation — covered by peierls-unified L1 elementwise and row-sum tests against the slab/sphere uniform-source closed forms.
.. vv-status: peierls-sphere-nystrom documented


.. vv-status: peierls-sphere-nystrom tested

with

.. math::

   K_{ij}
     \;=\; \frac{\Sigma_t(r_i)}{2}\!
       \sum_{k=1}^{n_\theta}\!\sum_{m=1}^{n_\rho}\!
         w_{\theta, k}\,\sin\theta_k\,w_{\rho, m}(r_i, \theta_k)\,
         e^{-\tau_{ikm}}\,L_j\!\bigl(r'_{ikm}\bigr).

The kernel matrix is assembled by ``build_volume_kernel`` in
:mod:`orpheus.derivations.continuous.peierls_nystrom.sphere`, which dispatches to
:func:`peierls_geometry.build_volume_kernel` with the sphere
geometry singleton pre-bound. The per-sample optical depth
:math:`\tau_{ikm}` is computed by the shared multi-annulus walker,
identical to the cylinder case — the 1-D radial annulus crossings
are geometry-agnostic.

Ray optical-depth walker
------------------------

The optical depth along the ray from :math:`r_i` in direction
:math:`\theta` over distance :math:`\rho`,

.. math::
   :label: peierls-sphere-ray-optical-depth

   \tau(r_i, \rho, \theta)
     \;=\; \int_{0}^{\rho}\!
       \Sigma_t\!\bigl(r'(r_i, s, \theta)\bigr)\,\mathrm ds,

.. (vv-status rationale) definition: Definitional / notation introduction. Definition of τ(r,ρ,θ) along a ray walker; used as a primitive by the assembled solver (no isolated identity to verify).
.. vv-status: peierls-sphere-ray-optical-depth documented


.. vv-status: peierls-sphere-ray-optical-depth tested

shares **the same walker** as the cylinder's
:eq:`peierls-cylinder-ray-optical-depth`: the boundary crossings
:math:`|\mathbf{r}(s)|^{2} = r_k^{2}` give the same quadratic in
:math:`s`, whose roots partition the ray into annular segments of
constant :math:`\Sigma_t`. Because the embedding ambient space
(2-D disc vs 3-D ball) enters only through the solid-angle weight
of :eq:`peierls-sphere-equation` and the kernel :math:`\kappa_d`,
the walker itself — which only sees the 1-D radial :math:`\Sigma_t`
profile and the 1-D chord algebra — is reusable verbatim.

L0-verified against closed-form traversals in
``tests/derivations/test_peierls_sphere_geometry.py``:

- ``TestSphereOpticalDepthAlongRay.test_homogeneous_1region_linear_in_rho``
  — short-circuit :math:`\tau = \Sigma_t\,\rho` for a bare sphere.
- ``test_scales_linearly_with_sig_t`` — doubling :math:`\Sigma_t`
  doubles :math:`\tau` at every :math:`(r_i, \theta, \rho)`.
- ``test_two_annulus_radial_transit`` — radial outward ray that
  crosses from the inner shell into the outer shell, pinning the
  per-segment accumulation against the analytic partitioning.
- ``test_two_annulus_through_centre_diameter`` — through-centre
  ray that traverses the inner shell twice (once on each side of
  the centre) and the outer shell twice, covering the four
  boundary-crossing algebra.

Row-sum identity — homogeneous and multi-region
-----------------------------------------------

The row-sum identity is the single most diagnostic consistency
check for the Peierls operator on any geometry. For the sphere the
structure is **identical** to the cylinder's
:ref:`peierls-cylinder-row-sum`, and the same
:math:`u = \tau(\rho)` change-of-variables argument carries
through: applying :math:`K` to a pure-scatter emission density
:math:`q \equiv \Sigma_t` reproduces the spatially uniform
:math:`\varphi \equiv 1` because the :math:`1/\Sigma_t` factor left
by the Jacobian is absorbed.

**Homogeneous sphere.** For a bare homogeneous sphere of radius
:math:`R`, the infinite-medium identity for the identity-LHS
kernel :eq:`peierls-sphere-nystrom` is

.. math::

   \sum_{j=1}^{N} K_{ij} \;=\; \Sigma_t(r_i) \qquad (R \to \infty).

The finite-sphere deficit :math:`\Sigma_t - \sum_j K_{ij}` equals
:math:`\Sigma_t \cdot P_{\rm esc}(r_i)` (the uncollided escape
probability weighted by :math:`\Sigma_t`). For :math:`R = 10` MFP,
:math:`\max_i |\Sigma_t - \sum_j K_{ij}| < 10^{-3}` at
:math:`r_i \le R/2`. Tested in
``TestSphereRowSumIdentity.test_interior_row_sum_equals_sigma_t``
in ``tests/derivations/test_peierls_sphere_prefactor.py``.

The deficit is **monotone increasing** from centre to surface
(``test_deficit_grows_toward_boundary``), and shrinks under
quadrature refinement at every interior observer
(``test_convergence_under_quadrature_refinement``).

**Multi-region sphere.** The naive "apply :math:`K` to
:math:`\mathbf 1`" identity fails when :math:`\Sigma_t` is
piecewise-constant across shells, for the same reason as the
cylinder: the change of variables :math:`u = \tau(\rho)` gives

.. math::

   \int_{0}^{\rho_{\max}}\!e^{-\tau(\rho)}\,\mathrm d\rho
     \;=\; \int_{0}^{\tau_{\max}}\!
       \frac{e^{-u}}{\Sigma_t\!\bigl(r'(u)\bigr)}\,\mathrm du,

and the :math:`1/\Sigma_t` factor depends on where along the ray
the source point sits. The correct multi-region identity is
obtained by applying :math:`K` to :math:`q = \Sigma_t`:

.. math::
   :label: peierls-sphere-row-sum-identity

   \sum_{j=1}^{N} K_{ij}\,\Sigma_t(r_j) \;=\; \Sigma_t(r_i)
   \qquad\text{(multi-region, } R \to \infty\text{)}.

.. vv-status: peierls-sphere-row-sum-identity documented

The :math:`\Sigma_t(r_j)` factor absorbs the :math:`1/\Sigma_t`
left behind by the change of variables, restoring
:math:`\int_0^\infty e^{-u}\,\mathrm du = 1` independently of
:math:`\Sigma_t` variation along the ray.

.. warning::

   As on the cylinder, a test that applied :math:`K` to
   :math:`\mathbf 1` and compared to :math:`\Sigma_t(r_i)` would
   silently fail for the multi-region case even when the
   implementation is correct. The row sum :math:`\sum_j K_{ij}`
   instead equals a ray-path-weighted average of
   :math:`1/\Sigma_t`, which is not a local quantity. The sphere
   test suite mirrors the cylinder pattern: the **homogeneous**
   identity uses :math:`\mathbf 1` because :math:`\Sigma_t` is
   constant; any future multi-region row-sum test must apply
   :math:`K` to :math:`\Sigma_t`, not to :math:`\mathbf 1`.

Surface-to-volume Green's function :math:`G_{\rm bc}`
------------------------------------------------------

White-BC closure needs the sphere's surface-to-volume Green's
function :math:`G_{\rm bc}(r_i)` — the scalar flux at interior
observer :math:`r_i` induced by a **unit uniform isotropic inward
partial current** :math:`J^{-}` on the spherical surface. This
section derives the compact observer-centred form used by the
implementation and documents the design choice to parametrise by
directions at the observer rather than area elements on the
surface.

For a uniform isotropic inward partial current :math:`J^{-}` on
the sphere, the inward angular flux on the surface is
:math:`\psi_{\rm in} = J^{-}/\pi`, since the partial current and
the isotropic angular flux are related by :math:`J^{-} =
\int_{\Omega \cdot \hat n < 0}|\Omega \cdot \hat n|\,\psi\,
\mathrm d\Omega = \pi\,\psi_{\rm in}` for an isotropic inward
hemisphere. The scalar flux at interior observer :math:`r_i` is
obtained by integrating the attenuated emission over **directions
at the observer**:

.. math::

   \varphi(r_i)
     \;=\; \psi_{\rm in}\!\int_{4\pi}\!
       e^{-\tau_{\rm surf}(r_i, \Omega)}\,\mathrm d\Omega
     \;=\; \frac{J^{-}}{\pi}\,\cdot\,2\pi\!\int_{0}^{\pi}\!
       \sin\theta\,e^{-\tau_{\rm surf}(r_i, \theta)}\,
       \mathrm d\theta,

where
:math:`\tau_{\rm surf}(r_i, \theta) = \int_0^{\rho_{\max}(r_i,\theta)}
\Sigma_t(r'(s))\,\mathrm ds` is the optical depth along the ray
from :math:`r_i` in direction :math:`\theta` to the surface.
Dividing by :math:`J^{-}`:

.. math::
   :label: peierls-sphere-G-bc

   G_{\rm bc}^{\rm sph}(r_i)
     \;=\; 2\!\int_{0}^{\pi}\!\sin\theta\,
       e^{-\tau_{\rm surf}(r_i, \theta)}\,\mathrm d\theta.

.. (vv-status rationale) derivation: Intermediate step in a derivation chain whose capstone identity is tested elsewhere. Surface-to-volume Green's function definition for the sphere; the assembled G-bc is verified by peierls-unified at L1 in TestSphereKernelRowSum and the per-face tests.
.. vv-status: peierls-sphere-G-bc documented


.. vv-status: peierls-sphere-G-bc tested

.. note::

   **Observer parametrisation vs surface parametrisation.** The
   standard textbook derivation writes :math:`G_{\rm bc}(r_i)` as
   an integral over the boundary surface area element
   :math:`\mathrm dA_{\rm surf} = R^{2}\sin\theta'\,\mathrm d\theta'\,
   \mathrm d\phi'` with a :math:`\cos\theta'` Lambertian weight
   (the projection of the inward normal onto the ray) and a
   :math:`1/d^{2}` geometric attenuation:

   .. math::

      G_{\rm bc}^{\rm surf}(r_i)
        \;=\; \frac{1}{\pi}\!\iint_{\rm surf}\!
          \frac{\cos\theta'\,e^{-\tau_{\rm surf}}}{d(r_i, \theta')^{2}}\,
            \mathrm dA_{\rm surf}.

   The **observer form** :eq:`peierls-sphere-G-bc` and the
   **surface form** are equivalent via change of variables: for
   every inward-pointing ray at the observer there is one and only
   one entry point on the surface, and the Jacobian of the mapping
   exactly cancels the :math:`1/d^{2}` attenuation and the
   :math:`\cos\theta'` Lambertian weight. The observer
   parametrisation is **structurally simpler** — no
   :math:`\cos\theta'`, no :math:`1/d^{2}`, no extra branch
   choice — and the angular range is the natural :math:`[0, \pi]`
   of polar-angle integration.

   This is the same Jacobian-cancellation principle that eliminates
   the :math:`1/\rho^{2}` volume singularity in the polar form of
   the volume kernel. Choosing the observer parametrisation is the
   design decision that makes :math:`G_{\rm bc}` a smooth integral
   that Gauss–Legendre quadrature handles spectrally.

**Vacuum limit (sanity check).** As :math:`\Sigma_t \to 0` the
exponential collapses to unity and

.. math::

   G_{\rm bc}^{\rm sph}(r_i)\,\big|_{\Sigma_t = 0}
     \;=\; 2\!\int_{0}^{\pi}\!\sin\theta\,\mathrm d\theta
     \;=\; 2 \cdot 2 \;=\; 4.

Physically: a uniform isotropic inward partial current of strength
:math:`J^{-}` on an empty ball fills the interior with scalar flux
:math:`4 J^{-}` (:math:`4\pi` sr of angular flux, each mode
:math:`\psi_{\rm in} = J^{-}/\pi`, integrated gives :math:`4\pi \cdot
J^{-}/\pi = 4 J^{-}`). This limit is tested in
``TestSphereGBCVacuumLimit.test_vacuum_G_bc_is_four`` in
``tests/derivations/test_peierls_sphere_prefactor.py``: a
:math:`\Sigma_t R = 10^{-8}` sphere gives
:math:`G_{\rm bc} = 4` to :math:`10^{-5}` at every interior
observer.

Rank-1 white-BC closure — geometry-aware surface divisor
----------------------------------------------------------

Under the **rank-1 Mark / isotropic closure**, the white-BC
correction to the volume kernel is of outer-product form
:math:`K_{\rm bc}[i, j] = u_i\,v_j` with

.. math::
   :label: peierls-rank1-white-bc-factors

   u_i \;=\; \frac{\Sigma_t(r_i)\,G_{\rm bc}(r_i)}{A_d},
   \qquad
   v_j \;=\; r_j^{d-1}\,w_j\,P_{\rm esc}(r_j),

where :math:`A_d` is the cell-surface measure, :math:`d \in
\{2, 3\}` is the ambient dimension, and :math:`P_{\rm esc}(r_j)`
is the uncollided escape probability.

For the **sphere** (:math:`d = 3`, :math:`A_d = 4\pi R^{2}`,
volume-element area :math:`A_j = 4\pi r_j^{2} w_j`), the
:math:`4\pi` azimuthal factor cancels between :math:`A_d` and
:math:`A_j`, leaving a ratio :math:`A_j / A_d = r_j^{2} w_j /
R^{2}`. The implementation therefore uses a surface divisor of
:math:`R^{2}`:

.. math::

   u_i^{\rm sph} = \frac{\Sigma_t(r_i)\,G_{\rm bc}(r_i)}{R^{2}},
   \qquad
   v_j^{\rm sph} = r_j^{2}\,w_j\,P_{\rm esc}(r_j).

For the **cylinder** the analogous ratio is :math:`r_j w_j / R`,
so the divisor is :math:`R`. These two cases are dispatched by
:meth:`CurvilinearGeometry.rank1_surface_divisor`, which returns
:math:`R` for the cylinder and :math:`R^{2}` for the sphere.

.. warning::

   **R-vs-R² gotcha.** A sphere Peierls implementation that
   re-uses cylinder scaffolding **without** updating the
   surface-divisor from :math:`R` to :math:`R^{2}` under-counts
   :math:`u_i` by a factor of :math:`R`. The symptom is an
   enormous overestimation of :math:`k_{\rm eff}`: the white-BC
   correction, having the wrong normalisation, feeds a spuriously
   large boundary source back into the fission eigenvalue. A
   previous attempt — see the historical retraction in
   :ref:`issue-100-retraction` and the full debate in
   `Issue #100 <https://github.com/deOliveira-R/ORPHEUS/issues/100>`_
   — hit exactly this wall and reported
   :math:`k_{\rm eff} \approx 6.7` for a 1-G bare sphere where
   :math:`k_\infty = 1.5`. The
   :meth:`~CurvilinearGeometry.rank1_surface_divisor` abstraction
   exists precisely to make this mistake impossible in new code.

The implementation is thin: :func:`build_white_bc_correction` in
:mod:`orpheus.derivations.continuous.peierls_nystrom.sphere` calls
:func:`compute_G_bc`, :func:`compute_P_esc`, and assembles the
rank-1 outer product with the geometry-aware divisor — all via the
unified :class:`~orpheus.derivations.continuous.peierls_nystrom.geometry.CurvilinearGeometry`
dispatch in :func:`peierls_geometry.build_white_bc_correction`.

.. _issue-100-retraction:

Rank-1 white-BC numerical evidence (sphere)
-------------------------------------------

The :math:`R^{2}`-divisor rank-1 closure is pinned by
``TestWhiteBCRank1ErrorScan`` against the bare-sphere
:math:`k_\infty = 1.5` reference at
:math:`(\Sigma_t, \Sigma_s, \nu\Sigma_f) = (1, 0.5, 0.75)` and
:math:`n_\theta = n_\rho = 20`, :math:`n_\phi = 32`, ``dps = 25``:

.. list-table:: Rank-1 white-BC :math:`k_{\rm eff}` scan, bare sphere
   :header-rows: 1
   :widths: 14 18 18 22

   * - :math:`R` / MFP
     - :math:`k_{\rm eff}` (sphere)
     - err vs :math:`k_\infty`
     - cylinder comparator
   * - 1.0
     - 1.0963
     - 26.9 %
     - 1.19 (21 %)
   * - 2.0
     - 1.3914
     - 7.2 %
     - 1.40 (7 %)
   * - 5.0
     - 1.4897
     - 0.7 %
     - 1.48 (2 %)
   * - 10.0
     - 1.4957
     - 0.3 %
     - 1.49 (1 %)
   * - 20.0
     - 1.4945
     - 0.4 %
     - —
   * - 30.0
     - 1.4946
     - 0.4 %
     - —

The sphere's rank-1 white-BC behaviour **parallels the
cylinder's** — large error at thin :math:`R`, monotone convergence
to :math:`k_\infty` as :math:`R \to \infty`. The
beyond-:math:`R = 10` MFP plateau at
:math:`|k_{\rm eff} - k_\infty|/k_\infty \approx 0.3\,\text{--}\,0.4\%`
is the quadrature-noise floor at the cited
:math:`(n_\theta, n_\rho, n_\phi)` order, not a rank-1 closure
defect. Under refinement the floor drops further
(``TestQuadratureConvergence.test_k_eff_converges_under_refinement``
monitors this at :math:`R = 4` MFP).

Row-sum residuals under the rank-1 white-BC correction
(:math:`\max_i |K_{\rm tot} \cdot \Sigma_t - \Sigma_t|` for the
homogeneous sphere at :math:`\Sigma_t = 1`):

.. list-table:: Row-sum residuals (homogeneous sphere)
   :header-rows: 1
   :widths: 20 30 30

   * - :math:`R`
     - vacuum :math:`K`
     - white :math:`K_{\rm tot}`
   * - 2.0 MFP
     - :math:`5.4\cdot10^{-1}`
     - :math:`6.2\cdot10^{-2}`
   * - 5.0 MFP
     - :math:`4.0\cdot10^{-1}`
     - :math:`9.4\cdot10^{-3}`
   * - 10.0 MFP
     - :math:`2.9\cdot10^{-1}`
     - :math:`5.8\cdot10^{-3}`

At :math:`R = 5` MFP the rank-1 closure recovers the row-sum
identity to better than 1 % — pinned in
``TestSphereWhiteBCRowSum.test_medium_sphere_residual_below_five_percent``
and ``test_thick_sphere_residual_below_two_percent``.

**Issue #100 historical retraction.** A pre-:math:`R^{2}`-divisor
sphere prototype reported :math:`k_{\rm eff} \approx 6.7` and
attributed the blow-up to a structural rank-1 failure on the
sphere (:math:`P_{\rm esc} / G_{\rm bc}` 40 %-variation argument).
That conclusion is **retracted** — the spurious factor of
:math:`R` injected by the cylinder-port divisor
(:math:`u_i / R` instead of :math:`u_i / R^{2}`) accounts for the
inflation at :math:`R \sim 1` MFP, and the
:math:`P_{\rm esc} / G_{\rm bc}` ratio argument conflated
:math:`u_i\,v_j` with the volume-to-volume coupling. The full
debate (pre-correction "ratio varies 40 % so rank-1 fails"
argument and the divisor-bug post-mortem) lives in
`Issue #100 <https://github.com/deOliveira-R/ORPHEUS/issues/100>`_
and the sister stub in :doc:`/theory/references/peierls_nystrom` §8. The
**residual** rank-1 deficit at :math:`R \to 0` is the genuine
flat-source-accuracy limit shared by sphere and cylinder; it is
tracked under
`Issue #103 <https://github.com/deOliveira-R/ORPHEUS/issues/103>`_
(higher-rank N1 closure), not as a sphere-specific bug.

Boundary conditions — rank-1 white and vacuum
----------------------------------------------

**Vacuum.** :math:`S_{\rm bc} \equiv 0`. The kernel :math:`K` is
the full operator; the eigenvalue problem is

.. math::

   \bigl[\mathrm{diag}(\Sigma_t) - K\,\mathrm{diag}(\Sigma_s)\bigr]
     \varphi \;=\; \frac{1}{k}\,K\,\mathrm{diag}(\nu\Sigma_f)\,\varphi,

solved by fission-source power iteration in
:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.solve_peierls_1g`
(with ``geometry=_pg.SPHERE_1D``) and ``boundary="vacuum"``. This
is the clean closure: no approximation enters beyond the quadrature
orders.

**Rank-1 white.** The unified :math:`K_{\rm vol} + K_{\rm bc}`
structure with :math:`K_{\rm bc}` the rank-1 outer product derived
above, solved by the same power iteration via
:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.solve_peierls_1g`
(with ``geometry=_pg.SPHERE_1D``) and ``boundary="white"``. Accuracy is governed by the cell
:term:`optical thickness` (``test_thin_sphere_rank1_error_bounded`` /
``test_medium_sphere_rank1_error_bounded`` /
``test_thick_sphere_rank1_near_k_inf``).

.. note::

   The rank-1 closure collapses because of radial symmetry, not
   because of the specific dimensionality. For a 1-D radial sphere
   with isotropic scattering, the re-entering partial current is
   scalar (:math:`J^{-}` has a single degree of freedom by
   rotational symmetry), and the scalar balance
   :math:`J^{-} = J^{+}` is exact at every surface point. What the
   rank-1 Mark closure approximates is the **angular shape** of
   :math:`J^{-}` — treating it as isotropic in the inward
   hemisphere. That approximation is exact in the thick-cell
   limit (where the angular flux at the boundary is the integrated
   average over the slab's emission density) and degrades as
   :math:`R \to 0` (where the angular dependence on the surface is
   Fourier-rich). Issue #103 (N1) is the planned higher-rank
   angular expansion that lifts this restriction.

Relationship to the CP flat-source sphere solver
--------------------------------------------------

The CP flat-source method for the sphere
(:func:`~orpheus.cp.solver.solve_cp` on ``sph1D`` meshes;
:mod:`orpheus.derivations.continuous.flat_source_cp.sphere`) integrates the native
:math:`e^{-\tau}` 3-D kernel analytically over each concentric
shell to produce the :math:`E_3` second-difference formula
(:eq:`second-diff-sph` and :eq:`self-sph` above). The Peierls
reference **bypasses that integration** entirely:

- the kernel is :math:`e^{-\tau}`, not :math:`E_3`,
- the spatial representation is a piecewise polynomial of degree
  :math:`p - 1` per panel, not a piecewise constant,
- the ray integration is performed numerically in observer-centred
  polar coordinates, not analytically over annular shell
  boundaries.

So the two methods share **almost nothing** except the underlying
point kernel they both derive from. A sign error, off-by-one, or
factor-of-two in the :math:`E_3(\tau_a) - E_3(\tau_b) - E_3(\tau_c)
+ E_3(\tau_d)` sphere CP second-difference formula — which would
cancel between the CP solver and a CP-self-verification test —
would be caught by the Peierls reference. Conversely, a systematic
error in the :math:`e^{-\tau}` evaluation (unlikely given that it
is ``numpy.exp``) would be caught by the CP eigenvalue tests. The
two together triangulate the spherical integral-transport stack.

The **kernel-level** relationship is exactly parallel to the
cylinder case:

.. list-table:: Peierls vs CP kernel pairing by geometry
   :header-rows: 1
   :widths: 20 30 30

   * - Geometry
     - Peierls (pointwise)
     - CP (flat-source)
   * - Slab
     - :math:`E_1`
     - :math:`E_3` (2nd diff.)
   * - Cylinder
     - :math:`\mathrm{Ki}_1`
     - :math:`\mathrm{Ki}_3` (2nd diff.)
   * - Sphere
     - :math:`e^{-\tau}`
     - :math:`E_3` (2nd diff.)

The Phase-A tests pin this relationship numerically:

- ``TestPeierlsSphereSelfConvergence.test_k_eff_cauchy_convergence_at_thick_R``
  and ``test_flux_cauchy_convergence_at_thick_R`` — Peierls is
  Cauchy-convergent in :math:`k_{\rm eff}` and flux profile under
  :math:`(n_\theta, n_\rho, n_\phi)` refinement.
- ``TestCPvsPeierlsSphereAtThickR.test_k_eff_agrees_at_thick_R`` —
  CP :math:`k_{\rm eff}` and Peierls :math:`k_{\rm eff}` agree to
  < 2 % at :math:`R = 10` MFP.
- ``TestCPvsPeierlsSphereAtThickR.test_flux_shape_agrees_at_thick_R``
  — volume-weighted normalised flux profiles agree to L2
  < 5 % at :math:`R = 10` MFP.

Phase B will unify ``cp_sphere.py`` and ``cp_cylinder.py`` under a
single ``cp_geometry.py`` module (Issue #107 / N6), mirroring the
already-completed Peierls unification; the rank-1 white-BC
closure parity between CP flat-source and Peierls rank-1 at the
thick-:math:`R` limit, verified by the tests above, is the
correctness gate for that unification.

Verification evidence
---------------------

Three classes of independent checks gate the Peierls sphere
implementation: geometry primitives (L0), prefactor / kernel
normalisation / row-sum (L0), and eigenvalue / flux-shape (L1).
All 35 sphere tests pass; cylinder regression (31 tests) passes
unchanged.

.. list-table:: Spherical Peierls verification summary
   :header-rows: 1
   :widths: 40 18 20 22

   * - Check
     - Level
     - Tolerance
     - Identity / eq.
   * - Geometry constants (prefactor, sin θ, r², R²)
     - L0
     - exact
     - :eq:`peierls-sphere-polar`
   * - :math:`\rho_{\max}` closed forms (5 cases)
     - L0
     - :math:`10^{-12}`
     - :eq:`peierls-sphere-rho-max`
   * - :math:`r'` closed forms (3 cases)
     - L0
     - :math:`10^{-12}`
     - :eq:`peierls-sphere-r-prime`
   * - Optical-depth walker (4 cases)
     - L0
     - :math:`10^{-12}`
     - :eq:`peierls-sphere-ray-optical-depth`
   * - Composite radial GL (3 integrals)
     - L0
     - :math:`10^{-12}`
     - :math:`\int_0^R 4\pi r^{2}\,\mathrm dr`
   * - :math:`G_{\rm bc}` vacuum limit
     - L0
     - :math:`10^{-5}`
     - :eq:`peierls-sphere-G-bc`
   * - Row sum (homogeneous, :math:`R = 10` MFP)
     - L0
     - :math:`<10^{-3}` interior
     - :eq:`peierls-sphere-nystrom`
   * - Row sum (white-BC, :math:`R = 10` MFP)
     - L0
     - :math:`<2\times 10^{-2}`
     - :eq:`peierls-sphere-row-sum-identity`
   * - Thick vacuum limit (:math:`R = 30` MFP)
     - L1
     - :math:`<10^{-2}` vs :math:`k_\infty`
     - vacuum fixed point
   * - Vacuum :math:`k_{\rm eff}(R)` monotonicity
     - L1
     - monotone on 5 :math:`R` values
     - vacuum fixed point
   * - Quadrature convergence (:math:`R = 4` MFP)
     - L1
     - :math:`|\Delta k|_{\rm next}` < prev
     - :eq:`peierls-sphere-nystrom`
   * - Rank-1 white-BC error scan
     - L1
     - :math:`<35\%` at :math:`R=1`, :math:`<1\%` at :math:`R=10`
     - :eq:`peierls-sphere-G-bc`
   * - CP-vs-Peierls :math:`k_{\rm eff}` (:math:`R = 10` MFP)
     - L1
     - :math:`<2\%`
     - :eq:`second-diff-sph`
   * - CP-vs-Peierls flux shape (:math:`R = 10` MFP)
     - L1
     - :math:`L^2 < 5\%`
     - :eq:`second-diff-sph`

:cite:`CaseZweifel1967` tabulates bare-sphere critical-radius
:math:`R_c` values as a function of :math:`c = \nu\Sigma_f /
\Sigma_a` (1-group) and offers a literature tie-point analogous to
the cylinder's Sanchez 1982 tie-point. The Peierls-sphere test
suite currently pins the solver empirically via the vacuum-BC
thick-limit and the monotone-:math:`R` scan rather than against
the Case–Zweifel table directly — Cardinal Rule L4 forbids
hand-transcription. Programmatic ingestion of the table (and the
parametrized ``TestCaseZweifelTiePoint`` test class it gates) is
tracked under `Issue #143
<https://github.com/deOliveira-R/ORPHEUS/issues/143>`_.

Numerical cost
--------------

The :math:`(\theta, \rho)` tensor-product quadrature dominates. For
each observer :math:`r_i` and each :math:`\theta_k`, the kernel
assembly evaluates :math:`e^{-\tau}` at :math:`n_\rho` points —
which is a single ``numpy.exp`` call per sample (no special-function
recurrence, unlike the cylinder's :math:`\mathrm{Ki}_1`). Dominant
cost: :math:`O(N \cdot n_\theta \cdot n_\rho)` exponential
evaluations per group. For :math:`N = 10` radial nodes,
:math:`(n_\theta, n_\rho) = (24, 24)`, ``dps = 20``, kernel
assembly takes :math:`\approx 1` s on current hardware (cheaper
than the cylinder by the :math:`\mathrm{Ki}_1`-vs-``exp`` speed
ratio); eigenvalue power iteration is a further :math:`O(N^{3})`
LU per iteration, typically converging in 20–30 iterations to
:math:`10^{-10}` eigenvalue tolerance.

Short-circuit: the homogeneous single-shell branch of
:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.compute_G_bc` bypasses
the multi-annulus walker and computes :math:`\tau_{\rm surf} =
\Sigma_t\,\rho_{\max}` directly, making the bare-sphere case
:math:`\sim 2\times` faster than the multi-region path.

.. seealso::

   :mod:`orpheus.derivations.continuous.peierls_nystrom.sphere` — thin facade; the
   sphere-specific API names, the
   ``_build_peierls_sphere_case`` registry builder, and the
   ``continuous_cases`` registration.

   :mod:`orpheus.derivations.continuous.peierls_nystrom.geometry` — unified polar-form
   Nyström infrastructure; the
   :class:`~orpheus.derivations.continuous.peierls_nystrom.geometry.CurvilinearGeometry`
   class with ``kind = "sphere-1d"`` and
   ``kind = "cylinder-1d"`` handles both geometries through one
   code path.

   :class:`~orpheus.derivations.continuous.peierls_nystrom.geometry.PeierlsSolution`
   — canonical result container; same dataclass shape for
   ``geometry_kind="sphere-1d"`` / ``"cylinder-1d"`` / ``"slab"``.

   :func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.solve_peierls_1g`
   — 1-group vacuum- or white-BC eigenvalue driver. Pass
   ``geometry=_pg.SPHERE_1D`` (or ``CYLINDER_1D`` / ``SLAB_POLAR_1D``)
   and ``boundary="vacuum"`` for the scaffold-level verification
   gate.

   ``tests/derivations/test_peierls_sphere_geometry.py`` — 17 L0
   tests for angular/radial geometry primitives, the composite
   Gauss–Legendre builder, the :math:`\rho_{\max}` closed forms,
   the :math:`r'` closed forms, and the multi-annulus
   optical-depth walker.

   ``tests/derivations/test_peierls_sphere_prefactor.py`` — 6 L0
   tests for the row-sum identity (homogeneous and
   white-BC-corrected), and the :math:`G_{\rm bc}` vacuum-limit
   sanity check.

   ``tests/derivations/test_peierls_sphere_eigenvalue.py`` — 4 L1
   tests: vacuum-BC thick limit, :math:`k_{\rm eff}(R)`
   monotonicity, quadrature convergence, white-BC thick-limit
   sanity.

   ``tests/derivations/test_peierls_sphere_white_bc.py`` — 4 L1
   tests pinning the rank-1 closure error at :math:`R \in \{1, 2,
   5, 10\}` MFP (Issue #103 bounds).

   ``tests/cp/test_peierls_sphere_flux.py`` — 4 L1 tests for
   Peierls self-convergence and CP-vs-Peierls flux / eigenvalue
   agreement at :math:`R = 10` MFP.

   :doc:`/theory/references/peierls_nystrom` — cross-cutting architectural page:
   §2 dimensionally-reduced kernels, §3 Jacobian cancellation,
   §8 white-BC rank-1 closure and Issue #100 historical record.


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
     - :ref:`theory-trajectory-resolvent` (Green's function)
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
       :func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.solve_peierls_mg`
       (Issue #104); shipped registry rows for slab + hollow
       cyl/sph 2G
     - Production via
       :func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere_mg`
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

:cite:`PomraningSiewert1982` Eq. (21) is the **structurally-independent** vacuum-sphere
reference for both families. Its derivation path
(integrate-over-:math:`\mu` then add half-spaces) is genuinely
different from the Sanchez 1986 cosh-even-extension path; PS-1982
itself confirmed via method-of-characteristics gives a third
independent confirmation of the same kernel. The local PDF copy is
``1982NSE80-481.pdf`` in the repo root (commercial title:
*JQSRT* 28(6), 503-506 (1982)).

ORPHEUS implementation:
:func:`orpheus.derivations.continuous.peierls_nystrom.ps1982_reference.solve_ps1982_vacuum_sphere`
— Nyström quadrature on the
:math:`[E_1(|r-x|) - E_1(r+x)]` kernel via
:func:`mpmath.expint`. Used as the L1 cross-check truth source for
Variant α vacuum BC; also serves as a candidate L1 reference for the
Nyström family's vacuum sphere closure when high precision is needed.

The Sanchez 1986 / 2002 architecture references
------------------------------------------------

:cite:`SanchezTTSP1986` is the load-bearing literature for both families.
Its Appendix gives the Green's function :math:`t(r',\Omega'\to r,
\Omega)` for sphere with the full
:math:`(\alpha,\beta)`-parametrised BC; the angle-integrated reduction
:math:`g_\alpha(\rho'\to\rho) = \int t\,\mathrm d\Omega` is what the
Nyström family discretises (Phase 4 matrix-Galerkin) and what the
Green's function family *avoids* discretising (Variant α works with
:math:`t` directly). Both perspectives are correct math; the
discretisation choice is what makes them different.

:cite:`Sanchez2002` extends the trajectory architecture to *lattice*
geometries via the periodic-trajectory closure
:math:`\psi = \psi_q(L)/(1-\psi_{bd}(L))\cdot\psi_{bd} + \psi_q`
(Eq. 15) — algebraically parallel to Variant α's
:math:`\psi_{\rm surf} = T(\mu_{\rm surf})\,B(\mu_{\rm surf})`
bounce sum, but for a different geometry. Cross-check on the
universality of the multi-bounce factor.

Hébert 2009 and Stamm'ler 1983
-------------------------------

:cite:`Hebert2020` Chapter 3 is the textbook reference for the
**rank-1** white-BC closure :math:`(1-P_{ss})^{-1}` (§3.8.5;
Eq. 3.323 = :math:numref:`hebert-3-323`) and the rank-2 per-face
F.4 closure (§3.8.4). :cite:`Stamm1983` Chapter 4 gives the same closure
in a different notation (Stamm'ler Eq. 34 = Hébert 3.323). Variant α's
V_α2 algebraic identity :math:`T_{00}^{\rm sphere} =
P_{ss}^{\rm sphere}` (Eq. :math:numref:`peierls-greens-V-alpha-2` in
:ref:`theory-trajectory-resolvent`) explains *why* Variant α at rank-1
agrees bit-for-bit with the existing Phase 4
``boundary="specular_multibounce"`` at :math:`N=1` and with
``boundary="white_hebert"`` rank-1: all three reduce to the same
closed-form geometric series.

:cite:`Stamm1983` Chapter 6 is the architectural reference for the
**multi-group** Peierls form — group-wise iteration with
cross-group + fission source. The same form is implemented in
both ORPHEUS families
(:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.solve_peierls_mg`
for the Nyström family;
:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere_mg`
for the Green's function family).

The Garcia 2018 / 2020 / 2021 stable :math:`P_N` family
--------------------------------------------------------

Garcia's stable-:math:`P_N` family of papers (:cite:`Garcia2020`,
:cite:`Garcia2021`; the 2018 JCTT predecessor that establishes the
exterior-of-sphere stability) is the only modern external numerical
reference in the published literature for sphere homogeneous +
multi-region transport with reflective BC. Garcia 2021 specifically
covers the multi-region sphere with internal sources and provides
4-significant-figure converged scalar-flux profiles — the L1
flux-shape reference for Variant α multi-region (Plan-(b) Option 1;
documented at :ref:`theory-trajectory-resolvent`). Garcia 2021 is
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

- ``Pomraning Siewert (1982) On the integral form of the equation of transfer for a homogeneous sphere.pdf`` — Pomraning-Siewert 1982 (vacuum sphere
  integral form). Local copy.
- ``Hebert(2009)Chapter3.pdf`` — Hébert 2009 Ch. 3 (CP, integral
  transport, white-BC closures). Local copy.
- ``Sanchez(2002) Treatment of Boundary Conditions in Trajectory-Based Deterministic Transport Methods.pdf`` — Sanchez 2002 (periodic-trajectory closure,
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
  Equation by Collision Probability Methods.pdf`` — Sanchez 1977
  (anisotropic 2D CP). Local copy.
- ``Collision Probabilities for Finite Cylinders and Cuboids.pdf`` —
  finite-axial-extent cylinder + cuboid CP. Local copy.
- ``Integral form of the equation of transfer for a homogeneous
  sphere with linearly anisotropic scattering.pdf`` — Sanchez 1986
  precursor (linearly anisotropic scattering). Local copy.
- Sanchez 1986 :cite:`SanchezTTSP1986` — primary; PDF in repo root as
  ``Sanchez(1986)TTSP14.pdf``.
- Garcia 2018 / 2020 / 2021 :cite:`Garcia2020` :cite:`Garcia2021` — paywalled;
  no local copy.
- Mitsis 1963 (ANL-6787) — pseudo-slab equivalence for vacuum sphere;
  numerical reference. No local copy.
- Sood-Forster-Parsons 2003 (LA-13511) — analytical benchmark
  test set (critical sphere c_crit(R) tables). No local copy.

Codebase pointers:

- :mod:`orpheus.derivations.continuous.peierls_nystrom.geometry` — Nyström
  family unified solver
  (:func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.solve_peierls_1g`,
  :func:`~orpheus.derivations.continuous.peierls_nystrom.geometry.solve_peierls_mg`).
- :mod:`orpheus.derivations.continuous.trajectory_resolvent.greens_function` —
  Green's function family
  (:func:`~orpheus.derivations.continuous.trajectory_resolvent.greens_function.solve_greens_function_sphere`,
  ``_mg``, ``_mr``, ``_mr_fixed_source``).
- :mod:`orpheus.derivations.continuous.peierls_nystrom.ps1982_reference` —
  PS-1982 reference solver (vacuum sphere only).
- :mod:`orpheus.derivations.continuous.trajectory_resolvent.origins.specular.greens_function`
  — SymPy derivations V_α1, V_α2, V_α3 for the Green's function
  family.
- :file:`docs/theory/_peierls_nystrom_capability_matrix.inc.rst` —
  auto-generated from
  :func:`orpheus.derivations.continuous.peierls_nystrom.cases.capability_rows`;
  rebuilt every Sphinx build by
  :mod:`tools.verification.generate_capability_matrices`.


.. |times| unicode:: U+00D7
