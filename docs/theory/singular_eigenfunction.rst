.. _theory-singular-eigenfunction:

================================================================================
Singular Eigenfunction Expansion (Case 1960 family)
================================================================================

.. contents:: Contents
   :local:
   :depth: 2


Key Facts
=========

**Read this before modifying any solver in
:mod:`orpheus.derivations.continuous.singular_eigenfunction`.**

- **What this is**: a Pillar-2 reference family realising the angular
  Green's function :math:`G(\tau, \tau'; \mu, \mu')` constructed via
  Case ν-spectrum + half-range completeness. Built originally by
  K.M. Case 1960 [Case1960]_; extended to spheres by Mitsis 1963
  [Mitsis1963]_; extended to cylinders by Westfall–Metcalf 1973
  [WestfallMetcalf1973]_; extended to linearly-anisotropic
  reflected slabs and spheres by Atalay 1997 [Atalay1997]_.
- **Pillar classification**:

  * **Closed-form** for the criticality determinants (Atalay Eqs 46
    / 54 reduce to scalar arctan-equation roots in :math:`d`).
  * **Semi-analytical** for the cylinder Mitsis-WM Fredholm
    coupling (WM-72 Eqs 30-32 with mpmath/scipy quadrature) and
    for the linear-anisotropy machinery (X-function Eq 40,
    extrapolated endpoint Eq 42, K_j moments).

  See ``vv-principles`` § "The three pillars of verification" for
  what each pillar can and cannot prove. **MMS does not enter this
  family** — every reference here is an analytical / semi-analytical
  claim about the Boltzmann eigenvalue or eigenmode.
- **Geometry × anisotropy coverage**: cylinder isotropic
  (Westfall-Metcalf 1973); slab + sphere linearly anisotropic with
  Atalay 1997 parity flip; foundational primitives (X-function, ν₀
  dispersion root, half-range projections) shared across geometries
  via :mod:`...singular_eigenfunction.core`.
- **Critical numerical facts** (see :ref:`theory-se-errata` and
  :ref:`theory-se-precision-floor`):

  * **ERR-037** (μ = tanh(t) endpoint substitution, fixed
    2026-05-03): Atalay Eq 42 :math:`z_0` evaluator went from 1.5–2 %
    error to 6–7 digits via a 1-line variable substitution. Pinned
    by :func:`tests.derivations.test_case_method_z0`.
  * **ERR-038** (Atalay paper precision floor, characterised
    2026-05-03): the R=0.99 perfect-reflector cases stall at ~5 %
    because Atalay's Eq 46 is itself a *first-order* Fredholm
    approximation. This is documented in Atalay's text but was
    misdiagnosed twice as a bug. Honest tolerance: 7e-2 at
    R=0.99.
  * **WM-72 Eq 17 + q-formula corrections** (2026-05-03): the
    SymPy algebra-of-record discipline caught two errata in the
    primary literature transcription that survived three solver
    iterations.
- **Cross-references**: :ref:`theory-fn-method` is the structurally-
  independent collocation cross-check; :ref:`theory-galerkin-spectral`
  is the matrix-Galerkin cross-check; :ref:`theory-trajectory-resolvent`
  realises the same physics on the same Sood family via bouncing
  characteristics. :ref:`theory-sood-registry` is the truth-value
  catalogue.


Why one consolidated package
================================================================================

This page (and the underlying
:mod:`orpheus.derivations.continuous.singular_eigenfunction` package)
merges what were historically separate ``case_method/`` (Atalay 1997
slab + sphere with linear anisotropy) and ``singular_eigenfunction/``
(Westfall–Metcalf 1973 cylinder, isotropic) folders/pages. The merge is
justified because **all primary authors of the constituent papers
explicitly call this technique "Case singular eigenfunction expansion"**
and cite Case 1960 [Case1960]_ as the foundational source:

* **Case 1960** [Case1960]_ — the foundational source. Introduces the
  *elementary solutions* (a discrete pair :math:`\pm\nu_0` plus a
  continuum on :math:`(-1, 1)`) that subsequent authors named the
  "Case singular eigenfunctions" after.

* **Mitsis 1963 (ANL-6787)** [Mitsis1963]_, abstract:
  *"Transport solutions to the monoenergetic plane, spherical, and
  cylindrical critical problems with isotropic scattering are developed
  by the method of singular expansion modes."*

* **Westfall–Metcalf 1973** [WestfallMetcalf1973]_, introduction:
  *"Since the introduction of the singular eigenfunction expansion
  technique by Case [Ref. 1] in 1960, a wide variety of transport
  problems have been treated by this method."*

* **Atalay 1997** [Atalay1997]_, abstract: *"Case's singular
  eigenfunction method is used to formulate the criticality conditions.
  In addition to available bi-orthogonality relations in the literature,
  some parallel relations are derived to obtain the solution."*

The two ORPHEUS folders represented **two parametric variations on the
same method**, not two methods:

* :mod:`...singular_eigenfunction.cylinder` (Westfall–Metcalf 1973
  family): cylinder, **isotropic scattering**, modified-Bessel-K
  radial kernels via the addition theorem + full-range completeness
  theorem.
* :mod:`...singular_eigenfunction.slab` and
  :mod:`...singular_eigenfunction.sphere` (Atalay 1997 family): slab
  + sphere, **linear anisotropy**, half-range bi-orthogonality +
  first-order Fredholm iteration.

Different geometry, different scattering anisotropy, different kernel
reduction (exponential :math:`E_n` 's vs Bessel :math:`K_0` 's), but
**identical mathematical machinery** above the trusted-library line:
the discrete Case eigenvalue :math:`\nu_0` (same dispersion function),
continuum modes on :math:`(-1, 1)`, half-range orthogonality, Fredholm
reduction of the boundary condition. The consolidation also cures the
previous ``case_method/`` folder's violation of the project no-author-
folders rule.

The dispersion function

.. math::
   :label: case-dispersion-function

   \Lambda(\nu) = 1 - c\nu\,\mathrm{atanh}(1/\nu) ,
   \qquad \nu \in \mathbb{C}

is shared with :mod:`...fn_method` via
:func:`...fn_method.core.dispersion.case_nu0` — a *medium* property,
acceptable cross-package reuse below the trusted-library line. Above
that line the F_N method and this package are structurally independent
(see :ref:`theory-se-vs-fn`).

.. (vv-status rationale) governing: Case dispersion function — defines the discrete eigenvalue ν_0 of the Boltzmann operator at multiplying media (c > 1). Verified by V_se-cyl.1 (test_v_se_cyl_1_dispersion_function) and re-derived inside both fn_method and singular_eigenfunction packages.
.. vv-status: case-dispersion-function documented

.. _theory-se-vs-fn:
.. _theory-case-vs-fn-method:

Why a separate package from the F_N method
================================================================================

The F_N method (Siewert–Benoist 1979 + Grandjean–Siewert 1979 for
slab; Siewert–Thomas 1986 for sphere; ``fn_method/``) and this
package's singular-eigenfunction methods (Atalay 1997 slab/sphere;
WM-72 cylinder) belong to the same broad "Case 1960 / Mitsis 1963 /
McCormick–Kušcer 1965" family, but the structural mathematics is
fundamentally different.

* **F_N (slab/sphere).** Imposes the exact half-space exit-distribution
  equation at :math:`(N+1)` collocation points and solves a
  determinantal equation :math:`\det M(R) = 0`. Slab and sphere are
  unified by a geometry sign :math:`s = \pm 1`; the scalar flux is
  recovered post-hoc via the KLL 1974 Fredholm iteration on
  :math:`A(\nu)` (see :ref:`theory-fn-method`).

* **Atalay (slab/sphere).** Uses the **full normal-mode expansion**
  (discrete + continuum) and reduces the boundary-condition closure
  to an arctan-equation criticality condition (Eqs 46 / 54). The two
  methods access the same Wiener-Hopf X-function but exercise it
  differently: F_N uses :math:`X(\xi_\beta)` at collocation points
  :math:`\xi_\beta`; Atalay uses :math:`X(\pm\nu_0)` and
  :math:`X(-\nu)` for :math:`\nu \in (0, 1)` integrated against
  half-range moments (the K_j / L_j integrals).

* **Westfall–Metcalf (cylinder).** WM 1973 explicitly notes:
  *"we found that the solution as formulated by Mitsis is not
  convergent"* for the bare cylinder. WM-72's reformulation uses
  Bessel-kernel expansion + iterative Fredholm scheme — different
  mathematical machinery again. The F_N method is **not applicable**
  to the cylinder geometry without major modification, and Mitsis-
  style Wiener-Hopf is documented as non-convergent there. Hence
  cylinder lives in singular_eigenfunction; the F_N cylinder slot
  in :ref:`theory-fn-method` is a permanent stub flagging the
  out-of-pillar status.

The two packages are **structurally independent** above the trusted-
library line (sharing only ``numpy``, ``scipy.integrate.quad``,
``mpmath``).

.. _spectrum-and-cases-expansion-theorem:

The Case spectrum and the expansion theorem
================================================================================

This section is the singular-eigenfunction pillar's math-rich centre.
It walks the reader from the transport equation's eigenvalue problem
to Case's expansion theorem, half-range completeness, and the
geometry-specific reductions that the :class:`Spectrum` class
encapsulates. The class
(:class:`orpheus.derivations.continuous.singular_eigenfunction.Spectrum`)
is the math-heart of this pillar — a single named *concept* that
covers slab + sphere + cylinder, isotropic + linearly-anisotropic, on
a unified factory contract. The exposition below ties the class's
docstring to the literature line by line; reading them together is
the intended workflow.

The transport operator's eigenvalue problem
--------------------------------------------------------------------------------

For the homogeneous one-speed transport equation in a multiplying
medium with isotropic scattering,

.. math::
   :label: spectrum-transport-equation

   \mu\,\partial_x \psi(x, \mu)
   + \Sigma_t \psi(x, \mu)
   = c\,\Sigma_t\,\bar{\psi}(x),
   \qquad
   \bar{\psi}(x) = \frac{1}{2}\int_{-1}^{1} \psi(x, \mu')\,d\mu',

seek separated solutions of the form
:math:`\psi(x, \mu) = \phi_\nu(\mu)\,e^{-\Sigma_t x / \nu}` with
unknown eigenvalue :math:`\nu`. Substituting and dividing by the
exponential gives the **Case-eigenfunction equation**:

.. math::
   :label: spectrum-case-eigenfunction-equation

   (\nu - \mu)\,\phi_\nu(\mu)
   = c\,\nu \int_{-1}^{1} \phi_\nu(\mu')\,d\mu'.

The function :math:`\phi_\nu(\mu)` is the **singular eigenfunction**
of the transport operator at eigenvalue :math:`\nu`. The right-hand
side has no :math:`\mu` dependence — only an integral that we
normalise without loss of generality to unity, leaving

.. math::
   :label: spectrum-case-eigenfunction-explicit

   \phi_\nu(\mu) = \frac{c\,\nu}{2 (\nu - \mu)}.

This formula is *almost* well-defined: it has a pole at
:math:`\mu = \nu`. For :math:`\nu \notin [-1, 1]` the pole is
outside the integration domain and the formula is well-behaved
everywhere. For :math:`\nu \in [-1, 1]` the pole is *inside* the
domain and the formula is singular — hence the name *singular
eigenfunction*. The correct interpretation in the singular case
involves principal-value distributions and a delta-function residue,
treated below.

The discrete spectrum: dispersion relation
--------------------------------------------------------------------------------

Outside the continuum :math:`\nu \notin [-1, 1]`, the formula
Eq. :eq:`spectrum-case-eigenfunction-explicit` integrates to

.. math::

   \int_{-1}^{1} \phi_\nu(\mu)\,d\mu
   = \frac{c\,\nu}{2}\,\ln\!\left(\frac{\nu + 1}{\nu - 1}\right)
   = c\,\nu\,\mathrm{atanh}(1/\nu).

Imposing the normalisation
:math:`\int_{-1}^{1} \phi_\nu(\mu')\,d\mu' = 1` (which our scaling
of Eq. :eq:`spectrum-case-eigenfunction-explicit` already encoded)
gives the **Case dispersion relation**:

.. math::
   :label: spectrum-dispersion-relation

   \Lambda(\nu) := 1 - c\,\nu\,\mathrm{atanh}(1/\nu) = 0.

This equation has solutions :math:`\nu = \pm\nu_0` (the
**discrete eigenvalues** — a positive/negative pair by the symmetry
of :math:`\Lambda`). The character of :math:`\nu_0` depends on
:math:`c`:

* **Sub-multiplying** (:math:`c < 1`): :math:`\nu_0 > 1` is
  *real* and the modes :math:`\phi_{\pm\nu_0}(\mu)\,
  e^{-x/\nu_0}` represent diffusing-and-decaying neutrons.
* **Super-multiplying** (:math:`c > 1`): :math:`\nu_0 = i u_0`
  is *purely imaginary* and the modes oscillate in
  :math:`x` instead of decaying, capturing the criticality
  configuration where neutron production balances leakage.
* **Critical** (:math:`c = 1`): :math:`\nu_0` diverges to
  :math:`\pm\infty`; the discrete pair coalesces with the
  continuum's edge.

The dispersion function
:math:`\Lambda(\nu) = 1 - c\nu\,\mathrm{atanh}(1/\nu)` is a
**medium property** — it depends on :math:`c` alone, not on the
geometry. The same :math:`\nu_0` appears in the slab, sphere,
and cylinder problems for the same medium. This is *the* load-bearing
fact behind the cross-package shared dispersion-root primitive
:func:`...fn_method.core.dispersion.case_nu0`: F_N and
singular_eigenfunction read :math:`\nu_0` from the same numerical
function below the trusted-library line because they are computing
the same mathematical object. Above the trusted-library line, the
two methods are structurally independent (different collocation /
projection / iteration schemes), but the *spectrum itself* is
shared. **The eigenvalue is a property of the medium, not the
geometry.**

V_se-cyl.1 verifies this geometry-independence claim algebraically
in SymPy.

The continuum spectrum: principal value plus delta
--------------------------------------------------------------------------------

For :math:`\nu \in (-1, 1)`, the formal expression
:math:`\phi_\nu(\mu) = c\nu/(2(\nu - \mu))` has a pole inside the
integration domain and the eigenfunction equation
:eq:`spectrum-case-eigenfunction-equation` cannot be satisfied by a
function — the formal solution is a **distribution**. The correct
treatment, due to Case (1960) §3, is

.. math::
   :label: spectrum-continuum-eigenfunction

   \phi_\nu(\mu) = \mathrm{P}\,\frac{c\,\nu}{2 (\nu - \mu)}
                  + \lambda(\nu)\,\delta(\nu - \mu),

where :math:`\mathrm{P}` denotes the Cauchy principal value (the
distribution that integrates against test functions by deleting an
:math:`\epsilon`-neighbourhood of the singularity and taking
:math:`\epsilon \to 0`), and the dispersion function on the
continuum,

.. math::

   \lambda(\nu) = 1 - c\,\nu\,\mathrm{atanh}(\nu),

is the residue weight that makes the eigenfunction equation
:eq:`spectrum-case-eigenfunction-equation` hold *as a
distribution*. Both pieces are essential — neither the principal
value alone nor the delta alone solves the equation; only the
combination does.

The continuum eigenvalues form a **continuous parameter family**
indexed by :math:`\nu \in [-1, 1]`. The full spectrum of the
transport operator is then

.. math::
   :label: spectrum-full-decomposition

   \Sigma_{\rm transport} = \{\pm\nu_0\} \cup [-1, 1].

This is *the spectrum*. Different geometries select different
*subsets* of this spectrum at boundaries, but all three
geometries (slab, sphere, cylinder) work on the *same* spectrum.
This is what justifies the unified :class:`Spectrum` class — it
encapsulates the spectrum :math:`\Sigma`, and dispatches to the
geometry-specific projection that pins the expansion coefficients.

The fact that the continuum mode has a non-zero weight
:math:`\lambda(\nu)` at the diagonal :math:`\nu = \mu` is what
makes naive Gauss-Legendre evaluation of the cylinder Fredholm
kernel **silently saturate at the wrong value** (this is exactly
the Numerical Bug Signatures §6 ERR-036 fingerprint applied to
the singular-eigenfunction pillar). The Mitsis-Zweifel
singular-subtraction identity (V_se-cyl.8, Eq.
:eq:`wm72-singular-subtraction`) absorbs the principal-value +
delta residue into a single regular integral plus a residue
term — the load-bearing structural identity behind the
cylinder's :func:`...cylinder.one_group._build_M_A_phi`
diagonal evaluation.

The expansion theorem (Case 1960; Mika 1961)
--------------------------------------------------------------------------------

The load-bearing claim that justifies *using* the Case spectrum is
the **expansion theorem**:

  **Theorem (Case 1960, Theorem 1; completeness in Mika 1961).**
  Any admissible angular flux on the homogeneous transport
  equation in a finite spatial domain admits a unique expansion

  .. math::
     :label: spectrum-expansion-theorem

     \psi(x, \mu) = a_+ \phi_{\nu_0}(\mu)\,e^{-x/\nu_0}
                  + a_- \phi_{-\nu_0}(\mu)\,e^{x/\nu_0}
                  + \int_{-1}^{1} A(\nu)\,\phi_\nu(\mu)\,
                    e^{-x/\nu}\,d\nu

  with three sets of expansion coefficients:

  * :math:`a_+` (dominant outgoing diffusion mode amplitude),
  * :math:`a_-` (incoming diffusion mode amplitude — vanishes
    on a half-space, non-zero on a finite slab),
  * :math:`A(\nu)` (continuum density on :math:`[-1, 1]` —
    the analog of a Fourier amplitude for the continuum
    eigenfunctions).

Mika 1961 proves the *completeness* claim — that no admissible
:math:`\psi(x, \mu)` is missed by this expansion. The combined
result (Case 1960 expansion + Mika 1961 completeness) gives:

* **Existence**: every angular flux has at least one expansion.
* **Uniqueness**: every angular flux has at most one expansion.

The expansion is the singular-eigenfunction-pillar counterpart of
the Galerkin / Legendre projection in the F_N method (and indeed
the Galerkin projection IS a *truncation* of this expansion onto a
finite Legendre basis). Truncation introduces an
:math:`O(N^{-p})` error; the singular-eigenfunction expansion is
*exact* (up to numerical evaluation of :math:`\nu_0` from the
dispersion relation and the X-function from its integral
definition).

What the expansion theorem reduces the criticality problem to
is: *given* that any :math:`\psi` admits the expansion
:eq:`spectrum-expansion-theorem`, what conditions do the boundary
constraints impose on the coefficients :math:`(a_+, a_-, A(\nu))`?
The answer is a half-range completeness theorem (next), which
provides the projection operator that pins the expansion
coefficients from boundary data.

Half-range completeness and the X-function (Inönü 1973)
--------------------------------------------------------------------------------

At a boundary surface, the angular flux's relevant function space
splits into incoming (:math:`\mu > 0`) and outgoing (:math:`\mu < 0`)
**half-ranges**. The half-range completeness theorem (Case-Zweifel
1967 Ch. 4; Inönü 1973 for finite media) states that any half-range
function admits a unique expansion onto

.. math::

   \{X(\mu)^{-1}\,\phi_{\nu_0}(\mu),\,
     X(\mu)^{-1}\,\phi_\nu(\mu) : \nu \in (0, 1)\},

where :math:`X(\mu)` is the **Wiener-Hopf X-function** of Inönü
1973:

.. math::
   :label: spectrum-x-function

   X(\mu) = \exp\!\left[
       \frac{c\mu}{2}
       \int_0^1 \frac{\mathrm{atanh}(\mu')\,d\mu'}
                       {(\mu' - \mu)\,(c\mu'\,\mathrm{atanh}(\mu') - 1)}
       \right].

The X-function is the **projection operator** that completes the
half-range basis: given boundary data on the half-range, the
expansion coefficients :math:`(a_+, a_-, A(\nu))` are determined
by a (linear) projection involving X. This is the structural
counterpart of inverting a Gram matrix in a Galerkin scheme — the
X-function is the singular-eigenfunction-pillar's Gram operator.

ORPHEUS implements :math:`X(\mu)` in
:func:`orpheus.derivations.continuous.fn_method.core.x_function.x_function_atalay`
with a critical numerical detail: the integrand has an algebraic
pole at :math:`\mu' = 1` cancelled by an opposing
:math:`(c\mu'\,\mathrm{atanh}(\mu') - 1)` factor, but ``mp.quad``
saturates at the wrong value if asked to compute the cancellation
directly (Numerical Bug Signatures §7, ERR-037 — the canonical
instance of the *quadrature endpoint pole-cancellation slow
convergence* fingerprint). The fix is the substitution
:math:`\mu' = \tanh(t)` mapping :math:`(0, 1) \to (0, \infty)` with
Jacobian :math:`\mathrm{sech}^2(t)` that *exactly* cancels the
pole under change-of-variables; the integrand becomes smooth at the
new endpoint :math:`t \to \infty` and standard quadrature
recovers 6-7 digits at ``mp.dps = 15``.

The X-function carries the medium dependence (the integrand depends
on :math:`c`); it does NOT depend on the geometry. Slab, sphere, and
cylinder all use the same X-function — only the *projection*
they apply at the boundary differs.

Geometry reductions: slab, sphere, cylinder
--------------------------------------------------------------------------------

The expansion theorem + half-range completeness give the **machinery**;
the geometry supplies the **boundary conditions** that pin the
expansion coefficients :math:`(a_+, a_-, A(\nu))`. Each
:class:`Spectrum` geometry instance is one such reduction. The
*spectrum itself* is the same; the *projection* differs.

**Slab** (Atalay 1997 Eq. 46). For a symmetric slab of half-thickness
:math:`d` with linearly-anisotropic scattering :math:`f_1` and
specular-style reflection :math:`R \in [0, 1)` at both faces, the
even-mode criticality condition reduces to a scalar arctan-equation
root in :math:`d`:

.. math::

   \tan\!\Big(\pm\frac{\pi}{2} - \theta_{\rm LHS}^{(46)}\Big)
   = \theta_{\rm RHS}^{(46)}

where :math:`\theta_{\rm LHS}^{(46)}` and :math:`\theta_{\rm RHS}^{(46)}`
are explicit functions of :math:`(c, R, f_1, d, \nu_0, \bar\nu, z_0,
K_0, K_1, K_2)` (see :ref:`theory-case-slab-eq46` below for the full
form). The K-moments :math:`K_j(c, R, d) = \int_0^1 \mu^j\,T(R, \mu,
d)\,H(\mu)\,d\mu` are half-range integrals against the Atalay
:math:`T(R, \mu, d)` kernel and the half-range projection
:math:`H(\mu)`.

**Sphere** (Atalay 1997 Eq. 54, derived via Mitsis 1963 parity flip).
The antisymmetric BC :math:`\psi(x, \mu) = -\psi(-x, -\mu)` reduces
the sphere problem to the *odd-mode* counterpart of the slab problem
on :math:`[-R, R]`. The structural changes from slab to sphere are
surgical:

* Kernel :math:`T(R, \mu, d) \to T_1(R, \mu, d)` — a sign flip in the
  second exponential of the T numerator and denominator.
* K-moments :math:`K_j \to L_j` — same kernel structure with the
  T → T_1 substitution.
* LHS criticality term: :math:`\sin \leftrightarrow \cos` shuffle and
  :math:`R \to -R` in the reflection-coefficient terms.

**The discrete eigenvalue :math:`\nu_0` and the continuum on
:math:`[-1, 1]` are identical** — the same :class:`Spectrum`
instance applies, only the boundary projection differs.

**Cylinder** (Westfall-Metcalf 1972, isotropic only). The cylindrical
addition theorem brings Bessel-K kernels into play, introducing
*additional* spectral structure beyond the basic Case eigenfunctions:
the Bessel modes :math:`I_0(R/\nu)` couple radial and angular
dimensions in a way slab/sphere geometries don't. Resolution needs
the **Mitsis-WM Fredholm iteration** (Eqs.
:eq:`wm72-eq30-bare` / :eq:`wm72-eq31` / :eq:`wm72-eq32`) with
**Mitsis-Zweifel singular subtraction** (V_se-cyl.8) handling the
principal-value + delta residue cleanly.

Linear anisotropy (Atalay 1997)
--------------------------------------------------------------------------------

Replacing isotropic scattering :math:`c\,\bar{\psi}(x)` by linearly
anisotropic :math:`c\int (1 + 3 f_1\,\mu\mu')\,\psi(x, \mu')\,d\mu'`
adds a term to the dispersion relation but **preserves the
spectrum's structure**: still one discrete pair :math:`\pm\nu_0`
plus a continuum on :math:`[-1, 1]`. The new dispersion function
includes :math:`f_1` corrections; the X-function picks up a
:math:`f_1`-dependent integrand. Atalay's K_j and L_j moments
integrate the new structure explicitly.

The validity range :math:`c \le 1 + 1/(3 f_1)` (Atalay Eq 5)
bounds the regime where the transport operator has *only one* pair
of discrete modes. Outside that band, complex-conjugate eigenvalue
pairs appear (Dahl-Sjöstrand 1979; Kohut 1993), and Atalay's
first-order Fredholm iteration fails to detect them. This validity
bound is the singular-eigenfunction pillar's structural envelope —
:class:`Spectrum` does not check the bound at construction (the
:math:`(c, f_1)` combination might be valid for slab but not
sphere, etc.), but the underlying Atalay solver raises if asked to
solve outside the bound.

The cylinder pillar (WM-72) is **isotropic only**. Linearly-
anisotropic cylinder is research-grade and not in the package.
:class:`Spectrum.from_problem` rejects cylinder + non-zero
:math:`f_1` at construction for this reason — surfacing the
out-of-pillar status at the facade boundary so callers know the
geometry/material combination is not shippable.

Where Spectrum sits in the three-meanings taxonomy
--------------------------------------------------------------------------------

The Sanchez-Chandrasekhar three-meanings taxonomy
(:ref:`reference-solvers-three-meanings`) classifies every
reference solver by *which mathematical object* it constructs and
*how* it constructs it. :class:`Spectrum` lives under meaning
**(γ): singular-eigenfunction angular Green's function** —
directly construct :math:`G(\tau, \tau'; \mu, \mu')` as a sum
over ν-spectrum eigenfunctions weighted by X-function residues.
It is the "spectrum-decomposing" attack on the same
boundary-value problem.

The other two meanings, on the same physical configuration,
are:

* **(α): trajectory resolvent** — :class:`Billiard` in
  :mod:`...trajectory_resolvent`. Trace bouncing characteristics
  through phase space, sum the multi-bounce series
  :math:`T = (I - S)^{-1}` (Birkhoff transfer-operator
  resolvent on the billiard table). The Green's function is the
  *path-integral* sum.

* **(β): spectral resolvent** — closed-form spectral
  μ-integration of the within-medium angular Green's function
  (Sanchez 1986 Eq. A1 / PS-1982 Eq. 21). Currently a stub in
  ORPHEUS (``spectral_resolvent/`` reserved); the closed-form
  spectral kernel is currently obtained indirectly through
  :class:`Billiard` rather than via the direct PS-1982 evaluator.

When all three constructions exist for the same problem,
agreement at all three points is L1-grade evidence per the
project's structural-independence rule: the three integrands are
*structurally* distinct (ν-spectrum vs ray-traced phase space vs
spectral μ-integration), so triple agreement is not coincidence.

Where the F_N method sits in this picture
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:class:`MomentSpace` (in :mod:`...fn_method`) is *also* under (γ):
it works in the Case eigenfunction representation but **truncates**
the basis to a finite Legendre moment expansion. The F_N method
projects boundary data onto :math:`(N+1)` collocation points and
solves :math:`\det M = 0` for the critical configuration. The
truncation is what makes F_N tractable at small N (typical
operating point :math:`N = 8`-:math:`12`); it is also what makes
F_N approximate (truncation error :math:`O(N^{-p})`) rather than
exact.

:class:`Spectrum` and :class:`MomentSpace` therefore live on the
**same pillar (γ)** but with different closures: Spectrum carries
the full :math:`(\pm\nu_0, [-1, 1])` spectrum; MomentSpace truncates
to :math:`(N+1)` Legendre coefficients. They cross-check each
other at the geometry boundary (slab and sphere — the cylinder is
not in F_N's pillar; that's why the package separation exists),
which is *meaningful* L1 evidence: Spectrum's full-expansion
exactness combined with MomentSpace's spectral-convergence
:math:`O(N^{-p})` agreement is a constructive proof that both
pillars are correctly identifying the same critical configuration.

The Sanchez-Chandrasekhar taxonomy is what tells us that
:class:`Billiard` agreement (a different pillar) is *structurally*
stronger than :class:`MomentSpace` agreement — Billiard exercises
a different integrand (ray-traced phase space) than Spectrum's
ν-spectrum integrand. In the verification matrix, the
``Spectrum`` ↔ ``Billiard`` cross-check is the strongest L1
gate; the ``Spectrum`` ↔ ``MomentSpace`` cross-check is L1 but
*sibling* (same pillar).

The Spectrum class as math-heart
--------------------------------------------------------------------------------

The :class:`Spectrum` class
(:class:`orpheus.derivations.continuous.singular_eigenfunction.Spectrum`)
encapsulates the geometry, the materials, the BC, and the
quadrature size for a single instance of the singular-eigenfunction
attack. It is the **3rd concrete instance of the math-heart pattern**
in the project, alongside :class:`Billiard` (trajectory_resolvent;
Birkhoff transfer-operator resolvent on a billiard) and
:class:`MomentSpace` (fn_method; Galerkin half-range Legendre
projection).

All three classes:

1. Are **frozen dataclasses** that own a :class:`GeometrySpec` plus
   method-specific configuration (``fn_order`` for MomentSpace,
   ``alpha_payload`` for Billiard, ``n_modes`` for Spectrum).
2. Are constructed via a **factory** ``from_problem(materials:
   dict[int, Mixture], geometry: GeometrySpec, ...)`` that accepts
   the production-protocol input shape (the same shape the discrete
   CP/SN/MOC solvers consume).
3. Expose a **Protocol-conforming surface** — ``materials``,
   ``geometry_spec``, ``method_name`` — that the unifying
   :class:`TransportSolver` Protocol consumes for cross-method
   dispatch.
4. Return a **shared cross-method result type** —
   :class:`CriticalSolution` from ``solve_critical``,
   :class:`FluxSolution` from ``solve_fixed_source`` — making
   results substitutable at the cross-method protocol boundary.

The "≥3 instances" threshold is what *empirically validates* the
unifying Protocol. With Spectrum landing the third sibling, the
Protocol's structural design (factory shape, property surface,
shared return types) is no longer a one-off pattern from
trajectory_resolvent or a two-off from fn_method + trajectory_resolvent;
it is a pattern that survived three independent mathematical
instances. The Protocol is not posited; it is observed.

Construction examples
~~~~~~~~~~~~~~~~~~~~~

Bare cylinder Sood ``Ua-1-0-CY`` benchmark (:math:`c = 1.30`,
isotropic):

.. code-block:: python

   import numpy as np
   from orpheus.derivations.common.geometry_spec import GeometrySpec
   from orpheus.derivations.common.xs_library import make_mixture
   from orpheus.derivations.continuous.singular_eigenfunction import (
       Spectrum,
   )
   from orpheus.geometry.mesh import BC

   mix = make_mixture(
       sig_t=np.array([1.0]),
       sig_c=np.array([0.0]),
       sig_f=np.array([0.3]),
       nu=np.array([2.0]),
       chi=np.array([1.0]),
       sig_s=np.array([[0.7]]),  # c = (0.7 + 0.6)/1.0 = 1.30
   )
   from orpheus.geometry.structured_geometry import (
       Region, StructuredGeometry,
   )
   geom = StructuredGeometry(
       geometry="CYL",
       regions=(Region(mat_id=0, outer_thickness_cm=1.72500292),),
       bcs=(BC.vacuum,),
   )
   spec = Spectrum(geometry=geom, materials={0: mix}, n_modes=24)
   sol = spec.solve_critical()
   # sol.parameter_value ≈ 1.7250035 mfp
   # sol.eigenvalue == 1.0 (k_eff at criticality)

Reflected slab with linear anisotropy (Atalay 1997 Table 2,
:math:`R = 0.50`, :math:`f_1 = 0.20`):

.. code-block:: python

   mix = make_mixture(
       sig_t=np.array([1.0]),
       sig_c=np.array([0.0]),
       sig_f=np.array([0.3]),
       nu=np.array([2.0]),
       chi=np.array([1.0]),
       sig_s=np.array([[0.7]]),
       sig_s1=np.array([[0.20 * 0.7]]),  # f1 = SigS[1]/SigS[0] = 0.20
   )
   geom = StructuredGeometry(
       geometry="SLB",
       regions=(Region(mat_id=0, outer_thickness_cm=2.0),),  # placeholder full width
       bcs=(BC.vacuum, BC("partial", {"albedo": 0.50})),  # R = 0.50 outer
   )
   spec = Spectrum(geometry=geom, materials={0: mix}, n_modes=8)
   sol = spec.solve_critical()
   # sol.parameter_value: critical half-thickness in mfp

Discipline boundary
~~~~~~~~~~~~~~~~~~~

:class:`Spectrum` is a **thin facade** over the existing
function-level API in :mod:`...slab.one_group`,
:mod:`...sphere.one_group`, and :mod:`...cylinder.one_group`. The
load-bearing implementation lives at the function level; the
class is the integration point for the cross-method protocol. The
foundation tests pin **bit-equality** between the class API and
the function API (verified via ``float.hex`` exact-bit comparison)
so calling through the class introduces zero numerical drift.

Above the trusted-library line, Spectrum and Billiard are
**structurally independent** — they share only ``numpy`` /
``scipy.special`` / ``mpmath``, and the dispersion-root primitive
:func:`...fn_method.core.dispersion.case_nu0` (a medium property,
not method machinery). The cross-check between
``Spectrum.solve_critical()`` and ``Billiard.solve_critical()`` on
the same Sood case is therefore L1-grade structurally-independent
evidence per the project's ``vv-principles`` skill § "structural
independence applies above the trusted-library line" rule.

Cylinder — Westfall–Metcalf 1973, bare radially-reflected, isotropic
================================================================================

.. _theory-se-wm72-derivation:

Bare-critical cylinder via WM-72 — overview
--------------------------------------------------------------------------------

For a bare infinite cylinder with monoenergetic, isotropic scattering
and :math:`c > 1` (multiplying medium so :math:`\nu_0 = i u_0` is
purely imaginary), the WM-72 derivation reduces the problem to the
coupled Fredholm system (WM-72 Eqs 30-32, **bare-cylinder limit**
:math:`a_0 = b_0 = 1`, :math:`d_0 = 0`, :math:`D(\nu) = 0`,
:math:`A(\nu) = B(\nu)`):

.. math::
   :label: wm72-eq30-bare

   \Phi'(\mu) = \Phi'_0(\mu)
              - c \int_0^1 \frac{A'(\nu)\,\nu^2\,H(\nu, \mu)}{\nu + \mu}\,d\nu

.. math::
   :label: wm72-eq31

   A'(\nu) = \frac{1}{N_2(\nu)}
             \int_0^1 \mu^2\,\eta_{2\nu}(\mu)\,\Phi'(\mu)\,d\mu

.. math::
   :label: wm72-eq32

   g(R) := c \int_0^1 \frac{\mu^2\,\Phi'(\mu)}{\mu^2 + u_0^2}\,d\mu = 0

.. (vv-status rationale) governing: WM-72 Eq 30 in the bare-cylinder reduction (a_0=b_0=1, d_0=0, D=0); the entire Mitsis-WM Fredholm iteration is the verification (test_singular_eigenfunction_cylinder + test_wm72_table_ii_six_configurations).
.. vv-status: wm72-eq30-bare documented

.. (vv-status rationale) governing: WM-72 Eq 31 — Mitsis-Zweifel singular-subtraction recovery of A'(ν) from Φ'(μ); verified by V_se-cyl.8 (test_v_se_cyl_8_singular_subtraction_eq31).
.. vv-status: wm72-eq31 documented

.. (vv-status rationale) governing: WM-72 Eq 32 — the bare-cylinder criticality condition; root-find on g(R) = 0 gives the critical radius. Verified at L1 against Sood Ua-1-0-CY (3e-7 relative).
.. vv-status: wm72-eq32 documented

where:

* :math:`\Phi'(\mu)` is the angular-flux-related profile;
  :math:`\Phi'_0(\mu) = -I_0(R/\mu)\,q(\nu_0, \mu)\,\eta_0(\mu)` is
  the inhomogeneous source term derived from the discrete Case mode.
* :math:`A'(\nu)` is the continuum amplitude.
* :math:`N_2(\nu) = \int_0^1 \mu^2 \eta_{2\nu}^2(\mu)\,d\mu` is the
  WM-72 Eq 21d normalisation integral.
* :math:`\eta_0(\mu) = c\nu_0^2/(\nu_0^2 - \mu^2)` is the
  discrete-mode pseudo-eigenfunction (V_se-cyl.2 — *with* the typo
  correction).
* :math:`\eta_{2\nu}(\mu)` is the continuum mode.
* :math:`q(\nu, \mu) = (R/\nu)\,K_0(R/\mu)\,I_1(R/\nu)
  + (R/\mu)\,K_1(R/\mu)\,I_0(R/\nu)` is the kernel function
  (V_se-cyl.5, *with* the q-formula correction).
* :math:`H(\nu, \mu) = (\nu - \mu)^{-1}\,
  [I_0(R/\mu)\,q(\nu, \mu)/I_0(R/\nu) - 1]` is the non-singular
  function appearing in Eq 30.

The criticality condition Eq 32 is one scalar nonlinear equation in
one scalar unknown :math:`R`; Brent root-find gives :math:`R_c`.

Production solver — full Mitsis-WM Fredholm method
--------------------------------------------------------------------------------

.. _theory-se-wm72-numerics:

The Branch-2 production solver
(:func:`orpheus.derivations.continuous.singular_eigenfunction.cylinder.solve_singular_eigenfunction_cylinder_bare_critical`)
implements the full Mitsis-WM method-of-record. Two algorithmic
choices distinguish it from the original 1973 implementation:

**1. Linear system instead of Jacobi iteration.** Substituting the
discretised Eq 31 (mapping :math:`\Phi' \to A'`) into discretised
Eq 30 gives

.. math::
   :label: wm72-coupled-linear-system

   (\mathbb{I} + c\,M_{A\phi}\,M_{\phi A})\,\mathbf{A}'
   = M_{A\phi}\,\boldsymbol{\Phi}'_0

.. (vv-status rationale) derivation: Combined linear system for the WM-72 bare-cylinder Fredholm coupling — substitutes Eq 31 into Eq 30 and lets numpy.linalg.solve replace WM-72's 1973-era Jacobi iteration. Verified by the L1 critical-radius gate at all six WM-72 Table II configurations.
.. vv-status: wm72-coupled-linear-system documented

where :math:`M_{A\phi}` and :math:`M_{\phi A}` are the discretised
integral operators. Solving by ``numpy.linalg.solve`` is faster and
more accurate than Jacobi iteration; WM-72 themselves used iteration
in 1973, but the equivalent linear-system formulation is what the
modern API permits.

**2. Mitsis-Zweifel singular subtraction** for :math:`M_{A\phi}`. The
:math:`\eta_{2\nu}(\mu)` continuum mode has a Cauchy P.V. + a
:math:`\lambda(\nu)\delta(\nu-\mu)` from the dispersion relation
(WM-72 Eq 19). The singular subtraction trick is

.. math::
   :label: wm72-singular-subtraction

   \int_0^1 \mu^2\,\eta_{2\nu}(\mu)\,\Phi'(\mu)\,d\mu
   = \int_0^1 \frac{c\,\nu^2\,[\mu^2 \Phi'(\mu) - \nu^2 \Phi'(\nu)]}
                    {\nu^2 - \mu^2}\,d\mu + \nu^2\,\Phi'(\nu) ,

.. (vv-status rationale) derivation: Mitsis-Zweifel singular-subtraction identity — collapses the Cauchy P.V. + λδ continuum-mode kernel into a regular GL-quadrable integrand plus a single residue. Verified by V_se-cyl.8 (test_v_se_cyl_8_singular_subtraction_eq31).
.. vv-status: wm72-singular-subtraction documented

absorbing both the Cauchy P.V. of the regular-η₂ν part and the
:math:`\lambda(\nu)\,\delta` from Eq 19 into a single regular integral
plus a residue. The diagonal point :math:`\mu_j = \nu_i` of the
regular integrand is evaluated by **Lagrangian-interpolation
differentiation matrix** on the GL nodes — exactly the technique
WM-72 cite on p. 7: *"evaluating the derivative term by Lagrangian
interpolation over all points."*

**3. Scaled Bessel functions** (``i0e``, ``i1e``, ``k0e``, ``k1e``)
throughout to avoid overflow at large :math:`R/\mu`. The exponential
factors in :math:`I_0(R/\mu)\,K_n(R/\mu)` cancel pairwise; the
:math:`I_n(R/\nu)\,K_n(R/\mu)` products with :math:`\nu \neq \mu`
carry exponential :math:`e^{R/\nu - R/\mu}` factors that we keep in
scaled form.

Convergence achieved
~~~~~~~~~~~~~~~~~~~~~

At :math:`n_{\rm grid} = 24` (matching WM-72's original 24-GL
quadrature order), the solver reproduces every WM-72 Table II
configuration to better than 5e-7 relative:

.. list-table:: WM-72 Table II benchmark agreement
   :header-rows: 1
   :widths: 15 30 30 25

   * - ``c``
     - ``R_c`` (mfp) truth
     - Solver result
     - Relative error
   * - 1.05
     - 5.411288 (WM-72)
     - 5.4112891
     - 4e-7
   * - 1.10
     - 3.577391 (WM-72)
     - 3.5773921
     - 3e-7
   * - 1.20
     - 2.287209 (WM-72)
     - 2.2872099
     - 4e-7
   * - 1.30
     - 1.72500292 (Sood)
     - 1.72500349
     - 3e-7
   * - 1.40
     - 1.396979 (WM-72)
     - 1.39697910
     - 5e-8
   * - 2.00
     - 0.668613 (WM-72)
     - 0.66861305
     - 8e-8

Wall-clock time per solve: ≤ 0.1 s on a typical container CPU.

Convergence rate is **near-spectral** for smooth integrands (the GL
quadrature + barycentric differentiation give the expected
exponential-rate convergence on smooth analytic functions).
Empirical: error at :math:`n=12` is ≤ 1e-5; at :math:`n=24` ≤ 1e-6;
at :math:`n=48` ≤ 1e-7. This is **4-6 orders of magnitude better**
than the original Phase B1 prototype (which used direct Nyström on
WM-72 Eq 6a + single-cell product integration on the log-singular
kernel diagonal, achieving only :math:`O(1/n)` algebraic convergence
with a 1e-3 floor at :math:`n=128`).

Test gates: :mod:`tests.derivations.test_singular_eigenfunction_cylinder`
(one foundation test per ``derive_*()`` in the cylinder origins
module, plus production-solver sanity tests, an L1 Sood
``Ua-1-0-CY`` reference-value gate at 1e-5, and a parametrized L1
gate over all six WM-72 Table II configurations).

V_se-cyl — eight SymPy verifications
-------------------------------------

The cylinder Branch-1 SymPy module
(:mod:`...origins.cylinder_derivations`) ships eight foundation-tagged
``derive_*()`` functions that pin the WM-72 algebra letter-for-letter.
Each is one verifiable structural identity. Together they form the
**fully-verified algebra-of-record** for the WM-72 derivation.

.. _theory-se-V-se-cyl-1:

V_se-cyl.1 — Dispersion function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`orpheus.derivations.continuous.singular_eigenfunction.origins.cylinder_derivations.derive_dispersion_function`.
**Test gate:**
:func:`tests.derivations.test_singular_eigenfunction_cylinder.test_v_se_cyl_1_dispersion_function`.

The dispersion function :eq:`case-dispersion-function` is identical
to slab/sphere — reflecting the *medium-property nature* of the
Case singular eigenfunctions. A discrete eigenvalue :math:`\nu_0`
satisfies :math:`\Lambda(\nu_0) = 0`; for :math:`c > 1` the root is
purely imaginary, :math:`\nu_0 = i u_0` with :math:`u_0 > 0`. SymPy
verifies the dispersion-function form against WM-72 Eq 18 and against
the Case 1960 form (relabeling :math:`c` ↔ Case's notation).

.. _theory-se-V-se-cyl-2:

V_se-cyl.2 — Discrete pseudo-eigenfunction (catches WM-72 Eq 17 typo)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`...origins.cylinder_derivations.derive_discrete_pseudo_eigenfunction`.
**Test gate:**
:func:`tests.derivations.test_singular_eigenfunction_cylinder.test_v_se_cyl_2_discrete_pseudo_eigenfunction`.

The discrete pseudo-eigenfunction is

.. math::

   \eta_0(\mu) = \frac{c\,\nu_0^2}{\nu_0^2 - \mu^2} ,

which satisfies WM-72 Eq 15 directly. **SymPy V_se-cyl.2 caught a
typo in the printed WM-72 Eq 17 on first run.** As printed,
Eq 17 reads :math:`\eta_{0l}(\mu) = c_l\,\nu_{0l}/(\nu_{0l}^2 -
\mu^2)` (single :math:`\nu_0` in numerator). Direct substitution
into Eq 15 gives:

* LHS reduces to :math:`c\nu_0\mu^2`;
* RHS evaluates to :math:`c\nu_0^2\mu^2`.

The mismatch is one power of :math:`\nu_0`. The corrected form
:math:`\eta_{0l}(\mu) = c_l\,\nu_{0l}^2/(\nu_{0l}^2 - \mu^2)`:

* Satisfies Eq 15 exactly.
* Reproduces Eq 21d's :math:`\nu_0^4` factor on
  :math:`N_0 = \int_0^1 \mu^2 \eta_0^2\,d\mu`.
* Closes the half-range normalisation Eq 14 under dispersion.

This is the canonical case for "re-derive every published equation
in SymPy" as a publication-grade discipline (see the
algebra-of-record skill).

.. _theory-se-V-se-cyl-3:

V_se-cyl.3 — Bessel-Wronskian identity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`...origins.cylinder_derivations.derive_bessel_wronskian_identity`.
**Test gate:**
:func:`tests.derivations.test_singular_eigenfunction_cylinder.test_v_se_cyl_3_bessel_wronskian_identity`.

The Bessel-Wronskian identity

.. math::
   :label: bessel-wronskian

   K_1(z)\,I_0(z) + I_1(z)\,K_0(z) = 1/z

is used in the WM-72 Eq 9 integrodifferential reduction. SymPy
verifies the identity by direct expansion (it is a textbook result
following from
:math:`I_0(z)\,K_0'(z) - I_0'(z)\,K_0(z) = -1/z` plus
:math:`I_0' = I_1`, :math:`K_0' = -K_1`).

.. (vv-status rationale) derivation: Bessel-Wronskian identity — standard mathematical identity used in WM-72 Eq 9 reduction; verified at L0 by V_se-cyl.3.
.. vv-status: bessel-wronskian documented

.. _theory-se-V-se-cyl-4:

V_se-cyl.4 — Bare-cylinder reduction (a_0 = b_0 = 1, d_0 = 0)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`...origins.cylinder_derivations.derive_bare_cylinder_reduction`.
**Test gate:**
:func:`tests.derivations.test_singular_eigenfunction_cylinder.test_v_se_cyl_4_bare_cylinder_reduction`.

For the bare cylinder (:math:`c_1 = c_2 = c`, no reflector) WM-72
Eq 27's source-prefactor term :math:`(c_2 - c_1) \to 0` cancels the
inhomogeneous part of the discrete amplitude :math:`d_0`, leaving:

* :math:`a_0 = b_0 = 1` (discrete-mode normalisation; one degree
  of freedom remains).
* :math:`d_0 = 0` (continuum-mode discrete amplitude vanishes).
* :math:`D(\nu) = 0` (continuum source vanishes).
* :math:`A(\nu) = B(\nu)` — **NOT zero**, but the two continuum
  amplitudes equal each other from Eq 33's middle-term reduction.

This last point is a **correction vs the original Phase B1 stylized
SymPy**, which assumed :math:`A = B = 0` for the bare cylinder. The
correct reduction preserves :math:`A(\nu) = B(\nu)` as the active
Fredholm unknown; setting it to zero would over-constrain the system
and give wrong critical radii.

.. _theory-se-V-se-cyl-5:

V_se-cyl.5 — Bare-cylinder criticality structure (catches q-formula typo)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`...origins.cylinder_derivations.derive_bare_cylinder_criticality_condition`.
**Test gate:**
:func:`tests.derivations.test_singular_eigenfunction_cylinder.test_v_se_cyl_5_bare_cylinder_criticality_condition`.

The corrected q-formula is

.. math::
   :label: wm72-q-formula

   q(\nu, \mu) = \frac{R}{\nu}\,K_0(R/\mu)\,I_1(R/\nu)
              + \frac{R}{\mu}\,K_1(R/\mu)\,I_0(R/\nu) .

.. (vv-status rationale) derivation: WM-72 q-formula — the second-term denominator must be μ (not R) to satisfy the Wronskian identity q(μ, μ) = 1; corrected from original Phase B1 SymPy. Verified by V_se-cyl.5 (test_v_se_cyl_5_bare_cylinder_criticality_condition).
.. vv-status: wm72-q-formula documented

The :math:`(R/\mu)` denominator in the second term is *forced* by
the Wronskian identity :math:`q(\mu, \mu) = 1` — which is itself
structurally required by Eq 29b's :math:`\nu \to \mu` limit. The
original Phase B1 SymPy module wrote :math:`q(\nu, \mu) =
(R/\nu)\,K_0\,I_1 + R\,K_1\,I_0` (no :math:`\mu` denominator on the
second term), giving :math:`q(\mu, \mu) \approx 0.72` numerically
— inconsistent with Eq 29b. The corrected form (now in V_se-cyl.5)
closes the Wronskian identity exactly:

.. math::

   q(\mu, \mu) = \frac{R}{\mu}\,[K_0(R/\mu)\,I_1(R/\mu)
                                 + K_1(R/\mu)\,I_0(R/\mu)]
              = \frac{R}{\mu}\cdot\frac{1}{R/\mu}
              = 1 .

This was caught by the post-hardening cross-check between V_se-cyl.5's
algebraic claim and the Wronskian identity V_se-cyl.3 — exactly the
discipline the algebra-of-record skill recommends.

.. _theory-se-V-se-cyl-6:

V_se-cyl.6 — Discrete eigenfunction normalisation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`...origins.cylinder_derivations.derive_discrete_eigenfunction_normalization`.
**Test gate:**
:func:`tests.derivations.test_singular_eigenfunction_cylinder.test_v_se_cyl_6_discrete_eigenfunction_normalization`.

The discrete normalisation :math:`N_0 = \int_0^1 \mu^2 \eta_0^2(\mu)
\,d\mu` matches WM-72 Eq 21d. With the corrected V_se-cyl.2
:math:`\eta_0(\mu) = c\nu_0^2/(\nu_0^2 - \mu^2)`, the integral
evaluates to

.. math::

   N_0 = \frac{c^2\,\nu_0^4}{2}\!\left[\frac{1}{1 - \nu_0^2}
                                       + \frac{1 - 2\nu_0^2}{2\nu_0^2}\,
                                         \log\!\frac{1 + \nu_0}{\nu_0 - 1}
                                       \right]

(WM-72 Eq 21d, verbatim) — and the :math:`\nu_0^4` prefactor is
the cross-check that V_se-cyl.2 carries the right power of
:math:`\nu_0`.

.. _theory-se-V-se-cyl-7:

V_se-cyl.7 — Bare-cylinder flux reconstruction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`...origins.cylinder_derivations.derive_flux_reconstruction_bare_cylinder`.
**Test gate:**
:func:`tests.derivations.test_singular_eigenfunction_cylinder.test_v_se_cyl_7_flux_reconstruction_bare_cylinder`.

The bare-cylinder neutron density profile is

.. math::
   :label: wm72-rho-bare-cylinder

   \rho(r) = b_0\,J_0(r/u_0)

.. (vv-status rationale) governing: Bare-cylinder neutron density profile — the dominant Case eigenfunction with imaginary ν_0 = i u_0 for c > 1; verified by V_se-cyl.7 (test_v_se_cyl_7_flux_reconstruction_bare_cylinder).
.. vv-status: wm72-rho-bare-cylinder documented

with :math:`b_0 = 1` per the bare-cylinder normalisation. SymPy
verifies that with :math:`\nu_0 = i u_0` (purely imaginary for
:math:`c > 1`), the modified-Bessel form :math:`I_0(r/(i u_0))`
reduces to the ordinary-Bessel form :math:`J_0(r/u_0)` — the
canonical "imaginary argument" Bessel-function relation.

.. _theory-se-V-se-cyl-8:

V_se-cyl.8 — Mitsis-Zweifel singular-subtraction structural identity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**SymPy derivation:**
:func:`...origins.cylinder_derivations.derive_singular_subtraction_eq31`.
**Test gate:**
:func:`tests.derivations.test_singular_eigenfunction_cylinder.test_v_se_cyl_8_singular_subtraction_eq31`.

V_se-cyl.8 is the load-bearing algebra behind the Branch-2 production
solver's diagonal handling. The Mitsis-Zweifel identity
:eq:`wm72-singular-subtraction` collapses the PV residue and
:math:`\lambda \delta` continuum-mode kernel via the dispersion
identity :math:`(1-\lambda) + \lambda = 1`. SymPy verifies the
identity by expanding both sides on a worked
:math:`\eta_{2\nu}(\mu) = c\nu/(\mu - \nu) + \lambda(\nu)\delta(\mu -
\nu)` decomposition.

This identity is what makes the Branch-2 production solver's diagonal
treatment work cleanly: the regular integrand
:math:`c\,\nu^2\,[\mu^2 \Phi'(\mu) - \nu^2 \Phi'(\nu)]/(\nu^2 - \mu^2)`
is finite and smooth at :math:`\mu = \nu` (L'Hôpital limit), and the
:math:`\nu^2 \Phi'(\nu)` residue is evaluated exactly. Without
V_se-cyl.8, the diagonal entry would have to be approximated by
some ad-hoc nearest-neighbour rule, which would degrade the
near-spectral convergence rate.

.. _theory-se-wm72-xverif:

Cross-check vs Variant α at Sood ``Ua-1-0-CY``
--------------------------------------------------------------------------------

The hardened WM-72 solver agrees with the Variant α cylinder solver
at the Sood ``Ua-1-0-CY`` configuration to ≤ 1e-5 relative. Both
solvers reproduce the published Sood :math:`r_c = 1.72500292` mfp:

* **Variant α** at 8.5e-6 (already shipped at
  :func:`tests.derivations.test_trajectory_resolvent_cylinder_xverif_sood2003`,
  via bouncing-characteristic integration with analytical
  bounce-period summation).
* **WM-72** at ≤ 3e-7 (this module, via singular-eigenfunction
  Fredholm coupling with Mitsis-Zweifel subtraction).

Both methods share **only** the dispersion-root primitive
(``case_nu0``) — a medium property independent of geometry. Above
the trusted-library line, the methods are entirely disjoint:

* **Variant α**: angle-resolved scalar transport in
  :math:`(r, \mu, \phi)` phase space, with bouncing characteristics
  and analytical bounce-period summation.
* **WM-72**: scalar-density integral transport equation reduced via
  singular-eigenfunction expansion, with Bessel-kernel matrix
  Fredholm coupling between :math:`\Phi'(\mu)` and :math:`A'(\nu)`.

Per ``algebra-of-record`` § "Structural independence applies above
the trusted-library line", the cross-check at the same precision
level (1e-5) is a true structurally-independent L1 anchor for the
Sood reference value. A third leg via ``peierls_nystrom``
(Bickley-Naylor :math:`\mathrm{Ki}_3` integrals) is available for
future expansion.

Test gate:
:mod:`tests.derivations.test_singular_eigenfunction_cylinder_xverif`.

Slab — Atalay 1997 reflected, linearly anisotropic
================================================================================

Atalay's 1997 paper extends the singular-eigenfunction expansion to
linearly anisotropic scattering with reflective BC. The paper's
**load-bearing technical contribution** is the four NEW parallel
half-range relations Eqs 28-31 that close the deficit blocking
previous reflected-slab attempts.

.. _theory-case-half-range:

The four new parallel half-range relations
--------------------------------------------------------------------------------

The McCormick–Kušcer 1965 bi-orthogonality relations
[McCormickKuscer1965]_ (Atalay Eqs 18-21) integrate the weight

.. math::

   [\phi_{0+}(\mu) + B c\nu_0/2]\,\gamma(\mu)\,(\nu_0 - \mu)

against the four half-range basis members
:math:`\{1, \mu, \mu^2, \mu^3\}`. These relations suffice for the
half-space Milne problem but **NOT** for the reflected-slab
boundary condition Atalay Eq 16, which requires a second weight

.. math::

   [\phi_\nu(\mu) + c\nu/(2\bar\nu)]\,\gamma(\mu) ,

evaluated against the same four basis members. Atalay's Eqs 28-31
fill the gap.

The structural parallelism is verified in
:func:`...origins.slab_sphere_derivations.derive_atalay_half_range_eqs28_to_31`
(V_case.3): Eqs 28-31 share the Wiener-Hopf X-function factor
structure (:math:`X(\pm\nu_0)`, :math:`X(\pm\nu')`, :math:`N(\nu)`)
of the McCormick–Kušcer set. The closing identity is:

.. list-table:: Half-range bi-orthogonality completion
   :header-rows: 1
   :widths: 30 35 35

   * - Weight × basis
     - McCormick-Kušcer set
     - Atalay parallel set
   * - Weight type
     - :math:`[\phi_{0+}(\mu) + Bc\nu_0/2]\gamma(\mu)(\nu_0 - \mu)`
     - :math:`[\phi_\nu(\mu) + c\nu/(2\bar\nu)]\gamma(\mu)`
   * - Atalay Eqs
     - 18, 19, 20, 21
     - 28, 29, 30, 31
   * - Closes
     - Half-space Milne
     - Reflected slab/sphere

V_case.3 verifies the structural parallelism — both sets reduce to
linear combinations of the same X-function evaluations.

.. _theory-case-slab-eq46:

Slab criticality (Eq 46)
--------------------------------------------------------------------------------

The slab criticality condition is the closed-form arctan equation
(**even-mode**, fundamental):

.. math::
   :label: singular-eigenfunction-eq46

   \pm\frac{\pi}{2} - \arctan\!\frac{R \sin[(d-z_0)/|\nu_0|] + \sin[(d+z_0)/|\nu_0|]}
                                       {R \cos[(d-z_0)/|\nu_0|] - \cos[(d+z_0)/|\nu_0|]}
   = \arctan\!\frac{\big(K_0 \bar\nu - 3 f_1 (1-c) \bar\nu
                       [K_1 \bar\nu d(\nu_0^2) - K_0 \nu_0^2 d(\bar\nu^2)]\big) |\nu_0|}
                  {(1+K_2) d(\nu_0 \bar\nu) d(-\nu_0 \bar\nu)
                   + K_1 \bar\nu d(\nu_0^2) - K_0 \nu_0^2 d(\bar\nu^2)}

**SymPy derivation:**
:func:`...origins.slab_sphere_derivations.derive_atalay_critical_slab_eq46`.
**Test gate:**
:func:`tests.derivations.test_case_method_symbolic.test_v_case_5_critical_slab_eq46`.

Eq 46 emerges from Atalay Eq 43 (the boundary-condition closure
written in :math:`(R e^{\pm i a_1}, e^{\pm i a_2})` form) by
observing that the numerator and denominator of the closure ratio
are complex conjugates — the log of their ratio is :math:`2i \arg(z)`,
which gives the arctan form.

The :math:`\pm \pi/2` ambiguity reflects multiple eigenvalue modes;
the fundamental mode is the smallest positive :math:`d`. The Branch-2
production solver brackets the mode-1 root via prominence-filtered
zero-crossing of :math:`\sin(\theta_{LHS} - \theta_{RHS})` (which
crosses zero smoothly at the critical thickness without :math:`\tan`
blowing up).

The :math:`K_j` symbols are half-range moments — see
:ref:`theory-case-x-function-eq40` for the X-function and
:ref:`theory-case-validity-eq5` for the validity bound. The
:math:`d(x) \equiv 1 - 3 f_1\,x` shorthand is Atalay's
linearly-anisotropic prefactor:

.. math::

   d(x) = 1 - 3 f_1\,x ,

which reduces to :math:`d \equiv 1` for isotropic scattering
(:math:`f_1 = 0`).

Atalay's tabulated Eq 46 critical thicknesses:

.. list-table:: Atalay 1997 Table 2 (vacuum + reflected slab,
   :math:`f_1 = 0`)
   :header-rows: 1
   :widths: 12 22 22 22 22

   * - :math:`c`
     - :math:`R = 0` (mfp)
     - :math:`R = 0.25` (mfp)
     - :math:`R = 0.50` (mfp)
     - :math:`R = 0.75` (mfp)
   * - 1.30
     - 1.87766
     - 1.40621
     - 0.89317
     - 0.40758
   * - 1.50
     - 1.30528
     - 0.97735
     - 0.62094
     - 0.28319
   * - 2.00
     - 0.66766
     - 0.50049
     - 0.31787
     - 0.14491

The Branch-2 solver
(:func:`...slab.solve_case_method_slab_critical`) reproduces these
to ≤ 5e-2 absolute, with the residual gap explained as a paper
precision floor at small thicknesses (see
:ref:`theory-se-precision-floor`).

Sphere — Atalay 1997 reflected, linearly anisotropic, parity flip
================================================================================

.. _theory-case-sphere-eq54:

Sphere criticality (Eq 54) via parity flip
--------------------------------------------------------------------------------

The sphere criticality condition is

.. math::
   :label: singular-eigenfunction-eq54

   \arctan\!\frac{\sin[(d+z_0)/|\nu_0|] - R \sin[(d-z_0)/|\nu_0|]}
                {\cos[(d+z_0)/|\nu_0|] + R \cos[(d-z_0)/|\nu_0|]}
   = \arctan\!\frac{\big(L_0 \bar\nu - 3 f_1 (1-c) \bar\nu
                       [L_1 \bar\nu d(\nu_0^2) - L_0 \nu_0^2 d(\bar\nu^2)]\big) |\nu_0|}
                  {(1+L_2) d(\nu_0 \bar\nu) d(-\nu_0 \bar\nu)
                   + L_1 \bar\nu d(\nu_0^2) - L_0 \nu_0^2 d(\bar\nu^2)}

**SymPy derivation:**
:func:`...origins.slab_sphere_derivations.derive_atalay_critical_sphere_eq54_via_parity_flip`.
**Test gate:**
:func:`tests.derivations.test_case_method_symbolic.test_v_case_6_critical_sphere_eq54_via_parity_flip`.

The sphere problem is treated as the **odd-mode** of the slab
problem on :math:`[-R, R]` via the antisymmetric BC

.. math::

   \psi(x, \mu) = -\psi(-x, -\mu)

(Atalay Eq 47 — the :math:`r\Phi(r) \to \Psi(x, \mu)` substitution
picks up an :math:`r`-sign factor, exactly as in the Mitsis 1963
sphere-as-slab derivation that underlies the Siewert-Thomas 1986
F_N sphere). The structural change reduces to:

* :math:`T(R, \mu) \to T_1(R, \mu)` — sign flip on the second
  exponential in the T-function.
* :math:`K_j \to L_j` — same integrand structure with
  :math:`T \to T_1`.
* sin↔cos shuffle in the LHS arctan argument.

**Numerical parity-flip equivalence at vacuum BC.** At :math:`R = 0`
(vacuum BC), both :math:`T(0,\mu) = T_1(0,\mu) = e^{-2d/\mu}`, so
:math:`K_j = L_j` bit-for-bit. Verified in the parity-flip test
(:mod:`tests.derivations.test_case_method_slab_sphere_parity_flip`).

For non-zero :math:`R`, the :math:`T_1` sign flip on the second
exponential makes :math:`L_j \neq K_j`; both Branch-1 (SymPy) and
Branch-2 (numpy) implementations parameterise the kernel on
:math:`{\rm geometry\_sign} \in \{+1, -1\}` exactly as F_N method's
slab/sphere unification does.

Where :math:`d` here is the sphere **radius** (interpretation: the
half-thickness of the equivalent slab equals the sphere radius),
the parity-flip is a clean :math:`s = +1 \to s = -1` change with
the corresponding :math:`L_j` rather than :math:`K_j` evaluations.

Anisotropy — linearly anisotropic scattering primitives
================================================================================

The Atalay framework is a rich apparatus of shared primitives — the
X-function, the extrapolated endpoint, the validity bound. Each is a
medium property (no geometry parameter) and is consumed by both
slab Eq 46 and sphere Eq 54.

.. _theory-case-x-function-eq40:

X-function (Atalay Eq 40)
--------------------------------------------------------------------------------

.. math::
   :label: singular-eigenfunction-eq40

   X(\mu) = \exp\!\Bigg\{ -\frac{c}{2} \int_0^1 d\nu\, g_1(c,\nu)\,
       \Big[d^2(\nu^2)\Big(1 + \frac{c\nu^2}{1-\nu^2}\Big)
            + 3 f_1 (1-c)^2 \nu^2 d(-\nu^2)\Big]\,\ln(\nu - \mu) \Bigg\}

**SymPy derivation:**
:func:`...origins.slab_sphere_derivations.derive_atalay_x_function_eq40`.
**Test gate:**
:func:`tests.derivations.test_case_method_symbolic.test_v_case_8_x_function_eq40`.

The Wiener-Hopf X-function is a medium-only quantity (depends on
:math:`c` and :math:`f_1`, no geometry parameter). SymPy verifies
the integrand structure against Atalay Eq 40 and against the
isotropic-limit closed form:

.. math::

   X^{\rm iso}(\mu) = \frac{1}{1 - \mu}\,
                       \exp\!\left[\frac{1}{\pi}\int_0^1
                            \arg \Lambda^{+}(\tau)\,
                            \frac{d\tau}{\tau - \mu}\right] .

For multi-line continued-Plemelj derivations connecting Eq 40 to
Eq 26 (the implicit X-function definition via Plemelj-Sokhotski),
see Case 1960 [Case1960]_ § IV and Case-Zweifel 1967 § 4.

.. note:: **X-function divergent integrand** (open investigation).
   Phase-2.2 of the ERR-038 cascade discovered that the Atalay Eq 40
   integrand as transcribed appears to be **logarithmically
   divergent** at :math:`\nu \to 1`. The bracket carries
   :math:`1/(1-\nu^2)` (simple pole at :math:`\nu = 1`); the
   supposed cancelling factor :math:`g_1(\nu)` only goes as
   :math:`1/\log^2(1-\nu)`, so the product
   :math:`g_1 \cdot \mathrm{bracket} \cdot \log(\nu - \mu)` decays
   only as :math:`\log(\nu - \mu)/[(1-\nu)\,\log^2(1-\nu)]`, whose
   primitive is :math:`\propto -1/\log(1-\nu)` — bounded but
   non-zero as truncation point :math:`\to` endpoint.

   Empirically, X-function values drift 6e-3 across
   ``mp.dps`` 15→60 (:math:`X(-0.55) = 0.866` at dps=15; 0.860 at
   dps=60). This is **not** convergence — it is the algorithm
   sampling closer and closer to the divergent endpoint at higher
   precision. Either:

   (a) Atalay's Eq 40 has a missing :math:`(1-\nu^2)` factor
       (typesetter loss),
   (b) the integral is interpreted as a Cauchy principal value (and
       Atalay's Tables are generated under that regularisation), or
   (c) the closed-form Case–Plazcek–Hofmann 1961 X-function for
       isotropic should be substituted.

   Phase-2.2 also showed that switching the X-function evaluator to
   a μ=tanh(t) substituted form shifts X(-0.99) by 1.2 % but only
   shifts the K_j moments by 0.3-0.6 % and 2d_crit by ≤ 0.1 %. The
   X-function tanh-fix is therefore **a robustness improvement, not
   a correctness fix** — pinned by
   :mod:`tests.derivations.test_case_method_x_function`. See
   ``.claude/agent-memory/numerics-investigator/kj_residual_xfunction_divergent_2026_05_03.md``
   for the full investigation.

.. _theory-case-extrapolated-endpoint-eq42:

Extrapolated endpoint :math:`z_0` (Atalay Eq 42) and the ERR-037 fix
---------------------------------------------------------------------

.. math::
   :label: singular-eigenfunction-eq42

   z_0 = -\frac{\nu_0}{2} \ln\!\frac{d(-\nu_0\bar\nu)}{d(\nu_0\bar\nu)}
        + \frac{c\,\nu_0}{4} \int_0^1 d\mu\, g_1(c,\mu)\,
                                \Big[d^2(\mu^2)\Big(1 + \frac{c\mu^2}{1-\mu^2}\Big)
                                     + 3 f_1 (1-c)^2 \mu^2 d(-\mu^2)\Big]\,
                                \ln\!\frac{\nu_0 + \mu}{\nu_0 - \mu}

**SymPy derivation:**
:func:`...origins.slab_sphere_derivations.derive_atalay_extrapolated_endpoint_eq42`.
**Test gate:**
:func:`tests.derivations.test_case_method_z0`.

Atalay Eq 42 is the **Milne extrapolated endpoint** for linearly-
anisotropic scattering — the distance beyond the physical surface
where the asymptotic-form scalar flux extrapolates to zero. Together
with the Eq 46 / Eq 54 criticality conditions, :math:`z_0` is a
load-bearing input.

.. _theory-se-err037:

ERR-037 — μ = tanh(t) endpoint substitution (fixed 2026-05-03)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Failure mode:** the Eq 42 integrand carries the bracket
:math:`[1 + c\mu^2/(1-\mu^2)]` with a simple pole at :math:`\mu = 1`
that is **algebraically cancelled** by the
:math:`\lambda^2(\mu) \sim \log^2(1-\mu)` growth in :math:`g_1(\mu)`.
The cancellation is mathematically exact, but **inefficient under
direct mp.quad over (0, 1)**: at ``mp.dps = 15`` the gap to Atalay
Table 1 is 3.3 %; at dps=25 it is 2.1 %; at dps=35 it is 1.6 %.
Slow monotone convergence to the wrong asymptote — the canonical
**Signature 7** fingerprint (see ``numerical-bug-signatures``
§ Signature 7).

**Fix:** the substitution :math:`\mu = \tanh(t)` maps
:math:`(0, 1) \to (0, \infty)` with Jacobian

.. math::

   1 - \mu^2 = \mathrm{sech}^2(t),\quad
   \frac{\mu^2}{1-\mu^2} = \sinh^2(t),\quad
   \tanh^{-1}(\mu) = t,\quad
   d\mu = \mathrm{sech}^2(t)\,dt ,

so the :math:`\mathrm{sech}^2(t)` Jacobian **cancels the
:math:`1/(1-\mu^2)` pole exactly**. The transformed integrand is
exponentially decaying at :math:`t \to \infty`
(:math:`g_1 \sim 1/t^2`) and ``mp.quad`` resolves it to **6-7 digits
at dps=25** in sub-millisecond wall-clock.

**Impact** (relative error on :math:`z_0`):

.. list-table:: ERR-037 ``z_0`` accuracy improvement
   :header-rows: 1
   :widths: 15 30 30 25

   * - :math:`c`
     - :math:`z_0` pre-fix
     - :math:`z_0` post-fix
     - Atalay Table 1 truth
   * - 1.10
     - 0.6373
     - 0.645971
     - 0.645971 (4e-7 rel)
   * - 1.30
     - 0.5356
     - 0.547144
     - 0.547144 (1.3e-8 rel)
   * - 1.50
     - 0.4664
     - 0.474869
     - 0.474869 (1e-6 rel)
   * - 2.00
     - 0.3520
     - 0.357551
     - 0.357551 (1e-6 rel)

Worst-case post-fix: **4e-7 relative** across all 10 Atalay Table 1
entries. Pinned at 1e-5 absolute by
:func:`tests.derivations.test_case_method_z0.test_atalay_z0_table1_isotropic`.

The fix propagates downstream: at vacuum BC :math:`R = 0`, the slab
critical thickness improves from ~1.1 % error to ~0.1 %; the sphere
:math:`R_c` improves from 0.48 % to **0.001 %** at Sood
``Ua-1-0-SP``.

**Anti-pattern caught — convention drift hypothesis falsified.** The
Wave 2-B closeout (2026-05-02) had misdiagnosed the gap as a
"Case-Zweifel completeness-sum normalisation discrepancy" with
Atalay's published form, framing the issue as a "multi-day fix"
requiring a convention-bridge investigation. The actual fix was a
1-line variable substitution. The diagnostic that resolves the
ambiguity:

* **Gap increases monotonically with dps but doesn't converge** →
  endpoint quadrature problem.
* **Gap is constant in dps** → structural / convention error.

The gap going from 3.3 % at dps=15 → 1.6 % at dps=35 was the
smoking gun pointing at the quadrature endpoint cause. See
``.claude/agent-memory/numerics-investigator/wave2_convention_drift_falsified_2026_05_03.md``
for the full investigation memo.

ERR-037 K_j extended-fix FALSIFIED
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Wave 2-B closeout speculated that the same μ=tanh(t)
substitution should fix the residual gap on K_j moments
(:func:`...core.half_range._atalay_K_or_L_moment_value`).
**Investigation tested directly: hypothesis FALSIFIED.**

* K_j outer ``scipy.quad`` converges fine (rel drift ~1e-11 between
  baseline ``epsabs=1e-10`` and tight ``epsabs=1e-13``). **NOT
  Signature 7 at the K_j level.**
* Bottleneck is X(-ν) accuracy: scipy backend has ~1e-3 relative
  error vs mpmath dps=30, propagated as :math:`X^2` into K_j.
* Atalay Eq 40's X-function integrand carries the same
  :math:`1/(1-\nu^2)` bracket pole and cancelling factor
  :math:`g_1 \sim 1/\log^2(1-\nu)`. BUT the product is
  **logarithmically divergent** at :math:`\nu \to 1` (primitive
  :math:`\propto -1/\log(1-\nu)`, bounded but truncation-dependent).
* mpmath dps 15→60: X(-0.55) drifts 0.866 → 0.860 (drift 6e-3,
  monotone, slow). Both the :math:`(0, 1)` and :math:`(0, \infty)`
  substituted forms give finite values that depend on the
  algorithm's stopping point near the endpoint.
* Implementing μ=tanh(t) in :func:`_atalay_X_function_scipy` AND
  :func:`_atalay_X_function_mpmath` shifted the gap from 4.40% →
  4.44% on R=0.75 (marginal degradation). Reverted.

**Conclusion:** this is **NOT** Signature 7. The residual gap
originates in the X-function's algorithm-dependent regularisation
of a divergent integrand. Atalay's Eq 40 is either (a) missing a
factor lost in typesetting, (b) interpreted under a specific
regularisation (Cauchy P.V. or implicit via the Eq 26 definition
path), or (c) needs the closed-form Case-Plazcek-Hofmann 1961
X-function for isotropic. Tracked as a separate investigation,
regression-pinned by
:mod:`tests.derivations.test_case_method_x_function`. See the
investigation memo
``.claude/agent-memory/numerics-investigator/kj_residual_xfunction_divergent_2026_05_03.md``.

The lesson: **don't apply Signature 7's substitution recipe blindly
to call sites that look similar without verifying the same pole
structure exists.** Phase-1 verdict on convergence-with-tolerance
fingerprint should be: "scipy.quad converges fine at the OUTER
level → look INSIDE the integrand for sub-routines whose own
convergence is the bottleneck."

.. _theory-case-validity-eq5:

Validity bound (Atalay Eq 5)
--------------------------------------------------------------------------------

.. math::
   :label: singular-eigenfunction-eq5

   c \le 1 + \frac{1}{3 f_1}

.. (vv-status rationale) governing: Atalay validity bound — for f_1 > 0, c must lie in a band; outside it, complex eigenvalue pairs appear (Dahl-Sjöstrand 1979) and Atalay's first-order Fredholm iteration breaks. Verified by V_case.9 (test_v_case_9_validity_bound_eq5).
.. vv-status: singular-eigenfunction-eq5 documented

**SymPy derivation:**
:func:`...origins.slab_sphere_derivations.derive_atalay_validity_bound_eq5`.
**Test gate:**
:func:`tests.derivations.test_case_method_symbolic.test_v_case_9_validity_bound_eq5`.

Atalay Eq 5 is the **one-pair-of-discrete-modes** range. For
:math:`f_1 = 0` (isotropic) the bound is trivial (all :math:`c`);
for :math:`f_1 = 0.30` the upper limit is :math:`c \le 19/9 \approx
2.111`. Outside this band complex-conjugate eigenvalue pairs appear
(Dahl–Sjöstrand 1979; Kohut 1993) and Atalay's first-order Fredholm
iteration breaks. The Branch-2 slab/sphere
solvers raise :class:`ValueError` if :math:`c` violates the bound,
via :func:`...anisotropy.linear.check_atalay_validity`.

.. _theory-se-precision-floor:

Atalay Eq 46 / Eq 54 first-order Fredholm precision floor (ERR-038)
====================================================================

**Status:** **PAPER LIMITATION CHARACTERISED 2026-05-03** — NOT a
code bug. This entry documents a *reference* limitation, not a
solver defect, and exists so future investigators don't re-chase
the gap as a bug.

**Mechanism:** Atalay 1997 § 2 (p. 236) derives Eq 46 from the full
Fredholm equation Eq 32 by the explicit step

   *"we here skip the zeroth order and proceed directly with the
   first order approximation. This provides us the required optimum
   accuracy. The first order approximation necessitates that we
   omit the integral term in Eq.(32)."*

Tables 2-5 are then computed from Eq 46 with first-order accuracy.
Atalay (p. 246) further states:

   *"as in the work of Kaper et al. (1974) for isotropic scattering,
   one may consider to iterate further until a better convergence
   is obtained… we expect some improvement in the accuracy
   especially for the small slab thicknesses."*

**Empirical fingerprint** (c=1.30 :math:`f_1 = 0`, cascade
data):

.. list-table:: ERR-038 — error scales with 1/d_crit
   :header-rows: 1
   :widths: 30 30 20 20

   * - 2d_atalay (mfp)
     - Our 2d (mfp)
     - Rel error
     - Regime
   * - 20.0
     - 19.99961
     - 0.002 %
     - Moderate-d
   * - 2.0
     - 2.00071
     - 0.036 %
     - Moderate-d
   * - 0.20
     - 0.20641
     - 3.2 %
     - Small-d
   * - 0.01456 (R=0.99)
     - 0.01529
     - 5.0 %
     - Small-d, perfect-reflector

The error scales monotonically with :math:`1/d_{\rm crit}` — the
signature of an omitted-higher-order term that vanishes at large
:math:`d` (where :math:`T(R, \mu) \to -1` saturates and the omitted
Fredholm integral contributes negligibly to the boundary residue)
and dominates at small :math:`d`.

**Cascade evidence ruling out alternative mechanisms:**

1. **Conditioning / cancellation.** Ruled out — K_j moments are
   dps-independent at bit-identical level across dps 15→40.
2. **Singular asymptotic** :math:`R \to 1`. Ruled out — error is
   smooth in :math:`R`; at :math:`R = 0.99`, :math:`c = 1.024` with
   :math:`2d = 0.20` (Atalay Table 6) the same 3.2 % error appears,
   demonstrating the gap is in :math:`1/d_{\rm crit}` not in
   :math:`1/(1-R)`.
3. **Different quadrature pole.** Ruled out — Phase 1.5 built a
   :math:`\mu = \tanh(t)` substituted X-function (mpmath) that gives
   1.2 % different X(-0.99) values, but Phase 2.2 confirmed this
   changes K_j by 0.3-0.6 % and 2d_crit by < 0.1 % across all 15
   Table 2 cells.
4. **X-function singular branch.** Ruled out by the same evidence.
5. **Atalay paper limitation.** Confirmed by Atalay's own text +
   the :math:`1/d_{\rm crit}` error scaling fingerprint + the
   moderate-d machine-precision self-consistency.

The diagnostic that resolves "code bug vs paper floor" is **scaling
the same physical problem to a regime where the paper's
approximation is exact** (here: large :math:`d`) and verifying
machine precision there. This is the structurally-independent
ground that lets the investigation close as a paper limitation
rather than open as a bug. Without this test, "uniform 5 % offset
everywhere" (likely solver bug) cannot be distinguished from "5 %
scaling with 1/d" (likely paper floor).

**L1 tests that pin the floor:**

* :func:`tests.derivations.test_case_method_slab.test_slab_atalay_table2_r099_first_order_floor`
  — pins the 7e-2 paper precision floor at R=0.99, with
  ``catches("ERR-038")`` so any future improvement (e.g. a
  higher-order Fredholm iteration that closes the gap) signals as
  an unexpected pass.
* :func:`tests.derivations.test_case_method_slab.test_slab_atalay_table6_moderate_d_consistency`
  — pins **1e-4 relative** agreement at :math:`2d \ge 2` mfp, the
  structurally-independent moderate-d ground.

**Lesson — reference contamination by under-reading the reference's
own caveats.** When investigating a numerical disagreement with a
published reference, **read the paper's stated approximation level
explicitly before assuming the gap is a code bug.** Atalay's text
twice explicitly states the published values are first-order
approximations with degraded precision at small slab thicknesses.
The Wave 2-B closeout had listed R=0.99 as "still 10%+ off, needs
careful analysis of the singular limit 2d→0", treating it as a
mathematical singular-limit problem requiring multi-day asymptotic
analysis. The actual issue is fully documented in Atalay's own
caveats. See
``.claude/agent-memory/numerics-investigator/atalay_r099_paper_floor_2026_05_03.md``
for the cascade memo.

V_case — eight SymPy verifications (Atalay 1997)
=================================================

The slab/sphere Branch-1 SymPy module
(:mod:`...origins.slab_sphere_derivations`) ships eight foundation-
tagged ``derive_*()`` functions covering the Atalay derivation chain.

.. _theory-se-V-case-1:

V_case.1 — Linearly-anisotropic dispersion reduction (Atalay Eqs 11→12)
------------------------------------------------------------------------

**SymPy derivation:**
:func:`...origins.slab_sphere_derivations.derive_atalay_dispersion_linear_anisotropic`.
**Test gate:**
:func:`tests.derivations.test_case_method_symbolic.test_v_case_1_atalay_dispersion`.

The linearly-anisotropic dispersion relation Eq 11 reduces
algebraically to Eq 12 — the 1-pair-of-modes form. SymPy verifies the
reduction via direct expansion + ``simplify()``.

.. _theory-se-V-case-2:

V_case.2 — Symmetry conditions Eqs 13/14 + 47-49
-------------------------------------------------

**SymPy derivation:**
:func:`...origins.slab_sphere_derivations.derive_atalay_symmetry_conditions_eq13_14_47_to_49`.
**Test gate:**
:func:`tests.derivations.test_case_method_symbolic.test_v_case_2_atalay_symmetry`.

V_case.2 verifies the slab-even / sphere-odd symmetry conditions
parametrising the parity flip. The slab problem on :math:`[-d, d]`
admits the symmetric reflection :math:`\psi(x, \mu) = \psi(-x, -\mu)`
(Atalay Eq 13/14); the sphere on :math:`[-R, R]` admits the
antisymmetric :math:`\psi(x, \mu) = -\psi(-x, -\mu)` (Atalay Eq 47).
Eqs 48-49 propagate the parity into the discrete + continuum mode
decomposition. The structural identity SymPy proves: the boundary
projection equations transform under the parity flip exactly by the
T → T_1 substitution that the Branch-2 implementation uses.

.. _theory-se-V-case-3:

V_case.3 — Half-range relations Eqs 28-31 (the load-bearing extension)
-----------------------------------------------------------------------

**SymPy derivation:**
:func:`...origins.slab_sphere_derivations.derive_atalay_half_range_eqs28_to_31`.
**Test gate:**
:func:`tests.derivations.test_case_method_symbolic.test_v_case_3_half_range`.

V_case.3 is described in :ref:`theory-case-half-range`. Atalay's
load-bearing technical contribution — the four parallel half-range
relations that close the deficit blocking previous reflected-slab
attempts. SymPy verifies the structural parallelism with the McCormick-
Kušcer 1965 set: both reduce to linear combinations of the same
X-function evaluations.

.. _theory-se-V-case-4:

V_case.4 — Fredholm form Eqs 27→32
-----------------------------------

**SymPy derivation:**
:func:`...origins.slab_sphere_derivations.derive_atalay_fredholm_form_eq27_to_eq32`.
**Test gate:**
:func:`tests.derivations.test_case_method_symbolic.test_v_case_4_fredholm_form`.

V_case.4 verifies the Fredholm-form prefactor reduction Eqs 27 →
32. The prefactor algebra is a load-bearing step in the Atalay
derivation: it converts the boundary projection (an integral
equation in the continuum amplitude :math:`A(\nu)`) into the form
that admits the first-order iteration that produces Eq 46.

.. _theory-se-V-case-5:

V_case.5 — Critical slab Eq 46
-------------------------------

**SymPy derivation:**
:func:`...origins.slab_sphere_derivations.derive_atalay_critical_slab_eq46`.
**Test gate:**
:func:`tests.derivations.test_case_method_symbolic.test_v_case_5_critical_slab_eq46`.

V_case.5 verifies the structural reduction of Atalay Eq 43 (the
boundary-condition closure in :math:`(R e^{\pm i a_1}, e^{\pm i
a_2})`-form) to the arctan form Eq 46. Numerator and denominator of
the closure ratio are complex conjugates; the log of their ratio is
:math:`2i \arg(z)`, which gives the arctan form. SymPy verifies the
algebraic identity.

.. _theory-se-V-case-6:

V_case.6 — Critical sphere Eq 54 via parity flip
-------------------------------------------------

**SymPy derivation:**
:func:`...origins.slab_sphere_derivations.derive_atalay_critical_sphere_eq54_via_parity_flip`.
**Test gate:**
:func:`tests.derivations.test_case_method_symbolic.test_v_case_6_critical_sphere_eq54_via_parity_flip`
plus the numerical parity-flip equivalence test
:mod:`tests.derivations.test_case_method_slab_sphere_parity_flip`.

V_case.6 verifies that the sphere Eq 54 is the slab Eq 46 under the
parity flip :math:`T \to T_1`, :math:`K_j \to L_j`, sin↔cos shuffle.
The numerical complement at vacuum BC (where :math:`T = T_1` so
:math:`K_j = L_j`) confirms bit-for-bit agreement between the slab
solver at :math:`s=+1` and the sphere solver at :math:`s=-1`.

.. _theory-se-V-case-7:

V_case.7 — Extrapolated endpoint Eq 42 (post-fix)
--------------------------------------------------

**SymPy derivation:**
:func:`...origins.slab_sphere_derivations.derive_atalay_extrapolated_endpoint_eq42`.
**Test gate:** see :ref:`theory-case-extrapolated-endpoint-eq42`
(the ERR-037 narrative).

V_case.7 verifies the **integrand structure** of Eq 42 against
Atalay's published form. The L1 numerical gate against Atalay Table
1 lives at :func:`tests.derivations.test_case_method_z0`; after
the ERR-037 fix the gate is at 1e-5 absolute. **V_case.7 alone was
insufficient as L1 verification** — the integrand structure was
correct but the quadrature evaluator was unconverged. This is the
canonical case for "L0 derivative-level verification + L1
value-level numerical verification" — both are needed.

.. _theory-se-V-case-8:

V_case.8 — X-function Eq 40
----------------------------

See :ref:`theory-case-x-function-eq40`.

.. _theory-se-V-case-9:

V_case.9 — Validity bound Eq 5
-------------------------------

See :ref:`theory-case-validity-eq5`.

Origins — Branch-1 SymPy algebra-of-record
================================================================================

The :mod:`...singular_eigenfunction.origins` subpackage hosts both
SymPy modules:

* :mod:`...singular_eigenfunction.origins.cylinder_derivations` —
  Westfall–Metcalf 1973 cylinder derivations: bare-cylinder integral-
  equation reduction, dispersion function (re-derivation matching
  Case 1960 / WM-72 Eq 18), Bessel-Wronskian identity used in the
  integrodifferential reduction (Eq 9), pseudo-eigenfunction
  structure (Eq 17 + 19, **with typo correction**), bare-cylinder
  closure (Eq 32 + 30 simplification, **with q-formula correction**),
  Mitsis-Zweifel singular subtraction (Eq 31).

* :mod:`...singular_eigenfunction.origins.slab_sphere_derivations` —
  Atalay 1997 slab + sphere derivations: dispersion relation
  reduction Eq 11 → 12 (linearly anisotropic), parity flip Eqs 47-49
  mapping slab-odd → sphere, half-range relations Eqs 28-31, Fredholm-
  form prefactor Eq 32, criticality conditions Eqs 46 (slab) / 54
  (sphere), X-function Eq 40, extrapolated endpoint Eq 42, validity
  bound Eq 5.

.. _theory-case-sood-coverage:

Sood case coverage (Atalay)
--------------------------------------------------------------------------------

The Atalay-anchored case catalogue lives in
:mod:`orpheus.derivations.continuous.sood_registry.atalay1997` and
covers the **reflected + linearly-anisotropic cross-product cases**
that lie outside both the Sood/Forster/Parsons LA-13511 truth set
(which focuses on bare configurations) and the Burkart-Ishiguro-
Siewert 1976 F_N reference [BurkartIshiguroSiewert1976]_ (vacuum-only). Specifically, Atalay
tabulates:

* :math:`(c, R, f_1)` triples for :math:`R \in \{0, 0.25, 0.50,
  0.75, 0.99\}` and :math:`f_1 \in \{0, 0.10, 0.20, 0.30\}` (slab,
  even modes, Tables 2-5).
* :math:`f_1 = 0.10` only (sphere, odd modes, Table 10).

The current ORPHEUS Atalay case catalogue ships 7 cases (6 slab + 1
sphere); see :doc:`sood_registry` for the full list.

.. _theory-se-errata:

Errata caught (V&V publication artifacts)
================================================================================

This implementation slice caught **two errata in the primary
literature**, both discovered by the algebra-of-record SymPy
re-derivation discipline. Documenting them here as a publication
artifact of the V&V process — both errata persisted across multiple
prior solver iterations until the SymPy discipline caught them.

WM-72 Eq 17 — :math:`\nu_0` should be :math:`\nu_0^2`
-------------------------------------------------------

WM-72 Eq 17 (printed): :math:`\eta_{0l}(\mu) = c_l\,\nu_{0l}/
(\nu_{0l}^2 - \mu^2)`. Direct substitution into Eq 15 fails:

* LHS reduces to :math:`c\nu_0\mu^2`;
* RHS reads :math:`c\nu_0^2\mu^2`.

Mismatched power of :math:`\nu_0`. The correct form is
:math:`\eta_{0l}(\mu) = c_l\,\nu_{0l}^2/(\nu_{0l}^2 - \mu^2)`,
which:

* Satisfies Eq 15 exactly.
* Reproduces Eq 21d's :math:`\nu_0^4` factor on
  :math:`N_0 = \int_0^1 \mu^2\eta_0^2\,d\mu`.
* Closes the half-range normalisation Eq 14 under dispersion.

V_se-cyl.2 caught the typo on first SymPy run.

q-formula in Phase B1 SymPy — :math:`R` should be :math:`R/\mu`
----------------------------------------------------------------

The Phase B1 SymPy module wrote :math:`q(\nu, \mu) = (R/\nu)\,K_0\,
I_1 + R\,K_1\,I_0`, but the WM-72 paper (Eq 28 footnote) actually
uses :math:`(R/\mu)\,K_1\,I_0` — the :math:`R` in the second term
should be divided by :math:`\mu`. The Wronskian identity
:math:`q(\mu, \mu) = 1` (which is structurally required by
Eq 29b's :math:`\nu \to \mu` limit) FORCES the :math:`(R/\mu)`
denominator. The original B1 form gave :math:`q(\mu, \mu) \approx
0.72`, inconsistent with Eq 29b. The corrected form (now in
V_se-cyl.5) closes the Wronskian identity exactly.

This was caught by the post-hardening cross-check between
V_se-cyl.5's algebraic claim and the Wronskian identity V_se-cyl.3
— a load-bearing demonstration that **cross-checking every derived
identity against an independent structural relation is mandatory**,
not optional.

The discipline of re-deriving every published equation in SymPy
and the discipline of cross-checking every derived identity against
the Wronskian (or any other independent structural relation) are
now both empirically validated as load-bearing for V&V correctness.

References
==========

.. [Case1960] Case, K. M. (1960). "Elementary Solutions of the
   Transport Equation and Their Applications." *Annals of Physics*
   **9**, 1-23.

.. [Mitsis1963] Mitsis, G. F. (1963). "Transport Solutions to the
   Monoenergetic Critical Problems." Argonne National Laboratory
   report ANL-6787.

.. [WestfallMetcalf1973] Westfall, R. M. & Metcalf, D. R. (1973).
   "Singular Eigenfunction Solution of the Monoenergetic Neutron
   Transport Equation for Finite Radially Reflected Critical
   Cylinders." *Nuclear Science and Engineering* **52**, 1-11.
   DOI 10.13182/NSE73-A23285.

.. [Atalay1997] Atalay, M.A. (1997).
   "The reflected slab and sphere criticality problem with anisotropic
   scattering in one-speed neutron transport theory."
   *Progress in Nuclear Energy* **31**\ (3), 229-252.
   DOI: 10.1016/0149-1970(95)00094-1.

.. [McCormickKuscer1965] McCormick, N.J., Kušcer, I. (1965).
   "Bi-orthogonality relations for solving half-space transport problems."
   *J. Math. Phys.* **6**, 1939.

.. [BurkartIshiguroSiewert1976] Burkart, A.R., Ishiguro, Y., Siewert, C.E. (1976).
   "Neutron transport in two dissimilar media with anisotropic scattering."
   *Nuclear Science and Engineering* **61**, 72-81.

* **Metcalf & Zweifel 1968** — *Nucl. Sci. Eng.* **33**, 318. — the
  singular-subtraction technique used in WM-72 Eqs 31 and 33.
* **Atkinson 1976** — SIAM, *A Survey of Numerical Methods for the
  Solution of Fredholm Integral Equations of the Second Kind*,
  Chapter 6 (product integration).
* **Berrut & Trefethen 2004** — *SIAM Review* **46**, 501-517,
  "Barycentric Lagrange Interpolation." — the differentiation
  matrix used for the Eq 31 diagonal handling.
* **Kaper, Lindeman, Leaf 1974** — *NSE* **54**, 94. — the
  Wiener-Hopf Fredholm iteration baseline; used in
  :ref:`theory-fn-method` for slab + sphere interior flux
  reconstruction.
* **Dahl, E.B. & Sjöstrand, N.G. 1979** — *Nuclear Science and
  Engineering* **69**, 114. — Eigenvalue spectrum of multiplying
  slabs and spheres for monoenergetic neutrons with anisotropic
  scattering. Identifies the regime outside Atalay Eq 5 where
  complex-conjugate eigenvalue pairs appear.

Internal references:

* Method-implementer closeouts:
  ``.claude/agent-memory/method-implementer/wave_2b_atalay_case_method.md``
  (Atalay slab/sphere shipment),
  ``.claude/agent-memory/method-implementer/wm72_cylinder_phase_b1.md``
  (cylinder shipment).
* Numerics-investigator memos:
  ``.claude/agent-memory/numerics-investigator/wave2_convention_drift_falsified_2026_05_03.md``
  (ERR-037),
  ``.claude/agent-memory/numerics-investigator/atalay_r099_paper_floor_2026_05_03.md``
  (ERR-038),
  ``.claude/agent-memory/numerics-investigator/kj_residual_xfunction_divergent_2026_05_03.md``
  (X-function divergent integrand).
* :doc:`fn_method` — companion F_N collocation reference, sharing
  the Wiener-Hopf X-function below the trusted-library line.
* :doc:`trajectory_resolvent` — Variant α reference family on the
  same Sood ``Ua-1-0-CY`` truth value (cylinder cross-check).
* :doc:`sood_registry` — Sood + Atalay case catalogue.
