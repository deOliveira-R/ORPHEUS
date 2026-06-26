.. _galerkin-projection:

==========================================================
Galerkin Projection — the discrete frame and its consumers
==========================================================

Every method in ORPHEUS that transitions a function between a
**fine** representation and a **coarse** one — discrete-ordinate
angular flux to spherical-harmonic moments, fine-energy fluxes to
broad-group cross-sections, regional flux to homogenised cross-
sections — does so via a **(reconstruction, analysis)** pair of
linear operators :math:`(R, M)`. In harmonic-analysis language
this pair is exactly the two operational **faces of a discrete
frame**: the analysis operator :math:`M = T` (sampled values →
coefficients) and the reconstruction :math:`R` (coefficients →
values, the canonical-dual synthesis). This page is the canonical
home for that pair, for the discrete-**frame** abstraction that
realises it (:class:`~orpheus.numerics.frame.FrameBase` and its
discipline subclasses), for the **discipline** (Galerkin vs
Petrov-Galerkin — test space equal to vs different from trial space)
that distinguishes its variants, and for the **cross-method consumer
catalog** of where the same primitives are used in the codebase.

The user's architectural rule that motivated this page:

  *"Galerkin projection is generally useful — same machinery is
  used in cross-section energy condensation. Catalog everything
  that needs to be implemented and where."*

Binding a :class:`~orpheus.numerics.basis.Basis` to a
:class:`~orpheus.numerics.measure.DiscreteMeasure` through a single
frame (Frame/Basis carve, ``refactor/operator-inverse-algebra``;
discipline-type hierarchy, Issue #268) puts one mechanism in front
of every consumer, instead of a separate projection / reconstruction
operator class per method.

.. contents::
   :local:
   :depth: 2


.. note:: **What changed (2026-06-24, Issue #268).** Two earlier
   framings on this page were **reversed** by the P1 discipline-type
   carve (``refactor/operator-inverse-algebra``), and the page now
   reflects the shipped architecture:

   1. **Discipline is a TYPE, not a property and not an operator
      marker.** The earlier draft alternately (a) carried the
      discipline as marker ABCs ``GalerkinProjection`` /
      ``PetrovGalerkinProjection`` *on the operator role* and (b)
      proposed collapsing it to a derived *property* of the frame
      (``measure == basis.canonical_measure``). Both are retired. The
      discipline is now a genuine **kind of frame**, carried by the
      frame TYPE:
      :class:`~orpheus.numerics.frame.FrameBase` →
      :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` →
      :class:`~orpheus.numerics.frame.GalerkinFrame`. The
      :mod:`orpheus.numerics.projection` module keeps only the two
      abstract operator **roles**
      (:class:`~orpheus.numerics.projection.AnalysisOperator`,
      :class:`~orpheus.numerics.projection.ReconstructionOperator`),
      which the frame faces subclass.

   2. **Homogenisation and condensation are Petrov-Galerkin.** An
      intermediate draft argued they were "Galerkin in the natural
      :math:`L^2(\phi V)` (resp. spectrum) metric". That reading folds
      the solution (the flux :math:`\phi`) into the *metric* — it is
      legitimate **only** for forward-flux, reaction-rate-only
      reduction, and breaks under the eigenvalue-consistent
      (adjoint-weighted) homogenisation reactor physics requires
      (where the test weight is the adjoint :math:`\varphi^*`, not the
      forward flux). The solution-weighting therefore lives on the
      **test side = a distinct test basis = the frame type**, never as
      a weight smuggled onto the measure: **the measure carries the
      axis and the fixed** :math:`L^2` **metric, never the
      discipline**.

   The homogenization derivation — the headline Petrov-Galerkin
   consumer of this hierarchy — lives in
   :ref:`sn-homogenization-petrov-galerkin-frame`
   (:doc:`discrete_ordinates`); it was rewritten to this same
   Petrov-Galerkin framing under Issue #268 (the earlier
   ":math:`L^2(\phi V)`-Galerkin" reading is retired there, with the
   forward-flux metric-fold shown to be the Galerkin *degenerate* of the
   eigenvalue-consistent adjoint-weighted case).


Key Facts
=========

- A **frame** binds a :class:`~orpheus.numerics.basis.Basis` (the
  synthesis / trial side — the functions :math:`\{e_k\}` and their
  convention) to a :class:`~orpheus.numerics.measure.DiscreteMeasure`
  (the domain — the sample points and the fixed quadrature-weight
  :math:`L^2` metric). The
  :class:`~orpheus.numerics.frame.FrameBase` so formed emits two
  operational faces:

  * the **analysis** face :math:`M = T` — sampled values →
    coefficients (``frame.analysis``), measured against the **test**
    basis;
  * the **reconstruction** face :math:`R` — coefficients → values,
    the canonical-dual synthesis (``frame.reconstruction``), purely
    trial-side.

  Together :math:`(R, M)` define a Galerkin-style discretisation of
  any :math:`A : V \to V` as :math:`A_h = M A R : W \to W`
  (Brenner & Scott 2008, §3.4).

- The **discipline** — whether the *test* functions equal the *trial*
  functions — is carried by the frame **type**, a genuine
  Liskov-correct hierarchy (Issue #268):

  .. code-block:: text

     FrameBase                 abstract; the discipline-FREE mechanics
     │                         (table, spaces, reconstruction face, the
     │                         analysis-face wiring — none depend on the
     │                         test side)
     └─ PetrovGalerkinFrame    explicit TEST basis (test ≠ trial); the
        │                      analysis measures against test functions
        │                      χ_k that need NOT equal the trial φ_k, so
        │                      M* ≠ R (an oblique dual)
        └─ GalerkinFrame       test IS trial — STRENGTHENS the promise to
                               Π* = R (a symmetric, here-diagonal Gram).
                               The angular spherical-harmonic projection
                               is the canonical pure-Galerkin frame.

  :class:`~orpheus.numerics.frame.GalerkinFrame` *is-a*
  :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` with
  ``test is trial``: it strengthens (never weakens) the base promise.

- **The measure never carries the discipline.** The
  :class:`~orpheus.numerics.measure.DiscreteMeasure` carries the axis
  and a fixed :math:`L^2` metric (the quadrature weights). The
  solution-weighting (forward flux :math:`\phi`, adjoint
  :math:`\varphi^*`) that distinguishes a Petrov-Galerkin instance is
  a first-class **test basis** — the test *space* — not a metric on
  the measure. This is the load-bearing rule the homogenisation /
  condensation consumers obey.

- The :mod:`orpheus.numerics.projection` module carries only the two
  abstract operator **roles**:
  :class:`~orpheus.numerics.projection.AnalysisOperator`
  (:math:`M : V \to W`) and
  :class:`~orpheus.numerics.projection.ReconstructionOperator`
  (:math:`R : W \to V`), which the two frame faces subclass. The
  discipline is the frame's type, never a marker on these roles.

- The single concrete frame shipping today is the
  **spherical-harmonic frame**
  :meth:`Quadrature.angular_frame(L)
  <orpheus.numerics.quadrature.Quadrature.angular_frame>` — the
  :class:`~orpheus.numerics.basis.SphericalHarmonicBasis` of order
  :math:`L` bound to an :math:`S^2` cubature. It is a
  :class:`~orpheus.numerics.frame.GalerkinFrame` (``test is trial``)
  and a **4π-tight frame**. The headline future
  :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` consumer is
  **eigenvalue-consistent spatial homogenisation** (adjoint-weighted;
  §18 of Grand Report v3); spectrum-weighted energy condensation
  (§17) is the second.

- Every concrete :class:`~orpheus.numerics.frame.GalerkinFrame`
  satisfies the **idempotency-on-coefficients** invariant on a
  sufficiently-exact quadrature:

  .. math::

     M \, R \;=\; c_{V}\,I_{W},

  where :math:`c_V` is a scalar that depends on the inner-product
  convention of :math:`V`. For the SN spherical-harmonic frame on a
  Lebedev quadrature, :math:`c_V = 4\pi` — this is the **L1
  idempotency** identity :eq:`pi-r-equals-4pi-i` verified at multiple
  :math:`L` against multiple Lebedev orders. (A 4π-tight frame is one
  whose frame operator :math:`S = T^*T` is :math:`4\pi` times the
  identity; the tightness constant IS this :math:`c_V`.)


The discrete frame — analysis, synthesis, and the frame operator
================================================================

The :math:`(R, M)` pair is the language of **frame theory**
(Christensen 2016, *An Introduction to Frames and Riesz Bases*). A
discrete frame is a countable family :math:`\{e_k\}` in a Hilbert
space :math:`V` for which two operators are defined:

* the **analysis operator** :math:`T : V \to W`,
  :math:`(T f)_k = \langle e_k, f \rangle_V` — it *analyses* a
  function into its coefficients against the frame elements;
* the **synthesis operator** :math:`T^* : W \to V`,
  :math:`T^* c = \sum_k c_k\,e_k` — the formal adjoint of
  :math:`T`, the bare expansion with NO weighting and NO dual
  factor.

Their composition is the **frame operator**
:math:`S = T^* T : V \to V`. A frame is **tight with constant
:math:`c`** when :math:`S = c\,I` — the frame elements then behave
like an orthonormal basis up to the scalar :math:`c`, and the
inversion is trivial: :math:`f = c^{-1} T^* T f`.

In the ORPHEUS algebra, the analysis operator IS the **analysis face**
:math:`M = T` (measured against the test basis), and the
**reconstruction** :math:`R` is the **canonical-dual synthesis** —
:math:`T^*` weighted by the dual frame's Gram-inverse so that
:math:`M R` recovers the band-limited identity (up to tightness). The
bare :math:`T^* = S_0` (the *naked synthesis*) is the shared
:meth:`~orpheus.numerics.basis.Basis.synthesize` primitive on the
:class:`~orpheus.numerics.basis.Basis`; the analysis face :math:`M`
and the reconstruction face :math:`R` are each :math:`S_0`
post-multiplied by exactly one diagonal weight family (the quadrature
weight :math:`w_n` for analysis; the addition-theorem factor
:math:`2\ell+1` for reconstruction).

Given a fine space :math:`V` (e.g. :math:`L^2(S^2)` for the angular
flux) and a coarse coefficient space :math:`W` (e.g. polynomials of
degree :math:`\le L` on :math:`S^2`), a Galerkin-style discretisation
splits as:

.. math::
   :label: galerkin-pair

   R \;:\; W \to V, \qquad
   M \;:\; V \to W,
   \qquad
   M \, R \;=\; c_V \, I_W

where :math:`c_V` is the frame's tightness constant
(:math:`c_V = 1` in the fully-orthonormal case;
:math:`c_V = 4\pi` for the no-prefactor real spherical harmonics).

.. vv-status: galerkin-pair documented

The frame is fully determined by **three** ingredients — and the
third, the discipline, is the type, not a fourth parameter:

1. The **domain** :math:`V` and its inner product
   :math:`\langle \cdot, \cdot \rangle_V` — fixed by the
   :class:`~orpheus.numerics.measure.DiscreteMeasure` (the sample
   nodes and the fixed quadrature-weight :math:`L^2` metric). For SN
   angular flux, :math:`V = L^2(S^2)` and the inner product is the
   W-weighted discrete sum
   :math:`\langle f, g \rangle_W = \sum_n w_n f_n g_n` on the
   angular cubature.
2. The **trial basis** of :math:`W` — fixed by the
   :class:`~orpheus.numerics.basis.Basis` (the synthesis / trial
   side; it owns the reconstruction). For SN scattering, the basis is
   the real spherical harmonics :math:`\{Y_\ell^m\}_{\ell \le L}`.
3. The **test basis** — the analysis (measured) side, fixed by the
   frame **type**. A :class:`~orpheus.numerics.frame.GalerkinFrame`
   uses the trial basis itself (``test is trial``); a
   :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` carries an
   explicit, generally different test basis (e.g. the indicator
   weighted by the within-group spectrum, or by the region flux).

Once the basis and the measure are bound and the frame type is
chosen, :math:`M` and :math:`R` are uniquely determined up to the
:math:`c_V` normalisation.


The Galerkin discipline
=======================

The Galerkin discipline is characterised by **test space equals
trial space** — the :class:`~orpheus.numerics.frame.GalerkinFrame`
case. The defining identity is

.. math::
   :label: galerkin-self-adjoint

   M^* \;=\; R
   \quad \text{(under the } V \text{ inner product, orthonormal basis)}.

.. vv-status: galerkin-self-adjoint documented

i.e. the analysis face's Hilbert adjoint is its reconstruction. This is
why a single basis :math:`\{e_k\}` produces both :math:`M` and
:math:`R`:

.. math::
   :label: galerkin-construction

   (M f)_k &\;=\; \langle e_k, f \rangle_V, \\
   R \, c     &\;=\; \sum_k c_k\,e_k.

.. vv-status: galerkin-construction documented

.. warning::

   The identity :math:`M^* = R` holds when the basis
   :math:`\{e_k\}` is orthonormal in :math:`V`. When the basis is
   only orthogonal — the case for the no-:math:`4\pi/(2\ell+1)`-
   prefactor real spherical harmonics ORPHEUS uses — the strict
   Hilbert adjoint :math:`M^*` and the addition-theorem
   reconstruction face :math:`R` differ by a **diagonal-in-:math:`\ell`
   scaling**. Specifically the strict adjoint is the *naked*
   synthesis (no :math:`(2\ell+1)` factor), while the frame's
   ``reconstruction`` face
   (:attr:`FrameBase.reconstruction
   <orpheus.numerics.frame.FrameBase.reconstruction>`) carries the
   :math:`(2\ell+1)` factor that the Pℓ scattering reconstruction
   needs:

   .. math::

      (M^* c)_n
      &\;=\; \sum_{\ell, m} Y_\ell^m(\hat\Omega_n)\,c_\ell^m
        \quad\text{(no factor — strict adjoint)}, \\
      (R c)_n
      &\;=\; \sum_\ell (2\ell+1)\,\sum_m Y_\ell^m(\hat\Omega_n)\,
             c_\ell^m
        \quad\text{(with factor — addition-theorem)}.

   The analysis face's representation transpose
   :meth:`frame.analysis.apply_transpose
   <orpheus.numerics.basis.Basis.analyze_transpose>` is
   :math:`w_n\,S_0` (the *naked* synthesis weighted by the
   quadrature weight); its metric-aware Hilbert adjoint
   ``frame.analysis.H`` is :math:`M^* = g_C\,S_0`; and
   :meth:`frame.reconstruction.apply
   <orpheus.numerics.basis.Basis.reconstruct>` returns :math:`R`
   (with the :math:`(2\ell+1)` factor). The capability dishonesty
   that conflated the bare transpose with the Hilbert adjoint was
   caught by QA review and corrected as ERR-039 (see the project's
   L0 error catalog) — :math:`R` and :math:`M^*` are both useful
   and coexist as distinct frame faces; they differ by exactly the
   Gram diagonal :math:`g_C = \mathrm{diag}(4\pi/(2\ell+1))`, so the
   docstrings and this page name them explicitly.

The Galerkin invariant :eq:`galerkin-pair` is then a consequence of
the basis being orthogonal in :math:`V`-inner-product. Concretely,
for the spherical-harmonic Galerkin pair:

.. math::

   (M R)_{\ell m, \ell' m'}
   &\;=\; \sum_n w_n\,Y_\ell^m(\hat\Omega_n)\,
                       Y_{\ell'}^{m'}(\hat\Omega_n) \\
   &\;=\; \frac{4\pi}{2\ell+1}\,
          \delta_{\ell\ell'}\,\delta_{mm'},

so :math:`M R = \mathrm{diag}(4\pi/(2\ell+1))`. Composing with the
reconstruction face's :math:`(2\ell+1)` factor (the addition-theorem
weight) yields :math:`M R = 4\pi I` — the L1 identity that the
test
``tests/numerics/test_spherical_harmonic_space.py``
verifies at :math:`L = 2,\,3,\,4` (see :eq:`pi-r-equals-4pi-i` in
:ref:`spherical-harmonics`). This :math:`4\pi` is precisely the
frame's **tightness constant** :math:`c_V`: the frame operator
:math:`S = T^*T` equals :math:`4\pi\,I`, so the spherical-harmonic
frame is a 4π-tight frame.

.. note::

   The identity :math:`M R = 4\pi I` is **not** identity-on-the-
   nose because the no-prefactor convention pushes the
   :math:`4\pi/(2\ell+1)` factor onto the orthogonality. A
   spherical-harmonic frame with a strict :math:`M R = I`
   invariant could be built by dividing the analysis weights by
   :math:`4\pi`, but the project chose to absorb the factor at the
   reconstruction face (the :math:`(2\ell+1)` weight) so the
   addition theorem reads cleanly. See :ref:`spherical-harmonics`
   for the convention rationale.


The Petrov-Galerkin discipline
==============================

The Petrov-Galerkin discipline is characterised by **test space
differs from trial space** — the
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame` case. The
:math:`(R, M)` pair is built from two distinct bases —
:math:`\{e_k\}` for the trial space (the reconstruction basis, owned
by the :class:`~orpheus.numerics.basis.Basis`) and :math:`\{f_k\}`
for the test space (the explicit ``test_basis`` carried by the
frame):

.. math::
   :label: petrov-galerkin-construction

   (M g)_k \;=\; \langle f_k, g \rangle_V, \qquad
   R \, c \;=\; \sum_k c_k\,e_k.

.. vv-status: petrov-galerkin-construction documented

The pair satisfies :math:`M R = I_W` (the coefficient extraction
:math:`G^{-1} M` uses the *cross* Gram
:math:`G_{kj} = \langle f_k, e_j \rangle`), but :math:`M^* \ne R` —
the distinct test space makes the Hilbert adjoint distinct from the
reconstruction (an oblique, not canonical, dual).

The canonical Petrov-Galerkin pairs in reactor physics:

.. list-table:: Petrov-Galerkin pairs
   :header-rows: 1
   :widths: 22 24 28 26

   * - Use
     - Trial basis :math:`\{e_k\}`
     - Test basis :math:`\{f_k\}`
     - Reference
   * - Energy condensation
     - Indicator on broad group :math:`G`
     - Within-group spectrum
       :math:`\phi_g \cdot \mathbf{1}_{g \in G}`
     - Hébert 2009, §6.2
   * - Spatial homogenisation (reaction-rate)
     - Indicator on region :math:`R`
     - Region forward flux
       :math:`\phi \cdot \mathbf{1}_{i \in R}`
     - Smith 1986; Hébert 2009 §13
   * - Spatial homogenisation (eigenvalue-consistent)
     - Indicator on region :math:`R`
     - Region **adjoint**
       :math:`\varphi^* \cdot \mathbf{1}_{i \in R}`
     - Hébert 2009 §13; Stamm'ler & Abbate 1983
   * - Stochastic Galerkin
     - Polynomial-chaos basis (Hermite, Legendre)
     - Same basis under PCE inner product
     - Xiu & Karniadakis 2002

In each case the test basis encodes a **physical weighting** — the
importance of the fine-space slot from the solver's perspective — so
that the coarse coefficients faithfully preserve reaction rates /
flux-volume integrals / variance moments.

.. _petrov-galerkin-not-weighted-metric:

.. note:: **Why this is Petrov-Galerkin and not "Galerkin in a
   weighted metric".** It is tempting to absorb the test weight
   (:math:`\phi`, :math:`\varphi^*`) into the *measure* and call the
   result an orthogonal (Galerkin) projection in an
   :math:`L^2(\phi V)` metric. That fold is legitimate only for the
   forward-flux, reaction-rate-only row of the table: there the test
   weight and the integrand multiplier coincide, and the two readings
   are the same map. It **breaks** for the eigenvalue-consistent row,
   where the test weight is the **adjoint** :math:`\varphi^*` and the
   integrand is the forward flux — these are different functions, so
   no single metric on the measure reproduces the bilinear functional
   :math:`\langle \varphi^*, \Sigma\,\phi\rangle`. Folding either one
   into the metric would mis-weight the other. The discipline must
   therefore live on the **test side** (the frame type), and the
   measure stays a fixed :math:`L^2` metric. The full adjoint-weighted
   derivation lands with the concrete homogenisation frame (a later
   phase of ``refactor/operator-inverse-algebra``); this page fixes
   the *architecture* of where the weighting lives, not the numerical
   derivation of the collapse.

No concrete Petrov-Galerkin frame ships in the P1 carve; the
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame` base is in place
so the first concrete instance — eigenvalue-consistent spatial
homogenisation — lands as a ``test_basis`` choice, not a new
mechanism.


.. _frame-eigenbasis-ownership:

Why the discipline splits — an operator owns its frame iff the frame is its eigenbasis
======================================================================================

The page so far has *asserted* the discipline of each consumer
axis-by-axis: the angular spherical-harmonic frame is Galerkin
(:eq:`galerkin-pair`), energy condensation and spatial homogenisation
are Petrov-Galerkin (the consumer table). This section supplies the
**root cause** that decides which axis lands in which discipline, and
with it the deeper architectural fact behind the whole Frame campaign:

  **An operator owns its frame if and only if the frame is its
  eigenbasis — i.e. the basis that diagonalises the operator by a
  symmetry of the phase space.**

For the angular axis the symmetry is the rotation group :math:`SO(3)`,
the eigenbasis is the spherical harmonics, the diagonalisation is a
*theorem* (Funk–Hecke + Schur), and the owner is the **scattering
operator** — hence a :class:`~orpheus.numerics.frame.GalerkinFrame`
(orthogonal eigenbasis, ``test is trial``). For the energy and spatial
axes there is **no such symmetry**, no eigenbasis, and therefore no
operator that owns the frame: the projection is solution-*weighted*
and lives on the test side — hence a
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame`. The two
disciplines are *one* structural cause read at two axes, not two
independent conventions.

The structural leg — Funk–Hecke makes the spherical harmonics scattering's eigenbasis
-------------------------------------------------------------------------------------

The anisotropic scattering source operator is an angular **integral
kernel** (:ref:`integral-kernel-category`, :doc:`operator_algebra`):
the source at ordinate :math:`\hat\Omega` reads the flux at *every*
ordinate, weighted by the scattering kernel

.. math::
   :label: scattering-zonal-kernel

   (S_{\rm aniso}\,\psi)(\hat\Omega)
   \;=\; \int_{4\pi}
         \Sigma_s(\hat\Omega \cdot \hat\Omega')\,
         \psi(\hat\Omega')\;d\hat\Omega',

where the kernel depends on the incoming and outgoing directions
**only through their cosine** :math:`\hat\Omega \cdot \hat\Omega'`. A
kernel of this form is called a **zonal** kernel on the sphere
:math:`S^2` (it is invariant under a simultaneous rotation of both
arguments). Two classical theorems pin its spectrum.

.. vv-status: scattering-zonal-kernel documented
   The zonal-kernel form of the anisotropic scattering source is the
   literature-standard transport definition (Bell & Glasstone 1970
   §1.6; Lewis & Miller 1993 §4.7); it is a transcription, not a
   solver claim. The implementing kernel R∘Λ∘M is pinned by the
   0-ULP windowed-vs-full crosscheck
   ``tests/sn/operators/test_scattering_kernel_crosscheck.py`` and the
   addition-theorem identity :eq:`real-sh-addition-theorem`.

**Funk–Hecke theorem.** For any zonal kernel
:math:`k(\hat\Omega \cdot \hat\Omega')` on :math:`S^2`, the spherical
harmonics are eigenfunctions of the integral operator
:math:`(T_k f)(\hat\Omega) = \int_{S^2} k(\hat\Omega\cdot\hat\Omega')\,
f(\hat\Omega')\,d\hat\Omega'`, with an eigenvalue that depends on
:math:`\ell` **only** (not on :math:`m`):

.. math::
   :label: funk-hecke-eigenvalue

   T_k\,Y_\ell^m \;=\; \lambda_\ell\,Y_\ell^m,
   \qquad
   \lambda_\ell \;=\; 2\pi \int_{-1}^{+1} k(t)\,P_\ell(t)\,dt.

.. vv-status: funk-hecke-eigenvalue documented
   The Funk–Hecke eigenvalue formula is a classical result (Müller
   1966, *Spherical Harmonics*, Lecture Notes in Mathematics 17,
   §"Funk-Hecke"); transcribed here as the structural ground for the
   ownership ruling. The eigenvalues realised in code (the per-ℓ
   Legendre moments Σ_{s,ℓ}) are the diagonal of
   :class:`~orpheus.transport.operators.scattering.LegendreMomentScattering` Λ.

Applied to the scattering kernel
:math:`k = \Sigma_s(\hat\Omega\cdot\hat\Omega')`, the eigenvalues are
exactly the **Legendre moments of the differential scattering cross
section**, :math:`\lambda_\ell = \Sigma_{s,\ell}` — which are
precisely the per-:math:`\ell` block entries of the diagonal operator
:math:`\Lambda` =
:class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`
(:eq:`scattering-as-tensor-product-sum`, :doc:`operator_algebra`). The
spherical harmonics are therefore not *a* convenient basis for
scattering — they are *the* eigenbasis, forced by the rotational
invariance of the kernel.

**The kernel factorisation is the spectral theorem written out.** A
self-adjoint operator :math:`A` on a finite-dimensional space has the
spectral decomposition :math:`A = U\,\Sigma\,U^*`, with :math:`U` the
unitary whose columns are the eigenvectors and :math:`\Sigma` the
diagonal of eigenvalues. The discrete ORPHEUS scattering kernel is
*literally* this decomposition:

.. math::
   :label: scattering-spectral-theorem

   S_{\rm aniso}
   \;=\;
   \underbrace{R}_{=\,U}\;\circ\;
   \underbrace{\Lambda}_{=\,\Sigma}\;\circ\;
   \underbrace{M}_{=\,U^*},

with

* :math:`M` (the frame's **analysis** face, ``frame.analysis``) =
  :math:`U^*` — the change of basis *into* the eigenbasis (project the
  flux onto its harmonic moments :math:`\phi_\ell^m`);
* :math:`\Lambda` (=
  :class:`~orpheus.transport.operators.scattering.LegendreMomentScattering`) =
  :math:`\Sigma` — the diagonal multiply by the spectrum
  :math:`\Sigma_{s,\ell}`, one scalar per :math:`\ell`-block;
* :math:`R` (the frame's **reconstruction** face,
  ``frame.reconstruction``) = :math:`U` — the synthesis *out of* the
  eigenbasis (rebuild the per-ordinate source).

The **addition theorem** :eq:`real-sh-addition-theorem` —
:math:`\sum_m Y_\ell^m(\hat\Omega)\,Y_\ell^m(\hat\Omega') =
P_\ell(\hat\Omega\cdot\hat\Omega')` — is exactly the *spectral
resolution* of the zonal kernel: it expresses the rank-:math:`(2\ell+1)`
projector onto the degree-:math:`\ell` eigenspace as an outer product
of harmonics. Reading :math:`S = R\circ\Lambda\circ M` as
:math:`U\Sigma U^*` is what makes the conjugation
:math:`S = \texttt{frame.conjugate}(\Lambda)` (the scattering **2-cell**
of the operator-algebra double category, :doc:`operator_algebra`) a
*spectral* statement and not merely a convenient bracketing.

**Schur's lemma fixes the block structure and the weights.** The
function space :math:`L^2(S^2)` decomposes into the
:math:`SO(3)`-irreducible subspaces
:math:`L^2(S^2) = \bigoplus_\ell V_\ell`, where
:math:`V_\ell = \mathrm{span}\{Y_\ell^m\}_{m=-\ell}^{\ell}` is the
degree-:math:`\ell` eigenspace of dimension :math:`2\ell+1`. The
scattering source operator commutes with every rotation (it is built
from a zonal kernel), so it lies in the :math:`SO(3)`-**commutant**.
By **Schur's lemma**, any operator in the commutant acts as a *scalar*
on each irreducible block — which is the :math:`m`-independence of the
Funk–Hecke eigenvalue, now derived from symmetry rather than computed
from an integral. The block dimensions :math:`\dim V_\ell = 2\ell+1`
are the origin of:

* the :math:`(2\ell+1)` reconstruction factor on the frame's
  reconstruction face (:eq:`sh-addition-theorem-reconstruction`,
  :ref:`spherical-harmonics`) — the irrep dimension; and
* the Gram diagonal :math:`g_C = \mathrm{diag}(4\pi/(2\ell+1))`
  (:eq:`sh-space-metric`) — the :math:`SO(3)` Plancherel weight on
  each block.

So the entire numerical apparatus of the spherical-harmonic frame —
the per-:math:`\ell` block structure, the :math:`(2\ell+1)` factor,
the :math:`4\pi`-tightness — is the representation theory of
:math:`SO(3)` acting on a rotationally-invariant kernel. The frame is
Galerkin **because** the eigenbasis of a self-adjoint (here
:math:`SO(3)`-invariant) operator is orthogonal: ``test is trial`` is
forced by symmetry, not chosen.

The asymmetry that fixes ownership — streaming does not diagonalise
-------------------------------------------------------------------

If both transport operators were diagonalised by the spherical
harmonics, "scattering owns the frame" would be a coin toss. They are
not, and the asymmetry is what assigns ownership.

The **streaming** operator :math:`\hat\Omega\cdot\nabla` carries the
direction :math:`\hat\Omega`, which is itself the :math:`\ell = 1`
vector irrep of :math:`SO(3)`. It does **not** commute with rotations
(rotating the frame rotates the gradient direction), so it is **not**
in the commutant and **not** diagonalised by the harmonics. By the
Clebsch–Gordan rule
:math:`V_1 \otimes V_\ell = V_{\ell-1}\oplus V_\ell\oplus V_{\ell+1}`,
multiplication by a component of :math:`\hat\Omega` couples
:math:`V_\ell` to :math:`V_{\ell\pm 1}`:

.. math::
   :label: streaming-pn-recurrence

   \mu\,Y_\ell^m
   \;=\;
   \frac{\ell+1}{2\ell+1}\,Y_{\ell+1}^{m}
   \;+\;
   \frac{\ell}{2\ell+1}\,Y_{\ell-1}^{m}

— the **Pℓ moment recurrence**, the block-**tridiagonal** coupling that
makes the PN coefficient matrix tridiagonal in :math:`\ell` rather
than diagonal (Fletcher 1983, Eq. 5; Hébert 2009, §3.6–3.7). Streaming
in the harmonic basis is *tolerated*, not *diagonalised*: the basis is
chosen to make collision diagonal, and streaming is then expressed
(awkwardly, tridiagonally) in those same coordinates.

.. vv-status: streaming-pn-recurrence documented
   The μ·Y_ℓ recurrence is the standard Legendre/spherical-harmonic
   recurrence (Fletcher 1983 NSE 84:33 Eq. 5; Hébert 2009 §3.6); a
   transcribed structural identity, not a solver claim. ORPHEUS does
   not yet ship a PN solver — the recurrence is documented here as the
   asymmetry that fixes frame ownership, not as a verified code path.

The ownership conclusion is then forced:

.. list-table:: Which operator the spherical-harmonic basis diagonalises
   :header-rows: 1
   :widths: 26 26 28 20

   * - Operator
     - Symmetry of the kernel
     - Action in the SH basis
     - Diagonalised?
   * - Scattering :math:`\Sigma_s(\hat\Omega\cdot\hat\Omega')`
     - :math:`SO(3)`-invariant (zonal)
     - **diagonal** per :math:`\ell`-block (Funk–Hecke + Schur)
     - **yes** — its eigenbasis
   * - Streaming :math:`\hat\Omega\cdot\nabla`
     - :math:`\ell=1` tensor; **not** invariant
     - block-**tridiagonal** :math:`\ell\!\leftrightarrow\!\ell\pm 1`
       (Clebsch–Gordan)
     - no — merely tolerated

Because the spherical harmonics are the eigenbasis of *scattering* and
nothing else in the transport operator, the frame is owned by the
scattering operator. In ORPHEUS this ownership is a concrete fact:
:class:`~orpheus.transport.operators.scattering.ScatteringOperator` holds the frame as
a cached property and binds its order to the scattering order,
``quadrature.angular_frame(self.scattering_order)`` (the canonical
constructor + the :math:`L`-binding). The frame *object* lives in the
method-agnostic :class:`~orpheus.numerics.frame.GalerkinFrame`
hierarchy; the *constructor ownership* sits on scattering.

Literature corroboration — no falsifier across six transport families
---------------------------------------------------------------------

The structural argument is corroborated by every documented transport
method: in SN, MoC, CP, PN, first-collision-source/ray-effect, and
random-ray, the flux→spherical-harmonic-moment projection
:math:`M = \int Y_\ell^m\,\psi\,d\hat\Omega` is invoked **solely** to
evaluate the anisotropic scattering source. A falsifier would be a
documented, structurally-independent *non-scattering* use of the flux
moment projection; none exists.

.. list-table:: Literature: the flux→SH-moment projection is anisotropic-scattering-rooted
   :header-rows: 1
   :widths: 22 30 48

   * - Reference
     - Equation / section
     - What it establishes
   * - Hébert 2009, *Applied Reactor Physics*
     - §3.3, Eq. (3.55) [the projection :math:`M`]; Eq. (3.54)
       [its sole use]; Eq. (3.42); Eq. (3.57)
     - :math:`M` (Eq. 3.55) is used **only** in the scattering source
       (Eq. 3.54). The integral / characteristic form natively wants
       **isotropic** sources (Eq. 3.42); fission is isotropic
       (Eq. 3.57). Only anisotropic scattering forces the moments.
   * - Brockmann 1981, NSE **77** (4), 377
     - Eq. (47) [the Legendre flux moment]
     - Introduces :math:`\Phi_\ell(x,E) = 2\pi\int P_\ell(\mu)\,\Phi\,d\mu`
       expressly "considering the problem of anisotropic neutron
       scattering"; the *same* projection is reused across SN, FEM,
       and orders-of-scattering.
   * - Fletcher 1983, NSE **84**, 33
     - Eq. (7) [moment equation]; Eq. (5) [streaming recurrence]
     - The moment equation is **diagonal** in :math:`\ell` "because of
       the orthogonality of spherical harmonics" (scattering); the
       streaming term produces the :math:`\ell\!\leftrightarrow\!\ell\pm1`
       **tridiagonal** coupling. PN's moments exist because the basis
       diagonalises scattering.
   * - Ahrens 2014, *Lagrange Discrete Ordinates*, arXiv:1405.3968
     - Eq. (7); abstract
     - The **negative-space proof**: LDO *removes* :math:`M` —
       "no spherical harmonic moments are needed" — precisely by
       reformulating the scattering source. An authority stating the
       only reason standard SN computes :math:`M` is the scattering
       source.
   * - External / boundary sources (Hébert 2009, Eq. 3.30/3.168)
     - —
     - An anisotropic external source is **specified** in spherical
       harmonics as input data :math:`Q_\ell^m`; it is never
       *projected* from the flux. Not a flux→moment projection.
   * - Anisotropic BCs (albedo / white / specular)
     - Hébert 2009, Eqs. 3.43–3.47
     - Handled in **ordinate (direction) space**, not via moments.
       Only PN expresses BCs in moments — because moments are PN's
       native unknowns, which exist to diagonalise scattering.

The single recurring exception, PN, is the method that makes the
moments the unknowns *in order to* diagonalise scattering — so even
there the root cause is scattering. The convergence of five
independent references on the same conclusion, plus Ahrens' explicit
removal of :math:`M`, leaves the claim "the spherical-harmonic angular
projection is intrinsically a scattering concern" with **zero
cross-validation against any non-scattering flux-moment use**.

The unifying principle — symmetry decides Galerkin vs Petrov-Galerkin
---------------------------------------------------------------------

The eigenbasis criterion is the *single structural cause* of the
discipline split this page documents. Reading it across all three
reduction axes:

.. list-table:: Symmetry decides the discipline
   :header-rows: 1
   :widths: 22 24 28 26

   * - Reduction axis
     - Phase-space symmetry
     - Eigenbasis?
     - Discipline (frame type)
   * - **Angular scattering**
     - :math:`SO(3)` rotational (zonal kernel)
     - **yes** — spherical harmonics, *orthogonal* (Funk–Hecke + Schur)
     - **Galerkin**, owned by the scattering operator
   * - **Energy condensation**
     - none (general :math:`G\times G` group-transfer matrix)
     - no — no Funk–Hecke, no clean spectrum
     - **Petrov-Galerkin**, solution-weighted, owned by no operator
   * - **Spatial homogenisation**
     - none (arbitrary heterogeneous geometry)
     - no
     - **Petrov-Galerkin**, solution-weighted, owned by no operator

When a symmetry of the phase space diagonalises the operator, the
eigenbasis is *forced* — and because the eigenbasis of a self-adjoint
operator is orthogonal, the projection is **Galerkin** (``test is
trial``, :math:`M^* = R`) and the operator **owns** the frame (the
frame's order is the operator's order). When there is no symmetry —
the energy group-transfer matrix and the spatial homogenisation kernel
have no rotational or other structure to exploit — there is no
eigenbasis, no operator that naturally owns the basis, and the
coarse-graining must instead **weight** the projection by the solution
(the within-group spectrum :math:`\phi_g`, the region flux
:math:`\phi`, or the adjoint :math:`\varphi^*`). That weighting is a
*test* basis distinct from the trial basis, so the projection is
**Petrov-Galerkin** (:math:`M^* \ne R`) and lives on the test side,
owned by the *projection verb*, never by an operator. This is also why
**fission does not own an angular frame**: fission's concern is the
energy axis (the :math:`\chi\nu\Sigma_f` group-transfer), which has no
eigenbasis — its angular emission is isotropic
(:math:`\ell = 0`-only), so there is no harmonic structure for it to
own.

The architectural payoff is that the Galerkin-vs-Petrov-Galerkin
type split (:class:`~orpheus.numerics.frame.GalerkinFrame` vs
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame`) is *derived*
from a single physical question — *does a symmetry diagonalise the
operator?* — rather than asserted axis-by-axis. The
:ref:`homogenisation note above <petrov-galerkin-not-weighted-metric>`
(why folding the solution into the metric breaks for the
eigenvalue-consistent case) is the *converse* of the same principle:
absent an eigenbasis, the solution-weighting cannot be hidden in a
fixed :math:`L^2` metric — it is irreducibly a distinct test space.

.. _frame-eigenbasis-relocation-tripwire:

The relocation tripwire — when scattering stops owning the constructor
----------------------------------------------------------------------

"Scattering owns the frame" is true **today** because scattering is
the *only* consumer of the spherical-harmonic frame whose truncation
order :math:`L` is physically meaningful. The constructor ownership
:meth:`ScatteringOperator.frame
<orpheus.transport.operators.scattering.ScatteringOperator.frame>` binds the frame
order to ``self.scattering_order``. This binding **relocates** to the
discipline-neutral factory
:meth:`Quadrature.angular_frame(L)
<orpheus.numerics.quadrature.Quadrature.angular_frame>` the moment a
**second** consumer arrives with an :math:`L` *independent of*
``scattering_order``:

* an **output detector / response functional** of anisotropic order
  :math:`L_d > L_{\rm scatter}` — a flux moment projection whose
  truncation is set by the *detector*, not the scattering kernel
  (structurally independent); or
* a **PN / SPN flux** expansion of order
  :math:`L_{\rm flux} > L_{\rm scatter}` — needing
  :math:`\max(L_{\rm flux}, L_{\rm scatter})`, not ``scattering_order``.

No such consumer exists today: the only output functional ORPHEUS
computes is the :math:`\ell = 0` scalar flux (via the
:class:`~orpheus.sn.solver.SNSolver`'s angular integration), which does
**not** use the frame. The factory
:meth:`~orpheus.numerics.quadrature.Quadrature.angular_frame` already
exists and already anticipates this relocation — its naming tripwire
names the permanent *angular axis*, not the contingent
spherical-harmonic basis, so a second consumer is a signature change
(``angular_frame(basis=...)``), not a rename. Until that second
consumer lands, the canonical constructor home is the scattering
operator, because scattering is the operator whose eigenbasis the
frame *is*.

A second, structurally distinct trigger is **cross-method use of**
:class:`~orpheus.transport.operators.scattering.ScatteringOperator` (`#261`). The
operator is method-agnostic in principle — every transport method
scatters — but a
:class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame` needs an
angular **measure** (a :class:`~orpheus.numerics.quadrature.Quadrature`)
to exist at all, and CP / MoC carry none (CP integrates angle away into
collision probabilities; MoC uses a track quadrature, not an
:math:`S^2` ordinate set). So the moment scattering is consumed by a
measure-free method, the frame cannot live *as a field on* the shared
operator. Two resolutions are open (deferred to `#261`; user,
2026-06-25): **(a) relocate** the frame to where the angular measure
lives (the
:meth:`~orpheus.numerics.quadrature.Quadrature.angular_frame` factory —
the original W-E idea), or **(b) specialize**
:class:`~orpheus.transport.operators.scattering.ScatteringOperator` per method — an SN
subclass that holds the frame (it carries the :math:`S^2` measure) over
a measure-free cross-method base. Where the independent-:math:`L`
triggers above make the *order* foreign, this one makes the *measure*
absent. The eigenbasis ruling still holds — scattering owns the frame
*wherever an angular measure exists to build it* — it merely stops being
expressible as a field on a single method-agnostic operator.


Cross-method consumer table
===========================

The frame pair is **infrastructure**, not SN-only. Every method that
lifts an angular / energy / spatial axis between fine and coarse
representations builds on these primitives, always through a frame of
the appropriate discipline type.

.. list-table:: Where the frame pair is consumed
   :header-rows: 1
   :widths: 22 22 22 16 18

   * - Consumer
     - Frame type
     - Pair
     - Status
     - Reference
   * - **SN aniso scattering**
     - GalerkinFrame
     - ``frame.analysis`` /
       ``frame.reconstruction``
       (``Quadrature.angular_frame(L)``)
     - **Live** (Wave 1)
     - §9 (Grand Report v3 line 1230)
   * - PN solver
     - GalerkinFrame
     - Same SH ``frame.analysis``
       on the moment-space side
     - Pending (PN solver not implemented)
     - §10 (lines 1299–1305)
   * - Energy condensation
     - PetrovGalerkinFrame
     - Within-group-spectrum test basis, indicator trial basis
     - Pending (concrete ``test_basis`` needed)
     - §17 (line 3935); Hébert 2009 §6.2
   * - Spatial homogenisation (eigenvalue-consistent)
     - PetrovGalerkinFrame
     - Adjoint-flux test basis, indicator trial basis
     - Pending — the **headline** PG consumer
     - §18; Hébert 2009 §13
   * - Stochastic Galerkin
     - GalerkinFrame
     - Polynomial-chaos basis on the random-input axis
     - Pending
     - §15A.7
   * - MC adjoint moments
     - GalerkinFrame
     - Same SH ``frame.analysis``
       used as a variance-reduction estimator
     - Pending
     - Lewis & Miller 1993 §10
   * - SN sensitivity
     - GalerkinFrame (adjoint)
     - Same pair, applied to the adjoint flux
     - Pending
     - Cacuci 2003

The two architectural payoffs:

* **One mechanism per discipline type**, not one per consumer. The
  spherical-harmonic :class:`~orpheus.numerics.frame.GalerkinFrame`
  emits the same ``analysis`` / ``reconstruction`` faces whether SN
  uses them for scattering or PN uses them for streaming — the
  difference is which axis the face is wrapped onto via the
  :class:`~orpheus.numerics.operator.TensorProductOperator`
  algebra (see :ref:`operator-algebra` and the tensor-product
  section there).
* **One V&V chain per discipline**. The Galerkin idempotency tests
  in :file:`tests/numerics/test_spherical_harmonic_space.py` cover
  every :class:`~orpheus.numerics.frame.GalerkinFrame` consumer, not
  just SN. The Petrov-Galerkin frames will inherit the
  cross-Gram / rate-preservation tests when their concrete
  ``test_basis`` lands.


Discipline as a type, not a property or an operator marker
==========================================================

The discipline — Galerkin vs Petrov-Galerkin — is a genuine **kind of
object**, so it is carried by the frame **type** (Issue #268). Two
rejected alternatives clarify why.

**Rejected (a): discipline as a marker ABC on the operator role.** An
earlier draft put the discipline on the analysis operator via marker
mixins ``GalerkinProjection`` / ``PetrovGalerkinProjection``:

.. code-block:: python

   # RETIRED: discipline declared on the operator role
   class GalerkinProjection(AnalysisOperator, ABC): ...
   class PetrovGalerkinProjection(AnalysisOperator, ABC): ...

This declares a discipline on the *role* (:math:`M`) when the
discipline is really a fact about the *frame* — the relationship
between the test and trial spaces, which an analysis operator in
isolation cannot express. The marker ABCs were retired; the
:mod:`orpheus.numerics.projection` module now carries only the two
discipline-free operator roles
(:class:`~orpheus.numerics.projection.AnalysisOperator`,
:class:`~orpheus.numerics.projection.ReconstructionOperator`).

**Rejected (b): discipline as a derived property of the frame.** An
intermediate draft proposed collapsing the distinction to a boolean
property (``measure == basis.canonical_measure``), on the theory that
homogenisation is "really Galerkin in a weighted metric". That reading
folds the solution into the metric and breaks for the
eigenvalue-consistent (adjoint-weighted) case (see the Petrov-Galerkin
note above). A property cannot encode a genuinely different object;
the discipline is a TYPE.

**Shipped: discipline as a Liskov-correct type hierarchy.** The frame
type names the discipline, and the user's pedantic naming rule is
satisfied — **a reader of a type name knows its properties without
reading the docstring**:

.. code-block:: python

   from orpheus.numerics.frame import (
       FrameBase, PetrovGalerkinFrame, GalerkinFrame,
   )

   # Galerkin: test IS trial, Π* = R (a canonical dual)
   sh_frame = quad.angular_frame(L)            # -> GalerkinFrame

   # Petrov-Galerkin: explicit test basis, M* ≠ R (an oblique dual)
   homog = PetrovGalerkinFrame(
       basis=mesh.indicator_basis(),
       measure=mesh.volume_measure,
       test_basis=adjoint_weighted_indicators,  # the discipline
   )

A reader of ``PetrovGalerkinFrame(...)`` immediately knows the
analysis measures against a distinct test basis, so :math:`M^* \ne R`;
a reader of :class:`~orpheus.numerics.frame.GalerkinFrame` knows the
strengthened :math:`M^* = R` promise holds. The hierarchy answers the
discipline question without reading prose, and a
:class:`~orpheus.numerics.frame.GalerkinFrame` that is handed a
distinct ``test_basis`` raises (the contradiction is unrepresentable).


Numerical evidence
==================

The L1-tagged tests in
:file:`tests/numerics/test_spherical_harmonic_space.py` verify the
Galerkin discipline's invariants on the spherical-harmonic
:class:`~orpheus.numerics.frame.GalerkinFrame`'s ``analysis`` /
``reconstruction`` faces:

1. **Idempotency** (4π-tightness):
   :math:`M R c = 4\pi c` on band-limited
   coefficient input, verified at :math:`L = 2,\,3,\,4` against
   Lebedev orders :math:`7,\,13,\,17`. See
   :eq:`pi-r-equals-4pi-i` in :ref:`spherical-harmonics`.
2. **Adjoint pairing**:
   :math:`\langle M \psi, c \rangle = \langle \psi, M^* c \rangle_W`
   with :math:`M^* = g_C\,S_0`, verified to ``rtol=1e-12`` on a
   Lebedev order-13 grid at :math:`L = 3`.

The tests verify mathematical identities of the operator algebra
(V&V level L1 — equation verification by analytical reference). The
companion L0/foundation shape and capability tests verify software
invariants (frame face spaces, capability sets) and are tagged
accordingly. No Petrov-Galerkin numerical evidence ships yet — the
cross-Gram / rate-preservation gates land with the first concrete
``test_basis``.


Implementation map
==================

* :class:`~orpheus.numerics.frame.FrameBase` — the abstract discrete
  frame: binds a :class:`~orpheus.numerics.basis.Basis` to a
  :class:`~orpheus.numerics.measure.DiscreteMeasure` and emits the
  ``analysis`` (:math:`M = T`) and ``reconstruction`` (:math:`R`)
  faces. Carries the discipline-free mechanics (table, the two
  spaces, the reconstruction face, the analysis-face wiring); the
  single mechanism for every choice-dependent change-of-basis.
* :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` — the general
  discipline: an explicit ``test_basis`` distinct from the trial
  basis, so :math:`M^* \ne R`. The base for homogenisation and
  condensation.
* :class:`~orpheus.numerics.frame.GalerkinFrame` — the Galerkin
  specialisation (``test is trial``, :math:`M^* = R`). The angular
  spherical-harmonic frame is the canonical instance.
* :class:`~orpheus.numerics.basis.Basis` — the synthesis (trial)
  side ABC: tabulate, naked synthesis :math:`S_0`, the three
  weighted contractions, and the discrete Gram.
* :class:`~orpheus.numerics.basis.SphericalHarmonicBasis` — the
  first concrete basis (real spherical harmonics); carries the
  no-prefactor convention and the
  :attr:`~orpheus.numerics.basis.SphericalHarmonicBasis.addition_theorem_factor`
  :math:`(2\ell+1)`.
* :class:`~orpheus.numerics.projection.AnalysisOperator` — the
  abstract fine→coarse operator role :math:`M : V \to W`; the
  ``analysis`` face subclasses it.
* :class:`~orpheus.numerics.projection.ReconstructionOperator` —
  the abstract coarse→fine operator role :math:`R : W \to V`; the
  ``reconstruction`` face subclasses it.
* :meth:`Quadrature.angular_frame(L)
  <orpheus.numerics.quadrature.Quadrature.angular_frame>` — builds
  the order-:math:`L` spherical-harmonic
  :class:`~orpheus.numerics.frame.GalerkinFrame` on a quadrature; the
  single home of the :math:`S^2` embedding.

The full-space projector — the operator that projects the SN
:math:`(N, n_x, n_y, n_g)` angular flux onto the
:math:`(L+1, 2L+1, n_x, n_y, n_g)` moment field — is built as a
**tensor product** of the angular-axis analysis face :math:`M`
and identity operators on the spatial / energy axes:

.. code-block:: python

   from orpheus.numerics.operator import IdentityOperator

   frame = quad.angular_frame(L)
   M = frame.analysis
   M_full = M & IdentityOperator() & IdentityOperator() & IdentityOperator()

The ``&`` dunder constructs the
:class:`~orpheus.numerics.operator.TensorProductOperator`. See
:ref:`operator-algebra` and the **Tensor product algebra** section
there for the relationship between this operator-algebra type and
the underlying numpy primitives (``np.einsum``, ``np.tensordot``,
``np.kron``).


History — from operator classes to the discipline-type frame
============================================================

The spherical-harmonic projection and reconstruction were first
shipped (Wave 0 of the SN performance plan) as standalone operator
classes ``HarmonicMomentProjection`` / ``HarmonicMomentReconstruction``
under a three-level inheritance
(``ProjectionOperator`` → ``GalerkinProjection`` → concrete). Two
naming-audit corrections then established the discipline-must-be-
typed pedagogy.

The Frame/Basis carve (``refactor/operator-inverse-algebra``)
took the next step: the projection :math:`M = Y^*W` and the
addition-theorem reconstruction :math:`R = (2\ell+1)\,S_0` are NOT
two unrelated operator classes — they are the **two faces of one
discrete frame** binding the SH basis to the angular measure. The
standalone operator classes were retired into the frame faces:

* ``HarmonicMomentProjection`` → ``frame.analysis``
  (:attr:`FrameBase.analysis <orpheus.numerics.frame.FrameBase.analysis>`,
  the analysis face :math:`M = T`);
* ``HarmonicMomentReconstruction`` → ``frame.reconstruction``
  (:attr:`FrameBase.reconstruction
  <orpheus.numerics.frame.FrameBase.reconstruction>`, the
  reconstruction face :math:`R`);
* the :math:`(2\ell+1)` factor moved onto
  :attr:`SphericalHarmonicBasis.addition_theorem_factor
  <orpheus.numerics.basis.SphericalHarmonicBasis.addition_theorem_factor>`
  (one home for the SH convention).

The P1 discipline-type carve (Issue #268) took the final step: the
discipline, which an earlier draft had carried as marker ABCs on the
operator role and a later draft proposed collapsing to a frame
*property*, became the frame **type**
(:class:`~orpheus.numerics.frame.FrameBase` →
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame` →
:class:`~orpheus.numerics.frame.GalerkinFrame`). The architectural
payoff: one mechanism (the frame), the discipline visible in the type,
and the eigenvalue-consistent homogenisation case correctly typed as
Petrov-Galerkin (test = adjoint-weighted indicator) rather than
mis-folded into a weighted measure.


References
==========

* Brenner, S. C. and Scott, L. R. (2008). *The Mathematical Theory
  of Finite Element Methods*, 3rd ed. Springer. §3.4 (Galerkin /
  Petrov-Galerkin general framework — test vs trial space).
* Christensen, O. (2016). *An Introduction to Frames and Riesz
  Bases*, 2nd ed. Birkhäuser. (The analysis operator :math:`T`, the
  synthesis operator :math:`T^*`, the frame operator
  :math:`S = T^*T`, tight frames, and the canonical dual — the
  harmonic-analysis foundation of the
  :class:`~orpheus.numerics.frame.FrameBase` abstraction.)
* Bell, G. I. and Glasstone, S. (1970). *Nuclear Reactor Theory*.
  Van Nostrand Reinhold. §1.6 (spherical-harmonic moment
  projection in transport).
* Lewis, E. E. and Miller, W. F. Jr. (1993). *Computational Methods
  of Neutron Transport*. ANS. §4.7 (Pℓ Galerkin reconstruction with
  the :math:`(2\ell+1)` factor).
* Müller, C. (1966). *Spherical Harmonics*. Lecture Notes in
  Mathematics **17**, Springer. (The Funk–Hecke theorem: spherical
  harmonics are the eigenfunctions of any zonal kernel on
  :math:`S^2`, with eigenvalue
  :math:`\lambda_\ell = 2\pi\int_{-1}^{1} k(t) P_\ell(t)\,dt` — the
  structural ground for "the SH frame is scattering's eigenbasis".)
* Hébert, A. (2009). *Applied Reactor Physics*. Polytechnique. §3.3
  (the flux→SH-moment projection :math:`M`, Eq. 3.55, used **only**
  in the scattering source Eq. 3.54; fission isotropic Eq. 3.57;
  integral form natively isotropic Eq. 3.42), §3.6–3.7 (the streaming
  :math:`\ell\!\leftrightarrow\!\ell\pm1` recurrence), §6.2 (energy
  condensation as a Petrov-Galerkin projection), §13
  (eigenvalue-consistent / adjoint-weighted spatial homogenisation).
* Brockmann, H. (1981). *Treatment of anisotropic scattering in
  numerical neutron transport theory*. Nucl. Sci. Eng. **77** (4),
  377–414. Eq. (47) — the Legendre flux moment
  :math:`\Phi_\ell = 2\pi\int P_\ell(\mu)\Phi\,d\mu` is introduced
  expressly for the anisotropic-scattering source and reused across
  SN / FEM / orders-of-scattering.
* Fletcher, J. K. (1983). *The solution of the multigroup neutron
  transport equation using spherical harmonics*. Nucl. Sci. Eng.
  **84**, 33–46. Eq. (7) — the moment equation is diagonal in
  :math:`\ell` "because of the orthogonality of spherical harmonics"
  (scattering); Eq. (5) — the streaming term produces the
  block-tridiagonal :math:`\ell\!\leftrightarrow\!\ell\pm1` coupling.
* Ahrens, C. D. (2014). *Lagrange Discrete Ordinates: a new angular
  discretization for the three-dimensional transport equation*.
  arXiv:1405.3968. Eq. (7) and abstract — the negative-space proof:
  LDO **removes** the SH moment projection ("no spherical harmonic
  moments are needed") precisely by reformulating the scattering
  source.
* Cacuci, D. G. (2003). *Sensitivity and Uncertainty Analysis,
  Volume I*. CRC Press. (Adjoint flux moments and the Galerkin
  pair on the adjoint side.)
* Xiu, D. and Karniadakis, G. E. (2002). *The Wiener-Askey
  polynomial chaos for stochastic differential equations*. SIAM J.
  Sci. Comput. 24 (2), 619–644. (Stochastic Galerkin on the random
  input axis.)
* Grand Report v3 §5.7 (line 664), §17 (line 3935), §32.7 — the
  catalog entries that drove the placement of these primitives in
  :mod:`orpheus.numerics.projection`.
