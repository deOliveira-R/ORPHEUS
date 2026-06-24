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

   The stale "Galerkin in :math:`L^2(\phi V)`" derivation still lives
   in :ref:`sn-homogenization-galerkin-frame`
   (:doc:`discrete_ordinates`); that section is flagged for the same
   reversal under Issue #268 and should be read through this
   correction until it is rewritten.


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
* Hébert, A. (2009). *Applied Reactor Physics*. Polytechnique. §6.2
  (energy condensation as a Petrov-Galerkin projection), §13
  (eigenvalue-consistent / adjoint-weighted spatial homogenisation).
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
