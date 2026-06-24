.. _galerkin-projection:

==========================================================
Galerkin Projection — the discrete frame and its consumers
==========================================================

Every method in ORPHEUS that transitions a function between a
**fine** representation and a **coarse** one — discrete-ordinate
angular flux to spherical-harmonic moments, fine-energy fluxes to
broad-group cross-sections, regional flux to homogenised cross-
sections — does so via a **(reconstruction, analysis)** pair of
linear operators :math:`(R, \Pi)`. In harmonic-analysis language
this pair is exactly the two operational **faces of a discrete
frame**: the analysis operator :math:`\Pi = T` (sampled values →
coefficients) and the reconstruction :math:`R` (coefficients →
values, the canonical-dual synthesis). This page is the canonical
home for that pair, for the discrete-**frame** abstraction that
realises it (:class:`~orpheus.numerics.frame.Frame`), the
**discipline** (Galerkin vs Petrov-Galerkin — canonical-dual vs
oblique-dual) that distinguishes its variants, and the
**cross-method consumer catalog** of where the same primitives are
used in the codebase.

The user's architectural rule that motivated this page:

  *"Galerkin projection is generally useful — same machinery is
  used in cross-section energy condensation. Catalog everything
  that needs to be implemented and where."*

Binding a :class:`~orpheus.numerics.basis.Basis` to a
:class:`~orpheus.numerics.measure.DiscreteMeasure` through a single
:class:`~orpheus.numerics.frame.Frame` (Frame/Basis carve,
``refactor/operator-inverse-algebra``) puts one mechanism in front
of every consumer, instead of a separate projection /
reconstruction operator class per method.

.. contents::
   :local:
   :depth: 2


.. warning:: **Superseded framing (2026-06-23, Issue #268).** This page
   names spatial homogenisation and energy condensation as the future
   **Petrov-Galerkin** instances (the consumer table, the
   "Petrov-Galerkin pairs" table, and
   ``class EnergyCondensation(PetrovGalerkinProjection)``). The spatial-
   homogenisation carve (``refactor/operator-inverse-algebra``) **derived
   that both are Galerkin in their natural weighted metric**, not
   Petrov-Galerkin: homogenisation is the :math:`L^2(\phi V)`-orthogonal
   projection and condensation is the :math:`L^2(\text{spectrum})`-
   orthogonal projection. The "oblique-dual / Petrov-Galerkin" reading is
   an artifact of insisting on the *unweighted* :math:`\mu_V` metric — in
   the native flux·volume (resp. spectrum) metric the test and trial
   spaces are the **same** indicator span and the projection is
   orthogonal. The two readings are the same map (the flux multiplier can
   live in the test function *or* in the measure); the discrete
   :class:`~orpheus.numerics.frame.Frame` reads the weighting off the
   **measure**, so it sees the Galerkin case.

   Consequence:
   :class:`~orpheus.numerics.projection.PetrovGalerkinProjection` has
   **no concrete instance**, and the Galerkin/Petrov-Galerkin *type*
   distinction collapses to a derived **property** of the frame
   (``measure == basis.canonical_measure``), not a class hierarchy. See
   the worked derivation in :ref:`sn-homogenization-galerkin-frame`
   (the first concrete instance proving spatial homogenisation is
   :math:`L^2(\phi V)`-Galerkin). The full reframe of this page plus the
   retirement of the ``GalerkinProjection`` / ``PetrovGalerkinProjection``
   ABCs is tracked under **Issue #268** — this note flags the page's own
   staleness; the prose below is **not yet rewritten** and should be read
   through this correction.


Key Facts
=========

- A **frame** binds a :class:`~orpheus.numerics.basis.Basis` (the
  synthesis / trial side — the functions :math:`\{e_k\}` and their
  convention) to a :class:`~orpheus.numerics.measure.DiscreteMeasure`
  (the analysis / test side — the sample points and quadrature
  weights). The :class:`~orpheus.numerics.frame.Frame` so formed
  emits two operational faces:

  * the **analysis** face :math:`\Pi = T` — sampled values →
    coefficients (``frame.analysis``);
  * the **reconstruction** face :math:`R` — coefficients → values,
    the canonical-dual synthesis (``frame.reconstruction``).

  Together :math:`(R, \Pi)` define a Galerkin discretisation of any
  :math:`A : V \to V` as :math:`A_h = \Pi A R : W \to W`
  (Brenner & Scott 2008, §3.4).

- Two **disciplines** distinguish how :math:`(R, \Pi)` is
  constructed — equivalently, which **dual frame** supplies the
  reconstruction:

  * **Galerkin** — test space equals trial space (the inner product
    used to define :math:`\Pi` is taken in :math:`V` itself); the
    reconstruction is the **canonical dual**. The canonical case for
    :math:`L^2`-orthogonal moment projection. Invariant:
    :math:`\Pi^* = R` under the :math:`V`-inner product (up to the
    Gram diagonal — see the warning below).
  * **Petrov-Galerkin** — test space differs from trial space (the
    inner product on :math:`V` and the inner product used to define
    :math:`\Pi` are different); the reconstruction is an **oblique
    dual**. The canonical case for energy condensation (within-group
    spectrum weighting) and spatial homogenisation (flux-volume
    weighting). Sibling discipline of Galerkin.

- The discipline is a **property of the frame** (which dual it
  uses), not a separate mechanism. The forward-looking projection
  ABCs in :mod:`orpheus.numerics.projection` —
  :class:`~orpheus.numerics.projection.AnalysisOperator` (most
  general) →
  :class:`~orpheus.numerics.projection.GalerkinProjection` /
  :class:`~orpheus.numerics.projection.PetrovGalerkinProjection`
  (discipline) — carry the discipline *vocabulary* as a typed
  hierarchy, so a future concrete analysis operator signals its
  dual-frame discipline in its type.

- The single concrete frame shipping today is the
  **spherical-harmonic frame**
  :meth:`Quadrature.angular_frame(L)
  <orpheus.numerics.quadrature.Quadrature.angular_frame>` — the
  :class:`~orpheus.numerics.basis.SphericalHarmonicBasis` of order
  :math:`L` bound to a :math:`S^2` cubature. It is a **4π-tight
  frame** (Galerkin discipline). Future Petrov-Galerkin frames land
  with energy condensation (§17 of Grand Report v3) and spatial
  homogenisation (§18). The discipline ABCs are in place so the
  symmetry of the two disciplines is visible to future readers.

- Every concrete Galerkin frame satisfies the
  **idempotency-on-coefficients** invariant on a sufficiently-exact
  quadrature:

  .. math::

     \Pi \, R \;=\; c_{V}\,I_{W},

  where :math:`c_V` is a scalar that depends on the inner-product
  convention of :math:`V`. For the SN spherical-harmonic frame on a
  Lebedev quadrature, :math:`c_V = 4\pi` — this is the **L1
  idempotency** identity :eq:`pi-r-equals-4pi-i` verified at multiple
  :math:`L` against multiple Lebedev orders. (A 4π-tight frame is one
  whose frame operator :math:`S = T^*T` is :math:`4\pi` times the
  identity; the tightness constant IS this :math:`c_V`.)


The discrete frame — analysis, synthesis, and the frame operator
================================================================

The :math:`(R, \Pi)` pair is the language of **frame theory**
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

In the ORPHEUS algebra, the analysis operator IS the **projection**
:math:`\Pi = T`, and the **reconstruction** :math:`R` is the
**canonical-dual synthesis** — :math:`T^*` weighted by the dual
frame's Gram-inverse so that :math:`\Pi R` recovers the
band-limited identity (up to tightness). The bare :math:`T^* = S_0`
(the *naked synthesis*) is the shared
:meth:`~orpheus.numerics.basis.Basis.synthesize` primitive on the
:class:`~orpheus.numerics.basis.Basis`; the analysis face
:math:`\Pi` and the reconstruction face :math:`R` are each
:math:`S_0` post-multiplied by exactly one diagonal weight family
(the quadrature weight :math:`w_n` for analysis; the
addition-theorem factor :math:`2\ell+1` for reconstruction).

Given a fine space :math:`V` (e.g. :math:`L^2(S^2)` for the angular
flux) and a coarse coefficient space :math:`W` (e.g. polynomials of
degree :math:`\le L` on :math:`S^2`), a Galerkin discretisation
splits as:

.. math::
   :label: galerkin-pair

   R \;:\; W \to V, \qquad
   \Pi \;:\; V \to W,
   \qquad
   \Pi \, R \;=\; c_V \, I_W

where :math:`c_V` is the frame's tightness constant
(:math:`c_V = 1` in the fully-orthonormal case;
:math:`c_V = 4\pi` for the no-prefactor real spherical harmonics).

.. vv-status: galerkin-pair documented

The frame is fully determined by:

1. The **fine space** :math:`V` and its inner product
   :math:`\langle \cdot, \cdot \rangle_V` — fixed by the
   :class:`~orpheus.numerics.measure.DiscreteMeasure` (the analysis /
   test side: the sample nodes and quadrature weights). For SN
   angular flux, :math:`V = L^2(S^2)` and the inner product is the
   W-weighted discrete sum
   :math:`\langle f, g \rangle_W = \sum_n w_n f_n g_n` on the
   angular cubature.
2. The **basis** of :math:`W` — fixed by the
   :class:`~orpheus.numerics.basis.Basis` (the synthesis / trial
   side). For SN scattering, the basis is the real spherical
   harmonics :math:`\{Y_\ell^m\}_{\ell \le L}`.
3. The **discipline** — a property of which **dual frame** supplies
   the reconstruction. Galerkin (canonical dual): the test space is
   the same basis. Petrov-Galerkin (oblique dual): the test space is
   different (e.g. weighted by the within-group spectrum
   :math:`\chi_g(\phi_g)`).

Once the basis and the measure are bound by a
:class:`~orpheus.numerics.frame.Frame`, :math:`\Pi` and :math:`R`
are uniquely determined up to the :math:`c_V` normalisation.


The Galerkin discipline
=======================

The Galerkin discipline is characterised by **test space equals
trial space**. The defining identity is

.. math::
   :label: galerkin-self-adjoint

   \Pi^* \;=\; R
   \quad \text{(under the } V \text{ inner product, orthonormal basis)}.

.. vv-status: galerkin-self-adjoint documented

i.e. the projection's Hilbert adjoint is its reconstruction. This is
why a single basis :math:`\{e_k\}` produces both :math:`\Pi` and
:math:`R`:

.. math::
   :label: galerkin-construction

   (\Pi f)_k &\;=\; \langle e_k, f \rangle_V, \\
   R \, c     &\;=\; \sum_k c_k\,e_k.

.. vv-status: galerkin-construction documented

.. warning::

   The identity :math:`\Pi^* = R` holds when the basis
   :math:`\{e_k\}` is orthonormal in :math:`V`. When the basis is
   only orthogonal — the case for the no-:math:`4\pi/(2\ell+1)`-
   prefactor real spherical harmonics ORPHEUS uses — the strict
   Hilbert adjoint :math:`\Pi^*` and the addition-theorem
   reconstruction face :math:`R` differ by a **diagonal-in-:math:`\ell`
   scaling**. Specifically the strict adjoint is the *naked*
   synthesis (no :math:`(2\ell+1)` factor), while the frame's
   ``reconstruction`` face
   (:attr:`Frame.reconstruction
   <orpheus.numerics.frame.Frame.reconstruction>`) carries the
   :math:`(2\ell+1)` factor that the Pℓ scattering reconstruction
   needs:

   .. math::

      (\Pi^* c)_n
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
   ``frame.analysis.H`` is :math:`\Pi^* = g_C\,S_0`; and
   :meth:`frame.reconstruction.apply
   <orpheus.numerics.basis.Basis.reconstruct>` returns :math:`R`
   (with the :math:`(2\ell+1)` factor). The capability dishonesty
   that conflated the bare transpose with the Hilbert adjoint was
   caught by QA review and corrected as ERR-039 (see the project's
   L0 error catalog) — :math:`R` and :math:`\Pi^*` are both useful
   and coexist as distinct frame faces; they differ by exactly the
   Gram diagonal :math:`g_C = \mathrm{diag}(4\pi/(2\ell+1))`, so the
   docstrings and this page name them explicitly.

The Galerkin invariant :eq:`galerkin-pair` is then a consequence of
the basis being orthogonal in :math:`V`-inner-product. Concretely,
for the spherical-harmonic Galerkin pair:

.. math::

   (\Pi R)_{\ell m, \ell' m'}
   &\;=\; \sum_n w_n\,Y_\ell^m(\hat\Omega_n)\,
                       Y_{\ell'}^{m'}(\hat\Omega_n) \\
   &\;=\; \frac{4\pi}{2\ell+1}\,
          \delta_{\ell\ell'}\,\delta_{mm'},

so :math:`\Pi R = \mathrm{diag}(4\pi/(2\ell+1))`. Composing with the
reconstruction face's :math:`(2\ell+1)` factor (the addition-theorem
weight) yields :math:`\Pi R = 4\pi I` — the L1 identity that the
test
``tests/numerics/test_spherical_harmonic_space.py``
verifies at :math:`L = 2,\,3,\,4` (see :eq:`pi-r-equals-4pi-i` in
:ref:`spherical-harmonics`). This :math:`4\pi` is precisely the
frame's **tightness constant** :math:`c_V`: the frame operator
:math:`S = T^*T` equals :math:`4\pi\,I`, so the spherical-harmonic
frame is a 4π-tight frame.

.. note::

   The identity :math:`\Pi R = 4\pi I` is **not** identity-on-the-
   nose because the no-prefactor convention pushes the
   :math:`4\pi/(2\ell+1)` factor onto the orthogonality. A
   spherical-harmonic frame with a strict :math:`\Pi R = I`
   invariant could be built by dividing the analysis weights by
   :math:`4\pi`, but the project chose to absorb the factor at the
   reconstruction face (the :math:`(2\ell+1)` weight) so the
   addition theorem reads cleanly. See :ref:`spherical-harmonics`
   for the convention rationale.


The Petrov-Galerkin discipline
==============================

The Petrov-Galerkin discipline is characterised by **test space
differs from trial space**. The :math:`(R, \Pi)` pair is built from
two distinct bases — :math:`\{e_k\}` for the trial space (the
reconstruction basis) and :math:`\{f_k\}` for the test space (the
projection inner product):

.. math::
   :label: petrov-galerkin-construction

   (\Pi g)_k \;=\; \langle f_k, g \rangle_V, \qquad
   R \, c \;=\; \sum_k c_k\,e_k.

.. vv-status: petrov-galerkin-construction documented

The pair satisfies :math:`\Pi R = I_W` (by definition of
"projection"), but :math:`\Pi^* \ne R` — the test inner product
makes the adjoint distinct from the reconstruction.

The canonical Petrov-Galerkin pairs in reactor physics:

.. list-table:: Petrov-Galerkin pairs
   :header-rows: 1
   :widths: 22 26 26 26

   * - Use
     - Trial basis :math:`\{e_k\}`
     - Test basis :math:`\{f_k\}`
     - Reference
   * - Energy condensation
     - Indicator on broad group :math:`G`
     - Within-group spectrum
       :math:`\phi_g \cdot \mathbf{1}_{g \in G}`
     - Hébert 2009, §6.2
   * - Spatial homogenisation
     - Indicator on region :math:`R`
     - Region flux × volume
       :math:`(V_i \phi_i) \cdot \mathbf{1}_{i \in R}`
     - Smith 1986; Hébert 2009 §13
   * - Stochastic Galerkin
     - Polynomial-chaos basis (Hermite, Legendre)
     - Same basis under PCE inner product
     - Xiu & Karniadakis 2002

In all three cases the test basis encodes a **physical
weighting** — the importance of the fine-space slot from the
solver's perspective — so that the coarse coefficients faithfully
preserve reaction rates / flux-volume integrals / variance
moments.

No concrete Petrov-Galerkin subclasses ship in Wave 0 of the SN
performance plan; the
:class:`~orpheus.numerics.projection.PetrovGalerkinProjection` ABC
is in place to flag the design symmetry. Energy condensation
(Issue #17, future) will install the first concrete instance.


Cross-method consumer table
===========================

The Galerkin / Petrov-Galerkin pair is **infrastructure**, not
SN-only. Every method that lifts an angular / energy / spatial
axis between fine and coarse representations builds on these
primitives.

.. list-table:: Where the Galerkin pair is consumed
   :header-rows: 1
   :widths: 22 22 22 16 18

   * - Consumer
     - Discipline
     - Pair
     - Status
     - Reference
   * - **SN aniso scattering**
     - Galerkin
     - ``frame.analysis`` /
       ``frame.reconstruction``
       (``Quadrature.angular_frame(L)``)
     - **Live** (this work, Wave 1)
     - §9 (Grand Report v3 line 1230)
   * - PN solver
     - Galerkin
     - Same SH ``frame.analysis``
       on the moment-space side
     - Pending (PN solver not implemented)
     - §10 (lines 1299–1305)
   * - Energy condensation
     - Petrov-Galerkin
     - Within-group-spectrum weighting on test side, identity on
       trial side
     - Pending (concrete subclass needed)
     - §17 (line 3935); Hébert 2009 §6.2
   * - Spatial homogenisation
     - Petrov-Galerkin
     - Flux-volume weighting on test side, identity on trial side
     - Pending
     - §18
   * - Stochastic Galerkin
     - Galerkin
     - Polynomial-chaos basis on the random-input axis
     - Pending
     - §15A.7
   * - MC adjoint moments
     - Galerkin
     - Same SH ``frame.analysis``
       used as a variance-reduction estimator
     - Pending
     - Lewis & Miller 1993 §10
   * - SN sensitivity
     - Galerkin (adjoint)
     - Same pair, applied to the adjoint flux
     - Pending
     - Cacuci 2003

The two architectural payoffs:

* **One mechanism per discipline**, not one per consumer. The
  spherical-harmonic :class:`~orpheus.numerics.frame.Frame` emits the
  same ``analysis`` / ``reconstruction`` faces whether SN uses them
  for scattering or PN uses them for streaming — the difference is
  which axis the face is wrapped onto via the
  :class:`~orpheus.numerics.operator.TensorProductOperator`
  algebra (see :ref:`operator-algebra` and the tensor-product
  section there).
* **One V&V chain per discipline**. The Galerkin idempotency tests
  in :file:`tests/numerics/test_spherical_harmonic_space.py` cover
  every Galerkin consumer, not just SN. Energy condensation will
  inherit the Petrov-Galerkin-specific tests when its concrete
  frame lands.


The discipline vocabulary — type signal vs documentation
========================================================

The spherical-harmonic frame realises the analysis / reconstruction
pair *generically* (one :class:`~orpheus.numerics.frame.Frame`
emits both faces; the discipline — canonical-dual vs oblique-dual —
is a property of the frame). Alongside that mechanism, the
:mod:`orpheus.numerics.projection` module keeps a **forward-looking
discipline vocabulary**: a typed ABC hierarchy that lets a *future*
concrete analysis operator (energy condensation, spatial
homogenisation) signal its dual-frame discipline in its type. The
user's pedantic naming rule drives it: **a reader of a type name
should know its properties without reading the docstring.**

.. code-block:: python

   class AnalysisOperator(LinearOperatorMixin, ABC):
       """Most general analysis (projection) ABC — Π = T."""

   class GalerkinProjection(AnalysisOperator, ABC):
       """Galerkin discipline: test space = trial space
          (canonical dual)."""

   class PetrovGalerkinProjection(AnalysisOperator, ABC):
       """Petrov-Galerkin discipline: test space ≠ trial space
          (oblique dual)."""

A reader of ``class EnergyCondensation(PetrovGalerkinProjection)``
immediately knows:

1. It's an **analysis operator** (inherits from
   :class:`AnalysisOperator`).
2. It follows **Petrov-Galerkin** discipline (inherits from
   :class:`PetrovGalerkinProjection`, so :math:`\Pi^* \ne R` — an
   oblique dual).
3. The remaining axis structure (which fine/coarse spaces) is the
   class's own concern.

No docstring needed for those facts. Compare to the alternative
(rejected) flattening:

.. code-block:: python

   # REJECTED: collapses Galerkin/Petrov-Galerkin distinction
   class EnergyCondensation(AnalysisOperator):
       """Petrov-Galerkin condensation on the energy axis.

       (must be read to know the discipline)
       """

The rejected form makes the discipline a docstring claim instead of
a type claim. When a future reader reaches the energy-condensation
implementation and sees ``class EnergyCondensation(AnalysisOperator)``
they have no way to know — without reading the docstring — whether
it satisfies the Galerkin :math:`\Pi^* = R` invariant or the
oblique-dual Petrov-Galerkin one. The type hierarchy answers this
without reading prose.

The same rule drives the rename of:

* ``ProjectionOperator`` → ``AnalysisOperator`` — the ``Analysis``
  prefix names the frame-theory *analysis operator* :math:`T` the
  ABC abstracts (parity with the
  :class:`~orpheus.numerics.projection.ReconstructionOperator`
  sibling that abstracts the synthesis-side :math:`R`), and the
  ``Operator`` suffix signals "carries the operator algebra"
  (parity with :class:`~orpheus.numerics.operator.LinearOperator`).

.. note::

   The spherical-harmonic analysis and reconstruction are now the
   :class:`~orpheus.numerics.frame.Frame`'s ``analysis`` /
   ``reconstruction`` faces, NOT standalone operator classes. The
   discipline ABCs above are retained as the typed vocabulary the
   future Petrov-Galerkin frames will declare against; the concrete
   ``GalerkinProjection`` / ``PetrovGalerkinProjection`` subclasses
   ship with energy condensation (§17) and spatial homogenisation
   (§18).


Numerical evidence
==================

The L1-tagged tests in
:file:`tests/numerics/test_spherical_harmonic_space.py` verify the
Galerkin discipline's invariants on the spherical-harmonic frame's
``analysis`` / ``reconstruction`` faces:

1. **Idempotency** (4π-tightness):
   :math:`\Pi R c = 4\pi c` on band-limited
   coefficient input, verified at :math:`L = 2,\,3,\,4` against
   Lebedev orders :math:`7,\,13,\,17`. See
   :eq:`pi-r-equals-4pi-i` in :ref:`spherical-harmonics`.
2. **Adjoint pairing**:
   :math:`\langle \Pi \psi, c \rangle = \langle \psi, \Pi^* c \rangle_W`
   with :math:`\Pi^* = g_C\,S_0`, verified to ``rtol=1e-12`` on a
   Lebedev order-13 grid at :math:`L = 3`.

The tests verify mathematical identities of the operator algebra
(V&V level L1 — equation verification by analytical reference). The
companion L0/foundation shape and capability tests verify software
invariants (frame face spaces, capability sets) and are tagged
accordingly.


Implementation map
==================

* :class:`~orpheus.numerics.frame.Frame` — the discrete frame that
  binds a :class:`~orpheus.numerics.basis.Basis` to a
  :class:`~orpheus.numerics.measure.DiscreteMeasure` and emits the
  ``analysis`` (:math:`\Pi = T`) and ``reconstruction`` (:math:`R`)
  faces. The single mechanism for every choice-dependent
  change-of-basis.
* :class:`~orpheus.numerics.basis.Basis` — the synthesis (trial)
  side ABC: tabulate, naked synthesis :math:`S_0`, the three
  weighted contractions, and the discrete Gram.
* :class:`~orpheus.numerics.basis.SphericalHarmonicBasis` — the
  first concrete basis (real spherical harmonics); carries the
  no-prefactor convention and the
  :attr:`~orpheus.numerics.basis.SphericalHarmonicBasis.addition_theorem_factor`
  :math:`(2\ell+1)`.
* :class:`~orpheus.numerics.projection.AnalysisOperator` — the
  most-general analysis (projection) ABC; the forward-looking
  discipline vocabulary.
* :class:`~orpheus.numerics.projection.GalerkinProjection` — the
  Galerkin discipline ABC (canonical dual).
* :class:`~orpheus.numerics.projection.PetrovGalerkinProjection`
  — the Petrov-Galerkin discipline ABC (oblique dual); sibling.
* :class:`~orpheus.numerics.projection.ReconstructionOperator` —
  the reconstruction-side ABC :math:`R : W \to V`, sibling of
  :class:`~orpheus.numerics.projection.AnalysisOperator`.
* :meth:`Quadrature.angular_frame(L)
  <orpheus.numerics.quadrature.Quadrature.angular_frame>` — builds
  the order-:math:`L` spherical-harmonic frame on a quadrature; the
  single home of the :math:`S^2` embedding.

The full-space projector — the operator that projects the SN
:math:`(N, n_x, n_y, n_g)` angular flux onto the
:math:`(L+1, 2L+1, n_x, n_y, n_g)` moment field — is built as a
**tensor product** of the angular-axis analysis face :math:`\Pi`
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


History — from operator classes to the discrete frame
=====================================================

The spherical-harmonic projection and reconstruction were first
shipped (Wave 0 of the SN performance plan) as standalone operator
classes ``HarmonicMomentProjection`` / ``HarmonicMomentReconstruction``
under a three-level inheritance
(``ProjectionOperator`` → ``GalerkinProjection`` → concrete). Two
naming-audit corrections then established the discipline-must-be-
typed pedagogy that survives in the discipline-ABC vocabulary above.

The Frame/Basis carve (``refactor/operator-inverse-algebra``)
took the next step: the projection :math:`M = Y^*W` and the
addition-theorem reconstruction :math:`R = (2\ell+1)\,S_0` are NOT
two unrelated operator classes — they are the **two faces of one
discrete frame** binding the SH basis to the angular measure. The
standalone operator classes were retired into
:class:`~orpheus.numerics.frame.Frame`:

* ``HarmonicMomentProjection`` → ``frame.analysis``
  (:attr:`Frame.analysis <orpheus.numerics.frame.Frame.analysis>`,
  the analysis face :math:`M = T`);
* ``HarmonicMomentReconstruction`` → ``frame.reconstruction``
  (:attr:`Frame.reconstruction
  <orpheus.numerics.frame.Frame.reconstruction>`, the
  reconstruction face :math:`R`);
* the :math:`(2\ell+1)` factor moved onto
  :attr:`SphericalHarmonicBasis.addition_theorem_factor
  <orpheus.numerics.basis.SphericalHarmonicBasis.addition_theorem_factor>`
  (one home for the SH convention).

The architectural payoff: one mechanism (the frame) instead of two
operator classes, and the **discipline becomes a property of the
frame** (which dual it uses) rather than a fact baked into a class
name. The discipline ABCs are kept as the typed vocabulary future
Petrov-Galerkin frames declare against.


References
==========

* Brenner, S. C. and Scott, L. R. (2008). *The Mathematical Theory
  of Finite Element Methods*, 3rd ed. Springer. §3.4 (Galerkin /
  Petrov-Galerkin general framework).
* Christensen, O. (2016). *An Introduction to Frames and Riesz
  Bases*, 2nd ed. Birkhäuser. (The analysis operator :math:`T`, the
  synthesis operator :math:`T^*`, the frame operator
  :math:`S = T^*T`, tight frames, and the canonical dual — the
  harmonic-analysis foundation of the
  :class:`~orpheus.numerics.frame.Frame` abstraction.)
* Bell, G. I. and Glasstone, S. (1970). *Nuclear Reactor Theory*.
  Van Nostrand Reinhold. §1.6 (spherical-harmonic moment
  projection in transport).
* Lewis, E. E. and Miller, W. F. Jr. (1993). *Computational Methods
  of Neutron Transport*. ANS. §4.7 (Pℓ Galerkin reconstruction with
  the :math:`(2\ell+1)` factor).
* Hébert, A. (2009). *Applied Reactor Physics*. Polytechnique. §6.2
  (energy condensation as a Petrov-Galerkin projection).
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
