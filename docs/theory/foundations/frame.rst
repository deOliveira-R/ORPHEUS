.. _galerkin-projection:

=========================================
The Discrete Frame — projection machinery
=========================================

Every method in ORPHEUS that transitions a function between a **fine**
representation and a **coarse** one — discrete-ordinate :term:`angular flux` to
spherical-harmonic moments, fine-energy fluxes to broad-group cross
sections, regional flux to homogenised cross sections — does so via a
**(reconstruction, analysis)** pair of linear operators :math:`(R, M)`.
In harmonic-analysis language this pair is exactly the two operational
**faces of a discrete frame**: the analysis operator :math:`M = T`
(sampled values → coefficients) and the reconstruction :math:`R`
(coefficients → values, the canonical-dual synthesis). This page is the
canonical home for that pair, for the discrete-**frame** abstraction that
realises it (:class:`~orpheus.numerics.frame.FrameBase` and its
discipline subclasses), and for the **discipline** — Galerkin vs
Petrov-Galerkin, test space equal to vs different from trial space — that
distinguishes its variants.

The page is organised **general case first**. The Petrov-Galerkin frame
(test :math:`\ne` trial) is the general object; the Galerkin frame (test
:math:`=` trial) is its symmetric specialisation — a structure the code
mirrors exactly, since :class:`~orpheus.numerics.frame.GalerkinFrame`
**is-a** :class:`~orpheus.numerics.frame.PetrovGalerkinFrame`
(:class:`~orpheus.numerics.frame.FrameBase` →
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame` →
:class:`~orpheus.numerics.frame.GalerkinFrame`). After the abstract frame
theory come the two disciplines and their concrete consumers — spatial
homogenisation and energy condensation (Petrov-Galerkin), spherical-
harmonic scattering projection (Galerkin) — and then the **advanced**
material: eigenbasis ownership, the cross-method consumer catalog, and
the adjoint-weighted seam. Binding a
:class:`~orpheus.numerics.basis.Basis` to a
:class:`~orpheus.numerics.measure.DiscreteMeasure` through a single
frame puts one mechanism in front of every consumer, instead of a
separate projection / reconstruction operator class per method.

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
   :ref:`sn-homogenization-petrov-galerkin-frame` (below); it was
   rewritten to this same
   Petrov-Galerkin framing under Issue #268 (the earlier
   ":math:`L^2(\phi V)`-Galerkin" reading is retired there, with the
   forward-flux metric-fold shown to be the Galerkin *degenerate* of the
   eigenvalue-consistent adjoint-weighted case).

.. note:: **What shipped since (P3 / P5 / P7).** The
   :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` base was empty of
   concrete consumers at the P1 carve; the **forward** (reaction-rate,
   :math:`\varphi^* = \varphi`) homogenisation (P3) and energy
   condensation (P5) have since shipped as concrete instances
   (:meth:`Solution.homogenize
   <orpheus.sn.solution.Solution.homogenize>`,
   :meth:`Solution.condense <orpheus.sn.solution.Solution.condense>`).
   This page (P7) is the capstone that ties the discipline-type
   hierarchy, the composed-operator verbs
   (:ref:`frame-composed-verbs`), the three-way Gram-structure gate
   (:ref:`frame-least-squares-discipline`), and the eigenbasis-ownership
   ruling (:ref:`frame-eigenbasis-ownership`) into one narrative, and
   reconciles its consumer table with the shipped frames. The one
   remaining frame discipline that is **theory-documented but not built**
   is the **eigenvalue-consistent** (adjoint-weighted) projection
   (:ref:`frame-adjoint-weighted-seam`, phase P6).

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
  and a fixed :math:`L^2` metric (the :term:`quadrature` weights). The
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

- **Two families of concrete frame ship today.** The
  **spherical-harmonic frame**
  :meth:`Quadrature.angular_frame(L)
  <orpheus.numerics.quadrature.Quadrature.angular_frame>` — the
  :class:`~orpheus.numerics.basis.SphericalHarmonicBasis` of order
  :math:`L` bound to an :math:`S^2` cubature — is a
  :class:`~orpheus.numerics.frame.GalerkinFrame` (``test is trial``)
  and a **4π-tight frame**. The forward (reaction-rate)
  **Petrov-Galerkin** consumers —
  :meth:`Solution.homogenize <orpheus.sn.solution.Solution.homogenize>`
  (space) and
  :meth:`Solution.condense <orpheus.sn.solution.Solution.condense>`
  (energy) — ship as concrete
  :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` instances with an
  explicit flux- / spectrum-weighted
  :class:`~orpheus.numerics.basis.WeightedIndicatorBasis` test basis
  (landed P3 / P5; full derivations in
  :ref:`sn-homogenization-petrov-galerkin-frame` and
  :ref:`sn-energy-condensation`, below). The
  **eigenvalue-consistent** (adjoint-weighted,
  :math:`\varphi^* \ne \varphi`) case is the **documented-future seam**:
  the theory is written (:eq:`sn-homogenization-bilinear`), the
  implementation (P6) is blocked on the adjoint flux
  :math:`\varphi^*` (:ref:`frame-adjoint-weighted-seam`).

- **The frame is also the production COMPOSER, not only the two faces.**
  Beyond ``analysis`` / ``reconstruction``, a
  :class:`~orpheus.numerics.frame.FrameBase` emits the composed-operator
  verbs a consumer uses directly — *define a frame, compose, done*:
  :meth:`conjugate <orpheus.numerics.frame.FrameBase.conjugate>`
  (:math:`R\circ A\circ M`, the scattering kernel) and
  :meth:`project <orpheus.numerics.frame.FrameBase.project>`
  (:math:`G^{-1}M`, the homogenise / condense verb). These are typed
  operator products, not hand-rolled numpy chains
  (:ref:`frame-composed-verbs`).

- **Three disciplines, gated by the trial basis's Gram structure.** The
  built frames cover the *row-sum-collapsible* Gram cases —
  :class:`~orpheus.numerics.frame.GalerkinFrame` (diagonal Gram,
  ``test is trial``) and the forward
  :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` (diagonal *or*
  partition-of-unity Gram). The third discipline — a least-squares
  frame over a **dense** cross-Gram needing the real
  :math:`(MR)^{-1}M` solve — is **designed but not built**:
  :meth:`FrameBase.project <orpheus.numerics.frame.FrameBase.project>`
  *refuses* a :class:`~orpheus.numerics.basis.GramStructure`
  ``DENSE`` trial (raising
  :class:`~orpheus.numerics.operator.NotInvertible`) rather than
  return a silently-wrong coarsening (:ref:`frame-least-squares-discipline`).

- Every concrete :class:`~orpheus.numerics.frame.GalerkinFrame`
  satisfies the **idempotency-on-coefficients** invariant on a
  sufficiently-exact quadrature:

  .. math::
     :label: galerkin-frame-idempotency

     M \, R \;=\; c_{V}\,I_{W},

  .. (vv-status rationale) Structural invariant: the general Galerkin
     idempotency-on-coefficients schema M R = c_V I. Its SN concrete instance
     (c_V = 4π) is :eq:`pi-r-equals-4pi-i`, the L1-verified form — the canonical
     pin ``tests/numerics/test_spherical_harmonic_space.py`` ``verifies("pi-r-equals-4pi-i")``
     constructs Π R = 4π I at multiple L / Lebedev orders. Not a separate solver claim.
  .. vv-status: galerkin-frame-idempotency documented

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


The Petrov-Galerkin frame
=========================

The **general** discrete frame is Petrov-Galerkin: the *test* functions
that measure the residual need not equal the *trial* functions that
reconstruct it. The code mirrors this generality — the Galerkin case is a
*subclass* (:class:`~orpheus.numerics.frame.GalerkinFrame` **is-a**
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame` with
``test is trial``) — so the Petrov-Galerkin frame is presented first, as
the general object, and the Galerkin frame (:ref:`the special case
<frame-galerkin-frame>`) as its symmetric specialisation.

In general
----------

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

Posing the adjoint
------------------

It is tempting to absorb the test weight (:math:`\phi`,
:math:`\varphi^*`) into the *measure* and call the result an orthogonal
(Galerkin) projection in an :math:`L^2(\phi V)` metric. That fold is
legitimate **only** for the forward-flux, reaction-rate-only row of the
Petrov-Galerkin table above: there the test weight and the integrand
multiplier coincide (both are the forward flux :math:`\phi`), and the two
readings are the *same* map — the **Galerkin degenerate**. It **breaks**
for the eigenvalue-consistent row, whose preserved functional is the
**bilinear** form :math:`\langle \varphi^*, \Sigma\,\phi\rangle`: the
test weight is the **adjoint** :math:`\varphi^*`, the integrand is the
**forward** flux :math:`\phi` — different functions — so no single metric
on the measure reproduces it (folding either one into the metric
mis-weights the other). The discipline must therefore live on the **test
side** (the frame type), and the measure stays a fixed :math:`L^2`
metric; only then does the adjoint fall out naturally — the adjoint
problem swaps the test basis to :math:`\varphi^*`, and the oblique dual
:math:`M^* \ne R` (not the canonical dual) is exactly what that swap
needs.

This is why the architecture is as it is. The frame was in fact first
posed as a *Galerkin* projection with the flux folded into the volume
measure; it was the adjoint-weighted requirement — the need to keep
:math:`\keff` stationary, which first-order perturbation theory ties to
the **adjoint**-weighted residual, not the forward-weighted one — that
forced the re-posing as **Petrov-Galerkin**, with the weighting on an
explicit test basis rather than smuggled onto the measure. Folding the
solution into the metric is precisely the mistake the #268 ruling
forbids: *the measure carries the axis and the fixed* :math:`L^2`
*metric, never the discipline.*

The forward (Galerkin-degenerate) Petrov-Galerkin frames have since
shipped (P3 / P5): :meth:`Solution.homogenize
<orpheus.sn.solution.Solution.homogenize>` and
:meth:`Solution.condense <orpheus.sn.solution.Solution.condense>` build a
concrete :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` with an
explicit flux- / spectrum-weighted
:class:`~orpheus.numerics.basis.WeightedIndicatorBasis` test basis — the
first concrete instances landing exactly as a ``test_basis`` choice on
the existing mechanism, not a new mechanism. The forward-flux derivation
is worked in full in :ref:`sn-homogenization-petrov-galerkin-frame`
(:ref:`§2c <sn-spatial-homogenization>`). The non-degenerate
(:math:`\varphi^* \ne \varphi`) eigenvalue-consistent case is the
**documented-future seam**: the bilinear derivation
(:eq:`sn-homogenization-bilinear`,
:ref:`sn-homogenization-why-petrov-galerkin`) is written, its
implementation (P6) blocked on the adjoint flux :math:`\varphi^*`
(:ref:`frame-adjoint-weighted-seam`). This subsection fixes the
*architecture* of where the weighting lives.

.. _sn-spatial-homogenization:

Applied to spatial homogenization
---------------------------------

Once a fine-mesh solution :math:`\phi_{i,g}` is in hand, a coarse-mesh
model that **reproduces every reaction rate** of the fine model can be
built by collapsing the fine cross sections onto the coarse cells. This
is *spatial homogenization* — a domain operation on the solution, not a
solver step. It is the spatial sibling of energy *condensation* (group
collapse); the two together are the classical "smear the detail you
have resolved into effective constants for a coarser calculation" move
(Hébert, *Applied Reactor Physics*, §13 for space, §6.2 for energy).

.. note::

   This is the **space-only** slice (the spatial sibling of energy
   condensation). It is **dimension-agnostic** — 1-D and 2-D fine meshes
   flow through the one frame body, because the coarse cell-indicator
   basis and the fine volume measure are n-D (see
   :ref:`sn-homogenization-petrov-galerkin-frame`). Energy is *not*
   condensed — the group structure (:math:`eg`) carries through unchanged
   — and the coarse mesh must share the fine mesh's outer boundary with
   internal coarse edges aligned to fine-cell edges (each coarse cell is
   a contiguous union of fine cells). The asymmetry between
   homogenization and condensation, and *why* they return different
   types, is the subject of :ref:`sn-condense-homogenize-asymmetry`.

The defining property: reaction-rate preservation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Homogenization is defined by *what it must preserve*, not by an
averaging recipe chosen for convenience. The physical quantity a
transport calculation actually consumes is the **volume-integrated
reaction rate** in each region and group,

.. math::
   :label: sn-homogenization-fine-rate

   r_{R,g} \;=\; \sum_{i \in R} V_i\,\Sigma_{i,g}\,\phi_{i,g},

.. (vv-status rationale) Definitional identity: the fine-mesh
   volume-integrated reaction rate. A derivation-decomposition step of
   the rate-preservation claim (:eq:`sn-homogenization-rate-preservation`),
   which is the verifiable solver claim (L0 gate
   ``tests.sn.test_homogenization``); this is its premise, not a separate
   claim.
.. vv-status: sn-homogenization-fine-rate documented

the sum over the fine cells :math:`i` contained in coarse cell
:math:`R`, with :math:`V_i` the fine-cell volume,
:math:`\phi_{i,g}` the converged fine :term:`scalar flux`, and
:math:`\Sigma_{i,g}` the fine cross section for whatever channel (total,
capture, fission, …) is in question. The coarse model carries one
effective cross section :math:`\Sigma_{R,g}` and one region flux
:math:`\Phi_{R,g}` per cell; for it to reproduce the fine rate it must
satisfy

.. math::
   :label: sn-homogenization-rate-preservation

   \Sigma_{R,g}\,\Phi_{R,g}
   \;=\;
   \sum_{i \in R} V_i\,\Sigma_{i,g}\,\phi_{i,g}.

Equation :eq:`sn-homogenization-rate-preservation` is the **only**
constraint homogenization imposes; everything below is derived from it.

The coarse region flux is fixed first, by the requirement that the
*production-free* particle inventory match — the region flux is the
flux integrated over the region,

.. math::
   :label: sn-homogenization-region-flux

   \Phi_{R,g} \;=\; \sum_{i \in R} V_i\,\phi_{i,g}.

.. (vv-status rationale) Representational identity: the coarse region
   flux is the flux integrated over the region. A definition consumed by
   the rate-preservation claim (:eq:`sn-homogenization-rate-preservation`),
   not an independent solver claim.
.. vv-status: sn-homogenization-region-flux documented

Substituting :eq:`sn-homogenization-region-flux` into
:eq:`sn-homogenization-rate-preservation` and solving for the effective
cross section gives the **flux·volume-weighted average**

.. math::
   :label: sn-homogenization-vector-collapse

   \Sigma_{R,g}
   \;=\;
   \frac{\sum_{i \in R} V_i\,\phi_{i,g}\,\Sigma_{i,g}}
        {\sum_{i \in R} V_i\,\phi_{i,g}}
   \;=\;
   \frac{\sum_{i \in R} w_{i,g}\,\Sigma_{i,g}}
        {\sum_{i \in R} w_{i,g}},
   \qquad
   w_{i,g} \equiv V_i\,\phi_{i,g}.

.. (vv-status rationale) Derivation-decomposition step: the
   flux·volume-weighted collapse obtained by solving the
   rate-preservation identity (:eq:`sn-homogenization-rate-preservation`)
   for the effective cross section. The verifiable content is the L0
   rate-preservation gate; this is the algebraic rearrangement, not a
   separate claim.
.. vv-status: sn-homogenization-vector-collapse documented

The weight :math:`w_{i,g} = V_i\,\phi_{i,g}` is the **flux·volume**
of the fine cell — the same quantity that appears in both the numerator
(rate) and the denominator (flux integral), so the average is a genuine
convex combination of the fine values: :math:`\Sigma_{R,g}` is bracketed
by the region's fine-cell extremes
:math:`\min_{i\in R}\Sigma_{i,g} \le \Sigma_{R,g} \le \max_{i\in R}\Sigma_{i,g}`.
This is *not* a separate design choice — it falls straight out of
:eq:`sn-homogenization-rate-preservation`. Choosing any other weight
(volume-only, unweighted) would break rate preservation at material
interfaces, which is exactly the regime homogenization exists to handle.

The flux·volume weight :math:`w_{i,g}` is the operation's whole signal,
and it is *not* a free parameter: it is the **test weighting** that rate
preservation forces. That is what the next subsection makes precise —
homogenization is the coefficient extraction of a **Petrov-Galerkin**
frame whose *test* basis is the flux-weighted cell indicator and whose
*trial* basis is the plain cell indicator, and ORPHEUS realises it by
routing through the one discrete
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame`, not a bespoke
membership matmul (see :ref:`sn-homogenization-petrov-galerkin-frame`).

The vector channels — total, capture, leakage-loss
(:math:`\Sigma_L`), fission, and production
(:math:`\nu\Sigma_f`) — each collapse through
:eq:`sn-homogenization-vector-collapse` with the *same* per-:math:`(R,g)`
weight. A group with zero region flux
(:math:`\Phi_{R,g} = 0`) has no reaction rate to preserve, so its
effective cross section is set to zero — the :math:`0/0` of
:eq:`sn-homogenization-vector-collapse` resolved by the only physically
meaningful value.

The matrix channels weight by the *source* group
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The scattering matrices :math:`\Sigma_{s,\ell}[g',g]` (one per Legendre
order :math:`\ell`) and the :math:`(n,2n)` matrix
:math:`\Sigma_{2n}[g',g]` carry **two** group indices, stored
``[g_from, g_to]`` (the ORPHEUS scattering convention — see
:ref:`theory-cross-section-data`). A naïve reuse of the scalar weight
would be wrong: the reaction rate that an out-scatter channel
:math:`g' \to g` actually produces is

.. math::
   :label: sn-homogenization-scatter-rate

   r_{R}^{\,g'\to g}
   \;=\;
   \sum_{i \in R} V_i\,\phi_{i,g'}\,\Sigma_{s,\ell,i}[g',g],

.. (vv-status rationale) Definitional identity: the out-scatter rate of
   a matrix channel, driven by the source-group flux. The premise that
   forces the source-group weighting; the verifiable claim is the L0
   matrix-channel rate-preservation gate
   (``test_rate_preservation_scattering_and_n2n``).
.. vv-status: sn-homogenization-scatter-rate documented

driven by the population of the **source** group :math:`g'` — the group
whose flux scatters *out*. The number of :math:`g'\to g` events scales
with how many particles are *in* :math:`g'`, i.e. with
:math:`\phi_{i,g'}`, not with the sink-group flux
:math:`\phi_{i,g}`. Rate preservation
:eq:`sn-homogenization-rate-preservation` therefore demands that the
effective matrix entry be weighted by the source-group flux·volume:

.. math::
   :label: sn-homogenization-matrix-collapse

   \Sigma_{s,\ell,R}[g',g]
   \;=\;
   \frac{\sum_{i \in R} V_i\,\phi_{i,g'}\,\Sigma_{s,\ell,i}[g',g]}
        {\sum_{i \in R} V_i\,\phi_{i,g'}}
   \;=\;
   \frac{\sum_{i \in R} w_{i,g'}\,\Sigma_{s,\ell,i}[g',g]}
        {\Phi_{R,g'}},

.. (vv-status rationale) Derivation-decomposition step: the
   source-group-weighted matrix collapse — the matrix-channel analogue
   of :eq:`sn-homogenization-vector-collapse`. The verifiable content is
   the L0 matrix-channel rate-preservation gate (which catches a
   g_from↔g_to swap, vv Mode 2); this is the algebraic form, not a
   separate claim.
.. vv-status: sn-homogenization-matrix-collapse documented

so the weight :math:`w_{i,g'} = V_i\,\phi_{i,g'}` rides the **first**
(``g_from``) matrix axis. In the code this falls out of the *test side*
for free: the :class:`~orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis`
carries the per-group flux :math:`\phi` as its test weight, and its
**leading-aligned broadcast** aligns that weight's group axis to the
*first* trailing (``g_from``) axis of whatever field it analyses — a
vector channel weights elementwise, a ``[g_from, g_to]`` matrix channel
weights by its source group — *before* the region integral. The
:math:`1/\Phi_{R,g'}` normalisation is the frame's diagonal Gram
(:meth:`FrameBase.project <orpheus.numerics.frame.FrameBase.project>`'s
:meth:`FunctionSpace.apply_inverse_metric
<orpheus.numerics.space.FunctionSpace.apply_inverse_metric>`), whose
:math:`\Phi_{R,g'}` rides the ``g_from`` axis and broadcasts over the
trailing ``g_to`` axis. The :math:`(n,2n)` channel collapses identically
— same source-group weighting on its ``[g_from, g_to]`` layout. Both ride
the *same* ``sigma_frame`` because
:meth:`MaterialXSField.project_through
<orpheus.transport.mesh.material_xs_field.MaterialXSField.project_through>`
routes every rate-bearing channel through it. The mechanism that carries
this — the discrete :class:`~orpheus.numerics.frame.PetrovGalerkinFrame`
— is derived in :ref:`sn-homogenization-petrov-galerkin-frame`.

.. warning::

   The source-group weighting is the subtle point of the whole
   operation and a textbook variable-swap trap (vv-principles failure
   **Mode 2**, ``SigS`` vs ``SigS^T``). Weighting the matrix collapse by
   the **sink** group :math:`g` instead of the **source** group
   :math:`g'` produces an effective scattering matrix that does *not*
   preserve the out-scatter rate — a bug that is invisible on a
   single-material or flat-flux region (where every group's weight is
   proportional) and only bites on a heterogeneous, multi-group region
   with an asymmetric flux spectrum. The regression gate
   (:ref:`sn-homogenization-verification`) catches it precisely because
   its reference loop weights by ``g_from`` *explicitly*.

The fission spectrum is production-weighted
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The emission spectrum :math:`\chi_g` is **not** a reaction rate — it is
a probability distribution (a simplex,
:math:`\sum_g \chi_g = 1`; see :eq:`emission-spectrum-simplex` in
:ref:`theory-cross-section-data`). Flux·volume-weighting it would not
preserve anything physical and could leave the simplex. The
rate-preserving choice is to weight :math:`\chi` by each fine cell's
**fission production rate**

.. math::
   :label: sn-homogenization-production-weight

   p_i \;=\; \sum_g \nu\Sigma_{f,i,g}\,\phi_{i,g}\,V_i,

.. (vv-status rationale) Definitional identity: the per-cell fission
   production rate used as the χ-mixing weight. A premise of the
   χ-collapse, not a separate solver claim; the simplex/null-law content
   it feeds is verified by the data-layer ``Mixture`` invariant tests.
.. vv-status: sn-homogenization-production-weight documented

so that the homogenized spectrum is the production-weighted convex
average

.. math::
   :label: sn-homogenization-chi-collapse

   \chi_{R,g}
   \;=\;
   \frac{\sum_{i \in R} p_i\,\chi_{i,g}}
        {\sum_{i \in R} p_i}.

.. (vv-status rationale) Representational identity: the
   production-weighted convex average of the fine emission spectra — the
   spatial analogue of the multi-fissile χ_mix
   (:eq:`emission-spectrum-chi-mix`). The simplex/null law it must
   satisfy is the data-layer invariant verified by
   ``test_homogenized_chi_is_simplex`` + ``Mixture.__post_init__``; this
   label is the mixing formula, not a separate solver claim.
.. vv-status: sn-homogenization-chi-collapse documented

Because each fine :math:`\chi_i` is a probability simplex and the
weights :math:`p_i \ge 0`, :math:`\chi_R` is a **convex combination of
simplices, hence itself a simplex** — it is exactly the spatial analogue
of the production-weighted multi-fissile mixing
:math:`\chi_{\rm mix}` of :eq:`emission-spectrum-chi-mix` (where the
weights are per-isotope production rather than per-cell). The simplex /
null law is re-validated when the homogenized
:class:`~orpheus.data.macro_xs.mixture.Mixture` is constructed
(:meth:`Mixture.__post_init__
<orpheus.data.macro_xs.mixture.Mixture.__post_init__>`); a coarse cell
with no fissile fine cells (:math:`\sum_i p_i = 0`) gets
:math:`\chi_R = 0`, the null-law branch.

Balance is preserved cell-by-cell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The definitional total-XS balance every :class:`Mixture` carries,

.. math::
   :label: sn-homogenization-balance

   \Sigma_t
   \;=\;
   \Sigma_c + \Sigma_L + \Sigma_f
   + \operatorname{rowsum}(\Sigma_{s0})
   + \operatorname{rowsum}(\Sigma_{2n}),

.. (vv-status rationale) Literature-transcribed definition: the total-XS
   balance identity every Mixture carries (the same identity as
   :eq:`sigT-computed`). Its preservation under homogenization is gated
   by ``test_homogenized_materials_balance`` (which calls
   ``Mixture.assert_balanced``); the identity itself is a data-layer
   definition, not an SN solver claim.
.. vv-status: sn-homogenization-balance documented

survives the collapse **cell-by-cell** when the fine materials balance.
The argument is one line once the weighting is understood: fix a coarse
cell :math:`R` and group :math:`g`. Every *removal* channel on the
left- and right-hand sides of :eq:`sn-homogenization-balance` —
:math:`\Sigma_t,\ \Sigma_c,\ \Sigma_L,\ \Sigma_f`, and the row-sums
:math:`\sum_{g'}\Sigma_{s0}[g,g']`, :math:`\sum_{g'}\Sigma_{2n}[g,g']`
— is a *removal from group* :math:`g`, so each collapses with the
**same** weight :math:`w_{i,g} = V_i\phi_{i,g}` (the row-sum of a
``[g_from, g_to]`` matrix over its sink index :math:`g'` is a removal
*from* the source group :math:`g`, weighted by :math:`g`'s flux — the
source-group weighting of :eq:`sn-homogenization-matrix-collapse`
restricted to a diagonal-of-the-source row). Because every term shares
the one weight, the homogenized balance is the *same convex average* of
the fine balances:

.. math::
   :label: sn-homogenization-balance-preservation

   \Sigma_{t,R,g}
   - \Big(\Sigma_{c,R,g} + \Sigma_{L,R,g} + \Sigma_{f,R,g}
   + \operatorname{rowsum}(\Sigma_{s0,R})_g
   + \operatorname{rowsum}(\Sigma_{2n,R})_g\Big)
   \;=\;
   \frac{\sum_{i\in R} w_{i,g}\,\big(\text{fine balance residual}_{i,g}\big)}
        {\sum_{i\in R} w_{i,g}}
   \;=\; 0,

since each fine residual is zero. No separate "rebalance the homogenized
total" step is needed — preservation is automatic, and the homogenized
``Mixture`` passes :meth:`Mixture.assert_balanced
<orpheus.data.macro_xs.mixture.Mixture.assert_balanced>`. (Had the
vector channels and the matrix row-sums collapsed with *different*
weights, the balance would break — which is another way of seeing why
the source-group weighting of :eq:`sn-homogenization-matrix-collapse`
is forced, not chosen.)

.. _sn-homogenization-petrov-galerkin-frame:

Homogenization is a Petrov-Galerkin projection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Everything above derives the flux·volume average from rate
preservation. This subsection takes the second view — the one that
fixes *what kind of operator* homogenization is, and therefore *how it
is implemented*. The answer is a single sentence that the rest of the
subsection unpacks:

  Homogenization is the coefficient extraction :math:`G^{-1} M` of a
  **Petrov-Galerkin** frame: the *trial* basis is the plain coarse-cell
  indicator :math:`\mathbf{1}_R`, the *test* basis is the
  flux-weighted indicator :math:`\chi_R = \phi\,\mathbf{1}_R`, and the
  measure is the bare geometric volume measure :math:`\mu_V`.

This is not decoration. It is the reason the production code routes
:meth:`Solution.homogenize <orpheus.sn.solution.Solution.homogenize>`
through the *same* discrete :class:`~orpheus.numerics.frame.FrameBase`
abstraction that carries SN anisotropic-scattering moment projection —
one mechanism for every fine→coarse change of representation (Cardinal
Rule 2, single source of truth), instead of a bespoke membership
matmul per method. It is the consumer the discrete-frame theory page
points at as the headline **Petrov-Galerkin** instance
(Issue #268); the test functions differ
from the trial functions (:math:`\chi_R = \phi\,\mathbf{1}_R \ne
\mathbf{1}_R`), so the discipline is genuinely Petrov-Galerkin, carried
by the frame **type**
(:class:`~orpheus.numerics.frame.PetrovGalerkinFrame`).

.. warning::

   **This corrects an earlier draft of this section.** A previous
   version argued homogenization was the ":math:`L^2(\phi V)`-orthogonal
   **Galerkin** projection" — that the flux multiplier could be folded
   into the *measure* (read :math:`\langle\Sigma,\phi\,\mathbf{1}_R
   \rangle_{\mathrm{d}V} = \langle\Sigma,\mathbf{1}_R\rangle_{\phi V}`),
   making test and trial the *same* span in a flux-weighted metric. That
   reading is **forward-flux, reaction-rate-only**, and it
   structurally breaks for the eigenvalue-consistent homogenization
   reactor physics actually requires (see
   :ref:`sn-homogenization-why-petrov-galerkin` below). Folding the
   solution into the metric is precisely the mistake the #268 ruling
   forbids: *the measure carries the axis and the fixed* :math:`L^2`
   *metric, never the discipline.* The flux is a **test-weighting the
   solution emits**, living on the test side — the frame type — not on
   the geometry's measure.

The trial space, the test space, and the cross-Gram
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Let :math:`\Sigma(x)` be the fine cross-section field — a function on
the spatial domain, piecewise-constant on the *fine* cells. The coarse
model can only carry one value per coarse cell :math:`R`, so the
**trial** space — where the answer lives — is

.. math::
   :label: sn-homogenization-coarse-space

   W \;=\; \operatorname{span}\{\mathbf{1}_R\}_R,
   \qquad
   \mathbf{1}_R(x) =
   \begin{cases} 1 & x \in R, \\ 0 & \text{otherwise,} \end{cases}

.. (vv-status rationale) Structural/representational identity: names the
   coarse trial space as the span of the coarse-cell indicators (the P0
   space). The implementing object is
   :class:`~orpheus.numerics.basis.IndicatorBasis`; the verifiable content
   is the membership-table / Gram bit-identity gated by
   ``tests.numerics.test_indicator_basis``, not a solver claim.
.. vv-status: sn-homogenization-coarse-space documented

the span of the **coarse-cell indicators** (the piecewise-constant /
P0 / box space; Brenner & Scott 2008 §3.4). A Galerkin projection would
*test* the residual against these same trial functions. Homogenization
does not: rate preservation forces the residual to be tested against the
**flux-weighted** indicators

.. math::
   :label: sn-homogenization-test-functions

   \chi_R(x) \;=\; \phi(x)\,\mathbf{1}_R(x),

.. (vv-status rationale) Structural/representational identity: names the
   Petrov-Galerkin TEST functions as the flux-weighted coarse-cell
   indicators χ_R = φ·1_R. The implementing object is
   :class:`~orpheus.numerics.basis.WeightedIndicatorBasis`; the verifiable
   content is the weighted-analysis bit-identity gated by
   ``tests.numerics.test_weighted_indicator_basis`` and the Mode-11
   routing sentinel, not a solver claim.
.. vv-status: sn-homogenization-test-functions documented

— a genuinely different basis from the trial :math:`\mathbf{1}_R`. With
test :math:`\ne` trial, the projection is **Petrov-Galerkin**: the
coarse coefficients are the solution of the Petrov-Galerkin normal
equations (test the residual against every test function, in the bare
geometric metric :math:`\mu_V` = weight :math:`V`),

.. math::
   :label: sn-homogenization-normal-equations

   \big\langle \chi_R,\; \Sigma - \Sigma_W \big\rangle_{V}
   \;=\; 0
   \quad \forall R
   \;\;\Longleftrightarrow\;\;
   c_R \;=\;
   \frac{\langle \chi_R,\, \Sigma \rangle_{V}}
        {\langle \chi_R,\, \mathbf{1}_R \rangle_{V}}
   \;=\;
   \frac{\sum_{i\in R} V_i\,\phi_{i,g}\,\Sigma_{i,g}}
        {\sum_{i\in R} V_i\,\phi_{i,g}},

.. (vv-status rationale) Derivation-decomposition step: the
   Petrov-Galerkin normal equations for the P0 projection with a
   flux-weighted test basis, whose solution IS the flux·volume collapse
   (:eq:`sn-homogenization-vector-collapse`). The verifiable content is
   the L0 rate-preservation gate plus the φV-vs-dV discriminator
   (``test_homogenization_is_flux_weighted_not_volume_weighted``); this
   is the projection-theoretic reading, not a separate claim.
.. vv-status: sn-homogenization-normal-equations documented

where :math:`\langle \chi_R, f\rangle_V = \int_R \phi\,f\,\mathrm{d}V`
is the flux-weighted region integral the test functions induce. Because
the indicators (trial *and* test) have **disjoint support**, the
**cross-Gram**

.. math::
   :label: sn-homogenization-cross-gram

   G_{RS} \;=\; \langle \chi_R,\, \mathbf{1}_S \rangle_{V}
   \;=\; \delta_{RS}\,\sum_{i\in R} V_i\,\phi_{i,g}
   \;=\; \delta_{RS}\,\Phi_{R,g}

.. (vv-status rationale) Structural identity: the cross-Gram of the
   homogenisation Petrov-Galerkin frame is diagonal (disjoint indicator
   supports), its diagonal being the region flux integral
   :eq:`sn-homogenization-region-flux`. The diagonality is exercised by the
   diagonal-Gram fast path in the L0 rate-preservation gate
   (``test_homogenization_is_flux_weighted_not_volume_weighted`` + the Mode-11
   ``apply_inverse_metric`` routing sentinel). A derivation-decomposition
   structural identity, not a separate solver claim.
.. vv-status: sn-homogenization-cross-gram documented

is **diagonal**, so the normal equations decouple cell-by-cell and each
coefficient is exactly the flux·volume average
:eq:`sn-homogenization-vector-collapse`. The denominator is the region
mass :math:`m_R = G_{RR} = \Phi_{R,g}` — the region flux integral
:eq:`sn-homogenization-region-flux` *is* the diagonal of the
cross-Gram. (Contrast the spherical-harmonic
:class:`~orpheus.numerics.frame.GalerkinFrame`, whose Gram is the
*symmetric* :math:`\langle Y_k, Y_j\rangle = \delta_{kj}/(2\ell+1)`
because there test :math:`=` trial; here the two factors of the Gram are
*different* bases, but disjoint support still diagonalizes it.)

**The test weighting is derived, not chosen.** Had the residual been
tested against the *plain* indicators :math:`\mathbf{1}_R` (the Galerkin
choice, test :math:`=` trial in the bare :math:`\mu_V` metric) the
projection would have been the **volume average**
:math:`\sum_i V_i \Sigma_i / \sum_i V_i`, which does *not* preserve the
reaction rate across a flux gradient. Matching rate preservation
:eq:`sn-homogenization-rate-preservation` is what *forces* the test
functions to be flux-weighted (:math:`\chi_R = \phi\,\mathbf{1}_R`)
rather than the plain :math:`\mathbf{1}_R`. This is the load-bearing
discriminator the regression gate
``test_homogenization_is_flux_weighted_not_volume_weighted`` pins: a
coarse region spanning a vacuum→reflective flux tilt over two materials
makes the flux-weighted and volume-only effective :math:`\Sigma_t`
numerically distinct, and production *must* match the flux-weighted one.

.. _sn-homogenization-why-petrov-galerkin:

Why Petrov-Galerkin and not Galerkin
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The flux-weighted projection *can* be written as a Galerkin projection
in a flux-weighted metric — fold :math:`\phi` from the test function
into the measure,

.. math::
   :label: sn-homogenization-metric-fold

   \big\langle \phi\,\mathbf{1}_R,\; \Sigma \big\rangle_{V}
   \;=\;
   \big\langle \mathbf{1}_R,\; \Sigma \big\rangle_{\phi V},

.. (vv-status rationale) Structural identity: the metric-fold that
   re-expresses the forward (φ*=φ) Petrov-Galerkin projection as a
   Galerkin projection in the L²(φV) metric. It is exact for the
   forward-flux case ONLY and is the convenience the #268 ruling rejects
   as the general framing (it folds the solution into the metric).
   Reframes the operator type; the numerical content is the same
   rate-preservation gate, not a new solver claim.
.. vv-status: sn-homogenization-metric-fold documented

making test and trial the *same* span :math:`\{\mathbf{1}_R\}` in the
:math:`L^2(\phi V)` metric. **This is a forward-only convenience, not
the structure.** It works here only because the test weight equals the
trial-side solution (:math:`\phi^* = \phi`) — the *forward* degenerate.
The homogenization reactor physics actually requires is
**eigenvalue-consistent**: the effective cross sections must keep the
multiplication factor :math:`\keff` stationary, and by first-order
perturbation theory :math:`\keff` is stationary with respect to the
**adjoint-weighted** residual. The functional that must be preserved is
therefore the **bilinear** form

.. math::
   :label: sn-homogenization-bilinear

   \big\langle \varphi^*,\, \Sigma\,\varphi \big\rangle,
   \qquad
   \Sigma_R \;=\;
   \frac{\int_R \varphi^*\,\Sigma\,\varphi\;\mathrm{d}V}
        {\int_R \varphi^*\,\varphi\;\mathrm{d}V},

.. (vv-status rationale) Literature-transcribed definition: the
   eigenvalue-consistent (adjoint-weighted) effective cross section, the
   bilinear ⟨φ*, Σφ⟩ that keeps k_eff stationary by first-order
   perturbation theory (Hébert §6; the generalized-equivalence /
   adjoint-weighting principle). The non-degenerate (φ*≠φ)
   Petrov-Galerkin case; implementation deferred to P6, blocked on the
   adjoint flux φ* (the #276 adjoint-transport campaign). Not yet an SN
   solver claim on this branch — documents the structure the forward
   case degenerates from.
.. vv-status: sn-homogenization-bilinear documented

with **test** functions :math:`\varphi^*\cdot\mathbf{1}_R` and
**trial** functions :math:`\varphi\cdot\mathbf{1}_R` that are now
genuinely distinct (:math:`\varphi^* \ne \varphi` away from a
self-adjoint problem). There is **no metric in which test equals
trial** — the map is irreducibly Petrov-Galerkin, :math:`M^* \ne R`.
The forward homogenization this slice ships is the **Galerkin
degenerate** of that map (:math:`\varphi^* = \varphi`, the flux is its
own adjoint weighting): it is a *legal* Galerkin reading because of the
coincidence :math:`\varphi^* = \varphi`, but the coincidence is *not*
the structure, so the honest framing — the one that survives the lift to
:math:`\varphi^* \ne \varphi` — is Petrov-Galerkin. ORPHEUS therefore
builds it as a :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` with
an explicit flux-weighted test basis, *not* a
:class:`~orpheus.numerics.frame.GalerkinFrame` with a flux-weighted
measure. The adjoint-weighted (:math:`\varphi^* \ne \varphi`) case
:eq:`sn-homogenization-bilinear` is **documented theory only** — its
implementation (phase P6 of the unified Frame-projection campaign) is
deferred, **blocked** on a converged adjoint flux :math:`\varphi^*` (the
#276 adjoint-transport campaign; see the capstone seam
:ref:`frame-adjoint-weighted-seam`). This section sets it up as the
non-degenerate sibling the forward case descends from.

.. note::

   The distinction is invisible numerically on *this* branch — forward
   homogenization with :math:`\varphi^* = \varphi` produces the same
   numbers either way (the metric-fold :eq:`sn-homogenization-metric-fold`
   is an exact identity for :math:`\varphi^* = \varphi`). It is an
   **architecture** distinction: writing the discipline on the *type*
   (an explicit test basis) rather than on the *measure* (a flux-folded
   metric) is what lets the adjoint-weighted case land as a no-op change
   of the test basis (:math:`\varphi \to \varphi^*`) rather than a
   re-derivation. The Mode-11 routing sentinel below pins that the flux
   genuinely lives on the test side, so a regression that silently
   re-folded it into the metric would be caught even though it would
   not move a single number.

The measure carries the axis, never the discipline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The implementation key is that the flux weight rides the **test basis**,
and the :class:`~orpheus.numerics.measure.DiscreteMeasure` carries only
the bare geometric volume :math:`\mu_V`. The frame factors the test
functional into a geometric measure and a solution multiplier,

.. math::
   :label: sn-homogenization-radon-nikodym

   \langle \chi_R, f\rangle_V \;=\; \int_R \phi\,f\;\mathrm{d}\mu_V,
   \qquad
   \chi_R \;=\; \phi\cdot\mathbf{1}_R,

.. (vv-status rationale) Structural identity: the test functional splits
   into the geometric base measure μ_V and the flux density φ carried on
   the test basis. It is the design rationale for carrying φ as the test
   weight (an integrand multiplier) rather than as ng separate measures
   or as a metric on the measure; the verifiable content is the measure /
   weighted-basis bit-identity gates and the Mode-11 routing sentinel,
   not a solver claim.
.. vv-status: sn-homogenization-radon-nikodym documented

i.e. the **geometric base measure** :math:`\mu_V` (group-independent,
the coarse/fine mesh's :attr:`volume_measure
<orpheus.transport.mesh.material_mesh.MaterialMesh.volume_measure>` — a
:class:`~orpheus.numerics.measure.DiscreteMeasure`) is multiplied at
integration time by the flux :math:`\phi` *carried on the test basis*.
The code carries exactly this split: the
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame` binds the trial
:class:`~orpheus.numerics.basis.IndicatorBasis` and the *bare*
group-independent :math:`\mu_V`, and the flux enters through the
**test basis**
:class:`~orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis`,
whose :meth:`analyze
<orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis.analyze>`
folds the per-group flux into the integrand on a trailing tensor axis —
``test.analyze(phi * channel_fine, …)`` — so the whole group structure
rides one frame.

This is *why* the test weight is **not** smuggled onto the measure.
:class:`~orpheus.numerics.measure.DiscreteMeasure`'s ``weights`` array
stays **1-D** (one mass per atom) and group-independent; a per-group
:math:`\mu_{\phi V}` would be :math:`n_g` distinct measures, and — worse
— a *measure*-borne flux weight is exactly the metric-fold the #268
ruling forbids as the general framing: it works for forward homogenization
and breaks under :math:`\varphi^* \ne \varphi`. Keeping :math:`\phi` on
the test basis instead of the measure forces the correct reading:
:math:`\phi` is a test-weighting the *solution* emits, not a property of
the geometry. The geometry (the mesh) owns one measure :math:`\mu_V`; the
solution owns the flux :math:`\phi`; the *frame type* (Petrov-Galerkin,
with its explicit test basis) carries the discipline. The
:class:`~orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis`
is **test-only** by construction — its :meth:`evaluate
<orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis.evaluate>`
is the *weight-free* geometric membership (the weight is an *analysis*
weight, not a tabulation), and its synthesis-side operations *raise*
(the Petrov-Galerkin reconstruction is purely trial-side; building a
weighted synthesis before a consumer exists would make a half-consistent
basis).

The mesh yields the basis; it does not inherit it
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The coarse trial space :eq:`sn-homogenization-coarse-space` is realised
by a **new** concrete numerics basis,
:class:`~orpheus.numerics.basis.IndicatorBasis` — the second concrete
:class:`~orpheus.numerics.basis.Basis` after
:class:`~orpheus.numerics.basis.SphericalHarmonicBasis`, and the
piecewise-constant (P0 / characteristic-function) analogue of it. The
coarse :class:`~orpheus.geometry.mesh.Mesh1D` **yields** this view via
:meth:`coarse.indicator_basis() <orpheus.geometry.mesh.Mesh1D.indicator_basis>`,
exactly symmetric with how it already yields
:meth:`coarse.volume_measure <orpheus.geometry.mesh.Mesh1D.volume_measure>`.

The mesh is **not** a :class:`~orpheus.numerics.basis.Basis` subclass,
and the reason is a clean role separation: a
:class:`~orpheus.numerics.basis.Basis` is the **measure-free** half of a
frame, while a mesh **carries the volume measure**. A mesh that *were* a
basis would conflate the two roles of the frame pair (the
discipline-free trial side and the measured test side) into one object.
So the mesh yields *both* views — its measure-free indicator basis and
its volume measure — and the
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame` binds them together
with the flux-weighted test basis. The yielded
:class:`~orpheus.numerics.basis.IndicatorBasis` is **geometry-free**: it
holds only the per-axis edge arrays, so :mod:`orpheus.numerics` carries
no dependency on :mod:`orpheus.geometry`. Its
:meth:`evaluate <orpheus.numerics.basis.IndicatorBasis.evaluate>` builds
the :math:`(n_{\rm fine} \times n_{\rm coarse})` one-hot **membership
table** :math:`T[i,R] = \mathbf{1}_R(x_i)` by a per-axis
``searchsorted`` followed by :func:`numpy.ravel_multi_index` in ``"ij"``
order — the *same* flat-cell ordering the volume measure uses for its
nodes, so the table column index and the measure node index agree by
construction in any dimension (no 1-D special case in the membership
machinery). This is what makes homogenization **dimension-agnostic**:
:meth:`Solution.homogenize <orpheus.sn.solution.Solution.homogenize>`
flattens its ``(ng, *spatial)`` flux to ``(n_fine, ng)`` in the same
``"ij"`` order and a 1-D or 2-D mesh flows through the one frame body
(pinned end-to-end by ``test_homogenize_2d_rate_preservation``).

The coefficient-extraction verb and its normalisation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The homogenization map is the frame's **coefficient-extraction verb**
:meth:`frame.project <orpheus.numerics.frame.FrameBase.project>`
= :math:`G^{-1} M`:

.. math::
   :label: sn-homogenization-frame-projector

   \Sigma_R \;=\; \big(G^{-1} \circ M\big)\,\Sigma,
   \qquad
   M = \text{analysis} \;\big(\textstyle\int_R \phi\,\cdot\;\mathrm{d}V\big),\;\;
   G^{-1} = \operatorname{diag}(1/\Phi_R),

.. (vv-status rationale) Structural/representational identity: the
   homogenization map as the inverse-Gram ∘ analysis coefficient
   extraction (``FrameBase.project``) of the Petrov-Galerkin frame.
   Each factor is a Frame primitive whose bit-identity is gated by
   ``tests.numerics.test_indicator_basis`` /
   ``tests.numerics.test_weighted_indicator_basis`` /
   ``tests.numerics.test_frame``; the Mode-11 sentinel
   ``test_homogenize_routes_through_the_petrov_galerkin_frame`` pins that
   ``homogenize`` actually calls them. Not a separate solver claim.
.. vv-status: sn-homogenization-frame-projector documented

read right-to-left:

#. **analysis** :math:`M` = ``frame.analysis.apply``, which delegates to
   the *test* basis's weighted analysis: the flux-weighted region
   integral :math:`(M\Sigma)_R = \sum_{i\in R} V_i\,\phi_i\,\Sigma_i =
   \int_R \phi\,\Sigma\,\mathrm{d}V` — the region reaction rate. The
   diagonal of the cross-Gram is recovered by the *same* face applied to
   the constant field, :math:`(M\,\mathbf 1)_R = \sum_{i\in R}
   V_i\,\phi_i = \Phi_R` (a single ``analysis ∘ reconstruction`` probe
   of the all-ones coefficient vector; the off-diagonals are
   structurally zero, so the row-sum IS the diagonal — see
   :meth:`frame.gram <orpheus.numerics.frame.FrameBase.gram>`).
#. **inverse Gram** :math:`G^{-1} = \operatorname{diag}(1/\Phi_R)` =
   :meth:`FunctionSpace.apply_inverse_metric
   <orpheus.numerics.space.FunctionSpace.apply_inverse_metric>` on a
   coarse coefficient space whose installed metric is the diagonal Gram
   :math:`\Phi_R = M\,\mathbf 1`. The normalisation :math:`1/\Phi_R` is
   **measure-dependent** (the region mass changes with the flux weight),
   *unlike* the spherical-harmonic :math:`2\ell+1` factor which is
   analytic and measure-free. A measure-dependent factor **cannot** live
   on the measure-free :class:`~orpheus.numerics.basis.Basis`, so the
   trial :meth:`reconstruct
   <orpheus.numerics.basis.IndicatorBasis.reconstruct>` stays the plain
   (identity-dual) broadcast and the :math:`1/\Phi_R` normalisation is
   applied **separately** by the coefficient space's metric. The metric
   is a **Moore–Penrose pseudo-inverse**: a coarse cell with zero region
   flux (:math:`\Phi_R = 0`) is in the metric's null space and gets
   effective :math:`\Sigma_R = 0` — the :math:`0/0` branch of
   :eq:`sn-homogenization-vector-collapse` resolved for free, with no
   special-casing in :meth:`Solution.homogenize`.

For the **matrix channels** the same verb runs with the source-group
flux as the test weight: :math:`\phi_{g'}` rides the ``g_from`` axis (the
leading axis the test weight aligns to — see
:meth:`WeightedIndicatorBasis._weighted
<orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis>`'s
leading-aligned broadcast), and ``apply_inverse_metric`` broadcasts the
per-region Gram :math:`\Phi_R[:, g_{\rm from}]` over the trailing
``g_to`` axis — so the source-group normalisation of
:eq:`sn-homogenization-matrix-collapse` falls out of the
**trailing-axis metric-broadcast** rather than needing its own code
path. The :math:`\chi` channel uses the *identical* machinery in a
*separate* frame with a *different* test weight — the per-cell production
density :math:`p_i = \sum_g \nu\Sigma_{f,i,g}\,\phi_{i,g}\,V_i`
(:eq:`sn-homogenization-production-weight`) — so its Gram becomes the
region production :math:`\sum_{i\in R} V_i\,p_i` and the projection is the
production-weighted convex average :eq:`sn-homogenization-chi-collapse`.

Two test weightings, two frames — one conserved rate each
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The vector and matrix :math:`\Sigma` channels and the emission spectrum
:math:`\chi` do **not** share a frame, because they preserve **two
different conserved rates**, and a Petrov-Galerkin frame carries exactly
one test weighting:

.. list-table::
   :header-rows: 1
   :widths: 22 30 26 22

   * - Frame
     - Channels
     - Conserved rate
     - Test weight :math:`w`
   * - ``sigma_frame``
     - :math:`\Sigma_t,\Sigma_c,\Sigma_L,\Sigma_f,\nu\Sigma_f` (vectors);
       :math:`\Sigma_{s,\ell},\Sigma_{2n}` (``[g_from, g_to]`` matrices)
     - **reaction rate** :math:`\sum_i V_i\phi_{i,g}\Sigma_{i,g}`
     - per-group flux :math:`\varphi` (matrices weight by the
       **source** group)
   * - ``emission_frame``
     - :math:`\chi` (emission spectrum)
     - **emission rate** :math:`\sum_i p_i\chi_{i,g}`
     - production density :math:`p=\sum_g\nu\Sigma_{f,g}\varphi_g`

Both frames bind the **same** trial
:class:`~orpheus.numerics.basis.IndicatorBasis` and the **same**
geometric measure :math:`\mu_V`; they differ *only* in the test basis
they carry. :meth:`Solution.homogenize` builds both and hands them to
:meth:`MaterialXSField.project_through
<orpheus.transport.mesh.material_xs_field.MaterialXSField.project_through>`,
which owns the **channel → frame** taxonomy: it routes every rate-bearing
:math:`\Sigma` channel through ``sigma_frame`` and :math:`\chi` through
``emission_frame``, projecting the *whole* cross-section field as one
object and returning one effective
:class:`~orpheus.data.macro_xs.mixture.Mixture` per coarse cell. The
caller owns the flux, so the caller builds the test weightings; the field
owns *which* weighting each channel needs.

.. note::

   That homogenization actually *executes* through these Frame readers —
   :meth:`IndicatorBasis.evaluate
   <orpheus.numerics.basis.IndicatorBasis.evaluate>` (the trial
   membership), :meth:`WeightedIndicatorBasis.analyze
   <orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis.analyze>`
   (the **test-side** flux-weighted reader),
   :meth:`FrameBase.project <orpheus.numerics.frame.FrameBase.project>`
   (the :math:`G^{-1}M` verb), and
   :meth:`FunctionSpace.apply_inverse_metric
   <orpheus.numerics.space.FunctionSpace.apply_inverse_metric>` — rather
   than a green rate gate riding a buggy refactor that recomputes
   membership inline or quietly re-folds :math:`\phi` into the metric, is
   pinned by the **Mode-11 sentinel**
   ``test_homogenize_routes_through_the_petrov_galerkin_frame``
   (vv-principles **Mode 11**, gate-never-executes-the-rewired-path): it
   monkeypatch-counts each reader and asserts the counter is positive
   after a ``homogenize`` run. The load-bearing count is
   ``WeightedIndicatorBasis.analyze`` — a bit-identity-preserving
   regression that kept the *old* Galerkin metric-fold (folding
   :math:`\phi` into the coefficient-space metric, test = plain trial
   indicator) would produce **identical numbers** yet never construct the
   weighted test basis, leaving that counter at zero. The
   rate-preservation identity :eq:`sn-homogenization-rate-preservation`
   remains THE correctness claim (the L0 gate); the sentinel makes the
   *Petrov-Galerkin structure* load-bearing for *this* implementation.

Why route through the Frame at all
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The membership-matmul the prior slice shipped was correct; the carve to
the :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` is an
**architecture** move, not a correctness fix (Cardinal Rule 2). Three
payoffs justify it:

* **One mechanism, not one per method.** The angular-flux →
  spherical-harmonic-moment projection of SN anisotropic scattering (a
  :class:`~orpheus.numerics.frame.GalerkinFrame`) and the fine → coarse
  cross-section projection of homogenization (a
  :class:`~orpheus.numerics.frame.PetrovGalerkinFrame`) are the *same*
  mechanism — a discrete frame's analysis/reconstruction pair —
  differing only in *which* :class:`~orpheus.numerics.basis.Basis` pair
  (trial / test) is bound and *which discipline type* carries them.
  Routing both through the :class:`~orpheus.numerics.frame.FrameBase`
  hierarchy collapses a twin path (coding-elegance anti-pattern 1) into
  one body.
* **Energy condensation becomes the same shape.** The deferred
  ``.condense`` sibling is the identical Petrov-Galerkin frame with the
  spatial :class:`~orpheus.numerics.basis.IndicatorBasis` replaced by a
  *spectral* group-indicator basis and the measure replaced by
  :math:`L^2(\text{spectrum})`. Establishing the frame routing here means
  condensation lands as a no-op extension through the same body, not a
  third arm.
* **The pseudo-inverse handles the empty-region edge case for free.** The
  :math:`0/0` of a flux-free coarse cell is the metric's null space, so
  it needs no guard in :meth:`Solution.homogenize` — the projection
  algebra absorbs it.

.. _sn-condense-homogenize-asymmetry:

The condense / homogenize asymmetry law
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Homogenization (space) and condensation (energy) are siblings, but they
are **not** symmetric — they return different types, and the asymmetry
is structural, not incidental:

.. list-table::
   :header-rows: 1
   :widths: 18 28 28 26

   * - Operation
     - Collapses
     - Mesh coupling
     - Return type
   * - **homogenize**
     - space (fine cells → coarse cells)
     - **mesh-COUPLED** — the effective materials *are* the coarse
       cells; geometry and materials are born together
     - :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`
       (geometry + materials)
   * - **condense**
     - energy (fine groups → coarse groups)
     - **mesh-DECOUPLED** — the condensed cross sections are
       group-structure data, independent of where they sit
     - ``dict[int, Mixture]`` (portable materials, *no* geometry)

Spatial homogenization is **mesh-coupled**: a homogenized material has
no meaning apart from *the coarse cell it homogenizes*, because the
flux·volume weights :math:`w_{i,g}` are tied to specific fine cells
inside a specific coarse region. The natural product is therefore the
coarse geometry *carrying* its materials — a
:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`, the
mesh+materials data carrier minted for exactly this purpose. (This is
why ``MaterialMesh`` exists as the middle type between a bare
:class:`~orpheus.geometry.mesh.Mesh1D` and a method phase space: a
homogenized model is *materials-and-geometry-together but not yet
method-specific* — it has no quadrature until
:meth:`SNMesh.from_material_mesh
<orpheus.sn.mesh.augmented_mesh.SNMesh.from_material_mesh>` promotes it.)

Energy condensation is **mesh-decoupled**: a condensed cross-section set
is just a coarser :class:`Mixture` — group-structure data that can be
attached to *any* geometry. Its natural product is a portable
``dict[int, Mixture]`` keyed by material id, with no geometry attached.

The asymmetry has a clean **frame-theoretic** reading once
homogenization is seen as the Petrov-Galerkin frame projection of
:ref:`sn-homogenization-petrov-galerkin-frame`. Both operations are the
*same* frame mechanism — analysis ∘ inverse-Gram — and differ *only* in
their trial basis :math:`K` (the
:class:`~orpheus.numerics.basis.Basis` bound into the frame):

* **homogenize** binds a **geometric** :math:`K` — the spatial
  :class:`~orpheus.numerics.basis.IndicatorBasis` (cell indicators
  :math:`\{\mathbf{1}_R\}`). Its coefficients *are* the coarse cells, so
  the result is **mesh-coupled** and the natural product is geometry +
  materials (:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`).
* **condense** binds a **spectral** :math:`K` — a group-indicator basis
  on the energy axis (broad-group indicators :math:`\{\mathbf{1}_G\}`)
  under the :math:`L^2(\text{spectrum})` measure. Its coefficients are
  *group-structure data*, carrying no spatial identity, so the result is
  **mesh-decoupled** and the natural product is a portable
  ``dict[int, Mixture]``.

The return-type asymmetry is therefore not incidental: it is the
:math:`K`-axis (space vs energy) of the *one* projection mechanism
showing through. A geometric trial basis births geometry; a spectral
trial basis births portable group data.

The condensation half is realised by :meth:`Solution.condense
<orpheus.sn.solution.Solution.condense>`; its theory — the per-material
rate-preserving collapse, the fractional-overlap re-binning that lifts
the nesting requirement, and the same Petrov-Galerkin discipline — is
:ref:`sn-energy-condensation` below, the energy-axis transpose of this
section.

.. _sn-homogenization-verification:

Verification
~~~~~~~~~~~~

The gate is :mod:`tests.sn.test_homogenization` (level **L0**, term
verification — it checks the defining identity term-by-term, not a
solver claim). Its load-bearing test asserts the rate-preservation
identity :eq:`sn-homogenization-rate-preservation` directly:

* **Vector channels** (``test_rate_preservation_vector_channels``) —
  for every channel
  (:math:`\Sigma_t,\ \Sigma_c,\ \Sigma_L,\ \Sigma_f,\ \nu\Sigma_f`),
  every coarse region, and every group, assert
  :math:`\Sigma_{R,g}\,\Phi_{R,g} = \sum_{i\in R} V_i\,\Sigma_{i,g}\,
  \phi_{i,g}` to machine precision.
* **Matrix channels** (``test_rate_preservation_scattering_and_n2n``) —
  the same identity for every Legendre order of
  :math:`\Sigma_{s,\ell}` and for :math:`\Sigma_{2n}`, with the
  reference loop weighting by the **source** group ``g_from``
  *explicitly* — which is what makes it catch a ``g_from``↔``g_to``
  swap (vv-principles **Mode 2**).
* **n-D** (``test_homogenize_2d_rate_preservation``) — the same
  rate-preservation identity end-to-end through a *real* 2-D
  ``solve_sn``, exercising the n-D membership
  (:func:`numpy.ravel_multi_index` ``"ij"``) and the
  ``(ng, nx, ny) → (n_fine, ng)`` flatten the dropped 1-D guard opens, a
  flux tilt keeping :math:`\phi` non-flat so the flux-weighting is
  genuinely activated.

The reference these are checked against is a **structurally-independent**
explicit per-region Python loop over the fine cells — *not* a re-call of
the production frame projection (vv-principles **L11**: a cross-check
must be structurally independent, not merely procedurally independent;
a frame-vs-loop comparison sharing the *same* region reduction would
share any bug in the reduction). Companion invariants pin the rest of
the contract:

.. list-table::
   :header-rows: 1
   :widths: 46 54

   * - Test
     - What it pins
   * - ``test_homogenized_materials_balance``
     - Balance :eq:`sn-homogenization-balance` survives the collapse —
       every removal channel shares the per-:math:`(R,g)` weight.
   * - ``test_homogenized_chi_is_simplex``
     - :math:`\chi_R` is a probability simplex (convex average of
       producing simplices, :eq:`sn-homogenization-chi-collapse`).
   * - ``test_chi_is_production_weighted``
     - :math:`\chi_R` uses the **production** weight
       :eq:`sn-homogenization-production-weight`, not a flux- or
       volume-weight — the simplex test is blind to *which* convex weight,
       so this pins the weight choice directly.
   * - ``test_homogenization_is_flux_weighted_not_volume_weighted``
     - The load-bearing **flux-weighted-test** guard: over a
       vacuum→reflective flux tilt the flux-weighted and volume-only
       effective :math:`\Sigma_t` are numerically distinct, and
       production MUST match the flux-weighted one
       (:eq:`sn-homogenization-normal-equations`) — reds a regression
       that drops :math:`\phi` from the test weight (Galerkin /
       volume-only averaging).
   * - ``test_homogenize_routes_through_the_petrov_galerkin_frame``
     - **Mode-11 sentinel**: ``homogenize`` actually calls
       :meth:`IndicatorBasis.evaluate
       <orpheus.numerics.basis.IndicatorBasis.evaluate>` (trial
       membership), :meth:`WeightedIndicatorBasis.analyze
       <orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis.analyze>`
       (the **test-side** flux-weighted reader),
       :meth:`FrameBase.project
       <orpheus.numerics.frame.FrameBase.project>`, and
       :meth:`FunctionSpace.apply_inverse_metric
       <orpheus.numerics.space.FunctionSpace.apply_inverse_metric>` — the
       Petrov-Galerkin routing
       (:eq:`sn-homogenization-frame-projector`) is on the gate's call
       graph, not bypassed by an inline recompute *or a silent re-fold of*
       :math:`\phi` *into the metric* (which would keep the numbers and
       leave ``WeightedIndicatorBasis.analyze`` at zero calls).
   * - ``test_effective_xs_bracketed_by_fine_extremes``
     - :math:`\Sigma_{t,R}` is bracketed by the region's fine-cell
       extremes — a physical sanity check independent of the rate gate.
   * - ``test_identity_homogenization_recovers_per_cell_materials``
     - Homogenizing onto the *same* fine mesh recovers each cell's
       original material (degenerate limit: one fine cell per coarse).
   * - ``test_single_material_region_recovers_material``
     - A coarse cell containing only material :math:`m` homogenizes
       back to :math:`m` (the flux weight cancels).
   * - ``test_outer_boundary_mismatch_raises``
     - The guard: a coarse mesh whose outer boundary differs from the
       fine mesh raises ``ValueError``.

The :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` data
contract itself — the ``ng``-consistency check, the volume measure, the
XS-field build, and the ``SNMesh(MaterialMesh)`` data/behavior split
(including a bit-identity check that ``SNMesh``'s inherited data block
matches a standalone ``MaterialMesh``) — is gated separately by
:mod:`tests.transport.test_material_mesh`.


.. _sn-energy-condensation:

Applied to energy condensation
------------------------------

Spatial homogenization (:ref:`sn-spatial-homogenization`) collapses the
*space* axis; **energy condensation** is its **energy-axis transpose** —
it collapses a fine-group cross-section set onto a coarser group
structure, **spectrum-weighted**, so that each coarse group reproduces
the fine reaction rate. The two are the classical pair of "smear the
detail you have resolved into effective constants for a coarser
calculation" moves (Hébert, *Applied Reactor Physics* :cite:`Hebert2009`,
§13 for space, §3.5 for energy).

In ORPHEUS condensation lives at two layers, mirroring how
homogenization splits between :meth:`Solution.homogenize` (the
orchestration) and :meth:`MaterialXSField.project_through
<orpheus.transport.mesh.material_xs_field.MaterialXSField.project_through>`
(the per-channel collapse):

* :meth:`Mixture.condense
  <orpheus.data.macro_xs.mixture.Mixture.condense>` — the per-material
  channel collapse (the data layer). Given a coarse target
  :class:`~orpheus.data.energy_grid.EnergyGrid` and a representative
  spectrum, ``mix.condense(target, spectrum)`` builds the fine→coarse
  fractional-overlap trial internally
  (:meth:`mix.energy_grid.overlap_to(target)
  <orpheus.data.energy_grid.EnergyGrid.overlap_to>`) and returns the
  condensed (coarse-group)
  :class:`~orpheus.data.macro_xs.mixture.Mixture`. It is **data-native** —
  every object it touches (the grid, the overlap factory, the
  Petrov-Galerkin frame) lives in ``data`` / ``numerics``, with **no**
  transport dependency.
* :meth:`Solution.condense <orpheus.sn.solution.Solution.condense>` —
  the orchestration (the SN layer). It derives each material's
  representative spectrum from the solved flux and returns a
  **portable** ``dict[int, Mixture]`` keyed by material id.

.. note::

   This is the **energy-only** slice (the energy sibling of spatial
   homogenization). Geometry is **not** touched — the result is portable
   few-group cross sections, not a mesh. The asymmetry between the two
   operations, and *why* they return different types
   (``dict[int, Mixture]`` vs
   :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`), is
   :ref:`sn-condense-homogenize-asymmetry` above. Energies obey the
   canonical fast-first convention (group ``0`` = fastest, descending
   boundaries; :ref:`canonical-group-convention`).

The defining property: reaction-rate preservation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Condensation, like homogenization, is defined by *what it must
preserve*, not by an averaging recipe. The quantity a transport
calculation consumes is the **reaction rate** in each group, and the
group-collapse cross section is *defined* so as to preserve it (Hébert
:cite:`Hebert2009` Eq. 3.103; Stamm'ler & Abbate, *Methods of Steady-State
Reactor Physics in Simplified Geometry* (1983), Eq. VI(6b) — two
independent authoritative textbooks state it identically). Fix a coarse
group :math:`G`, made up of fine groups :math:`g \in G`. The fine-group
reaction rate of any vector channel is

.. math::
   :label: energy-condensation-fine-rate

   r_G \;=\; \sum_{g \in G} \varphi_g\,\Sigma_g,

.. (vv-status rationale) Definitional identity: the fine-group reaction
   rate, a plain group sum because ORPHEUS φ_g is already
   group-integrated (the bin width is inside φ_g; see
   :eq:`energy-condensation-counting-measure`). A
   derivation-decomposition premise of the rate-preservation claim
   (:eq:`energy-condensation-rate-preservation`), which is the verifiable
   claim (L1 gate ``tests.data.test_mixture_condense``).
.. vv-status: energy-condensation-fine-rate documented

with :math:`\varphi_g` the per-material representative flux (the test
weight, fixed below) and :math:`\Sigma_g` the fine cross section for
whatever channel (total, capture, fission, …) is in question. The sum
has **no** :math:`\mathrm{d}E` or :term:`lethargy` width because ORPHEUS's
:math:`\varphi_g` is already the group-integrated flux
:math:`\int_g \phi\,\mathrm{d}E` — the bin width is baked into the flux
(:eq:`energy-condensation-counting-measure`). The coarse model carries
one effective cross section :math:`\Sigma_G` and one coarse-group flux
:math:`\Phi_G` per group; for it to reproduce the fine rate it must
satisfy

.. math::
   :label: energy-condensation-rate-preservation

   \Sigma_G\,\Phi_G
   \;=\;
   \sum_{g \in G} \varphi_g\,\Sigma_g,

This is the **only** constraint condensation imposes on the vector
channels; everything below is derived from it. It is the energy-axis
copy of :eq:`sn-homogenization-rate-preservation` (with the fine-cell
volume·flux weight :math:`V_i\phi_{i,g}` replaced by the fine-group flux
:math:`\varphi_g`, and the spatial region :math:`R` replaced by the
energy group :math:`G`). The coarse-group flux is fixed first, by the
production-free inventory match — the coarse flux is the fine flux
summed over the group,

.. math::
   :label: energy-condensation-coarse-flux

   \Phi_G \;=\; \sum_{g \in G} \varphi_g.

.. (vv-status rationale) Representational identity: the coarse-group flux
   is the fine flux summed over the group (the diagonal Gram Φ_G of the
   PG frame). A definition consumed by the rate-preservation claim
   (:eq:`energy-condensation-rate-preservation`), not an independent
   solver claim.
.. vv-status: energy-condensation-coarse-flux documented

Substituting :eq:`energy-condensation-coarse-flux` into
:eq:`energy-condensation-rate-preservation` and solving for the
effective cross section gives the **spectrum-weighted average**

.. math::
   :label: energy-condensation-vector-collapse

   \Sigma_G
   \;=\;
   \frac{\sum_{g \in G} \varphi_g\,\Sigma_g}
        {\sum_{g \in G} \varphi_g},

.. (vv-status rationale) Derivation-decomposition step: the
   spectrum-weighted collapse obtained by solving the rate-preservation
   identity (:eq:`energy-condensation-rate-preservation`) for the
   effective cross section — Hébert Eq. 3.103. The verifiable content is
   the L1 rate-preservation gate; this is the algebraic rearrangement,
   not a separate claim.
.. vv-status: energy-condensation-vector-collapse documented

the flux-weighted reaction-rate-preserving average (Hébert
:cite:`Hebert2009` Eq. 3.103 ≡ Stamm'ler VI(6b)). Because the weight
:math:`\varphi_g` appears in both the numerator (rate) and denominator
(flux), :math:`\Sigma_G` is a genuine convex combination of the fine
values: it is bracketed by the group's fine extremes
:math:`\min_{g\in G}\Sigma_g \le \Sigma_G \le \max_{g\in G}\Sigma_g`.
This is *not* a separate design choice — it falls straight out of
:eq:`energy-condensation-rate-preservation`. Choosing any other weight
(width-only, unweighted) would break rate preservation wherever the
spectrum varies across the coarse group, which is exactly the regime
condensation exists to handle. The vector channels — total, capture,
leakage-loss (:math:`\Sigma_L`), fission, and production
(:math:`\nu\Sigma_f`) — each collapse through
:eq:`energy-condensation-vector-collapse` with the *same* weight; a
coarse group with zero flux (:math:`\Phi_G = 0`) has no reaction rate to
preserve, so its effective cross section is zero — the :math:`0/0`
resolved by the only physically meaningful value (the frame's
Moore–Penrose Gram, below).

The counting measure: why the weight is :math:`\varphi_g`, not :math:`\Delta u_g\,\varphi_g`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A subtle but load-bearing point distinguishes the energy axis from the
spatial axis. In Hébert's continuous formulation the flux-weighted
average of a *distribution* (a reaction rate) is a plain lethargy
integral (Eq. 3.96, :math:`\langle X\rangle_g = \int_g X\,\mathrm{d}u`)
while the average of a *function* (the flux) carries a :math:`1/\Delta
u_g` lethargy-width normalisation (Eq. 3.97). The discrete weight that
preserves the rate is therefore

.. math::
   :label: energy-condensation-counting-measure

   w_g \;=\; 1
   \qquad\text{(counting), not}\qquad
   w_g \;=\; \Delta u_g \;=\; \ln(E_{g}^{\rm upper}/E_{g}^{\rm lower}),

.. (vv-status rationale) Structural identity: the energy-axis measure is
   COUNTING (w=1), because ORPHEUS φ_g is already group-integrated, so
   the discrete rate is a plain group sum. Design rationale for the
   measure choice; the verifiable content is the rate-preservation gate
   (a Δu weight breaks it) plus the
   :class:`~orpheus.numerics.measure.DiscreteMeasure` bit-identity, not a
   solver claim.
.. vv-status: energy-condensation-counting-measure documented

because ORPHEUS stores :math:`\varphi_g = \int_g \phi\,\mathrm{d}E`
already integrated over the bin ("eV-free"). The discrete rate is then a
plain group **sum** :math:`r_G = \sum_{g\in G}\varphi_g\Sigma_g`
(:eq:`energy-condensation-fine-rate`) — the bin width is already inside
:math:`\varphi_g`. Verified against the frame algebra:
:math:`\Sigma_G\cdot\Phi_G = \sum_g w_g\,\varphi_g\Sigma_g` equals the
physical rate **iff** :math:`w_g = 1`; introducing a :math:`\Delta u_g`
weight would double-count the width and break rate preservation. This is
the energy-axis analogue of the spatial case, where :math:`\phi_i` *is*
a density and therefore *does* need the geometric volume measure
:math:`V_i` (:eq:`sn-homogenization-fine-rate`); here the measure is
:math:`w_g = 1` because the integration is already done. Lethargy is the
node *coordinate*, never the *weight* (it reappears below as the
within-group flux model, :eq:`energy-condensation-overlap-fraction`,
which sets a fine group's *split* across coarse groups — a basis datum —
not the measure). The
:class:`~orpheus.numerics.measure.DiscreteMeasure` the condensation
frame binds is therefore a **counting** measure
(:attr:`weights = ones <orpheus.numerics.measure.DiscreteMeasure>`,
``support="energy"``).

The matrix channels: a two-axis collapse (sink summed, source averaged)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The scattering matrices :math:`\Sigma_{s,\ell}[g',g]` (one per Legendre
order :math:`\ell`) and the :math:`(n,2n)` matrix
:math:`\Sigma_{2n}[g',g]` carry **two** group indices, stored
``[g_from, g_to]`` (the ORPHEUS scattering convention — see
:ref:`theory-cross-section-data`). They collapse by a **two-axis** rule
that has *no spatial precedent* and is the energy-condensation analogue
that most differs from homogenization. The in-scatter rate from coarse
group :math:`G` into coarse group :math:`G'` that must be preserved is

.. math::
   :label: energy-condensation-scattering-collapse

   \Phi_G\,\Sigma_{s,\ell,G\to G'}
   \;=\;
   \sum_{g \in G}\sum_{g' \in G'}
   \varphi_g\,\Sigma_{s,\ell}[g, g'],

.. (vv-status rationale) Definitional identity: the in-scatter rate of a
   matrix channel that the 2-axis collapse must preserve — the SINK axis
   g' summed (every fine scatter into any fine group of G' is a scatter
   into G'), the SOURCE axis g flux-averaged. Hébert Eq. 3.104,
   Stamm'ler VI(6c). The verifiable claim is the L1 scattering-collapse
   gate (which catches a g_from↔g_to swap, vv Mode 2).

driven by the population of the **source** group :math:`G` (the group
whose flux scatters *out* — the number of :math:`G\to G'` events scales
with how many particles are *in* :math:`G`). Decompose the collapse into
its two axes:

#. **Sink axis** :math:`g'` (the destination group) is **summed**: any
   scatter into *any* fine group :math:`g'` of coarse :math:`G'` is a
   scatter into :math:`G'`. The destination has no rate to average — it
   is an accumulation. In matrix form the sink-sum is a right-multiply
   by the membership table :math:`T` (below):
   :math:`\Sigma^{\rm sink}[g, G'] = \sum_{g'} \Sigma_{s,\ell}[g, g']\,T[g', G'] = (\Sigma_{s,\ell}\,T)[g, G']`.
#. **Source axis** :math:`g` (the origin group) is **flux-averaged**, by
   the *same* :eq:`energy-condensation-vector-collapse` rule applied to
   the sink-summed matrix:
   :math:`\Sigma_{s,\ell,G\to G'} = \big(\text{project}\,(\Sigma_{s,\ell}\,T)\big)[G, G']`.

So the matrix collapse is

.. math::
   :label: energy-condensation-matrix-collapse

   \Sigma_{s,\ell,G\to G'}
   \;=\;
   \frac{\sum_{g \in G} \varphi_g
         \big(\sum_{g'\in G'} \Sigma_{s,\ell}[g, g']\big)}
        {\sum_{g \in G} \varphi_g}
   \;=\;
   \operatorname{project}\!\big(\Sigma_{s,\ell}\,T\big),

.. (vv-status rationale) Derivation-decomposition step: the two-axis
   matrix collapse — sink summed (@T), source flux-averaged (project) —
   the matrix-channel form of :eq:`energy-condensation-scattering-collapse`.
   The verifiable content is the L1 scattering-collapse gate plus its
   three mutation probes (swap axes / sum both / project both); this is
   the algebraic form, not a separate claim.
.. vv-status: energy-condensation-matrix-collapse documented

with :math:`\operatorname{project}` the source-group flux average
(:meth:`frame.project <orpheus.numerics.frame.FrameBase.project>`) and
:math:`T` the fine→coarse membership table. The :math:`(n,2n)` channel
collapses identically. In the code
(:meth:`Mixture.condense
<orpheus.data.macro_xs.mixture.Mixture.condense>`) this reads exactly
``frame.project(mat @ T)`` per matrix channel.

.. warning::

   The sink-sum / source-average asymmetry is the subtle point of the
   whole operation, and it is **the structural difference from spatial
   homogenization**, which flux-weights *both* matrix axes
   (:eq:`sn-homogenization-matrix-collapse` runs ``project`` on a single
   axis with the source weight, but the spatial collapse never
   *sums* a sink axis — there is no spatial sink-summation because
   homogenization keeps the group structure). Three wrong collapses each
   produce a numerically distinct — and rate-breaking — coarse matrix,
   and each is a textbook variable-swap / missing-factor trap
   (vv-principles failure **Mode 2** and **Mode 3**):

   * **swap the axes** (flux-weight the *sink*, sum the *source*,
     ``project(SigS.T @ T)``-style) — wrong source/sink roles;
   * **sum both axes** (``T.T @ SigS @ T``) — drops the source
     flux-weight entirely;
   * **project both axes** (flux-weight the sink too) — this is exactly
     the ``homogenize`` behaviour, which is *wrong* for condensation
     (it would guard against "the implementer copied ``homogenize``
     verbatim").

   The regression gate
   ``tests.data.test_mixture_condense::TestG3ScatteringTwoAxisCollapse``
   reds on all three because its in-scatter-rate reference loop sums the
   sink and flux-averages the source *explicitly* (a hand-coded double
   ``for`` over fine groups — structurally independent of the production
   ``project``).

The fission spectrum is a pure birth-group sum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The emission spectrum :math:`\chi_g` is **not** a reaction rate — it is
a probability distribution (a simplex, :math:`\sum_g \chi_g = 1`; see
:eq:`emission-spectrum-simplex` in :ref:`theory-cross-section-data`).
Flux-weighting it would not preserve anything physical and could leave
the simplex. The rate-preserving choice is the **pure birth-group sum**

.. math::
   :label: energy-condensation-chi-collapse

   \chi_G \;=\; \sum_{g \in G} \chi_g
   \;=\; \chi \,@\, T,

.. (vv-status rationale) Representational identity: the χ collapse is a
   pure birth-group sum (the @T contraction = the table's birth-group-sum
   role), NOT flux-weighted — χ is a probability mass, so the coarse
   probability is the sum of the fine masses landing in G (Hébert
   Eq. 3.112, Stamm'ler VI(6a)). The simplex/null law it must satisfy is
   the data-layer ``Mixture`` invariant; this label is the collapse
   formula, not a separate solver claim.
.. vv-status: energy-condensation-chi-collapse documented

— the probability mass of the fine birth groups landing in coarse group
:math:`G`, summed (Hébert :cite:`Hebert2009` Eq. 3.112; Stamm'ler VI(6a)).
This **differs from spatial homogenization**, whose :math:`\chi`
collapse is a *production-weighted convex average across cells*
(:eq:`sn-homogenization-chi-collapse`): there are many fine cells
contributing different spectra to one coarse cell, so they must be
*mixed*; here there is one material whose birth groups are merely
*re-binned*, so the coarse spectrum is the partial sum of the fine
probability mass. The sum **preserves the simplex**:

.. math::
   :label: energy-condensation-chi-simplex-preservation

   \sum_G \chi_G
   \;=\;
   \sum_G \sum_{g\in G} \chi_g
   \;=\;
   \sum_g \chi_g
   \;=\; 1,

because the partition is a partition of unity over coarse groups (every
fine group's mass is counted exactly once;
:eq:`energy-condensation-partition-of-unity`). A flux-weighted
projection would give :math:`\sum_G\chi_G \ne 1`, destroying the
simplex — which is why :math:`\chi` is routed through the *table*
contraction ``χ @ T``, **not** through :meth:`frame.project
<orpheus.numerics.frame.FrameBase.project>`. The simplex / null law is
re-validated when the condensed
:class:`~orpheus.data.macro_xs.mixture.Mixture` is constructed
(:meth:`Mixture.__post_init__
<orpheus.data.macro_xs.mixture.Mixture.__post_init__>`).

.. note::

   Post fast-first flip (:ref:`canonical-group-convention`), coarse
   group ``0`` is the fastest, so a fission spectrum — which peaks in the
   fast range — is **fast-peaked**: the bulk of :math:`\chi_G` sits in
   the low coarse-group indices. (On the production 421-group grid the
   χ peak energy ≈ 1.15 MeV lands a few coarse groups in, not at index
   0, because the 20-MeV grid ceiling puts a small high-energy tail
   above the fission peak — so the physically-correct invariant the
   real-data gate pins is *cumulative* fast-half mass
   :math:`\sum_{G<G_{\max}/2}\chi_G > 0.5`, not ``argmax == 0``.)

Balance is preserved group-by-group
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The definitional total-XS balance every
:class:`~orpheus.data.macro_xs.mixture.Mixture` carries (the same
identity :eq:`sigT-computed` / :eq:`sn-homogenization-balance`),

.. math::
   :label: energy-condensation-balance

   \Sigma_t
   \;=\;
   \Sigma_c + \Sigma_L + \Sigma_f
   + \operatorname{rowsum}(\Sigma_{s0})
   + \operatorname{rowsum}(\Sigma_{2n}),

.. (vv-status rationale) Literature-transcribed definition: the total-XS
   balance identity every Mixture carries. Its preservation under
   condensation is the energy-axis copy of the homogenization argument —
   every removal channel collapses with the SAME per-coarse-group flux
   weight Φ_G. Gated by ``Mixture.assert_balanced`` on the condensed
   mixture (foundation invariant); not an SN solver claim.
.. vv-status: energy-condensation-balance documented

survives the collapse **group-by-group** when the fine material
balances, by the *same* one-line argument as homogenization. Fix a
coarse group :math:`G`. Every *removal* channel —
:math:`\Sigma_t,\ \Sigma_c,\ \Sigma_L,\ \Sigma_f`, and the row-sums
:math:`\sum_{g'}\Sigma_{s0}[g,g']`, :math:`\sum_{g'}\Sigma_{2n}[g,g']` —
is a *removal from group* :math:`g`, so each collapses with the **same**
source-group weight :math:`\varphi_g`. Crucially the matrix row-sum
:math:`\sum_{G'}\Sigma_{s0,G\to G'}` equals the *source-flux-average of
the fine row-sum*: the sink-sum :math:`\Sigma_{s0}\,T` then a row-sum
over coarse :math:`G'` telescopes (partition of unity) back to the fine
total out-scatter
:math:`\sum_{g'}\Sigma_{s0}[g,g']`, which then source-averages by the
same :math:`\varphi_g` as the vector channels. Because every term shares
the one weight :math:`\varphi_g`, the condensed balance is the *same
flux-weighted average* of the fine balances:

.. math::
   :label: energy-condensation-balance-preservation

   \Sigma_{t,G}
   - \Big(\Sigma_{c,G} + \Sigma_{L,G} + \Sigma_{f,G}
   + \operatorname{rowsum}(\Sigma_{s0,G})
   + \operatorname{rowsum}(\Sigma_{2n,G})\Big)
   \;=\;
   \frac{\sum_{g\in G} \varphi_g\,\big(\text{fine balance residual}_g\big)}
        {\sum_{g\in G} \varphi_g}
   \;=\; 0,

since each fine residual is zero. No "rebalance the condensed total"
step is needed — preservation is automatic, and the condensed
``Mixture`` passes :meth:`Mixture.assert_balanced
<orpheus.data.macro_xs.mixture.Mixture.assert_balanced>` whenever the
fine one does. This is the operation's correctness **oracle**: the
balance identity is a free, structurally-independent regression guard on
every condense. (Had the vector channels and the matrix row-sums
collapsed with *different* weights — e.g. had the matrix been projected
on both axes — the balance would break, another way of seeing why the
sink-sum / source-average asymmetry of
:eq:`energy-condensation-matrix-collapse` is *forced*, not chosen.)

.. _sn-condensation-fractional-overlap:

The non-nested problem: fractional-overlap re-binning
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Everything above assumed a fine→coarse membership table :math:`T[g,G]`.
For a **nested** coarse structure — coarse boundaries a subset of the
fine boundaries — that table is one-hot (:math:`T[g,G]\in\{0,1\}`, each
fine group wholly inside one coarse group) and the collapse reduces to
the exact group-sum of :eq:`energy-condensation-vector-collapse`. But the
production case is **not nested**, and the table must generalise.

Why one-hot containing-interval fails
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ORPHEUS condenses a **421-group** library onto the WIMS-D 69- and
172-group structures (:cite:`WIMSD`). These structures were defined
*independently* of the 421-group grid, so their boundaries do **not**
align with the fine grid (the draft boundary-mismatch report flags 19
non-coincident boundaries for 172→69 alone), and — the harder part —
**the coarse grid is locally finer than the fine grid** in narrow
resonance and thermal bands. A naïve "assign each fine group wholly to
the coarse group its representative energy falls in" (a one-hot
``searchsorted`` containing-interval rule) then leaves **empty coarse
groups**: 3 empty for 421→69, 22 empty for 421→172 — a coarse group
narrower than the fine spacing receives no fine group's representative
energy at all. An empty coarse group has zero total cross section, which
is unphysical and breaks a downstream solve.

This is a well-known stage distinction in reactor-physics data
processing. There are **two** group-averaging stages, and only the
second has the nesting constraint:

.. list-table::
   :header-rows: 1
   :widths: 26 40 34

   * - Stage
     - What it averages
     - Boundary constraint
   * - **pointwise → multigroup**
       (NJOY/GROUPR, OpenMC ``mgxs``, MC²-3)
     - the *continuous-energy* :math:`\sigma(E)\phi(E)` directly over any
       boundaries
     - **none** — the actual cross-section shape inside each group is the
       truth, so any structure is trivially integrable
   * - **multigroup → fewer-group**
       (AMPX/MALOCS — *the stage ORPHEUS is in*)
     - *already-discretised* fine groups via a fine→coarse map
     - **nesting** — a fine group cannot be split, because the input
       carries only the fine-group *average*, not the within-group shape

The collapse-stage codes (AMPX's MALOCS module) **require** the coarse
boundaries to be a subset of the fine boundaries: their input is a
fine→coarse *correspondence array* (e.g. "the first 4 fine groups → broad
group 1"), which structurally cannot express a fine group straddling a
broad boundary. ORPHEUS deliberately goes **beyond** MALOCS — it lifts
the nesting requirement with conservative fractional re-binning, a
capability the production deterministic-library codes mostly lack (they
sidestep it by re-integrating the continuum, which ORPHEUS cannot do from
a pre-grouped 421-library).

The fix: a fractional partition of unity
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A fine group :math:`g` that *straddles* a coarse boundary apportions a
**fraction** of its rate to each coarse group it overlaps. The membership
table becomes **fractional**,

.. math::
   :label: energy-condensation-overlap-fraction

   T[g,G] \;=\; f_{g,G}
   \;=\; \frac{\int_{g \cap G} w(E)\,\mathrm{d}E}
              {\int_g w(E)\,\mathrm{d}E},

.. (vv-status rationale) Literature-transcribed definition: the
   flux-weighted overlap fraction of a straddling fine group — the
   fraction of its within-group flux falling in coarse group G (the
   standard reactor-physics conservative re-binning, Q2 of the P5
   literature pull). The verifiable content is the partition-of-unity
   bit-identity (:eq:`energy-condensation-partition-of-unity`) and the
   1/E lethargy-ratio gate (``TestF4WithinGroupModelOracle``); a
   definition, not a solver claim.
.. vv-status: energy-condensation-overlap-fraction documented

the fraction of fine group :math:`g`'s
(:math:`w`-weighted) interval lying in coarse group :math:`G`, with
:math:`w(E)` an assumed **within-group flux model**. By the integral's
additivity the table is a **partition of unity**:

.. math::
   :label: energy-condensation-partition-of-unity

   \sum_G T[g,G] \;=\;
   \frac{\int_g w(E)\,\mathrm{d}E}{\int_g w(E)\,\mathrm{d}E} \;=\; 1
   \qquad \forall g.

.. (vv-status rationale) Structural identity: the membership table's rows
   sum to 1 (every fine group's rate is partitioned, not duplicated or
   dropped) — the property that makes the collapse conservative and the χ
   sum simplex-preserving. Gated by ``TestF2PartitionOfUnity``
   (``table.sum(axis=1) == 1``); a representational invariant, not a
   solver claim.
.. vv-status: energy-condensation-partition-of-unity documented

so each fine group's rate is *partitioned* (counted once), never
duplicated or dropped — the collapse stays conservative, and the χ-sum
stays a simplex. The general fractional collapse

.. math::
   :label: energy-condensation-fractional-collapse

   \Sigma_G
   \;=\;
   \frac{\sum_g f_{g,G}\,\varphi_g\,\Sigma_g}
        {\sum_g f_{g,G}\,\varphi_g}

reduces **exactly** to the one-hot group-sum
:eq:`energy-condensation-vector-collapse` when the structure is nested
(:math:`f_{g,G}\in\{0,1\}`) — so the nested case is the *regression-safe
degenerate*, not a separate code path.

The within-group flux model :math:`w(E)` is **selectable**
(:class:`~orpheus.data.energy_grid.WithinGroupSpectrum`, a strategy
Protocol). The default — built first — is **1/E (flat in lethargy)**,
:class:`~orpheus.data.energy_grid.InverseEnergySpectrum`:

.. math::
   :label: energy-condensation-lethargy-overlap

   \int_{lo}^{hi} \frac{\mathrm{d}E}{E} \;=\; \ln\!\frac{hi}{lo},
   \qquad
   f_{g,G} \;=\;
   \frac{\ln(hi_{g\cap G}/lo_{g\cap G})}{\ln(hi_g/lo_g)},

.. (vv-status rationale) Literature-transcribed definition: the 1/E
   (flat-in-lethargy) overlap fraction is a lethargy ratio — the
   asymptotic slowing-down spectrum, NJOY IWT=3, the standard first
   choice for condensation. The verifiable content is the
   ``InverseEnergySpectrum.integrated_weight`` = ln(hi/lo) bit-identity
   and the ``TestF4`` 1/E-vs-flat-energy discriminator; a definition, not
   a solver claim.
.. vv-status: energy-condensation-lethargy-overlap documented

the lethargy-overlap ratio (the asymptotic slowing-down spectrum; the
standard first choice for condensation, NJOY ``IWT=3``). Flat-in-energy
(NJOY ``IWT=2``) and the library weighting spectrum (fission + 1/E +
Maxwellian, NJOY ``IWT=4``) are future options on the same strategy
seam. The model is the **only new numerics surface** the non-nested case
adds — everything else (the Petrov-Galerkin frame, :meth:`frame.project
<orpheus.numerics.frame.FrameBase.project>`, the diagonal Gram, rate
preservation) is unchanged.

.. _sn-condensation-petrov-galerkin-frame:

Condensation is a Petrov-Galerkin projection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Like homogenization (:ref:`sn-homogenization-petrov-galerkin-frame`),
condensation is the coefficient extraction :math:`G^{-1}M` of a
**Petrov-Galerkin** frame — the *energy-axis instance* of the *same*
mechanism, which is exactly why the numerics core is reused verbatim:

  Condensation is :meth:`frame.project
  <orpheus.numerics.frame.FrameBase.project>` of a
  :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` whose *trial*
  basis is the fractional group-overlap indicator
  :math:`\mathbf{1}_G` (carried by
  :class:`~orpheus.numerics.basis.OverlapBasis`), whose *test* basis is
  the spectrum-weighted indicator :math:`\varphi\,\mathbf{1}_G` (carried
  by :class:`~orpheus.numerics.basis.WeightedIndicatorBasis`), and whose
  measure is the **counting** measure :math:`\mu` (weight 1).

This is not decoration: it is why :meth:`Mixture.condense
<orpheus.data.macro_xs.mixture.Mixture.condense>` routes through the
*same* discrete :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` that
carries SN anisotropic-scattering moment projection and spatial
homogenization — one mechanism for every fine→coarse change of
representation (Cardinal Rule 2, single source of truth), not a bespoke
membership matmul per axis. The frame projection machinery
(:class:`~orpheus.numerics.frame.PetrovGalerkinFrame`,
:class:`~orpheus.numerics.basis.OverlapBasis`,
:class:`~orpheus.numerics.basis.WeightedIndicatorBasis`,
:meth:`FunctionSpace.apply_inverse_metric
<orpheus.numerics.space.FunctionSpace.apply_inverse_metric>`) was built
for homogenization anticipating energy as the second consumer; the data
layer reaches it because ``data → numerics`` is a permitted layer edge.

The trial / test / measure separation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Three distinct objects carry the three roles, cleanly separated — the
campaign's three-way trial/test/measure split holding verbatim
on the energy axis:

.. list-table::
   :header-rows: 1
   :widths: 20 30 28 22

   * - Role
     - Object
     - Carries
     - On the energy axis
   * - **trial** :math:`K`
     - :class:`~orpheus.numerics.basis.OverlapBasis`
     - the fractional membership :math:`T[g,G]` (partition geometry +
       within-group split)
     - the *only* new surface — the fractional table
   * - **test** :math:`M^*`
     - :class:`~orpheus.numerics.basis.WeightedIndicatorBasis`
     - the spectrum :math:`\varphi_g` as an analysis weight
     - the flux is the test weight, **never** the measure
   * - **measure** :math:`\mu`
     - :class:`~orpheus.numerics.measure.DiscreteMeasure`
     - the counting metric (:math:`w_g=1`), group-independent
     - fixed :math:`L^2`, **never** the discipline

The test functions :math:`\varphi_g\,\mathbf{1}_G` differ from the trial
functions :math:`\mathbf{1}_G` (:math:`\varphi_g\,\mathbf{1}_G \ne
\mathbf{1}_G`), so the projection is genuinely **Petrov-Galerkin**,
carried by the frame **type**
(:class:`~orpheus.numerics.frame.PetrovGalerkinFrame`). With trial and
test indicators of (fractional) partition-of-unity support, the
**cross-Gram** is **diagonal** — the all-ones probe through the frame
gives :math:`(M\,\mathbf 1)_G = \sum_g \varphi_g\,T[g,G] = \Phi_G`, the
coarse-group flux :eq:`energy-condensation-coarse-flux`, which *is* the
diagonal of the Gram — so the normal equations decouple group-by-group
and :meth:`frame.project <orpheus.numerics.frame.FrameBase.project>`
returns the spectrum-weighted average
:eq:`energy-condensation-vector-collapse` with a per-group
**reciprocal**, not a linear solve.

.. note::

   For a *fractional* (straddling) table two coarse columns share a fine
   row, so the **off-diagonal** cross-Gram
   :math:`G_{GG'} = \sum_g \varphi_g\,T[g,G]\,T[g,G']` is **not**
   structurally zero. This is correct and harmless:
   :meth:`frame.project <orpheus.numerics.frame.FrameBase.project>` uses
   only the **diagonal** :math:`G_{GG} = \Phi_G` (each coarse group's
   rate is its *own* functional — one DOF per group, the P0 / rank-0
   space), so it ignores the off-diagonals by construction. The
   :meth:`OverlapBasis.mass_matrix
   <orpheus.numerics.basis.indicator_basis.IndicatorBasis.mass_matrix>`
   inherits a docstring claiming a *diagonal* Gram (true for the one-hot
   parent, false for a fractional table) — that claim is **latent**: no
   consumer calls ``mass_matrix``, and the frame's row-sum probe never
   forms the full Gram. A *future* least-squares consumer needing the
   dense Gram (a non-indicator, richer coarse basis — which cross
   sections never want, a P1 coarse XS is not rate-meaningful) must
   compute it for the fractional case, not trust the inherited diagonal
   claim. See :ref:`frame-least-squares-discipline` for why this is
   **not** a ``LeastSquaresFrame`` (its trigger — test = :math:`A`·trial,
   a dense SPD Gram needing a real solve — is absent here; that
   discipline is designed-but-not-built, gated by
   :class:`~orpheus.numerics.basis.GramStructure` ``DENSE``).

Why the spectrum is the test weight, not a measure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The natural-seeming alternative — "treat the spectrum as a measure and
project the cross section onto a coarse basis in the
:math:`L^2(\text{spectrum})` metric" — is the energy-axis form of the
metric-fold (:eq:`sn-homogenization-metric-fold`), and it is **refused**
for the *same* reason it is refused for homogenization
(:ref:`sn-homogenization-why-petrov-galerkin`, the #268 ruling: *the
measure carries the axis and the fixed* :math:`L^2` *metric, never the
discipline*). Three structural grounds, identical to the spatial case:

#. **For a P0 (indicator) coarse basis the least-squares fit
   *coincides* with the flux-weighted average** — with disjoint /
   fractional partition-of-unity indicators the normal-equations Gram is
   diagonal and the least-squares solution is
   :math:`\Sigma_G = \sum_g w_g\Sigma_g / \sum_g w_g`, the flux-weighted
   average verbatim when :math:`w=\varphi`. So "least-squares" does not
   select a *new* frame — it re-derives the same Petrov-Galerkin average
   under a different name.
#. **Folding :math:`\varphi` into the measure breaks under
   adjoint-weighting.** The eigenvalue-consistent condensation reactor
   physics ultimately requires preserves the **bilinear**
   :math:`\langle\varphi^*,\Sigma\varphi\rangle` with
   :math:`\varphi^*\ne\varphi` (the same
   :eq:`sn-homogenization-bilinear` structure on the energy axis), where
   test :math:`= \varphi^*\mathbf{1}_G \ne` trial :math:`=
   \varphi\,\mathbf{1}_G` and **no single metric** reproduces the
   two-sided weighting. Forward-flux reaction-rate-only condensation is
   the degenerate :math:`\varphi^*=\varphi` case where the fold happens
   to work; the *type* (an explicit test basis) encodes the general
   case, so the adjoint-weighted lift (a later phase of the unified
   Frame-projection campaign, P6) lands as a no-op change of the test
   weight (:math:`\varphi \to \varphi^*`), not a re-derivation.
#. **The discipline is a property the** *type* **carries, never the
   measure.** Keeping :math:`\varphi` on the
   :class:`~orpheus.numerics.basis.WeightedIndicatorBasis` test side
   forces the correct reading: :math:`\varphi` is a test-weighting the
   *solution* emits, not a property of the energy axis. The energy axis
   (the grid) owns one counting measure; the solution owns the spectrum
   :math:`\varphi`; the frame *type* (Petrov-Galerkin, with its explicit
   test basis) carries the discipline.

The flux-weighted average is the **rank-0 moment** of **Generalized
Energy Condensation** (Rahnema, Douglass & Forget 2008 :cite:`Rahnema2008`):
GEC expands the within-coarse-group flux in orthogonal functions
:math:`\varphi(E)\approx\sum_n\varphi_{n,G}P_n(E)`, and the zeroth
moment (the constant / piecewise-constant basis function on :math:`G`)
*recovers the standard flux-weighted multigroup average exactly* — "the
zeroth moment generates the standard few-group equation". The higher
moments (:math:`n\ge1`) add the within-coarse-group spectral detail the
simple average discards; that is faithful within-group reconstruction
(honest upscaling), and it is **deferred** — it would be a richer trial
basis on the *same* frame (`GitHub #275
<https://github.com/deOliveira-R/ORPHEUS/issues/275>`_), no architectural
change.

.. _sn-condensation-downsampling:

Downsampling only: condensation loses information
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The load-bearing semantic — a design ruling, not an implementation
detail — is that **condensation is a one-way, lossy, downsampling
operation.** Collapsing fine groups into coarse groups *discards* the
within-coarse-group spectral structure (that is the rank-:math:`>0` GEC
content above). The continuous-projection view *could* fabricate detail
(a "64 → 200 group" upscaling would invent sub-group structure the data
never carried), but the group-collapse stage ORPHEUS is in **cannot**:
the input is only the fine-group *averages*. The asymmetry is encoded in
three places:

* **A global upscaling guard.**
  :meth:`EnergyGrid.overlap_to
  <orpheus.data.energy_grid.EnergyGrid.overlap_to>` (the binary
  mismatch factory) **refuses** a coarse target with *more* groups than
  the source — :math:`n_{\rm coarse} > n_{\rm fine}` raises
  ``ValueError`` ("condensation only DOWNSAMPLES; a finer target would
  fabricate sub-group structure the data does not contain"). Reconstructing
  a finer structure from group-integrated data is fabrication; the guard
  makes it unrepresentable.
* **The within-group model is the *explicit* assumption.** Where the
  coarse grid is *locally* finer than the fine grid (the resonance /
  thermal bands), the unavoidable *local* interpolation is done by the
  within-group flux model :math:`w(E)`
  (:eq:`energy-condensation-overlap-fraction`) — a bounded, named,
  selectable assumption, not a silent invention.
* **The provenance report.**
  :attr:`OverlapBasis.fractional_columns
  <orpheus.numerics.basis.overlap_basis.OverlapBasis.fractional_columns>`
  (carried on the trial basis :meth:`overlap_to
  <orpheus.data.energy_grid.EnergyGrid.overlap_to>` returns)
  lists the coarse-group indices whose data leaned on :math:`w(E)`
  (the columns that received a *fractional* — strictly between 0 and 1 —
  contribution). It is **empty** for a nested condensation (pure
  rate-preserving collapse, no assumption) and non-empty exactly where
  the coarse grid is locally finer than the fine grid. This is the
  data-vs-assumption provenance: a caller can see precisely which coarse
  groups are pure collapse and which lean on the spectral model. (The
  companion :attr:`OverlapBasis.dominant_column
  <orpheus.numerics.basis.overlap_basis.OverlapBasis.dominant_column>` —
  the ``argmax`` containing-coarse map — is the former
  ``GroupCondensation.coarse_of_fine``.)

.. warning::

   Faithful reconstruction / honest *upscaling* (recovering
   within-coarse-group detail via the rank-:math:`>0` GEC moments) is
   **not** what this slice ships, and the upscaling guard deliberately
   prevents accidentally posing it as a fine-target ``condense``. Do
   **not** read the within-group model :math:`w(E)` as faithful
   reconstruction — it is a bounded local-interpolation assumption used
   *only* to apportion a straddling fine group's already-known rate, not
   to invent new spectral structure. Honest upscaling is a future
   capability (`GitHub #275 <https://github.com/deOliveira-R/ORPHEUS/issues/275>`_).

.. _sn-condensation-grid-frame-axis:

The grid is a frame axis: dual measure / basis views
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

An :class:`~orpheus.data.energy_grid.EnergyGrid` is the energy analogue
of a coarse :class:`~orpheus.geometry.mesh.Mesh1D`, and — exactly like a
mesh — it is a **frame axis** that yields *both* halves of a discrete
frame, the two roles an axis plays in a projection:

.. list-table::
   :header-rows: 1
   :widths: 22 30 26 22

   * - View
     - Method
     - Role
     - Spatial twin
   * - **measure** :math:`\mu`
     - :meth:`EnergyGrid.as_measure
       <orpheus.data.energy_grid.EnergyGrid.as_measure>`
     - the **source** — the axis you project *from* (a counting
       :class:`~orpheus.numerics.measure.DiscreteMeasure`, :math:`w_g=1`,
       ``support="energy"``)
     - :meth:`Mesh1D.volume_measure
       <orpheus.geometry.mesh.Mesh1D.volume_measure>`
   * - **basis** :math:`\mathbf{1}`
     - :meth:`EnergyGrid.as_basis
       <orpheus.data.energy_grid.EnergyGrid.as_basis>`
     - the **target** — the axis you project *to* (the group-indicator
       :class:`~orpheus.numerics.basis.IndicatorBasis`, one-hot)
     - :meth:`Mesh1D.indicator_basis
       <orpheus.geometry.mesh.Mesh1D.indicator_basis>`

So a *nested* condensation is just ``fine.as_measure()`` →
``coarse.as_basis()`` — the two unary views suffice. The non-nested
production case needs **one more thing**, and it is irreducibly
**binary**: the fractional membership table reads *both* grids' edges at
once (a fine group straddles a *coarse* boundary — neither grid alone
knows the straddle), so it cannot be a unary view of either. That is
:meth:`EnergyGrid.overlap_to
<orpheus.data.energy_grid.EnergyGrid.overlap_to>`, the
``(fine, coarse) → OverlapBasis`` factory, with the containment

.. math::
   :label: energy-condensation-nested-subset

   \texttt{coarse.as\_basis()}\ \text{(nested, one-hot)}
   \ \subset\
   \texttt{fine.overlap\_to(coarse)}\ \text{(non-nested, fractional)}.

.. (vv-status rationale) Structural identity: the nested one-hot target
   view (``as_basis``) is the degenerate of the binary fractional
   ``overlap_to`` — the same containment the partition-of-unity table
   collapses to when every straddle fraction is 0 or 1. A
   representational relationship between the two views, gated by the
   ``TestF3NestedDegeneracy`` bit-identity (a nested ``overlap_to`` table
   equals the ``searchsorted`` one-hot); not a solver claim.
.. vv-status: energy-condensation-nested-subset documented

The returned trial is an :class:`~orpheus.numerics.basis.OverlapBasis`,
which **IS-A** :class:`~orpheus.numerics.basis.IndicatorBasis` carrying
the fractional table — so the nested one-hot view and the non-nested
fractional view are the *same type* of object, the degenerate and the
general case, never two code paths.

.. _sn-condensation-no-frame-subclass:

Why no ``CondensationFrame`` — the data-native shape
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Two design decisions shape where this machinery lives, and both are the
kind of "what was tried and rejected" a future session needs spelled
out so it does not re-litigate them.

**Condensation is data-native (no transport dependency).** An earlier
plan put a ``CondensationFrame`` in :mod:`orpheus.transport` (a
``transport/frames/`` package, symmetric with the angular
:class:`~orpheus.numerics.frame.GalerkinFrame`). That was **overturned**.
Condensation's carrier is the :class:`~orpheus.data.macro_xs.mixture.Mixture`
— a *data* type — and every object the collapse touches (the
:class:`~orpheus.data.energy_grid.EnergyGrid`, the
:meth:`overlap_to <orpheus.data.energy_grid.EnergyGrid.overlap_to>`
factory, and the :class:`~orpheus.numerics.frame.PetrovGalerkinFrame`)
lives in ``data`` / ``numerics``. The layering forbids the move (``data``
must **not** depend on ``transport``), and nothing in the operation needs
transport: it is a pure cross-section *re-binning*. So the collapse verb
is :meth:`Mixture.condense
<orpheus.data.macro_xs.mixture.Mixture.condense>` (data), reaching the
frame through the permitted ``data → numerics`` edge; only the SN-layer
*orchestration* (deriving the per-material spectrum from a solved flux)
lives above, in :meth:`Solution.condense
<orpheus.sn.solution.Solution.condense>`.

**There is no frame subclass at all** — no ``CondensationFrame``, no
``HomogenizationFrame``. The frame is a plain
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame`; the "condensation"
identity is **not** a new kind of frame, it lives in *two* ordinary
places:

#. the **binary overlap factory** :meth:`EnergyGrid.overlap_to
   <orpheus.data.energy_grid.EnergyGrid.overlap_to>` (which builds the
   energy-specific fractional trial), and
#. the **collapse verb** :meth:`Mixture.condense
   <orpheus.data.macro_xs.mixture.Mixture.condense>`, which shares its
   channel-assembly assembler :meth:`Mixture.from_dense_channels
   <orpheus.data.macro_xs.mixture.Mixture.from_dense_channels>` with
   spatial homogenization's :meth:`MaterialXSField.project_through
   <orpheus.transport.mesh.material_xs_field.MaterialXSField.project_through>`
   — one home for "assemble a ``Mixture`` from coarsened dense channels",
   not one per verb (Cardinal Rule 2).

Minting a ``CondensationFrame`` would be wrong on two counts. First it is
**false symmetry** with homogenization: homogenization is intrinsically a
*two-frame* operation (a flux-weighted frame for the rate channels **and**
a production-weighted frame for :math:`\chi` — :meth:`project_through
<orpheus.transport.mesh.material_xs_field.MaterialXSField.project_through>`
takes both), so even there a single "HomogenizationFrame" is the wrong
shape. Second it is **unjustified type-minting** (the project's
type-vs-property rule): a new frame *type* is earned only by a new frame
*morphism*, and condensation introduces none — its analysis,
reconstruction, Gram, and :meth:`project
<orpheus.numerics.frame.FrameBase.project>` are the *unchanged*
Petrov-Galerkin operations. The only genuinely new surface is the
fractional trial *basis* (:class:`~orpheus.numerics.basis.OverlapBasis`),
which is exactly where the novelty belongs — a basis, not a frame.

.. _sn-condensation-verification:

Verification
~~~~~~~~~~~~

The gates are :mod:`tests.data.test_energy_grid`,
:mod:`tests.data.test_mixture_condense`, and
:mod:`tests.sn.test_condensation`. The two solver-facing claims carry
``@pytest.mark.verifies`` markers tying them to the equations above:

* ``energy-condensation-rate-preservation``
  (:eq:`energy-condensation-rate-preservation`) — the vector-channel
  rate-preservation identity, asserted for every channel
  (:math:`\Sigma_t,\ \Sigma_c,\ \Sigma_L,\ \Sigma_f,\ \nu\Sigma_f`) and
  every coarse group: :math:`\Sigma_G\,\Phi_G = \sum_{g\in G}
  \varphi_g\,\Sigma_g` to one ULP.
* ``energy-condensation-scattering-collapse``
  (:eq:`energy-condensation-scattering-collapse`) — the two-axis matrix
  collapse, asserted as the preserved in-scatter rate :math:`\Phi_G\,
  \Sigma_{s,\ell,G\to G'} = \sum_{g\in G}\sum_{g'\in G'}\varphi_g\,
  \Sigma_s[g,g']` for every Legendre order and the :math:`(n,2n)` matrix.

These are **L1** (equation) claims against a **closed-form** reference,
*not* eigenvalue claims — condensation is a data-reduction operation, not
a solve, so there is no :math:`k` to verify; the pillar question is "does
the reduction preserve the rate functional", answered by closed-form
hand-summation. The correctness oracle is a **structurally-independent**
explicit per-group Python loop over the fine groups — *not* a re-call of
the production :meth:`frame.project
<orpheus.numerics.frame.FrameBase.project>` (vv-principles **L11**: a
cross-check must be structurally independent, not merely procedurally
independent; a frame-vs-frame comparison would share any reduction bug).
The scattering rate is one ULP, not bit-identical, because the
``@ T`` matmul reduction tree differs from the explicit group-by-group
sum — FP-non-associativity, principled-equivalent per the `vv-principles`
bit-identity-vs-principled-equivalence criteria (drift = reduction-depth
× ULP); the gate uses ``np.testing.assert_allclose(rtol=1e-12)``, never
exact ``==``.

The discriminator and the companion invariants:

.. list-table::
   :header-rows: 1
   :widths: 48 52

   * - Test
     - What it pins
   * - ``TestG1RatePreservationVectors`` /
       ``TestF1StraddleRatePreservation``
     - The **rate-preservation anchor**
       (:eq:`energy-condensation-rate-preservation`, the
       ``energy-condensation-rate-preservation`` verifies-target):
       :math:`\Sigma_G\,\Phi_G = \sum_{g\in G}\varphi_g\,\Sigma_g` for
       every vector channel and every coarse group, against the hand-sum
       oracle — nested (G1) and straddling-fractional (F1).
   * - ``TestG2WithinGroupVaryingDiscriminator`` /
       ``TestF4ModelDiscriminatorSUT``
     - The load-bearing **flux-weighting discriminator** (vv Mode 7): a
       fine spectrum that *varies within* each coarse group (e.g.
       :math:`\varphi=[1,4,2,0.5]`) makes the flux-weighted and
       arithmetic-average collapses numerically distinct, and rate
       preservation MUST match the flux-weighted one — reds a regression
       that drops :math:`\varphi`. A *flat* spectrum would null the
       weighting (the Mode-7 trap), so the fixture is asserted
       within-group-varying.
   * - ``TestG3ScatteringTwoAxisCollapse`` (three mutations)
     - The sink-sum / source-average asymmetry
       (:eq:`energy-condensation-matrix-collapse`) — each of swap-axes /
       sum-both / project-both produces a *numerically different* coarse
       matrix → each reds (vv Mode 2 / Mode 3; the project-both mutation
       guards against copying ``homogenize`` verbatim).
   * - ``TestChiBirthGroupSum`` (χ sum-vs-project guard)
     - :math:`\chi_G = \chi\,@\,T` preserves :math:`\sum\chi=1`
       (:eq:`energy-condensation-chi-collapse`); a flux-*projected* χ
       sums to :math:`\ne 1`, destroying the simplex — pinned separately
       from the projected channels.
   * - ``TestF2PartitionOfUnity`` / ``TestF3NestedDegeneracy``
     - The table is a partition of unity
       (:eq:`energy-condensation-partition-of-unity`,
       ``rows.sum == 1``), and the fractional table reduces
       **bit-identically** to the one-hot ``searchsorted`` table for a
       nested structure (the regression-safe degenerate).
   * - ``TestF4WithinGroupModelOracle``
     - The 1/E overlap fraction is the lethargy ratio
       :math:`\ln(hi_{\cap}/lo_{\cap})/\ln(hi_g/lo_g)`
       (:eq:`energy-condensation-lethargy-overlap`), and 1/E ≠
       flat-energy on a straddling group (the model is load-bearing, not
       cosmetic).
   * - ``TestF5UpscalingGuard``
     - :math:`n_{\rm coarse} > n_{\rm fine}` raises ``ValueError`` (the
       downsampling-only guard), with a positive control that a valid
       downsample succeeds.
   * - ``TestF6LocalInterpolationReport``
     - :attr:`OverlapBasis.fractional_columns
       <orpheus.numerics.basis.overlap_basis.OverlapBasis.fractional_columns>`
       is empty for a nested condensation and lists the straddle columns
       for a non-nested one (the data-vs-assumption provenance).
   * - ``TestG4BalanceRegression`` (positive + negative)
     - A condensed balanced ``Mixture`` passes
       :meth:`Mixture.assert_balanced` (:eq:`energy-condensation-balance`);
       a hand-built rate-broken condensed ``Mixture`` raises (vv #11:
       the negative leg pins the *invariant*, not merely the raising).
   * - ``TestG5WimsDerivationValidation`` (Table 11.3)
     - The containing-interval partition derived by the rule reproduces
       the published ``CONDENSE_172_TO_69`` (:cite:`WIMSD` Table 11.3) on the
       coincident-boundary groups, collecting the known 19 non-coincident
       boundaries as expected (failing only on a *new* mismatch).
   * - ``test_real_pwr_421_to_wims69_condensation_succeeds``
     - **L2 integration**: a *real* 421-group production mixture
       (:func:`~orpheus.data.macro_xs.recipes.pwr_like_mix`) condenses to
       WIMS-69 with **no empty coarse groups** (:math:`\Sigma_t>0` ∀
       coarse — the one-hot empty-group bug *gone*), balance preserved,
       χ fast-half-mass :math:`>0.5`, and a non-empty
       ``locally_interpolated``.
   * - Mode-11 routing sentinel
     - ``Mixture.condense`` actually calls :meth:`frame.project
       <orpheus.numerics.frame.FrameBase.project>` and
       :meth:`WeightedIndicatorBasis.analyze
       <orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis.analyze>`
       (the **test-side** spectrum reader) — the Petrov-Galerkin routing
       is on the gate's call graph, not bypassed by an inline matmul (vv
       **Mode 11**).

The intrinsic-property tests pin the new value objects' defining laws:
:class:`~orpheus.data.energy_grid.EnergyGrid` (strictly **descending**
boundaries — the #265 monotonicity slice — all-positive energies,
partition completeness, with positive *and* negative legs), and the
fractional-overlap trial
:class:`~orpheus.numerics.basis.OverlapBasis` (the partition-of-unity
law, ``rows.sum == 1``; and the containing-interval law on its
:attr:`~orpheus.numerics.basis.overlap_basis.OverlapBasis.dominant_column`:
every fine group → exactly one *dominant* coarse group, contiguous,
**fast-first** ordering — the orientation pin that catches the silent
descending-edge column-reversal trap, vv Mode 6).

≥2 groups throughout (69-, 172-, and 421-group cases — never the
degenerate 1-group case). Every sentinel / negative leg uses
``np.testing.*`` / ``pytest.raises`` / ``pytest.fail``, never a bare
``assert`` (vv **Mode 8**: ``-O`` strips bare asserts).

.. _frame-galerkin-frame:

The Galerkin frame
==================

The **Galerkin** frame is the special case ``test is trial`` — the
:class:`~orpheus.numerics.frame.GalerkinFrame` specialisation of the
Petrov-Galerkin base above. It strengthens the base promise
:math:`M R = I_W` (up to tightness) to the self-dual :math:`M^* = R`
(under an orthonormal trial basis), and its canonical instance is the
angular spherical-harmonic frame.

.. _frame-galerkin-in-general:

In general
----------

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
      :label: galerkin-strict-adjoint-vs-reconstruction

      (M^* c)_n
      &\;=\; \sum_{\ell, m} Y_\ell^m(\hat\Omega_n)\,c_\ell^m
        \quad\text{(no factor — strict adjoint)}, \\
      (R c)_n
      &\;=\; \sum_\ell (2\ell+1)\,\sum_m Y_\ell^m(\hat\Omega_n)\,
             c_\ell^m
        \quad\text{(with factor — addition-theorem)}.

   .. (vv-status rationale) Representational identity: distinguishes the strict
      Hilbert adjoint M* (naked synthesis, no 2ℓ+1) from the reconstruction
      face R (with 2ℓ+1) — the ERR-039 distinction. Each face is verified under
      its own label: M* = g_C·S_0 by :eq:`hilbert-adjoint-equals-metric-times-S0`,
      the (2ℓ+1) synthesis by :eq:`sh-addition-theorem-reconstruction`. A
      face-distinction framing, not a separate solver claim.
   .. vv-status: galerkin-strict-adjoint-vs-reconstruction documented

   The analysis face's representation transpose
   :meth:`frame.analysis.apply_transpose
   <orpheus.numerics.basis.Basis.analyze_transpose>` is
   :math:`w_n\,S_0` (the *naked* synthesis weighted by the
   quadrature weight); its metric-aware Hilbert adjoint
   ``frame.analysis.H`` is :math:`M^* = g_C\,S_0`; and
   :meth:`frame.reconstruction.apply
   <orpheus.numerics.basis.Basis.reconstruct>` returns :math:`R`
   (with the :math:`(2\ell+1)` factor). The adjoint-face dishonesty
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


.. _frame-spherical-harmonic-galerkin:

Applied to spherical-harmonics projection
-----------------------------------------

The spherical-harmonic angular frame is the canonical
:class:`~orpheus.numerics.frame.GalerkinFrame`. Its Galerkin discipline
is *forced*, not chosen: the anisotropic scattering source is a rotation-
invariant (zonal) integral kernel whose eigenbasis — by Funk–Hecke — is
the spherical harmonics, and the eigenbasis of a self-adjoint operator is
orthogonal. Every term of the frame is named below: the analysis face
:math:`M` (the quadrature-weight :math:`w_n` projection onto the harmonic
moments), the reconstruction face :math:`R` (the :math:`(2\ell+1)`
addition-theorem synthesis), the diagonal eigenvalue operator
:math:`\Lambda` (the Legendre moments :math:`\Sigma_{s,\ell}`), and the
tightness constant :math:`c_V = 4\pi` for which :math:`M R = 4\pi I`
(:eq:`pi-r-equals-4pi-i`, derived in the
:ref:`general treatment of the Galerkin frame <frame-galerkin-in-general>`
above). The whole kernel is the spectral theorem
:math:`S = R\circ\Lambda\circ M = U\Sigma U^*` written out.

The anisotropic scattering source operator is an angular **integral
kernel** (:ref:`integral-kernel-category`, :doc:`/theory/foundations/operator_algebra`):
the source at :term:`ordinate` :math:`\hat\Omega` reads the flux at *every*
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
(:eq:`scattering-as-tensor-product-sum`, :doc:`/theory/foundations/operator_algebra`). The
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

.. (vv-status rationale) Representational identity: the anisotropic scattering
   kernel written as the spectral theorem S = R∘Λ∘M = UΣU*. The implementing
   kernel R∘Λ∘M is pinned by the 0-ULP windowed-vs-full crosscheck
   ``tests/sn/operators/test_scattering_kernel_crosscheck.py`` and the
   addition-theorem identity :eq:`real-sh-addition-theorem` (same kernel as the
   sibling :eq:`scattering-zonal-kernel`). A spectral-theorem framing, not a
   separate solver claim.
.. vv-status: scattering-spectral-theorem documented

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
of the operator-algebra double category, :doc:`/theory/foundations/operator_algebra`) a
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

.. _frame-composed-verbs:

The frame's composed-operator verbs
===================================

A consumer does **not** hand-roll the analysis / reconstruction faces.
The point of binding a basis to a measure through a
:class:`~orpheus.numerics.frame.FrameBase` is that the frame then emits
the **composed operators** the method actually applies — *define a
frame, compose, done* (Cardinal Rule 2: the composition **is** the
production path, not a parallel "semantic" reading layered over a
hand-rolled numpy chain). Three composed verbs cover every consumer:

.. list-table:: The frame's composed-operator verbs
   :header-rows: 1
   :widths: 30 20 50

   * - Verb
     - Composition
     - Consumer
   * - :meth:`conjugate(A) <orpheus.numerics.frame.FrameBase.conjugate>`
     - :math:`R \circ A \circ M`
     - SN anisotropic scattering
       :math:`S_{\ell\ge 1} = R\,\Lambda\,M` (project to moments,
       multiply by the spectrum :math:`\Sigma_{s,\ell}`, reconstruct
       the per-ordinate source)
   * - :meth:`reconstruct_after(A) <orpheus.numerics.frame.FrameBase.reconstruct_after>`
     - :math:`R \circ A`
     - inputs **already** in coefficient space — the angular-windowed
       SN moment iterate, whose bulk is already :math:`M\psi`, so only
       :math:`R\,\Lambda` remains (wiring it to ``conjugate`` would
       double-project)
   * - :meth:`project(f) <orpheus.numerics.frame.FrameBase.project>`
     - :math:`G^{-1} M`
     - the **homogenise / condense** coefficient extraction —
       :meth:`Solution.homogenize
       <orpheus.sn.solution.Solution.homogenize>`,
       :meth:`Solution.condense <orpheus.sn.solution.Solution.condense>`

Each returns a **typed**
:class:`~orpheus.numerics.operator.OperatorProduct` (or, for
``project``, the inverse-Gram ∘ analysis chain), whose ``apply`` runs
exactly the numpy contraction a hand-rolled
``reconstruction.apply(A.apply(analysis.apply(x)))`` would — now as
**one named operator** with the
:class:`~orpheus.numerics.operator.OperatorProduct` space-compatibility
guard enforcing that :math:`A` composes between the faces (its
``domain`` is the analysis codomain, its ``codomain`` the
reconstruction domain).

:meth:`conjugate <orpheus.numerics.frame.FrameBase.conjugate>` is the
**2-cell** of the (Representation × Role) carrier double category
(:ref:`operator-algebra`): a coefficient-space Role-morphism :math:`A`
conjugated by the horizontal Representation-adjoint pair
:math:`(M, R)`. When the frame is the operator's *eigenbasis* — the SH
angular frame is the scattering kernel's, by Funk–Hecke —
:math:`R\circ\Lambda\circ M` **is** the spectral theorem
:math:`U\Sigma U^*` written out, and the frame is then *owned* by that
operator (:ref:`frame-eigenbasis-ownership`). The
coefficient-extraction verb :meth:`project
<orpheus.numerics.frame.FrameBase.project>` is the Petrov-Galerkin
:math:`G^{-1}M` derived term-by-term for the homogenisation consumer in
:ref:`sn-homogenization-petrov-galerkin-frame`; its diagonal-Gram normalisation
:attr:`gram <orpheus.numerics.frame.FrameBase.gram>` is the row-sum
probe of :ref:`frame-least-squares-discipline`.


.. _frame-least-squares-discipline:

The least-squares discipline — designed, not built
==================================================

The discipline split is carried one level deeper than the
Galerkin / Petrov-Galerkin *type* by the **trial basis's Gram
structure** — the declaration
:class:`~orpheus.numerics.basis.GramStructure` that decides whether the
coefficient extraction :meth:`project
<orpheus.numerics.frame.FrameBase.project>` can use a cheap row-sum
probe or needs a full dense solve. ``project`` normalises by the
cross-Gram :math:`G = MR`; the frame computes it with a single
``analysis(reconstruction(ones))`` **row-sum probe**, but that probe
equals the required normalisation only under one of two structural
conditions:

.. list-table:: The trial-basis Gram structure decides the projection machinery
   :header-rows: 1
   :widths: 22 30 26 22

   * - ``GramStructure``
     - Trial Gram :math:`MR`
     - Projection normalisation
     - Built?
   * - ``DIAGONAL``
     - diagonal (orthogonal harmonics; disjoint / nested cell / group
       indicators)
     - the row sum **is** the diagonal — a reciprocal
     - **yes** (Galerkin SH; forward homogenisation)
   * - ``PARTITION_OF_UNITY``
     - not diagonal, but membership rows sum to 1
       (:class:`~orpheus.numerics.basis.OverlapBasis`)
     - :math:`R\mathbf 1 = \mathbf 1` collapses the probe to the
       per-region weight — still a reciprocal
     - **yes** (forward condensation)
   * - ``DENSE``
     - neither (a tapered weight, a higher-rank GEC moment)
     - needs the real :math:`(MR)^{-1}M` least-squares solve
     - **no** — :meth:`project` *refuses* (#275)

The first two rows are the **built** frames: a
:class:`~orpheus.numerics.frame.GalerkinFrame` (diagonal Gram,
``test is trial``) and the forward
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame` (diagonal *or*
partition-of-unity). For both, :attr:`FrameBase.gram
<orpheus.numerics.frame.FrameBase.gram>` is a single row-sum probe and
:meth:`project <orpheus.numerics.frame.FrameBase.project>` is a per-cell
reciprocal (a Moore–Penrose pseudo-inverse, so an empty / zero-flux
region maps to :math:`\Sigma_R = 0` for free), **not** a linear solve.

The third row is the **third discipline** — a least-squares frame over a
**dense** cross-Gram. It is the natural sibling of
:class:`~orpheus.numerics.frame.GalerkinFrame` under the Petrov-Galerkin
base (the designed hierarchy is ``FrameBase → PetrovGalerkinFrame →
{GalerkinFrame, LeastSquaresFrame}``): its trigger is a trial basis
whose :math:`MR` is genuinely dense — ``test`` :math:`= A\cdot`\ ``trial``
for some non-identity :math:`A`, a dense SPD Gram needing a real solve —
for which the row-sum probe is **wrong**. It is **designed but not
built**: the base :class:`~orpheus.numerics.basis.Basis` defaults to
``GramStructure.DENSE`` (the safe refusal), and
:meth:`FrameBase.project <orpheus.numerics.frame.FrameBase.project>`
raises :class:`~orpheus.numerics.operator.NotInvertible` on a
``DENSE`` trial rather than return a silently-wrong coarsening. The
known future consumer is **higher-rank Generalized Energy Condensation**
(within-coarse-group spectral moments :math:`n \ge 1`; Rahnema,
Douglass & Forget 2008) — a richer coarse basis than the rank-0 P0
indicator — deferred to `GitHub #275
<https://github.com/deOliveira-R/ORPHEUS/issues/275>`_. No
``LeastSquaresFrame`` type exists today; the name marks the seam, not a
shipped class.

.. note::

   Cross sections never need this dense seam at rank 0: a P0 (indicator)
   coarse cross section is the only rate-meaningful one, and its
   partition-of-unity Gram is row-sum-collapsible. The dense
   :math:`(MR)^{-1}M` solve becomes load-bearing only for a
   non-indicator coarse basis — exactly the GEC :math:`n \ge 1`
   moments — so the refusal is a forward-looking guard, not a
   present-day gap.


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

The **structural leg** — that Funk–Hecke makes the spherical harmonics
scattering's eigenbasis, so the frame is Galerkin — is worked in
:ref:`frame-spherical-harmonic-galerkin` above; the subsections below
read its consequence for *ownership*: the asymmetry that assigns the
frame to scattering, the literature corroboration, and the unifying
principle.

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
   * - Spatial homogenisation (forward / reaction-rate)
     - PetrovGalerkinFrame
     - Flux-weighted test basis, indicator trial basis
       (:meth:`Solution.homogenize
       <orpheus.sn.solution.Solution.homogenize>`)
     - **Live** (P3)
     - :ref:`sn-homogenization-petrov-galerkin-frame`; Hébert 2009 §13
   * - Energy condensation (forward / reaction-rate)
     - PetrovGalerkinFrame
     - Spectrum-weighted test basis, fractional-overlap trial basis
       (:meth:`Solution.condense
       <orpheus.sn.solution.Solution.condense>`)
     - **Live** (P5)
     - :ref:`sn-energy-condensation`; Hébert 2009 §6.2
   * - Homogenisation / condensation (eigenvalue-consistent)
     - PetrovGalerkinFrame
     - **Adjoint**-flux test basis (:math:`\varphi^*`), indicator /
       overlap trial basis
     - Pending (P6) — blocked on :math:`\varphi^*`
       (:ref:`frame-adjoint-weighted-seam`)
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
  just SN. The forward Petrov-Galerkin frames now carry their own
  rate-preservation **L0** gates (:mod:`tests.sn.test_homogenization` —
  the per-channel rate identity, the φV-vs-dV discriminator, and the
  Mode-11 routing sentinel; :ref:`sn-homogenization-verification`). The adjoint-weighted
  (:math:`\varphi^* \ne \varphi`) cross-Gram gates land with P6.


.. _frame-adjoint-weighted-seam:

The adjoint-weighted seam — documented, not built
=================================================

The forward Petrov-Galerkin frames that ship today
(:meth:`Solution.homogenize <orpheus.sn.solution.Solution.homogenize>`,
:meth:`Solution.condense <orpheus.sn.solution.Solution.condense>`)
weight the test functions by the **forward** flux
:math:`\chi_R = \varphi\,\mathbf{1}_R`. That is the
**Galerkin-degenerate** (:math:`\varphi^* = \varphi`) case of the
projection reactor physics ultimately wants. The general,
**eigenvalue-consistent** projection weights the test functions by the
**adjoint** flux :math:`\varphi^*`, preserving the bilinear functional

.. math::
   :label: sn-homogenization-adjoint-weighted

   \Sigma_R \;=\;
   \frac{\int_R \varphi^*\,\Sigma\,\varphi\;\mathrm{d}V}
        {\int_R \varphi^*\,\varphi\;\mathrm{d}V},

.. (vv-status rationale) Documents code that does not exist yet: the
   eigenvalue-consistent (adjoint-weighted, φ*≠φ) homogenisation projection is
   NOT built (campaign phase P6, blocked on the #276 adjoint flux φ*). Only the
   forward (φ*=φ) degenerate ships. The specification the P6 slice implements —
   the sibling bilinear identity :eq:`sn-homogenization-bilinear` carries the
   same documented-only status. Not a verified solver claim.
.. vv-status: sn-homogenization-adjoint-weighted documented

so that the multiplication factor :math:`\keff` stays stationary under
the homogenisation (by first-order perturbation theory :math:`\keff` is
stationary with respect to the adjoint-weighted residual). The full
derivation — why this is *irreducibly* Petrov-Galerkin
(:math:`M^* \ne R`: there is **no** metric in which test
:math:`= \varphi^*\mathbf{1}_R` equals trial
:math:`= \varphi\,\mathbf{1}_R` when :math:`\varphi^* \ne \varphi`) — is
:eq:`sn-homogenization-bilinear` in
:ref:`sn-homogenization-why-petrov-galerkin`,
and it applies *verbatim* on the energy axis for condensation
(:ref:`sn-energy-condensation`).

**It generalises a degenerate that already ships.** The forward
reaction-rate functional ORPHEUS already computes —
:class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`,
:math:`\int_R \langle\Sigma_x, \varphi\rangle\;\mathrm{d}V` — is the
:math:`\varphi^* = 1` degenerate of the bilinear
:math:`\langle\varphi^*, \Sigma\varphi\rangle`. The eigenvalue-consistent
lift is the **same** Petrov-Galerkin frame with the implicit
:math:`\varphi^*` (which the forward frame takes to be the forward flux
:math:`\varphi`) replaced by a *real* adjoint flux on the test basis — a
no-op change of the ``test_basis`` (:math:`\varphi \to \varphi^*`),
**not** a re-derivation. Writing the discipline on the frame **type** (an
explicit test basis) rather than on the measure (a flux-folded metric) is
precisely what buys this: the adjoint case is a one-line test-basis swap.

.. important:: **Status — theory documented, implementation NOT built
   (P6).** The adjoint-weighted projection is **documented theory only**.
   No adjoint flux :math:`\varphi^*` is wired into
   :meth:`Solution.homogenize <orpheus.sn.solution.Solution.homogenize>`
   or :meth:`Solution.condense <orpheus.sn.solution.Solution.condense>`
   today; both run the forward (:math:`\varphi^* = \varphi`) degenerate.
   The implementation (campaign phase **P6**) is **blocked** on a
   converged adjoint flux — the scalar moment :math:`\varphi^*` of the
   adjoint angular flux :math:`\psi^*` solving
   :math:`(L + C - S)^{\mathsf T}\psi^* = q^*` (the SN adjoint-transport
   campaign, `#276`; the discrete scattering adjoint
   :math:`S^{\mathsf T}` it builds on is :ref:`sn-scattering-adjoint`).
   That campaign is in flight but has not yet delivered a
   :math:`\varphi^*` field a homogenisation consumer can read. Until it
   does, this section is the **specification** the P6 slice implements:
   the bilinear identity :eq:`sn-homogenization-bilinear` carries
   ``.. vv-status: documented`` — the structure the forward case
   degenerates from, **not** a verified solver claim. Do **not** read the
   forward homogenisation's green rate gates as evidence for the
   adjoint-weighted case.


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
companion L0/foundation shape and predicate tests verify software
invariants (frame face spaces, the ``is_invertible`` / ``is_adjointable``
predicates) and are tagged accordingly. The **forward** Petrov-Galerkin frames now ship their own
**L0** numerical evidence — the per-channel rate-preservation identity,
the φV-vs-dV (flux- vs volume-weighting) discriminator, the simplex /
production-weight :math:`\chi` gates, and the Mode-11 routing sentinel —
in :mod:`tests.sn.test_homogenization`
(:ref:`sn-homogenization-verification`),
together with the condensation gates of :ref:`sn-energy-condensation`.
The adjoint-weighted (:math:`\varphi^* \ne \varphi`) cross-Gram evidence
lands with P6 (:ref:`frame-adjoint-weighted-seam`).


Implementation map
==================

* :class:`~orpheus.numerics.frame.FrameBase` — the abstract discrete
  frame: binds a :class:`~orpheus.numerics.basis.Basis` to a
  :class:`~orpheus.numerics.measure.DiscreteMeasure` and emits the
  ``analysis`` (:math:`M = T`) and ``reconstruction`` (:math:`R`)
  faces. Carries the discipline-free mechanics (table, the two
  spaces, the reconstruction face, the analysis-face wiring); the
  single mechanism for every choice-dependent change-of-basis. Also
  emits the **composed-operator verbs**
  (:ref:`frame-composed-verbs`):
  :meth:`conjugate <orpheus.numerics.frame.FrameBase.conjugate>`
  (:math:`R\circ A\circ M`),
  :meth:`reconstruct_after <orpheus.numerics.frame.FrameBase.reconstruct_after>`
  (:math:`R\circ A`), and
  :meth:`project <orpheus.numerics.frame.FrameBase.project>`
  (:math:`G^{-1}M`) with its :attr:`gram
  <orpheus.numerics.frame.FrameBase.gram>` cross-Gram probe.
* :class:`~orpheus.numerics.basis.GramStructure` — the trial basis's
  projection-validity declaration (``DIAGONAL`` / ``PARTITION_OF_UNITY``
  / ``DENSE``) that decides whether :meth:`project
  <orpheus.numerics.frame.FrameBase.project>` uses the row-sum probe or
  refuses (:ref:`frame-least-squares-discipline`).
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
