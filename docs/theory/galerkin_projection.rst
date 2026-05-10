.. _galerkin-projection:

============================================================
Galerkin Projection — discipline, primitives, and consumers
============================================================

Every method in ORPHEUS that transitions a function between a
**fine** representation and a **coarse** one — discrete-ordinate
angular flux to spherical-harmonic moments, fine-energy fluxes to
broad-group cross-sections, regional flux to homogenised cross-
sections — does so via a **(reconstruction, projection)** pair of
linear operators :math:`(R, \Pi)`. This page is the canonical
home for that pair, the **discipline** (Galerkin vs Petrov-Galerkin)
that distinguishes its variants, and the **cross-method consumer
catalog** of where the same primitives are used in the codebase.

The user's architectural rule that motivated this page:

  *"Galerkin projection is generally useful — same machinery is
  used in cross-section energy condensation. Catalog everything
  that needs to be implemented and where."*

Lifting the Galerkin pair into :mod:`orpheus.numerics.projection`
(Wave 0, step C0.5 of the SN performance plan) puts one set of
primitives in front of every consumer, instead of a separate
implementation per method.

.. contents::
   :local:
   :depth: 2


Key Facts
=========

- A **projection** :math:`\Pi : V \to W` is a linear operator from a
  fine space :math:`V` to a coarse coefficient space :math:`W`. Its
  paired **reconstruction** :math:`R : W \to V` is the linear
  operator that lifts coefficients back to :math:`V`. Together
  :math:`(R, \Pi)` define a Galerkin discretisation of any
  :math:`A : V \to V` as
  :math:`A_h = \Pi A R : W \to W` (Brenner & Scott 2008, §3.4).

- Two **disciplines** distinguish how :math:`(R, \Pi)` is
  constructed:

  * **Galerkin** — test space equals trial space (the inner product
    used to define :math:`\Pi` is taken in :math:`V` itself). The
    canonical case for :math:`L^2`-orthogonal moment projection.
    Invariant: :math:`\Pi^* = R` under the :math:`V`-inner product.
  * **Petrov-Galerkin** — test space differs from trial space (the
    inner product on :math:`V` and the inner product used to define
    :math:`\Pi` are different). The canonical case for energy
    condensation (within-group spectrum weighting) and spatial
    homogenisation (flux-volume weighting). Sibling discipline of
    Galerkin.

- The naming hierarchy in :mod:`orpheus.numerics.projection` is
  **deliberately three-level** so the type itself signals which
  discipline a concrete projection follows:
  :class:`~orpheus.numerics.projection.ProjectionOperator` (most
  general) → :class:`~orpheus.numerics.projection.GalerkinProjection`
  / :class:`~orpheus.numerics.projection.PetrovGalerkinProjection`
  (discipline) →
  :class:`~orpheus.numerics.projection.HarmonicMomentProjection`
  (concrete basis).

- The single concrete pair shipping today is
  :class:`HarmonicMomentProjection` /
  :class:`HarmonicMomentReconstruction` for real spherical harmonics
  on a :math:`S^2` cubature. Future Petrov-Galerkin pairs land with
  energy condensation (§17 of Grand Report v3) and spatial
  homogenisation (§18). The ABCs are in place so the symmetry of
  the two disciplines is visible to future readers.

- Every concrete pair satisfies the **idempotency-on-coefficients**
  invariant on a sufficiently-exact quadrature:

  .. math::

     \Pi \, R \;=\; c_{V}\,I_{W},

  where :math:`c_V` is a scalar that depends on the inner-product
  convention of :math:`V`. For the SN spherical-harmonic Galerkin
  pair on a Lebedev quadrature, :math:`c_V = 4\pi` — this is the
  **L1 idempotency** identity :eq:`pi-r-equals-4pi-i` verified at
  multiple :math:`L` against multiple Lebedev orders.


The (reconstruction, projection) pair — anatomy
================================================

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

where :math:`c_V` is the inner-product normalisation constant
(:math:`c_V = 1` in the fully-orthonormal case;
:math:`c_V = 4\pi` for the no-prefactor real spherical harmonics).

.. vv-status: galerkin-pair documented

The pair is fully determined by:

1. The **fine space** :math:`V` and its inner product
   :math:`\langle \cdot, \cdot \rangle_V`. For SN angular flux,
   :math:`V = L^2(S^2)` and the inner product is the W-weighted
   discrete sum :math:`\langle f, g \rangle_W = \sum_n w_n f_n g_n`
   on the angular cubature.
2. The **basis** of :math:`W`. For SN scattering, the basis is the
   real spherical harmonics :math:`\{Y_\ell^m\}_{\ell \le L}`.
3. The **discipline**. Galerkin: the test space is the same basis.
   Petrov-Galerkin: the test space is different (e.g. weighted by
   the within-group spectrum :math:`\chi_g(\phi_g)`).

Once these three are chosen, :math:`\Pi` and :math:`R` are
uniquely determined up to the :math:`c_V` normalisation.


The Galerkin discipline
=======================

The Galerkin discipline is characterised by **test space equals
trial space**. The defining identity is

.. math::
   :label: galerkin-self-adjoint

   \Pi^* \;=\; R
   \quad \text{(under the } V \text{ inner product)}.

.. vv-status: galerkin-self-adjoint documented

i.e. the projection's Hilbert adjoint is its reconstruction. This is
why a single basis :math:`\{e_k\}` produces both :math:`\Pi` and
:math:`R`:

.. math::
   :label: galerkin-construction

   (\Pi f)_k &\;=\; \langle e_k, f \rangle_V, \\
   R \, c     &\;=\; \sum_k c_k\,e_k.

.. vv-status: galerkin-construction documented

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
reconstruction's :math:`(2\ell+1)` factor (the addition-theorem
weight) yields :math:`\Pi R = 4\pi I` — the L1 identity that the
test
``tests/numerics/test_projection_operators.py::
TestGalerkinIdempotencyOnLebedev::test_pi_R_is_identity_on_band_limited``
verifies at :math:`L = 2,\,3,\,4`.

.. note::

   The identity :math:`\Pi R = 4\pi I` is **not** identity-on-the-
   nose because the no-prefactor convention pushes the
   :math:`4\pi/(2\ell+1)` factor onto the orthogonality. A
   :class:`HarmonicMomentProjection` with a strict :math:`\Pi R =
   I` invariant could be built by dividing the projection weights by
   :math:`4\pi`, but the project chose to absorb the factor at the
   reconstruction (the :math:`(2\ell+1)` weight) so the addition
   theorem reads cleanly. See :ref:`spherical-harmonics` for the
   convention rationale.


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
     - ``HarmonicMomentProjection``
       /
       ``HarmonicMomentReconstruction``
     - **Live** (this work, Wave 1)
     - §9 (Grand Report v3 line 1230)
   * - PN solver
     - Galerkin
     - Same
       ``HarmonicMomentProjection``
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
     - Same
       ``HarmonicMomentProjection``
       used as a variance-reduction estimator
     - Pending
     - Lewis & Miller 1993 §10
   * - SN sensitivity
     - Galerkin (adjoint)
     - Same pair, applied to the adjoint flux
     - Pending
     - Cacuci 2003

The two architectural payoffs:

* **One implementation per discipline**, not one per consumer.
  ``HarmonicMomentProjection`` is the same class whether SN uses
  it for scattering or PN uses it for streaming — the difference
  is which axis it's wrapped onto via the
  :class:`~orpheus.numerics.operator.TensorProductOperator`
  algebra (see :ref:`operator-algebra` and the tensor-product
  section there).
* **One V&V chain per discipline**. The Galerkin idempotency tests
  in :file:`tests/numerics/test_projection_operators.py` cover
  every Galerkin consumer, not just SN. Energy condensation will
  inherit the Petrov-Galerkin-specific tests when its concrete
  subclass lands.


The naming hierarchy — type signal vs documentation
===================================================

The user's pedantic naming rule: **a reader of a type name should
know its properties without reading the docstring.** This drives
the three-level inheritance:

.. code-block:: python

   class ProjectionOperator(LinearOperatorMixin, ABC):
       """Most general projection ABC."""

   class GalerkinProjection(ProjectionOperator, ABC):
       """Galerkin discipline: test space = trial space."""

   class PetrovGalerkinProjection(ProjectionOperator, ABC):
       """Petrov-Galerkin discipline: test space ≠ trial space."""

   @dataclass(frozen=True)
   class HarmonicMomentProjection(GalerkinProjection):
       """Concrete: Galerkin on real spherical harmonics."""

A reader of ``class HarmonicMomentProjection(GalerkinProjection)``
immediately knows:

1. It's a **projection** (inherits from
   :class:`ProjectionOperator`).
2. It follows **Galerkin** discipline (inherits from
   :class:`GalerkinProjection`, which adds the
   :math:`\Pi^* = R` invariant).
3. The basis is **real spherical harmonics on the angular axis**
   (the class name).

No docstring needed for those facts. Compare to the alternative
(rejected) hierarchy:

.. code-block:: python

   # REJECTED: collapses Galerkin/Petrov-Galerkin distinction
   class HarmonicMomentProjection(ProjectionOperator):
       """Galerkin projection on real spherical harmonics.

       (must be read to know the discipline)
       """

The rejected form makes the discipline a docstring claim instead of
a type claim. When a future reader reaches the energy-condensation
implementation and sees ``class EnergyCondensation(ProjectionOperator)``
they have no way to know — without reading the docstring — whether
it satisfies the Galerkin :math:`\Pi^* = R` invariant or not. The
type hierarchy answers this without reading prose.

The same rule drives the rejection of:

* ``HarmonicReconstruction`` → ``HarmonicMomentReconstruction``
  (parity with ``HarmonicMomentProjection`` — the type signals
  that the same moment basis is at play).
* ``Projection`` → ``ProjectionOperator`` (parity with
  ``LinearOperator`` — the ``Operator`` suffix signals "carries
  the operator algebra").


Numerical evidence
==================

Two L1-tagged test classes in
:file:`tests/numerics/test_projection_operators.py` verify the
Galerkin discipline's invariants on the
``HarmonicMomentProjection`` / ``HarmonicMomentReconstruction``
pair:

1. **Idempotency**:
   :math:`\Pi R c = 4\pi c` on band-limited
   coefficient input, verified at :math:`L = 2,\,3,\,4` against
   Lebedev orders :math:`7,\,13,\,17`. See
   :eq:`pi-r-equals-4pi-i` in :ref:`spherical-harmonics`.
2. **Adjoint pairing**:
   :math:`\langle \Pi \psi, c \rangle =
    \langle \psi, R_{\text{no-factor}} c \rangle_W`,
   verified to ``rtol=1e-12`` on a Lebedev order-13 grid at
   :math:`L = 3`.

The tests are tagged ``@pytest.mark.l1`` because they verify
mathematical identities of the operator algebra (V&V level L1 —
equation verification by analytical reference). The L0-tagged
shape and capability tests (``TestHarmonicMomentProjectionShapes``,
``TestCapabilities``) verify software invariants and are tagged
``@pytest.mark.l0``.


Implementation map
==================

* :class:`orpheus.numerics.projection.ProjectionOperator` — the
  most-general ABC.
* :class:`orpheus.numerics.projection.GalerkinProjection` — the
  Galerkin discipline ABC, carrying the
  :meth:`assert_galerkin_idempotency` invariant test.
* :class:`orpheus.numerics.projection.PetrovGalerkinProjection`
  — the Petrov-Galerkin discipline ABC; sibling.
* :class:`orpheus.numerics.projection.HarmonicMomentProjection` —
  the concrete spherical-harmonic Galerkin projection.
* :class:`orpheus.numerics.projection.HarmonicMomentReconstruction`
  — the paired reconstruction (the addition-theorem
  :math:`(2\ell+1)`-weighted expansion).
* :func:`orpheus.numerics.spherical_harmonics.evaluate_real_sh` —
  the underlying :math:`Y_\ell^m` evaluator (see
  :ref:`spherical-harmonics`).

The full-space projector — the operator that projects the SN
:math:`(N, n_x, n_y, n_g)` angular flux onto the
:math:`(L+1, 2L+1, n_x, n_y, n_g)` moment field — is built as a
**tensor product** of the angular-axis :math:`\Pi` and identity
operators on the spatial / energy axes:

.. code-block:: python

   from orpheus.numerics.operator import IdentityOperator
   from orpheus.numerics.projection import HarmonicMomentProjection

   M = HarmonicMomentProjection.from_measure(quad.measure, L=L)
   M_full = M & IdentityOperator() & IdentityOperator() & IdentityOperator()

The ``&`` dunder constructs the
:class:`~orpheus.numerics.operator.TensorProductOperator`. See
:ref:`operator-algebra` and the **Tensor product algebra** section
there for the relationship between this operator-algebra type and
the underlying numpy primitives (``np.einsum``, ``np.tensordot``,
``np.kron``).


History — why three classes, not one
=====================================

The Wave 0 plan originally drafted a single class
``HarmonicMomentProjection`` with no ABC. Two corrections in the
naming-audit pass:

1. **The discipline must be visible in the type, not the docstring.**
   Without the ``GalerkinProjection`` ABC layer, a concrete
   ``EnergyCondensation`` Petrov-Galerkin class would be
   indistinguishable at the type level from a Galerkin pair — the
   invariants ``Π R = I`` (Petrov-Galerkin) and ``Π R = c_V I``
   with :math:`\Pi^* = R` (Galerkin) are different, but a single
   ``ProjectionOperator`` ABC cannot encode the difference.

2. **The future reader will see only the type signature.** When
   energy condensation lands, a maintainer reading
   ``class EnergyCondensation(PetrovGalerkinProjection)`` knows
   immediately that :math:`\Pi^* \ne R` without reading docs —
   they can write code (e.g. an adjoint sensitivity calculation)
   that respects the discipline. With a flat hierarchy they would
   have to grep the docstring.

The cost is one extra ABC layer; the payoff is a self-documenting
type system that propagates discipline across the inheritance graph
mechanically.


References
==========

* Brenner, S. C. and Scott, L. R. (2008). *The Mathematical Theory
  of Finite Element Methods*, 3rd ed. Springer. §3.4 (Galerkin /
  Petrov-Galerkin general framework).
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
