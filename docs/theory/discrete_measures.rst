.. _discrete-measures:

============================================
Discrete Measures and Quadrature Composition
============================================

.. note::

   **Stub page** — full narrative pending under Issue 18 of
   ``.claude/plans/sn_reshape.md``. This stub installs the
   load-bearing equation labels (:eq:`discrete-measure-integrate`,
   :eq:`discrete-measure-pushforward`,
   :eq:`bundle-measure-disintegration`) and the cross-references
   that downstream theory pages and tests rely on.

Key Facts
=========

- A **quadrature rule** is a discrete measure
  :math:`\mu = \sum_i w_i \, \delta_{x_i}` on a measurable space
  :math:`(\mathcal{X}, \Sigma)`. Lebesgue integration of a test
  function against :math:`\mu` is the familiar quadrature sum
  (:eq:`discrete-measure-integrate`).
- Four canonical composition operations are exposed by the project's
  primitive :class:`~orpheus.numerics.measure.DiscreteMeasure`:

  - **Tensor product** ``μ * ν`` — product measure
    :math:`\mu \otimes \nu` on
    :math:`\mathcal{X} \times \mathcal{Y}`.
  - **Direct sum** ``μ + ν`` — concatenation of atoms on a shared
    space; the foundation of composite (panelled) quadrature.
  - **Pushforward** ``μ.pushforward(φ)`` — image measure
    :math:`\varphi_* \mu` (:eq:`discrete-measure-pushforward`).
  - **Restriction** ``μ.restrict(E)`` — indicator-multiplication
    :math:`\mathbf{1}_E \cdot \mu`; supports half-range SN sweeps
    and bundle cuts.

- :class:`~orpheus.numerics.measure.BundleMeasure` carries the
  **disintegration** :eq:`bundle-measure-disintegration`: a base
  measure on :math:`\mathcal{B}` plus, for each base point, a fiber
  measure on :math:`\pi^{-1}(b)`. SN does not consume bundles in the
  Wave A campaign; MoC ray-bundle quadratures (Wave 2) will.
- The space tag is a runtime ``str`` — Python generics over a
  measurable-space type are not expressive enough without runtime
  overhead, so the project uses opaque string tags for sanity checks
  and documentation only.


Definitions
===========

A discrete (atomic) measure on a measurable space
:math:`(\mathcal{X}, \Sigma)` is a finite sum of weighted Dirac
deltas:

.. math::
   :label: discrete-measure-definition

   \mu = \sum_{i=1}^{N} w_i \, \delta_{x_i},
   \qquad x_i \in \mathcal{X},\quad w_i \in \mathbb{R}.

.. vv-status: discrete-measure-definition documented

The Lebesgue integral of an integrable test function :math:`f`
against :math:`\mu` collapses to a finite sum,

.. math::
   :label: discrete-measure-integrate

   \int_{\mathcal{X}} f \, d\mu \;=\; \sum_{i=1}^{N} w_i \, f(x_i),

which is the **quadrature** identity. Polynomial-exact rules
(Gauss-Legendre, Gauss-Chebyshev, …) are characterised by the
maximum degree :math:`p` such that
:eq:`discrete-measure-integrate` is exact for every polynomial of
degree :math:`\le p` — the field's
:attr:`~orpheus.numerics.measure.DiscreteMeasure.degree_of_exactness`
records this.

For a measurable map
:math:`\varphi : \mathcal{X} \to \mathcal{Y}`, the **pushforward**
or image measure :math:`\varphi_* \mu` satisfies the
change-of-variables identity

.. math::
   :label: discrete-measure-pushforward

   \int_{\mathcal{Y}} g \, d(\varphi_* \mu)
   \;=\; \int_{\mathcal{X}} (g \circ \varphi) \, d\mu
   \;=\; \sum_i w_i \, g(\varphi(x_i)).

.. (vv-status rationale) discrete-measure-pushforward: Verified
   transitively by the foundation test
   ``test_pushforward_invertible_map_change_of_variables`` —
   integrates :math:`x^2` against the φ-image of an equispaced
   measure on :math:`[0, 1]` under :math:`\varphi(x) = x^2` and
   confirms convergence to :math:`\int_0^1 x^4 \, dx = 1/5`,
   exactly the change-of-variables identity in this label. Tagged
   ``foundation`` rather than carrying ``verifies(...)`` because
   the identity is a software invariant of the discrete pushforward
   operation, not a physics-equation claim with an L0..L3 ladder
   slot.
.. vv-status: discrete-measure-pushforward documented

The pushforward in
:class:`~orpheus.numerics.measure.DiscreteMeasure` is
**φ-image semantics**: weights are preserved verbatim. If the
caller wants Radon-Nikodym semantics
:math:`g(y)\,dy = h(x)\,dx`, the Jacobian
:math:`|\det D\varphi|^{-1}` must be applied to the weights
explicitly.

For a fibered space :math:`\pi : \mathcal{X} \to \mathcal{B}`, the
**disintegration theorem** (Bourbaki 1969, Intégration VI §3) gives

.. math::
   :label: bundle-measure-disintegration

   \int_\mathcal{X} f \, d\mu
   \;=\; \int_\mathcal{B} \!\left[ \int_{\pi^{-1}(b)} f \, d\mu_b
   \right] d\mu_\mathcal{B}(b),

.. (vv-status rationale) bundle-measure-disintegration: Verified
   transitively by the foundation tests
   ``test_bundle_measure_separable_smoke`` (independent fiber, base
   :math:`\otimes` fiber agrees with explicit nested sum) and
   ``test_bundle_measure_varying_fiber`` (fiber depends on base
   point, dispatch path is exercised). MoC's ray-bundle quadratures
   in Wave 2 of the SN-reshape campaign will provide the physics-grade
   verification once the abstraction is consumed.
.. vv-status: bundle-measure-disintegration documented

implemented in
:class:`~orpheus.numerics.measure.BundleMeasure` as a base measure
plus a lazy fiber-measure factory.


Composition algebra — metadata propagation
==========================================

The four composition operations propagate the
:class:`~orpheus.numerics.measure.DiscreteMeasure` metadata as
follows. Where a field is marked **dropped**, the operation
generally cannot guarantee the field is still meaningful, and the
result carries ``None`` until the caller re-asserts it via
:meth:`~orpheus.numerics.measure.DiscreteMeasure.with_metadata`.

.. list-table::
   :header-rows: 1
   :widths: 18 18 16 16 16 16

   * - Operation
     - ``space``
     - ``n_points``
     - ``dim``
     - ``degree_of_exactness``
     - ``invariance_group``
   * - tensor product ``μ * ν``
     - ``f"{μ.space} × {ν.space}"``
     - :math:`N_\mu \cdot N_\nu`
     - :math:`d_\mu + d_\nu`
     - :math:`\min(p_\mu, p_\nu)` if both set, else dropped
     - dropped (Issue 3 may refine when both factors carry
       compatible groups)
   * - direct sum ``μ + ν``
     - shared (raises if mismatched)
     - :math:`N_\mu + N_\nu`
     - shared (raises if mismatched)
     - :math:`\min(p_\mu, p_\nu)` if both set, else dropped
     - dropped (concatenation generally breaks invariance)
   * - pushforward ``μ.pushforward(φ)``
     - explicit ``new_space=`` argument or ``f"φ_*({μ.space})"``
     - unchanged
     - inferred from ``φ``'s output shape
     - dropped (φ-image of a polynomial-exact rule is not
       polynomial-exact in general — undefined unless ``φ`` is
       affine)
     - dropped (φ generally does not preserve a group action)
   * - restriction ``μ.restrict(E)``
     - unchanged
     - ``mask.sum()``
     - unchanged
     - dropped (restricted rule is not Gauss on the cut domain)
     - dropped (caller may re-tag if :math:`E` is invariant)

Be honest about what is dropped: the field becoming ``None`` is the
correct behaviour, not a degradation. If a downstream consumer
needs the field, it must re-establish the claim and stamp it back
with :meth:`~orpheus.numerics.measure.DiscreteMeasure.with_metadata`.


1-D primitive constructors
==========================

The module ships three 1-D primitives:

- :func:`~orpheus.numerics.measure.gauss_legendre`
  ``(n)`` — Gauss-Legendre on :math:`[-1, 1]`,
  ``degree_of_exactness = 2n - 1``, weights sum to 2.
- :func:`~orpheus.numerics.measure.gauss_chebyshev`
  ``(n)`` — Gauss-Chebyshev (first kind) on :math:`[-1, 1]` with
  weight :math:`(1 - x^2)^{-1/2}`, ``degree_of_exactness = 2n - 1``
  (in the weighted sense), weights sum to :math:`\pi`.
- :func:`~orpheus.numerics.measure.equispaced`
  ``(a, b, n)`` — midpoint rule on :math:`[a, b]`,
  ``degree_of_exactness = 1``, weights sum to :math:`b - a`.

Higher-dimensional rules are built by composition:
``μ_S2 = gauss_legendre(n_mu) * equispaced(0.0, 2*np.pi, n_phi)``
gives a product quadrature on
:math:`\mu \in [-1, 1] \times \varphi \in [0, 2\pi)` —
the standard polar-azimuthal split for transport on the unit
sphere.


Symmetry groups for quadrature invariance
=========================================

A discrete measure :math:`\mu = \sum_i w_i \, \delta_{x_i}` is
**:math:`G`-invariant** under the action of a group :math:`G` iff
every :math:`g \in G` permutes the support points among themselves
in a weight-preserving way:

.. math::

   \forall g \in G,\;\;\;
   \exists \text{ a permutation } \pi_g \in S_N
   \text{ such that }\;
   g(x_i) = x_{\pi_g(i)}
   \;\text{ and }\;
   w_i = w_{\pi_g(i)}.

When this holds, integrating any :math:`G`-invariant integrand
:math:`f` (i.e. :math:`f \circ g = f` for all :math:`g \in G`) against
:math:`\mu` is automatically :math:`G`-symmetric in the result —
no spurious low-order asymmetry leaks from the quadrature rule into
moments where physics demands their absence.

Quadrature selection in ORPHEUS therefore reduces to a containment
check: given a geometry with natural symmetry group
:math:`G_{\text{geom}}` (slab :math:`\to SO(2) \times Z_2`,
sphere :math:`\to O(3)`, hexagonal lattice :math:`\to D_{6h}`, …)
and a candidate quadrature with invariance group
:math:`G_{\text{quad}}`, the rule is

.. math::
   :label: subgroup-of-o3-containment

   G_{\text{geom}} \subseteq G_{\text{quad}}
   \;\Longleftrightarrow\;
   \forall g \in G_{\text{geom}},\; g \in G_{\text{quad}}.

.. (vv-status rationale) subgroup-of-o3-containment: Verified
   transitively by the foundation tests in
   :file:`tests/numerics/test_symmetry.py` — every named relation in
   the static lattice (``Trivial ⊂ Z_2 ⊂ O_h ⊂ O(3)``,
   ``SO(2) ⊂ O(2) ⊂ O(3)``, the ``O_h \not\subset SO(3)`` improper-
   rotation test, the parameterised ``C_n ⊂ C_m \iff n | m`` rule,
   and reflexivity for every named entry) is asserted directly.
   Tagged ``foundation`` rather than carrying a verification ladder
   slot because the containment lattice is a software invariant
   (Hamermesh 1962 §2.5; Stiefel & Fässler 1979 Ch. 4) — there is
   no L0..L3 physics-equation claim to ladder against.
.. vv-status: subgroup-of-o3-containment documented

The "if and only if" is the standard set-theoretic definition of a
subgroup; the value here is the *application* — guaranteeing that
selecting a quadrature whose invariance group contains the geometry's
symmetry group is sufficient to preserve every symmetry the geometry
exhibits.

ORPHEUS's relevant sub-lattice of :math:`O(3)` is **finite and
small**. Slab, sphere, 1-D / 2-D Cartesian geometries combined with
the four shipped quadrature families (Gauss-Legendre, Lebedev,
level-symmetric, product) span exactly the eight named entries
encoded in :class:`~orpheus.numerics.symmetry.SubgroupOfO3` —
``Trivial``, ``Z_2``, ``SO(2)``, ``O(2)``, ``O_h``, ``I_h``,
``SO(3)``, ``O(3)`` — plus the parameterised families :math:`C_n`
and :math:`D_{nh}` reserved for forthcoming hex / triangular
lattices. The subgroup-containment relation is therefore captured
exhaustively by a static dict; no character-table or generator-based
machinery is implemented because none is yet needed.

Lebedev quadratures (Lebedev 1976) and level-symmetric :math:`S_N`
quadratures (Carlson & Lathrop 1968) are :math:`O_h`-invariant *by
construction*: their generating points are constrained to lie on
the orbits of the octahedral group through the choice of
free-parameter algebraic equations. The
:meth:`~orpheus.numerics.symmetry.SubgroupOfO3.is_invariant` check
in this module therefore *confirms* what the construction
guarantees rather than *discovering* invariance — its job is to
verify that the project's wrapper code preserves the property the
literature establishes, and to reject accidental loss (e.g., a
quadrature reshaped through a non-symmetric pushforward).

The non-trivial design choice for ``is_invariant`` is the
**fingerprint strategy**: rather than enumerating all 48 elements of
:math:`O_h` or all 120 elements of :math:`I_h`, the check applies a
generating set (8 sign flips × 6 coordinate permutations for
:math:`O_h`; 60 proper rotations + inversion for :math:`I_h`) and
verifies that every generator's image of the node array agrees with
the input under nearest-neighbour matching. The 1-D case
(:math:`SO(2)`-invariance of a measure on :math:`[-1, 1]`) is
trivial — there is no azimuthal coordinate to rotate — and the
:math:`Z_2` reflection :math:`x \to -x` is the only non-trivial
1-D check.

Tensor-product propagation of invariance is **deliberately dropped**
to ``None`` (see the metadata-propagation table above):
:math:`\mu \otimes \nu` is :math:`G_\mu \times G_\nu`-invariant, but
the diagonal subgroup of :math:`G_\mu \times G_\nu` is the natural
"join" only when both factors carry compatible group actions on a
shared space. Encoding the joint-group construction is left for a
later issue if and when a consumer needs it; today the field is
re-stamped on tensor-product results via
:meth:`~orpheus.numerics.measure.DiscreteMeasure.with_metadata`
when the caller knows the result is invariant.

Domain-specific adapters: AngularQuadrature
===========================================

The SN solver consumes angular quadratures through ~50
attribute-access sites (``quad.mu_x``, ``quad.weights``,
``quad.reflection_index('x')``, ``quad.spherical_harmonics(L)``, …)
in sweeps, BiCGSTAB operators, mesh constructors, and solvers.
Re-routing every site through
:meth:`DiscreteMeasure.integrate` would impose a per-site method-call
overhead on hot inner loops where the historical attribute access
costs nothing, and would needlessly couple the entire SN module to
the measure API.

The Issue 4 design (see ``.claude/plans/sn_reshape.md``) preserves
both worlds via a **bridge pattern**:

* The publication-grade quadrature *rules* live in
  :mod:`orpheus.numerics.quadrature` as free functions returning
  :class:`DiscreteMeasure` instances tagged with
  ``invariance_group`` and ``degree_of_exactness``:

  - :func:`~orpheus.numerics.quadrature.gauss_legendre_on_mu` —
    1-D rule on :math:`\mu \in [-1, 1]`,
    :math:`SO(2)`-invariant, ``degree_of_exactness = 2n - 1``.
  - :func:`~orpheus.numerics.quadrature.lebedev_sphere` —
    Lebedev rule on :math:`S^2`, :math:`O_h`-invariant,
    ``degree_of_exactness = order``.
  - :func:`~orpheus.numerics.quadrature.level_symmetric_sn` —
    Carlson-Lathrop level-symmetric :math:`S_N` rule on
    :math:`S^2`, :math:`O_h`-invariant.
  - :func:`~orpheus.numerics.quadrature.product_mu_phi` — the
    polar-times-azimuthal product rule used by the cylindrical SN
    sweep.

* The SN-side adapters at :mod:`orpheus.sn.quadrature` —
  :class:`~orpheus.sn.quadrature.GaussLegendre1D`,
  :class:`~orpheus.sn.quadrature.LebedevSphere`,
  :class:`~orpheus.sn.quadrature.LevelSymmetricSN`,
  :class:`~orpheus.sn.quadrature.ProductQuadrature` — wrap those
  rules and **cache numpy views** of the underlying measure's
  ``nodes`` / ``weights`` arrays as the legacy ``mu_x`` / ``mu_y``
  / ``mu_z`` / ``weights`` fields. The view aliasing is safe
  because :class:`DiscreteMeasure` is ``frozen=True``; the cached
  views never get out of sync with the measure.

The cached-view layer means the existing SN consumers see no API
change: ``quad.mu_x`` is still a numpy array, accessed at the same
:math:`O(1)` cost. The underlying :class:`DiscreteMeasure` is
exposed as ``adapter.measure`` for callers that want the
structurally-richer object — composability via
:meth:`~DiscreteMeasure.__mul__` /
:meth:`~DiscreteMeasure.pushforward` /
:meth:`~DiscreteMeasure.restrict`, plus the ``invariance_group``
that the upcoming ``select_quadrature`` registry (Issue 5) will
consume.

The bit-identical contract
--------------------------

The four legacy classes have ~50 active call sites and 11 frozen
regression snapshots at ``tests/sn/regression/snapshots/`` (one per
canonical SN configuration: slab / sphere / cylinder × 1G / 2G ×
homogeneous / multi-region × DD / aniso). The bridge refactor
preserves both:

* **Bit-identical foundation tests** at
  ``tests/numerics/test_rules_*.py`` use :func:`numpy.array_equal`
  (not :func:`numpy.allclose`) to verify that each rule function's
  nodes/weights are exact-bit equal to the legacy adapter's
  output.
* **Regression snapshots** rerun every SN configuration end-to-end
  through ``DDSweep`` / ``BiCGSTAB`` / power iteration and assert
  match against the frozen ``.npz`` baselines at the iterative
  solver's convergence floor (``rtol=1e-12``, ``atol=1e-13``).

Both layers are required: the bit-identical match at the rule level
guarantees no node-order or floating-point drift entered the
discretisation pipeline; the snapshot match guarantees that the
pipeline's downstream consumers (sweeps, operators, solvers) still
produce the same answer they always did.

Forward references
==================

This page is consumed progressively. The current stub installs the
labels and the composition algebra; the following pieces will land
in subsequent issues of the SN reshape campaign:

- **Issue 5** (``module:numerics``) will install a
  ``select_quadrature`` registry that consumes
  :meth:`~orpheus.numerics.symmetry.SubgroupOfO3.is_subgroup_of` to
  drive geometry → quadrature selection. The named entries and the
  containment lattice in :eq:`subgroup-of-o3-containment` are the
  data the registry reads.
- **Wave 2 (MoC migration)** will consume
  :class:`~orpheus.numerics.measure.BundleMeasure` for ray-bundle
  quadratures whose fibers (the parallel rays through the geometry)
  vary with the base direction. The bundle abstraction is shipped
  in Phase 0 because retrofitting it later would force a rewrite
  of MoC's quadrature handling.

The narrative in :ref:`theory-reference-solvers` and
:ref:`theory-discrete-ordinates` will progressively replace any
quadrature-specific wording with cross-references to
:eq:`discrete-measure-integrate` and the composition algebra
above.


References
==========

* Stoer, J. and Bulirsch, R. (2002). *Introduction to Numerical
  Analysis*, 3rd ed. Springer. §3.6 (Gauss quadrature);
  Theorem 3.6.20 (polynomial exactness).
* Trefethen, L.N. (2008). "Is Gauss quadrature better than
  Clenshaw-Curtis?" *SIAM Review* **50**, 67-87. The measure-theoretic
  perspective on spectral integration rules.
* Bourbaki, N. (1969). *Éléments de mathématique. Intégration*,
  Chapter VI §3 (disintegration of measures).
* Hamermesh, M. (1962). *Group Theory and Its Application to Physical
  Problems*. Addison-Wesley. §2.5 (finite point groups), §9.4
  (octahedral and icosahedral groups).
* Stiefel, E. and Fässler, A. (1979). *Group Theoretical Methods and
  Their Applications*. Birkhäuser. Chapters 4-5 (crystallographic
  point groups, invariant theory of finite groups).
* Lebedev, V.I. (1976). "Quadratures on a sphere." *USSR Computational
  Mathematics and Mathematical Physics* **16**, no. 2, 10-24. The
  octahedral-invariant construction validated by
  :meth:`~orpheus.numerics.symmetry.SubgroupOfO3.is_invariant`.
* Carlson, B.G. and Lathrop, K.D. (1968). "Transport theory: the
  method of discrete ordinates." In *Computing Methods in Reactor
  Physics*, Greenspan, Kelber, Okrent, eds., Gordon & Breach.
  Level-symmetric :math:`S_N` quadratures, also :math:`O_h`-invariant
  by construction.
