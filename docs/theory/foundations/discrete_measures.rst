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
    space; the foundation of composite (panelled) :term:`quadrature`.
  - **Pushforward** ``μ.pushforward(φ)`` — image measure
    :math:`\varphi_* \mu` (:eq:`discrete-measure-pushforward`).
  - **Restriction** ``μ.restrict(E)`` — indicator-multiplication
    :math:`\mathbf{1}_E \cdot \mu`; supports half-range SN :term:`sweeps <sweep>`
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


.. _discrete-measure-partition:

Partition by labelling predicate
================================

The fifth composition operation (added in Wave 0 of the SN
performance plan, step C0.4) realises the **inverse** of direct sum:
given a measure :math:`\mu` and a labelling map
:math:`\ell : \mathcal{X} \to L` on the support points, recover the
disjoint decomposition

.. (vv-status rationale) Verified by
   ``tests/numerics/test_measure_partition.py`` — disjoint-union
   coverage, weight conservation, generic-S² octant predicate on
   Lebedev orders 5/9/17, and the round-trip identity
   :math:`\mu = \bigoplus_\lambda \mu_\lambda`.
.. vv-status: discrete-measure-partition documented

.. math::
   :label: discrete-measure-partition

   \mu \;=\; \bigoplus_{\lambda \in L} \mu_\lambda,
   \qquad
   \mu_\lambda
   \;=\; \sum_{i : \ell(x_i) = \lambda} w_i \, \delta_{x_i}.

The decomposition is **disjoint** by construction (every node
appears in exactly one partition entry) and **preserves total
mass** (sum of partition weights equals total mass of :math:`\mu`).
The
:meth:`~orpheus.numerics.measure.DiscreteMeasure.partition_by`
method returns a tuple of
:class:`~orpheus.numerics.measure.DiscreteMeasurePartition` entries,
each carrying the partition label, the indices into the parent
measure, and the restricted measure :math:`\mu_\lambda` on the same
space as :math:`\mu`.

Round-trip identity
-------------------

For any partition predicate :math:`\ell`, the direct sum of
partition measures equals the parent modulo ordering:

.. math::
   :label: partition-round-trip

   \mu \;=\; \mu_{\lambda_1} \oplus \mu_{\lambda_2}
              \oplus \cdots \oplus \mu_{\lambda_K}.

.. vv-status: partition-round-trip documented

This is the **inverse-of-direct-sum** identity, verified by
``tests/numerics/test_measure_partition.py`` on a Lebedev grid
under the octant-sign predicate. The ``invariance_group`` and
``degree_of_exactness`` fields are dropped on each partition entry
— a partition typically breaks any global invariance, and the
polynomial-exactness of the parent rule does not survive
restriction.

Octant partition — the canonical use
------------------------------------

The motivating use is the **angular octant partition** of a
spherical cubature: with the predicate

.. math::
   :label: octant-sign-predicate

   \ell(\hat\Omega)
   \;=\; \bigl(\mathrm{sign}\,\mu_x,\;
                \mathrm{sign}\,\mu_y,\;
                \mathrm{sign}\,\mu_z\bigr),

.. vv-status: octant-sign-predicate documented

the partition produces the eight octants of :math:`S^2` (or four
in 2-D, where :math:`\mu_z = 0` is degenerate and contributes a
separate ``sign=0`` entry). The
:class:`~orpheus.numerics.quadrature.Quadrature` class exposes
``octants`` as a cached property delegating to ``partition_by`` —
the SN sweep then iterates **octants** (4 in 2-D — structural)
and **anti-diagonals per octant** (sweep-DAG topology —
structural), with the :term:`ordinate` axis (:math:`N_{\text{oct}}`)
**internal** to every :meth:`apply` call within an octant.

This is the operator-algebra-level realisation of the SN
performance plan's central principle: *the metric should know its
iterative structure*. The Python ``for n in range(N)`` over
ordinates that exists in the legacy 2-D wavefront sweep is
replaced by a tensor-product action vectorised across the
octant-restricted ordinate axis.

Cross-method consumers
----------------------

The partition primitive is **infrastructure** — it is consumed by
multiple solvers wherever a measure is naturally split by a
discrete labelling predicate (per :ref:`tensorial-framing`):

.. list-table:: Cross-method consumers of partition_by
   :header-rows: 1
   :widths: 24 38 38

   * - Solver
     - Predicate
     - Source
   * - **SN 2-D sweep**
     - Octant-sign predicate
       :eq:`octant-sign-predicate` — partitions the angular
       cubature into the 4 (2-D) or 8 (3-D) sweep octants.
     - Wave 2 of SN performance plan (closed Issue #4); the SN
       sweep then iterates octants and dispatches each octant to a
       per-octant
       :class:`~orpheus.sn.loss_representation.sweep_graph.SweepDependencyGraph` —
       see :ref:`sweep-octant-dependency-graph` for the full
       architecture.
   * - **MoC track-bundle direction grouping**
     - Polar-angle bin predicate — partitions tracks by polar
       angle :math:`\theta`.
     - Future MoC consumer.
   * - **MC boundary-current scoring by hemisphere**
     - Sign-of-outward-normal predicate — partitions track-
       boundary intersections into incoming / outgoing.
     - Future MC consumer.
   * - **SN boundary realiser** (§16A.5)
     - Per-boundary-face incoming-ordinate predicate —
       partitions ordinates by which face they enter.
     - Future SN boundary-condition consumer.

Test invariants
---------------

The partition primitive carries three :math:`\mathrm{L0}` invariants
verified by
``tests/numerics/test_measure_partition.py``:

1. **Disjoint union**: every parent index appears in exactly one
   partition entry's ``indices`` array.
2. **Weight conservation**: the sum of restricted measure weights
   over all partitions equals the parent's total mass.
3. **Per-partition well-formedness**: every partition's
   ``measure`` carries the same ``space`` tag as the parent and
   has consistent ``nodes.shape[0] == weights.shape[0]``.

The tests run on Lebedev orders 5, 9, and 17 under the octant-sign
predicate and on a randomised generic predicate.


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

Higher-dimensional rules are built by composition. The tensor product
``gauss_legendre(n_mu) * equispaced(0.0, 2*np.pi, n_phi)`` gives a rule
on the **parameter rectangle**
:math:`(\mu, \varphi) \in [-1,1] \times [0, 2\pi)` — the standard
polar-azimuthal split.

.. warning::

   That product is **not yet a rule on** :math:`S^2`, and this page
   claimed it was until 2026-08-02. `[M]` It carries ``nodes`` of shape
   ``(n, 2)`` on support ``"[-1,1] × [0,6.28…]"``, while the production
   sphere rule
   :func:`~orpheus.numerics.quadrature.rules_product.product_mu_phi`
   carries ``(n, 3)`` on ``"S^2"``. The missing step is the
   **pushforward** onto direction cosines,

   .. math::

      (\mu, \varphi) \;\longmapsto\;
      \bigl(\sqrt{1-\mu^2}\cos\varphi,\;
            \sqrt{1-\mu^2}\sin\varphi,\; \mu\bigr),

   which is a change of *space*, not merely of coordinates — and it is
   exactly the step under which
   :attr:`~orpheus.numerics.measure.DiscreteMeasure.degree_of_exactness`
   is dropped, because polynomial exactness in :math:`(\mu, \varphi)` is
   not polynomial exactness in the direction cosines. ``product_mu_phi``
   performs the map itself rather than composing these two primitives;
   the two differ further in the azimuthal convention
   (``equispaced`` places **midpoints**, the production rule
   **left-endpoints**).


Symmetry groups for quadrature invariance
=========================================

A discrete measure :math:`\mu = \sum_i w_i \, \delta_{x_i}` is
**:math:`G`-invariant** under the action of a group :math:`G` iff
every :math:`g \in G` permutes the support points among themselves
in a weight-preserving way:

.. math::
   :label: discrete-measure-g-invariance

   \forall g \in G,\;\;\;
   \exists \text{ a permutation } \pi_g \in S_N
   \text{ such that }\;
   g(x_i) = x_{\pi_g(i)}
   \;\text{ and }\;
   w_i = w_{\pi_g(i)}.

.. (vv-status rationale) discrete-measure-g-invariance: Notation definition —
   G-invariance of a discrete measure (every g permutes the support points in a
   weight-preserving way). The definition the quadrature-selection containment
   check :eq:`subgroup-of-o3-containment` builds on; its concrete instances are
   the foundation gates in :file:`tests/numerics/test_symmetry.py`
   (``test_lebedev_is_octahedral_invariant`` positive +
   ``test_lebedev_is_NOT_icosahedral_invariant`` negative, via
   ``SubgroupOfO3.is_invariant``). A definitional identity, not a solver claim.
.. vv-status: discrete-measure-g-invariance documented

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
level-symmetric, product) span the eight named entries encoded in
:class:`~orpheus.numerics.symmetry.SubgroupOfO3` — ``Trivial``,
``Z_2``, ``SO(2)``, :math:`D_{\infty h}`, ``O_h``, ``I_h``,
``SO(3)``, ``O(3)`` — plus the parameterised families :math:`C_n`
and :math:`D_{nh}`.

.. note::

   Three statements in this paragraph were falsified by later work and
   are corrected above; they are recorded because each was load-bearing
   when written.

   * It named :math:`O(2)`. That entry is **realized** as axial
     rotations plus :math:`\sigma_h`, which is :math:`C_{\infty h}`,
     not :math:`O(2)` — and under either reading the asserted relation
     :math:`D_{nh} \subseteq O(2)` is false, because :math:`D_{nh}`
     carries :math:`C_2` axes lying *in* the plane. Renamed to
     :math:`D_{\infty h}`, which does contain every :math:`D_{nh}` and
     is the group a cylinder actually carries.
   * It called :math:`C_n` / :math:`D_{nh}` "reserved for forthcoming
     hex / triangular lattices". They are in **active use**: the
     product rule on :math:`S^2` is tagged :math:`D_{n_\varphi h}`,
     computed from its nodes. (It previously advertised :math:`SO(2)`
     — a claim **no** finite point set on :math:`S^2` can satisfy.)
   * It said containment "is captured exhaustively by a static dict; no
     character-table or generator-based machinery is implemented because
     none is yet needed." Both halves are now false. Containment between
     two FINITE groups is decided by **computed matrix containment** on
     their realized operator sets; the static table is consulted only
     when one side is continuous. The lattice itself is a **computed**
     Hasse diagram of maximal-subgroup relations, walked downward from
     high symmetry to find the symmetry a node set actually has — the
     crystallographic construction. A *declared* invariance group is a
     claim with no construction behind it, and such claims shipped false
     here more than once; a computed one cannot lie about the object it
     was computed from.

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

Quadrature selection algorithm
==============================

Issue 5 of the SN reshape campaign installs a tagged registry on top of
the four quadrature-rule primitives, plus a selector that picks the
right :class:`~orpheus.numerics.measure.DiscreteMeasure` for a
geometry-and-target-degree pair. The selector lives at
:mod:`orpheus.numerics.quadrature.registry`.

Selection is fundamentally a **four-stage filter** in priority order:

1. **G compatibility (symmetry).** The rule's invariance group must
   contain every symmetry of the geometry —
   :math:`G_{\text{geom}} \subseteq G_Q`. A quadrature with **less**
   symmetry than the geometry imprints spurious low-order asymmetry
   on a symmetric problem (Lebedev 1976, §1; Stiefel & Fässler 1979
   Ch. 5). A quadrature with **more** symmetry is harmless — the
   extra symmetry is unused, never violated. The lattice
   :eq:`subgroup-of-o3-containment` drives this filter; see
   :class:`~orpheus.numerics.symmetry.SubgroupOfO3` for the
   containment table.

2. **V compatibility (polynomial exactness, Galerkin sense).** The
   rule's degree of exactness must reach the target: :math:`\deg(Q)
   \ge d`. Each rule's degree is parameter-dependent — :math:`2n - 1`
   for Gauss-Legendre, :math:`N - 1` for level-symmetric :math:`S_N`
   (conservative), ``order`` for Lebedev, :math:`\min(2 n_\mu - 1,
   n_\phi - 1)` for the product rule. The selector inverts each rule's
   formula and picks the smallest parameter set meeting the target.
   Lebedev's gap structure (no rules at orders 33, 37, 39, 43, 45, 49)
   is handled by rounding up to the next tabulated order; if the target
   exceeds the table's top end (47 in scipy's tabulation), the rule is
   rejected at this stage with a clear message.

3. **Structural compatibility.** The consumer can request boolean
   flags that the rule must satisfy:

   - ``positive_weights`` — physics consumers reject negative-weight
     rules (Newton-Cotes :math:`n \ge 9`) because they amplify
     quadrature noise.
   - ``axis_aligned`` — Cartesian-friendly: nodes lie on coordinate
     axes (or include axis-aligned orbits, in the Lebedev case).
   - ``level_structured`` — exposes per-:math:`\mu` polar levels for
     the cylindrical SN sweep's azimuthal redistribution coefficients
     (Bailey et al. 2009 Eq. 50).
   - ``half_range_clean`` — ``measure.restrict(lambda x: x > 0)``
     gives a valid quadrature without re-normalisation.

   Only flags passed with value ``True`` constrain the search; ``False``
   and missing keys are interpreted as "don't care."

4. **Minimum points (cost).** Among candidates passing 1-3, pick the
   smallest :math:`N`. Each ordinate is a sweep direction in the SN
   solver, and sweeps dominate the runtime — the cheapest valid rule
   wins.

Formally, the selection criterion is

.. math::
   :label: quadrature-selection-criterion

   Q^{\star} \;=\; \arg\min\Bigl\{\, n(Q) \;:\;\;
   G_{\text{geom}} \subseteq G_Q
   \;\wedge\; \deg(Q) \ge d
   \;\wedge\; F_{\text{req}} \subseteq F_Q
   \,\Bigr\},

.. (vv-status rationale) quadrature-selection-criterion: Verified
   transitively by the foundation tests in
   :file:`tests/numerics/test_registry.py` — every stage of the
   four-stage filter has a happy-path test
   (``test_select_slab_returns_gauss_legendre``,
   ``test_select_sphere_returns_gauss_legendre``,
   ``test_select_cylinder_with_level_structured_returns_product``,
   ``test_select_cartesian2d_prefers_lebedev_over_ls_sn``) and a
   negative path (``test_no_rule_fits_raises_with_log``,
   ``test_truly_incompatible_flags_raises``). Tagged ``foundation``
   rather than carrying a verification ladder slot because the
   selection chain is a software invariant (the predicate
   :math:`\mathcal{D}_Q = S^2/G^0_{\text{geom}} \wedge
   \Gamma_{\text{geom}} \subseteq \operatorname{Sym}(Q) \wedge
   \deg(Q) \ge d \wedge
   F_{\text{req}} \subseteq F_Q`), not a physics-equation claim with
   an L0..L3 ladder slot — the ladder lives on the rules themselves
   (``test_rules_*.py``).
.. vv-status: quadrature-selection-criterion documented

where :math:`n(Q)` is the number of nodes,
:math:`\mathcal{D}_Q` is the domain the rule's nodes live on,
:math:`\operatorname{Sym}(Q) \subseteq O(3)` is the group the rule's
nodes are **computed** to be invariant under,
:math:`\deg(Q)` is the polynomial-exactness degree, and
:math:`F_Q \subseteq \{\text{positive\_weights}, \text{axis\_aligned},
\text{level\_structured}, \text{half\_range\_clean}\}` is the rule's
structural-flag set.

Geometry → angular-symmetry assignment
--------------------------------------

The selector's static geometry table
(:data:`~orpheus.numerics.quadrature.registry.GEOMETRY_ANGULAR_SYMMETRY`)
records **both halves** of each geometry's angular symmetry, because
the two halves place two different demands on a quadrature.

A geometry's symmetry group :math:`G` splits by how the action is
used. The continuous part :math:`G^0` is **spent** by the dimensional
reduction — a slab integrates the azimuth out analytically — and what
it determines is the angular *domain*, the quotient :math:`S^2/G^0`.
The finite residual :math:`\Gamma = G/G^0` is still **owed**: it
cannot be integrated away, so a quadrature must realize it as a
permutation of the ordinates. This is the half a reflecting boundary
consumes — the face reflection maps ordinate :math:`m` to ordinate
:math:`m'` exactly only when the node set is closed under it.

.. list-table::
   :header-rows: 1
   :widths: 18 12 12 12 46

   * - Geometry
     - :math:`G^0` (spent)
     - Domain :math:`S^2/G^0`
     - :math:`\Gamma` (owed)
     - Rationale
   * - ``"slab"``
     - :math:`SO(2)`
     - :math:`[-1,1]`
     - :math:`Z_2`
     - 1-D in :math:`z`. Azimuthal rotation about the normal is
       integrated out, leaving :math:`\mu = \cos\theta` alone. Owed is
       :math:`\mu \to -\mu`, which pairs the two sweep senses;
       Gauss-Legendre nodes are symmetric (Stoer-Bulirsch §3.6), so it
       holds.
   * - ``"sphere"``
     - :math:`SO(2)`
     - :math:`[-1,1]`
     - :math:`Z_2`
     - 1-D **radial** spherical SN reduces to GL on :math:`\mu_r`, the
       cosine of the angle between ordinate and radial direction
       (Lewis & Miller 1993 §4.4). The continuous problem is
       :math:`O(3)`-symmetric; the radial reduction spends the azimuth
       about :math:`\hat r`. Here the spent half is not free — its
       fiber action reappears in the sweep as the
       angular-redistribution :math:`\alpha` term.
   * - ``"cylinder"``
     - trivial
     - :math:`S^2`
     - :math:`D_{2h}`
     - An axisymmetric cylinder is :math:`\phi`-independent in
       *space*, which does not reduce the *angular* domain: both
       angular degrees of freedom survive. The cylindrical SN sweep
       also requires per-:math:`\mu` polar-level structure; request it
       via ``level_structured=True``.
   * - ``"cartesian2d"``
     - trivial
     - :math:`S^2`
     - :math:`D_{2h}`
     - :math:`D_{2h} \cong (\mathbb{Z}_2)^3` is generated by the three
       coordinate-plane mirrors, and its chambers are exactly the
       octants the sweep decomposes into.

.. admonition:: What this table looked like before 2026-08-02, and why
                it could not work
   :class: warning

   The table held a single group per geometry, and it recorded the
   **spent** half: ``slab``, ``sphere`` and ``cylinder`` all read
   :math:`SO(2)`. That is a true statement about the *geometry* and a
   useless one for *selecting a quadrature*, because **no finite point
   set on** :math:`S^2` **is** :math:`SO(2)`-**closed**. The gate
   :math:`G_{\text{geom}} \subseteq G_Q` was therefore unsatisfiable
   by any discrete azimuthal rule, and could only ever pass on a false
   declaration — which is exactly what it did: ``product_mu_phi``
   advertised :math:`SO(2)`, and that single falsehood was the only
   reason the gate admitted the product rule for a cylinder.

   ``cartesian2d`` read :math:`O_h`, an over-claim by a factor of 6:
   :math:`O_h` demands the :math:`x \leftrightarrow z` exchange and
   the diagonal mirrors, which are symmetries of a *cube* and never of
   a z-uniform problem. It rejected the even-:math:`n_\phi` product
   rule, which genuinely does carry every symmetry that geometry
   needs.

   The two errors were compensating, which is why the selector looked
   healthy: the rule-side field was doing double duty as a *domain*
   tag and a *symmetry* tag, so a wrong domain answer and a wrong
   symmetry answer cancelled. Splitting spent from owed, and computing
   :math:`\operatorname{Sym}(Q)` from the nodes instead of reading a
   declaration, makes the gate both satisfiable and discriminating.

Worked examples
---------------

The four canonical happy-path selections:

* ``select_quadrature("slab", target_degree=15)``: the three
  :math:`S^2` rules fail stage 0 ("geometry 'slab' discretises
  :math:`[-1,1]` … but the rule's nodes live on :math:`S^2`"), leaving
  GL1D (n=8 → degree 15, 8 nodes).

* ``select_quadrature("sphere", target_degree=15)``: same path as
  slab — identical spent/owed split after the 1-D radial reduction.

* ``select_quadrature("cylinder", target_degree=5,
  level_structured=True)``: ``ProductQuadrature`` with
  ``n_mu=3, n_phi=6``. GL1D is rejected on *domain* (a
  :math:`\mu`-marginal cannot carry two angular degrees of freedom),
  which is a stronger and earlier objection than "no level structure".

* ``select_quadrature("cylinder", target_degree=4,
  level_structured=True)``: the inversion asks for ``n_phi=5``, and an
  **odd** azimuthal count is not :math:`D_{2h}`-invariant —
  :math:`\sigma_x` sends :math:`\phi \to 180^\circ - \phi`, mapping the
  :math:`0^\circ` node onto :math:`180^\circ`, which is not a node of
  a 5-fold grid. The selector refuses it and falls back to
  ``LevelSymmetricSN``. This is ERR-042 expressed structurally rather
  than as a hand-written guard, and it is the decisive behavioural
  difference from the old table.

* ``select_quadrature("cartesian2d", target_degree=5)``: Lebedev,
  LS_N and the even-:math:`n_\phi` product rule all satisfy
  :math:`D_{2h}`; Lebedev wins on cost (14 nodes vs 48 and 18). Only
  GL1D is rejected, on domain. The log shows the cost losers as
  *valid candidates that lost*, not rejections — their names are
  absent from ``log.rejected``.

Explainability log
------------------

Every :func:`~orpheus.numerics.quadrature.registry.select_quadrature`
call returns a
:class:`~orpheus.numerics.quadrature.registry.SelectionLog` carrying
the full decision provenance:

.. code-block:: python

   measure, log = select_quadrature("cartesian2d", target_degree=5)
   print(log.summary())
   # select_quadrature(geometry='cartesian2d', target_degree=5)
   # -> LebedevSphere({'order': 5}) [1 rejected]

   for name, reason in log.rejected:
       print(f"  - {name}: {reason}")
   # - GaussLegendre1D: domain mismatch: geometry 'cartesian2d'
   #   discretises S^2 (= S^2/Trivial), but the rule's nodes live
   #   on [-1,1]

The rejection names its *stage*, and the stage is the diagnosis. A
``domain mismatch`` says the rule is the wrong shape of object for
this geometry; a ``symmetry mismatch`` says it is the right shape but
breaks a symmetry the geometry needs, and names the group:

.. code-block:: python

   _, log = select_quadrature(
       "cylinder", target_degree=4, level_structured=True
   )
   print(dict(log.rejected)["ProductQuadrature"])
   # symmetry mismatch: geometry 'cylinder' owes D_2h, which the
   # rule's nodes at {'n_mu': 3, 'n_phi': 5} are not invariant under

When no rule fits the constraints, :func:`select_quadrature` raises
:class:`~orpheus.numerics.quadrature.registry.QuadratureSelectionError`
with a fully-populated ``exc.log`` so the caller can introspect the
failure without re-running the search.

Why a registry and not a hardcoded ladder
------------------------------------------

``solve_sn`` keeps explicit quadrature-passing as the canonical API:
explicit is better than implicit, and a quadrature is a load-bearing
modelling choice the user must own. The registry is **opt-in
convenience** — a default for prototyping and a documentation
artifact. The structural-flag tags double as Sphinx teaching content;
each :class:`~orpheus.numerics.quadrature.registry.QuadratureSpec`
docstring narrates the trade-off in the rule's design space.

Domain-specific quadrature: the Quadrature class
================================================

The SN solver consumes angular quadratures through ~50
attribute-access sites (``quad.mu_x``, ``quad.weights``,
``quad.reflection_index('x')``, ``quad.spherical_harmonics(L)``, …)
in sweeps, operators, mesh constructors, and solvers.
Re-routing every site through
:meth:`DiscreteMeasure.integrate` would impose a per-site method-call
overhead on hot inner loops where the attribute access
costs nothing, and would needlessly couple the entire SN module to
the measure API.

The design (see ``.claude/plans/sn_reshape.md``) preserves both
worlds by layering the SN-side derived data over the mathematical
rule:

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
  - :func:`~orpheus.numerics.quadrature.periodic_trapezoid` — the
    :math:`n`-point rule on the circle :math:`S^1`, exact for
    *trigonometric* polynomials of degree :math:`n - 1`. The
    azimuthal factor of the product rule below.
  - :func:`~orpheus.numerics.quadrature.product_mu_phi` — the
    polar-times-azimuthal product rule used by the cylindrical SN
    sweep.

  The circle rule is worth reading closely even though it is
  one-dimensional, because it is where two of this page's
  distinctions become concrete rather than notional.

  **A rule on a domain is one object, and the domain decides the
  claim.** Its nodes coincide numerically with
  :func:`~orpheus.numerics.measure.equispaced`'s on
  :math:`[0, 2\pi]` — and the two claims are incomparable: degree
  :math:`n - 1` against the Fourier basis on :math:`S^1` versus
  degree :math:`1` against polynomials on an interval. Composing a
  polar rule with an azimuthal one by taking :math:`\min` of those
  two bare integers gives **1**, where the sphere product rule's
  true degree is :math:`\min(2n_\mu - 1,\, n_\varphi - 1)`. That is
  the whole argument for
  :class:`~orpheus.numerics.exactness.ExactnessClaim`.

  **A parameter no gate can see still decides something.** The
  rule's ``shift`` — the rotation of the node lattice, in units of
  one step — provably cannot change the degree (it enters the
  aliasing sum only as a phase on an already-vanishing factor). What
  it decides instead is whether a node sits *on* a mirror axis, and
  hence whether the singular set :math:`\Sigma = \{\xi = 0\}` of a
  fold through that mirror is empty. Exactly two shifts in
  :math:`[0, 1)` leave the node set mirror-symmetric at all —
  :math:`0` and :math:`\tfrac{1}{2}` — so the pair is a
  classification rather than a choice of convention, and the shift
  is a parameter of the rule rather than a boolean flag on it.
  Because the nodes are generated as roots of unity by
  :func:`~orpheus.numerics.roots_of_unity.roots_of_unity`,
  :math:`\Sigma` is decided by ``sin == 0.0`` — an equality, not a
  tolerance whose ``eps`` would otherwise set a fold's
  well-posedness condition.

* The single :class:`~orpheus.numerics.quadrature.Quadrature` class
  wraps a :class:`DiscreteMeasure` and exposes named classmethod
  factories — :meth:`Quadrature.gauss_legendre(n)
  <orpheus.numerics.quadrature.Quadrature.gauss_legendre>`,
  :meth:`Quadrature.lebedev(order)
  <orpheus.numerics.quadrature.Quadrature.lebedev>`,
  :meth:`Quadrature.level_symmetric(sn_order)
  <orpheus.numerics.quadrature.Quadrature.level_symmetric>`,
  :meth:`Quadrature.product(n_mu, n_phi)
  <orpheus.numerics.quadrature.Quadrature.product>` — each returning
  a ``Quadrature`` around the canonical measure produced by the rule
  function above, plus the SN-side derived data (reflection
  partners, octant partition, level structure) cached at
  construction. The legacy ``mu_x`` / ``mu_y`` / ``mu_z`` /
  ``weights`` / ``N`` surface survives as ``@property`` views over
  the underlying ``measure.nodes`` / ``measure.weights``, so the ~50
  consumer sites see no API change.

.. note::

   **Superseded architecture (R-1 Phase A detour-C).** The original
   Issue-4 design shipped a *bridge pattern*: four SN-side adapter
   classes at ``orpheus.sn.quadrature`` (``GaussLegendre1D`` /
   ``LebedevSphere`` / ``LevelSymmetricSN`` / ``ProductQuadrature``)
   plus an ``AngularQuadrature`` Protocol, each caching denormalised
   ``mu_x`` / ``mu_y`` / ``mu_z`` / ``weights`` *fields* over the
   measure. R-1 Phase A detour-C retired all four classes and the
   Protocol, collapsing them into the named classmethod factories on
   the single :class:`Quadrature` class. The Pattern-7 hazard of the
   old design (four parallel adapter classes each with their own
   denormalised storage) is closed by construction: there is now
   exactly one producer of the ordinate data — the ``measure`` — and
   the named accessors are derived views with no separate storage.

The named accessors keep the same :math:`O(1)` cost as the old
cached fields. The underlying :class:`DiscreteMeasure` is exposed as
``quad.measure`` for callers that want the structurally-richer
object — composability via
:meth:`~DiscreteMeasure.__mul__` /
:meth:`~DiscreteMeasure.pushforward` /
:meth:`~DiscreteMeasure.restrict`, plus the ``invariance_group``
that the
:func:`~orpheus.numerics.quadrature.registry.select_quadrature`
registry (Issue 5; see "Quadrature selection algorithm" above)
consumes for automated rule selection.

The bit-identical contract
--------------------------

The consolidation preserved the exact node/weight arrays the four
legacy adapter classes produced, pinned by two test layers:

* **Bit-identical foundation tests** at
  ``tests/numerics/test_rules_1d.py`` /
  ``tests/numerics/test_rules_product.py`` /
  ``tests/numerics/test_rules_sphere.py`` use
  :func:`numpy.array_equal` (not :func:`numpy.allclose`) to verify
  that each rule function's nodes/weights are exact-bit equal to the
  legacy reference output.
* **Regression snapshots** rerun every SN configuration end-to-end
  through the sweep / operator / power-iteration pipeline and assert
  match against the frozen baselines at the iterative solver's
  convergence floor (``rtol=1e-12``, ``atol=1e-13``).

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

- **Wave 2 (MoC migration)** will consume
  :class:`~orpheus.numerics.measure.BundleMeasure` for ray-bundle
  quadratures whose fibers (the parallel rays through the geometry)
  vary with the base direction. The bundle abstraction is shipped
  in Phase 0 because retrofitting it later would force a rewrite
  of MoC's quadrature handling.

Issue 5's ``select_quadrature`` registry has landed and is documented
in the "Quadrature selection algorithm" section above
(:eq:`quadrature-selection-criterion`).

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
