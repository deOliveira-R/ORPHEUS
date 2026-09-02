.. _discrete-measures:

============================================
Discrete Measures and Quadrature Composition
============================================

.. note::

   **Stub page** — full narrative pending under Issue 18 of
   ``.claude/plans/sn_reshape.md``. This stub installs the
   load-bearing equation labels (:eq:`discrete-measure-integrate`,
   :eq:`discrete-measure-pushforward`,
   :eq:`discrete-measure-quotient`, :eq:`folded-level-arc`,
   :eq:`bundle-measure-disintegration`) and the cross-references
   that downstream theory pages and tests rely on.

Key Facts
=========

- A **quadrature rule** is a discrete measure
  :math:`\mu = \sum_i w_i \, \delta_{x_i}` on a measurable space
  :math:`(\mathcal{X}, \Sigma)`. Lebesgue integration of a test
  function against :math:`\mu` is the familiar quadrature sum
  (:eq:`discrete-measure-integrate`).
- The canonical composition operations exposed by the project's
  primitive :class:`~orpheus.numerics.measure.DiscreteMeasure`:

  - **Tensor product** ``μ * ν`` — product measure
    :math:`\mu \otimes \nu` on
    :math:`\mathcal{X} \times \mathcal{Y}`.
  - **Direct sum** ``μ + ν`` — concatenation of atoms on a shared
    space; the foundation of composite (panelled) :term:`quadrature`.
  - **Pushforward** ``μ.pushforward(φ)`` — image measure
    :math:`\varphi_* \mu` (:eq:`discrete-measure-pushforward`).
  - **Quotient** ``μ.quotient(G)`` — the fold
    :math:`\mu/G` onto the orbifold :math:`\mathcal{X}/G`
    (:eq:`discrete-measure-quotient`): pushforward onto orbit
    representatives, then :meth:`consolidate
    <orpheus.numerics.measure.DiscreteMeasure.consolidate>` of the
    collapsed atoms. Refuses unless the measure is *certified*
    :math:`G`-invariant. Mass is preserved — the discriminator
    against restriction, which drops it. Sidecar level structure
    *descends* along a fiberwise quotient by pure selection, and on
    a mirror fold each polar level becomes a single **arc** whose
    stored :math:`\eta`-order IS the fiber order in march
    orientation (:eq:`folded-level-arc`).
  - **Restriction** ``μ.restrict(E)`` — indicator-multiplication
    :math:`\mathbf{1}_E \cdot \mu`; supports half-range SN :term:`sweeps <sweep>`
    and bundle cuts.
  - **Orbit-space declaration** ``μ.on_orbit_space(M/H)`` — the *same
    atoms*, re-read as chart coordinates of an orbit space. It moves no
    node, changes no weight and preserves mass exactly, because it
    applies nothing: only what the measure KNOWS about its support
    changes (:ref:`measure-on-orbit-space`). It is how the slab's
    :math:`\mu`-rule says it lives on :math:`S^2/SO(2)_x` rather than on
    the interval a chart happens to map onto. Neither a
    ``pushforward`` (no map) nor a ``quotient`` (nothing is folded).

- :class:`~orpheus.numerics.measure.BundleMeasure` carries the
  **disintegration** :eq:`bundle-measure-disintegration`: a base
  measure on :math:`\mathcal{B}` plus, for each base point, a fiber
  measure on :math:`\pi^{-1}(b)`. SN does not consume bundles in the
  Wave A campaign; MoC ray-bundle quadratures (Wave 2) will.
- The support (``μ.support``) **is the point set the atoms live on**, a
  :class:`~orpheus.numerics.manifold.Manifold` — with ``dim``,
  ``contains``, ``__mul__`` and a catalogued ``quotient``. (Distinct from
  the derived ``μ.space``, the induced discrete-:math:`L^2`
  :class:`~orpheus.numerics.space.FunctionSpace`.)

  ⛔ ~~*It was a runtime* ``str`` *— "Python generics over a
  measurable-space type are not expressive enough without runtime
  overhead, so the project uses opaque string tags for sanity checks and
  documentation only."*~~ **RETIRED 2026-09-01 (tracker 2.0c).** ⭐ The
  ruling that sentence cites was *answered, not overruled*: it rejects a
  **phantom type parameter** — erased at runtime, buying nothing — and it
  is still correct. What ships instead is a closed sum of ordinary
  **values**, which costs no generics and no runtime overhead, and which
  the tag could not do. The three support rules in the propagation table
  below were string interpolation and are now the manifold's own algebra:
  :math:`\mu \otimes \nu` multiplies the supports, the fold takes
  ``support.quotient(G)``, and the direct sum compares objects.
- A measure **generates space factors, not only whole spaces**:
  :meth:`μ.axis(label) <orpheus.numerics.measure.DiscreteMeasure.axis>`
  is the axis-composed sibling of ``μ.space``, minting one
  :class:`~orpheus.numerics.axis.Axis` that keeps the weights, drops
  the nodes, and records ``μ`` itself as its generator so the
  forgetting is recoverable.
  :meth:`Quadrature.axis <orpheus.numerics.quadrature.directional.Quadrature.axis>`
  is the angular sibling and records the RULE, because
  ``level_indices`` lives on the
  :class:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure`
  side-channel and not among the measure's five fields. The doctrine, the identity
  exclusion and the seams live on
  :ref:`the space layer's page <spaces-axis-generator>`; ``kind`` is
  never a parameter of either mint — the generator's TYPE implies it.


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

Two derived operations build on it. **Consolidation**
(:meth:`~orpheus.numerics.measure.DiscreteMeasure.consolidate`)
merges coincident atoms by summing their weights — the reduction a
non-invertible :math:`\varphi` leaves undone. It moves no node and
changes no integral, so it is the one operation that *preserves*
both metadata claims. **Quotient**
(:meth:`~orpheus.numerics.measure.DiscreteMeasure.quotient`) is the
composite of the two, with the defining law that makes it a quotient
and not a convention:

.. math::
   :label: discrete-measure-quotient

   \mu / G \;=\; \mathrm{pushforward}(\varphi_{\mathrm{rep}})
   \texttt{.consolidate()},
   \qquad
   \int f \, d(\mu/G) \;=\; \int f \, d\mu
   \quad \text{for every } G\text{-invariant } f,

.. (vv-status rationale) discrete-measure-quotient: Verified by the
   foundation gates in ``tests/numerics/test_measure.py`` — two
   independent representative sections (the verb's first-appearing
   member vs the geometric ξ → |ξ| reference) realize the same
   quotient orbit-by-orbit; the G-invariant-integral law with a
   ξ-odd knob-dependent control; orbit-stabilizer weights read off
   the certificate plus the Burnside count; T24 free-action uniform
   2w bit-exact; refusal without a certificate; and the
   consumes-the-symmetry idempotence pair. Software invariants of
   the measure composite, not a physics-equation claim with an
   L0..L3 ladder slot.
.. vv-status: discrete-measure-quotient documented

where :math:`\varphi_{\mathrm{rep}}` sends each node to its orbit's
first-appearing member (the only *group-generic* section — a
geometric section such as :math:`\xi \to |\xi|` exists only for a
mirror). The summed weights are the orbit-stabilizer weights
:math:`W = w \cdot |G|/|\mathrm{Stab}(x)|`, derived rather than
chosen. The quotient is defined only on a measure whose
:math:`G`-invariance is *certified*
(:func:`~orpheus.numerics.symmetry.orbit_certificate`); when the
action is **free** — empty singular set :math:`\Sigma = \varnothing`
— every orbit has length :math:`|G|` and the weights collapse to a
uniform :math:`|G| \cdot w`, which is the fold's well-posedness
condition: a rule is admissible for a fold iff it places no node on
:math:`\Sigma`.

A quadrature rule often travels with *sidecar structure* — for the
SN product rules, the per-level metadata
(:class:`~orpheus.numerics.quadrature.LevelStructure`) the
cylindrical sweep marches by. That structure **descends along the
quotient by pure selection** (:meth:`LevelStructure.quotient
<orpheus.numerics.quadrature.rules_sphere.LevelStructure.quotient>`):
a quotient never moves a node — every folded atom is an orbit
representative, a bit-copy of a parent node — so the angular charts
descend as selections and each level's order descends as a
*restriction* (a subsequence of a sorted sequence is sorted; the
sort convention is spelled once, in the producer). The descent is
defined exactly when the action is **fiberwise**: it may permute a
level's circle but must not move weight between levels, certified
per level by mass conservation. A level-merging fold such as the
:math:`\mu`-mirror :math:`\sigma_z` is refused.

On the :math:`\sigma_y` mirror fold the descent has a sharp
geometric consequence — **a level becomes an arc**. Each level's
circle meets the singular set :math:`\Sigma = \{\xi = 0\}` at
:math:`\omega \in \{0, \pi\}`, and the fold keeps a single arc
between those points, on which

.. math::
   :label: folded-level-arc

   \frac{\partial \eta}{\partial \omega}
   \;=\; -\sin\theta \, \sin\omega \;<\; 0
   \quad \text{on } \omega \in (0, \pi)
   \qquad\Longrightarrow\qquad
   \operatorname{sort}_{\uparrow \eta}
   \;=\; \operatorname{traverse}_{\downarrow \omega} .

.. (vv-status rationale) folded-level-arc: Verified by the Q5.3
   foundation gates in ``tests/numerics/test_rules_sphere.py`` —
   ``test_a_folded_level_is_an_arc_in_march_order`` asserts, per
   level of four folded configs (staggered 4×8 / 2×4, staggered
   3×5 with one Σ endpoint on the arc, node-aligned 4×8 with both),
   strict η-injectivity AND strict ω-descent along the stored
   order; the selection/restriction mechanics, the fiberwise mass
   certificate, and both refusal arms are gated alongside. A
   geometry-of-the-rule invariant, not a physics-equation claim
   with an L0..L3 ladder slot.
.. vv-status: folded-level-arc documented

The :math:`\eta`-key — 2-to-1 on a full level, the mechanism behind
issue #326 — is *injective* on the arc, so the stored
:math:`\eta`-ascending order becomes a genuine ordering of the
fiber: the arc traversed in strictly decreasing :math:`\omega`,
which is the azimuthal march order. One order, seen through two
charts. The two ordering accessors that the adjudication of #326
had deliberately kept side by side merged on this theorem — the
:math:`\omega`-ascending ``fiber()`` accessor retired as the same
total order in the opposite orientation, with no consumer.

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

The composition operations propagate the
:class:`~orpheus.numerics.measure.DiscreteMeasure` metadata as
follows. Where a field is marked **dropped**, the operation
generally cannot guarantee the field is still meaningful, and the
result carries ``None`` until the caller re-asserts it via
:meth:`~orpheus.numerics.measure.DiscreteMeasure.with_metadata`.

.. list-table::
   :header-rows: 1
   :widths: 18 18 16 16 16 16

   * - Operation
     - ``support``
     - ``n_points``
     - ``dim``
     - ``degree_of_exactness``
     - ``invariance_group``
   * - tensor product ``μ * ν``
     - ``f"{μ.support} × {ν.support}"``
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
     - explicit ``new_space=`` argument or ``f"φ_*({μ.support})"``
     - unchanged
     - inferred from ``φ``'s output shape
     - dropped (φ-image of a polynomial-exact rule is not
       polynomial-exact in general — undefined unless ``φ`` is
       affine)
     - dropped (φ generally does not preserve a group action)
   * - consolidate ``μ.consolidate()``
     - unchanged
     - one atom per distinct position (:math:`\leq N`)
     - unchanged
     - **preserved** (merging by position changes no integral)
     - **preserved** (a group permuting the atoms still permutes
       the merged ones)
   * - quotient ``μ.quotient(G)``
     - ``f"{μ.support}/{G.name}"``
     - one atom per :math:`G`-orbit
     - unchanged
     - dropped (an honest claim would be against the pushforward
       reference :math:`\varphi_* \lambda`, a type with no consumer
       yet)
     - dropped (the fold CONSUMES the symmetry: the representative
       section is not :math:`G`-equivariant, so no residual
       invariance is guaranteed)
   * - restriction ``μ.restrict(E)``
     - unchanged
     - ``mask.sum()``
     - unchanged
     - dropped (restricted rule is not Gauss on the cut domain)
     - dropped (caller may re-tag if :math:`E` is invariant)
   * - orbit-space declaration ``μ.on_orbit_space(M/H)``
     - the orbit space itself, e.g. ``S^2/SO2_x`` (**raises** unless its
       ``realization`` is the current support)
     - unchanged
     - unchanged
     - **preserved** — the reference measure is a measure on the
       *chart*, and the chart is unchanged
     - dropped (a subgroup of :math:`O(3)` is a claim about an
       **embedding**, and the orbit space fixes one the chart did not;
       the adopter re-tags for the named axis)

Be honest about what is dropped: the field becoming ``None`` is the
correct behaviour, not a degradation. If a downstream consumer
needs the field, it must re-establish the claim and stamp it back
with :meth:`~orpheus.numerics.measure.DiscreteMeasure.with_metadata`.

.. _measure-on-orbit-space:

``on_orbit_space`` — declaring the point set, applying nothing
---------------------------------------------------------------

The last row of that table is the newest verb (tracker 2.4 of #429,
2026-09-01) and the odd one out: every other operation in this section
*computes* something. This one computes nothing at all.

A rule built on a manifold :math:`C` is a rule on **every** orbit space
:math:`M/H` whose invariant chart lands on :math:`C`. Declaring which
one is a statement about the measure's *type*, and there is no
arithmetic that can supply it:

.. math::

   \mu \;=\; \sum_i w_i\,\delta_{x_i} \ \text{on } C
   \qquad\longmapsto\qquad
   \mu \;=\; \sum_i w_i\,\delta_{x_i} \ \text{on } M/H,
   \qquad \text{realization}(M/H) = C .

Same :math:`x_i`, same :math:`w_i` — `[M]` literally the same arrays,
not copies (``adopted.nodes is original.nodes``). Same integral, same
mass, same exactness claim.

**Why it is neither of the two verbs it resembles.**

* Not a :meth:`~orpheus.numerics.measure.DiscreteMeasure.pushforward`:
  a pushforward applies a map :math:`\varphi` and therefore *moves the
  points*. Nothing is applied here.
* Not a :meth:`~orpheus.numerics.measure.DiscreteMeasure.quotient`: a
  fold starts from a measure on the **base** :math:`M`, identifies
  orbits, and drops the atom count. A :math:`\mu`-rule was never on
  :math:`S^2` to begin with — there is nothing to fold.

⟹ the discriminator worth remembering is *where the measure starts*.
``quotient`` walks **down** from the base; ``on_orbit_space`` stays put
and re-labels the chart it was already on.

**The one precondition.** ``on_orbit_space`` refuses unless the orbit
space's ``realization`` is this measure's current support. `[M]` handing
a :math:`\mu`-rule the mirror quotient ``S^2/sigma_y``, whose chart is
the disk :math:`D^2`:

.. code-block:: text

   a measure on '[-1,1]' cannot be read on 'S^2/sigma_y': that orbit
   space's chart is 'D^2'. A rule declares the orbit space whose CHART
   it was built on; to fold a rule on the base manifold, use quotient()

The message names the mismatch **and** the other verb, because "wrong
orbit space" and "you meant to fold" are the two ways a caller arrives
here.

**What it buys, measured.** The slab's polar rule is the first consumer.
`[M]` 2026-09-01, ``Quadrature.gauss_legendre(8)``: ``support.name``
moves ``'[-1,1]' → 'S^2/SO2_x'`` and ``space.name`` moves
``'L2[[-1,1]]' → 'L2[S^2/SO2_x]'``, with nodes and weights bit-identical
to ``gauss_legendre_on_mu(8)``. Before the declaration an 8-node slab
**angular** space and an 8-node **spatial** rule on :math:`[-1,1]` built
from those same arrays compared ``==`` *and* hash-equal; after it they
do not. The full argument, and why the axis in ``SO2_x`` is a parameter
rather than a convention, is on
:doc:`/theory/foundations/manifolds` —
:ref:`manifold-on-orbit-space` and
:ref:`manifold-so2-axis-is-a-parameter`. **Edited there, consumed
here.**


.. _measure-has-versus-spent:

``invariance_group`` and ``quotient_group`` — what a measure HAS and what it SPENT
-----------------------------------------------------------------------------------

One of the fields in the propagation table above names a subgroup of
:math:`O(3)`, and it has a sibling that does **not** appear there —
because the sibling is not a stored field, so there is nothing for an
operation to propagate. The two are almost complementary.
:attr:`invariance_group
<orpheus.numerics.measure.DiscreteMeasure.invariance_group>` is the
stored one: a recorded subgroup under which the atoms — nodes and
weights — are closed, i.e. what the measure **HAS**. It is a
*declaration*, not a computed stabiliser, and its ``None`` means
*unspecified*, never *trivial*.
:attr:`quotient_group
<orpheus.numerics.measure.DiscreteMeasure.quotient_group>` is *derived*
and never stored (tracker 2.0c of #429): the group the measure's support
was folded by, read straight off :attr:`Quotient.by
<orpheus.numerics.manifold.Quotient.by>`, so there is no second home for
it to drift from — the same argument
:meth:`~orpheus.numerics.measure.DiscreteMeasure.restrict`,
:meth:`~orpheus.numerics.measure.DiscreteMeasure.consolidate` and
:meth:`~orpheus.numerics.measure.DiscreteMeasure.partition_by` would
otherwise each have to honour by hand, every time they rebuild a
measure.

For a **point set** the two come apart in both directions, and the
shipped angular rules demonstrate it. The denominator is every
``classmethod`` factory on
:class:`~orpheus.numerics.quadrature.directional.Quadrature`, enumerated
by ``vars(Quadrature)`` — `[M]` **five of five**, so this is the family
and not a sample (``vv-principles`` #31's finite-roster corollary):

.. list-table:: `[M]` 2026-09-01 — all five shipped ``Quadrature`` factories
   :header-rows: 1
   :widths: 26 20 27 27

   * - Rule
     - ``support.name``
     - HAS (``invariance_group``)
     - SPENT (``quotient_group``)
   * - ``Quadrature.lebedev(17)``
     - ``'S^2'``
     - ``OctahedralOh``
     - ``None``
   * - ``Quadrature.level_symmetric(8)``
     - ``'S^2'``
     - ``OctahedralOh``
     - ``None``
   * - ``Quadrature.product(4, 8)``
     - ``'S^2'``
     - ``Dnh(8)``
     - ``None``
   * - ``Quadrature.gauss_legendre(8)``
     - ``'S^2/SO2_x'``
     - ``Mirror('x')``
     - ``SO2('x')``
   * - ``Quadrature.folded_product(4, 8)``
     - ``'S^2/sigma_y'``
     - ``None``
     - ``Mirror('y')``

The slab's polar rule carries **two different groups in the two slots**;
and the :math:`\sigma_y` fold HAS nothing *because* it spent
:math:`\sigma_y` — a fold keeps one representative per orbit, so the
surviving set is no longer closed under the mirror it was folded by.
**Spending a symmetry destroys having it**, which is why reading either
slot as the other is ``plan-authoring`` §3's ambiguous-name hazard.

⭐ **A BASIS carries only one of these, and that asymmetry is a theorem,
not an omission.** A function on an orbit space :math:`M/H`, pulled back
to :math:`M`, is exactly an :math:`H`-invariant function — being a
function on the quotient *is* being :math:`H`-invariant — so for
functions "has" and "spent" are one property.
:attr:`Basis.invariance_group
<orpheus.numerics.basis.base.Basis.invariance_group>` is consequently a
single ``@final`` slot, named for what the basis HAS and **derived** by a
``match`` on what its :attr:`~orpheus.numerics.basis.base.Basis.domain`
SPENT (#429 tracker 2.1b). The two ``None``\ s then mean different
things for the same reason: a full-sphere rule's ``quotient_group`` is
``None`` because it spent nothing, while the full-sphere harmonics'
``invariance_group`` is ``Trivial`` because they *have* the trivial
group. The derivation, the lower-bound caveat, and the pairing verdict
it makes computable are on
:doc:`/theory/foundations/manifolds` —
:ref:`manifold-has-versus-spent` and
:ref:`manifold-invariance-pairing`. **Edited there, consumed here.**


.. _discrete-measure-partition:

Partition by labelling predicate
================================

A further composition operation (added in Wave 0 of the SN
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


.. _1-d-primitive-constructors:

1-D primitive constructors
==========================

The module ships three 1-D primitives:

- :func:`~orpheus.numerics.measure.gauss_legendre`
  ``(n)`` — Gauss-Legendre on :math:`[-1, 1]`, exact to algebraic
  degree :math:`2n - 1` **against** ``legendre``; weights sum to 2.
- :func:`~orpheus.numerics.measure.gauss_chebyshev`
  ``(n)`` — Gauss-Chebyshev (first kind) on :math:`[-1, 1]`, exact to
  algebraic degree :math:`2n - 1` **against** ``chebyshev_t``, whose
  weight is :math:`(1 - x^2)^{-1/2}`; weights sum to :math:`\pi`.
- :func:`~orpheus.numerics.measure.equispaced`
  ``(a, b, n)`` — midpoint rule on :math:`[a, b]`, exact to algebraic
  degree 1 against ``uniform([a,b])``; weights sum to :math:`b - a`.

.. important::

   **The first two agree on every attribute a reader would think to
   compare, and integrate differently.** ``[M]`` at :math:`n = 4` both
   carry support :math:`[-1,1]`, orthogonal system ``ALGEBRAIC`` and
   degree 7 — and on the *unweighted* :math:`\int_{-1}^{1} x^6 dx`,
   Legendre gives :math:`0.285714` (exactly :math:`2/7`) while
   Chebyshev gives :math:`0.981748`. The **reference measure** is the
   only thing that distinguishes them, which is why a degree is stored
   as half of an :class:`~orpheus.numerics.exactness.ExactnessClaim`
   rather than as a bare integer, and why the quadrature selector's V
   conjunct compares claims rather than degrees
   (:ref:`quadrature-selection-algorithm`, stage 2). A degree is an
   index into the orthogonal system *of a measure*; the integer alone
   does not say which.

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

The groups themselves carry an order relation — containment in the
:math:`O(3)` lattice. For two subgroups :math:`H, K \subseteq O(3)`,

.. math::
   :label: subgroup-of-o3-containment

   H \subseteq K
   \;\Longleftrightarrow\;
   \forall g \in H,\; g \in K.

.. (vv-status rationale) subgroup-of-o3-containment: Verified
   transitively by the foundation tests in
   :file:`tests/numerics/test_symmetry.py` — every named relation in
   the static lattice (``Trivial ⊂ σ_z ⊂ O_h ⊂ O(3)``,
   ``SO(2)_z ⊂ D_∞h ⊂ O(3)`` with the two OTHER axes asserted NOT
   inside D_∞h and all three inside SO(3), ``C_n ⊂ SO(2)_z`` for every
   n, the ``O_h \not\subset SO(3)`` improper-rotation test, the
   parameterised ``C_n ⊂ C_m \iff n | m`` rule, and reflexivity for
   every named entry) is asserted directly. The axis-dependent rows
   arrived with the SO(2) parameterisation on 2026-09-01; the row this
   comment used to name, ``SO(2) ⊂ O(2) ⊂ O(3)``, predates the O(2) →
   D_∞h rename recorded three paragraphs below.
   Tagged ``foundation`` rather than carrying a verification ladder
   slot because the containment lattice is a software invariant
   (Hamermesh 1962 §2.5; Stiefel & Fässler 1979 Ch. 4) — there is
   no L0..L3 physics-equation claim to ladder against.
.. vv-status: subgroup-of-o3-containment documented

The "if and only if" is the standard set-theoretic definition of a
subgroup; what it buys is a *decidable* order over the named entries
and the parameterised families, implemented by
:meth:`~orpheus.numerics.symmetry.SubgroupOfO3.contains` and its
reverse-direction synonym ``is_subgroup_of``.

.. warning::

   **This relation is not the quadrature-selection gate**, and reading
   it as one is precisely the error the selector shipped until
   2026-08-02. Two things break when it is used that way. Whether a
   *rule* carries a symmetry is a question about its **nodes** —
   answered by :eq:`discrete-measure-g-invariance` through
   :meth:`~orpheus.numerics.symmetry.SubgroupOfO3.is_invariant` —
   whereas a lattice query can only consult a *declared* tag, which may
   be false or merely under-claimed. And the geometry side is not one
   group at all: it splits into the continuous half the dimensional
   reduction spends (which fixes the angular **domain**) and the
   discrete half still owed (which the nodes must realize as a
   permutation). Both corrections are set out under
   :ref:`quadrature-selection-algorithm`.

   ⛔ This paragraph read that "quadrature selection in ORPHEUS
   therefore reduces to a containment check", with the geometry side
   given as its whole natural symmetry group (slab
   :math:`\to SO(2) \times \sigma_x`, sphere :math:`\to O(3)`,
   hexagonal lattice :math:`\to D_{6h}`), and closed by claiming that
   containment "is sufficient to preserve every symmetry the geometry
   exhibits". Those readings were retired on 2026-08-02 and survived
   here until 2026-08-14, three screens above the section that
   describes what replaced them — which is the worse failure mode,
   because a page that contradicts itself can be cited for either
   sentence.

ORPHEUS's relevant sub-lattice of :math:`O(3)` is **finite and
small**. Slab, sphere, 1-D / 2-D Cartesian geometries combined with
the four shipped quadrature families (Gauss-Legendre, Lebedev,
level-symmetric, product) span the `[M]` **six** named entries encoded
in :class:`~orpheus.numerics.symmetry.SubgroupOfO3` — ``Trivial``,
:math:`D_{\infty h}`, ``O_h``, ``I_h``, ``SO(3)``, ``O(3)`` — plus
`[M]` **four** parameterised families :math:`C_n`, :math:`D_{nh}`,
:math:`\sigma_a` and :math:`SO(2)_a`.

⛔ This sentence read *"the **seven** named entries … ``Trivial``,
``SO(2)``, …"* until 2026-09-01. The axial rotation group left the enum
that day and became :class:`~orpheus.numerics.symmetry.SO2`, a frozen
dataclass parameterised by its axis — the same move the reflection had
made on 2026-08-02, for the same reason and after the same kind of
wrong answer. **The enum is now exactly the groups that have no axis to
name.** The argument is
:ref:`manifold-so2-axis-is-a-parameter`; the one-line version is that
the tree carries two polar axes at once — the real spherical-harmonic
basis and the slab's marginal are about :math:`x`, every product rule's
polar factor is about :math:`z` — so a group tag that did not name its
axis was a claim about the wrong pole in one of its two uses.

.. note::

   The reflection group joined the parameterised families on
   2026-08-02. It was a named, parameter-free ``Z_2`` entry realized as
   :math:`\sigma_z`, under a docstring calling the plane a convention.
   `[M]` It is not: ``product(4, 3)`` is :math:`\sigma_z`-closed and
   **not** :math:`\sigma_x`-closed, and the tag answered ``True``.
   Worse, the tag meant two different things by node shape — the 3-D
   arm tested :math:`\sigma_z`, the 1-D arm tested :math:`x \to -x`
   plane-free — so on an asymmetric :math:`\mu`-set the two arms
   disagreed, the 3-D one *certifying* a set that violates
   :math:`\mu \to -\mu`. Naming the plane is what makes the two arms
   answer one question.

   ⭐ **The axial rotation group followed on 2026-09-01, and the
   repetition is the point.** ``SO2`` was a named, parameter-free entry
   realized about :math:`z`: its exactness criterion asked whether every
   node had :math:`\rho = \sqrt{x^2 + y^2} = 0`. `[M]` run against the
   slab marginal's own embedded nodes — whose first row is
   :math:`(-0.96028\ldots,\,0,\,0)` — that criterion returns **False**,
   so the one shipped rule whose orbit space *is* :math:`S^2/SO(2)`
   reported that it was not invariant under the group it had been
   quotiented by. Same defect shape, one axis over, one month later.
   ⟹ **a group realized on a coordinate axis carries that axis as
   data.** ``ORPHEUS`` has now reached that ruling twice by shipping a
   wrong answer first; :math:`C_n` is the remaining candidate, and it is
   deliberately un-parameterised until a second consumer appears
   (:ref:`manifold-so2-axis-is-a-parameter`).

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
the input under nearest-neighbour matching.

A **continuous** group is not sampled at all — that is ERR-072, and the
criterion is exact. A finite point set is closed under :math:`SO(2)_a`
**iff every node lies ON axis** :math:`a`, because the orbit of any
point at distance :math:`\rho > 0` from that axis is a whole circle and
no finite set can contain one. So the predicate reads
:math:`\rho_a \le \texttt{atol}`, where :math:`\rho_a` is the norm of
the two coordinates *other* than :math:`a`'s.

⛔ **This paragraph read** *"the 1-D case (*:math:`SO(2)`*-invariance of
a measure on* :math:`[-1,1]`*) is trivial — there is no azimuthal
coordinate to rotate"* **until 2026-09-01, and that is exactly the
reading the axis parameter refutes.** It is not trivial and it is not
free: a 1-D measure's nodes are *embedded* before the predicate runs,
and the axis they land on decides the answer. `[M]` on
``Quadrature.gauss_legendre(8).measure``, whose support declares
:math:`S^2/SO(2)_x`, ``SO2('x').is_invariant`` is **True** while
``SO2('y')`` and ``SO2('z')`` are **False** — one true row of three, not
a vacuous yes. The embedding is read off the support: a rule on
:math:`S^2/SO(2)_a` embeds along :math:`a` (column 0 is the convention
only for a *bare* interval, which names no axis). `[M]`
``gauss_legendre_on_polar_orbit(8, "z")`` — the same eight floats,
declared about :math:`z` — is :math:`SO(2)_z`-invariant and **not**
:math:`SO(2)_x`-invariant, the exact mirror image.

The reflection :math:`x \to -x` remains the non-trivial *discrete* 1-D
check. That reflection is :math:`\sigma_x`, not "a reflection": under
the canonical :math:`(\mu, 0, 0)` embedding :math:`\sigma_y` and
:math:`\sigma_z` fix every node **pointwise** and hold trivially, while
:math:`\sigma_x` carries the whole content. The 1-D and 3-D checks are
the same question asked of the same embedding.

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

.. _quadrature-selection-algorithm:

Quadrature selection algorithm
==============================

Issue 5 of the SN reshape campaign installs a tagged registry over four
of the quadrature-rule primitives — Gauss-Legendre on :math:`\mu`,
Lebedev, level-symmetric :math:`S_N`, and the polar-times-azimuthal
product rule — plus a selector that picks the right
:class:`~orpheus.numerics.measure.DiscreteMeasure` for a
geometry-and-target-degree pair. The registry entries are the four
:class:`~orpheus.numerics.quadrature.registry.QuadratureSpec` values in
:data:`~orpheus.numerics.quadrature.registry.quadrature_registry`; the
selector lives at :mod:`orpheus.numerics.quadrature.registry`.

Selection is fundamentally a **five-stage filter** in priority order.
The stages are numbered here exactly as the selector numbers them, so a
rejection message read off a
:class:`~orpheus.numerics.quadrature.registry.SelectionLog` — which
always names its stage — points at the paragraph that explains it.

0. **Domain compatibility.** The rule's nodes must live on the angular
   domain the geometry's *dimensional reduction* left behind:
   :math:`\mathcal{D}_Q = S^2 / G^0_{\text{geom}}`. A slab integrates
   the azimuth out analytically, so its angular variable is
   :math:`\mu = \cos\theta` on :math:`[-1,1]`, and a rule whose nodes
   are points of :math:`S^2` is the **wrong shape of object** — not an
   over-resolved one, not an expensive one, but a rule for a different
   integral. A cylinder retains both angular degrees of freedom and
   wants :math:`S^2`.

   This stage is **not** a refinement of the symmetry stage below. It
   is the other half of the same decomposition of the geometry's
   angular symmetry (see
   :class:`~orpheus.numerics.quadrature.registry.AngularSymmetry`, and
   the subsection *Spent and owed* below), and the two conjuncts are
   logically independent: without stage 0, the symmetry stage alone
   admits a Lebedev rule for a slab, because an :math:`O_h`-invariant
   rule certainly satisfies the :math:`\sigma_x` a slab owes.

1. **G compatibility (symmetry).** The rule's node set must be closed
   under the geometry's *discrete residual* symmetry —
   :math:`\Gamma_{\text{geom}} \subseteq \operatorname{Sym}(Q)` — the
   half of the geometry's symmetry group that survives the reduction
   and must therefore be realized as a permutation of the ordinates. A
   quadrature with **less** symmetry imprints spurious low-order
   asymmetry on a symmetric problem (Lebedev 1976, §1; Stiefel &
   Fässler 1979 Ch. 5) and cannot represent a reflecting boundary
   exactly: the face reflection maps ordinate :math:`m` to ordinate
   :math:`m'` *exactly* only if the node set is closed under it. A
   quadrature with **more** symmetry is fine — the extra symmetry is
   unused, not violated.

   :math:`\operatorname{Sym}(Q)` is **computed from the rule's nodes**
   by :meth:`~orpheus.numerics.symmetry.SubgroupOfO3.is_invariant`
   applied to a generating set (:eq:`discrete-measure-g-invariance`),
   never read from a declared tag. The containment lattice
   :eq:`subgroup-of-o3-containment` is the order relation *on the
   groups*; it is not the gate, and routing the gate through it fails
   in **both** directions:

   - A declaration can be **false**, and then the gate passes on a lie.
     ``product_mu_phi`` advertised :math:`SO(2)`, which no finite point
     set on :math:`S^2` can satisfy, and that single falsehood was the
     only reason this gate admitted the product rule for a cylinder.
   - A declaration can be **true but not maximal**, and then the gate
     rejects a rule it must accept. Every measure is
     ``Trivial``-invariant, so tagging Gauss-Legendre ``Trivial`` is not
     a lie — and ``[M]`` a lattice query against that tag answers
     ``False`` for the slab's owed :math:`\sigma_x`, while
     ``admits_symmetry``, which asks the nodes, answers ``True``. The
     rule the slab must accept would be the one rejected.

   Asking the nodes cannot go wrong either way: a computed group cannot
   lie about the object it was computed from. This is pinned by
   ``test_selector_asks_the_nodes_not_the_declared_tag``, which injects
   the understated-but-true declaration precisely so the gate is not
   defended by a counter-example that a later bug fix could dissolve.
   All four shipped rules now declare their maximal group honestly —
   ``[M]`` :math:`\sigma_x`, :math:`O_h`, :math:`O_h` and
   :math:`D_{n_\varphi h}` — so on today's registry the two routes
   agree. That agreement is a property of today's tags, not of the
   design, which is why the principle is stated rather than measured.

   ⛔ This stage required :math:`G_{\text{geom}} \subseteq G_Q`, the
   geometry's whole group inside the rule's **declared** group, and
   named :eq:`subgroup-of-o3-containment` as the mechanism that drove
   it, until 2026-08-14. Every clause of that was retired on
   2026-08-02. :math:`G_{\text{geom}}` recorded the **spent**
   continuous half and became :math:`\Gamma_{\text{geom}}`, the
   **owed** discrete half; :math:`G_Q`, a declared parameter-free field
   on the spec, became :math:`\operatorname{Sym}(Q)`, computed from the
   instantiated nodes; and the lattice stopped being the mechanism. The
   retired gate was unsatisfiable by any discrete azimuthal rule and
   could only ever pass on a false declaration — see the admonition
   under *Geometry → angular-symmetry assignment* below for the
   compensating pair of errors that made it look healthy.

2. **V compatibility (exactness).** The rule's exactness **claim** must
   dominate the query's:
   :math:`\mathcal{E}(Q) \succeq (\lambda_{\text{geom}},\, d)`. One
   claim dominates another iff it is against the **same reference
   measure** *and* its degree is at least as large — the domination
   display under :eq:`quadrature-selection-criterion` states it
   formally.

   **The reference is not decoration.** A degree is an *index into the
   orthogonal system of a measure*, so the same integer means different
   things against different measures, and a rule can match on space, on
   orthogonal system **and** on degree while being exact against the
   wrong thing. That is not hypothetical, and this page's own two
   :ref:`1-D primitives <1-d-primitive-constructors>` are the witness:
   ``[M]`` :func:`~orpheus.numerics.measure.gauss_legendre` and
   :func:`~orpheus.numerics.measure.gauss_chebyshev` at :math:`n = 4`
   agree on all three (support :math:`[-1,1]`, system ``ALGEBRAIC``,
   degree 7) and differ on the **unweighted**
   :math:`\int_{-1}^{1} x^6 \, dx` — a monomial comfortably inside that
   degree — by **0.696**: Legendre returns :math:`0.285714`, which is
   exactly :math:`2/7`, while Chebyshev returns :math:`0.981748`, about
   :math:`3.4\times` the true value. Neither rule is broken; each is
   delivering its full advertised accuracy against the measure it was
   built for.

   Transport integrates :math:`\phi = \int \psi \, d\Omega`
   **unweighted**, so :math:`\lambda_{\text{geom}}` is Lebesgue measure
   on whatever angular domain the reduction left behind — ``legendre``
   on :math:`[-1,1]`, ``uniform(S^2)`` on the sphere. It is *derived
   from the spent group*, by :attr:`AngularSymmetry.reference
   <orpheus.numerics.quadrature.registry.AngularSymmetry.reference>`,
   exactly as the domain is (see *Spent and owed* below).

   Most rules have a parameter-dependent degree — :math:`2n - 1` for
   Gauss-Legendre, ``order`` for Lebedev,
   :math:`\min(2 n_\mu - 1,\, n_\phi - 1)` for the product rule — and
   the selector **inverts** it to pick the smallest parameter set
   reaching :math:`d`. ⚠ Level-symmetric :math:`S_N` has **no formula**:
   its realized degree is build-measured, and :math:`N - 1` is only the
   lower bound the inversion uses, which is why the inversion
   over-shoots (see the ⛔ blocks below for the measured table).
   Lebedev's gap structure (no rules at orders 33, 37, 39, 43, 45, 49)
   is handled by rounding up to the next tabulated order; if the target
   exceeds the table's top end (47 in scipy's tabulation), the rule is
   rejected at this stage with a clear message.

   **An inversion is a request, not evidence.** So the claim the *built*
   rule actually carries is verified, and there are three ways to fail
   it, each with its own gate in
   :file:`tests/numerics/test_registry.py`: the rule carries **no
   claim** at all (``test_a_rule_with_no_exactness_claim_at_all_is_refused``);
   its claim is against the **wrong reference**
   (``test_a_rule_exact_against_the_wrong_measure_is_refused``); or its
   degree falls **short of the target**, which is a spec whose own
   inversion over-promised and is therefore caught rather than trusted
   (``test_an_inversion_that_over_promises_is_caught_not_trusted``).
   Each arm is mutation-verified to bite — disabling it reddens exactly
   one gate (measured at the change, 2026-08-14; the unmutated control
   leaves 55 passing).

   ⛔ This conjunct read :math:`\deg(Q) \ge d` until 2026-08-14, and it
   was never *checked*. The stage only **inverted**: it asked each spec
   for parameters and trusted the answer, so nothing in the selector
   ever compared a built rule against the target at all. Both halves of
   that are now repaired, and they are separate defects — a bare integer
   is not a claim (it omits the measure that gives it meaning), and an
   unverified request is not a verification (a spec that over-promises
   was believed). The witness for the first half is
   ``test_gauss_chebyshev_clears_every_stage_except_the_reference``: on
   a slab query, Gauss-Chebyshev at :math:`n = 4` is **admitted by
   stages 0, 1, 2-as-degree and 3** and refused by the reference alone,
   so without this conjunct list position decides whether it wins.

   ⚠ **Level-symmetric has a top end too, and it is** :math:`S_{18}`:
   from :math:`S_{20}` the per-orbit weight solve has no *positive*
   solution **on this node seed**, so the rule does not exist and its
   inverter returns ``None`` — the same "cannot reach the target with
   any supported parameters" channel Lebedev uses past order 47. The
   frontier belongs to the seed, not to the level-symmetric shape,
   which is why it has moved once already and may move again. The
   measurements, and why positivity is not tradeable, are at
   :ref:`quadrature-ls-positivity` and are deliberately not restated
   here: the frontier is **discovered by attempting the construction**
   rather than compared against a constant, and a second copy of the
   limit — in this page, or in the inverter — is exactly the thing that
   drifts away from the node set. The ⛔ below is that drift, measured.

   ⛔ This bullet gave the level-symmetric degree as
   :math:`\max(3,\, N - 1)`, and the paragraph above it put the top end
   at :math:`S_{12}` with ``[M]`` a minimum weight of :math:`-0.027` at
   :math:`S_{14}`, until 2026-08-14. Both were true of the pre-#337
   node seed :math:`\mu_1^2 = 4/(N(N+2))` and were falsified on
   2026-08-08, when **#337** replaced it with the moment-matched root.
   ``[M]`` the measured degree is now :math:`3, 5, 7, 9, 11, 11, 15,
   15, 17` at :math:`S_2 \ldots S_{18}` — no clean formula in
   :math:`N` — so :math:`\max(3, N-1)` *under*-claims by 2 wherever it
   is wrong; and ``[M]`` the smallest weight at :math:`S_{14}` is
   :math:`+0.01299` — positive — with the first negative one appearing
   at :math:`S_{20}`. Neither number was ever mismeasured. Their
   **configuration** stopped being the shipped one, which is the
   failure mode a copied measurement has and a pointer to the producing
   gate does not.

   ⚠ **The retired formula is a spot-check trap**, which is worth
   knowing before reaching for one: :math:`\max(3, N-1)` happens to be
   *right* at :math:`S_2`, :math:`S_{12}`, :math:`S_{16}` and
   :math:`S_{18}`, and wrong at :math:`S_4`, :math:`S_6`, :math:`S_8`,
   :math:`S_{10}` and :math:`S_{14}`. Four of the nine buildable orders
   confirm it — including :math:`S_{12}`, the order the retired
   frontier made salient. Checking one order is not checking a formula;
   the gate that decides this sweeps every order against the
   closed-form monomial integral
   (``tests/numerics/test_advertised_degree_is_measured.py``).

   ⛔ Earlier still, the level-symmetric entry read :math:`N - 1`
   **(conservative)** until 2026-08-06.  It was neither: ``[M]`` the
   realized degree was **3 at every order**, so the number over-claimed
   from :math:`S_6` up (by 12 at :math:`S_{16}`) *and* under-claimed at
   :math:`S_2` — the signature of a formula describing a different
   construction rather than a cautious bound on this one.  Issue
   **#327** solved the weights per :math:`O_h` orbit; the degree became
   measured against the closed-form monomial integral, and the word
   "conservative" was retired because the figure is exact.

3. **Structural compatibility.** The consumer can request boolean
   flags that the rule must satisfy:

   - ``positive_weights`` — physics consumers reject negative-weight
     rules (Newton-Cotes :math:`n \ge 9`) because they amplify
     quadrature noise.
   - ``axis_aligned`` — Cartesian-friendly: nodes lie on coordinate
     axes (or include axis-aligned orbits, in the Lebedev case).
   - ``level_structured`` — exposes per-:math:`\mu` polar levels for
     the cylindrical SN sweep's azimuthal redistribution coefficients
     (Bailey-Morel-Chang 2010 Eq. (50)), which march per level in
     the order its Eq. (52) fixes.  ⛔ Credited to a **non-existent**
     "Bailey et al. 2009" *Annals of Nuclear Energy* record until
     2026-08-27; the equation numbers were right, the record was not
     (:ref:`bmc-equation-map`).
   - ``half_range_clean`` — ``measure.restrict(lambda x: x > 0)``
     gives a valid quadrature without re-normalisation.

   Only flags passed with value ``True`` constrain the search; ``False``
   and missing keys are interpreted as "don't care."

4. **Minimum points (cost).** Among candidates passing 0-3, pick the
   smallest :math:`N`. Each ordinate is a sweep direction in the SN
   solver, and sweeps dominate the runtime — the cheapest valid rule
   wins. The sort is stable, so an exact tie in node count is broken by
   position in :data:`~orpheus.numerics.quadrature.registry.quadrature_registry`.
   The count is read off the **built** measure that stages 0 and 1
   already forced into existence, never from a parallel formula that
   could drift from it.

Formally, the selection criterion is

.. math::
   :label: quadrature-selection-criterion

   Q^{\star} \;=\; \arg\min\Bigl\{\, n(Q) \;:\;\;
   \mathcal{D}_Q = S^2 / G^0_{\text{geom}}
   \;\wedge\; \Gamma_{\text{geom}} \subseteq \operatorname{Sym}(Q)
   \;\wedge\; \mathcal{E}(Q) \succeq (\lambda_{\text{geom}},\, d)
   \;\wedge\; F_{\text{req}} \subseteq F_Q
   \,\Bigr\},

.. (vv-status rationale) quadrature-selection-criterion: Verified
   transitively by the foundation tests in
   :file:`tests/numerics/test_registry.py` — every stage of the
   five-stage filter has a happy-path test
   (``test_select_slab_returns_gauss_legendre``,
   ``test_select_sphere_returns_gauss_legendre``,
   ``test_select_cylinder_with_level_structured_returns_product``,
   ``test_select_cartesian2d_prefers_lebedev_over_ls_sn``); the two
   conjuncts added on 2026-08-02 carry their own gates
   (``test_support_is_derived_from_the_spent_group_not_declared`` for
   the derived domain,
   ``test_the_two_stages_are_independent_and_both_load_bearing`` for
   their independence,
   ``test_owed_symmetry_selects_by_azimuthal_parity`` and
   ``test_select_cylinder_rejects_odd_azimuthal_product_rule`` for what
   the owed group discriminates, and
   ``test_selector_asks_the_nodes_not_the_declared_tag`` for
   computed-not-declared); and there is a negative path
   (``test_no_rule_fits_raises_with_log``,
   ``test_truly_incompatible_flags_raises``). Tagged ``foundation``
   rather than carrying a verification ladder slot because the
   selection chain is a software invariant (the predicate
   :math:`\mathcal{D}_Q = S^2/G^0_{\text{geom}} \wedge
   \Gamma_{\text{geom}} \subseteq \operatorname{Sym}(Q) \wedge
   \mathcal{E}(Q) \succeq (\lambda_{\text{geom}}, d) \wedge
   F_{\text{req}} \subseteq F_Q`), not a physics-equation claim with
   an L0..L3 ladder slot — the ladder lives on the rules themselves
   (``test_rules_*.py``). The exactness conjunct's own three rejection
   arms carry a gate each
   (``test_a_rule_with_no_exactness_claim_at_all_is_refused``,
   ``test_a_rule_exact_against_the_wrong_measure_is_refused``,
   ``test_an_inversion_that_over_promises_is_caught_not_trusted``),
   with ``test_gauss_chebyshev_clears_every_stage_except_the_reference``
   as the independence witness for the reference half and
   ``test_every_registered_rule_speaks_one_of_the_two_reference_measures``
   pinning the shipped registry against the two derived references.
.. vv-status: quadrature-selection-criterion documented

where :math:`n(Q)` is the number of nodes,
:math:`\mathcal{D}_Q` is the domain the rule's nodes live on,
:math:`G^0_{\text{geom}}` and :math:`\Gamma_{\text{geom}}` are the
continuous and discrete halves of the geometry's angular symmetry
(:class:`~orpheus.numerics.quadrature.registry.AngularSymmetry`),
:math:`\operatorname{Sym}(Q) \subseteq O(3)` is the group the rule's
nodes are **computed** to be invariant under,
:math:`\mathcal{E}(Q)` is the rule's exactness **claim**
(:class:`~orpheus.numerics.exactness.ExactnessClaim` — a reference
measure together with a degree against it),
:math:`\lambda_{\text{geom}}` is the measure the geometry integrates
against, and
:math:`F_Q \subseteq \{\text{positive\_weights}, \text{axis\_aligned},
\text{level\_structured}, \text{half\_range\_clean}\}` is the rule's
structural-flag set. **Domination** is

.. math::

   \mathcal{E}(Q) \succeq (\lambda, d)
   \quad\Longleftrightarrow\quad
   \operatorname{ref}\mathcal{E}(Q) = \lambda
   \;\wedge\;
   \deg\mathcal{E}(Q) \ge d ,

a **partial** order, not a total one: two claims against different
references are simply **incomparable**, which is the whole point.
:math:`\deg = 7` against ``chebyshev_t`` is neither better nor worse
than :math:`\deg = 7` against ``legendre`` — it is an answer to a
different question, and a selector that ranked them would be ranking
answers to two questions on one scale.

.. admonition:: The equation stated the retired predicate for twelve
                days after the prose around it had been corrected
   :class: warning

   Until 2026-08-14 the body of :eq:`quadrature-selection-criterion`
   read

   .. math::

      Q^{\star} \;=\; \arg\min\Bigl\{\, n(Q) \;:\;\;
      G_{\text{geom}} \subseteq G_Q
      \;\wedge\; \deg(Q) \ge d
      \;\wedge\; F_{\text{req}} \subseteq F_Q
      \,\Bigr\},

   that is, three conjuncts instead of four, with the domain conjunct
   absent entirely and the symmetry conjunct in its retired
   declared-tag form. The 2026-08-02 change that split spent from owed
   rewrote the geometry table, the worked examples, the rejection
   messages **and the predicate quoted inside this equation's own
   vv-status rationale**, and left the labelled equation alone.

   The V conjunct then moved later the same day, from
   :math:`\deg(Q) \ge d` to the claim-domination form above, for an
   unrelated reason — it was a bare integer compared across measures
   that give it different meanings, and it was never checked at all.
   That history is in the ⛔ under stage 2; the two corrections are
   separate events that happen to share a date.

   **The tell was legible without reading any code.** The "where" list
   directly beneath the equation defined :math:`\mathcal{D}_Q` and
   :math:`\operatorname{Sym}(Q)` — two symbols that did not occur in
   the equation it annotated — and did not define :math:`G_Q`, which
   did. A definition list that does not match its own equation is
   reporting a correction that stopped one line short.

   **No build could have caught it.** The label existed, so every
   :eq:`quadrature-selection-criterion` reference resolved and Sphinx
   was silent at every severity, nitpicky mode included; the V&V matrix
   listed the label as covered, because coverage is recorded against
   the *label*, not against what the label says. A labelled equation is
   an API — correcting the prose around it does not correct it, and the
   equation is where a reader in a hurry actually looks.

Spent and owed: why the domain and the symmetry are two stages
--------------------------------------------------------------

A geometry's angular symmetry group :math:`G \subseteq O(3)` does not
act on the angular variable as one undifferentiated thing. It splits by
**how the action is used**, and the two halves place two different
demands on a quadrature
(:class:`~orpheus.numerics.quadrature.registry.AngularSymmetry`):

* :math:`G^0`, the continuous part, is **spent** by the dimensional
  reduction. A slab is invariant under every rotation about :math:`z`,
  so :math:`\psi` depends on :math:`\Omega` only through
  :math:`\mu = \Omega\cdot\hat z`; the azimuth is integrated out
  analytically and never discretised. What this half determines is the
  *domain* — the angular variable lives on the quotient
  :math:`S^2/G^0`. In curvilinear geometry the spent half is not free:
  its non-trivial fiber action reappears in the sweep as the
  angular-redistribution (:math:`\alpha`) term, which is where the
  reduction is paid for.
* :math:`\Gamma = G/G^0`, the finite residual, is still **owed**. It
  cannot be integrated away, so a quadrature must realize it as a
  permutation of the ordinates. This is the half a reflecting boundary
  condition consumes.

The two demands are of different **kinds**, and that is the structural
reason one stage cannot carry both. Stage 0 is a question about the
*carrier*: is this node set even a set of points of the right space? It
compares a tag on the measure against the derived quotient, and its
answer is a type, not a tolerance. Stage 1 is a question about the
*arrangement of the nodes within that carrier*: given that they are
points of the right space, is the set closed under a finite group
action? Its answer is a permutation-existence check
(:eq:`discrete-measure-g-invariance`) run against the nodes. Folding
the two into one predicate over one group would require comparing
:math:`G^0`, which is continuous, against a finite point set — and that
comparison has no true instances, which is exactly the shape the
retired gate had.

**The domain is derived, never stored — and so is the measure on it.**
:attr:`AngularSymmetry.support
<orpheus.numerics.quadrature.registry.AngularSymmetry.support>`
computes :math:`S^2/G^0` from the spent group rather than holding a
second, independent column, so there is no state in which the domain
and the spent group disagree — they are one fact. An unmapped quotient
raises :exc:`NotImplementedError` rather than guessing, because a wrong
domain answer silently admits a rule of the wrong dimensionality. Both
halves are pinned by
``test_support_is_derived_from_the_spent_group_not_declared``.

:attr:`AngularSymmetry.reference
<orpheus.numerics.quadrature.registry.AngularSymmetry.reference>` is
the same construction one level down, and it answers the question the
support cannot: **the support says which space, the reference says
which measure on it.** ``[M]`` an **axial rotation** spent
:math:`\to` ``legendre``; nothing spent :math:`\to` ``uniform(S^2)``;
anything else raises :exc:`NotImplementedError`. ⭐ Note the predicate:
since 2026-09-01 the reference arm keys on
:attr:`rotation_axis
<orpheus.numerics.symmetry.SubgroupOfO3.rotation_axis>` being non-``None``,
not on equality with one group — and it is *right* to, whatever the
axis, because the pushforward of :math:`d\Omega` under
:math:`\mu = \hat\Omega\cdot\hat e_a` is :math:`2\pi\,d\mu` (Archimedes'
hat-box) for **every** axis. The measure on the chart does not know
which pole produced it; only the space does. Both are Lebesgue measure
on the domain the reduction left, because transport integrates
:math:`\phi = \int \psi \, d\Omega` unweighted — and ``legendre``
**is** Lebesgue on :math:`[-1,1]`: its weight is :math:`w(x) = 1` and
``[M]`` its mass is exactly :math:`2`. The name records the polynomial
family the measure generates, not a weighting.

A rule can get the space right and the measure wrong — that is exactly
the gap stage 2's reference conjunct exists to close, and the
Gauss-Chebyshev witness above is what it looks like. ``[M]`` on the
shipped registry the two facts happen to travel together
(``GaussLegendre1D`` claims ``legendre``; ``LebedevSphere``,
``LevelSymmetricSN`` and ``ProductQuadrature`` all claim
``uniform(S^2)``), which is *why* the witness has to be constructed
rather than found — and why a conjunct with no shipped counter-example
still earns its place.

**Neither conjunct implies the other**, measured on the shipped
registry — drop either stage and the gate admits something wrong:

.. list-table::
   :header-rows: 1
   :widths: 22 14 14 10 10 30

   * - Candidate
     - Geometry
     - Nodes live on
     - Stage 0
     - Stage 1
     - What the missing stage would have let through
   * - ``lebedev_sphere(order=5)``
     - ``"slab"``
     - :math:`S^2`
     - ✗
     - ✓
     - An :math:`S^2` cubature serving a :math:`\mu`-marginal: being
       :math:`O_h`-invariant, it is certainly
       :math:`\sigma_x`-closed.
   * - ``gauss_legendre_on_polar_orbit(n=8, axis="x")``
     - ``"cylinder"``
     - :math:`S^2/SO(2)_x`
     - ✗
     - ✓
     - A polar marginal serving a geometry whose two angular degrees
       of freedom both survive.
   * - ``gauss_legendre_on_mu(n=8)``
     - ``"slab"``
     - :math:`[-1,1]`
     - ✗
     - ✓
     - ⭐ **New row, 2026-09-01.** The *chart-level* rule, refused by
       its own geometry: it names an interval, and the slab discretises
       an orbit space. This is why the registry's ``GaussLegendre1D``
       spec is ``partial(gauss_legendre_on_polar_orbit, axis="x")`` and
       the chart rule is deliberately unregistered
       (:ref:`manifold-polar-orbit-rule`).
   * - ``product_mu_phi(3, 5)``
     - ``"cylinder"``
     - :math:`S^2`
     - ✓
     - ✗
     - An odd-:math:`n_\varphi` grid, on which no reflecting
       :math:`x` or :math:`y` face is representable as an ordinate
       permutation (ERR-042, stated structurally).
   * - ``product_mu_phi(3, 6)``
     - ``"cylinder"``
     - :math:`S^2`
     - ✓
     - ✓
     - Admitted — the control row.

`[M]` every cell above was re-measured 2026-09-01. Two of the rows are
the direct content of
``test_the_two_stages_are_independent_and_both_load_bearing`` — the
Lebedev-for-a-slab row (symmetry alone would admit it) and the
odd-:math:`n_\varphi`-for-a-cylinder row (domain alone would) — with
``test_owed_symmetry_selects_by_azimuthal_parity`` covering the parity
pair. They are why the split is not a refactor of one predicate into
two: a single "symmetry" stage is a strictly weaker gate than either of
these columns, admitting the union of both failure sets.

⚠ The two ``gauss_legendre`` rows carry no dedicated independence test
and are listed as *measurements*, not as gate content — the polar-orbit
row because it is the mirror image of the Lebedev one (a 1-D rule
offered to a 2-D geometry), and the chart-level row because the thing
it demonstrates is a **registration** decision rather than a selection
one: that rule is never presented to stage 0, because it is not in the
registry. Both reduce to the same fact, which is that stage 0 compares
*one* datum on both sides — the geometry names the group it spends, the
rule names the group its orbit space was quotiented by.

Why the stages do not evaluate in their own order
--------------------------------------------------

The conjunction in :eq:`quadrature-selection-criterion` is order-free —
a conjunction has no order — but the *evaluation* is not, and it splits
**one conjunct across two moments**:

* The V constraint runs **first, as an inversion**: ask each spec for
  the smallest parameters reaching :math:`d`. Then the rule is built.
  Then stages 0 and 1 interrogate its nodes.
* The V constraint runs **again, as a verification**, *after* stages 0
  and 1 — checking that the claim the built rule actually carries
  dominates :math:`(\lambda_{\text{geom}}, d)`. Then stage 3, then
  stage 4.

The first half is forced: :math:`\mathcal{E}(Q)` is a property of an
*instantiated* rule, stages 0 and 1 are questions about a rule's
*nodes*, and a rule has no nodes until its parameters are fixed — which
only the V constraint does. The order follows from what each predicate
must look at, not from a speed argument.

The second half is a *placement* choice, and it has a reason worth
stating because it looks arbitrary: :math:`\lambda_{\text{geom}}`
**encodes** :math:`\mathcal{D}_Q` — the reference and the support are
both derived from the same spent group — so a rule on the wrong space
fails the reference conjunct *and* the domain conjunct, and the domain
stage gives the better diagnosis ("your nodes live on :math:`S^2`, this
geometry discretises :math:`S^2/SO(2)_x`" beats "your measure is
``uniform(S^2)``, we wanted ``legendre``"). ⭐ The domain half of that
diagnosis got strictly sharper on 2026-09-01: it used to name the chart
:math:`[-1,1]`, which is the *same* interval for all three axes and so
could not distinguish a rule about the wrong pole. Running the
verification
after stages 0 and 1 therefore reserves it for the case it *uniquely*
catches: **a rule on the right space carrying the wrong measure.** That
case has no shipped instance, which is exactly why it needs a
constructed witness and a stated reason rather than a comment.

⛔ This subsection described the V constraint as a single moment — "V
first, then instantiate, then stages 0 and 1" — until 2026-08-14,
because until then the stage genuinely was one moment: it inverted and
never verified.

The deeper reason the symmetry conjunct cannot be hoisted above the
build is that **a rule's invariance group is parameter-dependent**. The
product rule's is :math:`D_{n_\varphi h}` — ``[M]`` :math:`D_{5h}` at
``(n_mu=3, n_phi=5)`` and :math:`D_{6h}` at ``(3, 6)``, which is
precisely why the cylinder refuses the first and admits the second. No
parameter-free field on
:class:`~orpheus.numerics.quadrature.registry.QuadratureSpec` can state
that truthfully, which is why the spec deliberately carries **no**
invariance group, and why the older design that tried to had to lie.
The rule's own measure keeps its ``invariance_group`` metadata as the
single source of that tag; duplicating it on the spec was a twin source
of truth.

The domain conjunct is a weaker case of the same argument — every
shipped rule's ``support`` happens to be parameter-independent, so it
*could* have been a field — but since stage 1 forces the build anyway,
the selector reads every such fact off the one constructed object: the
support, the computed invariance group, the exactness claim, **and the
node count stage 4 sorts on**. The winning candidate's measure is then
carried forward rather than rebuilt. There is therefore exactly one
construction of the measure the selector returns, and no second,
divergence-capable one — the same single-source argument that makes the
spec carry no invariance group and no parallel node-count formula.

The cost of building a candidate before rejecting it is real and
accepted: ``[M]``
:func:`~orpheus.numerics.quadrature.registry.select_quadrature` has no
production consumer today — ``solve_sn`` takes an explicit quadrature,
for the reasons under *Why a registry and not a hardcoded ladder*
below — so selection is not a hot path, and paying a construction to
ask an honest question is the right trade against caching a claim that
can go stale.

Geometry → angular-symmetry assignment
--------------------------------------

The selector's static geometry table
(:data:`~orpheus.numerics.quadrature.registry.GEOMETRY_ANGULAR_SYMMETRY`)
records **both halves** of each geometry's angular symmetry, for the
reasons set out under *Spent and owed* above. One row per supported
geometry; a new geometry is added here, never in the selector itself:

.. list-table::
   :header-rows: 1
   :widths: 18 12 12 12 46

   * - Geometry
     - :math:`G^0` (spent)
     - Domain :math:`S^2/G^0`
     - :math:`\Gamma` (owed)
     - Rationale
   * - ``"slab"``
     - :math:`SO(2)_x`
     - :math:`S^2/SO(2)_x`
     - :math:`\sigma_x`
     - 1-D. Azimuthal rotation about the streaming axis is integrated
       out, leaving :math:`\mu = \hat\Omega\cdot\hat e_x` alone. Owed is
       :math:`\mu \to -\mu`, which pairs the two sweep senses;
       Gauss-Legendre nodes are symmetric (Stoer-Bulirsch §3.6), so it
       holds. ⭐ Spent and owed now name the **same** axis, and so does
       the real spherical-harmonic pole — one axis, stated on every half.
   * - ``"sphere"``
     - :math:`SO(2)_x`
     - :math:`S^2/SO(2)_x`
     - :math:`\sigma_x`
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

.. note::

   ⛔ **The two 1-D rows read** :math:`SO(2)` **/ domain**
   :math:`[-1,1]` **until 2026-09-01.** Both halves changed at tracker
   2.4 and for one reason. The spent group is now :math:`SO(2)_x`,
   because a group realized on a coordinate axis carries that axis as
   data (:ref:`manifold-so2-axis-is-a-parameter`); and the domain is now
   the **orbit space** rather than its chart, because
   :attr:`AngularSymmetry.support
   <orpheus.numerics.quadrature.registry.AngularSymmetry.support>`
   derives it as ``SPHERE.quotient(spent)`` instead of returning a
   hand-written ``COSINE_INTERVAL``. `[M]` the two are not
   interchangeable: all three :math:`S^2/SO(2)_a` realize onto the
   **same** interval, so taking the realization is exactly the step that
   throws the axis away — and it is the step that let a slab angular
   rule and a spatial rule on :math:`[-1,1]` be the same value.

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
  :math:`S^2` rules fail stage 0, `[M]` verbatim —
  ``domain mismatch: geometry 'slab' discretises S^2/SO2_x, but the
  rule's nodes live on S^2`` — leaving GL1D (n=8 → degree 15, 8 nodes),
  whose own measure `[M]` declares ``support=S^2/SO2_x``.

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
artifact. What the registry carries per rule is the structural-flag
*values* in the four
:class:`~orpheus.numerics.quadrature.registry.QuadratureSpec` entries,
each annotated by a one-line comment at its assignment site; the prose
explaining what the flags *mean* is the ``Attributes`` section of the
single ``QuadratureSpec`` class docstring, and stage 3 above.

⛔ This paragraph said that each ``QuadratureSpec`` docstring narrates
the trade-off in the rule's design space, until 2026-08-14. It promised
a mechanism that cannot exist. ``QuadratureSpec`` is a dataclass, so
its four *instances* share one class docstring: ``[M]``
``spec.__doc__ is type(spec).__doc__`` is ``True`` for all four, and no
instance carries a ``__doc__`` of its own. Nor does anything in
:mod:`orpheus.numerics.quadrature.registry` become a Sphinx table row,
because the module is deliberately not ``automodule``-rendered — it
carries ``.. math:: :label:`` docstrings that would collide with this
page's labels. **This page is the rendered narration**, which is why a
drift like this has to be caught by reading, never by a build.

The promise was minted twice: the module's own docstring stated it more
concretely still, as "every entry's docstring becomes a Sphinx table
row". That duplication is itself the tell — a claim about a *rendering*
mechanism had never been checked against the rendering, and the module
asserting it is the one that is not rendered.

Domain-specific quadrature: the Quadrature class
================================================

The SN solver consumes angular quadratures through ~50
attribute-access sites (``quad.mu_x``, ``quad.weights``,
``quad.ordinate_permutation(motion)``, ``quad.spherical_harmonics(L)``, …)
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
  :class:`DiscreteMeasure` instances tagged with ``invariance_group``
  and ``exactness`` — the latter a whole
  :class:`~orpheus.numerics.exactness.ExactnessClaim` (a reference
  measure *and* a degree against it), of which ``degree_of_exactness``
  is a derived view:

  - :func:`~orpheus.numerics.quadrature.gauss_legendre_on_mu` —
    the **chart-level** 1-D rule on :math:`\mu \in [-1, 1]`,
    :math:`\sigma_x`-invariant, exact to algebraic degree
    :math:`2n - 1` against ``legendre``. It names **no** orbit space,
    on purpose: it is simultaneously the raw material of the slab's
    rule (:math:`\mu` about :math:`x`) and the polar factor of every
    product rule (:math:`\mu` about :math:`z`), so a factor declaring
    one axis inside a product about the other would be a false claim.
  - :func:`~orpheus.numerics.quadrature.gauss_legendre_on_polar_orbit`
    — the same atoms READ on :math:`S^2/SO(2)_a` for a named axis, via
    :meth:`~orpheus.numerics.measure.DiscreteMeasure.on_orbit_space`
    (:ref:`measure-on-orbit-space`). This is what the SN slab solver
    consumes through :meth:`Quadrature.gauss_legendre
    <orpheus.numerics.quadrature.directional.Quadrature.gauss_legendre>`
    and what the registry builds, as
    ``partial(gauss_legendre_on_polar_orbit, axis="x")``. `[M]` its
    nodes and weights are bit-identical to the chart rule's; its
    ``support`` is ``S^2/SO2_x``, its ``space`` ``L2[S^2/SO2_x]``, and
    its ``invariance_group`` is re-stamped ``Mirror(axis)`` for the
    named axis.
  - :func:`~orpheus.numerics.quadrature.lebedev_sphere` —
    Lebedev rule on :math:`S^2`, :math:`O_h`-invariant, exact to
    spherical-harmonic degree ``order`` against ``uniform(S^2)``.
  - :func:`~orpheus.numerics.quadrature.level_symmetric_sn` —
    Carlson-Lathrop level-symmetric :math:`S_N` rule on
    :math:`S^2`, :math:`O_h`-invariant, against ``uniform(S^2)``;
    its degree is build-measured, with no formula in :math:`N`.

  - :func:`~orpheus.numerics.quadrature.periodic_trapezoid` — the
    :math:`n`-point rule on the circle :math:`S^1`, exact for
    *trigonometric* polynomials of degree :math:`n - 1`. The
    azimuthal factor of the product rule below.
  - :func:`~orpheus.numerics.quadrature.product_mu_phi` — the
    polar-times-azimuthal product rule used by the cylindrical SN
    sweep.

  .. note::

     ⛔ The Gauss-Legendre entry read ":math:`SO(2)`-invariant" until
     2026-08-14 — the last surviving copy on this page of the claim
     retired on 2026-08-02, that a finite point set can be closed under
     a continuous group. ``[M]``
     ``gauss_legendre_on_mu(8).invariance_group`` is :math:`\sigma_x`,
     the mirror through the plane normal to :math:`x`, which is what
     :math:`\mu \to -\mu` becomes once the polar marginal is embedded
     as :math:`(\mu, 0, 0)`. And the bare integer
     ``degree_of_exactness`` was the whole exactness tag until the
     reference joined it; reading it alone is reading half a claim,
     which is precisely the state that let a Gauss-Legendre and a
     Gauss-Chebyshev degree look interchangeable.

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

.. admonition:: The generation principle, per family (#325)
   :class: note

   **A node set defined by a symmetry group is generated by the group
   action — sign flips, coordinate swaps, index arithmetic, all exact
   on floats — never by evaluating a parameterization of it.**
   Evaluating trig at symmetric angles destroys the symmetry to
   rounding, and no downstream tolerance recovers it.  Where each
   shipped family stands:

   * ``lebedev`` — fully algebraic (:math:`O_h` orbits of
     :math:`\pm(a,b,c)` permutations): mirrors and tangential zeros
     exact by construction.
   * ``periodic_trapezoid`` / ``product`` / ``folded_product`` — the
     azimuthal nodes are roots of unity (octant fold +
     integer-arithmetic fixed points); the symmetry CHECKER's own
     :math:`C_n` operators consume the same generator, so one
     trigonometric spelling exists in the whole chain.
   * MoC azimuthal — the same roots-of-unity call at
     :math:`(2k+1)/(4n)` of a turn; the mirror is the index map
     :math:`k \mapsto n-1-k`.
   * ``level_symmetric`` — a HYBRID by nature: the eight-octant sign
     replication is exact group action (a level's :math:`\eta` values
     are bit-copied), while the orbit representative is
     float-*solved* — since #337 from the :math:`\int\mu_z^N` moment
     condition (before it, from a formula).  The exactness that
     matters (48/48 :math:`O_h` operators land bit-exactly) comes
     from reading every cosine out of the ONE level array, so the
     symmetry is an integer permutation regardless of how the
     representative was produced.

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
  function above, plus the SN-side derived data (octant partition,
  level structure) cached at construction.  Reflection partners are
  no longer precomputed: the ``reflection_index`` table was retired
  (§7d.3), and a motion's ordinate permutation is derived on demand
  by :meth:`Quadrature.ordinate_permutation
  <orpheus.numerics.quadrature.Quadrature.ordinate_permutation>` —
  a certified match (bijection, equal weights, no bare
  nearest-neighbour) returning a typed permutation.  The legacy
  ``mu_x`` / ``mu_y`` / ``mu_z`` / ``weights`` / ``N`` surface
  survives as ``@property`` views over the underlying
  ``measure.nodes`` / ``measure.weights``, so the ~50 consumer sites
  see no API change.

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

.. _discrete-measures-quadrature-axis:

``quad.axis(label)`` — the rule as a space factor
--------------------------------------------------

Beside the ordinate accessors sits one mint, added at campaign-1 phase
CS5 (2026-08-29):
:meth:`Quadrature.axis(label)
<orpheus.numerics.quadrature.directional.Quadrature.axis>` returns the
angular :class:`~orpheus.numerics.axis.Axis` — the space factor whose
measure is this rule's weights — carrying **this rule** as its
:attr:`~orpheus.numerics.axis.Axis.generator`.

Three properties of the accessor are decisions rather than details, and
each one is load-bearing:

#. **It delegates the structural mint and upgrades only the
   provenance.** The shape/weights/``kind`` logic lives once, in
   :meth:`DiscreteMeasure.axis
   <orpheus.numerics.measure.DiscreteMeasure.axis>`; ``Quadrature.axis``
   calls it and then re-points the generator with
   ``dataclasses.replace``. There is no second copy of the mint to drift
   — the Pattern-2 single-source rule applied to a two-tier accessor.
   Because ``replace`` re-runs the axis's ``__post_init__``, the upgrade
   cannot bypass a construction invariant: the all-ones collapse, the
   :math:`-0.0` normalization and the read-only weights flag all survive
   it, and that survival is asserted rather than assumed.
#. **The generator is the rule, not the wrapped measure** — and that
   is the opposite of what the layering instinct suggests. The measure is the lower object, but it is also the
   quadrature *minus its side-channels*: ``level_indices`` is carried by
   :class:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure`,
   which is not one of the measure's five fields. A measure-typed
   generator would answer three of the four names a curvilinear consumer
   needs. The refutation is developed at
   :ref:`spaces-generator-why-quadrature`.
#. ``kind`` **is not a parameter.** A quadrature's components are point
   values on the sphere, so the generator's type implies ``NODAL`` and
   the mis-declaration "a nodal basis carrying no nodes" is unspellable
   on this path.

The consumer-facing consequence is that a caller holding an
axis-composed space no longer needs the rule handed to it separately:
``space.axis("angular").generator`` answers ``mu_x`` / ``eta`` /
``mu_z`` / ``level_indices`` and the rest of the ordinate surface. ⚠ The
generator is a **live reference** — a
:class:`~orpheus.numerics.quadrature.directional.Quadrature` is a
mutable dataclass — whereas the axis's own ``weights`` are a read-only
defensive copy. Read the generator for what the axis forgot; read the
axis for the factor measure of record
(:ref:`spaces-axis-generator`).

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
* Bailey, T.S., Morel, J.E. and Chang, J.H. (2010). "The Asymptotic
  Diffusion-Limit Accuracy of :math:`S_N` Angular Differencing
  Schemes." *Nuclear Science and Engineering* **165**, no. 2, 149-169,
  doi:10.13182/NSE08-66.  **Eq. (50)** is the cylindrical (R-Z)
  :math:`\alpha`-recursion that ``level_structured`` exists to feed;
  **Eq. (52)** fixes the order of a level.  ⛔ Both were credited to a
  **non-existent** "Bailey et al. 2009", *Annals of Nuclear Energy*
  **35**, 1929-1936, until 2026-08-27 — the equation numbers were
  right, the record was not (:ref:`bmc-equation-map`).
