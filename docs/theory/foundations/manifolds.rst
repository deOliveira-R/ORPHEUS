.. _manifolds:

=================================================================
Manifolds: the Point Set, the Orbit Space, and What a Basis Eats
=================================================================

.. contents:: Contents
   :local:
   :depth: 2


.. Machine header — the ``nexus-meta`` schema for this page (PROVISIONAL).
.. Seeded 2026-08-31 at tracker 2.0a of the angular-spaces campaign (#429),
.. under user ruling D0.7. This page owns LEVEL 1 of the three-level stack;
.. ``spaces`` owns levels 2 and 3 and ``discrete_measures`` owns the measure
.. that lives on level 1. The schema is provisional pending a full re-audit
.. of the corpus.

.. dropdown:: Machine header — ``nexus-meta`` schema (PROVISIONAL)
   :color: muted

   .. code-block:: yaml

      module: numerics
      concept: manifolds
      role: "the point-set layer — the manifold M a measure is supported on and a basis function is defined over, its algebra (product, orbit space, membership), the invariant-theoretic derivation that produces an orbit space, the TWO coordinate systems an orbit space honestly has (the invariant chart's codomain and a section's image), and the three-level separation (manifold / fields on it / coefficients) that keeps a FunctionSpace from being mistaken for a domain"
      depends_on: []
      related: [discrete_measures, spaces, frame, spherical_harmonics]
      status: "MINTED, gated, WIRED, and CONSUMED. Two catalogued derivations ship (S^2/SO(2)_a for the three ROTATION axes and S^2/<sigma_a> for the three MIRROR axes — six keys, two procedures — plus the derived identity quotient), and a Quotient carries BOTH coordinate systems after the 2026-08-31 two-slot ruling. `Space = str` and its six SPACE_* tags are RETIRED (tracker 2.0c, 2026-09-01): `DiscreteMeasure.support`, `GeneratingMeasure.support`, `UniformMeasure.support`, `ProductMeasure.support`, the `ReferenceMeasure` Protocol and `AngularSymmetry.support` all carry a Manifold, and `Basis.domain` does too (2.1). Tracker 2.1b (2026-09-01) read a SECOND answer off that same slot: `Basis.invariance_group` is DERIVED from `domain` by a match on its TYPE (a Quotient of the sphere -> its `by`; the sphere -> Trivial; anything else -> None), so a basis declares the symmetry its functions HAVE by naming the manifold they EAT. `[M]` 6 of 6 shipped bases answer, the property is @final, and it cost zero subclass edits and no new field. ERR-080's pairing therefore has BOTH operands and is a computable lattice verdict — `[M]` `Trivial contains SO2('x')` is False for the slab, while the shipped fold's two halves are literally ONE group object. Nothing CONSUMES that verdict yet: the frame's pairing gate is tracker 2.2, and a gate written on the FRAME's measure would be inert today because that measure still carries the forged S^2. Tracker 2.4 (2026-09-01) gave the axial rotation group its AXIS — `SO2(axis)` beside `Mirror(axis)` — and made the slab's polar rule DECLARE its orbit space: `Quadrature.gauss_legendre(8).measure.support.name == 'S^2/SO2_x'`, via the new verb `DiscreteMeasure.on_orbit_space`. That is this page's first PRODUCTION consumer, and it collapsed the registry twin (`AngularSymmetry.support` now calls `SPHERE.quotient`). Tracker 2.3 (2026-09-02) gave the category its ARROWS: `ManifoldMap(domain, codomain, apply)` is a frozen value type, composition `psi @ phi` is refused across mismatched endpoints, and `DiscreteMeasure.pushforward` now READS its target off `phi.codomain` (`new_space=` retired) and refuses a map out of the wrong point set — by manifold VALUE, so the slab's `S^2/SO2_x` rule and the chart rule on `[-1,1]`, whose nodes are `np.array_equal`, are told apart. Three arrows are typed: `archimedes(axis)` ([-1,1] x S^1 -> S^2, Archimedes' hat-box, `[M]` the product rule is bit-identical to its retired hand loop on 60 of 60 configurations and its support IS the chart's codomain); the orbit retraction inside `quotient()`; and `barycentre(orbit_space)` (S^2/SO(2)_a -> Ball(3), since 1 - norm(mu e_a)^2 = 1 - mu^2 = det P / 4). ERR-080 restated in that vocabulary: it is the barycentre map with a FORGED codomain — `[M]` the forgery's nodes are np.array_equal to the honest map's image and differ only in the type claimed. 2.3 is an ENABLER and repairs nothing: no membership check runs inside a map (that refusal is tracker 2.0b, at measure construction), the forgery arm stays a raw constructor BY DESIGN until 3.4, and `[M]` the ERR-080 gate still declares three xfail(strict=True) rows. Neither the entry's chart nor its section ships — `[M]` `Quotient.fundamental_domain` still has zero readers outside this module, and the pushforward REFERENCE measure is deferred to 3.1 because `[M]` a module-scope manifold -> exactness import closes a two-hop cycle that kills 5 of 5 fresh import orders. ERR-080 itself is still OPEN"


This page develops the **point-set layer** — the thing a measure is
supported *on*, the thing a basis function is defined *over*, and the
thing a quotient by a symmetry group produces. It is the layer *below*
:doc:`/theory/foundations/spaces` (which types the vector space a
discrete field lives in) and *below*
:doc:`/theory/foundations/discrete_measures` (which types the weighted
atom list). Until 2026-08 the corpus had no object for it at all: level
1 was a ``str``.

The organizing claim is one sentence: **a manifold is a first-class
value with an algebra — product, orbit space, membership — and it is a
different object from the function space of fields defined on it.**
Everything on this page follows from taking that seriously: why a
1-dimensional quadrature declaring its nodes to live on :math:`S^2` is
a forgery that a type can refuse, why ``f"{a} × {b}"`` was already
performing a product, why an orbit space is an **orbifold** and not a
quotient manifold, and why :math:`\det P = 4(1-\mu^2)` is the same
polynomial three times over.

.. warning::

   **"Manifold" is the second thing in this corpus with that name, and
   the older one is unrelated.** `[M]` counted on the tree as it stood
   before this page: **21** occurrences of the word in
   ``docs/theory/``, of which **14** are the
   S\ :sub:`N` *solution manifold* — the affine set
   :math:`\psi^\star + \ker A` that a singular loss operator admits
   instead of a unique solution
   (:doc:`/theory/foundations/field_algebra`,
   :doc:`/theory/methods/sn/cartesian_multid`). That object is a coset
   in a *vector* space, reached by a splitting's gauge freedom. *This*
   page's manifold is a point set in **Man**, reached by a chart. They
   share no machinery and no type; only the word.

   The discriminator when reading: a *solution* manifold is always
   qualified by the word "solution" or by :math:`\ker A`, and it is
   never an argument to anything. A :class:`Manifold
   <orpheus.numerics.manifold.Manifold>` is always the thing something
   else is *defined over*.

.. note::

   **Three words for the same level-1 object, all standard, all kept.**
   A measure calls it its **support**
   (:attr:`DiscreteMeasure.support
   <orpheus.numerics.measure.DiscreteMeasure.support>`); a basis
   function calls it its **domain**; the quadrature registry calls it
   the angular **domain** :math:`S^2/G^0`
   (:attr:`AngularSymmetry.support
   <orpheus.numerics.quadrature.registry.AngularSymmetry.support>`).
   These are not twins — *support of a measure* and *domain of a
   function* are both correct usage for the same manifold, and
   category-theoretically ``dom`` is the source of a morphism in both
   cases (in **Man** for a basis function, in **Vect** for an
   operator). What was missing was not a word; it was the object all
   three words name.

   ⚠ ``support`` is nonetheless a **misnomer of a different kind**, and
   the corpus already says so: `[M]`
   ``gauss_legendre(8).measure.support`` is ``'[-1,1]'`` while
   :math:`\operatorname{supp}(\mu)` is 8 points. The tag names the
   ambient manifold, not the support of the measure. Renaming it is
   tracked with the migration (:ref:`manifold-seams`); this page uses
   *manifold* for the object and quotes ``support`` when it means the
   slot.

.. admonition:: Key Facts (the point-set layer)
   :class: tip

   - **Three levels, and the tree named two.** The manifold :math:`M`;
     the fields on it :math:`L^2(M)`, which is what a
     :class:`~orpheus.numerics.space.FunctionSpace` is; and the
     coefficients :math:`\mathbb{R}^K`, which is what a
     :class:`~orpheus.numerics.basis.base.Basis`'s shape is. A basis
     function **eats a point of** :math:`M` — that is why a
     ``FunctionSpace`` cannot be a basis's domain
     (:ref:`manifold-three-levels`).
   - **The level-2 check passes; the level-1 check had nothing to
     check.** That is why :ref:`ERR-080 <manifold-err-080>` survived.
     On a slab the frame's arrow ``measure.space → basis.space``
     is between two well-formed spaces — `[M]` ``L2[S^2]`` of shape
     :math:`(8,)` into ``spherical_harmonic_space`` of shape
     :math:`(3,5)` — while the nodes it carries are not points of the
     manifold the basis needs.
   - **Membership is what makes a support claim falsifiable.**
     :meth:`Manifold.contains
     <orpheus.numerics.manifold.Manifold.contains>` is the defining
     equation of the member evaluated on the candidate nodes. `[M]` on
     ``gauss_legendre(8).angular_frame(2)`` the production measure
     declares ``support='S^2'`` over rows whose norms are
     :math:`0.1834 \ldots 0.9603` — **0 of 8** within :math:`10^{-12}`
     of 1 — and :class:`Sphere <orpheus.numerics.manifold.Sphere>`
     refuses them, three hops upstream of the wrong answer
     (:ref:`manifold-err-080`).
   - ✅ **The type is WIRED** (tracker 2.0c, 2026-09-01): ``Space = str``
     and its six ``SPACE_*`` tags are retired, and every measure in the
     tree carries a ``Manifold``. ⛔ **The membership PREDICATE is still
     not enforced at construction**, which is the half ERR-080 needs:
     :meth:`contains` ships and is gated, but nothing calls it on the way
     in, so the forged :math:`(\mu, 0, 0)` measure is still
     *constructible*. ERR-080 is **open** — its refusal is tracker 2.0b
     plus the fused fix step (:ref:`manifold-seams`).
   - ✅ **…and CONSUMED** (tracker 2.4, 2026-09-01): the slab's polar
     quadrature now *declares* the orbit space it lives on. `[M]`
     ``Quadrature.gauss_legendre(8).measure.support.name`` is
     ``'S^2/SO2_x'`` and its induced space is ``'L2[S^2/SO2_x]'``, with
     nodes and weights **bit-identical** to the chart-level rule
     (``np.array_equal`` on both). The declaration is a repair, not
     wiring: `[M]` before it, an 8-node slab **angular** space and an
     8-node **spatial** rule on :math:`[-1,1]` compared ``==`` *and*
     hash-equal (:ref:`manifold-orbit-space-declaration`).
   - ⭐ **A basis declares the symmetry its functions HAVE by naming the
     manifold they EAT** (tracker 2.1b, 2026-09-01). A function on
     :math:`M/H`, pulled back to :math:`M`, *is* an
     :math:`H`-invariant function — so for FUNCTIONS the group a basis
     HAS and the group its domain SPENT are one property, and
     :attr:`Basis.invariance_group
     <orpheus.numerics.basis.base.Basis.invariance_group>` is a
     ``match`` on ``domain``: stored nowhere, ``@final``, `[M]` **6 of
     6** shipped bases answering with **0** subclass edits. For a POINT
     SET they come apart, which is why a measure needs **two** slots:
     `[M]` the slab's rule HAS ``Mirror('x')`` while it SPENT
     ``SO2('x')``, and the :math:`\sigma_y` fold HAS **nothing**,
     because folding destroys the closure it spends
     (:ref:`manifold-has-versus-spent`).
   - ⛔ **The pairing is now a MEASUREMENT, and still not a refusal.**
     Both operands of the ERR-080 pairing now exist — `[M]`
     ``Trivial ⊇ SO2('x')`` is **False** for the slab's rule against the
     full-sphere harmonics its frame binds today, and ``Mirror('y') ⊇
     Mirror('y')`` is **True** for the shipped fold, where the two
     halves are the ``by`` of **one** memoised :class:`Quotient`
     (``has is spent``). Nothing consumes the verdict: the frame's
     pairing gate is tracker 2.2, and `[M]` the slab frame's own measure
     still carries the forged ``S^2``, so a gate written there would be
     inert today (:ref:`manifold-invariance-pairing`).
   - ⭐ **The rotation axis is a PARAMETER, not a convention — because
     the tree carries two poles.** :math:`SO(2)` left the parameter-free
     enum on 2026-09-01 and became ``SO2(axis)``, exactly as the
     reflection had on 2026-08-02. `[M]` the real spherical-harmonic
     basis takes ``cos θ = μ_x`` while a product rule's polar factor is
     :math:`\mu_z`, and **one** Gauss–Legendre rule serves both roles —
     so the group a marginal was quotiented by cannot be spelled without
     naming its axis (:ref:`manifold-so2-axis-is-a-parameter`).
   - ⭐ **A map carries its two point sets, so a codomain cannot be
     forged at the call site** (tracker 2.3, 2026-09-02).
     :class:`~orpheus.numerics.manifold.ManifoldMap` gives the
     category its **arrows**; ``pushforward`` reads its target off
     ``φ.codomain`` (``new_space=`` retired) and refuses a map whose
     ``domain`` is not the measure's support — by manifold VALUE, so
     `[M]` the slab's rule and the chart rule, whose nodes are
     ``np.array_equal``, are told apart. Three arrows ship: the
     Archimedes chart, whose codomain the product rule now READS
     (`[M]` bit-identical to the retired hand loop on **60 of 60**
     configurations); the orbit retraction inside ``quotient()``; and
     the orbit **barycentre**
     :math:`\mu \mapsto \mu\,\hat e_a`, whose honest codomain is
     ``Ball(3)`` because :math:`1-\lVert\mu\hat e_a\rVert^2 = 1-\mu^2
     = \tfrac14\det P`. ⟹ **ERR-080 is that map with a forged
     codomain** — `[M]` the forgery's nodes are ``np.array_equal`` to
     the honest map's image and differ only in the type claimed
     (:ref:`manifold-arrows`).
   - ⛔ **2.3 is an ENABLER: it repairs nothing.** No membership check
     runs inside a map (that refusal is tracker 2.0b, at measure
     construction); the forgery arm stays a raw constructor **by
     design**, because routing it through ``pushforward`` would force
     it to tell the truth; and `[M]` the ERR-080 gate still declares
     **three** ``xfail(strict=True)`` rows, untouched. Neither the
     entry's chart nor its section ships — `[M]`
     ``fundamental_domain`` still has zero readers outside the module
     (:ref:`manifold-arrows-not-built`).
   - **The algebra was already running, spelled as string
     interpolation.** `[M]` ``measure.py:588`` was
     :meth:`__mul__ <orpheus.numerics.manifold.Manifold.__mul__>`
     (``f"{self.support} × {other.support}"``), ``:1022`` was
     :meth:`quotient <orpheus.numerics.manifold.Manifold.quotient>`
     (``f"{self.support}/{group.name}"``) and ``:802`` was the
     pushforward's codomain (``f"φ_*({self.support})"``). The
     interpolation *was* the operation, unnamed
     (:ref:`manifold-string-tag`).
   - **The type is a CLOSED SUM split by TOTALITY**, not a polymorphic
     hierarchy. ``dim`` / ``name`` / ``contains`` / ``__mul__`` are
     answerable by every manifold and live on the abstract base; the
     derivation fields are answerable only by a quotient and live on
     :class:`~orpheus.numerics.manifold.Quotient` alone, so asking a
     sphere for a syzygy ideal is a type error rather than a ``None``
     (:ref:`manifold-closed-sum`).
   - **An orbit space is derived, by invariant theory, once per
     entry** — minimal invariants of :math:`\mathbb{R}[x]^H`, the
     syzygy ideal by elimination, then the Procesi–Schwarz condition
     :math:`P \succeq 0` with
     :math:`P_{ij} = \langle \nabla p_i, \nabla p_j\rangle`
     (:eq:`manifold-procesi-schwarz`). For :math:`S^2/SO(2)_a` this gives
     :math:`P = \operatorname{diag}(1, 4p_2)`,
     :math:`\det P = 4p_2 = 4(1-\mu^2)` and the orbit space
     :math:`[-1,1]` (:eq:`manifold-s2-mod-so2`); for the shipped
     cylindrical fold :math:`S^2/\langle\sigma_a\rangle` it gives
     :math:`P = \operatorname{diag}(1,1,4p_3)` and the **closed unit
     disk** :math:`D^2` (:eq:`manifold-s2-mod-mirror`). **Both families
     are one derivation reading the axis off the group**, so the
     catalogue holds `[M]` **six** keys served by **two** procedures.
   - ⭐ **An orbit space has TWO honest coordinate systems, and a
     Quotient carries both** (user ruling, 2026-08-31). ``realization``
     is the **invariant chart's codomain** — the same language as
     ``generators`` / ``gram`` / ``det_gram``; ``fundamental_domain`` is
     a **section's image**, in the BASE's coordinates, which is what
     :meth:`DiscreteMeasure.quotient
     <orpheus.numerics.measure.DiscreteMeasure.quotient>` actually
     emits. They answer different questions and neither subsumes the
     other: `[M]` the chart is **Mode-12 blind** to the ERR-080 forgery
     while the section refuses it (:ref:`manifold-two-coordinate-systems`).
   - ⭐ **ERR-080's level-1 half is a botched SECTION of**
     :math:`S^2/SO(2)_x`. The realization is a chart; a consumer needed a
     section; the tree fabricated one by zero-padding to
     :math:`(\mu,0,0)`, which is off :math:`S^2` — `[M]` norms
     :math:`0.183\ldots0.960`, while an honest :math:`\varphi = 0`
     half-meridian is on the sphere to :math:`0.0`. ⚠ That names the
     **level-1** half only; the level-2 repair is still the trivial
     isotypic sub-basis :math:`\{Y_\ell^0\}\cong\{P_\ell\}`
     (:ref:`manifold-err-080-is-a-section`). Since tracker 2.4 the slab
     at least *names* the orbit space it needs a section of; declaring
     the section itself is tracker 2.3, and it is a CHOICE
     (:ref:`manifold-the-axis-convention-for-a-section`).
   - ⚠ **An orbit space is an ORBIFOLD, not a quotient manifold.**
     Where the action is not free, the image of the fixed-point set is
     a **singular stratum**. For :math:`S^2/SO(2)_a` that locus is
     :math:`\mu = \pm 1`, the poles — exactly where
     :math:`\det P` vanishes, and exactly where the curvilinear
     S\ :sub:`N` :math:`\alpha`-dome closes. A design that assumes a
     quotient is a smooth submersion is wrong there and only there
     (:ref:`manifold-singular-stratum`). For the :math:`\sigma_a` fold
     the stratum is the disk's **boundary circle** — a locus, not a
     point set, which is why ``singular_stratum`` is a symbolic locus
     and not a ``tuple[float, ...]``
     (:ref:`manifold-stratum-is-a-locus`).
   - **One polynomial, three appearances, three epistemic statuses.**
     :math:`(1-\mu^2)` is the squared :math:`SO(2)_a`-orbit radius
     (**derived**, this page); the redundant harmonic
     :math:`Y_2^{+2}` that makes the slab Gram rank-deficient
     (**measured**: `[M]` ratio spread :math:`8.9\times10^{-16}`); and
     the curvilinear angular-redistribution coefficient
     (**an identity of polynomials, with the mechanism unproved** —
     the reduction has not been derived, and the ruling that must
     settle it is Phase 1.3 of #429).
   - **The catalogue is the engine's SEED, not its rival.** A general
     orbit-space engine is **deferred, not refused**; the binding
     requirement is on the DATA MODEL — a catalogue entry must *be* the
     derivation procedure's output, so an engine ships by *computing*
     these fields instead of reading them, introducing no new type. The
     falsifiable check, and `[M]` **7 of 9** of the procedure's outputs
     are slots today (:ref:`manifold-engine-seed`).
   - ⚠ **This module imports nothing from** :mod:`orpheus.numerics` **at
     runtime, and that is load-bearing** — *more* so since tracker 2.4,
     which added the reverse edge. `[M]` ``symmetry.py`` now imports
     **both** ``manifold`` (for :class:`Quotient`, to read a polar
     marginal's axis) and ``measure`` at module scope, and ``measure``
     imports ``manifold``. So a module-scope ``manifold → symmetry`` edge
     would close a **direct 2-cycle**, not merely the 3-cycle
     ``measure → manifold → symmetry → measure`` it closed before
     (:ref:`manifold-import-cycle`).


.. _manifold-three-levels:

Three levels, and why a function space is not a domain
======================================================

A basis function is a map, and the question *"what is its domain?"* is
a question about what it **eats**:

.. math::

   Y_\ell^m : S^2 \longrightarrow \mathbb{R},
   \qquad
   P_\ell : [-1,1] \longrightarrow \mathbb{R},
   \qquad
   \mathbf{1}_{C_i} : \mathbb{R}^d \longrightarrow \{0,1\} .

The argument is a **point** — a unit direction, a cosine, a position —
not a vector of degrees of freedom. `[M]`
:class:`~orpheus.numerics.space.FunctionSpace` is documented as *"a
finite-dimensional vector space of discrete fields"* and carries a
``shape``: the tensor shape of the DOFs. So :math:`L^2(S^2)` is the
space the harmonics are **elements of**, and it is never the space they
are **maps from**. ``FunctionSpace`` answers *"what do these live in"*;
the domain question is *"what do these consume"*. Different arrows.

Separating them gives three levels. The tree named two:

.. list-table:: The three levels
   :header-rows: 1
   :widths: 14 30 28 28

   * - Level
     - The object
     - What a ``Basis`` carried
     - What a ``DiscreteMeasure`` carried
   * - 1 — the manifold
     - :math:`M`: :math:`S^2`, :math:`[-1,1]`, :math:`M/H`,
       :math:`\mathbb{R}^d`, energy, an index set
     - ⛔ **nothing**
     - ``support`` — a bare ``str``
   * - 2 — fields on :math:`M`
     - :math:`L^2(M)`, discretized at :math:`N` nodes
     - —
     - ``.space`` (a
       :class:`~orpheus.numerics.space.FunctionSpace`)
   * - 3 — coefficients
     - :math:`\mathbb{R}^K`
     - ``.space`` (a
       :class:`~orpheus.numerics.space.FunctionSpace`)
     - —

Read the table's last two columns across row 3: a basis and a measure
both have a ``.space``, and they are spaces of **different levels**.
That is not a naming accident to be tidied — a frame's whole job is the
pair of arrows between them
(:doc:`/theory/foundations/frame`). What is missing is row 1, on both
sides.

⭐ **And the level-2 check already passes, which is why the defect
survived.** On a slab the frame's analysis arrow
``measure.space → basis.space`` is between two perfectly well-formed
spaces: `[M]` on ``gauss_legendre(8).angular_frame(2)`` the domain is
``L2[S^2]`` of shape :math:`(8,)` and the codomain is
``spherical_harmonic_space`` of shape :math:`(3,5)`, the metric
resolves, and the pairing computes. Nothing is wrong at level 2.
**The check that fails is one level down**, on a manifold where no
object existed — so there was nothing to compare and nothing to refuse.

⚠ And note what the level-2 name is: `[M]` ``L2[S^2]``. It is
*derived* — ``measure.py:331`` builds it as ``f"L2[{self.support}]"``
— so the forged level-1 tag propagates upward verbatim, and a reader
inspecting the space sees a confident, wrong statement about the
manifold. A derived name is only as true as what it derives from.

.. note::

   **Why the operator-vocabulary collision is not one.** The proposal
   to spell a basis's level-1 slot ``support`` — to dodge
   :class:`~orpheus.numerics.operator.LinearOperator`'s existing
   ``domain`` / ``codomain`` vocabulary — was raised and **refuted** on
   2026-08-31, on two independent grounds. First, ``support`` is
   *mathematically false for a basis*: :math:`\operatorname{supp}(f)`
   means *where* :math:`f` *is non-zero*, and for an indicator basis
   that is exactly ONE cell per function — so
   ``IndicatorBasis.support = "spatial_R1"`` would state something
   untrue. (For the spherical-harmonic basis it is accidentally
   near-right, almost all of :math:`S^2`, which is what let the
   proposal past.) Second, the collision is a *word*, not a type:
   ``dom`` is the source of a morphism in both readings, in **Man** for
   a basis function and in **Vect** for an operator. Same functor,
   different categories. The slot is therefore ``domain``, and it
   :ref:`landed 2026-09-01 <manifold-seams>` (tracker 2.1) as an
   abstract property of the ABC — `[M]` answered by **6 of 6** shipped
   subclasses, the denominator being ``Basis.__subclasses__()`` walked
   recursively at runtime. ⛔ This sentence read *"it is not yet built"*
   until that day.


.. _manifold-string-tag:

What the string tag could not do
================================

Level 1 **was** ``Space = str`` (``measure.py:111``, retired 2026-09-01 by
tracker 2.0c), with the module's own comment on the alias set reading,
verbatim:

   *"Common aliases used across the project. These are recommendations,
   not constraints; user code may pass arbitrary strings."*
   — ``measure.py:113-114``, as it stood

That is an honest description of a slot with no semantics. Three
consequences followed from it, each measured; the type answers all three,
and the argument for minting it was these sites, not the size of the
migration. **All three are now closed in the code** — the sections below are
kept in the present tense of the tag they describe, because they are the
*reason* the type exists and re-tensing them would erase the evidence.

.. _manifold-err-080:

(a) Nothing could be refused — the ERR-080 forgery
--------------------------------------------------

A 1-dimensional angular quadrature carries no azimuthal information.
:meth:`Quadrature.angular_frame
<orpheus.numerics.quadrature.directional.Quadrature.angular_frame>`
nonetheless builds its measure by ``column_stack``\ ing three
axis-cosine arrays — two of which are a zero *fallback*, not data — and
declared the result ``support=SPACE_SPHERE`` (today ``support=SPHERE`` —
the forgery survived the retype; only its spelling changed). The rows are then
:math:`(\mu, 0, 0)` with
:math:`\lVert\Omega\rVert = |\mu| \neq 1`: points of :math:`[-1,1]`,
not of :math:`S^2`.

.. note::

   **Scoped, 2026-09-01.** The paragraph above is now a statement about the
   **1-D arm alone**, and that is the whole of the change: the
   ``column_stack`` used to run for *every* rule, so a Lebedev or
   level-symmetric frame also rebuilt a measure it had been handed. Today a
   rule whose nodes already are three-component directions hands the frame
   **its own measure**, and the construction survives only where there is
   genuinely nothing honest to build — `[M]` 10 of the 12 shipped rules route,
   2 do not. It lives in :meth:`Quadrature._harmonic_frame_measure
   <orpheus.numerics.quadrature.directional.Quadrature._harmonic_frame_measure>`
   with its retirement trigger written beside it.

   ⭐ The repair also reversed two losses this page did not record, because
   they are not what the forgery is *about*: the rebuilt measure carried three
   of :class:`~orpheus.numerics.measure.DiscreteMeasure`'s five fields, so it
   dropped ``invariance_group`` and ``exactness`` as well as falsifying
   ``support``. `[M]` 10 of 12 rules carry a group, **0 of 12 frames did** —
   and at that moment the frame's forged ``support`` was still a *string*
   tag, so :attr:`DiscreteMeasure.phase
   <orpheus.numerics.measure.DiscreteMeasure.phase>` matched none of its
   manifold arms and fell through to the ``invariance_group`` fallback,
   which the rebuild had just dropped. *The angular frame's own measure
   could not say it was angular*: it raised ``NotImplementedError`` on all
   twelve. ⟹ the transferable form, worth more than the instance: **"the
   rebuild loses X" is a completeness claim over the source type's FIELD
   LIST**, and its denominator is ``dataclasses.fields(T)`` — not the
   concept you happen to be chasing.

   ⚠ Re-measured at HEAD after the 2.0c retype: that same frame measure now
   carries a real :class:`~orpheus.numerics.manifold.Sphere` support and
   `[M]` answers ``phase == 'angular'`` from the manifold arm, with
   ``invariance_group`` and ``exactness`` still ``None``. The **forgery is
   unchanged** — the nodes still are not on :math:`S^2` — so read this as
   the raise having moved, not the defect (:ref:`manifold-err-080`).

`[M]` reproduced 2026-08-31 on
``Quadrature.gauss_legendre(8).angular_frame(2)``, reading the
production measure's own nodes:

.. list-table:: The declared support against the nodes
   :header-rows: 1
   :widths: 46 54

   * - Reading
     - Value
   * - ``measure.support``
     - ``'S^2'``
   * - :math:`\lVert\Omega_n\rVert`, sorted
     - ``0.1834 0.1834 0.5255 0.5255 0.7967 0.7967 0.9603 0.9603``
   * - rows within :math:`10^{-12}` of the sphere
     - **0 of 8**
   * - ``Sphere().contains(nodes)``
     - ``False``
   * - ``Sphere().contains(nodes / ‖nodes‖)`` — positive leg
     - ``True``
   * - ``Interval(-1, 1).contains(nodes[:, 0])``
     - ``True`` — the manifold the nodes *actually* belong to

The positive leg is not decoration. A contract predicate exercised only
against a broken instance validates the *raising*, not the *claim*
(``vv-principles`` #11); the third and fourth rows together are what
make the refusal a statement about :math:`S^2` rather than about the
function.

That forgery is the first link of **ERR-080** — a live wrong-answer
defect in :math:`P_L` scattering, publicly reachable through
``solve_sn(..., scattering_order >= 2)`` on any 1-D chart, with a
second symptom (a crash whose message blames a layer three hops
downstream) at higher :math:`(N, L)`. The full mechanism, its
:math:`\ell = 1` masking, the census over 105 slab rows and the
proposed repair are in the
:doc:`error catalogue </theory/verification/error_catalog>` (search
``ERR-080``) and in
:doc:`/theory/foundations/spherical_harmonics`; this page carries only
the part that is a **manifold** claim.

.. warning::

   ⛔ **The predicate exists; it is not wired.** `[M]` 2026-08-31,
   ``grep`` over ``orpheus/`` and ``tests/``: the only importers of
   :mod:`orpheus.numerics.manifold` are its own test module.
   :meth:`Quadrature.angular_frame
   <orpheus.numerics.quadrature.directional.Quadrature.angular_frame>`
   still writes a string, ``DiscreteMeasure`` still takes a ``str``,
   and **ERR-080 is still open**, held by an ``xfail(strict=True)``
   regression gate. Nothing on this page has repaired it. The
   construction-time refusal is a *capability that now exists*; wiring
   it is tracker items 2.0b / 2.1 of #429
   (:ref:`manifold-seams`).

.. _manifold-string-algebra:

(b) The algebra ran as string concatenation
-------------------------------------------

`[M]` re-measured by AST over ``orpheus/`` + ``tests/`` (2026-08-31,
independently of the plan that first reported it, and agreeing with
it): of **62** ``.support`` attribute reads — 31 in ``orpheus/``, 31 in
``tests/`` — **18** feed a string operation (an f-string, a
``.startswith``, or a ``+``), and **all 18 are in** ``orpheus/``, none
in ``tests/``. Four of the 18 are the ones that matter here:

.. list-table:: Four sites, three of them verbs, spelled as interpolation
   :header-rows: 1
   :widths: 22 40 38

   * - Site (as of 2026-08-31)
     - What it computed, until the verb replaced it
     - The verb it *is*
   * - ``measure.py:588``
     - ``new_space = f"{self.support} × {other.support}"``
     - :meth:`Manifold.__mul__
       <orpheus.numerics.manifold.Manifold.__mul__>` — the product
       manifold :math:`M \times N`
   * - ``measure.py:1022``
     - ``new_space = f"{self.support}/{group.name}"``
     - :meth:`Manifold.quotient
       <orpheus.numerics.manifold.Manifold.quotient>` — the orbit
       space :math:`M/H`
   * - ``measure.py:802``
     - ``f"φ_*({self.support})"``
     - the codomain of a pushforward :math:`\varphi_* \mu`
   * - ``measure.py:331``
     - ``name = f"L2[{self.support}]"``
     - the *derived* level-2 name — the level-1 object smuggled inside
       a level-2 label

⛔ **Every cell of the middle column is now history, and each of the four
died on a different date.** The column header read *"What it computes
today"* until 2026-09-02; it was true when written and the campaign
repealed it row by row. Rows 1 and 2 became :meth:`Manifold.__mul__
<orpheus.numerics.manifold.Manifold.__mul__>` and
:meth:`Manifold.quotient
<orpheus.numerics.manifold.Manifold.quotient>` at tracker 2.0c
(2026-09-01); row 4's ``str`` became a typed
:class:`~orpheus.numerics.manifold.Manifold` the same day (the
interpolation survives, but of a point set rather than a tag — see the
✅ below); and row 3 is the last to go, at tracker 2.3 (2026-09-02),
where the pushforward's codomain stopped being *named at the call site*
and became a field of the map (:ref:`manifold-arrow-type`). The rows
stay because the *argument* is what they carry: the mint added no
algebra, and this is the evidence.

This is the strongest available form of the project's own type-minting
criterion (``coding-standards``, *Type vs property*: mint **iff** there
are :math:`\ge 2` non-isomorphic realizations **and** a non-identity
morphism is applied). The morphisms are not merely *applicable* — they
are *shipped*. The mint adds no algebra; it gives a name to one that
already runs.

⚠ **And the last row was not a naming quibble — it shipped as a
falsehood.** ``measure.py:331`` at least *derives* its name, from a
``str``; ``basis/indicator_basis.py`` **hard-coded** it as
``f"L2[coarse_cells_R{self.ndim}]"``. And
:meth:`EnergyGrid.as_basis <orpheus.data.energy_grid.EnergyGrid.as_basis>`
builds an :class:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis`
over an **energy index** partition
(``edges = arange(n_groups + 1) - 0.5``), so that basis named its own
coefficient space after a *spatial* manifold it has nothing to do
with. `[M]` reproduced 2026-08-31 on a two-group grid, and again
2026-09-01 immediately before the repair:

.. code-block:: python

   >>> eg = EnergyGrid(edges=np.array([1e6, 1.0, 1e-5]))    # 2 GROUPS
   >>> mesh = Mesh1D(edges=np.array([0.0, 0.5, 1.0]), ...)  # 2 CELLS
   >>> eg.as_basis().space.name              # BEFORE #429 tracker 2.1
   'L2[coarse_cells_R1]'
   >>> mesh.indicator_basis().space.name
   'L2[coarse_cells_R1]'                     # ...the very same value
   >>> eg.as_basis().space == mesh.indicator_basis().space
   True

The two compared ``==`` **and hash-equal**, so a 2-group energy space and
a 2-cell spatial space were one value:
:class:`~orpheus.numerics.space.FunctionSpace` identity is
``(name, shape)``, and a false name is therefore not cosmetic but an
illegal state that IS representable.

✅ **REMEDIED 2026-09-01 by #429 tracker 2.1.** The
:class:`~orpheus.numerics.basis.base.Basis` ABC now asks every basis what
its functions EAT —
:attr:`~orpheus.numerics.basis.base.Basis.domain`, a
:class:`Manifold` — and an
:class:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis` takes the
manifold it partitions as a required constructor field, so the name
derives:

.. code-block:: python

   >>> eg.as_basis().space.name
   'L2[coarse_cells(energy)]'
   >>> mesh.indicator_basis().space.name
   'L2[coarse_cells(spatial_R1)]'
   >>> eg.as_basis().space == mesh.indicator_basis().space
   False

⭐ **What made the defect invisible for as long as it lived is worth
recording, because it is not carelessness.** At four of the five
production sites the basis and its
:class:`~orpheus.numerics.measure.DiscreteMeasure` are built in the SAME
function, three to five lines apart — and the *measure* named the
manifold correctly the whole time (``support="energy"``,
``"spatial_R1"``, ``f"index({label})"``). The answer was never
unavailable, only unasked; a hard-coded f-string is exactly the shape
that cannot be contradicted by the object sitting beside it. The durable
gate is therefore not *"the name is right"* — which any self-consistent
lie satisfies — but *"the two halves of one frame name ONE manifold"*
(``tests/numerics/test_basis_domain.py::test_d6``).

⭐ And assigning the type was itself a census: `[M]` it immediately
separated two manifolds the string tag ``"energy"`` had conflated — the
continuous energy axis in eV (:class:`Interval`, partitioned by
``tests/data/test_energy_grid.py``) from the multigroup *index* axis
(:class:`EnergyGroups`, what production partitions). Both have ambient
dimension 1, so no dimensional check could have found it; only naming the
point set does.

⭐ **And the slot answers a SECOND question, which is what tracker 2.1b
collected the same day.** Once a basis names the manifold its functions
eat, it has already declared the symmetry those functions *have*: a
function on an orbit space :math:`M/H` is, pulled back to :math:`M`,
exactly an :math:`H`-invariant function. So
:attr:`Basis.invariance_group
<orpheus.numerics.basis.base.Basis.invariance_group>` is **derived from**
``domain`` and stored nowhere — the second operand of the ERR-080 pairing,
obtained for no new field and no subclass edit
(:ref:`manifold-basis-invariance-group`).

The remaining half is the *measure's* side, and since tracker 2.0c it is a
**name** rather than a type: ``support`` became a :class:`Manifold` that
same day, so ``measure.py:371`` now derives ``f"L2[{self.support.name}]"``
from a typed point set. What is still missing is that
:class:`~orpheus.numerics.space.FunctionSpace` records only the resulting
string; one that carried its own manifold would collapse both spellings
into one — the level-2 half of this repair, tracked at
:ref:`manifold-seams`. ⛔ This paragraph read *"``support`` is still a
``str``, so ``measure.py:331`` derives a correct name from an untyped
tag"* until 2026-09-01, and it was true when written: 2.1 and 2.0c landed
hours apart, and the second one repealed the first one's premise.

.. _manifold-string-drift:

(c) The vocabulary drifted, and a nonsense quotient is spellable
----------------------------------------------------------------

`[M]` **2026-08-31, on the tree as it then stood** — read this
subsection in the past tense throughout — both ``'S^2/<sigma_y>'`` (a
hand-written ``new_space=`` in ``tests/numerics/test_measure.py``) and
``'S^2/sigma_y'`` (two further sites, and what the production
:meth:`DiscreteMeasure.quotient
<orpheus.numerics.measure.DiscreteMeasure.quotient>` emitted, since
``SubgroupOfO3.Mirror("y").name == "sigma_y"``) shipped. They named
**one** quotient and were **unequal under** ``==``. Also shipped as
support literals: ``'img'``, ``'probe'``, ``'[-1,1]^slab'``.

.. note::

   ⛔ **The mechanism this subsection describes is retired, and its two
   halves died separately.** The *tag* went at tracker 2.0c
   (2026-09-01), when ``support`` became a
   :class:`~orpheus.numerics.manifold.Manifold` and two strings naming
   one quotient stopped being expressible. The *hand-written target*
   went at tracker 2.3 (2026-09-02): ``new_space=`` is retired and a
   pushforward reads its codomain off the map
   (:ref:`manifold-arrow-type`). The subsection stays because the
   demonstration below — a nonsense quotient accepted by one verb and
   refused by the other — is the argument that the mint was a repair
   and not a re-spelling. ⚠ Line numbers in the original `[M]` are
   deliberately dropped here: which line of a test file carried a
   spelling is not a durable fact, and three of the four cited files
   have since been edited by this campaign.

The sharper demonstration is what the quotient verb accepts. `[M]`
2026-08-31, folding a :math:`\sigma_y`-invariant 4-node set by
``Mirror("y")`` under three different tags:

.. list-table:: What the two quotient verbs accept
   :header-rows: 1
   :widths: 34 33 33

   * - Declared ``support``
     - ``DiscreteMeasure.quotient`` result
     - :meth:`Manifold.quotient
       <orpheus.numerics.manifold.Manifold.quotient>`
   * - ``'S^2'``
     - ``'S^2/sigma_y'`` — accepted
     - ``NotImplementedError``, naming the missing derivation
   * - ``'probe'``
     - ``'probe/sigma_y'`` — accepted
     - (not a manifold; unspellable)
   * - ``'not_a_manifold_at_all'``
     - ``'not_a_manifold_at_all/sigma_y'`` — accepted
     - (not a manifold; unspellable)

⚠ **State the scope of that exactly, because the obvious summary is
false.** :meth:`DiscreteMeasure.quotient
<orpheus.numerics.measure.DiscreteMeasure.quotient>` is **not**
unchecked: it calls
:func:`~orpheus.numerics.symmetry.orbit_certificate` and refuses a
measure that is not :math:`G`-invariant. What is unchecked is the
**tag** — the new support is minted by interpolation from whatever
string the old one held, for any group, with no claim that :math:`M/H`
is a manifold anyone has derived. Two different objects are being
gated and un-gated: the *nodes* are certified, the *manifold* is
asserted.

⭐ And the same run shows the constructor's other blind spot: the four
nodes used above have :math:`\lVert x\rVert = \sqrt2`, and
``DiscreteMeasure(nodes=..., support="S^2")`` accepted them without
complaint. That is ERR-080's mechanism in four lines, with no
quadrature involved.


.. _manifold-closed-sum:

The type: a closed sum, split by TOTALITY
=========================================

The ruled shape (user, 2026-08-31) is a **closed sum**, and the axis of
the split is **totality**:

- the operations **every** manifold answers — ``dim``, ``name``,
  ``contains``, ``__mul__`` — are abstract on the base
  :class:`~orpheus.numerics.manifold.Manifold`, so no variant may
  abstain;
- the operations only a **quotient** can answer — the invariant
  generators, the syzygy ideal, the Procesi–Schwarz matrix and its
  determinant, the singular stratum, the provenance — live on
  :class:`~orpheus.numerics.manifold.Quotient` alone.

A sphere has no syzygy ideal, and under this shape it cannot be asked
for one: the question is a :exc:`TypeError`, not a ``None``.

.. _manifold-shape-evidence:

The evidence, not the taste
---------------------------

Two shipped precedents were **read**, not recalled, and the deciding
measurement is a field census.

.. list-table:: The two precedents, measured 2026-08-31
   :header-rows: 1
   :widths: 26 74

   * - Precedent
     - What it actually is
   * - :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`
     - A **registered sibling hierarchy**: `[M]` 7 direct subclasses
       (``AlbedoBoundary``, ``PeriodicBoundary``, ``PrescribedInflow``,
       ``ReflectiveBoundary``, ``VacuumInflow``, ``WhiteBoundary``,
       ``ZeroFluxBoundary``) and 7 registry keys. The right shape when
       every member answers the same questions differently.
   * - :class:`~orpheus.numerics.symmetry.SubgroupOfO3`
     - **ONE class over a** ``_tag`` **sum**: `[M]` 13 methods
       dispatching with **10** ``isinstance`` calls and **0** ``match``
       statements; not a dataclass, and **not frozen** —
       ``g._tag = 'MUTATED'`` succeeds. The right *data model* for a
       small stable member set, with the dispatch it would get today.

The ruled shape is the second precedent's data model with the first's
class-per-variant realization — the two halves of :math:`M/H` kept
structurally parallel — and the reason it is not a polymorphic
hierarchy over a nullable base is measurable:

.. list-table:: Where the derivation fields belong
   :header-rows: 1
   :widths: 40 30 30

   * - Field
     - Answerable by a quotient?
     - Answerable by a sphere / interval?
   * - invariant generators :math:`p_1 \ldots p_k`
     - yes
     - **no**
   * - syzygy ideal :math:`I`
     - yes
     - **no**
   * - :math:`P_{ij} = \langle\nabla p_i, \nabla p_j\rangle`
     - yes
     - **no**
   * - :math:`\det P`
     - yes
     - **no**
   * - the singular stratum
     - yes
     - **no**
   * - ``derived_by`` provenance
     - yes
     - **no**

`[M]` **every** derivation field is in the left column. On a
polymorphic hierarchy they would therefore sit on the *base*, returning
``None`` for every non-quotient — which is exactly the tax
:attr:`SubgroupOfO3.mirror_axis
<orpheus.numerics.symmetry.SubgroupOfO3.mirror_axis>` already pays
(`[M]` ``None`` for ``SO2('x')``, ``Dinfh``, ``Oh`` and ``O3``, ``1``
for ``Mirror('y')``) and which ``directional.py:522`` already branches
on to raise :exc:`NotImplementedError`. Repeating that tax on a
brand-new type, with six fields instead of one, is the design the
closed sum refuses.

⭐ **And tracker 2.4 doubled the tax rather than repaying it**, which
is worth knowing before reading it as an argument against
axis-parameterisation. The axial rotation group's axis needed the same
accessor, so :attr:`SubgroupOfO3.rotation_axis
<orpheus.numerics.symmetry.SubgroupOfO3.rotation_axis>` joined it as
the **continuous dual**: `[M]` ``SO2('x') → 0``, ``SO2('z') → 2``, and
``None`` for everything else — *including* the groups that CONTAIN
axial rotations (:math:`D_{\infty h}`, :math:`SO(3)`, :math:`O(3)`),
because the accessor identifies the group whose elements are *exactly*
the rotations about one coordinate axis, not any group in which such
rotations sit. The two accessors are **mutually exclusive** on the
shipped family: `[M]` over **18** groups (the six named entries, three
mirrors, three axial rotations, :math:`C_{2,3,4,6}` and :math:`D_{2h},
D_{6h}`), ``mirror_axis`` is non-``None`` on exactly the three
:math:`\sigma_a` and ``rotation_axis`` on exactly the three
:math:`SO(2)_a` — **zero** groups answer non-``None`` to both.

.. note::

   **Why not a phantom type parameter.** ``.claude/plans/sn_reshape.md``
   Issue 2 — quoted verbatim in ``measure.py:106-109`` — rules *"don't
   try to enforce* ``Space`` *types via Python generics; not expressive
   enough without runtime overhead"*. That ruling stands and **does not
   cover this**: it rejects ``Generic[Tag]`` *phantom* parameters,
   which ``coding-standards`` also rejects, because they are erased at
   runtime and do not specialize dunders — one ``__mul__`` body would
   serve every instantiation. A first-class **value** with real methods
   is a different proposal, and the one the three shipped morphisms
   above require.

.. _manifold-members:

The members
-----------

`[M]` **ten** variants — ``Manifold.__subclasses__()`` restricted to
this module, which is also what the exhaustiveness gate compares
against. They fall in three families (curved, flat, discrete) plus
**four** constructors that take another manifold as an argument:
``Product``, ``Ball``, ``FundamentalDomain`` and ``Quotient``. Two
dimensions travel with each member and are easy to conflate: the
**topological** dimension ``dim`` (what the manifold *is*) and the
**ambient** coordinate count (how many columns ``contains`` consumes).
They differ for every curved member.

⛔ **Correction (2026-08-31).** The first version of this paragraph read
*"nine variants … plus two recursive constructors"*. Both halves were
wrong when written: the mint shipped **eight** concrete variants
(``git show b8c05d16:orpheus/numerics/manifold.py | grep -c '^class .*(Manifold)'``
:math:`\to` 8) and the table below listed eight rows, so the prose and
its own table disagreed. The two-slot ruling then added
:class:`~orpheus.numerics.manifold.Ball` and
:class:`~orpheus.numerics.manifold.FundamentalDomain`, taking the count
to ten. The count is now stated with the command that produces it,
because a member roster is a universal and the shipped
``test_every_variant_is_reachable_from_this_modules_list`` is the only
thing that keeps it honest.

.. list-table:: The shipped members
   :header-rows: 1
   :widths: 16 8 8 20 48

   * - Variant
     - ``dim``
     - ambient
     - ``name``
     - Membership predicate
   * - ``Sphere``
     - 2
     - 3
     - ``S^2``
     - :math:`\bigl|\lVert x\rVert - 1\bigr| \le \varepsilon`
   * - ``Circle``
     - 1
     - 2
     - ``S^1``
     - :math:`\bigl|\lVert x\rVert - 1\bigr| \le \varepsilon`
   * - ``Interval(a, b)``
     - 1
     - 1
     - ``[a,b]``
     - finite, and :math:`a - \varepsilon \le x \le b + \varepsilon`
   * - ``RealSpace(d)``
     - :math:`d`
     - :math:`d`
     - ``spatial_Rd``
     - finite
   * - ``EnergyGroups(ng)``
     - 0
     - 1
     - ``energy``
     - integral, and :math:`0 \le g < n_g`
   * - ``IndexSet(label, n)``
     - 0
     - 1
     - ``index(label)``
     - integral, and :math:`0 \le i < n`
   * - ``Product(L, R)``
     - :math:`\dim L + \dim R`
     - sum
     - ``L × R``
     - each factor admits its own coordinate block
   * - ``Ball(d)``
     - :math:`d`
     - :math:`d`
     - ``D^d``
     - finite, and :math:`\lVert p\rVert^2 \le 1 + \varepsilon`
   * - ``FundamentalDomain(M, n⃗, ℓ)``
     - :math:`\dim M -` #antipodal pairs
     - :math:`M`'s
     - ``M|ℓ``
     - :math:`M` admits it, **and**
       :math:`\langle p, n_i\rangle \ge -\varepsilon` for every normal
   * - ``Quotient(M, H, …)``
     - :math:`\dim` of the realization
     - realization's
     - ``M/H``
     - **either** coordinate system — the chart's, or the
       fundamental domain's, dispatched on the ambient width
       (:ref:`manifold-two-coordinate-systems`)

`[M]` the names reproduce the retired ``SPACE_*`` string tags
**verbatim** — ``S^2``, ``S^1``, ``[-1,1]``, ``[0,1]``, ``[0,inf)``,
``R``, ``energy``, ``spatial_R1``, ``index(angular)`` — so the
migration off ``Space = str`` could not silently re-word an error message
or an ``L2[...]`` space name that a test pins. That is a deliberate
constraint on the type, gated by
``tests/numerics/test_manifold.py::TestTypeLaws::test_the_names_reproduce_the_retired_string_tags``.

✅ **The constraint paid, and it is now a fact rather than a design
intent.** `[M]` 2026-09-01, immediately before tracker 2.0c touched a call
site: **10 of 10** tag constants reproduce exactly, and so do the two
*derived* spellings the tag vocabulary built by hand — ``S^2/sigma_y`` via
:attr:`Quotient.name` (against ``f"{support}/{group.name}"``) and
``[-1,1] × [-1,1]`` via :attr:`Product.name` (against ``f"{a} × {b}"``). The
single divergence in the whole migration is the affine remap
``LEGENDRE.on(0, 1)``, whose f-string over floats read ``"[0.0,1.0]"`` where
``Interval(0.0, 1.0).name`` normalises to ``"[0,1]"`` — a re-baseline of two
test pins, and an improvement: the object is the same interval whatever the
float repr of its endpoints, which the string could not say.

Two consolidations are visible in that table and were the reason the
member list in the plan's own row was incomplete:

- ``Interval`` **is ONE family, not three tags.** The retired
  ``SPACE_INTERVAL_M11``, ``SPACE_INTERVAL_01``, ``SPACE_HALF_LINE``
  and ``SPACE_R`` are four *members* of it — which is what
  ``generating_measure.py:420``'s ``support=f"[{a},{b}]"`` was already
  saying, one interpolation at a time.
- ``IndexSet`` **was minted twice, under incompatible spellings.**
  `[M]` ``frame.py:759`` builds ``f"index({axis_label})"`` and
  ``sn/operators/loss_kernel_gauge.py:1169`` builds
  ``f"sn_trace_orbit{orbit}_g{group}"``, whose "points" are trace DOF
  indices cast to float. One family; ``label`` is what distinguishes
  the instances.

⚠ ``EnergyGroups`` is deliberately **not** a bare ``IndexSet``. A group
index and any other integer-noded counting rule are indistinguishable
from their nodes alone, which is precisely why the tag had to supply
the physical identity in the first place; the measure's derived
``phase`` depends on it.

⭐ **The last two members were minted by a DERIVATION, not by a survey**
— which is the healthiest reason for a variant to exist, and worth
recording as such. Neither ``Ball`` nor ``FundamentalDomain`` was on
anyone's list of manifolds the tree might want. They arrived because
the second catalogue entry produced objects the shipped member set
could not name:

- ``Ball`` — :math:`S^2/\langle\sigma_a\rangle` **is** the closed
  2-disk in invariant coordinates (:ref:`manifold-s2-sigma-y` derives
  it), and the nearest shipped 2-D member,
  ``Product(COSINE_INTERVAL, COSINE_INTERVAL)``, is the bounding
  **square**. The discriminator is measured, not stylistic: `[M]`
  :math:`(0.9, 0.9)` is in the square and **not** in the disk, and it
  is the image of no direction at all, since
  :math:`\eta^2 + \mu^2 = 1.62 > 1` forces :math:`\xi^2 = -0.62 < 0`.
  Adopting the square because it already ships would have admitted
  :math:`(\eta,\mu)` pairs that no :math:`\Omega` maps to
  (:ref:`manifold-realization-refuted`).
- ``FundamentalDomain`` — the section's image, the *other* of a
  quotient's two coordinate systems, and the one every measure the
  tree emits through ``.quotient(...)`` actually speaks
  (:ref:`manifold-two-coordinate-systems`). One rule covers both a
  half-space and a hyperplane: an **antipodal pair** of normals
  :math:`\{n, -n\}` spells the equality :math:`\langle p, n\rangle = 0`,
  so `[M]` ``FundamentalDomain(SPHERE, (e_y,), …).dim == 2`` (the
  :math:`\sigma_y` hemisphere) while
  ``FundamentalDomain(SPHERE, (e_y, -e_y, e_x), …).dim == 1``
  (an :math:`SO(2)` half-meridian) — from the same field, with no
  second slot and no flag.

.. note::

   **The membership tolerance is a construction tolerance, not a
   physics one.** ``_MEMBERSHIP_ATOL`` is :math:`10^{-12}`, chosen to
   match the construction quality of the shipped quadrature rules
   (whose nodes are exact to a few ULP). `[M]` a node at
   :math:`1 + 10^{-9}` is refused; at :math:`1 + 10^{-13}` and
   :math:`1 + 10^{-14}` it is admitted. No caller should widen it to
   make a measure fit — the whole value of the predicate is that a
   forged support has nowhere to hide, and the ERR-080 forgery misses
   by :math:`4\times10^{-2}` to :math:`8\times10^{-1}`, not by ULPs.
   ``contains`` is a **universal**, not a mean: one bad row is enough.


.. _manifold-orbit-space:

The orbit space, derived
========================

An orbit space :math:`M/H` is not declared and not guessed. It is
**computed**, by a standard construction from real invariant theory,
once per (manifold, group) pair. This section states the procedure,
runs it in full on the **first** pair the tree catalogued, and draws
the four consequences that pay for it downstream. The **second** pair —
the shipped cylindrical fold :math:`S^2/\langle\sigma_y\rangle` — is
run in the section after it, because its answer forced a change to the
data model and so belongs after the consequences rather than beside
them (:ref:`manifold-second-entry`).

`[M]` the catalogue holds **six** keys today —
``Sphere/SO2_x``, ``Sphere/SO2_y``, ``Sphere/SO2_z``,
``Sphere/sigma_x``, ``Sphere/sigma_y``, ``Sphere/sigma_z`` — served by
**two** procedures, because each family shares one derivation that reads
the axis off the group. The identity quotient :math:`M/\{e\} = M` is a
seventh answer and is *not* a table row: it is derived on the spot, for
every manifold, because it is a theorem (:ref:`manifold-twin-lookup`).

⛔ This sentence read *"**four** keys … ``Sphere/SO2``"* until
2026-09-01, when tracker 2.4 gave the axial rotation group its axis and
the single ``SO2`` key became three. The **procedure** count did not
move, and that is the point: parameterising a family costs keys, never
derivations (:ref:`manifold-so2-axis-is-a-parameter`).

.. _manifold-derivation-procedure:

The procedure
-------------

Let :math:`G` be a compact group acting orthogonally on
:math:`\mathbb{R}^n`, and let :math:`X \subseteq \mathbb{R}^n` be a
:math:`G`-stable real algebraic set. Five steps:

**1. Minimal generators of the invariant ring.** Find
:math:`p_1, \ldots, p_k` generating :math:`\mathbb{R}[x]^G`. That a
finite generating set exists is Hilbert–Weyl; finding it is the step a
Gröbner engine would automate.

**2. The orbit map separates orbits.** The polynomial map
:math:`p = (p_1, \ldots, p_k) : \mathbb{R}^n \to \mathbb{R}^k` is
proper and separates :math:`G`-orbits, so

.. math::

   \mathbb{R}^n / G \;\cong\; p(\mathbb{R}^n) \;\subseteq\; \mathbb{R}^k ,

and the image is a **closed semialgebraic** set. The orbit space is
therefore describable by finitely many polynomial equalities and
inequalities — which is what makes the whole construction finite.

**3. The syzygy ideal gives the equalities.** The generators need not
be algebraically independent; the relations among them are

.. math::

   I \;=\; \ker\bigl(\mathbb{R}[y] \to \mathbb{R}[x],\;
                     y_i \mapsto p_i \bigr),

computed by elimination. :math:`V(I)` — the variety of that ideal — is
the algebraic set the image lies in. When the invariants *are*
algebraically independent, :math:`I = (0)` and there are no equalities.

**4. Procesi–Schwarz gives the inequalities.** The image is cut out of
:math:`V(I)` by one positive-semidefiniteness condition on the
**gradient Gram matrix** of the invariants, re-expressed in the
invariants themselves:

.. math::
   :label: manifold-procesi-schwarz

   p(\mathbb{R}^n) \;=\;
   \bigl\{\, y \in V(I) \;:\; P(y) \succeq 0 \,\bigr\},
   \qquad
   P_{ij} \;=\; \bigl\langle \nabla p_i,\, \nabla p_j \bigr\rangle .

This is the theorem the whole construction turns on: Procesi and
Schwarz, *Inequalities defining orbit spaces*, **Inventiones
Mathematicae 81** (1985), 539–554. :math:`V(I)` supplies the
equalities; :math:`P \succeq 0` supplies the inequalities, and together
they are a complete description.

**5. Intersect with the ideal of** :math:`X`. Steps 1–4 quotient the
whole ambient space. Adjoining the defining ideal of :math:`X` — for
:math:`S^2`, the single relation :math:`\lVert x\rVert^2 = 1` —
restricts the answer to :math:`X/G`.

.. (vv-status rationale) manifold-procesi-schwarz is a
   LITERATURE-TRANSCRIBED theorem statement (Procesi & Schwarz 1985):
   it has no implementing ORPHEUS function to verify against and makes
   no solver claim. What IS verifiable is its INSTANCE at
   :eq:`manifold-s2-mod-so2`, whose P-matrix, determinant and stratum
   are recomputed symbolically by the foundation gates in
   tests/numerics/test_manifold.py::TestQuotient
   (test_the_procesi_schwarz_matrix_matches_the_hand_derivation,
   test_det_P_vanishes_exactly_on_the_singular_stratum). Those tests
   carry @pytest.mark.foundation and deliberately NO verifies(...),
   per vv-principles' foundation-tier rule.
.. vv-status: manifold-procesi-schwarz documented

.. _manifold-s2-so2:

The worked entry: :math:`S^2 / SO(2)_a`
----------------------------------------

This is the **first** entry the tree catalogued, and it is the entry
every 1-dimensional angular discretisation is secretly using. Every
line below was **re-derived and re-run** in this session, independently
of the catalogue, for **all three axes**, and then compared against it.

Write :math:`a` for the **rotation axis** and :math:`b, c` for the other
two — the same convention the mirror entry uses for its mirrored axis
(:ref:`manifold-s2-sigma-y`). The shipped slab and sphere geometries
spend :math:`a = x`; a product rule's polar factor is about
:math:`a = z`. **Why the axis is a parameter at all, rather than a
convention fixed once, is** :ref:`manifold-so2-axis-is-a-parameter` —
read it before treating the three entries as redundant.

**Step 0 — the group.** :math:`SO(2)_a` is the group of proper
rotations fixing :math:`\hat e_a` and mixing the other two coordinates,
acting on :math:`\mathbb{R}^3` by (shown for :math:`a = z`, and
:math:`\det R_\theta = +1` for every :math:`a`)

.. math::

   R_\theta =
   \begin{pmatrix}
     \cos\theta & -\sin\theta & 0 \\
     \sin\theta & \phantom{-}\cos\theta & 0 \\
     0 & 0 & 1
   \end{pmatrix}.

It is compact, **connected**, and :math:`\dim = 1` — every structural
difference from the mirror entry follows from those three words
(:ref:`manifold-chart-section-asymmetry`).

**Step 1 — the invariants.** Two invariants generate:

.. math::

   p_1 = x_a, \qquad p_2 = x_b^2 + x_c^2 .

`[M]` verified symbolically for general :math:`\theta`, **on all three
axes** — both satisfy :math:`p(R_\theta x) - p(x) = 0` after
``simplify``, with the non-invariant control :math:`x_b` correctly
reported **not** invariant in each case. The control matters: a check
that passes on everything is not a check.

**Step 3 — the syzygy ideal is empty.** The Jacobian (shown for
:math:`a = z`)

.. math::

   \frac{\partial (p_1, p_2)}{\partial (x, y, z)}
   =
   \begin{pmatrix} 0 & 0 & 1 \\ 2x & 2y & 0 \end{pmatrix}

has `[M]` generic rank **2** on every axis, equal to the number of
invariants, so :math:`p_1` and :math:`p_2` are algebraically independent
and :math:`I = (0)`. There are no equalities; the orbit space is cut out
by inequalities alone.

**Step 4 — the Procesi–Schwarz matrix.** With
:math:`\nabla p_1 = \hat e_a` and
:math:`\nabla p_2 = 2(x_b \hat e_b + x_c \hat e_c)`,

.. math::
   :label: manifold-s2-mod-so2

   P \;=\;
   \begin{pmatrix} 1 & 0 \\ 0 & 4 p_2 \end{pmatrix},
   \qquad
   \det P \;=\; 4 p_2 ,
   \qquad\text{so}\qquad
   \mathbb{R}^3 / SO(2)_a \;=\; \{\, p_2 \ge 0 \,\}.

The two invariants have orthogonal gradients everywhere — the first
points along :math:`a`, the second lies in the :math:`bc`-plane — which
is why :math:`P` is diagonal and the condition collapses to a single
inequality. Note that :math:`P` itself carries **no** trace of which
axis was chosen: the axis lives entirely in the *generators*, which is
exactly why the three entries are one derivation and three keys.

**Step 5 — restrict to the sphere.** Adjoining
:math:`p_1^2 + p_2 = 1` and writing
:math:`\mu = p_1 = \hat\Omega \cdot \hat e_a`:

.. math::

   \det P \big|_{S^2} \;=\; 4\,(1 - \mu^2),
   \qquad
   S^2 / SO(2)_a \;=\; \{\, \mu \in \mathbb{R} : 1 - \mu^2 \ge 0 \,\}
   \;=\; [-1, 1].

`[M]` SymPy's ``solve_univariate_inequality`` on
:math:`4 - 4\mu^2 \ge 0` returns ``Interval(-1, 1)`` on every axis — the
orbit space is *computed*, not asserted — and the zero set is exactly
:math:`\{-1, +1\}`.

⚠ **The three orbit spaces are isometric and their realizations are
identical, and they are still three different quotients.** `[M]`
``SPHERE.quotient(SubgroupOfO3.SO2('x')).realization ==
SPHERE.quotient(SubgroupOfO3.SO2('z')).realization`` is ``True`` (both
are ``Interval(-1.0, 1.0)``) while the two ``Quotient`` values compare
``False``. That is the correct answer and it is the whole point:
:math:`\mu` is the cosine against a *different* direction in each, so a
rule on :math:`S^2/SO(2)_x` and a rule on :math:`S^2/SO(2)_z` can carry
byte-identical nodes and mean different functions of direction. A type
that identified them would re-admit exactly the confusion
:ref:`manifold-so2-axis-is-a-parameter` describes.

.. (vv-status note) manifold-s2-mod-so2 is the INSTANCE of
   :eq:`manifold-procesi-schwarz` for the axial-rotation catalogue
   family. Its content IS checked, and tightly: the P-matrix and its
   determinant are recomputed symbolically and compared with
   sp.simplify, and the singular stratum is SOLVED for rather than
   compared to a literal, by
   tests/numerics/test_manifold.py::TestQuotient::{test_the_procesi_schwarz_matrix_matches_the_hand_derivation,
   test_det_P_vanishes_exactly_on_the_singular_stratum}.
   .
   ⛔ UN-SENTINELED 2026-09-01 (tracker 2.4), and the reason is worth
   keeping because the sentinel carried its own exit condition and the
   condition FIRED. This block used to read `.. vv-status:
   manifold-s2-mod-so2 documented`, arguing that the gates above are
   @pytest.mark.foundation and that vv-principles' foundation tier
   carries no verifies(...) marker by rule -- so a verifies edge would
   assert a claim class the gates do not make. Two things falsified it.
   (a) [M] tests/numerics/test_slab_orbit_space.py:258,
   test_d1_three_axes_three_quotients_one_derivation, now carries
   @pytest.mark.verifies("manifold-s2-mod-so2") and asserts THIS
   equation's content per axis -- the invariants p1 = x_a and
   p2 = x_b^2 + x_c^2 against the shipped generators, det P = 4 p_2,
   the realization, and the pairwise distinctness of the three
   quotients. (b) [M] the foundation/verifies exclusion is not this
   project's practice: 65 tests tree-wide carry BOTH markers, the
   algebra-of-record SymPy-identity shape (tests/derivations/*_symbolic.py
   is full of them), so the combination produces a real edge here as it
   does there.
   .
   Keeping the sentinel would have made the generated matrix contradict
   itself -- the label listed with 1 verifying test AND in the
   Documented-only set, which is excluded from the orphan gate. The
   directive is therefore removed rather than re-argued, and this note
   stays so a future auditor does not re-add one. The equation is now
   verifies-covered; the two TestQuotient gates above remain its
   tightest checks and carry no marker, by the same 65-test convention.

.. list-table:: My re-derivation against the shipped catalogue entry, all three axes
   :header-rows: 1
   :widths: 44 56

   * - Check
     - Result, `[M]` 2026-09-01, for :math:`a \in \{x, y, z\}`
   * - :math:`P` (mine) :math:`-` ``entry.gram``
     - ``simplify`` :math:`\to` the zero :math:`2\times2` matrix, **3 of
       3 axes**
   * - :math:`\det P` (mine) :math:`-` ``entry.det_gram``
     - ``simplify`` :math:`\to 0`, **3 of 3**
   * - :math:`(p_1 - x_a,\; p_2 - (x_b^2+x_c^2))` (mine)
       :math:`-` ``entry.generators``
     - ``simplify`` :math:`\to (0, 0)`, **3 of 3** — and this is the row
       that *sees* the axis: the shipped generators read
       ``(p1 - x, p2 - y**2 - z**2)`` for :math:`a = x` and
       ``(p1 - z, p2 - x**2 - y**2)`` for :math:`a = z`
   * - ``entry.realization``
     - ``Interval(a=-1.0, b=1.0)`` — the **same value** on all three
   * - ``entry.dim`` / ``Sphere().dim``
     - ``1`` / ``2``
   * - ``entry.syzygy``
     - ``()`` — empty, as derived
   * - ``entry.singular_stratum`` / ``entry.is_free``
     - ``1 - u0**2`` / ``False``. The stratum is a **locus**, and
       solving it gives :math:`\{-1, +1\}` — the poles
       (:ref:`manifold-stratum-is-a-locus`)
   * - ``entry.fundamental_domain``
     - ``None`` — and that is an answer, not a gap: **no** section of
       :math:`S^2 \to S^2/SO(2)_a` is canonical, since every
       half-meridian is one and none is distinguished
       (:ref:`manifold-two-coordinate-systems`)
   * - ``entry.name`` / ``entry.derived_by``
     - ``'S^2/SO2_x'``, ``'S^2/SO2_y'``, ``'S^2/SO2_z'`` / ``'hand'``
   * - ``type(entry.gram)``
     - ``ImmutableDenseMatrix``, so the ``Quotient`` is **hashable** —
       required by the memo below, and by any ``set``/``dict`` keyed on
       a measure's support

⚠ One reproduction note, kept because it is the cheap trap in step 4.
Re-expressing :math:`P` in the invariants is a substitution, and
``sp.Matrix(...).subs(x**2 + y**2, p2)`` **silently fails** on
:math:`4x^2 + 4y^2`, which does not literally contain the node
:math:`x^2+y^2`. My first run therefore produced
:math:`\det P = 4x^2 + 4y^2`, an empty solution set for the stratum,
and a spurious disagreement with the catalogue. The failure was mine,
not the entry's; ``factor`` before substituting (or, as the shipped
builder does, substitute :math:`x_b \to \sqrt{p_2},\, x_c \to 0` —
legal because the expression is constant on an orbit) fixes it. A
disagreement with a reference is not a refutation until you have
diagnosed whose it is.

.. _manifold-so2-axis-is-a-parameter:

Why the rotation axis is a PARAMETER — the tree carries two poles
------------------------------------------------------------------

The derivation above is written in :math:`a`, and the catalogue holds
three keys where it held one. That is not generality for its own sake.
It is forced, and the forcing is a *measured* property of this codebase
rather than a mathematical nicety: **ORPHEUS has two polar axes in
simultaneous use, and one Gauss–Legendre rule serves both.**

Until 2026-09-01 the axial rotation group was a parameter-free enum
member ``SO2``, realized about :math:`z` — its exactness criterion asked
whether every node had :math:`\rho = \sqrt{x^2+y^2} = 0`. The three
facts that make that untenable are all in the tree today:

.. list-table:: The two poles, `[M]` 2026-09-01 on the live tree
   :header-rows: 1
   :widths: 30 12 58

   * - Site
     - Pole
     - What fixes it there
   * - The slab / sphere polar marginal
     - :math:`x`
     - The 1-D embedding is :math:`(\mu, 0, 0)`, so the residual mirror
       is :math:`\sigma_x` and the marginal's nodes lie on the
       :math:`x`-axis.
   * - ``_evaluate_real_sh``
       (:mod:`orpheus.numerics.basis.spherical_harmonic_basis`)
     - :math:`x`
     - ``cos_theta = mu_x``. The real spherical-harmonic pole **is**
       :math:`x`, so :math:`Y_\ell^0 = P_\ell(\mu_x)` and the Legendre
       polynomials are the :math:`m = 0` members *of that* basis.
   * - Every product rule's polar factor
       (:func:`~orpheus.numerics.quadrature.product_mu_phi`)
     - :math:`z`
     - `[M]` on ``Quadrature.product(4, 8)``, the :math:`z` column takes
       exactly **4** distinct values and they are the ``leggauss(4)``
       nodes, while :math:`x` and :math:`y` take **9** each. The
       azimuth winds about :math:`z`.
   * - :math:`C_n`, :math:`D_{nh}`, :math:`D_{\infty h}`
     - :math:`z`
     - The standard setting for the finite families, which the lattice
       already assumed and this page already recorded.

⟹ **the same function**,
:func:`~orpheus.numerics.quadrature.gauss_legendre_on_mu`, is the raw
material of *both* rows: it is the slab's rule on :math:`S^2/SO(2)_x`
and the polar factor of a product rule on :math:`S^2/SO(2)_z`. A group
tag on that rule that did not name its axis would be a claim about the
wrong pole in one of the two uses, whichever way it was fixed.

**The symptom the bare tag actually produced.** The retired criterion
was ``all(hypot(x, y) <= atol)`` — :math:`\rho` measured **about**
:math:`z`, unconditionally. `[M]` run verbatim against the slab
marginal's own embedded nodes, whose first row is
:math:`(-0.9602898564975362,\, 0,\, 0)`, it returns **False**: the
:math:`\rho` it computes is :math:`|\mu|`, not :math:`0`. So the one
shipped rule whose orbit space *is* :math:`S^2/SO(2)_x` reported that it
was **not** invariant under the group it had been quotiented by. The
campaign plan for #429 records this as its "Part IV obstacle 1". With
the axis named, the same criterion is exact *and* discriminating: `[M]`
on ``Quadrature.gauss_legendre(8).measure``, ``SO2('x')`` is **True**
while ``SO2('y')`` and ``SO2('z')`` are **False**.

**Two alternatives were available and both were refused.**

*Fix the axis by fiat* — pick one pole and standardise. That is exactly
what the retired bare tag did, and it re-ships the defect shape the
reflection family had already been cured of on 2026-08-02: an unnamed
group realized on one axis while a consumer needs another
(:doc:`/theory/foundations/discrete_measures` records the ``Z_2`` case,
where ``product(4, 3)`` is :math:`\sigma_z`-closed and **not**
:math:`\sigma_x`-closed while the tag answered ``True``). One curable
family per repair is not a pattern; two is.

*Move the slab to the* :math:`z` *pole* — make one pole true. There is
no single pole to standardise **on**: the move would have to touch the
sweep, the real spherical-harmonic basis, every ``Mirror('x')`` row in
the geometry table, and every frozen slab snapshot — to buy a
convention, not a capability. And it would still leave the product
rule's factor and the slab's rule as *different* objects that a shared
function returns, which is the thing the axis parameter states directly.

⭐ **The general rule this instance is the second witness for.** A group
realized on a coordinate axis carries that axis as **data**, not as a
setting: :class:`~orpheus.numerics.symmetry.Mirror` reached this ruling
for the plane on 2026-08-02 and
:class:`~orpheus.numerics.symmetry.SO2` for the rotation axis on
2026-09-01, both after shipping a wrong answer, both by leaving the
parameter-free enum for a frozen dataclass. The enum is now exactly the
groups that have **no** axis to name: `[M]` **six** members —
``Trivial``, :math:`D_{\infty h}`, :math:`O_h`, :math:`I_h`,
:math:`SO(3)`, :math:`O(3)` — beside **four** parameterised families
:math:`C_n`, :math:`D_{nh}`, :math:`\sigma_a`, :math:`SO(2)_a`.

.. note::

   ⚠ **Why** :math:`C_n` **did not follow.** It is cyclic about
   :math:`z` and is *not* axis-parameterised, and that asymmetry is
   deliberate rather than an oversight: no consumer has yet needed a
   :math:`C_n` about another axis. The project's own rule is to unify
   after the **second** instance, not before it
   (``coding-standards``), and both axis families crossed that line by
   shipping a wrong answer. When a :math:`C_n` about :math:`x` first
   appears, this is the paragraph that says what to do.

.. _manifold-so2-axis-lattice:

What the axis buys in the containment lattice
----------------------------------------------

The parameter is not decoration on a name — it changes the **order**.
:math:`SO(2)_a` is the one continuous family with a parameter, so its
relations are neither in the enum-to-enum table nor computable from a
finite realization, and they live in their own arm
(``symmetry._so2_contains``). `[M]` measured on the live tree, over all
three axes:

.. list-table:: The axial rotation group's edges
   :header-rows: 1
   :widths: 34 16 50

   * - Relation
     - Holds
     - Why
   * - :math:`SO(2)_a \subseteq SO(3),\, O(3)`
     - **every** :math:`a`
     - Proper rotations, inside the proper rotations.
   * - :math:`SO(2)_a \subseteq D_{\infty h}`
     - :math:`a = z` **only**
     - :math:`D_{\infty h}`'s continuous factor *is* :math:`SO(2)_z`; a
       rotation about :math:`x` does not preserve the :math:`z` axis.
   * - :math:`C_n \subseteq SO(2)_a`, every :math:`n`
     - :math:`a = z` **only**
     - :math:`SO(2)_z = \bigcup_n C_n` in the standard setting.
   * - :math:`\{e\} \subseteq SO(2)_a`
     - **every** :math:`a`
     - And nothing else finite: a mirror and every :math:`D_{nh}` carry
       :math:`\det = -1` elements.
   * - :math:`SO(2)_a \subseteq SO(2)_b`, :math:`a \ne b`
     - **never**
     - Two distinct axial rotation groups meet only in :math:`\{e\}`.

⭐ **The lattice and the invariance predicate were re-checked against
each other, which is the gate that catches a wrong edge.** For every
asserted edge :math:`A \subseteq B` and every measure :math:`m`,
:math:`B\text{-invariant}(m) \Rightarrow A\text{-invariant}(m)` must
hold — the compatibility law (``vv-principles`` #15), which is the loop
that exposed ERR-072. `[M]` re-run 2026-09-01 over **15** groups (the
six named entries, three mirrors, three axial rotations, :math:`C_2`,
:math:`C_4`, :math:`D_{2h}`) × **6** fixtures (the declared slab
marginal, the chart-level :math:`\mu`-rule, a marginal declared about
:math:`z`, ``product(4,8)``, ``level_symmetric(4)``,
``folded_product(4,8)``): **0 violations over 342 (edge × fixture)
pairs.**

The most informative row is the one that could not exist before: `[M]`
``gauss_legendre_on_polar_orbit(8, "z")`` — the same eight nodes,
declared about :math:`z` — is invariant under :math:`SO(2)_z`,
:math:`D_{\infty h}` **and** :math:`C_4`, while the slab's
:math:`x`-declared rule is invariant under :math:`SO(2)_x` and none of
those three. Same floats, different groups, because the axis says where
the nodes are embedded.

.. _manifold-quotient-is-memoised:

The derivation runs once — ``Manifold.quotient`` is memoised
-------------------------------------------------------------

An orbit space is **derived once and recorded**; that is the
catalogue's whole philosophy, and until 2026-09-01 nothing enforced it,
because nothing on a hot path asked for a quotient. Tracker 2.4 changed
that: *every slab quadrature now carries one*, so an unmemoised
:meth:`Manifold.quotient
<orpheus.numerics.manifold.Manifold.quotient>` would put a SymPy
derivation on the construction path of every slab solve.

`[M]` 2026-09-01, ``SPHERE.quotient(SubgroupOfO3.SO2('x'))`` on this
machine (host ``.venv``, CPython 3.14), cache cleared before the cold
reading:

.. list-table::
   :header-rows: 1
   :widths: 30 22 48

   * - Reading
     - Cost
     - Note
   * - cold (the SymPy derivation)
     - **6.43 ms**
     - Gradients, the ``simplify`` chain, the determinant.
   * - warm (mean over 1000 lookups)
     - **0.76 µs**
     - ``CacheInfo(hits=1000, misses=1)``.
   * - ratio
     - :math:`\approx 8500\times`
     -

Two properties make the memo *legal*, and both are worth stating
because a memo on the wrong object is a shared-mutable-state bug:
every :class:`~orpheus.numerics.manifold.Manifold` is a **frozen value
type** and hashable — which is why tracker 2.4 also retyped every
builder's ``gram`` to :class:`sympy.ImmutableMatrix`, since a mutable
``Matrix`` field makes the whole ``Quotient`` unhashable — and the
returned :class:`~orpheus.numerics.manifold.Quotient` is itself frozen,
so sharing one object across callers is safe for the same reason.

⭐ **A second memo landed in the same step, and it is the larger
number.** The containment machinery asks
:func:`~orpheus.numerics.symmetry.SubgroupOfO3.contains` for a *realized
operator set* before it consults any lattice table, so every question
involving :math:`I_h` rebuilt the icosahedral group's 120-element
closure. That was tolerable while the invariance walk offered one axial
rotation group; offering **three** multiplied the traffic. `[M]` on this
machine: ``_icosahedral_ops`` costs **155 ms** cold and returns 120
elements; a single ``maximal_invariance_groups`` walk calls it **33**
times on the slab marginal and **48** times on ``product(4, 8)`` — so
without the memo the slab walk alone would spend :math:`33 \times 155\
\text{ms} \approx 5.1\ \text{s}` rebuilding a constant. Memoised, `[M]`
the same walk is **5.3 ms** warm (197.5 ms on the first call of a fresh
process, which is the one cold build).

.. _manifold-dimension-drop:

Consequence 1 — the dimension drops by the group's
---------------------------------------------------

:math:`\dim S^2 = 2`, :math:`\dim SO(2) = 1`, and the quotient has
dimension 1. That is the generic count: the orbits are the
constant-latitude circles, each 1-dimensional, and the quotient records
one number per orbit. It is *only* generic — at the two poles the orbit
is a single point, and the drop there is 2, not 1. Which is the next
consequence.

.. _manifold-singular-stratum:

Consequence 2 — the action is not free, so the quotient is an ORBIFOLD
-----------------------------------------------------------------------

The poles :math:`\mu = \pm 1` lie on the rotation axis, so their
stabilizer is the whole of :math:`SO(2)`: the action there has a fixed
point, the orbit collapses from a circle to a point, and the
pushforward of the uniform measure vanishes. The image of the
fixed-point set is the **singular stratum**.

Concretely, :math:`[-1,1]` is a manifold *with boundary* — an orbifold
— not a quotient manifold, and

.. math::

   \det P \big|_{S^2} \;=\; 4\,(1-\mu^2) \;=\; 0
   \quad\Longleftrightarrow\quad \mu = \pm 1 :

**the stratum is exactly where** :math:`\det P` **vanishes**, which is
why the shipped entry *derives* it rather than declaring it (the
foundation gate solves :math:`\det P = 0` and compares against
``singular_stratum``). Anything designed on the assumption that a
quotient is a smooth submersion is wrong there — and *only* there,
which is what makes the stratum worth carrying as a field rather than
a caveat.

⭐ Two shipped objects already live on that stratum, from opposite
directions:

- the curvilinear S\ :sub:`N` :math:`\alpha`-dome **closes** at
  :math:`\mu = \pm 1`, because the redistribution coefficient
  :math:`(1-\mu^2)` vanishes there
  (:doc:`/theory/methods/sn/curvilinear_one_group`);
- the spherical-harmonic evaluator's ``on_axis`` guard fires when
  :math:`\sin\theta \approx 0`, i.e. on directions *along* the polar
  axis — the same locus, detected numerically and named nothing.

⭐ The second catalogued entry has an exact cylindrical analogue of the
first bullet, measured on production data: the fold's stratum is the
disk's boundary circle, the shipped quadrature nodes sit strictly
inside it, and the march seeds sit exactly **on** it, where the
:math:`\alpha`-dome closes (:ref:`manifold-orbifold-discretised`).

.. _manifold-one-polynomial:

Consequence 3 — one polynomial, three appearances
--------------------------------------------------

:math:`(1-\mu^2)` shows up three times in this corpus. They are the
same polynomial; whether they are the same *object* is three different
questions, with three different answers, and conflating them is the
error this subsection exists to prevent.

.. list-table:: :math:`(1-\mu^2)`, by epistemic status
   :header-rows: 1
   :widths: 26 44 30

   * - Appearance
     - What it is
     - Status
   * - the squared orbit radius
     - a point of :math:`S^2` at latitude :math:`\mu` sits at distance
       :math:`\sqrt{1-\mu^2}` from the rotation axis, so its
       :math:`SO(2)`-orbit is a circle of that radius and
       :math:`\det P = 4\,r_{\rm orbit}^2`
     - **DERIVED** on this page, from
       :eq:`manifold-procesi-schwarz`
   * - the redundant harmonic
     - on a 1-D rule the retained :math:`Y_2^{+2}` column is
       :math:`\propto (1-\mu^2)`, which is what makes the discrete
       Gram rank-deficient
     - **MEASURED** — see below
   * - the angular-redistribution coefficient
     - the :math:`(1-\mu^2)/r \cdot \partial_\mu \psi` term in the
       spherical streaming operator
     - ⚠ **an identity of polynomials only.** The *mechanism* —
       that the redistribution term is the connection of the
       phase-space quotient — is **unproved**

**The measured one.** `[M]` 2026-08-31, reading the live frame's own
basis and measure at :math:`L = 2`. The coefficient array has shape
:math:`(3,5) = 15` slots, of which **9** are actual
:math:`(\ell, m)` harmonics and 6 are padding. On a genuine 3-D rule
all 9 light up: ``lebedev(11)`` gives **9 live of 15**. On the slab
``gauss_legendre(8)`` only **5** do —
:math:`\{(0,0),(1,0),(2,0),(2,1),(2,2)\}` — and the discrete Gram
over those five has ``matrix_rank`` **4**, on live singular values
:math:`2.70755,\; 1.41922,\; 4.92450\times10^{-1},\;
4.74468\times10^{-2}`. So one of the five is a linear combination of
the others.

Which one is the point: dividing the :math:`(2,2)` column by
:math:`(1-\mu^2)` gives a **constant** :math:`0.866025`, with a spread
of :math:`8.9\times10^{-16}` across all eight nodes. The column *is*
that polynomial, to round-off — and :math:`(1-\mu^2)` is
:math:`\det P / 4`, which on the quotient is a function of
:math:`\mu` alone and therefore already spanned. The rank deficiency
is a theorem about :math:`S^2/SO(2)_x`, not a conditioning accident.

⚠ Do not quote a fifth singular value. It is a noise-floor reading and
`[M]` does not reproduce between runs; the reproducible statements are
the four live values, the counts **9 of 15** (3-D) versus **5 of 15**
(slab), and the rank **4**.

⚠ **And do not read the third row as settled.** The polynomial identity
is real and it is suggestive — the same expression appears as the
orbit-space boundary and as the coefficient that vanishes at the same
locus — but "the curvilinear redistribution term *is* the quotient's
connection" is a claim about a **reduction that has not been carried
out**. It is Phase 1.3 of #429, whose stated done-when admits exactly
two outcomes: a derivation, or an explicit ruling that the coincidence
is accidental in 1-D spherical geometry. Until one of those lands, the
honest statement is the one in the table: three occurrences of one
polynomial, one of them derived here, one of them measured, one of them
open.

⭐ The second catalogued entry produces the **cylindrical twin** of that
open row — :math:`1 - \eta^2 - \mu^2` is simultaneously the fold's
:math:`\det P` on :math:`S^2` and the locus on which the cylindrical
:math:`\alpha`-dome closes and the march seeds sit, measured exactly
(:ref:`manifold-orbifold-discretised`). It is the same *kind* of
statement as the third row here — an identity of loci with an unproved
mechanism — so closing one does not close the other, and neither may be
cited for the other.

.. _manifold-gelfand:

Consequence 4 — the quotient is a Gelfand pair, so :math:`\Lambda` is forced
-----------------------------------------------------------------------------

:doc:`/theory/foundations/spherical_harmonics` and
:doc:`/theory/foundations/frame` already own the statement that the
P\ :math:`_\ell` scattering kernel factors as :math:`R\,\Lambda\,M` with
:math:`\Lambda` diagonal, and the derivation (Funk–Hecke plus Schur's
lemma, read as the spectral theorem :math:`A = U\Sigma U^{*}`) lives at
:eq:`sh-funk-hecke-eigenvalue` and :ref:`frame-eigenbasis-ownership`.
**Edited there, consumed here.** What this page adds is the *quotient*
register — the reason the same factorization is forced by the orbit
space, stated in the vocabulary of :math:`S^2/SO(2)_a` rather than of the
scattering operator.

:math:`S^2` is itself a homogeneous space, :math:`S^2 = SO(3)/SO(2)`,
and :math:`(SO(3), SO(2))` is a **Gelfand pair** — the convolution
algebra of :math:`SO(2)`-bi-invariant functions on :math:`SO(3)` is
commutative. So the object under discussion is really the double coset
space

.. math::

   SO(2) \,\backslash\, SO(3) \,/\, SO(2),

whose bi-invariant functions are the **zonal spherical functions** of
the pair — for :math:`(SO(3), SO(2))` exactly the Legendre polynomials
:math:`P_\ell`. Commutativity of that algebra is what forces every
:math:`SO(3)`-equivariant zonal operator to be simultaneously
diagonalised with a multiplier depending on :math:`\ell` alone. Read
this way, :math:`R\,\Lambda\,M` with :math:`\Lambda` diagonal is a
**theorem of harmonic analysis on the quotient**, not a chosen
factorization that happens to work — which is the same conclusion the
two pages above reach from the operator side, by a different route.

⭐ The same lens names the repair ERR-080 needs. The sub-basis that
fixes a 1-D chart is not a "zonal special case" bolted on for slabs; it
is the **trivial isotypic component** of the :math:`SO(2)` action —
:math:`\{Y_\ell^0\} \cong \{P_\ell\}` — and that says *why those
slots and not others*. It is also the reason the repaired Gram is
expected to be exactly diagonal rather than merely better conditioned:
Gauss–Legendre integrates :math:`P_\ell P_{\ell'}` exactly to degree
:math:`2N-1`, so the Gram becomes exactly :math:`2/(2\ell+1)` on the
diagonal. The falsifiable form of that prediction is recorded with
ERR-080, not here.


.. _manifold-second-entry:

The second entry, and the two coordinate systems it forced
==========================================================

The first catalogued entry answered the question *what is*
:math:`M/H` with a single object, and
:class:`~orpheus.numerics.manifold.Quotient` stored it in a single
field, ``realization``, which both
:meth:`contains <orpheus.numerics.manifold.Quotient.contains>` and the
ambient-width helper read. The second entry — the **shipped cylindrical
fold** :math:`S^2/\langle\sigma_y\rangle`, which
:meth:`Quadrature.folded_product
<orpheus.numerics.quadrature.directional.Quadrature.folded_product>`
performs on every curvilinear rule — cannot be stored that way, and the
reason is not a detail of the entry. It is that an orbit space has
**two** honest coordinate systems, and the tree produces data in both.

⚠ **Why the first entry could not expose the fork — and it is not the
dimensions.** The chart codomain and a section have the *same* ``dim``
in **both** entries, so dimension cannot discriminate them: `[M]`
``Interval(-1,1).dim`` is :math:`1` and an :math:`SO(2)` half-meridian
written as ``FundamentalDomain(SPHERE, (e_y, -e_y, e_x), …).dim`` is
also :math:`1`; ``Ball(2).dim`` is :math:`2` and the
:math:`\sigma_y` hemisphere ``FundamentalDomain(SPHERE, (e_y,), …).dim``
is also :math:`2`. Indeed the agreement is now a **construction law**,
gated in
:meth:`Quotient.__post_init__ <orpheus.numerics.manifold.Quotient>` —
a quantity that must always agree cannot tell two cases apart. Two
measured facts hid the fork instead:

1. **No section of** :math:`S^2 \to S^2/SO(2)_a` **is canonical.** Every
   half-meridian is one and none is distinguished, so there was nothing
   to put in a second slot even had one existed. That is the normal
   situation for a positive-dimensional group.
2. **The tree's** :math:`SO(2)` **data is already in chart
   coordinates.** `[M]` ``gauss_legendre(8).measure.nodes`` has shape
   :math:`(8,)` and holds the invariant :math:`\mu` itself, so the
   realization and the data speak the same language and nobody had to
   choose. `[M]` ``folded_product(4,8).measure.nodes`` has shape
   :math:`(16,3)` — the base's ambient columns — so for the second
   entry the same slot cannot even *see* the tree's own nodes.

.. _manifold-s2-sigma-y:

The derivation — :math:`S^2/\langle\sigma_a\rangle` is the closed disk
-----------------------------------------------------------------------

The procedure of :ref:`manifold-derivation-procedure`, run in full.
Write :math:`a` for the mirrored axis and :math:`b, c` for the other
two; the shipped fold is :math:`a = y`. Every line below was re-derived
in this session, independently of the catalogue entry, and then
compared against it.

**Step 0 — the group.**
:math:`\sigma_a : x_a \mapsto -x_a`, with
:math:`H = \langle\sigma_a\rangle = \{e, \sigma_a\}` of order 2. `[M]`
:math:`\det\sigma_y = -1` and :math:`\sigma_y^2 = I`. That
determinant is not decoration: :math:`\sigma_a` is an **improper**
element, and specifically a *reflection* — it fixes a hyperplane
pointwise — which predicts steps 1 and 3 before either is run.

**Step 1 — the invariants.** A polynomial is :math:`\sigma_a`-invariant
iff it is **even in** :math:`x_a`, so

.. math::

   p_1 = x_b, \qquad p_2 = x_c, \qquad p_3 = x_a^2 .

`[M]` verified symbolically, with two non-invariant controls: for
:math:`x` and :math:`z` and :math:`y^2` the difference
:math:`p(\sigma_y x) - p(x)` is :math:`0`, while for the controls
:math:`y` and :math:`xy` it is :math:`-2y` and :math:`-2xy`. A check
that passes on everything is not a check.

*Completeness*, by Molien's series — not by eyeballing. Molien's
formula gives the Hilbert series of :math:`\mathbb{R}[x]^H` from the
group alone, and `[M]`

.. math::

   M(t) \;=\; \tfrac12 \sum_{g \in H} \frac{1}{\det(I - t g)}
         \;=\; \frac{1}{(1-t)^2 (1-t^2)}
         \;=\; 1 + 2t + 4t^2 + 6t^3 + 9t^4 + 12t^5 + O(t^6),

which is *exactly* the Hilbert series of a free polynomial algebra on
generators of degrees :math:`1, 1, 2`. `[M]` the difference of the two
series simplifies to :math:`0`. That single equality carries two facts
at once: **completeness** — the subalgebra
:math:`\mathbb{R}[x_b, x_c, x_a^2]` has the same graded dimension as
the whole invariant ring in every degree, so it *is* the whole ring —
and **freeness**, i.e. an empty syzygy ideal, which step 3 then
re-derives independently.

*Minimality*, by counting :math:`\dim(\mathfrak{m}/\mathfrak{m}^2)`
degree by degree, where :math:`\mathfrak{m}` is the ideal of
positive-degree invariants: `[M]` degree 1 contributes **2** new
generators, degree 2 contributes **1**, and degrees 3–5 contribute
**0** — so the minimal generating set has exactly three members.

.. warning::

   **A trap in that count, recorded because it produced a
   self-consistent wrong answer.** A "decomposable" is a product of
   **two or more** positive-degree invariants. Counting products of
   :math:`k \ge 1` factors instead includes the generators themselves,
   and then every degree reports *"0 new generators needed"* — i.e. the
   empty set generates the invariant ring. The output is internally
   consistent and completely wrong. The predicate is :math:`k \ge 2`.

**Step 3 — the syzygy ideal is empty, and predictably so.** `[M]` the
elimination
:math:`\langle u_1 - x_b,\, u_2 - x_c,\, u_3 - x_a^2\rangle \cap
\mathbb{R}[u]` returns a Gröbner basis with **no** :math:`u`-only
generator, so :math:`I = (0)`; the Jacobian
:math:`\partial(p_1,p_2,p_3)/\partial(x,y,z)` has determinant
:math:`\pm 2 x_a` — `[M]` exactly :math:`-2y` for the shipped
:math:`a = y` ordering, the sign being an artefact of which two
coordinates are kept first — and generic rank **3**, equal to the
number of invariants. It is the **rank** that carries the argument.

⭐ **But the theorem is better than the computation here.** A
reflection generates a *reflection group*, so by
**Chevalley–Shephard–Todd** its invariant ring is a polynomial ring —
hence free, hence :math:`I = (0)`. The elimination is the mechanical
route; the answer was never in doubt. This is a real structural
contrast with :math:`S^2/SO(2)_a`, whose syzygy ideal is *also* empty but
for the unrelated reason that its two invariants happen to be
algebraically independent: :math:`SO(2)` is not a reflection group and
no such theorem applies to it.

**Step 4 — Procesi–Schwarz.** The three gradients
:math:`\nabla p_1 = e_b`, :math:`\nabla p_2 = e_c`,
:math:`\nabla p_3 = 2 x_a e_a` are mutually orthogonal, so :math:`P` is
diagonal and :math:`P \succeq 0` collapses to a single inequality:

.. math::
   :label: manifold-s2-mod-mirror

   P \;=\; \operatorname{diag}\bigl(1,\, 1,\, 4 p_3\bigr),
   \qquad
   \det P \;=\; 4 p_3 ,
   \qquad\text{so}\qquad
   \mathbb{R}^3 / \langle\sigma_a\rangle \;=\; \{\, p_3 \ge 0 \,\}.

.. (vv-status rationale) manifold-s2-mod-mirror is the second INSTANCE
   of :eq:`manifold-procesi-schwarz`, for the S^2/<sigma_a> catalogue
   entry, and is classified exactly as its sibling
   :eq:`manifold-s2-mod-so2` is, for the same structural reason. Its
   content IS checked, and tightly: the P-matrix and its determinant are
   recomputed symbolically and compared with sp.simplify, the syzygy is
   asserted empty, and the stratum is SOLVED for (and shown to be a
   one-parameter family, i.e. a curve) rather than compared to a
   literal, by tests/numerics/test_manifold.py::
   TestTheSigmaYFoldIsExpressibleAndDiscriminating::{
   test_the_derivation_reproduces_procesi_schwarz,
   test_the_stratum_is_a_CIRCLE_not_a_point_set}. Those gates carry
   @pytest.mark.foundation and deliberately NO verifies(...): they
   assert an intrinsic law of a data type, with no flux, eigenvalue or
   convergence claim behind them, and vv-principles' foundation tier
   carries no verifies marker by rule. A verifies edge here would mint a
   coverage claim of a class the gates do not make.
.. vv-status: manifold-s2-mod-mirror documented

**Step 5 — restrict to the sphere.** Adjoining the sphere's ideal
:math:`p_1^2 + p_2^2 + p_3 = 1` makes :math:`p_3` *eliminable*, and
what remains is two-dimensional:

.. math::

   \det P \big|_{S^2} \;=\; 4\,\bigl(1 - p_1^2 - p_2^2\bigr),
   \qquad
   S^2/\langle\sigma_a\rangle \;\cong\;
   \bigl\{\, (p_1, p_2) : p_1^2 + p_2^2 \le 1 \,\bigr\} \;=\; D^2 .

`[M]` **the re-derivation agrees with the shipped entry exactly.** The
construction, so it regenerates from this page: form the three
gradients of :math:`(x, z, y^2)` symbolically, build
:math:`P_{ij} = \langle\nabla p_i, \nabla p_j\rangle`, substitute
:math:`y^2 \to p_3` to re-express it in the invariants, then compare
against ``SPHERE.quotient(SubgroupOfO3.Mirror("y"))`` with the entry's
own generator symbol mapped onto :math:`p_3`. ``sympy.simplify`` of the
difference is the zero :math:`3\times3` matrix, and of the
determinants, :math:`0`.

⚠ The substitution step has a trap the sibling entry records
(:ref:`manifold-s2-so2`), and it does **not** bite here: ``subs``
matches syntactic nodes, so ``subs(x**2 + y**2, p2)`` silently fails on
:math:`4x^2+4y^2` — whereas ``4*y**2`` literally contains the node
``y**2``, so ``subs(y**2, p3)`` succeeds. A derivation that "happens to
work" is still worth asserting: check that no free :math:`x, y, z`
remains in :math:`P` after the substitution.

⭐ **Step 5 supplies the equality that the syzygy ideal did not**, and
that is worth stating before an engine is written. In *both* catalogued
entries ``syzygy`` is honestly ``()``, and in both the real equality
arrives at step 5 — :math:`p_1^2 + p_2 = 1` for :math:`SO(2)`,
:math:`p_1^2 + p_2^2 + p_3 = 1` here. An engine that emits only
:math:`I` has emitted only half the equalities.

In transport coordinates the invariants are the direction cosines
themselves: :math:`p_1 = \mu_x = \eta` (the radial cosine),
:math:`p_2 = \mu_z = \mu` (the axial cosine), and the eliminated
:math:`p_3 = \mu_y^2 = \xi^2`. So the orbit space of the shipped
cylindrical fold is

.. math::

   S^2/\langle\sigma_y\rangle
   \;=\; \{\, (\eta, \mu) : \eta^2 + \mu^2 \le 1 \,\},
   \qquad \xi^2 = 1 - \eta^2 - \mu^2 \ \text{recovered from it}.

⚠ **The dimension does NOT drop**, and that one line is the source of
everything in the rest of this section. :math:`H` is **finite**, so
:math:`\dim H = 0` and :math:`\dim(S^2/H) = 2 - 0 = 2`; generic orbits
are two points, not curves. Contrast :math:`S^2/SO(2)_a`, where
:math:`2 - 1 = 1`. With no dimensional reduction the invariant chart
buys nothing as a *data* representation — :math:`3 \to 2` floats with
the third recoverable — while for :math:`SO(2)` it buys a genuine
:math:`3 \to 1` reduction (:ref:`manifold-chart-section-asymmetry`).

.. _manifold-stratum-is-a-locus:

The stratum is a LOCUS, and that retyped a field
-------------------------------------------------

:math:`\det P = 4 x_a^2` vanishes exactly on the mirror's own
fixed-point set — the great circle :math:`S^2 \cap \{x_a = 0\}`, which
in the realization's coordinates is the disk's **boundary**:

.. math::

   \operatorname{Fix}(\sigma_a) \cap S^2
   \;=\; \{\, \xi = 0 \,\}
   \;\longleftrightarrow\;
   \{\, \eta^2 + \mu^2 = 1 \,\} \;=\; \partial D^2 .

Every point of it is fixed by :math:`\sigma_a`, so its stabilizer is
all of :math:`H`, the orbit collapses from two points to one, and the
quotient is an **orbifold with boundary** — the same conclusion as
:ref:`manifold-singular-stratum` reaches for :math:`SO(2)`, by the same
route. In transport terms the stratum is :math:`\xi = \mu_y = 0`: the
purely **meridional** directions.

⛔ **And this is what retyped a field.** ``singular_stratum`` was
``tuple[float, ...]`` and the first entry stored ``(-1.0, 1.0)``. A
**circle is not a tuple of floats**. The first catalogued entry's
*shape* had silently become the field's *type*: a stratum is a locus,
and two poles are merely a locus that happens to be finite. The field
is now a SymPy expression in the realization's coordinates whose
vanishing set is the stratum, with ``None`` for a free action — so
`[M]` :math:`SO(2)` stores ``1 - u0**2`` (solving to
:math:`\{-1,+1\}`, unchanged in content) and :math:`\sigma_a` stores
``1 - u0**2 - u1**2``, whose solution set is a one-parameter family.
:attr:`is_free <orpheus.numerics.manifold.Quotient.is_free>` reads
``is None`` rather than ``== ()``.

⭐ **Why the stratum is STORED at all, when** ``det_gram`` **already
determines it.** This looks like a Pattern-2 twin and is not, and the
distinction is worth the paragraph because it is the criterion for
every future field: *a value that its owner cannot recompute from its
own state is derivation output, and storing it is right.* Recovering
the locus needs the **base's** defining ideal —
:math:`\det P = 4p_2` becomes :math:`4(1-\mu^2)` only after
substituting :math:`p_1^2 + p_2 = 1`, and :math:`\det P = 4p_3` becomes
:math:`4(1 - p_1^2 - p_2^2)` only after
:math:`p_1^2 + p_2^2 + p_3 = 1` — and a
:class:`~orpheus.numerics.manifold.Quotient` does not carry that ideal.
Step 5 of the procedure is exactly where that ideal enters, which is
why the stratum is a *fifth-step* output and not a property of
:math:`\det P` alone.

.. _manifold-two-coordinate-systems:

Two honest coordinate systems: a chart codomain AND a section's image
----------------------------------------------------------------------

The ruled shape (user, 2026-08-31) is **two slots**, and the two answer
different questions. State them that way when populating a new entry,
because the failure mode is putting the right object in the wrong slot:

.. list-table:: The two slots, and the question each answers
   :header-rows: 1
   :widths: 22 39 39

   * -
     - ``realization``
     - ``fundamental_domain``
   * - The question
     - *"What does the invariant chart of* :math:`M/H` *map ONTO?"*
     - *"Which points of* :math:`M` *are the chosen orbit
       representatives?"*
   * - Whose coordinates
     - the **invariants'** — the same language as ``generators``,
       ``gram`` and ``det_gram``
     - the **base's** ambient coordinates
   * - :math:`S^2/SO(2)_a`
     - ``Interval(-1, 1)`` — the polar cosine :math:`\mu`
     - ``None``: no section is canonical
   * - :math:`S^2/\langle\sigma_a\rangle`
     - ``Ball(2)`` — the disk :math:`(\eta, \mu)`
     - ``FundamentalDomain(SPHERE, ((0.0, 1.0, 0.0),), 'y>=0')`` — the
       closed hemisphere, named ``S^2|y>=0``
   * - Who produces data in it
     - a rule *born* in the chart, e.g.
       ``gauss_legendre(8).measure.nodes``, shape :math:`(8,)`
     - :meth:`DiscreteMeasure.quotient
       <orpheus.numerics.measure.DiscreteMeasure.quotient>`, **always** —
       e.g. ``folded_product(4,8).measure.nodes``, shape :math:`(16,3)`
   * - Can it see the ERR-080 forgery?
     - ⛔ **no** — Mode-12 blind, measured below
     - ✅ yes — the forged rows are not on :math:`S^2`

⭐ **The producer had already chosen, and it chose the section.** `[M]`
:meth:`DiscreteMeasure.quotient
<orpheus.numerics.measure.DiscreteMeasure.quotient>` computes orbit
representatives and then pushes forward along
``lambda nodes: nodes[representative]`` — a **selection** of parent
nodes, applying no chart. So *every* measure the tree emits through
``.quotient(...)`` carries the base's ambient columns, by construction;
``folded_product``'s :math:`(16,3)` is not a stylistic choice but the
only thing that method can produce. Under a chart-only reading,
``Quotient.contains`` could validate **none** of them.

**What each half of the type does.**

- ``Quotient.contains`` accepts **either** language and dispatches on
  the ambient width. This is the one place in the type where the
  distinction is a genuine local split rather than a repeated tag test,
  and it is the reason the fork is resolved in one method instead of at
  every call site.
- ``_ambient`` still reports the **realization's** width. That is
  deliberate and is not a compromise: a
  :class:`~orpheus.numerics.manifold.Product` factor must have one
  canonical width, or a product's coordinate split would be ambiguous.
  ``contains`` is deliberately the wider of the two.
- ``Quotient.__post_init__`` gates the pair: the two views describe one
  object, so their ``dim`` must agree. `[M]` this is a real check and
  not a tautology — the fundamental domain *derives* its ``dim`` from
  the base less one per antipodal normal pair, while the realization
  *states* its own, so a hemisphere offered against a 1-D realization
  is refused where it is written:

.. code-block:: text

   S^2/sigma_y: the fundamental domain 'S^2|y>=0' has dim 2 but the
   realization '[-1,1]' has dim 1 — the two must describe the same
   orbit space. Check the normals: an antipodal PAIR spells an
   equality and drops a dimension; a lone normal does not.

.. warning::

   ⛔ **The half-spaces must be CLOSED, and the witness is production
   data.** `[M]` the cylindrical march seeds each polar level at
   :math:`\xi = 0` **exactly** — the seed of level :math:`p` is
   :math:`(-\sqrt{1-\mu_p^2},\, 0,\, \mu_p)`, on :math:`S^2` to
   :math:`0.0` and on the stratum to :math:`0.0` — so a strict
   :math:`\langle p, n\rangle > 0` would refuse a direction the
   production sweep actually marches from. This is
   ``coding-elegance`` anti-pattern #18's half (ii) — *every legal
   value must be admitted*, which is a claim about the **producers**
   and is the half that gets skipped. Gated by
   ``test_the_half_space_is_CLOSED_because_production_marches_from_it``.

   ⚠ The shipped folded *rule* cannot witness this: `[M]` its 16 nodes
   have :math:`\mu_y \in [0.1945,\, 0.8688]`, strictly positive,
   because the even-:math:`n_\varphi` staggering makes the fold **free**
   on the nodes. The closure's edge data is the only witness available,
   which is why the gate is built on the seeds rather than on the
   quadrature.

.. _manifold-realization-refuted:

Five single-slot candidates, measured and refused
--------------------------------------------------

Before the two-slot ruling, five single-object candidates were put to
the shipped data. **All five fail**, and they fail in two disjoint
ways, which is itself the argument for two slots: the two that admit
the tree's nodes are blind to the chart, and the three that admit chart
points are blind to the nodes. The matrix below is **measured**, cell
by cell — and
"REFUSE (shape)" is a raised :exc:`ValueError` from the ambient-width
check, not a ``False`` return, which is a behavioural difference a
caller must handle (:ref:`manifold-gotcha-shape-vs-verdict`).

The five inputs: the shipped folded nodes :math:`(16,3)`; the **orbit
twins**, the same nodes with :math:`\mu_y \to -\mu_y`, i.e. the wrong
representative; the **ERR-080 forgery** :math:`(\mu, 0, 0)`,
:math:`(8,3)`, not unit-norm; and the chart images of the first and
third, :math:`(16,2)` and :math:`(8,2)`.

.. list-table:: Every candidate against every input, measured 2026-08-31
   :header-rows: 1
   :widths: 26 15 15 15 15 14

   * - Candidate
     - shipped
     - twins
     - forgery
     - shipped charted
     - forgery charted
   * - ``SPHERE``
     - ADMIT
     - ⛔ **ADMIT**
     - refuse
     - REFUSE (shape)
     - REFUSE (shape)
   * - ``RealSpace(2)``
     - REFUSE (shape)
     - REFUSE (shape)
     - REFUSE (shape)
     - ADMIT
     - ⛔ **ADMIT**
   * - ``COSINE_INTERVAL × COSINE_INTERVAL``
     - REFUSE (shape)
     - REFUSE (shape)
     - REFUSE (shape)
     - ADMIT
     - ⛔ **ADMIT**
   * - ``Ball(2)`` alone
     - REFUSE (shape)
     - REFUSE (shape)
     - REFUSE (shape)
     - ADMIT
     - ⛔ **ADMIT**
   * - the hemisphere alone
     - ADMIT
     - refuse
     - refuse
     - REFUSE (shape)
     - REFUSE (shape)
   * - ⭐ **SHIPPED: both slots**
     - ADMIT
     - refuse
     - refuse
     - ADMIT
     - ⚠ ADMIT

Reading the rows:

**The sphere itself.** ``realization = SPHERE`` is the convenient
placeholder, and the one to refuse hardest. It does *not* regress
ERR-080 (the forgery is still
refused, since :math:`\lVert(\mu,0,0)\rVert = |\mu| \ne 1`). What it
loses is **the fold itself**: the orbit twins are admitted, and more
sharply, `[M]` under this choice ``Quotient.contains`` becomes the
*same function* as ``SPHERE.contains`` — **no input exists that the
quotient refuses and its base admits**. A predicate that cannot
distinguish :math:`M/H` from :math:`M` is ``vv-principles`` #17's *arm
with no witness*, decidable at design time with no mutation needed.
⛔ And it is topologically false, not merely weak: :math:`D^2 \ncong
S^2` — the disk is contractible with boundary :math:`S^1` and
:math:`\chi = 1`, the sphere has no boundary and :math:`\chi = 2`. That
``dim`` happens to agree (both 2) carries no information here either,
precisely because :math:`H` is finite.

**The two shipped 2-D members**, ``RealSpace(2)`` and the square. Both
buy nothing the disk does not, and both are strictly weaker:
``RealSpace(2)`` drops the disk inequality entirely, and
``COSINE_INTERVAL × COSINE_INTERVAL`` is the bounding **square**, whose
discriminating witness is measured — `[M]` :math:`(0.9, 0.9)` is in the
square and not in the disk, and corresponds to **no direction at all**,
since :math:`\eta^2 + \mu^2 = 1.62 > 1` forces :math:`\xi^2 = -0.62`.
Reusing a shipped member because it ships is how a type acquires a
predicate that admits impossible points.

**The disk alone.** ``Ball(2)`` is what the *documented* meaning of
``realization`` requires, and the sharpest refusal of the five, because
it fails on the
very defect the type was minted for. ⛔ `[M]` **the chart is Mode-12
blind to ERR-080.** The chart :math:`(x,y,z) \mapsto (x,z)` discards
:math:`\mu_y`, and :math:`\mu_y` is precisely what the forgery
falsifies, so the forged row :math:`(\mu, 0)` is a **perfectly legal**
point of the disk — it is the orbit of the real direction pair
:math:`(\mu, \pm\sqrt{1-\mu^2}, 0)`. Measured:
:math:`\max |(\mu,0)|^2` over the eight forged rows is
:math:`0.9221566084920586 < 1`, so *every* forged row lands inside the
closed disk. The mechanism is exact, not statistical
(``vv-principles`` Mode 12 — the measured functional's stabiliser
contains the error class): no tolerance, refinement or fixture choice
can expose it.

**The hemisphere alone** — the only single object that admits the
shipped nodes and refuses both wrong inputs. Its cost is that it
**redefines** the field rather than adding to it: ``realization`` would
stop meaning *chart codomain*, the type would then answer in the base's
language for :math:`\sigma_a` and in the chart's for :math:`SO(2)` —
the exact vocabulary drift the mint exists to end
(:ref:`manifold-string-drift`) — and the derivation fields
(``generators``, ``gram``, ``det_gram``) would be speaking a coordinate
system the ``realization`` beside them no longer names.

⚠ **And read the shipped row honestly: it does not dominate every
cell.** The two-slot design still admits the *charted* forgery, and
that is correct rather than a residual defect — in chart coordinates
:math:`(\mu, 0)` genuinely **is** a point of the orbit space, and no
predicate over the chart can know it was built by zero-padding. What
the second slot buys is that the data the tree actually produces
arrives in *section* coordinates, where the predicate that can see the
forgery is the one that runs.

.. _manifold-err-080-is-a-section:

ERR-080's level-1 half is a botched section of :math:`S^2/SO(2)_x`
------------------------------------------------------------------

The chart-versus-section question is not new with :math:`\sigma_y`. It
arises for :math:`SO(2)` too — the moment any consumer of a
1-dimensional rule needs a 3-D direction — and the tree has been
answering it, silently and wrongly, for as long as :ref:`ERR-080
<manifold-err-080>` has existed.

The realization :math:`[-1,1]` is the **chart's** codomain,
unambiguously: a section of :math:`S^2 \to S^2/SO(2)_a` is a half-meridian
*inside* :math:`S^2 \subset \mathbb{R}^3`, ambient 3, and
``Interval(-1, 1)`` is ambient 1. So when
:meth:`Quadrature.angular_frame
<orpheus.numerics.quadrature.directional.Quadrature.angular_frame>`
needs three columns it is not asking for the chart at all — it is
asking for a **section**, which the tree never had. It fabricated one
by zero-padding:

.. list-table:: The fabricated section against an honest one
   :header-rows: 1
   :widths: 40 60

   * - Construction
     - `[M]` 2026-08-31, on ``gauss_legendre(8)``
   * - what the tree builds — ``column_stack`` of the three
       axis-cosine arrays, two of them a zero *fallback*
     - rows :math:`(\mu, 0, 0)`; norms
       :math:`0.1834 \ldots 0.9603`; ``Sphere().contains`` → ``False``
   * - an honest :math:`\varphi = 0` half-meridian. ⛔ this row spelled
       it :math:`\mu \mapsto (\sqrt{1-\mu^2},\, 0,\, \mu)` — a
       :math:`z`-pole section — until 2026-09-01. With the axis a
       parameter it is written in the axis's own language,
       :math:`\mu \mapsto \mu\,\hat e_a + \sqrt{1-\mu^2}\,\hat e_b`, and
       the slab's :math:`a` is :math:`x`:
       :math:`\mu \mapsto (\mu,\, \sqrt{1-\mu^2},\, 0)`
     - `[M]` on :math:`S^2` to :math:`0.0` (max
       :math:`\bigl|\lVert v\rVert - 1\bigr|`);
       ``Sphere().contains`` → ``True``. ⭐ Note what the comparison now
       shows: the fabrication is this map with the :math:`\hat e_b` term
       **dropped**, which is precisely why its rows fall short of the
       unit sphere
       (:ref:`manifold-the-axis-convention-for-a-section`)

⟹ **ERR-080's first link is not "a wrong tag". It is a section
fabricated where none was declared** — the realization is a chart, a
consumer needed a section, and zero-padding is what a codebase does
when the object it needs has no slot. With ``fundamental_domain`` in
the type, that padding has nowhere to live: an entry either declares a
section or declares that it has none, and :math:`S^2/SO(2)_a` honestly
declares ``None``.

.. warning::

   ⚠ **This names the level-1 half only. Do not read it as the
   repair.** Making the section land on :math:`S^2` makes the nodes
   *points of the manifold*; it does **not** fix the level-2 half. On
   any :math:`\varphi = 0` section every :math:`Y_\ell^{m \ne 0}` is
   evaluated at a *chosen* azimuth that carries no information, and the
   corpus's recorded repair for that is unchanged: the **trivial
   isotypic sub-basis** :math:`\{Y_\ell^0\} \cong \{P_\ell\}`
   (:ref:`manifold-gelfand`, and the falsifiable form with ERR-080 in
   the :doc:`error catalogue </theory/verification/error_catalog>`).
   Both halves are needed; this section establishes only the first, and
   **ERR-080 remains open** (:ref:`manifold-seams`).

   ⚠ Declaring a section for :math:`SO(2)` would also be a **choice**,
   not a derivation — the :math:`\varphi = 0` half-meridian is one of a
   continuum — so it belongs to the step that makes a slab declare its
   quotient, not to the derivation. The shipped entry therefore carries
   ``fundamental_domain=None`` on purpose.

.. _manifold-chart-section-asymmetry:

Why the two entries diverge: a structural asymmetry, not a style choice
-------------------------------------------------------------------------

This is the transferable half of the whole section, and it is what a
future entry's author needs before populating either slot. The two
shipped entries make **opposite** choices, and both are locally
correct; what could not serve both was the *type*.

.. list-table:: The asymmetry, term by term
   :header-rows: 1
   :widths: 26 37 37

   * -
     - :math:`S^2/SO(2)_a`
     - :math:`S^2/\langle\sigma_a\rangle`
   * - the group
     - compact **connected**, :math:`\dim = 1`
     - **finite**, :math:`\dim = 0`, and a *reflection*
   * - :math:`\dim(S^2/H)`
     - :math:`2 - 1 = 1` — **drops**
     - :math:`2 - 0 = 2` — **does not drop**
   * - the invariant chart's codomain
     - :math:`[-1,1] \subset \mathbb{R}^1`
     - the closed disk :math:`D^2 \subset \mathbb{R}^2`
   * - the chart *as data*
     - :math:`3 \to 1` floats: a genuine **reduction**
     - :math:`3 \to 2` floats, the third recoverable: **no reduction**
   * - a canonical section?
     - ⛔ **no** — every half-meridian is one; :math:`\varphi = 0` is an
       arbitrary pick
     - ✅ **yes** — the mirror determines the closed half-space, and
       being a *reflection* makes it **strict**
   * - what the tree ships as ``measure.nodes``
     - the **chart**, :math:`(8,)`
     - the **section**, :math:`(16,3)`
   * - `[M]` ``measure.support`` (2026-09-02)
     - ``'S^2/SO2_x'`` — the *quotient's* name. ⛔ This cell read
       ``'[-1,1]'`` — the *realization's* name — until 2026-09-02; it
       was true when written and tracker 2.4 repealed it on
       2026-09-01, when the slab's rule began *declaring* its orbit
       space (:ref:`manifold-orbit-space-declaration`). The row's
       point survives intact: what the rule ships as ``nodes`` is
       still the **chart**, one column, and only the *support* now
       says which orbit space those chart coordinates are for.
     - ``'S^2/sigma_y'`` — the *quotient's* name

⟹ **For a positive-dimensional group the chart is strictly cheaper and
no section is canonical, so the chart wins and the tree ships it. For a
finite reflection the chart is no cheaper and the section IS canonical,
so the section wins and the tree ships that.** Neither branch was
wrong; the single-slot type was.

.. warning::

   ⛔ **"Canonical because it is a reflection" does not generalise to a
   rotation — leave** ``fundamental_domain=None`` **for** :math:`C_n`.
   The closed half-space is a *strict* fundamental domain (it meets
   every orbit exactly once) only because :math:`\sigma_a`'s
   fixed-point set lies **in** it: a free orbit
   :math:`\{(x,y,z), (x,-y,z)\}` with :math:`y > 0` meets
   :math:`\{y \ge 0\}` once, and a stratum orbit :math:`\{(x,0,z)\}` is
   a single point that also lies in it. For a rotation :math:`C_n` the
   closed sector's two meridian edges are identified **with each
   other** by the group, so the closed sector maps 2-to-1 on its
   boundary and is *not* homeomorphic to the orbit space. A
   fundamental-domain slot filled for a :math:`C_n` entry would be
   stating something false, and the type cannot catch it — the ``dim``
   gate would pass.

.. note::

   **The hemisphere IS a legitimate realization set-theoretically, and
   is NOT a diffeomorphic one — both halves matter.** `[M]` the chart
   :math:`c : H^+ \to D^2`, :math:`(x,y,z) \mapsto (x,z)`, with inverse
   :math:`(p_1,p_2) \mapsto (p_1, \sqrt{1-p_1^2-p_2^2}, p_2)`, is a
   continuous bijection from a **compact** set onto a **Hausdorff**
   one, hence a homeomorphism (no separate inverse-continuity argument
   needed). It is *not* a diffeomorphism: :math:`\partial y/\partial
   p_i = -p_i/\sqrt{1-p_1^2-p_2^2}` blows up on the stratum, and from
   the forward side :math:`\mathrm{d}c` annihilates :math:`e_y` there —
   rank 1 on the boundary circle, 2 in the interior. That is a
   **Whitney fold**, and it shows up in the measure as an integrable
   singularity in :math:`1/\lvert y \rvert` — that is, in
   :math:`1/\lvert\xi\rvert`, the coordinate the fold quotients:
   `[M]` :math:`\int_{D^2} \mathrm{d}p_1\, \mathrm{d}p_2 /
   \lvert y\rvert = 2\pi`, the area of a hemisphere — finite.

   ⚠ **The fold does not bite** ``contains``. Membership is a level-1,
   set-theoretic question and the homeomorphism settles it; the fold
   bites at **levels 2 and 3** — what a basis function eats, and what a
   derivative operator differentiates. So it must not be cited as an
   argument against a disk realization *for membership purposes*. The
   arguments against the disk alone are the ones in
   :ref:`manifold-realization-refuted`, and they are entirely
   different.

.. _manifold-orbifold-discretised:

The orbifold is already realized in the shipped cylindrical sweep
------------------------------------------------------------------

:ref:`manifold-singular-stratum` records that two shipped objects live
on the :math:`SO(2)` stratum from opposite directions. The
:math:`\sigma_y` entry has an exact cylindrical analogue, and it is
measured:

.. list-table:: The fold's stratum against the shipped cylindrical data
   :header-rows: 1
   :widths: 46 54

   * - Object
     - `[M]` 2026-08-31, ``folded_product(4, 8)`` on a cylinder
   * - the 16 quadrature nodes
     - :math:`1 - \eta^2 - \mu^2 \in [0.0378,\, 0.7549]` — **strictly
       interior**; the fold is free on them
   * - the four march seeds — the starting angular edge per level,
       ``AngularRedistribution.mu_start_per_level``
     - the seed direction is
       :math:`(-\sqrt{1-\mu_p^2},\, 0,\, \mu_p)`, on :math:`S^2` to
       :math:`0.0` and with
       :math:`1 - \eta^2 - \mu^2 = 0.0` **exactly**, on all four —
       i.e. **on the stratum**
   * - the azimuthal cell edges of each level
     - the nodes sit at :math:`\omega/\pi \in \{0.125, 0.375, 0.625,
       0.875\}`, the staggered midpoints of four cells partitioning
       :math:`(0,\pi)`, so the edges are
       :math:`\omega/\pi \in \{0, \tfrac14, \tfrac12, \tfrac34, 1\}`
       and the two **outer** edges are :math:`\omega = 0, \pi`, where
       :math:`\xi = \sin\theta \sin\omega = 0`
   * - the :math:`\alpha`-dome per level
     - five edge values per level, closing at both ends. Levels 0 and 3
       read :math:`[0,\, 0.2566,\, 0.3629,\, 0.2566,\, 0]`; levels 1
       and 2 read
       :math:`[0,\, 0.8900,\, 1.2587,\, 0.8900,\, 1.1\times10^{-16}]`

⟹ **That is what an orbifold looks like when you discretise it:** the
*interior* of the fundamental domain carries the nodes, and its
*boundary* — the stratum — carries the closure's degenerate data, the
:math:`\alpha`-dome's zeros and the march seed.

⚠ Two naming traps in that table, both `[M]` and both worth knowing
before quoting it. ``mu_start_per_level`` holds a **radial** cosine
:math:`\eta = -\sin\theta_p`, not a polar :math:`\mu` — the name is the
half-angle thread's, not this page's; and the field's own docstring
spells the level's polar cosine :math:`\xi_p`, while :math:`\xi`
everywhere on this page is :math:`\mu_y`, the *azimuthal* cosine that
the fold quotients. The values are unambiguous — `[M]`
``mu_start_per_level`` equals :math:`-\sqrt{1-\mu_p^2}` on the level
cosines exactly — but the symbols are not.

⚠ **This is an identity of LOCI, established by arithmetic — the
mechanism is unproved.** Whether the cylindrical redistribution term
*is* the :math:`\sigma_y` quotient's connection is the exact
cylindrical twin of the open spherical question recorded in the third
row of :ref:`manifold-one-polynomial`, and it is not closed here. It
has the same two admissible outcomes: a derivation, or an explicit
ruling that the coincidence is geometric bookkeeping. Cite this
subsection for the measured coincidence, never for the mechanism.


.. _manifold-engine-seed:

The catalogue is the engine's SEED, not its rival
=================================================

An orbit space *can* be computed from scratch: steps 1–4 above are
mechanical, and a Gröbner-basis engine would run them. The project has
ruled that it will **not build that engine yet** — and the ruling is
worth quoting exactly, because the obvious paraphrase ("we refused the
engine") is the wrong one:

   *"We're not outright ruling out building the engine. We're ruling
   that we will not prematurely build the engine. The embryo should be
   such that if the day ever arises that we decide building the engine
   is the right step, it should be a development of the embryo, instead
   of having to do the entire engine from scratch for a code base that
   was not ready to receive it."* — user, 2026-08-31 (decision D0.1)

The groups that occur in transport number about a dozen. A Gröbner
engine for them is abstraction without a consumer, and its failure mode
is debugging elimination orderings instead of transport. So each entry
is derived once by the procedure and recorded — **deferred, not
refused.**

.. _manifold-engine-data-model:

The binding requirement is on the DATA MODEL, not the signature
----------------------------------------------------------------

⛔ The first version of this ruling said *"the catalogue and the engine
have the same interface — a second backend behind an unchanged
signature"*, and it was **rejected as too weak**, in the ruling's own
words as *"the twin-path risk wearing a compliment"*. A shared
signature guarantees only that the *call site* survives. The engine
would still arrive with its own representation of polynomials, ideals
and PSD conditions, plus a translation layer to whatever the catalogue
happened to store — a from-scratch build with a seam, which is exactly
what the ruling forbids.

⟹ **A catalogue entry must BE the derivation procedure's output, not a
human summary of its answer.** The procedure emits, per entry: the
invariant generators; the syzygy ideal; the matrix :math:`P` and
:math:`\det P`; the chart; a section, when one is canonical; the
pushforward measure; the stratum where :math:`\det P` vanishes.
**Those are the entry's fields.** An engine then ships by *computing*
them instead of reading them — a development, with no new vocabulary
and no seam.

⭐ **The list above grew by one, and the growth is the ruling working
rather than the ruling slipping.** The section was not on the
procedure's output list until the second entry produced one
(:ref:`manifold-two-coordinate-systems`); it was added to the
*procedure*, and a slot was added to match, which is exactly the
direction the check below permits. What it forbids is the reverse — a
field the procedure does not emit, or an output the procedure emits
that the entry has to summarise in prose.

The ruling comes with its own falsifiable check, and it is the question
to ask of any future edit here:

   *Given a catalogue entry, could an engine populate its fields
   without introducing a single new type?*

If the answer is no, the embryo has drifted from being a seed and the
ruling has been violated — however clean the interface looks.

.. list-table:: The procedure's outputs against the shipped slots
   :header-rows: 1
   :widths: 34 22 44

   * - Procedure output
     - Slot on :class:`~orpheus.numerics.manifold.Quotient`
     - Note
   * - invariant generators :math:`p_1 \ldots p_k`
     - ``generators``
     - SymPy expressions in the ambient coordinates
   * - syzygy ideal :math:`I`
     - ``syzygy``
     - ``()`` when the invariants are independent
   * - :math:`P_{ij}`
     - ``gram``
     - re-expressed in the invariants
   * - :math:`\det P`
     - ``det_gram``
     - its zero locus is the orbit-space boundary
   * - the singular stratum
     - ``singular_stratum``
     - a **locus** in the realization's coordinates — derivation
       output that a ``Quotient`` cannot recompute, because recovering
       it needs the base's own ideal
       (:ref:`manifold-stratum-is-a-locus`)
   * - provenance
     - ``derived_by``
     - ``"hand"`` / ``"engine"``
   * - a section of :math:`M \to M/H`, when canonical
     - ``fundamental_domain``
     - its **image**, in the base's coordinates; ``None`` is an
       answer, not a gap
       (:ref:`manifold-two-coordinate-systems`)
   * - the chart :math:`M/H \to N`
     - ⛔ **not a slot**
     - only its **codomain** ships, as ``realization`` — the mirror of
       the row above, where only the *image* ships. Neither **map** is
       a field. ⭐ **Narrowed 2026-09-02 (tracker 2.3):** a map is now
       a *type* (:class:`~orpheus.numerics.manifold.ManifoldMap`), so
       this row is a missing **slot** rather than a missing
       vocabulary. `[M]` none of the three arrows the tree ships is
       this chart or that section — the reasons are enumerated per
       arrow at :ref:`manifold-arrows-not-built`.
   * - the pushforward measure :math:`\pi_*\mu`
     - ⛔ **not a slot**
     - the measure descends today via
       :meth:`DiscreteMeasure.quotient
       <orpheus.numerics.measure.DiscreteMeasure.quotient>`, which
       knows nothing about this catalogue. ⚠ **And it cannot be added
       by importing** — `[M]` 2026-09-02, a module-scope
       ``manifold → exactness`` edge closes a two-hop cycle and kills
       5 of 5 fresh import orders; the viable mechanism is a
       function-scope import inside the derivation function
       (:ref:`manifold-arrows-not-built`, item 2). Tracker **3.1**.

`[M]` **7 of 9** — by ``dataclasses.fields``, ``Quotient`` declares ten
fields: ``base``, ``by``, ``realization``, ``fundamental_domain``,
``generators``, ``syzygy``, ``gram``, ``det_gram``, ``derived_by``,
``singular_stratum``, of which the first two are the entry's *inputs*.
So the seed is real for seven of the procedure's nine outputs and
**incomplete for two**, and those two are named as seams rather than
left to be discovered (:ref:`manifold-seams`). Stating the fraction is
the point: a ruling whose compliance is claimed but not counted is not
checkable. (`[M]` it read **6 of 8** until the two-slot ruling on
2026-08-31; both the numerator and the denominator moved, which is why
the two numbers are given together and not as a percentage.)

⭐ **Why the provenance field exists at all.** ``derived_by`` is read by
nothing today, and a reviewer could reasonably call it speculative. It
is not: a mixed hand/engine state must be *expressible*, or an
incremental engine rollout would have to be all-or-nothing — and an
incremental rollout is exactly what the ruling anticipates. The field
is the difference between a migration that can be staged and one that
cannot.

.. _manifold-refusal-names-the-work:

The refusal names the missing WORK, not the gap
------------------------------------------------

The engine's absence must be a work item a fresh session can pick up,
not a wall. `[M]` ``SPHERE.quotient(SubgroupOfO3.OctahedralOh)``
raises, verbatim:

.. code-block:: text

   no catalogue entry for S^2/Oh: derive it once (minimal invariants
   p_1..p_k of R[x]^H; syzygy ideal I by elimination; Procesi-Schwarz
   P_ij = <grad p_i, grad p_j> with P >= 0; intersect with the ideal of
   S^2) and register it in orpheus/numerics/manifold.py's
   _ORBIT_CATALOGUE, or implement the derivation engine. Catalogued
   today (manifold CLASS / group): ['Sphere/SO2_x', 'Sphere/SO2_y',
   'Sphere/SO2_z', 'Sphere/sigma_x', 'Sphere/sigma_y',
   'Sphere/sigma_z'].

Four things are in that message and all four are load-bearing: which
pair was asked for, the **procedure** to answer it, **where** to put
the answer, and what is catalogued already. The last is spelled
*manifold CLASS* on purpose — the request is named by manifold
*instance* name (``S^2``) while the catalogue is keyed by class, and a
message that silently switched vocabularies would send its reader
looking for a key that does not exist.

The catalogue is keyed on the **pair** ``(manifold class, group
name)``, because a quotient is binary dispatch: it is a property of
neither operand alone.

.. _manifold-tests-are-the-spec:

The regression tests are the engine's acceptance suite, written first
---------------------------------------------------------------------

Every catalogue entry ships a **symbolic** regression test that
reproduces its own derivation. Because the fields *are* the
procedure's outputs, those tests are not merely regression pins: they
are the engine's **specification**, and the engine ships on the day it
reproduces them by computation instead of by lookup. A specification
written before the implementation cannot be shaped to flatter it
(``vv-principles`` #17).

Concretely, for the two shipped entries, the assertions are on the
**symbolic value** rather than on a string or a float:

.. code-block:: python

   # tests/numerics/test_manifold.py::TestQuotient
   assert s2_mod_so2.syzygy == ()
   assert sp.simplify(s2_mod_so2.det_gram - 4 * p2) == 0
   assert sp.simplify(
       s2_mod_so2.gram - sp.Matrix([[1, 0], [0, 4 * p2]])
   ) == sp.zeros(2, 2)

   # ...and the stratum is DERIVED, not compared to a literal:
   det_on_sphere = sp.simplify(s2_mod_so2.det_gram.subs(p2, 1 - mu**2))
   roots = sorted(sp.solve(sp.Eq(det_on_sphere, 0), mu))
   assert [float(r) for r in roots] == [-1.0, 1.0]
   assert sp.simplify(s2_mod_so2.singular_stratum - (1 - u0**2)) == 0

   # tests/numerics/test_manifold.py
   #   ::TestTheSigmaYFoldIsExpressibleAndDiscriminating
   assert sp.simplify(fold.det_gram - 4 * u2) == 0
   assert sp.simplify(fold.gram - sp.diag(1, 1, 4 * u2)) == sp.zeros(3, 3)
   assert fold.realization == Ball(2)
   assert fold.dim == 2                     # H is FINITE: no drop

   # ...and the stratum is shown to be a CURVE, not a point set:
   assert sp.simplify(fold.singular_stratum - (1 - u0**2 - u1**2)) == 0
   sols = sp.solve(sp.Eq(fold.singular_stratum, 0), u1)
   assert len(sols) == 2 and all(u0 in s.free_symbols for s in sols)

⭐ **The last two lines are the shape of the assertion that a retyped
field needed.** Solving the locus and checking that the solutions
*retain a free symbol* is what distinguishes a curve from a point set,
and it is a claim no comparison against a literal could make. The
:math:`SO(2)` gate's own stratum assertion was likewise rewritten to
solve rather than compare — `[M]` it read
``s2_mod_so2.singular_stratum == (-1.0, 1.0)`` until 2026-08-31 — and
it survived the retyping **without weakening**, because it had already
been written to solve :math:`\det P = 0` rather than to trust the
stored value.

.. _manifold-twin-lookup:

⭐ The tree already performs this lookup, one level up
------------------------------------------------------

The mint was not introducing a new idea. `[M]`
:attr:`AngularSymmetry.support
<orpheus.numerics.quadrature.registry.AngularSymmetry.support>` — which
predates it — already computed :math:`S^2/G^0` from the *spent* group by
catalogue lookup, in the string vocabulary, and already raised
:exc:`NotImplementedError` on an uncatalogued quotient with the same
shape of message.

✅ **The twin is COLLAPSED as of tracker 2.4 (2026-09-01), in the
direction reading (iv) below predicted.** ``AngularSymmetry.support`` no
longer holds a table: it *calls* ``SPHERE.quotient(spent)``, so the two
lookups are one call and cannot disagree. `[M]` for a slab,
``GEOMETRY_ANGULAR_SYMMETRY['slab'].support is
SPHERE.quotient(SubgroupOfO3.SO2('x'))`` — object **identity**, not
merely equality, because the catalogue memoises
(:ref:`manifold-quotient-is-memoised`). The section is kept as the
argument that got there.

.. list-table:: Two lookups of :math:`S^2/H`, re-measured 2026-09-01 after the tracker-2.4 collapse
   :header-rows: 1
   :widths: 16 30 30 24

   * - :math:`H`
     - ``AngularSymmetry(...).support.name``
     - ``SPHERE.quotient(H).name``
     - Reading
   * - ``SO2('x')``
     - ``'S^2/SO2_x'``
     - ``'S^2/SO2_x'``
     - **Identical object.** ⛔ This row read ``'[-1,1]'`` on both sides
       until 2026-09-01 — the registry returned the CHART, which is the
       axis-blind spelling a slab rule could share with a spatial
       interval.
   * - ``SO2('z')``
     - ``'S^2/SO2_z'``
     - ``'S^2/SO2_z'``
     - Identical object. A row that could not be *spelled* before the
       axis was a parameter.
   * - ``Trivial``
     - ``'S^2'``
     - ``'S^2/Trivial'``
     - ⚠ **The one surviving divergence, and it is deliberate.** The
       registry short-circuits to the base, because a geometry that
       spends nothing discretises the sphere and ``'S^2'`` is the name
       every 2-D/3-D rule declares. The catalogue reports the
       *derivation's own output*, whose ``realization`` **is** that
       sphere. Same point set, two registers; a committed row pins them.
   * - ``sigma_y``
     - ``'S^2/sigma_y'``
     - ``'S^2/sigma_y'``
     - Now answered, by delegation — but see (iii): it is a row the
       registry has no *business* being asked, not one it was missing.
   * - ``Oh``
     - :exc:`NotImplementedError`
     - :exc:`NotImplementedError`
     - The refusal is now literally the same refusal
       (:ref:`manifold-refusal-names-the-work`).

Four readings, all useful — and all four survive the collapse, because
what collapsed is the *implementation*, not the distinction between the
two questions.

**(i)** On the overlapping rows the typed catalogue reproduced the
registry's answer exactly, which was the cheapest available evidence
that the type is a *re-typing* of an existing fact and not a rival one.
The pin is
``tests/numerics/test_manifold.py::TestQuotient::test_the_derived_orbit_space_agrees_with_the_hand_written_table``,
and reading its docstring is worth more than reading this paragraph,
because it records what the collapse did to its own claim class — twice,
in opposite directions.

* ⛔ **The axial row was DEMOTED by the collapse**, and correctly. Once
  the registry derives its domain *through* ``SPHERE.quotient``, asking
  whether the two agree is asking one call whether it equals itself —
  ``coding-standards``' single-sourcing clause exactly. The fix stays;
  the gate's *description* is what had to move. What survives on that
  row is a different and still-real claim: that the registry hands out
  the **orbit space**, carrying its spent group, rather than the chart.
* ⭐ **The same gate was PROMOTED by tracker 2.0c**, with no line of its
  body changing. It used to compare ``mine.name == theirs`` because
  ``theirs`` was a *string* — the strongest claim the string vocabulary
  admitted. Both sides are ``Manifold`` values now, so the assertion is
  **object equality**: not that two producers spell the orbit space the
  same way, but that they produce the same point set. A name gate is
  satisfied by any self-consistent lie; this one is not.
* ✅ **What is load-bearing on every row after the collapse** is the
  hand-written ``expected`` column, which is authored independently of
  both producers: it is a genuine external pin on the Procesi–Schwarz
  derivation, and the one input that could still redden the row is a
  wrong derivation.

The *do not re-baseline* note stands: if it reddens, one of the two is
wrong about a quotient, and which one is a mathematical question, not a
test-maintenance one.

**(ii)** ⛔ **The** ``Trivial`` **row read** :exc:`NotImplementedError`
**in the first version of this page, and it was already false when the
page landed.** Comparing the two lookups is what exposed the gap —
:math:`S^2/\{e\} = S^2` is legal and trivially derivable — and the same
commit that published the table (``fba4205a``) closed it, by
**deriving** the answer rather than tabulating it: the trivial group's
invariant ring is the whole polynomial ring, so :math:`p_i = x_i`,
:math:`P = I`, :math:`\det P = 1`, which vanishes nowhere, hence no
stratum, hence a free action — right vacuously, the only element being
the identity. That doubles as a **positive control on the machinery**:
the procedure reproduces a known answer. The row is corrected here as
history rather than deleted, because a gap reported into the corpus has
the shortest shelf life of anything on a page — the report is what
triggers the repair.

**(iii)** ⭐ **The** ``sigma_y`` **row is not a gap in the registry, and
reading it as one would send its repair in the wrong direction.** The
two lookups quotient by **different halves of the symmetry**.
:attr:`AngularSymmetry.support
<orpheus.numerics.quadrature.registry.AngularSymmetry.support>` is
defined as :math:`S^2/G^0` — the *continuous isotropy* a dimensional
reduction spends — while a mirror is a member of the *discrete
residual* :math:`\Gamma`, the finite half a quadrature must realize as
an ordinate permutation. `[M]` the shipped geometry table says so
directly: ``GEOMETRY_ANGULAR_SYMMETRY["cylinder"]`` is
``continuous_isotropy=Trivial, discrete_residual=Dnh(2)`` — nothing
continuous is spent on a cylinder, so its declared angular domain is
the whole sphere. :meth:`Manifold.quotient
<orpheus.numerics.manifold.Manifold.quotient>` has no such restriction:
it quotients by any subgroup, which is why it can answer
:math:`S^2/\langle\sigma_y\rangle` at all.

**(iv)** Two lookups of one fact is a Pattern-2 twin by construction,
and reading (iii) sharpens what the collapse *is*. The registry's
``support`` is not a rival catalogue; it is the **special case**
:math:`H = G^0`, so the collapse is
``support = base.quotient(G⁰)`` rather than a merge of two tables.

✅ **LANDED, tracker 2.4, 2026-09-01** — and the shipped form differs
from the one predicted here in one instructive respect. This paragraph
predicted ``base.quotient(G⁰).realization.name``, i.e. the *chart's
name*, because at the time ``support`` was a ``str``. What shipped is
``base.quotient(G⁰)`` — the **orbit space itself**, dropping both the
``.realization`` and the ``.name``. That is not a detail: taking the
realization is exactly the axis-blind step
(:ref:`manifold-so2-axis-is-a-parameter`), since all three
:math:`S^2/SO(2)_a` realize onto the *same* interval. The prediction
was correct about which object is the special case and wrong about how
much of it to keep, and the intervening retype (2.0c) is what made the
better answer expressible.

.. warning::

   ⚠ **A live consequence of (iii), re-measured 2026-09-01: still
   latent, and its stated CAUSE is now wrong.** ``admits_domain`` is
   ``measure.support == self.support``. `[M]` the cylinder declares
   ``S^2`` (a :class:`~orpheus.numerics.manifold.Sphere`) while the
   shipped cylindrical rule carries ``folded_product(4,8).measure.support
   == S^2/sigma_y`` (a :class:`~orpheus.numerics.manifold.Quotient`), so
   ``GEOMETRY_ANGULAR_SYMMETRY["cylinder"].admits_domain(...)`` is
   **False** — stage 0 would still reject the tree's own fold.

   ⛔ This block read *"the gate is a **string comparison** … two
   different quotients that the string vocabulary cannot tell apart"*
   until 2026-09-01. Both halves are void: tracker 2.0c made it a
   ``Manifold`` value comparison, and the two quotients are now
   perfectly distinguishable — which is the point. **The mismatch was
   never the vocabulary; it was the claim.** A rule folded by a member
   of :math:`\Gamma` lives on :math:`S^2/\Gamma'` while the geometry
   declares :math:`S^2/G^0`, and those are two genuinely different orbit
   spaces. Typing them made the disagreement *legible* instead of
   removing it, which is the correct outcome and the reason the fix is
   still not to loosen the comparison.

   It bites nothing today for one reason only: `[M]` ``folded_product``
   is **not in** ``quadrature_registry`` (four specs ship —
   ``GaussLegendre1D``, ``LebedevSphere``, ``LevelSymmetricSN``,
   ``ProductQuadrature``), so the selector never presents it to stage 0.
   The day it is registered, this is the first thing that fires.
   Recorded as a seam (:ref:`manifold-seams`).


.. _manifold-orbit-space-declaration:

The first production consumer: a quadrature DECLARES its orbit space
=====================================================================

Everything above is about the type. This section is about the day it
was first *used*, which is tracker 2.4 (2026-09-01): the slab's polar
quadrature stopped naming the interval a chart happens to map onto and
started naming the orbit space it lives on.

.. list-table:: The slab's polar rule, before and after
   :header-rows: 1
   :widths: 30 32 38

   * -
     - Before (the chart)
     - After, `[M]` 2026-09-01
   * - ``measure.support.name``
     - ``'[-1,1]'``
     - ``'S^2/SO2_x'``
   * - ``measure.space.name``
     - ``'L2[[-1,1]]'``
     - ``'L2[S^2/SO2_x]'``
   * - ``measure.quotient_group``
     - ``None``
     - ``SubgroupOfO3.SO2('x')``
   * - ``measure.phase``
     - ``'angular'``, via a **fallback** on the :math:`O(3)` tag
     - ``'angular'``, from the **manifold alone**
   * - nodes / weights
     - —
     - **bit-identical** to ``gauss_legendre_on_mu(8)``
       (``np.array_equal`` on both)

⭐ **This is a repair, not wiring, and the measurement says so.** Take
the eight nodes and weights of ``Quadrature.gauss_legendre(8)`` and
build a *spatial* rule from them on ``Interval(-1, 1)``. Before the
declaration the two induced function spaces were `[M]` ``==`` **and
hash-equal** — an 8-node slab angular space and an 8-node spatial rule
were the same value, so any cache, ``set`` or operator-domain check
keyed on the space would have conflated them. After it, `[M]` ``==`` is
``False`` and the hashes differ. That is the energy/spatial collision of
tracker 2.1 recurring one level up, and 2.0c could not close it: while
both supports were honestly ``Interval(-1, 1)``, there was nothing to
tell apart.

.. _manifold-on-orbit-space:

``on_orbit_space`` — the third verb, and why it is neither of the others
------------------------------------------------------------------------

The declaration is performed by a new measure verb,
:meth:`DiscreteMeasure.on_orbit_space
<orpheus.numerics.measure.DiscreteMeasure.on_orbit_space>`. Its
semantics are equation-free, which is the whole of its content:
**the same atoms, re-read as chart coordinates of an orbit space.**
Same nodes, same weights — `[M]` the *same array objects*, not copies —
and only what the measure KNOWS about its support changes, from "an
interval" to "the polar marginal of a sphere, about this axis".

It is easy to mistake for the two verbs the corpus already has, and it
is neither:

.. list-table::
   :header-rows: 1
   :widths: 20 26 26 28

   * -
     - ``pushforward(φ)``
     - ``quotient(G)``
     - ``on_orbit_space(M/H)``
   * - Starts from
     - a measure on :math:`\mathcal{X}`
     - a measure on the **base** :math:`M`
     - a measure on the **chart** :math:`C`
   * - Does it move a node?
     - yes, applies :math:`\varphi`
     - no — selects orbit representatives
     - **no — applies nothing at all**
   * - Node count
     - unchanged
     - one per orbit (drops)
     - unchanged
   * - Mass
     - preserved
     - preserved (orbit-stabilizer weights)
     - **untouched**
   * - What changes
     - the points
     - the points and the support
     - **only the support**

⟹ a :math:`\mu`-rule was never on :math:`S^2` to begin with, so there
is no fold to perform; and no map is applied, so there is no
pushforward. The verb exists because "this list of numbers is a
coordinate list *for* :math:`S^2/SO(2)_x`" is a fact about the measure's
type that no arithmetic can supply.

**The one precondition, refused where the declaration is written.**
``on_orbit_space`` raises unless the orbit space's ``realization`` **is**
this measure's current support — the chart it was built on. `[M]`
handing a :math:`\mu`-rule the mirror quotient
``SPHERE.quotient(SubgroupOfO3.Mirror('y'))``, whose chart is the disk
:math:`D^2`, raises verbatim:

.. code-block:: text

   a measure on '[-1,1]' cannot be read on 'S^2/sigma_y': that orbit
   space's chart is 'D^2'. A rule declares the orbit space whose CHART
   it was built on; to fold a rule on the base manifold, use quotient()

The message names the mismatch *and* the other verb, because the two
failure modes ("wrong orbit space" and "you meant to fold") are exactly
the two ways a caller gets here.

**What the metadata does, and why.**

* ``invariance_group`` is **DROPPED**. A subgroup of :math:`O(3)` is a
  claim about an *embedding*, and the orbit space fixes an embedding —
  its axis — that the chart did not. The adopter re-tags: `[M]`
  :func:`~orpheus.numerics.quadrature.gauss_legendre_on_polar_orbit`
  immediately stamps ``Mirror(axis)`` back, for the *named* axis, which
  is a strictly more specific claim than the chart's could be. This is
  the same discipline the metadata-propagation table on
  :doc:`/theory/foundations/discrete_measures` applies everywhere: the
  field becoming ``None`` is correct behaviour, and the caller who knows
  the residual re-establishes it.
* ``exactness`` **survives**. The reference measure is a measure on the
  chart and the chart is unchanged, so the claim is untouched: `[M]`
  ``exact to algebraic degree 15 against legendre`` before and after,
  on ``n = 8``.

.. _manifold-polar-orbit-rule:

Two objects on one interval: the chart rule and the orbit rule
---------------------------------------------------------------

The declaration could not simply be pushed into
:func:`~orpheus.numerics.quadrature.gauss_legendre_on_mu`, and the
reason is the two-poles fact of
:ref:`manifold-so2-axis-is-a-parameter` in its operational form. That
function serves **two** roles, and only one of them has an axis to name:

* it is the raw material of the slab's rule, whose :math:`\mu` is the
  cosine against :math:`x`;
* it is the **polar factor** of every product rule in
  :mod:`orpheus.numerics.quadrature.rules_product`, whose :math:`\mu` is
  the cosine against :math:`z`.

A factor that declared :math:`S^2/SO(2)_x` while sitting inside a
product about :math:`z` would be a false claim about the object it is
part of. So the tree carries two functions:

.. list-table::
   :header-rows: 1
   :widths: 34 33 33

   * -
     - ``gauss_legendre_on_mu(n)``
     - ``gauss_legendre_on_polar_orbit(n, axis)``
   * - Support
     - ``COSINE_INTERVAL`` — the chart, naming no axis
     - ``SPHERE.quotient(SO2(axis))``
   * - ``invariance_group``
     - ``Mirror('x')`` (the canonical :math:`(\mu,0,0)` embedding)
     - ``Mirror(axis)``
   * - Consumed by
     - the product rules' polar factor; the raw material below
     - :meth:`Quadrature.gauss_legendre
       <orpheus.numerics.quadrature.directional.Quadrature.gauss_legendre>`
       and the ``GaussLegendre1D`` registry spec
   * - In the registry?
     - **no**, deliberately
     - **yes**, as ``partial(..., axis="x")``

⭐ **The chart-level rule is deliberately NOT registered, and stage 0 is
what enforces it.** `[M]` on the slab's ``AngularSymmetry``:
``admits_domain`` is **False** for ``gauss_legendre_on_mu(8)`` (support
``[-1,1]``), **False** for a marginal declared about :math:`y` or
:math:`z`, and **True** only for the :math:`x`-declared rule. One fact
on both sides — the geometry names the group it spends, the rule names
the group its orbit space was quotiented by — so a rule about the wrong
pole is refused by the same comparison that refuses a sphere cubature.

`[M]` the live rejection text a slab selection now emits for the three
:math:`S^2` rules reads:

.. code-block:: text

   domain mismatch: geometry 'slab' discretises S^2/SO2_x, but the
   rule's nodes live on S^2

.. note::

   ⚠ **The** ``phase`` **fallback arm did NOT become unreachable, and a
   pre-flight predicted it would.** :attr:`DiscreteMeasure.phase
   <orpheus.numerics.measure.DiscreteMeasure.phase>` classifies a
   measure by the TYPE of its support manifold, with one fallback: a
   rule on a bare :class:`~orpheus.numerics.manifold.Interval` carrying
   an :math:`O(3)` invariance tag is angular. Tracker 2.4 was expected
   to close that arm by making the slab declare a sphere quotient. It
   did — for the slab. `[M]` over eight shipped rules
   (``gauss_legendre_on_mu``, ``gauss_legendre_on_polar_orbit``,
   ``Quadrature.gauss_legendre`` / ``product`` / ``folded_product`` /
   ``level_symmetric`` / ``lebedev``, ``periodic_trapezoid``) the arm is
   reached by exactly **one**: ``gauss_legendre_on_mu`` itself, which for
   the reason above must keep a bare interval. Every other rule now
   answers from its manifold's type. The honest statement is that the
   fallback is **scoped to the chart-level rule**, not retired — and it
   is a live example of the plan hazard where an "unreachable after this
   step" prediction is falsified by the same step's own design
   constraint.

.. _manifold-the-axis-convention-for-a-section:

What a section will have to choose, in the axis's language
------------------------------------------------------------

:ref:`manifold-err-080-is-a-section` establishes that ERR-080's level-1
half is a fabricated **section** of :math:`S^2 \to S^2/SO(2)_x`, and that
``fundamental_domain=None`` is the honest entry because no section is
canonical. Tracker 2.4 does not change that — it changes the *language*
the eventual choice has to be made in, and that is worth writing down
before someone makes it.

Now that the axis is a parameter, the :math:`\varphi = 0` half-meridian
has an axis-general spelling. Writing :math:`a` for the rotation axis
and :math:`b, c` for the other two, the candidate section is

.. math::

   \mu \;\longmapsto\;
   \mu\,\hat e_a \;+\; \sqrt{1-\mu^2}\,\hat e_b \;+\; 0\cdot\hat e_c ,

which is on :math:`S^2` by construction. `[M]` for :math:`a = x` this is
:math:`\mu \mapsto (\mu, \sqrt{1-\mu^2}, 0)` and it is the object the
tree fabricated as :math:`(\mu, 0, 0)` — the fabrication is the same map
with the :math:`\hat e_b` term dropped, which is exactly why its rows
have norms :math:`|\mu| < 1` rather than :math:`1`.

⚠ **It is still a CHOICE.** Every half-meridian is a section;
:math:`\varphi = 0` merely names one, and which one you pick is a
convention about where azimuth zero sits, not a derivation. The shipped
catalogue entry therefore keeps ``fundamental_domain=None`` on purpose,
and the choice belongs to whichever step declares a section — not to
the orbit-space derivation, which would then be smuggling a convention
into a theorem.

.. note::

   ⛔ **This paragraph read "and it belongs to tracker 2.3 … the step
   that mints the typed** ``Chart`` **— not to the orbit-space
   derivation" until 2026-09-02.** Tracker 2.3
   landed on that date and the prediction is half right, in the way
   worth recording (``coding-standards``, the *falsified prediction*
   tense class): the **phase** was right and a typed map did land, but
   it is :class:`~orpheus.numerics.manifold.ManifoldMap` rather than a
   ``Chart``, and **the choice was not made**. `[M]` 2026-09-02 the
   three arrows 2.3 types are a parametrisation of the *base*
   (``archimedes``), a per-measure orbit *retraction*, and the orbit
   *barycentre*, which lands off :math:`S^2` by construction — none of
   them is a section, and
   :attr:`Quotient.fundamental_domain
   <orpheus.numerics.manifold.Quotient.fundamental_domain>` is still
   read by nothing outside :mod:`orpheus.numerics.manifold` itself.
   The naming ruling is the reason: *a chart is* :math:`M \supset U \to
   \mathbb{R}^n`, and only the **inverse** of the Archimedes map is
   one, so a type called ``Chart`` would have mis-described two of its
   own three instances (:ref:`manifold-arrows`). The section is still
   owed, and it is still a choice.


.. _manifold-arrows:

Maps between manifolds: the ARROWS
==================================

Everything above this point is about **objects**. A category needs
arrows too, and the tree had been drawing three of them freehand:
whenever a construction wanted to move a point set somewhere else it
applied a callable and then *named the destination by hand*, at the
call site, in whatever vocabulary was locally convenient. That is the
shape a forged claim takes — the same shape :ref:`ERR-080
<manifold-err-080>` has — because a destination named at the call site
is a claim nobody else made and nothing can contradict.

Tracker 2.3 (2026-09-02) gives the arrows a type. It adds no
mathematics: every one of the three maps below was already being
computed, correctly, in production. What it adds is that the
**codomain travels with the map** instead of being supplied by the
caller, so *"apply this map and declare the result to live on*
:math:`S^2`\ *"* stops being a sentence anyone can write.


.. _manifold-arrow-type:

The type: a map carries its two point sets
-------------------------------------------

:class:`~orpheus.numerics.manifold.ManifoldMap` is a frozen value with
three fields — ``domain``, ``codomain``, ``apply`` — and it is the
point-level analogue of a
:class:`~orpheus.numerics.operator.LinearOperator`: where an operator
carries the two *function spaces* it maps between, a map carries the
two *point sets*. The design ruling (user, 2026-09-02) was for **one**
value type with named maps as factories — ``archimedes(axis)``,
``barycentre(orbit_space)`` — exactly as
:data:`~orpheus.numerics.manifold.SPHERE` and ``LEGENDRE`` are values
of their own types rather than subclasses of them.

Two arrows are *induced* by every :math:`\varphi : M \to N`, and they
are the reason the type exists rather than being decoration:

- the **pushforward of a measure**,
  :meth:`DiscreteMeasure.pushforward
  <orpheus.numerics.measure.DiscreteMeasure.pushforward>` — the image
  measure :math:`\varphi_*\mu` lives on :math:`N` *because*
  :math:`\varphi` says so (:eq:`discrete-measure-pushforward`);
- the **pullback of a function**, :math:`f \mapsto f \circ \varphi`,
  which is what the change-of-variables identity evaluates on the
  pushed measure. ⛔ Nothing in the tree applies a pullback through
  this type yet: the planned consumer is the map that restricts a
  basis to an orbit space (tracker 3.4b), and until it lands the
  second arrow is a statement about what the type *is*, not about
  what ships.

⭐ **The verb that changed is the pushforward, and its three states are
the campaign in miniature.** The same operation, spelled three ways
across five weeks:

.. list-table:: ``pushforward`` — who names the target
   :header-rows: 1
   :widths: 16 40 44

   * - When
     - Spelling
     - Who names :math:`\mathcal{Y}`, and what could contradict them
   * - until 2026-09-01
     - ``μ.pushforward(f)``
     - **Nobody.** The support was *fabricated* as
       ``f"φ_*({self.support})"`` — a manifold nobody has derived,
       wearing a name that makes it look like one
       (:ref:`manifold-string-algebra`).
   * - 2026-09-01 (2.0c)
     - ``μ.pushforward(f, new_space=Y)``
     - **The call site.** An improvement — only :math:`\varphi`'s
       author knows :math:`\mathcal{Y}`, and now somebody had to say
       it — but the caller is not always that author, and a caller who
       is wrong is unopposed.
   * - 2026-09-02 (2.3)
     - ``μ.pushforward(φ)``
     - **The map.** ``new_space=`` is retired; the target is *read*
       off ``φ.codomain``, and the verb additionally **refuses** a map
       whose ``domain`` is not this measure's support.

The refusal is by manifold **value**, not by array shape, and the
sharpest witness for that is a pair of measures carrying *literally the
same numbers*. `[M]` 2026-09-02:
``gauss_legendre_on_polar_orbit(4, "x").nodes`` is ``np.array_equal``
to ``gauss_legendre_on_mu(4).nodes`` — the slab's rule is the chart
rule with a declared orbit space (:ref:`manifold-polar-orbit-rule`) —
and handing the first to the product embedding raises where the second
is accepted:

.. code-block:: text

   ValueError: cannot push a measure on 'S^2/SO2_x × S^1' forward
   along a map out of '[-1,1] × S^1': the map's domain must be the
   measure's support. Build the map out of 'S^2/SO2_x × S^1', or hand
   this verb a measure on '[-1,1] × S^1' — the same numbers on a
   different manifold are a different measure.

Note where the difference surfaces: not on the polar factor but on the
**product**, because :meth:`Manifold.__mul__
<orpheus.numerics.manifold.Manifold.__mul__>` carried it there. The
algebra of objects is what makes the arrow's guard discriminating.

.. note::

   **No membership check runs inside the map**, deliberately. A map
   whose ``apply`` lands outside its declared codomain is a real
   defect, and its ruled home is
   :meth:`~orpheus.numerics.manifold.Manifold.contains` at *measure
   construction* (tracker 2.0b) — one refusal, on the object that
   actually escapes, rather than two half-refusals. ⛔ That is still
   **not built**: nothing calls ``contains`` on the way in, so a
   measure whose nodes are not points of its support remains
   constructible and **ERR-080 remains open**
   (:ref:`manifold-seams`).

   The **Jacobian** is likewise not this type's business. ``pushforward``
   is the :math:`\varphi`-image with weights preserved verbatim; a
   change of variables against a target *reference* measure is the
   caller's, and that asymmetry is documented where the identity lives
   (:doc:`/theory/foundations/discrete_measures`).


.. _manifold-three-arrows:

The three maps the tree was already spelling
---------------------------------------------

`[M]` 2026-09-02, at 2.3's opener the tree drew three arrows around the
quotient and none of them was typed. One of them it drew **twice** —
once honestly and once as ERR-080.

.. list-table:: Three arrows, four spellings
   :header-rows: 1
   :widths: 22 30 24 24

   * - Arrow
     - How it was spelled
     - Codomain, before
     - Codomain, after
   * - the orbit retraction
       :math:`M \to M/H`
     - a ``lambda`` plus a hand-written ``new_space=``, inside
       :meth:`DiscreteMeasure.quotient
       <orpheus.numerics.measure.DiscreteMeasure.quotient>`
     - named at the call site
     - ``self.support.quotient(group)`` — the **catalogue's own
       object**, by identity
   * - the Archimedes chart
       :math:`[-1,1]\times S^1 \to S^2`
     - a hand-written double loop plus the literal
       ``support=SPHERE``, inside ``spherical_product``
     - a literal
     - :data:`~orpheus.numerics.manifold.SPHERE`, **read** off
       ``archimedes("z").codomain``
   * - the orbit barycentre
       :math:`S^2/SO(2)_a \to D^3`, **honest** spelling
     - inline in ``symmetry._embedded_nodes``, which embeds a polar
       marginal for an invariance check
     - (untyped — a bare array)
     - :class:`~orpheus.numerics.manifold.Ball`\ ``(3)``, via
       :func:`~orpheus.numerics.manifold.barycentre`
   * - the **same** map, **forged** spelling
     - inline in ``Quadrature._harmonic_frame_measure``'s 1-D arm
     - ⛔ ``support=SPHERE`` — **a lie**
     - ⛔ unchanged, **by design** (see below)

⭐ **The fourth row is the whole argument.** Rows three and four compute
*the same image* — `[M]` 2026-09-02 on ``gauss_legendre(8)``, the
forgery's nodes are ``np.array_equal`` to
``barycentre(measure.support)(measure.nodes)`` — and differ only in
what they claim about where that image lives. The honest one says
:math:`D^3` and is right (`[M]` ``Ball(3).contains`` → ``True``); the
forged one says :math:`S^2` and is wrong (`[M]` ``Sphere().contains``
→ ``False``, norms :math:`0.1834 \ldots 0.9603`). A codomain that is a
**field of the map** cannot be forged at the call site; that is the
entire purchase of this type.

The forgery therefore **stays a raw constructor** until tracker 3.4
retires it, and that is not an oversight: it *cannot* be re-spelled
through ``pushforward`` without telling the truth about its codomain,
so re-spelling it would silently repair ERR-080's level-1 half in a
step whose subject is the type system. The arm carries a comment
naming the map it is a forgery of.


.. _manifold-archimedes:

Archimedes: :math:`[-1,1] \times S^1 \to S^2`
----------------------------------------------

Write :math:`a` for an axis and :math:`b, c` for its two cyclic
successors (:math:`z \to x, y`; :math:`x \to y, z`; :math:`y \to z, x`
— a right-handed frame in every case). The map is

.. math::

   \varphi_a(\mu, (\cos\varphi, \sin\varphi))
   \;=\; \mu\,\hat e_a
   \;+\; \sqrt{1-\mu^2}\,\bigl(\cos\varphi\,\hat e_b
                               + \sin\varphi\,\hat e_c\bigr),

which for :math:`a = z` is the direction-cosine triple every product
rule has always used —
:math:`\mu_x = \sin\theta\cos\varphi`,
:math:`\mu_y = \sin\theta\sin\varphi`,
:math:`\mu_z = \mu`, with :math:`\sin\theta = \sqrt{1-\mu^2}` — stated
as an equation in the module docstring of
``orpheus.numerics.quadrature.rules_product`` and pinned against the
map by hand in
``tests/numerics/test_manifold.py::test_archimedes_about_z_is_the_labelled_equation_verbatim``.

**Its relation to the orbit-space chart is exact, and it is a
projection.** The chart of :math:`S^2/SO(2)_a` is
:math:`\pi(\Omega) = \Omega\cdot\hat e_a`
(:ref:`manifold-s2-so2`, step 5), and composing gives

.. math::

   \pi \circ \varphi_a \;=\; \mathrm{pr}_1 ,

`[M]` bit-exactly: over 500 random :math:`(\mu, \varphi)` pairs per
axis, :math:`\max\lvert \pi(\varphi_a(\mu,\varphi)) - \mu \rvert =`
**0.0** for :math:`a \in \{x,y,z\}`, with
:math:`\max\bigl\lvert\lVert\varphi_a\rVert - 1\bigr\rvert \le`
**2.22e-16**. So the polar factor of a product rule *is* the orbit-space
coordinate, and the circle factor is the fibre the quotient forgets —
which is why one Gauss–Legendre rule can serve both a slab marginal and
a product rule's polar factor
(:ref:`manifold-so2-axis-is-a-parameter`).

.. admonition:: Why it is named for Archimedes
   :class: note

   **Archimedes' hat-box theorem** is the statement about this map that
   transport actually needs: the pushforward of the uniform measure
   :math:`d\Omega` along :math:`\mu = \Omega\cdot\hat e_a` is
   :math:`2\pi\,d\mu` — *uniform* on :math:`[-1,1]`, for every axis. A
   Gauss–Legendre rule in :math:`\mu` times any circle rule is therefore
   exact on the sphere against Lebesgue measure, which is the theorem
   :func:`~orpheus.numerics.quadrature.rules_product.spherical_product_claim`
   composes the two factors' claims through.

   ⚠ **Edited elsewhere, consumed here.** The corpus already owns that
   fact in the *registry* register — which reference measure the
   quadrature registry answers for a spent group, and why it keys on
   ``rotation_axis`` being non-``None`` rather than on one group
   (:doc:`/theory/foundations/discrete_measures`). This page owns it in
   the *map* register: the hat-box is a statement about
   :math:`\varphi_a`, and the pushforward reference measure it names is
   a field of the **catalogue entry**, not of the map
   (:ref:`manifold-arrows-not-built`).

⚠ **It is a parametrisation, not a chart in the strict sense**, and the
place it fails is not incidental — it is the stratum. The circle factor
**collapses** at :math:`\mu = \pm 1`: a whole fibre maps to one pole.
That is exactly the singular locus of :math:`S^2/SO(2)_a`
(:ref:`manifold-singular-stratum`), so the map is injective off the
stratum and the stratum is precisely where it is not. Its inverse on
:math:`S^2 \setminus \{\pm\hat e_a\}` is the :math:`(\mu, \varphi)`
chart. The collapse is gated directly: `[M]` on a
:math:`7 \times 8` grid the eight images at :math:`\mu = +1` are
``np.array_equal`` to :math:`\hat e_a` repeated, and likewise at
:math:`\mu = -1`.

**What the product rule now reads.**
:func:`~orpheus.numerics.quadrature.rules_product.spherical_product`
was a hand-written double loop that filled a
:math:`(n_\mu n_\varphi, 3)` array and then declared
``support=SPHERE``. It is now the algebra it always was:

.. code-block:: python

   measure = (polar * azimuthal).pushforward(archimedes("z")).with_metadata(
       invariance_group=group,   # DERIVED from the factors' generators
       exactness=claim,          # DERIVED through the product theorem
   )

— the tensor product is the measure's own
:meth:`~orpheus.numerics.measure.DiscreteMeasure.__mul__` on the
product manifold :math:`[-1,1]\times S^1`, the embedding is the typed
chart, and the support is the chart's codomain. `[M]` 2026-09-02,
re-run independently of the gate: over **60** configurations
(:math:`n_\mu \in \{2,3,4,5,6\}` :math:`\times` :math:`n_\varphi \in
\{6,8,10,16,24,32\}` :math:`\times` both circle shifts), **0** differ
from the transcribed pre-2.3 loop — ``np.array_equal`` on nodes and on
weights — and ``measure.support`` **is** ``archimedes("z").codomain``,
one object rather than a literal that happens to agree.

⭐ Two derivations of the same rule now meet: the *claim* side composes
the factors' exactness through the product theorem, and the *measure*
side composes their atoms through this chart. Neither reads the other,
and both refuse a mismatched factor pair — the claim side on the
exactness system, the measure side on the manifold.


.. _manifold-barycentre:

The orbit barycentre: :math:`S^2/SO(2)_a \to D^3`, and why it is not a section
-------------------------------------------------------------------------------

An orbit of :math:`SO(2)_a` is the circle
:math:`\{\Omega : \Omega\cdot\hat e_a = \mu\}`, of radius
:math:`\sqrt{1-\mu^2}` about the point :math:`\mu\,\hat e_a` on the
axis. That point is the orbit's **barycentre** — its mean under the
fibre's uniform measure — and the map

.. math::

   \beta_a : \mu \;\longmapsto\; \mu\,\hat e_a

is what a consumer wanting *one representative point per orbit* keeps
reaching for. Where it lands is a one-line computation, and it is the
whole story:

.. math::

   1 - \lVert \mu\,\hat e_a \rVert^2 \;=\; 1 - \mu^2
   \;=\; \tfrac14 \det P ,

the **squared orbit radius**, which is the quantity the catalogue
already records as :attr:`Quotient.det_gram
<orpheus.numerics.manifold.Quotient.det_gram>` — `[M]` the shipped
entry's ``det_gram`` is :math:`4 p_2` (:eq:`manifold-s2-mod-so2`),
which restricted to the sphere by :math:`p_1^2 + p_2 = 1` is
:math:`4(1-\mu^2)` (:ref:`manifold-s2-so2`, step 5); the identity
reproduces to **0.0** over nine :math:`\mu` values. So :math:`\beta_a` lands **on**
:math:`S^2` exactly where the orbit radius vanishes — the two poles,
i.e. the singular stratum — and strictly **inside** the ball
everywhere else. Its codomain is
:class:`~orpheus.numerics.manifold.Ball`\ ``(3)`` and can be nothing
else: `[M]` ``Ball(3).contains`` → ``True`` on the whole image,
``Sphere().contains`` → ``False`` on the interior and ``True`` on the
two poles.

⛔ **It is not a section, and the distinction is the one ERR-080 got
wrong.** A section of :math:`S^2 \to S^2/SO(2)_a` lands **on**
:math:`S^2` by picking a representative *direction*; for a
positive-dimensional group no such pick is canonical, which is why
every :math:`S^2/SO(2)_a` entry carries ``fundamental_domain=None`` on
purpose (:ref:`manifold-err-080-is-a-section`). The barycentre is
canonical *precisely because it is not a representative*: it is the
orbit's mean, and a mean of unit vectors is not a unit vector.

⟹ **ERR-080, restated in one sentence of this section's vocabulary:
the forgery is** :math:`\beta_a` **with its codomain declared as**
:math:`S^2` **instead of** :math:`D^3`. Everything else about the
computation is right. That is why the defect is invisible to any
arithmetic check and why no tolerance reaches it: the numbers are the
correct barycentres, and what is false is a *type*.

.. list-table:: The two spellings of :math:`\beta_a`, measured on ``gauss_legendre(8)``
   :header-rows: 1
   :widths: 30 35 35

   * -
     - ``symmetry._embedded_nodes``
     - ``Quadrature._harmonic_frame_measure`` (1-D arm)
   * - what it computes
     - :math:`\mu \mapsto \mu\,\hat e_a`, axis read off the support's
       group
     - :math:`\mu \mapsto (\mu, 0, 0)` by zero-padding
   * - `[M]` are the images equal?
     - **yes** — ``np.array_equal``, both directions
     - **yes**
   * - declared codomain
     - ``Ball(3)`` — via
       :func:`~orpheus.numerics.manifold.barycentre` since 2.3
     - ⛔ ``SPHERE``
   * - `[M]` is the declaration true?
     - ✅ ``Ball(3).contains`` → ``True``
     - ⛔ ``Sphere().contains`` → ``False``; norms
       :math:`0.1834\ldots0.9603`
   * - why it wants the barycentre
     - an invariance check *should* use it: a rotation about :math:`a`
       **fixes** the barycentre, so the point is genuinely
       :math:`SO(2)_a`-invariant
     - it wants a **direction**, which the barycentre is not — so it
       is the wrong map, honestly applied

⭐ **The honest spelling now reads the map** (Pattern 2 — one spelling
of one concept). ``symmetry._embedded_nodes``'s polar-marginal arm is
``barycentre(measure.support)(nodes)``, so the two cannot drift; `[M]`
bit-identical on **12** rows (``gauss_legendre_on_polar_orbit(n, axis)``
for :math:`n \in \{2,4,8,16\}` and all three axes). And the map
**refuses** anything that is not an axial orbit space — a mirror
quotient, the trivial quotient, a bare interval — because a mirror
orbit is a pair of points with no axis to lie on, and a point of an
interval is not an orbit at all.

.. note::

   ⚠ **A brief for this step reported that** :class:`Ball`
   **had zero production consumers before it.** `[M]` that is not what
   the tree says, and the correct statement is sharper. ``Ball(2)`` has
   been production since 2026-08-31 as the :math:`\sigma_y` entry's
   ``realization`` (``manifold.py``, the mirror derivation) and is
   matched in the ambient-dimension table. What had **never** been
   constructed anywhere — `[M]` ``git grep "Ball("`` over ``orpheus/``
   and ``tests/`` at the pre-2.3 commit returns **six** lines — the
   class definition, one ``match`` pattern, and **four** constructions,
   every one of them ``Ball(2)`` — is ``Ball(3)``, and what is new
   in kind is that a :class:`Ball` is now the **codomain of an arrow**
   rather than a field of a catalogue entry. The production docstring
   states the weaker claim; it is reported rather than edited here.


.. _manifold-arrow-composition:

Composition, functoriality, and the fold as a two-arrow chain
---------------------------------------------------------------

Composition is ``psi @ phi`` for :math:`\psi \circ \varphi`, refused
unless ``phi.codomain == psi.domain`` — the same guard, at the same
tier, as an operator product over unequal spaces:

.. code-block:: text

   ValueError: cannot compose: the inner map lands on 'S^2' but the
   outer map is defined on '[-1,1] × S^1'. A composition psi @ phi
   needs phi.codomain == psi.domain.

The one law worth stating — and, for a map of finite point sets whose
membership :meth:`~orpheus.numerics.manifold.Manifold.contains`
already governs, essentially the only intrinsic thing such a map can
get wrong — is that the pushforward is a **functor**:

.. math::
   :label: manifold-map-functoriality

   (\psi \circ \varphi)_* \mu \;=\; \psi_*\bigl(\varphi_* \mu\bigr).

.. (vv-status rationale) manifold-map-functoriality: A structural law
   of the arrow type, not a solver claim, so it carries no L0..L3
   ladder slot and no ``verifies(...)`` marker. It is gated by the
   ``foundation`` test
   ``tests/numerics/test_manifold.py::TestManifoldMap::test_functoriality_the_pushforward_of_a_composite_is_the_composite_of_pushforwards``
   (nodes and weights by ``np.array_equal``, and the support by
   value), and re-measured on production objects by the shipped-fold
   table in this section.
.. vv-status: manifold-map-functoriality documented

⭐ **This is not an abstract law here — the shipped cylindrical fold is
literally a chain of two of these arrows.**
:meth:`Quadrature.folded_product
<orpheus.numerics.quadrature.directional.Quadrature.folded_product>`
builds a product rule and then folds it, and after 2.3 both halves are
typed maps:

.. math::

   [-1,1] \times S^1
   \;\xrightarrow{\ \varphi_z\ }\; S^2
   \;\xrightarrow{\ \rho\ }\; S^2/\langle\sigma_y\rangle ,

with :math:`\rho` the orbit retraction that
:meth:`DiscreteMeasure.quotient
<orpheus.numerics.measure.DiscreteMeasure.quotient>` builds from the
invariance certificate. Because :math:`\varphi_z` lands on
:math:`S^2` and :math:`\rho` is defined there, ``rho @ chart``
type-checks, and :eq:`manifold-map-functoriality` says the one-shot
route must agree with the two-step one.

`[M]` 2026-09-02, measured — the composite built with ``@``, pushed in
one step and consolidated, against the two-step route, against the
shipped rule:

.. list-table:: The fold, three ways
   :header-rows: 1
   :widths: 18 14 24 22 22

   * - :math:`(n_\mu, n_\varphi)`
     - :math:`N` after the fold
     - one-shot ``ρ @ φ`` vs two-step
     - vs shipped ``folded_product``
     - ``support``
   * - (2, 8)
     - 8
     - ``array_equal`` ✅
     - ``array_equal`` ✅
     - ``'S^2/sigma_y'``, by **identity** with the catalogue entry
   * - (4, 8)
     - 16
     - ✅
     - ✅
     - ✅
   * - (4, 16)
     - 32
     - ✅
     - ✅
     - ✅
   * - (6, 10)
     - 30
     - ✅
     - ✅
     - ✅
   * - (3, 24)
     - 36
     - ✅
     - ✅
     - ✅

⚠ **Read the fixture, not just the ticks.** The reconstruction uses the
**staggered** circle rule, which is what ``folded_product`` selects
(:math:`\Sigma = \varnothing`, the fold's well-posedness condition);
with the node-aligned shift the same code agrees with itself on both
routes but does *not* reproduce the shipped rule, because it is then a
different rule. `[M]` at :math:`(2, 8)`: node-aligned puts **4** nodes
on :math:`\Sigma = \{\xi = 0\}`, so its 16 atoms fold into **10**
orbits — sizes ``[1,1,1,1,2,2,2,2,2,2]``, the four singletons being
the fixed points — while staggered puts **0** there and folds into
**8**, all of size 2. The functoriality half is fixture-independent;
the third column is a statement about which circle rule ships, and the
singleton orbits are exactly the well-posedness condition the fold's
own :math:`\Sigma = \varnothing` requirement names.


.. _manifold-arrows-not-built:

What 2.3 did NOT build
------------------------

Three things, stated so the next phase does not re-derive a decision
already taken and so no reader mistakes an arrow for a repair.

**(1) Neither of the catalogue entry's own two maps ships.** The
:ref:`engine data model <manifold-engine-data-model>` lists the
entry's *chart* :math:`M/H \to N` and its *section* as procedure
outputs that are not slots, and 2.3 does not change that — it changes
only whether such a thing could be *expressed*. None of the three
arrows above is either of them:

.. list-table::
   :header-rows: 1
   :widths: 30 32 38

   * - Arrow
     - Type
     - Why it is not the entry's chart or section
   * - :math:`\varphi_a`, ``archimedes``
     - :math:`[-1,1]\times S^1 \to S^2`
     - a parametrisation of the **base**, not a chart of the orbit
       space. Its first component *is* the chart, in the sense
       :math:`\pi \circ \varphi_a = \mathrm{pr}_1` measured above, but
       the chart :math:`\pi : S^2 \to [-1,1]` itself is still not a
       value anywhere.
   * - :math:`\rho`, the retraction
     - :math:`M \to M/H`
     - built **per measure** from an invariance certificate, so it
       depends on the atoms and not only on :math:`(M, H)` — it cannot
       be a field of an entry. Its image stays in the base's
       coordinates, which is the section's coordinate system, but it
       is not a section: it is the quotient map with a chosen
       representative per *realized* orbit.
   * - :math:`\beta_a`, ``barycentre``
     - :math:`S^2/SO(2)_a \to D^3`
     - a map **out of** the orbit space, landing off :math:`S^2`. A
       section is a map **into** the base. See above.

⟹ ``fundamental_domain=None`` on every :math:`S^2/SO(2)_a` entry is
still the honest answer, and :attr:`Quotient.fundamental_domain
<orpheus.numerics.manifold.Quotient.fundamental_domain>` still has
`[M]` **zero production readers**. The section remains a *choice*
(:ref:`manifold-the-axis-convention-for-a-section`), and 2.3 declined
to make it.

**(2) The pushforward reference measure is deferred to tracker 3.1 —
and the reason is a measured import cycle.** An orbit space's
pushforward reference (the :math:`2\pi\,d\mu` of the hat-box) is a
field of the catalogue **entry**, not of the map, because it is a
property of :math:`(M, H)` rather than of any one arrow. It is
answered today by a twin on the registry —
:attr:`AngularSymmetry.reference
<orpheus.numerics.quadrature.registry.AngularSymmetry.reference>`,
``LEGENDRE`` for any axial rotation and ``UNIFORM_ON_SPHERE`` for the
trivial group — and collapsing that twin onto the entry is the same
move tracker 2.4 made for ``support`` (:ref:`manifold-twin-lookup`).

⛔ It cannot be done by adding an import. `[M]` 2026-09-02 by AST with
relative imports resolved, :mod:`orpheus.numerics.exactness` imports
:mod:`orpheus.numerics.manifold` at **module scope, twice** — once for
:class:`~orpheus.numerics.manifold.Manifold` and once for
:data:`~orpheus.numerics.manifold.CIRCLE` and
:data:`~orpheus.numerics.manifold.SPHERE`. A module-scope
``manifold → exactness`` edge therefore closes a **two-hop** cycle,
with no import order that survives it: demonstrated on a throwaway
package carrying exactly that topology — no production file touched —
**5 of 5** entry points die with

.. code-block:: text

   ImportError: cannot import name 'Manifold' from partially
   initialized module 'pkg.manifold' (most likely due to a circular
   import)

and **5 of 5** import cleanly with the ``TYPE_CHECKING`` guard
restored (the positive control, without which a clean reading carries
no information). ⟹ the viable mechanism is a
:class:`~orpheus.numerics.manifold.Quotient` field populated **inside**
the derivation function through a function-scope import — the idiom
``_sphere_mod_so2`` already uses for SymPy — never at module scope.
This is the same guard as :ref:`manifold-import-cycle`, on a second
pair of modules.

**(3) ERR-080 is not repaired, and 2.3 moves neither of its gates.**
Nothing here calls
:meth:`~orpheus.numerics.manifold.Manifold.contains` on the way into a
measure, so the forged :math:`(\mu, 0, 0)` measure is still
constructible; the forgery arm is still a raw
:class:`~orpheus.numerics.measure.DiscreteMeasure` constructor by
design; and the level-2 half — the trivial isotypic sub-basis
:math:`\{Y_\ell^0\} \cong \{P_\ell\}` — is untouched. `[M]` by AST
over ``tests/sn/solve/test_pl_order_does_not_move_the_infinite_medium_flux.py``
the module still declares **three** ``@pytest.mark.xfail(strict=True)``
rows and 2.3 edits none of them. What 2.3 buys ERR-080 is a
*sentence*: the defect now has a name in the type system's own
vocabulary — :math:`\beta_a` with a forged codomain — and one honest
implementation of that map to point at.


.. _manifold-basis-invariance-group:

The second operand: a basis declares the symmetry its functions HAVE
====================================================================

The section above gave the pairing's **measure** side: a rule that says
which orbit space its atoms live on. This section is the **basis** side
— tracker 2.1b, 2026-09-01 — and the whole of its content is that the
basis side needed *no new field*. It follows from one elementary fact
about functions on a quotient:

.. math::

   \mathcal{F}(M/H) \;\;\xrightarrow[\;\cong\;]{\;\;\pi^*\;\;}\;\;
   \bigl\{\, f \in \mathcal{F}(M) \;:\; f \circ h = f \ \ \forall h \in H \,\bigr\}
   \;=\; \mathcal{F}(M)^H ,

with :math:`\pi : M \to M/H` the orbit projection and
:math:`\mathcal{F}` any function class the projection respects —
continuous, measurable, :math:`L^2`, smooth off the singular stratum.
Nothing here needs regularity; the statement is set-theoretic. Pulling a
function back along :math:`\pi` produces an :math:`H`-invariant function
on :math:`M`, and every :math:`H`-invariant function descends to the
quotient, because it is constant on orbits and the orbits *are* the
points of :math:`M/H`. Being a function on :math:`M/H` **is** being
:math:`H`-invariant, spelled two ways. So a basis that has
already named :math:`M/H` as its :attr:`domain
<orpheus.numerics.basis.base.Basis.domain>` has already declared its
group: the group is :attr:`Quotient.by
<orpheus.numerics.manifold.Quotient.by>`, sitting inside the slot
tracker 2.1 minted.

⭐ **The tracker asked for the wrong object, and the phase opener said
so.** Its row read, verbatim, *"``Basis.invariance_group`` — absent;
derivable for every shipped basis"*, and the plan's own census measured
*"0 of 6 subclasses answer it"*. Both true — and the design they invited
was a second abstract property with six overrides, kept in agreement
with ``domain`` by hand: exactly the two-tags-that-drift shape this page
exists to argue against. `[M]` **0 of 6** shipped
bases carried the name before this step — ``git show HEAD:`` over every
module of :mod:`orpheus.numerics.basis` and over
:mod:`orpheus.sn.operators.loss_kernel_gauge`, which is where
:class:`~orpheus.sn.operators.loss_kernel_gauge.LossKernelBasis` lives,
returns zero occurrences in each. After it, **6 of 6** answer (the
denominator is ``Basis.__subclasses__()`` walked recursively at runtime,
not a hand-list) and the basis-side diff is **one** file,
``orpheus/numerics/basis/base.py``: one concrete ``@final`` property on
the ABC, **zero** subclass edits. That is ``coding-standards``' *clean
before extending* landing as a no-op extension through a single generic
body, and it is the same dissolution tracker 2.0d's ``quotient_group``
**field** underwent at 2.0c, one level over: the fact was already in the
type, so the field would have been its second home.

.. _manifold-has-versus-spent:

HAS and SPENT: one slot for a function, two for a point set
------------------------------------------------------------

A measure carries **two** group slots and a basis carries **one**, and
that asymmetry is not an oversight on either side — it is the difference
between a point set and a function.

* :attr:`DiscreteMeasure.invariance_group
  <orpheus.numerics.measure.DiscreteMeasure.invariance_group>` — what the
  atom list **HAS**: a stored field recording a subgroup under which the
  nodes, weights included, are closed. It is a *declaration*, not a
  computed stabiliser, and ``None`` means unspecified rather than
  trivial.
* :attr:`DiscreteMeasure.quotient_group
  <orpheus.numerics.measure.DiscreteMeasure.quotient_group>` — what it
  **SPENT**: the group its support was folded by, derived from
  :attr:`Quotient.by <orpheus.numerics.manifold.Quotient.by>` and stored
  nowhere (tracker 2.0c).

For a POINT SET these come apart in both directions, and the shipped
rules realise three of the four combinations. The table is **exhaustive
over the family**, not a sample: its denominator is every
``classmethod`` factory on
:class:`~orpheus.numerics.quadrature.directional.Quadrature`,
enumerated by ``vars(Quadrature)`` — `[M]` **five of five**
(``vv-principles`` #31's finite-roster corollary: for an enumerable
shipped set, probe every member, because the one you skip is where the
counterexample lives).

.. list-table:: `[M]` 2026-09-01 — HAS and SPENT on all five shipped rules
   :header-rows: 1
   :widths: 26 20 27 27

   * - Rule
     - ``support.name``
     - HAS (``invariance_group``)
     - SPENT (``quotient_group``)
   * - ``lebedev(17)``
     - ``'S^2'``
     - ``OctahedralOh``
     - ``None``
   * - ``level_symmetric(8)``
     - ``'S^2'``
     - ``OctahedralOh``
     - ``None``
   * - ``product(4, 8)``
     - ``'S^2'``
     - ``Dnh(8)``
     - ``None``
   * - ``gauss_legendre(8)``
     - ``'S^2/SO2_x'``
     - ``Mirror('x')``
     - ``SO2('x')``
   * - ``folded_product(4, 8)``
     - ``'S^2/sigma_y'``
     - ``None``
     - ``Mirror('y')``

⚠ The missing combination is (**no** HAS, **no** SPENT) — an untagged
rule on the bare sphere. Nothing shipped is in that state, which is a
fact about the rules, not about the type.

Read the last two rows. The slab's polar rule HAS :math:`\sigma_x` and
SPENT :math:`SO(2)_x` — **two different groups in two slots on one
measure**, so no single field could carry both. And the
:math:`\sigma_y`-folded product rule HAS **nothing**, precisely because
it spent :math:`\sigma_y`: folding keeps one representative per orbit,
and a set with one point of each mirror pair is no longer closed under
the mirror. *Spending a symmetry destroys having it* — which is why
reading either slot as the other is ``plan-authoring`` §3's
ambiguous-name hazard, and why
:attr:`~orpheus.numerics.measure.DiscreteMeasure.quotient_group`'s own
docstring says so.

For FUNCTIONS the two collapse, by the isomorphism above. There is no
"folded away" to lose: a basis on :math:`M/H` *is* a set of
:math:`H`-invariant functions on :math:`M`, so what it has and what its
domain spent are one property. A basis therefore carries **one** slot,
named for what it HAS and read off what its domain SPENT:

.. list-table:: `[M]` 2026-09-01 — the one slot, over all six shipped bases
   :header-rows: 1
   :widths: 36 30 34

   * - Basis
     - ``domain.name``
     - ``invariance_group``
   * - ``SphericalHarmonicBasis(L)``, :math:`L \in \{0,1,3,7\}`
     - ``'S^2'``
     - ``Trivial``
   * - ``MirrorEvenSphericalHarmonicBasis(L=2, mirror_axis=a)``,
       :math:`a \in \{x,y,z\}`
     - ``'S^2/sigma_a'``
     - ``Mirror(a)`` — and ``is domain.by``
   * - ``IndicatorBasis`` from
       :meth:`EnergyGrid.as_basis
       <orpheus.data.energy_grid.EnergyGrid.as_basis>`
     - ``'energy'`` (an :class:`EnergyGroups`)
     - ``None``
   * - ``IndicatorBasis`` from
       :meth:`Mesh1D.indicator_basis
       <orpheus.geometry.mesh.Mesh1D.indicator_basis>`
     - ``'spatial_R1'`` (a :class:`RealSpace`)
     - ``None``
   * - ``WeightedIndicatorBasis``, ``OverlapBasis`` — both **delegate**
       ``domain`` to the basis they wrap
     - the wrapped basis's
     - ``None``, by delegation
   * - ``LossKernelBasis``
     - ``'index(sn_trace_orbit(...)_g)'`` (an :class:`IndexSet`)
     - ``None``

The mirror-even row is the one that carries the design: `[M]` the
answer is not merely *equal* to ``domain.by``, it **is** it —
``basis.invariance_group is basis.domain.by`` — so there is no second
object that could drift. A stored copy that happened to be right would
pass ``==`` and fail ``is``, which is what
``test_e2b_a_mirror_even_harmonic_HAS_its_mirror_read_off_its_domain``
asserts.

.. _manifold-invariance-three-arms:

Three arms, and why the answer off the sphere is ``None``, not ``Trivial``
---------------------------------------------------------------------------

The derivation is a ``match`` on the **type** of ``domain``, with three
arms:

.. list-table::
   :header-rows: 1
   :widths: 34 20 46

   * - ``domain``
     - Answer
     - Why
   * - ``Quotient(base=Sphere(), by=H)``
     - :math:`H`
     - the functions descend from :math:`M/H`, so they are exactly the
       :math:`H`-invariant ones
   * - ``Sphere()``
     - ``Trivial``
     - :math:`O(3)` **acts** on the domain and the basis has spent none
       of it — a domain of :math:`S^2` promises no invariance, whatever
       the individual functions happen to have (see the lower bound
       below)
   * - anything else
     - ``None``
     - no subgroup of :math:`O(3)` acts at all

⚠ **The third arm is a category answer, and** ``Trivial`` **would be a
lie.** No subgroup of :math:`O(3)` acts on a spatial mesh, an
energy-group index or a trace-DOF index set — there is no rotation of a
list of group boundaries. ``Trivial`` names the subgroup :math:`\{e\}`
**of** :math:`O(3)`, so writing it asserts that :math:`O(3)` acts on
this domain at all; ``None`` says the question does not arise. The
distinction is exactly the one
:attr:`DiscreteMeasure.phase <orpheus.numerics.measure.DiscreteMeasure.phase>`
already draws for the *same* manifolds when it refuses to classify a
non-angular support as angular, and it is why the two spellings of
"nothing" on the two sides mean opposite things: a full-sphere rule's
``quotient_group`` is ``None`` because it **spent nothing**, while
full-sphere harmonics' ``invariance_group`` is ``Trivial`` because they
**have** the trivial group. Same word in English, different lattice
elements — :math:`\{e\}` on one side, *no answer* on the other.

⭐ And the arms are decided by the domain and by nothing else, which is
a testable claim rather than a description:
``test_e4_the_group_is_decided_by_the_domain_and_by_nothing_else`` runs
all three on **one** class shape whose instances differ only in
``domain``, so an implementation keyed on the subclass — an
``isinstance`` on ``SphericalHarmonicBasis``, say — would give every
stub the same answer and fail.

.. _manifold-invariance-lower-bound:

The reading is a LOWER BOUND, and that is why the property is ``@final``
-------------------------------------------------------------------------

The domain gives the symmetry the basis is *guaranteed* to have, not the
largest one it happens to have. `[M]` ``SphericalHarmonicBasis(L=0)`` is
a single constant function, invariant under all of :math:`O(3)` — and it
answers ``Trivial``, at :math:`L \in \{0,1,3,7\}` alike, because its
domain says :math:`S^2` and a domain of :math:`S^2` promises nothing
more. The property is a **declaration read off a type**, never a
computed stabiliser.

Under-declaring is therefore *legal and lossy*: a basis invariant under
more than its domain shows will be refused pairings it could have
admitted, once the frame checks the two halves
(:ref:`manifold-invariance-pairing`). The remedy is to **declare the
finer domain** — a Legendre basis on :math:`S^2/SO(2)_x` rather than the
full harmonics on :math:`S^2` — which is tracker 3.4, and which is the
level-2 half of ERR-080's repair.

⛔ **The remedy is never to override the property**, and the type
enforces that: an override lets ``domain`` and ``invariance_group``
disagree, which is precisely the two-homes-for-one-fact state the
derivation exists to make unspellable. Hence ``@final`` — `[M]`
``Basis.__dict__['invariance_group'].fget.__final__`` is ``True``. It is
the same argument that keeps
:attr:`~orpheus.numerics.measure.DiscreteMeasure.quotient_group` derived
rather than stored, applied to the object on the other side of the
frame.

.. _manifold-invariance-pairing:

The pairing, measured: ERR-080 as a lattice verdict
-----------------------------------------------------

With both operands in hand, the check ERR-080 needs can finally be
*written down*. The rule is a containment in the subgroup lattice:

.. math::

   \text{admissible}
   \quad\Longleftrightarrow\quad
   \underbrace{G_{\text{spent}}}_{\texttt{measure.quotient\_group}}
   \;\subseteq\;
   \underbrace{G_{\text{have}}}_{\texttt{basis.invariance\_group}} ,

read as: *the symmetry a rule folded away must be one the basis's
functions are blind to.* If the rule kept one representative per
:math:`H`-orbit, then any function that distinguishes points within an
orbit has been handed a sample that cannot see the distinction, and the
pairing is a forgery whatever its shapes say.

Measured on the objects that ship, and on the pairing the tree
**actually forms**: each rule against the basis its own
``angular_frame(2)`` binds. The denominator is again all five
``Quadrature`` factories, and `[M]` **exactly one of the five fails**:

.. list-table:: `[M]` 2026-09-01 — ``rule.measure`` vs ``rule.angular_frame(2).basis``
   :header-rows: 1
   :widths: 26 22 18 34

   * - Rule
     - Basis its frame binds
     - SPENT / HAVE
     - ``have.contains(spent)``
   * - ``folded_product(4, 8)``
     - ``MirrorEvenSphericalHarmonicBasis`` (``mirror_axis=1``)
     - ``Mirror('y')`` / ``Mirror('y')``
     - ✅ **True** — and the two are the *same object*
       (``have is spent``), because ``basis.domain is measure.support``
       is one memoised :class:`Quotient`
   * - ``gauss_legendre(8)``
     - ``SphericalHarmonicBasis``
     - ``SO2('x')`` / ``Trivial``
     - ⛔ **False** — ERR-080, as a lattice verdict
   * - ``lebedev(17)``
     - ``SphericalHarmonicBasis``
     - ``None`` / ``Trivial``
     - spent nothing, so nothing to contain — admitted
   * - ``level_symmetric(8)``
     - ``SphericalHarmonicBasis``
     - ``None`` / ``Trivial``
     - the same
   * - ``product(4, 8)``
     - ``SphericalHarmonicBasis``
     - ``None`` / ``Trivial``
     - the same

The first row is the whole design in one line: the fold's two halves do
not merely *agree*, they read one group object out of one manifold, so
no drift between them is representable. The second row is ERR-080 —
stated, for the first time in this corpus, as a verdict a predicate
could return rather than as a story about zero-padded columns. And the
shape of the table matters as much as the verdict in it: the defect is
**not** a general weakness of the harmonic frame — `[M]` **1 of the 5**
shipped rules fails its own frame's pairing, and it is the 1-D one.
⚠ That is a *different* denominator from ERR-080's own scope census,
which counts ``(constructor, order)`` rows (`[M]` 7 of 15 non-zero, 5 of
them this defect); the two are not comparable and neither implies the
other.

.. note::

   ⛔ **Nothing refuses on this verdict yet, and the reason is worth
   stating precisely.** The frame's pairing gate — the plan's **G2**,
   fused with its siblings at tracker 2.2 — is not written. And a gate
   written naively on the *frame's* measure would be **inert on the very
   row it is for**: `[M]` 2026-09-01,
   ``Quadrature.gauss_legendre(8).angular_frame(2).measure.support.name``
   is still ``'S^2'`` (the surviving 1-D forgery,
   :ref:`manifold-err-080`) while
   ``Quadrature.gauss_legendre(8).measure.support.name`` is
   ``'S^2/SO2_x'``. That is why the gate below reads the verdict off the
   **quadrature's** measure — the object that knows what it spent — and
   why the negative leg is a *measurement made spellable*, not a
   refusal (``plan-authoring`` §6c: a step that adds a gate must land
   with the case the gate catches, and this step deliberately adds the
   operand rather than the gate).

   ⛔ **ERR-080 is OPEN.** It is held by the ``xfail(strict=True)`` gate
   at
   ``tests/sn/solve/test_pl_order_does_not_move_the_infinite_medium_flux.py``,
   three rows red by design. Nothing on this page repairs it.

**The campaign plan's "Part IV" lattice table, as a test.** That
four-row admissibility table — the same section of the #429 plan quoted
at :ref:`manifold-err-080-is-a-section` — was the *done-when* for this
step, and it now runs on real objects rather than on names
(``test_e5_part_IV_lattice_table_runs_on_the_objects_that_ship``):

.. list-table::
   :header-rows: 1
   :widths: 10 30 26 34

   * - Row
     - Basis space (HAVE)
     - Rule (SPENT)
     - Verdict
   * - 1
     - full :math:`S^2` harmonics — ``Trivial``
     - slab — ``SO2('x')``
     - ⛔ refused: this is ERR-080's pairing, and the refusal is
       categorical — no tolerance is involved
   * - 2
     - :math:`S^2/SO(2)_x` — ``SO2('x')``
     - slab — ``SO2('x')``
     - ✅ the repair
   * - 3
     - :math:`S^2/SO(2)_x` — ``SO2('x')``
     - full sphere — ``None``
     - ✅ a smaller space on a full rule is legal
   * - 4
     - :math:`S^2/\langle\sigma_y\rangle` — ``Mirror('y')``
     - fold — ``Mirror('y')``
     - ✅ the shipped fold

⚠ Two readings of that table are worth pinning. Row 3's measure side is
spelled ``None`` and not ``Trivial``, for the reason
:ref:`manifold-invariance-three-arms` gives: a full-sphere rule has SPENT
nothing, and the lattice element ``None`` stands for on that side is
:math:`\{e\}`, which every group contains. And rows 2 and 3 need a basis
on :math:`S^2/SO(2)_x` — tracker 3.4's Legendre basis, which `[M]` does
**not** ship (``class LegendreBasis`` is **0** hits over ``orpheus`` and
``tests``); the gate stands one in with a test-local stub declaring that
*domain*, carrying an explicit retirement trigger in its docstring. That
is a fixture with an expiry date, not a permanent one.

**The gates.** Section E of ``tests/numerics/test_basis_domain.py``,
`[M]` **six** functions and **eleven** collected rows — the module went
13 rows → 24, and the V&V matrix's ``numerics/test_basis_domain`` row
reads the same +11 independently — all ``@pytest.mark.foundation``
(the property is a type law, not a solver claim, and carries no theory
equation label):

.. list-table::
   :header-rows: 1
   :widths: 24 76

   * - Gate
     - What it pins
   * - ``test_e1``
     - ⭐⭐ the keystone, both legs: the fold's two halves read ONE
       group object (``==``, then ``is``), and the slab's pairing is
       **refusable** — asserted on the quadrature's measure, for the
       reason in the note above.
   * - ``test_e2``
     - the full-sphere harmonics HAVE ``Trivial`` at
       :math:`L \in \{0,1,3,7\}` — :math:`L = 0` included on purpose,
       as the lower-bound witness.
   * - ``test_e2b``
     - the mirror-even basis HAS its mirror **by identity**
       (``is domain.by``), over all three axes, with a negative leg
       showing a different mirror is incomparable — so the answer moves
       with the axis rather than being one constant.
   * - ``test_e3``
     - the category leg: ``None`` on every non-angular basis, including
       both delegating wrappers — with a positive control (the same
       class shape on a sphere answers ``Trivial``), so the arm is not
       "everything is ``None``" (``vv-principles`` #11).
   * - ``test_e4``
     - the three arms on one class shape differing only in ``domain``,
       with the quotient arm exercised through a **second group family**
       (``SO2``), since the shipped fold basis only ever brings a
       ``Mirror``.
   * - ``test_e5``
     - Part IV's four-row table above, on shipping objects.


.. _manifold-gotchas:

Gotchas
=======

.. _manifold-import-cycle:

The module imports nothing from ``numerics`` at runtime — on purpose
--------------------------------------------------------------------

:mod:`orpheus.numerics.manifold` references
:class:`~orpheus.numerics.symmetry.SubgroupOfO3` under
:data:`typing.TYPE_CHECKING` only. That is **load-bearing, not
tidiness**, and it is the kind of constraint that is invisible until it
is violated.

⚠ What makes the guard *affordable* is not that the module never touches
a group — it does, and increasingly: `[M]` the catalogue builders read
:attr:`group.name <orpheus.numerics.symmetry.SubgroupOfO3.name>`,
:attr:`group.mirror_axis
<orpheus.numerics.symmetry.SubgroupOfO3.mirror_axis>` and, since tracker
2.4, :attr:`group.rotation_axis
<orpheus.numerics.symmetry.SubgroupOfO3.rotation_axis>`. It is that the
group always arrives **as an argument**, never as a name this module has
to resolve — so every one of those reads is duck-typed at runtime and
needs no import. A design that instead *constructed* a group here (say,
to normalise a caller's tag) would force the import and close the cycle.

`[M]` by AST, **re-measured 2026-09-01 after tracker 2.4**, over the
three modules, with relative imports resolved and ``TYPE_CHECKING``
bodies separated:

.. list-table:: Every ``manifold`` / ``measure`` / ``symmetry`` edge among the three
   :header-rows: 1
   :widths: 26 24 20 30

   * - Site
     - Edge
     - Scope
     - Note
   * - ``manifold.py:71``
     - ``manifold → symmetry``
     - ``TYPE_CHECKING``
     - **The guard.** The only edge out of this module, and it is
       erased at runtime.
   * - ``measure.py:89``
     - ``measure → manifold``
     - module, **runtime**
     - Landed at tracker 2.0c, when ``support`` became a
       :class:`~orpheus.numerics.manifold.Manifold`.
   * - ``measure.py:108``
     - ``measure → symmetry``
     - ``TYPE_CHECKING``
     - Annotation only.
   * - ``measure.py:1141``
     - ``measure → symmetry``
     - **function** scope
     - Inside :meth:`DiscreteMeasure.quotient
       <orpheus.numerics.measure.DiscreteMeasure.quotient>`; its comment
       already says why.
   * - ``symmetry.py:102``
     - ``symmetry → manifold``
     - module, **runtime**
     - ⭐ **New at tracker 2.4** — ``_polar_axis_of`` needs
       :class:`~orpheus.numerics.manifold.Quotient` to read a polar
       marginal's axis off its support. `[M]` **re-measured
       2026-09-02**: the same line now imports **three** names,
       ``AXIS_INDEX``, ``Quotient`` and
       :func:`~orpheus.numerics.manifold.barycentre` — tracker 2.3
       moved the axis-index table here (:math:`\mathbb{R}^3`'s axes are
       a convention of *this* module's curved members) and made the
       polar embedding read the typed map (:ref:`manifold-barycentre`).
       The edge's direction and scope did not move; only its payload
       grew, which is the direction that keeps the guard affordable.
   * - ``symmetry.py:103``
     - ``symmetry → measure``
     - module, **runtime**
     -
   * - ``exactness.py:115``, ``:116``
     - ``exactness → manifold``
     - module, **runtime**, **twice**
     - ⭐ **A second pair, first measured 2026-09-02.**
       :mod:`orpheus.numerics.exactness` imports
       :class:`~orpheus.numerics.manifold.Manifold` and then
       :data:`~orpheus.numerics.manifold.CIRCLE` /
       :data:`~orpheus.numerics.manifold.SPHERE`, both at module scope.
       So a ``manifold → exactness`` edge — which is what carrying an
       orbit space's pushforward *reference measure* on the catalogue
       entry would want — closes a **two-hop** cycle of its own
       (:ref:`manifold-arrows-not-built`, item 2).

⟹ **the guard is now strictly more load-bearing than when it was
written.** Before tracker 2.4 the loop a runtime ``manifold → symmetry``
edge would have closed was the three-hop
``measure → manifold → symmetry → measure``. Now ``symmetry`` imports
``manifold`` **directly**, at module scope, so the same edge closes a
**two-hop** cycle ``manifold ⇄ symmetry`` as well. Two independent
cycles, one guard.

`[M]` demonstrated on a throwaway three-module package carrying exactly
the table's edge topology — no production file touched, which is what
makes the demonstration safe to run at all (mutating a module and
restoring it is not crash-safe). With the guard **removed**, i.e. a
module-scope ``manifold → symmetry`` added:

.. code-block:: text

   ImportError: cannot import name 'Quotient' from partially
   initialized module 'pkg.manifold' (most likely due to a circular
   import)

⭐ **and it fails in all three import orders** — ``measure`` first,
``manifold`` first, ``symmetry`` first. That is the sharpest change
tracker 2.4 made to this guard. The pre-2.4 three-hop cycle was
**order-dependent**: importing the right module first walked the loop in
an order that happened to work, so a smoke test could report green on a
broken façade. A two-hop cycle has no such order, so the failure is now
unconditional — worse to introduce, and cheaper to detect. With the
``TYPE_CHECKING`` guard restored — the shipped shape — all three
modules import cleanly in all three orders (the positive control,
without which a clean reading carries no information).

⚠ **Two ways this is easy to get wrong, both measured.**

1. **An AST import census that filters on**
   ``node.module.startswith("orpheus")`` **silently drops every
   relative import**, because an ``ImportFrom`` with ``level > 0``
   carries an *unqualified* ``.module``. A census written that way
   reports ``symmetry → measure`` as **absent** and concludes there is
   no cycle — wrong, and in the reassuring direction. Resolve relative
   imports, and give the census a positive control.
2. ⛔ **This item read "the cycle is not live today —** ``measure``
   **does not import** ``manifold`` **yet"** and was written on
   2026-08-31, when it was true and the guard was purely prophylactic.
   `[M]` it is false since tracker 2.0c: ``measure.py:89`` is a
   module-scope runtime ``measure → manifold``, and tracker 2.4 added
   ``symmetry.py:102`` on top. **Both cycles are live now.** A
   module-scope ``manifold → symmetry`` edge added today does not "pass
   every test and become fatal later" — it breaks
   ``import orpheus.numerics.measure`` immediately. The lesson survives
   its own correction, in the opposite direction: a prophylactic guard
   stops being prophylactic without anything editing it, so "is this
   still just precaution?" is a question with a shelf life.

.. _manifold-gotcha-ambient:

Topological dimension is not ambient dimension
-----------------------------------------------

``dim`` is what the manifold *is*; the ambient count is how many
columns :meth:`contains <orpheus.numerics.manifold.Manifold.contains>`
consumes. They differ for every curved member — a sphere is
``dim == 2`` in ``3`` ambient coordinates — and a
:class:`~orpheus.numerics.manifold.Product` needs the *ambient* count
to know where to split a point's coordinates, not the topological one.
The module keeps the ambient count in a deliberately **exhaustive**
``match`` with a raising fall-through, so a new member that forgets it
fails loudly rather than silently mis-splitting a product's
coordinates; a foundation gate walks every shipped variant through it.

⭐ **A** :class:`~orpheus.numerics.manifold.Quotient` **is where the two
counts come apart hardest, and the type answers them differently on
purpose.** Its ambient count — what ``_ambient`` reports, and therefore
what a :class:`~orpheus.numerics.manifold.Product` uses to split a
point's coordinates — is the **realization's**, because a product
factor must have one canonical width or the split is ambiguous. Its
``contains``, by contrast, accepts **either** coordinate system and
dispatches on the width it is handed
(:ref:`manifold-two-coordinate-systems`). So for the shipped fold `[M]`
``_ambient`` is :math:`2` (the disk) while ``contains`` also accepts
the :math:`(16,3)` section points the tree's own measures carry.

⚠ That asymmetry is deliberate and it is the one place on this page
where a single object answers two questions with two different numbers.
Read it as: *the canonical coordinate for composition is the chart's;
the predicate is as wide as the honest languages the object has.* An
earlier version of this paragraph stated only the first half — *"the
ambient count is the realization's, because membership is decided in
the coordinates the chart lands in"* — which was true of the
single-slot type and became false for ``contains`` on 2026-08-31.

.. _manifold-gotcha-shape-vs-verdict:

A wrong ambient dimension is a refusal, not a ``False``
--------------------------------------------------------

``SPHERE.contains(np.zeros((4, 2)))`` raises :exc:`ValueError` naming
the expected ambient dimension. It does not return ``False``. The
distinction is the difference between *"these are not points of this
manifold"* — a verdict about the data — and *"you handed me something
that cannot be points of anything with this ambient dimension"* — an
error in the caller. Collapsing the two would let a shape bug read as a
membership failure, and a membership predicate that returns ``False``
for malformed input is a predicate whose ``False`` means nothing.

.. _manifold-gotcha-not-a-manifold:

``EnergyGroups`` and ``IndexSet`` are 0-dimensional, and that is honest
-----------------------------------------------------------------------

A finite index set carries no metric structure and no smooth structure;
calling it a manifold is a stretch that the type makes deliberately.
The justification is that the *algebra* is what level 1 is for — an
energy axis composes with a spatial one under :math:`\times`, and a
measure on the pair is the tensor product — and a "manifold" that
refused to admit the counting factors would force every composite
support back into a string. ``dim`` is ``0`` for both, which is the
correct topological dimension of a finite discrete set, and neither
carries a chart.


.. _manifold-seams:

What is NOT built (the standing seams)
======================================

Stated explicitly so no reader mistakes a shipped *type* for a shipped
*migration*, and so the next phase does not re-derive a decision
already taken. The whole of the first row is what makes this page a
description of a capability rather than of a repair.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Not built
     - Where it lands, and what stands in for it today
   * - ⛔ ~~**Any production consumer at all**~~
     - ✅ **REMEDIED 2026-09-01 (tracker 2.0c).** *(Recorded as written —
       `[M]` 2026-08-31: "the only importers of
       :mod:`orpheus.numerics.manifold` are its own test module;
       ``Space = str`` is still live at ``measure.py:111`` with its six
       ``SPACE_*`` aliases, ``DiscreteMeasure.support`` is still a
       ``str``".)* The alias and all six tags are retired; ``support`` is
       a ``Manifold`` on all six implementors, the tensor product and the
       fold route through :meth:`Manifold.__mul__` and
       :meth:`Manifold.quotient`, and :attr:`DiscreteMeasure.phase`
       dispatches on the manifold's TYPE instead of on string prefixes.

       ⚠ **ERR-080 is still open** — held by the same
       ``xfail(strict=True)`` gate at
       ``tests/sn/solve/test_pl_order_does_not_move_the_infinite_medium_flux.py``.
       Retyping the slot is what makes the refusal *spellable*; it does
       not make it *fire*. Nothing calls :meth:`contains` on the way in,
       so the forged measure is still constructible — that is the row
       below, and the fused fix step.
   * - ``Basis.domain``
     - ✅ **LANDED 2026-09-01 (2.1).** No
       :class:`~orpheus.numerics.basis.base.Basis` could state the
       manifold its functions consume, which is why the ERR-080 pairing
       had nothing to check (:ref:`manifold-three-levels`). It is now an
       abstract property, so a basis that cannot say what it eats
       refuses to be constructed, and all six shipped subclasses answer.
       :class:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis`
       takes it as a **constructor field** (``partition_of``) rather
       than deriving it, and the reason was measurable in advance:
       `[M]` by AST, of 18 ``IndicatorBasis(...)`` construction sites
       tree-wide **4 are in** ``orpheus/``, and those four partition
       **three different manifold families** — a finite index set
       (``frame.py``, paired with ``support=f"index({axis_label})"``
       three lines below), :math:`\mathbb{R}^d` at two ranks
       (``geometry/mesh.py``), and the energy counting set
       (``data/energy_grid.py``). Any value the class *derived* from its
       own fields would hard-code one of the three. ⭐ The prediction
       held, and execution added a fourth family the string tag had
       hidden: `[M]` a partition by energy **VALUE** in eV is an
       :class:`Interval`, not the :class:`EnergyGroups` **index** axis
       production partitions — both ambient dimension 1, so only naming
       the point set separates them.
   * - ``Basis.invariance_group``
     - ✅ **LANDED 2026-09-01 (2.1b), and DERIVED — this row exists to
       record that no slot was added.** The pairing ERR-080 needs has two
       operands and the basis's was missing; the tracker recorded the
       property as *absent and derivable*, which invited a second
       abstract property with six overrides, and the phase opener found
       the answer already sitting in ``domain.by``.
       A function on :math:`M/H` *is* an :math:`H`-invariant function, so
       the group is read by a ``match`` on the domain's TYPE — `[M]`
       **6 of 6** shipped bases answer, with **0** subclass edits and one
       ``@final`` property on the ABC. This is tracker 2.0d's
       ``quotient_group`` FIELD dissolving into :attr:`Quotient.by
       <orpheus.numerics.manifold.Quotient.by>` at 2.0c, replayed one
       level over (:ref:`manifold-basis-invariance-group`).

       ⛔ **Still not built: the CONSUMER.** Both operands exist and
       `[M]` the verdict is computable — ``Trivial ⊇ SO2('x')`` is
       ``False`` for the slab — but nothing refuses on it. The frame's
       pairing gate is tracker **2.2**, fused with its siblings, and
       `[M]` a gate written on the *frame's* measure would be inert
       today because that measure still carries the forged ``S^2``
       (:ref:`manifold-invariance-pairing`). ERR-080 remains open.
   * - ``FunctionSpace.manifold``, and the derived ``L2[...]`` name
     - **Still open, and the reason narrowed.** Two sites build a level-2
       name by interpolating a level-1 tag, and both now interpolate a
       typed :class:`Manifold`: `[M]` ``measure.py:371`` is
       ``f"L2[{self.support.name}]"`` (2.0c retyped ``support``) and
       ``basis/indicator_basis.py:355`` is
       ``f"L2[coarse_cells({self.domain.name})]"`` (2.1 gave the basis a
       ``domain``). ⛔ This row read *"one of them
       (``basis/indicator_basis.py:284``) **hard-codes** it and `[M]` is
       already **false** for the energy-grid basis"* until 2026-09-01 —
       true when written, repaired by 2.1
       (:ref:`manifold-string-algebra`). What remains is the seam
       itself: a ``FunctionSpace`` still records the **string**, so the
       two producers agree by discipline rather than by construction, and
       a space that carried its own manifold would collapse both
       spellings into one.
   * - The two MAPS — the ``chart`` and the section — and the
       pushforward measure, as entry fields
     - **Phase 1.1 / 3.1, and NARROWED at 2.3.** `[M]` 7 of the
       derivation procedure's 9 outputs are slots on
       :class:`~orpheus.numerics.manifold.Quotient` today. What ships
       of each map is only its *end*: the chart's **codomain**
       (``realization``) and the section's **image**
       (``fundamental_domain``); neither map itself is a slot — which
       is why ``Quotient.contains`` must accept both languages rather
       than normalising to one
       (:ref:`manifold-two-coordinate-systems`).

       ⭐ **What tracker 2.3 changed (2026-09-02) is that a map is now
       expressible.** ⛔ This row read *"nothing can currently apply
       one"* until then; that is now false in general —
       :class:`~orpheus.numerics.manifold.ManifoldMap` is a value type
       and three arrows ship — and still true of *these two*. `[M]`
       none of the three is the entry's chart or its section, for
       three different structural reasons enumerated at
       :ref:`manifold-arrows-not-built`, and
       :attr:`~orpheus.numerics.manifold.Quotient.fundamental_domain`
       still has **zero** readers outside
       :mod:`orpheus.numerics.manifold` itself.

       The pushforward measure is not a slot either: a measure
       descends via :meth:`DiscreteMeasure.quotient
       <orpheus.numerics.measure.DiscreteMeasure.quotient>`, which
       knows nothing about the catalogue
       (:ref:`manifold-engine-data-model`). ⚠ `[M]` 2026-09-02 that one
       is blocked on an **import cycle**, not on a design question: a
       module-scope ``manifold → exactness`` edge kills 5 of 5 fresh
       import orders, so the field must be populated inside the
       derivation function.
   * - A ``ManifoldMap`` for the ERR-080 forgery arm
     - ⛔ **Deliberately NOT built, and it is the point.**
       ``Quadrature._harmonic_frame_measure``'s 1-D arm computes the
       orbit barycentre and declares it on :math:`S^2`; `[M]`
       2026-09-02 its nodes are ``np.array_equal`` to
       :func:`~orpheus.numerics.manifold.barycentre`'s image. It stays
       a **raw** :class:`~orpheus.numerics.measure.DiscreteMeasure`
       constructor because routing it through
       :meth:`~orpheus.numerics.measure.DiscreteMeasure.pushforward`
       would force it to name ``Ball(3)`` — i.e. would repair
       ERR-080's level-1 half inside a step whose subject is the type
       system. Retired at tracker **3.4**, together with the ``if``
       (:ref:`manifold-barycentre`).
   * - The remaining catalogue entries
     - **Phase 1.1.** `[M]` **six keys** ship — ``(Sphere,
       "SO2_x"/"SO2_y"/"SO2_z")`` and ``(Sphere,
       "sigma_x"/"sigma_y"/"sigma_z")`` — served by **two**
       procedures, since each family shares one derivation that reads
       the axis off the group. The identity quotient is a seventh
       answer, derived rather than tabulated. ⛔ This row read **four
       keys**, with a single ``(Sphere, "SO2")``, until tracker 2.4
       parameterised the axial rotation group on 2026-09-01; note the
       *procedure* count did not move
       (:ref:`manifold-so2-axis-is-a-parameter`). The expected
       remainder covers :math:`\mathbb{Z}_2` antipodal, :math:`C_n` /
       :math:`D_n` about an axis, the :math:`O_h` sublattice for octant
       symmetry, :math:`SO(3)`, and :math:`SO(2)\times\mathbb{R}_z` for
       the 1-D cylinder. ⚠ Whoever adds a :math:`C_n` entry must leave
       ``fundamental_domain=None``: a closed sector is **not** a
       strict fundamental domain, and the ``dim`` gate cannot catch a
       wrong one (:ref:`manifold-chart-section-asymmetry`).
   * - Collapsing the twin lookup
     - ✅ **DONE at tracker 2.4, 2026-09-01.**
       :attr:`AngularSymmetry.support
       <orpheus.numerics.quadrature.registry.AngularSymmetry.support>`
       no longer holds a table: it calls ``SPHERE.quotient(spent)``, so
       `[M]` for a slab its answer *is* the catalogue's object, by
       identity. The shipped collapse keeps the **orbit space**, where
       this row predicted ``…​.realization.name`` — a difference that
       matters, since taking the realization is the axis-blind step
       (:ref:`manifold-twin-lookup`, reading (iv)). The ``Trivial`` row
       remains two producers, deliberately.
   * - The ``support`` tag's own vocabulary split
     - ✅ **DONE at trackers 2.0c + 2.4.** `[M]` today
       ``gauss_legendre(8).measure.support.name`` is ``'S^2/SO2_x'``
       and ``folded_product(4,8).measure.support.name`` is
       ``'S^2/sigma_y'`` — **both** the quotient's name, both typed
       :class:`~orpheus.numerics.manifold.Quotient` values. ⛔ This row
       read *"``gauss_legendre(8).measure.support`` is ``'[-1,1]'`` —
       the realization's name … the registry's stage-0 gate compares it
       by string equality"*. The register split is closed and the gate
       is a value comparison. ⚠ What SURVIVES is the *disagreement it
       exposed*, which was never a vocabulary problem: `[M]`
       ``GEOMETRY_ANGULAR_SYMMETRY["cylinder"].admits_domain`` on the
       shipped fold is still **False**, because a rule folded by a
       member of :math:`\Gamma` genuinely does not live on the geometry's
       :math:`S^2/G^0`. Latent only because ``folded_product`` is not a
       registered spec (:ref:`manifold-twin-lookup`).
   * - The derivation ENGINE
     - **Deferred, not refused** — the ruling, the falsifiable
       compliance check and the acceptance suite that is already
       written are at :ref:`manifold-engine-seed`.
   * - Renaming ``DiscreteMeasure.support``
     - **With the migration.** The slot names the ambient manifold,
       not :math:`\operatorname{supp}(\mu)`; the corpus already
       records the misnomer
       (:doc:`/theory/foundations/discrete_measures`). Renaming a slot
       that `[M]` **87** ``support=`` keyword arguments pass — 29 in
       ``orpheus/``, 58 in ``tests/`` — is a migration act, not a docs
       one.
   * - An ``automodule`` for this module
     - **Not scheduled, deliberately.** `[M]` 2026-08-31: of 48
       ``automodule`` directives in the doc source, **6** are
       ``orpheus.numerics.*`` — ``axis``, ``convergence``,
       ``coupled_system``, ``eigenvalue``, ``field``, ``functional``.
       ``manifold``'s two closest siblings in the three-level stack,
       :mod:`orpheus.numerics.measure` and
       :mod:`orpheus.numerics.space`, are **not** among them, so
       surfacing level 1 alone would make
       :class:`~orpheus.numerics.manifold.Manifold` a live link while
       :class:`~orpheus.numerics.space.FunctionSpace` beside it in the
       same sentence stays plain text. Surfacing the package is its
       own task. ⚠ Consequence for editors: the Python-domain
       cross-references on this page render as **plain text with no
       warning at any severity**, so a stale one is invisible to
       ``sphinx -W`` and must be caught by an import-resolution grep.


.. _manifold-verification:

Verification
============

The gates live in ``tests/numerics/test_manifold.py``: `[M]` **50 test
functions, 70 collected rows**, run under the canonical
``python -O -m pytest`` invocation. Several functions are parametrized
— over the shipped-variant list, over three bases, and (since tracker
2.3) over the three rotation axes — and the two counts are given
separately because they move for different reasons: adding a *variant*
moves the second and not the first.

⚠ The row count here is the one the generated V&V matrix reports for
this module, which is a second, independent reading of the same tree
(``docs/theory/verification/matrix.rst``, ``numerics/test_manifold``
row). `[M]` it read **44** before the two-slot ruling, **56** after,
and **70** after tracker 2.3 added ``TestManifoldMap``. An earlier
version of this paragraph said *"30 test functions, 40 collected
rows"*; both numbers were wrong when written — the module had 32
functions and 44 rows — which is why the count is now stated with the
instrument that produces it.

**Every one carries** ``@pytest.mark.foundation`` **and none carries**
``verifies(...)``, and that is the correct tier rather than an
omission. ``foundation`` is the V&V ladder's *orthogonal* category —
software invariants with no theory-page equation label behind them:
data-structure laws, factory outputs, algebraic reduction invariants.
The assertions here are the intrinsic laws of the type (dimension
additivity, membership, the quotient's dimension drop, the recorded
derivation), not an L0–L3 claim about a solver, a flux or an
eigenvalue. Tagging them ``verifies`` would mint a coverage edge that
an audit would then trust.

Seven groups, and what each is for:

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Group
     - What it pins
   * - The type's own laws
     - The base is uninstantiable (a sum, not a member); every variant
       answers the three total operations; variants are frozen and
       compare by value; ``dataclasses.replace`` **re-runs the
       construction invariant** rather than being a hole in it; the
       names reproduce the retired string tags verbatim.
   * - The product algebra
     - :math:`\dim(M\times N) = \dim M + \dim N` over **every ordered
       pair** of the shipped members; the name reproduces the retired
       interpolation; multiplying a non-manifold is refused;
       membership splits the coordinate blocks.
   * - Membership, **both legs**
     - `[M]` **8 of the 9** tests in the ``TestMembership`` class
       assert a positive verdict beside their negative one
       (``vv-principles`` #11 — a contract predicate tested only
       against a broken instance validates the *raising*, not the
       *claim*); the ninth is
       ``test_a_wrong_ambient_dimension_is_a_typed_refusal``, which
       asserts a **raise** rather than a verdict and so has no positive
       leg to carry (:ref:`manifold-gotcha-shape-vs-verdict`). The
       load-bearing row is the ERR-080 forgery:
       the negative leg refuses :math:`(\mu,0,0)`, the positive leg
       admits the same nodes normalised, and a third assertion places
       them on :math:`[-1,1]` where they belong. (The count is scoped
       to that class on purpose — the fold's own membership gates,
       below, carry their legs differently.)
   * - The recorded derivation
     - The symbolic regression tests of
       :ref:`manifold-engine-seed` — for **both** entries, the
       :math:`P` matrix, the determinant, the empty syzygy, and the
       stratum **solved for** rather than compared to a literal.
   * - ⭐ The :math:`\sigma_y` fold, on production data
     - The load-bearing gate of the two-slot ruling carries **three**
       legs on the shipped
       ``Quadrature.folded_product(4, 8).measure.nodes``: the section
       ADMITS them, and REFUSES both wrong inputs — the orbit twins
       (which ``realization = SPHERE`` would have admitted) and the
       ERR-080 forgery. A single-leg gate could not tell those two
       candidate designs apart, since ``SPHERE`` also refuses the
       forgery. A companion row asserts the chart is **Mode-12 blind**
       to the same forgery while the section is not — i.e. the
       blindness is pinned as a *property*, not merely avoided; and
       ``test_the_half_space_is_CLOSED_because_production_marches_from_it``
       pins the non-strict inequality against the march seed, the only
       witness available since the shipped rule never populates the
       stratum.
   * - The two coordinate systems agree
     - ``__post_init__``'s dimension law, with both legs
       (``vv-principles`` #11): a hemisphere offered against a 1-D
       realization is **refused at construction**, and every shipped
       entry **satisfies** it. A third row pins the rule that makes one
       field express both a half-space and a hyperplane —
       :math:`\dim = 2` for a lone normal, :math:`\dim = 1` once its
       antipode joins it.
   * - ⭐ The ARROWS (``TestManifoldMap``, tracker 2.3)
     - The type's own laws first — a map is a frozen value with two
       endpoints; composition is refused across mismatched endpoints
       (**both** legs, plus the ``TypeError`` on a non-map); the
       pushforward is **functorial**
       (:eq:`manifold-map-functoriality`), nodes and weights by
       ``np.array_equal``. Then the two named maps, each with a
       positive leg *and* a refusal leg: ``archimedes`` lands on
       :math:`S^2` for **all three** axes and collapses the fibre onto
       the stratum, and equals the direction-cosine triple spelled by
       hand; ``barycentre`` lands **inside** the ball and on the sphere
       **only at the poles** — which is the ERR-080 discriminator
       stated as a property of a *map* rather than of a quadrature —
       satisfies :math:`1-\lVert b\rVert^2 = 1-\mu^2`, and refuses a
       mirror quotient, the trivial quotient and a bare interval with
       :math:`SO(2)` as the positive control. A final row pins the
       Pattern-2 collapse: ``symmetry._embedded_nodes`` is
       ``np.array_equal`` to the map it now reads
       (:ref:`manifold-barycentre`).

⭐ Two exhaustiveness gates are worth naming separately, because they
are what make "closed sum" a checkable claim rather than a description.
``test_every_variant_is_reachable_from_this_modules_list`` compares
``Manifold.__subclasses__()`` against the module's own exercised list,
so a member added to the module but not to the tests fails **there**;
``test_ambient_dimension_is_defined_for_every_variant`` walks every
shipped variant through the exhaustive ``match``. The benefit of a
closed sum is that an operation can be checked against every member,
and that benefit is only real if the member list is itself pinned.

.. note::

   **A second module carries the CONSUMER-side gates**, and it is where
   this page's basis-facing claims are pinned:
   ``tests/numerics/test_basis_domain.py``, `[M]` **24** collected rows
   (was 13 before tracker 2.1b; the count is the generated V&V matrix's
   ``numerics/test_basis_domain`` row, the same independent instrument
   used above). Section D pins ``domain`` — every shipped basis answers,
   a basis that cannot say what it eats **cannot be constructed**, and
   the flagship ``test_d6`` pins *"the two halves of one frame name ONE
   manifold"*. Section E pins ``invariance_group``
   (:ref:`manifold-basis-invariance-group`). Every row there is
   ``@pytest.mark.foundation`` for the same reason as this module's:
   these are the type's own laws, not an L0–L3 claim about a flux.


.. _manifold-development-history:

Development history
===================

Reverse-chronological changelog of the architectural milestones of
*this page's* subject — the point-set layer. The space layer's own
changelog is at :ref:`spaces-development-history`. Entries marked
*(in development)* live on an unmerged feature branch and have no
landed merge-to-``main`` hash yet; **trust** ``git`` **over this table
for merge status.**

.. list-table::
   :header-rows: 1
   :widths: 10 50 12 28

   * - When
     - Architectural milestone
     - Issue
     - Where
   * - 2026-09-02
     - **The category gets its ARROWS, and a codomain stops being
       something a caller can assert.** Every construction that moved
       a point set had been applying a callable and then *naming the
       destination by hand* — which is the exact shape of
       :ref:`ERR-080 <manifold-err-080>`, since a destination named at
       the call site is a claim nothing can contradict.
       :class:`~orpheus.numerics.manifold.ManifoldMap` makes it a
       field: one frozen value type (``domain``, ``codomain``,
       ``apply``) with named maps as **factories**, ruled that way
       (user, 2026-09-02) for the same reason
       :data:`~orpheus.numerics.manifold.SPHERE` is a value and not a
       subclass. :meth:`DiscreteMeasure.pushforward
       <orpheus.numerics.measure.DiscreteMeasure.pushforward>` retires
       ``new_space=`` and **reads** its target, and additionally
       refuses a map out of the wrong point set — by manifold VALUE,
       which is what makes it discriminating: `[M]` the slab's
       :math:`S^2/SO(2)_x` rule and the chart rule on :math:`[-1,1]`
       have ``np.array_equal`` nodes and only one of them is accepted.
       Three arrows are typed — the **Archimedes** chart
       :math:`[-1,1]\times S^1 \to S^2` (named for the hat-box
       theorem; `[M]` :math:`\pi \circ \varphi_a = \mathrm{pr}_1`
       bit-exactly, and the product rule is now
       ``(polar * azimuthal).pushforward(archimedes("z"))``,
       bit-identical to its retired hand loop on **60 of 60**
       configurations with its support **identical** to the chart's
       codomain); the orbit **retraction** inside ``quotient()``, now
       landing on the catalogue's own object; and the orbit
       **barycentre** :math:`\mu \mapsto \mu\,\hat e_a`, whose honest
       codomain is ``Ball(3)`` because
       :math:`1 - \lVert\mu\hat e_a\rVert^2 = 1-\mu^2 = \tfrac14\det P`
       — it lands ON the sphere only at the stratum, so it is
       **canonical precisely because it is not a section**. ⭐ That
       gives ERR-080 a one-sentence statement in the type system's own
       vocabulary — *the barycentre map with a forged codomain* —
       and `[M]` the forgery's nodes are ``np.array_equal`` to the
       honest map's image, so what is false about it is a **type** and
       nothing else. The honest spelling
       (``symmetry._embedded_nodes``) now reads the map, collapsing a
       Pattern-2 twin (`[M]` bit-identical on 12 rows).
       ⛔ **An enabler, not a repair**: no membership check runs inside
       a map, the forgery arm stays a raw constructor **by design**
       until tracker 3.4, and the gate still declares three
       ``xfail(strict=True)`` rows (:ref:`manifold-arrows`).
     - `#429 <https://github.com/deOliveira-R/ORPHEUS/issues/429>`_
     - *(in development)* ``fix/angular-phantom-support``; tracker 2.3.
       ⚠ The code was **uncommitted in the working tree** when this row
       was written — trust ``git log`` over this cell for its hash.
   * - 2026-09-01
     - **A basis learns what it EATS — and therefore what symmetry it
       HAS.** :class:`~orpheus.numerics.basis.base.Basis` gained the
       level-1 slot the three-level table had listed as ⛔ *nothing*
       (:ref:`manifold-three-levels`): ``domain``, a :class:`Manifold`,
       abstract on the ABC so a basis that cannot say what it consumes
       **cannot be constructed**. That closed a live falsehood — `[M]`
       ``basis/indicator_basis.py`` hard-coded its coefficient space's
       name as ``f"L2[coarse_cells_R{ndim}]"``, so a 2-group **energy**
       basis and a 2-cell **spatial** basis compared ``==`` *and*
       hash-equal; they do not now
       (:ref:`manifold-string-algebra`). ⭐ Assigning the type was itself
       a census: it separated the continuous energy axis in eV
       (:class:`Interval`) from the multigroup *index* axis
       (:class:`EnergyGroups`), which the tag ``"energy"`` had conflated
       at equal ambient dimension. **Then the same slot answered a
       second question for free.** A function on :math:`M/H` *is* an
       :math:`H`-invariant function, so
       :attr:`~orpheus.numerics.basis.base.Basis.invariance_group` is a
       ``match`` on ``domain`` — `[M]` **6 of 6** shipped bases answer,
       ``@final``, **0** subclass edits, no new field. The tracker had
       recorded the property as *absent and derivable*, which invited six
       overrides; the phase opener dissolved them, exactly as tracker
       2.0d's ``quotient_group`` **field** had dissolved into
       :attr:`Quotient.by <orpheus.numerics.manifold.Quotient.by>` one
       step earlier. With both operands in hand
       the ERR-080 pairing became a **lattice verdict** — `[M]`
       ``Trivial ⊇ SO2('x')`` is ``False`` for the slab, while the
       shipped fold's two halves are literally one group object
       (:ref:`manifold-basis-invariance-group`). ⛔ **Nothing refuses on
       that verdict**: the frame's pairing gate is tracker 2.2, and
       ERR-080 stays open, held by its ``xfail(strict=True)`` gate.
     - `#429 <https://github.com/deOliveira-R/ORPHEUS/issues/429>`_
     - *(in development)* ``fix/angular-phantom-support``; trackers 2.1
       and 2.1b. ⚠ 2.1b was **uncommitted in the working tree** when this
       row was written — trust ``git log`` over this cell for its hash.
   * - 2026-09-01
     - **The axial rotation group gets its AXIS, and the type gets its
       first production consumer.** :math:`SO(2)` left the parameter-free
       enum and became ``SO2(axis)``, beside ``Mirror(axis)`` and for
       the same reason one month later: `[M]` the tree carries **two
       poles** — the real spherical-harmonic basis and the slab's polar
       marginal are about :math:`x`, every product rule's polar factor
       and the finite families are about :math:`z` — and one
       Gauss–Legendre rule serves both roles, so the group a marginal
       was quotiented by cannot be spelled without its axis
       (:ref:`manifold-so2-axis-is-a-parameter`). The catalogue went
       from four keys to **six**, still two procedures, because the
       :math:`SO(2)` derivation now reads its axis off the group exactly
       as the mirror one does (:ref:`manifold-s2-so2`). Downstream, the
       slab's rule **declares** its orbit space through a new measure
       verb, :meth:`on_orbit_space
       <orpheus.numerics.measure.DiscreteMeasure.on_orbit_space>` — same
       atoms, new support — so `[M]`
       ``gauss_legendre(8).measure.support.name`` is ``'S^2/SO2_x'`` and
       an 8-node angular space no longer compares equal to an 8-node
       spatial rule on the same interval, which it did before
       (:ref:`manifold-orbit-space-declaration`). ⭐ That collapsed the
       registry twin this page had listed as a seam: ``AngularSymmetry``
       now *calls* ``SPHERE.quotient`` (:ref:`manifold-twin-lookup`).
       ⚠ Two memos landed with it, and they are not optimisation
       polish: the catalogue derivation is ~6 ms of SymPy and every slab
       quadrature now carries one, and the icosahedral operator set was
       being rebuilt tens of times per invariance walk once the axial
       family offered three axes
       (:ref:`manifold-quotient-is-memoised`). ⛔ **ERR-080 remains
       open**; what this bought it is the *vocabulary* for its section,
       not the section
       (:ref:`manifold-the-axis-convention-for-a-section`).
     - `#429 <https://github.com/deOliveira-R/ORPHEUS/issues/429>`_
     - *(in development)* ``fix/angular-phantom-support``; tracker 2.4.
       ⚠ The code was **uncommitted in the working tree** when this row
       was written — trust ``git log`` over this cell for its hash.
   * - 2026-08-31
     - **An orbit space gets its second coordinate system, and the
       catalogue its second derivation.** Deriving the shipped
       cylindrical fold :math:`S^2/\langle\sigma_y\rangle` produced an
       object the single-slot type could not hold: because :math:`H` is
       **finite** the dimension does not drop, so the invariant chart
       buys no reduction and a *section* is canonical — while every
       measure the tree emits through ``.quotient(...)`` already speaks
       the section's coordinates. Ruled: **two slots**. ``realization``
       keeps its documented meaning (the chart's codomain, in the
       invariants' language) and ``fundamental_domain`` carries the
       section's image, in the base's; ``contains`` accepts both and
       dispatches on ambient width, ``_ambient`` still reports the
       chart's, and ``__post_init__`` gates their dimensions
       (:ref:`manifold-two-coordinate-systems`). Two variants were
       minted *by the derivation*, ``Ball`` and ``FundamentalDomain``,
       and ``singular_stratum`` was retyped from ``tuple[float, ...]``
       to a symbolic **locus** — the :math:`\sigma_y` stratum is a
       circle, and the first entry's shape had become the field's type
       (:ref:`manifold-stratum-is-a-locus`). Four candidate single-slot
       realizations were measured and refused, including the disk
       alone, which is **Mode-12 blind** to ERR-080
       (:ref:`manifold-realization-refuted`). ⭐ The finding that pays
       for the rest: **ERR-080's level-1 half is a botched section of**
       :math:`S^2/SO(2)` — a chart where a section was needed, faked by
       zero-padding off the sphere
       (:ref:`manifold-err-080-is-a-section`). Three mirror keys ship,
       one procedure.
     - `#429 <https://github.com/deOliveira-R/ORPHEUS/issues/429>`_
     - *(in development)* ``fix/angular-phantom-support``
       (``b55bba56``); tracker 1.1 for the entry, user ruling of
       2026-08-31 for the two slots
   * - 2026-08-31
     - **The manifold becomes an OBJECT.** Level 1 of the three-level
       stack stops being ``Space = str`` — an opaque tag whose own
       comment called the entries *"recommendations, not
       constraints"* — and becomes a closed sum with an algebra:
       ``dim`` / ``name`` / ``contains`` / ``__mul__`` total on the
       base, and the invariant-theoretic derivation fields on
       :class:`~orpheus.numerics.manifold.Quotient` alone. `[M]`
       **eight** variants at this commit (ten after the ruling above),
       two of them consolidating families the string
       vocabulary had spelled several ways
       (:ref:`manifold-members`). The three morphisms the mint names
       were already running as string interpolation at
       ``measure.py:588`` / ``:1022`` / ``:802``
       (:ref:`manifold-string-algebra`). The first catalogue entry
       ships, ``S^2/SO(2) = [-1,1]``, derived by the Procesi–Schwarz
       procedure and carrying its own symbolic regression tests
       (:ref:`manifold-s2-so2`) — which, because the entry's fields
       *are* the procedure's outputs, are the deferred derivation
       engine's acceptance suite, written before it
       (:ref:`manifold-engine-seed`). The identity quotient
       :math:`M/\{e\} = M` is **derived, not tabulated**, and doubles as
       a positive control on the machinery
       (:ref:`manifold-twin-lookup`). ⛔ **Type only: zero production
       consumers, and ERR-080 remains open**
       (:ref:`manifold-seams`).
     - `#429 <https://github.com/deOliveira-R/ORPHEUS/issues/429>`_
     - *(in development)* ``fix/angular-phantom-support``
       (``b8c05d16``); tracker 2.0a, user ruling D0.7 for the mint and
       2.0a-R for the shape


References
==========

* Procesi, C. and Schwarz, G. (1985). "Inequalities defining orbit
  spaces." *Inventiones Mathematicae* **81**, no. 3, 539–554,
  doi:10.1007/BF01388581. The theorem :eq:`manifold-procesi-schwarz`
  transcribes: the image of the orbit map is cut out of the syzygy
  variety by the single condition
  :math:`P \succeq 0` on the gradient Gram matrix of the invariants.
  ⚠ **Not held locally.** `[M]` 2026-08-31 the paper is absent from
  ``scratch/literature/`` (78 items, all reactor-physics and transport
  literature) and has no OCR sidecar, so no page- or equation-level
  verification against the scan was possible. The statement above is
  the standard form of the theorem and the volume/year are
  over-determined and consistent (*Inventiones* **81** is 1985); a
  session that acquires the paper should verify the theorem's own
  numbering and record it here.
* Schwarz, G. (1975). "Smooth functions invariant under the action of a
  compact Lie group." *Topology* **14**, 63–68. The result behind step
  2 — that the orbit map is proper and separates orbits, so smooth
  invariants factor through it.
* Weyl, H. (1946). *The Classical Groups: Their Invariants and
  Representations*, 2nd ed. Princeton University Press. Finite
  generation of the invariant ring (step 1), for the classical groups
  this corpus uses.
* Satake, I. (1956). "On a generalization of the notion of manifold."
  *Proceedings of the National Academy of Sciences* **42**, 359–363.
  The original definition of a V-manifold — what is now called an
  **orbifold** — which is what :math:`S^2/SO(2)_a` is and a quotient
  manifold is not (:ref:`manifold-singular-stratum`).
* Helgason, S. (1984). *Groups and Geometric Analysis*. Academic Press.
  Chapter IV (spherical functions on a Gelfand pair; the Funk–Hecke
  theorem as the :math:`(SO(3), SO(2))` instance), which is the
  literature behind :ref:`manifold-gelfand`.
* Hamermesh, M. (1962). *Group Theory and Its Application to Physical
  Problems*. Addison-Wesley. §2.5 (finite point groups) — the source
  :doc:`/theory/foundations/discrete_measures` already cites for the
  group side of the same construction.
