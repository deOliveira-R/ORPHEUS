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
      status: "MINTED, gated, and WIRED. Two catalogued derivations ship (S^2/SO(2) and S^2/<sigma_a> for the three mirrors, plus the derived identity quotient), and a Quotient carries BOTH coordinate systems after the 2026-08-31 two-slot ruling. `Space = str` and its six SPACE_* tags are RETIRED (tracker 2.0c, 2026-09-01): `DiscreteMeasure.support`, `GeneratingMeasure.support`, `UniformMeasure.support`, `ProductMeasure.support`, the `ReferenceMeasure` Protocol and `AngularSymmetry.support` all carry a Manifold, and `Basis.domain` does too (2.1). ERR-080 itself is still OPEN — the forgery's REFUSAL at construction is tracker 2.0b + the fused fix step, not this page's capability"


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
     (:eq:`manifold-procesi-schwarz`). For :math:`S^2/SO(2)` this gives
     :math:`P = \operatorname{diag}(1, 4p_2)`,
     :math:`\det P = 4p_2 = 4(1-\mu^2)` and the orbit space
     :math:`[-1,1]` (:eq:`manifold-s2-mod-so2`); for the shipped
     cylindrical fold :math:`S^2/\langle\sigma_a\rangle` it gives
     :math:`P = \operatorname{diag}(1,1,4p_3)` and the **closed unit
     disk** :math:`D^2` (:eq:`manifold-s2-mod-mirror`).
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
     :math:`S^2/SO(2)`. The realization is a chart; a consumer needed a
     section; the tree fabricated one by zero-padding to
     :math:`(\mu,0,0)`, which is off :math:`S^2` — `[M]` norms
     :math:`0.183\ldots0.960`, while an honest :math:`\varphi = 0`
     half-meridian is on the sphere to :math:`0.0`. ⚠ That names the
     **level-1** half only; the level-2 repair is still the trivial
     isotypic sub-basis :math:`\{Y_\ell^0\}\cong\{P_\ell\}`
     (:ref:`manifold-err-080-is-a-section`).
   - ⚠ **An orbit space is an ORBIFOLD, not a quotient manifold.**
     Where the action is not free, the image of the fixed-point set is
     a **singular stratum**. For :math:`S^2/SO(2)` that locus is
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
     :math:`(1-\mu^2)` is the squared :math:`SO(2)`-orbit radius
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
     runtime, and that is load-bearing.** `[M]` ``symmetry.py:98``
     imports ``measure`` at module scope, so once ``measure`` imports
     ``manifold`` a module-scope ``manifold → symmetry`` edge closes the
     cycle ``measure → manifold → symmetry → measure``
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
   different categories. The slot is therefore ``domain``, and it is
   :ref:`not yet built <manifold-seams>`.


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
   and since :attr:`DiscreteMeasure.phase
   <orpheus.numerics.measure.DiscreteMeasure.phase>` keys on that group, *the
   angular frame's own measure could not say it was angular*: it raised
   ``NotImplementedError`` on all twelve. ⟹ the transferable form, worth more
   than the instance: **"the rebuild loses X" is a completeness claim over the
   source type's FIELD LIST**, and its denominator is
   ``dataclasses.fields(T)`` — not the concept you happen to be chasing.

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

   * - Site
     - What it computes today
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

The remaining half is the *measure's* side: ``support`` is still a
``str``, so ``measure.py:331`` derives a correct name from an untyped tag.
A :class:`~orpheus.numerics.space.FunctionSpace` that carried its own
manifold would collapse both spellings into one — the level-2 half of
this repair, tracked at :ref:`manifold-seams`.

.. _manifold-string-drift:

(c) The vocabulary drifted, and a nonsense quotient is spellable
----------------------------------------------------------------

`[M]` both ``'S^2/<sigma_y>'`` (``tests/numerics/test_measure.py:701``,
a hand-written ``new_space=``) and ``'S^2/sigma_y'``
(``:777``, ``:947``, and what the production
:meth:`DiscreteMeasure.quotient
<orpheus.numerics.measure.DiscreteMeasure.quotient>` emits, since
``SubgroupOfO3.Mirror("y").name == "sigma_y"``) ship. They name **one**
quotient and are **unequal under** ``==``. Also shipped as support
literals: ``'img'``, ``'probe'``, ``'[-1,1]^slab'``.

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
(`[M]` ``None`` for ``SO2`` and ``Oh``, ``1`` for ``Mirror('y')``) and
which ``directional.py:522`` already branches on to raise
:exc:`NotImplementedError`. Repeating that tax on a brand-new type,
with six fields instead of one, is the design the closed sum refuses.

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

`[M]` the catalogue holds **four** keys today —
``Sphere/SO2``, ``Sphere/sigma_x``, ``Sphere/sigma_y``,
``Sphere/sigma_z`` — served by **two** procedures, because all three
mirrors share one derivation that reads the axis off the group. The
identity quotient :math:`M/\{e\} = M` is a fifth answer and is *not* a
table row: it is derived on the spot, for every manifold, because it is
a theorem (:ref:`manifold-twin-lookup`).

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

The worked entry: :math:`S^2 / SO(2)`
--------------------------------------

This is the **first** entry the tree catalogued, and it is the entry
every 1-dimensional angular discretisation is secretly using. Every
line below was **re-derived and re-run** in this session, independently
of the catalogue, and then compared against it.

**Step 1 — the invariants.** :math:`SO(2)` here is the group of
rotations about :math:`\hat z`, acting on :math:`\mathbb{R}^3` by

.. math::

   R_\theta =
   \begin{pmatrix}
     \cos\theta & -\sin\theta & 0 \\
     \sin\theta & \phantom{-}\cos\theta & 0 \\
     0 & 0 & 1
   \end{pmatrix}.

Two invariants generate:

.. math::

   p_1 = z, \qquad p_2 = x^2 + y^2 .

`[M]` verified symbolically for general :math:`\theta` — both satisfy
:math:`p(R_\theta x) - p(x) = 0` after ``simplify``, with the
non-invariant control :math:`x` correctly reported **not** invariant.
The control matters: a check that passes on everything is not a check.

**Step 3 — the syzygy ideal is empty.** The Jacobian

.. math::

   \frac{\partial (p_1, p_2)}{\partial (x, y, z)}
   =
   \begin{pmatrix} 0 & 0 & 1 \\ 2x & 2y & 0 \end{pmatrix}

has `[M]` generic rank **2**, equal to the number of invariants, so
:math:`p_1` and :math:`p_2` are algebraically independent and
:math:`I = (0)`. There are no equalities; the orbit space is cut out by
inequalities alone.

**Step 4 — the Procesi–Schwarz matrix.** With
:math:`\nabla p_1 = (0,0,1)` and :math:`\nabla p_2 = (2x, 2y, 0)`,

.. math::
   :label: manifold-s2-mod-so2

   P \;=\;
   \begin{pmatrix} 1 & 0 \\ 0 & 4 p_2 \end{pmatrix},
   \qquad
   \det P \;=\; 4 p_2 ,
   \qquad\text{so}\qquad
   \mathbb{R}^3 / SO(2) \;=\; \{\, p_2 \ge 0 \,\}.

The two invariants have orthogonal gradients everywhere, which is why
:math:`P` is diagonal and the condition collapses to a single
inequality.

**Step 5 — restrict to the sphere.** Adjoining
:math:`p_1^2 + p_2 = 1` and writing :math:`\mu = p_1 = \hat\Omega \cdot
\hat z`:

.. math::

   \det P \big|_{S^2} \;=\; 4\,(1 - \mu^2),
   \qquad
   S^2 / SO(2) \;=\; \{\, \mu \in \mathbb{R} : 1 - \mu^2 \ge 0 \,\}
   \;=\; [-1, 1].

`[M]` SymPy's ``solve_univariate_inequality`` on
:math:`4 - 4\mu^2 \ge 0` returns ``Interval(-1, 1)`` — the orbit space
is *computed*, not asserted — and the zero set is exactly
:math:`\{-1, +1\}`.

.. (vv-status rationale) manifold-s2-mod-so2 is the INSTANCE of
   :eq:`manifold-procesi-schwarz` for the one shipped catalogue entry.
   Its content IS checked, and tightly: the P-matrix and its
   determinant are recomputed symbolically and compared with
   sp.simplify, and the singular stratum is SOLVED for rather than
   compared to a literal, by
   tests/numerics/test_manifold.py::TestQuotient::{test_the_procesi_schwarz_matrix_matches_the_hand_derivation,
   test_det_P_vanishes_exactly_on_the_singular_stratum}.
   It is nonetheless `documented`, not verifies-covered, for a reason
   that is structural rather than a gap: those gates are
   @pytest.mark.foundation -- an intrinsic law of a data type, with no
   flux, eigenvalue or convergence claim behind it -- and
   vv-principles' foundation tier carries NO verifies(...) marker by
   rule. A verifies edge here would need a test-side change and would
   assert a claim class the gates do not make. If the equation ever
   acquires a solver-facing consumer (e.g. the descent of a kernel
   along the quotient), that consumer's gate is where the marker
   belongs.
.. vv-status: manifold-s2-mod-so2 documented

.. list-table:: My re-derivation against the shipped catalogue entry
   :header-rows: 1
   :widths: 44 56

   * - Check
     - Result
   * - :math:`P` (mine) :math:`-` ``entry.gram``
     - `[M]` ``simplify`` :math:`\to` the zero :math:`2\times2` matrix
   * - :math:`\det P` (mine) :math:`-` ``entry.det_gram``
     - `[M]` ``simplify`` :math:`\to 0`
   * - ``entry.realization``
     - ``Interval(a=-1.0, b=1.0)``
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
       :math:`S^2 \to S^2/SO(2)` is canonical, since every half-meridian
       is one and none is distinguished
       (:ref:`manifold-two-coordinate-systems`)
   * - ``entry.name`` / ``entry.derived_by``
     - ``'S^2/SO2'`` / ``'hand'``

⚠ One reproduction note, kept because it is the cheap trap in step 4.
Re-expressing :math:`P` in the invariants is a substitution, and
``sp.Matrix(...).subs(x**2 + y**2, p2)`` **silently fails** on
:math:`4x^2 + 4y^2`, which does not literally contain the node
:math:`x^2+y^2`. My first run therefore produced
:math:`\det P = 4x^2 + 4y^2`, an empty solution set for the stratum,
and a spurious disagreement with the catalogue. The failure was mine,
not the entry's; ``factor`` before substituting (or, as the shipped
builder does, substitute :math:`x \to \sqrt{p_2},\, y \to 0` — legal
because the expression is constant on an orbit) fixes it. A
disagreement with a reference is not a refutation until you have
diagnosed whose it is.

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
is a theorem about :math:`S^2/SO(2)`, not a conditioning accident.

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
space, stated in the vocabulary of :math:`S^2/SO(2)` rather than of the
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

1. **No section of** :math:`S^2 \to S^2/SO(2)` **is canonical.** Every
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
contrast with :math:`S^2/SO(2)`, whose syzygy ideal is *also* empty but
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
are two points, not curves. Contrast :math:`S^2/SO(2)`, where
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
   * - :math:`S^2/SO(2)`
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

ERR-080's level-1 half is a botched section of :math:`S^2/SO(2)`
------------------------------------------------------------------

The chart-versus-section question is not new with :math:`\sigma_y`. It
arises for :math:`SO(2)` too — the moment any consumer of a
1-dimensional rule needs a 3-D direction — and the tree has been
answering it, silently and wrongly, for as long as :ref:`ERR-080
<manifold-err-080>` has existed.

The realization :math:`[-1,1]` is the **chart's** codomain,
unambiguously: a section of :math:`S^2 \to S^2/SO(2)` is a half-meridian
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
   * - an honest :math:`\varphi = 0` half-meridian,
       :math:`\mu \mapsto (\sqrt{1-\mu^2},\, 0,\, \mu)`
     - on :math:`S^2` to :math:`0.0`; ``Sphere().contains`` → ``True``

⟹ **ERR-080's first link is not "a wrong tag". It is a section
fabricated where none was declared** — the realization is a chart, a
consumer needed a section, and zero-padding is what a codebase does
when the object it needs has no slot. With ``fundamental_domain`` in
the type, that padding has nowhere to live: an entry either declares a
section or declares that it has none, and :math:`S^2/SO(2)` honestly
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
     - :math:`S^2/SO(2)`
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
   * - `[M]` ``measure.support``
     - ``'[-1,1]'`` — the *realization's* name
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
       the row above, where only the *image* ships. Neither **map**
       does
   * - the pushforward measure :math:`\pi_*\mu`
     - ⛔ **not a slot**
     - the measure descends today via
       :meth:`DiscreteMeasure.quotient
       <orpheus.numerics.measure.DiscreteMeasure.quotient>`, which
       knows nothing about this catalogue

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
   P_ij = <grad p_i, grad p_j> with P >= 0; intersect with the ideal
   of S^2) and register it in orpheus/numerics/manifold.py's
   _ORBIT_CATALOGUE, or implement the derivation engine. Catalogued
   today (manifold CLASS / group): ['Sphere/SO2', 'Sphere/sigma_x',
   'Sphere/sigma_y', 'Sphere/sigma_z'].

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

The mint is not introducing a new idea. `[M]`
:attr:`AngularSymmetry.support
<orpheus.numerics.quadrature.registry.AngularSymmetry.support>` — which
predates it — already computes :math:`S^2/G^0` from the *spent* group
by catalogue lookup, in the string vocabulary, and already raises
:exc:`NotImplementedError` on an uncatalogued quotient with the same
shape of message. Two independent orbit-space lookups now exist:

.. list-table:: Two lookups of :math:`S^2/H`, re-measured 2026-08-31 after the two-slot ruling
   :header-rows: 1
   :widths: 18 34 48

   * - :math:`H`
     - ``AngularSymmetry(...).support``
     - ``SPHERE.quotient(H).realization.name``
   * - ``SO2``
     - ``'[-1,1]'``
     - ``'[-1,1]'`` — **they agree**
   * - ``Trivial``
     - ``'S^2'``
     - ``'S^2'`` — **they agree**
   * - ``sigma_y``
     - :exc:`NotImplementedError`
     - ``'D^2'`` — ⭐ the catalogue answers a row the registry
       **structurally cannot**; see below
   * - ``Oh``
     - :exc:`NotImplementedError`
     - :exc:`NotImplementedError`

Four readings, all useful.

**(i)** On the two overlapping rows the typed catalogue reproduces the
registry's derived answer exactly, which is the cheapest available
evidence that the type is a *re-typing* of an existing fact and not a
rival one. A committed row pins them together
(``test_the_trivial_answer_agrees_with_the_shipped_string_twin``) with
a *do not re-baseline* note: if it ever reddens, one of the two is
wrong about a quotient, and which one is a mathematical question, not
a test-maintenance one.

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
and the migration is what collapses it — but reading (iii) sharpens
what the collapse *is*. The registry's ``support`` is not a rival
catalogue; it is the **special case** :math:`H = G^0`, and the collapse
is therefore ``support = base.quotient(G⁰).realization.name`` rather
than a merge of two tables. That collapse is tracker 2.0b, and it is
listed as a seam rather than performed here.

.. warning::

   ⛔ **A live consequence of (iii), measured, and latent rather than
   broken today.** The registry's stage-0 admission gate is a **string
   comparison**: ``admits_domain`` is
   ``measure.support == self.support``. `[M]` the cylinder declares
   ``'S^2'`` while the shipped cylindrical rule carries
   ``folded_product(4,8).measure.support == 'S^2/sigma_y'``, so
   ``GEOMETRY_ANGULAR_SYMMETRY["cylinder"].admits_domain(...)`` is
   **False** — stage 0 would reject the tree's own fold.

   It bites nothing today for one reason only: `[M]` ``folded_product``
   is **not in** ``quadrature_registry`` (four specs ship —
   ``GaussLegendre1D``, ``LebedevSphere``, ``LevelSymmetricSN``,
   ``ProductQuadrature``), so the selector never presents it to stage 0.
   The day it is registered, this is the first thing that fires, and
   the fix is not to loosen the comparison: it is that a rule folded by
   a member of :math:`\Gamma` lives on :math:`S^2/\Gamma'` while the
   geometry declares :math:`S^2/G^0`, and those are two different
   quotients that the string vocabulary cannot tell apart. Recorded as
   a seam (:ref:`manifold-seams`).


.. _manifold-gotchas:

Gotchas
=======

.. _manifold-import-cycle:

The module imports nothing from ``numerics`` at runtime — on purpose
--------------------------------------------------------------------

:mod:`orpheus.numerics.manifold` references
:class:`~orpheus.numerics.symmetry.SubgroupOfO3` under
:data:`typing.TYPE_CHECKING` only, and calls no method on it. That is
**load-bearing, not tidiness**, and it is the kind of constraint that
is invisible until it is violated.

`[M]` by AST, 2026-08-31, with relative imports resolved:

- ``symmetry.py:98`` is ``from .measure import DiscreteMeasure`` — a
  **module-scope, runtime** edge ``symmetry → measure``;
- ``measure.py``'s own ``symmetry`` imports are **both deferred**: line
  91 under ``if TYPE_CHECKING:``, and line 1005 inside
  :meth:`DiscreteMeasure.quotient
  <orpheus.numerics.measure.DiscreteMeasure.quotient>` at function
  scope, whose comment already says why.

So once the migration makes ``measure`` import ``manifold``, a
module-scope ``manifold → symmetry`` edge closes the loop
``measure → manifold → symmetry → measure``.

`[M]` demonstrated on a throwaway three-module package with exactly
that edge topology (no production file touched):

.. code-block:: text

   ImportError: cannot import name 'DiscreteMeasure' from partially
   initialized module 'pkg.measure' (most likely due to a circular
   import)

and with the ``TYPE_CHECKING`` guard in place — the shipped shape — all
three modules import cleanly in either order.

⚠ **Two ways this is easy to get wrong, both measured.**

1. **An AST import census that filters on**
   ``node.module.startswith("orpheus")`` **silently drops every
   relative import**, because an ``ImportFrom`` with ``level > 0``
   carries an *unqualified* ``.module``. A census written that way
   reports ``symmetry → measure`` as **absent** and concludes there is
   no cycle — wrong, and in the reassuring direction. Resolve relative
   imports, and give the census a positive control.
2. **The cycle is not live today.** ``measure`` does not import
   ``manifold`` yet, so a module-scope ``manifold → symmetry`` edge
   added right now would break nothing, pass every test, and become
   fatal at the commit that lands 2.0b. The guard is prophylactic;
   removing it as "unnecessary" is the failure mode it exists to
   prevent.

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
   * - ``FunctionSpace.manifold``, and the derived ``L2[...]`` name
     - **2.0c.** Two sites build a level-2 name by interpolating a
       level-1 tag; one of them
       (``basis/indicator_basis.py:284``) hard-codes it and `[M]` is
       already **false** for the energy-grid basis
       (:ref:`manifold-string-algebra`).
   * - The two MAPS — the ``chart`` and the section — and the
       pushforward measure, as entry fields
     - **Phase 1.1 / 3.1.** `[M]` 7 of the derivation procedure's 9
       outputs are slots on
       :class:`~orpheus.numerics.manifold.Quotient` today. What ships
       of each map is only its *end*: the chart's **codomain**
       (``realization``) and the section's **image**
       (``fundamental_domain``); neither map itself is a slot, so
       nothing can currently *apply* one — which is why
       ``Quotient.contains`` must accept both languages rather than
       normalising to one (:ref:`manifold-two-coordinate-systems`).
       The pushforward measure is not a slot either: a measure
       descends via :meth:`DiscreteMeasure.quotient
       <orpheus.numerics.measure.DiscreteMeasure.quotient>`, which
       knows nothing about the catalogue
       (:ref:`manifold-engine-data-model`).
   * - The remaining catalogue entries
     - **Phase 1.1.** `[M]` **four keys** ship — ``(Sphere, "SO2")``
       and ``(Sphere, "sigma_x"/"sigma_y"/"sigma_z")`` — served by
       **two** procedures, since all three mirrors share one
       derivation that reads the axis off the group. The identity
       quotient is a fifth answer, derived rather than tabulated. The
       expected remainder covers :math:`\mathbb{Z}_2` antipodal,
       :math:`C_n` / :math:`D_n` about an axis, the :math:`O_h`
       sublattice for octant symmetry, :math:`SO(3)`, and
       :math:`SO(2)\times\mathbb{R}_z` for the 1-D cylinder.
       ⚠ Whoever adds a :math:`C_n` entry must leave
       ``fundamental_domain=None``: a closed sector is **not** a
       strict fundamental domain, and the ``dim`` gate cannot catch a
       wrong one (:ref:`manifold-chart-section-asymmetry`).
   * - Collapsing the twin lookup
     - **2.0b.** :attr:`AngularSymmetry.support
       <orpheus.numerics.quadrature.registry.AngularSymmetry.support>`
       performs its own orbit-space lookup in the string vocabulary;
       `[M]` the two agree on both overlapping rows, and the registry's
       is the **special case** :math:`H = G^0`, so the collapse is
       ``support = base.quotient(G⁰).realization.name`` rather than a
       merge of two tables (:ref:`manifold-twin-lookup`).
   * - The ``support`` tag's own vocabulary split
     - **2.0b, and independent of any one entry.** `[M]`
       ``gauss_legendre(8).measure.support`` is ``'[-1,1]'`` — the
       *realization's* name — while
       ``folded_product(4,8).measure.support`` is ``'S^2/sigma_y'`` —
       the *quotient's* name. One slot, two registers, and the
       registry's stage-0 gate compares it by string equality:
       `[M]` ``GEOMETRY_ANGULAR_SYMMETRY["cylinder"].admits_domain``
       on the shipped fold is **False**. Latent only because
       ``folded_product`` is not a registered spec
       (:ref:`manifold-twin-lookup`). A migration that retypes the slot
       and leaves the two registers disagreeing has fixed the type and
       not the drift.
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

The gates live in ``tests/numerics/test_manifold.py``: `[M]` **42 test
functions, 56 collected rows**, run under the canonical
``python -O -m pytest`` invocation. Two functions are parametrized —
one over the 13-member shipped-variant list, one over three bases — and
the two counts are given separately because they move for different
reasons: adding a *variant* moves the second and not the first.

⚠ The row count here is the one the generated V&V matrix reports for
this module, which is a second, independent reading of the same tree
(``docs/theory/verification/matrix.rst``, ``numerics/test_manifold``
row). `[M]` it read **44** before the two-slot ruling and **56** after.
An earlier version of this paragraph said *"30 test functions, 40
collected rows"*; both numbers were wrong when written — the module had
32 functions and 44 rows — which is why the count is now stated with
the instrument that produces it.

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

Six groups, and what each is for:

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

⭐ Two exhaustiveness gates are worth naming separately, because they
are what make "closed sum" a checkable claim rather than a description.
``test_every_variant_is_reachable_from_this_modules_list`` compares
``Manifold.__subclasses__()`` against the module's own exercised list,
so a member added to the module but not to the tests fails **there**;
``test_ambient_dimension_is_defined_for_every_variant`` walks every
shipped variant through the exhaustive ``match``. The benefit of a
closed sum is that an operation can be checked against every member,
and that benefit is only real if the member list is itself pinned.


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
  **orbifold** — which is what :math:`S^2/SO(2)` is and a quotient
  manifold is not (:ref:`manifold-singular-stratum`).
* Helgason, S. (1984). *Groups and Geometric Analysis*. Academic Press.
  Chapter IV (spherical functions on a Gelfand pair; the Funk–Hecke
  theorem as the :math:`(SO(3), SO(2))` instance), which is the
  literature behind :ref:`manifold-gelfand`.
* Hamermesh, M. (1962). *Group Theory and Its Application to Physical
  Problems*. Addison-Wesley. §2.5 (finite point groups) — the source
  :doc:`/theory/foundations/discrete_measures` already cites for the
  group side of the same construction.
