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
      role: "the point-set layer — the manifold M a measure is supported on and a basis function is defined over, its algebra (product, orbit space, membership), the invariant-theoretic derivation that produces an orbit space, and the three-level separation (manifold / fields on it / coefficients) that keeps a FunctionSpace from being mistaken for a domain"
      depends_on: []
      related: [discrete_measures, spaces, frame, spherical_harmonics]
      status: "the type is MINTED and gated (tracker 2.0a, 2026-08-31); it has ZERO production consumers — the migration off `Space = str` is 2.0b/2.0c/2.1 and has not landed"


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
   - ⛔ **The predicate is NOT yet wired.** The type ships and is
     gated; `[M]` it has **zero** production consumers, and
     ``Space = str`` is still live at ``measure.py:111``. ERR-080 is
     open. This page documents a **capability**, not a fix
     (:ref:`manifold-seams`).
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
     :math:`[-1,1]` (:eq:`manifold-s2-mod-so2`).
   - ⚠ **An orbit space is an ORBIFOLD, not a quotient manifold.**
     Where the action is not free, the image of the fixed-point set is
     a **singular stratum**. For :math:`S^2/SO(2)` that locus is
     :math:`\mu = \pm 1`, the poles — exactly where
     :math:`\det P` vanishes, and exactly where the curvilinear
     S\ :sub:`N` :math:`\alpha`-dome closes. A design that assumes a
     quotient is a smooth submersion is wrong there and only there
     (:ref:`manifold-singular-stratum`).
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
     falsifiable check, and `[M]` **6 of 8** of the procedure's outputs
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

Level 1 was ``Space = str`` (``measure.py:111``), with the module's own
comment on the alias set reading, verbatim:

   *"Common aliases used across the project. These are recommendations,
   not constraints; user code may pass arbitrary strings."*
   — ``measure.py:113-114``

That is an honest description of a slot with no semantics, and it is
still what ships. Three consequences follow from it, each measured; the
type answers all three, and the argument for minting it is these
sites, not the size of the migration.

.. _manifold-err-080:

(a) Nothing could be refused — the ERR-080 forgery
--------------------------------------------------

A 1-dimensional angular quadrature carries no azimuthal information.
:meth:`Quadrature.angular_frame
<orpheus.numerics.quadrature.directional.Quadrature.angular_frame>`
nonetheless builds its measure by ``column_stack``\ ing three
axis-cosine arrays — two of which are a zero *fallback*, not data — and
declares the result ``support=SPACE_SPHERE``. The rows are then
:math:`(\mu, 0, 0)` with
:math:`\lVert\Omega\rVert = |\mu| \neq 1`: points of :math:`[-1,1]`,
not of :math:`S^2`.

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

⚠ **And the last row is not a naming quibble — it already ships as a
falsehood.** ``measure.py:331`` at least *derives* its name, from a
``str``; ``basis/indicator_basis.py:284`` **hard-codes** it as
``f"L2[coarse_cells_R{self.ndim}]"``. And
:meth:`EnergyGrid.as_basis <orpheus.data.energy_grid.EnergyGrid.as_basis>`
builds an :class:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis`
over an **energy index** partition
(``edges = arange(n_groups + 1) - 0.5``), so that basis names its own
coefficient space after a *spatial* manifold it has nothing to do
with. `[M]` reproduced 2026-08-31 on a two-group grid:

.. code-block:: python

   >>> eg = EnergyGrid(edges=np.array([1e6, 1.0, 1e-5]))   # 2 groups
   >>> eg.as_basis().space.name
   'L2[coarse_cells_R1]'

A :class:`~orpheus.numerics.space.FunctionSpace` that carried its
manifold would make both names derived and this one **unspellable** —
which is the level-2 half of the same repair, tracked at
:ref:`manifold-seams`.

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

Nine variants, in three families plus two recursive constructors. Two
dimensions travel with each member and are easy to conflate: the
**topological** dimension ``dim`` (what the manifold *is*) and the
**ambient** coordinate count (how many columns ``contains`` consumes).
They differ for every curved member.

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
   * - ``Quotient(M, H, …)``
     - :math:`\dim` of the realization
     - realization's
     - ``M/H``
     - decided in the **realization's** coordinates

`[M]` the names reproduce the retired ``SPACE_*`` string tags
**verbatim** — ``S^2``, ``S^1``, ``[-1,1]``, ``[0,1]``, ``[0,inf)``,
``R``, ``energy``, ``spatial_R1``, ``index(angular)`` — so the
migration off ``Space = str`` cannot silently re-word an error message
or an ``L2[...]`` space name that a test pins. That is a deliberate
constraint on the type, gated by
``tests/numerics/test_manifold.py::TestTypeLaws::test_the_names_reproduce_the_retired_string_tags``.

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
runs it in full on the one pair the tree ships, and draws the three
consequences that pay for it downstream.

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

This is the one entry the tree catalogues, and it is the entry every
1-dimensional angular discretisation is secretly using. Every line
below was **re-derived and re-run** in this session, independently of
the catalogue, and then compared against it.

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
     - ``(-1.0, 1.0)`` / ``False``
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
:math:`\det P`; the chart; the pushforward measure; the stratum where
:math:`\det P` vanishes. **Those are the entry's fields.** An engine
then ships by *computing* them instead of reading them — a development,
with no new vocabulary and no seam.

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
     - *derived* from ``det_gram``, and gated as such
   * - provenance
     - ``derived_by``
     - ``"hand"`` / ``"engine"``
   * - the chart :math:`M/H \to N`
     - ⛔ **not a slot**
     - only its **codomain** ships, as ``realization``
   * - the pushforward measure :math:`\pi_*\mu`
     - ⛔ **not a slot**
     - the measure descends today via
       :meth:`DiscreteMeasure.quotient
       <orpheus.numerics.measure.DiscreteMeasure.quotient>`, which
       knows nothing about this catalogue

`[M]` **6 of 8** — by ``dataclasses.fields``, ``Quotient`` declares
``base``, ``by``, ``realization``, ``generators``, ``syzygy``,
``gram``, ``det_gram``, ``derived_by``, ``singular_stratum``. So the
seed is real for six of the procedure's outputs and **incomplete for
two**, and those two are named as seams rather than left to be
discovered (:ref:`manifold-seams`). Stating the fraction is the point:
a ruling whose compliance is claimed but not counted is not checkable.

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
   today (manifold CLASS / group): ['Sphere/SO2'].

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

Concretely, for the one shipped entry, the assertions are on the
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
   assert s2_mod_so2.singular_stratum == (-1.0, 1.0)

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

.. list-table:: Two lookups of :math:`S^2/H`, measured 2026-08-31
   :header-rows: 1
   :widths: 20 40 40

   * - Spent :math:`H`
     - ``AngularSymmetry(...).support``
     - ``SPHERE.quotient(H).realization.name``
   * - ``SO2``
     - ``'[-1,1]'``
     - ``'[-1,1]'`` — **they agree**
   * - ``Trivial``
     - ``'S^2'``
     - ``NotImplementedError`` — ⛔ the catalogue lacks the identity
       quotient
   * - ``Oh``
     - ``NotImplementedError``
     - ``NotImplementedError``

Three readings, all useful. **(i)** On the one overlapping row the
typed catalogue reproduces the registry's derived answer exactly, which
is the cheapest available evidence that the type is a *re-typing* of an
existing fact and not a rival one. **(ii)** The ``Trivial`` row is a
real gap: :math:`S^2/\{e\} = S^2` is legal, trivially derivable
(:math:`p_i = x_i`, :math:`P = I`, :math:`\det P = 1`, no stratum) and
absent. **(iii)** Two lookups of one fact is a Pattern-2 twin by
construction, and the migration is what collapses it — the registry's
``support`` should *be* the catalogue's answer rather than a parallel
``if``/``if``/``raise``. That collapse is tracker 2.0b, and it is
listed as a seam rather than performed here.


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

A :class:`~orpheus.numerics.manifold.Quotient`'s ambient count is its
**realization's**, because membership in an orbit space is decided in
the coordinates the chart lands in — :math:`\mu`, not
:math:`(x, y, z)`.

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
   * - ⛔ **Any production consumer at all**
     - **#429 tracker 2.0b.** `[M]` 2026-08-31: the only importers of
       :mod:`orpheus.numerics.manifold` are its own test module.
       ``Space = str`` is still live at ``measure.py:111`` with its six
       ``SPACE_*`` aliases, ``DiscreteMeasure.support`` is still a
       ``str``, and **ERR-080 is still open** — held by an
       ``xfail(strict=True)`` gate at
       ``tests/sn/solve/test_pl_order_does_not_move_the_infinite_medium_flux.py``.
       Every refusal described on this page is a capability that
       *would* fire once the slot is retyped.
   * - ``Basis.domain``
     - **2.1.** No :class:`~orpheus.numerics.basis.base.Basis` can
       state the manifold its functions consume, which is why the
       ERR-080 pairing has nothing to check
       (:ref:`manifold-three-levels`).
       :class:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis`
       is expected to take it as a **constructor field** rather than
       derive it, and the reason is measurable: `[M]` by AST, of 18
       ``IndicatorBasis(...)`` construction sites tree-wide **4 are in**
       ``orpheus/``, and those four partition **three different
       manifold families** — a finite index set
       (``frame.py:755``, paired with ``support=f"index({axis_label})"``
       three lines below), :math:`\mathbb{R}^d` at two ranks
       (``geometry/mesh.py:444`` and ``:753``), and the energy counting
       set (``data/energy_grid.py:220``). Any value the class *derived*
       from its own fields would hard-code one of the three.
   * - ``FunctionSpace.manifold``, and the derived ``L2[...]`` name
     - **2.0c.** Two sites build a level-2 name by interpolating a
       level-1 tag; one of them
       (``basis/indicator_basis.py:284``) hard-codes it and `[M]` is
       already **false** for the energy-grid basis
       (:ref:`manifold-string-algebra`).
   * - The ``chart`` and the pushforward measure as entry fields
     - **Phase 1.1 / 3.1.** `[M]` 6 of the derivation procedure's 8
       outputs are slots on
       :class:`~orpheus.numerics.manifold.Quotient` today; the chart
       ships only as its **codomain** (``realization``), and the
       pushforward measure not at all — a measure descends via
       :meth:`DiscreteMeasure.quotient
       <orpheus.numerics.measure.DiscreteMeasure.quotient>`, which
       knows nothing about the catalogue
       (:ref:`manifold-engine-data-model`).
   * - The other ~11 catalogue entries
     - **Phase 1.1.** One entry ships: ``(Sphere, "SO2")``. The
       expected set covers :math:`\mathbb{Z}_2` antipodal,
       :math:`C_n` / :math:`D_n` about an axis, the :math:`O_h`
       sublattice for octant symmetry, :math:`SO(3)`, and
       :math:`SO(2)\times\mathbb{R}_z` for the 1-D cylinder. The
       identity quotient :math:`S^2/\{e\}` is a measured gap
       (:ref:`manifold-twin-lookup`).
   * - Collapsing the twin lookup
     - **2.0b.** :attr:`AngularSymmetry.support
       <orpheus.numerics.quadrature.registry.AngularSymmetry.support>`
       performs its own two-row orbit-space lookup in the string
       vocabulary; the two agree today `[M]` on the one overlapping
       row, and one of them should become a consumer of the other
       (:ref:`manifold-twin-lookup`).
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

The gates live in ``tests/numerics/test_manifold.py``: `[M]` **30 test
functions, 40 collected rows**, all passing under the canonical
``python -O -m pytest`` invocation.

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

Four groups, and what each is for:

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
     - `[M]` **8 of the 9** membership tests assert a positive
       verdict beside their negative one (``vv-principles`` #11 — a
       contract predicate tested only against a broken instance
       validates the *raising*, not the *claim*); the ninth is
       ``test_a_wrong_ambient_dimension_is_a_typed_refusal``, which
       asserts a **raise** rather than a verdict and so has no positive
       leg to carry (:ref:`manifold-gotcha-shape-vs-verdict`). The
       load-bearing row is the ERR-080 forgery:
       the negative leg refuses :math:`(\mu,0,0)`, the positive leg
       admits the same nodes normalised, and a third assertion places
       them on :math:`[-1,1]` where they belong.
   * - The recorded derivation
     - The symbolic regression tests of
       :ref:`manifold-engine-seed` — the :math:`P` matrix, the
       determinant, and the stratum **solved for** rather than
       compared to a literal.

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
     - **The manifold becomes an OBJECT.** Level 1 of the three-level
       stack stops being ``Space = str`` — an opaque tag whose own
       comment called the entries *"recommendations, not
       constraints"* — and becomes a closed sum with an algebra:
       ``dim`` / ``name`` / ``contains`` / ``__mul__`` total on the
       base, and the invariant-theoretic derivation fields on
       :class:`~orpheus.numerics.manifold.Quotient` alone. Nine
       variants, two of them consolidating families the string
       vocabulary had spelled several ways
       (:ref:`manifold-members`). The three morphisms the mint names
       were already running as string interpolation at
       ``measure.py:588`` / ``:1022`` / ``:802``
       (:ref:`manifold-string-algebra`). One catalogue entry ships,
       ``S^2/SO(2) = [-1,1]``, derived by the Procesi–Schwarz
       procedure and carrying its own symbolic regression tests
       (:ref:`manifold-s2-so2`) — which, because the entry's fields
       *are* the procedure's outputs, are the deferred derivation
       engine's acceptance suite, written before it
       (:ref:`manifold-engine-seed`). ⛔ **Type only: zero production
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
