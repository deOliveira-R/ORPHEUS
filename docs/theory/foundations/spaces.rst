.. _function-spaces:

==========================================================
Function Spaces: Axes, Measure, and the Collapse Doctrine
==========================================================

.. contents:: Contents
   :local:
   :depth: 2


.. Machine header — the ``nexus-meta`` schema for this page (PROVISIONAL).
.. Seeded at campaign-1 CS1 (2026-08-20) as ``field_algebra``'s sibling:
.. that page owns the ELEMENT algebra, this one owns the SPACES those
.. elements live in. The schema is provisional pending a full re-audit of
.. the corpus.

.. dropdown:: Machine header — ``nexus-meta`` schema (PROVISIONAL)
   :color: muted

   .. code-block:: yaml

      module: numerics
      concept: function_spaces
      role: "the space layer — a function space as the ordered product of its AXES (index shape, factor measure, basis kind, generator identity), the counting-measure theorem on the energy axis, and the collapse doctrine that decides which axes survive a degeneracy and why"
      depends_on: [field_algebra, frame]
      related: [discrete_measures, operator_algebra, operator_adjoint]
      status: "seeded at campaign-1 CS1 (the Energy axis); the Spatial / Quadrature / Harmonic axes are CS2"


This page develops the **space layer**: what a discrete function space
*is*, what it is made of, and what happens to it when a physical
degeneracy collapses one of its factors. It is the companion to
:doc:`/theory/foundations/field_algebra`, which types the **elements**
(a flux is a point of the positive cone :math:`K \subset V`); this page
types :math:`V` itself.

The organizing claim is one sentence: **a function space is the ordered
product of its axes, and an axis carries exactly four things — an index
shape, a factor measure, a basis kind, and the identity of the generator
that produced it.** Everything else on this page follows from taking
that seriously: why the energy metric is the identity as a *theorem*
rather than a default, why the homogeneous solver's spatial factor
survives its own collapse while the angular factor of a scalar space
does not, and why "two copies of :math:`\mathbb{R}^n` with different
inner products are the same space" is a claim this corpus has
**overturned**.

.. note::

   **Three unrelated things in this corpus are called an "axis". Keep
   them apart.**

   - A **space-factor axis** (:class:`~orpheus.numerics.axis.Axis`,
     :class:`~orpheus.numerics.axis.EnergyAxis`) — *this page's*
     subject. One tensor factor of a function space; it carries a
     measure and a basis kind, and it knows nothing about geometry.
   - A **geometric axis** (``Axis1D`` / ``AxisMesh`` /
     ``RadialAxisMesh`` in :mod:`orpheus.transport.mesh.axis`) — one
     coordinate DIRECTION of a structured mesh, carrying edges, face
     labels and a coordinate system. The name is a known misnomer
     (it declares per-axis geometry and creates no mesh); the rename is
     tracked as issue **#393**, and the two vocabularies are being kept
     deliberately distinct in the meantime. Where both are in play the
     code spells the space-factor one out: ``from orpheus.numerics.axis
     import Axis as SpaceFactorAxis``
     (:mod:`orpheus.transport.mesh.material_mesh`).
   - A **symmetry axis** — the invariant line of a rotation or mirror,
     as used throughout :doc:`/theory/foundations/discrete_measures`
     and the quadrature machinery.

   ⚠ The first two collide on an attribute NAME, on one object. `[M]`
   for a carrier ``mm``, ``mm.axes`` is the GEOMETRIC tuple
   (``(AxisMesh,)`` — the mesh's coordinate directions) while
   ``mm.bulk_space.axes`` is the SPACE-FACTOR tuple
   (``(EnergyAxis, Axis)``). They are different attributes of different
   objects that happen to share a spelling; neither is derived from the
   other.

.. note::

   **Two symbols on this page look alike.** :math:`V` (italic) is the
   **function space** — the object this page is about — and
   :math:`V_{\rm cell}` is a **cell volume**, a weight on the spatial
   axis. They meet constantly here (the whole point of clause 1 is that
   a spatial factor's weight is what distinguishes two spaces), so the
   volume is always written out in full. This follows
   :doc:`/theory/foundations/field_algebra`, which draws the same
   distinction for the same reason.

.. admonition:: Key Facts (the space layer)
   :class: tip

   - **A space IS its axis tuple.** Shape is the concatenation of the
     axes' shapes; the metric is the tensor product of the factor
     measures, stored PER AXIS and never densified
     (:eq:`spaces-axis-product`;
     :meth:`FunctionSpace.of_axes
     <orpheus.numerics.space.FunctionSpace.of_axes>`).
   - **An axis is (index shape, factor measure, basis kind, generator
     identity).** ``weights=None`` **is** the counting measure —
     deliberately and always; an axis has no "unbound" state, so the
     legacy two-state ambiguity of ``inner_product_weights`` cannot
     arise on this type.
   - **The two one-line discriminators of the collapse doctrine.**
     *(i)* **Can the admissible fields be integrated over the collapsed
     domain?** — symmetry-forced constancy on infinite-measure orbits
     means no, so the collapse must NORMALIZE. *(ii)* **Is the
     surviving convention consulted INSIDE the family, or only at
     re-embedding?** — inside ⟹ the axis persists; only at re-embedding
     ⟹ the convention lives on the arrow and the axis DROPS.
   - **A collapse that DROPS an axis is realized by two typed arrows,
     never one.** :meth:`FunctionSpace.retraction
     <orpheus.numerics.space.FunctionSpace.retraction>` is fiber
     integration :math:`R = \pi_*`; :meth:`~orpheus.numerics.space.FunctionSpace.section`
     is its right inverse :math:`E` with :math:`R\circ E = \mathrm{id}`.
     They differ by exactly the axis's total mass
     (:math:`R^\dagger = \Sigma w\,E`, the *plain* broadcast versus the
     normalized one), and they carry different TYPES so that scalar
     cannot be dropped at a call site. Both are induced by a rank-one
     indicator frame and memoized on the space
     (:ref:`spaces-collapse-pair`).
   - **EnergyGrid is a 1-D mesh in energy** — groups are its cells,
     group boundaries are its faces, and condensation is the
     mesh-overlap map
     (:meth:`EnergyGrid.overlap_to
     <orpheus.data.energy_grid.EnergyGrid.overlap_to>`). The one-group
     member is the one-CELL energy mesh; it keeps its edges because they
     define :math:`\bar\sigma`.
   - **The quotient point records the DENSITY CONVENTION.** The
     homogeneous solver's spatial factor is not absent and not
     measureless: it is a one-point axis whose unit weight *is* the
     normalized "per unit volume" convention, and the pairing consumes
     it. A genuine one-cell mesh with :math:`V_{\rm cell} \neq 1` is a
     **different** space.
   - **The energy metric is the identity as a THEOREM.** Multigroup
     flux components are group INTEGRALS (covariant, extensive) and
     cross sections are flux-weighted group AVERAGES (contravariant,
     intensive), so :math:`\int\sigma\varphi\,\mathrm{d}E =
     \sum_g \sigma_g\varphi_g` exactly and no group widths appear
     (:eq:`energy-condensation-counting-measure`). Consequence:
     :math:`V \cong V^*` isometrically along energy and the adjoint
     there is the plain transpose. A weighted
     :class:`~orpheus.numerics.axis.EnergyAxis` is REFUSED at
     construction.
   - **NODAL vs MODAL is the coordinate-cone question.** A nodal factor
     has one (components are cell/point values); a modal factor does not
     (components are expansion coefficients). ``has_coordinate_cone`` is
     three-valued — ``True`` / ``False`` / ``None`` for legacy spaces —
     and :meth:`Field.cone_violations
     <orpheus.numerics.field.Field.cone_violations>` consults it,
     REFUSING on ``False`` rather than manufacturing violations out of a
     basis choice.
   - ⚠ **Identity is a MIGRATION statement, both halves true.** The
     chartered doctrine is *identity = the axes' structural content*, so
     **metric differences imply space differences**. The current
     realization still compares ``(name, shape)``. The bridge is
     :meth:`of_axes <orpheus.numerics.space.FunctionSpace.of_axes>`,
     which derives the name deterministically and injectively from the
     axes — so for axis-built spaces the chartered identity already
     flows through the nominal comparison TODAY. The S3 flip retires the
     bridge.


.. _spaces-the-axis:

The axis: four slots, and what each one decides
===============================================

An **axis** is one tensor factor of a function space — the value object
recording *(index shape, factor measure, basis kind, generator
identity)*. It is the unit the composition machinery reasons about:
partitions, collapses, frames and (later) :math:`\oplus`-lifts act **per
axis**, never on an anonymous position of a monolithic shape tuple.

.. list-table:: The four slots
   :header-rows: 1
   :widths: 16 30 54

   * - Slot
     - Type
     - What it decides
   * - ``shape``
     - ``tuple[int, ...]``, rank :math:`\ge 1`
     - The factor's index set. Rank :math:`> 1` is admissible and
       deliberate: a spherical-harmonic axis is :math:`(L+1, 2L+1)`, and
       a rank-:math:`d` spatial axis is a legal CS2 design choice. Rank
       0 is refused — a factor with no index set has nothing to measure.
   * - ``weights``
     - ``NDArray | None`` over exactly ``shape``
     - The **factor measure**. ``None`` IS the counting measure. This is
       the slot the pairing consumes, and (through identity) the slot
       that makes a quotient point and a one-cell mesh different spaces.
   * - ``kind``
     - :class:`~orpheus.numerics.axis.BasisKind`
     - ``NODAL`` ⟹ the factor carries a **coordinate cone**;
       ``MODAL`` ⟹ it does not. Keyword-only with **no default** — the
       basis character is physics and must be spelled at every mint.
   * - identity
     - structural, **per subclass**
     - :math:`(\text{type}, \text{label}, \text{shape}, \text{kind},
       \text{weights bytes})` plus each subclass's own data. Two axes
       differing only in measure are DIFFERENT axes.

Composition, and the metric that comes with it
----------------------------------------------

A space is the ordered product of its axes. Two things are determined by
that sentence alone, and both are what
:meth:`FunctionSpace.of_axes
<orpheus.numerics.space.FunctionSpace.of_axes>` implements:

.. math::
   :label: spaces-axis-product

   V \;=\; \bigotimes_{a} V_a,
   \qquad
   \mathrm{shape}(V) \;=\; \mathrm{shape}(V_1) \frown \cdots \frown
   \mathrm{shape}(V_n),
   \qquad
   G_V \;=\; \bigotimes_{a} G_a,
   \quad
   G_a = \operatorname{diag}(w_a) \ \text{ or } \ I,

.. (vv-status rationale) Structural/representational identity: it STATES
   what ``FunctionSpace.of_axes`` constructs (shape = concatenation of
   the factors' shapes; metric = tensor product of the per-axis factor
   measures, applied factor-by-factor and never materialized). Not a
   solver claim — no flux, no eigenvalue, no discretization error. The
   verifiable content is the CS1 foundation battery
   (``tests/numerics/test_space_of_axes.py``: shape concatenation, the
   per-axis metric against an independently-built dense reference, the
   no-densification proof, and the derived name's determinism across
   processes) plus ``tests/numerics/test_axis.py`` for the factor laws.
.. vv-status: spaces-axis-product documented

where :math:`\frown` is tuple concatenation and :math:`w_a` is axis
:math:`a`'s weight array. The second half is the operational content:
the measure lives **per axis**, so composing two 2000-point weighted
axes stores :math:`2 \times 2000` weights and never the
:math:`4{,}000{,}000`-entry outer product. The metric is *applied*
factor by factor — each axis's weights broadcast into their own slot of
the element, with leading and trailing ones for the neighbouring
factors — so an interior weighted axis works exactly like a leading one.
That is a position the legacy prefix-broadcast convention could not
reach.

.. warning::

   **The legacy twin is still live, and CS2 retires it.** ``V * W`` (the
   ``*`` dunder →
   :class:`~orpheus.numerics.space.TensorProductSpace`) is the PRE-axis
   composition mechanism: it DENSIFIES the metric into an outer-product
   ``inner_product_weights`` and derives its name by joining the
   factors' names. CS1 keeps it — it threads the ``axes`` record when
   both sides carry one, and bridges axis-borne measures into its dense
   weights on mixed products, so no value is ever lost — and CS2
   collapses the live mints onto axis concatenation and retires the
   densifier. Until then the rule is: **new axis-aware code composes
   with** ``of_axes``; ``*`` **is the legacy surface.**

Canonical storage: one measure, one spelling, one identity
-----------------------------------------------------------

Because the measure is part of the identity, two spellings of the *same*
measure would be two identities of the same space — the exact twin the
axis layer exists to prevent. Three construction rulings close it
(2026-08-20), and each is enforced in
:meth:`Axis.__post_init__ <orpheus.numerics.axis.Axis>`:

#. **All-ones weights collapse to** ``None`` **at construction.** The
   counting measure has ONE spelling. Without this,
   ``weights=None`` and ``weights=np.ones(shape)`` would be the same
   measure with unequal identities. `[M]` the two axes compare equal and
   hash equal after canonicalization.
#. **Weights are canonicalized as** ``w + 0.0`` **and stored
   read-only.** :math:`-0.0` and :math:`+0.0` are one measure and must
   be one byte pattern — the identity key reads ``weights.tobytes()``,
   so an un-normalized :math:`-0.0` would mint a second name for one
   measure. The addition also forces a fresh allocation, so mutating the
   caller's array can never move an axis's hash after it has been used
   as a dictionary key.
#. **Non-finite weights are REFUSED.** A factor measure has finite
   weights.

There is deliberately **no non-negativity guard**. CS2's quadrature axes
legally carry signed weights (level-symmetric families with negative
weights are real, shipped objects), and the axis is the wrong layer to
outlaw them. `[M]` ``Axis("x", (2,), weights=[-1.0, 2.0], kind=NODAL)``
constructs.

Identity is structural, and it is per subclass
-----------------------------------------------

Equality and hashing read the structural content — never object
identity, and never a subset of the fields. Two consequences are
load-bearing and both are gated:

- **An** :class:`~orpheus.numerics.axis.EnergyAxis` **never equals a
  generic** :class:`~orpheus.numerics.axis.Axis` **carrying the same
  field tuple.** The identity of an axis is *what kind of generator
  produced this factor*, not a bag of fields. `[M]` the comparison is
  ``False``.
- **A synthetic axis never equals a** ``from_grid`` **axis of the same**
  ``ng``. Same index set, no partition data — a different axis. `[M]`
  ``EnergyAxis.synthetic(2) == EnergyAxis.from_grid(grid_2g)`` is
  ``False``.

.. _spaces-identity-bridge:

The derived name is the identity bridge (the S3 seam)
------------------------------------------------------

Space identity is nominally ``(name, shape)`` and stays that way until
the S3 flip. The bridge that makes the chartered doctrine *already true
for axis-built spaces* is the name derivation: ``of_axes`` computes the
name deterministically and injectively from the axes' structural content
— a length-prefixed, type-tagged content digest, never Python's
``hash()``, so it is stable across processes. Different axis tuples mint
different names, hence different spaces, **today**.

The mechanism is visible in one measurement. The homogeneous carrier's
quotient point and a genuine one-cell mesh with :math:`V_{\rm cell} = 2` have the
same shape :math:`(2, 1)` and the same readable prefix, and are
**unequal** spaces:

.. code-block:: text

   quotient point (volumes = [1.0])  ->  energy(2,)*spatial(1,)#<digest A>
   one-cell mesh  (volumes = [2.0])  ->  energy(2,)*spatial(1,)#<digest B>
   digest A != digest B   ->   the two spaces compare UNEQUAL

That is the "metric differences imply space differences" doctrine
flowing through the nominal identity with no flag day. The readable
prefix is for humans; the digest is the identity. S3 retires the bridge
by comparing the axes tuple directly.

.. note::

   **What was overturned.** Until this campaign the
   :class:`~orpheus.numerics.space.FunctionSpace` docstring taught that
   *two copies of* :math:`\mathbb{R}^n` *are "the same" space regardless
   of which inner product is installed*. That is false in the only sense
   the operator algebra cares about: the metric defines ``.H``, so two
   spaces differing only in metric are spaces where the same symbol
   denotes **different operators**, and composing across them must be
   refused rather than silently accepted. The docstring now carries both
   halves — the chartered doctrine and the current nominal realization —
   because stating only the target would lie forward and stating only
   the present would lie backward.

.. _spaces-nodal-modal:

NODAL and MODAL: which factors have a coordinate cone
------------------------------------------------------

:doc:`/theory/foundations/field_algebra` establishes that flux lives in
the positive cone :math:`K \subset V` and that cone membership is an
**element predicate**. The space layer supplies the missing half: *on
which spaces is that predicate even meaningful?*

- A **NODAL** factor has components that are point or cell VALUES (an
  indicator-like basis). Pointwise nonnegativity of the coefficients
  **is** nonnegativity of the function, so the coordinate cone
  :math:`K = \{x \ge 0\}` is the physical positive cone. The
  discrete-ordinates axis, the spatial cell axis and the energy group
  axis are all nodal.
- A **MODAL** factor has components that are expansion COEFFICIENTS. A
  positive function may have negative coefficients, so a per-component
  sign test answers a question about the *basis*, not about the
  function. The spherical-harmonic moment axis is modal.

This dichotomy is the ray-effect / negative-source dichotomy seen from
the algebra side: positivity is native to the quadrature axis,
rotational equivariance is native to the harmonic axis, and no angular
basis has both.

:attr:`FunctionSpace.has_coordinate_cone
<orpheus.numerics.space.FunctionSpace.has_coordinate_cone>` is
**three-valued**, and the third value is the point:

.. list-table::
   :header-rows: 1
   :widths: 14 30 56

   * - Value
     - When
     - What consumers do
   * - ``True``
     - axis-built, ALL factors ``NODAL``
     - Answer. The sign test is meaningful; the arithmetic is unchanged
       from the legacy path.
   * - ``False``
     - axis-built, ANY factor ``MODAL``
     - **Refuse**, with a typed error naming the space and the reason.
       Answering would manufacture violations out of a basis choice.
   * - ``None``
     - ``axes is None`` — legacy, not migrated
     - Pre-CS1 behavior, unchanged. The question cannot be answered
       structurally, and collapsing ``None`` into ``False`` would fire
       the refusal on every legacy space in the tree.

`[M]` **the refusal has no production witness yet, deliberately.** The
only axis mint in ``orpheus/`` today is
:attr:`MaterialMesh.bulk_space
<orpheus.transport.mesh.material_mesh.MaterialMesh.bulk_space>`, whose
factors are both ``NODAL``. Since :meth:`of_axes
<orpheus.numerics.space.FunctionSpace.of_axes>` is the only ROOT
producer of an ``axes`` record (``*`` and
:meth:`dual <orpheus.numerics.space.FunctionSpace.dual>` merely thread
one through, so both need an axis-built ancestor), it follows that
every harmonic-moment space in the tree is still legacy
(``axes is None``) and therefore takes the ``None`` arm.
The ``False`` arm's witness is test-constructed
(``tests/numerics/test_field.py``, gates E1/E2 — the refusal and the
same values answered on an all-nodal space, the positive-and-negative
pair ``vv-principles`` anti-pattern #11 requires). The arm becomes
production-reachable when CS2 mints the harmonic axis.


.. _spaces-counting-measure-theorem:

The counting-measure theorem on the energy axis
===============================================

The energy factor is the one axis whose measure is not a modelling
choice. It is forced — by the multigroup convention itself — to be the
counting measure, and *that* is why the energy metric is the identity.
This section derives it, states what the derivation buys at the space
layer, and says where the claim's single source of truth lives.

.. important::

   **Single source of truth.** The counting-measure claim
   (:math:`w_g = 1`, not :math:`w_g = \Delta u_g`) is stated, derived
   from Hébert's continuous formulation, and gated on
   :doc:`/theory/foundations/frame` — see
   :eq:`energy-condensation-counting-measure` inside
   :ref:`sn-energy-condensation`, whose verifiable content is the
   rate-preservation gate (a :math:`\Delta u_g` weight breaks it). This
   page does **not** restate that equation. What it owns is the claim's
   *space-layer consequence* — the energy metric is
   :math:`G_E = I`, hence :math:`V \cong V^*` isometrically along energy
   and the adjoint there is the plain transpose — and the fact that the
   consequence is now enforced at construction. Edited there, consumed
   here.

Covariant and contravariant: why no width appears
--------------------------------------------------

Write the continuous pairing that every reaction rate is:

.. math::

   r \;=\; \int_0^\infty \sigma(E)\,\varphi(E)\,\mathrm{d}E .

Discretizing splits into two *different* kinds of object, and the whole
theorem is that ORPHEUS discretizes each one in its own natural
variance:

- **The flux components are group INTEGRALS** — covariant, extensive
  quantities:

  .. math::

     \varphi_g \;=\; \int_{E_{g+1}}^{E_g} \varphi(E)\,\mathrm{d}E .

  The bin width is *inside* :math:`\varphi_g`; the stored number is
  "eV-free" and is a member of :math:`V`.

- **The cross sections are flux-weighted group AVERAGES** —
  contravariant, intensive quantities:

  .. math::

     \sigma_g \;=\;
     \frac{\displaystyle\int_{E_{g+1}}^{E_g} \sigma(E)\,\varphi(E)\,
     \mathrm{d}E}
          {\displaystyle\int_{E_{g+1}}^{E_g} \varphi(E)\,\mathrm{d}E} .

  This is a co-vector: it is a functional *on* fluxes, and it lives in
  :math:`V^*`.

Substituting the second definition into a group-by-group sum and
collapsing the denominator against :math:`\varphi_g` gives the pairing
back, **exactly** and with no measure factor left over:

.. math::

   \sum_g \sigma_g\,\varphi_g
   \;=\;
   \sum_g \frac{\int_g \sigma\varphi\,\mathrm{d}E}{\int_g
   \varphi\,\mathrm{d}E}\;\int_g \varphi\,\mathrm{d}E
   \;=\;
   \sum_g \int_g \sigma\varphi\,\mathrm{d}E
   \;=\;
   \int_0^\infty \sigma(E)\varphi(E)\,\mathrm{d}E
   \;=\; r .

That identity is the theorem. The two variances were *chosen* so that
the discrete pairing is the continuous integral with **weight one**;
introducing a lethargy width :math:`\Delta u_g` would double-count the
width and break rate preservation
(:eq:`energy-condensation-counting-measure`). Lethargy is the node
*coordinate*, never the *weight*.

.. note::

   **This is the spatial axis's exact opposite, and the contrast is the
   fastest way to remember which is which.** A spatial flux
   :math:`\phi_i` **is a density** — it was never integrated over the
   cell — so pairing it against a cross section requires the geometric
   volume measure :math:`V_i`
   (:eq:`sn-homogenization-fine-rate`). An energy flux
   :math:`\varphi_g` **is already an integral**, so pairing it requires
   nothing. Same equation, opposite measures, and the difference is
   entirely in which variance the discretization chose. Getting this
   backwards is the classical missing-width / double-counted-width bug
   in group-constant generation.

What the theorem buys at the space layer
-----------------------------------------

Three statements, each a direct consequence:

#. **The energy metric is** :math:`G_E = I_{n_g}`. Not "defaults to",
   not "is conventionally taken as" — *is*, by the argument above. In
   axis vocabulary: the energy factor's ``weights`` is ``None``, which
   IS the counting measure.
#. :math:`V \cong V^*` **isometrically along energy.** The Riesz map on
   that factor is the identity, so the distinction between a flux
   (:math:`V`) and a cross section (:math:`V^*`) is carried entirely by
   the *role*, never by a metric that would have to be applied to move
   between them.
#. **The adjoint along energy is the plain transpose.** Since
   :math:`A^\dagger = G^{-1} A^{\mathsf T} G`
   (:doc:`/theory/foundations/operator_adjoint`) and :math:`G_E = I`,
   the Hilbert adjoint and the Euclidean transpose coincide *on that
   factor*. This is why the energy-only operators of the homogeneous
   solver could be bound to a real space with **no value motion at
   all** — see :ref:`spaces-development-history`.

Construction enforces it
-------------------------

The theorem is not left as prose for a future contributor to violate.
:class:`~orpheus.numerics.axis.EnergyAxis` **refuses weights at
construction**, with a message that states the reason; a deliberately
non-physical weighted toy must use a generic
:class:`~orpheus.numerics.axis.Axis`, which is exactly what the
``.H``-sensitivity control in the CS1 battery does. Both constructors —
:meth:`EnergyAxis.from_grid <orpheus.numerics.axis.EnergyAxis.from_grid>`
and :meth:`EnergyAxis.synthetic
<orpheus.numerics.axis.EnergyAxis.synthetic>` — mint ``weights=None``.
The axis is also ``NODAL`` and rank-1 by construction, both refused
otherwise.

.. _spaces-energy-grid-is-a-mesh:

EnergyGrid is a 1-D mesh in energy
-----------------------------------

The reading that makes the energy axis a *member of the same family* as
the spatial one, rather than a special case:

.. list-table::
   :header-rows: 1
   :widths: 26 34 40

   * - Mesh concept
     - In energy
     - Carrier
   * - cells
     - the groups
     - ``shape = (ng,)``
   * - faces
     - the group boundaries
     - ``edges`` (``ng + 1`` values)
   * - the mesh-overlap map
     - condensation (fine → coarse)
     - :meth:`EnergyGrid.overlap_to
       <orpheus.data.energy_grid.EnergyGrid.overlap_to>`
   * - a one-cell mesh
     - the one-group member
     - :math:`\bar\sigma` still needs its edges and its weighting
       spectrum

Edges follow the canonical fast-first convention — strictly
DESCENDING, group ``0`` the fastest
(:ref:`canonical-group-convention`). The invariant is checked once, in
:class:`~orpheus.data.energy_grid.EnergyGrid`'s own construction, and
deliberately not re-checked on the axis.

Two consequences of the faces reading are worth stating explicitly,
because both were design forks:

- **Identity is** :math:`n_g` **plus the edges' CONTENT.** Not the
  ``EnergyGrid`` object. ⚠ ``Mixture.energy_grid`` mints a FRESH
  ``eq=False`` grid on every access, so ``is`` and ``==`` are both
  ``False`` for two reads of one mixture; an axis identity keyed on grid
  object identity would make two mints from one mixture disagree inside
  a single legitimate solve. The axis reads ``edges.tobytes()``.
- **A synthetic axis is NOT a from-grid axis.** ``synthetic(ng)`` is the
  honest spelling for a fixture or library that declares a group COUNT
  with no boundary energies — `[M]` every shipped ``get_mixture`` pair
  has ``eg is None``, so this is the common case, not the exotic one.
  Same index set, no partition data: a different axis.

.. _spaces-vv-collapse-hook:

The :math:`V` / :math:`V^*` collapse hook (declared, not yet built)
-------------------------------------------------------------------

The counting theorem has a sequel the axis is *recording data for*, and
this page declares it now so the eventual implementation is checked
against a stated contract rather than invented:

- **Condensation acts as a plain SUM on** :math:`V`. Integrals add:
  :math:`\Phi_G = \sum_{g \in G}\varphi_g`
  (:eq:`energy-condensation-coarse-flux`).
- **Condensation acts as a flux-weighted AVERAGE on** :math:`V^*`.
  Averages re-weight:
  :math:`\Sigma_G = \sum_{g\in G}\varphi_g\Sigma_g \big/
  \sum_{g\in G}\varphi_g`
  (:eq:`energy-condensation-vector-collapse`).
- **Collapse adjoint-consistency IS precisely that pair being mutually
  adjoint** under the counting pairing. Rate preservation
  (:eq:`energy-condensation-rate-preservation`) is what that
  adjointness says in physics vocabulary — which is why the corpus
  reads condensation as a Petrov-Galerkin projection rather than an
  averaging recipe (:ref:`sn-energy-condensation`).

.. warning::

   **This is a DECLARATION, not shipped machinery.** The pair above is
   what the group structure recorded on
   :class:`~orpheus.numerics.axis.EnergyAxis` exists to feed; the
   morphisms themselves are scheduled for a later phase / campaign 2
   (⚠ this row read "S7" until 2026-08-24: a bare plan-internal step
   number, and it COLLIDES with campaign 1 CS4b's own step S7, which
   landed that day and built none of this — trust the tree, not the
   number), and the condensation that ships today
   (:meth:`Mixture.condense
   <orpheus.data.macro_xs.mixture.Mixture.condense>`,
   :meth:`Solution.condense <orpheus.sn.solution.Solution.condense>`)
   does not route through an axis. Do not read this section as a
   description of a code path.


.. _spaces-collapse-doctrine:

The collapse doctrine: which axes survive a degeneracy, and why
===============================================================

This is the page's load-bearing content, and it is presented
**dialectically** — the question, then the two answers that were tried
and refuted *together with the questions that refuted them*, then the
doctrine that stands. That ordering is deliberate pedagogy, not
archaeology: the refuted versions are each *almost* right, and a reader
who meets only the final statement will re-derive one of them within a
week. The record was produced by an extended design dialogue on
2026-08-19/20 and is preserved at
``.claude/plans/cs1_energy_space_design.md`` Appendix A.


.. _spaces-collapse-the-question:

The question: two prior doctrines in tension
---------------------------------------------

The homogeneous infinite-medium solver, diffusion, S\ :sub:`N` and
P\ :sub:`N` differ in which tensor factors their spaces carry. The
architecture record that preceded this campaign stated two rules about
rank-one collapses, and they disagreed.

- **The retract rule.** Every canonical relation between spaces —
  the outflow half of a face slot, scalar inside angular, :math:`\ell=0`
  inside moments — has one shape: a projection :math:`\pi: V \to V'`, an
  embedding :math:`\iota: V' \to V`, with :math:`\pi\circ\iota =
  \mathrm{id}` and, under the right metrics, :math:`\iota = \pi^{H}` up
  to a scalar. On this reading a retract's codomain **DROPS** the
  collapsed axis: *scalar space is not angular space with a trivial
  axis*, and conflating the two is precisely the source of the classical
  :math:`4\pi` bookkeeping errors. What *is* a subspace of angular space
  is :math:`\iota(V_{\rm scal})`, the isotropic functions — a different
  object from :math:`V_{\rm scal}` itself.
- **The quotient rule.** The homogeneous solver is not measureless — it
  is maximally *quotiented*. Translation invariance quotients the
  spatial axis to a **one-point axis with unit measure**, and the "per
  unit volume" intensive convention *is* that normalized quotient
  measure. On this reading the collapsed axis **PERSISTS**.

Both collapses are rank-one. Both arguments are correct about their own
case. So the question the space layer has to answer is exactly:
**which rank-one collapses leave an axis behind, and what decides it?**


.. _spaces-collapse-version-1:

Version 1 — compactness (REFUTED)
----------------------------------

*Proposed.* A collapse over a **compact** domain integrates: the measure
is consumed by the integration, nothing problem-specific is left, and
the axis drops. That is angle — :math:`S^2` is compact, the total
:math:`4\pi` is universal. A collapse over a **non-compact** domain
cannot integrate (the integral diverges), so it must normalize, and the
normalization constant is problem data — the axis persists as the
quotient point. That is space — :math:`\mathbb{R}^d` under translations.

It reproduces both prior rules, it is a single criterion, and it is
wrong.

.. admonition:: ⛔ The refuting question — energy, and Bateman
   :class: error

   *Where does ENERGY sit?* Energy is collapsed by integration over
   intervals, so compactness would put it on the "integrate ⟹ drop"
   side. But the one-group flux manifestly **keeps** its axis — the
   shipped layout is :math:`(1, *\mathrm{spatial})`, not
   :math:`(*\mathrm{spatial})` — and it *must*: the group-averaged cross
   section :math:`\bar\sigma` is defined only relative to an interval
   and a weighting spectrum, and the Bateman / depletion pairing
   :math:`\langle\bar\sigma, \phi\rangle` consumes exactly that datum.
   Drop the axis and :math:`\bar\sigma` loses the thing that defines it.

   ⟹ **Compactness decides integrate-versus-normalize, at best. It
   cannot decide persist-versus-drop.** Those are two different
   questions, and Version 1 conflated them.


.. _spaces-collapse-version-2:

Version 2 — "energy is effectively finite-measure" (REFUTED)
-------------------------------------------------------------

*Proposed.* Patch energy in as a third case: it is *effectively*
finite-measure — the group structure is topped by a highest energy
:math:`E_0`, and the flux is integrable below it — so it behaves like a
compact domain for integration purposes while still carrying
problem-specific partition data.

.. admonition:: ⛔ The refuting question — what is the measure of :math:`(0, \infty)`?
   :class: error

   The energy domain has **infinite** Lebesgue measure. There is no
   upper limit on neutron energy; the grid top :math:`E_0` is a
   *library's* practicality, not a fact about the space. Calling energy
   "effectively finite-measure" smuggles a data-layer truncation into a
   space-layer definition, and it would make the doctrine depend on
   which library you loaded.

   The correct question is not about the domain at all. It is:
   **does the integral of that infinite tail CONVERGE?** — that is,
   is the admissible field class in :math:`L^1` over the collapsed
   domain?

   ⟹ **The discriminator is a property of the FIELDS the physics
   admits, never of the bare domain.** That reformulation is what
   unlocked the standing doctrine, and it is why the doctrine below
   talks about an *admissible field class* rather than about measures of
   sets.


.. _spaces-collapse-doctrine-standing:

The standing doctrine: two forks
---------------------------------

Version 1's real error was answering one question when there are two.
The doctrine separates them.

**Fork 1 — does the collapse INTEGRATE or NORMALIZE?** Decided by the
integrability of the **admissible field class** over the collapsed
domain. The mechanism is structural, and it runs through symmetry:

- A quotient **by a symmetry group** forces the fields to be CONSTANT
  along group orbits. If the group is **non-compact** (translations of
  :math:`\mathbb{R}^d`), its orbits carry infinite Haar measure, and a
  nonzero constant on infinite measure is never :math:`L^1`. Integration
  is impossible *structurally* — whatever the physics — and the only
  surviving functional is the normalized average, "per unit orbit
  measure".
- If the group is **compact** (rotations; :math:`S^2` is an orbit), Haar
  measure is finite, integration is always available, and the total
  (:math:`4\pi`) is canonical.
- With **no symmetry acting** — energy is the case — nothing forces
  constancy, and the *physics* decides integrability. For neutron
  spectra it holds: the fission spectrum :math:`\chi` decays
  super-exponentially above source energies and the thermal Maxwellian
  :math:`\to 0` as :math:`E \to 0`, so the tail integral converges. The
  practical grid top :math:`E_0` is the library's assertion that the
  neglected tail is below tolerance — a data-layer truncation
  statement, not a space-layer fact.

**Fork 2 — does the AXIS PERSIST?** Decided by **consultation**: the
axis survives iff the collapse leaves data that the surviving family
still consults — in its own pairing, or in its identity and guards. A
collapse that leaves only *re-embedding* conventions puts them on the
:math:`(\pi, \iota)` arrows and drops the axis; the family itself
changes.

Three clauses cover every collapse in the corpus:

.. list-table:: The three clauses
   :header-rows: 1
   :widths: 6 30 34 30

   * - #
     - Situation
     - Verdict
     - Instance
   * - **1**
     - collapse along **non-compact group orbits**
     - **NORMALIZE**; the axis **persists** as the quotient point with
       unit weight — the density convention is consulted by the
       member's OWN pairing
     - the homogeneous spatial slot (translations of
       :math:`\mathbb{R}^d`). In the modal branch the surviving datum is
       the mode parameter — buckling's :math:`B`.
   * - **2**
     - **partition-integration** of an :math:`L^1` class, no symmetry
       acting
     - **INTEGRATE** per cell; the **nodal mesh-axis persists** all the
       way down to its one-cell member — the partition (boundaries,
       weighting spectrum) is problem data consulted by identity, by
       guards and by the :math:`V`/:math:`V^*` pairing
     - energy. The one-group member is the one-CELL energy mesh, and
       :math:`\langle\bar\sigma, \phi\rangle` is exactly what consumes
       its edges. Likewise a genuine one-cell slab keeps its axis with
       weight :math:`V_{\rm cell} \neq 1`.
   * - **3**
     - **whole-domain integration** over a **compact canonical orbit**
     - the axis **DROPS** — the total is universal, so nothing
       problem-specific survives on the axis; the rebroadcast convention
       lives on the embedding :math:`\iota`, consulted only when LEAVING
       the family
     - angle. :math:`S^2` is a compact-group orbit and :math:`4\pi` is
       universal, so scalar spaces carry NO angular slot and are a
       different family.

.. admonition:: The two one-line tests
   :class: tip

   #. *Can the admissible fields be integrated over the collapsed
      domain?* — symmetry-forced constancy on infinite-measure orbits
      ⟹ **no** ⟹ normalize.
   #. *Is the surviving convention consulted INSIDE the family, or only
      at re-embedding?* — inside ⟹ **axis**; only at re-embedding ⟹
      **arrow**.

Notice what the doctrine does to the tension it was built to settle: it
does not pick a winner. **The retract rule and the quotient rule are
both right, about different clauses** — clause 3 is the retract, clause
1 is the quotient — and what was missing was the second fork, which
neither rule states.


.. _spaces-collapse-retrodictions:

Retrodictions: the doctrine against the shipped tree
------------------------------------------------------

A doctrine invented to settle one dispute is worth little if it only
settles that dispute. The rows below are layouts the doctrine was NOT
built from, and the clause column is its *prediction* of each.

⚠ **One row is not a retrodiction — read the third column.** The
buckling member is a design conclusion of the same dialogue, not
something the tree ships; it is listed because the doctrine's clause-1
modal branch is what predicts its shape, and a table of confirmations
that quietly includes an aspiration is the exact defect this corpus
otherwise catches in plans.

.. list-table::
   :header-rows: 1
   :widths: 62 14 24

   * - Fact
     - Clause
     - Status
   * - the homogeneous shape :math:`(n_g, 1)` — the spatial point is
       PRESENT, and its unit volume is consumed by the reaction-rate
       functionals inside the solve
     - 1
     - `[M]` **ships**
   * - the scalar family :math:`(n_g, *\mathrm{spatial})` — **no**
       angular slot; the :math:`4\pi` conventions live at the embeddings
     - 3
     - `[M]` **ships**
   * - the one-group layout :math:`(1, *\mathrm{spatial})` — the energy
       axis persists, with its edges
     - 2
     - `[M]` **ships**
   * - one-cell meshes keep weight :math:`V_{\rm cell} \neq 1`, and are
       therefore NOT the quotient point
     - 2
     - `[M]` **ships** (gated; see below)
   * - partial currents: a hemisphere collapse parameterized by
       :math:`\hat n` — the parameter survives on the FACE structure,
       the angular content on the arrow
     - 3 (+ a face axis)
     - `[M]` **ships**
   * - the buckling member: a size-1 **MODAL** spatial axis carrying
       :math:`B`, live angle, field :math:`\mathbb{F} = \mathbb{C}`
     - 1 (modal branch)
     - ⛔ **NOT built** — a prediction (campaign 2)

The fourth row is the one worth pausing on, because it is where the
doctrine becomes *mechanically* enforced rather than merely stated. A
quotient point and a genuine one-cell mesh have the same shape and
differ only in measure — and since the measure is part of axis identity
and axis identity derives the space name, they are **unequal spaces**
(:ref:`spaces-identity-bridge`). Nothing else in the tree can see the
difference: `[M]` a scalar metric commutes with every operator, so
:math:`A^\dagger = G^{-1}A^{\mathsf T}G = A^{\mathsf T}` whenever
:math:`G = cI` — the quotient-vs-one-cell distinction is provably
invisible to ``.H``, to every norm and to every value gate. **Space
identity is the only instrument that carries it**, which is why the
CS1 battery guards it with an identity gate and pairs that with an
explicit must-stay-green proof that no value gate can.


.. _spaces-quotient-family:

Consequences: the quotient family, and what a degenerate axis stores
---------------------------------------------------------------------

Clause 1 says the homogeneous solver's spatial axis persists. The
obvious follow-up is *so what does it store, and why is that useful?*
The answer is the **density convention** — the quotient's unit
normalization weight, made a first-class, composable object rather than
an unwritten agreement. Three payoffs, in increasing order of how much
they hurt when the convention is implicit:

#. **The pairing consumes it, so rates follow automatically.** If the
   convention ever changes — per lattice cell of volume
   :math:`V_{\rm cell}`
   instead of per unit volume, say — every functional follows by
   arithmetic, because the weight is the thing they integrate against.
   `[M]` **this is live, and it is measurable — since CS4a K2, of the
   SPACE's own measure.**
   :func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`
   normalises its flux through the space's pairing
   (``space.inner_product(Σx, φ)``), which contracts against the
   posing's per-axis weights — so the quotient point's weight is
   genuinely consumed, not decorative. The separated measurement (each
   measure varied ALONE, 2026-08-21): minting the point with weight 2.0
   moves the flux and both specific rates by exactly the measure ratio
   while :math:`k_\infty` is unchanged, and doubling the CARRIER's cell
   volumes moves **nothing** — the carrier supplies cross sections, not
   the measure. (Until K2 the rates read ``mesh.volume_measure``
   instead: the same experiment then read :math:`0.225` vs :math:`0.450`
   between the quotient carrier, weights ``[1.0]``, and a one-cell slab
   of width 2, weights ``[2.0]`` — a true measurement whose two measures
   were varied TOGETHER, so it could not distinguish which one the rate
   consumed; the pre-K2 value path was in fact bit-identically inert to
   the space weight.)
#. **Family coherence.** The :math:`B = 0` fiber and the :math:`B \neq
   0` members share the slot, so fiberwise machinery — Fourier
   convergence analysis, :math:`\rho(B)` — reads a uniform signature
   across the family instead of special-casing the degenerate member.
#. **Boundary-of-family maps get a home and a Jacobian.** Pushing
   homogenized constants into a meshed context is a map between two
   members of the family; with the measure on the axis, the conversion
   is the **measure ratio** — an object with a name — instead of an
   invisible unit-convention shift. That invisible shift is the
   classical missing-volume-factor bug class in homogenization.

The quotient family, and where the modes live
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The homogeneous solver is the **terminal** member of a family, not an
isolated special case:

- **The infinite-medium member.** Translation invariance quotients the
  spatial axis to the one-point axis with unit measure; isotropy
  quotients angle to the :math:`\ell = 0` retract (clause 3 — the
  angular axis is GONE, not trivial). The spectrum problem lives on the
  energy axis tensored with the quotient point, and everything transfers
  with no special cases: the energy axis is nodal, so the coordinate
  cone, the Hilbert-metric brackets, ray normalization and the
  irreducibility gate all apply verbatim
  (:doc:`/theory/foundations/infinite_medium`).
- **The buckling / B\ :sub:`1` member is the INTERMEDIATE quotient** — a
  one-dimensional **MODAL** spatial axis parameterized by :math:`B`, on
  which streaming is multiplication by :math:`iB\mu`. It is modal
  because a Fourier mode is a coefficient, not a cell value; the field
  is :math:`\mathbb{C}`; and Fourier convergence analysis is exactly
  *the solver diagonalized over this quotient family*, computed
  fiberwise. The :math:`B = 0` fiber is the infinite-medium member, and
  it is the same slot — which is payoff 2 above, stated concretely.
- **P**\ :sub:`N` **is the ANGULAR buckling.** The two hierarchies are
  the same construction on two different groups: irreps of translations
  give the buckling ladder on the spatial axis, irreps of rotations give
  the :math:`P_N` ladder on the angular axis. The parallel also
  *predicts* clause 3's asymmetry: the trivial angular slot stays
  ABSENT on scalar-family members, exactly as the trivial spatial slot
  would if space were compact.

.. note::

   **CS1 does not build the buckling member; it only refuses to
   foreclose it.** The complex field :math:`\mathbb{F} = \mathbb{C}` and
   the modal spatial axis are scheduled work (campaign 2). What CS1
   guarantees is that ``MODAL`` exists as a first-class basis kind from
   day one, so the member can be minted later without re-typing the
   axis, and that the consumers which cannot answer on a modal factor
   (the cone predicate) refuse rather than silently answer.

.. _spaces-symmetry-monotonicity:

The symmetry-monotonicity law
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Clause 1 quotients "by the symmetry group" — but *whose* symmetry?
Symmetry lives at the **geometry**, and every arrow down the modelling
lattice at best preserves it and usually reduces it. It never increases
it:

.. math::

   G_{\rm medium} = G_{\rm geom} \cap \operatorname{Stab}(\text{material
   assignment}),
   \qquad
   G_{\rm mesh} = G_{\rm geom} \cap \operatorname{Stab}(\text{cells}),
   \qquad
   G_{\rm pullback} \subseteq G_{\rm medium} \cap G_{\rm mesh}.

A material assignment breaks the geometry's symmetry unless the
assignment is itself invariant; a uniform mesh keeps discrete
translations while an unstructured one keeps nothing.

The consequence for the collapse doctrine is precise and it settles a
"where does this belong?" question that the mesh vocabulary could not
answer: **clause 1's quotient consumes the MEDIUM's surviving symmetry,
not raw geometry's.** The infinite homogeneous medium is exactly the
member whose assignment stabilizer is *everything* — a uniform
assignment breaks nothing — so the full translation group survives to be
quotiented. Buckling then restricts to an irrep of that surviving group.

.. warning::

   **The** ``Medium`` **layer this law argues for is CHARTERED, not
   shipped.** The lattice the law is stated over —
   geometry :math:`\to` {medium, mesh} :math:`\to` the pullback that
   carries materials on cells — is a design conclusion of the same
   dialogue, scheduled as its own phase (CS1.5) between CS1 and CS2.
   Today the homogeneous path still constructs a degenerate
   :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` through
   ``from_materials``, which is the pullback wearing a constructor story
   it does not have a mesh for. Read this subsection as *why the
   quotient is licensed one level above the mesh*, not as a description
   of a class that exists.



.. _spaces-collapse-pair:

The realized machinery: the axis collapse pair
==============================================

The doctrine above decides **which** axes survive a degeneracy. This
section is the machinery that performs the collapse when the verdict is
"drop", and it is where the doctrine stops being a rule and becomes two
typed arrows a call site can hold.

The bridge is the retract rule's own sentence
(:ref:`spaces-collapse-the-question`): *a projection*
:math:`\pi: V \to V'`, *an embedding* :math:`\iota: V' \to V`, *with*
:math:`\pi\circ\iota = \mathrm{id}` *and, under the right metrics,*
:math:`\iota = \pi^{H}` *up to a scalar*. Everything below is that
sentence made precise and made executable:

- :math:`\pi` is :class:`~orpheus.numerics.operator.AxisRetractionOperator`,
  minted by :meth:`FunctionSpace.retraction
  <orpheus.numerics.space.FunctionSpace.retraction>`;
- :math:`\iota` is :class:`~orpheus.numerics.operator.AxisSectionOperator`,
  minted by :meth:`FunctionSpace.section
  <orpheus.numerics.space.FunctionSpace.section>`;
- **the scalar is** :math:`\Sigma w`, the axis's total mass, and it is
  not a convention: it is the :math:`1\times 1` Gram of the rank-one
  frame that mints the pair (:ref:`spaces-collapse-pair-frame`);
- the two are **different types** precisely so that the "up to a
  scalar" cannot be silently dropped at a call site — which is the
  ERR-051 failure class made unspellable.

.. code-block:: python

   from orpheus.numerics.axis import Axis, BasisKind
   from orpheus.numerics.space import FunctionSpace

   V = FunctionSpace.of_axes(
       Axis("angular", (4,), weights=w, kind=BasisKind.NODAL),
       Axis("energy",  (2,),            kind=BasisKind.NODAL),
       Axis("spatial", (5,), weights=V_cell, kind=BasisKind.NODAL),
   )
   R = V.retraction("angular")   # V -> energy (x) spatial   (fiber integration)
   E = V.section("angular")      # energy (x) spatial -> V   (the section)

   R.apply(E.apply(phi))         # == phi
   R.H.apply(phi)                # the plain broadcast — NOT E

Both verbs are **memoized on the space**, one mint per space per axis
label, so a carrier that caches its spaces (every solver carrier does)
gets warm operators for free: ``sn.angular_bulk_space.retraction("angular")``
builds the generator once and returns the same object thereafter
(`[M]` identity, gated).


.. _spaces-collapse-pair-two-arrows:

Two arrows, not one — and the scalar between them
--------------------------------------------------

Write :math:`V = V_{\rm ax} \otimes V'` for the product whose first
factor is the axis being collapsed, with the product metric
:math:`G_V = \operatorname{diag}(w) \otimes G_{V'}`
(:eq:`spaces-axis-product`). The **retraction** is fiber integration
over that factor — the pushforward :math:`\pi_*` along the projection
that forgets the axis:

.. math::
   :label: spaces-collapse-retraction

   (R\,\psi)(\cdot) \;=\; \sum_{n} w_n\,\psi(n,\cdot),
   \qquad
   R : V \longrightarrow V' .

.. (vv-status rationale) Structural/representational identity: it STATES
   what ``AxisRetractionOperator.apply`` computes (the axis's factor
   measure contracted over the axis's ndarray dims). Not a solver claim
   — no flux, no eigenvalue, no discretization error. The verifiable
   content is the CS4b S6 foundation battery
   (``tests/numerics/test_axis_marginal.py``): the tightness row pins
   this contraction against the mint frame's own analysis content, and
   G6.5 pins it bit-identically against a hand-spelled einsum on the
   real S\ :sub:`N` carrier.
.. vv-status: spaces-collapse-retraction documented

.. implements:: spaces-collapse-retraction
   :by: py:method:orpheus.numerics.operator.AxisRetractionOperator.apply

   **Implemented by** 2 sites. The kernel is the operator's; the
   canonical field-level consumer re-spells nothing — since CS4b S6.2
   :meth:`AngularField._integrate_angular_values
   <orpheus.transport.fields._bases.AngularField._integrate_angular_values>`
   IS this ``apply``, so the tree has one realization of the angular
   reduction and not two.

.. implements:: spaces-collapse-retraction
   :by: py:method:orpheus.transport.fields._bases.AngularField._integrate_angular_values

Its Hilbert adjoint is where the two-arrow discipline is forced. Under
the product metric the axis weights **cancel exactly**:

.. math::
   :label: spaces-collapse-adjoint-is-pullback

   R^{\dagger}
   \;=\; G_V^{-1}\,R^{\mathsf T}\,G_{V'}
   \;=\; \bigl(\operatorname{diag}(w)\otimes G_{V'}\bigr)^{-1}
         \bigl(\operatorname{diag}(w)\otimes G_{V'}\bigr)\,\pi^{*}
   \;=\; \pi^{*},

.. (vv-status rationale) Structural identity of the metric sandwich: the
   Euclidean transpose of a weighted contraction is the WEIGHTED
   scatter, and the product metric's own axis block is exactly what
   removes those weights again — so the Hilbert adjoint of fiber
   integration is the UNWEIGHTED broadcast. Representational, not a
   solver claim (no physics enters; it is true for every axis measure).
   The verifiable content is the CS4b S6 foundation battery's G6.3 row
   in ``tests/numerics/test_axis_marginal.py``, which pins the
   adjunction on the physical metrics and carries the vv #19 NEGATIVE
   leg (the same pairing under a deliberately stripped spatial measure
   must break at O(1)).
.. vv-status: spaces-collapse-adjoint-is-pullback documented

.. no-implementation:: spaces-collapse-adjoint-is-pullback
   :kind: identity

   **Nothing implements this.** It is an identity between two things
   that are each computed elsewhere and never equated in production:
   the left side is produced by the generic metric-aware adjoint
   wrapper (``R.H``, which knows nothing about axes), the right side by
   :meth:`AxisSectionOperator.apply
   <orpheus.numerics.operator.AxisSectionOperator.apply>` scaled by
   :math:`\Sigma w`. No line forms the comparison — that is the point:
   the cancellation is what lets the adjoint be free rather than
   bespoke. It is *measured* by the G6.4 row of
   ``tests/numerics/test_axis_marginal.py``.

where :math:`\pi^{*}` is the **pullback** — the plain, unweighted
broadcast of :math:`\varphi` across the axis. So
:math:`(R, R^{\dagger}) = (\pi_*, \pi^{*})` is the discrete realization
of the fiber-integration / pullback adjunction, and the pairing

.. math::

   \langle R\psi,\ \varphi\rangle_{V'}
   \;=\;
   \langle \psi,\ \pi^{*}\varphi\rangle_{V}

holds on the *physical* metrics with no correction factor at all
(`[M]` bounded by :math:`3.7\times10^{-13}` relative over 200 draws on
the shipped S\ :sub:`N` carrier and :math:`1.6\times10^{-13}` on the
synthetic three-axis fixture — a scalar-valued identity whose relative
residual is set by cancellation in the inner products, not by the
operators; see :ref:`spaces-collapse-pair-evidence`).

**The pullback is not the section.** :math:`\pi^{*}` broadcasts
:math:`\varphi` unchanged, so :math:`R\pi^{*}\varphi = (\Sigma w)\varphi`
— it is a right inverse of :math:`R` only after division by the axis's
total mass. That division is the second arrow:

.. math::
   :label: spaces-collapse-section

   (E\,\varphi)(n,\cdot) \;=\; \frac{\varphi(\cdot)}{\Sigma w},
   \qquad
   E \;=\; \frac{\pi^{*}}{\Sigma w}
   \;=\; R^{\dagger}\,(R\,R^{\dagger})^{-1},
   \qquad
   R\circ E \;=\; \mathrm{id}_{V'} .

.. (vv-status rationale) Structural/representational identity: it STATES
   what ``AxisSectionOperator.apply`` computes (divide by the axis's
   total mass, then broadcast) and identifies it as the Moore-Penrose
   pseudo-inverse of the retraction in the two spaces' own metrics. Not
   a solver claim. The verifiable content is the CS4b S6 foundation
   battery in ``tests/numerics/test_axis_marginal.py``: G6.1 (the
   section law), G6.2 (idempotence of the composite projector), G6.6
   (bit-identity with the shipped isotropic-source kernel on the real
   S\ :sub:`N` carrier) and the gram-derivation row that pins the
   divisor against the mint frame's own Gram entry.
.. vv-status: spaces-collapse-section documented

.. implements:: spaces-collapse-section
   :by: py:method:orpheus.numerics.operator.AxisSectionOperator.apply

   **Implemented by** 3 sites. The kernel is the operator's; the mint
   supplies the divisor from the frame's Gram, and the canonical
   producer-side consumer re-spells nothing — since CS4b S6.2
   :meth:`AngularSourceSink.from_isotropic
   <orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink.from_isotropic>`
   IS this ``apply``.

.. implements:: spaces-collapse-section
   :by: py:function:orpheus.numerics.frame._collapse_pair

.. implements:: spaces-collapse-section
   :by: py:method:orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink.from_isotropic

Since :math:`R` is rank-deficient by construction, *many* right inverses
exist — putting all the mass on one ordinate, :math:`(\iota_0
\varphi)(n,\cdot) = \delta_{n0}\,\varphi/w_0`, is one. :math:`E` is
distinguished among them by being the **minimum-norm** one, which is
what :math:`R^{\dagger}(RR^{\dagger})^{-1}` says: it is the
Moore–Penrose pseudo-inverse in the two spaces' own metrics. `[M]` on
the synthetic fixture and one draw (``default_rng(21)``),
:math:`\lVert E\varphi\rVert_V = 1.700` against
:math:`\lVert \iota_0\varphi\rVert_V = 4.389` for the one-ordinate
right inverse — a factor 2.6 — with both satisfying
:math:`R\circ(\cdot) = \mathrm{id}` at the round-off floor.

The four arrows close into a square, and every entry is one of the two
operators or one of their adjoints:

.. list-table:: The collapse square (:math:`W \equiv \Sigma w`)
   :header-rows: 1
   :widths: 22 30 48

   * - Arrow
     - What it is
     - Identity
   * - :math:`R = \pi_*`
     - fiber integration (the weighted axis sum)
     - :math:`R\,R^{\dagger} = W\,\mathrm{id}`
   * - :math:`R^{\dagger} = \pi^{*}`
     - the pullback (the **plain** broadcast)
     - :math:`R^{\dagger} = W\,E` — the anti-ERR-051 scalar
   * - :math:`E`
     - the section (broadcast **after** dividing by :math:`W`)
     - :math:`R\circ E = \mathrm{id}`, and :math:`E` is minimum-norm
   * - :math:`E^{\dagger}`
     - the :math:`w`-mean (average over the axis)
     - :math:`E^{\dagger} = R/W`, and :math:`E^{\dagger}\circ R^{\dagger} = \mathrm{id}`

and the composite in the other order,

.. math::

   P \;\equiv\; E\circ R,
   \qquad
   (P\psi)(n,\cdot) \;=\; \frac{1}{\Sigma w}\sum_m w_m\,\psi(m,\cdot),

is the **conditional expectation onto axis-constant functions** — the
:math:`w`-mean projector. It is idempotent and, because
:math:`P^{\dagger} = R^{\dagger}E^{\dagger} = (WE)(R/W) = P`,
self-adjoint in :math:`G_V`: an *orthogonal* projector, not merely an
oblique one (`[M]` self-adjointness bounded by
:math:`4.0\times10^{-14}` relative over 200 draws — again a
cancellation-limited scalar identity, see
:ref:`spaces-collapse-pair-evidence`).
That is the precise sense in which "the isotropic part of
:math:`\psi`" is a well-defined object rather than a convention.

.. warning::

   **An axis whose SIGNED measure sums to zero has NO section, and the
   retraction over it is still legal.** The asymmetry is structural, not
   defensive: :math:`R` is a contraction and needs no division, while
   :math:`E` divides by :math:`\Sigma w` — and at :math:`\Sigma w = 0`
   the rank-one Gram is singular, so the mint frame has no canonical
   dual and no section EXISTS to hand back. The mint therefore leaves
   that arm unminted and :meth:`FunctionSpace.section
   <orpheus.numerics.space.FunctionSpace.section>` refuses at access,
   naming the cause. Signed axis weights are deliberately legal on
   :class:`~orpheus.numerics.axis.Axis` (a :math:`\sigma`-folded
   quadrature can carry them), so this is a reachable state and not a
   theoretical one.

.. note::

   **This is where the doctrine's "up to a scalar" gets its value.**
   The retract rule said :math:`\iota = \pi^{H}` *up to a scalar* and
   left the scalar unnamed, which is exactly the gap the classical
   :math:`4\pi` bookkeeping errors live in. Here the scalar is
   :math:`\Sigma w`, it is read off the mint frame's Gram, and the two
   arrows carry different **types** — so a call site that reaches for
   ``R.H`` where it wanted ``E`` does not silently rescale a source by
   :math:`4\pi`; it holds an object of the wrong class. `[M]` the gate
   asserts precisely that ``R.H`` is *not* an
   :class:`~orpheus.numerics.operator.AxisSectionOperator`.


.. _spaces-collapse-pair-naming:

Why "retraction" and "section" — and why "embedding" was rejected
------------------------------------------------------------------

The names are canonical, and the reason is worth stating because the
first spelling shipped and was replaced within a day.

:math:`R\circ E = \mathrm{id}` makes :math:`(R, E)` a **split
epi/mono pair** in the categorical sense (Mac Lane, *Categories for the
Working Mathematician*, §I.5): a morphism with a right inverse is a
**split epimorphism** and the right inverse is called a **section**;
dually the left inverse of a split monomorphism is called a
**retraction**. The collapse doctrine's own prose already used "the
retract rule", so :math:`R` inherited the right word immediately.

The first implementation (S6.0) called the other arrow
``AxisEmbeddingOperator``. That name was retired the same day for two
independent reasons:

#. **It cannot discriminate the pair it exists to discriminate.** Any
   injective structure-preserving map is an embedding — and
   :math:`\pi^{*} = R^{\dagger}` is one too. So "embedding" names a
   property both arrows have, in a design whose entire purpose is to
   keep :math:`R^{\dagger}` and :math:`E` apart. The
   :math:`\Sigma w` weld the ERR-051 class is made of was hiding in
   the *name*.
#. **The object is defined relative to** :math:`R`. :math:`E` is not "a
   map into :math:`V` that happens to be injective"; it is *the* right
   inverse of a specific retraction — :math:`R\circ E = \mathrm{id}`
   IS its definition. The categorical name for that is **section**, and
   nothing weaker carries the relation.

"Embedding" survives in this corpus only as a generic adjective (the
doctrine's :math:`\iota`, the :math:`S^2` embedding of
:doc:`/theory/foundations/spherical_harmonics`). It is not the name of
an operator.

.. note::

   **The naming rule this instantiates** is the same one
   :ref:`the frame hierarchy <frame-discipline-as-a-type>` follows: *a
   reader of a type name knows its properties without reading the
   docstring.* ``AxisSectionOperator`` tells a reader that composing it
   after the retraction is the identity; ``AxisEmbeddingOperator`` told
   them only that it was injective, which is true of the wrong arrow
   too.


.. _spaces-collapse-pair-frame:

The pair is frame-induced — a stage-2 generator at rank one
------------------------------------------------------------

The kernels above could be written by hand: an einsum and a broadcast,
six lines each. They are not, and the reason is Cardinal Rule 2. The
hand-written pair is a **concept-level twin** of machinery this corpus
already ships — the discrete frame
(:doc:`/theory/foundations/frame`) — specialized to rank one. Two
places would have had to agree, forever, about what "the axis measure"
means and what the normalization divisor is.

So the pair is *induced*. At the mint site
(:func:`orpheus.numerics.frame._collapse_pair`) a literal frame is
built over the axis's index set, read for its induced data, and
discarded:

.. code-block:: python

   frame = GalerkinFrame(
       basis=IndicatorBasis(edges_per_axis=(np.array([-0.5, n - 0.5]),)),
       measure=DiscreteMeasure(
           nodes=np.arange(n, dtype=float),
           weights=flat_weights,              # weights=None IS the counting measure
           support=f"index({axis_label})",
       ),
   )
   kernel_weights = frame.measure.weights          # the analysis face's content
   total_weight   = float(frame.discrete_gram[0, 0])   # the rank-one Parseval metric

The basis is a **single-region indicator** covering every index
:math:`\{0,\dots,n-1\}` — one region, so exactly one coefficient, so a
:math:`1\times1` Gram. Under that basis every table entry is
:math:`1`, and the frame's Gram
(:math:`G_{jk} = \sum_n w_n \phi_j(x_n)\phi_k(x_n)`) collapses to the
axis's total mass:

.. math::
   :label: spaces-collapse-rank-one-gram

   G \;=\; \bigl[\textstyle\sum_n w_n\bigr] \;=\; [\,\Sigma w\,],
   \qquad
   E \;=\; R_{\rm frame}\circ G^{-1},

.. (vv-status rationale) Literature-transcribed / structural: it states
   the rank-one instance of the Parseval theorem already derived at
   :ref:`frame-parseval-metric` (the frame's codomain metric is the
   INVERSE discrete Gram), specialized to a single-region indicator
   basis where the Gram is 1x1 and its entry is the measure's total
   mass. Not a solver claim. The verifiable content is the
   gram-derivation row of ``tests/numerics/test_axis_marginal.py``,
   which asserts the section's divisor IS the literal frame's
   ``discrete_gram[0, 0]``, and the tightness row, which pins the
   minted kernels against that frame's own face contents.
.. vv-status: spaces-collapse-rank-one-gram documented

.. implements:: spaces-collapse-rank-one-gram
   :by: py:function:orpheus.numerics.frame._collapse_pair

   **Implemented by** 1 site. The mint is where the identity is
   *used* — it reads ``frame.discrete_gram[0, 0]`` and stores it as the
   section's divisor. The Gram itself is the frame's generic
   :math:`O(NK^2)` einsum, not a rank-one special case, which is the
   whole point: nothing in the collapse pair re-derives what a frame
   already computes.

so the section's divisor **is** the Parseval metric of
:ref:`frame-parseval-metric` at :math:`K = 1`. "Divide by
:math:`4\pi`" is therefore not a convention this code chose; it is the
inverse Gram of a frame, obtained the same way every other metric in
the corpus is obtained.

The frame is discarded on return. That is deliberate, and it is the
**stage-2 generator discipline** — the ruling this section exists to
realize (user, 2026-08-24):

   A stage-2 generator induces structure on both the space and the
   operator, and the two inductions must be minted together, at one
   site. Frame: induces the HarmonicAxis metric (space side), mints
   Analysis/Synthesis (operator side) — consistency is the tightness
   gate. Scheme: induces the trace descriptor and basis kind (space
   side), mints the closure (operator side) — consistency is one
   closure serving both apply and solve, which is ERR-026's structural
   closing. Mesh and Quadrature are the degenerate cases (space side
   only). Forgetting = retaining the induced parts; accessors are
   provenance.

Three consequences, each visible in the code:

**Both inductions at one site.** The mint constructs the retraction and
the section together and returns them as a pair. There is no path that
produces one without the other, so they cannot disagree about the axis,
the dims, or the marginal space — `[M]` the two arrows share ONE
marginal-space instance (``R.codomain is E.domain``), gated.

**Forgetting means copying the induced parts out.** The operators
retain the bound spaces, the ndarray dims the axis occupies, the flat
weights, and the scalar divisor — and nothing else. In particular they
do **not** retain a frame *face*: a face is a view holding
``frame:`` (:class:`~orpheus.numerics.frame.FrameBase`'s
``_FrameAnalysis`` / ``_FrameReconstruction``), so keeping one would
keep the generator alive through it. `[M]` read the face
dataclasses — this is why the mint copies arrays rather than storing
``frame.analysis``.

**Consistency is a gate, not an instance.** Because the generator is
thrown away, nothing at runtime *forces* the minted kernels to agree
with the frame that produced them. The **tightness gate** supplies
that: it rebuilds the literal frame independently and pins all three
correspondences — :math:`R` against the analysis content,
:math:`R^{\mathsf T}` against ``analyze_transpose``, and :math:`E`
against reconstruction composed with :math:`G^{-1}`. `[M]` all three
are bit-exact on **200 of 200** draws, which is a *stronger* statement
than the section law two sections up: there the exactness is a property
of the draw, here it is a property of the construction, because the two
sides evaluate the same reduction in the same order. The operator
kernels are hand einsums and the frame path runs the basis's table
einsums, so these are two different float programs and agreement is a
real claim, not a tautology.

.. note::

   **The latent generalization, recorded and not built.** The mint is
   parameterized by exactly one choice — the basis. Swap the
   single-region ``IndicatorBasis`` for a single-region
   ``WeightedIndicatorBasis`` and the same site produces a *profiled*
   collapse: the Petrov-Galerkin test side of a
   :math:`\chi`-class emission collapse, where the axis is not
   averaged uniformly but against a spectrum. That is the same
   machinery, a different basis, and it is built when its consumer
   lands (CS4c) — not before. Recording it here is what stops a future
   session minting a second, parallel mechanism for it.


.. _spaces-collapse-pair-clause-gate:

The clause gate: which axes admit, and why energy refuses
-----------------------------------------------------------

The mint refuses an axis the doctrine says must persist. This is the
one place in the tree where the collapse doctrine is **enforced** rather
than described, so it is worth reading the admission table against
:ref:`spaces-collapse-doctrine-standing` clause by clause.

.. list-table:: Admission at the mint
   :header-rows: 1
   :widths: 26 12 62

   * - Axis
     - Clause
     - Verdict, and why
   * - **angle** (an untyped ``Axis`` today)
     - 3
     - **ADMIT.** Whole-domain integration over a compact canonical
       orbit: the total is universal, nothing problem-specific survives
       on the axis, so the drop-form marginal is exactly right and the
       re-broadcast convention lives on the arrows :math:`E` /
       :math:`\pi^{*}`.
   * - a typed :class:`~orpheus.numerics.axis.EnergyAxis`
     - 2
     - **REFUSE**, with a pointer. Partition-integration of an
       :math:`L^1` class: the energy axis PERSISTS at its one-cell
       member, because :math:`\langle\bar\sigma,\varphi\rangle` consumes
       the partition. A drop-form marginal here would be a second
       mechanism for a collapse the tree already implements as
       **condensation**.
   * - a ``MODAL`` axis
     - —
     - **REFUSE.** Contracting expansion COEFFICIENTS with the basis
       mass is not an integral of the represented function. The modal
       average is the coefficient at the average slot; slice it.
   * - an untyped generic axis, whatever its label
     - 3
     - **ADMIT.** The gate reads the axis's TYPE, never its label
       string.
   * - the only axis of a single-axis space
     - —
     - **REFUSE.** Its marginal would be a bare scalar, which is not a
       :class:`~orpheus.numerics.space.FunctionSpace`. Contract with
       the space's inner product instead.
   * - a space with ``axes is None``
     - —
     - **REFUSE.** A densified legacy product has no named factors to
       marginalise over.

The energy row is the load-bearing one, and it is the doctrine's
Version-1 refutation cashed out in code
(:ref:`spaces-collapse-version-1`). Energy is collapsed *by
integration*, so a compactness-flavoured reading puts it on the "drop"
side; the shipped one-group layout is :math:`(1, *\mathrm{spatial})`
and keeps its axis, because :math:`\bar\sigma` is defined only relative
to an interval and a weighting spectrum. The machinery that performs an
energy collapse correctly therefore already exists — it is
:meth:`EnergyGrid.overlap_to
<orpheus.data.energy_grid.EnergyGrid.overlap_to>` and the
Petrov-Galerkin condensation frames of
:ref:`sn-energy-condensation` — and the refusal message names it, so
a caller who reaches for the wrong tool is handed the right one rather
than a wrong answer.

The clause gate is not hypothetical on a shipped carrier. `[M]` on the
:math:`S_N` fixture above, whose ``angular_bulk_space`` axes are
``(Axis, EnergyAxis, Axis)``:

.. code-block:: text

   sn.angular_bulk_space.retraction("angular")  -> OK, codomain (2, 5)
   sn.angular_bulk_space.retraction("spatial")  -> OK, codomain (4, 2)
   sn.angular_bulk_space.retraction("energy")   -> TypeError: ... is a typed
       EnergyAxis, which PERSISTS at its one-cell member (collapse doctrine
       clause 2 ...). The energy collapse is condensation: use
       EnergyGrid.overlap_to / the Petrov-Galerkin condensation frame, not a
       drop-form marginal.

Two of the carrier's three factors marginalise; the one the doctrine
says must persist refuses, by TYPE, with the successor named in the
message.

.. warning::

   **The gate reads the TYPE, and today only energy has one.** A
   generic ``Axis(label="energy", ...)`` — a synthetic test factor — is
   ADMITTED, and correctly so: refusing on the label string would be
   stringly-typed dispatch, and it would refuse a legitimate synthetic
   fixture that carries none of energy's physics. The consequence is
   that the clause gate is **structural for energy and permissive for
   everything else** until CS2 lands the typed spatial / quadrature /
   harmonic axes; at that point the verdict becomes axis-family
   polymorphism and each family answers for itself. Until then, do not
   read "the mint admitted it" as "the doctrine says it drops".


.. _spaces-collapse-pair-refuted:

What was tried, and what refuted it
------------------------------------

**The hand-derived pair (S6.0) — superseded within a day.** The first
implementation minted the two operators from the space with
hand-spelled kernels and a hand-chosen divisor
(``weights.sum()``). It was correct and it was a twin: the design
dialogue that followed established the exact correspondence with the
rank-one frame, which is Cardinal Rule 2's stop condition. The
re-carve (S6.0b) kept the operator shells verbatim — admission, dims
bookkeeping, bound spaces, the einsum/broadcast kernels — and changed
only where the two retained numbers come from. That is why the re-carve
is bit-identity-safe by construction and why the equivalence gates
survive it unchanged.

**An** ``Axis`` **→ measure accessor — refused.** The mint needs a
:class:`~orpheus.numerics.measure.DiscreteMeasure` built from the axis.
The obvious move is to give :class:`~orpheus.numerics.axis.Axis` a
public accessor that produces one. It was rejected under the same
ruling that governs the frame: *accessors are provenance*. An accessor
would make the generator reachable from the axis forever, so the axis
would carry a permanent dependence on frame machinery it does not
need. The mint builds the measure with a **local** helper instead; the
axis stays four slots and nothing more.

**Caching the pair on the** ``Axis`` **— refused for the same reason,
and it would have been wrong anyway.** The pair is not a function of
the axis alone: its domain is the whole product and its codomain is the
product minus that factor, so two spaces sharing an axis have
*different* collapse pairs. The cache belongs to the space, and it is
there (memoized in the frozen dataclass's ``__dict__``, one entry per
axis label).

**Retaining the frame on the operators — refused.** Keeping the frame
would make consistency automatic instead of gated, which sounds
strictly better. It is not: it retains a whole generator (basis, table,
measure, two spaces) on every collapse operator in the tree to secure a
property that a two-line gate already secures, and it re-opens the
question of whether two operators built from *equal* frames are the
same operator. Consistency here is carried by content-determinism plus
the tightness gate, not by instance sharing. The rule that would flip
this: a second consumer needing the *identical frame instance* for
measure-consistency or anti-aliasing. None exists today.


.. _spaces-collapse-pair-evidence:

Numerical evidence
-------------------

`[M]` 2026-08-24, measured against the tree at HEAD. **The construction,
so the tables regenerate from this page.** The *synthetic* fixture is
the three-axis product ``angular(4, w=[0.3, 0.7, 0.5, 0.5]) ⊗
energy(2, counting) ⊗ spatial(5, V=[0.2, 0.3, 0.4, 0.7, 1.4])`` built
with :meth:`FunctionSpace.of_axes
<orpheus.numerics.space.FunctionSpace.of_axes>`; the weights are
non-uniform on purpose so that no cancellation flatters a law. The
:math:`S_N` fixture is the shipped carrier —
``SNMesh(Mesh1D(edges=[0, 0.2, 0.5, 0.9, 1.6, 3.0], cartesian, vacuum),
Quadrature.gauss_legendre(4), 2 groups)`` — and the operators come from
``sn.angular_bulk_space.retraction("angular")`` / ``.section("angular")``.
Inputs are ``numpy.random.default_rng(seed).standard_normal(shape)``.

Every entry below is a **bound over 200 independent draws**, not a
single reading: the residual of an exact-in-real-arithmetic identity is
a property of the numbers that happen to be involved, so one seed's
value is not reusable (:math:`\max_k \lVert a-b\rVert_\infty /
\lVert b\rVert_\infty` over ``default_rng(1000+k)``,
:math:`k = 0..199`).

.. list-table:: The square, measured (bound over 200 draws)
   :header-rows: 1
   :widths: 40 30 30

   * - Identity
     - synthetic 3-axis
     - :math:`S_N` carrier (GL4 slab)
   * - :math:`R\circ E = \mathrm{id}`
     - :math:`1.5\times10^{-16}`; ``array_equal`` on **123 of 200**
     - :math:`0.0` — ``array_equal`` on **200 of 200**
   * - :math:`P = E\circ R` idempotent
     - :math:`1.5\times10^{-16}`; ``array_equal`` on **130 of 200**
     - :math:`0.0` — ``array_equal`` on **200 of 200**
   * - :math:`R^{\dagger} = \pi^{*}` (the plain broadcast)
     - :math:`2.4\times10^{-16}`
     - :math:`2.3\times10^{-16}`
   * - :math:`R^{\dagger} = (\Sigma w)\,E`
     - :math:`2.4\times10^{-16}`
     - :math:`2.3\times10^{-16}`
   * - :math:`E^{\dagger} = R/\Sigma w`
     - :math:`3.6\times10^{-16}`
     - :math:`3.6\times10^{-16}`
   * - :math:`R\,R^{\dagger} = (\Sigma w)\,\mathrm{id}`
     - :math:`2.2\times10^{-16}`
     - :math:`2.2\times10^{-16}`
   * - the adjunction :math:`\langle R\psi,\varphi\rangle_{V'} = \langle\psi,\pi^{*}\varphi\rangle_V`
     - :math:`1.6\times10^{-13}`
     - :math:`3.7\times10^{-13}`
   * - :math:`P` self-adjoint in :math:`G_V`
     - :math:`4.0\times10^{-14}`
     - :math:`2.1\times10^{-14}`
   * - :math:`R \equiv` the shipped angular reduction
     - —
     - **bit-exact** (``np.array_equal``)
   * - :math:`E \equiv` the shipped isotropic-source kernel
     - —
     - **bit-exact** (``np.array_equal``)

⚠ The two **pairing** rows sit three orders above the others, and that
is arithmetic rather than a defect: an inner product of two random
fields can very nearly cancel, so the *relative* residual of a
scalar-valued identity is bounded by the conditioning of that
cancellation, not by the operators. Their gate row is written at
``rtol=1e-13`` for exactly this reason, and a *tighter* tolerance there
would be a latent false red (``vv-principles`` #16).

.. warning::

   **Do not read "bit-exact" as a law — it is a property of the draw.**
   :math:`R\circ E = \mathrm{id}` is exact in real arithmetic and holds
   at the round-off floor in IEEE-754; whether the floor is *zero*
   depends on how :math:`\sum_n w_n(\varphi/\Sigma w)` happens to
   re-associate for the particular numbers involved. `[M]` on the
   synthetic fixture — the one the gate uses — ``np.array_equal``
   FAILS on **844 of 2000** seeds (worst relative deviation
   :math:`1.5\times10^{-16}`, i.e. about one ULP), and the idempotence
   row fails on **57 of 200**. The originally-shipped G6.1/G6.2 rows
   pinned ``array_equal`` on seeds that happened to land in the exact
   set — seed-fragile — and were re-pinned at ``nulp=1`` the same day
   this audit measured the fragility (their docstrings carry the
   sweep). On the shipped :math:`S_N` carrier the identity
   is bit-exact on 200 of 200 seeds — there :math:`\Sigma w = 2`
   exactly *and* the symmetric Gauss–Legendre weights re-associate
   cleanly — which is why the production-facing rows can be pinned at
   ``np.array_equal`` honestly. Two consequences worth carrying: a
   fixture whose weights are chosen "non-uniform so no cancellation
   flatters a law" buys angular discrimination at the cost of exact
   re-association, so the general tier there is
   ``assert_array_almost_equal_nulp(..., nulp=4)``; and a multi-dim
   axis is one ULP by construction, because a flattened 2-D measure
   sums more terms.

**The divisor is the Gram entry, and that is stronger than
"the divisor is** ``weights.sum()``\ **".** The two agree on most
fixtures and not on all — the Gram is an einsum reduction
(``einsum("n,nj,nk->jk", w, table, table)``) and ``ndarray.sum`` is a
pairwise reduction, so they can differ by a ULP:

.. list-table:: Divisor vs. the naive total, on shipped quadratures
   :header-rows: 1
   :widths: 16 32 32 20

   * - rule
     - divisor (frame Gram entry)
     - ``quad.weights.sum()``
     - identical?
   * - ``gauss_legendre(4)``
     - ``2.0``
     - ``2.0``
     - yes
   * - ``gauss_legendre(8)``
     - ``1.9999999999999998``
     - ``2.0``
     - **no** (1 ULP — ``nextafter(2.0, 0)``; the gate's GL8 row pins
       the bound at ``nulp=1``)
   * - ``gauss_legendre(16)``
     - ``2.0``
     - ``2.0``
     - yes
   * - ``gauss_legendre(32)``
     - ``2.0``
     - ``2.0``
     - yes
   * - ``gauss_legendre(64)``
     - ``2.0000000000000004``
     - ``2.0000000000000004``
     - yes

At ``gauss_legendre(8)`` the shipped
:meth:`AngularSourceSink.from_isotropic
<orpheus.transport.source_sinks.angular_source_sink.AngularSourceSink.from_isotropic>`
therefore differs from a hand-written :math:`Q/\Sigma w` by
:math:`2.0\times10^{-16}` relative — one ULP, and the *induced* value
is the principled one. The lesson for a future gate: pin the divisor
against ``frame.discrete_gram[0, 0]`` (exact, always), never against
``weights.sum()`` (a fixture-dependent coincidence).

**The section is the harmonic frame's isotropic column — at the
measure level.** A natural question is whether the collapse pair
duplicates the spherical-harmonic frame's :math:`\ell = 0` channel.
Measured on the :math:`S_N` carrier above with
``HarmonicFrame.from_galerkin(sn.quad.angular_frame(L))``, feeding a
moment field that is zero except at :math:`(\ell, m) = (0,0)`:

.. list-table:: The isotropic column against the section
   :header-rows: 1
   :widths: 26 18 28 28

   * - frame
     - Gram max off-diagonal
     - :math:`\text{face}^{\dagger}(e_0\varphi)` vs :math:`E\varphi`
     - :math:`\text{reconstruction}(e_0\varphi)/W` vs :math:`E\varphi`
   * - slab, :math:`L=1`
     - :math:`5.6\times10^{-17}`
     - :math:`5.6\times10^{-17}`
     - **0.0** (``array_equal``)
   * - sphere, :math:`L=1`
     - :math:`5.6\times10^{-17}`
     - :math:`1.1\times10^{-16}`
     - **0.0** (``array_equal``)
   * - slab, :math:`L=2`
     - :math:`1.155`
     - :math:`16.17`
     - **0.0** (``array_equal``)
   * - sphere, :math:`L=2`
     - :math:`1.155`
     - :math:`16.17`
     - **0.0** (``array_equal``)

Two readings, and the second is the one to carry:

#. The **adjoint** correspondence — the section is the isotropic column
   of the harmonic frame's *physical* adjoint — holds exactly when the
   measured Gram is DIAGONAL, i.e. when the Parseval metric exists at
   all (:ref:`frame-parseval-metric`). ⚠ The discriminator is the
   **Gram**, not the geometry: the :math:`L=2` rows read identically on
   slab and sphere, because the angular frame is built from
   ``sn.quad`` and knows nothing about the spatial coordinate system.
   A 1-D polar Gauss–Legendre rule has no azimuthal nodes, so the
   :math:`m \ne 0` modes are not orthogonal under it and the Gram is
   dense at :math:`L\ge2` — `[M]` refining the polar order does not fix
   it (``gauss_legendre(8)`` at :math:`L=2` reads the same
   :math:`1.155` off-diagonal and the same :math:`16.17`).
#. The **measure-level** correspondence —
   :math:`E = \text{reconstruction}(e_0\,\cdot)/W` — is bit-exact in
   *every* configuration, dense Gram included, because it never touches
   a metric. That is the honest statement of the relationship, and it
   is also the reason the collapse pair is minted from an **indicator**
   frame rather than lifted out of the harmonic one: the collapse rides
   the measure, and it must keep working where the harmonic frame's
   metric does not exist.


.. _spaces-collapse-pair-gates:

Verification — cite the gate, never copy its numbers
------------------------------------------------------

The battery is ``tests/numerics/test_axis_marginal.py``, ``foundation``
-tagged throughout: these are software and mathematical invariants of a
*construction*, not equation claims, so no row carries
``verifies(...)``. Each row's docstring names the mutation that reddens
it.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Row family
     - What it pins
   * - ``TestSectionLaws``
     - :math:`R\circ E = \mathrm{id}`, idempotence of :math:`P`, and
       that the marginal space is the remaining axes **verbatim**
       (measures intact, so the marginal's metric stays physical).
   * - ``TestAdjointPairing``
     - the adjunction on the physical metrics, **with** its vv #19
       negative leg — the same pairing under a deliberately stripped
       spatial measure must break at O(1), because a positive reading
       alone cannot discriminate metric-loaded from metric-blind.
   * - ``TestTwoArrows``
     - :math:`R^{\dagger} = (\Sigma w)E` (the anti-ERR-051 row) and the
       *type* discrimination — ``R.H`` must not be an
       :class:`~orpheus.numerics.operator.AxisSectionOperator`.
   * - ``TestShippedKernelEquivalence``
     - bit-identity with the canonical angular reduction and the
       isotropic-source kernel on the real :math:`S_N` carrier, and
       that the angular marginal **is** the carrier's scalar bulk
       space. ⚠ Both equivalence rows are pinned against kernels
       hand-spelled **in the test**: since S6.2 the production targets
       route through the very operators under test, so a production
       comparison would be tautological.
   * - ``TestAxisGeneric``
     - the verbs are not angular-only — an untyped axis is admitted
       whatever its label, and a multi-dimensional axis contracts all
       of its dims with its own measure.
   * - ``TestAdmission``
     - every refusal in the clause table above, plus the shape guards
       in both directions, and the zero-total-weight asymmetry.
   * - ``TestFrameInduction``
     - the generator discipline itself: **tightness** (the minted
       kernels against an independently rebuilt literal frame's face
       contents), the **gram-derivation** of the divisor, the clause-2
       energy refusal, and that both verbs share one memoized mint.

.. _spaces-fences:

What is NOT built (the CS2 and CS4 seams)
=========================================

Stated explicitly so no reader mistakes a chartered design for shipped
machinery, and so the next phase does not re-derive a decision already
taken.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Not built
     - Where it lands, and what stands in for it today
   * - ``SpatialAxis`` / ``QuadratureAxis`` / ``HarmonicAxis``
     - **CS2.** There are no axis subclasses beyond
       :class:`~orpheus.numerics.axis.EnergyAxis`; the spatial factor of
       :attr:`MaterialMesh.bulk_space
       <orpheus.transport.mesh.material_mesh.MaterialMesh.bulk_space>`
       is a generic :class:`~orpheus.numerics.axis.Axis` labelled
       ``"spatial"``. The quotient point is a generic instance today and
       gets re-homed as ``SpatialAxis.quotient_point()`` when the
       subclass lands.
   * - Axis-built COMPOSITE and TRACE spaces
     - **CS2.** ⚠ Re-measured 2026-08-24: the *bulk* half of this
       fence has fallen. Campaign 1 CS4b moved the angular family onto
       axis-built carrier mints, so `[M]` on a shipped 1-D
       :class:`~orpheus.sn.mesh.augmented_mesh.SNMesh` the scalar bulk
       ``(energy, spatial)``, the angular bulk
       ``(angular, energy, spatial)`` and the scheme-widened
       ``angular_trial_space`` are ALL axis-built and all report
       ``has_coordinate_cone is True``. What is still legacy
       (``axes is None``, ``has_coordinate_cone is None``) is the
       **composite** :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
       and the flat **trace** buffers — the block/direct-sum structure
       the axis layer does not yet compose (see the
       :math:`\oplus` row below). This is also why the axis collapse
       pair refuses a non-axis-built space: it has no named factors to
       marginalise over.
   * - :math:`\oplus` composition (direct sums of spaces)
     - **CS2's opener.** The axis layer composes with :math:`\otimes`
       only; block/composite structure still rides
       :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
       and the coupled-block machinery.
   * - The identity flip (compare the axes tuple, not the name)
     - **S3.** Until then the derived name is the bridge
       (:ref:`spaces-identity-bridge`), and ``axes`` is declared
       ``compare=False``.
   * - Retiring the densifying ``__mul__``
     - **CS2.** The legacy ``*`` path is documented above and kept
       working; its own gates live in a separate test module so the
       retirement is a file-level move.
   * - The condensation morphisms on :math:`V` / :math:`V^*`
     - **Campaign 2.** Declared at
       :ref:`spaces-vv-collapse-hook`; the axis records the group
       structure they will consume. (This row said "S7" until
       2026-08-24 — a colliding plan-internal step number, see the
       warning at that anchor.)
   * - ``Medium`` and the mesh-conformity guard
     - **CS1.5.** See :ref:`spaces-symmetry-monotonicity`.
   * - Making the operators' ``space`` slot mandatory
     - **CS4.** Every operator's ``space`` is still
       ``FunctionSpace | None``; a ``None`` operand makes the
       composition guard skip rather than validate.


.. _spaces-development-history:

Development history
===================

Reverse-chronological (latest first) changelog of the architectural
milestones of *this page's* subject — the space layer. The field
algebra's own changelog is on
:doc:`/theory/foundations/field_algebra`, the operator algebra's on
:doc:`/theory/foundations/operator_algebra`, and the S\ :sub:`N`
solver's is :ref:`sn-development-history`. Entries marked *(in
development)* live on an unmerged feature branch and have no landed
merge-to-``main`` hash yet; trust ``git`` over this table for merge
status.

.. list-table::
   :header-rows: 1
   :widths: 10 50 12 28

   * - When
     - Architectural milestone
     - Issue
     - Where
   * - 2026-08-24
     - **The space becomes the construction key, and it mints the
       collapse pair** (campaign 1, phase CS4b, steps S5–S7).
       *Construction goes space-primary* (S5): the carrier gains
       :attr:`SNMesh.angular_trial_space
       <orpheus.sn.mesh.augmented_mesh.SNMesh.angular_trial_space>`
       (the scheme-widened angular mint), the composite allocators go
       space-keyed, and the mesh-keyed leaf **sugar tier is deleted** —
       every call site now names a space, not a carrier.
       *The axis collapse pair is minted on the space* (S6): the verbs
       :meth:`FunctionSpace.retraction
       <orpheus.numerics.space.FunctionSpace.retraction>` /
       :meth:`~orpheus.numerics.space.FunctionSpace.section` return the
       split epi/mono pair, memoized per axis label, and S6.0b re-carved
       their REALIZATION so both are the induced output of a
       single-region indicator frame built and discarded at one site —
       the **stage-2 generator discipline**, with the section's divisor
       read off the frame's :math:`1\times1` Parseval metric rather
       than chosen by hand (:ref:`spaces-collapse-pair`). The angular
       reduction and the isotropic-source projection are re-keyed onto
       that pair, so each has ONE realization tree-wide, and
       ``AxisEmbeddingOperator`` is renamed
       :class:`~orpheus.numerics.operator.AxisSectionOperator` on
       canonical-naming grounds (:ref:`spaces-collapse-pair-naming`).
       *The mesh-less carrier's two meanings un-weld* (S7): promoting
       the infinite-medium 1-cell carrier to an :math:`S_N` phase space
       raises a typed :class:`ValueError` (pre-repair a messageless
       bare ``assert`` that ``-O`` stripped into a deep
       ``AttributeError``), and
       :attr:`MaterialMesh.areas
       <orpheus.transport.mesh.material_mesh.MaterialMesh.areas>` names
       its own three cases instead of blaming a 2-D mesh for all of
       them. The homogeneous solver's reaction rates become the typed
       integrated co-vector, **re-posed onto the solver's own pose** so
       the pose stays the measure authority
       (:doc:`/theory/foundations/infinite_medium`).
     - —
     - merged @ ``55bb47b9`` —
       ``b00bf2d7`` … ``2690a434`` (space-primary construction),
       ``048144db`` (the pair), ``19b85775`` (the frame induction +
       the rename), ``ffb8f286`` (space-derived truncation),
       ``78925753`` / ``53e7d207`` (the re-keys and the packer
       re-home), ``1f8e0323`` (the S7 repairs), ``2e054bfc`` (the
       typed rate co-vector)
   * - 2026-08-20
     - **The space layer gains AXES, and the energy axis is the first
       one** (campaign 1, phase CS1). A new
       :mod:`orpheus.numerics.axis` mints
       :class:`~orpheus.numerics.axis.Axis` (frozen; structural
       per-subclass identity; canonical measure storage — all-ones
       collapses to ``None``, :math:`-0.0` normalized, non-finite
       refused, signed weights deliberately legal) and
       :class:`~orpheus.numerics.axis.EnergyAxis`
       (``from_grid`` / ``synthetic``; identity = :math:`n_g` + edges
       CONTENT; weighted axes refused by the counting-measure theorem).
       :meth:`FunctionSpace.of_axes
       <orpheus.numerics.space.FunctionSpace.of_axes>` composes a space
       as the ordered product of its axes with a **per-axis metric
       path** (no densification) and a deterministic, injective
       **derived name** — the identity bridge that makes "metric
       differences imply space differences" true today
       (:ref:`spaces-identity-bridge`).
       :attr:`FunctionSpace.has_coordinate_cone
       <orpheus.numerics.space.FunctionSpace.has_coordinate_cone>` makes
       the NODAL/MODAL dichotomy machine-readable, and
       :meth:`Field.cone_violations
       <orpheus.numerics.field.Field.cone_violations>` consults it.
       :attr:`MaterialMesh.bulk_space
       <orpheus.transport.mesh.material_mesh.MaterialMesh.bulk_space>`
       mints the scalar bulk through ONE uniform formula, so the
       degenerate carrier's quotient point, a genuine one-cell mesh and
       a meshed carrier all fall out of the same body; the homogeneous
       solver then poses :math:`A = C - K_{\rm iso}` and
       :math:`F` on that real space, retiring both hand-written
       **production** ``basis_shape=(ng, 1)`` spellings and turning the
       ``OperatorSum`` guard from *skipped* into *validating* — `[M]`
       zero ``basis_shape=(ng, 1)`` call sites remain anywhere in
       ``orpheus/``, though the keyword survives as an explicit
       override. The four
       ``test_monomorphic_leaves`` strict-xfail rows are deleted and
       succeeded by a positive floor. The **counting-measure theorem is
       why this moved no values**: identity metrics along both factors,
       and guards that compare spaces rather than values.
     - —
     - merged @ ``55bb47b9`` —
       ``1afff47b`` (the axis), ``f4876354`` (``of_axes`` + per-axis
       metric + cone metadata), ``e8769897`` / ``24a991ba`` (the
       operators' space slot renamed and widened), ``6bd782ab`` (the
       homogeneous posing), ``6da1b23c`` (the cone consult)
   * - 2026-08-19
     - **Flux is typed into the positive cone** :math:`K \subset V`,
       which is what made *this* page's question askable: once cone
       membership is an element predicate rather than a constructor
       invariant, "on which spaces is the predicate meaningful?" becomes
       a question about the SPACE, and the answer is the basis kind of
       its factors (:ref:`spaces-nodal-modal`). See
       :doc:`/theory/foundations/field_algebra`.
     - #331
     - merged ``f9d571b5``

.. admonition:: Verification — cite the gate, never copy its numbers
   :class: note

   The CS1 battery is ``foundation``-tagged throughout: these are
   software and mathematical invariants of a *type*, not equation
   claims, so no gate carries ``verifies(...)``.

   - ``tests/numerics/test_axis.py`` — the intrinsic laws of the axis
     concept: rank, measure canonicalization, the refusals, structural
     identity per subclass, and the ``synthetic`` / ``from_grid``
     inequality.
   - ``tests/numerics/test_space_of_axes.py`` — composition: shape
     concatenation, the per-axis metric against an independently built
     reference, the no-densification proof, the derived name's
     determinism across processes, and ``has_coordinate_cone``.
   - ``tests/numerics/test_field.py`` (gates E1/E2) — the cone
     consult's **positive and negative pair**: a MODAL space REFUSES
     with a typed error naming the space, and the same values on an
     all-NODAL space answer exactly what the legacy path answers.
   - ``tests/homogeneous/test_operator_spaces.py`` — the positive floor
     (all five homogeneous operators plus :math:`K` report the SAME
     space), the refusal witnesses (a 2g-vs-4g sum; :math:`M^{-1}(2g)
     \circ F(4g)`), the energy arm's ``from_grid``-vs-``synthetic``
     discrimination, and the ``.H`` **loaded/blind pair** required by
     ``vv-principles`` anti-pattern #19 — a bit-identity leg for the
     shipped scalar-metric case, paired with a deliberately
     non-physical per-group-weighted axis on which ``.H`` demonstrably
     MOVES.
   - ``tests/homogeneous/test_byte_stability.py`` — the migration gate
     that measured the theorem. It pins the homogeneous solve
     bit-exactly (``np.array_equal`` and exact ``==``, never
     ``allclose``) against a baseline captured immediately before the
     wiring, over every producing mixture the tree ships. `[M]` it held
     bit-exactly across the rewiring on 2026-08-20 — which is the
     evidence for the "no value motion" claim above. It is a CS1
     migration gate by design and retires after the merge cycle,
     subsumed by the L1 correctness anchor and the materialization byte
     pin.
