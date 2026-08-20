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
   morphisms themselves are scheduled for S7 / campaign 2, and the
   condensation that ships today
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
   `[M]` **this is live, and it is measurable.**
   :func:`~orpheus.homogeneous.solver.solve_homogeneous_infinite`
   normalises its flux through
   :class:`~orpheus.transport.reaction_rate_functional.IntegratedReactionRate`,
   whose ``evaluate`` contracts against ``mesh.volume_measure`` — so
   the quotient point's weight is genuinely consumed, not decorative.
   Feed the SAME flux :math:`\varphi = (1, 1)` and the same 2-group
   mixture to the quotient carrier (weights ``[1.0]``) and to a
   one-cell slab of width 2 (weights ``[2.0]``) and the production
   rates are :math:`0.225` and :math:`0.450` — exactly the measure
   ratio, with no other input changed.
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
   * - Axis-built S\ :sub:`N`, diffusion and P\ :sub:`N` spaces
     - **CS2.** Only the scalar bulk of a
       :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` is
       axis-built today, and its only consumer is the homogeneous
       solver. Every other space in the tree is legacy
       (``axes is None``) and therefore reports
       ``has_coordinate_cone is None``.
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
     - **S7 / campaign 2.** Declared at
       :ref:`spaces-vv-collapse-hook`; the axis records the group
       structure they will consume.
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
     - *(in development)* branch ``feature/cs1-energy-space`` —
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
