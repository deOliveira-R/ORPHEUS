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
      role: "the space layer — a function space as the ordered product of its AXES (label, index shape, factor measure, kind, support, atoms and coefficient layout, all derived from the axis's ORIGIN, which is provenance and excluded from identity), one identity key per space, the counting-measure theorem on the energy axis, and the collapse doctrine that decides which axes survive a degeneracy and why"
      depends_on: [field_algebra, frame]
      related: [manifolds, discrete_measures, operator_algebra, operator_adjoint]
      status: "seeded at campaign-1 CS1 (the Energy axis); the generator slot landed at CS5 (2026-08-29); the densifying half of the `*` product retired at CS4c step 6 item 6.2a (2026-09-07); the metric became a DERIVED object over the axes, with a dense form positioned on it as an overlay, at item 6.2c-i (2026-09-08); the Harmonic and Legendre axis subclasses landed at item 6.2c-ii (2026-09-08), so both moment heads are axis-built; the spatial-moment tail became the discretization scheme's own MODAL moment axis at item 6.2c-iii (2026-09-08), so every space that carries a tail is axis-built too; the Spatial / Quadrature subclasses stay CS2 and unscheduled; the axis gained its unit, and the marginal over a MODAL head became an integral, at S1 of the numerics untangling (2026-10-19); at S7 (2026-10-22) the generator slot became the Origin closed sum, the axis content gained support, atoms and layout, a space's identity became one key, the two head classes and the two family axis subclasses retired onto one generic coefficient axis carrying the basis's layout, and every basis's coefficient space became axis-built (FunctionSpace.of_basis)"


This page develops the **space layer**: what a discrete function space
*is*, what it is made of, and what happens to it when a physical
degeneracy collapses one of its factors. It is the companion to
:doc:`/theory/foundations/field_algebra`, which types the **elements**
(a flux is a point of the positive cone :math:`K \subset V`); this page
types :math:`V` itself.

The organizing claim is one sentence: **a function space is the ordered
product of its axes, and an axis is identified by what fixes its
factor's space — its label, its index shape, its factor measure, its
kind (with the coordinates of the constant function, its unit), the
manifold its coordinates live on, the atoms a nodal factor samples and
the coefficient layout of a basis's factor — all of it read off the
axis's ORIGIN, the object that minted it, which the identity
deliberately EXCLUDES.** Everything else on this page follows from taking that
seriously: why the energy metric is the identity as a *theorem* rather
than a default, why the homogeneous solver's spatial factor survives its
own collapse while the angular factor of a scalar space does not, why an
axis can hand a consumer back the direction cosines it forgot, and why
"two copies of :math:`\mathbb{R}^n` with different inner products are
the same space" is a claim this corpus has **overturned**.

.. warning::

   **Two things are easy to run together: what an axis IS and what
   MINTED it.** The axis's content is typed per subclass (an
   :class:`~orpheus.numerics.axis.EnergyAxis` is not a generic
   :class:`~orpheus.numerics.axis.Axis`) and holds everything its origin
   implies about the factor's space: the kind, the support, the atoms,
   the layout. The origin itself, the INSTANCE that minted the axis, is
   provenance and is never identity (:ref:`spaces-axis-origin`).

.. note::

   **Three unrelated things in this corpus are called an "axis". Keep
   them apart.**

   - A **space-factor axis** (:class:`~orpheus.numerics.axis.Axis`,
     :class:`~orpheus.numerics.axis.EnergyAxis`) — *this page's*
     subject. One tensor factor of a function space; it carries a
     measure, a kind, and, for a factor sampled from a measure, the
     atoms and the manifold they lie on. The rank-1 spatial factor's
     origin is the carrier's volume measure, so its atoms are the cell
     centres: the cell positions are part of the space (since S7,
     ruling 71), while the mesh's edges, faces and coordinate system are
     not (:ref:`spaces-axis-origin`).
   - A **geometric axis** (``Axis1D`` / ``AxisMesh`` /
     ``RadialAxisMesh`` in :mod:`orpheus.mesh.axis`) — one
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
   - ⚠ **…and a space is NOT a domain.** There are three levels — the
     manifold :math:`M`, the fields on it :math:`L^2(M)`, and the
     coefficients :math:`\mathbb{R}^K` — and a
     :class:`~orpheus.numerics.space.FunctionSpace` is the **second**.
     A basis function eats a *point* of :math:`M`, so this page's
     object can never be a basis's domain, and the frame's level-2
     arrow type-checking says nothing about the level-1 pairing (which
     is how :ref:`ERR-080 <manifold-err-080>` survived). Level 1 was a
     bare ``str`` until 2026-08-31 and now has its own type and its own
     page: :doc:`/theory/foundations/manifolds`.
   - **An axis is identified by its content** *(label, shape, weights,
     kind, support, nodes, layout)*; the kind, the support, the atoms and
     the layout are read off its origin, never passed
     (:ref:`spaces-axis-identity`).
     ``weights=None`` **is** the counting measure — deliberately and
     always; an axis has no "unbound" state, so the legacy two-state
     ambiguity of ``inner_product_weights`` cannot arise on this type.
   - **The metric is an OBJECT, and a space has exactly ONE source for
     it.** Since campaign 1 P7 a space may carry a typed
     :class:`~orpheus.numerics.metric.HilbertMetric` — the Hadamard
     weight is its diagonal special case — so a form with off-diagonal
     structure is expressible at last; resolution is
     ``metric`` > ``inner_product_weights`` > Euclidean, with axis-built
     spaces routing through their axes instead, and a construction guard
     refuses any two sources at once. ⚠ Consequently
     ``inner_product_weights is None`` **no longer implies Euclidean** —
     read a space's pairing through
     :meth:`~orpheus.numerics.space.FunctionSpace.inner_product`, never
     through the slot (:ref:`spaces-metric-object`).
   - **An axis records its ORIGIN, a closed sum, and the mint is a
     section.** :attr:`Axis.origin <orpheus.numerics.axis.Axis.origin>` is
     one of ``FromMeasure``, ``FromQuadrature``, ``FromBasis``,
     ``FromFrame`` or ``Declared`` (ruling 64), read by an exhaustive
     ``match``; re-minting from the recorded object gives the axis back
     (:eq:`spaces-axis-generator-section`), so a consumer holding the
     SPACE recovers ``mu_x`` / ``eta`` / ``mu_z`` / ``level_indices``
     through ``case FromQuadrature(quad)`` without being handed the
     quadrature separately, and one refusal,
     :func:`~orpheus.numerics.axis.quadrature_refusal`, names the axis, the
     asker and the arm when the origin is another (:ref:`spaces-axis-origin`).
   - ⛔ **The origin is provenance and stays out of identity — a
     RULING.** *Two axes with the same content are the same axis whatever
     object minted them*: ``origin`` is declared ``field(compare=False)``,
     the content encoder's one spelling of "not content"
     (:ref:`structured-geometry-content-identity`), while everything the
     origin implies (kind, support, atoms, layout) is stored on the axis
     and compared. ``Axis.of_quadrature(q)`` and
     ``Axis.of_measure(q.measure, "angular")`` are one axis
     (:ref:`spaces-generator-identity-exclusion`).
   - **The two one-line discriminators of the collapse doctrine.**
     *(i)* **Can the admissible fields be integrated over the collapsed
     domain?** — symmetry-forced constancy on infinite-measure orbits
     means no, so the collapse must NORMALIZE. *(ii)* **Is the
     surviving convention consulted INSIDE the family, or only at
     re-embedding?** — inside ⟹ the axis persists; only at re-embedding
     ⟹ the convention lives on the arrow and the axis DROPS.
   - **A collapse that DROPS an axis is realized by typed arrows built
     from the axis's unit, never by one verb.**
     :meth:`AxisRetractionOperator.of(space, label)
     <orpheus.numerics.axis_marginal.AxisRetractionOperator.of>` is fibre
     integration :math:`R = m^{\mathsf T}\otimes\mathrm{id}`, with the
     covector :math:`m = G_k u` (:math:`u` the coordinates of the constant
     function :math:`1` on the axis, :attr:`Axis.unit
     <orpheus.numerics.axis.Axis.unit>`), so a ``MODAL`` head marginalises
     to the integral of the function it represents;
     :meth:`AxisSectionOperator.of(space, label)
     <orpheus.numerics.axis_marginal.AxisSectionOperator.of>` is its right
     inverse :math:`E = (u/\langle m,u\rangle)\otimes\mathrm{id}` with
     :math:`R\circ E = \mathrm{id}`, and ``R.H`` is the pullback
     :math:`\pi^{*} = u\otimes\mathrm{id}`. The adjoint and the section
     differ by exactly the axis's mass (:math:`R^\dagger = \langle
     m,u\rangle\,E`), and they carry different TYPES so that scalar
     cannot be dropped at a call site. All three are minted at one site
     and memoized on the space (:ref:`spaces-collapse-pair`).
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
     it. A genuine one-cell mesh is a **different** space, with
     :math:`V_{\rm cell} \neq 1` by its measure and with
     :math:`V_{\rm cell} = 1` by its atom: the mesh's spatial axis keeps
     its cell centre on :math:`\mathbb{R}`, the counting point has none.
     The two are isometric and unequal, and the collapse
     :math:`\mathbb{R} \to \{*\}` is a map, a pushforward, never an
     equality (ruling 85, :ref:`spaces-axis-identity`).
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
   - **A space's identity is ONE key**,
     :meth:`FunctionSpace.identity_key()
     <orpheus.numerics.space.FunctionSpace.identity_key>`, and ``==`` and
     ``hash`` derive from it in the base class only (S7, #631). An
     axis-built space is its axes plus the content of every dense form
     positioned over them (the matrix; its inverse face is derived, not
     compared); an axes-less leaf is its class, name and shape; a dual, a
     non-axis tensor product, a full-field composite and a coupled system
     are their members' keys, tagged. So **metric differences, atom
     differences and domain differences all imply space differences**,
     and an axis-built space is never equal to a hand-named one wearing
     its label. `[M]` 2026-10-22, the six collisions #631 and its review
     measured each separate now: among 15 shipped rules on :math:`S^2`,
     0 of 105 measure-space pairs compare equal, where 21 did
     (:ref:`spaces-identity-key`).

.. _spaces-the-axis:

The axis: its slots, and what each one decides
==============================================

An **axis** is one tensor factor of a function space — the value object
recording what fixes that factor's coordinate space, together with its
**origin**, the object that minted it. It is the unit the composition
machinery reasons about: partitions, collapses, frames and (later)
:math:`\oplus`-lifts act **per axis**, never on an anonymous position of a
monolithic shape tuple. Three slots are passed (``label``, ``shape``,
``weights``), one is passed and is provenance (``origin``), and four are
read off the origin at construction (``kind``, ``support``, ``nodes``,
``layout``); the unit is a property of the kind.

.. list-table:: The slots of :class:`~orpheus.numerics.axis.Axis`
   :header-rows: 1
   :widths: 14 30 56

   * - Slot
     - Type
     - What it decides
   * - ``label``
     - ``str``
     - The factor's role (``"angular"``, ``"spatial"``, ``"energy"``, a
       coefficient layout's own label such as ``"harmonic"``). Content: two
       factors that agree on everything else and play different roles are
       different axes.
   * - ``shape``
     - ``tuple[int, ...]``, rank :math:`\ge 1`
     - The factor's index set. Rank :math:`> 1` is admissible and
       deliberate: a spherical-harmonic head is :math:`(L+1, 2L+1)`, and
       a rank-:math:`d` spatial axis is a legal CS2 design choice. Rank
       0 is refused — a factor with no index set has nothing to measure.
   * - ``weights``
     - ``NDArray | None`` over exactly ``shape``
     - The **factor measure**. ``None`` IS the counting measure. This is
       the slot the pairing consumes. A raw partition head carries
       ``None`` because a partition is measure-free until a frame binds a
       fine measure (ruling 88, :ref:`spaces-collapse-pair-frame`).
   * - ``origin``
     - ``FromMeasure | FromQuadrature | FromBasis | FromFrame | Declared``,
       keyword-only, no default
     - **Provenance, never identity** (ruling 64): the object that minted
       the factor, read by exhaustive ``match``. A consumer recovers through
       it what the axis does not keep (a quadrature's level fibration, the
       frame that can re-dress a head). :ref:`spaces-axis-origin` develops
       it.
   * - ``kind``
     - ``BasisKind.NODAL`` or :class:`Modal(unit) <orpheus.numerics.axis.Modal>`,
       a closed sum; read from the origin
     - ``NODAL`` ⟹ the factor carries a **coordinate cone**; the modal
       arm ⟹ it does not, and it CARRIES the factor's unit. Not a
       constructor parameter: ``Axis(kind=NODAL, generator=basis)`` cannot
       be spelled.
   * - ``unit``
     - ``NDArray`` over exactly ``shape``, read-only; a property read from
       the kind (``kind.unit_over(shape)``)
     - The **coordinates of the constant function** :math:`1` in this
       factor's representation (ruling 37 and ruling 41 of the numerics
       untangling): all ones on a ``NODAL`` axis; on a modal axis the unit
       the ``Modal`` arm carries: :math:`e_0` on a continuum head, the
       moments :math:`\Phi^{\mathsf T} w` on a head a frame has dressed,
       :math:`e_0` on the LD spatial-moment axis, and a declared zero on
       the loss-kernel gauge's modes (rulings 86 and 87). The marginal over
       the axis integrates against :math:`m = G_k u` (:ref:`spaces-axis-unit`).
   * - ``support``
     - a :class:`~orpheus.geometry.manifold.Manifold`, or ``None``; read
       from the origin
     - The manifold the factor's coordinates live on: a measure's support
       (:math:`S^2`, the orbit space :math:`S^2/O(2)_x`, a real line), a
       basis's domain. ``None`` on a ``Declared`` axis. Content (ruling 71).
   * - ``nodes``
     - read-only ``NDArray`` or ``None``; read from the origin
     - The **atoms** a nodal factor samples (a quadrature's directions, a
       carrier's cell centres). ``None`` on a coefficient or ``Declared``
       axis. Content (ruling 71).
   * - ``layout``
     - a coefficient layout or ``None``; read from the origin
     - How a basis indexes its coefficients
       (:mod:`orpheus.numerics.coefficient_layout`): the rectangular
       :class:`~orpheus.numerics.coefficient_layout.HarmonicLayout`, the
       flat :class:`~orpheus.numerics.coefficient_layout.LegendreLayout`,
       the partition's :class:`~orpheus.numerics.coefficient_layout.CellLayout`
       (its edges are content), or a basis computed from data's
       :class:`~orpheus.numerics.coefficient_layout.ModeLayout`. ``None`` on
       a sampled or declared axis. Content (ruling 81); an axis refuses a
       layout whose shape is not its own.
   * - identity
     - structural, **per subclass**
     - :math:`(\text{type}, \text{label}, \text{shape}, \text{weights},
       \text{kind}, \text{support}, \text{nodes}, \text{layout})` plus each
       subclass's own fields (an ``EnergyAxis``'s edges), digested by the
       one content encoder; never the origin. :ref:`spaces-axis-identity`
       states what each part decides and the collisions it removed.

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
   (``tests/gates/numerics/test_space_of_axes.py``: shape concatenation, the
   per-axis metric against an independently-built dense reference, the
   no-densification proof, and the derived name's determinism across
   processes) plus ``tests/gates/numerics/test_axis.py`` for the factor laws.
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

   **The legacy twin survives; its DENSIFYING half retired 2026-09-07**
   (CS4c step 6 item 6.2a). ``V * W`` (the ``*`` dunder →
   :class:`~orpheus.numerics.space.TensorProductSpace`) is still the
   PRE-axis composition mechanism — it takes whole *spaces* as factors
   and derives its name by joining theirs — but it no longer builds a
   dense weight array. It never populates ``inner_product_weights`` at
   all. Two arms:

   * **Every factor axis-built** → the product's ``axes`` is the
     concatenation of the factors' axes, and there is no metric object:
     the per-axis path of :eq:`spaces-axis-product` applies unchanged.
     (This arm is untouched by 6.2a.)
   * **Any factor axes-less** → the product's ``axes`` stays ``None``
     (an axis is never fabricated for a space that did not declare one)
     and the metric is a lazy
     :class:`~orpheus.numerics.metric.FactoredMetric` carrying ONE
     positioned entry per **axis** of an axis-built factor (a
     :class:`~orpheus.numerics.metric.DiagonalMetric` on that axis's
     block, or ``None`` for a counting-measure axis), one entry per
     dense-slot leaf factor, and a factor's own metric object riding
     verbatim — a nested ``FactoredMetric`` flattens, a
     :class:`~orpheus.numerics.metric.DenseMetric` positions on its
     block. Every block Euclidean ⟹ no metric at all.

   ⛔ **This block read** *"The legacy twin is still live, and CS2
   retires it … on all-diagonal factors it DENSIFIES the metric into an
   outer-product* ``inner_product_weights`` *… CS1 keeps it — it threads
   the* ``axes`` *record when both sides carry one, and bridges
   axis-borne measures into its dense weights on mixed products … and
   CS2 collapses the live mints onto axis concatenation and retires the
   densifier"* **until 2026-09-07**, and that was an accurate
   description of the tree until item 6.2a. The product formed
   :math:`w_V \otimes w_W` once and stored it; a *mixed* product (one
   axis-built factor, one not) therefore had to BRIDGE the axis-borne
   measure into that dense slot, because a weighted axis-built factor
   read through the dense slot alone would have been treated as
   Euclidean — a value bug, not a representation choice. The bridge and
   the outer-product builder retired together with the slot they fed.

   What survives of the old rule is the half that was never about the
   densifier: **new axis-aware code composes with** ``of_axes``, because
   ``of_axes`` is the only ROOT producer of an ``axes`` record — ``*``
   and :meth:`dual <orpheus.numerics.space.FunctionSpace.dual>` merely
   thread one through, so both need an axis-built ancestor (the same
   closure argument the nodal/modal refusal rests on, below).
   ``*`` is for a product whose factors are **not** all
   axis-built. `[M]` 2026-09-07, by AST over ``orpheus/**/*.py``, that is
   the harmonic/moment family and nothing else — **four** production
   sites: the angular head :math:`\otimes` the axis-built cell group, in
   :meth:`HarmonicFrame.moment_space_on
   <orpheus.transport.frames.harmonic_frame.HarmonicFrame.moment_space_on>`
   and in the SN carrier's content-equal field-side twin
   :meth:`SNProblem.moment_space
   <orpheus.sn.problem.SNProblem.moment_space>`, each of which
   ALSO appended, for a widened angular space, a second ``*`` factor —
   the Euclidean, axes-less ``SpatialMomentSpace`` (contributing a
   ``None`` metric entry). Item **6.2c** axis-ifies the head; until then
   see :ref:`spaces-metric-propagation`, where dropping a factor's
   metric is measured as a value bug rather than a representation
   choice.

   ⚠ **The field-side member of those four MOVED later the same day
   (2026-09-07, CS4c step 6 item 6.2b), and the count did not.** Until
   that item the field-side twin was
   :meth:`MomentField._space_for_mesh_and_L
   <orpheus.transport.fields._bases.MomentField._space_for_mesh_and_L>`,
   which formed the product on **every call**; the item gives the mint to
   the carrier as a cache keyed on ``(L, spatial_moments)``, and the field
   method survives as a pure READ of it. So the census was the same four —
   the hub's head-:math:`\otimes`-cell and its tail, the frame's
   head-:math:`\otimes`-cell and its tail — while the number of ``*``
   *evaluations* per windowed solve falls from `[M]`
   :math:`2\cdot\texttt{max\_inner} + 6` to a constant in the iteration
   budget (gated in
   ``tests/gates/sn/mesh/test_hub_owns_the_moment_space.py``).

   ⛔ **The count is now TWO, not four (2026-09-08, CS4c step 6 item
   6.2c-iii).** The two TAIL products are gone: a widened space no
   longer appends a separate Euclidean class beside the head, it carries
   the discretization scheme's own MODAL ``moment_axis`` as one more
   axis of its axis-built cell group. The hub composes that tail through
   the fields' own composer,
   :meth:`BulkField.compose_spatial_moments
   <orpheus.transport.fields._bases.BulkField.compose_spatial_moments>`
   (which returns :meth:`FunctionSpace.of_axes
   <orpheus.numerics.space.FunctionSpace.of_axes>`, not a ``*``), and
   the frame simply keeps the angular space's own axes. `[M]` 2026-09-08,
   the same AST census over ``orpheus/**/*.py``: **two** production
   ``*`` sites remain — ``head * bulk`` at
   ``sn/problem.py:1313`` and ``self.basis_space *
   FunctionSpace.of_axes(*axes[1:])`` at
   ``transport/frames/harmonic_frame.py:499``.

   ⭐ And with 6.2c-ii and 6.2c-iii both landed, **neither remaining
   site is the mixed case the rule above describes any more**: the head
   is axis-built, the cell group is axis-built, so ``*`` takes its
   *every-factor-axis-built* arm and the product's ``axes`` are the
   factors' concatenated. `[M]` 2026-09-08, on a 1-D LD carrier under
   ``gauss_legendre(2)``, the moment product's axes read
   ``[legendre, energy, spatial]`` at width 1 and
   ``[legendre, energy, spatial, spatial_moment]`` at width 2, at
   :math:`L = 1` (diagonal Gram) and :math:`L = 2` (DENSE Gram) alike —
   four rows, all axis-built. So ``*`` here is no longer buying the
   mixed-product bridge; it is spelling a concatenation ``of_axes``
   would spell, and the two sites are a **retirement candidate** rather
   than a live counter-example to the ``of_axes``-first rule.
   The composer serves EVERY tailed mint now — angular, scalar and the
   carrier's moment product alike — and its axes-less arm is a caller
   error rather than a case, because a space that carries a tail is
   axis-built by construction.

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
outlaw them. `[M]` ``Axis("x", (2,), weights=[-1.0, 2.0],
origin=Declared(BasisKind.NODAL))`` constructs.

.. _spaces-axis-identity:

An axis is identified by what fixes its space (ruling 71)
---------------------------------------------------------

Equality and hashing read the structural content, never object identity
and never a subset of the fields: the digest of *(type, label, shape,
weights, kind, support, nodes, layout)* plus each subclass's own fields,
through the one content encoder
(:class:`~orpheus.values.content.ContentIdentity`). Each part is there
because some pair of factors differs in it and nothing else, and each such
pair is two spaces:

- a **nodal** axis is fixed by its support, its atoms, its measure and its
  label. Atoms alone do not suffice: `[M]` the same Gauss–Legendre
  :math:`\mu`-rule about two poles, ``gauss_legendre_on_polar_orbit(8,
  "x")`` and ``(8, "z")``, has bit-equal nodes and weights and lies on two
  orbit spaces, :math:`S^2/O(2)_x` and :math:`S^2/O(2)_z`, and only the
  support separates the two axes;
- a **modal** (coefficient) axis is fixed by its domain, its layout, its
  kind with its unit, its measure and its label. The domain separates the
  real harmonics on :math:`S^2` from their mirror-even restriction on
  :math:`S^2/\sigma_y`, which share the rectangular layout and the
  continuum Gram, and the Legendre heads about two axes, which share the
  flat layout;
- a **partition** head is fixed by its cell edges, carried by its
  :class:`~orpheus.numerics.coefficient_layout.CellLayout`. Its masses are
  not enough: edges :math:`[0,1,2]` and :math:`[1,2,3]` have equal widths,
  equal masses and one support, and are two coefficient spaces because the
  cells are the partition's atoms (the main agent's ruling of 2026-10-21,
  an instance of ruling 71);
- an :class:`~orpheus.numerics.axis.EnergyAxis` never equals a generic
  :class:`~orpheus.numerics.axis.Axis` carrying the same field tuple (the
  class is a part), and ``EnergyAxis.synthetic(2)`` never equals
  ``EnergyAxis.of_edges(edges_2g)``: same index set, no partition data, a
  different axis.

The **origin is not a part** (:ref:`spaces-generator-identity-exclusion`):
``Axis.of_quadrature(q)`` and ``Axis.of_measure(q.measure, "angular")`` are
one axis, and two content-equal rule instances mint one axis.

**A unit cell and the quotient point are distinct spaces; collapse is a
map** (ruling 85, the user, 2026-10-21). The homogeneous pose's spatial
factor is a ``Declared`` one-point counting axis with no atom; a one-cell
mesh of unit width is a ``FromMeasure`` axis whose atom is its cell centre
on :math:`\mathbb{R}`. `[M]` 2026-10-22 (``Axis("spatial", (1,),
origin=Declared(NODAL))`` against ``Axis.of_measure`` of a one-atom
measure at :math:`0.5` with weight :math:`1`): both carry the counting
measure (``weights`` is ``None`` on both after canonicalization), and they
compare unequal, through the atom. They are isometric, and the collapse
:math:`\mathbb{R} \to \{*\}` that relates them is a pushforward, an
operator, never an identity of spaces (:ref:`spaces-collapse-doctrine`).
The bridge gate that compares the two asserts exactly that: they share
their arms (one energy axis, one counting measure) and differ by the atom
(``tests/gates/homogeneous/test_operator_spaces.py``
``::test_minted_space_and_a_unit_cell_carriers_bulk_space_share_their_arms``).

.. _spaces-identity-key:

A space's identity is one key (#631)
------------------------------------

:meth:`FunctionSpace.identity_key()
<orpheus.numerics.space.FunctionSpace.identity_key>` is the one
definition of what a space is, and ``==`` and ``hash`` derive from it in
the base class only (S7, finding F7 of the S6a review); no member of the
family re-states equality. The key is a tagged tuple per kind of space:

.. list-table:: The identity key
   :header-rows: 1
   :widths: 30 70

   * - space
     - key
   * - axis-built (``axes`` is not ``None``), including an axis-built
       tensor product
     - ``("axes", axes, overlay)``: the axis tuple, and for each axis the
       content digest of the dense form positioned over it (``None`` where
       there is none). A dense form's content is its matrix; its inverse
       face is derived from the matrix and is not compared.
   * - axes-less leaf
     - ``("named", class, name, shape)``. The factory-built leaves (the
       two trace spaces, the two radial-characteristic spaces) fold their
       content into a digest name, and the class keeps a hand-built space
       wearing such a name from passing for one.
   * - :class:`~orpheus.numerics.space.DualSpace`
     - ``("dual", primal.identity_key())``: covariance is identity (ruling
       75), so a dual never equals a primal.
   * - :class:`~orpheus.numerics.space.TensorProductSpace` with an
       axes-less factor
     - ``("tensor", factor keys)``.
   * - :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
     - ``("full_field", interior key, trace key)``.
   * - :class:`~orpheus.numerics.coupled_system.CoupledSpace`
     - ``("coupled", member keys)``.

The members of a composite are required fields: a full-field space without
its blocks, a coupled space without its systems and a tensor product
without its factors cannot be constructed.

**The collisions the key removes.** Issue #631 measured four pairs of
spaces that compared equal while being different spaces; the design added a
fifth and the S6a review a sixth. Each separates now, beside a control pair
that must stay equal:

.. list-table:: `[M]` before (``main`` at ``dd17f678``, the S7 census and design probes) and after (2026-10-22, this tree, ``archivist_probe.py`` beside the S7 gates)
   :header-rows: 1
   :widths: 30 30 14 26

   * - pair
     - why it is two spaces
     - before
     - after, and the gate
   * - 1. the measure spaces of two rules on :math:`S^2` with different
       weights
     - the weights are the metric
     - 21 of 105 pairs equal among 15 shipped rules, 20 of them with
       different weights
     - 0 of 105; ``test_axis_identity.py`` I1, I1b
   * - 2. ``level_symmetric(4)`` against ``product(2, 12)``: equal weights
       (all :math:`4\pi/24`), different directions
     - the atoms
     - equal axes
     - unequal; I2, and I2b for one rule about two poles
   * - 3. the real harmonics against their mirror-even restriction, at
       :math:`L = 1, 2, 3`
     - the domain
     - equal heads
     - unequal; ``test_heads_on_the_axis.py`` H1
   * - 4. indicator heads on edges :math:`[0,1,2]` and :math:`[0,0.5,2]`,
       raw and frame-dressed
     - the edges
     - equal, raw and dressed (dressed metrics :math:`[1, 1]` against
       :math:`[2, 0.667]`)
     - unequal; H2, and H2b for shifted edges of equal width
   * - 5. a dense-dressed head against the same axes with another Gram
     - the dense form
     - separated only by round-off in the unit's last bits, so two heads
       whose units came out bit-identical would compare equal whatever
       their metrics
     - unequal, by the overlay; I5, I5b, and I5c for one matrix with its
       inverse supplied against derived
   * - 6. ``from_blocks(a, b).dual()`` against ``from_blocks(a*, b*)``,
       hand-named ``a*``, ``b*``
     - covariance
     - equal, with equal hashes (the composite folded member NAMES)
     - unequal; I6

Two consequences of the key are worth stating. The measure space of a
frame over a quadrature IS the :math:`S_N` angular space over the same
rule (one :math:`L^2(\mu)`, not two): :meth:`FunctionSpace.of_measure
<orpheus.numerics.space.FunctionSpace.of_measure>` mints it through
:meth:`Axis.of_measure <orpheus.numerics.axis.Axis.of_measure>` under the
measure's :attr:`~orpheus.numerics.measure.DiscreteMeasure.axis_label`
(its phase, ``"angular"`` for a rule on the sphere or an orbit space), and
`[M]` ``measure_space(angular_frame(level_symmetric(6), 2)) ==
of_axes(Axis.of_quadrature(level_symmetric(6)))`` (gate LB-ii). Before S7
the measure space was a bare ``FunctionSpace("L2[…]",
inner_product_weights=w)`` and never equal to the angular space. And no
production number moved: `[M]` the S7 no-move baseline (100 cases:
solves, operators, heads, frames, the honesty registry) differs in 0 of
100, with its axis-weight control reddening 68 (S7b commit ``ac6c6558``).

.. _spaces-identity-bridge:

The derived name — once the identity bridge, now the readable label
-------------------------------------------------------------------

⛔ **The bridge reading was retired on 2026-09-07** by the identity flip
(structural ``__eq__``, CS4c step 6). The section keeps its label and its
history, because the bridge is why the chartered doctrine was already
true for axis-built spaces between CS1 (2026-08-20) and the flip, and
because the name derivation it describes still ships — what changed is
what that derivation is FOR.

:meth:`of_axes <orpheus.numerics.space.FunctionSpace.of_axes>` computes
the name deterministically and injectively from the axes' structural
content — a length-prefixed, type-tagged content digest, never Python's
``hash()``, so it is stable across processes. The digest is
:func:`~orpheus.values.content.name_digest` over the tuple of axes,
the first 8 bytes of the one content encoder's digest in hexadecimal
(:ref:`structured-geometry-content-identity`).

**Until the flip that injectivity WAS the identity.** Space identity was
``(name, shape)`` for every space, so the only way an axis-built space
could carry *metric differences imply space differences* was to mint a
different NAME for every different axis tuple: different axes → different
names → different spaces, with no flag day. That is the bridge.

**The identity reads the axes directly.** An axis-built space's identity
key is its axes and their overlay (:ref:`spaces-identity-key`), so the
doctrine holds without passing through a string. The derived name
survives as two things: the readable LABEL a human reads in a traceback,
and what keeps ``repr`` and the space-content guard messages
content-distinguishing. It no longer carries any composite's identity: a
full-field or coupled composite is keyed by its members' keys, because a
member's name does not carry covariance (#631 pair 6, finding F7 of the
S6a review; until S7 :meth:`FullFieldSpace.from_blocks
<orpheus.numerics.spaces.full_field_space.FullFieldSpace.from_blocks>`'s
digest of its members' names was that identity).

The mechanism is visible in one construction. The homogeneous problem's
quotient point and a genuine one-cell mesh with :math:`V_{\rm cell} = 2`
have the same shape :math:`(2, 1)` and the same readable prefix, and are
**unequal** spaces:

.. code-block:: text

   quotient point (volumes = [1.0])  ->  energy(2,)*spatial(1,)#<digest A>
   one-cell mesh  (volumes = [2.0])  ->  energy(2,)*spatial(1,)#<digest B>
   before the flip:  digest A != digest B  ->  UNEQUAL  (via the names)
   after  the flip:  axes differ in measure bytes  ->  UNEQUAL  (directly)

`[M]` 2026-09-07 on the landed carve, building both spaces as
``FunctionSpace.of_axes(EnergyAxis.synthetic(2), Axis("spatial", (1,),
weights=w, kind=BasisKind.NODAL))`` with ``w = [1.0]`` and ``w = [2.0]``
(the all-ones measure canonicalizes to ``None``, so the two axes differ
in their measure bytes):

.. list-table:: What the flip changed, on one pair of mints
   :header-rows: 1
   :widths: 46 18 18 18

   * - comparison
     - before
     - after
     - reads
   * - the two mints above (different measures)
     - ``False``
     - ``False``
     - unchanged — the doctrine held through the bridge, and now holds
       directly
   * - two independent mints of the SAME axes
     - ``True``
     - ``True``
     - unchanged; the hashes agree too
   * - an axis-built space vs a hand-named
       :class:`~orpheus.numerics.space.FunctionSpace` carrying its exact
       ``name`` and ``shape``
     - ``True``
     - ``False``
     - **flipped** — a label is no longer a way to spell somebody
       else's space
   * - ``A * B`` (both factors axis-built) vs
       ``FunctionSpace.of_axes(*A.axes, *B.axes)``
     - ``False``
     - ``True``
     - **flipped** — Q-T4 realized: an axis product is not a different
       *kind* of space, and the two mints derive DIFFERENT names
       (``…#a ⊗ …#b`` against one concatenated digest) while carrying
       the same axes

The last two rows are the whole content of the flip. The third says the
derived name stopped being a *credential* — before, anything that could
spell the string was that space. The fourth says the name stopped being
the *identity* — two spellings of one axis tuple were two spaces, and are
now one. `[R]` the "before" column is derivable rather than re-run: the
pre-flip body is exactly ``self.name == other.name and self.shape ==
other.shape`` (``git show 823f97dd:orpheus/numerics/space.py``), so it is
decided by the names printed above.

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

⛔ **This paragraph claimed** *"`[M]` the refusal has no production
witness yet, deliberately. The only axis mint in* ``orpheus/`` *today is*
``MaterialMesh.bulk_space``\ *, whose factors are both* ``NODAL``\ *…
The arm becomes production-reachable when CS2 mints the harmonic axis"*
**until 2026-09-07**, and both halves had gone false — not by item
6.2a, which touches no axis declaration, but by the CS4b axis migration
that preceded it and was never reconciled here.

`[M]` re-measured 2026-09-07, by AST over ``orpheus/**/*.py``:
:meth:`of_axes <orpheus.numerics.space.FunctionSpace.of_axes>` has
**seven** production call sites, not one. And the ``False`` row now has
a production occupant: :attr:`SNProblem.angular_trial_space
<orpheus.sn.problem.SNProblem.angular_trial_space>` appends the
scheme-owned ``moment_axis`` — declared modal (``Modal(e_0)``), because
a DG cell moment is an expansion coefficient — to an axis-built base
whenever the bound scheme is multi-moment, and `[M]`
``LinearDiscontinuous.is_multi_moment`` is ``True`` while
``DiamondDifference``'s is ``False``. So an LD carrier's trial space is
an axis-built space with a MODAL factor. ⚠ Scope: that measures the
*occupant*, not the *firing* — whether any consumer asks
``has_coordinate_cone`` of that space is a separate census this page
does not carry. The named test witness stands either way
(``tests/gates/numerics/test_field.py``, gates E1/E2 — the refusal and the
same values answered on an all-nodal space, the positive-and-negative
pair ``vv-principles`` anti-pattern #11 requires).

What survived intact is the CLOSURE ARGUMENT and its consequence for
the harmonic family. Since ``of_axes`` is the only ROOT producer of an
``axes`` record (``*`` and :meth:`dual
<orpheus.numerics.space.FunctionSpace.dual>` merely thread one through,
so both need an axis-built ancestor), and the angular head declares no
axes, every harmonic-moment space in the tree was still legacy on that
date (``axes is None``) and therefore took the ``None`` arm — `[M]`
2026-09-07 the moment product's metric is a
:class:`~orpheus.numerics.metric.FactoredMetric` while its ``axes`` is
``None``. Item **6.2a** did not move that: retiring the densifier
changed how a product carries its measure, not whether its factors
declare axes. Item **6.2c** does, by minting the harmonic axis (this
sentence said "CS2" until the 2026-09-07 ruling).


⛔ **The last reading above is dated 2026-09-07; since S7 every basis's
coefficient space is axis-built.** `[M]` 2026-10-22
(``archivist_probe.py``, :meth:`FunctionSpace.of_basis
<orpheus.numerics.space.FunctionSpace.of_basis>` of each basis):
``SphericalHarmonicBasis(L=1)``, its mirror-even restriction and
``LegendreBasis(L=2, axis="x")`` each span one ``Modal`` head axis
(``has_coordinate_cone`` ``False``), and an ``IndicatorBasis`` on edges
:math:`[0, 1, 3]` spans one ``NODAL`` cell axis with no measure
(``has_coordinate_cone`` ``True``); the indicator head a frame dresses is
``Modal(masses)`` (unit :math:`[1, 2]`, cone ``False``; ruling 83). The
production moment product reads ``False`` too: `[M]` 2026-09-22
(``scratch/_definitions/round2/moment_space_probe.py``),
``SNProblem.moment_space(1)`` is an axis-built ``TensorProductSpace``
with axes (``MODAL``, ``NODAL``, ``NODAL``) and cone ``False`` on a slab,
a sphere and a 2-D Cartesian problem, where each problem's all-nodal
``bulk_space`` reads ``True`` as the control. The ``None`` arm now holds
only the hand-built, axes-less spaces.

.. _spaces-basis:

The basis: the synthesis side of a frame, and what its coefficients are
------------------------------------------------------------------------

A **basis** is a finite family of functions
:math:`\{\phi_k\}_{k=1}^{K}` on a manifold :math:`M`, carried together
with the operations that turn a coefficient vector into function values.
The manifold is what the functions **eat**
(:attr:`Basis.domain <orpheus.numerics.basis.base.Basis.domain>`: a
direction on :math:`S^2`, an orbit of a stabiliser, a point of a
partitioned spatial or energy domain, an index); the coefficients live
in :math:`\mathbb R^K`, the basis's **coefficient space**
(:meth:`FunctionSpace.of_basis(basis)
<orpheus.numerics.space.FunctionSpace.of_basis>`). Those
are two different levels of the three-level picture, the manifold and
the coefficients, with the fields on the manifold between them
(:ref:`manifold-three-levels`). This section is the corpus's one
definition of the word; the abstract class is
:class:`~orpheus.numerics.basis.base.Basis`.

**What a basis computes.** Its central operation is the **naked
synthesis**

.. math::

   S_0 : \mathbb R^{K} \to \mathbb R^{N},
   \qquad
   (S_0\,c)_n \;=\; \sum_{k=1}^{K} \phi_k(x_n)\, c_k ,

the expansion evaluated at sample points :math:`x_n`, with no weights
and no dual factor; frame theory calls it the synthesis operator,
:math:`T^{*}` in the convention ORPHEUS follows, where :math:`T` is the
analysis operator (:cite:`CasazzaLynch2016` §4.1; Christensen 2016,
*An Introduction to Frames and Riesz Bases*, 2nd ed., calls the same
operator :math:`T`; :ref:`frame-discrete-frame-definition`). A basis
tabulates :math:`\Phi_{nk} = \phi_k(x_n)` once
(:meth:`~orpheus.numerics.basis.base.Basis.evaluate`, the only method that
takes points) and every other operation contracts that cached table, so
the frame that holds the table never re-tabulates on the hot path. The
four weighted contractions are each :math:`S_0` or its transpose times
**one** diagonal family:

.. list-table:: The basis's contractions (``basis/base.py``'s module docstring is the source)
   :header-rows: 1
   :widths: 24 30 46

   * - Verb
     - Action
     - Which weight, and why
   * - ``analyze``
     - :math:`M f = S_0^{\mathsf T}(w \odot f)`
     - the **measure's** weights :math:`w_n`: analysis is the measured,
       test side
   * - ``analyze_transpose``
     - :math:`M^{\mathsf T} c = w \odot S_0 c`
     - the same :math:`w_n`; the representation transpose, NOT the
       Hilbert adjoint
   * - ``reconstruct``
     - :math:`R c = S_0(d \odot c)`
     - the basis's own **canonical-dual factor** :math:`d_k`
       (:math:`2\ell + 1` for the real spherical harmonics); no measure
   * - ``reconstruct_transpose``
     - :math:`R^{\mathsf T} f = d \odot S_0^{\mathsf T} f`
     - the same :math:`d_k`; measure-free, like ``reconstruct``

The table is the conceptual unity, and the code keeps each row as its
own fused contraction rather than ``weight ⊙ synthesize``: floating-point
non-associativity would move the factored form at the ULP level, and the
scattering kernel that consumes these rows is pinned at 0 ULP. The
discrete Gram
:math:`\sum_n w_n\,\phi_j(x_n)\,\phi_k(x_n)`
(:attr:`~orpheus.numerics.frame.FrameBase.discrete_gram`) is the matrix
of the **Gram operator** :math:`T T^{*}` on the :math:`K` coefficients,
read with the Euclidean coefficient metric in which :math:`T^{*}` is the
naked synthesis :math:`S_0`. It is not the frame operator
:math:`S = T^{*} T`, which acts on the :math:`N` sample values; the two
share their nonzero spectrum, so frame bounds read off either agree, but
they are different objects (:ref:`frame-discrete-frame-definition`). It
equals the continuum Gram when the
quadrature is exact to the basis's degree. It is the frame's
:eq:`frame-discrete-gram`, computed from the cached table; the basis
alone does not mint it, because it pairs the basis with a measure and
that pairing is the frame (the former ``Basis.mass_matrix`` was retired
onto it).

**The basis is choice-free; the frame makes the choice.** A basis knows
its functions and their normalization convention and nothing about which
points sample them. The :class:`~orpheus.numerics.measure.DiscreteMeasure`
supplies the nodes, the weights and therefore the domain's inner product;
a :class:`~orpheus.numerics.frame.FrameBase` binds a basis to a measure
and mints the two faces from the pair (:eq:`galerkin-pair`). The two
sides split as follows:

- the **basis is the synthesis (trial) side**: it owns the reconstruction
  :math:`R`, and the frame's codomain is its coefficient space;
- the **analysis side belongs to the frame**: the analysis face
  :math:`M` is minted from a test basis and the measure's weights. A
  :class:`~orpheus.numerics.frame.GalerkinFrame` binds the test basis to
  the trial basis (``test is trial``); a
  :class:`~orpheus.numerics.frame.PetrovGalerkinFrame` holds a separate
  test basis, and the shipped one,
  :class:`~orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis`
  (cell indicators times a nodal weight such as the region flux), builds
  only ``evaluate`` and ``analyze`` and refuses the other four verbs
  (``synthesize``, ``reconstruct``, ``reconstruct_transpose``,
  ``analyze_transpose``) by name: none has a consumer,
  and a Petrov–Galerkin reconstruction is trial-side.

The separation is why one basis can serve several frames: the same
basis bound to two different quadratures is two frames with two
different discrete Grams, because the Gram belongs to the pairing and
not to either half. It is also why a pairing can be wrong while each
half is right, and why the basis must answer what it eats: ERR-080 was
the full-sphere real harmonics evaluated at the ordinates of a 1-D rule,
whose :math:`\mu_y = \mu_z = 0` meant "no azimuthal information" and
were read as "azimuth zero", a pairing the level-2 spaces accepted and
only the level-1 manifolds refute (:ref:`manifold-three-levels`). The symmetry a
basis HAS is read off its domain, never stored
(:attr:`Basis.invariance_group
<orpheus.numerics.basis.base.Basis.invariance_group>`, ``@final``): a
function on a quotient :math:`M/H` is exactly an :math:`H`-invariant
function on :math:`M`.

**What the coefficient space carries.** A basis carries two values its
coefficient space is read from, with no space import on the basis side
(the concept order puts ``basis`` below ``axis``, so a basis cannot mint
its own axis; ruling 81): its **layout**
(:attr:`Basis.layout <orpheus.numerics.basis.base.Basis.layout>`, a value
from :mod:`orpheus.numerics.coefficient_layout`) and its **continuum Gram**
diagonal over that layout (:attr:`Basis.coefficient_weights
<orpheus.numerics.basis.base.Basis.coefficient_weights>`; for the real
harmonics :math:`g_C = \operatorname{diag}\bigl(4\pi/(2\ell+1)\bigr)`
spread over each degree's live columns, zero on the padding).
:meth:`FunctionSpace.of_basis <orpheus.numerics.space.FunctionSpace.of_basis>`
is the one coefficient-space mint: one generic axis over
:meth:`Axis.of_basis <orpheus.numerics.axis.Axis.of_basis>`, labelled and
shaped by the layout, measured by ``coefficient_weights``, on the basis's
domain, with origin ``FromBasis(basis)``. The frame re-dresses that space
with the **discrete** Parseval metric (:ref:`frame-parseval-metric`), and
the dressed space, not the basis's, is the moment space the tree binds
(:ref:`frame-moment-space-single-home`). An indicator basis is
measure-free: its ``coefficient_weights`` is ``None``, because the cell
masses exist only against a bound measure (a spherical shell is
:math:`4\pi r^2\Delta r`, a coarse group counts its fine groups), so its
raw head's metric verbs stay Euclidean and its axis marginal refuses to
integrate, naming the frame (ruling 88, :ref:`spaces-collapse-pair-frame`).

**The Gram structure is a declaration the type enforces.** A trial
basis's frame normalises the projection by the cross Gram :math:`MR`,
and the frame computes that Gram with one row-sum probe,
``analysis(reconstruction(ones))``. The probe equals the required
normalisation under exactly two structural conditions, so each basis
declares which one holds
(:attr:`Basis.gram_structure <orpheus.numerics.basis.base.Basis.gram_structure>`,
a :class:`~orpheus.numerics.basis.base.GramStructure`):

- ``DIAGONAL`` — disjoint supports or orthogonal functions: :math:`MR` is
  diagonal and the row sum is the diagonal;
- ``PARTITION_OF_UNITY`` — overlapping functions whose rows sum to one:
  :math:`MR` is not diagonal, but :math:`R\mathbf 1 = \mathbf 1`
  collapses the probe to the per-region weight anyway;
- ``DENSE`` — neither: the true projection needs the dense
  :math:`(MR)^{-1}M` solve, which is not built, so
  :func:`project <orpheus.numerics.frame_faces.project>`
  REFUSES a dense trial rather than return a silently wrong coarsening.

The base class answers ``DENSE``, the safe refusal: a new basis earns the
shortcut only by declaring, having checked its Gram, that it may use it.

**The family, measured.** `[M]` 2026-09-22, every concrete ``Basis``
subclass in ``orpheus/`` (an AST census of class definitions whose
bases name a basis class, seven production members, which a runtime
``__subclasses__`` walk reproduces once
``orpheus.sn.operators.loss_kernel_gauge`` is imported; the three test
doubles the same census finds under ``tests/`` are excluded, one of them
nested inside a test function, where a header grep anchored at column 0
misses it). The coefficient-space column is `[M]` 2026-10-22
(``archivist_probe.py``) for the four members constructed there, and read
from the code (each class's ``layout`` and ``coefficient_weights``) for
the other three:

.. list-table:: The concrete bases
   :header-rows: 1
   :widths: 30 24 16 30

   * - Class
     - Domain (what it eats)
     - Gram structure
     - Coefficient space
   * - :class:`~orpheus.numerics.basis.indicator_basis.IndicatorBasis`
     - the partitioned manifold (spatial, energy, index)
     - ``DIAGONAL``
     - one ``NODAL`` cell axis, :class:`~orpheus.numerics.coefficient_layout.CellLayout`
       of the edges, no measure (cone ``True``)
   * - :class:`~orpheus.numerics.basis.overlap_basis.OverlapBasis`
       (an ``IndicatorBasis``)
     - the fine partition
     - ``PARTITION_OF_UNITY``
     - as ``IndicatorBasis``
   * - :class:`~orpheus.numerics.basis.weighted_indicator_basis.WeightedIndicatorBasis`
     - the wrapped indicator's
     - ``DENSE`` (the default)
     - the wrapped indicator's layout; a TEST basis only
   * - :class:`~orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis`
     - :math:`S^2`
     - ``DIAGONAL``
     - one ``Modal`` axis, :class:`~orpheus.numerics.coefficient_layout.HarmonicLayout`,
       measure :math:`g_C`, support :math:`S^2` (cone ``False``)
   * - :class:`~orpheus.numerics.basis.spherical_harmonic_basis.MirrorEvenSphericalHarmonicBasis`
     - :math:`S^2/\sigma_a`, a quotient
     - ``DIAGONAL``
     - as its parent, on support :math:`S^2/\sigma_y`: a different space
   * - :class:`~orpheus.numerics.basis.legendre_basis.LegendreBasis`
     - :math:`S^2/O(2)_a`, a quotient
     - ``DIAGONAL``
     - one ``Modal`` axis, :class:`~orpheus.numerics.coefficient_layout.LegendreLayout`,
       support :math:`S^2/O(2)_a` (cone ``False``)
   * - :class:`~orpheus.sn.operators.loss_kernel_gauge.LossKernelBasis`
     - an index set of trace degrees of freedom
     - ``DIAGONAL``, a theorem: the table is :math:`G`-orthonormal
     - one ``Modal`` axis, :class:`~orpheus.numerics.coefficient_layout.ModeLayout`
       with a declared zero unit, no measure, support the block's index
       set (rulings 86, 87)

**Nodal and modal belong to the axis, not to the basis.**
:class:`~orpheus.numerics.axis.BasisKind` is an axis slot
(:ref:`spaces-the-axis`), and it records what one coefficient of the
factor IS: ``NODAL`` when it is a point or cell **value** (an
indicator-like basis, so a coordinate cone exists and a sign test on the
coefficients is a sign test on the function), ``MODAL`` when it is an
**expansion coefficient** (a positive function can have negative ones,
:ref:`spaces-nodal-modal`). A basis enters that slot through its layout:
``FromBasis(basis).kind`` is one exhaustive ``match`` on
``basis.layout``, ``NODAL`` for a cell layout (unit ones) and
``Modal(layout.unit)`` for a degree or mode layout. So every one of the
7 members has a kind on its coefficient axis: the three indicator
members ``NODAL``, the three harmonic members and the loss-kernel basis
modal. The kind of a frame-dressed head is the frame's to say
(``FromFrame``): a dressed head holds moments, so it is
``Modal`` with unit :math:`\Phi^{\mathsf T} w`, and the dressed
indicator head is ``Modal(masses)`` until S6c (ruling 83, #633).

**Two related surfaces.** The harmonic family also satisfies the
structural protocol :class:`~orpheus.numerics.basis.base.TruncatedBasis`
(a truncation order :math:`L`, a degree layout, a domain, the continuum
Gram ``coefficient_weights``, and ``at_order``), which is
what a consumer asks for when it needs "the coefficient space of order
:math:`L`" without naming one class; and
:class:`~orpheus.numerics.basis.descent.Descent` pulls a basis back along
a quotient map (:ref:`manifold-descent`).

.. _spaces-axis-generator:

.. _spaces-axis-origin:

The origin: what minted an axis, and what its content keeps
-----------------------------------------------------------

Every axis records the object that minted it, in one slot,
:attr:`Axis.origin <orpheus.numerics.axis.Axis.origin>`. The slot's type is
a closed sum of five arms (ruling 64 of the numerics untangling,
:mod:`orpheus.numerics.axis`):

.. math::

   \texttt{Origin} \;=\; \texttt{FromMeasure}(\mu)
   \;\mid\; \texttt{FromQuadrature}(q)
   \;\mid\; \texttt{FromBasis}(b)
   \;\mid\; \texttt{FromFrame}(f)
   \;\mid\; \texttt{Declared}(\kappa) .

The axis does not take its kind, its support, its atoms or its coefficient
layout as constructor arguments. It reads all four off its origin, so the
generating object and the axis cannot disagree about them:

.. list-table:: The five arms (the module docstring of :mod:`orpheus.numerics.axis` is the source)
   :header-rows: 1
   :widths: 22 20 20 16 22

   * - arm, and the object it carries
     - kind it implies
     - support (content)
     - atoms (content)
     - coefficient layout (content)
   * - :class:`~orpheus.numerics.axis.FromMeasure` — a discrete measure :math:`\mu`
     - ``NODAL``
     - :math:`\mu`'s support
     - :math:`\mu`'s nodes
     - none
   * - :class:`~orpheus.numerics.axis.FromQuadrature` — a directional quadrature :math:`q`
     - ``NODAL``
     - the support of :math:`q`'s measure
     - the nodes of :math:`q`'s measure
     - none
   * - :class:`~orpheus.numerics.axis.FromBasis` — a basis :math:`b`
     - the layout's: ``NODAL`` for a cell layout, :class:`Modal(unit)
       <orpheus.numerics.axis.Modal>` for a degree or mode layout
     - :math:`b`'s domain
     - none
     - :math:`b`'s layout
   * - :class:`~orpheus.numerics.axis.FromFrame` — a frame :math:`f`
     - ``Modal``, unit :math:`\Phi^{\mathsf T} w` (a mode basis's
       declared unit instead)
     - the domain of :math:`f`'s trial basis
     - none
     - the trial basis's layout
   * - :class:`~orpheus.numerics.axis.Declared` — a kind :math:`\kappa`
     - :math:`\kappa`
     - none
     - none
     - none

``Axis(kind=NODAL, generator=basis)``, a nodal axis claiming a spectral
generator, was spellable before S7 and is not now: there is no ``kind``
parameter, and ``FromBasis(b).kind`` is one exhaustive ``match`` on
``b.layout`` that raises on a value that is not a coefficient layout.
``dataclasses.replace(axis, origin=...)`` re-runs ``__post_init__``, which
re-derives the kind, the support, the atoms and the layout, so a
re-dressed axis cannot keep a stale kind.

**Which arm each shipped axis takes.** Each row is pinned by a named gate,
cited rather than copied:

- the angular axis of every :math:`S_N` problem is ``FromQuadrature``
  (:meth:`Axis.of_quadrature <orpheus.numerics.axis.Axis.of_quadrature>`;
  ``tests/gates/numerics/test_axis_generator.py::TestG3MintFidelity``);
- the rank-1 spatial axis of a carrier is ``FromMeasure`` over the
  carrier's own volume measure (cell centres as atoms, volumes as weights;
  ``tests/gates/numerics/test_space_of_axes.py``
  ``::test_the_spatial_axis_is_minted_through_the_carriers_own_measure``,
  whose centres and volumes are derived by hand from the edge list);
- a basis's continuum coefficient axis is ``FromBasis``
  (:meth:`Axis.of_basis <orpheus.numerics.axis.Axis.of_basis>`), and the
  moment head a frame dresses is ``FromFrame``
  (:func:`~orpheus.numerics.frame_faces.basis_space`;
  ``tests/gates/numerics/test_heads_on_the_axis.py`` H0c, H2c);
- four axes are ``Declared`` because no generating object exists, or none
  may be imported into the numerics layer (ruling 82): every
  :class:`~orpheus.numerics.axis.EnergyAxis` (its content is its edges), the
  homogeneous pose's one-point counting axis, the rank-:math:`d` spatial
  axis of a :math:`d \ge 2` carrier, and the linear-discontinuous scheme's
  spatial-moment axis, which is ``Modal(e_0)`` and whose generator, the
  scheme, lives in the transport layer
  (``tests/gates/numerics/test_axis_generator.py``
  ``::TestG5TheQuadratureRefusalIsOneSpelling::test_G5_6a_the_shipped_atomless_axes_are_Declared``,
  ``tests/gates/numerics/test_space_of_axes.py``
  ``::test_the_rank_d_spatial_axis_is_generator_less_BY_CONTRACT``).

The rank-:math:`d` arm is ``Declared`` by contract (G6b): a
:class:`~orpheus.numerics.measure.DiscreteMeasure` is a flat atom list,
nodes :math:`(N, d)` and weights :math:`(N,)`, so
:meth:`Axis.of_measure <orpheus.numerics.axis.Axis.of_measure>` mints at
rank 1, and a rank-:math:`d` spatial axis has shape
:math:`(n_x, n_y, n_z)`. Minting it flat would be a different space for
every :math:`d \ge 2` carrier. The day a rank-:math:`d` measure-to-axis
pairing exists, the contract row is inverted deliberately. A consequence
the S7 review recorded (ruling 85): a 1-D slab shifted rigidly is a
different space, because its atoms move, while the same shift of a 2-D
mesh leaves the space equal, because the rank-:math:`d` arm carries no
atoms.

**The forgetful map, and what it forgets.** Write :math:`g` for a
measure-presenting origin (a measure, or the measure of a quadrature),
:math:`\mu_g = \sum_{i=1}^{n} w_i\,\delta_{x_i}` on a manifold
:math:`M_g`. The axis it mints keeps

.. math::
   :label: spaces-axis-forgetful-map

   \mu_g \;=\; \sum_{i=1}^{n} w_i\,\delta_{x_i},
   \qquad
   \mathcal{F}_\ell(g) \;=\;
   \bigl(\,\ell,\ (n,),\ w,\ \texttt{NODAL},\ M_g,\ x\,\bigr) ,

.. (vv-status rationale) Structural/representational identity: it STATES
   what ``Axis.of_measure`` / ``Axis.of_quadrature`` construct from a
   measure — the label, the rank-1 index shape ``(n_points,)``, the weight
   vector, the NODAL kind implied by the origin's arm, the support and the
   atoms (ruling 71, S7). Not a solver claim — no flux, no eigenvalue, no
   discretization error. The verifiable content is
   ``tests/gates/numerics/test_axis_generator.py`` G3 (the mint-fidelity
   roster over all five shipped ``Quadrature`` factories) and the S7 rows
   ``tests/gates/numerics/test_axis_identity.py`` I1, I1b, I2, I2b (the
   weights, the atoms and the support each separate two axes).
.. vv-status: spaces-axis-forgetful-map documented

where :math:`\ell` is the factor's label. Until S7 the content stopped at
the kind: the axis kept the weights and dropped the atoms and the
manifold. Since S7 (ruling 71) it keeps both, because they fix the factor's
space (:ref:`spaces-axis-identity`). What :math:`\mathcal{F}_\ell` still
forgets is only what is not a property of the space: the side channels the
generating object carries beside its measure. For a quadrature that is the
polar-level fibration (``level_indices``, carried by its
:class:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure`) and the
folding record; for a frame it is the frame itself, the one object that can
re-dress a moment head at another order. ``NODAL`` is not a parameter: a
discrete measure's components are point values with a coordinate cone
(:ref:`spaces-nodal-modal`), so the arm implies it.

**The section law.** The mint records its own input in the origin, and the
projection :math:`p(a)`, the object inside ``a.origin``, reads it back:

.. math::
   :label: spaces-axis-generator-section

   p\bigl(m_\ell(g)\bigr) \;=\; g \quad\text{for every generator } g,
   \qquad
   m_{\ell(a)}\bigl(p(a)\bigr) \;=\; a \quad\text{for every minted axis } a ,

.. (vv-status rationale) Structural identity (the SECTION law): it STATES
   that the origin-minted axis and the origin projection are a two-sided
   inverse pair — the mint records its own input (``a.origin.quadrature is
   g`` / ``a.origin.measure is g``, identity checks) and re-minting from the
   recorded object reproduces the axis (``Axis.of_quadrature(g, a.label) ==
   a``, a structural-equality check). Not a solver claim. The verifiable
   content is ``tests/gates/numerics/test_axis_generator.py`` G8
   (``TestG8TheMintIsASectionOfTheForgetfulMap``, the angular roster plus
   the measure leg) and the spatial leg in
   ``tests/gates/numerics/test_space_of_axes.py``.
.. vv-status: spaces-axis-generator-section documented

with :math:`m_\ell` one of :meth:`Axis.of_measure
<orpheus.numerics.axis.Axis.of_measure>` or :meth:`Axis.of_quadrature
<orpheus.numerics.axis.Axis.of_quadrature>`. The law's domain is the
measure- and quadrature-minted axes. A ``FromBasis`` axis satisfies the
same law through :meth:`Axis.of_basis <orpheus.numerics.axis.Axis.of_basis>`
(:ref:`spaces-moment-head-axis-built`); a ``Declared`` axis is outside the
domain, not a counterexample, because it has no generating object. A
consumer holding only the space therefore recovers ``mu_x``, ``eta``,
``mu_z`` and ``level_indices`` without being handed the quadrature
separately, and a consumer that needs the frame of a dressed head gets it
from ``FromFrame``.

.. _spaces-generator-why-quadrature:

Why the angular arm is the QUADRATURE and not its measure
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``FromQuadrature(q)`` and ``FromMeasure(q.measure)`` mint one axis: the
content is the measure's, and the arm is provenance. The angular arm exists
because the measure does not answer everything the angular consumers read
through the space. They read four names, and a bare measure answers three:

.. list-table:: What each candidate generator can answer
   :header-rows: 1
   :widths: 22 30 48

   * - Name
     - On a ``DiscreteMeasure``?
     - Where it actually lives
   * - ``mu_x`` / ``mu_z``
     - the data is there
     - columns of ``measure.nodes`` (shape :math:`(N, d)`), read by
       :meth:`Quadrature.axis_cosines
       <orpheus.numerics.quadrature.directional.Quadrature.axis_cosines>`.
       `[M]` 2026-08-29, ``q.mu_x`` is bit-equal to ``measure.nodes[:, 0]``
       and ``q.mu_z`` to ``measure.nodes[:, 2]``.
   * - ``eta``
     - the data is there
     - the same column 0 under its cylindrical-frame name
       (:math:`\eta = \hat\Omega\cdot\hat r`), an alias, not a second
       array.
   * - ``level_indices``
     - **no**
     - the :class:`~orpheus.numerics.quadrature.rules_sphere.LevelStructure`
       side channel, carried by the ``Quadrature``. A ``DiscreteMeasure``'s
       fields are ``nodes``, ``weights``, ``support``,
       ``invariance_group`` and ``exactness``.

The level fibration is structure the atom list does not contain: the
polar-level partition of the ordinate set is a datum the curvilinear
:math:`\alpha`-recursion needs. Recording the measure would have shipped a
provenance that answers three of the four names. The general lesson is to
choose a provenance type by what the consumer reads, not by which object
ranks lower: when the lower object is a projection of the higher one, the
layer argument narrows the contract silently. The refutation leg is gated:
``tests/gates/numerics/test_axis_generator.py``
``::TestG4TheFourNamesAnswerThroughTheSpace::test_the_bare_measure_cannot_answer_the_fourth_name``.

.. _spaces-generator-route-gate:

Reading the origin: one ``match``, one refusal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A consumer reads the origin by an exhaustive ``match``, never through a
narrowing accessor. `[M]` 2026-10-22, ``git grep "case FromQuadrature"``
over ``orpheus/``: five production match sites, in four functions, recover
the angular geometry through the space: the two angular closures in
:mod:`orpheus.sn.angular.closure`
(``MorelMontryAngularSweep``, ``IdentityAngularClosure``), the streaming
terms of :class:`~orpheus.sn.mesh.reduced_operator.ReducedStreamingOperator`,
and :meth:`HarmonicFrame.for_space
<orpheus.transport.frames.harmonic_frame.HarmonicFrame.for_space>` and
:meth:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame.moment_space_on`
(the second a presence guard). Each matches ``case FromQuadrature(quad)``
and, on any other arm, raises
:func:`~orpheus.numerics.axis.quadrature_refusal`, the one spelling of the
message. It names the axis, the asking consumer and the arm the axis was
minted through, so a refusal is diagnosable from its text alone. The match
is a type claim, not a ``None`` check: a measure-minted axis handed to a
consumer that needs a quadrature refuses, naming ``FromMeasure``. Two more
``match`` sites read the other arms:
:func:`~orpheus.numerics.spaces.moment_head.truncated_head` re-mints a head
through ``FromFrame`` or ``FromBasis`` (:ref:`spaces-moment-head-axis-built`),
and the axis marginal refuses a raw cell head through ``FromBasis``
(:ref:`spaces-collapse-pair-frame`). The census of the ``FromQuadrature`` sites is
gated: ``tests/gates/numerics/test_axis_identity.py`` O3.

**Proving that a consumer reads through the axis.** When a re-point changes
where a consumer gets a datum and leaves the datum alone, every
before-and-after value gate compares :math:`X` with :math:`X`: the bound
quadrature IS the object the consumer used to receive. The instrument that
discriminates is a **decoy**: an axis minted from a quadrature that carries
different angular data, so that any movement of the produced value can only
come from a read that went through the axis to its origin. Two decoys serve
two reads, and since S7 they differ in what space identity sees:

- the **node-scale** decoy (nodes :math:`\times\,0.9`, weights kept) moves
  ``mu_x``. Since S7 its axis is UNEQUAL to the true one, because the atoms
  are content (#631 pair 2); until S7 it was identity-equal, and that
  invisibility was a leg of the gate. The gate's route legs stand without
  it. `[M]` 2026-08-29: it moves 4 of 4 packets on a 4-cell slab, 4 of 4 on
  the sphere and 8 of 12 on the cylinder; the four cylinder packets that do
  not move are the :math:`\eta = 0` member of each :math:`\mu`-level, and
  :math:`0.9 \times 0 = 0`. It is the only node decoy both admission tiers
  accept on a curvilinear chart: a node roll breaks the dome's antisymmetry
  :math:`\sum_n w_n\mu_n = 0` and is refused by
  :func:`~orpheus.sn.angular.redistribution.angular_redistribution`;
  negation and reversal preserve the antisymmetry and are refused one tier
  later by the Morel–Montry closure's membership guard
  (:math:`\tau_m \in [0,1]`,
  :func:`~orpheus.sn.angular.closure.morel_montry_tau_per_level`).
- the **level-roll** decoy (the same measure, the level list rolled by one)
  leaves ``mu_x`` untouched and moves only the ``level_indices`` read. It is
  still invisible to space identity, because the level structure is a side
  channel and not content. `[M]` 2026-08-29: 4 of 12 packets move on
  ``folded_product(4, 6)``; the floor is 4 because the per-level radial
  cosines are a palindrome (:math:`0.440, 0.814, 0.814, 0.440`), so a roll
  of one fixes half the levels.

Gates: ``tests/gates/sn/mesh/test_reduced_operator.py``
``::…::test_K1_the_packet_reads_THROUGH_the_axis`` (per chart) and
``::…::test_K2_the_cylinder_index_read_is_a_separate_route``;
``tests/gates/sn/sweep/curvilinear/test_angular_closure.py``
``::TestP4RemTheClosureMintReadsThroughTheAxis``, where a weight decoy
witnesses the :math:`\Delta A/w` read. A weight decoy's admissibility
depends on the order: `[M]` 2026-08-29 on the sphere, weights
:math:`\times\,0.9` is admitted at ``gauss_legendre(2)`` and refused at
orders 4, 6 and 8 (:math:`\tau = 1.195, 1.047, 1.059`), because shrinking
every weight pushes the outermost :math:`\tau` past 1.

.. _spaces-generator-protocol:

The angular contract declares what its consumer reads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``AngularMeasure`` :class:`~typing.Protocol` in
:mod:`orpheus.sn.angular.redistribution`, the structural contract the
:math:`\alpha`-dome recursion consumes, declares ``level_structure:
LevelStructure | None`` and ``axis(label="angular") -> Axis``, with
:meth:`Axis.of_quadrature <orpheus.numerics.axis.Axis.of_quadrature>` as
its concrete implementer. ``None`` is the honest reading for a slab or 2-D
rule, which still answers ``level_indices`` with one degenerate level. The
member is declared and read directly because a defaulted
``getattr(x, "level_structure", None)`` inside a guard does not raise the
day the attribute is renamed: it takes the ``None`` branch and the
admission check disappears silently (``vv-principles`` #28). The label
defaults at the quadrature because a directional quadrature IS the angular
generator, so the literal ``"angular"`` appears at no mint site and a
label twin across mint sites cannot be spelled. Gate:
``tests/gates/sn/angular/test_redistribution.py``
``::TestG9TheProtocolDeclaresWhatItsConsumersRead``.

.. _spaces-generator-identity-exclusion:

.. _spaces-generator-identity-third-answer:

The origin is provenance; its content is identity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``origin`` is declared ``field(compare=False)`` on
:class:`~orpheus.numerics.axis.Axis`, the content encoder's one spelling of
"not content" (:ref:`structured-geometry-content-identity`), so equality,
hash and the ``of_axes`` name digest ignore it. What the origin implies
(the kind, the support, the atoms, the layout) is stored on the axis as
fields and IS content. The doctrine in one sentence: **two axes with the
same content are the same axis whatever object minted them**; the origin
records where the axis came from, the content records what it is.

Three facts keep the exclusion. First, it is the doctrine: an axis
re-weighted through a quadrature and the same content minted through that
quadrature's measure are one axis
(``tests/gates/numerics/test_axis_generator.py`` G1a, G1b), and the name
digest is blind to the arm (``tests/gates/numerics/test_space_of_axes.py``
``::test_of_axes_name_is_BLIND_to_the_generator``). Second, two carriers
built from equal inputs hold distinct rule instances, and their spaces are
equal because the instance is not content (issue #403's hazard never
reaches axis identity). Third, the arms that hold an object have no content
encoding of their own: a :class:`~orpheus.numerics.quadrature.directional.Quadrature`
is a mutable dataclass and a
:class:`~orpheus.numerics.measure.DiscreteMeasure` holds writeable arrays,
so admitting either as a part would make equality fall back to identity,
make ``hash`` raise ``TypeError`` and make the name mint raise
:class:`~orpheus.values.content.ContentlessError` (`[M]` 2026-10-02). The
ruling does not rest on the third fact: an origin with a content encoding
would still be excluded.

.. warning::

   **The origin's object is a LIVE REFERENCE; the axis's own copies are the
   authority.** A ``Quadrature`` is mutable, so the arrays of
   ``axis.origin.quadrature`` can in principle change after the axis was
   minted. The axis's ``weights`` and ``nodes`` cannot: ``__post_init__``
   stores ``w + 0.0`` and the atoms the same way, read-only and defensively
   copied. The factor measure and the atoms of record are ``axis.weights``
   and ``axis.nodes``; read the origin for what the axis does not keep (the
   level fibration, the frame).

.. _spaces-generator-not-a-reverse-accessor:

What the origin is NOT: an axis-to-measure accessor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CS4b phase S6.0b refused giving :class:`~orpheus.numerics.axis.Axis` a public
accessor that produces a :class:`~orpheus.numerics.measure.DiscreteMeasure`
(:ref:`spaces-collapse-pair-refuted`), and the refusal stands. Such an
accessor points axis to measure and must manufacture its output for every
axis that was not minted from one (a ``Declared`` axis, a coefficient head).
The origin points generator to axis and manufactures nothing: it records the
real object at the one moment it is in scope, the mint. The axis marginal
builds no measure either: it reads the unit :math:`u` and the metric block
:math:`G_k` the axis and its space carry (:ref:`spaces-axis-unit`).

.. _spaces-generator-none-inventory:

.. _spaces-generator-seams:

.. _spaces-generator-gates:

Verification — the gates, by name
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Cite these rather than copying their numbers; each re-measures itself.

.. list-table::
   :header-rows: 1
   :widths: 24 30 46

   * - Claim
     - Gate
     - What a mutation does to it
   * - the origin is provenance, the atoms are content
     - ``tests/gates/numerics/test_axis_generator.py``
       ``::TestG1GeneratorIsProvenanceNotIdentity``
     - ``G1a``: two content-equal rules mint one axis. ``G1b`` (re-posed
       at S7): ``FromQuadrature(q)`` and ``FromMeasure(q.measure)`` mint
       one axis, and the same weights ``Declared`` with no atoms are
       another axis. ``G1c``: a ``DiscreteMeasure`` still raises on ``==``.
   * - the digest is blind to the arm
     - ``tests/gates/numerics/test_space_of_axes.py``
       ``::test_of_axes_name_is_BLIND_to_the_generator``
     - admitting the origin into the content moves the derived name.
   * - the mint reproduces the literal
     - ``tests/gates/numerics/test_axis_generator.py``
       ``::TestG3MintFidelity``
     - pins the threading and the canonicalization; the two sides share
       the weight values (both read ``measure.weights``).
   * - the four names answer through the space
     - ``tests/gates/numerics/test_axis_generator.py``
       ``::TestG4TheFourNamesAnswerThroughTheSpace``
     - over all five ``Quadrature`` classmethod factories
       (``gauss_legendre``, ``level_symmetric``, ``product``,
       ``folded_product``, ``lebedev``), with the refutation leg that a
       bare measure does not answer ``level_indices``.
   * - the mint is a section
     - ``tests/gates/numerics/test_axis_generator.py``
       ``::TestG8TheMintIsASectionOfTheForgetfulMap``
     - minting at another shape, or from an object that did not produce
       the axis, reds it.
   * - the one refusal
     - ``tests/gates/numerics/test_axis_generator.py``
       ``::TestG5TheQuadratureRefusalIsOneSpelling``
     - the message carries the axis label, the asker, ``"minted through"``
       and the arm; its fragments are disjoint from the neighbouring
       ``"axis lookup:"`` refusal; the atomless shipped axes are
       ``Declared``.
   * - each arm implies its kind and exposes its content; no ``kind``
       beside an origin
     - ``tests/gates/numerics/test_axis_identity.py`` O1, O2, O5
     - an axis constructor accepting a ``kind``, or an arm reporting
       another kind, reds.
   * - the ``FromQuadrature`` consumers are the declared ones
     - ``tests/gates/numerics/test_axis_identity.py`` O3, O4
     - O3 is an AST census of the functions holding a
       ``case FromQuadrature`` pattern against the declared four; a new
       site reds it. O4: ``HarmonicFrame.for_space`` refuses a
       ``Declared`` and a ``FromMeasure`` axis, naming itself.
   * - the spatial chain, anchored independently
     - ``tests/gates/numerics/test_space_of_axes.py``
       ``::test_the_spatial_axis_is_minted_through_the_carriers_own_measure``
     - its volumes and centres are hand-derived from the edge list, not
       read back from the mesh.
   * - the rank-:math:`d` seam
     - ``tests/gates/numerics/test_space_of_axes.py``
       ``::test_the_rank_d_spatial_axis_is_generator_less_BY_CONTRACT``
     - the contract row; inverting it is a deliberate act.
   * - the consumers read through the axis
     - the K1 and K2 keystones and the closure-mint class named in
       :ref:`spaces-generator-route-gate`
     - a consumer that reads its datum by another path stays green on the
       true axis and fails to move under the decoy.
   * - the angular contract
     - ``tests/gates/sn/angular/test_redistribution.py``
       ``::TestG9TheProtocolDeclaresWhatItsConsumersRead``
     - removing the ``level_structure`` declaration reds it.

The design record of the provenance slot before S7 (the ``generator``
field, its narrowing accessor ``generator_as``, the three answers to the
identity exclusion, the honest-``None`` inventory and the decoy catalogue)
is this page's text at ``main`` commit ``dd17f678``, lines 1002 to 2033
(``git show dd17f678:docs/theory/foundations/spaces.rst``); the milestone
is a row of :ref:`spaces-development-history`.

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
:meth:`EnergyAxis.of_edges <orpheus.numerics.axis.EnergyAxis.of_edges>`
and :meth:`EnergyAxis.synthetic
<orpheus.numerics.axis.EnergyAxis.synthetic>` — mint ``weights=None``.
``of_edges`` takes the boundary energies alone, so the numerics axis names
nothing in the data layer; a group structure is read as an axis by its own
view, :meth:`EnergyGrid.as_axis <orpheus.data.energy_grid.EnergyGrid.as_axis>`,
which lives in data beside the grid's measure and basis views (layer_dag.md
P2, 2026-10).
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


.. _spaces-metric-object:

The metric is an OBJECT: three sources, one resolution
======================================================

Everything above treats a space's metric as a **measure** — a weight per
index, stored per axis, applied by broadcasting. That is the right object
for a *measure*, and it is not the general object. A discrete Hilbert
space's metric is a symmetric positive-semi-definite **bilinear form**

.. math::

   \langle x, y\rangle_G \;=\; y^{\mathsf T} G\, x ,
   \qquad G = G^{\mathsf T},\quad G \succeq 0 ,

and nothing in that definition makes :math:`G` diagonal. Until campaign 1
phase **P7** (2026-08-30) the tree could spell exactly two realizations,
both of them a thing that is *multiplied* into the element — a broadcast
weight array (the Hadamard product :math:`G\odot x`) and the per-axis
factor measures of an axis-built space — so a form with off-diagonal
structure was **unspellable at every level**, not merely unimplemented.
P7 makes the metric a thing that is **applied**: a typed
:class:`~orpheus.numerics.metric.HilbertMetric` family owned by the
space, of which the Hadamard weight is the diagonal special case it
always was.

.. admonition:: Why this is a foundations concern and not a frame detail
   :class: tip

   Two independent consumers demanded the same missing primitive, which
   is what promoted it from a loose end to a layer:

   - a **spherical-harmonic frame on a rule that does not integrate its
     products exactly**, whose measured discrete Gram couples the modes,
     so no diagonal metric satisfies Parseval on it
     (:ref:`frame-parseval-dense-arm`). `[M]` 2026-10-18 on
     ``frame_faces.angular_frame(folded_product(4, 4), 2)`` the live block's largest
     off-diagonal is 1.000 of the Cauchy–Schwarz scale. The frame that
     first demanded the primitive was the slab one, whose density was
     ERR-080's signature and went with its repair (#429); and
   - the **curvilinear multi-moment cell mass** :math:`M/V`, whose true
     value at a spherical pole cell is
     ``[[1, 0.5], [0.5, 0.4]]``
     — dense *and* cell-dependent.

   Only the first is closed by P7, and the reason is measured rather
   than asserted: see :ref:`spaces-metric-not-on-the-axis` below.

.. _spaces-metric-three-sources:

Three sources, and the order that resolves them
-----------------------------------------------

A :class:`~orpheus.numerics.space.FunctionSpace` now has three places a
metric can come from, and **exactly one of them may be occupied**:

.. list-table:: The metric sources, in resolution order
   :header-rows: 1
   :widths: 20 34 46

   * - Source
     - Realization
     - When it is the right one
   * - ``axes``
     - the axes-DERIVED object — one
       :class:`~orpheus.numerics.metric.FactoredMetric` with a
       :class:`~orpheus.numerics.metric.DiagonalMetric` per weighted axis,
       each factor's weight vector applied on its own index block (CS4c
       step 6 item 6.2c-i, 2026-09-08; until then an inline per-axis loop,
       ``FunctionSpace._apply_axes_weights``, spelled the same arithmetic a
       second time)
     - An axis-built space (:eq:`spaces-axis-product`). The axes **are**
       the metric source; the space never consults the weights slot, and
       consults ``metric`` only as the positioned OVERLAY of forms
       described in the next row.
   * - ``metric``
     - the object itself —
       :class:`~orpheus.numerics.metric.DiagonalMetric`,
       :class:`~orpheus.numerics.metric.DenseMetric`, or a
       :class:`~orpheus.numerics.metric.FactoredMetric` product of them.
       Beside ``axes`` (item 6.2c-i, ruling R-6.2c-2) it is admitted only
       as the positioned OVERLAY of forms — a ``FactoredMetric`` with one
       entry per axis, a form only on a block whose axis carries no
       measure, never a ``DiagonalMetric`` (a diagonal measure is the
       axis's own to carry) — and is MERGED into the axes-derived object,
       never substituted for it: every other axis keeps supplying its
       own block, and a product of axis-built factors carries its
       factors' overlays concatenated beside its axes.
     - A form no weight array can spell. The founding occupant is the
       frame's matrix Parseval metric :math:`G^{+}`, and since item
       6.2c-ii (2026-09-08) it rides an AXIS-BUILT space: on a
       dense-Gram frame the dressing empties the head axis's measure and
       positions :math:`G^{+}` on the space's derived object
       (:ref:`spaces-moment-head-axis-built`). `[M]` **17 of 75** shipped
       (rule, :math:`L`) rows measure ``DENSE``, **2 of them at**
       :math:`L \le 2` (``gauss_legendre(2)`` and ``folded_product(2,4)``,
       both at :math:`L = 2`) — so this is a shipped path, not a corner.
   * - ``inner_product_weights``
     - resolved to a
       :class:`~orpheus.numerics.metric.DiagonalMetric`, whose
       arithmetic is operation-for-operation the arms that used to live
       inline
     - The legacy dense weight array — every pre-P7 metric, and still
       the ordinary case.
   * - *(none)*
     - the verbs short-circuit without constructing an object
     - Euclidean, :math:`\sum_i x_i y_i`.

Resolution happens once per space, in
``FunctionSpace._resolved_metric``, and the three metric verbs —
:meth:`~orpheus.numerics.space.FunctionSpace.apply_metric`,
:meth:`~orpheus.numerics.space.FunctionSpace.apply_inverse_metric` and
:meth:`~orpheus.numerics.space.FunctionSpace.inner_product` — delegate
to whatever it returns. Composite spaces
(:class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`,
:class:`~orpheus.numerics.coupled_system.CoupledSpace`) override the
three verbs to dispatch **per block** and therefore inherit the family
for free: a dense-metric interior block flows through a composite with
no composite-level change.

**Exclusivity is a construction guard with three pairwise arms**
(``axes``/``weights``, ``axes``/``metric``, ``weights``/``metric``), and
the restructuring is not cosmetic. Before P7 the check sat *behind* an
``if self.axes is None: return`` early exit, so the
``(weights, metric)`` arm was **structurally unreachable** — a guard arm
with no possible witness (``vv-principles`` #17's granularity trap, in
its purest form: the guard would have been mutated as a unit and
certified by whichever arm the suite happened to reach). Each arm now
carries its own gate in ``tests/gates/numerics/test_space_of_axes.py``.

.. warning::

   **The flip found a whole subclass family standing outside the
   construction contract.** All three
   :class:`~orpheus.numerics.space.FunctionSpace` subclasses that
   defined their own ``__post_init__`` at P7 —
   ``SphericalHarmonicSpace`` (retired at S7), the since-retired
   ``SpatialMomentSpace``, and
   :class:`~orpheus.numerics.spaces.angular_trace_space.AngularFaceTraceSpace`
   — **overrode without chaining**, so *every* base construction guard
   (not only the new one) had been silently opted out of since each
   subclass was written. All three were made to call
   ``super().__post_init__()`` first, and the rule held as the family
   changed. `[M]` 2026-10-22, a runtime walk of
   ``FunctionSpace.__subclasses__`` after importing every ``orpheus``
   module: 3 of the 10 subclasses define ``__post_init__``
   (``DualSpace``, ``CoupledSpace``, ``AngularFaceTraceSpace``), and all
   3 call ``super().__post_init__()``. It surfaced only because the
   first dressing attempt left the
   spherical-harmonic basis's continuum weights in the field *beside*
   the new metric object and nothing refused it — the two-source illegal
   state was representable exactly where the guard was supposed to live.
   This is the standing argument for making a guard's own arms
   individually falsifiable: the arm that has no witness is the arm that
   is not running.

.. important::

   **⚠ An empty weights slot no longer means "Euclidean".**
   ``inner_product_weights is None`` describes the *diagonal source*
   only; a dense-metric space reads ``None`` there while carrying a
   real, non-Euclidean metric.
   Read the space's behaviour through
   :meth:`~orpheus.numerics.space.FunctionSpace.inner_product`, never
   through the slot. `[M]` re-measured 2026-09-07 after item 6.2a, and
   the production exposure is still small — worth stating so the hazard
   is not over-read: of the **221** lines across **57** files in
   ``orpheus/`` + ``tests/`` that name the slot, **28** carry an
   ``is None`` / ``is not None`` spelling (prose included) and only
   **four** are production BRANCHES — all four inside ``space.py``
   itself: the two arms of the exclusivity guard, the resolution in
   ``_resolved_metric``, and — fourth since 6.2a — the dense-slot-leaf
   arm of the tensor product's factored builder, which positions a
   hand-built factor's weights as a
   :class:`~orpheus.numerics.metric.DiagonalMetric` on that factor's
   block. (This paragraph read *"198 lines across 52 files … 20 are*
   ``is None`` / ``is not None`` *branches and only three are
   production"* until then. The fourth is a change of SPELLING, not of
   exposure: the pre-6.2a builder funnelled both a dense-slot leaf and
   an axis-built factor through one local ``w is not None`` test — the
   latter via the densifier — where the arms are now separate and the
   densifier is gone.)

.. _spaces-metric-moore-penrose:

The Moore–Penrose doctrine, extended from a reciprocal to a matrix
-------------------------------------------------------------------

:meth:`HilbertMetric.apply_inverse
<orpheus.numerics.metric.HilbertMetric.apply_inverse>` is the
**pseudo**-inverse :math:`G^{+}` — the reciprocal on the metric's range,
zero on its kernel — for every realization. On the diagonal arm this is
the doctrine the tree has carried since the trace-metric work (a
tangential ordinate with :math:`|\Omega\cdot n| = 0` has weight zero and
must map to zero, not to infinity). The dense arm extends it unchanged,
and there it is **not a convenience: it is the only thing that exists.**

`[M]` 2026-10-18 on ``frame_faces.angular_frame(Quadrature.folded_product(4, 4), 2)``,
the mirror-even spherical harmonics on :math:`S^2/\sigma_y` and the
smallest of the three dense angular heads production builds (the
frame-bridge census of 2026-10-18), the discrete Gram is :math:`15\times15` with **6 live diagonal slots and rank 5**. Its
singular values are :math:`4\pi = 12.566,\ 4.189,\ 4.189,\ 2.513,\
2.513` and then round-off (the largest :math:`2.7\times10^{-16}`, sixteen
orders below the smallest genuine mode). The dropped mode is structural,
not small: its singular vector, a combination of two :math:`\ell = 2`
columns, annihilates the evaluation table on the rule's 8 nodes
(:math:`\|Av\|_\infty = 2.8\times10^{-16}`), so the rule cannot tell those
two harmonics apart. :func:`numpy.linalg.inv` raises there;
:func:`numpy.linalg.pinv` returns the object Parseval needs. The other two
dense heads are singular too: ``frame_faces.angular_frame(folded_product(4, 4), 4)``,
:math:`45\times45` with 15 live slots and rank 8, and
``frame_faces.angular_frame(product(4, 8), 4)``, :math:`45\times45` with 25 live slots
and rank 23. (Probes ``scratch/metric_docstring/probe_flagship.py``,
``probe_kernel.py`` and ``probe_theory.py``, numpy 2.4.4. Until the
ERR-080 repair of 2026-09-02 the example here was the slab
``frame_faces.angular_frame(gauss_legendre(8), 2)``, 5 live slots and rank 4; that
frame's Gram is now diagonal, :math:`3\times3` and rank 3.)

Two construction choices in
:class:`~orpheus.numerics.metric.DenseMetric` are worth their reasons:

#. **The cutoff is PINNED, not left to the library default.** An
   implicit default is a silent dependency on a NumPy version, and a
   mutation arm that over-truncates is only meaningful against a fixed
   value. The shipped ``rcond`` is :math:`10^{-12}`, and the cliff it
   must sit below is computable rather than guessed: NumPy's ``rcond``
   is relative to the largest singular value, so truncation begins at
   :math:`\sigma_{\min}^{\rm live}/\sigma_{\max} = 2.513 / 12.566 = 0.2`.
   `[M]` scanned on that Gram, the draw-free range of the Parseval ratio
   over the band-limited subspace is :math:`[1, 1]` to nine digits at
   every scanned ``rcond`` from :math:`10^{-15}` to :math:`0.19`, as it
   must be, since no truncation can occur below the cliff. From
   :math:`0.21` upward it is :math:`[0, 1]`: the genuine :math:`2.513`
   pair is discarded and a field living in it reads zero.
   :math:`10^{-12}` sits 11.3 decades below the cliff and 3.8 above the
   largest round-off residue as ``eigh`` reports it
   (:math:`1.67\times10^{-16}` relative), and that lower margin moves
   with the LAPACK driver (the window and its census are the comment on
   ``_DENSE_METRIC_RCOND`` in :mod:`orpheus.numerics.metric`).
#. **The pseudo-inversion declares hermiticity.** `[M]` 2026-10-18 over
   six dense frames (``scratch/metric_docstring/probe_hermitian.py``),
   the same call *without* ``hermitian=True`` returns a matrix whose
   asymmetry :math:`\max|M - M^{\mathsf T}|` ranges from
   :math:`4.6\times10^{-17}` to :math:`1.75\times10^{-15}` (the largest
   on ``frame_faces.angular_frame(product(4, 4), 2)``), against at most
   :math:`6.9\times10^{-17}` with it: one to two orders, and far below
   the type's symmetry admission of :math:`10^{-12}`. The declaration is
   kept because it makes the pseudo-inverse cut on ``eigh``, the
   decomposition the cutoff window is measured on. (On the pre-repair
   slab Gram of 2026-08-30 the two read :math:`4.74\times10^{-14}` and
   :math:`1.1\times10^{-16}`.)

The type refuses two illegal states outright rather than repairing them:
an **asymmetric** matrix (an asymmetric form is not an inner product —
producers symmetrize, the type does not), and an explicitly supplied
inverse face that fails the **Penrose identity** (an inconsistent
:math:`(G, G^{+})` pair is not a representation choice). Both thresholds
are module constants, so a gate on the type quotes the type's own
number (``vv-principles`` #16).

.. _spaces-metric-parseval-theorem:

Parseval is a THEOREM for any Gram, singular or not
----------------------------------------------------

The pseudo-inverse is not a graceful degradation of the invertible
story; it is what makes the story *general*. For a band-limited field
:math:`\psi = S_0 c` the analysis face returns the **covariant** moments
:math:`\varphi = M\psi = Gc` identically
(:eq:`frame-analysis-is-the-gram`, pure algebra). Measuring those
moments under :math:`G^{+}` and using the first Penrose identity
:math:`G G^{+} G = G`:

.. math::
   :label: spaces-pseudo-inverse-parseval

   \|M\psi\|^2_{G^{+}}
   \;=\; (Gc)^{\mathsf T} G^{+} (Gc)
   \;=\; c^{\mathsf T}\bigl(G G^{+} G\bigr) c
   \;=\; c^{\mathsf T} G\, c
   \;=\; \|S_0 c\|^2_W .

.. (vv-status rationale) A representational identity: it derives the
   Parseval isometry (:eq:`frame-parseval-isometry`, whose SSOT is
   :ref:`frame-parseval-metric` on the frame page) for a SINGULAR Gram
   by substituting the first Penrose identity. No flux, no eigenvalue,
   no discretization error — the metric-object layer's statement of a
   claim the frame page owns in the frame register. The verifiable
   content is the DenseMetric law battery
   (``tests/gates/numerics/test_dense_metric.py``, hand-derived literals in
   exact binary fractions plus the range-projector leg) and, on the
   frame side, the six-mechanism dense dressing gate, the isometry
   gate's three dense rows, and the wrong-metric discriminator
   (``tests/gates/numerics/test_frame.py``).
.. vv-status: spaces-pseudo-inverse-parseval documented

.. no-implementation:: spaces-pseudo-inverse-parseval
   :kind: identity

   **Nothing implements this.** Both sides are computed — the left by
   :meth:`FunctionSpace.inner_product
   <orpheus.numerics.space.FunctionSpace.inner_product>` on the frame's
   dressed coefficient space, the right by the same verb on the measure
   space — and no line in production forms the comparison. That is the
   point: the equality is what lets the dressing be a one-line
   ``replace`` rather than a bespoke correction. It is *measured* by the
   Parseval isometry gate.

The identity holds because :math:`G G^{+} G = G` is one of the four
Penrose conditions, true by definition of the pseudo-inverse — no
hypothesis about rank, conditioning or invertibility enters. `[M]`
2026-10-18 on ``frame_faces.angular_frame(folded_product(4, 4), 2)``
(``scratch/metric_docstring/probe_theory.py``) the residual is
:math:`\max|G G^{+} G - G| = 1.78\times10^{-15}` absolute, which is
:math:`1.41\times10^{-16}` relative to :math:`\max|G| = 4\pi`, and
:math:`1.71\times10^{-16}` in the Frobenius ratio. Name the norm when
quoting it: the absolute figure is an order larger only because
:math:`\max|G|` is :math:`4\pi`.

.. _spaces-metric-one-spelling:

The pairing has ONE spelling, and that is a bit-identity contract
------------------------------------------------------------------

:meth:`HilbertMetric.pairing
<orpheus.numerics.metric.HilbertMetric.pairing>` is
``float(np.sum(self.apply(x) * y))`` — :math:`\sum (Gx)\odot y` — for
**every** realization, inherited rather than overridden. Two reasons,
both load-bearing:

**Bit-identity with the legacy path.** The pre-P7 diagonal spelling was
``np.sum(w * x * y)``, which Python evaluates left-to-right as
``(w*x)*y``; ``DiagonalMetric.apply(x)`` *is* ``w_b * x``, so routing
the shipped diagonal path through the family preserves the reduction
tree exactly. The obvious alternative — densify and matmul,
:math:`y^{\mathsf T}(\operatorname{diag}(w)\,x)` — is **not** the same
program. `[M]` 40 seeds :math:`\times` 500 draws at :math:`n = 15`, the
two spellings disagree bitwise on **60.4 %–69.8 %** of draws, with the
worst gap ranging over **46–16 384 ULP** (relative
:math:`9.2\times10^{-15}`–:math:`2.1\times10^{-12}`) across seeds.
Note which half of that is publishable: the *fraction* is a stable
property of the arithmetic (two different reduction trees disagree on
about two draws in three), while the *ULP gap* is unbounded wherever the
sum nearly cancels and must never be frozen into a docstring or a gate
band. What the measurement licenses is the structural claim: a matmul
pairing would move pinned numbers tree-wide, so it is deliberately
unspellable here.

**One source of truth for the adjoint.**
``AdjointOperator`` builds :math:`A^{\dagger} = G^{-1}A^{\mathsf T}G`
from ``apply_metric``/``apply_inverse_metric`` while the pairing that
*judges* it comes from ``inner_product``. Deriving the pairing from
``apply`` makes those two agree by construction — the ERR-067 family
(two spellings of one metric diverging silently) becomes unspellable
rather than merely untested.

.. _spaces-riesz-legs:

The Riesz legs: the metric as two arrows, and the adjoint as their composition
------------------------------------------------------------------------------

The metric object is applied through three verbs on the space
(:meth:`~orpheus.numerics.space.FunctionSpace.apply_metric`,
:meth:`~orpheus.numerics.space.FunctionSpace.apply_inverse_metric`,
:meth:`~orpheus.numerics.space.FunctionSpace.inner_product`). Two of them
are also **arrows**: operators with a declared domain and codomain, which
the operator algebra can compose like any other operator. They are the
**Riesz legs**. This section is their one definition.

**The mathematics.** A finite-dimensional real Hilbert space :math:`V`
with symmetric positive-definite metric :math:`G` is isomorphic to its
dual :math:`V^{*}` (the space of linear functionals on :math:`V`), and
the isomorphism is fixed by the metric. The Riesz representation theorem
states it: every functional :math:`f \in V^{*}` is the pairing with
exactly one vector, :math:`f(x) = \langle v_f, x\rangle_G` for all
:math:`x` (:cite:`Kreyszig1978` §3.8, Theorem 3.8-1, p. 188). In coordinates the map from a
vector to the functional it represents multiplies by :math:`G`, and the
inverse map multiplies by :math:`G^{-1}`. Differential geometry names the
same pair the **musical isomorphisms** of a metric, :math:`\flat`
("flat", which lowers an index) and :math:`\sharp` ("sharp", which raises
one) (:cite:`Lee2018` Ch. 2, "Raising and Lowering Indices", p. 26), and
ORPHEUS adopts those names. Lee's :math:`\sharp` is the inverse of
:math:`\flat` and needs a nonsingular metric; the raising leg below
extends it to a singular one through the pseudo-inverse, which is
ORPHEUS's choice and not Lee's:

.. (V&V scope note) Definitional: it STATES what the two leg classes
   apply (RieszLowerOperator.apply delegates to FunctionSpace.apply_metric,
   RieszRaiseOperator.apply to FunctionSpace.apply_inverse_metric), plus
   the arrow bookkeeping. Not a solver claim: no flux, no eigenvalue, no
   discretization error. Wired, L0, to the test a mutation of the
   lowering leg reddens: tests/gates/numerics/test_riesz_legs.py
   ::test_lowering_realizes_the_metric_pairing (G-H2, the defining
   pairing <flat x, y>_Euclid == <x, y>_G against a hand-spelled weight
   array, bitwise). The raising half's Moore-Penrose behaviour is pinned
   by G-H1's singular row (a 1/G raising leg reddens it and nothing else
   in the module), which is the witness of spaces-riesz-round-trip below
   and carries no marker. The dual-space admission row (G-H4,
   test_the_legs_admit_a_dual_space) reads the dual metric, not a leg's
   arithmetic, so it carries no marker either.

.. math::
   :label: spaces-riesz-lower-raise

   \flat_V : V \to V^{*}, \quad \flat_V\,x = G\,x ,
   \qquad\qquad
   \sharp_V : V^{*} \to V, \quad \sharp_V\,f = G^{+} f ,

.. implements:: spaces-riesz-lower-raise
   :by: orpheus.numerics.operator.RieszLowerOperator.apply

   **Implemented by** the lowering leg, whose ``apply`` is
   ``space.apply_metric(x)``: the leg contributes the ARROW (``domain`` the
   primal :math:`V`, ``codomain`` ``V.dual()``) and borrows the arithmetic
   from the space, so the metric keeps one home.

.. implements:: spaces-riesz-lower-raise
   :by: orpheus.numerics.operator.RieszRaiseOperator.apply

   **Implemented by** the raising leg, whose ``apply`` is
   ``space.apply_inverse_metric(f)``: the Moore–Penrose pseudo-inverse of
   :ref:`spaces-metric-moore-penrose`, on every metric realization.

where :math:`G^{+}` is the Moore–Penrose pseudo-inverse, the unique
matrix satisfying the four Penrose conditions :math:`G G^{+} G = G`,
:math:`G^{+} G G^{+} = G^{+}`, :math:`(G G^{+})^{\mathsf T} = G G^{+}`
and :math:`(G^{+} G)^{\mathsf T} = G^{+} G` (:cite:`Penrose1955` Theorem 1,
p. 406, stated there for the conjugate transpose; the inverse was first
defined by :cite:`Moore1920`, through conditions on column and row
spaces rather than these four equations). On an
invertible metric :math:`G^{+} = G^{-1}` and the pair is the Riesz map and
its inverse exactly. The defining property of the lowering leg is that
the Euclidean pairing of :math:`\flat x` with :math:`y` IS the metric
pairing, :math:`\langle \flat x, y\rangle = y^{\mathsf T} G x =
\langle x, y\rangle_G`, which is what the test G-H2
(``test_lowering_realizes_the_metric_pairing``) asserts bitwise.

Three facts about the dual end are load-bearing:

- **The dual space carries the inverse form.** :math:`V^{*}` has
  :math:`V`'s coordinates and the inverse Riesz form :math:`G^{+}` as its
  metric, :math:`\langle f, g\rangle_{V^{*}} = f^{\mathsf T} G^{+} g`
  (:ref:`spaces-dual-space`). That is the pairing the Riesz map induces on
  the dual: :math:`\langle \flat x, \flat y\rangle_{V^{*}} = x^{\mathsf
  T} G G^{+} G y = \langle x, y\rangle_V`. It is a distinct space,
  ``V.dual() != V``, so a composition chain tells a vector from a
  functional.
- **The legs are posed on every space, primal or dual.** On :math:`V^{*}`
  the lowering leg :math:`\flat_{V^{*}} : V^{*} \to V^{**} = V` applies
  :math:`G^{+}` and the raising leg :math:`\sharp_{V^{*}}` applies
  :math:`G`, so the double lowering is the identity on the metric's range:
  for :math:`w = [0.5, 2, 4]` and :math:`x = \mathbf 1`,
  ``RieszLowerOperator(V.dual()).apply(RieszLowerOperator(V).apply(x))``
  returns :math:`x` exactly (`[M]` 2026-10-21, the archivist's probe on
  the S6 worktree; gate G-H4, ``test_the_legs_admit_a_dual_space``, and
  ``test_dual_metric.py`` L1 and L2).
- **Each leg is its own transpose, and the adjoint of the lowering leg is
  the raising leg.** Every shipped metric realization is symmetric (a
  diagonal weight, or a :class:`~orpheus.numerics.metric.DenseMetric`
  admitted only through its symmetry guard, or a
  :class:`~orpheus.numerics.metric.FactoredMetric` built from such
  factors), so :math:`\flat^{\mathsf T} = \flat` and
  :math:`\sharp^{\mathsf T} = \sharp` under the reflexive identification
  :math:`V^{**} = V`, and both legs report
  :attr:`~orpheus.numerics.operator.LinearOperator.is_adjointable`, so
  their ``.dual()`` runs. Both ends of a leg are declared spaces, so its
  adjoint is :func:`~orpheus.numerics.operator.posed`, and the Hilbert
  adjoint is the law :math:`\flat_V^{*} = \sharp_V`. With
  :math:`W = V^{*}` in :eq:`spaces-adjoint-riesz-composition`,
  :math:`\flat_V^{*} = \sharp_V \circ \flat_V^{\mathsf T} \circ
  \flat_{V^{*}} = G^{+}\,G\,G^{+} = G^{+}` (the second Penrose
  condition), which is what ``RieszRaiseOperator(V)`` applies (`[M]` 2026-10-21: on the weights
  above ``RieszLowerOperator(V).H.apply(y)`` equals
  ``RieszRaiseOperator(V).apply(y)`` to the printed digits; gate
  ``test_dual_metric.py`` L3, by reciprocity).

The round trip, and the Moore–Penrose choice on the kernel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A metric that is only positive **semi**-definite has a kernel, and the
S\ :sub:`N` boundary trace carries one: the trace metric is
:math:`|\Omega\cdot\hat n|\,w_n`, which vanishes on the tangential
ordinates (:ref:`bc-trace-structure`; on one face the kernel is exactly
the rows neither :ref:`half-trace <bc-half-trace>` carries). Composing the
legs then does not return the identity. From the third and fourth Penrose conditions, with
:math:`G` symmetric:

.. (vv-status rationale) An identity of the pseudo-inverse, not a solver
   claim: G⁺G is the orthogonal projector onto range(G) for any symmetric
   G (Penrose conditions 1, 2 and 4). Implemented by
   FunctionSpace.project_live (S6b, ruling 74), which forms exactly
   apply_inverse_metric(apply_metric(x)). Its measurement is
   tests/gates/numerics/test_riesz_legs.py G-H1: identity on a strictly
   positive metric (the positivity precondition asserted), the
   kernel-zeroing projector on a singular diagonal one, and non-finite
   values refused. [M] 2026-09-22 a 1/G raising leg reddens exactly the
   singular row, 1 of the module's 11 tests. The projector's own gates are
   tests/gates/numerics/test_quotient_metric.py B3, B3b and B4 (♯∘♭ = P
   bitwise).
.. vv-status: spaces-riesz-round-trip documented

.. math::
   :label: spaces-riesz-round-trip

   \sharp_V \circ \flat_V \;=\; G^{+} G \;=\; P_{\operatorname{range}(G)} ,

.. implements:: spaces-riesz-round-trip
   :by: orpheus.numerics.space.FunctionSpace.project_live

   **Implemented by** the live projector, whose body is
   ``self.apply_inverse_metric(self.apply_metric(x))``: the round trip of
   the two legs' arithmetic, read through the space so a composite applies
   it block by block (:ref:`spaces-quotient`). The composite adjoint uses
   the legs in the other order around an operator (below). The identity
   is what makes the legs *honest* on a singular metric, and it is
   measured by G-H1.

the orthogonal projector onto the range of :math:`G`. It is the identity
exactly when :math:`G` is strictly positive. On a diagonal metric with
zero weights it zeroes those slots and keeps the others, which on a trace
block is the **tangential-zeroing projector**. On a dense singular metric
it is not diagonal at all: `[M]` 2026-09-22, a rank-2 symmetric
:math:`G = BB^{\mathsf T}` on :math:`\mathbb R^4` (``B`` a seeded
:math:`4\times 2` Gaussian draw) through
:class:`~orpheus.numerics.metric.DenseMetric`, the round trip built
column by column is symmetric, idempotent and of rank 2, and has
off-diagonal entries of order :math:`10^{-1}` (the probe is
``scratch/_definitions/riesz_probe.py``). It also equals
``pinv(G) @ G`` to the last bit, which is not an independent check:
:class:`~orpheus.numerics.metric.DenseMetric` stores
``pinv(G, hermitian=True)`` as its inverse face. A gate asserting
``♯∘♭ == id`` would therefore be blind to a ``1/G`` raising leg on every
strictly positive fixture and a false red on a legal trace space, which
is why G-H1 carries a singular fixture beside the positive one: `[M]`
2026-09-22, a ``1/G`` raising leg installed in-process reddens exactly
that singular row, 1 of the module's 11 tests, and the unmutated run
reads 11 of 11 green.

**Why the pseudo-inverse and not some other inverse.** The kernel is
exactly where the metric cannot see: a vector :math:`z \in \ker G`
pairs to zero with everything. So wherever :math:`\sharp` is used to
solve for a vector from its pairings, that vector is determined only up
to :math:`\ker G`, and :math:`G^{+}` picks one representative, the one
with no component in :math:`\ker G`. Because :math:`G` is symmetric,
:math:`\operatorname{range}(G)` and :math:`\ker G` are Euclidean-orthogonal,
so this is also the representative of least Euclidean norm. The next
subsection shows the choice is exact for the adjoint whenever the
operator does not carry kernel content into what the codomain metric
measures. :math:`G^{-1}` does not exist there (:func:`numpy.linalg.inv`
raises on every dense angular Gram the tree builds,
:ref:`spaces-metric-moore-penrose`),
and a reciprocal :math:`1/w` would put infinities in the tangential
slots.

The adjoint is the legs composed around the transpose
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For an operator :math:`A : V \to W` between spaces with metrics
:math:`G_V` and :math:`G_W`, the Hilbert adjoint :math:`A^{*} : W \to V`
is defined by reciprocity (:cite:`Kreyszig1978` §3.9, Definition 3.9-1,
p. 196, for genuine Hilbert spaces, where it exists and is unique):

.. math::

   \langle A x, y\rangle_W \;=\; \langle x, A^{*} y\rangle_V
   \qquad \text{for all } x \in V,\; y \in W .

Here :math:`W` is the codomain *space*, the spelling the code uses; the
subscript :math:`W` on a norm elsewhere on this page, as in
:eq:`spaces-pseudo-inverse-parseval`, names the quadrature-weighted
pairing instead. In coordinates the left side is :math:`(Ax)^{\mathsf T} G_W y =
x^{\mathsf T}\,(A^{\mathsf T} G_W y)` and the right side is
:math:`x^{\mathsf T} G_V\,(A^{*}y)`. Equal for every :math:`x` means
:math:`G_V\,(A^{*}y) = A^{\mathsf T} G_W\, y`. Solving with the
pseudo-inverse and reading the product right to left gives three arrows:

.. (V&V scope note) The arrow-typed factorization of the Hilbert
   adjoint: with V = W and G⁻¹ read as the pseudo-inverse it IS
   g-adjoint-definition, one quantity in two spellings (X4), not a
   second definition. Not a solver claim. Wired, L0, to the two tests a
   mutation of either leg reddens, both in
   tests/gates/numerics/test_riesz_legs.py (G-A1):
   ::test_adjoint_apply_is_bit_identical_to_the_inline_formula (the
   composition bit-identical to G_V⁺ Aᵀ(G_W y)) and
   ::test_adjoint_transpose_is_the_leg_theorem ((A*)ᵀ = ♭_W A ♯_V
   bit-identical to a hand-built right side, with a wrong-order negative
   leg). G-H3 (counting spies proving the composite ROUTES through the
   legs) and the dagger-dual square stay green under a leg mutation and
   carry no marker. The per-leg mutation battery
   tests/gates/sn/architecture/test_monomorphic_leaves.py
   ::test_each_riesz_leg_is_individually_load_bearing exercises the same
   legs on the SN leaves.

.. math::
   :label: spaces-adjoint-riesz-composition

   A^{*} \;=\; \sharp_V \circ A^{\mathsf T} \circ \flat_W ,
   \qquad
   W \xrightarrow{\;\flat_W\;} W^{*}
     \xrightarrow{\;A^{\mathsf T}\;} V^{*}
     \xrightarrow{\;\sharp_V\;} V ,

.. implements:: spaces-adjoint-riesz-composition
   :by: orpheus.numerics.operator.AdjointOperator

   **Implemented by** the dagger arrow, which builds the three factors AT
   CONSTRUCTION from the inner operator's two ends: the codomain's
   lowering leg, the metric-free dual arrow (the object ``inner.dual()``
   returns, constructed directly), and the domain's raising leg;
   ``apply`` is their composition and nothing else.

.. implements:: spaces-adjoint-riesz-composition
   :by: orpheus.numerics.operator.LinearOperator.dual

   **Implemented by** the middle factor's public mint: the dual arrow
   :math:`A^{\mathsf T} : W^{*} \to V^{*}` is the representation
   transpose carried between the dual spaces, and every metric lives in
   the legs, never here. The adjoint does not call this verb; it
   constructs the same class directly.

read as: lower the output-side vector to a functional on :math:`W`, pull
that functional back along :math:`A` (the metric-free transpose, between
the dual spaces), and raise the result to a vector of :math:`V`. With
:math:`V = W` and :math:`G^{-1}` read as :math:`G^{+}` this is
:eq:`g-adjoint-definition`, :math:`A^{\dagger} = G^{-1}A^{\mathsf T}G`;
the two labels are one quantity, the composite page's endomorphic
spelling and this page's arrow-typed one, and the second is the one that
generalises to :math:`V \ne W`.

.. _spaces-adjoint-kernel-condition:

**When the kernel is non-trivial the reciprocity is exact under one
condition.** Substituting the pseudo-inverse and using
:math:`G_V G_V^{+} = G_V^{+} G_V = P_{\operatorname{range}(G_V)}`,

.. math::

   \langle A x, y\rangle_W - \langle x, A^{*} y\rangle_V
   \;=\; x^{\mathsf T}\bigl(I - G_V G_V^{+}\bigr) A^{\mathsf T} G_W\, y
   \;=\; \bigl\langle A\,P_{\ker G_V}\,x,\; y\bigr\rangle_W ,

so the identity holds for every :math:`x` and :math:`y` if and only if
:math:`A` maps :math:`\ker G_V` into :math:`\ker G_W`, and it holds for
every :math:`x \in \operatorname{range}(G_V)` unconditionally. `[M]`
2026-09-22, twice. In plain numpy (``scratch/_definitions/riesz_probe.py``,
one draw), a random :math:`A : \mathbb R^4 \to \mathbb R^3` against the
rank-2 :math:`G_V` above: the defect equals the predicted
:math:`\langle A P_{\ker} x, y\rangle_W` to a relative
:math:`1.8\times10^{-15}` and vanishes to :math:`4\times10^{-15}` on
:math:`x \in \operatorname{range}(G_V)`. Through the production arrow
(``A.H`` on a space carrying a rank-2
:class:`~orpheus.numerics.metric.DenseMetric`,
``scratch/_definitions/riesz_production_probe.py``), over 200 seeds: the
defect matches the prediction to a median relative
:math:`3.5\times10^{-16}` and at most :math:`2.1\times10^{-13}`, and on
:math:`\operatorname{range}(G_V)` it is at most
:math:`7.4\times10^{-14}` of the pairing; a positive-definite control
metric reads a defect at round-off.

On the S\ :sub:`N` trace block the kernel is the tangential rows, and
the condition holds there exactly, though not only through the trace
selectors. `[M]` 2026-09-22 (#493) on a 2-D Cartesian
``Quadrature.product(4, 4)`` problem (a non-square :math:`3\times4`
mesh, two groups, two materials, posed through
:class:`~orpheus.sn.problem.SNProblem`), where 224 of the composite's
832 degrees of freedom are tangential (on every face, 8 of the 16
ordinates graze it): the loss :math:`A = (L + C) - S - N_{2n} - B` maps
each of the 224 tangential unit vectors to exactly zero and writes
nothing into a tangential row, so the predicted defect
:math:`\langle A P_{\ker} x, y\rangle_W` is bitwise zero on 20 of 20
seeds under vacuum, reflective and mixed boundaries, and the measured
reciprocity defect is at most :math:`1.5\times10^{-14}` of
:math:`|\langle A x, y\rangle_W|` (the linear-discontinuous closure
reads the same: 448 tangential columns exactly zero, defect at most
:math:`2.6\times10^{-14}`). Two mechanisms carry it. The boundary law
:math:`B` reads and writes the trace only through the selectors of
:eq:`bc-trace-restriction-pair`, whose index sets exclude the
tangential rows. The streaming operator :math:`L` writes the trace only
through the same selectors, but it does **read** the tangential rows:
the multi-dimensional sweep walk carries a grazing ordinate (one with
:math:`\mu_a = 0` along an axis :math:`a`) as if it moved in the
positive :math:`a` direction, and takes its inflow from the whole slot
of the lower :math:`a` face rather than through a selector. That read
never reaches a weighted row, because it is multiplied by the ordinate's
streaming coefficient :math:`|\mu_a|/\Delta_a`, which is exactly
``0.0``: every shipped rule's tangential cosines are exact zeros, and
the trace metric and the streaming stencil read the same cosine. A
coefficient :math:`\delta` (relative to the largest one) in its place
leaves the streaming operator a defect of about :math:`15\,\delta` of
the pairing (measured at :math:`\delta = 10^{-3}` and
:math:`10^{-13}`), so the exactness depends on that cosine being
exactly zero. Writing into a tangential row, by contrast, can never
break the identity, because the codomain metric gives those rows zero
weight: that every output of the loss is zero on the tangential rows is
true (:ref:`g-adjoint`, Key Facts), but it is not the reason for the
exactness. Two gates in
``tests/gates/sn/operators/test_g_adjoint_reciprocity.py`` pin it on the same
``product(4, 4)`` box (with mixtures of their own), under vacuum and
reflective boundaries. The
claim-level gate is the two ``cart2d_product44`` rows of
``test_g_adjoint_reciprocity_full_block``, which read the reciprocity
functional itself and, their docstring records, redden when the grazing
streaming coefficient is :math:`10^{-8}` of the largest but not at
:math:`10^{-13}`. The structural gate is
``test_tangential_trace_slots_are_a_zero_summand_of_the_loss``, which
asserts bitwise that the loss maps every tangential unit vector to zero
and writes nothing into a tangential row; it resolves the coefficient
down to :math:`10^{-13}` and also sees the two defects the functional
is blind to, a write into a tangential row and a tangential ordinate
misclassified as inflow or outflow.
`[REFUTED 2026-09-22]` for the question "what makes the
S\ :sub:`N` adjoint exact": this paragraph said until then that an
operator reaching the trace block only through the selectors "can
neither read nor write those rows"; the fact the measurement
establishes is that :math:`L` reads them and multiplies the read by an
exact zero.

**What the arrow form buys over the arithmetic.** The composition
:math:`G_V^{+}\,A^{\mathsf T}(G_W\,y)` could be, and until CS4c step 1
(2026-08-30) was, written inline. Promoting the factors to arrows is what
the operator algebra needs from them, for four reasons:

#. **The metric keeps one home.** Each leg delegates to the space's
   resolved :class:`~orpheus.numerics.metric.HilbertMetric`, so the
   adjoint serves a diagonal weight, a dense Gram and a composite
   bulk :math:`\oplus` trace metric (a singular trace block included)
   through the same three calls; G-A1 pins the leg composition
   bit-identical to the pre-leg inline formula (`[M]` 2026-09-22, the
   eleven tests of ``test_riesz_legs.py`` and the 40 rows of the per-leg
   battery pass under ``python -O -m pytest``).
#. **The dagger laws become structure.** ``A.H.H is A`` is an object
   identity (the adjoint of the adjoint returns the inner operator), and
   the transpose of the adjoint is a theorem of the legs,
   :math:`(A^{*})^{\mathsf T} = \flat_W \circ A \circ \sharp_V`, since
   the metrics are symmetric. That theorem is what made ``A.H``
   adjointable and closed #375.
#. **Each factor is a separately mutable seam.** A mutation that drops
   the metric from both legs at once is a similarity transformation, and
   a similarity is invisible wherever the operator commutes with the
   metric; a mutation that drops one leg is not. Splitting the paired
   mutation into single legs moved the ledger battery from 9 of 20 to
   20 of 20 rows red (the SN development history's Riesz-legs entry,
   :ref:`sn-development-history`; the battery re-measures itself in
   ``test_monomorphic_leaves.py::test_each_riesz_leg_is_individually_load_bearing``).
#. **The adjoint needs both ends, so an operator without them has none.**
   The legs are built from ``inner.domain`` and ``inner.codomain``. An
   operator that has not declared its spaces therefore has no Hilbert
   adjoint to take, and :class:`~orpheus.numerics.operator.AdjointOperator`
   refuses it at construction with
   :class:`~orpheus.numerics.operator.MissingAdjoint`, rather than
   silently applying the Euclidean transpose under the adjoint's name.
   The only exemption is the metric-free stratum (a constant multiplier
   :math:`cI` commutes with every metric,
   :class:`~orpheus.numerics.operator.ConstantMultiplier`), whose missing
   end becomes an :class:`~orpheus.numerics.operator.IdentityOperator`
   leg. A field multiplier commutes only with a diagonal metric and is not
   exempt (`#625 <https://github.com/deOliveira-R/ORPHEUS/issues/625>`_). This is why
   the adjoint requires a **bound operator** (:ref:`bound-operator`).

**Where it is realised.** The two legs are operators built from a space,
``RieszLowerOperator(space)`` and ``RieszRaiseOperator(space)``
(:class:`~orpheus.numerics.operator.RieszLowerOperator`,
:class:`~orpheus.numerics.operator.RieszRaiseOperator`), and
:meth:`LinearOperator.dual <orpheus.numerics.operator.LinearOperator.dual>`
mints the dual arrow. :class:`~orpheus.numerics.operator.AdjointOperator`,
reached as ``A.H``, builds the three factors from the operator's two
ends, ``RieszRaiseOperator(A.domain) ∘ A.dual() ∘
RieszLowerOperator(A.codomain)``, and composes them. The space carries no
leg property: the ``FunctionSpace.riesz_lower`` and
``FunctionSpace.riesz_raise`` properties, which no production code read
(`[M]` 2026-09-22, an AST census of ``orpheus/``), retired in S1 of the
numerics untangling (2026-10-19), when the space lost its operator
imports. The landing of the legs (CS4c steps 1–2, commits ``68a9c9f3``
and ``733d96f3``) is recorded in the S\ :sub:`N` development history
(:ref:`sn-development-history`).

.. _spaces-dual-space:

The dual space: the covectors, carrying the inverse form
--------------------------------------------------------

:math:`V^{*}` is the space of linear functionals on :math:`V`. In
coordinates a functional is an array of :math:`V`'s shape, so the dual has
the primal's coordinates; what makes it a different space is its metric.
The Riesz map :math:`\flat = G` is an isometry from :math:`V` onto its
dual exactly when the dual pairs functionals through the inverse form,

.. math::

   \langle f, g\rangle_{V^{*}} \;=\; f^{\mathsf T}\,G^{+}\,g ,
   \qquad
   \langle \flat x, \flat y\rangle_{V^{*}}
   \;=\; x^{\mathsf T} G\,G^{+} G\, y \;=\; \langle x, y\rangle_V ,

with :math:`G^{+}` the Moore–Penrose pseudo-inverse
(:ref:`spaces-riesz-legs`; :math:`G G^{+} G = G` is the first Penrose
condition, so the isometry holds on a singular metric too). ORPHEUS spells
this as two objects, one per side of the metric/space split.

**The metric of the dual.** :meth:`HilbertMetric.dual
<orpheus.numerics.metric.HilbertMetric.dual>` returns a
:class:`~orpheus.numerics.metric.DualMetric`, the primal metric with its
two faces swapped: its ``apply`` is the primal's ``apply_inverse`` and its
``apply_inverse`` the primal's ``apply``. Every realization already stores
both faces (a :class:`~orpheus.numerics.metric.DiagonalMetric` its weights
and their masked reciprocal, a :class:`~orpheus.numerics.metric.DenseMetric`
:math:`G` and the stored ``pinv(G, hermitian=True)``), so the dual
computes no new inverse and its arithmetic is the primal's, bit for bit.
That is why the dual is one wrapper over every realization rather than a
``DiagonalMetric(1/w)``: the reciprocal re-inverted differs from the
masked reciprocal in the last bit, and :math:`1/(1/w) \ne w` in floating
point. The dual is cached per metric object and is an involution by object
identity, ``m.dual().dual() is m``, because :math:`(G^{+})^{+} = G` for a
symmetric positive semi-definite form. A
:class:`~orpheus.numerics.metric.FactoredMetric` duals each block,
:math:`(G_1\otimes G_2)^{+} = G_1^{+}\otimes G_2^{+}`, so the dual of a
product metric stays a positioned product whose entries are the blocks'
:class:`~orpheus.numerics.metric.DualMetric`\ s. Gates:
``tests/gates/numerics/test_dual_metric.py`` M1, M1b (the face swap and
its positioned block form, bitwise), M2 (the involution), M3 (a dense dual
computes no inverse: a spy on ``pinv`` counts zero calls) and M4 (the
pairing keeps its one spelling, :ref:`spaces-metric-one-spelling`).

**The dual of a leaf space.** :meth:`FunctionSpace.dual
<orpheus.numerics.space.FunctionSpace.dual>` of a leaf space is a
:class:`~orpheus.numerics.space.DualSpace` whose content IS its primal.
It copies no axes, weights or metric object: its resolved metric is the
primal's resolved metric's ``.dual()``, read per call through the primal.
Its constructor refuses a shape that is not the primal's and any metric
source of its own (``ValueError``), since such a dual would disagree with
its primal. `[M]` 2026-10-21, the archivist's probe on the S6 worktree: on
the dense fixture of ``test_dense_metric.py`` (the exact-binary
:math:`G = [[2, 0.5, 0], [0.5, 1, 0.25], [0, 0.25, 4]]`) and
:math:`x = \mathbf 1`, the primal pairs :math:`\langle x, x\rangle_V = 8.5`
and the dual :math:`\langle x, x\rangle_{V^{*}} = 1.3`, which is
:math:`x^{\mathsf T}G^{+}x` to the last two bits. Gates: S1 (the dual
reads the inverse metric bitwise), S2, S2b and S2c (the dual pairing is
the pairing of the raises, against an exact rational inverse on the dense
fixture), S3 and R4 (the primal is the content; the refusals).

**V\* is not V** (ruling 75 of ``.claude/plans/numerics_untangling.md``).
A dual space equals another dual space iff their primals are equal, and it
never equals a primal space, in either order of the comparison:
``V.dual() == V`` is ``False`` and ``V == V.dual()`` is ``False``, and
``hash(V.dual())`` is ``hash(("dual", V))``. Covariance is therefore part
of a space's identity, and the composition guard refuses a covector where
a vector is expected: ``R @ R.dual()`` and a sum joining :math:`V` and
:math:`V^{*}` raise
:class:`~orpheus.numerics.operator.IncompatibleOperatorComposition`
(gates S4, S5, G1, G2).

**One dual per space.** The dual is a function of the space alone, so
``V.dual()`` is memoised on the space: ``V.dual() is V.dual()``, and
``V.dual().dual() is V`` (the involution :math:`V^{**} = V`, by object
identity). ``DualSpace.of(V)`` is the same object as ``V.dual()``. The
memo is what bounds the dual arrows' cost: an adjoint solve rebuilds its
dual arrows on every matvec (``OperatorProduct.dual`` rebuilds a word's
dual word each call, `#632 <https://github.com/deOliveira-R/ORPHEUS/issues/632>`_),
and each arrow asks for its ends' duals. `[M]` 2026-10-21, the S6a
verification spec's probe p1 (``solve_sn_adjoint`` on the 3-region 2-group
slab, GL-8, 40 cells): without the memo one solve minted 5 645
``DualSpace``\ s, and 10 965 at tight tolerances, every one from a dual
arrow's ``domain`` or ``codomain``; with it, 2 (the duals of the interior
and the trace block). The gate is ``test_dual_metric_sn.py`` K2: the mint
count is the same at two tolerance settings that run different numbers of
outer iterations, and is bounded by the number of distinct primals.

**The dual of a composite is the composite of the duals.** Duality is
additive on a direct sum and multiplicative on a tensor product, so a
composite space's dual is built from its members' duals rather than
wrapped as a :class:`~orpheus.numerics.space.DualSpace` (which would have
no member metrics to read, and would pair Euclidean):

.. math::

   (V_{\rm bulk} \oplus V_{\rm trace})^{*} = V_{\rm bulk}^{*} \oplus V_{\rm trace}^{*},
   \qquad
   (V \otimes W)^{*} = V^{*} \otimes W^{*} .

:class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace` duals
its interior and trace blocks;
:class:`~orpheus.numerics.coupled_system.CoupledSpace` duals its systems
and swaps its two zero factories, since this space's cotangent zero is
the zero of its dual; :class:`~orpheus.numerics.space.TensorProductSpace`
duals its factors. Each override is the one hook ``_make_dual``, and the
memo and the involution are the base's. The composite's metric verbs then
dispatch to the member duals' inverse forms. `[M]` 2026-10-21, the archivist's
probe: for :math:`V` with weights :math:`[0.5, 2, 4]`, :math:`W` with
weights :math:`[3, 5]` and :math:`x = \mathbf 1`,
``(V * W).dual() == V.dual() * W.dual()`` and its pairing reads
:math:`1.4\overline{6} = \sum_{ij} w_{V,i}^{-1} w_{W,j}^{-1}`. Gates:
``test_dual_metric_sn.py`` C1 (a composite dual reads the member inverse
metrics), C2 (it is the biproduct of the member duals), R2b
(``DualSpace.of`` on the coupled carrier and its full-field member is the
composite dual).

The tensor product resolves each factor's metric once, through the factor's
own resolution (an axis-built factor's per-axis entries, a metric object,
a dense-slot leaf's weights, a dual's inverse form). A dual factor carries
no ``axes``, ``metric`` or weights of its own, so a builder that read
those raw fields would see it as Euclidean: `[M]` 2026-10-21, the S6a
elegance review's probe P3, on ``V.dual() * W`` the pairing read
``385.0`` where :math:`G^{+}\otimes w` gives ``163.25``, before the
builder was changed to read each factor's resolved metric (gate
``test_dual_metric.py`` R1).

**The dagger–dual square commutes, as a law.** With every dual space
carrying its inverse form, the dual arrow
:math:`A^{\mathsf T} : W^{*} \to V^{*}` is a bound operator whose ends
give Riesz legs, and its Hilbert adjoint is the base composition
:eq:`spaces-adjoint-riesz-composition` over those ends,
:math:`(A^{\mathsf T})^{*} = \sharp_{W^{*}} A \,\flat_{V^{*}} = G_W A G_V^{+}`.
The transpose of the primal adjoint is the same matrix,
:math:`(A^{*})^{\mathsf T} = (G_V^{+} A^{\mathsf T} G_W)^{\mathsf T} =
G_W A G_V^{+}`, since both metrics are symmetric. The two are computed by
different routes (the base sandwich on the dual arrow, and the primal
adjoint's ``.dual()``) and agree bitwise, because a dual metric is its
primal with the faces swapped and the two routes run the same three
delegations in the same order. Gates: ``test_dual_metric.py`` D2 and
``test_dual_metric_sn.py`` C3 (bitwise on the S\ :sub:`N` carrier), R5
(reciprocity of the dual of a sum in the dual end metrics), and
``test_riesz_legs.py``'s dagger–dual row. D1 asserts that neither the dual
arrow nor the adjoint arrow defines a route of its own (``adjoint`` is not
in ``vars(_DualOperator)``, ``dual`` is not in ``vars(AdjointOperator)``).

.. warning::

   **A dual arrow does not swap its carrier's role.** For
   :math:`S` : flux :math:`\to` source, ``S.dual()`` is typed
   :math:`W^{*} \to V^{*}` but its action still accepts a flux-classed
   field and returns a source-classed one, so ``S.dual().H`` cannot take
   ``S.dual()``'s own output. Covariance is typed on spaces, not yet on the
   carriers (`#633 <https://github.com/deOliveira-R/ORPHEUS/issues/633>`_).
   A test of a dual arrow samples the class its ``apply`` admits.

.. _spaces-quotient:

A degenerate metric is a quotient: rank and the live projector
--------------------------------------------------------------

A metric with a kernel (the S\ :sub:`N` trace's tangential ordinates,
where :math:`|\Omega\cdot\hat n| = 0`; dead moment-head slots; a
rank-deficient Parseval head) is an inner product on the quotient
:math:`V / \ker G`, not on :math:`V`. ORPHEUS keeps the coordinates and
lets the metric carry its kernel (ruling 74): the tangential values stay
where the sweep and the reflections read them, and no shape changes at any
trace, full-field or coupled carrier. `[M]` the S6 census
(``scratch/numerics_untangling/packages/S6_census/report.md``): 3 380 of
62 432 spaces built over the numerics and ``sn/operators`` suites have a
degenerate metric, in those three families.

Two readings make the quotient explicit:

- :meth:`HilbertMetric.rank(shape) <orpheus.numerics.metric.HilbertMetric.rank>`
  is the dimension of :math:`V/\ker G` on elements of ``shape``. A diagonal
  metric counts its nonzero broadcast weights; a dense metric uses
  ``matrix_rank`` at the same cutoff its stored pseudo-inverse uses, so the
  rank counts exactly the directions :math:`G^{+}` inverts; a factored
  metric multiplies its block ranks (a Euclidean block is full rank); a
  dual metric has its primal's rank, since :math:`G^{+}` and :math:`G`
  share their range. :attr:`FunctionSpace.rank
  <orpheus.numerics.space.FunctionSpace.rank>` reads its metric (the
  element size when the space is Euclidean), and a direct sum or a coupled
  space sums its members' ranks.
- :meth:`FunctionSpace.project_live <orpheus.numerics.space.FunctionSpace.project_live>`
  is :math:`P = G^{+}G` (:eq:`spaces-riesz-round-trip`), the orthogonal
  projector onto :math:`\operatorname{range} G = (\ker G)^{\perp}`. It is
  the identity on a non-degenerate space and zeroes the kernel
  coordinates of a diagonal degenerate one: for weights :math:`[1, 0, 2]`
  the rank is 2 and :math:`P\mathbf 1 = [1, 0, 1]` (`[M]` 2026-10-21, the
  archivist's probe). Two elements are the same element of the quotient
  iff their live projections agree, which is how a law that holds modulo
  the kernel is stated and gated.

**The adjoint is the quotient adjoint.** By the kernel condition
(:ref:`the kernel condition <spaces-adjoint-kernel-condition>`), reciprocity holds on every
:math:`x` iff the operator maps :math:`\ker G_V` into :math:`\ker G_W`.
`[M]` 2026-10-21, ruling 80 (``scratch/numerics_untangling/packages/S6bc_census/``,
probes ``p2_kernel.py`` and ``p3_leaves.py``): on the cylinder
``folded_product(2, 6)``, the sphere GL-7 and the slab GL-7, the
kernel-to-live block of every production S\ :sub:`N` operator is exactly
``0.0``, on 87 of 87 rows (the reflective and vacuum boundary, streaming,
scattering, :math:`N_{2n}`, fission, the radial-characteristic blocks and
the sweep); a control boundary operator that feeds a tangential slot into
a live one reads ``0.70``. Production ``.H`` equals
:math:`P_V A^{*} P_W` to :math:`2.1\times10^{-16}`. The adjoint laws are
therefore stated modulo the kernel:

- ``A.H.H is A`` is a choice of representative on :math:`V/\ker G`. On the
  cylinder sweep the double adjoint computed out differs from :math:`A` by
  0.519, all of it on the kernel-to-kernel block (the sweep passes the
  tangential slots through as the identity), so the two agree after
  projection.
- Reciprocity, :math:`\flat^{*} = \sharp` and ``M @ M.inverse()`` are
  gated in the G-seminorm, never component by component on the kernel:
  :math:`M M^{-1} q = P q`, and :math:`M M^{-1} e_k = 0` on every dead
  slot :math:`k`.

Gates: ``tests/gates/numerics/test_quotient_metric.py`` (B1 the metric
rank against hand counts, B2 a Parseval head's rank against the sampled
harmonic rank, B3 and B3b the projector, B4 the round trip, B5 the
retraction closed form against the sandwich modulo the kernel) and
``tests/gates/sn/operators/test_quotient_adjoint_sn.py`` (T1 the trace
rank is the non-tangential count, K80 every degenerate-ended leaf maps
the kernel into the kernel, A1 to A3 the sweep's laws on the quotient).

.. _spaces-metric-propagation:

Composition: a dropped metric is a VALUE bug wearing a representation costume
-------------------------------------------------------------------------------

A new representation is invisible to every consumer that propagates the
old one by copying an array. Three such sites exist, and all three were
measured *silently wrong* on a dense-metric space before P7 taught them
— each failing in the flattering direction (reverting to Euclidean, or
to a plausible-looking value, with no error and no warning):

.. list-table:: The three propagation sites
   :header-rows: 1
   :widths: 26 37 37

   * - Site
     - Pre-P7 behaviour on a dense-metric operand
     - What it does now
   * - :meth:`DualSpace.of
       <orpheus.numerics.space.DualSpace.of>`
     - copied ``inner_product_weights`` only, so the dual of a
       dense-metric space read the **plain Euclidean** pairing.
       `[M]` on the ``test_dense_metric`` fixture the dual read
       :math:`4.5` where the primal reads :math:`23.25`
     - reads the primal's metric object through the primal, dualised:
       the dual carries the inverse form :math:`G^{+}`
       (:ref:`spaces-dual-space`), so it reads
       :math:`x^{\mathsf T}G^{+}x`, never a copy and never Euclidean
   * - :meth:`TensorProductSpace.from_factors
       <orpheus.numerics.space.TensorProductSpace.from_factors>`
       (the legacy ``*``)
     - densified the factor measures into one outer-product weight
       array — which has no slot for a matrix factor, so the dense
       factor was **dropped** and its block went Euclidean
     - **since 2026-09-07 (CS4c step 6 item 6.2a) the factored arm is
       the ONLY non-axis arm**: a lazy
       :class:`~orpheus.numerics.metric.FactoredMetric` positioned per
       **axis** of an axis-built factor and per dense-slot leaf factor,
       Kronecker never materialized, and ``inner_product_weights`` never
       populated by a product at all.
       ⛔ This cell read *"grows a dense-factor arm: a lazy*
       ``FactoredMetric``\ *, one positioned entry per factor"* until
       then. That was true of P7 and it described an arm that was
       **inert**: `[M]` the step-6 activation census counted **0**
       entries across 11 SN runs (~450 product mints), because no factor
       reaching ``*`` on an SN path carried a metric OBJECT — the only
       production installs are the frame's
       ``DenseMetric.inverse_of(discrete_gram)`` and the two ψ½ ray
       spaces, none of which enters a product. What 6.2a changed is the
       arm's REACHABILITY, not its existence: the dense arm it shared the
       dispatch with is gone, so every axes-less product now takes it
   * - ``FrameBase.gram`` (retired 2026-09-08 → :meth:`gram_inverse <orpheus.numerics.frame_faces.CrossGramInverse.of>`)
     - ``replace(test_space, inner_product_weights=diagonal)`` on a
       dressed test space handed the **row-sum probe** a space whose
       ``apply_inverse_metric`` ignores that slot and applies its own
       matrix
     - **stripped** the metric object while installing the probe
       diagonal — the cross-Gram machinery must never inherit the
       trial-side Parseval dressing. ⭐ Item 6.2c-ii went further and made
       the state unspellable: the normalisation is now the typed arrow
       :class:`~orpheus.numerics.frame_faces.CrossGramInverse`, whose action
       reads the probe diagonal and NO space's metric
       (:ref:`frame-gram-inverse-arrow`)

The third is the sharpest, because it is not a missing feature but a
live projection error on a production path
(:func:`project <orpheus.numerics.frame_faces.project>`
is what homogenisation and condensation call). `[M]` on the overlap
frame — Gram
``[[1.25, 0.25], [0.25, 1.25]]``,
condition number :math:`1.50` — the pre-P7 spelling computes
:math:`G\,(Mf)` instead of :math:`\operatorname{diag}^{+}(Mf)` and
returns ``[7.0, 11.0]`` where the true projection is
:math:`[8/3,\,16/3]`: **162.5 %** wrong on the first component,
**106.3 %** on the second. Post-P7 that spelling is not merely wrong,
it is **unrepresentable** — the exclusivity guard refuses the
two-source space at construction, `[M]` with
``ValueError: space 'L2[coarse_cells(spatial_R1)]' carries BOTH dense
inner_product_weights and a metric object`` (the space name gained its
manifold at #429 tracker 2.1; the guard is unchanged). The illegal state that
produced a 162 % value error is now a state that cannot be built.

.. note::

   **The tensor product's laziness is the same discipline the axis path
   already follows.** :class:`~orpheus.numerics.metric.FactoredMetric`
   stores ``(block_shape, factor)`` pairs and applies them in sequence,
   each to its own index block; the Kronecker product
   :math:`G_1\otimes G_2\otimes\cdots` is never formed. The pairing is
   nonetheless exact against a densified reference: `[M]` a
   :math:`(3,)` dense-metric factor times a :math:`(2,)` weighted
   factor reads ``141.5``, bit-equal to
   :math:`y^{\mathsf T}\bigl(G\otimes\operatorname{diag}(w)\bigr)x`
   built with :func:`numpy.kron`, where dropping the dense factor would
   have read ``45.0``. On separable probes the pairing factorizes
   exactly, which is what makes a hand-derived pin possible:
   :math:`23.25 \times 2.5 = 58.125` in exact binary fractions.

.. warning::

   **This does NOT mean an axis-built product can carry a matrix.**
   :eq:`spaces-axis-product` says the metric of an ``of_axes`` product
   is :math:`\bigotimes_a G_a` with
   :math:`G_a = \operatorname{diag}(w_a)` or :math:`I`, and that
   remains exactly true — see :ref:`spaces-metric-not-on-the-axis`. The
   :class:`~orpheus.numerics.metric.FactoredMetric` arm belongs to the
   **legacy** ``*`` composition, which takes whole *spaces* as factors
   and therefore can meet one that carries a metric object. ⛔ This
   paragraph closed *"The two paths were not unified by P7 and CS2 still
   retires the second"* until **2026-09-07**. The retirement it promised
   landed as **CS4c step 6 item 6.2a** — and what retired is the
   *densifier*, not ``*``: the product's dense weights array and the
   bridge that fed it are gone, so both paths now apply their measure
   factor-by-factor. They are still not one path, and the residue is
   narrow and named: an axes-less factor makes the product axes-less, so
   ``*`` cannot answer an axis query. `[M]` 2026-09-07 the only production
   occupants are the four harmonic/moment mints; item **6.2c** axis-ifies
   the angular head and closes it.

.. _spaces-metric-dressing-evidence:

The dense dressing, and the only evidence that can adjudicate a metric
------------------------------------------------------------------------

The founding consumer installs it:
:func:`basis_space <orpheus.numerics.frame_faces.basis_space>` on a ``DENSE`` verdict
now returns the basis's space dressed with
``DenseMetric.inverse_of(discrete_gram)`` — the pseudo-inverse at the
pinned cutoff, with the exact symmetrized Gram kept as the inverse face
(:math:`(G^{+})^{+} = G` for a symmetric PSD form, so this is exact and
strictly better conditioned than a second pseudo-inversion) — and
**strips** the basis's continuum weights, because the dressing
*replaces* the metric on that arm exactly as the diagonal arm
overwrites it. Full account, including the pre-repair slab Gram's table:
:ref:`frame-parseval-dense-arm`.

⛔ **Reciprocity cannot adjudicate a metric, and it never will.** The
Hilbert-adjoint identity :math:`A^{\dagger} \equiv G^{-1}A^{\mathsf T}G`
holds to :math:`1.4\times10^{-16}` for **every** invertible :math:`G` —
Euclidean, random, :math:`(V\!\cdot\! w)^3` — because ``.H`` is *built
from* the stored metric and a sandwich always reproduces the pairing it
was assembled from (ERR-039's first shield, at this layer). A
reciprocity gate therefore proves *loadedness* and carries **zero**
information about *which* metric is installed. The instrument that can
fail is one that compares the metric against a quantity defined without
it — Parseval, i.e. the field's own norm — plus hand-derived literals
with no solver, frame or quadrature in the chain.

.. list-table:: The wrong-metric discriminator — `[M]` 2026-10-18, the draw-free range of the Parseval ratio over the band-limited subspace of ``frame_faces.angular_frame(folded_product(4, 4), 2)``
   :header-rows: 1
   :widths: 40 24 36

   * - Metric on the coefficient space
     - Parseval ratio range
     - Verdict
   * - the basis's **continuum** Gram (pre-F-0, undressed)
     - :math:`[6.317,\ 157.9]`
     - the ERR-039 defect
   * - the best **diagonal** candidate, :math:`1/\operatorname{diag}(G)`
     - :math:`[1.000,\ 2.000]`
     - a diagonal metric is *provably insufficient* here, not merely
       unavailable
   * - the **matrix** pseudo-inverse :math:`G^{+}` (shipped)
     - :math:`[1.000,\ 1.000]`
     - the theorem, :eq:`spaces-pseudo-inverse-parseval`

The middle row is the whole justification for the phase, and it is the
row a diagonal-only design could never have produced. The range is the
generalized-eigenvalue range of the Rayleigh quotient
:math:`(Gc)^{\mathsf T} M (Gc) / c^{\mathsf T} G c`, so it is a property of
the frame and cannot move with a random draw
(``scratch/metric_docstring/probe_theory.py``). ⚠ It is **frame-dependent**:
`[M]` the diagonal candidate ranges over :math:`[0.707,\ 1.842]` on
``frame_faces.angular_frame(product(4, 4), 2)`` and :math:`[0.889,\ 1.111]` on
``frame_faces.angular_frame(level_symmetric(4), 3)``, so a single draw there can land
near 1. The committed discriminator,
``test_no_diagonal_metric_can_satisfy_parseval_on_a_dense_frame``, reads
``frame_faces.angular_frame(folded_product(2, 4), 3)``, where the diagonal range is
:math:`[1.000,\ 3.000]` and never drops below 1. (Until the ERR-080
repair this table read one seeded :math:`\psi` on the slab
``frame_faces.angular_frame(gauss_legendre(8), 2)``: 25.53, 1.806 and
0.999999999999999.)

**The dressing moves a production adjoint, and nothing observed it.**
This is not plumbing. The scattering operator builds
``frame_faces.angular_frame(quadrature, scattering_order)`` in production, so
dressing a ``DENSE`` frame changes the spaces every ``.H`` through the
adjoint sandwich reads. `[M]` comparing the dressed and undressed
analysis adjoints **as operators** (column by column on the coefficient
basis, so the number carries no random draw):

.. list-table:: Movement of ``FrameAnalysis.of(frame).H`` under the dressing
   :header-rows: 1
   :widths: 34 22 22 22

   * - Frame
     - :math:`\max|\Delta|`
     - relative (max-norm)
     - relative (Frobenius)
   * - ``product(4,4)``, :math:`L=2`
     - :math:`12.49`
     - :math:`0.994`
     - :math:`0.985`
   * - ``gauss_legendre(8)``, :math:`L=2`
     - :math:`12.39`
     - :math:`0.986`
     - :math:`0.980`
   * - ``level_symmetric(4)``, :math:`L=3`
     - :math:`12.49`
     - :math:`0.994`
     - :math:`0.980`

The two operators are not a small correction apart; they are
essentially unrelated, which is the correct reading of a repair whose
"before" state the tree's own docstring called *"the stored-metric
sandwich, NOT the physical Hilbert adjoint"*. ⚠ Do not freeze a
single-vector figure for this: a random-probe relative movement on
``product(4,4)`` is `[M]` **0.879–0.986** over 200 seeds, so a quoted
point value inside that band is one draw's reading, not a property of
the frame.

.. _spaces-metric-frame-square:

What does NOT ride along: the frame square's scalar closure
-------------------------------------------------------------

Installing the right metric does **not** make every identity that
mentions a metric come true, and the distinction is worth a subsection
because it is exactly where a reader will over-read the repair.

The spherical-harmonic frame square closes on one scalar,
:math:`M^{*} = R/W` and :math:`R^{*} = W M`
(:eq:`frame-square-closure-sh`), and that identity needs the per-degree
Gram :math:`G_\ell` to be **a single number for each** :math:`\ell`.
Writing :math:`d` for the addition-theorem factors
(:math:`d_\ell = 2\ell+1`) and :math:`Y` for the synthesis table, the
closure is the operator statement

.. math::

   Y\Bigl(G^{+} - \tfrac{1}{W}\operatorname{diag}(d)\Bigr) \;=\; 0 ,

i.e. the metric and the reconstruction weights need agree only **modulo
the kernel of** :math:`Y`. That makes the relation decidable and it
splits the shipped frames three ways rather than two:

.. list-table:: `[M]` 2026-08-30 — the closure residual is not implied by the metric being right
   :header-rows: 1
   :widths: 30 14 22 34

   * - Frame
     - Verdict
     - rel :math:`\|M^{*}y - Ry/W\|`
     - Why
   * - the six ``DIAGONAL`` sphere frames — LS\ :sub:`4` at
       :math:`L{=}1` and :math:`L{=}2`, and LS\ :sub:`8`, product
       :math:`8\times8`, folded :math:`8\times8` and Lebedev-13 at
       :math:`L{=}2`
     - ``DIAGONAL``
     - :math:`\le 1\times10^{-15}`
     - each live :math:`\ell` block is one constant; per-degree
       diagonal spread :math:`\le 4\times10^{-15}`
   * - ⛔ ``gauss_legendre(8)``, :math:`L=2` — **this row is HISTORY**
     - was ``DENSE``
     - was :math:`0.30`–:math:`10.2` (200 seeds)
     - the live :math:`\ell{=}2` diagonal was
       :math:`[0.4,\,0.8,\,0.8]` — no :math:`G_\ell` existed, at
       **any** metric. ⭐ Both :math:`0.8`\ s were ERR-080's fabricated
       :math:`m \ne 0` slots; since #429's fused commit (2026-09-02) a
       1-D rule binds a FLAT Legendre head and this frame is
       ``DIAGONAL`` with :math:`G_\ell = 2/(2\ell+1)` — one number per
       degree — so it moves into the **first** row of this table.
       `[M]` 2026-09-02, 200 seeds: rel :math:`\le 5.1\times10^{-16}`,
       and the same at ``gauss_legendre(4)``, ``(16)`` and at
       :math:`L = 3`
   * - ``product(4,4)`` :math:`L{=}2`; ``level_symmetric(4)``
       :math:`L{=}3`; ``folded_product(2,4)`` :math:`L{=}3`
     - ``DENSE``
     - `[M]` 2026-09-02, 200 seeds:
       :math:`7.6\times10^{-4}`–:math:`0.47`,
       :math:`4.3\times10^{-2}`–:math:`0.15`,
       :math:`0.26`–:math:`1.65`
     - same cause, milder: only the top degree's block is non-constant.
       (The third is new here — it takes over the slab's role as the
       loudest breaker, on a SHIPPED rule and a basis the repair does
       not touch)
   * - ``folded_product(4,6)``, :math:`L=3`; ``gauss_legendre(2)``,
       :math:`L=2`
     - ``DENSE``
     - :math:`\le 3\times10^{-15}` / `[M]`
       :math:`\le 1.6\times10^{-15}` (200 seeds each)
     - their coupling is a pure **rank deficiency**, so the
       disagreement lives entirely in :math:`\ker Y`. ⭐ The second is
       a 2026-09-02 addition and it is the stronger exhibit, because it
       is a **Legendre** frame and its rank deficiency has a closed
       form: :math:`P_2` vanishes identically at ``GL_2``'s two nodes,
       which ARE its roots (the dead-slot theorem —
       :math:`G_{22} \sim 10^{-33}`). Two families, one mechanism

The last row is the one that makes the statement a theorem rather than a
correlation, and it was found by measuring rather than by reasoning
(and it now has a second, independently-caused member — see the
``gauss_legendre(2)`` note in that row).
That frame's only live off-diagonal couples two :math:`\ell = 3` slots
whose :math:`2\times2` block is
``[[0.6732, 0.8691], [0.8691, 1.1220]]``
— `[M]` determinant :math:`-8.7\times10^{-17}`, **rank 1**: the two
harmonics are linearly *dependent* on that folded node set. So
:math:`G^{+}` and :math:`\operatorname{diag}(d)/W` differ by an
:math:`O(1)` amount that :math:`Y` annihilates: `[M]`
:math:`\|Y D\|_\infty = 4.4\times10^{-16}` with
:math:`\|D\|_\infty = 0.557`, against
:math:`\|Y D\|_\infty = 6.30` on the slab.

⟹ **the** ``DIAGONAL`` **verdict is sufficient for the scalar
closure, and** ``DENSE`` **does not decide it.** That is why the
frame-square gate keeps only the
diagonal parameters — not because a dense frame always breaks the
closure (one shipped frame does not), but because the verdict does not
*imply* it, and a gate must not assert what its population cannot
guarantee. The discriminator is Gram structure, never geometry: a
*sphere* rule (``product(4,4)``) is among the frames that break it.

.. _spaces-metric-not-on-the-axis:

Where the metric may NOT live, and what stays refused
-------------------------------------------------------

Three rulings bound the phase. Each is a decision with a reason, not a
scheduling accident.

**1. The axis keeps a MEASURE; it does not grow a form.**
:attr:`Axis.weights <orpheus.numerics.axis.Axis.weights>` stays a 1-D
weight vector, and :eq:`spaces-axis-product` stays literally true. The
distinction is not notation: *a measure is diagonal by nature* — it
assigns a number to each atom of an index set — whereas a Gram is a
**bilinear form** on the space the atoms span, and the two coincide only
when the basis is orthogonal on that measure. The generator induces, the
space holds: a frame that measures its own Gram dresses the *space* it
mints, and no axis is asked to carry a matrix it has no way to produce.

**2. Identity stays metric-blind** — ⛔ **OVERTURNED 2026-09-08, twice
over; the paragraph is preserved because its REASON survived while its
claim did not.** As written (P7, 2026-08-30): *"Space identity is*
:math:`(\text{name}, \text{shape})`\ *, so a dressed and an undressed
space of the same name compare equal and hash equal, and the frame's*
``basis_space == basis.space`` *invariant survives the dressing
untouched. … (The chartered doctrine that* metric differences imply
space differences *is unaffected — it flows through the axis-derived
name,* :ref:`spaces-identity-bridge`\ *, and an axis-built space cannot
carry a metric object at all.)"*

Two later steps falsified both of its clauses, in opposite directions:

* the **identity flip** (CS4c step 6, 2026-09-07) made an axis-built
  space compare by its ``axes`` tuple DIRECTLY, so the doctrine stopped
  flowing through the derived name and started being the definition; and
* item **6.2c-ii** (2026-09-08) made the moment heads axis-built, so
  ``basis_space`` and the basis's own space (``basis.space`` until S7,
  ``FunctionSpace.of_basis(basis)`` since) differ in their head axis's
  measure and are therefore **different spaces**. `[M]` re-measured on the same
  33 shipped (rule, :math:`L`) rows: ``basis_space(frame) ==
  frame.basis.space`` is **0 of 33** where it was 33 of 33.

Item **6.2c-i** (2026-09-08) also repealed the parenthetical's last
clause: an axis-built space CAN carry a metric object, as the positioned
overlay of forms merged into its axes-derived metric
(:ref:`spaces-metric-three-sources`) — which is how a dense Gram lives on
an axis-built head at all.

What survives untouched is the **reason** the ``metric`` field is
declared ``compare=False``, and it is worth separating from the claim it
was attached to: it is **structurally mandatory**, not a statement about
identity. An ndarray inside a dataclass-generated ``__eq__`` makes the
comparison return an array and ``hash()`` raise, so a *subclass's*
generated ``__eq__`` must never reach the field — exactly the reason
``inner_product_weights`` and ``axes`` carry the same flag. The manual
``__eq__`` reads ``axes`` regardless; ``compare=False`` governs what the
dataclass machinery generates, not what identity IS.

⚠ **The consequence for any GATE over these two spellings, and it
INVERTED.** *(2026-09-02:)* a ``==`` assertion was structurally blind to
which spelling a producer bound, so it could not adjudicate between them;
`[M]` over 33 shipped (rule, :math:`L`) rows the two were ``(name,
shape)``-equal on **33 of 33** and metric-different on **33 of 33**, and
the gate that pinned #429 tracker 2.5's choice therefore had to assert
the metric ARRAY with the other spelling as its negative control
(``vv-principles`` #19). Since 6.2c-ii the metric IS the identity on
these spaces, so ``==`` is the discriminating instrument — which is what
lets the hub/frame agreement gate assert ``==`` and mean it
(``tests/gates/sn/mesh/test_hub_and_frame_agree_on_the_moment_space.py``). ⭐
The transferable half: *whether a* ``==`` *gate is blind or
discriminating is a property of the identity relation, not of the gate* —
so a re-typing that moves a field into the identity silently upgrades
every equality assertion in the corpus, and a docstring calling one of
them "metric-blind" goes stale without its file being touched. Full
account: :ref:`frame-moment-space-single-home`.

**3. Expressible is not known — the curvilinear moment mass stays
refused.** The multi-moment cell mass on a curved chart was blocked on
**two** independent things, deliberately named because they need
opposite repairs: the *machinery* (a non-Hadamard space metric,
`GitHub #409 <https://github.com/deOliveira-R/ORPHEUS/issues/409>`_) and
the *value* (which :math:`G` is right is pinned by physical functionals
and needs the cell-solve consumer of
`#158 <https://github.com/deOliveira-R/ORPHEUS/issues/158>`_). P7 lands
the machinery only. The refusal in
``DiscretizationSchemeBase._assert_moment_mass_is_expressible`` therefore
re-derives to name **#158 alone**, and nothing installs a guess: a
metric that can now be *spelled* is not a metric anyone has *chosen*.

.. note::

   **Why the axis ruling and the curvilinear refusal are the same
   measurement read twice.** The two consumers need different strengths,
   and that asymmetry is the strongest argument for an
   operator-valued metric over a "matrix weights slot on ``Axis``":

   - a dense angular frame (today ``frame_faces.angular_frame(folded_product(4, 4), 2)``;
     the slab frame before the ERR-080 repair) needs **one
     cell-independent** dense :math:`K\times K` coefficient Gram — a
     per-axis matrix could carry that;
   - the curvilinear mass is
     :math:`G_{(i,k),(i,j)} = V_i\,M_{kj}(r_i)` — dense in the moment
     index **and** varying along the spatial index, i.e.
     **non-separable** across (spatial :math:`\times` moment). No
     per-axis object of any shape, vector or matrix, can carry it.

   A design that closed the first consumer by widening ``Axis`` would
   have closed exactly one of the two and made the second harder, by
   putting the metric on the factor that cannot express it. The
   space-held object closes the first and leaves the second's *machinery*
   available for whatever chooses its value.


.. _spaces-moment-head:

The angular HEAD of a moment space — two families, one door
-----------------------------------------------------------

The metric sections above are about *what a factor measures with*. This
one is about a question that had only one answer until #429 (2026-09-02):
**what SHAPE does the leading factor of a moment space have, and who is
allowed to know?**

A moment field's space is ``<head> ⊗ cells``: the head is the
coefficient space of the basis the quadrature's frame binds, one axis
whose layout is a degree layout. Two families ship. The real spherical
harmonics' head is the rectangular :math:`(L+1, 2L+1)` table, a
:class:`~orpheus.numerics.coefficient_layout.HarmonicLayout`, with the
addition-theorem-shifted :math:`[\ell + m]` column and zero padding
outside :math:`|m| \le \ell`. The Legendre basis
:math:`\{P_\ell(\mu)\}_{\ell \le L}` on the orbit space
:math:`S^2/O(2)_a` has a **flat** head, :math:`(L+1,)`, a
:class:`~orpheus.numerics.coefficient_layout.LegendreLayout`: one
coefficient per degree, because the trivial isotypic component of the
:math:`O(2)_a` action — equivalently of its rotation half, which has the
same orbits (:ref:`manifold-orbit-space-stabiliser`) — is
one-dimensional in every degree (:ref:`manifold-descending-slots`).

Two families with different ranks make the layout a *variable*, so it
needs an owner, and the owner is the **layout value the basis carries**
(ruling 81). The basis ranks below the axis in the numerics concept order,
so it cannot mint its axis; it carries
:attr:`Basis.layout <orpheus.numerics.basis.base.Basis.layout>`, a value
from :mod:`orpheus.numerics.coefficient_layout`, and the one generic
coefficient axis carries that value as its ``layout`` slot. The family is
the layout, not a class: there is no harmonic or Legendre axis subclass and
no head space class. A consumer reads the layout through **one door**,
:func:`~orpheus.numerics.spaces.moment_head.moment_head_axis`, which
returns the head's single axis when its layout has degree structure and
otherwise raises a ``TypeError`` naming the space and the asking owner;
:func:`~orpheus.numerics.spaces.moment_head.head_layout` is that axis's
layout. `[M]` 2026-10-22, ``git grep`` over ``orpheus/``: the door's
production readers are the moment field's head check and its
``head_layout`` (:mod:`orpheus.transport.fields._bases`), the fission
moment operator (:mod:`orpheus.transport.operators.fission`) and the
transfer operators' head read (:mod:`orpheus.transport.operators.transfer`),
which spells the material field's per-degree contraction.

.. list-table:: What the layout answers, per family
   :header-rows: 1
   :widths: 30 34 36

   * - Question
     - ``HarmonicLayout(L)``
     - ``LegendreLayout(L)``
   * - ``shape``
     - :math:`(L+1,\ 2L+1)`
     - :math:`(L+1,)`
   * - ``isotropic_slot``
     - ``(0, 0)``
     - ``(0,)``
   * - ``degree_block(l)``
     - ``(l, 0:2l+1)``
     - ``(l,)``
   * - ``unit`` (the coefficients of the constant function :math:`1`)
     - :math:`e_0` at the isotropic slot (:math:`Y_{00} = 1`)
     - :math:`e_0` (:math:`P_0 = 1`)
   * - ``spread(per_degree)``
     - a per-degree value over each degree's :math:`2\ell+1` live
       columns, zero on the padding
     - the row itself
   * - rank ⟹ where the group axis is
     - 2 ⟹ ``values.shape[2]``
     - 1 ⟹ ``values.shape[1]``
   * - continuum metric (``coefficient_weights``)
     - :math:`4\pi/(2\ell+1)` per degree, spread over the live columns
     - :math:`4\pi/(2\ell+1)` per degree

Two further layouts index coefficient spaces that are not moment heads,
so the door refuses them: :class:`~orpheus.numerics.coefficient_layout.CellLayout`,
one coefficient per cell of a partition, which carries the partition's
per-axis edges because the cells are its atoms, and
:class:`~orpheus.numerics.coefficient_layout.ModeLayout`, a row of mode
amplitudes computed from data (the loss-kernel gauge's modes), whose unit
the basis declares.

.. _spaces-moment-head-axis-built:

The head is AXIS-BUILT, and the head axis's measure IS the metric
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Every head carries exactly one axis, and that axis's content carries the
whole identity of the head:

.. list-table:: The head axis, on the two spellings of a head
   :header-rows: 1
   :widths: 22 39 39

   * -
     - the basis's continuum head, ``FunctionSpace.of_basis(basis)``
     - the frame-dressed head, ``basis_space(frame)``
   * - origin
     - ``FromBasis(basis)``
     - ``FromFrame(frame)``
   * - ``label``, ``shape``, ``layout``
     - the layout's (``"harmonic"`` or ``"legendre"``)
     - the frame's trial basis's layout's
   * - ``weights``
     - the continuum Gram, ``basis.coefficient_weights``
     - the frame's Parseval inverse of its discrete Gram, or ``None`` with
       :math:`G^{+}` positioned as a dense overlay where the discrete Gram
       is dense (:ref:`spaces-metric-three-sources`)
   * - ``kind``
     - ``Modal(e_0)``
     - ``Modal(Φᵀw)``, the moments of the constant function
   * - ``support``
     - the basis's domain: :math:`S^2`, :math:`S^2/\sigma_a`, or
       :math:`S^2/O(2)_a`
     - the trial basis's domain

``Modal`` because a spectral coefficient may be negative for a positive
function, so neither spelling has a coordinate cone.

⭐ **Which pole a Legendre head was folded about is its domain.** The tree
carries two poles, so the Legendre heads about :math:`x` and about
:math:`z` are two physically different spaces, and their continuum
weights are ``array_equal``. The domain :math:`S^2/O(2)_a` is a content
part of the axis (ruling 71), so the two heads separate through it (`[M]`
2026-10-22; gate ``tests/gates/numerics/test_heads_on_the_axis.py`` H6).
The same part separates the harmonics on :math:`S^2` from their mirror-even
restriction (#631 pair 3, H1). No field is added for the pole: the
domain already says it, which is why ruling 81 dissolved the separate
``spent_axis`` field the Legendre head carried before S7.

⚠ **The measure enters the identity, and that is the whole point.** Two
heads of one order that differ only in their metric are two spaces —
which is what lets the corpus's chartered doctrine (*metric differences
imply space differences*) reach the moment layer, and what forced the
choice of WHICH metric the one moment space carries. That choice was
ruled Parseval; the full account, with both ends priced, is
:ref:`frame-the-one-moment-space`.

⭐ **Provenance is the object that spans the head.** ``FromBasis`` carries
the basis that was passed to ``of_basis`` and ``FromFrame`` the frame that
dressed it, never a fresh family instance built from the order: a
σ-even restriction's head must record the σ-even basis, or the origin
answers a question about some other object. `[M]` 2026-09-08, a first
spelling that minted a fresh ``SphericalHarmonicBasis(L=L)`` as the
provenance reddened 48 gate rows asserting that the head's recorded
object is ``frame.basis``.

**Truncation re-mints through the origin.**
:func:`~orpheus.numerics.spaces.moment_head.truncated_head` lowers a head
to ``L_new`` by matching its axis's origin: ``FromFrame`` re-poses the
frame at the lower order and returns its dressed head
(:meth:`GalerkinFrame.at_order
<orpheus.numerics.frame.GalerkinFrame.at_order>`), ``FromBasis``
re-spans the basis (``basis.at_order(L_new)``) and returns its
``of_basis``; any other origin refuses. Never a slice of the parent's
dressing: the discrete Gram's verdict can flip with :math:`L` (`[M]`
``folded_product(2, 4)`` is DENSE at :math:`L = 2` and DIAGONAL at
:math:`L = 1`), so a slice is undefined there. `[M]` 2026-10-22,
``truncated_head`` of the :math:`L = 3` dressed head of
``level_symmetric(6)`` equals the frame's own :math:`L = 1` dressed head
(gate H7).

**A zero modal unit is legal** (ruling 86, amending ruling 41's
"nonzero"). It is the honest unit of a factor whose span is orthogonal to
the constants, and then the marginal over the factor integrates every
field to 0. A family whose span contains :math:`1` mints its unit nonzero
by construction. The loss-kernel gauge's basis DECLARES its unit zero
(ruling 87): `[M]` 2026-10-21, :math:`|\Phi^{\mathsf T} G\,1|` is at most
1.8e-16 of its scale on each of 102 blocks over six reflective boxes and
exactly 0 on only 12 of them, so a computed unit would be round-off and
the coefficient space's identity would hang on arithmetic order. The
witness is ``tests/gates/sn/operators/test_loss_kernel_gauge.py``
``::test_the_kernel_modes_are_orthogonal_to_the_constants``; the proof is
owed (#636).

⭐ **The two continuum metrics agreeing is not a coincidence and it is the
reason the descent is an ISOMETRY rather than merely an isomorphism.** The
Legendre Gram is taken against the **pushforward** of :math:`d\Omega`
along the quotient map — :math:`\pi_*\,d\Omega = 2\pi\,d\mu` by
Archimedes' hat-box, which is what
:func:`~orpheus.numerics.generating_measure.pushforward_reference`
carries — so

.. math::
   :label: spaces-legendre-pushforward-gram

   \int_{S^2/O(2)_a} P_\ell^2 \; \mathrm{d}(\pi_*\Omega)
   \;=\; \int_{-1}^{1} P_\ell(\mu)^2 \, 2\pi \,\mathrm{d}\mu
   \;=\; \frac{4\pi}{2\ell+1},

exactly the harmonics'
:attr:`metric_per_ell <orpheus.numerics.basis.spherical_harmonic_basis.SphericalHarmonicBasis.metric_per_ell>`.
`[M]` 2026-10-22, ``LegendreBasis(L=2, axis="x")``'s coefficient axis
carries :math:`[12.5664, 4.1888, 2.5133] = 4\pi\,[1, 1/3, 1/5]`.

.. warning::

   ⚠ **NOT the bare Legendre mass-2 normalisation**
   :math:`2/(2\ell+1)`, which is the Gram against :math:`\mathrm{d}\mu`
   with no :math:`2\pi`. The two differ by a factor :math:`2\pi` and the
   wrong one would move **every operator end's metric** on every 1-D
   solve. The discriminator is which measure the orbit space carries,
   and the orbit space carries the pushforward of the sphere's — a
   1-D angular rule is a rule on :math:`S^2/O(2)_a`, not on an
   abstract interval.

   ⚠ And do not confuse either with the **discrete** Gram a frame
   measures, which is against the RULE's weights: `[M]` a
   Gauss–Legendre rule sums to :math:`W = 2`, so its discrete Legendre
   Gram is :math:`2/(2\ell+1)` and its Parseval-dressed ``basis_space``
   carries the inverse, :math:`(2\ell+1)/2`. Continuum, discrete and
   Parseval-dressed are three different arrays on this space exactly as
   they are on the harmonic one (:ref:`spaces-metric-three-sources`).

.. (vv-status rationale) A normalisation identity: it states which Gram
   the Legendre head carries and why, not what a solver computes. Its
   verifiable content is ``LegendreBasis.coefficient_weights``, sourced
   from ``LegendreBasis.metric_per_ell`` so the formula has one home, and
   the bit-identity of the descent, which would fail under the mass-2
   normalisation.
.. vv-status: spaces-legendre-pushforward-gram documented

⛔ **Why the door is a repair and not decoration: on a flat head a read
of the rectangular layout returns the wrong array and raises NOTHING.** A
consumer indexing ``values[0, 0]`` on an :math:`(L+1, n_g, n_x)` tensor
gets *group 0's spatial slice* — well-shaped, silently wrong. The sites
that read the rectangular layout as if it were the contract until #429
(2026-09-02): ``scalar_flux``, ``isotropic_part``, ``anisotropic_part``,
``l_block``, the fission :math:`\ell = 0` dyad, ``ng`` (which located the
group axis at a hard-coded index 2), ``zeros_for_problem_and_L``, and the
material field's per-degree group contraction, whose ``einsum`` spelled
the :math:`m` axis into its subscripts and would have contracted the GROUP
axis as if it were :math:`m`. All of them read the layout now, and a head
whose layout has no degree structure is refused by name rather than
contracted wrongly. The frame-side account is :ref:`frame-g0-descent-arrow`.

The spatial-moment TAIL of a widened moment space (LD,
``spatial_moments > 1``) is the discretization scheme's own ``Declared``
``Modal(e_0)`` :meth:`moment_axis
<orpheus.transport.spatial.scheme.DiscretizationSchemeBase.moment_axis>`
(:ref:`the account in the SN chapter <spatial-moment-space>`). `[M]`
2026-09-08, 1-D LD under ``gauss_legendre(2)``: the widened product's
axes read ``[legendre, energy, spatial, spatial_moment]`` at :math:`L = 1`
and :math:`L = 2` alike.

The design record of the two head classes and the ``MomentHead``
protocol they satisfied until S7 is this page's text at ``main`` commit
``dd17f678``, lines 3703 to 3938; the milestone is a row of
:ref:`spaces-development-history`.


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
differ only in measure — the measure is part of axis identity, and since
the identity flip (structural ``__eq__``, CS4c step 6, 2026-09-07) an
axis-built space IS its axis tuple, so they are **unequal spaces**.
Until the flip the same conclusion arrived one step later, through the
axis-derived name (:ref:`spaces-identity-bridge`). Nothing else in the
tree can see the difference: `[M]` a scalar metric commutes with every
operator, so
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
   Read this subsection as *why the quotient is licensed one level above
   the mesh*, not as a description of a class that exists.

   ⛔ **Corrected 2026-09-08 (CS4c coda).** This warning used to add
   that "today the homogeneous path still constructs a degenerate
   :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` through
   ``from_materials``, which is the pullback wearing a constructor story
   it does not have a mesh for." That is now history: the factory is
   deleted and the homogeneous path constructs no carrier at all — its
   problem hub
   (:class:`~orpheus.homogeneous.solver.HomogeneousProblem`) mints the
   quotient pose and every field on it straight from the
   :class:`~orpheus.data.macro_xs.mixture.Mixture`
   (:doc:`infinite_medium`). So the *misplacement* the sentence
   complained about is gone, and the ``Medium`` layer's own charter is
   untouched by that: what the infinite path has today is an
   **aggregate** over a mixture, not the stage-1 medium this law argues
   for — the lattice above still describes a class that does not
   exist.



.. _spaces-collapse-pair:

The realized machinery: the axis marginal
=========================================

The doctrine above decides **which** axes survive a degeneracy. This
section is the machinery that performs the collapse when the verdict is
"drop", and it is where the doctrine stops being a rule and becomes three
typed arrows a call site can hold.

The bridge is the retract rule's own sentence
(:ref:`spaces-collapse-the-question`): *a projection*
:math:`\pi: V \to V'`, *an embedding* :math:`\iota: V' \to V`, *with*
:math:`\pi\circ\iota = \mathrm{id}` *and, under the right metrics,*
:math:`\iota = \pi^{H}` *up to a scalar*. Everything below is that
sentence made precise and made executable:

- :math:`\pi` is the retraction
  :class:`~orpheus.numerics.axis_marginal.AxisRetractionOperator`, built by
  :meth:`AxisRetractionOperator.of(space, label)
  <orpheus.numerics.axis_marginal.AxisRetractionOperator.of>`;
- :math:`\iota` is the section
  :class:`~orpheus.numerics.axis_marginal.AxisSectionOperator`, built by
  :meth:`AxisSectionOperator.of(space, label)
  <orpheus.numerics.axis_marginal.AxisSectionOperator.of>`;
- :math:`\pi^{H}`, the Hilbert adjoint of the retraction, is the pullback
  :class:`~orpheus.numerics.axis_marginal.AxisPullbackOperator`, returned by
  ``R.H`` and minted with :math:`R` (:ref:`spaces-collapse-pair-pullback`);
- all three are built from one datum the axis owns, its **unit**
  :math:`u`, the coordinates of the constant function :math:`1` in the
  axis's own representation (:attr:`Axis.unit
  <orpheus.numerics.axis.Axis.unit>`; :ref:`spaces-axis-unit`), and from
  the axis's metric block :math:`G_k`;
- **the scalar is** :math:`\langle m, u\rangle`, the mass of the measure
  the axis carries (:math:`\Sigma w` on a nodal axis, :math:`4\pi` on a
  harmonic head), and it is not a convention: it is the integral of the
  constant function, :eq:`spaces-axis-marginal-integral` at :math:`c = u`;
- the section and the adjoint are **different types** precisely so
  that the "up to a scalar" cannot be silently dropped at a call site,
  which is the ERR-051 failure class made unspellable.

The three arrows live in :mod:`orpheus.numerics.axis_marginal`, a module
ranked after the operator module (the concept order of the numerics
layer); the space itself carries no lift verb. The constructors are
downstream of the space: ``.of(space, label)`` takes the space and the
axis label.

.. code-block:: python

   import numpy as np

   from orpheus.numerics.axis import Axis, BasisKind
   from orpheus.numerics.axis_marginal import (
       AxisPullbackOperator,
       AxisRetractionOperator,
       AxisSectionOperator,
   )
   from orpheus.numerics.space import FunctionSpace

   w = np.array([0.3, 0.7, 0.5, 0.5])
   V_cell = np.array([0.2, 0.3, 0.4, 0.7, 1.4])
   V = FunctionSpace.of_axes(
       Axis("angular", (4,), weights=w, kind=BasisKind.NODAL),
       Axis("energy", (2,), kind=BasisKind.NODAL),
       Axis("spatial", (5,), weights=V_cell, kind=BasisKind.NODAL),
   )
   R = AxisRetractionOperator.of(V, "angular")  # V -> energy (x) spatial (fibre integration)
   E = AxisSectionOperator.of(V, "angular")     # energy (x) spatial -> V (the section)

   phi = np.random.default_rng(0).standard_normal(R.codomain.shape)
   np.testing.assert_array_almost_equal_nulp(R.apply(E.apply(phi)), phi, nulp=4)
   assert isinstance(R.H, AxisPullbackOperator)    # the pullback, NOT E
   assert R.H.H is R
   assert E.mass == 2.0                     # the mass <m, u> = sum(w)

The three arrows are **memoized on the space**, one mint per space per
axis label, so a carrier that caches its spaces (every solver carrier
does) gets warm operators: ``AxisRetractionOperator.of(sn.angular_bulk_space,
"angular")`` mints once and returns the same object thereafter, and the
retraction, the section and the pullback share ONE marginal-space
instance (``R.codomain is E.domain``; `[M]` identity, gated by
``test_the_constructors_memoise_one_marginal``). The memo lives in the
space instance's private slot ``"_axis_marginal"``
(:func:`orpheus.numerics.instance_cache.cached_on`), not on the axis: the
marginal is not a function of the axis alone, since its domain is the
whole product and its codomain the product minus that factor.


.. _spaces-axis-unit:

The unit, and the covector of integration
-----------------------------------------

Marginalising over an axis is integrating along it, and integrating a
function is pairing it with the constant function :math:`1`. So the one
datum the marginal needs from the axis is **the coordinates of** :math:`1`
**in the axis's own representation**, which the axis owns as its unit
:math:`u` (:attr:`Axis.unit <orpheus.numerics.axis.Axis.unit>`). Write the
axis's coordinates as a vector :math:`c \in \mathbb R^{N_k}` (one entry per
slot of the axis), the function they represent on the axis's domain
:math:`D_k` as :math:`f_c`, the measure the axis carries as
:math:`\mu_k`, and the axis's metric block as :math:`G_k`, so that
:math:`c^{\mathsf T} G_k c' = \langle c, c'\rangle_{G_k} = \int_{D_k} f_c\,
f_{c'}\,\mathrm d\mu_k`. The unit is defined by :math:`f_u = 1`. Then for
every :math:`c`,

.. math::
   :label: spaces-axis-marginal-integral

   \langle m, c\rangle
   \;=\; m^{\mathsf T} c
   \;=\; u^{\mathsf T} G_k\, c
   \;=\; \langle u, c\rangle_{G_k}
   \;=\; \int_{D_k} 1\cdot f_c \,\mathrm d\mu_k
   \;=\; \int_{D_k} f_c \,\mathrm d\mu_k ,
   \qquad
   m \;\equiv\; G_k\, u .

.. implements:: spaces-axis-marginal-integral
   :by: py:method:orpheus.numerics.axis_marginal.AxisRetractionOperator.apply

   **Implemented by** 2 sites. The retraction's ``apply`` is the pairing
   :math:`\langle m, c\rangle` along the axis; the mint forms the covector
   :math:`m = G_k u` from the axis's metric block,
   ``FunctionSpace._axis_block(k)``, the same entry the space's derived
   metric composes (a positioned form on the space's overlay when the axis
   is measure-less, else the axis's weights, else ``None``, the counting
   measure), so the covector and the space's metric cannot disagree.

.. implements:: spaces-axis-marginal-integral
   :by: py:function:orpheus.numerics.axis_marginal._mint

The vector :math:`m` is the **covector of integration**: the axis's metric
applied to its unit. Each step of :eq:`spaces-axis-marginal-integral` is
one fact: the definition of :math:`m`; the definition of the
:math:`G_k`-inner product; the hypothesis that :math:`G_k` is the inner
product the axis's measure induces on the represented functions; and
:math:`f_u = 1`. At :math:`c = u` it reads :math:`\langle m, u\rangle =
\mu_k(D_k)`, the **mass** of the axis's measure, which is the section's
divisor below. The equation is the **modal marginal**: before the unit
existed, an axis of ``MODAL`` kind was refused at the mint, because
contracting its coefficients with its weights is not an integral of the
function they represent (:ref:`spaces-collapse-pair-refuted`). With the
unit, the contraction is :math:`m^{\mathsf T} c` with :math:`m = G_k u`,
and that is the integral on every axis whose unit is right. The
verifiable content is the integration family of
``tests/gates/numerics/test_axis_marginal_lift.py`` (G2,
:ref:`spaces-collapse-pair-gates`): on the coefficient heads, the
retraction against an integral computed in the test with a
:math:`40\times80` product rule on the sphere and an independent
evaluation of the basis; on the moment heads, the commuting square
below.

**The unit on every kind of axis.** The unit is set by whatever mints the
axis, like the weights (ruling 37 of the numerics-untangling plan): the
forgetful image of a stage-2 generator keeps the induced data, and the
unit is one of them (:ref:`spaces-collapse-pair-frame`).

.. list-table:: The unit, the covector and the mass, per axis kind
   :header-rows: 1
   :widths: 22 22 14 14 14 14

   * - Axis (where it is minted)
     - Coordinates
     - :math:`G_k`
     - :math:`u`
     - :math:`m = G_k u`
     - mass :math:`\langle m, u\rangle`
   * - nodal, with a measure :math:`w` (an angular rule, cell volumes)
     - point values :math:`f(x_n)`
     - :math:`\operatorname{diag}(w)`
     - :math:`\mathbf 1`
     - :math:`w`
     - :math:`\Sigma w`
   * - nodal, counting measure (``weights=None``)
     - point values
     - :math:`I`
     - :math:`\mathbf 1`
     - :math:`\mathbf 1`
     - :math:`N_k`
   * - continuum Legendre head
       (:meth:`FunctionSpace.of_basis
       <orpheus.numerics.space.FunctionSpace.of_basis>` of a
       ``LegendreBasis``)
     - coefficients :math:`c_\ell` of :math:`\sum_\ell c_\ell P_\ell`
     - :math:`\operatorname{diag}\bigl(4\pi/(2\ell+1)\bigr)`
     - :math:`e_0`
     - :math:`4\pi\,e_0`
     - :math:`4\pi`
   * - continuum harmonic head
       (:meth:`FunctionSpace.of_basis
       <orpheus.numerics.space.FunctionSpace.of_basis>` of a
       ``SphericalHarmonicBasis``)
     - coefficients :math:`c_{\ell m}` of :math:`\sum c_{\ell m} Y_{\ell m}`
     - :math:`4\pi/(2\ell+1)` per :math:`\ell`
     - :math:`e_0` (slot :math:`(0,0)`)
     - :math:`4\pi\,e_0`
     - :math:`4\pi`
   * - a head a frame has dressed
       (``frame_faces._dressed``; every ``SNProblem.moment_space(L)``)
     - moments :math:`\varphi = \Phi^{\mathsf T} W \psi`
     - :math:`B^{+}`, :math:`B = \Phi^{\mathsf T} W \Phi`
     - :math:`\Phi^{\mathsf T} w = B\,e_0`
     - :math:`B^{+} B\, e_0`
     - :math:`\Sigma w` (to rounding)
   * - the LD scheme's spatial-moment axis
       (``DiscretizationSchemeBase.moment_axis``)
     - Legendre coefficients on the reference cell
     - :math:`\operatorname{diag}(1, \tfrac13, \dots)`
     - :math:`e_0`
     - :math:`e_0`
     - :math:`1` (the cell average)

The continuum heads and the moment heads use the same harmonics in the
convention where :math:`P_0 = Y_{00} = 1` (no orthonormalising prefactor;
the :math:`\ell = 0` metric entry is :math:`4\pi`,
:doc:`/theory/foundations/spherical_harmonics`), so :math:`e_0` represents
the constant function exactly. The continuum Legendre head carries the
sphere's measure pushed to :math:`\mu` (`[M]` its weights are
:math:`4\pi/(2\ell+1)`: :math:`12.566, 4.189, 2.513, 1.795` at
:math:`L = 3`), so its mass is :math:`4\pi` and not :math:`2`.

**The two variances of a modal head.** Two heads of the same order and the
same harmonics can hold different coordinates of the same function, and
the unit is what tells them apart.

- A **coefficient head** (minted from a basis) holds the *contravariant*
  expansion coefficients :math:`c` of :math:`f_c = \sum_i c_i \phi_i`, and
  its metric is the continuum Gram :math:`\int \phi_i \phi_j\,\mathrm d\mu`.
  The constant is :math:`\phi_0`, so :math:`u = e_0`.
- A **moment head** (a head a frame has dressed) holds the *covariant*
  moments :math:`\varphi_i = \sum_n w_n\,\phi_i(\Omega_n)\,\psi_n`, the
  pairings of the ordinate field :math:`\psi` with the basis, and its
  metric is the pseudo-inverse of the discrete Gram
  :math:`B = \Phi^{\mathsf T} W \Phi` (:math:`\Phi_{ni} =
  \phi_i(\Omega_n)` the frame's table, :math:`W = \operatorname{diag}(w)`;
  :ref:`frame-parseval-metric`). The coordinates of the constant on such a
  head are the moments of :math:`1`, :math:`u = \Phi^{\mathsf T} W
  \mathbf 1 = \Phi^{\mathsf T} w`, and because column 0 of :math:`\Phi`
  is :math:`\mathbf 1`, :math:`u = B\,e_0`. That vector is
  :math:`(\Sigma w)\,e_0` only on a rule that integrates every live
  harmonic of the head exactly; on ``folded_product(4, 4)`` at
  :math:`L = 4` its other entries reach 3.67 (`[M]` 2026-10-19, the unit table of
  :ref:`spaces-collapse-pair-evidence`).

The weights alone do not fix which of the two a head holds: two dense
frame-dressed heads of the same order, with ``weights=None`` and the form
positioned on the space, would be EQUAL axes with different units
(:math:`B e_0` on different rules) if the unit were not content. So the
unit lives inside the kind: the kind is the closed sum
``BasisKind.NODAL | Modal(unit)`` (ruling 41 of the numerics-untangling
plan), :class:`~orpheus.numerics.axis.Modal` is a content-identified
value carrying the unit, and the kind is part of the axis's content. On
a modal axis the unit therefore enters ``==``, ``hash`` and the digest;
on a nodal axis, whose kind is the bare enum member, there is nothing to
enter, so every nodal axis keeps the digest it had before the unit
existed (ruling 39, which an override of ``content_parts`` enforced
until ruling 41 made it a property of the type).

**Why the moment head integrates.** On a moment head,
:eq:`spaces-axis-marginal-integral` has to be checked directly, because
:math:`G_k = B^{+}` is the inner product of the frame's *analysis*, not a
Gram of represented functions. Let :math:`P = B^{+}B`, the orthogonal
projector onto :math:`\operatorname{range}(B)`. For a stored moment vector
:math:`\varphi = \Phi^{\mathsf T} W \psi`,

.. math::

   \langle m, \varphi\rangle
   \;=\; e_0^{\mathsf T} B\, B^{+} \varphi
   \;=\; e_0^{\mathsf T} P\, \varphi
   \;=\; e_0^{\mathsf T} \varphi
   \;=\; \sum_n w_n\,\phi_0(\Omega_n)\,\psi_n
   \;=\; \sum_n w_n\,\psi_n ,

where the first step is :math:`m = B^{+} B e_0` with :math:`B` and
:math:`B^{+}` symmetric, and the third needs :math:`P\varphi = \varphi`,
that is :math:`\operatorname{range}(\Phi^{\mathsf T} W) \subseteq
\operatorname{range}(B)`. With positive weights this always holds:
:math:`B = (W^{1/2}\Phi)^{\mathsf T}(W^{1/2}\Phi)`, so
:math:`\operatorname{range}(B) = \operatorname{range}(\Phi^{\mathsf T})
\supseteq \operatorname{range}(\Phi^{\mathsf T} W)`. 117 of the 120
shipped rules of the population below have only positive weights; the
other three (``lebedev(13)``, ``lebedev(25)``, ``lebedev(27)``, with 8, 6
and 18 negative weights) still satisfy it at :math:`L = 2, 4, 6`, where
`[M]` the marginal of the moments matches :math:`\sum_n w_n \psi_n` to at
most :math:`2.0\times10^{-15}` of its largest entry (the unit table of
:ref:`spaces-collapse-pair-evidence`). So the marginal over a moment
head is the marginal over the ordinates it was analysed from, the
**commuting square** :math:`R_{\rm head}(\Phi^{\mathsf T} W \psi) =
R_{\rm angular}(\psi)` that gate G2b pins with the moments computed in
the test. The mass is :math:`\langle m, u\rangle = e_0^{\mathsf T} P B
e_0 = B_{00} = \Sigma w` in exact arithmetic, and to rounding in the code
(:ref:`spaces-collapse-pair-evidence`).

.. dropdown:: First got wrong: the unit read through the generator basis
   :color: muted

   The first design (2026-10-19, before any code) delegated the unit of a
   modal axis to its generator basis, a new verb
   ``Basis.constant_coefficients()`` returning :math:`e_0` per family.
   The test-architect's specification measured it before it was built.
   On a moment head :math:`u = e_0` gives :math:`\langle m, \varphi\rangle
   = e_0^{\mathsf T} B^{+}\varphi`, which on a rule with
   :math:`B e_0 = (\Sigma w) e_0` is :math:`(1/\Sigma w)\int\psi\,\mathrm
   d\mu`: `[M]` 0.5 on Gauss–Legendre and :math:`1/(4\pi) = 0.0796` on
   the sphere rules, and on ``folded_product(4, 4)`` at :math:`L = 4`
   not a multiple of the integral at all (the ratio ranges from
   :math:`-0.11` to :math:`0.18` over 20 random columns, `[M]` 2026-10-19). The basis knows its coefficients of :math:`1`;
   it does not know whether a head built on it holds coefficients or
   moments. That fact is decided where the head is dressed, so the unit
   is minted there (ruling 37), and the reciprocity and section gates
   could never have seen the error, because they hold for every nonzero
   :math:`u` (:ref:`spaces-collapse-pair-gates`).

**A modal axis is minted with its unit, or not at all** (rulings 40 and
41). The modal arm of the kind is ``Modal(unit)``, so a modal factor
without its unit cannot be spelled. What remains to refuse is refused at
construction: a kind that is neither arm (``TypeError``, *"kind must be
BasisKind.NODAL or Modal(unit)"*; since S7 the kind is read off the
origin, so this arm is reached only through ``Declared``), a unit that is
not finite (``ValueError``, *"a modal unit must be finite"*, in
``Modal``), and a unit over another shape than the axis's (``ValueError``,
*"must live over exactly the axis shape"*). So a modal factor that cannot
integrate is unrepresentable. A ZERO unit is legal (ruling 86, amending
ruling 41's "nonzero"): it is the honest unit of a factor whose span is
orthogonal to the constants, the marginal over it integrates every field
to 0, and its section does not exist, because the mass
:math:`\langle m, u\rangle` is 0. The mint treats a mass as zero when it
is round-off of its own terms' scale (at most :math:`64\,\varepsilon`
times :math:`\sum |m_i u_i|`), so a section never divides by a round-off
mass (the private mint ``axis_marginal._mint``). The one
refusal the marginal itself makes on the origin is the raw partition head,
which is measure-free (ruling 88, :ref:`spaces-collapse-pair-frame`). A planned refusal for a basis
that cannot represent the constant was dropped: `[M]` (the S1
specification, by reading the families) 0 of the 5 shipped
families with a synthesis lacks the constant (Legendre, the spherical
and the mirror-even harmonics carry degree 0, an indicator partition
covers its support, an overlap basis is a partition of unity).


.. _spaces-collapse-pair-two-arrows:

Three arrows, and the mass between them
---------------------------------------

Write :math:`V = V_k \otimes V'` for the product whose first factor is
the axis being collapsed, with the product metric
:math:`G_V = G_k \otimes G_{V'}` (:eq:`spaces-axis-product`), and a field
as its slices :math:`\psi = (\psi_n)_n`, :math:`\psi_n \in V'`. The
**retraction** is fibre integration over that factor, the pushforward
:math:`\pi_*` along the projection that forgets the axis:

.. math::
   :label: spaces-collapse-retraction

   (R\,\psi)(\cdot) \;=\; \sum_{n} m_n\,\psi(n,\cdot),
   \qquad
   R \;=\; m^{\mathsf T}\otimes \mathrm{id}_{V'} ,
   \qquad
   R : V \longrightarrow V' .

On a nodal axis :math:`m = w` and this is the weighted sum over the
axis: the scalar flux :math:`\phi = \sum_n w_n \psi_n` on the angular
axis, the volume integral on a spatial one.

.. (vv-status rationale) Structural/representational identity: it STATES
   what ``AxisRetractionOperator.apply`` computes (the covector m = G_k u
   contracted over the axis's ndarray dims). Not a solver claim: no
   flux, no eigenvalue, no discretization error. The verifiable content
   is the nodal battery ``tests/gates/numerics/test_axis_marginal.py``
   (the tightness row against the single-region indicator frame's own
   analysis content; G6.5 bit-identically against a hand-spelled einsum
   on the real S_N carrier) and the G5 rows of
   ``tests/gates/numerics/test_axis_marginal_lift.py`` (array_equal with
   the hand einsum over the raw weights, 16 nodal spaces). The modal
   content of the same contraction is the separate, testable label
   spaces-axis-marginal-integral.
.. vv-status: spaces-collapse-retraction documented

.. implements:: spaces-collapse-retraction
   :by: py:method:orpheus.numerics.axis_marginal.AxisRetractionOperator.apply

   **Implemented by** 2 sites. The kernel is the operator's; the
   canonical field-level consumer re-spells nothing: since CS4b S6.2
   :meth:`AngularField._integrate_angular_values
   <orpheus.transport.fields._bases.AngularField._integrate_angular_values>`
   IS this ``apply`` on the angular axis, so the tree has one
   realization of the angular reduction and not two.

.. implements:: spaces-collapse-retraction
   :by: py:method:orpheus.transport.fields._bases.AngularField._integrate_angular_values

Its Hilbert adjoint is where the three-arrow discipline is forced. Under
the product metric the adjoint is the **pullback**, the constant function
along the axis times the marginal field, with no metric in it:

.. math::
   :label: spaces-collapse-adjoint-is-pullback

   R^{\dagger} \;=\; \pi^{*} \;=\; u \otimes \mathrm{id}_{V'},
   \qquad
   (\pi^{*}\varphi)(n,\cdot) \;=\; u_n\,\varphi(\cdot) .

.. (vv-status rationale) Structural identity: under the product metric
   G_V = G_k (x) G_V', <R psi, phi> = <psi, u (x) phi>_G for every psi
   and phi (the derivation in spaces-collapse-pair-pullback), so the
   Hilbert adjoint of fibre integration is the pullback u (x) id, with no
   inverse of G. Representational, not a solver claim (no physics
   enters; it holds for every axis on which the product condition holds,
   and it is the unique adjoint when G_V is invertible). The verifiable
   content is ``tests/gates/numerics/test_retraction_adjoint_is_the_pullback.py``
   (S6.1: the type, the involution, array_equal with the broadcast; S6.2:
   the product condition and reciprocity, the weights read from the raw
   axes; S6.4: the generic sandwich as the X4 witness), the G3 family of
   ``tests/gates/numerics/test_axis_marginal_lift.py`` (reciprocity over
   38 (space, axis) rows, modal heads with dense forms included), and the
   G6.3 row of ``tests/gates/numerics/test_axis_marginal.py`` with its vv
   #19 NEGATIVE leg (the same pairing under a deliberately stripped
   spatial measure must break at O(1)).
.. vv-status: spaces-collapse-adjoint-is-pullback documented

.. implements:: spaces-collapse-adjoint-is-pullback
   :by: py:method:orpheus.numerics.axis_marginal.AxisPullbackOperator.apply

   **Implemented by** 3 sites. The pullback's ``apply`` is the right-hand
   side, :math:`u \otimes \varphi`; the retraction's ``adjoint`` is the
   equality itself, returning that operator as :math:`R^{\dagger}` (so
   ``R.H`` is the closed form, not the generic wrapper). The third site
   is a hand-written twin on the angular axis, where :math:`u = \mathbf
   1`: ``AngularSourceSink.__add__`` with an isotropic partner adds
   ``other.values[None]``, the same broadcast spelled outside its type
   (:ref:`spaces-collapse-pair-pullback`, the warning).

.. implements:: spaces-collapse-adjoint-is-pullback
   :by: py:method:orpheus.numerics.axis_marginal.AxisRetractionOperator.adjoint

.. implements:: spaces-collapse-adjoint-is-pullback
   :by: py:method:orpheus.transport.fields.angular_source_sink.AngularSourceSink.__add__

So :math:`(R, R^{\dagger}) = (\pi_*, \pi^{*})` is the discrete
realization of the fibre-integration / pullback adjunction, and the
pairing

.. math::

   \langle R\psi,\ \varphi\rangle_{V'}
   \;=\;
   \langle \psi,\ \pi^{*}\varphi\rangle_{V}

holds on the *physical* metrics with no correction factor at all
(`[M]` 2026-08-24 on nodal axes: bounded by :math:`3.7\times10^{-13}`
relative over 200 draws on the shipped S\ :sub:`N` carrier and
:math:`1.6\times10^{-13}` on the synthetic three-axis fixture, a
scalar-valued identity whose relative residual is set by cancellation in
the inner products, not by the operators; see
:ref:`spaces-collapse-pair-evidence`. Gate G3 normalises the residual by
:math:`\lVert G\psi\rVert\,\lVert\pi^{*}\varphi\rVert` instead and bounds
it at :math:`16\varepsilon` over 38 rows, modal heads included).

**The pullback is not the section.** :math:`\pi^{*}` places
:math:`\varphi` along the unit unchanged, so :math:`R\pi^{*}\varphi =
(m^{\mathsf T}u)\,\varphi = \langle m, u\rangle\,\varphi`: it is a right
inverse of :math:`R` only after division by the axis's mass. That
division is the third arrow:

.. math::
   :label: spaces-collapse-section

   (E\,\varphi)(n,\cdot) \;=\; \frac{u_n\,\varphi(\cdot)}{\langle m, u\rangle},
   \qquad
   E \;=\; \frac{u}{\langle m, u\rangle}\otimes\mathrm{id}_{V'}
   \;=\; \frac{\pi^{*}}{\langle m, u\rangle}
   \;=\; R^{\dagger}\,(R\,R^{\dagger})^{-1},
   \qquad
   R\circ E \;=\; \mathrm{id}_{V'} .

.. (vv-status rationale) Structural/representational identity: it STATES
   what ``AxisSectionOperator.apply`` computes (divide by the axis's mass
   <m, u>, then place along the unit) and identifies it as the
   Moore-Penrose right inverse of the retraction in the two spaces' own
   metrics. Not a solver claim. The verifiable content is the G4 family
   of ``tests/gates/numerics/test_axis_marginal_lift.py`` (R o E = id at
   nulp max(4, N_k) and E o R idempotent, over 38 rows), the G5 rows
   (bit identity with the nodal pair; the mass array_equal to math.fsum
   of the weights over 120 shipped rules), and the nodal battery in
   ``tests/gates/numerics/test_axis_marginal.py``: G6.1, G6.2, G6.6
   (bit-identity with the shipped isotropic-source kernel on the real
   S_N carrier).
.. vv-status: spaces-collapse-section documented

.. implements:: spaces-collapse-section
   :by: py:method:orpheus.numerics.axis_marginal.AxisSectionOperator.apply

   **Implemented by** 3 sites. The kernel is the operator's; the mint
   supplies the mass, summed correctly rounded; the canonical
   producer-side consumer re-spells nothing: since CS4b S6.2
   :meth:`AngularSourceSink.from_isotropic
   <orpheus.transport.fields.angular_source_sink.AngularSourceSink.from_isotropic>`
   IS this ``apply`` on the angular axis.

.. implements:: spaces-collapse-section
   :by: py:function:orpheus.numerics.axis_marginal._mint

.. implements:: spaces-collapse-section
   :by: py:method:orpheus.transport.fields.angular_source_sink.AngularSourceSink.from_isotropic

The three identities in :eq:`spaces-collapse-section` follow in one line
each: :math:`R E \varphi = (m^{\mathsf T}u)\,\varphi/\langle m,u\rangle =
\varphi`; :math:`R R^{\dagger} = R\pi^{*} = \langle m, u\rangle\,
\mathrm{id}`, a nonzero scalar, so :math:`(RR^{\dagger})^{-1} =
1/\langle m, u\rangle`; and :math:`E = \pi^{*}/\langle m, u\rangle`.

Since :math:`R` is rank-deficient by construction, *many* right inverses
exist. On a nodal axis, putting all the mass on one slot,
:math:`(\iota_0\varphi)(n,\cdot) = \delta_{n0}\,\varphi/w_0`, is one.
:math:`E` is distinguished among them by being the **minimum-norm** one,
which is what :math:`R^{\dagger}(RR^{\dagger})^{-1}` says: it is the
Moore–Penrose pseudo-inverse in the two spaces' own metrics. `[M]`
2026-08-24 on the synthetic fixture and one draw (``default_rng(21)``),
:math:`\lVert E\varphi\rVert_V = 1.700` against :math:`\lVert
\iota_0\varphi\rVert_V = 4.389` for the one-slot right inverse, a factor
2.6, with both satisfying :math:`R\circ(\cdot) = \mathrm{id}` at the
round-off floor.

The four arrows close into a square, and every entry is one of the
operators or one of their adjoints:

.. list-table:: The collapse square (:math:`M \equiv \langle m, u\rangle`, the mass)
   :header-rows: 1
   :widths: 22 30 48

   * - Arrow
     - What it is
     - Identity
   * - :math:`R = \pi_*`
     - fibre integration (the pairing with :math:`m` over the axis)
     - :math:`R\,R^{\dagger} = M\,\mathrm{id}`
   * - :math:`R^{\dagger} = \pi^{*}`
     - the pullback (:math:`u\otimes\varphi`; on a nodal axis the
       **plain** broadcast)
     - :math:`R^{\dagger} = M\,E`, the anti-ERR-051 scalar
   * - :math:`E`
     - the section (:math:`u\otimes\varphi` **after** dividing by
       :math:`M`)
     - :math:`R\circ E = \mathrm{id}`, and :math:`E` is minimum-norm
   * - :math:`E^{\dagger}`
     - the :math:`m`-mean (on a nodal axis the :math:`w`-average)
     - :math:`E^{\dagger} = R/M`, and :math:`E^{\dagger}\circ R^{\dagger} = \mathrm{id}`

and the composite in the other order,

.. math::

   P \;\equiv\; E\circ R,
   \qquad
   (P\psi)(n,\cdot) \;=\; \frac{u_n}{M}\sum_{n'} m_{n'}\,\psi(n',\cdot),

is the **conditional expectation onto the functions constant along the
axis**: it is idempotent and, because :math:`P^{\dagger} =
R^{\dagger}E^{\dagger} = (ME)(R/M) = P`, self-adjoint in :math:`G_V`, an
*orthogonal* projector, not merely an oblique one (`[M]` 2026-08-24 on
nodal axes, self-adjointness bounded by :math:`4.0\times10^{-14}`
relative over 200 draws, again a cancellation-limited scalar identity,
see :ref:`spaces-collapse-pair-evidence`). That is the precise sense in
which "the isotropic part of :math:`\psi`" is a well-defined object
rather than a convention.

**The mass is summed correctly rounded.** The divisor is
:math:`\langle m, u\rangle`, computed as ``math.fsum(m * u)`` over the
flattened axis (ruling 38): the correctly rounded sum of the products
:math:`m_n u_n`. On a nodal axis each product is :math:`w_n` exactly, so
the divisor is the correctly rounded :math:`\Sigma w`, a named quantity
with a platform-free value; `[M]` ``array_equal`` to ``math.fsum(w)`` on
120 of 120 shipped rules (gate G5-mass). Every other spelling of the sum
rounds differently on some rule, and the spelling the pair used before
this ruling was the least accurate of them
(:ref:`spaces-collapse-pair-evidence`, the divisor table). On a moment
head the products :math:`m_n u_n` are themselves rounded, so the divisor
is :math:`\Sigma w` to a few ULP, not exactly.

**And** :math:`R\circ E` **is not bit-exact.** :math:`R E \varphi =
\sum_n m_n\,(u_n\varphi/M)` re-associates a sum of :math:`N_k` terms, so
it reproduces :math:`\varphi` to the reduction depth and no better:
`[M]` 2026-10-19 at ``7baf444c``, at most :math:`N_k/3` ULP on each of the 120 shipped
rules (the largest count 67 ULP on ``lebedev(41)``, :math:`N_k = 590`),
and bit-exact on all 20 draws on only 2 of the 120 rules
(:ref:`spaces-collapse-pair-evidence`). The gate's tolerance is therefore
``nulp = max(4, N_k)`` (G4), and a consumer that needs :math:`R \circ E =
\mathrm{id}` (the DSA low-order map) asserts it at :math:`N_k\varepsilon`
relative.

.. warning::

   **An axis whose measure has zero mass has NO section, and the
   retraction over it is still legal.** The asymmetry is structural, not
   defensive: :math:`R` is a contraction and needs no division, while
   :math:`E` divides by :math:`\langle m, u\rangle`, and at zero mass
   (a signed measure summing to zero) no right inverse of the form
   :math:`R^{\dagger}(RR^{\dagger})^{-1}` exists. The mint therefore
   leaves that arm unminted, and :meth:`AxisSectionOperator.of
   <orpheus.numerics.axis_marginal.AxisSectionOperator.of>` refuses at
   access with a ``ValueError`` naming the cause (*"zero total weight"*).
   Signed axis weights are deliberately legal on
   :class:`~orpheus.numerics.axis.Axis` (a :math:`\sigma`-folded
   quadrature can carry them, and three shipped Lebedev rules do), so
   this is a reachable state and not a theoretical one.

.. note::

   **This is where the doctrine's "up to a scalar" gets its value.**
   The retract rule said :math:`\iota = \pi^{H}` *up to a scalar* and
   left the scalar unnamed, which is exactly the gap the classical
   :math:`4\pi` bookkeeping errors live in. Here the scalar is the mass
   :math:`\langle m, u\rangle`, the integral of the constant function,
   and the arrows carry different **types**, so a call site that reaches
   for ``R.H`` where it wanted ``E`` does not silently rescale a source
   by :math:`4\pi`; it holds an object of the wrong class. `[M]` the
   gate asserts precisely that ``R.H`` is *not* an
   :class:`~orpheus.numerics.axis_marginal.AxisSectionOperator`, and that
   it *is* an :class:`~orpheus.numerics.axis_marginal.AxisPullbackOperator`.


.. _spaces-collapse-pair-pullback:

The third arrow: the adjoint is minted as the pullback
------------------------------------------------------

``R.H`` returns :class:`~orpheus.numerics.axis_marginal.AxisPullbackOperator`,
minted from the same induced data as :math:`R` and :math:`E` (the bound
spaces, the ndarray dims of the axis, the unit). Its ``apply`` is
:math:`u\otimes\varphi` (on a nodal axis the plain broadcast) and its
``apply_transpose`` the pairing with :math:`u` over the axis;
``R.H.H is R``. The identity :eq:`spaces-collapse-adjoint-is-pullback` is
the implementation, not a property the generic adjoint happens to have.

**The derivation, with the step each hypothesis enters.** Write the full
space as :math:`V = V_k\otimes V'`, a field as its slices
:math:`\psi = (\psi_n)_n` with :math:`\psi_n \in V'`, the retraction as
:math:`R\psi = \sum_n m_n\psi_n` with :math:`m = G_k u`, and let the full
metric be the **product**

.. math::

   G_V \;=\; G_k\otimes G_{V'} ,

the *product condition*. For every :math:`\psi\in V` and
:math:`\varphi\in V'`,

.. math::

   \langle R\psi,\varphi\rangle_{V'}
   &= \Bigl(\sum_n m_n\,\psi_n\Bigr)^{\!\mathsf T} G_{V'}\,\varphi
      \qquad \text{(the definition of } R\text{)} \\
   &= \sum_{n,n'} \psi_n^{\mathsf T}\,\bigl((G_k)_{nn'}\,u_{n'}\,G_{V'}\bigr)\,\varphi
      \qquad \text{(}m_n = \textstyle\sum_{n'} (G_k)_{nn'}u_{n'}\text{; linearity)} \\
   &= \psi^{\mathsf T}\,\bigl(G_k\otimes G_{V'}\bigr)
      \,(u\otimes\varphi)
      \qquad \text{(the double sum in block form)} \\
   &= \langle \psi,\ \pi^{*}\varphi\rangle_V
      \qquad \text{(the product condition; } \pi^{*}\varphi = u\otimes\varphi\text{)} .

The Hilbert adjoint is *defined* by
:math:`\langle R\psi,\varphi\rangle_{V'} = \langle\psi,
R^{\dagger}\varphi\rangle_V` for every :math:`\psi`, and the computation
above *exhibits* a vector with that property, :math:`u\otimes\varphi`.
Nothing in it inverts a metric: the covector :math:`m` was built as
:math:`G_k u`, so the :math:`G_k` the full metric carries is already the
one :math:`R` contracted with. The generic metric sandwich
:math:`\sharp_V\circ R^{\mathsf T}\circ\flat_{V'}` reaches the same
vector only by applying :math:`G_V^{+}` to :math:`m\otimes G_{V'}\varphi`,
which returns :math:`(G_k^{+}G_k u)\otimes\varphi`, equal to
:math:`u\otimes\varphi` exactly when :math:`u\in\operatorname{range}(G_k)`.
The adjoint is unique when :math:`G_V` is invertible, which under the
product condition means :math:`G_k` and :math:`G_{V'}` invertible (on a
nodal axis, every weight :math:`w_n\neq 0`); then
:math:`R^{\dagger} = \pi^{*}` exactly, with no metric in the result.

Each hypothesis has a failure, and each failure is a different object:

- **The full metric is not a product.** Then the third step fails, and
  no closed form exists. The space layer cannot spell this state on an
  axis-built space: :meth:`FunctionSpace.of_axes
  <orpheus.numerics.space.FunctionSpace.of_axes>` derives the metric
  from the axes, and a metric object is admitted only as a
  ``FactoredMetric`` with one entry per axis, a positioned form only on
  a measure-less axis's own block (:ref:`spaces-metric-object`). `[M]`
  2026-10-19, the S1 specification's census: on 112 (space, axis) pairs
  of 39 multi-axis spaces (the angular bulk, angular trial, bulk and
  moment spaces of six S\ :sub:`N` problems, seven heads ⊗ cells, two
  synthetic spaces), :math:`\lVert G - G_k\otimes G_{\rm rest}\rVert /
  \lVert G\rVert` is 0 on 110 and :math:`6.0\times10^{-17}` on the other
  2; gate G3 builds a full-shape ``DenseMetric`` on an axis-built space
  and requires the refusal (*"DERIVED from its axes"*). The non-product
  spaces of the tree (``FullFieldSpace``, bulk ⊕ trace, and the coupled
  system's space) are axes-less, and the by-label resolver refuses them
  (*"not axis-built"*).
- **The marginal metric is not the full metric with the block
  removed.** The marginal is :meth:`FunctionSpace.without_axis
  <orpheus.numerics.space.FunctionSpace.without_axis>`, which keeps every
  other axis's measure and positioned form (its inverse is
  :meth:`FunctionSpace.with_axis
  <orpheus.numerics.space.FunctionSpace.with_axis>`, the product widened
  by one axis, which the moment-tail composer
  :meth:`BulkField.compose_spatial_moments
  <orpheus.transport.fields._bases.BulkField.compose_spatial_moments>`
  uses); both read the overlay through one private reader,
  ``FunctionSpace._overlay_forms``. A marginal rebuilt with
  ``FunctionSpace.of_axes(*other_axes)`` would drop the forms and make
  :math:`G_{V'}` wrong whatever the adjoint is.
- **The axis block is singular.** A zero weight, or a frame-dressed head
  of deficient rank (``folded_product(4, 4)`` at :math:`L = 2` has a
  form of rank 5 on 15 slots), makes :math:`G_V` singular. The adjoint
  is then not unique: the pullback and the sandwich are two adjoints
  that both satisfy reciprocity, and they differ by a vector in
  :math:`\ker G_k\otimes V'`. On a moment head they agree, because its
  unit :math:`B e_0` lies in :math:`\operatorname{range}(B) =
  \operatorname{range}(B^{+})`. Reciprocity is blind to a write into the
  kernel, so the reciprocity rows cannot tell two adjoints apart there;
  what pins the agreement on the moment heads is G3's sandwich row, which
  compares the closed form with ``AdjointOperator(R)`` directly on the 19
  modal rows, the singular heads among them.

The closed form is therefore the Hilbert adjoint of every retraction the
mint can produce, and ``.H`` carries no guard.

A **form on the collapsed axis's own block is** :math:`G_k`. A dense form
positioned on a measure-less axis (the one form an axis measure cannot
carry) is the metric of that block, and the covector is computed with it,
:math:`m = G_k u`. Gate G3 carries such spaces (the S6 fixture ``D2``,
and the dense frame-dressed heads, rank 5 of 15, 8 of 45 and 23 of 45)
in its reciprocity population.

**Why the closed form, and not the sandwich.** Where every weight is
non-zero the two are one operator in exact arithmetic, so the choice is
a ruling about which spelling is the law, and three measurements decide
it (`[M]` 2026-10-02, on nodal axes, before the unit existed; on a
nodal axis :math:`u = \mathbf 1`, and gate G5 pins the arrows of today
``array_equal`` to the ones measured).

#. *Re-association.* On a diagonal metric the sandwich computes
   :math:`(w\,G\,\varphi)/(w\,G)` element by element, which rounds away
   from :math:`\varphi`: ``np.array_equal`` with the broadcast failed on
   **1 178 of 1 200** draws over six spaces (worst 3 ULP, on the
   linear-discontinuous trial space), and on 297 of 300 random
   S\ :sub:`N` detectors (two slab fixtures × Gauss–Legendre 5, 8, 16 ×
   50 draws; worst relative difference :math:`2.6\times10^{-16}`). Over
   a census of 1 891 suite mints the two agreed bit for bit on 46.
#. *A dense form.* The sandwich's :math:`\sharp` applies a dense form's
   pseudo-inverse, whose rounding is :math:`\kappa(G)\,\varepsilon`
   relative to the result (:math:`\kappa` the form's condition number):
   worst :math:`2.3\,\kappa\varepsilon` on ``D1`` (:math:`\kappa = 2.18`)
   and :math:`0.8\,\kappa\varepsilon` on ``D3`` (:math:`\kappa = 5.27`)
   over 200 draws each; gate S6.4 bounds it at
   :math:`8\,\kappa\varepsilon`. The pullback reads no metric and has no
   such term.
#. *A degenerate measure.* On a three-node axis with weights
   :math:`(0.5, 0, 1)` beside a weighted spatial axis, the pullback and
   the sandwich differ by 0.64 (maximum absolute, one draw) while both
   satisfy reciprocity (residuals 0.0 and :math:`1.1\times10^{-16}`):
   the adjoint is not unique there, and the closed form is the one the
   adjunction names.

The ruling (the orchestrator, 2026-10-02, design (B) of the step-6
specification, ``.claude/plans/reference_p1_spec.md`` §1.6): the
retraction gets its closed-form adjoint, because :math:`\pi^{*}` is the
operator algebra's own law for the adjoint of fibre integration under a
product measure, the :math:`(\pi_*, \pi^*)` adjunction, and a sandwich
that rounds differently from that law is the rounding of a
re-association, not a principled choice. Rejected, design (A): re-spell
the S\ :sub:`N` detector lift (below) as the sandwich and accept the ULP
change. The sandwich keeps a job, as the independent witness of the
closed form (gate S6.4).

.. dropdown:: First got wrong: the adjoint was the sandwich, and the mint dropped kept forms
   :color: muted

   Until step 6 of the first phase of the reference-solution campaign
   (#405, 2026-10-02) the pair was two operators, and ``R.H`` returned
   ``AdjointOperator`` wrapped around :math:`R`, the sandwich. The
   identity said that the sandwich *is* the plain broadcast; nothing in
   the tree computed the broadcast as the adjoint. The retraction's
   docstring had read ``R† = π*, the plain broadcast, [M]
   np.array_equal`` since 2026-08-24, and it was false: the
   "measurement" was of a different expression. The test-architect found
   it while asking whether the S\ :sub:`N` adjoint's hand-written
   broadcast was already ``R.H``; it was not, by up to 3 ULP.

   The same step found that the mint built the marginal as
   ``FunctionSpace.of_axes(*other_axes)``, dropping every positioned form
   of a kept axis. `[M]` 2026-10-02, the test-architect's mint census (a
   pytest plugin wrapping the mint while ``tests/gates/numerics``,
   ``transport``, ``sn``, ``diffusion`` and ``homogeneous`` ran: 9 113
   passed, 4 skipped, 45 xfailed; artefacts under
   ``scratch/reference_architecture/p1step6/ta/``): 1 891 mints, the
   product condition held on 1 891 of 1 891, and 0 of 1 891 full spaces
   carried a positioned form, so no shipped result was wrong through it.
   The hand-built witnesses ``D1`` (a dense form on a kept, measure-less
   axis) and ``D3`` (the composition route ``(angular ⊗ spatial) * dense
   head``) failed the product condition by 0.80 and 0.79 relative, and
   the broadcast violated reciprocity on them by 4.05 and 3.67. Step 6
   added ``without_axis`` and ``with_axis``.

   Step 6 also **refused** ``D2``, a dense form on the collapsed axis
   itself, as "two measures on one block": the retraction then
   contracted with the axis's weights while the space's metric on that
   block was the form. That refusal retired with the unit (S1,
   2026-10-19): the covector is :math:`m = G_k u`, so a form on the
   block is the measure the retraction integrates with, and ``D2`` is a
   member of the reciprocity population.

**One object, minted once.** The pullback is constructed inside the
retraction's own constructor, so the two refer to each other and
``R.H.H is R``. That satisfies the two-inductions clause
(:ref:`spaces-collapse-pair-frame`): :math:`R` is constructed only at
the mint, so all three arrows come from one mint. A second
``AxisPullbackOperator(R)`` is refused with a ``TypeError`` (*"minted by
its retraction, once"*), so the involution holds for every pullback
that can exist, not only for the one ``R.H`` returns.

.. list-table:: The three arrows of one axis marginal, three types
   :header-rows: 1
   :widths: 18 30 22 30

   * - Arrow
     - Type
     - What ``apply`` does
     - Law, with its measured rounding
   * - :math:`R = \pi_*`
     - :class:`~orpheus.numerics.axis_marginal.AxisRetractionOperator`
     - the pairing with :math:`m = G_k u` over the axis
     - :math:`R\,R^{\dagger} = \langle m, u\rangle\,\mathrm{id}`
   * - :math:`E`
     - :class:`~orpheus.numerics.axis_marginal.AxisSectionOperator`
     - divide by :math:`\langle m, u\rangle`, then place along :math:`u`
     - :math:`R\circ E = \mathrm{id}` to at most :math:`N_k/3` ULP
       (`[M]` 2026-10-19, 120 shipped rules × 20 draws)
   * - :math:`R^{\dagger} = \pi^{*}`
     - :class:`~orpheus.numerics.axis_marginal.AxisPullbackOperator`
     - :math:`u\otimes\varphi`
     - ``np.array_equal`` to the broadcast on nodal axes (`[M]`
       2026-10-02, 1 600 of 1 600 draws over eight spaces, gate S6.1);
       :math:`R^{\dagger} = \langle m, u\rangle\,E` to at most 1 ULP
       there (the cylinder; 0 ULP on the 7 spaces whose
       :math:`\Sigma w = 2`)

The population of the last column is the one gate S6.1 had when it was
measured: the synthetic
three-axis product, the S\ :sub:`N` bulk spaces of a three-interval
non-uniform slab at Gauss–Legendre 5, 8 and 16, that slab's
linear-discontinuous trial space, a two-interval cylinder's bulk space
at ``folded_product(4, 8)``, and ``D1`` and ``D3``; draws
``standard_normal`` scaled by :math:`10^{k}`, :math:`k\in[-6, 6]`, seed
11. Since S1 the gate's population also carries ``D2``, a dense form on
the collapsed axis, nine spaces in all. Gate G5 pins the same three arrows ``array_equal`` against the hand
einsum, the broadcast and ``broadcast(φ / E.mass)`` on 16 nodal
spaces.

**Three arrows, three types, and ERR-051.** The section and the adjoint
differ by exactly the mass, so a single undiscriminated verb ("lift",
"embed", "broadcast") would have had to choose one convention, and a
call site re-pointed from one use to the other would have rescaled a
quantity by the mass in silence. That is the failure class of ERR-051,
where an idempotency check asserted the wrong one of two identities a
factor :math:`4\pi` apart. Each arrow is its own type, so a call site
holds the object its use needs or an object of the wrong class; gate
``test_the_two_arrows_are_different_types`` asserts that ``R.H`` is an
``AxisPullbackOperator`` and not an ``AxisSectionOperator``.

.. _spaces-collapse-pair-two-lifts:

The two lifts of an angle-integrated quantity
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

On the angular axis the two arrows into the full space are the two
ways an angle-integrated quantity enters phase space, and the quantity's
**role** picks the arrow (the user's ruling of 2026-10-02, recorded in
``.claude/plans/reference_cache.md``, "P1 step 6 opened"):

- **A source rate** :math:`Q(\vec r, g)` enters the right-hand side
  through the section, :math:`q = E\,Q`, so that the rate is kept:
  :math:`R\,q = Q`. Every bulk isotropic source in S\ :sub:`N` already
  does this
  (:meth:`AngularSourceSink.from_isotropic
  <orpheus.transport.fields.angular_source_sink.AngularSourceSink.from_isotropic>`
  is :math:`E`'s ``apply``).
- **A detector response** :math:`\Sigma_d(\vec r, g)` is the functional
  :math:`\psi \mapsto \langle\Sigma_d, R\psi\rangle_{V'}`, the reaction
  rate the detector reads from the scalar flux. Its Riesz representative
  in :math:`V` is :math:`R^{\dagger}\Sigma_d = \pi^{*}\Sigma_d`, by the
  derivation above: the detector enters phase space through the
  adjoint, with no division by anything.

The two lifts differ by :math:`R\circ R^{\dagger} = \Sigma w`, the mass
of the measure the angular axis carries, and that mass is **never
typed**: it is :math:`4\pi` for a rule on the whole sphere
(``folded_product``, ``product``, ``level_symmetric``, ``lebedev``) and
2 for a one-dimensional Gauss–Legendre rule, whose ordinates are points
of the orbit space of :math:`\mu` alone, the sphere quotiented by the
rotations and reflections about the polar axis
(:ref:`manifold-orbit-space`). A hand-written :math:`1/4\pi` is
therefore wrong on the slab and the sphere, where :math:`\Sigma w = 2`,
and a detector lifted by the source's arrow reads :math:`1/\Sigma w`
of its true response. The user's question that settled it (2026-10-02):
*"integrated over all directions is measure related, and divided by 4π
is also measure related … a defective operator"*. Both are the measure,
the operators already existed, and the convention became two typed
arrows.

The continuous counterpart, on :math:`S^2` with
:math:`\mathrm d\Omega = \mathrm d\mu\,\mathrm d\varphi`, is Branch 1's
:mod:`orpheus.derivations.common.angular_measure`, which derives the
mass :math:`4\pi` by integration and proves both lifts in SymPy
(:ref:`structured-geometry-two-lifts-branch-1`).

**The first consumer: the S**\ :sub:`N` **adjoint entry.**
:func:`~orpheus.sn.solver.solve_sn_adjoint_fixed_source` lifted its
detector with a hand-written ``np.broadcast_to(sigma_d[None], …)``:
:math:`R^{\dagger}`'s content, but not :math:`R^{\dagger}`. Since step 6
it reads ``AxisRetractionOperator.of(problem.angular_bulk_space, "angular").H.apply(sigma_d)``
(:ref:`sn-adjoint-dual-lift`). The change was bit-exact: `[M]` the qa
review's capture before and after, 36 of 36 arrays ``np.array_equal``
(12 scalar adjoint fluxes and 24 angular-flux leaves; the
three-interval vacuum | reflective slab, 2 and 4 groups, diamond and
linear-discontinuous, Gauss–Legendre 6, 3 random sparse detectors
each; ``scratch/reference_architecture/p1step6/qa/p4_sn_bitexact.py``).
Bit-exactness is a consequence of the closed form: the sandwich
re-spelling would have moved 297 of 300 detectors by up to 3 ULP.

.. warning::

   **The adjoint has a second spelling still, and a hand broadcast
   survives.** A composite operator that holds a retraction (a sum or
   a product with :math:`R` as a leaf) daggers through the generic
   sandwich, not through the leaf's pullback, so :math:`R^{\dagger}`
   has two spellings in the tree (#558). And
   ``AngularSourceSink.__add__`` with a ``ScalarSourceSink`` partner
   adds ``other.values[None]``, a hand-written :math:`\pi^{*}`
   (the dunder note in :doc:`/theory/foundations/operator_algebra`
   records that it is the pullback and not the section). Neither is
   wrong in value; both are the concept spelled outside its type.

**The gates.** ``tests/gates/numerics/test_retraction_adjoint_is_the_pullback.py``
(``foundation``; the specification is
``.claude/plans/reference_p1_spec.md`` §1.6, rows S6.1 to S6.4):

.. list-table::
   :header-rows: 1
   :widths: 22 78

   * - Row
     - What it pins, and the mutation that reddens it
   * - S6.1
     - ``R.H`` is an ``AxisPullbackOperator``, ``R.H.H is R``, the two
       spaces swap, and ``R.H.apply`` is ``np.array_equal`` to the
       broadcast, over nine spaces (the eight above and ``D2``) × 200
       draws. Red when ``.adjoint()`` returns the sandwich (8 of 8 in
       step 6's battery), and the involution
       leg red when the pullback's adjoint is a fresh retraction. A
       companion row asserts that a second pullback cannot be minted.
   * - S6.2
     - the product condition, with the weights read from the raw
       ``Axis`` objects and never from :math:`R`, and reciprocity at
       ``rtol=1e-13``. Red when the pullback is scaled by 2, and on
       ``D1`` and ``D3`` when the marginal drops the kept forms.
   * - S6.3
     - the marginal of ``D1`` carries ``D1``'s dense form. Red when the
       marginal drops the forms. Until S1 the row also required the mint
       to refuse ``D2``; that refusal retired with the unit, and ``D2``
       joined the S6.1, S6.2 and S6.4 population.
   * - S6.4
     - the X4 witness: the closed form against an explicitly built
       sandwich, ``nulp=4`` on a diagonal metric and
       :math:`8\,\kappa\varepsilon` on a dense form. Red when
       :math:`\pi^{*}` is replaced by :math:`E` or by
       :math:`R^{\mathsf T}`. Designed green, and said so in its
       docstring: :math:`\pi^{*} := \Sigma w\,E` is the same operator to
       0 or 1 ULP, so this row cannot tell them apart (S6.1 can).
   * - verbs
     - ``without_axis`` and ``with_axis`` keep every other axis's form,
       with ``of_axes`` over the same axes as the control that drops it.

At step 6 (2026-10-02) two rows of
``tests/gates/numerics/test_axis_marginal.py`` changed meaning
(``retirement-audit`` D.14 and D.15). The G6.3 pairing row is
**promoted**: it compared the sandwich, whose reciprocity is its own
definition, so it could not fail on any :math:`R`; it now compares the
closed form, which reciprocates only because the mint admits product
metrics. The G6.4 row's rounding sentence now describes the broadcast
against :math:`\Sigma w\,(\varphi/\Sigma w)`. Step 6's mutation battery
(``scratch/reference_architecture/p1step6/impl/mut_6a.py``, in process):
the sandwich as the adjoint reddened S6.1 on 8 of 8 spaces and the
two-arrows row; the pullback × 2 reddened S6.1, S6.2 and S6.4 on 8 of 8
each and the G6.3 and G6.4 rows; the refusal removed reddened S6.3; the
marginal dropping kept forms reddened S6.2 and S6.4 on ``D1`` and ``D3``
and S6.3's positive leg.


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

The marginal is induced — the axis carries the stage-2 datum
------------------------------------------------------------

The kernels above could be written by hand: an einsum and a broadcast,
six lines each. What keeps them from being a hand-written twin of the
frame machinery (:doc:`/theory/foundations/frame`) is where their two
retained numbers come from: the unit :math:`u` and the metric block
:math:`G_k`, which the axis and the space already carry, and nothing the
marginal chooses for itself.

The rule is the **stage-2 generator discipline** (user, 2026-08-24):

   A stage-2 generator induces structure on both the space and the
   operator, and the two inductions must be minted together, at one
   site. Frame: induces the harmonic head axis's metric (space side), mints
   Analysis/Synthesis (operator side) — consistency is the tightness
   gate. Scheme: induces the trace descriptor and basis kind (space
   side), mints the closure (operator side) — consistency is one
   closure serving both apply and solve, which is ERR-026's structural
   closing. Mesh and Quadrature are the degenerate cases (space side
   only). Forgetting = retaining the induced parts; accessors are
   provenance.

The unit is one of the induced parts, so it is minted where the axis is
minted, by whatever mints it (ruling 37 of the numerics-untangling plan):

- a **measure** mints a nodal axis, whose unit is all ones by the kind
  (the axis fills it in at construction);
- a **basis** spans a continuum head, and its layout carries the
  coefficients of :math:`1` in that basis, :math:`e_0` for the harmonic
  and Legendre layouts (``FromBasis(basis).kind`` is
  ``Modal(basis.layout.unit)``; :meth:`FunctionSpace.of_basis
  <orpheus.numerics.space.FunctionSpace.of_basis>`); a partition's cell
  layout is ``NODAL``, unit ones, and a mode basis declares its unit,
  zero for the loss-kernel gauge (rulings 86 and 87);
- a **frame**, when it dresses a head (``frame_faces._dressed``), writes
  the analysis of :math:`1` from its own data, :math:`u =
  \Phi^{\mathsf T} w`, as it writes the head's weights or its positioned
  form: the trial basis's ``analyze`` applied to the all-ones ordinate
  field with the frame's table and weights;
- the **LD scheme** writes :math:`e_0` when it mints the spatial-moment
  axis (``DiscretizationSchemeBase.moment_axis``).

The axis keeps the induced datum as content, and the object that minted
it stays provenance (:attr:`Axis.origin <orpheus.numerics.axis.Axis.origin>`).
The marginal reads the origin in one place only: a ``FromBasis`` axis
whose layout is a :class:`~orpheus.numerics.coefficient_layout.CellLayout`
is a partition's raw coefficient head, which is **measure-free**, so the
marginal over it refuses with a ``TypeError`` naming the frame that would
supply the fine measure (ruling 88: a cell carries volume only against
the measure a frame binds, :math:`4\pi r^2\Delta r` on a spherical shell,
the fine-group count on a coarse group). The frame-dressed head,
``Modal(masses)`` (ruling 83), integrates: `[M]` 2026-10-22, the raw head
on edges :math:`[0, 1, 3]` refuses and the dressed head maps the constant
function to :math:`3` (gate ``tests/gates/numerics/test_heads_on_the_axis.py``
H2d).

Three consequences, each visible in the code
(:mod:`orpheus.numerics.axis_marginal`):

**All three arrows at one site.** The private mint ``_mint`` constructs
the retraction (and, inside it, the pullback) and the section together
and returns them as one record, memoised per (space, label). There is no
path that produces one without the others, so they cannot disagree about
the axis, the dims or the marginal space: `[M]` the arrows share ONE
marginal-space instance (``R.codomain is E.domain``), gated.

**Forgetting means copying the induced parts out.** The operators retain
the bound spaces, the ndarray dims the axis occupies, the axis's shape,
one flat vector over the axis (the covector :math:`m` on the retraction,
the unit :math:`u` on the pullback and the section) and, on the section,
the scalar mass, and nothing else: no axis, no frame, no generator.

**Consistency is a gate, not an instance.** Nothing at runtime forces
the minted kernels to agree with any frame. On a nodal axis the
correspondence the corpus relies on is with the rank-one frame (below),
and the **tightness gate** supplies it by building that frame
independently in the test and pinning :math:`R` against its analysis
content, :math:`R^{\mathsf T}` against ``analyze_transpose``, and
:math:`E` against reconstruction composed with the inverse Gram, each
``np.array_equal``. The operator kernels are einsums over the covector
and the unit, and the frame path runs the basis's table einsums, so
these are two different float programs and the agreement is a
comparison of two implementations.

**The rank-one frame, as a mathematical fact.** On a nodal axis, the
single-region indicator basis covering the index set
:math:`\{0,\dots,n-1\}`, bound to the axis's measure in a Galerkin frame,
has one coefficient, a table of ones, and a :math:`1\times1` Gram
(:math:`G_{jk} = \sum_n w_n \phi_j(x_n)\phi_k(x_n)`) equal to the axis's
mass:

.. math::
   :label: spaces-collapse-rank-one-gram

   G \;=\; \bigl[\textstyle\sum_n w_n\bigr] \;=\; [\,\langle m, u\rangle\,],
   \qquad
   E \;=\; R_{\rm frame}\circ G^{-1},

.. (vv-status rationale) Structural: the rank-one instance of the
   Parseval theorem derived at frame-parseval-metric (the frame's
   codomain metric is the INVERSE discrete Gram), specialised to a
   single-region indicator basis on a nodal axis, where the Gram is 1x1
   and its entry is the mass. Not a solver claim, and no production code
   builds this frame: the marginal is built from the unit (u = 1, m = w
   on a nodal axis). The verifiable content is the tightness row of
   ``tests/gates/numerics/test_axis_marginal.py``
   (``TestFrameInduction``), which builds the literal frame in the test
   and pins the three arrows against its faces.
.. vv-status: spaces-collapse-rank-one-gram documented

so the section's divisor is the Parseval metric of
:ref:`frame-parseval-metric` at :math:`K = 1`, and "divide by
:math:`4\pi`" is the inverse Gram of a frame, not a convention this code
chose. The correspondence is exact in real arithmetic; in floating point
the frame's Gram entry is a three-operand einsum and the section's
divisor is the correctly rounded sum, which differ on 68 of 120 shipped
rules (:ref:`spaces-collapse-pair-evidence`), so the tightness row is
pinned on a fixture whose mass is exactly 2.0.

.. note::

   **A latent generalization, recorded and not built.** On a nodal axis
   the marginal is the pairing with :math:`m = w`, a uniform average
   against the axis's measure. A *profiled* collapse (the
   Petrov-Galerkin test side of a :math:`\chi`-class emission collapse,
   where an axis is averaged against a spectrum rather than uniformly)
   is the same construction with a different test weight. It is built
   when its consumer lands, not before; recording it here is what stops
   a future session minting a second, parallel mechanism for it.


.. _spaces-collapse-pair-clause-gate:

The clause gate: which axes admit, and why energy refuses
-----------------------------------------------------------

The mint refuses an axis the doctrine says must persist. This is the
one place in the tree where the collapse doctrine is **enforced** rather
than described, so it is worth reading the admission table against
:ref:`spaces-collapse-doctrine-standing` clause by clause. Every refusal
is raised at one site, the mint in :mod:`orpheus.numerics.axis_marginal`,
or in the space's by-label resolver it calls first.

.. list-table:: Admission at the mint
   :header-rows: 1
   :widths: 26 12 62

   * - Axis
     - Clause
     - Verdict, and why
   * - **angle** (an untyped ``NODAL`` ``Axis``)
     - 3
     - **ADMIT.** Whole-domain integration over a compact canonical
       orbit: the total is universal, nothing problem-specific survives
       on the axis, so the drop-form marginal is exactly right and the
       re-broadcast convention lives on the arrows :math:`E` /
       :math:`\pi^{*}`.
   * - a ``MODAL`` head (Legendre, harmonic, a frame-dressed moment head,
       the LD spatial-moment axis)
     - 3
     - **ADMIT.** The marginal is the integral of the represented
       function, :eq:`spaces-axis-marginal-integral`, through the axis's
       unit; on a moment head it is the integral over the ordinates the
       moments were analysed from.
   * - an axis whose block carries a positioned form (a dense form on a
       measure-less axis)
     - —
     - **ADMIT.** The form is :math:`G_k`, and the covector is
       :math:`m = G_k u`.
   * - a typed :class:`~orpheus.numerics.axis.EnergyAxis`
     - 2
     - **REFUSE** (``TypeError``, *"condensation"*), with a pointer.
       Partition-integration of an :math:`L^1` class: the energy axis
       PERSISTS at its one-cell member, because
       :math:`\langle\bar\sigma,\varphi\rangle` consumes the partition. A
       drop-form marginal here would be a second mechanism for a
       collapse the tree already implements as **condensation**.
   * - an untyped generic axis, whatever its label
     - 3
     - **ADMIT.** The gate reads the axis's TYPE, never its label
       string.
   * - the only axis of a single-axis space
     - —
     - **REFUSE** (``ValueError``, *"bare scalar"*). Its marginal would
       be a bare scalar, which is not a
       :class:`~orpheus.numerics.space.FunctionSpace`. Contract with the
       space's inner product instead.
   * - a space with ``axes is None``
     - —
     - **REFUSE** (``TypeError``, *"not axis-built"*, in the by-label
       resolver ``FunctionSpace._axis_index``). An axes-less space (a
       hand-named space, or a product with an axes-less factor) has no
       named factors to marginalise over.
   * - a label naming zero or several axes
     - —
     - **REFUSE** (``ValueError``, *"names 0 axes"*, in the same
       resolver), naming the space's axis inventory.
   * - the **section** over an axis of zero mass
     - —
     - **REFUSE** at :meth:`AxisSectionOperator.of
       <orpheus.numerics.axis_marginal.AxisSectionOperator.of>`
       (``ValueError``, *"zero total weight"*); the retraction over the
       same axis is legal (:ref:`spaces-collapse-pair-two-arrows`).

The energy row is the load-bearing one, and it is the doctrine's
Version-1 refutation cashed out in code
(:ref:`spaces-collapse-version-1`). Energy is collapsed *by
integration*, so a compactness-flavoured reading puts it on the "drop"
side; the shipped one-group layout is :math:`(1, *\mathrm{spatial})`
and keeps its axis, because :math:`\bar\sigma` is defined only relative
to an interval and a weighting spectrum. The machinery that performs an
energy collapse correctly therefore already exists, it is
:meth:`EnergyGrid.overlap_to
<orpheus.data.energy_grid.EnergyGrid.overlap_to>` and the
Petrov-Galerkin condensation frames of
:ref:`sn-energy-condensation`, and the refusal message names it, so
a caller who reaches for the wrong tool is handed the right one rather
than a wrong answer. The energy clause is checked before the
single-axis clause, so a lone energy axis is refused toward
condensation and not as a bare scalar (gate G6, the wiring-order row).

The clause gate is not hypothetical on a shipped carrier. `[M]` 2026-10-19 at
``7baf444c``, on the :math:`S_N` fixture of :ref:`spaces-collapse-pair-evidence`
(a five-cell slab, ``gauss_legendre(4)``, 2 groups), whose
``angular_bulk_space`` axes are ``(Axis, EnergyAxis, Axis)``, and on its
moment space at :math:`L = 1`, whose head is a frame-dressed Legendre
head:

.. code-block:: text

   AxisRetractionOperator.of(sn.angular_bulk_space, "angular")  -> OK, codomain (2, 5)
   AxisRetractionOperator.of(sn.angular_bulk_space, "spatial")  -> OK, codomain (4, 2)
   AxisRetractionOperator.of(sn.angular_bulk_space, "energy")   -> TypeError: axis
       marginal: axis 'energy' is a typed EnergyAxis, which PERSISTS at its
       one-cell member (collapse doctrine clause 2 ...). The energy collapse is
       condensation: use EnergyGrid.overlap_to / the Petrov-Galerkin
       condensation frame, not a drop-form marginal.
   AxisRetractionOperator.of(sn.moment_space(1), "legendre")    -> OK, codomain (2, 5);
       the section's mass is 2.0

Three of the carrier's factors marginalise, the modal head among them;
the one the doctrine says must persist refuses, by TYPE, with the
successor named in the message.

.. warning::

   **The gate reads the TYPE, and today only energy has one.** A
   generic ``Axis(label="energy", ...)``, a synthetic test factor, is
   ADMITTED, and correctly so: refusing on the label string would be
   stringly-typed dispatch, and it would refuse a legitimate synthetic
   fixture that carries none of energy's physics. The consequence is
   that the clause gate is **structural for energy and permissive for
   everything else** until the typed spatial and quadrature axes land
   (CS2); at that point the verdict becomes axis-family polymorphism and
   each family answers for itself. Until then, do not read "the mint
   admitted it" as "the doctrine says it drops".


.. _spaces-collapse-pair-refuted:

What was tried, and what refuted it
------------------------------------

Each entry is a design the marginal does not have, with the structural
reason it does not.

**The hand-derived pair (S6.0, 2026-08-24), superseded within a day.**
The first implementation minted the two operators from the space with
hand-spelled kernels and a hand-chosen divisor (``weights.sum()``). It
was correct and it was a twin of the rank-one frame
(:eq:`spaces-collapse-rank-one-gram`), which is Cardinal Rule 2's stop
condition. The re-carve (S6.0b) kept the operator shells verbatim and
changed only where the two retained numbers came from: a single-region
indicator frame was built over the axis's index set at the mint, read
for its measure weights and its :math:`1\times1` Gram entry, and
discarded.

**The rank-one indicator frame as the mint (S6.0b), retired by the unit
(S1, 2026-10-19).** The frame read two numbers, the weights and their
sum, and the axis now supplies both directly, as the metric block and
the unit. Three reasons retired it. *(i)* It could not express a modal
marginal: built over the index set with the axis's weights, its analysis
contracts a modal axis's coefficients with its weights, which is not an
integral of the represented function, so the mint refused every
``MODAL`` axis ("the modal average is the coefficient at the average
slot; slice it"). With the unit the modal marginal is
:eq:`spaces-axis-marginal-integral`. *(ii)* It could not read a form on
the collapsed axis's block, so the mint refused that state as "two
measures on one block"; with :math:`m = G_k u` the form is the measure.
*(iii)* Its divisor, the Gram entry, a three-operand einsum, was the
least accurate spelling of the mass (up to 39 ULP off the correctly
rounded sum, :ref:`spaces-collapse-pair-evidence`). The stage-2
discipline survives the retirement, because the axis, the forgetful
image of its generator, now carries the induced datum
(:ref:`spaces-collapse-pair-frame`).

**The unit read through the generator basis, refuted before it was
built (2026-10-19).** A verb ``Basis.constant_coefficients()`` would have
returned :math:`e_0` per family; on a frame-dressed head the right unit
is :math:`B e_0`, and :math:`e_0` integrates to :math:`1/\Sigma w` of the
answer, or to no multiple of it (the dropdown in
:ref:`spaces-axis-unit`). The variance of a head's coordinates is decided
where the head is dressed, so the unit is minted there.

**Keeping the divisor's old bits, rejected (ruling 38).** Spelling the
mass as the retired frame's einsum would have kept every downstream
array bit for bit. It was rejected as an API spelling whose only purpose
is to reproduce a legacy reduction tree (``vv-principles``, bit-identity
versus principled equivalence): the correctly rounded mass is a named
quantity, its reference is platform-free, and the drift is the reduction
depth times the ULP. The moved pins were re-baselined in the S1 commit
(the section divisor on 68 of 120 shipped rules, up to 39 ULP on
``lebedev(29)``; one S\ :sub:`N` regression case,
``2d_2g_LS4_dd_8x4_het_si``, moved its :math:`k_{\rm eff}` by 1 ULP,
inside its :math:`10^{-11}` tolerance).

**An** ``Axis`` **→ measure accessor, refused (S6.0b).** The indicator
mint needed a :class:`~orpheus.numerics.measure.DiscreteMeasure` built
from the axis, and the obvious move was a public accessor on
:class:`~orpheus.numerics.axis.Axis` producing one. It was rejected under
the same ruling that governs the frame: *accessors are provenance*. An
accessor would have made a generator reachable from the axis forever,
and a pre-CS5 axis had dropped the nodes, so the only node set it could
synthesise was the index set (``nodes = arange(n)``): a synthetic
measure under a name readers would take for the generating one. The
mint built that measure with a local helper instead. Since S1 the mint
builds no measure at all; it reads the unit and the metric block. The
refusal stands, and the :attr:`~orpheus.numerics.axis.Axis.origin`
slot does not overturn it, because the origin points **generator → axis** and
manufactures nothing: it records the real object at the one moment it
is in scope (:ref:`spaces-generator-not-a-reverse-accessor`).

**Caching the marginal on the** ``Axis``, **refused.** The marginal is
not a function of the axis alone: its domain is the whole product and
its codomain the product minus that factor, so two spaces sharing an
axis have *different* marginals. The cache belongs to the space instance
(the private slot ``"_axis_marginal"``, one entry per axis label).

**Retaining the frame on the operators, refused (S6.0b).** Keeping the
indicator frame would have made consistency automatic instead of gated.
It retained a whole generator (basis, table, measure, two spaces) on
every collapse operator to secure a property a gate already secured, and
it re-opened the question of whether two operators built from *equal*
frames are the same operator. Since S1 no frame exists to retain.


.. _spaces-collapse-pair-evidence:

Numerical evidence
-------------------

`[M]` 2026-08-24, measured on the nodal pair of that day. On a nodal axis
:math:`u = \mathbf 1` and :math:`m = w`, and gate G5 pins the arrows of
today ``array_equal`` to the hand einsum and the broadcast those
measurements used, so the rows below describe today's arrows; the
divisor table and the unit table further down were measured on
2026-10-19 at ``7baf444c``. **The construction,
so the tables regenerate from this page.** The *synthetic* fixture is
the three-axis product ``angular(4, w=[0.3, 0.7, 0.5, 0.5]) ⊗
energy(2, counting) ⊗ spatial(5, V=[0.2, 0.3, 0.4, 0.7, 1.4])`` built
with :meth:`FunctionSpace.of_axes
<orpheus.numerics.space.FunctionSpace.of_axes>`; the weights are
non-uniform on purpose so that no cancellation flatters a law. The
:math:`S_N` fixture is the shipped carrier —
``SNProblem(Mesh1D(edges=[0, 0.2, 0.5, 0.9, 1.6, 3.0], cartesian, vacuum),
Quadrature.gauss_legendre(4), 2 groups)`` — and the operators come from
``AxisRetractionOperator.of(sn.angular_bulk_space, "angular")`` /
``AxisSectionOperator.of(sn.angular_bulk_space, "angular")``.
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

.. note::

   **Read the two** :math:`R^{\dagger}` **rows with their date.** They
   were measured on 2026-08-24, when ``R.H`` was the generic metric
   sandwich, so they record the sandwich's rounding. Since #405 P1 step 6
   (2026-10-02) ``R.H`` is the closed form
   :class:`~orpheus.numerics.axis_marginal.AxisPullbackOperator`: the row
   :math:`R^{\dagger} = \pi^{*}` is ``np.array_equal`` by construction
   (`[M]` 1 600 of 1 600 draws over eight spaces, gate S6.1), and
   :math:`R^{\dagger} = (\Sigma w)\,E` is at most 1 ULP
   (:ref:`spaces-collapse-pair-pullback`). The other rows are
   unchanged by that step, which touched neither :math:`R` nor :math:`E`.

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

**The divisor is the correctly rounded mass, and every other spelling
rounds differently somewhere.** `[M]` 2026-10-19, re-run at
``7baf444c``
(``scratch/numerics_untangling/packages/S1_doc_probes/s1_numbers.py``, run
from outside the repository as ``PYTHONPATH=<repo> <repo>/.venv/bin/python
-O <probe>``), over the 120 shipped rules of
gate G5-mass: Gauss–Legendre 2 to 64, level-symmetric 2 to 16,
``product`` and ``folded_product`` over :math:`(2, 4, 6, 8)\times(4, 8,
12, 16)`, and 17 Lebedev orders from 3 to 41. Each spelling of
:math:`\Sigma w` is compared with ``math.fsum(w)``:

.. list-table:: Spellings of the mass against the correctly rounded sum (120 rules)
   :header-rows: 1
   :widths: 46 18 36

   * - spelling
     - rules where it differs
     - worst, in ULP of the sum
   * - ``einsum("n,nR,nS->RS", w, 1, 1)[0, 0]`` (the retired frame's
       :math:`1\times1` Gram entry, the divisor until S1)
     - 68
     - 39 (``lebedev(29)``)
   * - ``einsum("n,n->", w, 1)``
     - 41
     - 11 (``lebedev(27)``)
   * - ``w @ 1``
     - 21
     - 4 (``lebedev(35)``)
   * - ``w.sum()``
     - 29
     - 3 (``lebedev(25)``)
   * - ``AxisSectionOperator.of(V, "angular").mass``
     - 0
     - 0 (``array_equal`` on 120 of 120)

So no spelling but the correctly rounded one is the same number on every
rule, and the one the pair shipped with was the furthest from it. The
lesson for a gate: pin the divisor against ``math.fsum`` of the stored
weights, a platform-free reference, never against another reduction of
the same weights. On the rules the S\ :sub:`N` regression corpus uses,
the change of spelling moved the divisor on three: by 2 ULP on
``level_symmetric(4)``, 1 ULP on ``gauss_legendre(16)`` and 4 ULP on
``product(4, 8)`` (the S1 specification, §7; one regression case's
:math:`k_{\rm eff}` moved by 1 ULP).

**The unit on the frame-dressed heads.** `[M]` 2026-10-19 at ``7baf444c``
(``scratch/numerics_untangling/packages/S1_doc_probes/s1_heads.py``, the same invocation): each head is
``frame_faces.basis_space(frame_faces.angular_frame(quad, L))`` ⊗ a counting axis of
20 slots, :math:`\psi` is ``default_rng(7).standard_normal((N, 20))``,
the moments are the frame's own analysis :math:`\Phi^{\mathsf T}W\psi`,
and the reference is :math:`\sum_n w_n \psi_n`. The column "with
:math:`e_0`" is the marginal the refuted design would have computed,
:math:`(G_k e_0)^{\mathsf T}\varphi` divided by the reference, over the
20 columns.

.. list-table:: The unit :math:`u = B e_0` against :math:`e_0`, on moment heads
   :header-rows: 1
   :widths: 24 8 16 18 18 16

   * - rule, order
     - Gram
     - :math:`\max\lvert u - \Sigma w\,e_0\rvert`
     - :math:`\max\lvert R\varphi - \int\psi\rvert / \max\lvert\int\psi\rvert`
     - ratio with :math:`e_0`
     - mass :math:`-\ \Sigma w`
   * - ``gauss_legendre(8)``, :math:`L = 3`
     - diagonal
     - :math:`2.8\times10^{-17}`
     - :math:`2.1\times10^{-17}`
     - 0.5000
     - 0 ULP
   * - ``lebedev(7)``, :math:`L = 3`
     - diagonal
     - :math:`1.8\times10^{-15}`
     - :math:`8.7\times10^{-17}`
     - 0.0796
     - +1 ULP
   * - ``level_symmetric(4)``, :math:`L = 3`
     - dense
     - :math:`3.6\times10^{-15}`
     - :math:`5.6\times10^{-17}`
     - 0.0796
     - −2 ULP
   * - ``folded_product(4, 4)``, :math:`L = 2`
     - dense
     - :math:`6.7\times10^{-16}`
     - :math:`7.0\times10^{-17}`
     - 0.0796
     - 0 ULP
   * - ``folded_product(4, 4)``, :math:`L = 4`
     - dense
     - **3.67**
     - :math:`2.5\times10^{-16}`
     - :math:`-0.11` to :math:`0.18`
     - +1 ULP
   * - ``product(4, 8)``, :math:`L = 4`
     - dense
     - :math:`7.1\times10^{-15}`
     - :math:`7.7\times10^{-17}`
     - 0.0796
     - −4 ULP
   * - ``lebedev(13)``, :math:`L = 2, 4, 6` (8 negative weights)
     - diagonal
     - :math:`\le 2.9\times10^{-15}`
     - :math:`\le 2.0\times10^{-15}`
     - 0.0796
     - 0 ULP
   * - ``lebedev(25)``, :math:`L = 2, 4, 6` (6 negative weights)
     - diagonal
     - :math:`2.8\times10^{-14}`
     - :math:`\le 1.9\times10^{-15}`
     - 0.0796
     - −16 ULP
   * - ``lebedev(27)``, :math:`L = 2, 4, 6` (18 negative weights)
     - diagonal
     - :math:`5.7\times10^{-14}`
     - :math:`\le 2.4\times10^{-16}`
     - 0.0796
     - −32 ULP

Three readings. *(i)* With the unit the head mints, the marginal of the
moments is the ordinate integral at round-off on every row, the
signed-weight rules included. *(ii)* With :math:`e_0` it is
:math:`1/\Sigma w` of the integral (0.5 on the line, :math:`1/(4\pi) =
0.0796` on the sphere) wherever the rule integrates the head's live
harmonics exactly, and no multiple of it where it does not
(``folded_product(4, 4)`` at :math:`L = 4`, the row whose unit differs
from :math:`\Sigma w\,e_0` at O(1)). *(iii)* The mass on a moment head is
:math:`\Sigma w` to a few ULP, not exactly: the unit is the frame's
analysis of :math:`\mathbf 1`, an einsum, so :math:`u_0` already carries
that reduction's rounding (up to 32 ULP on ``lebedev(27)``); gate G2c
bounds the mass at :math:`8\varepsilon` relative on its five heads.

**The section is not bit-exact against the retraction.** `[M]` 2026-10-19 at
``7baf444c``, the 120 rules × 20 draws ``standard_normal(3)`` scaled by
:math:`10^{k}`, :math:`k\in[-6, 6]`, seed 7: the worst ULP of
:math:`R(E\varphi)` against :math:`\varphi` is at most :math:`N_k/3` on
every rule. The ratio is attained on ``lebedev(3)`` (2 ULP,
:math:`N_k = 6`), the largest count is 67 ULP on ``lebedev(41)``
(:math:`N_k = 590`), and 2 of the 120 rules are bit-exact on all 20
draws. The bound grows
with the reduction depth :math:`N_k`, which is why gate G4 is written at
``nulp = max(4, N_k)`` and not as ``array_equal``.

**The section is the harmonic frame's isotropic column — at the
measure level.** A natural question is whether the angular marginal
duplicates the spherical-harmonic frame's :math:`\ell = 0` channel.
`[M]` 2026-10-18 (``scratch/metric_docstring/probe_iso_column.py``, run
as ``PYTHONPATH=. .venv/bin/python -O <probe>``): on the :math:`S_N`
carrier above, on its spherical twin (the same edges and the same
``gauss_legendre(4)``), and on the slab carrier with dense-Gram
sphere rules, build
``HarmonicFrame.from_galerkin(frame_faces.angular_frame(sn.quad, L))``, feed a moment
field that is zero except at the :math:`(0,0)` slot, and compare both
faces with :math:`E\varphi`, :math:`E` the angular section
``AxisSectionOperator.of(sn.angular_bulk_space, "angular")``. Each entry is the worst
relative deviation over 20 draws ``default_rng(1000+k)``; "row-0
coupling" is :math:`\max_{j>0}|G_{0j}|/G_{00}`. ``SNProblem`` on the
spherical mesh refuses the 3-D rules (``ValueError``: an ordinate's
``tau_raw`` lies outside :math:`[0, 1]`), so the dense rows run on the
slab mesh; the frame reads only
``sn.quad``, so the mesh does not enter.

.. list-table:: The isotropic column against the section
   :header-rows: 1
   :widths: 24 14 14 24 24

   * - frame
     - Gram verdict, max off-diagonal
     - row-0 coupling
     - :math:`\text{face}^{\dagger}(e_0\varphi)` vs :math:`E\varphi`
     - :math:`\text{reconstruction}(e_0\varphi)/W` vs :math:`E\varphi`
   * - slab ``gauss_legendre(4)``, :math:`L=1`
     - DIAGONAL, :math:`1.4\times10^{-17}`
     - :math:`6.8\times10^{-18}`
     - :math:`2.1\times10^{-16}`
     - **0.0** (``array_equal`` 20 of 20)
   * - slab ``gauss_legendre(4)``, :math:`L=2`
     - DIAGONAL, :math:`5.6\times10^{-17}`
     - :math:`2.3\times10^{-17}`
     - :math:`2.1\times10^{-16}`
     - **0.0** (20 of 20)
   * - sphere ``gauss_legendre(4)``, :math:`L=1`
     - DIAGONAL, :math:`1.4\times10^{-17}`
     - :math:`6.8\times10^{-18}`
     - :math:`2.2\times10^{-16}`
     - **0.0** (20 of 20)
   * - sphere ``gauss_legendre(4)``, :math:`L=2`
     - DIAGONAL, :math:`5.6\times10^{-17}`
     - :math:`2.3\times10^{-17}`
     - :math:`2.2\times10^{-16}`
     - **0.0** (20 of 20)
   * - ``product(4, 4)``, :math:`L=2`
     - DENSE, :math:`1.088`
     - :math:`2.2\times10^{-17}`
     - :math:`1.7\times10^{-16}`
     - **0.0** (20 of 20)
   * - ``folded_product(4, 4)``, :math:`L=2`
     - DENSE, :math:`1.088`
     - :math:`5.3\times10^{-17}`
     - :math:`3.5\times10^{-16}`
     - **0.0** (20 of 20)
   * - ``folded_product(2, 4)``, :math:`L=3`
     - DENSE, :math:`2.793`
     - :math:`1.1\times10^{-16}`
     - :math:`2.1\times10^{-16}`
     - **0.0** (20 of 20)
   * - ``level_symmetric(4)``, :math:`L=3`
     - DENSE, :math:`0.193`
     - :math:`1.2\times10^{-16}`
     - :math:`9.0\times10^{-16}`
     - :math:`4.2\times10^{-16}` (0 of 20)
   * - ``product(4, 8)``, :math:`L=4`
     - DENSE, :math:`1.112`
     - :math:`4.7\times10^{-17}`
     - :math:`5.3\times10^{-16}`
     - :math:`6.0\times10^{-16}` (0 of 20)
   * - ``folded_product(4, 4)``, :math:`L=4`
     - DENSE, :math:`3.665`
     - :math:`0.292`
     - :math:`0.356`
     - **0.0** (20 of 20)

Two readings, and the second is the one to carry:

#. The **adjoint** correspondence — the section is the isotropic column
   of the harmonic frame's *physical* adjoint — holds exactly when the
   constant mode :math:`e_0` is an eigenvector of the discrete Gram, that
   is, when the rule integrates every other live basis function to zero
   (row-0 coupling zero). `[R]` Then :math:`G^{+}e_0 = e_0/G_{00}`, and
   the adjoint face, which synthesizes from :math:`G^{+}` applied to the
   moments, returns the constant :math:`\varphi\,Y_{00}/G_{00} =
   \varphi/W`. A DIAGONAL Gram satisfies this, but a DENSE Gram does too
   whenever only the higher modes couple: `[M]` 5 of the 6 dense rows
   above read round-off, and the one that breaks,
   ``folded_product(4, 4)`` at :math:`L=4`, is the one whose rule
   integrates a non-constant mirror-even harmonic to a nonzero value
   (row-0 coupling :math:`0.292`, deviation :math:`0.356`). ⚠ The
   discriminator is that one row of the Gram, not the geometry and not
   the Gram's diagonality: the slab and sphere rows read alike because
   the frame is built from ``sn.quad`` and knows nothing about the
   spatial coordinate system.
#. The **measure-level** correspondence —
   :math:`E = \text{reconstruction}(e_0\,\cdot)/W` — never touches a
   metric, so it holds in *every* configuration, the broken dense row
   included. It is bit-exact (``array_equal`` 20 of 20) on 8 of the 10
   rows and two to three ULP off on the other two
   (``level_symmetric(4)`` and ``product(4, 8)``). That is
   the honest statement of the relationship, and it is also the reason
   the angular marginal is built from the nodal axis's own unit and
   measure rather than lifted out of the harmonic frame: the marginal
   rides the measure, and it must keep working where the harmonic
   frame's adjoint column does not match. (The table was measured on
   2026-10-18, when the section's divisor was the rank-one frame's Gram
   entry; on ``gauss_legendre(4)`` every spelling of the mass is 2.0, and
   on the dense rows the divisor has since moved by at most a few ULP to
   the correctly rounded sum.)

Until the ERR-080 repair (#429, 2026-09-02) this table measured a 1-D
rule bound to the full spherical harmonics on a fabricated azimuth:
`[M]` 2026-08-24, the :math:`L=2` rows read a Gram off-diagonal of
:math:`1.155` and an adjoint deviation of :math:`16.17` on slab and
sphere alike (``gauss_legendre(8)`` read the same), and the reading
drawn from them, that the adjoint correspondence needs a DIAGONAL Gram,
was the sufficient condition taken for the exact one.

.. _spaces-collapse-pair-gates:

Verification — cite the gate, never copy its numbers
------------------------------------------------------

Three files, ``foundation``-tagged throughout. Each row's docstring names
the mutation that reddens it. The software invariants of the
construction carry no ``verifies(...)``; the modal marginal,
:eq:`spaces-axis-marginal-integral`, is an equation claim, and the
integration family G2 is its gate.

``tests/gates/numerics/test_axis_marginal_lift.py`` holds the six
families of the S1 specification
(``scratch/numerics_untangling/packages/S1_gate_spec.md``) and the rows
the S1 review added, 314 rows,
over a registry built from production factories (nodal axes of six
S\ :sub:`N` problems, the twelve production moment spaces, five
frame-dressed heads ⊗ energy ⊗ spatial, two continuum heads ⊗ cells, the
dense-form spaces ``D1``–``D3``, and synthetic, zero-weight and
zero-mass axes):

.. list-table::
   :header-rows: 1
   :widths: 22 10 68

   * - Family
     - Rows
     - What it pins, and what it cannot see
   * - G1, the unit is the constant
     - 37
     - nodal units are ones (``array_equal``); a coefficient head's unit
       is :math:`e_0` and the function it represents, evaluated in the
       test with numpy and scipy, is 1; a moment head's unit equals the
       moments of 1 computed in the test with an independent evaluation
       of the basis and ``math.fsum``, and represents 1 at the nodes; the
       LD spatial-moment axis's unit is :math:`e_0` and its marginal is
       slot 0, the cell average, read from the array directly.
   * - G2, integration
     - 33
     - the modal marginal :eq:`spaces-axis-marginal-integral`: on
       coefficient heads against a :math:`40\times80` product rule on the
       sphere built in the test; on the twelve production moment spaces
       the commuting square :math:`R_{\rm head}(\Phi^{\mathsf T}W\psi) =
       R_{\rm angular}(\psi)`; and the mass of a modal head, :math:`4\pi`
       or :math:`\Sigma w`. With G1, the only family that sees a wrong
       unit.
   * - G3, reciprocity
     - 58
     - the product condition from the raw axes and forms, then
       :math:`\langle R\psi,\varphi\rangle = \langle\psi,
       \pi^{*}\varphi\rangle`, over 38 (space, axis) rows, modal heads
       with dense forms included; the closed form against the generic
       sandwich ``AdjointOperator(R)`` on the 19 modal rows (to
       :math:`64\varepsilon` of the largest entry), the X4 witness that
       on a singular head the two adjoints still agree because
       :math:`u = B e_0 \in \operatorname{range}(G^{+})`; and the refusal
       of a non-product metric on an axis-built space. ⛔ Blind to the
       unit by construction: the algebra is self-consistent for every
       nonzero :math:`u`.
   * - G4, the section laws
     - 38
     - :math:`R\circ E = \mathrm{id}` at ``nulp = max(4, N_k)`` and
       :math:`E\circ R` idempotent. ⛔ Blind to the unit and the metric,
       because the divisor is built from the same :math:`m`.
   * - G5, the nodal pair and the mass
     - 136
     - on 16 nodal spaces, :math:`R`, :math:`\pi^{*}` and :math:`E`
       ``array_equal`` against the hand einsum over the raw weights, the
       broadcast and ``broadcast(φ / E.mass)``; on 120 shipped
       rules, the mass ``array_equal`` to ``math.fsum(w)``.
   * - G6, the refusals
     - 12
     - every row of the admission table above by its shortest distinctive
       message fragment, the wiring order (energy before single-axis), the
       kind's refusals (neither arm; a zero unit; a unit over another
       shape), the shape guards, the second
       pullback, and the memo (one marginal instance per (space, label)).

``tests/gates/numerics/test_retraction_adjoint_is_the_pullback.py`` holds
the closed-form adjoint, S6.1 to S6.4 (:ref:`spaces-collapse-pair-pullback`
and its gate table below the two lifts), and the two axis verbs.

``tests/gates/numerics/test_axis_marginal.py`` holds the nodal battery:

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
       negative leg: the same pairing under a deliberately stripped
       spatial measure must break at O(1), because a positive reading
       alone cannot discriminate metric-loaded from metric-blind.
   * - ``TestTwoArrows``
     - :math:`R^{\dagger} = (\Sigma w)E` (the anti-ERR-051 row) and the
       *type* discrimination: ``R.H`` must not be an
       :class:`~orpheus.numerics.axis_marginal.AxisSectionOperator`, and
       must be an :class:`~orpheus.numerics.axis_marginal.AxisPullbackOperator`.
   * - ``TestShippedKernelEquivalence``
     - bit-identity with the canonical angular reduction and the
       isotropic-source kernel on the real :math:`S_N` carrier, and
       that the angular marginal **is** the carrier's scalar bulk
       space. ⚠ Both equivalence rows are pinned against kernels
       hand-spelled **in the test**: the production targets route
       through the very operators under test, so a production
       comparison would be tautological.
   * - ``TestAxisGeneric``
     - the arrows are not angular-only: an untyped axis is admitted
       whatever its label, and a multi-dimensional axis contracts all
       of its dims with its own measure.
   * - ``TestAdmission``
     - the by-label refusals, the single-axis refusal, the shape guards
       in both directions, and the zero-mass asymmetry.
   * - ``TestFrameInduction``
     - the **tightness** row: the arrows against a single-region
       indicator frame built independently in the test
       (:eq:`spaces-collapse-rank-one-gram`), a comparison of two
       implementations since no production code builds that frame; the
       divisor ``array_equal`` to ``math.fsum``; the GL8 divisor as a
       1-ULP bound; the clause-2 energy refusal; and the one-mint memo.


.. _spaces-fences:

What is NOT built (the standing seams)
======================================

Stated explicitly so no reader mistakes a chartered design for shipped
machinery, and so the next phase does not re-derive a decision already
taken.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Not built
     - Where it lands, and what stands in for it today
   * - ``SpatialAxis`` / ``QuadratureAxis``
     - **CS2, unscheduled.** The coefficient families need no subclass:
       since S7 one generic coefficient axis carries the basis's layout
       (ruling 81). The only axis subclass is
       :class:`~orpheus.numerics.axis.EnergyAxis`; the spatial factor of
       :attr:`MaterialMesh.bulk_space
       <orpheus.transport.mesh.material_mesh.MaterialMesh.bulk_space>` is
       a generic :class:`~orpheus.numerics.axis.Axis` labelled
       ``"spatial"``, and the quotient point is a generic ``Declared``
       instance.
   * - Axis-built COMPOSITE and TRACE spaces
     - **CS2.** ⚠ Re-measured 2026-08-24: the *bulk* half of this
       fence has fallen. Campaign 1 CS4b moved the angular family onto
       axis-built carrier mints, so `[M]` on a shipped 1-D
       :class:`~orpheus.sn.problem.SNProblem` the scalar bulk
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
   * - The rank-:math:`d` measure-to-axis pairing
     - **CS2.** A :class:`~orpheus.numerics.measure.DiscreteMeasure` is
       a flat atom list, so :meth:`Axis.of_measure
       <orpheus.numerics.axis.Axis.of_measure>` mints at rank 1;
       the rank-:math:`d` arm of :attr:`MaterialMesh.bulk_space
       <orpheus.transport.mesh.material_mesh.MaterialMesh.bulk_space>`
       therefore is ``Declared``, with no atoms, **by contract** (G6b),
       gated as such. Inverting that row is a deliberate CS2 act. Until
       then a rigid shift of a :math:`d \ge 2` mesh is the same space
       while a shifted 1-D slab is not (ruling 85,
       :ref:`spaces-axis-origin`).
   * - ``FunctionSpace.manifold`` — the level-1 slot
     - **#429 tracker 2.0c.** A space records the index shape of its
       DOFs and not the point set those DOFs discretise, so the
       manifold is smuggled through a NAME STRING. ✅ **Half discharged
       2026-09-01 by tracker 2.1**: ``basis/indicator_basis.py`` used to
       hard-code ``f"L2[coarse_cells_R{self.ndim}]"`` — **false** for the
       energy-grid basis, an index partition calling itself spatial, and
       ``==``/hash-equal to a same-sized spatial space — and now derives
       ``f"L2[coarse_cells({self.domain.name})]"`` from a
       :class:`~orpheus.geometry.manifold.Manifold` the caller declares.
       ⚠ **What remains is that the SPACE records only the string.**
       `[M]` both producers now interpolate a typed
       :class:`~orpheus.geometry.manifold.Manifold` — ``measure.py:371``
       is ``f"L2[{self.support.name}]"`` since tracker 2.0c retyped
       ``support``, and the indicator basis is the line above — but a
       ``FunctionSpace`` still holds the resulting *name*, so the two
       producers agree by discipline rather than by construction.
       ``tests/gates/numerics/test_basis_domain.py::test_d6`` pins that
       agreement, and pins the one pair that does NOT yet agree in
       spelling (``LossKernelBasis``'s bare label against ``IndexSet``'s
       ``index(...)``). ⛔ This clause read *"`[M]` ``measure.py:331``
       still derives ``f"L2[{self.support}]"`` from a ``str``"* until
       2026-09-01: true when written, and repealed hours later by 2.0c,
       which is the campaign's own step. The level-1 view of this seam,
       and the rest of that migration, is at :ref:`manifold-seams`.
   * - The identity flip (an axis-built space IS its axis tuple)
     - ✅ **LANDED 2026-09-07 — the identity flip, CS4c step 6**
       (structural ``__eq__``), and made one key at S7
       (:ref:`spaces-identity-key`): an axis-built space is its axes and
       their dense overlay; an axes-less leaf is its class, name and
       shape, which for the digest-named leaf classes IS content
       identity; a composite is its members' keys; one of each is never
       the same space.
       ⛔ This row read *"S3. Until then the derived name is the bridge,
       and axes is declared compare=False"* until 2026-09-07. Its first
       half is now history
       (:ref:`spaces-identity-bridge`); its second was never the whole
       mechanism — the field stays ``compare=False`` so a subclass's
       dataclass-generated ``__eq__`` never reaches an ndarray, while the
       manual ``__eq__`` reads it. ⚠ The "S3" of that retired sentence is
       this plan-internal step, NOT the landed **CS4b S3** content-digest
       re-key that gave the composites and trace spaces their names —
       which is why the landed thing is named *the identity flip* here
       and everywhere after it.
   * - Retiring the densifying ``__mul__``
     - ✅ **LANDED 2026-09-07 — CS4c step 6 item 6.2a.** ``*`` survives
       and stopped densifying: the dense outer-product weights builder
       and the mixed-product bridge that fed it are gone, and a product
       never populates ``inner_product_weights``. What remains of the
       twin is *axes-lessness*, not densification — ✅ item **6.2c-ii**
       (2026-09-08) axis-ified both angular heads and closed that for the
       un-widened product, and ✅ item **6.2c-iii** (2026-09-08) closed
       it for the WIDENED one too, by making the tail the discretization
       scheme's own MODAL ``moment_axis`` and retiring the axes-less
       ``SpatialMomentSpace`` (:ref:`spaces-moment-head-axis-built`).
       ⛔ This clause read *"a WIDENED moment product is still axes-less,
       because its* ``SpatialMomentSpace`` *tail is"* until 6.2c-iii
       landed the same day.
       ⛔ And this row read *"CS2. The legacy* ``*`` *path is documented
       above and kept working; its own gates live in a separate test
       module so the retirement is a file-level move"* until
       **2026-09-07**. The
       first clause is now history; the **second was wrong about the
       mechanism**, and the landing says so. A gate that pins
       *behaviour* migrates with the behaviour, it does not travel with
       a file: ``test_space_algebra.py``'s two dense-slot rows
       (``test_inner_product_factorises_weighted``,
       ``test_inner_product_mixed_euclidean_and_weighted``) and
       ``test_space_of_axes.py``'s mixed-product third leg were
       **re-keyed in place** onto the factored metric's applied VALUES —
       each still asserts the same outer product, now as
       ``tp.apply_metric(...)`` rather than as a stored tensor — and the
       new arm-agreement band lives in
       ``tests/gates/numerics/test_tensor_product_metric_is_factored.py``. No
       file moved.
   * - The condensation morphisms on :math:`V` / :math:`V^*`
     - **Campaign 2.** Declared at
       :ref:`spaces-vv-collapse-hook`; the axis records the group
       structure they will consume. (This row said "S7" until
       2026-08-24 — a colliding plan-internal step number, see the
       warning at that anchor.)
   * - ``Medium`` and the mesh-conformity guard
     - **CS1.5.** See :ref:`spaces-symmetry-monotonicity`.
   * - Making the operators' ``space`` slot mandatory
     - ✅ **LANDED, in two halves.**  The *carrier* half — a
       constructor that will not build an unbound operator — landed
       for :math:`C` / :math:`S` / :math:`F` and the energy bindings
       at CS4c K2b / step 2 / step 3 / step 5 (two mandatory
       keyword-only ends on
       :class:`~orpheus.transport.operators.bound_operator.BoundOperator`).
       The
       *declaration* half — the ends a leaf ANNOTATES — landed for the
       :math:`L` and :math:`B` family at **CS4c step 6 item 6.4**
       (2026-09-07, ``550be030``): the eight ``domain`` / ``codomain``
       properties on
       :class:`~orpheus.sn.operators.streaming.StreamingOperator`,
       :class:`~orpheus.sn.operators.boundary.SNBoundaryOperator`,
       :class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`
       and
       :class:`~orpheus.sn.operators.boundary.SNMaskedBoundaryOperator`
       read ``FunctionSpace``, not ``Optional[FunctionSpace]``.  A
       bound leaf's ends are never ``None``.
       ⛔ This row read *"CS4. Every operator's* ``space`` *is still*
       ``FunctionSpace | None``\ *; a* ``None`` *operand makes the
       composition guard skip rather than validate"* until 2026-09-07.
       Its first clause is now history.  Its second was a claim about a
       mechanism that **had already been repaired elsewhere**: the
       composition guard's skip is keyed on the runtime VALUE, not on
       the annotation, and no L/B binding has answered ``None`` since
       the ends became mesh-derived properties; the sibling story that
       an unbound end made the ``.H`` metric branch conditional stopped
       being true at CS4c **step 1**, when
       :class:`~orpheus.numerics.operator.AdjointOperator` began
       REFUSING an unbound inner outright.  Item 6.4 is therefore
       TYPING-ONLY — no value moves — and the honest reading of the
       flip is that it makes the type system state a law the runtime
       already enforced.  ⚠ One deliberate exception survives, and it
       is not an oversight:
       :class:`~orpheus.sn.operators.boundary.RadialCharacteristicBoundaryOperator`\ 's
       CONSTRUCTOR parameter stays ``FullFieldSpace | None``, because
       the constructor **is** the narrow — its guard ("the pose carries
       no ψ½ ray") has two witnesses that pass a runtime ``None`` on
       the slab, and flipping the parameter would make a
       runtime-reachable guard statically unreachable (the harmful-stub
       trap).  The row is about the ENDS a leaf declares, not about
       every parameter that mentions a space.
   * - A CHOSEN curvilinear moment mass
     - **Not scheduled here — and the reason is a split, not a delay.**
       P7 landed the *machinery* for a non-Hadamard metric; which
       :math:`G` a curved chart's multi-moment cell wants is a
       **value** question, pinned by physical functionals and needing
       the cell-solve consumer of `#158
       <https://github.com/deOliveira-R/ORPHEUS/issues/158>`_. The
       scheme family still refuses, now naming that one blocker.
       Expressible is not known
       (:ref:`spaces-metric-not-on-the-axis`).
   * - The Riesz legs (♭ / ♯)
     - **Landed after this table was written.** The legs are operators
       built from a space, ``RieszLowerOperator(space)`` (♭, applies
       :math:`G`) and ``RieszRaiseOperator(space)`` (♯, applies
       :math:`G^{+}`), the classes
       :class:`~orpheus.numerics.operator.RieszLowerOperator` and
       :class:`~orpheus.numerics.operator.RieszRaiseOperator`, and
       :class:`~orpheus.numerics.operator.AdjointOperator` IS
       :math:`A^{*} = \sharp_{A.\mathrm{domain}}\circ A.\mathrm{dual}()
       \circ \flat_{A.\mathrm{codomain}}` — P7's metric family is
       the one metric arithmetic they wrap, as this row anticipated, so
       there is no third spelling of it. Their definition is
       :ref:`spaces-riesz-legs`; the landing is recorded in the SN
       development history (the Riesz-legs entry).


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
   * - 2026-10-22
     - **An axis is identified by what fixes its space, and a space by
       one key** (numerics untangling, package S7, halves S7a and S7b;
       rulings 58, 63, 64, 71, 81 to 88; #631). The provenance slot became
       the closed sum ``Axis.origin`` (``FromMeasure``, ``FromQuadrature``,
       ``FromBasis``, ``FromFrame``, ``Declared``), from which the kind,
       support, atoms and layout are derived; the axis content gained the
       support, the atoms and the coefficient layout; ``FunctionSpace``
       gained ``identity_key()``, the one definition of ``==`` and
       ``hash``, which separates all six measured collisions;
       ``FunctionSpace.of_measure`` and ``FunctionSpace.of_basis`` became
       the two coefficient- and sample-space mints; the basis carries its
       layout (``orpheus.numerics.coefficient_layout``), so the two head
       space classes, the two family axis subclasses, the ``MomentHead``
       protocol and ``TensorProductSpace.find_factor`` retired onto one
       generic coefficient axis and the door ``moment_head_axis``. See
       :ref:`spaces-axis-origin`, :ref:`spaces-axis-identity`,
       :ref:`spaces-identity-key` and :ref:`spaces-moment-head`; the full
       record is the row of :ref:`operator-algebra-history`.
     - #631
     - ``[LANDED e80ec738]`` (S7a ``47e7bfaa``, S7b ``ac6c6558``, the inverse-face fix ``96ec89fd``, docs ``56b3ef92``)
   * - 2026-10-19
     - **The axis owns its unit, and the marginal is built from it**
       (numerics untangling, package S1; rulings 36 to 41). Every axis
       carries :attr:`~orpheus.numerics.axis.Axis.unit`, the coordinates
       of the constant function :math:`1` in its own representation,
       minted with the axis: all ones on a ``NODAL`` axis; :math:`e_0` on
       the continuum Legendre and harmonic heads; the moments
       :math:`\Phi^{\mathsf T} w = B e_0` on a head a frame has dressed;
       :math:`e_0` on the LD spatial-moment axis. Ruling 41 makes the
       kind a closed sum, ``BasisKind.NODAL | Modal(unit)``: the modal arm
       is the value class :class:`~orpheus.numerics.axis.Modal`, which
       carries the unit, ``Axis.unit`` is read from the kind, and a modal
       factor without a unit cannot be spelled; the unit is content on
       modal axes only by that type, so every nodal digest is unchanged.
       The metric block :math:`G_k` has one answer,
       ``FunctionSpace._axis_block(k)``, shared by the space's metric and
       the covector. The three arrows
       move to :mod:`orpheus.numerics.axis_marginal` and are built as
       :math:`R = (G_k u)^{\mathsf T}\otimes\mathrm{id}`,
       :math:`\pi^{*} = u\otimes\mathrm{id}` and :math:`E = (u/\langle
       m,u\rangle)\otimes\mathrm{id}` by
       ``AxisRetractionOperator.of(space, label)`` /
       ``AxisSectionOperator.of(space, label)``, with the correctly
       rounded mass (``math.fsum``, ``AxisSectionOperator.mass``) as the
       divisor. The modal marginal
       becomes an integral (:eq:`spaces-axis-marginal-integral`), so the
       MODAL refusal and the "two measures on one block" refusal retire;
       the single-region indicator frame the mint built and discarded
       (``frame._collapse_pair``) is deleted, and so are the space's
       ``retraction``, ``section``, ``riesz_lower`` and ``riesz_raise``.
       The divisor moved on 68 of 120 shipped rules (up to 39 ULP).
       Record: :ref:`spaces-axis-unit`, :ref:`spaces-collapse-pair`.
     - —
     - ``5f850a69``, ``71896877`` on ``refactor/eprime-wave3`` *(in development)*
   * - 2026-10-02
     - **The retraction's adjoint is minted as the pullback** (#405 P1
       step 6). :class:`~orpheus.numerics.axis_marginal.AxisPullbackOperator`,
       the plain broadcast :math:`\pi^{*}`, is minted by the retraction
       from the same induced data, ``R.H`` returns it (it was the generic
       metric sandwich, 1 to 3 ULP away), ``R.H.H is R``, and a second
       one cannot be minted. The mint admits only product metrics: the
       marginal keeps every positioned form of the axes that remain
       (it dropped them, a defect reached by 0 of 1 891 suite mints), and
       a collapsed axis carrying a form is refused. Two space verbs,
       :meth:`~orpheus.numerics.space.FunctionSpace.without_axis` and
       :meth:`~orpheus.numerics.space.FunctionSpace.with_axis`, over one
       overlay reader. The S\ :sub:`N` adjoint entry lifts its detector
       by :math:`R^{\dagger}`, bit-exact. Record:
       :ref:`spaces-collapse-pair-pullback`.
     - #405
     - ``61f82a18`` on ``main``
   * - 2026-10-02
     - **The axis and every space-name digest move onto the one content
       encoder** (#405 P1 step 5). :class:`~orpheus.numerics.axis.Axis`
       takes its ``==`` and ``hash`` from
       :class:`~orpheus.values.content.ContentIdentity` over its
       dataclass fields (``label``, ``shape``, ``weights``, ``kind``, the
       generator declared ``compare=False``), which replace
       the hand-written ``_identity_key`` tuple and the private
       ``_structural_bytes`` encoder (which encoded a scalar by its
       ``repr``, so ``1`` and ``1.0`` differed, and tagged the class
       without its module). ``FunctionSpace.of_axes``, the three trace
       spaces' name mints and ``FullFieldSpace.from_blocks``, four
       encoders with no type tags or length prefixes, call
       :func:`~orpheus.values.content.name_digest`, the one 8-byte
       space-name digest, and the trace spaces read their layout's
       content through ``FaceLayout.structure``. The
       generator exclusion is answered a third time: an inclusion now
       fails for both generator types
       (:ref:`spaces-generator-identity-third-answer`). Record of the
       encoder: :ref:`structured-geometry-content-identity`.
     - #405
     - ``a5113ac0`` on ``main``
   * - 2026-09-12
     - **A** ``Quadrature`` **is a VALUE, so the generator exclusion
       becomes a RULING** (the consumers campaign, step 1 / S1b; ruling
       R-cc3). A quadrature is a generating datum of every SN problem,
       and a save state must compare by CONTENT — so
       :class:`~orpheus.numerics.quadrature.directional.Quadrature`
       became ``eq=False`` with a hand-written ``_identity_key`` over
       nodes, weights, support, invariance group, exactness claim, level
       structure and folding (the derived frame cache is excluded — it
       is derived, not generating). ⛔ That **inverts one leg** of this
       page's own argument for keeping
       ``Axis.generator`` out
       of ``_identity_key``: `[M]` a key holding a ``Quadrature`` no
       longer raises on ``==`` or ``hash`` and no longer splits the
       ``of_axes`` digest. The exclusion is unchanged and is now carried
       by the doctrine — *two axes with identical structural content are
       the same axis whatever instance produced them* — with a
       :class:`~orpheus.numerics.measure.DiscreteMeasure` generator as
       the surviving structural leg. The ``G1c`` sub-case was re-posed in
       the landing commit to assert both halves.
       See :ref:`spaces-generator-identity-exclusion`.
     - `#459 <https://github.com/deOliveira-R/ORPHEUS/issues/459>`_
     - ``2c1667b0`` on ``main``
   * - 2026-09-08
     - **The metric becomes a DERIVED object over the axes, and both
       moment heads become AXIS-BUILT** (CS4c step 6 item 6.2c, rulings
       R-6.2c-1 / R-6.2c-2). Two commits. **6.2c-i** replaced the
       "axes win, the object is ignored" short-circuit with one
       :class:`~orpheus.numerics.metric.FactoredMetric` DERIVED from the
       axes — a :class:`~orpheus.numerics.metric.DiagonalMetric` per
       weighted axis on its own block, nothing on a counting axis — and
       admits an explicit metric object beside the axes only as the
       positioned **OVERLAY** of forms: one entry per axis, in order, a
       form ONLY on a block whose axis carries no measure, MERGED into
       the derived entries rather than substituted for them. That is the
       home a dense Gram needs (*a Gram is a FORM, never on*
       :attr:`Axis.weights <orpheus.numerics.axis.Axis.weights>`), and the
       guard refuses the three ways of spelling a second source: a bare
       object beside axes, blocks that do not follow the axes, a
       positioned ``DiagonalMetric``, and a form beside a weighted axis.
       ⭐ The overlay reading is the carve's own correction — the first
       spelling let an explicit object WIN wholesale, so a caller had to
       restate every other axis's measure inside it (a twin spelling) and
       a product would have dropped a factor's form. The inline
       ``FunctionSpace._apply_axes_weights`` loop retired into
       :meth:`DiagonalMetric.apply_block
       <orpheus.numerics.metric.DiagonalMetric.apply_block>` —
       operation-for-operation the same reshape-and-multiply, so
       bit-identical on every axis-built space. **6.2c-ii** then made
       ``SphericalHarmonicSpace``
       and
       ``LegendreSpace``
       axis-built — one ``MODAL``
       ``HarmonicAxis`` /
       ``LegendreAxis``, the latter carrying
       its spent :math:`O(2)_a` axis in the identity — so the metric
       entered the identity and the corpus's chartered doctrine reached
       the moment layer (:ref:`spaces-moment-head-axis-built`). `[M]` 17
       of 75 shipped (rule, :math:`L`) rows measure a ``DENSE`` Gram and
       2 of them sit at :math:`L \le 2`, so the positioned-form path is
       exercised, not hypothetical. Consequences elsewhere: the ONE
       moment space the tree binds became the frame's Parseval-dressed
       head (:ref:`frame-the-one-moment-space`); the ``MODAL`` arm of the
       ``Axis.generator`` slot, declared and
       unbuilt since CS5, is live (:ref:`spaces-generator-seams`); and
       ``has_coordinate_cone`` got its first production ``False``.
     - —
     - ``numerics/space.py``, ``numerics/metric.py``,
       ``numerics/axis.py``, ``numerics/spaces/``
   * - 2026-08-30
     - **The metric becomes an OBJECT — a space can carry a
       non-diagonal** :math:`G` (campaign 1, phase P7). The metric stops
       being a thing that is *multiplied* into the element (a broadcast
       weight array, or the per-axis factor measures) and becomes a
       thing that is **applied**: a typed
       :class:`~orpheus.numerics.metric.HilbertMetric` family —
       :class:`~orpheus.numerics.metric.DiagonalMetric` (the Hadamard
       realization, bit-identical to the arithmetic it replaces),
       :class:`~orpheus.numerics.metric.DenseMetric` (a symmetric matrix
       on the flattened leading block, Moore–Penrose inverse face,
       refusing an asymmetric matrix or a Penrose-inconsistent pair) and
       :class:`~orpheus.numerics.metric.FactoredMetric` (a lazy
       per-block tensor product; the Kronecker product is never
       materialized). The space gains a third metric SOURCE and a
       three-arm exclusivity guard — of whose arms the
       ``(weights, metric)`` one had been structurally unreachable — and
       the three propagation sites that used to copy a weight array are
       taught: :meth:`DualSpace.of
       <orpheus.numerics.space.DualSpace.of>` threads the object,
       ``TensorProductSpace.from_factors`` grows a dense-factor arm, and
       ``FrameBase.gram`` (retired at item 6.2c-ii in favour of the
       :meth:`gram_inverse <orpheus.numerics.frame_faces.CrossGramInverse.of>` arrow)
       **strips** it while installing the row-sum probe (`[M]` the
       pre-P7 spelling returned ``[7.0, 11.0]`` for a projection whose
       true value is :math:`[8/3,\,16/3]` — a silent value error, and
       now an unconstructible state). The founding consumer is the
       frame: a ``DENSE`` verdict is no longer a refusal but the matrix
       pseudo-inverse dressing, so **Parseval is a theorem on every
       frame** (:eq:`spaces-pseudo-inverse-parseval`), the recorded F-0
       limitation is repaired, and `[M]` a production analysis adjoint
       moves by :math:`98\,\%` in Frobenius relative — a correctness
       repair, not plumbing. What did NOT change:
       :attr:`Axis.weights <orpheus.numerics.axis.Axis.weights>` stays a
       1-D measure (a measure is diagonal by nature; a Gram is a form),
       identity stays metric-blind, and the curvilinear moment mass
       stays refused on its VALUE alone
       (:ref:`spaces-metric-not-on-the-axis`). Found by the flip: all
       three :class:`~orpheus.numerics.space.FunctionSpace` subclasses
       that override ``__post_init__`` did so **without chaining**, and
       had been outside every base construction guard since they were
       written.
     - `#409 <https://github.com/deOliveira-R/ORPHEUS/issues/409>`_
     - ``6a0e0473`` (the family), ``bae73fa7`` (exclusivity + the
       factored product), ``f1f30cea`` (the dense dressing),
       ``af9f95f1`` (the curvilinear refusal names one blocker) —
       *(in development,* branch ``feature/p7-nondiagonal-metric``\ *;*
       ``git`` *is the merge-status authority)*
   * - 2026-08-29
     - **A producer binds the space factor, and reads the angular
       geometry through it** (the streaming campaign's P4-remainder) —
       the first solve-time consumer of
       ``Axis.generator``, which discharges
       CS5's third seam one day after it was declared.
       :class:`~orpheus.sn.mesh.reduced_operator.ReducedStreamingOperator`
       gains an
       :attr:`~orpheus.sn.mesh.reduced_operator.ReducedStreamingOperator.angular_axis`
       field, minted by each of the three streaming factories from the
       measure they were already handed, and
       :meth:`~orpheus.sn.mesh.reduced_operator.ReducedStreamingOperator.streaming_terms`
       recovers the radial cosine and the level fibration **through** it.
       The recovery goes via one accessor,
       ``Axis.generator_as``
       — a typed narrow whose refusal names both the axis and the asking
       consumer, and which is load-bearing rather than defensive because
       the declared union's other arm cannot answer either read. The
       angular-closure family's construction contract widens to
       ``cls(angular, pairing, angular_axis)``, so both concrete members
       mint through the same narrow. In the other direction the
       **courier dies** — the field whose docstring said it was *"held
       so a consumer that needs the weights, the level partition or the
       cosines does not have to be handed the quadrature separately"*,
       which is a weld's own confession:
       :class:`~orpheus.sn.angular.redistribution.AngularRedistribution`
       sheds its ``quadrature`` field and is pure :math:`\alpha` data —
       ``coord``, ``alpha_per_level``, ``mu_start_per_level`` — pinned
       structurally by field-set equality so a re-addition reds. Two
       supporting single-source moves: ``Axis.of_quadrature``'s label
       **defaults**, making a label twin across mint sites unspellable,
       and the cylinder admission probe's defaulted ``getattr`` on
       ``level_structure`` hardens to a direct read now that the contract
       declares the member. Every value comparison over this re-point is
       ``X == X``, so the acceptance evidence is a **route** gate with a
       decoy generator, not a value gate
       (:ref:`spaces-generator-route-gate`).
     - —
     - ``ac485104`` (the dead ``_weight_of`` retires), ``ad04e236`` (the
       binding, the courier's death, and the gates), ``1fb70c15`` (the
       admission probe reads the declared contract) — *(in development,*
       branch ``feature/p4rem-producer-binds-axis``\ *;* ``git`` *is the
       merge-status authority)*
   * - 2026-08-29
     - **An axis can name the generator that made it** (campaign 1,
       phase CS5). The axis gains a fifth slot,
       ``Axis.generator`` — **provenance,
       never identity** — and the mint routes THROUGH the generator:
       :meth:`Axis.of_measure
       <orpheus.numerics.axis.Axis.of_measure>` (the
       axis-composed sibling of ``DiscreteMeasure.space``,
       ``NODAL`` by construction) and :meth:`Axis.of_quadrature
       <orpheus.numerics.axis.Axis.of_quadrature>`
       (delegating the structural mint down and upgrading only the
       provenance). An axis is a **forgetful map** from its generator —
       it keeps the weights and drops the nodes — and the mint is that
       map's **section**, so the forgetting is recoverable
       (:eq:`spaces-axis-generator-section`): a consumer holding only
       the SPACE recovers ``mu_x`` / ``eta`` / ``mu_z`` /
       ``level_indices``. The generator is EXCLUDED from
       ``_identity_key`` — recorded here as *structurally mandatory, not
       taste, since an inclusion makes* ``Axis.__eq__`` *and*
       ``hash(Axis)`` *RAISE*, which was true of the tree on this row's
       date. ⛔ **The quadrature half of that reason inverted on
       2026-09-12** (S1b gave ``Quadrature`` content identity); the
       exclusion survives as a RULING, with the ``DiscreteMeasure`` leg
       still raising
       (:ref:`spaces-generator-identity-exclusion`). **Two** of the
       tree's three nodal mint sites consume it, bit-identically:
       :attr:`SNProblem.angular_bulk_space
       <orpheus.sn.problem.SNProblem.angular_bulk_space>`
       collapses to ``Axis.of_quadrature(self.quad)``, and the rank-1 arm of
       :attr:`MaterialMesh.bulk_space
       <orpheus.transport.mesh.material_mesh.MaterialMesh.bulk_space>`
       mints through the carrier's own ``volume_measure``; the third,
       the homogeneous pose's counting point, keeps its literal and is
       documented as honestly generator-less. The rank-:math:`d` spatial
       arm likewise stays literal, by contract. The chartered ``DiscreteMeasure``-typed accessor
       was **refuted** by the phase's opening ground measurement — a
       bare measure answers three of the done-when's four names
       (:ref:`spaces-generator-why-quadrature`). Sibling repair: the
       ``AngularMeasure`` Protocol in
       :mod:`orpheus.sn.angular.redistribution` now declares the
       ``level_structure`` member its cylinder-factory consumer was
       already reading past the contract
       (:ref:`spaces-generator-protocol`).
     - —
     - merged @ ``cb3cd15b`` — ``4e7b8977`` (the slot, the two mints,
       the three consumers, the gates), ``b0bfc06c`` (the Protocol
       declaration), ``cb3cd15b`` (the gate-quality repairs the docs
       pass's own census found: the roster's fifth factory, the
       unmarked gate class, the module head still saying "four slots")
   * - 2026-08-24
     - **The space becomes the construction key, and it mints the
       collapse pair** (campaign 1, phase CS4b, steps S5–S7).
       *Construction goes space-primary* (S5): the carrier gains
       :attr:`SNProblem.angular_trial_space
       <orpheus.sn.problem.SNProblem.angular_trial_space>`
       (the scheme-widened angular mint), the composite allocators go
       space-keyed, and the mesh-keyed leaf **sugar tier is deleted** —
       every call site now names a space, not a carrier.
       *The axis collapse pair is minted on the space* (S6): the verbs
       ``FunctionSpace.retraction`` / ``FunctionSpace.section`` return the
       split epi/mono pair, memoized per axis label, and S6.0b re-carved
       their REALIZATION so both are the induced output of a
       single-region indicator frame built and discarded at one site —
       the **stage-2 generator discipline**, with the section's divisor
       read off the frame's :math:`1\times1` Parseval metric rather
       than chosen by hand (:ref:`spaces-collapse-pair`). The angular
       reduction and the isotropic-source projection are re-keyed onto
       that pair, so each has ONE realization tree-wide, and
       ``AxisEmbeddingOperator`` is renamed
       :class:`~orpheus.numerics.axis_marginal.AxisSectionOperator` on
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
       ⛔ **Both S7 halves were superseded at the CS4c coda
       (2026-09-08)** — do not go looking for either in the tree. The
       mesh-less 1-cell carrier retired with its ``from_materials``
       factory, taking the :math:`S_N`-promotion refusal and one of
       ``areas``' three arms with it (two remain, and ``mesh is None``
       carries one meaning again, so there is nothing left for S7's
       ``ndim`` discrimination to discriminate); and the rate
       co-vectors' cross-section fields are now *born* on the pose, so
       the re-posing step does not exist rather than being done
       correctly. See :ref:`homogeneous-development-history`.
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
       (``from_grid`` / ``synthetic``, since layer_dag.md P2 ``of_edges`` /
       ``EnergyGrid.as_axis``; identity = :math:`n_g` + edges
       CONTENT; weighted axes refused by the counting-measure theorem).
       :meth:`FunctionSpace.of_axes
       <orpheus.numerics.space.FunctionSpace.of_axes>` composes a space
       as the ordered product of its axes with a **per-axis metric
       path** (no densification) and a deterministic, injective
       **derived name** — the identity bridge that made "metric
       differences imply space differences" true from that day, three
       weeks before the identity flip made it direct
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
       a meshed carrier all fell out of the same body (the degenerate
       carrier retired at the CS4c coda, 2026-09-08; the uniformity
       claim is unchanged for the two survivors, and a genuine
       unit-width one-cell mesh's ``bulk_space`` is still ``==`` the
       homogeneous problem's pose — the identity-bridge gate pins it);
       the homogeneous
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

   - ``tests/gates/numerics/test_axis.py`` — the intrinsic laws of the axis
     concept: rank, measure canonicalization, the refusals, structural
     identity per subclass, and the ``synthetic`` / ``of_edges`` (until
     layer_dag.md P2, ``from_grid``) inequality.
   - ``tests/gates/numerics/test_space_of_axes.py`` — composition: shape
     concatenation, the per-axis metric against an independently built
     reference, the no-densification proof, the derived name's
     determinism across processes, and ``has_coordinate_cone``.
   - ``tests/gates/numerics/test_field.py`` (gates E1/E2) — the cone
     consult's **positive and negative pair**: a MODAL space REFUSES
     with a typed error naming the space, and the same values on an
     all-NODAL space answer exactly what the legacy path answers.
   - ``tests/gates/homogeneous/test_operator_spaces.py`` — the positive floor
     (all five homogeneous operators plus :math:`K` report the SAME
     space), the refusal witnesses (a 2g-vs-4g sum; :math:`M^{-1}(2g)
     \circ F(4g)`), the energy arm's gridded-vs-``synthetic``
     discrimination, and the ``.H`` **loaded/blind pair** required by
     ``vv-principles`` anti-pattern #19 — a bit-identity leg for the
     shipped scalar-metric case, paired with a deliberately
     non-physical per-group-weighted axis on which ``.H`` demonstrably
     MOVES.
   - ``tests/gates/homogeneous/test_kinf_exact_reference.py`` — the
     homogeneous solve's value gate: :math:`k_\infty`, the gauged flux and
     both condensed cross sections against the exact rational answer for
     their float inputs, within a forward-error bound derived per case
     (:ref:`homogeneous-exact-reference`), over every producing mixture the
     tree ships.  The evidence for the "no value motion" claim above is its
     predecessor's: the migration gate ``test_byte_stability.py`` pinned
     the solve bit-exactly (``np.array_equal`` and exact ``==``) against a
     baseline captured immediately before the wiring, and `[M]` it held
     across the rewiring on 2026-08-20.  It retired on 2026-10-01, when a
     macOS update moved the LAPACK eigen-solve under it by 1 ULP with no
     ORPHEUS change: a byte pin on a platform library's output pins the
     platform.
