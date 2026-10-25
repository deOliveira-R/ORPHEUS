.. _architecture-conceptual-view:

The conceptual view — a Problem, a Solution, and the mathematics between them
=============================================================================

This page maps ORPHEUS as mathematics: what the code is made of, in the
order a reader meets it, with every concept pointing at the theory page
that defines it and at the class that realises it. It defines nothing
itself. The layering page (:ref:`architecture-layering`) says where each
package lives; this page says what the packages are, together. A concept
whose definition the corpus does not yet carry is said to be owed, with
the nearest page, rather than defined here.

The tree it describes is ``main`` at ``b4865b5e`` (2026-09-21), except
that every row and paragraph on the boundary was rewritten at the tip of
the wave-2 boundary carve, ``471c4737`` (2026-10-24), which replaced the
trace-law architecture by the face role, the wall condition on a cone, the
binding and the per-method realizer table (the present architecture is
:ref:`bc-architecture-of-record`, and each retired name maps to its
successor in :ref:`bc-wave2-crosswalk`). Every
class name is a live role, so a rename without an update here is a dead
reference; the gate is the Nexus ``dead_references`` sweep, not the Sphinx
build, which renders an unresolved role as plain text without a warning
for the many modules here that no autodoc page renders. Every count
carries its command in the three memos the page was written from
(``scratch/_claude_md/explorer_{problem_side,solution_side,theory_map}.md``,
untracked).

.. contents::
   :local:
   :depth: 1


Everything is a map on spaces
-----------------------------

Everything in the problem layer is a map on spaces. A **field** is a map
*space → values*: cross sections map the space to a cross-section value,
the flux maps it to a flux value. An **operator** is a map *space →
space*. The organising object of the codebase is therefore the **space
with labelled axes**, and nothing fundamental requires a mesh object. A
space is the ordered product of its axes (:eq:`spaces-axis-product`,
:ref:`spaces-the-axis`); an axis is identified by what fixes its
factor's space, *(label, shape, factor measure, kind, support, atoms,
coefficient layout)*, all of it derived from its **origin**, the closed
sum ``FromMeasure | FromQuadrature | FromBasis | FromFrame | Declared``
that records what minted it and is provenance, never identity
(:ref:`spaces-axis-generator`). A space's identity is one key,
:meth:`~orpheus.numerics.space.FunctionSpace.identity_key`, read off its
axes (and any dense form positioned over them) or off its members
(:ref:`spaces-identity-key`), so metric differences, atom differences and
domain differences all imply space differences. Realising a law or a kernel
means binding it to a space, a domain and a codomain
(:ref:`bound-operator`): that binding is what makes its properties
concrete, and it is why the adjoint below is a composition and not a
second implementation.

Classes: :class:`~orpheus.numerics.space.FunctionSpace` and
:meth:`~orpheus.numerics.space.FunctionSpace.of_axes`;
:class:`~orpheus.numerics.axis.Axis` with
:class:`~orpheus.numerics.axis.BasisKind`;
:class:`~orpheus.numerics.field.Field`;
:class:`~orpheus.numerics.operator.LinearOperator`.


The Problem half — posing is a filtration
-----------------------------------------

The Problem is the material mesh augmented, stage by stage, until the
phase space is fully defined on every axis, and it carries the posed
question. The chain of commitments is a **filtration**, a monotone
refinement of partitions of phase space: materials commit the
per-material partition and the energy axis's data; the geometry overlay
commits regions, identifications and boundary data; the mesh commits
cells and the spatial measure; state fields put data on the cells; the
method head commits the terminal refinement of every axis, and only
there does every axis formally construct. Coarsening is the same chain
walked backward: the axis marginal is the projection onto a coarser stage
(:ref:`spaces-collapse-pair`), and homogenisation and condensation are its
Petrov–Galerkin form (:ref:`sn-homogenization-petrov-galerkin-frame`,
:ref:`sn-condensation-petrov-galerkin-frame`).

In the tree the augmentation is two-staged, data then behaviour.
:class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` is the
method-agnostic data (the axes, the region-to-material map, the group
count, the cell volumes and the cross-section field
:class:`~orpheus.transport.mesh.material_xs_field.MaterialXSField`,
which is a Problem datum, :ref:`sn-sigma-is-a-problem-datum`), and a
method's problem *is* a material mesh that adds one method's behaviour:
:class:`~orpheus.sn.problem.SNProblem` adds the quadrature, the spatial
scheme, the angular closure, the scattering order and the declared boundary
conditions (a map from each integer tag the faces carry to a wall condition);
:class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh` adds the scalar
trace and the realised albedo operators. The two are the witnesses of the
:class:`~orpheus.transport.method.TransportMethod` Protocol, whose one
body :func:`~orpheus.transport.method.realize_boundary` walks the faces of
the problem's binding to the method's realizer table
(:ref:`bc-realizer-layer`). The third hub,
:class:`~orpheus.homogeneous.solver.HomogeneousProblem`, is posed from
one :class:`~orpheus.data.macro_xs.mixture.Mixture` on the energy axis
alone. "Problem" is thus a concept with three realisations and no shared
type, which is the first entry of the debt list below.

The hub (:ref:`the-problem-hub`) owns what is consumed: its spaces as cached mints
(:attr:`~orpheus.transport.mesh.material_mesh.MaterialMesh.bulk_space`,
the SN angular bulk, trial, trace and full-field spaces, the moment space
of the angular frame, :ref:`frame-moment-space-single-home`), its fields,
and its bound operators, of which there is one fission operator per
Problem (:ref:`sn-one-fission-per-problem`). Its last step is the
**operator pencil**, on a transport hub :math:`(A, F)`, the loss and the
fission production, with the family :math:`A - \sigma F`
(:ref:`the-operator-pencil`, :eq:`pencil-family`, where the general pencil
is written :math:`(A, M)`; :class:`~orpheus.numerics.pencil.OperatorPencil`),
on one space, with no inverse and no resolvent of its own. The two questions are typed on the
pencil (:ref:`eigenvalue-posing`;
:ref:`sn-the-problem-poses-its-pencil`):
:class:`~orpheus.numerics.posing.EigenPosing` is the homogeneous question
:math:`A\psi = \mu F\psi` read through a spectral map
(:data:`~orpheus.numerics.posing.K_MAP`, :math:`k = 1/\mu`), whose
unknown is a ray in the cone plus a scalar;
:class:`~orpheus.numerics.posing.SourcePosing` is the affine question
:math:`A\psi = q`, whose unknown is a coset. A subcritical multiplying
source is the composition ``SourcePosing(pencil.at(1), q)``, never a
third type. The affine form itself is the wall condition's,
:eq:`affine-bc-form`, and there is no affine operator by ruling: every
affine condition is a linear operator plus a typed source.

The posings are the question *bound* to a pencil. What is asked exists
before any pencil, as a physics-free value
(:ref:`structured-geometry-question-values`):
:class:`~orpheus.numerics.question.Eigen` names one direction of the
system's parameter space by an opaque key, a base point by its offsets
from the physical point, and a mode (the pole wanted);
:class:`~orpheus.numerics.question.FixedSource` holds a source and
:class:`~orpheus.numerics.question.Response` a detector, so the role is
the type and no value carries an adjoint flag. The k-eigenvalue, the
classical c-eigenvalue, a boron search and a critical size are four keys
of one ``Eigen``; the keys that exist are a set of cells,
:class:`~orpheus.data.cells.CellCoefficient` (k is every fission-emission
cell), and the width of one interval,
:class:`~orpheus.geometry.extent.GeometryExtent`. A question is written
down with the system it is asked of as a **reference specification**
(:ref:`structured-geometry-specification`), the key a reference cache
stores an answer under. The layer of the filtration it is posed at is its
type: :class:`~orpheus.specification.specification.InfiniteMediumSpecification`
holds one material and is posed on energy alone, the infinite medium being
the point in position and in direction (:ref:`infinite-medium-definition`),
and :class:`~orpheus.specification.specification.GeometrySpecification`
holds the materials and a geometry with its face roles. Either is admitted in a
canonical form: spectators dropped, every key resolved to a coordinate of
the problem, the datum fitted to its groups, regions and readable
coordinates. Nothing binds a question to a system's operators yet: deriving
the pencil and the spectral map from a parameter, and the mode law, are the
posing sequence's unit 6 (#529).


The objects on an axis — measure, basis, frame, cone
----------------------------------------------------

A **discrete measure** :math:`\mu = \sum_i w_i \delta_{x_i}`
(:eq:`discrete-measure-definition`) carries nodes, weights, a typed
support :class:`~orpheus.geometry.manifold.Manifold` and an invariance
group; tensor product, direct sum, pushforward along a typed
:class:`~orpheus.geometry.manifold.ManifoldMap`, restriction, partition
and quotient are its verbs (:class:`~orpheus.numerics.measure.DiscreteMeasure`).
A quadrature wraps one such measure and mints the angular axis and the
angular frame from it (:class:`~orpheus.numerics.quadrature.directional.Quadrature`).
The space a measure samples is :math:`L^2(\mu)`,
:meth:`FunctionSpace.of_measure(measure, label)
<orpheus.numerics.space.FunctionSpace.of_measure>`: one nodal axis that
keeps the weights, the atoms and the manifold they lie on. A **basis**
spans a coefficient space, :meth:`FunctionSpace.of_basis(basis)
<orpheus.numerics.space.FunctionSpace.of_basis>`: one axis on the basis's
domain, indexed by the coefficient layout the basis carries
(:mod:`orpheus.numerics.coefficient_layout`: the rectangular harmonic
table, the flat Legendre row, one coefficient per cell of a partition, or
a row of computed modes) and measured by its continuum Gram, modal for a
spectral family and nodal for a partition
(:ref:`spaces-moment-head-axis-built`). There is one generic coefficient
axis; the family is its layout, not a subclass. The energy axis is a
one-dimensional mesh in energy whose one-cell limit persists
(:ref:`spaces-counting-measure-theorem`, :ref:`spaces-energy-grid-is-a-mesh`;
:class:`~orpheus.numerics.axis.EnergyAxis`). A basis, the choice-free
synthesis side of a frame, is defined at :ref:`spaces-basis`, within the
three-level picture of the manifold page (:ref:`manifold-three-levels`);
the class is :class:`~orpheus.numerics.basis.base.Basis`.

The flux lives in the **positive cone**
(:ref:`cone-ordered-vector-space`, :eq:`positive-cone-definition`). The
cone is a predicate, never a class, by ruling
(:ref:`cone-membership-is-a-predicate`): a field answers
:meth:`~orpheus.numerics.field.Field.cone_violations`, a space answers
whether its coordinates carry a cone at all (nodal yes, modal no), and a
scheme declares whether it preserves it.


Stage-2 generators — frames and schemes
---------------------------------------

A **frame** is a *(basis, measure)* pairing (:ref:`galerkin-projection`;
the definition sits at the frame page's "The discrete frame" section,
:eq:`galerkin-pair`). It induces the space *and* the operators at one
site: the basis fixes the codomain, dressed with the frame's Parseval
metric, the inverse discrete Gram (:ref:`frame-parseval-metric`,
:eq:`frame-discrete-gram`), and the measure fixes the domain; it emits the
**analysis** face :math:`M` (samples to coefficients) and the
**reconstruction** face :math:`R` (coefficients to values), and the two
close on coefficients, :math:`M R = c_V I` (:eq:`galerkin-frame-idempotency`),
so :math:`R \circ M` is a projector up to the frame's scalar. The
projection discipline is the type:
:class:`~orpheus.numerics.frame.PetrovGalerkinFrame` holds an explicit
test basis (:eq:`petrov-galerkin-construction`),
:class:`~orpheus.numerics.frame.GalerkinFrame` binds test to trial
(:eq:`galerkin-construction`), and then the two faces are each other's
adjoint up to one scalar, :math:`M^{*} = R/W` under the Parseval metric on
a diagonal Gram (:ref:`frame-square-closure-section`; the bare
:math:`M^{*} = R` needs an orthonormal basis, which the real harmonics
are not), and
:class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame` is the
transport realisation minting carrier-typed faces. The spatial
**scheme** (:ref:`discretization-closures`;
:class:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase`)
plays the same role for the spatial axis: it mints the modal moment axis
with its own mass diagonal, so the trial space and the sweep kernels
cannot spell the mass twice.

Both are **stage-2 generators**: they induce structure on the space and
on the operator together, consistency between the two inductions is the
gate (tightness for a frame, one closure serving apply and solve for a
scheme), and after minting the generator is forgotten, the operators
retaining only the induced data. A mesh and a quadrature are the
degenerate, space-side-only cases. The axis is the forgetful image that
keeps the induced data, and one of them is its **unit**, the coordinates
of the constant function :math:`1` in its own representation
(:attr:`~orpheus.numerics.axis.Axis.unit`), written by whatever mints the
axis: ones from a measure, :math:`e_0` from a basis, the moments
:math:`\Phi^{\mathsf T} w` from a frame that dresses a head
(:ref:`spaces-collapse-pair-frame`, :ref:`spaces-axis-unit`).


Between spaces — retraction, section, restriction, lift, trace
--------------------------------------------------------------

Four typed families move a field between a space and a subspace, a
boundary or a member, each a split pair :math:`(r, e)` with
:math:`r \circ e = \mathrm{id}`. Along an axis, the **retraction**
:math:`R = \pi_{*} = m^{\mathsf T}\otimes\mathrm{id}` is integration
along the axis, the pairing with the covector :math:`m = G_k u` built from
the axis's unit and metric block (angular integration on the angular axis,
the volume integral on the spatial one, the integral of the represented
function on a modal head; :eq:`spaces-axis-marginal-integral`), and the
**section** :math:`E = (u/\langle m,u\rangle)\otimes\mathrm{id}` is the
mass-normalised constant field (:ref:`spaces-collapse-pair-naming`;
"embedding" is not an operator name in this corpus, by that ruling);
:class:`~orpheus.numerics.axis_marginal.AxisRetractionOperator`,
:class:`~orpheus.numerics.axis_marginal.AxisSectionOperator`, built by
``.of(space, label)``. The retraction's Hilbert adjoint is a third arrow
and a third type, the **pullback** :math:`R^{\dagger} = \pi^{*} =
u\otimes\mathrm{id}` (on a nodal axis the plain broadcast), minted with
the retraction and returned by ``R.H``
(:class:`~orpheus.numerics.axis_marginal.AxisPullbackOperator`;
:ref:`spaces-collapse-pair-pullback`); it differs from the section by the
axis's mass :math:`\langle m, u\rangle` (:math:`\Sigma w` on a nodal
axis), which is why an angle-integrated source
enters phase space through :math:`E` and a detector through
:math:`R^{\dagger}` (:ref:`spaces-collapse-pair-two-lifts`). Along an index
subset, the **trace restriction** :math:`\gamma_S` gathers and its
transpose :math:`\iota_S` scatters (:eq:`bc-trace-restriction-pair`,
:ref:`bc-domain-narrowing`;
:class:`~orpheus.numerics.elementary_operators.TraceRestrictionOperator`); the
boundary trace is one whole-boundary space whose inflow and outflow
halves are selectors over the sign of :math:`\Omega \cdot \hat n`
(:ref:`bc-trace-structure`;
:class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`),
and one face's inflow or outflow half is a space of its own, a
:ref:`half-trace <bc-half-trace>`
(:class:`~orpheus.numerics.spaces.angular_trace_space.AngularFaceTraceSpace`).
Along a system, the **system restriction** selects a member of a
coupled field and its transpose is extension by zero
(:ref:`coupled-block-system-restriction`,
:eq:`coupled-block-system-restriction-pair`;
:class:`~orpheus.numerics.coupled_system.SystemRestrictionOperator`), on
the coupled space, the direct sum of the members
(:ref:`coupled-block-n-general-machinery`;
:class:`~orpheus.numerics.coupled_system.CoupledSpace`). Into the composite, the **lift**
carries a bulk action onto :math:`\text{bulk} \oplus \text{trace}` by
extension by zero on the trace (:ref:`cs4c-ends-select-the-body`;
:class:`~orpheus.transport.operators.lift.BulkLift`;
:class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`).

A face of the boundary has one **role**: a *deck*, which says the domain
continues across the face (a rigid motion that glues it to a partner face,
or to itself for a mirror), or an integer *tag*, which says the domain ends
there and which only names a colour of the true boundary, the exact twin of
a material id (:class:`~orpheus.geometry.deck.MirrorDeck`,
:class:`~orpheus.geometry.deck.TranslationDeck`, :ref:`geometry-face-roles`).
The geometry says which faces share an exterior model and never what the
model is. A **wall condition** is what a tag is bound to: the affine pair
of a response network :math:`R` and an inflow :math:`q`
(:class:`~orpheus.transport.boundary_conditions.WallCondition`,
:eq:`affine-bc-form`), a descriptor and not an operator. :math:`R = sP +
d\,(E\circ Rt)` is a point :math:`(s, d)` of the convex cone
:math:`s, d \ge 0,\ s + d \le 1`, so vacuum, white, specular and a
prescribed inflow are facades over one type, and :math:`R = 0` is one value
however it is reached; a combination that leaves the cone is a
:class:`~orpheus.transport.boundary_conditions.WallDisplacement`, which a
problem cannot bind. A **binding**
(:func:`~orpheus.transport.boundary_conditions.bind`) joins the declared
tag-to-condition map to the faces once, as an exact cover of the tags, and
is lawful only when every condition is invariant under the **colouring
stabiliser** of the problem, the group of isometries that map the material
runs onto themselves and each face role onto the same tag or the same deck
type (:func:`~orpheus.geometry.colouring.colouring_symmetry`, computed from
the material runs, so refining a mesh cannot change it: the three producers'
raw partitions, a geometry's intervals, a mesh's regions and a material
mesh's cells, are merged into runs at one site). A method meets the binding
through its **realizer table**, one per method, which answers admission (is
the declaration realized exactly, approximated, or refused on a named axis;
a plumbing refusal names the issue that tracks the missing machinery) and realization together,
through one door, so a query and a realization cannot disagree
(:ref:`bc-method-realizability`;
:class:`~orpheus.transport.boundary_conditions.RealizerTable`). This is why
the sweep reads the inflow as a given and never re-applies the condition.


The metric is an object, and the adjoint falls out
--------------------------------------------------

A space resolves its metric into one object
(:ref:`spaces-metric-object`;
:class:`~orpheus.numerics.metric.HilbertMetric`, factored per weighted
axis or dense on a Gram head), with the Moore–Penrose inverse on the
kernel that the tangential trace slots occupy, and one spelling of the
pairing. The **Riesz legs** are first-class arrows: ♭ lowers,
:math:`V \to V^{*}`, applying :math:`G`
(:class:`~orpheus.numerics.operator.RieszLowerOperator`); ♯ raises,
:math:`V^{*} \to V`, applying :math:`G^{+}`
(:class:`~orpheus.numerics.operator.RieszRaiseOperator`); they are
defined beside the metric object (:ref:`spaces-riesz-legs`,
:eq:`spaces-riesz-lower-raise`), and the SN development history records
their landing. The **adjoint** (:ref:`g-adjoint`,
:eq:`g-adjoint-definition`, :eq:`spaces-adjoint-riesz-composition`) is
the three-factor composition
:math:`A^{*} = \sharp_V \circ A^{\mathsf T} \circ \flat_W`, built at
construction from the bound spaces
(:class:`~orpheus.numerics.operator.AdjointOperator`, reached as ``A.H``);
the metric-free transpose is ``A.dual()`` between the dual spaces
(:ref:`spaces-metric-propagation`), and ``apply_transpose`` is its raw
array verb. There is no ``.T`` property anywhere in the tree. Composers
propagate the transpose by law, sums, products, tensor products
(:ref:`tensor-product-spaces`) and block operators, and
``(A.H).inverse()`` is ``A.inverse().H`` as an object identity, so the
adjoint sweep is the reverse scan reached through the sweep operator's
transpose, and no driver module spells a transpose (0 of 24
``apply_transpose`` call sites, 2026-09-21). The adjoint needs both
ends, which is why it exists only for a **bound operator**, one
constructed with its domain and codomain declared
(:ref:`bound-operator`), with one exemption: the constant multipliers
:math:`cI`, which commute with every metric, so their transpose is their
adjoint on any space. The dual space :math:`V^{*}` carries the inverse
Riesz form :math:`G^{+}` and is distinct from :math:`V`
(:ref:`spaces-dual-space`), and a degenerate metric is an inner product on
the quotient :math:`V/\ker G` (:ref:`spaces-quotient`).


Symmetry and quotient
---------------------

A **symmetry group** is realised, never tabulated
(:ref:`manifold-realization`;
:class:`~orpheus.geometry.symmetry.SubgroupOfO3` with its computed
:class:`~orpheus.geometry.symmetry.Realization`, the identity component
plus one representative per component): containment, the normaliser and
the connected part's fixed set are each one computation. Whether a
measure is invariant is the measure's question
(:eq:`discrete-measure-g-invariance`;
:meth:`~orpheus.numerics.measure.DiscreteMeasure.is_invariant_under`),
asked on the measure's own support, so that a folded rule's permutation
and its invariance cannot disagree. The support must be a **G-space**,
one of the families a rigid motion acts on, and the action is one law:
:math:`g` acts iff :math:`g(M) = M` (:ref:`manifold-gspace`,
:eq:`manifold-gspace-action-law`;
:class:`~orpheus.geometry.manifold.GSpace`); energy groups, an index set
and a product are not G-spaces and refuse the question. Every invariance
verb is a composition of the measure's pushforward with that action
through one group-free relabelling (:eq:`manifold-invariance-composition`,
:ref:`discrete-measure-invariance-module`;
:meth:`~orpheus.numerics.measure.DiscreteMeasure.relabelling_to`). An **orbit space** is named by its stabiliser, one spelling
(:ref:`manifold-orbit-space`, :ref:`manifold-orbit-space-stabiliser`,
:ref:`manifold-quotient-map`; :class:`~orpheus.geometry.manifold.Quotient`,
:meth:`~orpheus.numerics.measure.DiscreteMeasure.quotient`), and the lift
from the orbit space back to the ambient measure is the Reynolds projector
:math:`P_H`, the orbit barycentre, which is not a section
(:ref:`manifold-lift`, :eq:`manifold-reynolds-projector`,
:ref:`manifold-reynolds-projector-section`). **Descent** pulls a basis back along the
quotient map (:ref:`manifold-descent`;
:class:`~orpheus.numerics.basis.descent.Descent`). The group's exact role
in posing: the orbit partition bounds how coarse an admissible pose may
be, and among admissible refinements the good ones respect it; refinement
is the flow, symmetry the admissibility bound.

The same construction acts on **positions**. A 1-D coordinate system is a
**chart**: its coordinate :math:`c` (the slab's :math:`x`, the distance to
the cylinder's axis or the sphere's centre) is the quotient map by its
symmetry group :math:`G_c = \{g \in E(3) : c \circ g = c\}`, whose
generic isotropy is the angular symmetry the problem spends and whose
singular strata are the centre and the axis
(:ref:`chart-and-chord-chart`; :class:`~orpheus.geometry.chart.Chart`).
A straight line is a **line** in Plücker coordinates
(:class:`~orpheus.geometry.line.Line`), and its **chord** through the
level sets of :math:`c` is solved once in the orbit space, every 3-D
length being an orbit-space length times the obliquity
:math:`1/|P\Omega|` (:ref:`chart-and-chord-chord`;
:class:`~orpheus.geometry.chord.ConcentricPartition`).


The Solution half — splitting, resolvent, iteration, outcome
------------------------------------------------------------

The Problem poses; the Strategy labels; the driver applies. The
**Strategy** is a value and not the Problem's
(:ref:`sn-splitting-is-a-strategy-value`;
:class:`~orpheus.sn.splitting.Splitting`, minted at one site from a
schedule): it labels each loss term implicit or explicit, so
:math:`A = M - N` with :math:`M` and :math:`N` derived through the
algebra's own sums, and it certifies itself against the Problem by the
law :math:`M - N = A` (this :math:`M` is the splitting's implicit part;
the pencil's second member is written :math:`F` above). The Strategy
inverts only the implicit part: ``implicit.inverse()`` is the sweep
(:class:`~orpheus.sn.operators.sweep_operator.SweepOperator`) or the block
back-substitution on a carrying mesh
(:class:`~orpheus.numerics.coupled_system.CoupledSubstitutionOperator`),
and the driver only applies it. The inverse of the split itself,
:math:`(M - N)^{-1}`, is minted by the Strategy as a
:class:`~orpheus.numerics.green_operator.GreenOperator` with the
realisation it chose (``Splitting.green``; :ref:`green-strategy-mints`):
its relaxation face builds the driver, so
:class:`~orpheus.numerics.iteration.SourceIteration` runs
:math:`\psi \leftarrow M^{-1}(q + N\psi)` with the residual stop under a
:class:`~orpheus.numerics.green_operator.Richardson` realisation, and
:class:`~orpheus.numerics.iteration.KrylovAcceleration` hands the matvec to
GMRES under a :class:`~orpheus.numerics.green_operator.Krylov` one. The eigen resolvent :math:`A^{-1}F` of :eq:`eigen-resolvent` is
never formed; the iteration realises it. One outer loop serves every
iterative family, :func:`~orpheus.numerics.eigenvalue.power_iteration` over
the :class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` Protocol; the
0-D baseline solves its pencil directly.

Every level records itself and nothing stores a verdict: each driver
mints an :class:`~orpheus.numerics.convergence.IterationRecord` whose
``converged`` is derived from co-indexed criteria and a budget. The
convergence contract is two-sided: a best-effort exit is legal and made
audible once at the public entry
(:func:`~orpheus.numerics.convergence.warn_if_unconverged`), and a claimed
convergence is re-measured through a real forward apply and refused
beyond its tolerance. The answer is fused with its question and its gauge
into an **outcome** (:class:`~orpheus.numerics.outcome.EigenOutcome`,
:class:`~orpheus.numerics.outcome.SourceOutcome`; the gauge is the
section that picked the representative,
:class:`~orpheus.numerics.gauge.ScaleGauge`), reported on
(:class:`~orpheus.numerics.outcome.ExitReport`), and packaged once.

The **Solution** (:ref:`the-solution-outcome`;
:ref:`sn-solution-carries-its-posing`) is the frozen five-tuple
*problem, outcome, strategy, exit report, record*
(:class:`~orpheus.sn.solution.SolutionBase`). Its kind, eigen or source,
is the outcome's type parameter; its role is the class:
:class:`~orpheus.sn.solution.Solution` carries the forward verbs
(reaction rates, homogenisation, condensation) and
:class:`~orpheus.sn.solution.AdjointSolution` carries importance and
structurally lacks them. The flux members are derived readers of the
outcome's state. The adjoint entries dagger the same objects, the
implicit part, the gains and the posing, through ``.H``
(:ref:`sn-adjoint`), which is how the adjoint falls out of the algebra
rather than out of a second solver.


The concept table
-----------------

.. list-table::
   :header-rows: 1
   :widths: 22 34 44

   * - Concept
     - Class or function
     - Defined at
   * - Space, axis, basis kind, origin, coefficient layout
     - :class:`~orpheus.numerics.space.FunctionSpace` (its :meth:`~orpheus.numerics.space.FunctionSpace.identity_key`, :meth:`~orpheus.numerics.space.FunctionSpace.of_measure`, :meth:`~orpheus.numerics.space.FunctionSpace.of_basis`), :class:`~orpheus.numerics.axis.Axis`, :class:`~orpheus.numerics.axis.BasisKind`, :class:`~orpheus.numerics.axis.Modal`, the origin arms :class:`~orpheus.numerics.axis.FromMeasure`, :class:`~orpheus.numerics.axis.FromQuadrature`, :class:`~orpheus.numerics.axis.FromBasis`, :class:`~orpheus.numerics.axis.FromFrame`, :class:`~orpheus.numerics.axis.Declared`; :mod:`orpheus.numerics.coefficient_layout`
     - :ref:`spaces-the-axis`, :eq:`spaces-axis-product`, :ref:`spaces-axis-origin`, :ref:`spaces-axis-identity`, :ref:`spaces-identity-key`, :ref:`spaces-moment-head`
   * - Discrete measure; manifold; quadrature
     - :class:`~orpheus.numerics.measure.DiscreteMeasure`, :class:`~orpheus.geometry.manifold.Manifold`, :class:`~orpheus.numerics.quadrature.directional.Quadrature`
     - :eq:`discrete-measure-definition`
   * - Basis
     - :class:`~orpheus.numerics.basis.base.Basis`, :class:`~orpheus.numerics.basis.base.GramStructure`. On the reference side of the branch line, which imports none of these (:ref:`architecture-reference-insulation`): :class:`~orpheus.derivations.continuous.characteristic.basis.PanelBasis`, discontinuous nodal Lagrange panels over a body's orbit coordinate, graded toward walls and interfaces and even (polynomials in :math:`c^2`) on the panel touching a singular stratum, with its Gram matrix in the chart's volume measure, whose density is the kernel's :meth:`Chart.measure_density <orpheus.geometry.chart.Chart.measure_density>`, the derivative of :meth:`CoordSystem.measure <orpheus.geometry.coord.CoordSystem.measure>`'s one definition. No production basis is nodal Lagrange, so it has no twin.
     - :ref:`spaces-basis`, :ref:`manifold-three-levels`; :ref:`characteristic-panel-basis`
   * - Cone
     - :meth:`~orpheus.numerics.field.Field.cone_violations`
     - :ref:`cone-ordered-vector-space`, :eq:`positive-cone-definition`, :ref:`cone-membership-is-a-predicate`
   * - Frame; analysis; reconstruction
     - :class:`~orpheus.numerics.frame.FrameBase`, :class:`~orpheus.numerics.frame.PetrovGalerkinFrame`, :class:`~orpheus.numerics.frame.GalerkinFrame`, :class:`~orpheus.transport.frames.harmonic_frame.HarmonicFrame`; :class:`~orpheus.numerics.projection.AnalysisOperator`, :class:`~orpheus.numerics.projection.ReconstructionOperator`
     - :ref:`galerkin-projection`, :eq:`galerkin-pair`, :eq:`galerkin-frame-idempotency`, :eq:`galerkin-construction`, :eq:`petrov-galerkin-construction`
   * - Gram; Parseval metric
     - :attr:`~orpheus.numerics.frame.FrameBase.discrete_gram`
     - :ref:`frame-analysis-is-the-gram-section`, :eq:`frame-discrete-gram`, :ref:`frame-parseval-metric`
   * - Scheme
     - :class:`~orpheus.transport.spatial.scheme.DiscretizationSchemeBase`
     - :ref:`discretization-closures`
   * - Unit of an axis (the constant function's coordinates); covector of integration
     - :attr:`~orpheus.numerics.axis.Axis.unit`, :attr:`~orpheus.numerics.axis_marginal.AxisRetractionOperator.covector`
     - :ref:`spaces-axis-unit`, :eq:`spaces-axis-marginal-integral`
   * - Retraction; section; pullback (the retraction's adjoint)
     - :class:`~orpheus.numerics.axis_marginal.AxisRetractionOperator`, :class:`~orpheus.numerics.axis_marginal.AxisSectionOperator`, :class:`~orpheus.numerics.axis_marginal.AxisPullbackOperator`
     - :ref:`spaces-collapse-pair`, :ref:`spaces-collapse-pair-naming`, :ref:`spaces-collapse-pair-pullback`, :ref:`spaces-collapse-pair-two-lifts`
   * - Trace restriction; trace space; half-trace
     - :class:`~orpheus.numerics.elementary_operators.TraceRestrictionOperator`, :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`, :class:`~orpheus.numerics.spaces.angular_trace_space.AngularFaceTraceSpace`
     - :eq:`bc-trace-restriction-pair`, :ref:`bc-trace-structure`, :ref:`half-trace <bc-half-trace>`
   * - System restriction; coupled space
     - :class:`~orpheus.numerics.coupled_system.SystemRestrictionOperator`, :class:`~orpheus.numerics.coupled_system.CoupledSpace`
     - :ref:`coupled-block-system-restriction`, :eq:`coupled-block-system-restriction-pair`, :eq:`coupled-block-system-restriction-laws`; the coupled space: :ref:`coupled-block-n-general-machinery`
   * - Lift; full-field space
     - :class:`~orpheus.transport.operators.lift.BulkLift`, :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
     - :ref:`cs4c-ends-select-the-body`
   * - Face role; deck; boundary tag
     - ``FaceRole`` (a deck or an ``int``), :class:`~orpheus.geometry.deck.MirrorDeck`, :class:`~orpheus.geometry.deck.TranslationDeck`, :func:`~orpheus.geometry.deck.pairing_of` (the :class:`~orpheus.geometry.deck.FacePairing` a realizer reads), :func:`~orpheus.geometry.deck.parse_face_roles`; the value both meshes store, ``orpheus.mesh.face_roles.FaceRoles`` (``Mesh1D.face_roles``, ``outer_role``)
     - :ref:`geometry-face-roles`, :ref:`bc-architecture-of-record`
   * - Wall condition on the cone; response network; inflow; displacement
     - :class:`~orpheus.transport.boundary_conditions.WallCondition` (the facades ``vacuum``, ``specular``, ``white``, ``uniform``, ``framed`` and ``from_amplitudes``), :class:`~orpheus.transport.boundary_conditions.ResponseNetwork`, :class:`~orpheus.transport.boundary_conditions.InflowTerms`, :class:`~orpheus.transport.boundary_conditions.WallDisplacement`; the method-only conditions share the base :class:`~orpheus.transport.boundary_conditions.BoundaryCondition` (zero flux is ``orpheus.diffusion.zero_flux.ZeroFluxCondition``)
     - :eq:`affine-bc-form`, :ref:`bc-architecture-of-record`
   * - Binding; declared conditions
     - :class:`~orpheus.transport.boundary_conditions.BoundaryConditions` (the twin of ``Materials``), :class:`~orpheus.transport.boundary_conditions.Binding`, :func:`~orpheus.transport.boundary_conditions.bind`
     - :ref:`bc-architecture-of-record`
   * - Colouring stabiliser
     - :func:`~orpheus.geometry.colouring.colouring_symmetry`, :class:`~orpheus.geometry.colouring.ColouringSymmetry`, :class:`~orpheus.geometry.colouring.UncomputedSymmetries` (the typed pessimistic answer of a carrier with no group machinery)
     - :ref:`bc-architecture-of-record`
   * - Admission; realizer table; realization
     - :class:`~orpheus.transport.boundary_conditions.RealizerTable`, :class:`~orpheus.transport.boundary_conditions.Admission` (:class:`~orpheus.transport.boundary_conditions.Realized`, :class:`~orpheus.transport.boundary_conditions.Approximated`, :class:`~orpheus.transport.boundary_conditions.Refused`), :func:`~orpheus.transport.method.realize_boundary`; the tables ``SN_REALIZERS``, ``DIFFUSION_REALIZERS``, ``CP_REALIZERS``, ``MOC_REALIZERS``, ``MC_REALIZERS``; the typed SN operators ``RealizedBoundary`` and its arms ``ZeroMap``, ``ScaledOrdinatePermutation``, ``RankOneReemission``, ``RealizedSum`` (``orpheus.sn.boundary.realized_operators``); the certification :func:`~orpheus.numerics.quadrature.ordinate_pairing.ordinate_pairing`
     - :ref:`bc-realizer-layer`, :ref:`bc-method-realizability`, :ref:`bc-refusal-axes`
   * - Retired 2026-10-24 (wave-2 carve): boundary law; realizer
     - ``BoundaryTraceLaw`` and its concretes, ``LawSum`` and ``LawScaled``, the ``BC`` tag, ``SNBoundaryRealizer`` and the ``BoundaryRealizer`` Protocol, the string registries, ``FaceLaws`` and ``face_laws``. Successors: the rows above.
     - :ref:`bc-wave2-crosswalk`
   * - Metric
     - :class:`~orpheus.numerics.metric.HilbertMetric`
     - :ref:`spaces-metric-object`
   * - Riesz legs
     - :class:`~orpheus.numerics.operator.RieszLowerOperator`, :class:`~orpheus.numerics.operator.RieszRaiseOperator`
     - :ref:`spaces-riesz-legs`, :eq:`spaces-riesz-lower-raise`, :eq:`spaces-riesz-round-trip`
   * - Bound operator
     - :class:`~orpheus.numerics.operator.LinearOperator` (its ``domain`` and ``codomain``), :class:`~orpheus.transport.operators.bound_operator.BoundOperator`
     - :ref:`bound-operator`
   * - Adjoint; dual
     - :class:`~orpheus.numerics.operator.AdjointOperator`, :meth:`~orpheus.numerics.operator.LinearOperator.dual`
     - :ref:`g-adjoint`, :eq:`g-adjoint-definition`, :eq:`spaces-adjoint-riesz-composition`, :ref:`spaces-metric-propagation`
   * - Symmetry group; invariance
     - :class:`~orpheus.geometry.symmetry.SubgroupOfO3`; :class:`~orpheus.geometry.manifold.GSpace` (its :meth:`~orpheus.geometry.manifold.GSpace.action`); :meth:`~orpheus.numerics.measure.DiscreteMeasure.relabelling_to` and the invariance verbs on :class:`~orpheus.numerics.measure.DiscreteMeasure`
     - :ref:`manifold-realization`, :eq:`discrete-measure-g-invariance`, :ref:`manifold-gspace`, :eq:`manifold-invariance-composition`
   * - Orbit space; quotient; descent
     - :class:`~orpheus.geometry.manifold.Quotient`, :class:`~orpheus.numerics.basis.descent.Descent`
     - :ref:`manifold-orbit-space`, :ref:`manifold-orbit-space-stabiliser`, :ref:`manifold-quotient-map`, :ref:`manifold-reynolds-projector-section`, :ref:`manifold-descent`
   * - Chart (a coordinate system's orbit map and symmetry group); singular stratum; the measure's density (the area of a level set)
     - :class:`~orpheus.geometry.chart.Chart` (derived from its kept columns and its linear group), :class:`~orpheus.geometry.chart.SingularStratum`; a line's image :class:`~orpheus.geometry.chart.RadialImage`, :class:`~orpheus.geometry.chart.AxialImage`; :meth:`CoordSystem.measure_density <orpheus.geometry.coord.CoordSystem.measure_density>` (read through :meth:`Chart.measure_density <orpheus.geometry.chart.Chart.measure_density>`)
     - :ref:`chart-and-chord-chart`, :eq:`geometry-radial-coordinate`, :ref:`chart-and-chord-isotropy`, :ref:`chart-and-chord-strata`, :eq:`geometry-measure-density`
   * - Line (Plücker coordinates); ray as a line with a start
     - :class:`~orpheus.geometry.line.Line`
     - :ref:`chart-and-chord-lines`
   * - Concentric partition; chord; crossings; measure on lines; the orbit space of oriented lines (the line domain)
     - :class:`~orpheus.geometry.chord.ConcentricPartition`, :class:`~orpheus.geometry.chord.Chord`, :class:`~orpheus.geometry.chord.Crossings`; :meth:`~orpheus.geometry.chart.Chart.beam_density`; :class:`~orpheus.geometry.chart.LineDomain` (:meth:`Chart.line_domain <orpheus.geometry.chart.Chart.line_domain>`)
     - :ref:`chart-and-chord-chord`, :eq:`geometry-line-crossing-law`, :eq:`geometry-crossing-order`, :eq:`geometry-chord-segment-lengths`, :eq:`geometry-cylinder-axial-factor`, :ref:`chart-and-chord-location`, :eq:`geometry-measure-on-lines`, :eq:`geometry-cauchy-mean-chord`, :eq:`geometry-line-domain`
   * - Mesh-free function; angular chart
     - :class:`~orpheus.numerics.mesh_free_function.RegionwiseConstant`, :class:`~orpheus.numerics.mesh_free_function.Symbolic`; :class:`~orpheus.geometry.coord.AngularChart` (:attr:`CoordSystem.angular_chart <orpheus.geometry.coord.CoordSystem.angular_chart>`)
     - :ref:`structured-geometry-mesh-free-functions`, :ref:`structured-geometry-angular-chart`, :ref:`structured-geometry-two-lifts-branch-1`
   * - Material mesh; method hubs
     - :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`, :class:`~orpheus.sn.problem.SNProblem`, :class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh`, :class:`~orpheus.homogeneous.solver.HomogeneousProblem`
     - :ref:`architecture-layering`; hub: :ref:`the-problem-hub`
   * - Pencil; posings
     - :class:`~orpheus.numerics.pencil.OperatorPencil`, :class:`~orpheus.numerics.posing.EigenPosing`, :class:`~orpheus.numerics.posing.SourcePosing`
     - :ref:`the-operator-pencil`, :eq:`pencil-family`, :ref:`eigenvalue-posing`, :ref:`sn-the-problem-poses-its-pencil`
   * - Reference kernel: dense pencil in weak form, fundamental mode, least solution; and its deliberate twin across the branch line
     - :class:`~orpheus.derivations.common.dense_pencil.DensePencil` (:meth:`~orpheus.derivations.common.dense_pencil.DensePencil.fundamental`, :meth:`~orpheus.derivations.common.dense_pencil.DensePencil.spectrum`, :meth:`~orpheus.derivations.common.dense_pencil.DensePencil.adjoint`, :meth:`~orpheus.derivations.common.dense_pencil.DensePencil.least_solution`, :meth:`~orpheus.derivations.common.dense_pencil.DensePencil.reach`), :func:`~orpheus.derivations.common.quadrature.composite_gauss_legendre`. **Twins, never to be merged:** :meth:`DensePencil.fundamental <orpheus.derivations.common.dense_pencil.DensePencil.fundamental>` (the references' Perron–Frobenius extraction: refuses a complex, non-positive or not strictly dominant eigenvalue and a sign-changing vector) and :func:`~orpheus.numerics.eigenvalue.dominant_eigenpair` (production's: takes the largest real part and refuses only a complex one; its docstring claims the full contract, `#580 <https://github.com/deOliveira-R/ORPHEUS/issues/580>`_). They differ in strength by design, and they exist twice so that a reference agreeing with production shares no extraction code with it.
     - :ref:`verification-reference-kernel`, :ref:`architecture-reference-insulation`
   * - Face role to wall, in the closed references; and the one translation to production, test-side
     - :class:`~orpheus.specification.exterior.WallAmplitudes` (a wall's cone coordinates, with the vertices ``VACUUM``, ``MIRROR``, ``WHITE``), :func:`~orpheus.specification.exterior.admit_exterior` (the exterior as an exact cover of the geometry's tags, admitted once and called by every reference that holds the pair), :func:`~orpheus.derivations.common.reference_body.amplitudes_at` (the one reader of a face role's amplitudes; read by :meth:`Walls.of <orpheus.derivations.continuous.characteristic.walls.Walls.of>` and by :func:`~orpheus.derivations.common.reference_body.specular_albedos`). The reference imports no production response and no typed law (the user's ruling of 2026-10-06). The only translation from a reference's amplitudes to production's conditions is the test-side loop ``tests/gates/_reference_pairing.py``, built on the primitives-taking :meth:`WallCondition.from_amplitudes <orpheus.transport.boundary_conditions.WallCondition.from_amplitudes>`; the string-tag registry and its twin parse ``_law_from_tag`` are retired with the typed law family, so no pair of twins remains across the branch line at this row.
     - :ref:`characteristic-walls-roles`, :ref:`architecture-reference-insulation`
   * - Question; mode; point
     - :class:`~orpheus.numerics.question.Eigen`, :class:`~orpheus.numerics.question.FixedSource`, :class:`~orpheus.numerics.question.Response`; :class:`~orpheus.numerics.question.Fundamental`, :class:`~orpheus.numerics.question.Nearest`; the point a :class:`~orpheus.values.content.FrozenMapping`
     - :ref:`structured-geometry-question-values`, :ref:`structured-geometry-question-values-point`, :ref:`structured-geometry-question-values-role`
   * - Reference specification; coordinate; channel
     - :class:`~orpheus.specification.specification.InfiniteMediumSpecification`, :class:`~orpheus.specification.specification.GeometrySpecification`; :class:`~orpheus.data.cells.CellCoefficient`, :class:`~orpheus.geometry.extent.GeometryExtent`; :class:`~orpheus.data.cells.Channel`
     - :ref:`structured-geometry-specification`, :ref:`structured-geometry-specification-layer`, :ref:`structured-geometry-specification-canonical`, :ref:`structured-geometry-specification-coordinates`
   * - Strategy; splitting; schedule
     - :class:`~orpheus.sn.splitting.Splitting`, :func:`~orpheus.sn.splitting.resolve_schedule`
     - :ref:`sn-splitting-is-a-strategy-value`
   * - Resolvent
     - :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`, :class:`~orpheus.numerics.coupled_system.CoupledSubstitutionOperator`
     - :eq:`eigen-resolvent`
   * - Green of a split; realisation; representation
     - :class:`~orpheus.numerics.green_operator.GreenOperator`, :class:`~orpheus.numerics.green_operator.Richardson`, :class:`~orpheus.numerics.green_operator.Krylov`, :class:`~orpheus.numerics.iteration.Representation`
     - :ref:`green-strategy-mints`, :ref:`green-realisations`
   * - Iteration; outer loop
     - :class:`~orpheus.numerics.iteration.SourceIteration`, :class:`~orpheus.numerics.iteration.KrylovAcceleration`, :func:`~orpheus.numerics.eigenvalue.power_iteration`
     - :ref:`eigenvalue-posing`
   * - Record; exit report; gauge
     - :class:`~orpheus.numerics.convergence.IterationRecord`, :class:`~orpheus.numerics.outcome.ExitReport`, :class:`~orpheus.numerics.gauge.ScaleGauge`
     - :ref:`the-solution-outcome`
   * - Outcome; Solution
     - :class:`~orpheus.numerics.outcome.EigenOutcome`, :class:`~orpheus.numerics.outcome.SourceOutcome`, :class:`~orpheus.sn.solution.SolutionBase`, :class:`~orpheus.sn.solution.Solution`, :class:`~orpheus.sn.solution.AdjointSolution`
     - :ref:`the-solution-outcome`, :ref:`sn-solution-carries-its-posing`


What is still hand-rolled — the debt list
-----------------------------------------

Anything hand-rolled in production is a debt that needs refinement; a
guard is the same debt in another spelling. Measured 2026-09-21 on
``b4865b5e``; the commands are in the memos named at the top.

- **No shared Problem type.** Three hubs (:ref:`the-problem-hub`) realise the concept
  (:class:`~orpheus.sn.problem.SNProblem`,
  :class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh`,
  :class:`~orpheus.homogeneous.solver.HomogeneousProblem`; only two are
  named ``*Problem``) and no ``Problem`` ABC or Protocol exists; :class:`~orpheus.transport.method.TransportMethod`
  covers the two method-meshes only. The Problem → Solution carve, a
  standalone Problem module with a thin solver per family, is the
  consumers campaign's, and :class:`~orpheus.homogeneous.solver.HomogeneousProblem`
  lives in the solver module until then.
- **The diffusion hub mints neither pencil nor posing**: ``def pencil``
  and ``def eigen_posing`` exist on 2 of 3 hubs, ``def source_posing`` on
  1; the diffusion solver assembles its loss, its fission and an exact LU
  resolvent itself. Its result, like the CP and MoC results, is not on
  the Solution shape (per-family fields: Diffusion, CP and MoC carry a
  record only, Homogeneous an outcome only).
- **The α posing is stated and not minted**:
  :data:`~orpheus.numerics.posing.ALPHA_MAP` is defined and no hub mints an
  α :class:`~orpheus.numerics.posing.EigenPosing`.
- **The inner algorithm is a string.** Source iteration versus Krylov is
  dispatched on ``inner_solver: str`` by the SN coordinator, and the
  Solution records it only as the record's label; the schedule string
  survives on the entry signatures, resolved at one site. The SN outer
  iterate is a bare array; the 2-D angular window and the DSA corrector
  travel outside the Strategy.
- **The adjoint fixed-source posing is spelled by hand, and no iteration
  is derived from its posing** (#484). ``solve_sn_adjoint_fixed_source``
  names its operator directly (``pencil.H.lhs``) where
  :meth:`~orpheus.numerics.posing.SourcePosing.H` would dagger the
  forward posing, so the σ = 0 member is spelled twice; every source
  entry builds its iteration from the Splitting with the σ = 1 production
  bolted on as an extra gain while the posing is only recorded, and the
  convergence-claim check reads the Splitting's residual, so a
  posing-versus-iteration mismatch is a diagnostic, never a refusal; the
  adjoint Strategy is a separate Jacobi mint daggered piecewise. That is
  why a multiplying-source adjoint, ``source_posing(q).H(q*)`` with
  ``F.H``, does not fall out today. The eigen pair is tied
  (``eigen_posing.H()``, ``k† = k`` pinned) and the pure-transport pair is
  witnessed by the duality gate.
- **Adjoints the eager gate refuses**: the boundary Gauss–Seidel reverse
  scan, the linear-discontinuous transpose kernel, and 11 of 67 concrete
  :class:`~orpheus.numerics.operator.LinearOperator` subclasses with no
  ``apply_transpose`` (48 own one, 8 inherit one); there is no diffusion
  adjoint entry.
- **The angular measure's mass is still typed by hand in three places,
  and the retraction's adjoint has a second spelling** (`[M]` 2026-10-02,
  #405 P1 step 6's census and review): MoC lifts an isotropic source by
  a hand-written :math:`1/(4\pi)` instead of the angular section (#556);
  the manufactured-solution builders divide by a hand-read weight sum
  (12 of 12 S\ :sub:`N` cases) and the Green's-function sphere by a typed
  :math:`4\pi`, where one Branch-1 measure object should derive it
  (#557); a composite holding a retraction daggers through the generic
  sandwich rather than the leaf pullback (#558), and
  ``ScalarSourceSink.__add__`` adds a hand-written broadcast,
  :math:`\pi^{*}` outside its type.
- **Chords, point location and the measure on lines are spelled outside
  their kernel** (`[M]` 2026-10-05, the geometry census at ``a336bde4``,
  :ref:`chart-and-chord-deferred`): 13 chord square roots and 14
  geometric discriminants outside the SymPy origins, 8 point-location
  spellings with two boundary conventions, and 3 realisations of the
  measure on lines, none calling
  :class:`~orpheus.geometry.chord.ConcentricPartition`. Its one consumer
  is the characteristic reference (:ref:`theory-characteristic-reference`),
  which chords its lines through a body's partition and through the
  panel partition refining it, and is built to replace the trajectory
  resolvents that hold those spellings (#405).
- **The measure's density is spelled twice in the input layer**:
  :meth:`CoordSystem.measure_density <orpheus.geometry.coord.CoordSystem.measure_density>`,
  the derivative of the one measure, and production's
  ``compute_areas_1d``, the face areas of a 1-D mesh, which are the same
  level-set areas. Retiring the second onto the first moves bits on the
  sphere in S\ :sub:`N`, so it is its own step, `#584
  <https://github.com/deOliveira-R/ORPHEUS/issues/584>`_
  (:ref:`chart-and-chord-measure-density`).
- **Exit-report members still listed as not yet built** (the outcome
  module's ``NotYet``):
  the carrying eigen exit's balance (#354, open) and the daggered eigen
  exit (#353, open); the list also names the linear-discontinuous residual
  (#310), whose issue is closed.
- **Restriction siblings share no base** by ruling ("no consumer treats
  any restriction generically yet"); the rank-d spatial axis is
  ``Declared`` and carries no atoms (a gated contract, G6b); there is no ``Cone`` class, no
  embedding operator and no affine operator, each by ruling and so not
  debt.
- **The tagged guards**: ``ELEGANCE-DEBT[guard]`` occurs 9 times in
  9 sites under ``orpheus/`` (`[M]` 2026-10-24 at ``471c4737``, a Python
  scan of every line of ``orpheus/**/*.py`` for the tag; positive control:
  the scan finds the ``#457`` site in ``orpheus/transport/full_field.py``).
  The sites: the full-field carrier (#457); the lock on withdrawn
  reference generators in ``orpheus/derivations/common/withdrawal.py``
  (#506, retired when each reference family's generator returns a
  :class:`~orpheus.reference.certificate.ReferenceCertificate` whose
  standing carries the
  :class:`~orpheus.reference.withdrawal.Withdrawal`;
  :ref:`vv-withdrawn-generators`); the Sood registry's refusal of a case
  with no citation (#405, retired when a published solution replaces the
  registry's case class); the 2-D mesh's cell integrals, refused where the code
  would go (#569, retired when ``Mesh2D`` carries region labels and the
  mesh-free function owns its second coordinate); the characteristic
  reference's refusal of a traversed slot whose 3-D length overflowed
  (#582, two sites, retired when the kernel's chord refuses or resolves a line within
  an underflow of parallel to the level sets) and its refusal of a slab
  direction within the grazing bound, whose crossing parameters are
  differences of large absolute positions (#585, retired when a slab line's
  crossings are formed relative to a point inside the body); and the measure's
  G-space guards ``DiscreteMeasure._gspace`` and ``_acted_on_in`` (#624,
  two sites, retired when a measure is typed by its support,
  :ref:`manifold-gspace`). The boundary composition's check on its direct
  children (#551) retired with ``LawSum`` and ``LawScaled`` at the wave-2
  carve: a combination of wall conditions is arithmetic on a cone. A wall condition
  does not add to a non-condition (its ``__add__`` returns ``NotImplemented``,
  so ``w + 3`` is Python's ``TypeError``), and an operator operand is refused
  by the operand check of ``OperatorSum`` in ``numerics/operator.py``, which
  raises ``TypeError``.
  ``# TODO`` occurs once; ``raise NotImplementedError`` 72
  times in 27 files under ``orpheus/`` excluding ``derivations/`` (119 in
  39 with it), the population a retirement audit walks (``[M]``
  2026-10-03, a line count over ``orpheus/**/*.py``).
