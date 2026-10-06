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

The tree it describes is ``main`` at ``b4865b5e`` (2026-09-21). Every
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
:ref:`spaces-the-axis`); an axis is the value object *(shape, factor
measure, basis kind)* plus a ``generator`` slot that is provenance and
deliberately excluded from identity (:ref:`spaces-axis-generator`), so
metric differences imply space differences. Realising a law or a kernel
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
walked backward: the collapse pair is the projection onto a coarser stage
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
scheme, the angular closure, the scattering order and the boundary laws;
:class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh` adds the scalar
trace and the realised albedo laws. The two are the witnesses of the
:class:`~orpheus.transport.method.TransportMethod` Protocol, whose one
body :func:`~orpheus.transport.method.resolve_boundary_conditions` walks
the faces to a typed law and the method's realizer
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
third type. The affine form itself is the boundary law's,
:eq:`affine-bc-form`, and there is no affine operator by ruling: every
affine law is a linear operator plus a typed source.

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
holds the materials and a geometry with its laws. Either is admitted in a
canonical form: spectators dropped, every key resolved to a coordinate of
the problem, the datum fitted to its groups, regions and readable
coordinates. Nothing binds a question to a system's operators yet: deriving
the pencil and the spectral map from a parameter, and the mode law, are the
posing sequence's unit 6 (#529).


The objects on an axis — measure, basis, frame, cone
----------------------------------------------------

A **discrete measure** :math:`\mu = \sum_i w_i \delta_{x_i}`
(:eq:`discrete-measure-definition`) carries nodes, weights, a typed
support :class:`~orpheus.numerics.manifold.Manifold` and an invariance
group; tensor product, direct sum, pushforward along a typed
:class:`~orpheus.numerics.manifold.ManifoldMap`, restriction, partition
and quotient are its verbs (:class:`~orpheus.numerics.measure.DiscreteMeasure`).
A quadrature wraps one such measure and mints the angular axis and the
angular frame from it (:class:`~orpheus.numerics.quadrature.directional.Quadrature`).
An axis is the measure's forgetful image: the weights are kept, the nodes
dropped, and the axis is nodal; a **basis** of the harmonic family
mints modal axes (:class:`~orpheus.numerics.axis.HarmonicAxis`,
:class:`~orpheus.numerics.axis.LegendreAxis`;
:ref:`spaces-moment-head-axis-built`). The energy axis is a
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
degenerate, space-side-only cases. The discipline is stated at rank one
by the collapse pair (:ref:`spaces-collapse-pair-frame`): a single-region
indicator basis bound to the axis measure is built, read for its induced
retraction and section, and discarded.


Between spaces — retraction, section, restriction, lift, trace
--------------------------------------------------------------

Four typed families move a field between a space and a subspace, a
boundary or a member, each a split pair :math:`(r, e)` with
:math:`r \circ e = \mathrm{id}`. Along an axis, the **retraction**
:math:`R = \pi_{*}` is the measure contraction (angular integration on the
angular axis, the volume integral on the spatial one) and the **section**
:math:`E` is the measure-normalised constant field
(:ref:`spaces-collapse-pair-naming`; "embedding" is not an operator name
in this corpus, by that ruling);
:class:`~orpheus.numerics.operator.AxisRetractionOperator`,
:class:`~orpheus.numerics.operator.AxisSectionOperator`. The retraction's
Hilbert adjoint is a third arrow and a third type, the **pullback**
:math:`R^{\dagger} = \pi^{*}`, the plain broadcast, minted with the pair
and returned by ``R.H``
(:class:`~orpheus.numerics.operator.AxisPullbackOperator`;
:ref:`spaces-collapse-pair-pullback`); it differs from the section by the
axis's mass :math:`\Sigma w`, which is why an angle-integrated source
enters phase space through :math:`E` and a detector through
:math:`R^{\dagger}` (:ref:`spaces-collapse-pair-two-lifts`). Along an index
subset, the **trace restriction** :math:`\gamma_S` gathers and its
transpose :math:`\iota_S` scatters (:eq:`bc-trace-restriction-pair`,
:ref:`bc-domain-narrowing`;
:class:`~orpheus.numerics.operator.TraceRestrictionOperator`); the
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

A **boundary law** is a descriptor, not an operator: it carries the typed
factors :math:`G` and :math:`R` of the affine form
(:class:`~orpheus.geometry.boundary.BoundaryTraceLaw`) and becomes an
operator only through the method's realizer
(:ref:`bc-method-realizability`;
:class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`), which is why
the sweep reads the inflow as a given and never re-applies the law.


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
(:ref:`bound-operator`), with one exemption: the metric-free pointwise
stratum, whose transpose needs no metric.


Symmetry and quotient
---------------------

A **symmetry group** is realised, never tabulated
(:ref:`manifold-realization`;
:class:`~orpheus.numerics.symmetry.SubgroupOfO3` with its computed
:class:`~orpheus.numerics.symmetry.Realization`, the identity component
plus one representative per component): containment, the normaliser and
the connected part's fixed set are each one computation. Whether a
measure is invariant is the measure's question
(:eq:`discrete-measure-g-invariance`;
:meth:`~orpheus.numerics.measure.DiscreteMeasure.is_invariant_under`,
delegating to :mod:`orpheus.numerics.invariance`), asked on the measure's
orbit space so that a folded rule's permutation and its invariance cannot
disagree. An **orbit space** is named by its stabiliser, one spelling
(:ref:`manifold-orbit-space`, :ref:`manifold-orbit-space-stabiliser`,
:ref:`manifold-quotient-map`; :class:`~orpheus.numerics.manifold.Quotient`,
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
and the driver only applies it:
:class:`~orpheus.numerics.iteration.SourceIteration` runs
:math:`\psi \leftarrow M^{-1}(q + N\psi)` with the residual stop, and
:class:`~orpheus.numerics.iteration.KrylovAcceleration` hands the matvec to
GMRES. The eigen resolvent :math:`A^{-1}F` of :eq:`eigen-resolvent` is
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
   * - Space, axis, basis kind
     - :class:`~orpheus.numerics.space.FunctionSpace`, :class:`~orpheus.numerics.axis.Axis`, :class:`~orpheus.numerics.axis.BasisKind`
     - :ref:`spaces-the-axis`, :eq:`spaces-axis-product`, :ref:`spaces-axis-generator`
   * - Discrete measure; manifold; quadrature
     - :class:`~orpheus.numerics.measure.DiscreteMeasure`, :class:`~orpheus.numerics.manifold.Manifold`, :class:`~orpheus.numerics.quadrature.directional.Quadrature`
     - :eq:`discrete-measure-definition`
   * - Basis
     - :class:`~orpheus.numerics.basis.base.Basis`, :class:`~orpheus.numerics.basis.base.GramStructure`
     - :ref:`spaces-basis`, :ref:`manifold-three-levels`
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
   * - Retraction; section; pullback (the retraction's adjoint)
     - :class:`~orpheus.numerics.operator.AxisRetractionOperator`, :class:`~orpheus.numerics.operator.AxisSectionOperator`, :class:`~orpheus.numerics.operator.AxisPullbackOperator`
     - :ref:`spaces-collapse-pair`, :ref:`spaces-collapse-pair-naming`, :ref:`spaces-collapse-pair-pullback`, :ref:`spaces-collapse-pair-two-lifts`
   * - Trace restriction; trace space; half-trace
     - :class:`~orpheus.numerics.operator.TraceRestrictionOperator`, :class:`~orpheus.numerics.spaces.angular_trace_space.AngularTraceSpace`, :class:`~orpheus.numerics.spaces.angular_trace_space.AngularFaceTraceSpace`
     - :eq:`bc-trace-restriction-pair`, :ref:`bc-trace-structure`, :ref:`half-trace <bc-half-trace>`
   * - System restriction; coupled space
     - :class:`~orpheus.numerics.coupled_system.SystemRestrictionOperator`, :class:`~orpheus.numerics.coupled_system.CoupledSpace`
     - :ref:`coupled-block-system-restriction`, :eq:`coupled-block-system-restriction-pair`, :eq:`coupled-block-system-restriction-laws`; the coupled space: :ref:`coupled-block-n-general-machinery`
   * - Lift; full-field space
     - :class:`~orpheus.transport.operators.lift.BulkLift`, :class:`~orpheus.numerics.spaces.full_field_space.FullFieldSpace`
     - :ref:`cs4c-ends-select-the-body`
   * - Boundary law; realizer
     - :class:`~orpheus.geometry.boundary.BoundaryTraceLaw`, :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`
     - :eq:`affine-bc-form`, :ref:`bc-realizer-layer`, :ref:`bc-method-realizability`
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
     - :class:`~orpheus.numerics.symmetry.SubgroupOfO3`, :mod:`orpheus.numerics.invariance`
     - :ref:`manifold-realization`, :eq:`discrete-measure-g-invariance`
   * - Orbit space; quotient; descent
     - :class:`~orpheus.numerics.manifold.Quotient`, :class:`~orpheus.numerics.basis.descent.Descent`
     - :ref:`manifold-orbit-space`, :ref:`manifold-orbit-space-stabiliser`, :ref:`manifold-quotient-map`, :ref:`manifold-reynolds-projector-section`, :ref:`manifold-descent`
   * - Chart (a coordinate system's orbit map and symmetry group); singular stratum
     - :class:`~orpheus.geometry.chart.Chart` (derived from its kept columns and its linear group), :class:`~orpheus.geometry.chart.SingularStratum`; a line's image :class:`~orpheus.geometry.chart.RadialImage`, :class:`~orpheus.geometry.chart.AxialImage`
     - :ref:`chart-and-chord-chart`, :eq:`geometry-radial-coordinate`, :ref:`chart-and-chord-isotropy`, :ref:`chart-and-chord-strata`
   * - Line (Plücker coordinates); ray as a line with a start
     - :class:`~orpheus.geometry.line.Line`
     - :ref:`chart-and-chord-lines`
   * - Concentric partition; chord; crossings; measure on lines
     - :class:`~orpheus.geometry.chord.ConcentricPartition`, :class:`~orpheus.geometry.chord.Chord`, :class:`~orpheus.geometry.chord.Crossings`; :meth:`~orpheus.geometry.chart.Chart.beam_density`
     - :ref:`chart-and-chord-chord`, :eq:`geometry-line-crossing-law`, :eq:`geometry-crossing-order`, :eq:`geometry-chord-segment-lengths`, :eq:`geometry-cylinder-axial-factor`, :ref:`chart-and-chord-location`, :eq:`geometry-measure-on-lines`, :eq:`geometry-cauchy-mean-chord`
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
   * - Boundary tag to typed law, in the closed references; and its deliberate twin across the branch line
     - :data:`~orpheus.derivations.continuous.characteristic.walls.TAG_REGISTRY` (each admitted kind, with the wall's outward sign, to the typed law it names; read by :meth:`Walls.of <orpheus.derivations.continuous.characteristic.walls.Walls.of>` through the law's two factors). **Twins, never to be merged:** the reference's registry and production's parse ``orpheus.transport.method._law_from_tag``. Production's is versatile: it consults each method's admission table, reads the axis and the sign from a multi-dimensional face label, admits the ``albedo`` kind with no re-emission shape, and drops an undeclared parameter (``BC("white", {"albedo": a})`` parses as albedo 1, `#583 <https://github.com/deOliveira-R/ORPHEUS/issues/583>`_). The reference's is closed: five kinds, ``partial`` among them, each taking exactly its declared parameters and refusing any other. A shared parse would make every production change to a tag's meaning a change to the references, which the insulation principle forbids (the user's ruling of 2026-10-06).
     - :ref:`characteristic-walls-registry`, :ref:`architecture-reference-insulation`
   * - Question; mode; point
     - :class:`~orpheus.numerics.question.Eigen`, :class:`~orpheus.numerics.question.FixedSource`, :class:`~orpheus.numerics.question.Response`; :class:`~orpheus.numerics.question.Fundamental`, :class:`~orpheus.numerics.question.Nearest`; the point a :class:`~orpheus.numerics.content.FrozenMapping`
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
  :class:`~orpheus.geometry.chord.ConcentricPartition`, which has no
  consumer yet (#405).
- **Exit-report members still listed as not yet built** (the outcome
  module's ``NotYet``):
  the carrying eigen exit's balance (#354, open) and the daggered eigen
  exit (#353, open); the list also names the linear-discontinuous residual
  (#310), whose issue is closed.
- **Restriction siblings share no base** by ruling ("no consumer treats
  any restriction generically yet"); the rank-d spatial axis is
  generator-less (a gated contract); there is no ``Cone`` class, no
  embedding operator and no affine operator, each by ruling and so not
  debt.
- **The tagged guards**: ``ELEGANCE-DEBT[guard]`` occurs four times
  under ``orpheus/``: the full-field carrier (#457); the lock on withdrawn
  reference generators in ``orpheus/derivations/common/withdrawal.py``
  (#506, retired when each reference family's generator returns a
  :class:`~orpheus.reference.certificate.ReferenceCertificate` whose
  standing carries the
  :class:`~orpheus.reference.withdrawal.Withdrawal`;
  :ref:`vv-withdrawn-generators`); the Sood registry's refusal of a case
  with no citation (#405, retired when a published solution replaces the
  registry's case class); and the boundary composition's check on its
  direct children (#551, retired when the composition's operands are typed
  as responses). ``# TODO`` occurs once; ``raise NotImplementedError`` 72
  times in 27 files under ``orpheus/`` excluding ``derivations/`` (119 in
  39 with it), the population a retirement audit walks (``[M]``
  2026-10-03, a line count over ``orpheus/**/*.py``).
