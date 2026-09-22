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
:class:`~orpheus.numerics.operator.AxisSectionOperator`. Along an index
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
:class:`~orpheus.numerics.gauge.ScaleGauge`), certified
(:class:`~orpheus.numerics.outcome.ExitCertificate`), and packaged once.

The **Solution** (:ref:`the-solution-outcome`;
:ref:`sn-solution-carries-its-posing`) is the frozen five-tuple
*problem, outcome, strategy, certificate, record*
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
   * - Retraction; section
     - :class:`~orpheus.numerics.operator.AxisRetractionOperator`, :class:`~orpheus.numerics.operator.AxisSectionOperator`
     - :ref:`spaces-collapse-pair`, :ref:`spaces-collapse-pair-naming`
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
   * - Material mesh; method hubs
     - :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh`, :class:`~orpheus.sn.problem.SNProblem`, :class:`~orpheus.diffusion.augmented_mesh.DiffusionMesh`, :class:`~orpheus.homogeneous.solver.HomogeneousProblem`
     - :ref:`architecture-layering`; hub: :ref:`the-problem-hub`
   * - Pencil; posings
     - :class:`~orpheus.numerics.pencil.OperatorPencil`, :class:`~orpheus.numerics.posing.EigenPosing`, :class:`~orpheus.numerics.posing.SourcePosing`
     - :ref:`the-operator-pencil`, :eq:`pencil-family`, :ref:`eigenvalue-posing`, :ref:`sn-the-problem-poses-its-pencil`
   * - Strategy; splitting; schedule
     - :class:`~orpheus.sn.splitting.Splitting`, :func:`~orpheus.sn.splitting.resolve_schedule`
     - :ref:`sn-splitting-is-a-strategy-value`
   * - Resolvent
     - :class:`~orpheus.sn.operators.sweep_operator.SweepOperator`, :class:`~orpheus.numerics.coupled_system.CoupledSubstitutionOperator`
     - :eq:`eigen-resolvent`
   * - Iteration; outer loop
     - :class:`~orpheus.numerics.iteration.SourceIteration`, :class:`~orpheus.numerics.iteration.KrylovAcceleration`, :func:`~orpheus.numerics.eigenvalue.power_iteration`
     - :ref:`eigenvalue-posing`
   * - Record; certificate; gauge
     - :class:`~orpheus.numerics.convergence.IterationRecord`, :class:`~orpheus.numerics.outcome.ExitCertificate`, :class:`~orpheus.numerics.gauge.ScaleGauge`
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
  convergence certificate reads the Splitting's residual, so a
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
- **Certificate members still listed as not yet built** (the outcome
  module's ``NotYet``):
  the carrying eigen exit's balance (#354, open) and the daggered eigen
  exit (#353, open); the list also names the linear-discontinuous residual
  (#310), whose issue is closed.
- **Restriction siblings share no base** by ruling ("no consumer treats
  any restriction generically yet"); the rank-d spatial axis is
  generator-less (a gated contract); there is no ``Cone`` class, no
  embedding operator and no affine operator, each by ruling and so not
  debt.
- **The one tagged guard**: ``ELEGANCE-DEBT[guard]`` occurs once under
  ``orpheus/`` (the full-field carrier, #457); ``# TODO`` once;
  ``raise NotImplementedError`` 66 times in 24 files under ``orpheus/``
  excluding ``derivations/`` (112 in 35 with it), the population a
  retirement audit walks.
