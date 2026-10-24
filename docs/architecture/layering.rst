.. _architecture-layering:

The Layer Contract
==================

This page records the **layering criterion** that organizes the ORPHEUS
package tree, the package-to-layer assignment that follows from it, the
import-linter test that enforces it (:file:`tests/gates/test_layer_imports.py`),
the transitional exemptions captured in its ``WHITELIST``, and the second
contract the same test enforces inside L0: the closed references import no
production machinery, neither the mesh, the transport layer nor a method
package, and of :mod:`orpheus.numerics` only its interface vocabulary,
never its mathematics (:ref:`architecture-reference-insulation`).

The contract is load-bearing. The whole point of organizing code by
mathematical knowledge layer is to make a class of bugs — bugs of
*coupling* — impossible by construction. If :mod:`orpheus.numerics` could
import :mod:`orpheus.sn`, then a numerics primitive could come to depend
on an SN-specific convention; if :mod:`orpheus.transport` could import a
method package, a method's idiosyncrasy could leak into the
transport-vocabulary layer; if a method could import a sibling method,
the two would silently couple through shared types. The criterion below
forbids each of those edges; the linter makes the forbidding executable.


The criterion
-------------

   **A module's home is the lowest-knowledge layer whose vocabulary
   suffices to define it. Imports flow only from more-knowledge to
   less-knowledge.**

Two clauses, both load-bearing:

1. *Lowest-knowledge layer whose vocabulary suffices.* The math-layer
   primitive :class:`~orpheus.numerics.operator.LinearOperator` is defined
   without any neutron-physics vocabulary; its home is therefore
   :mod:`orpheus.numerics`, not :mod:`orpheus.transport` or
   :mod:`orpheus.sn`. The transport-vocabulary type
   :class:`AngularFlux` is defined using "ordinate" and "moment" — concepts
   from transport theory but not specific to any discretization; its home
   is therefore :mod:`orpheus.transport`, not :mod:`orpheus.sn`. The
   SN boundary realizer ``SNBoundaryRealizer``
   is defined using SN-specific face-coordinate decoding (the flat
   ``from_flat_with_traces`` codec that once illustrated this point has
   since moved up to the transport-layer :class:`~orpheus.transport.timed_full_field.TimedFullField`);
   its home is therefore :mod:`orpheus.sn`.

2. *Imports flow only from more-knowledge to less-knowledge.* An L3
   method package may import an L1 primitive (the method *uses* the math);
   an L1 primitive must not import from an L3 package (the math does not
   *know about* any method). The arrows point downward in the layer
   diagram below.

A useful test for whether a candidate module is in the right layer: *if
this module's docstring uses vocabulary from layer N, can it be lifted to
layer N-1 by simply removing the layer-N words?* If yes, the module
belongs in N-1 (you accidentally specialized something general). If no,
the module belongs in N (the layer-N vocabulary is load-bearing).


The layer table
---------------

The layers, top-to-bottom in the import order (each layer may import only
the layers below it):

.. list-table::
   :header-rows: 1
   :widths: 12 38 50

   * - Layer
     - Knows
     - Packages
   * - **L4** orchestration
     - wiring a run; driver / entry point
     - thin scripts; ``plotting.py``
   * - **L3** discretization
     - one method's machinery
     - :mod:`orpheus.sn`, ``orpheus.pn`` (planned — no such package
       exists yet),
       :mod:`orpheus.moc`, :mod:`orpheus.cp`, :mod:`orpheus.mc`,
       :mod:`orpheus.diffusion`, :mod:`orpheus.kinetics` (transitional —
       dissolves under P3.6), :mod:`orpheus.fuel`,
       :mod:`orpheus.thermal_hydraulics`, :mod:`orpheus.homogeneous`
   * - **L2** transport vocabulary
     - the transport equation's objects; method-agnostic
     - :mod:`orpheus.transport` (created by P3.3)
   * - **(input)** mesh
     - the discretisation overlay on the geometry: cells, subdivision,
       per-axis primitives
     - :mod:`orpheus.mesh`
   * - **(input)** specification
     - a question written with the materials and geometry it is asked of:
       the reference specification, the key a reference cache stores under
     - :mod:`orpheus.specification`
   * - **(input)** geometry + data
     - shapes, coordinate systems and boundary laws; nuclear data
     - :mod:`orpheus.geometry`, :mod:`orpheus.data`
   * - **L1** mathematics
     - functional analysis, linear algebra, measure theory; no neutrons
     - :mod:`orpheus.numerics`
   * - **(values)**
     - what a value is (the admitted reals) and when two values are the
       same (content identity); imports nothing from ``orpheus``
     - :mod:`orpheus.values`
   * - **L0** references
     - Branch-1 analytical / SymPy / mpmath references; the one exception
       to the top-to-bottom order: it may import L1 and the input layer
       (see the notes below)
     - :mod:`orpheus.derivations`

A few notes the table is too compact to capture:

* The **input layer** is not a strict member of the L0/L1/L2/L3 stack;
  it provides primitive types that the layers above it consume.
  Geometries, meshes and nuclear-data structures are inputs in the same
  sense that a function argument is an input — they cross every layer
  boundary but carry no algorithmic knowledge. The linter forbids
  :mod:`orpheus.numerics` an edge to :mod:`orpheus.data`, guarded or not
  (layer_dag.md P2, 2026-10: the energy axis is minted from its edges by
  :meth:`~orpheus.numerics.axis.EnergyAxis.of_edges`, and the views that
  read nuclear data as an axis, :meth:`EnergyGrid.as_axis
  <orpheus.data.energy_grid.EnergyGrid.as_axis>` and
  :func:`~orpheus.data.materials.energy_axis_of`, live in data), and one
  to :mod:`orpheus.mesh`. It admits the edge to :mod:`orpheus.geometry`,
  and that edge is the intended direction (layer_dag.md ruling 1): the
  point groups (:mod:`orpheus.geometry.symmetry`), the point sets they act
  on (:mod:`orpheus.geometry.manifold`) and the face naming
  (:mod:`orpheus.geometry.faces`) are geometry, and numerics' measures,
  quadratures and bases import them. Until layer DAG WP-B (2026-10-22) the
  first two were numerics modules importing
  :mod:`orpheus.geometry.transformation`. The reverse edges, geometry's
  boundary laws importing numerics, are what keeps
  {geometry, numerics} a package cycle, pending layer_dag.md P1's boundary
  carve.

* The input layer has one internal order: :mod:`orpheus.mesh` sits
  **above** :mod:`orpheus.geometry`. A mesh is an overlay on a geometry
  (it divides the geometry's intervals into cells), so
  :mod:`orpheus.mesh` may import :mod:`orpheus.geometry`,
  :mod:`orpheus.numerics` and :mod:`orpheus.data`, and never
  :mod:`orpheus.transport` or a method package; :mod:`orpheus.geometry`,
  :mod:`orpheus.data` and :mod:`orpheus.numerics` never import
  :mod:`orpheus.mesh`. Binding cross sections to the cells is not an
  input-layer job: the material mesh
  :class:`~orpheus.transport.mesh.material_mesh.MaterialMesh` stays at
  L2 in :mod:`orpheus.transport.mesh`. The package and its reason are
  on :doc:`/api/mesh`.

* The input layer has a second package above :mod:`orpheus.geometry`:
  :mod:`orpheus.specification` composes materials (:mod:`orpheus.data`), a
  geometry (:mod:`orpheus.geometry`) and a question
  (:mod:`orpheus.numerics.question`), and is neither data nor geometry, so
  it has a package of its own (#405 P1 step 8). It imports
  :mod:`orpheus.data`, :mod:`orpheus.geometry` and :mod:`orpheus.numerics`
  and never :mod:`orpheus.mesh`, :mod:`orpheus.transport`, a method
  package or :mod:`orpheus.derivations`; :mod:`orpheus.data`,
  :mod:`orpheus.geometry`, :mod:`orpheus.numerics` and :mod:`orpheus.mesh`
  never import it. The mesh and the specification are therefore siblings:
  neither imports the other. :mod:`orpheus.derivations` may import it,
  because the reference registry will hold specifications. Its two
  coordinates live one level down, each in the package whose vocabulary
  defines it: :class:`~orpheus.data.cells.CellCoefficient` names a
  ``Mixture``'s channels and lives in :mod:`orpheus.data`;
  :class:`~orpheus.geometry.extent.GeometryExtent` names an interval and
  lives in :mod:`orpheus.geometry`
  (:ref:`structured-geometry-specification-coordinates`).

* **L0** sits below **L2**, beside **L1** and the input layer, not below
  **L1**. The linter forbids :mod:`orpheus.derivations` exactly
  :mod:`orpheus.transport` and every L3 package, the same set it forbids
  :mod:`orpheus.mesh`; it forbids :mod:`orpheus.numerics`,
  :mod:`orpheus.geometry` and :mod:`orpheus.data` that set plus
  :mod:`orpheus.mesh`. So a reference may use the mathematics layer and
  describe its problem in input-layer vocabulary (a geometry, a mesh, a
  ``Mixture``, a boundary condition), and it may never name a transport object such as
  ``MaterialMesh``; lifting a reference's problem to a method's problem is
  a production verb at L2 and above. The derivations ship reference
  solvers built from SymPy, ``mpmath``, or pure analytical closed forms.
  Production code that needs a structurally independent reference
  imports L0; the L3-uses-L0 pattern is documented in
  :doc:`/theory/verification/index`. Inside L0, the closed references
  (``derivations/continuous/`` and ``derivations/common/``) are held to a
  stricter rule than the layer table: they never import
  :mod:`orpheus.mesh`, and of :mod:`orpheus.numerics` they may import only
  five interface modules (:ref:`architecture-reference-insulation`).

* **L4** is permissible to import everything. It is the only layer
  where wiring a run can pull in transport types, method-specific
  problems, and the math layer simultaneously. The single-file
  ``plotting.py`` is an L4 example; entry-point scripts in
  ``examples/`` are L4.


.. _architecture-problem-and-solver:

Problem and Solver are not a layer
----------------------------------

.. important::

   **The type names in the design table below are the design's
   vocabulary, not the tree's.** The Problem side has since been
   reified under other names: the posed question is
   :class:`~orpheus.numerics.posing.EigenPosing` or
   :class:`~orpheus.numerics.posing.SourcePosing` over an
   :class:`~orpheus.numerics.pencil.OperatorPencil`, and the Problem is
   the per-method hub (:class:`~orpheus.sn.problem.SNProblem`,
   :class:`~orpheus.homogeneous.solver.HomogeneousProblem`); the map of
   the concepts is :ref:`architecture-conceptual-view`. ``Eigenproblem``,
   ``Arnoldi``, ``TimeStepper``, ``CriticalityProblem``,
   ``AlphaEigenproblem``, ``FixedSourceProblem``, ``InitialValueProblem``
   and ``SweepPreconditionedSolver`` remain design names with no class of
   that name, written as literals so that no role asserts a class the
   interpreter cannot produce. What ORPHEUS ships in each role is
   tabulated in :ref:`architecture-problem-solver-today` immediately
   below.

The ``Problem`` and ``Solver`` families are NOT layers. They
are math-object families (like :class:`~orpheus.numerics.field.Field` and
:class:`~orpheus.numerics.operator.LinearOperator`) that recur at every
layer with a layer-appropriate vocabulary:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Layer
     - Problem (declarative)
     - Solver (iterative)
   * - L1 (math)
     - ``Eigenproblem`` (generic, ``Ax = λx``)
     - ``PowerIteration``, ``Arnoldi``, ``TimeStepper``
   * - L2 (transport)
     - ``CriticalityProblem``, ``AlphaEigenproblem``,
       ``FixedSourceProblem``, ``InitialValueProblem``
     - (transport-vocabulary scheduler; method-agnostic)
   * - L3 (method)
     - (method-specific problem types if any)
     - ``SweepPreconditionedSolver``, DSA, TSA, JFNK

A consumer at L3 would construct an L2 ``Problem``
(``CriticalityProblem(loss, fission)``) and an L1 solver
(``PowerIteration``) and compose them. The Problem is the declarative
description; the Solver is the algorithmic iteration. They are
orthogonal axes, not layers, and they recur at each layer with
appropriate vocabulary.


.. _architecture-problem-solver-today:

What fills each role today
~~~~~~~~~~~~~~~~~~~~~~~~~~

The declarative/iterative split above is settled and load-bearing, and
the Problem side is reified since 2026-09: the pencil and the two posings
are types in :mod:`orpheus.numerics`, and each method's hub mints them.
Every row below is a live cross-reference, so the table doubles as the
gap measure — a row with no live role is a genuine hole.

.. list-table::
   :header-rows: 1
   :widths: 32 68

   * - Design name (above)
     - What ORPHEUS ships
   * - L1 ``Eigenproblem`` + ``PowerIteration``
     - The method-agnostic
       :class:`~orpheus.numerics.eigenvalue.EigenvalueSolver` Protocol is
       the boundary, and
       :func:`~orpheus.numerics.eigenvalue.power_iteration` is the single
       power-iteration loop in the codebase. The *problem* is a type
       since 2026-09: :class:`~orpheus.numerics.posing.EigenPosing`, a
       pencil with a spectral map, minted by the hub; the loop consumes
       the pair of methods a solver exposes across that Protocol.
   * - L2 ``CriticalityProblem``
     - :class:`~orpheus.numerics.iteration.KEigenvalue` — the
       operator-triple realization of that same boundary, carrying the
       k-posing :math:`A_{\rm loss} = A - S`, :math:`M = F`,
       :math:`k = \mu` (the full posing table is at
       :ref:`eigenvalue-posing`). The declarative type is
       :class:`~orpheus.numerics.posing.EigenPosing` with
       :data:`~orpheus.numerics.posing.K_MAP`, minted by
       :attr:`~orpheus.sn.problem.SNProblem.eigen_posing` and by the
       homogeneous hub.
   * - L2 ``FixedSourceProblem``
     - :class:`~orpheus.numerics.iteration.SourceIteration`, with
       :class:`~orpheus.numerics.iteration.KrylovAcceleration` as the
       accelerated arm; the declarative type is
       :class:`~orpheus.numerics.posing.SourcePosing`, minted by
       :meth:`~orpheus.sn.problem.SNProblem.source_posing` as the pencil's
       member at :math:`\sigma = 1`.
   * - L2 ``AlphaEigenproblem``
     - Not built. The :math:`\alpha`-eigenvalue row
       (:math:`A_{\rm loss} = L+C-S-N_{2n}-F-B`, :math:`M = 1/v`,
       :math:`\alpha = -1/\mu` with :math:`\mu` the eigenvalue of
       :math:`A^{-1}M`) is a documented seam in
       :mod:`orpheus.numerics.eigenvalue`'s package header — a posing the
       existing loop would accept, with no constructor yet.
   * - L2 ``InitialValueProblem``
     - Not built at this boundary. The coupled point-kinetics /
       thermal-hydraulics transient in :mod:`orpheus.kinetics` integrates
       its own ODE state vector through ``scipy.integrate.solve_ivp``; it
       never poses a transport operator, so it is not an instance of this
       row.
   * - L1 ``Arnoldi`` / ``TimeStepper``
     - Reserved at the ``eigenvalue_method`` constructor selector on
       :class:`~orpheus.numerics.iteration.KEigenvalue`: only ``"power"``
       is implemented, and any other value raises at construction rather
       than failing later.
   * - L3 ``SweepPreconditionedSolver``
     - Diffusion-synthetic acceleration ships as
       :class:`~orpheus.sn.acceleration.dsa.DSACorrection` over
       :class:`~orpheus.sn.acceleration.dsa.DSALowOrderSystem` (see
       :ref:`sn-acceleration`); TSA and JFNK are not built.

Read the two tables together: the *vocabulary* of this section is
settled, and its reification landed as the pencil, the posings and the
per-method hubs rather than as an ``orpheus.transport.problems``
package; what remains open (a shared Problem type over the three hubs,
the diffusion hub's pencil, the α posing) is the debt list of
:ref:`architecture-conceptual-view`.


The import-linter test
----------------------

The criterion is enforced by :file:`tests/gates/test_layer_imports.py`. The
test walks every Python module under :file:`orpheus/`, parses its
imports via Python's ``ast`` module (NOT regex — regex misses
``TYPE_CHECKING`` blocks, multi-line imports, and function-body lazy
imports), and reports every edge that violates the layer contract.

The forbidden-edge dictionary is:

.. code-block:: python

   FORBIDDEN_EDGES: dict[str, frozenset[str]] = {
       # The values package is the bottom of the graph: it imports
       # nothing from orpheus (layer_dag.md ruling 6).
       "values": ORPHEUS_TOP_LEVEL - VALUES_PACKAGES,

       # L1 imports nothing above itself; data sits strictly above it
       # (layer_dag.md P2). Geometry is still admitted, pending P1.
       "numerics": DATA_PACKAGES | MESH_PACKAGES | SPECIFICATION_PACKAGES
                   | REFERENCE_PACKAGES | L2_PACKAGES | L3_PACKAGES,

       # Geometry and data never import the mesh overlay, the
       # specification or the reference package above them.
       "geometry": MESH_PACKAGES | SPECIFICATION_PACKAGES | REFERENCE_PACKAGES | L2_PACKAGES | L3_PACKAGES,
       "data":     MESH_PACKAGES | SPECIFICATION_PACKAGES | REFERENCE_PACKAGES | L2_PACKAGES | L3_PACKAGES,

       # The mesh imports geometry, data and L1, never the specification,
       # the reference package, L2 or L3.
       "mesh": SPECIFICATION_PACKAGES | REFERENCE_PACKAGES | L2_PACKAGES | L3_PACKAGES,

       # The specification imports data, geometry and L1 only.
       "specification": MESH_PACKAGES | REFERENCE_PACKAGES | L2_PACKAGES | L3_PACKAGES | L0_PACKAGES,

       # The reference package sits above the specification, below
       # transport, L3 and derivations.
       "reference": MESH_PACKAGES | L2_PACKAGES | L3_PACKAGES | L0_PACKAGES,

       # L2 imports L1 + inputs only.
       "transport": L3_PACKAGES,

       # L3 methods cannot import sibling L3 packages.
       "sn":         L3_PACKAGES - {"sn"},
       "pn":         L3_PACKAGES - {"pn"},
       # ... etc for every L3 package ...

       # L0 (derivations) imports L1 + inputs only, as the input layers do.
       "derivations": L2_PACKAGES | L3_PACKAGES,
   }

A failing test names the offending module in the parametrised test ID,
so the bug-finding signal is module-local. The test is tagged
``@pytest.mark.foundation`` — a software contract, not a theory
claim.


Tolerances
----------

The linter ships with one tolerance, the WHITELIST below. An import inside
an ``if TYPE_CHECKING:`` block is an edge like any other.

**The retired TYPE_CHECKING exemption (layer-DAG P3).** Until P3 the
linter skipped a ``TYPE_CHECKING``-guarded import when the source layer
was L1 or L2 (layer-DAG P2 first narrowed it to an L2 source and an L3
target, when it removed the two guarded ``numerics`` → ``data`` imports), on the argument that a typing-only import creates no runtime
edge. The census at retirement (`[M]` 2026-10-22, over the 393 modules of
``orpheus/``, with the gate's own import walker) found 7 tolerated imports
in 6 files: all 7 were ``transport`` naming ``sn`` or ``diffusion`` types,
3 of them were dead (imported and never used), and 0 came from
``numerics``. The other 4 typed a transport factory's carrier as
``SNProblem`` where it read only transport vocabulary (a ``MaterialMesh``
member, a surface a transport Protocol names, the ``angular_trace`` door),
or sat on an SN algorithm living in L2 (the per-ordinate q½ source fold,
moved to :mod:`orpheus.sn.operators.radial_characteristic`). The runtime
argument was true and beside the point: a layer is a set of concepts, and
an L2 annotation naming a method type says that the L2 concept depends on
the method's. The exemption hid exactly those misplaced concepts, so it is
retired, and a leg in ``_FORBIDDEN_LEGS`` witnesses that a
``TYPE_CHECKING`` import from ``transport`` into ``sn`` (and from
``numerics`` into ``transport``) is refused.

**WHITELIST.**

  An explicit ``frozenset[tuple[str, str]]`` of
  ``(module_relative_path, target_top_level_package)`` pairs that the
  linter MUST pass even though :data:`FORBIDDEN_EDGES` would reject
  them. Every entry carries a ``RETIRE_IN_P3_FOLLOWUP`` comment naming
  its retirement trigger.

  The whitelist is empty since layer DAG P4b and P4c (2026-10-22):

  .. code-block:: python

     WHITELIST: frozenset[tuple[str, str]] = frozenset()

  Each of its last two rows was a reference module that posed a production
  problem (a manufactured-solution harness for a method): it used production
  as a black box, NOT to share algebra. That is categorically
  different from an algebra-sharing import (which would be structurally
  contaminating per :doc:`/theory/verification/index`). The retirement
  trigger was the module's migration to a method-side test
  or to an external benchmark harness.

  A third row, ``("derivations/continuous/sood_registry/builders.py",
  "cp")``, was retired by layer DAG P4 (2026-10-22) with the module it
  exempted: its lazy import of ``CPParams`` was the only edge closing the
  package cycle between ``cp`` and ``derivations``, and the function it
  served had no caller. A Sood case's interface to production is its
  ``to_geometry()`` and its ``materials``, and a consumer meshes the
  geometry itself. ``test_every_exemption_row_names_a_live_edge`` (below)
  keeps a row from outliving its edge.

  The row ``("derivations/continuous/mms/moc.py", "moc")``, the second
  of the two rows the whitelist held until then, was retired by layer DAG P4b (2026-10-22). The MoC manufactured-solution
  module built the production track mesh (``MOCMesh`` over an
  ``MOCQuadrature``) for its one test; that builder, ``build_moc_mesh``,
  moved into the test, ``tests/gates/moc/test_mms.py``, and the reference
  keeps only the case's declarations in plain data (the FSR partition
  ``fsr_edges(n_annuli)``, which owns the equal-area annuli and the
  square-border cell out to the Wigner-Seitz radius; the outer law; the
  quadrature orders and ray spacing), the closed ansatz and the
  standalone sweep, whose track-mesh argument it reads without naming
  its type. That last read is a duck type no gate sees: a rename in
  ``MOCMesh`` is caught only by the MMS test. Since the
  continuous-reference walk imports every ``orpheus.derivations`` module,
  the import had also loaded ``orpheus.moc`` into every consumer of the
  registry; ``test_loading_the_reference_registry_loads_no_moc_module``
  now reads ``sys.modules`` in a fresh interpreter to keep it out.

  The other row, ``("derivations/continuous/mms/sn.py", "transport")``, was
  retired by layer DAG P4c (2026-10-22). The non-vacuum S\ :sub:`N`
  manufactured-solution cases had wrapped their closed inflow trace into
  an :class:`~orpheus.transport.fields.AngularBoundarySourceSink` inside
  function bodies, and the module carried a bundler
  (``build_nonvacuum_fixed_source``) into a
  :class:`~orpheus.transport.timed_full_field.TimedFullField`. Each case
  now declares the trace in plain data
  (``inflow_face_values(mesh, n_face_moments)``, a ``{face: values}``
  dict), and the translation into production lives with the tests that
  pair a case with a solve (``mms_prescribed_inflow`` and
  ``build_nonvacuum_fixed_source`` in ``tests/gates/sn/_test_helpers.py``).
  Moving the bundler into ``transport`` instead was refused: transport
  would then name a manufactured-solution case protocol, a verification
  concept in production.


.. _architecture-reference-insulation:

The closed references import no production machinery
----------------------------------------------------

The layer table lets L0 import L1 whole and the mesh overlay. For the
closed references both are narrowed, by the user's rulings of 2026-10-06,
whose principle the user gave in these words:

   "the reference methods (like Fn or trajectory resolvent) are highly
   closed, purpose built and limited scope. They are fundamentally
   different than production, which is versatile, generalist and large in
   scope. So they ask for different things from their machinery and the
   churn should be concentrated on production, whereas once references
   reach a good architecture, churn should be extremely limited."

**What the principle implies.** A reference is valuable because it stays
put: a value it produced last month is comparable with the value it
produces today. Production numerics is the opposite kind of code. The
operator algebra, the pencil and the iteration family in
:mod:`orpheus.numerics` are built to be general (any operator on any space,
the adjoint derived through ``.H``) and they change whenever production
grows. If a reference computed with that machinery, every production
refactor would be a change to every reference, and agreement between a
reference and a production solver would rest partly on shared code. So a
reference depends only on slow-moving upstreams: numpy, scipy, mpmath, and
the reference kernel of its own, which lives in
``orpheus/derivations/common/`` and is small, closed and changed only by
ruling. Its dense linear algebra, the module
:mod:`orpheus.derivations.common.dense_pencil`, is derived on
:ref:`verification-reference-kernel`; its composite Gauss rule is
:func:`~orpheus.derivations.common.quadrature.composite_gauss_legendre`.

**What a closed reference may import.** Two rules, one per kind of
upstream:

* **No production machinery.** A closed reference never imports
  :mod:`orpheus.mesh`, :mod:`orpheus.transport` or any method package
  (the set ``REFERENCE_FORBIDDEN_PACKAGES`` in
  :file:`tests/gates/test_layer_imports.py`, the mesh, L2 and L3 packages
  of the layer table; the user's ruling "Add mesh and methods"). The layer
  table already forbids L0 the transport layer and the methods; the new
  member is the mesh, the discretisation overlay a reference has no use
  for, since it poses its problem on the geometry and evaluates on its own
  points. Its input vocabulary stays importable: :mod:`orpheus.data`,
  :mod:`orpheus.geometry`, :mod:`orpheus.specification` and
  :mod:`orpheus.reference`.
* **Of numerics, only the interface.** Below.

**Interface versus mathematics.** A reference still has to speak to its
consumers: it receives a question, it is read through observables, it is
keyed and memoised by content. That vocabulary is shared, so it is
imported; what is computed is not. The allowlist is the set
``REFERENCE_INTERFACE_NUMERICS`` in :file:`tests/gates/test_layer_imports.py`,
four submodules of :mod:`orpheus.numerics`. Content identity and the scalar
parsers were members until layer DAG ruling 6 (2026-10-22) moved them to
:mod:`orpheus.values`, below every package, which every reference imports
freely; a reference still normalises a memo key's argument at its definition
site with :func:`~orpheus.values.scalars.canonical_reals` or
:func:`~orpheus.values.scalars.parse_real`, the rule the content encoder
applies, so the key and the stored value agree (#405 P3 step 2c, C1).

.. list-table::
   :header-rows: 1
   :widths: 26 74

   * - Submodule
     - What a reference takes from it
   * - :mod:`orpheus.numerics.question`
     - the question asked of a system (:class:`~orpheus.numerics.question.Eigen`,
       :class:`~orpheus.numerics.question.FixedSource`,
       :class:`~orpheus.numerics.question.Response`) and the mode selector
   * - :mod:`orpheus.numerics.observable`
     - what is read off an answer (an eigenvalue, a flux integral, a point
       value, a ratio)
   * - :mod:`orpheus.numerics.mesh_free_function`
     - the mesh-free functions a specification states (a per-region table,
       a symbolic function)
   * - :mod:`orpheus.numerics.traced_memo`
     - the traced memo, a pure function of content memoised on disk

The test for membership is whether the module computes anything a
reference's value depends on. The posed-question types
:class:`~orpheus.numerics.pencil.OperatorPencil`,
:class:`~orpheus.numerics.posing.EigenPosing` and
:class:`~orpheus.numerics.posing.SourcePosing` fail it: they are built on
:class:`~orpheus.numerics.operator.LinearOperator`, the production operator
algebra, so they are mathematics, and the reference spells its own
:class:`~orpheus.derivations.common.dense_pencil.DensePencil` instead.

**The exemptions, and why.** The rule binds ``derivations/continuous/``
and ``derivations/common/`` (the user's ruling on the gate's scope,
2026-10-06: "Continuous references only"). Outside it, by design:

* ``derivations/discrete/`` is outside the two roots. Its modules are
  algebras of record of a production discretization, so their subject is
  a production object and they follow production's changes by nature:
  :mod:`orpheus.derivations.discrete.sn.balance` imports the roots of
  unity on purpose (since layer DAG WP-B, 2026-10-22, the module is
  :mod:`orpheus.geometry.roots_of_unity`; the move does not widen what a
  closed reference may import, so the insulation gate refuses
  :mod:`orpheus.geometry.roots_of_unity`, :mod:`orpheus.geometry.symmetry`
  and :mod:`orpheus.geometry.manifold` under the two roots by name).
* One module under the roots builds production objects and is exempt
  by name, its entry of ``_REFERENCE_INSULATION_EXEMPT`` carrying its
  reason in the gate. The S\ :sub:`N` manufactured-solution harness poses
  production meshes: ``derivations/continuous/mms/sn.py`` imports
  :mod:`orpheus.mesh` and the production
  :mod:`orpheus.numerics.quadrature` and
  :mod:`orpheus.numerics.moment_layout` to pose S\ :sub:`N` meshes (it
  imports nothing from :mod:`orpheus.transport` since layer DAG P4c,
  2026-10-22). ``derivations/continuous/mms/moc.py`` was the second exempt
  module until layer DAG P4b (2026-10-22), when its MoC problem posing
  moved to its test; since then it imports nothing from :mod:`orpheus.mesh`
  or :mod:`orpheus.moc`.

**The deliberate twins.** A primitive both branches need exists once on
each side of the line: the Perron–Frobenius extraction is
:meth:`DensePencil.fundamental <orpheus.derivations.common.dense_pencil.DensePencil.fundamental>`
for the references and :func:`~orpheus.numerics.eigenvalue.dominant_eigenpair`
for production. The duplicate is the point (each copy is verified on its
own, so a reference agreeing with a production solver shares no
eigen-extraction code with it), and it is recorded by name in the concept
table of :ref:`architecture-conceptual-view` so that no later review merges
the two.

**The gate.** :file:`tests/gates/test_layer_imports.py`, section "The
closed references import only the interface vocabulary of numerics",
four ``foundation`` tests, and a fifth in the section "Every exemption
row names a live edge":

* ``test_closed_references_import_only_the_numerics_interface`` walks every
  module under ``orpheus/derivations/``, keeps those under the two roots
  and not exempt, and fails on any import of a package in
  ``REFERENCE_FORBIDDEN_PACKAGES``, and on any import of an
  :mod:`orpheus.numerics` submodule outside the allowlist or of the
  package itself; the two refusals carry different messages. It reads
  imports with the same ``ast`` walker as the layer linter, plus the
  string-literal argument of ``importlib.import_module`` and
  ``__import__``. As in the layer linter since layer-DAG P3, an import
  inside ``if TYPE_CHECKING:`` is refused here: a reference has no reason
  to name a production type even for a type checker.
* ``test_reference_insulation_refuses_each_import_shape`` feeds nine
  synthetic sources, one per import shape (absolute, ``import a.b.c``,
  relative, function-body, ``TYPE_CHECKING``, the bare package by
  ``from orpheus.numerics import``, ``importlib.import_module``, and the
  package as a name, absolute and relative), and requires exactly one
  numerics violation from each.
* ``test_reference_insulation_refuses_production_machinery`` feeds six:
  an absolute import of the mesh, a function-body import from the
  transport layer, a method package, the mesh package as a name, a
  relative import of a method package, and a method import from the path
  of the retired Sood builder (re-adding its exemption turns that leg
  red), and requires exactly one production-machinery violation from
  each.
* ``test_reference_insulation_admits_the_interface_and_the_exempt``
  requires no violation from nine sources: three interface imports, an
  exempt harness importing a production quadrature, a module under
  ``derivations/discrete/``, numpy with scipy, an exempt MMS harness
  importing a method, the geometry with the data, and the specification
  with the reference package.
* ``test_every_exemption_row_names_a_live_edge`` removes each
  ``WHITELIST`` row and each ``_REFERENCE_INSULATION_EXEMPT`` entry in
  turn and requires the named file to exist and to redden the predicate
  the row exempts (2 of 2 rows and 2 of 2 entries live after P4). Its
  first red is the tree of P4 with the Sood builder deleted and its two
  rows not yet retired: it reports both as naming a file that does not
  exist.

Its first red is the tree before the reference kernel landed. `[M]`
2026-10-06: the gate's own predicate, run on
``orpheus/derivations/common/eigenvalue.py`` as committed at ``a15cb1b7``,
returns one violation, the function-body import of
:mod:`orpheus.numerics.eigenvalue` by the infinite-medium adjoint spectrum,
which that module now takes from the dense pencil; on the reference
kernel's tree it returns none. Every other non-interface import of
:mod:`orpheus.numerics` under ``orpheus/derivations/`` is exempt (`[M]` the
same day, ``git grep`` over the import lines of ``orpheus/derivations``,
4 import lines in 2 files): ``Quadrature`` and ``moment_layout`` in
``continuous/mms/sn.py``, and ``roots_of_unity`` in
``discrete/sn/balance.py``. The production-machinery rule's first red is
the Sood builder without its exemption (`[M]` 2026-10-06, the gate's
predicate run with the entry removed): two violations for
:mod:`orpheus.mesh` and one for ``orpheus.cp.solver``; on the tree with the
exemption the whole gate returns none. P4 (2026-10-22) deleted the
builder and retired both rows; the first red is the recreated module's
function-body ``CPParams`` import, refused by both tree-wide tests
(``test_no_forbidden_imports`` with "forbidden: derivations → cp", and
``test_closed_references_import_only_the_numerics_interface`` with "never
imports production machinery").

**A change to the allowlist is a contract change.** The interface
vocabulary is itself a channel through which production churn can reach a
reference, so the allowlist is kept minimal. Adding a submodule to
``REFERENCE_INTERFACE_NUMERICS``, removing a package from
``REFERENCE_FORBIDDEN_PACKAGES``, or adding a module to
``_REFERENCE_INSULATION_EXEMPT``, is reviewed as a change to this contract:
the commit says why the module computes nothing a reference's value depends
on (or, for an exemption, why its subject is a production object), and the
same reasoning is written here. The geometric kernel both branches import
(:mod:`orpheus.geometry.chart`, :mod:`orpheus.geometry.line`,
:mod:`orpheus.geometry.chord`) is held to the same standard of stability by
the same ruling; no gate enforces that, so it rests on review.


When to break the rule
----------------------

The criterion is a constraint, not a moral imperative. If a real
engineering need requires a transgression, the procedure is:

1. **Make the import explicit and local.** Use a function-body lazy
   import rather than a module-level import; this localizes the
   coupling to a single function rather than to the whole module.

2. **Add a WHITELIST entry** in :file:`tests/gates/test_layer_imports.py`
   with a ``RETIRE_IN_<phase-or-issue>`` comment naming the
   retirement trigger. The whitelist makes the exemption visible and
   gives a future contributor a place to start when refactoring.

3. **Open an issue** describing why the layering needs adjustment, OR
   why the module needs to move to a different layer, OR why the
   criterion itself needs revision.

The linter is a tool for catching unintentional coupling — not a tool
for forbidding intentional coupling. The discipline is "every
exemption is named and justified", not "no exemptions exist."


Historical context
------------------

The layer contract was formalized in Phase 3 of the
``moment-space-and-layering`` plan (2026-05). The packages had largely
converged on the contract *before* the linter landed (the discipline
had been enforced by earlier waves of refactoring); the P3.1 commit
only made the contract executable. Three ``derivations/`` whitelist
entries were the entirety of the violations across 243 Python modules:
the MoC harness's, still listed; the Sood builder's, retired by layer
DAG P4 (2026-10-22); and the diffusion cases' benchmark cross-check,
since retired; the S\ :sub:`N`
harness's entry came later.

The earlier waves that converged on the contract:

* Wave 0 → Wave 11 (the typed-field / boundary-realizer cascade) moved
  shape contracts into :mod:`orpheus.numerics` and method-specific
  realizers into :mod:`orpheus.sn`, retiring the legacy shared types
  that crossed L1/L3 boundaries.

* Phase 1 of moment-space-and-layering (2026-05) added the typed
  spherical-harmonic space at L1 and split the SN-specific
  ``apply_traced`` into a generic moment-projection primitive at L1
  + a thin SN consumer at L3. The Frame/Basis carve later re-homed
  that primitive as the spherical-harmonic
  :class:`~orpheus.numerics.frame.GalerkinFrame`'s ``analysis`` face.

The Phase 3 refactor packages the convergence as an enforced contract;
subsequent Phase 3 steps (P3.2 through P3.6) make further structural
moves under the protection of the linter.

The closed-reference contract (:ref:`architecture-reference-insulation`)
was added on 2026-10-06 with the reference kernel, as the precursor of the
characteristic references' rebuild (#405); its first red was the
infinite-medium adjoint spectrum's import of the production
eigen-extraction.
