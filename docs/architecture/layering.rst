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
   SN boundary realizer :class:`~orpheus.sn.boundary.realizer.SNBoundaryRealizer`
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
  boundary but carry no algorithmic knowledge. The linter does not
  forbid :mod:`orpheus.numerics` an edge to :mod:`orpheus.geometry` or
  :mod:`orpheus.data`; it does forbid one to :mod:`orpheus.mesh`.

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
       # L1 imports nothing above itself.
       "numerics": MESH_PACKAGES | SPECIFICATION_PACKAGES | L2_PACKAGES | L3_PACKAGES,

       # Geometry and data never import the mesh overlay or the
       # specification above them.
       "geometry": MESH_PACKAGES | SPECIFICATION_PACKAGES | L2_PACKAGES | L3_PACKAGES,
       "data":     MESH_PACKAGES | SPECIFICATION_PACKAGES | L2_PACKAGES | L3_PACKAGES,

       # The mesh imports geometry, data and L1, never the specification,
       # L2 or L3.
       "mesh": SPECIFICATION_PACKAGES | L2_PACKAGES | L3_PACKAGES,

       # The specification imports data, geometry and L1 only.
       "specification": MESH_PACKAGES | L2_PACKAGES | L3_PACKAGES | L0_PACKAGES,

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

The linter ships with two tolerances:

**TYPE_CHECKING exemption.**

  Imports inside an ``if TYPE_CHECKING:`` block do not create a runtime
  edge — they exist only for static type checkers (mypy, pyright). An L1
  or L2 module may legitimately import an L3 type *inside* a
  ``TYPE_CHECKING`` block when the type appears only in a string-quoted
  annotation. The linter's ``ast`` walker recognizes ``TYPE_CHECKING``
  guards and skips imports inside them when the source layer is L1 or
  L2.

**WHITELIST.**

  An explicit ``frozenset[tuple[str, str]]`` of
  ``(module_relative_path, target_top_level_package)`` pairs that the
  linter MUST pass even though :data:`FORBIDDEN_EDGES` would reject
  them. Every entry carries a ``RETIRE_IN_P3_FOLLOWUP`` comment naming
  its retirement trigger.

  The whitelist holds three ``derivations/`` edges into L2 and L3:

  .. code-block:: python

     WHITELIST: frozenset[tuple[str, str]] = frozenset({
         # RETIRE_IN_P3_FOLLOWUP — MMS source uses MOCMesh / MOCQuadrature
         ("derivations/continuous/mms/moc.py", "moc"),
         # RETIRE_IN_P3_FOLLOWUP — sood_registry lazy-imports CPParams
         ("derivations/continuous/sood_registry/builders.py", "cp"),
         # RETIRE_IN_P3_FOLLOWUP — the non-vacuum MMS reference lazily builds
         # its prescribed-inflow source from transport vocabulary
         ("derivations/continuous/mms/sn.py", "transport"),
     })

  Each entry is a reference module that poses or builds a production
  problem (a manufactured-solution harness for a method, the Sood
  registry's collision-probability builder): it uses production as a
  black box, NOT to share algebra. These are categorically
  different from algebra-sharing imports (which would be structurally
  contaminating per :doc:`/theory/verification/index`). The retirement
  trigger for each is the module's migration to a method-side test
  or to an external benchmark harness.


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
five submodules of :mod:`orpheus.numerics`:

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
   * - ``orpheus.numerics.content``
     - content identity, the one encoder of a persistent key
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
  :mod:`orpheus.derivations.discrete.sn.balance` imports
  :mod:`orpheus.numerics.roots_of_unity` on purpose.
* Three modules under the roots build production objects and are exempt
  by name, each entry of ``_REFERENCE_INSULATION_EXEMPT`` carrying its
  reason in the gate. The two manufactured-solution harnesses pose a
  production method: ``derivations/continuous/mms/sn.py`` imports
  :mod:`orpheus.mesh`, :mod:`orpheus.transport` and the production
  :mod:`orpheus.numerics.quadrature` and
  :mod:`orpheus.numerics.moment_layout` to pose S\ :sub:`N` problems, and
  ``derivations/continuous/mms/moc.py`` imports :mod:`orpheus.mesh` and
  :mod:`orpheus.moc`. ``derivations/continuous/sood_registry/builders.py``
  builds production collision-probability problems, importing
  :mod:`orpheus.mesh` twice and ``orpheus.cp.solver`` once. The transport
  and method imports of these three are also entries of the layer
  linter's ``WHITELIST`` above.

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
four ``foundation`` tests:

* ``test_closed_references_import_only_the_numerics_interface`` walks every
  module under ``orpheus/derivations/``, keeps those under the two roots
  and not exempt, and fails on any import of a package in
  ``REFERENCE_FORBIDDEN_PACKAGES``, and on any import of an
  :mod:`orpheus.numerics` submodule outside the allowlist or of the
  package itself; the two refusals carry different messages. It reads
  imports with the same ``ast`` walker as the layer linter, plus the
  string-literal argument of ``importlib.import_module`` and
  ``__import__``. Unlike the layer linter's L1/L2 tolerance, an import
  inside ``if TYPE_CHECKING:`` is refused here: a reference has no reason
  to name a production type even for a type checker.
* ``test_reference_insulation_refuses_each_import_shape`` feeds seven
  synthetic sources, one per import shape (absolute, ``import a.b.c``,
  relative, function-body, ``TYPE_CHECKING``, the bare package, and
  ``importlib.import_module``), and requires exactly one numerics
  violation from each.
* ``test_reference_insulation_refuses_production_machinery`` feeds four:
  an absolute import of the mesh, a function-body import from the
  transport layer, a method package, and a relative import of a method
  package, and requires exactly one production-machinery violation from
  each.
* ``test_reference_insulation_admits_the_interface_and_the_exempt``
  requires no violation from eight sources: two interface imports, an
  exempt harness, a module under ``derivations/discrete/``, numpy with
  scipy, the exempt Sood builder importing a method, the geometry with the
  data, and the specification with the reference package.

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
exemption the whole gate returns none.

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
the two still listed for the MoC harness and the Sood builder, and the
diffusion cases' benchmark cross-check, since retired; the S\ :sub:`N`
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
