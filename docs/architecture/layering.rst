.. _architecture-layering:

The Layer Contract
==================

This page records the **layering criterion** that organizes the ORPHEUS
package tree, the package-to-layer assignment that follows from it, the
import-linter test that enforces it (:file:`tests/test_layer_imports.py`),
and the transitional exemptions captured in its ``WHITELIST``.

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
   is therefore :mod:`orpheus.transport`, not :mod:`orpheus.sn`. The SN
   B1''-eq-map factory ``AngularFlux.from_flat_with_traces`` is defined
   using SN-specific face-coordinate decoding; its home is therefore
   :mod:`orpheus.sn`.

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
     - :mod:`orpheus.sn`, :mod:`orpheus.pn` (planned),
       :mod:`orpheus.moc`, :mod:`orpheus.cp`, :mod:`orpheus.mc`,
       :mod:`orpheus.diffusion`, :mod:`orpheus.kinetics` (transitional —
       dissolves under P3.6), :mod:`orpheus.fuel`,
       :mod:`orpheus.thermal_hydraulics`, :mod:`orpheus.homogeneous`
   * - **L2** transport vocabulary
     - the transport equation's objects; method-agnostic
     - :mod:`orpheus.transport` (created by P3.3)
   * - **(input)** geometry + data
     - mesh geometry; nuclear data
     - :mod:`orpheus.geometry`, :mod:`orpheus.data`
   * - **L1** mathematics
     - functional analysis, linear algebra, measure theory; no neutrons
     - :mod:`orpheus.numerics`
   * - **L0** references
     - Branch-1 analytical / SymPy / mpmath references
     - :mod:`orpheus.derivations`

A few notes the table is too compact to capture:

* The **input layer** is not a strict member of the L0/L1/L2/L3 stack;
  it provides primitive types that every layer (including L1) may
  consume. Geometry meshes and nuclear-data structures are inputs in the
  same sense that a function argument is an input — they cross every
  layer boundary but carry no algorithmic knowledge.

* **L0** sits BELOW **L1** in the import hierarchy. The derivations live
  in :mod:`orpheus.derivations` and ship reference solvers built from
  SymPy, ``mpmath``, or pure analytical closed forms. They have *less*
  algorithmic knowledge than the production primitives in
  :mod:`orpheus.numerics` (a SymPy expression is structurally simpler
  than a numpy iteration). Production code that needs a structurally
  independent reference imports L0 — the L3-uses-L0 pattern is
  documented in :doc:`/verification/index`.

* **L4** is permissible to import everything. It is the only layer
  where wiring a run can pull in transport types, method-specific
  problems, and the math layer simultaneously. The single-file
  ``plotting.py`` is an L4 example; entry-point scripts in
  ``examples/`` are L4.


Problem and Solver are not a layer
----------------------------------

The :class:`Problem` and :class:`Solver` families are NOT layers. They
are math-object families (like :class:`Field` and
:class:`~orpheus.numerics.operator.LinearOperator`) that recur at every
layer with a layer-appropriate vocabulary:

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Layer
     - Problem (declarative)
     - Solver (iterative)
   * - L1 (math)
     - :class:`Eigenproblem` (generic, ``Ax = λx``)
     - :class:`PowerIteration`, :class:`Arnoldi`, :class:`TimeStepper`
   * - L2 (transport)
     - :class:`CriticalityProblem`, :class:`AlphaEigenproblem`,
       :class:`FixedSourceProblem`, :class:`InitialValueProblem`
     - (transport-vocabulary scheduler; method-agnostic)
   * - L3 (method)
     - (method-specific problem types if any)
     - :class:`SweepPreconditionedSolver`, DSA, TSA, JFNK

A consumer at L3 typically constructs an L2 :class:`Problem`
(``CriticalityProblem(loss, fission)``) and an L1 solver
(``PowerIteration``) and composes them. The Problem is the declarative
description; the Solver is the algorithmic iteration. They are
orthogonal axes, not layers, and they recur at each layer with
appropriate vocabulary.


The import-linter test
----------------------

The criterion is enforced by :file:`tests/test_layer_imports.py`. The
test walks every Python module under :file:`orpheus/`, parses its
imports via Python's ``ast`` module (NOT regex — regex misses
``TYPE_CHECKING`` blocks, multi-line imports, and function-body lazy
imports), and reports every edge that violates the layer contract.

The forbidden-edge dictionary is:

.. code-block:: python

   FORBIDDEN_EDGES: dict[str, frozenset[str]] = {
       # L1 imports nothing above itself.
       "numerics": L2_PACKAGES | L3_PACKAGES,

       # Input layers import L1 only.
       "geometry": L2_PACKAGES | L3_PACKAGES,
       "data":     L2_PACKAGES | L3_PACKAGES,

       # L2 imports L1 + inputs only.
       "transport": L3_PACKAGES,

       # L3 methods cannot import sibling L3 packages.
       "sn":         L3_PACKAGES - {"sn"},
       "pn":         L3_PACKAGES - {"pn"},
       # ... etc for every L3 package ...

       # L0 (derivations) sits below L1.
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

  At the time of P3.1's landing, the whitelist contains three
  ``derivations/`` → ``L3`` edges:

  .. code-block:: python

     WHITELIST: frozenset[tuple[str, str]] = frozenset({
         # RETIRE_IN_P3_FOLLOWUP — inline-import benchmark cross-check
         ("derivations/continuous/cases/diffusion.py", "diffusion"),
         # RETIRE_IN_P3_FOLLOWUP — MMS source uses MOCMesh / MOCQuadrature
         ("derivations/continuous/mms/moc.py", "moc"),
         # RETIRE_IN_P3_FOLLOWUP — sood_registry lazy-imports CPParams
         ("derivations/continuous/sood_registry/builders.py", "cp"),
     })

  Each entry is a Branch-1-uses-production-as-a-black-box benchmark —
  the reference imports a production solver to *cross-check* a
  reference value, NOT to share algebra. These are categorically
  different from algebra-sharing imports (which would be structurally
  contaminating per :doc:`/verification/index`). The retirement
  trigger for each is the module's migration to a method-side test
  or to an external benchmark harness.


When to break the rule
----------------------

The criterion is a constraint, not a moral imperative. If a real
engineering need requires a transgression, the procedure is:

1. **Make the import explicit and local.** Use a function-body lazy
   import rather than a module-level import; this localizes the
   coupling to a single function rather than to the whole module.

2. **Add a WHITELIST entry** in :file:`tests/test_layer_imports.py`
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
only made the contract executable. The 3 ``derivations/`` whitelist
entries above were the entirety of the violations across 243 Python
modules.

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
