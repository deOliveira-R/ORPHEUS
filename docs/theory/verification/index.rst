.. _theory-verification:

============
Verification
============

How ORPHEUS establishes that its solvers are *correct*. This part
collects the machinery and the evidence in one place: the
verification principles and evidence taxonomy, the test-harness
tagging contract, the cross-method regression protocol, the binding
reference-solution contract, the per-method evidence chapters (one
per solver family), the cross-suite evidence summary, and the
auto-generated per-equation V&V matrix.

ORPHEUS enforces a four-level verification & validation ladder —
**L0** term verification (hand calculation vs code), **L1** equation
verification (analytical / manufactured solutions with asserted
convergence order), **L2** integration (multi-group, heterogeneous),
and **L3** validation (experimental comparison) — plus an orthogonal
**foundation** bucket for software invariants outside the physics
ladder. Every test declares its rung via ``pytest.mark.l0`` /
``l1`` / ``l2`` / ``l3`` / ``foundation`` and, for equation-level
tests, the equation labels it exercises via
``pytest.mark.verifies``. Solver-vs-solver agreement is **L4
(informational)** — recorded by the cross-method protocol, never a
stand-in for L0–L3 correctness evidence. The normative ladder
definition is :ref:`vv-level-ladder` in :doc:`principles`; see
:ref:`vv-vocabulary` for the binding vocabulary.

The **error catalog** (``ERR-NNN`` entries in
``.claude/skills/vv-principles/error_catalog.md``) records every
L0-caught bug; each entry is pinned by a regression test carrying
``@pytest.mark.catches("ERR-NNN")``, and the matrix flags any
catalog entry without a catcher as a publication-blocker.

The chapters
============

.. list-table::
   :header-rows: 1
   :widths: 24 46 30

   * - Chapter
     - What it documents
     - Read it when
   * - :doc:`principles`
     - The evidence doctrine: verification vs validation,
       structural independence, the claim-layer taxonomy, the
       normative L0–L3 ladder (+ foundation + the L4 ruling), the
       three reference pillars, how the classification systems
       relate, the retired-Richardson record, the failure-mode
       catalogue, and the error-catalog contract.
     - Designing a verification strategy; judging what evidence a
       claim needs and what a reference can prove.
   * - :doc:`harness`
     - The test-harness contract: the ladder definition table,
       marker conventions, the tagging-precedence chain, the audit
       CLI, the ``.. vv-status:`` sentinel schema, and foundation
       tests.
     - Authoring or tagging any test; running or extending the
       audit.
   * - :doc:`cross_method`
     - The L4 cross-method regression protocol: solver adapters,
       agreement tolerances, reference-contamination discipline,
       and truth traceability.
     - Adding a solver or case to the cross-method net; citing an
       agreement result.
   * - :doc:`reference_solutions`
     - The binding reference contract: V&V vocabulary discipline,
       the operator-form taxonomy, the
       ``ContinuousReferenceSolution`` API, kernel primitives, and
       the reference-tier audit.
     - Writing or hardening a reference solution.
   * - :doc:`homogeneous`
     - The infinite-medium matrix eigenvalue — the shared analytical
       L1 anchor every solver family verifies against.
     - Tracing any solver's homogeneous case to its reference.
   * - :doc:`sn`
     - The discrete-ordinates evidence: the MMS ladder (1D/2D,
       heterogeneous, anisotropic, curvilinear), the eigenvalue
       cases, and the numerical-sensitivity record.
     - Working on any SN solver component.
   * - :doc:`collision_probability`
     - The CP evidence: E₃/Ki₄ semi-analytical eigenvalue grids, the
       9-case grid design, property gates, and the extended QA suite.
     - Working on the CP solver or its kernels.
   * - :doc:`method_of_characteristics`
     - The MOC evidence: the four-level suite, eigenvalue and
       convergence cases, and the per-segment MMS gate.
     - Working on the MOC solver.
   * - :doc:`monte_carlo`
     - The MC evidence: analytical and CP-referenced cases, the
       four-level suite, and the failure-mode coverage map.
     - Working on the MC solver.
   * - :doc:`diffusion`
     - The diffusion evidence: buckling anchor, continuous-reference
       cases with measured orders, and the retired-Richardson record.
     - Working on the diffusion solver.
   * - :doc:`summary`
     - The cross-suite evidence compilation: the per-method evidence
       map, the registry-generated reference-case table, cross-method
       coverage, the structural property-test inventory, convergence
       studies, and the run-book.
     - Surveying the whole suite's results; running the tests.
   * - :doc:`matrix`
     - **Auto-generated** on every Sphinx build from the test
       registry: level × module grid, per-equation coverage, orphan
       equations, documented-only labels, phantom verifies-targets,
       ERR-catalog coverage, and unmarked tests.
     - Auditing coverage; asking "what pins this equation?"

The continuous reference solvers that supply the analytical and
semi-analytical truth values are catalogued in the
:ref:`references part <theory-reference-solvers>`.

.. toctree::
   :maxdepth: 2

   principles
   harness
   cross_method
   reference_solutions
   homogeneous
   sn
   collision_probability
   method_of_characteristics
   monte_carlo
   diffusion
   summary
   matrix
