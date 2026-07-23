.. _theory-verification:

============
Verification
============

How ORPHEUS establishes that its solvers are *correct*. This part
collects the machinery and the evidence in one place: the
verification principles and evidence taxonomy, the test-harness
tagging contract, the cross-method regression protocol, the binding
reference-solution contract, and the auto-generated per-equation
V&V matrix.

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
stand-in for L0–L3 correctness evidence. See :ref:`vv-vocabulary`
for the binding definitions.

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
     - The evidence taxonomy: the three reference pillars and their
       boundaries, error-catalog vs paper-floor evidence classes,
       and the three-meanings Green's-function taxonomy the suite
       consumes.
     - Designing a verification strategy; judging what a reference
       can and cannot prove.
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
   matrix
