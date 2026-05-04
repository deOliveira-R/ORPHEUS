.. _verification-index:

==========================================
Verification — V&V infrastructure & matrix
==========================================

.. contents:: Contents
   :local:
   :depth: 2


Key Facts
=========

**Read this page first when investigating a verification gap, designing
a new test, or auditing the V&V status of a solver.**

- ORPHEUS verifies every solver against analytical or semi-analytical
  references derived from the solver's *own* governing equation. No
  cross-solver "benchmarking" stands in for verification — that is
  level L4 (informational), strictly distinct from L0–L3 (correctness).
  See :ref:`vv-vocabulary` for the binding definitions.
- The verification suite carries three artefacts on every Sphinx build:

  * :doc:`reference_solutions` — vocabulary discipline,
    operator-form taxonomy, kernel primitives, and the
    :class:`~orpheus.derivations.ContinuousReferenceSolution`
    contract.
  * :doc:`matrix` — auto-generated V&V matrix from the test registry
    (level × module × equation grid), equation coverage table,
    L0 error-catalog coverage, and the orphan-equation list.
  * :doc:`/theory/reference_solvers` — the Pillar-2 reference-solver
    catalogue with the three-meanings Green's-function taxonomy and
    per-solver pillar classification.

- The **error catalog** (``ERR-NNN`` entries in
  ``.claude/skills/vv-principles/error_catalog.md``) records every
  L0-caught bug. Each entry names the failure mode (1–6 from the AI
  failure-modes table), the evidence-class that previously hid it,
  and the test that catches it via ``@pytest.mark.catches("ERR-NNN")``.
  The matrix's "L0 error-catalog coverage" section flags any catalog
  entry without a catcher as a **publication-blocker**.


.. _verification-evidence-classes:

Evidence classes — the three pillars and their boundaries
=========================================================

Every verification reference belongs to exactly one of three pillars
described in detail in the project's ``vv-principles`` skill. The
table below summarises the boundaries each pillar can support; the
verification matrix in :doc:`matrix` resolves which pillar each
test claims, and :doc:`reference_solutions` gives the operator-form
taxonomy each pillar must commit to.

.. list-table::
   :header-rows: 1
   :widths: 22 18 18 22 20

   * - Pillar
     - Convergence-order
     - Flux-shape
     - Eigenvalue
     - When applicable
   * - Closed-form
     - ✓ (against exact)
     - ✓ (under assumptions)
     - ✓ (exact)
     - Limited regimes (homogeneous medium, simple geometry).
   * - **MMS** — Method of Manufactured Solutions
     - ✓ (great flexibility)
     - ✓ (any imposed shape)
     - **✗ (source-driven)**
     - Any operator that admits a non-vanishing trial solution. MMS
       does NOT prove eigenvalues by construction.
   * - Semi-analytical
     - ✓ (against arb-precision integral)
     - ✓
     - ✓
     - Hard cases with no closed form; the most common pillar in
       this project.

The semi-analytical pillar rests on a **two-step correctness
ladder**: integrator correctness (assumed for trusted upstream
libraries like ``scipy.integrate.quad`` / ``mpmath.quad``) and
**reduction correctness** (the math that takes the equation to a
single integral). If the reduction is wrong, the integral evaluates
exactly the wrong equation — a *reference contamination* instance.

ERR-038 (Atalay 1997 paper-precision floor, characterised in commit
``4c83e09``) is the project's canonical example of an evidence-class
distinct from a code bug: the published reference numbers are
themselves at a 1e-3-grade precision floor, so cross-checks against
them must respect that floor. ORPHEUS classifies this as
*reference-precision-floor* evidence — neither correct ("they
agree") nor incorrect ("the code is wrong"), but a publication-grade
ceiling that constrains the L1 claim. See ERR-038 in the error
catalog for the full diagnostic.


.. _verification-error-catalog-vocabulary:

Error catalog vs. paper-floor evidence
======================================

The ``ERR-NNN`` catalog mixes two evidence classes that the matrix
distinguishes:

- **Code bugs** (most ERR entries). A bug was committed, caught by an
  L0 / L1 test, fixed, and the test pinned with
  ``@pytest.mark.catches("ERR-NNN")`` to prevent regression.
- **Reference-precision floors** (ERR-038 family). The published
  numerical reference is the bottleneck — agreement at the published
  precision is achievable; tighter cross-checks would need a
  higher-precision reference, which does not exist in the literature
  at this time.

Both classes require an entry in the catalog. The distinction matters
when reading the matrix: a "verified" tag against a paper-floor
reference is a *bounded* claim, not an absolute one. Hardening the
claim requires either a structurally-independent reference at higher
precision (the canonical ORPHEUS path: Atkinson product Nyström
hardening of Atalay reference values, ERR-036) or a different
problem class.


.. _verification-three-meanings-cross-link:

Three-meanings taxonomy (where this verification suite consumes it)
====================================================================

The reference solvers in :doc:`/theory/reference_solvers` realise
three structurally-independent constructions of the Green's function
(Meanings α / β / γ). The verification matrix exploits this:

- (α) **Trajectory resolvent** — slab / cylinder / sphere / annulus /
  hollow sphere via :ref:`theory-trajectory-resolvent`.
- (β) **Spectral resolvent** — sphere reserved (gap; PS-1982 Eq. 21
  direct evaluator pending).
- (γ) **Singular-eigenfunction angular Green's** — slab / sphere /
  cylinder via :ref:`theory-singular-eigenfunction` (criticality)
  and :mod:`orpheus.derivations.continuous.fn_method` (interior flux
  reconstruction via KLL 1974).

When a problem is realised by all three constructions, **triple
agreement** (α) ≈ (β) ≈ (γ) is the highest-confidence L1 evidence
the project produces — agreement across three structurally-distinct
integrands. Production of (β) for sphere is the headline
infrastructure gap (literature memo:
``.claude/scratch/sanchez_chandrasekhar_gap.md``).


Pages
=====

.. toctree::
   :maxdepth: 2

   reference_solutions
   matrix
