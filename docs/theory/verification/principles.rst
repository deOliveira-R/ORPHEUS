.. _theory-verification-principles:

=======================
Verification Principles
=======================

.. This page is the designated single owner of the verification
   principles — the L0..L3 ladder, the evidence pillars, the
   reference tiers, the operator-form taxonomy, and how the four
   classification axes relate. The full treatment is authored in
   task #10 stage V5. The sections below were moved verbatim from
   the dissolved ``docs/verification/index.rst`` front page so the
   part restructure (stage V2) loses nothing.


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

The reference solvers in :doc:`/theory/references/index` realise
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


.. The three sections below were moved verbatim at task #10 stage V3
   from the Overview / Verification Methodology / Reference Case
   Types chapters of the dissolved ``docs/theory/verification.rst``.
   V5 integrates them into the authored principles treatment.


Overview
========

ORPHEUS verifies every transport solver against analytical reference solutions
derived from each method's own mathematical equations using SymPy.  Each
derivation is **self-contained**: it starts from the solver's formulation and
derives the expected eigenvalue independently.  No cross-verification (comparing
one solver against another) **stands in for** this evidence — solver-vs-solver
agreement is level **L4 (informational)**, strictly distinct from the L0–L3
correctness ladder.  Each solver's L0–L3 verification stands on its own, as if
every other solver were deleted; the separate cross-method regression protocol
(:doc:`/theory/verification/cross_method`) then adds L4 agreement gates on top,
never in place, of that ladder.

The same code that produces the LaTeX equations in this part also produces
the reference values consumed by the ``pytest`` test suite.  This is the
**single source of truth**: equations in the documentation cannot drift from
the values in the tests because both come from the same ``orpheus/derivations/``
package.


Verification Methodology
========================

Each verification case defines three things:

1. **Cross sections** — synthetic macroscopic data for abstract regions
   (not specific materials).  This isolates the numerical method from
   the cross-section processing pipeline.

2. **Analytical eigenvalue** — derived from the solver's own equations
   using SymPy (symbolic where possible, numerical for special-function
   kernels like E₃ and Ki₄).

3. **Tolerance** — principled, not arbitrary.  The tolerance is bounded
   by the dominant error source:

.. list-table:: Tolerance Rationale
   :header-rows: 1
   :widths: 20 15 20 45

   * - Method
     - Tolerance
     - Error source
     - Rationale
   * - Homogeneous
     - < 10⁻¹²
     - FP arithmetic
     - Direct ``numpy.eigvals`` of small dense matrix
   * - CP slab
     - < 10⁻⁶
     - Power iteration
     - Solver keff_tol=10⁻⁷ bounds the error
   * - CP cylinder
     - < 10⁻⁵
     - Power iteration + Ki₄ interpolation
     - Solver keff_tol=10⁻⁶ plus Ki₄ table resolution
   * - SN (homogeneous)
     - < 10⁻⁸
     - Power iteration
     - Flat flux is exact in DD; only iteration error
   * - SN (heterogeneous)
     - O(h²)
     - Spatial discretisation
     - Richardson extrapolation from own mesh convergence
   * - MOC (homogeneous)
     - < 10⁻⁴
     - Ray spacing + iteration
     - Flat-source exact; convergence limited by ray density
   * - MC
     - z < 5σ
     - Statistical
     - Central Limit Theorem; σ ~ 1/√N_active
   * - Diffusion
     - O(h²)
     - FD spatial discretisation
     - Analytical buckling eigenvalue (bare) or Richardson (reflected)


Reference Case Types
====================

Analytical
----------

The eigenvalue is computed in closed form or as the eigenvalue of a
finite-dimensional matrix derived symbolically by SymPy.  No numerical
solver is involved.

**Examples**: homogeneous infinite medium (matrix eigenvalue), diffusion
bare slab (buckling eigenvalue), SN/MOC/MC homogeneous (each derived from
the solver's own equations showing that the infinite-medium eigenvalue
is the exact solution).

Semi-Analytical
---------------

The reference involves a special function (E₃, Ki₃/Ki₄) evaluated to
integrator precision, followed by a finite matrix eigenvalue.  The only
approximation is the numerical :term:`quadrature` for the special function, which
is controlled to machine precision via ``scipy.special.expn`` (E₃) or
high-resolution lookup tables (Ki₃/Ki₄ with 20,000 points).

**Examples**: all CP slab and CP cylinder cases.  The CP matrix is
exact for the collision probability formulation; the eigenvalue is
a finite matrix problem.

.. _richardson-extrapolation:

Richardson-Extrapolated
-----------------------

For discretised solvers on heterogeneous problems, the reference is the
converged limit of the solver itself, estimated by running at 4 mesh
refinement levels and extrapolating assuming O(h²) convergence:

.. math::
   :label: richardson-extrapolation-formula

   k_{\rm ref} \approx k_h + \frac{k_h - k_{2h}}{2^p - 1}

where :math:`p = 2` for :term:`diamond-difference <diamond difference>` and finite-difference schemes.

This is legitimate formal verification of the **convergence rate** — it
tests that the implementation converges at the theoretically expected order,
even though the reference value is self-generated.

**Examples**: SN heterogeneous slab, MOC heterogeneous pin cell, diffusion
fuel+reflector.
