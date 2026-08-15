.. _theory-verification-principles:

=======================
Verification Principles
=======================

.. This page is the designated single owner of the verification
   doctrine — the claim-layer taxonomy, the normative L0..L3 ladder
   (+ foundation + the L4 ruling), the three reference pillars, and
   the map relating the classification systems. Authored at task #10
   stage V5 from the ``vv-principles`` skill + its reference.md,
   integrating the sections parked here at stages V2/V3. The
   operational twin of this page is the ``vv-principles`` skill
   (``.claude/skills/vv-principles/``): same doctrine, packaged as
   the agent-side decision instrument (NEVER→instead redirects,
   review recipes, the live failure-mode registry). New failure
   modes and anti-patterns land in the skill first; this page
   carries the doctrine and its rationale.

Verification is the discipline that separates *"the code produces
plausible numbers"* from *"the code solves its governing equations
correctly."* In a codebase where numerical kernels are written and
reviewed by AI agents, that separation cannot rest on care and
plausibility — the characteristic failure modes are exactly the
*plausible* substitutions (a sign flip, a swapped subscript, a
missing :math:`2\pi`) that survive every "looks right" review. It
rests instead on a small set of principles about **what a claim
says**, **what evidence can support it**, and **where a reference
solution is allowed to come from**.

This page is the normative home of those principles. The
:doc:`harness` chapter owns the *mechanics* (how a test declares
its level, how declarations are audited); :doc:`reference_solutions`
owns the *reference contract* (vocabulary, operator forms, tiers,
the ``ContinuousReferenceSolution`` API); the per-method chapters
own the *evidence*; :doc:`matrix` renders the coverage. This page
owns the doctrine they all instantiate.

.. contents::
   :local:
   :depth: 1


Verification and validation — which screw is being turned
=========================================================

Following Oberkampf & Roache: **verification** proves *code matches
mathematics*; **validation** proves *mathematics matches nature*.
Two distinct screws.

The distinction is load-bearing in practice. If the equation itself
is wrong — a missing term in the discretisation, a wrong recursion
in the angular sweep — verification can pass cleanly: the code
faithfully matches the (wrong) mathematics. Only validation (L3)
catches it, as a systematic error that refuses to close against
experiment no matter how the verification cases are tightened.
Conflating the two produces false confidence in both directions:
"all my verification tests pass" does not mean "my equations are
right," and "my code agrees with the experiment" does not mean "my
code is verified." Each screw needs its own evidence.

The binding vocabulary contract — what the words *verification*,
*validation*, *reference solution*, and *benchmark* are allowed to
mean in this repository — lives at :ref:`vv-vocabulary`.


.. _verification-structural-independence:

Structural independence — the root principle
============================================

Everything else on this page inherits from one insight: **a
verification test proves correctness only if it exercises a
different mathematical path than the code under test.** When test
and code share an upstream identity, an integrand, or a primitive,
an error in that shared element is invisible — both sides agree,
and both are wrong.

The canonical instance is ERR-032 (all ``ERR-NNN`` entries live in
the error catalog,
``.claude/skills/vv-principles/error_catalog.md``). Two
"independent" cross-checks of a Peierls kernel both applied the
identity :math:`\int E_2 = 1 - E_3`. The correct identity is
:math:`\int E_2 = \tfrac{1}{2} - E_3`. The two implementations
agreed to :math:`10^{-39}` — not because either was right, but
because they shared the upstream identity. The bug surfaced only
when a third reference, derived from a structurally different angle
(a particle-balance row-sum), refused to close.

The distinction that names the phenomenon:

- Two derivations are **procedurally independent** if they use
  different code paths, different functions, different files.
  Procedural independence catches typos, copy-paste corruption, and
  index drift *between the two paths*.
- They are **structurally independent** only if they exercise a
  different integrand, a different identity, or a different
  recursion. Structural independence catches errors in the
  *content* — wrong identity, wrong integrand, wrong constant.

Procedural independence is not structural independence. Two
derivations of the same closed form are procedural independence in
a costume.

**Operational corollary — ship references in pairs.** When shipping
a new analytical reference, force the cross-check to come from a
different structural angle: a *kernel check* (row-sum
:math:`\sum_j P_{ij} = 1`, particle balance, positivity bounds)
**and** a *closed form* (eigenvalue, asymptotic limit, degenerate
parameter regime) — two structurally distinct angles on the same
physics, never two derivations of the same expression.

The principle applies to *all* test design, not just numerical
cross-checks. ERR-051 is the test-design instance: a
contract-validation helper asserted the wrong invariant
(:math:`\Pi R = I` instead of :math:`\Pi R = 4\pi I`), and its sole
test fed it a deliberately-broken input constructed precisely so
the wrong assertion succeeded at raising. A negative-only test of
an ``assert_X`` method validates the *raising behaviour*, not the
*invariant claim* — every contract validator needs at least one
positive test (correct instance, must not raise) and one negative
test (broken instance, must raise), with the broken instance built
independently of the assertion under test.


.. _verification-claim-layers:

What is being claimed — the claim-layer taxonomy
================================================

Claims about a solver output are layered, and each layer adds
dependencies. Verify lower layers before higher ones, and match
evidence to the *claim's* layer::

                 ┌────────────────────────────────┐
                 │  Eigenvalue claim              │  depends on eigenvalue solver
                 │  (k_eff, k_inf)                │  + flux shape + discretisation
                 └────────────────────────────────┘
                               ↑ depends on
                 ┌────────────────────────────────┐
                 │  Flux-shape claim              │  depends on the discrete model
                 │  (ψ(x, μ, E), φ(x))            │  + boundary conditions
                 └────────────────────────────────┘
                               ↑ depends on
                 ┌────────────────────────────────┐
                 │  Convergence-order claim       │  pure math; lowest dependency
                 │  (O(h^p), MMS slope)           │  verifies parts AND whole
                 └────────────────────────────────┘

Three reclassifications to apply whenever reading a claim:

- **Convergence-order results are math claims, not solver claims.**
  They prove the discretisation is *consistent* — nothing about the
  solved value being correct. A solver can converge at exactly the
  theoretical order to the wrong limit (ERR-006 did; see the
  :ref:`failure-mode catalogue <verification-failure-modes>`).
  MMS lives at this layer.
- **Flux-shape results are model claims, not eigenvalue claims.**
  They depend on the equation and the boundary conditions, not on
  the eigenvalue iteration. MMS reaches this layer when the
  manufactured source is structurally independent of the code's
  primitives.
- **Eigenvalue results are solver claims.** They bring the
  iteration scheme, the normalisation, and the convergence test
  into scope. MMS does *not* reach this layer — k-eigenvalue
  verification needs an analytical eigenvalue or a
  structurally-independent semi-analytical reference (see
  :ref:`the pillars <verification-evidence-classes>`).


.. _vv-level-ladder:

The V&V ladder — how deep the evidence goes
===========================================

ORPHEUS classifies correctness evidence on a four-level ladder,
with one orthogonal bucket and one parallel track::

   VERIFICATION — "Are we solving the equations right?"
     L0  Term verification        hand calculation vs code, per term
     L1  Equation verification    analytical solutions, MMS, convergence order
     L2  Integration              multi-group + heterogeneous, self-convergence

   VALIDATION — "Are we solving the right equations?"
     L3  Validation               experimental data (ICSBEP, IRPhE, SINBAD)

   INFORMATIONAL — parallel to the ladder, never a rung
     L4  Benchmarking             code-to-code agreement

   ORTHOGONAL TO THE LADDER
     foundation                   software invariants with no physics equation

Every test in ``tests/`` declares its rung via ``pytest`` markers;
the declaration mechanics, the tagging-precedence chain, and the
audit that enforces them are the :doc:`harness` contract. What
follows is the *meaning* of each rung — what it proves, and just as
importantly, what it deliberately does not prove.

L0 — term verification
----------------------

Every term of the discretised equation, isolated and verified by
hand calculation against the code's output — right sign *and* right
magnitude. Evidence at this level zeroes out every term except the
one under test (via boundary condition, material choice, or
geometry), checks both polarities of any sign-bearing term, and
uses non-uniform profiles (uniform profiles hide index errors).

L0 does **not** prove that the discretised equation is the right
equation, that the assembled terms produce the right global result,
or that the iteration converges. It proves the bricks, not the
house — which is precisely why it is the only defense against the
plausible-substitution failure modes catalogued
:ref:`below <verification-failure-modes>`.

L1 — equation verification
--------------------------

The full discretised operator matches an analytical solution, an
MMS-manufactured solution, or a high-precision semi-analytical
reference — with the convergence order *measured* against the
theoretical order of the discretisation, not merely "passes a
tolerance."

L1 does **not** prove that heterogeneous problems are correct, that
multi-group coupling is correct, or that the solver handles
realistic configurations. Those are L2's job.

L2 — integration
----------------

Multi-group, heterogeneous, mesh-refined problems: the full solver,
in a realistic configuration, self-converges to a stable answer as
spatial and angular resolution increase, and its converged values
agree with structurally-independent references where they exist
(the flat-flux degenerate limit, a transfer-matrix eigenvalue, an
MMS-pinned operator). Self-convergence at the expected order is L2
evidence *about stability and consistency*; the converged **value**
still inherits its correctness from the L1 references below it —
a self-generated limit is never its own reference (see the
:ref:`retired Richardson class <richardson-extrapolation>`).

L2 does **not** prove the solver matches reality, nor that the
cross-section data is right.

L3 — validation
---------------

The solver, with realistic cross sections, matches experimental
measurements (ICSBEP criticality benchmarks, IRPhE reactor-physics
experiments, SINBAD shielding benchmarks) within acceptance
criteria stated **before** the comparison. L3 is sequenced, not
aspirational: it starts after L1/L2 maturity, because L3 agreement
over an unverified L2 base is accidental agreement.

The necessity chain
-------------------

Each rung requires the rungs below it:

- L1 without L0 = compensating errors (two wrong terms cancelling
  in the case you happened to test).
- L2 without L1 = masked components (the integration test averages
  over the component error).
- L3 without L2 = accidental agreement (experiment matched through
  cancelling model errors or tuned data).
- L4 without L0–L2 = proves nothing at all.

Foundation tests — orthogonal, not a rung
-----------------------------------------

Before placing a test on the ladder, ask: *is this
physics-equation-tied or software-invariant?* Foundation tests
verify software invariants — data-structure contracts, factory
outputs, algebraic-reduction identities of pre-physics building
blocks. They have no equation ``:label:`` in the theory corpus to
point at, so they carry ``@pytest.mark.foundation``, never a
``verifies(...)`` marker, and the audit reports them in their own
bucket. Promoting a physics test to "foundation" to dodge the
ladder is the failure mode the bucket's tagging rules exist to
prevent — the taxonomy, tiebreak rules, and anti-patterns live at
:ref:`vv-foundation-tests`.

.. _verification-l4-ruling:

L4 — benchmarking is parallel, never a rung
-------------------------------------------

Two ORPHEUS solvers agreeing is *cross-implementation agreement*,
not *correctness evidence* — both could share a textbook error, a
convention, or a copy-pasted lineage. L4 therefore sits **parallel
to** the ladder, never on it: it produces information about whether
two implementations agree and zero information about whether either
is correct. Every L4 claim must name its L0–L3 backing.
"ORPHEUS-CP agrees with ORPHEUS-MC" is not a verification
statement; "ORPHEUS-CP agrees with ORPHEUS-MC, both independently
verified at L1 against analytical references" is.

One refinement, from the :doc:`cross-method protocol
<cross_method>`: when two methods are each independently
L1-verified *and* their reference constructions are structurally
independent, their pairwise agreement gate is tagged **L1** — not
because agreement became correctness evidence, but because the gate
is a *consumer* of the two existing L1 references: a drift in
either method off its own L1-pinned truth reddens the pair gate, so
the gate carries L1-strength **regression** power. The *primary* L1
evidence remains each method's own truth gate; the agreement is
never a stand-in for it. That is the whole L4 ruling in one
sentence: agreement gates may *inherit* the strength of the
references they consume, but they never *generate* correctness
evidence of their own.

.. _verification-1g-degeneracy:

1-group degeneracy — the canonical warning
------------------------------------------

:math:`k = \nu\Sigma_f / \Sigma_a` is flux-shape independent: a
1-group eigenvalue is a material-property ratio, computable without
solving the transport equation at all. A 1-group eigenvalue test
therefore cannot detect *any* error in the spatial, angular, or
scattering operators. The ORPHEUS cylindrical diamond-difference
recurrence bug (ERR-006) survived twenty passing tests before a
multi-group run caught it. **Multi-group (≥ 2 groups) is mandatory
for any test claiming to verify transport.** This is an instance of
a general lens — the measured functional's invariance group
contains the error class (Mode 12,
:ref:`below <verification-failure-modes>`) — degenerate enough to
deserve its own named rule.


.. _verification-evidence-classes:

Where references come from — the three pillars
==============================================

Every verification reference belongs to exactly one of three
pillars. The pillar structure is not bureaucracy: it enumerates the
only three ways a reference can be derived *independently of the
discretisation under test*, captured by a duality plus a fallback:

- **"Given an equation, find the solution"** → **closed-form**
  analytical. The classical task: governing equation + boundary
  conditions + assumptions → an exact expression; the eigenvalue
  falls out of the boundary problem.
- **"Given a solution, find the equation's source"** → **MMS**
  (Method of Manufactured Solutions). The inverse task: impose
  :math:`\psi_{\text{chosen}}`, substitute into the operator,
  derive the source :math:`Q^{\text{ext}}` that makes it exact.
  You buy flexibility (any imposed shape); you sell the eigenvalue
  (the problem is source-driven by construction).
- **When neither question closes algebraically:** reduce the
  equation to a single integral and evaluate to arbitrary
  precision → **semi-analytical**.

Closed-form and MMS are both *analytical* (exact by construction);
semi-analytical is *exact via arbitrary-precision numerics*. The
distinction matters because each pillar has a different evidence
boundary:

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

Closed-form analytical
----------------------

The reference is derived from the governing equation along a
different path — symbolically in SymPy or by hand. Symbolic ≠
discrete by construction, so a correctly-derived closed form is
ground truth in the limit. ORPHEUS instances that close: the
homogeneous infinite medium (matrix eigenvalue
:math:`k_\infty = \lambda_{\max}(\mathbf A^{-1}\mathbf F)` — the
shared L1 anchor, :doc:`homogeneous`), bare-slab diffusion (sine
eigenfunction with buckling), flat-flux transport limits derived
from each solver's own equations. Heterogeneous and curvilinear
closed forms are rare — which is why the other two pillars exist.

The derivation itself must be verified: ERR-032 (above) is a
closed-form reference that was *wrong*, because its derivation
shared an upstream identity with the code it was meant to check.

MMS — and why its rules are what they are
-----------------------------------------

MMS is constructive: choose :math:`\psi_{\text{chosen}}`, apply the
continuous operator symbolically, obtain :math:`Q^{\text{ext}}`,
run the code with that source, and require the code to reproduce
:math:`\psi_{\text{chosen}}` under refinement at the theoretical
order. Each operational rule exists for a mechanical reason:

- **The trial solution must not vanish under derivatives.**
  Polynomials die at finite derivative order; whatever the operator
  does beyond that order becomes invisible. Trigonometric and
  exponential trials self-reproduce under differentiation and keep
  every order observable.
- **The trial solution must be non-trivial at boundaries.** A
  :math:`\psi` that vanishes at the boundary satisfies the BC
  equation trivially and tests nothing about the BC
  implementation.
- **The manufactured source must be structurally independent of
  the code's primitives.** If :math:`Q^{\text{ext}}` is computed by
  the same numerical kernels the solver uses, the test is a
  tautology. Derive it symbolically.
- **The trial must stress the discretisation, not minimise source
  complexity.** Human-written MMS trends toward simple trials
  because hand-deriving :math:`Q^{\text{ext}}` is error-prone —
  a heuristic that protects the human, not the verification. An
  agent deriving the source in SymPy has no such constraint and
  must override the inherited bias deliberately: mixed scales,
  higher-frequency modes, non-trivial group coupling that
  activates the scattering off-diagonals, angularly non-trivial
  shapes in curvilinear contexts. A simple trial function tests
  less (Mode 7, :ref:`below <verification-failure-modes>`).
- **MMS proves no eigenvalue — mechanically, not as a
  limitation.** The problem is source-driven: you imposed the
  solution, the derived source makes it exact, and the eigenvalue
  is whatever you started with. There is no eigenvalue information
  in the construction. Eigenvalue claims must be matched to a
  closed-form or semi-analytical reference.

Semi-analytical — the two-step correctness ladder
-------------------------------------------------

The reduction from governing equation to single-integral form is
analytical; the evaluation is numerical at arbitrary precision.
Correctness rests on a two-step chain:

1. **Integrator correctness** — assumed for well-tested upstream
   libraries (``scipy.integrate.quad``, ``mpmath.quad``); a
   verification obligation in its own right for custom integrators
   (e.g. the Bickley–Naylor kernels in
   :mod:`orpheus.derivations.common.kernels`).
2. **Reduction correctness** — the load-bearing math. If the
   reduction is wrong, the integrator answers the wrong question
   exactly; convergence at the right order to the wrong limit is
   the signature (ERR-006).

ORPHEUS instances: the E₃-based collision-probability eigenvalue
grids and the Peierls integral references
(:mod:`orpheus.derivations.continuous.peierls_nystrom`, Nyström
collocation at 30+ digits) — the most common pillar in this
project.

Ancillary references — never pillars
------------------------------------

Three familiar reference *uses* are not evidence sources:

- **Independent re-derivation** — a different mathematical path to
  the same closed form. Strong only when the paths are
  *structurally* independent (different identity, different
  integrand); otherwise it is ERR-032 waiting to happen.
- **Code-to-code (L4)** — cross-implementation agreement,
  exclusively. Two codes agreeing share whatever they share —
  possibly the error. See the :ref:`L4 ruling
  <verification-l4-ruling>`.
- **Monte Carlo** — itself a numerical method whose geometry
  tracker, free-flight sampler, and tally estimators need prior
  verification. Once verified, MC is a *consumer* of references,
  never a source; a deterministic solver compared against an
  unverified MC tally is an L4 claim dressed as L1 evidence.


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


Reference tiers and operator forms — the other two axes
=======================================================

Two further classification systems, both owned by the
:doc:`reference contract page <reference_solutions>`, grade and
scope individual references:

**Tiers** grade a reference value's *provenance* — how structurally
independent it is from any solver: T1 (closed-form symbolic), T1.5
(analytical MMS — exact by construction, source-driven), T2
(semi-analytical to root-finder tolerance), T2.5 (semi-analytical
but sharing a discretisation family with the solver under test —
tolerated, scheduled for replacement), and the two **banned**
grades: T3 (the solver under test used as its own reference) and T4
(a different solver used as the reference). The full definitions
and the per-module audit ledger live at
:ref:`verification-campaign-migration`. The tier system is
:ref:`structural independence <verification-structural-independence>`
made enforceable for references: the T3/T4 bans are what keep L2
self-consistency from being dressed up as L1 equation verification.

**Operator forms** scope which equation a reference actually
solves. Every :class:`~orpheus.derivations.ContinuousReferenceSolution`
commits to exactly one operator-form tag (``"differential-sn"``,
``"integral-peierls"``, ``"diffusion"``, …) and consuming tests
assert that the target solver discretises the same form — a
reference for one equation has no business verifying a solver for
another. The taxonomy lives at :ref:`operator-form-taxonomy`. This
is the type-system defense against the reference-contamination
sub-case "the reference solved a different problem than you
thought" (:ref:`below <verification-reference-contamination>`).


.. _verification-classification-map:

How the classification systems relate
=====================================

The corpus carries five classification systems — claim layers, the
ladder, pillars, tiers, operator forms. They are not rivals; they
answer five different questions about one verification claim:

.. list-table::
   :header-rows: 1
   :widths: 22 38 40

   * - System
     - Question it answers
     - Owner
   * - Claim layer
     - *What* is being claimed? (convergence-order / flux-shape /
       eigenvalue)
     - :ref:`this page <verification-claim-layers>`
   * - Ladder rung
     - *How deep* is the evidence? (L0–L3; foundation orthogonal;
       L4 parallel)
     - :ref:`this page <vv-level-ladder>`
   * - Pillar
     - *Where does the reference come from?* (closed-form / MMS /
       semi-analytical)
     - :ref:`this page <verification-evidence-classes>`
   * - Tier
     - *How trustworthy is the reference value?* (T1 … T2.5;
       T3/T4 banned)
     - :ref:`reference contract <verification-campaign-migration>`
   * - Operator form
     - *Which equation does the reference solve?*
     - :ref:`reference contract <operator-form-taxonomy>`

A complete verification claim fixes all five coordinates, and the
systems constrain each other:

- The **pillar bounds the claim layer**: MMS supports
  convergence-order and flux-shape claims, never eigenvalue claims
  (the capability table above).
- The **tier enforces the ladder's necessity chain**: a T3
  reference (solver as its own truth) collapses L1 into L2
  self-consistency; the ban is what makes an "L1" tag mean
  something.
- The **operator form guards the pillar's premise**: a reference
  is only ground truth *for the equation it solves*; the form tag
  makes consuming a mismatched reference unspellable.
- The **claim layer selects the ladder evidence**: an eigenvalue
  claim needs L1 evidence from an eigenvalue-capable pillar
  (closed-form or semi-analytical); an order claim can rest on MMS
  alone.

**Worked example.** The diffusion 2-region fuel+reflector case
(:doc:`diffusion`): the claim is an *eigenvalue* (+ flux shape) —
the top claim layer; the evidence is *L1* (equation verification
against an external truth); the reference is *semi-analytical*
(transcendental transfer-matrix root, solved by a double-precision
Brent root-find) at tier *T2*, committed to operator form
``"diffusion"``; and the consuming
tests in :mod:`tests.diffusion.test_continuous_reference` assert
the form match and carry the L1 marker. Five coordinates, one
claim, no ambiguity about what has been proven.


.. _richardson-extrapolation:

Richardson extrapolation — a retired reference class
====================================================

For a discretised solver with a known convergence order :math:`p`,
running at mesh spacings :math:`h` and :math:`2h` lets the
converged limit be *estimated* from the solver's own refinement
sequence:

.. math::
   :label: richardson-extrapolation-formula

   k_{\rm ref} \approx k_h + \frac{k_h - k_{2h}}{2^p - 1}

.. (vv-status rationale) historical-record definition: the Richardson
.. extrapolation formula, a DELIBERATELY RETIRED reference class (tier T3,
.. methodologically banned — the solver-as-its-own-reference). No live test
.. exists by design; the label documents the technique and why it was
.. replaced.
.. vv-status: richardson-extrapolation-formula documented

with :math:`p = 2` for :term:`diamond-difference <diamond
difference>` and finite-difference schemes. The formula is sound
mathematics for what it says: an error-cancellation estimate of the
limit of *this solver's own sequence*, and the backbone of measured
convergence-order assertions (an L2 self-convergence mechanic).

What it can never be is a **reference**. A Richardson-extrapolated
"truth" is the solver under test used as its own reference — tier
**T3, methodologically banned** (:ref:`the tier ledger
<verification-campaign-migration>`). It verifies the convergence
*rate* while assuming precisely what verification must prove: that
the limit being converged to is the right one. ERR-006 is the
standing demonstration — a curvilinear sweep converging at exactly
:math:`O(h^2)` to a wrong limit, invisible to any self-referential
gate. Early ORPHEUS used Richardson-extrapolated heterogeneous
references; the verification campaign retired every one of them in
favour of structurally-independent pillars:

- **SN heterogeneous** — retired at campaign Phase 2.1a; replaced
  by analytical MMS continuous references (T1.5, six references)
  plus the T2 transfer-matrix eigenvalue
  (:doc:`sn`).
- **MOC heterogeneous** — retired at campaign Phase 2.2a; replaced
  by the pin-cell MMS continuous reference (T1.5;
  :ref:`moc-mms-verification`).
  :mod:`orpheus.derivations.continuous.cases.moc` now derives
  homogeneous cases only.
- **Diffusion 2-region** — the Richardson reference served through
  the migration window, then retired for the transcendental
  transfer-matrix reference (T2) once the operator-algebra solver
  landed; the historical record, including the agreement of the
  retired value with its successor at the Richardson floor, is
  preserved at :ref:`the diffusion chapter's retirement record
  <diffusion-2region-richardson>`.

The equation above is kept for the record and for its live use in
measured-order assertions; no ORPHEUS reference value is
Richardson-extrapolated today.


The verification-case contract
==============================

Every verification case defines three things:

1. **Cross sections** — synthetic macroscopic data for abstract
   regions (not specific materials). This isolates the numerical
   method from the cross-section processing pipeline.

2. **Reference value** — derived from the solver's own governing
   equations by pure mathematics (SymPy where symbolic, mpmath for
   special-function kernels), through one of the three pillars.
   Each derivation is **self-contained**: it starts from the
   method's formulation and derives the expected value
   independently, as if every other solver were deleted.

3. **Tolerance** — principled, not arbitrary: bounded by the
   dominant error source of the case.

.. list-table:: Tolerance rationale by dominant error source
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
     - The gates pass ``keff_tol`` **explicitly** at 10⁻⁷–10⁻⁸; the
       shipped default is 10⁻⁶ (``CPParams.keff_tol``). ⛔ An increment
       test does not *bound* the error — see the cross-family note below
   * - CP cylinder
     - < 10⁻⁵
     - Power iteration + cylinder-kernel interpolation
     - Same explicit gate tolerances, plus the interpolant's fixed
       precision (⛔ the 20 000-point Ki lookup table was superseded by a
       fixed-precision Chebyshev interpolant at Phase B.4, #94 —
       ``orpheus/cp/solver.py:68-70``)
   * - SN (homogeneous)
     - < 10⁻⁸
     - Power iteration
     - Flat flux is exact in DD; only iteration error
   * - SN (heterogeneous)
     - case-specific
     - Spatial discretisation
     - MMS measured-order gates (T1.5) + transfer-matrix eigenvalue
       (T2); tolerance set by the discretisation error at the
       tested refinement, never by a self-generated limit
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
     - Analytical buckling eigenvalue (bare, T1) or transcendental
       transfer-matrix reference (2-region, T2)

The same code that produces the LaTeX equations in this part also
produces the reference values consumed by the ``pytest`` suite —
the **single source of truth**: equations in the documentation
cannot drift from the values in the tests because both come from
the same :mod:`orpheus.derivations` package.


.. _verification-user-path:

Route through the user's path — and treat a bypass as a gap
===========================================================

A verification case must reach the system under test **the way a user
reaches it**. Stated as the ruling that established it (2026-08-05):

   Tests must route through the machinery that a user would exercise
   without bypassing code functionality. Or else it's not testing the
   path the users go through.

This is not a style preference about test ergonomics. A test that
assembles an internal object by hand and injects it past the layers a
user traverses has verified **that object**, and has said nothing about
the declaration, resolution, and wiring that were supposed to produce
it. Every one of those skipped steps is production code that ships
unexercised — and the failure is silent in the worst way, because the
suite is green and the coverage report counts the injected object's
lines.

The canonical instance in this codebase is the boundary source. The
affine boundary law is

.. math::

   \gamma_-\psi \;=\; L\,\gamma_+\psi \;+\; q ,

and a manufactured-solution driver needs a non-zero :math:`q`. Supplying
it as :meth:`~orpheus.transport.source_sinks.AngularBoundarySourceSink.prescribed_inflow`
straight into a composite source verifies the *channel* — the RHS carries
it, the sweep consumes it — while saying nothing about whether **declaring**
a :class:`~orpheus.geometry.boundary.PrescribedInflow` produces it. It did
not: until the mesh-BC bridge landed, a declared inflow realized into an
operator nothing consumed, so the declaration was inert and no green test
could tell. The honest driver declares the boundary condition and lets the
machinery materialise :math:`q`.

Test-owned machinery is a GAP SIGNAL
------------------------------------

The rule has a second half, and it is the more useful one.

A test may legitimately implement a production **Protocol** — an
:class:`~orpheus.geometry.boundary._source.InflowSourceSpec` evaluating a
manufactured solution on a face is using the machinery, not bypassing it.
But when a test has to *build* something in order to exercise production,
that is evidence the production shape is **missing**, not evidence the
test is clever:

   If the MMS is exercising custom machinery that is not part of
   production, then the proper shape of the machinery should be
   implemented so that the MMS can use it. It is a sign of a gap.

So a test-authored implementation is a **stopgap with an owner**: it
ships, and it ships with a phase that promotes the shape into production.
The diagnostic question at review time is not *"is this test allowed to
define that class?"* but *"why does production not already offer it?"* —
and the answer is nearly always that a real capability was never built,
with the test's private version quietly standing in for it.

Two corollaries worth stating, because both have bitten:

* **A bypass and a stopgap look identical in the diff.** Both add a class
  to a test module. The discriminator is whether the object is *reached
  through* the production path (stopgap — the Protocol is public and the
  wiring is exercised) or *injected past* it (bypass). Ask which
  production lines run because of it.
* **A convenience factory can BE the bypass.** ``prescribed_inflow`` is
  production code, public and correct; using it in a driver still skips
  the declaration tier. "It is a production symbol" does not establish
  that the user path was travelled.

.. _verification-principled-equivalence:

Bit-identity vs principled equivalence
======================================

Bit-identity is an *implementation* property, not a math property.
An ``np.array_equal`` regression contract is a strong gate while
the implementation is unchanged — verification by inheritance from
a previously-verified value. The same gate becomes the wrong gate
when a refactor deliberately changes the floating-point reduction
tree: the two implementations compute the same value in real
arithmetic and disagree at IEEE-754 ULP level because addition is
not associative.

A non-bit-exact change is acceptable **only when all three** of the
following hold; reject if any fails:

1. **The new formulation is principled at every step** — each
   intermediate is a named, inspectable domain quantity (a
   per-group production rate), not "whatever the reduction order
   happened to produce."
2. **The new value is verified against a structurally-independent
   reference.** Old-vs-new ULP distance is necessary but never
   sufficient — both could be wrong by the same offset.
3. **The drift is FP-non-associativity, dimensionally
   explainable** — bounded by (reduction depth) × ULP, or for
   iterative solvers (iteration count) × (condition number) × ULP.
   Drift beyond these bounds is an algorithmic change masquerading
   as FP noise.

When all three hold, narrow the regression contract for the touched
primitive (``assert_array_almost_equal_nulp``) and preserve
bit-identity elsewhere; document the relaxation.

A companion discipline guards the *diagnosis* step: an error
measured **in isolation** (an offline residual, a per-kernel
discrepancy) is not automatically "the dominant floor" or "the
improvement this change buys." Internal self-consistency is
necessary, never sufficient — a matvec≡sweep round-trip at
:math:`10^{-16}` proves two paths solve the *same* operator, not
that the operator's fixed point is correct (ERR-061). Before
crediting an isolated error as the floor, it must survive three
end-to-end checks: an **end-to-end swap** (wire the piece into the
full solver; the claimed effect persists in the converged answer),
a **term-silent control** (zero the term; the answer is
byte-identical where the term should not matter), and
**amplification** — the sharpest disproof: grow the term 3–10× and
require the converged answer to get strictly *worse* against a
structurally-independent reference. A term whose amplification
moves nothing was never the floor.


.. _verification-failure-modes:

How verification fails — the failure-mode catalogue
===================================================

Solver-bug modes (1–6) — the code is wrong
------------------------------------------

The characteristic failure modes of AI-generated numerical code are
mechanically explainable, not arbitrary. Sub-word tokenizers chunk
text into reusable units; symbols that look similar —
:math:`\Sigma_a` and :math:`\Sigma_f`, :math:`\mu_x` and
:math:`\mu_y`, :math:`E_2` and :math:`E_3` — share or co-locate
tokens, sitting adjacent in embedding space, and probabilistic
decoding lands on a wrong-but-plausible substitution at rates far
above human typo rates. The six modes are the observable signature:

.. list-table::
   :header-rows: 1
   :widths: 5 20 33 42

   * - #
     - Mode
     - Example
     - Detection (L0 strategy)
   * - 1
     - Sign flip
     - :math:`(a-b)` vs :math:`(b-a)`
     - Heterogeneous eigenvalue diverges under refinement
   * - 2
     - Variable swap
     - ``mu_x`` vs ``mu_y``; ``SigS`` vs ``SigS^T``
     - Per-ordinate flat-flux residual; asymmetric 2G inputs
   * - 3
     - Missing factor
     - Missing ``ΔA/w``, ``2π``, volume
     - Fixed-source flux spike at r=0 vs ``Q/Σ_t`` analytic
   * - 4
     - Wrong recursion
     - :math:`\alpha_{m+1/2}` index drift
     - Per-ordinate flat-flux residual
   * - 5
     - Index error
     - ``face[i]`` vs ``face[i+1]``
     - Non-uniform mesh produces detectably different keff
   * - 6
     - Convention drift
     - Definition site vs usage site disagree
     - 2G heterogeneous with asymmetric SigS — wrong group ratio

The same defense — L0 term verification — also protects against
every other tokenized-error generator: human subscript typos, OCR
drift (:math:`\Sigma` vs :math:`E`), copy-paste index corruption
when adapting a sibling block, and convention drift between files
written months apart. Every caught instance is catalogued as an
``ERR-NNN`` entry (:ref:`the error-catalog contract
<verification-error-catalog-vocabulary>`).

.. _verification-test-design-modes:

Test-design modes (7–12) — the test cannot see the bug
------------------------------------------------------

The second family is subtler: the solver bug is real, but the test
is *structured* so it cannot mathematically observe it. The defense
is test review at design time, not more L0. Each mode below names
its mechanism, its canonical ORPHEUS instance, and its defense; the
operational recipes (mutation protocols, sentinel patterns) live in
the ``vv-principles`` skill.

**Mode 7 — MMS simplification bias.** The ansatz *nulls* the
hardest term by design: an isotropic-in-:math:`\mu` curvilinear
trial makes the angular-redistribution term — the sweep's hardest
math — cancel identically, so the code path never runs (ERR-026
was invisible to exactly such an MMS). Defense: every
multi-dimensional MMS declares which terms its ansatz **activates**
and which it **nulls**; if the nulled set includes a term with a
catalogued ERR, redesign; always pair an isotropic ansatz with an
angularly non-trivial companion.

**Mode 8 — compiled-out assertion.** A bare ``assert`` under the
canonical ``python -O`` invocation is stripped to a no-op: the test
collects, passes, and asserts nothing — a tripwire that cannot
trip. Defense: per gate, decide the runtime mode explicitly;
sentinel gates use ``np.testing.assert_*`` / ``pytest.fail``
(function calls survive ``-O``) or run without ``-O``.

**Mode 9 — acceleration verified only in a degenerate regime.** A
splitting or accelerator must not change the converged fixed point
— but the invariance is often gated on a regime where the wrong
formulation is *accidentally exact* (a σ_r-fold exact on the
isotropic reflective box yet 46–56 % wrong on vacuum/heterogeneous;
an octant Gauss–Seidel correct on an axis-aligned quadrature yet
converging to a wrong fixed point on a diagonal cubature, ERR-056).
Defense: gate fixed-point invariance on a configuration that
*breaks* the degeneracy (anisotropic flux; diagonal cubature),
asserting equality with the un-accelerated fixed point separately
from any rate claim.

⭐ **And the premise itself can fail.** Every clause above guards a
false GREEN; this one guards a false RED and an ill-posed claim.
"A splitting must not change the converged fixed POINT"
presupposes that a fixed point *exists*. If the operator is
**singular** there is a fixed **manifold**: :math:`\ker A` is
splitting-invariant, but the complementary invariant subspace is
not — so the oblique projector whose range the iteration freezes
differs by splitting, and two perfectly correct splittings
legitimately return **different members**. The gate then reds with
no bug present, and the received wisdom is exactly what makes that
red look like a defect worth chasing. ⚠ The discriminator is *not*
"the bulk did not move": manifold selection and a genuinely
**incoherent** schedule (:math:`M - N \neq A`, the ERR-056 family)
present the same signature from a distance. Three checks separate
them — (a) :math:`M - N \equiv A` bit-exactly for both splittings;
(b) with the kernel **removed**, do the schedules agree on the
*boundary* as well as the bulk; (c) is
:math:`\psi_A - \psi_B` actually in :math:`\ker A`. An incoherent
splitting moves the bulk too. ⟹ **before writing any FP-invariance
gate, ask whether the operator is singular on that gate's own
fixture**; if it is, gate a functional that is
:math:`\perp \ker A` (equivalently: gate the quotient), or pin the
gauge explicitly — never the raw state. The worked ORPHEUS instance
is #344, the SN within-group loss operator on a closed reflective
diamond box (:ref:`sn-loss-kernel-gauge`), where two further traps
travel with it: :math:`\ker A` is **pure-trace**, so every bulk gate
is silent by construction; and the excitation is a **parity** effect,
so a :math:`4, 8, 16, 32` refinement ladder is a single congruence
class and reports "no effect" — a refinement ladder must break the
arithmetic pattern of its own step before any universal is claimed
over it (:ref:`verification-anti-patterns`).

**Mode 10 — activated but unconstrained.** The term's code path
genuinely runs, yet flipping its sign moves nothing: the term
enters the measured quantity as a higher-order-small forcing
absorbed below the convergence floor. "The code path runs" is not
"a sign error is caught" — only a red mutation proves that.
Defense: mutation-verify every activated term that carries a
sign/convention trap; if no isolating regime exists in which the
term dominates, verify the structural pair instead
(machine-precision threading proof + a consumed-sign mutation
moving the answer O(1)) and say so honestly — in the test's scope
note *and* as a prophylactic ``.. warning::`` on the theory page,
so the next reader cannot mint the over-claim the test never
supported.

**Mode 11 — the gate never executes the rewired path.** A refactor
moves logic onto a new production reader, and the closeout credits
an existing green gate as its evidence — but the gate's call graph
never reaches the rewired line (the consumer routes around it). The
green proves the unchanged siblings are unchanged, not that the new
path is correct (the #236 c-fold: a file-write sentinel proved the
credited 640-second twin never called the new reader at all).
Defense: sentinel-instrument the exact rewired line and confirm the
credited gate executes it; then mutation the stamp and confirm the
gate reds.

**Mode 12 — the invariant-functional gate.** The gate measures a
derived functional — an eigenvalue, a balance sum, a normalised
shape — whose **invariance group contains the error class**, so the
error is annihilated algebraically before any assertion sees it:
not sub-floor but *exactly* zero, at every tolerance, in every
regime. The homogeneous :math:`K = A^{-1}F` carve is the canonical
instance: the factor-swap :math:`FA^{-1}` is *similar* to the true
product, so every k-level gate — however tight, however independent
its reference — was designed-green while :math:`\lVert \Delta K
\rVert = O(1)`. Defense, at gate-design time and analytically:
enumerate the measured functional's stabiliser (similarity +
transpose for spectra; per-term cancellation for balance sums;
global scaling for normalised shapes), intersect it with the threat
model, and for any mutation class inside it, gate the constructed
**object itself** against an independently-posed reference — pin
the operator, not just its spectrum. Anti-pattern instances of the
same lens: a 1-group eigenvalue (flux-shape independent —
:ref:`verification-1g-degeneracy`) and "particle balance holds"
(telescoping sums annihilate per-ordinate errors that cancel).
When the invariance is the artefact of a *wrong metric* (a
zero-weight block placing the error inside the stabiliser), the
remedy can be to repair the metric itself so the same functional
becomes a real catcher (ERR-067) — with a control leg proving the
unmutated baseline still passes on the newly-exercised input.

.. _verification-reference-contamination:

Reference contamination — the trusted reference was the problem
---------------------------------------------------------------

The most seductive failure mode in numerical V&V: convergence
checks pass, plots look sensible, two solvers agree — and the
agreement is meaningless because the *reference* was unverified, or
solved a different problem than you thought. Recurring sub-cases:

1. **MC-vs-MC** — two samples from the same generator dressed as
   cross-validation.
2. **Deterministic-vs-unverified-MC** — the solver's error bounded
   by an unknown error; an L4 claim dressed as L1 evidence (the T4
   ban).
3. **Reference converged under the wrong boundary condition** —
   precise, reproducible, and wrong (the recurring Davison
   method-of-images instance, issue #132).
4. **L4 dressed as validation** — "we match the validated code,"
   where the reference code was itself only benchmarked.

The defense is a concrete tracing exercise: *your reference → what
proves it? → what proves that? → … → a structurally-independent
analytical or symbolic ground.* If the chain terminates in "another
code," "the textbook table," or "the previous version of this
solver," you have shown consistency, not correctness — and
consistency is necessary, never sufficient. A related discipline
binds agreement tolerances: a cross-method gate pairing two
verified references must declare its tolerance at the *larger* of
the two truth floors — tighter would mean one method was calibrated
against the other (:doc:`cross_method`).


.. _verification-anti-patterns:

Anti-pattern redirects
======================

The doctrine above compresses into a set of *never → instead*
redirects — the review checklist form of this page. Each names the
seductive move and the principled replacement:

- **Never claim verification from L4 agreement alone** — require an
  L0–L2 chain to a structurally-independent reference; two ORPHEUS
  solvers agreeing is cross-implementation agreement.
- **Never assert closeness to another in-repo solver as a truth
  gate** — match the claim to a reference at the right level and
  pillar.
- **Never accept a 1-group eigenvalue as transport evidence** —
  demand ≥ 2 groups (:ref:`verification-1g-degeneracy`).
- **Never accept homogeneous-only verification** — flat flux nulls
  every redistribution and weight-cancellation term; demand at
  least one heterogeneous, mesh-refined, multi-group case.
- **Never read "convergence rate is correct" as "result is
  correct"** — verify the converged-to value; O(h²) to the wrong
  limit is still O(h²).
- **Never trust an untraced reference** — treat it as contamination
  until the chain to a structurally-independent ground is shown.
- **Never treat "two derivations agree" as proof** — check
  *structural* independence first (ERR-032 agreed to
  :math:`10^{-39}`).
- **Never accept "particle balance holds" as L0 evidence** —
  telescoping sums hold by construction; require per-ordinate
  flat-flux residuals.
- **Never conflate validation with verification** — state which
  screw is being turned.
- **Never accept "it produces reasonable numbers"** — enumerate
  every term, isolate it, verify sign *and* magnitude; sign-flipped
  small terms look reasonable.
- **Never test a contract validator only negatively** — one
  positive and one negative case, with the broken instance built
  independently of the assertion (ERR-051).
- **Never credit a "behaviour-neutral" claim from a fast proxy** —
  re-prove neutrality per consumer with direct value comparison;
  a neutrality proven for one consumer's contract can be false for
  another's (ERR-063).


.. _verification-error-catalog-vocabulary:

The error-catalog contract
==========================

Every bug caught during development is logged as an ``ERR-NNN``
entry in ``.claude/skills/vv-principles/error_catalog.md``, with
its failure mode (1–6), how it hid (which evidence class fooled the
previous tests), which test now catches it, and a one-sentence
lesson. Each entry is pinned by a regression test carrying
``@pytest.mark.catches("ERR-NNN")``, and the :doc:`matrix` flags
any entry without a catcher as a publication blocker.

A ``catches`` marker is a **coverage claim, not a topic tag**: it
asserts that *this* test goes red when the exact documented bug is
re-introduced. A test can live in the same module as the bug and be
structurally blind to it (a matrix-assembly pin that never sees the
dropped inflow factor). Every new marker is mutation-verified —
re-drop the documented bug, confirm this specific test fails under
the canonical ``-O`` invocation; if a different test catches it,
the marker belongs there.

The catalog mixes two evidence classes that the matrix
distinguishes:

- **Code bugs** (most entries) — committed, caught by an L0/L1
  test, fixed, and pinned against regression.
- **Reference-precision floors** (the ERR-038 family) — the
  published reference is itself the bottleneck. ERR-038 (the Atalay
  1997 paper-precision floor, characterised in commit ``4c83e09``)
  is the canonical example: the published numbers sit at a
  1e-3-grade floor, so cross-checks against them must respect that
  floor. This is neither "they agree" nor "the code is wrong" but a
  publication-grade ceiling that *bounds* the claim: a "verified"
  tag against a paper-floor reference is a bounded claim, hardened
  only by a structurally-independent higher-precision reference
  (the canonical path: Atkinson product-Nyström hardening of the
  Atalay values, ERR-036) or a different problem class.

Both classes require an entry; the class determines what "fixing"
means. A future session reading the catalog must attend to the
class before attempting a "fix."


Worked case studies
===================

Six canonical epistemic-failure studies — each answering *what
evidence class would have caught this* — live alongside the error
catalog in ``.claude/skills/vv-principles/scripts/``, ordered by
abstraction: ERR-006 (right order, wrong limit), ERR-032
(procedural ≠ structural independence), issue #100 (falsification
is not localisation), issue #123 (single-quadrature signal is not
closure quality), issue #132 (agreement with the reference is not
agreement with the physics), and issue #226 (agreement in the
eigenvalue is not agreement in the operator — Mode 12). They are
the long-form companions to the failure-mode catalogue above.
