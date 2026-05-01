# V&V Principles — reference

This document is the philosophy. [SKILL.md](SKILL.md) is the decision
instrument. Read this when you need to *understand why* the rules in
SKILL.md are what they are, when you are training a new agent on V&V
practice, or when you face an edge case the SKILL.md taxonomy does
not directly cover.

The order of sections is load-bearing. Structural independence
(§1) is the deepest insight; everything else inherits from it.

---

## §1 Structural independence — verification fails when test and code share an error source

The L11 lesson: a verification test proves correctness only if it
exercises a different mathematical path than the code. When test and
code share an upstream identity, an integrand, or a primitive, an
error in that shared element is invisible to the test. Both sides
agree — and both are wrong.

**Concrete: ERR-032 (the ∫E_2 antiderivative bug).** Two "independent"
cross-checks of a Peierls kernel both applied
`∫E_2 = 1 − E_3`. The correct identity is `∫E_2 = ½ − E_3`. The two
implementations agreed to 1e-39 — not because either was right, but
because they shared the upstream identity. The bug surfaced only
when a third reference, derived from a structurally different angle
(particle-balance row-sum), refused to close.

**Definition.** Two derivations are *procedurally independent* if
they use different code paths, different functions, different files.
They are *structurally independent* only if they exercise a different
integrand, a different identity, or a different recursion. The
distinction matters because:

- Procedural independence catches typos, copy-paste corruption, and
  index drift between the two paths.
- Structural independence catches errors in the *content* — wrong
  identity, wrong integrand, wrong constant.

Procedurally-independent ≠ structurally-independent.

**Operational corollary.** When shipping a new analytical reference,
force the cross-check to come from a different structural angle.
Ship in pairs:

- **Kernel check** — row-sum (∑ P_ij = 1), particle balance
  (production = absorption + leakage), positivity bounds.
- **Closed form** — eigenvalue, asymptotic limit, known solution
  in a degenerate parameter regime.

**NEVER** ship two derivations of the same closed form as
"independent verification." That is procedural independence in a
costume. The kernel check and the closed form must be two structurally
distinct angles on the same physics.

---

## §2 Oberkampf–Roache formal frame + tokenization-grounded L0

### The formal frame

Verification proves *code matches mathematics*. Validation proves
*mathematics matches nature*. Two distinct screws.

This distinction is load-bearing in ORPHEUS practice. If the equation
itself is wrong (deeper screwup — a missing term in the discretisation,
a wrong recursion in the angular sweep), verification can pass cleanly.
The code matches the (wrong) mathematics. Validation (L3) catches it
as systematic error that won't close against experiment, no matter
how the verification cases are tightened.

Conflating the two frames produces false confidence. "All my
verification tests pass" does not mean "my equations are right."
"My code agrees with the experiment" does not mean "my code is
verified." Each screw needs its own evidence.

### Why L0 verification at all? The tokenization grounding

L0 verification — every term, hand-checked, isolated, sign and
magnitude verified — exists because the failure structure of
AI-generated numerical code is mechanically distinct from the
failure structure of hand-written code.

Sub-word tokenizers (BPE) chunk text into reusable units. Symbols
that *look* similar to humans — `Σ_a` and `Σ_f`, `μ_x` and `μ_y`,
`E_2` and `E_3` — may share or co-locate tokens, putting them in
adjacent embedding space. The model's probabilistic decoding lands
on a wrong-but-plausible substitution at rates orders of magnitude
higher than human typo rates.

**Therefore:** plausible-substitution errors (sign flip, variable
swap, missing factor, wrong recursion, index error, convention
drift — the 6 modes catalogued in SKILL.md) are NOT random mistakes.
They are the *observable signature* of a mechanically-explainable
failure mode. **L0 verification — every term, hand-checked — is the
only defense.**

### AI-targeted, not AI-exclusive

The same principle protects against any tokenized-error generator:

- **Human typos** — particularly subscript/superscript errors.
- **OCR drift** — `Σ` vs `E`, `1` vs `l` vs `I`.
- **Copy-paste corruption** — index drift when adapting one block of
  code to a sibling block.
- **Convention drift across files** — definition site and usage site
  written months apart, with the convention silently changing.

L0 verification is the canonical defense for all of them. The
tokenization argument explains *why* the defense matters now more
than ever, but the defense itself is older than the failure mode it
mitigates today.

---

## §3 The V&V ladder elaborated

### L0 — Term verification

**Definition.** Every term in the discretised equation, isolated and
verified by hand calculation against the code's output.

**What it proves.** Each operator term computes the value the
mathematics says it should compute, with the right sign and the
right magnitude.

**What evidence counts.** Hand calculation on a problem where every
term except the one under test has been zeroed out (via boundary
condition, material choice, or geometry). Both polarities of any
term that can change sign. Index ordering verified with a
non-uniform profile (uniform profiles hide index errors).

**What it does NOT prove.** That the discretised equation is the
right equation. That the assembly of terms produces the right global
result. That the iteration scheme converges.

### L1 — Equation verification

**Definition.** The full discretised equation matches an analytical
solution, an MMS-derived solution, or a high-precision semi-analytical
reference.

**What it proves.** The assembly of terms — discretised operator,
boundary conditions, source — produces the value the continuous
equation predicts (in the limit, with the expected convergence rate).

**What evidence counts.** Closed-form solutions (homogeneous infinite
medium, slab with constant source, Peierls integral). MMS solutions
where the manufactured source is structurally independent of the
code's primitives. Convergence rate verified against the theoretical
order of the discretisation.

**What it does NOT prove.** That heterogeneous problems are correct.
That multi-group coupling is correct. That the solver handles
realistic configurations.

### L2 — Integration testing

**Definition.** Multi-group + heterogeneous + self-convergence
test cases, with mesh refinement in space and angle.

**What it proves.** The full solver, in a realistic configuration,
converges to a stable answer as resolution increases.

**What evidence counts.** Multi-group eigenvalue convergence on a
heterogeneous geometry (fuel + moderator). Per-ordinate flat-flux
residual on a curvilinear sweep. Self-convergence: the answer at
fine resolution agrees with the answer at coarse resolution within
the expected discretisation error.

**What it does NOT prove.** That the solver matches reality. That
the cross-section data is correct. That the equation-of-state
assumptions hold.

### L3 — Validation

**Definition.** The solver, with realistic cross sections, matches
experimental measurements within stated tolerances.

**What it proves.** The mathematics matches nature for the class of
problems represented by the experiment.

**What evidence counts.** ICSBEP critical-mass benchmarks. IRPhE
reactor physics experiments. SINBAD shielding benchmarks.
Acceptance criteria stated **before** comparison.

**What it does NOT prove.** That the solver is correct outside the
experimental range. That a different cross-section library would
give the same agreement. That a deeper equation error is not being
masked by cross-section tuning.

### L4 — Benchmarking (parallel to the ladder)

L4 is a different *act*: confirmation between codes presupposes both
are independently verified. Two unverified codes agreeing proves
nothing — they may share a textbook error, a numerical convention, or
a copy-pasted lineage. Useful for confidence-building and publication;
**NEVER** as evidence of correctness.

Every L4 claim MUST name its L0–L3 backing. "ORPHEUS-CP agrees with
ORPHEUS-MC" is not a verification statement; "ORPHEUS-CP agrees with
ORPHEUS-MC, both of which are independently verified at L1 against
analytical references" is.

### Foundation tests — orthogonal to the ladder

Before applying the ladder, ask: *is this physics-equation-tied or
software-invariant?*

Foundation tests (`@pytest.mark.foundation`) verify software
invariants — data structures, factory outputs, algebraic-reduction
invariants. They have no `:label:` to point at because there is no
equation in `docs/theory/` they correspond to. They never carry
`verifies(...)`. They sit outside the ladder.

Examples: `test_geometry` volumes summing to the total; `Mesh1D`
frozen-immutability; subdivision producing equal volumes. These are
correctness claims about the software substrate, not about the
physics it implements. Reserve `@pytest.mark.foundation` exclusively
for these cases — promoting a physics test to "foundation" to dodge
the V&V ladder is the failure mode this category is designed to
prevent.

---

## §4 The three pillars elaborated — philosophy and genesis

The pillar structure is not bureaucracy. It emerges from the only
three ways a *reference* solution can be derived independently of the
discretisation under test. SKILL.md states the pillars operationally;
this section argues *why* there are exactly three.

### §4.1 The duality at the centre

A reference must come from somewhere structurally independent of the
discretisation. There are exactly three such places, captured by two
questions plus a fallback:

- **"Given an equation, find the solution"** → closed-form analytical.
  Classical task: governing equation + BCs + assumptions → exact
  expression for ψ(x, μ, E). Eigenvalue falls out of the boundary
  problem; flux shape is the eigenfunction.

- **"Given a solution, find the equation source"** → MMS. Inverse
  task: impose ψ_chosen, substitute into the operator, derive the
  source Q^ext that produces it. You buy flexibility (any ψ_chosen);
  you sell the eigenvalue (problem is source-driven by construction).

- **Fallback: reduce to a single integral, evaluate to arbitrary
  precision** → semi-analytical. When neither question closes
  algebraically but reduction to integral form does. Any well-tested
  numerical integrator gives arbitrary precision.

The duality ("equation→solution" ↔ "solution→equation") is the
conceptual hook that makes the rule structure sticky. Without it,
the constraints feel arbitrary; with it, they are mechanical
consequences of the inversion.

### §4.2 Pillar 1 — Closed-form analytical

The reference is derived from the governing equation along a different
path (typically symbolic in SymPy or hand-derived). If structurally
independent of the code's discretisation — and symbolic ≠ discrete by
construction — it is ground truth in the limit. ORPHEUS examples that
close: infinite homogeneous medium (k_inf = ν Σ_f / Σ_a from the
matrix eigenvalue), bare slab diffusion (sine eigenfunction with
buckling), bare slab transport at flat-flux. Heterogeneous +
curvilinear closed forms are rare.

The derivation itself must be verified. **ERR-032** is the canonical
counter-example: an analytical reference that was wrong because its
derivation shared an upstream identity (`∫E_2 = 1 − E_3`) with the
code's discretisation. See §1 (structural independence) for the
operational corollary.

### §4.3 Pillar 2 — MMS, with the *why* of its rules

MMS is constructive: choose ψ_chosen, plug into the operator, derive
Q^ext, run the code with that source, check that the code reproduces
ψ_chosen under refinement.

**Why "non-vanishing under derivatives" matters.** Polynomials vanish
at finite derivative order. If the operator computes `∂²ψ/∂x²` and
`ψ = x³`, both discrete and continuous evaluations yield `6x`; the
test passes regardless of any error in the second-derivative
evaluation itself. Trigonometric and exponential trial functions
self-reproduce under derivatives (up to factors), so they generate
non-trivial residuals at every order. Operationally: the test only
sees what the trial function makes visible; a trial function that
vanishes at some derivative order makes that order invisible.

**Why "non-trivial at boundaries" matters.** If `ψ_chosen = sin(πx/L)`
on `x ∈ [0, L]`, it vanishes at both boundaries — the BC equation is
trivially satisfied. The code's BC implementation does nothing
non-trivial; the test exercises only the interior operator. To verify
BC handling, ψ_chosen must produce non-zero boundary values or
non-zero boundary flux that the code's BC implementation must enforce.

**Why MMS gives no eigenvalue.** MMS is source-driven by construction.
You imposed ψ_chosen; the source Q^ext you derived produces it; the
eigenvalue is whatever k you started with. There is no external
eigenvalue truth to check against — only the operator's faithfulness
to its own discretisation. The eigenvalue claim **MUST** be matched to
a closed-form or semi-analytical reference, never to MMS.

**Source-independence rule.** Q^ext must be generated symbolically
(SymPy) or by hand. If Q^ext is computed by the same numerical
primitives the code under test uses, the test becomes a tautology.

**Why AI agents inherit a simplification bias — and MUST override it.**
Training corpora are written by humans, who simplify trial functions
because hand-derivation of Q^ext from a complex ψ_chosen is
error-prone. The training signal teaches the model that "good MMS
uses simple trial functions." The simplification heuristic protects
humans from arithmetic errors; it does not serve the verification
objective. AI agents using SymPy derive Q^ext programmatically —
the constraint that justifies simplicity for humans does not bind.

The MMS purpose is to **stress-test** the discretisation. A simple
trial function tests less. Override the inherited bias: pick
ψ_chosen for the operator's failure modes, not for the ease of
deriving its source. Concrete moves:

- Mix scales: `ψ = sin(πx/L) · exp(-α x)` with `α L ≈ 5` exercises
  both gradient and decay regimes simultaneously.
- High-frequency oscillation: `ψ = sin(k π x/L)` with `k ≥ 3`
  forces the discretisation to resolve modes the textbook `k = 1`
  case does not exercise.
- Near-singular boundary behaviour (when BCs admit it): `ψ` with
  `dψ/dx` large near the boundary stresses BC handling beyond what
  smooth trial functions catch.
- Group-coupling for multi-group transport: `ψ_g(x)` with
  group-dependent shapes that activate the scattering matrix
  off-diagonals, NOT a uniform shape across groups.

The simplification bias is **AI-targeted** in the same sense as the
tokenization grounding for L0 (§2): a mechanically-explainable
inherited heuristic that the human constraint justifies but the AI
constraint does not. Override deliberately.

### §4.4 Pillar 3 — Semi-analytical, with the 2-step ladder

The reduction from governing equation to single-integral form is
analytical; the integral evaluation is numerical at arbitrary
precision. The Peierls reference solver in
`orpheus.derivations.continuous.peierls` is the canonical ORPHEUS
instance.

The 2-step correctness ladder:

1. **Integrator correctness.** `scipy.integrate.quad`, `mpmath.quad`
   are well-tested upstream — correctness is commonly assumed.
   Custom integrators (e.g. the Bickley-Naylor `ki_n_mp` in
   `orpheus.derivations.common.kernels`) require their own
   verification before the pillar applies.

2. **Reduction correctness.** The equation→integral reduction is the
   load-bearing math. If the reduction is wrong, the integrator gives
   an exact answer to the wrong question. **ERR-006** is the canonical
   instance: convergence at the right order to the wrong continuous
   limit.

If both steps hold, the pillar's chain of trust is:
integrator (assumed/verified) → reduction (verified at L0/L1) →
arbitrary-precision evaluation.

### §4.5 Ancillary references — NOT pillars

These are uses of existing references, not new evidence sources.

- **Independent re-derivation** is a strong cross-check if the paths
  are *structurally* independent (different identity, different
  integrand); weak if only procedurally independent.
- **Code-to-code (L4)** is *cross-implementation agreement*, NOT
  correctness evidence. Reserve **exclusively** for confidence-building.
  Two codes agreeing share whatever they share — possibly an error.
- **Monte Carlo** is itself a numerical method needing prior
  verification. CP-vs-MC is L4 benchmarking until MC is verified
  against an analytical / probability-chain reference. Once verified,
  MC is a *consumer* of references, **NEVER** a source.

---

## §5 Worked case studies — index

The case studies live as individual files under `scripts/`,
indexable by symptom. Each follows the 4-bullet epistemic-failure
schema (the bug; what evidence existed; why it didn't catch; what
evidence class would have). Add a new case by copying
`scripts/_template.md`. The bar: every entry MUST answer "what
evidence class would have caught this" sharply — entries that
cannot answer that question are catalog material, not skill
material.

The five canonical cases, ordered by epistemic abstraction:

1. **[err006_convergence_to_wrong_limit.md](scripts/err006_convergence_to_wrong_limit.md)** — convergence at the right order is necessary, not sufficient (cylindrical/spherical SN, redistribution + α-recursion).
2. **[err032_shared_algebra_crosscheck.md](scripts/err032_shared_algebra_crosscheck.md)** — procedural independence is not structural independence (slab white-BC ∫E₂ antiderivative).
3. **[issue100_routing_localisation.md](scripts/issue100_routing_localisation.md)** — falsification is not localisation; mode-routing bugs require an internal algebraic invariant (Class-B sphere normalisation).
4. **[issue123_quadrature_crossing.md](scripts/issue123_quadrature_crossing.md)** — single-quadrature signal is not closure quality; require ≥2-quadrature signed-error stability (Direction-C/Q/N falsifications).
5. **[issue132_reference_contamination.md](scripts/issue132_reference_contamination.md)** — agreement with the reference is not agreement with the physics; audit the BC translation step (Davison image series solving the wrong BC).

---

## §6 Reference contamination — named failure mode

**Definition.** Your trusted reference itself was unverified, or
solved a different problem than you thought. The most seductive
failure mode in numerical V&V: convergence checks pass, plots look
sensible, two solvers "agree" — and the agreement is meaningless
because the reference was contamination.

Reference contamination is what §1 (structural independence) is
designed to prevent. Naming it as a failure mode in its own right is
necessary because it appears in distinct sub-cases that each need
their own counter-pattern.

### Sub-cases

1. **MC-vs-MC.** Two Monte Carlo runs (same code or sibling codes)
   agreeing within stated uncertainty. Looks like cross-validation;
   is actually two samples from the same generator. Resolves nothing
   if the geometry tracker, free-flight sampler, or tally estimator
   has a shared bug.

2. **CP-vs-unverified-MC.** A CP solver compared against an MC
   reference that has no L0–L2 backing. The CP error is bounded by
   the (unknown) MC error. An L4 claim dressed as L1 evidence.

3. **Method-of-images converged to the wrong BC.** A reference
   solution computed under a boundary condition that subtly differs
   from the problem being solved. The converged-to value is precise,
   reproducible, and wrong. Issue #132 (Davison) is the recurring
   instance in this codebase.

4. **L4 dressed as validation.** Code-to-code agreement presented as
   "we match the experiment-validated reference," when the reference
   was itself only L4-validated. Two unverified codes agreeing
   produces no validation information.

### How to avoid

Every reference you trust **MUST** be traceable back to a
structurally-independent analytical or symbolic ground (rank 1 in
the reference hierarchy). If the trace breaks — if at any point in
the chain you cannot point at a closed-form derivation that exercises
a different identity than the code under test — the reference is
contamination-prone.

The tracing exercise is concrete:

```
Your reference  →  what proves it?  →  what proves that?  →  ...  →  rank 1 ground?
```

If the chain terminates in "another code" or "the textbook table" or
"the previous version of this solver," you have not yet shown
correctness. You have shown consistency.

Consistency is necessary; it is **NEVER** sufficient.
