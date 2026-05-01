---
name: vv-principles
description: PROACTIVELY use when reviewing claims of correctness, designing verification plans, or evaluating whether evidence supports a claim. Provides the V&V hierarchy (L0–L3 + foundation), the 6 AI failure modes catalogue, the reference hierarchy by structural independence, anti-patterns, and the hierarchical claim taxonomy. Preloaded by qa, test-architect, numerics-investigator, and archivist.
---

# V&V Principles — claim taxonomy, evidence hierarchy, anti-patterns

This skill is the **decision instrument**. The pedagogy lives in
[reference.md](reference.md). Open this file during reviews,
verification-plan design, and bug triage. Open `reference.md` when
you need the philosophy or the worked case studies.

---

## CRITICAL: Anti-patterns to flag immediately

Each line below is a redirect: **NEVER** do X — **instead** do Y. If
you see the left-hand pattern in a PR, claim, or doc, raise it before
any other review work.

1. **NEVER** claim verification on the basis of L4 agreement alone —
   **instead** require an L0–L2 evidence chain pointing at a
   structurally-independent reference. Two ORPHEUS solvers agreeing is
   _cross-implementation agreement_, NOT _correctness evidence_.
2. **NEVER** assert `np.allclose` against another solver in this
   codebase — **instead** match the claim to a reference at the right
   level (analytical for L1, MMS for spatial convergence, MC tally for
   L4 cross-check only after MC itself is verified).
3. **NEVER** accept a 1-group eigenvalue test as evidence of solver
   correctness — **instead** demand ≥2 groups. k = νΣ_f/Σ_a is
   flux-shape independent; 1G is degenerate.
4. **NEVER** accept a homogeneous-only verification — **instead**
   demand at least one heterogeneous, mesh-refined, multi-group case.
   Flat flux nulls every redistribution and weight-cancellation term.
5. **NEVER** read "convergence rate is correct" as "result is correct"
   — **instead** verify the converged-to value. O(h²) to the wrong
   limit is still O(h²).
6. **NEVER** trust a reference that has not been traced back to a
   structurally-independent analytical or symbolic ground —
   **instead** treat it as **reference contamination** until the
   trace is shown. The most seductive failure mode: MC vs MC,
   CP vs unverified MC, method-of-images converged to the wrong BC.
7. **NEVER** treat "two derivations agree" as proof — **instead**
   check whether they are _structurally_ independent. ERR-032 (two
   antiderivatives both using `∫E_2 = 1 − E_3` instead of `½ − E_3`)
   agreed at 1e-39 because they shared the upstream identity, not
   because either was right.
8. **NEVER** accept "particle balance holds" as L0 evidence —
   **instead** require per-ordinate flat-flux residual. Telescoping
   sums hold by construction even with wrong per-ordinate balance.
9. **NEVER** conflate validation with verification — **instead** state
   which screw is being turned. If the equation itself is wrong,
   verification can pass cleanly; only L3 catches it.
10. **NEVER** accept "it produces reasonable numbers" — **instead**
    enumerate every term, isolate it, and verify sign AND magnitude.
    Sign-flipped small terms look reasonable.

---

## The 6 AI failure modes — mechanism and detection

These failure modes are mechanically explainable, NOT arbitrary —
they are the observable signature of sub-word tokenizer co-location.
**L0 verification is the only defense.** See reference.md §2 for the
mechanism (tokenization grounding, AI-targeted but not AI-exclusive).

| #   | Mode             | Example                                | Detection (L0 strategy)                                   |
| --- | ---------------- | -------------------------------------- | --------------------------------------------------------- |
| 1   | Sign flip        | `(a − b)` vs `(b − a)`                 | Heterogeneous eigenvalue diverges under refinement        |
| 2   | Variable swap    | `mu_x` vs `mu_y`; `SigS` vs `SigS^T`   | Per-ordinate flat-flux residual; asymmetric 2G inputs     |
| 3   | Missing factor   | Missing `ΔA/w`, `2π`, volume           | Fixed-source flux spike at r=0 vs `Q/Σ_t` analytic        |
| 4   | Wrong recursion  | `α_{m+1/2}` index drift                | Per-ordinate flat-flux residual                           |
| 5   | Index error      | `face[i]` vs `face[i+1]`               | Non-uniform mesh produces detectably different keff       |
| 6   | Convention drift | Definition site vs usage site disagree | 2G heterogeneous with asymmetric SigS — wrong group ratio |

The catalogued instances live in `error_catalog.md` (ERR-NNN entries).

---

## Hierarchical claim taxonomy — verify the lower layers first

Claims are layered. Each layer adds dependencies. Verify lower layers
before higher ones, and match evidence to the _claim's_ layer.

```
              ┌────────────────────────────────┐
              │  Eigenvalue claim              │  depends on eigenvalue solver
              │  (k_eff, k_inf)                │  + flux shape + discretisation
              └────────────────────────────────┘
                            ↑ depends on
              ┌────────────────────────────────┐
              │  Flux-shape claim              │  depends on the discrete model
              │  (ψ(r,μ,E), φ(r))              │  + boundary conditions
              └────────────────────────────────┘
                            ↑ depends on
              ┌────────────────────────────────┐
              │  Convergence-order claim       │  pure math; lowest dependency
              │  (O(h^p), MMS slope)           │  verifies parts AND whole
              └────────────────────────────────┘
```

Layer reclassifications to apply when reading a claim:

- **Convergence-order results are _math claims_, NOT _solver claims_.**
  They prove the discretisation is consistent — nothing about the
  solved value being correct. MMS lives at this layer.
- **Flux-shape results are _model claims_, NOT _eigenvalue claims_.**
  They depend on the equation and the BC, not on the eigenvalue
  iteration. MMS reaches this layer when the source is structurally
  independent of the code's primitives.
- **Eigenvalue results are _solver claims_.** They bring the iteration
  scheme, normalisation, and convergence test into consideration. MMS
  does NOT directly reach this layer — k-eigenvalue verification needs
  an analytical eigenvalue (homogeneous infinite medium, transfer
  matrix) or a structurally-independent semi-analytical reference.

---

## CRITICAL: The three pillars of verification

Every verification reference is one of three kinds. Each kind proves a
different thing. **NEVER** name a reference vaguely as "analytical" —
**instead** identify which pillar it belongs to, because each pillar
has a different evidence boundary.

### The duality at the centre

Two questions reveal the pillar split:

- **"Given an equation, find the solution"** → **closed-form** analytical solutions
- **"Given a solution, find the equation source"** → **MMS** (Method of Manufactured Solutions)

When neither question closes algebraically:

- **"Reduce the equation to a single integral, evaluate to arbitrary precision"** → **semi-analytical**

Closed-form and MMS are both *analytical* (exact by construction).
Semi-analytical is *exact via arbitrary-precision numerics*. The
distinction matters when judging what claims a pillar can support.

### What each pillar proves

| Pillar          | Convergence-order                  | Flux-shape            | Eigenvalue            | When it applies                                         |
| --------------- | :--------------------------------: | :-------------------: | :-------------------: | ------------------------------------------------------- |
| Closed-form     | ✓ (against exact)                  | ✓ (under assumptions) | ✓ (exact)             | Limited regimes (homogeneous, simple geometry)          |
| **MMS**         | ✓ (great flexibility)              | ✓ (any imposed shape) | **✗** (source-driven) | Any operator that admits a non-vanishing trial solution |
| Semi-analytical | ✓ (against arb-precision integral) | ✓                     | ✓                     | Hard cases with no closed form                          |

**MMS does NOT prove eigenvalues.** This is mechanical, not a
limitation. By construction MMS is a *source-driven* problem — you
imposed the solution, derived the source that makes it true, and the
eigenvalue is whatever k you started with. There is no eigenvalue
information in MMS to verify against. **NEVER** make eigenvalue claims
on the basis of MMS evidence — **instead** match the eigenvalue claim
to a closed-form or semi-analytical reference.

### MMS operational rules

- **Trial solution MUST NOT vanish under derivatives.** Trigonometric
  and exponential functions are the canonical candidates. Polynomials
  vanish at finite derivative order and produce trivial residuals.
- **Trial solution MUST be non-trivial at boundaries** to verify
  boundary-condition handling. A solution that vanishes at the
  boundary by construction tests nothing about the BC.
- **Trial solution MUST stress-test the numerical method, NOT
  minimise source complexity.** Human MMS designs trend toward simple
  sources because hand-derivation of Q^ext is error-prone. AI agents
  using SymPy have no such constraint — the source is derived
  programmatically. **NEVER** pick "the simplest trig that satisfies
  the BCs" when stronger trial functions exist — **instead** pick
  ψ_chosen for stress-test value: high-frequency oscillation, mixed
  scales, near-singular boundary behaviour, non-trivial group-coupling
  for multi-group transport. The simplification heuristic that
  protects humans from arithmetic errors does not serve verification.
  See reference.md §4.3 for the mechanism.
- **Manufactured source MUST be structurally independent of the
  code's primitives.** If the source is generated by the same
  numerical primitives the code uses, MMS becomes a tautology.

### Semi-analytical correctness ladder

Semi-analytical correctness rests on a two-step chain:

1. **Integrator correctness.** For `scipy.integrate`, `mpmath.quad`,
   etc., correctness is commonly assumed (well-tested upstream). For
   custom integrators, integrator correctness is itself a
   verification requirement before this pillar applies.
2. **Reduction correctness.** The reduction from equation to single
   integral is the pillar's load-bearing math. If the reduction is
   wrong, the integral is exact for the wrong equation — a reference
   contamination instance (see anti-patterns).

If both steps hold, the integral evaluation gives the solution to
arbitrary precision. The Peierls reference solver in
`orpheus.derivations.continuous.peierls` is the canonical ORPHEUS
instance.

### Structural independence — applies across all three pillars

Whichever pillar you use, the chain of trust **MUST** terminate in a
structurally-independent ground. **Procedurally-independent ≠
structurally-independent.** Two derivations that use different code
paths but exercise the same integrand or identity are *procedurally*
independent only. When shipping a new reference, force the cross-check
to come from a different *structural* angle — a kernel check (row-sum,
particle balance) AND a closed-form check (eigenvalue, asymptotic
limit) — **NEVER** two derivations of the same closed form.

### Ancillary references — NEVER pillars

These are NOT pillars; they are ancillary uses of references that
already exist:

- **Independent re-derivation** — a different mathematical path to the
  same closed form. Strong cross-check if the paths are *structurally*
  independent (different identity / different integrand). Weak if only
  procedurally independent.
- **Code-to-code (L4)** — Reserve **exclusively** for cross-implementation
  agreement. **NEVER** proves correctness — both implementations could
  be wrong. Every L4 claim **MUST** name its L0–L3 backing.
- **Monte Carlo** — itself a numerical method that needs verification
  (geometry tracker, free-flight sampler, collision physics, tally
  estimators). Useful as a *consumer* of references; **NEVER** a
  *source* of them. Comparing CP-vs-MC is L4 benchmarking, not
  verification, until MC itself has been verified against an
  analytical / probability-chain reference.

---

## V&V level taxonomy — the ladder

```
VERIFICATION — "Are we solving the equations right?"
  L0  Term verification        hand calc vs code, per term
  L1  Equation verification    analytical solutions, MMS, convergence order
  L2  Integration testing      multi-group + heterogeneous, self-convergence

VALIDATION — "Are we solving the right equations?"
  L3  Validation               experimental data (ICSBEP, IRPhE, SINBAD)

INFORMATIONAL — parallel to the ladder
  L4  Benchmarking             code-to-code — produces zero correctness info

ORTHOGONAL TO THE LADDER
  foundation                   software invariants — no theory-page :label:
                               (data structures, factory outputs, algebraic
                               reduction invariants). Use @pytest.mark.foundation;
                               NEVER carry verifies(...).
```

- **L4 is parallel to the correctness ladder, not part of it.**
  L4 produces information about whether two implementations agree —
  it produces zero information about whether either is correct.
  Every L4 claim **MUST** name its L0–L3 backing.
- **L3 is sequenced, not aspirational.** ICSBEP / IRPhE / SINBAD data
  exists; L3 starts after L1 maturity (when the verification matrix
  has populated, verified entries below it). L3 without L2 is
  accidental agreement.
- **Necessity chain.** L1 without L0 = compensating errors. L2 without
  L1 = masked components. L3 without L2 = accidental agreement. L4
  without L0–L2 = proves nothing.

---

## CRITICAL: 1-group degeneracy — canonical statement

**k = νΣ_f / Σ_a is flux-shape independent.** A 1-group eigenvalue
test cannot detect any error in the spatial, angular, or scattering
operators — the result is a material-property ratio, computable
without solving the transport equation. **Multi-group (≥2G) is MUST
for any verification claim.** This statement is the canonical
reference cited by `CLAUDE.md` Cardinal Rule 6, `qa/AGENT.md`, and
`test-architect/AGENT.md`.

---

## CRITICAL: Log every caught bug

This skill owns the L0 error catalog (`error_catalog.md` in this
skill directory). Every agent that loads `vv-principles` is bound by
the following directive.

**MUST** log every bug caught during development → `error_catalog.md`
with:

- **ERR-NNN** (next sequential ID)
- **Failure mode** (1–6 from the AI failure modes table)
- **How it hid** — what evidence-class fooled the previous tests
- **Which test catches it** — linked via `@pytest.mark.catches("ERR-NNN")`
- **Lesson** — one sentence

The catalog is a QA publication artifact and the skill's primary
self-improvement vehicle. **NEVER** close a numerical-bug
investigation without an ERR entry. The catalog grows the skill;
gaps in the catalog mean lessons did not propagate.

---

## Sign-pattern + magnitude fingerprint diagnostic

Sign-pattern + magnitude scaling form a 2-D fingerprint that pins bug
class before debugger steps. The full fingerprint catalog lives in the
adjacent skill — see
[../numerical-bug-signatures/SKILL.md](../numerical-bug-signatures/SKILL.md).
**Read fingerprints before opening mpmath.**

---

## Pointers

- **Catalogued bugs (ERR-NNN):** `error_catalog.md` in this skill
  directory. Every L0-caught bug carries: failure mode (1–6), how it
  hid, which test catches it, lesson.
- **Worked case studies:** `scripts/` in this skill directory, one
  file per epistemic-failure case. See `scripts/_template.md` to
  add a new one.
- **Adjacent skills:** [`numerical-bug-signatures`](../numerical-bug-signatures/SKILL.md)
  (recognition catalog), [`probe-cascade`](../probe-cascade/SKILL.md)
  (factor isolation), [`nexus-verification`](../nexus-verification/SKILL.md)
  (graph-based coverage audit — invoke its tools during a V&V review).

For the philosophy (structural independence, Oberkampf–Roache frame,
tokenization grounding, reference contamination), read
[reference.md](reference.md).
