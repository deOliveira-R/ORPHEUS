---
name: nexus82-monte-carlo-implementers
description: The MC theory page's equation→code truth (nexus#82) — 22 of 22 DECLARABLE, and the transferable prior that a page of ALGORITHMIC RULES has ~100% implementer density where a page of ALGEBRAIC LAWS has ~50%; plus the measured coverage defect (0 of 4 MC private kernels is imported by any test).
metadata:
  type: project
---

Sibling of [[nexus82-operator-algebra-implementers]] and
[[nexus82-loss-representation-implementers]], same campaign (the licensing-grade
V&V ledger: convert `implements` guesses into declarations). Findings file:
`.claude/inventories/declare_mc.md`. Measured at HEAD `58e46c6f`, 2026-08-18.

**Why:** the ledger cannot adjudicate a coverage claim whose equation has no
declared implementer. The MC page's 22 in-scope equations carry 195 claims.

**How to apply:** read the KIND PRIOR below before starting the next page in the
fanout — it predicts how much of the page will be declarable and where the
judgement calls sit, which is the expensive part of the estimate.

## ⭐ The transferable half — the KIND PRIOR

The three pages audited so far split cleanly by what the page IS:

| page kind | example | declarable |
|---|---|---|
| **algebraic laws / typing rules / canonical forms** | `operator_algebra.rst` | 21 of 40 (19 NONE) |
| **representation laws** (leaf-sum, removal-σ, facewise-separable) | `sn/loss_representation.rst` | mixed; several genuine NONE |
| **algorithmic rules — sample / estimate / wrap / split / branch** | `monte_carlo.rst` | **22 of 22** |

⟹ A chapter that describes an ALGORITHM has ~100 % implementer density: a
sampling rule, an estimator, a population-control rule and a boundary wrap all
have a site by construction. Do not spend the operator-algebra page's
`{identity, law, canonical-form}` triage effort on such a page — spend it on the
two shapes below, which are the only judgement calls that arose.

## The two equation shapes that are genuine judgement calls

1. **An EXPECTATION IDENTITY of a procedure** (`roulette-conservation`,
   `splitting-weight-conservation` — "the proof" subsections). Nothing computes
   `E[w_after]` or `E[N]`; the identity is a *property of* a unique procedure.
   DECLARABLE against that procedure only if the ledger accepts "the unique
   realisation whose expectation this is" as an implementer. ⚠ Check for a guard
   first and say so: `orpheus/mc/solver.py` has **zero** `assert` statements, so
   there is no enforcement site — and under `-O` there could not be one anyway.
   **Rule all such rows the same way; they are structurally identical.**
2. **A TEST-TOLERANCE equation** (`hetero-tolerance`: `|k_MC − k_ref| < 5σ + 0.06 k_ref`).
   This KIND is not in the campaign's `{identity, law, canonical-form, definition}`
   list. Its only computational sites are the two test functions that ALSO claim
   it, so declaring them makes the claim self-adjudicating (a vacuous green).
   Alternative verdict `NOTHING:definition`. Expect one of these per verification-
   heavy chapter.

## Durable facts about the MC subsystem (survive line drift)

- **The whole chapter lives in ONE module**, `orpheus/mc/solver.py`, and in
  **five functions plus three `NeutronBank` methods**: `_precompute_xs`
  (majorant, χ-cumsum, lethargy widths), `_random_walk` (free flight, decompose,
  virtual-collision probability, direction sampling, periodic wrap, branching,
  scattering CDF, fission weight, χ rebirth, collision estimator),
  `_russian_roulette` (all three roulette labels), `_split_heavy` (both splitting
  labels), `solve_monte_carlo` (keff-cycle, keff-mean, sigma-keff, |Δu|), and
  `NeutronBank.{initialize, normalize_weights, save_start_weights}`.
- The only equation implemented OUTSIDE `orpheus/mc/` is `ws-pitch`
  (`StructuredGeometry.wigner_seitz_pin_cell`, plus MC's own
  `ConcentricPinCell.default_pwr`). ⚠ The same equal-area conversion is
  separately LABELLED on the MoC page (`moc-wigner-seitz`, `pitch-recovery`) and
  the CP page (`wigner-seitz`) — four labels, one identity. Any declaration here
  is a page-ownership question, not a code question.
- ⛔ `EnergyGrid.lethargy_widths` is the NEGATION of `mc-lethargy-width-sign`
  (`log(edges[:-1]/edges[1:])`, positive by construction) — same physics,
  opposite sign convention. Never declare it against the MC label.

## Measured defects found while doing this (file as gaps, they are real)

- **`[M]` 0 of 4 MC private kernels is imported by any test.**
  `grep -rn "_random_walk|_russian_roulette|_split_heavy|_precompute_xs" tests/`
  returns **0 hits across the whole test tree**. Four L0 gates instead
  REPLICATE the solver logic inline, saying so in their own comments
  (`tests/mc/test_properties.py:397, 460, 524, 620`;
  `tests/mc/test_gaps.py:640-646`). ⟹ once declarations land, those `verifies`
  claims adjudicate as REFUTED — correctly. Do not soften a declaration to
  green a replicating gate (this is L-013's frozen-RHS family: the test and the
  SUT are two copies, so the gate cannot see the SUT drift).
- **`collision-estimator` is implemented only as the SUM.** The equation's
  `1/(N_act·V)` normalisation exists nowhere; `MCResult.tally` is returned raw.
- **`ws-pitch`'s `p = R√π` direction is duplicated 5× across test modules**
  (`tests/mc/test_monte_carlo.py`, `test_cross_verification.py`,
  `test_convergence.py` ×2, `test_gaps.py`) rather than routed through a helper —
  and that expression IS the ERR-017 fix.
- **All 22 equations had ZERO `implements` edges** before this pass — not even a
  name-token guess, against a corpus of 13 206 such edges. So on this page the
  "a partial answer stands the guess down" hazard does not bite.
