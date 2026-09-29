---
name: coproduct-direct-sum-two-level-frames
description: Verdicts of the 2026-09-28 structural attack on the coproduct prototype vs FullFieldSpace/CoupledSpace — the sum lives at TWO levels (CoproductManifold at points, DirectSumSpace = its L² image), labels are essential iff the summand family has automorphisms, the discrete trace map reads the INFLOW (A_tt not diagonal), the trace block is singular under all-reflective BCs (null space = reflective cycles), Schur reduction direction is a rank guard, the role grading is two typed arrows (ι_! and Restriction.dual()), and the role enums have 0 production readers.
metadata:
  type: project
---

# Coproduct / direct sum — the two-level frame (attack memo
`scratch/posing_sequence/open_items/coproduct/attack_structure/memo.md`, probes beside it;
re-ground every number against the live tree before acting)

**Why:** the user asked for BOTH forms attacked before ruling A2 (one `FunctionSpace` over the
coproduct vs the two production carriers). Recommendation given: ADOPT WITH CHANGES.

**How to apply (the durable verdicts, each with its number):**
- **Two levels, no category error.** `⊔` at the point level is a coproduct ONLY (no projection;
  point count adds: 160/24/180 vs 4096/128/5184 for a product); its L² image is the biproduct
  `⊕` with both `ι*` and `ι_!`. Name each by its perfect match: `CoproductManifold` (0 hits,
  unspent) and `DirectSumSpace` (the tree's word; `CoproductSpace` is a misnomer).
- **Labels are essential iff `Aut(summand family) ≠ 1`.** `[M]` on two ISOMORPHIC summands the
  reversed twin passes 4/4 gates and only the label leg separates; 0 of 3 fixtures are isomorphic
  today; the first consumer is the 6-group precursor system (0 carriers in `kinetics/`). Both
  production digests and both adapter grammars are label-blind (D18: two homes).
- **The discrete outflow map reads the INFLOW.** `[M]` SN slab: `A_tt` of `L+C` has 16 of 32
  off-diagonal `+1` entries (the DD chain's `(−1)^n ψ_in`); diffusion 4 of 8 (`J⁺` rows read `J⁻`).
  So γ is not a map bulk → trace; it is the loss row `[A_tb | A_tt(out←in)]`, and the ruled
  `R_V` has domain `bulk ⊕ trace_in`. The trace-inflow memo's "`L_tt = diag(+I, −I)`" is REFUTED.
- **The trace block is singular under all-reflective BCs**: rank(`A_tt`) = 24 of 32 (n = 8 and
  n = 7), null-space dimension = number of reflective cycles = ordinate pairs × groups = 8
  (predicted by the label digraph, then measured). Full `A` stays full rank. So the ontology's
  REDUCE-the-trace clause needs a rank guard; the bulk-eliminated RESPONSE FORM
  (`A_tt − A_tb A_bb⁻¹ A_bt`, the CP transmission operator, `[R]`) is available 4 of 4.
- **Schur both ways from four `block()` calls**: 1e-15 with red controls ≈ 1; the response form
  densifies the trace block (56 → 232 of 1024).
- **Role grading = two typed arrows.** `ι_!` (primal) and `Restriction.dual()` (cotangent; the
  tree's `LinearOperator.dual()`); the prototype's `role_of(value)` (1 site) is the smell, and
  the `ι_bulk*.H` failure is the SAME defect (the leaf Riesz leg returns an ndarray; the arrow
  knew the role). `DualSpace.of` forgets sum structure in 3 of 3 carriers (control: axes kept)
  → `DirectSumSpace.dual()` = sum of duals is owed.
- **Role enums**: 0 production readers, 0 branches, `OperatorSum.apply = a(x)+b(x)` (the
  `operator.py:200` comment is stale); 93 test sites pin them. The key set of the block grid is
  the structural pattern; a rectangular coupling block's pattern (`labels(cod) × labels(dom)`)
  is unspellable by an endomorphism enum (`[M]` sphere: seed lands on A/bulk only, emission reads
  A/bulk only, `A_BB` FULL with `block_role None` — 3 of 4 blocks). Declarations survive on 2
  monolithic FULL leaves only.
- **Point set lives on the MEASURE**, not the space: 0 of 8 leaves reach a manifold through
  every axis; `DiscreteMeasure.__add__` (`measure.py:718`) is the coproduct measure restricted to
  ONE support — widen it, do not mint.
- **Carrier move**: `CoupledField` gains labels (14 positional sites); `Composite` stays (230
  name reads) as the 2-ary labelled element — respects the provisional ruling.
- **MC seam**: level-1 (manifold) code has 0 finite-basis references; the bank = tagged points =
  `Injection.__call__`'s output.

Related: [[pullback-pair-frobenius-perron-frames]] (coproduct rows unspellable then),
[[coupled-system-field-bc-frames]] (D5: the 2×2 is a free re-association),
[[posing-ontology-clean-attack-frames]] (trace = pullback along the coproduct injection),
[[field-role-typing-faceflux-frames]], [[mc-seam-feynman-kac-particle-frames]].
