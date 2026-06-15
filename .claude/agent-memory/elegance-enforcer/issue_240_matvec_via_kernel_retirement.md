---
name: issue-240-matvec-via-kernel-retirement
description: #240 Phase-1 PASS — retires the matvec_via_kernel bit-id compromise; uniform Cartesian matvec via polymorphic residual_kernel_batch; assert_regression re-baseline
metadata:
  type: project
---

# #240 Phase 1 — retire `matvec_via_kernel` (commit `001c836`, branch `feature/sn-space-angle-tier2`)

PASS clean. The REVERSAL of the Increment-B (#158) compromise the user flagged as
"never compromise architecture for bit-identity" [[feedback-principled-over-bit-identical]].
This is the textbook "what good looks like" for a compromise-retirement carve.

**What landed:** the `matvec_via_kernel: ClassVar[bool]` capability flag (on `CellUpdateBase`)
+ its `LinearDiscontinuous` `=True` override DELETED. `_OneDimScanWalk._apply_walk`
Cartesian arm (`loss_representation.py:2023-2045`) now calls
`sn_mesh.cell_update.residual_kernel_batch(...)` UNIFORMLY for DD AND LD — zero
isinstance/flag/scheme-name branch. Pure virtual dispatch:
`DiamondDifference.residual_kernel_batch` (`diamond.py:343`, reconstructs diamond march
`psi_out = 2ψ̄ − in_a` at `:376`) vs `LinearDiscontinuous.residual_kernel_batch`
(`linear_discontinuous.py:462`, halves `s` to `g = 0.5*s_axes[0]` at `:429`).

**RULINGS (all PASS):**
1. **DD+LD unification = genuine Pattern-2 SSOT, no flag.** `s_axes = 2|μ|/Δ` is the
   scheme-agnostic shared interface; each scheme owns its closure factor inside its
   own kernel override. Anti-pattern-3 (boolean-flag dispatch) RETIRED. Grep: zero
   live `matvec_via_kernel` refs; 2 surviving mentions are historical "retiring the
   ... special case" docstrings (legit).
2. **`if curvature == "cartesian"` = legitimate GEOMETRY branch, NOT scheme special-case.**
   Curvilinear arm rides `cell_balance_for_streaming` (M-M angular thread = structurally
   NOT a `(a, inverse_denom, w)` coefficient — real math difference). Single predicate
   reused across walk (`:2092`,`:2104`) → dispatcher+guards cannot drift (Smell-#7 strength).
3. **Curvilinear inline `2ψ̄ − ψ_in` (`:2079`) = HONEST single-occupant geometry.**
   `git blame` → predates this commit (`3ea348c` Inc-B). Curvilinear SN provably DD-only:
   `LinearDiscontinuous.affine_scan_coefficients` RAISES on non-neutral curvature
   (`linear_discontinuous.py:538`) + `_require_slab` (`:188`). A per-scheme branch in
   curvilinear would be unrepresentable (one occupant). Face-recon co-located with the
   M-M density it depends on — extracting to share DD's kernel = unify-over-difference
   (Pattern-6 violation). NOT a flag-worthy twin.
4. **`assert_regression(kind="direct")` re-baseline = anti-pattern-#17 done RIGHT.**
   Gate at `nulp=reduction_depth` (principled FP-reduction-depth bound, `_regression_assert.py:162`);
   bit-identity RETAINED as escalatable `DriftWarning` re-pinnable via `-W error::DriftWarning`.
   NOT tolerance-loosening-to-hide-a-bug: structurally-independent refs green (LD MMS O(h²),
   DD analytical k∞), gap documented in docstring+commit. `_assert_arm` consolidation of 3
   slab arms = clean unify-after-two. Scope correct: boundary-slot stays strict
   `assert_array_equal` (`test_bc_extraction_matvec.py:353`, faces unchanged); 2-D Cartesian
   sentinel stays strict (`:365`, rides `loss_action` walk, untouched by 1-D carve).

**Phase-2 deferral CONFIRMED correct (do NOT flag):** `s_axes = 2|μ|/Δ` bakes DD's
`2 = 1/w_DD`, forcing LD's `0.5*s` un-bake. Principled fix `SNMesh.streaming → g`
couples to #239 (2-D ScanMarch consumes same `streaming`) → land together = unify-after-two.

**Gate I RAN (.venv/bin/python):** `TestT4bPreT4RegressionSnapshot` + `TestVacuumMatvecBitIdentity`
→ SLB arms 6/6 PASS; the only reds are 3 **SPH** arms at ~1e15 ULP = the pre-existing
curvilinear structural reds #195/#209 the commit DELIBERATELY keeps flagged (curvilinear
baselines UNTOUCHED; nULP gate hard-fails them rather than masking). Working-as-designed.

**2 NITS (non-blocking):** NIT-1 latent Pattern-6 — diamond-march `2ψ̄−ψ_in` in 2 sites
(`diamond.py:376` + `loss_representation.py:2079`); disjoint geometries, predates commit,
revisit when published LD-curvilinear closure lands (the 2nd occupant). NIT-2 — commit cites
"max ~64 ULP at near-zero cancellation" but gate uses `nulp=nx`; suite green so not firing,
add a docstring line confirming near-zero outliers are within `nx` ULP (or occur at exact-zero
entries where ULP distance is exact).
