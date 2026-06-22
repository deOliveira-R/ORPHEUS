---
name: issue-158-ld-spatial-occupant
description: #158 Tier-2a Step-1 LinearDiscontinuous cell-update occupant — PASS-WITH-NITS; the (2,ng) moment-source + Schur-scalar-contract rulings Step 2 builds on
metadata:
  type: project
---

# #158 Tier 2a Step 1 — `LinearDiscontinuous` (LD) spatial occupant

First non-DD occupant of the swappable per-cell `CellUpdate` contract
(`orpheus/sn/spatial/`). Slab/Cartesian only, 2-moment (ψ̄ + ψ̂), O(h²),
diffusion-limit-consistent. Reviewed PASS-WITH-NITS on `main` (untracked files).

**Why:** Step 2 wires the per-cell walk to pass LD's source; the two
load-bearing design decisions (A, B below) had to be ratified BEFORE Step 2.

**How to apply:** when reviewing Step 2, the rulings below are SETTLED — do not
re-litigate A or B; check that Step 2 honors the trigger conditions named.

## Decision A — the `(2, ng)` moment-source carrier: ACCEPTABLE staged carrier, KEEP
- NOT shape/stringly-typed dispatch (anti-pattern #3/#4): no runtime shape-branch
  inside one function. The SCHEME is the dispatch (`create("linear_discontinuous")`),
  statically knows its source is `(2,ng)`. Shape = occupant CONTRACT, not discriminator.
  Pattern 5 (scheme is the primitive) working correctly.
- Base contract types `source: np.ndarray` shape-scheme-defined (DD `(ng,)`, LD `(2,ng)`) — extensibility honored, not violated.
- Pattern 6 correctly applied: typed `MomentSource`/slope-`CellResult` from ONE
  consumer = instance-with-ceremony. DEFER to ≥2 multi-moment consumers.
- Illegal state caught loud at boundary (`_ld_source_moments` ValueError names the
  DD-`(ng,)`-drops-slope failure = missing-factor bug class made unspellable).
- **TRIGGER for Step 2:** the `source[0]`/`source[1]` POSITIONAL packing becomes a
  Pattern-7 convention replicated at every PRODUCER once the global moment-contract
  wires a producer. THEN a named carrier (`LDMomentSource(average,slope)` or named
  axis) earns its keep. `_moment_source` test helper is the seed. Record in Step-2
  plan; do NOT build now.

## Decision B — Schur-complement scalar contract: CORRECT, the elegant call, KEEP
- Structural asymmetry exploited: slope ψ̂ is cell-LOCAL (no cross-cell coupling);
  only the face closure ψ_out=ψ̄+ψ̂ crosses the boundary. Eliminating ψ̂ locally
  RECOGNIZES it's not part of the inter-cell contract, not hiding it.
- `residual` returns scalar `(ng,)` defect in ψ̄ = byte-shape-identical to DD's
  `denom·ψ̄−(source+numer)` → NO contract change, no ripple to contract consumers.
- `_schur_terms` = clean Pattern-2 single-source (solve+apply consume it), mirrors
  DD's `cell_balance_terms`/`cell_balance_for_streaming`. Symmetry-in-math→code.
- Slope recoverable from (avg,outflow) pair (ψ̂=ψ_out−ψ̄), no stranded CellResult field.
- Keeping the 2-component system EXPLICIT would be LESS elegant (represents a
  cross-cell-meaningless quantity, invites wrong discontinuous-closure read).

## The 2 nits (approval conditions, neither blocks architecture)
- **NIT-1 (twin arithmetic, do-now):** `update:276-284` re-fetches s_hat (2nd
  `_ld_source_moments` call) + re-spells the slope-moment row by hand, while
  `_schur_terms` already folds that row into `eff_source=s_bar−s_hat·mu·θ/d2p`.
  Slope-row algebra in TWO places. Bug habitat = the DOCUMENTED slope-row SIGN trap
  (LM-1989 §1.4/§6) spelled twice; round-trip CANNOT catch a slope-only drift
  (residual consumes only the scalar ψ̄ path; slope verified by the single-config
  linear-exactness oracle ONLY). Fix = `_schur_terms` returns a named `_LDCellTerms`
  bundle incl. s_hat (also kills the `_,_` discard at `residual:310`); slope
  reconstruction gets ONE home.
- **NIT-2 (doc SSOT, Cardinal Rule 3 close-out):** the regenerated 2×2 + Schur forms
  (`S`,`D₂'`,`eff_source`,`eff_numer_upstream`) live in docstring `:39-74` AND code
  `:251-254` AND (refs) theory page + SymPy script. Published LM-1989 box was
  INTERNALLY INCONSISTENT → this codebase's regen is the SSOT, now in 3-4 places.
  anti-pattern-#11 exception OK for this step (SymPy-pinned) IFF canonical home =
  theory page w/ docstring deferring. Dispatch archivist at close-out; feature not
  DONE until Sphinx thorough.

## Deliberate divergences NOT to flag (recorded so next review doesn't re-raise)
- LD `_schur_terms:240` adds explicit `abs_mu/volume is None` fail-loud guard;
  DD `residual` (`diamond.py:227`) + `cell_balance.py:303` read them BARE (trust the
  Step-2.5 populate-all invariant). LD's guard is the BETTER discipline (-O-safe,
  Pattern 4) and diverges in the SAFE direction → KEEP. The bare-access in DD/
  cell_balance is the latent gap (pre-existing, out of scope; only flag if touching DD).
- `is_affine_scannable=False` redundant w/ ABC default but appropriate explicit
  restatement of a documented trait → not a smell.

## Test gate = exemplary (no changes)
Linear-exactness oracle = structurally-independent (physical property, not algebra
re-statement) → catches the slope-row sign trap the round-trip can't. Guards both
illegal states (curvilinear refuse + DD-shape refuse). Residual mode#2 (linearity)
+ mode#3 (∂r/∂source=−1). n_groups{1,2,4} = 1G control + 2G/4G real gate (corner-probe).

## Increment A (group-2 batched kernel + solver guard) — PASS-WITH-NITS, reviewed separately
Group-1 NIT-1 CONFIRMED FIXED (the `_LDCellTerms` bundle exists; slope has ONE home).
Increment A = the DAG-family group-2 (`_kernel_terms`÷V → `cell_kernel_batch`/
`residual_kernel_batch`) + the solver scan-cache guard + the `cell_update` kwarg.

- **RULING (a) group1(×V)≡group2(÷V) SEPARATE = ACCEPTABLE Pattern-2 stance, PASS.**
  NOT a forbidden twin: they live in different CONVENTION LAYERS (per-cell contract
  ×V w/ ∂r/∂source=−1 vs batch d-generic Cartesian `s=2|μ|/Δ` ÷V). Unifying = unify-
  over-the-difference (the ×V/÷V split IS what each layer expresses). SAME shape as
  DD's per-cell `residual`(via `cell_balance_for_streaming`) vs `*_kernel_batch`
  (inlined fold) split. Structural anchor = `test_group1_equals_group2_flat` machine-
  precision pin (the aggressive-retirement-exception shape: keep both spellings IFF a
  foundation equiv test pins them). LD's group-2 is TIGHTER than DD: solve+apply share
  `_kernel_terms`; DD's batch pair inlines the fold TWICE (no shared helper). DD carries
  3 spellings (scalar helper + 2 inlined folds); LD group-2 = 1.
- **RULING (b) solver guard `reduced is not None and cell_update.is_affine_scannable`
  (solver.py:907) = PRINCIPLED, no scheme-knowledge leak, PASS.** Reads a CAPABILITY
  TRAIT off the polymorphic scheme, NOT identity (no `isinstance(DiamondDifference)`,
  no `.key==`). SAME `is_affine_scannable` trait the scan-family `supports` gates read
  (CumprodScan:703/ScanMarch:1220) → predicate single-sourced as a CONCEPT across
  cache-build + selection, cannot drift (strength case, not the second-spelling smell).
  Cache built from `affine_scan_coefficients` (ABC default RAISES for non-affine) → guard
  makes "scan cache for a scheme w/ no scan recurrence" unrepresentable (Pattern 4).
- **VERIFIED the dispatch end-to-end:** `SNMesh.cell_update` ALWAYS populated (DD default
  geometry.py:271-273 → guard can't AttributeError). `default_for` (loss_repr:1468) routes
  non-affine 1-D past all 3 scan-gated reps → `FullFieldWavefront` (supports==is_cartesian
  :995, the d-generic spine) which reads `cell_update=self.mesh.cell_update` + calls the
  kernel polymorphically (:892/:954). LD reaches prod w/ ZERO new strategy. The rejected
  `OneDimPerCellWalk` twin correctly ABSENT (it'd be the s_y=0 split-over-no-difference).
- **NIT-1 (do-now, CONCERN):** comment `linear_discontinuous.py:360-363` "the two are the
  SAME LD, pinned consistent by the group1≡group2 gate" OVERSTATES the pin. Test
  (`test_group1_equals_group2_flat`:285) pins ONE config, `Q̂=0` flat ONLY. group-1
  `eff_source` folds slope-source `s_hat`; group-2 `_kernel_terms.rhs` has NO Q̂ term →
  they coincide only on the s_hat=0 slice = exactly what's tested. The UNTESTED term is the
  slope-source path, which carries the DOCUMENTED LM-1989 §1.4/§6 slope-row SIGN TRAP — the
  highest-risk term, skipped by the pin. Bug habitat: Inc-C wires Q̂≠0, maintainer trusts a
  "general equivalence" never tested for the slope term. FIX = scope the comment to the
  Q̂=0 slice + name slope-Q̂ equiv as Inc-C-deferred; optional parametrize the pin over ≥2
  (mu,h) (one-point→functional). anti-pattern-#17 / Smell-#7 (comment asserts a property the
  code doesn't establish).
- **NIT-2 (record only):** group-2 positional `psi_in[0]`/`s_axes[0]` = same Pattern-7
  producer-convention trigger as the (2,ng) source — earns a named per-axis carrier when
  Inc-D wires bilinear multi-D LD. Record in Inc-D plan; do NOT build now.
- **ARCH watch (no issue):** DD's `*_kernel_batch` pair is now the LESS elegant sibling
  (inlines fold twice, no `_kernel_terms`). Do NOT change in this PR (fold-order = d=2
  bit-id contract, out of scope). Single-source IF a future edit must touch both folds.
- **Do-NOT-flag (staging, ratified):** d≥2 raises (Inc-D bilinear); Q̂=0 flat (Inc-C);
  curvilinear raises (group-1 guard). FullFieldWavefront docstring "d=1 slab + d=2" now
  also = LD's d=1 non-affine spine — doc imprecision, out of changed-files scope.
