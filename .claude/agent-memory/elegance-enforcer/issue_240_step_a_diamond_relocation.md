---
name: issue-240-step-a-diamond-relocation
description: #240 Step A — relocate DD's diamond closure factor 2=1/w_DD from the producer (SNMesh.streaming) into the scheme kernels; streaming now emits RAW g=|μ|/Δ. PASS-WITH-NITS.
metadata:
  type: project
---

# #240 Step A — diamond-2 producer→scheme relocation (PASS-WITH-NITS)

Branch `feature/sn-space-angle-tier2`, working tree (UNCOMMITTED at review). The producer
`SNMesh.streaming(axis)` / `_setup_cartesian` dropped the `2.0` and now emits the RAW
down-face streaming `g = |μ|/Δ` (was `2|μ|/Δ`); each scheme applies its OWN closure factor.
This is the architectural correction the user's `2|μ|/Δ` question pointed at: the `2` is
DD's diamond `1/w_DD` (from `ψ_out=2ψ̄−ψ_in`), NOT geometry, so it was polluting the
scheme-agnostic accessor (LD had to un-bake `0.5*s`; 2-D ScanMarch inlined it).

## RULINGS (settled — the relocation is the RIGHT move; do not re-litigate)
- **Q1 producer→scheme relocation = CORRECT.** The streaming accessor now reads as
  scheme-agnostic geometry (`g = |μ|·A_down/V`, a real named unit-bearing quantity, Pattern
  3). SSOT verified at producer: grep shows ZERO remaining `2.0*np.abs(axis_cosines)` /
  pre-scaled streaming producers; every prod consumer reads `sn_mesh.streaming(a)` and
  applies its own factor (loss_representation.py:633, :2944 feed str_axes → the DAG walks).
- **DD kernels apply the 2 UNIFORMLY in the left fold** (`diamond.py` `cell_kernel_batch`
  :339-343 + `residual_kernel_batch` :379-383): `denom += 2.0*s_a` AND `numer += 2.0*s_a*in_a`.
  Both gain the 2 → the "denom-only-2 = non-uniform-2× bug → wrong ψ̄" hazard is explicitly
  documented in the docstring (Q4 satisfied). The closure `2ψ̄−in` is a SEPARATE 2 (diamond
  MEAN), correctly distinguished in prose.
- **Both DAG walks single-source the 2 through the kernel.** `walk_full` (:709) +
  `walk_windowed` (:751) pass raw `str_axes_octant` straight into `level_op.cell` →
  `cell_kernel_batch`/`residual_kernel_batch`. The `2` lives in exactly ONE place per scheme
  (the kernel fold). No double-2 / missing-2.
- **LD `g = s_axes[0]`** (linear_discontinuous.py:429, was `0.5*s_axes[0]`) — LD never wanted
  DD's factor; now reads the raw g directly. Correct.
- **Q2 the 2-D ScanMarch `# NOTE(#239)` inline-DD-stays = ACCEPTABLE DOCUMENTED DEFERRAL,
  NOT a twin-path violation.** The 2-D row-march (`loss_representation.py` ScanMarch
  `_sweep_interior` :1311 + `_loss_action_interior` :1430; `scan.py:_scanmarch_row`) does NOT
  route through `cell_kernel_batch` — it is a genuinely different algorithm (DD-only ÷V
  prefix-march, no coefficient cache; the apply/solve _CellResidual/_CellSolve asymmetry).
  This is the same "select narrow / specialise on measured cost" ruling as the standing
  Inc-B 2-D-inline seam. DD is the ONLY 2-D consumer today (bilinear 2-D LD deferred to #158
  Inc D). Issue **#239 EXISTS + OPEN** ("2-D ScanMarch: lift to the group-3 coefficient
  model") → the removal trigger is tracked (anti-pattern #11 satisfied). NOT premature
  unification (Pattern 6). The diamond mean `out_y=2ψ̄−ψ_y_in` (:1448) correctly LEFT (it's
  the cell-average recon, not the streaming factor — prose now says so).
- **Q5 test re-baselines HONEST hand calcs.** Design-A hand-calc: denom 18 = 2+2·3+2·5,
  numer 120 = 16+2·3·4+2·5·8, ψ̄=120/18 (verified by hand). nd-admission: 4.7 =
  0.5+2·(0.3+0.7+1.1), 1.9 = 0.5+2·0.7 (both verified). The sweep_regression cross-check
  (test_streaming_stencil_reproduces_dd_denominator) keeps the OLD hardcoded `2|μ|/Δ` formula
  as the structurally-INDEPENDENT oracle and asserts `2*streaming(0)==old` — structural
  independence preserved (the old formula is hand-written, not transcribed from the new
  accessor). Per-ordinate-loop reference + stencil-value tests + the kernel sha256 goldens
  all correctly carry the new convention. sha256 EXPECTED dict re-hashed with a #240 comment
  explaining why the fold stays bit-identical (2·g == former 2|μ|/Δ, power-of-2 commute).
- **Bit-identity is real**: `(2|μ|)/Δ == 2·(|μ|/Δ)` bit-for-bit (IEEE power-of-2 scaling
  commutes with rounding); the change only relocates an exact ×2/×0.5. Hence the 505/1/4
  green under DriftWarning-error with NO re-baseline of the regression goldens (only the
  hand-calc + kernel-sha goldens move, because their INPUTS changed convention).

## NITS (non-blocking — Q3 naming)
- **NIT-1 (do-now cosmetic, Smell-#7 naming):** the `sx2 = 2.0*s_x` / `sy2 = 2.0*s_y` named
  intermediates in BOTH 2-D arms (loss_representation.py:1317-1318, :1432-1433) are clear
  enough WITH the inline `# DD x denom-term 2g_x` comment, BUT the follow-on line
  `alpha = 2.0 * sx2 / D_row - 1.0` (:1330) reads `2.0 * (2.0*s_x)` = "diamond-closure-2 ×
  denom-term-2g" — TWO distinct 2s multiplied, the trailing comment says so but the
  expression alone is opaque. Optional: name the closure-2 product (`two_x_recon = 2.0*sx2`)
  or keep the comment. Zero correctness impact; the value is right (verified bit-id). Q4's
  "both gain the 2" invariant IS clearly expressed in the diamond.py docstring; the 2-D arms
  rely on the comment rather than a structural guarantee (acceptable — they're the deferred
  #239 seam).
- **NIT-2 (record, the 1-D matvec Pattern-2 dup — ALREADY FLAGGED in-code):**
  `loss_representation.py:2058` `s_axes=((abs_mu / V[i])[:, None, None],)` RE-DERIVES the raw
  g that `streaming(0)` already produces. The `# NOTE(#240)` at :2050-2053 honestly flags it
  as a deferred single-source follow-up ("pending a widths-vs-volumes bit-id check"). This is
  the standing Inc-B-era `abs_mu/V[i]` dup, now correctly emitting raw g (dropped the stray
  `2.0*`). Acceptable-for-now: it is a SECOND spelling of one quantity, but the producer is
  the SSOT and the dup is tracked. The live hazard (the institutional twin-matvec pattern): a
  future edit to the `streaming` convention (e.g. A_down/V for non-tensor-product geometry)
  lands on the accessor but NOT this re-derivation. The note names the collapse destination
  (`streaming(0)[global_dir, i]`). Should be a filed follow-up if not already under #239/#240.

## Idiom catalogue
- The diamond `2` now joins `is_affine_scannable` / `matvec_via_kernel` as scheme-owned
  closure: the PRODUCER emits raw geometry (`g`), the SCHEME owns its closure factor in its
  kernel fold. Future schemes apply their own w. This is the clean end-state the Inc-B
  coefficient-model carve was building toward.
- The 2-D ScanMarch remains the ONE bounded DD-only seam where the `2` is inline (not
  kernel-routed); #239 is its tracked collapse.

## Gates (per dispatch brief — NOT re-run this review)
- `tests/sn/sweep/core tests/sn/solve` = 505 passed/1 skip/4 xfail under DriftWarning-error,
  NO drift. operators+cartesian_2d = 532 passed, 7 pre-existing reds (#195/#209 SPH, #214
  mu_y). MMS 2-D/LD-slab/het + k∞ ≥2G = 41 passed.
