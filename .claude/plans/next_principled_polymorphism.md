# #240 — Principled-over-bit-identical + complete polymorphism (post-#158-B)

> **Origin (user, 2026-06-15):** *"I do not favour bit-identical. I favour
> principled changes... You seem to have made compromises to architecture to
> preserve bit identical behaviour, and that was the wrong choice... any
> polymorphism must offer equivalent capabilities... All sweep strategies must
> accept polymorphic discretization without inline behaviour."*
> Durable principle: [[feedback-principled-over-bit-identical]].
> Branch `feature/sn-space-angle-tier2`.

## STATUS (2026-06-15)

- **Phase 0 (favoritism removal): memory DONE; skill-text directive → HAND-OFF
  (below), NOT applied on this branch** (the branch's `vv-principles`/`coding-elegance`
  are STALE — they predate the `.claude/rules/` migration that lives on main; the
  parallel session owns those files; never-stage constraint).
- **Phase 1 (retire `matvec_via_kernel`): ✅ DONE + COMMITTED `fde76ac` (feat)
  + elegance PASS + qa SUPPORTED.** Cartesian matvec routes DD+LD uniformly through
  `residual_kernel_batch` (no scheme flag); DD slab snapshots re-baselined to
  principled-equiv (`assert_regression`, bit-id kept as escalatable `DriftWarning`);
  curvilinear keeps `cell_balance` (M-M, DD-only single-occupant). Verified:
  scan-solve byte-id 505/1/4; operators 7 pre-existing reds only; sweep/spatial 544;
  analytical+MMS+regression 118; LD matvec≡sweep + two-paths green; SPH reds NOT
  masked (hard-fail ~1e15 ULP).
- **Phase 2 (s_axes=g + #239 + audit): NEXT.** test-architect FIRST.

## Phase 0 — favoritism removal: the HAND-OFF (apply on MAIN, not this branch)

The durable principle is captured in [[feedback-principled-over-bit-identical]]
(survives). The skill-text directive below must be applied to **main's CURRENT**
`vv-principles`/`coding-elegance` (NOT this branch's stale copies) by the parallel
session / at the branch's main-merge:

1. **`vv-principles` SKILL.md §"Bit-identity vs principled-equivalence"** — add:
   > *"NEVER compromise the ARCHITECTURE (add a capability flag, special-case a
   > scheme, leave a polymorphic path incomplete, inline scheme behaviour) to
   > preserve bit-identity. When the 3 criteria hold, RE-BASELINE the snapshot —
   > the architecture is the compounding asset; a regenerated ~1-ULP golden is
   > nothing. Bit-identity is welcome ONLY as a free bonus when the implementation
   > is genuinely unchanged; prefer a gate that keeps it as an ESCALATABLE warning
   > (e.g. `DriftWarning` under `-W error`) over an always-on strict byte-pin."*
   (This is the directive whose absence let #158-B add `matvec_via_kernel`.)
2. **`coding-elegance` SKILL.md (Pattern 2 / architecture)** — one-line xref:
   "architecture > bit-identity; never grow a flag/branch/inline-special-case to
   keep old bytes — re-baseline."
3. **`.claude/rules/{coding-standards,vv-testing}.md` (MAIN only — absent here)** —
   soften "bit-identical where possible" / "never relax a tolerance" →
   "principled-equivalent; re-baseline freely; never compromise architecture."
4. Reframe the memories that *celebrate* "bit-identical" (project_issue_206_matvec_carve,
   project_wave_o_operator_algebra) to "principled-equiv is the standard; bit-id was a
   free bonus" when next touched (non-blocking).

(qa note 2026-06-15: the vv-principles structural-independence/masking framework is
already adequate — the gap is specifically the ARCHITECTURE-compromise directive above,
not the masking check.)

## Phase 1 — DONE (record)

Reverted the `matvec_via_kernel` compromise. Commit `fde76ac`:
- `cell_update.py`: DELETED `matvec_via_kernel` ClassVar; fixed `is_affine_scannable` doc.
- `linear_discontinuous.py`: DELETED the override.
- `loss_representation.py` `_apply_walk`: `if curvature=="cartesian"` → uniform
  `residual_kernel_batch` (DD+LD); curvilinear arm unchanged (M-M + inline DD diamond
  march — single-occupant geometry, documented).
- `affine_closure.py`: docstring embodies principled-over-bit-identical.
- `test_streaming_operator.py` + `test_bc_extraction_matvec.py`: slab snapshots
  re-baselined; **bulk** gate → `assert_regression(kind="direct", reduction_depth=nx)`;
  **boundary** gate STRICT `assert_array_equal` (0 ULP — free bonus). SPH/CYL/cart2d
  baselines UNTOUCHED.

Drift measured: bulk mean ~1 ULP, max relΔ ~1e-14 (the 64-ULP max is a ULP-metric
artifact at a near-zero cancellation value); boundary byte-identical.

## Phase 2 — s_axes=g + #239 (2-D ScanMarch) + complete-polymorphism audit

**test-architect FIRST** (convention-crossing operator-algebra carve, L17).

### The s_axes convention finding (from the user's `2|μ|/Δ` question + the audit)

- **Single producer:** `SNMesh.streaming(axis)` returns `2|μ_axis|/Δ_axis`
  (`geometry.py:738`, stores `_streaming_axes`). The DAG walks (`sweep_graph.py`
  walk_full:712 / walk_windowed:757) just SLICE `operands.str_axes` (built from
  `sn_mesh.streaming(a)`, `loss_representation.py:633/2924`). The 1-D `_apply_walk`
  (`loss_representation.py:2038`) DUPLICATES the construction as `2.0*abs_mu/V[i]`
  (a Pattern-2 smell — should consume `streaming` too).
- **Consumers:** DD kernel (`diamond.py:336-338,372-374`) uses `s` directly
  (`denom += s_a`); LD kernel (`linear_discontinuous.py:429`) does `g = 0.5*s_axes[0]`.
- **The `2`** = DD's diamond closure factor `1/w_DD` (provable: `ψ_out = 2ψ̄ − ψ_in`
  substituted into the cell balance gives `denom = Σ_t + 2g`, `g = |μ|A_down/V`). It
  coincides numerically with the geometric `A_total/A_down = 2` (cartesian), but the
  PRINCIPLED primitive is the raw down-face `g = |μ|·A_down/V`.
- **The fix:** `SNMesh.streaming → g` (raw down-face streaming). DD kernel applies its
  diamond `2` internally (`denom += 2*s_a`, `numer += 2*s_a*in_a`); LD drops the `0.5`
  (`g = s_axes[0]`); `_apply_walk:2038` drops the `2.0`. No scheme constant in the
  shared interface or the scheme-agnostic producers.
- **COUPLING to #239:** the 2-D ScanMarch (`loss_representation.py:1318`
  `alpha = 2.0*s_x/D - 1.0`, `:1320` `beta`) ALSO consumes `streaming` with inline DD
  coefficients that assume `s = 2g`. Changing `streaming → g` forces rewriting those
  inline coefficients = #239's body. So s_axes=g and #239 land TOGETHER (user-confirmed
  2026-06-15: "Commit P1 now; s_axes+#239 = P2").
- **Re-baselines** 1-D AND 2-D DD (the d=2 legacy `sigt+sx+sy` fold-order shifts ~1 ULP
  when the `2*s` moves inside the kernel) → principled-equiv re-baseline of the cart2d
  snapshots + the window-equivalence/d=2-bit-id tests.

### #239 — 2-D ScanMarch onto the coefficient model

`scan.py` `_scanmarch_row` (`:341-342` inline `0.5*(in+out)` / `2*avg-in`) +
`loss_representation.py` `_sweep_interior` (`:1318-1320` inline DD `alpha`/`beta`;
`:1419,1434` `2.0*psi_bar_row`) — remove ALL inline DD discretization; route the 2-D
row-march through `affine_scan_coefficients` + `affine_closure` so it consumes ANY
affine scheme. Pin: DD 2-D ≡ its FullFieldWavefront oracle through the genericized path
(the full 2-D LD two-paths gate needs the bilinear 2-D LD kernel — Increment D).

### Complete-polymorphism audit

grep `isinstance(.*DiamondDifference)`, scheme-name branches, `if .*cell_update`
special-cases, inline discretization (`0.5*`, `2.0*psi`, `2*s`) in EVERY sweep strategy
(`CumprodScan`, `ScanMarch`, `FullFieldWavefront`, `MovingFrontierWindow`,
`_OneDimScanWalk`, `_OctantWalk`). Each hit = real structural difference (document) or
compromise (fix). Also: the `_apply_walk:2038` duplicate `streaming` construction
(Pattern-2 — route through `SNMesh.streaming`).

## Phase 3 — resume the LD roadmap

Increment C (#37, diffusion limit: φ̂ in the iterate; xfail tripwire flips green) ·
Step-4 close-out (#36: theory-page LD coefficient-model section + SymPy derivation into
`derivations/`) · Increment D (#38, multi-D bilinear LD `cell_kernel_batch` → unlocks the
2-D LD two-paths gate).

## Verification (every phase)
- `python -O -m pytest tests/sn/sweep/core tests/sn/solve -W "error::...DriftWarning"`
  (re-baselined goldens) + LD suite + `tests/sn/operators` + the sweep/spatial sweep.
  Pre-existing reds to route around: #195/#209 (curvilinear SPH structural — the
  `pre_t4_snapshots` sphere arms + `bc_extraction` `[*-SPH]`), #214 (`mu_y` 2-D),
  #212 (`continuous_get` hang — never run all of `tests/sn`).
- elegance-enforcer + qa per phase. NO bit-identity compromise — re-baseline.

## Tracking
- **#240** (parent) — favoritism removal (P0) + `matvec_via_kernel` revert (P1, DONE) +
  s_axes=g + complete-polymorphism audit (P2).
- **#239** (2-D ScanMarch → coefficient model) — P2, lands WITH s_axes=g (shared producer).

## Follow-up filed/needed
- **Stale curvilinear bit-id snapshots** (the SPH `pre_t4_snapshots` arms +
  `bc_extraction` `[*-SPH]`): strict bit-id baselines red since #206/ERR-058 — need a
  principled re-baseline-to-current OR confirm a real bug (tied to #195/#209). File/track
  as the curvilinear snapshot cleanup (NOT this carve).
