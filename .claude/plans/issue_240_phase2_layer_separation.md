# #240 Phase 2 — separate DiscretizationScheme (agnostic) · SN loss-representation · InvertibleOperator

> **Durable in-repo recovery anchor** (the project rule: plans live in ORPHEUS/.claude/, not ~/.claude).
> Parent: `.claude/plans/next_principled_polymorphism.md`. Branch `feature/sn-space-angle-tier2`.
> Approved 2026-06-15. **STATUS: Step A DONE + committed `d717d4d`. Steps B/C/D NEXT.**

## Context (the user's redirect, 2026-06-15)

A `CellUpdate → DiscretizationScheme` rename surfaced that SN-specific loss-rep realization has
LEAKED into the "discretization" concept. Two investigations (diffusion-mirror on the real
`orpheus/diffusion/solver.py`; an InvertibleOperator Σ-flow trace) confirmed the user's three
diagnoses:

1. **`affine_closure.py` is NOT agnostic.** Its `cell_average`/`source_emission`/(proposed)
   `outgoing_face_from_average` encode the SN sweep's first-order-upwind affine recurrence
   `ψ_out = a·ψ_in + b` + convex face blend. The diffusion solver is a 2nd-order **symmetric**
   stencil (no upwind, no `a`, no `b`, no face-blend cell-average) — consumes NONE of them. These
   are **SN loss-rep realization**; being free functions over one consumer's data they should be
   **methods** of the SN loss representation.
2. **The `Σ_t` in `cell_balance` is a leak.** `InvertibleOperator(L, C)` already puts the Σ-choice
   in `C` (σ_t transport / σ_r removal) and `.solve` honors it (threads `C.sigma` into the sweep,
   `operator.py:1533`). But the **matvec** path reaches for `L.sigma_t` (`loss_representation.py:626,
   1824, 2268`) and feeds `cell_balance` — so `InvertibleOperator(L, C(σ_r))` would sweep with σ_r
   but **matvec with σ_t** (a silent L21 "matvec≡sweep" violation; invisible only because production
   builds L and C from the same σ_t and σ_r has no callers). `cell_balance` itself is principled (Σ
   is a *parameter*); the fix is the caller sourcing it from `C`, not `L`.
3. Therefore the agnostic scheme is THIN: **basis only** (LD `θ`, moment count/projection, DG order,
   traits). Everything Σ/streaming-bearing is SN loss-rep.

## Target architecture (the 3 layers)

- **L4 DiscretizationScheme (agnostic, thin):** basis identity (DD/LD/Step), LD `θ`, moment
  projection (`_ld_source_moments`), traits. NO Σ, NO streaming, NO closure ops.
- **L6 SN loss representation (owns its realization; consumes the operator):** the affine-scan
  closure ops (as METHODS) + the cell-coefficient assembly (`cell_kernel_batch`/`affine_scan_coefficients`)
  + the sweep/matvec walks. Computes the **action of the InvertibleOperator it is given**, reading
  the diagonal Σ from **C** (not `L.sigma_t`).
- **L5 InvertibleOperator:** `L + C`; the Σ-choice lives in `C`; the loss-rep is its
  `.loss_representation`. (Exists; `.solve` honors C — only matvec leaks.)
- **L2 SN method-augmented mesh:** `streaming = |μ|·A_down/V` (raw `g`; the diamond `2` is the
  scheme's). **L1 geometry** (`Mesh` V/A/Δ) already shared (CP+SN).

Mirror-validated role map: agnostic = basis+traits; SN-leaked = `cell_kernel_batch`/
`residual_kernel_batch`/`affine_scan_coefficients`(a,inv)/`update`/`residual`/LD Schur (take Σ/
streaming); SN-sweep-infra = `CellVisit`/`UpstreamState`/`CellResult` (the `StreamingTerms` in
`CellVisit` is the leak vector). `affine_closure.py` = embryonic SN-loss-rep closure (becomes methods).
`cell_balance.py` = the established SN-loss-rep balance home (the relocation destination).

## Staged execution

### Step A — `s_axes = g` (L2 de-pollution) — ✅ DONE + committed `d717d4d`
Producer `geometry.py:1436` → raw `g = |μ|/Δ` (dropped 2.0; accessor+docstrings); DD kernels
`diamond.py` internalise `2·s_a` in **denom AND numer** (closure `2ψ̄−in` unchanged); LD
`linear_discontinuous.py:429` drops the `0.5`; 2-D ScanMarch (`loss_representation.py`
`_sweep_interior`/`_loss_action_interior`) `2·s` consistency with inline DD STAYING (flagged
`#239`); 1-D matvec `:2053` dropped the literal 2.0 (kernel applies it; `abs_mu/V[i]` dedup tracked
on #240). **BIT-IDENTICAL** (only relocates an exact power-of-2; IEEE rounding commutes —
`(2|μ|)/Δ == 2·(|μ|/Δ)`). Strict gate `tests/sn/sweep/core tests/sn/solve -W
"error::tests.sn.regression._regression_assert.DriftWarning"` = 505/1/4, zero drift, NO snapshot
regenerated. (⚠ the test-architect spec's `-W` path `orpheus.sn...` is WRONG — the real module is
`tests.sn.regression._regression_assert`.) Re-baselined convention-encoding unit tests only:
`test_cell_kernel_batch` hand-calc (Design A: param = raw g; denom=Σ_t+Σ2g; hand-derived 18, 120/18)
+ source-of-record sha256 + `_per_ordinate_loop_reference` + `test_sweep_regression` stencil tests +
`test_sweep_graph` legacy-inlined reference + `test_sweep_graph_nd_admission` oracles (denom 4.7, 1.9)
+ `test_linear_discontinuous` 3 LD-kernel sites (`2*mu/h`→`mu/h`, qa-caught). Verified §3 frozen
positive (slab+cart2d) + §4 polymorphic (ScanMarch≡full-field) + §5 references (MMS 2-D/LD-slab/het,
k∞≥2G = 41p); only the 7 pre-existing reds (#195/#209 SPH ×5, #214 mu_y ×2). elegance PASS-WITH-NITS
(NIT-2 = the `abs_mu/V[i]` dedup, tracked on #240; NIT-1 = optional naming in the #239 seam, left),
qa SUPPORTED, Sphinx clean.

### Step B — close the Σ matvec leak — NEXT
`loss_action`/`loss_action_transpose` receive the diagonal from **C** (pass the `InvertibleOperator`,
or thread σ explicitly) instead of `operator.sigma_t`; fix the 3 leak sites (`loss_representation.py:626,
1824, 2268`) + keep the `StreamingOperator` Resolution-A `−σ⊙` subtraction consistent with the matvec
σ. Leaf kernels + `cell_balance` need NO change (parameter-driven). NEW capability gate: build
`InvertibleOperator(L, C(σ_r))` with σ_r ≠ σ_t and assert **matvec ≡ sweep** (true L21) — the test
that would have caught the leak. Behavior-preserving for production (σ_t=σ_t today). Minimal-seam
detail in the explorer Σ-flow memo (this session). `CollisionOperator.sigma` is the diagonal handle;
`InvertibleOperator` exposes `.streaming`(L)/`.diagonal`(C)/`.sigma`.

### Step C — affine-closure free functions → SN loss-rep methods — after B
Move `cell_average`/`source_emission`/`outgoing_face_from_average` (add the missing inverse op
`(ψ̄−(1−w)·in)/w`) from `affine_closure.py` to methods on the SN loss-representation base; remove the
≥4 inlined `2ψ̄−in` (`diamond.py:340,376`, `loss_representation.py` curvilinear matvec, per-cell
`update`). Bit-identical at `w=½`.

### Step D — DiscretizationScheme basis-only + relocate assembly (+rename, +#239) — needs its own design pass
Make the scheme thin (basis+traits); relocate the polymorphic cell-coefficient assembly
(`cell_kernel_batch`/`affine_scan_coefficients`/LD Schur) to the SN loss-rep; **#239** (2-D ScanMarch
onto the coefficient model — needs the NEW generic `outgoing_face_from_average`) folds in; rename
`CellUpdate→DiscretizationScheme` + `CellUpdateBase→DiscretizationSchemeBase` + `cell_update`→`scheme`
attribute. Nexus-rename map: MAIN tree only (the dry-run also lists the sibling worktree
`.claude/worktrees/nexus-workspace-wiring/` — NEVER touch it) and **Base-first** (the `CellUpdate`
substring lives inside `CellUpdateBase`); apply with Edit/sed on the verified site list, NOT
`rename(dry_run=False)`. OPEN design question: home of the *scheme-polymorphic* SN assembly
(generic-parameterized-by-`w`/basis vs an SN-side per-scheme realization hierarchy). Dispatch
test-architect (L17 convention carve) + a design pass FIRST. ⚠ honest caveat: a genuinely BILINEAR
2-D LD (y-slope row coupling, #158 Inc D) does NOT fit "x-scan + transverse-direct-y"; the 2-D
abstraction is correct for the DD/Step/slopeless class and Inc D extends it.

## Verification discipline (every step)
Strict gate `python -O -m pytest tests/sn/sweep/core tests/sn/solve -W "error::tests.sn.regression._regression_assert.DriftWarning"`
(CORRECT path — not the spec's `orpheus.sn`). Route around the 7 pre-existing reds: SPH #195/#209
(`test_vacuum_bulk_bit_identical_1d[*-SPH]` + `test_sphere_{1,2}g_apply_bit_identical`), #214
(`test_2d_mesh_resolution`, `two_d_cartesian_loss_action`). NEVER run all `tests/sn` (#212
`continuous_get` hang). Run `tests/sn/spatial` too (the LD kernel tests — Step A's qa caught the gap
there). Per step: elegance-enforcer + qa on the diff; Sphinx clean; NO bit-identity compromise —
re-baseline if principled.

## Commits (explicit paths only; NEVER the 3 forbidden untracked or the parallel-session set)
One commit per step (feat/refactor) + a `chore(claude)` records commit. Trailer
`Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. NOT pushed/merged.
NEVER stage: `.claude/plans/r1_phase_a_dim_agnostic_ultraplan.md`, `derivations/diagnostics/
diag_s69_scanmarch_vs_window_bench.py`, `scratch/literature/`, the parallel-session set
(`CLAUDE.md`, `.claude/skills/`, `.claude/rules/`, `.claude/lessons.md`, `.claude/hooks/`),
`docs/_build/`. The `.claude/agent-memory/{elegance-enforcer,qa}/` writes ARE this session's review
agents' legitimate memory → chore commit.

## EXIT (when #240 complete — per `next_principled_polymorphism.md`)
Resume `.claude/plans/issue_158_spatial_cellupdate_carve.md` §"Increment C" (#37, LD diffusion
limit); then Step-4 (#36) + Increment D (#38); ff-merge `feature/sn-space-angle-tier2` → `main` +
apply the Phase-0 hand-off to main's CURRENT vv-principles/coding-elegance.
