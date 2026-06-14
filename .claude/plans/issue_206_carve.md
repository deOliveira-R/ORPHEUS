# #206 — Unify the SN matvec/sweep WDD recurrence (single source of denom/a/b) + relocate the 1-D matvec off the operator

## ⭐ START HERE (post-compaction, cold-start)

**✅ PHASES A + B + C ALL DONE (2026-06-14). #206 is CODE-COMPLETE on `refactor/sn-cellupdate-seam-slab`.
Next action = merge ff-only to `main` (carries the #206 close) → then resume Tier 2
(`.claude/plans/sn_space_angle_discretization_plan.md`: #158/#233 spatial ∥ #235 angular).**
See the "⭐ PHASE A + B + C DONE" block below for the full record. The deferred optional nULP
denom-fold (SCAN option) is documented there — NOT part of #206 (the denom is already structurally
single-sourced; both specialists recommended bit-identical-first, the fold's elegance win is marginal
and it forfeits the strict 0-ULP gate).

- **Branch:** `refactor/sn-cellupdate-seam-slab` (off `main @ eab05ab`). **Steps 1–2 of the seam
  (the `cell_update` capability + cache delegation) are DONE + verified but UNCOMMITTED** — they
  live as working-tree modifications (compaction does not touch the tree; they survive). The
  steps-1–2 code files: `orpheus/sn/{geometry,loss_representation,operator,solver}.py`,
  `orpheus/sn/spatial/{cell_update,diamond,sweep_cache}.py`,
  `tests/sn/operators/test_streaming_operator.py`,
  `tests/sn/sweep/core/{test_cell_update_protocol,test_sweep_cache}.py`. (Planning artifacts also
  uncommitted: this file + `issue_236_s2_seam_carve.md` + `sn_space_angle_discretization_plan.md` +
  the agent-memory memos.)
- **What's DONE (don't redo):** `is_affine_scannable` trait + DD's `affine_scan_coefficients` /
  `cell_average_from_faces` / `outgoing_face_from_average`; `CollisionCache.from_geometry(geom,sig_t,
  cell_update)` delegates; `SNMesh.cell_update` retyped `CellUpdate`→`CellUpdateBase` (L18). Bit-identical
  (probe + `test_dd_regression` + sweep/core 437 + cartesian_2d + curvilinear-streaming-equil green).
  Detail = `.claude/plans/issue_236_s2_seam_carve.md` DONE block.
- **Recovery artifacts (READ):** (1) this plan (phases + gates + hazards + DO-NOT); (2) test-architect
  verification spec = `.claude/agent-memory/test-architect/issue_206_cellupdate_seam_verification.md`;
  (3) explorer blast-radius memo = `.claude/agent-memory/explorer/issue_206_matvec_blast_radius.md`
  (its conclusions are also in "Key facts the audit established" below). The §236/space-angle/cluster
  memory all point here.
- **First moves:** (Phase 0) capture the A-NEW heterogeneous-σ_t non-flat single-sweep+matvec baseline
  on truly-clean `eab05ab` (stash the steps-1–2 mods → capture → restore — same stash trick used to
  confirm the pre-existing sphere reds). (Phase A1) route the SLAB closures: matvec `_compute_LpC`
  face-advance + `_run_1d_sweep:2004` → `cell_update.outgoing_face_from_average`/`cell_average_from_faces`;
  gate `CumprodScan`/`ScanMarch.supports` on `is_affine_scannable`. ⚠ line numbers DRIFTED from the
  issue text — re-grep the closure sites; the audit's numbers in this file were current at `eab05ab`.
- **Pre-existing reds (NOT ours — route around, do not fix here):** `test_sphere_{1g,2g}_apply_bit_identical`
  (stale post-ERR-058) + **`TestVacuumMatvecBitIdentity::test_vacuum_bulk_bit_identical_1d[*-SPH]`**
  (3 seeds; the SAME stale-sphere-snapshot class — the committed `bc_extraction_baseline/` sphere `.npz`
  predates the ERR-058 closure-seed fix; PROVEN pre-existing 2026-06-14: `_compute_LpC`/the matvec denom
  helper/the `.npz` are byte-identical `eab05ab`→foundation, and `test_g_adjoint_reciprocity[SPH]` +
  `test_phase_c_crosscheck` sphere-kinf PASS → current sphere matvec is CORRECT, the `.npz` is stale)
  + #212 `test_keff_slab::test_heterogeneous_absolute_keff` hang. See the gate recipe.
- **NEVER `git add`:** `.claude/plans/r1_phase_a_dim_agnostic_ultraplan.md`,
  `derivations/diagnostics/diag_s69_scanmarch_vs_window_bench.py`, `scratch/literature/`, `docs/_build/`.
  Stage explicit paths only. Commit trailer: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.
  Host env → `.venv/bin/python`; canonical test `python -O -m pytest`; bg tmp = `$CLAUDE_JOB_DIR/tmp`.
- **⚠ Steps 1–2 are uncommitted.** Decide with the user whether to commit them as the clean foundation
  (`feat(sn): …`) before Phase A, OR fold them into the Phase-A1 commit. Do not commit unprompted.

## Context

**Why this change.** The 1-D SN within-group solve has two applications of one operator
(L21): the **sweep** `(L+C)⁻¹q` (forward substitution) and the **matvec** `(L+C)ψ`. Today
they are computed by **separate, duplicated machinery**:

- **sweep** = free functions `_run_1d_sweep`/`_sweep_1d_unified` (`loss_representation.py`),
  driven by the two-stratum cache (`affine_scan_coefficients` — the σ_t-epoch `(a, 1/denom)`
  I just landed on `DiamondDifference`);
- **matvec** = `_MSpatialOperatorSum._compute_LpC` ON THE OPERATOR (`operator.py:343`),
  hand-rolling the SAME WDD recurrence via `cell_balance_for_streaming`.

Both are a **shared 1-D scan core** consumed by BOTH `CumprodScan` AND `ScanMarch`'s `s_y=0`
1-D degeneration (`ScanMarch.supports` = `is_1d or 2-D-cartesian`). So the WDD denominator
`denom = 2|μ|·A_down + dA_w·c_out + σ_t·V` is spelled at **three sites** (sweep cache, matvec,
matvec-transpose), and the diamond closure `ψ̄=½(in+out) ⟺ out=2ψ̄−in` is inlined at ~10 sites.
This is the **ERR-026 twin-path debt class**: the moment one site is "improved" without the
others, the sweep and matvec silently diverge on non-flat ψ.

The 2-D wavefront family already fixed its half of this in **S6.3b (`3a79ab3`)**: it moved the
`(L+C)` walk OFF the operator INTO `LossRepresentation.loss_action`, with `apply = loss_action(ψ)
− σ_t·ψ` (Resolution-A `L=(L+C)−C`) as the only operator glue, and a SHARED `_OctantWalk` frame
that forks only at the cell-kernel + emit objects (never a bool flag). **#206 is the unfinished
1-D analogue.** The intended outcome: the WDD recurrence + closure live in ONE place
(`DiamondDifference`, the `cell_update` seam), the 1-D matvec lives in the representation layer
(not the operator), and sweep≡matvec is a code fact.

**Status entering this carve.** The `cell_update` seam capability + cache delegation (the §2
steps 1–2) are DONE + verified, uncommitted on `refactor/sn-cellupdate-seam-slab`. See
`.claude/plans/issue_236_s2_seam_carve.md` (its DONE block). Branch is off `main @ eab05ab`.

**Specialist inputs (read for detail):**
- explorer blast-radius memo (in this session's transcript) — anatomy + consumer map + hazards.
- test-architect spec: `.claude/agent-memory/test-architect/issue_206_cellupdate_seam_verification.md`.

**⚠ Exit:** once #206 lands (all phases merged ff-only to `main`), RESUME the parent thrust at
`.claude/plans/sn_space_angle_discretization_plan.md` **Tier 2** (#158/#233 spatial accuracy ∥
#235 angular). #206 is the FOUNDATIONAL blocker that unblocks an independently-selectable
spatial scheme reaching production through ONE seamed recurrence.

---

## Key facts the audit established (load-bearing)

1. **`_compute_LpC` (operator.py:343) ≡ `_compute_decomposition` (operator.py:772) are
   byte-twins with COPY-PASTED walks** — `_compute_LpC` emits full `(L+C)ψ`;
   `_compute_decomposition` emits the `(M_spatial, M_angular_redist)` split via the SAME walk.
   `apply` (`:1154`) recombines the pair. This is a Pattern-2 duplication *inside* the operator.
2. **Production surface (L20):** `_compute_LpC`/`_compute_LpC_transpose` — exactly **2
   delegations each** (`CumprodScan` + `ScanMarch`-1D `loss_action`/`loss_action_transpose`,
   loss_representation.py:738/752/1362/1428). **No test calls them by path** (privatised behind
   `apply`/`loss_action` in Wave-T T.5) → **the L20 test-migration is docstring-only.**
   `_compute_decomposition` is the entangled one: **3 production consumers** (transient
   orchestrator `:269`, `M_spatial.apply` `:1154`, `M_angular_redist.apply` `:1272`) — do NOT
   casually delete it; the `M_spatial`/`M_angular_redist` split is a finer concern than `(L+C)`.
3. **The matvec IS a scan.** `_compute_LpC`'s per-cell forward loop
   (`psi_face_in = 2·psi_cell − psi_face_in`, `:469`) is a first-order recurrence — the SAME
   shape `ScanMarch._loss_action_interior` reconstructs faces with (α=−1, β=2ψ̄) via
   `scan._x_scan_faces`. So the 1-D matvec can be re-expressed as: consume
   `affine_scan_coefficients` + reconstruct faces via `outgoing_face_from_average`. The
   face-advance `2·psi_cell−psi_face_in` is `outgoing_face_from_average(psi_cell, psi_face_in)`
   un-routed (bit-identical — same algebra).
4. **Bit-identity caveat:** the matvec's density-form `(denom·ψ−numer)/V` vs the sweep's
   `b = 2·QV·inv_denom` scan are algebraically equal but may differ at ULP across the
   V-multiply/÷. Face-closure routing is bit-identical; denom/coefficient routing for the matvec
   is **principled-equivalent** (nULP, pinned by the dual-view contract) unless the matvec is
   kept in its exact density grouping.
5. **Coupling hazards (Phase C):** transient orchestrator (`:269`); the `M_spatial`/
   `M_angular_redist` split (`_compute_decomposition`'s 3 consumers); the curvilinear angular
   adjoint `closure.angular_adjoint` (`_compute_LpC_transpose:752`, pinned by
   `test_g_adjoint_reciprocity`); the Carlson coupled-pole seed (ERR-058, ALREADY single-sourced
   through `pole_angular_closure` — keep routing, don't re-inline); the **degenerate cylinder
   ordinate branch** (`:503`, NO slab/scan analogue — the hard part); `cell_balance_for_streaming`
   is multi-consumer (also `DiamondDifference.residual`) — stop the matvec calling it, don't delete it.

---

## Plan — 4 phases, each an independently-landable PR (bit-identical/principled-equiv + green + ff-merge)

Slab-first WITHIN each phase (curvilinear is ERR-058-adjacent → second). Main-agent writes
inline, turn-by-turn (surgical SN carve — NOT method-implementer). The face/denom seam already
has its capability built; these phases ROUTE consumers through it and relocate the bodies.

### Phase 0 — baseline + gate wiring  (prereq, cheap)
- Capture the **A-NEW** baseline at `eab05ab` BEFORE any change: raw single-sweep
  (`angular_flux`+`scalar_flux`) AND single-matvec on a fixed-seed random Q with
  **heterogeneous σ_t, ≥2G, slab+sphere+cyl** (`assert_regression(kind="direct",
  reduction_depth=nx)`). Heterogeneous+non-flat is mandatory — `dA_w·c_out` is dead on flat ψ.
- Wire the test-architect gate recipe (below). Confirm `TestVacuumMatvecBitIdentity` +
  `--capture-baseline` (ROOT `conftest.py`) is the reusable random-ψ bit-id gate.

### Phase A — single-source the WDD recurrence + closure (the #206 headline: "one denom/a/b")
Route BOTH the matvec (`_compute_LpC`/`_compute_LpC_transpose`) and the sweep (`_run_1d_sweep`)
through the `DiamondDifference` seam. Collapses 3 denom sites → 1, ~10 closure sites → 1.
- **A1 (bit-identical):** the diamond face closure. `_compute_LpC:469`/`:960` +
  `_run_1d_sweep:2004` (slab) → `outgoing_face_from_average`/`cell_average_from_faces`. Gate
  `CumprodScan.supports`/`ScanMarch.supports` on `is_affine_scannable`.
- **A2 (principled-equiv):** the denom/coefficients. `_compute_LpC`'s `cell_balance_for_streaming`
  denom → `affine_scan_coefficients` (`inverse_denom`); the matvec emits `(denom·ψ̄ − numer)/V`
  from `(1/inverse_denom)` + `|μ|·A_total·ψ_in`. `cell_balance_for_streaming` stays (still used by
  `DiamondDifference.residual`); the matvec just stops calling it.
- Order within A: slab → curvilinear (`:2167` + the M-M thread) → ScanMarch-2D (`:1406`/`:1418` +
  `scan.py:331-332`).
- Gate: A-NEW + `test_dd_regression` under `-W error::DriftWarning` (A1 strict) +
  `test_cache_populator_matches_cell_balance_terms` (the dual-view contract, now literal).

### Phase B — extract the shared 1-D-scan frame
Recognize `_run_1d_sweep` (solve) + `_compute_LpC` (apply) as **two directions of ONE 1-D scan
walk** and build a shared frame analogous to `_OctantWalk`: forks ONLY at the cell-kernel
(scan-solve via `ordinate_scan` / scan-apply via the α=−1,β=2ψ̄ face reconstruction) + the emit
policy (`_SweepEmit` angular/moment vs `(L+C)ψ` bulk + boundary). NEVER a bool `is_solve` flag.
Relocate the free helpers (`_run_1d_sweep`/`_sweep_1d_unified`/`_ensure_*_cache`) into it —
shared by CumprodScan.sweep + ScanMarch-1D (NOT folded into CumprodScan). Pure relocation →
bit-identical (whole 1-D suite under `-W error::DriftWarning`).

### Phase C — move the 1-D matvec off the operator (the 1-D half of S6.3b)
- `CumprodScan.loss_action`/`ScanMarch`-1D own `(L+C)ψ` via the shared 1-D frame; retire the
  `operator.M_spatial._compute_LpC` delegations. `apply = loss_action(ψ) − σ_t·ψ` (Resolution-A,
  the only glue — mirror `apply_transpose`).
- Collapse the `_compute_LpC`/`_compute_decomposition` byte-twin: re-express
  `_compute_decomposition` (the M_spatial/M_angular_redist split) by reusing the relocated
  coefficients with the angular term isolated — keep its 3 consumers working.
- PRESERVE: `closure.angular_adjoint` (transpose), the Carlson seed via `pole_angular_closure`,
  the degenerate-cylinder branch, the 2-D `NotImplementedError` geometry dispatch.
- This is the highest-risk phase; do slab → sphere → cylinder; lean on the four-leg gate set.

---

## Verification (test-architect spec — full memo in agent-memory)

**Four-leg standoff (L14), 1-D slab+sphere+cyl:**
- leg 1 (sweep≡independent ref): `test_phase_c_crosscheck.py` (trajectory_resolvent Variant-α
  Green's fn + `kinf_homogeneous`).
- leg 2 (matvec≡ref): `test_bc_extraction_matvec.py::TestStreamingEquilibriumValue` +
  `diag_p42_adjoint_oracle.py` (transpose).
- leg 3 (sweep≡matvec twin): `test_loss_action_convention.py` +
  `test_streaming_operator_decomposition.py` + `test_sweep_vs_apply_consistency.py`.
- leg 4 (refinement): `test_phase_c_crosscheck` flux-shape rows.
- adjoint: `test_g_adjoint_reciprocity.py` + its L11 negative control
  `test_wrong_trace_metric_breaks_reciprocity` (MUST stay green; only -O-firing curvilinear
  `angular_adjoint` gate).
- Mode-9 (L27): heterogeneous σ_t + anisotropic P1 ≥2G + **per-ordinate** (not weight-summed),
  non-flat. `test_affine_carve_bit_identity::si_slab_2g_het` (sha256 anchor) + a sphere aniso companion.

**Gate-run recipe** (route around pre-existing reds — confirmed NOT ours):
⚠ `DriftWarning` MUST use the QUALIFIED path — bare `error::DriftWarning` raises
`AttributeError: module 'builtins' has no attribute 'DriftWarning'` (pytest can't resolve it).
1. Phase A/B strict bit-id (1-D, -O): `pytest -O tests/sn/regression tests/sn/sweep
   tests/sn/solve/test_affine_carve_bit_identity.py
   -W "error::tests.sn.regression._regression_assert.DriftWarning"
   -k "not (sphere_1g_apply or sphere_2g_apply)"`  (now ALSO runs the A-NEW gate
   `tests/sn/sweep/core/test_affine_carve_baseline.py`)
2. A/B+C value (1-D, -O): `pytest -O tests/sn/verification/analytical/test_phase_c_crosscheck.py
   tests/sn/operators/test_g_adjoint_reciprocity.py
   tests/sn/operators/test_loss_action_convention.py
   tests/sn/operators/test_streaming_operator_decomposition.py
   tests/sn/operators/test_bc_extraction_matvec.py
   -k "not (vacuum_bulk_bit_identical_1d and SPH)"`  ⟵ routes around the stale sphere `.npz`
3. Mode-8 bare-assert gates (NO -O):
   `pytest tests/sn/operators/test_streaming_operator.py
   "tests/sn/operators/test_bc_extraction_matvec.py::TestStreamingEquilibriumValue"
   -k "not (sphere_1g_apply or sphere_2g_apply)"`
4. `--deselect tests/sn/eigenvalue/test_keff_slab.py::test_heterogeneous_absolute_keff` (#212 hang).

**Pre-existing reds (NOT this carve; PROVEN pre-existing at `eab05ab` 2026-06-14):**
`test_sphere_{1g,2g}_apply_bit_identical` + `TestVacuumMatvecBitIdentity::test_vacuum_bulk_bit_identical_1d[*-SPH]`
(BOTH stale post-ERR-058 #195 sphere snapshots — current sphere matvec is CORRECT per the passing
`test_g_adjoint_reciprocity[SPH]` + `test_phase_c_crosscheck` references; the committed `.npz` baselines
are the stale party, a #195 re-capture task, NOT this carve) + #212 het-keff hang.
FOUNDATION baseline (recipes, 2026-06-14): recipe 2 = 70 pass / 3 SPH-stale fail / 3 skip; recipe 3 = 65 pass.
Per-phase: qa + elegance-enforcer; `python -m tests._harness.audit` exit 0; Sphinx clean.

---

## Critical files
- `orpheus/sn/operator.py` — `_MSpatialOperatorSum` (`:304`), `_compute_LpC` (`:343`),
  `_compute_LpC_transpose` (`:595`), `_compute_decomposition` (`:772`), `apply`/`apply_transpose`,
  the 3 `_compute_decomposition` consumers (`:269`/`:1154`/`:1272`), `M_spatial`/`M_angular_redist`.
- `orpheus/sn/loss_representation.py` — `CumprodScan`/`ScanMarch` `sweep`+`loss_action`(+transpose)
  (`:704`/`:727`/`:752`/`:1221`/`:1332`/`:1362`/`:1428`); `_OctantWalk` (`:445`, the template);
  `_run_1d_sweep` (`:1810`), `_sweep_1d_unified` (`:1712`), `_ensure_*_cache`.
- `orpheus/sn/spatial/diamond.py` — the seam (`affine_scan_coefficients`,
  `cell_average_from_faces`, `outgoing_face_from_average`); `residual` (`:189`, keeps
  `cell_balance_for_streaming`).
- `orpheus/sn/spatial/scan.py` — `ordinate_scan`, `_x_scan_faces`/`_scanmarch_row` (`:331-332`).
- `orpheus/sn/spatial/cell_balance.py:120` — `cell_balance_for_streaming` (multi-consumer, KEEP).
- New: the shared 1-D-scan frame (Phase B); A-NEW baseline test under `tests/sn/sweep/core/`.

## ⭐ PHASE A + B + C DONE — #206 ready to merge (2026-06-14)

**✅ PHASE C COMPLETE + committed** (`refactor/sn-cellupdate-seam-slab`, 5 commits
`eaafbe1`→`f55b90c`): the 1-D `(L+C)` matvec moved OFF the operator INTO `_OneDimScanWalk`.
- `eaafbe1` forward (`_compute_LpC`→`_OneDimScanWalk.loss_action`, verbatim, bit-identical).
- `7300f3e` transpose (`_compute_LpC_transpose`→`.loss_action_transpose`, verbatim).
- `acbcddc` collapse the `_compute_decomposition` byte-twin: shared `_apply_walk` core
  (emits `m_full` + optionally `m_ang` via `emit_angular`), `loss_action` / `loss_action_decomposed`
  thin assemblers; `_compute_decomposition` → one-line delegation.
- `c1edaf7` retire `_compute_LpC`/`_compute_LpC_transpose` (431 lines) + type/doc cleanup
  (`sigma_t` passed directly — no operator duck-typing; buffer-presence guards; stale-ref fixes).
- `f55b90c` review-finding doc fixes (NIT-1 stale transpose msg + NIT-2 phantom ψ-keyed cache narrative).

**DECISIONS (load-bearing):** (1) **DENSITY** option — relocate VERBATIM (keep `cell_balance_for_streaming`
density grouping) → strict **0-ULP bit-identical**, NOT the nULP scan-cache denom-fold. The denom is
ALREADY structurally single-sourced (`cell_balance_for_streaming` ≡ `affine_scan_coefficients` denom,
byte-for-byte). (2) **`emit_angular` fork is PERF-LOAD-BEARING:** always-emitting `m_ang` was a MEASURED
+14–19% on the matvec hot path (the production Krylov `apply` discards `m_ang`); the hot `loss_action`
(`emit_angular=False`) skips it → re-benchmarked back to pre-C4 cost. (3) the SPLIT stays on the operator
(operator algebra); `_compute_decomposition` sources the full from the relocated walk → ONE spatial walk.
Operator glue `apply = loss_action ± σ_t·ψ` (Resolution-A) UNCHANGED.

**GATES (ALL GREEN):** strict bit-id A-NEW sweep+matvec 9 @ 0 ULP (`-W error::...DriftWarning` in
sweep/core+solve, the genuine strict floor — qa proved it LIVE via a 1-ULP perturbation that FAILED);
value+blast (operators+sweep/core) 928; Mode-8 no-O 108; g-adjoint reciprocity 12 (+L11 control
quantified: true-metric 1.4e-16 vs wrong-metric 0.16); decomposition invariant T4b/T4c; leg-1 kinf
closed-form + decomposition direct-call 20; broad solve 60 + eigenvalue 57. **Reviews:** elegance-enforcer
PASS-WITH-NITS (both NITs fixed in `f55b90c`); qa all-6-claims SUPPORTED (mechanical AST-diff: forward
+ transpose + degenerate-branch byte-identical mod the 2 substitutions). **Sphinx** build exit 0 (no new
warnings; Nexus refreshed); **V&V audit** exit 0.

**FOLLOW-UPS filed:** #237 (degenerate-cylinder matvec/sweep branch has NO test stressor —
`level_symmetric` has no `|μ_x|<eps` ordinate; PRE-EXISTING hole, branch is correct-by-AST-diff),
#238 (do `M_spatial`/`M_angular_redist` separately-applicable leaves still earn their keep now that
`_compute_decomposition` is cold? — likely real consumer is #2 DSA / #200 precond). The optional nULP
denom-fold (fork sweep+apply at a shared cell-kernel via the scan cache) is DEFERRED — a clean future
step, oracle-pinnable, NOT needed for #206.

**PRE-EXISTING REDS (route-around, NOT this carve, all tracked):** #232 (2 2-D Cartesian `mu_y` BC),
#195 (3 stale SPH vacuum/apply `.npz` snapshots), #212 (`test_keff_slab::test_heterogeneous_absolute_keff`
hang). Recipes + `-k`/`--deselect` in the "Verification" section above.

---

### (historical) PHASE A + B DONE — Phase C was NEXT (2026-06-14, post specialist plans)

**Phase A COMPLETE + committed** (`refactor/sn-cellupdate-seam-slab`, 8 commits `3ca8562`→`f61e1b0`):
the diamond face closure is single-sourced through the `cell_update` seam across CumprodScan +
ScanMarch-1D + ScanMarch-2D + matvec. A2-denom FOLDED INTO Phase C (user decision). Curvilinear
sweep closure landed in A1. scan.py carries ZERO inline diamond arithmetic. Gates: A-NEW strict 6,
sweep/core 443, recipe2 70, recipe3 65, 2-D gate 102 — all strict where escalation works; elegance PASS×2.

**Phase B COMPLETE + committed** (`1facf96` refactor + `8faa08c` elegance + `9b1f316` matrix regen):
the 1-D sweep free helpers (`_sweep_1d_unified`/`_run_1d_sweep`/`_ensure_*_cache`) RELOCATED into
`_OneDimScanWalk` (`@dataclass(frozen=True) mesh`, the 1-D analogue of `_OctantWalk`). `_sweep_1d_unified`
→`.sweep`, `_run_1d_sweep`→`._run`, cache ensures→methods; `_initial_guess_values` KEPT module-level.
2 prod callers route through `_OneDimScanWalk(self.mesh).sweep(...)`. SOLVE-only (fork at GEOMETRY
`is_slab`, NOT solve/apply; matvec stays on operator). B2 gate `TestOneDimScanWalkFrame` (frozen
mesh-holder + no-`is_solve` AST tripwire). PURE RELOCATION → bit-identical: A-NEW strict 6,
sweep/core+affine_carve strict 446 (=baseline), regression drift-profile reproduced byte-for-byte
(12 DriftWarnings, identical ULP), dispatch+B2 14; elegance PASS-WITH-NITS (both nits fixed: 3
theory-doc xrefs + 4 blank lines); Sphinx exit 0 no new warnings (Nexus graph REFRESHED to branch);
audit exit 0. ⚠ test-architect FINDING 1 persists: `-W error::DriftWarning` is INERT in
`tests/sn/regression/` (conftest forces `always`) — the strict gate is sweep/core+solve, regression
is tolerance-only.

**⏭ PHASE C (NEXT — the highest-risk phase, the four-leg standoff):** move the 1-D matvec OFF the
operator (`_MSpatialOperatorSum._compute_LpC`/`_compute_LpC_transpose`) INTO `_OneDimScanWalk` as the
APPLY-direction kernel, + fold the A2-denom (`affine_scan_coefficients`). The explorer's "real Phase-C
design question": the matvec is a per-cell recurrence on a SEPARATE object — Phase C must re-express
it as the α=−1,β=2ψ̄ scan (via `_x_scan_faces`) AND unify with `M_spatial`'s spatial sum, then
`apply = loss_action(ψ) − σ_t·ψ` (Resolution-A, the only operator glue). PRESERVE: `_compute_decomposition`'s
3 consumers (M_spatial/M_angular_redist split), `closure.angular_adjoint` (transpose), the Carlson seed
via `pole_angular_closure`, the degenerate-cylinder branch. START Phase C with: test-architect
(four-leg standoff plan) + explorer (apply-kernel unification design) dispatches — proactive triggers
for a high-risk operator-algebra carve crossing the operator↔loss_representation boundary. Likely
principled-equiv (nULP) on the denom, bit-identical on the closure (per "Key facts" #4). Best done
with FRESH context (correctness-critical).

### ⚠ test-architect FINDING 1 — `-W error::DriftWarning` is INERT in `tests/sn/regression/`
`tests/sn/regression/conftest.py:25` forces `always::...DriftWarning`, OVERRIDING the CLI `-W error`.
So "strict bit-id via `-W error`" is a NO-OP for the regression suite. The GENUINE strict gate is
`tests/sn/sweep/core` (no conftest override — the A-NEW gate lives here) + `tests/sn/solve/test_affine_carve_bit_identity.py`.
The regression suite is a TOLERANCE gate (diff the DriftWarning summary, don't expect zero).
**Pre-existing 1-D regression drift profile (Phase B/C MUST reproduce verbatim, NOT zero):**
`cyl_1g_homogeneous_{LS4,product}_dd_n20` (scalar ~272k/298k ULP), `sphere_2g_homogeneous_dd_n20`
(k_eff 3032 / scalar 30580 ULP), `sphere_2g_p1_aniso_dd_n20` (scalar 278 ULP). ALL slab rows strict-clean.

### Phase B verification (test-architect)
- **Strict gate (the REAL one):** `pytest -O tests/sn/sweep/core tests/sn/solve/test_affine_carve_bit_identity.py -W "error::tests.sn.regression._regression_assert.DriftWarning" -p no:cacheprovider`. Foundation-for-B baseline: sweep/core **443 pass / 1 skip / 4 xfail**; +affine_carve+regression **16 pass**; A-NEW **6 strict**; dispatch **12**.
- **Drift-profile reproduce:** `pytest -O tests/sn/regression/test_dd_regression.py -k "not 2d"` → diff DriftWarning summary vs above.
- **B2 import-surface (NEW gate to add):** extend `tests/sn/sweep/core/test_unified_sweep_dispatch.py` (or new `test_one_dim_walk.py` mirroring `test_one_octant_walk.py`): relocated frame resolves from new home; CumprodScan/ScanMarch-1D `.sweep` still bind through the frame; NO `is_solve`/`is_apply` bool param (signature inspection — the anti-degradation tripwire).
- **B3 cumprod≡wavefront:** `tests/sn/sweep/core/test_wavefront_cumprod_equivalence.py` EXISTS + green (CumprodScan≡FullFieldWavefront + analytical k_inf). Re-use.
- No new frame oracle needed (A-NEW covers it). No SPH route-around needed (those reds are in `operators/`, not Phase-B paths).

### Phase B design (explorer blast-radius + synthesis) — `_OneDimScanWalk`
- **Relocation is CLEAN:** `_sweep_1d_unified:1719`, `_ensure_geom_cache:1782`, `_ensure_coll_cache:1791`,
  `_run_1d_sweep:1817` (+`_initial_guess_values:1675`) have ZERO external callers; the 2 prod callers
  (`CumprodScan.sweep:726`, `ScanMarch.sweep:1246`) call `_sweep_1d_unified` IDENTICALLY; tests mention
  names in docstrings only. Pure intra-module relocation → bit-identical, no import/test churn.
- **Template:** `_OctantWalk` (`:444-682`) = `@dataclass(frozen=True)` sole field `mesh: SNMesh`, no `__init__`;
  `_interior_walk` shared seam; forks via injected OBJECTS (cell-kernel callable + `_SweepEmit`), NEVER a bool
  (tripwire `:466-468`, `test_one_octant_walk.py`). Build the 1-D analogue `_OneDimScanWalk(mesh)`.
- **What the 1-D frame hosts:** the two-stratum cache ensure/stash (preserve the `sn_mesh._{geom,coll}_cache`
  contract — solver pre-stashes at `solver.py:902-914`, `_ensure_*` is the lazy fallback; keep build-count=1),
  the slab joint-batch (`:1961-2033`) + curvilinear per-ordinate (`:2040-2197`) bodies, `_initial_guess_values`,
  boundary-face I/O (slab xmin+xmax, curvilinear xmax-only), the moment-mode `ValueError`.
- **Phase-B fork = GEOMETRY (`is_slab`), NOT solve/apply** — the 1-D scan's only variation is geometry.
  Phase B wires SOLVE only; leave the seam where the Phase-C apply-kernel attaches.
- **DO NOT retrofit `_SweepEmit`** — the 1-D body writes ndarrays directly at chain/ordinate granularity
  (≠ `_SweepEmit`'s per-octant granularity); keep direct writes (simplest + bit-identical).
- **Coupling (preserve):** pole_angular_closure/Carlson seed (curvilinear, `:2049`/`:2080-2088`, Pattern-2 with
  matvec — route, don't re-inline), boundary_flux, `initial_guess` passthrough (slab ignores, curvilinear consumes).
- **⚠ Phase-C real design question (explorer gap):** the 1-D matvec is `M_spatial._compute_LpC` (separate object,
  per-cell recurrence), NOT a scan in a shared frame. Phase C must pull the per-ordinate recurrence into the
  frame as the apply-kernel (re-express matvec as the α=−1,β=2ψ̄ scan via `_x_scan_faces`) AND unify with
  `M_spatial`'s spatial sum + fold A2-denom (`affine_scan_coefficients`). The four-leg standoff lives here.
- Memos: `.claude/agent-memory/explorer/issue_206_phase_b_frame_blast_radius.md` (blast radius),
  `.claude/agent-memory/test-architect/issue_206_cellupdate_seam_verification.md` (refined Phase-B gates).

### Phase B implementation note
The `_run_1d_sweep` body is ~380 lines; relocating it into a method re-indents the whole body + churns
`sn_mesh`→`self.mesh`. Best done with FRESH context (a 380-line bit-identity-critical move in a bloated
context is an L16-class error risk). Gate at every step with the strict A-NEW + sweep/core gates.

## DO NOT
- Do NOT fold the shared 1-D core into `CumprodScan` — it's shared with `ScanMarch`-1D (build a
  shared frame, mirroring `_OctantWalk`).
- Do NOT delete `_compute_decomposition` or `cell_balance_for_streaming` (live consumers).
- Do NOT re-inline the Carlson seed or the angular closure (route through `pole_angular_closure`).
- Do NOT let the relocated matvec return `Lψ` (it returns `(L+C)ψ`; `−C` is the operator glue once).

## Post-approval step 0 (durable home)
Mirror this plan to `.claude/plans/issue_206_carve.md` (the repo path the cluster memory +
parent + §2 plan now point to; `~/.claude/` is ephemeral per CLAUDE.md). Then file a one-line
#206 comment: the issue text's locations are stale (`operator.py:904`/`sweep.py:371` gone); the
live duplication is `cell_balance.py:120` (matvec) vs `diamond.py` `affine_scan_coefficients` (sweep).
