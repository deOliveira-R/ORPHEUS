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
- **Phase 2 — SCOPE EXPANDED (user redirect 2026-06-15) to a 3-LAYER SEPARATION:
  agnostic DiscretizationScheme (basis only) · SN loss-representation (owns its
  realization, consumes the InvertibleOperator's Σ — not `L.sigma_t`) · the
  InvertibleOperator (Σ-choice lives in C). AUTHORITATIVE detail + staging now in
  `.claude/plans/issue_240_phase2_layer_separation.md` (the §"Phase 2" below is the
  SUPERSEDED narrower s_axes=g+#239 framing — kept for the convention analysis).**
  Staged A→B→C→D: **Step A (s_axes=g) ✅ DONE + committed `d717d4d`** (bit-identical;
  elegance PASS-WITH-NITS, qa SUPPORTED, Sphinx clean). **Step B (close the Σ matvec
  leak) NEXT**, then C (affine-closure free funcs → SN loss-rep methods) + D
  (scheme basis-only + relocate assembly + rename + #239; needs a design pass).
  Two investigations grounded it: diffusion-mirror (`affine_closure` is SN-specific,
  not agnostic) + InvertibleOperator Σ-flow trace (`.solve` honors C, matvec leaks
  `L.sigma_t` — silent L21 break under a removal-form operator).

## Architecture model — the diffusion-mirror factoring (Phase 2 rationale + the deferred target)

Derived by mirroring the SN solver against a *hypothetical* diffusion solver: assemble both
explicit LHS matrices and demand "what's shared" is **coherent from both seats** (the strong
signal). SN matrix row (per ordinate): `[Σ_t·V + streaming-diag]·ψ̄_i − (streaming coupling)·ψ̄_{i-1}`
— lower-triangular (directional, 1st-order → *sweep*). Diffusion row:
`−D̃_{i-½}·φ_{i-1} + [Σ_r·V + D̃_{i-½}+D̃_{i+½}]·φ_i − D̃_{i+½}·φ_{i+1}` — symmetric tridiagonal
(2nd-order → *SPD solve*), `D̃ = D·A/Δ`.

### The six layers (refined 2026-06-15 with the InvertibleOperator ↔ loss-rep separation)
1. **Spatial geometry** (`A_face`, `V`, `Δ`, incidence) — method-agnostic. ALREADY shared:
   `orpheus/geometry/{coord,mesh}.py` (`compute_volumes_1d`/`compute_areas_1d`, `Mesh1D`/`Mesh2D`);
   CP reads `mesh.areas` directly (`cp/solver.py:369`). Smells (low-pri): `SNMesh.volumes/.areas`
   re-export `mesh.*` (obscures origin — read `Mesh` like CP); d≥3 volume stranded in SNMesh
   (`geometry.py:373`, no `Mesh3D`).
2. **Method-augmented mesh** — geometry × the method's per-face coefficient. SN
   `streaming = |μ_n|·A_down/V` (the ONE legit SN augmentation, `_setup_cartesian` `geometry.py:1435`);
   diffusion `leakage = D·A/Δ`. Curvilinear "reduced" geometry (`geometry/reduced_operator.py`,
   shared with MoC/CP): pure part (`face_areas`/`delta_A` from `Mesh`) + M-M augmented
   (`alpha_half`, `redist_dAw = ΔA/w`, `tau_mm`).
3. **Material data** (`Σ_t/Σ_s/νΣ_f`; `D` from `Σ_tr`) — shared xs library (DATA only).
4. **DiscretizationScheme** (DD/LD/Step) — the method-AGNOSTIC closure: flux representation
   (DOFs/cell) + reference-element overlaps (mass `∫P_aP_b`, derivative `∫P_aP_b'`, stiffness
   `∫P_a'P_b'`, face evals `P_a(±1)`). **NO `|μ|`, NO `D`, NO `Σ`.** ⚠ NOT what the code holds today.
5. **InvertibleOperator** — defines WHAT'S ON THE LHS: which operator terms (`L+C−S−F` …) AND which
   `Σ` is on the reaction diagonal (transport `Σ_t` with within-group scatter in the SOURCE, vs
   removal `Σ_r=Σ_t−Σ_s0` folded in). **The `Σ`-choice (transport vs removal) lives HERE — it is an
   OPERATOR concern, NOT a scheme or loss-rep concern** (user, 2026-06-15). Reads Layers 2+3.
6. **Loss representation** — REALIZES the operator's action: assemble the explicit matrix, OR
   forward-substitution (sweep), OR matvec (Krylov subspace). **Agnostic to the operator's contents**
   — it realizes whatever `InvertibleOperator` it is given. To do so it pulls the operator's
   terms+`Σ` (L5), the scheme's closure (L4), and `g` (L2), assembles the cell coefficients, and
   schedules them.

### Coherence result (the strong signal)
Both loss reps need from the scheme the SAME core — flux representation + **mass** overlap (→ the
`Σ·V` reaction diagonal) + **face** map (→ currents) — differing ONLY in which derivative-order
overlap (SN: 1st-derivative = streaming; diffusion: stiffness = leakage), both from the same basis.
And both pull `Σ` from THEIR operator (L5), not the scheme. Shared seam = **the scheme exposes the
basis/overlaps; the InvertibleOperator supplies `Σ`; each loss rep assembles + schedules its own
realization.** Coherent from both seats. ✓

### What's mislayered TODAY (the sharpened point-2 finding)
The current `CellUpdate`/"DiscretizationScheme" CONFLATES Layers 4+5+6: it holds the **loss-rep's
realization work** (assemble the SN cell coefficients) in three schedule-forms, each TAKING `sig_t`
(the OPERATOR's `Σ`, L5) and `|μ|/g` (L2):
- group 1 `update`/`residual` (per-cell, DAG); group 2 `cell_kernel_batch`/`residual_kernel_batch`
  (÷V batched, DAG → `sweep_graph.py:849/890`); group 3 `affine_scan_coefficients` → `(a,
  inverse_denom, w)` (×V scan → `CollisionCache` `sweep_cache.py:465` + the 2-D ScanMarch).
- groups 2 & 3 are **convention-twins** (÷V kernel vs ×V scan) of the SAME relation; the
  method-agnostic closure (DD's `w`, LD's `θ`) is buried + replicated across all three.
→ `affine_scan_coefficients` is **doubly mislocated**: it is the loss-rep's REALIZATION of the
  operator (not a scheme concern), AND it takes the OPERATOR's `Σ` (`sig_t`, L5, not the scheme's).
  The scheme should see ONLY the basis/closure — no `Σ`, no `g`.

### Target vs Phase-2 scope (unify-after-two)
- **TARGET (with #2 consistent-DSA diffusion — the second consumer that validates the seam):**
  scheme = method-agnostic basis/overlaps (L4 only); the cell-coefficient assembly (groups 1/2/3 /
  `affine_scan_coefficients`) becomes the SN **loss representation** realizing its `InvertibleOperator`
  — pulling `Σ` from the operator, `g` from the mesh, the closure from the scheme, assembling
  GENERICALLY (no scheme-branch — the FEM weak-form pattern), and scheduling (matrix/sweep/matvec).
  A future diffusion loss rep reuses the SAME scheme + its own operator. This seam IS consistent DSA
  (diffusion discretization = angular moment of the SN scheme).
- **PHASE 2 (now) — non-speculative cleanup ONLY:** s_axes=g (L2 de-pollution), ScanMarch → scheme
  coefficients (remove inline DD; consistent with the target — the coefficients it consumes are
  realization work that later migrates), `CellUpdate → DiscretizationScheme` rename (correct OBJECT
  name; document it CURRENTLY also carries the loss-rep realization + the operator's `Σ`, both
  migrating with #2). Name `affine_scan_coefficients` accurately + note the DAG twin + that it takes
  the operator's `Σ`. **DO NOT** extract the overlap interface, move assembly to the loss rep, or
  strip `Σ` from the scheme now — that needs the second consumer (DSA) to validate (unify-after-two).
  Record the geometry smells (V/A re-export, Mesh3D gap) low-priority.

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
- **Geometry low-pri** (from the diffusion-mirror map): `SNMesh.volumes/.areas` re-export
  `mesh.*` (diffusion should read `Mesh` directly like CP `cp/solver.py:369`); no `Mesh3D`
  (d≥3 tensor-product volume stranded in `SNMesh.geometry.py:373`). Promote when a shared
  spatial-geometry base is needed (i.e. with #2 diffusion).

## EXIT — when all #240 phases are complete

When Phase 2 lands (s_axes=g + #239 + `DiscretizationScheme` rename + audit; verified +
elegance/qa reviewed + committed) and #240 is closed, resume here:

1. **IMMEDIATE → `.claude/plans/issue_158_spatial_cellupdate_carve.md` §"Increment C" (#37).**
   The LD diffusion limit: thread `φ̂` (flux slope) into the within-group iterate → the
   scattering source `(Q̄, Q̂=Σ_s·φ̂)`; the `test_ld_thick_diffusive_limit_xfail` tripwire
   flips to xpass; #233 pole-cell lift. Then Step-4 close-out (#36: theory-page LD
   coefficient-model section + SymPy derivation into `derivations/`) and Increment D (#38:
   multi-D bilinear LD `cell_kernel_batch` → the 2-D LD two-paths gate). State/roadmap =
   `[[project-issue-158-ld-dag]]`; parent thrust = `.claude/plans/sn_space_angle_discretization_plan.md`
   (Tier 2a).
2. **DEFERRED ARCHITECTURAL → with #2 (consistent DSA):** realize the TARGET in the
   Architecture-model section — scheme = basis/overlaps (L4); the SN cell-coefficient
   assembly → SN loss representation (pulling `Σ` from the `InvertibleOperator`); a diffusion
   loss rep reuses the scheme. The diffusion-mirror exploration (this plan) is the design
   record; it lands when the second method makes the overlap interface concrete.
3. **MERGE:** `feature/sn-space-angle-tier2` → `main` (ff-only). Apply the Phase-0 hand-off
   to main's CURRENT `vv-principles`/`coding-elegance` + soften `.claude/rules/`; coordinate
   with the parallel tool-routing session. Then delete the branch.
