# #240 Phase 2 Step D — home every method into its native base (DiscretizationScheme vs LossRepresentation)

> **Durable in-repo recovery anchor** (project rule: plans live in ORPHEUS/.claude/, not ~/.claude).
> Parent: `.claude/plans/next_principled_polymorphism.md` / `issue_240_phase2_layer_separation.md`.
> Branch `feature/sn-space-angle-tier2`. Approved 2026-06-16.
> **STATUS: D1–D4 DONE + committed. D5 (RE-SCOPED — a CAMPAIGN: #239 N-D DD scan-march +
> N-dim LD / complete polymorphism, subsumes #158 Inc D) + D6 (docs) NEXT. Pick up at D5's
> design pass (proactive test-architect FIRST). See the D5 section + [[project-issue-158-ld-dag]].**
> Commits: A+B `f0d68c3`/`4937c3a` · D1 `8bc1a49` · D2 `784edeb` · D3 `4f04126` · D4 `c40a341`
> (+ chore records each). D4 finding: scheme was ALREADY Σ-stateless → D4 = the
> diffusion-readiness contract gate + the Base interface note (no code change). The deferred
> model-agnostic param rename (cross-domain-attacker) is filed as **#241**
> (`total_xs`→`reaction_xs`, `cell_average_weight`→`face_blend_weight`; `streaming` KEPT).
> All bit-identical/principled; gates green throughout (strict 505/1/4, full 1083→1093).

## Context

Steps A + B are committed. Step B established the loss-rep reads Σ from the operator's C
(single-sourced). Step C (affine-closure → methods) was INVERTED to run after D, then the
inversion DISSOLVED once a literature evaluation (Larsen-Morel 1989, Larsen-Morel-Miller
1987, Bailey-Morel-Chang 2010, Hébert + CFD box/κ-scheme refs) **confirmed the user's
structural claim**:

* **Step ↔ first-order upwind** (donor cell, `ψ̄=ψ_out`, w=1).
* **DD ↔ central — precisely the Keller box scheme** (central advection + trapezoidal
  reaction; box ≡ diamond proven). DD's negative-flux pathology IS the central `Pe>2` wiggle.
* **LD ↔ DG-P1 with upwind flux** (linear/2nd-order upwind; θ=⅓ = P1 Legendre normalization).
* The coefficient triple `(a, 1/denom, w)` is a **generic advection–reaction discretization**
  (swap `|μ|→a`, `Σ_t→r`; formulas unchanged). `w(Σ)` ½→1 is the generic **Péclet/κ-scheme
  blend** (SymPy-verified `w→½` as Σ→0), NOT an SN artifact.
* Genuinely transport-specific = the **ANGULAR axis**: quadrature ordinates + the curvilinear
  angular redistribution (`(1−μ²)/r·∂_μψ`, α-recursion, Morel-Montry τ). Literature splits
  the diffusion-limit analysis into a *spatial* paper (LMM-1987) and an *angular* paper
  (BMC-2010) — space ⊗ angle are factorizable axes.

**Decisive consequence (refutes the original Step C/D framing in
`issue_240_phase2_layer_separation.md`):** the assembly + closure ops are NOT SN-loss-rep
realization that leaked into the scheme — they are a **genuinely generic advection–reaction
*spatial* discretization** that correctly belongs ON the scheme. The assembly does NOT
relocate to the loss-rep; the closure ops stay in the spatial layer → no circular import, no
C-before-D inversion needed.

**User directive (2026-06-16):** task unchanged (home the dangling functions/leaked methods
into their native base) — only the *how*. Each member → **DiscretizationScheme Base**
(generic advection–reaction spatial) or **LossRepresentation Base** (SN angular/sweep); the
cut is the proposal below (up for discussion). The generic advection–reaction form **must be
REALIZED** (Σ-stateless, parameterized by wave-speed + reaction-rate) so SN consumes it
cleanly AND **in preparation for diffusion** — the CONFIRMED next model (standalone diffusion
solver AND consistent-DSA preconditioner, #2). The full **Spatial × Angular tensor product**
is the NEXT campaign, from this clean base.

## The homing (proposed — "what belongs where")

Native-base test: *function of (wave-speed, reaction-rate, source, basis) → generic
advection–reaction → DiscretizationScheme Base; angular/ordinate/sweep-structure →
SN-specific → LossRepresentation Base.*

**→ DiscretizationScheme Base (generic advection–reaction SPATIAL — diffusion-consumable):**
* Reconstruction ops — `cell_average` (`(1−w)·in+w·out`), `source_emission` (`QV·inv/w`),
  NEW `outgoing_face_from_average` (`(ψ̄−(1−w)·in)/w`). Currently DANGLING free functions in
  `spatial/affine_closure.py` → become Base methods/staticmethods.
* Coefficient assembly — `affine_scan_coefficients`, `cell_kernel_batch`,
  `residual_kernel_batch`, the LD Schur (`_schur_terms`/`_kernel_terms`). Already on the
  scheme (`diamond.py`/`linear_discontinuous.py`) — CONFIRMED correctly placed; Σ-stateless (D4).
* Basis — `θ` (LD), the `w`-rule, the moment layout `_ld_source_moments`, traits.

**→ LossRepresentation Base (SN-specific angular + sweep):**
* Sweep/matvec walks (octant/scan/wavefront) — already there.
* The **angular redistribution** (Morel-Montry α-recursion/τ) — FUSED in `diamond.py`'s
  curvilinear `update` + loss-rep `_apply_walk` curvilinear thread. *Cartesian* inlined
  `2ψ̄−in` homes NOW via `outgoing_face_from_average`; the deep **curvilinear angular
  extraction is the NEXT campaign** (the spatial⊗angular fusion already guarded with
  `NotImplementedError`).
* Ordinate/octant structure, per-cell sweep vocab (`CellVisit`/`UpstreamState`/`CellResult`).

## Sequenced sub-steps (each verifiable + committable; route around the 7 reds; NEVER all `tests/sn` — #212)

**D1 — `outgoing_face_from_average` + collapse the inlined closures.** Add the generic
inverse `(ψ̄−(1−w)·in)/w`; route the inlined sites through it: DD `diamond.py:169`/`:346`/
`:384`, LD `linear_discontinuous.py:459`/`:484`, `scan.py:344` + `loss_representation.py:1431`/
`1452`. BIT-IDENTICAL at each scheme's w. Kills the realization-leak duplication. [old-Step-C core]

**D2 — home the reconstruction ops onto the DiscretizationScheme Base.** Move
`cell_average`/`source_emission`/`outgoing_face_from_average` from `affine_closure.py` free
functions → Base methods/staticmethods. Migrate callers (loss-rep scan-solve
`loss_representation.py:2632/2647/2800/2818` + D1 sites). No circular import (spatial layer
owns them; loss-rep already imports the scheme).

**D3 — the rename.** `CellUpdateBase→DiscretizationSchemeBase` FIRST (substring collision),
then `CellUpdate→DiscretizationScheme`, then `cell_update→scheme`. Nexus `rename`, MAIN tree
only (NEVER the worktree). Blast radius: ~29 prod+16 test (`*Base`), ~26 prod+12 test
(`CellUpdate`), ~78 prod+65 test (`cell_update` attr), + docs. Registry `key=` strings
UNCHANGED. Mechanical, bit-identical. (Module-file rename OPTIONAL.)

**D4 — realize the Σ-stateless generic advection–reaction form.** Coefficient methods take
reaction-rate (Σ) as an EXPLICIT param sourced from C (extends Step B); scheme holds NO Σ
state. Interface → `(wave-speed, reaction-rate, source, geometry) → coefficients/cell-update`,
diffusion-ready. `cell_balance.total_xs` already supports this.

**D5 — complete spatial-closure polymorphism in N-D (#239 + N-dim LD).** ⭐ RE-SCOPED by the
user (2026-06-16): *"LD, like any other spatial closure, needs to be able to run on all sweep
strategies."* Complete-polymorphism principle ([[feedback-principled-over-bit-identical]],
"equivalent capabilities across polymorphism"): every spatial closure (DD/Step/LD) must have a
VALID path on every sweep strategy in every dimension; the only permitted exclusions are
STRUCTURAL (declared via `supports()`/`Compatibility(reason)`), never arbitrary code gaps.
**CURRENT GAP: LD is 1-D-ONLY** — it rides `CumprodScan` (1-D scan) but its multi-D kernel
(`_kernel_terms`/`cell_kernel_batch`) RAISES `NotImplementedError` for `len(s_axes) != 1`. DD runs
1-D + N-D (scan + wavefront); LD's N-D is the hole. D5 is now a CAMPAIGN (own design pass + a
proactive test-architect FIRST), two threads:

* **D5a — #239: the N-D DD scan-march onto the coefficient model.** The 2-D `ScanMarch`
  (`loss_representation.py` `_sweep_interior`/`_loss_action_interior`) inlines DD `2g`/`alpha`/`beta`;
  fold onto the coefficient model via `outgoing_face_from_average` (D1) + an **N-D coefficient cache**
  (`CollisionCache`/`from_geometry` is 1-D-only) + the transverse `s_y·ψ_y` folded into
  `source_emission`. Principled ~1-ULP re-baseline (#158-B1; two-paths oracle + cart2d apply snapshot,
  Mode-9 het/non-square/diagonal cubature).
* **D5b — N-dim LD (bilinear multi-axis slope) — SUBSUMES #158 Increment D (#38).** LD's multi-D
  closure is BILINEAR (independent slope per axis) → does NOT fit the scan-march's "x-scan +
  transverse-DIRECT-y" (transverse coupling is a SLOPE, not a face value), so LD-in-N-D rides the
  **DAG wavefront** (`cell_kernel_batch`/`residual_kernel_batch`, the full per-cell bilinear Schur),
  NOT the scan-march. Close the arbitrary code gap: LD's `_kernel_terms`/`cell_kernel_batch` handle
  `len(s_axes) > 1` (the multi-axis Schur). Then LD runs N-D on the wavefront; the scan-march stays a
  DD/Step (slopeless) OPTIMIZATION whose `supports()` returns `Compatibility(False, "bilinear slope
  coupling needs the wavefront")` for LD-in-N-D — a STRUCTURAL exclusion (mathematically true, like
  "can't invert streaming without collision"), NOT a polymorphism violation. Verify: each scheme has
  ≥1 valid strategy per dim; routing-flip (DD≠LD same config); LD's absolute anchor (multi-D LD MMS
  O(h²) — the LM-1989 slope-sign trap, the recurring `Q̂≠0` slope-source coverage gap the elegance/qa
  reviews flagged).

**Connections:** the `SweepStrategy` first-class abstraction (`.claude/plans/sn_sweep_strategy.md` —
protocol `sweep`+`residual`+`supports`; `Compatibility(ok,reason)` selection on `is_cartesian`+`ndim`)
is the natural home for the `supports()` exclusions; #158 Increment D (#38) is SUBSUMED. `is_affine_scannable`
likely needs to be dimension-aware (LD: True 1-D / False multi-D) OR `supports()` carries it.

**D6 — docs + next-campaign issue.** Archivist: add Step/DD/LD ↔ advection-scheme
correspondence to `docs/theory/discrete_ordinates.rst` + scheme docstrings (upwind/box-central/
DG-P1-upwind, box≡diamond nuance, `w(Σ)` Péclet coupling, spatial⊗angular factorization,
lumped-step-vs-step-characteristic). File a NEW issue for the full **Spatial × Angular tensor
product** (extract curvilinear angular redistribution → distinct AngularScheme).

## Verification (test-architect spec)

* **D1/D2/D3: BIT-IDENTICAL** — strict DriftWarning gate (CORRECT path
  `tests.sn.regression._regression_assert.DriftWarning`):
  `python -O -m pytest tests/sn/sweep/core tests/sn/solve -W "error::tests.sn.regression._regression_assert.DriftWarning"`
  (baseline 505/1/4). Rename: migrate name/attr-pinned tests; oracle/solver tests survive.
* **D4: bit-identical** prod + a NEW capability gate (explicit reaction-rate ≠ mesh Σ → matches
  direct build), mirroring Step B's removal-form teeth.
* **D5a (#239 DD scan-march): ~1-ULP** — two-paths oracle `test_scanmarch_sweep_equals_oracle`/
  `_residual_equals_oracle` + a cart2d apply arm on `TestT4bPreT4RegressionSnapshot` (Mode-9:
  het/2G-asym/non-square, diagonal cubature).
* **D5b (N-dim LD): NEW capability** (not bit-id) — multi-D LD MMS O(h²) on the wavefront (the
  LM-1989 slope-sign trap; MUST include a non-vanishing slope source `Q̂≠0`, the recurring coverage
  gap); a per-(scheme×strategy×dim) coverage matrix (every scheme has ≥1 valid strategy per dim;
  `supports()` rejects LD-on-scan-march-in-N-D with the structural reason, asserted as a NEGATIVE
  gate); routing-flip DD≠LD on the same N-D config. ⚠ proactive **test-architect FIRST** — this is a
  cross-strategy/cross-dim carve (the L17 trigger) AND a new L0/L1 reference (multi-D LD MMS).
* **Basis/contract gate (vv L11):** thin scheme exposes basis/traits; recon ops reachable as
  Base methods; Mode-8 negative arms use `pytest.fail`.
* **Routing-flip (#158-B3 analogue):** DD vs LD same 1-D het 2G → DIFFERENT + each scheme's
  absolute anchor (LD `TestLDLinearExactness`/`test_mms_ld_slab` O(h²); DD k∞).
* **Route-arounds:** `-k "not (vacuum_bulk_bit_identical_1d and SPH) and not (sphere_1g_apply_bit_identical or sphere_2g_apply_bit_identical) and not test_2d_mesh_resolution and not two_d_cartesian_loss_action"`.
  Dirs: `tests/sn/operators tests/sn/spatial tests/sn/sweep/core tests/sn/sweep/cartesian_2d tests/sn/solve`.
* Per sub-step: elegance-enforcer + qa; Sphinx clean; `tests._harness.audit` exit 0.

## Execution discipline
method-implementer for D1/D2/D4/D5 against the spec; elegance + qa before each commit;
archivist for D6. One commit per sub-step + a `chore(claude)` records commit; trailer
`Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`; explicit paths only
— NEVER the 3 forbidden untracked / parallel-session set / `docs/_build/`. `.claude/agent-memory/*`
→ chore commit. Branch NOT pushed/merged. **Step C (#48) is SUBSUMED by D1+D2.**

## EXIT
The full Spatial × Angular tensor-product campaign (the D6 issue). Then resume
`issue_158_spatial_cellupdate_carve.md` §"Increment C" (#37) → Step-4 (#36) + Inc D (#38);
ff-merge → `main` + apply the Phase-0 hand-off.
