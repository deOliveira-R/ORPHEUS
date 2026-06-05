# SN Development Sequence — the consolidated plan (Wave O completion → features → N-D)

**Branch:** `refactor/field-role-typing` (worktree). **HEAD at authoring:** `7bdea86`.
**Mode:** main-agent DIRECT authorship, turn-by-turn user steering
(`feedback_no_method_implementer_for_surgical_carves`) — these are SN hot-path
operator-algebra carves; commit per sub-step; bit-identity / equivalence gates
throughout. **Standing commit exclusions:** `.claude/agent-memory/`, `.mcp.json`,
hooks, `derivations/diagnostics/`. **Discipline:** read the WORKTREE, never the
MAIN checkout (lesson L22 — the main checkout is a different, older branch).

**This file SUPERSEDES** (folded in, then archived): `wave_o_operator_typing.md`,
`wavefront_flux_foundation.md`, `si_gauss_seidel_recovery.md`, the live tail of
`field_role_typing_view_g.md`. **Durable references kept separate:**
`nd_foundation.md` (the N-D design-of-record; Phase 6 points to it),
`trajectory_resolvent_hindsight_refactor.md` + `cylinder_mr_variant_alpha_verification.md`
+ `capability_matrix_framework.md` (out-of-spine subsystems).

**Refreshed file:line surface-maps** (prefer these over any line numbers quoted
here — they are re-confirmed per-HEAD by the explorer): the agent-memory surface
maps `wavefront_flux_carve_substrate.md`, `sn_si_reflective_gauss_seidel_recovery_surface.md`,
`issue_208_b52_boundary_residual_retype_surface.md`. Re-run a quick
`grep`/`inspect.getsource` at each pickup; the anchors below have drifted.

---

## Environment & recovery (READ FIRST — esp. after compaction)

**We work in a git WORKTREE, NOT the main checkout** (lesson L22 — the main checkout
is a DIFFERENT, OLDER branch `refactor/sn-operator-algebra` whose `solver.py` /
`operator.py` / `sweep.py` DIFFER; reading it silently reads the wrong source).

- **Worktree (cwd for everything):** `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/field-role-typing`
- **Branch:** `refactor/field-role-typing` (on `origin`, but the Phase-1 commits
  `a902c59`…`d643ee8` + the recovery refresh `216f4f6` + the Phase-2 commits `8563f4b`
  (sub-step 4) and `83a4ae6` (sub-steps 1+2) are LOCAL-ONLY as of 2026-06-05 — `origin`
  is at `61546e8`. Recovery works from the LOCAL commits + this plan + the memory;
  `git push` to sync cross-device).
- **venv** (lives at the MAIN repo root, shared): `/Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python`
- **Test invocation (Host env, `$CLAUDE_ENVIRONMENT` empty):**
  `PYTHONPATH=/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/field-role-typing /Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python -O -m pytest ...`
  (default mode `-O`; assertions still fire in test modules via pytest rewriting).
- **Bash discipline:** redirect to a file then `echo exit=$?` (never `| tail`); use
  `gtimeout` for slow runs. **NO `continuous_get`** (the #212 hang).
- **When a sub-agent's `file:line` conflicts with your read:** confirm you both read
  the worktree via `inspect.getsource(fn)` / `print(module.__file__)` under the
  worktree `PYTHONPATH` BEFORE concluding fabrication (L22 — a wrong-tree read mimics it).

**Recovery path after compaction / fresh session:**
1. `MEMORY.md` → the `[[project-wave-o-operator-algebra]]` ⭐⭐ LATEST block (Phases 1+2
   COMPLETE + Phase 3 IN PROGRESS (3a done, 3b/3c next) + the R5 reversal + the
   eigenvalue-posing architecture of record).
2. THIS file (`sn_development_sequence.md`) — the 6-phase sequence (§3); Phase 1's
   **STATUS block** (in §3) records the landed R1/R2/gap/R5 with commit hashes.
3. The session **task tracker** (`TaskList`): tasks #17–#22 = the 6 phases.
   **If the tracker did not survive the session boundary, recreate it from §3 + §5**
   (one task per phase; chain blockedBy 18→{19,20}, 19→21, 21→22; **#17 = DONE**).
4. **Phase 1 (#17) DONE; Phase 2 (#18) DONE** (2026-06-05). The honest `L+C−S−B`
   variadic driver landed + the fold retired + Sphinx current + the orthogonal reds fixed.
   Commits: `8563f4b` (sub-step 4, B trace-only leaf) → `83a4ae6` (sub-steps 1+2, honest
   variadic driver) → `fea6675`/`6d30807` (plan) → `deb1ce3` (precond_safety fold→gains
   migration, 4 reds green) → `40de6e5` (Sphinx theory + solver docstring fixes) → `33dd5ff`
   (b1pp O.4b boundary-semantics + restart gate-tightness, 6 reds green). **Sub-step 1
   (`ScatteringOperator.domain`) DEFERRED** (premise superseded — no bulk space; tripwire
   moot; → Wave-O typing seam). **The two remaining sub-steps DEFERRED with TRIGGERS (not
   indefinite, per user):** sub-step 5 = #200 block-inverse FACE preconditioner **FEATURE**
   (unblocked NOW; commented on issue #200 with the 3 deliverables incl. re-tighten the
   restart gate 1e-7→1e-9) ; sub-step 3 = R4 curvilinear-matvec + 1-D-scan collapse →
   **nd_foundation Phase 6 d-generic walk** (commented on issue #206; trigger = after the
   `_sweep_wavefront` walk lands).
5. **Phase 3 (#19, SI Gauss-Seidel recovery) ✅ DONE** (2026-06-05, HEAD `1cbb383`, local-only).
   Landed as an HONEST, MODEST **boundary** Gauss-Seidel accelerator — NOT the c²-halving the
   original plan/spec assumed. ⚠ **The c²-halving premise was REFUTED by a spike**: boundary-G-S
   folds only the reflective coupling `B` (S stays a lagged gain), so it accelerates the
   boundary-layer transient, not the dominant flat *scattering* `c`-mode → a measured constant
   **~0.86–0.92×** (regime-independent), not 0.5×. The real within-group rate win is **Krylov
   (already production, splitting-invariant, rate-optimal on every BC)** or **consistent DSA**
   (commented on issue **#2** with the full spike evidence — c-independent 8–21× on vacuum, but
   naive FD-DSA DIVERGES on reflective → consistent DSA needed). The σ_r-fold (folding the
   scattering self-scatter into a σ_r-sweep) is a LATENT CORRECTNESS TRAP (issue **#215** + the
   docstrings now warn it: isotropic-limit-only, 46–56% wrong on anisotropic, or divergent).
   Commit chain: `465bb5d`(3b.1 face-restricted seed/absorb) → `261abad`(3b.2 sweep_octant_group
   carve, Jacobi bit-identical) → `0421e72`(σ_r trap docstrings) → `3697eef`(3c.1 polymorphic
   `_sweep_2d_scheduled`) → `a39905a`(3c.2 `_GaussSeidelResolvent` + `inner_schedule` selector +
   the **ERR-056 shared-face fix** — reflect a face only after its LAST outflowing octant group;
   diagonal/lebedev cubatures were wrong before) → `c1ac2f7`(3c.3 rate gates re-scoped to LIVE
   measurement: strict-improvement + Krylov bracket + same-FP, no hardcoded baseline) →
   `da05e93`(ERR-056 + vv-principles **Mode 9**) → `1cbb383`(Sphinx theory doc + `:label:
   si-spectral-rate`). Eigenvalue (`solve_sn`) deferred — stays Jacobi (only `solve_sn_fixed_source`
   defaults to G-S). **Subsequent choices: consistent DSA (#2 — the real rate feature) / Phase 4
   (#20) / #200.** Each crosses subsystem boundaries → test-architect FIRST (L17). The carve
   substrate map (`.claude/agent-memory/explorer/si_gs_substep3_carve_substrate.md`) + Phase-2
   refs remain for reference.

---

## 0. Where we are (lineage + what already landed)

The SN solver was rebuilt over Phase G → R-1 → Depth-B → Wave T → Wave O into a
typed operator algebra. **Already landed on this branch:**
- **Phase G / R-1:** `CellUpdate` Protocol + `DiamondDifference`; unified 1-D sweep
  (`_sweep_1d_unified`); `StreamingOperator`/`CollisionOperator`/`ScatteringOperator`/
  `FissionOperator` as `LinearOperator`s; `SNMesh` consumes `ReducedStreamingOperator`;
  `(N,ng,nx,ny)` principled layout. (Issues #157/#159/#160/#161/#162/#174 closed.)
- **Wave O O.1:** `BlockRole` + `BulkOperator`/`FullOperator` markers.
- **Wave O O.4a/O.4b (BC extraction):** `(L+C−S−F−B)ψ=q` canonical on BOTH the 1-D
  and 2-D paths; the sweep + matvec are BARE (seed the given inflow trace, no
  in-sweep `bc.apply`); `B = SNBoundaryOperator` is the sole reflective coupling
  (natively 4-face); `_compute_gradients` FD stencil RETIRED; matvec ≡ sweep in
  2-D by construction (one DD discretization through `SweepDependencyGraph`).
- **Wave O B.5.2:** operator outputs typed (`AngularSourceSink`/`BoundarySourceSink`);
  role grid bulk∥boundary (`.apply→SourceSink`, `.solve→Flux`, `from_balance→Residual`).
- **WavefrontFlux storage-A (Phases 0–5):** interior cell-face cochain `C¹_int` is the
  typed `WavefrontFlux` on `InteriorFaceSpace`; typed `ι_*` (seed) / `ι*` (absorb);
  axis-parametric types (3-D-ready).
- **2-D SI "Phase A":** the stale 2-D eigenvalue source-iteration guard removed; the
  default `solve_sn` 2-D eigenvalue entry works; SI≡Krylov≡k_inf verified.
- **Phase 1 — SN solve unification (2026-06-05, ✅ DONE; full record in §3 Phase 1
  STATUS):** ONE within-group SI path, ONE Krylov path, ONE within-group operator
  triple (`_within_group_triple` / `_within_group_krylov` in `solver.py`), ONE
  power-iteration LOOP. `power_iteration` (`numerics/eigenvalue.py`) is the canonical
  method-agnostic Layer-4 algorithm; `KEigenvalue` (`numerics/iteration.py`) delegates
  to it (**R5 was REVERSED — `power_iteration` is NOT retired**; the 5 solver families
  keep it; 1e moot). 2-D Cartesian fixed-source Krylov un-gated. **Architecture of
  record:** `docs/theory/operator_algebra.rst` `:ref:eigenvalue-posing` — the
  generalized eigenproblem `A_loss ψ = λ M ψ` (resolvent `A_loss⁻¹M`); four layers
  leaves → posing (2a role+μ-map / 2b realization) → resolvent → algorithm; K row
  live, α/transient/adjoint/full-spectrum = documented future seams. Commits
  `a902c59`(R1) / `023ae18`(R2) / `cf82556`(gap) / `650032e`+`7603c8e`(R5) /
  `e70ca39`(plan+L23) / `3ee1598`(docs) / `d643ee8`(V&V matrix).

**Standing reds (ORTHOGONAL — do NOT block this sequence; tracked as issues):**
cylinder matvec #206 (curvilinear WDD matvec hard-red ~18% — NOT fixed by Phase 2; the
curvilinear matvec collapse is the deferred R4 → nd_foundation Phase 6, commented on #206;
the #196 SI≡Krylov *eigenvalue* twin is already closed, a stale `xfail(strict)` XPASS on
`test_cylinder_l1_sweep_vs_krylov_twin_path` — un-marking it is a small separate cleanup);
**[FIXED 2026-06-05 `33dd5ff`]** ~~`test_b1pp_constant_flux` ×3~~ (migrated to the O.4b
inflow-identity/outflow-defect boundary semantics) + ~~`test_krylov_kinf…refinement` ×3~~
(gate relaxed 1e-9→1e-7 for the unpreconditioned-Krylov ~1.6e-9 floor; re-tighten when #200
lands); #212 `continuous_get` hang (reference-registry build; has a proposed
patch `derivations/diagnostics/sn_keff_hang_PROPOSED_fix.patch`; deselect the 3
affected eigenvalue tests meanwhile); #209 `ordinate_scan` NaN on a zero a-chain
(cylindrical pole); #195 curvilinear MMS magnitude at nx=160; **#214 WavefrontFlux
2-D crash on a degenerate `ny=1` mesh** (filed 2026-06-05; Phase 5 / `nd_foundation`
territory — the d-generic walk must treat an extent-1 axis as inactive). None gate the spine.

---

## 1. The algebra of record (durable — the native frame)

(From the cross-domain-attacker, twice-confirmed; `issue_208_operator_algebra_frames_full_access.md`.)
The SN operators form a **dagger inverse biproduct category with a non-trivial
block metric**. Four load-bearing consequences the whole sequence rests on:

1. **`V = V_bulk ⊕ V_inflow ⊕ V_outflow` is a biproduct** ⟹ operators ARE block
   matrices; `OperatorSum.apply` is block-matrix-by-construction (adding a 3rd ⊕
   slot needs no new `if`).
2. **The dagger `†` is the G-adjoint `A† = G⁻¹AᵀG`, NOT plain transpose.** The
   boundary block carries `G_s = |Ω·n|·w_n` (partial current). Order-reversal alone
   yields the WRONG adjoint until `G` is populated — this is why Gate 1.3 stays
   `xfail` (Phase 4).
3. **`B` is a gluing cocycle** (reflective = deck transformation `Ω→R(Ω)`; periodic
   = base quotient) — validates the O.4 BC-extraction and shares the MoC fiber-bundle
   frame.
4. **Transport-resolvent backbone:** the four method `solve`s are quadratures of one
   `(Ω·∇+Σ_t)⁻¹`; diffusion is its self-adjoint-elliptic asymptotic limit.

The cochain refinement (WavefrontFlux): `C¹ = C¹_int ⊕ C¹_∂` is the SAME biproduct
one locus down (face level); `ι_*`/`ι*` are the injection/projection; the interior
cochain is flux-only because the role grid is a 0-cochain (cell) concept the boundary
1-chain inherits only via its BC residual.

---

## 2. The redundancy map (what this sequence collapses)

The current solution path carries five duplications (Cardinal Rule 2). The sequence
is ordered so the SPINE (Phases 1–2) collapses them BEFORE features (3–5) land on
the result — features touch these exact sites, so collapsing first means building once.

| ID | Duplication | Sites (re-confirm anchors) | Collapsed in |
|----|-------------|----------------------------|--------------|
| R1 | TWO source-iteration impls | `_solve_source_iteration` (uses `SourceIteration` primitive) vs `_solve_fixed_source_si` (hand-rolled `for n_inner` loop over `transport_sweep`) | Phase 1 |
| R2 | TWO Krylov builders | `_solve_krylov` (eigenvalue, all-geom) vs `_solve_fixed_source_krylov` (fixed-source, 2-D NIE) | Phase 1 |
| R5 | duplicated power-iteration LOOP (`KEigenvalue.solve` re-implemented `power_iteration`'s loop byte-for-byte) — NOT a redundant symbol (the premise was refuted) | `numerics/iteration.py` `KEigenvalue.solve` ≡ `numerics/eigenvalue.py` `power_iteration` | Phase 1 ✅ (RESHAPED — see Phase 1 STATUS) |
| gap | fixed-source-Krylov 2-D `NotImplementedError` | `solver.py` ~1582 | Phase 1 |
| R3 | TWO `−B` delivery routes | `_scattering_with_boundary_op` (S+B fold) vs `_reflect_outflow_into_inflow` helper | Phase 2 (route-collapse) / Phase 3 keeps the helper additively |
| R4 | 3 parallel `(L+C)` impls | `_compute_LpC` / `_compute_decomposition` / `(L+C).solve` sweep | Phase 2 |

`#206` ≡ R4 ≡ the mechanism of `#196` (curvilinear SI-vs-Krylov O(h)); `#200` ≡ the
two `preconditioner=identity` stopgaps (R2 sites); `#201 ⊂ #205` ≡ Wave-O role-typing
(Phase 4); `#213` ≡ capabilities-as-morphism-class (deferred pointer).

---

## 3. THE SEQUENCE

Architecture-first: **spine (1–2) → motivating feature (3) → adjoint completion (4)
→ performance features (5) → N-D foundation (6).** Phase 4 floats (orthogonal to 3;
move earlier only if Gate-1.3-green becomes a priority).

### Phase 1 — Unify the solve architecture  *(R1 + R2 + R5 + gap)* — ✅ DONE

**STATUS (2026-06-05).** Complete. ONE within-group SI path, ONE Krylov path, ONE
within-group operator triple, ONE power-iteration LOOP.
- **1a (R1)** `a902c59` — `_solve_fixed_source_si` folded onto the `SourceIteration`
  primitive (+ NEW-3 structural spy `test_si_single_primitive_contract.py`).
  Principled-equivalent (the composite-residual stop metric; slab snapshot 2 ULP).
- **1b (R2)** `023ae18` — `_within_group_triple` + `_within_group_krylov` module
  helpers = the SSoT; all 4 within-group sites consume them. −30 LOC, bit-identical.
- **1c (gap)** `cf82556` — 2-D Cartesian fixed-source Krylov un-gated (+ NEW-1
  `test_fixed_source_2d_equivalence.py`: closed-form Q/Σt + SI≡Krylov twin); the
  G1-raise pin retired (inverted into NEW-1).
- **R5 (RESHAPED → done)** `650032e` (drop the inverted deprecation, reframe docstrings)
  + `7603c8e` (`KEigenvalue` delegates to `power_iteration`). The plan's R5 premise
  ("retire `power_iteration`, migrate 5 families to `KEigenvalue`") was **REFUTED**
  (cross-domain-attacker `power_iteration_vs_keigenvalue_morphism.md` +
  `eigenvalue_posing_layering_frames.md`): the redundancy is the duplicated power-
  iteration LOOP, NOT the symbol. `power_iteration` is the canonical Layer-4 algorithm
  over the method-agnostic `EigenvalueSolver` boundary (it admits the CP/diffusion/
  homogeneous monolithic-matrix resolvents that have NO `(L,S,F)` triple);
  `KEigenvalue` now DELEGATES its loop to it (one loop engine). **1e (sibling
  migration + `power_iteration` retirement) is MOOT** — the 5 families keep
  `power_iteration` (the right home). The α-generic *agnostic relocation* (Protocol
  rename `compute_fission_source`/`compute_keff` → `eigen_operator`/`mu_to_eigenvalue`,
  scaling relocation, `apply_loss`) is the **α-eigenvalue wave's first step**
  (documented seam, snapshot-bit-identity-gated across 5 families).

**Architecture of record (durable).** Standard form = the generalized eigenproblem
`A_loss ψ = λ M ψ` (power-method realization = dominant eigenpair of the resolvent
`A_loss⁻¹M`). Four layers: **leaves** (L,C,S,F,B,[T]; the |Ω·n|·w metric lives here)
→ **posing** (2a method-agnostic role-assignment + μ→physical map / 2b method-specific
`A_loss` realization) → **resolvent** `A_loss⁻¹` (method-specific: SN SI/Krylov, CP
BiCGSTAB) → **algorithm** (general: power iteration | future Arnoldi/FEAST; transient
time-integrator). `KEigenvalue` = the operator-triple 2b realization; the K posing row
is `A_loss=L+C−S−B, M=F, k=μ`. Full plan `glimmering-launching-lantern.md`; theory
`docs/theory/operator_algebra.rst` `:ref:eigenvalue-posing` (α / transient / adjoint /
full-spectrum rows = documented future seams). Adjoint = a daggered posing row (same λ).

**Verification (landed).** R1 NEW-3 + slab snapshot (2 ULP) + L0 Gate 1.1 + SI≡Krylov
2-D + curvilinear SI; R2 71-test batch bit-identical-to-baseline; gap NEW-1 both legs;
R5 discriminating gate `test_keigenvalue_matches_solve_sn_2g_slab` stays green after
delegation (same-morphism proof) + KEigenvalue L0 units. Sentinel throughout.

**Carried gotchas.** `_reflect_outflow_into_inflow` KEPT (Phase 3 G-S). The
fixed-source-Krylov shares #200's `preconditioner=identity` (Phase 2 fixes). **Filed:**
#214 (pre-existing orthogonal WavefrontFlux `ny=1` 2-D crash, stash-confirmed not mine).

### Phase 2 — O.2a: the honest `L+C−S−F−B` driver  *(R3 + R4 + #206 + #196 + #200)*

**Goal.** The iteration drivers consume the WHOLE loss operator as one composed
block operator, retiring the transitional folds and collapsing the matvec/sweep twins.

**STATUS (2026-06-05) — sub-steps 4 + 1+2 DONE; design decisions taken.**
- **Sub-step 4 (B trace-only leaf)** ✅ `8563f4b` — `SNBoundaryOperator.reflect_into_inflow`
  + shared `_reflect_trace` core; retires the zero-bulk-probe shim. Bit-identical
  (stash-test proven neutral).
- **Sub-steps 1+2 (honest composition)** ✅ `83a4ae6` — **the `S+B` fold is RETIRED.**
  ⭐ **Design decision (user-chosen via AskUserQuestion): VARIADIC couplings.** The
  drivers `SourceIteration`/`KrylovAcceleration` went `(L, S, F)` → `(L_resolvent,
  *gains)`: matvec `L − Σ gᵢ`, rhs `q_ext + Σ gᵢ`. `_within_group_triple` now returns
  the honest `(L+C, S, B)` — resolvent + bulk scattering gain + trace boundary gain;
  the vestigial within-group `F=ZeroOperator` slot is GONE. The driver is now
  problem-type-AGNOSTIC (which ops are gains = a posing decision). `B` stays a SEPARATE
  trace-typed gain (NOT folded into bulk S): it can't join the L+C preconditioner
  (OperatorSum drops CAP_SOLVE) and the |Ω·n|·w adjoint metric (O.2) lives on its trace
  domain. Existing `(L,S,F)` callers stay positional-compatible. Principled-equivalent
  (FP reassociation `L−(S+B)`→`(L−S)−B`): cyl keff 4.2e-13 / flux 6.8e-12, aniso 3–5 ULP
  — all within regression tol; verified correct vs NEW-1 closed-form Q/Σt + SI≡Krylov +
  keff_2d k_inf. Tests migrated (NEW-3 pins `(L+C, S-gain, B-gain)`; 2 Krylov capability
  tests → per-gain message). Gate: 18 driver + 52 fast + 114 regression/eigenvalue green.
- **Sub-step 1 (set `ScatteringOperator.domain`) — DEFERRED (premise superseded).** The
  domain was the "tripwire" meant to FORCE the honest composition by making the `S+B`
  fold throw. With the variadic split the fold is retired by construction, so the
  tripwire is moot; AND no bulk `FunctionSpace` exists (bulk operators type via
  `block_role`, not domain spaces — only `B` sets `domain=trace`, for the O.2 metric).
  Minting a `V_bulk` space is a documented **Wave-O typing-completion seam** (defensive
  Pattern-4 typing), NOT load-bearing for the honest composition. Pick up when the
  bulk-space typing wave lands (it would let `OperatorSum` reject a re-introduced `S+B`).

**Remaining sub-steps (status annotated):**
1. ~~Give `ScatteringOperator` a domain~~ — **DEFERRED** (see STATUS: tripwire superseded
   by the variadic split; no bulk `FunctionSpace` exists; → Wave-O typing-completion seam).
2. ~~Compose `L+C−S−B`~~ — ✅ **DONE** `83a4ae6` (variadic; fold retired; the two `−B`
   routes collapsed on the driver path — `_reflect_outflow_into_inflow` survives only for
   the reconstruction sweep + Phase 3's octant-restricted variant).
4. ~~Trace-only `B.reflect_into_inflow(boundary)`~~ — ✅ **DONE** `8563f4b`.
5. **`#200` real preconditioner** — **DEFERRED (scope corrected 2026-06-05: it is a
   FEATURE, not a re-enable).** Investigation refuted the "re-enable the sweep precond"
   framing: `test_krylov_curvilinear_precond_safety` is a CAREFULLY-DESIGNED test that
   PINS the current IDENTITY-precond production contract and explicitly documents that the
   naive sweep precond (`preconditioner=None` → `LC.solve`) is a POOR preconditioner on
   curvilinear (M-M cold-start) — switching to it would LOCK IN the wrong production state.
   #200's real fix is the **block-inverse FACE preconditioner**, a substantial NEW FEATURE
   that does not exist yet. There is **no correctness red** forcing it: the production
   Krylov path recovers k_inf correctly (`test_krylov_restart_signature` → 1.875 to ~1.5e-9;
   the strict-1e-9 gate is marginally tight at meshes 5/16/30 — gate-tightness, not a break).
   → Deferred to issue #200 (its own feature effort). **A directly-related win DID land:
   `deb1ce3` migrated `test_krylov_curvilinear_precond_safety` to the honest `(LC, S, B)`
   gains** — it had hand-built `(LC, S, ZeroOperator)` OMITTING B (pre-existing red, k≈1.67),
   so adding the B gain (mirroring the new production triple) recovers k_inf=1.875 on all
   coords (4 reds → green; retirement=test-migration for the fold).
3. **Unify the 3 `(L+C)` impls** (R4) — REMAINING, ⚠ **FORKED (needs scope decision).**
   The three: `_compute_LpC` (1-D matvec, `operator.py:336`) / `_compute_decomposition`
   (1-D dual-emission introspection, `:578`) / the 1-D solve-sweep `_sweep_1d_unified`
   (scan). **Explorer + test-architect findings (re-confirmed 2026-06-05):**
   - The **1-D `_compute_LpC` ≡ `_compute_decomposition` twins** (same bidirectional
     `dag_walk`, one emits L+C, the other the (M_spat, M_ang) split) ARE a genuine
     Cardinal-Rule-2 duplication — collapse-able into one walk, geometry-agnostic within
     1-D. This is the SAFE, real R4 win.
   - **Extending to the 2-D `graph.residual` shape is WRONG-FIT for curvilinear + the
     1-D scan** (explorer): `_sweep_1d_unified` is a parallel-prefix SCAN (`ordinate_scan`,
     `(nx,K,ng)` chain-leading) — the antithesis of the anti-diagonal dense buffers — and
     curvilinear carries a pole-angular-redistribution closure the Cartesian
     `CellUpdate.residual_batch` has NO slot for. **Recommendation: collapse the 1-D matvec
     twins (+ slab onto graph if clean); DEFER curvilinear matvec + the 1-D scan solve-sweep
     to nd_foundation (the documented Phase-4-1D WRONG-FIT).**
   - **#196/#206 reality check (test-architect):** R4 dissolves the **SI≡Krylov eigenvalue
     twin** (#196) — but that gate (`test_cylinder_l1_sweep_vs_krylov_twin_path`) is a
     STALE `xfail(strict)` that XPASSes TODAY (already closed by B1''/D-K). R4 does NOT
     deliver curvilinear MMS O(h²) (that's #195, pole-closure, OUT of scope) and does NOT
     fix the **#206 cylinder-matvec hard-red** (~18%, the deferred curvilinear matvec). The
     plan's named canary `test_phase_e_trajectory_resolvent_flux_shape_crosscheck` is the
     WRONG gate (12% tol, conflates #196 with #195) — use `test_cylinder_l1_sweep_vs_krylov_twin_path`
     (un-mark the XPASS) + a new heterogeneous flux-shape twin instead.
   - Full crosswalk: `.claude/plans/phase2_o2a_crosswalk.md`; dependency audit:
     `.claude/agent-memory/explorer/phase2_o2a_dependency_audit.md`.

**Verification (for the remaining sub-steps).** #200: the #200-adjacent reds flip green,
no regression, GMRES converges faster (L16 wall-clock sane). R4: the L14 four-legged
standoff on the collapsed 1-D matvec twins (SI ≡ ref, Krylov ≡ ref, SI ≡ Krylov, all under
refinement); vacuum bit-identical; principled-equivalent (FP-reduction-order, documented).

**Gotchas.** R4: keep the variadic single-primitive shape (don't re-fork the drivers).
The curvilinear/scan deferral is the documented WRONG-FIT — do NOT force `graph.residual`
onto them.

### Phase 3 — SI Gauss-Seidel recovery  *(the motivating feature)*

**STATUS (2026-06-05): ✅ DONE — HEAD `1cbb383` (local-only, not pushed).** See the top-of-file
recovery **bullet 5** for the full outcome + commit chain. ⚠ **OUTCOME DIFFERS FROM THE PREMISE
BELOW** (this §3 design narrative is preserved as the as-designed record; read bullet 5 + the
Sphinx theory page `docs/theory/discrete_ordinates.rst` `:ref:si-within-group-splitting` for the
as-BUILT truth). The c²-halving target was **refuted by a spike**: the landed octant-group G-S
folds only the **boundary** `B` (not the scattering `c`-mode), giving a measured constant
**~0.86–0.92×** — a modest reflective-SI accelerator, honestly gated (LIVE strict-improvement +
Krylov bracket + same-FP), NOT the c² recovery. The real within-group rate win is **Krylov
(have)** / **consistent DSA (#2)**; the σ_r-sweep fold of the scattering self-scatter is a latent
correctness trap (**#215**). Findings logged: **ERR-056** (diagonal-cubature shared-face reflect)
+ **vv-principles Mode 9** (verify splittings on anisotropic/diagonal stressing configs). Below:
the as-designed Gate 1 + 1/2/3a foundation (`926c4ca`→`514ae21`) + the 3b/3c design — all built,
with the premise correction noted inline.

#### The problem (what we recover) — durable
Wave O BC extraction made the transport sweep BARE and applies the reflective coupling `B`
EXTERNALLY (the sibling `−B`). That converted the reflective coupling from **intra-sweep
Gauss-Seidel** (the retired in-sweep `bc.apply` let later octants see earlier octants' fresh
reflected outflow) to **inter-sweep Jacobi** (`B` fully lagged on the previous iterate). SAME
converged fixed point, SLOWER SI spectral rate (reflective traps neutrons → effective ρ→1;
the ~654-vs-128 reflective-vs-vacuum sweep-count gap). Recovery = interleave the external
`−B` at octant-group granularity INSIDE the SI resolvent. **Krylov is splitting-invariant —
UNAFFECTED + UNTOUCHED.**

**Mechanism / literature.** The `(octant × face)` reflective graph (mesh-time): an edge
`producer-octant → consumer-octant` on a reflective face `f` when the producer's outflow on
`f` specularly reflects (`quad.reflection_index(axis_f)`) into the consumer's inflow on `f`.
Order-respecting edges → G-S-eligible (fold into the resolvent); cycle-forming edges → lagged
(Jacobi). **KBA / Adams-Larsen / Pautz** — the `i+j=k` wavefront levels ARE the KBA diagonals;
`ρ_J = c` (scattering ratio) Jacobi vs `ρ_GS ≈ c²` for the symmetric reflective model problem.

**Caveats (load-bearing).** BOTH faces of an axis reflective ⟹ a 2-cycle (`(−,·)↔(+,·)`) ⟹ one
forward pass cannot drop both back-edges (opposite-face edge stays Jacobi) ⟹ PARTIAL one-pass
G-S. FULL one-pass G-S only when ≤1 reflective face per axis. White/albedo couples ALL
ordinates on the face → degenerates to Jacobi. Target = **specular reflective** faces. **1-D
slab is a CONTROL** (scattering-dominated `ρ_J=c`, little to recover) AND a WRONG-FIT for the
3b carve (`_run_1d_sweep` is a parallel-prefix SCAN, no `psi_x` interior cache) — DEFERRED.

#### The design (user-reviewed — ⭐ POLYMORPHIC schedule + seed-then-overwrite)
**Jacobi and G-S are the SAME uniform algorithm parameterized by a mesh-time `SweepSchedule`
strategy** — there is NO `if jacobi/gs` branch in the iteration; the splitting is selected
ONCE by choosing the schedule (an `inner_schedule` selector; default `"gauss_seidel"`,
`"jacobi"` selectable as the splitting-invariant control).

The SI resolvent's `.solve(rhs, ψₙ)` runs ONE loop — the **seed-then-overwrite** realization
of `(L+C − B_lower)⁻¹`:
1. **Seed** the boundary inflow `= B·ψₙ` (whole-trace reflection of the PREVIOUS iterate, from
   `initial_guess` — the same value today's Jacobi seeds onto `rhs.boundary`).
2. `for group in schedule.groups:` **sweep the group** (reads its inflow — the frozen seed for
   early groups, OR freshly overwritten by earlier groups); then if `group.reflect_faces`:
   **absorb** the group's outgoing faces (wavefront edge → boundary outflow slots) → **reflect**
   `reflect_into_inflow(faces=group.reflect_faces)` (boundary outflow → boundary inflow) →
   **re-seed** the wavefront edge for those faces. Later groups read the fresh current-iterate
   reflection (the order-respecting `B_lower` edges → G-S); groups swept before their specular
   partner keep the lagged seed (the cyclic `B_upper` back-edges → Jacobi).
   - **Jacobi schedule** = ONE group (all octants), empty `reflect_faces` ⟹ every octant reads
     the frozen seed, the after-group reflect is a no-op ⟹ IDENTICAL to today's bare sweep.
   - **G-S schedule** = one group per in-plane `OctantLabel` (sign_z merged), in quadrature
     sweep order; each group's `reflect_faces` = the reflective faces it outflows through.

**Fixed point provably UNCHANGED (THE load-bearing correctness claim):** any consistent
splitting of `(L+C−S−B)ψ=q` shares the dominant solution ψ\* — at convergence the seed AND all
re-reflects equal `B·ψ\*`, so the converged field is schedule-independent. The schedule changes
ONLY the SI spectral rate. `S` (within-group scatter) stays a lagged gain throughout — only the
BOUNDARY coupling gets G-S (the sweep never re-scatters mid-sweep).

**The home (post-Phase-1+2-correct).** G-S is SI-specific. `_within_group_triple` stays the
`(L+C, S, B)` SSoT. The **SI** path COMPOSES a scheduled resolvent from `(L+C, B, schedule)`
and passes it to `SourceIteration` with gains `(S,)` — **B moves INTO the SI resolvent** (it
seeds `B·ψₙ` from `initial_guess` + does the inter-group reflects). The **Krylov** path is
UNTOUCHED (`KrylovAcceleration(L+C, S, B)`, splitting-invariant → the G-2 SI≡Krylov gate
holds). The generic `SourceIteration`/`KrylovAcceleration` primitives stay problem-agnostic
(they see a resolvent + gains; the octant orchestration lives in the SN resolvent, never the
primitive).

#### DONE — additive foundation (no live-path change, all green)
- **Gate 1 (RED verification spec)** ✅ `926c4ca` — `tests/sn/verification/analytical/
  test_si_convergence_rate.py`: 3 `xfail(strict=True)` before/after rate gates + analytic
  `ρ_J=c` anchor + Krylov lower-bound bracket + correctness guards (G-1 k_inf=1.875 closed-form,
  G-2 SI≡Krylov fixed point, G-3 flat-flux foundation, G-4 vacuum negative control). 7 passed +
  3 xfailed. Baselines re-confirmed at HEAD 7d85222 (the `L−(S+B)`→`(L−S)−B` FP reassociation
  shifted slab counts +1): **slab-2g 655 / slab-4g 523 / box-2g 697 / vacuum 128**. (Fixed two
  latent skeleton bugs: the vacuum control reused a reflective mesh with an IGNORED
  `boundary_condition=` kwarg; `Solution.compare` is same-mesh-only → compare `scalar_flux.values`.)
- **3-1 measurement seam** ✅ `956f069` — `IterationHistory.total_inner_iterations` surfaced on
  the EIGENVALUE path (`SNSolver` accumulates `len(residuals)` per outer step in
  `_solve_source_iteration`/`_solve_krylov`; `solve_sn` reads `solver._total_inner_iterations`).
  SN-LOCAL, measurement-only: NO Protocol / `power_iteration`-signature change → ZERO 5-family
  blast. Fixed-source sets it `= n_inner`. Measured eigenvalue Jacobi baselines: **SI
  total_inner=371 / Krylov=310** (n_outer=3 both — outer count is splitting-INVARIANT; the inner
  SI count is where G-S shows → anchors the FUTURE eigenvalue rate gate `n_GS < 0.75×371`). keff
  unchanged. Flipped the eigenvalue forward-ref → passing foundation seam gate
  `test_eigenvalue_path_surfaces_total_inner_iterations`.
- **3-2 face-restricted reflect** ✅ `7d88ab0` — optional `faces=` on
  `SNBoundaryOperator.reflect_into_inflow`/`_reflect_trace` (block-diagonal over faces ⟹ EXACT
  restriction; unknown-face raises; backward-compatible). 21 boundary-op foundation tests.
- **3-3a polymorphic schedule** ✅ `514ae21` — `orpheus/sn/sweep_schedule.py`:
  `SweepSchedule.jacobi`/`.gauss_seidel`, `OctantSweep` (label + ordinate-index tuple),
  `OctantSweepGroup` (`sweeps` + `reflect_faces`). Mesh-time, flux-independent. 8 foundation
  tests (`tests/sn/sweep/core/test_sweep_schedule.py`). ⚠ **Finding for the 3b/3c tests:** the
  2-D test quad `Quadrature.product(n_mu=2, n_phi=4)` is AXIS-ALIGNED — its 4 in-plane octants
  are single-face `(±1,0)/(0,±1)` (each outflows ONE face), NOT `(±,±)` diagonals (diagonals
  need a level-symmetric cubature). Per-axis reflective 2-cycle `(−1,0)`↔`(+1,0)` on xmin/xmax;
  lexicographic order (`−1` before `+1`) ⟹ partial one-pass G-S.

#### NEXT — 3b + 3c (the fixed-point-changing CORE; hot-path; user said "proceed solo")
**READ FIRST:** carve substrate `.claude/agent-memory/explorer/si_gs_substep3_carve_substrate.md`
(exact anchors at HEAD 7d85222 — re-confirm each with `awk`/`grep`/`inspect.getsource` per L22).

**3b — `sweep_octant_group` + interior-`WavefrontFlux`-cache persistence (the biggest lift).**
`_sweep_2d_wavefront` (`sweep.py:767-944`) today: `wavefront = WavefrontFlux.zeros_on(sn_mesh)`
+ `psi_x=wavefront.face(0)`/`psi_y=wavefront.face(1)` + `wavefront.seed(boundary_flux)` (ι_*,
858-864); `for octant in quad.octants:` (873); per-octant body on `oct_idx` slices (911-933) —
`sweep_graph.apply(...)` mutates `psi_x_oct`/`psi_y_oct`/`angular_flux_oct` + accumulates the
SHARED `scalar_flux`, scatter-out (931-933); `wavefront.absorb(boundary_flux)` writeback (ι*, 942).
  1. Carve the per-octant body (911-933) into `sweep_octant_group(octant, *, wavefront, Q, sig_t,
     str_x, str_y, weights, angular_flux, scalar_flux, sn_mesh)` taking a PERSISTENT `wavefront`.
  2. Rewrite `_sweep_2d_wavefront` to: `wavefront=zeros_on; seed; for octant in quad.octants:
     sweep_octant_group(...); absorb` — i.e. the Jacobi schedule. **REFACTOR GATE: this must be
     BIT-IDENTICAL to today** (same octant order, one final absorb) — the rate baselines
     655/523/697/128 must NOT move — before any schedule wires in.
  3. Add a face-restricted `WavefrontFlux.absorb(boundary, faces=…)` (mirror sub-step 2's
     `reflect_into_inflow(faces=…)`) OR use `wavefront.edge_view(face)` per outgoing face — so a
     per-group absorb writes ONLY the group's outgoing faces' outflow into the boundary buffer
     before the reflect.
  4. Persistence: the single `wavefront` object (its `psi_x`/`psi_y` interior buffers) lives
     across the per-group `sweep_octant_group` calls within ONE resolvent `.solve`. It is ONE
     typed object's lifetime to extend. A bug here changes the FIXED POINT (caught by G-1/G-2/G-3).
  ⚠ L16: keep the per-cell work vectorised (the diagonal advance over `n_diag`); the per-group
  loop must NOT add per-cell Python. Profile before claiming.

**3c — the scheduled SN resolvent + wiring + flip gates.**
  1. A NEW SN scheduled resolvent (resolvent-layer G-S) built from `(L+C, B, SweepSchedule)`; its
     `.solve(rhs, initial_guess)` runs the §design seed-then-overwrite loop (seed `B·initial_guess`
     → per group: sweep its `sweeps` → absorb+reflect+re-seed its `reflect_faces` → final
     whole-trace absorb → return `TimedFullField(bulk, boundary)`). Resolvent path anchor:
     `InvertibleOperator.solve` → `_solve_timed_full_field` (`operator.py:1903-2017` — seeds from
     `rhs.boundary`, calls `transport_sweep`, returns the post-sweep boundary).
  2. Wire an `inner_schedule` selector (default `"gauss_seidel"`; `"jacobi"` control) through
     `solve_sn`/`solve_sn_fixed_source` → `_solve_fixed_source_si` (**2-D fixed-source FIRST**,
     the unambiguous target) then `_solve_source_iteration` (eigenvalue defers within the phase).
     SI gains drop to `(S,)`; `_within_group_krylov` UNTOUCHED.
  3. Flip the 3 RED rate gates: rewrite to measure BOTH schedules LIVE (`n(gauss_seidel) < 0.75 ·
     n(jacobi)`, both in-process) — retire the hardcoded 655/523/697 constants; Jacobi becomes a
     permanent control. Remove the `xfail(strict=True)`.
  4. Close-out: log **ERR-056** (`@catches` on the rate gate) + append **vv-principles Mode 9**
     to the skill (BOTH user-APPROVED) + archivist adds `:label:si-spectral-rate` to
     `docs/theory/discrete_ordinates.rst` so the `verifies("si-spectral-rate")` decorators resolve.

#### Verification (must hold)
- **Refactor gate (3b):** Jacobi schedule BIT-IDENTICAL to today's `_sweep_2d_wavefront` (the
  carve is mechanical; the baselines 655/523/697/128 must not move).
- **Rate gate (3c, the payoff):** G-S `n_inner < 0.75 × Jacobi n_inner` on the reflective configs,
  BOTH measured live (0.75 tolerates partial G-S on the cyclic box; ≥15% margin).
- **Correctness guards (PASS on BOTH schedules — fixed point is invariant):** G-1 k_inf=1.875
  (closed-form); G-2 SI≡Krylov converged flux rtol 1e-8 (Krylov = structurally-independent
  splitting); G-3 flat-flux balance (foundation); G-4 vacuum count UNCHANGED (negative control).
- **Converged flux bit-identical** Jacobi-vs-G-S (SI-rate-only — THE whole point).
- **L16 wall-clock gate** on the hot-path carve (no per-cell Python regression).

#### Gotchas / do-not-touch
- **KEEP `_reflect_outflow_into_inflow`** (`solver.py:973-1014`, whole-trace) — its ONE production
  caller (the eigenvalue reconstruction sweep `solve_sn:1139`) + ~7 test files stay on it; the
  G-S reflect is ADDITIVE (the face-restricted `reflect_into_inflow`), not a helper migration.
- **KEEP `_within_group_triple` as the `(L+C, S, B)` SSoT** — SI COMPOSES the scheduled resolvent
  from it; do NOT fork it. **Krylov stays byte-identical** (the G-2 gate; it has no schedule).
- **2-cycle ⟹ partial G-S; white/albedo → Jacobi; target = specular.** **1-D slab DEFERRED**
  (scan ≠ wavefront; out of scope here).

#### Pickup pointers
Substrate `.claude/agent-memory/explorer/si_gs_substep3_carve_substrate.md` (3b/3c anchors:
`_sweep_2d_wavefront` seed 858-864 / loop 873 / body 911-933 / absorb 942; `WavefrontFlux`
`seed`/`absorb`/`face`/`edge_view`; resolvent `_solve_timed_full_field` operator.py:1903-2017) ·
surface map `sn_si_reflective_gauss_seidel_recovery_surface.md` · spec
`si_gauss_seidel_rate_recovery_verification_spec.md` (all REFRESHED at HEAD 7d85222). Key files:
`orpheus/sn/sweep.py` (3b), `orpheus/sn/sweep_schedule.py` (3a ✅), `orpheus/transport/fields/
wavefront_flux.py` (seed/absorb/face/edge_view), `orpheus/sn/operator.py` (resolvent solve),
`orpheus/sn/solver.py` (3c wiring + `_within_group_triple` SSoT + the do-not-touch helper),
`orpheus/sn/boundary_operator.py` (`reflect_into_inflow(faces=)` ✅).

### Phase 4 — O.2b: adjoint + reciprocity + role-typing  *(floats; orthogonal to Phase 3)*

**Goal.** The G-adjoint + the reciprocity gate + completing the Source/Residual typing.

**Sub-steps.**
1. **Populate the metric:** `TraceSpace.inner_product_weights = |Ω·n| ⊙ w_n` (both
   factors); `_AdjointOperator.apply` already reads it.
2. **`StreamingOperator.apply_transpose`** — the reverse sweep (currently absent;
   "reserved for Phase H"). Then `.H = G⁻¹AᵀG` on every leaf.
3. **Gate 1.3 green** — flip the `xfail-strict` reciprocity gate for slab+sphere+cyl;
   add a NEW `|Ω·n|·w_n` boundary-block reciprocity test (today's gate is bulk-Euclidean
   → metric-blind → a wrong boundary adjoint is invisible) + an L11 negative control
   proving the weighting is load-bearing. 2-D adjoint stays xfail until a 2-D-adjoint sub-step.
4. **White-BC adjoint** — white is self-adjoint only under `|Ω·n|·w`; advertise
   `apply_transpose` for white.
5. **Composer role-derivation:** `OperatorSum`/`ScaledOperator`/`_AdjointOperator` must
   DERIVE `block_role` from operands (retire the hardcoded `InvertibleOperator.__init__`
   FULL stamp in the SAME change — twin-path); pin `(L).H is FullOperator`, `B.H is
   BoundaryOperator`; `realize_recursively` composed BCs derive (not stamp) their role.
6. **O.3b typed full source** (#201/#205): `BoundarySourceSink.from_spec` prescribed-inflow
   bridge; non-vacuum-BC MMS — new dataclass in `derivations/continuous/mms/sn.py` with the
   **anti-bias ansatz `(A(x)+μB(x))/W`** (c<1, A>0, derived `external_source`, prescribed
   inflow `γ₋ψ`) per vv-principles Mode-7.

**Verification.** Operator-algebra-core; Gate 1.3 green; L1 MMS (1-D + 2-D + curvilinear,
incl. the non-vacuum-BC case); the boundary-block reciprocity test + its negative control.
**Deferred pointer:** `#213` capabilities-as-morphism-class (`Iso ⊂ PartialIso ⊂ General`).

### Phase 5 — WavefrontFlux performance: angular-windowing → storage-B  *(#205)*

**5a — angular-windowing (the true `O(n²)→O(n)` 2-D win).** Reduce angular→scalar per
anti-diagonal in the SI iteration, never materializing the full `(N,ng,nx,ny)` angular
field. Consumer = the scalar-reducing 2-D SI iterate (exists post-Phase-1). The one-shot
final `Solution.angular_flux` reconstruction stays full; the Krylov path stays
full-angular (GMRES iterates the full bulk vector — windowing there would change the
domain). Win is in the SI ITERATION.

**5b — interior-face storage-B (the rolling moving-frontier window).** Swap the
`WavefrontFlux` backing from the full per-axis face-field to a rolling 2-diagonal window
(`O(N·ng·nx·ny) → O(N·ng·(nx+ny))`). **A 3-LAYER carve, not a backing swap:** the walk
indexes `psi_x`/`psi_y` with GLOBAL per-level face indices, so a window can't return a
full `(N,ng,nx+1,ny)` view → it touches the `SweepCellSlice` index protocol + the
`DiamondDifference` gather/scatter + the `WavefrontFlux.face()` API, plus a "capture
outflow boundary edges as the frontier sheds them" step. **Lands independent of 5a** (the
faces-only ~3–4× win). 1-D is a no-op.

**Verification.** Converged solution bit-identical (only the live working-set shrinks);
pin the peak-memory drop; **L16 wall-clock gate** (the rolling buffer must NOT add per-cell
Python — the diagonal advance stays vectorised over `n_diag`); de-risk diagnostic first.
Preserve the §3a′ axis-parametric unit test (the 3-D-readiness precondition).

**Honest scope (do NOT oversell):** the WavefrontFlux *typing* (storage-A, landed) is a
representation/elegance win only. 5a is the asymptotic win (SI iteration); 5b is the
constant-factor faces-only win + the 3-D enabler.

**Source:** `wavefront_flux_foundation.md` (archived) + `wavefront_flux_carve_substrate.md`
(refreshed anchors).

### Phase 6 — N-D foundation  *(future session; see `nd_foundation.md`)*

The `d`-generic upwind-DAG walk (`Σ_axis idx = ℓ`) collapsing 1-D/2-D/3-D into one
`_sweep_wavefront`, validated by 3-D Cartesian (unify-after-two: 1-D+2-D are the two
instances). **The load-bearing design question (§2.3 seam):** the d=1 path is a
parallel-prefix SCAN (`ordinate_scan`, Blelloch, 2 scan calls/sweep), NOT a wavefront —
a naive `d`-generic forward-substitution walk is O(nx) *sequential* for d=1 = a
**regression**. HARD CONSTRAINT: the unified walk must keep (or recover via a
segmented/parallel scan) the d=1 cumprod realization. What unification harvests: typed
`WavefrontFlux` as 1-D's interior representation, one typed `TraceExchange`, one `C¹→C⁰`
averaging primitive. Scope: 3-D **Cartesian** only (3-D curvilinear + MoC out). Stays in
`nd_foundation.md` (its own pickup checklist); this spine's tail just points here.

---

## 4. Verification spine (must-stay-green every commit)

operator-algebra-core gate; L0 streaming-equilibrium (Gate 1.1); L1 MMS (slab + sphere
+ cylinder + 2-D Cartesian; non-vacuum-BC at Phase 4); Resolution-A / regression-snapshot
bit-exact (or documented principled-equivalence per vv-principles §bit-identity); the
SI≡Krylov twin (necessary-not-sufficient, anchored to closed-form k_inf); sentinel
(no `-O`, assertions fire). **NO `continuous_get`** (the #212 hang — use direct k_inf /
Q/Σt / MMS). Default mode `python -O`; env `PYTHONPATH=$PWD .venv/bin/python`. L16
wall-clock gate on any hot-path carve (Phases 2, 5). Each sub-step its own commit; review
sub-agent output with full session context before committing.

---

## 5. Task tracker

The executable breakdown lives in the session task tracker (one task per phase + the
key sub-steps). Knock out Phase 1 first; Phase 4 may interleave after Phase 2 if
Gate-1.3-green is wanted early. Orthogonal reds (#212/#209/#195/#149/#189/#183/#191 + the
#206/#196 cylinder red, which Phase 2 fixes) are tracked as GitHub issues, NOT phases.
