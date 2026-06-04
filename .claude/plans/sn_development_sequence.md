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
- **Branch:** `refactor/field-role-typing` (pushed to `origin`).
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
1. `MEMORY.md` → the `[[project-wave-o-operator-algebra]]` ⭐⭐ LATEST block (the
   reshaped dependency chain + 2-D SI Phase A landing).
2. THIS file (`sn_development_sequence.md`) — the 6-phase sequence (§3).
3. The session **task tracker** (`TaskList`): tasks #17–#22 = the 6 phases, chained.
   **If the tracker did not survive the session boundary, recreate it from §3 + §5**
   (one task per phase; chain blockedBy 17→18→{19,20}, 19→21, 21→22).
4. Phase 1 (#17) is the only unblocked phase — start there (test-architect FIRST).

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

**Standing reds (ORTHOGONAL — do NOT block this sequence; tracked as issues):**
cylinder matvec #206/#196 (curvilinear WDD O(h) SI-vs-Krylov asymmetry — *fixed by
Phase 2*); #212 `continuous_get` hang (reference-registry build; has a proposed
patch `derivations/diagnostics/sn_keff_hang_PROPOSED_fix.patch`; deselect the 3
affected eigenvalue tests meanwhile); #209 `ordinate_scan` NaN on a zero a-chain
(cylindrical pole); #195 curvilinear MMS magnitude at nx=160. None gate the spine.

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
| R5 | deprecated `power_iteration` (5 solvers) vs unused `KEigenvalue` + ~8 inline test copies | `numerics/eigenvalue.py` vs `numerics/iteration.py` | Phase 1 |
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

### Phase 1 — Unify the solve architecture  *(R1 + R2 + R5 + gap)*

**Goal.** ONE within-group solve primitive, selectable SI/Krylov, consumed by ONE
`KEigenvalue` outer AND the standalone fixed-source entry. After this: exactly one
source-iteration path, one Krylov path, one power iteration.

**Sub-steps.**
1. **Migrate `_solve_fixed_source_si` onto `SourceIteration`** (R1) — retire its bespoke
   `for n_inner` loop; build the same `SourceIteration(LC, S+B, F=0)` triple the
   eigenvalue inner uses, differing only in `q_ext` (external Q vs fission source) and
   the returned scalar contract. (Keep `_reflect_outflow_into_inflow` available — Phase 3
   needs the additive octant-restricted variant; the whole-trace helper stays for the
   final reconstruction sweep.)
2. **One Krylov builder** (R2) — factor the byte-identical `KrylovAcceleration(LC, S+B,
   F=0, precond=identity, restart=N·ng·nx·ny)` construction shared by `_solve_krylov`
   and `_solve_fixed_source_krylov` into one site parametrised by `q_ext`.
3. **Un-gate fixed-source-Krylov-2D** (gap) — same recipe as 2-D SI Phase A: de-risk
   (it's the twin of the working eigenvalue Krylov), add an SI≡Krylov 2-D fixed-source
   equivalence pin, delete the `solver.py:~1582` guard.
4. **One power iteration** (R5) — migrate `solve_sn` (and the 4 sibling solvers: cp /
   diffusion / homogeneous / moc) off the deprecated `numerics.eigenvalue.power_iteration`
   onto `numerics.iteration.KEigenvalue`; retire `power_iteration` + `EigenvalueSolver`;
   rewire the ~8 inline `for _ in range(N): compute_fission_source → solve_fixed_source
   → compute_keff` loops in `test_keff_2d.py` to the primitive (retirement = test migration).

**Verification.** Bit-identical / principled-equivalent k_eff + flux on the existing
regression snapshots (all geometries); SI≡Krylov on the 2-D fixed-source case; the
`KEigenvalue` migration must reproduce `power_iteration`'s eigenvalue trajectory
(pin the keff history). Sentinel (no `-O`). Each sibling-solver migration its own commit.

**Gotchas.** The 5 cross-solver `power_iteration` callers make R5 the widest blast
radius — do the dependency audit (L20) before deleting. The fixed-source-Krylov un-gate
shares #200's `preconditioner=identity` stopgap — keep identity here, fix in Phase 2.

**Proactive:** test-architect FIRST (carve crosses eigenvalue↔fixed-source and the
deprecated↔new outer boundary; L17 crosswalk).

### Phase 2 — O.2a: the honest `L+C−S−F−B` driver  *(R3 + R4 + #206 + #196 + #200)*

**Goal.** The iteration drivers consume the WHOLE loss operator as one composed
block operator, retiring the transitional folds and collapsing the matvec/sweep twins.

**Sub-steps.**
1. **Give `ScatteringOperator` a domain** — the `S+B` fold (`_scattering_with_boundary_op`)
   currently type-checks ONLY because `S.domain is None` (the forcing tripwire). Setting
   the domain makes the fold throw, forcing the honest composition.
2. **Compose `L+C−S−F−B`** on the direct-sum carrier; the driver inverts/applies it
   directly. Retires the `S+B` fold AND collapses the two `−B` routes (R3): the driver
   path no longer needs `_reflect_outflow_into_inflow` (Phase 3 will have added the
   octant-restricted variant; the whole-trace helper survives only for the final
   reconstruction sweep).
3. **Unify the 3 `(L+C)` impls** (R4 = #206): `_compute_LpC` (production) /
   `_compute_decomposition` (dual-emission introspection) / the solve-sweep into one
   geometry-agnostic action (the 2-D `_apply_2d_cartesian` already shares the DD closure
   via `graph.residual` — extend that shape to 1-D). **This fixes #196** (the curvilinear
   SI-vs-Krylov O(h) asymmetry is the symptom of the twin) by construction.
4. **Trace-only `B.reflect_into_inflow(boundary)`** entry on `SNBoundaryOperator` —
   retires the zero-bulk-probe shim + its ~6 call sites.
5. **`#200` real preconditioner** — replace `preconditioner=lambda q: q` (both R2 sites,
   now one) with the block-inverse / sweep-as-preconditioner (the silent-fallback bug is
   gone post-Phase-1.2; safe to re-enable).

**Verification.** Curvilinear MMS SI-vs-Krylov flux-shape now O(h²) both paths (the
#196 canary `test_phase_e_trajectory_resolvent_flux_shape_crosscheck` flips); vacuum
bit-identical; reflective convergence-equivalence; the 3-impl unification is bit-identical
or principled-equivalent (FP-reduction-order, documented). Operator-algebra-core gate.

**Gotchas.** Operation-order discipline (vv-principles bit-identity) on the unified
`(L+C)` action. The composed-operator restructure touches `SourceIteration`/
`KrylovAcceleration` signatures — keep Phase 1's single-primitive shape.

### Phase 3 — SI Gauss-Seidel recovery  *(the motivating feature)*

**Goal.** Recover the intra-sweep Gauss-Seidel reflective coupling lost to BC
extraction. Bare sweep + external `−B` is inter-sweep **Jacobi** (B fully lagged in N);
legacy intra-sweep `bc.apply` was **G-S** (later octants saw earlier octants' fresh
reflected outflow). Same fixed point, slower SI rate. **Krylov is splitting-invariant
(UNAFFECTED).**

**Mechanism.** The **`(octant × face)` reflective graph** (mesh-time): edges
`(producer octant, face) → (consumer octant, face)` via the specular permutation
`quad.reflection_index`. Order-respecting edges → G-S-eligible (fold into the inverted
triangular factor `L⁻¹_oct`); cycle-forming edges → lagged (Jacobi). **Literature home:
KBA / Adams-Larsen / Pautz** — the `i+j=k` wavefront levels ARE the KBA diagonals.

**Caveats (load-bearing).** 2-D with BOTH xmin+xmax reflective ⟹ `(−,+)↔(+,+)` is a
2-cycle ⟹ one forward pass can't drop both back-edges (opposite-face edge stays Jacobi).
Full one-pass G-S only when ≤1 reflective face per axis. White/albedo couples all
ordinates → degenerates to Jacobi. Target = **specular reflective** faces. **1-D slab is
a CONTROL** (scattering-dominated ρ_J = c; little to recover).

**Sub-steps.**
1. **Measurement seam:** add `IterationHistory.total_inner_iterations` (the eigenvalue-path
   SI rate is currently unobservable — `n_inner is None`; fixed-source already surfaces it).
2. **Octant/face-restricted reflect:** ADD a new restricted-reflect via the canonical
   `SNBoundaryOperator` (face-restrictable) — **additive**, NOT a migration; the whole-trace
   `_reflect_outflow_into_inflow` stays for Jacobi + final reconstruction (~8 test files
   import it — don't touch).
3. **Keep the interior edge cache alive across octant-groups** within one sweep (the
   biggest structural lift): the ephemeral `WavefrontFlux` `psi_x`/`psi_y` must persist
   across octant-group calls and the boundary writeback relocate to after each group, so
   the interleaved `−B` sees fresh outflow. A bug here changes the FIXED POINT (caught by
   correctness guards).
4. **Target order:** 2-D fixed-source SI (the unambiguous target) FIRST; eigenvalue SI
   defers within the phase.

**Verification.** Phase-1 RED gate first: `xfail(strict=True)` rate gate (Mixture B
fully-reflective 2-D box, ε=1e-8, ratio pins ≥15% margin). Phase-2 flips green:
`n_inner ≤ 0.75 × Jacobi_baseline` (the 0.75 tolerates partial G-S on the cyclic box);
1-D + vacuum unchanged; **converged flux bit-identical** (SI-rate-only — this is the
whole point). vv-principles **Mode 9** (value-preserving solver-quality regression on an
UNMEASURED cost axis — this bug's own class) → log **ERR-056** when it lands (user
APPROVED). Optional literature-researcher for the `ρ_GS = c²` model-problem factor.

**Source:** `si_gauss_seidel_recovery.md` (archived) + `sn_si_reflective_gauss_seidel_recovery_surface.md`
(refreshed anchors; NOTE its Q1 "bare scattering_op at solver.py:571" drift-flag is
STALE — all 3 inner drivers use the `S+B` fold).

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
