# SI Gauss–Seidel recovery — interleaved `−B` reflect schedule (Wave O, #208-adjacent)

**Branch:** `refactor/field-role-typing` (worktree). **Base HEAD:** `266fcf5`
(B.5.2 fully landed). **Mode:** main-agent DIRECT authorship, turn-by-turn user
steering (`feedback_no_method_implementer_for_surgical_carves`) — touches the
green SI drivers + the bare 2-D sweep; commit per phase.
**Standing exclusions from commits:** `.claude/agent-memory/`, `.mcp.json`,
hooks, `derivations/diagnostics/`.

**⭐ SUBSTRATE STATUS (2026-06-04):** the typed `WavefrontFlux` substrate this
plan was waiting on is **LANDED** (storage-A Phases 0–3, HEAD `888775c`+): the
2-D sweep + matvec interior cochain is `WavefrontFlux` with typed `seed`=ι_* /
`absorb`=ι* / `edge_view`. **Sequencing: this G-S work is picked up AFTER Phase 6
(storage-B, the moving-frontier window) per the user directive** — so re-confirm
the substrate (incl. whatever Phase 6 changed the `WavefrontFlux` backing to)
before starting. The `(octant×face)` reflective graph composes as a typed
`sweep octant → ι* absorb → −B reflect → ι_* seed next octant` on this substrate.

**Design inputs (read these):**
- test-architect: `.claude/agent-memory/test-architect/si_gauss_seidel_rate_recovery_verification_spec.md`
- explorer: `.claude/agent-memory/explorer/sn_si_reflective_gauss_seidel_recovery_surface.md`

---

## 0. The diagnosis (what we're recovering)

BC extraction (O.4a.2/O.4b) made the sweep BARE (reads `ψ.boundary.inflow` as
fixed-for-the-sweep) and moved the reflective coupling to the external sibling
`−B` (`SNBoundaryOperator`), applied **once per full sweep** → inter-sweep
**Jacobi** (B fully lagged in `N`). The retired `bc.apply`-in-sweep read the
live `psi_x` buffer, so later octants saw earlier octants' fresh reflected
outflow → intra-sweep **Gauss–Seidel**. Same fixed point, **slower SI rate**.
In splitting terms we want `M = (L+C−B_lower)` (forward-substitution over the
octant order) instead of `M = (L+C)` with all of B lagged. Krylov is UNAFFECTED
(GMRES sees the whole `(L+C−S−F−B)`; splitting-invariant).

**The win is essentially 2-D** (multiple coupled octants per axis). 1-D slab has
2 direction groups; its SI rate is scattering-dominated (ρ_J = c), so G-S on the
reflective coupling barely moves it — treat 1-D as a CONTROL (count unchanged),
not a beneficiary. **Empirical, not a priori** — the rate harness (Phase 1)
decides.

---

## 1. Scope + the key constraints (from the two design passes)

- **Target path = 2-D fixed-source SI** (`_solve_fixed_source_si`, solver.py:1341):
  LIVE, and its inner count `n_inner` IS surfaced → measurable/testable. The
  reflect is `_reflect_outflow_into_inflow(solver._boundary_flux, sn_mesh)`
  (solver.py:1424) run once immediately before `transport_sweep` (1425). THE
  FIRST + UNAMBIGUOUS TARGET.
- **Eigenvalue SI** (`_solve_source_iteration`:534 → `SourceIteration(LC,
  _scattering_with_boundary_op, …)`:653) couples via the `S+B` fold (B·ψ.outflow
  in `rhs.boundary`, Jacobi). Its inner count is **NOT surfaced** (`n_inner is
  None`; the outer keff count is splitting-invariant) → needs a measurement
  seam. **DEFER to a later phase** (after the seam + the fixed-source win).
- **Measurement seam (test-architect requirement):** add
  `IterationHistory.total_inner_iterations` so the eigenvalue-path SI rate
  becomes observable. Fixed-source already surfaces `n_inner`.
- **Sweep factorability (explorer Q3):** `transport_sweep`/`_sweep_2d_wavefront`
  (sweep.py:756-940) is monolithic (`for octant in quad.octants`:852, per-octant
  `sweep_graph.apply`:903-925, seed 4 faces up front:836-843, writeback after
  loop:934-937). Octant bodies ARE cleanly factorable (each reads frozen seed,
  touches only its `oct_idx` slots; `scalar_flux` is an order-free additive
  accumulator). **Biggest lift:** the ephemeral `psi_x`/`psi_y` interior edge
  cache dies at function exit — the G-S driver must keep it ALIVE across
  octant-group calls within one sweep and relocate the boundary writeback to
  after each group so the interleaved `−B` sees fresh outflow.
- **Cyclic caveat (explorer Q4):** specular pairing = `quad.reflection_index`
  (directional.py:365-387); octants = `quad.octants` (directional.py:413-431).
  2-D with BOTH xmin+xmax reflective ⇒ `(−,+)↔(+,+)` is a 2-cycle ⇒ a single
  forward pass cannot drop both edges (the opposite-face back-edge stays
  Jacobi-lagged). Full one-pass G-S only when ≤1 reflective face per axis
  (half-symmetric domain). White/albedo couples ALL ordinates on the face →
  degenerates to Jacobi. **The recovery targets specular reflective faces; the
  back-edge staying Jacobi is expected (the legacy was also partial G-S by its
  fixed order).**
- **Additive, not migrating (explorer Q6):** KEEP the whole-trace
  `_reflect_outflow_into_inflow` for the Jacobi path + the final reconstruction
  sweep (solver.py:1146). ADD a new octant/face-restricted reflect + schedule
  for the G-S path. Express via the canonical `SNBoundaryOperator` (face-
  restrictable, boundary_operator.py:140-200) so it survives the O.2 collapse
  that retires the helper. ~8 test files import the helper — do not touch it.

---

## 2. Architecture — G-S as a selectable splitting

The G-S schedule is a SELECTABLE iteration-splitting for the SI driver, NOT a
sweep change. The sweep primitive stays bare; `−B` stays single-sourced; the
driver orchestrates `(bare octant-group sweep) → (face-restricted −B reflect) →
(next group)`. This is exactly the driver-splitting flexibility O.2 will own —
the recovery is a down-payment on it. The schedule (the specular octant DAG +
its topological order) is a **mesh-time derived object** (mirror
`SweepDependencyGraph`).

Pieces:
1. `IterationHistory.total_inner_iterations` (measurement seam).
2. `sweep_octant_group(octant, psi_x, psi_y, scalar_acc, …)` — factor from
   `_sweep_2d_wavefront`:880-925; bare; reads current seed + persistent
   `psi_x`/`psi_y`; writes its octant's outflow to the shared edge cache.
3. The specular octant SCHEDULE (mesh-time): the **`(octant × face)` reflective
   graph** (cross-domain-attacker, frame-validated — see
   `.claude/agent-memory/cross-domain-attacker/field_role_typing_faceflux_frames.md`,
   Frame 2 = sparse triangular factorization). A reflective BC is a directed
   edge `(producer octant, face) → (consumer octant, face)` via the specular
   permutation `quad.reflection_index`. **Edges that respect the octant
   topological order (producer swept before consumer) are Gauss-Seidel-ELIGIBLE
   → fold into the inverted triangular factor `L⁻¹_oct`; cycle-forming edges
   MUST be lagged (Jacobi).** "Which back-edges are in the factor" reduces to
   ONE topological comparison on this graph — no per-BC branching. (2-D
   both-faces-reflective = a 2-cycle → opposite-face edge stays Jacobi, §1.)
   **Literature home: KBA / Adams-Larsen / Pautz parallel-SN sweep scheduling —
   our `i+j=k` levels ARE the KBA wavefront diagonals**; this graph is exactly
   what KBA pipelining / Pautz optimal-scheduling operate on.
4. A face-restricted reflect (`SNBoundaryOperator` applied to one face's
   outflow→inflow) interleaved between octant groups by the schedule.
5. A G-S sweep driver that allocates `psi_x`/`psi_y` once, runs groups in
   schedule order with reflects between, accumulates `scalar_flux`.
6. Wire it as a SELECTABLE path in `_solve_fixed_source_si` (Phase 2); extend to
   `_solve_source_iteration` (Phase 3).

**Substrate note (build order):** the attacker confirmed the bare-sweep +
external `−B` is **pure-Jacobi by construction** (all reflective coupling
lagged); the G-S end is reachable by folding the order-respecting back-edges.
If the **WavefrontFlux** foundation (`.claude/plans/wavefront_flux_foundation.md`)
lands FIRST (the chosen order), this schedule becomes an explicit TYPED
composition — `sweep octant → ι* absorb → −B reflect → ι_* seed next octant` —
on `WavefrontFlux`/`BoundaryFlux`, instead of a raw-buffer-timing dance. Same
construction either way; typed substrate makes an illegal schedule
unrepresentable (Pattern 4).

---

## 3. Phasing (de-risk first; turn-by-turn)

### Phase 0 — DE-RISK (diagnostic, excluded)
A throwaway in `derivations/diagnostics/`: prototype the octant-interleaved
reflect on a small 2-D reflective box; confirm (a) same converged fixed point as
the current Jacobi path (bit-level on the converged flux), (b) fewer SI sweeps to
converge. STOP if the fixed point moves or the count doesn't drop. Cheap proof
before touching production.

### Phase 1 — measurement seam + rate harness (RED gate)
- Add `IterationHistory.total_inner_iterations`.
- Land the test-architect skeleton `tests/sn/verification/analytical/test_si_convergence_rate.py`:
  the `xfail(strict=True)` rate gate (RED on today's Jacobi — documents the
  regression), the correctness guards (k_inf=1.875, SI≡Krylov fixed point,
  flat-balance, vacuum negative control), config = Mixture B fully-reflective
  (slab control + 2-D box beneficiary), ε=1e-8, ratio/inequality pins (≥15%
  margin, never exact counts). Establish the Jacobi baseline empirically.
- **Gate:** the rate test is RED (xfail-strict), the correctness guards GREEN.

### Phase 2 — the G-S recovery on 2-D fixed-source SI (flip the gate GREEN)
- `sweep_octant_group` factor + the specular schedule + the interleaved
  face-restricted `−B` + the G-S sweep driver; wire as the selectable path in
  `_solve_fixed_source_si`. Keep the whole-trace helper for the Jacobi path +
  reconstruction.
- **Gate:** the rate test flips GREEN (`n_inner ≤ 0.75 × Jacobi_baseline` on the
  2-D box; 1-D slab + vacuum unchanged), correctness guards STILL green
  (fixed point unchanged — k_inf=1.875, SI≡Krylov, flat-balance), the octant
  equivalence snapshots (`test_2d_octant_sweep_equivalence.py`) + sentinel hold.
  elegance-enforcer review.

### Phase 3 — extend to eigenvalue SI (optional, after the seam)
- Apply the same schedule to the `S+B`-fold path in `_solve_source_iteration`;
  use `total_inner_iterations` to flip the eigenvalue-rate row (currently
  `xfail strict=False`) green.

### Phase 4 — docs + lessons
- archivist: the G-S/Jacobi splitting section in `docs/theory/discrete_ordinates.rst`
  (the L7-trap section already discusses traversal-order G-S, :1670-1800) +
  `operator_algebra.rst` (the splitting-as-selectable-`M` framing); the
  `verifies("si-spectral-rate")` theory label.
- vv-principles **Mode 9** (value-preserving solver-quality regression on an
  UNMEASURED cost axis — this bug's class; test-architect parked it pending
  approval) → ERR-056 when this lands. **APPROVE.**
- Optional: literature-researcher for the ρ_GS = c² model-problem factor
  (Adams & Larsen 2002; Hageman & Young consistently-ordered SOR) to tighten the
  rate gate 0.75 → ~0.55.

---

## 4. Verification (must-stay-green) + references
- Rate: `test_si_convergence_rate.py` (the before/after gate + Krylov lower
  bracket). Correctness: k_inf=1.875 (closed-form), SI≡Krylov fixed point
  (structurally-independent splitting cross-check), flat-balance (foundation),
  vacuum negative control. Snapshots: `test_2d_octant_sweep_equivalence.py`
  (octant order/G-S; case-3 `@catches("ERR-003")`). Sentinel 36/36 (no `-O`).
- Default `python -O`; env `PYTHONPATH=$PWD .venv/bin/python`. NO `continuous_get`
  (#212). Bash: `> file 2>&1; echo exit=$?`.
- PRE-EXISTING REDS (orthogonal): cylinder matvec #206/#196; #212 hang;
  `test_krylov_curvilinear_precond_safety`/`test_b1pp`/`test_krylov_restart`
  (#200-adjacent Krylov budget — a DIFFERENT convergence axis, do not conflate
  with this SI-rate work).

## 5. Risk ranking
1. **`psi_x`/`psi_y` lifetime across octant-group calls (HIGH/structural).** The
   factor must keep the interior edge cache alive + relocate writeback; a bug
   here changes the FIXED POINT (caught by the correctness guards, not just the
   rate gate). Phase 0 de-risk targets exactly this.
2. **Cyclic back-edge (MEDIUM).** 2-D both-reflective can't be full G-S in one
   pass; the gate's 0.75 factor must tolerate partial G-S (the test-architect
   set it conservatively). Don't over-promise full c² speedup on the cyclic box.
3. **Schedule correctness (MEDIUM).** Wrong specular pairing → wrong inflow →
   wrong fixed point (correctness guards catch it). Derive from
   `quad.reflection_index`, verify against the legacy octant order.
4. **Eigenvalue-path drift (LOW, deferred).** Phase 3; the seam + the fold's
   reflect interleave; gated behind the fixed-source win.
