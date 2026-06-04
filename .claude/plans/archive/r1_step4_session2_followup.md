# R-1 Step 4 — Session 2 follow-up plan

**Branch**: `refactor/sn-operator-algebra` (continues from HEAD `d1c1c9e`)
**Origin**: session 1 retrospective (commits A–F landed; G + H deferred; 8 process/architecture issues surfaced).

.. note::

   **2026-05-21 — IN FLIGHT**.  Phase 0 landed (commits ``58ce7c6``,
   ``a05edcd``, ``ff6c57b``).  Phase 1.1 is mid-flight: Groups A, B, C done
   uncommitted on top of ``ff6c57b``; Group D ~90% done; Group E pending.
   The IN-FLIGHT state, all decisions made, and the resuming-fresh
   checklist live in
   ``.claude/plans/r1_step4_session2_state.md`` — read that FIRST when
   picking this up.  The user direction during Phase 1.1 expanded the
   scope (singledispatch on producers + sweep API change + legacy
   migration); see the state doc for the architecture decisions log.

   GitHub issues filed this session:

   * **#202** — A1 producer normalisation (closed by Phase 1.1 commit `de8822d`).
   * **#203** — A2 KrylovAcceleration ``CAP_STATELESS_INVERSE`` (CLOSED — superseded by Phase 1.2 structural fix in commit `c93355c`; see §1.2 below).
   * **#204** — A5 promote diagnostics to permanent L1 (Phase 1.3 — next).
   * **#205** — Cross-method field architecture (deferred; the dimensional
     sin is acknowledged in code comments).
   * **#206** — WDD recurrence machinery unification (deferred).

## Why this plan exists

Session 1 shipped a working carve but with ~3× the debugging cost a principled session would have needed. The retrospective surfaced eight failure modes, ranked by leverage. Without containment, the same patterns will resurface in Step G (the deferred retirement + sister carve) AND in every future operator-algebra refactor.

Containment must happen BEFORE Step G touches code, because:

- Step G repeats the carve pattern (G1 / G2 = same shape as D / E).
- Step G3 retirement compounds whatever convention drift is already in place.
- Step H Sphinx narrative is wrong by construction if Pattern 7 violations remain on the production path.

This plan sequences the fixes into 7 phases. **Phases 0–2 are process + architectural cleanup that DOES NOT touch the deferred carve.** Phase 3 onwards is the deferred work, now standing on a cleaner foundation.

Each phase has a single commit, an acceptance gate, and explicit reversibility. The commits are independently reviewable — if any phase reveals a deeper issue we land what's principled and pause.

---

## Phase 0 — Capture the lessons in durable artifacts

**Goal**: prevent the eight failure modes from recurring on the NEXT operator-algebra carve (Step G1/G2 directly, MoC and CP later).

### 0.1 Update `.claude/lessons.md` with the four session-1 lessons

Add to `.claude/lessons.md`:

- **L17 — Convention crosswalk before carve**: a carve that crosses subsystem boundaries (operator algebra ↔ sweep, scalar ↔ per-ordinate, packed ↔ typed) MUST start with a written crosswalk table enumerating each subsystem's input/internal/output convention. The crosswalk is the architecture; the code is its transcription. Skipping it costs ~3× debug time (session 1 evidence: three convention bugs each took 1+ hour to diagnose; the crosswalk would have been ~15 minutes).
- **L18 — Pattern 7 at the producer, not the consumer**: when a producer (e.g. `ScatteringOperator.apply`) outputs a value in convention A and consumers (e.g. `InvertibleOperator.solve`) expect convention B, fix the producer. Consumer-side bridges multiply with every new consumer; producer-side normalisation costs once. Session 1's `* sum_w` bridge in `InvertibleOperator.solve` is the textbook example to retire.
- **L19 — `None` defaults that depend on unstated invariants are dangerous**: `KrylovAcceleration(preconditioner=None)` silently routed through `L.solve` assuming statelessness. `L.solve` was NOT stateless (read `rhs(1)` history). The silent fallback was the bug. Default values for behavioural parameters MUST either advertise their preconditions in the type system OR require explicit caller choice. Same lesson applies in MoC, CP, kinetics whenever a primitive auto-selects behaviour based on another primitive's capability advertisement.
- **L20 — Retirement requires an upstream dependency audit**: a retirement step that says "retire X" without enumerating "who calls X" is incomplete. Session 1's Step G under-estimated scope by ~2× because the plan didn't audit `solve_sn_fixed_source`'s callgraph through the legacy symbols. The audit is `grep -rn "<symbol>"` plus a callgraph walk for each retirement target. Cost: 10 minutes per retirement target. Skipped cost: full re-plan mid-session.

**Acceptance**: four lessons added to `.claude/lessons.md`, each with the failure mode + the structural fix + a generalisation pointer (where this lesson applies beyond SN).

**Commit**: `docs(lessons): R-1 Step 4 session-1 lessons — convention crosswalk, Pattern 7, None-default, retirement audit`

### 0.2 Update `coding-elegance` skill with the carve crosswalk template

Add a new section "Convention crosswalk template" to `.claude/skills/coding-elegance/SKILL.md` under Pattern 7. Template:

```
## Convention crosswalk template (use before any operator-algebra carve)

For each subsystem you'll cross:

| Subsystem | Input convention | Internal convention | Output convention |
|-----------|------------------|---------------------|-------------------|
| Producer A| ...              | ...                 | ...               |
| Consumer B| ...              | ...                 | ...               |
| Bridge    | (which way)      | (which transform)   | (which value)     |

The "Bridge" row is where Pattern 7 demands action: either move the
bridge to the producer (definition site), or document why it must
stay at the consumer (rare; load-bearing reason required).

Apply to: per-ordinate vs iso scalar, /W normalisation,
sign conventions (μ axis direction), packed vs typed layout.
```

**Acceptance**: skill file updated; main-agent + sub-agents that preload `coding-elegance` (method-implementer, qa, test-architect, numerics-investigator) see it on dispatch.

**Commit**: `docs(skills): coding-elegance — convention crosswalk template (R-1 Step 4 lesson)`

### 0.3 Update `subagent-handoff-protocol` skill with the proactive-dispatch trigger

Add: **before** any operator-algebra carve that crosses subsystem boundaries (defined as: typed ↔ packed, per-ordinate ↔ iso, normal ↔ adjoint, scalar ↔ angular), the main agent MUST dispatch `test-architect` to design the verification plan. The dispatch is REACTIVE today (after bugs); making it PROACTIVE prevents the bug.

**Acceptance**: skill file updated with the new trigger condition explicit in the agent table.

**Commit**: `docs(skills): subagent-handoff-protocol — proactive test-architect trigger for operator-algebra carves`

---

## Phase 1 — Architectural cleanup (A1, A2, A5)

**Goal**: retire the two technical debts session 1 created so Step G stands on a clean foundation.

### 1.1 A1 — `ScatteringOperator.apply` typed branch applies `/sum_w` at the producer

**File**: `orpheus/sn/scattering.py`.

Change the typed-AngularFlux branch of `ScatteringOperator.apply` (~line 836):

```python
# Before:
combined: PerOrdinateSource = iso + aniso
return AngularFlux(combined.values, mesh)

# After:
# R-1 follow-up (A1) — per-ordinate normalisation lives at the
# producer: iso source rescaled by 1/sum_w so the result is already
# the per-ordinate value the operator-algebra consumer expects.
sum_w = float(mesh.quad.weights.sum())
combined: PerOrdinateSource = iso / sum_w + aniso
return AngularFlux(combined.values, mesh)
```

Then drop the corresponding `* sum_w` bridge from `InvertibleOperator.solve` and the `/ sum_w` rescaling from `_solve_krylov` / `_solve_source_iteration` in `solver.py`. These three rescaling sites all dissolve.

**Migration of `test_S_typed_matches_raw`**: this test pinned `S.apply(typed).values == S.apply(bare)`. After A1 the typed branch is `/sum_w`-already-applied; the bare branch is not. Rewrite the test to pin the NEW contract: `S.apply(typed).values * sum_w == S.apply(bare)` (the consistent algebraic relationship).

**Verify**:
- `tests/sn/test_scattering_operator.py` — all pass (the `foldable_part + residual_part = full` identity is unaffected).
- `tests/sn/test_operators_apply_typed.py` — update `test_S_typed_matches_raw`.
- `tests/sn/l1_analytical/test_kinf_homogeneous.py` — 35 pass + 5 xfailed (unchanged).
- `tests/sn/test_invertible_operator.py` — 18 pass (unchanged).
- New L0 test: `tests/sn/test_scattering_operator.py::TestProducerSideNormalisation::test_typed_apply_returns_per_ordinate_already_normalised`.

**Acceptance**: the three rescaling sites collapse to ONE producer-side normalisation. Sphinx narrative for `ScatteringOperator` documents the convention; `coding-elegance` Pattern 7 example list updated to cite this as the concrete fix.

**Commit**: `refactor(sn): ScatteringOperator typed apply normalises iso at producer (Pattern 7; retires consumer-side /sum_w bridge)`

**Reversibility**: revert is mechanical (paste the bridge back in three sites; revert the producer change).

### 1.2 SHIPPED — sweep+matvec unification through M-M's psi_half_seed (#203 superseded)

**Status**: COMPLETE — commit `c93355c` (5 files, net **−35 LOC**).
**Closes**: #203 (A2 — superseded by the structural fix).
**Reframing**: the A2 plan (CAP_STATELESS_INVERSE capability advertisement)
was a SYMPTOM-FIX for the silent-fallback bug.  The architectural read on
2026-05-22 (user-led discussion) identified the deeper cause: the 1D sweep
was duplicating M-M's Carlson coupled-pole seed inline (`carlson_inward_sweep_from_source`
direct call + inline `Q_bar` derivation), while the matvec already routed
through `MorelMontryAngularSweep.precompute_psi_state → psi_half_seed`.
Two paths to the same recurrence — Pattern 2 (single source of truth)
violation and concept leakage.

**The unification (one strategy, two consumers)**:

* `_run_1d_sweep` (`orpheus/sn/sweep.py`) drops the inline Q_bar derivation
  AND the direct `carlson_inward_sweep_from_source` call.  Builds the
  per-level `CarlsonSweepContext` and calls
  `closure.psi_half_seed(psi_level, context)` — the SAME strategy the
  matvec invokes.  Per-level normalisation (`sum_w_level`) replaces the
  latent global-weights bug that affected cylindrical multi-level quadratures.

* `InvertibleOperator.solve(rhs, *, initial_guess=None)`
  (`orpheus/sn/operator.py`) takes the previous iterate as an explicit
  kwarg.  The Carlson seed travels through `initial_guess`, not through
  `rhs(1)`.  AngularFlux's history machinery (`history_depth`, `.stash`,
  `(1)` indexing) preserved unchanged — reserved for future time-derivative
  tracking per user direction.

* `SourceIteration._solve_with_seed` (`orpheus/numerics/iteration.py`)
  dispatches via `inspect.signature` introspection: typed
  `InvertibleOperator.solve` receives `initial_guess=psi_prev`; generic
  `MatrixOperator.solve(b)` (foundation L0 fixtures) receives positional
  rhs only.  `_attach_previous` helper retired.

* `apply_sweep_1d` retired (zero callers).

**The silent-fallback bug class evaporates structurally** — no
`CAP_STATELESS_INVERSE` capability needed.  `InvertibleOperator.solve` is
now a pure function of `(rhs, initial_guess, boundary)`.  GMRES residual
calls hit the preconditioner with `initial_guess=None` → cold-start
Carlson seed (deterministic, no garbage).

**Verification** (matches A1 baseline exactly):
- L1 analytical: 35 pass + 5 xfailed in 4:55
- Curvilinear streaming-equilibrium + dispatch + snstreaming: 60 pass + 1 skip in 15:31
- Regression: 5 failed (ALL PRE-EXISTING on baseline `ff6c57b` — slab-2g stale snapshots, sphere-P1-aniso ZeroDivision, 2D NotImpl)
- Iteration: 18/18 · Invertible operator: 19/19

**Key lesson** (captured in L21): when sweep and matvec are different
applications of the same operator, they should share ONE strategy.
Reduce strategies; don't add alternatives.  The Carlson seed is M-M's
internal business; the only difference between matvec and sweep is which
ψ they feed it (current target vs previous iterate).

### 1.3 A5 — Promote diagnostics to permanent L1 tests

**Files**:
- `derivations/diagnostics/diag_r1_step_d_probe_b_identity_precond.py` → `tests/sn/test_krylov_curvilinear_precond_safety.py`
- `derivations/diagnostics/diag_r1_step_e_invertible_solve_w_bridge.py` → fold into `tests/sn/test_invertible_operator.py`

For each: add `@pytest.mark.l1` (or `@pytest.mark.l0` for the InvertibleOperator unit tests), `@pytest.mark.verifies(...)` linking to the relevant Sphinx equation labels, and `@pytest.mark.catches("ERR-NNN")` after logging an ERR-catalog entry.

After A1 lands, the W-bridge test asserts the NEW producer-side convention (the bridge is gone; the test verifies its absence is correct).

**Acceptance**: 1.2s permanent L1 cross-checks replace the full L1 sweep as the fast gate. Full L1 sweep stays as the final gate before merge to main.

**Commit**: `test(sn): promote R-1 Step 4 diagnostics to permanent L1 (curvilinear-precond + W-bridge regression catches)`

### 1.4 Log ERR catalog entries

Update `.claude/skills/vv-principles/error_catalog.md` with:

- **ERR-NNN — Convention drift: per-ordinate vs iso magnitudes on the typed-AngularFlux operator algebra**. Failure mode 6 (convention drift). How it hid: legacy SI used IsotropicSource which auto-rescales; new typed algebra used per-ordinate everywhere; bridge was missing. Caught by: slab-2g eigenvalue keff=5.0 vs ref=1.875 (Step D) and slab-2g SI keff=1.484 vs 1.875 (Step E). Lesson: convention crosswalk before carve.
- **ERR-NNN — Silent preconditioner fallback breaks stateful-inverse contract**. Failure mode 4 (wrong recursion / state). How it hid: `preconditioner=None` defaulted to `L.solve` which read `rhs(1)` history; GMRES residual vectors have no history; sweep silently used in-iteration default. Caught by: sphere-2g-krylov 470× slowdown + keff oscillation.

**Acceptance**: error catalog has two new entries; corresponding pytest marks (`catches="ERR-NNN"`) cross-link.

**Commit**: `docs(vv): error catalog — ERR-NNN convention-drift and silent-precond-fallback (R-1 Step 4)`

---

## Phase 2 — Step G pre-work (dependency audit + crosswalk)

**Goal**: NEVER touch retirement code until the audit is written down.

### 2.1 Dependency audit for each retirement target — SHIPPED (2026-05-22)

**Status**: COMPLETE — `.claude/plans/r1_step4_g_dependency_audit.md`
(976 lines) is the authoritative artifact.  Dispatched the `explorer`
agent against `a638798`; audit covers all 12 listed symbols + the
augmentation questions (A: typed vs legacy `apply` equivalence, B:
2-D Cartesian retirement scope, C: `_with_traces` relocation cost, D:
GitHub issues against the scope).

**Key findings that amend §3 / §4 of this plan** (the audit
supersedes the original table here):

- **SURPRISE-1 — `EquationMap` is a KEEPER, not a retirement target.**
  Consumed by the `_with_traces` family (also keepers) AND by the
  typed `StreamingOperator._ensure_eq_map` / `CollisionOperator._ensure_eq_map`.
  Phase 5 (§5) relocates it alongside the `_with_traces` helpers
  into `orpheus/sn/angular_flux.py`.  **§4.4 below is wrong as
  originally written** — `EquationMap` does not retire.
- **SURPRISE-2 — `build_equation_map` (Cartesian) CANNOT retire in
  Step G.**  The typed 2-D Cartesian path in
  `StreamingOperator._ensure_eq_map` (`operator.py:1887`) and
  `CollisionOperator._ensure_eq_map` (`operator.py:2222`) consumes it.
  Defer Cartesian-`build_equation_map` retirement to Phase A
  alongside 2-D Cartesian absorption per #199 item 5.
- **SURPRISE-3 — `solution_to_angular_flux*` are consumed by
  `transport_operator_matvec_unified` helpers** (the unified matvec
  engine).  Retirement either rewires those helpers or preserves
  the decoders as `_with_traces` wrappers.
- **SURPRISE-4 + SURPRISE-5 — 2-D Cartesian fixed-source Krylov has
  no typed target post-G1; the principled answer is `_solve_fixed_source_si`
  (geometry-agnostic via `transport_sweep`).**  G1 should
  `NotImplementedError` on 2-D Krylov; G2's SI carve handles 2-D
  Cartesian naturally.
- **SURPRISE-6 — `derivations/diagnostics/diag_krylov_iter_breakdown.py:118`
  calls `_make_sweep_preconditioner` with a stale 3-arg signature**
  (`(eq_map, n_unknowns, sum_w)` vs production 2-arg).  Diagnostic
  is broken on current tree; will simply break harder at G3a.  No
  action required.
- **SURPRISE-7 — `tests/sn/test_phase_c_gates.py` has 11
  `SNStreamingOperator` instantiation sites NOT mentioned in Issue
  #199.**  These need migration / retirement / reframe before G3d
  can land.  **Expand #199 body before G3d starts**.
- **SURPRISE-8 — `docs/theory/discrete_ordinates.rst` has 10+
  `SNStreamingOperator` references** plus equation-map narrative
  sections.  Phase 6 / H Sphinx narrative is ~3-4 hours minimum.

### Legacy vs path-forward — user-directed scope correction (2026-05-22)

**SURPRISE-1's original conclusion was wrong** ("`EquationMap` is a
keeper").  Per user direction:

> "The legacy architecture needs to be retired and equation map is
> going with it.  Differentiate what is legacy architecture from
> the path forward.  The path forward is operator algebra, unified
> indices that do not require adapters."

* **Legacy architecture (RETIRES)**: packed-1D vector layout with
  `EquationMap` slot-dispatch — includes `EquationMap` class, all
  `build_equation_map*` factories, all `solution_to_angular_flux*`
  decoders, **AND** the `_with_traces` family (`pack_with_traces`,
  `solution_to_angular_flux_with_traces`,
  `build_equation_map_with_traces`), **AND** the typed operators'
  `_ensure_eq_map` machinery (the slot lookup is itself an
  adapter), **AND** `AngularFlux.to_flat_with_traces` /
  `from_flat_with_traces` (typed ↔ legacy bridges), **AND**
  `SNStreamingOperator`, **AND** `transport_operator_matvec_unified`
  in its current packed-face-slot signature.  See the audit's
  CORRECTION section for the full inventory.

* **Path forward**: typed `AngularFlux.values` (shape
  `(N, ng, nx, ny)`) + `BoundaryFlux.xmax_face` / `xmin_face`
  (shape `(N, ng)`) — natively indexed by ordinate + group; NO
  `eq_map.face_outer_ordinate` slot lookup.  Operators consume
  `AngularFlux` and emit `AngularFlux`; the face-inflow ordinate
  mask is derived from the quadrature, not from a precomputed
  slot map.  Phase 5 of this plan (relocate `_with_traces` family
  to `angular_flux.py`) is **hollowed out** — every former
  relocation target retires instead.

**Retirement sequence (canonical — supersedes the older §4.1-§4.4
plan AND the prior 8-step list)**:

1. **G0 (NEW pre-work)** — write a native-shape unified matvec that
   consumes `AngularFlux.values` + `BoundaryFlux.xmax_face` /
   `xmin_face` directly; derives inflow ordinates from
   `sn_mesh.quad`; returns cell + face residuals in matching
   native shapes.  Pin equivalence (bit-identity or principled-
   equivalent per `vv-principles`) to the current
   `transport_operator_matvec_unified` packed-face-slot path on
   slab + sphere + cylinder.
2. **G1** — carve `_solve_fixed_source_krylov` onto
   `KrylovAcceleration` + typed `AngularFlux` (1-D only;
   `NotImplementedError` on 2-D).  Orphans
   `_make_sweep_preconditioner` and `_build_rhs_*`.
3. **G2** — carve `_solve_fixed_source_si` onto `SourceIteration` +
   typed `AngularFlux` (geometry-agnostic — handles 2-D Cartesian;
   SURPRISE-5).
4. **G3a** — delete `_make_sweep_preconditioner` +
   `_build_rhs_{cartesian, spherical, cylindrical}` (mechanical;
   closes #174).
5. **G3b** — migrate `StreamingOperator._apply_typed` +
   `CollisionOperator._apply_typed` to the native matvec from G0;
   drop `_ensure_eq_map`, `_eq_map` cache, `n_unknowns`.
6. **G3c** — retire the legacy `transport_operator_matvec_unified`
   (packed-face-slot signature).
7. **G3d** — retire `pack_with_traces`,
   `solution_to_angular_flux_with_traces`,
   `build_equation_map_with_traces`, **`EquationMap`** (the user-
   corrected scope: keeper status revoked).
8. **G3e** — retire `solution_to_angular_flux_{spherical, cylindrical}`
   + `build_equation_map_{spherical, cylindrical}`.
9. **G3f** — retire `SNStreamingOperator` (closes #199 items 1-4).
   Migrate the 11 sites in `test_phase_c_gates.py` (SURPRISE-7) and
   the typed `(L+C)` apply comparison in
   `test_keigenvalue_matches_solve_sn_2g_slab`; delete
   `test_snstreamingoperator.py`; rewrite or retire monkey-patch
   fixtures in `test_l1_standoff_slab_cylinder.py` +
   `test_unified_matvec_*.py`.
10. **G3g** — retire `AngularFlux.to_flat_with_traces` +
    `AngularFlux.from_flat_with_traces` (no typed ↔ legacy bridge
    left to feed).
11. **G3h — DEFERRED to Phase A** — retire Cartesian
    `build_equation_map` + Cartesian-1-D fallback in
    `solution_to_angular_flux`.  Blocked by 2-D Cartesian
    absorption (#199 item 5).
12. **Phase 6 / H** — Sphinx narrative (~3-4 h minimum per
    SURPRISE-8).

**GitHub issue cross-references**:

- **#199** (Step G omnibus) — Step G ships items 1-4; defers item 5
  (2-D Cartesian absorption) to Phase A.  Body should be amended
  to mention `test_phase_c_gates.py` migration (SURPRISE-7) AND
  the G0 native-matvec pre-work.
- **#174** (`_build_rhs_cartesian` softer refactor) — subsumed by
  G3a; close with reference to the G3a commit when it lands.
- **#160** (originating issue for `SNStreamingOperator`) — open
  since pre-R-1; G3f retirement is the terminal close-out.

### 2.2 Convention crosswalk for the legacy fixed-source carve

Same template as Phase 0.2, applied to `_solve_fixed_source_krylov` and `_solve_fixed_source_si`. Specifically enumerate:

- `external_source` shape and convention (per-ordinate? iso? `/W`-applied?)
- `_build_rhs_*` output shape and convention (packed-1D, `/W`-applied iso + per-ordinate aniso? what?)
- `solution` shape (packed-1D with what eq_map?)
- `phi` shape and convention (scalar flux `(ng, nx, ny)`)

This is the same artifact as Phase 0.2 but for the fixed-source path. **It must exist before any Step G1/G2 code lands.**

**Acceptance**: crosswalk written to `.claude/plans/r1_step4_g_convention_crosswalk.md`.

**Commit**: `docs(plan): R-1 Step G convention crosswalk for solve_sn_fixed_source`

### 2.3 Dispatch `test-architect` to design Step G1/G2 verification plan

**Brief**: design the L0 + L1 verification plan for the deferred carves. Specifically: which existing tests pin the legacy `solve_sn_fixed_source` behaviour, which new tests are needed for the typed-flux carve, which L1 MMS gates must continue passing.

**Acceptance**: `test-architect` returns a verification plan that lands as `.claude/plans/r1_step4_g_verification_plan.md`. Session 2 implementation follows this plan.

---

## Phase 3 — G1 / G2: carve `solve_sn_fixed_source` onto typed contract

**Goal**: replicate Steps D / E for the fixed-source path. After this phase, `solve_sn_fixed_source` consumes the typed operator algebra. The legacy symbols are unused but not yet deleted.

### 3.1 G1 — carve `_solve_fixed_source_krylov` onto `KrylovAcceleration + AngularFlux`

Structurally identical to Step D's carve of `SNSolver._solve_krylov`:

- Build typed `q_ext` from `external_source` (per-ordinate, post-A1 no `/sum_w` rescaling needed for S).
- Build operator triple `(LC, S, ZeroOperator)`.
- `preconditioner=lambda q: q` (explicit identity — kept for L1 anchor bit-identity until #200 ships the production preconditioner choice; the silent-fallback bug class that originally motivated the explicit identity is GONE post-Phase-1.2, see §1.2).
- One `KrylovAcceleration.solve(q_ext)` call replacing the outer source-iteration + GMRES inner combo.

The legacy outer SI wrap (line 1464 `for n_outer in range(max_inner)`) collapses — KrylovAcceleration solves `(L+C-S)·ψ = q_ext` directly, no outer needed.

**Verify**: MMS curvilinear convergence rate (`tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`) and the slab MMS path.

**Commit**: `refactor(sn): _solve_fixed_source_krylov onto KrylovAcceleration + typed AngularFlux (R-1 Step 4 G1)`

### 3.2 G2 — carve `_solve_fixed_source_si` onto `SourceIteration + AngularFlux`

Structurally identical to Step E's carve. Same operator triple. SI loop's `_solve_with_seed` forwards the previous iterate via `initial_guess` explicitly (post-Phase-1.2 — `_attach_previous` retired in `c93355c`; AngularFlux history machinery preserved for time derivatives, no longer load-bearing for the Carlson seed path).

**Verify**: same MMS gates as G1.

**Commit**: `refactor(sn): _solve_fixed_source_si onto SourceIteration + typed AngularFlux (R-1 Step 4 G2)`

---

## Phase 4 — G3: retire the legacy symbols

**Goal**: ALL legacy operators / decoders / RHS builders / preconditioners are gone. The retirement is principled because Phase 3 leaves no production callers.

### 4.1 G3a — retire `_make_sweep_preconditioner` + `_build_rhs_*`

Free functions and SNSolver method. Zero production callers after Phase 3. Tests that reference these are deleted (none currently — these are internal helpers).

**Commit**: `refactor(sn): retire _make_sweep_preconditioner + _build_rhs_* helpers (R-1 Step 4 G3a)`

### 4.2 G3b — retire `solution_to_angular_flux*` and `build_equation_map*` legacy decoders

The `_with_traces` versions stay (used by `AngularFlux.from_flat_with_traces`). The non-`_with_traces` versions retire.

Delete corresponding test files / sites that exercised the legacy decoder directly.

**Commit**: `refactor(sn): retire legacy decoders + equation_map builders (R-1 Step 4 G3b)`

### 4.3 G3c — retire `SNStreamingOperator`

Now nothing in production reads `solver.L`. The class deletes.

`tests/sn/test_snstreamingoperator.py` deletes entirely (superseded by typed-operator tests).
`tests/sn/test_streaming_operator.py` `_packed_psi` migrates to typed AngularFlux (same pattern as Step F for collision).

**Commit**: `refactor(sn): retire SNStreamingOperator (R-1 Step 4 G3c)`

### 4.4 G3d — retire `EquationMap` class

After G3a/b/c, only `EquationMap` is left from the legacy module. Delete the class and its helpers.

**Commit**: `refactor(sn): retire EquationMap class — legacy packed-1D eq_map gone (R-1 Step 4 G3d)`

---

## Phase 5 — G4: relocate the keepers

`pack_with_traces`, `solution_to_angular_flux_with_traces`, `build_equation_map_with_traces` are STILL used (by `AngularFlux.from_flat_with_traces` / `to_flat_with_traces`). Move them to be private helpers of `AngularFlux`.

**File**: `orpheus/sn/angular_flux.py` (gain three private static methods or module-private functions). `orpheus/sn/operator.py` loses them.

**Verify**: `AngularFlux.from_flat_with_traces` round-trip identity test continues to pass.

**Commit**: `refactor(sn): relocate _with_traces helpers into angular_flux module (R-1 Step 4 G4)`

---

## Phase 6 — H: Sphinx narrative + closeout

**Goal**: durable documentation of the operator algebra, the typed contract, the convention crosswalk, the four AI-failure-mode catches.

### 6.1 Sphinx page for the four-operator algebra

`docs/theory/sn_operator_algebra.rst` (new or extended) with:

- The full algebra `(L + C - S - F) ψ = q_ext` derivation.
- The typed contract: `AngularFlux` as the only shape; `InvertibleOperator` as the algebraic identity; the convention crosswalk.
- The Carlson seed / partner-flux state plumbing via `rhs(1)`.
- The /W convention (now at the producer).
- The known limitations (#200 sphere-4g-krylov, #201 typed sources).
- References: Lewis & Miller §3.2, Adams & Larsen 2002 §III, Hébert §3.9.4.

**Acceptance**: `sphinx-build docs docs/_build/html -W` clean.

**Commit**: `docs(sn): operator-algebra narrative for the typed AngularFlux carve (R-1 Step 4 H)`

### 6.2 Memory entries

Update sub-agent memories for `method-implementer`, `numerics-investigator`, `archivist` with the lessons (most live in `.claude/lessons.md` already from Phase 0; this step is just the cross-references).

**Acceptance**: agents have pointers to the lessons + plan + crosswalk artefacts.

**Commit**: `docs(agent-memory): R-1 Step 4 closeout + lesson cross-references`

### 6.3 Final L1 sweep + retirement audit

Run full `tests/sn/l1_analytical/` and `tests/sn/regression/`. Grep for any remaining references to retired symbols (only OK in agent memory and historical commits).

**Acceptance**: full SN suite green; `grep -rn "SNStreamingOperator\|EquationMap\|build_equation_map[^_]\|solution_to_angular_flux[^_]"` returns only doc/memory/history matches.

---

## Issue tracking — file these before starting Phase 1

Cross-session containment requires GitHub issues per Cardinal Rule 4:

- **A1 issue**: "ScatteringOperator typed apply should normalise iso at producer (Pattern 7)" — body cites session-1 diagnosis; closes when Phase 1.1 lands.
- **A2 issue**: "KrylovAcceleration default precond should require CAP_STATELESS_INVERSE" — body cites Step D bug; closes when Phase 1.2 lands.
- **A3 issue (already exists as #201)** — leave open; this plan does NOT close it. The typed Source distinction is bigger scope (Phase A).
- **A5 issue**: "Promote R-1 Step 4 diagnostics to permanent L1 tests" — closes when Phase 1.3 lands.
- **Step G omnibus issue**: "R-1 Step 4 G + H — retire legacy SN operators + Sphinx closeout" — body links this plan; closes when Phase 6 lands.

---

## Scope NOT in this plan

- **A3 (typed `Source` distinction)**: stays at issue #201; Phase A absorbs it.
- **MoC / CP analogous carves**: out of scope; the lessons in `.claude/lessons.md` shape those when they happen.
- **The pre-existing 39 packed-vector failures in non-`test_collision_operator.py` files**: those tests exercise `SNStreamingOperator` directly; they retire with Step G3c (Phase 4.3), no separate migration.

## Reversibility / safety

Every phase is independently revertible. Phases 0–2 don't touch production. Phases 3–5 each ship a self-contained carve / retirement; reverts cap at one phase boundary.

Phase 3 (G1 / G2) is the highest-risk phase: it touches `solve_sn_fixed_source` (public MMS API). The MMS L1 convergence-rate tests are the primary gate. If MMS convergence rate degrades on curvilinear, the convention crosswalk from Phase 2.2 is the diagnostic. If it doesn't degrade but a specific MMS case fails, follow the Step D / Step E diagnostic pattern (manual outer loop, print boundary state, check sum_w bridge).

## Estimated effort (full session 2)

- Phase 0: 1 hour (mostly text)
- Phase 1: 3-4 hours (three small refactors + test migrations)
- Phase 2: 2 hours (audit + crosswalk + test-architect dispatch)
- Phase 3 (G1/G2): 4-6 hours (carves + L1 verification)
- Phase 4 (G3): 2 hours (mechanical retirement; tests already migrated)
- Phase 5 (G4): 1 hour
- Phase 6 (H): 3-4 hours (Sphinx narrative + final checks)

---

## 2026-05-21 reality check

Phase 1.1's "3-4 hours" estimate was for the original (conservative) A1
plan.  Mid-session the scope expanded to "AGGRESSIVE A1" per user direction:

- ``@singledispatchmethod`` on ``ScatteringOperator.apply`` AND
  ``FissionOperator.apply`` (vs the original "just rescale at the producer").
- ``PerOrdinateSource.from_isotropic`` factory.
- ``transport_sweep`` signature change (drop ``iso_source``).
- Both 1D unified AND 2D wavefront sweep paths migrate.
- ``InvertibleOperator.solve`` bridge drop + docstring rewrite.
- ``SNStreamingOperator.solve`` signature change.
- Legacy ``_solve_fixed_source_si`` + ``_solve_fixed_source_krylov`` +
  ``_make_sweep_preconditioner`` migrate.
- 8 test files migrate.
- All decided in surgical mode (no method-implementer; turn-by-turn user
  steering).

The actual Phase 1.1 scope expanded to ~30 sites / ~6 production files /
~8 test files.  Mid-flight at context-handoff (~30 commits-worth of
changes in 1 atomic commit per Path B).

The IN-FLIGHT state document
(``.claude/plans/r1_step4_session2_state.md``) is the resuming-fresh
anchor.  Re-estimating remaining Phase 1.1 effort: ~1-2 hours
(wide-gate triage + Group E new L0 tests + commit).  Phase 1.2 / 1.3 /
1.4 / 2 estimates above remain valid.

Total: ~16-20 hours = 2-3 working sessions. Recommend splitting:
- **Session 2a**: Phase 0 + Phase 1 + Phase 2 (cleanup + audit).  Ships a cleaner foundation.
- **Session 2b**: Phase 3 (G1/G2 carves).
- **Session 2c**: Phase 4 + Phase 5 + Phase 6 (retire + relocate + narrative).
