# B.5.2 — Holistic operator-output `BoundaryResidual` retype (Wave O / #208)

**Branch:** `refactor/field-role-typing` (worktree
`.claude/worktrees/field-role-typing/`). **Base HEAD when this plan was
written:** `cd29d69` (O.4b E3 verification gates committed).
**Mode:** main-agent DIRECT authorship, turn-by-turn user steering
(`feedback_no_method_implementer_for_surgical_carves`) — this carve touches the
GREEN 1-D SI/Krylov drivers, so it is high-stakes surgical work. Commit per
sub-step.
**Standing exclusions from commits:** `.claude/agent-memory/`, `.mcp.json`,
hooks, `derivations/diagnostics/`.

**Pre-flight reads (do these FIRST):**
- `.claude/agent-memory/explorer/issue_208_b52_boundary_residual_retype_surface.md`
  — THE ground-truth surface map (HEAD `02f8f0a`). Q1 (the 12 emission sites
  table), Q2 (the two-hat tension), Q3 (the type-agnostic flat round-trip
  escape hatch), Q5 (the test-migration list), Q6 (the emission pattern), Q7 +
  Synthesis (the O.2 boundary + the resolution shapes). **This plan is the
  execution synthesis of that map — read the map for every `file:line`.**
- `.claude/agent-memory/test-architect/issue_208_o4b_2d_verification_plan.md`
  — the §bit-identity-vs-principled-equivalence framing (the retype is a TYPE
  change with IDENTICAL `.values`, so value-based tests pass; only output-TYPE
  asserts migrate).
- `orpheus/transport/residuals/boundary_residual.py` (the target leaf, its
  `from_balance`), `orpheus/transport/source_sinks/boundary_source_sink.py`
  (the source-role leaf for the SI rhs), `orpheus/transport/fields/_bases.py`
  (`BoundaryField` arithmetic / class gate).

---

## 0. What this carve is (and is NOT)

**The "dimensional sin", boundary half.** Every operator's `.apply` output
`.boundary` is currently mistyped `BoundaryFlux`. The bulk half of the fix
(`.bulk` → `AngularSourceSink`) ALREADY LANDED (commit `f400743`, follow-up
`81289b6` — every `# B.5.2:` comment in the code is that bulk work). This carve
is the **boundary half**: retype the operator-output `.boundary` to its honest
role leaf.

The honest types (the role grid, plan `wave_o_operator_typing.md:652-685`):
```
matvec/composed OUTPUT boundary  →  BoundaryResidual   (a defect r_Γ = γ₋ψ − (R·G·γ₊ψ + q))
SI/Krylov source rhs boundary    →  BoundarySourceSink (the inflow source q.boundary + B·ψ.outflow)
SOLVE output boundary (sweep)    →  BoundaryFlux       (the solution trace — UNCHANGED, leave alone)
```

**Why it can land NOW (the blocker is gone).** A prior naive attempt
(`wave_o_operator_typing.md:176-197`) failed with `TimedFullField boundary type
mismatch: BoundaryResidual vs BoundaryFlux` because the PRE-extraction boundary
slot was OVERLOADED (matvec defect AND SI inflow-flux seed at once). **BC
extraction has SINCE landed** (1-D O.4a.2 `4c0ff96`/`2bdc66d`; 2-D O.4b
`7638d1f`/`dfeb604`/E3 `cd29d69`). Extraction de-overloads the slot by
construction → the honest typing is now a consequence that CAN land.

**NOT in this carve (defers to O.2 proper):** the `|Ω·n|·w` G-metric adjoint
(`.H = G⁻¹AᵀG`), Gate-1.3 green, the full `L+C−S−F−B` driver consuming the loss
operator (retiring the `_scattering_with_boundary_op` `S+B` fold + the
`_reflect_outflow_into_inflow` helper twin), the white-BC adjoint. The retype
needs NONE of these.

---

## 1. THE CRUX — the two-hat tension and its resolution

`SNBoundaryOperator.apply` (`B`) emits a **non-zero** boundary `B·ψ.outflow`
(the reflective inflow). It is consumed TWO ways:

- **Krylov matvec** (`numerics/iteration.py:676-683`):
  `out = L.apply(ψ) − (S+B).apply(ψ) − F.apply(ψ)` — typed TimedFullField
  arithmetic. The result boundary is the **residual** → wants `BoundaryResidual`.
  The `BoundaryField` class gate (`timed_full_field.py:280`, strict
  `type(a) is type(b)`) forces L, (S+B), F output boundaries to SHARE a class →
  all `BoundaryResidual`.
- **SI rhs** (`numerics/iteration.py:455`):
  `rhs = F.apply(ψ) + (S+B).apply(ψ) + q_ext` — the SOURCE of `L·ψ = rhs`. The
  sweep CONSUMES `rhs.boundary` as the inflow SEED
  (`operator.py:1968-1986`, verified). So `rhs.boundary` = `B·ψ.outflow +
  q_ext.boundary` is a **source** → wants `BoundarySourceSink`.

**One operator cannot emit `BoundaryResidual` (Krylov) AND `BoundarySourceSink`
(SI) at once.** That is the two-hat problem — the same source/residual overload
that #208 exists to fix, moved up one level.

### The resolution (the chosen shape — Shape B-minimal)

**`B.apply` → `BoundaryResidual` (one hat, the Krylov/residual role). The SI
driver STOPS sourcing its boundary seed from `(S+B).apply.boundary`; it composes
the inflow source separately as a `BoundarySourceSink`.**

Concretely:
- Krylov: `L.apply − (S+B).apply − F.apply` composes
  `BoundaryResidual − BoundaryResidual − BoundaryResidual` (closed). Ravels via
  `to_flat` (type-agnostic, Q3); the GMRES iterate reconstructs as `BoundaryFlux`
  off the flux `solution_template` (Q3 escape hatch) — the solution side is
  UNTOUCHED.
- SI: the driver builds `rhs.boundary` as a `BoundarySourceSink` =
  `q_inflow_source ⊕ (B·ψ.outflow as a source)`. It does NOT add
  `(S+B).apply.boundary` (now a residual) to `q_ext`. The `B·ψ.outflow` inflow
  source comes from a **source-role reflection** — either a new
  `BoundarySourceSink`-typed entry point on `SNBoundaryOperator` (e.g.
  `reflect_into_inflow_source(ψ) → BoundarySourceSink`, the source sibling of
  `apply`), or by re-typing the values produced by the existing
  `_reflect_outflow_into_inflow` helper (`solver.py:959`) as `BoundarySourceSink`.

This is the START of O.2's driver split (the two `−B` routes begin to separate
by role: the matvec route stays residual; the SI source route becomes a typed
source). O.2 finishes it (the `L+C−S−F−B` driver + fold retirement).

**Why this is the minimal honest landing (per the map's verdict §e):** the
matvec-residual typing and the SI-source typing are COUPLED through B and cannot
be honestly separated. The smallest green-keeping cut is exactly the three
pieces above. There is NO "matvec-only retype" that keeps green — the moment
`L.apply` → `BoundaryResidual`, the Krylov sum drags B in (class gate), and B in
the SI sum drags the SI-source typing in.

---

## B0 — DE-RISK FIRST (the prophylactic; mirror E0/D0)

Before ANY production edit, a throwaway diagnostic in
`derivations/diagnostics/` (excluded). Prove on small 1-D slab + 2-D Cartesian
meshes, reflective + vacuum:
1. **No cross-class throw.** Build the SI rhs the NEW way (bulk source ⊕
   `BoundarySourceSink` inflow source) and the Krylov matvec the NEW way (all
   `BoundaryResidual` output) — confirm neither hits the
   `TimedFullField._check_partner` gate (`timed_full_field.py:280`). This is the
   HIGH risk (§Risk 1) — the exact prior failure.
2. **Values unchanged.** The retype is a TYPE change only: assert the SI rhs
   `.boundary.values` and the Krylov matvec `.boundary.values` are BIT-IDENTICAL
   to the current (BoundaryFlux) path. The reflection `R·G·ψ.outflow` is the same
   arithmetic; only the wrapping class differs.
3. **Round-trip.** Confirm the Krylov `to_flat`/`from_flat` round-trips a
   `BoundaryResidual` output back to a `BoundaryFlux` iterate (Q3 escape hatch),
   bit-identical to today.
If any fails, STOP and re-examine the resolution before touching production.
(The `except TypeError` un-mask at `iteration.py:737` means a cross-class throw
surfaces LOUDLY now — good, it is the tripwire.)

---

## 2. Execution steps (ordered; all `file:line` at the map's HEAD `02f8f0a`)

The matvec-side group (steps 1–4) is **all-or-nothing** — every leaf's output
boundary in the Krylov sum must flip together. Land steps 1–5 + the SI-source
restructure (step 6) in ONE commit (the green gate is the full composed path);
the test migrations (step 7) co-land.

**Step 1 — L's 1-D matvec residual → `BoundaryResidual`.**
- `operator.py:543` build / `:570` return — `_MSpatialOperatorSum._compute_LpC`
  `m_boundary` (THE 1-D `(L+C)` boundary keystone; outflow=defect
  `streamed−given`, inflow=identity `given`).
- `operator.py:847` build / `:874` return — `_compute_decomposition`
  `m_spat_boundary` (the dual-emission TWIN of `_compute_LpC` — **MUST flip
  IDENTICALLY** or the Resolution-A bit-exact decomposition test breaks; §Risk 2).
- `operator.py:278` build / `:287` return — `_SpatialSweepDirection.apply`
  `masked_boundary` (per-direction summand, the slow standalone path).
- Mechanical: swap `BoundaryFlux.zeros_on(...)` → `BoundaryResidual.zeros_on(...)`
  + the `from orpheus.transport.fields.boundary_flux import BoundaryFlux` →
  add `from orpheus.transport.residuals import BoundaryResidual` in each method's
  import block. The face_view writes (`:546-565`) are UNCHANGED (same values).
- `operator.py:1274` (`StreamingOperator.apply` 1-D) inherits via
  `boundary=LpC_result.boundary` — no own alloc, no edit.

**Step 2 — L's 2-D matvec residual → `BoundaryResidual`.**
- `operator.py:1460` build / `:1474` return — `_apply_2d_cartesian` `out_boundary`
  (the O.4b Phase E bare emission; 4 faces, outflow=defect, inflow=identity).
  Swap the class. **WILL trip the A2D-1 source-hash pin** (§Risk 3) — refresh it.

**Step 3 — M_angular zero → `BoundaryResidual.zeros_on`.**
- `operator.py:880` — `_compute_decomposition` `m_ang_tff` (M_angular_redist,
  zero boundary, BulkOperator per MA-Q4).

**Step 4 — C/S/F zero boundaries → `BoundaryResidual.zeros_on`.**
- `operator.py:1592` — `CollisionOperator.apply`.
- `scattering.py:1150` — `ScatteringOperator.apply` (TimedFullField branch).
- `fission.py:321` — `FissionOperator.apply` (TimedFullField branch).
- **DO NOT touch** `operator.py:1623` (`CollisionOperator.solve` — flux codomain)
  or `operator.py:2010` (`InvertibleOperator._solve_timed_full_field`
  `boundary_buf` — the SWEEP solution trace, a `BoundaryFlux`). These are SOLVE
  outputs (fluxes), not residuals.

**Step 5 — B's output → `BoundaryResidual`.**
- `boundary_operator.py:183` build / `:197` return — `SNBoundaryOperator._apply_faces`
  (both `apply` and `apply_transpose`). Swap to `BoundaryResidual`. Now `B.apply`
  wears ONE hat (the Krylov/residual role).

**Step 6 — the SI-source restructure (THE crux; do with B0 in hand).**
- `q_ext.boundary` → `BoundarySourceSink.zeros_on` at the source-build sites:
  `solver.py:616` (SI `q_ext_composite`), and the fixed-source builds
  `solver.py:748` / `:1576` (verify exact lines at execution HEAD).
- The SI rhs boundary composition: the driver must compose `rhs.boundary` as a
  `BoundarySourceSink` (`q_inflow_source ⊕ B·ψ.outflow-as-source`), NOT via
  `(S+B).apply.boundary`. Add the source-role reflection entry point
  (`SNBoundaryOperator.reflect_into_inflow_source(ψ) → BoundarySourceSink`, or
  re-type `_reflect_outflow_into_inflow`'s written values as `BoundarySourceSink`).
  Touch points: `_solve_source_iteration` (eigenvalue SI, the `S+B` fold route,
  `solver.py:643`), `_solve_fixed_source_si` (`solver.py:~1409` direct-loop
  `_reflect_outflow_into_inflow` route), the final reconstruction sweep
  (`solver.py:~1131`). The SWEEP consumes `rhs.boundary` type-agnostically
  (`operator.py:1968-1986` copies `.face_view` into a fresh buffer) — so the seed
  still reaches the sweep once the class is consistent.
- **Both `−B` routes** (the `S+B` fold for eigenvalue SI; the
  `_reflect_outflow_into_inflow` helper for the direct loops) must be
  source-typed consistently. (They collapse fully at O.2.)

**Step 7 — test migrations (co-land with the production commit).**
- `tests/sn/operators/test_native_matvec.py:320` — MIGRATE
  `assert isinstance(result.boundary, BoundaryFlux)` → `BoundaryResidual`
  (`result = (L+C).apply`). The ONE hard matvec-output-TYPE assert.
- `tests/sn/operators/test_bc_extraction_2d.py` (`TestBoundaryResidual2D*`,
  `:362/:442`) — ADD a positive type assert that the 2-D matvec output boundary
  IS `BoundaryResidual` (the natural home; the E3 gates currently read values
  type-agnostically — now pin the type too).
- `tests/sn/operators/test_sn_boundary_operator.py` (`:140/:167`) — ADD a type
  assert that `B.apply` output is `BoundaryResidual`.
- Refresh the **A2D-1 source-hash pin** (`test_streaming_operator.py`
  `TestT4dApply2DCartesianSourceHashPin.EXPECTED_SHA256` + a History entry) —
  editing `_apply_2d_cartesian` trips it (intended tripwire).
- `tests/transport/test_timed_full_field.py:87` — SURVIVES (it asserts the
  zeros-FACTORY arg, not an operator output). Do not touch.

---

## 3. Verification — must-stay-green (1-D AND 2-D, SI AND Krylov)

Run at EVERY sub-commit (the carve changes the matvec TYPE flow on both paths):
- **1-D SI + Krylov:** `tests/sn/solve/test_fixed_source_g1.py`,
  `test_b1pp_verification.py` (L+C composition `:137/:224/:288`; null-space
  `.boundary.values==0` `:237` — value-based, survives).
- **1-D matvec + the TWIN:** `test_native_matvec.py` (`:320` migrated),
  `test_streaming_operator.py`, `test_streaming_operator_decomposition.py`
  (Resolution-A bit-exact — the `_compute_LpC`/`_compute_decomposition` twin MUST
  flip identically, §Risk 2).
- **2-D:** `test_bc_extraction_2d.py` (E3 gates + the new type assert),
  `test_2d_l2_matvec_correctness.py` (Gate K k_inf=1.875),
  `test_streaming_operator.py` A2D-1 (refresh), the octant snapshots
  (value-based, survive).
- **Operator leaf zeros:** `test_operators_apply_typed.py` (F/S/C zero-boundary
  VALUE asserts — survive), `test_collision_operator.py`,
  `test_scattering_operator.py`, `test_fission_operator.py`.
- **B:** `test_sn_boundary_operator.py` (face_view value asserts survive + the
  new type assert).
- **Type system:** `test_field_units.py` (BoundaryResidual units),
  `test_typed_residuals.py` (`from_balance`),
  `test_boundary_source_sink_residual.py` (cross-class raise pins).
- **Sentinel** 36/36 (NO `-O`); homogeneous-2G keff; sweep L1 MMS (slow —
  the matvec TYPE flow changed; run them).
- Default mode `python -O`; env `PYTHONPATH=$PWD
  /Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python`. **NO `continuous_get`**
  (#212 hang). Bash gotcha: `> file 2>&1; echo exit=$?` (never `| tail; echo $?`).

**KNOWN PRE-EXISTING REDS — NOT caused by this carve (do not chase):**
1. Cylinder matvec `tests/sn/sweep/curvilinear/test_unified_matvec_cylinder.py`
   (27 cases, 1-D, #206/#196, bisected to `7072f9b`).
2. **#212** `continuous_get` hang (`test_heterogeneous_transport`,
   `test_keff_slab`, `test_mms_heterogeneous`). Fix patch:
   `derivations/diagnostics/sn_keff_hang_PROPOSED_fix.patch`.

---

## 4. Risk ranking (most → least likely to break green 1-D)

1. **The two-hat / SI-source cross-class throw (HIGH).** The moment B →
   `BoundaryResidual`, the OLD SI rhs `F + (S+B) + q_ext` throws
   `TimedFullField boundary type mismatch` (the exact prior failure,
   `wave_o_operator_typing.md:178-179`). The `except TypeError` un-mask
   (`iteration.py:737`) surfaces it LOUDLY. **Mitigation = the §1 resolution
   (step 6: source-type the SI rhs, don't add operator-output boundaries to the
   source). De-risk with B0 first.**
2. **The 1-D twin (MEDIUM).** `_compute_LpC` (`:543`) and `_compute_decomposition`
   (`:847`) are Cardinal-Rule-2 twins — BOTH must flip identically or the
   Resolution-A bit-exact test (`test_streaming_operator_decomposition.py`) breaks.
3. **A2D-1 source-hash pin (MEDIUM, deterministic).** Editing `:1460` trips it —
   refresh the hash + History entry (expected, not a bug).
4. **`test_native_matvec.py:320` (LOW, mechanical)** — one type assert to migrate.
5. **Solution build (NONE).** Verified type-agnostic (map Q5); `final_boundary`
   stays `BoundaryFlux` (it is the SOLVE trace, not a residual).

---

## 5. Pickup checklist (fresh post-compaction session)

1. Confirm branch `refactor/field-role-typing`, `git log --oneline -3` shows
   `cd29d69` (E3 gates) at/near HEAD, tree clean of tracked non-excluded paths.
2. Read the explorer map (Q1 table = the exact sites; Synthesis (b) = the SI-rhs
   resolution shapes; (d) = risk ranking) and this plan §1 (the crux).
3. **B0 de-risk** (no cross-class throw + values bit-identical + round-trip) —
   STOP if it fails.
4. Steps 1–6 in ONE production commit (the green gate is the composed path) +
   step 7 test migrations co-landed. Turn-by-turn with the user (high-stakes,
   touches green 1-D drivers).
5. Verify §3 must-stay-green (1-D SI+Krylov, 2-D, twin, A2D-1, sentinel).
6. Update `wave_o_operator_typing.md` + the auto-memory: B.5.2-boundary DONE;
   what O.2 still owns (G-metric adjoint, fold retirement, the two-route
   collapse, white-BC adjoint, Gate-1.3).
7. Sphinx: the `operator_algebra.rst` boundary-typing section (archivist) — the
   honest role grid (BoundaryFlux=solution / BoundaryResidual=defect /
   BoundarySourceSink=source) now realized on the boundary.

---

## 6. Why this is the honest minimum (verdict, from the map §e)

A "matvec-only retype that keeps green WITHOUT touching the SI driver" is **NOT
possible** — the two-hat coupling through B forces at least the SI-source
source-typing. The coherent standalone landing is exactly:
`matvec output → BoundaryResidual` + `SI/Krylov source → BoundarySourceSink` +
`B's SI contribution source-typed`. It completes the `f400743` bulk retype on
the boundary and STARTS (does not finish) O.2's driver split. Everything else in
O.2 (G-metric adjoint, `.H`, Gate 1.3, the `L+C−S−F−B` driver + fold/helper
retirement, white-BC adjoint) is genuinely separable and defers.
