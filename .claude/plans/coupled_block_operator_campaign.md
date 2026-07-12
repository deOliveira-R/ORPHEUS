# Coupled block operator — the ψ½ Ray-Characteristic system realized in operator algebra

**Campaign plan — 2026-07-07.** Posing user-confirmed. **Supersedes** the FaceField carve's
remaining phases (the guard-retirement C3, the `StartingDirection→RadialCharacteristic`
rename C4, and the docs C5 all fold in here) and **subsumes** the coupled-system machinery.
Substrate: the **C2 FaceField carve LANDED** (`c637407` C2a / `4081c0d` C2b / `be5e7f8` C2c
on branch `refactor/sn-walk-unification`; the ψ½ field is `StartingDirectionField(FaceField
[tuple[int,int,str]])`, sibling of the spatial `BoundaryField`, `G_sd = V_cell` metric fix).

**Binding lens:** [[feedback-build-the-machinery-operator-algebra]] — the goal is the
MACHINERY (LEAN but LARGE), realize the operator-algebra formalism ALWAYS (the four
questions: *what is it in operator algebra / how is it posed / what invariant tests it / how
is it solved*). Risk/effort set the SEQUENCE (small steps so a bug stays findable), never the
WHETHER. The welded un-named `A_BA` is a defect in the machinery (a coupling block never
posed in the algebra), not an optional cleanup.

---

## The object (posing — confirmed 2026-07-07)

The within-group augmented SN system is a **2×2 coupled block operator** over two systems
(the composite biproduct re-partitioned — `Mat₂(Mat₂(𝒞)) ≅ Mat₄(𝒞)`; NOT a new type):

```
[ A_AA   A_AB ] [ transport ]      A_AA = L + C − S − B     (System A operator)
[ A_BA   A_BB ] [ ray       ]      A_BB = RayOp             (System B operator)
```

- **System A (transport):** state = `FullField` (`BulkField` bulk angular flux ⊕ spatial
  `BoundaryField` as its BC). `A_AA = L+C−S−B`.
- **System B (ray):** state = `RadialCharacteristicFlux` (a codim-1 `FaceField` — the ψ½
  cells at each radial cell) carrying **TWO boundary conditions**: r=R Dirichlet (the seed /
  given inflow data) **and** r=0 pole-reflection (a two-point BVP). `A_BB = RayOp` (the banded
  radial straight-characteristic march μ∂_r+σ_t, Hébert 3.434–3.435, + the two BCs).
- **`A_AB`** (ray→transport): the ray seeds the bulk angular recurrence — **posed as an operator**.
- **`A_BA`** (transport→ray): the bulk emission folded onto the ray source — the **un-named
  Schur elimination** currently hand-rolled into S/F — **posed as an operator** (the fold =
  named machinery, single-sourced via `fold_moments_to_starting_direction`).

**"Deciding where the block goes IS operator algebra":** once all four are posed as
operators, whether `A_BA` sits in the resolvent (folded), the explicit block operator, or a
DSA preconditioner is a *composition choice* the machinery supports. That flexibility is the
realization. The block-SOLVE stays spectral-character-keyed (the existing `InvertibleOperator`
/ `SourceIteration` / `MatrixInverseOperator` family — resolvent block-triangular-direct +
scattering iterative), which is already correctly realized.

## Specialist-verified facts (the measured algebra — 3 dispatches, 2026-07-07)

- **Within `(L+C)` the ray↔bulk coupling is block-triangular, solved DIRECTLY** (ρ=0, one
  pass): measured `A_sb = 0` EXACT (seed←bulk within L+C), `A_bs = 7.5` (`A_AB` M-M feed),
  `A_ss = 5.0` (`A_BB` march); `solve(apply(ψ))−ψ = 3.5e-16` (#284 forward-substitution).
- **Only the scattering coupling iterates:** `A_BA` lives ENTIRELY in the lagged gain S/F —
  measured `S_sb = 0.183` (bulk→ray source), `S_bs = 0` EXACT (ψ½ zero moment weight),
  `ρ(M⁻¹N) = 0.371 < c` (outer SI, Adams–Larsen). So the API reflects BOTH a block-triangular
  DIRECT solve (ray/bulk, in the resolvent) AND a block-ITERATIVE solve (scattering).
- **Extract is principled-equivalent, not bit-identical** (5.5e-16, different reduction tree);
  a NAIVE extraction returns O(1) garbage (the sweep treats inflow/seed rows as given data —
  the row-contract MUST be preserved). Principled-equivalence + the invariant test is the bar
  ([[feedback-principled-over-bit-identical]]).
- **The corrections that stand:** System B is a **two-point BVP** (r=R + r=0-pole → 2 BC data);
  the block-solve is **spectral-keyed** (reuse existing solvers); the ρ-honest **block residual**
  `r = Aψ − q` (via `evaluate_residual`, `solver.py:231`) is the convergence metric, NOT
  `‖Δψ‖`; **DSA** is a future CONSUMER (its R/P operators are a separate feature — not a
  risk-deferral); **RQI** (KKT/indefinite) and **multiphysics** (nonlinear/Newton) live in
  DIFFERENT machinery — out of scope.
- **Pre-carve baselines banked:** `derivations/diagnostics/diag_coupled_0{1,2}_*.py` (the
  block-triangular certificate + welded=exact-inverse) — promote to
  `tests/sn/operators/test_psi_half_coupling.py`.
- **Stale doc to fix:** the `"dormant until d3"` comments (`scattering.py:1562`,
  `fission.py:508`) — d3 LANDED @ `a29ab2d`; the A_BA arms are LIVE.

## The sequence — small de-risking steps, each a COMPLETE operator-algebra realization + its invariant

Each step answers the four questions and lands its invariant test (test-architect designs the
suite). Principled-equivalent, verified — not necessarily bit-identical.

1. **Pose System B.** `RadialCharacteristicFlux` (the `StartingDirection→RadialCharacteristic`
   rename lands here) + its two BCs (r=R Dirichlet, r=0 pole-reflection) + `RayOp = A_BB`.
   *Invariant:* the radial two-point BVP is well-posed; the straight-char march IS its solve
   (the r=R + r=0 data close it). *Corner→BC:* the corner is already separately-keyed;
   `_reflect_starting_direction` (`boundary.py:300`) is already its BC arm.
2. **Name `A_BA`.** Un-weld the Schur fold out of S (`scattering.py:1571`) and F
   (`fission.py:513`) into a named coupling operator (reusing `fold_moments_to_starting_
   direction`); fix the stale `"dormant"` comments. *Invariant:* it reproduces the S/F seed
   arms (single-source fold; `A_BA` = Fold ∘ bulk-emission). Low-risk (lagged gain OUTSIDE the
   resolvent — cannot touch #284).
3. **Name `A_AB`.** Pose the ray→bulk seed injection (`precompute_psi_state` reading
   `cells(p,−1)`) as an operator. *Invariant:* matches the in-sweep injection; the resolvent
   realizes it as forward-substitution; the adjoint `seed_cells_bar` is `A_AB.H`.
4. **Assemble `CoupledOperator [[A_AA,A_AB],[A_BA,A_BB]]`** + widen the role lattice
   (`SystemRole {A,B,COUPLED}` cloning the set-union `_join_block_roles`; System B first-class
   — it is a LIVE 3rd block absent from today's `BlockRole`). *Invariant:* block matvec ≡ the
   current augmented apply; block `.H` ≡ the current adjoint (the `(A.H).inverse()=(A.inverse
   ()).H` swap law extends). **Highest-risk step** if it touches the interleaved (L+C) walk —
   gate hard; may WRAP the sweep as the resolvent's block-triangular solve of the posed
   operators (naming, not re-implementing) as the first landing, extraction second.
   **INCLUDES THE LIFT** (user ruling 2026-07-08): S/F → pure bulk, `A_BA = Fold ∘ K_iso` composed
   by the driver, the fold migrates to sn beside A_BB — see "Step-4 LIFT deliverable" below.
5. **The solve.** Pose the resolvent (block-triangular DIRECT) + scattering (ITERATIVE) via
   the spectral-keyed solvers; the ρ-honest block residual as the stopping test. *Invariant:*
   coupled solve ≡ current solve (principled-equivalent, 5.5e-16 vs the dense-LU oracle).
6. **Retire the guards + `Optional` + the two-channel collapse.** Presence = the mesh-
   enumerated block-list / the operator's arity (structural, mismatch unconstructable); retire
   the 7 `_require/_refuse_starting_direction` guards + `FullField.Optional[starting_direction]`;
   retire the `sweep()` `starting_direction_source/_flux` kwargs (the seed rides only the
   composite — user ruling 2026-07-06). *Invariant:* a carrying-mesh composite missing the ray
   block is unconstructable; a non-carrying one carrying it is unconstructable.
7. **Verify + docs + close-out.** Promote `diag_coupled_0{1,2}`; the operator-algebra theory
   page (the 2×2 coupled system, the four blocks, the seed-as-BC); full SN+numerics+transport
   wall; pyright ≤ 1; sphinx -W. Comment #280/#282; the DSA seam (#2) named as a consumer.

## ⏸ Compaction points

Any step boundary is valid; recommended defaults ⏸ after **step 3** (System B + both couplings
named — the algebra is posed, the hard assembly remains) and ⏸ after **step 5** (the coupled
operator + solve landed — only retirement + docs remain). Before each: commit everything +
append a `## CHECKPOINT` block (step statuses done@hash / next, deviations, rulings) + tell the
user it is safe to `/compact`. Re-anchor from this plan + `git log` (trust git).

## Standing constraints

- **Surgical mode** — main agent writes directly, user steers; NO `method-implementer`;
  `test-architect` (proactive — this is the numerics↔transport↔sn triple-boundary operator-
  algebra carve), `explorer`, `qa`, `elegance-enforcer`, `archivist` encouraged.
- **Tests** — canonical `.venv/bin/python -O -m pytest -p no:xdist --timeout=300
  -p no:cacheprovider` SERIAL (xdist unstable for tests/sn+tests/numerics).
- **pyright** — oracle is `npx pyright` / the ratchet CLI (baseline `transport:1` = accepted
  #288); streamed `<new-diagnostics>` = documented #226 LSP artifact, IGNORE.
- **Pushes HELD**; `.claude/agent-memory/*` + `scratch/*` + `.claude/skills/*` NOT the main
  agent's to commit; never `git checkout`/`restore` on uncommitted files (monkeypatch to
  revert mutations); commits stamp the SESSION model per [[feedback-truthful-model-attribution]]
  (**Opus 4.8** through 2026-07-10 B.1; **Fable 5** from the 2026-07-10 B.2 resume on — the
  user's /model switch); sphinx -W must stay clean.

## Out of scope / deferred AS CONSUMERS (not risk-deferrals)

- **DSA (#2)** — SN⊗diffusion via R/P; a CONSUMER of this coupled-block machinery once its
  restriction/prolongation operators exist (a separate feature). The seam is named at step 7.
- **RQI** — KKT/indefinite saddle-point (`SaddlePointOperator`, trigger #294 mixed-form) —
  different machinery.
- **Multiphysics** — nonlinear (Newton/JFNK) — different machinery.

---

## Invariant suite (test-architect blueprint, 2026-07-07) — "what invariant tests it?"

**Home:** NEW `tests/sn/operators/test_psi_half_coupling.py` (`pytest.mark.foundation`), reusing
anchor helpers from `test_starting_direction_metric.py` (`_recip_defect`, `_v_cell_seed`,
`_build_sphere`, `_composite`, `_dense`, `_blocks`) + `test_loss_transpose_solve.py`
(`_probe_augmented_matrix_one_group`). Extends `test_inverse_adjoint_coherence.py` (the swap
law on the coupled operator). **Carrying member = sphere-GL S4 ONLY**; cylinder/slab are the
non-carrying **control** (step 6). **≥2G mandatory** on every value row.

**Step 0 (LAND FIRST — the regression floor):** promote `derivations/diagnostics/
diag_coupled_0{1,2}_*.py` verbatim → `TestRegressionFloor` (6 tests: block-triangular #284
certificate `A_sb=0`; bulk→ray coupling lives in the lagged S-gain `S_sb≠0`/`S_bs=0`;
`ρ(M⁻¹N)≤c`; folded-seed nilpotent ρ=0; welded-sweep = exact inverse 3.5e-16; **extract =
principled-equiv vs dense-LU of assembled (L+C), 5.5e-16** — the oracle every EXTRACT step
pins against). Land BEFORE touching production.

**Per-step invariant classes** (each: defining law + in-process monkeypatch tooth + oracle):
- **Step 1 `TestA_BB_RadialBVP`:** march `solve∘apply=id`; ≡ closed-form ODE `φ=q/σ+(φ_R−q/σ)
  e^{−σ(R−r)}` O(Δr²) **on a GRADED mesh**; r=R Dirichlet data propagates (nonzero inflow); r=0
  pole-continuation `ψ½⁺(0)=ψ½⁻(0)`; reflective corner swap; fixed-source `Q/Σ_t`. WRAP→bit-id.
- **Step 2 `TestA_BA_SchurFold`:** the fold `Σ_ℓ (2ℓ+1)/2·(±1)^ℓ` on a **MANUFACTURED
  anisotropic ≥2-moment** input (iso is blind to `P_ℓ(±1)`); ℓ=0 ⇒ ½Q₀; A_BA ≡ the S arm AND
  the F arm (fissile `A/2g`) single-sourced (Mode-11 wrap-counter on the fold helper).
  WRAP→bit-id.
- **Step 3 `TestA_AB_SeedInjection`:** ≡ in-sweep injection; `A_AB.H`=`seed_cells_bar`
  (Euclidean transpose, per-group); `A_sb=0` triangular; seed-is-consumed asymmetry. WRAP→bit-id.
- **Step 4 `TestCoupledOperator`:** block matvec ≡ current apply (0-ULP view); **`.H` V_cell
  reciprocity stays Mode-12-closed** (the HIGHEST-value row — a Euclidean block adjoint on
  System B silently REOPENS Mode-12); swap law extends; `SystemRole` set-union join.
- **Step 5 `TestCoupledSolve`:** ≡ current fixed-source (`Q/Σ_t`, non-fissile) AND k_eff
  (fissile `A/2g`); EXTRACT ≡ dense-LU (principled); **ρ-honest block residual `r=Aψ−q` via
  `evaluate_residual` as the STOP test, not `‖Δψ‖`**; cold-residual lag-death classifier.
- **Step 6 `TestPresenceStructural`:** carrying-mesh-missing-seed → unconstructable; non-
  carrying-carrying-seed → unconstructable (`match=` the specific message, + positive control
  — the 7 guards have ZERO negative tests today, teeth NET-NEW); mixed-presence no-silent-drop.

**The 6 load-bearing refutations (fire before the ink dries):** (1) matvec≡current is bit-id
only for WRAP; EXTRACT is principled-equiv (split the gate). (2) block `.H` does NOT inherit
the metric automatically — a Euclidean block adjoint on System B reopens Mode-12 (ERR-067);
the V_cell-reciprocity gate on the assembled `.H` is mandatory. (3) the fold gate MUST
manufacture an anisotropic case (both S/F arms feed ℓ=0 only → iso blind to `P_ℓ(±1)`). (4)
non-fissile configs → k_inf is a `nan` dead gate; use fixed-source `Q/Σ_t` + a SEPARATE
fissile `A/2g` for the F arm. (5) presence guards have ZERO negative tests → net-new teeth,
`match=` the specific message (a downstream crash false-greens). (6) cylinder/slab are the
non-carrying CONTROL, not "other geometries."

**Mutation targets (monkeypatch-only, never `git checkout`):** `psi_half_angle_seed.py:148/150/
154`, `starting_direction_space.py:471` (fold sign), `loss_representation/__init__.py:2846`
(`_seed_rows_transpose`), `boundary.py:298` (`_reflect_starting_direction`), `{scattering,
fission}.py` (the `sd_out` arms).

---

## ⏸ CHECKPOINT — design phase complete (2026-07-07)

**Committed:** the FaceField C2 carve (`c637407`/`4081c0d`/`be5e7f8`), the C2-done design-note
checkpoint (`cc5a809`), this campaign plan (`e5713a3`) — branch `refactor/sn-walk-unification`,
**40 ahead of main, pushes HELD, tree clean.** The philosophy [[feedback-build-the-machinery-
operator-algebra]] + the #295 LayoutBearingSpace follow-up are saved. The 3 coupled-system
specialist memos persist in `.claude/agent-memory/{cross-domain-attacker,numerics-investigator}/`;
the numerics-investigator's `diag_coupled_0{1,2}` are on disk in `derivations/diagnostics/`.
The invariant suite above is the durable verification blueprint.

**NEXT (execution phase):** Step 0 — promote `diag_coupled_0{1,2}` → `test_psi_half_coupling.py::
TestRegressionFloor` (the floor, land first) → then Step 1 (pose System B: `RadialCharacteristic
Flux` rename + the 2-BC `RayOp`). Re-anchor after `/compact` from THIS plan + `git log`.

---

## Step 1 design — SETTLED + EXECUTING (2026-07-07)

**Step 0 DONE @ `2834b90`** (regression floor 6/6 green; the 2 scratch diagnostics retired).
**Step 1 IN PROGRESS.** Specialist memos: explorer `a17685c8…` (march + rename map),
cross-domain-attacker `ad3b7c3f…` (composition frame + precedent), test-architect `a5e2bdff…`
(12-gate suite + mutation ledger).

**Decomposition — 3 bit-identical commits, RENAME FIRST ("clean before extending": name the
substrate, then extend it with the new operators):**

- **1a — Rename `StartingDirection → RadialCharacteristic`.** 550 refs / 24 prod + 27 test + 17
  doc. Atomicity-pinned on the `FullField.starting_direction` dataclass field (+39 constructor
  sites) + `FullFieldSpace`. 7 CamelCase symbols (`…Field/Flux/SourceSink/Residual/Displacement/
  Space` + bare `StartingDirection` role-grid) + snake `starting_direction*`. 5 prod module
  `git mv` + 2 test-module `git mv`. **18 silent-render doc `:class:`/`:func:` cross-refs**
  (`discrete_ordinates.rst`, `loss_representations.rst`, …) — `-W` does NOT catch these; rename
  in lockstep + grep-verify. Space-separated physics prose ("starting direction" μ=±1, Hébert
  citations) is NOT renamed (sed hits only underscore/CamelCase symbols). Mechanical, bit-id.
- **1b — Un-weld `SNBoundaryOperator → B_a + B_b`.** B_a = the trace arm (`_reflect_trace`),
  B_b = the ray-corner arm (`_reflect_starting_direction` extracted into its own operator).
  `B = B_a + B_b` via `OperatorSum`, bit-identical (disjoint composite blocks). **`split()`
  relocates onto B_a** (B_b schedule-atomic — closes the LATENT BUG below); **per-leaf
  `is_adjointable`** (B_b: reflective/vacuum → yes, white/albedo/periodic → no); **exact-zero
  sibling blocks** (NOT `None` — the d1 mixed-presence law raises under `FullField.__add__`).
  Gates: `TestBoundaryUnweld` (per-block byte-id, vacuum + reflective), `TestB_b_RayBoundary`
  (+ the G_sd reciprocity gate), `TestSplitInteraction`.
- **1c — Pose `RayOp = A_BB`.** The radial straight-characteristic march — a WRAP of the
  standalone `carlson_inward_sweep_from_source` (+ `…_transpose`) at `psi_half_angle_seed.py:
  87-207`; production sites `loss_representation/__init__.py:4092-4117` (solve, 2 legs) +
  `:4614-4643` (transpose). Reads r=R corner as GIVEN data; r=0 pole internal (inward-exit-face
  → outward-entry-face leg-chaining). Gates: `TestA_BB_RadialBVP` (foundation: solve∘apply=id,
  WRAP bit-id via a Mode-11 call-counter sentinel, pole-thread, Dirichlet propagation, fixed-
  source Q/Σ) + the L1 closed-form-ODE convergence gate `φ=q/σ+(φ_R−q/σ)e^{−σ(R−r)}` on a GRADED
  mesh (measured orders 1.71→1.85→1.92) in a sibling `tests/sn/operators/test_ray_operator.py`
  (`pytest.mark.l1` — L9: don't conflate foundation + verifies). ⚠ `RayOp.apply` (banded matvec)
  vs the composite `A_ss` block is the ONE non-bit-id Step-1 spot — gate at `rtol`, not
  `array_equal`, unless it calls literally the same code.

### RULING P1 — the boundary-composition precedent (cross-domain-attacker, DURABLE)

A block-composed system's boundary is the **direct sum of per-system boundary blocks over the
composite biproduct carrier**, `B = Σ_x ι_x∘B_x∘π_x`. Each `B_x` is produced by a single
**carrier-indexed realization** of the one physical face law (`realize(law, carrier)`), never
hand-coded per carrier. The **off-diagonal structure is keyed to the face physics, not fixed**:
reflection/albedo/vacuum are per-system endomorphisms ⇒ block-diagonal `diag(B_a, B_b)`; a face
that **transmits/converts** one system's outflow into another's inflow (CP/MoC interface
currents, Schwarz transmission, mode-coupling albedo) is an off-diagonal `B_ab` ⇒ **never assume
diagonal**. **DSA is block-diagonal** (its SN↔diffusion coupling is interior volumetric R/P, not
boundary) — so this precedent + the carrier-indexed realizer pays for DSA directly. **Corollary
— grading:** any solver-internal grading of a boundary (Gauss-Seidel schedule, octant, group)
is a refinement of exactly ONE system's block `B_x` and MUST live on that leaf, never on the
composite `Σ_x B_x` (a grading on the composite silently replicates every un-graded sibling).

### RULING P2 — B_b.H is Mode-12-SAFE: Euclidean = Hilbert (test-architect, measured)

The ray corner gauge is **symmetric** (`g₊ = g₋ = V(R)` — both corners at r=R; measured
130.82 exact), so `B_b.H = G_sd⁻¹ B_bᵀ G_sd = B_bᵀ`. B_b advertises the **Euclidean
`apply_transpose` ONLY** — no per-leaf `.H` (`.H` is realized once at the composite via
`G⁺·apply_transpose·G`, L19). It is a metric PROPERTY, not a triviality → land a **B_b-isolated
G_sd-reciprocity gate** (control = 0 + two teeth: transpose-flip → 0.96 RED, corner-gauge
asymmetry → 0.33 RED). This CLOSES the §16.A A4 honest-scope gap ("G3 cannot see this arm") the
welded corner arm always carried.

### LATENT BUG (closed by 1b, not a separate fix)

`split()` doubles the ray corner: `SNMaskedBoundaryOperator.apply → _apply_faces` restricts only
`_reflect_trace(rows=…)` but calls `_reflect_starting_direction` UNRESTRICTED, so
`B_lower.apply + B_upper.apply ≠ B.apply` on System B's block (corner emitted by BOTH halves).
Latent today — `split()`'s only consumer (`ScheduledInvertibleOperator`) is seedless multi-D
Cartesian G-S; the seed-carrying sphere runs the Jacobi degenerate (`B_lower=0`). Goes LIVE on
curvilinear multi-D (a mesh with BOTH a ψ½ seed AND an octant schedule — **#22**, this
campaign's path). The un-weld fixes it structurally (`split()`→B_a; B_b schedule-atomic);
`TestSplitInteraction` pins it.

### Step-1 EXECUTION STATUS (2026-07-07)

- **Step 0 DONE @ `2834b90`** — regression floor 6/6.
- **1a DONE @ `f9e4c4a`** — the `StartingDirection → RadialCharacteristic` rename (550 refs /
  55 files, bit-identical; plan capture @ `b015e36`, 1b map @ `edd2954`).
- **1b DONE @ `604c890`** — the boundary un-weld `B → B_a + B_b`. B_b =
  `RadialCharacteristicBoundaryOperator` (System B ray boundary); the latent split-double-count
  closed structurally; `_within_group_triple` returns the composite; the reconstruction reflects
  per-system. Bit-identical (ratchet transport:1, floor byte-identical, 932+589 passed, sphinx
  -W). NEW gates `TestBoundaryUnweld`/`TestB_b_RayBoundary`/`TestSplitInteraction` with
  mutation-proven orthogonal teeth (reciprocity control=0 + wrong-transpose=1.00 /
  gauge-asymmetry=0.33). elegance-enforcer PASS (no violations; all fixes applied). 2 tests
  rewired to the successor API (reflective-sphere Q/Σ recovery + SI/eig single-primitive).
- **1c DONE @ `77e7b72` (2026-07-07) — `RadialCharacteristicOperator` built + 15 gates
  landed, all green + mutation-verified. ⟹ STEP 1 COMPLETE (1a rename + 1b boundary
  un-weld + 1c RayOp).** Verify: ratchet `transport:1`, full `tests/sn/operators` 719
  passed, sphinx -W exit 0, elegance-enforcer PASS. NEW
  `orpheus/sn/operators/radial_characteristic.py`:
  `RadialCharacteristicOperator(sn_mesh, total_cross_section)`, `domain=codomain=
  radial_characteristic_space`. Realizes the **resolvent action** `.solve` (the two-leg
  Carlson march WRAP, keyed by `space.levels` = the p_idx positions) + its **adjoint**
  `.solve_transpose` (the ISOLATED reverse march — the in-sweep reverse MINUS the
  `psi_angle_bar` bulk-coupling term, which is A_AB's transpose, steps 2–3). Smoke: exact
  Euclidean adjoint `⟨solve(u),v⟩=⟨u,solve_transpose(v)⟩` at **3.3e-16**; WRAP bit-identical
  to direct engine calls. **DESIGN DECISIONS LOCKED (elegance-review + user may overrule):**
  (1) forward `.apply` (μ∂_r+σ_t matvec) **loud-deferred to step 4** — it is woven into
  `(L+C).apply` (R13); a standalone spelling now = a twin or a premature hot-path carve.
  (2) `is_invertible=False` + `is_adjointable=False` in 1c (base semantics: `is_invertible`
  advertises the `inverse()` OPERATOR-returning act, which lands in step 4 WITH the forward
  `.apply` — the involution web `inverse().solve == forward apply` needs both directions;
  shipping `inverse()` now = a resolvent whose `solve` raises = a false capability). The
  math-invertibility is expressed by the working `.solve`. (3) `block_role=None` (the
  SystemRole lattice is step 4). (4) **σ_t is a typed `CrossSectionField`** (USER
  RULING 2026-07-07, elegance-review-flagged — deviates from the plan's earlier
  bare-ndarray ratification): mesh-bound, so the constructor guards mesh-identity
  (closes the Pattern-4 illegal-state — a σ_t from a foreign mesh would silently
  march the wrong Δr; a bare ndarray couldn't catch it), aligns A_BB with the
  inverse-family sibling `InvertibleOperator`'s typed coefficient, and reads like
  `C_ray = M[σ_t]`. The field's own space invariant makes the (ng,nx) shape correct
  by construction (the explicit shape check is retired as redundant). Step 4 sources
  it as `mat_xs.total_cross_section_field`. Gates: `TestA_BB_RadialBVP` foundation (9, in
  `test_psi_half_coupling.py`: adjoint-consistency+tooth, Mode-11 WRAP sentinel, pole
  continuation, 2.5a reversal on graded, Dirichlet propagation+tooth, Q/Σ equilibrium ≥2G,
  constructor CONTROL cylinder/slab/malformed) + L1 convergence (6, NEW
  `tests/sn/operators/test_ray_operator.py`, `pytest.mark.l1`: closed-form ODE
  `φ=q/σ+(φ_R−q/σ)e^{−σ(R−r)}` uniform+graded, orders 1.82→1.95; 4 DD-coefficient +
  Mode-5-index-drift teeth). **⚠ TRANSIENT TWIN (Cardinal Rule 2, TRACKED):** `.solve`'s
  two-leg orchestration mirrors the in-sweep block `loss_representation:4104-4119` (the
  engine is single-sourced; the orchestration is not yet). Retired at step 4/5 when the
  coupled resolvent routes the production `(L+C)` ray solve through this operator — see the
  RETIREMENT LIST below. Both sides behaviour-pinned meanwhile (regression floor + sweep
  suite ⟂ the Mode-11 WRAP gate). Recon (the durable build spec) retained below.

  **Home + name:** NEW `orpheus/sn/operators/radial_characteristic.py`; class
  `RadialCharacteristicOperator` (A_BB interior; the boundary sibling `RadialCharacteristic
  BoundaryOperator` = B_b landed in 1b). `domain = codomain = sn_mesh.radial_characteristic_space`.

  **The WRAP source (verbatim, lift into `RayOp.solve`):** the two-leg march + pole-chaining is
  `loss_representation/__init__.py:4092-4117` — read `q_minus=cells_view(src,·,-1)`,
  `q_plus=cells_view(src,·,+1)`, `corner_in=corner_view(src,·,-1)`; inward
  `cells_minus,pole_face = carlson_inward_sweep_from_source(q_minus, σ_t, dr, corner_in)`;
  outward (reversed data) `cells_plus_rev,corner_out = carlson_inward_sweep_from_source(
  q_plus[:,::-1], σ_t[:,::-1], dr[::-1], pole_face)`; write cells(-1)=cells_minus,
  corner(-1)=corner_in, cells(+1)=cells_plus_rev[:,::-1], corner(+1)=corner_out. The engine
  `carlson_inward_sweep_from_source` (+ transpose `…_transpose`) is `psi_half_angle_seed.py:
  87-207` (already tested in `test_psi_half_angle_seed.py`; Hébert 3.434/3.435, `:label:
  hebert-3-434`). RayOp WRAPS it — single source, does NOT re-implement.

  **σ_t / dr sourcing (confirmed):** σ_t = a CONSTRUCTOR PARAM (ng,nx) — the C_ray collision
  coefficient, matching the sweep's `sig_t` param + `MultiplicationOperator(coefficient)`; the
  test provides it, Step 4 provides `solver.mat_xs.total_cross_section` (shape (ng,nx) verified).
  dr = `sn_mesh.axis_widths[0]` (nx,) — mesh geometry. So `RadialCharacteristicOperator(sn_mesh,
  total_cross_section)`.

  **✅ P-IDX-VS-LEVEL — RESOLVED @ `bf07922` (numerics-investigator verdict: NOT A BUG; key by
  `p_idx`).** `RadialCharacteristicSpace.levels` are level POSITIONS (`enumerate(raw)`), NOT
  sparse μ-indices; the gate `p_idx in seed_levels` and the space key share ONE source
  (`mesh.radial_characteristic_levels`), so a passing `p_idx` is ALWAYS a valid slot key — no
  crash, no wrong-slot, for any config. Multi-level-sparse carrying meshes DON'T occur (sphere
  always 1 level `(0,)`; cylinder never carries; verified 16-quadrature battery). **RayOp MUST
  key by `p_idx` (level POSITION), NOT `level`** — `level` is the walk's `mu_level_idx`, `None`
  on the sphere (the only carrying geometry), so passing it would raise. The invariant is pinned
  by `tests/sn/mesh/test_radial_characteristic_slot_coordination.py` (34 tests) + an inoculating
  comment at the solve site.

  **Forward `.apply` (banded matvec μ∂_r+σ_t) — DEFER to Step 4 as a documented seam:** the
  forward action of A_BB is currently woven into `(L+C).apply` (the ray residual on the +1
  outflow corner, R13); extracting it standalone is the Step-4 CoupledOperator's job (it needs
  the (B,B) matvec block). RayOp in 1c realizes the RESOLVENT (`.solve` = the march) + its
  adjoint (via `carlson_inward_sweep_transpose`). For the `solve∘apply=id` invariant, use the
  test-architect's fallback: the dense `(L+C).apply` restricted to the ray block as the forward
  reference (NOT a from-scratch matvec).

  **Gates** (test-architect design, home = sibling `tests/sn/operators/test_ray_operator.py`
  `pytest.mark.l1` for the convergence gate; foundation gates extend `test_psi_half_coupling.py`):
  `TestA_BB_RadialBVP` — the L1 closed-form-ODE convergence `φ=q/σ+(φ_R−q/σ)e^{−σ(R−r)}` on a
  GRADED mesh (measured orders 1.71→1.85→1.92, already validated; mutations at `psi_half_angle_
  seed.py:150/151/154` + `dr[::-1]` graded-vs-uniform Mode-5 control); WRAP bit-id vs the
  in-sweep `:4092-4117` (Mode-11 call-counter sentinel on `carlson_inward_sweep_from_source`);
  pole continuation `ψ½⁺(0)=ψ½⁻(0)` (sentinel-captured leg args); r=R Dirichlet propagation
  (φ_R=0 vs 1 → interior differs by the exp envelope); fixed-source Q/Σ equilibrium; ≥2G. Build
  q½ DIRECTLY (uniform), NOT via `fold_moments_to_radial_characteristic` (that fold is Step 2).

### RETIREMENT LIST (tracked debt — retired at the named step, "retire as you go")

Enumerated so no superseded path falls silently out of view (coding-standards
"retire as you go"). Each entry: what, `file:line`, retired-at, and the guard holding it
until then.

1. **The in-sweep two-leg ψ½ orchestration** — `loss_representation/__init__.py:4104-4119`
   (the direct-seed march inside `(L+C).solve`; reads/marches/writes the ray blocks). A
   transient twin of `RadialCharacteristicOperator.solve` (the engine is single-sourced; the
   orchestration is duplicated). **Retire at step 4/5** — when the coupled block-triangular
   resolvent routes the production `(L+C)` ray solve through `A_BB` (the CoupledOperator's
   `(B,B)` block), the inline block collapses to that call (the `phi_aux = cells_minus` M-M
   seed use at `:4123` reads the operator's returned inward cells). Guard until then: the
   in-sweep by the regression floor + `tests/sn/sweep/**`; `A_BB.solve` by the Mode-11 WRAP
   bit-identity gate. (Step 5's `EXTRACT ≡ dense-LU` gate is the equivalence that licenses
   the collapse.)
2. **The in-sweep reverse ψ½ orchestration** — `loss_representation/__init__.py:4621-4649`
   (inside `sweep_transpose`; the reverse march + the `psi_angle_bar` bulk term). Its pure
   ray-block half is `A_BB.solve_transpose`; the `psi_angle_bar` term is `A_BA`'s transpose
   (steps 2–3). Retire alongside entry 1 (the coupled adjoint-solve routes through
   `A_BB.solve_transpose` + the named `A_AB`/`A_BA` transposes).
3. **The in-sweep A_AB forward placement** — `loss_representation/__init__.py:3168-3186`
   (inside `(L+C).apply` / `_OneDimScanWalk._apply_walk`; the seed's `angular_numer` term
   `−(ΔA/w)·c_in·ψ_{m-1/2}/V` fused into `m_full = (denom·ψ − numer)/V`). A transient twin
   of `RadialCharacteristicSeeding.apply` (the M-M *kernel* `precompute_psi_state` /
   `cell_contribution` is single-sourced; only the `∓numer/V` orchestration is duplicated —
   a DIFFERENT production entry point from entries 1–2, which are the resolvent/solve side).
   **Retire at step 4/5** — when the CoupledOperator routes the `(L+C)` bulk rows through
   `A_AA + A_AB` (the block matvec), the inline seed placement collapses to `A_AB.apply`.
   Guard until then: the in-sweep by the regression floor + `tests/sn/sweep/**`; `A_AB.apply`
   by the `TestA_AB_SeedInjection` bit-identity gate + Mode-11 WRAP sentinel.
4. **The in-sweep A_AB transpose placement** — `loss_representation/__init__.py:~3475-3590`
   (inside `(L+C).apply_transpose`; the `numer_bar = −ob/V` gather at `:3498` → `angular_adjoint`
   at `:3576` → `seed_cells_bar` landed on `cells(p,-1)` via `_seed_rows_transpose:2846`). Its
   isolated projection is `RadialCharacteristicSeeding.apply_transpose`. Retire alongside
   entry 3 (the coupled adjoint block matvec routes the bulk cotangent → ray through the named
   `A_AB.apply_transpose`).

### OPEN DECISION D1 — the B_a naming asymmetry (elegance-enforcer, user's call)

`SNBoundaryOperator` (B_a) keeps the generic name while `RadialCharacteristicBoundaryOperator`
(B_b) is descriptive — so B_a reads as *primary*, B_b as *auxiliary*, though the algebra frames
them co-equal (System A / System B boundaries). Honest today (B_b IS subsidiary — block-
triangular in the ray, lagged-scatter coupled), but when DSA adds a 3rd sibling
`DSABoundaryOperator`, "which is THE SN boundary?" becomes concrete. Two options: rename B_a →
`SNTraceBoundaryOperator` (symmetric, greppable — but a wide rename, and "trace" slightly
undersells the |Ω·n|·w metric it carries), OR keep the asymmetry (the docstrings already say
"B_a — System A's (trace) boundary", honest). **Current call: KEEP `SNBoundaryOperator`** — the
rename is premature (DSA forces the symmetry question with a concrete 3rd instance; defer-until-
the-shape-is-clear). Revisit when DSA lands. User may overrule.

### SEAM — carrier-indexed realizer (DEFERRED, user-ruled 2026-07-07)

The precedent's *production* half — `realize(law, carrier)` producing B_a (trace) + B_b (ray
corner) from one face law, extending the existing `SNBoundaryRealizer` to be carrier-aware.
**Built when DSA's diffusion-albedo adds the 3rd carrier** (≥3 makes `realize(law, carrier)` pay;
today the ray realization is a ~2-line specular-swap kind-fact). B_b relocates the hand-coded
kind-dispatch for now. NOT a risk-defer — the abstraction genuinely doesn't yet pay at 2
carriers; the RULING (P1) captures the shape so DSA builds it.

### 1b implementation map — SCOPED, ready to execute (2026-07-07, `f9e4c4a`)

The variadic-gains driver makes this a clean slot-in (B_b is a lagged gain, exactly as B was
separated from S in Wave O #208). Consumption audit (post-rename): the boundary op is consumed
three ways — (1) as a `.apply` gain (matvec `−Σgᵢ.apply` + rhs `+Σgᵢ.apply`); (2) `.split()`
ONLY in the seedless multi-D-Cartesian G-S branch (`_select_si_resolvent`, `solver.py:558-564`
— `sn_mesh.is_cartesian and not is_1d`, never a seed-carrying mesh); (3) `reflect_inflow_inplace`
at the reconstruction (`solver.py:1719`) on a FRESH `SNBoundaryOperator(sn_mesh)`, decoupled from
the triple.

**boundary.py:**
- `SNBoundaryOperator` STAYS = **B_a** (System A trace boundary). Drop `_reflect_radial_
  characteristic` from `_apply_faces` (`:293-295`); B_a emits a **present-but-ZERO** ray corner
  (`RadialCharacteristicSourceSink.zeros_on(mesh)` when the input carries a seed, `None` when
  not) — NOT bare `None` on a seed-carrying input (the d1 mixed-presence law raises under
  `FullField.__add__`; cross-domain-attacker refutation). Remove the `radial_characteristic`
  kwarg from `reflect_inflow_inplace` (`:401-448`). `split()`, `reflect_into_inflow`,
  `is_adjointable` (already trace-only) unchanged.
- NEW `RadialCharacteristicBoundaryOperator` = **B_b** (System B ray boundary): `apply` = zero
  bulk + **present-zero trace** (`AngularBoundarySourceSink.zeros_on`) + `_reflect_radial_
  characteristic(seed)` (the method MOVES here); `apply_transpose` = the Euclidean mirror;
  `is_adjointable` per-leaf (reflective/vacuum → True, white/albedo/periodic → False);
  `domain=codomain=full_field_space`; `block_role=BOUNDARY`. Add `reflect_corner_inplace(seed)`
  (the in-place ray seeding the reconstruction needs).
- `B = B_a + B_b` via `OperatorSum` reproduces today's welded `SNBoundaryOperator.apply`
  BIT-IDENTICALLY (disjoint present-zero blocks sum exactly).

**solver.py:**
- `_within_group_triple` (`:224-228`): `B_a = SNBoundaryOperator(sn_mesh)`; return
  `B = B_a + RadialCharacteristicBoundaryOperator(sn_mesh)` iff
  `sn_mesh.radial_characteristic_space is not None`, else `B_a`. Widen the return annotation to
  `SNBoundaryOperator | OperatorSum` (the split-only-seedless invariant guarantees `.split()`
  sees B_a; seed-carrying ⟹ Jacobi ⟹ B is the composite gain, never split).
- Reconstruction ray seeding (`:1719`): split into `B_a.reflect_inflow_inplace(bf)` (trace) +
  `B_b.reflect_corner_inplace(seed)` (ray).

**LATENT BUG auto-fixed:** B_a's masked `split()` halves now emit ZERO ray corner (the arm is
gone from `_apply_faces`), so `B_lower + B_upper + B_b = B` with no ray double-count.

**Gates (extend `test_psi_half_coupling.py`):** `TestBoundaryUnweld` (per-block byte-id vacuum +
reflective), `TestB_b_RayBoundary` (corner swap ± transpose, vacuum-zero, white/albedo
`NotImplementedError` `match=`, block-locality, **the G_sd-reciprocity gate** control=0 + 2 teeth
0.96/0.33), `TestSplitInteraction` (split lives on B_a; B_b unaffected). Mutation ledger + the
measured teeth are in the test-architect memo `a5e2bdff`.

---

## Step 2 — EXECUTION STATUS (DONE, 2026-07-08)

**DONE — the ψ½ Schur fold un-welded into ONE named operator.** All green + mutation-verified +
elegance PASS; bit-identical (731 `tests/sn/operators` pass unchanged incl. the curvilinear
aniso-P1 MMS end-to-end fold), ratchet **`transport:1`**, sphinx **-W exit 0**.

**The realization.** `A_BA = Reconstruction ∘ Emission`; the shared factor is the **Reconstruction**
(the fold), NOT `A_BA` itself (the emission — `K_iso·φ₀` for S, `χνΣf·φ/k` for F — is
operator-specific; the FULL A_BA block assembles at step 4). The user rejected the initial
`RadialCharacteristicFold` name (it conflated the fold FACTOR with the A_BA block) — the fold IS
the angular reconstruction `R` sampled at the closed μ=±1 rays, so:

- **NEW `orpheus/transport/operators/radial_characteristic_reconstruction.py`** — class
  `RadialCharacteristicReconstruction(LinearOperator)` (name USER-chosen 2026-07-08 over
  `…Emission`/`…Restriction`; "Restriction" would collide with the DSA multigrid R, #2).
  `apply(moments (n_moments,ng,nx)) → RadialCharacteristicSourceSink` (broadcast fold onto every
  carried level); `apply_transpose(ray cotangent) → (n_moments,ng,nx)` (the Euclidean transpose,
  sum over levels/signs). `is_adjointable=True`, `codomain=ray space`, `domain=None` (untyped
  bulk-moment intermediate, K_iso-precedent), `n_moments` FIXED at construction (Pattern 4 — the
  apply/apply_transpose adjoint pair agrees on the moment dim by construction; default 1 = iso
  production reach). `block_role=None` (SystemRole lattice is step 4).
- **HOME = transport — but a step-2 TRANSIENT, NOT native (USER RULING 2026-07-08).** The fold's
  *native* home is beside A_BB in `sn/operators/radial_characteristic.py` (both are blocks of the
  ONE ψ½ coupled operator). It sits in transport at step 2 ONLY because S/F (transport) still
  CONSUME it (transport ⊥ sn at runtime — grep-verified: the only `from orpheus.sn` in transport/
  are `TYPE_CHECKING`). **That consumption IS the monolithic posing this campaign un-welds** — the
  ψ½ starting-direction ray is a curvilinear-SN augmentation a model-generic gain should not own;
  S/F owning a ray arm is the welded coupling. **Step 4 reverses the arrow** (see the lift
  deliverable below). The earlier "native, dependency-forced split" framing was the main agent
  going along with the inherited posing — the user caught it; the reconstruction's module docstring
  now labels the transport home as a transient.
- **Single-sourced the weights** (Pattern 2/7): `_radial_characteristic_reconstruction_weights(n,
  sign) = (2ℓ+1)/2·(±1)^ℓ` in `radial_characteristic_space.py`, shared by the forward
  `fold_moments_to_radial_characteristic` AND the NEW `fold_moments_to_radial_characteristic_transpose`
  — the P₁(−1) sign spelled ONCE; the **S-adjoint's hand-rolled `0.5` is RETIRED** to it (the
  ERR-022 forward/adjoint-coefficient-drift class dissolved by construction).
- **Routed** S-fwd + F-fwd (`.apply(emission[None])`) + **S-adjoint** (`.apply_transpose(chi_seed)`
  → `K_isoᵀ(m_bar[0])`); the 3 hand-rolled fold loops **DELETED** (aggressive retirement, complete
  — no `sd_values`/`cells_bar_sum` tells survive). Stale `"dormant until d3"` comments fixed (d3
  landed @ `a29ab2d`; the arms are LIVE).
- **`from_angular_source` NOT routed through** (the one deviation from the AskUserQuestion preview,
  user-flagged): it is the per-ordinate birth factory doing per-level Legendre **analysis** (a
  DISTINCT operation, not the S/F broadcast-fold twin); it already shares the fold math (helper +
  weights). Forcing it through the broadcast operator would be correct only via the R12a 1-level
  coincidence (fragile). elegance-enforcer ACCEPT: different physics, summed into the same q½ block.

**Gates** (`TestA_BA_SchurFold`, `test_psi_half_coupling.py`, extended by test-architect memo
`a651e6e9`): the fold contract on a MANUFACTURED anisotropic ≥2-moment input (refutation #3), the
Mode-11 wrap-counter (S+F both route through the reconstruction — re-pointed to `_rcr_mod`, the
reconstruction module's fold binding, after the module-level import moved the patch target), the
½·emission bit-identity vs the old loop, the fold-transpose Euclidean contract, the S fwd↔adj
consistency, the non-carrying cyl/slab CONTROL. The 3 formerly-xfail direct-surface gates are now
LIVE (BIND flipped to `RadialCharacteristicReconstruction.apply/.apply_transpose`; xfail markers
removed). Teeth mutation-verified in-process (M1 weights→fold RED, M2 transpose shares weights, M3
bypass→count 0 RED).

**FOLLOW-UPS (out of scope — for the user to file / step 7):**
1. `from_angular_source = Reconstruction ∘ Analysis` — name the per-level Legendre `Analysis`
   companion (currently unnamed inline `legvander`+`einsum`; anti-#5/#6). `module:sn`/`type:improvement`.
2. **A_BA theory prose** — the operator-algebra theory page's A_BA section (the Schur fold, the
   `(2ℓ+1)/2·(±1)^ℓ` reconstruction, the broadcast/sum adjoint) → campaign **step 7** docs.

### Step-4 LIFT deliverable — reverse the ψ½ coupling arrow (USER RULING 2026-07-08)

The step-4 CoupledOperator assembly MUST include **the lift** that reverses the S/F→fold dependency
and gives the fold its native sn home (the completion of the un-weld step 2 began — step 2 named the
fold; step 4 lifts it out of S/F). Locked deliverables:

1. **S/F → pure bulk.** Drop the `(ray, bulk)` seed arms from `ScatteringOperator`/`FissionOperator`
   (scattering.py S-fwd + S-adj, fission.py F-fwd) — the model-generic gain STOPS producing ψ½ ray
   sources; its composite `apply` returns bulk (⊕ boundary) only.
2. **`A_BA = Fold ∘ K_iso ∘ integrate` (scattering) / `Fold ∘ FissionKernel` (fission)** — a NAMED
   sn coupled-block operator composing the SHARED transport emission kernels (`IsotropicScattering`
   = K_iso, already standalone; the fission dyad) with the fold. NO dependency on the full S/F.
3. **The SI driver applies `A_BA` as its own lagged gain** (the Wave-O #208 pattern that separated
   `B` from `S`): the ray seed is born driver-side, not inside S/F.
4. **Move the fold → `orpheus/sn/operators/radial_characteristic.py`** (beside A_BB). Once S/F no
   longer consume it, the runtime `transport → sn` ban lifts and the fold joins its self-block
   sibling; **retire `orpheus/transport/operators/radial_characteristic_reconstruction.py`**.
5. **The ψ½ DATA types STAY** at the transport/numerics data layer (`RadialCharacteristicSourceSink`,
   the carrier space, the residual) — generic transport ops (`M[σ]`, the residual) act on them; only
   the *operator* migrates. NO data-type re-homing.

*Invariant:* the driver's `S_bulk + A_BA` ≡ today's monolithic `S.apply` composite (bit-identical on
the sphere); the regression floor + the `TestA_BA_SchurFold` bit-identity gate pin it.

---

## Step 3 — EXECUTION STATUS (DONE, 2026-07-08)

**DONE — `A_AB` posed as ONE named operator, BOTH directions realized.** All green +
mutation-verified + elegance PASS-with-nits (all nits addressed); forward AND transpose
bit-identical to the in-sweep injection (0-ULP `array_equal`), Euclidean adjoint defect
**0.0**; ratchet **`transport:1`** (sn 0), sphinx **-W exit 0**, `tests/sn/operators` 735 +
`test_psi_half_coupling` 42, full `tests/sn -m "not slow"` 0 reds.

**The realization — CELL-LOCAL ANGULAR, both directions clean WRAPs (refines the pre-compaction
framing).** The pre-compaction analysis said A_AB "defers its forward like A_BB". The deeper
trace (explorer-verified, 5/5 claims ACCEPT) corrected this: the M-M recurrence
`ψ_half[m+1] = (ψ[m]−(1−τ)ψ_half[m])/τ` runs over ORDINATES at a FIXED cell, so **A_AB is
cell-local angular** — the seed at cell i couples ONLY to cell i's bulk ordinates, NO spatial
weaving (the radial march is A_BB's job). A_BB's defer-the-forward rationale (spatially-woven
matvec, no standalone engine) simply does NOT apply. So both directions realize as thin WRAPs
of the single-sourced closure methods (`precompute_psi_state`/`cell_contribution`/
`angular_adjoint` on `MorelMontryAngularSweep`), and `is_adjointable=True`.

- **NEW `RadialCharacteristicSeeding(LinearOperator["RadialCharacteristicField","AngularField"])`**
  in `orpheus/sn/operators/radial_characteristic.py` — BESIDE A_BB (native sn home; NO transport
  transient, NO lift — A_AB's consumer is the sn `(L+C)` walk, not model-generic S/F, so nothing
  forces it below sn). Name **USER-chosen 2026-07-08** (`…Seeding` over `…SeedInjection`/`…Injection`;
  "Injection" would collide with the DSA multigrid transfer term, the same concern that rejected
  "Restriction" for A_BA). Scope **USER-chosen: realize BOTH directions** ([[feedback-build-the-
  machinery-operator-algebra]] — the full honest posing, not the minimal defer-forward).
  - `apply(seed) → AngularSourceSink`: WRAP `precompute_psi_state(np.zeros((N,ng,nx)),
    radial_characteristic=seed)` (bulk ZEROED → isolates A_AB from A_AA's redistribution by
    LINEARITY) + per (carrying level, cell) `cell_contribution` → `−angular_numer/V`. Reads ONLY
    the inward `cells(p,-1)` leg.
  - `apply_transpose(cotangent) → RadialCharacteristicSourceSink`: the local `numer_bar[p]=−ob/V`
    gather (exact transpose of the forward `−·/V`) → `angular_adjoint` → `seed_cells_bar` on
    `cells(p,-1)`; discards `psi_ang_bar` (that's A_AA^T). `+1` leg + corners EXACTLY 0.
  - **σ-INDEPENDENT** — constructor takes `sn_mesh` ONLY (contrast A_BB's σ_t): with ψ_cell=0 the
    collision/streaming terms drop out, only the angular numerator survives (positively asserted
    in the gate: reference identical at two σ slopes).
- **The shared-kernel projection (step-4 seed):** `precompute_psi_state`/`cell_contribution`/
  `angular_adjoint` each serve BOTH A_AA (bulk redistribution) and A_AB (seed); A_AB projects out
  its block (zero the bulk / discard `psi_ang_bar`). Step-4 CoupledOperator calls the ONE kernel
  and routes into A_AA + A_AB — never a twin kernel.
- **Gates** (`TestA_AB_SeedInjection`, `test_psi_half_coupling.py`, test-architect memo
  `a_ab_seed_injection_verification.md`): forward bit-id + Mode-11 sentinel (precompute called
  once with all-zero psi_view + seed passed; cell_contribution ≥ nx) + σ-independence; transpose
  bit-id + zero-slots + Mode-11 (angular_adjoint once); **Euclidean adjoint consistency < 1e-11**
  (THE correctness gate — compares apply↔transpose, catches a shared-method bug the L13 bit-id
  gates are blind to) + committed tooth (cell_contribution sign flip → defect O(1)); reads-only-
  inward-leg asymmetry; non-carrying cyl/slab CONTROL + mesh-identity. Teeth all mutation-verified
  RED in-process (cell_contribution flip, angular_adjoint seed_bar flip, precompute-reads-+1).
- **Retirement entries 3-4 ADDED** (the A_AB in-sweep matvec placement twins at
  `loss_representation:3168-3186` fwd + `~3475-3590` transpose — a DIFFERENT entry point from the
  A_BB resolvent twins 1-2; retire at step 4/5 when the CoupledOperator routes the bulk rows
  through A_AA+A_AB). Docstring carries a symmetric "Tracked transient twin" note.
- **Docs (Cardinal Rule 3):** `.. automodule:: orpheus.sn.operators.radial_characteristic` added
  to `docs/api/discrete_ordinates.rst` (was MISSING since 1c — A_BB + A_AB docstrings rendered
  nowhere); this SURFACED + FIXED a latent RST bug in the A_BB module docstring (a paragraph split
  a bullet list). `sn/operators/__init__.py` inventory + V&V matrix regenerated.

**FOLLOW-UP (step 4 — the block-operator theory-page PROSE for A_AB, absent as it was for A_BA;
natural home = the CoupledOperator assembly docs.)**

## ⏸ CHECKPOINT — step 3 complete (2026-07-08)

**Committed on branch `refactor/sn-walk-unification`** (pushes HELD): Step 3 = `RadialCharacteristic
Seeding` (A_AB) + `TestA_AB_SeedInjection` + docs. Steps 0/1(1a/1b/1c)/2/3 DONE; the algebra is
POSED (System B + all THREE couplings A_BA/A_BB/A_AB named + realized). **Deviation from this
plan's step-3 line:** the plan said "matches the in-sweep injection; the adjoint seed_cells_bar is
A_AB.H" — realized as intended, PLUS the forward `.apply` is ALSO realized (the cell-local
discovery made it a clean WRAP, not a defer) — a STRENGTHENING, user-ratified. **Rulings:** name
`…Seeding` (collision-avoidance); realize-both (build-the-machinery); native sn home (consumer-
layer forces the home, not the data types — A_AB's consumer is sn, unlike A_BA's transport S/F).

**NEXT (step 4 — the HARD one): assemble `CoupledOperator [[A_AA,A_AB],[A_BA,A_BB]]`** + widen the
role lattice (`SystemRole {A,B,COUPLED}`) + **THE LIFT** (S/F→pure bulk, A_BA=Fold∘K_iso by the
driver, fold→sn) + **THE WALK UN-WEAVE** (retirement entries 1-4: extract A_AB/A_BB forward matvecs
+ the solve-side orchestrations from the interleaved `(L+C)` walk into the explicit block matvec).
Highest-risk step (touches the interleaved walk) — gate hard, WRAP-first then extract. Re-anchor
from THIS plan + `git log` (trust git). ⏸ COMPACTION recommended NOW (after step 3).

---

## Step 4 — DECOMPOSITION + EXECUTION (the HARD step; WRAP-first ruling, 2026-07-08)

**Ground truth (explorer + test-architect dispatches, 2026-07-08).** The KEY finding: **all four
extraction-target operators already exist + are tested** (A_BB.solve/solve_transpose @ 1c;
A_AB.apply/apply_transpose @ step 3). So the walk un-weave is *route the production walk through the
built operator, then DELETE the inline twin* — NOT "build the operator." The ONE genuinely-unbuilt
piece is **A_BB.apply** (the forward μ∂_r+σ_t matvec, deferred from 1c) — the CoupledOperator's (B,B)
block + `inverse()` need it. Memos: explorer `ac4355b6…` (blast-radius map, current `file:line`),
test-architect `.claude/agent-memory/test-architect/coupled_operator_step4_verification.md` (the
extended gate spec), elegance-enforcer `.claude/agent-memory/elegance-enforcer/coupled_block_boundary_
unweld_rulings.md`.

**Verification bar SPLITS three ways** (one gate can't serve all three — refutation #1 operational):
4d WRAP = `array_equal` 0-ULP vs current augmented apply (**TRANSIENT** — retired at 4e); 4c LIFT =
bit-id vs captured baseline; 4e EXTRACT = principled-equiv **5.5e-16** vs the dense-LU-of-assembled-
(L+C) oracle (regression floor row 6). The DURABLE gate is the principled A1-EXTRACT (green at 0-ULP
through 4a-4d, stays green at 5.5e-16 after 4e — never churns; the WRAP bit-id proof is a separate
transient characterization, honestly retired at 4e). Load-bearing row: **A2 — the assembled `.H`
V_cell reciprocity stays Mode-12-closed** (a Euclidean block adjoint on System B reopens ERR-067;
reds on ALL geometries — the "slab stays green" instinct is FALSE, burned on 2.5c).

**The 5 risk-ordered sub-commits** (WRAP-first, user-ruled 2026-07-08):
- **4a — SystemRole {A,B,COUPLED} lattice.** DONE @ `4113d27` (see below).
- **4b — A_BB.apply forward + inverse() + is_adjointable/is_invertible flips.** The one unbuilt piece.
  ⚠ predicate-flip contract blast (test-architect flag 2): grep EVERY `is_adjointable is False` /
  `is_invertible is False` for `RadialCharacteristicOperator` (`test_capability_survival.py`,
  `test_inverse_operator_equivalence.py`, the contract table) + update in the SAME landing. Gate:
  forward matvec ≡ dense `(L+C).apply` ray-block; `solve∘apply=id` both directions.
- **4c — THE LIFT.** Move the fold `RadialCharacteristicReconstruction` transport→sn (beside A_BB);
  build `A_BA = Fold∘K_iso∘integrate` (S) / `Fold∘FissionKernel` (F); drop the ray arms — S-fwd
  (`scattering.py:1572-1597`) + S-adj (`:1806-1854`, **ASYMMETRIC** — ray-out present-ZERO, bulk-out
  gets `K_isoᵀ∘Reconstructionᵀ` pullback; a symmetric "drop" LOSES the bulk term) + F-fwd
  (`fission.py:493-534`) → pure bulk. **HAZARD 3: the LIFT + driver-gain-add MUST land in ONE commit**
  (presence law — S/F dropping the ray while the iterate still carries the 3rd block trips
  `FullField.__sub__` mixed-presence). **HAZARD 5: S and F lift onto DIFFERENT seams** — S = within-
  group lagged gain (the `(S,B)` tuple); F = the OUTER `q_ext` (within-group fission=0; the F ray
  source rides `compute_fission_source`, `solver.py:951`, NOT the gain tuple). Retire
  `transport/operators/radial_characteristic_reconstruction.py`. **DRIVER-WIRING FORK** (→ AskUserQuestion
  at 4c start): A_BA as its own gain slot (widen the `(L+C,S,B)` triple — ~13 unpack sites) vs composed
  into the scatter gain (`S_bulk+A_BA` in one slot — less blast; ⚠ the shared-`K_iso`-emission elegance:
  S_bulk and A_BA share the `K_iso∘integrate` emission — do NOT recompute it twice). Gates L1-L5 (L4 =
  the DECISIVE Mode-11 "gate-never-executes-the-rewired-path" driver-routing sentinel).
- **4d — Assemble CoupledOperator (WRAP-first).** Name the 2×2 exposing the four blocks; `.apply`/
  `.solve`/`.H` delegate to the current fused walk (bit-id). Stamp `system_role=COUPLED` explicitly;
  **C-fwd (elegance ruling): explicit-stamp the A_AA block = `SystemRole.A`, do NOT join-derive it**
  (transport-None poisons the join). Gates: the 4d-WRAP transient characterization (`array_equal`,
  RETIRED at 4e); **A2a** composite-G forward reciprocity (reds ALL geoms — tooth M-ADJ-metric at
  `_AdjointOperator.apply`) + **A2b** seed-isolated V_cell; **A3** swap law `(A.H).inverse()=
  (A.inverse()).H` extends `test_inverse_adjoint_coherence.py` (⚠ grep every is_adjointable/
  is_invertible False assertion when the coupled predicates land); **E4** two-anchor (φ=Q/Σ_t
  non-fissile + k_inf fissile, ≥2G).
- **4e — THE WALK UN-WEAVE (EXTRACT).** Route the walk through the built operators + delete the inline
  twins (retirement entries 1-4, current lines: entry1 `_run:4104-4131` +A_AB feed `:4135-4136`; entry2
  `_run_transpose:4633-4661`; entry3 `_apply_walk:3168-3186` +precompute `:3078-3084`; entry4
  `loss_action_transpose:3410→3498→3576→3580-3590` +`_seed_rows_transpose:2791-2846`). **HAZARD 1:
  A_AA+A_AB are FUSED in the shared kernel** (`precompute_psi_state`/`cell_contribution`/
  `angular_adjoint` compute them summed) — isolate A_AB by zeroing psi_view (as `A_AB.apply` already
  does); expect principled-equiv, NOT bit-id. **HAZARD 2: the solve is a within-(L+C) block-triangular
  forward-substitution** — ray marched FIRST (A_BB⁻¹), then `phi_aux=cells_minus` seeds the bulk loop
  (`psi_angle` is shared mutable ray→bulk state); preserve the ordering. **HAZARD 8: non-carrying
  cyl/slab share the kernel** — must stay bit-exact (the CONTROL). Retire the 4d-WRAP transient gate;
  flip A1-EXTRACT + E1-E3 live. HIGHEST-RISK.

### 4a — DONE @ `4113d27` (SystemRole {A,B,COUPLED} lattice)

NEW `SystemRole {A,B,COUPLED}` enum + `_join_system_roles` (union: same→same, different→COUPLED, None
propagates) in `numerics/operator.py`, **CLONING block_role's mechanism exactly** (base-Protocol
`system_role=None` attr + the 3 composer propagations: OperatorSum join, _AdjointOperator/ScaledOperator
passthrough). Stamped the 4 ψ½ operators: A_BB→B, B_b→B, A_AB→COUPLED, A_BA→COUPLED; model-generic
System-A ops (C/S/F/B_a) stay None (honest — no layering violation; transport-None would poison the join
anyway). `SystemRole.A` defined but UNSTAMPED in 4a (the honest 3rd of a complete partition; becomes the
CoupledOperator's A_AA-block label at 4d — **C-fwd**). Gates: `TestSystemRoleLattice`
(`test_psi_half_coupling.py`) + `TestSystemRoleEnum` (`test_operator_protocols.py` — the exact-3-member
pin mirroring BlockRole's). Verified: ratchet transport:1, 62+46 tests, both teeth mutation-RED, sphinx
-W exit 0. **elegance-enforcer PASS-WITH-NITS** (0 violations; all 3 nits fixed: SystemRole exact-member
pin, ScaledOperator propagation case, bidirectional join cross-ref + collapse trigger). **RULING:** the
parallel enum is **structurally FORCED** — B_b sits non-None on BOTH axes (block_role=BOUNDARY AND
system_role=B), which a merged enum cannot represent; the generic `RoleAxis` abstraction is premature at
2 axes (**collapse trigger = a 3rd parallel axis**: DSA/multiphysics).

### 4b — DONE @ `6b4d584` (A_BB forward via shared-kernel extraction)

**USER RULING (2026-07-08): "extract the shared kernel now", NOT a tracked twin** (chosen over the
WRAP-first default when the operator docstring flagged the twin-vs-extract fork). The forward is
single-sourced: 3 shared functions in `psi_half_angle_seed.py` beside `carlson_inward_sweep_from_source`
— `radial_characteristic_residual_march` (per-leg kernel, RELOCATED from the walk's `_seed_residual_
march`), `radial_characteristic_forward_residual` (two-leg orchestration), `..._transpose` (PURE A_BB
transpose). BOTH the `(L+C)` walk (`_seed_rows_forward`/`_seed_rows_transpose` → thin wrappers; the
transpose wrapper ADDS the A_AB `seed_cells_bar` coupling, NOT part of A_BB) AND `A_BB.apply`/
`apply_transpose` route through them. `inverse()` = the generic `InverseOperator` (apply→solve,
solve→apply); is_invertible/is_adjointable→True.

**test-architect refutations (R1-R4 — corrected my "bit-exact by construction" assumption):** solve∘apply
=id is ~FP ULP for the CELLS (the forward's 2/Δr and the march's Δr·σ+2 reassociate, L7), bit-exact 0
ONLY on the +1 outflow corner → test on the CONSISTENT subspace ψ0=solve(q0) (the +1 corner is a free
datum apply MEASURES / solve OVERWRITES). The transpose split is correct block decomposition (pure
A_BBᵀ shared; A_ABᵀ in the walk). Gate 4 is EUCLIDEAN (the metric `.H` is the composite's job, 4d).
**The predicate-flip contract was EMPTY** (A_BB not composed in production, not in
`test_capability_survival`) — the flip is additive-safe.

**Verified:** walk BIT-IDENTICAL after extraction (all blocks, fwd+transpose, array_equal); `TestA_BB_
Forward` (5 gates incl. the Mode-11 both-namespace anti-twin spy) + both teeth (fwd recurrence σ×3,
transpose sign-flip) mutation-RED w/ green controls; `tests/sn -m "not slow"` **2073 passed / 0 reds**;
ratchet transport:1; sphinx -W exit 0; **elegance-enforcer PASS** (0 violations; all 5 crux decisions
principled; 1 stale-class-docstring NIT fixed). Retirement: `_seed_residual_march` retired (relocated);
the SOLVE orchestration twin (A_BB.solve vs in-sweep `:4104-4131`) REMAINS → step 4e.

**FOLLOW-UP (optional, step 7):** a canonical positive A_BB row in `test_capability_survival.py`
(gate 5c already catches a predicate flip-back; the canonical table is a nice-to-have durable pin).

**NEXT: 4c — THE LIFT** (S/F→pure bulk, A_BA=Fold∘K_iso driver-born, fold→sn, retire the transport
reconstruction module). ⚠ HAZARD 3 (the LIFT + driver-gain-add land in ONE commit — presence law) +
HAZARD 5 (S=within-group gain, F=outer q_ext — different seams) + the **DRIVER-WIRING FORK** (→
AskUserQuestion at 4c start: A_BA as its own gain slot vs composed with the scatter gain; ⚠ the
shared-K_iso-emission elegance — S_bulk and A_BA share the `K_iso∘integrate` emission, don't recompute
twice). Re-anchor from THIS block + `git log` (trust git).

### 4c GROUNDING + OPEN QUESTIONS (2026-07-08, pre-implementation, `22a96cc`)

**DRIVER-WIRING FORK RESOLVED (user, 2026-07-08): OWN GAIN SLOT.** Gains widen
`(S,B)→(S,A_BA,B)`; A_BA a peer coupling (mirrors #208's B-from-S un-weld — `_within_group_
triple:204` "the transitional `S + B` fold is RETIRED"). S stays pure model-generic bulk.
Cost accepted: S_bulk + A_BA each call `K_iso∘integrate` — cheap, single-sourced (the SAME
`IsotropicScattering` object, not a twin). Rejected: composed `S⊕A_BA` gain (would re-weld
the very fold #208 retired).

**Own-slot driver blast (precise, `orpheus/sn/solver.py`):** 1 producer `_within_group_
triple:168` (returns `(L+C,S,B)`) + 3 helper sigs (`_within_group_si:603`, `_within_group_
krylov:372`, `_select_si_resolvent:538`) + 4 caller unpacks (`LC,S,B = _within_group_triple`
@ `:1374` eig-SI / `:1565` eig-Krylov / `:2461` fix-SI / `:2663` fix-Krylov). **A_BA is
Jacobi-ONLY** — it exists iff `radial_characteristic_space is not None` (seed-carrying 1-D
curvilinear), which ALWAYS falls to Jacobi (G-S needs `is_cartesian and not is_1d`), so the
G-S split arm NEVER sees A_BA ⟹ append A_BA to the Jacobi gains only; G-S arm untouched.

**Scatter A_BA design (clean):** `A_BA_scatter = Fold ∘ K_iso ∘ integrate_angular`, wrapping
the migrated Fold (`RadialCharacteristicReconstruction`) + the SHARED
`ScatteringOperator.isotropic_kernel` (`= IsotropicScattering(mat_xs)+IsotropicN2N(mat_xs)`,
`scattering.py:911`). House style = A_AB `RadialCharacteristicSeeding` (bespoke
`LinearOperator` subclass wrapping shared kernels). Candidate: ONE generic sn class
parametrized by the emission kernel (K_iso for scatter, the fission kernel for F) — name TBD
(user historically steers these; `RadialCharacteristicEmission`?). Fold MATH kernel
(`fold_moments_to_radial_characteristic{,_transpose}`, `numerics/spaces/radial_characteristic_
space.py:447/510`) STAYS (data layer, deliverable 5).

**⚠ OPEN — F ray-seed is a TWO-FOLD-FACTORY question (HAZARD 5, subtler than "compute_fission_
source"):**
- The eigenvalue F ray seed does NOT ride F-fwd (`fission.py:522`, `RadialCharacteristic
  Reconstruction`): `compute_fission_source:951` passes a bare scalar ndarray → fission's
  scalar arm (no seed). It rides `_radial_characteristic_source_from_per_ordinate`
  (`solver.py:1364`, eig-SI q_ext assembly) = `RadialCharacteristicSourceSink.from_angular_
  source(q_ext_per_ord)` (`:526-535`).
- So there are TWO ψ½-fold factories: (a) `RadialCharacteristicReconstruction` (moments-fold,
  step-2's named A_BA fold, used by S/F `.apply` arms) and (b) `from_angular_source`
  (per-ordinate-fold; docstring `:529` claims "the ONE q½-fold factory" — solver +
  fixed-source rhs + `transport_sweep`). They AGREE at ℓ=0 (production reach):
  `from_angular_source ≈ Fold∘integrate`, `Reconstruction = Fold`.
- **QUESTION (test-architect map + user):** does 4c UNIFY these (A_BA = the single ray-seed
  path, retiring `from_angular_source`'s driver route) or leave `from_angular_source` for its
  OTHER consumers (fixed-source rhs, `transport_sweep`)? The invariant "A_BA ≡ today's ray
  seed bit-identical" must pin against whichever is TODAY's LIVE path (eig = `from_angular_
  source`, NOT `Reconstruction`). test-architect dispatched (bg, own-slot gate spec L1–L5 +
  S-adj asymmetry ruling + refutations).

### 4c GATE SPEC + OBLIGATIONS (test-architect memo, own-slot, 2026-07-08)

Full table: `.claude/agent-memory/test-architect/coupled_operator_step4_verification.md`. Durable
load-bearing findings:

**THE reframe — NO adjoint SI driver in production** (`grep 'apply_transpose|\.H' solver.py` =
EMPTY). Forward own-slot gains call only `A_BA.apply`. The S-adj `K_isoᵀ∘Reconstructionᵀ` bulk
pullback is exercised ONLY by test-side `.H` reciprocity gates, which feed a **present-ZERO seed
cotangent** ⟹ a lost pullback is a **FALSE GREEN** in every existing adjoint gate (R2).

**THREE obligations on the ONE commit (HAZARD 3):**
1. LIFT (S/F ray-arm drop) + own-slot gain-add land TOGETHER — the intermediate mixed-presence
   state trips `evaluate_residual:300` / `FullField.__sub__`; acceptance = full `tests/sn -m "not
   slow"` green in the single commit.
2. **S-adj reciprocity-leaf:** A_BA MUST realize `apply_transpose` at 4c (NOT deferred to 4d), AND
   the reciprocity gate's loss `(L+C−S−B).H` MUST gain the A_BA leaf → `(L+C−S−A_BA−B).H` in the
   SAME commit — else `.H` silently loses the pullback (hidden by present-zero seeds). **Sole
   catcher = L1-ADJ with a NONZERO seed cotangent** (A2a/A2b are present-zero-blind).
3. **Pinned-row FLIP (R3):** `TestRegressionFloor::test_bulk_to_ray_coupling_lives_in_the_lagged_
   scattering_gain` (test:206) asserts `S_sb>1e-6` on the MONOLITHIC `S.apply`; the LIFT makes S
   pure-bulk ⟹ `S_sb→0` ⟹ REDDENS. Migrate it in the LIFT commit (re-point the claim to the A_BA
   gain). Do NOT let it fail as "expected".

**S-adj asymmetry ruling:** the pullback `bulk_bar += w·K_isoᵀ(Reconstructionᵀ χ_seed)`
(`scattering.py:1837-1844`) = `A_BA.apply_transpose(χ_seed)` (Foldᵀ→K_isoᵀ→integrateᵀ=×w). After
the lift: `S.apply_transpose` pure-bulk; `A_BA.apply_transpose` carries the pullback; sum = monolith.

**F outer-seam (HAZARD 5) + the two-factory resolution:** the eigenvalue F ray seed is TODAY built
at `solver.py:1364` via `from_angular_source` (NOT F-fwd's `Reconstruction` arm — `compute_fission_
source:951` passes a bare scalar ndarray → fission's SCALAR arm, no seed). 4c migrates the `:1364`
route → `A_BA_fission.apply` (L4-F sentinel). `from_angular_source` STAYS for its non-fission
consumers (fixed-source rhs, `transport_sweep`). ⚠ the migration REMOVES a `from_isotropic∘integrate`
round-trip ⟹ F ray seed may be principled-equiv (ULP), NOT bit-id — MEASURE (`array_equal` if exact,
else `assert_array_almost_equal_nulp`). ⚠ CHECK F-fwd's ray arm (`fission.py:508-534`) production
reachability in the retirement audit (may be DEAD ⟹ pure retirement).

**Gates** (all `-O`): L1-FWD (S/F pure-bulk fwd, `array_equal`) · **L1-ADJ** (nonzero-χ_seed pullback
catcher — THE decisive one) · L2 (A_BA=Fold∘K_iso∘integrate bit-id vs `_ba_oldloop_reference` +
Mode-11 fold-spy `n==2·n_levels`) · L3 (`S_bulk⊕A_BA≡monolith`, `array_equal` the OBJECT — Mode-12) ·
**L4-S** (own-slot driver-routing sentinel: wrap `A_BA_scatter.apply`, `n>0` in a real solve +
`A_BA in gains`) · **L4-F** (F outer-seam: wrap `A_BA_fission.apply`, `n>0` in outer iter +
`A_BA_fission NOT in within-group gains`) · L5 (retirement + `_rcr_mod` repoint). Driver mechanics:
SI `rhs+=g.apply` (`iteration.py:619`), Krylov `out-=g.apply` (`:868`); widen `_within_group_triple:
249`→`(L+C,S,A_BA,B)`, `_select_si_resolvent:599`, `*gains` spread `:658`. R1: `array_equal` legit
(fold home-move, no re-assoc) EXCEPT the F round-trip + any summed-rhs reduction-order (→ nulp).

### 4c EXECUTION FINDINGS (main agent, 2026-07-08) — scatter DONE, F reframed

**Scatter LIFT built + verified (production edits landed, uncommitted):**
- `RadialCharacteristicEmission` (A_BA, generic over `emission_kernel`) + migrated
  `RadialCharacteristicReconstruction` (Fold) in `sn/operators/radial_characteristic.py`.
  SMOKE-VERIFIED: scatter `apply` bit-id vs the old fold (maxabs **0.0**); `apply_transpose`
  an EXACT Euclidean adjoint (defect **0.0**); `F.kernel ≡ F.apply` emission bit-for-bit.
- S-fwd/S-adj/F-fwd → PURE BULK (present-zero ray). Scatter A_BA wired as its OWN gain via the
  new single-source `_lagged_gains(S, B, sn_mesh)` (SI Jacobi arm + both Krylov calls) — this
  AVOIDS the `_within_group_triple` rename (99-ref blast) the plan's "widen the triple" implied;
  A_BA injected at the gain-assembly seam instead. **The scatter LIFT is BIT-IDENTICAL end-to-end**
  (the scattering seed just moved S.apply → A_BA_scatter; q_ext's fission seed untouched).
- `tests/sn/operators/test_psi_half_coupling.py`: **6 reds, all expected** — the R3 flip
  (`TestRegressionFloor::test_bulk_to_ray_coupling_lives_in_the_lagged_scattering_gain`) + 5
  step-2 `TestA_BA_SchurFold` gates that asserted S/F EMIT the ray (now A_BA's job). These are
  the migrate-to-step-4-architecture list. 45 pass (bit-id elsewhere).

**⚠ F REFRAME — the two-factory premise was FALSE (no twin):** `from_angular_source`
(`radial_characteristic_source_sink.py:157`) computes Legendre moments (`legvander`+einsum) then
folds through the SHARED `fold_moments_to_radial_characteristic` kernel — the SAME kernel
`Reconstruction`/A_BA use. So the fission ray seed ALREADY routes through the one fold (no
Cardinal-Rule-2 twin); `from_angular_source` is the PER-ORDINATE typed entry, `Reconstruction`
the MOMENTS entry, over one kernel. Further: **F-fwd's ray arm is DEAD in production**
(`fission_op.apply` only ever gets a scalar — the drop is pure retirement), and the live fission
seed rides the OUTER `q_ext` (`solver.py:1401/1570`), which is a valid FACTORED A_BA_fission
(F.kernel in `compute_fission_source` + `from_angular_source`'s fold). So migrating it to
`Reconstruction(fission_source[None])` is NOT a twin-fix — it is a ULP re-baseline + shape-fiddling
(fission_source is `(ng,nx,ny)`; the fold wants `(n_moments,ng,nx)`) that merely removes a
`from_isotropic→re-project` round-trip. Site `:1956` folds the TOTAL source (fission+scatter), NOT
pure fission — MUST keep `from_angular_source`.

**⟹ RECONSULTING the user** (the "migrate fission only" ruling rested on the false twin premise):
recommend **scatter LIFT bit-identical THIS commit; leave the fission seed on `from_angular_source`**
(already folds through the shared kernel — no twin; more general = handles anisotropic); file the
round-trip removal as a minor-optimization follow-up if wanted. RadialCharacteristicEmission stays
generic (machinery; F.kernel smoke-verified) even though production instantiates it for scatter.

**USER RULING (2026-07-08): TWO immediate back-to-back commits, NOT a follow-up issue.**
- **Commit 1 (4c) = scatter LIFT — DONE @ `fbcb5aa` (2026-07-08).** S/F ops→pure bulk;
  `RadialCharacteristicEmission` (A_BA) + migrated Fold in sn; A_BA_scatter own gain
  (`_lagged_gains`); the `_EmissionKernel` Protocol (typing the apply+apply_transpose contract
  the bare `LinearOperator` doesn't surface); reciprocity leaf `(L+C−S−A_BA−B).H`; gates
  `TestCoupledLift` {L1-FWD/L1-ADJ/L2/L3/L4-S} + R3 row re-point + the 3-gain single-primitive
  contract; retired the transport `radial_characteristic_reconstruction.py`. Fission seed UNCHANGED
  (from_angular_source). BIT-IDENTICAL (0-ULP smoke + curvilinear keff). VERIFIED: tests/sn -m "not
  slow" 2078/0, ratchet transport:1, sphinx -W 0, elegance PASS (all 3 conditions met). test-arch
  refutation banked: deleting the WHOLE A_BA leaf keeps reciprocity green (leaf-invariant, Mode-12)
  ⟹ L1-ADJ (nonzero seed) catches a WRONG transpose, L4-S sentinel catches a DROPPED gain.
- **Commit 2 (immediately after, this session) = F outer-seam migration — PRODUCTION DONE
  (uncommitted).** NEW `_radial_characteristic_fission_seed(fission_source, sn_mesh)` helper
  (`solver.py`) = the direct moments-fold `Reconstruction(fission_source[None])`; the eig-SI
  (`:1401`) + eig-Krylov (`:1570`) fission seed routes through it (`from_isotropic→from_angular_
  source` round-trip GONE). `:1956` (TOTAL source) + `:2172` (external) KEEP `from_angular_source`
  (the per-ordinate typed entry — docstrings now clarify the two typed entries over the ONE kernel,
  no twin). VERIFIED principled-equiv (maxabs **1.2e-15** vs the old route, NOT bit-id — the removed
  `·w` round-trip reassociates); curvilinear/homogeneous eigenvalue suite ABSORBS it (31 passed,
  unchanged); ratchet transport:1 (sn 0). **DONE @ `88d94b4`.** Gates (test-architect): L1F fission
  value gate (bit-id vs moments-fold + nulp=32 vs the retired route) + L4F outer-seam sentinel +
  teeth. **test-architect refutation banked (memo L4-F row):** a GLOBAL `Reconstruction.apply`
  counter is Mode-11-BLIND here — the scatter A_BA folds through the SAME reader every SI iteration
  (`global_fold=322` even with the fission seed reverted), so the sentinel must measure the
  SEAM-SPECIFIC fold DELTA (`seam_fold_delta=0` on revert → RED). VERIFIED: tests/sn -m "not slow"
  **2082 passed / 0 reds**, ratchet transport:1, sphinx -W 0.

**⟹ STEP 4c (THE LIFT) COMPLETE @ `fbcb5aa` (scatter, bit-id) + `88d94b4` (fission, ~ULP).** S/F are
pure bulk; A_BA = `RadialCharacteristicEmission` (scatter = own within-group gain; fission = the
direct moments-fold at the outer q_ext seam); the Fold migrated to sn; the transport reconstruction
module retired; the reciprocity `.H` carries the A_BA leaf. NEXT = **Step 4d** (assemble the 2×2
`CoupledOperator`, WRAP-first) — note the 4d typing reconciliation elegance-enforcer flagged: A_BB/
A_AB are sub-space block-typed, A_BA/B_b are FullField-embedded; 4d chooses how to reconcile (wrap
A_BB/A_AB to embedded, unwrap A_BA/B_b to blocks, or embed-at-assembly). Owed docs (step 4/5): the
A_BA/A_AB block-operator theory-page prose.
Rationale (user): follow-up work we CAN do immediately should be done immediately, not filed —
but bit-id and ULP-re-baseline stay in SEPARATE commits (bisectable, each cleanly verifiable).

---

## Step 4 RE-SCOPED (2026-07-10) — the N-general block machinery (user design dialogue)

**The 4d "WRAP-first name the 2×2" framing is SUPERSEDED.** A design dialogue (2026-07-09/10)
re-posed 4d around the operational criterion "*if we assembled the explicit matrix with all
blocks in it, it must work*" and then generalized. Settled design + rulings below; the earlier
4d gate memo (`coupled_operator_step4_verification.md`) is superseded on structure (its
invariants — type-safe matvec, `.H` Mode-12 reciprocity, two-anchor value — carry forward).

### Why the OperatorSum/present-zero route was REJECTED (user, decisive)

The flat `OperatorSum` + `FullField→FullField` present-zero padding **keeps wrong
multiplications representable** — a padded block accepts any `FullField`; nothing at the type
level stops a ray operator from receiving a bulk field or an emission from landing in the trace
slot. Present-zero padding *is* the loss of Pattern-4 safety; `system_role` tags are runtime
metadata, not type constraints. The honest object is a **typed block operator over a typed
block vector**, where the block matvec is the only spelling and type-checks per block.

### The settled design — semantics-agnostic, N-general (the machinery, not the ψ½ special case)

- **`CoupledField` = N systems**, each a complete `(interior ⊕ boundary)` composite — i.e.
  each system is a `FullField` in its own right.
- **`CoupledOperator` = the N×N block grid** `A_ij : System_j → System_i` (diagonal = each
  system's self-operator `L+C−S−B`; off-diagonals = couplings; missing block = zero map, so
  coupling sparsity is explicit).
- **Block matvec** `y_i = Σ_j A_ij x_j`, per-block TYPE-CHECKED (`A_AB @ x_A` is a type error —
  illegal states unrepresentable, Pattern 1 ∘ Pattern 4).
- **Three consumption modes at the block level** (the stencil-assembly `Op → Mat` functor):
  `assemble` places each block at its `(row_i, col_j)` OFFSET → one flat Mat₄ → direct factor
  (small); `apply` = block matvec (Krylov); `solve` = block-triangular / block-Jacobi /
  block-G-S (large). **The block structure PROVIDES the offsets → this is the scoped
  realization of the deferred `LocalToGlobalMap`** (the "structured consumer" the 2-P0
  ruling waited for; `assembled_operator.py:48`). Block granularity = a solve-strategy + size
  knob, not physics.
- **The co-producing mechanism** `build_coupled_system(sn_mesh, mat_xs) → (CoupledOperator,
  CoupledSpace)` emits the operator AND the matching block-field space together (aligned by
  construction); keyed on the mesh's ray-carrying status → non-carrying degenerates to 1×1
  (no System B) → applying a System-B block is unconstructable. **This SUBSUMES step 6**
  (structural presence: `Optional[ray]` + the 7 guards dissolve — presence = block existence).
- **The ray LEAVES the augmented `FullField`** and becomes System B's own `(interior ⊕
  boundary)` composite. Folds step 6's ray-extraction into this build.
- **Within System A stays flat-sum + `block_role` tags** — `L` is irreducibly FULL (couples
  bulk↔trace), so bulk/trace is NOT cleanly block-separable; typed-2×2 belongs only where the
  decomposition is clean (the A/B split). `SystemAField` remains a composite. Right tool per
  level (nested: typed blocks at A/B, flat-sum within A).

### The ψ½ instance (#1) — a 2×2 special case of the above

System A = the entire SN transport (`FullField` bulk⊕boundary), `A_AA = L+C−S−B_a`; System B =
the radial-characteristic closure (its own interior⊕boundary), `A_BB = RayOp − B_b`; off-diag
`A_AB` = ray→bulk seed (`RadialCharacteristicSeeding`), `A_BA` = bulk→ray emission
(`RadialCharacteristicEmission`, RE-TYPED from FullField-embedded to the true `SystemA→SystemB`
block — the elegance memo's anticipated "re-type at 4d").

### FOUNDATIONAL insight → GH #296 (the FORWARD quest, deferred)

SN is ALREADY a coupled block system, solved with implicit per-axis block strategies (space =
sweep = block-triangular direct #284; angle = source iteration = block-Jacobi over scattering;
energy = multigroup = block-G-S #273). Monolithic treatment was a CHOICE. Reifying "block"
unifies solve-strategy posing across SN axes, other transport methods, domain decomposition
(sparse inter-domain coupling), and CP-finer-than-interface-current (volumetric coupling).
**#296 tracks reframing the SN stack (and beyond) onto blocks — a FUTURE quest.** THIS campaign
builds the machinery N-general in SHAPE but migrates ONLY ψ½ (existing axes keep their
correctness-gated representations). User: "leaving any SN stack reformulation to a future quest."

### RULED sub-design (2026-07-10, USER: path (i))

**RULED path (i)** (user: "we're building already, so let's make it right"): generalize
`FullField` → carrier-generic `System[Interior,
Boundary]` (System A = `System[AngularFlux, AngularBoundaryFlux]` bit-id; System B =
`System[RadialCharacteristicFlux, RadialCharacteristicBoundary]`) — the honest uniform object,
but touches the widely-used `FullField`; vs (ii) keep `FullField` for System A, give System B a
bespoke `(interior⊕boundary)` composite — lighter, but "each system is a FullField" stays
conceptual not typed. LEAN = (i) (user framing "System B would also have a FullField").

### RULED: B1 rename + the corrected A/B/C sequencing (2026-07-10)

**B1 RULED (user): rename the carrier INTERIOR ACCESSOR `.bulk` → `.interior` NOW** ("this
rename will never be as easy as today"; sub-agents execute, main agent owns scope + verifies).
SCOPE = the **accessor family ONLY**: `.bulk`→`.interior` (564), `bulk=`→`interior=` (280),
`bulk_space`→`interior_space` (36), `n_bulk`/`bulk_shape` internals, `CompositeField.bulk`
protocol, docs `.bulk`/`:attr:` refs. **CONSCIOUS KEEP** (distinct axes — NOT off-pattern):
`BulkField` (the interior-field TYPE — `.interior` *holds* a `BulkField`),
`BlockRole.BULK`/`BulkOperator` (the operator ACTION-FOOTPRINT classification), compound names
(`BulkAnalysisOperator`). Each family stays internally consistent ⟹ greppability holds; a
full-lexicon `bulk`→`interior` is an OPTIONAL later pass (additive, no rework).

**ORDERING CORRECTION (explorer, load-bearing):** `System[Interior, Boundary]` is 2-block;
today's `FullField` is 3-block (`interior ⊕ boundary ⊕ Optional[ray]`). You CANNOT make
FullField a PURE 2-block `System` while the ray lives in it → the ray must LEAVE first; pure-2-
block-ness FALLS OUT of the eviction, it is NOT a separate first step. (Supersedes the earlier
4d.0-first framing.) FullField blast-radius (explorer): 564 `.bulk` readers, 280 `bulk=`,
382-degree hub; the numerics drivers (`iteration.py` gain-fold) + `OperatorSum`/`_AdjointOperator`
are the ONLY composite-generic seams (stay bit-id); everything else is shape-specific.

**Corrected phases:**

- **Phase A (bit-identical) — rename + carrier-generic.**
  - **A1 (mechanical, sub-agents): DONE @ `42f1a34`** (93 files; tests/sn **2082/0 = baseline**,
    diffusion+homogeneous+cp 276, ratchet transport:1; verified structurally clean — the streamed
    `bulk`/`starting_direction` diagnostics were STALE LSP #226 noise, grep+ratchet confirmed).
    KEPT (conscious): `BulkField` type, `BlockRole.BULK`/`BulkOperator`, compound names, prose
    domain-noun "bulk". **DOCS-DEBT:** the docs `:attr:`bulk``/`psi.bulk` accessor refs are stale
    (silent-render, no `-W` break) — DEFERRED to ONE consolidated archivist pass AFTER A2 (covers
    both the accessor rename AND the new `System[Interior,Boundary]` carrier docs).
  - **A2 (surgical, main agent):** generalize `FullField` → carrier-generic `System[Interior,
    Boundary]`, System A = `System[AngularFlux, AngularBoundaryFlux]`, ray KEPT Optional for now;
    diffusion/CP (already pure 2-block `FullField(interior, boundary)`, zero ray) ADOPT the
    generic base — the living bit-id proof. Gate: per-consumer bit-id (ERR-063 discipline — direct
    value compare, not a proxy) + the multi-instantiation crux `System[ScalarFlux,
    ScalarBoundaryFlux]` (a hardcoded-`AngularFlux` bug in the generic body reds scalar / greens
    angular — L20 sharpened).
- **Phase B (principled-equiv 5.5e-16) — System B + machinery + ray eviction.** Build System B =
  `System[ray-interior, ray-boundary]` + `CoupledField`/`CoupledOperator` (semantics-agnostic,
  N-general: type-checked block matvec, block `assemble` = offset placement → flat Mat₄, block
  `.H`); RE-TYPE `A_BA`/`B_b` to true `SystemA↔SystemB` cross-space blocks; block-aware driver
  `[ψ_A, ψ_B]`; **evict the ray from FullField** → System A falls out pure-2-block; retire the
  mixed-presence law (SUBSUMES step 6). ⚠ **ERR-053**: re-size GMRES `restart = n_dof` =
  Σ both systems' `to_flat().size` or GMRES silently truncates. Gates: block matvec type-safety;
  **assemble ≡ probe** (principled-equiv `rtol=1e-11` via `_dense(.apply)` NOT `.as_matrix()` —
  a tautology trap; offset-swap tooth on an ASYMMETRIC toy 2×2 — symmetric is offset-transpose-
  blind); block `.H` Mode-12 reciprocity (reds ALL geoms, "slab stays green" is FALSE).
- **Phase C — block solve (step 5) + walk un-weave (4e).**

Gate spec (re-written for the re-scope):
`.claude/agent-memory/test-architect/coupled_operator_step4_verification.md`.

### A2 design — extract the generic `System[Interior, Boundary]` base (2026-07-10)

**Key finding (grounding `full_field_space.py`):** `FullFieldSpace` + the `CompositeField`
protocol (numerics) are ALREADY carrier-generic — family-blind, duck-typed on
`.interior`/`.boundary`/(optional `.radial_characteristic`), used by SN (angular trace) AND
diffusion/CP (scalar trace). So the SPACE side is DONE; A2 extracts the generic CONCRETE composite
BASE from the `FullField` DATACLASS (whose 2-block ALGEBRA — `_recombine`/dunders/`to_flat`/
`from_flat`/mesh-identity — is already carrier-agnostic in structure; only the leaf TYPES
`interior: BulkField`, `boundary: BoundaryField` are specialized).

**Approach (A2-i):** a generic 2-block frozen-dataclass base `System[Interior, Boundary]`
(`kw_only`, holds the 2-block algebra) + `FullField(System[AngularFlux/BulkField,
AngularBoundaryFlux/BoundaryField])` subclass that ADDS `Optional[radial_characteristic]` and
overrides the hooks to thread it — the ray is transport-SN-specific, NOT part of the generic
2-block System, a TEMPORARY subclass extension Phase B removes (collapsing FullField → pure
`System[·,·]`). diffusion/CP adopt the PURE base `System[ScalarFlux, ScalarBoundaryFlux]` (the
living bit-id proof + the multi-instantiation gate — a hardcoded-`AngularFlux` bug in the generic
body reds scalar). `TimedFullField(FullField)` unchanged. `FullFieldSpace` already `getattr`s the
ray → works for both the 2-block base and 3-block FullField.

**Ray-threading:** base `_recombine(*, interior, boundary)` is 2-block; FullField overrides to add
`radial_characteristic`; `FullFieldSpace._rebuild` threads the ray only when the field carries it
(already field-driven via `getattr`/`_seed_space_for`). **NAME RULED: `Composite[Interior, Boundary]`**
(user, 2026-07-10 — structural not domain-role; "System" misleads under domain decomposition; the
principle is saved to [[feedback-high-signal-names]]).

### A2.1 EXECUTION — Composite extraction (2026-07-10, in verification)

**Realized:** `full_field.py` now holds `Composite[Interior, Boundary]` (generic 2-block base,
`Interior`/`Boundary` bound to `BulkField`/`BoundaryField`) + `FullField(Composite[BulkField,
BoundaryField])` subclass. The algebra lives ONCE on `Composite`: the six dunders + `copy` route
through **two per-shape hooks** `_map_binary` / `_map_unary` (+ `_recombine`), and the flat protocol
through `_flat_parts` / `_from_flat` (INSTANCE hooks, so `to_flat`/`from_flat` live once with NO
classmethod-override narrowing — the from_flat Liskov fix). FullField overrides ONLY those hooks to
thread the ψ½ block. `TimedFullField` UNCHANGED (inherits FullField's ray-aware hooks; overrides
`_recombine` alone). Error messages preserved VERBATIM (the shared `"same-class partner"` Field-base
vocabulary kept — a divergence to `"Composite partner"` was reverted).

**Gate (bit-identical):** ratchet `transport:1`; composite carrier `test_full_field`+`test_timed_full_field`
**51 passed**; `transport`+`numerics`+`diffusion` **1355 passed**; `tests/sn` **2082/0 = baseline**. The
streamed `full_field.py` import diagnostics are #226 LSP noise (ratchet CLI clean).

## ⏸ PHASE A COMPLETE (2026-07-10) — recommended /compact point

Three commits, all verified, branch `refactor/sn-walk-unification` **69 ahead of main, pushes HELD**:
- **A1 `42f1a34`** — `.bulk`→`.interior` accessor rename (bit-identical).
- **A2.1 `c439bf4`** — extract the generic `Composite[Interior, Boundary]` base; `FullField(Composite
  [BulkField, BoundaryField])` subclass with the temporary ψ½ block; algebra ONCE via the
  `_map_binary`/`_map_unary`/`_recombine` + `_flat_parts`/`_from_flat` hooks; `TimedFullField` unchanged
  (bit-identical, sn 2082/0).
- **A2.2 `d6b1490`** — `test_composite.py`: the pure `Composite` base's intrinsic affine-torsor laws on a
  SCALAR carrier (the base hooks are dead-until-tested since FullField overrides them; doubles as the N2
  multi-instantiation gate). 14 passed.

**NEXT = Phase B** (principled-equiv 5.5e-16) — the HARD arc. Build System B = `Composite[ray-interior,
ray-boundary]` + `CoupledField`/`CoupledOperator` (semantics-agnostic N-general: type-checked block
matvec, block `assemble` = offset placement → Mat₄, block `.H`); RE-TYPE `A_BA`/`B_b` to true
`SystemA↔SystemB` cross-space blocks; block-aware driver `[ψ_A, ψ_B]`; **evict the ray from FullField**
→ FullField collapses to pure `Composite`; **relax the `CompositeField` protocol (numerics) to 2-block**
(drop the required `radial_characteristic` member — this is what unblocks diffusion adopting the pure
`Composite` base); retire the mixed-presence law (SUBSUMES step 6). ⚠ **ERR-053**: re-size GMRES
`restart=n_dof` = Σ both systems' `to_flat().size`. Gates (test-architect memo
`coupled_operator_step4_verification.md`): block matvec type-safety; **assemble ≡ probe** (principled-
equiv rtol=1e-11 via `_dense(.apply)`, offset-swap tooth on an ASYMMETRIC toy 2×2); block `.H` Mode-12
reciprocity (reds ALL geoms). Re-anchor from THIS block + `git log` (trust git).

---

## ⏸ PHASE B.1 COMPLETE (2026-07-10) — System B IS an independent composite; recommended /compact point

**"Pose System B as an independent composite" (the user's explicit Phase-B entry) is DONE** in four
verified, ADDITIVE commits on `refactor/sn-walk-unification` (**74 ahead of main, pushes HELD**). The
unified `RadialCharacteristicField` path + all its consumers are UNTOUCHED — the split coexists until the
operator-retyping arc migrates + retires it. Ratchet stayed `transport:1` every commit.

- **B.1a `030e550`** — relax the `Composite` base bounds `BulkField`/`BoundaryField` → `Field`; the
  System-A locus guards RELOCATE to `FullField.__post_init__` (verbatim messages, BEFORE super), the
  `mesh` property reads via `getattr` (generic `Field` carries no mesh). USER-RULED: relax-to-Field +
  subclasses-guard (over locus-marker ABCs). **BIT-IDENTICAL** (full not-slow wall 3592/0).
- **B.1b `df05587`** — the two split spaces `RadialCharacteristicInteriorSpace` (cells, `G_sd=V_cell`) /
  `RadialCharacteristicBoundarySpace` (corner, `G=V(R)`) on a shared `_RadialCharacteristicSubSpace`
  base, + `SNMesh.radial_characteristic_{interior,boundary}_space` cached props. CLEAN-BEFORE-EXTEND: a
  NEW `_radial_characteristic_legs` (`_RayLeg` NamedTuple) + `_validate_ray_space_inputs` single-source
  the leg order + metric; the unified `for_levels` + mesh prop REFACTORED through them (bit-id). Gate
  (20): arange-fingerprint split-fidelity (NOT split==unified — self-referential to the shared walk),
  the MANDATORY multi-level `for_levels((0,2,5))` fixture (mutation-verified: reverse-levels is
  sphere-GREEN / multi-RED — the config-blindness proof), metric partition + conservation, presence.
- **B.1c `59b2ad0`** — the split leaf types: two INDEPENDENT `FaceField[tuple[int,int]]` loci
  `RadialCharacteristic{Interior,Boundary}Field` (like Bulk vs Boundary, NOT a shared base) + their
  `…Flux` + `…Displacement` role leaves, registered (displacements EAGERLY for `_BY_REP`). The two flux
  carry DISTINCT `_carrier_rep`s (their distinct Field bases) ⟹ `flux ⊖ flux` mints the CORRECT per-locus
  displacement (no collision — the test-architect's flagged hazard, closed by construction). Gate (19).
- **B.1d `fa61797`** — `RadialCharacteristicComposite(Composite[…InteriorFlux, …BoundaryFlux])` — a
  TRIVIAL subclass (NO hook overrides — the payoff of B.1a; System B adds no 3rd block) + `from_mesh`
  (presence-gated) + the split-fidelity BRIDGE `from_unified`/`to_unified` (`to_unified∘from_unified==id`
  bitwise = the retirement licence). Gate (16) incl. the multi-instantiation crux — mutation-verified:
  re-tightening the base bound REDs the ray composite while System-A FullField stays GREEN (System B IS
  the catcher for a B.1a regression).

**Design rulings locked this phase (durable):** (1) relax-to-Field + each specialization guards its
concrete leaf types on the FIELD base (admits both flux + displacement composites). (2) interior/boundary
are DIFFERENT LOCI (independent field bases) but the SPACES share a base (same structural kind). (3)
single-source the `(level,sign)` leg walk across all three ray spaces (Pattern 2 — the split can't drift;
the fidelity gate is a safety net). (4) the split-fidelity oracle is a HAND-KNOWN arange fingerprint +
the mandatory multi-level fixture (single-level sphere is offset-blind — L20). (5) `RadialCharacteristic
Composite` name chosen by the main agent (flagged for the user; could converge to the freed
`RadialCharacteristicField` at unified-leaf retirement). Specialist memo: test-architect
`coupled_operator_b1_split_verification.md` (G-A/G-B/G-C + the config-blindness ledger).

**NEXT — the HARD principled-equiv arc (the operator retyping, still Phase B):** RE-TYPE `A_BA`
(`RadialCharacteristicEmission`) + `B_b` (`RadialCharacteristicBoundaryOperator`) from FullField-embedded
→ true `SystemA↔SystemB` cross-space blocks acting on `RadialCharacteristicComposite`; build the typed
`CoupledField`/`CoupledOperator` machinery (test-architect `coupled_operator_step4_verification.md` M1-M5
+ A2 + the assemble≡probe centrepiece); block-aware driver `[ψ_A, ψ_B]`; **evict the ray from FullField**
(→ it collapses to a pure 2-block `Composite`; migrate diffusion/CP onto the pure base — the living
bit-id proof); relax the numerics `CompositeField` protocol to 2-block (drop the required
`radial_characteristic`); retire the mixed-presence law + the unified `RadialCharacteristicField` (the
`from_unified`/`to_unified` bridge licenses it). ⚠ ERR-053 GMRES `restart=n_dof` re-size. Re-anchor from
THIS block + `git log` (trust git).

---

## Phase B.2 — RULED decomposition (2026-07-10, user; session resumed on Fable 5)

Two AskUserQuestion rulings at arc start:

1. **Additive spine, 4 sub-steps** (over one atomic landing):
   - **B.2a** — the numerics `CoupledSpace`/`CoupledField`/`CoupledOperator` machinery
     (semantics-agnostic, N-general; gates M1–M5 on synthetic asymmetric toys per the
     test-architect step-4 memo). ADDITIVE, zero production touch.
   - **B.2b** — re-type `A_BA` (codomain) + `B_b` (domain/codomain) onto
     `RadialCharacteristicComposite` (P3 value-unmoved gates; re-point `TestA_BA_SchurFold`
     / `TestB_b_RayBoundary` probes), with **thin transient FullField-gain adapters** at the
     gain seam (`_lagged_gains` / the `B_a + B_b` sum) so production stays bit-identical.
     Mints the split-locus SOURCE-SINK role leaves the block codomains need.
   - **B.2c** — `build_coupled_system` co-producing builder (P1–P2 presence-structural;
     SUBSUMES step 6 / task #34).
   - **B.2d** — driver iterate → `CoupledField [ψ_A, ψ_B]` (W1 Mode-11 sentinel) + **ray
     eviction from FullField** (pure 2-block collapse; diffusion/CP onto the pure base) +
     `CompositeField` protocol → 2-block + mixed-presence-law retirement + the B.2b adapters
     retired + E4 anchors + ⚠ ERR-053 `restart=n_dof` re-size.
2. **The unified leaf DEMOTES in B.2, it does not retire** — block-level carriers
   (CoupledField, block signatures, driver, A_BA/B_b) go split-composite; the A_BB/A_AB
   ENGINES + the fused `(L+C)` walk keep the unified layout internally behind
   `from_unified`/`to_unified` bridges at their block boundaries (transient, cheap). Full
   unified-leaf retirement lands at **Phase C (4e)** with the walk-slot rewrite. (CORRECTS
   the B.1 checkpoint's "retire the unified `RadialCharacteristicField` in B.2" — that is
   unreachable before 4e: the walk's `cells_view`/`corner_view` internals march the unified
   buffer. The bridge is the licence; 4e is the executioner.)

Grounding dispatched at arc start: explorer blast-radius memo
`.claude/agent-memory/explorer/b2_ray_eviction_blast_radius.md` (the
`FullField.radial_characteristic` consumer map, mixed-presence/present-zero sites,
`CompositeField` protocol consumers, unified role-leaf inventory, the walk's narrowest
bridging boundary, ERR-053 sites) — feeds B.2b/B.2d scoping. **Its four structural
verdicts:** V1 — the (d) eviction is ATOMIC across carrier + space + protocol
(`_rebuild`/`_recombine` kwarg coupling forces one commit); V2 — the fused-walk bridge is
exactly SIX signatures (`InvertibleOperator.solve`/`solve_transpose` streaming.py:1049/1145,
`loss_action`/`_transpose` loss_rep:2860/3319 + packs) and everything below marches
`values + unified slot views` — the demote ruling is exactly right; V3 — machinery must
precede the re-type (the gain sum + OperatorSum guard need the typed carrier to exist);
V4 — A_BB is ALREADY leaf-typed (its re-type is a leaf swap, no padding to shed). B.2 must
mint the split **SourceSink** pair (A_BA codomain + B_b emission); the split **Residual**
pair only iff `evaluate_residual` re-types at (d). ERR-053 sites auto-track iff
`CoupledField` satisfies the ravellable protocol (it does — pinned by the interop gate).

### B.2a — DONE @ `81cfd25` (2026-07-10)

**The N-system machinery landed in numerics** (`orpheus/numerics/coupled_system.py`):
`SystemField` (member contract = the iteration.py ravellable pair + copy; arithmetic
deliberately duck-typed — the affine-torsor member signatures fit no simple Protocol
spelling), `CoupledField` (member-wise algebra delegation; the flat protocol packs in
system order ⟹ the Krylov boundary + every `restart = n_dof` site count both systems
automatically — ERR-053 closed by conformance), `CoupledSpace(FunctionSpace)`
(member-wise metric dispatch; `system_slices` = **the scoped LocalToGlobalMap** the
assembly layer's 2-P0 ruling deferred), `CoupledOperator(LinearOperator)` (typed N×N
grid; mis-placed space-declaring block UNCONSTRUCTABLE; None block = structural zero;
`assemble()` scatters block emissions via `scipy.sparse.block_array` →
`SparseAssembledOperator`; NO solve/inverse — step 5). **The Hilbert adjoint comes free
and Mode-12-closed**: zero adjoint code here — `.H` = the existing `_AdjointOperator`
over the Euclidean transposed grid + the member-wise metrics (measured reciprocity
1.3e-16; the M-ADJ-metric tooth reds O(1)).

Gates: `tests/numerics/test_coupled_operator.py` (39, foundation — the step-4 memo's
M1–M5 on synthetic toys; teeth are PERMANENT in-process monkeypatch tests: offset-swap +
dropped-block red the assemble≡probe centrepiece on a same-size ASYMMETRIC grid,
M-ADJ-metric reds the reciprocity). VERIFIED: numerics+transport wall **1312/0**, ratchet
`transport:1` (numerics 0), sphinx -W exit 0, **elegance-enforcer PASS-WITH-NITS (0
violations; all 3 nits fixed pre-commit** — the SystemField docstring over-claim
corrected, collapse-trigger notes added to the two deliberately-coextensive helpers).
Key elegance verdicts banked: `CoupledField` NESTS the transport `Composite` (inter- vs
intra-system layers; both fan `op` to members, the math single-sources on the Field
leaves — collapse trigger = a 3rd direct-sum fan-out container); `eq=False` follows the
numerics `Field` convention (the `Composite` divergence needs no reconciling).

**FOLLOW-UP — FILED as #297 (2026-07-10, user-confirmed):** numerics carries THREE
member-contract concepts (`Vector` named / `_is_ravellable` ad-hoc private /
`SystemField`) — converge spellings 2+3 onto a named `Ravellable` Protocol home (the
`Vector` promotion precedent; `Vector` itself stays — the algebra/serialization split
is principled). Full scope + the ERR-053 interop-gate constraint on the issue.

**⚠ B.2b–d carried reminder (elegance-enforcer):** M2 (assemble≡probe) is really the
field-layout ≡ space-layout coextensiveness gate — on the toys `to_flat().size ==
prod(space.shape)` holds by contract; on the REAL members (`FullField`,
`RadialCharacteristicComposite`, multi-axis shapes) the three offset spellings
(`to_flat`, `system_slices`, `block_array` inference) are a LIVE thing to re-pin, not a
contract-guaranteed identity.

**NEXT = B.2b:** re-type `A_BA` (codomain) + `B_b` (domain/codomain) onto
`RadialCharacteristicComposite` behind thin transient FullField-gain adapters at the gain
seam (production bit-identical); mint the split-locus SourceSink leaves
(`RadialCharacteristicInteriorSourceSink` + boundary sibling, ANGULAR_RATE_UNITS, over
the locus field bases — the unified `RadialCharacteristicSourceSink(RadialCharacteristicField)`
is the pattern); P3 value-unmoved gates (re-point the `TestA_BA_SchurFold` /
`TestB_b_RayBoundary` probes). Open B.2b design points (main-agent proposes, user
steers): System B's member SPACE for the composite (a family-blind `FullFieldSpace.
from_blocks(interior=ray_interior_space, trace=ray_boundary_space)` vs a bespoke sibling)
and whether the source-role composite rides the SAME `RadialCharacteristicComposite`
class (role-erased slots, the FullField precedent / #289 F2 erasure) or a role-typed
sibling. Re-anchor from THIS block + `git log` (trust git).

### B.2b — RULED design points (2026-07-10, user; both = the recommended options)

- **DP1 (System B's member space): reuse `FullFieldSpace` family-blind.**
  `SNMesh.radial_characteristic_composite_space` (cached, presence-gated None) =
  `FullFieldSpace.from_blocks(interior_space=ray_interior, trace_space=ray_boundary,
  name="radial_characteristic")`. `from_blocks` gains the `name` param (identity stays
  honest — the space name signals the instance); `_rebuild` gains PRESENCE-DISPATCH
  (pass the `radial_characteristic` kwarg to `_recombine` only when the carrier exposes
  the slot — today it ALWAYS passes it, which would TypeError on a pure 2-block
  `Composite`; FullField behavior byte-unchanged). Zero new space classes; the
  member-wise metric dispatch stays single-sourced; this IS the post-eviction end-state
  (one composite-space class, instances differ in members), so (d) simplifies in place.
- **DP2 (source-role composite): same class, role-erased.** Re-bind
  `RadialCharacteristicComposite`'s static params from the flux leaves to the FIELD
  BASES (`…InteriorField`/`…BoundaryField` — the FullField precedent; the runtime guard
  already admits them; `⊖` already mints displacement members onto the same class).
  Role parses at consumption (B_b's #289-F2 parse moves onto the boundary member).
  `from_unified`/`to_unified` become ROLE-PRESERVING via an exact-class leaf table
  (Flux ↔ flux pair, SourceSink ↔ source pair, Displacement ↔ displacement pair) —
  one bridge body, no twins.
- **Units refinement over the brief's shorthand (Cardinal Rule 1):** the split
  DISSOLVES the unified SourceSink's documented corner-units deviation into two honest
  per-locus declarations — `RadialCharacteristicInteriorSourceSink.UNITS =
  ANGULAR_RATE_UNITS` (volumetric folded emission, like `AngularSourceSink`) and
  `RadialCharacteristicBoundarySourceSink.UNITS = ANGULAR_FLUX_UNITS` (trace-like
  corner datum, like `AngularBoundarySourceSink` — the "on the trace a source does not
  pick up the volumetric cm⁻³" convention).
- **Execution spine (each commit green):** b1 = split SourceSink leaves + composite
  re-bind + role-generic bridge; b2 = the space instance (`from_blocks` name +
  `_rebuild` presence-dispatch + the SNMesh cached prop); b3 = the A_BA/B_b re-type +
  solver-private transient gain adapters (`_lagged_gains` solver.py:619 + the
  `B_a + B_b` sum :251) + P3 re-points (`TestCoupledLift` isinstance gates → the
  adapter; the L4-S sentinel spy pins the wrapped CLASS method, so it fires through
  the adapter unchanged). Test-architect delta memo:
  `coupled_operator_b2b_retype_verification.md` (dispatched at arc start).

### B.2b — DONE @ `81c1995` (2026-07-10); recommended /compact point

Three green commits (`42b596d` b1 / `b15bf2e` b2 / `81c1995` b3 + the two plan/matrix
docs commits `427c159`, `c426b33`):

- **b1** — `RadialCharacteristicInteriorSourceSink` (ANGULAR_RATE) +
  `…BoundarySourceSink` (ANGULAR_FLUX — the trace convention; the split DISSOLVED the
  unified leaf's documented corner-units deviation); the composite re-bound to the
  FIELD BASES (honestly role-erased); `from_unified`/`to_unified` role-preserving via
  the exact-class `_UNIFIED_TO_SPLIT` table (+ its inverse; mixed-role and off-table
  refusals loud). +16 gates (G-b1.1–1.4 incl. the role⊕values split rows).
- **b2** — `FullFieldSpace.from_blocks(name=…)` + `_rebuild` PRESENCE-DISPATCH (the
  seed kwarg passed iff the carrier exposes the slot — FullField byte-unchanged; a
  pure 2-block carrier no longer AttributeErrors; ONE documented transitional cast,
  retires at (d)) + cached `SNMesh.radial_characteristic_composite_space`
  (identity `("radial_characteristic", (n_i + n_b,))`). +7 gates (G-b2.1/2.2: metric
  trio ≡ direct split-space application on the REAL composite, BOTH branch arms,
  role-preserving rebuild, the coextensiveness + flat-order pin).
- **b3** — **A_BA**: `apply → RadialCharacteristicComposite` (SOURCE members; no
  present-zero padding — "writes bulk" UNSPELLABLE), `apply_transpose` takes the
  composite cotangent → 2-block System-A FullField (ray=None); domain/codomain
  DECLARED (full_field_space / the composite space); fold engine unified behind the
  bridge (demote ruling). **B_b**: domain=codomain= the composite space; composite →
  composite (the #289-F2 role parse relocated onto the boundary member;
  `_reflect_corner` re-typed onto the split member — ONE law body);
  **ctor guard** (unconstructable seedless — the old None-pass-through + the dead
  seedless is_adjointable arm retired); `reflect_corner_inplace` non-None, bridges
  internally (the solver call site `_reflect_boundary_inplace` presence-guards).
  **Adapters** (TRANSIENT, retire at (d), each beside its block): `_RayEmission
  FullFieldGain` (radial_characteristic.py) + `_RayBoundaryFullFieldGain`
  (boundary.py) — byte-identical FullField faces; wired at `_lagged_gains` + the
  `_within_group_triple` sum; DELEGATION proven (the new G-b3.3 spy-count==1 gate
  keeps the L4-S sentinel non-vacuous). Re-points landed per the memo + two extra
  same-class sites the wall surfaced (test_invertible_operator's hand fold;
  test_si_single_primitive_contract's gain-1 isinstance → wraps-predicate);
  `_dense_seed` retired (successor `_dense_ray`, unified-layout probes).

**VERIFIED:** psi_half+g_adjoint 84/84; **FULL sn+transport+numerics not-slow wall
3457/0**; ratchet `transport:1` (numerics 0); sphinx -W exit 0 (verifies-marker
skips = the parked phantom). Every value row `array_equal` per the memo's §0 ruling
(pure re-labeling — no tolerance anywhere).

**NEXT = B.2c:** `build_coupled_system(sn_mesh, mat_xs) → (CoupledOperator,
CoupledSpace)` — the typed 2×2 grid over (full_field_space ⊕
radial_characteristic_composite_space) with A_AA/A_AB/A_BA/A_BB placed
by-construction (P1 alignment; P2 presence-structural: the seedless grid is 1×1);
SUBSUMES step 6 / task #34's presence collapse for the grid arm. The M2
assemble≡probe re-pin on REAL members (the elegance carried reminder) lands HERE —
CoupledSpace.system_slices vs FullField/composite to_flat on multi-axis members.
A_AB/A_BB stay unified-internal behind their block boundaries (V2/V4: A_BB is
already leaf-typed — its grid entry is direct; A_AB's re-type is a leaf swap).
Then B.2d (driver CoupledField + ray eviction + adapter retirement + protocol
2-block + ERR-053 re-size + E4 anchors). Re-anchor from THIS block + `git log`.

### B.2c — DONE @ 8b06052 (2026-07-10); recommended /compact point

Two green commits (`1097647` c1 / `8b06052` c2 + the plan/matrix docs commits). Delta gate
spec: test-architect `coupled_operator_b2c_builder_verification.md` (G-c1.1–1.3 +
G-c2.1–2.6; findings F1–F4; rulings R1–R4).

- **c1 `1097647`** — the grid-entry re-types. **A_BB** (`RadialCharacteristicOperator`):
  ALL FOUR action surfaces (apply/apply_transpose/solve/solve_transpose) re-typed
  composite → composite (domain=codomain=`radial_characteristic_composite_space`);
  solve gains the #289-F2 SOURCE-role parse; `inverse()` untouched (R4 — delegates,
  speaks composite for free). **A_AB** (`RadialCharacteristicSeeding`): domain → the
  composite space; codomain → `full_field_space` (apply emits the seed's bulk term as
  the FullField INTERIOR member over present-zero trace/ψ½ — System A IS its honest
  row space; apply_transpose reads `.interior`, discards trace/ψ½ structurally, emits
  the composite cotangent). Engines stay UNIFIED behind the bitwise bridge (demote
  ruling; dissolves at 4e). **Shared parse minted** (rule-of-three):
  `RadialCharacteristicComposite.require_member(x, mesh=, context=)` guards all THREE
  System-B block boundaries (A_BB + A_AB new, B_b converged — its transpose-path
  message no longer hardcodes ".apply"). Every value row `array_equal` (§0 re-label);
  re-points per the memo F4 ledger (TestA_BB_RadialBVP/Forward, TestA_AB_SeedInjection,
  test_ray_operator's L1 helper — orders untouched). NEW G-c1 gates: the four-surface
  container/identity-pin row (+ bridge-drop & value-corruption teeth) and the A_AB
  grid-entry container row.
- **c2 `8b06052`** — `orpheus/sn/coupled_system.py::build_coupled_system(sn_mesh,
  mat_xs, *, scattering_order=0) → (CoupledOperator, CoupledSpace)` (exported from
  `orpheus.sn`). Blocks: (A,A) `LC − S − B_a` stamped `SystemRole.A` (C-fwd);
  (A,B) `+Seeding` (loss sign internal); (B,A) `−Emission` (gain negated);
  (B,B) `A_BB − B_b`. P1 aligned-by-construction; **P2 presence predicate = the
  System-B member-space read itself** (the Optional narrow IS the branch; ratchet
  stays clean by construction, no cast). The dead-slot hazard DOCUMENTED + witnessed,
  not guarded (R3). The builder is a tracked transient construction twin of
  `_within_group_triple`/`_lagged_gains` (collapses at (d)). **PRESENCE FIX forced by
  `grid.H`**: `Emission.apply_transpose`'s cotangent ψ½ slot flips ray=None →
  PRESENT-ZERO (the transposed row sums it with A_AAᵀ's 3-block output — the
  mixed-presence law correctly refused the b3 2-block shape; G-b3.2 re-pinned
  None→present-zero; the gain adapter is immune, byte-identity preserved). B.2d flips
  ALL System-A emissions to 2-block together.
- **Gates landed** (`TestCoupledBuilder`, 7): G-c2.1 P1 + identity-pinned member
  spaces + the **F2 runtime-apply proof** (the unified and composite spaces COLLIDE on
  `(name, shape) = ("radial_characteristic", (ni+nb,))` ⟹ `FunctionSpace.__eq__` and
  hence the grid's per-block check CANNOT see a still-unified-typed block — runtime is
  the sufficient catcher; transitional aliasing, dissolves at 4e with the unified
  space); G-c2.2 P2 shapes + the four ctor-refusal bypass-proof; **G-c2.3 THE
  centrepiece** grid ≡ the COMPLETE fused loss (local helper WITH the B_b adapter on a
  REFLECTIVE sphere — F1: `_full_loss_case` omits B_b, vacuum masks it) at per-row
  bars (bulk rtol=1e-11 = the intrinsic M-M block split, the campaign's FIRST
  principled-equiv row in the block realization; trace/ray rtol=1e-12; dead slot
  array_equal 0) + teeth (off-diagonal swap UNCONSTRUCTABLE; dropped −B_b reds on
  reflective); G-c2.4 M2-on-real (`to_flat` ≡ `system_slices` + coverage + round-trip
  on real multi-axis members; block_array arm pinned UNAVAILABLE — F3/R2, re-enters
  only if an SN walk assembler wires streaming assembly); G-c2.5 forward block-`.H`
  reciprocity < 1e-12 on the real 2×2 + the M-ADJ-metric tooth (Mode-12/ERR-067);
  G-c2.6 the live-ray double-count POSITIVE witness (**REMOVE at B.2d**).

**VERIFIED:** psi_half 70/70 (63 + 7); operator neighborhood 148/148; FULL
sn+transport+numerics not-slow wall **3467/0** (5 skipped, 36 xfailed, 9:37); ratchet `transport:1`; sphinx -W exit 0
(verifies-marker skips = the parked phantom).

**DEFERRED to B.2d per R1/R2:** the A2a swap-law/inverse `coupled` fixture
(`grid.is_invertible` is False until step 5), E4 solve-anchors (φ=Q/Σ_t, k_inf — the
centrepiece's forward-equivalence to the anchored fused loss is the B.2c tie), the
dead-slot guard (structural at the eviction).

**NEXT = B.2d:** driver iterate → `CoupledField [ψ_A, ψ_B]` (W1 Mode-11 sentinel:
wrap the driver's route into `CoupledOperator.apply`, count>0) + **ray eviction from
FullField** (ATOMIC across carrier + space + protocol — explorer V1; `full_field_space`
collapses to 2-block, the coupled DOF count goes honest, the dead-slot hazard
dissolves, G-c2.6 witness REMOVED, G-b3.2/A_AB presence pins flip to 2-block) +
`CompositeField` protocol → 2-block + mixed-presence-law retirement + **the B.2b/B.2c
adapters + the builder's construction twin retired** (driver block-native; the triple
re-founds on the grid) + ERR-053 `restart=n_dof` re-size verification (CoupledField
conformance auto-tracks — pinned by the interop gate) + E4 anchors + the A2a coupled
fixture (with step-5's block solve if is_invertible lands there). Bridge the fused
walk at the SIX V2 signatures. Archivist docs follow-up: `discrete_ordinates.rst` (18
refs) + `loss_representations.rst` (8) describe the 3-block FullField as current.

### B.2d — RULED design points (2026-07-10, user; all four = the recommended options)

- **Spine d1→d2→d3 (each commit green):** d1 = the block-native driver on the B.2c
  dead-slot convention (iterate → `CoupledField[ψ_A 3-block, ψ_B]`; M = the fused walk
  behind a pack/split bridge; N = the coupled gain grid; the B.2b/B.2c adapters + the
  builder's construction twin RETIRE; W1 Mode-11 sentinel; ZERO touch to
  streaming/loss_rep — principled-equiv at ~ULP: rhs-assembly reassociation + GMRES
  dead-slot padding; the walk sees bit-identical inputs). d2 = THE ATOMIC EVICTION
  (FullField → 2-block; space + `CompositeField` protocol + mixed-presence law retire;
  ALL presence pins flip together; the SIX walk signatures re-type to leaf kwargs;
  `evaluate_residual` → coupled + the split Residual mint; G-c2.6 REMOVED; DOF count
  honest). d3 = E4 anchors + the A2a coupled fixture + the ERR-053 gate + archivist
  dispatch + estate.
- **DP-splitting: the System record.** ONE builder `build_within_group_system(sn_mesh,
  mat_xs, *, scattering_op=None, scattering_order=0)` in `orpheus/sn/coupled_system.py`
  returns frozen `WithinGroupSystem(loss, space, resolvent, gains)` — the NAMED
  `A = M − N` splitting (Hackbusch): `.loss` = the typed grid; `.resolvent` = M (the
  coupled-bridged fused walk on carrying meshes — pack `[ψ_A, ψ_B]` → 3-block via
  `to_unified`, fused `(L+C)` acts, split back; plain LC seedless); `.gains` = N (ONE
  `CoupledOperator [[S+B_a, ∅],[Emission, B_b]]` carrying — the ∅ (A,B) slot is
  structural: Seeding lives in M, the walk's welded feed, so `N = M − A` has no
  ray→bulk gain; `(S, B_a)` tuple seedless, B_a LAST per the existing boundary-gain
  convention the G-S arm parses). `_within_group_triple` + `_lagged_gains` RETIRE into
  it; the four solve sites consume the record; the solver's cached `scattering_op`
  injects via kwarg (a cache seam, not a flag). The fused `LC.apply` on a LIVE-ray
  3-block IS `M.apply` (welded Seeding feed + RC rows), so `matvec = M.apply − N.apply
  = A.apply` — consistent with the grid's complementary dead-slot arrangement.
- **DP-seedless: carrying only.** The coupled iterate appears exactly where System B
  exists (P2's spirit — a 1-member coupled wrapper on a seedless mesh is ceremony
  without algebra). Seedless paths (multi-D Cartesian G-S split, 2-D harmonic
  windowing, `_maybe_window`) are ZERO-TOUCH; carrying meshes are 1-D curvilinear ⟹
  always Jacobi, never windowed — the coupled arm bypasses both machineries
  structurally.
- **DP-Solution: own member.** Post-eviction `Solution.angular_flux` = the honest
  2-block System-A composite; `Solution.radial_characteristic:
  RadialCharacteristicComposite | None` carries System B's converged state as its own
  typed member (downstream System-A readers untouched; ray readers re-point).

### B.2d d1 — DONE @ `c0f23f6` (2026-07-11); recommended /compact point

One green feature commit + the matrix regen (`fa3ec11`). Delta gate spec:
test-architect `coupled_operator_b2d_driver_eviction_verification.md`
(G-d1.1–1.8 / G-d2.x / G-d3.x; findings F1–F5; rulings R1–R6 — F1–F5 all fired
and were resolved as ruled).

- **The record**: `orpheus/sn/coupled_system.py::build_within_group_system(...)
  → WithinGroupSystem(loss, space, resolvent, gains)` — THE single construction
  site; `_within_group_triple` + `_lagged_gains` retired into it (the builder's
  construction-twin status dissolved). Carrying: `resolvent =
  CoupledInvertibleOperator` (M — the fused walk behind the exact
  `_split_fused_state`/`_fuse_coupled_state` bridge; all FOUR surfaces incl.
  the reverse-scan pair; `inverse()` = `CoupledSweepOperator`, an
  InverseWrapMixin sibling, seed accepted-and-dropped per #282/2.5d);
  `gains = (N,)`, `N = [[S+B_a, ∅],[Emission, B_b]]` (∅ structural — Seeding
  lives in M; (A,A) stamped SystemRole.A). Seedless: `(L+C, (S, B_a))` pure
  re-package, G-S/windowing ZERO-touch. `build_coupled_system` = the loss-grid
  VIEW (delegation). Elegance verdict banked: the M−N split closes on SHARED
  piece objects (three views of one object set — cannot drift).
- **The four solve sites** consume the record; carrying iterate/q_ext =
  `CoupledField[ψ_A 3-block dead-ray, ψ_B]` split at the transitional birth
  seams (all marked retire-at-d2). The B.2b adapters DELETED — the fused
  spelling survives as TEST ORACLES (`FusedRayEmissionGain`/`FusedRayBoundaryGain`
  in `tests/sn/_test_helpers.py`; dissolve at d2 with the 3-block). R3 ruled:
  NO `.interior` face on CoupledField (the N-general container stays
  vocabulary-clean) — F1 (g1 spy) + F5 (`_flux_displacement_leaf` walks
  `systems[0]`) both re-pointed. Two #289-F2 role parses minted (compute_keff
  + finalize trace reads); the windowed arm narrows structurally.
- **The GMRES exact-breakdown carve-out** (`numerics/iteration.py`): the
  transitional dead-pad makes the coupled matvec SINGULAR; warm-started late
  outers break down AT the exact solution (final pr_norm literal 0.0) and
  scipy stamps info>0 — the warn branch now recognizes exact breakdown as
  convergence (a PERMANENT invariant; the ψ½ pad was the first caller, per
  the elegance NIT-2 reword). Root cause (the pad) dissolves at d2.
- **Gates landed**: `TestWithinGroupSystem` (11 — G-d1.1 W1 wrap-counter
  sentinel ×2 inner_solvers + the dead-slot rider + the fused-bypass tooth;
  G-d1.2 identity pins; G-d1.4 bridge round-trip + teeth; G-d1.5 N ≡ fused
  gains `array_equal` control + teeth, the dropped-B_b mutation UPGRADED to
  unconstructable by the all-None-column guard; G-d1.6 Mode-9 het-VACUUM
  2-region sphere same-fixed-point at SAFETY×tol ×2 solvers + the
  `array_equal` slab control; G-d1.7 the F5 diagnostic gate; G-d1.8 the
  ERR-053 d1 sizing pin). G-d1.3 executed as the 3-search audit (grep ZERO
  live references). Re-point ledger: 17 files (the ~12 F2 import re-points;
  the L4-S/L4-F structural pairs + the SI single-primitive contract re-pinned
  onto the record's block shapes; G-b3.3 deleted with its referent; the
  g-adjoint full-loss + psi_half centrepiece/ρ probes onto the oracles).

**VERIFIED:** FULL sn+transport+numerics not-slow wall **3477/0** (5 skipped,
36 xfailed, 19:09); psi_half 80/80; ratchet `transport: 1` (sn 0 — no casts,
no ignores); sphinx -W exit 0; smoke sphere k(SI) ≡ k(Krylov) @ 3.2e-11 with
the Solution ray intact; elegance-enforcer **PASS-WITH-NITS (0 violations)**,
NIT-2/NIT-4 fixed pre-commit, NIT-3 waived (defensible asymmetry, retires at
d2).

**CARRIED to d2/d3 (from the elegance + memo ledgers):**
- **NIT-1 (MUST, d3 archivist)** — 13 stale `_within_group_triple` doc refs
  (broken `:func:` roles render silently, no `-W` catch):
  `operator_algebra.rst` :1750 :5025 :5133 :5237 :5329 :5722 :6736 +
  `discrete_ordinates.rst` :12233 :13820 :14746 :14820 :14842 :17125 —
  re-point to `build_within_group_system` + fix tense; prioritize the
  present-tense ones. Plus the pre-listed 3-block-FullField sections
  (`discrete_ordinates.rst` 18 refs / `loss_representations.rst` 8).
- **Estate (d3)** — the pre-existing LC-spelling triplication: the builder +
  `SNSolver.self.L` (:849 init + :929 rebind_sigma_t) — collapse the two
  `self.L` legs onto the builder or a shared `_build_lc` primitive
  (`operator_algebra.rst:1748` documents it; the 3rd spelling IS the L-002
  collapse trigger).
- **NIT-2 anchor (SHOULD, test-architect at d2/d3)** — a minimal general
  singular-consistent GMRES test so the exact-breakdown guard stays exercised
  after the d2 pad removal.

**NEXT = d2 (THE ATOMIC EVICTION)** — one commit per explorer V1 + the memo's
G-d2.1–2.7: `FullField`/`TimedFullField` → pure 2-block (carrier hooks +
`zeros` keying + `__post_init__` seed checks); `FullFieldSpace` seed arms +
`CompositeField` protocol → 2-block; the mixed-presence law retires (6 raise
sites) + ALL present-zero producer pads (S/Sᵀ/F/B_a, A_BAᵀ ray-out,
`_radial_characteristic_scaled`); the SIX walk signatures re-type to leaf
kwargs (the CoupledInvertibleOperator bridge passes
`radial_characteristic_source/flux`/`seed_cot` explicitly — F4: the d2
re-point of test_282/walk-baselines/native-matvec IS the six-signature
catcher, values stay bit-identical); `evaluate_residual` → coupled + the
split Residual mint (R4 — rides the eviction atomically) +
`boundary_vs_interior_split` seed-gap closes; `Solution.angular_flux` →
2-block + NEW `Solution.radial_characteristic` member (DP-Solution);
G-b3.2/A_AB presence pins flip 2-block; **G-c2.6 witness REMOVED** (G-d2.7);
the timed presence law retires (test-migration only); the ERR-053 file
MIGRATES (F3 — coupled `to_flat` count; the honest end-to-end gate lands at
d3); honest-DOF gate G-d2.3 (`Δ == n_seed`); diffusion/CP walls = the
untouched bit-id oracle (G-d2.6); the birth seams + oracles + `_rebuild`
transitional cast dissolve. Then d3 (E4 anchors + A2a forward arm + the
honest ERR-053 gate + the archivist dispatch + estate). Re-anchor from THIS
block + `git log` (trust git).
Re-anchor from THIS block + `git log` (trust git).

### B.2d d2 — DONE @ `e5d1acf` (2026-07-11); recommended /compact point

THE ATOMIC EVICTION landed as ONE feature commit + the matrix regen
(`1bbd997`). FullField/TimedFullField/FullFieldSpace/CompositeField are pure
2-block; a live-ray ψ_A is UNREPRESENTABLE (G-c2.6 witness REMOVED — G-d2.7;
memo R3 discharged exactly as ruled: the type system is the guard).

- **The six walk signatures** speak EXPLICIT leaf kwargs: forward pair
  `(radial_characteristic_flux, radial_characteristic_source)` (read/filled
  swaps with direction — the (source, flux) pair is the fused factor's
  (rhs-side, state-side) ray legs), transpose pair `(seed_cot, seed_cot_out)`.
  No-legs on carrying = the ray-decoupled (A,A) block action via internal
  zero-substitution (bit-identical to the retired dead-slot arithmetic — the
  grid's OperatorSum members need no kwargs); joint solves REQUIRE the pair
  (no block-inverse spelling exists until step 5); seedless refuses. Guards:
  `_require_/_refuse_radial_characteristic` + `_require_leg_pair` (one
  family). The walk BELOW is untouched — the frozen `walk_matvec_*` +
  affine-carve baselines reproduce 0-ULP through the re-type (F4 discharged).
- **The bridge dissolved**: `_split_fused_state`/`_fuse_coupled_state`
  deleted; `CoupledInvertibleOperator`'s four surfaces = `_require_coupled_
  pair` parse + exact `to_unified`/`from_unified` round trips + walk-filled
  buffers (matvec/transposes → SourceSink buffer, solve → Flux buffer).
  `_system_b_member`'s live-ray refusal arm dissolved (unspellable).
- **Producer pads retired**: S/Sᵀ/F/B_a present-zero emissions, A_AB/A_BAᵀ
  codomain pads, `_zero_radial_characteristic_like`, and the C-operator's
  `_radial_characteristic_scaled` (+ its `invert=` boolean-flag smell) —
  the ray's σ term is A_BB's alone.
- **Driver native**: `_coupled_flux_state`/`_coupled_source_state` births at
  all four solve sites; `_build_fixed_source_rhs` returns the driver-ready
  state (accepts an explicit-q½ CoupledField); `evaluate_residual` coupled
  arm mints the SPLIT residual pair (`RadialCharacteristicInterior/
  BoundaryResidual`; unified leaf RETIRED) and the bare arm REFUSES a
  carrying mesh (the elegance CONCERN — the eviction's signature footgun, a
  DSA caller passing `solution.angular_flux` alone would silently drop r_B —
  closed loudly in-commit); `boundary_vs_interior_split`'s √(b²+i²)=‖r_A‖
  EXACT (the closed diagnostic gap); `Solution.angular_flux` 2-block + NEW
  `Solution.radial_characteristic` member with the presence biconditional
  (DP-Solution). Honest DOF: Δ == n_seed pinned (G-d2.3, merged into the
  psi_half G-d1.8 successor); ERR-053 file migrated (F3).
- **NEW production predicate beyond the ruled brief** (flagged): `SweepOperator.
  is_adjointable` gains the B.2d THIRD factor — a carrying-mesh fused wrap
  cannot thread `seed_cot`, so it refuses EAGERLY; `CoupledSweepOperator` is
  the joint swap-law home. The coherence file's sphere rows re-expressed as
  the COUPLED swap-law arm (G1c/G2c/G4c on M + the predicate row); its fused
  rows went seedless-only.
- **Gate re-shapes vs the d1-era spellings** (each with the successor's
  rationale in-file): G-c2.3 → `grid ≡ M − N` (two-path: explicit blocks vs
  welded legs; same per-row bars); G-d1.5 → `N ≡ pieces`; G-d1.4 → the M
  leg-plumbing gate (zero-leg ≡ no-leg; live-ψ_B non-vacuity; joint round
  trip on cells legs — corners are given-data slots); G-d1.6 het → SI ≡
  Krylov cross-driver (the fused reference is unrepresentable; the d1 proof
  vs pre-d1 stands as history); the d1.1 bypass tooth → the DP-seedless
  negative control (a fused bypass cannot be CONSTRUCTED); Mode-12 seed-flip
  closure re-expressed on the coupled M; the FusedRay*Gain oracles deleted;
  the carrier file REWRITTEN to presence-facts (System-B existence) + the
  V_cell metric rows (composite-member ≡ unified gauge).

**VERIFIED:** FULL sn+transport+numerics+diffusion+cp not-slow wall
**3722/0** (5 skipped, 36 xfailed, 17:31) — tests/diffusion + tests/cp GREEN
UNTOUCHED (G-d2.6, zero re-points); psi_half 78/78; ratchet `transport: 1`
(sn 0; no casts; ONE relocated `# type: ignore[override]` marker
FullField→TimedFullField, net ledger 0); sphinx -W exit 0 ×2; smoke: sphere-SI
+ both slab arms k BIT-IDENTICAL to pre-d2 HEAD (12 digits, worktree-probed),
sphere-Krylov Δ1.8e-10 = the honest-DOF Krylov-subspace change (anticipated
principled-equiv arm); elegance-enforcer **PASS-WITH-NITS, 0 violations**
(CONCERN + NIT-2 fixed in-commit; NIT-1 4× repack acceptable — collapse
trigger banked; NIT-3 seed_cot naming asymmetry noted).

**CARRIED to d3:**
- E4 anchors G-d3.1 (φ=Q/Σ_t reflective pure-absorber carrying sphere ≥2G,
  rtol≤1e-10, both drivers) + G-d3.2 (k_inf homogeneous reflective vs
  derivations closed form, paired with the het cross-driver row).
- G-d3.3 the honest END-TO-END ERR-053 gate (restart from the coupled ravel
  at both Krylov sites + the stall-deficit tooth) — the d2 migration landed
  the count pin; the end-to-end proof is d3.
- G-d3.4 A2a forward arm in the coherence file (grid.H reciprocity in the
  swap-law home) + the `xfail(strict=False)` swap-law inverse arm
  (`grid.is_invertible` False until step 5) — R5/R11.
- Archivist dispatch: the 13 `_within_group_triple` refs (d1-carried) + the
  3-block-FullField sections (`discrete_ordinates.rst` ~18 refs — the
  role-quadruple table row was minimally re-pointed to the split pair at d2,
  the SECTION rewrite is d3 — + `loss_representations.rst` ~8 presence-law
  refs) + the six-signature/leaf-kwarg protocol documentation.
- Estate: the LC triplication (`SNSolver.self.L` :~849/:~929 vs the builder)
  + the NIT-2 singular-consistent GMRES anchor test (the exact-breakdown
  guard's first caller — the dead pad — is GONE; the guard is now
  caller-less-but-permanent, needs its own minimal anchor).
- Pre-existing (NOT d2 regressions): audit MISSING ERR-042/045/047/051/055
  markers in untouched geometry/sweep files.

**NEXT = d3** — E4 anchors + A2a + honest ERR-053 + archivist + estate; then
Phase C (#47 4e walk un-weave, #41 step-5 block solve, #35/#34). Re-anchor
from THIS block + `git log` (trust git).

### B.2d d3 — DONE (2026-07-11); **B.2d COMPLETE**; recommended /compact point

Six commits: `1de9592` (the B_aᵀ fix) → `4841694` (LC collapse) → `d3fbfab`
(the d3 gates) → `f9d50f4` (theory docs) → `66d2244` (matrix) → `82a9db1`
(catches-markers) [+ the marker-driven matrix delta + this checkpoint].

- **THE FIND — the B_aᵀ vacuum diagonal leak (fixed @ `1de9592`).** The new
  A2a grid-reciprocity arm CAUGHT a live production bug on its FIRST run
  (defect 2.6e-4, het-VACUUM coherence sphere). Localization: the (A,A)
  quadrant only; dense discriminator → exactly 4 diagonal trace entries,
  expected (G⁻¹AᵀG) −1 vs grid.H −2 = LC's honest −1 PLUS a spurious +1 from
  ``B_aᵀ``. Root cause: ``_reflect_trace``'s transpose arm OUTPUT-projected
  ``lawᵀ`` onto the outflow rows — extracting the vacuum mask law's
  (zero-on-inflow ⊕ identity-on-outflow) harmless-in-the-forward identity
  block. The honest transpose of the row-projected forward ``B_face =
  P_inflow ∘ law`` is ``lawᵀ ∘ P_inflow`` (INPUT restriction, full image
  written). Off-diagonal permutation laws (reflective/albedo) are
  bit-identical under either spelling — why every reflective-fixture
  reciprocity gate was blind (the ERR-063 masking family), and why
  ``test_adjointable_when_all_faces_support`` had PINNED the masked regime
  (vv anti-pattern #12) — rewired onto the defining law ``got =
  lawᵀ(P_inflow·ψ)``. Leaf pins minted (psi_half TestBoundaryUnweld: vacuum
  ``Fᵀ ≡ T ≡ 0``; reflective dense ``T ≡ F.T`` bit-equal). The raw verb is
  now licensed by the ``adjointable()`` TypeGuard + a loud per-face
  ``MissingAdjoint`` (spec §39.1; ratchet-clean where the old ``getattr``
  was pyright-invisible), bite-tested via the stub-face fixture (the
  elegance do-now). Forward solves never call ``B_aᵀ`` ⟹ k-values untouched
  BY CONSTRUCTION (the wall's kinf gates confirm; no worktree probe needed).
  B_b audited clean (KIND-dispatched, vacuum zero both directions).
  **PROPOSED (user rules): catalog as a new ERR entry** — the vv
  error_catalog lives under `.claude/skills/` (not the main agent's to
  edit unprompted).
- **The gates (@ `d3fbfab`)**: G-d3.1 (psi_half
  ``TestWithinGroupSystemAnchors``) — flat-flux equilibrium END-TO-END,
  reflective pure-absorber carrying sphere, 2G group-graded (Mode-6), BOTH
  drivers, rtol 1e-10, per-ord ψ + φ + **the Solution.radial_characteristic
  member at the same equilibrium** (layout-free set check — the fold
  convention is exact on flat); teeth ride G-d1.5 per the memo. **G-d3.2
  discharged by CROSS-REFERENCE, not a twin** (Cardinal 2):
  ``test_kinf_homogeneous`` already runs production ``solve_sn`` on the
  carrying sphere × {1,2,4}G × both drivers vs the derivations closed form
  at 1e-10 (tighter than the memo's 1e-8) — the class docstring records the
  Mode-12 pairing (k_inf never credited against shape-class mutations).
  G-d3.3 — the END-TO-END site proof in the migrated ERR-053 file: a
  ``KrylovAcceleration.__init__`` spy through BOTH production drivers
  asserts every constructed ``restart`` == the SPACE-derived coupled sum
  (bulk + trace + System-B flat), Mode-11 fired-check + the bulk-only
  mutation tooth through ``_within_group_krylov``. G-d3.4 — the A2a
  loss-GRID forward ``.H`` reciprocity arm in the coherence file (<1e-12 +
  the metric-strip tooth) + the ``xfail(strict=False)`` swap-law inverse arm
  (dormant until step 5 flips ``grid.is_invertible`` — R5/R11; the wall's
  37th xfail).
- **Estate**: the LC triplication collapsed @ `4841694` —
  ``build_streaming_collision(sn_mesh, mat_xs)`` in ``coupled_system.py`` is
  THE one ``L + C`` spelling (the builder + both solver ``self.L`` legs call
  it; two dead solver imports retired; byte-identical extract-function;
  elegance-confirmed home + name). The NIT-2 GMRES exact-breakdown anchors
  (numerics ``test_iteration.py``): the real-scipy singular-consistent arm
  (exact solution, literal-0.0 tail, NO ERR-053 warning) + the
  stubbed-gmres guard-branch pair (info>0 + 0.0-tail must NOT warn /
  nonzero tail MUST) — **probe fact banked: current scipy's convergence
  test is ``<=``, so info>0-with-0.0-tail is unreachable via real solves on
  small systems; the stub pins the branch deterministically,
  version-proof.**
- **ERR-marker audit (qa @ `82a9db1`, every marker mutation-verified)**:
  5 closed — ERR-020 (the Phase-F marker-migration find: catcher re-homed,
  marker left behind — L-022 applied to markers), 031, 040 (trace-space
  layer; the envisioned BC-layer second defense never built), 051 (the
  catalog-authorized ΠR=4πI extension), 055 (file-level). **4 remain
  HONESTLY missing (041/042/045/047)** — each needs an UNBUILT BC-layer
  production invariant (never-raised error classes
  ``VacuumAppliedToOutgoingTraceError``/``ReflectionDidNotMapInflowToOutflowError``,
  the no-op ``assert_source_lives_on_incoming_trace`` ABC default, the
  weight-measure check that delegates to involution — ERR-042 would be a
  BLIND marker today). **PROPOSED follow-up issue** (module:geometry +
  module:tests, type:improvement) — user confirms before filing (GitHub =
  outward action).
- **Docs @ `f9d50f4`** (archivist, main-agent reviewed): the 13
  ``_within_group_triple`` refs re-pointed (2 deliberate historical
  literals); the route-(a) section reframed post-eviction; the
  loss_representations presence-law section → "How the walk sees ψ½ — the
  six-signature leaf-kwarg protocol" (read-vs-fill table verified against
  the code, incl. the code's own ``D_BB`` seed-diagonal notation); the
  FaceField admonition PLANNED → PARTIALLY-built (parent landed
  `4081c0d`); dev-history d1+d2 entries. In-review correction: the LC note
  names ``build_streaming_collision`` (the archivist wrote pre-collapse
  prose). **Archivist FLAG banked (untouched, follow-up)**: the
  dev-history route-(a) entry (~:20650) still says "zero-metric
  RadialCharacteristicSpace" — the refuted ERR-067 ghost-metric wording vs
  the corrected ``G_sd = V_cell`` body; an L-015 dev-history
  reconciliation pass.

**VERIFIED:** FULL sn+transport+numerics+diffusion+cp not-slow wall
**3733/0** (5 skipped, 37 xfailed, 11:21); ratchet `transport: 1` (sn 0 —
the TypeGuard bridge is the typed spelling of the fix); sphinx -W exit 0
(twice: archivist state + post-review state); audit MISSING 9 → 4;
elegance-enforcer **PASS-WITH-NITS, nothing blocking** (the one do-now —
the MissingAdjoint bite-test — landed in `1de9592`).

**NEXT = Phase C**: #47 (4e walk un-weave), #41 (step-5 block solve —
flips the A2a inverse arm live + retires the joint-solve pair requirement),
#35/#34 (docs + guard retirement). Local task #51 (B.2d) CLOSED. Re-anchor
from THIS block + `git log` (trust git).

## Phase C 4e — RULED design points (2026-07-11, user; P0 memos absorbed)

P0 grounding @ `9fdbd47`: explorer `campaign_4e_walk_unweave_map.md` (the six
signatures + the four welded regions + the overlap matrix + the unified-family
blast radius, all current file:line) + test-architect
`coupled_operator_4e_unweave_verification.md` (survivors ledger, S1/S2
single-source proofs, per-surface equivalence architecture, mutation ledger).
KEY P0 finding: the EXTRACT is smaller than the #47 task card — the A_BB
MATVEC rows are already single-sourced (`radial_characteristic_forward_
residual(+_transpose)`, no twin); the genuine welded twins are ONLY the two
SOLVE orchestrations (`_run:4207-4237` ≡ `A_BB.solve:498-522`,
`_run_transpose:4739-4762` ≡ `solve_transpose:564-597`, both WRAP-pinned
byte-identical today) + the FUSED A_AB composition (H1).

- **R-4e.1 — H1-NARROW scope.** Extract ONLY the two welded A_BB solve
  orchestrations (M.solve/M.solve_transpose ray legs route through the
  `RadialCharacteristicOperator` engine). The A_AB feed STAYS M's fused
  internal — consistent with the B.2d DP-splitting ruling ("Seeding lives in
  M"; the N grid's (A,B)=∅ is structural). Consequences: `_seed_rows_forward`/
  `_seed_rows_transpose` SURVIVE (matvec-side, shared-kernel-routed) ⟹ the
  TA's M1 headline re-aim becomes VERIFY-still-bites, not re-aim; the sphere
  forward AND transpose matvec stay bit-exact; e2's solve legs are expected
  bit-id-or-nearly (the WRAP gates pin welded ≡ engine byte-identical TODAY)
  — measure under `-W error::DriftWarning`, never assume; any drift is
  localized and explained or it is a bug. S2 tripwire form under narrow:
  post-e2 `loss_representation` has ZERO `carlson_inward_sweep_*` references
  (the engine + psi_half keep them; the walk's only carlson callers are the
  two welded marches). The cylinder non-carrying direct-seed FOLD is NOT
  A_BB — stays in the walk.
- **R-4e.2 — spine e1 → e1b → e2, each commit green.**
  e1 = THE WALK-SLOT REWRITE (atomic, like d2, forced by cells_view/
  corner_view coupling): the walk marches split-native (composite members;
  closure kwarg re-types unified→Interior leaf — same `.cells(p,-1)`
  accessor), the six signatures + 4 thin wrappers + `transport_sweep`
  re-type their leaf kwargs to the split composite, M's 4×(buffer +
  to_unified + from_unified) dissolve, driver factories mint split directly,
  `B_b.reflect_corner_inplace` takes the split boundary member, engines +
  Emission/Reconstruction go split-native (the ×9 + ×2 internal bridges),
  THEN the unified family retires (leaf base + Flux/SourceSink/Displacement,
  unified `RadialCharacteristicSpace`, `SNMesh.radial_characteristic_space`
  — presence readers re-point to `radial_characteristic_composite_space`/
  `_levels`, presence-coextensive by construction — `from_unified`/
  `to_unified` + tables, `from_angular_source` re-homed split-side).
  BIT-IDENTICAL GATE: every frozen baseline (`walk_matvec_{slab,sphere,cyl}`
  nulp=1, affine-carve ALL geoms, test_282 C(i)/C(ii), G2c array_equal)
  under `-W error::DriftWarning` — ANY drift at e1 is a representation bug,
  not an accepted delta (TA refutation 4).
  e1b = the FREED-NAME rename (own mechanical commit).
  e2 = the solve-leg extraction (route through the A_BB engine; delete the
  two welded marches; the WRAP twin-coexistence gates retire BY DESIGN —
  S1 engine-execution sentinel both drivers + S2 carlson-absence tripwire +
  row-6 oracle + test_282 C(i) cold <1e-11 carry the coverage; preserve H2
  ray-first ordering — `phi_aux = cells_minus` handoff reads the engine's
  returned flux member).
- **R-4e.3 — the freed name: `RadialCharacteristicComposite` → 
  `RadialCharacteristicField` at e1b** (user-ruled — the FullField mirror;
  the B.1d flagged convergence FIRES). The name REBINDS from the retired
  unified leaf to the composite ⟹ the retirement audit is STRICTER: at e1,
  every existing `RadialCharacteristicField` reference (code, tests, docs)
  is explicitly dispositioned dead-vs-composite BEFORE e1b reminst the name;
  isinstance checks that meant the leaf must not silently re-bind. Docs
  re-point at the archivist pass; archaeology collision accepted.
- **Step-5 boundary (do NOT cross):** no `grid.is_invertible` flip, the A2a
  inverse arm STAYS xfail, the joint-solve pair law survives; M's deeper
  re-pose as honest block-triangular substitution over the grid is #41.

### 4e e1 — DONE @ `63702e7` (2026-07-11); recommended /compact point

THE WALK-SLOT REWRITE landed as ONE commit (production + the 32-file test
migration), **bit-identical by measurement**:

- **Production**: the six walk signatures + 4 thin wrappers +
  `transport_sweep` take `RadialCharacteristicComposite` legs (same kwarg
  names); `_seed_rows_forward/_transpose` + the psi_half kernels re-typed
  to `(interior, boundary)` buffer pairs (same arithmetic, two view
  targets); `_run`/`_run_transpose`/`_apply_walk`/`loss_action_transpose`
  march member views; the closure kwarg → the Interior leaf; M's four
  surfaces pass members straight (the 4×(buffer+to_unified+from_unified)
  dissolved); driver factories mint composites
  (`_radial_characteristic_zeros`→`from_mesh`, `…source_from_per_ordinate`
  →`source_from_angular`, fission seed via the re-typed Reconstruction);
  `B_b.reflect_corner_inplace(composite)`; the A_BB/A_AB/A_BA/fold engines
  split-native (A_BB ctor guards on the interior space — the pyright
  narrowing fix). RETIRED: the unified leaf base + 3 role leaves (3
  modules deleted), the unified `RadialCharacteristicSpace` class,
  `SNMesh.radial_characteristic_space`, the bridge + `_UNIFIED_TO_SPLIT`
  tables, `_RayLeg.unified_key`. NEW composite factories:
  `source_zeros_on` + `source_from_angular` (the ONE fold factory,
  re-homed). The F2 space-identity aliasing dissolved structurally.
- **The migration agent's load-bearing finding (RNG faithfulness)**: the
  unified buffer interleaved cells-then-corner PER LEG; the shared test
  helper `random_radial_characteristic_composite` draws in that exact
  order, reproducing the retired RNG stream bit-for-bit — without it the
  frozen baselines would have "drifted" falsely. Plus one real catch: the
  `_install_forward_spy` 4-arg delegation would have TypeError'd on the
  6-arg kernel (the retired-symbol grep is blind to surviving NAMES with
  changed signatures — census lesson).
- **Retired test rows (each with named successor, per the migration
  report)**: the `TestSplitFidelityBridge` class (8 — the licence, now
  discharged by bit-identity), 2 split-space unified-interleave rows, the
  unified-gauge carrier row (re-posed structurally independent), the
  `_seed_to_coupled_layout` re-label, the b2c bridge-mechanism tooth
  (re-aimed on kernel output). Mode-12 ERR-067 seed-flip tooth re-aimed
  onto the pair return — VERIFIED still bites (>1e-6).

**VERIFIED:** frozen `walk_matvec_{slab,sphere,cyl}` + affine-carve at
nulp=1 ZERO drift (the e1 acceptance — TA refutation 4 satisfied: no
arithmetic fold); **k BIT-IDENTICAL pre-vs-post** on sphere/slab ×
SI/Krylov (probe vs a `009b4dd` worktree, 15 digits); FULL not-slow wall
**3721/0** (5 skipped, 37 xfailed; the −12 vs 3733 = the honest
retirements above); ratchet `transport: 1` (sn 0); sphinx -W exit 0
(theory-page prose refs to retired symbols render silently — the
archivist 4e pass is BANKED, not yet run).

**NEXT = e1b** (the freed-name rename `RadialCharacteristicComposite` →
`RadialCharacteristicField`, R-4e.3 — pure mechanical, own commit, strict
dead-vs-composite disposition of every existing reference BEFORE the
name reminted; grep-driven, ~code+tests+docs) **then e2** (the solve-leg
extraction, R-4e.1 H1-narrow: route `_run`/`_run_transpose`'s welded ray
marches through the `A_BB` engine, delete the twins; expected
bit-id-or-nearly — the WRAP gates pin byte-identity today — measure
under DriftWarning, never assume; S1 engine-execution sentinel both
drivers + S2 carlson-absence tripwire in loss_representation; the WRAP
twin-coexistence gates retire WITH the twins; test_282 C(i) cold
<1e-11 + row-6 + G-d3.1 carry the coverage). Then the archivist 4e docs
pass + matrix regen + the #47 close-out. Re-anchor from THIS block +
`git log` (trust git).

### 4e e1b — DONE @ `ea7f919c` (2026-07-11)

The freed-name remint landed exactly per R-4e.3: class →
`RadialCharacteristicField`, module → `transport/radial_characteristic_field.py`
(the `full_field.py` mirror, git-mv), `SNMesh.radial_characteristic_field_space`,
test helper `random_radial_characteristic_field`. Pre-remint disposition: the
freed name was fully DEAD in live code (13 prose hits — 6 docs + 3 code
comments + 4 test comments; zero isinstance / live bindings; the two `:class:`
xrefs carry the dead `fields._bases` path so cannot re-bind); 7 non-docs
archaeology sites clarified, 6 docs hits banked for the archivist pass.
Verified: targeted 150/0 under DriftWarning-as-error; **FULL tree wall 6338/0**
(18 skipped, 55 xfailed — the first whole-tree wall this campaign, strictly
broader than e1's scoped 3721); ratchet `transport: 1`; sphinx -W exit 0.
NEXT = e2 per the block above.

### 4e e2 — DONE @ `98fe2e36` (2026-07-11) — THE UN-WEAVE IS COMPLETE

The solve-leg extraction landed per R-4e.1 (H1-narrow): both welded ray
orchestrations DELETED from the walk; the production `(L+C)` path routes
System B through the NAMED resolvent.

- **Forward**: `_run` solves System B ONCE, up front, via
  `RadialCharacteristicOperator.solve` constructed over the walk's OWN σ_t
  (`CrossSectionField.from_mesh(sigma_t_gx, self.mesh)` — byte-faithful, the
  values array IS the walk's); the carrying branch reads the marched inward
  cells off the carrier as the M-M seed (ray-first preserved). **Transpose**:
  the reversed ordinate loop accumulates `psi_angle_bar` onto a COPY of the
  seed cotangent (the fused A_ABᵀ feed made visible — H1-narrow keeps it M's
  internal), then ONE `A_BB.solve_transpose(aug)` after the loop IS `m_seed`.
  Late (function-local) imports — module-level would cycle via
  `operators.streaming → loss_representation`.
- **carlson refs in loss_representation 8 → 0**; dead `dr` locals removed;
  the operator module's SOLVE-twin notes flipped terminal (the :712 A_AB twin
  note STAYS — step-5 scope under H1-narrow); test 369's prose refreshed.
- **Gates**: test_282 Mode-11 RE-AIMED onto `RadialCharacteristicOperator.
  solve` (CLASS-level wrap) = S1-fwd; NEW S1-transpose sentinel; NEW S2
  carlson-absence tripwire (all in test_282 §16.F). The engine's own
  two-leg-reference/WRAP/pole-continuation gates survive unchanged (they pin
  the engine, not the twin). Row-6 oracle (389) + round-trip (369) +
  C(i)/G-d3.1 carry the value coverage, all green through the engine route.
- **VERIFIED**: **k BIT-IDENTICAL pre-vs-post, 15 digits, sphere/slab ×
  SI/Krylov** (probe vs an `ea7f919c` worktree, confirmed importing the
  twin-bearing tree — 8 carlson refs); frozen walk/affine baselines nulp=1
  ZERO drift; **FULL tree wall 6340/0** (+2 = the new gates); ratchet
  `transport: 1`; sphinx -W exit 0; matrix 6600 → 6602.

**NEXT**: the archivist 4e docs pass (discrete_ordinates.rst ~20 unified-era
refs incl. the role-quadruple table + the 2 dead `:class:` xrefs to
`fields._bases.RadialCharacteristicField`; loss_representations.rst FaceField
diagram + the six-signature section; silent `:class:` refs to retired
symbols) + elegance-enforcer review of e1/e1b/e2 + matrix regen + the #47
close-out. Then #41 (step 5), #34/#35 per the Phase C order.

### 4e CLOSED @ `9c66fbe6` (2026-07-11) — task #47 done, STEP 4 (#40) COMPLETE

The docs pass + review landed: archivist re-pointed 13 dead xrefs across 5
retired module paths to ZERO (sources + built HTML — `-W` is blind to
silently-unresolved refs), redrew the role-quadruple table + FaceField
hierarchy for the split/composite architecture, disambiguated the freed-name
collision everywhere, recorded the un-weave in 6 places + the Phase C 4e
Development-history entry, and fixed one adjacent refuted-framing survivor
(the 2026-07-04 changelog's "zero-metric" → SPD `V_cell`; L-015 sweep gap).
Elegance-enforcer verdict on e1/e1b/e2: **ZERO violations, ship-quality**;
its 3 findings fixed in the same commit (the stale "bridge to the unified
layout" class docstring; 2 broken source-sink xrefs de-linkified).
**BANKED FOR STEP 5 (#41)**: the `RadialCharacteristicOperator` construction
triplication (`coupled_system.py:680` + the two walk legs) is coextensive
today — give it a reciprocal cross-ref or a shared construction seam when
the block solve lands. Agent IDs: archivist `a7b8d65faa208f705`, enforcer
`a946468abd9ec5742`. Sequence this segment: `ea7f919c` (e1b) → `eb3a2168` →
`98fe2e36` (e2) → `ffbd9761` → `9c66fbe6` (docs+close). NEXT = **#41 step 5**
(the block solve), then #34/#35. Recommended /compact point.

### The two d3 HELD decisions — RESOLVED (user rulings, 2026-07-11)

1. **B_aᵀ vacuum-transpose find → ERR-068 MINTED** in the vv error catalog
   (appended after the pending ERR-067; Mode-2 composition-order swap on a
   row-projected transpose, hidden by the ERR-063 masking family + an
   anti-pattern-#12 pin). Marker `@pytest.mark.catches("ERR-068")` landed on
   `TestBoundaryUnweld::test_b_a_vacuum_transpose_is_the_honest_zero` (born
   RED pre-fix @ `1de9592`, the teeth evidence). Per ownership, the
   `.claude/skills/` files (ERR-067 + ERR-068 + SKILL.md edits) remain
   UNCOMMITTED — the user commits the skills substrate.
2. **ERR-041/042/045/047 marker gaps → NO GitHub issue; user chose
   BUILD-THE-INVARIANTS-NOW**, queued as local task **#52** blocked by #41
   (step 5): mint the never-raised BC-layer error classes
   (`VacuumAppliedToOutgoingTraceError`/`ReflectionDidNotMapInflowToOutflowError`),
   the real `assert_source_lives_on_incoming_trace`, the independent
   weight-measure check; then catchers (mutation-verified RED) + markers.

### ⏸ COMPACTION — pre-step-5 checkpoint (2026-07-12)

**Everything committed through `69621ffd`** (the ERR-068 marker + d3-rulings
record); branch `refactor/sn-walk-unification` **154 ahead of origin/main,
pushes HELD**; scoped tree CLEAN. Standing exception on disk (survives
compaction, NOT the main agent's to commit): `.claude/skills/vv-principles/
{error_catalog.md, SKILL.md}` carry the user's pending ERR-067 + ERR-068
entries. Last verification state: full tree wall **6340/0** @ e2 (`98fe2e36`;
only docstrings/tests/docs touched since, targeted suites 130/0 + 116/0 +
82/0), ratchet `transport: 1`, sphinx -W exit 0, matrix current @ 6602 rows
(+1 ERR-068 catches row).

**NEXT = #41 (STEP 5 — THE BLOCK SOLVE).** Re-anchor from:

- **The gate spec summary** — THIS plan :174-176 (`TestCoupledSolve`):
  block solve ≡ current fixed-source (`Q/Σ_t`, non-fissile) AND k_eff
  (fissile `A/2g`); EXTRACT ≡ dense-LU (principled-equiv, NOT bit-id — the
  step-0 row-6 oracle `test_extract_to_dense_is_principled_equivalent_not_
  bit_identical` is the pinned reference); **the ρ-honest block residual
  `r = Aψ − q` via `evaluate_residual` as the STOP test, not `‖Δψ‖`**;
  cold-residual lag-death classifier. Fuller spec rows: the test-architect
  memo `.claude/agent-memory/test-architect/coupled_operator_step4_
  verification.md` (+ the 6 refutations at :182-190).
- **What step 5 flips**: `grid.is_invertible` False→True ⟹ the DORMANT
  `xfail(strict=False)` A2a swap-law inverse arm (G-d3.4, the wall's 37th
  xfail — R5/R11: `A.H.inverse() ≡ A.inverse().H`) goes LIVE; the
  joint-solve explicit-pair requirement (`radial_characteristic_source/
  _flux` kwargs threaded by every driver) can retire onto the block
  resolvent surface.
- **Banked seams to resolve IN step 5**: (i) the `RadialCharacteristicOperator`
  construction TRIPLICATION (`coupled_system.py:680` + the two walk legs in
  `loss_representation/__init__.py` `_run`/`_run_transpose` — coextensive
  today; give it ONE construction seam or a reciprocal cross-ref; enforcer
  4e note); (ii) the :712 A_AB "tracked transient twin" note in
  `operators/radial_characteristic.py` — its retirement decision (route the
  bulk rows through `A_AA + A_AB` or keep fused) is step-5 scope under
  H1-narrow.
- **After #41**: task **#52** (the BC-layer invariants for ERR-041/042/045/
  047, user-ruled build-now) unblocks; then #34 (step 6 — guard retirement,
  presence structural) / #35 (step 7 — verify + docs + close-out).
- **Standing constraints** (unchanged): surgical mode (main agent writes
  production; test-architect/explorer/qa/archivist/elegance-enforcer
  dispatchable); pushes HELD; canonical SERIAL pytest invocation
  (`.venv/bin/python -O -m pytest -p no:xdist --timeout=300
  -p no:cacheprovider`); ratchet oracle `python -m tests._harness.
  pyright_ratchet` baseline `transport: 1`; streamed LSP diagnostics =
  the #226 artifact (CLI is the arbiter); never `git checkout`/`restore`
  on uncommitted paths; commits stamp the session model.

Re-anchor from THIS block + `git log` (trust git, not the summary).
