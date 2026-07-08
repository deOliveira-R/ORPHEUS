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
  revert mutations); commits stamp **Opus 4.8**; sphinx -W must stay clean.

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
