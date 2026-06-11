# SN Sweep-Strategy Carve (C3.4 + C3.5) — polymorphic sweep/matvec dispatch

## ⭐ STATUS / RECOVERY (read first)

Part of the N-D layout campaign. **Authoritative index = `.claude/plans/nd_layout_foundation.md`**;
this file is the self-contained design for the sweep-strategy carve, which **re-scopes** the
original C3.4 ("wire a `_wavefront_1d_sweep` adapter") + C3.5 ("orchestration d-generic") around a
first-class `SweepStrategy` abstraction.

- **Designed 2026-06-10** in a multi-turn design conversation; every decision below is LOCKED
  (see §"Decisions locked").
- **⭐⭐⭐ LATEST (2026-06-11, S6.5 @ `1d86030` + S5.2 @ `f193b34` DONE; S6.9 MEASUREMENT
  DONE — window-fate decision PENDING USER).**
  **S5.2 (the post-S6.5 G4/G6 gates) @ `f193b34`:** NEW
  `tests/sn/solve/test_scan_march_end_to_end.py` — G4.a fixed-source Mode-9 FP-invariance
  (aniso+het+vacuum, non-flat guard) + G4.b all-reflective shed-order pin (ERR-056 class;
  source confined to the fuel half — the uniform source measured max/min=1.067 and TRIPPED
  the guard, config fixed not the guard; d=2 limitation stated, full diagonal stressor =
  d=3) + G6 eigenvalue with ScanMarch forced (closed-form k_inf 2G + SI≡Krylov het
  flux-shape). Forcing = ONE `default_for` patch (post-S6.5 all doors read it at call
  time), as a CONTEXT MANAGER (a fixture would force the reference leg too — caught at
  design), with explicit NON-VACUITY counters on `_sweep_interior`/`loss_action`. The
  config builders MIRROR (not import) the affine golden's private `_build_2d` —
  elegance-enforcer ruling: the principled Pattern-2 EXCEPTION (independent transcriptions
  keep config drift observable; a shared builder would let one edit silently move both the
  sha256 reference and the FP gate). 4/4 green; elegance PASS zero conditions.
  **S6.9/S5.3 MEASUREMENT (benchmark `derivations/diagnostics/diag_s69_scanmarch_vs_window_bench.py`,
  uncommitted per the diagnostics exclusion; numbers ALSO on #222): ScanMarch WINS
  everywhere measured.** Bare sweep sm/win = 0.57–0.84× (best at LS4 large grids 0.57–0.61;
  narrows with angular order: LS8 0.69–0.71, LS16 0.84); matvec 0.55–0.78×; peak memory
  IDENTICAL 0.98–0.99× (the rolling frontier has NO memory edge over the row-march at d=2;
  full-field is 1.3–1.4× both — the ORACLE niche is memory-irrelevant); end-to-end
  `solve_sn_fixed_source` 48×48 LS8 2G het = 0.82× (10.5 s vs 12.8 s — the kernel win
  amortized by scattering/moment/reflect overhead). **DECISION PENDING (Fork B2, USER
  call): flip the 2-D Cartesian default window→ScanMarch and retire `MovingFrontierWindow`?**
  Retirement surface if YES: `MovingFrontierWindow` + `_MovingFrontier` +
  `SweepDependencyGraph.walk_windowed` (the windowed storage walk — FullFieldWavefront
  uses `walk_full`) + the window≡full oracles re-target ScanMarch≡full; Fork-B2 discipline
  = regenerate the affine sha256 golden IN the flip commit ([[scan-march-verification]]
  §G5) + the G4 Mode-9 gates (ALREADY GREEN @ `f193b34`). CAVEAT pinned at S0: the d=2
  zero-copy contiguity argument was the window's edge and it did NOT materialize vs the
  row-march; the d=3 surface-layout question is moot until a 3-D compute path exists
  (neither rep has one today; ScanMarch generalizes scan(x)∘march(y,z), the window
  generalizes frontier_dim=d−1 — S4).
  **⏭ NEXT = the user's window-fate decision → if B2: the flip+retirement commit; then
  S5.5 (Sphinx architecture page for the FINAL state — wait for the fate decision so the
  page documents the survivor set), S6.6 (deferred), C3.6/C4/C5.**
- **PRIOR (2026-06-11, S6.5 DONE + COMMITTED @ `1d86030` — ONE representation
  instance; L21 is a TYPE FACT).** The four doors collapsed: (1) the APPLY door
  (`StreamingOperator.loss_representation` cached_property) IS the one instance;
  (2) the SOLVE door — `_solve_timed_full_field` now runs
  `self.loss_representation.sweep(rhs.bulk.values, ...)` DIRECTLY (NEW
  `InvertibleOperator.loss_representation` plain property → the streaming leaf's cached
  instance — deliberately UNcached, caching the forward would mint a second handle; the
  `AngularSourceSink` wrap-unwrap dropped — `from_mesh` stores BY REFERENCE +
  `_unwrap_source` returns `.values`, so the object flow is identical, bit-identical by
  construction); (3) the G-S resolvent threads
  `interior=self._invertible.loss_representation._sweep_interior`; (4) `_sweep_2d_wavefront`'s
  `interior` is REQUIRED (the fresh-`MovingFrontierWindow` default removed; the 2 direct
  test callers migrated to first-class `MovingFrontierWindow(mesh).sweep`). **DELIBERATE
  scope boundary:** the module-level `transport_sweep` REMAINS the operator-FREE functional
  entry (its one production caller = the `solve_sn` post-convergence reconstruction,
  no operator in scope, the `from_isotropic` /W projection load-bearing) — exactly TWO
  production `default_for` sites remain (the operator's cached property + transport_sweep).
  **NEW `tests/sn/operators/test_one_representation_instance.py`** (gate memo §4 adapted
  post-(f), spy-based call-time-self capture — both tests FAILED at pre-S6.5 HEAD
  (distinct ids / AttributeError), flipped PASS in the commit). Test migrations: 3
  `test_invertible_operator.py` spies re-pointed `transport_sweep`→`patch.object(CumprodScan,
  "sweep")` (their 3 reds caught the relocation); `cart2d_2g_nonsquare`/`het_operands`
  promoted to `tests/sn/_test_helpers.py` (2nd consumer). Gates: one-instance 5/5; bit-id
  set 49; eigenvalue+operators+solve 585 (SI≡Krylov≡k_inf); broad not-slow 2368/4/35
  (+2 = the new tests, xfail set unchanged); collection 5104; retirement audit clean;
  Sphinx baseline-only; elegance PASS zero conditions.
  **⏭ NEXT = S5.2 (task #17): the G4/G6 end-to-end scan-march gates**
  ([[scan-march-verification]] §4: G4.a FP-invariance on `si_2d_p1_aniso_het`-class
  aniso+het+vacuum config, G4.b reflective shed-order w/ level-symmetric + the d=2
  limitation statement, G6 SI≡Krylov≡k_inf with ScanMarch FORCED — post-S6.5 forcing =
  monkeypatch `loss_representation.default_for`, seen at call time by ALL doors).
  THEN S6.9/S5.3 (measure window-vs-scanmarch + the window-fate decision = USER decision
  on measured numbers), S5.5 (Sphinx architecture page, archivist), S6.6 (ExplicitMatrix,
  deferred), C3.6/C4/C5.
- **PRIOR (2026-06-11, S6.4 COMPLETE — ALL SIX SUB-STEPS (a)–(f) DONE + COMMITTED):
  the unified walk is REAL and the geography matches the algebra.** **(f) landed at `6e959b0`:**
  `sweep.py` DISSOLVED (byte-verbatim relocations, elegance diff-confirmed) — scan primitives →
  `spatial/scan.py`; orchestration (`transport_sweep` + `_sweep_1d_unified` + the schedule loop
  RENAMED `_sweep_scheduled` + `_sweep_2d_wavefront`) → `loss_representation.py` ORCHESTRATION
  section; the lazy-import CYCLE IS GONE; 16 importers + 4 spies re-pointed (3 spy failures
  caught = the false-green-after-relocation class); Sphinx roles repointed/de-roled +
  api automodule moved. **WavefrontFlux + InteriorFaceSpace RETIRED (user-approved)** with the
  succession recorded everywhere (the concept lives as `_MovingFrontier` per-level seed/shed +
  `FullFieldWavefront._octant_face_cochain`/`_edge_outflow`; the `wavefront-flux-cochain` theory
  anchor KEPT, derivation preserved as history). Gates: collection 5102 clean; broad not-slow
  sn+numerics+transport 2366 passed; fresh `-E` Sphinx zero non-baseline warnings;
  elegance PASS zero conditions. **⏭ NEXT = S6.5 — unify the doors → ONE representation
  instance. The POST-(f) DOOR INVENTORY (pinned 2026-06-11, HEAD `5112bab`):**
  (1) the APPLY door — `operator.py:1547` `StreamingOperator.loss_representation`
  cached_property → `default_for(self.sn_mesh)`;
  (2) the SOLVE door — `loss_representation.py:1613` `transport_sweep` →
  `default_for(sn_mesh).sweep(...)` (a SECOND fresh instance every sweep call);
  (3) the G-S resolvent — `solver.py:380` `interior=default_for(sn_mesh)._sweep_interior`
  (a THIRD instance per `_solve_scheduled` call);
  (4) the window default — `loss_representation.py:2369` `_sweep_2d_wavefront`
  `interior is None → MovingFrontierWindow(sn_mesh)._sweep_interior` (a fresh construction,
  bypassing `default_for` — the documented (b)-era seam).
  S6.5 = the operator HOLDS the one instance; `solve` (via `InvertibleOperator` →
  `transport_sweep` or threading the instance directly) and `apply` consume the SAME object;
  the resolvent + window-default plumbing collapse onto it. The discriminating one-instance
  test ([[s6-relayering-verification]] §4, designed at S6.0-prime — re-validate its door line
  refs against THIS inventory, the memo predates (f)) flips xfail→xpass. NOTE: representations
  are stateless frozen dataclasses (mesh-only field), so "one instance" is a STRUCTURAL
  type-fact goal (the L21 invariant becomes construction-enforced), not a perf fix — the
  per-shape DAG cache (c) already de-duplicated the heavy state. THEN S6.9/S5.3 (measure +
  window-fate), S5.5 (Sphinx architecture page), S6.6 (ExplicitMatrix, deferred), C3.6/C4/C5.
- **PRIOR ((a)–(e) detail):** Commit chain (NOT pushed): `2958aa1` (the S6.4 AMENDMENT locked: d-generic
  `_OctantWalk`, kernel-parameterized schedule loop, matvec→schedule, A2D-1 retire, NEW (e)+(f))
  → `7abad0e` **(a)** `_OctantWalk` + both matvec frames through it; A2D-1 RETIRED (successor =
  the window≡full MATVEC output oracle); NEW `tests/sn/operators/test_one_octant_walk.py` (SPY
  xfail-strict + AST anti-boolean tripwire + both-matvec-variants pin)
  → `1b4b0c0` **(b)** the sweep frames in: `_SolveOperands`+`_SweepEmit` (guarded
  angular-XOR-moment TYPE) + `_OctantWalk.sweep_group`; `_sweep_2d_scheduled` now
  kernel-parameterized (required `interior`; **ScanMarch gained G-S composability free**);
  `sweep_octant_group` + `_sweep_2d_scanmarch` RETIRED; **the one-walk SPY FLIPPED xfail→PASS**
  (`L.apply` and `A.solve` provably hit ONE `_OctantWalk._interior_walk`)
  → `b3b7a7d` **(c)** DAG ownership mesh→family: `SweepDependencyGraph.for_shape(shape)`
  (lru_cache(8), returns a cached `MappingProxyType` — mutation unrepresentable, identity
  stable); geometry.py SCRUBBED (the curvilinear None-slots + build site + import GONE);
  `_DAGWavefront.sweep_graphs` property; tests migrated
  `test_snmesh_sweep_graphs.py`→`tests/sn/primitives/test_dag_ownership.py`
  → `d1b2f02` **(d)** the full-field oracle folded: `FullFieldWavefront` kernels-only
  (`_octant_face_cochain` in-edge-only seeding ≡ the whole-trace ι_* seed BYTE-identically —
  upwind never-read-before-write); `_sweep_full_field` RETIRED; `_sweep_2d_scheduled` buffers
  d-GENERIC; ⚠ **WavefrontFlux + InteriorFaceSpace now have NO production consumer — fate =
  USER decision at (f)** (the typed ι boundary algebra dissolved into the shared frame)
  → `3438131` the HEALED cylinder twin-path xfails un-marked (pre-existing D-K heal, surfaced
  by the first full slow-suite run since; 3× green incl. a stashed-HEAD reproduction; bare
  asserts in test modules ARE live under `-O` — pytest rewrites them, probe-verified)
  → `27ed995` **(e)** the graph's 4 walks → `walk_full`+`walk_windowed` × `_CellSolve`(guarded)
  /`_CellResidual` level ops; diamond = PURE cell algebra (storage adapters + `SweepCellSlice`
  RETIRED; the CellUpdate extension point = the storage-free kernel pair); tests
  `test_cell_update_batch.py`→`test_cell_kernel_batch.py` + the kernel sha256
  source-of-record pin (the GENUINE hash exception) + two-walks-not-four; §15A.2 theory
  rewritten (archivist).
  EVERY sub-step byte-identity-gated (window≡full SWEEP+MATVEC `array_equal`, affine sha256
  golden, G2.c nulp, SI≡Krylov≡k_inf) + broad-suite green + elegance-enforcer PASS + Sphinx
  clean. ZERO "edit both in lockstep" markers remain; ONE octant frame; ONE O.4b block.
  **⏭ NEXT = S6.4(f)** (module geography; see the AMENDMENT bullet §S6.2–S6.6): `sweep.py`
  dissolves; `_x_scan_faces`/`_scanmarch_row`→`spatial/scan.py`; walk+schedule
  driver+`_sweep_1d_unified`+`transport_sweep`→`loss_representation.py`; RENAME
  `_sweep_2d_scheduled`→`_sweep_scheduled`; ~30 test files re-point + ~14 `transport_sweep`
  Sphinx roles (gate-memo ADDENDUM §(f) has the checklist); ⚠ ASK USER: retire or keep
  `WavefrontFlux`/`InteriorFaceSpace`. THEN S6.5 (one representation instance; the solver.py
  `default_for` plumbing + `_sweep_2d_wavefront`'s window default are its collapse targets).
- **PRIOR (2026-06-10, S6.3 DONE — the walk moved OFF the operator): S6.3a + S6.3b.**
  The architectural heart of S6.  **S6.3a** (Mode-8 prereq, commit `ee91db4`): migrated
  `test_g_adjoint_reciprocity` off 5 bare `assert`s → `pytest.fail` (the curvilinear
  sphere/cyl angular-second-triangular-factor reciprocity was a FALSE GREEN under `-O`;
  now fires, 12/12).  **S6.3b** (the walk-move): each `LossRepresentation.loss_action` now
  OWNS the `(L+C)ψ` walk — the 5 `_apply_*` bodies MOVED off `StreamingOperator` into
  `loss_representation.py` (1-D → `operator.M_spatial._compute_LpC`; 2-D
  window/scanmarch/full-field → the moved octant walks) and all 5 `_apply_*` DELETED (456
  lines).  `operator.apply = loss_action − σ_t·ψ.bulk` (Resolution A `L = (L+C) − C`; the
  `−C` glue COLLAPSED 5×→1×, a Pattern-2 win); `apply_transpose` mirrors (CumprodScan
  carries `closure.angular_adjoint`, reciprocity green).  **OUTPUT BYTE-IDENTICAL** (the
  `−σ_t` relocation is the same arithmetic, verified at a 150-test pre-deletion gate with
  the dead `_apply_*` still present + 182-test post-deletion + affine-carve golden +
  window≡full MATVEC oracle migrated to `loss_action` + scan-march G2.c + keff).  **A2D-1
  RETARGETED** (`inspect.getsource` `_apply_2d_cartesian` → `MovingFrontierWindow.loss_action`,
  hash `d18135c8…`→`b455b2de…`, provenance line added; output-identity via the oracle, NOT
  source).  **NEW convention pin** `test_loss_action_convention.py` (4✓): the
  NON-tautological structural anchor — flat-reflective `L·ψ_flat=0` ⟹ `loss_action=σ_t·ψ`
  (proves it is the FULL `(L+C)` loss, not bare `L`) + the `−C` glue cross-checked vs an
  INDEPENDENT `CollisionOperator`.  ⚠ caught + fixed an **S6.2 docstring corruption** (the
  file-internal `replace_all residual→loss_action` also hit docstring cross-refs to the
  SIBLING `graph.residual`/`residual_windowed` — 3 cosmetic spots, all fixed; lesson L25).
  elegance-enforcer **PASS-with-nits** (1 fixed: the Protocol/base `loss_action` docstrings
  still said "applies L" — a contract-surface bug-habitat for the future S6.6
  `ExplicitMatrix` leaf; now state `(L+C)ψ`).  Sphinx **build succeeded**; stale `_apply_*`
  doc-refs cleaned across 10 test files (only the A2D-1 historical provenance retains them).
  Commits: `ee91db4` (S6.3a) + `3a79ab3` (S6.3b).  **NEXT = S6.4 — RESCOPED to the unified
  kernel-parameterized octant walk (sweep≡matvec) + DAG ownership (committed `7bbeb5f`; see
  §S6.4 below — the two "S6.4" bullets).** The `test-architect` S6.4 gate plan **LANDED** at
  `.claude/agent-memory/test-architect/s6_4_unified_walk_verification.md` — **READ IT FIRST** at
  pickup. Key calls: **A2D-1 → RETIRE** (the source-hash is subsumed by the `window≡full` MATVEC
  output oracle + brittle on a shared walk; fallback = retarget to `_OctantWalk2D._interior_walk`);
  the **one-walk discriminating test** = a SPY on `_OctantWalk2D._interior_walk` hit by BOTH
  `L.apply` + `A.solve` (xfail now → xpass at staging-b); the **DAG-ownership test** =
  `not hasattr(SNMesh, "sweep_graphs")` + DAG-free reps grep-gate + the `geometry.py:343/350`
  None-slot GONE; ⚠ kernel = a KERNEL OBJECT + EMIT policy, NEVER a boolean `is_solve` flag (the
  memo §6 ships a structural tripwire); NON-SQUARE shapes (12×7, 5×9) are load-bearing (x↔y moat).
  **STAGING** (each bit-id-gated): (a) extract `_OctantWalk2D` for the 2 MATVEC variants → (b) bring
  `sweep_octant_group` + `_sweep_2d_scanmarch` IN [SPY flips] → (c) DAG-ownership → (d) fold
  `FullFieldWavefront` last. THEN S6.5 (unify the two `default_for` doors → the one-instance
  discriminating test flips xfail→xpass).
  ⭐ **AMENDED 2026-06-11 (user-approved, native-place architecture audit): the walk is d-GENERIC
  from birth (`_OctantWalk`, NOT `_OctantWalk2D`); `_sweep_2d_scheduled` becomes
  kernel-parameterized (ScanMarch gains G-S free); the matvec iterates `SweepSchedule.jacobi`;
  A2D-1 RETIRES at (a); TWO NEW SUB-STEPS (e) collapse the graph's 4 level-walks → 2
  direction-parameterized walks, (f) `sweep.py` DISSOLVES (module geography; the lazy-import
  cycle goes away). Full amendment text = the "S6.4 — AMENDMENT" bullet in §S6.2–S6.6 below.**
- **PRIOR (2026-06-10, S6 EXECUTION STARTED): S6.0-prime + S6.2 DONE.** The
  `test-architect` S6 verification plan landed (memo
  `.claude/agent-memory/test-architect/s6_relayering_verification.md`): the anchor set
  (A2D-1 hash byte-stable THROUGH S6.2 — its regen is an S6.3 event), the "one
  representation instance" discriminating test (fail-now/pass-at-S6.5; the two-doors
  `default_for` sites VERIFIED DISTINCT — `operator.py:1584` apply-door vs
  `sweep.py:241` solve-door, each constructs a fresh instance), the
  `loss_action`-returns-`(L+C)ψ` pin, and a NEW Mode-8 prerequisite it caught:
  `test_g_adjoint_reciprocity` asserts via bare `assert` (INERT under `python -O`), so
  S6.3's curvilinear angular-transpose coverage must be migrated to
  `np.testing`/`pytest.fail` BEFORE S6.3 leans on it. **S6.2 (rename) COMMITTED,
  bit-identical:** name CONFIRMED BY USER = **`LossRepresentation`** (`residual →
  loss_action`, `residual_transpose → loss_action_transpose`, handle
  `op.loss_representation`); module `sweep_strategy.py → loss_representation.py` (git-mv,
  history kept); `IncompatibleStrategy → IncompatibleRepresentation`, `SWEEP_STRATEGIES →
  LOSS_REPRESENTATIONS`. Bodies STILL delegate to `operator._apply_*` (NO behavior
  change; the walk-off-operator + the `(L+C)ψ` convention flip are S6.3). Gates GREEN
  `-O`: 108 directly-affected (A2D-1 hash + dispatch contract + scan-march G2.c +
  cumprod≡spine + window≡full) + 1172 broad (sweep/operators/solve/keff_2d) — 0 failed.
  `grep SweepStrategy` clean (2 historical phase-labels lowercased to "sweep-strategy
  carve"; the type is fully gone). Cross-domain-attacker memo: agentId
  `af87c1668054aaae7` (S6 plan) + `a8826036ef5b5e945` (S6 design). **NEXT = S6.3** (move
  the walk OFF the operator; `apply = loss_action − σ_t·ψ`; delete the 4 `_apply_*`;
  A2D-1 REGEN output-byte-identical via the window≡full oracle; do the Mode-8 reciprocity
  fix FIRST; dispatch elegance-enforcer — it earns its keep on the logic move).
- **PRIOR (2026-06-10): S5.0 + S5.1 DONE + committed; S6 CAPSTONE design LOCKED.** Read
  the **`# ⭐⭐⭐ S6 — THE OPERATOR / REPRESENTATION RE-LAYERING`** section at the BOTTOM of this
  file FIRST — it is the architectural keystone (the user's "why does the operator hold a matvec
  per strategy / apply IS a sweep action / where does the explicit matrix live" questions,
  cross-domain-attacker-validated across 4 expert frames). Commits since `954ddf4` (NOT pushed):
  **`4ca7a05`** S5.0 (ERR-057 denormal `ordinate_scan` fix) · **`8913229`** S5.1a (the `ScanMarch`
  sweep, `scan(x)∘march(y)`, oracle-verified) · **`0613298`** S5.1b (the row-march matvec twin +
  the `_x_scan_faces` Pattern-2 primitive). `ScanMarch` is a complete self-consistent strategy
  (sweep + matvec, both row-march; G2 + G2.c green; 2 elegance-enforcer PASS). **NEXT = S6** (the
  re-layering: operator = pure algebra, rename `SweepStrategy → SpatialRepresentation`, one
  instance both doors, `_OctantWalk2D` frame collapses the 2 Fork-B1 twins; the explicit matrix
  becomes a peer representation). S6.0-prime = dispatch `test-architect` FIRST. **The S5.3
  measurement (does ScanMarch beat the window?) folds into S6.9** (decide the window's fate on
  clean representations). Cross-domain-attacker memo: agentId `a8826036ef5b5e945`.
- **S0 ✅ DONE (2026-06-10):** `test-architect` verification plan delivered (memo
  `.claude/agent-memory/test-architect/sweep_strategy_carve_verification.md`). Key findings:
  `reduced is not None` ⟺ `is_1d` for ALL 4 constructible meshes (carve is behavior-preserving);
  `is_cartesian` is ORTHOGONAL → `supports` keys on BOTH; the anchor set (L2 `test_affine_carve_bit_identity`
  PRIMARY, L3 A2D-1 source-hash WRAP-not-RELOCATE, L5/L6 window≡full, L7/L8 value grounds); the Mode-8
  dispatch-pin hazard; the 2-D-adjoint deferral gate; the synthetic-d3 idiom.
- **S1 ✅ DONE (2026-06-10), bit-identical:** new `orpheus/sn/sweep_strategy.py` (Protocol + `_SweepStrategy`
  guard base + `_DAGWavefront` family + `CumprodScan`/`MovingFrontierWindow`/`FullFieldWavefront` thin
  wrappers + `Compatibility`/`IncompatibleStrategy` + `SWEEP_STRATEGIES` + `default_for`); `SNMesh.is_cartesian`
  (`curvature is None`); `transport_sweep` rewired to `default_for(sn_mesh).sweep(...)` (scattered branch retired;
  moment-on-1-D guard moved into `CumprodScan`); `test_unified_sweep_dispatch.py` MIGRATED (spy→selection+routing,
  `-O`-safe). Gates GREEN under `python -O`: **1622 passed / 0 failed** (86 anchor + 1164 broad sweep/operators/
  solve/regression/spatial + 372 eigenvalue/primitives/verification); the cumprod-spine equivalence stays
  `xfail` (S3 wires it). Regression drift = documented baseline (6920 ULP on 2d-aniso, within tol). elegance-enforcer
  PASS-with-nits (2 inline nits applied). Docs: inbound xrefs to `sweep_strategy` kept as plain literals
  (S5 adds the automodule + theory page; S1 stays code-only per the phase boundary). Committed `f6b4ad5`.
- **S2 DONE (2026-06-10), bit-identical:** the MATVEC twin now routes through the same polymorphic
  `sweep_strategy`. `SweepStrategy.residual(operator, psi)` + `residual_transpose(operator, phi)` added to the
  Protocol + base + 3 strategies (CumprodScan -> `operator._apply_1d`/`_apply_1d_transpose`; MovingFrontierWindow ->
  `_apply_2d_cartesian`; FullFieldWavefront -> `_apply_2d_cartesian_full_field`; the shared 2-D adjoint deferral
  lives on `_DAGWavefront.residual_transpose`). `StreamingOperator` gained a `sweep_strategy` cached_property
  (`default_for(self.sn_mesh)`); `apply`->`strategy.residual(self,psi)`, `apply_transpose`->`strategy.residual_transpose`;
  the 1-D bodies EXTRACTED verbatim into `_apply_1d`/`_apply_1d_transpose`; `_apply_2d_cartesian[_full_field]`
  UNTOUCHED (A2D-1 WRAP-not-relocate, hash green). TWO EVIDENCE-BASED DEVIATIONS from the plan's letter
  (both elegance-ruled JUSTIFIED): (1) KEPT the 3 `_MSpatial` defensive raise-guards
  (`_compute_LpC`/`_compute_LpC_transpose`/`_compute_decomposition`), NOT "remove all 4 raises"; L20 caller audit
  shows `_compute_decomposition` has non-strategy callers (`M_spatial.apply`/`M_angular_redist.apply`/transient
  orchestrator) so it is load-bearing; the others guard internal helpers; none are dispatch points. (2) the 2-D
  adjoint keeps raising `NotImplementedError` (NOT `IncompatibleStrategy`) - the mesh IS compatible (forward works),
  only the adjoint FEATURE is deferred, so `IncompatibleStrategy` would be the wrong concept. Gates GREEN under
  `python -O`: **1572 passed / 0 failed** (610 matvec-heavy operators/cartesian_2d/solve/eigenvalue + 962 broad
  sweep/regression/spatial/primitives/verification); PRIMARY `test_affine_carve_bit_identity` + A2D-1 source-hash +
  `test_sweep_vs_apply_consistency` + `test_g_adjoint_reciprocity` + `_compute_decomposition` paths all green.
  elegance-enforcer PASS (clean, no nits; twin-delivery smell did NOT materialize, matvec single-sourced). Sphinx
  clean. Committed `e08573e`.
- **S3 DONE (2026-06-10):** `FullFieldWavefront` is now the genuine d-generic verification SPINE, and the
  d=1 cumprod optimization is verified against it. **S3a (committed `37ce528`, d=2 bit-identical):**
  `SNMesh.streaming(axis)` (the d-uniform `2|μ|/Δ` accessor, pulled forward from S5);
  `_sweep_2d_full_field`->`_sweep_full_field` (sweep.py) + `_apply_2d_cartesian_full_field`->`_apply_full_field`
  (operator.py) generalized IN PLACE to any Cartesian d (`signs[:ndim]` in-plane projection of the
  full-angular octant label, `streaming(a)` over `range(ndim)`, `*spatial` shapes, `not any(signs)` pure-z
  guard dormant at d=1, ellipsis einsum); `FullFieldWavefront.sweep`/`.residual` call the d-generic bodies +
  `supports` widened to any-d Cartesian (overriding `_DAGWavefront`'s d=2-only). d=2 bit-identity anchored by
  the UNCHANGED window (`window ≡ full` oracle stays `np.array_equal` green -> the generalization preserved d=2;
  no twin). **S3b (uncommitted in tree -> committing now): the equivalence + adapter retirement.**
  `test_wavefront_cumprod_equivalence.py` REWRITTEN to strategy-vs-strategy (`CumprodScan(mesh).sweep` vs
  `FullFieldWavefront(mesh).sweep` at d=1, `assert_array_almost_equal_nulp` @ `_NULP_BOUND=128`, Mode-9
  anisotropic/het config); the hand adapters (`_wavefront_1d_sweep`/`_cumprod_1d_sweep`/`_SpineNotLanded`/
  `_spine_landed`/`test_diamond_difference_importable`) RETIRED (zero dangling refs); xfail->pass. ⭐ KEY: the
  d=1 cumprod ≡ spine equivalence HOLDS within nulp (the spine's whole-trace BC seed/absorb matches the cumprod
  chain seed at FP-association) -> the spine transitively inherits the analytical `k_inf=1.875` anchor.
  ⚠ The explorer's `(N,ng,nx,1)` phantom-y note was STALE (the C-phases already removed it); both strategies
  are rank-1 `(N,ng,nx)`, no layout bridge needed. Gates GREEN `-O`: S3a broad 1197 passed/0 failed + d=2
  window≡full oracle 8 + d=1 equivalence 4. elegance-enforcer **PASS (zero nits)** — streaming(axis)
  pull-forward JUSTIFIED, supports widening HONEST, two-full-field-bodies = the aggressive-retirement
  fuller-view-oracle exception (now DOUBLY pinned: d=2 window≡full bit-id + d=1 cumprod≡spine nulp). Sphinx
  clean.
- **S4 DONE (2026-06-10), d=2 bit-identical + d=1/d=3 admitted:** the rolling moving-frontier window — the
  storage-B PRODUCTION sweep optimization — generalized from the hardcoded 2-diagonal to the general
  `frontier_dim = d-1` rolling slab (a *point* at d=1, *line* at d=2, *surface* at d=3). **All d-dependent index
  arithmetic moved into a mesh-time `_FrontierPlan`** (`sweep_graph.py`: per-level slab `read`/`write`
  selectors + domain-edge `seed`/`shed` index maps; L16 zero per-sweep recompute), so the
  `apply_windowed`/`residual_windowed` walks are now DIMENSION-AGNOSTIC and read like the cochain trace algebra
  (`seed` ι_* / `incoming` gather / `emit` advance / `shed` ι*). `_MovingFrontier` → a `d`-tuple of slabs (free
  axes ghosted +1 on their own coord, the determined/last axis parity-rolled). The graph field
  `window_slots`/`window_edges` → single `window_plan`; `_window_metadata`/`_bounds`/`seed_x`/`seed_y` RETIRED
  (zero live refs). The walk SIGNATURES took per-axis TUPLES (`inflow`/`str_axes_octant`/`capture`) replacing
  the `inflow_x`/`inflow_y`/... pairs; the two orchestrators (`sweep.sweep_octant_group`,
  `operator._apply_2d_cartesian`) migrated to them (A2D-1 source-hash REGENERATED with a history line).
  ⭐ **d=2 BIT-IDENTITY preserved**: the per-level `read` selector degenerates to the SAME contiguous
  `slice(i0,i1+1)` (box-contiguous ⟺ `det <= 1` ⟺ `d <= 2`), so the generalized walk reproduces the legacy
  2-diagonal byte-for-byte AND keeps the measured contiguity speedup (PROFILED: window = **0.909× the
  full-field walk** at S8/64²/4g — still a NET SPEEDUP, win RETAINED, no d=2 fast-path needed). At d≥3 the
  anti-hyperplane is a SIMPLEX (not a box) → fancy index: the `(d-1)`-slab memory win holds; the d=3
  contiguity/SPEED is the one measured-cost question DEFERRED to profiling (no 3-D quadrature yet; OUT of the
  correctness gate). **Decision A resolved CONSERVATIVELY (a refinement of the initial widen-supports
  proposal): `MovingFrontierWindow.supports` STAYS `d==2` — its strategy entry (`_sweep_2d_wavefront` /
  `_apply_2d_cartesian`) is still the 2-D orchestrator (d-generic orchestration is C3.6, needs a 3-D
  mesh/quadrature). Widening would let `default_for` build it on a 1-D/3-D mesh whose `.sweep` then crashes in
  the 2-D orchestrator — a latent illegal-state. So NO `sweep_strategy.py` change: the WALK is general, the
  STRATEGY selects narrow honestly (the governing principle).** VERIFICATION:
  `test_sweep_graph_window_equivalence.py` REWRITTEN d-generically — `window ≡ full` at d=1 (frontier_dim==0
  point) / d=2 (`np.array_equal` bit-id anchor) / synthetic d=3 (B7 admission idiom, no quadrature); + a
  `(d-1)`-slab backing-size invariant. Gates GREEN `-O`: window≡full 28 + end-to-end affine-carve golden + A2D-1
  hash + nd-admission + cumprod-spine equivalence; **broad suite: sweep 617 passed/0 genuine-fail (the lone
  curvilinear-cylinder timeout is the pre-existing #206-area slow test, untouched by S4 — passes without the
  artificial 45s cap) + operators/solve 524 passed/0 fail**. elegance-enforcer **PASS-with-nits** (both
  forward-looking/non-blocking: the d=1 `is_point` base is the more-honest data structure — phantom-axis route
  REJECTED; `det` is the SSOT of the free/determined partition — one-line forward-note added). diamond.py
  kernel docstring de-staled ("d=2 PRODUCTION path" → "(d−1)-frontier").
- **S4 POLISH DONE + committed (`dc5d5c9`, bit-identical):** review flagged the `for a in range(d)` frontier
  loops as a latent-iterable smell. Profiling (cProfile, 64²/S8/4g) localized the regression delta to the
  per-level abstraction layer (NOT the shared kernel — which is the floor and CANNOT be vectorized: its
  left-fold `((σ_t+s₀)+s₁)` is load-bearing for d=2 bit-identity AND a `np.stack(...).sum(0)` is MEASURED 36%
  SLOWER, the bit-identity is FREE). Made the latent iterables explicit: `_MovingFrontier.emit` zips the three
  parallel per-axis collections (slabs ⨯ write-selectors ⨯ out_faces); `incoming` iterates the slabs directly;
  `seed`/`shed` store ONE axis-tagged record PER PRESENT EDGE (no `None`-padding — the walk iterates the edges,
  not the axes); the walk hoists the `(slice(None),*cell_idx)` selectors + builds `s_axes` via `zip`. BIT-
  IDENTICAL (only the iteration shape changed); measured d=2 window/full-field ratio **0.909×→0.862×**, ~192k
  fewer calls/walk; 529-test sweep/cartesian_2d/solve/keff_2d gate green. ⭐ **DEEP-DIVE → ISSUE #222 (the S5
  HEADLINE):** conceptually `apply_windowed` is forward-substitution on a lower-triangular operator; the
  anti-diagonal wavefront is ONE valid schedule, **row-march + x-scan is another** (VERIFIED ≡ wavefront @
  2.2e-15, **1.75× faster** single-octant) — and the within-row recurrence IS the first-order linear scan the
  1-D `CumprodScan` already uses (`ordinate_scan`), so it **UNIFIES 1-D `CumprodScan` + the 2-D window into one
  scan-march primitive**, adopts the flux-independent `a_attenuation` two-stratum cache the wavefront lacks
  (subsumes #206), and INHERITS the closed-form/pair-monoid conditioning robustness for free. The schedule is
  TWO orthogonal axes — **schedule** (anti-diagonal / row-march) × **backend** (forward-sub / closed-form-scan /
  division-free pair-monoid, the latter two dispatched by cumprod-underflow conditioning) — with
  `FullFieldWavefront` as the unconditionally-stable ORACLE the fast scan is pinned against (the S1–S4 pattern).
  ⚠ underflow-freedom is ALGORITHMIC not geometric: the "no internal boundary" geom change was BC-LAYER (the
  r=0 zero-area face + the `a=0` pole reset SURVIVE; ERR-054/#209 still live) and gradual cumprod underflow is
  intrinsic to the contractive recurrence — the pair-monoid backend already handles both. **NEXT = S5** (frontend
  `Compatibility` finalize; retire the d=2 orchestrators' hand-listed `str_x`/`str_y` onto a `streaming(axis)`
  axes-map + the `OctantLabel.sign_x/sign_y/streams_in_2d` 2-D shims — DEFERRED here because the orchestrators
  stay 2-D until C3.6; Sphinx theory page via archivist) — and **fold #222 in as the S5 design headline**
  (schedule × backend, wavefront-oracle-pinned), starting with a `test-architect` principled-equivalence plan.
- **Depends on C3.0–C3.3 (all DONE):** the wavefront spine is already dimension-generic — the
  per-octant DAG `SweepDependencyGraph.from_cartesian(shape)` (C3.1), the diamond cell kernel
  `cell_kernel_batch` (C3.2a), the full-field walk `graph.apply`/`graph.residual` (C3.2b), and the
  geometry octant build + `is_1d`/5-gate dispatch (C3.3). The spine is B7-verified at d=1 **and**
  d=3 (`test_sweep_graph_nd_admission`).
- **First implementation step is a `test-architect` dispatch** — this carve routes the *matvec*
  from one shape to another across a subsystem boundary, which is the proactive trigger in
  `subagent-handoff-protocol` ("operator-algebra carve crossing subsystem boundaries"). Do NOT
  start S2 before the verification plan lands.
- **Branch:** `worktree-sn-nd-layout` (off `main`); 33 commits ahead of origin, NOT pushed.

---

## THE PROBLEM THIS SOLVES

The 1-D-vs-multi-D sweep dispatch is scattered and procedural, and it is the *same decision*
spelled three different ways:

1. **`transport_sweep`** (`orpheus/sn/sweep.py:104`) branches on dimensionality:
   1-D → `_sweep_1d_unified` (the Blelloch `ordinate_scan`), 2-D → `_sweep_2d_wavefront`
   (the `_MovingFrontier` window).
2. **The matvec** branches in **5 operator gates** (C3.3: `not sn_mesh.is_1d` in
   `_compute_LpC`, `_compute_decomposition`, `_compute_LpC_transpose`, `apply`,
   `apply_transpose`).
3. **The verification oracle** (the full-field spine `graph.apply`) is reached only through
   **hand-built TEST adapters** — `_wavefront_1d_sweep` and `_cumprod_1d_sweep` in
   `tests/sn/sweep/core/test_wavefront_cumprod_equivalence.py` — because the spine is not a
   *selectable* path. There is no way to ask "which sweep methods are valid for this mesh?".

Adding a method, a dimensionality, or a frontend means touching all three. That is an
enum-style branch repeated at every call site — **cyclomatic complexity, not abstraction**
(`coding-elegance` anti-pattern: stringly-typed / special-case dispatch). An enum parameter
threaded into `transport_sweep` would only make it worse (a second branch axis).

---

## GOVERNING PRINCIPLE (the rule for every strategy)

> **Construct each strategy as general as its algorithm naturally allows. Select narrow.
> Specialize the implementation only on a *measured internal* performance cost.**

Three separable layers, never conflated:

- **Construct general (capability).** If the algorithm is *naturally* d-general, the CODE handles
  any d. A strategy that is d-specific in code *only because its algorithm is intrinsically
  d-specific* is legitimate (a prefix scan needs a total order → a chain → 1-D; there is no
  "2-D chain scan"). A strategy that is d-specific *because we wrote a narrow crutch* is a smell —
  it means the general structure did not surface.
- **Select narrow (policy).** Whether we OFFER / recommend / default a strategy at a given
  `(geometry, ndim)` is a SEPARATE layer (`supports` / `default_for`), independent of the code's
  capability. "Don't pick the window at d=1, pick cumprod" is a recommendation, *not* a reason to
  leave the window's code unable to express d=1.
- **Specialize on measured cost only.** The sole justification to restrict an implementation's
  d-range is a *measured* hot-path regression where the general construction makes the *kept*
  cases slower than a narrower construction would. Even then, keep the general path as the pinned
  fallback/oracle (the `feedback_aggressive_retirement` "fuller-view oracle" exception).

This principle is the acceptance lens for every phase below.

---

## THE ARCHITECTURE

### The strategy protocol (the algebra)

```python
class SweepStrategy(Protocol):
    def sweep(self, Q, sig_t, boundary, *, initial_guess=None) -> tuple[...]:  ...  # (L+C)⁻¹ q
    def residual(self, psi, ...) -> ...: ...                                        # (L+C) ψ   (matvec twin)
    @classmethod
    def supports(cls, mesh) -> "Compatibility": ...                                 # selection layer
```

Each strategy carries **both** the forward sweep AND the matvec twin, so the operator's `apply`
routes through `strategy.residual(...)` and the 5 C3.3 gates collapse.

### The hierarchy (what is shared vs not)

```
SweepStrategy (Protocol: sweep, residual, supports)
├── _DAGWavefront            ── reads mesh.sweep_graphs (the anti-hyperplane DAG) + the DD kernel
│   ├── FullFieldWavefront     buffer = full field     · Cartesian, any d · the ORACLE
│   └── MovingFrontierWindow   buffer = rolling (d−1)-frontier · Cartesian, d≥2 (built general) · prod opt
└── CumprodScan             ── reads the two-stratum scan cache · 1-D, any geometry · prod opt
```

**Substrate lives on the mesh, not the strategy** (strategies are lightweight consumers):
- `sweep_graphs` — the per-octant anti-hyperplane DAG (`from_cartesian(spatial_shape)`), built only
  for Cartesian meshes (C3.3); curvilinear sets it `None`.
- the two-stratum scan cache (`geom`+`coll` → `ordinate_scan`) — the chain-recurrence substrate.

**Is it the same DAG?** Yes for two of three: `FullFieldWavefront` and `MovingFrontierWindow`
consume the **same** `sweep_graphs` — they are two *buffer policies* over one anti-hyperplane walk
(full field vs rolling `(d−1)`-frontier), already pinned bit-identical by the C3.2b
`window ≡ full` oracle. `CumprodScan` builds no DAG — 1-D is a chain, the Blelloch closed form
needs no graph. The DD physics (closure, attenuation) is shared *conceptually* across all three but
expressed two ways (scan recurrence vs explicit kernel), pinned equivalent by the equivalence tests
— the accepted dual-view.

### Why these three (per the governing principle)

- **`CumprodScan`** — intrinsically 1-D (chain prefix scan needs a total order). Geometry-blind:
  slab, sphere, cylinder share one body via `_sweep_1d_unified` → `ordinate_scan` (the angular
  redistribution folds into the scan's affine source `b`, ordinates processed in angular order).
  Legitimately 1-D *by the algorithm's nature*, not by our hand.
- **`FullFieldWavefront`** — naturally d-general; ALREADY d∈{1,2,3} (the DAG + kernel + walk are
  d-generic, B7-verified d=1/d=3). The slow, readable spine; the verification oracle.
- **`MovingFrontierWindow`** — naturally d-general (the moving frontier is the `(d−1)`-dim rolling
  slab: a point at d=1, a line at d=2, a surface at d=3). **Must be constructed on
  `frontier_dim = d−1`**, capable of any d (S4). The 2-diagonal is the `frontier_dim == 1`
  instance, not the thing itself.

---

## SELECTION / COMPATIBILITY MODEL (the frontend-checkable layer)

Applicability is a **declared, queryable capability** — "make illegal states unrepresentable"
applied to method selection. The compatibility signal is the genuine criterion (the coordinate
system), NOT the `sweep_graphs is None` substrate proxy:

```python
# mesh.is_cartesian  ::=  curvature is None   (new one-line property, parallel to is_1d)

class CumprodScan:
    @classmethod
    def supports(cls, mesh): return Compatibility(mesh.is_1d, "requires a 1-D mesh")

class FullFieldWavefront:
    @classmethod
    def supports(cls, mesh): return Compatibility(mesh.is_cartesian, "requires Cartesian geometry")

class MovingFrontierWindow:
    @classmethod
    def supports(cls, mesh):
        return Compatibility(mesh.is_cartesian and mesh.ndim >= 2,
                             "requires Cartesian geometry, d ≥ 2")
```

`Compatibility(ok: bool, reason: str)` — the `reason` lets a teaching frontend gray-out a method
*and explain why* ("Moving-frontier window — requires Cartesian geometry, d ≥ 2"), which is
pedagogically load-bearing (ORPHEUS teaches reactor physics).

**Three consumers, ONE predicate (single source of truth):**
1. **Frontend** — `[S for S in SWEEP_STRATEGIES if S.supports(mesh).ok]`. A cylinder
   (non-Cartesian) → only `CumprodScan` → the dropdown shows only Blelloch.
2. **Factory default** — `SweepStrategy.default_for(mesh)` picks the best *available* production
   optimization, **falling back to the spine** when no optimization exists:

   | mesh | applicable | `default_for` |
   |---|---|---|
   | Cart-1D | `{FullField(oracle), CumprodScan}` | `CumprodScan` |
   | Cart-2D | `{FullField(oracle), MovingFrontierWindow}` | `MovingFrontierWindow` |
   | Cart-3D | `{FullField}` (window not built yet) | `FullField` ← never stuck |
   | Cyl/Sph-1D | `{CumprodScan}` | `CumprodScan` |

   The day a d=3 window lands (S4 + a 3-D mesh), Cart-3D's default flips from `FullField` to the
   window automatically — one predicate, no caller change.
3. **Construction guard** — `Strategy(mesh)` raises `IncompatibleStrategy` if `not supports.ok`, so
   even a bypassed UI cannot build an illegal pairing.

That `supports` predicate **is** the `is_1d`/`curvature` dispatch scattered across `transport_sweep`
+ the 5 operator gates today — now declared once per strategy. The whole point: "add 3-D window
support" becomes "extend one strategy, widen one predicate," not a hunt through every call site.

---

## MATVEC UNIFICATION (the C3.3-gate collapse)

The 5 C3.3 gates each spell `not sn_mesh.is_1d`. With `strategy.residual`, they collapse to
"ask the mesh for its strategy and delegate":

- `StreamingOperator.apply` (the live 2-D dispatch, `operator.py:~1458`) → `strategy.residual(...)`.
- The 4 raise-gates (`_compute_LpC`, `_compute_decomposition`, `_compute_LpC_transpose`,
  `apply_transpose`) → the strategy either implements the path or its absence is the
  `IncompatibleStrategy`/not-applicable signal (no hand-written "multi-D not wired" raise — the
  strategy's `supports` carries that).
- The operator holds `strategy = SweepStrategy.default_for(sn_mesh)` (selected once at
  construction) and the hot path is branchless `strategy.residual(...)`.

This is the operator-algebra carve crossing a subsystem boundary → **`test-architect` proactive
dispatch is mandatory before S2.**

---

## VERIFICATION STRATEGY

- **Bit-identical default gate (S1–S2).** `default_for(mesh)` must pick the SAME path the old
  dispatch did, so `transport_sweep` and the operator matvec are byte-for-byte unchanged on every
  existing mesh. The strategies in S1/S2 are *thin wrappers over the existing code* — pure
  dispatch refactor, zero algorithm change. Gate: the C3.2b/C3.3 bit-identity anchors
  (`test_affine_carve_bit_identity`, the DD regression snapshot, A2D-1 source hash) stay green.
- **Solve-vs-solve equivalence (S3).** Replace the hand-built adapters with strategy-vs-strategy
  through the real API: `CumprodScan(mesh).sweep(...)` vs `FullFieldWavefront(mesh).sweep(...)` at
  d=1; `MovingFrontierWindow(mesh)` vs `FullFieldWavefront(mesh)` at d=2; folding in
  `test_2d_full_field_oracle`. ONE pattern, parametrized by mesh ndim.
- **vv-principles Mode 9 (splitting/acceleration FP-invariance).** The equivalence tests run on an
  **anisotropic / heterogeneous** config (NOT the degenerate isotropic-reflective box), asserting
  the converged value is strategy-invariant to solver tolerance — a swapped optimization MUST NOT
  move ψ*. The analytical anchor (`k_inf = 1.875`, `test_cumprod_path_hits_analytical_kinf`) pins
  the converged value so the spine inherits it.
- **Principled-equivalence, not bit-identity, across strategies.** Cumprod vs spine differ at
  FP-association ULP (`assert_array_almost_equal_nulp`) — the existing `_NULP_BOUND` gate.
- **Synthetic d=3 admission (S4).** The general `MovingFrontierWindow` is verified `window ≡ full`
  on a synthetic 3-D shape (B7-style, no 3-D quadrature) to PROVE the `frontier_dim` construction is
  genuinely general — correctness, separate from the perf-quality question.

---

## PHASES (each independently bit-identical-gated + committable)

**S0 — `test-architect` dispatch (PROACTIVE, mandatory).** Verification plan for the carve:
which existing tests pin the legacy dispatch, the solve-vs-solve equivalence specs, the
bit-identity default gate, the vv-Mode-9 config, the compatibility/registry/guard tests. Output
shapes S1–S5.

**S1 — Strategy skeleton + 3 strategies as thin wrappers, SWEEP side.** `SweepStrategy` protocol +
`_DAGWavefront` base + `FullFieldWavefront`/`MovingFrontierWindow`/`CumprodScan` wrapping the
EXISTING sweep code (no behavior change) + `is_cartesian` property + `supports`/`default_for` +
`SWEEP_STRATEGIES` registry. Rewire `transport_sweep` → `default_for(mesh).sweep(...)`.
**Bit-identical** (each strategy wraps the path the old branch chose).

**S2 — Matvec side: `strategy.residual` + collapse the 5 gates.** Each strategy gains `.residual`
(wrapping `residual_windowed` / `graph.residual` / the 1-D matvec). Route `StreamingOperator.apply`
+ the 4 raise-gates through the strategy. **Bit-identical** (wraps existing matvec paths). Operator
holds `default_for(sn_mesh)`. *(Operator-algebra carve — S0 test-architect plan governs.)*

**S3 — Solve-vs-solve equivalence + retire adapters.** Rewrite `test_wavefront_cumprod_equivalence`
+ fold `test_2d_full_field_oracle` as strategy-vs-strategy, parametrized by ndim, on an
anisotropic/heterogeneous config (Mode 9). **Retire** `_wavefront_1d_sweep`, `_cumprod_1d_sweep`,
and the per-d oracle drivers `_sweep_2d_full_field` (`feedback_aggressive_retirement` +
`feedback_retirement_means_test_migration`). Measure the cumprod-vs-spine speedup through the
selector (the original C3.4 perf ask).

**S4 — `MovingFrontierWindow` `frontier_dim = d−1` generalization (the principle, embodied).**
Refactor `_MovingFrontier` / `_window_metadata` from hardcoded 2-diagonal to the general
`(d−1)`-frontier rolling slab. Verify d=1 (the trivial `frontier_dim == 0` base — free, no cost to
d≥2: a 0-dim frontier is a clean degenerate base, no hot-loop branch) + d=2 (bit-identical) +
**synthetic d=3 admission** (`window ≡ full` on a synthetic shape). **FLAG, do not assume:** the
d=2 speedup is a *zero-copy contiguity* win (the frontier *line* stored compactly is basic-slice
addressable; the full-field grid-diagonal needs fancy-indexed copies). Whether the general
`frontier_dim` layout preserves that contiguity at the d=3 *surface* is the ONE place the governing
principle's "measured cost" exception may bite — settled by profiling the d=3 frontier, not assumed.
If it loses the speedup, the fix is a measured d=2 contiguous *fast-path kept alongside* the general
frontier (pinned equivalent), never a d=1 exclusion.

**S5 — Frontend compatibility + cleanup + docs + the SCAN-MARCH headline (#222).** `Compatibility(ok,
reason)` finalized for the frontend. Retire the C3.2b elegance CONCERN — the `str_axes` hand-listed axis tuple
at the orchestrators → the strategy holds an `axes`-keyed `sn_mesh.streaming(a)` map (ONE axis-order
source). Retire the `OctantLabel.sign_x/y` / `streams_in_2d` 2-D shims. **Sphinx docs** for the
strategy architecture (dispatch the `archivist`): the protocol, the hierarchy, the compatibility
grid, the governing principle, the spine-as-always-available-capability story.

⭐ **HEADLINE = ISSUE #222 (scan-march unification — READ THE ISSUE).** The S4 deep-dive proved the d-D DD
sweep is forward-substitution on a lower-triangular operator; the anti-diagonal wavefront is ONE valid
schedule, **row-march + x-scan** another (VERIFIED ≡ wavefront @ 2.2e-15, **1.75× faster** single-octant,
reuses the 1-D `ordinate_scan` per line). This reframes S5 from "window vs scan as rival strategies" into
**two orthogonal axes**: (1) **schedule** = anti-diagonal-wavefront / row-march; (2) **backend** =
forward-sub / closed-form-scan / division-free pair-monoid (the latter two already dispatched by
cumprod-underflow conditioning in `ordinate_scan`). `FullFieldWavefront` stays the unconditionally-stable
ORACLE the fast scan-march is pinned against (the S1–S4 oracle-plus-fast-path pattern). The scan-march
UNIFIES `CumprodScan` + the 2-D window into one primitive, adopts the flux-independent `a_attenuation`
two-stratum cache the wavefront lacks (subsumes #206), and inherits conditioning robustness for free.
⚠ underflow-freedom is ALGORITHMIC not geometric (the "no internal boundary" change was BC-layer; the r=0
zero-area face + `a=0` pole reset SURVIVE — ERR-054/#209 live; gradual cumprod underflow is intrinsic; the
pair-monoid handles both). START with a `test-architect` principled-equivalence plan (the wavefront oracle
pins it at nulp; snapshots regenerate; a thick/long-chain conditioning test). SCOPE (a phase, not a polish):
boundary-outflow shedding + moment-output mode + the matvec twin (`residual`) + d-generalization
(`scan(x)∘march(y)`; the (y,z)-plane at 3-D) + per-line coefficients to keep the (d−1)-slab memory win.

*Mapping to the original plan:* S1+S2 = the "C3.4/C3.5 strategy + matvec together" the user chose;
S3–S5 subsume the original C3.4 (cumprod oracle + speedup) and C3.5 (orchestration d-generic,
str_axes fix, shim retirement). The original C3.6 (end-to-end 3-D — needs a 3-D mesh/quadrature)
remains downstream.

---

## DECISIONS LOCKED (2026-06-10 design conversation)

1. **Two selectable wavefront strategies** (`FullFieldWavefront` + `MovingFrontierWindow`) sharing
   a `_DAGWavefront` base — NOT one strategy with a buffer-policy flag. Rationale: the frontend
   should let a user pick the oracle explicitly for verification, and they have genuinely different
   memory profiles (the 3× peak-memory split C3.2b measured).
2. **The wavefront oracle is Cartesian-only.** Curvilinear is not an anti-hyperplane lattice; its
   verification reference is the per-cell `cell_update` march the scan is already pinned against.
   No apples-to-oranges oracle. Curvilinear meshes keep `sweep_graphs = None`.
3. **Selection keys on geometry first, via `is_cartesian`** (the genuine criterion), not the
   `sweep_graphs is None` substrate proxy. Honestly a 2-axis `(coord × ndim)` compatibility; the
   DAG family is gated on Cartesian, `CumprodScan` on the orthogonal `ndim == 1`.
4. **Strategies are constructed as general as their algorithm naturally allows** (the governing
   principle). `FullFieldWavefront` already d∈{1,2,3}; `MovingFrontierWindow` built on
   `frontier_dim = d−1` (any d); `CumprodScan` 1-D by the algorithm's nature, not by our hand.
5. **Selection ≠ implementation generality.** "Don't recommend the window at d=1" lives in
   `supports`/`default_for`; it does NOT justify a d=2-hardcoded window. The only implementation
   d-restriction allowed is a *measured* internal performance cost.
6. **The general path is never absent.** `FullFieldWavefront` (the spine) is the applicable
   fallback for any Cartesian d, so `default_for` is never stuck (Cart-3D defaults to the spine
   until a window is built). The code is never "incapable" at any Cartesian dimensionality.

---

## RISKS / OPEN

- **d=3 window contiguity (S4).** The one genuine "measured performance price" candidate — the
  zero-copy speedup may be d=2-line-contiguity-specific. Decided by profiling, not assumed. Memory
  win generalizes for free; speed win is the open question.
- **Bit-identity discipline (S1–S2).** The whole value of S1/S2 is that they are pure dispatch
  refactors. If any strategy's wrapper diverges byte-for-byte from the path it replaces, the carve
  has a hidden behavior change — the bit-identity anchors are the tripwire.
- **`residual` surface area (S2).** The matvec twin must cover forward + adjoint
  (`apply_transpose`) for each strategy. The curvilinear adjoint + the 2-D adjoint are partly
  deferred today (the C3.3 raise-gates) — the strategy must preserve those deferral boundaries
  (an unsupported adjoint is `IncompatibleStrategy`/not-applicable, not a silent wrong answer).
- **Frontend is hypothetical.** No frontend exists yet; `supports`/`Compatibility` are built for
  the factory + tests now, frontend-ready by construction. Do not build UI here.

---

## VERIFICATION GATES / ACCEPTANCE

- S1–S2: every existing sweep + operator + eigenvalue + solve + curvilinear suite green under
  `python -O`; the bit-identity anchors (`test_affine_carve_bit_identity`, DD regression snapshot,
  A2D-1 source hash) UNCHANGED; elegance-enforcer PASS (no enum/branch dispatch survives in a hot
  path; selection is single-sourced in `supports`/`default_for`).
- S3: solve-vs-solve equivalence green at the `_NULP_BOUND` on an anisotropic/heterogeneous config;
  `k_inf = 1.875` anchor intact; zero remaining `_wavefront_1d_sweep`/`_cumprod_1d_sweep`/
  `_sweep_2d_full_field` references (retirement audit).
- S4: `window ≡ full` bit-identical at d=2 (unchanged); d=1 + synthetic d=3 `window ≡ full` green;
  the speedup measurement recorded (d=2 confirmed; d=3 flagged pending or measured).
- S5: Sphinx builds clean; the `str_axes` + `OctantLabel` shims retired (retirement audit);
  no orchestrator hand-lists a positional per-axis tuple.

---

# ⭐⭐⭐ S6 — THE OPERATOR / REPRESENTATION RE-LAYERING (CAPSTONE; design LOCKED 2026-06-10)

> **READ THIS FIRST when picking up S6.** This is the architectural keystone of the SN
> operator algebra, surfaced by the user's two questions (2026-06-10): *"Why does the operator
> need `_apply_2d`? Why doesn't `apply` work for any dimensionality? Maybe the application is
> a sweep action (sweep using the operator)? If we had the explicit matrix, where would it
> live?"* The design below is **validated by the cross-domain-attacker across four independent
> expert frames** (numerical linear algebra, graph theory, category theory, algebra/representation
> — all four name the abstraction **"representation," not "strategy"**). It is LOCKED. The
> cross-domain-attacker memo (agentId `a8826036ef5b5e945`) is the structural backing.

## S6.0 — THE NATIVE FRAME (the "why")

**`(L+C)` (the within-group loss operator: streaming `L` + collision `C = σ_t⊙`) is a
lower-triangular `LinearOperator`.** Everything follows:

- **Lower-triangular under the sweep partial order.** Each cell depends only on *upstream*
  cells along an ordinate (`dag_walk_cell_indices` / `SweepDependencyGraph`). The per-octant
  DAG **IS the sparsity pattern** that triangularizes `(L+C)`. (Block-diagonal in the ordinate
  `n`; a *second* triangular factor in the angular index for curvilinear — the M-M/Carlson
  recurrence, `operator.py` `_compute_LpC_transpose` ~:583.)
- **The sweep is forward-substitution; the matvec is the row-action.** `solve = (L+C)⁻¹q`
  (sweep) and `apply = (L+C)ψ` (matvec) are the **two actions of ONE triangular operator**.
  They share the gather/scatter walk byte-for-byte and differ ONLY in the cell kernel
  (`cell_kernel_batch` solve vs `residual_kernel_batch` apply) — `diamond.py:367-376` is the
  proof, already true at the kernel layer (the L21 "one walk, two kernels").
- **The four walks are four *representations of one operator*, not four operators.**
  `CumprodScan` (1-D chain), `MovingFrontierWindow` (anti-diagonal, rolling frontier),
  `ScanMarch` (row-march), `FullFieldWavefront` (full-field oracle) are *schedules/buffer
  policies* over the SAME lower-triangular solve. **"apply the operator" IS "walk the
  dependency graph"** — a representation concern, not an operator concern.

### The smell (what's structurally wrong today)

Two defects, both instances of `coding-elegance` Smell #16 (two paths claiming to be one operator):

1. **Half-done dispatch.** `operator.apply(ψ) → strategy.residual(operator, ψ) →
   operator._apply_<specific>(ψ)`. The matvec *implementation* lives back ON the operator as
   one method per strategy (`_apply_1d`, `_apply_2d_cartesian`, `_apply_2d_cartesian_scanmarch`,
   `_apply_full_field`). The strategy merely forwards. **The dispatch was relocated, not
   eliminated — enum-dispatch wearing a Protocol's clothes.**
2. **TWO DOORS to the strategy.** `apply` reaches it via `StreamingOperator.sweep_strategy →
   default_for`; `solve` reaches it via `InvertibleOperator.solve → transport_sweep →
   default_for` — a SECOND, independent selection. So **"matvec ≡ sweep (L21)" is a coincidence
   maintained by `test_2d_full_field_oracle.py`, NOT a theorem enforced by construction.** The
   two halves of one operator's representation are selected on opposite sides of the
   operator/strategy boundary (capability asymmetry: `CAP_APPLY` on `StreamingOperator`,
   `CAP_SOLVE` on the composite `InvertibleOperator`).

## S6.1 — THE TARGET ARCHITECTURE (three layers)

```
Layer 1  OPERATOR ALGEBRA  (StreamingOperator / InvertibleOperator / OperatorSum)
  · L = (L+C) − C  factoring; composition (L+C−S−F/k); .H + G-metric (domain/codomain)
  · apply(ψ) = representation.loss_action(self, ψ) − σ_t·ψ.bulk      # the −C is the ONLY glue
  · apply_transpose(φ) = representation.loss_action_transpose(self, φ) − σ_t·φ.bulk
  · solve(q) = representation.sweep(self, q)                          # inverts the FULL (L+C)
  · HOLDS NO _apply_* methods.  ONE representation instance, reached by BOTH solve and apply.
        │  dispatches spatial application to ▼
Layer 2  SPATIAL REPRESENTATION  (rename of SweepStrategy)
  · The (L+C) operator's REPRESENTATION.  Protocol: sweep / loss_action / loss_action_transpose
    / supports.  TWO families:
      – matrix-FREE traversal:  CumprodScan · MovingFrontierWindow · ScanMarch · FullFieldWavefront
      – assembled SPARSE (S6.6, deferred):  ExplicitMatrix  (sweep = spsolve_triangular,
        loss_action = M @ ψ)
  · loss_action returns (L+C)ψ  (NOT L·ψ — see S6.3 confirmation); the operator subtracts C.
        ├── _OctantWalk2D (NEW base): the SHARED 2-D octant frame + the O.4b boundary block;
        │     MovingFrontierWindow + ScanMarch provide only `_interior_walk` (anti-diag | row-march)
        ├── CumprodScan          1-D total-order chain  (separate; geometry-aware via the cache)
        └── FullFieldWavefront   whole-trace WavefrontFlux oracle  (separate mechanism)
        │  reads ▼
Layer 3  STRUCTURE  (on the mesh)
  · SweepDependencyGraph = the sparsity structure (explicit DAG; carries update_batch / residual
    / residual_windowed); the matrix-free walks traverse it (or the implicit scan order)
  · DiamondDifference (cell_update) = the per-cell closure, two kernels (solve / apply)
```

### Locked decisions (cross-domain-attacker-validated, with the code grounding)

1. **`operator = pure algebra`; remove the four `_apply_*` methods.** `L = (L+C) − C`
   (`operator.py:1761`), `__add__`, `.H` + G-metric all stay on the operator.
2. **`loss_action → (L+C)ψ`, operator does `−C` — CONFIRMED right (not `L·ψ`).** Three code
   facts: (a) `residual_kernel_batch` natively yields `(Σ_t+Σsₐ)ψ̄ − (Q+Σsₐ·inₐ) = (L+C)ψ` at
   zero source (`diamond.py:363`); (b) every matvec subtracts `σ_t·probe` AT THE END
   (`operator.py:1762`); (c) the pure-z octant sets `LpC=σ_t·ψ` so `L·ψ=0` after subtraction
   (`operator.py:1731`). Returning `L·ψ` would force `−C` INTO every walk, re-coupling the
   collision diagonal into the streaming representation — the opposite of the factoring.
3. **The O.4b typed boundary-residual block lives in `_OctantWalk2D.loss_action`, NOT the
   operator algebra.** It is *trace bookkeeping of the walk* (reads `trace.outflow_indices`,
   writes `streamed − given`). The `−B` off-diagonal that IS algebra already lives as a sibling
   operator — leave it there.
4. **ONE representation instance reached by BOTH `solve` and `apply`** (closes the two-doors
   smell). Discriminating test: `assert StreamingOperator.sweep_strategy IS the instance
   InvertibleOperator.solve uses`. Today they construct two; after S6, one.
5. **RENAME `SweepStrategy → SpatialRepresentation`** (`residual → loss_action`,
   `residual_transpose → loss_action_transpose`). NOT a unify-after-two-instances defer: the
   benefit is ESTABLISHED (4 expert frames agree) and the 2nd instance already exists
   (`CumprodScan` + `MovingFrontierWindow` are already two representations); `SweepStrategy`
   actively mis-describes `ScanMarch.residual` (a matvec) and the future matrix (a solve).
   Rename the ABSTRACTION now (the carve is open; bolting it on later = re-touch every call
   site twice); defer only the matrix *implementation* (no consumer yet). Alt name considered:
   `LossRepresentation`. **Confirm the final name with the user at S6.2 execution.**

## S6.2 – S6.6 — THE STAGED REFACTOR (each independently gated + committable)

> ⚠ **PROACTIVE: S6.0-prime = dispatch `test-architect` FIRST.** This is an operator-algebra
> carve crossing the operator↔representation boundary (the `subagent-handoff-protocol` proactive
> trigger). The plan it returns pins the bit-identity gates (A2D-1 source-hash, `window≡full`
> oracle, `test_affine_carve_bit_identity`, `SI≡Krylov≡k_inf`) + the NEW "one representation
> instance" discriminating test BEFORE any code moves. Per the A2D-1 pitfall: the carve WILL
> regenerate the source hash — refresh it in the SAME commit and assert *output* byte-identity
> via the `window≡full` oracle (not source identity).

- **S6.2 — Rename + the protocol reshape.** `SweepStrategy → SpatialRepresentation`;
  `residual → loss_action` (return `(L+C)ψ`); `residual_transpose → loss_action_transpose`.
  Update the Protocol, the 4 concrete classes, the registry, and all call sites. The
  `loss_action` bodies still delegate to the operator's `_apply_*` for now (pure rename, no
  behavior change — bit-identical). **Confirm the name with the user.**
- **S6.3 — Move the walk OFF the operator.** `operator.apply(ψ) = representation.loss_action(
  self, ψ) − σ_t·ψ.bulk`. Move each `_apply_*` body INTO its representation's `loss_action`
  (returning `(L+C)ψ`); the operator's `−C` subtraction is the only remaining glue. Delete the
  four `_apply_*` methods from `StreamingOperator`. Regenerate A2D-1 (output-byte-identical via
  oracle). Same for `apply_transpose` → `loss_action_transpose` (1-D wired; 2-D deferred raise
  preserved; ⚠ the curvilinear angular second-triangular-factor — `closure.angular_adjoint`,
  `operator.py:727` — MUST ride with `CumprodScan.loss_action_transpose`, else a spatial-only
  reverse silently drops the angular transpose).
- **S6.4 — ONE kernel-parameterized octant walk shared by `sweep` + `loss_action`
  (RESCOPED 2026-06-10 — user insight: "matvec and sweep are ONE walk, differing only in the
  cell kernel + the target; long-due debt").** The DEEPEST L21 realization. matvec and sweep
  walk the SAME octant DAG / SAME frontier / SAME boundary seed-shed — the face recurrence
  `out = 2ψ̄ − in` is LITERALLY identical. They fork ONLY at:
  (i) the **CELL KERNEL** — `diamond.py` `cell_kernel_batch` (SOLVE: ψ̄ UNKNOWN, divide by
  `D = σ_t + Σs`, uses source `q`; :273) vs `residual_kernel_batch` (APPLY: ψ̄ GIVEN,
  pure-reflection `α = −1` no divide, evaluate the residual; :331); and
  (ii) the **EMIT policy** — sweep → ψ + scalar/moment projection; matvec → residual bulk + the
  O.4b boundary defect.
  The graph ALREADY proves it: PAIRED methods `graph.apply`/`apply_windowed` (sweep, :682/:849)
  ↔ `graph.residual`/`residual_windowed` (matvec, :770/:929) walk ONE graph, forked at the
  kernel. ⭐ `ScanMarch` ALREADY embodies the unification — its `sweep` + `loss_action` both call
  the ONE primitive `_x_scan_faces` (sweep.py:1186), differing only in `(α, β)`. The DAG family
  does NOT yet: the octant orchestration is written THREE times in lockstep
  (`sweep.py::sweep_octant_group` sweep + `MovingFrontierWindow.loss_action` +
  `ScanMarch.loss_action` matvec); the codebase ADMITS it (`sweep.py:1328` "edit both in
  lockstep" + the 2 Fork-B1 IOU notes).
  **THE CARVE: extract `_OctantWalk2D` — the ONE octant traversal (octant projection + pure-z +
  frontier + boundary seed/shed + the O.4b boundary block), parameterized by a CELL-KERNEL object
  (`apply_windowed` | `residual_windowed`) + an EMIT policy — NOT a boolean flag.** `sweep` and
  `loss_action` STAY the two NAMED public faces (the q-vs-ψ applications) but both delegate to the
  ONE walk. `sweep_octant_group` + the matvec `loss_action` bodies COLLAPSE; ALL THREE "edit both
  in lockstep" twins DISSOLVE; "matvec ≡ sweep" becomes STRUCTURAL, not test-maintained. Also
  HARMONIZE the signatures (today `loss_action(operator, ψ)` reaches the operator for σ_t while
  `sweep(Q, sig_t, boundary)` takes raw arrays — both should take the SAME problem data).
- **S6.4 — the CLASS STRUCTURE (three orthogonal concerns, ONE design; AMENDS S6.1 Layer 2+3).**
  (1) `_OctantWalk2D` = the shared sweep+matvec traversal (kernel + emit params, above) + the
  O.4b boundary block, which is duplicated 3× today (window + scan-march + full-field) → it lives
  HERE, shared by all three.
  (2) **DAG ownership** (user-raised; AMENDS S6.1 LOCKED Layer 3): move the DAG OFF the mesh
  (`geometry.py:350 self.sweep_graphs = None` — the illegal-state smell) INTO the `_DAGWavefront`
  family {window, full-field}, lazily built + **cached per mesh-SHAPE** by an accessor the family
  OWNS (not the mesh), so DAG-free reps (`CumprodScan`, `ScanMarch`) NEVER mention
  `sweep_graphs`/`OctantLabel` and the mesh reverts to pure geometry. VERIFIED only `_DAGWavefront`
  uses it (`sweep.py:910/1533` + `loss_representation.py:524/672`). S6.3 set this up (the
  window/full-field `loss_action` now live in the family's module; `sn_mesh.sweep_graphs` →
  `self.sweep_graphs`; the sweep `_sweep_2d_wavefront`/`_sweep_full_field` fold into the shared
  walk or take the DAG as a param).
  (3) **The composition**: `MovingFrontierWindow` ∈ BOTH `_OctantWalk2D` AND `_DAGWavefront`
  → mixin/composition, NOT single inheritance; `ScanMarch` ∈ `_OctantWalk2D` only (DAG-free
  row-march); `FullFieldWavefront` ∈ both (DAG + the shared boundary block, different interior
  `WavefrontFlux` walk).
  Each sub-step independently bit-identity-gated (window≡full sweep + matvec, scan-march G2.c,
  affine-carve golden, keff). **AMENDS S6.1: Layer 2 WALK becomes kernel-parameterized (shared by
  sweep + matvec); Layer 3 STRUCTURE moves from "on the mesh" → "owned by the DAG-using family."**
- **S6.4 — AMENDMENT (LOCKED 2026-06-11, native-place architecture audit; user-approved).**
  The audit confirmed the direction fork (solve vs apply) bottoms out in the DISCRETIZATION —
  `diamond.py::cell_kernel_batch` (solve) / `residual_kernel_batch` (apply), pure + storage-free
  + d-generic — and is then RE-SPELLED at three altitudes: the octant frame (6× = 3 reps × 2
  directions; the original S6.4 target), the graph level-walks (4× = `apply`/`residual` +
  `apply_windowed`/`residual_windowed` — the frontier seed/incoming/emit/shed loop written
  twice, the full-field level loop twice), and the cell storage adapter (2× =
  `update_batch`/`residual_batch`, already thin via the factored gather/scatter). The NATIVE
  architecture spells the fork ONCE — at the representation surface (`sweep` vs `loss_action`)
  — threaded down as a (cell-kernel, emit-policy) pair; everything between surface and cell
  algebra is direction-agnostic traversal. (ScanMarch spells the same fork as scan
  COEFFICIENTS — solve `α = 2s/D−1` vs reflection `α = −1` — the Blelloch closed-form
  evaluation of the same WDD closure, pinned by G2.c; a second legitimate spelling, in
  `spatial/`.) The LOAD-TIME IMPORT CYCLE (loss_representation imports sweep's bodies; sweep
  lazily imports `default_for` back) is the module-geography smell: bodies live apart from
  their owner. Decisions LOCKED:
  (1) the walk is **d-GENERIC from birth** — `_OctantWalk` (NOT `_OctantWalk2D`): signs/faces/
  inflow/captures are per-axis tuples over `mesh.ndim` (folding the d-generic oracle into a
  2-D-hardcoded walk at (d) would regress C3); byte-identical at d=2.
  (2) `_sweep_2d_scheduled` becomes **kernel-parameterized** (representation supplies the solve
  interior) → `_sweep_2d_scanmarch` dissolves into "Jacobi schedule × scanmarch kernel" and
  **ScanMarch gains Gauss-Seidel for free** (the inter-group reflect is kernel-agnostic);
  the G-S resolvent keeps the window default.
  (3) the matvec octant iteration moves from `quad.octants` onto `SweepSchedule.jacobi`
  (single-sources the in-plane projection `_octant_sweep`; bit-safe — the matvec has no
  cross-ordinate reductions).
  (4) **A2D-1 RETIRES** at sub-step (a) in favour of the `window≡full` MATVEC `array_equal`
  oracle (gate memo §2; a source-hash on a shared frame trips on every legitimate refactor
  with no behavior signal).
  (5) TWO NEW SUB-STEPS appended to the staged (a)–(d):
  **(e) collapse the graph's 4 level-walks → 2 direction-parameterized walks** (full,
  windowed), reusing the (a) emit-policy objects — same anchors (window≡full both directions,
  G2.c, affine golden); the `CellUpdate` override point for future Step/LD becomes the pure
  kernel pair (they supply cell algebra; storage handled once above);
  **(f) module geography — `sweep.py` DISSOLVES**: `_x_scan_faces`/`_scanmarch_row` →
  `spatial/scan.py` (the scan family's "graph+walk together"); the octant walk + schedule
  driver + `_sweep_1d_unified` + `transport_sweep` → `loss_representation.py`;
  DAG construction cache = `SweepDependencyGraph.for_shape(shape)` (construction with the
  graph) accessed via the `_DAGWavefront` family accessor (ownership with the family; sub-step
  (c) unchanged, better-homed); `solver.py` re-points (`_reflect_outflow_into_inflow` already
  native there — a boundary concern); consumers = 3 production files + ~30 test files
  (retirement = test migration). The lazy-import cycle DISSOLVES — no workaround left.
  ⭐ (f)-SCOPE ADDITIONS from the (d) execution (2026-06-11): (i) RENAME
  `_sweep_2d_scheduled` → `_sweep_scheduled` when it relocates (its buffer setup went
  d-generic at (d) — the name lies; one rename instead of two); (ii) DECIDE the fate of
  `WavefrontFlux` + `InteriorFaceSpace` — the (d) fold left the typed cochain with NO
  production consumer (its ι_*/ι* boundary algebra dissolved into the shared `_OctantWalk`
  frame; the fuller view survives as the full-cochain kernel's raw per-axis buffers); its own
  test file remains. Retire (aggressive-retirement: the structural-reference role transferred)
  or keep as a typed teaching/MoC-adjacent concept — USER call at (f).
  Each of (e)/(f) is pure relocation, independently bit-identity-gated; (f) additionally gated
  by clean collection + the full anchor set. Gate-memo addendum for (e)/(f): test-architect
  (dispatched at (a) start).
- **S6.5 — Unify the two doors.** `InvertibleOperator.solve` + `StreamingOperator.apply` reach
  ONE `SpatialRepresentation` instance (the `InvertibleOperator` holds/derives the same instance
  as its `streaming` operand). `solve` calls `representation.sweep`; `apply` calls
  `representation.loss_action`. The "assert same instance" test goes green → the L21 invariant
  becomes a TYPE FACT, not a test coincidence. ⚠ Keep the asymmetry: `loss_action → (L+C)ψ`
  + operator `−C`, but `sweep → (L+C)⁻¹q` DIRECTLY (the inverse is of the full composite — no
  `+C` to add back; `apply` is on the leaf `L`, `solve` is on the composite).
- **S6.6 — (FORWARD / deferred) the explicit-matrix representation.** `ExplicitMatrix`
  representation: assemble `(L+C)` as a sparse lower-triangular CSR (block-diagonal per ordinate;
  the DD stencil → sparse) in `__post_init__`; `sweep = spsolve_triangular`; `loss_action = M@ψ`.
  Reuse CP's `numerics.operator` assembly idiom (`cp/solver.py` already assembles+inverts at the
  SAME `numerics.operator.LinearOperator` layer — `operator.py:18-23` names both as peer
  representations). VALUE: a second, *construction-exact* solve/apply oracle (`spsolve_triangular`
  IS the inverse of `M@ψ`), `scipy.sparse.tril(M)` makes triangularity runnable, and it unlocks
  direct/ILU preconditioners + a sparsity view for DSA. Build when a consumer wants it; the
  renamed Protocol already admits it.

## S6.7 — FORWARD DESIGN / CROSS-METHOD (pollination + scope fences)

- **CP (collision probability):** already lives at the explicit-matrix end (assembles dense
  matrices + `np.linalg.solve`, `cp/solver.py:547`). SN's `ExplicitMatrix` is the SAME
  representational pattern at the SAME `numerics.operator` layer — build it with CP's assembly
  primitives ⇒ ~free.
- **Diffusion:** wraps as `scipy.sparse.linalg.LinearOperator` + `bicgstab`
  (`diffusion/solver.py:192`). SN's representation should round-trip through
  `numerics.operator.as_scipy_linop` (already exists) so the SN Krylov inner loop + diffusion's
  BiCGSTAB share one adapter. **S6 gate:** confirm `as_scipy_linop(representation)` still drives
  the SN Krylov inner after the rename.
- **MoC — SCOPE FENCE (do NOT cross).** MoC uses fiber bundles + characteristic rays, NOT a
  cell-DAG ([[project_moc_structure]]). Keep `SpatialRepresentation` the **SN loss operator's**
  representation; do NOT generalize it to a transport-wide concept. The `CellVisit`/`dag_walk`
  abstraction is SN-only.

## S6.8 — VERIFICATION (the gates)

- **Bit-identity inheritance:** A2D-1 source-hash regenerates (output byte-identical via the
  `window≡full` oracle); `test_affine_carve_bit_identity`, the DD regression snapshots,
  `SI≡Krylov≡k_inf` all green under `python -O`.
- **The discriminating test (the structural payoff):** `StreamingOperator.sweep_strategy IS`
  the instance `InvertibleOperator.solve` uses — proves "one operator, two actions, one
  representation" by construction (Smell #16 closed).
- **`ScanMarch.loss_action ≡ FullFieldWavefront.loss_action`** (the renamed G2.c) stays green
  through the move.
- **(S6.6) the explicit-matrix oracle:** `spsolve_triangular(M, q) == CumprodScan.sweep(q)` and
  `M @ ψ == CumprodScan.loss_action(ψ)` to nulp on a small slab; `scipy.sparse.tril(M)` is
  exactly zero above the visit-order diagonal (triangularity made runnable).
- **elegance-enforcer PASS:** zero `_apply_*` survive on the operator; one octant frame; the
  two IOU twins retired (retirement audit).

## S6.9 — RELATIONSHIP TO S5.3 / sequencing

S6 **is** the consolidation the S5.3 IOUs point at, pulled forward (the cross-domain-attacker:
"this carve IS the S5.3 trigger; do it now"). After S6 makes the architecture right (operator =
algebra; representation = walk; one instance; `_OctantWalk2D` frame), the S5.3 **measurement**
(does `ScanMarch` beat `MovingFrontierWindow` end-to-end?) decides whether the window
*representation* retires (Fork B2) — operating now on CLEAN representations. The re-layering is
valid regardless of the flip outcome (it makes both window + scan-march clean representations
sharing `_OctantWalk2D`; retirement, if it comes, just deletes one `_interior_walk`).

## S6 — RECOVERY (read order at pickup)

1. This S6 section (the locked design).
2. The cross-domain-attacker memo (agentId `a8826036ef5b5e945`) — the 4-frame structural backing,
   the code-grounded confirmations, the pitfalls.
3. The S5.1 commits `8913229` (sweep) + `0613298` (matvec) — the `ScanMarch` strategy + the
   `_x_scan_faces` Pattern-2 primitive + the two Fork-B1 IOU notes S6.4 retires.
4. `diamond.py:367-376` (the gather/scatter-shared, kernel-forked proof) + `operator.py:1641`
   (`_apply_2d_cartesian`) + `operator.py:1791` (`_apply_2d_cartesian_scanmarch`) + the renamed
   `sweep_strategy.py`.
5. FIRST ACTION = dispatch `test-architect` (S6.0-prime) for the bit-identity + "one instance"
   gate plan; do NOT move code before it lands.
