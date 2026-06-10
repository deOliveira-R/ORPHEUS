# N-D layout foundation — `(N, ng, *spatial)` array convention (eliminate the phantom `ny=1`)

**Branch/worktree:** `worktree-sn-nd-layout` (fresh off `main` @ `f0ffb30`).
**Mode:** main-agent direct authorship, turn-by-turn user steering
(`feedback_no_method_implementer_for_surgical_carves`). Commit per stage; gates
throughout. Standing commit exclusions: `.claude/agent-memory/`, `.mcp.json`,
hooks, `derivations/diagnostics/`.

---

## ⭐⭐⭐ POST-COMPACTION RECOVERY — READ THIS FIRST (2026-06-09)

**One line:** the production `(N, ng, *spatial)` rank-d carve **AND** the full
per-file test migration (C2c) **AND** the C2d integrated gate are **DONE** — full
suite **2345 passed / 0 genuine failures**, branch pushed to `origin`. This session
also CLOSED #214 (degenerate-`ny=1`, via an SO(2) product quadrature) and ran a
Lebedev-quadrature audit (~5 min suite speedup). **C3 IN PROGRESS — see the C3
SESSION UPDATE immediately below (architecture PIVOT: dimension-agnostic wavefront).**

### ⭐ SESSION UPDATE 2026-06-09 (C3 START — ARCHITECTURE PIVOT: dimension-agnostic wavefront)
**USER DIRECTIVE (re-shapes C3):** *"The wavefront should work for 1, 2 and 3D,
dimension agnostic. The cumprod should be an optimized option, and we can check that
wavefront matches cumprod for 1D (and even measure the improvement)."* This OVERRIDES
the `nd_foundation.md §2.3` deferral ("1-D is a parallel-prefix SCAN not a wavefront →
WRONG-FIT, defer the fold"). The hard constraint it flagged (a naive d-generic
forward-sub walk is O(nx) SEQUENTIAL for d=1 → perf regression) is RESOLVED: cumprod
stays the SELECTABLE d=1 optimization (1-D production keeps Blelloch speed); the
wavefront becomes the general mechanism + the 1-D verification oracle.

**ARCHITECTURE (spine + pluggable optimizations — collapses 3 sweep + 2 matvec bodies
into ONE spine + oracle-pinned optimizations):**
- WAVEFRONT DAG full-field walk (`SweepDependencyGraph.from_cartesian(shape)` +
  d-generic `DiamondDifference.cell_kernel_batch` summing `denom=σ_t+Σ_axis s_axis`) =
  the dimension-agnostic SPINE + verification ORACLE (d=1,2,3, ONE code path).
- cumprod parallel-prefix scan (`ordinate_scan`/`_sweep_1d_unified`) = the d=1
  OPTIMIZATION (default-selected; ZERO perf regression). Pinned wavefront(1D)≡cumprod
  = principled-equiv @ nulp (DIFFERENT FP association orders, anchored on analytical
  k_inf=1.875 NOT old-vs-new).
- `_MovingFrontier` window = the d≥2 OPTIMIZATION (pinned ≡ full-field by
  `tests/sn/sweep/cartesian_2d/test_2d_full_field_oracle.py`).

**C3 SUB-STAGING:** C3.0 scaffold (test-first) · C3.1 d-generic DAG spine · C3.2
d-generic diamond kernel + d-generic full-field walk · C3.3 geometry d-generic
(`itertools.product` + build the 1-D graph + fold the 5 `ny>1` gates + `is_1d`) ·
C3.4 wavefront-1-D wiring + cumprod oracle + perf measurement · C3.5 orchestration
d-generic (`xmin/...`→per-axis `FaceLabel`) + retire the sweep/matvec twins · C3.6
synthetic-3-D admission + docs (folds C5).

**DONE THIS SESSION:**
- ✅ **C3.0** (`5b941eb`): test-first scaffold — `test_wavefront_cumprod_equivalence.py`
  (wavefront≡cumprod oracle + cumprod→analytical-k_inf anchor [live] + perf measurement)
  + `test_sweep_graph_nd_admission.py` (B1–B6 synthetic d=1/d=3 structural pins + the
  d=2 `from_cartesian`==`from_cartesian_2d` legacy-equivalence). xfail-gated green.
- ✅ **C3.1** (`9b75374`): `OctantLabel(signs:tuple)` d-generic (compat props
  `sign_x`/`sign_y`/`streams_in_2d`/`ndim`/`streams`) + `SweepDependencyGraph.from_cartesian(shape,*,label)`
  d-generic (levels = the `Σ idx=k` anti-hyperplane, `Σ(n−1)+1` levels; C-MAJOR within
  a level ⟹ d=2 BIT-IDENTICAL to legacy `from_cartesian_2d`, now a thin deprecated
  alias). Fields → `face_in`/`face_out`/`spatial_shape` tuples (+ compat props
  `face_in_x`/`nx`/…). `_MovingFrontier` window = d=2-ONLY optimization
  (`window_slots`/`window_edges` None otherwise; `apply_windowed`/`residual_windowed`
  fail-loud guard). 6 production + 4 test-file `OctantLabel(a,b)`→`OctantLabel((a,b))`;
  2 validation `match=` strings updated. ⭐ GATES ALL GREEN: smoke (d=2 bit-id vs legacy
  incl. `ny=1` degenerate + d=1 chain + d=3 admission); scaffold d-gen pins flipped
  xfail→pass (49 passed/3 xfailed); 2-D moat bit-identical; sweep/core+cartesian_2d+
  primitives+slab+dd_regression **754 passed -O** (DD drift **6920 ULP == pre-C3.1
  baseline → ZERO new drift**); `test_sweep_graph` 76 passed no-`-O`; k_inf=1.875.

- ✅ **C3.2a** (`9063487`): `DiamondDifference.cell_kernel_batch`/`residual_kernel_batch`
  dimension-generic via a per-axis TUPLE API (`psi_in: tuple`, `s: tuple` → `(psi_avg,
  psi_out: tuple)`). The axis reduction is an EXPLICIT LEFT FOLD (NOT `sum()`) so the
  `((sigt+s_0)+s_1)` / `((Q+s_0·in_0)+s_1·in_1)` order reproduces legacy `sigt+sx+sy` /
  `Q+sx·in_x+sy·in_y` BIT-FOR-BIT at d=2 (IEEE-754 non-associativity is load-bearing).
  4 callers (`update_batch`/`residual_batch` + `apply_windowed`/`residual_windowed`) pass
  the (x,y) pair as d-tuples + unpack. B5(d=3)+B6c(d=1) admission pins WIRED to the REAL
  kernel (was hand-oracle). Retired 2 dead imports. GATES: scaffold+window-equiv+2-D
  octant+full-field oracle+cell_update_batch+dd_regression **112 passed -O**; DD drift
  **6920 ULP == baseline → ZERO new drift**.

- ✅ **C3.2b** (d-generic full-field WALK — so d=1/d=3 CAN walk): `SweepCellSlice`
  (`cell_update.py`) 2-axis fields → per-axis TUPLES (`cell_idx`/`face_in_idx`/
  `face_out_idx`/`psi_faces`/`str_axes`); `DiamondDifference.update_batch`/`residual_batch`
  gather/scatter → d-generic via shared helpers `_cell_face_selector` (the
  `(:, :, *cell_idx)` advanced index with axis-a's cell idx replaced by its face idx) +
  `_gather_cell_inputs` + `_scatter_outgoing_faces` — a Pattern-2 extraction collapsing the
  REAL twin the two batch methods carried (gather/scatter identical, only the cell math
  differs); `SweepDependencyGraph.apply`/`residual` walk `for ii,jj in levels` →
  `for cell_idx in levels` + d-generic `_make_slice` + `buf[(slice(None),…,*cell_idx)]`
  scatter; signatures `psi_x/y_octant`→`psi_faces_octant: tuple`, `str_x/y_octant`→
  `str_axes_octant: tuple`. The KERNEL (`cell_kernel_batch`) was already d-generic (C3.2a),
  UNTOUCHED → bit-identity left-fold intact. ⭐ **WavefrontFlux-typing acceptance item
  RESOLVED (elegance-enforcer PASS):** the buffers handed to the walk are OCTANT-RESTRICTED
  working copies (`face(a)[oct_idx].copy()`, `(N_oct,…)` not the full mesh `(N,…)`) — they
  CANNOT be a mesh-bound `WavefrontFlux` without violating its `_check_partner` binding (a
  type-lie in the OTHER direction; Pattern 4 reversed). So the typed axis↔buffer map lands
  at the ORCHESTRATOR (`tuple(wavefront.face(a) for a in wavefront.axes)` — the typed
  `WavefrontFlux.axes` map, replacing hardcoded `psi_x=face(0); psi_y=face(1)`); the slice +
  kernel carry the raw octant-projected tuple in that axis order (Pattern 7 — normalise at
  the definition site). NO fictitious octant-submesh (single-consumer; defer per Pattern 6).
  ⭐ NEW **B7 walk-admission test** (`test_sweep_graph_nd_admission.py`): the d-generic WALK
  (not just the kernel) verified at d=1 AND d=3 via the geometry-agnostic matvec==sweep
  round-trip (`residual(apply(x))≈0`, 10 passed) — makes "d=1/d=3 CAN walk" a TESTED fact.
  Migrated `test_cell_update_batch.py`+`test_sweep_graph.py`+`test_sweep_graph_window_equivalence.py`
  (the SweepCellSlice/apply/residual external readers; windowed `str_x/y_octant` kwargs KEPT —
  the window path stays d=2). ⭐ **A2D-1 source-hash pin REFRESHED** (`2d0c9f…`→`b17d69…`): it
  had been STALE since C3.1 (`9b75374`'s `OctantLabel((sx,sy))` migration changed
  `_apply_2d_cartesian`'s text), uncaught for 3 commits because C3.1/C3.2a gates skipped the
  operators suite — caught + fixed now (the change is value-preserving; matvec full-field
  oracle equivalence confirms production correct). GATES ALL GREEN: sweep core+cartesian_2d
  **475 passed** (incl. window≡full bit-id oracle + end-to-end full-field oracle + B7);
  operators+solve **524 passed** (incl. sha256 bit-identity anchor `test_affine_carve_bit_identity`
  = sub-ULP proof of ZERO byte change in the 2-D path + A2D-1); eigenvalue **53 passed** (2-D
  k_eff); primitives 271 / spatial 4 / regression 13. ⚠ The 8 verification "failures" are
  pytest-TIMEOUTS (>60s) on PRE-EXISTING-slow paths C3.2b NEVER touches (1-D slab cumprod +
  cylinder curvilinear — diff + import audit prove independence; slab confirmed passes at
  generous timeout; cylinder `cyl_2g_3reg…n40` is heavy/non-converging #206-adjacent). DD
  drift stays **6920 ULP == baseline → ZERO new drift**. elegance-enforcer **PASS-WITH-NITS**
  (no violations; committable). Commit chain: `5b941eb`(C3.0)→`9b75374`(C3.1)→`9063487`(C3.2a)
  →`10e2587`(C3.2a elegance)→`f88ef7d`(plan)→**C3.2b** (this commit).

**⏭ NEXT = C3.3 (geometry d-generic) — file:line targets:**
1. **Octant build → `itertools.product`.** `geometry.py:1311-1313`: the dict comp
   `{OctantLabel((sx,sy)): SweepDependencyGraph.from_cartesian_2d(nx=…, ny=…, label=…)}`
   (the `(sx,sy)` come from a hardcoded 2-D octant list above 1311) → `for signs in
   itertools.product((-1,+1), repeat=ndim): OctantLabel(signs): from_cartesian(spatial_shape,
   label=…)`. Retires `from_cartesian_2d` (the C3.1 alias) + its legacy-equivalence pin
   `test_d2_from_cartesian_matches_legacy` becomes a frozen golden (see its NOTE).
2. **Build the 1-D sweep graph (None today).** `geometry.py:342,349`: `self.sweep_graphs =
   None` for the non-2-D mesh → build the d=1 graph dict (`product((-1,+1), repeat=1)` = 2
   octants on `from_cartesian((nx,))`). This is what C3.4's `_wavefront_1d_sweep` consumes.
3. **Fold the 5 `ny>1` dispatch gates → a dimensionality test.** `operator.py:373,618,831`
   (`if curvature == "cartesian" and ny > 1:`) + `operator.py:1460,1530` (`if curv is None
   and sn_mesh.ny > 1:`) select the 2-D-Cartesian wavefront path vs 1-D/curvilinear. Replace
   the `ny > 1` phantom-axis test with `ndim == 2` / `not is_1d` (genuine dimensionality).
4. **Genericize `is_1d`.** `geometry.py:628-630`: `return self.ny == 1` → `return self.ndim
   == 1` (or `len(spatial_shape) == 1`). Add a dimensionality property if missing.
   Production octant read sites STAY (`operator.py:1697,1812`, `sweep.py:902,1250` —
   `sweep_graphs[OctantLabel((sx_eff,sy_eff))]` lookup by label).
THEN C3.4 (wire d=1 wavefront via the d-generic walk in `_wavefront_1d_sweep` — see
`test_wavefront_cumprod_equivalence.py` `_wavefront_1d_sweep` docstring sketch + the
`_SpineNotLanded` raise to delete — + flip the cumprod oracle xfail→pass + MEASURE the
speedup; window_slots/window_edges → `SweepOptimization` sum type now it's the 2nd
optimization) + C3.5 (orchestration d-generic — incl. the `str_axes` → `axes`-keyed
`sn_mesh.streaming(a)` map fix from the elegance CONCERN — + retire twins + the
`OctantLabel.sign_x/y`/`streams_in_2d` shims) + C3.6 (3-D admission + docs).

⭐ **ELEGANCE-REVIEW ACCEPTANCE ITEMS (carry forward):**
- ✅ **C3.2b WavefrontFlux-typing — RESOLVED** at the orchestrator boundary (see C3.2b above;
  elegance-enforcer confirmed the orchestrator-not-slice landing is the CORRECT call, not a
  shortcut — the octant transient cannot validly be a mesh-bound cochain).
- **C3.4/C3.5 (elegance-enforcer CONCERN, from C3.2b):** `str_axes = (streaming_x, streaming_y)`
  at both orchestrators is a SECOND hand-listed axis-order tuple NOT derived from
  `wavefront.axes` (latent-only at d=2 — both tuples length-2, order-trivial; the habitat opens
  at the d≥3 orchestrator). When d≥3 orchestration lands, expose an `axes`-keyed
  `sn_mesh.streaming(a)` / `streaming_axes` map so `str_axes = tuple(sn_mesh.streaming(a) for a
  in wavefront.axes)` — ONE axis-order source for BOTH per-axis tuples; acceptance criterion
  "no orchestrator hand-lists a positional per-axis tuple". A call-site contract comment pins
  it for now (`sweep.py`/`operator.py` orchestrators).
- **C3.4:** revisit `window_slots`/`window_edges` (`tuple|None`) as a
  `SweepOptimization|None` sum type ONLY once the d=1 cumprod accelerator makes it the 2nd
  optimization instance (`feedback_unify_after_two_instances`) — NOT before.
- **ACCEPTED as-is:** batch kernel ≠ Pattern-2 twin of scalar `cell_balance_for_streaming`
  (Cartesian multi-axis streaming vs curvilinear single-axis + Morel–Montry angular);
  explicit-left-fold-over-`sum()` exemplary (bit-identity why-at-the-what); the partial
  `OctantLabel.sign_x/y` / `streams_in_2d` shims retire in C3.5; the `(slice(None), *cell_idx)`
  idiom stays inline (2 distinct leading-axis shapes — a helper would need a flag).

Details below are the PRE-pivot C3 plan (superseded by the staging above).

### ⭐ SESSION UPDATE 2026-06-09 (cont., post-compaction) — production gaps found+fixed; test migration ~⅔ done
**The test migration surfaced 3 PRODUCTION rank-2 remnants the k_inf/DD
validation never exercised (so the earlier "NO production bugs" was
incomplete).** ALL FIXED + re-validated (commit `94fcae5`):
1. `solution.py` `reaction_rate_density` `einsum("gxy,gxy->gxy")` → elementwise
   `xs * scalar_flux.values` (σ⊙φ, rank-generic). Never on the k_inf path.
2. `solver.py` `__debug__` cell-flatten invariant `total_cross_section.transpose(1,2,0)`
   → `np.moveaxis(...,0,-1)` + `reshape(*spatial_shape, ng)`. Fired only WITHOUT `-O`
   (default mode is `-O` → the gate was silently inert post-carve).
3. **HarmonicMomentField was a rank-2 ISLAND** (its `_phase_space_shape` never
   flipped). Flipped HMF + the `MomentDisplacement` mirror to `(L+1,2L+1,ng,*spatial)`
   — 5 shape-tuple sites, NO einsum cascade (moment einsums already ellipsis-generic;
   2-D wavefront writes stay legit 2-D). HMF now rank-d; the 2-D-only GATING of the
   windowed-moment path (Morel–Montry seed needs per-ordinate ψ) is untouched.
   Plus HARMLESS `n_dof=N*ng*nx*ny`→`N*ng*prod(spatial_shape)` + retired dead
   `_scalar_flux_from_angular`. **An explorer swept ALL of `orpheus/`: NO remaining
   LIVE rank-2 bugs**; residual `ny>1` dispatch gates (`operator.py` `_apply_2d_cartesian`
   selector) + `geometry.py:629 is_1d` are HARMLESS today → genericize in **C3**;
   stale `(ng,nx,ny)` docstrings → later sweep.
RE-VALIDATED: k_inf=1.875 slab/sphere/cyl (scalar_flux now genuinely rank-1
`(ng,nx)`); 2-D octant moat 7/7; DD regression 13/13; 9 bc_extraction `.npy`
baselines re-captured + INDEPENDENTLY verified value-preserving (`new==old[...,0]`).

**TEST MIGRATION STATUS: C2c COMPLETE — ALL regions migrated + per-region green**
(modulo pre-existing reds #206/#195/#212/#214). ~38 files across 4 waves of
parallel general-purpose sub-agents (brief = the validated method; reference diff
`git show 1774586`). Per-region gates: operators+primitives+solve+numerics+transport
**1650 passed/0 failed**; sweep core+slab **397**; cartesian_2d **38** (genuine-2-D,
LEFT UNTOUCHED); curvilinear/verification/eigenvalue green. Genuine-2-D files
(`cartesian_2d/*`, `test_keff_2d`, `test_mms_2d`, `test_sweep_graph`, `test_cell_update_batch`)
gate-confirmed green + untouched (rank-d-correct already). #206 cylinder hand-ref now
CORRECTLY xfails (migration made it COMPUTE — was dying at IndexError). #195 stays an
expected xfail(strict). 9 bc_extraction `.npy` + phase_c_crosscheck `.npz` snapshots
load value-preserving (verified `new==old[...,0]`).
**Commit chain (this session):** `1774586`(ref batch) → `94fcae5`(PRODUCTION carve-completion:
3 gaps + HMF island) → `9870495`(wave-1) → `598a256`(plan) → `af38cfb`(wave-2/3) →
`ce6442f`(wave-4). **C2d DONE — full-suite integrated gate GREEN: 130 files, 2345 passed, 0 genuine
failures, 13 xfailed (#206/#195), 7 skipped (all legit).** ⭐ **#214 CLOSED** (the
degenerate-`ny=1` 2-D-wavefront crash — the ONLY test broken *specifically* by
dimensionality handling; a full skip/xfail audit confirms no others). The rank-d
carve incidentally fixed it (`ny=1` Mesh2D eigenvalue → k_inf, no more
"no 'y' interior face-field"); `test_sweep_1d_2d_consistency` un-skipped + green in
1.4s, with `Lebedev(17)`→`product(n_mu=8,n_phi=4)` chosen by the quadrature objects'
`invariance_group`/`degree_of_exactness` tags (#154): SO(2) matches the slab-in-x
degenerate-`ny=1` symmetry + the 1-D GL reference; Lebedev was an O_h S² SO(3) cubature,
500× slower. Commit `Closes #214` + issue comment. **NEXT: C3 (sweep DAG d-generic;
FOLD IN the deferred `ny>1` dispatch-gate genericization in `operator.py`
`_apply_2d_cartesian` selector + `geometry.py:629 is_1d` → dimensionality test, per the
explorer audit) → C4 (#220 — boundary dict, still OPEN/forward) → C5 (3-D admission
pins + the `(ng,nx,ny)` docstring sweep). Other dim-adjacent OPEN issues are FORWARD
work, NOT closeable: #220 (C4), #219 (MethodSpace arch); #210 (2-D MMS pin) is deferred
for PERFORMANCE not dimensionality.**

⭐ **LEBEDEV QUADRATURE AUDIT (this session, post-#214) — committed `e9b123a`+`0b8cbba`.**
Prompted by #214's fix: since ORPHEUS has NO 3-D geometry, Lebedev (an `O_h` `S²`
cubature for SO(3) spherical-harmonic MOMENT integration, N=110 @ order 17) is the
WRONG tool for a plain 2-D transport sweep. Audited ALL ~140 Lebedev usages (26 files).
**Choice framework via the objects' OWN tags** (`q.measure.invariance_group`,
`.degree_of_exactness`, `q.N`): genuine 2-D Cartesian (`O_h`) → `level_symmetric(4)`
(same `O_h`, N=24); degenerate/axial (SO2) → `product` (the #214 case); Lebedev ONLY
for high-L moment integration or testing quadrature-independence. P1 anisotropy needs
just the degree-2 moment → `level_symmetric(4)` (doe=3) integrates it exactly, so Lebedev
is overkill even for P1. **Sped up the 4 slow ones** (the only ones worth it; the rest are
<2s — left): `test_keff_2d` 322→90s (was timing out the 200s cap), `test_discrete_ordinates_2d`
115→35s, `test_fixed_source_2d_equivalence` 6.4→1.5s, `test_solver_components` one solve.
≈5 min suite speedup. Every swap justified by a discriminator (flux-shape-independent
k_inf / convergence-trend / two-sided equivalence — swap BOTH legs); NO assertion/tolerance
weakened; green under -O AND without. **KEPT Lebedev (load-bearing):** `test_gl_and_lebedev_agree`
(Lebedev IS the subject), `test_z_ordinates_contribute` (`level_symmetric` has ZERO pure-z
ordinates → a swap would make it VACUOUS — the sharp catch), the spherical-harmonic
orthogonality/addition-theorem + scattering moment tests. Legit Lebedev (cubature-property
+ moment-integration tests) untouched.

### Environment (Host env, `$CLAUDE_ENVIRONMENT` empty)
- **Worktree (cwd for everything):** `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/sn-nd-layout`
- **Branch:** `worktree-sn-nd-layout` (off `main` @ `f0ffb30`; `main` is the
  up-to-date trunk — Wave O + Phases 1–5 + #208 all merged; the old
  `field-role-typing` worktree refs in other memories are STALE). NOT pushed.
- **venv:** `/Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python`
- **Test invocation:** `PYTHONPATH=/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/sn-nd-layout /Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python -O -m pytest -p no:cacheprovider ...`
  Default mode `-O`. **Sentinels run WITHOUT `-O`** (`pytest -m sentinel`, Mode 8).
  **Deselect the #212 hang:** `--deselect "tests/sn/eigenvalue/test_keff_slab.py::test_heterogeneous_absolute_keff"`.
  Wrap slow runs in `gtimeout 600` and redirect to a file (full `tests/sn/` times
  out ~25 min; gate in slices). NO `continuous_get` (#212).

### Commit chain (all on the branch; clean tree at HEAD)
`646960d` plan → `05545ba` C1 (rank-adaptive local-op broadcasts, bit-identical)
→ `ee04870` C2 lever+builders (WIP) → `d6e401b` C2b structured-1-D matvec
(swapaxes(0,1), `(ng,N,nx)`, pole closure) → `6ae3da8` rank-d MMS sources +
`_build_fixed_source_rhs` validation → `79f49e6` test literal-`1` migration +
DD-snapshot regen → `6ca0154` `legacy_proxy_matvec` helper → `074372a`+`9fff2ee`
plan/recalibration.

### VALIDATED (production correctness — re-run any to confirm post-compaction)
- slab/sphere/cyl `solve_sn(get_mixture("A","2g"), Mesh1D, …)` → **k_inf=1.875**
  (closed-form anchor; smoke script pattern is in this plan's history).
- **2-D bit-identical**: `tests/sn/sweep/cartesian_2d/test_2d_octant_sweep_equivalence.py`
  green; `_apply_2d_cartesian` UNTOUCHED.
- **DD regression 13/13**: `tests/sn/regression/test_dd_regression.py` (10 1-D
  snapshots regenerated + VERIFIED value-preserving: 6 squeeze-only, 4 FP-drift
  1e-16..6.8e-12 within tol; 2-D untouched).
- MMS 1-D path runs (sources are rank-d). **No production bugs** — every
  remaining suite failure is TEST-side phantom.

### THE REMAINING WORK — do these IN ORDER
**(1) C2c test migration — ✅ DONE.** ~38 files migrated across 4 sub-agent waves;
production carve-completed (3 late-found rank-2 gaps + HMF island, `94fcae5`). See the
SESSION UPDATE block above for the full record + method (still valid if more phantom
surfaces: build via the in-scope mesh's `*mesh.spatial_shape`; drop phantom extraction
indices; gate per file; NEVER blanket-script — the 3 reverted approaches are recorded).
**(2) C2d full-suite gate — ✅ DONE.** 130 files, **2345 passed, 0 genuine failures**,
13 xfailed (#206/#195), 7 skipped (legit). Per-file loop (each file isolated, 200s cap;
deselect #212) is the robust gating recipe — a monolithic `pytest tests/sn/` buffers
output + times out. ⭐ ALSO this session: #214 CLOSED + a Lebedev-quadrature audit
(see SESSION UPDATE above).
**(3) C3 sweep DAG d-generic (NEXT) ⭐** (`OctantLabel(signs:tuple)`, `from_cartesian_2d`→
   `from_cartesian(dims)`, `diamond.cell_kernel_batch` `denom=σt+Σ_a s[a]`,
   `geometry` octant build `itertools.product((-1,+1), repeat=ndim)`; 2-D
   bit-identical sum-in-x,y-order; synthetic-3-D shape admission). See `nd_foundation.md` §2.
   **FOLD IN** the explorer-audited HARMLESS `ny>1` dispatch-gate genericization:
   `operator.py` `_apply_2d_cartesian` selectors (`sn_mesh.ny > 1` ×5) +
   `geometry.py:629 is_1d` (`self.ny == 1`) → a dimensionality test (`reduced is None` /
   `len(spatial_shape)`). These are correct-today via ny=1 but carry the phantom.
**(4) C4 boundary inventory #220** (`_resolve_bcs`→`dict[FaceLabel,op]`; still OPEN/forward,
   nd_foundation §2.6 needs it for 3-D).
**(5) C5 3-D shape-admission pins + docs** (extend `docs/theory/index_convention.rst`
   `(N,ng,nx,ny)`→`(N,ng,*spatial)`; sweep the stale `(ng,nx,ny)` docstrings the
   explorer flagged across `solver.py`/`operator.py`/`material_xs_field.py`/etc.).

### Task tracker: #1 C1 ✅ · #2 C2a ✅ · #3 C2b ✅ · #4 C2c ✅ · #5 C2d ✅ ·
**#6 C3 (NEXT)** · #7 C4 (#220) · #8 C5. Chain blockedBy 6→7→8.

---

## STATUS (live — 2026-06-09)

### ⭐⭐ UPDATE 2026-06-09 (latest) — TEST MIGRATION IS A DEEP PER-FILE REWRITE (recalibrated)
The remaining test migration is NOT a scriptable "mechanical tail" — it is a
**deep per-file rewrite** of ~47 files (esp. the operator-algebra suites:
`test_streaming_operator` 35 fails, `test_bc_extraction_matvec` 30,
`test_streaming_operator_decomposition` 19, `test_native_matvec` 16,
`test_invertible_operator`/`test_g_adjoint_reciprocity`/`test_collision_operator`
/`test_operators_apply_typed` ~12 each). Each thoroughly assumes rank-2 spatial
in CONSTRUCTIONS + EXTRACTIONS (`[..., 0]`, `[:, :, :, 0]`) + slices + assertions.
**Three scripting approaches were tried and REVERTED** (all clean now):
1. blanket `(N,ng,nx,ny)`→`(N,ng,nx)` — REGRESSED hidden-2-D cases (e.g.
   `test_streaming_operator_decomposition` 5×5 CART). Mixed files unsafe.
2. `(N,ng,nx,ny)`→`(N,ng,*<meshvar>.spatial_shape)` — 87 NameErrors: the mesh var
   VARIES per construction (helper functions take it as a param with other names;
   `sn_mesh`/`sn`/`op`/`local_sn_mesh`). Not a single-var-per-file problem.
3. rank-conditional `(N,ng,nx) if ny==1 else (N,ng,nx,ny)` — robust (no meshvar,
   no NameError) but only +7 passes for heavy `if ny==1` clutter; doesn't touch
   the extraction class. Low gain / high clutter.
**Correct path = methodical PER-FILE (per-function) manual migration**: in each
test/helper, use the in-scope mesh's `spatial_shape` (or build the field via
`AngularFlux.zeros_on(mesh)` / `.from_mesh`) for constructions, and drop the
phantom index on extractions. Use the suite as the guide. NOT scriptable safely.
**Production carve is DONE + VALIDATED + COMMITTED — unaffected by this** (the
test suite is verification scaffolding; the branch is not merged). Commit chain
ends at `074372a`; clean working tree.

### ⭐ UPDATE 2026-06-09 (late) — PRODUCTION CARVE DONE+VALIDATED; test tail remains
**PRODUCTION `(N,ng,*spatial)` COMPLETE + VALIDATED + COMMITTED** (commit chain
`05545ba` C1 → `ee04870`+`d6e401b` lever+builders+structured-1-D → `6ae3da8`
MMS sources+fixed-source validation → `79f49e6` test shape-tuples/extractions +
snapshot regen → `6ca0154` legacy_proxy_matvec helper). VALIDATION: slab/sphere/
cyl `solve_sn` → k_inf=1.875; 2-D bit-identity moat green; DD regression 13
passed (10 1-D snapshots regenerated + VERIFIED value-preserving: 6 squeeze-only,
4 FP-drift 1e-16..6.8e-12 within tol). NO production bugs — every full-suite
failure is TEST-side.

**REMAINING = test migration tail (C2c, mechanical, NOT blanket-safe).** ~44 test
files build `(N,ng,nx,ny)` / `(ng,nx,ny)` with a local **`ny=1` VARIABLE** (a
literal-`1` regex can't see it; ~22 files are MIXED 1-D+2-D so `(N,ng,nx,ny)`→
`(N,ng,nx)` BREAKS the 2-D cases — confirmed: a blanket attempt regressed
`test_streaming_operator_decomposition`'s 5×5 CART cases, REVERTED). **Correct
universal fix = per-file `(N,ng,nx,ny)`→`(N,ng,*<meshvar>.spatial_shape)`**
(rank-correct for both 1-D and 2-D). Also remaining: `[..., 0]`/`[:, :, X, 0]`
extractions on now-rank-1 flux/source that the literal-`1` regex missed; some
`np.full((nx,1),…)` per-group stacks. Pattern catalog + the safe-vs-mixed file
split are in this session's history. The DONE-and-committed shape-tuple/extraction
migration (79f49e6) + helper (6ca0154) handled the literal-`1` subset.

**Task tracker** (recreate from this §STATUS + §Staging if it didn't survive a
session boundary): #1 C1 ✅ · #2 C2a builder layer ✅ · #3 C2b structured-1-D
surgery (IN PROGRESS) · #4 C2c 1-D test migration · #5 C2d C2 gate · #6 C3 sweep
DAG d-generic · #7 C4 boundary #220 · #8 C5 3-D admission pins + docs. Chain
blockedBy 3→4→5→6→7→8.


**C1 ✅ COMMITTED** (2 commits: `646960d` plan, `05545ba` rank-adaptive local-op
prep). 922 passed, A2D-1 re-blessed, sentinels green without `-O`. Bit-identical.

**C2 ⏳ IN PROGRESS (the lever + rank change) — branch is MID-CARVE (NOT green;
expected for an atomic rank change). UNCOMMITTED WIP in the worktree.**

DONE (the builder layer — all bit-identical on 2-D where reachable; verified the
1-D fail-loud advances past each):
- `_bases.py` `_shape_for_mesh` (Angular + Scalar) → `(N,ng,*spatial)` / `(ng,*spatial)` — THE LEVER, flipped.
- `source_sinks/angular_source_sink.py` `from_isotropic` (expected + broadcast → spatial_shape).
- `source_sinks/scalar_source_sink.py` `__add__` (160) + `as_per_ordinate` (198) `[None,:,:,:]`→`[None]`.
- `solver.py` `initial_flux_distribution` (886), `compute_group_production_rate` (943/955/958), `compute_group_absorption_rate` (982/985/988) `"g...,g...->g..."` + `moveaxis(0,-1).reshape(-1,ng)`, `n_dof` (1295), `_setup` rebind sig_t_1d (823).
- `material_xs_field.py` `_ensure_cell_views` (reshape→spatial), the 4 per-material verbs `for mid,idx` + `(slice(None),*idx)` (487/509/558/597) — bit-identical on 2-D (np.where already ndim).
- `geometry.py` `__init__` mat_map/volumes rank-d (283-284), `full_field_space` g_bulk rank-generic (824-836).

REMAINING (the structured-1-D surgery — careful per-FUNCTION reads; 2-D path
`_apply_2d_cartesian` is UNTOUCHED & stays bit-identical). Pattern: drop
`[:, :, 0]`/`[:, :, :, 0]` (sig_t/Q now rank-1); `transpose(1,0,2,3)`→`transpose(1,0,2)`
(and the `[...,0]` variant → `transpose(1,0,2)`); `out_g_first` zeros
`(ng,N,nx,ny)`→`(ng,N,nx)` + writes `[..., i, 0]`→`[..., i]`; pole seed
`psi_view[:, :, 0, 0]`→`[:, :, 0]`; drop `[..., None]` / `[:, :, None]` re-pad.
⭐ KEY SIMPLIFICATION: `transpose(1, 0, 2, 3)` (the N↔ng swap) → **`swapaxes(0, 1)`**
— rank-GENERIC (works rank-1 AND rank-2, no per-rank branch; reads as "swap the
ordinate and group axes"). Do the `transpose(1, 0, 2, 3)[..., 0]` → `swapaxes(0, 1)`
replace FIRST (the `[...,0]` phantom-drop vanishes at rank-1), then the bare
`transpose(1, 0, 2, 3)` → `swapaxes(0, 1)`. Same spirit: the rate `transpose(1,2,0)`
→ `moveaxis(0, -1)` (already applied in solver). These transpose sites are ALL
1-D-only (2-D uses `_apply_2d_cartesian`), so swapaxes is correct + 2-D-irrelevant.
SITES (re-confirm line #s):
- `operator.py` `_compute_LpC` (385-528: 389,390,392,393,404, the `[..., i, 0]` writes, 528),
  `_compute_decomposition` (625-740: 628,629,633, writes, 740 re-pad),
  `_compute_LpC_transpose` / 3rd variant (845-1007: 847,852,864,1006,1007), 1102.
- `sweep.py` `_run_1d_sweep`/`_sweep_1d_unified` (400,464,465 drops; 644 transpose; 773 `[:, :, None]` re-pad).
- `spatial/pole_angular_closure.py` 852 `transpose(1,0,2,3)[...,0]`→`transpose(1,0,2)` (1193 is a DOCSTRING — leave).
- `solver.py` 875 (rebind sig_t_1d, same as 823) + the `__debug__` cell-flattening check
  720/745/748 (`reshape(nx,ny,ng)`→`reshape(*spatial,ng)`, `transpose(1,2,0)`→`moveaxis(0,-1)`) — fires only WITHOUT `-O`.
- TEST MIGRATION: slab/sphere/cyl tests hardcoding `(N,ng,nx,1)` → `(N,ng,nx)` (e.g.
  `tests/sn/sweep/slab/test_unified_matvec_slab.py` `np.zeros((N,ng,nx,1))` / `np.full((ng,nx,1),…)`).

GATE LADDER for C2: (1) slab smoke `solve_sn(get_mixture("A","2g"), Mesh1D slab)` → k_inf
≈1.2659 (2g A); (2) sphere+cyl smoke; (3) full `tests/sn/` + 2-D bit-identity moat
(`test_2d_octant_sweep_equivalence.py` array_equal) + sentinels w/o `-O`; (4) 1-D
re-validation via k_inf + MMS (NOT old bytes). Then C3 (sweep DAG d-generic), C4
(#220 boundary), C5 (3-D admission pins + docs).

## Goal (user, 2026-06-09)
Make EVERY array genuinely `(N, ng, *spatial)`, spatial rank == ndim: `(nx,)`
1-D, `(nx,ny)` 2-D, `(nx,ny,nz)` future 3-D. **NOT adding 3-D compute** — a
dimension-agnostic FRAMEWORK that takes 3-D later. **Full rank-d, all paths incl.
curvilinear** ("how we get to a good architecture").

## The crux (empirically confirmed)
The pack convention (`SNMesh.spatial_shape`/`n_unknowns_flat`, from the ultraplan
C1 already on `main`) is genuine rank-d, but the typed FIELD arrays
(`_bases.py:189/255 _shape_for_mesh` → `(N,ng,nx,ny)`) are ALWAYS rank-2: 1-D
carries a phantom `ny=1` (`AngularFlux.values.shape==(4,1,5,1)` for a 5-cell
slab). **2-D is already `(N,ng,*spatial)`-correct → stays bit-identical.** The
phantom lives entirely on the 1-D path; curvilinear is also `ndim=1` → slab +
curvilinear move to rank-d TOGETHER (shared `_shape_for_mesh`).

## Operator-machinery audit (explorer, the partition the staging rests on)
**THE LEVER:** `_bases.py` `BulkField._shape_for_mesh`→`(N,ng,*spatial_shape)`,
`ScalarField._shape_for_mesh`→`(ng,*spatial_shape)`. `Field.__post_init__` enforces
`values.shape==space.shape` STRICTLY → flipping the lever makes every stale rank-2
producer **fail LOUD**. No load-bearing `ny=1` broadcast exists — the phantom is
vestigial positional indexing.

**SURVIVE UNTOUCHED (rank-agnostic already):** `FissionOperator`+`RankOneOperator`
(reduce axis=0, broadcast *spatial); all `MomentProjection`/`HarmonicMomentReconstruction`
einsums (ellipsis `n...`); `ordinate_scan`/`_pair_monoid_scan`; the
`MorelMontryAngularSweep` pole closure INTERNALS (already rank-1 — drops y on
entry: `psi_view.transpose(1,0,2,3)[...,0]`); all iteration drivers
(`SourceIteration`/`KrylovAcceleration`/`KEigenvalue`, via the ravellable
protocol); `TimedFullField.to_flat`/`from_flat` (ravel/reshape vs template);
the 1-D face layout (already `(N,ng)`); the ENTIRE 2-D Cartesian path.

**MUST CHANGE:**
- *Rank-ADAPTIVE-safe (works for current rank-2 AND future rank-1, bit-identical
  on both → prep commit C1):* `operator.py:1963/1994` CollisionOperator
  `sigma[None,:,:,:]`→`sigma[None]`; `angular_flux.py:131` `"n,ngij->gij"`→
  `"n,ng...->g..."`; any sibling `[None,:,:,:]`/`ngij`/`gxy` in scattering/sources.
- *Coupled to the lever (sigma rank must match ψ rank → C2):*
  `material_xs_field.py:369-375` `(ng,nx,ny)` reshape → `(ng,*spatial_shape)` +
  the per-material `(ix,iy)` verbs + `cells_by_material` → ndim-generic
  (`np.nonzero(mat_map)` → ndim index tuple); `geometry.py __init__`
  `(nx,ny)`/`ny=1` normalization (`volumes`/`mat_map` reshape, `is_1d`).
- *Structured 1-D (the rank change body → C2):* the `operator.py` 1-D matvec
  triple `_compute_LpC`/`_compute_decomposition`/`_compute_LpC_transpose`
  (`np.zeros((ng,N,nx,ny))` ny=1, `psi_view[:,:,0,0]` pole seed →`[:,:,0]`,
  `[:,global_dir,i,0]` writes, `sigma_t[:,:,0]`, `volumes[:,0]`,
  `.transpose(1,0,2,3)[...,0]`, re-pad `[...,None]`); `sweep.py
  _sweep_1d_unified` entry/exit (`[:,:,:,0]` drop + `[:,:,None]` re-add).
- *Solver (C2):* `solver.py` initial-guess `np.ones((ng,nx,ny))`, shape validators
  `expected=(N,ng,nx,ny)`, rate reductions `reshape(nx*ny,ng)`→`prod(spatial_shape)`.
- *Tests (C2):* slab/sphere/cylinder tests hardcoding `(N,ng,nx,1)`.

**Latent (NOT the phantom, NOT load-bearing, leave for a separate cleanup):**
`to_flat.size != n_unknowns_flat` (boundary block: `n_unknowns_flat` counts
outflow-only ordinates, `BoundaryFlux` stores all N). `n_unknowns_flat` has ZERO
production callers. Bulk parts coincide → rank change neither causes nor fixes it.

## Staging (commits, each gated)
- **C1 — local-operator rank-adaptive prep** (bit-identical on 2-D AND current
  rank-2 1-D; ALL tests green). `sigma[None]`, ellipsis einsums. Checkpoint.
- **C2 — THE LEVER + coupled + structured 1-D + tests** (the rank change). Flip
  `_shape_for_mesh`; fix material_xs, geometry `__init__`, the 1-D matvec triple,
  `_sweep_1d_unified`, solver factories/validators; migrate 1-D tests. Curvilinear
  sequenced LAST within C2, most-gated. 2-D bit-identical; 1-D re-validated via
  `k_inf=νΣ_f/Σ_a` + MMS (NOT old-vs-new bytes).
- **C3 — sweep DAG d-generic** (`OctantLabel(signs:tuple)`, `from_cartesian_2d`→
  `from_cartesian`, `diamond.cell_kernel_batch` `denom=σt+Σ_a s[a]`,
  `SweepCellSlice` per-axis). 2-D bit-identical (sum in x,y order); 3-D shape-admit.
- **C4 — boundary inventory #220** (`_resolve_bcs`→`dict[FaceLabel,op]` loop,
  `boundary_face_layout`→`face_labels`; retire named `bc_xmin/...`).
- **C5 — 3-D shape-admission pins + docs** (synthetic-3-D pack/level pins, no
  compute; extend `docs/theory/index_convention.rst` `(N,ng,nx,ny)`→`(N,ng,*spatial)`).

## Gates (test-architect, paths confirmed on `main`)
- **2-D bit-identity moat (every commit):** `tests/sn/sweep/cartesian_2d/test_2d_octant_sweep_equivalence.py`
  (nulp=64) + `tests/sn/regression/test_dd_regression.py` `-W error::DriftWarning`.
- **Sentinels WITHOUT `-O`** (Mode 8): `pytest -m sentinel`. Incl. per-ordinate
  flat-flux `(L+C)·1=σ_t·1` (ERR-026), closed-form `test_kinf_homogeneous.py`.
- **1-D re-validation (C2):** `k_inf=νΣ_f/Σ_a` (`tests/sn/eigenvalue/`), MMS
  (`tests/sn/verification/mms/` slab+sphere+cyl), SI≡Krylov twin.
- **Fast per-commit (`-O`):** `tests/sn/primitives/ tests/sn/sweep/core/
  tests/sn/sweep/slab/ tests/sn/sweep/curvilinear/ tests/sn/solve/`.
- **L16 wall-clock:** the `Σ_axis` stays a `d≤3` vectorized loop, NO per-cell Python.
- Non-square `(5,7)` for Mode-6 `order=` pins; extent-1 `(1,7)`/`(nx,1)` for #214.
- Default `python -O`; env `PYTHONPATH=$worktree .venv/bin/python`. NO
  `continuous_get` (#212). Standing reds (don't block, confirm red at clean HEAD):
  #206 cyl-matvec, #195 MMS@160, #212, #214.

## References
- `.claude/plans/nd_foundation.md` (Phase-6 design-of-record §1–§2),
  `.claude/plans/sn_development_sequence.md` (Phase 6 + verification spine).
- explorer memo `sn_phantom_axis_rank_change_audit.md`; `[[lessons-L17]]`
  (crosswalk before carve), `[[lessons-L24]]` (re-characterise at pickup),
  `[[feedback-aggressive-retirement]]`.
