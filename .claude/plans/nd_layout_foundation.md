# N-D layout foundation — `(N, ng, *spatial)` array convention (eliminate the phantom `ny=1`)

**Branch/worktree:** `worktree-sn-nd-layout` (fresh off `main` @ `f0ffb30`).
**Mode:** main-agent direct authorship, turn-by-turn user steering
(`feedback_no_method_implementer_for_surgical_carves`). Commit per stage; gates
throughout. Standing commit exclusions: `.claude/agent-memory/`, `.mcp.json`,
hooks, `derivations/diagnostics/`.

---

## ⭐⭐⭐ POST-COMPACTION RECOVERY — READ THIS FIRST (2026-06-09)

**One line:** the production `(N, ng, *spatial)` rank-d carve (eliminate the
phantom `ny=1`) is **DONE + VALIDATED + COMMITTED**; what remains is a **deep,
manual, per-file test migration (~47 files)** + then phases C3/C4/C5.

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

**TEST MIGRATION STATUS: operators+primitives+solve+numerics+transport REGION
= 1650 passed, 0 failed (DONE).** Migrated via 3 waves of parallel general-purpose
sub-agents (brief = the validated method; reference diff `git show 1774586`). Files
done (~22): the 2 reference + wave-1 (typed_residuals, angular_flux, field_units,
timed_full_field, solution, typed_source_sinks, harmonic_moment_field, streaming_operator,
native_matvec, iteration, affine_flux_algebra) + wave-2/3 (collision, bc_extraction_matvec,
solver_components, invertible, operators_apply_typed, g_adjoint_reciprocity, fixed_source_g1,
flux_displacement_diagnostics, si_single_primitive_contract). **REMAINING = the
sweep/verification/eigenvalue region** (curvilinear + slab + 1-D verification +
sweep/core-1-D; the `cartesian_2d/*`, `test_keff_2d`, `test_mms_2d` are GENUINE 2-D →
already green, LEAVE). Method is PROVEN; gate the region (`gtimeout 1200 ... tests/sn/sweep/
tests/sn/verification/ tests/sn/eigenvalue/ --deselect ...test_heterogeneous_absolute_keff`,
buffered so use a file), fan out per failing file, verify-then-regen any `.npy`/snapshot.
**Commit chain (this session):** `1774586`(ref batch) → `94fcae5`(PRODUCTION fixes) →
`9870495`(wave-1 tests) → `<wave-2/3 tests>`. Then C2d full-suite gate, then C3/C4/C5.

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
**(1) C2c test migration (NEXT) — deep, manual, PER-FILE. NOT scriptable.**
~47 files reference rank-2 phantom. Worst offenders (operators dir gate:
160 fail / 319 pass): `test_streaming_operator` 35, `test_bc_extraction_matvec` 30,
`test_streaming_operator_decomposition` 19, `test_native_matvec` 16,
`test_operators_apply_typed` 12, `test_invertible_operator` 12,
`test_g_adjoint_reciprocity` 12, `test_collision_operator` 12,
`test_typed_residual_evaluation` 4, `test_solver_components` 1. Plus the
sweep/solve/eigenvalue/primitives/verification files (not yet gated). Find them:
`grep -rlE "\((N, ng|ng|quad\.N, ng|N_ord, ng), (nx|n_cells|nr), ny\)|\.values\[[^]]*, 0\]|:, :, :, 0\]" tests/sn/`.
  - **METHOD (per file/function):** for CONSTRUCTIONS, build via the in-scope
    mesh — `AngularFlux.zeros_on(mesh)` / `.from_mesh(arr, mesh)` or
    `np.zeros((N, ng, *<the in-scope mesh>.spatial_shape))`. For EXTRACTIONS on
    the now-rank-1 `(ng,nx)` flux / `(N,ng,nx)` field: `.values[0, :, 0]`→`[0, :]`,
    `.values[:, :, 0].T`→`.values.T`, `.values[:, 0, 0]`→`[:, 0]`,
    `[:, :, :, 0]`→drop, `.bulk.values[:, 0, 0, 0]`→`[:, 0, 0]`. `from_isotropic`
    now expects `(ng, *spatial)`. Gate per file: it's fast (operators ~5s).
  - **DO NOT REPEAT (all tried + reverted):** (a) blanket `→(N,ng,nx)` regresses
    hidden-2-D cases; (b) scripted `→*<meshvar>.spatial_shape` gives 87 NameErrors
    (mesh var varies per construction, esp. in helper fns); (c) rank-conditional
    `if ny==1 else` is low-gain/high-clutter + misses extractions. **Go per-file.**
  - `legacy_proxy_matvec` helper ALREADY fixed (`6ca0154`); the DD-snapshot
    generator `_generate_snapshots.py` already rank-d (`79f49e6`). The genuine
    2-D sweep tests (`test_sweep_graph*`, `test_cell_update_batch`) use real `ny>1`
    — they are NOT phantom; leave them.
**(2) C2d full-suite gate** once C2c green: `tests/sn/ tests/transport/ tests/numerics/`
   slices (deselect #212; sentinels w/o -O; `-W error::DriftWarning` on regression).
**(3) C3 sweep DAG d-generic** (`OctantLabel(signs:tuple)`, `from_cartesian_2d`→
   `from_cartesian(dims)`, `diamond.cell_kernel_batch` `denom=σt+Σ_a s[a]`,
   `geometry` octant build `itertools.product((-1,+1), repeat=ndim)`; 2-D
   bit-identical sum-in-x,y-order; synthetic-3-D shape admission). See `nd_foundation.md` §2.
**(4) C4 boundary inventory #220** (`_resolve_bcs`→`dict[FaceLabel,op]`).
**(5) C5 3-D shape-admission pins + docs** (extend `docs/theory/index_convention.rst`).

### Task tracker (recreate if lost): #1 C1 ✅ · #2 C2a ✅ · #3 C2b ✅ ·
#4 C2c (IN PROGRESS — the deep per-file test migration above) · #5 C2d ·
#6 C3 · #7 C4 · #8 C5. Chain blockedBy 4→5→6→7→8.

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
