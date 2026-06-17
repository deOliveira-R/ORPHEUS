---
name: issue-240-d5b-s3-unified-matvec-closeout
description: #240/#38/#37 D5b-S3 (merged A+B) — the UNIFIED moment matvec. Landed the full architecture (branch removal + d-generic moment matvec + carrier threading): cell_kernel_batch/residual_kernel_batch collapse to ONE d-generic dense path (d=1 scalar arm + d≥2 raise BOTH retired); the len(s_axes)>1 moment gate at _CellSolve/_CellResidual is GONE (pure scheme trait); ZERO scheme-isinstance; every iterate/source/probe/scattering carrier threads the spatial-moment axis shape-agnostically via face_moment_tail (DD byte-identical). DD/Step GATE 4 = 513P/1skip/4xf (preserved). 2-D LD fixed-source RUNS end-to-end + 2-D LD MMS O(h²) PASSES (GATE 2 2-D). BUT two coupled DEEP issues remain → numerics-investigator: (1) the thick-cell DIFFUSION LIMIT is NOT recovered (LD converges to DD only with refinement: nx=4→39%, nx=16→8%, nx=64→0.8% — the flat-source signature persists despite the slope source threading + a verified matvec≡sweep consistency); (2) the d=1 SCAN (CumprodScan) moment threading is NOT done (1-D LD MMS crashes on the moment source). The matvec architecture is structurally sound (round-trip verified to 1e-16); the diffusion-limit physics needs the UBLD asymptotic analysis.
metadata:
  type: project
---

# #240 D5b-S3 (merged A+B) — the UNIFIED moment matvec

**Branch** `feature/sn-space-angle-tier2`. **NOT committed** (main agent reviews + commits).
Host env, `.venv/bin/python -O`. Canonical `python -O -m pytest`; NEVER all `tests/sn` (#212).

## STATUS — the ARCHITECTURE is DONE; the diffusion-limit PHYSICS + the d=1 scan are OWED

The USER-settled decision (unify the matvec across all d; net-remove branches) is
IMPLEMENTED IN FULL at the architecture level. The branch-removal headline is
achieved, DD byte-identity is preserved, the 2-D LD path works end-to-end, and the
unified moment matvec is verified consistent (matvec≡sweep round-trip to 1e-16).
BUT the thick-cell diffusion limit is NOT recovered and the d=1 scan is NOT threaded
— both are coupled DEEP numerical issues (the slope-source physics) that need
numerics-investigator. I stopped at the wrong-answer cascade per the brief's
discipline ("If isolation runs longer than ~30 minutes, dispatch
numerics-investigator").

## THE BRANCH-REMOVAL GREP (the headline deliverable — CLEAN)

```
$ grep -rnE "len\(s_axes\) *> *1|isinstance.*(Discontinuous|Diamond)" \
    orpheus/sn/sweep_graph.py orpheus/sn/loss_representation.py \
    orpheus/sn/spatial/linear_discontinuous.py orpheus/sn/spatial/diamond.py \
    orpheus/sn/spatial/scheme.py orpheus/sn/operator.py
  (no output — ZERO matches)
```

The `len(s_axes) > 1` moment gate at `_CellSolve.cell` (sweep_graph:883) AND
`_CellResidual.cell` (:929) is GONE — the emit/probe path is now the pure scheme
trait `spatial_basis_per_axis > 1` (the dimension conjunct that admitted d=1 LD to
the scalar path is removed). The d=1 scalar matvec kernel twin (`_kernel_terms`) is
retired. The d≥2 `_CellResidual` raise is retired. The only remaining `len(s_axes)`
uses are `d = len(s_axes)` (computing the dimension, NOT a branch). ZERO
`isinstance(scheme, ...)` in the sweep/matvec layer.

## WHAT LANDED (the complete architecture)

### 1. The unified d-generic moment kernel (`linear_discontinuous.py`)
- `cell_kernel_batch` / `residual_kernel_batch` COLLAPSED to ONE d-generic dense
  UBLD path (`_ubld_system` → `per_cell_solve` / `A·ψ⃗ − R`) for EVERY d. The d=1
  arm now returns the FULL `(N_oct, ng, n_diag, 2^d=2)` moment vector (slot 0 =
  average, slot 1 = slope). The `_kernel_terms` scalar Q̂=0 helper is DELETED.
- `_ubld_system` is SHAPE-AGNOSTIC over the source: a genuine `(…, 2^d)` moment
  source (its trailing axis IS the moment axis — discriminated by RANK vs the cell
  batch, NOT trailing-axis size, which would be fragile when n_diag==2^d) flows
  through verbatim; a flat scalar source lifts onto slot 0.
- `_ubld_inflow` appends the explicit `2^{d-1}` face-moment axis at d=1 (where the
  walk's cochain has NO moment axis, `face_moment_tail(1)==()`), keyed on
  `face_moment_tail(2^{d-1})==()` (deterministic on d, NOT a shape probe — a scalar
  face's trailing n_diag==1 is indistinguishable from a 2^{d-1}==1 moment by shape).
- `_ubld_outgoing_faces` drops the d=1 face-moment axis (tail `()`) to match the
  walk's scalar d=1 cochain.
- **⭐ THE MASS-NORMALIZATION (the matvec/sweep moment-source consistency, NEW):**
  `residual_kernel_batch` divides the residual `A·ψ⃗ − R` by the diagonal Legendre
  mass `M_ii = θ^{|i|}`. The UBLD RHS folds the cell source MASS-WEIGHTED
  (`R_source = M·S⃗`, the test-function projection), but `A = (L+C) − S` subtracts
  `S.apply(ψ)` RAW at the OperatorSum level — so `(L+C)ψ` must be put in raw
  per-moment units or the slope rows disagree by θ. VERIFIED: `(L+C)ψ` (M-normed)
  used as a raw moment source recovers ψ exactly (round-trip 1e-16); the d=1 closed
  form still matches the dense solve (`test_production_kernel_equals_dense`).

### 2. The cell-emit φ̂ accumulation (`sweep_graph.py`)
`_CellSolve.cell` keeps the FULL `2^d` axis (no drop) and emits via the
spatial-moment-AGNOSTIC einsums `"ng...,n->g..."` / `"nlm,ng...,n->lmg..."` (the
SAME convention `_SweepEmit.pure_z` already used) — the trailing moment axis rides
through with NO per-axis branch. `_CellResidual.cell` raise retired; the probe
carries the moment axis via `psi_avg_probe_octant[(:,:,*cell_idx)]`.

### 3. The buffer widening (`loss_representation.py`) — `_spatial_moment_tail`
NEW `_LossRepresentation._spatial_moment_tail` property (= `face_moment_tail(per_axis^d)`,
the bulk-field sibling of `_n_face_moments`). Every bulk buffer appends it:
`_OctantWalk.loss_action` (LpC + Q_zero), `MovingFrontierWindow`/`FullFieldWavefront`
`_sweep_interior`/`_loss_action_interior` (angular_oct + LpC_oct + n_face_moments to
walk_windowed), `_sweep_scheduled` (angular/scalar/moment buffers), the d=1
`_apply_walk` Cartesian matvec (`out_g_first` + the moment probe via
`swapaxes(0,1)[:,:,None]`). DD/Step (per_axis==1) → `()` tail, byte-identical.

### 4. The carrier threading (shape-agnostic, NOT `if LD-d≥2`)
- **Iterate carriers SELECT the axis:** `_unwindowed_cold_start` (NEW, the 1-D/curvilinear
  AngularFlux sibling of `_windowed_cold_start`) + `_windowed_cold_start` both pass
  `spatial_moments=scheme.spatial_basis_per_axis`. Both fixed-source + eigenvalue SI
  paths use them.
- **The typed output wraps:** `InvertibleOperator._solve_timed_full_field`,
  `_GaussSeidelResolvent.solve_moments`, `_OneDimScanWalk.loss_action`/`loss_action_decomposed`
  all pass `spatial_moments=per_axis` to `AngularFlux`/`HarmonicMomentField.from_mesh*`.
- **The external source moment-lift:** `_build_fixed_source_rhs` + `_lift_external_source_to_moments`
  (NEW helper) lift the flat external source onto slot 0 (Q̂=0); `_average_moment_scalar`
  (NEW) reduces the Solution scalar flux to the cell-average moment.
- **The Krylov path:** `n_dof *= per_axis^d`; template = `_unwindowed_cold_start`;
  `_average_moment_scalar` reduction. The flat ravel (`to_flat`/`from_flat`) tracks the
  template shape automatically (VERIFIED — `replace(template.bulk, values=...)` preserves
  the moment space).

### 5. The scattering Σ_s⊗I flow (`scattering.py`, `_bases.py`, `angular_flux.py`, `harmonic_moment_field.py`)
- `BulkField.spatial_moments_per_axis` (NEW property) reads the width OFF the field's
  space (`find_factor(SpatialMomentSpace).per_axis`) — the single source of truth.
- `AngularFlux.integrate_angular` + `HarmonicMomentField.scalar_flux` propagate the width
  to the child `ScalarFlux` (the `ng...` einsum is already agnostic).
- `_assemble_per_ordinate_source` allocates the iso/aniso accumulators with
  `spatial_moments=_spatial_moments_of(phi)` so the in-place `Σ_s^T·φ` `+=` doesn't
  broadcast-fail. `add_iso_source`/`add_n2n_source` preserve Q's width on re-wrap.
  `build_aniso_source` + the HarmonicMomentField aniso arm wrap with the iterate's width.
- `SNBoundaryOperator._apply_faces` zero bulk carries the input's moment width (so
  `(L+C)ψ − B.apply(ψ)` composes element-wise in the matvec).
- The 2-D `_sweep_scheduled` pure-z `Q/Σ_t` broadcasts Σ_t over the trailing moment axis.

## GATES (verbatim paste-back, L12)

### GATE 4 (negative control — DD/Step byte-identity) — PASSES, IDENTICAL pre/post
```
$ .venv/bin/python -O -m pytest tests/sn/sweep/core tests/sn/solve \
    -W "error::tests.sn.regression._regression_assert.DriftWarning" -q
513 passed, 1 skipped, 4 xfailed, 2 warnings in 84.04s
```
The S2 baseline `513/1/4` is PRESERVED pre==post, NO golden `.npy` moved. The
`test_2d_ld_sweep_runs_genuine_ld_not_dd` test (in this gate's set) now PASSES
(2-D LD runs + differs from DD). DD byte-identity falls out of the
`face_moment_tail(per_axis^d)==()` formula, NOT an `if`.

### GATE 2 (2-D leg) — PASSES (the 2-D LD MMS O(h²) + DD≠LD discrimination)
```
$ .venv/bin/python -O -m pytest tests/sn/verification/mms/test_mms_ld_2d.py -q
3 passed, 1 warning in 324.54s
```
`test_ld_2d_converges_second_order_smoke` (O(h²) with the slope source active) +
`test_ld_2d_two_paths_ffw_equals_mfw` + `test_dd_and_ld_2d_converge_to_different_values`
ALL green. (⚠ 324s — the d≥2 dense per-cell solve is slow; #227 perf.)

### Pre-existing reds (spatial+operators) — EXACTLY 7, ZERO new (git-stash confirmed)
```
$ .venv/bin/python -O -m pytest tests/sn/spatial tests/sn/operators -q
7 failed, 573 passed, 4 skipped, 1 xfailed
```
The 7: sphere 1-D matvec SPH ×3 (`test_vacuum_bulk_bit_identical_1d[*-SPH]`) +
`Face 'ymin' mu_y` ×2 (`test_2d_mesh_resolution`, `test_two_d_cartesian_loss_action_returns_result`)
+ sphere curvilinear apply ×2 (`test_sphere_{1,2}g_apply_bit_identical`). `git stash`
confirmed identical at clean tree — ZERO introduced by this carve.

### LD foundation + dispatch + numerics/transport — ALL GREEN
```
$ .venv/bin/python -O -m pytest tests/sn/spatial/test_ld_ubld_symbolic.py \
    tests/sn/spatial/test_ld_ubld_primitive.py tests/sn/spatial/test_linear_discontinuous.py \
    tests/sn/sweep/core/test_unified_sweep_dispatch.py -q
59 passed
$ .venv/bin/python -O -m pytest tests/numerics tests/transport -q
848 passed
```

### GATE 3 (the genuine Mode-9 — SI≡Krylov on (L+C−S_full)) — PASSES (smoke), strong positive
A 2-D LD fixed-source smoke (4×3, level_symmetric S4, c=0.7 self-scatter, vacuum):
```
2-D LD SI vs Krylov max rel diff = 4.99e-11   (≈ inner_tol — PASSES)
```
SI-with-lag and Krylov (the now-UNRAISED d≥2 LD matvec) converge to the SAME flux on
the identical `(L+C−S_full)` operator — two structurally-different code paths agreeing
to solver tol. This is STRONG evidence the unified moment matvec is structurally
CORRECT (the matvec twin of the sweep). ⚠ Per vv Mode-9: this proves SI≡Krylov on the
SAME operator — it does NOT prove that operator is the diffusion-limit-consistent one
(that is GATE 1 / OWED-1, a SEPARATE claim). The test-architect's gate config
(anisotropic + diagonal cubature + non-zero scatter + 2G-asym + non-square) should be
landed as the formal `test_ld_2d_krylov_matches_si_full_operator` once OWED-1 is fixed.

### GATE 1 (the thick-diffusion tripwire) — STILL XFAILED (the diffusion limit is NOT recovered)
```
$ .venv/bin/python -O -m pytest \
    "tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_thick_diffusive_limit_xfail" -q
1 xfailed, 1 warning
```
I deliberately LEFT the xfail (did NOT flip to a passing gate) because flipping it
while the physics is broken would be a FALSE GREEN. The tripwire correctly stays
xfailed until the diffusion-limit physics is fixed.

### GATE 2 (1-D leg) + the d=1 two-paths gate — FAILING ON THE SCAN MOMENT THREADING
```
$ .venv/bin/python -O -m pytest tests/sn/verification/mms/test_mms_ld_slab.py -q
3 failed, 2 passed, 1 xfailed
  FAILED test_sn_1d_slab_ld_mms_converges_second_order — broadcast (16,1,20,2) (the scan QV·V)
  FAILED test_sn_1d_slab_ld_mms_krylov_matches_si — same (the SI leg crashes)
  FAILED test_ld_two_paths_scan_equals_dag_oracle — scan path not moment-threaded
```
All 3 are the d=1 SCAN (CumprodScan) not yet threading the moment source (item 5).

## ⭐⭐ THE TWO OWED ITEMS (coupled — both need the diffusion-limit slope-source physics)

### OWED-1 (LOAD-BEARING): the thick-cell DIFFUSION LIMIT is NOT recovered
The Krylov 1-D LD matvec RUNS, carries the φ̂ moment, and the matvec≡sweep is
VERIFIED CONSISTENT (round-trip 1e-16; `(L+C)ψ` M-normed == the raw source the sweep
needs). The scattering slope source `Σ_s·φ̂` IS active (φ̂ ≠ 0 in the converged
iterate). BUT the thick-cell limit is NOT recovered — LD converges to DD only with
MESH REFINEMENT, the EXACT flat-source signature Inc C is supposed to eliminate:
```
nx=4:  DD_mid=2.41  LD_mid=1.47  rel=0.39    (thick cells — should match, doesn't)
nx=16: DD_mid=2.36  LD_mid=2.18  rel=0.08
nx=64: DD_mid=2.36  LD_mid=2.34  rel=0.008   (thin cells — converges, as flat LD does)
```
σ_t=40, c=0.99, so σ_t·h=10/cell at nx=4. The converged φ̂≈-0.02 is TINY relative to
φ̄≈1.5 — for a diffusion profile the within-cell slope should be MUCH larger. So the
slope is being computed but the slope-scattering coupling does NOT produce the
diffusion-limit-consistent operator. **This is a wrong-answer cascade in the UBLD
diffusion-limit physics** — the architecture is sound (verified), the bug is in the
slope-source coupling magnitude/sign/ordering. Needs the algebra-of-record (Adams-2001
UBLD asymptotics, LMM-1987 four-limits) + probe-cascade. Candidate suspects:
  (a) the moment ORDERING — does `integrate_angular`'s slope moment align with the
      UBLD slope row? (the Kronecker `[bar, ŷ, x̂, x̂y]` order; d=1 is `[bar, hat]`).
  (b) the Σ_s⊗I vs Σ_s⊗M convention — I chose Σ_s⊗I (raw) + matvec M-normalization;
      maybe the diffusion limit needs a different moment-source weighting.
  (c) a missing cross-cell slope-continuity / a sign on the slope-row scattering.

### OWED-2: the d=1 SCAN (CumprodScan) moment threading (item 5 of the directive)
The d=1 production sweep (CumprodScan `_run`) is NOT threaded for moments — it
produces a scalar `angular_flux` and `QV_per_ord = Q*V` broadcast-fails on the
`(N,ng,nx,2)` source. Item 5 requires `scan_xV` to gain an `s_hat` slope-source arg,
single-sourced with `schur_xV`'s fold through `D1ClosedForm`, and `_run` to
reconstruct ψ̂ per cell + write the moment iterate. This is COUPLED to OWED-1 (the
slope physics must be right first; the scan is the perf path for the SAME operator).
Until done, the 1-D LD SI path crashes (the 1-D MMS uses the scan).

## FILES CHANGED (mine, this session)
- `orpheus/sn/spatial/linear_discontinuous.py` — kernel collapse (one d-generic moment
  path; `_kernel_terms` deleted) + `_ubld_system` shape-agnostic source + `_ubld_inflow`
  d=1 face-axis append + `_ubld_outgoing_faces` d=1 tail + `residual_kernel_batch`
  M-normalization. `face_moment_tail` import.
- `orpheus/sn/sweep_graph.py` — `_CellSolve.cell` agnostic-einsum emit (gate retired) +
  `_CellResidual.cell` raise retired. Dropped `AVERAGE_MOMENT` import (now unused).
- `orpheus/sn/loss_representation.py` — `_spatial_moment_tail` property + all bulk-buffer
  widenings (octant-walk LpC/Q_zero, both `_sweep_interior`/`_loss_action_interior`,
  `_sweep_scheduled`, the d=1 `_apply_walk` matvec + its 3 output wraps) + pure-z Σ_t
  moment broadcast.
- `orpheus/sn/operator.py` — `_solve_timed_full_field` typed-wrap `spatial_moments`.
- `orpheus/sn/solver.py` — `_unwindowed_cold_start` + `_lift_external_source_to_moments`
  + `_average_moment_scalar` (NEW helpers) + `_GaussSeidelResolvent.solve_moments` wrap +
  `_build_fixed_source_rhs` lift + the 4 cold-start sites + the Krylov n_dof/template/reduce
  + the windowed-reconstruction scalar reduce.
- `orpheus/sn/scattering.py` — `_spatial_moments_of` helper + `_assemble_per_ordinate_source`
  accumulator widening + `add_iso_source`/`add_n2n_source` re-wrap width + `build_aniso_source`
  + HarmonicMomentField aniso arm wraps.
- `orpheus/sn/boundary_operator.py` — `_apply_faces` zero-bulk moment width.
- `orpheus/transport/fields/_bases.py` — `BulkField.spatial_moments_per_axis` property.
- `orpheus/transport/fields/angular_flux.py` — `integrate_angular` width propagation.
- `orpheus/transport/fields/harmonic_moment_field.py` — `scalar_flux` width propagation.
- `orpheus/sn/spatial/_ubld.py` — docstring updates (kernel is dense, scan/per-cell use the
  closed form).
- `docs/theory/discrete_ordinates.rst` — `ld-ubld-unified-moment-matvec` stub (labels
  `ld-ubld-unified-moment-residual` + `:mod:` + archivist TODO; build exit 0, no new warnings).
- `tests/sn/spatial/test_linear_discontinuous.py` — 2 tests updated (d=1 kernel returns the
  moment vector now: `psi_avg[..., AVERAGE_MOMENT]` for the scalar avg; the d=1 face is
  scalar, no moment axis).
- `tests/sn/spatial/test_ld_ubld_primitive.py` — `test_production_kernel_equals_dense` same.
NOT mine (pre-existing uncommitted): `.claude/skills/vv-principles/*`, `.claude/plans/*`,
`docs/verification/matrix.rst` (Sphinx auto-regen), `.claude/agent-memory/elegance-enforcer/*`,
the 3 forbidden untracked.

## LESSON (candidate for the crosswalk / coding-elegance)
A carve that "unifies the matvec" by collapsing a Schur-reduced scalar path onto a
d-generic dense path inherits a SUBTLE units mismatch: the dense UBLD RHS is
mass-weighted (`M·S⃗`, the test-function projection), but `A = (L+C) − S` subtracts the
scattering source RAW at the OperatorSum. The matvec residual MUST be M-normalized
(÷ the diagonal Legendre mass) to match — else the slope rows disagree by `θ^{|i|}`
and the round-trip silently fails on the slope (the AVERAGE row coincidentally agrees
because `M_00=1`). Map the SOURCE-NORMALIZATION as an explicit Pattern-7 crosswalk row
whenever a carve routes a source through BOTH a mass-weighted cell system (sweep) AND
an element-wise OperatorSum subtract (matvec). AND: the diffusion-limit recovery is a
SEPARATE correctness claim from the matvec-consistency — a verified matvec≡sweep
round-trip does NOT imply the operator is the diffusion-limit-consistent one (the
slope-source coupling can be self-consistent but physically wrong). Gate the limit
against the structurally-independent diffusion reference, NOT just the round-trip.
