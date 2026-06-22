---
name: issue-240-d5b-s3-owed2-scan-closeout
description: #240/#38/#37 D5b-S3 OWED-2 — the d=1 LD CumprodScan moment-source fix (the SWEEP-side analog of the ERR-061 matvec frame fix). The d=1 LD source-iteration scan crashed on the now-active spatial-moment source (broadcast (N,ng,nx,2) vs scalar (1,1,nx)) — a REGRESSION reds-ing 3 tests. FIX = (1) scan_xV's slope fold single-sourced through a new D1ClosedForm._slope_fold (shared with schur_xV) + new scan_slope_face_source (the face-chain b slope contribution) + scan_reconstruct (per-cell ψ̄,ψ̂ — the scan twin of the per-cell Schur); (2) the SAME octant_moment_frame_signs involution the DAG _CellSolve uses, applied on the scan via _reframe (source global→sweep IN, moment ψ̄/ψ̂ sweep→global OUT, scalar face stays sweep-frame). 3 red tests GREEN; scan≡DAG to 1e-16; diffusion limit recovered on BOTH SI(scan) AND Krylov(matvec) nx=4 rel=4.11%<5%; GATE 4 DD/Step byte-id 513/1/4 IDENTICAL; ZERO new failures; branch-free grep clean.
metadata:
  type: project
---

# #240 D5b-S3 OWED-2 — the d=1 LD moment SCAN (the sweep-side ERR-061 analog)

**Branch** `feature/sn-space-angle-tier2`. **NOT committed** (main agent commits).
Host env `.venv/bin/python`; canonical `python -O -m pytest`; NEVER all `tests/sn` (#212).

## STATUS — DONE + VERIFIED. The last D5b-S3 blocker is closed.

The d=1 LD CumprodScan source-iteration / scan path now carries the
spatial-moment iterate φ̂ + the threaded scattering-slope source `Σ_s·φ̂`,
exactly mirroring the DAG/matvec frame fix (ERR-061). The 3 red tests are GREEN,
the diffusion limit is recovered on the SI/scan path (not just Krylov), and the
DD/Step negative control is byte-identical.

## ROOT CAUSE (regression, NOT a new feature)

Pre-S3 the LD iterate was scalar, so the d=1 SI path (selected by
`CumprodScan.supports` = `is_1d and is_affine_scannable`) worked with a flat
source. S3 made the iterate moment-valued (`(N,ng,nx,2)`) → the scan source
seam (`_OneDimScanWalk._run`, the slab joint-batch `source_emission(QV,...)`)
broadcast-crashed `(N,ng,nx,2)` vs the scalar `ordinate_scan` recurrence. A
REGRESSION on a previously-valid path → a commit blocker.

## THE FIX (single-sourced; the SWEEP analog of the matvec frame fix)

Two coupled pieces, both single-sourced (Pattern 2):

### 1. The slope fold — single-sourced through `D1ClosedForm` (`orpheus/sn/spatial/_ubld.py`)

NEW private `D1ClosedForm._slope_fold(V, s_hat)` → `(mu_Adown, d2p,
eff_source_shift, slope_source)` = THE ONE site of the `s_hat` slope-row algebra.
`schur_xV` (matvec/per-cell) was REFACTORED to call it; the scan calls it too.

NEW `D1ClosedForm.scan_slope_face_source(V, s_hat, inverse_denom, w)` → the
SLOPE source's contribution to the face-chain affine source `b`:
`b = source_emission(s_bar, inv, w) + scan_slope_face_source(...)`
where the slope term = `θ·s_hat/D₂' − eff_source_shift·inv/w`. This is a pure
ADDITION on top of the existing flat-source `source_emission` (DD/Step never
reach it → byte-identical).

NEW `D1ClosedForm.scan_reconstruct(V, s_bar, s_hat, psi_in)` → `(ψ̄, ψ̂)` per
cell from the chained upstream face — the scan twin of the per-cell Schur (it
just calls `schur_xV` + the slope reconstruction; ONE algebra source). The
KEY math fact (verified): with a slope source, `cell_average(ψ_in, ψ_out, w) ≠
ψ̄` (the convex blend is only the flat-source closure). So the scan can NOT use
`cell_average` for ψ̄ — it MUST reconstruct ψ̄,ψ̂ from `schur_xV` using the
chained ψ_in. The FACE chain (`a·ψ_in + b`) propagates `ψ_out = ψ̄+ψ̂`
correctly via the slope-augmented `b`; the cell moments come from the Schur.

NEW LD scheme method `moment_scan_closure(abs_mu, A_down, V, sig_t)` →
`D1ClosedForm` (the LD-only opt-in handle; rebuilds `g = |μ|A_down/V` from the
SAME geometry the cached `(a, inverse_denom, w)` come from, so the scan's flat
`b` and its slope correction agree to FP).

### 2. The frame involution — the SAME map the DAG `_CellSolve` uses (`orpheus/sn/loss_representation.py` `_run` slab branch)

The cell kernel works in the per-ordinate SWEEP frame; the iterate φ̂ + the
scattering source `Σ_s·φ̂` live in the GLOBAL frame. The `_run` slab moment
branch applies `octant_moment_frame_signs((direction_sign,), per_axis)` via the
SHARED `_reframe` helper (NOT a scan twin — the exact same primitive the d≥2
`_CellSolve` and the d=1 `_apply_walk` matvec ride):
- source moments `QV_full_chain` global→sweep IN (`_reframe(QV, frame_signs)`);
- reconstructed `(ψ̄, ψ̂)` `mom_sweep` sweep→global OUT (`_reframe(mom, frame_signs)`);
- the scalar OUTGOING FACE stays sweep-frame (it propagates along the chain,
  never crossing into the global iterate). The d=1 face cochain is 2^{d-1}=1
  (scalar) → no reframe.
Backward sweep flips the slope (`frame_signs=[1,-1]`) so `φ̂=Σ_n w_n ψ̂_n` does
not cancel forward against backward ordinates → the diffusion limit holds.

### Entry (`_run`): moment-axis threading + flat-source lift
- `moment_tail = face_moment_tail(per_axis**ndim)`; `is_moment = moment_tail != ()`.
- `angular_flux`/`scalar_flux` allocate the trailing `2^d` axis; `QV` broadcasts
  V over it.
- A FLAT scalar source (`Q.ndim == 3`, the two-paths oracle / a manufactured Q̄)
  is LIFTED onto slot 0 (slope 0), discriminated by RANK — the SAME contract the
  DAG `_ubld_system` uses. A genuine moment source (the SI scattering iterate)
  rides through.
- The slab loop BRANCHES on `is_moment`: the flat (DD/Step) arm is the existing
  body VERBATIM (byte-identical); the moment (LD) arm is the new slope-source scan.

## WHY THE FACE/CELL SPLIT (the load-bearing math insight)

For flat-source LD, `ψ̄ = (1−w)ψ_in + w·ψ_out` (convex blend, the closure). With
a slope source, ψ̄ and ψ_out DECOUPLE. The correct chain:
1. FACE chain `ψ_out = a·ψ_in + b_face` with `b_face` = flat emission + slope
   term (so the next cell's ψ_in is the correct dense `ψ̄+ψ̂`).
2. CELL moments `(ψ̄, ψ̂)` from the per-cell Schur using the chained ψ_in (NOT
   `cell_average`).
Verified against a from-scratch dense d=1 chain (FACE/PBAR/PHAT all match to
1e-12) AND against the live DAG (scan≡DAG to 1e-16 on a 2G-het non-flat config).

## GATES (L12 paste-back)

3 red tests → GREEN (full LD slab file):
```
$ .venv/bin/python -O -m pytest tests/sn/verification/mms/test_mms_ld_slab.py -q
7 passed, 1 warning in 16.67s
```
(`test_sn_1d_slab_ld_mms_converges_second_order`, `test_sn_1d_slab_ld_mms_krylov_matches_si`,
`test_ld_two_paths_scan_equals_dag_oracle` were the 3 reds — now green; +4 already-green.)

GATE 4 (DD/Step byte-identity — the negative control):
```
$ .venv/bin/python -O -m pytest tests/sn/sweep/core tests/sn/solve \
    -W "error::tests.sn.regression._regression_assert.DriftWarning" -q
513 passed, 1 skipped, 4 xfailed, 2 warnings in 85.22s
```
IDENTICAL pre/post, no DriftWarning, no golden moved (the fix is a NO-OP at
per_axis=1: `face_moment_tail(1)==()` → `is_moment=False` → the existing flat
slab body runs verbatim; `octant_moment_frame_signs(_, 1)` returns None).

Diffusion repro on BOTH solvers (nx=4, σ_t=40, c=0.99):
```
source_iteration  : DD=2.4069 LD=2.3080 rel=0.0411 -> PASS  (the scan/SI path)
krylov            : DD=2.4069 LD=2.3080 rel=0.0411 -> PASS  (the matvec path)
```
The SI scan and Krylov matvec converge to the SAME (4-dp identical) diffusion-
limit-consistent value, rel<0.05. The scan recovers the diffusion limit too.

GATE 2 (2-D LD MMS — NOT touched, the negative control for the 2-D path):
```
$ .venv/bin/python -O -m pytest tests/sn/verification/mms/test_mms_ld_2d.py -q
3 passed, 1 warning in 332.40s
```

spatial dir (the LD foundation + primitive):
```
$ .venv/bin/python -O -m pytest tests/sn/spatial -q
70 passed, 1 warning in 66.70s
```
spatial + operators (the documented PRE-EXISTING reds, ZERO new):
```
$ .venv/bin/python -O -m pytest tests/sn/spatial tests/sn/operators -q
7 failed, 575 passed, 4 skipped, 1 xfailed
```
The 7: sphere 1-D matvec SPH ×3 (`test_vacuum_bulk_bit_identical_1d[0/1/2-SPH]`)
+ `Face 'ymin' mu_y` ×2 (`test_2d_mesh_resolution`,
`test_two_d_cartesian_loss_action_returns_result`) + sphere apply ×2
(`test_sphere_{1,2}g_apply_bit_identical`) — the EXACT set the prior closeouts
(unified-matvec, diffusion-limit) pinned as pre-existing (git-stash confirmed
there). My changes touch ONLY the d=1 slab scan path (none of the 7 exercise it:
they are sphere matvec + 2-D mu_y BC). ZERO new failures.

numerics + transport + sweep cache:
```
$ .venv/bin/python -O -m pytest tests/numerics tests/transport tests/sn/sweep/core/test_sweep_cache.py -q
875 passed, 1 skipped, 1 warning in 1.47s
```
eigenvalue + l1_analytical + primitives (DD eigenvalue rides the SI scan):
```
$ .venv/bin/python -O -m pytest tests/sn/eigenvalue tests/sn/l1_analytical tests/sn/primitives -q
354 passed, 6 warnings in 164.07s
```
LD foundation + _ubld refactor pins:
```
$ .venv/bin/python -O -m pytest tests/sn/spatial/test_ld_ubld_symbolic.py \
    tests/sn/spatial/test_ld_ubld_primitive.py tests/sn/spatial/test_linear_discontinuous.py \
    tests/sn/spatial/test_ld_slope_frame.py -q
6 + 11 + 21 + 2 = 40 passed
```

Branch-freeness (single-source discipline — empty = clean):
```
$ grep -rnE "len\(s_axes\) *> *1|isinstance.*(Discontinuous|Diamond)" \
    orpheus/sn/loss_representation.py orpheus/sn/spatial/_ubld.py \
    orpheus/sn/spatial/linear_discontinuous.py orpheus/sn/spatial/scheme.py \
    orpheus/sn/spatial/scan.py orpheus/sn/sweep_graph.py
  (no output — ZERO matches)
```
The moment gate is the established trait `face_moment_tail(per_axis**ndim) != ()`,
NOT isinstance, NOT `len(s_axes)>1`.

## HOW THE INVOLUTION IS SINGLE-SOURCED BETWEEN DAG + SCAN

The scan does NOT define its own frame map. It calls the SAME
`octant_moment_frame_signs((direction_sign,), per_axis)` (`_ubld.py`) +
`_reframe(arr, frame_signs)` (`sweep_graph.py`) the d≥2 `_CellSolve`/`_CellResidual`
and the d=1 `_apply_walk` matvec ride. The d=1 `direction_sign` IS the
`signs_eff` the d≥2 sites pass. `_reframe`'s guard
(`arr.shape[-1] != frame_signs.shape[0]`) makes a scalar face / DD source pass
untouched. So the sweep⇄global involution lives in exactly ONE place; the scan
is a new CONSUMER, not a twin. The slope-fold `s_hat` algebra similarly lives in
ONE place (`D1ClosedForm._slope_fold`), shared by `schur_xV` (matvec/per-cell),
`scan_slope_face_source`, and `scan_reconstruct`.

## FILES CHANGED (mine, this session)
- `orpheus/sn/spatial/_ubld.py` — `D1ClosedForm._slope_fold` (NEW, the shared
  `s_hat` fold); `schur_xV` REFACTORED to call it (byte-equivalent — 38 existing
  primitive/symbolic/LD tests stay green); `scan_xV` docstring (slope note);
  `scan_slope_face_source` (NEW); `scan_reconstruct` (NEW).
- `orpheus/sn/spatial/linear_discontinuous.py` — `moment_scan_closure` (NEW, the
  LD-only opt-in d=1 closed-form handle for the moment scan).
- `orpheus/sn/loss_representation.py` — `_OneDimScanWalk._run`: entry moment-tail
  + flat-source rank-lift; the slab joint-batch loop BRANCHED into a flat
  (byte-identical) arm + a moment (LD slope-source) arm with the frame involution.
- `docs/theory/discrete_ordinates.rst` — `ld-ubld-moment-scan` stub (label +
  `:eq: ld-ubld-moment-scan-source` + `:mod:`/`:meth:` cross-refs + archivist
  TODO; Sphinx build exit 0, label/eq/TODO confirmed in HTML + Nexus graph).
NOT mine (prior uncommitted D5b-S3 sessions): the rest of the diff stat
(scattering, solver, operator, _bases, angular_flux, harmonic_moment_field,
boundary_operator, sweep_graph `_reframe`/`_CellSolve` frame field — all from the
unified-matvec + diffusion-limit sessions).

## EDGE NOTE (acceptable, guarded)
The curvilinear `else` branch in `_run` still writes scalar `angular_flux` — it
would break if `is_moment=True`. But curvilinear LD is UNPUBLISHED and REJECTED
at the cache build (`affine_scan_coefficients` raises;
`test_ld_curvilinear_scan_rejected` passes), so the moment-curvilinear pairing is
unrepresentable — the guard makes the illegal state unreachable. If a future
curvilinear LD lands, the curvilinear `_run` arm must gain the same moment
threading.

## OWED (the rest of D5b-S3 / D5b)
- archivist DISPATCH emitted (rich narrative for `ld-ubld-moment-scan`) —
  followup:false.
- The non-vacuum domain-inflow moment trace (S4 boundary widening) is unrelated.
- #246 (the typed `SpatialMomentSpace` predicate replacing the 4 value-based
  shape probes) is the standing follow-up; my `is_moment = moment_tail != ()`
  adds a 5th `face_moment_tail`-derived gate (consistent with the existing 4,
  closed by #246 together).

## LESSON (candidate for algebra-of-record / coding-elegance)
A moment-carrying parallel-prefix SCAN is NOT a drop-in widening of the scalar
scan: with a slope SOURCE the convex face-blend closure `ψ̄=(1−w)ψ_in+w·ψ_out`
DECOUPLES from ψ_out (it holds only for the flat-source LD), so the scan must (a)
propagate the FACE chain `ψ_out=ψ̄+ψ̂` via a slope-augmented affine `b`, and (b)
reconstruct `(ψ̄,ψ̂)` per cell from the chained ψ_in via the per-cell Schur —
NOT via `cell_average`. Conflating the two (using `cell_average` for ψ̄ on the
moment scan) gives the WRONG cell average while the FACE chain looks right. Map
the closure-validity-regime (flat-source-only convex blend) as an explicit
crosswalk row whenever a carve threads a SOURCE moment through a recurrence whose
reconstruction was derived flat-source. The single-source win: the same
`D1ClosedForm._slope_fold` powers the matvec per-cell Schur AND the scan, and the
same `octant_moment_frame_signs`/`_reframe` involution powers the DAG AND the
scan — the scan is a CONSUMER of the matvec's machinery, never a twin.
