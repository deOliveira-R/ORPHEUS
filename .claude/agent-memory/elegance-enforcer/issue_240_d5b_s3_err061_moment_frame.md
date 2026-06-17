---
name: issue-240-d5b-s3-err061-moment-frame
description: ERR-061 LD diffusion-limit fix review — octant_moment_frame_signs involution + _reframe shape-guard CONCERN + the d=1 _frame_signs twin DO-NOW.
metadata:
  type: project
---

# #240 D5b-S3 (ERR-061) moment-frame involution — PASS-WITH-NITS

The diffusion-limit fix: LD slope moment ψ̂_n stored in per-ordinate SWEEP frame but
summed by the angular reduction φ̂=Σ_n w_n ψ̂_n as if GLOBAL-x; backward ordinates (μ<0)
have sign-flipped sweep frame → forward+backward slopes cancelled → diffusion limit lost.
Fix = re-sign the moment vector at the octant boundary (global→sweep IN, sweep→global OUT).
GATE-4 DD byte-id 513/1/4 zero drift; diffusion limit nx=4 rel 38.9%→4.1%.

## What the fix did (3 files, all uncommitted on `feature/sn-space-angle-tier2`)
- `_ubld.py`: NEW `octant_moment_frame_signs(octant_signs, spatial_basis_per_axis)` =
  ∏_a (octant_sign_a)^{o_a} in tensor-Legendre Kronecker layout. `[1.0]` for per_axis==1.
- `sweep_graph.py`: NEW module-level `_reframe(arr, frame_signs)` (shape-guarded involution);
  `moment_frame_signs` field on `_CellSolve`/`_CellResidual`; routes source/probe (IN) +
  moment/residual (OUT) through `_reframe`. Outgoing FACE deliberately NOT re-framed.
- `loss_representation.py`: `_moment_frame_signs(signs_eff)` method (SSOT wrapper); wired at
  4 d≥2 sites; d=1 matvec uses a SEPARATE `_frame_signs(direction_sign)` closure.

## VERDICTS
**[PASS] the primitive** — math verified NUMERICALLY: involution v*v==1 all octants; d=2
`(−1,−1)`→`[1,−1,−1,+1]` (both slopes flip, cross-moment x̂y flips twice → invariant).
Textbook Pattern-7 producer-owns-the-sign. `signs_eff` maps grazing-0→+1 at the call site
(loss_rep:606) BEFORE the primitive sees it → "zero octant never reached" is STRUCTURAL not
probabilistic. Reads like the physics.

**[PASS] branch-free preserved + TIGHTENED** — grep-clean: zero `len(s_axes)>1`, zero
scheme-isinstance. The fix REMOVES the old d≥2 dimension gate (`psi_avg[...,AVERAGE_MOMENT]`
slice + the `_CellResidual` `len(s_axes)>1` NotImplementedError raise) in favor of
spatial-moment-AGNOSTIC einsums (`ng...` / `nlm,ng...`, same as `_SweepEmit.pure_z`). This
dissolves the D5b-S2 `len(s_axes)>1 and per_axis>1` vs `_n_face_moments!=1` coextensivity smell
on the solve side. Genuine architectural-forwardness GAIN.

**[PASS] 4× `_moment_frame_signs(signs_eff)` NOT a dup** — one cheap call per octant, signs
differ; the (windowed|full)×(solve|apply) grid is pre-existing structure.

**[DO-NOW NIT] the d=1 `_frame_signs` closure IS a Pattern-2 twin of `_moment_frame_signs`**
(loss_rep:2243-2246 vs :300-318). BOTH encode `per_axis==1→None else octant_moment_frame_signs`.
Only diff = `signs_eff` (d-tuple) vs `(direction_sign,)` (1-tuple) — but a d=1 octant's signs_eff
IS `(direction_sign,)`. Same decision, two spellings. Bug habitat = the EXACT ERR-061 / Phase-F#6
class this fix closes; S4 convention refinement (metric-weight / θ-mass into the frame) lands on the
method, d=1 closure keeps old convention, and 1-D LD moment axis is the LEAST byte-id-tested.
COLLAPSE: `self._moment_frame_signs((direction_sign,))` (verify reachability from `_OneDimScanWalk`;
if not a holder, hoist to shared base or module-level free func).

**[CONCERN] `_reframe` shape-guard `arr.shape[-1] != frame_signs.shape[0]` = correctness-by-
CONVENTION not by-construction.** Value-based moment-axis discriminator (the recurring
`_n_face_moments!=1` / `len(s_axes)>1` / `q.ndim>sig.ndim+1` family). Safe for S3 because:
(1) short-circuits on `None` before touching `.shape` → probe only runs under LD;
(2) under LD `psi_avg`/`psi_bar`/`residual` carry trailing `2^d` BY KERNEL CONTRACT;
(3) only non-moment array is `Q_cells` — d=1 flat-source trailing = `n_diag==1 ≠ 2`; d≥2 LD `Q`
ALWAYS built with `moment_tail` (the ~6 construction sites bake it). COLLISION (`n_diag==2^d`,
reachable on a 4×4 d=2 mesh) is UNREACHABLE only because no flat-source d≥2 LD path exists yet.
S4 (non-vacuum moment trace / external moment source) is when it first becomes reachable → silent
sign-flip (inverted ERR-061), NO exception. Required: file `module:sn` issue "unify spatial-moment-
axis discrimination behind ONE typed predicate" (`SpatialMomentSpace` factor minted by D5b-S3-A0
#240), trigger = `n_diag==2^d`, deadline = S4. The honest fix passes an explicit `is_moment_valued`
flag from the caller-who-knows, or queries the field's space factor once arrays are typed.

## Approval conditions returned
1. DO-NOW: collapse `_frame_signs` → `_moment_frame_signs` (one frame-convention home).
2. RECORD: `module:sn`/`type:improvement` issue for the shape-probe→typed-predicate debt, S4 trigger.

## WHOLE-DELIVERABLE FOLLOW-UP (2026-06-17) — the COMPLETE D5b-S3 reviewed = PASS-WITH-NITS
The full unified-moment-matvec + ERR-061 + OWED-2 scan reviewed together (23-file working tree).
**Verdict PASS-WITH-NITS, zero BLOCK.** The 5 acceptance bars all hold:
- **Bar 1 (branch-freeness):** GENUINE. grep clean (zero `len(s_axes)>1`, zero scheme-isinstance,
  zero `if d==1` dispatch — the one `d==1` hit is a COMMENT asserting absence; `key=` is the registry
  key). Dimensional dispatch is the pure scheme trait `spatial_basis_per_axis` + `face_moment_tail`
  formula + geometry-keyed supports. Curvilinear LD stays a DECLARED raise.
- **Bar 2 (single-sourcing — THE key check):** PROVABLY ONE source each. `_slope_fold` (`_ubld.py:400`)
  is the ONE `s_hat` site; all 3 consumers delegate (schur_xV:462, scan_slope_face_source:522,
  scan_reconstruct→schur_xV:543). NO twin of the fold scan-vs-DAG. `_reframe` defined ONCE
  (sweep_graph:87), imported by loss_rep. `octant_moment_frame_signs` owns the None-for-DD decision.
- **Bar 3 (M-norm legibility):** STRONG PASS. `residual_kernel_batch` docstring :639-649 NAMES the
  convention (M·S⃗ test-fn projection vs raw OperatorSum subtract), `M_ii=θ^{|i|}` named w/ meaning,
  inline `M=⊗diag(1,θ)` at :655. NOT a magic divide. (Σ_s⊗I+M⁻¹ vs Σ_s⊗M = documented-at-site choice,
  investigator ruled correct; not re-litigated.)
- **Bar 4 (retired scalar arm):** `_kernel_terms` DELETED from production (only comment + derivations/
  + symbolic tests reference it). `kernel_rhs` survives as a flat-source primitive VIEW on the live
  D1ClosedForm w/ ZERO production callers, foundation-pinned (`test_divV_kernel_view_equals_dense`) —
  the ÷V flat-source reference-oracle, NOT a retired matvec arm. Acceptable (watch-item, not nit).
- **Bar 5 (reads-like-the-math):** PASS. Named carriers (s_bar/s_hat/mom_sweep/mom_global/QV_chain_sweep);
  φ̂ reduction `einsum("k,kgxp->gxp", w_ords, mom)` reads as Σ_n w_n ψ̂_n; `_CellSolve.__post_init__`
  exactly-ONE-mode guard = textbook Pattern-4.
**DO-NOW from the ERR-061 review (1 above) = RESOLVED:** the `_frame_signs` closure is GONE; d=1 sites
(loss_rep:2246/:2950) call the primitive directly because `_OneDimScanWalk` does NOT inherit
`_LossRepresentation` (can't reach the method). What remains is a 2-line BINDING dup (`per_axis=...;
octant_moment_frame_signs((sign,),per_axis)`), NOT a decision twin (the None-decision lives once in the
primitive). Downgraded ERR-061 DO-NOW → NIT (offer: hoist the binding to a module-level
`frame_signs_for(scheme, signs)` free func both arms call). CONCERN (2 above) = #246 (already filed,
confirmed in brief; `_reframe` guard byte-identical to prior review, did not regress; +2 more inline
`face_moment_tail(...)==()` probes in solver helpers = same family, same #246).
**GATES re-run by me:** LD foundation+frame+primitive 34✓ (after a ONE-TIME non-reproducible flake on
`test_ld_slope_moment_global_frame_consistency` under -O — probed directly: production stores GLOBAL-frame
slope +0.04596 under BOTH -O and non-O, both solvers; 3× re-run of exact command 34✓; transient, NOT a
defect). solve dir under DriftWarning-error 60✓. LD slab MMS 7✓ (diffusion tripwire + scan≡dag + krylov≡si).
**Nits returned (all non-blocking):** (a) d=1 frame-sign binding 2-line dup → optional free-func hoist;
(b) `kernel_rhs` orphan primitive view → keep as oracle, no action; (c) the 1-time flake → if it recurs,
hunt cross-test state leakage (a scheme/quad cache) — log only.
