---
name: issue-257-s8-streaming-pure-L-verification
description: #257 S8 gate spec — StreamingOperator (L+C)−C fold → pure L (σ-free) + C=M[σ_t] + loss-rep=L+C + apply(FullField) re-point + S/F fibration. PRINCIPLED-equivalence (not bit-id, the recomposition re-associates FP). PRE-IMPL, branch feature/field-typed-operator-algebra.
metadata:
  type: project
---

# #257 S8 — StreamingOperator → pure L (BEHAVIORAL composition carve) — GATE SPEC

PRE-IMPL. Branch `feature/field-typed-operator-algebra` (HEAD `90ea1cd`). The
operator-composition-seam proactive trigger (MUST dispatch test-architect first).
Authoritative plan = `.claude/plans/issue_257_coefficient_field_promotion.md`
(⭐ CURRENT STATUS + §S8). Live anchors RE-CONFIRMED this session (the census
line numbers had drifted — `as_scipy_linop` already RETIRED in S7, not at :1785).

## The system-under-test (4 bundled changes, ONE seam)
1. **Drop `(L+C)−C`**: `StreamingOperator.apply` (`operator.py:419-420`) today =
   `loss_action(σ_t,ψ).bulk − σ_t[None]·ψ.bulk` (fused full loss minus C). S8 →
   pure streaming `Ω·∇ψ`, σ-FREE. σ leaves `StreamingOperator`'s surface
   (`sigma_t` field `:308`); `C` = the shared `MultiplicationOperator[σ_t]`
   (S3b, `transport/multiplication_operator.py:115`; `CollisionOperator:512` is its
   subclass). Loss-rep (`loss_representation.py`) stays `L+C` = the sweep = σ lives
   HERE legitimately (cell optical depth). `InvertibleOperator.apply` (`:839/865`)
   single-sources σ from `self.diagonal.sigma` (= `self.sigma`, `:833-835`) — ALREADY
   does this since #240 Step B.
2. **apply(FullField) re-point**: leaves typed `LinearOperatorMixin["TimedFullField"]`,
   stamp `history_depth=psi.history_depth` on outputs (`:425/469`, MultOp `:212`,
   Fission `:429`). S8 → consume/emit timeless `FullField` (`transport/full_field.py`;
   `TimedFullField=Cofree(FullField,d)`). Driver `advance`s timeless result into
   timed state (`timed_full_field.py:274`).
3. **Helper widening** (`loss_representation.py`): `_require_typed_composite` (1
   helper, 4 sites `operator.py:410/455/864/881`) + ~14 `loss_action`/`_transpose`
   defs + ~15 `psi:"TimedFullField"` param sites. Invariance seam `operator.py:641/794`
   RESOLVES when both ends read `FullField`.
4. **S/F fibration**: `ScatteringOperator` (`scattering.py:286`, `np.ndarray|HarmonicMomentField`)
   + `FissionOperator` (`fission.py:327` `@singledispatchmethod` over {TimedFullField,
   ScalarFlux, ndarray}) branch on carrier TYPE — single-`V` is a partial lie for these two.

## A. Equivalence classification — PRINCIPLED-EQUIVALENT (NOT bit-identical)
**Verdict: principled-equivalence.** The math `(L+C)−C=L` is exact, but the
recomposition changes the FP reduction tree: TODAY `streaming_action(ψ)+σ·ψ` is
computed FUSED inside `loss_action` then `−σ·ψ` subtracted (the σ·ψ appears, cancels,
re-appears); S8 computes pure `streaming_action(ψ)` once + `C.apply=σ⊙ψ` once and ADDS
at the OperatorSum level. Different associativity ⇒ ≤ few-ULP drift on the field-level
matvec. ⭐ PROBE-GROUNDED (the `test_removal_form_matvec_sweep.py` premise correction,
verified this session): the WDD matvec is AFFINE in σ in the forward direction
(`M(σ)ψ=streaming_action(ψ)+σ·ψ`), so `L_pure.apply+C.apply ≡ M(σ_t)ψ` to ≤2 ULP — the
S8 recomposition is the SAME affine algebra, just with `streaming_action` now exposed
as the leaf instead of `loss_action−σ·ψ`.
- Per-output tolerance: field-level `(L+C−S−B)·ψ` matvec → `assert_regression(kind=direct,
  reduction_depth=nx)` (single-step FP-non-assoc, the EXISTING gate convention for these
  snapshots). SI rhs / converged keff & flux → iterative `SAFETY(10)×conv_tol`. Boundary
  trace `L.boundary` → expected STRICT 0-ULP (pure-L still emits the same outflow defect
  `ψ_out=2ψ̄−ψ_in`; C/S/F never touch boundary).
- Gate = vv 3-criteria: (1) named principled intermediate = `streaming_action` (pure
  `Ω·∇ψ`, the σ-free leaf — a genuine domain quantity, not a reduction artefact) +
  `C=M[σ_t]` (the multiplier); (2) struct-indep ref = k_inf closed-form + Q/Σ_t (NOT
  old-vs-new ULP proximity); (3) drift dimensionally explained (`reduction_depth×ULP`).
  ⚠ Do NOT force a 0-ULP snapshot on the bulk — that is the WRONG gate (the S6.3/S7
  lesson; `vv-principles` §bit-id-vs-principled).

## B. Legacy-pinning inventory (what currently pins `(L+C)−C` / σ-on-L / timed codomain)
**BREAK BY DESIGN — re-baseline/migrate (the L20 retirement scope):**
- ⭐ `tests/sn/operators/test_streaming_operator_decomposition.py` — THE file most
  affected. `TestResolutionADecomposition` (`(L+C).apply≡M` bit-exact) STAYS GREEN by
  construction (the composite still = M). But `TestSubtractiveDefinition` (pins
  `L.apply.bulk == M.bulk − σ_t·ψ.bulk` — the `(L+C)−C` subtractive form, `array_equal`)
  is the OLD definition of L → BREAKS BY DESIGN; re-baseline to `L_pure.apply.bulk ==
  streaming_action(ψ).bulk` (= `M − σ·ψ` to nULP, no longer the *defining* equation).
  `TestResolutionADifferentFromPriorWrong` (ERR-058 σ_t=0 coincidence) — RE-EXAMINE: with
  pure-L, σ_t=0 is the SAME pure streaming → the "differs/coincides at σ_t=0" framing
  dissolves (there is no σ on L to set to 0). RETIRE or re-point to the loss-rep.
- ⭐ `test_streaming_operator.py::TestT4bPreT4RegressionSnapshot` +
  `::TestT4cPreT4RegressionSnapshotCurvilinear` — these snapshot `StreamingOperator.apply`
  output (the `(L+C)−C` value). Under S8 `L.apply` returns PURE streaming (≠ the old
  `M−σ·ψ` only by the affine algebra to ULP, BUT the snapshot is of the OLD operator's
  full action). DECISION: the snapshots are of `L`'s apply; pure-L `streaming_action`
  numerically == `M−σ·ψ` to ≤2 ULP (affine), so these can stay as `assert_regression`
  (the kind=direct nULP gate ALREADY absorbs ~1 ULP) — BUT VERIFY: re-run pre-impl and
  confirm the live pure-L value reproduces the frozen snapshot within `reduction_depth×ULP`.
  If the boundary 0-ULP strict gate trips, it's a real bug (boundary must not move).
  ⚠ The 5 SPHERE snapshots are ALREADY baseline-red (#250, NOT ours) — `-k "not
  (sphere_1g_apply or sphere_2g_apply)"`.
- `test_loss_action_convention.py` — STAYS GREEN (it pins the loss-rep `loss_action`
  returns `(L+C)ψ` + the `−C` glue). Under S8 the `−C` glue MOVES from inside `L.apply`
  to the OperatorSum (`L_pure + C`). The structural anchor `test_loss_action_is_full_loss_LpC_flat_reflective`
  (`loss_action(flat)=σ_t·ψ`, `apply(flat)≈0`) — `apply` now = `streaming_action(flat)`
  which is STILL ≈0 on the reflective fundamental (pure streaming annihilates flat).
  RE-POINT `test_apply_equals_loss_action_minus_independent_collision_het` to assert the
  NEW glue: `(L_pure+C).apply == loss_action` (the composite recovers the full loss).
**STAY BYTE-GREEN (must not move):**
- `test_collision_operator.py` (C unchanged — already `M[σ_t]` subclass), `test_multiplication_operator.py`,
  `test_invertible_operator.py` (composite still = M(σ)).
- `test_removal_form_matvec_sweep.py` — its STRUCTURAL teeth (`op.apply==M(σ_r)` array_equal)
  + the σ_C==σ_t bit-id arms STAY (S8 makes them MORE true — `L_pure+C(σ_r)` is exactly
  `streaming_action+σ_r·ψ=M(σ_r)ψ`, no longer relies on the affine-cancellation accident).
- `test_g_adjoint_reciprocity.py` (the curvilinear `Lᵀ` G-adjoint; `−C` glue moves but
  the reciprocity oracle `op.H==G⁻¹opᵀG` is composition-invariant — VERIFY stays green).

## C. New-convention catchers (ADD — the teeth)
- **C1 — pure-L intrinsic σ-freedom (THE defining property of the carve).** Construct
  `L = StreamingOperator(sn_mesh)` (no σ). Mutate σ on the loss-rep / build C(σ_a) vs
  C(σ_b) → `L.apply(ψ)` must be BYTE-IDENTICAL (`array_equal`) regardless of σ. Only
  `C`/the loss-rep see σ. This is the user's intrinsic-property directive + the Mode-2
  catcher (a σ leak back onto L). ⚠ TEETH: mutation-verify a stub that re-reads σ into
  `streaming_action` REDDENS this. Per-geometry (slab/sphere/cyl); ≥2G het σ.
- **C2 — per-ordinate flat-flux residual** (Modes 1/2/4 — sign flip / μ_x↔μ_y swap /
  α-recursion drift in pure-L). On uniform Q, uniform σ_t, reflective box: flat ψ=Q/σ_t
  → `streaming_action(ψ_flat)=0` PER ORDINATE (not just the summed balance — anti-pattern
  #8). Reuse `test_per_ordinate_flat_flux_consistency`'s structure; this is now a DIRECT
  test of the pure-L leaf (it was previously folded inside loss_action). Curvilinear
  MANDATORY (Signature 1; flat-ψ nulls the redistribution → must check per-ordinate).
- **C3 — heterogeneous ≥2G eigenvalue** (Mode 6 — σ single-source convention drift).
  `test_kinf_homogeneous` 2eg/4eg × {slab,sphere,cyl} × {krylov,SI} — STAYS GREEN
  (k_inf closed-form, struct-indep). ADD a HETEROGENEOUS 2G k_eff/standoff arm
  (`test_l1_standoff_slab_cylinder.py` is the home) — homogeneous nulls redistribution
  (#4); 1G is degenerate (#3, Cardinal Rule). The σ that feeds C vs the σ in the
  loss-rep sweep MUST be the SAME single source — a drift shows as wrong group ratio.
- **C4 — Mode-9 FP-invariance (the recomposition-graph trap).** Converged flux ==
  pre-S8 converged flux to solver tol, on a config that BREAKS the degenerate
  coincidence: ANISOTROPIC flux (VACUUM / heterogeneous / streaming — NOT the
  fully-reflective isotropic box, which would be blind to a recomposition bug) AND a
  DIAGONAL cubature (`level_symmetric`/`lebedev`) for the 2-D/curvilinear arm (the
  shared-face ERR-056 hazard). Because the composition graph CHANGES (L+C now adds two
  leaves vs the fused subtract), an isotropic-reflective-only FP check is BLIND. Use
  `test_fixed_source_2d_equivalence.py`'s Leg-1 Q/Σ_t (≥2G, vacuum) + a curvilinear
  het arm. Assert against the LIVE pre-S8 value captured in-process OR the Q/Σ_t
  closed form (the struct-indep ground), not a hardcoded baseline.
- **C5 — apply(FullField) timeless round-trip.** (a) Operator OUTPUT carries no history:
  `L.apply(FullField)` returns a `FullField` (not `TimedFullField`) — `isinstance` +
  `not hasattr history`/`history_length==0`. (b) Driver re-attach: the SI/Krylov driver
  `advance`s the timeless result into the timed state; converged TIMED state unchanged
  vs pre-S8 (the iteration sees the comonad, the operator does not). (c) `advance`
  type-guard still fires (`timed_full_field.py:306-317`). Per-geometry; assert the
  converged `.scalar_flux` == pre-S8 to `SAFETY×conv_tol`.
- **C6 — fibration dispatch parity after re-typing.** `ScatteringOperator.apply` +
  `FissionOperator.apply` (`@singledispatchmethod`) still dispatch correctly across
  {FullField/TimedFullField, ScalarFlux, ndarray} after the carrier re-point. Pin each
  arm returns the right type + value (`test_operators_apply_typed.py` already covers
  the typed lifts — EXTEND for FullField input; assert `F.apply(FullField)` and
  `F.apply(ScalarFlux)` agree on the bulk via the shared ScalarFlux branch, the
  single-source-of-truth `fission.py:413`).

## D. L1/MMS structural backstop (must STAY GREEN — the real correctness ground)
The composition change must NOT move the converged-to limit (anti-pattern #5). Exact paths:
- `tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py` + `test_curvilinear_aniso_convergence.py`
  (`catches("ERR-026")`) — the curvilinear streaming closure MMS, the load-bearing
  redistribution backstop. ⭐ The streaming carve is closest to THIS math.
- `test_mms_ld_slab.py`, `test_mms_ld_2d.py`, `test_mms.py`, `test_mms_aniso.py`,
  `test_mms_curvilinear.py`, `test_mms_heterogeneous.py` — slab/2-D LD + aniso MMS (O(h²) order).
- `test_kinf_homogeneous.py` + `test_kinf_homogeneous_tolerance.py` — closed-form k_inf
  (the eigenvalue ground; MMS does NOT prove eigenvalues).
- `test_si_convergence_rate.py` — the RATE gate (`history.n_inner`, ρ_J=c scattering
  ratio). S8 is value-preserving so the rate must be UNCHANGED; this gate stays green.
- `test_l1_standoff_slab_cylinder.py`, `test_phase_c_crosscheck.py` (trajectory_resolvent
  struct-indep), `test_prescribed_inflow_consistency.py`.

## E. Mode-11 (the gate must EXECUTE the rewired path)
- C1 → the NEW pure-L `StreamingOperator.apply` (`streaming_action`, σ-free). Sentinel:
  file-write in the new pure-L body; confirm C1's run fires it. ⚠ `StreamingOperator.apply`
  has ZERO graph callers (Nexus confirmed) — it's reached only via `OperatorSum`/the
  driver. So a gate that calls `L.apply(ψ)` DIRECTLY (C1, C2) is REQUIRED — a solve-only
  gate (C4, D) routes through `(L+C).solve`=the SWEEP, which is `loss_representation` and
  does NOT call `L.apply` (the matvec leaf). The SWEEP path and the MATVEC path are
  SEPARATE production lines — C4/D pin the sweep; C1/C2/C-snapshot pin the matvec leaf.
- C5 → the timeless-codomain emit + driver `advance`. Sentinel: confirm the operator
  output `type is FullField` (not TimedFullField) AND the driver's `advance` call fires
  (counter on `advance`).
- C6 → each `@singledispatchmethod` arm. Sentinel: dispatch-table hit per type.
- ⚠ The legacy SNAPSHOT gates (`TestT4b/c`) — VERIFY they still EXECUTE `StreamingOperator.apply`
  after the re-typing (they build `L` + call `L.apply(state)` directly — they DO reach
  the matvec leaf; this is the one place the matvec leaf is exercised by a regression gate).

## F. Sequencing — ORDERED SUB-SEQUENCE, not one atomic carve
RECOMMEND landing as an ordered sequence, each independently gated:
1. **S8a — apply(FullField) re-point + helper widening** (mechanical, behavioral-NEUTRAL
   to VALUES). Gate: bit-identical converged state (`advance` decouple is FP-neutral);
   C5 + the FullField discriminating tests. LAND FIRST (it's the typing seam; the
   `(L+C)−C` drop in S8b is cleaner once codomain is FullField).
2. **S8b — drop `(L+C)−C` → pure-L + C=M[σ_t]** (THE behavioral composition change,
   principled-equiv). Gate: C1/C2/C3/C4 + re-baseline B's snapshots + the D backstop.
3. **S8c — S/F fibration resolve** (`@singledispatchmethod` re-type). Gate: C6.
⭐ RISKIEST = S8b (the recomposition graph change — Mode-9 trap, the FP-tree re-association).
S8a is the second-riskiest (the timeless decouple touches every operator + the driver;
a botched `advance` silently drops history → C5's round-trip catches it). Splitting lets
S8a land bit-identical (a clean checkpoint) before the value-moving S8b.

## G. Failure-mode exposure (most-exposed → catcher)
- **Mode 6 (convention drift — σ single-source)** MOST exposed: σ now lives in TWO
  conceptual homes (C=M[σ_t] surface vs loss-rep sweep cell-depth) → a drift between
  them is the headline risk. Catcher: C3 (2G het k_eff/standoff) + C1 (σ-free L proves
  L does NOT read σ).
- **Mode 9 (splitting/degenerate-regime)**: the recomposition graph change. Catcher: C4
  (anisotropic vacuum/het + diagonal cubature FP-invariance).
- **Mode 2 (variable swap / σ leak onto L)**: C1 (intrinsic σ-freedom) + C2 (per-ordinate).
- **Mode 11 (gate-never-executes-rewired-path)**: §E — C1/C2 call `L.apply` DIRECTLY
  (the matvec leaf has zero callers + the sweep path routes around it).
- **Mode 1/4 (sign/recursion in pure streaming)**: C2 per-ordinate flat-flux (curvilinear).
- NO new failure mode (no vv table row). NO new ERR until a real bug (next free ERR-063
  per the LD cluster — VERIFY at impl time, may have advanced).

## Route-around reds (baseline, NOT ours — every gate)
`-k "not (sphere_1g_apply or sphere_2g_apply)"` (#250 5 stale SPHERE snapshots) + the
#232 mu_y 2 + `--deselect …test_keff_slab::test_heterogeneous_absolute_keff` (#212 hang).
NEVER run all `tests/sn` (#212). 7 baseline reds total; every stage keeps exactly these.

Extends [[issue-257-s3-multiplication-operator-verification]] (C=M[σ_t] promotion) +
[[issue-257-s6-integral-kernel-category-verification]] (S/F as typed Kernels) +
[[issue-206-phase-c-verification]] (the matvec-leaf-vs-sweep separation, DriftWarning
escalation) + [[feedback-regression-tolerance-design]] (direct/iterative tol) +
[[feedback-vv-tagging]] (foundation no-verifies) + [[feedback-principled-over-bit-identical]].
