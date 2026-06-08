---
name: issue-197-typed-cylinder-lessons
description: Consolidated durable lessons from the Issue #197 typed-cylinder sequence (PR-TYPED-2 typed fields, R-1 Step 4b Krylov carve, PR-TYPED-6/6c unified matvec). The gap-fill kernels NOT captured by the committed 6_5/6b/6c-family closeouts. All work LANDED on origin/main.
metadata:
  type: project
---

# Issue #197 typed-cylinder sequence — consolidated durable lessons

All work in this note is **LANDED on origin/main** (BC realizer + typed fields +
unified matvec + SNStreamingOperator retirement all merged). This is the
distillation of the durable lessons from the five uncommitted gap-fill diary
entries (`issue_197_pr_typed_2`, `issue_197_r1_step4b`,
`issue_197_pr_typed_6_situational_analysis`,
`issue_197_pr_typed_6_v3_situational_analysis`,
`issue_197_pr_typed_6c_step6_equivalence_audit`, `issue_197_pr_typed_6c_closeout`)
that are NOT already preserved in the committed
`issue_197_pr_typed_{0,1,3,4,5,6_5,6_foundation,6b}_closeout.md` family. The
bulk of the 6c narrative (B1'' face-state, cylinder twin-path, ERR-049) lives in
the committed `issue_197_pr_typed_6_5_closeout.md`; the Option-β /
`compute_psi_half_per_level` foundation lives in the committed
`issue_197_pr_typed_6b_closeout.md`. This note holds the remaining kernels.

Companions (committed, in this worktree): [[issue_197_pr_typed_6_5_closeout]],
[[issue_197_pr_typed_6b_closeout]], [[issue_197_pr_typed_5_closeout]].

---

## Lesson 1 — the discrete curvilinear WDD/M-M operator is NOT additively separable (the v6 architectural finding)

**The single durable kernel salvaged from the two `*_situational_analysis`
diary entries** (both were point-in-time STOP-and-report memos of an in-flight
state now fully landed via Option β — the analysis stands, the in-flight status
is moot):

The continuous SN operator splits cleanly:
`L_continuous = Ω·∇ψ + (1−μ²)/r ∂ψ/∂μ` (streaming + angular redistribution).
A brief asked for the **discrete** matvec to be **positively constructed** from
three additive peer-leaves `(L_spatial + L_angular + C).apply(ψ)`, explicitly
forbidding the subtractive form.

**Finding: the discrete decomposition is NOT additive at the per-cell-residual
level.** With the WDD spatial closure + Morel-Montry angular closure, the
discrete streaming part and discrete angular-redistribution part **share the DD
coupling constants `c_in`, `c_out`** (computed from `α_in`, `α_out`, `τ`, `ΔA/w`).
Setting `σ_t=0`, `source=0`, `angular_upstream=None` in `cell_update.residual`
does NOT yield pure streaming for sphere/cylinder — the `(ΔA/w)·c_out·ψ̄` term
(the discrete angular-redistribution coupling) survives. The brief's "set σ_t=0
and the residual is streaming-only" is true ONLY for Cartesian where `ΔA/w = 0`.

**Resolution (what landed): Resolution A's subtractive leaf is the principled
form.** `StreamingOperator.apply(ψ) := M(ψ; σ_t) − σ_t⊙ψ`, with `(L + C).apply`
reproducing the full matvec `M(ψ; σ_t)` bit-exactly. The "one algebra, two
consumers" goal (sweep AND matvec route through `cell_balance_for_streaming`,
the Pattern-2 single source of truth) was met WITHOUT peer-leaf separation. The
key architectural identity, verified algebraically:
`cell_update.residual / V_i = (L + C).apply(ψ)[n, i]` — the cell-balance residual
EQUALS the matvec output per (n,i) modulo a per-cell `V_i` factor.

**Durable rule**: when a brief asks to decompose a discrete operator additively,
first check whether the discretization's closure constants couple the
"separate" continuous terms. If they do, the additive peer-leaf split is
unspellable without subtraction; the subtractive `M(ψ) − diag` form (Resolution
A) is the principled resolution, and "one algebra body, shared by both
consumers" is the achievable elegance goal, not "additive peer leaves".

A secondary salvaged kernel: the v3 analysis correctly predicted the
performance trap — a literal per-ordinate `cell_update.residual` call inside the
matvec inner loop would be ~40× slower than the ordinate-vectorized body. Option
β (factor `cell_balance_for_streaming`, vectorize over the ordinate-mask)
preserves both performance AND single-source-of-truth. This is why the landed
implementation vectorizes over `n_mask` inside `StreamingOperator.apply` rather
than looping per ordinate.

---

## Lesson 2 — the face-block linearity anti-pattern (R-1 Step 4b; candidate vv-principles signature)

**The load-bearing discovery of the Krylov carve.** When carving
`SNSolver._solve_krylov` onto the typed `AngularFlux` algebra with a sweep
preconditioner on the B1''-packed flat vector (cells + outer-face + inner-face
slots), GMRES converged to a numerical residual of ~1e-31 but to the **WRONG ψ**
(2G slab k_eff = 1.634 vs expected 1.875 — a 24% drift). The **cell-block**
residual was machine zero; the **face-block** residual was 0.18.

**Root cause**: the sweep's outflow at `boundary.xmax_face` is a NONLINEAR
functional of `rhs` (the WDD outflow at the converged face, not a linear
operator on the face-slot input). Encoding this into the preconditioner's face
block breaks Krylov convergence theory: `M(rhs_1 + rhs_2) ≠ M(rhs_1) + M(rhs_2)`
on the face slots.

**Fix**: identity on the face block — the inverter copies the rhs's face-block
values through unchanged (preconditioner acts as identity on faces); GMRES then
drives the face residual to zero via its own iterations. Converges in ~36
iterations to the correct ψ.

| inverter face block | iters | residual | result |
|---|---|---|---|
| sweep outflow | 9953 | 2.2e-31 | wrong |
| zero | 9931 | 1.1e-31 | wrong |
| **identity on rhs** | **36** | **6.7e-14** | **correct** |

**The diagnostic fingerprint** (worth a `vv-principles` numerical-bug signature):
GMRES converges to wrong ψ with the **cell block at machine precision AND the
face residual non-zero**. Detection: print `(A·ψ_solve − b).boundary` after
convergence — non-zero face residual + machine-zero cell residual ⟹ the
preconditioner is non-linear on the face block. The anti-pattern in one line:
**treating a sweep's output face state as a linear-operator face value.**

---

## Lesson 3 — the SI 4c Carlson previous-iterate hook (why the SI carve was deferred)

The Krylov path carved cleanly (Lesson 2); the **SI carve (4c) was reverted**.
A naive `SourceIteration(LC, S, ZeroOp, inverter=sweep).solve(...)` converged to
a fixed point but the WRONG one (2G slab k_eff = 1.25 vs 1.875; 33% drift),
regressing 6 `test_solver_components` tests.

**Root cause**: the legacy SI threads the **previous iteration's angular flux**
as `initial_guess` to `transport_sweep`, so the curvilinear Carlson coupled-pole
seed receives `σ_t·φ_0(prev)/W` rather than the in-iteration source `Q_p/W`. The
`SourceIteration` primitive's `inverter(rhs)` contract passes only
`rhs = S·ψ + F·ψ + q_ext` — not the previous iterate. The Carlson seed then
falls back to the fresh-start convention, which changes the curvilinear
iteration's fixed point on heterogeneous / multi-group cases.

**Durable rule**: a generic `SourceIteration` primitive that passes only `rhs`
to its inverter cannot reproduce a sweep whose seed depends on the previous
iterate. Either (a) extend the primitive with an SN-specific previous-iterate
hook, or (b) restructure the Carlson seed so the in-iteration-source fallback
gives the same fixed point. The Krylov path does NOT need this hook (it operates
on the full packed vector, not iterate-to-iterate), which is why the carve split:
Krylov landed, SI stayed on the legacy bit-identical loop.

---

## Lesson 4 — per-geometry equivalence audit before retirement (6c Step 6)

Before retiring `SNStreamingOperator`, a per-geometry audit quantified where the
legacy bundle equals the new `(L+C)` algebra at ULP and where it does NOT.
On 4-cell 2G fixtures (`scratch/step6_equivalence_quantify.py`):

| Geometry | rel diff | Status | Cause |
|----------|---------:|--------|-------|
| Sphere | ~1.7e-16 | bit-identical at ULP | (no divergence — same algebra, same code) |
| Slab | ~2.3–3.5 | **structural** divergence | FD (1st-order, legacy) vs WDD (2nd-order, unified) |
| Cylinder | ~0.6–1.0 | **structural** divergence | legacy ERR-049 routing bug; unified is correct |

**The pattern: encode the EXPECTED divergence as strict-xfail guards.** The
divergence was pinned in `TestCompositionEquivalence` with `xfail(strict=True)`:
sphere PASSES (bit-identical), slab+cylinder XFAIL (expected structural gap). If
someone later "fixes" the legacy matvec to match the unified body, the xfail
spontaneously passes → the suite FAILS → flags the change for review. This is
the structural artifact that documents a divergence permanently in the test
suite, distinct from a numeric-tolerance contract. (The guards retired WITH the
bundle at PR-TYPED-7; the record lives here.)

**Sub-lesson (L14)**: the cylinder twin-path bug (`SI` vs `Krylov` reaching
different discrete fixed points at rel ≈ 4e-3) was invisible to every
single-algorithm L1 anchor — each path passed `trajectory_resolvent` at < 3%
rel; only the four-leg standoff (composed matvec ≡ SI sweep ≡ legacy matvec ≡
convergence-under-refinement) exposed the inter-algorithm gap. **Pre-retirement
L14 four-leg standoffs catch what single-algorithm L1 anchors cannot** — now
project verification discipline for any matvec/sweep duality.

---

## Lesson 5 — PR-TYPED-2 typed-field landing decisions (the durable choices)

The three SN flux types (`AngularFlux`, `ScalarFlux`, `BoundaryFlux`) landed
simultaneously (`AngularFlux.integrate_angular() → ScalarFlux` couples them;
`BoundaryFlux.pole_phi_prev: ScalarFlux`). The decisions worth keeping:

- **`BoundaryFlux` is MUTABLE; `AngularFlux`/`ScalarFlux` are FROZEN.** The
  sweep's persistent-BC contract is a write-through cache (reflective-BC partners
  read previous-sweep outgoing-face writes); frozen-and-return-new would force
  per-sweep allocation churn the hot path can't afford. `BoundaryFlux` is a
  STATEFUL accumulator; the value fluxes flow through dunder arithmetic and
  return fresh instances. The mutable/frozen split is principled per the type's
  ROLE (accumulator vs value).
- **Dropped the `Face` enum in favour of named attributes** (`xmin_face`,
  `xmax_face`, `xmin_xmax_buf`, `ymin_ymax_buf`). Decision rule surfaced:
  **enum vs named-attribute is "is the face an INDEX (named attribute) or a
  DISPATCH KEY (enum)?"** Here the 4 face slots are never branching dispatch
  keys, so named attributes win. (Candidate `coding-elegance` Pattern 4
  decision-flow addition.)
- **Operator `apply` accepts `np.ndarray | TypedField → matching`** via a
  4-line isinstance dispatch (typed-in → typed-out; raw-in → raw-out preserves
  the legacy packed-vector / FD-matvec consumers). `SNStreamingOperator.apply`
  KEPT on raw ndarray — the packed vector has BC equations interleaved
  (`EquationMap.ordinate`), so it is NOT an `AngularFlux`.
- **Pattern 4 confirmation in the wild**: after retiring the 7 SNSolver
  `@property` shims (`sig_t`/`sig_a`/.../`_cells_by_mat`), every test site that
  still read `solver.sig_t` surfaced as an immediate `AttributeError` at the test
  boundary — not a wrong number, not a silent buffer drift. Stringly-typed dict
  keys (`psi_bc["bc_xmim"]` typo) had previously lazy-initialised a fresh wrong
  buffer silently; the named-attribute `BoundaryFlux` makes the typo unspellable.

---

## Cross-references to landed code

- BC realizer waves (0,1,5,7,8,9,10,11) + C188 — see the `wave_*` notes in this dir.
- Unified matvec + B1'' face-state — committed `issue_197_pr_typed_6_5_closeout.md`.
- `compute_psi_half_per_level` Option-β foundation — committed `issue_197_pr_typed_6b_closeout.md`.
- ERR-049 (legacy cylinder matvec routing) — `.claude/skills/vv-principles/error_catalog.md`.
- Scattering S = S_foldable + S_residual split — `issue_196_phase_g_step3_4a_closeout.md`.
