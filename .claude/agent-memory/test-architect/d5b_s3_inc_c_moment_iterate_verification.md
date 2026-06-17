---
name: d5b-s3-inc-c-moment-iterate-verification
description: Verification spec for #240 D5b-S3 — folding Increment C (the spatial-moment iterate φ̂ + the scattering-slope source Σ_s·φ̂ + the genuine (2^d,ng) moment-source carrier + closing the d≥2 Krylov 2-D LD matvec raise) into the SN solver. THE LOAD-BEARING RESOLUTION: Inc C is a PHYSICS-COMPLETION, not a splitting/accelerator — the converged FP CHANGES vs S2 (S2 converges to (L+C−S_flat)ψ=Q which drops the slope-scattering source = diffusion-limit-INCONSISTENT; S3 converges to (L+C−S_full)ψ=Q = the diffusion-limit-CONSISTENT operator). So the durable plan's "Mode-9 FP-invariance vs the S2 flat-source FP" wording is WRONG and would ship a correctness bug behind a green gate. The genuine Mode-9 invariant is SI-with-lagged-moment-iterate ≡ direct/Krylov solve of the SAME full UBLD operator (L+C−S_full). PRIMARY correctness gates = thick-diffusion-limit (1-D xfail→PASS + 2-D analog) + transport MMS still O(h²). Branch feature/sn-space-angle-tier2, HEAD f45a219. Extends d5-nd-polymorphism-verification §D5b.
metadata:
  type: project
---

# #240 D5b-S3 — Inc C moment-iterate fold: verification spec (L17, PRE-IMPL)

**Status:** PRE-IMPLEMENTATION. Branch `feature/sn-space-angle-tier2`,
HEAD `f45a219` (S2 done + committed `495af60`). Host env, canonical
`python -O -m pytest`. NEVER run all `tests/sn` (#212 hang).

This memo is the S3-specific spec. It EXTENDS
[[d5-nd-polymorphism-verification]] §D5b (the S2/D5b.0–.5 gates) — the
gates here are NEW (the S2 gates pin S2's flat-source closure; S3
changes the closure, so S3 owns its own correctness chain). Read the
S2 crosswalk `.claude/plans/issue_240_d5b_s2_crosswalk.md` for the
moment ordering `[bar, ŷ, x̂, x̂y]` and the Q-lift S3 replaces.

---

# 0. THE LOAD-BEARING RESOLUTION — physics-completion, NOT iteration-only

**VERDICT: the brief's analysis is CORRECT. Inc C is a
PHYSICS-COMPLETION; the converged fixed point CHANGES vs S2.** The
durable plan's S3 gate wording — *"the moment-source splitting MUST
NOT move the converged FP ... the un-accelerated flat-source FP to
solver tol"* — is **IMPRECISE AND WRONG as a gate**, and writing it
would ship a correctness bug behind a green gate. Here is the precise
operator statement, probe-grounded against the live code:

## What S2 converges to vs what S3 converges to

The scattering source `S·φ` is computed **outside** the sweep (in the
SI loop as the lagged gain `S`, or as the Krylov matvec's `−S` block),
then fed INTO the sweep as `Q_cells`. The sweep kernel inverts
`(L+C)`; the iteration folds in `S`. So the discrete operator the
iteration converges to is `(L + C − S)ψ = Q_ext`, and the question
"what does S2 vs S3 converge to" is entirely "what is the discrete
`S` they use".

- **S2 (today, `495af60`):** `ScatteringOperator.apply` reads ONLY the
  cell-average scalar flux `φ̄` (the ℓ=0 *spatial-average* moment — see
  `scattering.py:1093` `HarmonicMomentField` arm and `:1067`
  `AngularFlux` arm, both reduce to a scalar-per-cell `φ`). The kernel
  then **internally lifts** that scalar into slot 0 (`AVERAGE_MOMENT`)
  of the `2^d` moment-source vector and leaves the slope rows
  (`Ŝ_x, Ŝ_y, Ŝ_xy`) **ZERO** (S2 crosswalk §"Q_cells carrier"). Call
  this operator `S_flat`: it is `Σ_s ⊗ e₀e₀ᵀ` over the spatial-moment
  axis (scatter only the average moment; zero into the slopes).

  ⇒ **S2 converges to `(L + C − S_flat)ψ = Q_ext`.** This is O(h²)
  (the LD streaming/collision discretization is second-order) but
  **NOT diffusion-limit-consistent** — the literature memo
  (`multi_d_ld_closure.md:316`) states it verbatim: *"the flat-source
  multi-D UBLD is O(h²) but NOT diffusion-limit-consistent until the
  moment source iterate is threaded."* The 1-D xfail tripwire's
  `reason=` says the same: *"flat-source LD (Q̂=0) loses the diffusion
  limit; needs the canonical slope source Σ_s·φ̂."*

- **S3 (this step):** `ScatteringOperator.apply` reads `φ̂` too (the
  per-axis *spatial slopes* of the scalar flux) and produces a genuine
  `(2^d, ng)` moment source whose slope rows carry `Σ_s·φ̂`. Call this
  `S_full`: `Σ_s ⊗ I_{spatial-moment}` (scatter every spatial moment,
  not just the average). DD/Step at `n=1` (per_axis=1) reduce to
  `Σ_s ⊗ I_1` ≡ `S_flat` ≡ the existing flat scattering — **so S_full
  ≡ S_flat for DD/Step by construction** (the negative-control
  bit-identity, gate 4).

  ⇒ **S3 converges to `(L + C − S_full)ψ = Q_ext`.** This is the
  **diffusion-limit-CONSISTENT** operator. It is a DIFFERENT discrete
  operator from S2's for the multi-moment schemes (LD-d1 and LD-d≥2),
  and its converged solution is DIFFERENT — *that is the entire point*:
  it is what makes the `ld-thick-diffusive` tripwire flip xfail→PASS,
  *because the FP moves to the correct (diffusion) one*.

## Why "S3-FP == S2-flat-source-FP" is the WRONG gate

`S_flat ≠ S_full` whenever any spatial slope moment of `φ` couples
through scattering — i.e. whenever the converged flux has a non-zero
per-cell slope AND `Σ_s ≠ 0`. That is the GENERIC scattering case (any
heterogeneous / thick / non-flat problem). Asserting `S3-FP == S2-FP`
to solver tol would:
1. PASS on the degenerate regimes where the slope-scattering source is
   ~0 (optically thin, weak scatter, near-flat flux) — exactly the
   regimes that hide the bug (vv H2 — flat flux nulls redistribution);
2. FAIL (correctly!) on the thick-diffusive regime — but the plan
   wording frames that failure as a BUG, so an implementer chasing a
   green gate would "fix" it by NOT threading the slope source, i.e.
   by un-doing Inc C. **The gate would actively reward shipping the
   un-completed physics.**

This is a Mode-9 *mis-application*: Mode-9 (verify a SPLITTING/
acceleration to the same FP) presupposes the change is FP-preserving
by definition. Inc C is NOT a splitting — it changes the operator. The
plan conflated "thread the moment source" (a physics change) with
"split the iteration" (an FP-preserving change). They are different
categories.

## The GENUINE Mode-9 invariant (the brief's reframing — CONFIRMED)

The real iteration-only invariant — the legitimate Mode-9 gate, and
the natural sibling of D5b.4 — is:

> **SI-with-lagged-moment-iterate ≡ direct/Krylov solve of the SAME
> full UBLD operator `(L + C − S_full)ψ = Q`, to solver tol.**

i.e. lagging `S_full·φ̂` across sweeps (SI) converges to the SAME flux
as the simultaneous (Krylov) solve of the full operator. This verifies
the *iteration scheme* (the lag) does not change the operator's true
solution — only the rate. It is the within-group analog of "SI ≡
Krylov on the same operator" (D5b.4, which today is the kernel-level
matvec twin; S3 LIFTS it to the full `(L+C−S_full)` operator because
S3 closes the d≥2 Krylov matvec raise). This gate is gate 3 below.

**Both S2-flat-source-FP and S3-full-source-FP are LEGITIMATE fixed
points of DIFFERENT operators.** The verification job is NOT to make
them equal — it is to (a) prove S3's operator is the CORRECT one
(thick-diffusion limit, gate 1), (b) prove the threading didn't break
the existing O(h²) transport consistency (gate 2), and (c) prove the
iteration (SI lag) finds the SAME solution as the direct solve of S3's
operator (gate 3, the genuine Mode-9).

---

# 1. GATE LIST (PRIMARY correctness gates clearly separated from the FP-invariance gate)

Config discipline applied throughout (vv Cardinal Rule 6 + H1/H2):
≥2G + heterogeneous + non-flat for every scattering claim; NON-SQUARE
2-D (x↔y swap defence); `level_symmetric` quad in 2-D (genuine `mu_y`,
#214-safe); Mode-8 `-O`-safe (`np.testing.assert_*` / `pytest.fail` /
`assert_regression`, NEVER bare `assert`).

## ════ PRIMARY CORRECTNESS GATES (these PROVE Inc C is right) ════

### GATE 1 — THE thick-diffusion limit (the headline correctness gate)

This is the single gate that proves Inc C is correct, not just runs.
It has THREE legs (1-D flip, 2-D analog, basis-discrimination).

**Reference (the structurally-independent ground):** the interior
**diffusion solution** in the asymptotically-thick scattering-dominated
limit (Larsen-Morel-Miller 1987 = LMM-1987, the four-limits criterion;
Adams-2001 the multi-D quadrilateral verdict; BLA-1992 the 2-D
rectangular asymptotic analysis). The diffusion solution is NOT
another ORPHEUS SN solver — it is the analytic limit of the *scaled
transport equation* as `ε → 0` (`σ_t = O(1/ε)`, `σ_a = O(ε)`,
`Q = O(ε)`, `c = σ_s/σ_t → 1`). **This is structurally independent of
the LD kernel** (it is the continuous diffusion PDE limit, not any
discrete operator).

Two acceptable reference realizations (pick per leg):
- **(a) DD-as-anchor (1-D leg, cheap):** DD IS interior-diffusion-
  consistent in 1-D (LMM-1987 Table I) — the EXISTING 1-D tripwire
  uses `LD ≈ DD interior` as the operational reference. This is
  defensible BECAUSE DD's diffusion-limit consistency is itself a
  published, separately-verified property (NOT a circular SN-vs-SN
  comparison — it is "compare against a discretization PROVEN to hit
  the diffusion limit"). Document the chain: LD→diffusion via DD-as-a-
  proven-diffusion-limit-proxy. ⚠ This is the WEAKEST link of gate 1
  (it is a discretization proxy, not the continuous limit) — leg (b)
  is the stronger anchor; keep BOTH.
- **(b) continuous-diffusion-solution-as-anchor (the strong leg):** for
  a homogeneous thick slab with the scaled XS + a known source, the
  interior diffusion solution `φ_diff(x)` is closed-form (a
  cosh/quadratic profile depending on BC). Assert the LD interior
  flux → `φ_diff` within `O(ε)` + the discretization error. This is
  the structurally-independent continuous reference. RECOMMEND for the
  1-D leg's value anchor (not just `LD ≈ DD`).

**Leg 1a — 1-D tripwire flips xfail→PASS.**
- File: `tests/sn/verification/mms/test_mms_ld_slab.py::test_ld_thick_diffusive_limit_xfail` (`:235`).
- Action: it is `@pytest.mark.xfail(strict=False, reason="#158 Increment C ...")`.
  When S3 lands, **flip `strict=False` → `strict=True`** (or remove the
  xfail and rename `test_ld_thick_diffusive_limit`). Per
  [[feedback_vv_tagging]] an xfail that NOW HOLDS but is left
  `strict=False` is a silent false-pass habitat — the gate must
  ASSERT the pass once Inc C lands. The body already asserts
  `rel < 0.05` (LD midpoint ≈ DD midpoint at `σ_t·h=10, c=0.99`).
  ⭐ **DECISION FLAG for the implementer:** removing the xfail vs
  flipping to strict-True. RECOMMEND: rename to
  `test_ld_thick_diffusive_limit` (drop the `_xfail` suffix), remove
  the `xfail` marker entirely, KEEP the `rel < 0.05` assert. A
  strict-xpass that you then have to remember to un-strict is churn;
  if Inc C makes it pass, it is a PASS, full stop.
  ⚠ **1G CAVEAT:** this test is 1-GROUP (`_make_1g_mixture`). Per
  Cardinal Rule 6, a 1G test is degenerate for *eigenvalue* claims —
  but this is a FIXED-SOURCE flux-SHAPE claim (the diffusion-limit
  flux profile), and the diffusion limit is a real flux-shape
  phenomenon visible in 1G (the slope-scattering source `Σ_s·φ̂` is
  active in 1G — `Σ_s` self-scatter `g→g` couples the slope). So 1G is
  ACCEPTABLE for the *existence* of the limit, BUT — per gate 5 below —
  ADD a 2G-heterogeneous companion (leg 1c) so the group-coupled slope
  source (the `Σ_s^T` convention on the slope rows, Mode #6) is
  exercised. Do NOT ship gate 1 on 1G alone.

**Leg 1b — the 2-D Cartesian analog (NEW).**
- File: NEW `tests/sn/verification/mms/test_mms_ld_2d.py::test_ld_2d_thick_diffusive_limit`
  (the 2-D home already exists from S2).
- Config: a 2-D Cartesian homogeneous thick slab/box, `σ_t·h ≈ 10`/cell,
  `c → 1` (`c=0.99`), uniform isotropic source, vacuum edges. The
  reference: `LD-2D interior ≈ DD-2D interior` (DD-2D is also
  diffusion-limit-consistent on Cartesian cells — it is the WDD diamond,
  which preserves the limit) AND/OR the continuous 2-D diffusion
  solution. ⭐ **The basis-discrimination is the load-bearing 2-D
  content:** the UBLD (`2^d` with the `xy` cross moment) PASSES; the
  naïve simplex `1+d` FAILS (Adams-2001) — see leg 1c.
- Reference realization: RECOMMEND `LD-2D ≈ DD-2D` (cheap, DD-2D is the
  proven proxy) PAIRED with a coarse continuous-diffusion sanity band.
  The tolerance `rel < 0.05` mirrors the 1-D leg (the limit is reached
  to O(ε) + O(h²); `σ_t·h=10` puts ε≈0.1, so 5% is the right band —
  document the `ε` budget).
- ⚠ **2G companion:** ALSO run a 2G-heterogeneous thick config (leg 1c
  below subsumes this) so the gate is not 1G-degenerate.
- Markers: `@pytest.mark.l1 @pytest.mark.slow @pytest.mark.catches("ERR-NNN")`
  IF a real basis-choice bug is caught (see Mode-7/basis hazard §3) —
  otherwise `@l1 @slow` + a `verifies("ld-cartesian-2d")` once D6 mints
  the label. NO `verifies` until D6 mints it; link `transport-cartesian-2d`
  (EXISTS) so it is not a total orphan.

**Leg 1c — the basis-discrimination (the silent-physics-bug catcher).**
- Claim: the UBLD (`2^d`, with `xy`) PASSES the thick-diffusion limit;
  the simplex `1+d` (no `xy`) FAILS it (Adams-2001). Since S2 already
  built the UBLD `2^d` object (not the simplex), this is a
  CHARACTERIZATION/regression guard: a future "optimization" that
  drops the `xy` cross moment (to save a moment) would silently break
  the diffusion limit while passing every non-diffusion test.
- File: `test_mms_ld_2d.py::test_ld_2d_thick_diffusive_needs_cross_moment`
  (or fold the rationale into leg 1b's docstring + a `catches` if a bug
  is found). RECOMMEND: a docstring note in leg 1b that the PASS is
  basis-dependent + the Adams-2001 citation; a standalone discrimination
  test is OPTIONAL (the `2^d` object is already built — there is no
  simplex code path to discriminate against unless someone adds one).
  ⭐ If the implementer is tempted to add a "fast" simplex path, leg 1c
  becomes MANDATORY. Flag: track as a follow-up note, not a blocking
  gate, until a simplex path exists to discriminate.
- Config: 2G het, NON-SQUARE, thick. The catcher is "drop the `xy`
  moment row from the source → the limit fails" (mutation test).

**Failure modes caught:** the silent diffusion-limit-inconsistency
(the WHOLE POINT of Inc C); Mode #3 missing factor (the slope-source
`Σ_s·φ̂` magnitude); the basis-choice trap (simplex vs bilinear).

### GATE 2 — transport-regime MMS still O(h²) with the FULL moment source threaded

**Claim:** convergence-order (math, L1) + flux-shape. Threading the
genuine `(2^d, ng)` moment source (replacing S2's internal scalar
lift) must NOT break the existing transport-regime O(h²) — and the
1-D LD slab MMS + the 2-D smoke must still hold.

**What stays / what changes:**
- **1-D LD slab MMS** (`test_mms_ld_slab.py::test_sn_1d_slab_ld_mms_converges_second_order`, `:82`):
  STAYS GREEN. The manufactured 1-D case is a fixed-source MMS in a
  transport (non-diffusion) regime; threading the moment source must
  not perturb its O(h²) (the source IS the MMS-derived source — now
  potentially with a non-zero slope row). ⚠ **SUBTLETY:** the 1-D MMS
  `external_source` builds a `(2, ng)` source today (`_ld_source_moments`).
  S2 threads it through. If S3 changes how the scattering slope feeds
  (the scattering source's slope row now non-zero), the WITHIN-GROUP
  scattering of the MMS case changes — but the MMS `Q_ext` is derived
  to make `ψ_chosen` exact REGARDLESS (the total RHS `Σ_s^T φ + Q_ext`
  reproduces `ψ_chosen`). So as long as the MMS derives `Q_ext` against
  the SAME within-group transport PDE the solver implements, O(h²)
  holds. RE-RUN to confirm; if it drifts, the scattering-slope wiring
  disagrees with the MMS's PDE (a real bug — investigate, don't relax).
- **2-D LD smoke** (`test_mms_ld_2d.py::test_ld_2d_converges_second_order_smoke`, `:90`):
  STAYS GREEN (its ansatz vanishes on the edges, vacuum auto). Same
  argument. The smoke uses the EXISTING isotropic het MMS
  (`build_2d_cartesian_heterogeneous_mms_case`) which is angularly flat
  → its scalar flux has a spatial slope per cell, so the slope-
  scattering source `Σ_s·φ̂` IS now active in the smoke (S3 changes the
  within-group source). RE-RUN; O(h²) must hold (the MMS source absorbs
  the change). If it drifts → bug.

**Claims S3 UNLOCKS vs which STAY S4 (the brief's explicit ask):**
- **S3 UNLOCKS:** (i) the d≥2 Krylov 2-D LD matvec (the
  `_CellResidual.cell` raise at `sweep_graph.py:929` closes — the
  matvec now has the `2^d` moment probe it needs); (ii) the
  thick-diffusion limit (gate 1); (iii) the genuine moment-source
  CARRIER in `_SolveOperands.Q` (the `(2^d, ng)` source replaces the
  internal scalar lift). With the moment source threaded, the 2-D MMS
  CAN now run with a `Q̂≠0` slope source per axis — meaning a
  STRENGTHENED 2-D MMS could in principle make a flux-shape claim.
- **STAYS S4:** the `@verifies("ld-cartesian-2d")` flux-SHAPE claim
  with the **strengthened Mode-7 stress ansatz** (`SN2DCartesianLDStressMMSCase`,
  `μ_x·B + μ_y·C` slope drivers, a0>0 non-vanishing boundary) + the
  **non-vanishing-boundary `BoundaryFlux` moment trace**. S3 unlocks
  the *capability* (Q̂≠0 moment source threads), but the FULL flux-shape
  verification needs (a) the strengthened ansatz (which S3 does not
  mint) AND (b) the non-vanishing boundary trace (which needs the
  `BoundaryFlux` widening S4 owns — the S2 crosswalk §"Domain-edge
  inflow" defers any non-zero domain inflow to S4). **So S3's MMS gate
  stays a CONVERGENCE smoke + the thick-diffusion VALUE anchor; the
  flux-shape `verifies` claim is S4.** ⭐ State this in the S3 PR: "S3
  threads the moment source and unlocks the limit + the Krylov matvec;
  the flux-shape verification (strengthened ansatz + non-vanishing
  boundary trace) is S4."

**Q̂ posture for S3 (resolving the durable plan's open question
`§"The Q̂≠0 slope-SOURCE posture"`):** S3 DOES thread the moment source
(that is Inc C). So the durable-plan branch "if D5b threads the moment
source, the MMS MUST supply Q̂≠0" is the active one. BUT the
strengthened-ansatz MMS that supplies Q̂≠0 per axis is S4 — so at S3
the moment-source CARRIER exists and is exercised by the
SCATTERING-slope path (`Σ_s·φ̂`, internally generated), NOT by an
external MMS Q̂. The S3 gate for the carrier is gate 3 (FP-invariance
on the full operator, which exercises `Σ_s·φ̂` slope rows) + gate 1
(the limit, which IS the slope-scattering source doing its job). The
external-Q̂ MMS slope-SOURCE sign trap stays S4.

### ════ THE FP-INVARIANCE GATE (the genuine Mode-9 — CLEARLY SEPARATE) ════

### GATE 3 — SI-with-lagged-moment-iterate ≡ direct/Krylov solve of the SAME full UBLD operator

**This is the ONLY legitimate FP-invariance gate. It does NOT assert
S3-FP == S2-FP (that would be wrong, §0). It asserts the ITERATION
SCHEME (lagging the moment source) converges to the SAME flux as the
simultaneous solve of the identical operator `(L+C−S_full)`.**

**Claim:** software invariant (foundation). On a NON-DEGENERATE config,
the within-group SI (which lags `S_full·φ̂` across sweeps) and the
Krylov solve (which solves `(L+C−S_full)ψ=Q` simultaneously via the
now-closeable d≥2 LD matvec) converge to the SAME scalar flux to
solver tol. This is where the closed `_CellResidual.cell` raise
(`sweep_graph.py:929`) gets exercised end-to-end.

- File: `tests/sn/verification/mms/test_mms_ld_2d.py::test_ld_2d_krylov_matches_si_full_operator`
  (mirror the 1-D `test_sn_1d_slab_ld_mms_krylov_matches_si` `:123`,
  but now WITH a non-zero `scattering_order=0` P0 self-scatter so the
  slope-scattering source is active — see config).
- **Config (Mode-9 degeneracy-break — MANDATORY, vv Mode-9):** ANISOTROPIC
  flux (NOT the fully-reflective isotropic box) + a DIAGONAL cubature.
  Concretely:
  * **Geometry/BC:** 2-D Cartesian, NON-SQUARE (`nx≠ny`), VACUUM edges
    (vacuum makes the flux streaming/anisotropic — the isotropic-
    reflective box is FORBIDDEN per Mode-9). Het `Σ_t(x,y)`.
  * **Quadrature:** `level_symmetric` (a DIAGONAL cubature — NOT an
    axis-aligned product quad; shared faces exercise the multi-D
    coupling; #214-safe genuine `mu_y`). This is the Mode-9
    requirement for "angular-schedule changes" — and the matvec's
    face reconstruction per axis is exactly such a schedule.
  * **Groups/scatter:** ≥2G ASYMMETRIC downscatter (vv H1; the
    group-coupled `Σ_s^T` on the slope rows, Mode #6) + a NON-ZERO P0
    within-group self-scatter `Σ_s0^{g→g}` (so `S_full`'s slope rows
    are GENUINELY non-zero — a zero-scatter config makes `S_full ≡
    S_flat` and the gate degenerate-blind, the exact Mode-9 trap). ⭐
    **The self-scatter must be strong enough that `Σ_s·φ̂` is a
    non-negligible fraction of the slope source** — otherwise SI and
    Krylov agree trivially because the slope source ≈ 0. RECOMMEND
    `c = Σ_s/Σ_t ≈ 0.5–0.8` (moderate scatter — not the thick limit,
    that is gate 1) so the slope-scattering feedback is active but the
    iteration still converges in a sane count.
  * **Flux non-flatness:** non-flat per-ordinate / het source so the
    per-cell spatial slope `φ̂` is genuinely non-zero (flat `φ̂` nulls
    `Σ_s·φ̂` — vv H2).
- **Tolerance:** `rtol = SAFETY × inner_tol ≈ 1e-9, atol=1e-11`
  (iterative → SAFETY×conv_tol, [[feedback_regression_tolerance_design]]).
  NOT `array_equal` (SI and Krylov are different reduction trees +
  different iteration counts).
- **⚠ The matvec-raise-closure is the SUT here.** S3 closes
  `sweep_graph.py:929`. This gate is the END-TO-END exercise of that
  closure (the kernel-level matvec twin is D5b.4 at S2; this lifts it
  to the full `(L+C−S_full)` operator with the moment iterate). A
  sweep without a verified matvec is half a carve (L14 leg-3 standoff).
- **PAIR with gate 1** (gate 3 proves SI ≡ Krylov on the SAME
  operator; gate 1 proves that operator is the CORRECT one — gate 3
  alone is necessary-not-sufficient: both SI and Krylov could solve
  the SAME WRONG operator and agree).
- Markers: `@pytest.mark.foundation` (a software invariant — SI ≡
  Krylov on one operator; no equation claim).
- **Failure modes caught:** #6 convention drift (`Lψ` vs `(L+C−S)ψ` in
  the matvec); #4 wrong recursion in the apply-direction `2^d`-moment
  face reconstruction; the Mode-9 degenerate-blindness (broken by the
  anisotropic + diagonal-cubature + non-zero-scatter config).

## ════ NEGATIVE CONTROL + DEGENERACY GUARD ════

### GATE 4 — DD/Step bit-identity (the negative control)

**Claim:** at `n=1` (per_axis=1), `S_full = Σ_s ⊗ I_1 ≡ S_flat ≡ the
existing flat scattering`. DD/Step are byte-identical pre/post S3 — the
strict DriftWarning gate stays at the S2 baseline.

- **Exact invocation:**
  ```
  python -O -m pytest tests/sn/sweep/core tests/sn/solve \
    -W "error::tests.sn.regression._regression_assert.DriftWarning"
  ```
  MUST stay at the **S2 baseline** (the S2 record says `562 / 2 skip /
  4 xfail`; the prior-memo DD strict gate was `505/1/4` — **RE-CONFIRM
  the live count at pickup**; the gate is "unchanged vs the immediately
  preceding HEAD `f45a219`", not a hard-coded number). The
  `tests.sn.regression...` path is load-bearing (`orpheus.sn...`
  silently fails to escalate — FINDING-1, [[issue_206_phase_c_verification]]).
- **Why it works:** S3's scattering change is gated on the
  spatial-moment axis. DD/Step have `spatial_basis_per_axis=1` → the
  `_CellSolve.cell` emit gate `len(s_axes) > 1 and spatial_basis_per_axis > 1`
  (`sweep_graph.py:883`) stays FALSE → no trailing moment axis → the
  scattering source stays a scalar-per-cell → `ScatteringOperator.apply`
  on a scalar/average-only flux is UNCHANGED. The S3 wiring MUST
  preserve this: the `(2^d, ng)` moment-source carrier and the `φ̂`
  read MUST be NO-OPs when `per_axis == 1`. ⭐ **This is the
  backward-compat invariant gating the whole carve** (S2 crosswalk
  §"n_spatial_moments lattice": DD/Step at per_axis==1 stay EXACTLY).
- **If RED on a SWEEP/SOLVE snapshot** → the S3 scattering change
  leaked into the DD/Step (n=1) path (a bug — investigate, don't
  relax). The gate DISCRIMINATES: DD/Step must be inert to S3.
- Markers: the regression machinery is `@pytest.mark.foundation`.
- **Failure mode caught:** the carve leaking the moment-source change
  into the single-moment path (Mode #3 / #6).

### GATE 5 — 1-group degeneracy guard

**Claim (vv Cardinal Rule 6 / H1):** every scattering/limit claim in
gates 1+3 uses a 2G-ASYMMETRIC heterogeneous config (in ADDITION to
any 1G smoke). The slope-scattering source `Σ_s·φ̂` has a group-coupled
convention (`Σ_s^T` on the slope rows — Mode #6 convention drift in the
group axis). A 1G test cannot see a `Σ_s^T` transpose bug on the slope
rows (1G `Σ_s` is symmetric/scalar). **Enforcement, not a separate
test:** gate 1 leg 1c (2G het thick) + gate 3 (2G asym scatter) carry
this. The 1-D tripwire (leg 1a) is 1G — ACCEPTABLE only because leg 1c
adds the 2G companion. ⭐ **Reject the S3 plan if the ONLY
limit/scatter gate is 1G.** The 1-D xfail flip is necessary (it is the
#37 deliverable) but NOT sufficient — the 2G companion is mandatory.

---

# 2. THE CROSSWALK (L17 — the moment-iterate threading convention)

Pattern 7 convention crosswalk for the S3 carve (the subsystems
crossed: scalar-flux iterate ↔ spatial-moment iterate ↔ scattering ↔
moment-source carrier ↔ kernel RHS ↔ matvec probe). Written BEFORE
code. Sites verified live @ HEAD `f45a219` where line-numbered.

## The spatial-moment axis convention

Per S2 crosswalk: per-cell moments = `per_axis^d` (DD=1; LD-d1=2;
LD-d2=4; LD-d3=8), ordered tensor-Legendre `[bar, ŷ, x̂, x̂y]` in 2-D
(index `2*o_x + o_y`). The SPATIAL-MOMENT axis is a NEW trailing axis
on the iterate flux field. **It is ORTHOGONAL to the HARMONIC-moment
(angular ℓ,m) axis** — do NOT conflate. `φ̄` = (spatial-average,
angular-ℓ0) = the scalar flux today. `φ̂` = (spatial-slope rows,
angular-ℓ0). For S3 (P0/isotropic-dominant), the scattering reads the
spatial moments of the *scalar* flux (the angular-ℓ0 / 0th harmonic
moment), so the spatial-moment axis rides on the scalar flux.

| Subsystem | Input convention | Internal convention | Output convention |
|-----------|------------------|---------------------|-------------------|
| **Cell solve** (`_CellSolve.cell`, `sweep_graph.py:858`) | `Q_cells` `(N_oct or 1, ng, n_diag, 2^d)` moment source | `cell_kernel_batch` returns `psi_avg` with trailing `2^d` (d≥2 LD) | emit ACCUMULATES the FULL `2^d` spatial moments into the iterate (TODAY: reduces to slot-0 average `psi_avg[..., AVERAGE_MOMENT]` `:884`) |
| **Iterate flux field** (between-sweep) | scalar `φ̄` `(ng, *spatial)` TODAY | — | NEW: `(n_spatial_moments, ng, *spatial)` — `φ̄` is slot 0, `φ̂` the slope rows |
| **`ScatteringOperator.apply`** (`scattering.py:925`) | scalar flux / `HarmonicMomentField` (spatial-average only) TODAY | `Σ_s` matmul on the group axis | NEW: lift `Σ_s ⊗ I_{spatial-moment}` — apply `Σ_s` to EVERY spatial moment row independently |
| **Moment-source carrier** (`_SolveOperands.Q`) | scalar avg lifted to slot 0 inside kernel (S2) | — | NEW: genuine `(…, 2^d)` carrier; slot 0 = `Σ_s·φ̄ + Q_ext`, slope rows = `Σ_s·φ̂` |
| **Kernel RHS** (`_ubld_system`, `linear_discontinuous.py`) | scalar lift to `[..., 0]` (S2) | `M·S⃗_moments + F_in·traces` | NEW: read the full `(2^d, ng)` source (replacing the internal scalar lift) |
| **Matvec probe** (`_CellResidual.cell`, `sweep_graph.py:913`, RAISES `:929`) | RAISES today (needs `2^d` probe) | NEW: full `2^d`-moment probe | residual on the full moment vector; closes the raise |

## The Bridge rows (Pattern 7 — where to act)

1. **The `Σ_s ⊗ I_{spatial-moment}` lift (the LOAD-BEARING bridge).**
   `ScatteringOperator.apply` today maps `φ̄ → Σ_s·φ̄` (one spatial
   moment). Inc C maps `(φ̄, φ̂_x, φ̂_y, φ̂_xy) → (Σ_s·φ̄, Σ_s·φ̂_x,
   Σ_s·φ̂_y, Σ_s·φ̂_xy)` — `Σ_s` applied to EACH spatial moment row
   independently (the spatial-moment axis is a SPECTATOR to the
   group-matmul, exactly as the cell axis is). ⭐ **Pattern-7
   placement:** the lift belongs at the PRODUCER (`ScatteringOperator`),
   NOT at each consumer. The cleanest realization: the spatial-moment
   axis is a leading/trailing broadcast axis the existing per-material
   group matmul (`apply_p0_in_scatter` / `apply_legendre_scattering_moments`)
   is already layout-agnostic over (the docstring at `scattering.py:726`
   says `MomentProjection`/`HarmonicMomentReconstruction` are
   "layout-agnostic in their trailing axes"). So `Σ_s ⊗ I_spatial` is
   "the group matmul, vectorized over one MORE spectator axis". DD/Step
   pass a length-1 spatial-moment axis → `Σ_s ⊗ I_1` ≡ today
   (gate 4). ⚠ **VERIFY** `apply_p0_in_scatter` broadcasts over an
   extra trailing axis cleanly (probe: explorer is mapping this — flag
   the assumption).
2. **The `φ̂` accumulation in the emit** (`_CellSolve.cell:884`). Today
   the emit reduces `psi_avg` to slot 0. S3: accumulate the FULL `2^d`
   spatial moments into the `(n_spatial_moments, ng, *spatial)` iterate.
   The harmonic-moment emit branch (`:893`) and the angular branch
   (`:888`) both currently write the average; S3 widens the iterate
   write. ⚠ The scalar-flux iterate that DRIVES scattering is the
   angular-ℓ0 component; for the windowed (moment) path the spatial
   moments ride on the ℓ=0 harmonic moment. Keep the spatial-moment and
   harmonic-moment axes distinct (crosswalk §"spatial-moment axis").
3. **The matvec probe** (`_CellResidual.cell:929` raise). S3 supplies
   the `2^d`-moment probe (the iterate now carries it) → `residual_kernel_batch`
   gets its full probe → close the raise. This is gate 3's SUT.

## Mode-7/8/9 hazards specific to S3

- **Mode-9 (the BIG one — the brief's load-bearing concern).** The
  durable plan's "FP-invariance vs S2-flat-source-FP" is a
  MIS-APPLICATION of Mode-9 (§0). Inc C is NOT a splitting — it changes
  the operator. **Do NOT write a gate asserting S3-FP == S2-FP.** The
  legitimate Mode-9 gate is gate 3 (SI ≡ Krylov on the SAME `S_full`
  operator), and IT must run on the anisotropic + diagonal-cubature +
  non-zero-scatter config (else the FP-coincidence makes SI≡Krylov
  trivially — the slope source ≈ 0). The TWO Mode-9 sub-instances now
  in the D5 area (D5b.3 schedule-equiv, gate-3 SI-vs-Krylov-with-lag)
  both want the same degeneracy-break config.
- **Mode-8 (`-O` strips bare assert).** ALL S3 gates `np.testing.*` /
  `pytest.fail` / `assert_regression`. The 1-D tripwire body uses bare
  `assert rel < 0.05` (`test_mms_ld_slab.py:270`) — ⚠ **it is INERT
  under `-O`**. When flipping it xfail→PASS (leg 1a), MIGRATE the bare
  `assert` to `np.testing.assert_array_less` / `pytest.fail` (it was
  xfail so the inert-assert never mattered; once it is a LIVE PASS gate
  it MUST fire under `-O`). This is a Mode-8 catch the flip introduces.
- **Mode-7 (MMS simplification bias).** Gate 2's transport MMS reuses
  the EXISTING isotropic cases (1-D slab + 2-D het smoke) — those are
  Mode-7-acceptable for S3 (they verify the moment source didn't break
  O(h²); the smoke is honestly NOT a flux-shape claim). The
  STRENGTHENED μ-non-trivial stress ansatz (which WOULD make a
  flux-shape claim) is S4 — so S3 does NOT need the Mode-7 override.
  ⭐ But the thick-diffusion gate (gate 1) IS the real correctness
  content S3 ships — it is a VALUE anchor against the diffusion limit,
  not a rate claim, so it sidesteps the "rate ≠ correctness" trap (vv
  §5) by construction.

---

# 3. EXISTING-TEST DISPOSITION (retire / invert / migrate / stay-green)

| Test | file:line | Verdict | Reason |
|------|-----------|---------|--------|
| `test_ld_thick_diffusive_limit_xfail` | `test_mms_ld_slab.py:235` | **FLIP xfail→PASS + Mode-8 migrate** | #37/Inc C deliverable. Remove xfail, rename `..._limit`, migrate bare `assert:270`→`np.testing`. THE 1-D leg of gate 1. |
| `test_sn_1d_slab_ld_mms_converges_second_order` | `test_mms_ld_slab.py:82` | **STAY GREEN** | gate 2 — transport O(h²) survives the moment-source thread. RE-RUN; drift = bug. |
| `test_ld_2d_converges_second_order_smoke` | `test_mms_ld_2d.py:90` | **STAY GREEN** | gate 2 — 2-D smoke survives. Its scalar flux has spatial slopes → `Σ_s·φ̂` now active; O(h²) must hold. |
| `test_ld_2d_two_paths_ffw_equals_mfw` | `test_mms_ld_2d.py:142` | **STAY GREEN** | the FFW≡MFW schedule-equiv (D5b.3); S3's moment iterate runs on both schedules → must still agree. RE-RUN. |
| `test_dd_and_ld_2d_converge_to_different_values` | `test_mms_ld_2d.py:191` | **STAY GREEN** | DD≠LD discrimination (D5b.5); S3 makes LD MORE different from DD (the slope source) → the gap widens, stays green. |
| `test_cell_kernel_batch_admits_multi_d` | `test_linear_discontinuous.py` (S2 INVERTED) | **STAY GREEN** | the kernel admits d≥2; S3 doesn't re-touch the kernel admission. |
| D5b.4 kernel matvec twin (`residual_kernel_batch` at solved state) | `test_linear_discontinuous.py::TestLDKernel` | **STAY GREEN + LIFT** | S2 verifies the KERNEL matvec twin; S3 LIFTS it to the full `(L+C−S_full)` operator (gate 3). The kernel-level twin stays; gate 3 adds the operator-level end-to-end. |
| DD strict gate (`sweep/core + solve -W error::DriftWarning`) | (invocation) | **STAY at S2 baseline** (gate 4) | DD/Step bit-id negative control. RE-CONFIRM live count at pickup. |
| `_CellResidual.cell` raise | `sweep_graph.py:929` | **CLOSE** | S3 supplies the `2^d` probe; the raise goes. Gate 3 is the end-to-end exercise. |

**Green floor (MUST stay green throughout S3):** the gate-2 transport
MMS (1-D + 2-D smoke) + the D5b.3/D5b.5 foundation gates + the DD
strict gate (S2 baseline) + `test_kinf_homogeneous`(≥2G, the
structurally-indep eigenvalue anchor — UNTOUCHED by S3) + `test_mms_2d`
(DD). Route around the pre-existing reds (the d5 memo's 7-red `-k`
list); NEVER all `tests/sn` (#212).

---

# 4. SELF-IMPROVEMENT — the Mode-9 mis-application is a NEW failure shape

**FLAG AT DELIVERY (candidate vv-principles addition, NOT yet a
skill promotion):** the S3 design surfaced a NEW test-design failure
shape distinct from the catalogued Mode-9:

> **Mode-9-MISCLASSIFICATION:** treating a PHYSICS-COMPLETION (a change
> that deliberately moves the converged fixed point to a more-correct
> one) AS a splitting/acceleration (which by definition must NOT move
> the FP), and therefore writing an "FP-invariance" gate that asserts
> the NEW FP equals the OLD (incomplete) FP. The gate then PASSES on
> the degenerate regime where the completion is inactive (flux flat /
> weak scatter) and FAILS on the regime where the completion matters —
> and the failure is mis-read as a bug, rewarding REMOVAL of the
> completion. This is the inverse risk of canonical Mode-9: there the
> gate is too WEAK (degenerate-blind); here the gate is too STRONG (it
> forbids the correct FP change) AND degenerate-blind simultaneously.

The catalogued Mode-9 says "verify a splitting to the SAME FP". The
prerequisite — *first classify whether the change IS FP-preserving* —
is implicit and was dropped by the durable plan. **The
classification step is the missing discipline:** before writing any
FP-invariance gate, ask "is this change FP-preserving (splitting/
accel/schedule/representation) or FP-changing (physics completion /
operator change)?" Only the former gets a same-FP gate; the latter
gets a CORRECTNESS gate against a structurally-independent reference
(here: the diffusion limit).

This is the SECOND Inc-C-area instance of "the iteration-only vs
operator-change distinction is load-bearing" (the first: the 1-D Inc-B
flat-source caveat, [[issue_158_ld_spatial_verification]]). **RECOMMEND:
when delivering S3, add a one-line note to the vv-principles Mode-9 row
— "before gating FP-invariance, CLASSIFY: a physics-completion that
moves the FP gets a correctness gate against an independent reference,
NOT a same-FP gate."** Promote to a full Mode-10 only if a real
production bug is later shown to have shipped behind a same-FP gate on
a physics-completion. For now: documented here.

---

# 5. DELIVERY CHECKLIST (the S3 implementer's gate sequence)

1. **Re-confirm the FP semantics** (§0): S3 converges to
   `(L+C−S_full)`, S2 to `(L+C−S_flat)`. Do NOT write a S3-FP==S2-FP
   gate.
2. **Gate 4 FIRST** (the negative control): wire `Σ_s ⊗ I_spatial`
   such that `per_axis==1` (DD/Step) is a NO-OP. Re-run the DD strict
   gate → MUST stay at the S2 baseline. (Land the backward-compat
   invariant before the new physics.)
3. **Gate 2** (transport O(h²) survives): RE-RUN the 1-D LD slab MMS +
   the 2-D smoke + D5b.3/D5b.5. Drift = bug (the moment-source wiring
   disagrees with the MMS PDE — investigate, don't relax).
4. **Gate 3** (the genuine Mode-9): close the `:929` matvec raise; the
   end-to-end SI ≡ Krylov on `(L+C−S_full)` on the anisotropic +
   diagonal-cubature + non-zero-scatter NON-SQUARE 2G config.
5. **Gate 1** (THE correctness gate): flip the 1-D tripwire xfail→PASS
   (Mode-8-migrate the bare assert) + the 2-D analog leg + the 2G
   companion (gate 5). Anchor against the diffusion limit (DD-proxy +
   continuous-diffusion value band).
6. **Gate 5** (1G guard): confirm gates 1+3 carry a 2G het config.
   REJECT if the only scatter/limit gate is 1G.
7. **Claims ledger** in the PR: S3 UNLOCKS {Krylov d≥2 matvec, the
   limit, the `(2^d,ng)` carrier}; STAYS S4 {the strengthened
   flux-shape ansatz + non-vanishing boundary trace + the
   `verifies("ld-cartesian-2d")` claim}.
8. **Self-improvement** (§4): note the Mode-9-misclassification to the
   vv-principles Mode-9 row at delivery.
9. Route around the 7 pre-existing reds; NEVER all `tests/sn` (#212).

---

## Cross-links
- Extends [[d5-nd-polymorphism-verification]] §D5b (S2/D5b.0–.5; this
  memo owns S3, the moment-iterate fold). The S2 crosswalk is
  `.claude/plans/issue_240_d5b_s2_crosswalk.md` (moment ordering, the
  Q-lift S3 replaces, the n_spatial_moments lattice / backward-compat
  invariant).
- [[issue_158_ld_spatial_verification]] (the 1-D Inc-B flat-source
  caveat = the 1-D precursor of this FP-changing-vs-iteration-only
  distinction; the thick-diffusive tripwire originated there).
- [[feedback_regression_tolerance_design]] (iterative → SAFETY×conv_tol
  for gates 2/3; DriftWarning strict for gate 4).
- [[feedback_vv_tagging]] (the xfail strict=False → flip-to-PASS
  discipline for gate 1 leg 1a).
- [[issue_206_phase_c_verification]] (FINDING-1: the
  `tests.sn.regression...` DriftWarning escalation path for gate 4).
- vv-principles: Mode 9 (the mis-application is §0/§4); hierarchical
  taxonomy (gate 1 = flux-shape/value claim, NOT eigenvalue — MMS
  can't prove eigenvalues, and the diffusion limit is the structurally-
  independent value reference); Cardinal Rule 6 / H1 (gate 5);
  bit-identity-vs-principled (gate 4 = the negative-control bit-id).
- Literature: Adams-2001 (the multi-D thick-diffusion verdict, the
  basis-discrimination), BLA-1992 (2-D LD asymptotic), LMM-1987 (the
  four-limits criterion + DD's 1-D diffusion-limit consistency = the
  gate-1 DD-proxy ground). `.claude/agent-memory/literature-researcher/multi_d_ld_closure.md`
  §5 (the diffusion limit — "flat-source UBLD is O(h²) but NOT
  diffusion-limit-consistent until the moment source is threaded" =
  the verbatim FP-semantics confirmation).
