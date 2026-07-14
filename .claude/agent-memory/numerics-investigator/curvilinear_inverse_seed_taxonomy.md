---
name: curvilinear-inverse-seed-taxonomy
description: Ruling + recipe — is the 1-D curvilinear (L+C).solve a seed-independent direct inverse (SweepOperator) or a seed-lagged iteration? Slab/cyl YES, sphere NO. Seed-entry map, coupling structure, feasibility, MMS-blindness.
metadata:
  type: reference
---

# Curvilinear `(L+C).solve`: SweepOperator or lagged iteration? (#226 taxonomy ruling)

Diagnosed on `refactor/inverse-as-operator` (2026-07). READ-ONLY probe of whether the
curvilinear WDD sweep is an honest direct triangular inverse. Backs `lessons.md` L14.

## Verdict (per geometry)
- **SLAB** = honest `SweepOperator` (direct inverse): seedΔ **0.0** (bit), residual **8e-16**.
- **CYLINDER = QUADRATURE-DEPENDENT** (the old "CYLINDER YES" here was over-generalized
  from a LEVEL-SYMMETRIC probe; CORRECTED 2026-07-05, #280 Phase 2.5b, product-cyl repro
  cold err **0.575**):
  - **level-symmetric** = trivially direct (cold err **7.8e-16**) because the first-swept
    ordinate has `c_in[m0]=0` (raw τ=1 → (1−τ)/τ=0) — a **DEAD** seed, not "α-dome
    telescoping" (that mechanism attribution was WRONG).
  - **product** = **SEED-LAGGED** (cold err **0.575**, threaded→M⁻¹ at 6.7e-16): the #229
    τ-clamp maps raw τ=0 → ½, so `(1−τ)/τ=1` → `c_in[m0]≠0` → the seed = ψ_{m0} (t=0,
    first-swept ordinate's own flux) is a LIVE self-coupling read off the ITERATE.
    **RETIRABLE by a direct single-pass fold** — see the #280 §below.
- **SPHERE** = WAS seed-DEPENDENT lagged (seedΔ **4.6e-2**, cold residual **5e5**), **FIXED
  by route (a)** — the direct `carlson_inward_sweep_from_source` starting-direction march
  (LANDED `a29ab2d`, L15). No longer a lag.

## #280 Phase 2.5b — product-cyl seed lag is FEASIBLE to retire by a single-pass fold
Ruling (2026-07-05, read-only feasibility): **FEASIBLE (clean pure-diagonal fold).**
- **Root cause** — the ONLY lagged element is the seed. Seeding the cold solve with the
  CONVERGED flux → single-pass exact (6.7e-16). Seed = ψ_{m0}, t=0 EXACTLY, m0 =
  first-swept ordinate (SAME ordinate it feeds) → a per-ordinate self-reference.
- **The self-coupling is a PURE DIAGONAL fold.** The augmented `(L+C)` matrix M is EXACTLY
  walk-order triangular (triu=0); the seed ordinate m0's output depends on NOTHING but
  itself (other-ordinate coupling = 0.0), and its nx×nx self-block is triangular. Zeroing
  the seed in the matvec changes ONLY the m0 self-block **diagonal** (off-diag change = 0.0
  exact) — so the seed self-coupling folds into the m0 cell diagonal with coefficient
  **`κ = dA_w[m0]·c_in[m0]`**, exactly as self-scatter folds. The "seed telescopes away"
  story (this file's old text, the `test_282` docstring line ~629, ERR-026 docs) is a
  MIS-ATTRIBUTION for product — the seed genuinely contributes O(1) to the diagonal; M is
  triangular because that contribution is ON the block-diagonal (forward-sub handles it).
- **POC** — inject the local triangular solve of each m0 self-block as the seed → cold
  solve = M⁻¹ single-pass (5e-16). The direct fold reproduces M_aug⁻¹ bit-close.
- **Change site** — `_run` non-carrying `else` branch, `loss_representation/__init__.py`
  **4091–4101** (`phi_aux = closure.edge_extrapolated_seed(psi_level, p_idx)` reads the
  iterate). The **MATVEC needs NO change** (`precompute_psi_state`/`cell_contribution`
  build the seed from their input = already the correct triangular operator; the solve just
  must reproduce its forward-substitution). Fix = solve m0 in-place with κ folded into its
  cell diagonal (NOT read from the iterate) — structurally the cylinder analogue of the
  sphere's route (a) starting-direction solve, but with a live self-redistribution fold.
- **Re-baseline** — the fixed point is **MACHINE-IDENTICAL** (lagged-SI-converged ≡ folded ≡
  M⁻¹ to ~5e-16, different FP paths, far under any keff/MMS rtol). UNLIKE sphere route (a)
  (which changed the seed VALUE/closure → moved keff(N) at
  fixed N, L15), the cyl is NON-carrying: the converged seed IS ψ_{m0} either way, so the
  fold changes only HOW ψ_{m0} is obtained (iterate-lag → in-place). keff / MMS / frozen
  `walk_matvec_cyl_2g` matvec baselines / converged snapshots **DO NOT MOVE**; only a
  cold/intermediate single-within-group-sweep value moves (principled: now exact).
- **LS + degenerate** — LS is INERT (c_in[m0]=0, dead). Degenerate pure-azimuthal ordinates
  (μ_r≈0, is_degenerate=[2,6,10,…]) sit MID-level DOWNSTREAM of m0 (m0 = most-negative μ_x,
  never degenerate) → they inherit the correct M-M thread once m0 is folded; no special
  handling. τ-clamp #229 is NOT an obstruction (it DEFINES κ via (1−τ)/τ=1; triu=0 stands).
- Diagnostics: `derivations/diagnostics/diag_280_cyl_product_seed_lag.py` (characterization
  + strict-xfail TARGET gate) + `diag_280_cyl_fold_feasibility.py` (structural proofs +
  POC). Promote both to `tests/sn/sweep/` (the feasibility proofs are permanent
  routing/triangularity gates; the xfail flips green when the fold lands).

## Where the seed enters (file:line)
SWEEP (`orpheus/sn/loss_representation/__init__.py`):
- `_initial_guess_values` (2098–2139) extracts `initial_guess.bulk.values`; `None`→`None`.
- `_run` curvilinear arm (3149–3323): `ig_values=_initial_guess_values(...)` (3162);
  `psi_level = psi_g_first[:,ords,:]` or **zeros** on cold start (3185–3188);
  `phi_aux = closure.psi_half_seed(psi_level, carlson_ctx)` (3197) — **THE lagged seed**;
  `psi_angle = phi_aux.copy()` (3198). The M-M recurrence (3300–3304) and the pole r=0
  continuity capture/consume (`pole_outflow`, 3228/3286) are **feed-forward within the
  current sweep** — NOT lagged.
MATVEC (`orpheus/sn/spatial/pole_angular_closure.py`): `precompute_psi_state` (838–886)
builds the WHOLE half-grid from the apply INPUT ψ via the same `psi_half_seed` (881).
`InvertibleOperator.apply=loss_action(self.sigma,psi)` (streaming.py:758);
`.solve→loss_representation.sweep(...,initial_guess=)` (streaming.py:858→__init__.py:1037).
Default seed = **`AngularEdgeExtrapolation`** (streaming/pole closure `__init__` 594–600),
NOT `CarlsonInwardSweep` (ERR-058 §195 superseded it); reads the two most-inward ordinate
cell-AVERAGES (psi_half_angle_seed.py 789–811), exact on flat AND **linear-in-μ**.

## Coupling structure = one LOCAL CYCLE (not global, not feed-forward)
The extrapolation seed reads ordinate-0's cell-average; ordinate-0's redistribution reads
the seed → **seed ↔ first-ordinate cycle** (per level). The sweep breaks it by LAGGING
(reads previous iterate); a matvec builds it from its input (consistent). ⇒ `solve` inverts
`apply` ONLY at the fixed point (probe6: `solve(apply(ψ0),ig=ψ0)`→5e-13; cold→24–3900×).
The seed-iteration `ψ←solve(b,ig=ψ)` on the BARE operator **DIVERGES** ×~70/iter.

## Feasibility — a direct sphere sweep IS implementable (seed-dep is a CHOICE)
The μ=±1 starting-direction equation is CLOSED: `(1−μ²)=0` kills redistribution → pure
streaming+collision. `carlson_inward_sweep_from_source` (psi_half_angle_seed.py:428) ALREADY
implements it (Hébert 3.434–3.435). Direct options: (a) carry ψ(·,μ=−1) as EXPLICIT
per-level state (cleanest — the hidden state made explicit); (b) small per-level block-solve
over {seed, first 1–2 ordinates × cells}; (c) source-driven seed from the RHS moment (breaks
matvec/sweep symmetry — a pure operator has no "source"). This is exactly **issue #200**'s
"block-inverse face preconditioner".

## Why the MMS suite is blind (Mode-7) + how production dodges it
Every curvilinear MMS ansatz is ≤ linear-in-μ (`sin(πr/R)`, `(A(r)+B(r)μ)/W`) = the seed's
EXACT regime → SI converges O(h²) on the whole ladder, seed-lag INVISIBLE (`solve_sn_fixed_
source` docstring claim "SI≡Krylov on MMS ladders" holds ONLY there). A higher-order-in-μ
field (uniform-source sphere) → SI **NaN** at every c∈[0,0.99]. **But config-dependent**
(2026-07-01 falsifier re-run): the uniform-source sphere `solve_sn_fixed_source` gives SI
**NaN** + Krylov **negative flux** (min −163) on a COARSE `level_symmetric`(S4)/16-cell vacuum
body, yet PHYSICAL flux (∈[0.84,1.99], both SI & Krylov) on a `gauss_legendre`(S16)/40-cell
body — slab/cyl clean in BOTH; the existing `test_unclamped_sphere_flux_strictly_positive`
(GL, `solve_sn` eigenvalue) passes. So "sphere production is broken" is an over-claim: the
NaN/negativity is the DOWNSTREAM consequence of the seed-lag on a stress quadrature/mesh, not
a blanket failure — quad+resolution is a real axis. Production dodges: eigenvalue
keff is shape-independent (SI ok, "flux shape closure-bias drift" admitted, solver.py:1707);
curvilinear GMRES ships the **IDENTITY** precond (`_within_group_krylov` solver.py:332,
"cold-seed sweep destabilises curvilinear precond" per `test_krylov_curvilinear_precond_
safety.py`), with #200 tracking the real one.

## Diagnostics (tmp — promote candidates)
`/Users/rodrigo/.claude/jobs/84fd66f8/tmp/diag_curvilinear_seed_sensitivity.py` (slab/cyl
seed-independence + machine residual = regression gate; sphere seed-dep = xfail-until-#200)
· `diag_sphere_fixedpoint_consistency.py` (fixed point is CORRECT). Natural home
`tests/sn/operators/` beside `test_removal_form_matvec_sweep.py` (whose `_ROUNDTRIP_CASES`
already EXCLUDES sphere for this exact reason, lines 337–346).

## Falsifier recipes (#226 plan §16 — re-run 2026-07-01, all CONFIRMED)
- **F3 round-trip discriminator** (`falsifier_f3.py`): `inv.apply(A.apply(x))==x` cleanly
  separates the honest `(L+C).inverse()` SweepOperator (slab 2G vacuum: **3.84e-16**) from the
  confused `_GaussSeidelResolvent` (2-D reflective: **2.667** O(1) RED — its `apply`=(L+C) and
  `solve`=(L+C−B_lower)⁻¹ invert DIFFERENT operators, the plan's W2). Control = clean (L+C)
  round-trip on the same 2-D mesh **7.02e-16**. Reach the G-S resolvent via
  `_select_si_resolvent(LC, None, B, sn2, "gauss_seidel")` (2-D Cartesian ONLY; 1-D falls back
  to Jacobi). NO existing gate runs the round-trip on the G-S resolvent (net-new). Clean side
  IS gated: `test_removal_form_matvec_sweep_roundtrip` (via .solve) + `test_inverse_operator_
  equivalence` (inverse≡solve).
- **F5 G-Neumann constructibility** (`falsifier_f5.py`): `Σ_k ((L+C)⁻¹S)^k (L+C)⁻¹ q` built
  with `inv=LC.inverse()`∘`S.apply` decays geometrically ratio **0.4995≈c=0.5** and converges
  to a structurally-independent **dense LU** solve of (L+C−S) at **4.19e-13** (τ=40 mfp slab,
  1G). The SI iterate from a zero start IS this partial sum; the dense LU (direct elimination)
  is the independent oracle since the sweep-inverse is exercised ONLY on the Neumann side.
