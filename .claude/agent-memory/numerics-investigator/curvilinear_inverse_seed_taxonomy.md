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
- **SLAB / CYLINDER** = honest `SweepOperator` (direct inverse): seedΔ **0.0** (bit),
  `‖Aψ−b‖∞/‖b‖∞` **~8e-16 / 7e-16**. Cylinder is seed-blind because each μ-level's
  α-dome telescopes to 0 at the level edges (absorbs the seed) — same mechanism the
  ERR-026 docs cite for "cylinder MMS passes despite wrong seed".
- **SPHERE** = seed-DEPENDENT lagged/preconditioned iteration, exact ONLY at the SI
  fixed point. seedΔ(X1,X2)=**4.6e-2**, cold-seed round-trip residual **5e5** (abs 1.8e6),
  localized to the outer |μ|=0.816 level. Calling it a machine-precision one-shot inverse
  is DIShonest.

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
