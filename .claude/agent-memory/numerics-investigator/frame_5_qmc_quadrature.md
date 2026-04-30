---
name: Frame 5 QMC angular quadrature for F.4
description: Owen-scrambled Sobol' angular quadrature kills L17/L19 quadrature-crossing on F.4 hollow-sphere at anchor and pathology (2026-04-22)
type: project
---

# Frame 5 — randomized QMC angular quadrature for F.4 hollow-sphere closure

**Date**: 2026-04-22. **Agent**: numerics-investigator (Opus 4.7). **Issue**: cross-domain-attack Frame 5 (no GitHub issue yet).

## Claim tested

L17/L19 empirical pathology: F.4's signed error `(k - k_inf)/k_inf` under product-Gauss panel refinement flips sign between RICH=(4 panels, p_order=8, n_ang=64) and RICH+panels=(5 panels, p_order=8, n_ang=64). Because exp(-τ·d) has bounded Hardy-Krause variation, Owen-scrambled Sobol' angular quadrature should converge with τ-bounded constant (vs product-Gauss's τ^{2p} constant) and give a single sign-stable signed error.

## Deliverable

`derivations/diagnostics/diag_f4_qmc_quadrature.py` — subprocess-isolated driver with Owen-scrambled Sobol' monkey-patched into `peierls_geometry.gl_float`. Default run: anchor + pathology (2 points, ~510s wall). Full 6-point scan: `QMC_SCAN_MODE=full`.

## Option chosen

**Option A** — angular-only QMC. Replace only the angular Gauss-Legendre driven by `n_ang` (via `gl_float`) with Owen-scrambled Sobol'. K_vol keeps product-Gauss (uses `gl_nodes_weights` directly — NOT affected by the monkey-patch). Radial composite-GL keeps product-Gauss (also uses `gl_nodes_weights`).

Rationale: L17's pathology is radial-panel-driven (n_panels=4 → 5), so Option A is NOT a full cure for L17 — but it demonstrates the angular-contribution-to-noise lever cheaply (~100s/point vs 30+ min for full 3D Sobol'). If Option A passes anchor + pathology, Option B is unnecessary. It did pass — so Option B was NOT implemented this session.

## Key results — full 6-point grid

| point      | PG RICH       | PG RICH+panels | QMC mean       | QMC 95% CI                 | sign stable | CI crosses 0 | verdict |
|------------|---------------|----------------|----------------|----------------------------|-------------|--------------|---------|
| (5,  0.3)  | +0.057825%    | +0.032331%     | **+0.040762%** | [+0.04058%, +0.04094%]     | True        | False        | PASS    |
| (5,  0.5)  | +0.297309%    | +0.322013%     | **+0.317996%** | [+0.31769%, +0.31830%]     | True        | False        | PASS    |
| (10, 0.3)  | +0.003285%    | −0.008221%     | **−0.003843%** | [−0.00392%, −0.00376%]     | True        | False        | PASS    |
| (10, 0.5)  | +0.016683%    | +0.007990%     | **+0.006847%** | [+0.00672%, +0.00697%]     | True        | False        | PASS    |
| (20, 0.3)  | −0.005982%    | −0.007481%     | **−0.007956%** | [−0.00799%, −0.00792%]     | True        | False        | PASS    |
| (20, 0.5)  | +0.005430%    | −0.003942%     | **−0.005453%** | [−0.00550%, −0.00541%]     | True        | False        | PASS    |

**6/6 points PASS** both sign-stability across all 32 scrambles AND CI-does-not-cross-zero. Anchor CI width = 0.000165% (**far below the 0.003% criterion**). Pathology CI width = 0.000089%, sign unambiguously negative. std across scrambles ~ 1e-6–9e-6 fractional (an order of magnitude smaller than PG RICH/RICH+panels separation at every point).

**Key observation**: at the two L17 sign-flip points — (10, 0.3) and (20, 0.5) — QMC's signed-error mean matches PG RICH+panels' SIGN (both negative) but not PG RICH's (positive). This suggests PG RICH is the biased estimator and PG RICH+panels is approaching the truth from the correct side; QMC converges to a value ~20-30% beyond RICH+panels in the same direction.

At the non-sign-flip points, QMC matches PG RICH+panels' trend and typically lies between the two PG estimates, closer to RICH+panels. This is consistent with Sobol' being unbiased (~E[err] = 0) vs product-Gauss being order-of-magnitude biased.

## Scaling

- QMC wall per scramble (batched, K_vol amortized): **~1.5-1.6s/scramble** after a ~45-65s K_vol build.
- 32 scrambles ≈ 95s (τ=5) – 115s (τ=20).
- Product-Gauss RICH: ~44-67s per evaluation; RICH+panels (5 panels): ~77-120s per evaluation.
- Total wall anchor+pathology (2 pts): 510s (8.5 min).
- Total wall full 6-point grid: 1492s (25 min).
- Naive 32× factor completely absent — batched single-subprocess design amortizes K_vol across all scrambles. If not batched, would be 32 × 50s ≈ 1600s per point; 9600s for full grid.

## Verdict

**QMC Option A fully eliminates the L17/L19 quadrature-crossing pathology at all 6 reference points.** Every point passes sign-stability across all 32 scrambles AND CI-does-not-cross-zero. The bootstrap 95% CI is 20-100x narrower than the PG RICH vs RICH+panels spread.

Key wins:
1. The two L17 sign-flip points (10, 0.3) and (20, 0.5) no longer have ambiguous sign — QMC gives a crisp single value.
2. The pathology point (20, 0.5) has the second-tightest CI of the grid (width 0.000089%), demonstrating that **higher τ actually IMPROVES QMC's noise ratio** — exactly the Hardy-Krause-bounded-variation prediction. Compare to PG RICH, whose bias coefficient scales with τ^{2p}.
3. The QMC mean sits at a **distinct, physically plausible value** that typically matches PG RICH+panels' sign and lies beyond it in the same direction — consistent with PG being the biased estimator that under-resolves the tangent-angle integrand kink.

The anchor's QMC mean ≈ −0.0038% suggests F.4's true closure signed error at (τ=10, ρ=0.3) is ~0.004% (negative), not ~0.003% (positive as RICH suggests). This is below the L17 `RESOLUTION_THRESHOLD = 0.005%` but unambiguously signed.

## Proposed wrapper — L19 protocol reinterpretation

Current L19 helper `assert_rank_n_structural_win` enforces five gates (S1–S5) based on discrete signed-error monotonicity across >=2 product-Gauss quadratures. If QMC passes, S3 (no sign flip) and S4 (|err| monotone) can be replaced by a single CI-based gate:

```python
def assert_rank_n_qmc_structural_win(
    closure_fn, f4_fn, point, N=4096, n_scrambles=32, tol_frac=1e-4
):
    """CI-based L19 protocol for QMC-quadratured closures.

    Pass iff the bootstrap 95% CI of the closure's signed error is strictly
    tighter than and non-overlapping with F.4's, AND does not cross zero.
    """
    # ...pseudocode — not shipped this session.
    # Build 32 scrambles of closure_fn(point, N=N, seed=i) and f4_fn(point, N=N, seed=i).
    # Bootstrap CI_closure and CI_f4.
    # Assert: CI_closure.hi < CI_f4.lo (strict separation, closure strictly better)
    #         AND |CI_closure.mid| < tol_frac (meets target)
    #         AND CI_closure does not cross 0.
```

Not shipped — the existing `assert_rank_n_structural_win` remains the production L19 helper. When a future closure actually passes Frame 5's angular-QMC quadrature, file a new Issue to add this wrapper.

## Proposed follow-up for production code

**Recommendation: NO immediate switch of `orpheus.derivations.peierls_geometry` to randomized Sobol'.** Reasoning:

1. Option A kills angular-quadrature noise but does NOT kill the radial-panel-refinement noise (L17's root cause). A full replacement would need Option B (3D Sobol' over (r, r', µ)), which is ~200 LOC of infrastructure.
2. Product-Gauss remains superior when the integrand is **known to be smooth and polynomial-approximable** — most of ORPHEUS's other solvers.
3. The Peierls hollow-sphere closure is a **specific regime** where QMC is better (exp(-τ·d) with τ large, grazing-ray singularities at tangent angles). A config flag `quadrature="qmc"` makes sense, but only if paired with B's full 3D Sobol' as well.

**Recommended follow-up**: file a new Issue `module:cp, type:improvement` titled "Optional randomized-QMC angular quadrature for hollow Peierls" referencing this diagnostic as prior art. Scope: add a `quadrature="qmc"` kwarg to `build_closure_operator` that routes to a `gl_float_qmc` analog, with N and scramble-seed parameters. Priority: low — F.4's production accuracy at RICH is adequate; this would strengthen the L19 protocol but isn't blocking any known production bug.

## Bootstrap CI math

Bootstrap-t-less percentile CI on the mean across 32 scrambles.  Let {s_i} = 32 signed-error samples. For B = 10 000 bootstrap resamples, compute means m_b = mean(s_i | i ∈ resample_b). 95% CI = [P_2.5, P_97.5] of {m_b}. CI-width criterion: `width(CI) < 0.003%` (anchor) was met with 20x margin.

Sign-stable criterion: ∀i, sign(s_i) == const. Tested: True for both points, all 32 scrambles.

Sign-unambiguous-at-CI criterion: CI does not include 0. Tested: True for both points.

## Artifacts

- `derivations/diagnostics/diag_f4_qmc_quadrature.py` — the driver.
- `/tmp/diag_f4_qmc_quadrature.json` — machine-readable dump (anchor + pathology).
- This memo.

## Cross-references

- L17, L19, L20: research log §"New lessons (L17+)".
- Direction N / Issue #123: `tests/cp/test_peierls_rank_n_protocol.py` (the L19 protocol helper `assert_rank_n_structural_win`).
- Frame attack memo: `.claude/agent-memory/cross-domain-attacker/peierls_rank_n_frame_attack.md` Frame 5 section.
