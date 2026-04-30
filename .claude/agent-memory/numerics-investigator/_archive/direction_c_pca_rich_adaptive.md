---
name: Direction C — PCA sectors with scale calibration at RICH quadrature
description: Sanchez-Santandrea 2002 PCA sectors + per-sector scale DOF tested at RICH quad on the hollow-sphere rank-N closure. Full 6-point grid complete. Hypothesis RH-C falsified — RICH-wins at 3/6 points are all L17 quadrature-crossing artifacts under RICH+panels verification.
type: project
---

# Direction C — PCA sectors at RICH quadrature (Issue #121, 2026-04-22)

## STATUS: FALSIFIED (complete)

6-point grid + M-curve + 3-point crossing-verification complete. PCA M=3
uniform-Brent @ RICH produces 3 apparent wins (5.0/0.3, 5.0/0.5, 20.0/0.5)
with ratios 4.5×–10× vs F.4. All 3 were probed with a 3-quadrature
stability check (RICH, RICH+panels=(5,8,64), RICH+pp=(5,10,64)) and all 3
FAIL the stability test:

- Two of three show PCA *signed* error flipping sign between RICH and
  RICH+panels — textbook L17 quadrature-crossing signature.
- The third (5.0/0.3) shows non-monotone PCA error that does not
  converge under refinement.

F.4 does the same thing at (20.0/0.5) — its own structural floor is
below its quadrature noise. That was already known from the pre-kill
anchor probe (L17).

**Conclusion**: PCA M=3 with uniform-scale calibration has no
structural-floor advantage over F.4. The RICH-win ratios are
quadrature-residual cancellation, not truncation-residual reduction.
**RH-C refuted.**

## Protocol (as run)

- Optimization at RICH (4, 8, 64). Brent 1D on uniform α = (s, ..., s).
- Reporting at RICH (same quadrature).
- Stability check at 3 richer quadratures (RICH vs RICH+panels vs
  RICH+pp) on the 3 surprise wins. ULTRA (5, 10, 96) attempted but
  single evals exceed 600s wall; dropped in favor of the more feasible
  RICH+panels and RICH+pp.
- Subprocess per-point with `timeout 600`–1200 s shell wall. Originally
  one point ran >700 s (Brent fallback to `bounded` on σ_t·R=20 due to
  mid-bracket not <right); rerun at 1200 s succeeded.

## Results

### Six-point reference grid — PCA M=3 uniform-scale Brent @ RICH

| σ_t·R | ρ   | F.4 RICH   | PCA α=1 RICH | α*     | PCA α* RICH | PCA/F.4 |
|-------|-----|------------|--------------|--------|-------------|---------|
| 5.0   | 0.3 | 0.05782%   | 1.14708%     | 1.1770 | 0.01279%    | 0.22 *  |
| 5.0   | 0.5 | 0.29731%   | 1.06849%     | 1.0671 | 0.02867%    | 0.10 *  |
| 10.0  | 0.3 | 0.00328%   | 0.93325%     | 1.2294 | 0.00833%    | 2.54    |
| 10.0  | 0.5 | 0.01668%   | 1.01288%     | 1.1197 | 0.02633%    | 1.58    |
| 20.0  | 0.3 | 0.00598%   | 0.43402%     | 1.1502 | 0.01123%    | 1.88    |
| 20.0  | 0.5 | 0.00543%   | 0.49227%     | 1.0909 | 0.00055%    | 0.10 *  |

Starred (*) = PCA apparently beats F.4 at RICH. All 3 shown below to be
quadrature-crossing artifacts.

### M-curve @ anchor (σ_t·R=10, ρ=0.3) RICH

F.4 @ anchor RICH: **0.00328%**.

| M | α=1 RICH | α*     | PCA α* RICH |
|---|----------|--------|-------------|
| 2 | 0.93325% | 1.2385 | 0.00718%    |
| 3 | 0.89123% | 1.2294 | 0.00833%    |
| 5 | 0.89123% | 1.2179 | 0.00803%    |

PCA plateaus at ~0.008% for M ∈ {2, 3, 5}: **M-independence confirms
the PCA truncation floor at the anchor ≈ 2.5× F.4 RICH**. No benefit
to more sectors.

### Stability check — signed error at 3 quadratures

At each PCA α*, evaluated signed k-residual at three progressively
richer quadratures. **A structural win keeps sign and magnitude.**

| σ_t·R | ρ   | α*     | quad         | F.4 signed   | PCA signed   | comment |
|-------|-----|--------|--------------|--------------|--------------|---------|
| 5.0   | 0.3 | 1.1770 | RICH         | +0.05782%    | +0.01296%    | PCA wins |
| 5.0   | 0.3 | 1.1770 | RICH+panels  | +0.03233%    | +0.05086%    | PCA LOSES (magnitude grew 4×) |
| 5.0   | 0.3 | 1.1770 | RICH+pp      | +0.05830%    | +0.02325%    | PCA wins again (unconverged) |
| 5.0   | 0.5 | 1.0671 | RICH         | +0.29731%    | +0.02846%    | PCA wins (small) |
| 5.0   | 0.5 | 1.0671 | RICH+panels  | +0.32201%    | **-0.07598%**| SIGN FLIP — L17 crossing |
| 5.0   | 0.5 | 1.0671 | RICH+pp      | +0.38548%    | -0.05668%    | still negative |
| 20.0  | 0.5 | 1.0909 | RICH         | +0.00543%    | +0.00026%    | both near-zero |
| 20.0  | 0.5 | 1.0909 | RICH+panels  | **-0.00394%**| **-0.05822%**| SIGN FLIP (both!) |
| 20.0  | 0.5 | 1.0909 | RICH+pp      | +0.00974%    | -0.04643%    | F.4 flips again; PCA still neg |

(5.0, 0.3): PCA signed err non-monotone 0.013 → 0.051 → 0.023 — doesn't
converge. Both F.4 and PCA noise-floor-limited at ≲ 0.06%.

(5.0, 0.5): F.4 monotone increasing (truncation-residual-dominated), PCA
sign-flips — crossing artifact. Apparent win is real in |·| but PCA is
unconverged.

(20.0, 0.5): Both F.4 and PCA in the noise floor of the discretization.
F.4 signed err swings 0.005% → -0.004% → +0.010%. PCA signed err swings
0.0003% → -0.058% → -0.046%. The RICH 0.0003% was a coincidental zero
crossing.

### Quadrature-crossing at anchor (pre-computed perturbation scan)

Retained from pre-kill session. At uniform α*=1.2294, M=3, (σ_t·R=10,
ρ=0.3):

| α                          | RICH err  | notes |
|----------------------------|-----------|-------|
| [1.229, 1.229, 1.229]      | 0.0083%   | uniform Brent optimum |
| [1.329, 1.229, 1.229]      | 0.0002%   | α_0 perturbed +0.1    |
| [1.320, 1.229, 1.229]      | ~0%       | sign flip near 1.30–1.35 |

| α                          | RICH err  | ULTRA err | Δ    |
|----------------------------|-----------|-----------|------|
| [1.229, 1.229, 1.229]      | 0.00833%  | 0.01347%  | RICH < ULTRA (crossing nearby) |
| [1.320, 1.229, 1.229]      |           | 0.00662%  | RICH "minimum" is ULTRA local |

F.4 reference: RICH=0.00328%, ULTRA=0.00102%.

## Per-sector perturbation — full table from anchor

(Pre-computed, uniform α*=1.2294 baseline.)

| α vector                   | RICH err   |
|----------------------------|------------|
| [1.229, 1.229, 1.229]      | 0.00833%   |
| [1.129, 1.229, 1.229]      | 0.01620%   |
| [1.179, 1.229, 1.229]      | 0.01229%   |
| [1.279, 1.229, 1.229]      | 0.00409%   |
| [1.329, 1.229, 1.229]      | 0.00019%   |
| [1.229, 1.129, 1.229]      | 0.12565%   |
| [1.229, 1.179, 1.229]      | 0.06950%   |
| [1.229, 1.279, 1.229]      | 0.05843%   |
| [1.229, 1.329, 1.229]      | 0.13121%   |
| [1.229, 1.229, 1.129]      | 0.42129%   |
| [1.229, 1.229, 1.179]      | 0.23353%   |
| [1.229, 1.229, 1.279]      | 0.26322%   |
| [1.229, 1.229, 1.329]      | 0.59250%   |

α_2 (grazing sector) 10× more noise-sensitive than α_0. Consistent with
the RICH quadrature being strained in the grazing direction.

## Uniform-scale extended scan @ anchor RICH

| s (uniform α) | RICH err |
|----|------|
| 0.90 | 1.121% |
| 0.95 | 1.028% |
| 1.00 | 0.918% |
| 1.05 | 0.785% |
| 1.10 | 0.624% |
| 1.20 | 0.181% |
| 1.2294 (Brent optimum) | 0.008% |
| 1.30 | 0.521% |
| 1.40 | 1.668% |
| 1.60 | 6.267% |

Steep minimum at s ≈ 1.23. Similar shape to L14's 1.8066 finding for
rank-(1,1,1) Legendre (numerical values differ: different basis → 
different "gauge" reference).

## Final verdict

**Does Direction C beat F.4 at any point in the 6-point grid at RICH?**

At face value, yes — PCA-RICH beats F.4-RICH at 3 of 6 points
(ratios 0.10, 0.22, 0.10). **But** all three apparent wins fail the
stability test: PCA signed-error flips sign under quadrature refinement
(2 of 3) or fails to converge (1 of 3). The wins are quadrature-noise
cancellation artifacts, identical in mechanism to L17. No structural
advantage.

F.4 itself shows the same pathology at (20, 0.5) — its signed error
flips sign between RICH and RICH+panels. F.4's true structural floor
is below RICH's quadrature noise at high σ_t·R. The lesson L17 applies
symmetrically: to claim any rank-N closure beats another at high σ_t·R,
the comparison must be at quadrature that RESOLVES both closures, not
merely RICH.

**Verdict: PCA M=3 uniform-Brent has no structural-floor advantage over
F.4. At best it matches F.4 within their shared quadrature noise.
RH-C refuted. Issue #121 should be closed with "falsified, converges to
the L17+L14 story".**

## Lessons to promote

- **L17 (confirmed, complete)**: F.4 RICH = 0.003% is quadrature-noise-
  limited. ANY comparison at RICH where the alternative is also tuned at
  RICH is NOT a structural-floor comparison. The L16 retraction
  principle extends: every claimed "RICH-win" below F.4-RICH must be
  verified at ≥ 1 richer quadrature. If the signed error flips sign or
  fails to decrease monotonically, it is a crossing artifact.

- **L18 (confirmed)**: Per-sector α tuning inherits the quadrature-noise
  coupling of each sector. The grazing sector c ∈ [2/3, 1] has ~10×
  the RICH-noise-sensitivity of the polar sector. This is NOT a
  structural-basis property but a quadrature-mesh interaction.

- **L19 (new)**: Single-point RICH-numbers are now shown to be
  misleading on BOTH sides of the rank-N comparison (F.4 AND PCA).
  Closures must be compared at quadrature resolving the smaller of
  their structural floors — not just the larger. For Peierls/F.4 at
  σ_t·R ≥ 10 the floor is ≲ 0.001%, which is NOT resolved by RICH.
  All future rank-N closure work must report the signed-error stability
  table across ≥ 2 quadratures as standard practice.

## Files

- `derivations/diagnostics/diag_pca_sectors_rich_adaptive.py` (main)
- `derivations/diagnostics/diag_pca_sectors_hollow_sph.py` (E6 original,
  BASE-quadrature reference — frozen artifact)
- `/tmp/dc_*.log` (per-point and crossing-check logs — not committed)
- `/tmp/direction_c_grid_runner.py`, `/tmp/direction_c_crossing_check.py`
  (one-shot runners used for this session — not committed)

## Recommendations

1. Close Issue #121 with "falsified at RICH, ratios are L17/L19 crossing
   artifacts, no structural advantage over F.4".
2. Do NOT lift NotImplementedError guards on `boundary="white_rank2"`.
   F.4 remains production.
3. Promote `test_pca_uniform_scale_does_not_beat_f4_at_rich_anchor` to
   `tests/cp/test_peierls.py` as a canary regression test.
4. Land L17+L19 in the research log: RICH is NOT the final word on
   structural floors at σ_t·R ≥ 10. For any future closure claim, a
   ≥ 2-quadrature signed-error stability table is mandatory.
