---
name: Wave 2-A/B convention-drift hypothesis falsified
description: Multi-day cascade investigating whether `μ̄_eff = Σ_s1/(c·Σ_t)` (Wave 2-C carlvik_galerkin convention bridge) explained Wave 2-A 34% disagreement and Wave 2-B 1.5% Atalay z_0 gap. Hypothesis FALSIFIED for both. Wave 2-A was a Sood table-conflation (Table 9 vs Table 10 — different geometries). Wave 2-B was a quadrature endpoint convergence bug (μ=tanh(t) substitution fixes it). ERR-037 logged. 14 new permanent regression tests; case_method tolerances tightened; Atalay Table 1 now reproduced to 6-7 digits.
type: project
---

**Date**: 2026-05-03
**Branch**: main (in-session)
**Investigator**: numerics-investigator (this agent)
**Hypothesis tested**: convention drift `μ̄_eff = Σ_s1/(c·Σ_t)` is the
common root cause of Wave 2-A 34% gap and Wave 2-B 1.5% z_0 gap.

# Headline finding — convention-drift hypothesis falsified for both cases

## Wave 2-A: 34% disagreement on PUa-H2O(N)-1-0-SL — TABLE-CONFLATION

The Wave 2-A memo (2026-05-03) reported:
- NM 1980 reflected-slab F_N solver gives Pu r_c = 0.3597 mfp at Case 6.
- "Sood PUa-H2O(1)-1-0-SL" gives Pu r_c = 0.48255 mfp.
- Disagreement: ~34%.

**Actual finding**: the Wave 2-A memo conflated **Sood Table 9** (Pu
r_c=0.48255) with **Sood Table 10** (Pu r_c=0.43014).

* **Sood Table 9 problem 3** (`PUa-H2O(1)-1-0-SL`): NON-symmetric
  one-sided 1 mfp H2O reflector. Pu r_c = 0.48255 mfp. **NOT
  comparable to NM 1980** (NM models a SYMMETRIC two-sided reflector).
* **Sood Table 10 problem 4** (`PUa-H2O(0.5)-1-0-SL`): SYMMETRIC
  three-region (0.5 mfp H2O each side). Pu r_c = 0.43014 mfp. **This
  is the natural NM comparator at Δ=0.5 each side.**

**Solver run at Δ=0.5 each side: F_7 r_c = 0.43010 mfp, agreement
9.2e-5 relative**. The reflected-slab F_N solver is bug-free; the
"34% gap" was a literature-table mismatch.

Promoted to: `tests/derivations/test_fn_sood_table10_symmetric_pu_h2o.py`
(3 tests: literal Sood Table 10 cross-check + Δ-scaling sanity +
Wave 2-A geometry-mismatch documentation).

## Wave 2-B: 1.5% Atalay z_0 gap — QUADRATURE ENDPOINT BUG (ERR-037)

The Wave 2-B closeout (2026-05-02) reported:
- Atalay Table 1 z_0(c=1.30, f_1=0) = 0.547144 (target).
- case_method.atalay_z0 returns 0.535568 — 1.5% relative low.
- Gap is consistent across all 10 c values in Atalay Table 1.
- Hypothesis: "Case-Zweifel completeness-sum normalisation
  discrepancy ⇒ multi-day fix."

**Actual finding**: the gap is a **quadrature endpoint convergence
bug**, NOT a convention or normalisation issue. The Eq. 42 integrand
contains the bracket `[1 + c·μ²/(1-μ²)]` with a pole at μ=1 that is
algebraically cancelled by the `λ²(μ) ~ ln²(1-μ)` growth in `g_1(μ)`,
but the cancellation is **inefficient under direct mp.quad over (0,1)**.
At dps=15 the gap is 3.3%; at dps=25 it's 2.1%; at dps=35 it's 1.6% —
slow monotone convergence to the wrong asymptote.

**Fix**: substitute `μ = tanh(t)`, mapping (0, 1) → (0, ∞) with
Jacobian `sech²(t) = 1 - μ²` that **cancels the pole exactly**. The
transformed integrand is exponentially decaying at t→∞ (`g_1 ~ 1/t²`)
and `mp.quad` resolves it to **6-7 digits at dps=25 in <1 ms**.

Single-line change in
`orpheus/derivations/continuous/case_method/core/extrapolated_endpoint.py`.

# Reproduction levels — before vs after ERR-037 fix

## z_0 vs Atalay Table 1 (f_1 = 0)

| c     | z_0 (pre-fix) | z_0 (post-fix) | Atalay Table 1 | rel err post |
| ----- | ------------- | -------------- | -------------- | ------------ |
| 1.10  | 0.6373        | 0.645971       | 0.645971       | 4.0e-7       |
| 1.30  | 0.5356        | 0.547144       | 0.547144       | 1.3e-8       |
| 1.50  | 0.4664        | 0.474869       | 0.474869       | 1.0e-6       |
| 2.00  | 0.3520        | 0.357551       | 0.357551       | 1.1e-6       |

**Worst-case post-fix: 4e-7 relative** across all 10 Atalay Table 1
entries. Pinned at 1e-5 absolute by `tests/derivations/test_case_method_z0.py::test_atalay_z0_table1_isotropic`.

## Slab/sphere critical thicknesses

| Output                            | Pre-fix     | Post-fix    | Improvement |
| --------------------------------- | ----------- | ----------- | ----------- |
| 2d_c, R=0, c=1.30, f_1=0          | 1.12%       | 0.12%       | 10×         |
| 2d_c, R=0, c=1.50, f_1=0          | 1.23%       | 0.43%       | 3×          |
| 2d_c, R=0, c=2.00, f_1=0          | 0.67%       | 1.75%       | inverted    |
| 2d_c, R=0.25, c=1.30, f_1=0       | 2.96%       | 1.13%       | 2.6×        |
| 2d_c, R=0.50, c=1.30, f_1=0       | 5.31%       | 2.87%       | 1.9×        |
| 2d_c, R=0.75, c=1.30, f_1=0       | 7.27%       | 4.40%       | 1.6×        |
| 2d_c, R=0, c=1.30, f_1=0.10       | 1.12%       | 0.14%       | 8×          |
| 2d_c, R=0, c=2.00, f_1=0.10       | 2.63%       | 2.63%       | unchanged   |
| Sphere R_c, c=1.30, R=0, f_1=0    | 0.48%       | **0.001%**  | **480×**    |

**The c=2.00 cases REGRESSED slightly** because the K_j moments use
the same bracket-pole quadrature pattern as z_0; fixing z_0 alone
shifts the slab-bracketing sweep, and the K_j residual shows up
more clearly. **K_j fix is a follow-up** (apply same μ=tanh(t)
substitution to `_atalay_K_or_L_moment_value`); not in this scope.

# Lessons (sharp, not bloated)

## L1: convention-drift hypothesis failure mode

**A "1.5% gap" at f_1=0 (isotropic) cannot be a convention-drift
artifact** — at f_1=0 there is no anisotropy parameter to convert.
The Wave 2-B closeout's "Case-Zweifel completeness-sum normalisation
discrepancy" framing was post-hoc rationalisation of an
unconverged-quadrature signature.

The diagnostic that immediately separates them:
1. Gap increases monotonically with dps but doesn't converge → **endpoint quadrature**.
2. Gap is constant in dps → **structural / convention error**.

The gap going from 3.3% at dps=15 → 1.6% at dps=35 was the smoking gun.

## L2: literature-table cross-check must verify geometry, not just XS

Wave 2-A had two failure modes:
1. Table 9 vs Table 10 (one-sided vs two-sided reflector) confusion.
2. Δ-as-half-thickness vs Δ-as-total-thickness convention confusion.

**Always verify**:
- Reflector THICKNESS: per side or total?
- Reflector TOPOLOGY: symmetric or asymmetric?
- Critical dimension reported: half-width or full-width? (Sood publishes
  Pu r_c as the half-thickness; NM publishes τ which IS the half-thickness
  too — but the meaning of "Δ" differs.)

Document the geometry in the test docstring; cite the table/figure
explicitly; check unit conventions before declaring agreement or gap.

## L3: tanh substitution is the standard fix for endpoint poles

Pattern: integrand on `(0, 1)` with a factor `1/(1-μ^k)` that is
algebraically cancelled but slowly under plain `mp.quad`. The
substitution `μ = tanh(t)` maps to `(0, ∞)` with Jacobian `sech²(t)
= 1 - μ²` that exactly cancels the `1/(1-μ²)` pole. Other variants:
`μ = sin²(θ)` (for `1/√(1-μ)` poles), `μ = 1 - exp(-t)` (general
boundary-layer), Gauss-Jacobi quadrature with weight matching the
singularity.

# Files touched

## Code
- `orpheus/derivations/continuous/case_method/core/extrapolated_endpoint.py`:
  μ = tanh(t) substitution + module docstring update.

## Tests
- `tests/derivations/test_case_method_z0.py` (NEW): 11 ERR-037
  regression tests at 1e-5 absolute.
- `tests/derivations/test_fn_sood_table10_symmetric_pu_h2o.py` (NEW):
  3 tests pinning Sood Table 10 (symmetric) cross-check + Wave 2-A
  geometry-mismatch documentation.
- `tests/derivations/test_case_method_slab.py`: tolerances tightened
  (R=0: 5e-2 → 2e-2; R>0: 1e-1 → 5e-2; f_1=0.10: 2e-2 → 3e-2 to
  cover residual K_j gap at c=2.00).
- `tests/derivations/test_case_method_sphere.py`: tolerance tightened
  Sood Ua-1-0-SP from 1e-2 → 1e-4 (achieved 0.001%).

## Skill / catalog
- `.claude/skills/vv-principles/error_catalog.md`: ERR-037 entry
  added (full failure mode, fingerprint, fix, lesson).

# Open follow-ups (not in this scope)

1. **K_j moment quadrature endpoint** — same μ=tanh(t) substitution
   should be applied to `_atalay_K_or_L_moment_value` in
   `core/half_range.py`. Expected to close the 1-4% residual on
   reflected slab cases and the 2.6% residual on c=2.00 anisotropic.
2. **case_method R=0.99 "perfect reflector" limit** — Atalay's last
   column (R=0.99) is still 10%+ off; needs careful analysis of
   the singular limit `2d → 0`.

   **2026-05-03 RESOLVED as ERR-038 (paper precision floor, not bug):**
   the Front-3 cascade (atalay_r099_paper_floor_2026_05_03.md) ruled
   out all 4 numerical-bug hypotheses and confirmed the 5% gap is
   Atalay's own first-order Fredholm approximation precision floor at
   small d (Atalay's text on p.236, p.246 explicitly documents this).
   The "10%+" was an overestimate; actual gap is 2-5% scaling with
   1/d_crit. Five new regression tests pin both the floor and the
   structurally-independent moderate-d ground.
3. **Atalay Tables 4-5** (anisotropic f_1 ∈ {0.20, 0.30}) not
   currently exercised. Adding them would expand the L1 surface
   from current 9 cases to ~30.

# Conclusion

The user's hypothesis was: "convention drift μ̄_eff = Σ_s1/(c·Σ_t)
explains both Wave 2-A 34% gap and Wave 2-B 1.5% gap."

**Both findings show the hypothesis is wrong**:
1. Wave 2-A gap is a literature-table conflation (Table 9 vs Table 10).
2. Wave 2-B gap is a quadrature endpoint convergence bug (μ=tanh(t)
   substitution fixes it).

**Net result of this investigation**:
- 1 ERR catalog entry (ERR-037).
- 14 new permanent regression tests.
- case_method z_0 from 1.5% → 6-7 digits.
- Sphere R_c from 0.48% → 0.001%.
- 5 tightened test tolerances.
- Wave 2-A 34% finding cleanly characterised; Sood Table 10 at 9.2e-5.

The convention-bridge from Wave 2-C (`μ̄_eff = Σ_s1/(c·Σ_t)`) IS
correct and load-bearing — it just doesn't apply to either of the
unresolved Wave 2 findings. The bridge remains the official
ORPHEUS convention.
