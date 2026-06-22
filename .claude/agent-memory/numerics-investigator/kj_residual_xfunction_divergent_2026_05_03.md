---
name: K_j residual gap is X-function divergent integrand, NOT Signature 7
description: Wave 2-B follow-up to apply ERR-037's μ=tanh(t) substitution from z_0 to K_j moment quadrature. Phase 1 fingerprint FALSIFIED the hypothesis. K_j outer scipy.quad converges fine (~1e-11). Bottleneck is X(-ν) which has a logarithmically-divergent integrand (Atalay Eq 40 as transcribed). μ=tanh(t) substitution does NOT close the gap (4.40% → 4.44%, marginal degradation). The X-function gives finite values via algorithm-dependent regularisation; production scipy regularisation differs from Atalay's implicit regularisation. Open as separate investigation.
type: project
---

**Date**: 2026-05-03
**Branch**: main (in-session)
**Investigator**: numerics-investigator (this agent)
**Hypothesis tested**: extending the μ=tanh(t) substitution that fixed
ERR-037 (z_0 evaluator) to `_atalay_K_or_L_moment_value` (K_j/L_j moment
quadrature) closes the residual 1-4 % gap on Atalay Table 2 reflected
slabs.

# Headline finding — hypothesis FALSIFIED

The K_j outer integrand does NOT have the same algebraic-pole structure
as z_0. Phase-1 convergence-with-quadrature-tolerance fingerprint shows:

- K_j scipy.quad converges to ~1e-11 between baseline (epsabs=1e-10) and
  tight (epsabs=1e-13). NOT Signature 7 at the K_j call site.
- Bottleneck is X(-ν) accuracy: scipy backend at default tolerance
  has ~1e-3 relative error vs mpmath dps=30, propagated into K_j as X².

# Why X(-ν) is hard — the Atalay Eq 40 integrand is divergent

Atalay's Eq 40 X-function:

```
X(μ) = exp{-c/2 ∫₀¹ g_1(c,ν) [d²(ν²)(1+cν²/(1-ν²)) + 3f₁(1-c)²ν²d(-ν²)] ln(ν-μ) dν}
```

The bracket carries a `1/(1-ν²)` simple pole at ν=1. The supposed
cancelling factor `g_1 = 1/[λ²(ν) + (πcν/2 d(ν²))²]` has λ ~ -c·atanh(ν)
~ -(c/2) log[2/(1-ν)] at ν → 1, so `g_1 ~ 4/(c² log²(1-ν))`. The product
`g_1 · bracket · log(ν-μ)` decays only as `~log(ν-μ)/[(1-ν) log²(1-ν)]`,
whose primitive is `∝ -1/log(1-ν)` — bounded but **non-zero** as
truncation point → endpoint.

Empirically (mpmath dps 15-60, X(-0.55) at c=1.30, f1=0):

```
dps=15: X = 0.86648
dps=25: X = 0.86349
dps=30: X = 0.86263
dps=40: X = 0.86154
dps=60: X = 0.86037
```

Drift ≈ 6e-3 across dps 15→60. This is **NOT** convergence — this is the
algorithm sampling closer and closer to the divergent endpoint at higher
dps, integrating a wider portion of the divergent tail.

The `μ=tanh(t)` substitution maps `(0,1)→(0,∞)` exactly:

```
∫₀¹ g_1·bracket·log(ν-μ) dν = ∫₀^∞ g_1·[d²(sech² + cν²) + 3f₁(1-c)²ν²d(-ν²)sech²]·log(ν-μ) dt
```

Bit-perfect substitution at intermediate ν=0.7 confirmed (diff=0.0).
But the substituted form is also divergent — `g_1 ~ 1/t²` and
`bracket·dν ~ d²·c·1` (constant) at large t, so the integrand decays
only as `1/t²·log(ν-μ)` ~ const/t². ∫₁^T const/t² dt = const(1 - 1/T).
**Bounded** but every truncation T introduces O(1/T) error.

# What the production code does

`_X_integrand_real` in `core/x_function.py` integrates `(0,1)` directly
via scipy.quad with `points=[0.5]`. scipy.quad in float64 uses adaptive
Gauss-Kronrod and stops sampling when the contribution to the absolute
error falls below `epsabs=1e-12`. The "stopping point" is determined by
float-precision (last sample at ν ≈ 1 - O(ε_machine)), giving X(-0.55) =
0.86152 — close to the dps≈40 mpmath value, but NOT the same.

# What I tried, and why it didn't help

Implemented the μ=tanh(t) substitution in BOTH `_atalay_X_function_scipy`
and `_atalay_X_function_mpmath` over `(0, T_MAX=20)` and over
`[0, mp.inf]`. Net effect on Atalay Table 2:

| R    | pre-tanh-X gap | post-tanh-X gap |
| ---- | --------------:| ---------------:|
| 0.25 |          1.13% |           1.135% |
| 0.50 |          2.87% |           2.892% |
| 0.75 |          4.40% |           4.436% |

**Marginal degradation**, not improvement. Reverted.

The parallel investigation (r099 perfect-reflector front) had ALREADY
reached this conclusion earlier in the session — see
`derivations/diagnostics/diag_atalay_r099_06_pin_xfunction_fix.py`
which shows R=0.99 still 5.03% off after a tanh-X-function patch.

# What's actually wrong — provenance question

The integrand of Atalay Eq 40 as transcribed appears to be
logarithmically divergent. Either:

1. Atalay's Eq 40 has a missing `(1-ν²)` factor (typesetter loss).
2. The integral is interpreted as a Cauchy principal value or otherwise
   regularised (and Atalay's Table 2 is generated under that
   regularisation).
3. The closed-form Case-Plazcek-Hofmann 1961 X-function for isotropic
   should be substituted (see e.g. McCormick 1964 or Case-Zweifel 1967).

The X-function defined IMPLICITLY by Eq 26 (`X(μ) = ∫γ(ν)d(ν²)/(ν-μ)dν`)
has a different convergence behaviour. The Eq 40 representation as a
calculation formula must connect via Plemelj-Sokhotski-style
integration; the divergence is presumably absorbed in the
exponentiation, but only under a specific regularisation.

# Open as a separate investigation

This is **NOT** Signature 7. The probe-cascade label "extends ERR-037"
is INCORRECT for this investigation. To close the residual gap requires:

- (a) Find Atalay's actual regularisation (probably via direct Eq 26
  `X(μ) = ∫γd(ν²)/(ν-μ)dν` with explicit closed-form reduction of the
  divergent log-tail in the IMPLICIT exp(...) representation). Expensive.
- (b) Switch to closed-form Case-Plazcek-Hofmann X-function for
  isotropic c-cases. Very different code path, but possibly tractable.
- (c) Accept the residual 1-4 % gap as an analytical-vs-implementation
  regularisation mismatch. Document it in the test tolerance and ship.

Outside scope of this session per the task brief.

# Files touched (this investigation)

- `derivations/diagnostics/diag_kj_phase1_fingerprint.py` — K_j outer
  scipy.quad convergence fingerprint (rules out Signature 7 at K_j level).
- `derivations/diagnostics/diag_kj_phase1b_xfunction.py` — X(-ν)
  scipy-vs-mpmath agreement (proved the bottleneck is X, not K_j outer).
- `derivations/diagnostics/diag_kj_x_function_divergent_integrand.py`
  — pinned regression: X integrand drifts 6e-3 across dps 15→60
  (logarithmic divergence fingerprint), AND pinned residual K_j gap
  at the post-ERR-037 baseline values (1.13%, 2.87%, 4.40%).

# Lessons sharpened

## L (sharpening of existing wave2 lesson)

**Both Wave 2-A and Wave 2-B had a "looks like Signature 7" pattern at
first glance** (slow monotone convergence with dps), but the underlying
mechanism was different:

- Wave 2-A: literature-table conflation (Sood Table 9 vs Table 10).
- Wave 2-B z_0: **TRUE Signature 7** (algebraic pole cancelled by log²,
  μ=tanh(t) regularises). Fixed.
- K_j residual: NOT Signature 7. Algebraically-divergent integrand,
  algorithm-dependent regularisation. μ=tanh(t) doesn't help because
  the integral itself is divergent.

The Phase-1 verdict on convergence-with-tolerance fingerprint should
be: "scipy.quad converges fine at the OUTER level → look INSIDE the
integrand for sub-routines whose own convergence is the bottleneck."
Don't apply Signature 7's substitution recipe blindly to call sites
that look similar without verifying the SAME pole structure exists.

# Conclusion

The Wave 2-B follow-up scope ("apply μ=tanh(t) to K_j") is closed as
**hypothesis FALSIFIED**. Production code unchanged from baseline. The
residual ~1-4 % gap on Atalay Table 2 reflected slabs is pinned by
`derivations/diagnostics/diag_kj_x_function_divergent_integrand.py`
and tracked as an open investigation requiring deeper analytical work
on the X-function regularisation.
