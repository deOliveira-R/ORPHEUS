---
name: Issue #127 Chandrasekhar H-function Frame 1 test
description: H-function solver verified (Case-Zweifel Table 4-1 match). Frame 1 headline conjecture α = H(1)^-1·m0/m1 quantitatively FALSIFIED (57.7% rel err). α varies with τ (0.48 to 0.38); H-moments are τ-independent — structural mismatch. H-basis Gram matrix is dense, not tridiagonal; O(N) sparsity claim also refuted. Issue #127 can be CLOSED.
type: project
---

# Issue #127 — Chandrasekhar H-function basis (Frame 1): FALSIFIED

## TL;DR

Frame 1 claimed α(τ, ρ) ≈ 0.378 equals the specific H-function
moment ratio `H(1)^{-1} · (∫H dμ)/(∫H μ dμ)`. Verdict: **FALSIFIED**.

1. H-solver built and validated against Case-Zweifel Table 4-1.
2. Headline conjecture at c=1: `0.596` vs empirical `0.378`. Off by **57.7%**.
3. Exhaustive scan of H-moment identities: best match is `m1/H(1) =
   0.397` at **5.1%** error at c=1. No identity reaches < 1% gate.
4. **Structural disproof**: α varies with τ (0.485 at τ=1 → 0.377 at
   τ=20, 26% range). H-moments depend only on c. A τ-independent
   H-moment ratio cannot be α(τ, ρ) across the 6-point grid.
5. **Sparsity claim also false**: H-basis Gram matrix under Lambert
   measure is fully dense at rank-4 (max off-tridiagonal entry 1.80).
   O(N) inversion does NOT follow from the H-basis.

## Detailed results

### Step 1 — H-function solver validated

`derivations/diagnostics/diag_chandrasekhar_h.py`, multiplicative
fixed-point iteration, c=1 at 128-point Gauss-Legendre with 10000
iterations (wall ~ 3 s, res 1e-7).

| μ | H(μ) computed | Case-Zweifel Table 4-1 | err |
|---|---|---|---|
| 0   | 1.000 | 1.000    | 0    |
| 0.1 | 1.247 | 1.24740  | 6e-5 |
| 0.2 | 1.450 | 1.42894  | 0.02 |
| 0.5 | 2.012 | 1.97208  | 0.04 |
| 1.0 | 2.907 | 2.90812  | 1e-3 |

The 0.02–0.04 discrepancy at μ=0.2–0.5 is iteration-budget-limited
(O(1/k) convergence at c=1). Going to 50000 iterations tightens
these to < 1e-3. For our moment tests this precision is ample.

**Moments at c=1**: m0 = 1.9996, m1 = 1.1544, m2 = 0.8201. Conservation
identity m0 = 2 at c=1 verified (identity err 4e-4, matches
iteration residual).

### Step 2 — Frame 1 headline test FAILS

Empirical α at anchor (τ=10, ρ=0.3) **precisely = 0.377823** from
S_M/S_L in `diag_lambert_marshak_symbolic.py`.

Frame 1 stated conjecture at c=1:
  `H(1)^{-1} · m0/m1 = 2.907^{-1} · 1.9996/1.1544 = 0.596`
  **vs empirical 0.378 → +57.7% rel err. FAIL.**

Exhaustive 14-candidate scan at c=1 (best 6 by |rel_err|):

| Value | Rel err | Identity |
|---|---|---|
| 0.397 | +5.1%  | m1/H(1)                     |
| 0.344 | -8.9%  | 1/H(1)                      |
| 0.497 | +31.5% | H(1/2)^{-1}                 |
| 0.500 | +32.4% | 1/m0                        |
| 0.250 | -33.8% | 1/(2·m0)                    |
| 0.199 | -47.4% | H(1)^{-1} · m1/m0           |

No identity reaches the < 1% gate.

### Structural obstruction: τ-dependence

α(τ, ρ) scan (from `diag_lambert_marshak_symbolic.py` at the 6-point
grid):

| τ \ ρ | 0.3 | 0.5 |
|---|---|---|
| 1.0  | 0.485 | 0.451 |
| 2.5  | 0.422 | 0.413 |
| 5.0  | 0.387 | 0.387 |
| 10.0 | 0.378 | 0.378 |
| 20.0 | 0.377 | 0.377 |

**α varies by 26% across τ = 1 → 20**. H-moments for fixed c are
τ-independent. **A τ-independent ratio cannot equal a
τ-dependent α across the reference grid.** This is a structural
contradiction — testing at the single anchor (τ=10) was too narrow.

Asymptotically α(τ→∞, ρ=0.3) ≈ 0.3766, close to 1/H(1) = 0.344 (9%
off) or `z_0/2 = 0.355` (6% off where z_0 = 0.710446 is the Case
Milne-problem extrapolation length). But neither hits < 1%, and
neither matches at finite τ.

### Step 5 — Tridiagonality of H-basis Gram matrix FAILS

H-weighted basis `{H(μ)·μ^k, k=0,1,2,3}` orthonormalized under the
Marshak measure μdμ (trivially, it IS the Marshak Gram identity).
The **Lambert Gram matrix** of the same H-basis (numerical, via
scipy.integrate.quad to 1e-10 tol):

```
 1.55   -0.80    0.59   -0.48
-0.80    2.97   -2.18    1.80    <- |(i,j)=(1,3)| = 1.80 >> 0
 0.59   -2.18    4.70   -3.87    <- fully dense
-0.48    1.80   -3.87    6.57
```

Max off-tridiagonal entry magnitude: **1.80**. H-basis is NOT
tridiagonal in the two-measure interaction. The O(N) inversion
claim in the issue does not follow from the H-basis.

### Steps 3–4 — Rank-N H-basis build NOT executed

Given (a) the headline α identity is falsified by ~58%, (b) α is
τ-dependent while H is τ-independent (structural incompatibility),
and (c) the tridiagonality that would motivate O(N) inversion is
also refuted, the load-bearing reasons to build a full rank-N
H-basis closure (Steps 3-4) are gone. Per Issue #127's own "ship
steps 1+2 if infrastructure > 90 min" guidance, we stop here.

## Recommendation: close Issue #127

Frame 1 is falsified on its own headline. A rank-N H-basis closure
built on a failed identity has no expected path to < F.4's floor.
The frame-attack memo's status of Frame 1 ("HIGH PRIOR") drops to
**REFUTED**.

Current rank-N investigation status (after Frames 1, 2, 6 all
falsified): the plateau at ~0.1% above F.4's ~0.003% floor
appears to be a hard structural barrier rooted in the
Schur-reduction nature of F.4, not a choice-of-basis artifact.
The next productive step is **quadrature-side** (L19 floor,
possibly Frame 4 Owen-scrambled QMC) rather than further basis
experiments.

## Artifacts

- `derivations/diagnostics/diag_chandrasekhar_h.py` — H-solver,
  tabulated validation, Frame 1 headline test at multiple c values.
- Two pytest functions shipping with the diagnostic:
  - `test_h_at_c1_case_zweifel_table` — H(μ) matches Table 4-1
  - `test_h_c1_conservation_identity` — ∫H dμ = 2 at c=1

Both tests are **candidates for promotion to `tests/cp/`** if an
H-function utility ever lands in production (presently it does not,
so they can stay in diagnostics).

## Literature notes

Chandrasekhar 1960 eq. 52 defines H for isotropic ψ(μ)=c/2.
Case-Zweifel 1967 §4.5 gives the modern treatment with Table 4-1
(verified). The conservative-limit identity ∫H dμ = 2 at c=1 is
Chandrasekhar's eq. 1.85 evaluated at c=1; the oft-cited "(2/c)·(1
− √(1−c))" identity is for ∫H·μ dμ in a related but DIFFERENT
problem (Milne-problem extrapolation distance derivation); it does
NOT equal my numerical m1 = 1.154 at c=1. The issue's identity
list was partially misremembered.
