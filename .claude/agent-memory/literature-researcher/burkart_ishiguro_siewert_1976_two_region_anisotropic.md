---
name: Burkart-Ishiguro-Siewert 1976 — two-region linearly-anisotropic Milne+critical
description: Full PDF read (NSE 61, 72-77). Critical reflected slab + Milne with linear anisotropy. NOT F_N method — Case singular eigenfunctions + Chandrasekhar H/S functions. Tables I/II provide reference truth values.
type: reference
---

# Burkart-Ishiguro-Siewert 1976 NSE 61, 72-77 — what's actually in the paper

DOI 10.13182/NSE76-A28462. PDF in scratch/literature/, fully read.

## Verdict on method (CRITICAL framing fix)

This paper is **NOT F_N method**. It is the **Case singular-eigenfunction
method + Chandrasekhar H/S-function regularization**, applied to:
1. **Section II** — Milne problem for two adjoining half-spaces (each with
   its own c_α, b_α; α=1,2 means "core/reflector").
2. **Section III** — critical reflected slab (core c_1>1, reflector c_2<1,
   half-thickness `a`).

Linear anisotropy only: scattering kernel is `(1/2) c_α (1 + b_α μ μ')`,
NO P_2 / P_3 / P_L extension. This is a hard ceiling — the paper does
NOT generalize to general L.

## The equations (verbatim transcription)

### Transport equation (Eq. 1)
```
μ ∂Ψ_α/∂x + Ψ_α(x,μ) = (1/2) c_α ∫_{-1}^{1} Ψ_α(x,μ') (1 + b_α μ μ') dμ'
```
- `c_α` = mean secondaries per collision in region α
- `b_α` = linear-anisotropy coefficient
- `l_α := b_α (1 - c_α)` (notation, used everywhere below)

### Dispersion function (Eq. 5)
```
Λ_α(z) = 1 - c_α z R_α(z,z) tanh^{-1}(1/z) + c_α l_α z²
```
where `R_α(x,y) := 1 + l_α x y` (Eq. 4). Discrete eigenvalue ξ_α ∉ (-1,1)
is the positive zero of Λ_α. Notation: `ξ_1 = ν_0` (core), `ξ_2 = η_0` (reflector).

### Continuum eigenfunctions (Eqs. 4, 6, 7)
```
φ_α(ξ_α, μ) = (c_α ξ_α / 2) R_α(ξ_α, μ) [1/(ξ_α - μ)]    (4) discrete
φ_α(ξ, μ)   = (c_α ξ / 2)  R_α(ξ, μ) P/(ξ-μ) + λ_α(ξ) δ(ξ-μ)   (6) continuum
λ_α(ξ)      = 1 - c_α ξ R_α(ξ,ξ) tanh^{-1}ξ + c_α l_α ξ²      (7)
```
P = principal value.

### S-function (the regularization kernel)
```
S_2(μ',μ) = (c_2 μ μ')/(μ + μ') [1 - ĉ_2 (μ + μ') - l_2 μ μ'] H_2(μ') H_2(μ)
```
(Eq. 10). H_2 is **Chandrasekhar's H-function for region 2**. Convention:
```
ĉ_α := c_α l_α α_{α,1} / (2 - c_α α_{α,0})
q̂_α := 2(1 - c_α) / (2 - c_α α_{α,0})
α_{α,β} := ∫_0^1 H_α(μ) μ^β dμ      (Eq. 11)
```

### Critical condition (Section III)
Two coupled regular integral equations (Eqs. 36, 38) in unknowns
`Â(ν_0)` (discrete coefficient, core) and `D(ν) = Â(ν) e^{a/ν} e^{a/ν_0} / Â(ν_0)`
(Eq. 37). Solved iteratively. Two analytical approximations available:

- **Approximation B** (set D(ν)=0, ignore Eq. 38): yields
  ```
  a_{0B} = (1/2) π |ν_0| - z_{0B}              (39)
  ```
  with `z_{0B}` given by Eq. 40 (extrapolated endpoint, B-form).
- **Approximation M** (drop integral terms only): yields
  ```
  a_{0M} = (1/2) π |ν_0| - z_{0M}              (42)
  ```
  with z_{0M} given by Eq. 43.

Both reduce to identical when `b_1=b_2=0` (isotropic) **OR** when `c_2=0`
(bare reactor — single-medium Milne with anisotropy).

## Reference truth values (CRITICAL for V&V)

### TABLE I — extrapolated endpoints `z_{0B}, z_{0M}` (Milne)
Indexed by `(c_1, b_1, c_2, b_2)` with c_1∈{1.01, 1.06, 1.20, 1.50, 1.60},
b_1∈{0, 1}, c_2∈{0.4, 0.9}, b_2∈{0, 1}. **40 rows total.** Sample:
```
c_1=1.20, b_1=0.0, c_2=0.4, b_2=0.0:  z_0B = z_0M = 0.6992548
c_1=1.20, b_1=1.0, c_2=0.4, b_2=1.0:  z_0B = 1.037642, z_0M = 1.040783
```

### TABLE II — critical half-thicknesses `a_{0B}, a_{0M}, a_{exact}` (reflected slab)
Same indexing. **40 rows total.** Sample:
```
c_1=1.20, b_1=0.0, c_2=0.4, b_2=0.0:  a_0B=1.182975, a_0M=1.182479, a_exact=1.182419
c_1=1.50, b_1=1.0, c_2=0.4, b_2=1.0:  a_0B=0.5183184, a_0M=0.5167079, a_exact=0.5051470
```
Author claim (p. 77): exact `a` accurate to 3 sig figs for c_1≤1.2, only
1-2 sig figs for c_1≥1.5 (iterative solution truncation, not method limit).
**Match against Caroll-Aronson NSE 51 (1973)** confirmed in text.

## Mapping to Sood LA-13511 cases

The c_α∈(0,1) reflector is **NOT** the canonical Sood family (which uses
α=0 vacuum BC for bare configs and explicit reflector materials by name —
Pu/H2O/Fe — for reflected configs with TABULATED MATERIAL DATA, not c_α).

**Direct unlocks (anisotropic):**
- `*-1-1-SL` family (linearly-anisotropic slab, b≠0) — TABLE I provides
  Milne extrapolated endpoints (40 anisotropic Milne values).
- BIS-76 does NOT cover sphere or cylinder anisotropy.

**Indirect unlocks (two-region slab):**
- Reflected-slab family: TABLE II's 40 (c_1, b_1, c_2, b_2, a_exact)
  rows are independent reference values for any two-region slab
  reflected-critical solver. Slab-only.

**Estimate**: ~5-10 Sood slab cases get a second independent truth source
(BIS-76 vs Sood LA-13511 cross-check). The bulk value of BIS-76 is the
**Milne table** (Sood doesn't tabulate Milne, only critical) — opens a
NEW V&V axis (linear-anisotropy Milne) that Sood doesn't cover.

## Implementation feasibility — HEAVY (re-scoped)

**SymPy walls (severe):**
- Chandrasekhar H_2(μ) for **anisotropic** region — NOT closed form.
  Defined implicitly via nonlinear integral equation
  `H_α(μ) = 1 + c_α μ H_α(μ) ∫_0^1 [Ψ_α(μ')/(μ+μ')] H_α(μ') dμ'`
  with Ψ_α(μ) = (1/2)[1 + (something with l_α)]. **mpmath fallback mandatory.**
- α_{α,β} moment integrals (Eq. 11) — only available numerically.
- Coupled iterative solution of Eqs. 36 + 38 — Nyström quadrature.

**Architectural seams (per existing fn_method/ tree):**
- BIS-76 is **NOT F_N**. If integrating, must NOT live under `fn_method/`.
  Suggested package: `case_h_function/` or `singular_eigenfunction/`.
- Linear-anisotropy slab Milne is one solver; reflected-slab critical is
  a second solver sharing the Milne primitives. Factor as:
  ```
  case_h/
    core/dispersion.py        # Λ_α(z), eigenvalue ξ_α
    core/h_function.py        # iterative H_α(μ) solver
    core/s_function.py        # S_α(μ',μ) kernel
    milne/two_halfspace.py    # Section II
    critical/reflected_slab.py # Section III
  ```
- Reusing `fn_method/core/dispersion.py` is **structurally wrong**: F_N
  uses Marshak-moment matrix collocation; BIS-76 uses Case eigenfunction
  expansion. Different mathematical machinery.

## Errata check

Reviewed Eqs. 5, 7, 10, 14, 17-24, 36-43 by inspection. **Two minor
ambiguities** (not errata):
- Eq. 14 `W(ξ)`: `R_2` argument is `(ξ,ξ)`. Consistent with Eq. 4 if
  Region-2 functions are evaluated at Region-1 eigenvalues — confirm
  numerically before implementing.
- Eq. 35 sign of `ν(l_1 - l_2) q̂_2`: cross-medium term. Worth re-deriving
  symbolically.

No outright sign or factor errors found.

## Supplementary papers BIS-76 cites (acquisition queue)

- **McCormick & Kuščer, J. Math. Phys. 6, 1939 (1965)** — half-range
  orthogonality (Ref. 10). REQUIRED for the regularization step.
- **Pahor & Zweifel, J. Math. Phys. 10, 581 (1969)** — Case+Chandrasekhar
  coupling (Ref. 6). REQUIRED for the H-function machinery.
- **Caroll & Aronson, NSE 51, 166 (1973)** — independent reference
  values for cross-validation (Ref. 11).
- **Siewert & Burkart, NSE 58, 253 (1975)** — isotropic predecessor
  (Ref. 8). Should be read first — reduces to BIS-76 when b_1=b_2=0.
- **Burkart PhD Thesis, NCSU (1975)** — additional numerical results
  (Ref. 13). Possibly more tabulated values.

## Bottom line for ORPHEUS

BIS-76 is a **deep reference set** (80 tabulated values) but a **heavy
implementation**. If the goal is just unlocking Sood anisotropic slab
truth values, BIS-76's Table II provides ~40 independent reflected-slab
critical-thicknesses — value comes from the table, not from re-implementing
the method. Recommend: **use the published tables for V&V; do NOT
implement BIS-76 method until isotropic Case-eigenfunction infrastructure
exists** (Siewert-Burkart 1975 first).
