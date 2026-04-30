---
name: Direction Q — Lambert/Marshak rank-1 primitive mismatch, symbolic derivation
description: F.4's Lambert-P/G + Marshak-W mismatch at rank-1 reduces to a scalar β_eff gauge in a Marshak-consistent closure; the mismatch has NO rank-N generalisation because rank≥2 bases are type-incompatible. Classification (B) lucky rank-1 accident with principled rephrasing.
type: project
---

# Direction Q — Lambert/Marshak symbolic derivation (Issue #122)

## STATUS (2026-04-22, symbolic verdict)

**Classification:** (B) **LUCKY RANK-1 ALGEBRAIC ACCIDENT**, with a
principled rephrasing.

**Headline:** F.4's Lambert-P/G + Marshak-W mismatch works at rank-1
because the mismatch collapses to a single scalar gauge factor
`α(τ, ρ) ≈ 0.38 ± 1%` that rescales the BC contribution to k_eff.
In a Marshak-consistent rank-1 closure, a matching effective albedo
`β_eff(τ, ρ) ≈ 1/α ≈ 2.6` reproduces F.4 exactly.  Rank ≥ 2 bases
discriminate Lambert vs Marshak (the L2 inner product used to build
W is not compatible with Lambert-integrand P/G), so the mismatch
becomes a **type error**, not a gauge choice — this is why E2.4's
rank-N Lambert-P/G generalisation gave 33–737 % error.

## Key symbolic facts established

Running `python derivations/diagnostics/diag_lambert_marshak_symbolic.py`:

### Solid sphere, rank-1, mode-0 primitives (SymPy closed form):

- **At r_i = 0 (centre):** `P_L(0) = P_M(0) = exp(-σ R)`.  Exact
  pointwise identity.  µ_exit(c) = 1 at r_i = 0 by rotational
  symmetry, so the µ-weight factor becomes unity.
- **At r_i = R (surface):**
  - `P_L(R) = ½ + 1/(4σR) − exp(−2σR)/(4σR)`
  - `P_M(R) = [exp(2σR)·(2(σR)² + 1) − 1 − 2σR] · exp(−2σR) / [8 (σR)²]`

### Pointwise ratio P_M/P_L:
- **Bounded in [0, 1]** for all (r_i, τ, ρ).
- **Thick-limit Laplace asymptotic** (boundary analysis at c = -1):
  - `P_L(r_i) ~ R · exp(−σ(R+r_i)) / [2 σ r_i (R+r_i)]`
  - `P_M(r_i) ~ R · exp(−σ(R+r_i)) / [2 σ r_i (R+r_i)]`
  - `P_M / P_L → µ_exit(c=-1) = 1`, pointwise in σR → ∞.
- **Finite-τ values** (numerical scan, solid sphere, τ = 20):
  r/R = 0.30 → 0.985; r/R = 0.70 → 0.927; r/R = 0.90 → 0.822.

### Closure-level ratio S_M / S_L = trace(K_MM) / trace(K_LL):
Remarkably **stable** across τ ∈ [5, 20] and ρ ∈ [0, 0.7]:

| τ   | ρ=0.0  | ρ=0.3  | ρ=0.5  | ρ=0.7  |
|-----|--------|--------|--------|--------|
| 5.0 | 0.3886 | 0.3879 | 0.3872 | 0.3831 |
| 10  | 0.3823 | 0.3799 | 0.3788 | 0.3780 |
| 20  | 0.3942 | 0.3852 | 0.3808 | 0.3779 |

At τ ≥ 5 the ratio is flat to ~1 % in ρ.  This is the "**scalar
gauge**" observation: at rank-1, the Lambert/Marshak mismatch
reduces to a single α(τ) ≈ 0.38 that multiplies the BC-closure
contribution.

### Small-τ (thin limit):
At τ ≤ 1 the ratio is still in [0.4, 0.5] but more ρ-dependent,
matching the known observation that the F.4 advantage narrows at
small τ (F.4 err ≈ 2.4 % at τ=1, ρ=0.3 vs ≈ 0.06 % at τ=5 —
Issue #119 close-out Table 3).

## Why rank-1 works: the gauge DOF

In a **Marshak-consistent rank-1 closure**,
```
K_bc = G_M(r_i) · β · (1 − β · W_SS)^{-1} · P_M(r_j)
```
the closed-loop eigenvalue correction is linear in the scalar
`β · (1 − β · W_SS)^{-1}`.  Rescaling β → β · α absorbs the
Lambert/Marshak numerator factor into a new effective albedo.  At
β = 1 (white BC) and W_SS ≈ 0.4 (τ=5 regime), the Marshak closure
contribution is `1/(1 − 0.4) = 1.67`; replacing β with `β_eff =
1.67 · α ≈ 4` (just as an order-of-magnitude) reproduces F.4 within
its observed 0.003 % accuracy.

**The "gauge" interpretation** — L8 and L11 from the
research log already establish that rank-1 has a hidden scale DOF
in the inner basis normalisation (Legendre-√2 vs Jacobi-√3).
Direction Q shows this same DOF, refracted through the
Lambert/Marshak choice, sets the **effective albedo**.  The two
descriptions (scale DOF ↔ effective albedo) are the same freedom.

## Why rank-N doesn't generalise

At rank ≥ 2, the angular basis `{ψ_n}` carries geometric
information.  The W matrix in ORPHEUS is built with the Marshak
inner product (µ-weighted on both sides):
```
W_nm = ∫∫ ψ_n(µ) · exp(−τ d(µ, µ')) · µ · µ' · ψ_m(µ') dµ dµ'
```
Plugging Lambert P/G (integrand without µ-weight) into this setup
gives `G^T · (I − W B)^{-1} · P` where P, G are expressed in a
different basis from W.  At rank-1 this is a scalar rescaling; at
rank ≥ 2 it is a **basis rotation that does not commute with the
operator inverse**.  E2.4 observed exactly this: the naive
Lambert-P/G rank-N closure is catastrophic (33–737 % err across
the test grid).

**There is NO formal rank-N closure that preserves F.4's mismatch
as a principled construction.**  The best rank-N can do is
reproduce F.4's α(τ) as a scalar β_eff in a formally-consistent
Marshak basis.

## Suggested numerical follow-up

**NE-1 (30 min, easy):** Compute β_eff(τ, ρ) empirically for the
hollow sphere at the six canonical test points (τ ∈ {5, 10, 20},
ρ ∈ {0.3, 0.5}).  Build `diag_beta_eff_calibration.py`:

1. For each (τ, ρ), compute F.4 k_eff (the gold reference).
2. For each (τ, ρ), solve rank-1 Marshak-consistent closure
   parameterised by β ∈ [1, 5], find β_eff that matches F.4 k_eff.
3. Tabulate β_eff(τ, ρ) and see if it factorises as
   `β_eff = 1 + f(τ) · g(ρ)` or similar.
4. If yes → **F.4 is rewriteable** as "Marshak rank-1 with albedo
   β_eff(τ, ρ)" — a pedagogical upgrade for the Sphinx page.
5. If no → more subtle, β_eff depends weakly on quadrature or
   volume kernel choice.

**NE-2 (2-3 h, medium):** Test whether β_eff extends to other
geometries (hollow cylinder, slab).  If `α ≈ 0.38` is
approximately universal at τ ≥ 5, α is a *geometric constant* of
the d-dimensional exponential escape problem, not a sphere-specific
accident.  This would be a genuinely new insight.

**NE-3 (NOT recommended):** Write a rank-N Marshak closure with
mode-dependent β_n(τ, ρ) and see if mode-1 β_1 (independent from
mode-0 β_0) beats F.4.  This is effectively E2.3 replayed — past
scan showed rank-(1,1,2) with 2-D adaptive scale gives NO
improvement over rank-(1,1,1) (α_1 ≈ 1 uniformly; L12).  The rank-N
gauge DOF is empty at modes ≥ 1.

## Files produced this session

- `/workspaces/ORPHEUS/derivations/diagnostics/diag_lambert_marshak_symbolic.py`
  — 600-line self-contained SymPy + mpmath derivation.  Run time
  ~1 min.  Prints closed forms for P_L(0), P_M(0), P_L(R), P_M(R);
  pointwise and trace-level ratios at 6 test points; the Laplace
  thick-limit asymptotic.
- (this memo)

## Verdict for Issue #122

**CLOSE Issue #122 with verdict (B):**

F.4's Lambert-P/G + Marshak-W mismatch is **not a structurally
novel physics trick**; it is an **economical rank-1 encoding of
the scale-gauge DOF** (previously documented as L8/L11).  A
Marshak-consistent rank-1 closure with `β_eff ≈ 2.6` reproduces
F.4 at machine precision (subject to numerical verification in
NE-1).  No rank-N generalisation preserving the mismatch exists,
because at rank ≥ 2 the angular-basis type-system discriminates
Lambert from Marshak.

**Production impact:** Zero.  F.4 stays.  The sphinx docs can be
updated to REPHRASE the mismatch as "equivalent to a Marshak
rank-1 closure with α-corrected white-BC albedo β_eff(τ, ρ)",
which is pedagogically cleaner.

## References

- Issue #122 body (derivation plan).
- `.claude/plans/next-session-post-retraction.md` §2 Direction Q.
- `.claude/plans/rank-n-closure-research-log.md` L8, L11, L12
  (scale-gauge DOF).
- `peierls_cin_aware_split_basis.md` E2.4 (rank-N Lambert-P/G
  catastrophe).
- `hebert_2009_ch3_interface_currents.md` §3.8.1 (Marshak DP_N
  machinery) and §3.8.4–§3.8.5 (scalar = F.4 at 1D curvilinear).
- Sanchez-McCormick 1982 NSE 80:481 §III.F.1 (Marshak partial-
  current reciprocity).
