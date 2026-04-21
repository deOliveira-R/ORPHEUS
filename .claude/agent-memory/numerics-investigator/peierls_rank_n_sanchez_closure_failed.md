---
name: Rank-N Sanchez-McCormick §III.F.1 recipe does NOT close to 0.1%
description: Systematic scan of 60+ variants of the Sanchez per-face recipe for hollow sphere rank-N: best is 1.43% at N=2 (plateau), not the 0.1% target. Neither the µ-weighted orthonormal basis with consistent µ-weights on P, G, W, nor F.4/Sanchez hybrid works.
type: project
---

# Rank-N hollow-sphere: Sanchez-McCormick recipe plateaus at 1.43%

## Investigation date
2026-04-21. Target: ≤ 0.1 % k_eff residual at R=5, r_0/R=0.3, homogeneous Σ_t=1, Σ_s=0.5, νΣ_f=0.75 (k_inf=1.5). Current shipped N=1: 0.077 %. Prior N=2 best: 3.87 %.

## µ-weighted orthonormal basis (Gram-Schmidt on [0,1] under weight µ)

Symbolically verified at `derivations/diagnostics/derive_mu_weighted_basis.py`:

- f^0(µ) = √2
- f^1(µ) = 6µ − 4
- f^2(µ) = √6 · (10µ² − 12µ + 3)
- f^3(µ) = √2 · (70µ³ − 120µ² + 60µ − 8)

These are proportional to Jacobi P^{(0,1)}_n shifted to [0,1]. Promote to
`tests/derivations/` as a general mathematical utility test.

## Negative results — systematic scan (all on hollow sphere at target params)

### Recipe matrix (f^n basis × µ weights × Jacobian × Gelbard × Model A/B)

Best over ALL 60+ combinations:

| Basis | α_P | α_G | α_W | Jac | Gelbard | N=1 err | N=2 err | N=3 err |
|---|---|---|---|---|---|---|---|---|
| µ-ortho | 1 | 1 | µ²both | 0 | none | 2.55 % | **1.43 %** | 1.42 % |
| µ-ortho | 1 | 1 | 0 | 0 | none | 1.55 % | 1.60 % | — |
| µ-ortho | 0 | 1 | 0 | 0 | none | 13.3 % | 13.3 % | — |
| shifted Leg | 1 | 1 | 0 | 0 | none | 11.0 % | 10.9 % | — |
| shifted Leg | 1 | 1 | 0 | 0 | yes | 4.9 % | 4.9 % | — |
| F.4 mode 0 + µ-ortho n≥1 | — | — | — | — | — | **0.077 %** | **3.88 %** | 3.88 % |

The µ-ortho basis with consistent µ-weights plateaus at ~1.43 % for N ≥ 2
(tested to N=4). This is a self-consistent closure that solves a
slightly different equation — NOT the original physics problem.

### Findings

1. **N=1 bit-exactness with F.4 and Sanchez closure are MUTUALLY INCOMPATIBLE.**
   F.4's scalar mode-0 uses no µ weight (Lambert basis). Sanchez uses f^0=√2
   with µ weight. Any hybrid (F.4 mode 0 + Sanchez n≥1) gives bit-exact
   N=1 but DEGRADES at N=2 because the cross-mode W entries couple
   incompatible inner products.

2. **No constant prefactor fix works.** Scanning P and G prefactors in
   {0.5, 1.0, 2.0} with αW ∈ {0,1}: minimum N=2 error is 1.43 %, always
   with N=1 error of 2.54 %. The residual does NOT collapse to zero under
   any scalar rescaling.

3. **Gelbard (2n+1) factor cannot be placed.** Applied as D·(I-W·D)⁻¹,
   or D·(I-W)⁻¹·D, or D·(I-W)⁻¹, or (I-W)⁻¹·D: all WORSEN the N=2 error
   (5-40 %). Consistent with the memo's note that (2n+1) lives in the
   reflection only for the *unweighted* Legendre basis.

4. **The (ρ/R)² Jacobian is confirmed spurious.** Any recipe variant
   including (ρ/R)² gives 18-48 % error. Phase F.5 was right to drop it.

5. **Model A vs Model B split has negligible effect.** P=A/G=B gives 1.43 %
   vs P=A/G=A at 1.42 % — within quadrature noise. Not a structural
   variable for rank-N.

## Root cause analysis

The ORPHEUS F.4 scalar closure uses **Lambert angular-flux basis** for mode 0:
    basis = {f^0 = 1}, integrand measure = sin θ dθ (no µ)

Sanchez-McCormick §III.F.1 uses **µ-weighted partial-current basis** for mode 0:
    basis = {f^0 = (π A)⁻¹ or equivalently √2 / (π A)}, integrand measure = µ sin θ dθ

These are OBJECTIVELY DIFFERENT mathematical closures producing
different k_eff values. F.4 happens to give 0.077 % for this problem;
Sanchez mode-0 alone gives 2.55 %.

The memo Section 6c states: "The σ_t-dependent residual is the
fingerprint of a missing (Ω·n) = cos θ weight in the partial-current
moment." That hypothesis was tested and found incomplete — adding the
µ weight to Sanchez's primitives does give *improvement* (from 5.59 %
shipped to 1.43 % best), but does not close to 0.1 %.

**Probable interpretation**: for the observer-centered Peierls-Nyström
framework, the Sanchez per-face formalism is not the correct closure.
The correct closure must respect:
  - The Peierls kernel's fundamental observer-centered integration
  - The F.4 scalar mode-0 convention (0.077 %)
  - A rank-N extension that adds basis functions orthogonal to f^0 in
    the SAME inner product space

The closest candidate would be the plain shifted Legendre basis with
**no µ weight anywhere**, matching the (2n+1) factor in the reflection
(this is the unweighted-Lambert convention from memo §5). But that gives
10.9 % at N=2 per our scan — also doesn't close.

## Recipe NOT delivered — memo guidance was incomplete

The literature memo at `sanchez_mccormick_rank_n_per_face.md` was
correct about the structural fix (µ-weight in primitives, not `(ρ/R)²`
Jacobian), and the fix does reduce error from 5.59 % to 1.43 %. But
closing to 0.1 % requires MORE than the µ-weight fix — the memo's
"Bottom-line correction recipe" steps 1-4 produce the 1.43 % result, not
the ≤ 0.1 % the task specification demanded.

## Code changes

**None shipped.** `_build_closure_operator_rank_n_white` continues to
raise `NotImplementedError`. The shipped Phase F.5 `_marshak` variants
are preserved but unreachable. Diagnostic scripts committed as-is.

## Scripts (all in `derivations/diagnostics/`)

- `derive_mu_weighted_basis.py` — Gram-Schmidt derivation of µ-weighted
  orthonormal polynomials on [0,1]. Verified to 1e-12. **Promote to
  `tests/derivations/test_mu_weighted_basis.py` as general utility.**
- `diag_sanchez_recipe_scan.py` — 13-variant initial recipe scan.
- `diag_sanchez_fractional_scan.py` — 16-combination α^{P,G,We,Wa} scan.
- `diag_sanchez_modeB_and_gelbard.py` — Model A/B + Gelbard variants.
- `diag_sanchez_legendre_variants.py` — plain Legendre, √(2n+1)P̃_n.
- `diag_sanchez_N_convergence.py` — plateau at ~1.43 % proof (N=1..4).
- `diag_hybrid_f4_plus_sanchez.py` — F.4 mode-0 + Sanchez n≥1 hybrid.
- `diag_G_prefactor_scan.py` — P/G prefactor scaling scan.

## Suggested next steps

1. **Re-examine the memo rigorously.** The memo's §4 claim "`G_bc = 4·P_esc`
   matches the -4π V^{-1} / (π A V^{-1}) ratio" should imply the Sanchez
   closure reduces to F.4's N=1 — but empirically it doesn't (F.4 gives
   0.077 %, Sanchez gives 2.55 % at N=1). Either the memo's algebra has
   a factor error, or the F.4 scalar itself embeds a non-Sanchez
   approximation that accidentally gives very low error.

2. **Check F.4 as a rank-1 Mark closure, NOT Sanchez-McCormick.** F.4's
   W uses Lambert emission (cos·sin measure = µdµ once), not Sanchez's
   µ²·dµ double-µ convention. F.4 closure may be Marshak DP_0 (single
   mode) rather than Sanchez §III.F.1's IC formulation. These are
   different limits at higher N.

3. **Consider a fundamentally different closure paradigm.** The Stepanek
   form, or the Hébert 2020 §3 rank-N, or an explicit DP_N integration
   in the volume kernel may close where Sanchez §III.F.1 doesn't map
   cleanly onto ORPHEUS's observer-centered Nystrom form.

4. **N=1 bit-exactness is a HARD constraint** imposed by the Phase F.4
   regression. Any alternative rank-N closure must include F.4 as its
   N=1 reduction. This eliminates most candidates.

## Constraint held

Per the task spec: "If the recipe does NOT close and you've exhausted
systematic variants, revert all code changes." No code changes were
made to `orpheus/derivations/peierls_geometry.py` — the
`NotImplementedError` guard at
`_build_closure_operator_rank_n_white` remains in place. Only
diagnostic scripts and memory entries were added.
