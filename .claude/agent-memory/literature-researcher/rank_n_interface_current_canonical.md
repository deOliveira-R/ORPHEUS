---
name: Rank-N per-face interface-current canonical references
description: Which paper/book contains the authoritative DP_{N-1} surface-to-surface W matrix, normalization convention, and reciprocity form for Peierls solvers.
type: reference
---

**Canonical source for rank-N Marshak / DP_{N-1} per-face interface currents**:
Sanchez & McCormick (1982) "A Review of Neutron Transport Approximations,"
*NSE* 80(4), 481-535, DOI `10.13182/nse80-04-481`. §III.F is the
**only identified source** in the standard reactor-physics textbook
corpus (Ligou 1982, Sanchez 2002, Stamm'ler 1983 Ch.IV, Stacey 2007
Ch.9 all use scalar/DP-0). See
`rank_n_closure_four_references_synthesis.md` for the full
cross-check.

**Warning**: the §III.F.1 formulation is unique to the 1982 review
and has no independent cross-validation in the four textbook
references. The ORPHEUS implementation plateaus at ~1.42 % k_eff
error (N = 2-4, σ_t·R = 5) and is consistently 3-13× worse than the
scalar F.4 closure at σ_t·R ≤ 5. Treat §III.F.1 k_eff claims with
scepticism; validate any rank-N closure against the σ_t → 0
conservation identity `W_oo[m,n] + W_io[m,n] = δ_{mn}` before
trusting it.

Still untested candidates for a second rank-N derivation:
**Hébert (2009/2020)** *Applied Reactor Physics* Ch. 3 §3.4-3.5
(DOI `10.1515/9782553017445`, 3rd ed.) and Stepanek primary papers.
Stamm'ler-Abbate 1983 is NOT a source for rank-N — Ch.IV uses only
scalar CP.

**Normalization rule** (load-bearing, easy to get wrong):
- `(2n+1)` Gelbard factors belong in the **reflection operator R**,
  not in the **W transmission matrix**. So `R = diag(1,3,5,...,2N-1)`
  and `W` integrand has only `cos θ · sin θ · P_tilde_m · P_tilde_n
  · exp(-Σ_t ℓ)`.
- Physical surface response operator is `W · R`, not `W`.

**Generalized reciprocity at N ≥ 1** (Sanchez-McCormick eq. III.F.20):
`A_i · W_{ji}^{mn} = A_j · W_{ij}^{nm}` — note the **transposed mode
indices**. Scalar reciprocity (W_oi = (R/r_0)²·W_io for a sphere) is
the m=n=0 special case. Use as a built-in test: off-diagonal
reciprocity must hold *across* the transpose.

**Cross-mode chord-angle at exit face** (spherical case): if emission
angle at R is θ, the cosine at r_0 is `μ' = sqrt(1 - (R/r_0)² sin²θ)`
— EXACT Snell-like geometry, not an approximation. Sanchez-McCormick
eq. III.F.12.

**No published rank-N k_inf benchmark exists for a hollow spherical
annulus.** Closest truth sets (Sood-Forster-Parsons 2003, Ganapol
2024 *Foundations* 4(3) 27) are slab-only. For curvilinear,
validation has always been code-to-code. If ORPHEUS produces a
machine-precision convergence table at rank ≥ 8, that is itself a
publishable result.
