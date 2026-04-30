---
status: validated
project: ORPHEUS
modules: [cp, peierls, sn]
date_validated: 2026-04-30
primary_criteria: [structurally-simpler, algorithmic-advantage]
secondary_criteria: [structure-exposing, expressive]
requires_validation: false
---

# Precedent: Hilbert-Schmidt separable kernel — rank-axis is the integration variable

## Problem (current / minimal formulation)

A boundary or response operator with continuous structure in
an auxiliary variable µ (direction cosine, energy, time-of-
flight) of the form

```
K(x, y) = ∫ G(x, µ) · F(y, µ) · c(µ) dµ
```

is discretised by projecting onto a finite basis in (x, y) —
the matrix-Galerkin form

```
K_N = G_N · R · (I − T · R)^{-1} · P_N
```

where R, T are matrices of basis-overlap integrals in (x, y)
and rank-N is the truncation order in the spatial / response
indices.

ORPHEUS instance: Phase 5 continuous-µ specular multi-bounce
kernel for the sphere

```
K_bc^mb(r_i, r_j) = 2 ∫_0^1 G_in(r_i, µ) · F_out(r_j, µ) · f_mb(µ) dµ
f_mb(µ) = µ / (1 − exp(−σ·2R·µ))
```

with the matrix-Galerkin form blowing up as N grows: rank-2
gives +1.61% error and the operator norm of `(I − T·R)^{-1}`
diverges with N (verified in
`diag_specular_overshoot_02_TR_spectrum.py`).

Symptoms of the wrong-frame formulation (elegance smells
firing):

- Smell 2: dense matrix `(I − T·R)^{-1}` explicitly formed and
  inverted
- Smell 13: convergence-by-checking — no spectral-gap argument
  for the matrix-Galerkin form
- Smell 15: rank-N non-monotone — adding modes worsens the
  error on the sphere
- Additionally: the structurally distinct slab/cylinder/sphere
  behaviour (slab converges, cylinder is R-conditioning-bound,
  sphere diverges) has no named theorem — Smell 6 (iterative
  without fixed-point structure analyzed) fires through the
  Neumann series interpretation.

## Structural trigger

From reference.md Part A:

- A.7 Hilbert-Schmidt / separable kernels: kernel of form
  `K(x, y, µ) = a(x, µ) · b(y, µ) · c(µ)` integrated over µ
- A.3 Spectral theory of multiplication operators: continuous-
  limit `T·R` is multiplication by `m(µ) = exp(−σ·2R·µ)`;
  `ess_range(m) ∋ 1` for sphere (at µ → 0+) — matrix-Galerkin
  diverges

From reference.md Part C:

- Smell 2 fires
- Smell 13 fires
- Smell 15 fires (the structural tell on the sphere)

## Frame

Hilbert-Schmidt separable-kernel theory. Core objects:

- The kernel slice at fixed µ has rank 1: `K(·, ·, µ) =
  a(·, µ) ⊗ b(·, µ)`
- The total kernel rank is bounded by the cardinality of the
  µ-quadrature
- The natural finite-dimensional axis of the operator is µ,
  NOT the spatial indices x, y
- The structural mistake of matrix-Galerkin is treating
  (x, y) as the rank axis when µ is the rank axis

## Reformulation (sketch)

Discretise the integration variable directly with a
quadrature `{(µ_q, w_q)}_{q=1..Q}`:

```
K[i, j] = Σ_q w_q · a(x_i, µ_q) · b(y_j, µ_q) · c(µ_q)
```

No matrix inverse. No rank-N basis truncation. Storage
O(I·Q + J·Q); cost O(I·J·Q).

For the Phase 5 multi-bounce ORPHEUS instance: M1 reads

```
K_bc^mb[i, j] = Σ_q w_q · G_in(r_i, µ_q) · F_out(r_j, µ_q) · f_mb(µ_q)
```

with Q = 64 Gauss-Legendre on µ ∈ (0, 1] for homogeneous, or
Q = 32 piecewise-GL with knots at impact parameters for
multi-region.

## Elegance payoff

- **Structurally-simpler**: no matrix inverse;
  `(I − T·R)^{-1}` is replaced by a single quadrature sum.
  The conditioning hazard (`R = (1/2)M^{-1}`) disappears.
- **Algorithmic-advantage**: O(I·J·Q) deterministic cost
  vs O(I·J + N³) for matrix-Galerkin (where N is the
  rank-truncation order that fails to converge on the
  sphere). No iterative refinement; no condition-number
  blow-up at large σR.
- **Structure-exposing**: the rank-1-in-µ structure of the
  kernel becomes the primary object. The frame names which
  axis is the rank axis, which the matrix-Galerkin form
  hides.
- **Expressive**: the same form composes across slab,
  cylinder, sphere by changing only `G_in`, `F_out`, `f_mb`
  — the geometry-specific pieces — while the discretisation
  remains uniform.

## Concrete first test

The validation was done as part of the Phase 5 frame attack
on 2026-04-28; M1 is the production reference for the sphere
multi-bounce kernel.

Verification triad — reusable for any rank-truncation method
validation in this kernel class:

- **M1 (production reference)**: separable-quadrature form
  `K[i, j] = Σ_q w_q · a(x_i, µ_q) · b(y_j, µ_q) · c(µ_q)`.
  This is what production code calls.
- **M2 (verification cross-check)**: bounce-resolved
  generating-function expansion. Use
  `1/(1 − e^(−x)) = Σ_k e^(−kx)` to write
  `K = Σ_k K^(k)` where each `K^(k) = ∫ a · b · µ · e^(−k·σ·2R·µ) dµ`
  is non-singular. Truncation at K_max gives an explicit
  residual bound. M1 = M2 to 1e-5 at K_max ≈ 10 for
  τ_R = 2.5.
- **M4 (ground truth)**: independent reference. For the
  Phase 5 instance, MC k_eff on the surface Markov chain
  `Σ_k (T·R)^k` (specular preserves µ, so the chain
  factorises by µ and the closed form `1/(1 − e^(−σ·2Rµ))`
  is recovered). For other instances, closed-form
  references exist for restricted parameter regimes (e.g.,
  Sanchez 1986 TTSP, Trefethen-Embree §15 for
  multiplication-operator examples).

Pass condition for the Phase 5 instance:

- M1 vs M4 agree to < 0.1% at τ_R = 2.5 with Q = 64 GL on µ
- M1 vs M2 agree to 1e-5 at K_max = 10
- Matrix-Galerkin at N = 8 gives +1.61% (the discriminator —
  the wrong-frame method is structurally separated from M1
  by an order of magnitude)

## Literature path

- Tricomi, "Integral Equations" (1957) Ch. 2 — degenerate-
  kernel theory, separable kernels, finite-rank
  approximation
- Trefethen & Embree, "Spectra and Pseudospectra" (2005),
  §15 — discrete approximations of unbounded multiplication
  operators do not converge in operator norm; the named
  theorem behind matrix-Galerkin failure on the sphere
- Sanchez, "Approximate solutions of the two-dimensional
  integral transport equation by collision probability
  methods" (1986) TTSP — restricted closed-form references
  used for cross-validation

## Transferable pattern

Look for this reformulation when you see:

1. An operator with continuous structure in an auxiliary
   variable µ (direction cosine, energy group, time-of-flight,
   collision number) of the form
   `K = ∫ a(·, µ) ⊗ b(·, µ) · c(µ) dµ`
2. Discretisation that projects onto a basis in the spatial
   / response indices (matrix-Galerkin) and inverts a matrix
   that depends on the µ-structure
3. Rank-N non-monotone convergence (Smell 15) when basis
   order is increased — especially worse on geometries where
   the µ-integrand has a singular endpoint
4. A multiplication-operator interpretation where
   `ess_range(m) ∋ 1` over a non-zero-measure set (A.3
   sharpening)

The reformulation moves the discretisation from the spatial
indices to µ. Expected payoff: matrix inverse disappears,
conditioning hazard disappears, geometry-specific pathology
becomes a property of `f_mb` rather than the discretisation.

Apply to (candidate targets):

- Energy-group response kernels of the form
  `R(g → g') = ∫ ψ(E, g) · χ(E, g') · σ_t(E)^{-1} dE` —
  rank axis is E, not (g, g')
- Time-of-flight kernels in pulsed-source problems
- Collision-number-resolved escape probabilities (the
  M2 form is itself this pattern with the bounce count k as
  the rank axis)
- Any rank-N-closure investigation: before adding modes,
  test whether the integration variable is the rank axis
  the existing method is missing

The verification triad (M1 production / M2 bounce-resolved /
M4 ground-truth) generalises: production = direct quadrature
on the integration variable, verification = a different
expansion of the same kernel, ground-truth = a
discretisation-independent reference (Monte Carlo, closed
form, or external code).
