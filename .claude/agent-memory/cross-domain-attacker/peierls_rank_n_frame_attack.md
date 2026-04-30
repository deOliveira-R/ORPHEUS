---
name: Peierls rank-N white-BC closure — frame attack matches
description: Six frames attacked against the "rank-N plateaus 10–100× above F.4" conclusion. Highest-prior matches: Chandrasekhar H-function basis (two-measure harmonic analysis), Rayleigh-Ritz variational principle (F.4 IS rank-1 Ritz), and group-theoretic subgroup-trivial-irrep (white BC + SO(3) forces rank-1). Feynman-Kac surface Markov chain and differential-geometric connection-coefficient are secondary matches. Koksma-Hlawka QMC is a cross-cutting quadrature fix.
type: project
---

# Peierls rank-N white-BC — frame attack inventory

## Context

Production F.4 scalar closure reaches ~0.003% at RICH. Every
formally-consistent rank-N generalization plateaus at 0.07–1.36%.
Prior conclusion: "universal structural barrier at F.4's accuracy."
Frame attack targets: name the structural reason, or surface a
reformulation that breaks F.4.

## Matches with concrete first-test payoff

### 1. Two-measure harmonic analysis → Chandrasekhar H-function basis
- Trigger: upper-bidiagonal change-of-basis matrix M_nm between
  Lambert(dµ) and Jacobi(µdµ) ONBs.
- Insight: M_nm IS the Christoffel kernel transform for
  multiplication by µ (a three-term recurrence between two
  orthogonal polynomial families on the same interval).
  Chandrasekhar 1960 Ch. V derived the half-range polynomials
  that diagonalize the transport kernel at τ = ∞.
- First test: does α(τ) = 0.38 equal the specific
  H-function-weighted moment ratio? If yes, F.4's scalar gauge
  has a named analytic origin.
- Status: HIGH PRIOR — literature precedent exists, no ORPHEUS
  implementation.

### 2. Rayleigh-Ritz / Krein-Rutman variational principle
- Trigger: rank-2 Marshak is WORSE than rank-1 F.4 (1.36% vs
  0.003%) — a Galerkin-without-variational-principle tell.
- Insight: if F.4 IS rank-1 Rayleigh-Ritz on the Peierls-boundary
  operator, then rank-N failure is because the chosen rank-N
  bases are not nested Ritz subspaces (Marshak ⊄ Lambert).
  Monotone convergence from above is recovered with nested
  subspaces.
- First test: evaluate R[ψ_b = 1] at six grid points and check
  R[1] equals F.4 k_eff to 10⁻⁵ relative. Then build a rank-2
  nested Ritz subspace and verify monotone improvement.
- Status: HIGH PRIOR — this is the "non-linear eigenvalue"
  question 6 from the task; explicit variational framing.

### 3. Group theory — subgroup-trivial-irrep on white BC
- Trigger: white BC + isotropic scattering + SO(3) hollow sphere.
- Insight: white BC commutes only with the SO(2) × Z₂ subgroup
  fixing the surface normal. Its eigenspace on the half-range
  [0, 1] is 1-dimensional (trivial irrep). Rank-N adds modes
  that are exactly orthogonal to the BC operator's image — they
  carry no information the BC can constrain.
- First test: swap white BC for anisotropic (cosine²) BC. If
  rank-N then beats rank-1: symmetry was the barrier. If not:
  another frame needed.
- Status: HIGH PRIOR — falsifiable in one diagnostic.

### 4. Feynman-Kac surface Markov chain
- Trigger: Neumann series `(I - W)^{-1}` is a path sum; W is
  compact Markov kernel on outer surface; cavity = ballistic jump.
- Insight: stationary µ-distribution is exponential (Laplace) in
  thick limit; polynomial truncation has geometric convergence
  with τ-dependent constant — explains the 10–100× plateau.
- First test: MC-sample the stationary surface distribution,
  fit exponential tail, check ρ-independence of the decay rate.
- Status: MEDIUM PRIOR — may reproduce rather than supersede F.4.

### 5. Differential geometry connection coefficient
- Trigger: M_nm is upper-bidiagonal (three-term recurrence
  structure of a covariant derivative).
- Insight: "Lambert / Marshak mismatch" is parallel transport
  between two sections of the half-range bundle. Frame-covariant
  closure should exist.
- First test: verify µ · p_n^(Leg) = Σ c_nm p_m^(Jac) matches
  observed M_nm entry-by-entry; rewrite F.4 with conjugated W.
- Status: MEDIUM PRIOR — precedent in
  candidate_cylindrical_connection.md (not yet validated in
  ORPHEUS).

### 6. Koksma-Hlawka / randomized QMC
- Trigger: L19 quadrature-crossing pathology — signed errors
  flip sign under panel refinement.
- Insight: exp(-τ·d) is bounded Hardy-Krause variation in µ;
  QMC gives O(N⁻¹(log N)²) vs product-Gauss's τ-dependent
  constant. Scrambled Sobol' gives unbiased CIs.
- First test: 32 Owen-scrambled randomized-QMC runs at F.4
  anchor; check CI < 0.003% and sign-stability.
- Status: HIGH PRIOR on fixing the L19 protocol; orthogonal to
  the structural floor question.

## Matches checked and rejected (UNEXPLORED for this problem)

- Category theory, symplectic, de Rham/FEEC, Clebsch-Gordan,
  Bloch/crystallographic, Sobolev traces, complex analysis
  contour, dynamical systems, control theory, homogenization,
  H-matrix, graph theory.
- MaxEnt: only applies if rank-N is retained; subsumed by
  group-theory match.

## Key library insight for next session

The "rank-N non-monotone" behaviour (rank-2 worse than rank-1) is
a **diagnostic** for missing variational principle. Frame attacks
on eigenvalue-closure methods should always probe Rayleigh-Ritz
first when monotone rank improvement fails.
