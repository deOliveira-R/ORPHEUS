---
name: Elegance smell — rank-N non-monotone convergence
description: When a rank-N truncation ladder shows higher ranks performing WORSE than lower ranks, the underlying method is Galerkin projection without a variational principle. Rayleigh-Ritz with nested subspaces guarantees monotone convergence from above; non-monotone behaviour is a structural tell that the Ritz frame was abandoned.
type: feedback
---

# Elegance Smell 15 — Rank-N non-monotone convergence

## The smell

A basis-truncation ladder (rank-1, rank-2, rank-3, ...) where the
error does NOT decrease monotonically. Classic example:

- ORPHEUS Peierls rank-N white BC: rank-1 F.4 gives 0.003%,
  rank-2 Marshak gives 1.36%, rank-3 Marshak gives 1.36%.
- PCA-sector rank-N: M = 2 gives 0.007%, M = 3 gives 0.008%,
  M = 5 gives 0.008%.

**Why:** Galerkin projection onto an arbitrarily chosen basis
has no monotonicity theorem. Rayleigh-Ritz with NESTED subspaces
has the Courant-Fischer min-max theorem guaranteeing monotone
convergence from above.

**How to apply:** When you see rank-N non-monotone, probe for
the variational principle the method is supposed to minimize /
maximize. If rank-N is doing Galerkin without a Ritz structure
(non-nested bases, projection in the wrong Gram matrix), the
fix is to switch to Ritz — not to add more modes.

## Matching frame

A.5 Variational calculus / optimization — Rayleigh quotient.

## Diagnostic first test

Compute the Rayleigh quotient R[ψ_b] of the converged rank-1
ansatz directly. If R[rank-1-ansatz] matches the production
eigenvalue to 10⁻⁵ relative: the method IS variational, and the
Ritz frame applies. Build rank-2 as a NESTED subspace (Gram-
Schmidt orthogonalize the rank-2 basis against rank-1) and
confirm monotone decrease.

## Literature

- Courant & Hilbert 1953, Vol. I, §VI — min-max theorem.
- Wendroff 1961 NSE — variational discrete ordinates.
- Case & Zweifel 1967 Ch. 6 — variational methods for
  Boltzmann.

---

**Promoted to skill on 2026-04-30** — see `cross-domain-frames`
reference.md Smell #15. This memory file remains as evidence /
precedent until the live rank-N work threads close.
