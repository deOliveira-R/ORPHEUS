---
name: Inner-surface metric is the rank-N frontier (post-E2 summary)
description: F.4's 0.12% is quadrature-limited (true floor <0.01%); rank-(1,1,N) plateau at 0.99% is structural in basis METRIC not resolution; Jacobi-c² inner basis achieves 0.002-0.07% at σ_t·R=5 but catastrophic (247%) at σ_t·R=1 — adaptive inner-product weight is the unresolved frontier.
type: project
---

# Inner-surface metric is the rank-N frontier

## Core finding (post-Experiment-2, 2026-04-21)

Three empirical facts that reshape the rank-N problem:

1. **F.4's 0.12% residual at standard quadrature is pure quadrature error.**
   Refining (n_panels, p_order, n_ang) to (4, 8, 96) drops F.4 err to
   0.0025%. True F.4 structural floor < 0.01%. NON-monotone under
   refinement — 0.12% at baseline is a cancellation coincidence.

2. **Rank-(1,1,N) plateau at 0.99% IS structural** (not quadrature).
   Refining 2× makes it SLIGHTLY WORSE (1.08%). Inner mode energies at
   self-consistent rank-(1,1,8) solution: mode-0 = 89.8%, mode-1 = 10.2%,
   modes 2-7 = <0.04% combined. So rank-(1,1,2) captures all significant
   inner mode information but still err = 0.99%. **The 0.99% is not
   about resolution.**

3. **Jacobi-c² inner basis achieves 0.002% at σ_t·R=5 (rich quadrature)**,
   0.072% at base quadrature. BUT catastrophic at σ_t·R=1 (247% err). So
   the basis METRIC matters hugely — adaptive weight is needed but
   unresolved.

## Implication for the closure problem

The rank-N per-face white-BC closure is not about finding more angular
modes. It's about finding the RIGHT INNER-SURFACE METRIC (= weight in the
orthonormality inner product on c ∈ [0,1]). Legendre (c-weighted) and
Jacobi (c²-weighted) represent the same rank-2 information space BUT
give k_eff residuals that differ by 100-1000× depending on optical
thickness.

This explains why F.4's mode-0 closure at N=1 works (0.003% floor):
- it uses a Lambert-convention (no µ_exit weight) P/G with Marshak W.
- This mismatch is equivalent to a specific metric choice at N=1.
- The N=1 coefficient algebra absorbs the mismatch structure.

It also explains why basis rotations (split-basis, per-face, V-S) don't
break 1%: they're all Legendre-inner with different OUTER bases. Outer
basis doesn't matter (rank-1 BC bottleneck, L1). Inner metric is the
knob that moves.

## Rejected directions (post-E2)

- RH4: "Lambert P/G generalizes to split basis" — refuted (33-737% err).
- RH5: "Adding inner modes breaks plateau" — refuted (modes 2-7 are
  <0.04% of inner flux energy).

## Open frontier

- Derive the asymptotic inner ψ^+(c) at large σ_t·R. Likely form is
  some P_2(c) proxy — suggests α=0 Jacobi (i.e., c⁰ weight?) or
  P̃(c) = c√3 specifically.
- Find α(σ_t·R, ρ) that minimizes err universally. Cross-σ_t scan shows
  Legendre wins at σ_t·R<5, Jacobi-c² wins at σ_t·R≥5.
- Asymptotic matching between Legendre (thin) and Jacobi (thick) limits
  is a PROMISING next-session direction.

## Files

- `derivations/diagnostics/diag_cin_f4_quadrature_floor.py`: E2.1.
- `derivations/diagnostics/diag_cin_split_inner_enrichment.py`: E2.2 + E2.3.
- `derivations/diagnostics/diag_cin_split_lambert_pg.py`: E2.4 + E2.5.
- `derivations/diagnostics/diag_cin_split_jacobi_inner.py`: E2.6.
- Research log: `.claude/plans/rank-n-closure-research-log.md` (updated
  with Experiment 2, L5-L7, RH4-RH5, Directions I and J).
