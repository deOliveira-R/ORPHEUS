---
name: coupled-block-operator-numerics
description: Numerics assessment of coupled 2x2 block systems [[A_AA,A_AB],[A_BA,A_BB]] for the planned CoupledOperator re-architecture (SN psi-half + DSA + multiphysics). The method menu mapped to ORPHEUS, the psi-half convergence structure (direct-in-sweep / iterated-in-scattering), wrap-vs-extract verdict, and the "substrate already exists" ruling.
metadata:
  type: reference
---

# Coupled block operators `[[A_AA,A_AB],[A_BA,A_BB]]` — numerics assessment (#2 DSA / CoupledOperator)

DESIGN investigation (2026-07-07, branch `refactor/sn-walk-unification`, the FaceField
carve C2 landed, pre-DSA). NOT a bug hunt. The "coupled system being minted" = the SN
augmented ψ½ system viewed as System A (transport bulk⊕trace) × System B (ψ½ ray/RayOp).

## The ψ½ system has TWO different 2×2 partitions — don't conflate them
1. **Ray/bulk split, WITHIN the loss operator `(L+C)`** — block-TRIANGULAR, `A_BA=0`,
   solved DIRECTLY by forward substitution (the sweep, #284). MEASURED (dense unit-probe
   over `FullField.to_flat`, GL-S4 vacuum sphere nx5/2g, dims bulk40/trace8/seed24):
   `A_sb=A_st=0.0 EXACT` (seed rows self-contained), `A_bs=7.5` (seed→bulk M-M feed = A_AB),
   `A_ss=5.0` (ray self-block = A_BB, banded straight-char). **ρ=0, no iteration.**
2. **Scattering/fission split, the FULL operator `A=(L+C)−S−B`** — the OUTER source
   iteration is block-Gauss-Seidel on the regular splitting `A=M−N`, `M=(L+C)`, `N=S+B`.
   The bulk→ray coupling `A_BA` lives ENTIRELY in the LAGGED gain `S` (`S_sb=0.183` = the ray
   scattering source; `S_bs=0` EXACT = ψ½ has zero moment weight, can't feed the scalar flux —
   the ghost-metric physics). MEASURED `ρ(M⁻¹N)=0.371 < c=0.4` (Adams–Larsen bound, below c
   for vacuum leakage). **This is the ONE genuinely-iterated coupling.**

So: ψ½ coupled solve is **DIRECT at the ray/bulk level, ITERATIVE at the scattering level**.
Answer to "is A_AB triangular-in-sweep but A_BA at SI convergence?" = YES (option b).

## Wrap-vs-extract verdict (NUMERICAL): WRAP-first
- **Welded sweep is the EXACT direct inverse**: `(L+C).solve((L+C).apply(ψ))−ψ = 3.5e-16`
  (route-a works; a WRAP inherits it bit-for-bit).
- **Extract-to-dense is PRINCIPLED-equivalent, NOT bit-identical**: welded WDD forward-sub
  vs LAPACK LU of the SAME assembled `(L+C)` agree to `5.5e-16` (a few ULP, different
  reduction trees). Naive dense `M⁻¹` on an ARBITRARY rhs (random boundary/seed) gives O(1)
  garbage (1.3) — the sweep's `.solve` treats inflow/seed rows as GIVEN DATA, not equation
  rows; an extraction MUST preserve that row-contract.
- **Recommendation:** WRAP (name A_AB/A_BA as views over the sweep) — bit-identical, ZERO
  risk to #284 exactness or ERR-053 Krylov `restart` sizing. EXTRACT to materialized blocks
  ONLY when the DSA Schur complement / preconditioner needs them; accept ~1e-15 there.

## Folded vs lagged (Q5): folded strictly better; lagged is NOT a viable oracle
- Folded (route a): direct, ρ=0, 1 pass, 4e-16.
- Lagging route-a's DIRECT seed feed: NILPOTENT (ρ=0, 2 steps) — benign but wasteful
  (because `A_sb=0`, no back-edge). MEASURED.
- Historical #282 EDGE-EXTRAPOLATION seed (ψ½=E·ψ_bulk, a genuine cycle A_AA⁻¹·A_bs·E):
  **DIVERGES ρ≈70** for sphere (documented L14/[[curvilinear-inverse-seed-taxonomy]]). This is
  why route-a was needed. NO regime where lagged-edge-extrap beats the fold (cylinder-product
  the seed is a FREE bulk edge-ordinate τ_raw=0 → pure-diagonal fold; slab no seed).
- **Rebuilding lagged as a comparison oracle = NOT worthwhile**: the edge-extrap variant
  diverges (useless), and the **dense LU of the assembled `(L+C)` is the correct
  structurally-independent oracle** (5.5e-16; F5 precedent 4e-13).

## The general machinery is LARGELY ALREADY BUILT — do NOT mint a monolithic CoupledOperator
The dagger-biproduct block substrate exists: `FullFieldSpace` (bulk⊕trace⊕ψ½ direct sum,
per-block metric dispatch) + `BlockRole` (a leaf is CLASSIFIED by which blocks it touches;
off-diagonals are NOT materialized) + `OperatorSum` (`.H` propagates to leaves via
`apply_transpose`; `(A.H).inverse()=(A.inverse()).H` #280 swap law) + `GreenOperator` (the
`A=M−N` splitting-solve as a standalone operator). The SOLVE STRATEGY is a FAMILY keyed by
spectral character, already in code:
- triangular-direct (hyperbolic transport sweep) — `InvertibleOperator`/`ScheduledInvertibleOperator`
  forward substitution.
- splitting-iterative — `SourceIteration` (`iteration.py:611-661`, `ψ=M⁻¹(q+Σgᵢ.apply(ψ))`).
- preconditioned-Krylov — `KrylovAcceleration(preconditioner=…)` (`iteration.py:777,807-812`;
  default = sweep-as-preconditioner; SN ships IDENTITY #200; DSA #2 passes the diffusion-
  corrected precond).
- dense-direct (ELLIPTIC diffusion) — `MatrixInverseOperator(FlattenedOperator(A,tpl))`
  (`diffusion/operators.py:115`); the splitting/Neumann series DIVERGES for fine-mesh elliptic,
  so NEVER `A.inverse()` there. **Spectral character picks the strategy.**
Strongest existing coupled precedent = **diffusion IS a bulk⊕boundary block solve** with an
honest un-condensed `(J⁺,J⁻)` trace; eliminating it = **Schur complement onto the bulk**
(`operators.py:98-105`). A CoupledOperator should GENERALIZE that un-condensed spelling.
**Three block-maths stay distinct** (facefield_codim1_design.md §4): saddle-point/KKT (RQI,
diffusion mixed-form #294 — deferred `SaddlePointOperator`), dagger-biproduct (the composite —
ψ½+DSA land here), lossy Petrov-Galerkin (condensation). A unifying `BorderedOperator` was
REJECTED. ⟹ CoupledOperator = a VIEW/strategy over the existing biproduct, making A_AB/A_BA
first-class for DSA — NOT a new monolith, NOT the saddle-point family.

## ρ-honest stop is a coupled-system requirement (L11/Signature 9)
`SourceIteration` stops on the INCREMENT `‖Δψ‖` (`iteration.py:637-643`), which understates the
true error by 1/(1−ρ) as the dominance ratio ρ→1 (near-critical, DSA-accelerated). The
ρ-honest residual `r=Aψ−q` ALREADY EXISTS as `evaluate_residual` (`solver.py:231-333`, the DSA
restriction substrate) but is NOT the convergence metric. A CoupledOperator's block-solve
convergence test SHOULD measure the block residual, not the increment.

Diagnostics (self-contained, promotable to `tests/sn/operators/test_psi_half_coupling.py`):
`derivations/diagnostics/diag_coupled_01_psi_half_block_structure.py` (triangular certificate +
S-carries-lagged-coupling + ρ≤c + folded-nilpotent, 4 gates) and
`diag_coupled_02_wrap_vs_extract.py` (welded=exact-inverse + extract=principled-not-bit-identical,
2 gates). All 6 green under `-O`. Backs [[curvilinear-inverse-seed-taxonomy]], [[lessons-L14]].
