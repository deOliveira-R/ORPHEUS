---
name: p6-adjoint-weighted-homogenization-rulings
description: Durable design rulings from the P6 #281 B1+B2 adjoint-weighted homogenization/condensation review (the five-morphism taxonomy, B&G energy convention, is_same_phase_space strength) — reusable for the B3 C-gates and the anisotropic-moment seam.
metadata:
  type: project
---

# P6 (#281) adjoint-weighted homogenization/condensation — design rulings

Reviewed B1 (`f6ddfb17`) + B1b (`ed7bad7c`) + B2 (`f8fe7b1b`) on
`feature/sn-adjoint-weighted-homogenization`, 2026-07-26. Verdict COMMIT-READY: 0 MUST-FIX,
2 SHOULD (T4 unwelded inline twin → [[lessons-L19]]; the `Solution.condense` representative-
spectrum twin), 2 NIT. These are the rulings I VERIFIED so a B3 (C1/C2/C3 gates) or an
anisotropic-seam review can reuse them without re-deriving.

**Why:** P6 lands the eigenvalue-consistent collapse; B3 formalizes the C-gates (spec
`.claude/plans/p6_adjoint_verification_spec.md` §4) and the moment-resolved Σ_s,ℓ refinement is
a documented seam. This note is the design baseline those reviews build on.

**How to apply:** when the next P6 slice lands, check the new code against these rulings first;
they are the "what HELD" so a divergence is the signal.

## The five-morphism taxonomy (spatial, `project_through_bilinear`)
The field owns the channel→morphism map (two morphisms forward, FIVE bilinear). Each verified
numerically correct (T1/T1b/T2/T3) and welded to the derivation builders:
- **Σ_c,Σ_L,Σ_f** (response vectors) → **pair frame**, test weight `φ*⊙φ` (T1).
- **Σ_t** (collision channel) → **collision frame**, test weight `ρ_{i,g}=Σ_n w_n ψ*_{i,g,n}ψ_{i,g,n}`
  (T1b, USER-RULED option 2 — production implements the EXACT ANGULAR pairing, not the scalar P0
  truncation). VERIFIED the ρ reduction axis: `AngularFlux.values` is `(N,ng,nx,ny)` ordinate-first
  (angular_flux.py:129 `einsum("n,ng...->g...")`), so `rho = einsum("n,n...->...", w, psi_star*psi)`
  contracts axis-0 correctly.
- **Σ_s,ℓ, Σ_2n** (matrices) → **per-pair sink×source**, `φ*_g·φ_{g'}` per (g',g) (T2). Correct
  einsum: `phi`↔from-group, `phi_star`↔to-group, channel `[from,to]`.
- **νΣ_f** → **mixed-fold** (T3): num folds `ι_i=Σ_g φ*χ` (fine), den folds `ι̃_i=Σ_g φ*χ_R`
  (collapsed). V explicit in the einsum.
- **χ** → **canonical convex average** (T3): weight `ι·p`, V rides the measure ⟹ `q_i=V_i ι_i p_i`.

`Σ_f` (fission XS, a response vector, pair-frame) vs `νΣ_f` (fission SOURCE, mixed-fold dyad) take
DIFFERENT morphisms — that split is the physics, not an inconsistency.

## The B&G Ch. 6 energy convention (`Mixture.condense` bilinear arm)
Literature-adjudicated (memo `p6_literature_memo.md` Source E). Reads like B&G with named
intermediates (NOT procedural transcription): plain condensed-flux carrier (6.125); adjoint carrier
`Ψ†_G=⟨φ*φ⟩_G/Φ_G` (6.126-8); vectors `Σ=⟨φ*Σφ⟩/⟨φ*φ⟩` (6.135); matrices sink-fold-φ*-÷-Ψ†-then-
source-average (6.136); fission FACTORED with the rank-1 simplex rescale `s` moved into νΣf so χ_C
stays a simplex (k-invariant). T6 proves exactness at true spectra + O(ε²) spectrum-error gap.

## `SNMesh.is_same_phase_space` — the pairing-guard (augmented_mesh.py)
Strength = constituent `is`-identity (mesh/quad/materials) + scheme **TYPE**. RULINGS I confirmed:
- Right HOME (SNMesh IS the phase-space carrier), right STRENGTH (the two-entry solve_sn +
  solve_sn_adjoint pattern builds two wrappers over the same constituents; raw `is` wrongly rejects).
- **Pole-closure DELIBERATELY excluded is CORRECT** — the closure is a solve-time sweep strategy,
  orthogonal to field CONTRACTIBILITY (it doesn't change ψ's layout or the quadrature ρ uses). The
  guard asks "can these fields be contracted?", not "were they solved identically." (NIT raised: the
  docstring should NOTE the exclusion so a future maintainer doesn't "fix" it by adding the closure —
  which would over-strengthen and break the two-entry acceptance.)
- Scheme is TYPE-compared (not `is`) because each wrapper builds its own scheme instance; correct
  today (schemes are parameter-free strategy types) — collapse trigger is a future LAYOUT-parametrized
  scheme (`type()` would then be too weak).

## Layering (both directions correct)
- transport `IntegratedReactionRate.evaluate(phi, *, adjoint: ScalarFlux|None)` takes a `ScalarFlux`
  (NOT `AdjointSolution`) — transport can't see sn. The bilinear folds φ* into the weight and routes
  through the ONE `InnerProductFunctional` body (no parallel reduction). Degenerate `adjoint=None` /
  flat φ*=1 is bit-identical (`1.0*x==x` + shared body).
- data `Mixture.condense(..., adjoint_spectrum: np.ndarray)` takes a bare `(ng,)` array — data can't
  see transport, so bare-numpy is FORCED here (the sn-layer `Solution.condense` extracts the typed
  fields' representative spectra and passes arrays down). Correct anti-#13 exception.

## The Moore–Penrose zero convention is consistent (verified against live code)
Every new explicit einsum uses `np.divide(num, den, out=np.zeros_like(num), where=den!=0.0)`; the
frame's `apply_inverse_metric` (numerics/space.py) is `np.where(nonzero, x/wb, 0.0)`, `nonzero=wb!=0`.
Same "zero where zero weight" — verified against `_diagonal_apply_inverse_metric`, not the docstring.

## The derivation/impl split is NOT a twin
The SymPy builders (`derivations/common/homogenization.py`) and the production numpy einsums are two
languages for one math, each independently pinned (builder by worth-zeroing theorem; production by a
hand-reference smoke). That is the algebra-of-record structural-independence, the intended design —
do NOT flag it as Pattern-2 duplication. (The genuine twin is INTERNAL — [[lessons-L19]].)
