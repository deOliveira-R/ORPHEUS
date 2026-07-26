---
name: dsa-rp-angular-frame
description: VERDICT — DSA's restriction/prolongation is the ℓ=0 GalerkinFrame (angular_frame(0)), consistency = Schur of a two-moment Galerkin triple product; mint nothing
metadata:
  type: project
---

# DSA (R,P) = the ℓ=0 GalerkinFrame — issue #2 Phase 3, 3-P0 verdict

Frame-attack on the DSA-for-SN restriction/prolongation pair (roadmap
`.claude/plans/stencil_assembly_dsa_roadmap.md` §Phase-3 3b; deliverable
`.claude/plans/dsa_rp_frame_analysis.md`). Branch-verified worktree facts (L-005).
The five verdicts the plan-of-record consumes:

**Why:** DSA needed to reuse the landed frame machinery, not mint an ad-hoc R/P — and
the moment-0 projector already has ≥3 spellings (Smell #16 fires as a would-be 4th).
**How to apply:** when the DSA build (or any SN acceleration) reaches 3b, hand it these
verdicts; they resolve the property-vs-type / Galerkin-vs-PG / anti-mint questions.

1. **(R,P) is GALERKIN, = `angular_frame(0)`** (the ℓ=0 sub-block of the SH
   `GalerkinFrame`; `directional.py:417`, `spherical_harmonic_basis.py` L=0 branch:
   `Y₀=1`, factor `2·0+1=1`, Gram `4π`). R = `.analysis` (moment-0 = scalar flux
   `Σ w_n ψ_n`); P = normalized reconstruction `R∘G⁻¹` = the section `M⁺` (so `M∘P=I`).
   Π = P∘M is **W-self-adjoint** (worked at S₂: `Π=(1/4π)[[w₁,w₂],[w₁,w₂]]`, `WΠ`
   symmetric) ⟹ orthogonal projector ⟹ Galerkin. **NOT PG** — NO solution weighting on
   the test side (contrast spatial homogenization's flux-weighted `φ·1_R` test ⟹ M*≠R).
   The self-adjointizing measure IS the plain quadrature `w` that `angular_frame`
   already binds. This is L-009 on the V₀ (trivial-SO(3)-irrep) sub-block: Σ_s zonal ⟹
   Funk–Hecke ⟹ scattering's eigenbasis ⟹ Galerkin, inherited by ℓ=0.
2. **Consistency = `A_low = Schur_{ℓ=1}(R₁ A_high P₁)`** — a two-moment (ℓ≤1, =
   `angular_frame(1)`) Galerkin triple product with the ℓ=1 (current J) block
   Schur-eliminated; that elimination **IS the Fick/P1 closure** (`A_φJ A_JJ⁻¹ A_Jφ`).
   "Consistent" (Alcouffe) = the triple product on the **ASSEMBLED DD operator**
   (reduce-the-discrete ≠ discretize-the-reduced — they don't commute; Phase-2 assembly
   is what enables it). 3a's matrix compare `Schur(R₁ A_high P₁) ≟ LeakageOperator` is a
   COMPUTED consistency proof.
3. **Boundary arm = SAME Galerkin frame, DIFFERENT measure** `|Ω·n|w` (half-range trace,
   measure-is-metric). Analysis face = partial currents `(J⁺,J⁻)` = two-sign ℓ=0
   half-range moments (carrier `AngularBoundaryFlux.net_current` already minted, Phase
   1). **Marshak = the boundary Fick** (an incoming-moment Schur elimination via the
   albedo `J⁻=𝒜J⁺`), NOT a reconstruction. Vacuum→Marshak(𝒜=0) is exactly the DSA
   error-problem BC (zero inflow).
4. **Foreign frames fire, all reduce to existing machinery**: two-grid (angular
   coarsening, R/P = grid transfer, A_low = Galerkin coarse operator); deflation (the
   c→1 slow mode is the ℓ=0 near-null space — explains why ℓ=0 suffices + when it fails,
   Morel anisotropy → bump to `angular_frame(1)`); preconditioned Richardson
   (`M_DSA⁻¹ = I + A_low⁻¹Σ_s`, ρ≤0.2247c) = the `KrylovAcceleration.preconditioner`
   Callable + a `SourceIteration` correction-wrap (NO new two-grid/deflation class);
   Fourier (0.2247c = the non-deflatable intermediate-mode remainder, the gate target).
5. **ANTI-MINT: instantiate `angular_frame(0)`, mint NOTHING.** No `IsotropicBasis` (the
   constant IS `SphericalHarmonicBasis(L=0)`; minting forks R∘M + ULP-drifts — the
   `iso_source_frame_conjugation_unification` verdict). No ad-hoc R/P (one frame's two
   faces). Shared primitive = `Quadrature.angular_frame` (already "single source" for
   scattering §5.6 + in-sweep moments, `directional.py:434`). Discriminating anti-twin
   test = `array_equal` (0 ULP) `dsa_R.apply(r) ≡ angular_frame(0).analysis.apply(r)` —
   a hand-rolled einsum ULP-drifts. Scattering-P0-fastpath unification
   (`add_iso_source`/`apply_p0_in_scatter`) = a SEPARATE flagged out-of-scope campaign.

**New structural observation (reusable):** Fick and Marshak are the SAME move — a Schur
elimination of an odd (interior current) / incoming (boundary partial-current) moment.
An acceleration scheme's low-order operator is the Galerkin coarse operator `R A_high P`
of a symmetry-sub-block frame, post-composed with a Schur closure of the
retained-but-eliminated moments; "consistent" = the triple product on the discrete
(assembled) high-order operator. See [[iso_source_frame_conjugation_unification]] (the
ℓ=0-is-angular_frame(0)-not-a-new-basis prior) and lessons L-009 (symmetry ⟹ Galerkin).
