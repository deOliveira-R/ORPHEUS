---
name: harmonic-frame-ownership-funk-hecke
description: Falsifiable representation-theory verdict on WHO OWNS the angular spherical-harmonic frame (M:AngularFlux→HarmonicMomentFlux, R:reconstruct) — scattering (its eigenbasis) vs the angular phase-space (a general L²(S²) tool). VERDICT the claim "M is intrinsically a scattering concern" is CONFIRMED on structural grounds (NOT non-falsifiable): M is LITERALLY the change-of-basis into the scattering operator's eigenspace by Funk–Hecke (Σ_s(Ω·Ω') zonal ⟹ {Y_ℓ^m} diagonalize it, eigenvalues = Legendre moments Σ_{s,ℓ} = the diagonal of Λ) + Schur (Σ_s· ∈ SO(3)-commutant ⟹ scalar-per-ℓ-block with (2ℓ+1) degeneracy). Streaming Ω·∇ does NOT diagonalize (it is the ℓ=1 tensor operator ⟹ ℓ↔ℓ±1 PN recurrence by Clebsch–Gordan) — the asymmetry IS the crux and it OWNS the basis to scattering. The four falsification candidates ADJUDICATED: (A) "ψ∈L²(S²)" licenses an INFINITE basis but NOT the TRUNCATED ℓ≤L M/R the code has (L = the scattering spectrum's support Σ_{s,ℓ}=0 for ℓ>L) → scattering-derived; (B) anisotropic external source SH-expansion → not live (sources are iso, ℓ=0-embedded) + structurally weak (basis still dictated by the operator it feeds); (C) angular detector-response functional of order L_d>L_scatter → THE genuine falsifier-in-principle (M-truncation set by the DETECTOR not the kernel) but ABSENT from ORPHEUS (only ℓ=0 scalar flux via integrate_angular, not the frame); (D) PN streaming closure → scattering-rooted (moments exist only because the basis diagonalizes collision). Architecture ALREADY CORRECT: HarmonicFrame in transport/frames/ (thin GalerkinFrame), ScatteringOperator HOLDS it (scattering.py:510 cached_property), Λ-on-S / M,R-on-frame spectral-theorem split (A=R·Λ·M, U=R/M generic diagonalizing unitary REUSABLE, Λ=operator-specific spectrum), L=scattering_order binding (scattering.py:545). RELOCATION TRIPWIRE = 2nd frame consumer with an L independent of scattering_order (candidate-C detector, OR a PN/SPN flux order L_flux>L_scatter) ⟹ the CONSTRUCTOR ownership scattering_op.frame relocates off S to the neutral Quadrature.angular_frame(L) factory (which already exists + anticipates it). CROSS-METHOD DISANALOGY (sharpens it): energy condensation is NOT an eigenbasis (group-transfer matrix has NO rotational symmetry, no Funk–Hecke) → its frame is solution-WEIGHTED PetrovGalerkin, NOT operator-owned Galerkin. ONE cause explains the Galerkin-vs-PG split across the whole campaign: angular=Galerkin BECAUSE SO(3)-eigenbasis-orthogonal; spatial/energy=PG BECAUSE no symmetry. Read on main, clean tree — file:line ground truth (L-005 N/A).
metadata:
  type: project
---

# Who owns the angular SH frame — the Funk–Hecke / Schur ownership verdict

DESIGN VERDICT answering a falsifiable architectural claim posed head-on:
does the angular flux→moment projection M (and its reconstruction R) belong
ON the scattering operator (scattering's eigenbasis) or ON a general angular
phase-space tool? Read on `main`, clean tree — `scattering.py` /
`frame.py` / `harmonic_frame.py` / `harmonic_moment_flux.py` /
`harmonic_moment_source_sink.py` / `directional.py` / `solver.py` /
`loss_representation.py` file:line is GROUND TRUTH (no Nexus/worktree
staleness; L-005 N/A). This is the FIRST of the six M/R memos to answer
OWNERSHIP (the prior five — `projection_reconstruction_frame_pair`,
`projection_discipline_hierarchy`, `unified_frame_api_design`,
`coefficient_field_promotion`, `homogenization_measure_derivation` — all
work WITHIN the M/R pair: naming, weight-homes, type-vs-property, API verbs,
the Galerkin derivation). This one is one level UP: the FRAME's HOME.

## The load-bearing theorem (the claim is CONFIRMED, not non-falsifiable)

**M is LITERALLY the change-of-basis into the scattering operator's
eigenspace.** This is a theorem, not an analogy:

- **Funk–Hecke.** For a zonal kernel k(Ω·Ω') on S², the spherical harmonics
  are eigenfunctions of (T_k f)(Ω)=∫k(Ω·Ω')f(Ω')dΩ', eigenvalue
  λ_ℓ=2π∫₋₁¹k(t)P_ℓ(t)dt depending ONLY on ℓ. The scattering source operator
  Σ_s· IS such a T_k; its eigenvalues are exactly the Legendre moments
  Σ_{s,ℓ} = the diagonal of `LegendreMomentScattering` Λ (`scattering.py:211`,
  domain==codomain==SphericalHarmonicSpace `:349-368`). So M = the unitary U
  diagonalizing Σ_s·; Λ = mult by the eigenvalue spectrum; R = synthesis;
  S_aniso = R∘Λ∘M = the operator in its OWN eigenbasis (spectral theorem
  A=UΣU*). The addition theorem `Σ_m Y_ℓ^m(Ω)Y_ℓ^m(Ω')=(2ℓ+1)/4π·P_ℓ`
  (`scattering.py:79`) IS the spectral resolution of the zonal kernel.
- **Schur's lemma.** L²(S²)=⊕_ℓ V_ℓ is the isotypic decomposition into SO(3)
  irreps (dim V_ℓ=2ℓ+1). Any SO(3)-invariant operator (commutant) acts as a
  SCALAR per V_ℓ — forcing the block-diagonal-per-ℓ + (2ℓ+1) degeneracy. The
  (2ℓ+1) factor and g_C=4π/(2ℓ+1) are the irrep dims / SO(3) Plancherel
  weights; the addition theorem is Peter–Weyl.

**The asymmetry IS the crux.** Streaming Ω·∇ does NOT commute with SO(3) (Ω
is the ℓ=1 vector irrep); by Clebsch–Gordan V_1⊗V_ℓ=V_{ℓ-1}⊕V_ℓ⊕V_{ℓ+1} it
connects V_ℓ↔V_{ℓ±1} — exactly the PN moment recurrence (a tridiagonal, NOT
block-diagonal). The basis is chosen to diagonalize COLLISION; streaming is
then expressed (awkwardly) in it. ⟹ the basis is OWNED by scattering.

## The four falsification candidates — ADJUDICATED

| Candidate | Verdict | Structural reason |
|---|---|---|
| (A) ψ∈L²(S²) ⟹ any angular fn SH-expandable | SCATTERING-DERIVED | The infinite expansion is basis-agnostic but USELESS; the code's TRUNCATED ℓ≤L M/R (`HarmonicMomentFlux` (L+1,2L+1,…)) has L = the scattering spectrum's support (Σ_{s,ℓ}=0 for ℓ>L). "ψ∈L²" licenses an INFINITE basis, not the FINITE M/R. R∘M=band-limited projector≠I (`frame.py:74-77`), band CHOSEN to match the kernel. |
| (B) anisotropic external source in SH | NOT LIVE + weak | Every external source is ISO, ℓ=0-embedded via `from_isotropic`/`_lift_external_source_to_moments` (`loss_representation.py:1966,:456`) — a trivial injection, NOT an analysis M. Every `HarmonicMomentSourceSink` is BORN from Λ (`harmonic_moment_source_sink.py:16-23`). Even an aniso source's reason to be in {Y_ℓ^m} is to ADD to the scattering source in the SAME spectral coords — basis still operator-dictated. |
| (C) angular detector-response functional, order L_d | **THE GENUINE FALSIFIER-IN-PRINCIPLE** | ⟨R_d(Ω),ψ⟩ with R_d aniso to order L_d>L_scatter ⟹ an M-projection whose truncation L_d is set by the DETECTOR, not Σ_s — structurally INDEPENDENT. But ABSENT from ORPHEUS: the only output functional is the ℓ=0 scalar flux via `integrate_angular` (`solver.py:1349`), which does NOT use the frame. STRUCTURALLY genuine, EMPIRICALLY latent. |
| (D) PN streaming closure | SCATTERING-ROOTED | The closure truncates the ℓ↔ℓ±1 STREAMING recurrence, but the moments exist only because the basis diagonalizes COLLISION. The order L is scattering_order (Pℓ-scatter ⟹ PL flux). "Streaming in scattering's eigenbasis, then truncated." |

**Verdict on the claim.** CONFIRMED on structural grounds, with ONE named
escape hatch. The user's exact phrasing — "every method that needs
HarmonicFrame needs it because it's giving special attention to scattering" —
is structurally almost-non-falsifiable: non-falsifiable for the FLUX
projection (the flux's only reason to be in SH is that Σ_s, the operator
coupling it to itself, is diagonal there — A+D), falsifiable ONLY by an
OUTPUT functional (C) the codebase does not have. There is no
structurally-independent FLUX-side consumer; the structurally-independent
consumer is OUTPUT-side and latent.

## Architecture: ALREADY CORRECT (do not move the frame either way)

The placement is already in the code and right: `HarmonicFrame` in
`transport/frames/` as a thin `GalerkinFrame` specialization;
`ScatteringOperator` HOLDS one (`scattering.py:510` `cached_property frame`);
the windowed-SI resolvent sources its frame FROM `scattering_op.frame`
(`solver.py:617`, comment: "M sourced from the scattering operator's own
quadrature ⇒ bit-identical to S's internal projection"); the fuller-view
oracle (`operator.py:980`) IS the scattering apply's view.

- **Do NOT bury the frame as a private field of S.** The GalerkinFrame
  machinery (analysis/reconstruction/conjugate, `frame.py`) is method-agnostic
  and SHARED with the indicator-homogenization PetrovGalerkinFrame — the SAME
  spectral-decomposition abstraction (project→act-in-eigencoords→reconstruct)
  serves scattering's R∘Λ∘M AND homogenization's G⁻¹∘M. Burying it re-forks it.
- **The correct split is in place:** Λ (role-changing eigenvalue-mult edge) on
  the scattering operator (`LegendreMomentScattering`); M,R (role-preserving
  change-of-basis) on the frame (`harmonic_frame.py:48-55` documents this
  explicitly). Spectral-theorem factoring A=R·Λ·M: U=R/M is the diagonalizing
  unitary (generic, REUSABLE); Λ=the spectrum (operator-specific). Λ-on-frame
  or M,R-on-S would each be the error; the code makes NEITHER.
- **"Scattering owns it" = is the canonical CONSTRUCTOR + the L-binding,** which
  is exactly `angular_frame(self.scattering_order)` (`scattering.py:545`). That
  binding "frame's L = scattering's L" IS the architectural ownership statement.

## THE RELOCATION TRIPWIRE (the flip condition — record it)

The moment a SECOND frame consumer arrives with a truncation L INDEPENDENT of
scattering_order:
- candidate-C: an angular detector response R_d(Ω) of order L_d, OR
- a PN/SPN solver whose flux expansion order L_flux > L_scatter

…M/R stop being scattering-owned (≥2 structurally-independent consumers at
L_scatter and L_d). The CONSTRUCTOR ownership `scattering_op.frame` relocates
off S to the neutral `Quadrature.angular_frame(L)` factory — which ALREADY
EXISTS (`directional.py:417`) and ALREADY anticipates this (the
`angular_frame` naming tripwire `:439-447` names the AXIS not the basis,
exactly so a 2nd consumer is a signature change not a rename). The frame
object itself stays in `transport/frames/`; only the CONSTRUCTOR home moves
(both S and the detector then bind their own L off the neutral factory).
First test that REDs the current binding: a PN closure with L_flux>L_scatter
needs `max(L_flux,L_scatter)`, NOT `scattering_order` (`scattering.py:545`).

## CROSS-METHOD DISANALOGY (sharpens the verdict — energy is NOT an eigenbasis)

The SAME eigenbasis logic does NOT transfer to energy condensation, and the
disanalogy EXPLAINS the campaign's Galerkin-vs-PG split:
- Energy group-transfer Σ_s[g'→g] / χνΣf is a GENERAL G×G matrix — NO
  rotational structure, NO Funk–Hecke, NO clean spectrum. So the energy
  condensation basis (coarse-group indicators {1_G}) is NOT the eigenbasis of
  the transfer operator. Its frame is a flux-WEIGHTED projection (the
  `homogenization_measure_derivation` Galerkin-in-L²(φV) / `unified_frame_api`
  PetrovGalerkin), a MEASURE/test-side choice, NOT operator-owned.
- ⟹ ONE structural cause for the whole campaign's discipline split: **angular
  SH = Galerkin BECAUSE it is the SO(3)-eigenbasis (orthogonal, M*=R up to
  g_C); spatial/energy = PetrovGalerkin BECAUSE no symmetry (solution-weighted
  test≠trial).** Do NOT over-generalize "operator owns its eigenbasis frame"
  to energy — energy has no eigenbasis, only a weighted projection. This
  unifies and EXPLAINS the prior memos' Galerkin/PG verdicts from the
  representation-theory root (they asserted the split; this derives WHY).
- First test (discriminates the two disciplines from one cause): assert the
  angular HarmonicFrame is `GalerkinFrame` (test is trial, Gram = SO(3)
  Plancherel g_C=4π/(2ℓ+1) diagonal by Schur) while energy condensation is a
  genuine `PetrovGalerkinFrame` (test=φ·1_G≠trial=1_G, Gram=1/Φ_G diagonal by
  disjoint support but PG ∵ solution-weighted). M*=R for the former, M*≠R for
  the latter.

## Refuted frames (durable UNEXPLORED for this problem class)

- **Differential geometry / connection coefficients** — no curvature/1/r term
  in the algebraic SH contraction on fixed S²; the curvilinear angular
  redistribution (1−µ²)/r∂_µ is a SEPARATE (streaming-side) operator, not the
  scattering eigenbasis question. (L-001/L-008.)
- **Category theory / adjunction (M⊣R)** — duality fully captured by the frame
  analysis/synthesis pair + the named SO(3) Plancherel metric g_C; Funk–Hecke
  + Schur ARE the concrete content. "Adjunction" adds vocabulary, not a lever.
  (L-001.)
- **Tensor networks / MPO** — Y dense rank-full per ℓ; the per-ℓ block
  structure is a direct sum (isotypic), not a tensor-train. No bond dim.
  (L-001.)
- **Homology** — no ∂²=0; R∘M is the band-limited projector (idempotent), not
  a boundary differential. (L-001.)
- **Wiener–Hopf / Chandrasekhar H-function** — wrong solver family (half-space
  spectral), orthogonal to the eigenbasis-ownership question. (L-001.)
- **Optimization / variational** — the SH expansion minimizes the L²(S²)
  truncation residual (Rayleigh–Ritz), but adds no lever the spectral frame
  lacks; the eigenbasis is FORCED by kernel symmetry, not chosen by a
  minimization. (L-001.)

## Elegance assessment

- **Structure-exposing (strong):** names the moment field as the SPECTRAL
  COORDINATE of Σ_s (its eigenbasis representation), not the
  `harmonic_moment_flux.py:6` "natural data carrier of the Galerkin pipeline"
  ("pipeline" = operational language for "kernel in its own eigenbasis", the
  Smell-#15-adjacent unnamed-native-frame tell). The Galerkin-ness is DERIVED
  (self-adjoint zonal kernel ⟹ orthogonal eigenbasis ⟹ Galerkin), not an
  independent discipline choice (`frame.py:42`); the PN recurrence is DERIVED
  (Clebsch–Gordan, Ω=ℓ=1), not asserted.
- **Structurally-simpler:** ONE principle ("which operators are in the SO(3)
  commutant") classifies every transport operator's SH-basis structure
  (collision/iso-scatter=scalar-per-block; aniso Σ_{s,ℓ}=ℓ-dependent scalar;
  streaming=ℓ=1 ladder) AND explains the cross-campaign Galerkin-vs-PG split.
- **Confirms, does not refactor:** the rare attack whose structural verdict
  VALIDATES the current placement (S holds the frame, Λ-on-S/M,R-on-frame
  split, L=scattering_order). No Smell #16 on ownership — the frame is already
  single-sourced. The deliverable is the THEOREM (M=diagonalizing unitary) +
  the named relocation tripwire, not a code change.
