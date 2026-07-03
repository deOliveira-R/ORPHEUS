---
name: SH angular flux→moment projection M is rooted in anisotropic scattering (claim STANDS)
description: Falsification hunt for "the SH angular flux-moment projection M (psi→phi_l) exists BECAUSE of anisotropic scattering." Verdict CLAIM STANDS. In every documented transport method (SN, MoC, CP, PN, FCS/ray-effect, random-ray), the FLUX→SH-moment projection M=∫Y_l psi dΩ is invoked solely to evaluate the anisotropic scattering source Q=Σ(2l+1)/4π σ_{s,l} φ_l P_l. Fission is isotropic (φ_0 only); external/boundary sources are SPECIFIED in SH (input, not a flux projection); BCs are handled in ordinate space (albedo/white/specular) — only PN expresses BCs in moments because moments are PN's native unknowns. Architectural implication: HarmonicFrame belongs on the Scattering operator, NOT the angular phase-space.
type: reference
---

# SH angular flux→moment projection M: root cause = anisotropic scattering

**Task (2026-06-25):** falsify "every transport method needing the SH
angular flux→moment projection M (ψ(Ω) → φ_ℓ^m = ∫ Y_ℓ^m ψ dΩ, and
reconstruction R: moments → angular source) needs it BECAUSE of
(anisotropic) scattering." Hunt a documented, non-marginal,
structurally-independent NON-scattering use of the FLUX→SH-moment
projection.

**VERDICT: CLAIM STANDS.** No falsifier found. Every documented use of
the FLUX→SH-moment projection M across SN, MoC, CP, PN, FCS/ray-effect,
random-ray is anisotropic-scattering-rooted.

## The three load-bearing primary witnesses (all LOCAL or OA)

1. **Hébert 2009 Ch.3** (`scratch/literature/Hebert(2009)Chapter3.pdf`,
   clean text layer). The decisive equation chain:
   - **Eq. (3.55)**: `φ_ℓ^m(r,E) = ∫_{4π} dΩ R_ℓ^m(Ω) ψ(r,E,Ω)` — THE M
     operator. (Also Eq. (3.14) one-group.)
   - **Eq. (3.54)**: the SH rewrite of the **scattering source** of Eq.
     (3.51); this is the SOLE place Eq. (3.55) is invoked.
   - **Eq. (3.42)**: "The integral form [MoC/CP] is generally limited to
     **isotropic sources** in the LAB." → the angular dependence of the
     source (needing R: moments→Q(Ω)) is forced ONLY by anisotropic
     scattering.
   - **Eq. (3.57)**: fission source is **isotropic** (φ_0 / scalar flux
     only).
   - **Eq. (3.30)/(3.168)**: external/boundary source Q expanded in SH
     is *specified as input data* Q_ℓ^m — NOT a flux projection.
   - **Eqs. (3.43)–(3.47)**: albedo / specular / white / periodic BCs in
     SN/CP handled in **ordinate (direction) space**, NOT via moments.

2. **Brockmann 1981 NSE 77(4):377-414** (DOI 10.13182/NSE81-3, LOCAL,
   the canonical anisotropic-scattering review). **Eq. (47)**:
   `Φ_ℓ(x,E) = 2π ∫_{-1}^{+1} Φ(x,E,μ) P_ℓ(μ) dμ` — the M operator; the
   text: this projection "is of main interest considering the problem
   of anisotropic neutron scattering." Line ~916: the SAME Legendre
   flux-moment technique "is a natural choice in the spherical harmonics
   method, [but] is also used when the discrete ordinates method, the
   finite elements method, or the orders-of-scattering method is
   applied" — one scattering-driven projection reused across methods.

3. **Ahrens 2014 arXiv:1405.3968** ("Lagrange Discrete Ordinates", OA).
   The negative-space proof. **Eq. (7)**: `φ_n^m(r) = ∫_{S²} Ȳ_n^m(Ω)
   ψ(r,Ω) dΩ` (= M). Abstract, verbatim: "**To calculate the scattering
   source** in the LDO equations, **no spherical harmonic moments are
   needed** — only values of the angular flux." LDO removes M *precisely
   by* handling the scattering source differently; explicitly ties the
   moment-savings to "problems with strong anisotropies." An authority
   stating the only reason standard SN computes M is the scattering
   source.

## PN basis-choice question (was it chosen FOR scattering?)

YES — the SH basis DIAGONALIZES scattering; streaming merely tolerates
it. Evidence:
- **Fletcher 1983 NSE 84:33-46** (DOI 10.13182/NSE83-A17455, LOCAL).
  Eq. (7) moment equation: `(2ℓ+1) a_ℓ γ_ℓm = S_ℓm` with
  `a_ℓ = σ_t − σ_{s,ℓ}` — "the right side results because of the
  **orthogonality of spherical harmonics**" → scattering is DIAGONAL in
  ℓ. The streaming term Ω·∇ (Eq. 5 recurrences) produces the
  ℓ↔ℓ±1 TRIDIAGONAL coupling (block-tridiagonal coefficient matrix,
  Fletcher lines ~330/382).
- **Hébert §3.6–3.7**: streaming `Ω·∇=μ∂_x+η∂_y+ξ∂_z` (Eq. 3.139) under
  the Legendre recurrence `μP_ℓ=[(ℓ+1)P_{ℓ+1}+ℓP_{ℓ-1}]/(2ℓ+1)` couples
  φ_ℓ↔φ_{ℓ±1}; the within-group scattering `Σ_{w,ℓ}φ_ℓ^m` (Eq. 3.169) is
  diagonal in (ℓ,m). So: scattering diagonal (the WHY of the basis),
  streaming tridiagonal (tolerated).

## Candidate falsifiers — each examined, each NON-falsifying

| Candidate | Genuinely a FLUX→SH-moment projection? | Non-scattering? | Verdict |
|---|---|---|---|
| **MoC anisotropic** (NET 2022, DOI 10.1016/j.net.2022.03.041, OA; random-ray Tramm 2024 NSE P3, DOI 10.1080/00295639.2024.2394729) | YES — "angular flux, scattering source, cross section expanded in surface SH" | NO — the SH expansion IS the anisotropic scattering source | scattering-rooted |
| **First-collision source / ray-effect** (INL report `7246957.pdf`, OCR-readable) | YES — Eq.(6) `q_1st = HΨ_u = Σ_ℓ (2ℓ+1)/4π σ_{s,ℓ} Σ_m Φ_{u,ℓ,m} Y_ℓm`. Uncollided flux Ψ_u solved by RAY-TRACING in ordinate space (Eq.4a/5, NO moments); moments appear ONLY when feeding Ψ_u into scattering operator H | NO — q_1st IS the scattering operator H; needs φ_ℓ (ℓ≥1) iff σ_{s,ℓ}≠0 i.e. iff scattering anisotropic | scattering-rooted (ray-effect ≠ the projection's cause) |
| **Sharma 2001 ANE** "SH moments of angular flux for spherically symmetric systems" (DOI 10.1016/S0306-4549(00)00080-3; precursor S0306454999000626) | Marginal — SN solves in ordinate space, φ_ℓ EXTRACTED as OUTPUT diagnostic vs singular-eigenfunction benchmark | partly (geometry framing) | NON-falsifier: marginal (7 cites, 0 influential), output/post-processing, not a method-closing projection |
| **Functional-expansion / scattering-moment tallies** (MC, Griesheimer 2005; Legendre scattering moment tallies) | The "scattering moments" ARE σ_{s,ℓ} | NO — explicitly scattering moments | scattering-rooted by definition |
| **Anisotropic external/boundary source** (beam, surface) | NO — the SOURCE Q_ℓ^m is SPECIFIED as input (Hébert 3.30/3.168); the flux is not projected | n/a | not a flux projection |
| **Anisotropic BCs** (albedo/white/specular) in SN/CP | NO — handled in ordinate space (Hébert 3.43–3.47) | n/a | not a flux projection |
| **PN Marshak / white BC in moments** (Hébert 3.176–3.180; line 2429) | moments ARE PN's native unknowns; BC necessarily in moments | streaming/closure | NON-falsifier: PN-internal closure, not M added to another method; the moments exist for the PN solve which itself exists to diagonalize scattering |

## Architectural implication for ORPHEUS

The literature supports **"SH angular projection is intrinsically a
scattering concern"** → a HarmonicFrame belongs on the **Scattering
operator**, NOT as a general property of the angular phase-space. The
flux→moment M and source-reconstruction R are the analysis/synthesis
faces of the scattering operator's harmonic representation. The angular
*quadrature* (ordinate set) is the general angular-discretization tool;
the *SH moment projection* is specifically the scattering operator's
internal basis. Streaming, fission, external sources, and BCs all live
happily in ordinate/scalar space. The one method where moments are NOT
scattering-bound (PN) is the method that makes the moments the unknowns
*in order to diagonalize scattering* — so even there the root cause is
scattering.

## Local-folder gaps (NOT acquired, would need user)

- Askew 1972 + Halsall/CACTUS 1980 PDFs are scanned-image (no OCR text
  layer) — not extractable; Mazumdar 2015 (LOCAL, text OK) confirms MoC
  moments needed "in case of anisotropic scattering" (line 141-143).
- Sanchez 1982 "Review of neutron transport approximations" (LOCAL) is
  scanned-image — no text layer.
- Modern named anisotropic-MoC refs (Yamamoto, Sugimura, OpenMOC/MPACT
  treatment) are NOT local; OA abstracts (NET 2022, Tramm 2024)
  sufficed for the verdict.
