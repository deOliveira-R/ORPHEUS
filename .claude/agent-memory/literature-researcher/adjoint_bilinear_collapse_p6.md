---
name: adjoint-bilinear-collapse-p6
description: P6 (#281) adjoint/bilinear-weighted group-constant collapse — what the LOCAL corpus holds (Hebert Ch3 adjoint objects + Rayleigh stationarity; Dorning Ch8 bilinear kinetics params exactly + section 8.5.4 forward-times-adjoint homogenization from asymptotics; Roy Ch4 industrial flux-weight+equivalence-factors), the verified primary chain (Larsen 75/76, Zhang-Rizwan-uddin-Dorning 95/97/97, Ussachoff/Henry), and the five missing classics that carry the per-channel formulas.
metadata:
  type: project
---

# Adjoint-weighted (bilinear) few-group collapse — extraction record (P6 #281)

Full memo: `.claude/plans/archive/p6_literature_memo.md` (2026-07-25; all
equations spot-verified on rendered pages). Complements
[[energy-condensation-collapse-formulas]] (the flux-weighted baseline).

## The local-corpus partition (no conflicts)

- **Hébert 2009 Ch. 3** [LOCAL]: adjoint transport (3.60)/(3.61) —
  kernel transposed, fission dyad χ⊗νΣ_f → νΣ_f⊗χ; φ* is a FUNCTION of
  E (not a distribution) ⟹ multigroup adjoint = lethargy AVERAGE
  (3.118)/(3.119) while forward = group INTEGRAL; multigroup adjoint
  system (3.120)/(3.121) = TRANSPOSE of the forward multigroup operator
  with the SAME flux-weighted constants (identification with ⟨φ*⟩_g
  glossed, no consistency condition stated); Rayleigh ratio (3.63)
  written POINT-WISE (no ∫d³r, no leakage; numerator = Σ_j
  ⟨χ_jφ*⟩·⟨νΣ_{f,j}φ⟩ rank-1 factored) + p. 77 stationarity statement
  (δK_eff from δΣ without δφ, δφ*) = the local perturbation grounding.
  ALL of Hébert's condensation (3.103)-(3.112), kinetics (3.125)-(3.135)
  incl. 1/V and Λ, is FLUX-weighted — φ* never enters a collapse.
- **Dorning, Ch. 8 of Azmy-Sartori *Century in Review* (Springer 2010,
  DOI 10.1007/978-90-481-3411-3_8)** [LOCAL]: THE bilinear prescriptions.
  (8.39) F_t=(N₀†,M_tΨ); (8.44) ρ=(N₀†,(H_t−H₀)Ψ)/F_t — the δk formula
  exactly; (8.45) β̄ᵢ=(N₀†,M_tⁱΨ)/F_t (a worked adjoint-weighted
  spectrum-channel collapse); (8.46) Λ̄=(N₀†,Ψ)/F_t. §8.5.4 p. 436
  verbatim: homogenized constants "weighted by both the transport theory
  flux and the transport theory adjoint flux (or importance)" — falls
  out of Fredholm solvability in multiscale asymptotics, NOT a modeling
  choice. NO stationarity/"second-order" sentence anywhere in the
  chapter (grep-verified) — use Hébert p. 77 for that.
- **Roy, Ch. 4 same volume** [LOCAL]: industrial spatial homogenization
  = flux-volume average (4.44) × μ correction factors + GET
  discontinuity factors (4.47); adjoint appears ONLY as response weight
  R_d=⟨Φ†,Q⟩ (p. 191). Bilinear collapse absent from the production
  chain.

## Verified primary chain for the bilinear-homogenization result

Larsen 1975 JMP 16:1421 (10.1063/1.522714); Larsen 1976 NSE 60:357
(10.13182/NSE76-A26897); Chiang-Dorning 1980 (ANS proc.);
Zhang-Rizwan-uddin-Dorning: NSE 121:226 1995 (10.13182/NSE95-A28560),
TTSP 26:433 1997 (10.1080/00411459708017925), TTSP 26:765 1997
(10.1080/00411459708224422); Smith 1986 PNE 17:303
(10.1016/0149-1970(86)90035-1); Ussachoff 1955 Geneva P/656; Henry
WAPD-124 1955; Henry 1958 NSE 3:52 (10.13182/NSE58-1). All
CrossRef-resolved this session.

## B&G Ch. 6 — THE per-channel source (acquired + extracted 2026-07-26)

`Bell-Glasstone(1970)Nuclear_reactor_theory.pdf` [LOCAL]; printed ≈
PDF−18; all equations verified on rendered pp. 290/297/323-326.
**§6.4h** (printed 305-308) = the bilinear multigroup derivation
(variational, P1, plane, within-group separability ansatz):
- Vector: **(6.135)** σ_{i,g}=∫_g σψ†_{i,g}ψ_{i,g}dE, i=0,1 (per
  Legendre MOMENT; ≡ ⟨φ*Σφ⟩/⟨φ*φ⟩ via unit-overlap norm (6.126)).
- Scattering: **(6.136)** σ_{i,g′→g}=∫_g dE∫_{g′}dE′ σ_i(E′→E)
  ψ_{i,g′}(E′)ψ†_{i,g}(E) — **per-pair SINK-adjoint × SOURCE-flux
  CONFIRMED**, both moments.
- Carriers: (6.125) ∫ψ_0=1 ⟹ forward carrier = PLAIN group-integral
  flux; (6.126) ∫ψ_0ψ†_0=1 ⟹ coarse adjoint = FLUX-WEIGHTED group
  average ⟹ ⟨φ*φ⟩_G = Φ*_G·Φ_G exactly (plain-carrier rows carrying
  bilinear-valued constants; the (a)/(b) ORPHEUS conventions COINCIDE
  under separability). Hébert (3.118) plain-average = flat-flux approx
  of B&G's carrier. External source collapses ADJOINT-weighted
  (6.133)/(6.134).
- Fission: (4.38) fission-in-transfer ⟹ collapses per-pair by (6.136);
  separable kernel ⟹ νσ_{f,g′→g}=χ†_g·(νσ_f)_{g′}, χ†_g=∫_gχψ†_0dE
  (Σχ†_g≠1) — dyad factored (synthesis, flagged).
- Grounding: δk = **(6.71)** p. 279 (⟨φ*,δAφ⟩/⟨φ*,Fφ⟩ exactly);
  second-order THEOREM **(6.90)** J=J₀+(δΦ†,LδΦ) p. 293 (flux-only
  weighting errs at first order (Q†,δΦ)); eigenvalue functional (6.92).
- Consistency: §6.2c p. 272 — flux-weighted multigroup adjoint ≠
  energy-integrated continuous adjoint (the Hébert-gap closer); §6.4h
  adjoint system (6.137)/(6.138) shares (6.135)/(6.136) constants,
  "clearly adjoint"; §6.4g: variational route guarantees dual
  consistency. B&G Ch. 9 point kinetics consumes §6.4h constants.
- Practical p. 308: bilinear "superior … few groups"; many groups ⟹
  correction fades. Refs [35] Pitterle-Maynard TANS 8:205 (1965),
  Little-Hardie NSE 29:402, Buslik NSE 32:233, counterpoint
  Yasinsky-Kaplan NSE 31:354; [34] Henry NSE 27:493 (time-dep).

Remaining acquisitions (corroboration only): Stacey Ch. 13, Williams
CRC (Ronen 1986, OSTI 5707826), Henry 1975, D&H.

## Process note

Zotero MCP tools were not exposed in this session's tool surface at all
(distinct from the L-007 dead-server signature) — annotations
unchecked; re-query when the surface returns.
