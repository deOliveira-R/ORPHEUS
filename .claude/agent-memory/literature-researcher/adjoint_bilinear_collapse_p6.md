---
name: adjoint-bilinear-collapse-p6
description: P6 (#281) adjoint/bilinear-weighted group-constant collapse — what the LOCAL corpus holds (Hebert Ch3 adjoint objects + Rayleigh stationarity; Dorning Ch8 bilinear kinetics params exactly + section 8.5.4 forward-times-adjoint homogenization from asymptotics; Roy Ch4 industrial flux-weight+equivalence-factors), the verified primary chain (Larsen 75/76, Zhang-Rizwan-uddin-Dorning 95/97/97, Ussachoff/Henry), and the five missing classics that carry the per-channel formulas.
metadata:
  type: project
---

# Adjoint-weighted (bilinear) few-group collapse — extraction record (P6 #281)

Full memo: `.claude/plans/p6_literature_memo.md` (2026-07-25; all
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

## The gap (why acquisition is needed)

NO local source writes the per-channel fine→coarse bilinear collapse
(vector ⟨φ*Σφ⟩/⟨φ*φ⟩; scattering per-pair sink-φ*×source-φ; explicit
χ_G) — the local corpus gives the PAIRING STRUCTURE (Hébert 3.61/3.121,
Dorning 8.30/8.37/8.39) but not the collapse rule. Classics that carry
it: B&G 1970 Ch. 6 (top), Stacey Ch. 13, Williams CRC-Handbook chapter
(Ronen ed. 1986, handbook verified OSTI 5707826), Duderstadt-Hamilton,
Henry 1975. NOT in `scratch/literature/`; extraction stopped per the
no-unilateral-pivot rule.

## Process note

Zotero MCP tools were not exposed in this session's tool surface at all
(distinct from the L-007 dead-server signature) — annotations
unchecked; re-query when the surface returns.
