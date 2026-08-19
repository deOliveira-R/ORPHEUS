---
name: larsen-morel-2010-review-extraction
description: Full survey extraction of Larsen & Morel 2010 "Advances in Discrete-Ordinates Methodology" (Azmy-Sartori book ch.1, pp.1-84, DOI 10.1007/978-90-481-3411-3_1) — corpus §3.2 table verdicts (NO cylinder SN / NO adjoint construction / NO k-eigenvalue posing / ZERO verification), the M4b convergence-theory equation map, β≡BMC-α crosswalk, and the book-PDF workflow (page offset +12; pdftotext mangles math).
metadata:
  type: reference
---

# Larsen & Morel 2010 review — extraction record (2026-07-21)

**Cite when**: corpus-campaign §3.2/§3.3 gap tables; any "does the modern
review cover X?" question; convergence-theory (M4b) sourcing; the
symptom→cause (M2) table; β/α convention crosswalks.

**Full deliverable** (page/eq/verbatim evidence for every claim):
`.claude/plans/archive/phase_i_survey_larsen_morel_2010.md`. Source = LOCAL
`scratch/literature/Nuclear Computational Science - A Century in Review.pdf`
(whole Azmy-Sartori book, 476 PDF pages). CrossRef-verified:
DOI 10.1007/978-90-481-3411-3_1, pp. 1-84, online 2009-12-24.

## Book-PDF workflow (reusable for OTHER chapters of this same file)

- **Book page N = PDF page N+12** (ch.1 = PDF 13-96). Odd pages carry the
  running head "1 Advances in Discrete-Ordinates Methodology <p>".
- pdftotext extracts prose + references cleanly but **mangles all math**
  (Springer custom encoding: µ→\x07 etc.) — transcribe equations from
  visual `Read` passes on PDF pages; grep the text layer for locations.

## §3.2 table verdicts (all confirmed by whole-chapter greps + reads)

- **84 pp**: main text 1-75, 126 refs 76-82 (chronological bands), bios 83-84.
- **Cylinder: NO** — sphere-only curvilinear (pp.7-9 Eqs 1.20-1.26 +
  pp.36-37 flux dip 1.92-1.93, ~5 pp). 2-D cyl = ONE sentence
  - ⭐ **MINED 2026-08-12 (Q68): those pp.36-37 are §1.5.1 "Angular
    Derivatives", and their VERDICT PROSE is the load-bearing content this
    entry was missing** — *"Very few practical improvements beyond the
    elimination of the flux dip have been made"*; angular DFEM *"has not
    happened"* (hard in multi-D); the LD-in-angle scheme was *"found to be
    LESS ACCURATE than the weighted-diamond scheme"* because it abandons the
    accurate starting-direction flux, and the **hybrid** (quadratic-continuous
    first cell + LD after) is what wins. Eq. **(1.93)** = Morel-Montry's τ
    verbatim. Ref **[26] = Dudziak-O'Dell-Alcouffe LA-7911-PR (1979)** is the
    **starting-direction-flux primary**. `[SCAN PDF p.49]`. Full treatment in
    [[lathrop-2000-angular-scheme-taxonomy]].
  - ⚠ **Lesson for this file's own use:** an entry that logs a section's page
    range and equation numbers reads as covered. It is not — the *prose
    verdict* of a review section is usually the reason to have read it.
  ("generalized to 2-D cylindrical geometry [38]", p.37).
- **Adjoint of S_N: NO** — adjoint appears only as MC variance-reduction
  input (pp.21-23) + "self-adjoint" operator property (pp.49, 71 SAPD).
  No adjoint-of-the-discretization construction (B&G §6.2b style).
- **k-eigenvalue: NEVER POSED** — Eq.(1.1) fission inline, fixed-source;
  "criticality" absent; k exists only as biblio [126] (Warsa 2004 NSE 147:26).
- **Verification: ZERO** — verification/verify/manufactured/validat* = 0
  hits; "benchmark" once (p.22, MCNP-vs-ATTILA application anecdote).
- Review is application-pulled (radiative transfer / charged particles /
  well logging) — explains the reactor-side gaps. Scope claim p.1: prior
  books/reviews [32,46,71] "none of these covers the advanced work done
  during the past 20 years".

## M4b convergence-theory map (the review's unique value)

- ρ definitions (1.33a,b — b is the computable estimator); **ρ_SI = c**
  (1.34, p.12; derived 1.165-1.176, ω(λ)=(c/λ)tan⁻¹λ, sup at λ≈0; slow
  modes "depend nearly linearly on µ" p.57). Outer-iteration ρ₀ for
  radiative transfer (1.46) can be ≈0.9999 (p.17).
- **DSA**: (1.177)-(1.180); ρ ≈ 0.23c (p.58); mode-complementarity
  narrative; Alcouffe CONSISTENCY history p.59 ([6],[11],[23],[34],[35]);
  "unconditionally effective" vs "unconditionally efficient" definitions
  (p.59); partially-consistent [64] / approximate-consistent [61,67];
  Azmy [87] multi-D heterogeneity deficiency (§1.8.4) → Krylov rescue
  [124]; quasidiffusion [8] ≠ acceleration (extra truncation error).
- **Thick diffusion limit** (§1.4.4, 1.74-1.91): ε-scaling (1.78);
  hierarchy (1.82); quadrature admissibility Σw=2 (1.84), Σµw=0 (1.87),
  Σµ²w=2/3 (1.90); solvability condition = diffusion eq (1.91); scheme
  verdicts p.34 (SC→0; DD fails w/ anisotropic incident BC; LD ok 1-D +
  simplex, FAILS quads/hexes; bilinear/trilinear ok — matches
  [[multi-d-ld-closure]]); unresolved-boundary-layer open problem p.35.
- **Krylov** (§1.9): eigenvalue clustering > κ for GMRES (p.67);
  sweep-preconditioning works because integral operator = compact
  perturbation of identity (p.70); (1.218) pos-def "observed (but not
  rigorously proven [124])"; Peierls (1.219) SAPD on orthogonal meshes,
  LOST on unstructured tets [124] (p.71); DSA-preconditioned (1.220).
- **NO angular diffusion-limit asymptotics** — BMC 2010 NSE 165 is ABSENT
  from the 126-ref bibliography (chapter finalized ~2009). The review has
  ONLY the spatial half ([49],[52],[66],[105]); angular half lives solely
  in BMC ([[space-angle-discretization-separability]]).

## Convention crosswalk (extends the α-table in [[curvilinear-sweep-directness-ruling]])

- L&M-2010 sphere angular coefficients are named **β**:
  β_{n+1/2} = −2Σ_{n'≤n}µ_{n'}w_{n'} (1.23b) ⇔ β_{n+1/2}=β_{n−1/2}−2µ_n w_n,
  β_{1/2}=0 — **IDENTICAL to BMC sphere α (Eq.11)** incl. ×2 and Σw=2;
  minus sign absorbed into β (β≥0), redistribution enters (1.23a) with +,
  divisor r/w_n. M-M weighted diamond = their (1.93) ≡ BMC 15+42 ≡ Bailey 74.
- ⚠ **α is REPURPOSED as the SPATIAL WD weight** (1.30):
  ψ_{n,j}=((1+α)/2)ψ_{j+1/2}+((1−α)/2)ψ_{j−1/2}; τ=(1+α)/2. Never map
  "L&M α" to the angular recursion.
- Galerkin operator algebra (1.102-1.108): S=MΣD, D_{m,n}=P_m(µ_n)w_n
  (discrete→moment = ORPHEUS analysis face), M_{n,m}=((2m+1)/2)P_m(µ_n)
  (moment→discrete = reconstruction); Gauss ⇒ M=D⁻¹ (frame invertibility);
  S_N scattering = similarity transform of Legendre Σ (S_N/P_{N−1} slab
  equivalence). Straight-ahead invariance (1.124) → extended transport
  correction α=Σ_N leaves solution invariant, cuts c.

## Reference-list hygiene (silently correct when citing)

[87] Azmy 2002 JCP 182:213 misfiled under the 1995-1999 band; [88]
Adams-Wareing-Walters and [51] Baker-Koch printed "(1988)" but both are
1998 (NSE 130:18, NSE 128:312).
