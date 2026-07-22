---
name: adams-larsen-2002-iteration-theory
description: Full extraction of Adams & Larsen 2002 PNE 40 iteration/acceleration review — where every convergence-theory result lives, page map, notation traps
metadata:
  type: project
---

Adams & Larsen (2002), PNE 40(1):3-159 (PII S0149-1970(01)00023-3) is FULLY
extracted → deliverable `.claude/plans/phase_i_survey_adams_larsen_2002.md`
(Phase-I corpus survey, gaps M4b+M2). Local PDF + page-marked text dump
`scratch/literature/adams_larsen_2002_paged.txt` (journal p.N = PDF page N−2;
marker `=== p.N ===`).

**Headline results with locations**: ρ_SI = c continuous (2.17, p.26) AND
discrete-DD any N/h (3.28)-(3.30, p.51, Λ=(2/Σth)tan(Σthλ/2)); heterogeneous
bound min c ≤ ρ ≤ max c (3.31); DSA ω(λ) closed form (2.50, p.29),
σ_DSA ≤ 0.2247c (2.51); inconsistent DSA: ω monotone-increasing in Σth,
→ ω_SI (3.41/3.43/3.44, p.53), cell-centered variants DIVERGE [51=Reed 1971]
(p.54); consistent four-step (3.45)-(3.64) ⇒ (3.65) = Σth-independent
< 0.2247 (p.57); multi-D S2-DSA ρ=c/2 → 0.2247c as N→∞ (4.65)-(4.67, p.83);
sphere continuous Fourier = slab EXACTLY via integral-eq trick (4.35)-(4.43,
pp.78-79, "this analytic result is new"), NO discrete curvilinear Fourier
possible (p.79); cylinder = ⅓ page, zero theory (p.80). Krylov: Tφ=η̂ with
T=I−X_SI (6.4)-(6.5) IS ORPHEUS's GMRES-on-sweep-preconditioned system;
Derstine-Gelbard Table 1 (p.110) + Ashby Table 2 (p.113): GMRES rescues
inconsistent DSA but "nothing takes the place of a good preconditioner";
p.140 = the CG/GMRES-outside-DSA design-guidance quote (license for #200+#2
architecture). k-eigenvalue Ch.VIII: dominance ratio r ↔ c analogy (p.127),
SPI (8.18)-(8.25). False convergence (1.15)-(1.20, pp.9-10): error ≤ εσ/(1−σ),
fix = test ε(1−σ̂). M2 symptom→cause→fix table (14 entries) in the deliverable.

**Why:** the textbook corpus has ZERO convergence theory (corpus gap M4b);
this review is the canonical source for the docs' convergence/admissibility
section and the #2 DSA reading map (primaries Alcouffe [63], Larsen [83],
McCoy-Larsen [85], Morel [86] — all LOCAL in scratch/literature/).

**How to apply / traps:** (1) OCR renders ⅓Λ² as "½A2" in (3.41)/(3.43) —
equations in the deliverable were transcribed from PAGE IMAGES, trust those;
(2) their λ is in MEAN-FREE-PATH units (modes e^{iΣ_tλx}); (3) their cylinder
"μ" (p.80) is an ANGLE not a cosine (their ξ = ORPHEUS mu_z); (4) their
A = I−L⁻¹S and T = I−X_SI are sweep-preconditioned operators, NOT the ORPHEUS
honest matvec A = L+C−S−B; (5) printed citation slip: §VIII.B cites SPI as
"[131]" (Surzhikov) but the Wielandt ref is [132] Sutton NSE 98:169; (6) three
unrelated β's: SI-eigenvalue β(λ), TSA parameter β, CG scalar β.
