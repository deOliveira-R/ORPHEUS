---
name: verification-error-ontology-sources
description: #405 P2 error-ontology pull (2026-10-02) — which V&V / numerical-error primaries are LOCAL, their page offsets, key equation locators, and which are paywalled
metadata:
  type: reference
---

Memo: `scratch/reference_architecture/p2/error_ontology_literature.md` (8 sections + NEEDS).

LOCAL now (all untracked under scratch/literature, NO OCR sidecars yet — Mistral 429 that day):
- Oberkampf-Trucano 2002 SAND2002-0529 (OSTI 793406; printed p = PDF p). p.19 five error sources; p.29 Eq.(1)-(2) split at the h→0 CONSISTENCY LIMIT (finite h + I + round-off lumped); p.59 E1-E4; p.83 "numerical error is bias".
- Oberkampf-Trucano 2007 (printed = PDF+60): p.82-83 never tie model scales to mesh scales; p.83/86 series-truncation of Type-1 benchmarks.
- Arioli-Liesen-Miedlar-Strakos 2013 GAMM 36:102-129 author preprint (preprint pages): Eq.12-13 total = discretisation + algebraic (Pythagoras needs Galerkin+SPD); Eq.31-35 DWR split; Eq.37 η_NC+η_O+η_AE; Eq.40 eigenvalue needs dual residual.
- Atkinson 1997 chapters (printed pages, rendered OK): Thm 3.1.1 (3.1.35) projection error ≍ ||x−P_n x||; Nyström interp (4.1.6), Thm 4.1.2 (4.1.33/35); Thm A.1 p.516 (NO printed Neumann-tail bound); (6.5.190) p.303 names.
- Trefethen 2008 SIREV 50:67-87 (printed = PDF+66): Thm 4.2 (4.7) |a_n|≤2Mρ^-n p.75; Thm 4.3 (4.11) p.76; Thm 4.1 p.73. ATAP free = ch.1-6 ONLY (Bernstein ch.8 paywalled).
- Aurentz-Trefethen 2017 (arXiv 1512.01803): chopping = "rounding" for functions; T128 aliasing p.5.
- Lathrop 1968 NSE 32:357 ray effects (copied from NSE archive) — "not numerical but ... the discrete ordinates formulation itself".
- Already local, used: Adams-Larsen 2002 (printed = PDF+2) §I.C false convergence (1.18)-(1.19) p.10, inconsistent acceleration moves fixed point p.70; Azmy-Sartori book: Larsen-Morel p.6 "histogram σ ⇒ multigroup exact" (printed = PDF−12); Spanier ch.3 (printed = PDF−11) CLT p.134, Koksma-Hlawka p.137.

PAYWALLED (asked user): Oberkampf-Roy 2010 (DOI 10.1017/CBO9780511760396; 2nd ed 2025 10.1017/9781009031004), Roy 2005 JCP, Roache 1998, Babuška-Oden 2004, Babuška-Strouboulis 2001, Lewis-Miller 1984, Hoeffding 1963. Cambridge Core = Cloudflare wall. OSTI WAS reachable 2026-10-02 (contradicts [[verification-case-anatomy-sources]] note — retry before assuming down).
