---
name: math-foundations-citations-riesz-penrose
description: Verified units for Kreyszig 1978 (Riesz 3.8-1 p.188, adjoint Def 3.9-1 p.196), Lee 2018 IRM musical isomorphisms (Ch.2 p.26), Penrose 1955 Thm 1 p.406, Moore 1920 Bull AMS 26:394-395 — and where the scans are online
metadata:
  type: reference
---

Pure-math citations for `docs/theory/foundations/spaces.rst` (the Riesz legs,
2026-09-22 W7 pass; deliverable `scratch/_definitions/literature_check.md`).
None is in `scratch/literature/` or Zotero. All four were CONFIRMED against the
sources themselves:

- **Kreyszig 1978**, *Introductory Functional Analysis with Applications*, Wiley
  NY, ISBN 0-471-50731-8 (no DOI). §3.8 "Representation of Functionals on
  Hilbert Spaces": **3.8-1 Riesz's Theorem, printed p.188** (scan PDF p.203).
  §3.9 "Hilbert-Adjoint Operator" opens p.195; **Def. 3.9-1, p.196** (PDF
  p.211) `<Tx,y>=<x,T*y>`; existence and uniqueness = Thm 3.9-2, same page. His
  inner product is linear in the FIRST slot, so f(x)=<x,z>. Scan: physics.bme.hu
  course server (703 pp).
- **Lee 2018** IRM 2nd ed., GTM 176, DOI 10.1007/978-3-319-91755-9; Ch.2
  "Riemannian Metrics" pp.9-54 (`_2`). **Flat/sharp named on p.26**,
  subsection "Raising and Lowering Indices" (inside "Basic Constructions…",
  p.25); Eq.(2.13); the indexes give "sharp 26-28" and "musical isomorphisms 26".
  Scan offset: printed = PDF − 14 (lps.elte.hu, a corrected printing with the
  same pagination). Author's front matter and errata: sites.math.washington.edu/~lee/Books/RM/.
  Lee's sharp needs g NONSINGULAR, so a G⁺ sharp is an extension of his.
- **Penrose 1955** DOI 10.1017/S0305004100030401, 51(3):406-413; **Theorem 1
  p.406 = Eqs (3)-(6)**; does NOT cite Moore. The journal was *Proc. Camb. Phil.
  Soc.* until 1975 (the publisher's cover pages say "formerly Proceedings"). The
  PDF is fetchable from Cambridge Core's `aop-cambridge-core/content/view/...`
  URL.
- **Moore 1920**: abstract no.18 in Dresden, "The fourteenth western meeting of
  the AMS", Bull AMS 26(9) (June 1920), DOI 10.1090/s0002-9904-1920-03322-7;
  the abstract is on **pp.394-395**, and the title "On the reciprocal of the
  general algebraic matrix" is on p.389. Moore's conditions are RANGE conditions
  plus κλκ=κ, not Penrose's four equations. ams.org blocks curl behind
  Cloudflare; the open scan is archive.org `sim_american-mathematical-society-bulletin_1920-06_26`.

**How to apply:** Look here before re-verifying any of these. For another
pre-1990 Bull AMS item, use archive.org `sim_*` serial scans, which are open and
come with djvu.txt. Textbook scans turn up on university course servers, which
count as the source itself; a web page that only repeats the numbers does not.

**Frame-operator LETTERS (2026-09-22, `scratch/_definitions/frame_convention_check.md`).**
- Christensen 2016 (DOI 10.1007/978-3-319-25613-9): **T = SYNTHESIS**
  (`T{c_k}=Σc_k f_k`, ℓ²→H), read on p.558 §22.1 and p.166 §7.1. By
  derivation, T* is analysis and S = TT*. His definitions (§3.2 Thm 3.2.3,
  §5.1) are UNREAD, and the book is paywalled. ⭐ The free fallback is
  `page-one.springer.com/pdf/preview/<chapterDOI>`, which returns the first 2
  printed pages of ANY Springer chapter as a PDF and gets past the JS challenge
  that blocks `link.springer.com`.
- Casazza 2000 "Art of frame theory" (arXiv math/9910168, Type-3 fonts, so
  RENDER the pages): T = "preframe operator" `Te_n=f_n` (synthesis), T* =
  "frame transform", **S = TT*** (preprint pp.13-14). The reviewer's claim that
  it is the T=analysis source is FALSE.
- Casazza-Lynch 2016 PSAPM 73:1-51 (DOI 10.1090/psapm/073/00627; arXiv
  1509.07347, §4.1, pp.20-21): **T = ANALYSIS, T* = synthesis, S = T*T**. This
  is the ORPHEUS convention's real source.
