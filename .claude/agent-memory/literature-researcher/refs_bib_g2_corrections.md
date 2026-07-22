---
name: refs-bib-g2-corrections
description: "#231 Phase G2: docs/refs.bib complete (59 entries covering 59 rst + 2 docstring-only keys, 2 ruled aliases); ledger of in-text-definition errors CrossRef/local-PDF corrected (fields fixed, keys immutable)"
metadata:
  type: project
---

Phase G2 (#231, 2026-07-22, seeded @ 30a073fe + docstring follow-up)
built the complete `docs/refs.bib`: every docutils `.. [Key]` definition
in `docs/theory/**/*.rst` (59) AND in `orpheus/**/*.py` docstrings
(sweep found 11; only 2 not already theory keys), keys verbatim, every
DOI CrossRef-resolved. Final: **59 entries**, pybtex-clean. Gate:
(rst ∪ py defs) − bib == the 2 ruled aliases exactly; no orphans.

**Why:** the in-text definitions drifted from the published record — 10
of 59 carried wrong pages/volumes/authors/years. The BIB fields follow
the published record; the RST definition texts still carry the errors
until the downstream in-text swap deletes them.

**How to apply:** when any of these citations is touched again, trust
refs.bib, not the page's bracket definition. Corrections ledger:

| Key | In-text error → verified truth (source) |
|---|---|
| WuXieFischer1999 | "G. Wu, Z. Wu, B. Fischer" → **Yican Wu, Zhongsheng Xie, Ulrich Fischer**, NSE 133(3):350–357 (CrossRef; key was right, def+old-seed wrong) |
| Sanchez2002 | sole author → **Sanchez, Mao, Santandrea** (CrossRef) |
| SanchezTTSP1986 | TTSP vol 14 → **15**(3):333–343 (CrossRef) |
| Stepanek1981 | NSE 78 pp 171–179 → **53–65** (CrossRef) |
| Knyazev1993 | pp 385–389 → **368–374** (CrossRef; confirms [[knyazev-1993-cylinder-anisotropic-ic]]) |
| Garcia2020 | JCP vol 393, "in a sphere" → **405**:109139, "…in a **spherical shell**" (CrossRef) |
| Garcia2021 | JCP vol 433 → **424**:109856 (CrossRef; my own older memory also said 433 — both came from the def) |
| BurkartIshiguroSiewert1976 | pp 72–81 → **72–77** (CrossRef; memory was right) |
| Carlvik1966 (key) | def year 1966 → Geneva proc. **1965** (conf. 1964), paper P/681, Vol 2 p 225 (local PDF + memory) |
| Mitsis1963 | "G. F." → **George J.** Mitsis (local PDF title page) |
| McCormickKuscer1965 | def CONFLATED: title from the 1966 JMP 7:2036 bi-orthogonality paper + numbers from the 1965 JMP 6:1939 paper. Encoded the **1965** "Half-Space Neutron Transport with Linearly Anisotropic Scattering" (10.1063/1.1704744) — matches key year, def numbers, and the usage context (Atalay Eqs 18–21, linearly anisotropic) |
| ENDF102 | def's primary text = Kellett JEFF Report 20 (decay/FY) but key + usage ("ENDF-6 format … punched-card era") mean the **ENDF-6 Formats Manual** (ENDF-102/BNL-NCS-44945) — manual encoded primary, Kellett in note; page should split into two keys |

**Ruled consolidations (the swap MUST map these; no bib entry exists):**
`[PS1982]_ → :cite:`PomraningSiewert1982`` (same JQSRT 28:503–506 work;
PS1982 entry deleted per one-concept-one-spelling) and
`[Sood1999]_ → :cite:`SoodLA13511_1999`` (docstring-only key,
`singular_eigenfunction/__init__.py:139`, same LA-13511 work).

**MetcalfZweifel1968** (60th key, docstring-only:
`singular_eigenfunction/cylinder/one_group.py:147` +
`origins/cylinder_derivations.py:1176`): thin "NSE 33, 318" resolved to
**Part II**, "Solution of the Two-Group Neutron Transport Equation—II",
NSE **33(3):318–326** (1968), DOI **10.13182/NSE68-A19240** — the
docstring page 318 = Part II's first page; companion Part I is 307–317,
DOI 10.13182/NSE68-A19239 (in the entry note — same-duo same-volume
adjacent parts, the McCormick–Kuscer confusion pattern). Provenance:
local folder holds NEITHER the paper NOR NSE volumes (per-paper PDFs
only) → CrossRef fallback, coordinator-authorized. The WM-72 citing
context (singular-subtraction for η_ν-weighted half-range integrals,
Eqs 31/33) is Metcalf citing his own earlier technique. NB: the py
docstrings repeat rst errors (e.g. one_group.py:143 "Mitsis, G. F.") —
the bib is the truth there too.

Key-vs-year
oddities accepted as-is: CarlsonLathrop1965 (chapter published 1968),
Bickley (no year in key). New DOIs minted into the corpus:
BaileyMorelChang2010=10.13182/NSE08-66, MorelMontry1984=
10.1080/00411458408211661 (pp 615–633), Pautz2002=10.13182/NSE02-1,
Pomraning1989=10.13182/NSE89-A23622, Carlvik1968=10.13182/NSE31-295
(pp 295–303, title recovered), Bickley=10.1080/14786443508561483.
