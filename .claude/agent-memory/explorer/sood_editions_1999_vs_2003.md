---
name: sood-editions-1999-vs-2003
description: Where Sood LA-13511 (1999) and the 2003 journal edition differ in values, which the primary source backs, and the OCR sidecar trap
metadata:
  type: project
---

Measured 2026-09-29 (P1 step 2b census; memo `scratch/reference_architecture/p1step2b/edition_census.md`).
The user ruled that day that the 2003 edition is the reference.

- **The value differences are confined to problems 3, 4, 50-52, 67, 71, 72 and 73.** The rest is renumbering: 1999 tables = 2003 + 3, and 1999 Eq 19/20/28/29/32 = 2003 A.2/A.3/A.11/A.12/A.15. The literature-researcher memory `sood_2003_vs_1999_extraction.md` says "same radii, same ratios", which is false.
- **The 2003 XS are the primary's.** Siewert & Thomas 1986 (OCR in `scratch/literature_ocr/`, Tables I-II) gives exactly the 2003 values. The 1999 U-Al set is a 5-significant-figure rounding of them. The one inversion is the U-Al flux ratio: the primary's exact 3.125 equals 1999's, while 2003's 3.124951 comes from its own rounding of Σ12s.
- **Trap:** the OCR `.md` holds only `[tbl-N.md]` links. The table cells live in the `.mocr.json` sidecar (`pages[i].tables[j].content`). A grep of the `.md` for a tabulated value returns 0. Inline the tables first.
- **Instrument that worked:** a decimal-token multiset diff between the two editions' problem sections gives the complete delta in one pass. Control it with one known token per side.

- **2003 appendix equations print as "(.N)", with no "A"** (page image p. 99). "A.N" is our normalisation. 1999 Eq N = 2003 (.N−17) for the whole appendix; only Eq 52 and Eq 76 changed bodies.
- **References do not shift uniformly**: 1999 [1]→[7]; [2]-[25] +7; [26]-[42] +9; [43]-[44] +10; [45]→[55]. Map them by title, never by offset.
- **The registry's author labels are unreliable beside correct numbers** (2026-09-29: 10 of 24 Ref tokens mislabelled, e.g. "Mitsis (Ref. 15)" where [15] is Dahl-Sjostrand). Check a label against the reference list, not against the registry.
- The PDFs are in `scratch/literature/`; `pdftoppm -f P -l P -r 110` renders a page for a label check.

Related: [[reference-producer-landscape]].
