---
name: reference-local-literature
description: Where the transport textbooks live locally for verifying equation numbers (OCR markdown in scratch, image-only PDFs elsewhere) and their page offsets
metadata:
  type: reference
---

Equation numbers are verified against the source, never recalled. The `literature-researcher`
(W7) owns acquisition; this is only where an archivist finds a book already on disk.

- **Bell & Glasstone 1970** (`BellGlasstone1970`): OCR markdown, greppable,
  `scratch/literature_ocr/Bell-Glasstone(1970)Nuclear_reactor_theory.md`. Each page carries
  `## p. N` (PDF page) and `*[printed header: <printed page> ...]*`; cite the PRINTED page.
- **Stamm'ler & Abbate 1983** (`Stamm1983`): chapters 4 and 6 only, OCR markdown beside it.
- **Duderstadt & Hamilton 1976** (`Duderstadt1976`): image-only PDF (pdftotext returns nothing),
  `~/Downloads/Books/1976_Nuclear_Reactor_Analysis.pdf`; read pages with Read `pages=`.
  `[M]` 2026-10-02: PDF page = printed page + 17 (Ch. 4 §IV, the P1/Σ_tr derivation, is printed
  pp. 124–139).
- **Not on disk** `[M]` 2026-10-02: Case & Zweifel 1967 (in `refs.bib`, unverifiable here),
  Cercignani (not in `refs.bib`). Name them in `NEEDS:`, never cite an unverified equation number.
- `scratch/literature/` holds ~100 more PDFs (papers); list it before asking for a source.
