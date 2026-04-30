---
name: Phase 2d table-preservation constraint
description: When the plan says "preserve numerical tables" inside a section earmarked for relocation, the actual LoC delta is bounded by the table footprint, not by the plan's headline cut number
type: feedback
---

When a relocation plan reads "cut to ~25 LoC stub" but ALSO says
"preserve the numerical tables at lines X-Y", the actual achievable
delta is `(original LoC) − (table LoC) − (stub LoC)`, which is
much smaller than the headline number. The Phase 2d
collision_probability.rst Issue #100 cut was budgeted at −105 LoC
(130 → 25), but the preserved tables alone occupy ~70 LoC, so the
real delta was −30 LoC (130 → 100). This is correct, not a miss.

**Why:** Production numerical tables are load-bearing for current
test gates (`TestWhiteBCRank1ErrorScan`, `TestSphereWhiteBCRowSum`)
and the rationale-section directive ranks **table preservation**
above headline-LoC targets. The plan's envelope number was an
estimate written without table-LoC subtraction.

**How to apply:**

1. Before starting a relocation, count the LoC of the
   to-be-preserved tables (`awk 'NR==X,NR==Y' file | wc -l`).
2. Subtract from the section size to get the achievable narrative
   cut.
3. If the resulting delta is meaningfully smaller than the plan's
   target, state this in the commit message and final report
   ("plan target −X; actual −Y because tables preserved per
   plan's R3 risk-register clause") rather than treating it as a
   deviation requiring further cuts.
4. Reframe the section heading from a debate-narrative
   ("Issue #N — retraction of the X claim") to a
   production-evidence framing ("Rank-1 white-BC numerical
   evidence — sphere") so the table reads as canonical reference
   data rather than failed-experiment baggage. The retraction
   note becomes a paragraph at the end, pointing at the GH issue
   for the full debate.
5. Preserve the original `:ref:` anchor verbatim — even if the
   section title flips, callers from sister documents
   (`peierls_unified.rst:5332`, `:5449`) still cite the anchor
   and the new framing reads cleanly when their text is
   "see ... for the numerical evidence".
