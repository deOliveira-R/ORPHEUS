---
name: Large relocation cuts (>900 LoC) preserve anchors as stub
description: How to execute a >900-LoC single-section cut from a doc when the section's labels are externally referenced; stub pattern with explicit anchor preservation
type: feedback
---

**For very large relocation cuts (≥ 900 LoC, often a whole H1 section), preserve every `:label:` anchor inside the cut block as part of the replacement stub.** A 25-LoC stub is enough if it carries: (a) all preserved labels, (b) production decision, (c) "where the code lives" file pointers, (d) "where the close-out narrative lives" — the GitHub issue close-out comment URL, (e) git-archaeology hint (`git log --grep=...`).

**Why:** Task 2026-04-30 Phase 2a peierls_unified.rst cleanup. Issue #138 close-out cut 935 LoC of moment-form Nyström architectural narrative (the archived τ-Laguerre + closed-form moment-recursion derivations). The labels `theory-peierls-moment-form` and `theory-peierls-moment-form-failed-polar` were referenced from FOUR other sites in the same doc. Wholesale deletion would have produced 4 broken `:ref:`s, each requiring an inline rewrite. Instead a 25-LoC stub preserved both labels AND the load-bearing pointers, so the cross-references continued to resolve to a useful target.

**How to apply:**

1. **Inventory references to the cut content's labels FIRST.** Before any edit:
   ```bash
   grep -n "theory-peierls-moment-form\|theory-peierls-moment-form-failed-polar" docs/theory/peierls_unified.rst
   ```
   Anything inside the cut block is fine (will be deleted with the rest); anything OUTSIDE the cut block needs the label preserved in the stub OR the citing site rewritten.

2. **Stub structure (~25 LoC):**
   ```rst
   .. _label-A:

   .. _label-B:

   Section title — ARCHIVED
   ========================

   <1-paragraph reason-for-archive: "no longer needed because X
   replaces it; preserved for future Y">

   **Production decision.** <what is the active path now>

   **Where the code lives.**
   - :file:`derivations/archive/foo.py` — short purpose
   - :file:`derivations/archive/bar.py` — short purpose

   **Where the close-out narrative lives.** See
   `Issue #NNN <https://github.com/.../issues/NNN#issuecomment-...>`_.

   <1-line note about what's preserved at each anchor>
   ```

3. **External `:ref:`s often imply the cut content is "below" or "in this doc"** — phrasing like "see :ref:`X` for the full defect analysis and the diagnostic tables". After the cut, the diagnostic tables are in the GitHub issue, not in the doc. The citing site needs to be rewritten to say so:
   - Before: `See :ref:`X` for the full defect analysis and the diagnostic tables.`
   - After: `The full defect analysis and diagnostic tables now live in [Issue #NNN](URL); the in-tree pointer is :ref:`X`.`

4. **Validate `Sphinx -W` after every sub-step**, not just at the end. A broken `:ref:` from a missed citing site shows up as a build warning that you can fix in the same sub-step's commit (or as a follow-up commit `docs(peierls): #NNN — fix dangling :ref: from sub-step N`, never via `git amend`).

5. **`grep -c "^\.\. _label:"` before/after** is the quick sanity check on anchor preservation. Phase 2a started at 72 labels, ended at 65 — the 7 lost labels were all `:eq:` labels INSIDE the cut moment-form section (peierls-moment-K-source-form, peierls-moment-segment, peierls-moment-contraction, peierls-moment-J-E1, peierls-moment-J-Ki1, peierls-moment-vandermonde, peierls-moment-K-assembly), and ALL their `:eq:` references were also inside the cut block. The auto-generated `verification/matrix.rst` lists them as orphan equations — fine, that file is regenerated on every build and the orphan list shrinks accordingly.

**Counter-pattern (avoid):** "Just delete the section, the issue comment has it." — leaves dangling `:ref:`s and invalidates external code/test references that may still link to the labels. Cost of a 25-LoC stub is far less than the cleanup cost.
