---
name: Issue close-out relocation comments (failed-experiment narrative)
description: Pattern for drafting GitHub issue comments that relocate failed-experiment narrative out of Sphinx into issues; covers RST→GitHub-Markdown math, anchor preservation, OPEN vs CLOSED framing, line-range drift handling
type: feedback
---

When relocating ~50–500-LoC RST blocks from theory pages to GitHub issue comments (post-#138-style cleanup), use this structure verbatim. It maximises reconstruction-from-issue-alone for a future agent.

**Why:** Cardinal Rule 1 — Sphinx is the LLM's brain — only holds if the relocated knowledge lives somewhere durable. A future agent reading the issue alone must reconstruct what was tried, why it failed, and what the production decision became, without needing the doc.

**How to apply:** Every comment file follows this exact 5-block template.

## Template

```
# Issue #NNN <verb> — <topic>
**Date:** YYYY-MM-DD
**Doc cleanup commit:** <!-- COMMIT-HASH -->
**Source relocated from:** <path>:<L1>–<L2> [+ secondary range].

[For OPEN issues: "Status: OPEN. This is the active research record — not a post-mortem." + frame as hand-off to next agent.]
[For CLOSED issues: "Status: closed wontfix" / "closed-falsified" / etc.]

---

## Summary
- 3 bullets: what was tried, what killed it, what the production decision became.

Surviving doc anchors:
- :label: <eq-labels that survive elsewhere>
- :ref: <section-labels that survive elsewhere>
- citations preserved as plain text (Sanchez1982, Hébert3-323, etc.)

[For OPEN: "Hard evidence pin: ERR-NNN in tests/l0_error_catalog.md" if applicable.]

---

[Optional: ## Production decision (preserved in Sphinx) — verbatim copy of the decision that stays in the doc, with table.]

<details><summary>Full investigation record (verbatim from <path>:<L1>–<L2>)</summary>

[Verbatim RST content with these conversions:
 - .. math:: :label: foo  →  ```math ... ```\n(label: `foo`)
 - :math:`x`  →  $x$
 - :func:`x.y`  →  `x.y`
 - .. list-table:: with :header-rows: 1  →  GitHub markdown pipe table
 - **bold** stays
 - DO NOT paraphrase. Move; do not rewrite.]

</details>

---

## [Hand-off / Synthesis section]
- For OPEN issues: numbered list of what the next agent must read first, what tooling exists, what is still open.
- For CLOSED issues: closing posture stating the production verdict and any forward-pointing tracker.

References:
- Plain-text citation list. No RST :cite: directives — they don't render in GitHub.
```

## Pitfalls

1. **Block-math conversion.** `.. math:: :label: foo` blocks become triple-backtick-`math` fenced blocks. The label is preserved as a separate `(label: \`foo\`)` line *outside* the fence, because GitHub's KaTeX renderer ignores `:label:`-style attributes. Inline `:math:\`x\`` becomes `$x$`. Display math without a label works as `$$x$$` or as a `math`-fenced block; prefer the fence for multi-line math because GitHub's `$$ ... $$` is finicky inside lists.

2. **Multi-row list-tables.** RST `.. list-table::` with multi-line `* -` rows must collapse to single-line markdown rows. Multi-paragraph cells (with sub-bullets) must be flattened to a single cell with `<br>` or a `*notes*` style adjunct row.

3. **Cross-doc references.** RST `:ref:\`some-anchor\`` and `:eq:\`label\`` cannot resolve in GitHub markdown. Convert to backticked anchor names with a sentence describing what they point to. Cite the file path explicitly: `\`docs/theory/peierls_nystrom.rst\` ≈ lines 5410–5476`.

4. **Citation entries.** RST `[CitationKey] ...` reference blocks must be repeated as plain-text lines under "References:" — the GitHub viewer does not back-resolve them.

5. **OPEN vs CLOSED framing matters.** For OPEN issues, use "investigation log" not "post-mortem" in the title; lead with `Status: OPEN. This is the active research record — not a post-mortem.`; close with a `## Hand-off for the next agent` section listing what they should read first. For CLOSED issues, lead with "close-out" and end with a "Closing posture" paragraph.

6. **Line-range drift.** Plan handed line ranges captured pre-cleanup; they will hold within ±5 lines unless an unrelated commit landed between sessions. If the section title is named in the plan (e.g. "Villarino-Stamm'ler per-mode extension (novel, falsified ...)"), grep that title to re-derive the range. The titles are stable; the line numbers are not.

7. **Don't paraphrase.** The wording in the doc was chosen carefully at the time of the experiment (often by the literature-researcher or numerics-investigator who ran the falsification). Paraphrasing risks degrading semantic precision. Move; do not rewrite.

8. **Preserve numerical tables verbatim.** Pre-fix vs post-fix tables, signed-error scans, sphere/cyl convergence patterns are PRODUCTION DATA, not narrative. Even when the surrounding prose is "wrong-attempt history", the tables are usually load-bearing evidence for the OPEN/CLOSED status of the issue.

9. **Collapsibles aggressively for >50-LoC investigation records.** GitHub renders `<details><summary>...</summary>` cleanly. The 3-line summary at top is what every reader sees first; the verbatim block is what an agent reconstructing context reads.

10. **One-comment-per-issue, not one-comment-per-section.** Even when the source spans multiple sections (e.g. Comment 5 — #132 spans both the Probe G cascade narrative AND the follow-up directions), keep them in a single comment with separate collapsibles. The issue is the unit of memory; multi-comment fragmentation defeats the relocation goal.

11. **Cross-link sister issues.** Issues #119 / #122 / #123 / #132 form a connected web (all trace back to the same Lambert/Marshak basis-rotation algebra at $N \ge 2$). Each comment should explicitly link to its siblings — that is how a future agent reconstructs the *family* of obstructions, not just one of them.
