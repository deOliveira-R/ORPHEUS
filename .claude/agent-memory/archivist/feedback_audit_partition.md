---
name: Audit-then-edit doc cleanup partitions
description: Pattern for read-only audit reports that classify doc sections as KEEP/RELOCATE/TRIM/REMOVE before any editing
type: feedback
---

When the user asks for a doc-cleanup audit (typically before pushing a refactor like #138 that retires symbols), produce a **partition table** structured exactly like:

```
### Peierls section inventory       (one-line characterization per section, with line ranges)
### KEEP                            (table: lines | reason)
### RELOCATE                        (table: lines → issue # | summary)
### TRIM / REMOVE                   (table: lines | action | reason)
### Stale cross-refs to retired symbols  (table: line | stale text | replacement)
### Estimated total LoC reduction   (Before/After/Δ; per-section breakdown)
### Cross-cutting follow-ups        (issues to open if not already tracked)
```

**Why:** The user's directive is "doc reflects CODE in detail, but NOT extensive failed-experiment narratives." Failed-experiment narratives belong in the GitHub issue that started them — the issue title plus the agent's git log already houses the why. The doc should be tight: derivation, current code refs, current numerical evidence, structural gotchas that future readers must see.

**How to apply:**

1. **Always read the current code's canonical API** (e.g., `solve_peierls_1g(geometry=…, boundary=…)`) so stale cross-refs are found and the suggested replacement is concrete, not speculative.
2. **Do not edit during an audit** — the user explicitly asks for a report. The parent agent (or a follow-up edit pass) executes the changes. Returning a clean partition lets the user merge/reject by section.
3. **Distinguish three categories of "narrative"**:
   - **Load-bearing teaching** (R-vs-R² gotcha, "why polar form not chord form", four-strategy taxonomy for log singularities) → KEEP. These prevent future readers from making the same architectural mistake.
   - **Tracked-work narrative** ("Phase-4.2 item C8 is a prerequisite for…", "white-BC closure not yet implemented") → trim to one sentence + GitHub issue ref.
   - **Session-history / wrong-attempt narrative** ("Phase-4.1 debugging insight", "the Zotero MCP server was unreachable", "earlier conclusion is retracted") → relocate to issue or `peierls_unified` (the unified page is the natural home for retraction-in-context). Doc retains a 1-3 line pointer.
4. **Numerical-evidence tables (k_eff scans, row-sum residuals, convergence tables) are Cardinal Rule 2** ("numerical evidence: before/after with multiple cases"). Always KEEP them, even when the surrounding narrative is the wrong-attempt story — move the tables into the **current** subsection if necessary.
5. **Cross-doc cross-refs to retraction sections work** — `:doc:`peierls_unified` §8` is enough; do not re-state the wrong claim in the new page just so you can retract it. The retraction lives once, where the deep architectural narrative lives.
6. **Phase tags ("Phase-4.2", "Phase B.4", commit hashes in parens) are session-internal noise** in a theory page. Strip them; the cross-ref to a labelled section/issue is enough. Git log preserves commit attribution.
7. **Estimate LoC reduction per section** so the user knows where the bulk of the cleanup lives. The bulk is almost always in the most recent retraction/postmortem section, because that's where session-history accumulates fastest.

**Pitfall to avoid:** Do not confuse "retired symbol cross-ref" (must be replaced with the canonical post-refactor symbol) with "retired-attempt narrative" (move to issue). The first is mechanical; the second is editorial.

**Pitfall observed in this audit:** `seealso` blocks at the end of long sections often have stale wording ("thin facade") that pre-dates the registry-only refactor — the `:mod:` directive itself resolves correctly because the module still exists, so Sphinx does not warn. Grep for "thin facade", "wrapper", "trivial helper" specifically when the upstream refactor was a façade-collapse.
