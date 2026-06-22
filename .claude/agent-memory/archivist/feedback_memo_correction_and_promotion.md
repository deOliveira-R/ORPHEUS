---
name: memo-correction-and-promotion
description: Pattern for correcting a stale planning memo (in .claude/agent-memory/) and promoting its catalog content to a durable Sphinx page. Banner-at-top discipline; replace_all gotcha when old shape appears in correction text; preserve historical quotes verbatim with context notes; design recommendations untouched; cross-link the Sphinx promotion back to the memo as design-side detail.
metadata:
  type: feedback
---

# Pattern: correct a stale planning memo + promote its catalog to Sphinx

When a planning memo (typically in
`.claude/agent-memory/explorer/` or
`.claude/agent-memory/method-implementer/`) predates an
architectural discovery that flipped a load-bearing decision
(e.g. Issue #196 storage-layout flip), the docs-side cleanup PR
has two coupled deliverables:

1. **Correct the memo in place** — preserve the narrative,
   replace shape arithmetic / boolean directions / etc. The memo
   stays in agent-memory; it's the design-side detail the next
   refactor session reads.
2. **Promote the memo's catalog to Sphinx** — typically as a
   new H1 section in a related theory page. The Sphinx version
   is the durable, build-versioned reference; the memo retains
   the design rationale + open questions + Issue # linkages.

**Why both:** memo-only is fragile (lives in agent-memory,
not auto-discovered by Nexus on Sphinx build); Sphinx-only loses
the design-history nuance (open questions, risk register,
rejected alternatives). Coupling the two with a cross-link
chain solves both.

**How to apply (the correction half):**

1. **Banner at the top of the memo body**, immediately after the
   H1 title.  3 paragraphs: (a) what predated, (b) the
   correction summary, (c) what is the single source of truth
   now. Date-stamped.  Example structure:

   ```markdown
   # <H1 title>

   > **Correction (YYYY-MM-DD).** This memo predates the layout
   > discovery that drove Issue #NN...  The original proposal
   > <verbatim wrong choice> was wrong; the principled <decision>
   > is <correct choice>...
   >
   > Design recommendations (<list>) are conceptually
   > unchanged — only the <thing that changed> moved.
   >
   > The current production state under the principled layout is
   > the single source of truth — see `<sphinx page>` "<promoted
   > section>" and "<related section>".
   ```

2. **Sweep shape / direction / sign mentions** with `replace_all`.
   **Gotcha:** the banner you just wrote ALSO contains the
   old-shape mention (in the "original was wrong" sentence).
   `replace_all` will overwrite it, leaving a tautology
   ("the original `(N, ng, nx, ny)` was wrong; the principled
   layout is `(N, ng, nx, ny)`"). **Always re-read the banner
   after the replace_all sweep** and explicitly preserve the
   old-shape mention via a separate Edit (or stage the banner
   AFTER the sweep). Bug discovered Issue #196 PR-CLEANUP-DOCS
   when initial banner read "original `(N, ng, nx, ny)` was
   wrong; principled is `(N, ng, nx, ny)`" — caught at the
   verification re-read.

3. **Rewrite dataclass / code snippet arithmetic** for the new
   axes.  `__post_init__` validation, property getters that
   index into `.values.shape[…]`, `np.einsum` strings, slice
   expressions — all touched. Don't half-do it; an honest memo
   has internally-consistent code snippets.

4. **Preserve verbatim quotes from external issues** (e.g.
   Issue #197 NewType definitions). Add a context note above:
   "quoted verbatim; note the shapes here predate Issue #NN's
   layout flip and are kept for historical accuracy." Future
   sessions can then cross-check what the original issue said
   vs the migration that intervened.

5. **DO NOT rewrite design recommendations.** The dataclass APIs,
   operator-algebra coupling, partial-close-of-Issue-N logic,
   risk register — all stay. Only the shape / direction / sign
   ARITHMETIC moves. If a design recommendation became
   structurally wrong after the discovery, that's a separate
   refactor PR, not a doc-correction PR.

**How to apply (the promotion half):**

1. **New H1 section** in the related Sphinx theory page. Pick
   the page that's the canonical home for the architectural
   decision (here: `docs/theory/index_convention.rst`). NOT a
   new page (per typical PR's "no new pages" anti-rec).
2. **Insert position matters.** Place the new vocabulary section
   BEFORE the per-array reference table that consumes the
   vocabulary, so a reader sees the conceptual map first.
3. **Target LoC range** (200-400 for a typical 7-table catalog
   per directive).  Each table gets a 2-3 sentence prose
   introduction explaining its epistemic role.
4. **Cross-link "Existing counterpart" cells** to actual code via
   `:class:`, `:meth:`, `:func:`, `:mod:`. Aspirational entries
   (not-yet-implemented) get plain text "None" or "Future".
5. **Add ``.. _section-anchor:`` BEFORE the H1 header** to give
   the section a referenceable name. Wire the corrected memo's
   "see Sphinx" pointer to this anchor.
6. **Cross-link the Sphinx promotion back to the memo** from the
   Future Work / typed-field-resume subsection. Bidirectional:
   memo points at Sphinx as canonical; Sphinx points at memo as
   design-side detail.

**Anti-patterns to avoid:**

- DO NOT mass-replace shapes WITHOUT re-reading the banner — see
  banner-tautology gotcha above.
- DO NOT promote with `automodule` if the memo content is
  hand-written prose tables. The promotion is a manual
  authored section, not a code-extracted reference.
- DO NOT delete the memo after promotion. The memo still owns
  the design rationale + Issue # linkages + risk register; the
  Sphinx promotion is a strict subset (the catalog only).

**Verification gates:**

- `grep -c "<old shape>" memo` must equal 0 after the sweep.
- The new H1 section's header must be greppable in the Sphinx
  page.
- Sphinx `-W` must exit 0 — the new section's cross-refs all
  resolve, no malformed table markup.

Worked example: Issue #196 PR-CLEANUP-DOCS corrected
`typed_field_contracts_for_phase_g.md` for the
`(N, ng, nx, ny)` layout and promoted its §1 inventory (7
catalog tables) to the new "SN Field Vocabulary" H1 section in
`docs/theory/index_convention.rst` (362 LoC). The
"Future Work / Typed-field contract resume" subsection was
rewired to point at the new section + the corrected memo as
the design-side detail.
