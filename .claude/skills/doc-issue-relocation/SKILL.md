---
name: doc-issue-relocation
description: 'Use when relocating a failed-experiment / dead-end / archived investigation narrative out of a Sphinx theory page into a GitHub issue thread. Covers the 5-block comment template, RST→GitHub-Markdown conversion, OPEN-investigation-log vs CLOSED-post-mortem framing, and the ≥900-LoC stub pattern that preserves externally-cited :label: anchors when the cut spans a whole H1 section. Examples: "move §X investigation to issue #N", "this 700-line wrong-attempt narrative belongs in the close-out comment, not in Sphinx", "preserve the labels but cut the section".'
---

# doc-issue-relocation — moving narrative from Sphinx into GitHub issue threads

When a Sphinx section records a failed experiment, a falsified extension,
or an archived architectural attempt, the body usually does not belong in
the production theory page. But the *content* is load-bearing for any
future agent who might re-attempt the same investigation. This skill
prescribes the exact workflow for **moving** that content (never
paraphrasing) into GitHub issue comments while keeping cross-references
inside the doc tree intact.

## When to use

- A Sphinx section is tagged ARCHIVED, FALSIFIED, RETRACTED, or carries
  a "wrong-attempt history" framing AND its cleanup is gated on closing
  or progressing a GitHub issue.
- The cut block is large enough (≥ 50 LoC) that inline deletion would
  lose context the issue thread cannot reconstruct on its own.
- The cut block contains numerical evidence, derivations, or
  literature-synthesis tables that must remain reproducible by a future
  agent.
- A doc cleanup commit (e.g. post-#138-style sweeps) is restructuring
  multiple sections and you need a durable home for the relocated
  material.

## When NOT to use

- The section is **production decision** content that belongs in
  Sphinx — never relocate active theory.
- The cut is < 50 LoC and the issue body already captures the
  context — a one-line tombstone in the issue is enough; full
  relocation is overkill.
- The narrative is fully captured in `tests/l0_error_catalog.md`
  (ERR-NNN entries with failure-mode, hide-mechanism, fix). The
  catalog IS the post-mortem; relocating to the issue duplicates it.
- The work is still active and the doc section is the active research
  record — keep it in Sphinx until the investigation closes or the
  issue is the better home.

## The 5-block GitHub comment template

Every relocation comment follows this exact structure. Save the draft to
`.claude/scratch/relocation/<issue>_<slug>.md`, post via
`gh issue comment <N> --body-file …`, then commit the doc cut.

```
# Issue #NNN <verb> — <topic>
**Date:** YYYY-MM-DD
**Doc cleanup commit:** <!-- COMMIT-HASH (filled in after the doc commit lands) -->
**Source relocated from:** <path>:<L1>–<L2> [+ secondary range].

[Banner — pick one:
 OPEN  → "Status: OPEN. This is the active research record — not a post-mortem."
 CLOSED → "Status: closed wontfix" / "closed-falsified" / "closed-superseded by #M".]

---

## Summary
- 3 bullets: what was tried, what killed it, what the production decision became.

Surviving doc anchors:
- :label: <eq-labels that survive elsewhere in the doc>
- :ref:   <section-labels that survive elsewhere in the doc>
- citations preserved as plain text (Sanchez1982, Hébert3-323, …)

[For OPEN issues only: "Hard evidence pin: ERR-NNN in tests/l0_error_catalog.md" if applicable.]

---

[Optional block: ## Production decision (preserved in Sphinx) — verbatim
copy of the decision that stays in the doc, with table. Repeat verbatim
to anchor the comment to what's still live.]

<details><summary>Full investigation record (verbatim from <path>:<L1>–<L2>)</summary>

[Verbatim RST content with the conversions in the table below.
 DO NOT paraphrase. Move; do not rewrite.]

</details>

---

## Hand-off / Closing posture
- For OPEN issues: numbered list — what the next agent must read first,
  what tooling exists, what is still unanswered.
- For CLOSED issues: closing posture stating the production verdict
  and any forward-pointing tracker.

References:
- Plain-text citation list. No RST `:cite:` directives — they don't render.

<!-- Source: <doc-path>:<L1>–<L2>, <plan-path> if relocating from a plan file -->
```

## RST → GitHub-Markdown conversion table

GitHub renders KaTeX and pipe tables but does NOT understand RST
directives. Convert the verbatim block element-by-element.

| RST construct                                    | GitHub-Markdown rendering                                                                                         |
| ------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------- |
| `.. math:: :label: foo` (block)                  | Triple-backtick fence with language `math`. Put `(label: \`foo\`)` on the line *after* the fence — KaTeX ignores `:label:`. |
| `:math:\`x\`` (inline)                           | `$x$`. Avoid `$$x$$` for display math inside lists — flaky; prefer the fence.                                     |
| `:func:\`module.fn\``                            | `` `module.fn` ``. Flag if the function moved/renamed in the cut.                                                 |
| `:class:\`Foo\``, `:meth:\`Foo.bar\``            | `` `Foo` ``, `` `Foo.bar` ``                                                                                      |
| `.. list-table::` w/ `:header-rows: 1`           | GitHub pipe table. Multi-paragraph cells flatten with `<br>` or a sibling `*notes*` row.                          |
| Multi-line `* -` rows                            | Collapse to single-line; sub-bullets become `<br>·` separated runs.                                               |
| `[CitationKey]_` references                      | Plain-text `[CitationKey]` (no underscore). Repeat the bibliography under `## References:`.                       |
| `:ref:\`some-anchor\``, `:eq:\`label\``          | Backticked anchor name + a sentence describing what the anchor names; cite the file path: `` `docs/theory/peierls_unified.rst` ≈ lines 5410–5476 ``. |
| `.. note::`, `.. warning::` admonitions          | Markdown blockquote `> **Note.**` / `> **Warning.**` followed by the body indented under the quote.               |
| `**bold**` / `*italic*`                          | Stay as-is — both renderers agree.                                                                                |
| `.. code-block:: python`                         | Fenced ``` ```python ``` block.                                                                                    |
| `:sup:\`-3\``, `:sub:\`f\``                      | `<sup>-3</sup>`, `<sub>f</sub>`. KaTeX-internal exponents stay inside `$…$`.                                      |
| `\Sigt{}`, `\nSigf{}`, `\keff` (LaTeX macros)    | These are project macros from `docs/conf.py` — KaTeX in GitHub does NOT load them. Inline-expand to `\Sigma_t`, `\nu\Sigma_f`, `k_\text{eff}`. |

### Worked conversion example

RST source:

```rst
.. math:: \int_0^1 e^{-\tau/\mu} d\mu = E_2(\tau)
   :label: peierls-slab-E2

The escape kernel reduces to :math:`\frac{1}{2} E_2(\tau)` — see
:func:`compute_P_esc_outer` and the surviving anchor
:ref:`theory-peierls-slab-polar`.
```

GitHub-Markdown:

````markdown
```math
\int_0^1 e^{-\tau/\mu} d\mu = E_2(\tau)
```
(label: `peierls-slab-E2`)

The escape kernel reduces to $\frac{1}{2} E_2(\tau)$ — see
`compute_P_esc_outer` and the surviving doc anchor
`theory-peierls-slab-polar` in `docs/theory/peierls_unified.rst`.
````

## Collapsibles for ≥ 50-LoC verbatim blocks

GitHub renders `<details><summary>…</summary> … </details>` cleanly.
Use it whenever the verbatim record exceeds ~50 LoC. The `<summary>`
line is the only thing every reader sees first; the body is what an
agent reconstructing context unfolds.

```html
<details><summary>Full investigation record (verbatim from peierls_unified.rst:5410–5476)</summary>

[…verbatim block, converted per the table above…]

</details>
```

Multi-section cuts go in **a single comment** with multiple
collapsibles, not multiple comments. The issue is the unit of memory;
fragmenting across comments defeats the relocation goal. The single
exception is the parent/child split for multi-issue audit tables (see
the close-out narrative arc in the archivist agent).

## OPEN-investigation-log vs CLOSED-post-mortem framing

The framing is load-bearing — it tells the next agent whether the
issue is a frozen archaeological record or an active hand-off.

| Aspect             | OPEN (investigation log)                                                                | CLOSED (post-mortem)                                              |
| ------------------ | --------------------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| Title verb         | "Investigation log", "active research record"                                           | "close-out", "post-mortem", "wontfix"                             |
| Banner             | `Status: OPEN. This is the active research record — not a post-mortem.`                 | `Status: CLOSED — falsified` / `wontfix` / `superseded by #M`     |
| Hard evidence pin  | Reference an `ERR-NNN` entry if one exists; the catalog is the regression gate.         | Optional — the ERR-NNN may have already migrated to a fix commit. |
| Tail section       | `## Hand-off for the next agent` — numbered list of what to read first and what tooling | `## Closing posture` — single paragraph, production verdict       |
| xfail-strict pin   | Explicitly mention any xfail-strict test that flips on fix-landing                      | Not applicable                                                    |
| Sibling cross-link | If the issue is a partial close-out of a sibling-closed issue, lead with that link      | Optional cross-link to follow-up tracker                          |

## The ≥ 900-LoC stub pattern

When the cut block is a whole H1 section (typically ≥ 900 LoC, e.g.
"Moment-form Nyström assembly — ARCHIVED" in `peierls_unified.rst`)
**and** its `:label:` anchors are referenced from elsewhere in the
doc tree, do NOT delete the section. Leave a ~25-LoC stub that
preserves the anchors and points at the durable archive.

### Step 1 — inventory references BEFORE editing

```bash
grep -rn "label-A\|label-B" docs/theory/
```

Anything inside the cut block: fine, will be deleted with the rest.
Anything OUTSIDE the cut block: the citing site needs the label
preserved in the stub OR rewritten inline. The label-survival decision
is what determines whether you stub or inline-rewrite.

### Step 2 — write the stub

```rst
.. _label-A:

.. _label-B:

Section title — ARCHIVED
========================

<1-paragraph reason-for-archive: "no longer needed because X
replaces it; preserved for future Y if Z changes">

**Production decision.** <what is the active path now>

**Where the code lives.**

- :file:`derivations/archive/foo.py` — short purpose
- :file:`derivations/archive/bar.py` — short purpose

**Where the close-out narrative lives.** See
`Issue #NNN <https://github.com/.../issues/NNN#issuecomment-XXXX>`_.

<1-line note about what's preserved at each anchor for future searchers>
```

### Step 3 — rewrite citing sites

External `:ref:`s often imply the cut content is "below" or "in this
doc" — phrasing like "see :ref:\`X\` for the full defect analysis and
the diagnostic tables." After the cut, the diagnostic tables live in
the GitHub issue, not in the doc. Rewrite the citing site:

| Before                                                                                  | After                                                                                                                                  |
| --------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------- |
| `See :ref:\`X\` for the full defect analysis and the diagnostic tables.`                | `The full defect analysis and diagnostic tables now live in [Issue #NNN](URL); the in-tree pointer is :ref:\`X\`.`                     |
| `See Phase F.5 investigation`                                                           | `See Phase F.5 close-out for the X-reference synthesis, structural obstruction, and production decision`                               |

### Step 4 — verify per sub-step, not at the end

Validate `python -m sphinx -W -b html docs docs/_build/html` after
*every* commit in the relocation series, not just at the end. A broken
`:ref:` from a missed citing site shows up as a build warning fixable
in the same sub-step's commit (or as a follow-up commit
`docs(<scope>): #NNN — fix dangling :ref: from sub-step N`, never via
`git amend`).

### Step 5 — anchor sanity check

```bash
grep -c "^\.\. _" docs/theory/<file>.rst
```

Run before/after the cut. The delta equals the number of `:eq:`
labels INSIDE the cut block that you intentionally let go (their
`:eq:` references must also have been inside the cut block — confirm
with grep). The auto-generated `verification/matrix.rst` lists them
as orphan equations; fine, that file regenerates on every build and
the orphan list shrinks accordingly.

## Pitfalls

1. **Block-math conversion.** `.. math:: :label: foo` blocks become
   triple-backtick-`math` fenced blocks. Preserve the label as a
   separate `(label: \`foo\`)` line *outside* the fence.

2. **Multi-row list-tables.** RST `.. list-table::` with multi-line
   `* -` rows must collapse to single-line markdown rows.

3. **Cross-doc references.** RST `:ref:\`anchor\`` and `:eq:\`label\``
   cannot resolve in GitHub markdown. Convert to backticked anchor
   names with a sentence describing what they point to. Cite the file
   path explicitly: `` `docs/theory/peierls_unified.rst` ≈ lines 5410–5476 ``.

4. **Citation entries.** RST `[CitationKey]_` reference blocks must be
   repeated as plain-text lines under `## References:`. The GitHub
   viewer does not back-resolve them.

5. **OPEN vs CLOSED framing matters.** Use "investigation log" for
   open, "close-out" / "post-mortem" for closed. The banner and tail
   section change accordingly.

6. **Line-range drift.** Plan-handed line ranges captured pre-cleanup
   hold within ±5 lines unless an unrelated commit landed between
   sessions. If the section title is named in the plan, grep that
   title to re-derive the range. Titles are stable; line numbers are
   not.

7. **Don't paraphrase.** The wording in the doc was chosen carefully
   at the time of the experiment — often by the literature-researcher
   or numerics-investigator who ran the falsification. Move; do not
   rewrite.

8. **Preserve numerical tables verbatim.** Pre-fix vs post-fix tables,
   signed-error scans, sphere/cyl convergence patterns are PRODUCTION
   DATA, not narrative. Even when the surrounding prose is
   wrong-attempt history, the tables are usually load-bearing
   evidence for the OPEN/CLOSED status of the issue.

9. **Collapsibles aggressively for >50-LoC blocks.** The 3-line
   summary at top is what every reader sees first; the verbatim
   body is what an agent reconstructing context reads.

10. **One comment per issue, not one per section.** Even when the
    source spans multiple sections, keep them in a single comment
    with separate collapsibles. The issue is the unit of memory;
    multi-comment fragmentation defeats the relocation goal.

11. **Cross-link sister issues.** Issues that share an obstruction
    (e.g. the same Lambert/Marshak basis-rotation algebra at N≥2)
    form a connected web. Each comment must explicitly link to its
    siblings — that is how a future agent reconstructs the *family*
    of obstructions, not just one of them.

12. **Counter-pattern (avoid):** "Just delete the section, the issue
    comment has it." — leaves dangling `:ref:`s and invalidates
    external code/test references that may still link to the labels.
    Cost of a 25-LoC stub is far less than the cleanup cost.

## Output checklist

After the relocation lands:

- [ ] Comment posted to the right issue with the 5-block structure.
- [ ] Doc cleanup commit references the comment URL in the commit body.
- [ ] All externally-referenced `:label:`s preserved in the stub OR
      every citing site rewritten.
- [ ] `python -m sphinx -W -b html docs docs/_build/html` passes (or
      warning count unchanged from baseline).
- [ ] `grep -c "^\.\. _" docs/theory/<file>.rst` delta accounted for.
- [ ] Citing sites that said "see X below for tables" rewritten to
      point at the GitHub issue.
- [ ] Sibling issues cross-linked when the obstruction is shared.
- [ ] HTML comment marker `<!-- Source: <path> -->` left in the
      comment for future grep-archaeology.
