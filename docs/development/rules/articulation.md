---
harness:
  kind: rule
  budget_tokens: 800
  brief: >-
    a report is complete sentences; `[M]` measured with its command, `[R]` reasoned, `[HYPOTHESIS]` proposed; a bare number is read as measured; `[M]` on an inherited claim certifies that some measurement answered some question, so a `[M]` on a negative (absent, discarded, no consumers) is re-measured against the question at hand before it is built on.
---

# Articulation — the writing standard

Articulate means: take a concept apart so that a reader with none of your
context reassembles it losslessly. The loss is the measure. For a plan the loss
appears as a surprise — a later session believes something the tree
contradicts. For a doc it appears as a wrong reconstruction or a dead
reference. For a brief it appears as a sub-agent repairing the instruction
silently.

1. **Say what you mean.** When a literal phrase is available, use it. Mannered
   prose — metaphor and flourish in place of direct statement, "a dial worth
   turning" for "a parameter worth varying" — drags in connotations you did not
   choose and costs context for nothing.
2. **Straight to the point is not clipped.** One idea per sentence. Complete
   sentences in anything a reader keeps: a plan, a doc, a report, an issue.
   Shorthand only between tool calls.
3. **Write in the reader's vocabulary, not the session's.** A name you coined
   while working is defined at first use or not used. Expand an uncommon
   acronym once.
4. **Every named object is defined in the text or linked.** A symbol whose
   spelling is ambiguous is named by its structure — what it ranges over and
   what it is paired with. A file, function or flag is named only when the
   reader must go there.
5. **Markers carry epistemic status, one meaning per spelling:** `[M]`
   measured, with the command or date; `[R]` reasoned, not yet measured;
   `[HYPOTHESIS]`; `[REFUTED YYYY-MM-DD]`; `[LANDED <hash>]`;
   `[REMEDIED YYYY-MM-DD @<hash>]`. A bare number with no marker is read as
   measured. No decorative glyphs; if a glyph is used anywhere, it has one
   meaning everywhere.
6. **The two hottest cases.** A *plan* is a message to yourself or a sub-agent
   after context is gone: `plan-authoring` governs its claims, this rule governs
   its prose. *Documentation* is the project's brain: maximal effort on
   completeness and derivation, never on brevity, and a present-tense-false
   sentence in a doc is a bug fixed on sight.

- check: show the text to a colleague with none of the task's context — if they
  would be confused, so will the model. For a plan: can every step be executed
  from the text alone, without asking the author?
- tell: the reader meets an arrow chain, a compound stacked from hyphens, a
  label invented mid-session and used as if it were shared, or a paragraph
  that must be read twice.
