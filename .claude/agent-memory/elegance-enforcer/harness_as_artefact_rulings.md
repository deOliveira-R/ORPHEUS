---
name: harness-as-artefact-rulings
description: Reviewing a PROSE corpus (rules/skills/CLAUDE.md) as code under Pattern 2 — the five probes, the forced-copy finding, and what counts as a duplicate vs a citation.
metadata:
  type: project
---

Reviewing the agent harness itself (13 rule/skill/onboarding pages, W5/R1 of the
`harness_context_budget` campaign, 2026-09-21). Memo:
`scratch/_harness_eval/reeval/elegance_memo.md`. The durable part is the METHOD —
a prose corpus obeys Pattern 2 exactly like code, and these probes transfer to any
future harness, docs-architecture or rule-consolidation review.

**Why:** the user ruled consolidation-before-retirement, verbatim: *"look first at
things that become stronger as they are together."* That reframes a duplicate
census: the deliverable is not "which copy dies" but "what does the merged
definition carry that no copy carries today" — and for 18 of 25 rows the answer
was substantial, so a retirement-first pass would have deleted load-bearing halves.

**How to apply — the five probes for a prose-corpus Pattern-2 review:**

1. **State the STATED-vs-CITED predicate before counting.** A site counts only if it
   carries the imperative, `check:`, `tell:` or mechanism *in its own words*; a
   `[case]` link, a "see X" or a gate-list entry is a citation. Without this the
   regex census reports 5–7 sites where 2 are real. `[M]` here: the regex named
   bit-identity in 5 cores; reading found the statement in exactly 1.
2. **The FORCED COPY is the finding, not the duplicate.** Look for the one site whose
   reader *cannot load the definition* — here the brief template's "Rules that apply
   to you" block, read by three `omitClaudeMd` agents. A forced copy is a build
   product written by hand: the fix is a GENERATED include (the harness already
   injects `<!-- harness: error_index -->`), not a citation. It was already drifting
   (the page restates its own §W7 twenty lines later).
3. **A page whose preamble promises not to restate, and does.** `cardinal`'s preamble
   says "this page does not restate"; three of its five rules restate. That is X3
   (prose is not enforcement) applied to the harness — a prose claim with no line
   that fails when it is false. Check every page's self-description against its body.
4. **Intra-file duplicates hide behind cross-file ones.** `omitClaudeMd` is stated
   three times inside ONE always-on file (roles table, invariant 4, "Dispatch
   facts"). Count sites within a file, not only across files.
5. **Boilerplate preambles are free tokens.** 5 of 13 pages restate the page FORMAT
   (`check:`/`tell:`/`[case]`) and 5 restate the `[M]`/`[R]` gloss. A reader OBEYING
   a rule never needs the rule's own format; a reader WRITING one does — so the
   format belongs on the harness page, on demand. ≈1 100 always-on chars, zero agent
   actions depend on them.

**Do NOT collapse a genuine multi-face law.** Single-source-of-truth has four faces
by design: architecture (Cardinal Rule 2), code (`coding-elegance` Pattern 2),
evidence (`instrument-doctrine` X4), V&V (the three pillars). Each is load-bearing
in its own tier. What must stop is each face RESTATING the shared mechanism (the
α-normalised-AST decider; "name the mechanism that kept the copies equal and grep
its message for a witness"). Cite the mechanism, keep the face. Same ruling for the
declared floor/ceiling pairs (`coding-standards` § "A guard is elegance debt" vs
Pattern 4) — those are the MODEL, not a defect.

**Tiering heuristic that emerged:** a clause bites in a PHASE (a retirement, a
rename, a census, a battery, a commit) or under a PATH. Move the body out of
always-on and leave a one-line TRIGGER resident ("a delete, a rename or a re-home
loads the retirement audit"). The trigger is the only part that must be resident,
because the failure it prevents — a deleted symbol whose docstring references
nobody re-checks — surfaces silently, with no build warning at any severity.

See also [[lessons]] L13 (doc-carve certification) and the `coding-elegance` skill's
Pattern 2 / "a repeated conditional is a missing type" — here a repeated CLAUSE is a
missing DEFINITION, and the fix is one definition cited by ID.
