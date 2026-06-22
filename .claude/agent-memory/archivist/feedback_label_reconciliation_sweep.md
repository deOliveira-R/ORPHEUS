---
name: label-reconciliation-sweep
description: Pattern for resolving Sphinx :verifies: missing-label warnings on a docs-side cleanup PR. Three resolution paths (add label to existing eq, retarget decorator, add new labelled eq block); section-vs-equation label disambiguation via *-section suffix; in-scope vs out-of-scope triage by module label; Sphinx info-vs-warning severity distinction matters for the -W gate.
metadata:
  type: feedback
---

# Pattern: ``@pytest.mark.verifies`` missing-label reconciliation sweep

When a test suite carries many ``@pytest.mark.verifies(LABEL)``
decorators that point at non-existent ``.. math:: :label: LABEL``
blocks, the Sphinx Nexus extension emits one
"``no matching equation node math:equation:LABEL — skipping``"
info line per (decorator × test) site.  These appear in the
Sphinx output and (without -W) look like warnings — but per the
Nexus extension they emit at info severity, NOT warning severity,
so ``sphinx-build -W`` accepts them.

**Why:** the typical PR-CLEANUP-DOCS deliverable for an
architectural migration (e.g. Issue #196 PR-INDEX-1..7) accumulates
these warnings as test decorators outrun docs growth.  The
reconciliation sweep is the dedicated docs-PR that closes the
gap.

**How to apply:**

1. **Get the unique label list** via
   ``sphinx-build -q docs ... 2>&1 | grep "no matching equation node"
   | sed -E "s/.*verifies\('([^']+)'\).*/\1/" | sort -u``.
2. **Triage IN-SCOPE vs OUT-OF-SCOPE by module label** of the
   carrying test file.  An SN-architecture PR closes SN-label
   warnings; CP / FN / derivations labels are out-of-scope and
   recorded as such in the closeout (typical ratio: 1:1 in-scope
   vs deferred).
3. **For each in-scope label, pick a resolution path:**
   - **Path 1 — ADD ``:label:`` to existing eq.**  Use when the
     equation is already in prose / a code block in a theory page
     but missing the ``:label:`` directive.  Cheapest fix; the
     equation is already documented.
   - **Path 2 — UPDATE the decorator** to the existing label name.
     Use only when the decorator is a typo / rename and the eq
     exists under a different name.  Verify via the V&V audit
     output: ``python -m tests._harness.audit``.
   - **Path 3 — ADD new ``.. math:: :label:`` block**.  Use when
     the equation is genuinely missing from the theory; the
     decorator pointed at a future-planned label.  This is the
     most write-intensive path — surround the new eq with the
     prose context (motivation, citation, code pointer) per
     Cardinal Rule 3 (Sphinx as brain).

4. **Section-label vs equation-label disambiguation.**  If a label
   exists as a section anchor (``.. _foo:``) but the decorator
   expects an equation node, the right move is to **rename the
   section anchor to ``foo-section``** (add the suffix) and
   **create a new ``.. math:: :label: foo`` block** at the head of
   the renamed section.  ``vv-principles`` and the Nexus extension
   both treat ``:eq:`` (equation) and ``:ref:`` (section) as
   distinct node types — the same name can sit on both with the
   suffix discipline.  Always grep the codebase for the old anchor
   first and update every ``:ref:`` consumer to the ``-section``
   variant in the SAME commit.

5. **Severity discipline for the -W gate.**  Run
   ``sphinx-build -W -q docs ...; echo exit=$?``.  Exit 0 with
   "skipping" info lines is ACCEPTABLE — the lines are info
   severity, NOT warnings.  -W only escalates lines emitted at
   warning level by docutils / sphinx, not arbitrary stdout.

6. **Closeout table.**  In §5 of the closeout memo, list every
   in-scope label with: LABEL / resolution path (1/2/3) / source
   file or new addition.  Add an OUT-OF-SCOPE table immediately
   below listing every deferred label with its rationale (module
   ownership, pre-existing-not-introduced-by-this-PR).

**Anti-patterns to avoid:**

- DO NOT use ``:noindex:`` / ``:nowarn:`` to suppress the info
  lines.  These are NOT warnings; suppression is unnecessary and
  obscures real Sphinx state.
- DO NOT add an unlabelled ``.. math::`` block "to silence the
  warning" without surrounding context.  The Sphinx-as-brain rule
  requires every new equation to be rich with motivation,
  citation, and code pointer — minimum effort here yields
  maximum future confusion.
- DO NOT promote a typo'd label by adding both old and new — pick
  one and update the decorator if needed (Path 2).

**Pre-edit baseline gate.**  Always grep the pre-edit Sphinx output
for the baseline label count.  An OUT-OF-SCOPE label that appears
the same number of times pre- and post-edit is unchanged and
defensible in the closeout; one that suddenly grew is a regression
to investigate.

Worked example: Issue #196 PR-CLEANUP-DOCS resolved 9 SN labels
via this pattern.  6 used Path 3 (new eq block), 2 used Path 1
(label-on-existing) via the section-vs-eq disambiguation, 1
needed both (Phase-C section anchors + new labelled eq blocks
inside them).  9 atalay / nm-1980 labels deferred to module:cp
PRs.
