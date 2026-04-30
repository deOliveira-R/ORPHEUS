---
name: Misc-cluster plan triage
description: Patterns for triaging mixed-topic .claude/plans/ files where each plan ties to a different open or closed GitHub issue
type: feedback
---

When triaging a small heterogeneous batch of plans (each tied to its own issue), the action selection collapses to four common patterns. Document each so future triage skips re-deriving the rubric:

1. **Plan with explicit RETRACTED/OBSOLETE banner at top** → DELETE outright. The retraction notice is the durable record. (Example: `verification-spec-split-adaptive.md` had "RETRACTED 2026-04-22" at line 1.)

2. **Bug investigation plan whose findings are now in `tests/l0_error_catalog.md`** → DELETE. The catalog entries (with ERR-NNN, failure mode, "How it hid", "L1 test that catches", "Fix" with issue ref) ARE the postmortem. Cross-reference plan→ERR by symptom keywords (e.g., "panel boundary", "ray crossing", "cross-panel").

3. **Milestone plan whose deliverables fully shipped to code+tests** → DELETE without comment. (Example: `post-vacuum-bc-machine-precision-milestone.md` — verified by `grep` showing the planned functions and test classes both exist.) Don't leave a "shipped" comment when there's nothing the post would tell a future reader that the code+tests don't already.

4. **Phase-tracking plan with multi-stage status (some shipped, some superseded, some redefined)** → POST a close-out narrative comment to the parent issue. Use the 9-step structure from `feedback_closeout_docs.md`. Critical: include the **load-bearing rationale** (e.g., the Gauss-Laguerre rejection in the slab plan) — this is the intellectual content most likely to be re-derived if not preserved.

**The "appended-cleanup-comment" sub-pattern**: when the post-#138 cleanup ALREADY relocated a section to GitHub but the plan adds substantial NEW evidence (empirical scans, derivations, redefined test specs), use `gh api -X PATCH` on the existing comment ID rather than posting a sibling. Add a clear "Updated YYYY-MM-DD (source)" header line so the chronology stays legible. (Example: #129 comment 4348745701 grew from 5.8 KB to 15.6 KB by appending Stage 1 measurement tables + clean-limit candidate analysis A-D + Test 4.1/4.2 specs.)

**The "issue OPEN despite closeout-comment" pattern**: when an issue's last comment says "Closes #NNN" but the issue is still OPEN, the maintainer either (a) chose to keep it open as a tracker for follow-ups, or (b) forgot. Don't presume — leave a recommendation in the new comment ("Closure decision depends on whether 'Phase X' means the machinery vs the universal claim") and let the maintainer decide. Do NOT close the issue from the agent.

**Citation hygiene for relocated content**: when a plan archives findings to a GitHub comment, include the plan filename + size in an HTML comment header (`<!-- Source: .claude/plans/X.md (N KB, retired with this comment) -->`). This lets future agents grep for "retired with this comment" to find all plan-archive comments without indexing the deletions.
