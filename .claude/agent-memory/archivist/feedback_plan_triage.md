---
name: Plan-file triage decision framework
description: Five-action decision rubric (DELETE/UPDATE/POST/CONSOLIDATE/SUPERSEDED) for triaging cluster-owned plan files; archive-vs-delete heuristic and content-preservation discipline
type: feedback
---

When a parent agent dispatches a plan-cluster triage task, apply this rubric to each plan in the cluster:

**Decision rubric (pick exactly one):**

- **A — DELETE/ARCHIVE.** Plan content is fully captured in (a) issue comments, (b) Sphinx, (c) shipped code. Move to `.claude/plans/archive/<name>.md`, do NOT `rm` outright; archive preserves history. Rule of thumb: every plan worth tracking goes to archive, not the bit bucket. **Exception**: when post-cleanup close-out comments are *exhaustive* (e.g., post-#138 close-outs are 60–200 lines of falsification record per issue) and a plan is purely a pre-falsification roadmap or post-session snapshot whose *every* claim is now in the close-out comment, hard-delete is acceptable — archiving 9 plans for a single cluster bloats the archive. Reserve archive for the *master research log* (cited by name in existing comments; archiving requires patching those references via `sed -i 's|.claude/plans/X|.claude/plans/archive/X|g'` + `gh api PATCH .../comments/<id>`).
- **B — UPDATE EXISTING COMMENT** via `gh api -X PATCH repos/.../issues/comments/<id>` when augmenting a recent close-out comment with material from the plan that the close-out missed.
- **C — POST NEW COMMENT** when an OPEN issue's body reflects pre-decision state and a status update is needed (e.g., V-of-V carve-out from a sibling issue).
- **D — CONSOLIDATE** when 2+ plans describe the same forward-looking work; merge to `.claude/plans/consolidated/<topic>.md`.
- **E — SUPERSEDED.** Replaced by a later plan; record which.

**Why:** plans accumulate as scratchpads; without a triage discipline they become a corpus archaeology problem. Each plan must end either as deletion-with-content-preserved or as an active forward-looking artifact. "Just leave it" is never a valid action.

**How to apply:**

1. Read plan top-to-bottom. Note its scope, what it claims to deliver, and what's been landed (look for §"Landed deliverables" or §"Session trail" sections).
2. Identify all related issues by content + git log of cited commits + label match.
3. `gh issue view <N> --comments` to read current state of each related issue.
4. Cross-check Sphinx for content preservation: `grep -n` for plan section titles + key equation labels in `docs/theory/*.rst`.
5. **Verify drift before action**: a plan can claim "Stage 2 landed" while subsequent decisions (Cardinal Rule 2 invocations, V-of-V carve-outs) make the next stages obsolete. Read the recent CLOSED issues' close-out comments for these decisions.
6. Pick action; if B/C/D, draft in `.claude/scratch/plans-triage/<plan-name>.md` first.
7. Apply via `gh issue comment` or `gh api PATCH`.
8. `mv` plan to `.claude/plans/archive/` after action lands.
9. Final report = markdown table | Plan | Decision | Target | Action summary |.

**Anti-pattern to avoid:** writing a comment that summarises the plan's *aspirations* rather than the actual decisions/landed-state. The comment should read as a status update, not a plan summary. If the plan is fully shipped, the comment goes on the close-out issue stating "plan archived; X/Y/Z landed via commits A/B/C". If the plan is partially superseded by a Cardinal Rule decision, the comment goes on the OPEN issue stating "plan supersession reason; what's still tracked here vs what's now out-of-scope."

**Concrete example (2026-04-30 MG+topology+#138 cluster):**

- 4 of 6 plans were straight DELETE-to-archive (content shipped + already in Sphinx + close-out comments cite commits): `issue-104-multi-group-peierls`, `peierls-138-full-collapse`, `continuous-references-simplification`, `topology-based-consolidation`, `post-cp-topology-and-coordinate-transforms`.
- 1 plan triggered POST NEW COMMENT on TWO open issues (#115 + #116) because a Cardinal Rule 2 V-of-V carve-out from a sibling-issue (#138) close-out invalidated the plan's slab-retirement scope: `post-peierls-unification-close-out`. Both comments stated supersession reason + what's still tracked + plan archive location.
- Each comment ended with a "Plan housekeeping" section pointing readers to `.claude/plans/archive/<name>.md` so future searches don't dead-end.

**Concrete example 2 (2026-04-30 rank-N closure cluster, 9 plans):**

- 7 plans hard-deleted as SUPERSEDED — their entire content was already covered by post-#138 close-out comments on #119/#122/#132 (each comment is 60–206 lines of falsification record with verbatim cross-σ_t·R tables, Probe G localisation, V-S falsification, etc.). Hard delete (not archive) was correct because (a) the close-outs are exhaustive, (b) archiving 7 redundant plans would bloat the archive, (c) every cited equation label and commit is preserved in Sphinx and the comments. Plans deleted: `post-issue-119-rank-n-marshak-closure`, `post-marshak-primitives-landed-closure-open`, `post-rank-n-hollow-sphere-closure-characterization`, `issue-100-103-rank-n-class-b-multi-region`, `post-four-reference-synthesis-close-out`, `next-session-rank-n-hebert-and-beyond`, `next-session-post-retraction`.
- 1 plan archived (master research log `rank-n-closure-research-log.md`) because it is *cited by name* in #122 and #123 comment bodies; archiving required `sed -i 's|.claude/plans/X|.claude/plans/archive/X|g'` on saved comment bodies + `gh api -X PATCH .../comments/<id> -f body=@/tmp/...` to update both references.
- 1 plan distilled into a fresh POST NEW COMMENT on OPEN #112 (`post-rank-N-stepanek-calibration` → 240-line Stepanek-anchored 8-stage roadmap). The pre-existing #112 investigation-log comment carried the V1/V2/V4/V6 *what fails* data; the new comment carries the *what to try next* roadmap that the close-out lacked. Plan deleted after the comment landed.

The asymmetry between the two examples: when close-outs are skeletal, plans archive (cluster 1); when close-outs are exhaustive falsification records, plans hard-delete (cluster 2). The master research log is the universal exception — always preserve, since lessons logs accumulate generation-cost that close-outs only synthesise.
