---
name: plan-cluster-triage
description: 'Use when post-milestone cleanup requires triaging a batch of `.claude/plans/*.md` files — usually after an umbrella issue close-out lands and several sibling plans are now redundant, superseded, or partially shipped. Provides the 5-action rubric (DELETE/UPDATE/POST/CONSOLIDATE/SUPERSEDED), archive-vs-rm discipline, the master-research-log exception, the 4 sub-patterns for heterogeneous batches, and the gh-api-PATCH-existing-comment recipe. Examples: "9 plans tied to #119/#122/#132 — what to do?", "post-#138 cluster cleanup", "this plan was retracted in 2026-04-22, can I delete?".'
---

# plan-cluster-triage — post-milestone cleanup of `.claude/plans/*.md`

Plans accumulate as scratchpads. Without a triage discipline they
become a corpus-archaeology problem: a future session must distinguish
"shipped and recorded elsewhere" from "still load-bearing" by reading
each one top-to-bottom. This skill prescribes the rubric. Every plan
must end as either deletion-with-content-preserved or as an active
forward-looking artifact. **"Just leave it" is never a valid action.**

## When to use

- An umbrella issue (or a small cluster of related issues) just
  received a close-out comment and several plans in `.claude/plans/`
  refer to the same work.
- A milestone PR landed and the pre-implementation plan is no longer
  driving anything.
- A plan is explicitly banner-marked RETRACTED, OBSOLETE, or
  SUPERSEDED at line 1.
- A multi-stage phase-tracking plan has some stages shipped, some
  superseded, and the parent issue does not yet reflect the
  end-state.

## When NOT to use

- The plan is the active forward-looking artifact for an OPEN
  investigation — leave it, sharpen it instead.
- The plan is a master research log cited by name in shipped issue
  comments — always archive, never delete (see exception below).
- You haven't read the umbrella issue's close-out comment yet — the
  triage decision pivots on what the close-out already captures.

## The 5-action rubric

Pick exactly one action per plan.

| Action            | When                                                                                                     | Mechanics                                                                                                      |
| ----------------- | -------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| **A — DELETE/ARCHIVE** | Plan content fully captured in (a) issue comments, (b) Sphinx, (c) shipped code.                    | Default: `git mv` to `.claude/plans/archive/<name>.md`. Hard `rm` only when the umbrella close-out is exhaustive AND no comment cites the plan by name. See archive-vs-rm below. |
| **B — UPDATE EXISTING COMMENT** | Augmenting a recent close-out comment with material the close-out missed (new evidence, new test specs, derivation appendix). | `gh api -X PATCH repos/<owner>/<repo>/issues/comments/<id> -f body=@/tmp/new_body.md`. Add a `Updated YYYY-MM-DD (source)` header so chronology stays legible. |
| **C — POST NEW COMMENT** | An OPEN issue's body reflects pre-decision state and a status update is needed (e.g., V-of-V carve-out from a sibling issue invalidated this plan's scope). | `gh issue comment <N> --body-file <draft>`. Frame as **status update**, not plan summary (see anti-pattern below). |
| **D — CONSOLIDATE** | Two or more plans describe the same forward-looking work.                                                  | Merge into `.claude/plans/consolidated/<topic>.md`. Cite the original filenames in a "Provenance" footer.       |
| **E — SUPERSEDED** | Plan replaced by a later plan or by a shipped Cardinal Rule decision.                                       | Record which one. Move to `.claude/plans/archive/` with a 1-line `SUPERSEDED by …` banner at the top.           |

## Defaults

### Archive-not-rm by default

Move plans to `.claude/plans/archive/<name>.md` rather than `rm`. The
git history alone is not enough — a future agent searching for "what
was tried before" needs the file to exist on disk to read it. The
archive is cheap; the search-time cost of a missing record is high.

### Hard-delete is only correct when:

1. The umbrella close-out comment is **exhaustive** (60–200 lines of
   falsification record per issue with verbatim tables, probe
   localisation, V&V cross-references), AND
2. The plan is **not cited by name** in any comment (grep for the
   plan filename across `gh issue view` output for related issues),
   AND
3. Archiving would **bloat the archive** with redundant copies (e.g.
   9 plans for one cluster all saying the same thing).

When all three hold, hard-delete is acceptable. When even one fails,
archive.

### Master research-log plans always archive

Lessons logs accumulate generation-cost that close-outs only
synthesise. Even when the close-out is exhaustive, the master log is
the deepest-context source — keep it. Examples:

- `rank-n-closure-research-log.md`
- `peierls-investigation-log.md`

If the master log is **cited by name** in shipped comment bodies,
archiving requires patching those references. Workflow:

```bash
# Find citations
gh issue view <N> --json comments --jq '.comments[] | select(.body | contains(".claude/plans/X.md"))'

# Update saved comment bodies
sed -i 's|.claude/plans/X.md|.claude/plans/archive/X.md|g' /tmp/comment_*.md

# Patch via gh API
gh api -X PATCH repos/<owner>/<repo>/issues/comments/<id> -f body=@/tmp/comment_<id>.md
```

## Status-update vs plan-summary anti-pattern

**A new comment must read as a status update, not a plan summary.**
The plan's *aspirations* are already in the plan; the comment must
tell a future reader what *actually decided*, what *actually shipped*,
and what is *still in scope*.

| Anti-pattern (plan summary)                                         | Pattern (status update)                                                                              |
| ------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| "Plan #X proposes a 7-stage rollout: Stage 1 …, Stage 2 …, …"       | "Plan archived; Stages 1–4 landed via commits A/B/C/D. Stage 5 superseded by Cardinal Rule 2 V-of-V carve-out from #138 close-out — see [link]." |
| "The plan recommends Gauss-Laguerre quadrature for the polar axis." | "Gauss-Laguerre rejected post-investigation — failure mode documented in §22.9 audit comment, row 3. Production path is composite GL on `[0, R]`." |

When the plan is fully shipped, the status comment goes on the
close-out issue stating "plan archived; X/Y/Z landed via commits
A/B/C". When the plan is partially superseded by a Cardinal Rule
decision, the comment goes on the OPEN issue stating "plan
supersession reason; what's still tracked here vs what's now
out-of-scope."

## The 4 sub-patterns for heterogeneous batches

For mixed-topic plan batches (each plan ties to a different open or
closed issue), the action selection collapses to four common cases:

### 1. Explicit RETRACTED/OBSOLETE banner at line 1 → DELETE silent

The retraction notice is the durable record. Don't post a comment;
the banner already says everything a future reader needs.

> Example: `verification-spec-split-adaptive.md` had `RETRACTED 2026-04-22` at line 1 → hard-delete.

### 2. Bug investigation plan whose findings are now in `tests/l0_error_catalog.md` → DELETE

The catalog entries (with ERR-NNN, failure mode, "How it hid", "L1
test that catches", "Fix" with issue ref) ARE the post-mortem.
Cross-reference plan→ERR by symptom keywords (e.g., "panel boundary",
"ray crossing", "cross-panel").

### 3. Milestone plan whose deliverables fully shipped to code+tests → DELETE silent

Don't leave a "shipped" comment when there's nothing the post would
tell a future reader that the code+tests don't already.

> Verify: `grep -rn <plan-key-concept> orpheus/ tests/` confirms the
> planned functions and test classes both exist. If yes → DELETE.

> Example: `post-vacuum-bc-machine-precision-milestone.md` —
> verified by `grep` showing the planned functions and test classes
> exist → silent DELETE.

### 4. Phase-tracking plan with multi-stage status (some shipped, some superseded, some redefined) → POST close-out comment

Use the 9-step close-out narrative arc (see the archivist agent's
own playbook section in `AGENT.md`). Critical: include the
**load-bearing rationale** (e.g., a Gauss-Laguerre rejection, a
V-of-V carve-out reason). This is the intellectual content most
likely to be re-derived if not preserved.

## The "PATCH existing comment" sub-pattern

When a post-cleanup commit ALREADY relocated a section to GitHub but
the plan adds substantial NEW evidence (empirical scans, derivations,
redefined test specs), use `gh api -X PATCH` on the existing comment
ID rather than posting a sibling.

```bash
# Fetch existing body
gh api repos/<owner>/<repo>/issues/comments/<id> --jq .body > /tmp/existing.md

# Edit /tmp/existing.md — append new section under "Updated YYYY-MM-DD (source)" header

# Patch
gh api -X PATCH repos/<owner>/<repo>/issues/comments/<id> -f body=@/tmp/existing.md
```

The "Updated YYYY-MM-DD (source)" header preserves chronology and
tells the reader exactly when the new evidence landed and where it
came from.

> Example: #129 comment 4348745701 grew from 5.8 KB to 15.6 KB by
> appending Stage 1 measurement tables + clean-limit candidate
> analysis A-D + Test 4.1/4.2 specs.

## The "issue OPEN despite Closes-#NNN comment" pattern

When an issue's last comment says "Closes #NNN" but the issue is
still OPEN, the maintainer either (a) chose to keep it open as a
tracker for follow-ups, or (b) forgot. Don't presume — leave a
recommendation in the new comment ("Closure decision depends on
whether 'Phase X' means the machinery vs the universal claim") and
let the maintainer decide.

**Do NOT close the issue from the agent.** Closure is an explicit
maintainer action; agents recommend, maintainers act.

## Citation hygiene for relocated content

When a plan's content is archived to a GitHub comment, include the
plan filename + size in an HTML comment header:

```html
<!-- Source: .claude/plans/X.md (12.4 KB, retired with this comment) -->
```

This lets future agents grep for "retired with this comment" to find
all plan-archive comments without indexing the deletions in git.

## Cross-check before deletion

Before any DELETE action:

```bash
# Has the work shipped?
grep -rn <plan-key-concept> orpheus/ docs/ tests/

# Is the plan cited by name in any open or closed issue comment?
gh issue list --state all --search "in:comments .claude/plans/<filename>"
```

If nothing implements the concept and no issue exists, **file the
issue first** (Cardinal Rule 3: never let an improvement opportunity
pass undocumented), then delete the plan with the new issue as its
durable home.

> Example: `precision-floor-tool.md` (T3.4 / G-P2.1) had no
> implementation in code, no issue had been filed. New issue #145
> created with the plan content as the body, plan deleted, issue
> body cites the plan-file path in a `Provenance:` footer for
> lineage tracing.

## Workflow

1. Read each plan top-to-bottom. Note its scope, what it claims to
   deliver, and what's been landed (look for §"Landed deliverables"
   or §"Session trail" sections).
2. Identify all related issues by content + git log of cited commits
   + label match.
3. `gh issue view <N> --comments` to read current state of each
   related issue.
4. Cross-check Sphinx for content preservation: `grep -n` for plan
   section titles + key equation labels in `docs/theory/*.rst`.
5. **Verify drift before action**: a plan can claim "Stage 2 landed"
   while subsequent decisions (Cardinal Rule 2 invocations, V-of-V
   carve-outs) made the next stages obsolete. Read the recent
   CLOSED issues' close-out comments for these decisions.
6. Pick action; if B/C/D, draft in
   `.claude/scratch/plans-triage/<plan-name>.md` first.
7. Apply via `gh issue comment` or `gh api PATCH`.
8. `git mv` plan to `.claude/plans/archive/` (or `rm` per defaults
   above) after action lands.
9. Final report = markdown table | Plan | Decision | Target | Action summary |.

## Worked example — post-#138 quadrature cluster (2026-04-30)

The umbrella issue #133 received a §22.9 audit comment as part of the
post-#138 cleanup. The audit had per-primitive landing rows mapping
each acceptance criterion to a shipped commit. Six plans were tied to
the cluster: 4 quadrature, 2 specular.

### Triage walk-through

For each plan: read the §22.9 audit table top-to-bottom. Does every
acceptance criterion in the plan have a `**DONE — commit X**` row?

| Plan                                          | Audit-table coverage           | Action  | Reason                                                                              |
| --------------------------------------------- | ------------------------------ | ------- | ----------------------------------------------------------------------------------- |
| `quadrature-architecture.md`                  | All criteria → DONE rows       | DELETE  | Audit captures everything                                                           |
| `surface-centred-quadrature.md`               | All criteria → DONE rows       | DELETE  | Audit captures everything                                                           |
| `visibility-cone-substitution-rollout.md`     | Plan even has its own `Status: COMPLETE` header | DELETE | Plan self-confirms                                                |
| `slab-into-curvilinear.md`                    | All criteria → DONE rows       | DELETE  | Audit captures everything                                                           |
| `specular-bc-method-of-images.md`             | Captured in #133 retreat comment `Round 1 Front A` | DELETE | Image-series approach for sphere documented as never-attempted in close-out         |
| `specular-bc-phase4-multibounce-rollout.md`   | Captured in #133 retreat comment | DELETE | Mode-space `R_specular = (1/2)M⁻¹` shipped in §peierls-specular-bc                  |
| `precision-floor-tool.md` (T3.4 / G-P2.1)     | NO match — no implementation, no issue | NEW ISSUE then DELETE | Cardinal Rule 3: file new issue (#145) with plan body, delete plan after |

### Outcome

Six DELETEs (audit-captured), one new-issue-then-DELETE (precision-floor).
This skewed result (zero archives, zero comments) is correct because
the §22.9 audit was unusually thorough — every acceptance criterion
was already row-mapped to a commit. Compare to the post-#138 peierls
4-plan cluster, which had several UPDATE EXISTING COMMENT actions
because its close-outs were less exhaustive.

### Asymmetry rule

When close-outs are skeletal → plans archive (action A).
When close-outs are exhaustive falsification records → plans
hard-delete (action A with `rm`).
The master research log is the universal exception — always
preserve, since lessons logs accumulate generation-cost that
close-outs only synthesise.

## Output checklist

After triage:

- [ ] Each plan has exactly one action assigned.
- [ ] DELETE actions verified against shipped code + close-out audit.
- [ ] POST/UPDATE comments drafted as status updates, not plan
      summaries.
- [ ] Master research logs archived, not deleted.
- [ ] Plan-name citations in shipped comments patched if archive paths
      changed.
- [ ] No issue closed from agent — recommendations only.
- [ ] HTML `<!-- Source: -->` markers in archive-comments.
- [ ] Final markdown table produced for the parent agent.
