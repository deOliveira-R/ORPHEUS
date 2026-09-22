# Brief — topic-file blast-radius audit, owner: numerics-investigator

**Workflow:** the harness campaign's close-out (`.claude/plans/harness_context_budget.md`, "Memory distillation — CLOSE-OUT (2026-09-21)", owed item (1)). **Phase:** AUDIT — read and propose. You are the OWNER of `.claude/agent-memory/numerics-investigator/`; you propose, the orchestrator evaluates and applies, one commit per owner. **This dispatch is READ-ONLY on tracked files:** you write exactly one file, `scratch/_blast_radius/numerics-investigator/proposal.md` (untracked), and edit nothing else.

## Artefacts from the previous phase

- Your distillation table, whose "topic files judged archaeology" section named these candidates BY NAME without reading them: `.claude/plans/archive/memory_distillation_2026-09-21/numerics-investigator/table.md`.
- The candidates (3 files):
- `.claude/agent-memory/numerics-investigator/phase5a_moment_consuming_scatter_derisk.md`
- `.claude/agent-memory/numerics-investigator/r1_step_d_sphere_preconditioner_oscillation.md`
- `.claude/agent-memory/numerics-investigator/sn_sig_t_layout_drift_indexerror.md`
- Their inbound references, measured 2026-09-22 by the orchestrator: `scratch/_blast_radius/numerics-investigator/referrers.md` (a claim; re-run any line you build on).
- Your live digest and index: `.claude/agent-memory/numerics-investigator/lessons.md`, `.claude/agent-memory/numerics-investigator/MEMORY.md` (and `lessons_archive.md` / `_archive/` where you have one).

## Ask

Read every candidate file in full. For each, give ONE verdict with its reason:

- **RETIRE** — campaign narrative whose every behavioural lesson is already in your digest, your archive, a rule or a skill (name where, by section), so the file brings nothing forward; and its blast radius is clean once the referrer edits below are made.
- **SALVAGE, then RETIRE** — the file holds a reusable recipe, convention or correction that is nowhere else: quote it (the lines to carry, verbatim or distilled to failure → correction), name its destination (a digest section, a named topic file, or an archive section), then retire the file.
- **KEEP** — a live pointer, a standing convention, or a load-bearing referent (a `repo` referrer that consumes it); say which.

For every RETIRE and SALVAGE, list every referrer line from `referrers.md` with the exact edit that keeps the referrer true: delete the line, re-point it to the destination, or rewrite it as history. A referrer that is itself a candidate is edited only if it survives. Your own `MEMORY.md` index lines for retired files are deleted.

The standard's default posture (the main memory's `feedback_memory_distillation_standard.md`, which your distillation brief pasted whole): a cold topic file costs nothing per dispatch, so the win of deletion is corpus hygiene, not context; **delete only when the blast radius is clean AND the content brings nothing forward.** A candidate you cannot decide without git archaeology is a RETIRE: that is the definition of archaeology.

Two lessons every brief of this campaign carries: write `proposal.md` as soon as its first verdicts are ready and append as you go (a dispatch can die mid-flight); and verify every status line you cite ("#NNN closed", "landed at <hash>") against git or GitHub (`gh issue view N --json state`, `git merge-base --is-ancestor <hash> HEAD`), never against a memory's claim. "No recent commits touch this" is not evidence.

## Rules that apply to you

The project rules load with you; the items below are the task-specific ones.

- Read-only on tracked files; `proposal.md` is your one deliverable.
- Verify each referrer line by re-running its grep (`grep -n "<stem>" <file>`); a discrepancy is a finding.
- Nexus: agent-memory files are graph nodes. If `mcp__nexus__*` tools are available, `context("<node id>")` on one candidate tells you its inbound edges; if they are not available to you, say so in `NEEDS:` and rely on the grep census.
- Never `git checkout`, `git restore` or `git stash` anything.

## Return contract

Write `scratch/_blast_radius/numerics-investigator/proposal.md`: a table `| file | verdict | reason (one line) | referrer edits |`, then a "Salvaged" section with the quoted lines and their destinations, then a `NEEDS:` block for anything you could not obtain. Reply to the orchestrator in at most 120 words: the counts per verdict, anything surprising, and the path. No prose beyond that; the file carries the detail.
