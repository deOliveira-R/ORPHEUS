# Brief — topic-file blast-radius audit, owner: explorer

**Workflow:** the harness campaign's close-out (`.claude/plans/harness_context_budget.md`, "Memory distillation — CLOSE-OUT (2026-09-21)", owed item (1)). **Phase:** AUDIT — read and propose. You are the OWNER of `.claude/agent-memory/explorer/`; you propose, the orchestrator evaluates and applies, one commit per owner. **This dispatch is READ-ONLY on tracked files:** you write exactly one file, `scratch/_blast_radius/explorer/proposal.md` (untracked), and edit nothing else.

## Artefacts from the previous phase

- Your distillation table, whose "topic files judged archaeology" section named these candidates BY NAME without reading them: `.claude/plans/archive/memory_distillation_2026-09-21/explorer/table.md`.
- The candidates (5 files):
- `.claude/agent-memory/explorer/flux_torsor_vs_cone_inventory.md`
- `.claude/agent-memory/explorer/identity_step_census_durable.md`
- `.claude/agent-memory/explorer/problem_solution_split_census.md`
- `.claude/agent-memory/explorer/sn_solve_exit_and_reflective_default.md`
- `.claude/agent-memory/explorer/phase5_mu_resolved_primitive_inventory.md`
- Their inbound references, measured 2026-09-22 by the orchestrator: `scratch/_blast_radius/explorer/referrers.md` (a claim; re-run any line you build on).
- Your live digest and index: `.claude/agent-memory/explorer/lessons.md`, `.claude/agent-memory/explorer/MEMORY.md` (and `lessons_archive.md` / `_archive/` where you have one).

## Ask

Read every candidate file in full. For each, give ONE verdict with its reason:

- **RETIRE** — campaign narrative whose every behavioural lesson is already in your digest, your archive, a rule or a skill (name where, by section), so the file brings nothing forward; and its blast radius is clean once the referrer edits below are made.
- **SALVAGE, then RETIRE** — the file holds a reusable recipe, convention or correction that is nowhere else: quote it (the lines to carry, verbatim or distilled to failure → correction), name its destination (a digest section, a named topic file, or an archive section), then retire the file.
- **KEEP** — a live pointer, a standing convention, or a load-bearing referent (a `repo` referrer that consumes it); say which.

For every RETIRE and SALVAGE, list every referrer line from `referrers.md` with the exact edit that keeps the referrer true: delete the line, re-point it to the destination, or rewrite it as history. A referrer that is itself a candidate is edited only if it survives. Your own `MEMORY.md` index lines for retired files are deleted.

The standard's default posture (the main memory's `feedback_memory_distillation_standard.md`, which your distillation brief pasted whole): a cold topic file costs nothing per dispatch, so the win of deletion is corpus hygiene, not context; **delete only when the blast radius is clean AND the content brings nothing forward.** A candidate you cannot decide without git archaeology is a RETIRE: that is the definition of archaeology.

Two lessons every brief of this campaign carries: write `proposal.md` as soon as its first verdicts are ready and append as you go (a dispatch can die mid-flight); and verify every status line you cite ("#NNN closed", "landed at <hash>") against git or GitHub (`gh issue view N --json state`, `git merge-base --is-ancestor <hash> HEAD`), never against a memory's claim. "No recent commits touch this" is not evidence.

## Rules that apply to you

*Every datum this brief states (a hash's date, a count, an exemplar's behaviour, a `Class.attr (file:line)`) is a claim; verify it by one command before building on it, and report the discrepancy as a finding.*

- `articulation`: a report is complete sentences; `[M]` measured with its command, `[R]` reasoned, `[HYPOTHESIS]` proposed; a bare number is read as measured; `[M]` on an inherited claim certifies that some measurement answered some question, so a `[M]` on a negative (absent, discarded, no consumers) is re-measured against the question at hand before it is built on.
- `code-search`: `grep` is ugrep, and an anchor inside an alternation group matches nothing, silently: use `\b…\b` or `-P` with a lookbehind; a sub-agent has no `ToolSearch`, so if Nexus is missing say so in NEEDS: and fall back to Bash.
- `coding-standards`: tests run as `python -O -m pytest`; a bare `assert` outside a collected test module is stripped under `-O`, so a contract is a `raise` and a test-side check is `np.testing.assert_*`.
- `instrument-doctrine`: a zero from a filter is evidence only after a positive control of each shape it must find (X1); every count states its predicate, its tree and its exclusions, a universal is `k of N`, and a completeness claim is re-run in Python (`re` + `pathlib.rglob`) so its denominator is stated (X2); a membership or consumer question is answered by an AST pass, a line grep reported only as the prose column, since a public name's docstring fame reads as consumption.
- `plan-authoring`: a cited precedent's adjectives are verified by reading it, at the layer (data or binder) that preserves arity; a proposed name is grepped in the prose corpus with the verb it would own; a step that adds a gate lands with the case it catches, an `⟺` checked in both directions and every symbol on a law's RHS checked against the datum's methods.
- `process-discipline`: never `git checkout`, `git restore` or `git stash` a path that carries uncommitted edits: they revert to HEAD and destroy the work (L28); revert a mutation by monkeypatching in-process or by mutating a copy; a brief's framing is a claim, so take the count it presumes (branches, fibre, consumers by AST, fused attributes) before arguing its scope and report the smaller answer plainly; a rejected candidate carries its structural reason AND the question it was refuted for.

- Read-only on tracked files; `proposal.md` is your one deliverable.
- Verify each referrer line by re-running its grep (`grep -n "<stem>" <file>`); a discrepancy is a finding.
- Nexus: agent-memory files are graph nodes. If `mcp__nexus__*` tools are available, `context("<node id>")` on one candidate tells you its inbound edges; if they are not available to you, say so in `NEEDS:` and rely on the grep census.
- Never `git checkout`, `git restore` or `git stash` anything.

## Return contract

Write `scratch/_blast_radius/explorer/proposal.md`: a table `| file | verdict | reason (one line) | referrer edits |`, then a "Salvaged" section with the quoted lines and their destinations, then a `NEEDS:` block for anything you could not obtain. Reply to the orchestrator in at most 120 words: the counts per verdict, anything surprising, and the path. No prose beyond that; the file carries the detail.
