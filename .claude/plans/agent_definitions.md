# Agent definitions a context-free dispatch can trust — current, single-sourced, and enforced where the harness can enforce

**Status: LIVING PLAN (`plan-authoring` §0), opened 2026-09-22 at the first exchange. No implementation until the user rules the plan polished.** The seven issues waiting at `harness_context_budget.md` ⏸ COMPACTION POINT #13 (#488, #484, #491, #492, #494, #495, #496) wait until this is done.

**The instruction (the user, 2026-09-22, verbatim):** *"Before we tackle those I want to improve agent definitions."* and, after compaction, *"We want to review and improve agent definitions."* The review below is the first pass; the user has not yet said what they themselves want changed, so their own list is the first thing this plan asks for.

## 1. Goal, in the harness's own terms

A dispatched agent starts with no session context. Its definition (front matter, role block, body, preloaded skills, its own memory) plus the brief are all it knows. The goal: every definition tells its agent **exactly what it may do, what its role is, and the method no other page carries** — true against today's tree, with no statement that contradicts another page, and with every capability claim enforced by the harness where the harness can enforce it.

Done-when (a hypothesis until the shape is ruled): `python -m tools.harness --check` reads every line of every agent body (today it reads 0 of the hand bodies' lines), and a dispatch of each agent with a fixed probe brief reports tools, write scope and memory policy that match its definition (measured by dispatching, the harness page's rule).

## 2. What an agent definition is (the ontology question)

An agent definition is a **projection of the project's knowledge onto one role, plus a capability grant**. Its parts, and who owns each today:

| part | what it is | owner today | owner it should have `[R]` |
|---|---|---|---|
| capability | model, tools, memory scope, `omitClaudeMd`, and write scope | AGENT.md front matter (hand) + brief prose | front matter, and hooks for what front matter cannot express |
| role | phases, whom it calls, return contract | `docs/development/agents/<name>.md` → generated role block | unchanged |
| method | the procedure only this role runs (the qa L0 interrogatives, the investigator's cascade, the attacker's trigger procedure) | AGENT.md body (hand, unchecked) | the agent's source page, generated and checked |
| knowledge | domain snapshots (SN operator shape, project layout) and restated doctrine (V&V rules, the build gate) | AGENT.md body (hand, unchecked) | nowhere in the body: a pointer to the owning page, or a generated excerpt of it |
| memory | what the agent learned | `.claude/agent-memory/<name>/`, injected by the harness | unchanged; its write policy is one stated rule |

Every defect found below is either a part held by the wrong owner (knowledge copied into a body, where it drifts) or a capability stated in prose that the harness does not grant or does not enforce.

## 3. The review, 2026-09-22 (first pass, by the orchestrator)

Instruments: each of the nine AGENT.md files and role sources read in full; paths and symbols checked by `scratch/_agent_defs/agent_refs.py` (an AST pass over `orpheus/` for definitions, with `IterationRecord` as the positive control); defect patterns counted by `scratch/_agent_defs/agent_census2.py` (a regex census, positive control qa's known skill-edit line; its counts are read as leads and every instance cited below was read in its body); the Claude Code sub-agent documentation (`code.claude.com/docs/en/sub-agents`, fetched 2026-09-22); description and pasted-block sizes by `scratch/_agent_defs/agent_desc.py` (each probe runs as `ROOT=$PWD .venv/bin/python <probe>`). "Body" means the hand-maintained text after the generated role block.

### D1. The harness checks none of the hand-written bodies

`[M]` `tools/harness/targets/claude_code.py`: an agent's generated target is its role block only. `tools/harness/ids.py` `pages_under_check`: the ID and link checks read `docs/development/agents/*.md` (9 to 11 lines each), not the bodies (125 to 654 lines each). So a dead path, a retired section number or a stale symbol in a body is invisible to `--check` and to CI. D3 is what that invisibility produced.

### D2. Capability and mandate disagree

- **The cross-domain-attacker has no text search at all.** Its tools are Read, Grep, Glob, WebSearch, WebFetch; the standalone Grep and Glob tools no longer exist (the `code-search` rule, case 2026-06-14; `[R]` not re-measured in a sub-agent), and it has no Bash. Today's audit brief asked it to verify with `gh`/`git` and it could not (CP#14 observation 1).
- **Grep and Glob are listed in 8 of 9 allowlists** (all but the explorer) `[M]`: dead configuration that tells a reader of the front matter the agent can search.
- **Write and Edit arrive by side effect.** `[M]` the documentation: `memory:` "automatically enables Read, Write, and Edit" so the agent can manage its memory. This is why the roster shows Write and Edit for the explorer, the attacker and the literature-researcher. The literature-researcher's deliverable (a memo written incrementally) works only through this side effect. The explorer's body says "You NEVER modify files", yet the brief template asks it for a listing file and it runs probes from files.
- **The method-implementer is told to resume a sibling "by name with SendMessage"**, and 9 of 9 allowlists lack `SendMessage` `[M]`. The documentation: a sub-agent holds `SendMessage` only when its `tools:` list names it.
- **Edit scope and agent memory collide** (CP#14 observation 2). 6 of 9 bodies mandate a memory update after every task (archivist, attacker, elegance-enforcer, literature-researcher, numerics-investigator, qa; the census), while briefs say "read-only on tracked files" and agent memory is tracked. Four agents resolved it four ways today.
- **"Read-only" was read as "no mutation"** (CP#14 observation 4). qa's body item 11 already requires in-process mutation and forbids only edits to tracked production files; the brief's "read-only" overrode it, and qa proposed six `verifies` witnesses unmutated, three of which stayed green under the defect.

### D3. Stale knowledge in bodies (each `[M]` against today's tree)

- archivist: its project tree cites `docs/theory/discrete_ordinates.rst`, `docs/theory/collision_probability.rst` and `orpheus/derivations/{sn_contamination,sn_heterogeneous,cp_sphere}.py`, none of which exists (the theory tree was restructured into parts on 2026-07-15).
- test-architect: `derivations._xs_library.get_mixture` lives at `orpheus/derivations/common/xs_library.py`; `homogeneous_1d` has no definition anywhere in the repository.
- elegance-enforcer: `SNMesh` and `from_setup` have no definition in `orpheus/`.
- explorer: `TraceSpace` and `BoundaryResidual` have no class of those names (the tree has `AngularTraceSpace`, `ScalarTraceSpace`, `AngularBoundaryResidual`); its project layout omits `transport/` (the L2 layer) and `plotting.py`, and calls `numerics/` "Shared numerics (eigenvalue protocol)" where CLAUDE.md's layer table says "mathematics only: spaces, measures, quadrature, operators". The explorer runs with `omitClaudeMd`, so this stale map is the only map it has.
- qa, test-architect, archivist: cite `error_catalog.md`; the file is `docs/theory/verification/error_catalog.rst`.
- numerics-investigator, test-architect: cite `vv-principles` by numbered sections (`§1`, `§4`, `§6`, `§H2`); the skill's 11 `##` headings are named, not numbered.
- elegance-enforcer: "Read `CLAUDE.md` Cardinal Rules 1 and 2" (they live in `.claude/rules/cardinal.md`, already loaded for this agent); three "(See memory note: …)" pointers name main-agent memory notes the agent cannot read, whose content now lives in the `coding-standards` and `process-discipline` rules.
- elegance-enforcer: its body carries a pasted copy of the harness's memory instructions, 13 072 characters ≈ 3 268 tokens, ending "Your MEMORY.md is currently empty"; its MEMORY.md has 145 lines. The harness injects the real block at every dispatch (the documentation), so the agent reads the instructions twice and one copy is false.
- archivist: "Follow ALL four directives", followed by five.

### D4. One concept, two definitions (X4)

- **A VIOLATION's "three legs"**: the workflows page (W1-P3) says what, which pattern, the remedy; the elegance-enforcer body says bug-habitat argument, coextensiveness check, live-tree verification.
- **The docs build gate**: the archivist body says the gate is the warning SET diff, never a count (Quality checklist item 2; the build-gate section), and elsewhere "Verify the count is unchanged pre/post-edit, not the content" (the close-out quality gates).
- **Restated doctrine**: qa item 11 and the numerics-investigator's Rules 3, 4, 7 restate `vv-principles`; the test-architect's "Cardinal Rule" section restates the 1-group degeneracy under a name the cardinal page does not use (it has five rules). Each restatement is a copy that drifts when the skill moves.

### D5. Self-improvement directives aim at generated files and grow rule restatements

4 of 9 bodies direct the agent to edit `vv-principles`' SKILL.md (qa: "add it to the skill BEFORE completing the review"; numerics-investigator: "update `vv-principles` SKILL.md §Anti-patterns in the same commit"; test-architect: "append the row to the skill's table BEFORE delivering"; archivist Directive 5), and the method-implementer proposes "an edit to the relevant skill". `.claude/skills/vv-principles/SKILL.md` is GENERATED: a direct edit is drift that `--check` and CI reject. None of the five tells the agent to check the rules before proposing, which is CP#14 observation 6 (qa and the archivist proposed lessons that restate existing clauses).

### D6. Sub-agents told to converse mid-dispatch

The elegance-enforcer "Ask for clarification if the scope is ambiguous" and "ask the main agent before proceeding"; the method-implementer "report back to the user"; the literature-researcher "Before searching, clarify". A dispatched agent's only channel back is its return, and the return contract's `NEEDS:` block is that channel.

### D7. Role identity disagrees with the role block

The method-implementer body: "this agent BUILDS new code; numerics-investigator FIXES existing code." Its role block: "W1-P2 (build), W2 (the fix)", and W2 says "the fix lands (implementer or main agent)". The body's deliverable manifest is a Branch-1 SymPy module, a Branch-2 solver and an L1 cross-check, the shape of a published-formulation reference solver, while the project's builds this year were operator-algebra carves.

### D8. Context cost carried for nothing

- Descriptions: every holder of the Agent tool pays the roster, 5 977 characters ≈ 1 494 tokens `[M]`, of which the elegance-enforcer's three `<example>` blocks are 2 615 characters ≈ 653 tokens (44 %).
- The elegance-enforcer's pasted memory block (D3), ≈ 3 268 tokens per dispatch.
- Preloads are chosen per agent with no stated criterion; `vv-principles` + `numerical-bug-signatures` + `coding-elegance` ≈ 30K tokens ride every qa, test-architect, numerics-investigator and method-implementer dispatch (CP#14 census).

### D9. Naming

The archivist's front matter spells `name: Archivist`; its memory directory is `archivist/`, and its body carries a paragraph working around the case mismatch ("the worktree-isolation guard rejects writes to the non-canonical casing"). The root fix is the name. `[HYPOTHESIS]` the refused write of CP#14 observation 3 is the same mismatch or the permission layer; `.claude/settings.json` has no hook that guards Write (its Write hook is the PostToolUse Nexus brief), so it is not a project hook.

### Also found, needing a measurement before anything relies on it

- The documentation says a non-fork sub-agent never receives the main conversation's auto memory; the harness page records, `[M]` 2026-09-21, that the project memory index is inherited by every dispatch without `omitClaudeMd`. One of the two is out of date; re-measure by dispatching before a definition leans on either.
- Harness features no agent uses (the documentation): per-agent `hooks`, `disallowedTools`, `effort`, `maxTurns`, `permissionMode`, `isolation: worktree`. `isolation: worktree` branches from the default branch, not the parent's HEAD, so it cannot review uncommitted work (refuted FOR mutating uncommitted work under review; the FACT: it isolates an agent that edits committed code).

### What worked (keep)

The two-pass review (qa withdrew attacks with reasons); probes run as scripts from `/tmp` with `module.__file__` printed; the "every datum is a claim" sentence for the `omitClaudeMd` agents; the attacker's trigger procedure and its UNEXPLORED block; the literature-researcher's local-folder-first procedure.

## 4. Candidate shapes

**S1. Patch in place.** Fix every D-item in the hand bodies. Makes nothing unspellable: D1 stands, so the next rename re-stales a body silently.

**S2. Generate the whole body.** Move each hand body into `docs/development/agents/<name>.md` (role + method), so AGENT.md is front matter + one generated block. `--check` then reads every body line: dead links, dead IDs (`§6`, `error_catalog.md`), budgets. Knowledge leaves the bodies for pointers; for the `omitClaudeMd` explorer, the layer table is generated from its one source (CLAUDE.md's on-boarding page) instead of hand-copied. Makes D1 and D3 detectable and most of D4 visible; leaves D2 (capability in prose) and D5 (the learning channel) as review questions.

**S3. S2 plus a capability contract the harness enforces** `[HYPOTHESIS]`, recommended.
- tools that match the mandate: Bash for the attacker; Write and Edit listed wherever the role writes; Grep and Glob removed; `SendMessage` where a body resumes by name;
- a PreToolUse hook on Edit and Write that refuses any file carrying the GENERATED stamp, for every agent and the main agent (makes D5's direct skill edit unspellable at write time, before CI);
- a stated write scope per agent, enforced by a hook in its front matter where it is narrower than the tree (reviewers and the explorer: `scratch/`, `/tmp` and their own memory);
- one definition of "read-only" in the brief template (no edits to tracked files; in-process mutation and scratch writes allowed; the agent's own memory handled per the memory ruling);
- one learning channel (the memory ruling below): an agent never edits a skill or a rule; it proposes, and each proposal names the clause that does not already cover it.

Makes D2, D5 and the D4 restatements unspellable or checked; leaves domain judgement (which method, which preloads) to the rulings.

**Refuted so far:** granting agents the `Skill` tool so preloads can shrink (refuted FOR context cost by the 2026-09-21 ruling: the roster costs 2 613 first-turn tokens per dispatch and grows with plugins the project never uses; the FACT: a preload is the only reproducible way a sub-agent receives a skill).

## 5. Questions for the user

1. **What did you want improved?** The review above is mine; your own list decides the scope.
2. **Target shape:** S3 (recommended), S2, or S1.
3. **Agent memory:** (a) an agent writes its own memory freely and the orchestrator reviews and commits it at close-out; (b) an agent never writes memory mid-dispatch and returns a `LESSONS:` block the orchestrator distills; (c) (a), with each lesson naming the rule clause that does not already cover it.
4. **The method-implementer's identity:** the general W1-P2/W2 builder over the operator algebra (body rewritten to its role block), or the published-formulation reference-solver builder its body describes.
5. **Preload criterion** (proposed `[R]`): a skill is preloaded when the role applies it at every dispatch; otherwise the brief names the page. Applied, it would review, for instance, `coding-elegance` in qa (the enforcer's axis) and `numerical-bug-signatures` in the test-architect.

## 6. Rulings ledger

(none yet)

## 7. Implementation starts when

The user rules this plan polished. Before then: no edit to any AGENT.md, role source, hook or setting. The order, once ruled `[HYPOTHESIS]`: the measurements of §3's last block; the body move (S2) with D3 and D4 fixed as each body moves; the capability contract; the learning channel; then a probe dispatch of each agent confirming what it received.

## Resume surface

This file; `docs/development/harness.md` ("Adding or changing", "What loads, and what it costs"); `docs/development/workflows.md` ("The brief"); `.claude/agents/*/AGENT.md`; `harness_context_budget.md` ⏸ COMPACTION POINT #14 (the census and the eight dispatch observations this review extends).
