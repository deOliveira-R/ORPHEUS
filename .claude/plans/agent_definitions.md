# Agent definitions a context-free dispatch can trust — current, single-sourced, enforced where the harness can enforce, and each read as one refined mandate

**Status: LIVING PLAN (`plan-authoring` §0), opened 2026-09-22 at the first exchange. No implementation until the user rules the plan polished.** The seven issues waiting at `harness_context_budget.md` ⏸ COMPACTION POINT #13 (#488, #484, #491, #492, #494, #495, #496) wait until this is done.

**The instruction (the user, 2026-09-22, verbatim):** *"Before we tackle those I want to improve agent definitions."* and, after compaction, *"We want to review and improve agent definitions."* The review below (§3) is the first pass, by the orchestrator.

**The user's answer (2026-09-22):** the review is endorsed as the foundation ("an excellent foundation before what I meant"); S3 is the shape; D6 is corrected (now measured, §3); and what the user meant is five further items (§5), which follow the S3 foundation. The rulings are in §7.

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

**The boundary-of-failure loop (the user, 2026-09-22, ruling R9).** Three agents are one loop over the boundaries of failure of a capability:
- the **test-architect is foresight**: it predicts every boundary of failure before the capability exists, and its ladder of tests is those boundaries, ordered so that each rests on the verified ones below it;
- the **numerics-investigator is search**: it is dispatched when (1) the boundaries were identified but the implementation failed at one, or (2) a boundary was not identified and must be found; its probe cascade is the ladder built after the fact;
- the **qa is hindsight**: with the capability and its tests both existing, it checks that the boundaries were right and that each is tested. When hindsight finds a new requirement (a boundary nobody predicted) or a vacuous test (a duplicate or invalid boundary), the loop repeats from foresight.

Each of the three definitions states its place in this loop, and the loop is the reason the ladder principle (R7) is shared doctrine rather than one agent's method.

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

`[REFUTED 2026-09-22 as stated, the user's correction]` The harness has a mid-dispatch channel. `[M]` 2026-09-22, a background `general-purpose` haiku probe holding `SendMessage`: it sent a message to `main` at 03:29:00, the orchestrator replied, and the reply reached the probe at its sixth tool call (03:29:06); it held no `AskUserQuestion` (the documentation: that tool is removed from every sub-agent). The facts that survive: (i) a sub-agent can reach the orchestrator mid-dispatch, never the user; (ii) a reply is delivered at the agent's next tool call, so an agent waiting for one must keep working; (iii) 9 of 9 of our allowlists omit `SendMessage` `[M]`, so for our agents the channel is closed by our own configuration, and the four "ask" instructions above name no mechanism the agent holds. Whether to open it, and for which agents, is question Q1 (§6).

### D7. Role identity disagrees with the role block

The method-implementer body: "this agent BUILDS new code; numerics-investigator FIXES existing code." Its role block: "W1-P2 (build), W2 (the fix)", and W2 says "the fix lands (implementer or main agent)". The body's deliverable manifest is a Branch-1 SymPy module, a Branch-2 solver and an L1 cross-check, the shape of a published-formulation reference solver, while the project's builds this year were operator-algebra carves.

### D8. Context cost carried for nothing

- Descriptions: every holder of the Agent tool pays the roster, 5 977 characters ≈ 1 494 tokens `[M]`, of which the elegance-enforcer's three `<example>` blocks are 2 615 characters ≈ 653 tokens (44 %).
- The elegance-enforcer's pasted memory block (D3), ≈ 3 268 tokens per dispatch.
- Preloads are chosen per agent with no stated criterion; `vv-principles` + `numerical-bug-signatures` + `coding-elegance` ≈ 30K tokens ride every qa, test-architect, numerics-investigator and method-implementer dispatch (CP#14 census).

### D9. Naming

The archivist's front matter spells `name: Archivist`; its memory directory is `archivist/`, and its body carries a paragraph working around the case mismatch ("the worktree-isolation guard rejects writes to the non-canonical casing"). The root fix is the name. `[HYPOTHESIS]` the refused write of CP#14 observation 3 is the same mismatch or the permission layer; `.claude/settings.json` has no hook that guards Write (its Write hook is the PostToolUse Nexus brief), so it is not a project hook.

### Also found, needing a measurement before anything relies on it

- `[M]` 2026-09-22 (the D6 probe): a background sub-agent holding `SendMessage` messages `main` and receives the reply at its next tool call.
- The documentation says a non-fork sub-agent never receives the main conversation's auto memory; the harness page records, `[M]` 2026-09-21, that the project memory index is inherited by every dispatch without `omitClaudeMd`. `[M]` 2026-09-22, §8 step 1 (one background `general-purpose` haiku probe, no `omitClaudeMd`, answering from its context before any tool call): the index IS inherited, and it is the SESSION-START snapshot (the probe quoted the index line as it stood before this session edited it); the rules load (control: the cardinal page's heading present), and the path-scoped `plan-authoring` does not until triggered. The documentation is out of date for this harness; the harness page stands.
- `[M]` 2026-09-22, the same probe: a WRITE of a new file under `.claude/plans/` did not load the path-scoped `plan-authoring` rule (no system reminder after the write), where a READ under a scoped path does (the harness page, `[M]` 2026-09-21). So a path-scoped rule reaches whoever reads a matching file, and every Edit is preceded by a Read, but creating a new page loads nothing: the archivist's definition says to read the `documentation` rule before creating a page (R11).
- Harness features no agent uses (the documentation): per-agent `hooks`, `disallowedTools`, `effort`, `maxTurns`, `permissionMode`, `isolation: worktree`. `isolation: worktree` branches from the default branch, not the parent's HEAD, so it cannot review uncommitted work (refuted FOR mutating uncommitted work under review; the FACT: it isolates an agent that edits committed code).

### What worked (keep)

The two-pass review (qa withdrew attacks with reasons); probes run as scripts from `/tmp` with `module.__file__` printed; the "every datum is a claim" sentence for the `omitClaudeMd` agents; the attacker's trigger procedure and its UNEXPLORED block; the literature-researcher's local-folder-first procedure.

## 4. Candidate shapes

**S1. Patch in place.** Fix every D-item in the hand bodies. Makes nothing unspellable: D1 stands, so the next rename re-stales a body silently.

**S2. Generate the whole body.** Move each hand body into `docs/development/agents/<name>.md` (role + method), so AGENT.md is front matter + one generated block. `--check` then reads every body line: dead links, dead IDs (`§6`, `error_catalog.md`), budgets. `[REFUTED 2026-09-22 in part, at implementation]` The ID check skips code spans by design and has no `§N` kind, so neither example would have reddened: D3's stale references are paths in code spans. The instrument that sees them is a new one, `tools/harness/paths.py` (a code-span file check over every generated page), which landed with the move; its first run found 1 dead file in the skills and 17 in the moved bodies. A numbered section reference (`§6`) and a bare identifier (`SNMesh`) stay unchecked, and the rewrite removes them. Knowledge leaves the bodies for pointers; for the `omitClaudeMd` explorer, the layer table is generated from its one source (CLAUDE.md's on-boarding page) instead of hand-copied. Makes D1 and D3 detectable and most of D4 visible; leaves D2 (capability in prose) and D5 (the learning channel) as review questions.

**S3. S2 plus a capability contract the harness enforces** `[HYPOTHESIS]`, recommended.
- tools that match the mandate: Bash for the attacker; Write and Edit listed wherever the role writes; Grep and Glob removed; `SendMessage` where a body resumes by name;
- a PreToolUse hook on Edit and Write that refuses any file carrying the GENERATED stamp, for every agent and the main agent (makes D5's direct skill edit unspellable at write time, before CI);
- a stated write scope per agent, enforced by a hook in its front matter where it is narrower than the tree (reviewers and the explorer: `scratch/`, `/tmp` and their own memory);
- one definition of "read-only" in the brief template (no edits to tracked files; in-process mutation and scratch writes allowed; the agent's own memory handled per the memory ruling);
- one learning channel (the memory ruling below): an agent never edits a skill or a rule; it proposes, and each proposal names the clause that does not already cover it.

Makes D2, D5 and the D4 restatements unspellable or checked; leaves domain judgement (which method, which preloads) to the rulings.

**Refuted so far:** granting agents the `Skill` tool so preloads can shrink (refuted FOR context cost by the 2026-09-21 ruling: the roster costs 2 613 first-turn tokens per dispatch and grows with plugins the project never uses; the FACT: a preload is the only reproducible way a sub-agent receives a skill).

## 5. The user's scope — five items that follow the S3 foundation (2026-09-22)

The user, on what they meant: *"any of these should happen only after your proposal, which as I said, is an excellent foundation."* Each item: the user's intent, what the tree holds today (measured), and the design sketch (a hypothesis until ruled).

### 5.1 Each agent's memory audited against its definition

**Intent (the user):** find where an agent's lessons and surprises exist *"because the role is not better specified"*, and solve them by improving the definition.

**Today:** nine memory directories (index and digest sizes in the CP#14 census). The 2026-09-21 distillation asked a different question of the same files: which lessons restate a RULE (those retired to the rules). This audit asks which lessons fill a gap in the agent's own DEFINITION.

**Sketch `[HYPOTHESIS]`:** every memory entry is classed as (i) restates a rule or a skill: retire it (the distillation law); (ii) a gap in the definition, something the role should have told the agent before it had to learn it: fold it into the definition and retire the entry; (iii) experience local to this agent's work: keep it. The output per agent is the list of class-(ii) entries with the definition sentence each implies; it feeds the rewrite of §5.5. Who classifies: each agent on its own memory, with a fixed output schema, reviewed by the orchestrator (the agent knows why it wrote each entry; the orchestrator holds the definition's intent).

### 5.2 The test-architect builds on the tests that exist, as a ladder by complexity

**Intent (the user):** before designing, find the tests that already exist for the capability and order them by complexity ("this tests X, this other tests X and Y"), building a hierarchy until the full capability is tested; improve an existing test before creating a new one. *"Fewer tests that are well thought out, properly organized by difficulty, and build one upon another in a clearly structured way, are much superior to adding 10 half-thought, potentially vacuous tests."* The worked example, a boundary condition like the albedo: the foundations first (vacuum works; reflective works), then the edges (albedo 0 equals vacuum; albedo 1 equals reflective), then the interior (0 < α < 1).

**Today:** the test-architect's body opens at "Identify the feature being verified" and designs a matrix from scratch; it has no survey step and no "extend before adding" rule. `[M]` 2026-09-22: no clause of `vv-testing`, `vv-principles` or `coding-standards` says to extend an existing test before adding one (one `grep -iE` over the three sources for existing-test, duplicate-test, test-hierarchy and builds-on phrasings: 0 hits there, while the same pattern hits 5 lines elsewhere in `docs/development/`, the positive control). The nearest doctrine is `vv-principles`' necessity chain (L1 without L0 is compensating errors), which orders V&V LEVELS, not the tests of one capability.

**Sketch `[HYPOTHESIS]`:** the survey is the first step of every spec: the existing tests of the capability found by Nexus (the equations' `verifies` edges, the runtime exercisers of the touched symbols) and by grep, each placed on a ladder of rungs (foundation, edge, interior, composition), each rung naming the rungs it rests on; the spec's deliverable is the ladder with its gaps, each gap filled by improving a test when one sits on that rung and by a new test only when none does. An edge rung asserts equality to a foundation rung (albedo 0 against vacuum), a structurally independent reference by construction. The principle (a capability's tests form a ladder; a test not placed on it is a finding) may belong in `vv-principles` so that qa reviews against it too; the procedure belongs in the test-architect's definition (question Q2).

**The user's refinement (2026-09-22, ruling R7):** the test-architect still delivers a test matrix, and its rows follow the order **reuse > improve > new**. The reason for the hierarchy is diagnosis: a flat suite says *that* something failed, not *where the boundary of the failure is*. The expected shape is a prediction the matrix can be checked against. A new capability built on an established foundation (a new boundary condition) should be mostly reuse, with some improvement (a test generalised, for instance) and some new tests for the capability's own behaviour and limits; a genuinely new capability must have behaviour no existing test covers, by definition. A refinement, merge or generalisation that lets one capability express several specialised cases is the exception that tests little that is new, and even it typically opens new capabilities that need new tests. `[R]` A consequence worth stating in the clause: a ladder whose rungs each rest on verified lower rungs is the numerics-investigator's probe cascade built in advance, so a red reads its own diagnosis (the lowest red rung bounds the defect).

### 5.3 The archivist writes the present; the past goes to the page's end, the future to issues and plans

**Intent (the user):** documentation reflects the current state of the code. Older things go to the page's history section and its gotchas, or to a collapsible box showing something important that was first got wrong and how it was got right. *"The main documentation body should not mix past, present and future. It should be exclusively about the present. Past goes to auxiliary sections at the end of the page and future is the scope of GitHub issues and plans files."*

**Today `[M]` 2026-09-22:**
- The procedure exists: the theory-page template of #231 (its settled design, recorded in `.claude/plans/archive/sn_doc_architecture_231.md` §"The 9-section template") puts Gotchas at section 8 and History at section 9 as ONE collapsed changelog, and relocates narrative essays to issues. `sphinx_design` is loaded (`docs/conf.py`), and 10 theory pages already carry a Development history section.
- The procedure has no home in `docs/development/`: its authority is an archived plan and an issue comment, so no agent definition can point at it.
- The archivist's body teaches the opposite. Its "Close-Out Narrative Arc", which it calls its most-used playbook, keeps the motivation in the body with its tenses flipped, puts retraction tombstones above the content they retract, keeps falsified tables in the body, and lists a session trail. Its Directive 3 rubric scores "Failed approaches: full history with rationale" as excellent.

**Sketch `[HYPOTHESIS]`:** the documentation procedure gets one source in `docs/development/` (its form is question Q3), stating the user's ruling with the template. The archivist's definition then carries only its method: how to rewrite a page to the present tense, and where each piece of the past goes. The options are the History changelog row, the Gotchas section, a collapsed "first got wrong" box at the content it concerns, or the issue that closed the work. The close-out narrative becomes an issue comment plus one History row.

### 5.4 The method-implementer refuses a vague plan

**Intent (the user):** if the implementer finds itself working around the plan or reading vague instructions, it refuses the implementation. A vague plan, typically one that needs more attention to ontology, is work for the orchestrator and the user together. *"The method implementer needs an excellent plan to not have surprises since it cannot directly communicate with me."* Implementing published formulations is useless if the foundation has not been laid, when that foundation needs significant ontological exploration and dialogue.

**Today:** the body has no readiness check; its procedure starts from "Read the plan + cited literature" and its manifest assumes a published-formulation build (D7).

**Sketch `[HYPOTHESIS]`:** the definition opens with a readiness test the plan must pass before any code:
- every object the build creates or changes is named in the plan, with its home (module and layer);
- every convention crossing a subsystem boundary is in a crosswalk;
- the gates exist as a test-architect spec;
- the done-when is a predicate;
- no step leaves a choice open ("decide", "figure out", "as appropriate", "TBD", two candidate shapes);
- the plan carries the user's ruling that it is polished (`plan-authoring` §0).

A failure is returned as `REFUSED:` with the specific questions, and the same holds mid-build: the moment the agent would have to work around the plan, it stops and returns rather than improvise. An ontological question is never settled between the implementer and the orchestrator alone; it returns to the orchestrator and the user. The identity follows (ruling R3): a builder that executes a polished plan, of which a published formulation is one kind.

### 5.5 Each definition rewritten to read as one refined mandate

**Intent (the user):** no definition reads as the original plus amendments. The amendments are used to sharpen the definition itself, which then reads seamlessly as a refined version.

**Today `[M]` 2026-09-22** (`scratch/_agent_defs/amend_tells.py`, a regex census of amendment phrasing: dates, "promoted from", "RE-POSED", "retired", "used to", "no longer", "sharpens", "does not contradict", case narratives; positive control the test-architect's "promoted from the lessons digest … 2026-09-21"; the counts are leads, each hit to be read): 9 of 9 bodies carry at least one tell, 45 in all (archivist 14, elegance-enforcer 11, explorer 8, cross-domain-attacker 5, test-architect 3, the other four 1 each).

**Sketch `[HYPOTHESIS]`:** each body is rewritten once, from its role outward (identity, capability, method, return), with the S3 fixes, the §5.1 gaps and the §5.2–§5.4 mandates folded in. Dates, cases and narratives leave the definition for the evidence page, as the harness's distillation law already requires of rules ("every clause keeps its imperative, its `check:` and its `tell:`"). Done-when: the census reads 0 on every generated body, each surviving hit on a stated exception list.

## 6. Open questions

- **Q1 → ruled R6.**
- **Q2 → ruled R7.**
- **Q3 → ruled R11: a path-scoped rule.** The suggestion as made: The user's framing: the archivist is the only sub-agent that documents, so the choice is its definition or a preloaded skill, and a skill's description would reach the main agent's roster every session. The suggestion is a third form: a **path-scoped rule** `documentation`, with `paths:` on `docs/theory/**` and `docs/architecture/**`, carrying the present-only ruling and the page template. The reasons:
  - the main agent needs it too: it writes theory pages itself in a surgical carve (W3 sends only the changelog to the archivist) and fixes stale docs on sight (Cardinal Rule 3);
  - qa needs it to verify documentation claims in W4;
  - a path-scoped rule loads for whichever agent touches a matching file, main or sub-agent, at no cost until then, with no roster line and no preload (the harness page, `[M]` 2026-09-21, the `vv-testing` probe; `plan-authoring` and `coding-standards` load this way today);
  - the archivist's definition keeps only its own method (its build gate, its cross-reference grep, how it rewrites a page to the present) and points at the rule.

  One measurement is owed in Phase 1: that a WRITE of a new file under a scoped path loads the rule as a read does (`[R]`: only the read is measured).
- **Q4 → ruled R10**: decided per agent at the rewrite, where the orchestrator gives its best suggestion from several perspectives (what the role applies at every dispatch, context cost, what the brief can carry instead, what a missing preload has cost before).
- **Q5 → ruled R8.**

## 7. Rulings ledger

- **R1** (the user, 2026-09-22): the review is endorsed as the foundation; the shape is **S3**.
- **R2** (the user, 2026-09-22): agent memory is not option (b), because curating the agents' lessons mixes the orchestrator's role with curation; (a) or (c), the orchestrator's pick. **Picked: (c).** An agent writes its own memory, and each lesson names the rule or skill clause that does not already cover it (an `omitClaudeMd` agent reads `.claude/rules/` on demand to check); the orchestrator commits the memory diff at close-out and does not curate it. Why (c) over (a): it puts the "a rule already says this" check at the writer, where CP#14 observation 6 arose, and a lesson whose uncovered clause is the agent's own definition is marked at birth as a §5.1 class-(ii) gap. An agent's own memory is always inside its edit scope, whatever the brief says of the tree.
- **R3** (the user, 2026-09-22): the method-implementer's first requirement is a well-specified plan; it refuses a vague one (§5.4). Its identity is a builder executing a polished plan, not a published-formulation specialist.
- **R4** (the user, 2026-09-22): D6 is corrected; a sub-agent can reach the orchestrator before it returns (measured, §3).
- **R5** (the user, 2026-09-22): the five items of §5 are the scope that follows the S3 foundation.
- **R6** (the user, 2026-09-22, Q1): an agent that may need clarification from the orchestrator holds `SendMessage` and asks, rather than working around what it thinks was meant; the Key agents are the first candidates. **The orchestrator's application:** all nine, since each Support agent has its own clarification case (the literature-researcher's "not in the local folder" question is W7's own; the explorer's question scope; the attacker's artefact and question). **The protocol, in every definition `[R]`:** send the question and continue the work that does not depend on the answer (the reply arrives at the next tool call); if nothing is independent of it, return with the question in `NEEDS:` and be resumed by name with the answer (invariant 3), never spin tool calls waiting. An ontological gap in a plan is not a clarification: the method-implementer refuses (§5.4).
- **R7** (the user, 2026-09-22, Q2): tests are hierarchical and well structured, so that a failure's boundary is known; the principle goes into `vv-principles`, the procedure into the test-architect's definition; the test matrix orders its rows reuse > improve > new (§5.2).
- **R9** (the user, 2026-09-22): the boundary-of-failure loop of §2 — test-architect foresight, numerics-investigator search, qa hindsight, repeating when hindsight finds a new or a vacuous boundary.
- **R10** (the user, 2026-09-22, Q4): the preload criterion is decided per agent at the rewrite, on the orchestrator's multi-perspective suggestion.
- **R11** (the user, 2026-09-22, Q3): the documentation procedure is a path-scoped rule `documentation`.
- **R12** (the user, 2026-09-22): *"The plan is scoped enough to start. Begin working on it."* Implementation opens at §8 step 1.
- **R13** (the user, 2026-09-22, Q6): the literature-researcher downloads a freely published resource itself, into the scratch literature tree, and may fetch from OSTI; a paywalled paper or a secondary substitution is still a question.
- **R14** (the user, 2026-09-22, Q7): `~/Downloads/NSE/` is NOT a sanctioned location. The user will create a place holding the NSE archive and the other literature, and give the agent access; until it exists, the definition names only `scratch/literature/` and the OCR sidecars.
- **R15** (the user, 2026-09-22, Q8): the page template gets a home in `docs/development/`; the template and the page that follows it most closely are sharpened first, then the other pages are rewritten, in a separate documentation session: #498. The archivist's rewrite here points at the `documentation` rule and at #498 until the template page exists.
- **R16** (the user, 2026-09-22, Q9): each agent's memory retirements (classes (i) and (iv) of its audit, and the archivist's 29 unindexed topic files) land in that agent's rewrite commit.
- **R17** (the user, 2026-09-22): the review stage is option 2. qa classes every landed attack NUMERICAL, ARCHITECTURAL or BOTH; the elegance-enforcer, in parallel, constructs the illegal states the types still allow; the parent routes the architectural findings to the enforcer and the constructed states to qa or the test-architect (the workflows page, W1-P3 and W2-P3).
- **R18** (the user, 2026-09-22): a guard is debt or a declared scope boundary (`SCOPE-BOUNDARY[guard]`), and either way it names the machinery that would derive its property, so that building it can be ruled on. `[LANDED 1f18f20e]` `coding-standards`, the cardinal page, a `coding-elegance` Pattern 4 corollary, the ledger test, the orbit catalogue's door as the first instance.
- **R19** (the user, 2026-09-22): probes live only in `scratch/derivations/diagnostics/`. The tracked `derivations/diagnostics/` retires: each of its files moves to scratch and survives for a stated reason, is promoted into `tests/`, or is retired; and the scratch folder is then triaged and trimmed (`tests/derivations/_promotion_policy.md`). This is the campaign's last step: `[M]` 2026-09-22 the directory holds 24 tracked files and is cited by 11 theory pages, 35 test files and one module, and `tests/derivations/test_diagnostics_resolve_their_imports.py` scans it.
- **R8** (the user, 2026-09-22, Q5): the orchestrator writes every rewritten definition and the user reviews each. The prose is direct, without mannered speech, straight to the point: maximum effect with minimum context. Its instrument: each body is generated with a `budget_tokens` set at its measured size (the harness's budget law), so growth is a red, and the §5.5 amendment census reads 0.

## 8. Implementation order and its start condition

Implementation started 2026-09-22 (R12). The order `[HYPOTHESIS]`:

1. **Measure.** Whether a Key dispatch receives the main memory index (§3's open conflict), by dispatching. `[LANDED 2026-09-22]` it does, as the session-start snapshot; and a write does not load a path-scoped rule (§3).
2. **The S3 foundation.** `[LANDED 2026-09-22]` in five commits: the whole definition generated and a code-span file check (`ff1e99d7`; the check's CI-only red on an absent ignored directory fixed in `d95a81ce`); the two write guards (`1d5d12d1`); tools matched to mandates, `SendMessage` for all nine, the archivist's name, invariants 5 and 6 and the brief template's read-only (`7ae2396c`); the `documentation` rule and the `vv-principles` ladder (`86d4529a`). The front-matter changes (tools, `SendMessage`, the write-scope hooks, the name) reach a dispatch only after Claude Code restarts; step 5's probes measure them then.
   - The generator reads the whole body. The move is verbatim first, and `--check`'s first red on the moved bodies is its positive control; D3's dead references are then fixed so the move lands green.
   - The front matter matches the mandates (tools; `SendMessage` for all nine per R6; the archivist's name).
   - A PreToolUse hook refuses writes to GENERATED files.
   - Each agent's write scope is stated and enforced.
   - The brief template defines "read-only" and states the memory scope (R2).
   - The pages the rewritten bodies will point to are written: the documentation page in the form Q3 settles, and the test-ladder clause in `vv-principles` (R7).
3. **The §5.1 audit,** per agent, producing the class-(ii) gap lists. `[HYPOTHESIS]` dispatched after the restart, so that each audit dispatch also measures what the agent now holds (step 5's probe folded in): each agent classifies its own memory with a fixed schema (entry; class (i) restates a clause, naming it; (ii) a definition gap, with the definition sentence it implies; (iii) its own experience), writes the table under `scratch/_agent_defs/audit/<name>.md`, and reports its tool list, and whether a write outside its scope was refused, in the first line of its return.
   `[LANDED 2026-09-22]` All nine audits returned; the tables are `scratch/_agent_defs/audit/<name>.md`, the orchestrator's reading of them `scratch/_agent_defs/audit/_orchestrator_notes.md`. Each dispatch doubled as the step-5 probe.

   | agent | population | (i) restates | (ii) gap | (iii) own | (iv) stale | tools as intended | write probe |
   |---|---:|---:|---:|---:|---:|---|---|
   | archivist | 213 | 54 (+20 restating its own definition) | 21 → 18 sentences | 113 | 5 | yes | generated-guard REFUSED |
   | cross-domain-attacker | 98 | 28 | 4 | 64 | 2 | yes; no zotero tool | write-scope REFUSED |
   | elegance-enforcer | 62 | 19 | 10 | 33 | 8 | yes | write-scope REFUSED |
   | explorer | 77 | 35 | 13 | 27 | 2 | yes | write-scope REFUSED |
   | literature-researcher | 81 (93 rows) | 10 | 20 → 18 sentences | 57 | 6 | yes; no zotero tool | write-scope REFUSED |
   | method-implementer | 29 | 15 | 1 (+2 index-wide) | 9 | 4 | yes | generated-guard REFUSED |
   | numerics-investigator | 63 | 22 | 10 | 31 | 0 | yes | generated-guard REFUSED |
   | qa | 159 | 67 | 7 | 78 | 7 | yes | write-scope REFUSED |
   | test-architect | 170 | 60 | 28 → 17 sentences | 76 | 6 | yes | generated-guard REFUSED |

   `[M]` 9 of 9 dispatches held Read, Write, Edit, Bash and SendMessage, and every Key agent Agent; none held Grep, Glob, ToolSearch or AskUserQuestion; 9 of 9 write probes were refused by the intended hook, and the main agent's own edit of a generated rule was refused too (2026-09-22). Population and class counts are each agent's own, over its own stated population.

   **Cross-cutting findings, for the rewrite:**
   - **The orchestrator's briefs abbreviated the generated rules list** (2 of 2 Support agents measured it: the explorer got 5 of 7 items with `instrument-doctrine` cut at X2; the attacker lost `process-discipline`'s premise clause). The template says paste the list; the rewrite makes the definitions of briefing agents say "verbatim", and the orchestrator pastes the generated list rather than retyping it.
   - **Every "grow the skill yourself" directive contradicts the new guards** (qa's trigger into the generated `vv-principles`; the attacker's Growth Protocol into `cross-domain-frames`, refused by write-scope). The rewrite: an agent proposes a skill or rule edit in its return, naming the clause it extends; the orchestrator applies it.
   - **No Zotero tool reaches a sub-agent** (2 of 2 agents that list the server; port 23119 refused). The literature tier routes around Zotero until it is measured working.
   - **Path-scoped rules reach an `omitClaudeMd` agent on a Read** (the attacker, 2 of 2).
   - **The write-scope hook sees Edit, Write and MultiEdit only**; `mcp__nexus__rename` (applied) and `ingest`/`runtime_ingest` write files (`[R]`, the explorer). Add them to the read-only agents' matcher at the rewrite.
   - **Correctness defects inside definitions**, each fixed at its agent's rewrite: the literature-researcher's direction-cosine mapping is inverted for both sources it names; the numerics-investigator's Step 3 (uniform reflective box) is exactly singular for DD at d ≥ 2; the elegance-enforcer's institutional-knowledge #1 puts fission inside the loss operator and #5 is refuted on Python 3.14.3; the test-architect's `n_inner=None` gap and `strict=False` xfail are stale, and `pyproject.toml:33`'s "run sentinels WITHOUT -O" is wrong for collected tests; the archivist's `vv-status` line implies a value (`implemented`) that is not legal.
   - **Identities narrower than the work**: the elegance-enforcer's largest workload is documentation, harness prose and generator tools, not code; qa's is prose artefacts in 9 of 19 reviews; the test-architect's is operator-algebra carves, not solvers; the method-implementer adds two readiness conditions to §5.4 (a plan built on an unmerged branch is refused; the implementer never changes git state in the shared tree).
   - **Filed**: #497, the xref gate certifies `:mod:` targets only (a dead `:class:` target is DECLINED).
   - **Held, outside this campaign**: the method-implementer's L-013 (a separability claim possibly refuted by `R = R_spatial ⊗ A_angular`) is re-measured before it moves anywhere.

   **Questions for the user** (collected from the audits): (Q6) may the literature-researcher fetch a freely published primary standard (ENDF-102, the NJOY manual) directly, asking only for a paywalled paper or a secondary substitution? (Q7) is `~/Downloads/NSE/` (760 NSE volume zips) a sanctioned second tier-0 location beside `scratch/literature/`? (Q8) should the #231 page template get a home in `docs/development/`, or does the `documentation` rule's placement law suffice? (Q9) the memory retirements (classes (i) and (iv), and the archivist's 29 unindexed topic files) — apply them in each agent's rewrite commit?
   **The rewrite** `[LANDED 2026-09-22]`, 9 of 9, each approved by the user, then merged: method-implementer `62ab8f6e`; test-architect `324ba447` (+ `c2a7418b` the qa-attack standard, `69bd323e` the spec as a declared DAG); qa `ea2e5c7d`; elegance-enforcer `af479efb` (+ `1f18f20e` the guard classes); numerics-investigator `54809381`; archivist `ba3e0dab`; explorer `43defe8d` (+ the write-scope hook refusing Nexus's writers); literature-researcher `5ac5ffe8`; cross-domain-attacker `69bfc44e`. Every definition is generated, budgeted at its measured size, and reads 0 on the amendment census bar present-tense uses. One red on `main` from the qa rewrite (a test anchored on a removed heading), fixed at `a49ec558`; CI is watched after every push.
4. **The rewrite (§5.5),** one agent per commit, folding in D4, D5, D7, the gap list and the §5.2–§5.4 mandates, with the amendment census at 0.
5. **Retire `derivations/diagnostics/` (R19)**: triage each tracked probe (move to scratch with its reason, promote, or retire), repoint the 11 pages, 35 tests and one module that cite the directory, retire or re-aim its import-resolution test, then triage and trim the scratch folder.
6. **Confirm by dispatching:** each agent is probed for what it received (tools, write scope, memory policy).

## Resume surface

This file; `docs/development/harness.md` ("Adding or changing", "What loads, and what it costs"); `docs/development/workflows.md` ("The brief"); `.claude/agents/*/AGENT.md`; `harness_context_budget.md` ⏸ COMPACTION POINT #14 (the census and the eight dispatch observations this review extends).


# ⏸ COMPACTION POINT — 2026-09-22. The nine definitions are rewritten; next is §8 step 5, retiring `derivations/diagnostics/`

**The instruction (the user, verbatim):** *"I approve. Merge and prepare for context compaction before the last step."* The last step is R19 (§8 step 5).

**Where things stand.** `main` = `origin/main` at this commit's parent `69bfc44e` plus this plan update; CI green on `69bfc44e`. Landed this campaign (§8): the S3 foundation (whole definitions generated; the code-span file check; the generated-guard and write-scope hooks, the latter also refusing Nexus's applied `rename` and its `ingest` tools; tools matched to mandates; `SendMessage` for all nine; invariants 5 and 6; the `documentation` rule; the `vv-principles` ladder); the §5.1 audit; the nine rewrites; rulings R1–R19. Beside it: `SCOPE-BOUNDARY[guard]` (R18), the P3 review routing (R17), `rests_on` recorded in #358's plan, #497 (xref checker blind to `:class:` targets) and #498 (the page template, a separate documentation session) filed.

**Next: §8 step 5 (R19).** Probes live only in `scratch/derivations/diagnostics/`; the tracked `derivations/diagnostics/` retires. `[M]` 2026-09-22: 24 tracked files; cited by 11 theory pages, 35 test files and one module (`orpheus/derivations/continuous/peierls_nystrom/geometry.py`); `tests/derivations/test_diagnostics_resolve_their_imports.py` scans the directory; `.claude/plans/**` and `.claude/scratch/**` cite it as archaeology (not repointed). The procedure `[HYPOTHESIS]`, per `tests/derivations/_promotion_policy.md` (DELETE / PROMOTE / LEAVE):
1. Re-take the census (the counts above are from one `git grep`; re-run in Python with a positive control, X2).
2. Triage each tracked probe: promote into `tests/` (test-architect), move to `scratch/derivations/diagnostics/` with a stated reason, or retire (git history keeps it). A likely delegate: the numerics-investigator, which owns probes, with the orchestrator reviewing.
3. Repoint each citing page and test (a citation of a retired probe becomes its commit or its issue); retire or re-aim the import-resolution test.
4. Triage and trim `scratch/derivations/diagnostics/` itself (33 files).
5. Then §8 step 6: nothing further to probe (the audit dispatches measured it), so the campaign closes: the close-out, the plan archived, the memory index updated.

**Open, outside this step:** the literature location the user will create (R14); the method-implementer's 112 unindexed closeout files (asked, not yet ruled; R16 covered only classes (i) and (iv)); L-013's separability claim, owed a re-measure.
