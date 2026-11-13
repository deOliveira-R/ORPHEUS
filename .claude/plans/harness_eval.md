# What each model still needs from the harness: evaluate, sharpen or retire (living plan, opened 2026-11-13)

## Goal

Every instruction the harness loads into an agent changes an outcome for the model that reads it, and nothing else is loaded. An agent in a workflow then starts from the smallest context that still produces work to the project's standard, and each workflow stage is given to the model that does it well at the lowest cost.

The harness here means: CLAUDE.md, the rules (`.claude/rules/`), the skills (`.claude/skills/`), the agent definitions (`.claude/agents/*/AGENT.md`), the agents' memory, the session-start protocol and the hooks. Their sources are `docs/development/`; the `.claude/` copies are generated and are never edited by hand.

## Rulings (the user, 2026-11-13)

- **H1, the premise.** The harness was built for models since Opus 4.6. "Now you're all much smarter/capable than you were when we created this harness, so some of the skills, agents definition and other parts of the harness might contain totally unnecessary instructions."
- **H2, the method.** "A pure eval design, evaluate, sharpen or retire of certain parts of the harness to keep only what is really important." Not a rewrite by judgement.
- **H3, the reason.** Unnecessary text costs tokens, may hold a model back, and multiplies in a workflow.
- **H4, per model.** "We need to figure out what Sonnet, Opus and Fable needs, and what they don't. This might also help us figure out how to assign tasks to models for the dynamic workflow."
- **H5, order.** This runs before the ontology chart's build workflows (`ontology_chart.md`, ruling O11).

## Baseline `[M]` 2026-11-13

**What an agent holds before it does anything** (first-turn input tokens, read from the transcripts of the two workflow runs of this session, `wf_04ea90f3-46f` and `wf_99d2fbf4-c0a`):

| agent, as dispatched | model | n | first-turn context |
|---|---|---|---|
| `explorer` (a Support agent: `omitClaudeMd`, a narrow tool allowlist) | Sonnet | 16 | 55.5K (one outlier at 125K) |
| `elegance-enforcer` (a Key agent: CLAUDE.md, rules, three preloaded skills) | Opus | 16 | 98.7K |
| the default workflow agent (every tool, CLAUDE.md, rules) | Sonnet | 12 | 124.3K |
| the same | Opus | 1 | 123.6K |

**What the project's own text weighs** (characters / 3.6):

| block | tokens |
|---|---|
| always-on: CLAUDE.md and seven generated rules, plus the installed Nexus rule (`python -m tools.harness --check`) | 15.6K |
| path-scoped rules, when a matching file is read: `plan-authoring` 9.7K, `coding-standards` 3.3K, `vv-testing` 1.7K, `documentation` 0.6K | up to 15.3K |
| the session-start batch for the main agent: `vv-principles` 13.8K, `coding-elegance` 9.6K, `instrument-doctrine` 3.6K, the lessons index 3.4K | 30.4K |
| all 21 project skills | 88.6K |
| the memory index (`MEMORY.md`) | 3.4K |

| agent | AGENT.md | preloaded skills | own memory index and lessons |
|---|---|---|---|
| qa | 2.5K | 46.8K (8 skills) | 14.6K |
| test-architect | 3.2K | 44.2K (6) | 20.6K |
| archivist | 2.2K | 41.9K (7) | 14.6K |
| numerics-investigator | 2.9K | 41.5K (7) | 16.7K |
| method-implementer | 2.0K | 41.3K (5) | 3.5K |
| elegance-enforcer | 3.1K | 19.6K (3) | 21.5K |
| explorer | 1.7K | 9.3K (3) | 6.8K |
| literature-researcher | 2.2K | 6.8K (1) | 5.9K |
| cross-domain-attacker | 2.1K | 4.8K (2) | 18.8K |

Four skills carry most of the preload: `vv-principles` 13.8K (in 6 agents), `numerical-bug-signatures` 11.3K (2), `algebra-of-record` 9.8K (3), `coding-elegance` 9.6K (5). Whether an agent's memory index and lessons load on dispatch is not measured here.

**What a run costs** (medians per agent in the census run): an Opus reviewer took 120 turns, grew from 98.7K to 261K of context, wrote 90K tokens, and read 24.2M tokens from cache; a Sonnet census agent took 49 turns, 173K final context, 48K written, 5.9M read from cache. The cost of an agent is therefore close to (context) x (turns): every 1K removed from the first turn is removed from every later turn.

**What these numbers already say** `[R]`:
- The project's always-on text (15.6K) is a small part of a 124K first turn. The default workflow agent starts 69K heavier than `explorer`; the differences between the two are the tool allowlist and CLAUDE.md with the rules, and the rules are at most 15.6K of that. So most of the gap is `[HYPOTHESIS]` tool schemas and listings (the MCP servers, the skill and agent listings), not the project's prose. This is the first thing to measure, because if it holds, the largest saving is to dispatch every workflow stage as a typed agent with a narrow allowlist, and no instruction needs to be evaluated to get it.
- For the Key agents the preloaded skills (41K to 47K) are the largest project-owned block, three times the always-on text.
- The previous harness campaign (`harness_context_budget.md`, closed 2026-09-22) measured the same fixture across harness versions and found the movement inside its own noise, at one draw per cell. An evaluation with one run per arm will again decide nothing.

## The design `[R]` (proposed 2026-11-13; nothing below is ruled)

**E0, accounting.** No judgement involved. For each agent type and each model, a zero-tool probe reads the first-turn usage, and the floor is decomposed by removing one block at a time: the MCP servers, the skill listing, the tool allowlist, CLAUDE.md with the rules, each preloaded skill, the AGENT.md body, the agent's memory. Output: a table of tokens per block per agent, and which blocks the project controls. (CLAUDE.md, the rules and the memory index are a session-start snapshot, so removing them needs a separate checkout and a headless session, as the T6 fixture did; an AGENT.md body and a skill preload change per dispatch.)

**E1, classification.** Every unit (a rule clause, a skill section, an AGENT.md section) is sorted by kind, because only one kind needs a behavioural evaluation:

| kind | example | disposition |
|---|---|---|
| a fact the model cannot infer | `grep` here is ugrep and an anchor inside a group matches nothing; the layer table; a ruling of the user; a path | kept, stated once, in the fewest words |
| prose that a hook or a gate already enforces | "never `git add -A`" (the hook refuses it); generated files (a hook refuses the edit) | reduced to what the refusal does not say |
| routing | which agent, which tool, which workflow | kept where a wrong route was measured; else retired |
| history and evidence | a founding case told in the rule | already on the evidence pages; removed from the loaded text |
| behaviour coaching | "ask of every gate: what input would make it red?"; "name the intermediate" | evaluated (E2) |

**E2, the behavioural evaluation**, for coaching units only.
- *Tasks* come from the project's own record: each rule clause links the founding case on an evidence page (`docs/development/evidence/`), which states the situation, the wrong action a model took and the right one. A task re-creates the situation on a pinned commit. The T6 review fixture (ten known findings, scored blind, `scratch/_harness_eval/`) is reused as the first review task.
- *Arms*, per model (Sonnet, Opus, Fable): bare (the unit absent), with the unit, and with a one-sentence version of the unit. A unit is retired for a model when the bare arm does as well; sharpened when the one-sentence version does as well as the full text; kept otherwise.
- *Repeats*: at least 5 per cell, because one draw decides nothing (the T6 record).
- *Scoring*: a deterministic check wherever the outcome is a fact (did the agent take the wrong action; did the gate it wrote redden under the mutation), and a blind Opus reader with a rubric and a quoted line per cell elsewhere.
- *Controls*: one task where an instruction is known to matter (the ugrep anchor: a fact, so the bare arm must fail), which proves the evaluation can see a difference; and one where it is known not to.
- *The decision is per model*, so a unit can be dropped from a Sonnet executor's brief and kept for nothing, or kept for Sonnet and dropped for Opus. The per-model table is also the assignment table for the workflow stages.

**E3, the trim.** Edits go to the sources under `docs/development/`; `python -m tools.harness` regenerates; the Nexus rule and skills are trimmed upstream in sphinxcontrib-nexus. The budget check (`--check`) and the harness gates stay green. E0 is re-run and the before-and-after table is published with one instrument on both sides.

**E4, the assignment.** From E2's per-model results and E0's costs: which model each stage of the four workflow templates is given, and which typed agent (with which allowlist and which preload) it runs as.

**Sizing** `[R]`: E0 is one small workflow of zero-tool probes. E1 is one census workflow over about 150K tokens of harness text. A pilot of E2 with 6 tasks, 3 arms, 2 models and 5 repeats is 180 short runs; the full matrix is set after the pilot shows the noise.

## What is not evidence

The orchestrator's own sense of which instructions it no longer needs. A model cannot observe what it would have done without an instruction it has already read.

## Open questions for the user

1. **Scope of the environment.** If E0 confirms that tool schemas and listings dominate the first turn, part of the saving is outside the repository: the MCP servers and plugins configured at the user level (for example the Nexus tools every agent loads, the browser tools, and the skill listings of other projects). Is trimming that environment for this project in scope?
2. **Fable's role in the matrix.** Evaluate Fable as an executor and reviewer too, or only as orchestrator and advisor?
3. **The scale of E2**, after the pilot.

## Implementation starts when

The user rules the design. E0 changes nothing and can run first.
