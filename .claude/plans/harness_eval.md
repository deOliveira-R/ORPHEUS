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
- **H6, the environment is in scope.** "If we find tool schemas and listings dominate then we shall try to trim it. We will figure out exactly how to do that."
- **H7, Fable's place.** Fable is tested only on the most cognitively loaded tasks; its main function is orchestrating and advising. "There is probably no point in evaluating Fable if Opus or Sonnet are passing. We will have to figure out what hard means to Sonnet and Opus, and hill-climb them with coaching behaviour." The orchestrator's reading `[R]`: E2 is a ladder: a task goes to Sonnet first, to Opus where Sonnet fails, to Fable where Opus fails, and coaching is evaluated as what lifts a model over a task it fails bare.
- **H8.** The accounting stage (E0) may run. It ran; results below.

## Baseline `[M]` 2026-11-13

**What an agent holds before it does anything** (first-turn input tokens, read from the transcripts of the two workflow runs of this session, `wf_04ea90f3-46f` and `wf_99d2fbf4-c0a`):

| agent, as dispatched | model | n | first-turn context |
|---|---|---|---|
| `explorer` (a Support agent: `omitClaudeMd`, a narrow tool allowlist) | Sonnet | 16 | 55.5K (one outlier at 125K) |
| `elegance-enforcer` (a Key agent: CLAUDE.md, rules, three preloaded skills) | Opus | 16 | 98.7K |
| the default workflow agent (every tool, CLAUDE.md, rules) | Sonnet | 12 | 124.3K |
| the same | Opus | 1 | 123.6K |

**What the project's own text weighs** (characters / 3.6; `[REFUTED 2026-11-13]` as a token count by E0 below, which measured 2.6 characters per token; the character proportions stand):

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

Four skills carry most of the preload: `vv-principles` 13.8K (in 6 agents), `numerical-bug-signatures` 11.3K (2), `algebra-of-record` 9.8K (3), `coding-elegance` 9.6K (5). Whether an agent's memory loads on dispatch: measured in E0 below (it does, from the repository root).

**What a run costs** (medians per agent in the census run): an Opus reviewer took 120 turns, grew from 98.7K to 261K of context, wrote 90K tokens, and read 24.2M tokens from cache; a Sonnet census agent took 49 turns, 173K final context, 48K written, 5.9M read from cache. The cost of an agent is therefore close to (context) x (turns): every 1K removed from the first turn is removed from every later turn.

**What these numbers already say** `[R]`:
- The project's always-on text (15.6K) is a small part of a 124K first turn. The default workflow agent starts 69K heavier than `explorer`; the differences between the two are the tool allowlist and CLAUDE.md with the rules, and the rules are at most 15.6K of that. So most of the gap is `[HYPOTHESIS]` tool schemas and listings (the MCP servers, the skill and agent listings), not the project's prose. This is the first thing to measure, because if it holds, the largest saving is to dispatch every workflow stage as a typed agent with a narrow allowlist, and no instruction needs to be evaluated to get it.
- For the Key agents the preloaded skills (41K to 47K) are the largest project-owned block, three times the always-on text.
- The previous harness campaign (`harness_context_budget.md`, closed 2026-09-22) measured the same fixture across harness versions and found the movement inside its own noise, at one draw per cell. An evaluation with one run per arm will again decide nothing.

## E0 results `[M]` 2026-11-13 (accounting; ruled to run by the user the same day)

**Instruments.** (1) A workflow of 24 one-turn probes (`wf_f6402e54-a5b`): each of the nine project agents, `general-purpose` and the default workflow agent, on Sonnet and on Opus, with one fixed prompt and no tool call; first-turn usage read from each transcript as input + cache-creation + cache-read tokens. (2) Ten throwaway probe agents passed to a headless Sonnet session by `--agents`, each varying one block (the tool list, the Nexus server, the browser server, CLAUDE.md with the rules, one preloaded skill, the agent memory). (3) Twelve headless main sessions with one flag changed each. The transcripts record the full request (every tool schema, every instruction file, every preloaded skill), so each block's character size is read, not estimated. Scripts and data: `scratch/_harness_eval/e0/`. Nothing tracked and no memory file changed (`git status`, the memory directories' modification times).

**Controls.** Three repeats of one cell differ by 14 tokens. Sonnet and Opus differ by 6 tokens in each of the 11 agent types, once both are dispatched from the repository root. A block of known size (the 49 582-character `vv-principles` file appended to a main session) moved the reading by 18 671 tokens.

**The earlier estimates were low by a third.** The baseline tables above convert characters to tokens at 3.6. Measured: 2.6 characters per token for this project's prose, 2.6 to 2.9 for tool schemas. Every "tokens" figure in the two project-text tables above is superseded by the table below.

**The first turn of each agent, dispatched by a workflow, by block** (thousands of tokens, Sonnet; the block columns sum 5K to 6K short of each total, the remainder being the brief, the dispatch wrappers, the session context and the fit's constant; the total is measured, the split is the recorded character size of each block converted by a two-rate fit over all 34 sub-agent probes, largest residual 2.4K):

| agent | total | Nexus tools | other MCP | built-in tools | CLAUDE.md and rules | MEMORY.md | preloaded skills | AGENT.md and own memory | listings and environment |
|---|---|---|---|---|---|---|---|---|---|
| qa | 136.7 | 17.9 | 0 | 4.4 | 22.2 | 4.6 | 69.3 | 10.9 | 1.4 |
| test-architect | 135.4 | 17.9 | 0 | 4.4 | 22.2 | 4.6 | 65.8 | 14.0 | 1.4 |
| numerics-investigator | 130.5 | 17.9 | 0 | 4.4 | 22.2 | 4.6 | 61.9 | 11.9 | 1.4 |
| method-implementer | 125.0 | 17.9 | 0 | 4.4 | 22.2 | 4.6 | 61.6 | 8.6 | 1.4 |
| archivist | 128.9 | 17.9 | 0 | 4.4 | 22.2 | 4.6 | 62.6 | 11.3 | 1.4 |
| general-purpose | 122.9 | 17.9 | 13.3 | 45.9 | 22.2 | 4.6 | 0 | 1.4 | 12.9 |
| the default workflow agent | 122.6 | 17.9 | 13.3 | 45.9 | 22.2 | 4.6 | 0 | 0.9 | 12.9 |
| elegance-enforcer | 96.4 | 17.9 | 0 | 4.4 | 22.2 | 4.6 | 26.6 | 14.9 | 1.4 |
| explorer | 52.6 | 17.9 | 0 | 4.4 | 0 | 0 | 12.6 | 11.8 | 1.4 |
| cross-domain-attacker | 50.9 | 17.9 | 0 | 5.3 | 0 | 0 | 6.3 | 14.9 | 1.4 |
| literature-researcher | 38.6 | 0 | 0 | 5.3 | 0 | 0 | 9.1 | 16.1 | 1.4 |

**Direct differences from the one-block probes** (tokens; these are measured, and where they differ from the fitted columns above they are the better number):

| block | tokens |
|---|---|
| an agent holding one tool (`Read`), no CLAUDE.md | 5 271 |
| four more built-in tools (`Write`, `Edit`, `Bash`, `SendMessage`) | +4 140 |
| the Nexus server's 45 tools | +20 071 |
| the browser server's 25 tools | +6 989 |
| CLAUDE.md (project and user), the seven rules and MEMORY.md | +26 528 |
| `vv-principles` preloaded (64 070 characters as injected) | +25 590 |
| `coding-elegance` preloaded | +12 856 |
| agent memory enabled, on an agent with no memory files | +4 051 |
| every tool a headless session holds (87), with the skill and agent listings that come with them, against five built-ins | +47 815 |
| the same all-tools agent inside this interactive session (126 tools: the artifact, design, scheduling and messaging tools and the claude.ai connectors are added) | +39 181 more |

Preloaded skills by weight: `vv-principles` 24.6K (in 5 agents), `numerical-bug-signatures` 15.3K (2), `algebra-of-record` 13.1K (3), `coding-elegance` 13.0K (5), `research` 9.1K (1), `retirement-audit` 8.0K (6), `doc-issue-relocation` 6.4K (1), `nexus-elegance` 5.6K (1), `instrument-doctrine` 4.8K (4); the remaining seven are under 3.5K each. The largest single tools of the default agent: `Artifact` 12.6K, `ArtifactData` 4.3K, `DesignSync` 3.7K, `Monitor` 2.7K, `ArtifactComments` 2.6K. Instruction files: `process-discipline` 4.7K, MEMORY.md 4.6K, CLAUDE.md 4.1K, `workflows` 3.0K, `code-search` 2.4K, `cardinal` 2.2K, `nexus-tools` 2.1K, `instrument-doctrine` 1.8K, the user's CLAUDE.md 1.1K, `articulation` 0.9K.

**The orchestrator's own floor** (a headless main session, first turn; an interactive session such as this one holds about 39K more of tools): Sonnet 103 471, Opus 103 590, Fable 104 899. Of the Sonnet figure: 27 built-in tools 42.0K, Nexus 20.9K, the browser server 7.0K, the skill listing (90 entries) with the `Skill` tool 4.9K; with no tools, no MCP server and no skills 34.1K remain (the system prompt, CLAUDE.md, the rules, MEMORY.md, the agent listing, the session-start hook's text).

**What E0 settles.**
1. *The hypothesis holds for an untyped agent and fails for the Key agents.* In the default workflow agent the tool schemas are 77K of 122.6K and the project's always-on text is 26.8K. In the four heaviest Key agents the preloaded skills are 62K to 69K, half of the first turn, and the tools are 22K.
2. *The model does not change the weight.* What differs between Sonnet, Opus and Fable is what each does with the text, which is E2's question.
3. *One saving needs no evaluation, and two need a routing check.*
   - Free: the artifact, design, scheduling, messaging, connector and browser tools of an untyped agent, about 55K of its 77K of tools. No rule or skill of the project refers to them. A workflow stage dispatched as a typed agent with a short tool list does not carry them.
   - Routing check: the Nexus server costs 20.1K in each of the ten agent types that name it, whether or not the stage asks a graph question. The always-on `nexus-tools` and `code-search` rules and seven skills tell an agent to use those tools, so removing or narrowing them changes a capability the text relies on. The `Skill` tool is the same case: the measured 35.9K probe (five built-in tools with CLAUDE.md and the rules) holds neither, while the rules it inherits tell it to load skills.
   - Routing check: MEMORY.md (4.6K) is the orchestrator's index of campaigns; whether a sub-agent needs it is not a coaching question.
4. *What E2 must decide is now sized.* The preloads: `vv-principles` is 24.6K in five agents, and the four largest skills are 61K of qa's 69K. The always-on prose: 22.2K. These are the units whose removal needs evidence.
5. *The dispatch route changes the tools* (read from the recorded requests, 0 of 27 workflow probes and 2 of 2 `Agent`-tool probes). A workflow-dispatched agent does not hold the `Agent` tool, including the six Key agents whose definitions list it: inside a workflow a Key agent cannot spawn, which invariant 2 of the `workflows` rule does not say. Dispatched by the `Agent` tool, qa holds it (137.8K, 53 tools). An untyped agent holds 126 tools by the workflow route and 114 by the `Agent` tool (110.1K against 122.9K). No agent on either route holds `ToolSearch`; every MCP schema is loaded eagerly.
6. *The working directory changes what an agent loads, silently.* Five probes were dispatched while this session's working directory was a scratch sub-folder. `numerics-investigator` then read 120.7K against 130.2K from the repository root, and `method-implementer` 117.9K against 125.0K: the preloaded `vv-principles` arrived 15 273 characters shorter and the agent's own memory was absent. The untyped agent was unaffected (14 tokens). The build workflows dispatch into worktrees, so E4 must check what a worktree-rooted agent loads.
7. *The budget gate converts at the old rate.* `tools/harness/budget.py` uses 3.6 characters per token, measured on 2026-09-20 at 3.51 to 3.53. The same kind of text now measures 2.6, so `--check` under-reads by about a third `[R]` (the likeliest cause is the newer models' tokenizer; not verified).
8. *The probes' own reports of their tool counts were wrong in 2 of 24 cells* (36 and 46 reported, 53 and 51 held). Inventories are read from the recorded request only.

**Residue.** Fourteen headless sessions left transcripts in the user's Claude project folder (they hold the ablation data). The probe agent with memory enabled created an empty `.claude/agent-memory/p5-five-memory/`, removed.

## The design `[R]` (proposed 2026-11-13; nothing below is ruled)

**E0, accounting.** No judgement involved. For each agent type and each model, a zero-tool probe reads the first-turn usage, and the floor is decomposed by removing one block at a time: the MCP servers, the skill listing, the tool allowlist, CLAUDE.md with the rules, each preloaded skill, the AGENT.md body, the agent's memory. Output: a table of tokens per block per agent, and which blocks the project controls. (As run: throwaway agents passed to a headless session by `--agents` vary every block, CLAUDE.md included through `omitClaudeMd`, with no separate checkout.)

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
- *Arms*, per model (Sonnet first, Opus where Sonnet fails, Fable where Opus fails: ruling H7): bare (the unit absent), with the unit, and with a one-sentence version of the unit. A unit is retired for a model when the bare arm does as well; sharpened when the one-sentence version does as well as the full text; kept otherwise.
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

1. **The free saving (E0, item 3).** Dispatch every workflow stage as a typed agent with a short tool list now, before E1 and E2? It changes no instruction. The Nexus and `Skill` questions wait for the routing check.
2. **How the Nexus tools reach an agent.** A per-role allowlist of graph tools, or the command line through `Bash` with the skills rewritten to it. The second is a change upstream in sphinxcontrib-nexus.
3. **The scale of E2**, after the pilot.

## Implementation starts when

The user rules the design. E0 has run.
