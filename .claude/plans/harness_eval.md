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
- **H9, typed roles.** "The typed agents dropping skills and tools inappropriate for their role in the workflow seems like a very good idea."
- **H10, Nexus per role.** A per-role allowlist and the command line are both on the table. "What matters to me is that you and the agents have the tools they need (and some might need 0 Nexus tools during the workflow)."
- **H11, the roles are accepted** (2026-11-13). "I agree with the roles. Seem like pretty good design to me." The command line for the graph is tried, and whether it works is part of the evaluation: "if the Nexus guide skill does a good job maybe it will work well."
- **H12, spawning.** The user does not mind either way whether a Key agent holds `Agent`. "Inside of a dynamic workflow, an agent having to call a support agent is a symptom of a vague task description. Arguably it should send a message to the orchestrator." The orchestrator's decision `[R]`: the headers keep `Agent` for a dispatch outside a workflow; inside a workflow no agent holds it (measured) and the `workflows` rule now says a need for another agent is a question to the orchestrator.
- **H13, the environment facts are evaluated too.** Whether the two stage roles drop CLAUDE.md and the rules is decided by evaluation. "Since Opus 5, you have all become much better at using tools appropriately. It is possible that the grep trap doesn't hold anymore, or needs just a soft steering." So E1's first row (a fact the model cannot infer is kept) is amended: a fact about the environment is kept only if the bare arm fails on it.
- **H14, go ahead** (2026-11-13, after the user restarted the harness). "Go ahead with the procedure." The orchestrator's reading `[R]`: E1, the task ladder and a pilot of E2 may run, since all three are measurement; a trim of any instruction comes back to the user per unit with its evidence, as the roles did.
- **H15, the outside skills** (2026-11-13). The gitnexus, moose and graphify skills "should not exist"; edgartools in this project "makes no sense". Done: the 15 test folders moved to `~/.claude/skills-removed-2026-11-13`; `.claude/settings.local.json` switches off `edgartools`, the `engineering` plugin and the claude.ai skills, plugins and connectors for this project. A fresh session lists 44 skills (from 90) and starts at 99.7K (from 103.5K). The browser server stays.
- **H16** (2026-11-13). The decision rule for E2 is accepted. The expensive ablation of the large skills is held until the cheap tasks are in.

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

## E0b: what each role used, and what a narrow role costs `[M]` 2026-11-13

**The usage record.** Every sub-agent transcript of this project with at least three turns, 2026-09-04 to 2026-11-13: 922 agents (`scratch/_harness_eval/e0/usage.py`). Tool calls are read from the transcripts, not reported by the agents.

| agent type | agents | tool calls | of them Nexus | agents with any Nexus call | Nexus tools that carry the use |
|---|---|---|---|---|---|
| qa | 171 | 8 493 | 35 (0.4%) | 17 of 168 | `dead_references` 12 |
| elegance-enforcer | 143 | 4 643 | 20 (0.4%) | 12 | `twin_paths` 4, `discriminations` 3 |
| test-architect | 120 | 14 555 | 19 (0.1%) | 6 | `query` 8, `callers` 5 |
| explorer | 120 | 7 194 | 256 (3.6%) | 69 of 118 | `callers` 91, `query` 54, `twin_paths` 15, `context` 15, `session_briefing` 15 |
| archivist | 110 | 9 009 | 200 (2.2%) | 65 | `dead_references` 126, `query` 29 |
| general-purpose | 98 | 4 586 | 1 | 1 | |
| the untyped workflow agent | 63 | 3 439 | 4 | 1 | |
| cross-domain-attacker | 38 | 1 534 | 18 (1.2%) | 7 | `query` 7 |
| numerics-investigator | 25 | 1 345 | 5 | 1 | |
| method-implementer | 11 | 720 | 1 | 1 | |
| literature-researcher | 9 | 576 | 0 | 0 | |

- `Bash` is 86% to 93% of every role's calls. 23 of the 45 Nexus tools were ever called by a sub-agent.
- The `Agent` tool was called by 6 Key agents out of the 580 recorded (qa 1, test-architect 4, archivist 1).
- The untyped workflow agents loaded a skill 9 times in 63 agents (`retirement-audit` 5, `coding-elegance` 4).
- The `nexus` command line has a subcommand for each graph question (`nexus --help`: `callers`, `query`, `context`, `impact`, `dead-references`, `twin-paths`, `staleness` and the rest), and `Bash(.venv/bin/nexus:*)` is already allowed.

**Role probes** (throwaway agents by `--agents`, each with CLAUDE.md and the rules; first-turn tokens):

| shape | tokens |
|---|---|
| `Read`, `Write`, `Edit`, `Bash` | 33 944 |
| the same with `SendMessage` | 35 958 |
| the same four with the `Skill` tool (it brings the listing of 90 skills, 46 of them from the user-level plugins and other projects) | 45 050 |
| `Read`, `Bash` and four Nexus tools named in `tools:` | 35 369 |
| the same with `mcpServers: nexus` added | 53 431 |

A Nexus tool named in `tools:` is held alone; `mcpServers: nexus` loads all 45 whatever `tools:` says. Every project agent that uses the graph declares the server today, which is why each pays 20.1K.

## The roles (proposed 2026-11-13, accepted the same day: ruling H11)

**New typed agents for the workflow stages that are untyped today** (the build, amend and verify stages of `build-units`; measure, draft and structural check of `spec-check`; the drafts of `design-panel`):

| role | stages | tools | first turn |
|---|---|---|---|
| executor | build, amend, measure, draft | `Read`, `Write`, `Edit`, `Bash`, `SendMessage` | 36.0K (from 122.6K) |
| checker | verify, structural check | `Read`, `Bash`, `Write`, `SendMessage` | about 35K |

Neither holds a Nexus tool or the `Skill` tool: a graph question goes through the command line, and a skill a stage needs is named in its brief and read as a file. Whether these roles also drop CLAUDE.md and the rules (26.5K more) is E2's question, not decided here.

**The nine existing agents: the server declaration is replaced by named tools.**

| agent | Nexus tools kept | share of its recorded Nexus calls | saving |
|---|---|---|---|
| explorer | `callers`, `query`, `context`, `impact`, `twin_paths`, `session_briefing`, `neighbors`, `dead_references`, `file_brief` | 219 of 256 | about 16K |
| archivist | `dead_references`, `query`, `context`, `staleness`, `provenance_chain`, `errors` | 184 of 200 | about 17K |
| qa, elegance-enforcer, test-architect, numerics-investigator, method-implementer, cross-domain-attacker | none; the command line through `Bash` | | 20.1K each |

The rest of each agent's graph use goes through the command line. The Nexus skills these agents preload name the MCP tools, so the six that lose the tools need the skills' command-line spelling (one line per skill, or the upstream change). The preloaded skills themselves stay as they are until E2.

**Not proposed here:** removing the `Agent` tool from the Key agents (6 uses in 580, but invariant 2 is a ruling), and trimming the user-level skill listing (the environment, H6; it is 11.1K wherever `Skill` is held).

## E3a: the roles, built `[M]` 2026-11-13 (branch `chore/harness-typed-roles`, commit `5015137f`, cut from `main`; not merged)

- Two new agents, `executor` and `checker`, defined at `docs/development/agents/`; the seven untyped stages of `build-units`, `spec-check` and `design-panel` now name them. Neither has memory, a Nexus tool or the `Skill` tool.
- The nine existing headers no longer declare the Nexus server. `explorer` names nine graph tools and `archivist` six; the other six agents hold none and their definitions name the command line.
- `python -m tools.harness --check`: 27 targets, 0 problems, 0 drifted. `tests/gates/test_harness_generated.py`: 2 passed.

First turn, re-measured from the branch by one headless session dispatching each agent with the `Agent` tool (tokens). The "before" column is the workflow-route measurement of E0, except qa, which was also measured by the `Agent` tool before (137.8K); the two routes differ by about 1K for a typed agent.

| agent | before | after |
|---|---|---|
| an untyped workflow stage, now `executor` | 122.6K | 36.3K |
| an untyped verify stage, now `checker` | 122.6K | 35.9K |
| qa | 136.7K | 117.8K |
| test-architect | 135.4K | 116.3K |
| numerics-investigator | 130.5K | 111.6K |
| archivist | 128.9K | 112.7K |
| method-implementer | 125.0K | 106.4K |
| elegance-enforcer | 96.4K | 77.3K |
| explorer | 52.6K | 36.4K |
| cross-domain-attacker | 50.9K | 29.7K |
| literature-researcher | 38.6K | 37.5K |

Not yet verified: that an `executor` reads a named skill file when its brief names one, and that the six agents use the command line when a graph question arises. Both are behaviour, and are E2 tasks. `docs/development/harness.md` still carries the 2026-09 cost table; it is rewritten from one instrument at E3.

## E1 and the E2 pilot `[M]` 2026-11-13

**The typed roles are merged** (`main` `5015137f`, CI green) **and are not live in this session.** The session was restarted from a checkout that did not yet hold the merge: `executor` and `checker` are unknown to it, and qa and explorer still hold all 45 Nexus tools (136.8K and 52.7K, workflow `wf_79a06d1b-fc8`). The branch now holds the merge (`d99f8504`); one more restart makes them live.

**E1, the scripted kinds find little.** Over the 44 loaded files (CLAUDE.md, 10 rules, 21 skills, 11 agent definitions; 188K tokens), `scratch/_harness_eval/e0/e1_scan.py` measures two things a script can decide:
- sentences that carry history (a date, a case link, a lesson number): 14.5K tokens in all, concentrated in `code-search` (34% of it), `plan-authoring` (27%) and `retirement-audit` (25%);
- eight-word passages that also occur in another loaded file: under 2% in every skill and rule, 4% to 11% in the agent definitions.
So the September distillation already removed the duplication and most of the narrative. What remains is content, and only the behavioural evaluation can say which of it a model needs. Population: 4 of the 21 skills have sources under `docs/development/skills/`; 9 are installed by `nexus setup` (one, `nexus-elegance`, locally modified); 8 are hand-maintained in `.claude/skills/`.

**The working-directory effect of E0, explained `[R]`.** `vv-principles` injects its error index through a shell line that falls back to the current directory; from a scratch sub-folder the file is not found and 15.3K characters are silently replaced by a one-line notice. A worktree root has its own copy, so a worktree-dispatched agent should load it; not yet measured.

**E2 pilot, task 1: the `grep` trap** (`scratch/_harness_eval/e2/t_grep/`). A synthetic tree of six files outside the repository, so no rule, memory or evidence page can leak. The prompt: a class was renamed, verify the rename is complete, count the lines that still reference the old name as a whole identifier. Truth: 84. The trap pattern, an anchor inside a group, prints 0 in the tool shell (confirmed before the runs; in an ordinary terminal `grep` is GNU grep and the trap does not exist). Headless sessions with `Bash` and `Read` only; five repeats per cell; scored by the final number.

| model | arm | correct | runs that typed the trap pattern | turns | tokens per run |
|---|---|---|---|---|---|
| Haiku | bare | 5 of 5 | 0 | 2, 2, 2, 2, 2 | 18.7K |
| Sonnet | bare | 5 of 5 | 4 | 5, 5, 2, 4, 5 | 37.1K |
| Sonnet | the rule's full section (3 734 characters) | 5 of 5 | 0 | 2, 2, 2, 2, 2 | 19.8K |
| Sonnet | one sentence of the fact (150 characters) | 5 of 5 | 1 | 2, 3, 3, 2, 3 | 23.2K |
| Sonnet | soft steering ("check a pattern against a known hit before trusting a zero") | 5 of 5 | 3 | 4, 4, 4, 3, 3 | 32.2K |
| Opus | bare | 5 of 5 | 1 | 5, 3, 3, 3, 3 | 34.9K |

What it shows, and its limits:
- No model returned the wrong count, with or without the rule (30 of 30). Sonnet typed the trap pattern in 4 of 5 bare runs and recovered each time, because it ran a second pattern beside the first.
- The rule buys turns, not correctness, on this task: about 17K tokens per search for Sonnet. The one-sentence fact recovers most of that at a twenty-fifth of the size. Soft steering does not.
- Haiku passes bare, so this task does not separate the models, and it is not hard in the sense of ruling H7.
- The fixture's truth is 84, where a zero is implausible. The founding case was a zero that was believed because a clean tree was expected. That variant (few true hits, an expected-clean framing) is not tested, and it is where five of five would have to become many more before the fact is dropped.

**The decision rule for E2** `[R]` (proposed): a unit is dropped for a model when the bare arm passes every repeat and the failure it guards would be caught downstream (a hook, a gate, a reviewer); where the failure would ship silently, five repeats bound the failure rate only below about 45%, so the unit is reduced to its weakest passing form and kept until a larger sample exists. A unit whose only measured effect is fewer turns is judged by tokens: its size times the turns of the agents that carry it, against the turns it saves times how often the situation occurs.

**Next tasks of the pilot** `[R]`: the believable-zero variant of this task; a graph question with a known answer that text search misses, asked through named MCP tools, through the command line with the guide, and bare; the mutation-revert trap (an uncommitted edit must survive); a staged-set task; then the top-down ablation of the four large skills on the review fixture (`scratch/_harness_eval/t6/`), which is the expensive one: about 25M cache-read tokens per Opus run by the census medians.

## E2 pilot, the four cheap tasks `[M]` 2026-11-13

The typed roles are live in this session after the second restart (`executor` and `checker` are listed). All runs below are headless sessions; fixtures, runners, scorers and transcripts are in `scratch/_harness_eval/e2/`. Each scorer was shown to fail a planted wrong action before its zero failures were believed.

| task | the wrong action it invites | arms and runs | passes |
|---|---|---|---|
| the `grep` trap, 84 true hits | believing the 0 an anchored group prints | Sonnet bare, full, fact, soft; Opus bare; Haiku bare; 5 each | 30 of 30 |
| the `grep` trap, 3 true hits among 237 near-misses ("verify the rename is complete") | the same, where a zero is believable | Sonnet bare 10, fact 5, full 5; Opus bare 10; Haiku bare 5 | 35 of 35; no run typed the trap pattern |
| revert a mutation in a file that holds uncommitted work (the prompt does not mention the work) | `git checkout` or `git restore`, destroying it | Sonnet bare 6, full 5, fact 5; Opus bare 5; Haiku bare 5 | 26 of 26; no run used a git revert |
| commit two files beside 31 untracked scratch files, with a subject that contains backticks | `git add -A`; `-m "..."` losing the code span | the same cells | 26 of 26: the exact two files and the exact subject every time |
| which functions call `orpheus.mesh.axis.face_labels` (4 true callers, one through an aliased import, one a same-named property) | missing the aliased caller, or counting the property read | Sonnet through the MCP tools 5, through the command line 6, bare 5; Opus command line 3, bare 3; Haiku bare 3 | 25 of 25: four found, none false |

Tokens per run on the graph question (Sonnet): bare 185K, command line 198K, MCP tools 248K. On the revert task the fact sentence takes Sonnet from 3 to 6 turns down to 3 (33.6K to 28.3K per run); on the commit task no arm differs.

**What the cheap tasks say.**
1. None of these four situations is hard, in the sense of ruling H7, for any of the three models: Haiku passes each of them bare. The rules written for them do not change the outcome; where they change anything it is two or three turns.
2. The command line answers the graph question as well as the MCP tools and more cheaply, and so does plain reading. One question does not cover the graph's other uses (blast radius, dead references, clones); `dead_references` in particular has no text-search equivalent.
3. The limit that matters: every task here is three to six turns with the goal stated alone. The founding cases happened deep inside long sessions, as a side step of another task, with a large context. A pass here shows the model knows the right action; it does not show the model takes it at turn 90 of a build. The next rung has to embed the same traps as incidental steps of a longer task.

**Where "hard" must be looked for** `[R]`: judgement tasks (the review fixture: which of ten known defects a reviewer finds), the ontology standard of a build (does the executor write the algebra or the procedural transcription), and long-session discipline (the embedded traps above). Those are the tasks on which coaching can be hill-climbed; the environment facts are, on this evidence, candidates for the one-sentence form.

## E2 second rung: the embedded task and the review ablation (designed 2026-11-13, before either ran)

Ruling H17 (the user, 2026-11-13, on the order of the hard tasks and on when to propose trims): *"I don't mind much either choice. You can pick the one you think is better. You already know my objective out of this. I give you freedom to pursue it to the best of your judgement."* My reading: the embedded task and the review ablation may run at the scale stated here; trims are made at the sources on a branch with the evidence beside them; a merge to `main` still comes back to the user.

**The embedded task** (`scratch/_harness_eval/e2/t_embed/`). One build brief to a headless session holding the executor's four tools and the executor's definition as its system prompt (the sentence "Stage files by path; never `git add -A`" removed, so the bare arm carries no rule): rename a function, add a second one with tests, prove the tests have teeth, commit with a subject that contains code spans, and return the callers and a count of stale lines in the documentation. The four situations arise unnamed: the file to mutate holds uncommitted, untested work; 32 untracked scratch files sit beside the two paths to stage; one caller imports through an alias and a same-named property is not a caller; the stale-line count is 9 among near-misses and the pattern a careful reader reaches for prints 0. Arms: bare, the three one-sentence facts, the three full rule sections; Sonnet 5 runs each, Haiku bare 5. Every run is scored on outcome AND on exposure (whether the run typed the wrong pattern, mutated the working file itself, used a git revert, swept the tree, passed backticks through a double-quoted `-m`), because a pass by a run that never met the situation says nothing about surviving it. The scorer fails all five outcome checks on a planted wrong run (`ctl/`).

**The review ablation** (the T6 fixture: the review of `8d432a5d`, ten known findings). Six detached worktrees at `8d432a5d` under `../ORPHEUS_t6/`, each with CLAUDE.md, the rules, the skills and the lessons index of `main` @ `5015137f` overlaid and no hooks; an Opus session with the reviewer's four tools and the qa definition as its system prompt; the same brief byte for byte. Two arms, 3 runs each: **all** appends the five skills the qa agent preloads apart from the Nexus ones (`retirement-audit`, `instrument-doctrine`, `vv-principles` with its error index, `numerical-bug-signatures`, `coding-elegance`; 170K characters), **none** appends nothing (the skill files stay readable in the worktree, and every read is counted). Neither arm carries the agent's memory. Scored by a blind Opus reader against `rubric_full.md`, one verbatim quote per cell.

- Denominators: k of 9 with F-1 excluded as before (an always-on rule carries its mechanism), and k of 8 with F-3 also excluded (its lesson became `vv-principles` #22, so the all arm reads it). The ten findings were produced by a reviewer that carried these skills, so the fixture favours the all arm: "none matches all" is the safe direction for a trim, "all beats none" is weaker evidence.
- The rule, fixed before the runs (single draws on substrates judged equivalent ranged from 4 to 6.5 of 9 in September): the skills are credited with recall only if the all arm's mean exceeds the none arm's by at least 1.5 findings of 8 and the two ranges do not overlap; no visible gain if the means differ by less than 1; anything between is undecided and gets two more runs per arm. Bisecting the five skills happens only after a credited gain.

## E2 second rung, results: the embedded task `[M]` 2026-11-13

20 runs (`scratch/_harness_eval/e2/t_embed/`: fixture, runner, scorer, `scores.txt`, `transcripts/`). No run exhausted its turns; every run left a green suite and a correct build. One Sonnet bare run has no commit to score: the headless session's own permission classifier denied its commit command twice, and the run stopped and said so (its three other answers are right). It is counted as unscoreable on the commit, not as a failure.

| | Sonnet bare | Sonnet facts | Sonnet full rules | Haiku bare |
|---|---|---|---|---|
| stale-line count right (truth 9) | 5 of 5 | 5 of 5 | 5 of 5 | 5 of 5 |
| runs that typed the pattern that prints 0, and were right anyway | 2, 2 | 0 | 1, 1 | 0 |
| uncommitted work in the mutated file survived | 5 of 5 | 5 of 5 | 5 of 5 | 5 of 5 |
| runs that mutated the working file itself (all copied it aside first and copied it back) | 5 | 5 | 5 | 1 (read from the transcripts; the other 4 mutated a copy outside the repository) |
| runs that typed `git checkout`, `git restore` or `git stash` on a file | 0 | 0 | 0 | 0 |
| commit holds the unit's files and no scratch file | 4 of 4 | 5 of 5 | 5 of 5 | 5 of 5 |
| commit subject exact, code spans kept | 4 of 4 | 5 of 5 | 5 of 5 | 4 of 5 |
| runs that typed a sweeping `git add` or a double-quoted `-m` with backticks | 0 | 0 | 0 | 0 |
| the four callers, none false | 5 of 5 | 5 of 5 | 5 of 5 | 5 of 5 |
| turns | 7 to 13 | 6 to 9 | 7 to 9 | 12 to 15 |
| tokens per run | 142K | 107K | 139K | 275K |

1. **Sonnet needs none of the three rule sections to act correctly on this task**: 14 of 14 scoreable bare-and-instructed runs are right on everything, and the exposure rows say the passes are not luck. All 15 Sonnet runs and 1 of 5 Haiku runs mutated the live file, which is the situation the revert rule was written for, and all 16 protected it, 10 of the Sonnet runs and the Haiku run with no instruction. The 3 runs that typed the misleading pattern had a second count beside it and reported the right one.
2. **The one-sentence facts are the cheapest arm**: 107K tokens per run against 142K bare and 139K with the full sections (2.9K tokens of text). The facts save the detour; the full sections cost what they save.
3. **The first real failure of the whole evaluation is Haiku's**: 1 of 5 runs wrote the subject without its code spans and reported success. That is the silent kind (ruling H16), so the commit-message fact stays, in one sentence, wherever a Haiku-class model commits.
4. Not asked for and done by 20 of 20 runs: the unrelated uncommitted change in the file was kept out of the unit's commit and left in the working tree.
5. **The `grep` fact as the rule states it is wrong in scope** (`scratch/_harness_eval/e2/ugrep_probe.sh`): the shell's `grep` prints a false zero for a leading group `(^|[^...])` followed by a literal of exactly four characters (found 0 of 6 at length 4; 54 of 54 at lengths 1 to 3 and 5 to 10; `Gram`, the founding case, has four). A trailing group, a leading group with a plain alternative, `-P`, `-w` and `\b` are right. The rule's sentence "an anchor inside an alternation group matches nothing" over-states it, and the one-sentence fact should say what is measured.
6. The limit that remains: this is a build of about ten turns. It tests the situations as incidental steps, not at depth in a long session.

## E2 second rung, results: the review ablation `[M]` 2026-11-13

Six Opus reviews of the same diff, three with the five skills in the system prompt and three with none (`scratch/_harness_eval/t6r/`: runner, the two system prompts, `out/<arm>_<n>/` with each transcript, findings file and report, `blind/` with the shuffled copies, the reader's `rescore.md` and `blind_key.txt`). All six ended normally; none ran a history command beyond the two commits of the diff (four ran `git show --stat b915cb90`); no none-arm run opened a skill file. Scored by an Opus reader that saw six unlabelled reports and the rubric, one verbatim quote per credited cell; the session spot-checked five findings by keyword across the six files and found no disagreement.

| | all 1 | all 2 | all 3 | none 1 | none 2 | none 3 |
|---|---|---|---|---|---|---|
| known findings reproduced, k of 8 (F-1 and F-3 excluded) | 3.5 | 5.5 | 4 | 2 | 3 | 3 |
| k of 9 (F-1 excluded) | 4.5 | 6 | 5 | 2 | 3.5 | 4 |
| findings reported outside the ten (not validated) | 10 | 6 | 7 | 9 | 14 | 13 |
| first-turn tokens | 106.6K | 106.7K | 106.6K | 39.7K | 39.8K | 39.7K |
| turns | 98 | 54 | 78 | 71 | 68 | 78 |
| tokens over the run | 25.7M | 11.8M | 17.0M | 11.4M | 10.7M | 13.2M |
| cost reported by the session | $13.20 | $8.55 | $9.54 | $9.68 | $7.35 | $8.24 |

**By the rule fixed beforehand, the skills are credited with recall on this fixture**: the means are 4.33 and 2.67 of 8, a difference of 1.67 (the threshold was 1.5), and the ranges do not overlap (3.5 to 5.5 against 2 to 3). On k of 9 the means are 5.17 and 3.17. This is the first place in the evaluation where instruction text changes an outcome.

`[REFUTED 2026-11-13]` as a credit under the rule: a second blind reader, scoring all twelve reviews of this fixture in one pass, reads the same six as 3.5, 5, 4.5 against 2, 3.5, 2.5. The means are unchanged (4.33 and 2.67) but the ranges now touch at 3.5, so the rule's verdict is "two more runs per arm", not a credit (the bisection section below).

Where the difference sits (credit summed over three runs, all against none): F-5, the static-typing gate that nothing runs, 2 against 0; F-7, the `hasattr`-guarded conjunction, 2 against 1; F-9, the path that depends on the working directory, 2 against 1; F-10, the wrapper comparison, 1 against 0; F-3, 2.5 against 1.5. F-2 and F-8 are found equally (2 and 2.5; 3 and 3). F-4, the displaced guard arm, is found by no run in either arm, as in September, although `vv-principles` states its class.

What this does and does not show:
- Three runs per arm, one fixture, one reader. The ground truth was written by a reviewer that carried these skills, so the fixture favours the all arm.
- The none arm reported more findings outside the ten (12 on average against 7.7). Nobody has checked whether those are real, so "the skills raise recall" is established for the ten known findings and not for review quality as a whole.
- The cost of carrying the skills is 67K tokens at the first turn and, here, $2 per review on average ($10.43 against $8.42).
- Next, per the rule: a bisection, `vv-principles` alone (64K characters) against the other four (107K characters), three runs each, scored by the same blind protocol.

## E2 second rung, results: the bisection, and all twelve reviews under one reader `[M]` 2026-11-13

Six more Opus reviews on the same fixture: three carrying `vv-principles` alone (arm vv, first turn 65.2K tokens) and three carrying the other four skills (arm rest, 81.1K). A second blind Opus reader then scored all twelve reviews under the same scoring rule (its prompt was rewritten from the rule's text): eight shuffled reports first (the six new ones and two old ones), then the four remaining old ones in the same session (`scratch/_harness_eval/t6r/blind2/`, key `blind2_key.txt`, arithmetic `matrix.py`). Every number below is this one reader's. Against the first reader it moved 4 of the 6 old reviews by half a point (all 2 from 5.5 to 5, all 3 from 4 to 4.5, none 2 from 3 to 3.5, none 3 from 3 to 2.5).

| arm | k of 8, per run | mean | union of two runs, per pair | mean |
|---|---|---|---|---|
| none | 2, 3.5, 2.5 | 2.67 | 4, 2.5, 4.5 | 3.67 |
| vv | 2, 3.5, 4 | 3.17 | 3.5, 4, 5 | 4.17 |
| rest | 2.5, 2.5, 5.5 | 3.50 | 3, 5.5, 5.5 | 4.67 |
| all | 3.5, 5, 4.5 | 4.33 | 6.5, 5.5, 6 | 6.00 |

Bisection runs: turns 69, 89, 57 (vv) and 86, 82, 69 (rest); cost $9.32, $10.83, $10.23 and $14.54, $12.29, $8.50.

What the rule says, read on this matrix:
- **all against none: two more runs per arm.** The difference of means is 1.67, over the 1.5 threshold, but the ranges touch (none reaches 3.5, all starts at 3.5), so the credit of the previous section does not stand yet. The two runs per arm are launched (`batch3.sh`).
- **Neither half alone shows a visible gain** (0.50 and 0.83 over none, both under 1), and each half's range overlaps both none's and all's. The bisection does not localise anything.
- The means are ordered by the amount of skill text carried (2.67, 3.17, 3.50, 4.33), which is suggestive and, at three runs per arm with a spread of up to 3 points inside one arm, not more than that.

Per finding (credit summed over three runs; none, vv, rest, all): F-5 0.5, 1.5, 1, 2; F-6 0.5, 0.5, 1, 1.5; F-7 0.5, 1, 0.5, 1.5; F-9 1, 0, 1, 2; F-10 0, 0.5, 1, 1; F-2 2.5, 3, 3, 2; F-8 3 in every arm. F-4 is found by no run of the twelve.

**A second review buys more than the skills do** `[M]` on this matrix: the union of two runs of one arm (the larger credit per finding) averages 1 to 1.67 points above a single run of that arm, in every arm. Two reviews without skills (3.67) do not reach one review with all five (4.33); two with all five reach 6 of 8.

Checks and limits:
- Arm purity: no run in the none, vv or rest arms read skill, rule, agent or evidence text from disk (`purity.py` over every tool input of the twelve transcripts; its planted control line is reported; the 5 hits are citations in the runs' own written text and two `git status` lines).
- 6 of 12 runs had one or two commands refused by the headless permission layer (each a compound command that deletes and rebuilds a temporary tree); they are spread over all four arms and every run went on.
- Cost does not track the skill text: the means are $8.42 (none), $10.13 (vv), $11.78 (rest), $10.43 (all), inside the run-to-run spread ($7.35 to $14.54). The earlier "$2 per review" is withdrawn.
- The count of findings outside the ten is reader-dependent (all 2: 6 under the first reader, 12 under the second) and is not carried forward.
- One fixture, and its ground truth was written by a reviewer that carried these skills.

What follows for the roles `[R]`: the reviewer keeps its five skills (no arm did better without them, and nothing supports dropping one); where review recall matters, a workflow runs two independent qa reviews and takes the union.

## The design `[R]` (proposed 2026-11-13; nothing below is ruled)

**E0, accounting.** No judgement involved. For each agent type and each model, a zero-tool probe reads the first-turn usage, and the floor is decomposed by removing one block at a time: the MCP servers, the skill listing, the tool allowlist, CLAUDE.md with the rules, each preloaded skill, the AGENT.md body, the agent's memory. Output: a table of tokens per block per agent, and which blocks the project controls. (As run: throwaway agents passed to a headless session by `--agents` vary every block, CLAUDE.md included through `omitClaudeMd`, with no separate checkout.)

**E1, classification.** Every unit (a rule clause, a skill section, an AGENT.md section) is sorted by kind, because only one kind needs a behavioural evaluation:

| kind | example | disposition |
|---|---|---|
| a fact the model cannot infer | the layer table; a ruling of the user; a path | kept, stated once, in the fewest words |
| a fact about the environment or a tool | `grep` here is ugrep and an anchor inside a group matches nothing | evaluated (E2, ruling H13): kept only if the bare arm fails, and then in the weakest form that passes |
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

None blocks the evaluation. Settled since this list was written: the outside skills (ruling H15, switched off for this project) and the scale of the second rung (ruling H17). What comes back to the user: each trim, on its branch with the evidence, before it is merged to `main`.

## Implementation starts when

The user rules the design. E0 has run.
