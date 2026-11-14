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
- **all against none: two more runs per arm (run; see the subsection below).** The difference of means is 1.67, over the 1.5 threshold, but the ranges touch (none reaches 3.5, all starts at 3.5), so the credit of the previous section does not stand yet. The two runs per arm are launched (`batch3.sh`).
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

### The two more runs per arm `[M]` 2026-11-13

Four more reviews (all 4, all 5, none 4, none 5), scored by the same reader in the same session (S13 to S16), with the same purity and history checks (clean).

| arm | k of 8, five runs | mean | union of two runs, mean over the 10 pairs | cost, mean |
|---|---|---|---|---|
| none | 2, 3.5, 2.5, 2, 3 | 2.60 | 3.40 | $7.52 |
| all | 3.5, 5, 4.5, 6, 4 | 4.60 | 5.85 | $10.31 |

- The difference of means is 2.00 of 8. In 24 of the 25 all-against-none pairs the all run scores higher, and the 25th is a tie at 3.5 (all 1 and none 2). The rule's letter (disjoint ranges) is missed by that one tie; the rule named no terminal test after the extra runs, so the verdict here is the session's judgement and is stated as such: **the five skills together raise the reviewer's recall of the known findings on this fixture by about 2 of 8.**
- Where it sits (credit summed over five runs; all, none): F-5 4, 0.5; F-7 2.5, 0.5; F-6 2.5, 1; F-10 2, 0; F-9 3, 2; F-3 2.5, 2; F-2 4, 4; F-8 5, 5; F-4 0, 0.
- One review with the skills (4.60) beats the union of two without (3.40); two with the skills reach 5.85.
- Still one fixture, with ground truth written by a reviewer that carried these skills; which skill carries the gain is not known (the bisection above).

What follows for the roles `[R]`: the reviewer keeps its five skills, as a set; where review recall matters, a workflow runs two independent qa reviews and takes the union (about 1.25 more of 8 for a second review's cost).

## What the evidence supports so far `[R]` (2026-11-13; proposals, each comes back to the user before it reaches `main`)

| unit | evidence | proposal |
|---|---|---|
| the `grep` trap (`code-search`) | Sonnet right in 20 of 20 embedded builds and 3 of 3 that typed the trap; the bug's scope is narrower than the rule said | keep, one sentence, with the corrected scope `[LANDED 29cdb08e]` (merged to `main` 2026-11-14 on the user's go-ahead); a false zero ships silently |
| never revert a mutation with git (`process-discipline`) | 0 of 20 runs reverted with git; 16 mutated the live file and all copied it aside first | keep, one sentence; the loss is silent |
| commit message with backticks by `-F` | Sonnet 14 of 14 right bare; Haiku lost the code spans in 1 of 5 and reported success | keep, one sentence, in any brief a Haiku-class agent commits under |
| never `git add -A` | 0 of 20 runs; the hook refuses it and prints the reason | drop the prose from stage-role bodies (done in the embedded task's arms); the rule keeps only what the refusal does not say |
| the reviewer's five skills | 4.60 against 2.60 of 8 over five runs per arm; no single skill isolated | keep all five on qa |
| a second independent review | the union of two runs adds 0.80 (none) to 1.25 (all) of 8 | where recall matters, the template runs two qa reviews |
| CLAUDE.md and the rules on a stage role (H13) | the embedded builds ran with no project rule loaded and were green 20 of 20, on a toy package | executor and checker stay without them for mechanical stages; NOT yet shown for a build held to the ontology standard |

Assignment, as far as measured: Sonnet for build, amend, measure and draft stages (bare body plus three fact sentences); Haiku for fixed-schema categorising, not for committing; Opus for review, with the five skills.

**Not yet measured, and the next rung:** whether a Sonnet executor writes the algebra and not the procedural transcription (Cardinal Rule 1's warning), with and without `coding-elegance`, on a real ORPHEUS unit at a pinned commit, scored blind by an Opus elegance reader. This is the task that says what "hard" means for Sonnet, and it decides the build stages' preload.

## E2 third rung: a build held to the ontology standard (designed 2026-11-14, before any run)

Ruling H18 (the user, 2026-11-14): "Yes, merge and go ahead with the Sonnet evaluation."

**The fixtures are the project's own record.** On 2026-10-23 the `build-units` workflow built wave 2 unit by unit: a Sonnet stage wrote the gates, a Sonnet stage built (`[M]` the template at those commits names `model: 'sonnet'` for both), Opus qa and elegance reviewed, and the amend commit's body lists each finding with its disposition. So each unit carries its literal brief (in `wave2_boundary_carve.md` at the gates commit), its red gates, one historical Sonnet build made under that day's full harness, the findings a review raised on it, and the amended build.

| fixture | kind | base (gates red) | historical build | amended |
|---|---|---|---|---|
| L4, `ordinate_pairing` and the `TraceRole` tier | additive, numerics | `2f1d04e2` | `21a5e0e7` | `afc8ed6d` |
| S3, the diffusion realizer table and binding | carve, diffusion | `232fccf5` | `9c0e8fde` | `508c6478` |

`[M]` scorer anchors on L4 (`scratch/_harness_eval/sb/anchor.sh`): the original gates are 23 of 23 red at the base and 23 of 23 green on the historical build; the amended gate file reads 1 red of 25 on the historical build and 25 green on the amended one.

**Rungs** are how much of the ontology the brief hands over, because that is a larger variable than any skill: (1) the literal brief with the plan sections and the red gates, as the workflow gives it; (2) the plan sections without the gates; (3) the outcome alone. Rung 1 runs first; a later rung runs only where the earlier one passes.

**Arms.** Bare first: Sonnet with the executor body and the three environment facts, no CLAUDE.md, no rule, no skill; five runs per fixture. An arm is added only for a finding that recurs: `coding-elegance` in the system prompt (Sonnet), and Opus bare.

**Isolation, each with a control.** The run directory is a `git archive` export of the base with a fresh one-commit history (no later commit is reachable). `.claude/` is removed except the two plan files the brief cites. The interpreter is a separate virtual environment that sees the project's dependencies and not its editable install, so `orpheus` resolves only to the export; control: the original gates are red in the export under a plain `python -m pytest`. Every transcript is scanned for a read under the main checkout.

**Outcomes, scored by the session and never taken from the run:** the original gates green; the files touched inside the brief's list; the unit's neighbouring tests and the layer gate green; the amended gate file's reds; and, per fixture, a rubric of the historical findings (written from the amend commit's body and diff before any run), each judged present or absent in the run's diff by one blind Opus reader in one session, with the historical build and the amended build among the unlabelled diffs as anchors (the reader must find the findings in the first and not in the second, or its cells are not used).

**The rule, fixed now.** Bare passes rung 1 on a fixture when at least 4 of 5 runs have the original gates green inside the file list AND the mean count of rubric findings present is no higher than the historical build's. A finding present in at least 3 of 5 bare runs is recorded as hard for Sonnet at this rung; for those, each added arm gets five runs, and it is credited with the finding when the rate falls to at most 1 of 5. No runs beyond that: a finding no arm removes stays the reviewer's, which is where the workflow already puts it.

## E2 third rung, results: bare Sonnet at rung 1 `[M]` 2026-11-14

Five bare Sonnet builds per fixture (`scratch/_harness_eval/sb/`: `mkrun.sh`, `run.sh`, `score.py`, the briefs, the rubrics, `runs/<fixture>_<model>_<arm>_<n>/` with the transcript, the build's diff and `score.json`), each from the literal brief with the gates already committed. The scorer ran every test itself. Anchors: the historical build and the amended build of each unit, exported and scored the same way.

Changes to the design, made before the batch: the second fixture is unit L1 (the deck family and the colouring symmetry, geometry; base `34c9c9fd`, historical build `487d2d19`, amended `8a3e6c16`), not slice S3, whose build spans three commits; the export keeps every tracked plan file of the base, and for L1 the one untracked memo its brief cites (dated 2026-10-22).

| | L4, 5 bare runs | L4 historical build | L1, 5 bare runs | L1 historical build |
|---|---|---|---|---|
| original gates green | 5 of 5 | yes | 5 of 5 | yes |
| production files outside the brief's list | 0 of 5 | 1 file (`tangential.py`) | 0 of 5 | none |
| amended gates red | the same 1 in 5 of 5 | that 1 | the same 4 in 5 of 5 | those 4 |
| layer and neighbouring gates | the package-cycle gate red in 4 of 5 | green | green in 5 of 5 | green |
| cost per build | $0.28 to $0.41 | | $0.54 to $0.73 | |

- **Bare Sonnet reproduces the historical build's defects exactly as far as the amended gates see them**: every finding of the 2026-10-23 review that became a gate row is red in 10 of 10 bare builds, and it was red in the build made under that day's full harness. On these findings the harness text bought the Sonnet builder nothing, and its absence cost nothing.
- **What "hard" means for Sonnet at this rung is everything the gates do not demand.** The five gated findings are: an error that does not name its real cause (the folded rule), a value whose constructor does not state its invariants, a raw input not parsed before use, a law checked on part of its statement (the pairing's linear part only), and one convention with two addresses. The builds turn the given gates green and stop there.
- **One thing the bare builder does worse:** the L4 builds put the new function where the brief said and thereby closed an import cycle inside `numerics` (the cycle gate is red in 4 of 5, and in the untranscribed pilot). The historical builder avoided it by moving a constant to a new module outside its file list. The brief names neither the cycle gate nor the wider suite, and the 2026-10-23 builder had CLAUDE.md, which does. This is a fact a stage brief must carry: which gates beyond the unit's own to run before committing.
- All ten runs: no read under the main checkout, no denial, one or two commits, nothing uncommitted. First turn 10.5K to 10.7K tokens.

By the rule, bare does not pass rung 1 on L4 (the cycle) and matches the historical build on L1; all five gated findings recur in 5 of 5, so both added arms run on both fixtures (`batch2.sh`): Sonnet with `coding-elegance` in the system prompt, and bare Opus. The blind reader's pass over the full rubrics covers every arm in one session, after those runs.

## E2 third rung, results: the two added arms and the blind rubric pass `[M]` 2026-11-14

Twenty more builds (Sonnet with `coding-elegance` in the system prompt, first turn 22.9K tokens against 10.5K bare; bare Opus), five per fixture per arm. All 30 builds of the three arms: original gates green, inside the brief's file list, no read under the main checkout, no denial. One blind Opus reader per fixture judged every build's diff against the fixture's rubric, with the historical and the amended build among them unlabelled (`scratch/_harness_eval/sb/blind_<fixture>/scores.md`, keys `key_<fixture>.txt`, arithmetic `tab.py`).

**The gated findings (deterministic; builds in which the amended gate row is red, of 5).**

| | Sonnet bare | Sonnet + `coding-elegance` | Opus bare |
|---|---|---|---|
| L4, the folded rule's error | 5 | 5 | 5 |
| L4, the package-cycle gate | 4 | 5 | 4 |
| L1, the pairing's whole law | 5 | 5 | 5 |
| L1, one agreement convention | 5 | 5 | 5 |
| L1, the raw partition parsed | 5 | 5 | 0 |
| L1, values state their invariants | 5 | 5 | 3 |
| cost per build | $0.28 to $0.73 | $0.29 to $0.82 | $1.96 to $6.17 |

**The rubric (blind reader; findings of the 2026-10-23 review present, mean per build).** The anchors read as they must on L4 (historical 9 of 9 present; amended 0 of 8, and the ninth item is the extra file the review accepted). On L1 the amended build reads 0 of 11 and the historical build 8 of 11: items B2, B7 and B9 do not separate the anchors (B2 was posed against a test the builds do not contain) and are dropped, leaving 8.

| | historical build (full harness of that day, one draw) | Sonnet bare | Sonnet + `coding-elegance` | Opus bare |
|---|---|---|---|---|
| L4, of 9 | 9 | 4.67 (6 builds, the pilot included) | 4.80 | 5.00 |
| L1, of the 8 valid | 8 | 7.80 | 8.00 | 4.60 |

- **`coding-elegance` in the builder's prompt changes nothing** on either fixture, gated or read (4.80 against 4.67; 8.00 against 7.80). By the rule it is credited with no finding.
- **The bare Sonnet build is no worse than the historical one**, and on L4 better (4.67 against 9 of 9: no bare build has the redundant guard, the unread read surface or the extra file). The harness text of 2026-10-23 did not help its builder on these findings.
- **Opus as the builder is credited with three L1 findings** (present in 0 of 5: values that state their invariants, the raw input parsed, the sign of the scale under `np.spacing`), and with none on L4. These are completeness of a single function or type. It does not remove the findings about structure: the law stated in part (B1), the union re-spelled by hand (B3), the unrecorded twin (B4, A7), a convention with two addresses (B8), an illegal role left spellable (A4), an undeclared scope boundary (A5). Each is present in every build of every arm, except B8, absent in 1 of the 5 Opus builds.
- The reader and the gates disagree in one place and it is the gate's shape: on L1 "values state their invariants" is red in 3 of 5 Opus builds and read absent in 5 of 5 (the gate asks for a specific refusal); on L4 the folded rule is red in 5 of 5 Opus builds and read present in 3.

**What follows `[R]`.**
- The build stage runs on Sonnet with the bare executor body. No skill preload is supported for it, and Opus buys the removal of the completeness findings at seven to ten times the cost; a reviewer and an amend stage already exist for those.
- The structural findings are not reached by any builder arm at this rung. They are reached by the review (which is why the reviewer's skills and a second review matter), or by a gate: every one of them that was written as a gate row after the review would have been demanded of the builder had the gates stage written it. The leverage is upstream of the build, in the specification and the gates, which is the next thing to measure: the gates stage, by model and by preload, scored on whether its gates redden on the historical build's known defects.
- A stage brief names the gates beyond the unit's own that the builder runs before committing (the package-cycle gate was red in 13 of 15 L4 builds and no brief named it).
- Rungs 2 and 3 (less of the ontology handed over) were not run: rung 1 already leaves the structural findings in every build, so a thinner brief can only add to them.

## E2 fourth rung: the gates stage (designed 2026-11-14, before any run)

Ruling H19 (the user, 2026-11-14): "Yes. Go ahead with the gates-stage measurement."

**A fact found while designing, which shapes the rung.** `[M]` the plan the gates stage received (`wave2_boundary_carve.md` at `34c9c9fd`, sections 5.1 and 5.4) enumerates each unit's done-when gates one by one, with the inputs and the mutation that reddens each (ten for L1, four for L4). The historical gates stage transcribed a list; the defects the review later found are the ones that list did not name. So the stage is measured at two rungs:
- **G1, the literal brief** (the list given): does a gates writer add anything to the list that reddens on a known defect?
- **G2, the list withheld**: the export's plan has the "Done-when gates" bullets of the unit's section removed and the brief says to write the gates the unit's types and laws require. This measures the writer's own gate design.

**Fixtures:** the same two units, from the commit before the gates (L1 `d9d638e6`, L4 `f6690aeb`), exported and isolated as in the third rung. The writer commits tests only.

**Scoring, deterministic.** The run's test files are laid over three trees and run by the scorer: the base (every gate must be red or erroring), the historical build, and the amended build. A test that fails on the historical build and passes on the amended one is a **catch**; a test that fails on the amended build is a **wrong gate** (it demands something the specification does not, or is broken). Catches are mapped to the rubric's findings by reading the failing assertion. Anchors: the historical gate files must read 0 catches, the amended gate files 1 (L4) and 4 (L1).

**Arms and rule.** Sonnet bare at G1, five runs per fixture. The historical gates stage is the reference: 0 catches. An arm is credited at a rung when at least 4 of 5 runs have at least one catch and no more than one wrong gate. G2 runs for Sonnet bare, Sonnet with `vv-principles` and `instrument-doctrine` in the system prompt, and bare Opus, five runs per fixture each; G1 gets the added arms only if Sonnet bare is not credited there and G2 shows an arm that is.

## E2 fourth rung, results: the gates stage `[M]` 2026-11-14

Forty runs, all ended `success`: 0 of 40 touched a production file, left uncommitted work, read the main checkout or the removed harness files, or met a permission denial. First turn: about 10.5K tokens bare, 28.8K with the two skills. Scorer: `gscore.py` (catch: fails on the historical build, passes on the amended one); `gcat.py` splits the gates that fail on the amended build by their outcome on the historical one. Archive: `scratch/_harness_eval/sb/` (`gruns/<run>/{run.jsonl,gates.diff,gscore.json,gcat.json}`).

| cell | L4: runs with a catch | L1: runs with a catch | L4: runs with at most 1 gate failing on amended | L1: same | mean cost L4 / L1 |
|---|---|---|---|---|---|
| G1 (list given), Sonnet bare | 0 of 5 | 0 of 5 | 5 of 5 | 2 of 5 | $0.39 / $0.64 |
| G2 (list withheld), Sonnet bare | 0 of 5 | 3 of 4 measurable (run 5 does not collect on the amended build, so it could register no catch) | 5 of 5 | 0 of 4 measurable | $0.40 / $1.30 |
| G2, Sonnet with `vv-principles` and `instrument-doctrine` | 0 of 5 | 4 of 5 | 2 of 5 | 0 of 5 | $0.57 / $1.53 |
| G2, Opus bare | 0 of 5 | 2 of 5 | 4 of 5 | 4 of 5 | $3.06 / $5.46 |

**By the rule written before the runs, no arm is credited at either rung on either unit.** The rule is read per unit. L4 has 0 catches in 20 runs. On L1 only one run of 19 measurable (Opus, run 2) has both a catch and at most one gate failing on the amended build.

What the L1 catches are, by the finding they redden on (runs of 5): the raw partition is parsed (a malformed partition is refused): Sonnet bare 2, Sonnet with skills 4, Opus 1; the values state their invariants (the identity is first): Opus 2, both Sonnet arms 0; a self-paired face refuses a translation: Sonnet bare 1. The amended gate file catches four findings; no run of 20 reddens on "one agreement convention" or on the pairing law and its adjoint.

What this supports:

1. **With the gate list in the brief the gates writer adds nothing that catches a known defect** (0 of 10 runs), the same reading as the historical gates stage (0 catches on both units). The stage transcribes the list.
2. **Without the list, the writer reaches some of L1's findings (9 of 14 measurable runs) and none of L4's (0 of 15).** The L4 writers did test folded rules: 15 of 15 runs hold between 2 and 7 gates that name one (a text match on `fold|spent|orbit|Quotient|gauss_legendre|octant` per test body), and in 13 of 15 runs every one of them passes on both builds; in the other two they fail on both for one import cause each. They assert the plan's diagnosis (the pairing "does not cross"), which the historical build satisfies; the amended diagnosis (the rule has spent the mirror) came from the review. The L1 catches all sit where a plain refusal suffices (`ValueError` on a malformed input). `[R]` A catch needs the gate to pass on the amended build, so a finding whose fix is a new name or a new diagnosis is harder to register than one whose fix is a refusal; the L4 zero is a zero of this instrument, supported by the reading above and not by the count alone.
3. **The two verification skills are not credited for the gates writer.** 4 of 5 against 3 of 5 on L1 is inside what five runs can separate. (That arm's two L4 runs with 22 and 31 gates failing on both builds each have one cause, a module accessor that resolves to the re-exported function, `[M]` the first traceback of each; they are not counted against it.)
4. **Opus is not better at this stage** (2 of 5 on L1, 0 of 5 on L4) at about four to eight times the cost. Its gates are the closest to what either build spells (4 of 5 runs with at most one failing gate on both units).
5. **Gates fix the plan's decisions in place, including the ones the review later overturned.** Gates that pass on the historical build and fail on the amended one: all 5 Opus L4 runs assert that `TraceRole` names four tiers (the plan's decision 49; the amend replaced it with a separate type), and 14 of 20 L1 runs hold such a gate (the agreement convention, the self-pairing invariant, the digest rows). A gate written from the plan cannot catch an error in the plan.

`[R]` A limit of the second criterion, stated after the fact: "does not pass on the amended build" mixes three things, a gate that guessed a spelling neither build used (the package-export gate fails on both builds in 13 of 15 L1 G2 runs), a gate that pins a decision the amend changed, and a wrong gate. The instrument cannot tell the first from the third, because neither build was written against these gates. With the second criterion set aside, one cell passes the first: Sonnet with the two skills on L1 (4 of 5), beside 3 of 4 for bare Sonnet, which is no separation. A scorer limit: the outcome parser cuts a parametrised test id at its first space, so node counts are approximate; run-level verdicts do not rest on them.

Assignment this supports `[R]`: no added skill and no stronger model is credited at the gates stage, and the gate list belongs in the specification. Whether the stage can run bare is not shown (the section on the trims: 2 of 5 bare L1 runs wrote float-rule gates with a wrong expected value). The leverage is the specification (which gates it names, and whether its decisions are right) and the review.

## E3 and E4: the trims drafted and the assignment table `[R]` (2026-11-14; each branch waits for the user before `main`)

Two branches off `main` (`29cdb08e`), one unit each, not merged and not pushed. Both edit `build-units.js` on different lines, so the second to land is rebased. A first draft of each reached further than its evidence and was withdrawn after review: `omitClaudeMd` on the executor itself (it is also the drafter of `spec-check` and `design-panel` and the amender, none measured bare, and those two templates carry no rules block), and the gates stage as a bare executor (see below).

| branch | commit | what changes | evidence | check run |
|---|---|---|---|---|
| `chore/builder-agent` | `3ba92923` | a new stage agent `builder` with `omitClaudeMd`, used by the build stage of `build-units` only; its body keeps two environment facts whose loss is silent (the `grep` false zero, the git revert); the executor keeps CLAUDE.md and the rules and drops the `git add -A` prose; the texts that count "the three `omitClaudeMd` agents" say "Support agents"; a History paragraph on the workflows page | the build rung: 10 of 10 bare Sonnet builds green and inside the file list, no more review findings than the full-harness build | `python -m tools.harness --check`: 0 problems; `test_harness_generated.py` and `test_docstring_xrefs.py`: 48 passed; a probe: a dispatched builder reported no CLAUDE.md and 0 rules in context, a dispatched executor CLAUDE.md and 7 rules. The header is read when the harness starts, so it reaches a dispatch after a restart. |
| `chore/build-units-qa-passes` | `3a81ebb7` | a new argument `qa_passes` (default 1) | the review rung (a second review adds 0.80 to 1.25 of 8) | the template parses; a stub run dispatches two qa passes with `qa_passes: 2` and one by default |

Not changed, with the reason:
- The gates stage stays on the test-architect. `[M]` 2026-11-14: in 2 of 5 bare L1 runs with the list given, gates on the float rule assert an expected value that both builds contradict (run 1: five rows expect `(0,.1,.2,.30000001)` and its like refused as asymmetric; run 5: nine rows expect `(0,.1,.2,.3)` and its like symmetric), which is a wrong gate and not a spelling. The historical test-architect gates are one draw per unit with none, so the comparison is 2 of 5 against 0 of 1 and decides nothing; the bare writer adds no catch and is not shown to be as safe.
- The reviewer's five skills stay (credited as a set). The checker and the executor keep CLAUDE.md and the rules (their stages were not measured bare).
- The tree-wide ledger line in the template's rules block stays: the literal briefs carried none, and 13 of 15 L4 builds left `test_every_package_module_graph_is_acyclic_but_the_listed_cycles` (`tests/gates/test_layer_imports.py`, present at the base and named in that line) red. The trim "the brief template names the wider gates" is therefore already met by the template.
- The `-F` commit fact is not in the builder's body: Sonnet was right 14 of 14 bare, and the template's rules block states it.

The assignment for the workflow stages, as measured:

| stage | model | agent | loads | basis |
|---|---|---|---|---|
| gates | Sonnet | test-architect (unchanged) | as today | the gates rung: no arm credited; the bare writer not shown safe |
| build | Sonnet, a ruling for the user (below) | builder | the brief and the template's rules block | the build rung |
| amend, measure, draft | Sonnet | executor (unchanged) | as today | not measured bare |
| review (qa) | Opus | qa | its five skills; two independent passes where recall matters | the review rung |
| review (elegance) | Opus | elegance-enforcer | unchanged | not measured apart from qa |
| verify, structural check | Opus | checker | unchanged | not measured |
| fixed-schema categorising | Haiku | `general-purpose` | a schema; never commits | the first rung |
| specification | not measured | | | the gates and build rungs put the leverage here: the gate list and the plan's decisions are fixed upstream of every stage measured |

What "hard" means, as far as measured `[M]`: at the build stage Opus removes completeness findings that Sonnet leaves on the second unit (blind means 4.60 against 7.80 of 8 findings present; the raw partition left unparsed in 0 of 5 Opus builds against 5 of 5 Sonnet builds) and none on the first (5.00 against 4.67 of 9), at about five to ten times the cost per build ($1.96 to $6.17 against $0.28 to $0.82). The structural findings survive every builder arm and every gates arm, and were reached only by Opus review with its skills. So Sonnet for the build is a trade: the review and amend stages then carry the completeness findings Opus would not have left. The specification stage and the Fable rung of ruling H7 are not measured.

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
