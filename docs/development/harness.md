# The agent harness — what loads when, and how to change it

ORPHEUS is developed by a user steering Claude Code agents. The knowledge those
agents work from lives in this section of the docs; the `.claude/` directory is
one harness's **generated** view of it. This page records what loads into an
agent's context, what it costs, and the procedure for adding a rule, a skill or
an agent. Measured numbers are marked `[M]` (the 2026-09 evaluation,
`.claude/plans/harness_context_budget.md`, Parts I, VII and VIII, and the
2026-11 evaluation, `.claude/plans/harness_eval.md`); a budget is a budget, not
a measurement.

## What loads, and what it costs

Measured 2026-11-13 and 2026-11-16 by the harness evaluation
(`.claude/plans/harness_eval.md`, closed 2026-11-16; scripts and data in
`scratch/_harness_eval/e0/`). The instrument is the first-turn usage of a
zero-tool probe, read from its transcript as input + cache-creation +
cache-read tokens; three repeats of one cell differ by 14 tokens, and Sonnet
and Opus differ by 6 tokens per agent type, so the model does not change what
an agent holds.

**The first turn of each agent** (Sonnet, dispatched with the `Agent` tool
from the repository root, `[M]` 2026-11-13 after the typed stage roles landed):

| agent | first turn | what dominates it |
|---|---|---|
| qa | 117.8K | its preloaded skills (69.3K) |
| test-architect | 116.3K | its preloaded skills (65.8K) |
| archivist | 112.7K | its preloaded skills (62.6K) |
| numerics-investigator | 111.6K | its preloaded skills (61.9K) |
| method-implementer | 106.4K | its preloaded skills (61.6K) |
| elegance-enforcer | 77.3K | its preloaded skills (26.6K), CLAUDE.md and the rules |
| literature-researcher | 37.5K | its definition and own memory |
| explorer | 36.4K | nine Nexus tools, its definition and own memory |
| executor (the build, amend, measure and draft stages) | 36.3K | CLAUDE.md and the rules |
| checker (the verify and structural-check stages) | 35.9K | CLAUDE.md and the rules |
| cross-domain-attacker | 29.7K | its definition and own memory |
| an untyped agent (`general-purpose`, or a workflow stage that names no agent), by a workflow dispatch | 122.6K | tool schemas (77K), none of them the project's text |

The `builder` (the build stage, `omitClaudeMd`) was not measured as a
dispatch; it is the executor's tool list without CLAUDE.md, the rules and the
memory index. An agent dispatched into a git worktree starts within 1%
of the same agent at the repository root (`[M]` 2026-11-16, test-architect,
one probe per arm, totals only: 115 841 against 116 866 tokens); an agent dispatched while the session's
working directory is a folder with no `.claude/` loses its own memory and
part of its preload without a message (`[M]` 2026-11-13: 120.7K against
130.2K).

**What one block costs** (tokens, `[M]` 2026-11-13, each a difference between
two probes that differ by that block alone):

| block | tokens |
|---|---|
| an agent holding one tool (`Read`), no CLAUDE.md | 5 271 |
| four more built-in tools (`Write`, `Edit`, `Bash`, `SendMessage`) | +4 140 |
| the Nexus server's 45 tools (`mcpServers: nexus`; a Nexus tool named in `tools:` is held alone) | +20 071 |
| the browser server's 25 tools | +6 989 |
| CLAUDE.md (the project's and the user's), the seven always-on rules and the memory index | +26 528 |
| `vv-principles` preloaded | +25 590 |
| `coding-elegance` preloaded | +12 856 |
| agent memory enabled, on an agent with no memory files | +4 051 |
| every tool a headless session holds (87) against five built-ins | +47 815 |

The orchestrator's own first turn is 103.5K in a headless session (Sonnet
103 471, Opus 103 590, Fable 104 899), and about 39K more in an interactive
one, which holds more tools.

**What the always-on block is.** CLAUDE.md, the rules in `.claude/rules/` with
no `paths:` front matter, and the auto-memory index are inherited by every
dispatch that does not set `omitClaudeMd`; `omitClaudeMd: true` drops all
three together and nothing else. Every always-on file is generated or
installed (`tests/gates/test_harness_generated.py` fails on a rule that is
neither): one rule, the Nexus routing rule `nexus-tools`, is installed by
`nexus setup` from the extension and recorded in
`.claude/nexus-install-manifest.json`. Four rules are path-scoped and load
when a matching file is read, in the main agent and in a sub-agent alike:
`coding-standards` under `orpheus/**`, `tests/**`, `derivations/**` and
`tools/**`; `plan-authoring` under `.claude/plans/**`; `vv-testing` under
`tests/**` (ruled 2026-09-21); `documentation` under `docs/theory/**` and
`docs/architecture/**` (ruled 2026-09-22). The session-start batch
(`.claude/hooks/session-start.txt`: three skill cores, the lessons index, the
Nexus briefing) is paid by the main agent only.

**The budget and the measurement.** Each source page declares
`harness: {kind, budget_tokens}`, and `python -m tools.harness --check` holds
each generated file to it by an estimate, characters / 2.6
(`tools/harness/budget.py`). The rate is measured: a 49 582 character file
moved a Sonnet first turn by 18 671 tokens on 2026-11-13, and tool schemas
read 2.6 to 2.9. Until 2026-11-16 the rate was 3.6 (3.51 to 3.53 on a Haiku
probe, 2026-09-20), which read 28% low on the current models; on that day
the rate was changed and every `budget_tokens` multiplied by 3.6 / 2.6 and
rounded up to 50, so each file's length limit stayed where it was and the
printed numbers became tokens. `--check` prints the always-on block (CLAUDE.md,
six generated rules and the installed Nexus rule) at ≈21.9K; the probe
reads 26 528 for that block with the user's own CLAUDE.md (1.1K) and the
memory index (4.6K). An estimate is still not a measurement: a budget bounds
a file's length.

The 2026-09 target for the always-on block, at most 16K by the keep-minus-omit
probe, was met on a Haiku probe at the old rate (15 571 on 2026-09-21; 70 966
before that restructure). At the measured rate the same block is 26.5K, so
the target is not met and is open for a new ruling (issue #714).

## The one-way flow

```
docs/development/{rules,skills,agents}/*.md + lessons.md + onboarding.md     (SOURCE, MyST; each
        │                                        declares `harness: {kind, budget_tokens}` in its front matter)
        │  tools/harness/  (`python -m tools.harness`, a `_GENERATORS` row in docs/conf.py; `--check`, read by
        │  `tests/gates/test_harness_generated.py`: drift, budget, dead links, orphans. `targets/claude_code.py` is the
        │  one harness; the other modules name none, and adding a harness is one module under `targets/`)
        ▼
.claude/rules/*.md   .claude/skills/*/SKILL.md   .claude/lessons.md   AGENT.md definitions   CLAUDE.md   (GENERATED, committed)
```

- Evidence pages (`evidence/*.md`) have no `.claude/` copy: an agent reads them
  on demand when a core's link names them; the link target is the heading's
  MyST anchor, and a broken one fails the docs build.
- The lessons index (`.claude/lessons.md`) is rendered like a rule but is not
  auto-loaded; the session-start batch reads it.
- Harness-specific text is added by the generator, never written in the docs:
  the `!cat` line that injects the error-catalogue index into `vv-principles`,
  the `GENERATED` stamp, the AGENT.md markers.
- An AGENT.md's header (tools, model, memory, `omitClaudeMd`) is
  harness-specific and hand-maintained; everything after it is the agent's
  definition (role, method, return contract), generated whole between the
  markers, so every check below reads every line of it (since 2026-09-22;
  until then only a ten-line role block was generated, and the hand-written
  bodies, 125 to 654 lines each, were read by no check,
  `.claude/plans/archive/agent_definitions.md` D1).
- A hand edit of generated text is refused as it is made:
  `.claude/hooks/generated-guard.py`, a PreToolUse hook on Edit, Write and
  MultiEdit registered in `.claude/settings.json`, fires in the main agent and
  in every sub-agent, and refuses an edit of a stamped file, an edit that
  changes the text between `GENERATED` markers, and an edit of a file an
  installer recorded in its manifest; the refusal names the source to edit and
  the command that regenerates. A write through Bash is outside its reach, and
  the drift check of `--check` in CI is the backstop
  (`tests/gates/tools/test_write_guards.py`, one red per arm).
- Installed files: `nexus setup` writes the extension's own skills, hooks and
  its always-on routing rule `nexus-tools` into `.claude/` and records each in
  `.claude/nexus-install-manifest.json`. They are tracked, never generated and
  never edited here: an edit is a fork upstream cannot update, and
  `nexus setup --check` reports it as locally modified. The project's own
  search and graph clauses are the `code-search` rule, generated like any other.
- Two sources name `.claude/rules/` as where the rules live (the onboarding
  page's "Where the rules are", the cardinal page's "the floor"). With one
  harness the path is the reader's fact; when a second harness exists the
  sources stop naming it and the generator writes it, as it writes the stamp
  and the AGENT.md markers today (`[R]` ruled 2026-09-21; `[M]` the same day,
  over the onboarding page and the rule and skill sources, `lessons.md`
  excluded: 8 mentions of `.claude/` in 6 sources, the other 6 naming
  directories no harness generates).

## Adding or changing

- **A rule**: write `docs/development/rules/<name>.md` (an imperative with its
  `check:` and `tell:` wherever it names a failure; a link to its evidence once a
  founding case exists, since a rule distilled from analysis rather than from a
  surprise starts without one), declare `harness: {kind: rule, budget_tokens: N}`
  in its front matter (a round figure above the measured size, never more than
  400 above it: `tools/harness/budget.py`), add the page to `index.rst`, run the
  generator. A rule with no `paths:` is always-on for every Key agent and the
  main agent: it earns that only if it applies to every artefact an agent
  writes. A rule whose every clause bites only under some paths declares them
  (`paths:`), and loads when a matching file is READ, in the main agent and
  in a sub-agent alike (`[M]` 2026-09-21: a haiku sub-agent that read
  `tests/conftest.py` received `vv-testing` as a system reminder after the
  read, none before); a WRITE of a new file under the scope does not load it
  (`[M]` 2026-09-22: a haiku sub-agent that wrote a new file under
  `.claude/plans/` received no `plan-authoring`), so an edit, which the Edit
  tool makes only after a read, gets the rule, and a page created from nothing
  does not unless its author reads the rule first. When a
  Support brief must carry one sentence of the rule (the three Support agents
  load no rule), that sentence is the block's `brief:`; `tools/harness/brief.py`
  assembles every rule's brief into the generated list on
  [the workflows page](workflows.md#the-brief), so the copy the brief template
  needs is a build product and a hand edit to it is drift.
- **A citation by ID is checked.** `--check` resolves every plain-text ID a
  core cites (X1–X4, Cardinal Rule N, Pattern N, `vv-principles` or
  `coding-elegance` #N, mode N, ERR-NNN, Lnn, a coding-standards item such as
  B.4 of the retirement-audit skill, a plan-authoring tag with two or more hyphens) against the page that
  defines it (`tools/harness/ids.py`), the way the build checks a link's anchor;
  a renamed or retired definition reddens at every site that still cites it.
  The evidence pages are under check too (since 2026-09-21: a citation to a
  retired numbering is dead for today's reader, whatever the age of the
  sentence around it; the first run over the eight pages found two, both
  re-pointed). Not checked, so a clean run is read for what it is: a bare
  `#N` with no page named before it in the paragraph (it reads as an issue
  number and is a problem in its own right, except on a page whose subject is
  the numbered list, where it resolves there, and after `task`, `issue`,
  `nexus` or `PR`, where it is one), `L1`–`L4` (also the V&V levels), a tag
  with fewer than two hyphens, and `mode-0` (ERR-030's normalisation sense;
  modes are numbered from 1).
- **A file named in a code span is checked.** `--check` resolves every
  repository path, dotted module and bare file name a generated page names
  in a code span against the tracked tree (`tools/harness/paths.py`); a
  renamed or deleted file reddens at every page that still names it, which
  no build warning does, since a code span is not a link and the ID check
  skips code. A file the reader is to create, and an example, carry a
  placeholder (`diag_<topic>_01_characterize.py`); a span with a placeholder,
  a glob or whitespace is not read, nor is a fence, a bare identifier or a
  path git ignores (local by design). `[M]` 2026-09-22, its first run: 1
  dead file in the skills and rules (a skill citing an untracked scratch
  report; re-pointed to its evidence case), then 17 in the agent bodies once
  they moved under generation (`error_catalog.md` for the `.rst`, a moved
  module, a test file that never existed, output paths written as if they
  existed); 126 resolved after the repair.
- **A skill**: `docs/development/skills/<name>.md` with YAML front matter
  (`name`, `description` — the harness reads them) and the `harness:` block
  beside them, the page in `index.rst`; preload it from an agent's `skills:` list when that
  agent needs it at every dispatch.
- **An agent**: `docs/development/agents/<name>.md` (role, phases, supports,
  method, return contract: the whole definition, budgeted at its measured
  size like any page) with its `harness:` block, the page in `agents/index.rst`. An
  agent whose role writes nothing tracked (a reviewer, a Support agent)
  carries a `write-scope` hook in its front matter
  (`python3 "$CLAUDE_PROJECT_DIR/.claude/hooks/write-scope.py" <name>` on Edit, Write and
  MultiEdit), which refuses a write outside `scratch/`, the temporary
  directory and its own memory, so "read-only" is enforced rather than asked
  for (qa, elegance-enforcer, explorer, cross-domain-attacker,
  literature-researcher). The
  hand-maintained AGENT.md header decides the model, the tools (a Support agent
  omits `Agent`; no agent lists `Skill`), the memory scope and
  `omitClaudeMd`; an edit to `tools:` is live from the next harness start, an
  edit to `skills:` or the definition from the next dispatch.
- **Distilling a page into a core**: every clause keeps its imperative, its
  `check:` and its `tell:`; a clause's check must reach every mechanism its text
  names, so two mechanisms one check cannot reach are two clauses; a count in an
  appendix states its convention and is re-run after the edit. A lesson retires
  only when a clause carries each of its body's own numbered rules, matched by
  its check and never by its title.
- **A surprise, a lesson, a founding case**: append to the evidence page as a
  `###` heading (date + words) and link it from the clause it instances; a
  mechanism that recurs is a signal for a tool, not a paragraph.
- **A pointer to a definition does not restate any part of it.** A gloss beside
  a link ("whose census clause is X and Y") asserts the target's scope, reads
  as navigation so nobody audits it, and goes stale when the target grows
  (`[M]` 2026-09-20: the explorer role block named two of the census clause's
  four requirements, staled by the commit that added the fourth). Point, and
  let the definition speak.
- **A plan** stays in `.claude/plans/` (transient: executed, triaged, archived);
  its close-out record moves here.

- **A claim about what a dispatch receives is measured by dispatching**, never
  read off the frontmatter or a plan. `[M]` 2026-09-20: a zero-tool explorer
  probe, asked what it held, showed that `omitClaudeMd` drops CLAUDE.md, the
  rules and the PROJECT memory index but not the agent's own memory index or
  its preloaded skills; six role blocks had said "no memory index" from
  reading the F10 measurement, which had looked only at the project index.
  The same probe fixture (haiku, one-line system prompt, no tools) carries no
  role block, so a role-block change cannot move the T4 dispatch floor: a Key
  dispatch costs the fixed harness block, the instruction files, its AGENT.md,
  its preloaded skills and its own memory index, and each part is measured
  by a probe that carries it. `[M]` 2026-09-21: an agent's `tools:` allowlist is
  read when the harness starts, and its `skills:` list and definition at each
  dispatch: qa and method-implementer, with `Skill` added to their lists, held
  no `Skill` tool until Claude Code was restarted, then held it and loaded
  `nexus-guide` on demand, while the `retirement-audit` preload added to qa's
  `skills:` reached the next dispatch of the same session. No Key agent lists
  `Skill` (ruling 2026-09-21: a preload is specific and reproducible per
  dispatch, and the tool carries its own context cost), so a sub-agent's
  skills are its `skills:` preload and a rule that says "load skill X"
  reaches a sub-agent only through that list. The cost, `[M]` 2026-09-21 by
  the two-arm probe (haiku, `omitClaudeMd`, `tools: [Read]` against
  `tools: [Read, Skill]`, a fresh headless parent): listing `Skill` injects
  the roster of every skill the environment holds with its description, 77
  entries including every installed plugin's, for **2 613** first-turn tokens
  (5 822 against 3 209); the roster is the environment's, not the project's,
  so the cost grows with plugins the project never uses.

## Session start

`.claude/hooks/session-start.txt` (printed by a SessionStart hook after the
environment-health gates) names the batch the main agent loads: the
`vv-principles`, `coding-elegance` and `instrument-doctrine` cores, the lessons
index, and the Nexus session briefing. Rules and CLAUDE.md are a session-start
snapshot: an edit to a rule core is live from the NEXT session; an edit to a
skill is live at its next invocation; an edit to an AGENT.md is live at its
next dispatch.

## Acceptance instruments

`/context` lists the loaded memory files and their size; the `InstructionsLoaded`
hook logs which instruction files loaded and why; a zero-tool haiku probe
dispatched with a fixed question set (`scratch/_harness_eval/probe_questions.md`)
measures a dispatch's floor from its transcript's first-turn usage. The `gates`
workflow (`.github/workflows/gates.yml`) runs `--check`, the harness tests and
the strict build on every push, so a hand edit inside a generated file reds in
CI as well as locally ([the git workflow page](git_workflow.rst), "Continuous integration").
