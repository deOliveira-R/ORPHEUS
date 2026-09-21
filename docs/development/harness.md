# The agent harness — what loads when, and how to change it

ORPHEUS is developed by a user steering Claude Code agents. The knowledge those
agents work from lives in this section of the docs; the `.claude/` directory is
one harness's **generated** view of it. This page records what loads into an
agent's context, what it costs, and the procedure for adding a rule, a skill or
an agent. Measured numbers are marked `[M]` (the 2026-09 evaluation,
`.claude/plans/harness_context_budget.md`, Parts I, VII and VIII); a budget is a
budget, not a measurement.

## What loads, and what it costs

| block | who pays | before (2026-09-04) | after this restructure |
|---|---|---|---|
| CLAUDE.md + `.claude/rules/*.md` + the auto-memory index — inherited by every dispatch that does not set `omitClaudeMd` | main agent and every Key agent | ≈64.6K tokens | budgeted in each page's `harness:` front matter; every always-on file is generated (`tests/test_harness_generated.py` fails on a hand-maintained rule), so the sum `python -m tools.harness --check` prints is the whole block but the memory index: eight rule cores + CLAUDE.md ≈21.4K (`[M]` 2026-09-21 chars/3.6, after the re-evaluation; `vv-testing` ≈1.0K is path-scoped under `tests/**` and loads when a matching file is touched, in the main agent and in a sub-agent alike), with the memory index ≈2.5K on top. The keep − omit probe read **13 850** tokens with `coding-standards` and `plan-authoring` path-scoped as well (`[M]` 2026-09-21, haiku first turn on `f5f2cf2d`, 49 288 chars over 8 files at 3.56 chars/token); the user then ruled both always-on for now, which adds their ≈9.8K back (a re-run is owed); 28 705 after K3b, 30 655 after P1, 70 966 before the restructure |
| the session-start batch (`.claude/hooks/session-start.txt`): skill cores, the lessons index, the Nexus briefing | main agent only | ≈125K | ≈30K |
| a Support agent with `omitClaudeMd: true` | explorer, literature-researcher, cross-domain-attacker, haiku categorisers | ≈107K per dispatch (haiku, zero tools) | ≈36.5K — the harness-fixed block only (tool schemas, skill index, roster); `[M]` 2026-09-21 the omit arm reads 36 536 |

Measured 2026-09-20 on Claude Code 2.1.278: `omitClaudeMd: true` drops CLAUDE.md,
every rule file and the memory index together; the Nexus tool schemas (≈13K)
remain, since the harness loads them eagerly.

## The one-way flow

```
docs/development/{rules,skills,agents}/*.md + lessons.md + onboarding.md     (SOURCE, MyST; each
        │                                        declares `harness: {kind, budget_tokens}` in its front matter)
        │  tools/harness/  (`python -m tools.harness`, a `_GENERATORS` row in docs/conf.py; `--check`, read by
        │  `tests/test_harness_generated.py`: drift, budget, dead links, orphans. `targets/claude_code.py` is the
        │  one harness; the other modules name none, and adding a harness is one module under `targets/`)
        ▼
.claude/rules/*.md   .claude/skills/*/SKILL.md   .claude/lessons.md   AGENT.md role blocks   CLAUDE.md   (GENERATED, committed)
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
  harness-specific and hand-maintained; only the role block between the
  markers is generated.

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
  (`paths:`), and loads when a matching file is touched, in the main agent and
  in a sub-agent alike (`[M]` 2026-09-21: a haiku sub-agent that read
  `tests/conftest.py` received `vv-testing` as a system reminder after the
  read, none before). When a
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
  Not checked, so a clean run is read for what it is: a bare `#N` with no page
  named before it in the paragraph (it reads as an issue number and is a
  problem in its own right), `L1`–`L4` (also the V&V levels), a tag with fewer
  than two hyphens, and the evidence pages.
- **A skill**: `docs/development/skills/<name>.md` with YAML front matter
  (`name`, `description` — the harness reads them) and the `harness:` block
  beside them, the page in `index.rst`; preload it from an agent's `skills:` list when that
  agent needs it at every dispatch.
- **An agent**: `docs/development/agents/<name>.md` (role, phases, supports,
  return contract) with its `harness:` block, the page in `agents/index.rst`; the
  hand-maintained AGENT.md header decides the model, the tools (a Support agent
  omits `Agent`; no agent lists `Skill`), the memory scope and
  `omitClaudeMd`; an edit to `tools:` is live from the next harness start, an
  edit to `skills:` or the role block from the next dispatch.
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
  read when the harness starts, and its `skills:` list and role block at each
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
measures a dispatch's floor from its transcript's first-turn usage.
