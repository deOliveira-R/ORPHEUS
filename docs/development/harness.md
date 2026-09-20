# The agent harness — what loads when, and how to change it

ORPHEUS is developed by a user steering Claude Code agents. The knowledge those
agents work from lives in this section of the docs; the `.claude/` directory is
one harness's **generated** view of it. This page records what loads into an
agent's context, what it costs, and the procedure for adding a rule, a skill or
an agent. Numbers are `[M]` from the 2026-09 evaluation
(`.claude/plans/harness_context_budget.md`, Parts I and VII).

## What loads, and what it costs

| block | who pays | before (2026-09-04) | after this restructure |
|---|---|---|---|
| CLAUDE.md + `.claude/rules/*.md` + the auto-memory index — inherited by every dispatch that does not set `omitClaudeMd` | main agent and every Key agent | ≈64.6K tokens | budgeted in `tools/docs/harness_manifest.toml` (rule cores ≈17K + the untouched rules ≈7K + memory ≈2.3K) |
| the session-start batch (`.claude/hooks/session-start.txt`): skill cores, the lessons index, the Nexus briefing | main agent only | ≈125K | ≈30K |
| a Support agent with `omitClaudeMd: true` | explorer, literature-researcher, cross-domain-attacker, haiku categorisers | ≈107K per dispatch (haiku, zero tools) | ≈36.5K — the harness-fixed block only (tool schemas, skill index, roster) |

Measured 2026-09-20 on Claude Code 2.1.278: `omitClaudeMd: true` drops CLAUDE.md,
every rule file and the memory index together; the Nexus tool schemas (≈13K)
remain, since the harness loads them eagerly.

## The one-way flow

```
docs/development/{rules,skills,agents}/*.md  +  lessons.md        (SOURCE, MyST)
        │  tools/docs/generate_harness.py  (a `_GENERATORS` row in docs/conf.py;
        │  `--check` in CI: drift, budget, dead heading links)
        ▼
.claude/rules/*.md   .claude/skills/*/SKILL.md   .claude/lessons.md   AGENT.md role blocks   (GENERATED, committed)
```

- Evidence pages (`evidence/*.md`) have no `.claude/` copy: an agent reads them
  on demand when a core's link names them; the link target is the heading's
  MyST anchor, and a broken one fails the docs build.
- Harness-specific text is added by the generator, never written in the docs:
  the `!cat` line that injects the error-catalogue index into `vv-principles`,
  the `GENERATED` stamp, the AGENT.md markers.
- An AGENT.md's header (tools, model, memory, `omitClaudeMd`) is
  harness-specific and hand-maintained; only the role block between the
  markers is generated.

## Adding or changing

- **A rule**: write `docs/development/rules/<name>.md` (imperative + `check:` +
  `tell:` + a link to its evidence), add a `[[rule]]` entry with a
  `budget_tokens` to the manifest, add the page to `index.rst`, run the
  generator. A rule is always-on for every Key agent and the main agent: it
  earns that only if it applies to every artefact an agent writes.
- **A skill**: `docs/development/skills/<name>.md` with YAML front matter
  (`name`, `description` — the harness reads them), a `[[skill]]` entry, the
  page in `index.rst`; preload it from an agent's `skills:` list when that
  agent needs it at every dispatch.
- **An agent**: `docs/development/agents/<name>.md` (role, phases, supports,
  return contract), an `[[agent]]` entry, the page in `agents/index.rst`; the
  hand-maintained AGENT.md header decides the model, the tools (a Support agent
  omits `Agent`), the memory scope and `omitClaudeMd`.
- **A surprise, a lesson, a founding case**: append to the evidence page as a
  `###` heading (date + words) and link it from the clause it instances; a
  mechanism that recurs is a signal for a tool, not a paragraph.
- **A plan** stays in `.claude/plans/` (transient: executed, triaged, archived);
  its close-out record moves here.

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
