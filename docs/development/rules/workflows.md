# Workflows — roles, invariants, and the routes work takes

Routing guidance, not a restriction: the seven workflows below are the routes work
usually takes; an agent dispatches whenever the reason is worthwhile. The full
phase descriptions, the brief template and the return contract are in
[workflows](../workflows.md).

## Roles

| role | mandate | may spawn | agents |
|---|---|---|---|
| **Orchestrator** | owns phase transitions, the review stage, user rulings, issues, commits | yes | the main agent |
| **Key** | owns one phase of a workflow; calls Support agents freely, and any other Key agent when the reason is worthwhile (an implementer spawning a numerics-investigator keeps its own context clean) | yes | test-architect, method-implementer, numerics-investigator, archivist, qa, elegance-enforcer |
| **Support** | answers one question from its brief and returns; never spawns; the three `omitClaudeMd` agents are launched without the project rules, so the brief is the only place a project rule reaches them | no | explorer (available to any agent), literature-researcher, cross-domain-attacker — all `omitClaudeMd`; haiku categorisers (`general-purpose` with a fixed output schema; they inherit the rules) |

## Invariants

1. **The parent's review is independent of the child's.** A child may spawn
   any agent it needs, including qa; the parent still dispatches its own qa and
   elegance-enforcer on the result. Independence is guaranteed at the parent,
   where it is cheap.
2. **Depth 3 is the ceiling, 2 the norm**: the orchestrator dispatches a Key
   agent, which dispatches a Support agent. The harness default
   (`CLAUDE_CODE_MAX_SUBAGENT_SPAWN_DEPTH`) is 3; Support allowlists omit
   `Agent`.
3. **Continuity by name.** The test-architect that wrote a verification spec is
   resumed (`SendMessage`) at review time to confirm its gates landed; the
   implementer is resumed with the review findings. No re-briefing.
4. **Every brief carries** the workflow ID, the phase, the previous phase's
   artefact paths, and the return contract: a word cap (Opus runs long), files
   carry the detail, and a `NEEDS:` block for anything the agent could not
   obtain. A brief to one of the three `omitClaudeMd` agents also carries the rules
   that apply to it (the template's "Rules that apply to you" line): it sees
   no project rule.

## The workflows

- **W1 Build a capability** — P0 context: explorer, plus literature-researcher
  for a published formulation. P1 verification design: **test-architect** (a
  spec whose every gate names its first red). P2 build: **method-implementer**,
  or the main agent for a surgical carve; supports explorer,
  literature-researcher, cross-domain-attacker after the first pass. P3 review,
  in parallel, dispatched by the parent: **qa** and **elegance-enforcer**; any
  red resumes the implementer by name. P4 documentation: **archivist** (theory
  page, changelog; `sphinx -W` and `dead_references` clean). P5 close-out:
  issues, retirement audit, commit.
- **W2 Wrong answer** — P1 **numerics-investigator** (the probe cascade; may
  call test-architect for the permanent test). P2 the fix. P3 **qa** and
  **elegance-enforcer** in parallel, dispatched by the parent (qa with a mutation
  that re-introduces the defect). P4 **archivist** (the ERR
  entry). P5 close-out.
- **W3 Surgical carve** — the main agent writes with the user steering;
  explorer for the blast set; **test-architect** for gates and re-baselines;
  review as W1; **archivist** for the changelog.
- **W4 Documentation campaign** — **archivist**; qa for claim verification;
  gates `-W`, `dead_references`, `staleness`.
- **W5 Design review** — **cross-domain-attacker** + **elegance-enforcer** on a
  first-pass design; the output feeds W1's build phase.
- **W6 Tree-wide census** — haiku categorisers with a fixed output schema; the
  orchestrator aggregates; no nesting.
- **W7 Literature acquisition** — **literature-researcher**: the local folder
  first, then OCR sidecars; "not in the local folder" is a question to the
  user, never a pivot to a secondary source.

## Dispatch facts (measured 2026-09-20, Claude Code 2.1.278)

Nested dispatch works to three layers; a sub-agent's `tools:` allowlist decides
whether it holds `Agent`. Sub-agents carry `SendMessage` and can address named
siblings; a finished agent resumes on message with its full history. A dispatch
inherits CLAUDE.md, every rule file and the memory index (`[M]` 2026-09-20, the
keep − omit haiku probe at the first turn: 30 655 tokens after this restructure,
70 966 before it) unless the agent sets `omitClaudeMd: true`.
