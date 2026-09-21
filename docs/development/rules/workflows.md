---
harness:
  kind: rule
  budget_tokens: 1900
---

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
| **Support** | answers one question from its brief and returns; never spawns; its brief carries its rules (invariant 4) | no | explorer (available to any agent), literature-researcher, cross-domain-attacker — all `omitClaudeMd`; haiku categorisers (`general-purpose` with a fixed output schema; they inherit the rules) |

## Posture

Dispatch freely; never ask permission. Sub-agents (explorer,
numerics-investigator, literature-researcher, qa, archivist, test-architect,
…) are dispatched during investigation and implementation without pausing to
ask "shall I dispatch X?", and independent investigations run in parallel: the
user values momentum and parallelism, and the fleet's existence is already
approved. Every agent's output is still reviewed with full session context
before committing; the agent had none. The built-in Explore agent is denied;
`explorer` is the project's.

**Exception — surgical, high-correctness carves: the main agent writes
directly.** For operator-algebra carves, convention changes crossing three or
more subsystems, and anything in the `refactor/sn-operator-algebra` family,
`method-implementer` is NOT dispatched: the main agent writes the code with the
user steering step by step and `AskUserQuestion` checkpoints. The reason (`[R]`, the user's
ruling of 2026-05-20): the user's since-inception knowledge of the codebase
corrects the implementation in real time, while `method-implementer` runs a
brief to completion, so for surgical work the loss of turn-by-turn correction
costs more than the parallelism gains. The constraint is on
`method-implementer` alone (batch code production); test-architect, explorer,
qa and archivist remain available and encouraged. Any "surgical" framing by
the user selects this mode, and `method-implementer` re-enters scope at the
user's signal, when routine refactor cycles resume.

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
   that apply to it (the template's "Rules that apply to you" line):
   `omitClaudeMd` drops CLAUDE.md, every rule and the project memory index,
   not the agent's own memory or its preloaded skills (measured on the
   [harness page](../harness.md)).

## The workflows

One routing line each; the phases, the gates and the return contract are on
[the workflows page](../workflows.md).

- **W1 Build a capability** — a new solver, operator or method: context from
  explorer (and literature-researcher for a published formulation), the
  verification spec from **test-architect**, the build by
  **method-implementer** or the main agent for a surgical carve, **qa** and
  **elegance-enforcer** in parallel, **archivist**, close-out.
  [W1](../workflows.md#w1--build-a-capability)
- **W2 Wrong answer** — **numerics-investigator**'s probe cascade, the fix,
  **qa** with a mutation that re-introduces the defect and
  **elegance-enforcer** in parallel, **archivist** for the ERR entry.
  [W2](../workflows.md#w2--wrong-answer)
- **W3 Surgical carve** — the main agent writes with the user steering;
  explorer for the blast set, **test-architect** for gates and re-baselines,
  review as W1. [W3](../workflows.md#w3--surgical-carve)
- **W4 Documentation campaign** — **archivist**, qa for claim verification;
  gates `-W`, `dead_references`, `staleness`.
  [W4](../workflows.md#w4--documentation-campaign)
- **W5 Design review** — **cross-domain-attacker** and **elegance-enforcer**
  on a first-pass design; the output feeds W1's build.
  [W5](../workflows.md#w5--design-review)
- **W6 Tree-wide census** — haiku categorisers with a fixed output schema; the
  orchestrator aggregates; no nesting. [W6](../workflows.md#w6--tree-wide-census)
- **W7 Literature acquisition** — **literature-researcher**, the local folder
  first; "not in the local folder" is a question to the user, never a pivot.
  [W7](../workflows.md#w7--literature-acquisition)

## Dispatch facts

Nested dispatch works to three layers and a sub-agent's `tools:` allowlist
decides whether it holds `Agent`; sub-agents address named siblings with
`SendMessage`, and a finished agent resumes on message with its full history
(measured: [the page's history](../workflows.md#history)). What a dispatch
inherits and what it costs, measured by probe, is
[the harness page](../harness.md#what-loads-and-what-it-costs).
