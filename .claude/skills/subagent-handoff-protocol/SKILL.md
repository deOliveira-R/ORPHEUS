---
name: subagent-handoff-protocol
description: 'PROACTIVELY load whenever a sub-agent expects to coordinate with another sub-agent during its task. Documents the orchestrated-dispatch pattern that bridges around the Anthropic constraint that sub-agents cannot call the Agent tool themselves (no nested dispatch). Defines the bidirectional contract: sub-agents return DISPATCH_REQUEST blocks; the main agent dispatches/resumes on their behalf and returns DISPATCH_RESULT blocks with the dispatched-agent ID so the requester can address that specific instance for follow-up. Read by both sub-agents (request side) and the main agent (orchestration side). Examples: "I need a literature pull but I''m a sub-agent", "How do I follow up with the same numerics-investigator that ran earlier", "How does the main agent route between sub-agents", "Why can''t my agent dispatch another agent directly".'
---

# Subagent Handoff Protocol — orchestrated dispatch via the main agent

This skill institutionalizes the bridge that lets sub-agents coordinate
with other sub-agents through the main agent. It is read by **both
audiences**:

- **Sub-agents** read it to learn how to request a dispatch and how to
  follow up with a specific previously-dispatched instance.
- **The main agent** reads it (transitively, via CLAUDE.md) to learn
  the reciprocal orchestration responsibilities.

The protocol exists because of an Anthropic platform constraint, not
a project preference. Quoted directly from the Claude Code subagents
docs (`https://code.claude.com/docs/en/sub-agents`):

> This prevents infinite nesting (subagents cannot spawn other
> subagents) while still gathering necessary context.

The Agent tool is unavailable inside a sub-agent's context. Any
attempt to dispatch directly from a sub-agent will fail. All
cross-sub-agent coordination MUST flow through the main agent.

---

## CRITICAL: One-line summary

**Sub-agent**: When you need help, return a structured
`DISPATCH_REQUEST` block. The main agent will dispatch on your
behalf and resume you (via SendMessage) with a `DISPATCH_RESULT`
block that includes the dispatched agent's `agent_id` — capture
that ID if you want to address that specific instance for follow-up.

**Main agent**: When a sub-agent's return contains one or more
`DISPATCH_REQUEST` blocks, parse them, dispatch (or resume by ID),
and return each result wrapped in a `DISPATCH_RESULT` block carrying
the dispatched agent's ID. You may APPEND a `MAIN_AGENT_AUGMENTATION`
sub-block to a brief before dispatch to leverage your wider context;
NEVER modify the requesting sub-agent's original brief content. You
are a transparent bridge with optional augmentation — depth from
sub-agents, breadth from you.

---

## CRITICAL: The depth-vs-breadth principle

**Sub-agents have depth, the main agent has breadth.**

A sub-agent operates with deep focus on a narrow task — it knows the
specifics of its current investigation but lacks visibility into:

- What other sub-agents are working on simultaneously.
- The full conversation history with the user.
- Recent commits, plans, or decisions outside its briefing.
- Side discussions about strategy or scope that didn't reach its brief.

The main agent has the opposite profile: shallow per-thread but with
full lateral visibility — every sub-agent's prior result, the user's
goals, the broader project state.

The augmentation mechanism (defined later in this skill) lets the main
agent piggyback on a dispatched sub-agent's investigation to ask
ADDITIONAL questions that benefit from that sub-agent's depth, while
keeping the requesting sub-agent's brief and answer intact and
uncontaminated. This is how the system gets both depth (per
investigation) AND breadth (across investigations) without the
requesting sub-agent having to anticipate every question the main
agent might want answered.

---

## Agent table — which sub-agent to request

The project defines specialized sub-agents in `.claude/agents/`. Use
this table to choose the right `agent` field value when emitting a
`DISPATCH_REQUEST`. The main agent reads the same table when
fielding ad-hoc dispatches.

Each agent has a focus, a "proactively invoke when" trigger, and
preloaded skills that shape its behavior. Sub-agents whose preload
list includes `subagent-handoff-protocol` (THIS skill) can themselves
emit `DISPATCH_REQUEST` blocks; sub-agents without it must complete
in one shot or return to the user.

| Agent                     | Proactively invoke when                          | Preloaded skills                                                                                                            |
| ------------------------- | ------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------- |
| **explorer**              | Understanding code, tracing dependencies         | `nexus-exploring`, `nexus-guide`                                                                                            |
| **archivist**             | Writing/reviewing Sphinx docs                    | `nexus-verification`, `nexus-exploring`, `vv-principles`, `subagent-handoff-protocol`                                       |
| **qa**                    | Reviewing code, validating claims                | `nexus-verification`, `nexus-impact`, `nexus-debugging`, `vv-principles`, `numerical-bug-signatures`, `subagent-handoff-protocol` |
| **numerics-investigator** | Solver gives wrong answers                       | `nexus-debugging`, `nexus-impact`, `probe-cascade`, `vv-principles`, `numerical-bug-signatures`, `subagent-handoff-protocol` |
| **literature-researcher** | Need equations from papers / external references | `research`, `subagent-handoff-protocol`                                                                                     |
| **test-architect**        | Planning verification BEFORE implementation      | `nexus-verification`, `nexus-impact`, `vv-principles`, `numerical-bug-signatures`, `subagent-handoff-protocol`              |
| **cross-domain-attacker** | Detecting structural patterns                    | `cross-domain-frames`, `subagent-handoff-protocol`                                                                          |

(The `subagent-handoff-protocol` entry on each row is the target
state after rollout — existing AGENT.md frontmatter may not yet list
it. Adopting this skill across the agent fleet is the rollout step
that makes nested-coordination requests work for every agent.)

### Choosing the right agent

- **Need code understanding / tracing / file location** → `explorer`
- **Need a Sphinx page written or rewritten with full math context**
  → `archivist`
- **Reviewing correctness / V&V evidence claims** → `qa`
- **A solver gives a wrong-but-plausible answer** →
  `numerics-investigator`
- **Need an equation, citation, or reference value from a paper /
  literature corpus** → `literature-researcher`
- **Designing the verification plan BEFORE writing the
  system-under-test** → `test-architect`
- **Sense that the formulation might not match the problem's native
  mathematical structure** → `cross-domain-attacker`
- **Anything else / general** → `general-purpose`

If you can't find a clean match, prefer `general-purpose` with a very
explicit brief over forcing a poor match.

---

## The DISPATCH_REQUEST block (sub-agent → main)

When a sub-agent needs another sub-agent's help, it appends one or
more blocks of this shape to its final return:

```
--- DISPATCH_REQUEST ---
agent: <agent-type>
agent_id: <previously-returned-ID — OPTIONAL, see below>
brief: |
  <self-contained prompt for the dispatched agent.
   Assume zero context from your work. Include all
   relevant file paths, prior results, exact ask,
   expected output format.>
wait_for: <one-line description of what you need back>
followup: <true | false — whether you want to be resumed after the dispatch>
--- END DISPATCH_REQUEST ---
```

### Field semantics

| Field      | Required | Meaning                                                                                                                                                                                                                                         |
| ---------- | :------: | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `agent`    |   Yes    | The agent-type to invoke. Examples: `literature-researcher`, `numerics-investigator`, `archivist`, `cross-domain-attacker`, `test-architect`, `general-purpose`.                                                                                |
| `agent_id` |    No    | A specific previously-dispatched sub-agent's ID. **If provided**, the main agent resumes that exact instance via SendMessage (preserving its full prior context). **If absent**, the main agent dispatches a FRESH instance via the Agent tool. |
| `brief`    |   Yes    | Self-contained prompt. Include file paths, prior results, exact ask, expected output format. The dispatched agent has no prior context from your work.                                                                                          |
| `wait_for` |   Yes    | One-line summary of what result you need. Used by the main agent for traceability and to filter or summarise long results.                                                                                                                      |
| `followup` |   Yes    | `true` = main agent must resume you (via SendMessage) with the result. `false` = fire-and-forget; you've completed your task by delegating.                                                                                                     |

### Multiple blocks per return

A sub-agent MAY emit multiple `DISPATCH_REQUEST` blocks in one return.
The main agent will dispatch them in parallel (one Agent tool call per
block, sent in a single message for efficiency). All `DISPATCH_RESULT`
blocks come back together in the resume message.

If you intend a SEQUENCE of dispatches (B's brief depends on A's
result), do them one at a time across multiple resume cycles, not as
parallel requests.

---

## Main-agent augmentation (depth-vs-breadth in action)

When the main agent receives a `DISPATCH_REQUEST` from sub-agent A
asking for help from agent type B, the main agent MAY append an
additional sub-block to the brief that asks B for ADDITIONAL information
that benefits from B's depth — questions A didn't ask but the main
agent (with breadth across the whole session) knows are valuable.

The dispatched B receives a brief with this two-part structure:

```
--- ORIGINAL_BRIEF ---
<verbatim brief from A — UNMODIFIED>
--- END ORIGINAL_BRIEF ---

--- MAIN_AGENT_AUGMENTATION ---
<additional questions from the main agent, leveraging session breadth>
--- END MAIN_AGENT_AUGMENTATION ---
```

B then returns answers to BOTH parts in matching blocks:

```
--- ORIGINAL_RESPONSE ---
<B's answer to A's brief>
--- END ORIGINAL_RESPONSE ---

--- AUGMENTATION_RESPONSE ---
<B's answer to the main agent's augmented questions>
--- END AUGMENTATION_RESPONSE ---
```

The main agent then resumes A with **both** responses inside the
`DISPATCH_RESULT` block. A reads the `ORIGINAL_RESPONSE` to continue
its task. A MAY also read the `AUGMENTATION_RESPONSE` if relevant —
it sees what the main agent asked and what B answered. This keeps
the system transparent: A knows what extra context the main agent
extracted via its dispatch.

### Rules for augmentation

- **The main agent MUST NOT modify the `ORIGINAL_BRIEF`.** A's brief
  goes through verbatim. Augmentation is APPEND-only.
- **The main agent SHOULD augment when it has lateral context A lacks**
  — e.g., "another sub-agent X is also looking at this; please flag
  any inconsistency", or "the user's stated goal is Y; please mention
  if your finding affects Y".
- **The main agent SHOULD NOT augment with questions that just
  duplicate A's brief** in different words. The augmentation is for
  questions A could not have asked because A doesn't have the breadth.
- **Augmentation is OPTIONAL.** If the main agent has nothing extra
  to ask, B receives only the ORIGINAL_BRIEF (without delimiters in
  this case — just the brief content alone, as in the simple
  protocol). B's response then has no AUGMENTATION_RESPONSE block.
- **B MUST keep responses separate** when given an augmented brief.
  Mixing the two contaminates the per-audience context. If B can't
  answer the augmentation but answered the original, B should still
  emit an empty `--- AUGMENTATION_RESPONSE ---` block with a single
  line explaining why.
- **A SHOULD treat the `AUGMENTATION_RESPONSE` as informational, not
  authoritative for its own task.** The augmentation is the main
  agent's investigation, not A's; A's deliverables are based on the
  `ORIGINAL_RESPONSE`. A may flag if the AUGMENTATION_RESPONSE reveals
  something that contradicts A's brief, but A is not obligated to
  re-scope based on it.

### Why this design

The pattern preserves the requesting sub-agent's epistemic boundary
(its task, its evidence, its conclusions) while giving the main
agent a token-efficient way to harvest depth from a sub-agent's
investigation that's already happening. Without augmentation, the
main agent would either:

1. **Skip the questions** — losing useful information.
2. **Dispatch a separate parallel agent** with a similar brief — paying
   double the token cost for overlapping work.
3. **Modify A's brief** — contaminating A's task with side concerns A
   didn't ask for and might not be equipped to handle.

Augmentation gives a fourth, cleaner path: piggyback transparently,
keep audiences separated, return both responses to A so A understands
what extra context was extracted on the side.

---

## The DISPATCH_RESULT block (main → sub-agent)

When the main agent resumes you (via SendMessage), the message
contains one block per dispatched agent:

```
--- DISPATCH_RESULT ---
agent: <agent-type that ran>
agent_id: <ID of the sub-agent instance that ran — STORE THIS>
status: <completed | failed | timed_out>
result: |
  <the dispatched agent's full return, verbatim — including
   ORIGINAL_RESPONSE and AUGMENTATION_RESPONSE sub-blocks
   if augmentation was used>
--- END DISPATCH_RESULT ---
```

If the main agent augmented the brief, the `result` field contains
both the `ORIGINAL_RESPONSE` and `AUGMENTATION_RESPONSE` sub-blocks
verbatim. You read the `ORIGINAL_RESPONSE` for your task; the
`AUGMENTATION_RESPONSE` is informational about what extra context
the main agent extracted.

### MUST capture `agent_id`

If you (the requester) want any follow-up — clarification, refinement,
re-run with adjusted parameters — you can resume that exact instance
by including the `agent_id` in your next `DISPATCH_REQUEST`. Resuming
preserves the dispatched agent's prior context (it remembers what it
did and what it found); fresh dispatch loses that context and pays
the token cost to re-establish it.

Capture every `agent_id` you receive. The cost of storing the ID is
trivial; the cost of re-dispatching from scratch can be substantial.

### Failure handling

If a dispatched agent fails or times out, `status` will reflect that
and `result` will contain the error context. The requester decides
whether to:

- Retry with a refined brief (fresh dispatch — omit `agent_id`).
- Resume the same instance with a clarifying message (provide
  `agent_id`).
- Pivot to a different agent type entirely.
- Abort and return to the user with an error report.

---

## When to resume vs. fresh dispatch

The `agent_id`-or-not decision is the highest-leverage choice in this
protocol. Use the table below.

| Situation                                                                                                              | Use                                                |
| ---------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------- |
| Follow-up question on the SAME topic the agent already explored                                                        | **Resume** (provide `agent_id`)                    |
| Need elaboration on a specific point from their prior result                                                           | **Resume**                                         |
| Re-run the SAME check with adjusted parameters (the agent already loaded the relevant code/refs)                       | **Resume**                                         |
| Their prior context is load-bearing and re-establishing it costs many tokens                                           | **Resume**                                         |
| Different topic — even if same agent type                                                                              | **Fresh dispatch** (omit `agent_id`)               |
| Want a fresh perspective without anchoring on the prior conclusion                                                     | **Fresh dispatch**                                 |
| Prior dispatch finished cleanly and the next task is genuinely independent                                             | **Fresh dispatch**                                 |
| You're orchestrating a parallel comparison (two instances of the same agent type running independent investigations)   | **Fresh dispatch** (twice, possibly in one return) |
| Prior dispatch timed out or its context was cleaned up (`status: timed_out` or main agent reports SendMessage failure) | **Fresh dispatch** with synthetic recap            |

---

## Main-agent contract (the orchestrator's responsibilities)

This section is **for the main agent** and is referenced from
CLAUDE.md (the project-wide orchestration glue).

When the main agent receive a sub-agent's return:

1. **Scan the return for `--- DISPATCH_REQUEST ---` blocks** between
   their delimiters. Parse each block (YAML-style fields).

2. **For each block**, decide dispatch path:
   - If `agent_id` is provided → use **SendMessage** with `to:
<agent_id>` and the `brief` as the message body.
   - If `agent_id` is absent → use the **Agent tool** with
     `subagent_type: <agent>` and `prompt: <brief>`.

3. **Capture the dispatched agent's ID** from the dispatch result.
   For Agent tool calls this is the `agentId` returned in the result
   metadata. For SendMessage resumes the agent_id is the same ID
   that was used to address the message.

4. **If multiple blocks** are present and they're independent,
   dispatch them in parallel (one tool call per block, all sent in
   one message). If they're sequential or there's any reason to
   serialize, dispatch one at a time.

5. **Resume the original sub-agent** (via SendMessage to the
   requester) with one `DISPATCH_RESULT` block per dispatched agent,
   each carrying:
   - `agent` — the type that ran
   - `agent_id` — the ID of the instance that ran (this is the
     critical bit; the requester needs it for follow-up)
   - `status` — completed / failed / timed_out
   - `result` — the full return, verbatim

6. **If `followup: false`** on a DISPATCH_REQUEST, do NOT resume the
   original requester for that block. The dispatched agent's result
   is fire-and-forget; route it to the user if appropriate.

7. **Be a transparent bridge — with optional augmentation.** Do NOT:
   - Edit, summarize, or filter the requesting sub-agent's brief.
   - Edit, summarize, or filter results before resume.
   - Skip the `agent_id` field on resume.
   - Combine multiple results into a single block (one DISPATCH_RESULT
     per dispatched agent).

   You MAY:
   - **Append a `MAIN_AGENT_AUGMENTATION` sub-block** to the brief
     before dispatching, asking the dispatched agent for ADDITIONAL
     information that benefits from your wider context (other
     sub-agents' results, user intent, project state).
   - The dispatched brief then has structure: `--- ORIGINAL_BRIEF ---`
     (verbatim from requester) + `--- MAIN_AGENT_AUGMENTATION ---`
     (your additions). The dispatched agent will return matched
     `--- ORIGINAL_RESPONSE ---` and `--- AUGMENTATION_RESPONSE ---`
     sub-blocks. Pass BOTH back to the requester inside the
     `result` field of the `DISPATCH_RESULT`.
   - When augmenting, ensure the augmentation questions are
     genuinely additive (use your breadth) and not duplicative of
     the requester's brief (which would just confuse the dispatched
     agent).

8. **Surface errors transparently.** If dispatch or resume fails at
   the platform level (e.g., SendMessage to an expired agent_id
   returns failure), return a `DISPATCH_RESULT` block with `status:
   failed` and the platform error message in `result`. Do NOT
   silently retry or fall back without explicit instruction.

---

## Examples

### Example 1: Sub-agent A needs a literature pull, then continues

**A returns:**

```
[A's normal report content...]

--- DISPATCH_REQUEST ---
agent: literature-researcher
brief: |
  Find and extract the exact closed-form expression for
  the multi-region critical-sphere k_eff from Williams 2005
  Annals of Nuclear Energy 32, p. 215. PDF likely paywalled;
  if so, return the citation graph + best-found alternative.
  Output format:
    - equation number + verbatim transcription
    - ORPHEUS notation map
    - V&V pillar classification (closed-form / MMS / semi-analytical)
wait_for: equation + notation map + pillar
followup: true
--- END DISPATCH_REQUEST ---
```

**Main agent** dispatches `literature-researcher`. The Agent tool
returns the result with metadata including `agentId: lit_a1b2c3`.

**Main agent resumes A:**

```
--- DISPATCH_RESULT ---
agent: literature-researcher
agent_id: lit_a1b2c3
status: completed
result: |
  [literature-researcher's full return verbatim — equations,
   notation map, pillar classification, etc.]
--- END DISPATCH_RESULT ---
```

A captures `lit_a1b2c3` for any follow-up.

### Example 2: Sub-agent A wants follow-up from the SAME literature-researcher

**A returns** (after some further work):

```
--- DISPATCH_REQUEST ---
agent: literature-researcher
agent_id: lit_a1b2c3
brief: |
  Follow-up to your earlier extraction of Williams 2005:
  was the kernel in the equation you transcribed defined for
  vacuum BC, white BC, or both? If both, please separate the
  two cases and indicate which the multi-region tabulated
  k_eff values use.
wait_for: BC type per equation + which BC the table values use
followup: true
--- END DISPATCH_REQUEST ---
```

**Main agent** sees `agent_id: lit_a1b2c3` → uses SendMessage with
`to: lit_a1b2c3` and the brief as message. The literature-researcher
instance still has its prior context (knows what it already
extracted), so the answer is fast and contextful. Result returned
to A wrapped in a `DISPATCH_RESULT` block with the same `agent_id:
lit_a1b2c3`.

### Example 3: Sub-agent A delegates and exits (fire-and-forget)

```
--- DISPATCH_REQUEST ---
agent: archivist
brief: |
  Now that the prototype is shipped + verified (53 tests green,
  see commits aaaa..bbbb), please write the Sphinx theory page
  for the new functionality. Source artifacts:
    - Code: /workspaces/.../greens_function.py
    - SymPy: /workspaces/.../origins/specular/greens_function.py
    - Tests: tests/derivations/test_..._symbolic.py + ...
    - Closeout memo: .claude/agent-memory/...
  Output: Sphinx page docs/theory/peierls_greens.rst extended
  with the new sections; clean -W build.
wait_for: clean Sphinx build + section list
followup: false
--- END DISPATCH_REQUEST ---
```

**Main agent** dispatches the archivist; result goes to the user
(not back to A). A is done.

### Example 4: Sub-agent A wants two independent dispatches in parallel

```
--- DISPATCH_REQUEST ---
agent: literature-researcher
brief: |
  Extract Garcia 2021 J. Comp. Phys. 433 fixed-source benchmark
  table 5 numerical values + convention map.
wait_for: table values + convention notes
followup: true
--- END DISPATCH_REQUEST ---

--- DISPATCH_REQUEST ---
agent: cross-domain-attacker
brief: |
  Examine the multi-region operator I just sketched (path:
  /workspaces/.../greens_function.py § solve_..._mr) and
  check whether the discretisation choice matches the problem's
  native structure, or whether a foreign mathematical frame
  (BEM, spectral, etc.) would be cleaner.
wait_for: structural assessment + alternative frames (if any)
followup: true
--- END DISPATCH_REQUEST ---
```

**Main agent** dispatches both in one message (parallel Agent calls),
captures both agent_ids, returns BOTH `DISPATCH_RESULT` blocks in the
single resume message to A.

### Example 5: Main-agent augmentation in flight

Sub-agent A (numerics-investigator on a wrong-answer cascade)
returns:

```
[A's investigation report so far...]

--- DISPATCH_REQUEST ---
agent: literature-researcher
brief: |
  Find the closed-form expression for the cylinder Bickley-
  Naylor 3-D Ki_n integral series in Knyazev & Selivanov 2014.
  Specifically need Ki_2, Ki_3, Ki_4. Output: closed form +
  recursion relation + ORPHEUS notation map.
wait_for: closed forms + recursion
followup: true
--- END DISPATCH_REQUEST ---
```

Main agent has BREADTH the sub-agent lacks: it also has, from a
parallel session, the user's note that ORPHEUS is considering a
Plan-(c) cylinder Variant α extension. Main agent appends an
augmentation block before dispatch:

**Brief sent to literature-researcher** (with augmentation):

```
--- ORIGINAL_BRIEF ---
Find the closed-form expression for the cylinder Bickley-
Naylor 3-D Ki_n integral series in Knyazev & Selivanov 2014.
Specifically need Ki_2, Ki_3, Ki_4. Output: closed form +
recursion relation + ORPHEUS notation map.
--- END ORIGINAL_BRIEF ---

--- MAIN_AGENT_AUGMENTATION ---
While you have Knyazev open, please ALSO note: does the paper
provide the analogous expressions for the spherical analog
(if any), and does it give numerical-stability guidance for
high-n_quad evaluation? This information would inform a
parallel Plan-(c) cylinder Green's-function extension.
Answer separately in AUGMENTATION_RESPONSE — do not mix with
the original task.
--- END MAIN_AGENT_AUGMENTATION ---
```

literature-researcher returns:

```
--- ORIGINAL_RESPONSE ---
[Closed forms for Ki_2/Ki_3/Ki_4 + recursion + notation map,
 answering A's request]
--- END ORIGINAL_RESPONSE ---

--- AUGMENTATION_RESPONSE ---
Knyazev does not give explicit spherical analogs — the paper
is cylinder-only. For spherical Ki_n analogs, see Carlvik 1965
which the bibliography cites at Ref. 12. Numerical stability:
for high n_quad (>96) use mpmath with dps≥30; the recursion
becomes ill-conditioned via Σ-summation but stable via the
Bessel-K reformulation given in §4.2.
--- END AUGMENTATION_RESPONSE ---
```

Main agent resumes A with the FULL result (both responses) inside
`DISPATCH_RESULT`:

```
--- DISPATCH_RESULT ---
agent: literature-researcher
agent_id: lit_xy789
status: completed
result: |
  --- ORIGINAL_RESPONSE ---
  [Ki_n closed forms + recursion + notation map]
  --- END ORIGINAL_RESPONSE ---

  --- AUGMENTATION_RESPONSE ---
  [spherical analogs + numerical stability notes]
  --- END AUGMENTATION_RESPONSE ---
--- END DISPATCH_RESULT ---
```

A reads `ORIGINAL_RESPONSE` for its task. A may glance at the
`AUGMENTATION_RESPONSE` to know what extra context the main agent
extracted, but A does not need to act on it — that's the main
agent's lateral concern, not A's investigation.

---

## Anti-patterns

### For sub-agents

- **NEVER attempt to call the Agent tool directly.** It's not
  available inside a sub-agent's context; the call will fail. Use
  `DISPATCH_REQUEST` instead.
- **NEVER write a brief that assumes the dispatched agent has your
  context.** Every brief must be self-contained — file paths, prior
  results, exact ask, expected output format. The dispatched agent
  starts with zero history.
- **NEVER omit `agent_id`** when you specifically want to address
  the SAME instance you got an ID from earlier. Without it the main
  agent dispatches fresh, doubles the token cost, and loses the
  prior context.
- **NEVER include a `DISPATCH_REQUEST` in your prose body** without
  the `--- DISPATCH_REQUEST ---` and `--- END DISPATCH_REQUEST ---`
  delimiters. The main agent parses by delimiter; informal mentions
  ("the main agent should also dispatch X") will be missed.
- **NEVER request to dispatch yourself** (the same agent type as you
  are) just to "wait" for something. If you need to wait for an
  external event, return to the user. If you need to coordinate
  with another instance of your own type, that's a fresh dispatch
  by another agent type, not a self-dispatch.
- **NEVER request a dispatch as a way to avoid doing your own job.**
  Dispatch only when you genuinely need a different specialist's
  context, tools, or perspective. Dispatch is not a shortcut for
  effort.

### For the main agent

- **NEVER modify the requesting sub-agent's brief content** before
  dispatching. Augmentation is APPEND-only via a separate
  `MAIN_AGENT_AUGMENTATION` sub-block; the original brief is
  verbatim.
- **NEVER omit the dispatched agent's `agent_id`** when resuming the
  requester. The ID is the protocol's load-bearing piece — it
  enables follow-up.
- **NEVER combine multiple results into a single `DISPATCH_RESULT`**.
  One block per dispatched agent.
- **NEVER silently retry** a failed dispatch or resume. Surface the
  failure transparently in `status: failed`; let the requester
  decide.
- **NEVER act on a `DISPATCH_REQUEST` block without the proper
  delimiters.** If a sub-agent mentions dispatching in prose without
  the delimiter format, ask them to re-emit using the protocol
  format rather than guessing.
- **NEVER augment with questions duplicative of the requester's
  brief**. The augmentation must use your lateral breadth; rewording
  the requester's questions just confuses the dispatched agent and
  wastes tokens.
- **NEVER strip the `AUGMENTATION_RESPONSE`** from the result before
  returning to the requester. The requester benefits from seeing
  what extra context you extracted — they may notice an
  inconsistency you missed.

---

## Why the `agent_id` design matters

Without resume-by-ID, every cross-sub-agent question would force a
fresh Agent dispatch, re-loading the agent's tools, system prompt,
and any context it had built. For an agent like `literature-researcher`
that may have spent a session pulling and indexing 5 PDFs, losing
that context to a fresh dispatch wastes ~30k tokens and ~2 minutes
on every follow-up.

By exposing the `agent_id` and routing follow-ups through SendMessage
resume, the protocol lets sub-agents amortize the setup cost across
multiple questions to the same specialist instance. The cost is one
extra field in two block types; the benefit is order-of-magnitude
fewer redundant dispatches.

---

## Pointers

- **Constraint source**: Anthropic Claude Code docs,
  `https://code.claude.com/docs/en/sub-agents` — under the
  Plan-subagent description: _"This prevents infinite nesting
  (subagents cannot spawn other subagents) while still gathering
  necessary context."_
- **Main-agent reciprocal procedure**: should be added to project
  `CLAUDE.md` as the orchestration glue. Every sub-agent's
  `AGENT.md` should preload this skill.
- **Adjacent skills**:
  - `vv-principles` — when V&V verification is needed before
    declaring a sub-agent's deliverable done.
  - `behavioral-auto-regression` — for diagnosing a sub-agent that
    keeps trying to call Agent directly instead of using the
    DISPATCH_REQUEST protocol.
- **Existing usage** (after adoption): every project AGENT.md should
  list this skill in its `skills:` frontmatter. Search:
  `grep -l 'subagent-handoff-protocol' .claude/agents/`.
