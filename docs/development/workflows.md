# Workflows — phases, briefs, and the return contract

The always-on core ([rules/workflows](rules/workflows.md)) carries the roles,
the four invariants and the seven workflows as one-liners. This page carries
what an agent needs when it is inside one: the phases in detail, the brief
template, and the return contract. It replaces the retired
`subagent-handoff-protocol` skill, which was written for a harness in which a
sub-agent could not dispatch another; since at least Claude Code 2.1.259 (measured; see History) a sub-agent
can, up to three layers deep, so inter-agent coordination is a tool grant, not
a relay protocol. The relay's block formats survive only as the `NEEDS:` return
contract below.

## W1 — Build a capability

| phase | key agent | supports | gate | on failure |
|---|---|---|---|---|
| P0 context | orchestrator | explorer; literature-researcher when the formulation is published | the theory page's Key Facts read; the Nexus briefing run | — |
| P1 verification design | **test-architect** | explorer | a spec whose every gate names the input in today's tree that it rejects (its first red) | redesign before P2 |
| P2 build | **method-implementer**, or the main agent for a surgical carve | explorer, literature-researcher, numerics-investigator when a probe is needed; cross-domain-attacker after the first pass | scope suite green; the spec's gates land with the code | the implementer fixes; a second probe is a new dispatch, not a re-brief |
| P3 review, in parallel | **qa** and **elegance-enforcer**, dispatched by the orchestrator on the artefact | explorer | qa: term-level correctness, coverage, mutation; elegance: structure, three-leg verdicts | any red resumes the implementer by name with the findings |
| P4 documentation | **archivist** | explorer | theory page, derivations, changelog; `sphinx -W` clean; `dead_references` 0 | archivist iterates |
| P5 close-out | orchestrator | — | issues closed or filed; retirement audit run; one commit per landed unit | — |

## W2 — Wrong answer

**numerics-investigator** runs the probe cascade (drop one complication at a
time to the minimal reproducer); it may call test-architect for the permanent
test and literature-researcher for the reference formulation. The fix lands
(implementer or main agent). **qa** reviews the fix with a mutation that
re-introduces the defect and requires a red. **archivist** writes the ERR entry
and the theory-page note. Close-out as W1.

## W3 — Surgical carve

The main agent writes, with the user steering step by step; `method-implementer`
is not dispatched. explorer enumerates the blast set by contract (every site
that depends on the interface, most of which spell no symbol). **test-architect**
designs the gates and re-baselines per `vv-principles`' three criteria. Review
as W1-P3; **archivist** for the changelog.

## W4 — Documentation campaign

**archivist** is the key agent; qa verifies claims against the tree; explorer
answers structural questions. Gates: `sphinx -W`, `dead_references`, `staleness`.

## W5 — Design review

**cross-domain-attacker** (does the formulation match the problem's native
structure?) and **elegance-enforcer** (does the design meet the elegance
standard?) on a first-pass design, in parallel; the output feeds W1-P2.

## W6 — Tree-wide census

`general-purpose` haiku categorisers, each with a fixed output schema and a
bounded slice of the tree; the orchestrator aggregates. Low-freedom briefs: the
schema, the predicate, the slice, the positive control. No nesting.

## W7 — Literature acquisition

**literature-researcher**: `scratch/literature/` first, then the OCR sidecars
(`scratch/literature_ocr/`), then online. "Not in the local folder" is a
question to the user, never a unilateral pivot to a secondary source.

## The brief

Every brief carries, in this order:

```
Workflow: W<n> — phase P<k>
Artefacts from the previous phase: <paths>
Ask: <one paragraph, self-contained; assume zero context from this session>
Rules that apply to you: <only for Support agents, which see no project rules:
  the two or three that matter — e.g. the ugrep silent-zero hazard for a
  census; local-folder-first for literature>
Return contract: report under <N> words; the file(s) at <paths> carry the
  detail; end with a NEEDS: block (see below), empty if nothing is missing.
```

For an Opus-pinned agent the word cap is not optional: Opus 5 writes longer
responses and longer files by default, and every unbounded report lands in the
parent's context. Never pre-filter a reviewer's severity ("only report
high-severity issues" is followed literally and under-reports); ask for
everything and filter at the parent.

## The return contract

```
<report, under the cap>

NEEDS:
- <thing the agent could not obtain, and from whom: a file, a ruling, a
   dispatch it lacks the tool for> — or "none"
```

The parent acts on `NEEDS:` and resumes the agent by name with the answer. This
is the whole of what remains of the former `DISPATCH_REQUEST` /
`DISPATCH_RESULT` relay: an agent that holds `Agent` dispatches for itself; an
agent that does not (Support) asks in `NEEDS:`.

## History

Until 2026-09 the `subagent-handoff-protocol` skill (retired 2026-09-20) carried
a relay protocol built on the belief that "subagents cannot spawn other
subagents" was an Anthropic platform constraint. Measured 2026-09-05 and again
2026-09-20 (Claude Code 2.1.259 and 2.1.278): a depth-1 agent spawned a depth-2
agent which spawned a depth-3 agent that answered; the depth-3 agent had no
`Agent` tool, exactly the documented cap. The constraint was this project's own
`tools:` allowlists. The plan that made the change is
`.claude/plans/harness_context_budget.md` (Groups E and H).
