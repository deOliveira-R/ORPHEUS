# Workflows — phases, briefs, and the return contract

The always-on core ([rules/workflows](rules/workflows.md)) carries the roles,
the four invariants and the seven workflows as one-liners. This page carries
what an agent needs when it is inside one: the phases in detail, the brief
template, and the return contract. It replaces the retired
`subagent-handoff-protocol` skill, which was written for a harness in which a
sub-agent could not dispatch another; since at least Claude Code 2.1.261 (measured; see History) a sub-agent
can, up to three layers deep, so inter-agent coordination is a tool grant, not
a relay protocol. The relay's block formats survive only as the `NEEDS:` return
contract below.

## W1 — Build a capability

| phase | key agent | supports | gate | on failure |
|---|---|---|---|---|
| P0 context | orchestrator | explorer; literature-researcher when the formulation is published | the theory page's Key Facts read; the Nexus briefing run | — |
| P1 verification design | **test-architect** | explorer | a spec whose every gate names the input in today's tree that it rejects (its first red) | redesign before P2 |
| P2 build | **method-implementer**, or the main agent for a surgical carve | explorer, literature-researcher, numerics-investigator when a probe is needed; cross-domain-attacker after the first pass | scope suite green; the spec's gates land with the code | the implementer fixes; a second probe is a new dispatch, not a re-brief |
| P3 review, in parallel | **qa** and **elegance-enforcer**, dispatched by the parent on the artefact | explorer | qa: term-level correctness, coverage, mutation; elegance: structure, every violation with its three legs (what, which pattern, the remedy) | any red resumes the implementer by name with the findings |
| P4 documentation | **archivist** | explorer | theory page, derivations, changelog; `sphinx -W` clean; `dead_references` 0 | archivist iterates |
| P5 close-out | orchestrator | — | issues closed or filed; retirement audit run; one commit per landed unit | — |

## W2 — Wrong answer

**numerics-investigator** runs the probe cascade (drop one complication at a
time to the minimal reproducer); it may call test-architect for the permanent
test and literature-researcher for the reference formulation. The fix lands
(implementer or main agent). **qa** and **elegance-enforcer** review the fix in
parallel, dispatched by the parent: qa with a mutation that re-introduces the
defect and requires a red; the enforcer because a fix is where a patch replaces
the structural repair. **archivist** writes the ERR entry
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
(`scratch/literature_ocr/`), then online. The brief spells out the exact path
and, when the paper is known to be local, is phrased as "extract equations from
the local PDF at `<full path>`" rather than "find and acquire paper X": the
first is unambiguous and cheap. The user maintains the folder actively and has
every Nuclear Science & Engineering volume locally. "Not in the local folder"
is the agent's FIRST response, as a question to the user ("acquire it, or
will you add it?"), never a unilateral pivot to a secondary source: a pivot is a structural decision (a
different math path, possibly a weaker verification claim) that needs the
user's approval, not agent autonomy. The sidecar-first search, the scan as the
source of truth for every load-bearing equation, the paraphrase-and-page-cite
output discipline and the incremental write of the memo (the 2026-07-22
content-filter deaths: two agents killed mid-generation lost everything,
zero bytes written) are the agent's own procedure, carried by its AGENT.md,
not by the brief.

**Zotero liveness: 0 hits is not "empty".** A Zotero MCP server returning 0
hits on a known-present item together with connection-refused on port 23119 is
BROKEN, not a library that lacks the paper; the misreading is silent and burns
turns on repeated 0-hit queries. Stop querying, record "Zotero down — no
annotations checked", and proceed on the local folder and the web tier. This
applies to any agent briefed to consult Zotero, not only literature-researcher.

## The brief

Every brief carries, in this order:

```
Workflow: W<n> — phase P<k>
Artefacts from the previous phase: <paths>
Ask: <one paragraph, self-contained; assume zero context from this session>
Rules that apply to you: <for the three `omitClaudeMd` agents, which see no
  project rule and no project memory index — only their AGENT.md, their own
  agent memory and their preloaded skills — this line is the only place a
  project rule reaches them; a `general-purpose` categoriser inherits the rules
  and needs only its schema. Paste the generated list below, then the
  task-specific items: whether the agent may edit tracked files at all (a
  census is read-only); if another agent is editing the tree meanwhile, what,
  where and until when (L38); any negative you assert about the tree ("X has
  no gate") marked `[R]` for the agent to re-verify (L50). For literature: W7
  above, in full — `scratch/literature/` first, spelled out, then the OCR
  sidecars; "not in the local folder" is a question to the user, never a
  pivot; Zotero at 0 hits plus connection refused is down, not empty. For a
  design review (cross-domain-attacker): the artefact's path; the return is
  structural detection, no critique. If the task depends on Nexus: what to do
  when it is missing (a sub-agent has no `ToolSearch`). Any campaign pointer
  the agent needs, verbatim — it has no project memory index.>
Return contract: report under <N> words; the file(s) at <paths> carry the
  detail; a verification claim pastes the pytest summary line verbatim, in a code fence (L12);
  end with a NEEDS: block (see below), empty if nothing is missing.
```

The generated list, one item per rule that declares a `brief:` in its
`harness:` block (`tools/harness/brief.py` assembles it; a hand edit here is
drift):

<!-- BEGIN GENERATED brief rules — source: the harness.brief of every page under docs/development/rules/; edit the source, not this block -->
- `articulation`: a report is complete sentences; `[M]` measured with its command, `[R]` reasoned, `[HYPOTHESIS]` proposed; a bare number is read as measured.
- `coding-standards`: tests run as `python -O -m pytest`; a bare `assert` outside a collected test module is stripped under `-O`, so a contract is a `raise` and a test-side check is `np.testing.assert_*`.
- `instrument-doctrine`: a zero from a filter is evidence only after a positive control of each shape it must find (X1); every count states its predicate, its tree and its exclusions, a universal is `k of N`, and a completeness claim is re-run in Python (`re` + `pathlib.rglob`) so its denominator is stated (X2).
- `nexus-tools`: `grep` is ugrep, and an anchor inside an alternation group matches nothing, silently: use `\b…\b` or `-P` with a lookbehind; a sub-agent has no `ToolSearch`, so if Nexus is missing say so in NEEDS: and fall back to Bash.
- `process-discipline`: never `git checkout`, `git restore` or `git stash` a path that carries uncommitted edits: they revert to HEAD and destroy the work (L28); revert a mutation by monkeypatching in-process or by mutating a copy.
<!-- END GENERATED brief rules -->

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
2026-09-20 (Claude Code 2.1.261 and 2.1.278): a depth-1 agent spawned a depth-2
agent which spawned a depth-3 agent that answered; the depth-3 agent had no
`Agent` tool, exactly the documented cap. The constraint was this project's own
`tools:` allowlists. The plan that made the change is
`.claude/plans/harness_context_budget.md` (Groups E and H). Measured 2026-09-20
and 2026-09-21 (2.1.278): a finished agent resumes on a message from its parent
with its full history, so a reviewer is resumed by name rather than re-briefed;
a sub-agent holds `SendMessage`, and addresses a named sibling with it, only
when its `tools:` allowlist admits it (the `general-purpose` agent does, `qa`
with its explicit list does not), and a sub-agent with an explicit list holds
no `Skill` tool even when `Skill` is listed, so a project agent's skills are
its `skills:` preload and nothing else.
