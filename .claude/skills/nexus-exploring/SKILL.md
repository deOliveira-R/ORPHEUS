---
name: nexus-exploring
description: "Use when the user asks how code works, wants to understand architecture, trace execution flows, explore unfamiliar parts of the codebase, or hunt structural smells. Examples: \"How does X work?\", \"What calls this function?\", \"Show me the auth flow\", \"How does this equation connect to code?\", \"Find dead code / copy-paste twins / missing abstractions\", \"Which classes share an implicit interface?\""
---

# Exploring with Nexus

IMPORTANT: This skill is the dedicated tool for code exploration. It
complements Grep — use Nexus for structural queries, Grep for text search.

## Workflow

```
1. query({text: "<concept>"})                        → Find symbols
2. context({node_id: "<symbol>"})                    → 360-degree view
3. provenance_chain({node_id: "<symbol>"})            → Citation → equation → code
4. shortest_path({source: "<A>", target: "<B>"})      → How concepts connect
5. Read source files for implementation details
```

## Checklist

- [ ] `query` for the concept you want to understand
- [ ] `context` on key symbols for callers/callees/docs
- [ ] `provenance_chain` to trace mathematical origins
- [ ] `communities` to see functional groupings
- [ ] Read source files for implementation details

## Key Tools

**query** — find symbols by keyword:
```
query({text: "collision probability"})
→ Functions, classes, equations matching the search, sorted by connectivity
```

**context** — 360-degree view of a symbol:
```
context({node_id: "py:function:orpheus.sn.solver.solve_sn"})
→ All incoming/outgoing edges grouped by type
```

**provenance_chain** — mathematical traceability:
```
provenance_chain({node_id: "py:function:orpheus.sn.sweep.transport_sweep"})
→ Bailey2009 → Eq.transport-cartesian → transport_sweep
```

**shortest_path** — how concepts connect:
```
shortest_path({source: "doc:theory/collision_probability", target: "py:class:numpy.ndarray"})
→ Theory page → function → numpy dependency
```

## Symptom → tool (don't hand-roll these with query/grep)

Each of these questions has a DEDICATED tool that answers in one call
what generic `query`/`context` exploration approximates in ten:

| The user suspects… | Call |
|---|---|
| dead code / "does anything call this?" at scale | `dead_functions` |
| copy-paste twins, parallel implementations that drifted | `twin_paths` |
| the same string/enum tag branched on everywhere (missing type) | `discriminations` |
| classes sharing an implicit interface without a base class | `protocol_conformers` |
| a helper living in the wrong module (feature envy) | `native_place` |
| docs/docstrings referencing deleted symbols or equations | `dead_references` |
| "map the functional areas" / load-bearing hubs | `communities`, `god_nodes` |
| the few nodes holding separate areas together | `bridges` |
| what actually RAN (hotspots, fired-vs-static edges, coverage) | `runtime_runs` → `runtime_hotspots` / `runtime_edges` / `runtime_branches` |
| which TESTS executed a node (the only relation that can REFUTE a claim) | `runtime_exercisers(run, node)` |
| which marker pytest actually RESOLVED at collection | `runtime_markers` |

### ⛔ The caveat on each answer — skipping it is how a confident wrong answer ships

| the tool | what its answer does NOT contain |
|---|---|
| `impact` / `callers` | the cone follows `calls`/`type_uses`/`inherits`. `[M]` 2026-08-18 on a real corpus: **12–15 % recall** against execution evidence, and **0 of 300** proven test↔symbol pairs have ANY path over it — properties, dunders, callbacks and polymorphic dispatch mint no edge. And `callers` **empty ≠ dead**: read the `unresolved` block |
| `dead_functions` | **candidates, not verdicts.** `unresolved_calls > 0` means it is probably called and the resolver lost the edge |
| `provenance_chain` | an entry marked `inferred` was minted from a shared name token; `via` names the tokens |
| `verification_coverage` / `verification_audit` | read `code_evidence` FIRST. `[M]` 2026-08-19: **13 121 of 13 508** `implements` edges on one corpus are still guesses. With `run=`, split verdicts on `TestReference.source == "declared"` — the aggregate sums authored claims with two heuristic BFS tiers, and it INFLATES as declarations land |
| `runtime_exercisers` and every `run=` consumer | needs a capture with `dynamic_context = test_function` **and an unbroken `__init__.py` chain from the rootdir**. `[M]` 2026-08-19: one missing `__init__.py` made a 3330-test tree unattributable — `exercised_by` **0**, silently, with `rc=0` and a plausible bind count. **Check `exercised_by > 0` before trusting any capture** |
| `retest` | with `run=` each row carries `warrant: executed` — a capture PROVES it ran. Without one it is the same weak cone as `impact`, labelled `reachable` |
| `errors` | `total_entries: 0` means nothing is DECLARED, not that the project is clean |

## What the user actually says

Users describe **symptoms**, not tools. Measured: phrasings like these
reach the right tool far less often than the table above, because
nothing connects the complaint to the capability. Translate:

| What they say | What it means | Call |
|---|---|---|
| "we keep having to change several classes in lockstep and we forget one" | one concept spread across classes | `protocol_conformers`, `twin_paths`, `discriminations` |
| "two people built these independently and didn't talk" | duplicated computation | `twin_paths` |
| "this module feels bigger than it should be" | dead weight | `dead_functions` |
| "things live in surprising places" / "boundaries got fuzzy" | feature envy | `native_place`, `bridges` |
| "adding a case means touching five files" | missing type | `discriminations` |
| "I just deleted / renamed / moved X — what's left?" | orphaned references | `dead_references`, `callers`, `retest` |
| "the docs feel out of date" | two different drifts | `dead_references` (gone) + `staleness` (old) |

## Sweeps you run WITHOUT being asked

Some requests are a job, not a question, and the job includes checks the
user did not think to ask for. Do these unprompted:

- **"Health check this codebase" / onboarding / "what will confuse
  people"** → the whole smell family in one pass: `twin_paths`,
  `discriminations`, `native_place`, `protocol_conformers`,
  `dead_functions`, plus `dead_references` for the docs. Then
  `communities` / `god_nodes` for the map. A health check that reports
  only structure and misses the smells is not a health check.
- **"We're cutting a release" / "is the doc set publishable?"** →
  `dead_references` FIRST (a dead reference draws no build warning at
  any severity, so nothing else will catch it), then `staleness`, then
  `verification_audit`.
- **After any deletion or rename you performed** → `dead_references`
  before declaring the work done. Green tests do not cover prose.

If the project ships the `/doc-health` command or the dead-references
hook, that finding may already be in your context — act on it rather
than re-deriving it.

## The position bridge (LSP ↔ graph)

The language server and the graph are complementary: LSP resolves
precisely (definitions, references, live unsaved state, alias-aware);
the graph sees what LSP cannot (equations, tests, doc pages, V&V
chains). Bridge in BOTH directions:

- **Position → node**: any (file, line) — an LSP result, a stack
  trace, an editor cursor — feeds `node_at({file, line})` and returns
  the innermost enclosing graph node. Continue with `context` /
  `impact` / `provenance_chain` / `callers` for the cross-domain
  picture.
- **Node → position**: every AST-derived node result carries
  `file_path` + `lineno` — feed them straight to Read / the editor /
  an LSP request. No text-search round-trip.
- **Honesty**: the graph is a build-time snapshot. `node_at` WARNS
  when the queried file changed since the graph's stamped commit —
  positions may then map to the wrong symbol. Rebuild
  (sphinx-build / nexus analyze) and re-ask; do not trust a warned
  mapping for surgical work.

## Worktrees and ambient context

- A session inside a git worktree must query THAT worktree's graph:
  auto-alignment handles sessions launched there; after a mid-session
  EnterWorktree, call `use_workspace(<worktree name>)`. The
  `session_briefing` workspace block warns on branch mismatch.
- Language-server identity errors mentioning `.claude.worktrees.*`
  module paths are wrong-rooted-server noise, not code bugs —
  discount them.
- Projects may wire `nexus file-brief` into an edit-time hook: a few
  lines of graph context (callers, equations, tests, docs) appear
  automatically after you edit a file. The node IDs in a brief are
  copy-pasteable entry points — follow the hub ID with `context` or
  `impact` before reshaping a file. (The graph is always a snapshot:
  after substantial edits, rebuild before trusting positions.)

See [reference.md](reference.md) for full tool/schema/CLI reference.
