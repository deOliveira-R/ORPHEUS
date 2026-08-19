---
name: nexus-exploring
description: "Use when the user asks how code works, wants to understand architecture, trace execution flows, or explore unfamiliar parts of the codebase. Examples: \"How does X work?\", \"What calls this function?\", \"Show me the auth flow\", \"How does this equation connect to code?\""
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

## Symptom → tool

⛔ **This is NOT the instrument-choice table.** *Nexus vs grep vs Read* is
settled by the always-loaded `.claude/rules/nexus-tools.md` "Question → Tool"
table — do not re-derive it here, and do not restate it. This table covers the
two families that one does not: the **architecture smells** and the **runtime
overlay**, plus the position bridge — and for every family, the caveat that
makes the answer honest.

⭐ Consult it before hand-rolling a multi-query search. The third column is the
load-bearing one: skipping it is how a confident wrong answer gets reported.

### I have a position, not a name (the bridge into the graph)

| you have | tool | the caveat |
|---|---|---|
| `file:line` from a traceback, LSP result, or diff | `node_at(file, line)` | every other tool takes a node id; this is the only way in from a position |
| only a path | `file_brief(file)` | `null` means the file is not in the graph at all — new, excluded, or a stale build |

### What smells? (the missing-abstraction family — no rows for this in the rule)

| smell | tool | the caveat |
|---|---|---|
| the same computation, twice | `twin_paths` | symmetric-by-design pairs (`apply`/`apply_transpose`) legitimately resemble each other |
| a repeated conditional on a tag | `discriminations` | *a repeated conditional is a missing type* — but a genuinely open set may stay a tag |
| a function that belongs in a class | `native_place` | a pure, independently-tested free function is fine as-is; weigh `excluded_callers` |
| unreachable code | `dead_functions` | ⛔ **candidates, not verdicts.** `unresolved_calls > 0` means it is probably called and the resolver lost the edge |
| structural hotspots | `god_nodes`, `bridges`, `communities` | `god_nodes` excludes stdlib by default, or `float`/`int` win |
| undeclared Protocol conformance | `protocol_conformers` | matches method NAMES, not signatures — a lead, not a proof |

### What actually RAN? (the runtime overlay — the only relation that can REFUTE)

| question | tool | the caveat |
|---|---|---|
| which tests executed this node | `runtime_exercisers(run, node)` | ⛔ needs a capture with `dynamic_context = test_function` **and an unbroken `__init__.py` chain from the rootdir**. `[M]` 2026-08-19: one missing `tests/sn/__init__.py` made the whole SN tree unattributable — `exercised_by` **0**, silently, with `rc=0` and a plausible bind count. Check `exercised_by > 0` before trusting any capture |
| dispatch the static graph cannot see | `runtime_edges(mode="dynamic_only")` | annotation-mediated dispatch, and the resolved face of polymorphism |
| a conditional never taken both ways | `runtime_branches` | the dynamic twin of `discriminations` |
| a recompute smell | `runtime_hotspots(by="ncalls")` | a property called 10 000× is a caching opportunity |
| the observed stage ORDER | `runtime_timeline` | `viztracer` only — cProfile has counts, not order |
| which marker pytest RESOLVED | `runtime_markers` | ⭐ a module-level `pytestmark` is **invisible to AST extraction**; only a `--collect-only` capture sees it |
| what captures exist | `runtime_runs` | read `families` — a run without `exercised_by` cannot adjudicate a claim |

### The caveats on the rule's OWN rows

The rule routes these correctly; what it does not say is what the answer omits.

| the rule sends you to | what the answer does not contain |
|---|---|
| `impact` / `callers` for blast radius | ⛔ the cone follows `calls`/`type_uses`/`inherits`. `[M]` 2026-08-18: **12–15 % recall** against execution, and **0 of 300** proven test↔symbol pairs have ANY path over it — properties, dunders, callbacks and polymorphic dispatch mint no edge. And `callers` **empty ≠ dead**: read the `unresolved` block |
| `provenance_chain` for traceability | ⚠ an entry marked `inferred` was minted from a shared name token; `via` names the tokens |
| `verification_audit` for coverage | ⛔ read `code_evidence` first. `[M]` 2026-08-19: **13 121 of 13 508** `implements` edges are still guesses. With `run=`, split verdicts on `TestReference.source == "declared"` — the aggregate sums authored claims with two heuristic BFS tiers |

⭐ And the one the rule does not mention: `retest(scope, run=…)`. With a capture,
each row carries `warrant: executed` — a capture PROVES it ran. Without one it
is the same weak cone as `impact`, labelled `reachable`.

### Am I answering from the right checkout?

`workspaces` lists every checkout and its graph; `use_workspace(root)` switches.
⛔ Spawned inside `.claude/worktrees/*`, every query answers from the MAIN
checkout until you switch. Silence from `session_briefing` means the INDEXED
SOURCES match — not that the branch does.

See [reference.md](reference.md) for schemas, node-id grammar, edge types, the
`graph_query` pattern syntax, and the CLI.
