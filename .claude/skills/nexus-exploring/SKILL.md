---
name: nexus-exploring
description: "Use when the user asks how code works, wants to understand architecture, trace execution flows, or explore unfamiliar parts of the codebase. Examples: \"How does X work?\", \"What calls this function?\", \"Show me the auth flow\", \"How does this equation connect to code?\""
---

# Exploring with Nexus

IMPORTANT: This skill is the dedicated tool for code exploration. It
replaces Grep for all architecture, dependency, and flow questions.

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

See [reference.md](reference.md) for full tool/schema/CLI reference.
