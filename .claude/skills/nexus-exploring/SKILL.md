---
name: nexus-exploring
description: "Use when the user asks how code works, wants to understand architecture, trace execution flows, or explore unfamiliar parts of the codebase. Examples: \"How does X work?\", \"What calls this function?\", \"Show me the auth flow\", \"How does this equation connect to code?\""
---

# Exploring with Nexus

## When to Use

- "How does the SN solver work?"
- "What calls this function?"
- "Show me the main components"
- "How does this equation connect to code?"
- "What's the mathematical provenance of this function?"
- Understanding code or theory you haven't seen before

## Workflow

```
1. READ nexus://graph/stats                          → Graph overview
2. nexus query({text: "<concept>"})                   → Find related symbols
3. nexus context({node_id: "<symbol>"})               → 360-degree view
4. nexus shortest_path({source: "<A>", target: "<B>"}) → How concepts connect
5. nexus provenance_chain({node_id: "<symbol>"})       → Citation → equation → code chain
```

## Checklist

```
- [ ] nexus stats() — overview of graph size and types
- [ ] nexus query() for the concept you want to understand
- [ ] nexus context() on key symbols for callers/callees/docs
- [ ] nexus provenance_chain() to trace mathematical origins
- [ ] nexus communities() to see functional groupings
- [ ] Read source files for implementation details
```

## Tools

**query** — find symbols by keyword:
```
nexus query({text: "collision probability"})
→ Functions, classes, equations matching the search, sorted by connectivity
```

**context** — 360-degree view of a symbol:
```
nexus context({node_id: "py:function:sn_solver.solve_sn"})
→ All incoming/outgoing edges grouped by type (calls, documents, implements, etc.)
```

**provenance_chain** — mathematical traceability:
```
nexus provenance_chain({node_id: "py:function:sn_sweep.sweep_spherical"})
→ Bailey2009 → Eq.alpha-recursion → sweep_spherical → helper functions
```

**shortest_path** — how concepts connect:
```
nexus shortest_path({source: "doc:theory/collision_probability", target: "py:class:numpy.ndarray"})
→ Theory page → function → numpy dependency
```

## Resources

| Resource | What you get |
|----------|-------------|
| `nexus://graph/stats` | Node/edge counts by type, density |
| `nexus://graph/communities` | Functional areas with top members |
| `nexus://graph/schema` | Available node types, edge types, ID format |
| `nexus://briefing` | Session briefing: stale docs, coverage gaps, recent changes |

## Node ID Format

All node IDs follow `<domain>:<type>:<qualified_name>`:
- `py:function:sn_solver.solve_sn`
- `py:class:collision_probability.CPMesh`
- `math:equation:alpha-recursion`
- `doc:theory/discrete_ordinates`
- `std:label:theory-collision-probability`
