---
name: nexus-debugging
description: "Use when the user is debugging a bug, tracing an error, or asking why something fails. Especially powerful for numerical errors — traces from failing test to the equations that might be wrong. Examples: \"Why is X failing?\", \"Which equation is wrong?\", \"Trace this bug\""
---

# Debugging with Nexus

IMPORTANT: This skill is the dedicated tool for bug tracing and
equation diagnosis. It replaces Grep for all debugging exploration.

## The Key Insight

When a numerical test fails, the question isn't "what function is broken"
but "which equation is wrong or misimplemented." Nexus traces
test → call graph → equations → citations. No other tool can do this.

## Workflow

```
1. trace_error({test_node_id: "py:function:tests.test_X.test_Y"})  → Equations on path
2. provenance_chain({node_id: "<suspect equation>"})                 → Citation trail
3. context({node_id: "<suspect function>"})                          → All callers/callees
4. Read the theory page for the flagged equation
5. Compare code against the equation in the citation
```

## Key Tools

**trace_error** — the primary debugging tool:
```
trace_error({test_node_id: "py:function:tests.test_cp_slab.test_slab_2eg_2rg"})
→ Call chain: solve_cp → CPSolver.power_iteration → CPMesh.compute_pinf_group
→ Equations on path: collision-rate, e3-def, reciprocity
→ Citations: Hebert2009, Stamm1983
```

**provenance_chain** — citation trail:
```
provenance_chain({node_id: "py:function:orpheus.sn.sweep.transport_sweep"})
→ Implements: math:equation:transport-cartesian
→ From: Bailey2009 (Eq. 50)
```

See [../nexus-exploring/reference.md](../nexus-exploring/reference.md) for full reference.
