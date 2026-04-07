---
name: nexus-debugging
description: "Use when the user is debugging a bug, tracing an error, or asking why something fails. Especially powerful for numerical errors — traces from failing test to the equations that might be wrong. Examples: \"Why is X failing?\", \"Which equation is wrong?\", \"Trace this bug\""
---

# Debugging with Nexus

## When to Use

- "Why is this test failing?"
- "Which equation might be wrong?"
- "Trace where this error comes from"
- "keff is wrong — what should I check first?"
- Investigating numerical errors, bugs, or unexpected behavior

## Workflow

```
1. nexus query({text: "<error or symptom>"})                    → Find related symbols
2. nexus trace_error({test_node_id: "py:function:test_X"})       → Trace test → equations
3. nexus context({node_id: "<suspect function>"})                → See all callers/callees
4. nexus provenance_chain({node_id: "<suspect>"})                → Get the citation trail
5. Read source + theory page for the flagged equation
```

## Checklist

```
- [ ] Identify the failing test or symptom
- [ ] nexus trace_error() to find equations on the call path
- [ ] Review equations ranked by centrality (most likely root cause first)
- [ ] nexus provenance_chain() to find the original citation
- [ ] Read the theory page for the flagged equation
- [ ] Read the implementing code
- [ ] Compare code against the equation in the citation
```

## The Key Insight

When a numerical test fails, the question isn't "what function is broken" (the
debugger shows that) but "which equation is wrong or misimplemented." Nexus is
the only tool that can trace from a test → through the call graph → to the
equations on that path → to the citations those equations come from.

## Tools

**trace_error** — the primary debugging tool:
```
nexus trace_error({test_node_id: "py:function:test_cp_slab.test_slab_2eg_2rg"})
→ Call chain: solve_cp → CPSolver.power_iteration → CPMesh.compute_pinf_group
→ Equations on path: collision-rate, e3-def, reciprocity
→ Citations: Hebert2009, Stamm1983
→ Check collision-rate first (highest centrality)
```

**provenance_chain** — citation trail for a suspect function:
```
nexus provenance_chain({node_id: "py:function:sn_sweep.sweep_spherical"})
→ Implements: math:equation:alpha-recursion
→ From: Bailey2009 (Eq. 50)
→ Read: docs/theory/discrete_ordinates.rst
```

**context** — 360-degree view to understand call structure:
```
nexus context({node_id: "py:function:sweep.sweep_cylindrical"})
→ Called by: solve_sn, test_sn_cylindrical
→ Calls: compute_alpha_dome, numpy.cumsum
→ Implements: transport-cylindrical equation
```
