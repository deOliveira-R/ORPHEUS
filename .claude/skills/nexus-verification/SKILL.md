---
name: nexus-verification
description: "Use when the user wants to check verification status, documentation coverage, or doc-code drift. Examples: \"What's verified?\", \"Which docs are stale?\", \"What equations have no tests?\", \"Documentation coverage report\""
---

# Verification & Documentation Quality with Nexus

IMPORTANT: This skill is the dedicated tool for V&V assessment and
documentation drift detection. It replaces Grep for all coverage questions.

## Workflow

```
1. verification_coverage()                                → Full V&V map
2. verification_coverage({status_filter: "implemented"})  → Gaps: code but no tests
3. staleness()                                             → Docs that drifted
4. session_briefing()                                      → Combined overview
```

## Coverage Status Values

| Status | Meaning | Action |
|--------|---------|--------|
| **verified** | Equation + code + test | Fully traced |
| **tested** | Code + test, no equation link | Add IMPLEMENTS doc |
| **implemented** | Equation + code, no test | Write a test |
| **documented** | Equation only, no code | Implement or mark future |
| **orphan_code** | Code with no equation | Document the theory |

## Test Inventory Queries

Tests are indexed in the graph. Test node IDs: `py:function:tests.<file>.<function>`.

```
# List all test modules
query({text: "tests.test_", node_types: "module", limit: 50})

# List test functions in a module
neighbors({node_id: "py:module:tests.test_sn_1d", direction: "out", edge_types: "contains"})

# Find tests covering a function
impact({target: "py:function:orpheus.sn.solver.solve_sn", direction: "upstream"})

# Trace test → equations
trace_error({test_node_id: "py:function:tests.test_cp_slab.test_slab_cp_eigenvalue"})
```

## Checklist

- [ ] `verification_coverage` for full V&V status
- [ ] Review "implemented" entries (verification gaps)
- [ ] Review "documented" entries (unimplemented equations)
- [ ] Review "orphan_code" entries (undocumented theory)
- [ ] `staleness` to find docs needing updates
- [ ] Create GitHub Issues for each gap found

See [../nexus-exploring/reference.md](../nexus-exploring/reference.md) for full reference.
