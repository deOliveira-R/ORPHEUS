---
name: nexus-verification
description: "Use when the user wants to check verification status, documentation coverage, or doc-code drift. Examples: \"What's verified?\", \"Which docs are stale?\", \"What equations have no tests?\", \"Documentation coverage report\""
---

# Verification & Documentation Quality with Nexus

## When to Use

- "Which equations are verified?"
- "What's the documentation coverage?"
- "Which docs are stale?"
- "What equations have no implementing code?"
- "What code has no tests?"
- V&V status assessment
- Documentation quality review

## Workflow

```
1. nexus verification_coverage()                    → Full equation → code → test map
2. nexus verification_coverage({status_filter: "implemented"})  → Gaps: code but no tests
3. nexus staleness()                                → Docs that drifted from code
4. nexus session_briefing()                         → Combined overview
```

## Checklist

```
- [ ] nexus verification_coverage() for full V&V status
- [ ] Review "implemented" entries — equations with code but no tests (verification gaps)
- [ ] Review "documented" entries — equations with no implementing code
- [ ] Review "orphan_code" entries — code with no equation (undocumented theory)
- [ ] nexus staleness() to find docs needing updates
- [ ] Create GitHub Issues for each gap found
```

## Coverage Status Values

| Status | Meaning | Action |
|--------|---------|--------|
| **verified** | Equation + code + test | Good — fully traced |
| **tested** | Code + test, no equation link | Add IMPLEMENTS documentation |
| **implemented** | Equation + code, no test | Write a test — verification gap |
| **documented** | Equation only, no code | Either implement or mark as future work |
| **orphan_code** | Code with no equation | Document the theory behind it |

## Test Inventory Queries

Tests are indexed in the Nexus graph (via `nexus_extra_source_dirs`).
Test node IDs use the format `py:function:tests.<file>.<function>`,
e.g. `py:function:tests.test_sn_1d.test_homogeneous_exact`.
Test modules use `py:module:tests.<file>`.

**List all test files (modules):**
```
nexus query({text: "tests.test_", node_types: "module", limit: 50})
→ tests.test_sn_1d (degree 12), tests.test_cp_slab (degree 8), ...
```

**List all test functions in a specific test file:**
```
nexus neighbors({node_id: "py:module:tests.test_sn_1d", direction: "out", edge_types: "contains"})
→ test_homogeneous_exact, test_heterogeneous_convergence, test_angular_convergence, ...
```

**List all orpheus modules (to cross-reference with test coverage):**
```
nexus query({text: "orpheus.", node_types: "module", limit: 50})
→ orpheus.sn.solver, orpheus.cp.solver, orpheus.fuel.solver, ...
```

**Find which tests cover a specific solver function:**
```
nexus impact({target: "py:function:orpheus.sn.solver.solve_sn", direction: "upstream"})
→ d=1: tests.test_sn_1d.test_homogeneous_exact, tests.test_sn_solver_components.TestHomogeneousExact, ...
```

**Find which modules have NO tests (untested code):**
```
nexus impact({target: "py:module:orpheus.fuel.solver", direction: "upstream"})
→ If no tests.* nodes appear at any depth → module is untested
```

**Map modules to test coverage (replaces `ls` + `diff` shell commands):**
Run impact upstream on each solver's main function. If `tests.*` nodes
appear, it's tested. If none appear, it's untested:
```
nexus impact({target: "py:function:orpheus.thermal_hydraulics.solver.solve_thermal_hydraulics", direction: "upstream"})
→ No tests.* nodes → untested
nexus impact({target: "py:function:orpheus.sn.solver.solve_sn", direction: "upstream"})
→ 6 test functions at d=1 → well tested
```

**Find what a test function exercises:**
```
nexus context({node_id: "py:function:tests.test_sn_1d.test_homogeneous_exact"})
→ Calls: solve_sn, get (derivations), GaussLegendre1D, homogeneous_1d
```

**Trace from a test to the equations it verifies:**
```
nexus trace_error({test_node_id: "py:function:tests.test_cp_slab.test_slab_cp_eigenvalue"})
→ Call chain → equations on path → citations
```

## Module Coverage Audit (complete workflow)

To assess which modules are tested and which aren't, run this sequence:

**Step 1 — List all solver modules:**
```
nexus query({text: "orpheus.", node_types: "module", limit: 50})
→ orpheus.sn.solver, orpheus.cp.solver, orpheus.fuel.solver, ...
```

**Step 2 — List all test modules:**
```
nexus query({text: "tests.test_", node_types: "module", limit: 50})
→ tests.test_sn_1d, tests.test_cp_slab, ... (missing: test_fuel, test_thermal, test_kinetics)
```

**Step 3 — For each solver, check upstream for test coverage:**
```
nexus impact({target: "py:function:orpheus.sn.solver.solve_sn", direction: "upstream"})
→ d=1 includes tests.test_sn_1d, tests.test_sn_solver_components → TESTED

nexus impact({target: "py:function:orpheus.fuel.solver.solve_fuel_behaviour", direction: "upstream"})
→ No tests.* nodes at any depth → UNTESTED
```

**Step 4 — For tested modules, count test functions:**
```
nexus neighbors({node_id: "py:module:tests.test_sn_1d", direction: "out", edge_types: "contains"})
→ 5 test functions

nexus neighbors({node_id: "py:module:tests.test_sn_solver_components", direction: "out", edge_types: "contains"})
→ 18 test functions (classes contain methods — count those too)
```

This replaces `ls orpheus/ && ls tests/ | diff` shell patterns. The
graph knows the relationships; shell only knows filenames.

## Tools

**verification_coverage** — the V&V map:
```
nexus verification_coverage({status_filter: "implemented"})
→ math:equation:alpha-cylindrical — has code (sweep_cylindrical) but no test
→ math:equation:surface-to-surface — has code (CPMesh._compute_slab_rcp) but no test
```

**staleness** — doc-code drift:
```
nexus staleness()
→ STALE: doc:theory/discrete_ordinates
→   sn_sweep.sweep_cylindrical modified 3 days ago, doc unchanged for 2 weeks
→   Affected symbols: sweep_cylindrical, compute_alpha_dome
```

**session_briefing** — combined overview:
```
nexus session_briefing()
→ 5 stale docs, 12 verification gaps, 452 unresolved references
→ Priority: update discrete_ordinates.rst
```
