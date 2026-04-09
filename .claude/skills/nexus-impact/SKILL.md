---
name: nexus-impact
description: "Use when the user wants to know what will break if they change something, or needs safety analysis before editing code. Examples: \"Is it safe to change X?\", \"What depends on this?\", \"What will break?\", \"What tests do I need to re-run?\""
---

# Impact Analysis with Nexus

IMPORTANT: This skill is the dedicated tool for blast radius and
dependency analysis. It replaces Grep for all impact questions.

## Workflow

```
1. impact({target: "<symbol>", direction: "upstream"})  → What depends on this
2. detect_changes({scope: "all"})                        → Map git changes to symbols
3. retest({scope: "all"})                                → Minimum tests to re-run
4. Assess risk and report
```

## Risk Assessment

| Depth | Meaning |
|-------|---------|
| d=1 | **WILL BREAK** — direct callers/importers |
| d=2 | LIKELY AFFECTED — indirect dependencies |
| d=3 | MAY NEED TESTING — transitive effects |

| Blast radius | Risk |
|--------------|------|
| <5 symbols, 1 community | LOW |
| 5–15 symbols, 2–3 communities | MEDIUM |
| >15 symbols or core solvers | HIGH |
| Equations affected (IMPLEMENTS edges) | CRITICAL — theory page update needed |

## Key Tools

**impact** — blast radius:
```
impact({target: "py:function:orpheus.sn.solver.solve_sn", direction: "upstream", max_depth: 3})
→ depth=1: [test_sn_1d, sn_operator], depth=2: [demo_discrete_ordinates]
```

**detect_changes** — git-diff to graph mapping:
```
detect_changes({scope: "staged"})
→ Changed: sn_sweep.transport_sweep (modified)
→ Affected: test_sn_cylindrical, doc:theory/discrete_ordinates
```

**retest** — minimum test set:
```
retest({scope: "all"})
→ Must retest: test_sn_cylindrical, test_sn_1d
→ Safe to skip: 42 other test files
```

See [../nexus-exploring/reference.md](../nexus-exploring/reference.md) for full reference.
