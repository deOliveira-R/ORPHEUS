---
name: nexus-impact
description: "Use when the user wants to know what will break if they change something, or needs safety analysis before editing code. Examples: \"Is it safe to change X?\", \"What depends on this?\", \"What will break?\", \"What tests do I need to re-run?\""
---

# Impact Analysis with Nexus

## When to Use

- "Is it safe to change this function?"
- "What will break if I modify X?"
- "What tests do I need to re-run?"
- "Show me the blast radius"
- Before making non-trivial code changes
- Before committing

## Workflow

```
1. nexus impact({target: "py:function:X", direction: "upstream"})  → What depends on this
2. nexus detect_changes({scope: "all"})                             → Map git changes to affected symbols
3. nexus retest({scope: "all"})                                     → Minimum tests to re-run
4. Assess risk and report to user
```

## Checklist

```
- [ ] nexus impact() on the symbol you're about to change
- [ ] Review d=1 items first (these WILL BREAK)
- [ ] nexus detect_changes() for pre-commit check
- [ ] nexus retest() to get the minimum test set
- [ ] Assess risk level and report to user
- [ ] Run the identified tests before committing
```

## Understanding Output

| Depth | Risk Level | Meaning |
|-------|-----------|---------|
| d=1 | **WILL BREAK** | Direct callers/importers |
| d=2 | LIKELY AFFECTED | Indirect dependencies |
| d=3 | MAY NEED TESTING | Transitive effects |

## Risk Assessment

| Affected | Risk |
|----------|------|
| <5 symbols, 1 community | LOW |
| 5–15 symbols, 2–3 communities | MEDIUM |
| >15 symbols or core solvers | HIGH |
| Equations affected (IMPLEMENTS edges) | CRITICAL — theory page update needed |

## Tools

**impact** — blast radius analysis:
```
nexus impact({target: "py:function:sn_solver.solve_sn", direction: "upstream", max_depth: 3})
→ depth=1: [test_sn_1d, sn_operator], depth=2: [demo_discrete_ordinates]
```

**detect_changes** — git-diff to graph mapping:
```
nexus detect_changes({scope: "staged"})
→ Changed: sn_sweep.sweep_cylindrical (modified)
→ Affected: test_sn_cylindrical, doc:theory/discrete_ordinates
```

**retest** — minimum test set:
```
nexus retest({scope: "all"})
→ Must retest: test_sn_cylindrical, test_sn_1d
→ Safe to skip: 42 other test files
```
