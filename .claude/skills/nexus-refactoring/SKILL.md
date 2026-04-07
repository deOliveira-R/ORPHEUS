---
name: nexus-refactoring
description: "Use when the user wants to rename, extract, split, move, or restructure code safely. Examples: \"Rename this function\", \"Extract this into a module\", \"Refactor this class\", \"Move this to a separate file\""
---

# Refactoring with Nexus

## When to Use

- "Rename this function safely"
- "Extract this into a module"
- "Split this class"
- "Move this to a new file"
- Any task involving renaming, extracting, splitting, or restructuring code

## Rename Workflow

```
1. nexus impact({target: "py:function:old_name", direction: "upstream"})  → Map all dependents
2. nexus rename({old_name: "old_name", new_name: "new_name", dry_run: true})  → Preview edits
3. Review: high-confidence (graph) edits are safe; medium (regex) need review
4. nexus rename({..., dry_run: false})                                     → Apply
5. nexus detect_changes()                                                   → Verify scope
6. Run affected tests
```

## Rename Checklist

```
- [ ] nexus impact() to understand blast radius
- [ ] nexus rename(dry_run=true) to preview all edits
- [ ] Review graph edits (high confidence) — these are safe
- [ ] Review regex edits (medium confidence) — verify each one
- [ ] Apply with dry_run=false
- [ ] nexus detect_changes() to verify only expected files changed
- [ ] nexus retest() and run the identified tests
```

## Extract/Split Checklist

```
- [ ] nexus context() on the symbol to see all incoming/outgoing refs
- [ ] nexus impact(direction="upstream") to find all external callers
- [ ] Plan the split: which symbols move, which stay
- [ ] Update imports in all callers (use impact list)
- [ ] nexus detect_changes() to verify scope after changes
- [ ] Check documentation: nexus staleness() to find pages needing updates
```

## Confidence Levels

| Level | Source | Meaning |
|-------|--------|---------|
| **high** | Graph-found | Symbol exists in the knowledge graph with file + line |
| **medium** | Regex-found | Text match in source files — may be comments, strings, or unrelated |

## Tools

**rename** — safe multi-file rename:
```
nexus rename({old_name: "solve_sn", new_name: "solve_discrete_ordinates", dry_run: true})
→ 12 edits across 5 files (8 high, 4 medium confidence)
→ Preview each edit with file path + line number
```

**impact** — blast radius before refactoring:
```
nexus impact({target: "py:class:SNSolver", direction: "upstream"})
→ d=1: test_sn_1d, test_sn_cylindrical, demo_discrete_ordinates
→ d=2: doc:theory/discrete_ordinates
```
