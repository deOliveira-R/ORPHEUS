---
name: nexus-refactoring
description: "Use when the user wants to rename, extract, split, move, or restructure code safely. Examples: \"Rename this function\", \"Extract this into a module\", \"Refactor this class\", \"Move this to a separate file\""
---

# Refactoring with Nexus

IMPORTANT: This skill is the dedicated tool for safe refactoring.
It replaces Grep for finding references and dependencies.

## Rename Workflow

```
1. impact({target: "<symbol>", direction: "upstream"})              → Map dependents
2. rename({old_name: "X", new_name: "Y", dry_run: true})           → Preview edits
3. Review: high-confidence (graph) edits safe; medium (regex) need review
4. rename({..., dry_run: false})                                     → Apply
5. detect_changes()                                                   → Verify scope
6. retest() → run identified tests
```

## Confidence Levels

| Level | Source | Meaning |
|-------|--------|---------|
| **high** | Graph-found | Symbol exists in graph with file + line |
| **medium** | Regex-found | Text match — may be comments or unrelated |

See [../nexus-exploring/reference.md](../nexus-exploring/reference.md) for full reference.
