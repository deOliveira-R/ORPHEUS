---
description: Documentation health — dead references and drift, computed now
allowed-tools: Bash
---

The following was computed from the knowledge graph at the moment you were
invoked. It is a finding, not a suggestion to go look for one.

## Dead documentation references

!`N="${CLAUDE_PROJECT_DIR:-.}/.venv/bin/nexus"; [ -x "$N" ] || N="$(command -v nexus 2>/dev/null)"; D="${NEXUS_DB:-${CLAUDE_PROJECT_DIR:-.}/docs/_build/html/_nexus/graph.db}"; if [ -n "$N" ] && [ -f "$D" ]; then "$N" dead-references --db "$D" --format text --limit 25; else echo "(nexus or graph not found — run 'nexus setup' and build the docs; set NEXUS_DB to override the graph path)"; fi`

## Timestamp drift

!`N="${CLAUDE_PROJECT_DIR:-.}/.venv/bin/nexus"; [ -x "$N" ] || N="$(command -v nexus 2>/dev/null)"; D="${NEXUS_DB:-${CLAUDE_PROJECT_DIR:-.}/docs/_build/html/_nexus/graph.db}"; if [ -n "$N" ] && [ -f "$D" ]; then "$N" staleness --db "$D" 2>/dev/null | head -40; fi`

---

Address the dead references above: each one is a doc, docstring, or quoted
type annotation citing code or an equation label that **no longer exists**.
Sphinx renders them as plain text and emits no warning at any severity, so
they will not surface any other way.

For each target, decide and act:

- the symbol was **renamed** → update the reference to the new name;
- the symbol was **deleted** and the prose is now wrong → rewrite or remove
  the passage, don't just unlink it;
- the reference is to something **never implemented** → say so explicitly in
  the prose, or drop the claim.

Report what you changed and what you deliberately left, with reasons. If a
reported target actually resolves at runtime through `__getattr__` or
metaclass machinery, say so — static analysis cannot see that, and the
finding is a false positive worth recording.
