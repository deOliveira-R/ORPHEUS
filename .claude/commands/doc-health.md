---
description: Documentation health — dead references and drift, computed now
allowed-tools: Bash
---

The following was computed from the knowledge graph at the moment you were
invoked. It is a finding, not a suggestion to go look for one.

## Dead documentation references

!`cd "${CLAUDE_PROJECT_DIR:-.}" || exit 0; N=".venv/bin/nexus"; [ -x "$N" ] || N="$(command -v nexus 2>/dev/null)"; if [ -z "$N" ]; then echo "(nexus not installed — run: bash scripts/setup.sh --no-docs)"; else D="${NEXUS_DB:-$("$N" config db 2>/dev/null)}"; if [ -f "$D" ]; then "$N" dead-references --db "$D" --format text --limit 25; else echo "(no graph at ${D:-<unresolved — check .nexus/config.toml>} — build it: sphinx-build docs docs/_build/html)"; fi; fi`

## Timestamp drift

!`cd "${CLAUDE_PROJECT_DIR:-.}" || exit 0; N=".venv/bin/nexus"; [ -x "$N" ] || N="$(command -v nexus 2>/dev/null)"; [ -n "$N" ] || exit 0; D="${NEXUS_DB:-$("$N" config db 2>/dev/null)}"; [ -f "$D" ] || exit 0; "$N" staleness --db "$D" --project-root . 2>/dev/null | head -40`

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
