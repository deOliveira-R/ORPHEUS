# ORPHEUS — Open Reactor Physics Educational University System

## Session Start Protocol

Every session, before doing any work:

1. Read `.claude/lessons.md` — behavioral corrections from past sessions
2. Identify which module(s) the task involves
3. Dispatch the **explorer** agent on that module (reads Sphinx theory + Nexus graph)
4. Check open GitHub Issues for that module: `gh issue list -l module:<name>`

---

## Python Environment

All Python: `.venv/bin/python` (Python 3.14). Bare `python`/`python3`/`pip` are blocked by a PreToolUse hook.

---

## Cardinal Rules

These are non-negotiable. Violation of any cardinal rule is a session failure.

### 1. Track every improvement

NEVER let an improvement opportunity pass undocumented. When you notice
something that could be better — during implementation, debugging, code
review, or analysis — create a GitHub Issue immediately with the
appropriate `module:` label.

### 2. Sphinx IS the LLM's brain

Sphinx documentation is NOT a concise summary. It is the **specialized
knowledge base** that makes future sessions sharp. After implementing a
feature, the Sphinx docs must include: full derivations, design rationale
(not just what — WHY), what was tried and failed, gotchas, literature
references with equation numbers, numerical evidence. A feature is not
DONE until Sphinx is this thorough.

### 3. Log every caught bug

Every bug caught during development → `tests/l0_error_catalog.md` with:
ERR-NNN, failure mode (1–6), how it hid, which test catches it, lesson.
This is a QA publication artifact.

### 4. V&V levels are not interchangeable

```
VERIFICATION — "Are we solving the equations right?"
  L0  Term Verification        hand calc vs code, per term
  L1  Equation Verification    analytical solutions, MMS, convergence order
  L2  Integration Testing      multi-group + heterogeneous, self-convergence
VALIDATION — "Are we solving the right equations?"
  L3  Validation               comparison against experiment (ICSBEP, IRPhE)
INFORMATIONAL
  L4  Benchmarking             code-to-code — never proves correctness
```

Each level requires the levels below it. 1-group tests are DEGENERATE
(k = νΣ_f/Σ_a regardless of flux shape). Multi-group (≥2G) is mandatory.
Synthetic XS at L0–L2; real data only at L3+.

### 5. Absolute discipline every session

At session start: check open Issues for context.
During work: tag every new improvement immediately.
At session end: verify no orphan TODOs exist outside GitHub Issues.

---

## Knowledge Architecture

### Nexus Knowledge Graph

Nexus (`sphinxcontrib-nexus`) is the single knowledge graph for ORPHEUS.
It unifies **code structure** (call graphs, imports, inheritance, type
annotations) and **documentation structure** (equations, cross-references,
citations, theory pages) in one queryable graph. It runs as an MCP server
with 20 tools and 4 resources.

The graph is rebuilt automatically during every `sphinx-build`. Use the
MCP tools (`mcp__nexus__*`) for all code and documentation exploration.

**Before modifying a solver**: run `mcp__nexus__impact` for blast radius.
**Before committing**: run `mcp__nexus__detect_changes` to verify scope.
**When debugging**: run `mcp__nexus__trace_error` to trace from test to equations.
**When exploring**: run `mcp__nexus__context` for 360-degree symbol view.

For detailed workflows, read the nexus skills:

| Task | Skill |
|------|-------|
| "How does X work?" | `nexus-exploring` |
| "What breaks if I change X?" | `nexus-impact` |
| "Why is X failing?" / "Which equation is wrong?" | `nexus-debugging` |
| Rename / extract / refactor | `nexus-refactoring` |
| V&V status / "Which docs are stale?" | `nexus-verification` |
| Dependency migration | `nexus-migration` |
| Full tool/resource reference | `nexus-guide` |

### Sphinx Documentation (`docs/theory/`)

Theory, derivations, and equations behind each solver. Each theory page
has a **Key Facts** header — the essential equations, gotchas, and design
decisions. Read Key Facts before modifying any solver.

- **Before modifying a solver**: dispatch the **explorer** agent (reads theory + Nexus graph)
- **After modifying equations**: update the theory page and rebuild Sphinx
- **Documentation tasks**: use the **archivist** agent

### GitHub Issues

Single source of truth for improvements, bugs, and feature requests.
Labels: `module:sn`, `module:cp`, `module:moc`, `module:mc`, `module:diffusion`,
`module:geometry`, `module:data`, `module:tests`, `module:docs`,
`level:L0`, `level:L1`, `level:L2`, `type:bug`, `type:improvement`, `type:feature`.

---

## Specialized Agent Fleet

Six agents in `.claude/agents/` — use them instead of generic subagents.
Each has a `lessons.md` that sharpens across sessions.

| Agent | Invoke when | Key rule |
|-------|-------------|----------|
| **explorer** | Understanding code (mandatory first step) | Uses Nexus MCP + Sphinx. Replaces built-in Explore. |
| **archivist** | Writing/reviewing Sphinx docs | Sphinx-as-brain: full derivations. Uses Nexus for cross-ref audits and staleness detection. |
| **qa** | Reviewing code, validating claims | V&V hierarchy. 6 AI failure modes. Multi-group mandatory. |
| **numerics-investigator** | Solver gives wrong answers | 7-step diagnostic cascade. Uses Nexus trace_error for equation tracing. |
| **literature-researcher** | Need equations from papers | Source priority by topic. Maps notation to ORPHEUS. |
| **test-architect** | Planning verification BEFORE implementation | Analytical solution catalog. 1-group is degenerate. |

**After every agent invocation**: review the output with full session
context before committing. Sub-agents lack conversation history.

---

## Working Principles

- **Plan first**: enter plan mode for any non-trivial task (3+ steps)
- **Verify before done**: run tests, compare expected vs actual
- **Demand elegance**: is there a simpler way? No hacky fixes.
- **Fix root causes**: trace bugs to root cause, not symptoms
- **Code style**: Pythonic (dataclasses, type hints, scipy). No 1:1 MATLAB copy. Never transcribe values manually.

---

## Nexus — Quick Reference

### Always Do

- **MUST run impact analysis before editing any symbol.** `mcp__nexus__impact({target: "py:function:X", direction: "upstream"})` — report blast radius to the user.
- **MUST run detect_changes before committing.** `mcp__nexus__detect_changes({scope: "all"})` — verify only expected symbols affected.
- **MUST warn the user** if impact shows >15 affected symbols or equations on the path.
- When exploring, prefer `mcp__nexus__query` and `mcp__nexus__context` over grep.

### Impact Risk Levels

| Depth | Meaning | Action |
|-------|---------|--------|
| d=1 | WILL BREAK — direct callers/importers | MUST update these |
| d=2 | LIKELY AFFECTED — indirect deps | Should test |
| d=3 | MAY NEED TESTING — transitive | Test if critical path |

### Self-Check Before Finishing

1. `mcp__nexus__impact` was run for all modified symbols
2. No high-risk warnings were ignored
3. `mcp__nexus__detect_changes` confirms changes match expected scope
4. All d=1 dependents were updated

### Resources

| Resource | Use for |
|----------|---------|
| `nexus://graph/stats` | Graph overview |
| `nexus://graph/communities` | Functional areas |
| `nexus://graph/schema` | Node types, edge types, ID format |
| `nexus://briefing` | Session briefing: stale docs, coverage gaps, changes |
