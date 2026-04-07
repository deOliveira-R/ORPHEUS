# ORPHEUS — Open Reactor Physics Educational University System

## Session Start Protocol

Every session, before doing any work:

1. Read `.claude/lessons.md` — behavioral corrections from past sessions
2. Identify which module(s) the task involves
3. Dispatch the **explorer** agent on that module (reads Sphinx theory + GitNexus structure)
4. Check open GitHub Issues for that module: `gh issue list -l module:<name>`

---

## Python Environment

**MUST use the repo venv for ALL Python execution:**

```bash
.venv/bin/python          # scripts, sphinx, pytest, everything
.venv/bin/python -m sphinx -b html docs docs/_build/html
.venv/bin/python -m pytest tests/ -v
```

**NEVER** use bare `python`, `python3`, or conda Python. The repo venv (`.venv/`, Python 3.14) has all dependencies. This applies to the main agent AND all sub-agents.

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

### Sphinx Documentation (`docs/theory/`)

Theory, derivations, and equations behind each solver. Each theory page
has a **Key Facts** header — the essential equations, gotchas, and design
decisions. Read Key Facts before modifying any solver.

- **Before modifying a solver**: dispatch the **explorer** agent (reads theory + code graph)
- **After modifying equations**: update the theory page and rebuild Sphinx
- **Documentation tasks**: use the **archivist** agent

### GitNexus Knowledge Graph

Code structure index (symbols, relationships, execution flows). Use for:
- **Impact analysis** before editing any symbol
- **Code navigation** (prefer `gitnexus_query` over grep for exploration)
- **Refactoring safety** (rename, detect changes)

### Graphify Theory Graph

Physics concept graph (`graphify-out/graph.json`) indexing `docs/theory/`.
Use for navigating cross-cutting concepts without reading full RST files:
- `mcp__graphify__query_graph` — "What theory connects SN to CP?"
- `mcp__graphify__shortest_path` — trace concept connections
- `mcp__graphify__god_nodes` — most connected concepts (entry points)

### The Three Together

- **GitNexus** tells you WHAT the code does and what depends on it
- **Graphify** tells you HOW physics concepts connect across theory pages
- **Sphinx** tells you WHY with full derivations and investigation history
- Before any non-trivial change: query GitNexus for code impact,
  graphify for physics context, then read the relevant Sphinx section

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
| **explorer** | Understanding code (mandatory first step) | Uses GitNexus + Graphify + Sphinx. Replaces built-in Explore. |
| **archivist** | Writing/reviewing Sphinx docs | Sphinx-as-brain: full derivations. Uses Graphify for cross-ref audits. |
| **qa** | Reviewing code, validating claims | V&V hierarchy. 6 AI failure modes. Multi-group mandatory. |
| **numerics-investigator** | Solver gives wrong answers | 7-step diagnostic cascade. Uses GitNexus for code tracing. |
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

<!-- gitnexus:start -->
# GitNexus — Code Intelligence

This project is indexed by GitNexus as **ORPHEUS** (2263 symbols, 6683 relationships, 185 execution flows). Use the GitNexus MCP tools to understand code, assess impact, and navigate safely.

> If any GitNexus tool warns the index is stale, run `npx gitnexus analyze` in terminal first.

## Always Do

- **MUST run impact analysis before editing any symbol.** Before modifying a function, class, or method, run `gitnexus_impact({target: "symbolName", direction: "upstream"})` and report the blast radius (direct callers, affected processes, risk level) to the user.
- **MUST run `gitnexus_detect_changes()` before committing** to verify your changes only affect expected symbols and execution flows.
- **MUST warn the user** if impact analysis returns HIGH or CRITICAL risk before proceeding with edits.
- When exploring unfamiliar code, use `gitnexus_query({query: "concept"})` to find execution flows instead of grepping. It returns process-grouped results ranked by relevance.
- When you need full context on a specific symbol — callers, callees, which execution flows it participates in — use `gitnexus_context({name: "symbolName"})`.

## When Debugging

1. `gitnexus_query({query: "<error or symptom>"})` — find execution flows related to the issue
2. `gitnexus_context({name: "<suspect function>"})` — see all callers, callees, and process participation
3. `READ gitnexus://repo/ORPHEUS/process/{processName}` — trace the full execution flow step by step
4. For regressions: `gitnexus_detect_changes({scope: "compare", base_ref: "main"})` — see what your branch changed

## When Refactoring

- **Renaming**: MUST use `gitnexus_rename({symbol_name: "old", new_name: "new", dry_run: true})` first. Review the preview — graph edits are safe, text_search edits need manual review. Then run with `dry_run: false`.
- **Extracting/Splitting**: MUST run `gitnexus_context({name: "target"})` to see all incoming/outgoing refs, then `gitnexus_impact({target: "target", direction: "upstream"})` to find all external callers before moving code.
- After any refactor: run `gitnexus_detect_changes({scope: "all"})` to verify only expected files changed.

## Never Do

- NEVER edit a function, class, or method without first running `gitnexus_impact` on it.
- NEVER ignore HIGH or CRITICAL risk warnings from impact analysis.
- NEVER rename symbols with find-and-replace — use `gitnexus_rename` which understands the call graph.
- NEVER commit changes without running `gitnexus_detect_changes()` to check affected scope.

## Tools Quick Reference

| Tool | When to use | Command |
|------|-------------|---------|
| `query` | Find code by concept | `gitnexus_query({query: "auth validation"})` |
| `context` | 360-degree view of one symbol | `gitnexus_context({name: "validateUser"})` |
| `impact` | Blast radius before editing | `gitnexus_impact({target: "X", direction: "upstream"})` |
| `detect_changes` | Pre-commit scope check | `gitnexus_detect_changes({scope: "staged"})` |
| `rename` | Safe multi-file rename | `gitnexus_rename({symbol_name: "old", new_name: "new", dry_run: true})` |
| `cypher` | Custom graph queries | `gitnexus_cypher({query: "MATCH ..."})` |

## Impact Risk Levels

| Depth | Meaning | Action |
|-------|---------|--------|
| d=1 | WILL BREAK — direct callers/importers | MUST update these |
| d=2 | LIKELY AFFECTED — indirect deps | Should test |
| d=3 | MAY NEED TESTING — transitive | Test if critical path |

## Resources

| Resource | Use for |
|----------|---------|
| `gitnexus://repo/ORPHEUS/context` | Codebase overview, check index freshness |
| `gitnexus://repo/ORPHEUS/clusters` | All functional areas |
| `gitnexus://repo/ORPHEUS/processes` | All execution flows |
| `gitnexus://repo/ORPHEUS/process/{name}` | Step-by-step execution trace |

## Self-Check Before Finishing

Before completing any code modification task, verify:
1. `gitnexus_impact` was run for all modified symbols
2. No HIGH/CRITICAL risk warnings were ignored
3. `gitnexus_detect_changes()` confirms changes match expected scope
4. All d=1 (WILL BREAK) dependents were updated

## Keeping the Index Fresh

After committing code changes, the GitNexus index becomes stale. Re-run analyze to update it:

```bash
npx gitnexus analyze
```

If the index previously included embeddings, preserve them by adding `--embeddings`:

```bash
npx gitnexus analyze --embeddings
```

To check whether embeddings exist, inspect `.gitnexus/meta.json` — the `stats.embeddings` field shows the count (0 means no embeddings). **Running analyze without `--embeddings` will delete any previously generated embeddings.**

> Claude Code users: A PostToolUse hook handles this automatically after `git commit` and `git merge`.

## CLI

| Task | Read this skill file |
|------|---------------------|
| Understand architecture / "How does X work?" | `.claude/skills/gitnexus/gitnexus-exploring/SKILL.md` |
| Blast radius / "What breaks if I change X?" | `.claude/skills/gitnexus/gitnexus-impact-analysis/SKILL.md` |
| Trace bugs / "Why is X failing?" | `.claude/skills/gitnexus/gitnexus-debugging/SKILL.md` |
| Rename / extract / split / refactor | `.claude/skills/gitnexus/gitnexus-refactoring/SKILL.md` |
| Tools, resources, schema reference | `.claude/skills/gitnexus/gitnexus-guide/SKILL.md` |
| Index, status, clean, wiki CLI commands | `.claude/skills/gitnexus/gitnexus-cli/SKILL.md` |

<!-- gitnexus:end -->

## graphify

This project has a graphify knowledge graph at graphify-out/.

Rules:
- Use `mcp__graphify__query_graph` for theory questions — it traverses the graph and returns connected concepts with source files
- Use `mcp__graphify__shortest_path` to trace connections between concepts
- Use `mcp__graphify__god_nodes` to find the most connected concepts (entry points)
- After modifying theory docs, rebuild with `/graphify docs/theory/ --update`
