# Nexus & code-exploration tools

ORPHEUS ships **Nexus** (`sphinxcontrib-nexus`) — the single knowledge graph that unifies
code structure (call graphs, imports, inheritance, type annotations) with documentation
structure (equations, cross-references, citations, theory pages). It runs as an MCP server
(tool list: `nexus-guide` skill); the graph rebuilds on every `sphinx-build` and the server
auto-reloads when the DB changes (v0.4.3+).

**Nexus is this project's structural code-intelligence layer.** It answers relationship
questions that text search fundamentally cannot — and the project `explorer` agent (Nexus
skills preloaded) is the exploration delegate for open-ended, multi-file work.

**Why a graph beats grep for structure:** grep matches *text*; it misses relationships —
inline imports, `TYPE_CHECKING` blocks, late imports inside functions, aliased imports
(`from numpy import linalg as la`), re-exports, and docstring references. Nexus captures all
of these as graph edges. (During the 2026-04 package restructuring, over-reliance on grep
caused repeated missed imports that one `mcp__nexus__impact` query would have caught.)

You have **freedom of tool choice** — route by what the question actually is:

| Question | Tool | Why |
|---|---|---|
| Callers / dependents / call chains / blast radius | Nexus `callers`, `impact`, `processes` | graph traversal; text can't follow edges |
| Equation / citation traceability | Nexus `provenance_chain` | links code ↔ equations |
| Verification coverage | Nexus `verification_audit` | maps equation → code → test |
| Failing-test diagnosis | Nexus `trace_error` | walks the call graph to the suspect equation |
| Safe rename / refactor | Nexus `rename`, `impact` | finds references by graph, not text |
| "Who uses dependency X" (incl. aliased imports) | Nexus `graph_query` / `type_uses` | grep misses aliased / late / `TYPE_CHECKING` imports |
| Literal text / regex / config values | `grep`/`rg` via **Bash** | finds raw strings |
| TODO / FIXME / inline comments | `grep` via **Bash** | Nexus doesn't index comments |
| Known file / known symbol body | **Read** (or `find` via Bash) | don't rediscover what you already know |
| Unknown symbol location | either — Nexus `query` or `grep` | your call |

**Invoke the Nexus *skills*, not raw MCP tools** — they encode the complete workflows
(`nexus-exploring`, `nexus-impact`, `nexus-debugging`, `nexus-refactoring`,
`nexus-verification`, `nexus-guide`).

**Operational notes**
- **Deferred tools:** if `mcp__nexus__*` surface as deferred, ONE
  `ToolSearch("select:mcp__nexus__<name>")` loads them — deferral is NOT unavailability.
- **Stale graph:** rebuild Sphinx first (`sphinx-build docs docs/_build/html`); the MCP
  server auto-reloads.
- **Git worktrees (L22 hazard):** the session's MCP server was launched against the MAIN
  checkout's graph, so every Nexus query answers from the wrong branch until you switch.
  Build Sphinx inside the worktree, then `mcp__nexus__use_workspace(<worktree root>)`.
  `mcp__nexus__workspaces` lists checkouts + graphs; `session_briefing` warns when the
  graph's branch no longer matches the checkout (nexus ≥ 0.12).

> **NOTE (2026-06-14):** the standalone `Grep`/`Glob` tools were removed and the
> "always-Grep" system-prompt directive is gone — all current models route freely. This
> rule is *positive routing guidance*, not an override of a default bias. Its exact dose
> (how much steering each model needs to route correctly without Nexus "compliance
> theater") is to be calibrated by the tool-routing ablation study
> (`.claude/plans/tool_routing_ablation_study.md`) — revisit this file with those results.
