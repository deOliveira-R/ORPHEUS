# Nexus & code-exploration tools

This project ships **Nexus** (`sphinxcontrib-nexus`) — one knowledge graph unifying
code structure (call graphs, imports, inheritance, type annotations) with documentation
structure (equations, cross-references, citations, theory pages). It runs as an MCP
server; the graph rebuilds on every docs build and the server auto-reloads when the
database changes.

**Nexus is the structural code-intelligence layer.** It answers relationship questions
that text search fundamentally cannot.

**Why a graph beats grep for structure:** grep matches *text*; it misses relationships —
inline imports, `TYPE_CHECKING` blocks, late imports inside functions, aliased imports
(`from numpy import linalg as la`), re-exports, and docstring references. Nexus captures
all of these as graph edges.

You have **freedom of tool choice** — route by what the question actually is, and
**invoke the Nexus *skills*, not raw MCP tools**: each skill encodes the complete
workflow and names its tools, and `nexus-guide` lists them all.

| Question | Route | Why |
|---|---|---|
| Callers, dependents, blast radius, call chains — "what breaks if I change X" | Nexus: `nexus-impact` | text cannot follow an edge |
| "Who uses X", including aliased, late and `TYPE_CHECKING` imports | Nexus: `nexus-exploring` | grep misses exactly those |
| A structural smell: clones (`twin_paths`), dead code (`dead_functions`), one tag branched on everywhere (`discriminations`), an implicit interface (`protocol_conformers`), a helper in the wrong module (`native_place`) | Nexus: `nexus-exploring` | whole-graph sweeps; the symptom names its tool |
| The map of a codebase: its areas (`communities`, `god_nodes`), the few nodes holding areas together (`bridges`); how a subsystem works; what actually RAN | Nexus: `nexus-exploring` | hubs, bridges and the runtime overlay |
| Docs ↔ code ↔ test traceability, verification coverage, doc drift, dead documentation references, which tests must re-run | Nexus: `nexus-verification` | with a coverage capture (`run=`) it answers from what EXECUTED |
| A wrong answer, a failing test | Nexus: `nexus-debugging` | walks the call graph to the suspect equation |
| A rename, an extraction, a retirement | Nexus: `nexus-refactoring` | references by graph, not by text |
| Literal text, a regex, a config value, a TODO/FIXME or any comment | `grep`/`rg` via **Bash** | raw strings; Nexus does not index comments |
| A file or a symbol body you already know | **Read** (or `find` via Bash) | don't rediscover what you already know |
| Where an unknown symbol lives | either — Nexus `query` or `grep` | your call |

Over-using Nexus where a plain `Read` or `grep` was correct is as much a misselection as
grepping for a relationship question. Do not perform compliance theater.

**Users describe symptoms, not tools.** "We keep changing these classes in lockstep",
"two people built this separately", "things live in surprising places", "the docs feel
out of date" are all graph questions — route them to `protocol_conformers` / `twin_paths`
/ `native_place` / `dead_references` respectively rather than reading files until a
pattern appears (`nexus-exploring` carries the full table, "What the user actually says").

**Some checks are part of the job, not a request.** After you delete or rename anything,
run `dead_references` before calling it done — green tests do not cover prose, and a dead
documentation reference produces no build warning at any severity, so nothing else will
catch it. Before a release, and for any "health check" or onboarding review, run the
sweeps `nexus-exploring` lists under "Sweeps you run WITHOUT being asked" (the smell
family, `dead_references` and `staleness`).

**Operational notes**

- **Deferred tools — ⛔ and the escape hatch is MAIN-AGENT-ONLY.** If `mcp__nexus__*`
  surface as deferred, ONE `ToolSearch("select:mcp__nexus__<name>")` loads them —
  deferral is NOT unavailability. ⚠ **A sub-agent has no `ToolSearch` tool**, so this recovery path
  does not exist for it — `[M]` 2026-08-19, a sub-agent probe reported 45 `mcp__nexus__*`
  tools loaded eagerly and **no `ToolSearch` at all**. A sub-agent that finds Nexus
  genuinely absent cannot recover; it must say so and fall back to `Bash` (grep, or
  `python -c "from sphinxcontrib.nexus.export import load_sqlite"` against the graph DB).
  ⟹ **when a dispatch depends on Nexus, say in the brief what to do if it is missing** —
  otherwise the agent improvises silently, and its report cannot be told apart from a
  grep-derived one.
  This is the most common cause of an agent silently avoiding the graph.
- **Stale graph:** rebuild the docs first; the MCP server auto-reloads.
- **Git worktrees:** the session's MCP server may have been launched against the MAIN
  checkout's graph, so every query answers from the wrong branch until you switch. Build
  inside the worktree, then `use_workspace(<worktree root>)`; `workspaces` lists every
  checkout and its graph.
  ⚠ **`session_briefing` warns when files the graph INDEXES have changed — not when the
  branch differs.** Those are different questions, and reading the second into the first
  makes the warning look broken: an ordinary ff-merge-and-delete leaves the graph
  describing the checkout exactly while the branch name has moved on (`[M]` 2026-08-16:
  25 files differed from the build commit, **0 of them indexed**, and the briefing was
  right to stay quiet). ⟹ *silence means the indexed sources match*, not "same branch".
