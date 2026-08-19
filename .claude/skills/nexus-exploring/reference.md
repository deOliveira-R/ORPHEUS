# Nexus Reference

Full tool, resource, and schema reference for the Nexus knowledge graph.
This file is shared across all nexus-* skills.

## Tools (45)

### Exploration
| Tool | What it answers | Key args |
|------|----------------|----------|
| `query` | Find symbols by keyword | `text`, `node_types`, `limit` |
| `file_brief` | What the graph knows about one FILE — module node, hub, equations it implements, doc pages owed an update, and for a test file what its gates verify plus the pytest command. Start here when all you have is a path | `file` |
| `node_at` | Map a file position (LSP result, stack trace) to the innermost enclosing node; warns when the file changed since the graph was built | `file`, `line` |
| `context` | 360-degree view of a symbol | `node_id` |
| `neighbors` | Direct connections | `node_id`, `direction`, `edge_types` |
| `callers` | Functions that call this symbol | `node_id`, `transitive`, `max_depth` |
| `callees` | Functions this symbol calls | `node_id`, `transitive`, `max_depth` |
| `shortest_path` | How two concepts connect | `source`, `target`, `max_hops` |
| `god_nodes` | Most connected symbols | `top_n` |
| `stats` | Graph summary | — |
| `communities` | Functional groupings | `min_size` |
| `bridges` | Nodes connecting communities | `top_n` |
| `processes` | Execution flows from entry points | `min_length` |
| `graph_query` | Structured traversal (Cypher-like) | `pattern`, `limit` |

### Safety & Refactoring
| Tool | What it answers | Key args |
|------|----------------|----------|
| `impact` | Blast radius analysis | `target`, `direction`, `max_depth`, `edge_types` |
| `detect_changes` | Git diff → graph mapping | `scope` |
| `retest` | Minimum test set after changes — pass `run` and a covered symbol is answered from EXECUTION, not from the 12-15 %-recall call cone | `scope`, `run`, `limit` |
| `doc_impact` | the static-cone dual of `retest` — documented claims a change to this symbol puts in question, with `page:line#anchor` and a verified flag | `node_id`, `limit` |
| `rename` | Safe multi-file rename | `old_name`, `new_name`, `dry_run` |

### Architecture Smells (missing-abstraction family)
| Tool | What it answers | Key args |
|------|----------------|----------|
| `native_place` | Functions that may belong inside a class (Feature-Envy) | `min_callers`, `exclude`, `limit` |
| `twin_paths` | Independent implementations of the same computation (Type-2/3 clones) | `min_similarity`, `min_tokens`, `exclude`, `limit` |
| `discriminations` | Tags discriminated at many sites (candidate missing types) | `min_sites`, `exclude`, `limit` |
| `dead_functions` | Functions/methods with no static callers (dead-code candidates) | `exclude`, `limit` |
| `protocol_conformers` | Classes satisfying a Protocol's method-set without declaring it | `min_methods`, `exclude`, `limit` |

### Runtime Overlay (dynamic execution-flow — the smell family's dynamic counterpart)
The static graph is *what can run*; a runtime overlay is *what actually ran*. Capture is consumer-side (run a canonical workload under a tracer), then `runtime_ingest` joins the artifact onto node-IDs and stores it in a sidecar (`<project root>/.nexus/traces/<run>.json`) — never in `graph.db`, which is rebuilt on every `sphinx-build`. The sidecar sits beside the database, *outside* the Sphinx build output, because it is durable state: a profiled run costs minutes to reproduce, and a directory a clean build deletes would destroy it. The query tools take `run` as one name OR comma-separated names to **union the canonical suite**.
| Tool | What it answers | Key args |
|------|----------------|----------|
| `runtime_ingest` | Overlay a `cProfile` / `coverage --branch` / `viztracer` trace on the graph | `artifact`, `kind`, `run`, `source_prefix` (list), `root` |
| `runtime_runs` | List ingested runs | — |
| `runtime_hotspots` | Hot path / iteration counts (the dynamic stage DAG) | `run`, `by` (cumtime/ncalls/tottime), `limit` |
| `runtime_edges` | Fired-vs-static edges: `dynamic_only` (dispatch the static graph missed), `fired`, `dead` | `run`, `mode`, `node`, `substantive_only`, `limit` |
| `runtime_markers` | Tests by marker, **as pytest resolved it** (module-level `pytestmark`, class marks, conftest hooks — invisible to a decorator walk); carries runnable pytest ids | `run`, `marker`, `node`, `limit` |
| `runtime_branches` | Partial-branch nodes; discriminators ranked first (missing-type suspects) | `run`, `node`, `partial_only`, `limit` |
| `runtime_exercisers` | Which tests EXECUTED a node — the falsifier for a coverage claim (needs contexts) | `run`, `node`, `limit` |
| `runtime_timeline` | Observed execution sequence (a viztracer run): nodes by first entry | `run`, `max_depth`, `limit` |

### Code+Doc Fusion
| Tool | What it answers | Key args |
|------|----------------|----------|
| `provenance_chain` | Citation → equation → code chain, plus the math-to-math spine in `relations` (what a discrete form discretizes / derives from / approximates) | `node_id` |
| `verification_coverage` | V&V status map | `status_filter` |
| `verification_audit` | Complete V&V audit (single call) | — |
| `verification_gaps` | Untagged tests, unverified equations, missing err catchers | `module`, `level` |
| `errors` | Catalogued failure modes and the tests that catch them, UNCAUGHT FIRST; `total_entries: 0` means nothing is declared, not nothing is wrong | `limit` |
| `staleness` | Doc-code drift (git timestamps) + dead-reference summary | — |
| `dead_references` | Docs/docstrings citing symbols or equation labels that NO LONGER EXIST (Sphinx renders these as plain text with no warning) | `limit` |
| `session_briefing` | Session overview | — |
| `trace_error` | Failing test → equations on path | `test_node_id` |
| `migration_plan` | Dependency migration phases | `from_dep`, `to_dep` |

### Workspaces (git worktrees)
A graph is a snapshot of ONE checkout. A session working in a worktree must query THAT worktree's graph.

| Tool | What it answers | Key args |
|------|----------------|----------|
| `workspaces` | Every checkout of the project, which graph each carries, and how fresh | — |
| `use_workspace` | Switch the active graph to another checkout | `ref` (name or root path) |

### Ingestion
| Tool | What it answers | Key args |
|------|----------------|----------|
| `ingest` | Add a document to the graph | `file_path`, `llm_command` |

## Resources (4)

| Resource | Content |
|----------|---------|
| `nexus://graph/stats` | Node/edge counts by type, density |
| `nexus://graph/communities` | Functional areas with top members |
| `nexus://graph/schema` | Node types, edge types, ID format |
| `nexus://briefing` | Session briefing: stale docs, gaps, changes |

## Node ID Format

```
<domain>:<type>:<qualified_name>

py:function:orpheus.sn.solver.solve_sn
py:class:orpheus.cp.solver.CPMesh
py:method:orpheus.cp.solver.CPMesh.compute_pinf_group
py:module:orpheus.sn.solver
py:tag:geometry
math:equation:alpha-recursion
prf:algorithm:transport-sweep
doc:theory/discrete_ordinates
std:label:theory-collision-probability
```

## Edge Types (16)

| Edge | Meaning | Source |
|------|---------|--------|
| `contains` | Parent → child (module→function, class→method) | Sphinx + AST |
| `references` | Cross-reference (`:ref:`, `:term:`) | Sphinx |
| `documents` | Doc page → code symbol (`:func:`, `:class:`) | Sphinx |
| `equation_ref` | Doc → equation (`:eq:`) | Sphinx |
| `cites` | Doc → citation | Sphinx |
| `implements` | Code → equation (inferred) | Merge |
| `calls` | Function → function | AST |
| `imports` | Module → module | AST |
| `inherits` | Class → parent class | AST |
| `type_uses` | Function → type (from annotations) | AST |
| `tests` | Test function → tested function | AST |
| `derives` | Derivation → equation | AST |
| `discriminates_on` | Function → tag it branches on (`if x == "..."`, `match`) | AST |
| `discretizes` | Discrete statement → the continuous one it discretizes | Directive |
| `derives_from` | Specialization → the parent it was reduced from | Directive |
| `approximates` | Closure/truncation → the exact form it stands in for | Directive |

## graph_query Pattern Syntax

```
source_type -edge_type-> target_type [WHERE field=value]

Examples:
  function -calls-> function           # all function-to-function calls
  file -contains-> equation            # all equations in doc pages
  * -implements-> equation             # code implementing equations
  function -type_uses-> external WHERE name=numpy*   # numpy usage
  * -cites-> *                         # all citation edges
```

Wildcards: `*` matches any type. `name=prefix*` for prefix match.

## CLI Commands (for ! injection)

All output JSON to stdout. `--db` is optional: each command derives the
project's graph at `<project root>/.nexus/graph.db`, which `nexus config db`
prints. Pass `--db` only to open a different graph.

```bash
nexus callers <node_id> --db <path> [--transitive] [--max-depth 3]
nexus callees <node_id> --db <path> [--transitive] [--max-depth 3]
nexus audit --db <path> [--project-root .]
nexus briefing --db <path>
nexus context <node_id> --db <path>
nexus neighbors <node_id> --db <path> [--direction in|out|both] [--edge-types calls,imports]
nexus god-nodes --db <path> [--top-n 10]
nexus communities --db <path> [--min-size 3]
nexus bridges --db <path> [--top-n 10]
nexus processes --db <path> [--min-length 3]
nexus shortest-path <source> <target> --db <path> [--max-hops 8]
nexus graph-query "<pattern>" --db <path> [--limit 50]
nexus trace <test_node_id> --db <path>
nexus runtime-ingest <artifact> [--kind cprofile|coverage|viztracer] [--run NAME] [--source-prefix PFX ...] [--root DIR] [--note TEXT]
nexus runtime-runs --db <path>
nexus runtime-hotspots --db <path> [--run NAME[,NAME...]] [--by cumtime|ncalls|tottime] [--limit 20]
nexus runtime-edges --db <path> [--run NAME[,NAME...]] [--mode dynamic_only|fired|dead] [--node SUBSTR] [--substantive-only] [--limit 50]
nexus runtime-branches --db <path> [--run NAME[,NAME...]] [--node SUBSTR] [--all] [--limit 50]
nexus runtime-exercisers --db <path> [--run NAME[,NAME...]] [--node SUBSTR] [--limit 50]
nexus runtime-timeline --db <path> [--run NAME] [--max-depth N] [--limit 50]
nexus retest --db <path> [--project-root .] [--scope all|staged|unstaged|branch]
               [--run <cov-run>[,<cov-run>]] [--limit N]
nexus changes --db <path> [--project-root .] [--scope all|staged|unstaged|branch]
nexus rename <old> <new> --db <path> [--project-root .] [--apply]
```
