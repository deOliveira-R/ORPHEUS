# Nexus Reference

Full tool, resource, and schema reference for the Nexus knowledge graph.
This file is shared across all nexus-* skills.

## Tools (21)

### Exploration
| Tool | What it answers | Key args |
|------|----------------|----------|
| `query` | Find symbols by keyword | `text`, `node_types`, `limit` |
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
| `retest` | Minimum test set after changes | `scope` |
| `rename` | Safe multi-file rename | `old_name`, `new_name`, `dry_run` |

### Code+Doc Fusion
| Tool | What it answers | Key args |
|------|----------------|----------|
| `provenance_chain` | Citation → equation → code chain | `node_id` |
| `verification_coverage` | V&V status map | `status_filter` |
| `verification_audit` | Complete V&V audit (single call) | — |
| `staleness` | Doc-code drift | — |
| `session_briefing` | Session overview | — |
| `trace_error` | Failing test → equations on path | `test_node_id` |
| `migration_plan` | Dependency migration phases | `from_dep`, `to_dep` |

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
math:equation:alpha-recursion
doc:theory/discrete_ordinates
std:label:theory-collision-probability
```

## Edge Types (12)

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

All output JSON to stdout. Default db: `_nexus/graph.db`.

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
nexus retest --db <path> [--project-root .] [--scope all|staged|unstaged|branch]
nexus changes --db <path> [--project-root .] [--scope all|staged|unstaged|branch]
nexus rename <old> <new> --db <path> [--project-root .] [--apply]
```
