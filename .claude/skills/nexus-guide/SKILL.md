---
name: nexus-guide
description: "Use when the user asks about Nexus itself — available tools, how to query the knowledge graph, MCP resources, graph schema, or workflow reference. Examples: \"What Nexus tools are available?\", \"How do I use Nexus?\", \"What can I query?\""
---

# Nexus Guide

Quick reference for all Nexus MCP tools, resources, and the knowledge graph schema.

## Always Start Here

For any task involving code understanding, debugging, impact analysis, or verification:

1. **Read `nexus://briefing`** — session overview: stale docs, coverage gaps, recent changes
2. **Match your task to a skill below**
3. **Follow that skill's workflow and checklist**

## Skills

| Task | Skill |
|------|-------|
| Understand architecture / "How does X work?" | `nexus-exploring` |
| Blast radius / "What breaks if I change X?" | `nexus-impact` |
| Trace bugs / "Why is X failing?" / "Which equation is wrong?" | `nexus-debugging` |
| Rename / extract / split / refactor | `nexus-refactoring` |
| V&V status / "Which docs are stale?" | `nexus-verification` |
| Dependency migration / "Plan numpy → jax" | `nexus-migration` |
| CLI commands (analyze, serve) | `nexus-cli` |
| Tools, resources, schema reference | `nexus-guide` (this file) |

## Tools Reference (16 tools)

### Standard
| Tool | What it answers |
|------|----------------|
| `query` | "Find code/docs by keyword" |
| `context` | "360-degree view of a symbol" |
| `neighbors` | "What's directly connected?" |
| `impact` | "What breaks if I change this?" |
| `shortest_path` | "How do these concepts connect?" |
| `god_nodes` | "What are the central concepts?" |
| `stats` | "Graph summary" |

### Safety & Refactoring
| Tool | What it answers |
|------|----------------|
| `communities` | "What are the functional modules?" |
| `detect_changes` | "What did my git changes affect?" |
| `rename` | "Safe multi-file rename with preview" |
| `retest` | "Minimum tests to re-run after changes" |

### Code+Doc Fusion (unique to Nexus)
| Tool | What it answers |
|------|----------------|
| `provenance_chain` | "Citation → equation → code chain" |
| `verification_coverage` | "What's verified, what's not?" |
| `staleness` | "Which docs drifted from code?" |
| `session_briefing` | "What do I need to know right now?" |
| `trace_error` | "Which equation is wrong?" |
| `migration_plan` | "Plan my dependency migration" |

## Resources (4 resources)

| Resource | What you get |
|----------|-------------|
| `nexus://graph/stats` | Node/edge counts by type |
| `nexus://graph/communities` | Functional area summaries |
| `nexus://graph/schema` | Node types, edge types, ID format examples |
| `nexus://briefing` | Session briefing (stale docs, coverage gaps, changes) |

## Node ID Format

```
<domain>:<type>:<qualified_name>

py:function:sn_solver.solve_sn
py:class:collision_probability.CPMesh
py:method:CPMesh.compute_pinf_group
py:module:sn_solver
math:equation:alpha-recursion
doc:theory/discrete_ordinates
std:label:theory-collision-probability
std:term:widget
```

## Edge Types (13)

| Edge | Meaning | Source |
|------|---------|--------|
| `contains` | Parent → child (toctree, module→function, class→method) | Sphinx + AST |
| `references` | Cross-reference (`:ref:`, `:term:`) | Sphinx |
| `documents` | Doc page → code symbol (`:func:`, `:class:`) | Sphinx |
| `equation_ref` | Doc → equation (`:eq:`) | Sphinx |
| `cites` | Doc → citation | Sphinx |
| `implements` | Code → equation (inferred from co-occurrence) | Merge |
| `calls` | Function → function | AST |
| `imports` | Module → module | AST |
| `inherits` | Class → parent class | AST |
| `type_uses` | Function → type (from annotations) | AST |
| `tests` | Test function → tested function | AST |
| `derives` | Derivation → equation | AST |
