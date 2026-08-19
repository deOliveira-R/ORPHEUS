---
name: nexus-guide
description: "Use when the user asks about Nexus itself — available tools, how to query the knowledge graph, MCP resources, graph schema, or workflow reference. Examples: \"What Nexus tools are available?\", \"How do I use Nexus?\", \"What can I query?\""
---

# Nexus Guide

## Skills

| Task | Skill |
|------|-------|
| "How does X work?" | `nexus-exploring` |
| Structural smells (dead code, clones, missing types) | `nexus-exploring` |
| "What actually ran?" (runtime overlay) | `nexus-exploring` |
| Review a diff for architectural decay | `nexus-elegance` |
| "What breaks if I change X?" | `nexus-impact` |
| "Why is X failing?" | `nexus-debugging` |
| Rename / extract / refactor | `nexus-refactoring` |
| V&V status / stale docs / dead doc references | `nexus-verification` |
| Dependency migration | `nexus-migration` |
| CLI commands (analyze, serve) | `nexus-cli` |

## Tool choice: route by question shape

You have freedom of tool choice. Nexus is not a replacement for text
search — they answer different question shapes:

| Question | Tool |
|---|---|
| Relationships: callers, dependents, blast radius, call chains | Nexus |
| Doc↔code↔test traceability, coverage, drift | Nexus |
| "Who uses X" incl. aliased / late / `TYPE_CHECKING` imports | Nexus (grep misses these) |
| Literal text, regex, config values, TODO/FIXME comments | `grep`/`rg` via Bash |
| A file or symbol body you already know | `Read` |

Over-using Nexus where a plain `Read` was correct is as much a
misselection as grepping for a relationship question.

**Deferred tools:** if `mcp__nexus__*` surface as deferred, ONE
`ToolSearch("select:mcp__nexus__<name>")` loads them. Deferral is not
unavailability — treating it as such is the most common cause of an
agent silently avoiding the graph.

## Quick Start

```
1. READ nexus://briefing                    → Session overview
2. Match task to skill above
3. Follow that skill's workflow
```

## Bridges into the graph

- **From a position** (LSP result, stack trace, editor line):
  `node_at({file, line})` → graph node → `context`/`impact`.
- **From a node result**: AST-derived results carry
  `file_path`/`lineno` — open the source directly.
- **From an edit**: projects may inject `nexus file-brief` output
  via an edit-time hook — the brief's node IDs are entry points.
- **From a worktree**: after EnterWorktree, `use_workspace(<name>)`
  so queries answer from that checkout's graph.

## Full Reference

See [../nexus-exploring/reference.md](../nexus-exploring/reference.md) for
complete tool, resource, edge type, node ID format, graph_query syntax,
and CLI command reference.
