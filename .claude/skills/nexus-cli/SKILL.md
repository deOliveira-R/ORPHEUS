---
name: nexus-cli
description: "Use when the user needs to run Nexus CLI commands like analyze source code, start the MCP server, or check graph status. Examples: \"Index this repo\", \"Analyze the codebase\", \"Start the nexus server\""
---

# Nexus CLI

## Where the graph lives — ask, never hardcode

A project declares its graph location ONCE, in `.nexus/config.toml`
under `[graph]`. The CLI, the MCP server and the Sphinx extension all
read it from there, so none of the commands below needs to be told the
path. Ask for it instead of writing it out:

```bash
# The graph's path, bare, for shell capture
nexus config db

# Anchored to a specific checkout (a worktree, or running from elsewhere)
nexus config db --project-root /path/to/checkout

# Every resolved setting, as `key = value`
nexus config
```

Exit 0 = the value is known; exit 1 = the key is unset or unknown.

**Precedence for `--db`, in order:** an explicit `--db` flag > `[graph].db`
in the nearest `.nexus/config.toml` > the unconfigured fallback
`_nexus/graph.db` (relative to the working directory, and almost never
right for a real project — the artefacts live under the Sphinx output
directory, which only the build knows).

The config file is found by walking UP from `--project-root` on the
subcommands that take that flag, and from the **current working
directory** on those that do not (`status`, `query`, and most read-only
queries). So run `nexus` from inside the project tree, pass
`--project-root` where it is accepted, or pass `--db` explicitly.

⟹ **Pass `--db` only to override the declaration deliberately** — a
scratch graph, a second checkout, a graph you are diffing against.
Retyping the configured path by hand creates a second declaration, and
the second declaration is the one that goes stale.

## Commands

### Analyze Source Code

```bash
# Analyze the current directory into the configured graph.
# Merges with the graph already there — this is how AST detail gets
# added to a Sphinx-built graph.
nexus analyze .

# A different source tree, same configured graph
nexus analyze src/

# With specific sys.path entries (for non-standard project layouts)
nexus analyze . --sys-path 01.Discrete.Ordinates 02.Collision.Probability

# Auto-detect numbered directories as sys.path entries
nexus analyze . --auto-sys-path

# Also write a JSON copy to a path you name
nexus analyze . --json /tmp/graph.json

# Deliberate override — index into a scratch graph, leaving the
# project's own graph untouched
nexus analyze . --db /tmp/scratch-graph.db
```

### Start MCP Server

```bash
# Serves the configured graph
nexus serve

# From outside the checkout, or to pick one worktree of several
nexus serve --project-root /path/to/project
```

### MCP Configuration

Add to Claude Code's MCP config. Give the **project root**, not the
database — the server resolves the graph through that root's
`.nexus/config.toml`, so the entry keeps working when the graph moves:

```json
{
  "mcpServers": {
    "nexus": {
      "command": "nexus",
      "args": ["serve", "--project-root", "/path/to/project"]
    }
  }
}
```

This repository's own `.mcp.json` is the working example, with both
values expressed relative to `${CLAUDE_PROJECT_DIR}` so the entry is
checkout-independent.

## Sphinx Integration

Add to `docs/conf.py`:
```python
extensions = ['sphinxcontrib.nexus']
```

Settings belong in `.nexus/config.toml`, not in `conf.py` — one
declaration, read by the extension, the CLI and the server alike:

```toml
[graph]
output = "graph"                          # directory under the Sphinx outdir
db = "docs/_build/html/graph/graph.db"    # where consumers should look
```

`conf.py` still accepts the older `nexus_*` values (`nexus_output`,
`nexus_ast_analyze`, …) as a LOWER-precedence tier, but settings split
across two files are settings that drift. Prefer the TOML; this project
has retired its `nexus_*` conf.py values entirely.

After `sphinx-build`, the graph artefacts land in the configured output
directory: `graph.db` (SQLite), `graph.json`, and `graph.html` (the
interactive explorer). `nexus config db` prints the resolved database
path.

## Graph Freshness

The graph is automatically rebuilt during every `sphinx-build`. For
standalone use:

```bash
# Re-analyze after code changes
nexus analyze .

# The MCP server loads from the database — restart it after re-analysis
```
