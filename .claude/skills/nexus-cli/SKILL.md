---
name: nexus-cli
description: "Use when the user needs to run Nexus CLI commands like analyze source code, start the MCP server, or check graph status. Examples: \"Index this repo\", \"Analyze the codebase\", \"Start the nexus server\""
---

# Nexus CLI

## Commands

### Analyze Source Code

```bash
# Basic: analyze current directory
nexus analyze .

# With SQLite output
nexus analyze . --db _nexus/graph.db

# With specific sys.path entries (for non-standard project layouts)
nexus analyze . --sys-path 01.Discrete.Ordinates 02.Collision.Probability

# Auto-detect numbered directories as sys.path entries
nexus analyze . --auto-sys-path

# Also write JSON
nexus analyze . --db _nexus/graph.db --json _nexus/graph.json

# Merge with existing Sphinx-generated graph
nexus analyze src/ --db docs/_build/html/_nexus/graph.db
```

### Start MCP Server

```bash
# Start the MCP server
nexus serve --db _nexus/graph.db

# With project root for git operations
nexus serve --db _nexus/graph.db --project-root /path/to/project
```

### MCP Configuration

Add to Claude Code's MCP config:
```json
{
  "mcpServers": {
    "nexus": {
      "command": "nexus",
      "args": ["serve", "--db", "/path/to/graph.db", "--project-root", "/path/to/project"]
    }
  }
}
```

## Sphinx Integration

Add to `docs/conf.py`:
```python
extensions = ['sphinxcontrib.nexus']

# Optional configuration
nexus_output = '_nexus'           # Output directory (default)
nexus_ast_analyze = True          # Run AST analysis during build (default: True)
```

After `sphinx-build`, both `_nexus/graph.db` (SQLite) and `_nexus/graph.json` are generated.

## Graph Freshness

The graph is automatically rebuilt during every `sphinx-build`. For standalone use:

```bash
# Re-analyze after code changes
nexus analyze . --db _nexus/graph.db

# The MCP server loads from the database — restart it after re-analysis
```
