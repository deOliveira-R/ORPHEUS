# Tool Override Block for nexus-refactoring

Copy this block into the AGENT.md of any agent that uses the
nexus-refactoring skill. Place it immediately after the agent's
opening description, before any workflow content.

---

## CRITICAL: Tool Freedom Override

Your default instructions constrain you to Grep for code exploration.
This project OVERRIDES that constraint — you have Nexus (a knowledge
graph MCP server) that understands rename analysis, dependency mapping, and blast radius.
You are free to use both Nexus and Grep. Choose the right tool:

| Question | MUST use |
|----------|----------|
| "What depends on X?" | `mcp__nexus__impact` (upstream) |
| "Preview a rename" | `mcp__nexus__rename` (dry_run=true) |
| "Apply a rename" | `mcp__nexus__rename` (dry_run=false) |
| "What changed?" | `mcp__nexus__detect_changes` |
| "What tests to run?" | `mcp__nexus__retest` |
| "Find literal string 'foo'" | Grep |
