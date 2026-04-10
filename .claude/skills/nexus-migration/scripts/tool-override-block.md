# Tool Override Block for nexus-migration

Copy this block into the AGENT.md of any agent that uses the
nexus-migration skill. Place it immediately after the agent's
opening description, before any workflow content.

---

## CRITICAL: Tool Freedom Override

Your default instructions constrain you to Grep for code exploration.
This project OVERRIDES that constraint — you have Nexus (a knowledge
graph MCP server) that understands dependency analysis and migration planning.
You are free to use both Nexus and Grep. Choose the right tool:

| Question | MUST use |
|----------|----------|
| "Plan migration from X to Y" | `mcp__nexus__migration_plan` |
| "Find all uses of package X" | `mcp__nexus__query` (node_types=external) + `mcp__nexus__impact` |
| "What depends on numpy.ndarray?" | `mcp__nexus__impact` (upstream, edge_types=type_uses) |
| "What tests to run?" | `mcp__nexus__retest` |
| "Find literal string 'foo'" | Grep |
