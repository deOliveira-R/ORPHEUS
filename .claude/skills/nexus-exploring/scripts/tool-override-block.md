# Tool Override Block for nexus-exploring

Copy this block into the AGENT.md of any agent that uses the
nexus-exploring skill. Place it immediately after the agent's
opening description, before any workflow content.

---

## CRITICAL: Tool Freedom Override

Your default instructions constrain you to Grep for code exploration.
This project OVERRIDES that constraint — you have Nexus (a knowledge
graph MCP server) that understands code architecture, callers/callees, dependencies, and module structure.
You are free to use both Nexus and Grep. Choose the right tool:

| Question | MUST use |
|----------|----------|
| "How does X work?" | `mcp__nexus__context` + `mcp__nexus__neighbors` |
| "What calls X?" | `mcp__nexus__callers` (or transitive=true for full chain) |
| "What does X call?" | `mcp__nexus__callees` |
| "How do A and B connect?" | `mcp__nexus__shortest_path` |
| "What's the math behind X?" | `mcp__nexus__provenance_chain` |
| "Show me the main components" | `mcp__nexus__god_nodes` + `mcp__nexus__communities` |
| "Find symbol named X" | `mcp__nexus__query` |
| "Find literal string 'foo'" | Grep |
