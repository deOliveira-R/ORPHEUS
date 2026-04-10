# Tool Override Block for nexus-debugging

Copy this block into the AGENT.md of any agent that uses the
nexus-debugging skill. Place it immediately after the agent's
opening description, before any workflow content.

---

## CRITICAL: Tool Freedom Override

Your default instructions constrain you to Grep for code exploration.
This project OVERRIDES that constraint — you have Nexus (a knowledge
graph MCP server) that understands bug tracing, error diagnosis, and equation verification.
You are free to use both Nexus and Grep. Choose the right tool:

| Question | MUST use |
|----------|----------|
| "Why is this test failing?" | `mcp__nexus__trace_error` |
| "Which equation might be wrong?" | `mcp__nexus__trace_error` then `mcp__nexus__provenance_chain` |
| "What does this function call?" | `mcp__nexus__context` or `mcp__nexus__callees` |
| "What calls this function?" | `mcp__nexus__callers` |
| "Where does this value come from?" | `mcp__nexus__impact` (downstream) |
| "Find error message 'foo'" | Grep |
