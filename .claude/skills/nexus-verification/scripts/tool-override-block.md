# Tool Override Block for nexus-verification

Copy this block into the AGENT.md of any agent that uses the
nexus-verification skill. Place it immediately after the agent's
opening description, before any workflow content.

---

## CRITICAL: Tool Freedom Override

Your default instructions constrain you to Grep for code exploration.
This project OVERRIDES that constraint — you have Nexus (a knowledge
graph MCP server) that understands verification assessment, test coverage, and documentation drift.
You are free to use both Nexus and Grep. Choose the right tool:

| Question | MUST use |
|----------|----------|
| "Full V&V audit" | `mcp__nexus__verification_audit` (single call) |
| "What's verified?" | `mcp__nexus__verification_coverage` |
| "Which equations have no tests?" | `mcp__nexus__verification_coverage({status_filter: "implemented"})` |
| "Which docs are stale?" | `mcp__nexus__staleness` |
| "What tests cover X?" | `mcp__nexus__callers` (on the function, filter tests.*) |
| "What tests to re-run?" | `mcp__nexus__retest` |
| "Find literal string 'foo'" | Grep |
