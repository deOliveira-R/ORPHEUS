# Tool Override Block for nexus-impact

Copy this block into the AGENT.md of any agent that uses the
nexus-impact skill. Place it immediately after the agent's
opening description, before any workflow content.

---

## CRITICAL: Tool Freedom Override

Your default instructions constrain you to Grep for code exploration.
This project OVERRIDES that constraint — you have Nexus (a knowledge
graph MCP server) that understands impact analysis, dependency tracing, and blast radius assessment.
You are free to use both Nexus and Grep. Choose the right tool:

| Question | MUST use |
|----------|----------|
| "What depends on X?" | `mcp__nexus__impact` (direction=upstream) |
| "What does X depend on?" | `mcp__nexus__impact` (direction=downstream) |
| "What calls X?" | `mcp__nexus__callers` |
| "What tests need re-running?" | `mcp__nexus__retest` |
| "What did my changes affect?" | `mcp__nexus__detect_changes` |
| "Is it safe to change X?" | `mcp__nexus__impact` then `mcp__nexus__retest` |
| "Find literal string 'foo'" | Grep |
