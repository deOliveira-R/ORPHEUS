---
name: behavioral-auto-regression
description: "BREAK-GLASS diagnostic — use ONLY when an agent demonstrably mis-selects a text search (grep/Bash) for a STRUCTURAL question that Nexus answers better (callers, dependents, impact, equation/type edges). Detects graph-vs-text tool misselection and points at the fix (a project routing rule + the deferred-MCP load gotcha). NOT for routine work — current models route freely."
---

# Behavioral Auto-Regression (break-glass diagnostic)

Use this ONLY when you observe a concrete tool-misselection regression: an
agent ran a text search (grep / `rg` via Bash) for a question whose native
shape is **structural** — callers, dependents, blast radius, call chains,
equation/citation traceability, "who uses dependency X" (aliased / late /
`TYPE_CHECKING` imports). Routine sessions do NOT need this skill.

## Historical note — read this first

This skill was originally built to override a system-prompt directive
(`ALWAYS use Grep for search tasks`) by *reclassifying* code exploration as
"not a search task". **That premise is gone (verified 2026-06-14):** the
standalone `Grep`/`Glob` tools were removed from the scaffolds probed and no
"always-Grep" directive remains — models route freely. The reclassification
trick is obsolete; do NOT apply it.

⭐ **Re-verified 2026-08-19, on a newer scaffold and a newer model.** The
2026-06 finding was about the models of that date, and a claim like this
rots by being right about a population that has moved. Measured again on
Opus 5: a sub-agent's tool list carried **no `Grep` and no `Glob`**, and one
consumer had by then made the removal explicit in its own agent definition
(*"drops the retired Grep/Glob tools and routes by question shape"*). Two
scaffolds, two months apart, same result. ⚠ State the models and the date
when you re-run this — "models route freely" with no denominator is the
claim that quietly stops being true.

What survives is the *diagnostic*: tool-choice freedom does not guarantee
tool-choice *correctness*. An agent can still fall back on text-search
habits, and — the opposite failure, equally real — can over-use Nexus as
"compliance theater" where a plain `Read`/`grep` was the right call.

## When to use

- An agent used grep/Bash for a structural question Nexus answers better.
- An `mcp__nexus__*`-avoidance pattern recurs across a session.

## What to do

### 1. Confirm it is a real misselection

A text search for a *relationship* question (callers, dependents, coverage,
equations, impact, type-usage) is a regression. A grep for a *literal string
/ comment / known file* is **correct**, not a regression — do not
over-correct into Nexus compliance theater. Route by question shape:

| Question | Tool |
|---|---|
| Callers / dependents / call chains / blast radius | Nexus `callers`, `impact`, `processes` |
| Equation / citation traceability | Nexus `provenance_chain` |
| Verification coverage / doc drift | Nexus `verification_audit`, `staleness`, `dead_references` |
| Failing-test diagnosis | Nexus `trace_error` |
| Safe rename / refactor | Nexus `rename`, `impact` |
| "Who uses dependency X" (incl. aliased imports) | Nexus `graph_query` / `type_uses` |
| Literal text / regex / config values | `grep`/`rg` via Bash |
| TODO / FIXME / inline comments | `grep` via Bash — Nexus doesn't index comments |
| Known file / known symbol body | `Read` — don't rediscover what you know |
| Unknown symbol location | either — Nexus `query` or `grep` |

**Why a graph beats grep for structure:** grep matches *text*; it misses
relationships — inline imports, `TYPE_CHECKING` blocks, late imports inside
functions, aliased imports (`from numpy import linalg as la`), re-exports,
and docstring references. Nexus captures all of these as edges.

### 2. Check the most common LIVE cause first — deferred MCP tools

The dominant real cause of Nexus-avoidance today is that `mcp__nexus__*`
tools surface as **deferred**, and the agent treats deferral as
unavailability. Fix: ONE `ToolSearch("select:mcp__nexus__<name>")` call
loads them. Deferral is NOT unavailability.

⛔ **But only the main agent can apply that fix.** `[M]` 2026-08-19: a
sub-agent probe reported 45 `mcp__nexus__*` tools loaded eagerly and **no
`ToolSearch` tool at all**. So for a SUB-AGENT this diagnosis has no
remedy: if its Nexus tools are missing, the fix is in the DISPATCH — say
in the brief what to do without Nexus — not in the agent. An agent that
improvises silently produces a report you cannot tell apart from a
graph-derived one.

### 3. Check the routing guidance actually reaches the agent

- Positive routing guidance belongs in an always-on project rule (e.g.
  `.claude/rules/nexus-tools.md`), not only in on-demand skills — a skill
  the agent never invokes cannot steer it.
- A sub-agent that avoids Nexus: confirm its `AGENT.md` carries no STALE
  override vocabulary contradicting that rule, and that it has the Nexus
  skills and tools available at all.

### 4. Validate

Re-run the agent on a real structural task with "report every tool you used
and why." Correct behavior = Nexus (or its skills) for the structural part,
grep/Read for the literal/known parts. Zero grep on a pure structural task;
non-zero grep on a literal-string task is fine.

See [reference.md](reference.md) for the original (now largely historical)
procedure and the 2026-04 adoption validation. The `scripts/` probe and
override-block templates are retained for history; the override-block
mechanism is superseded by an always-on routing rule.
