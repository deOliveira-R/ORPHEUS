---
name: behavioral-auto-regression
description: "BREAK-GLASS diagnostic — use ONLY when an agent demonstrably mis-selects a text search (grep/Bash) for a STRUCTURAL question that Nexus answers better (callers, dependents, impact, equation/type edges). Detects graph-vs-text tool misselection and points at the fix (the always-on `.claude/rules/nexus-tools.md` rule + the deferred-MCP load gotcha). NOT for routine work — current models route freely."
---

# Behavioral Auto-Regression (break-glass diagnostic)

Use this ONLY when you observe a concrete tool-misselection regression: an agent ran a text
search (grep / `rg` via Bash) for a question whose native shape is **structural** — callers,
dependents, blast radius, call chains, equation/citation traceability, "who uses dependency
X" (aliased / late / `TYPE_CHECKING` imports) — i.e. a question `nexus-tools.md` routes to
Nexus. Routine sessions do NOT need this skill; the positive routing guidance lives in the
always-on rule `.claude/rules/nexus-tools.md`.

## Historical note — read this first

This skill was originally built to override a system-prompt directive (`ALWAYS use Grep for
search tasks`) by *reclassifying* code exploration as "not a search task". **That premise is
gone (verified 2026-06-14):** the standalone `Grep`/`Glob` tools were removed and no
"always-Grep" directive remains for Opus 4.8 / Sonnet 4.6 / Haiku 4.5 — all route freely. So
the reclassification trick is obsolete; do NOT apply it. What survives is the *diagnostic*:
tool-choice freedom does not guarantee tool-choice *correctness* — an agent can still fall
back on text-search habits, and (the opposite failure) can over-use Nexus as "compliance
theater" where a known `Read`/`grep` was correct.

## When to use

- An agent used grep/Bash for a structural question Nexus answers better (under-Nexus).
- A `mcp__nexus__*`-avoidance pattern recurs across a session.

## What to do

### 1. Confirm it is a real misselection
Any grep/Bash text search for a *relationship* question (callers, dependents, coverage,
equations, impact, type-usage) is a regression. A grep for a *literal string / comment /
known file* is **correct**, not a regression — don't over-correct into Nexus compliance
theater. Use the routing table in `.claude/rules/nexus-tools.md` as the gold standard.

### 2. Check the most common LIVE cause first — deferred MCP tools
The current real cause of Nexus-avoidance is that `mcp__nexus__*` tools surface as
**deferred**, and the agent treats deferral as unavailability. Fix: ONE
`ToolSearch("select:mcp__nexus__<name>")` call loads them. Confirm the agent knows this — it
is stated in `nexus-tools.md`.

### 3. Check the rule actually reaches the agent
- Main agent / any in-repo session: `.claude/rules/nexus-tools.md` auto-loads (always-on).
  If it isn't being followed, confirm the session is rooted in the repo and the rule isn't
  shadowed.
- A sub-agent that avoids Nexus: confirm its `AGENT.md` doesn't carry STALE override
  vocabulary that conflicts with the rule, and that it has the Nexus skills/tools available.
  (Stale AGENT.md Nexus-steering predates the rule and may now be redundant or contradictory
  — a known follow-up.)

### 4. Validate
Re-run the agent on a real structural task with "report every tool you used and why." Correct
behavior = Nexus (or its skills) for the structural part, grep/Read for the literal/known
parts. Zero grep on a pure structural task; non-zero grep on a literal-string task is fine.

See [reference.md](reference.md) for the original (now largely historical) procedure and the
validation results from the 2026-04 adoption. The `scripts/` probe and override-block
templates are retained for history; the override-block mechanism is superseded by the
`nexus-tools.md` rule.
