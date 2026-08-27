# Nexus & code-exploration tools

ORPHEUS ships **Nexus** (`sphinxcontrib-nexus`) — the single knowledge graph unifying
code structure (call graphs, imports, inheritance, type annotations) with documentation
structure (equations, cross-references, citations, theory pages). It runs as an MCP
server (tool list: `nexus-guide` skill); the graph rebuilds on every `sphinx-build` and
the server auto-reloads when the DB changes (v0.4.3+).

**Nexus is this project's structural code-intelligence layer.** It answers relationship
questions that text search fundamentally cannot — and the project `explorer` agent (Nexus
skills preloaded) is the exploration delegate for open-ended, multi-file work.

**Why a graph beats grep for structure:** grep matches *text*; it misses relationships —
inline imports, `TYPE_CHECKING` blocks, late imports inside functions, aliased imports
(`from numpy import linalg as la`), re-exports, and docstring references. Nexus captures
all of these as graph edges. (During the 2026-04 package restructuring, over-reliance on
grep caused repeated missed imports that one `mcp__nexus__impact` query would have caught.)

You have **freedom of tool choice** — route by what the question actually is:

| Question | Tool | Why |
|---|---|---|
| Callers / dependents / call chains / blast radius | Nexus `callers`, `impact`, `processes` | graph traversal; text can't follow edges |
| Equation / citation traceability | Nexus `provenance_chain` | links code ↔ equations |
| Verification coverage | Nexus `verification_audit` | maps equation → code → test |
| Docs referencing symbols that no longer exist | Nexus `dead_references` | renders as plain text; no build warning |
| Which tests must re-run after this change | Nexus `retest` (pass `run=`) | with a coverage capture it answers from what EXECUTED; the static cone alone has 12-15 % recall |
| Does the test that CLAIMS to verify this actually run it | Nexus `verification_audit` (pass `run=`) | a `verifies` marker is authored and unfalsifiable until evidence adjudicates it |
| Failing-test diagnosis | Nexus `trace_error` | walks the call graph to the suspect equation |
| Safe rename / refactor | Nexus `rename`, `impact` | finds references by graph, not text |
| "Who uses dependency X" (incl. aliased imports) | Nexus `graph_query` / `type_uses` | grep misses aliased / late / `TYPE_CHECKING` imports |
| Structural smells (clones, dead code, missing types) | Nexus `twin_paths`, `dead_functions`, `discriminations`, `native_place`, `protocol_conformers` | whole-graph sweeps |
| What actually RAN (hotspots, real dispatch) | Nexus `runtime_*` | dynamic overlay; static graph can't see it |
| A file position (LSP result, stack trace) | Nexus `node_at` | position → graph node |
| Literal text / regex / config values | `grep`/`rg` via **Bash** | finds raw strings |
| TODO / FIXME / inline comments | `grep` via **Bash** | Nexus doesn't index comments |
| Known file / known symbol body | **Read** (or `find` via Bash) | don't rediscover what you already know |
| Unknown symbol location | either — Nexus `query` or `grep` | your call |

Over-using Nexus where a plain `Read` or `grep` was correct is as much a misselection as
grepping for a relationship question. Do not perform compliance theater.

**Users describe symptoms, not tools.** "We keep changing these classes in lockstep",
"two people built this separately", "things live in surprising places", "the docs feel
out of date" are all graph questions — route them to `protocol_conformers` / `twin_paths`
/ `native_place` / `dead_references` respectively rather than reading files until a
pattern appears.

**Some checks are part of the job, not a request.** After you delete or rename anything,
run `dead_references` before calling it done — green tests do not cover prose, and a dead
documentation reference produces no build warning at any severity, so nothing else will
catch it. Before a release, and for any "health check" or onboarding review, sweep the
smell family (`twin_paths`, `discriminations`, `native_place`, `protocol_conformers`,
`dead_functions`) alongside `dead_references` and `staleness`.

**Invoke the Nexus *skills*, not raw MCP tools** — they encode the complete workflows
(`nexus-exploring`, `nexus-impact`, `nexus-debugging`, `nexus-refactoring`,
`nexus-verification`, `nexus-elegance`, `nexus-guide`).

## ⛔ `grep` here is **ugrep**, and one common construction fails SILENTLY

`[M]` 2026-08-26. `grep` in this environment is a shell function wrapping
**ugrep 7.5.0** (`ARGV0=ugrep … -G --ignore-files --hidden -I --exclude-dir=…`),
not GNU or BSD grep. Its regex dialect differs in at least one way that matters,
and the failure mode is the worst possible one: **zero matches, exit 1, no error
message** — indistinguishable from a clean tree.

**The failing construction: an anchor (`^` or `$`) INSIDE an alternation group.**
Isolated on a 3-line fixture containing `square Gram, while`:

| pattern | matches | |
|---|---|---|
| `grep -E 'Gram'` | 1 | ✅ |
| `grep -E '(^\|[^a-z_])Gram'` | **0** | ⛔ **silent false negative** |
| `grep -E '([^a-z_])Gram'` | 1 | ✅ (anchor removed) |
| `grep -P '(?<![A-Za-z_])Gram'` | 1 | ✅ (PCRE lookbehind, needs `-P`) |
| `grep -E '\bGram\b'` | 1 | ✅ |

⚠ **This is exactly the idiom a retirement audit reaches for.** *"the symbol, but
not preceded by a letter or underscore"* is how you separate `gram` from
`programs`, or a retired module name from a surviving attribute of the same
spelling — and writing it the natural way returns a confident, empty, wrong
answer. `coding-standards`' three-search audit is built on greps like this.

⟹ **Use `\b…\b`, or `-P` with a lookbehind, or drop the anchor. Never put `^`/`$`
inside a group.** And for any *completeness* claim — a residual check, a "no
consumers left" verdict, a done-when — **re-run it in Python** (`re` +
`pathlib.rglob`) rather than the shell: the pattern is then unambiguous, and it
is the only way to state a denominator you can trust.

⭐ **Validate the filter against a POSITIVE CONTROL before trusting a negative.**
One line — assert the pattern finds a member you already know exists. Without it,
a broken filter and a clean tree print the same thing, and the broken one reads
as *"nothing to do"*. See `.claude/lessons.md` L61 (six false negatives in one
session, two distinct mechanisms: this one, and zsh eating quotes/backticks out
of double-quoted patterns — that second one at least prints `(eval): bad math
expression` on a channel nobody reads).

**Operational notes**

- **Deferred tools — ⛔ and the escape hatch is MAIN-AGENT-ONLY.** If `mcp__nexus__*`
  surface as deferred, ONE `ToolSearch("select:mcp__nexus__<name>")` loads them —
  deferral is NOT unavailability. ⚠ **A sub-agent has no `ToolSearch` tool**, so this recovery path
  does not exist for it — `[M]` 2026-08-19, a sub-agent probe reported 45 `mcp__nexus__*`
  tools loaded eagerly and **no `ToolSearch` at all**. A sub-agent that finds Nexus
  genuinely absent cannot recover; it must say so and fall back to `Bash` (grep, or
  `python -c "from sphinxcontrib.nexus.export import load_sqlite"` against the graph DB).
  ⟹ **when a dispatch depends on Nexus, say in the brief what to do if it is missing** —
  otherwise the agent improvises silently, and its report cannot be told apart from a
  grep-derived one.
  This is the most common cause of an agent silently avoiding the graph.
- **Stale graph:** rebuild Sphinx first (`sphinx-build docs docs/_build/html`); the MCP
  server auto-reloads.
- **Git worktrees (L22 hazard):** the session's MCP server was launched against the MAIN
  checkout's graph, so every query answers from the wrong branch until you switch. Build
  Sphinx inside the worktree, then `mcp__nexus__use_workspace(<worktree root>)`;
  `mcp__nexus__workspaces` lists every checkout and its graph.
  ⚠ **`session_briefing` warns when files the graph INDEXES have changed — not when the
  branch differs.** Those are different questions, and reading the second into the first
  makes the warning look broken: an ordinary ff-merge-and-delete leaves the graph
  describing the checkout exactly while the branch name has moved on (`[M]` 2026-08-16:
  25 files differed from the build commit, **0 of them indexed**, and the briefing was
  right to stay quiet). ⟹ *silence means the indexed sources match*, not "same branch".

> **NOTE (2026-06-14, re-verified 2026-08-19):** the standalone `Grep`/`Glob` tools were
> removed and the "always-Grep" system-prompt directive is gone — models route freely.
> `[M]` re-measured on Opus 5: a sub-agent's tool list carries no `Grep` and no `Glob`.
> This rule is *positive routing guidance*, not an override of a default bias. Its exact
> dose (how much steering each model needs to route correctly without Nexus "compliance
> theater") is to be calibrated by the tool-routing ablation study
> (`.claude/plans/tool_routing_ablation_study.md`) — revisit this file with those results.
