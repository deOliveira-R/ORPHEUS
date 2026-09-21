# Nexus & code-exploration tools — evidence

Founding cases of [the nexus-tools rule](../rules/nexus-tools.md), moved (2026-09-20, when the rule came under generation) from the dated paragraphs, the fixture table and the notes of the hand-maintained `.claude/rules/nexus-tools.md`. Each entry's first line names the clause it belongs to; the text below it is the original measurement, unchanged (its glyphs stay — it is history). The rule links here by heading anchor.

## Cases

### 2026-04 restructuring grep misses

Clause: the preamble (why a graph beats grep for structure)

During the 2026-04 package restructuring, over-reliance on grep caused repeated missed
imports that one `mcp__nexus__impact` query would have caught. grep matches *text*; it
misses relationships — inline imports, `TYPE_CHECKING` blocks, late imports inside
functions, aliased imports (`from numpy import linalg as la`), re-exports, and docstring
references. Nexus captures all of these as graph edges.

### 2026-08-26 the ugrep fixture

Clause: `grep` here is ugrep, and an anchor inside an alternation group matches nothing, silently

`[M]` 2026-08-26. `grep` in this environment is a shell function wrapping
**ugrep 7.5.0** (`ARGV0=ugrep … -G --ignore-files --hidden -I --exclude-dir=…`),
not GNU or BSD grep (`[M]` 2026-09-20: the same fixture on ugrep 7.8.4 reproduces the
silent zero). Its regex dialect differs in at least one way that matters, and the failure
mode is the worst possible one: **zero matches, exit 1, no error message** —
indistinguishable from a clean tree.

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

### 2026-09-17 git grep has no word boundary

Clause: the ugrep clause's first check (the mirror in `git grep -E`)

The mirror in `git grep -E`: there `\b` is the thing that does not exist (POSIX ERE has no
word boundary), so `git grep -nE '\.(n_inner|n_outer)\b'` over a tree with nine such reads
printed NOTHING (`[M]` 2026-09-17, step 3 U6) — the same silent-zero shape, one binary
over; the countermeasure (re-run any completeness claim in Python with `re` +
`pathlib.rglob`, and validate the filter against a positive control) covers both.

### 2026-08-19 a sub-agent has no ToolSearch

Clause: Operational notes — deferred tools

If `mcp__nexus__*` surface as deferred, ONE `ToolSearch("select:mcp__nexus__<name>")` loads
them — deferral is NOT unavailability. ⚠ **A sub-agent has no `ToolSearch` tool**, so this
recovery path does not exist for it — `[M]` 2026-08-19, a sub-agent probe reported 45
`mcp__nexus__*` tools loaded eagerly and **no `ToolSearch` at all**. A sub-agent that finds
Nexus genuinely absent cannot recover; it must say so and fall back to `Bash` (grep, or
`python -c "from sphinxcontrib.nexus.export import load_sqlite"` against the graph DB).
⟹ **when a dispatch depends on Nexus, say in the brief what to do if it is missing** —
otherwise the agent improvises silently, and its report cannot be told apart from a
grep-derived one. This is the most common cause of an agent silently avoiding the graph.

`[M]` 2026-09-21, two dispatches (qa, method-implementer): no `ToolSearch`, corroborated; 40 `mcp__nexus__*` tools enumerated from the definitions, against the 45 above, and which half moved (the extension's version or the earlier count's predicate) is unmeasured.

### 2026-08-16 the briefing's silence

Clause: Operational notes — `session_briefing` warns about indexed files, not the branch

`session_briefing` warns when files the graph INDEXES have changed — not when the branch
differs. Those are different questions, and reading the second into the first makes the
warning look broken: an ordinary ff-merge-and-delete leaves the graph describing the
checkout exactly while the branch name has moved on (`[M]` 2026-08-16: 25 files differed
from the build commit, **0 of them indexed**, and the briefing was right to stay quiet).
⟹ *silence means the indexed sources match*, not "same branch".

### 2026-06-14 the removed Grep and Glob tools

Clause: Operational notes — the standalone `Grep` and `Glob` tools are gone

NOTE (2026-06-14, re-verified 2026-08-19): the standalone `Grep`/`Glob` tools were removed
and the "always-Grep" system-prompt directive is gone — models route freely. `[M]`
re-measured on Opus 5: a sub-agent's tool list carries no `Grep` and no `Glob`. This rule
is *positive routing guidance*, not an override of a default bias.

The note went on to defer the rule's exact "dose" (how much steering each model needs to
route correctly without Nexus compliance theatre) to a tool-routing ablation study
(`.claude/plans/tool_routing_ablation_study.md`, design v0.1 of 2026-06-14, never
executed). Retired 2026-09-20 as obsolete by user ruling: tool routing has since been
verified working in practice, so the study has no question left to answer; the plan file
was deleted in the same commit that brought this rule under generation.
