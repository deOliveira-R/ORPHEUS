---
harness:
  kind: rule
  budget_tokens: 1700
  brief: >-
    `grep` is ugrep, and an anchor inside an alternation group matches nothing, silently: use `\b…\b` or `-P` with a lookbehind; a sub-agent has no `ToolSearch`, so if Nexus is missing say so in NEEDS: and fall back to Bash.
---

# Code search in ORPHEUS — the shell's `grep`, the harness's tools and the project's graph cases, beside the bundled `nexus-tools` rule

What Nexus is, which question goes to the graph and which to `grep` or
`Read`, the skills to invoke instead of raw tools, the checks that are part
of the job, and the operational notes are the routing rule `nexus-tools`
beside this one: sphinxcontrib-nexus bundles it, `nexus setup` installs it,
`nexus setup --check` is its drift instrument, and it is never edited here
([the harness page](../harness.md), "The one-way flow"); its tool list is
the `nexus-guide` skill. This page holds what is ORPHEUS's: its cases, its
agents, its shell. Founding cases: [the evidence page](../evidence/code-search.md).

## The graph here

- The `explorer` agent (Nexus skills preloaded) is the exploration delegate
  for open-ended, multi-file work.
- The kinds of reference grep misses are enumerated with their instruments
  in `retirement-audit` A.2, and `dead_references` after a delete or a
  rename has its reason and its catcher in `retirement-audit` B.4
  ([case](../evidence/code-search.md#2026-04-restructuring-grep-misses)).
- **A stale graph announces itself** as a query returning an unexpected
  result: zero changes, an old module name. After a major file move or
  restructuring, assume it. check: rebuild first
  (`sphinx-build docs docs/_build/html`); the MCP server auto-reloads.

## `grep` here is ugrep, and an anchor inside an alternation group matches nothing, silently

`grep` in this environment is a shell function wrapping **ugrep**, not GNU or
BSD grep, and one construction fails in the worst possible way:
zero matches, exit 1, no message, indistinguishable from a clean tree. The
failing construction is an anchor (`^` or `$`) INSIDE an alternation group:
on a line containing `square Gram, while`, `grep -E '(^|[^a-z_])Gram'` returns
0, while `\bGram\b`, `-P '(?<![A-Za-z_])Gram'` and the anchor-free
`([^a-z_])Gram` return 1
([case](../evidence/code-search.md#2026-08-26-the-ugrep-fixture)). This is
exactly the idiom a retirement audit reaches for ("the symbol, but not preceded
by a letter or underscore"), so `coding-standards`' three-search audit is built
on greps like it.

- check: use `\b…\b`, or `-P` with a lookbehind, or drop the anchor; never put
  `^`/`$` inside a group. The mirror in `git grep -E`: there `\b` is the thing
  that does not exist (POSIX ERE has no word boundary), so
  `git grep -nE '\.(n_inner|n_outer)\b'` printed nothing over a tree with nine
  such reads ([case](../evidence/code-search.md#2026-09-17-git-grep-has-no-word-boundary)).
- check: for any completeness claim (a residual check, a "no consumers left"
  verdict, a done-when) re-run the pattern in Python (`re` + `pathlib.rglob`):
  the pattern is then unambiguous and the denominator can be stated (X2).
- check: a positive control before any negative (X1): one line asserting the
  pattern finds a member you already know exists
  ([L61](../evidence/lessons.md#l61-unvalidated-filter-clean): six false
  negatives in one session, two mechanisms — this one, and zsh eating quotes
  and backticks out of a double-quoted pattern, which at least prints
  `(eval): bad math expression` on a channel nobody reads).
- check: a probe is written to a FILE under a QUOTED heredoc (`<<'PY'`), its
  paths passed by environment variable, a backtick spelled `chr(96)`, and each
  pattern asserted against its own witness before it runs: an UNQUOTED heredoc
  (`<<PY`, chosen so a path interpolates) runs command substitution on every
  backtick pair in its body, so a pattern with backticks is silently rewritten
  before the interpreter sees it (`[M]` 2026-09-21: four markup patterns
  collapsed to "match any bold", 120 hits on clean prose; the same collapse
  the other way prints a clean 0). zsh also does not word-split an unquoted
  `$var`, and an unquoted word beginning with `=` (`echo ===`) is a command
  lookup that fails and aborts the whole compound, silently losing every grep
  sequenced after it: quote separators.
- check, for a probe that reads another program's output or runs from
  another tree: pytest's colour codes precede `FAILED`, so `grep -cE "^FAILED"`
  reads 0 on a red run (pass `--color=no`, or drop the anchor); and `python -c`
  puts the CWD at `sys.path[0]` ahead of `PYTHONPATH`, so a worktree probe run
  from the main tree imports the MAIN tree and prints HEAD's values for every
  commit (run a probe SCRIPT from outside the repository and print
  `module.__file__` first). `[M]` 2026-09, a nine-commit bisect: both, both
  flattering.
- tell: a confident, empty, wrong answer.

## The harness's search tools

- **The standalone `Grep` and `Glob` tools are gone**, so text search is
  `grep`/`rg` through Bash for every agent, and the routing rule is positive
  guidance, not an override of a default bias
  ([case](../evidence/code-search.md#2026-06-14-the-removed-grep-and-glob-tools)).
- **A sub-agent has no `ToolSearch`**
  ([case](../evidence/code-search.md#2026-08-19-a-sub-agent-has-no-toolsearch)):
  the routing rule states the fact and the check (a brief that depends on
  Nexus says what to do if it is missing); this page's `brief:` line is how
  that sentence reaches the three `omitClaudeMd` agents, which load no rule.
- **A worktree's graph** ([L22](../evidence/lessons.md#l22-worktree-not-main))
  and **the briefing's silence**
  ([case](../evidence/code-search.md#2026-08-16-the-briefings-silence)): the
  routing rule carries both checks; the cases are the project's.
