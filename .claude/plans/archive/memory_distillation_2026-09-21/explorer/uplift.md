# Uplift candidates — explorer distillation, 2026-09-21

Each candidate is a clause the digest keeps because no source that reaches
this agent carries it (`AGENT.md`, the two preloaded skills, the five brief
lines). Verdicts were made by reading the rule file named; the section is
quoted by its heading. Wording is proposed in the rule's own `check:` /
`tell:` form; the rule owner decides.

## U-1 `code-search` — zsh aborts a compound command on an unquoted `=`-word

- **Rule read:** `.claude/rules/code-search.md` § "`grep` here is ugrep …". Its
  zsh clauses cover quote/backtick eating inside a double-quoted pattern, the
  unquoted heredoc, and the unsplit `$var`; none covers `=cmd` expansion.
- **Proposed clause:** check: quote every separator in a batched Bash call
  (`printf 'NAME\n'`, `echo "---"`); an unquoted word beginning with `=`
  (`echo ===`) triggers zsh `=cmd` lookup, fails, and aborts the WHOLE
  compound command. tell: the greps after a separator print nothing and the
  only diagnostic is `== not found`, on a channel nobody reads.
- **Founding case:** explorer L-008 (archive §L-008).
- **Brief line:** the rule's `brief:` line need not change; the clause is
  environment-specific and rarely fatal, but it is silent.

## U-2 `instrument-doctrine` brief line — a membership question is parsed, not grepped

- **Rule read:** `.claude/rules/instrument-doctrine.md` § "X2. Every claim
  carries its population and its instrument", check: "A membership question
  is parsed with an AST, never grepped through a line window or a `head`."
  The brief line pasted to this agent carries X2's predicate / tree /
  exclusions / `k of N` / re-run-in-Python sentence and NOT this clause.
- **Proposed brief-line addition:** "a membership or consumer question is
  answered by an AST pass (calls, receivers, keyword sites), and a line-grep
  is reported only as the prose column — a public name's docstring fame
  reads as consumption."
- **Founding cases:** explorer L-038 (3 of 6 public `symmetry.py` functions:
  0 executable callers behind 4–5 `:func:` lines), L-043 (148 AST keyword
  sites vs 183 grep lines), L-045 (17 calls + 1 docstring line; the real
  member was a `def sweep` surrogate), L-037 (`Assign` ≠ `AnnAssign`).
- **Where the digest keeps it:** L-038, marked `[uplift]`.

## U-3 `articulation` brief line — `[M]` certifies a measurement, not the sentence

- **Rule read:** `.claude/rules/plan-authoring.md` § "§2 Mark the epistemic
  status of every claim", item [M]-SCOPE: "`[M]` certifies that a
  measurement happened; it does not certify that the measurement answers
  the sentence it sits in." `plan-authoring` is path-scoped
  (`.claude/plans/**`) and has no `brief:` line; `articulation`'s brief line
  carries the marker vocabulary and could carry this one sentence.
- **Proposed brief-line addition:** "`[M]` on an inherited claim certifies
  that some measurement answered some question; a `[M]` on a NEGATIVE
  (absent, discarded, no consumers) is re-measured against the question at
  hand before it is built on."
- **Founding case:** explorer L-025 (`[M] "solve_sn discards the solver it
  builds"` — true of the returned `Solution`, false at the call site the
  plan reused it for).
- **Where the digest keeps it:** L-025, marked `[uplift]`.

## U-4 `workflows` invariant 4 (the brief) — a sub-agent dispatched into a moving tree

- **Rule read:** `.claude/rules/workflows.md` § "Invariants", item 4 lists
  what every brief carries; `.claude/rules/process-discipline.md` § "Trust
  git for merge status" covers cross-session staleness only. Nothing that
  reaches a Support agent says what to do when the main session edits the
  audited files DURING the dispatch.
- **Proposed clause (for the brief template or invariant 4):** check: a
  Support agent auditing a subsystem with uncommitted edits opens with
  `git status --short` and `git diff --stat`, re-runs both at close, re-runs
  verbatim every search whose emptiness is a finding, runs `git ls-files
  --error-unmatch` on each file it calls landed, and tags every cited file
  "(at HEAD)" or "(untracked, in-flight)". tell: a "zero consumers" verdict
  or a "not yet built" premise in a report written while the carve ran.
- **Founding cases:** explorer L-007 (a census moved 1532 → 1552 between
  two greps), L-012 (eight files clean → modified during one dispatch; a
  "zero consumers" verdict flipped by a rule landed mid-dispatch).
- **Generalises to:** qa and test-architect dispatched mid-carve.

## U-5 The brief template — the brief's own data are claims

- **Rule read:** `plan-authoring` §1 PRECEDENT, §2 [M]-SCOPE, CHECKLIST-HALO,
  A-BRIEF'S-METHOD-IS-A-CLAIM cover the AUTHOR's duties; no brief line tells
  the RECIPIENT that a brief's timeline, count, exemplar or `Class.attr
  (file:line)` citation is verified before it is built on.
- **Proposed sentence for the "Rules that apply to you" generator (or
  `articulation`'s brief line):** "every datum the brief states — a hash's
  date, a count, an exemplar's behaviour, a `Class.attr (file:line)` — is a
  claim; verify it by one command before building on it, and report the
  discrepancy as a finding."
- **Founding cases:** explorer L-020 (the brief's "landed AFTER the audit"
  was backwards), L-023 (the "deliberate `=50`" helper defaulted to 4000),
  L-027 ("37 sites, all in tests/" reproduced under no convention),
  L-045 (the cited line sat in a different class).
- **Where the digest keeps it:** meta-lesson M-1; also offered for
  `AGENT.md` OP5 in `agent_md.md` item 1, which would make U-5 unnecessary
  for this agent.
