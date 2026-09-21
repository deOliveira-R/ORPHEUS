# Process discipline — evidence

Founding cases of [the process-discipline rule](../rules/process-discipline.md), moved (2026-09-20, when the rule came under generation) from the `> [M]` blockquotes and dated parentheticals of the hand-maintained `.claude/rules/process-discipline.md`. Each entry's first line names the clause it belongs to; the text below it is the original measurement, un-quoted and unchanged (its glyphs stay — it is history). The rule links here by heading anchor.

## Cases

### 2026-06-01 latent consumer

Clause: Bias toward completion

B.3, 2026-06-01: `BoundaryResidual` had no "consumer" only because the SN matvec computed the
boundary defect but mistyped it.

### 2026-06-21 seven campaigns in flight

Clause: Trust git for merge status

2026-06-21: ≈7 SN campaigns were mislabeled in-flight across the memory substrate; only one
branch was genuinely open — multiple sub-agents rediscovered this independently. Why it
matters: "resuming" merged work, or re-deriving a landed decision, is wasted effort, and a
stale "in-flight" claim makes a clean tree look unfinished.

### 2026-08-27 the 212-file sweep

Clause: The two git prohibitions are enforced by a hook; the riders are not

`[M]` 2026-08-27, un-weld P4.1a. `git add -A` swept **212** untracked `scratch/` files —
**74 670 lines** — into a commit about retiring a `coord` field, and it was pushed. The
rewrite that removed them then deleted all 212 from disk; they came back only from a
`cp -a` taken beforehand (745 files, restored byte-identical). Of the 212, **19 `.md`
memos are cited by path from tracked plans** — i.e. the sweep was not uniformly wrong,
which is precisely why reading the staged set matters rather than trusting a rule of
thumb about what `scratch/` is for. The streaming plan's own banner warns *"⚠ `scratch/`
is UNTRACKED — a `git clean` destroys it"*.

The sweep itself is refused since 2026-09-20 by `.claude/hooks/git-guard.py` (K3a of the
harness-context-budget campaign): `git add -A|--all|.` and a `git commit` while HEAD is
`main` exit 2 with the reason. The riders (read the staged set; `cp -a` aside before a
rewrite; content-exact check; re-point the hashes) are not enforced and stay in the rule.

### 2026-08-18 the battery the timeout killed

Clause: Mutation-testing an uncommitted file — the restore must be crash-safe

`[M]` 2026-08-18, nexus#88. An 8-arm battery over the full suite (~40 s each) hit the
2-minute limit; the `finally` did not run and `ast_analyzer.py` was left carrying a
deliberately-broken guard. `diff` against the copy-aside caught it on the next command.
Re-scoped to the two files that could redden, the same battery ran in ~30 s.

### 2026-08-18 twenty-nine of thirty red

Clause: After pushing, look at CI — and baseline a red one before adding to it

`[M]` 2026-08-18, nexus. CI had failed **29 of its last 30 runs**, across two
jobs, for over a day. I pushed **six** times into it without looking. The
baseline at `6e469e9` was **11** pyright errors; when I finally checked it
read **13** — my two, invisible in a number that had never been green. Had
the baseline been read at the first push, they would have been caught at the
commit that introduced them, where the fix is one narrowing.

The cost of the delay was not the fix (minutes) but the attribution: it took
a CI-log archaeology pass over an old run to establish which 11 were not
mine. ⚠ And note the direction — a long-red CI makes ADOPTING someone else's
failure the path of least resistance, because "it was already failing" is
true and exculpatory right up until it is your code.
