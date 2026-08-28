# Process discipline — bugs, completion, and issue hygiene

Workflow rules that operationalize Cardinal Rules 1 (Correctness) and 4 (GitHub issues are
the cross-session plan/log).

## Fix bugs immediately — never accumulate, never bury

When a bug is identified (test failure, code review, downstream cascade, L0 audit), the
**next code-touching action MUST be the fix**, before any new feature work. Do not defer
under "ship now, clean up later."

- Why (user, Phase 3B debrief): a known bug downstream means you can't distinguish "new
  bug" from "old bug propagated"; two simultaneous bugs can partially cancel so fixing one
  surfaces the other; every test below a known bug has reduced epistemic value; refactors
  over a buggy base propagate the bug into the new abstraction.
- Acceptable in-flight: bug known + fix in progress on the same branch. NOT acceptable: bug
  known + feature work proceeding in parallel. If the fix is genuinely large and the user
  explicitly chose ship-now-fix-later, that's deliberate — otherwise default to fix-now.
  When unsure, stop the feature work, ask, fix, resume. It is never too late to fix.

## Bias toward completion — finish known work while you hold the context

Don't defer work you KNOW you'll have to do to a future issue when the current session has
the **context advantage**. Deferral is for genuine cross-session/device handoff, not for
punting finishable work.

- Corollary: **"no current consumer" often means "consumer not yet wired," not
  "speculative."** Before concluding a type/feature is unused, hunt for LATENT consumers —
  quantities the code already computes but hasn't typed/routed yet. (B.3, 2026-06-01:
  `BoundaryResidual` had no "consumer" only because the SN matvec computed the boundary
  defect but mistyped it.)
- Separate the cheap, low-risk part (mint the type, do the rename — do it now) from the
  correctness-critical carve (wire it in — gate with test-architect/elegance review, still
  in-session if context allows).

## Trust git for merge-status — never a frozen memory claim

A memory or plan note that says a campaign is "NOT pushed / in-flight / on branch X / not
merged" is a **point-in-time snapshot**, not a standing fact. It is written mid-flight (often
in surgical-carve sessions) when it is true; the work merges in a later session and the note
is never updated — so the snapshot lies forward.

- **Before resuming, citing, or acting on** any "unmerged/in-flight" claim, reconcile against
  git: `git merge-base --is-ancestor <hash> HEAD` (ancestor ⇒ merged ⇒ it is done, not
  active). NEVER trust the frozen claim over git.
- **Update at merge-time:** a campaign note's terminal state is "merged @ `<hash>`", the same
  way a `Closes #NN` trailer closes an issue. After that it is archaeology — its lesson goes
  to `lessons.md`, its architectural milestone to the relevant theory page's Development-
  history changelog, and the note retires.
- **A plan's internal task numbers COLLIDE with real GitHub issue numbers.** A plan reading
  "Tasks #46–#52" may be that session's INTERNAL numbering while #46–#52 are unrelated live
  issues in another module. Never resolve "which issue tracks this?" from a plan's bare
  `#N` — read each candidate's title (`gh issue view N`). Same rule before *citing* an issue
  in code, docs, or a commit trailer: a docs pass that redirects a retired claim to "#183
  tracks this" is minting a new claim, and a closed or mis-numbered issue is the same
  defect class the pass was removing.
- Why it matters: "resuming" merged work, or re-deriving a landed decision, is wasted effort,
  and a stale "in-flight" claim makes a clean tree look unfinished. (2026-06-21: ≈7 SN
  campaigns were mislabeled in-flight across the memory substrate; only one branch was
  genuinely open — multiple sub-agents rediscovered this independently.)

## Never `git add -A` in this repo — stage explicitly, and READ the staged set

`git add -A` stages every untracked file. This repo carries a large, deliberately
**untracked** working directory (`scratch/` — probes, audit memos, gate logs; the streaming
plan's own banner warns *"⚠ `scratch/` is UNTRACKED — a `git clean` destroys it"*), and
`.gitignore` covers only two of its subdirectories. So the sweep is silent, enormous, and
lands inside a commit whose message describes something else entirely — which is what makes
it survive review: the subject line reads `refactor(geometry): …`, and nobody re-reads the
file list of a commit whose subject they just wrote.

⟹ **Stage by path** (`git add <paths>`), or `git add -u` when the change is confined to
tracked files. Before every commit, run `git status --porcelain` (not `-uno`) and read it.
The check is one command and it is the only thing standing between a working directory and
the history.

⭐ **And the sharper half, which the cleanup teaches: un-tracking a file by REWRITING
history DELETES it from the working tree.** `filter-branch`/`filter-repo` finish by checking
out the rewritten HEAD, and a file that is no longer in the new tree is removed from disk —
so the fix for an accidental commit destroys the very working files the accident preserved.
⟹ **`cp -a` the directory aside BEFORE the rewrite**, verify with `diff -rq` after, and
restore with `cp -Rn`. Same family as the `git checkout` clause below: any git operation
that writes the working tree is irrecoverable for untracked state.

⚠ Two riders for the rewrite itself, both cheap and both worth stating in the report:
- **Verify the rewrite is content-exact before pushing** —
  `git diff --quiet <old-tip> <new-tip> -- . ':(exclude)<removed-paths>'`. If that is
  empty, an already-run test gate **carries over without a re-run**, which is the
  difference between a 5-minute fix and a 60-minute one.
- **Every hash in the rewritten range changes**, so plans, memories and issue comments
  citing one now dangle (`plan-authoring` §7.1). Grep the tracked corpus for the old
  short hashes, re-point them, and leave a banner saying a dangling hash in that range
  predates the rewrite and should be re-found by commit SUBJECT.

> `[M]` 2026-08-27, un-weld P4.1a. `git add -A` swept **212** untracked `scratch/` files —
> **74 670 lines** — into a commit about retiring a `coord` field, and it was pushed. The
> rewrite that removed them then deleted all 212 from disk; they came back only from a
> `cp -a` taken beforehand (745 files, restored byte-identical). Of the 212, **19 `.md`
> memos are cited by path from tracked plans** — i.e. the sweep was not uniformly wrong,
> which is precisely why reading the staged set matters rather than trusting a rule of
> thumb about what `scratch/` is for.

## Mutation-testing an uncommitted file — never `git checkout` to revert

To prove a gate's teeth bite you mutate production code, run the gate (expect RED), then
revert. If the file has **uncommitted edits**, `git checkout -- <file>` / `git restore` /
`git stash` reverts to **HEAD — destroying your uncommitted work**, not to your pre-mutation
state.

- Revert a mutation by **monkeypatching in-process** (cleanest), or copy the file to a tmp
  path and mutate the copy. Reserve `git checkout`/`restore` for files you have **not** touched.
- Same data-loss family as the `.claude/*` checkout hazard (lessons L28): a `git checkout` on
  any path carrying uncommitted state is irrecoverable.
- ⭐ **The restore must be CRASH-safe, not merely exception-safe — a `finally` is not
  enough.** A battery that mutates, runs, and restores in a `try/finally` is safe against a
  failing test and unsafe against the thing that actually happens: the harness kills it. A
  `SIGTERM` at the 2-minute Bash timeout, a `pkill`, an interrupt — none of these let the
  `finally` complete, and the file is left MUTATED on disk, silently, in a working tree the
  next command will happily commit.
  ⟹ **copy the pristine file aside BEFORE the first mutation and verify against that copy
  after the run** (`diff -q`), rather than trusting the unwind. The copy is the only artefact
  that survives a kill. And budget the battery: N mutations × a full-suite run overruns the
  timeout by construction — scope the run to the files that can redden, and check the file's
  integrity FIRST after any battery that did not print its own completion line.
  > `[M]` 2026-08-18, nexus#88. An 8-arm battery over the full suite (~40 s each) hit the
  > 2-minute limit; the `finally` did not run and `ast_analyzer.py` was left carrying a
  > deliberately-broken guard. `diff` against the copy-aside caught it on the next command.
  > Re-scoped to the two files that could redden, the same battery ran in ~30 s.

## After pushing, look at CI — and a RED one must be baselined before you add to it

"`main` is always green" is a claim about CI, so it is only true if someone
reads CI. Pushing and moving on treats the gate as decoration.

⭐ **The non-obvious half, and the reason this is a rule rather than a
reminder: a CI that is ALREADY red cannot tell your regression from the
inherited one.** The run said `failure` before your push and `failure` after
it — identical output, zero information. This is the same defect the
`vv-principles` instrument rules describe (a reading that cannot change
carries nothing), except the instrument is the build.

⟹ On finding a red CI, **establish the baseline before fixing or adding
anything**: read the last run from before your work and count what failed.
Then your own contribution is a subtraction, and both halves are reportable —
*"the baseline was N, I added M, here is why"* — instead of a single number
that lets an inherited failure be quietly adopted or a new one quietly hidden.

> `[M]` 2026-08-18, nexus. CI had failed **29 of its last 30 runs**, across two
> jobs, for over a day. I pushed **six** times into it without looking. The
> baseline at `6e469e9` was **11** pyright errors; when I finally checked it
> read **13** — my two, invisible in a number that had never been green. Had
> the baseline been read at the first push, they would have been caught at the
> commit that introduced them, where the fix is one narrowing.
>
> The cost of the delay was not the fix (minutes) but the attribution: it took
> a CI-log archaeology pass over an old run to establish which 11 were not
> mine. ⚠ And note the direction — a long-red CI makes ADOPTING someone else's
> failure the path of least resistance, because "it was already failing" is
> true and exculpatory right up until it is your code.

## A refuted candidate is first-class output — record the structural REASON

When an investigation rejects a hypothesis, frame, design, or root-cause candidate, the
deliverable lists the rejected candidates **with the one-line structural reason each
failed** — never a bare "considered and rejected".

- **Why:** the reason is what stops the next session re-attacking a dead candidate. A
  rejection without its reason must be re-derived from scratch, at full cost — and the
  re-derivation usually costs more than the original, because the second attacker has less
  context about why it looked promising.
- **Who owes this:** `qa` (findings excluded from a report), `numerics-investigator`
  (hypotheses eliminated in the probe cascade), `explorer` (paths that dead-ended),
  `cross-domain-attacker` (the UNEXPLORED block), and the main agent in any plan that
  narrows a design space.

## Don't file issues for what you'll fix this session

GitHub issues are the cross-session/device handoff medium (Cardinal Rule 4). A finding that
will be remediated this session should be **fixed inline**, not filed — a back-to-back
open/close is just noise.

- After a QA/review batch: split findings into inline-fixable-now (fix them) vs
  cross-session followups (file with proper `module:`/`level:`/`type:` labels).
- Sub-agents (`qa`, `numerics-investigator`) sometimes default to "file an issue for
  tracking" — override that when in-session remediation is the plan.
