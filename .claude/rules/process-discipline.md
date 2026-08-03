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

## Mutation-testing an uncommitted file — never `git checkout` to revert

To prove a gate's teeth bite you mutate production code, run the gate (expect RED), then
revert. If the file has **uncommitted edits**, `git checkout -- <file>` / `git restore` /
`git stash` reverts to **HEAD — destroying your uncommitted work**, not to your pre-mutation
state.

- Revert a mutation by **monkeypatching in-process** (cleanest), or copy the file to a tmp
  path and mutate the copy. Reserve `git checkout`/`restore` for files you have **not** touched.
- Same data-loss family as the `.claude/*` checkout hazard (lessons L28): a `git checkout` on
  any path carrying uncommitted state is irrecoverable.

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
