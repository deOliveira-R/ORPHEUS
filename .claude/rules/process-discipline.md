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

## Don't file issues for what you'll fix this session

GitHub issues are the cross-session/device handoff medium (Cardinal Rule 4). A finding that
will be remediated this session should be **fixed inline**, not filed — a back-to-back
open/close is just noise.

- After a QA/review batch: split findings into inline-fixable-now (fix them) vs
  cross-session followups (file with proper `module:`/`level:`/`type:` labels).
- Sub-agents (`qa`, `numerics-investigator`) sometimes default to "file an issue for
  tracking" — override that when in-session remediation is the plan.
