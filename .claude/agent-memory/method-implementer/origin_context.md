---
name: Origin context
description: The empirical case study that motivated minting the method-implementer agent on 2026-05-02
type: project
---

The method-implementer agent was minted on 2026-05-02 from GitHub
issue #147. The motivating empirical evidence is the Plan 2 follow-on
A1+A2 work on `feature/peierls-greens-function` (commits
`4754993..6a61ed8`), where the recurring task pattern was:

> Plan → SymPy-prove load-bearing algebra → implement prototype →
> build structurally-independent reference (or dispatch
> literature-researcher) → cross-verify → if disagreement, isolate
> the bug → iterate → update memos → dispatch archivist for Sphinx

This loop ran ad-hoc in the main thread before the mint and produced
four diagnostic events that a specialised agent with the right
preloaded skills would have caught without user intervention:

1. Forgot to load the `vv-principles` skill at session start despite
   CLAUDE.md instructing it. (Caught by user.)
2. Calibrated L0 thresholds (`std/mean > 0.5`) against the buggy
   prototype's own output instead of against physics expectations.
3. Took the V_α1 numerical pass (closed sphere k_eff = k_inf to
   machine precision) as evidence the operator was correct. It
   wasn't — V_α1's algebraic closure (ψ_surf solving its own
   fixed-point) cancels the L_first variable identically, masking a
   forward-vs-backward trajectory bug. Only the PS-1982
   structurally-independent cross-check exposed it (6 % disagreement
   → bug fix → < 1e-4).
4. Assumed µ → −µ symmetry without making the assumption visible.
   The symmetry holds for the rank-1 closed-sphere isotropic
   eigenmode but breaks for vacuum BC.

These events are the canonical training signal for the agent. They
inform the **algebra-of-record** discipline (closed-form 1A and
semi-analytical 1B should both be implemented when the problem admits
them — 1A catches algebraic errors, 1B catches discretization errors
that 1A's algebraic cancellation might mask) and the **bias steering
to dispatch literature-researcher EARLY** (don't reconstruct what
you can read).

## Original 8 bias-steering lines — institutional homes

The original issue listed 8 bias-steering rules. After session
2026-05-02 they were absorbed into project-level skills, leaving only
one as procedural AGENT.md material:

| Original                                                | Lives in (after mint)                           |
| ------------------------------------------------------- | ----------------------------------------------- |
| 1. Never claim verification on closure-by-construction  | `vv-principles` (anti-patterns)                 |
| 2. Threshold calibration from physics, not prototype    | `vv-principles` (anti-patterns)                 |
| 3. Forward-vs-backward trajectory is Failure Mode 6     | `numerical-bug-signatures` (signature catalog)  |
| 4. Symmetry assumptions must be visible                 | `vv-principles` (claim-scope discipline)        |
| 5. Log-singular kernels need named recipes              | `algebra-of-record` (or future skill)           |
| 6. Cross-check disagreement → hypothesis ordering       | `vv-principles` / `algebra-of-record`           |
| 7. Structural independence beats agreement              | `vv-principles` + `algebra-of-record`           |
| **8. Dispatch literature-researcher EARLY**             | **AGENT.md** — procedural workflow rule (kept)  |

Lines 1, 2, 3 are flagged as "outstanding before mint" in the
follow-up issue (search for `method-implementer skill follow-up` in
GitHub). The mint proceeded without waiting because the AGENT.md is
lean and the skills auto-deliver new content on next preload.

## Pointers

- Originating issue: GitHub #147 (CLOSED 2026-05-02)
- Closeout memo (peierls A1-A3 + Plan-(b)):
  `.claude/agent-memory/numerics-investigator/peierls_greens_phase1_closeout.md`
- Variant α decision memo:
  `.claude/agent-memory/numerics-investigator/peierls_greens_variant_alpha_decision.md`
- Algebra-of-record skill:
  `.claude/skills/algebra-of-record/SKILL.md`
- V&V skill:
  `.claude/skills/vv-principles/SKILL.md`
- Branch: `feature/method-implementer-agent` (mint commit)
