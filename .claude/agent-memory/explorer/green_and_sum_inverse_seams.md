---
name: green-and-sum-inverse-seams
description: Post-S2 shape of the sum inverse / GreenOperator / seeded_inverse seams — who builds Green, where the seed is dropped, the .H tautology trap (2026-10-20, main@2375992d)
metadata:
  type: project
---

Recorded at main@2375992d (post-S2, pre-S4); re-verify before repeating.

- Production builds NO GreenOperator for a solve: 1 of 28 `GreenOperator(` calls is in orpheus/, inside
  `GreenOperator.dual`. The Strategy realises (M−N)⁻¹ at the driver layer directly (`_within_group_si`,
  `_within_group_krylov`, `KEigenvalue`, `solve_sn_adjoint_fixed_source` = `SourceIteration(seeded_inverse(M.H), N.H)`).
  `Splitting.explicit` is already N-as-gains-tuple, the same shape as `GreenOperator._gains`.
- Multi-term `OperatorSum.inverse()` always raises; 0 production `except NotInvertible`, so every runtime
  multi-term receiver is a refusal test. 9 tests pin the message fragment `GreenOperator(M, N)`.
- `OperatorSum` has no `solve`, so `InverseOperator(sum)` cannot realise itself (its apply = inner.solve).
- Seed paths of (αA)⁻¹: direct `_apply_scaled` seeds α·x₀; `seeded_inverse(αA)` wraps in `_SeededExactApply`
  and drops it — a family-membership effect, discriminated only with a wrap-member recorder (a plain test
  double is dropped on BOTH arms, so it cannot show the disagreement).
- Trap: making `G.H` = Green of the `.H` legs re-creates the override removed 2026-10-18; the honesty-gate
  witness `_sn_green` is exactly that formula (D.14 tautology). `GreenOperator(self.inner.H)` is not even
  constructible (ctor takes M, N).
- G0 gate: operator rank < iteration < green_operator; a realisation imported by operator.py is a violation,
  and PENDING may not grow.
- Memory-write path: the hook accepts the MAIN checkout's `.claude/agent-memory/explorer/`, not the worktree copy.

Related: [[operator-composite-algebra-seam]], [[adjointability-flag-semantics]].
