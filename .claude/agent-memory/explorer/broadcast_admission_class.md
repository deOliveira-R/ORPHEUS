---
name: broadcast-admission-class
description: ERR-122 is one instance of a class — numpy-broadcast fibrewise products admit non-members; where the every-run guard is and is not (measured 2026-10-30 at 6aa26d00)
metadata:
  type: project
---

Measured 2026-10-30 at `6aa26d00` (re-verify before repeating). Detail and scripts:
`scratch/membership_covector/broadcast/census.md` (main checkout).

- **The class.** A stored coefficient multiplied into an element by numpy `*`/`einsum` takes the broadcast closure of shapes as its domain. It has two faces:
  - operand-side: ERR-122, ERR-055;
  - coefficient-side, axes placed by position: ERR-062, ERR-087, ERR-119, and the 2026-08-04 note on trailing vs leading.
- **Census.** 48 of 61 in-scope element verbs have no every-run guard.
  - Every numerics space, metric or functional verb is guarded only under `__debug__`, or not at all.
- **Where the every-run guards are.**
  - The only every-run exact parse is transport's `admit_array`, on the plain bindings.
  - `DiagonalOperator._check_shape` checks only rank, so `(8,2,1)` for `(8,2,12)` is silent in both modes.
- **Silent in both modes (asserts on and `-O`).** `compute_group_absorption_rate`, `Solution.reaction_rate_density`, `DiscreteMeasure.integrate`, the numerics `FrameAnalysis`.
- **Probe trap.** Compare against the INTENDED input, never a pre-expanded one; against the expanded reference, "right by accident" is a tautology.
  - Class silent outcomes as (i) the outer product, (ii) the implicit pullback x⊗1 (a different member), (iii) a non-member out.
  - `[M]` over 32 rows: 5 / 11 / 4 respectively.
- **The 48 split by ruling.** 39 of the 48 are outside the #664 raw-tier ruling; only the 9 debug-only rows are ruled.
- **Already solved once.** The by-label solution is ERR-119's `_broadcast_axes`. `_broadcast_leading` is its positional twin.

See [[space-membership-and-covector-census]].
