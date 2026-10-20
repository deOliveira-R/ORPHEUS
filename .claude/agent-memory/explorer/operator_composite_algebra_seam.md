---
name: operator-composite-algebra-seam
description: Where OperatorSum/Product/Scaled are minted and read (2026-10-20, 7621fc5f) — dunders mint ~99%, 6 production leg readers, fusion is operand-type dispatch, (A+B).H wraps the sum
metadata:
  type: project
---

Shape of the composite operator algebra at `7621fc5f` (S2 census, full tables in
`scratch/numerics_untangling/packages/S2_census.md` while it lives).

- Minting: >= 99 % of production sums/scaled operators come from the `LinearOperator` dunders,
  not constructor calls (runtime spy). Direct production ctor calls: 9.
- Fusion is OPERAND-TYPE dispatch in three dunders: `StreamingOperator.__add__` (L + C only),
  `StreamingCollisionOperator.__sub__` (- B_lower), `BulkAnalysisOperator.__matmul__` (@ SweepOperator).
  The fused results are the only 3 subclasses (2 sums, 1 product); `_left_spine_terms` treats them as
  leaves via `type(node) is OperatorSum`.
- Production leg readers: 6 sites (green `_left_spine_terms`/`_negated`; named accessors on the 3
  subclasses). No production walker over the operator tree; roles are joined at construction.
- `(A+B).H` = `AdjointOperator(A+B)`; nothing distributes the adjoint. `A.H.H is A` is pinned through the pencil.
- `StreamingOperator` is a non-frozen dataclass: field `__eq__`, `__hash__ = None` (unhashable term).
- A string-spelled PRIVATE-slot pin exists: `getattr(kernel, "_a", None)` in `test_psi_half_coupling.py`
  (slow test) — invisible to a property spy and to a `.a` grep.

**Why:** S2 (ruling 23) turns the sum into an unordered formal linear combination; these are the seams.
**How to apply:** for any later composite change, start from the dunders and the fusion dispatch, not from
`OperatorSum(` grep; re-verify line numbers (operator.py moves every wave). Instrument: a `-p` pytest
plugin wrapping the leg properties with a caller-frame counter ([[census-predicates-bound-reference-and-activation-traceback]]).
