---
name: issue-196-phase-g-step3-4a-closeout
description: Phase G Step 3+4.a — ScatteringOperator.foldable_part / residual_part / foldable_sigma; the S = S_foldable + S_residual data split (pure algebra, no wrapper class). Predecessor to the committed step3_4bi_v2 (leaves + is_foldable_into_sigma_r predicate). LANDED on origin/main.
metadata:
  type: project
---

# Issue #196 Phase G Step 3+4.a — ScatteringOperator foldable/residual split

`refactor/sn-operator-algebra` commit `1fcd34c`, building on `2df0a61`. **LANDED
on origin/main.** Scattering-layer-only extension: three new methods on the
existing `ScatteringOperator` exposing `S = S_foldable + S_residual`. Pure
additive — `apply()` byte-for-byte unchanged. The committed
[[issue_196_phase_g_step3_4bi_v2_closeout]] is the consumer (the `StreamingOperator` +
`CollisionOperator` leaves + `is_foldable_into_sigma_r()` predicate); this note
is the data-API predecessor it builds on.

## What landed (`orpheus/sn/scattering.py`, +166 LoC)

- **`foldable_part() -> ScatteringOperator`** — sibling with
  `scattering_order=0`, `Y=None`, diagonal-only `sig_s[mid][0]`
  (`np.diag(np.diag(...))`, length-1 list), `sig_s0[mid]` = that diagonal,
  `sig2[mid] = zeros`. The within-group self-scatter only.
- **`residual_part() -> ScatteringOperator`** — sibling carrying the P0
  off-diagonal (`p0 - np.diag(np.diag(p0))`), every Pℓ≥1 block verbatim, `sig2`
  verbatim (n,2n is unconditionally residual), `Y is self.Y`.
- **`foldable_sigma() -> dict[int, np.ndarray]`** — `{mid: diag(sig_s[mid][0]).copy()}`,
  the per-material `(ng,)` self-scatter XS. The cheaper data view that
  substep 3+4.b's cell-balance fusion (`σ_r = σ_t − foldable_sigma()[mid][g]`)
  actually consumes.

**Foldability criterion (the math)**: only P0 within-group self-scatter is
foldable into the sweep's removal cross-section, because `Y_ℓ^m(Ω_n)` makes
Pℓ≥1 self-scatter direction-dependent — it cannot collapse to a scalar diagonal
in the cell-balance denominator. The algebraic identity
`S.apply ≈ foldable_part().apply + residual_part().apply` holds at `rtol=1e-14`
(FP-non-associativity bounded; sum-of-two-applies vs single-apply).

## Durable lesson — the split is DATA, not TYPE (Pattern 6 + Pattern 4 in tandem)

The temptation to introduce `ScatteringSplit` / `WithinGroupScattering` /
`ResidualScattering` wrapper classes was real (neat typed objects expressing
"this is the foldable half"). Pattern 6 (defer abstraction) + the user-locked
decision "use pure algebra, no `WithinGroupOperator` class" closed it. The same
`ScatteringOperator` class with different per-material arrays IS the split. The
caller `A_wg = L + C − S.foldable_part()` consumes the SAME operator type as the
parent `S`; the consumer cares about `apply()`, not which half it holds. **An
`is_foldable` boolean attribute would be a Pattern-4 violation**: there is no
value of such a flag that expresses "the sibling produced by `foldable_part()`"
better than the type itself does. The committed step3_4bi_v2 followed through —
its `is_foldable_into_sigma_r()` is a PREDICATE on the operator (a question you
ask), not stored state.

## Anti-mutation discipline (the purity contract)

Every new array in `foldable_part`/`residual_part` is a fresh numpy expression
that does NOT alias the parent (`np.diag(np.diag(...))`, `np.zeros_like`,
`p0 - np.diag(np.diag(p0))` all allocate). Cross-instance sharing is limited to
read-only `apply()` inputs (`Y`, `weights`, `cells_by_mat`, residual's `sig2`
and Pℓ≥1 blocks). `foldable_sigma` returns `.copy()` so consumers can mutate
freely. Purity tested by snapshot-before / call / assert-array-equal-after. The
purity tests verify `np.array_equal` (not `is`-identity), leaving future
implementations free to add `lru_cache` / `cached_property` without breaking the
contract.

## Verification

31 new tests (6 classes: `TestFoldablePart` 9, `TestResidualPart` 8,
`TestFoldableSigma` 4, `TestAlgebraicIdentity` 6, `TestPurity` 3, + 1 fixture
`solver_2g_p1_n2n` stressing cross-group P0 + non-zero P1 + non-zero n,2n).
48/48 in `test_scattering_operator.py`. 11 regression snapshots bit-identical at
rtol=1e-12 BY CONSTRUCTION (`apply` unchanged; the 3 new methods have ZERO
production consumers outside `scattering.py` — grep-verified). No new ERR-NNN.

**File-path note**: brief specified `tests/sn/scattering/test_scattering_operator.py`;
the on-disk path is `tests/sn/test_scattering_operator.py` (no `scattering/`
subdir). Tests added to the existing file. A brief-text-vs-disk-state mismatch,
not a coverage gap — worth checking the actual on-disk path before trusting a
brief's stated test location.
