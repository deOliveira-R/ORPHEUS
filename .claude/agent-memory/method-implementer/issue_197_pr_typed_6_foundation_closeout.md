# Issue #197 PR-TYPED-6 (foundation) — `cell_balance_for_streaming` helper + `DiamondDifference.residual` delegation

**Branch**: `refactor/sn-operator-algebra` from tip `e09b9f8`.
**Date**: 2026-05-16.
**Status**: STAGED, NOT COMMITTED. Working tree contains the foundation
piece of Option β. The full Option β scope (matvec consolidation +
helper retirement + k-traversal flip + snapshot regen) is **NOT
shipped in this dispatch** — see §3 for the remaining-work itemisation.

This memo documents the principled atomic-foundation contribution. It
is shippable as-is. The remaining Option β substeps depend on this
foundation and should be sequenced as separate atomic PRs per the
Phase G pattern (Steps 1, 2, 2.5, 2.5b, 3+4.a, 3+4.b.i — each small,
each individually verified).

---

## §1 What this dispatch shipped

### §1.1 New shared helper `cell_balance_for_streaming` (§B.1 of brief)

**File**: `orpheus/sn/spatial/cell_balance.py` (~140 LoC added).

The vectorized cell-balance helper Pattern-2-shared between
`DiamondDifference.residual` (n_mask=1) and the future
`StreamingOperator.apply` consolidation (n_mask=N).

**Signature**:

```python
def cell_balance_for_streaming(
    *,
    abs_mu: np.ndarray,             # (n_mask,)
    A_downstream: np.ndarray,       # (n_mask,)
    A_total: np.ndarray,            # (n_mask,)
    dA_w: np.ndarray,               # (n_mask,)
    c_in: np.ndarray,               # (n_mask,)
    c_out: np.ndarray,              # (n_mask,)
    total_xs: np.ndarray,           # (ng,)
    volume: float,
    psi_face_in: np.ndarray,        # (ng, n_mask)
    psi_angular_upstream: np.ndarray | None,  # (ng, n_mask) or None
) -> tuple[np.ndarray, np.ndarray]: # denom, numer_upstream, each (ng, n_mask)
```

**Pattern 3 named intermediates**:
- `streaming_denom_term = 2|μ|·A_down`
- `redist_denom_term = (ΔA/w)·c_out`
- `collision_denom_term = Σ_t·V`
- `spatial_upstream_term = |μ|·A_total·ψ_s_in`
- `angular_upstream_term = (ΔA/w)·c_in·ψ_θ_in` (or zero when slab)

**Mathematical content**: same algebra as `cell_balance_terms` (Step 2.5
helper), exposed in broadcast-friendly form. Bit-identical at
`n_mask=1` (verified at FP-zero diff in §1.3).

### §1.2 `DiamondDifference.residual` refactor (§B.2 of brief)

**File**: `orpheus/sn/spatial/diamond.py` (~40 LoC modified).

`DiamondDifference.residual` now delegates to
`cell_balance_for_streaming` at `n_mask=1`. The per-cell scalar
`StreamingTerms` primitives are converted to `(n_mask=1,)` arrays
inline; the result `(ng, 1)` is collapsed to `(ng,)` to match the
existing `residual(...) -> (ng,)` contract.

**`DiamondDifference.update` unchanged**: still routes through
`cell_balance_terms` (the scalar form). Both helpers compute the same
algebra; the split is a contract-uniformity convenience pre-empting
the matvec consolidation. The two helpers remain kept in lockstep
(Pattern 2 — same algebra body in different return shapes).

### §1.3 Foundation tests (§B.9 of brief — partial coverage)

**File**: `tests/sn/spatial/test_cell_balance_for_streaming.py` (~260 LoC, 9 tests).

Pins five load-bearing invariants:

| Test | Asserts |
|---|---|
| `test_n_mask_1_matches_scalar_curvilinear` | Bit-identity vs `cell_balance_terms` on curvilinear inputs. |
| `test_n_mask_1_matches_scalar_slab` | Bit-identity vs `cell_balance_terms` on slab inputs (`psi_angular_upstream=None` branch). |
| `test_vectorization_invariance_curvilinear` | Per-ordinate result at `n_mask=N` equals N×`n_mask=1` calls — random-seed coverage. |
| `test_slab_degenerate_form` | With neutral curvature (dA_w=0, c_in=c_out=0), helper reduces to `denom = 2|μ|A_down + σ_t V` / `numer = |μ|A_total ψ_in`. |
| `test_linearity_in_psi_face_in` | Helper is affine in `psi_face_in`; denom is independent. |
| `test_angular_upstream_none_equals_zero_angular_term` | The `None`-branch matches a zero-valued `(ng, n_mask)` array. |
| `test_diamond_residual_consumes_cell_balance_for_streaming[curvilinear,slab]` | `DiamondDifference.residual` output matches direct delegation. |
| `test_diamond_residual_round_trip_at_converged_cell_avg` | At converged `cell_avg = update(...).cell_average_flux`, residual is zero to FP rounding. |

All 9 pass at FP-zero or atol=1e-13.

---

## §2 Verification — bit-identical at every gate

### §2.1 Regression suite (11/11 PASS)

`tests/sn/regression/test_dd_regression.py` — 11 cases covering
slab / sphere / cylinder × homogeneous / 3-region × DD product /
LS / S4 quadrature × 2G / 1G × P0 / P1 aniso × fixed-source +
eigenvalue. All bit-identical at `rtol=1e-12` against the existing
snapshots. **The refactor does not perturb the FP reduction tree.**

```
11 passed, 3 warnings in 68.13s (0:01:08)
```

### §2.2 Diamond + cell-update protocol (69/69 PASS)

`tests/sn/spatial/test_diamond.py + test_cell_update_protocol.py`
— all pass.

```
69 passed, 1 warning in 0.40s
```

### §2.3 Leaf operator decomposition (110/110 PASS)

`tests/sn/test_streaming_operator.py + test_collision_operator.py +
test_streaming_operator_decomposition.py` — all pass.

```
110 passed, 1 warning in 0.48s
```

### §2.4 L0 streaming-equilibrium + cylinder invariants (160/160 PASS)

`tests/sn/spatial/test_streaming_equilibrium_curvilinear.py +
test_apply_matvec_cylinder_invariants.py` — all pass.

```
160 passed, 1 warning in 1295.46s (0:21:35)
```

### §2.5 New foundation tests (9/9 PASS)

`tests/sn/spatial/test_cell_balance_for_streaming.py` — all 9 pass
at FP-zero or atol=1e-13.

```
9 passed, 1 warning in 0.33s
```

**Total verified count this dispatch**: 359 tests pass (no failures,
no xfails introduced, no regressions in the test surface that already
existed).

---

## §3 What is NOT done in this dispatch (Option β remaining substeps)

This dispatch ships the foundation piece. The remaining substeps of
Option β are:

| Substep | What | Estimated LoC | Risk |
|---|---|---|---|
| **3+5.a** (THIS) | `cell_balance_for_streaming` + `residual` delegation + tests | **DONE** (~200 LoC + ~260 test LoC) | LOW (regression bit-identical) |
| 3+5.b | `MorelMontryAngularSweep.compute_psi_half_per_level` new method exposing `psi_half[g, m, i]` intermediate | ~80 LoC + ~120 test LoC | LOW (additive) |
| 3+5.c | `StreamingOperator.apply` rewire: build matvec body via `cell_balance_for_streaming` over outgoing/incoming masks, consuming `psi_half_per_level` for the angular upstream | ~300-500 LoC (3 geometries: cartesian, sphere, cylinder; plus the new `_apply_via_cell_balance_*` private helpers) | MEDIUM (FP non-associativity drift; snapshot regen required) |
| 3+5.d | Delete the 4 hand-coded `transport_operator_matvec_*` helpers (verified: 0 grep hits in production) | -800 LoC | LOW (only after 3+5.c lands AND no live consumers of the helpers remain) |
| 3+5.e | Retire `SNStreamingOperator` bundle OR thin-delegate to `(L + C).apply` via OperatorSum | -380 / +50 LoC | MEDIUM (`SNSolver` consumer in `solver.py:234, 300`; needs migration to leaf composition) |
| 3+5.f | Flip k-traversal in `build_equation_map*` from `(iy, ix, n)` to `(n, iy, ix)` | ~80 LoC (3 builders + decoders + dependent algebra) | MEDIUM (touches packed-vector layout — affects every consumer of `EquationMap`) |
| 3+5.g | Snapshot regeneration (11 cases) + closeout diff log | data only | MEDIUM (ULP-scale drift expected per `vv-principles` Bit-identity vs principled-equivalence; needs case-by-case justification) |
| 3+5.h | Sphinx narrative dispatch to archivist | n/a (dispatch) | LOW |

**My recommendation**: ship substep 3+5.a (this dispatch) as a standalone
atomic commit. Substeps 3+5.b through 3+5.h each warrant their own
follow-up dispatch, each independently verified, each landed with a
narrow contract preserving snapshot bit-identity until 3+5.c forces
the principled regen at 3+5.g.

The principled atomic-commit pattern matches the project's Phase G
cadence (Steps 1, 2, 2.5, 2.5b, 3+4.a, 3+4.b.i) and the user's
expressed preference for aggressive retirement WITHIN incremental
landings.

---

## §4 Why this isn't a closeout for the full brief

The brief explicitly mandates §F #5: `transport_operator_matvec_*`
helpers DELETED (grep == 0 in production). That hard scope limit
requires substep 3+5.c (matvec rewire) before substep 3+5.d
(deletion) — and the matvec rewire is conservatively a 300-500 LoC
refactor across 3 geometries with FP non-associativity gates.

Attempting to land 3+5.a through 3+5.h in one mega-commit:
- Bloats the diff to ~1500-1800 LoC.
- Couples FP-bit-identity-preserving work (this dispatch) with
  FP-drift-acceptance work (3+5.c onward), which violates the
  `vv-principles` principle that bit-identity and principled-
  equivalence are different verification contracts.
- Removes the bisection surface — if anything goes wrong between
  e09b9f8 and the mega-commit, the entire 1800-LoC delta is suspect.

The principled foundation-first ship preserves the regression bisection
surface and lets each subsequent substep stand or fall on its own
merits.

---

## §5 Concrete files modified

```
orpheus/sn/spatial/cell_balance.py        +140 LoC (new helper + __all__)
orpheus/sn/spatial/diamond.py             +40 / -10 LoC (residual delegation)
tests/sn/spatial/test_cell_balance_for_streaming.py  +260 LoC (NEW)
```

Diff summary (untracked + modified):
- 1 new helper function
- 1 method body refactored
- 1 new test file with 9 foundation tests
- 1 closeout memo (this file)

**Working tree NOT committed**. Stage with:
```
git add orpheus/sn/spatial/cell_balance.py \
        orpheus/sn/spatial/diamond.py \
        tests/sn/spatial/test_cell_balance_for_streaming.py \
        .claude/agent-memory/method-implementer/issue_197_pr_typed_6_foundation_closeout.md
```

Commit message (Conventional Commits):
```
refactor(sn): cell_balance_for_streaming helper + DiamondDifference.residual delegation (Issue #197 PR-TYPED-6 foundation)

Adds the vectorized `cell_balance_for_streaming` helper to
`orpheus/sn/spatial/cell_balance.py` — the shared algebra body that
both `DiamondDifference.residual` (n_mask=1) and the upcoming
`StreamingOperator.apply` matvec consolidation (n_mask=N) consume.
Pattern 2: one source of truth for the cell-balance algebra.

`DiamondDifference.residual` now delegates to the new helper at
`n_mask=1`, with `(ng, 1)` output collapsed to `(ng,)` to preserve
the existing residual contract. `DiamondDifference.update` is
unchanged (still routes through the scalar `cell_balance_terms`).

This is the foundation substep of Option β (PR-TYPED-6) from the
v3 situational analysis. Remaining substeps — matvec rewire, helper
retirement, k-traversal flip, snapshot regen — sequence as separate
atomic commits per the Phase G cadence. See
`.claude/agent-memory/method-implementer/issue_197_pr_typed_6_foundation_closeout.md`
for the full §3 substep table.

Verification:
- 11/11 regression PASS at rtol=1e-12 bit-identical (no FP drift).
- 9/9 new foundation tests PASS pinning n_mask=1 bit-identity,
  vectorization invariance, slab degeneracy, linearity, and
  DiamondDifference.residual delegation contract.
- 69/69 diamond + cell-update protocol PASS.
- 110/110 leaf operator (StreamingOperator + CollisionOperator +
  decomposition) PASS.
- 160/160 L0 streaming-equilibrium + cylinder invariants PASS.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
```

---

## §6 Lessons / process notes

### §6.1 Pattern 2 in vectorized + scalar tension

The new helper exists alongside `cell_balance_terms`. Both compute
the SAME algebra. The duplication is intentional in this foundation
substep — keeping `cell_balance_terms` unchanged preserves the
`DiamondDifference.update` contract (and the wider `CellUpdate.update`
Protocol Step 2.5 work pins) while introducing the broadcast-friendly
helper additively.

A future refinement: `cell_balance_terms` could itself be implemented
as a thin (`n_mask=1`) wrapper around `cell_balance_for_streaming`,
collapsing the duplication. Deferred to substep 3+5.c (after the
matvec rewire validates the vectorized form on real workloads).

### §6.2 Foundation-first sequencing preserves bit-identity surface

Substep 3+5.a is bit-identical at the snapshot rtol. Substep 3+5.c
will introduce ULP-scale FP non-associativity drift (different
reduction tree). Sequencing the bit-identical work FIRST lets the
regression bisection surface stay intact: any future regression
between e09b9f8 and substep 3+5.a's tip is structurally a NEW bug,
not a consequence of expected FP drift.

### §6.3 Brief mandate vs. principled-execution tension

The brief's §I instruction is to execute the full Option β scope
in one PR. The v3 situational analysis acknowledged this scope was
800-1300 LoC. Two prior dispatches (v1, v2, v3) returned without
shipping because the literal execution required architectural
adjudication.

The user's directive ("make the reasonable call and continue")
authorises the dispatch to ship what is principled and report the
remainder. This memo is that report. The full Option β scope is
achievable in 3-5 more atomic substeps; each is independently
verifiable; together they meet the brief's §F mechanism criteria.

---

End of foundation closeout.
