---
name: issue-196-phase-g-step1-closeout
description: Phase G Step 1 (Issue #196) closeout — type-system promotion of DiamondDifference → SNCellOperator(LinearOperator) and MorelMontryAngularSweep → AngularRedistribution(LinearOperator). Wrapper-only; bit-identical regression snapshots; sets up Step 2's call-site unification.
metadata:
  type: project
  branch: refactor/sn-operator-algebra
  step: Phase G Step 1
  date: 2026-05-12
---

# Phase G Step 1 — Closeout

## What shipped

**Production code**

- NEW `orpheus/sn/spatial/operators.py` (~440 lines) carrying two
  LinearOperator promotions:
  - `SNCellOperator(LinearOperatorMixin)` — wraps `DiamondDifference`.
    `solve(visit, total_xs, source, upstream_state) → CellResult`
    delegates to `DiamondDifference.update`.  `apply(cell_avg, *,
    visit, total_xs, upstream_state, source) → residual` computes
    the per-cell discretised operator residual `L_cell · ψ̄ − q` via
    per-branch helpers `_apply_curvilinear_residual` and
    `_apply_cylindrical_degenerate_residual` (the slab branch's
    residual lives inline because slab has no curvilinear bundle).
    Capabilities: `{CAP_APPLY, CAP_SOLVE}`.
  - `AngularRedistribution(LinearOperatorMixin)` — wraps
    `MorelMontryAngularSweep`.  `apply(psi_cells, *, alpha_half,
    redist_dAw, tau_mm, volume, level_indices=None, carlson_context=None) → R`
    delegates to `MorelMontryAngularSweep.__call__`.
    Capabilities: `{CAP_APPLY}` only (per Pattern 4 decision in
    open-question Q4 below).

- `orpheus/sn/spatial/__init__.py` updated to export both new
  classes.

**Tests**

- NEW `tests/sn/spatial/test_sncell_operator.py` — 83 tests filled
  in.  All pass.  Covers Gates 1 (bit-identity slab + sphere +
  cylinder + cyl-degenerate), 2 (apply-solve round-trip at
  `rtol=1e-12`), 3 (capability surface), 4 (geometry coverage as
  parametrize on Gates 1-3).

- NEW `tests/sn/spatial/test_angular_redistribution.py` — 85 tests
  filled in + 4 skipped.  Covers Gate 3 (capabilities), Gate 6 #3
  (Carlson seed equivalence — 54 tests), Gate 6 #4 (flat-flux
  closure — 16 tests), Gate 6 #5 (linearity in input — 12 tests).
  The 4 skipped tests are the round-trip class deferred until
  `CAP_SOLVE` ships in a later step.

- `tests/sn/spatial/test_sweep_vs_apply_consistency.py` extended —
  Phase F twin-path defense test body filled in (currently
  xfail-strict=False under Step 1 because the call sites haven't
  been unified yet; transitions to xpass at Step 2).  Step 2
  acceptance gate `test_sncell_operator_consumed_by_both_call_sites_at_step_2`
  marked skip with rationale until Step 2.

## Empirical evidence

| Gate | Pre-Step-1 | Post-Step-1 |
|---|---|---|
| Gate 1 — bit-identity slab | n/a (test new) | **83 pass** in `test_sncell_operator.py` |
| Gate 2 — apply-solve round-trip | n/a | residual ~7e-18 (slab); 0.0 (sphere); 0.0 (cyl-degen); passes `rtol=1e-12` on all 72 parametrize combos |
| Gate 3 — capability surface | n/a | 4 pass (SNCellOp) + 3 pass (AngularRedist) |
| Gate 4 — geometry coverage | n/a | covered as parametrize of Gates 1-3 |
| Gate 5 — Phase F twin-path defense | n/a | **xfail-strict=False** at Step 1 (expected — call-site unification is Step 2) |
| Gate 6 — M-M algebraic identities | n/a | 82 pass (54 Carlson equivalence + 16 flat-flux + 12 linearity) |
| Gate 7 — regression snapshots | 11 pass (baseline) | **11 pass** at `rtol=1e-12` — BIT-IDENTICAL |

**Full sn/spatial/ suite**: 306 pass, 5 skipped, 1 xfailed.

**Regression confirmation**: 11/11 snapshots stay bit-identical
(`np.array_equal` regression contract preserved — Step 1 is pure
plumbing, no float-reduction-tree change, no algorithmic change).

## Decisions on the 5 open questions

| # | Test-architect's recommendation | Method-implementer's decision | Rationale |
|---|---|---|---|
| Q1 — apply signature | Bound-per-cell pattern preferred for elegance | **Keyword-args pattern** (`apply(cell_avg, *, visit, total_xs, upstream_state, source)`) | The test skeleton's Gate 2 contract (line 422-426 of the original skeleton) explicitly used the keyword-args pattern.  Per the dispatch instructions ("read the test skeletons first to determine which pattern they encode, then implement to match"), I followed the test skeleton's contract.  The bound-per-cell pattern is a future refactor candidate but doesn't fit Step 1's test-contract-as-spec discipline. |
| Q2 — solve returns | Full `CellResult` | **Full `CellResult`** (as recommended) | Bit-identity gate compares every CellResult field via `np.array_equal`; returning the full result keeps the bit-identity contract field-by-field rather than collapsing to just `cell_average_flux`. |
| Q3 — apply input shape | `(ng, nx)` φ_0 (the M-M recurrence's natural input; seed folds ψ → φ_0 already) | **`(ng, M, nx)` psi_cells** | The M-M recurrence's actual input is `psi_cells` of shape `(ng, M, nx)` — the per-ordinate ψ on this level.  The seed folds ψ → φ_0 internally, but the recurrence ITSELF still consumes `(ng, M, nx)` at the slice `psi_level[:, m, :]` per ordinate.  Passing `(ng, nx)` would require an artificial broadcast inside the operator that doesn't reflect what the wrapped primitive consumes.  Pattern 5 (build the right primitive): the operator's `apply` shape MUST match the wrapped primitive's natural input shape. |
| Q4 — AngularRedist capabilities | Drop `CAP_SOLVE` at Step 1 | **Drop `CAP_SOLVE`** (as recommended) | Per Pattern 4 (advertise only what works).  No back-substitution code ships at Step 1; advertising `CAP_SOLVE` without the implementation would be the harmful-stub anti-pattern this module is designed against. |
| Q5 — MissingCapability trigger | Use `compose_requiring_adjoint(op)` or assert `not in caps` | **`op.H.apply(np.array([1.0]))` raises** | The Wave-0 `_AdjointOperator.apply` checks for `CAP_APPLY_TRANSPOSE` on the inner operator and raises `MissingCapability` at call time when missing.  This is the composer-level trigger that exercises Pattern 4's runtime enforcement.  Also asserted at surface level via `assert CAP_APPLY_TRANSPOSE not in op.capabilities`. |

## Divergences from the test skeleton (intentional)

1. **Round-trip assertion target**.  The Gate 2 test docstring said
   `assert_allclose(residual, q)`.  But mathematically the residual
   `L_cell · ψ̄ − q` at the solved value `ψ̄ = solve(q)` is `0`, not
   `q`.  I wrote the assertion against zero
   (`assert_allclose(residual, np.zeros_like(residual))`), which is
   the meaningful contract.  The docstring's `≡ q` wording reflects
   an older interpretation where `apply` is the action `L·x` (not the
   residual); under that interpretation the round-trip would be
   `apply(solve(q)) ≡ q`.  My choice (apply = residual) is more
   consistent with the per-cell discretised operator's algebra:
   `apply` returns "how far am I from solving the cell balance",
   which is zero at the solved value.

2. **`AngularRedistribution` input shape `(ng, M, nx)`** (Q3 above)
   instead of `(ng, nx)`.  Rationale documented above.

3. **`test_round_trip_at_q_zero_trivial`** uses `np.testing.assert_array_equal`
   rather than `assert_allclose(rtol=1e-12, atol=1e-13)`.  At
   `q = 0` and zero upstream, the residual is exactly `0` by
   construction (no FP rounding from non-trivial combinations).
   The exact `array_equal` contract is tighter and catches a wider
   class of regressions.

## Files touched

- **Created**: `orpheus/sn/spatial/operators.py` (~440 lines)
- **Updated**: `orpheus/sn/spatial/__init__.py` (added 2 exports)
- **Created**: `tests/sn/spatial/test_sncell_operator.py` (~500 lines, 83 tests + 0 skip)
- **Created**: `tests/sn/spatial/test_angular_redistribution.py` (~370 lines, 85 tests + 4 skip)
- **Updated**: `tests/sn/spatial/test_sweep_vs_apply_consistency.py` (Phase F twin-path defense body filled, Step 2 gate marked skip)

## What this does NOT close

- **ERR-026 manifestation #7** (SI-vs-Krylov O(h) WDD asymmetry on
  heterogeneous MR) — Step 1 is type-system promotion only.  The
  call sites still live in twin functions
  (`transport_operator_matvec_*` and `_sweep_1d_*`).
  Manifestation #7 closes at Step 2 by construction when both call
  sites route through the SAME `SNCellOperator` instance with the
  SAME closure config.
- **Issue #196** stays OPEN until Step 2 lands.
- **Phase E flux-shape sentinel xfail** — same trigger as
  manifestation #7; transitions to xpass at Step 2.

## Type-checking

`pyright` and `mypy` are not installed in the project's virtualenv;
the type-check step was skipped.  The new module uses standard type
annotations (`np.ndarray`, `LinearOperatorMixin` inheritance,
`frozen=True, slots=True` dataclasses) that should be straightforward
for any type-checker; if pyright is later installed, a clean run is
expected.

## Sphinx narrative

No new equation labels at Step 1.  The existing labels
(`dd-curvilinear-scalar`, `dd-mm-closure-constants`,
`phase-f-carlson-seed-source-driven`, `phase-f-q-bar-twin-forms`,
`mm-weights`, `pole-mm-recurrence`) cover the Step 1 claims via
`@pytest.mark.verifies` decorators on the new tests.  Bulk narrative
for Phase G's four-operator unification (touching all Steps 1-5)
will be picked up by the parallel archivist dispatch — Step 1 by
itself doesn't earn a narrative section.

## Next steps

1. **Step 2** — numerics-investigator pre-dispatch to characterise
   the per-cell SI-vs-Krylov drift on `sphere_2g_3reg` n=40 at the
   converged fixed point.  Then:
   - Construct `StreamingOperator(L)` wiring at
     `orpheus/sn/operator.py` or a new `orpheus/sn/streaming.py`.
   - Both `transport_operator_matvec_*` and `_sweep_1d_*` route
     through ONE `SNCellOperator` instance with the SAME closure
     config.
   - Phase E flux-shape sentinel xfail → xpass.
   - Phase F twin-path defense test xfail → xpass (this Step 1
     deliverable's xfail marker self-removes).
   - Step 2 acceptance gate
     `test_sncell_operator_consumed_by_both_call_sites_at_step_2`
     skip → pass.
2. **Steps 3-5** continue per the Phase G plan.

## Commit pair

Two atomic commits per Step 1's deliverable spec:

1. `feat(sn): SNCellOperator + AngularRedistribution LinearOperator promotion (Issue #196 Step 1)` — production code: `orpheus/sn/spatial/operators.py` + `__init__.py`.
2. `test(sn): apply-vs-solve invariants for cell + angular ops (Issue #196 Step 1)` — test files: 3 test files (2 new + 1 extended).

## Lessons / observations for future steps

1. **Residual ≠ action** for per-cell operators.  The "apply"
   method on a discretised per-cell operator naturally returns the
   residual `L·ψ̄ − q` (because q is a per-cell input, not a
   field-level concept).  The forward action `L·ψ̄` standalone
   doesn't have a clean per-cell interpretation because the
   "what flows in" is the same input that `solve` consumes.  This
   may motivate a different signature for the L (streaming
   operator) at Step 2 where the source is field-level.

2. **`MorelMontryAngularSweep` natural input is `(ng, M, nx)`** —
   the per-level slice.  The Carlson seed strategy on it folds
   ψ → φ_0 internally, but the recurrence itself reads per-ordinate
   slices.  The Step 2 unification will need to decide whether L's
   public `apply` takes `(ng, M, nx)` or the global `(N, nx, ny, ng)`
   shape from `solve_sn`.

3. **The bound-per-cell pattern** for `SNCellOperator` was the more
   elegant choice per Pattern 5 (visit/total_xs/upstream bound at
   construction; apply takes only cell_avg).  But the test skeleton
   pinned the keyword-args pattern.  Worth revisiting at Step 2 if
   per-iteration allocation cost surfaces as a concern.

## Pointers

- Plan: `.claude/plans/issue_196_phase_g_four_operator_unification.md`
- Reconciliation: `.claude/agent-memory/general-purpose/phase_g_four_operator_architecture_reconciliation.md`
- Step 1 gate design: `.claude/agent-memory/test-architect/issue_196_phase_g_step1_verification_gates.md`
- Phase F closeout: `.claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md`
- Production: `orpheus/sn/spatial/operators.py`
- Tests: `tests/sn/spatial/test_sncell_operator.py`, `tests/sn/spatial/test_angular_redistribution.py`
- Twin-path defense: `tests/sn/spatial/test_sweep_vs_apply_consistency.py`

## Linked memories

- `[[issue-168-phase-f-closeout]]` — the failure-mode profile this
  Step defends against.
- `[[issue-168-phase-d-closeout]]` — the original Carlson seed
  introduction that AngularRedistribution wraps.
