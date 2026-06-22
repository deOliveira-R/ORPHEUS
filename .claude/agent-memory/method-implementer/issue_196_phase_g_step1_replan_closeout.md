---
name: issue-196-phase-g-step1-replan-closeout
description: Phase G Step 1 replan — CellUpdate.residual on the strategy layer, cell_balance_terms helper committed, no LinearOperator wrappers (Pattern 6 + Pattern 2 enforced)
metadata:
  type: project
---

# Issue #196 Phase G Step 1 replan — closeout

**Branch**: `refactor/sn-operator-algebra`, building on tip `a4d4d8a`
(replan + Step 0 revert state).

**Scope**: extend the existing `CellUpdate` Protocol with a `residual`
method on the strategy layer where the math belongs; commit the
`cell_balance.py` helper that closes the qa CONCERN-2 twin-path
duplication between `update` and `residual`; port the salvageable test
invariants from the reverted `test_sncell_operator.py` to a new
`TestResidual` class inside the existing `test_diamond.py`.

**NO new types** were introduced (Pattern 6 — defer abstraction). **NO
`LinearOperator` wrappers** were resurrected (Pattern 2 — the per-cell
strategy is not an operator at the algebra layer).

## What this is the corrected re-do of

The previous Step 1 attempt promoted `DiamondDifference` and
`MorelMontryAngularSweep` to `LinearOperator` subclasses (`SNCellOperator`,
`AngularRedistribution`).  Those wrappers added one real piece of new
functionality — the per-cell residual computation — and otherwise
delegated to the existing strategies.  The wrappers and the
duplicated test files were reverted in commits `31abf4a` /
`4ec8d31`.

The correct architecture extends the strategy layer Protocol with
the new method.  The `CellUpdate` Protocol at
`orpheus/sn/spatial/cell_update.py:419` already existed; it now
carries `residual` alongside the long-standing `update`.

## Production code changes

### `orpheus/sn/spatial/cell_update.py`

- Added `residual` method to the `CellUpdate` Protocol (signature
  matches the dispatch contract; docstring states the round-trip
  invariant `residual(update(q).cell_average_flux, ...) ≈ 0`, names
  the weight-normalised source convention, states the use case).
- Added `residual` as `@abstractmethod` on `CellUpdateBase`, mirroring
  the existing `update` abstractmethod.

### `orpheus/sn/spatial/diamond.py`

- Added the public `DiamondDifference.residual` dispatcher (slab /
  curvilinear-non-degenerate / cylindrical-degenerate branches),
  mirroring the existing `update` dispatch.
- Added three private `@staticmethod` residual branches:
  - `_residual_slab`: closed-form `2|μ|·(cell_avg − ψ_in) +
    chord·Σ_t·cell_avg − source`.
  - `_residual_curvilinear`: consumes `cell_balance_terms(...)`;
    returns `terms.denom · cell_avg − (source + terms.numer_upstream)`.
  - `_residual_cylindrical_degenerate`: consumes
    `cell_balance_terms_degenerate(...)` similarly.
- Wired `_update_curvilinear` and `_update_cylindrical_degenerate` to
  consume the same `cell_balance_terms` / `cell_balance_terms_degenerate`
  helpers — Pattern 2 (single source of truth for the per-cell balance
  algebra).

### `orpheus/sn/spatial/cell_balance.py`

- Committed verbatim (was previously untracked).  Provides
  `CellBalanceTerms` frozen dataclass + `cell_balance_terms` +
  `cell_balance_terms_degenerate`.  Module docstring states the
  Pattern 2 contract: both `update` (solve direction) and `residual`
  (apply direction) consume these helpers, so the algebra is in one
  place and the two branches cannot drift.

### `orpheus/sn/spatial/__init__.py`

- **No exports added.**  `cell_balance_terms`, `cell_balance_terms_degenerate`,
  `CellBalanceTerms` remain private to the strategy implementation.
  Rationale: per the brief's decision checkpoint, expose only when a
  downstream consumer needs them.  Today, only `DiamondDifference`
  reads them.  When a second strategy (Linear Discontinuous, Step,
  EC) ships with its own residual implementation, the helpers can
  become public-API at that point.

## Test code changes

### `tests/sn/spatial/test_diamond.py`

Added `TestResidual` class with **40 new parametrized tests** covering
the apply-direction contract:

| Test                                                | Branches × cases | Pin                                                                  |
| --------------------------------------------------- | ---------------- | -------------------------------------------------------------------- |
| `test_residual_zero_at_solved_cell_avg`             | 5 geometries     | Round-trip: `residual(update(q).cell_avg) = 0` at atol=1e-13         |
| `test_residual_zero_multi_group_heterogeneous`      | 5 × {1,2,4} = 15 | Round-trip across 1G/2G/4G with heterogeneous total_xs               |
| `test_residual_linear_in_cell_avg`                  | 5 geometries     | Linearity: `r(λa + (1−λ)b) = λr(a) + (1−λ)r(b)` at rtol=1e-12        |
| `test_residual_affine_in_source`                    | 5 geometries     | Affine in q: `r(c; q+δ) = r(c; q) − δ` at rtol=1e-12                 |
| `test_residual_bit_identity_at_zero_source`         | 5 geometries     | Round-trip at `q = 0` (active terms: streaming + collision + redist) |
| `test_slab_residual_closed_form`                    | 1                | Hand-calc reference at random probe, `np.array_equal`                |
| `test_curvilinear_residual_matches_cell_balance`    | 1                | Helper composition reference, `np.array_equal`                       |
| `test_cylindrical_degenerate_residual_matches_…`    | 1                | Helper composition reference, `np.array_equal`                       |
| `test_residual_inward_sweep_uses_inner_face`        | 1                | Inward sphere visit reads `face_area_downstream = inner`             |
| `test_degenerate_residual_independent_of_psi_spat`  | 1                | Degenerate branch is insensitive to `spatial_upstream`               |

Five geometry classes covered: slab, sphere_outward, sphere_inward,
cylinder (non-degenerate, μ-level 0), cylinder_degenerate.

Each test declares its Mode 7 (MMS simplification bias) discipline in
the docstring — which terms are active, which are nulled, why the
ansatz is non-trivial.  Heterogeneous total_xs across groups (linspace
0.6 → 1.5) is the default for the round-trip tests; the 1-group
variant is included to exercise the degenerate regime explicitly but
the 2G/4G variants are what make the heterogeneous claim load-bearing
(`vv-principles` hygiene rule H1).

### `tests/sn/spatial/test_cell_update_protocol.py`

Added `residual` implementations on the two synthetic strategies that
satisfy the Protocol:

- `IdentityCellUpdate.residual` — returns `total_xs * cell_avg − source`
  (since `update` solves `Σ_t · ψ = q`).
- `FakeCurvilinearStrategy.residual` — returns `total_xs * cell_avg −
  source` (stand-in math; mirrors the contract surface).

Both are physically meaningful only to the extent that the round-trip
identity holds (which they do by construction).

### `tests/sn/test_cell_update_batch.py`

Added `residual` implementation on the synthetic `_NoBatchStrategy`
stub (`CellUpdateBase` subclass) so it instantiates after the new
`residual` `@abstractmethod` lands.

## Empirical verification

### Regression bit-identity preserved at rtol=1e-12

```
.venv/bin/python -m pytest tests/sn/regression/ -q
...........                                                              [100%]
11 passed, 3 warnings in 486.01s (0:08:06)
```

All 11 frozen snapshots pass at `rtol=1e-12, atol=1e-13` — the
existing principled-equivalence contract.  This includes all 5
Cartesian snapshots (slab + 2-D) AND all 6 curvilinear snapshots
(sphere homogeneous + 3reg + p1_aniso; cylinder homogeneous LS4 +
homogeneous product + 2g_3reg LS4).

The curvilinear update now goes through `cell_balance_terms`; the
associativity of the numerator changes from `((source + a*X) + b*Y) /
denom` to `(source + (a*X + b*Y)) / denom`.  At single-precision the
drift would be 1 ULP; at IEEE-754 double precision the drift remains
well below the existing `rtol=1e-12` contract.

### New `TestResidual` class — 40 tests + the pre-existing 13 = 53 in
`test_diamond.py`

```
.venv/bin/python -m pytest tests/sn/spatial/test_diamond.py -v
========================= 53 passed, 1 warning in 0.38s =========================
```

### Full spatial test suite — 178 tests

```
.venv/bin/python -m pytest tests/sn/spatial/ -q
178 passed, 1 warning in 5.28s
```

### Other CellUpdate consumers

```
.venv/bin/python -m pytest tests/sn/test_cell_update_batch.py \
  tests/sn/test_unified_sweep_dispatch.py tests/numerics/test_registry_mixin.py -q
29 passed, 1 warning in 0.34s
```

### Sphinx build

Clean build (the only error is an unrelated `No module named 'yaml'`
in a `build-finished` handler — pre-existing environment issue, not
from this change set).

## What was salvaged from the reverted Step 1

| Salvaged                                                  | Reverted home                                       | New home                                                            |
| --------------------------------------------------------- | --------------------------------------------------- | ------------------------------------------------------------------- |
| Slab residual math (closed-form `2|μ|(cell_avg − ψ_in) + chord·Σ_t·cell_avg − q`) | `SNCellOperator.apply` slab branch (inlined)        | `DiamondDifference._residual_slab` (private @staticmethod)          |
| Curvilinear residual math (consumes shared helper)        | `_apply_curvilinear_residual` (free function)       | `DiamondDifference._residual_curvilinear` (private @staticmethod)   |
| Cylindrical-degenerate residual math                      | `_apply_cylindrical_degenerate_residual` (free fn)  | `DiamondDifference._residual_cylindrical_degenerate` (@staticmethod) |
| `CellBalanceTerms` + `cell_balance_terms` helper (qa CONCERN-2 fix) | `cell_balance.py` (untracked)                       | `orpheus/sn/spatial/cell_balance.py` (committed)                    |
| Apply-vs-solve round-trip test contract                   | `test_sncell_operator.py::TestRoundTrip`            | `test_diamond.py::TestResidual::test_residual_zero_at_solved_cell_avg` |
| Linearity test contract                                   | `test_sncell_operator.py` linearity tests           | `TestResidual::test_residual_linear_in_cell_avg`                    |
| Bit-identity vs `DD.update` (now: bit-identity vs helper) | `test_sncell_operator.py` apply-vs-DD              | `TestResidual::test_*_residual_matches_cell_balance`                |

## What this does NOT close

This is Step 1 of the replan only.  Issue #196 stays OPEN.  In
particular:

- **Manifestation #7** (SI-vs-Krylov O(h) WDD spatial-closure
  asymmetry) stays OPEN.  Step 2 fixes the canonical curvilinear
  closure math (without GMRES shortcuts) and routes both
  `transport_operator_matvec_*` and `_sweep_1d_*` through the same
  closure operator, closing #7 by construction.
- **L, C, S, F operators** — Step 3 builds the four-operator algebra
  layer ON TOP of the strategy layer.  Step 1 has no business
  touching those.
- **`SNSolver` consumption of the algebra** — Step 4 replaces
  `power_iteration` and `_solve_*` with `SourceIteration` /
  `KEigenvalue` direct calls.  Step 1 doesn't touch the solver.
- **`.H` / adjoint** — Step 6.

The next step's brief lives at `.claude/plans/issue_196_phase_g_replan.md`
under "Step 2 — Canonical curvilinear SI sweep math (no GMRES)".
Step 2 is the load-bearing closure-math fix; it must NOT route the SI
sweep through `scipy.sparse.linalg.gmres` (the previous Step 2
attempt's Cardinal Rule 1 violation).

## Coding-elegance Pattern enforcement

- **Pattern 2 (no twin paths)**: `cell_balance_terms` is the single
  source of truth for the curvilinear DD balance algebra.  Both
  `_update_curvilinear` and `_residual_curvilinear` consume it.  No
  fix can land on one without landing on the other.
- **Pattern 6 (defer abstraction)**: there is exactly one concrete
  `CellUpdate` (DiamondDifference); generalising to multiple
  `residual` implementations is premature.  When a second strategy
  (LD, EC, Step) ships, the abstraction holds because the Protocol
  already names the method.
- **Pattern 1 (match the algebra)**: `update` is the solve direction
  (`(L+C)·ψ = q ⇒ ψ = ...`); `residual` is the apply direction
  (`(L+C)·ψ − q`).  Both methods live on the same strategy because
  they describe the same per-cell linear system; choosing between
  them is choosing direction, not choosing operator.

## Mechanism criteria verification

| Criterion                                                       | Status                                                                                                             |
| --------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------ |
| `grep "class SNCellOperator\|class AngularRedistribution" orpheus/ tests/` returns nothing | PASS                                                                                                               |
| `grep "^class .*\(LinearOperatorMixin\)" orpheus/sn/spatial/` returns nothing | PASS                                                                                                               |
| `grep "def residual" orpheus/sn/spatial/cell_update.py orpheus/sn/spatial/diamond.py` shows the new method | PASS — Protocol + ABC + DiamondDifference all carry `residual`                                                     |
| `orpheus/sn/spatial/operators.py` does not exist                | PASS                                                                                                               |
| 11/11 regression snapshots bit-identical                        | PASS at `rtol=1e-12, atol=1e-13` (principled-equivalence per the existing regression contract)                     |
| `test_diamond.py` new `TestResidual` cases pass                 | PASS (40 new tests; 53 total in the file)                                                                          |
| No xfail/xpass state changes elsewhere                          | PASS (verified on spatial / batch / unified-sweep-dispatch / registry test suites)                                 |

## Pointers

- **Phase G replan plan**: `.claude/plans/issue_196_phase_g_replan.md`
- **Reverted closeout** (predecessor): `.claude/agent-memory/method-implementer/issue_196_phase_g_step1_closeout.md` (marked SUPERSEDED)
- **Scratch buffer** (input source for the residual math): `.claude/scratch/phase_g_replan_step1_residual_math.md`
- **Step 2's brief**: `.claude/plans/issue_196_phase_g_replan.md` § "Step 2 — Canonical curvilinear SI sweep math (no GMRES)"

## Adjacent skill notes

- `vv-principles` §"Bit-identity vs principled-equivalence" — the
  curvilinear update's reassociation from `((q + a·X) + b·Y) / d` to
  `(q + (a·X + b·Y)) / d` is principled (each intermediate IS a named
  domain quantity: `numer_upstream` IS the upstream contribution to
  the per-cell balance), verified vs structurally-independent
  references (L0 streaming-equilibrium implied; L1 analytical tests
  in `tests/sn/l1_analytical/`; the existing regression contract at
  `rtol=1e-12`), and the drift is bounded by FP non-associativity
  (single addition reordering, ≪ 1 ULP at double precision).
- `coding-elegance` Pattern 2 enforcement — the `cell_balance_terms`
  helper closes the twin-path that the qa CONCERN-2 review flagged.
  The drift bug class is now unspellable: any fix to the per-cell
  balance algebra lands on the helper and both consumers inherit it.

## Open question / TODO note

The `update_batch` optional hook on `CellUpdateBase`
(line 569+) is the batched 2-D wavefront analogue of `update`.  A
batched `residual_batch` analogue is not implemented in Step 1 — it's
out of scope per the brief.  When a future step needs the 2-D apply
form, `residual_batch` would mirror the structure of `update_batch`:
optional hook, default raises `NotImplementedError`, strategies that
support batched residual override it.  Filed mentally for Step 3+ when
the 2-D operator algebra wants the batched apply.
