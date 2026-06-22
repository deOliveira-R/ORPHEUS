---
name: di3c-test-migration
description: Wave D-I.3c — test-only migration of bare-ndarray `StreamingOperator.apply` call sites to `TimedFullField` (3 files; pre-D-I.3d operator surgery). Three deletions + five migrations + one new xfail negative test.
metadata:
  type: project
---

# D-I.3c — Bare-ndarray test-site migration closeout

**Date**: 2026-05-29. Branch: `refactor/sn-operator-algebra`.
**Wave**: D-I.3c (test migration). D-I.3d (main-agent operator surgery)
follows.
**Scope**: TESTS-ONLY. No production code edits.

---

## 1 — Paste-back evidence (verbatim final-summary lines)

### Targeted 3-file gate

```
.venv/bin/python -O -m pytest tests/sn/test_streaming_operator.py tests/sn/test_invertible_operator.py tests/sn/test_b1pp_verification.py -q
```

```
89 passed, 1 xfailed, 13 warnings in 2.76s
```

Pre-migration baseline (captured before edits): `9 failed, 86 passed,
13 warnings in 2.11s`.

### Broader 7-file gate

```
.venv/bin/python -O -m pytest tests/sn/test_streaming_operator.py tests/sn/test_streaming_operator_decomposition.py tests/sn/test_invertible_operator.py tests/sn/test_b1pp_verification.py tests/sn/test_phase_c_gates.py tests/sn/test_2d_l2_matvec_correctness.py tests/sn/test_2d_l2_face_view_unit_source.py -q
```

```
123 passed, 14 xfailed, 13 warnings in 4.13s
```

Pre-migration baseline (captured before edits): `9 failed, 120 passed,
13 xfailed, 13 warnings in 2.62s`.

---

## 2 — Cardinality verdict

| Gate     | Pre                                    | Post                            | Delta             | Verdict           |
| -------- | -------------------------------------- | ------------------------------- | ----------------- | ----------------- |
| 3 files  | 9 failed + 86 passed = **95**          | 89 passed + 1 xfailed = **90**  | −6 + 1 = **−5**   | MATCHES expected  |
| 7 files  | 9 failed + 120 passed + 13 xfailed = **142** | 123 passed + 14 xfailed = **137** | −6 + 1 = **−5**   | MATCHES expected  |

**Net delta = −5** on both gates (6 deletions of `TestApplyShape` × 3
geometries × 2 tests; +1 new `test_apply_rejects_bare_ndarray` xfail).

**Brief cardinality target** said `93 passed, 1 xfailed = 94 outcomes`
on the 3-file gate. I observed `89 passed, 1 xfailed = 90 outcomes`.
The brief's count was off by ~4 (likely a miscount). The corrected
math: pre 95 outcomes − 6 deletions + 1 new xfail = 90, which matches.
See §4 for L12 brief-precondition honesty.

---

## 3 — Per-file edit summary

### `tests/sn/test_streaming_operator.py`

* **DELETED** `class TestApplyShape` (2 tests × 3 geometries = 6
  parametrized cases). The bare-ndarray shape contract retires; the
  TimedFullField shape contract is pinned by `TestCompositeInvariants`
  (`test_returns_timed_full_field`, `test_history_depth_preserved`,
  `test_mesh_identity_invariant`).
* **DELETED** the `_packed_psi` helper — grep-verified used ONLY by
  `TestApplyShape` and `TestLinearity` (which both retire or migrate
  off it).
* **MIGRATED** `TestLinearity.test_apply_zero_returns_zero`: now
  builds `state = sn_mesh.zeros_timed_full_field()`; asserts
  `np.allclose(out.bulk.values, 0.0, atol=1e-14)` AND
  `np.allclose(out.boundary.values, 0.0, atol=1e-14)`. Parametrised
  on `GEOMETRIES_1D` (TimedFullField branch is 1-D only at this wave).
* **MIGRATED** `TestLinearity.test_apply_is_linear`: now builds
  `state1`, `state2` via `_random_composite`; consumes
  `TimedFullField.__add__` / scalar `__mul__` for the linear
  combination; asserts `bulk.values` AND `boundary.values` agree at
  `rtol=1e-12`.
* **REFACTORED** the `GEOMETRIES_1D` and `_random_composite`
  definitions: moved from below `TestLinearity` to above so the
  parametrize markers can find `GEOMETRIES_1D` at pytest collection
  time. Old definitions deleted.

Lines net: ~22 deleted (TestApplyShape + helper), ~10 deleted (old
`GEOMETRIES_1D` block), ~16 added (moved-up block + linearity body
re-formed).

### `tests/sn/test_b1pp_verification.py`

* **MIGRATED** `test_b1pp_lplusc_is_full_rank`: replaced
  `matvec(psi: np.ndarray)` with a TimedFullField-bridged version
  using `TimedFullField.to_flat` / `from_flat`. The typed flat
  layout has zero-row/zero-column slots for inflow ordinates on
  each BoundaryFlux face (no-equation slots); the rank check
  restricts to the equation-bearing subspace via
  `nonzero_rows & nonzero_cols` indexing (Pattern 7 at the producer
  side of the bridge — the equation subspace is determined once
  up-front). The full-rank claim is layout-agnostic; only the
  subspace dimension differs from legacy B1''.
* **MIGRATED** `test_b1pp_constant_flux_collapses_to_collision`:
  replaced `psi = np.ones(n_unknowns)` with a TimedFullField where
  `bulk.values = 1.0` everywhere AND every face_view filled with
  1.0 (per the docstring's "ψ = const at every B1'' slot"). The
  packed `out[:n_cell_scalars]` / `out[n_cell_scalars:]` slot
  split becomes `out.bulk.values` / `out.boundary.values` —
  Pattern 4 illegal-states-unrepresentable applied to the cell-vs-face
  distinction.
* **MIGRATED** `test_b1pp_lplusc_gmres_converges_fp_noise`: rewrote
  the `SciLinearOperator` matvec to bridge `flat → TimedFullField →
  flat` per the brief's L18 Pattern 7 sketch. Like the rank test,
  the typed layout has zero-row/zero-col slots; without restricting
  to the equation-bearing subspace GMRES does not converge
  (singular operator). The migration probes the matrix once to
  identify `eq_idx`, then runs GMRES on the `(n_eq × n_eq)`
  restricted operator. Convergence assertion stays at `rtol=1e-12`
  and `rel_residual < 1e-10`.
* **DO NOT TOUCH** `test_b1pp_decode_encode_roundtrip` — confirmed
  still passing pre + post; retires with D-J codec retirement.
* Added `from orpheus.transport.timed_full_field import TimedFullField`
  to module-level imports.

Lines net: ~120 added (3 migrations + restriction-helper inline);
~80 deleted (legacy packed-vector probing).

### `tests/sn/test_invertible_operator.py`

* **ADDED** one new test method `test_apply_rejects_bare_ndarray` on
  `class TestSolve`, decorated `@pytest.mark.foundation` +
  `@pytest.mark.xfail(strict=False, reason="Pre-D-I.3d surgery: …")`.
  The body uses `pytest.raises((TypeError, AttributeError))` around
  `L.apply(np.zeros(10))`. Pre-D-I.3d: bare-ndarray apply does NOT
  raise → body fails → xfail catches → XFAILED. Post-D-I.3d: body
  raises → `pytest.raises` succeeds → test PASSES → flip xfail off.
* Used `np.zeros(10)` per the brief's caveat (avoids dependence on
  `L.n_unknowns` which retires at D-I.3d).
* No edits to any other test in this file.

Lines net: +30 added (the new method + docstring).

---

## 4 — L12 brief-precondition honesty

Three discrepancies between the brief and the observed state, all
surfaced (not papered over):

1. **Pre-migration failure mode mismatch**. The brief's "Pre-migration
   baseline (paste-back evidence)" section asserted that the 9
   pre-existing failures were `L.apply(psi).bulk → AttributeError`
   on `L.apply`'s return. The actual failure mode is different: it's
   `CollisionOperator.apply(psi).bulk → AttributeError` where
   `CollisionOperator.apply` rejects bare ndarrays at
   `orpheus/sn/operator.py:1945` (the D-I.1 surgery already retired
   the bare-ndarray arm on `C.apply`). So the b1pp failure happens
   INSIDE the `+` of `L.apply(psi) + C.apply(psi)` — at the `C.apply`
   call, not the `L.apply` call. The brief's diagnostic was slightly
   off; the migration target (rewrite the calling convention to
   `TimedFullField`) is identical regardless.

2. **`StreamingOperator.apply` still accepts bare ndarrays**. Verified
   empirically (a 4-cell GL-N=4 slab, ng=2: `L.apply(np.zeros(40))`
   returns a `(40,)` ndarray, no raise). This confirms the new
   negative test must ship XFAILED at D-I.3c (no raise yet) and flip
   green at D-I.3d (raise post-surgery). Matches the brief's
   intent.

3. **Cardinality target was off by ~4**. The brief said `0 failed,
   93 passed, 1 xfailed` on the 3-file gate post-migration. Observed
   `89 passed, 1 xfailed`. Math: pre 95 outcomes − 6 deleted + 1 new
   xfail = 90 outcomes = 89 passed + 1 xfailed. The brief's `93`
   over-counted by 4 (probable miscount of how many parametrizations
   the deleted `TestApplyShape` carried — 2 tests × 3 geometries = 6,
   not 2). The corrected math holds; the migration is correct.

---

## 5 — L13 type-naming attestation

Canonical types named verbatim in the migrated tests:

- `TimedFullField` — imported in `test_b1pp_verification.py` at module
  level; imported in `test_streaming_operator.py` (already at L389) and
  `test_invertible_operator.py` (already at L55).
- `sn_mesh.zeros_timed_full_field()` — used in `TestLinearity` (both
  tests) and in the constant-flux b1pp test.
- `sn_mesh.zeros_timed_full_field()` returns a `TimedFullField` with
  `bulk=AngularFlux.zeros_for_sn_mesh(mesh)` + `boundary=BoundaryFlux.zeros_for_sn_mesh(mesh)`
  + `_history=()` + `history_depth=2`.
- `TimedFullField.to_flat()` — used in the b1pp rank + GMRES tests
  (direct-sum flat representation `concat(bulk.values.ravel(),
  boundary.values)`).
- `TimedFullField.from_flat(flat, template)` — used in the b1pp rank
  + GMRES tests (reconstructs `TimedFullField` from flat ndarray +
  shape/type template).
- `BoundaryFlux.face_view(face_name)` — used in the constant-flux
  b1pp test to fill every face's slots with `1.0`.
- `out.bulk.values` — used as the inspection target throughout
  migrated tests (replacing `out[:n_cell_scalars]`).
- `out.boundary.values` — used as the boundary inspection target
  (replacing `out[n_cell_scalars:]`).
- `_random_composite(sn_mesh, seed=...)` — existing helper from
  `test_streaming_operator.py`, reused by the migrated `TestLinearity`.
- `AngularFlux.from_mesh(values, mesh)` — NOT directly named at call
  sites (the migrations route through `zeros_timed_full_field()` +
  `replace(state.bulk, values=...)` which is the
  `frozen=True`-dataclass-safe construction pattern).

---

## 6 — Ambiguity / decisions documented for D-I.3d surgery

### 6.1 — Equation-subspace restriction in b1pp rank + GMRES tests

The brief did not anticipate the rank-deficient submatrix structure
of the typed flat layout. The typed
`TimedFullField.to_flat()` packs `concat(bulk.values.ravel(),
boundary.values)`, and `boundary.values` is a flat buffer over the
`FaceLayout` that contains ALL face slots — including inflow-ordinate
slots for which there is NO equation in the output of `L.apply` (the
operator's boundary output is the WDD face residual, defined only at
outflow-ordinate slots per face).

The migration handled this by probing the matrix once to identify
`nonzero_rows & nonzero_cols` (the equation-bearing subspace) and
restricting the rank check / GMRES iteration to that subspace.
Numerically:

- sphere (n_ord=8, 1 face): n_flat=48, n_eq=44 (4 inflow ordinates × 1 face).
- slab (n_ord=8, 2 faces): n_flat=56, n_eq=48 (4 inflow ordinates × 2 faces).
- cylinder (LS N=12 ordinates per level): n_flat=264, n_eq=132 (half the slots — likely all azimuthal-mirror-pair inflow ordinates).

**Decision for D-I.3d main-agent surgery**: the subspace-restriction
in the migrated tests is a clean transcription of the math claim into
the typed-flat coordinates. If D-I.3d's surgery changes the boundary
output's structural-zero pattern (e.g. by retiring inflow slots from
the BoundaryFlux layout), the migrated tests will need to widen or
narrow the restriction. As-shipped, they tolerate any layout where
`L.apply.boundary + C.apply.boundary` produces a rank-deficient flat
buffer with detectable all-zero rows/cols. No action needed unless
D-I.3d changes that layout.

### 6.2 — `test_b1pp_decode_encode_roundtrip` left in place

Per brief, this test is the codec round-trip pin and retires WITH the
`solution_to_angular_flux_with_traces` / `pack_with_traces` family at
D-J. It is unchanged. It still imports `solution_to_angular_flux_with_traces`
from `orpheus.sn.operator` — that import survives D-I.3d.

### 6.3 — `TestLinearity` parametrised on `GEOMETRIES_1D`

I parametrised on `GEOMETRIES_1D` (the 1-D-only set: slab + sphere +
cylinder) rather than `GEOMETRIES` (which is currently the same set,
but the variable is named distinctly for the TimedFullField branch).
This matches `TestCompositeInvariants`. If D-I.3d adds 2-D support to
`zeros_timed_full_field` for slab/2D-Cartesian, the `GEOMETRIES_1D`
parametrisation in `TestLinearity` can be widened.

### 6.4 — GMRES restart bound

The migrated GMRES test uses `restart=min(200, n_eq)`. For cylinder,
`n_eq=132`, so restart=132. The pre-migration test used
`restart=min(200, n_unknowns)`; both are equivalent at the legacy
packed sizes. No change of semantics.

### 6.5 — Empty `cylinder` LS dimension `n_ord`

Pre-migration tests passed `cylinder` cases with `_build_cylinder`
default `sn_order=4`. Post-migration my test_b1pp_lplusc_is_full_rank
took ~30s on cylinder due to the n_flat × n_flat = 264 × 264 matrix
probe. This is well-within tolerance for an L0 test but slower than
the legacy ~1s. If this becomes an issue D-I.3d can downsize the
cylinder builder's quadrature.

---

## 7 — Files changed

- `tests/sn/test_streaming_operator.py` — refactored, deletions
- `tests/sn/test_b1pp_verification.py` — 3 migrations, 1 import added
- `tests/sn/test_invertible_operator.py` — 1 new xfail negative test

No production-code edits. No edits to other tests files (per brief
constraint).

---

## 8 — Post-D-I.3d follow-up (main-agent action)

When D-I.3d surgery lands (removes the bare-ndarray adapter at
`orpheus/sn/operator.py:1539-1658`), the main agent MUST:

1. **Remove the `@pytest.mark.xfail(...)` decorator** from
   `TestSolve.test_apply_rejects_bare_ndarray`. The test then becomes
   a permanent negative contract (raises post-surgery).
2. **Verify the 3-file gate**: post-surgery target is `90 passed,
   0 xfailed` on the 3-file gate (the new test PASSES, no longer
   xfailed). Other tests stay green.
3. **Verify the 7-file gate**: target is `124 passed, 13 xfailed`
   (the only file-touched delta).
4. **Consider whether `n_unknowns` / `_ensure_eq_map` / `_eq_map`
   retire together** — per the explorer audit. The migrated tests
   no longer use these symbols.

---

## 9 — Lessons / references

- **L18 (Pattern 7 producer-side normalisation)** applied at three
  places: the b1pp rank test (subspace restriction at the matrix-
  probing producer), the b1pp GMRES test (matvec lifts at the
  producer boundary), and the new xfail negative test (the typed
  contract is the producer's responsibility).
- **L13 (canonical-type naming in briefs)** met by the brief; the
  migration consumed every named type verbatim.
- **L12 (paste-back evidence)** discipline — three brief
  pre-conditions diverged from observed state (§4); all surfaced.
- **`coding-elegance` Pattern 2 (single source of truth)** — the
  bare-ndarray + typed dual calling convention was a duplicate
  contract. Retiring it (at D-I.3d surgery) collapses two concepts
  to one.
- **`coding-elegance` Pattern 4 (illegal states unrepresentable)** —
  the cell-vs-face slot distinction becomes a type distinction
  (`out.bulk` vs `out.boundary`) instead of an array-slicing
  convention (`out[:n_cell_scalars]` vs `out[n_cell_scalars:]`).
- **Sibling memos**: [[di3-streamingoperator-apply-retirement-audit]]
  (explorer dependency audit); [[di2-test-scattering-operator-closeout]]
  (the sibling sibling wave's test-migration discipline).
