---
name: di3-streamingoperator-apply-retirement-audit
description: Wave D-I.3 dependency audit — retiring the bare-ndarray adapter arm of StreamingOperator.apply (L1539-L1658) in orpheus/sn/operator.py. Inventories every production + test caller, the internal-orphan symbols that come along, and the readiness verdict.
metadata:
  type: project
---

# D-I.3 — `StreamingOperator.apply` bare-ndarray adapter retirement audit

**Surgery target**: `orpheus/sn/operator.py` L1539-L1658 — the
`if isinstance(psi, TimedFullField): return self._apply_timed_full_field(psi)`
KEEP-path (L1537-1538) survives; everything below it (decode-via-`solution_to_angular_flux*`,
wrap-into-`TimedFullField`, call `_apply_timed_full_field`, re-pack-via-`pack_with_traces`,
post-subtract `σ_t ⊙ ψ_cell`) RETIRES.

**Date**: 2026-05-29. Sibling commits: D-I.1 `c97897d`, D-I.2 `8a8ddbf`,
D-K.1 `400ca33`, D-K.2 `dadf4e8`.

---

## Section 1 — Production callers of `StreamingOperator.apply`

Inventory grouped by entry point. Routing column: T = `TimedFullField` (KEEP),
B = bare-ndarray (RETIREMENT-target), A = legacy `AngularFlux` (KEEP — the
docstring still lists this overload but no production caller uses it; absorbed
by the typed path).

| File:line | Context | Arg shape | Routes to |
|---|---|---|---|
| `orpheus/numerics/iteration.py:262` | `_default_keff_estimator` — `den = np.sum(L.apply(psi))` | Whatever `KEigenvalue.solve` was constructed against | **Not invoked from SN production** (SNSolver does NOT use `KEigenvalue`; `solver.py` has no `KEigenvalue(` construction site). Only test consumers. |
| `orpheus/numerics/iteration.py:454` | `SourceIteration._solve.<...>` — `rhs = self.F.apply(psi) + self.S.apply(psi) + q_ext` | n/a — calls `F.apply`/`S.apply`, NOT `L.apply` | n/a (SI calls `L.solve`, NOT `L.apply`) |
| `orpheus/numerics/iteration.py:668` | `KrylovAcceleration.solve` → `A_matvec` — `out = self.L.apply(psi) - self.S.apply(psi) - self.F.apply(psi)` | `_unravel_like(q_ext, psi_flat)` — when `q_ext` is `TimedFullField` (which it IS in `_solve_krylov`), this returns `TimedFullField` via the ravellable `from_flat` dispatch at L186-194 | **T** (KEEP) |
| `orpheus/numerics/iteration.py:937` | `KEigenvalue.solve` — `q_fission = self.F.apply(psi) / keff` | n/a — calls `F.apply`, NOT `L.apply` | n/a |
| `orpheus/sn/solver.py:_solve_krylov` (L599-735) | Wires `LC = L_leaf + C_t` → `KrylovAcceleration`; `q_ext_composite` is `TimedFullField` (L664-671) | `TimedFullField` at the boundary; `_unravel_like` preserves it | **T** (KEEP) — verified via the `_is_ravellable(TimedFullField)` protocol detection (`TimedFullField.from_flat` exists at `orpheus/transport/timed_full_field.py:422`) |
| `orpheus/sn/solver.py:_solve_si` (L470-595) | Wires `LC = L_leaf + C_t` → `SourceIteration`; passes `q_ext_composite` (`TimedFullField`, L542-552) | n/a — `SourceIteration` calls `L.solve`, NOT `L.apply` | n/a |
| `orpheus/sn/solver.py:1377-1410` | A debug/dispatch helper (`solver.py:1340-1410` block tracks "ReducedStreamingOperator retires G3f") — `L_leaf = StreamingOperator(...)` constructed but only as constructor sanity, no `.apply` invocation in the block above | n/a | n/a |

**Production verdict**: every production call site that ends up at
`StreamingOperator.apply` routes via `KrylovAcceleration.solve →
A_matvec` with `psi = _unravel_like(q_ext_composite, psi_flat) →
TimedFullField`. The bare-ndarray adapter at L1539-L1658 is
**zero-cited from production code**.

`compute_keff_from_residuals` at iteration.py:262 is unreachable from
SNSolver (which does not consume `KEigenvalue`); its tests
(`tests/numerics/test_iteration.py:414, 645`) build synthetic operator
triples — not `StreamingOperator` — so they don't drive the
bare-ndarray arm either.

---

## Section 2 — Test callers of `StreamingOperator.apply`

Tabulated per file. Argument-shape column tracks what is actually passed
to `L.apply(...)`.

| File | apply-call count | Arg shapes | Migration cost |
|---|---|---|---|
| `tests/sn/test_streaming_operator.py` | 34 raw matches; **8 actual `L.apply(...)` invocations** (L234, 243, 280, 292, 293 × 2, 396, 419, 440, 452, 481, 496) | Bare ndarray at L234, L243, L280, L292-293 (`_packed_psi`, `np.zeros_like(psi)`, linear-combo on bare ndarray); TimedFullField at L396, L419, L440, L452, L481, L496 (`_random_composite`, `zeros_timed_full_field`) | **MIGRATE 5 sites** (`TestApplyShape.test_apply_preserves_shape`, `test_apply_returns_packed_vector`, `TestLinearity.test_apply_zero_returns_zero`, `test_apply_is_linear`). Each is a `_packed_psi` swap to `_random_composite` + assertion rewrite from `out.ndim == 1` / `out.shape == psi.shape` to `isinstance(out, TimedFullField) and out.bulk.values.shape == state.bulk.values.shape`. Linearity test: same — use TimedFullField `__add__`/`__mul__` already exercised at L290-293 of `test_streaming_operator_decomposition.py`. |
| `tests/sn/test_invertible_operator.py` | 28 `StreamingOperator(...)` constructions; **`L.apply` direct on TimedFullField** at L268 only (`sum_out = L.apply(state) + C.apply(state)`, `state = _random_state` = TimedFullField). Other constructions feed `L+C → InvertibleOperator.solve` or `__add__`. | TimedFullField everywhere | **NO migration** — already typed |
| `tests/sn/test_b1pp_verification.py` | 7 `L.apply` calls; constructions at L126, 179, 234, 286 | **Bare ndarray** — L130-131 `matvec(psi: np.ndarray)`; L183-184 `psi = np.ones(n_unknowns); L.apply(psi)+C.apply(psi)`; L238-239 similar; L290 `psi = rng.standard_normal(eq_map.n_unknowns)`; L296 `pack_with_traces` decode/encode round-trip | **MIGRATE 3 tests** (`test_b1pp_lplusc_is_full_rank` L116, `test_b1pp_constant_flux_collapses_to_collision` L167, `test_b1pp_lplusc_gmres_converges_fp_noise` L224) **+ DECIDE on 1 test** (`test_b1pp_decode_encode_roundtrip` L279 — this is a *pure decoder-pair identity test* that uses `solution_to_angular_flux_with_traces` + `pack_with_traces` directly, NOT `L.apply`; it retires WITH the codec symbols at D-J). The 3 migration sites trade `np.ones(n_unknowns)` → `mesh.zeros_timed_full_field()` with bulk set to 1.0 and the cell-vs-face slot inspection becomes `out.bulk.values` and `out.boundary.face_view(name)`. GMRES test: switch matvec from packed to ravellable via `TimedFullField.to_flat`/`from_flat`. |
| `tests/sn/test_phase_c_gates.py` | 7 `StreamingOperator(...)`; `L.apply(...)` at L392 + indirect via `op.apply` where `op` is the operator being tested | TimedFullField everywhere (`_build_composite`, `_flat_psi_composite`, `_random_bulk`) | **NO migration** — already typed (D-K.5 migration 2026-05-29 documented in source) |
| `tests/sn/test_streaming_operator_decomposition.py` | 17 matches; `L.apply` at L173, 245, 316 | TimedFullField (`state = TimedFullField(...)`, lines 157-168, 234-241, 303-313) | **NO migration** — already typed (D-I.1 migration) |
| `tests/sn/test_2d_l2_matvec_correctness.py` | 5 matches; `L.apply(state).bulk.values` at L200 | TimedFullField (`mesh.zeros_timed_full_field()` + `rng.standard_normal` on `.bulk.values` and `.boundary.values`) | **NO migration** — already typed |
| `tests/sn/test_2d_l2_face_view_unit_source.py` | 8 matches; `L.apply(state)` at L124, 155, 191, 228 | TimedFullField (`_zero_state_with_unit_face` returns `TimedFullField`) | **NO migration** — already typed |
| `tests/sn/test_operators_apply_typed.py` | 5 matches; `Lpsi = L.apply(state)` at L215 | TimedFullField (`_random_state`) | **NO migration** — already typed |
| `tests/sn/test_unified_matvec_slab.py` | 0 `L.apply`; imports `build_equation_map`, `solution_to_angular_flux` at L53-54 | Tests `transport_operator_matvec_unified` directly (uses legacy codec only for fixture wiring) | **Codec consumers** — flagged for D-J, not D-I.3 |
| `tests/sn/test_unified_matvec_cylinder.py` | 0 `L.apply`; imports `build_equation_map_cylindrical`, `solution_to_angular_flux_cylindrical` at L47-48; uses `build_equation_map_cylindrical(...)` at L255, 303 | Tests `transport_operator_matvec_unified` directly | **Codec consumers** — flagged for D-J |
| `tests/sn/test_unified_matvec_sphere.py` | 0 `L.apply`; imports `build_equation_map_spherical` at L38; uses it at L137 | Tests `transport_operator_matvec_unified` directly | **Codec consumers** — flagged for D-J |
| `tests/sn/test_native_matvec.py` | 0 `L.apply`; imports `build_equation_map_spherical`, `build_equation_map_with_traces` at L48-49; uses `build_equation_map_with_traces` at L444, 460 | Tests `transport_operator_matvec_unified` directly | **Codec consumers** — flagged for D-J |
| `tests/sn/spatial/test_apply_matvec_cylinder_invariants.py` | 0 `L.apply`; imports `build_equation_map_cylindrical` at L56; uses at L102 | Direct codec consumer | **Codec consumer** — flagged for D-J |
| `tests/sn/test_fixed_source_g1.py` | 0 `L.apply`; uses **legacy codec symbols only as spy targets** (L356-384) to assert ZERO invocation along G1 Krylov path | n/a — sentinel test | **DELETE the spy block** at D-J (the assertion becomes trivially true when the symbols are gone) |
| `tests/sn/test_l1_standoff_slab_cylinder.py` | 0 `L.apply` calls (file uses `solve_sn` end-to-end); docstring mentions `StreamingOperator + CollisionOperator` for context only | n/a | **NO migration** |
| `tests/sn/l1_analytical/test_kinf_homogeneous_tolerance.py` | 0 `L.apply`; docstring mentions "GMRES on `L.apply` — no SI dominance ratio" for context only | n/a (uses `solve_sn`) | **NO migration** |
| `tests/sn/spatial/test_psi_half_angle_seed.py` | 0 `L.apply`; docstring mentions `SNStreamingOperator.apply_transpose` for context only | n/a | **NO migration** |
| `tests/numerics/test_iteration.py` | `KEigenvalue(` at L414, 645 (synthetic operator triples — no `StreamingOperator`) | n/a | **NO migration** |

**Net D-I.3 test-migration scope**: 5 sites in
`test_streaming_operator.py` + 3 tests in `test_b1pp_verification.py`
= **8 tests** require packed→TimedFullField rewrite. The other 7
test files that touch `L.apply` are already typed (D-I.1 / D-K.5
migration history).

Cardinality pin (post-migration verifier): `pytest tests/sn -k
"streaming_operator or invertible_operator or b1pp or phase_c or
streaming_operator_decomposition or 2d_l2 or operators_apply_typed"
--collect-only | wc -l` should match the pre-migration count.

---

## Section 3 — Internal-to-orphan-after-retirement (D-J targets)

Per `__all__` at `orpheus/sn/operator.py:95-106`, the candidates are:

| Symbol | Source location | Production callers after D-I.3 surgery | Test callers after D-I.3 surgery | Internal-orphan? |
|---|---|---|---|---|
| `EquationMap` (dataclass) | operator.py:116 | NONE (`StreamingOperator._ensure_eq_map` retires WITH the bare arm) | `test_b1pp_verification.py` decode/encode round-trip (L279); `test_unified_matvec_*.py` fixture wiring; `test_native_matvec.py`; `test_apply_matvec_cylinder_invariants.py`; `test_dag_walk.py` (`EquationMap.unknowns_at_cell_for_mask`) | **YES** — D-J retires from operator.py; tests either rewire (codec round-trip → BoundaryFlux+AngularFlux round-trip) or delete |
| `build_equation_map` | operator.py:215 | NONE | `test_unified_matvec_slab.py:53`; `test_fixed_source_g1.py:358,374` (sentinel spy only) | **YES** |
| `build_equation_map_spherical` | operator.py:470 | NONE | `test_unified_matvec_sphere.py:38,137`; `test_native_matvec.py:48`; `test_fixed_source_g1.py:356,368` (sentinel spy) | **YES** |
| `build_equation_map_cylindrical` (alias of `_spherical`) | operator.py:598 | NONE | `test_unified_matvec_cylinder.py:47,255,303`; `test_apply_matvec_cylinder_invariants.py:56,102`; `test_fixed_source_g1.py:357,371` (sentinel spy) | **YES** |
| `build_equation_map_with_traces` | operator.py:632 | NONE (only used INSIDE `_ensure_eq_map`) | `test_native_matvec.py:49,444,460` | **YES** |
| `solution_to_angular_flux` | operator.py:259 | **ONE — `StreamingOperator.apply` L1571** (inside the 2-D Cartesian bare-ndarray arm, L1563-1588). This call site DIES with D-I.3. | `test_unified_matvec_slab.py:54`; `test_fixed_source_g1.py:361,383` (sentinel spy) | **YES (after D-I.3 surgery)** |
| `solution_to_angular_flux_spherical` | operator.py:496 | NONE | `test_unified_matvec_sphere.py:20,76` (mentioned in docstring); `test_fixed_source_g1.py:359,377` (sentinel spy) | **YES** |
| `solution_to_angular_flux_cylindrical` (alias of `_spherical`) | operator.py:599 | NONE | `test_unified_matvec_cylinder.py:48` (imported, but the actual use is `solution_to_angular_flux_spherical` at L137 via the alias); `test_fixed_source_g1.py:360,380` (sentinel spy) | **YES** |
| `solution_to_angular_flux_with_traces` | operator.py:711 | **ONE — `StreamingOperator.apply` L1605** (inside the 1-D bare-ndarray arm, L1604-1608). DIES with D-I.3. | `test_b1pp_verification.py:43,292` (decode/encode round-trip) | **YES (after D-I.3 surgery)** |
| `pack_with_traces` | operator.py:782 | **ONE — `StreamingOperator.apply` L1644** (1-D bare-ndarray arm re-pack). DIES with D-I.3. | `test_b1pp_verification.py:281,296` (round-trip) | **YES (after D-I.3 surgery)** |
| `StreamingOperator._ensure_eq_map` | operator.py:1403 | Inside `StreamingOperator.apply` bare-ndarray arm at L1542 | `test_b1pp_verification.py:287` (`L._ensure_eq_map(ng=ng)`); `test_streaming_operator.py` does not appear to call it directly | **YES** — retire WITH the bare arm at D-I.3 |
| `StreamingOperator.n_unknowns` (cached_property) | operator.py:1433 | NONE | `test_b1pp_verification.py:128,181,236` (matvec size); `test_streaming_operator.py` via `_packed_psi` (L145) | **YES** — retires WITH the bare arm |
| `StreamingOperator._eq_map` (field) | operator.py:1401 | Internal | n/a | **YES** — retires WITH the bare arm |

**D-I.3 surgery removes**: the bare-ndarray adapter body L1539-L1658,
the `_ensure_eq_map` method (L1403-1431), the `_eq_map` field
(L1401), and the `n_unknowns` cached_property (L1433-1440). The
1-D 2-D-Cartesian dispatch branches both vanish; only the
`isinstance(psi, TimedFullField)` → `_apply_timed_full_field`
delegation survives.

**D-J retirement** (after D-I.3) removes the entire codec family from
`__all__` (operator.py:95-106), the dataclass `EquationMap` and the
two factory + decoder + encoder triples. The dependent tests
(`test_b1pp_decode_encode_roundtrip`, `test_native_matvec.py` codec
tests, `test_unified_matvec_*.py` codec-fixture wiring,
`test_fixed_source_g1.py` sentinel spy block,
`test_apply_matvec_cylinder_invariants.py`,
`test_dag_walk.py`) either rewire to TimedFullField or delete.

---

## Section 4 — `_apply_timed_full_field` body integration recommendation

After surgery the `apply` method becomes:

```python
def apply(self, psi: TimedFullField) -> TimedFullField:
    return self._apply_timed_full_field(psi)
```

**Recommendation: option (b) — merge inline**. Rename
`_apply_timed_full_field` to `apply`, paste its body inline, delete
the old `apply` entirely. Reasons (file:line evidence):

1. **No external caller** invokes `_apply_timed_full_field` directly
   — `grep -rn "_apply_timed_full_field" --include="*.py" orpheus/
   tests/` returns only the three internal references at
   operator.py:1538, 1554-1557, 1584 (all inside the soon-to-die
   bare-ndarray arm) plus `test_2d_l2_matvec_correctness.py:6` and
   `test_2d_l2_face_view_unit_source.py:38` in docstrings only.

2. **No subclass override** — `StreamingOperator` has no subclasses
   (verified — `grep -rn "class.*StreamingOperator" --include="*.py"
   orpheus/` returns only the one class definition at L1307 plus the
   retired-pre-D-K.2 `SNStreamingOperator` references in comments).

3. **No sister `_apply_*_cartesian` pattern requiring symmetry** —
   `_apply_2d_cartesian` is a 2-D-specific kernel called from inside
   the merged `apply` (via the `curv is None and ny > 1` branch at
   operator.py:1694-1696). It STAYS as a private helper because it
   IS a per-geometry kernel — the parallel public structure (1-D
   unified vs 2-D Cartesian) lives inside the merged `apply`'s
   geometry dispatch.

4. **Pattern 1 (single source of truth)** prefers one named entry
   point over a public/private twin where the public is a trivial
   delegator.

After the merge, the geometry dispatch (L1692-1699 of
`_apply_timed_full_field`) is the only branching point in `apply`,
which matches the codebase's "build primitives not products"
elegance standard.

---

## Section 5 — Audit table summary (ONE-PAGE)

| Symbol | Prod callers | Test callers | Internal-orphan after D-I.3? | D-J? |
|---|---|---|---|---|
| `StreamingOperator.apply` bare-ndarray arm (L1539-L1658) | 0 | 8 (5 in `test_streaming_operator.py`, 3 in `test_b1pp_verification.py`) | n/a (the retirement target) | n/a |
| `StreamingOperator._apply_timed_full_field` | 0 (self only) | 0 (docstring mentions only) | Rename → `apply` at D-I.3 | n/a |
| `StreamingOperator._apply_2d_cartesian` | self (`_apply_timed_full_field`) | 0 | KEEP as private kernel | n/a |
| `StreamingOperator._ensure_eq_map` | self (bare arm) | 1 (`test_b1pp_verification.py:287`) | YES — retire at D-I.3 | n/a |
| `StreamingOperator.n_unknowns` | 0 | `test_streaming_operator.py:145` (via `_packed_psi`); `test_b1pp_verification.py:128,181,236` | YES — retire at D-I.3 | n/a |
| `StreamingOperator._eq_map` field | self | 0 | YES — retire at D-I.3 | n/a |
| `EquationMap` | 0 (after D-I.3) | ~10 sites (codec + spies + dag_walk + b1pp) | YES | **D-J target** |
| `build_equation_map` | 0 | 3 sites | YES | **D-J target** |
| `build_equation_map_spherical` | 0 | 5 sites | YES | **D-J target** |
| `build_equation_map_cylindrical` (alias) | 0 | 7 sites | YES | **D-J target** |
| `build_equation_map_with_traces` | 0 | 3 sites (`test_native_matvec.py`) | YES | **D-J target** |
| `solution_to_angular_flux` | 1 — operator.py:1571 (2-D bare arm) | 3 sites | YES (after D-I.3) | **D-J target** |
| `solution_to_angular_flux_spherical` | 0 | 3 sites | YES | **D-J target** |
| `solution_to_angular_flux_cylindrical` (alias) | 0 | 3 sites | YES | **D-J target** |
| `solution_to_angular_flux_with_traces` | 1 — operator.py:1605 (1-D bare arm) | 2 sites (`test_b1pp_verification.py:43,292`) | YES (after D-I.3) | **D-J target** |
| `pack_with_traces` | 1 — operator.py:1644 (1-D bare arm) | 2 sites (`test_b1pp_verification.py:281,296`) | YES (after D-I.3) | **D-J target** |

---

## Section 6 — Surgery readiness verdict

### **READY** for D-I.3d surgery — zero production-code blockers.

**Evidence**:

1. Every production caller path to `L.apply` runs through
   `_solve_krylov`'s `KrylovAcceleration` (`orpheus/sn/solver.py:702`)
   with `q_ext_composite: TimedFullField` (L664-671). The
   `_unravel_like(q_ext_composite, psi_flat)` call at
   `iteration.py:667` dispatches to `TimedFullField.from_flat`
   (verified at `orpheus/transport/timed_full_field.py:422`) — the
   matvec at L668 receives `TimedFullField`, so the `isinstance(psi,
   TimedFullField)` guard at operator.py:1537 fires, and the
   bare-ndarray adapter is dead code.

2. `_solve_si` calls `L.solve`, not `L.apply`, so it does not exercise
   `StreamingOperator.apply` at all (the matvec lives inside the
   sweep on the `InvertibleOperator.solve` side).

3. `SNSolver` does NOT consume `KEigenvalue`, so `iteration.py:262`'s
   `den = np.sum(L.apply(psi))` never reaches `StreamingOperator`
   from production code.

4. The remaining `solver.py:1377-1410` `StreamingOperator(...)`
   construction is for a debug/dispatch scaffold and does not invoke
   `.apply`.

### Recommended D-I.3 surgery sequence

1. **D-I.3c (test migration, in this commit or sibling)** — rewrite
   the 8 bare-ndarray test sites:
   - `test_streaming_operator.py`: replace `_packed_psi` with
     `_random_composite` in `TestApplyShape` (2 tests) and
     `TestLinearity` (2 tests); rewrite assertions from
     `out.shape == psi.shape` / `out.ndim == 1` to TimedFullField
     invariants matching the sibling
     `TestCompositeInvariants.test_returns_timed_full_field` (L385).
   - `test_b1pp_verification.py`: rewrite
     `test_b1pp_lplusc_is_full_rank` (L116),
     `test_b1pp_constant_flux_collapses_to_collision` (L167), and
     `test_b1pp_lplusc_gmres_converges_fp_noise` (L224) to:
     - replace `n_unknowns = L.n_unknowns` with
       `state = sn_mesh.zeros_timed_full_field()`;
     - matvec returns `TimedFullField`; cell vs face inspection via
       `out.bulk.values` and `out.boundary.face_view(name)`;
     - GMRES test: bridge through `TimedFullField.to_flat`/`from_flat`
       (or call the matvec on `TimedFullField` directly via
       `KrylovAcceleration`).
   - `test_b1pp_decode_encode_roundtrip` (L279) — DEFER to D-J
     (it's a codec round-trip pin, not a `L.apply` test).

2. **D-I.3d (surgery)** — at `orpheus/sn/operator.py`:
   - Delete L1539-L1658 (the bare-ndarray adapter body after the
     `isinstance` guard);
   - Delete `_ensure_eq_map` (L1403-1431), `n_unknowns`
     (L1433-1440), `_eq_map` field (L1401);
   - Merge `_apply_timed_full_field` (L1660-1713) inline as `apply`,
     drop the `isinstance` guard;
   - Tighten the type hint on `apply(psi: TimedFullField) ->
     TimedFullField`;
   - Update the docstring to drop the `np.ndarray | AngularFlux`
     overload language;
   - Leave `_apply_2d_cartesian` (L1715+) untouched.

3. **D-J (subsequent commit)** — retire the codec family
   (`EquationMap`, `build_equation_map*`, `solution_to_angular_flux*`,
   `pack_with_traces`); migrate the remaining test sites
   (`test_b1pp_decode_encode_roundtrip`, `test_unified_matvec_*.py`,
   `test_native_matvec.py`, `test_apply_matvec_cylinder_invariants.py`,
   `test_dag_walk.py`) to typed-field equivalents or delete; remove
   the `test_fixed_source_g1.py:333-411` sentinel spy block (becomes
   trivially-true).

### Notes for the implementer

- The `legacy AngularFlux` (mesh-bound, single-class) overload
  language in the apply docstring (operator.py:1463-1472) describes
  a third arm that the current source does NOT actually carry —
  there's no `if isinstance(psi, AngularFlux):` branch in the
  L1442-L1658 body. The docstring is stale; the surgery is a clean
  two-case → one-case collapse, NOT a three-case → one-case.

- `test_b1pp_verification.py:287` calls `L._ensure_eq_map(ng=ng)`
  inside the constant-flux test (L167-207) to compute
  `n_cell_scalars = eq_map.n_eq * ng` for slot inspection. The
  D-I.3c rewrite reads slot semantics off `BoundaryFlux.layout` and
  `out.bulk.values` directly — no `eq_map` lookup needed.

- Related memo: [[dh1b-legacy-angular-flux-consumer-audit]] —
  describes the parallel Phase A retirement of
  `orpheus.sn.angular_flux.AngularFlux` (legacy class). D-I.3 is
  Wave D-I (operator algebra leaves), Phase A is separate.

- Related memo: [[issue-196-phase-g-replan-algebra]] — establishes
  the four-operator algebra contract (`L`, `C`, `S`, `F` all return
  `TimedFullField`) that D-I.3 closes for `L`. D-I.1 closed it for
  `C`; D-I.2 closed it for `S`. After D-I.3, all four leaves are
  `TimedFullField → TimedFullField` and `OperatorSum.apply`'s
  dunder distribution is type-safe end-to-end.
