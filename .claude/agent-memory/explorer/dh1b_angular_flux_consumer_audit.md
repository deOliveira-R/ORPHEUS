---
name: dh1b-angular-flux-consumer-audit
description: D-H.1b consumer migration audit — LEGACY orpheus.sn.angular_flux.AngularFlux usage sites; constructors, .boundary/.values/.history reads, Krylov flat-trace methods; migration mapping to new L2 AngularFlux + TimedFullField; sequence ordering with Phase A overlap risk.
metadata:
  type: project
---

# D-H.1b LEGACY AngularFlux Consumer Audit

Pre-implementation reconnaissance for Depth B step **D-H.1b**: migrate every
consumer of `orpheus.sn.angular_flux.AngularFlux` (legacy) to use:

- `orpheus.transport.fields.angular_flux.AngularFlux` (NEW pure-Field, no boundary, no history)
- `orpheus.transport.timed_full_field.TimedFullField` (NEW composite: bulk + boundary + history)

Source-of-truth files (already shipped):

- `orpheus/transport/fields/angular_flux.py` — pure-Field; 191 LOC; `values: NDArray`, `space: FunctionSpace`, `mesh: SNMesh`; `from_mesh`, `from_ndarray`, `zeros_for_sn_mesh`, `integrate_angular`, `N/ng/nx/ny` props; algebra inherited from `Field`.
- `orpheus/transport/timed_full_field.py` — composite container; 481 LOC; `bulk: Field`, `boundary: Field`, `_history`, `history_depth`; `advance(new_bulk, new_boundary)`, `at_lag(k)`, `history_length`, `to_flat`, `from_flat`, `copy`; algebra propagates to bulk + boundary.
- `orpheus/sn/geometry.py:809-847` — `SNMesh.zeros_timed_full_field(history_depth=2)` factory.
- `orpheus/numerics/iteration.py:163-231` — `_is_ravellable` / `_ravel` / `_unravel_like` / `_zeros_like` / `_l2_norm` already detect BOTH the new (`to_flat`/`from_flat`) and the legacy (`to_flat_with_traces`/`from_flat_with_traces`+`mesh`) protocols.

---

## 1. Legacy AngularFlux surface inventory

Legacy class at `orpheus/sn/angular_flux.py` (513 LOC). Public surface:

| Attribute / method | Purpose | Pure-Field carries? | Composite carries? |
|---|---|---|---|
| `values: NDArray` | Bulk flux values | YES (`new.values`) | YES via `state.bulk.values` |
| `mesh: SNMesh` | Phase-space carrier | YES (`new.mesh`) | YES via `state.bulk.mesh` |
| `space: FunctionSpace` | (not on legacy) | YES (new attribute) | n/a (per-member) |
| `boundary: BoundaryFlux` | Face state | NO — REMOVED | YES (`state.boundary`) |
| `_history: tuple` | Shift register | NO — REMOVED | YES (`state._history`) |
| `history_depth: int` | Buffer cap | NO — REMOVED | YES (`state.history_depth`) |
| `__call__(lag=0)` | Read lag-k frame | NO | YES (`state.at_lag(lag)`) — returns `TimedFullField`, not `AngularFlux` |
| `stash(new)` / `__lshift__` | Push new head | NO | YES (`state.advance(new_bulk, new_boundary)`) |
| `__len__` | 1 + history count | NO | YES via `1 + state.history_length` |
| `integrate_angular() → ScalarFlux` | Σ wₙψₙ | YES (`new.integrate_angular()`) | call `state.bulk.integrate_angular()` |
| `at_ordinate(n) → (ng, nx, ny)` | psi.values[n] | NO (currently) | n/a |
| `linf() → float` | max abs | YES (inherited Field method) | via `state.bulk.linf()` |
| `copy()` | Deep copy | YES (inherited) | YES (`state.copy()` — drops history) |
| `from_flat_with_traces(flat, mesh)` | Krylov decode | NO | YES (`TimedFullField.from_flat(flat, template)`) — DIFFERENT signature (template, not mesh) |
| `to_flat_with_traces()` | Krylov encode | NO | YES (`state.to_flat()`) |
| `N/ng/nx/ny` properties | Metadata | YES | YES via `state.bulk.N` etc. |
| `_validate_partner` | Layer 1 isinstance | YES (inherited) | YES (composite-level check) |
| Dunders `+/-/* / / -` | Algebra | YES (inherited Field) | YES (composite distribute) |

Two non-trivial signature deltas to watch:

1. **`from_flat_with_traces(flat, mesh)` vs `TimedFullField.from_flat(flat, template)`** — the new factory needs an entire `TimedFullField` template (carries shape AND mesh AND history_depth), not just an SNMesh. The Krylov adapter at `iteration.py:212-216` already handles this for `to_flat`/`from_flat` flavor.

2. **`__call__(lag)` returns `AngularFlux | None` vs `at_lag(lag)` returns `TimedFullField`** — type semantics shift. Pre-D-H consumers that did `psi(1).values` now must do `state.at_lag(1).bulk.values`. Pre-D-H `psi(1) is None` becomes `lag > state.history_length` — semantic check, not None-return. NB: the new `at_lag` raises `IndexError` instead of returning `None` (see `timed_full_field.py:376-382`).

---

## 2. Production consumers (`orpheus/`)

### 2a. Constructor call sites `AngularFlux(...)` — producer pattern

| File:line | Pattern | Pre-D-H args | Post-D-H replacement |
|---|---|---|---|
| `orpheus/sn/scattering.py:950` | `AngularFlux(combined.values, mesh)` | values + mesh | `AngularFlux.from_mesh(combined.values, mesh)` — pure-Field result (no boundary; this is a SOURCE return per the #205 dimensional sin) |
| `orpheus/sn/fission.py:236` | `AngularFlux(per_ord.values, psi.mesh)` | values + mesh | same as above |
| `orpheus/sn/solver.py:529, 645, 1368` | `AngularFlux(q_ext_per_ord.values, sn_mesh)` | values + mesh — used as RHS for SI/Krylov | new-L2 `AngularFlux.from_mesh(values, mesh)` for the BULK; the composite `TimedFullField` is built where the BOUNDARY state attaches (currently `q_ext_typed.boundary` mutated at solver.py:539) |
| `orpheus/sn/solver.py:1015-1017` | `AngularFlux(angular_flux, sn_mesh, boundary=solver._boundary_flux)` | values + mesh + boundary | `TimedFullField(bulk=L2AngularFlux.from_mesh(angular_flux, sn_mesh), boundary=solver._boundary_flux, history_depth=2)` — this is the eigenvalue `Solution` construction site |
| `orpheus/sn/solver.py:1233` | `AngularFlux(angular, sn_mesh)` (initial_guess in SI loop) | values + mesh | `L2AngularFlux.from_mesh(angular, sn_mesh)` — but caller is sweep `initial_guess` so see operator.py:2659 |
| `orpheus/sn/solver.py:1262-1264` | `AngularFlux(angular, sn_mesh, boundary=solver._boundary_flux)` | values + mesh + boundary | `TimedFullField(...)` — fixed-source SI Solution |
| `orpheus/sn/operator.py:1321` | `AngularFlux(m_cell, sn_mesh, boundary=m_boundary)` | matvec RETURN | replace return type with `TimedFullField` for the residual; SEE Phase A overlap |
| `orpheus/sn/operator.py:1596` | `AngularFlux(psi_view, sn_mesh, boundary=boundary_legacy)` | proxy bundle for legacy `SNStreamingOperator.apply` | Phase A retires `SNStreamingOperator`; AVOID |
| `orpheus/sn/operator.py:2069` | `AngularFlux(psi_cell, sn_mesh, boundary=boundary_in)` | packed→typed adapter | Phase A overlap (`StreamingOperator.apply` bare-ndarray packed path) |
| `orpheus/sn/operator.py:2151` | `AngularFlux(cell_values, sn_mesh, boundary=result.boundary)` | `_apply_typed` return | the central typed-flux return of `StreamingOperator.apply` — Phase A WILL rewrite this |
| `orpheus/sn/operator.py:2308-2310, 2359-2361` | `AngularFlux(self.sigma[None,:,:,:] * psi.values, psi.mesh)` | `CollisionOperator.apply` typed branch | NEW: return `TimedFullField` (zero boundary) OR keep pure `L2AngularFlux` if caller distributes externally. Decision pin: under `TimedFullField.__add__/sub`, `L.apply + C.apply` must both be `TimedFullField` for the dunder to compose. **Either both operators return `TimedFullField`, or the algebra fails the type check at composition.** |
| `orpheus/sn/operator.py:2672-2677` | `AngularFlux(angular, sn_mesh, boundary=boundary_buf, history_depth=rhs.history_depth)` | `InvertibleOperator.solve` return | `TimedFullField(bulk=L2AngularFlux.from_mesh(angular, sn_mesh), boundary=boundary_buf, history_depth=rhs.history_depth, _history=())` — the single sweep-return site |
| `orpheus/sn/geometry.py:779` | `AngularFlux(np.zeros(...), self)` | `zeros_angular_flux()` factory | retire factory; consumers call `mesh.zeros_timed_full_field()` (composite) or `L2AngularFlux.zeros_for_sn_mesh(mesh)` (pure bulk) |

### 2b. `psi.boundary` access — production

All in `orpheus/sn/`:

- **operator.py:1066-1068** — `transport_operator_matvec_unified` reads `boundary = psi.boundary; face_outer = boundary.xmax_face; face_inner = boundary.xmin_face`. This is THE 1-D matvec body; Phase A WILL rewrite. **AVOID.**
- **operator.py:1085** — `raise ValueError("Slab geometry requires psi.boundary.xmin_face to be populated …")`. Same matvec body, Phase A. AVOID.
- **operator.py:2075** — `result.boundary.xmax_face[eq_map.face_outer_ordinate, :]` (packed `_apply_typed` scatter). Phase A. AVOID.
- **operator.py:2080-2081** — same packed scatter, `xmin_face`. AVOID.
- **operator.py:2660** — `_copy_boundary_face_state(initial_guess.boundary, boundary_buf)`. This is `InvertibleOperator.solve` and IS in D-H.1b scope: `initial_guess` becomes a `TimedFullField`, so the line becomes `_copy_boundary_face_state(initial_guess.boundary, boundary_buf)` (no change in syntax since `state.boundary` is also a `BoundaryFlux`). **MIGRATE.**
- **operator.py:2662** — same call on `rhs.boundary`. `rhs` becomes a `TimedFullField`; `rhs.boundary` still resolves. MIGRATE (signature change only).
- **solution.py:230** — `return self.angular_flux.boundary`. This is the `Solution.boundary_flux` delegate property. Becomes `return self.angular_flux.boundary` where `self.angular_flux` is now a `TimedFullField` (its `.boundary` is the L2 `BoundaryFlux`). Property body unchanged; type signature changes. **MIGRATE (type-annot only).**
- **operator.py:2680-2697** — `_copy_boundary_face_state(src, dst)` helper. Operates on raw `BoundaryFlux`, no AngularFlux dependency. Inert.

### 2c. `psi.values` access — production

Sites that read `psi.values` where `psi` is a legacy AngularFlux (NOT a per-ordinate source nor a scalar flux nor a moment field):

- **scattering.py:619** — `values = angular_flux.values if is_typed else angular_flux` (in `build_aniso_source`). MIGRATE: `is_typed` check distinguishes; if `angular_flux` becomes `TimedFullField`, replace with `values = angular_flux.bulk.values if is_typed else angular_flux`.
- **scattering.py:652-654** — `mesh = angular_flux.mesh; moments_values = M.apply(values); …`. Same: `mesh = angular_flux.bulk.mesh`.
- **scattering.py:922** — `phi = psi.integrate_angular()` — `psi` is the `AngularFlux` `singledispatch` argument. The whole `@apply.register def _(self, psi: AngularFlux) → AngularFlux` BLOCK has to redispatch on `TimedFullField`. **DECISION POINT**: either (a) keep singledispatch on `AngularFlux` and operate on `psi.bulk` internally, or (b) rewrite to dispatch on `TimedFullField`. Recommended (b) for read-as-the-math.
- **solver.py:579, 695, 1399** — `psi_typed.integrate_angular().values`. MIGRATE to `state.bulk.integrate_angular().values`.
- **solver.py:1225** — `total_values = iso_per_ord.values + external_source`. `iso_per_ord` is a `PerOrdinateSource`, NOT an AngularFlux. Inert (only mentioned because grep flagged it).
- **operator.py:946** — docstring `psi.values` reference. Update.
- **operator.py:975-978** — docstring. Update.
- **operator.py:1009** — `psi_view = psi.values` in `transport_operator_matvec_unified`. Phase A. AVOID.
- **operator.py:1598, 2071** — packed adapters. Phase A. AVOID.
- **operator.py:2146** — `result.values - self.sigma_t[None,:,:,:] * psi.values` (`_apply_typed`). Phase A. AVOID.
- **operator.py:2309** — `self.sigma[None,:,:,:] * psi.values` (`CollisionOperator.apply` typed). MIGRATE (this is in D-H.1b scope — `CollisionOperator` is type-clean for the rewrite; under composite typing it becomes `self.sigma[None,:,:,:] * psi.bulk.values`).
- **operator.py:2360** — `q.values / self.sigma[None,:,:,:]` (`CollisionOperator.solve` typed). MIGRATE — same shape change.
- **operator.py:2631** — `rhs.values IS per-ordinate density` (docstring + `PerOrdinateSource.from_mesh(rhs.values, sn_mesh)`). MIGRATE: `rhs.bulk.values`.
- **operator.py:2639** — `source = PerOrdinateSource.from_mesh(rhs.values, sn_mesh)` — inside `InvertibleOperator.solve`. MIGRATE.
- **fission.py:225, 234, 236, 257** — `psi.integrate_angular()`, `fission_iso.values`, `psi.mesh`. MIGRATE (parallel to scattering.py).
- **sweep.py:554** — `initial_guess.values.transpose(1, 0, 2, 3)[..., 0]`. Phase A owns `sweep.py`. AVOID.
- **sweep.py:598** — `initial_guess.values[global_n, :, 0, 0]`. Phase A. AVOID.

### 2d. `psi(lag)` / `psi << new` / `psi.stash(new)` / `len(psi)` — production

Only ONE production reference and it's a DOCSTRING ghost — the runtime stash/lag plumbing was retired in Phase 1.2 (per the comment in `orpheus/sn/operator.py:2594-2602`):

- `operator.py:2474` — docstring `rhs(1)` reference to retired plumbing. **Docs-only edit.**
- `operator.py:2595` — comment about retired `rhs(1)`. Docs-only.
- `iteration.py:82-86` — docstring about retired `stash`. Docs-only.

**No production runtime call sites for stash / __lshift__ / __call__(lag) / __len__.** This is the cleanest possible D-H.1b finding: the history-shift functionality currently has ZERO production consumers. The new `TimedFullField.advance` / `.at_lag` will land with ONLY test consumers until BDF kinetics work needs them.

### 2e. `psi.integrate_angular()` / `psi.copy()` / `psi.linf()` — production

- `psi.integrate_angular()` — 4 production call sites: `scattering.py:922`, `fission.py:225`, `solver.py:579, 695, 1399`. ALL feed `.values` directly; consumer pattern is `psi_typed.integrate_angular().values`. Under composite typing: `state.bulk.integrate_angular().values`.
- `psi.copy()` — ZERO production call sites. (The `.copy()` matches in operator/sweep are all on `np.ndarray`, not AngularFlux.)
- `psi.linf()` — ZERO production call sites. Only test-side.

### 2f. `psi.from_flat_with_traces(...)` / `psi.to_flat_with_traces()` — production

Two production sites; BOTH inside Phase-A territory:

- `operator.py:2130` — `flat_in = psi.to_flat_with_traces()` (`StreamingOperator._apply_typed` 2-D Cartesian flat round-trip). Phase A absorbs the 2-D path. AVOID.
- `operator.py:2132` — `return AngularFlux.from_flat_with_traces(flat_out, sn_mesh)`. Same Phase-A territory. AVOID.

The Krylov adapter (`iteration.py:_ravel`, `_unravel_like`, `_zeros_like`, `_l2_norm`) consumes the protocol abstractly — it already supports BOTH `to_flat_with_traces`/`from_flat_with_traces` AND the new `to_flat`/`from_flat`. No edit needed there; **migration happens by the new-L2 producer (the composite `TimedFullField`) showing up in the typed slot.** The legacy detector branch retires when no `AngularFlux` instances remain in flight.

### 2g. `psi.at_ordinate(n)` — production

ZERO production call sites. Only test-side (one site).

---

## 3. Test consumers (`tests/...`)

13 test files import legacy `AngularFlux` (5,300 LOC combined):

| File | LOC | Pre-D-H surface used | Migration burden |
|---|---|---|---|
| `tests/sn/test_angular_flux_with_boundary.py` | 612 | stash/__lshift__/__call__(lag)/__len__/history_depth/from_flat_with_traces/to_flat_with_traces/copy/boundary keyword | **HIGHEST** — these are exactly the legacy-only methods; tests should retire wholesale and be re-implemented against `TimedFullField` |
| `tests/sn/test_native_matvec.py` | 511 | boundary= kwarg constructor, .values, .mesh, .boundary | Phase A overlap (matvec internals) — AVOID until Phase A lands its rewrite |
| `tests/sn/test_invertible_operator.py` | 768 | AngularFlux constructor (no boundary), .mesh, .values, history_depth keyword, `rhs(1)` retired-test reference, .copy | HIGH — InvertibleOperator.solve return type changes; per-test ~10 sites |
| `tests/sn/test_collision_operator.py` | 395 | AngularFlux constructor, no boundary, .values, .mesh | LOW — operator returns typed flux; rewire to TimedFullField |
| `tests/sn/test_operators_apply_typed.py` | 355 | AngularFlux constructor, .values, .mesh, to_flat_with_traces/from_flat_with_traces, .integrate_angular | MEDIUM — Krylov round-trip tests; rewire flat-methods to TimedFullField.from_flat/to_flat |
| `tests/sn/test_typed_fields.py` | 293 | AngularFlux constructor, .values, .mesh, .at_ordinate | MEDIUM — `at_ordinate` retires; replace with direct `psi.values[n]` indexing or keep as derived method on new L2 AngularFlux if needed |
| `tests/sn/test_harmonic_moment_field.py` | 558 | AngularFlux constructor, .values, .integrate_angular | LOW — moment-field tests; AngularFlux is a fixture builder |
| `tests/sn/test_fixed_source_g1.py` | 405 | AngularFlux constructor (no boundary) | LOW — fixture builder |
| `tests/sn/test_solution.py` | 424 | AngularFlux constructor (sol.angular_flux field), .values | MEDIUM — Solution.angular_flux field type changes to `TimedFullField`; rewire field-access patterns |
| `tests/sn/test_streaming_operator_decomposition.py` | 339 | AngularFlux(boundary=) constructor | Phase A — AVOID |
| `tests/sn/test_krylov_curvilinear_precond_safety.py` | 255 | AngularFlux constructor, .integrate_angular, retired rhs(1) ref | MEDIUM |
| `tests/sn/_test_helpers.py` | 131 | `legacy_proxy_matvec` builds an AngularFlux(boundary=) | Phase A — AVOID (this helper exists for legacy `transport_operator_matvec_unified`; Phase A retires it) |
| `tests/numerics/test_iteration_angular_flux.py` | 254 | AngularFlux constructor, to_flat_with_traces / from_flat_with_traces, .mesh | **HIGHEST in numerics**: tests the ravellable protocol on the legacy class. Either (a) retire entirely now that `tests/transport/test_timed_full_field.py` covers the new protocol, or (b) keep until D-H.2 retirement of legacy class. Recommendation: tag with `pytest.mark.skip` post-D-H.1b consumers complete, then delete in D-H.2 alongside `orpheus/sn/angular_flux.py`. |

---

## 4. Migration mapping table (consolidated)

| Pre-D-H pattern | Post-D-H equivalent | Notes |
|---|---|---|
| `AngularFlux(arr, mesh)` (no boundary) | `L2AngularFlux.from_mesh(arr, mesh)` | pure bulk |
| `AngularFlux(arr, mesh, boundary=bf)` | `TimedFullField(bulk=L2AngularFlux.from_mesh(arr, mesh), boundary=bf)` | composite |
| `AngularFlux(arr, mesh, history_depth=d)` | `TimedFullField(bulk=…, boundary=mesh.zeros_boundary_flux(), history_depth=d)` | composite, default boundary |
| `psi.values` | `psi.values` if pure; `state.bulk.values` if composite | type-check call site |
| `psi.mesh` | `psi.mesh` (pure) or `state.bulk.mesh` (composite) | symmetric to .values |
| `psi.boundary` | `state.boundary` | composite-only |
| `psi.boundary.xmax_face` | `state.boundary.xmax_face` | `BoundaryFlux` API unchanged |
| `psi(lag)` | `state.at_lag(lag)` | returns `TimedFullField`, NOT `AngularFlux`; raises `IndexError` not returns `None` past depth |
| `psi.stash(new)` / `psi << new` | `state.advance(new_bulk, new_boundary)` | callers MUST split new into bulk+boundary |
| `len(psi)` | `1 + state.history_length` | head counts as 1 |
| `psi.integrate_angular()` | `state.bulk.integrate_angular()` (composite) or `psi.integrate_angular()` (pure) | both return `ScalarFlux` |
| `psi.copy()` | `state.copy()` (composite, drops history) or `psi.copy()` (pure Field) | dunders inherited |
| `psi.linf()` | `state.bulk.linf()` (composite) or `psi.linf()` (pure) | inherited from Field |
| `psi.at_ordinate(n)` | `psi.values[n]` (direct indexing) | not migrated; consumers retire to indexing |
| `psi.N/ng/nx/ny` | `psi.N/ng/nx/ny` (pure) or `state.bulk.N` (composite) | new L2 AngularFlux already provides |
| `psi.from_flat_with_traces(flat, mesh)` | `TimedFullField.from_flat(flat, template)` | SIGNATURE CHANGE: template not mesh |
| `psi.to_flat_with_traces()` | `state.to_flat()` | semantics equivalent |
| `psi.history_depth` | `state.history_depth` | direct field on composite |

---

## 5. Highest-risk migration site

**`orpheus/sn/operator.py:2502-2677` — `InvertibleOperator.solve`**.

Why this is the single hardest site:

1. **Cross-type signature change**: `rhs: AngularFlux` → `rhs: TimedFullField`; `initial_guess: AngularFlux | None` → `initial_guess: TimedFullField | None`; return type `AngularFlux` → `TimedFullField`. This is the central SI/Krylov sweep API — every consumer (SourceIteration, KrylovAcceleration, solver.py SI/Krylov inner loops, fixed-source paths) sees the type change.
2. **Boundary plumbing is load-bearing**: lines 2658-2663 thread `initial_guess.boundary` → `boundary_buf` → sweep mutation → return. The `_copy_boundary_face_state` helper preserves the persistent partner-flux state for reflective BCs. If the new `state.boundary` shape diverges from the old `psi.boundary`, the reflective-BC fixed point shifts.
3. **History contract preservation**: line 2676 reads `history_depth=rhs.history_depth`. Under new semantics this stays as `history_depth=rhs.history_depth` BUT the new return is a composite WITH `_history=()`. **The current `InvertibleOperator.solve` does NOT call `.advance()` — it returns a NEW state, NOT an advanced state.** The SI/Krylov outer loop is responsible for the advance call (or not — under the new design, the outer iteration may just construct a fresh `TimedFullField` and discard history). Decision pin for the implementer: does the outer loop need history threading, or is the per-iterate composite always history-empty?
4. **Cross-operator algebra requirement**: `L.apply(state) + C.apply(state)` must compose under `TimedFullField.__add__`. That requires `StreamingOperator.apply` and `CollisionOperator.apply` to BOTH return `TimedFullField`. If only `StreamingOperator` returns composite, the `+` raises TypeError. **EITHER all four operators (L, C, S, F) get rewired to composite return, OR none.** Phase A owns StreamingOperator, so D-H.1b can ONLY safely rewire CollisionOperator if Phase A coordinates. This is the cross-worktree pinch point.
5. **Phase A overlap with sweep.py**: `transport_sweep(source, sig_t, sn_mesh, boundary_buf, initial_guess=initial_guess)` at line 2664 passes the typed `initial_guess` through to sweep code that Phase A is rewriting. If `initial_guess` becomes a `TimedFullField`, sweep.py:552 `initial_guess.values.transpose(...)` breaks (becomes `initial_guess.bulk.values.transpose`). Phase A will fix sweep.py independently; D-H.1b should avoid changing the InvertibleOperator.solve call signature until Phase A consensus.

**Recommendation**: D-H.1b should NOT migrate `InvertibleOperator.solve` in isolation. Either (a) coordinate with Phase A so both worktrees rewire together, or (b) D-H.1b leaves InvertibleOperator.solve on legacy `AngularFlux` and migrates only the Solution dataclass field + scattering/fission apply variants; the InvertibleOperator/sweep migration becomes a D-H.1c follow-up after Phase A lands its sweep rewrite.

---

## 6. Recommended migration sequence

Given the Phase-A overlap, the safest order is **leaf-first** — files with no downstream consumers (only consumed by Solution/tests) first; "trunk" files (consumed by every solver path) last.

### Leaf files (independent — migrate first)

1. **`orpheus/sn/solution.py`** (sol-level): change `Solution.angular_flux: AngularFlux` → `Solution.angular_flux: TimedFullField`. Update `boundary_flux` property (line 230) annotation. Update `compare` (line 326-327) to read `self.angular_flux.bulk.values`. This is the SINGLE field-type swap that every consumer of `Solution` (tests, examples) is forced to honor.

2. **`orpheus/sn/fission.py`** (operator leaf): rewire `FissionOperator.apply(AngularFlux)` register (line 204-236) to dispatch on `TimedFullField` instead. Internal compute: `psi.integrate_angular()` → `psi.bulk.integrate_angular()`; return `TimedFullField(bulk=L2AngularFlux.from_mesh(per_ord.values, psi.bulk.mesh), boundary=psi.bulk.mesh.zeros_boundary_flux())`. Single registered branch, ~30 LOC.

3. **`orpheus/sn/scattering.py`** (operator leaf): same pattern as fission, rewire `ScatteringOperator.apply(AngularFlux)` register (lines 886-950) and `build_aniso_source` (line 619). ~60 LOC.

4. **`orpheus/sn/operator.py` — CollisionOperator only** (lines 2276-2388 — the typed branches of `apply`/`solve`/`apply_transpose`): rewire from `AngularFlux` to `TimedFullField`. Two lines change semantically: 2308-2310 (`apply`) and 2359-2361 (`solve`). **NOTE the cross-worktree pinch**: Phase A is NOT rewriting CollisionOperator (per cross-worktree note); this is safe to migrate. But the COMPOSITE algebra `L + C` requires Phase A to coordinate. Recommended: rewire CollisionOperator typed branches to return `TimedFullField` ONLY if Phase A's StreamingOperator typed branches do the same in the same merge. Otherwise migrate CollisionOperator only to accept `TimedFullField` while still returning `L2AngularFlux` for now (then bump to composite return when Phase A aligns).

5. **`orpheus/sn/geometry.py:771-782`** — `zeros_angular_flux()` legacy factory. Retire (replaced by `mesh.zeros_timed_full_field()` and `L2AngularFlux.zeros_for_sn_mesh(mesh)`). All callers are tests.

### Mid files

6. **`orpheus/sn/solver.py`** — pure type-annotation churn at 5 constructor sites (529, 645, 1015, 1233, 1262, 1368). Plus 3 sites reading `psi_typed.integrate_angular().values` (579, 695, 1399). The solver dispatches into `InvertibleOperator.solve` (which is the trunk site under Phase A) — solver migration depends on InvertibleOperator migration. **Defer until trunk is settled.**

### Trunk files (CROSS-WORKTREE — coordinate with Phase A)

7. **`orpheus/sn/operator.py` — StreamingOperator + InvertibleOperator** (1280-2677). Phase A owns the 1-D / 2-D matvec internals (`transport_operator_matvec_unified`, `_apply_typed`, `solution_to_angular_flux*`). D-H.1b should NOT migrate this independently.

8. **`orpheus/sn/sweep.py`** — Phase A owns. **AVOID.**

9. **`orpheus/numerics/iteration.py`** — `_is_ravellable` detector already supports both protocols (lines 163-194). The detector retires when no legacy AngularFlux instances remain in flight (D-H.2, not D-H.1b).

### Test files

10. After production sites #1-#6 land, migrate the LOW-burden test files (test_collision_operator, test_harmonic_moment_field, test_fixed_source_g1). These are fixture-pattern migrations: replace `AngularFlux(arr, mesh)` with `L2AngularFlux.from_mesh(arr, mesh)` or `mesh.zeros_timed_full_field()` per test intent.

11. Migrate MEDIUM-burden tests (test_typed_fields, test_solution, test_operators_apply_typed, test_krylov_curvilinear_precond_safety) once `Solution` type signature is final.

12. The HIGH-burden tests (test_angular_flux_with_boundary, test_invertible_operator, test_iteration_angular_flux) wait until D-H.2 to RETIRE — these test legacy-specific behavior (`stash`, `__call__(lag)`, `to_flat_with_traces`, `boundary=` kwarg, `history_depth` on `AngularFlux`). The new `TimedFullField` equivalents are already covered by `tests/transport/test_timed_full_field.py`; the legacy tests delete in D-H.2.

13. Phase-A overlap tests (test_native_matvec, test_streaming_operator_decomposition, _test_helpers.legacy_proxy_matvec) — AVOID. Phase A will handle.

---

## 7. Phase A overlap — files D-H.1b MUST AVOID

Per cross-worktree note: Phase A in `.claude/worktrees/r1-phase-a-dim-agnostic/` is rewriting these files. D-H.1b SHOULD NOT touch them:

- `orpheus/sn/sweep.py` — entirely. Phase A rewrites `transport_sweep`, `_sweep_2d_wavefront`.
- `orpheus/sn/operator.py` lines 980-1321 — `transport_operator_matvec_unified` body and the `_sweep_direction` nested closure. Phase A territory.
- `orpheus/sn/operator.py` lines 1454-1700 — `SNStreamingOperator` (entire class). Phase A retires.
- `orpheus/sn/operator.py` lines 1900-2150 — `StreamingOperator.apply`, `_apply_typed`, and the packed-format bare-ndarray branches. Phase A overlap.
- `orpheus/sn/operator.py` lines 2128-2132 — 2-D Cartesian flat round-trip with `to_flat_with_traces` / `from_flat_with_traces`. Phase A.
- `tests/sn/test_native_matvec.py` — exercises `transport_operator_matvec_unified` directly.
- `tests/sn/test_streaming_operator_decomposition.py` — exercises StreamingOperator typed apply.
- `tests/sn/_test_helpers.py` `legacy_proxy_matvec` — Phase A retires this helper.

**Safe-to-touch surface in operator.py for D-H.1b**:

- Lines 2169-2406 (`CollisionOperator` class — typed branches only).
- Lines 2412-2677 (`InvertibleOperator` class — but only with Phase A coordination per §5).
- Lines 2680-2697 (`_copy_boundary_face_state` helper — no AngularFlux dependency).

---

## 8. Gaps / open questions

1. **Cross-operator composite return type**: do `L.apply`, `C.apply`, `S.apply`, `F.apply` all switch to returning `TimedFullField`, or do they stay as `L2AngularFlux` returns and only the outer SI/Krylov loop wraps into `TimedFullField`? The `TimedFullField.__add__` Layer-1 check at `timed_full_field.py:225` rejects bare `L2AngularFlux + TimedFullField` — so the four operators must agree on return type. **DECISION REQUIRED** before D-H.1b starts the operator carves.

2. **Solution.angular_flux semantic**: the eigenvalue Solution carries the final converged angular flux. Does that flux carry boundary state (load-bearing for downstream consumers like reaction-rate computation against material boundary)? Currently `sol.boundary_flux` delegates to `sol.angular_flux.boundary` — this implies YES. So `Solution.angular_flux` becomes a `TimedFullField` (composite). Confirm with the user before migrating solution.py.

3. **History semantics on the composite return path**: currently `InvertibleOperator.solve` reads `rhs.history_depth` but returns a fresh AngularFlux. The new `TimedFullField` flow: does `solve` call `.advance()` on `rhs` to push the result, or return a fresh composite with `_history=()`? The cleanest semantic is **return a fresh composite with empty history** — history is iteration metadata and the outer loop is responsible for `.advance()`-ing.

4. **Krylov adapter retirement**: `iteration.py:163-194` supports both protocols. The legacy branch retires in D-H.2 when no `to_flat_with_traces`-bearing types remain. The audit confirms no production code path needs the legacy branch after D-H.1b lands; only `tests/numerics/test_iteration_angular_flux.py` does. Recommend marking that test file `pytest.mark.skip(reason="D-H.2 retirement")` in D-H.1b commit, deleting in D-H.2.

5. **No `at_ordinate` in new L2 AngularFlux**: the new pure-Field does NOT carry `at_ordinate`. One test uses it (`test_typed_fields.py:202`). Either (a) add `at_ordinate` to new L2 AngularFlux for ergonomics, or (b) migrate the test to direct `psi.values[n]` indexing. Lean toward (b) — `at_ordinate` is a 1-liner indexing convenience; explicit indexing reads as the math.
