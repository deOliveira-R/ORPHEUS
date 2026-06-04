# Principled Index Migration — `(N, nx, ny, ng) → (N, ng, nx, ny)`

**Tracking issue**: Issue #196 (Phase G replan, prerequisite for Step 3+4.b type refactor).
**Branch**: `refactor/sn-operator-algebra` from tip `c21c2ef`.
**Start date**: 2026-05-15.

---

## §1 Context

The SN four-operator algebra `(L+C-S-F/k)ψ = q` requires one canonical shape for ψ across L, C, S, F. Today the codebase has three incompatible shapes (packed-vector for L/C leaves, `(N, nx, ny, ng)` for S, `(nx, ny, ng)` for F). The fix is two-stage: (1) pick the principled axis order; (2) introduce typed `AngularFlux` / `ScalarFlux` dataclasses on that order.

The principled order is **`(N, ng, nx, ny)`** — energy g second, NOT trailing. Derivation:

| Index | Coupling | Justification |
|---|---|---|
| `n` (ordinate) | None for within-group (P0 only) | Block-diagonal at top; sweep iterates outermost; Krylov outer batch |
| `g` (group) | None for within-group | Block-diagonal at second level (within-group means no cross-group coupling) |
| `i, j` (cells) | Streaming dependency chain | Sequential dag_walk required (the only "must-iterate" axis) |

The cross-section convention flips alongside: `σ_t : (nx, ny, ng) → (ng, nx, ny)`.

This migration must land **before** the typed dataclasses (next phase) so the dataclasses are born under the principled layout.

---

## §2 Target

- `ψ` (angular flux): `(N, ng, nx, ny)` for 2D; `(N, ng, nx)` for 1D.
- `φ` (scalar flux): `(ng, nx, ny)` for 2D; `(ng, nx)` for 1D.
- `σ_t`, `σ_a`, `νΣ_f`, etc.: `(ng, nx, ny)` for 2D; `(ng, nx)` for 1D.
- All sweep / matvec / operator internals: principled layout throughout.
- All public APIs: principled layout.
- All tests: principled layout fixtures.

---

## §3 Performance analysis

### Win: joint-batch ordinate_scan

Current call site: `for global_n in chain_def.ordinates: ordinate_scan((nx, ng), …)` — one scan per ordinate.

For SLAB: ordinates within a chain are not coupled (no M-M angular thread). Joint-batch over `(chain_size, ng, nx)`: one scan call per chain (2 chains per slab problem). Saves ~N Python invocations.

For CURVILINEAR (sphere/cylinder): the M-M angular thread couples ordinates sequentially within a level. The per-cell scan stays per-ordinate. Future work (research-level): reformulate the M-M recurrence as a parallel-prefix scan over ordinates; would unlock joint-batch for curvilinear too. Not in this migration — flag for separate investigation.

For ALL geometries: joint-batch over `(ng, nx)` within each ordinate's scan is the current behavior — verified preserved.

### Cost: per-cell `(ng,)` slice becomes strided

Under the principled layout, `ψ[n, :, i, j]` is non-contiguous (stride = `nx·ny`). This affects:
- Fallback degenerate-cell strategy path (rare, e.g., cylindrical |η|<10⁻¹⁵). Negligible total cost (~0.3 ms/sweep at nx=160).
- Per-cell unit tests of `CellUpdate.residual` / `.update`. Not performance-critical.

The hot path uses `ordinate_scan` (batched), which never takes per-cell `(ng,)` slices. No hot-path penalty.

### Net projection

Wall-clock impact: ±1% on the full test suite. Possibly slight win for high-N curvilinear (LS-Sn pin cells) due to better cell-axis contiguity per `(n, g)`. Possibly slight loss for problems heavy in the degenerate fallback. Within noise.

---

## §4 Migration strategy

**Bottom-up with bridge transposes.** Each PR converts one layer; outer layers stay on the old layout via transposes at the new layer's boundary. The transpose burden shrinks as the boundary moves outward. After the final PR, no bridges remain.

This is the standard pattern for big shape migrations. It keeps every intermediate commit green.

Rule: each PR leaves `pytest tests/sn/regression/ -q` green at `rtol=1e-12`. The regression snapshots stay valid throughout; only the FINAL PR (where the public sweep API flips) regenerates them under principled layout (verified by a one-shot transpose check against the old snapshots).

---

## §5 Stages

### PR-INDEX-1 — ordinate_scan joint-batch + sweep internals

**Scope**: `orpheus/sn/sweep.py` (`_run_1d_sweep` body). Internal layout flips to `(N, ng, nx, ny=1)`. Public `transport_sweep` signature unchanged (transposes at entry + exit).

**Changes**:
- `_run_1d_sweep` builds internal `angular_flux: (N, ng, nx, 1)`, `scalar_flux: (ng, nx, 1)`.
- Source arrays `QV_iso: (ng, nx)`, `Q_aniso_1d: (N, ng, nx)`.
- CollisionCache fields are READ as `(N, ng, nx)` (transposed at access if cache internal stays old).
- ordinate_scan call for SLAB: joint-batched over `(chain_size, ng, nx)` — one scan per chain (2 chains, vs N/2 today).
- ordinate_scan call for CURVILINEAR: per-ordinate `(ng, nx)` batch (unchanged). Document the M-M coupling as the obstruction to joint-batch.
- Entry transpose: `Q (nx, ny, ng) → (ng, nx, ny)`, `Q_aniso (N, nx, ny, ng) → (N, ng, nx, ny)`.
- Exit transpose: `angular_flux (N, ng, nx, ny) → (N, nx, ny, ng)`, `scalar_flux (ng, nx, ny) → (nx, ny, ng)`.

**Tests**: existing regression snapshots stay bit-identical (entry/exit transposes round-trip). Add a microbench test pinning the joint-batch behavior for slab.

**Performance gate**: Per-sweep Python overhead reduction ~N µs for slab. Measure on `slab_2g_3reg_dd_n40` benchmark.

**LoC estimate**: ~100-150 in `sweep.py`.

**Sub-agent**: method-implementer.

### PR-INDEX-2 — CollisionCache + GeometryCoefficients internal layout

**Scope**: `orpheus/sn/spatial/sweep_cache.py`.

**Changes**:
- `CollisionCache` fields flip from `(N, nx, ng)` to `(N, ng, nx)` (`a_attenuation`, `inverse_denom`, etc.).
- `CollisionCache.from_geometry` constructor: σ argument input shape changes from `(nx, ng)` to `(ng, nx)` for 1D.
- `GeometryCoefficients`: most fields are `(N, nx)` or `(N, nx, ny)` — no group axis. Unchanged.
- `_run_1d_sweep` reads cache fields directly (no transpose needed since PR-INDEX-1 already aligns).
- Bridge transpose at `_ensure_coll_cache`: input `sig_t (nx, ny, ng)` from caller; transpose to `(ng, nx, ny)` before passing to `from_geometry`.

**Tests**: existing cache-invariance tests (`test_sweep_cache.py`). Update assertion shapes.

**LoC**: ~80-120 in `sweep_cache.py` + adjustments to consumers.

**Sub-agent**: method-implementer.

### PR-INDEX-3 — Cross-section layout flip (`Mixture`, `assemble_cell_xs`)

**Scope**: `orpheus/data/macro_xs/mixture.py`, `orpheus/data/macro_xs/assemble.py` (or wherever `assemble_cell_xs` lives), all consumers.

**Changes**:
- `Mixture.sigma_t`, `Mixture.nu_sigma_f`, etc.: stored shape changes from `(ng,)` per-material — no change (already group-only for per-material). The per-CELL arrays in `SNSolver.sig_t`, `nu_sig_f` etc. flip from `(nx, ny, ng)` to `(ng, nx, ny)`.
- `assemble_cell_xs` output: `(ng, nx, ny)`.
- Bridge transpose at the boundary between SNSolver attributes (new layout) and any external caller still on the old layout. After PR-INDEX-5, no external consumers — all internal.
- Remove the bridge transpose added in PR-INDEX-2 (since cache now consumes the new layout natively).

**Tests**: assertions on `solver.sig_t.shape`, fixture assertions.

**LoC**: ~150-200.

**Sub-agent**: method-implementer.

### PR-INDEX-4 — Operators flip (S, F, legacy `SNStreamingOperator`, Resolution-A leaves)

**Scope**: `orpheus/sn/operator.py`, `orpheus/sn/scattering.py`, `orpheus/sn/fission.py`.

**Changes**:
- `ScatteringOperator.apply`: internal layout flips from `(N, nx, ny, ng)` to `(N, ng, nx, ny)`. The Galerkin pipeline `R · Λ · M` einsums update their axis labels. Public signature changes to accept and return `(N, ng, nx, ny)`.
- `FissionOperator.apply`: internal `phi` shape flips from `(nx, ny, ng)` to `(ng, nx, ny)`. Public signature ditto.
- `SNStreamingOperator.apply` (packed-vector form): the `EquationMap` packing order updates so the packed vector corresponds to `(N, ng, nx, ny)` traversal. Equivalent — just relabel the eq_map axis convention.
- `SNStreamingOperator.solve`: bridge transpose at boundary (in/out flip).
- `StreamingOperator` (Resolution A): internal apply layout flips. `eq_map.ix`, `eq_map.iy` retain semantics; the packed `sigma_packed` is built from the new `(ng, nx, ny)` σ_t shape.
- `CollisionOperator` (Resolution A): same.

**Tests**: `test_scattering_operator.py`, `test_fission_operator.py`, `test_streaming_operator.py`, `test_collision_operator.py`, `test_streaming_operator_decomposition.py`. Update fixture shapes.

**LoC**: ~300-400.

**Sub-agent**: method-implementer.

### PR-INDEX-5 — `SNSolver` + scalar-flux consumers + regression snapshots

**Scope**: `orpheus/sn/solver.py`, `tests/sn/regression/`.

**Changes**:
- `SNSolver.scalar_flux`, `angular_flux` attributes flip to principled layout.
- `compute_keff`, `compute_group_production_rate`, `compute_group_absorption_rate` consume `(ng, nx, ny)` φ.
- `KEigenvalue` outer loop reads/writes φ in new shape.
- `solve_sn`, `solve_sn_fixed_source` return new shapes.
- Regression snapshots: load → transpose `(nx, ny, ng) → (ng, nx, ny)` for φ and `(N, nx, ny, ng) → (N, ng, nx, ny)` for ψ. Verify bit-identity with old snapshots after transpose. If exact, save the new snapshots and update load code to skip the transpose.
- Public `transport_sweep` signature: input `Q (ng, nx, ny)`, `Q_aniso (N, ng, nx, ny)`; output `(angular_flux, scalar_flux)` with new shapes.
- Remove all entry/exit transposes added in earlier PRs.

**Tests**: full regression suite at `rtol=1e-12`. The snapshot test compares principled-layout values; the migration claim is "principled-equivalence via pure transpose" per `vv-principles` §"Bit-identity vs principled-equivalence".

**LoC**: ~250-350.

**Sub-agent**: method-implementer.

### PR-INDEX-6 — Test fixtures cleanup + docs

**Scope**: All `tests/sn/` fixtures; `docs/theory/operator_algebra.rst`, `docs/theory/discrete_ordinates.rst`.

**Changes**:
- Test fixtures construct ψ, φ, σ in principled layout directly (no transposes).
- Documentation explains the principled order with the §1 derivation table.
- Add a `docs/theory/index_convention.rst` page with the canonical statement.

**LoC**: ~200 (mostly mechanical fixture updates).

**Sub-agent**: archivist.

---

## §6 Acceptance criteria

After all 6 PRs land:

1. **No bridge transposes remain.** `grep -rn "transpose\b.*\([1-3], 0\|0, [1-3]" orpheus/sn/` returns no migration-shaped transposes (some legitimate domain-specific transposes may remain — case-by-case review).
2. **All public APIs use principled layout.** `transport_sweep`, `solve_sn`, operator `apply`s all consume/return `(N, ng, nx, ny)` shapes.
3. **All 11 regression snapshots pass at `rtol=1e-12`** in principled layout.
4. **The L0 streaming-equilibrium gate stays green** (26/26 sphere + cylinder).
5. **Cross-section convention is consistent**: every `σ_t`, `σ_a`, `νΣ_f` is `(ng, nx, ny)` end-to-end.
6. **Performance gate**: full SN suite wall-clock within ±5% of pre-migration baseline.
7. **Cross-domain audit**: no SN-shape conventions leaked into `numerics/`, `geometry/`, `data/` (cross-section storage stays in `data/` but uses principled order natively).
8. **Documentation**: `docs/theory/index_convention.rst` exists and is cross-referenced from the operator-algebra and discrete-ordinates pages.

After acceptance, the type refactor (Issue #197 partial close + structural `AngularFlux` / `ScalarFlux` dataclasses) can land on the principled foundation.

---

## §7 Deferred work (not in this migration)

- **Joint-batch ordinate_scan for curvilinear** via parallel-prefix reformulation of the M-M angular thread. Research-level algorithm work; separate dispatch to numerics-investigator with algebra-of-record discipline. Estimated win: 3-10× sweep speedup on cylindrical pin-cell problems.
- **JAX migration** of the principled layout. The `(N, ng)` joint batch becomes `jax.vmap(scan, axes=(0, 0))` natively. Issue #197 Option D trigger.
- **GPU port**: `(N, ng)` leading batch maps to GPU grid dimension; `(nx, ny)` maps to block dimension. Natural fit; future Wave.

---

## §8 Sub-agent dispatch sequence

1. PR-INDEX-1: method-implementer.
2. PR-INDEX-2: method-implementer (after PR-INDEX-1 verified).
3. PR-INDEX-3: method-implementer.
4. PR-INDEX-4: method-implementer.
5. PR-INDEX-5: method-implementer.
6. PR-INDEX-6: archivist.

Each dispatch is gate-keeper-reviewed before the next.

---

## §9 Risk register

| Risk | Mitigation |
|---|---|
| Bridge transposes accumulate; cleanup forgotten | PR-INDEX-6 explicitly removes all bridges; grep-gate enforces |
| Regression snapshots drift on transpose | One-shot transpose check in PR-INDEX-5 verifies bit-identity |
| Performance regression on curvilinear | Performance gate per PR; back out if >5% slowdown |
| The 3+4.b type refactor needs more shape changes than anticipated | This migration locks the layout; type refactor adds dunder methods on top, no further shape changes needed |
| Eq_map (packed vector) ordering ambiguity | PR-INDEX-4 explicitly documents the ng→nx→ny traversal order for the packed vector |

End.
