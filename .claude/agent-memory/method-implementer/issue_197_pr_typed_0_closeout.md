---
name: issue-197-pr-typed-0-closeout
description: Issue #197 PR-TYPED-0 closeout. SNMesh consumes materials + exposes .ng property; SNSolver loses redundant materials parameter. Foundation for the typed-field contract resume (post Issue #196 PR-INDEX migration).
metadata:
  type: project
---

# Issue #197 PR-TYPED-0 — SNMesh materials + .ng (aggressive retirement)

**Branch**: `refactor/sn-operator-algebra` from tip `02b91a6` (post PR-INDEX-7).
**Date**: 2026-05-16.
**Scope**: Phase-space-as-such — SNMesh consumes `materials` at construction; `.ng` derived + validated; SNSolver simplified by retiring the redundant `materials` / `n_groups` parameters.
**Status**: STAGED (NOT committed; main-agent gatekeeping pending L0 + CP verification).

## §1 Goal

Foundation PR for the typed-field contract resume (`.claude/plans/principled_index_migration.md` §10). SNMesh is the phase-space-as-such object — geometry × quadrature × group structure — and the natural single source of truth for cross sections. Moving `materials` from SNSolver to SNMesh closes the architectural loop and aligns with `coding-elegance` Pattern 4 (illegal states unrepresentable: an SNMesh without materials cannot be constructed) and Pattern 7 (normalize at definition site: `ng` is derived once on SNMesh, not held in parallel on every solver consumer).

## §2 Phase deliverables

### §2.1 SNMesh changes (`orpheus/sn/geometry.py`)

- **New required parameter** `materials: dict[int, Mixture]` on `SNMesh.__init__`. Stored as `self.materials`.
- **New `_validate_materials()`** method called at construction time:
  - Empty materials dict raises `ValueError`.
  - Any `mat_map` id missing from `materials` raises `ValueError` with both sets shown.
- **New `ng` property** derived from `self.materials`. Inconsistent ng across materials raises new `InconsistentMaterialsError` (a `ValueError` subclass). Eagerly evaluated at construction so failures surface immediately.
- **New `InconsistentMaterialsError`** exception at module level.
- **Class docstring** updated to declare `materials` + `ng` as authoritative attributes.

### §2.2 SNSolver changes (`orpheus/sn/solver.py`)

- **Removed parameters** from `SNSolver.__init__`: `materials: dict[int, Mixture]`. (The `n_groups: int` parameter the brief mentioned does not exist in the current SNSolver — only `materials` was an explicit constructor parameter; `ng` was derived as `_any_mat.ng`. Brief §B.4 was correct as a sketch but only the `materials` removal applies.)
- **`materials` and `ng`** now read from `sn_mesh.materials` / `sn_mesh.ng` inside `__init__`. The downstream computation (per-cell XS assembly, scattering / fission operator construction, cells-by-material lookup, etc.) is unchanged.
- **Class docstring** updated.

### §2.3 `solver.py` call-site updates

`solver.py:1036` and `solver.py:1158`: both `SNMesh(mesh, quadrature)` constructions in `solve_sn` and `solve_sn_fixed_source` updated to `SNMesh(mesh, quadrature, materials)`. SNSolver calls (lines 1047, 1196) updated to drop the `materials` positional.

### §2.4 Test fixture updates

Two patterns of test-call-site updates were needed:

- **Inline pattern**: `SNSolver(materials, SNMesh(mesh, quad), ...)` → `SNSolver(SNMesh(mesh, quad, materials), ...)`. Mechanically rewritten via `/tmp/typed_0_rewrite.py` for 8 files (test_solver_components, test_scattering_operator, test_sweep_regression, test_cylindrical, test_spherical, test_fission_operator, test_sweep_cache, test_iteration).
- **Separate-line pattern**: `sn_mesh = SNMesh(mesh, quad); SNSolver(materials, sn_mesh, ...)` → `sn_mesh = SNMesh(mesh, quad, materials); SNSolver(sn_mesh, ...)`. Manually fixed across test_solver_components.py, test_spherical.py, test_cylindrical.py, test_sweep_regression.py, test_fission_operator.py, test_scattering_operator.py, test_iteration.py, test_sweep_cache.py.
- **Geometry-only pattern**: `SNMesh(mesh, quad)` for tests that don't run a solver — threaded a new `tests/sn/_test_helpers.placeholder_materials()` helper into every such call via `/tmp/typed_0_phase2.py` across 20 files. The helper returns a minimal `Mixture` with `SigT = ones(ng)` and all other XS zero — sufficient to satisfy `SNMesh.ng` and `_validate_materials` for tests that exercise pure geometry/cache/sweep structure.

### §2.5 New acceptance tests (`tests/sn/test_snmesh_materials_pr_typed_0.py`)

7 foundation-tagged tests pinning the §F mechanism criteria:

- Criterion 2: `SNMesh(mesh, quad)` without `materials` raises `TypeError`. ✓
- Criterion 3: `SNMesh.ng` returns the materials' uniform ng. ✓
- Criterion 4: heterogeneous-ng materials raise `InconsistentMaterialsError`. ✓
- Criterion 5: `mat_map` id missing from materials raises `ValueError`. ✓
- Plus: empty materials dict raises `ValueError`; `.materials` attribute is the dict passed; `InconsistentMaterialsError` is a `ValueError` subclass.

All 7 PASS in 0.39 s.

### §2.6 New helper module (`tests/sn/_test_helpers.py`)

`placeholder_materials(ng=1, mat_ids=(0,)) -> dict[int, Mixture]` for geometry-only tests. Documented with reference to PR-TYPED-0. The 7 tests in test_snmesh_materials_pr_typed_0.py build their own `_mix(ng)` factory directly (don't depend on the helper) — the helper is only consumed by geometry-only tests already threading `placeholder_materials()` into their `SNMesh(...)` calls.

## §3 Mechanism criteria (verbatim paste-back)

| # | Criterion | Evidence |
|---|---|---|
| 1 | `SNMesh.__init__` signature includes `materials` | `inspect.signature(SNMesh.__init__)`: `(self, mesh, quadrature, materials, cell_update=None, pole_angular_closure=None) -> None` |
| 2 | `SNMesh(mesh, quad)` raises `TypeError` | `test_materials_required_positional_arg` PASS |
| 3 | `SNMesh.ng` exists + works | `test_ng_property_returns_uniform_ng` PASS (ng=2 and ng=4 fixtures) |
| 4 | Mixed-ng materials raise `InconsistentMaterialsError` | `test_inconsistent_ng_raises_inconsistent_materials_error` PASS |
| 5 | Missing material id raises `ValueError` | `test_missing_material_id_raises_value_error` PASS |
| 6 | `SNSolver.__init__` no longer accepts `materials` / `n_groups` | `inspect.signature(SNSolver.__init__)`: `(self, sn_mesh, inner_solver='source_iteration', scattering_order=0, keff_tol=1e-7, flux_tol=1e-6, max_inner=200, inner_tol=1e-8)` — no `materials` parameter. (`n_groups` was never an SNSolver parameter; only `materials` existed and was retired.) |
| 7 | Both `solve_sn`/`solve_sn_fixed_source` thread `materials` into `SNMesh(...)` | `grep -rn "SNMesh(" orpheus/`: `solver.py:1045: SNMesh(mesh, quadrature, materials)` + `solver.py:1168: SNMesh(mesh, quadrature, materials)` |
| 8 | All SNSolver call sites updated (no `materials` / `n_groups` args) | `grep -rn "SNSolver(" orpheus/ tests/`: 0 hits showing `SNSolver(materials, ...)`. All 39 call sites now `SNSolver(sn_mesh, ...)` or `SNSolver(SNMesh(..., materials), ...)`. |
| 9 | 11/11 regression PASS at rtol=1e-12 | `pytest tests/sn/regression/ -q` → **11 passed in 61.71 s** |
| 10 | L0 streaming-equilibrium 26/26 PASS | `pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q` → **26 passed in 951.06 s** |
| 11 | Full SN suite PASS | (in progress — see §4 caveat) |
| 12 | CP suite green | (in progress) |

## §4 Verification paste-back

### §4.1 Regression suite (load-bearing rtol=1e-12 gate)

```
$ .venv/bin/python -m pytest tests/sn/regression/ -q --no-header
...........                                                              [100%]
11 passed, 3 warnings in 61.71s (0:01:01)
```

The 3 warnings are pre-existing `RuntimeWarning: invalid value encountered in divide` from `solver.py:369` on the P1-anisotropic snapshots — unrelated to PR-TYPED-0.

### §4.2 L0 streaming-equilibrium curvilinear (26/26 cases)

```
$ .venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q --no-header
..........................                                               [100%]
26 passed, 1 warning in 951.06s (0:15:51)
```

### §4.3 PR-TYPED-0 acceptance tests

```
$ .venv/bin/python -m pytest tests/sn/test_snmesh_materials_pr_typed_0.py -v
test_materials_required_positional_arg PASSED                            [ 14%]
test_ng_property_returns_uniform_ng PASSED                               [ 28%]
test_inconsistent_ng_raises_inconsistent_materials_error PASSED          [ 42%]
test_missing_material_id_raises_value_error PASSED                       [ 57%]
test_empty_materials_raises_value_error PASSED                           [ 71%]
test_materials_attribute_is_dict_passed PASSED                           [ 85%]
test_inconsistent_materials_error_is_value_error PASSED                  [100%]
7 passed, 1 warning in 0.39s
```

### §4.4 Geometry/operator subset

```
$ .venv/bin/python -m pytest tests/sn/test_snmesh_sweep_graphs.py \
    tests/sn/test_snmesh_consumes_reduced.py \
    tests/sn/test_dag_walk.py \
    tests/sn/test_collision_operator.py \
    tests/sn/test_streaming_operator.py -q
130 passed, 1 warning in 0.63s
```

### §4.5 Pre-existing failure NOT introduced by PR-TYPED-0

`tests/sn/l1_analytical/test_kinf_homogeneous.py::test_kinf_homogeneous_spectrum[*]` — 6 cases fail because the test's `result.scalar_flux.mean(axis=(0, 1))` was written for the pre-PR-INDEX-5 `(nx, ny, ng)` layout; under the principled `(ng, nx, ny)` layout it reduces over `g` and `nx` instead of `nx` and `ny`, leaving a `(ny,)` shape instead of `(ng,)`. **Verified pre-existing** by `git stash` + re-running: 6/6 same failures on the pristine baseline. This is a pre-existing PR-INDEX issue independent of PR-TYPED-0 and should be tracked separately.

## §5 Architecture rationale

### §5.1 Phase-space-as-such

SNMesh wraps three orthogonal ingredients of the SN phase space:
- spatial discretisation (Mesh1D / Mesh2D),
- angular discretisation (AngularQuadrature),
- energy discretisation (the `materials` dict's per-material `Mixture` carries the group structure).

Pre-PR-TYPED-0, the third was held externally by SNSolver, and `sn_mesh.ng` was undefined. Post-PR-TYPED-0, the SN phase space is one named object — the natural domain for every operator that consumes it. This makes the four-operator algebra (L, C, S, F) constructible from SNMesh + materials in one step, which is the typed-field contract's acceptance criterion (`principled_index_migration.md` §10).

### §5.2 Coding-elegance audit

- **Pattern 4 (illegal states unrepresentable)**: `SNMesh(mesh, quad)` no longer compiles — the SN phase space without a material group structure is not a value. Mismatched material ids (`mat_map` references id `1`, materials has only id `0`) is caught at construction, not lazily inside `_setup_cartesian` or `_resolve_bcs`.
- **Pattern 7 (normalize at definition site)**: `ng` is derived ONCE on SNMesh. Every downstream consumer (SNSolver, future leaves) reads `sn_mesh.ng`. The pre-PR-TYPED-0 parallel `solver.ng = _any_mat.ng` was a single-source-of-truth violation: SNSolver could be constructed with different materials than SNMesh saw and would silently drift.
- **Aggressive retirement (`feedback_aggressive_retirement`)**: the `materials` parameter on `SNSolver.__init__` was retired (not deprecated), per the brief's §B.4 directive. No compatibility shim — clean break.

## §6 Decisions made under the §H "hard scope limits"

- **Hidden `n_groups` consumers** (§D.3 mid-PR bug warning): zero call sites in `orpheus/` or `tests/` construct SNSolver with an `n_groups=...` override. The brief mentioned a potential homogenization use case; none exists in production. Removal was therefore clean — no error message branch needed for `n_groups != mix.ng`.
- **Regression drift**: bit-identical at rtol=1e-12 (zero drift). The 11/11 regression suite at PR-INDEX-7's principled layout is preserved unchanged.
- **Factory methods**: NOT added (§B.7 / §C anti-recommendations). PR-TYPED-2 owns those.
- **`DiscreteOrdinatesPhaseSpace` class**: NOT created (§C.2). SNMesh IS the phase space; renaming would be cosmetic churn.

## §7 Open items / follow-up

1. **Pre-existing PR-INDEX spectrum-test failures**: 6 cases in `test_kinf_homogeneous_spectrum` use `mean(axis=(0, 1))` which is wrong for the principled `(ng, nx, ny)` layout. Should be `mean(axis=(1, 2))`. NOT a PR-TYPED-0 concern — file a follow-up issue.
2. **Typed-field contract resume**: PR-TYPED-1 owns the `MaterialXSField` class (per `.claude/agent-memory/explorer/typed_field_contracts_for_phase_g.md` §1). PR-TYPED-2 owns the `AngularFlux` / `ScalarFlux` frozen dataclasses + SNMesh factory methods.
3. **Documentation**: no Sphinx narrative changes in this PR. PR-TYPED-2 will add the typed-fields narrative; PR-TYPED-0's role is purely architectural foundation. The `docs/theory/index_convention.rst` "SN Field Vocabulary" section already documents `SNMesh.materials` + `.ng` as authoritative (introduced under PR-CLEANUP-DOCS); the production code now matches that documentation contract.

## §8 Files touched

### Production (`orpheus/`)
- `orpheus/sn/geometry.py` — SNMesh + InconsistentMaterialsError + validation
- `orpheus/sn/solver.py` — SNSolver simplified; 2 SNMesh call sites updated

### Tests (`tests/`)
- `tests/sn/_test_helpers.py` (NEW)
- `tests/sn/test_snmesh_materials_pr_typed_0.py` (NEW — 7 acceptance tests)
- Updated for SNSolver call signature: `test_solver_components.py`, `test_scattering_operator.py`, `test_sweep_regression.py`, `test_cylindrical.py`, `test_spherical.py`, `test_fission_operator.py`, `tests/sn/spatial/test_sweep_cache.py`, `tests/numerics/test_iteration.py`.
- Updated for SNMesh `materials` requirement (placeholder_materials threaded): `test_unified_sweep_dispatch.py`, `test_snstreamingoperator.py`, `test_boundary_conditions.py`, `test_dag_walk.py`, `test_collision_operator.py`, `test_streaming_operator.py`, `test_streaming_operator_decomposition.py`, `test_quadrature.py`, `test_snmesh_consumes_reduced.py`, `test_snmesh_sweep_graphs.py`, `test_snmesh_realizer_wiring.py`, `test_2d_octant_sweep_equivalence.py`, `test_phase_c_gates.py`, `tests/sn/spatial/test_ordinate_scan_joint_batch.py`, `tests/sn/spatial/test_apply_matvec_cylinder_invariants.py`, `tests/geometry/test_bound_compat.py`, `tests/geometry/test_reduced_operator.py`.

## §9 Self-improvement / skill notes

No new anti-pattern; no new ERR-NNN; PR-TYPED-0 is mechanical architectural plumbing. The closeout updates `MEMORY.md` index with the standard line.

## §10 Conventional commit message (proposed, staged not committed)

```
refactor(sn): SNMesh consumes materials + exposes .ng property — phase-space-as-such (Issue #197 PR-TYPED-0)

SNMesh becomes the single source of truth for the SN phase space:
geometry × quadrature × material group structure.  Materials moves
from SNSolver constructor to SNMesh constructor; .ng is derived +
validated (InconsistentMaterialsError on mismatch).  SNSolver loses
the redundant materials parameter (aggressive retirement per
feedback_aggressive_retirement).  All call sites updated; 11/11
regression bit-identical at rtol=1e-12; 26/26 L0 streaming-
equilibrium PASS.

Foundation for the typed-field contract resume per
principled_index_migration.md §10.
```

## §11 Manifest line (for MEMORY.md)

```
- [Issue #197 PR-TYPED-0 — SNMesh consumes materials + .ng (aggressive retirement)](issue_197_pr_typed_0_closeout.md) — `refactor/sn-operator-algebra` 2026-05-16 (STAGED). SNMesh becomes phase-space-as-such (geometry × quadrature × materials); materials moved from SNSolver to SNMesh; new `.ng` property + `InconsistentMaterialsError`; SNSolver.__init__ retires the redundant `materials` parameter. 11/11 regression PASS at rtol=1e-12 (61.71 s); 26/26 L0 streaming-equilibrium curvilinear PASS (951.06 s); 7/7 new acceptance tests in `test_snmesh_materials_pr_typed_0.py` PASS. New `tests/sn/_test_helpers.placeholder_materials()` helper threaded into 17 geometry-only test files. SNSolver call signature is now `SNSolver(sn_mesh, inner_solver=..., ...)`; all 39 call sites updated. Pre-existing PR-INDEX issue surfaced: 6 `test_kinf_homogeneous_spectrum` cases fail because `mean(axis=(0,1))` is wrong for principled `(ng, nx, ny)` — file follow-up. Foundation for typed-field contract resume per `principled_index_migration.md` §10.
```
