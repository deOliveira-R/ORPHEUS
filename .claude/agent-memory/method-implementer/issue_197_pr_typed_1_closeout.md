---
name: issue-197-pr-typed-1-closeout
description: Issue #197 PR-TYPED-1 closeout. MaterialXSField introduced as the single source of truth for SN cross sections; 8 leaked per-material dispatch loops collapsed into typed verbs. Foundation for the typed-field contract resume (PR-TYPED-2).
metadata:
  type: project
---

# Issue #197 PR-TYPED-1 — `MaterialXSField` + 8-loop collapse

**Branch**: `refactor/sn-operator-algebra` from tip `c4c269d` (post PR-TYPED-0).
**Date**: 2026-05-16.
**Scope**: Introduce `MaterialXSField` typed wrapper; collapse 8 leaked per-material dispatch loops; rewire `SNSolver`/`ScatteringOperator`/`FissionOperator` to consume the typed XS field through one handle.
**Status**: STAGED (NOT committed; main-agent gate-keeping pending).

## §1 Goal

`SNSolver.__init__` previously held SEVEN separate XS attributes — `sig_t`, `sig_a`, `sig_p`, `chi` (per-cell broadcast) + `sig_s`, `sig2`, `sig_s0` (per-material dicts) — plus the parallel topology dict `_cells_by_mat` and the (unused) `_sig2_sum`. Every per-material operation reached for all four handles by name; the dispatch pattern `for mid, (ix, iy) in cells_by_mat.items(): ...` leaked out of `SNSolver.__init__` into 7 sites in `scattering.py` plus 1 in `solver.py:compute_group_production_rate` — eight identical dispatch loops, each an independent opportunity for the per-material structure to drift from the source-of-truth.

PR-TYPED-1 introduces `orpheus.sn.material_xs_field.MaterialXSField` — the typed wrapper that:
- Owns the per-material `Mixture` dict + the spatial `mat_map` as the SINGLE source of truth.
- Exposes per-cell views (`total_cross_section`, `absorption_cross_section`, `fission_production`, `emission_spectrum`) lazily-cached for the operators that need cell-grid layout.
- Exposes per-material accessors (`sig_s_legendre(mid)`, `n2n_matrix(mid)`, `fission_production_per_material(mid)`, `chi_per_material(mid)`) for operations that exploit per-material structure.
- Exposes typed verbs (`apply_p0_in_scatter`, `apply_n2n`, `apply_legendre_scattering_moments`, `add_n2n_to_group_rate`, `foldable_sig_s`, `residual_sig_s`, `is_p0_diagonal_with_zero_n2n`, `foldable_sigma`) — the lifted-up forms of the formerly-leaked dispatch loops. Each verb captures one piece of math; the per-material loop lives ONLY inside the verb's body.

Per `coding-elegance` Pattern 7 (normalize at definition site): the per-material dispatch loop pattern lives in **ONE place** now (the `material_xs_field.py` module). Per Pattern 3 (named intermediates): each verb names the operation in domain language (`apply_p0_in_scatter` reads as "the in-scatter source"). Per Anti-pattern 2 (no twin paths): the seven separate XS attributes on `SNSolver` collapse to one `self.mat_xs`.

## §2 Phase deliverables

### §2.1 New module `orpheus/sn/material_xs_field.py` (∼730 LoC)

- `MaterialXSField` dataclass wrapping `(materials, mesh)`.
- `from_mesh(sn_mesh)` standard factory.
- `_synthetic_for_tests(...)` factory for foundation tests that exercise the per-material dispatch in isolation (used by `test_legendre_moment_scattering.py` and `test_scattering_operator.py`).
- `with_overridden_sig_s_and_n2n(sig_s_dense, n2n_dense)` sibling constructor used by `ScatteringOperator.foldable_part()` / `residual_part()` to build derived operators with modified scattering data.
- 4 per-cell views (lazy, cached via `_ensure_cell_views`).
- 5 per-material accessors (`cells_by_material` cached, `sig_s_legendre`/`n2n_matrix` route through cached `_sig_s_dense`/`_n2n_dense` populated by `_build_dense_caches`, `fission_production_per_material`/`chi_per_material` direct read-throughs).
- 8 typed verbs encapsulating per-material dispatch.

### §2.2 `SNMesh.material_xs_field()` factory (`orpheus/sn/geometry.py`)

Lazy-import method returning `MaterialXSField.from_mesh(self)`. The lazy import avoids the circular dependency `geometry → material_xs_field → geometry`.

### §2.3 `SNSolver` collapses 7 XS attrs → 1 `mat_xs` (`orpheus/sn/solver.py`)

- New `self.mat_xs = sn_mesh.material_xs_field()` is the single attribute holding all XS state.
- The retired attributes (`self.sig_t`, `sig_a`, `sig_p`, `chi`, `sig_s`, `sig2`, `sig_s0`, `_cells_by_mat`) become THIN read-through `@property` definitions routing to `self.mat_xs.*` — marked TRANSIENT in their docstrings. `_sig2_sum` is dropped entirely (no consumer).
- `assemble_cell_xs` import retained — it's invoked once at `__init__` to run the `__debug__` cell-flattening invariant check.
- Both operator constructors (`ScatteringOperator.from_solver_data`, `FissionOperator.from_solver_data`) called with `mat_xs=self.mat_xs`.
- `rebind_cross_sections(new_sig_t)` updated to override `self.mat_xs._sig_t_cell` directly (so the cell-flattening invariant remains pinned).
- `compute_group_production_rate`'s (n,2n) per-material loop replaced by `self.mat_xs.add_n2n_to_group_rate(rate, flux_distribution, self.volume)`.

### §2.4 `ScatteringOperator` consumes `MaterialXSField` (`orpheus/sn/scattering.py`)

- `LegendreMomentScattering` constructor: `sig_s` + `cells_by_mat` → single `mat_xs: MaterialXSField` parameter. `apply` body becomes one line: `return self.mat_xs.apply_legendre_scattering_moments(moments, L=self.L, skip_l0=self.skip_l0)`.
- `ScatteringOperator` constructor: collapses 10 separate handles (`n_ordinates`, `nx`, `ny`, `ng`, `sig_s`, `sig2`, `sig_s0`, `Y`, `weights`, `cells_by_mat`) to 2 (`mat_xs`, `quadrature`).
- `n_ordinates` / `nx` / `ny` / `ng` / `weights` / `Y` / `sig_s` / `sig2` / `sig_s0` / `cells_by_mat` retained as `@property` read-throughs (TRANSIENT — kept for one cycle for the four `_build_rhs_*` helpers in `solver.py` that still consume `solver.sig_s` / `solver.sig2` directly; PR-TYPED-2 retires).
- `add_iso_source`/`add_n2n_source` become 1-line delegators to the typed verbs.
- `build_aniso_source` updated for the new `LegendreMomentScattering` signature.
- `foldable_part` / `residual_part` / `is_foldable_into_sigma_r` / `foldable_sigma` all delegate to `mat_xs.foldable_sig_s` / `residual_sig_s` / `is_p0_diagonal_with_zero_n2n` / `foldable_sigma` — encapsulating the dispatch loops.
- `foldable_part` / `residual_part` build derived `ScatteringOperator` siblings via `mat_xs.with_overridden_sig_s_and_n2n(...)`.
- `Y` is a lazy property — only built on first access when `scattering_order > 0`, freeing the P0-only path from paying the spherical-harmonics cost.

### §2.5 `FissionOperator` consumes `MaterialXSField` (`orpheus/sn/fission.py`)

- Constructor: `chi` + `sig_p` → single `mat_xs: MaterialXSField` parameter.
- `chi` / `sig_p` retained as `@property` read-throughs onto `mat_xs.emission_spectrum` / `mat_xs.fission_production`.
- `apply` body unchanged in semantics — reads `self.chi` / `self.sig_p` through the new properties.

### §2.6 Test fixture updates

- `tests/sn/test_legendre_moment_scattering.py`: 4 direct `LegendreMomentScattering(sig_s=...)` constructor calls migrated to use `MaterialXSField._synthetic_for_tests(...)` factory + new `LegendreMomentScattering(mat_xs=...)` API.
- `tests/sn/test_scattering_operator.py`: 4 direct `ScatteringOperator(n_ordinates=..., sig_s=..., ...)` constructor calls migrated to `ScatteringOperator(mat_xs=..., quadrature=_StubQuad(...), scattering_order=...)`. New `_StubQuad` class added inline (minimal AngularQuadrature stand-in supplying `.N`, `.weights`, `.spherical_harmonics(L)`).
- Other tests pass unchanged via the read-through TRANSIENT shims on `SNSolver` / `ScatteringOperator`.

## §3 Mechanism criteria (verbatim paste-back)

| # | Criterion | Evidence |
|---|---|---|
| 1 | `MaterialXSField` exists in `orpheus/sn/material_xs_field.py` | `grep -c "class MaterialXSField" orpheus/sn/material_xs_field.py` → **1** |
| 2 | `SNMesh.material_xs_field()` factory exists | `grep -c "material_xs_field" orpheus/sn/geometry.py` → **4** (TYPE_CHECKING import + factory method + signature + docstring reference) |
| 3 | `SNSolver` uses single `self.mat_xs` instead of 7 separate attrs | `grep -c "self.mat_xs" orpheus/sn/solver.py` → **28** sites consuming the single handle; the seven separate XS attributes are now `@property` read-throughs |
| 4 | All 8 per-material loops in `scattering.py` + `solver.py:422` REPLACED by `mat_xs.*` calls | `grep "for mid.*cells_by_mat\|for mid.*items()" orpheus/sn/scattering.py orpheus/sn/solver.py` returns **0 hits**. All 8 sites now route through typed verbs in `material_xs_field.py` (where the 8 sites collapse into 9 dispatch-loop call sites, but ALL inside the producer module per Pattern 7). |
| 5 | `FissionOperator` constructor takes `mat_xs` | `grep "mat_xs:" orpheus/sn/fission.py` → 1 dataclass field declaration |
| 6 | `ScatteringOperator` constructor takes `mat_xs` | `grep "mat_xs:" orpheus/sn/scattering.py` → 2 dataclass field declarations (`LegendreMomentScattering` + `ScatteringOperator`) |
| 7 | 11/11 regression PASS at rtol=1e-12 | `pytest tests/sn/regression/ -q` → **11 passed in 62.86 s** |
| 8 | L0 streaming-equilibrium 26/26 PASS | `pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q` → **26 passed in 958.88 s** |
| 9 | Operator suites PASS | `pytest tests/sn/test_scattering_operator.py tests/sn/test_fission_operator.py tests/sn/test_collision_operator.py tests/sn/test_legendre_moment_scattering.py tests/sn/test_snstreamingoperator.py tests/sn/test_streaming_operator.py tests/sn/test_streaming_operator_decomposition.py -q` → **214 passed in 1.06 s** |
| 10 | Full SN suite PASS (minus pre-existing) | (running — see §4 for details) |
| 11 | CP suite green | (running) |

## §4 Verification paste-back

### §4.1 Regression suite (load-bearing rtol=1e-12 gate)

```
$ .venv/bin/python -m pytest tests/sn/regression/ -q --no-header
...........                                                              [100%]
11 passed, 3 warnings in 62.86s (0:01:02)
```

The 3 warnings are pre-existing `RuntimeWarning: invalid value encountered in divide` on the P1-anisotropic snapshots (unchanged from PR-TYPED-0).

### §4.2 L0 streaming-equilibrium curvilinear (26/26 cases)

```
$ .venv/bin/python -m pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q --no-header
..........................                                               [100%]
26 passed, 1 warning in 958.88s (0:15:58)
```

### §4.3 Operator suites (full)

```
$ .venv/bin/python -m pytest tests/sn/test_scattering_operator.py \
    tests/sn/test_fission_operator.py \
    tests/sn/test_collision_operator.py \
    tests/sn/test_legendre_moment_scattering.py \
    tests/sn/test_snstreamingoperator.py \
    tests/sn/test_streaming_operator.py \
    tests/sn/test_streaming_operator_decomposition.py -q --no-header
214 passed, 1 warning in 1.06s
```

### §4.4 Pre-existing failures (NOT introduced by PR-TYPED-1)

Confirmed via `git stash` + re-run that the following failures are pre-existing on the PR-TYPED-0 baseline (tip `c4c269d`):

- `tests/sn/l1_analytical/test_kinf_homogeneous.py::test_kinf_homogeneous_spectrum[*]` — 6 cases fail. Per `principled_index_migration.md` and the PR-TYPED-0 closeout, the test's `mean(axis=(0, 1))` is wrong for principled `(ng, nx, ny)` layout. NOT a PR-TYPED-1 concern.
- `tests/sn/test_solver_components.py::TestTransportSweep::test_matches_saved_reference` — shape `(2, 6, 4)` vs saved reference `(6, 4, 2)`. The saved npy file is in legacy `(nx, ny, ng)` layout; the test was never regenerated after PR-INDEX-5. Verified pre-existing on PR-TYPED-0 baseline. NOT a PR-TYPED-1 concern.

### §4.5 Full SN suite (in flight)

Running in background while writing the closeout — will update §4.5 with verbatim paste-back upon completion.

### §4.6 CP suite (in flight)

Running in background while writing the closeout — will update §4.6 with verbatim paste-back upon completion.

## §5 Architecture rationale

### §5.1 Why two access modes (per-cell + per-material)

The XS field admits two access patterns at different operator layers:

- **Per-cell views** are what L, C, and F want — they index by `(g, ix, iy)` (or `(g, ix, iy, n)` for the streaming operator's matvec). `total_cross_section`, `fission_production`, `emission_spectrum` are the natural reads.
- **Per-material accessors** are what S wants — every Legendre order's scattering matrix is `(ng, ng)` per material, and the per-material decomposition into foldable/residual parts is a per-material algebraic operation. Forcing S through a per-cell view would require building a `(ng, ng, nx, ny)` array per Legendre order — a 4-D blow-up that wastes memory and obscures the per-material structure.

`MaterialXSField` provides both; the lazy caches mean a consumer that only needs per-cell views pays nothing for the dense Legendre cache, and vice versa.

### §5.2 Coding-elegance audit

- **Pattern 3 (named intermediates)**: `apply_p0_in_scatter`, `apply_n2n`, `apply_legendre_scattering_moments`, `add_n2n_to_group_rate`, `foldable_sig_s`, `residual_sig_s`, `is_p0_diagonal_with_zero_n2n`, `foldable_sigma` — each named operation captures one piece of math. The reader of `S.foldable_part()` sees `self.mat_xs.foldable_sig_s()` and immediately knows the operation's intent without descending into a loop body.
- **Pattern 4 (illegal states unrepresentable)**: every operator (L, C, S, F) consumes `mat_xs` from a single source. Constructing F with a `chi` that mismatched the materials' actual fission spectrum was previously possible (the constructor just took bare arrays); now it routes through `mat_xs.emission_spectrum` which is derived from `materials` — divergence is structurally unspellable.
- **Pattern 7 (normalize at definition site)**: the per-material dispatch loop pattern lives in ONE module (`material_xs_field.py`). The 8 sites in `scattering.py`/`solver.py` that previously each carried their own copy of the loop are now 1-line delegations.
- **Anti-pattern 2 (no twin paths)**: SNSolver previously carried 7 redundant XS attrs (4 per-cell + 3 per-material dicts + 2 cell-by-mat structures = 9 if you count). They collapse to one `mat_xs`. The transient `@property` read-throughs keep the call-site API working for one cycle (`solver.sig_t`, `solver.sig_s`, etc.).
- **Aggressive retirement (`feedback_aggressive_retirement`)**: the dict-based attributes (`sig_s`/`sig2`/`sig_s0`/`_cells_by_mat`/`_sig2_sum`) are RETIRED from `__init__` — they are no longer initialized as instance dicts. They survive only as `@property` accessors that rebuild on demand from `mat_xs`. PR-TYPED-2 will retire those properties entirely once the `_build_rhs_*` helpers consume `mat_xs` directly.

### §5.3 Composability framing

The user's framing of composability was load-bearing in the design:

- **Mixing** (weighted volume-average of two MaterialXSFields → one homogenised one): NOT implemented in this PR. The homogenisation step lives outside SN today; a future Wave would add `MaterialXSField.homogenise(weights)`.
- **Restriction** (per-region subset → MaterialXSField on the smaller domain): NOT implemented in this PR. Hook point for CP/MoC consumer patterns; future Wave.
- **Action** (`mat_xs.fission_production * scalar_flux` as a typed dunder): DEFERRED to PR-TYPED-2 which introduces `AngularFlux` / `ScalarFlux` frozen dataclasses with `__mul__` dunders.

The current PR delivers the foundation that makes all three composability operations buildable; the dunders themselves wait for the typed flux types.

## §6 Decisions made

- **Non-frozen dataclass**: `MaterialXSField` is a regular `@dataclass`, not `@dataclass(frozen=True)`. The lazy caches (`_sig_t_cell`, `_cells_by_mat`, `_sig_s_dense`, `_n2n_dense`) require mutation. The frozen-with-`object.__setattr__` workaround would obscure the lazy-cache idiom and complicate testing. The "frozenness" we want is content-immutability of `materials` and `mesh` — Python doesn't enforce that structurally, but every consumer respects it by convention.
- **Lazy-import of `MaterialXSField` from `SNMesh.material_xs_field()`**: necessary to avoid the circular dependency `geometry.py → material_xs_field.py → (TYPE_CHECKING) geometry.py`. The lazy import is colocated with the factory method.
- **`MaterialXSField._synthetic_for_tests`**: deliberately public-by-convention (single-underscore prefix). Foundation tests that exercise per-material dispatch in isolation NEED to build synthetic MaterialXSFields without paying the full `SNMesh`+`Mixture` construction cost. The factory bypasses the per-cell view machinery (which requires real `Mixture` data) and only populates the dense caches. The single-underscore prefix signals "test-only, do not consume from production code".
- **Transient `@property` shims on `SNSolver` + `ScatteringOperator`**: kept for ONE cycle per the brief's recommendation. PR-TYPED-2 retires them when the four `_build_rhs_*` helpers (which currently read `solver.sig_s` / `solver.sig2` directly) are rewired to read `solver.mat_xs.sig_s_legendre(mid)` / `solver.mat_xs.n2n_matrix(mid)` directly. The shim cost is one extra dict construction per access — negligible on the test-driven access pattern; the `_build_rhs_*` helpers don't hit it in a hot loop because they index per-cell, not per-material.
- **`Y` becomes a lazy property on `ScatteringOperator`**: under the old constructor, `Y` was a constructor parameter (computed eagerly in `SNSolver.__init__` only if `L > 0`). Under the new constructor, `Y` is built on first read via `self.quadrature.spherical_harmonics(L)` — moved into the property. This frees the P0-only path from constructing `Y` at all (a tiny optimization, but consistent with the lazy-cache discipline).
- **`with_overridden_sig_s_and_n2n` sibling factory**: chosen over per-Mixture-copying for `foldable_part`/`residual_part`. Building modified Mixtures would require constructing new `csr_matrix`-backed Mixtures and re-deriving every per-cell view — overkill for the algebraic foldable/residual split which only touches the dense Legendre matrices.

## §7 Open items / follow-up

1. **PR-TYPED-2 owns the typed-field contract**: `AngularFlux` / `ScalarFlux` / `BoundaryFlux` frozen dataclasses with `__mul__` dunders. The transient `@property` read-throughs on `SNSolver` + `ScatteringOperator` retire as part of that PR (consumers rewired to read `mat_xs.*` and the typed flux fields directly).
2. **PR-TYPED-4 owns `HarmonicMomentField`**: the typed wrapper for the moment tensor passed through `LegendreMomentScattering.apply`. The current bare-ndarray contract works fine for now.
3. **`_build_rhs_*` helpers in `solver.py`**: still consume `solver.sig_s` / `solver.sig2` via the transient read-throughs. PR-TYPED-2's typed-field migration will rewire these to consume `solver.mat_xs.sig_s_legendre(mid)` / `solver.mat_xs.n2n_matrix(mid)` directly. Then the read-through properties on `ScatteringOperator` and `SNSolver` can retire.
4. **`compute_macro_xs` producer side**: unchanged (per §C anti-recommendation #4). The producer still emits `(N_cells, ng)` flat XS arrays via `assemble_cell_xs`; PR-TYPED-1 is consumer-side. CP/diffusion/MoC also untouched.
5. **No new ERR-NNN**: PR-TYPED-1 is mechanical architectural plumbing. No bugs caught during the migration (the regression suite stayed bit-identical at rtol=1e-12 throughout).
6. **No 9th per-material loop discovered**: §H asked me to flag any 9th loop. Confirmed: 8 sites in `scattering.py` + `solver.py` were the complete inventory.

## §8 Files touched

### Production (`orpheus/`)
- `orpheus/sn/material_xs_field.py` (NEW, ∼730 LoC) — `MaterialXSField` typed wrapper.
- `orpheus/sn/geometry.py` — `SNMesh.material_xs_field()` factory + TYPE_CHECKING import.
- `orpheus/sn/solver.py` — `SNSolver.__init__` collapses 7 XS attrs → 1 `mat_xs`; transient `@property` shims for `sig_t`/`sig_a`/`sig_p`/`chi`/`sig_s`/`sig2`/`sig_s0`/`_cells_by_mat`; `compute_group_production_rate` uses `add_n2n_to_group_rate` verb; `rebind_cross_sections` updated.
- `orpheus/sn/scattering.py` — `LegendreMomentScattering` and `ScatteringOperator` consume `mat_xs` + `quadrature`; 4 internal methods delegate to typed verbs; `add_iso_source`/`add_n2n_source` become 1-line delegators.
- `orpheus/sn/fission.py` — `FissionOperator` consumes `mat_xs`; `chi`/`sig_p` become read-through properties.

### Tests (`tests/`)
- `tests/sn/test_legendre_moment_scattering.py` — 4 constructor sites migrated to `MaterialXSField._synthetic_for_tests` + new API.
- `tests/sn/test_scattering_operator.py` — 4 constructor sites migrated; `_StubQuad` helper added inline.

## §9 Self-improvement / skill notes

No new anti-pattern, no new ERR-NNN. PR-TYPED-1 is architectural plumbing per the brief. The closeout memo's manifest line will follow this template:

```
- [Issue #197 PR-TYPED-1 — MaterialXSField + 8-loop collapse](issue_197_pr_typed_1_closeout.md) — refactor/sn-operator-algebra 2026-05-16 (STAGED). New MaterialXSField type with per-material AND per-cell views; SNSolver collapses 7 XS attrs to one mat_xs; 8 leaked per-material loops replaced by typed verbs (apply_p0_in_scatter, apply_n2n, apply_legendre_scattering_moments, add_n2n_to_group_rate, foldable_sig_s, residual_sig_s, is_p0_diagonal_with_zero_n2n, foldable_sigma); thin sig_t/sig_a/sig_p/chi/sig_s/sig2/sig_s0/_cells_by_mat read-through properties on SNSolver + transient sig_s/sig2/sig_s0/cells_by_mat/Y on ScatteringOperator as one-cycle shims (retired in PR-TYPED-2). 11/11 regression PASS at rtol=1e-12 (62.86 s); 26/26 L0 streaming-equilibrium curvilinear PASS (958.88 s); 214/214 operator suites PASS. Foundation for typed-field contract resume per principled_index_migration.md §10.
```

## §10 Conventional commit message (proposed, staged not committed)

```
refactor(sn): MaterialXSField + collapse 8 per-material dispatch loops (Issue #197 PR-TYPED-1)

Introduces orpheus/sn/material_xs_field.py — MaterialXSField is the
typed wrapper over per-material Mixture data + spatial mat_map.
Exposes lazy per-cell views (total_cross_section, etc.) AND
per-material accessors AND typed verbs that encapsulate the
formerly-leaked per-material dispatch loops.

Eight leaked per-material loops collapse into typed verb calls:
- scattering.py:234 LegendreMomentScattering.apply →
  mat_xs.apply_legendre_scattering_moments
- scattering.py:405 add_iso_source → mat_xs.apply_p0_in_scatter
- scattering.py:426 add_n2n_source → mat_xs.apply_n2n
- scattering.py:578 foldable_part → mat_xs.foldable_sig_s
- scattering.py:628 residual_part → mat_xs.residual_sig_s
- scattering.py:702 is_foldable_into_sigma_r →
  mat_xs.is_p0_diagonal_with_zero_n2n
- scattering.py:737 foldable_sigma → mat_xs.foldable_sigma
- solver.py:429 (n,2n) in compute_group_production_rate →
  mat_xs.add_n2n_to_group_rate

SNSolver collapses 7 XS attributes (sig_t/sig_a/sig_p/chi +
sig_s/sig2/sig_s0) plus _cells_by_mat plus _sig2_sum to ONE
mat_xs attribute.  Eight read-through @property shims preserve
the old API surface for one cycle (retired in PR-TYPED-2).

ScatteringOperator + FissionOperator constructors collapse 10+
per-material handles to two typed parameters (mat_xs, quadrature).

Foundation for the typed-field contract resume per
principled_index_migration.md §10.  No production semantics
change: 11/11 regression PASS at rtol=1e-12; 26/26 L0
streaming-equilibrium curvilinear PASS; 214/214 operator suites
PASS.
```

## §11 Manifest line (for MEMORY.md)

```
- [Issue #197 PR-TYPED-1 — MaterialXSField + 8-loop collapse](issue_197_pr_typed_1_closeout.md) — `refactor/sn-operator-algebra` 2026-05-16 (STAGED). New `MaterialXSField` type with per-material AND per-cell views; `SNSolver` collapses 7 XS attrs to one `mat_xs`; 8 leaked per-material loops replaced by typed verbs (`apply_p0_in_scatter`, `apply_n2n`, `apply_legendre_scattering_moments`, `add_n2n_to_group_rate`, `foldable_sig_s`, `residual_sig_s`, `is_p0_diagonal_with_zero_n2n`, `foldable_sigma`); thin `sig_t/sig_a/sig_p/chi/sig_s/sig2/sig_s0/_cells_by_mat` read-through properties on SNSolver + transient `sig_s/sig2/sig_s0/cells_by_mat/Y` on ScatteringOperator as one-cycle shims (retired in PR-TYPED-2). 11/11 regression PASS at rtol=1e-12 (62.86 s); 26/26 L0 streaming-equilibrium curvilinear PASS (958.88 s); 214/214 operator suites PASS. `LegendreMomentScattering` constructor signature `(sig_s, cells_by_mat, L, skip_l0) → (mat_xs, L, skip_l0)`. `ScatteringOperator` constructor collapses 10 handles (`n_ordinates`, `nx`, `ny`, `ng`, `sig_s`, `sig2`, `sig_s0`, `Y`, `weights`, `cells_by_mat`) to 2 (`mat_xs`, `quadrature`). New `MaterialXSField._synthetic_for_tests` factory for foundation tests. Foundation for typed-field contract resume per `principled_index_migration.md` §10.
```
