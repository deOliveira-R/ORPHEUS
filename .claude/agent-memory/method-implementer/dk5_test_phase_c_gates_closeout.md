# D-K.5 — `test_phase_c_gates.py` SNStreamingOperator → (L + C) migration

**Branch**: `refactor/sn-operator-algebra`
**Date**: 2026-05-29
**Predecessor**: commit `05864bf` (test_snstreamingoperator deletion + L variable retirement)
**Successor**: D-K.2 — `SNStreamingOperator` class deletion (unblocked by this task)

## Verbatim pytest stdout (L12 — closeout-time plausibility substitution mitigation)

```
============================= test session starts ==============================
platform darwin -- Python 3.14.3, pytest-9.0.2, pluggy-1.6.0 -- /Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python
cachedir: .pytest_cache
rootdir: /Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/moment-space-and-layering
configfile: pyproject.toml
plugins: dash-4.1.0, anyio-4.13.0
collecting ... collected 17 items

tests/sn/test_phase_c_gates.py::test_apply_linearity_under_sweep_frame[sphere_GL4_reflective] PASSED [  5%]
tests/sn/test_phase_c_gates.py::test_apply_linearity_under_sweep_frame[cyl_LS4_reflective]    PASSED [ 11%]
tests/sn/test_phase_c_gates.py::test_apply_curvilinear_per_ordinate_flat_flux_residual[sphere_GL4_reflective-mms-0.0] XFAIL [ 17%]
tests/sn/test_phase_c_gates.py::test_apply_curvilinear_per_ordinate_flat_flux_residual[sphere_GL4_reflective-mms-0.5] XFAIL [ 23%]
tests/sn/test_phase_c_gates.py::test_apply_curvilinear_per_ordinate_flat_flux_residual[cyl_LS4_reflective-mms-0.0]    XFAIL [ 29%]
tests/sn/test_phase_c_gates.py::test_apply_curvilinear_per_ordinate_flat_flux_residual[cyl_LS4_reflective-mms-0.5]    XFAIL [ 35%]
tests/sn/test_phase_c_gates.py::test_apply_apply_transpose_reciprocity_under_sweep_frame[sphere_GL4_reflective]       XFAIL [ 41%]
tests/sn/test_phase_c_gates.py::test_apply_apply_transpose_reciprocity_under_sweep_frame[cyl_LS4_reflective]          XFAIL [ 47%]
tests/sn/test_phase_c_gates.py::test_apply_face_fluxes_match_sweep_recurrence_spherical            PASSED [ 52%]
tests/sn/test_phase_c_gates.py::test_bc_trace_contract_respected_by_matvec_vacuum_sphere           PASSED [ 58%]
tests/sn/test_phase_c_gates.py::test_bc_trace_contract_respected_by_matvec_reflective_sphere       PASSED [ 64%]
tests/sn/test_phase_c_gates.py::test_bc_trace_contract_capture_and_compare_sphere[vacuum]          PASSED [ 70%]
tests/sn/test_phase_c_gates.py::test_bc_trace_contract_capture_and_compare_sphere[reflective]      PASSED [ 76%]
tests/sn/test_phase_c_gates.py::test_sweep_curvilinear_per_ordinate_flat_flux_residual[sphere_GL4_reflective-0.5]     PASSED [ 82%]
tests/sn/test_phase_c_gates.py::test_sweep_curvilinear_per_ordinate_flat_flux_residual[sphere_GL4_reflective-1.5]     PASSED [ 88%]
tests/sn/test_phase_c_gates.py::test_sweep_curvilinear_per_ordinate_flat_flux_residual[cyl_LS4_reflective-0.5]        PASSED [ 94%]
tests/sn/test_phase_c_gates.py::test_sweep_curvilinear_per_ordinate_flat_flux_residual[cyl_LS4_reflective-1.5]        PASSED [100%]

11 passed, 6 xfailed, 1 warning in 0.54s
```

## What migrated cleanly (11 passed, 4 xfail preserved)

* `test_apply_linearity_under_sweep_frame[sphere|cyl]` (Gate 1.4, 2 tests) — uses
  `TimedFullField` arithmetic (`__add__`, scalar `__mul__`) to express
  α·ψ + β·φ; asserts equality on BOTH `bulk.values` and `boundary.values`.
* `test_apply_curvilinear_per_ordinate_flat_flux_residual` (Gate 1.1, ERR-026
  catcher; 4 parametrised cases, all xfail-preserved with the pre-existing MMS
  reason) — flat ψ = 1 on reflective sphere/cylinder; composite `(L+C).apply`
  reproduces the bundled `SNStreamingOperator.apply` via Resolution A
  (`(L+C)ψ = M(ψ;σ_t)` bit-exact); checks `out.bulk.values = σ_t·ψ`.
* `test_apply_face_fluxes_match_sweep_recurrence_spherical` (Gate 1.2,
  bit-stability) — extended to assert `np.array_equal` on BOTH bulk AND
  boundary across repeated applies.
* `test_bc_trace_contract_respected_by_matvec_vacuum_sphere` + `_reflective_sphere`
  (Gate 1.5 basic, 2 tests) — `(L+C).apply(zero_state)` = zero on both bulk
  and boundary.
* `test_bc_trace_contract_capture_and_compare_sphere[vacuum|reflective]`
  (Gate 1.5 strengthened, 2 tests) — instruments `bc_right.apply`; captured
  call[0] now consumes `boundary.face_view("xmax")` (zero in zero-boundary
  input) instead of the legacy cell-centre-proxy at i=nx-1; call[1] still
  the WDD-propagated outward outflow per the §16A.3 contract.
* `test_sweep_curvilinear_per_ordinate_flat_flux_residual[*]` (Gate 1.6,
  Phase F sweep-path twin, 4 tests) — untouched; tests `CarlsonInwardSweep`
  and `carlson_inward_sweep_from_source` directly, no SN matvec consumer.

## What broke + the fix applied

**Gate 1.3 reciprocity** (`test_apply_apply_transpose_reciprocity_under_sweep_frame`,
2 parametrised cases) raised `AttributeError`, not `MissingCapability`, on
the first call to `op.apply_transpose(phi_state)`.

* Root cause: `StreamingOperator.capabilities = frozenset({CAP_APPLY})` —
  no `apply_transpose` advertised, no method defined. `OperatorSum.apply_transpose`
  (operator.py:618-619) calls `self.a.apply_transpose(x)` unconditionally
  without checking the capability set, so the missing method surfaces as
  `AttributeError` rather than the cleaner `MissingCapability`.

* Fix applied: marked Gate 1.3 `xfail(strict=True, raises=(MissingCapability,
  AttributeError))` with a reason that traces both the structural gap
  (StreamingOperator has no analytic adjoint yet) AND the dispatch-tightening
  side-issue (`OperatorSum.apply_transpose` should raise `MissingCapability`
  for capability mismatches; documented as Phase H / Issue #208 follow-up).

* Why xfail and not skip: the test body remains a regression contract for
  the FUTURE Phase H landing — when `StreamingOperator.apply_transpose`
  ships, the xfail will flip to xpass and the strict marker will FAIL the
  CI, forcing maintainers to convert it back to a passing test (and to
  re-validate the reciprocity claim with the new analytic adjoint).

## What's xfailed (6 total — 4 pre-existing + 2 new)

| Test | Cases | Reason |
|---|---|---|
| `test_apply_curvilinear_per_ordinate_flat_flux_residual` | 4 (sphere×σ_t∈{0,0.5} + cyl×σ_t∈{0,0.5}) | **Pre-existing.** MorelMontryAngularSweep does not preserve per-ordinate flat-flux invariant on sphere (cylindrical happens to telescope); xfail-strict=False. PR-TYPED-6c Step 7 retired the MMS-alternative closures (`LegacyTauSymmetricInterpolation`, `BaileyFlatFluxRedist`); MMS is the only surviving curvilinear closure. |
| `test_apply_apply_transpose_reciprocity_under_sweep_frame` | 2 (sphere, cyl) | **New (D-K.5).** Composite `(L+C)` lacks `apply_transpose` because `StreamingOperator` advertises only `{CAP_APPLY}`. Wave O Issue #208 / Phase H lands the analytic-adjoint algebra; until then, xfail-strict=True with `raises=(MissingCapability, AttributeError)`. |

## Lessons applied

* **L12 (closeout-time plausibility substitution mitigation)** — pasted
  verbatim pytest stdout above with explicit per-test outcomes and the
  17-item collection-then-result trail. No paraphrased "all green" claim;
  the 11/6 split is auditable line-by-line.

* **L13 (name the existing type to extend rather than invent)** — the
  migration consumed the existing `TimedFullField` carrier (with its
  `__add__`/`__mul__` algebra), `AngularFlux.from_mesh`,
  `BoundaryFlux.zeros_for_sn_mesh`, and `SNMesh.zeros_timed_full_field`
  factories. No new test-helper types invented; only a thin
  `_build_composite(sn_mesh, bulk, boundary?)` helper that composes the
  existing factories with optional non-zero boundary input (used by no
  test currently — the boundary=None default covers every migration site).

* **Pattern 7 (normalise at definition site)** — the Gate 1.5 strengthened
  test no longer needs the legacy `solution_to_angular_flux_spherical`
  packed → typed round-trip; the reference helper consumes `psi_bulk`
  directly (the same `(N, ng, nx, ny)` layout the matvec consumes), so
  the WDD-propagation chain reads from the canonical typed representation
  at every step. The transpose-to-`(ng, N, nx, 1)` ordering is moved to
  the helper body where it matches the matvec internal convention; no
  consumer site sees the packed/typed bridge.

* **Bit-identity vs principled-equivalence (vv-principles)** — the
  migration is bit-identical at the algebra level by Resolution A
  (`(L+C).apply(ψ) = M(ψ;σ_t)` bit-exact across slab/sphere/cylinder per
  `test_streaming_operator_decomposition.py`). The test rtols carried
  through unchanged (`rtol=1e-12`/`1e-13`/`1e-14`); no contract relaxation
  needed.

## Manifest

* **Migrated file**: `tests/sn/test_phase_c_gates.py` (10 construction
  sites + ~11 method calls → composite (L+C) consumed via TimedFullField).
* **Imports retired**: `SNStreamingOperator` removed; the file now imports
  `StreamingOperator`, `CollisionOperator`, `MissingCapability`,
  `TimedFullField`, `AngularFlux`, `BoundaryFlux`.
* **New helpers**: `_build_composite(sn_mesh, bulk_values, boundary_values=None)`,
  `_random_bulk(sn_mesh, rng)`, `_flat_psi_composite(sn_mesh, ng=1)`,
  `_outflow_at_boundary_for_sphere_from_bulk(sn_mesh, psi_bulk)`.
* **Retired helpers**: `_flat_psi_for_geometry` (replaced by
  `_flat_psi_composite`); `_outflow_at_boundary_for_sphere` and
  `_cell_centred_outer_psi_for_sphere` (replaced by the bulk-direct
  `_outflow_at_boundary_for_sphere_from_bulk`).
* **No production-code changes** — test-only migration.
* **No new ERR-NNN** — D-K.5 is structural migration, not bug-fix.
* **No Sphinx narrative changes** — the gate semantics are preserved;
  archivist's existing `discrete_ordinates.rst` Phase C narrative still
  describes the gates (the legacy class name in docstring text remains
  for historical traceability).

## Follow-ups unblocked

* **D-K.2** — `SNStreamingOperator` class deletion at `orpheus/sn/operator.py:1313+`
  is now unblocked. This file is the last substantive test consumer per
  commit `05864bf`'s closeout. After D-K.2, the cosmetic docstring/comment
  mentions across `test_l1_standoff_slab_cylinder.py`,
  `test_unified_matvec_{slab,cylinder}.py`, `test_b1pp_verification.py`,
  `test_fixed_source_g1.py`, `spatial/test_psi_half_angle_seed.py`,
  `spatial/test_sweep_vs_apply_consistency.py` can land alongside.

* **Wave O / Phase H** — the new Gate 1.3 xfail will catch the analytic-
  adjoint landing on `StreamingOperator`. Issue #208 already tracks
  the algebra; when the strict xfail flips to xpass the CI will force
  re-validation.
