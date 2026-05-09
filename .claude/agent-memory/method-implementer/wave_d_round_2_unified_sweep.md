---
name: Wave D Round 2 — unified sweep
description: SN reshape Wave D Round 2 (Issue #161) closeout — orpheus/sn/sweep.py rewritten as unified algorithm parameterized by cell_update; bit-identical to 11 frozen regression snapshots when cell_update=DiamondDifference (default).
type: project
---

`refactor/sn/unified-sweep` 2026-05-06 (descends from
`refactor/sn-operator-algebra` HEAD `2665ea3`). Wave D Round 2 rewrote
`orpheus/sn/sweep.py` (730 LOC, 4 dispatch paths) into a unified
algorithm parameterized by `cell_update: CellUpdate`.

**Why:** Cardinal Rule 2 architecture — collapse 4 string-equality
dispatch paths into one boolean-driven dispatch on
`sn_mesh.reduced.requires_upstream_angular_state`; expose
`SNMesh.cell_update` as the strategy slot for Wave C-extension's
LD/EC/Step rollout.

**How to apply:** When Wave C-extension lands LD/EC/Step strategies,
the unified sweep is ready — pass `cell_update=LinearDiscontinuous()`
at SNMesh construction. Two open design points remain:

1. **2-D wavefront cell-update parameterization** — kept as inlined DD
   for Wave D bit-identity. The Protocol's per-cell `(ng,)` shape
   doesn't accept `(n_diag, ng)` slices used by anti-diagonal
   vectorisation. Wave C-extension must extend the Protocol with a
   `update_batch(...)` method that accepts `(n_cells, ng)` shapes
   while preserving 1 ULP equality vs the per-cell `update`.

2. **`ReducedStreamingOperator.streaming_terms` cylindrical `abs_mu`
   bug** — extraction returns
   `abs_mu = abs(quadrature.mu_x[direction_idx])` for cylindrical
   geometry, but `direction_idx` for cylindrical is the WITHIN-LEVEL
   azimuthal index `m_local`, not the global ordinate index. The
   correct global ordinate is `level_idx[m_local]`. Worked around in
   the sweep via `_streaming_terms_with_abs_mu()`. Fix at the
   geometry layer needs to thread `level_indices` through the
   factory or accept the global ordinate as a separate argument.
   Pre-existing latent bug — masked by the Wave C DD test using
   synthetic streaming-terms with hand-set `abs_mu`.

**Bit-identity contract held:** 11/11 regression snapshots stay
`np.array_equal`. The default `SNMesh.cell_update = DiamondDifference()`
reproduces the inlined sweep math verbatim (Wave C verified via
hand-calc tests). The dispatch consolidation preserves 1 ULP equality.

**ERR-026 NOT closed by Wave D:** the 2 xfail-strict tripwires at
`tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
remain xfail. The curvilinear sweep's one-directional WDD closure
(the bug) is preserved bit-identically through DD's curvilinear
branch; Wave E (Issue #15) closes ERR-026 via solver-path migration
to Krylov-on-`apply` with the symmetric closure.

**Cylinder workaround scope:** The `_streaming_terms_for_inward_sweep`
swap helper handles the inward-vs-outward face-area convention. The
`_streaming_terms_with_abs_mu` helper handles the cylindrical
within-level vs global ordinate index. Both helpers are inline in
`orpheus/sn/sweep.py`; both should be retired when the geometry-layer
APIs are extended (Wave C-extension or Wave E).

**Files:**
- `orpheus/sn/sweep.py` — rewrite. ~730 LOC inlined → ~770 LOC dispatched
  (small growth from per-cell `streaming_terms()` + dispatch
  scaffolding; main code reduction is the consolidated dispatch).
- `orpheus/sn/geometry.py` — added `cell_update: CellUpdate` constructor
  argument with default `DiamondDifference()`.
- `tests/sn/test_unified_sweep_dispatch.py` — new. 9 foundation tests
  covering: dispatch routing (slab/sphere/cyl/2D), 1D cumprod fast-path
  preconditions (DD + GL1D + isotropic), default cell_update.
- `docs/theory/discrete_ordinates.rst` — extended with
  `unified-sweep-dispatch` section: dispatch rationale, cell-update
  Protocol parameter, 1D cumprod fast-path preconditions, 2D wavefront
  inlined-DD note, ERR-026 deferred-to-Wave-E note.
- `docs/verification/matrix.rst` — auto-regenerated (test count
  bumped 2758 → 2767).

**Tests passed (Wave D R2 verification gates):**
- `tests/sn/test_unified_sweep_dispatch.py` — 9/9 foundation (new)
- `tests/sn/spatial/` — 24/24 foundation (Wave C protocol + diamond)
- `tests/sn/test_sweep_regression.py` — 12/12 (geometry stencil reads)
- `tests/sn/regression/test_dd_regression.py` — **11/11 bit-identical**
  (gating contract held)
- `tests/sn/l1_analytical/` — 27/29 + 2 xfail intact (ERR-026 tripwires)
- `tests/derivations/test_sn_mms_anisotropic_symbolic.py` — 12/12
- `tests/sn/test_mms.py` + `test_mms_aniso.py` + `test_mms_2d.py` —
  4 + 1 + 3 (slab + 2D Cartesian MMS, including 2g hetero ~86s)
- `tests/sn/test_quadrature.py` + `test_boundary_conditions.py` +
  `test_properties.py` + `test_snmesh_consumes_reduced.py` — 80
  passing (Wave D R1 read-site tests + BC/quad tests)
- `tests/sn/test_scattering_operator.py` + `test_fission_operator.py`
  — 27/27 (Wave D R1 operator tests)
- `tests/sn/test_sweep_operator_inconsistency.py` — 4/4 (D3 staying
  as-is, Wave E rewriting)
- `sphinx-build -W -q` — exit 0
- `python -m tests._harness.audit` — orphan/ERR coverage unchanged
  (23 orphan equations + ERR-020 / ERR-031 missing — pre-existing)

**Pre-existing failures NOT introduced by this round:**
- `tests/sn/test_mms_curvilinear.py` — 2 failures (legacy isotropic
  ansatz, fail in baseline too — pre-existing ERR-026-adjacent issue;
  per the plan, file is not xfail-marked but L1-tagged; left alone
  since baseline state).

**Convention discipline learned:**
- The `streaming_terms` extraction has TWO directional conventions
  baked into the API: face_area_in/out is canonicalised to
  outward-sweep semantics; abs_mu is computed from `mu_x[direction_idx]`
  which is correct only for slab and sphere (where `direction_idx` is
  the global ordinate). Cylindrical needs special handling at the
  caller.
- Wave C's DD test for inward sweep was bit-identical-to-itself but
  did NOT match the inlined sweep's inward convention — the reference
  formula in the test mirrors the DD's call to `streaming_terms`
  (which uses `face_area_out = A[i+1]`), but the inlined sweep's
  inward branch reads `A_out = A[i]`. The test passed because both
  the DD and the reference use the same (incorrect-for-inward) face
  semantics; integration with the actual sweep surfaces the
  discrepancy. The fix lives in the sweep wrapper, not in the DD test
  (the DD test pins the per-cell math; the directional convention is
  the sweep's responsibility).

**Performance note:** Sphere 2g homogeneous regression test ran in
1.07s (baseline) → 1.80s (post-rewrite). Cylinder LS4 1g homogeneous
ran in 1.47s → 2.98s. Sphere 2g 3reg ran in 4.0s. The ~2× slowdown
on curvilinear paths is expected from the per-cell Python-level
dispatch through the Protocol (vs inlined math); affects iteration-
heavy eigenvalue cases. This is acceptable for Wave D — performance
optimisation (e.g., a vectorised `update_batch`) is Wave C-extension
work after LD/EC/Step land. 1-D cumprod fast path preserved at full
speed for slab+GL+isotropic; 2-D wavefront also preserved.
