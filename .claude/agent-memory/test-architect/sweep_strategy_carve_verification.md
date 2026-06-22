---
name: sweep-strategy-carve-verification
description: SweepStrategy carve (C3.4+C3.5) V&V plan. Polymorphic sweep/matvec dispatch (CumprodScan/FullFieldWavefront/MovingFrontierWindow). reduced⟺is_1d coincidence pin; A2D-1-wrap-not-relocate; Mode-8 dispatch-pin hazard; 2-D adjoint deferral gate; synthetic-d3 frontier_dim=d-1 admission.
metadata:
  type: project
---

S0 verification plan for the SweepStrategy carve (worktree `worktree-sn-nd-layout`,
HEAD `0f2ab34`). Plan = `.claude/plans/sn_sweep_strategy.md`. Replaces TWO scattered
predicates with ONE polymorphic strategy: sweep keys on `reduced is not None`
(`sweep.py:222`); 5 matvec gates key on `not is_1d` (`operator.py` lines 372/617/829/
1458/1528). Strategies: `CumprodScan`(1-D any geom), `FullFieldWavefront`(Cartesian any-d,
ORACLE), `MovingFrontierWindow`(Cartesian d≥2). NONE landed yet (pre-impl).

## ⭐ Structural findings (VERIFIED this session, load-bearing)

- **`reduced is not None` ⟺ `is_1d` for ALL 4 constructible meshes** (slab/sphere/cyl
  T/T/T; cart2d F/F). The two scattered predicates DO coincide today → the carve is
  behavior-preserving. PIN this (`test_reduced_iff_is_1d_coincidence`).
- **`is_cartesian` (= `curvature is None`) is ORTHOGONAL**: slab+cart2d both cartesian
  but route differently → strategy `supports`/`default_for` MUST key on BOTH
  `is_cartesian` AND `ndim`, NOT is_cartesian alone. (slab=cart `curvature=None`;
  sphere/cyl `'spherical'`/`'cylindrical'`.) Guard against a "select on is_cartesian
  only" regression.
- **Cylinder REQUIRES level-structured quadrature** (`level_symmetric`/`product`); GL
  raises `ValueError` (no LevelStructure side-channel). Use LS-4 for cyl/cart2d configs.
- **`from_cartesian(shape, label=OctantLabel(signs))` IS LANDED (C3.1).** The cumprod
  stub's `_SpineNotLanded`/`_wavefront_1d_sweep` adapter assumption is STALE — spine
  exists; only the d=1 *driver* was unwired. S3 drops the adapter, wires the real
  `CumprodScan.sweep`/`FullFieldWavefront.sweep` API.

## The anchor set (all GREEN @ 0f2ab34, batch = `47 passed, 3 xfailed` -O)

- **L2 PRIMARY S1-S2 tripwire** `tests/sn/solve/test_affine_carve_bit_identity.py` —
  sha256 golden of converged psi/phi via `solve_sn_fixed_source` end-to-end (routes
  sweep+matvec); 3 cases (2D-SI-windowed/2D-Krylov/1D-slab-SI); `-O`-safe `raise
  AssertionError`; MUST stay byte-identical.
- **L3 A2D-1** `test_streaming_operator.py::TestT4dApply2DCartesianSourceHashPin::
  test_apply_2d_cartesian_source_hash_unchanged` — sha256 of `inspect.getsource(
  _apply_2d_cartesian)`. ⚠ **S2 must WRAP not RELOCATE** `_apply_2d_cartesian` so the
  hash stays free-green; if relocated, regenerate hash + history line. Plan
  UNDER-SPECIFIES this — flag to implementer.
- **L4 DD-regression** non-square 8×4 2D cases catch x↔y swap (`assert_regression`
  SAFETY×conv_tol≈1e-11; 2d-aniso pre-drifts ~6920 ULP by design).
- **L5/L6 window≡full** `cartesian_2d/test_2d_full_field_oracle.py` +
  `core/test_sweep_graph_window_equivalence.py` (`assert_array_equal` bit-id) — FOLD
  into S3 d=2 strategy-vs-strategy; `_sweep_2d_full_field`+`_apply_2d_cartesian_full_field`
  RETIRE→wrapped by FullFieldWavefront (fuller-view-oracle exception, now selectable).
- **L7 φ=Q/Σ_t** `test_2d_octant_sweep_equivalence.py::*closed_form_anchor` (2×2
  numpy.linalg.solve) = d=2 converged-VALUE structural ground.
- **L8 k_inf=1.875** `test_wavefront_cumprod_equivalence.py::test_cumprod_path_hits_
  analytical_kinf` (transfer-matrix) = d=1 structural ground. Both L7/L8 STAY (vv §1.5:
  ULP-distance necessary-never-sufficient).
- **L9 hand-loop** `test_sweep_graph.py::TestApplyMatchesLegacyInlined::
  test_per_cell_loop_equivalence` non-square (3,4)/(4,3) = the x↔y-swap moat.

## ⚠ Mode-8 LIVE HAZARD (demonstrated)

`tests/sn/sweep/core/test_unified_sweep_dispatch.py` (the legacy-dispatch pin,
`TestDispatchByReducedProperty`) uses bare `assert` — passes under `-O` (7 passed) WHILE
asserting NOTHING (pytest warns "assertions ignored"). S1 MIGRATES it to strategy-
selection + converts bare assert → `np.testing`/`pytest.fail`. Sibling
`test_sweep_graph.py::TestAssertCellCoverage` also bare-assert (run WITHOUT -O or note
inert).

## Deferral boundary (RISK, now gated)

2-D Cartesian ADJOINT raises NotImplementedError TODAY (`apply_transpose` op.py:1528;
`_compute_LpC_transpose` op.py:617). S2 strategy MUST preserve → `IncompatibleStrategy`/
not-applicable, NEVER silent wrong answer. NEW gate `test_2d_cartesian_adjoint_raises_
incompatible`.

## S4 synthetic-d3 idiom (model)

`tests/sn/sweep/core/test_sweep_graph_nd_admission.py::_walk_roundtrip_residual` +
`test_d3_walk_residual_vanishes_at_apply_solution` — apply↔residual round-trip on a
synthetic `d`-tuple (NO real quadrature). S4's `MovingFrontierWindow` frontier_dim=d−1
`window≡full` mirrors it: d=1(frontier_dim=0 base)/d=2(bit-id)/synthetic-d3
(`assert_array_equal`). d=3 CONTIGUITY/speedup is OUT of correctness gate (profiling).

## Deselect / standing reds
Deselect `tests/sn/eigenvalue/test_keff_slab.py::test_heterogeneous_absolute_keff`
(#212 hang). Orthogonal held reds: #206 cyl-matvec (deferred adjoint), #195 MMS@160 —
both untouched by the 2-D-gated carve.
