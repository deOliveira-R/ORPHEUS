# #340 — the solves whose INNER exhausts its budget while the outer reports success

`[M]` 2026-08-10, HEAD `52650a86`. Instrument:
`$CLAUDE_JOB_DIR/tmp/count_silent_truncations.py` — a pytest plugin that wraps
`_warn_if_unconverged`, records every call where `converged and not fully_converged`,
and then delegates to the original UNCHANGED. Nothing warned that did not warn before.

Command (SERIAL, canonical `-O`):
```
PYTHONPATH=$CLAUDE_JOB_DIR/tmp .venv/bin/python -O -m pytest \
  tests/sn tests/numerics tests/cp tests/moc tests/diffusion \
  -m "not slow" -p count_silent_truncations -q
```

Run outcome: **9 failed / 5501 passed / 1 skipped / 130 deselected / 61 xfailed**
in 1492.92 s — the 9 are the known task-#51 baseline reds, so the instrument
perturbed nothing.

**25 calls across 20 tests.** EVERY one is `inner(source-iteration)`;
zero `inner(gmres)`, zero CP/MoC/diffusion (those three carry no record yet — N4).

These do NOT warn today and must NOT start warning in N6 commit 1 (the guard
stays on `converged`). They are the population commit 2 will interrogate with
N5's residual certificate — user ruling 2026-08-10: *certify them, then decide*.

⚠ `test_inner_tol_bias_collapses_at_1e_12` is a **deliberate** truncation — the
test exists to study inner-tolerance bias, so its starved inner is the fixture,
not a defect. That makes **R2** ('declare the deliberate truncations in-test')
and N6 the SAME population viewed twice: R2 declares the intentional ones so
commit 2 can stay silent about them on purpose rather than by accident.

| test | calls | level | ran/budget | rho |
|---|---|---|---|---|
| `tests/sn/eigenvalue/test_keff_2d.py::TestSIKrylov2DEquivalence::test_eigenvalue_jacobi_gauss_seidel_equivalence` | 2 | `inner(source-iteration)` | 500/500 | 0.9912 |
| `tests/sn/eigenvalue/test_keff_curvilinear.py::TestCylinderMultiGroupMultiRegion::test_2g_heterogeneous_fuel_moderator` | 1 | `inner(source-iteration)` | 500/500 | 0.9729 |
| `tests/sn/eigenvalue/test_keff_curvilinear.py::TestCylinderMultiGroupMultiRegion::test_2g_heterogeneous_product_different_resolutions` | 2 | `inner(source-iteration)` | 500/500 | 0.9729 |
| `tests/sn/eigenvalue/test_keff_curvilinear.py::TestCylinderMultiGroupMultiRegion::test_multigroup_eigenvector_not_flat` | 1 | `inner(source-iteration)` | 500/500 | 0.9729 |
| `tests/sn/eigenvalue/test_keff_curvilinear.py::TestMultiGroupMultiRegionSpherical::test_2g_heterogeneous_converges` | 1 | `inner(source-iteration)` | 500/500 | 0.9788 |
| `tests/sn/eigenvalue/test_keff_curvilinear.py::TestMultiGroupMultiRegionSpherical::test_multigroup_eigenvector_not_flat` | 1 | `inner(source-iteration)` | 500/500 | 0.9788 |
| `tests/sn/regression/test_dd_regression.py::test_dd_regression[2d_2g_LS4_dd_8x4_het_si]` | 1 | `inner(source-iteration)` | 300/300 | 0.9705 |
| `tests/sn/regression/test_dd_regression.py::test_dd_regression[cyl_2g_3reg_folded_4x8_dd_n40]` | 1 | `inner(source-iteration)` | 300/300 | 0.9578 |
| `tests/sn/regression/test_dd_regression.py::test_dd_regression[slab_2g_3reg_dd_n40]` | 1 | `inner(source-iteration)` | 300/300 | 0.9577 |
| `tests/sn/regression/test_dd_regression.py::test_dd_regression[sphere_2g_3reg_dd_n40]` | 1 | `inner(source-iteration)` | 300/300 | 0.9560 |
| `tests/sn/regression/test_dd_regression.py::test_dd_regression[sphere_2g_homogeneous_dd_n20]` | 1 | `inner(source-iteration)` | 300/300 | 0.9281 |
| `tests/sn/solve/test_sn_adjoint_certification.py::TestP13KEquality::test_heterogeneous_slab_k_equality` | 1 | `inner(source-iteration)` | 500/500 | 0.9558 |
| `tests/sn/solve/test_sn_adjoint_certification.py::TestP13Mutations::test_streaming_no_reversal_shifts_k_heterogeneous` | 1 | `inner(source-iteration)` | 500/500 | 0.9558 |
| `tests/sn/sweep/core/test_cache.py::test_collision_cache_invariance_under_source_iteration` | 1 | `inner(source-iteration)` | 50/50 | 0.9581 |
| `tests/sn/sweep/curvilinear/test_w1_clamp_silent_on_flat.py::test_homogeneous_reflective_sphere_iso_unchanged` | 1 | `inner(source-iteration)` | 300/300 | 0.9281 |
| `tests/sn/sweep/curvilinear/test_w1_clamp_silent_on_flat.py::test_unclamped_sphere_flux_strictly_positive[2g_R2_S64]` | 1 | `inner(source-iteration)` | 800/800 | 0.9927 |
| `tests/sn/verification/analytical/test_kinf_homogeneous_tolerance.py::test_inner_tol_bias_collapses_at_1e_12[2eg-cylinder]` | 2 | `inner(source-iteration)` | 100/100 | 0.9228 |
| `tests/sn/verification/analytical/test_kinf_homogeneous_tolerance.py::test_inner_tol_bias_collapses_at_1e_12[2eg-sphere]` | 2 | `inner(source-iteration)` | 100/100 | 0.9292 |
| `tests/sn/verification/analytical/test_kinf_homogeneous_tolerance.py::test_inner_tol_bias_collapses_at_1e_12[4eg-cylinder]` | 1 | `inner(source-iteration)` | 100/100 | 0.8892 |
| `tests/sn/verification/analytical/test_kinf_homogeneous_tolerance.py::test_inner_tol_bias_collapses_at_1e_12[4eg-sphere]` | 2 | `inner(source-iteration)` | 100/100 | 0.9047 |

**25 calls across 20 tests.**
