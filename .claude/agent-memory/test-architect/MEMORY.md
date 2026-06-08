# Test-Architect Memory Index

Read `lessons.md` at the START of every invocation. Durable verification
principles live in the `feedback_*` notes + the carve-lesson notes. The
committed suite records WHAT is verified and these notes keep WHY — and
several plan notes are themselves cited as the verification-rationale of
record by committed tests + `docs/theory/operator_algebra.rst`, so they are
load-bearing **provenance**, NOT redundant scaffolding. Do not delete
committed brain; hygiene the uncommitted working notes.

## Durable verification principles (read when designing a plan)

- [Test Architect lessons](lessons.md) — L1 homogeneous-blind-to-curvilinear; L2 walk the analytical-references checklist; L3 no mesh-independent transport eigenvalue reference (Issue #8).
- [Regression tolerance design](feedback_regression_tolerance_design.md) — `assert_regression`: iterative→`SAFETY(10)×conv_tol` read off run-config (SoT), direct→`nulp(reduction_depth)`; `DriftWarning` tripwire; `-O`-safe (no bare assert).
- [Eigen on non-fissile mixture is malformed](feedback_eigen_on_nonfissile_mixture.md) — k=prod/abs=0/0→nan dead gate; reformulate as fixed-source (same Pℓ/streaming path); corroborate vs flat inf-medium `(diagΣ_t−Σ_s0ᵀ)⁻¹Q`.
- [Diagnostic→test promotion](feedback_diagnostic_promotion.md) — verify-diag-runs-first; reproduce pins via public API; 3 foundation classes (invariant/negative-regression/math-origin); warning suppression; delete-after-pass.
- [V&V tagging idioms](feedback_vv_tagging.md) — module `pytestmark` vs per-test `verifies()`; foundation NO `verifies()`; xfail `strict=False`+`reason=` pointing at unlocking API/label.
- [Cross-method protocol design](feedback_cross_method_protocol.md) — reuse registry schema (thin `CrossMethodCase`); `max(tol_a,tol_b)` agreement; L1-not-L4 convention; verify truth values vs literature memos first.

## SN operator-algebra carve lessons (LANDED #208/Wave-O/Wave-T/Phase 4-5; tests live)

- [Issue #208 Wave-O carve lessons](issue_208_wave_o_carve_lessons.md) — snapshot-on-random-ψ + Q/Σt two-anchor; `--capture-baseline` ROOT-conftest gotcha; defect-vs-raw zero-input precision boundary; affine-torsor gate; sha256-golden bit-id; residual-before-adjoint; metric-blind reciprocity.
- [Wavefront-flux carve lessons](wavefront_flux_carve_lessons.md) — storage-A/B; inherited-octant-snapshot (interior edges ephemeral→outputs+boundary pinned); peak-pin Leg-1(size-ratio sharp) vs Leg-2(tracemalloc loose); nulp-einsum vs array_equal-typed-wrapper; SHED/CAPTURE highest-risk.
- [Phase 5a/5c angular-windowing lessons](phase5a_5c_angular_windowing_lessons.md) — THE TRAP (all 2-D snapshots P0-isotropic→a φ_ℓ≥1-dropping carve passes the whole suite); moment-tensor fuller-view oracle catches ℓ/m drift (Mode-5) scalar is blind to; nulp(n_diag) relaxation.
- [Phase 4.6 non-vacuum MMS ansatz](phase4_46_nonvacuum_mms_ansatz.md) — prior MMS all VANISH at boundary→q.boundary≠0 untested; `(A+μB)/W` with A=a0+a1·sin, **a0>0 load-bearing**; converged-VALUE-not-just-O(h²)-rate; sphere companion mandatory; G-4 reflective-mesh-not-kwarg dud.
- [SOTP separability verification](sotp_separability_verification.md) — master condition: SOTP needs Cartesian-product per-axis; coupled physics (per-material XS, M-M curvilinear recurrence) NOT SOTP-able→OperatorSum fallback; Route A array_equal vs Route B nulp; slab degenerate-separable hides curvilinear bugs.
- [SI convergence-rate verification](si_convergence_rate_verification.md) — rate≠value test-design; measurand `history.n_inner` (eigenvalue path =None, OPEN gap); target = analytic ρ_J=c; gate `<0.75×baseline`+value guards; 1G-OK for rate (flux-shape-indep by design).
- [Snapshot migration when production goes bare](snapshot_migration_when_production_goes_bare.md) — reusable recipe: shared-driver SoT; schema=persisted∩compared; vacuum-bit-id correctness gate; inheritance-from-bare needs closed-form anchor; retire false `@catches`; re-verify term activation.

## Test infrastructure & harness (LANDED)

- [SN sentinel (canary) harness](sn_sentinel_harness.md) — `@pytest.mark.sentinel` on the capability DAG; **vv Mode-8 `-O`-strips-bare-assert hazard** (gate runs WITHOUT -O); cosmic-ray mutation recipe + `git checkout` gotcha; thin-wrapper pattern; known gaps.
- [SN taxonomy reorg mapping](sn_taxonomy_reorg_mapping.md) — tests/sn capability-taxonomy reorg; file-destination map; cap-marker-via-dir-conftest; the xdist 3.8.0 + Py3.14.3 loadscope DEADLOCK (gate sequential per-dir, whole-tree -n auto unachievable).
- [Peierls rank-N / rank-2 closure](project_peierls_rank_n.md) — Marshak/DP_N rank-N + Phase F rank-2 per-face BC; regime-A bit-exact recovery (rtol 1e-14); §9.1 dual-route discipline (GL vs mpmath.quad of same integrand) for any new analytical formula; F.1/F.2 test layout.

## Landed-work verification plans (committed brain — KEEP; ★ = cited as provenance)

These pre-implementation plans/specs drove now-landed carves. The committed
tests are the live record of WHAT passed; these keep the gate design + WHY.
Several are cited by file-path from committed tests/docs — deleting them would
dangle that provenance, so they stay.

- [D-I.3 verification plan](D-I.3_verification_plan.md) — StreamingOperator bare-`ndarray` retirement (Depth B D-I.3): both 1-D B1''-packed + 2-D `ravel('F')` arms retire → `TimedFullField`-only.
- ★ [D-G BoundaryFlux pure-Field verification](dg_boundary_flux_pure_field_verification.md) — D-G BoundaryFlux pure-Field carve + SweepScratch split; pre-D-G snapshot + 3 specs + ULP-bound equivalence. Cited by `tests/transport/fields/test_boundary_flux.py` (+2 tests).
- [Issue #196 Phase G Step 1 gates](issue_196_phase_g_step1_verification_gates.md) — DiamondDifference→`SNCellOperator` + MorelMontry→`AngularRedistribution` type promotion; 7 gates + Phase-F twin-path defense.
- [Phase 1 verification plan (moment-space)](phase1_verification_plan.md) — ERR-039 fix via `SphericalHarmonicBasis/Space`/`MomentProjection` + `.T`-vs-`.H`; pin legacy + new-convention tests + L1 MMS gates.
- [Phase 3 verification plan (moment-space)](phase3_verification_plan.md) — layered re-architecture: P3.1 import-linter `FORBIDDEN_EDGES` + per-step gates P3.2–P3.6 + AngularFlux split.
- ★ [Wave T.3 scattering SOTP spec](wave_t_t3_scattering_verification_spec.md) — `ScatteringOperator.apply`→`SumOfTensorProductsOperator` pre-impl plan. Cited by `tests/sn/_fixtures/wave_t_t3/`.
- ★ [Wave T.4 streaming SOTP spec](wave_t_t4_streaming_verification_spec.md) — `StreamingOperator.apply`→`L_spatial+L_angular_redist` tensor-network pre-impl plan. Cited by `tests/sn/operators/test_streaming_operator.py` + `docs/theory/operator_algebra.rst`.
