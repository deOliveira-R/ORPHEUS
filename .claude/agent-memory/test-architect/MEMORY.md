# Test-Architect Memory Index

Three sections: (1) the lessons pointer — READ FIRST every dispatch;
(2) active/in-flight state — git-true; (3) durable reference (reusable
verification-design files cited by AGENT.md or holding a meta-recipe).
The failure-mode taxonomy lives in the `vv-principles` skill; the
reference inventory + XS mixtures live in `AGENT.md`. Do NOT re-add
campaign play-by-play here — it is merged archaeology.

## 1. Lessons (durable verification discipline — read at START)

- [Test Architect lessons](lessons.md) — L1 config-blindness (flat/1G/homog/slab/iso-snapshot); L2 reference↔claim-layer + structural independence; L3 MMS simplification-bias override + activate/null declaration; L4 prove every gate's teeth bite (Mode-8/xfail/mutation/call-count); L5 characterization-gate without calcifying; L6 Mode-10 structural teeth; L7 regression-tolerance = the claim (SAFETY×conv_tol, non-fissile→fixed-source, IEEE-vs-exact); L8 cross-method infra; L9 V&V tagging; L10 proactive-refute-the-premise payoff.

## 2. Active / in-flight verification work

- **None.** Every SN campaign whose closeout this memory once tracked
  (#206/#208/#236/#240/#247/#251/#257/#18/#19/#20) is MERGED to `main`
  (verified: #257 branch deleted, its hashes ancestors of HEAD; the #236
  separability gate `test_space_angle_separability.py` is on main). The
  `feature/sn-spatial-angular-product` branch is 30 commits ahead with
  its own "#236 complete" closeout — finished work awaiting an ff-merge,
  NOT open verification work. Before resuming any branch claim,
  re-verify with `git merge-base --is-ancestor <hash> HEAD`.

## 3. Durable reference (reusable verification-design recipes)

These hold a reusable RECIPE or are cited by name from `AGENT.md` — not
campaign logs. Their core lessons are in `lessons.md`; these files keep
the full worked method.

- [Convergence-RATE verification](si_convergence_rate_verification.md) — cited by AGENT.md §5. Verify an iterations-to-converge / acceleration claim vs the analytic SI spectral radius ρ=c; measurand `history.n_inner`; the still-OPEN eigenvalue-path `n_inner=None` API gap; rate-claims are flux-shape-independent → 1G-OK (declare the claim layer).
- [Snapshot migration when production goes BARE](snapshot_migration_when_production_goes_bare.md) — cited by AGENT.md §7. Shared-driver SoT, schema=persisted∩compared, VACUUM-bit-id correctness gate, snapshot-inheritance-needs-anchor, false-`@catches` retirement, term-activation re-verify after an inject.
- [SN sentinel (canary) harness](sn_sentinel_harness.md) — `@pytest.mark.sentinel` one-cheap-test-per-capability-node design; cosmic-ray mutation-validation recipe (`local` distributor, `git checkout` after every run); the per-NODE-sentinel-leaves-module-interior-uncovered gap-closer.
- [SOTP separability verification](sotp_separability_verification.md) — the SOTP master condition (separable ⟺ Cartesian-product per-axis; coupled physics → OperatorSum fallback); how to gate an `assert_separable` claim; Route-A array_equal vs Route-B nulp bound; slab is the degenerate case.

## 4. Durable design idioms (feedback)

- [Regression tolerance design](feedback_regression_tolerance_design.md) — `assert_regression`: iterative→`SAFETY(10)×conv_tol` off the run-config SoT, direct→`nulp(reduction_depth)`; `DriftWarning` tripwire; `-O`-safe.
- [Eigen on non-fissile mixture is malformed](feedback_eigen_on_nonfissile_mixture.md) — k=0/abs→nan dead gate; reformulate as fixed-source; corroborate vs flat inf-medium `(diagΣ_t−Σ_s0ᵀ)⁻¹Q`.
- [Diagnostic→test promotion](feedback_diagnostic_promotion.md) — verify-diag-runs-first; reproduce pins via public API; 3 foundation classes; warning suppression; delete-after-pass. (Policy SoT: `tests/derivations/_promotion_policy.md`.)
- [V&V tagging idioms](feedback_vv_tagging.md) — module `pytestmark` vs per-test `verifies()`; foundation carries NO `verifies()`; xfail `strict=False`+`reason=`.
- [Cross-method protocol design](feedback_cross_method_protocol.md) — reuse registry schema (thin wrapper); `max(tol_a,tol_b)` agreement; L1-not-L4 convention; verify truth values vs literature memos first.
