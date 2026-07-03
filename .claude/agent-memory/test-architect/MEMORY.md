# Test-Architect Memory Index

Three sections: (1) the lessons pointer — READ FIRST every dispatch;
(2) active/in-flight state — git-true; (3) durable reference (reusable
verification-design files cited by AGENT.md or holding a meta-recipe).
The failure-mode taxonomy lives in the `vv-principles` skill; the
reference inventory + XS mixtures live in `AGENT.md`. Do NOT re-add
campaign play-by-play here — it is merged archaeology.

## 1. Lessons (durable verification discipline — read at START)

- [Test Architect lessons](lessons.md) — L1 config-blindness (flat/1G/homog/slab/iso-snapshot); L2 reference↔claim-layer + structural independence; L3 MMS simplification-bias override + activate/null declaration; L4 prove every gate's teeth bite (Mode-8/xfail/mutation/call-count); L5 characterization-gate without calcifying; L6 Mode-10 structural teeth; L7 regression-tolerance = the claim (SAFETY×conv_tol, non-fissile→fixed-source, IEEE-vs-exact); L8 cross-method infra; L9 V&V tagging; L10 proactive-refute-the-premise payoff; L11 axis-transpose-of-shipped-reduction (gate the FLIPPING rule, not the shared machinery); L12 fast-path→composed FOLD = principled-equiv-not-0-ULP + unchanged-sibling-stays-0-ULP + "transpose-falls-out-free" gated on the ONE leaf transpose + Euclidean-not-`.H` reciprocity + Mode-11 leaf-wrap + outside-the-fold channel; L13 capability-string→typed-operator carve (`.solve`/`CAP_SOLVE`→`.inverse()`, `CAP_APPLY_TRANSPOSE`→`.H`): equiv-keystone=bit-id INHERITANCE+closed-form anchor+seed Mode-11 wrap; the WHOLE frozenset (BOTH axes) retires SAFELY (the "iff both" law is in the METHOD BODY not the set; half-modernized axis=twin-path HAZARD) — licensed by a FAITHFULNESS scaffold gate (derived `is_invertible`/`is_adjointable` ≡ old frozenset ∀ operator, asymmetry fixtures mandatory, then deleted; recursive pins permanent); runtime-query=derived `@property` NOT `isinstance(Supports*)` (class-uniform→blind to half-adjointable composites+value-dep singular edge; Protocols=static contract only); LITERAL-string `"apply_transpose"` precondition in 2 S† canaries breaks on deletion (grep both); M-ADJ-FORGE downstream=RAISE not wrong-value; adjoint-of-inverse via RECIPROCITY not round-trip (metric≠Euclidean activates curvilinear). L14 the FIRST *iterative* inverse (`GreenOperator`=`OperatorSum.inverse()` wrapping the SI driver): NO legacy-solve to inherit→keystone=dense-LU-anchor+G-Neumann NOT bit-id, `inverse().apply≡solve` is a TAUTOLOGY (solve:=inverse().apply)→not evidence; `ConvergenceFailure` MUST check TRUE residual `‖Aψ−q‖/‖q‖` NOT the driver INCREMENT (Signature-9 ρ-blind, c=0.99 falsifier); `is_invertible=a.is_invertible` FLIPS frozen `is_invertible is False` pins False→True by design (grep-audit WITHOUT adjoint-filter which HIDES compound-line hits; leaf/product STAY, plain-sum-invertible-leading FLIPS; "no general (A+B)⁻¹" retires) + faithfulness green ONLY if conditional CAP_SOLVE is LOCKSTEP + BOTH asymmetry fixtures (inv+ao & ao+inv; (T,F) uniquely IDs the a-rule); 3 ordering edges (L+C→Sweep via MRO / C+L→divergent-Green-loud-fail / (−S)+A→refuse-at-ctor-canonical-msg) + flatten-EXACT-OperatorSum-only; NEW iterative op needs its OWN Mode-11 seed pin (landed per-iterate spy one level below is BLIND + drop value-invisible); name-earning triage (G-Neumann+dense-anchor LOAD-BEARING / G-reciprocity INCLUDE-lighter non-symmetric-A Euclidean-not-`.H` / G-kernel FOLD-into-anchor; L1 Mode-9=het≥2G VACUUM slab); wrap-delegate MIXIN extraction is no-op⇒existing-suites-green-UNCHANGED + pyright abstract-sig teeth (#285 STRUCTURAL). L15 the RETIREMENT terminal step (frozenset→typed-predicate, EXECUTES L13's design): the coexistence scaffold (`is_X≡CAP_X∈caps`) DELETES→its PERMANENT successor is keystone-v2, a structural/value CONTRACT (`is_invertible False`⟹`_STRUCTURAL_ABSENT` not-hasattr OR `_VALUE_RAISE` eager-`NotInvertible`; pin WHICH per type) + eager-`.H`-or-`MissingAdjoint` mirror + bridge-consistency (`invertible(op)==is_invertible` pins the TypeGuard body); migration is a MECHANICAL-RULE + completeness-REGREP (returns-ZERO) NOT a table (surface WIDER than any map—127 reads/33 files, re-grep yourself; bare-string `"apply_transpose"` preconditions missed by CAP-grep); land-order W1-additive/W2-ATOMIC/W3-teeth (frozenset+its exception-class retire together—a W1 `MissingCapability→NotInvertible` flip breaks `pytest.raises(MissingCapability)` mid-flight; eager-`.H` is the ONE W1 behavior change; ratchet re-baselines DOWN); the pruned-body `OperatorProduct.solve` re-route = bit-id-INHERITANCE(`--capture-baseline`)+dense-anchor+Mode-11-sentinel, DIRECT-factors `array_equal` but the ITERATIVE Green-factor row `assert_regression(SAFETY×driver_tol)` NOT array_equal (L7 iterated-drift), non-commuting Diag@3-cycle-Perm activating.

## 2. Active / in-flight verification work

- **Diffusion integration #290 (pre-carve plan DELIVERED @ `d2a2a0c`).**
  [Verification plan](diffusion_integration_290_verification_plan.md) — the 7-deliverable spec for
  replacing `orpheus/diffusion/` with an in-algebra operator on MaterialMesh+ScalarFlux. CRUX
  (arithmetic-verified): a `SigS1=0, SigT=transport` Mixture reproduces legacy D BIT-IDENTICALLY →
  references DON'T move → all L1 anchors survive. 3 USER RULINGS: vacuum-BC semantics
  (keep zero-flux Dirichlet vs Marshak #182), data-seam bit-id vs re-baseline, Richardson retirement
  (H4 self-referencing → DELETE). Mode-12: k=λ_max(A⁻¹F) is similarity+transpose blind → object gate
  `A_diff.as_matrix()≡FD-stencil` + MMS flux-shape (#93, override the single-group/constant-D trap).
  Awaiting carve.
- **None (SN).** Every SN campaign whose closeout this memory once tracked
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
- [Operator space-guard only bites OperatorSum](operator_space_guard_only_bites_operatorsum.md) — the `OperatorSum`/`Product` domain/codomain guard is INVISIBLE to the SI/Krylov matvec (`L.apply−Σgᵢ.apply`, no OperatorSum; drivers validate only capabilities); it bites only on actually-composed sums (`InvertibleOperator`, `evaluate_residual`'s `LC-S-B`). `FunctionSpace.__eq__` by `(name,shape)`; activation-gate targets the composed sum, blocker = a new space name. (Banked from the P4.5b/c verification plan.)
- [Cross-layer relocation carve verification](cross_layer_relocation_carve_verification.md) — verifying "relocate a method-agnostic unit DOWN a layer + route hardcoded dispatch through a registry". H1 registration-timing regression MASKED by process-global registry state → fresh-process subprocess gate is mandatory (same-process = vacuous, Mode-11-adjacent). H2 verbatim move TRIPS `test_layer_imports.py` via a `TYPE_CHECKING` sn import (INPUT pkgs are NOT TC-exempt) → drop it, retype param `Any`. The "layer inversion" is usually documentation-only at runtime (`sys.modules` diff empty) — refute at design time. (Banked from the `realize_recursively` move spec.)

## 4. Durable design idioms (feedback)

- [Regression tolerance design](feedback_regression_tolerance_design.md) — `assert_regression`: iterative→`SAFETY(10)×conv_tol` off the run-config SoT, direct→`nulp(reduction_depth)`; `DriftWarning` tripwire; `-O`-safe.
- [Eigen on non-fissile mixture is malformed](feedback_eigen_on_nonfissile_mixture.md) — k=0/abs→nan dead gate; reformulate as fixed-source; corroborate vs flat inf-medium `(diagΣ_t−Σ_s0ᵀ)⁻¹Q`.
- [Diagnostic→test promotion](feedback_diagnostic_promotion.md) — verify-diag-runs-first; reproduce pins via public API; 3 foundation classes; warning suppression; delete-after-pass. (Policy SoT: `tests/derivations/_promotion_policy.md`.)
- [V&V tagging idioms](feedback_vv_tagging.md) — module `pytestmark` vs per-test `verifies()`; foundation carries NO `verifies()`; xfail `strict=False`+`reason=`.
- [Cross-method protocol design](feedback_cross_method_protocol.md) — reuse registry schema (thin wrapper); `max(tol_a,tol_b)` agreement; L1-not-L4 convention; verify truth values vs literature memos first.
