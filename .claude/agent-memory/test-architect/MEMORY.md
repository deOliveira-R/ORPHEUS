# Test-Architect Memory Index

One line per entry — detail lives in the linked file, NEVER inlined here (the
index is loaded whole every dispatch; keep it small). Four sections: (1) lessons
— READ `lessons.md` FIRST every dispatch; (2) active/in-flight state — git-true
(reconcile "unmerged" claims against git before acting); (3) durable reference
recipes; (4) design idioms. The failure-mode taxonomy lives in `vv-principles`;
the reference inventory + XS mixtures in `AGENT.md`. No campaign play-by-play
here — it is merged archaeology.

## 1. Lessons — a HOT digest over a COLD archive (READ the digest at START)

- **[Lessons — hot digest](lessons.md)** — 1013 lines (`[M]` 2026-09-21). One imperative rule per
  entry, grouped by behavioral family (eight meta-lessons · gates that cannot red ·
  harness discipline · config blindness · reference & claim layer · tolerance ·
  carve archetypes · snapshots & exactness · pure-math primitives). **Read this
  file whole, every dispatch.** Every entry ends with a `→ LNN` pointer into the
  archive; families 6 (carve archetypes) and 8 (pure-math primitives) are
  reference lookup tables — read only their one-paragraph meta-rule per
  dispatch, open `lessons_archive.md` §L87/§L88 when a carve matches a shape.
- **[Lessons — cold archive](lessons_archive.md)** — ~10 900 lines, sections L1–L88,
  append-ordered. The war stories, measured numbers, `file:line` detail, the
  carve-archetype lookup table (§L87) and the pure-math primitive long form
  (§L88). **Open ONE section at a time, only when the digest's pointer says the
  detail matters.** Never read it whole — that is ~55K tokens.
- NO lesson content is inlined here. The digest is the index over the archive;
  this file is the index over everything else. New lessons: add the RULE to the
  digest (with its `→ LNN`) and the war story as a new archive section.

## 2. Active / in-flight verification work

**Detail → [active campaigns](active_campaigns.md) and `lessons_archive.md` §LNN.**
ONE line each here — name, terminal status, pointer. Merge status comes from git
and GitHub, never from this list (`process-discipline`); a landed campaign is
archaeology and lives only in the archive.

- **Consumers step 3 — the Solution carries its POSING** — plan + anchors delivered 2026-09-17, PRE-carve (`[M]` 2026-09-21: `tests/sn/architecture/test_step3_solution_anchors.py` still carries 7 strict-xfail rows; #484, the adjoint posing, follows it). → **`L86`**
- **#432 — the axis-parameterised O(2) member** — issue OPEN; the stabiliser gates shipped. → **`L70`**
- **#235 — the 2-D angular closure's ranking INSTRUMENT** — design delivered, issue OPEN. → **`L48`**
- **#358 — the test-dependence DAG** — memo delivered, issue OPEN (the main memory's TEST-DAG thread). → **`L55`**
- **Everything else this list carried is MERGED** (`[M]` 2026-09-21, `gh issue view`: #459, #448, #426, #434, #429, #325, #337, #2, #280, #340, #344, #290 all CLOSED; the CS ladder CS1–CS5 with its P4 remainder and the consumers steps 1–2 COMPLETE 2026-09-18 @ `c27373b9`; the boundary machinery, G2/G5/G6, the three-DOF separation and the prior SN campaigns before that). The record is the SN theory page's development history and the archive §L17–§L85; a "#41" once listed here was a plan-internal number, not issue #41.

## 3. Durable reference (reusable verification-design recipes)

Reusable RECIPEs / cited by `AGENT.md`. Core lessons in `lessons.md`; these keep the worked method.

- [Convergence-RATE verification](si_convergence_rate_verification.md) — AGENT.md §5. Iterations-to-converge vs analytic SI ρ=c; measurand `history.n_inner`; the OPEN eigenvalue-path `n_inner=None` gap; rate-claims flux-shape-independent → 1G-OK.
- [Snapshot migration when production goes BARE](snapshot_migration_when_production_goes_bare.md) — AGENT.md §7. Shared-driver SoT; schema=persisted∩compared; VACUUM-bit-id gate; snapshot-inheritance-needs-anchor; false-`@catches` retirement; term-activation re-verify.
- [SN sentinel harness](sn_sentinel_harness.md) — `@pytest.mark.sentinel` one-cheap-test-per-capability-node; cosmic-ray mutation-validation (`git checkout` after each run); per-NODE-sentinel-leaves-interior-uncovered gap.
- [SOTP separability verification](sotp_separability_verification.md) — separable ⟺ Cartesian-product per-axis; coupled physics → OperatorSum fallback; Route-A array_equal vs Route-B nulp; slab degenerate.
- [Operator space-guard only bites OperatorSum](operator_space_guard_only_bites_operatorsum.md) — the domain/codomain guard is INVISIBLE to SI/Krylov matvec; bites only actually-composed sums; `FunctionSpace.__eq__` by `(name,shape)`; activation-gate the composed sum.
- [Cross-layer relocation carve](cross_layer_relocation_carve_verification.md) — relocate-down + registry-dispatch. H1 registration-timing MASKED by process-global state → fresh-process subprocess gate mandatory. H2 `TYPE_CHECKING` sn import trips `test_layer_imports`. Layer-inversion usually doc-only at runtime.
- [A3/#280 reverse-scan transpose-solve](a3_reverse_scan_transpose_verification.md) — reverse-DAG `apply_transpose`; retired-CAP→typed-predicate reconciliation; assembled-Mᵀ (Cartesian-only) vs dense-apply SPHERE keystone; 1-D loop spy + orientation-OBJECT AST tripwire. §7 CYLINDER arm: mandatory `product(n_mu=4,n_phi=8)` (LS nulls both hard terms=control); G1/G2-dense-Mᵀ-keystone/G3-full-field-recip(#284)/G4/G5; ERR-066 degenerate-drop tooth.
- [A_BA ψ½ Schur-fold un-weld](aba_schur_fold_unweld_verification.md) — lessons L22. Welded-fold un-weld (N sites→ONE source). 7 gate types: manufactured-anisotropic fold contract, Mode-11 wrap-counter EXACT `2·n_levels`, bit-id INHERITS + independent `½·emission`, two transpose gates, F-non-vacuity, cyl/slab None-ray control.
- [CoupledOperator Step-4 verification](coupled_operator_step4_verification.md) — N-general block machinery (ψ½=instance #1). 4d.0 `FullField`→`System[I,B]` structure-only (multi-instantiation synthetic CRUX). 4d.1 assemble≡probe principled-equiv + block-`.H` Mode-12. 4d.2 presence=block-existence. 4d.3 block-apply WRAPs fused walk + two-anchor.
- [CoupledOperator B.2b re-type](coupled_operator_b2b_retype_verification.md) — pure re-labeling → `array_equal` EVERY row (any rtol/nulp = RED FLAG). b1 split SourceSink + role-preserving bridge (role⊕values split-blind). b2 family-blind `from_blocks` + presence-dispatch. b3 A_BA/B_b onto ray composite + adapter-delegation sentinel.
- [A_AB seed-injection](a_ab_seed_injection_verification.md) — cell-local rectangular coupling (ray→bulk, σ-indep) = `A_bs` block. Equivalence gates SHARE closure methods → blind; the ONE catcher = gate-3 Euclidean fwd↔transpose adjoint-consistency + `test_radial_characteristic_metric` anchor. Sphere ONE level → multi-level untestable.
- [A_BB forward shared-kernel EXTRACT](radial_characteristic_forward_extract_verification.md) — Step 4b. Round-trip PRINCIPLED-EQUIV ~3 ULP not 0-ULP; `solve∘apply=id` only on CONSISTENT subspace; transpose seam adds seed_cells_bar A_AB term; EUCLIDEAN not V_cell metric; Mode-11 anti-twin routing sentinel.

## 4. Durable design idioms (feedback)

- [Regression tolerance design](feedback_regression_tolerance_design.md) — iterative→`SAFETY(10)×conv_tol` off run-config SoT, direct→`nulp(reduction_depth)`; `DriftWarning` tripwire; `-O`-safe.
- [Eigen on non-fissile mixture is malformed](feedback_eigen_on_nonfissile_mixture.md) — k=0/abs→nan dead gate; reformulate fixed-source; corroborate vs `(diagΣ_t−Σ_s0ᵀ)⁻¹Q`.
- [Diagnostic→test promotion](feedback_diagnostic_promotion.md) — verify-diag-runs-first; reproduce via public API; 3 foundation classes; delete-after-pass. (SoT: `tests/derivations/_promotion_policy.md`.)
- [V&V tagging idioms](feedback_vv_tagging.md) — module `pytestmark` vs per-test `verifies()`; foundation carries NO `verifies()`; xfail `strict=False`+`reason=`.
- [Cross-method protocol design](feedback_cross_method_protocol.md) — reuse registry schema; `max(tol_a,tol_b)` agreement; L1-not-L4; verify truth values vs literature memos first.
