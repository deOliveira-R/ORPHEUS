# Test-Architect Memory Index

One line per entry — detail lives in the linked file, NEVER inlined here (the
index is loaded whole every dispatch; keep it small). Four sections: (1) lessons
— READ `lessons.md` FIRST every dispatch; (2) active/in-flight state — git-true
(reconcile "unmerged" claims against git before acting); (3) durable reference
recipes; (4) design idioms. The failure-mode taxonomy lives in `vv-principles`;
the reference inventory in `AGENT.md` §2. No campaign play-by-play
here — it is merged archaeology.

## 1. Lessons

- **[Lessons — hot digest](lessons.md)** — 533 lines. The entries no rule, skill or definition clause carries (the workflows rule, invariant 6), each with a `→ LNN` pointer into the archive. Read it whole, every dispatch. Pruned 2026-09-22 by the agent-definitions audit: 73 entries restated a clause or are now carried by the definition.
- **[Lessons — cold archive](lessons_archive.md)** — 11573 lines, sections L1–L97, append-ordered: the war stories, the carve-archetype lookup (§L87) and the pure-math primitive long form (§L88). Open one section at a time, when the digest's pointer says the detail matters.

## 2. Active verification work

Merge status comes from git and GitHub, never from this list (`process-discipline`).

- **#432 — the axis-parameterised O(2) member** — issue OPEN; the stabiliser gates shipped. → **`L70`**
- **#235 — the 2-D angular closure's ranking instrument** — design delivered, issue OPEN. → **`L48`**
- **#358 — the test-dependence DAG** — memo delivered, issue OPEN. → **`L55`**
- **#405 P1 steps 3-8 (W1/W3)** — specs `.claude/plans/reference_p1_spec.md` §1.3a-§1.8; batteries `scratch/reference_architecture/p1step{3,5,6,7,8}/`. Step 8 landed `eaa74163`. On resume: confirm landed; re-run each step's owed arms on the real modules. → `L95`, `L99`-`L101`
- **Reflective cleanup (W3, `refactor/reflective-is-a-mirror`)** — gates delivered 2026-10-01: spec `scratch/boundary_ontology/reflective_cleanup_gates.md`; carry fixture + batteries `scratch/boundary_ontology/reflective_gates/`. On resume: re-run the carry gate post-carve.
- **Platform drift W2 (`fix/platform-independent-quadrature`)** — gates (a) CR Gauss rules, (b) GL fingerprint, (c) exact k∞ + bound, spec `scratch/platform_drift/gates_spec.md` (2026-10-01). On resume: confirm landed, re-run `cr_sim`/`r1_sim` arms. → **`L98`**
- **#405 P2, P3 and follow-ups (W3)** — specs `.claude/plans/reference_p2_spec.md`, `reference_p3_spec.md`; batteries `scratch/reference_architecture/p2/ta/`, `p3/ta/`, `p3/gates_repair/`. On resume: confirm landed; build the 7b.2 battery; re-run the P3 arms on the real module. → `L102`-`L106`
- **#200 Krylov gates (W3)** — evidence `scratch/reference_architecture/p3/krylov200/`, `gates200b/`. On resume: confirm landed; tag the record rows with the new ERR id. → `L107`
- **Characteristic-reference campaign (seed, P1 steps a–e, P2, P3 steps)** — CLOSED 2026-10-18 (`8ffafc6a`); each step's spec + battery under `scratch/characteristic_architecture/`; lessons `L108`–`L115`.
- **Under-resolved adjoint rung (W2 follow-up, `docs/frame-square-closure`)** — 2026-10-18, uncommitted: `TestAdjointAtAnUnderResolvedOrder` (LS4 L=2/3/5) + frame edge row; battery `scratch/adjoint_rung/arms.py` (`ARM=drop|mixed|short`). On resume: confirm landed.
- **Adjointability honesty gate (W2, ruling 27 of `numerics_untangling.md`)** — 2026-10-18, uncommitted on `main`: `tests/gates/numerics/test_adjointability_is_honest.py` (leg 1 flag⟺.H applies, leg 2 witness routes); battery `scratch/adjoint_honesty/battery/run.sh <arms>`. On resume: re-run after the fixes; 17 reds expected to clear.
- **Numerics untangling S1-S7 (W3)** — specs `scratch/numerics_untangling/packages/S*_gate_spec.md`, batteries `S*_gates/`. S7 landed `47e7bfaa`/`ac6c6558`. On resume: re-run each step's arms and bit baselines on the landed tree. → `L116`-`L118`
- **Wave 2 S2b/S3/S4 boundary carve gates (W3)** — prototypes + batteries `scratch/wave2_s2b/`, `wave2_s3/`, `wave2_s4/`. On resume: confirm landed; re-run the arms on the real realizer tables. → `L119`
- **Sub-cone measurements r2/r3 (rulings 23, 30)** — `scratch/verification_tree/subcones_r2/report.md`, `scratch/layer_dag/refiling/subcones_r3/`. On resume: only on the orchestrator's rulings.
- **#423 moment basis S0 (W3, `refactor/moment-basis`)** — spec 2026-10-25 `scratch/layer_dag/moment_basis/gate_spec.md`; drafts, plugins, batteries `.../ta/`. On resume: confirm the 3 drafts landed; re-run `ta/run_sites.sh` arms on the carved tree. → **`L120`**
- **Layer DAG 1b, Phase 2, Phase 4 (W3)** — specs `scratch/layer_dag/refiling/specs_1b/*.md`, `specs_2r/numerics_laws.md`, `scratch/layer_dag/phase4/spec.md`. On resume: main's rulings on each spec's NEEDS; confirm landed. → `L121`, `L124`-`L126`
- **Space membership U0-U3b (W3)** — specs `scratch/membership_covector/{spec,u0,u1b,u3b}/`; merged `1c1d8d1c`. On resume: only the owed multi-D arms (#616). → `L122`
- **Cone store B1-B5 (W1/W3)** — specs `scratch/verification_tree/cone_store/{spec_b1.md,b2/,b3/,b4/,b5/}`. On resume: the owed post-capture arms per spec. → `L127`-`L131`
- **#683 matrix carriers act on members (W3, `fix/683-assembled-operator-membership`)** — 2026-11-02: gates `tests/gates/numerics/test_matrix_carrier_membership.py` (17 rows) + 3 re-posed rows + Phase 2 fixture moved to Energy ⊗ point; spec `scratch/683/gate_spec.md`, battery `scratch/683/battery/`, census `scratch/683/census/`. On resume: SN census result; confirm landed.
- **Cone store list 3 marker repair (W3, rulings 51-53)** — CONFIRMED landed in the tree 2026-11-08 (uncommitted, `fix/rests-on-repair`), spec `scratch/verification_tree/cone_store/b5/census/repair/spec.md` + `table.json` (248 judged after ruling 53: 18 KEEP, 7 RE-POINT, 223 REMOVE). On resume: re-run `place` after the edits; expect 0 failing outside `traced_memo`.
- **Numerics reading pass design (W6, ruling 54)** — 2026-11-09: schema, `check.py`, controls, 28 packages, pilot P06 green; `scratch/verification_tree/cone_store/b5/reading/`. On resume: main's rulings on NEEDS (a)-(e); hand-read subject controls for the 24 packages without one. → `L132`
- Everything else is merged; the record is the SN theory page's development history and the archive.

## 3. Durable reference (reusable verification-design recipes)

Reusable RECIPEs / cited by `AGENT.md`. Core lessons in `lessons.md`; these keep the worked method.

- [SN sentinel harness](sn_sentinel_harness.md) — `@pytest.mark.sentinel` one-cheap-test-per-capability-node; cosmic-ray mutation-validation (copy the module aside first); per-NODE-sentinel-leaves-interior-uncovered gap.
- [SOTP separability verification](sotp_separability_verification.md) — separable ⟺ Cartesian-product per-axis; coupled physics → OperatorSum fallback; Route-A array_equal vs Route-B nulp; slab degenerate.
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
- [V&V tagging idioms](feedback_vv_tagging.md) — module `pytestmark` vs per-test `verifies()`; foundation carries NO `verifies()`; xfail `strict=True`+`reason=`.
- [Cross-method protocol design](feedback_cross_method_protocol.md) — reuse registry schema; `max(tol_a,tol_b)` agreement; L1-not-L4; verify truth values vs literature memos first.
