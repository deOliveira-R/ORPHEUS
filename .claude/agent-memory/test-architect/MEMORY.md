# Test-Architect Memory Index

One line per entry — detail lives in the linked file, NEVER inlined here (the
index is loaded whole every dispatch; keep it small). Four sections: (1) lessons
— READ `lessons.md` FIRST every dispatch; (2) active/in-flight state — git-true
(reconcile "unmerged" claims against git before acting); (3) durable reference
recipes; (4) design idioms. The failure-mode taxonomy lives in `vv-principles`;
the reference inventory in `AGENT.md` §2. No campaign play-by-play
here — it is merged archaeology.

## 1. Lessons

- **[Lessons — hot digest](lessons.md)** — 527 lines. The entries no rule, skill or definition clause carries (the workflows rule, invariant 6), each with a `→ LNN` pointer into the archive. Read it whole, every dispatch. Pruned 2026-09-22 by the agent-definitions audit: 73 entries restated a clause or are now carried by the definition.
- **[Lessons — cold archive](lessons_archive.md)** — 11573 lines, sections L1–L97, append-ordered: the war stories, the carve-archetype lookup (§L87) and the pure-math primitive long form (§L88). Open one section at a time, when the digest's pointer says the detail matters.

## 2. Active verification work

Merge status comes from git and GitHub, never from this list (`process-discipline`).

- **#432 — the axis-parameterised O(2) member** — issue OPEN; the stabiliser gates shipped. → **`L70`**
- **#235 — the 2-D angular closure's ranking instrument** — design delivered, issue OPEN. → **`L48`**
- **#358 — the test-dependence DAG** — memo delivered, issue OPEN. → **`L55`**
- **#405 P1 — the specification/question carve** — spec delivered 2026-09-25 (`.claude/plans/reference_p1_spec.md`), gates S1.1–S8.5; step 3 split 3a/3b/3c in its §1.3a (2026-09-29, 3a gates + battery written; battery at `scratch/reference_architecture/p1step3/battery3a/`). → **`L95`** Step 3c gates written 2026-09-29 (`test_mesh2d_face_laws.py`, `TestTheDefaultsAreGone`); battery owed on resume.
- **Reflective cleanup (W3, `refactor/reflective-is-a-mirror`)** — gates delivered 2026-10-01: spec `scratch/boundary_ontology/reflective_cleanup_gates.md`; carry fixture + batteries `scratch/boundary_ontology/reflective_gates/`. On resume: re-run the carry gate post-carve.
- **Platform drift W2 (`fix/platform-independent-quadrature`)** — gates (a) CR Gauss rules, (b) GL fingerprint, (c) exact k∞ + bound, spec `scratch/platform_drift/gates_spec.md` (2026-10-01). On resume: confirm landed, re-run `cr_sim`/`r1_sim` arms. → **`L98`**
- **#405 P1 step 5 (content identity, W3, `refactor/content-identity`)** — gates re-specified 2026-10-02 (spec §1.5; 4 files `tests/gates/*/test_content_identity*.py` + `_content_identity_helpers.py`; battery `scratch/reference_architecture/p1step5/gates_ta/battery/`). On resume: the owed post-carve arms (per-type, S5.6 per-space, S5.9, S5.10) and the green run. → **`L99`**
- **#405 P1 step 6 (phase-space functions, W3, `feature/phase-space-functions`)** — §1.6 re-specified 2026-10-02 (S6.1–S6.21, design (B): retraction `.H` = closed-form pullback π*); probes + plugins `scratch/reference_architecture/p1step6/ta/`. On resume: confirm landed; run S6.6 route arms post-carve.
- **#405 P1 step 7 (question values, W1)** — §1.7 re-specified 2026-10-02 (S7.1–S7.13); drafts + prototype + battery `scratch/reference_architecture/p1step7/ta/`. On resume: confirm landed, re-run the 14 arms on the real module, pin S7.12. → **`L100`**
- **#405 P1 step 8 (the specification, W1)** — gates landed `eaa74163` (two types: `InfiniteMediumSpecification | GeometrySpecification`); final battery `scratch/reference_architecture/p1step8/ta/battery/final/` (40 arms). On resume: only if the re-review reopens it. → **`L101`**
- **#405 P2 (solutions and certificates, W3)** — spec `.claude/plans/reference_p2_spec.md` (2026-10-03, gates R1.x-R8.x); step-1 drafts + prototype + 21-arm battery `scratch/reference_architecture/p2/ta/`. Steps 1–2 merged `55870ddc`; step-3 battery on `9f8d1f7f` (`ta/step3/`); §1.4 re-specified, gates + 14-arm battery on a prototype (`ta/step4/`). Step 4 WITHDRAWN to P4 (ruling 2026-10-03); §1.5 re-specified, step-5 gates in `ta/step5/` (worktree of `99862745`). Step 5 merged `a3ff64d0`; step 6 gates landed `4b724f04`, battery in `ta/step6/` (E3 repaired by a new row). Next: step 7a. → **`L102`**, `L103`
- **#405 P2 step 7b.1 (the uncertified reading, W3)** — gates + 32-arm battery delivered 2026-10-03 (spec §1.7b.1; battery `scratch/reference_architecture/p2/ta/step7b1/battery/`). Review round re-run 2026-10-03: 36 arms, 0 blind; #568 xfail row. → **`L104`**
- **#405 P2 step 7b.2 (certify_agreement migration, W3)** — §1.7b.2 re-specified 2026-10-03 (R7b2.1–13; 7b.2.0 region labels gates `tests/gates/mesh/test_mesh1d_regions.py` + 2 re-posed pins; probes `scratch/reference_architecture/p2/ta/step7b2/`). On resume: confirm landed; build the battery; measure the cylinder E2 reading's cost.
- **#405 P3 (the traced memo, W3)** — spec `.claude/plans/reference_p3_spec.md` (2026-10-04, gates M1.x–M5.x); prototype + gates + battery in `scratch/reference_architecture/p3/ta/` (`wt` prototype, `wt0` pristine first reds, `wt2` acceptance run). On resume: confirm landed; re-run the battery arms on the real module. → **`L105`**
- **#405 follow-up, the two slow files' gates (W3, `test/slow-reference-gates`)** — repaired 2026-10-04, evidence `scratch/reference_architecture/p3/gates_repair/`. On resume: confirm landed; after #200 re-time the slab order row. → **`L106`**
- **#200 sweep-preconditioned Krylov gates (W3, `fix/krylov-sweep-preconditioner`)** — delivered 2026-10-04, evidence `scratch/reference_architecture/p3/krylov200/`. On resume: confirm landed.
- **#200 follow-up gates (W3, same branch)** — ERR-053 re-homed + value row, inner-record rows (L0 + SN), 2-D fixture row; delivered 2026-10-05, evidence `scratch/reference_architecture/p3/gates200b/`. On resume: confirm landed; tag the record rows with the new ERR id. → **`L107`**
- **Chord-oracle hoist gates (W3, `refactor/chord-oracle-axial-lift`)** — delivered 2026-10-05: cylinder row retired, axial-cosine column law added; evidence `scratch/reference_architecture/p3/oracle_hoist/gates/`. On resume: confirm landed.
- **Geometric kernel seed (W5→W1, `characteristic_reference_architecture.md`)** — verification spec delivered 2026-10-05: `scratch/characteristic_architecture/seed_verification_spec.md`, probes `seed_spec_probes/`. Gates LANDED (uncommitted) 2026-10-05 on `feature/geometric-kernel-seed`: 5 files `tests/gates/geometry/test_{line,chart,chord,line_measure,kernel_corroboration}.py`, 34-arm battery + README `scratch/characteristic_architecture/seed_gates/`. On resume: confirm landed.
- **Characteristic reference P1 + P0.5 (W3)** — verification spec delivered 2026-10-06, revised the same day to the Krein/Galerkin/panels rulings: `scratch/characteristic_architecture/p1_verification_spec.md` (49 new gates, roster 122/99/77 of 298), probes `p1_gates/`. On resume: confirm landed; run the §12 battery.
- **Characteristic P1 step (a), the kernel verbs (W3, `feature/characteristic-kernel-verbs`)** — gates written 2026-10-06: `tests/gates/geometry/test_chord_transits.py`, `test_chart_directions.py`; 18-arm battery `scratch/characteristic_architecture/p1_step_a/battery/` (`run.sh <arms>`). On resume: the catalogue row (S^2/D_1h) reds until the entry lands; then run its owed quarter arm.
- **Characteristic P1 step (b) first rung, walls + line closure (W3, `feature/characteristic-walls-closure`)** — gates written 2026-10-06: `tests/gates/derivations/test_characteristic_{walls,closure}.py` (164 rows); 30-arm battery `scratch/characteristic_architecture/p1_step_b1/battery/` (`run.sh <arms>`). On resume: confirm landed; move the listed rows to `l0` once the archivist mints the labels.
- **Characteristic P1 step (b) rung 2, basis + line transport (W3, `feature/characteristic-basis-transport`)** — gates written 2026-10-06: `tests/gates/derivations/test_characteristic_{basis,transport}.py` + `_characteristic_mp.py` (mpmath refs); spec + battery `scratch/characteristic_architecture/p1_step_b2/` (`battery/run.sh <arms>`). On resume: confirm landed; move T1-T3/T7-T9 to `l0` once `characteristic-traversal-integrals` exists.
- **Characteristic P1 step (b) rung 3 gates (W3, `feature/characteristic-rung3`)** — written 2026-10-07: `tests/gates/derivations/test_characteristic_assembly.py` (33 fns), `tests/gates/geometry/test_{measure_density,line_domain}.py`; spec `scratch/characteristic_architecture/p1_step_b3/spec.md`; battery `.../gates/battery/` (`run.sh <arms>`, `battery_table.md`, 54 arms). On resume: confirm landed; mint labels then move planned levels. → **`L108`**
- **#586 cylinder rows restored (W3, `refactor/characteristic-line-cost`)** — 2026-10-07: thick cylinder escape legs tau 30/100 at 16 (+ERR-101 via A16), `_CLOSED_SLOW` both groups; probes `scratch/characteristic_architecture/p586/`. On resume: confirm landed.
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
