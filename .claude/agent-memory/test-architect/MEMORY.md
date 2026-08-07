# Test-Architect Memory Index

One line per entry — detail lives in the linked file, NEVER inlined here (the
index is loaded whole every dispatch; keep it small). Four sections: (1) lessons
— READ `lessons.md` FIRST every dispatch; (2) active/in-flight state — git-true
(reconcile "unmerged" claims against git before acting); (3) durable reference
recipes; (4) design idioms. The failure-mode taxonomy lives in `vv-principles`;
the reference inventory + XS mixtures in `AGENT.md`. No campaign play-by-play
here — it is merged archaeology.

## 1. Lessons — a HOT digest over a COLD archive (READ the digest at START)

- **[Lessons — hot digest](lessons.md)** — ~650 lines. One imperative rule per
  entry, grouped by behavioral family (gates that cannot red · harness discipline ·
  config blindness · reference & claim layer · tolerance · carve archetypes ·
  snapshots & exactness · pure-math primitives). **Read this file whole, every
  dispatch.** Every entry ends with a `→ LNN` pointer into the archive.
  ⚠ **A digest entry recording a GAP needs a LANDING note the moment the gap
  closes** — `L3`'s "no SN MMS exercises `q.boundary ≠ 0`" was true when written,
  the §4.6 fix landed, and the stale entry then generated a whole phase brief for
  work already done (`L40a`).
- **[Lessons — cold archive](lessons_archive.md)** — 3719 lines, sections L1–L40,
  append-ordered. The war stories, measured numbers, `file:line` detail and
  per-fixture tables. **Open ONE section at a time, only when the digest's pointer
  says the detail matters.** Never read it whole — that is ~48K tokens.
- NO lesson content is inlined here. The digest is the index over the archive;
  this file is the index over everything else. New lessons: add the RULE to the
  digest (with its `→ LNN`) and the war story as a new archive section.

## 2. Active / in-flight verification work

**Detail → [active campaigns](active_campaigns.md).** One line each here; open that
file for headlines, open findings and measured numbers. Reconcile against git first.

- **G2 geometric-transformation machinery** (`orpheus/geometry/transformation.py`) — **SHIPPED 2026-08-03**: `tests/geometry/test_transformation.py`, 42 gates / 96 cases, 8.9 s, pyright 0, **32/32 mutations caught, 0 blind**. Plan + closeout `scratch/g2_verification_plan.md` §12; lessons L35. G1 landed mid-plan and the coordinator then discharged A14/A15/A4/A1 (`cc5703ed`: `embedded_in`, the `Permutation` type, `NotAFinitePointGroupError`, the `SO(1)` message, `_ORTHOGONALITY_ATOL` 1e-10→1e-12). **A6 (typed failure reason) deliberately DEFERRED** ⟹ O4's input-side guard isolation is permanent. NEXT = G3.
- **G5 self-paired deck collapse** (`IdentityMap`+`SpecularMirror` → one `RigidMotion`-parameterised type) — PRE-carve plan `scratch/g5_verification_plan.md` (10 gates G1–G10, 16-mutation battery + 2 positive controls), lessons L36. Baseline `1 failed, 743 passed in 9.40s` (the 1 = the known task-#33 red). **Three refutations of the brief**: involution does NOT imply self-paired (`inversion`/half-turns admitted); the offset is bit-identically invisible; the cost estimate was 35× high. **Two existing gates DECAY to tautologies** under the collapse. `_orbit_closure` already takes `RigidMotion`, so the feared lossy axis-inverse is avoidable. OPEN user question: guard A vs guard B (§8.3).
- **G6.3 step 5 — bind the specular deck permutation `Γ₊→Γ₋`** — plan `scratch/g6_step5_verification_plan.md` (10 gates G1–G10, measured battery, §7 has 4 user rulings), lessons **L38**. **All four brief items landed MID-DISPATCH (uncommitted)** ⟹ POST-carve plan. ⛔ **Measured: the binding is gated by NOTHING** — drop / SWAP / collapse-to-one-space each = **0 new reds** over 1668 tests, while the identity-perm positive control = **+23**. ⛔ **And it evaporates at `TensorProductOperator`** (realizer output `domain=None`), so **step 8 cannot fire as written** — step 3 has the same hole. `is_involution` is quadrature-dependent unbound (lebedev already False ⟹ a lebedev-only gate is green on a no-op); `.H` is metric-blind EXACTLY (`G₋ == G₊∘π` bit-exact) ⟹ 2-nulp + a ×3 negative leg (`2/3`, `2`). Fast slice **24 s**, not the 5m45s the ladder implied. OPEN Q1: retire `is_involution` (zero consumers + self-conceding docstring)?
- **G6.3 step 7 — the deck-transformation uplift** (`SpatialWrap(axis:str)` + `_specular_kernel` → `PairedDeck(motion)` + one `_deck_kernel`) — spec `scratch/step7_deck_uplift_verification_plan.md`, lessons **L41**. **SHIPPED: `tests/geometry/test_paired_deck.py` (63 rows) + `tests/sn/operators/test_deck_kernel.py` (78 rows), 1.3 s, pyright 0, both batteries measured.** ⛔⛔ **The briefed contract was UNRUNNABLE**: `RigidMotion.permutes` is the AFFINE action, so it returns `None` for every deck element carrying a translation — the wrap included (main agent fixed it independently via `Quadrature.ordinate_permutation` → `linear_part.preserves`). ⛔ The brief's keystone fixture `product(4,4)` is DEGENERATE (local perm `== arange`); the keystone runs on `product(4,8)`. ⛔ π-vs-π⁻¹ reds only by RAISING ⟹ a second in-range mutation is mandatory. ⭐ Three EXCLUSIVE catchers: swap→binding rows (11), duplicate-body→the wrap counter (1), delete-`Permutation`'s-clause→the ERR-073 row (1). §8 migration debt DISCHARGED at `e13313a8`; residual = `reflection_index` (33 refs) which shelf-lifes one gate, and the periodic METRIC row (§6.1, out of scope).
- **G6 every operator knows its spaces** (`LinearOperator.domain/.codomain` `Optional=None` ⟹ `.H` silently degrades to the bare Euclidean transpose) — PRE-carve plan `scratch/g6_verification_plan.md` (999 lines; 7 gates, 8-mutation battery + 2 controls, MEASURED positive control), lessons L37. **The ⭐ metric gate IS buildable** — the Lambertian `AngularAverageOperator` is the ONLY non-weight-preserving boundary operator (4 of 5 shipped laws are Mode-12 blind); law holds 2.2e-16, degradation 5.5e-1, `domain:=None` reds it. **6 acceptance corrections C1–C6**: the shipped Lambertian is `is_adjointable=False` (needs B5's `apply_transpose`); "bit-identical under a layout permutation" is arithmetically IMPOSSIBLE for a reduction (25 % only) → use the FACE-PACKING permutation instead; **G6.4 cannot be boundary-scoped** — the `None` default is on the shared base and reaches `homogeneous/solver.py`. **DECAY LIST = 12 strict xfails + ~8 gates in 3 flavours (DIE/DECAY/INVERT), BLOCKING for the G6.4 commit.** OPEN user rulings: C1 (absorb B5?), C3 (boundary-only vs tree-wide: ~20 vs ~150 tests), C4 (space-generic escape type).
- **#325 symmetry-exact circle nodes** — plan + SHIPPED `tests/numerics/test_roots_of_unity.py` (17 gates, all mutation-proven), lessons L34. ⛔ **BLOCKED on a USER physics ruling**: exactness manufactures `argsort` ties worth **1.008 %** flux.
- **B3.4c periodic → partner face** — plan + dry-run gate module (95/95), lessons L33. All seven production steps landed mid-plan; **the TEST surface does not exist** and the two live strict xfails cannot detect the change.
- **P3 affine boundary source channel** (retire `IncomingSourceOperator`; `realize(PrescribedInflow)` → the zero morphism) — plan `scratch/p3_verification_plan.md` (8 sections, 4-candidate evaluation, 6-mutation battery + 2 controls), lessons **L39**. ⛔ **The campaign's "prophylactic, fenced twice" premise is FALSE**: the `block_role=None` stamp fences nothing (`|B(0)|=2.5`) and Krylov HARD-RAISES (`‖Aψ−q‖/‖q‖=1.718`) on BOTH sides of P2′ ⟹ P3 is a **fix**. ⭐ Keystone = the boundary law evaluated on the CONVERGED answer (`γ₋ψ|_f` reads bit-exact `2q`/`q`/`0`); superposition-in-`q` is a Mode-12 non-catcher. **The production carve + 4 test migrations landed DURING the dispatch** (8 reds → `555 passed`, 24.5 s) ⟹ §8 is the live reconciliation and carries the **11-item residual gap list** — Gates A/B/D all still ABSENT, `|B(0)|=2.5` lives only in a COMMENT, `DEAD TARGETS 2`, a label-less `:ref:`, battery un-run.
- **P4 non-trivial MMS through the DECLARED inflow channel** — design `scratch/p4_mms_design.md` (13 sections, 9-gate matrix, 8-mutation activation battery, 15 refuted candidates), lessons **L40**. ⛔⛔ **The brief's central risk was a LANDED fix**: 4 of 13 MMS builders are already non-vacuum + anisotropic (`build_slab_2g_nonvacuum_mms_case`), so P4 collapses from "build a reference" to **RE-ROUTE one** — declared ≡ hand-supplied `q_∂` measured BIT-IDENTICAL. ⭐ Orders `[2.041, 2.010, 2.003]`; the 8-mutation battery reds ×52–×4376 (`qiso` ×725 ⟹ the anisotropy is CONSTRAINED, not merely activated). ⛔ Three further refutations: the spec is evaluated at **TWO** shapes (`(N,)` guard probe + `(|Γ₋|,ng)`); the sphere is **NOT** refused (only its transpose is); **no public entry point can be handed a declared law** (`solve_sn_fixed_source` rebuilds the `SNMesh`) ⟹ P6's first item. ⛔ The brief's `B(0)=0` matvec row is **tautological** on the MMS fixture (`|B(x)|=0.0`) ⟹ moved to P5 on prescribed+reflective (`|B(x)|=1.320`).
- **Boundary machinery B3.2 → B3.4b + #21** — ALL DELIVERED on `refactor/operator-strategy-layers`; lessons L29–L32. Six open follow-ups incl. a LATENT `_law_from_tag` bug and a real orientation mismatch.
- **Three-DOF separation (operator ∥ splitting ∥ realization)** — full spec + **P0 SHIPPED** (`tests/sn/architecture/`), lessons L24–L26. The proposed AC was signature-tautological; the real red is the WELD.
- **#280 RESIDUE** — binding spec `.claude/plans/residue_verification_spec.md`, PRE-carve. Three contract findings SHAPE the carve; the LD `.H` moment-mass metric is load-bearing.
- **DSA for SN (#2)** — 13-gate battery `.claude/plans/dsa_verification_spec.md`, PRE-carve, runs AFTER #280. Lessons L23: the correction→0 PARTITION (7 of 8 canonical errors are FP-invisible).
- **Coupled Block Operator step-5 (#41)** — [spec](coupled_operator_step5_solve_verification.md); awaiting user rulings D2a/D3/D5.
- **Diffusion #290** — [plan](diffusion_integration_290_verification_plan.md); **MERGED @ `3a19133`**, so the DSA consumable is built.
- **A3/#280 Phase-2.5 walk-unification** — gate files authored, surface on `main`; recipe [a3-reverse-scan](a3_reverse_scan_transpose_verification.md), lessons L17/L18/L19.
- **Prior SN campaigns** (#206/#208/#236/#240/#247/#251/#257/#18/#19/#20) MERGED to `main` — NOT open work.

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
