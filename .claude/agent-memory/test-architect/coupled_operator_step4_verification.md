---
name: coupled-operator-step4-verification
description: The RE-SCOPED Step-4 gate spec for the Coupled Block Operator campaign (2026-07-10) — N-general, semantics-agnostic block machinery with ψ½ as instance #1 (SUPERSEDES the bespoke-2×2 4d). 4d.0 = generalize FullField→System[Interior,Boundary] (L20 structure-only carrier, System A bit-identical for EVERY consumer, the multi-instantiation crux); 4d.1 = CoupledField/CoupledOperator machinery (the LOAD-BEARING assemble≡probe gate front-and-centre + per-block type-safety NET-NEW teeth + block-.H Mode-12 reciprocity); 4d.2 = build_coupled_system co-producing builder (presence=block-existence, subsumes step 6); 4d.3 = wire block apply + E4 two-anchor. Carries the durable invariants (assemble≡probe, A2 composite/seed reciprocity, two-anchor). Sibling of aba_schur_fold + a_ab_seed_injection + inverse_adjoint_coherence (L19).
metadata:
  type: reference
---

# Verifying Step 4 (RE-SCOPED 2026-07-10) — the N-general block machinery

**When this recipe applies.** A campaign has NAMED all block operators of a
coupled system (Steps 1-3: A_BB resolvent, A_BA fold, A_AB seeding — all landed
on `refactor/sn-walk-unification`; 4a/4b/4c DONE @ `19d9831`) and now builds
**N-general, semantics-agnostic block machinery** — a `CoupledField` (N systems)
+ `CoupledOperator` (N×N typed block grid) with three consumption modes
(`assemble`/`apply`/`solve`) — landing the ψ½ 2×2 as **instance #1**. The
operational criterion the user set: *"if we assembled the explicit matrix with
all blocks in it, it must work."* This SUPERSEDES the earlier bespoke-2×2 4d
(`OperatorSum` + present-zero padding was REJECTED — it keeps wrong
multiplications representable; `system_role` tags are runtime metadata, not type
constraints; the honest object is a typed block operator over a typed block
vector where the block matvec is the ONLY spelling and type-checks per block —
Pattern-1 ∘ Pattern-4). The prior memo's INVARIANTS carry forward (assemble≡probe,
Mode-12 `.H` reciprocity, two-anchor value); its STRUCTURE (bespoke 2×2 WRAP) is
replaced. Sibling of `aba_schur_fold_unweld_verification.md` (Step 2),
`a_ab_seed_injection_verification.md` (Step 3); A3 EXTENDS `test_inverse_adjoint_
coherence.py` (L19).

**The object (worked instance).** `System[Interior, Boundary]` = the carrier-generic
generalization of today's `FullField` (`transport/full_field.py`). System A =
`System[AngularFlux, AngularBoundaryFlux]` (≡ today's FullField, MUST stay bit-id);
System B = `System[RadialCharacteristicFlux, RadialCharacteristicBoundary]` (the ψ½
ray, its own interior⊕boundary composite). `CoupledField` = `[SystemA, SystemB]`;
`CoupledOperator` = `[[A_AA, A_AB],[A_BA, A_BB]]` with A_AA=L+C−S−B_a, A_BB=RayOp−B_b,
A_AB=`RadialCharacteristicSeeding` (already leaf-typed `LinearOperator[Radial…Field,
AngularField]`), A_BA=`RadialCharacteristicEmission` RE-TYPED FullField-embedded →
true `SystemA→SystemB` block. Carrying member = **sphere-GL S4 ONLY**; cyl/slab =
non-carrying CONTROL (degenerate 1×1); **≥2G every value row**.

## The FOUR landings — the equivalence bar SPLITS by sub-commit

| sub-commit | what lands | equivalence bar | oracle | RED-if |
|---|---|---|---|---|
| **4d.0 carrier** | `FullField`→`System[Interior,Boundary]` generic; System A ≡ today | **per-consumer bit-id (`array_equal`/nulp)** | old-vs-new VALUE for EACH of to_flat/from_flat/space.shape/metric/±/zeros/copy | a hardcoded-leaf-type in the generic body (L20 multi-instantiation) / a proxy passed for a broken consumer (ERR-063) |
| **4d.1 machinery** | `CoupledField`+`CoupledOperator` (numerics, semantics-agnostic) | **assemble≡probe principled-equiv (nulp/`rtol=1e-11`)** + type-safety RAISE | `_dense(CoupledOperator.apply, tpl)` block matvec | offset-swap / wrong-block placement / Euclidean block-`.H` (Mode-12) / wrong-block matvec accepted |
| **4d.2 builder** | `build_coupled_system` co-produces (op, space); presence=block-existence | **structural (constructability)** | the mesh R12a predicate | carrying-mesh missing System B constructs / non-carrying carrying it constructs |
| **4d.3 wire+anchor** | block apply WRAPs fused walk (bit-id); extract@4e | **bit-id vs `_dense(LC.apply,tpl)` (5.5e-16 floor) + two-anchor** | #284 row-6 oracle + φ=Q/Σ_t + k_inf | value off floor / anchor drift |

**The DURABLE load-bearing gate is 4d.1's `assemble≡probe`** (the RE-SCOPE's
centrepiece — "if we assembled the matrix it must work" IS this gate). It is
semantics-agnostic (synthetic block fixture), principled-equiv (NOT 0-ULP —
L16: sparse COO→CSR scatter order ≠ block-matvec apply order), and the ψ½
instance rides it via the existing `_dense`/`_loss` #284 machinery. The 4d.3
WRAP bit-id is a TRANSIENT characterization (retired at 4e when the un-weave
deliberately breaks it — L5 characterize-without-calcify + retire-as-you-go).

## 4d.0 — the `System[Interior,Boundary]` carrier generalization (L20 recipe)

**This is the L20 pattern EXACTLY** (structure-only ABC generalization + stringly→
typed generic, intended bit-identical). Apply the ERR-063/Mode-12 discipline: prove
neutrality PER CONSUMER with a direct old-vs-new VALUE comparison, NEVER a proxy
(snapshots-didn't-move / type-checks-pass / no-guard-raised are Mode-12-blind).

- **N1 — per-consumer neutrality, System A byte-identical (`TestSystemCarrierNeutrality`).**
  For a seed-carrying `_sphere()` AND a non-carrying `_cyl_product()`/`_slab()`,
  assert byte-identity for EVERY consumer of the generalized carrier: `to_flat()`
  (`array_equal`), `from_flat(flat, tpl)` round-trip (`array_equal` bulk⊕boundary⊕
  seed), `full_field_space.shape` (==), `full_field_space.inner_product(x,y)` +
  `apply_metric`/`apply_inverse_metric` (`assert_array_equal`), `__add__`/`__sub__`/
  `__neg__`/`__mul__`/`__truediv__` (`array_equal`), `zeros(...)`, `copy()`, the
  `mesh` property, `_recombine`. Capture the OLD values pre-carve via the root
  conftest `--capture-baseline` (or an in-repo golden built from `HEAD~1`); do NOT
  compare against a re-derivation on the new code (self-referential). 16 files
  import `FullField` — the neutrality scope is that consumer set (grep-audit).
- **N2 — THE CRUX: the multi-instantiation synthetic fixture (config-blindness).**
  Production exercises the generic at exactly ONE leaf-type instantiation —
  `System[AngularFlux, AngularBoundaryFlux]` (System A). **A `System[Interior,
  Boundary]` bug in the type-parametrized body — a hardcoded `AngularFlux`/
  `AngularBoundaryFlux` leaf in `zeros`/`from_flat`/`_recombine` — is INVISIBLE to
  a System-A-only suite** (the leaf-type "index" is DEAD at the single production
  instantiation; sharpens L1/L11). MANUFACTURE a SECOND instantiation:
  `System[ScalarFlux, ScalarBoundaryFlux]` (both confirmed importable
  `BulkField`/`BoundaryField` subclasses — the diffusion/CP leaves, NO ray third
  block). The LOAD-BEARING mutation is a hardcoded `AngularFlux` slipping into the
  generic body: it REDs the scalar instantiation (`type(sys.bulk)` wrong) while
  System A STAYS green — **the asymmetry IS the config-blindness evidence** (same
  shape as L20's `2*pos→pos` / L11's product-RED/LS-GREEN). Assert `type(sys_scalar
  .bulk) is ScalarFlux` AND `type(sys_a.bulk) is AngularFlux` from ONE generic path.
- **N3 — metric is a REGRESSION-GUARD + the ERR-067 tripwire (L20).** "Structure-only,
  metric descends per-leaf" ⟹ the carve MUST NOT touch the per-block metrics
  (`G_bulk=V_cell·w_n`, trace `|Ω·n|·w`, seed `V_cell`). Keep GREEN + UNMODIFIED:
  `test_starting_direction_metric.py::test_shipped_metric_block_values` (atol=0.0),
  the `FullFieldSpace` composite reciprocity gates (`test_full_field_space.py`,
  `test_g_adjoint_reciprocity.py`). ADD ONE tripwire: a mutation giving the generic
  carrier a UNIFORM/absent block metric (the ERR-067 ghost `G_sd≡0` at the seed) MUST
  RED the byte-exact `V_cell` gate. That guardrail is "NO metric lives on the
  `System` ABC — it descends per-leaf." Mode-12-SAFE (asserts the weight ARRAY,
  never a scalar functional; feeds a NONZERO seed + a control leg — L18).
- **N4 — MRO sibling-not-child (L20).** `TimedFullField(FullField)` (`timed_full_
  field.py:147`) MUST stay a correct subclass after `FullField` becomes/aliases
  `System[AngularFlux, AngularBoundaryFlux]`. Gate: `System in TimedFullField.__mro__`,
  `TimedFullField._recombine` returns a `TimedFullField` (empty history, preserved
  `history_depth`), `__post_init__` flat-buffer validation fires EXACTLY ONCE up the
  `super()` chain. **The break mode is collapsing an intermediate name** — every
  existing `isinstance(x, FullField)` / `isinstance(x, TimedFullField)` consumer MUST
  still resolve (blast radius: `streaming.py:211/1003`, `scattering.py:1796`,
  `isotropic_scattering.py:272/286/382/395`, `solver.py:2140/2218/2577/2772`,
  `windowing.py:133` — grep-audit + assert each holds). ⚠ **alias-vs-subclass is an
  IMPLEMENTER decision** (is `FullField = System[Angular…]` an alias, or `class
  FullField(System)`?); the gate is agnostic — "every `isinstance(_, FullField)`
  consumer resolves correctly" — but FLAG that System B is NOT a FullField
  (`System[Radial…]≠System[Angular…]`), so per-block leaf-type guards do the
  discrimination, not `isinstance(_, FullField)`.

## 4d.1 — the CoupledField/CoupledOperator machinery (semantics-agnostic, numerics)

Home: NEW `tests/numerics/test_coupled_operator.py` (semantics-agnostic, SYNTHETIC
toy block operators — the machinery is N-general, NOT SN). Build a 2×2 `[[D11,D12],
[D21,D22]]` of small asymmetric dense/diagonal `LinearOperator`s over two DISTINCT
toy System types (`System1`/`System2` with different leaf shapes) so offsets AND
type-safety are both activated and hand-verifiable.

- **M1 — block matvec value `y_i = Σ_j A_ij x_j`.** `CoupledOperator.apply(x)` vs a
  hand-summed reference on the synthetic 2×2 (`array_equal` on the toy dense blocks
  — the machinery adds no reduction). Non-carrying degenerates to 1×1 (identity to
  the single block). Config: **asymmetric blocks** `D12≠D21` (else an offset-swap is
  Mode-12-transpose-blind, L16); x = non-flat fixed-seed random (§0.6).
- **M2 — THE LOAD-BEARING assemble≡probe gate (FRONT-AND-CENTRE).**
  `CoupledOperator.assemble().as_matrix() ≡ _dense(CoupledOperator.apply, tpl)`.
  The `assemble` mode places each block at its `(row_i,col_j)` OFFSET into one flat
  Mat (the block structure PROVIDES the offsets — the scoped `LocalToGlobalMap`);
  `apply` is the block matvec — **two structurally-different computations of the same
  matrix**. **Principled-equiv `assert_allclose(rtol=1e-11)` / `nulp`, NEVER
  `array_equal`** (L16: `SparseAssembledOperator` COO→CSR scatter-sum order ≠ block-
  matvec apply order). ⚠ **TAUTOLOGY TRAP (L16):** the probe MUST be `_dense(.apply)`
  (the block matvec) — if `CoupledOperator.as_matrix()` DELEGATES to `assemble().
  to_dense()`, a gate `assemble().as_matrix() ≡ .as_matrix()` is vacuous; force the
  RETAINED apply-probed view. **Teeth:** an offset-swap (place A_AB at A_BA's
  `(row,col)`) / a wrong-block placement / a dropped block reds it O(1) — on the
  ASYMMETRIC synthetic (a symmetric block-grid is offset-transpose-blind, Mode-12).
- **M3 — per-block TYPE-SAFETY, NET-NEW teeth (L4/anti-#11).** `A_AB @ x_A` (wrong-
  component block applied to the wrong System) is (i) a RUNTIME raise (`match=` the
  specific "expected System_j got System_i" message) AND (ii) a pyright error
  (`assert_type` on the block signature `A_ij : System_j → System_i` + a
  `# type: ignore`-free negative that reveals `reportArgumentType`). **The block-
  type + presence guards have ZERO negative tests today → teeth are NET-NEW, not
  migrated** — write the POSITIVE control (`A_AB @ x_B` succeeds, correct value) AND
  the negative (`A_AB @ x_A` raises the specific message). A bare `pytest.raises(...)`
  is FALSE-green if a DOWNSTREAM crash (shape mismatch) satisfies it — `match=` the
  block-guard's own message.
- **M4 — block `.H` = block transpose, Mode-12-closed (THE HIGHEST-VALUE ROW).**
  `CoupledOperator.H` = the block-transpose `(A.H)_ij = (A_ji).H` with the COMPOSITE
  metric conjugation preserved. See the Mode-12 deep-dive below: A2a (composite-G
  forward reciprocity, extends `test_inverse_adjoint_coherence.py` with a `coupled`
  fixture) + A2b (seed-isolated V_cell, extends `test_g_adjoint_reciprocity.py`). A
  EUCLIDEAN block adjoint (skip the metric conjugation) reopens ERR-067 → reds ALL
  geoms (refutation #3).
- **M5 — SystemRole {A,B,COUPLED} join carries to CoupledOperator (foundation, NO
  `verifies`).** `CoupledOperator.system_role = COUPLED`; the join `_join_system_roles`
  (`operator.py:318`) already landed (4a). Explicit-STAMP the A_AA block = `SystemRole.A`
  (C-fwd elegance ruling — do NOT join-derive; transport-None poisons the join). Tooth:
  a join returning `a` unconditionally → A⊔B=A RED. ⚠ **MUST NOT lock in Optional[ray]**
  — no "non-carrying → ray role None" pin (step-6 subsumed by 4d.2 presence).

## 4d.2 — the co-producing builder (presence=block-existence, SUBSUMES step 6)

`build_coupled_system(sn_mesh, mat_xs) → (CoupledOperator, CoupledSpace)` emits the
operator AND the matching block-field space, ALIGNED by construction. Here A_BA
RE-TYPES FullField-embedded → true `SystemA→SystemB`, B_b→System B's boundary block,
and the ray LEAVES System A's FullField into System B's composite (folds step-6 ray-
extraction). Home: `test_psi_half_coupling.py::TestCoupledBuilder`.

- **P1 — alignment by construction.** `op, space = build_coupled_system(sn, mat_xs)`:
  `space` = `CoupledSpace[SystemASpace, SystemBSpace]`, and `op.apply(space.zeros())`
  type-checks + shapes align (`op.codomain == space` per block). A mismatched
  op/space pairing is unconstructable (the builder is the ONLY constructor).
- **P2 — presence-STRUCTURAL, NET-NEW teeth (subsumes step 6).** Keyed on the mesh
  R12a predicate (`sn.radial_characteristic_space is not None`): a carrying `_sphere`
  → 2×2 (System B EXISTS); a non-carrying `_cyl_product`/`_slab` → 1×1 (System B
  UNCONSTRUCTABLE). **The 7 presence guards have ZERO negative tests today** — write
  positive (carrying builds 2×2, non-carrying builds 1×1) + negative (a carrying
  `CoupledField` MISSING System B raises `match=`-specific; a non-carrying one
  CARRYING System B raises `match=`-specific), each with its positive control. This
  IS "applying a System-B block on a non-carrying mesh is unconstructable" —
  Pattern-4 illegal-states, mesh-keyed.
- **P3 — A_BA / B_b re-type does not move the VALUE.** After A_BA re-types to
  `SystemA→SystemB`, its block matvec on System A's bulk ≡ the pre-re-type
  `RadialCharacteristicEmission.apply` restricted to the ray codomain (`array_equal`
  — bit-id inheritance, the value did not move, only the typed carrier). Same for B_b
  (System-B boundary block ≡ the pre-re-type present-zero-padded corner). ⚠ this is
  a RE-POINT of the existing `TestA_BA_SchurFold` / `TestB_b_RayBoundary` probes from
  the FullField-embedded layout to the block domain/codomain.

## 4d.3 — wire block apply + the two-anchor keystone

- **W1 — block apply WRAPs the fused walk (bit-id, TRANSIENT).** `CoupledOperator.apply`
  delegates to the current fused `(L+C)` walk → `array_equal` vs `_dense(LC.apply,tpl)`
  (the #284 row-6 machinery). **Mode-11 sentinel (deliverable):** wrap the driver's
  route into `CoupledOperator.apply` (count>0) — the ONLY proof the iterate rewired to
  `CoupledField`; a green bulk value gate measures the UNCHANGED sibling (Mode 11).
  RETIRED at 4e when the un-weave breaks bit-id → migrate to the principled-equiv
  `assemble≡probe` (M2) which stays green at 5.5e-16.
- **E4 — the two-anchor keystone (L2, closed-form structural).** All equivalence gates
  (N1, M1-M2, W1) are bit-id INHERITANCE = necessary-NOT-sufficient. Anchor the FORWARD
  value structurally-independently: non-fissile reflective `_sphere(bc="reflective")` →
  `φ=Q/Σ_t` (fixed-source flat infinite-medium, the single most powerful curvilinear
  diagnostic); fissile A/2g+chi → k_inf (≥2G, Cardinal Rule). ⚠ **non-fissile → k_inf is
  a nan dead gate** (refutation #6) — φ=Q/Σ_t for the non-fissile, SEPARATE fissile for
  k_inf. Mode-12 lens: k_inf anchors the FORWARD only (`eig(K)`); A2 is its adjoint
  complement (zero adjoint coverage in k_inf, L19).

## The Mode-12 reciprocity deep-dive (M4 / A2 — the highest-value row)

Assembling `A.H` for the block grid and doing a EUCLIDEAN block adjoint (skipping the
composite-metric conjugation) silently REOPENS Mode-12/ERR-067. TWO complementary
gates, TWO regimes (do NOT collapse them):

- **A2a — composite-G forward reciprocity (L19 Gate-1, extends `test_inverse_adjoint_
  coherence.py`).** Add a `coupled` fixture (`build_coupled_system(sn,mat_xs)[0]`)
  alongside the existing `_loss(sn)` parametrization. `⟨A.apply(ψ),x⟩_G = ⟨ψ,b⟩_G` for
  `x=A.H.inverse().apply(b)`; `b` bulk-only source-carried (a random FULL b falsely
  reds ~1e-1 — off the free-DOF/null-space slots per the file's own `_rand_bulk`
  doctrine); ψ FULL (bulk⊕trace random, seed present-ZERO — the well-posed regime);
  G = composite `V·w`(bulk) ⊕ `|Ω·n|·w`(trace) ⊕ `V_cell`(seed). STRUCTURALLY
  INDEPENDENT of the reverse-scan (arithmetic = forward apply + metric only). Tooth =
  **M-ADJ-metric** (skip the G⁺/G wrap in `_AdjointOperator.apply`, `operator.py:1127`)
  → RED **ALL 3 geoms** ~0.13-0.33. ⚠ **the "slab stays GREEN" prediction is FALSE**
  (L19, burned on 2.5c — the composite bulk⊕trace metric is non-trivial on EVERY
  geometry: non-uniform mesh + bulk-vs-trace weight mismatch): assert RED on ALL
  three, NEVER design it as a slab-green discriminator. This is refutation #3's
  "Euclidean block adjoint reopens ERR-067."
- **A2b — seed-isolated V_cell reciprocity (L18/B_b, the ERR-067-specific catcher).**
  `⟨A_ss x,y⟩_{G_sd}=⟨x,A_ssᵀ y⟩_{G_sd}` under `G_sd=V_cell` (production
  `inner_product_weights`, `radial_characteristic_space.py:346`) on System B's
  interior self-block, with a NONZERO seed (activate the previously-nulled input,
  L18) + CONTROL leg (unmutated holds <1e-12) + teeth: (i) ghost `G_sd→0` reds under
  a nonzero seed; (ii) seed-row forward flip (`_seed_rows_forward` but NOT
  `_seed_rows_transpose`) → SPD metric sees the apply-vs-transpose inconsistency, RED
  >1e-6 (0.107 measured); (iii) asymmetric corner gauge V[+1]→2V[+1] RED 0.33. **BOTH
  legs mandatory** — mutated-alone is fooled by a broken baseline (also off-tol on the
  new input, mimicking "caught" → false green). Extend `test_g_adjoint_reciprocity.py`.

## Carry-forward vs RE-POINT ledger (the ray leaves FullField at 4d.2)

**CARRY FORWARD unchanged (the ray is still System A's optional third block through
4d.0-4d.1; the fused walk still computes it through 4d.3 WRAP):**
- `TestRegressionFloor` #284 certificate + the `_dense`/`_blocks`/`_bn`/`_template`/
  `_loss`/`_sphere` helpers + the row-6 principled-equiv oracle (5.5e-16) — the
  oracle VALUE becomes the 4d.1 `assemble≡probe` floor.
- `TestCoupledLift` (4c LIFT — A_BA-as-own-gain, S/F pure-bulk) — orthogonal to the
  carrier generalization, landed, stays.
- `TestSystemRoleLattice` (4a) — the role enum + join; extends to `CoupledOperator`
  (M5), no re-point.
- `TestA_AB_SeedInjection` — A_AB is ALREADY leaf-typed (`LinearOperator[Radial…Field,
  AngularField]`); its direct-leaf probes are unaffected by the carrier generalization.

**RE-POINT when the ray leaves FullField (4d.2) — the FullField-embedded-seed layout
→ CoupledField `[SystemA ⊕ SystemB]`:**
- `_template`/`_blocks`/`_dense`/`_seed_composite`/`_random_composite`/`_dense_seed`/
  `_v_cell_seed` — all build/probe the FullField `radial_characteristic` third block;
  after 4d.2 the seed is System B, so the seed-slice `slice(nb+nt, nb+nt+ns)` becomes
  System B's sub-block of the CoupledField flat layout. RE-POINT to the
  `CoupledSpace` layout (SystemA = bulk⊕trace; SystemB = ray⊕ray-boundary).
- `TestRegressionFloor::test_loss_operator_is_block_triangular_in_the_ray` — its
  `A_sb`/`A_bs` measurement on the augmented `(L+C)` re-expresses as the
  `CoupledOperator` block structure (A_AA self-block triangular; A_AB the seed→bulk
  feed) **at 4e (the un-weave)**, NOT at 4d.0. ⚠ Do NOT re-point the floor
  prematurely (refutation #7): through 4d.3 WRAP the fused walk keeps it valid.
- `TestA_BA_SchurFold` + the `RadialCharacteristicEmission` (A_BA) `_dense(A_BA.apply,
  tpl)` probes on FullField — RE-POINT to the `SystemA→SystemB` block domain/codomain
  when A_BA re-types (4d.2, P3).
- `TestB_b_RayBoundary` — from "present-zero trace on augmented FullField" → "System
  B's boundary block" (4d.2, P3).

## Mutation ledger (monkeypatch-only; NEVER `git checkout` an uncommitted file — L4)

| gate | tooth | target `file:line` (⚠ best-estimate, may drift) | expected RED |
|---|---|---|---|
| N1 | proxy for a broken consumer | any `System` generic-body edit | per-consumer value diff |
| N2 (CRUX) | hardcoded `AngularFlux` leaf in generic body | NEW `System.zeros`/`_recombine`/`from_flat` | scalar-instantiation RED / System-A GREEN |
| N3 | ABC uniform/absent block metric (ERR-067 ghost) | NEW `System`/`CoupledSpace` metric | byte-exact V_cell RED + control holds |
| N4 | collapse the `FullField` intermediate name | NEW `System` re-parent | `isinstance(_,FullField)` consumer RED |
| M2 | offset-swap A_AB↔A_BA / wrong-block placement | NEW `CoupledOperator.assemble` | O(1) matrix diff (asymmetric blocks) |
| M2 tautology | probe via `.as_matrix()` not `.apply` | NEW `CoupledOperator.as_matrix` | vacuous-green (flag, don't ship) |
| M3 | `A_AB @ x_A` accepted (drop the block-type guard) | NEW `CoupledOperator.apply` block dispatch | no raise → RED (+ pyright reportArgumentType) |
| M4/A2a | skip G⁺/G wrap (Euclidean block-`.H`) | `numerics/operator.py:1127` (`_AdjointOperator.apply`) | 0.13-0.33 ALL geoms |
| A2b-i | ghost `G_sd→0` (metric→zeros) | `numerics/spaces/radial_characteristic_space.py:346` | seed-recip RED + control holds |
| A2b-ii | `_seed_rows_forward` flip (NOT transpose) | `sn/loss_representation/__init__.py:2738` | >1e-6 (~0.107) |
| M5 | join returns `a` (drops differ→COUPLED) | `numerics/operator.py:318` (`_join_system_roles`) | A⊔B=A RED |
| P2 | carrying builds 1×1 / non-carrying builds 2×2 | NEW `build_coupled_system` | constructability RED `match=` msg |
| W1 | driver forgot to route iterate → CoupledField | NEW driver wire | Mode-11 count 0 |
| E4 | sign-flip C+K_iso / drop a block | `homogeneous`/`solve_sn` | k_inf O(1) / φ≠Q/Σ_t |

## Refutations (fire before the ink dries)

1. **The 4d.0 CRUX — System A is exercised at ONE leaf-type instantiation.** A
   `System[Interior,Boundary]` generic bug is invisible to a System-A-only suite;
   MANUFACTURE `System[ScalarFlux, ScalarBoundaryFlux]` (N2) — the hardcoded-leaf
   mutation reds scalar / greens angular. Without it, "bit-identical for every
   consumer" is a FullField-only proxy (Mode-12/L11). THE highest-value 4d.0 row.
2. **Per-consumer neutrality is a VALUE claim, not a proxy (ERR-063/Mode-12).**
   snapshots-didn't-move / type-checks-pass / no-guard-raised are all Mode-12-blind
   to a per-consumer divergence — assert old-vs-new VALUE (`array_equal`/nulp) on
   EACH consumer's actual output.
3. **The block `.H` does NOT inherit the metric** — a Euclidean block adjoint reopens
   Mode-12/ERR-067 (M-ADJ-metric reds ALL 3 geoms, NOT slab-green — the "slab stays
   green" prediction is FALSE, L19/burned-on-2.5c). A2a+A2b mandatory.
4. **assemble≡probe is principled-equiv, NOT 0-ULP** (L16: COO→CSR scatter order ≠
   apply order) — `rtol=1e-11`/nulp, never `array_equal`. And the probe MUST be
   `_dense(.apply)` (block matvec), NOT `.as_matrix()` if that delegates to assemble
   (tautology trap).
5. **type-safety + presence guards are NET-NEW teeth** (L4) — the block-type + 7
   presence guards have ZERO negative tests; write positive control + negative
   `match=`-specific-message (a downstream shape-crash false-greens a bare
   `pytest.raises`).
6. **non-fissile → k_inf is a nan dead gate** — E4 fixed-source `φ=Q/Σ_t` + SEPARATE
   fissile A/2g for k_inf, ≥2G.
7. **the ray-leaves-FullField RE-POINT is 4d.2, NOT 4d.0** — System A stays seed-
   carrying (bit-id) through 4d.0; the #284 floor certificate stays "within (L+C)"
   through 4d.3 WRAP; only 4e un-weaves. Don't re-point the floor / `_template`
   prematurely.
8. **do NOT lock in Optional[ray]** — presence = block EXISTENCE (carrying→2×2,
   non-carrying→1×1); NO "non-carrying System A → ray=None" pin (P2 subsumes step 6).
9. **cyl/slab are the non-carrying CONTROL** — CoupledOperator degenerates to 1×1
   (System A only); the offset-swap/type-safety mutations run on the CARRYING sphere.
10. **the MRO is sibling-not-child** (L20) — TimedFullField stays a correct
    FullField/System subclass; the break mode is collapsing the `FullField`
    intermediate name; System B is NOT a FullField (leaf-type guards discriminate).

## Result contract

Land in four sub-commits (4d.0 carrier / 4d.1 machinery / 4d.2 builder / 4d.3 wire),
4e = the walk un-weave (EXTRACT, unchanged in spirit — retire the 4d.3 WRAP transient,
flip `assemble≡probe` to the durable floor). NEW files: `tests/numerics/test_coupled_
operator.py` (M1-M5, semantics-agnostic synthetic blocks) + `tests/sn/operators/
test_system_carrier_neutrality.py` (N1-N4, the L20 carrier). Extensions:
`test_inverse_adjoint_coherence.py` (A2a `coupled` fixture) + `test_g_adjoint_
reciprocity.py` (A2b) + `test_psi_half_coupling.py` (P1-P3, W1, E4; re-point
`_template`/`_blocks`/A_BA/B_b probes at 4d.2). Every tooth mutation-verified
in-process under `-O` (`np.testing`/`pytest.fail`, never bare assert; monkeypatch to
revert). Full `tests/sn -m "not slow"` + `tests/numerics` + ratchet transport:1 +
sphinx -W are the end-to-end acceptance (the eigenvalue/fixed-source wall the
operator-level gates don't own). The load-bearing deliverable is **M2 assemble≡probe**
("if we assembled the matrix it must work") + **M4/A2 Mode-12 reciprocity** (the
highest-value row) + **N2 the multi-instantiation crux** (the 4d.0 config-blindness).
