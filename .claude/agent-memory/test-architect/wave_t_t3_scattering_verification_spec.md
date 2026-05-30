---
name: wave-t-t3-scattering-verification-spec
description: Verification SPEC for Wave T step T.3 — lifting ScatteringOperator.apply into a SumOfTensorProductsOperator. Pre-implementation test plan to drive turn-by-turn T.3 execution.
metadata:
  type: project
---

# Wave T T.3 — Scattering as `SumOfTensorProductsOperator` — verification spec

**Branch.** `refactor/sn-operator-algebra` (worktree
`moment-space-and-layering`). T.3 follows the
`feedback_no_method_implementer_for_surgical_carves` rule: main agent
direct authorship with turn-by-turn user steering. This spec is the
test driver.

**Plan source.** `.claude/plans/wave_t_tensor_network.md` §6 step
T.3 + §5 verification strategy + §7 risk row 4 (highest leverage —
scattering is in the Krylov hot path of every SN solve).

**Grand report.** §15.2 lines 2046-2086. The §15.2-named identity is

```
S = SumOfTensorProductsOperator([
    SigmaMoment(xs.scatter, ell) & AngularMomentOperator(ell) & GroupScatteringMatrix(xs.scatter, ell)
    for ell in range(xs.scatter.order + 1)
])
```

i.e. `Σ_ℓ` of three-factor `TensorProductOperator` summands, one per
Legendre order. **NOTE: This naming is grand-report aspirational.**
This spec resolves the concrete factor design in §6 below.

**Bifurcation point (V&V pillars).** T.3 is a code-internal refactor.
The verification chain runs against pre-T.3 numerics (regression /
snapshot — see §3 for the principled-equivalence three-criteria gate)
backed by structurally-independent ground truth at L1 — homogeneous
`k_∞ = νΣ_f / Σ_a` (closed-form pillar) and P1 MMS (MMS pillar). No
eigenvalue claim is asserted from MMS evidence.

---

## §1 Scope + boundaries

### In scope (the lift)

1. **`ScatteringOperator.apply`** at `orpheus/sn/scattering.py:805-988` —
   the four `singledispatchmethod` arms (`TimedFullField`,
   `AngularFlux`, `ScalarFlux`, fallback `raise TypeError`). The
   bare-`np.ndarray` arm RETIRED in D-I.2 (commit per scattering.py:632)
   — there are NOW only three dispatch arms, plus the unsupported-type
   `raise TypeError` default.
2. **`build_aniso_source`** at `scattering.py:527-651` — the inner
   `R · Λ · M` pipeline. T.3's lift target.
3. **`LegendreMomentScattering.apply`** at `scattering.py:227-270` —
   the Λ inner numerics. Already algebraically the §15.2 inner
   summand (per the class docstring); T.3 may keep this as the
   `GroupScatteringMatrix` factor body OR wrap it in a thinner typed
   factor. Decision deferred to §6 architectural questions.

### Out of scope (deferred)

1. **Per-material per-ℓ einsum** at
   `orpheus/sn/material_xs_field.py:515-572`
   (`MaterialXSField.apply_legendre_scattering_moments`) — the inner
   numerics. Bit-identity preserved per plan §6 T.3 paragraph 3 ("The
   per-material per-ℓ einsum at material_xs_field.py:515-572 is the
   inner numerics. Bit-identical preserved."). T.3 wraps it; does
   NOT change it.
2. **`MomentProjection` and `HarmonicMomentReconstruction`** —
   `orpheus/numerics/projection.py:240-599`. Stay as-is.
3. **`foldable_part` / `residual_part`** — those sibling operators
   construct new `ScatteringOperator` instances via different XS
   data; they will inherit T.3's lift automatically because they
   route through the same `.apply`.
4. **`SumOfTensorProductsOperator._build` fusion** — out per plan §4.
5. **`apply_transpose` advertisement** on the lifted kernel —
   scattering currently has `{CAP_APPLY}` only. T.3 keeps it that
   way; adjoint surfaces in a later wave.
6. **PrescribedInflow / affine source compositions** — not relevant
   to scattering.

### Bit-identity discipline (the contract gate)

T.3 introduces a `Σ_ℓ` outer sum across summand `apply` results. The
existing inline path's outer sum is the per-ℓ Python for-loop inside
`MaterialXSField.apply_legendre_scattering_moments` (lines 558-571)
which accumulates `out[l, :n_m]` per-ℓ block-by-block — NOT a single
flat reduction.

**The reduction order of the `R · Λ · M` pipeline DOES NOT change**
under T.3 if the lift is `R ∘ Λ ∘ M` (one composition) rather than
`Σ_ℓ (R_ℓ ∘ Λ_ℓ ∘ M_ℓ)` (per-ℓ summand sum). See §6 Q1 for the
factor-shape decision.

**Default expectation.** T.3 keeps the inner `R · Λ · M` reduction
verbatim and lifts only the EXTERNAL apply shape. Result: **strict
bit-identity at the matvec level** for all four dispatch arms. The
nulp-relaxation route (per `vv-principles` §"Bit-identity vs
principled-equivalence") is available IFF the §6 Q1 decision elects
the per-ℓ summand form. The spec covers both routes (see §3 "Pass
criterion" column).

---

## §2 Pre-T.3 existing-test inventory

Two scattering-specific test files exist; one fission test class
already pioneered the T.2 lift pattern. The L1 MMS file is the
load-bearing downstream gate. All counts below are at the spec time
(2026-05-30).

### `tests/sn/test_scattering_operator.py` — 13 classes, ~50 tests

| Class                                | Line  | What it pins                                                                                | T.3 disposition       | Reason                                                                                                                                                                                                                                                                                          |
| ------------------------------------ | ----- | ------------------------------------------------------------------------------------------- | --------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `TestProtocolCompliance`             | 125   | `isinstance(op, LinearOperator)`; `capabilities == {CAP_APPLY}`; `apply` accepts AngularFlux | **KEEP**              | Capability set unchanged.                                                                                                                                                                                                                                                                       |
| `TestBitIdenticalExtractionP0`       | 154   | `add_iso_source` / `add_n2n_source` / delegators vs per-cell reference                      | **KEEP**              | T.3 does NOT touch P0 in-scatter / (n,2n) — those stay in `MaterialXSField.apply_p0_in_scatter` / `apply_n2n`. Bit-identity preserved.                                                                                                                                                          |
| `TestAnisotropicScatteringExtraction`| 240   | `build_aniso_source`: P0 returns None; iso flux → 0; delegator matches operator             | **KEEP / RELAX-nulp** | If T.3 keeps the inner `R∘Λ∘M` reduction verbatim (§6 Q1 Route A), KEEP at `rtol=1e-13`. If T.3 elects per-ℓ summand form (Route B), relax `test_delegator_matches_operator` to `nulp=4*L` per principled-equivalence gate. Current `assert_array_equal` at line 301 needs nulp swap for Route B. |
| `TestApplySemantics`                 | 309   | `apply(iso ψ)` matches manual Q_iso/sum_w broadcast; linearity; zero-input                  | **KEEP**              | Linearity is structural (must hold for any linear lift). Iso ψ test is `rtol=1e-13`. Both routes pass.                                                                                                                                                                                          |
| `TestProducerSideNormalisation`      | 390   | Pattern 7 / R-1 Step 4 A1 `/sum_w` invariant at apply boundary                              | **KEEP**              | T.3 does NOT change the producer-side normalisation. The `/sum_w` factor stays at the apply boundary (`scattering.py:929-930, 986-987`).                                                                                                                                                        |
| `TestCompositeInvariants`            | 462   | TimedFullField composite → composite output; implicit-zero BC; history depth                | **KEEP**              | Composite contract unchanged. T.3's lift is internal to bulk apply.                                                                                                                                                                                                                             |
| `TestP0AlgebraicIdentities`          | 524   | P0 uniform-flux hand-calc; (n,2n) doubling factor                                           | **KEEP**              | P0 + n2n unchanged.                                                                                                                                                                                                                                                                             |
| `TestFoldablePart`                   | 633   | `foldable_part` structure (scattering_order=0; diagonal P0; zero n2n)                       | **KEEP**              | Independent of T.3 lift.                                                                                                                                                                                                                                                                        |
| `TestResidualPart`                   | 710   | `residual_part` structure (zero diagonal P0; carries Pℓ + n2n verbatim)                     | **KEEP**              | Independent of T.3 lift.                                                                                                                                                                                                                                                                        |
| `TestFoldableSigma`                  | 783   | `foldable_sigma()` shape + diagonal extraction                                              | **KEEP**              | Independent.                                                                                                                                                                                                                                                                                    |
| `TestAlgebraicIdentity`              | 825   | `S.apply ≈ S.foldable_part().apply + S.residual_part().apply` at `rtol=1e-14`               | **KEEP / RELAX-nulp** | Currently `rtol=1e-14, atol=1e-15`. If T.3 elects per-ℓ summand reduction (Route B), the algebraic identity still holds at FP-non-associativity (rtol=1e-13 should pass for L=1, P1 typical). KEEP rtol=1e-14 first; relax only if it actually fails on the route chosen.                       |
| `TestPurity`                         | 935   | `foldable_part() / residual_part()` are pure functions                                      | **KEEP**              | Independent.                                                                                                                                                                                                                                                                                    |
| `TestIsFoldableIntoSigmaR`           | 974   | `is_foldable_into_sigma_r` structural predicate                                             | **KEEP**              | Independent.                                                                                                                                                                                                                                                                                    |

### `tests/sn/test_legendre_moment_scattering.py` — 6 classes

| Class                                   | Line  | What it pins                                                | T.3 disposition       | Reason                                                                                                                                                                                                                                                                  |
| --------------------------------------- | ----- | ----------------------------------------------------------- | --------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `TestCapabilities`                      | 79    | `Lam.capabilities == {CAP_APPLY}`                           | **KEEP**              | Λ stays apply-only.                                                                                                                                                                                                                                                     |
| `TestPerEllBlockDiagonal`               | 93    | `skip_l0`; per-ℓ block isolation                            | **KEEP**              | Λ block structure unchanged.                                                                                                                                                                                                                                            |
| `TestPerMaterialPartition`              | 133   | Per-material dispatch via `cells_by_material`               | **KEEP**              | Per-material dispatch lives in `MaterialXSField`. Unchanged.                                                                                                                                                                                                            |
| `TestEnergyContractionDirection`        | 169   | `Σ_{g_from}` energy axis contraction direction              | **KEEP**              | Convention unchanged.                                                                                                                                                                                                                                                   |
| `TestBitIdenticalToLegacyInlinedMath`   | 240   | Λ matches legacy inline `moment @ sig_s_l[l]` at `rtol=1e-15` | **KEEP**            | Λ inner numerics unchanged (assuming T.3 wraps `LegendreMomentScattering` rather than rewriting it).                                                                                                                                                                    |
| `TestComposesUnderOperatorAlgebra`      | 290   | Λ composes via `@` with `IdentityOperator`                  | **KEEP**              | Λ stays a `LinearOperator`.                                                                                                                                                                                                                                             |

### `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py` — anisotropic curvilinear MMS gate

| Test                                              | What it verifies                                       | T.3 disposition       | Reason                                                                                                                                       |
| ------------------------------------------------- | ------------------------------------------------------ | --------------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| All tests in file                                 | O(h²) DD convergence on curvilinear P1 anisotropic MMS | **KEEP (must stay green)** | Per plan §5.2 the L1 MMS gates are ground truth. 5 pre-existing DD-regression xfails per plan §6 line 295; the set must be unchanged. |

### `tests/sn/test_mms_aniso.py` — 1D slab P1 anisotropic MMS

| Test                                                  | What it verifies                                | T.3 disposition           | Reason                                                                                                              |
| ----------------------------------------------------- | ----------------------------------------------- | ------------------------- | ------------------------------------------------------------------------------------------------------------------- |
| `test_sn_p1_aniso_mms_converges_second_order`         | O(h²) DD convergence on 1D slab P1 MMS          | **KEEP (must stay green)** | The load-bearing L1 gate that catches Pℓ scattering bugs (per the test docstring: "drops the P1 contribution, uses the wrong (2ℓ+1) factor, or transposes the scattering matrix index"). |
| `test_sn_p1_aniso_mms_source_degrades_to_p0`          | α=0 → P1 source equals P0 source                | **KEEP**                  | Pre-solve MMS source check; independent of T.3 lift.                                                                |

### `tests/sn/l1_analytical/test_kinf_homogeneous.py` — k_∞ closed-form

| Test       | What it verifies                                     | T.3 disposition           | Reason                                                                                                                                                                                                                                  |
| ---------- | ---------------------------------------------------- | ------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| all       | `k_∞ = νΣ_f/Σ_a` on homogeneous reflective           | **KEEP (must stay green)** | The structurally-independent closed-form L1 reference for eigenvalue claims. **This is the principled-equivalence pillar #2 reference** if T.3 elects nulp-relaxation. Must hold to closed-form accuracy regardless of FP reduction order. |

### Total disposition summary

- **KEEP unchanged**: 17 of 19 classes / files.
- **KEEP / RELAX-nulp conditional**: 2 of 19 (`TestAnisotropicScatteringExtraction.test_delegator_matches_operator`, `TestAlgebraicIdentity` — only if Route B in §6 Q1 is elected).
- **DELETE**: 0. T.3 is purely additive.
- **UPDATE-type**: 0. T.3 does NOT change dispatch return types.

---

## §3 New tests required

| #  | Test name                                                                     | File                                                              | Class                              | Level                | What it verifies                                                                                                                                                                                                                                                                                                                                                                                                          | Reference                                                | Pass criterion                                                                                                                              |
| -- | ----------------------------------------------------------------------------- | ----------------------------------------------------------------- | ---------------------------------- | -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | -------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| L2-1  | `test_kernel_is_sum_of_tensor_products`                                       | `tests/sn/test_scattering_operator.py`                            | `TestSumOfTensorProductsKernel`    | foundation           | T.3 verification gate. `scattering_op.kernel` is a `SumOfTensorProductsOperator` instance with `len(summands) == xs.scatter.order + 1` (one per Legendre order including ℓ=0). Mirrors T.2's `TestRankOneTensorProductKernel` pattern.                                                                                                                                                                                  | Type-introspection                                       | `isinstance(kernel, SumOfTensorProductsOperator)` and `len(kernel.summands) == op.scattering_order + 1`                                       |
| L2-2  | `test_each_summand_is_3_factor_tensor_product`                                | same                                                              | same                               | foundation           | Each summand is a 3-factor `TensorProductOperator`. The §15.2 form `Σ_{s,ℓ} & A_ℓ & G_ℓ`. **Note**: §6 Q1 may elect a 2-factor or 1-factor summand instead; this test pins whatever shape is chosen (parametrise by `expected_factor_count`).                                                                                                                                                                              | Type-introspection                                       | `len(summand.ops) == EXPECTED_FACTOR_COUNT` for all summands; each factor is a `LinearOperator`                                              |
| L2-3  | `test_assert_separable_passes`                                                | same                                                              | same                               | foundation           | The `SumOfTensorProductsOperator.assert_separable()` invariant fires without raising. First production consumer of this method — important architectural milestone per plan §5.1 line 158.                                                                                                                                                                                                                                 | Built-in invariant                                       | `kernel.assert_separable()` returns None (no raise)                                                                                          |
| L2-4  | `test_kernel_capabilities_apply_only`                                         | same                                                              | same                               | foundation           | `kernel.capabilities == {CAP_APPLY}` — intersection across summands and across factors; rank-1 in some factor of every summand precludes CAP_SOLVE; transpose deferred                                                                                                                                                                                                                                                  | Capability set algebra                                   | `kernel.capabilities == frozenset({CAP_APPLY})`                                                                                              |
| L2-5  | `test_disjoint_axes_per_summand`                                              | same                                                              | same                               | foundation           | The factors of each summand act on disjoint axes. The disjoint-axes contract is the TP-factor precondition. **CRITICAL**: §6 Q2 may surface that the per-material per-ℓ einsum is NOT disjoint-axis (touches both group AND spatial). If so, the spatial-factor is `IdentityOperator` and the group/spatial structure is in a `RankOneOperator`-style non-disjoint primitive. This test pins the chosen factorisation.    | Type-introspection on `axis` attribute                  | Each factor's `axis` attribute (or equivalent) is distinct across the 3-factor product                                                       |
| L1-1  | `test_apply_bit_identical_to_pre_t3_snapshot_angular_flux`                    | same                                                              | `TestPreT3RegressionSnapshot`      | foundation           | Per-dispatch-arm regression snapshot. Take the pre-T.3 output (captured to a fixture file BEFORE T.3 commits) on a fixed (seed, mesh, mat) triple and compare post-T.3 output. AngularFlux arm.                                                                                                                                                                                                                          | Pre-T.3 captured numerics (regression-snapshot)          | **Route A** (factor-wrap, inner R∘Λ∘M unchanged): `np.testing.assert_array_equal`. **Route B** (per-ℓ summand): `nulp=4*(L+1)` per principled-equivalence gate (see plan §7 risk 3). |
| L1-2  | `test_apply_bit_identical_to_pre_t3_snapshot_scalar_flux`                     | same                                                              | same                               | foundation           | ScalarFlux dispatch arm — P0 + (n,2n) only path. Snapshot regression.                                                                                                                                                                                                                                                                                                                                                    | Pre-T.3 captured numerics                                | Strict `np.testing.assert_array_equal` (scalar-flux path is P0+n2n which T.3 doesn't touch)                                                  |
| L1-3  | `test_apply_bit_identical_to_pre_t3_snapshot_timed_full_field`                | same                                                              | same                               | foundation           | TimedFullField composite arm. Snapshot regression on bulk + boundary parts.                                                                                                                                                                                                                                                                                                                                              | Pre-T.3 captured numerics                                | Same as L1-1 (composite delegates to AngularFlux arm internally)                                                                              |
| L1-4  | `test_build_aniso_source_bit_identical_to_pre_t3`                             | same                                                              | same                               | foundation           | The lifted `R · Λ · M` pipeline (the actual carve target) output matches pre-T.3 on fixed (P_ℓ=2, ng=4, nx×ny=4×3, asymmetric SigS) input. The single highest-leverage bit-identity check.                                                                                                                                                                                                                                | Pre-T.3 captured numerics                                | Same as L1-1                                                                                                                                  |
| L3-1  | `test_kernel_apply_linearity_in_psi`                                          | same                                                              | `TestSumOfTensorProductsKernel`    | foundation           | Linearity of the lifted kernel `apply`: `kernel.apply(α·ψ + β·φ) ≈ α·kernel.apply(ψ) + β·kernel.apply(φ)` at FP precision. Linearity SHOULD be inherited from `SumOfTensorProductsOperator.apply` automatically; pin it here as the structural guard.                                                                                                                                                                       | Algebraic identity (structural)                          | `rtol=1e-12, atol=1e-13`                                                                                                                      |
| L3-2  | `test_kernel_apply_matches_dispatch_arms`                                     | same                                                              | same                               | foundation           | `op.kernel.apply(psi.values) == op.apply(psi).values` for each of the three dispatch arms (AngularFlux, ScalarFlux, TimedFullField bulk). Single source of truth for the underlying math — mirrors T.2's `test_kernel_apply_matches_apply_dispatch`. Catches dispatch-arm drift.                                                                                                                                          | Internal consistency (kernel = single source of truth)  | `np.testing.assert_array_equal`                                                                                                              |
| L3-3  | `test_per_summand_independence`                                               | same                                                              | same                               | foundation           | Each individual summand `kernel.summands[ell].apply(...)` produces the per-ℓ in-scatter contribution AND the sum across summands equals `kernel.apply(...)`. **The per-summand verification step required by the brief.** Verify each summand individually before trusting the sum.                                                                                                                                       | Decomposition consistency                               | `rtol=1e-14, atol=1e-15` for `kernel.apply(psi) == sum(s.apply(psi) for s in kernel.summands)`                                                |
| L3-4  | `test_kernel_zero_input_zero_output`                                          | same                                                              | same                               | foundation           | `kernel.apply(zeros) = zeros` — linearity zero-guard. Trivial but pins regression where the lift accidentally introduces a constant term (e.g., a misplaced bias).                                                                                                                                                                                                                                                       | Structural (linearity)                                  | `np.testing.assert_array_equal`                                                                                                              |
| L6-1  | `test_per_material_per_ell_einsum_invariance_p1`                              | `tests/sn/test_material_xs_field.py` (new file IF not existing)   | `TestApplyLegendreScatteringMoments` | foundation           | The `MaterialXSField.apply_legendre_scattering_moments` einsum at `material_xs_field.py:515-572` is bit-identical pre/post T.3. This is the inner numerics primitive that MUST stay unchanged per plan §6. Pin on fixed (P_ℓ=1, asymmetric SigS) input.                                                                                                                                                                  | Pre-T.3 captured numerics + hand-coded einsum reference | `np.testing.assert_array_equal` (no FP reduction reorder; this primitive is the leaf)                                                        |
| L6-2  | `test_per_material_per_ell_einsum_invariance_p3`                              | same                                                              | same                               | foundation           | Higher-order variant (P_ℓ=3) of L6-1. Catches ℓ-loop drift.                                                                                                                                                                                                                                                                                                                                                              | Pre-T.3 captured numerics                                | `np.testing.assert_array_equal`                                                                                                              |
| L5-1  | `test_wave_t_t3_apply_overhead_below_5pct`                                    | `tests/sn/performance/test_wave_t_overhead.py` (extend or create) | `TestT3PerformanceRegression`      | (not vv-tagged)      | Per plan §5.4. `ScatteringOperator.apply` is in the Krylov hot path. Measure walltime delta on the 1-D slab P1 case (≥1000 iterations to amortise warmup). Pre-T.3 baseline captured to fixture; post-T.3 must be ≤ 5% slower. **CRITICAL**: this is the §7 risk-2 gate. If exceeded, investigate fusion before shipping.                                                                                                | Pre-T.3 wallclock baseline (fixture)                    | `post_t3_walltime / pre_t3_baseline <= 1.05`                                                                                                  |
| L4-1  | `test_mms_p1_aniso_l1_gate` (delegated to existing)                           | `tests/sn/test_mms_aniso.py` (existing)                           | (n/a — module-level)               | l1                   | The pre-existing P1 aniso MMS convergence test. **Must stay green at the same xfail set**. Plan §5.2 — the L1 MMS gates are ground truth.                                                                                                                                                                                                                                                                                  | MMS pillar (structurally independent)                   | All pre-existing passing tests still pass; xfail set unchanged                                                                                |
| L4-2  | `test_mms_curvilinear_aniso_l1_gate` (delegated to existing)                  | `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py` | (existing)                       | l1                   | Curvilinear P1 aniso MMS. Same as L4-1 but on curvilinear geometry. Plan §5.2 + §6 line 295: 5 pre-existing DD-regression xfails must stay at the same set.                                                                                                                                                                                                                                                              | MMS pillar (curvilinear)                                | Same xfail set                                                                                                                                |
| L4-3  | `test_kinf_homogeneous_l1_gate` (delegated to existing)                       | `tests/sn/l1_analytical/test_kinf_homogeneous.py`                 | (existing)                         | l1                   | `k_∞ = νΣ_f / Σ_a` closed-form on homogeneous reflective. **The structurally-independent closed-form reference for the principled-equivalence three-criteria gate** (per `vv-principles` §"Bit-identity vs principled-equivalence"). If T.3 elects nulp-relaxation, this is the pillar-#2 evidence that the new value is verified against an independent reference, NOT just close to the old value.                  | Closed-form pillar (homogeneous infinite medium)        | Existing assertion tolerance                                                                                                                  |

### Test count summary

- **New foundation tests (L0 / L2 / L3 / L6 type-and-structure)**: 12.
- **New performance gates (L5)**: 1.
- **Existing L1 gates that must stay green (L4)**: 3 files (P1 aniso slab MMS; curvilinear P1 aniso MMS; k_∞ homogeneous).
- **Total new test files**: 0-1 (extend existing `test_scattering_operator.py`; possibly add `tests/sn/test_material_xs_field.py` if not already present).
- **Pre-T.3 snapshot fixtures required**: 1 file under `tests/sn/_fixtures/wave_t_t3/pre_t3_snapshots.npz` carrying captures for L1-1, L1-2, L1-3, L1-4, L6-1, L6-2.

---

## §4 Order of T.3 implementation gates

T.3 is a one-shot lift of a single operator family, but the **substep
ordering matters for catching bugs early**. Each substep gate must
pass before the next substep starts.

### Substep T.3a — Capture pre-T.3 snapshots

**Action.** Before any T.3 code change, run a one-shot script that
exercises each dispatch arm on fixed (seed, mesh, mat) triples and
writes outputs to `tests/sn/_fixtures/wave_t_t3/pre_t3_snapshots.npz`.

**Gate.** Snapshot file exists; `tests/sn/_fixtures/wave_t_t3/`
schema validated (one ndarray per dispatch arm × P_ℓ_order).

**Why first.** Per plan §5.1 — bit-identity is the contract. Cannot
demonstrate bit-identity post-T.3 without the pre-T.3 numbers in
hand. ALSO: the L1-1..L1-4 tests are useless WITHOUT this fixture.

### Substep T.3b — Add SOTP kernel construction (NO call-site rewire)

**Action.** Add a `kernel` cached_property on `ScatteringOperator`
that constructs the §15.2 `SumOfTensorProductsOperator` with
`xs.scatter.order + 1` summands. **Do not rewire** any `.apply`
arm yet.

**Gate.**
- L2-1: kernel structure (SOTP with correct summand count).
- L2-2: each summand is the chosen factor count.
- L2-3: `assert_separable` passes.
- L2-4: capabilities = {CAP_APPLY}.
- L2-5: disjoint axes per summand.

**Why this gate before rewire.** Catches §6 Q1 / Q2 design failures
EARLY. If the disjoint-axes contract fails, the factor design is
wrong; back out and redesign without contaminating dispatch arms.

### Substep T.3c — Wire `build_aniso_source` to use the kernel

**Action.** Rewire `build_aniso_source` to call `self.kernel.apply(...)`
instead of the inline `R · Λ · M` chain at `scattering.py:625-651`.
Keep the producer-side `/sum_w` projection at the apply boundary
(Pattern 7 invariant) — that lives OUTSIDE the kernel.

**Gate.**
- L1-4: `build_aniso_source` bit-identical (Route A) or nulp=4*(L+1) (Route B).
- L3-1: linearity holds.
- L3-3: per-summand independence (each summand produces a coherent per-ℓ block).
- L3-4: zero input → zero output.

**Why this gate before dispatch rewire.** `build_aniso_source` is
the inner numerics; if its post-T.3 output drifts unacceptably, all
three dispatch arms inherit the drift.

### Substep T.3d — Wire each dispatch arm to consume the kernel

**Action.** AngularFlux arm first (the legacy hot path). Verify
gate. Then ScalarFlux arm. Verify gate. Then TimedFullField arm.
Verify gate. **One arm at a time** — preserves the ability to
bisect drift to the offending arm.

**Gate per arm.**
- L1-1 / L1-2 / L1-3: snapshot regression bit-identical (Route A) or nulp.
- L3-2: `kernel.apply == apply(typed_input).values` for the arm.
- Existing `TestApplySemantics` + `TestProducerSideNormalisation` + `TestCompositeInvariants` arm-relevant tests stay green.

### Substep T.3e — Per-material einsum invariance

**Action.** No code change here — this is a verification step.

**Gate.**
- L6-1, L6-2: `MaterialXSField.apply_legendre_scattering_moments`
  output unchanged.

**Why this is a separate substep.** Per plan §6 T.3 paragraph 3 the
per-material per-ℓ einsum is "Bit-identical preserved". This gate
explicitly verifies that contract — defense against a T.3 author
accidentally modernising `material_xs_field.py` while they're in the
file.

### Substep T.3f — L1 MMS + closed-form L1 gates

**Action.** Run the full L1 suite.

**Gate.**
- L4-1: 1D slab P1 aniso MMS O(h²) convergence holds.
- L4-2: curvilinear P1 aniso MMS DD convergence — same xfail set.
- L4-3: k_∞ homogeneous matches closed form within existing tolerance.

**Why last (but mandatory).** L1 MMS is the structurally-independent
ground. If L1 MMS fails after L2 + L1 snapshot tests pass, the
snapshot itself was wrong (a coincidence — both old and new code
agree on a wrong value). Per plan §5.2 line 162: "the L1 MMS gates
are the ground truth. They must stay green at every Wave T commit.
If a substep breaks an L1 gate, the substep is wrong (regardless
of bit-identity claims — bit-identity is necessary but not
sufficient when reductions reorder)."

### Substep T.3g — Performance regression

**Action.** Run the 1D slab Krylov benchmark.

**Gate.**
- L5-1: ≤ 5% slowdown.

**Why last.** Optimisation is a separate concern from correctness;
shipping a correct but slow T.3 is a regression but not a
disqualifier (per plan §7 risk 2 mitigation).

### Master gate: COMMIT only when ALL of {T.3a-T.3g} pass

Per `feedback_no_method_implementer_for_surgical_carves`: each
substep is a separate commit IFF the previous gate passes
cleanly. If a gate fails, fix in place; do not commit a failing
substep.

---

## §5 Risks and how each test catches them

Cross-referenced from plan §7 risk register row 4 (T.3 is the
highest-leverage scattering rewire) and the test-architect
`numerical-bug-signatures` mapping.

### Risk R1 — Reduction-order drift exceeds principled-equivalence bound

**Mechanism.** T.3 Route B (per-ℓ summand reduction) reorders the
`Σ_ℓ` outer reduction. Drift bounded by `(order+1) × ULP ≈ 4 ULP`
for typical L=3.

**Caught by.** L1-1 / L1-4 snapshot tests fail under
`np.testing.assert_array_equal` (Route A) or `nulp ≤ 4*(L+1)`
(Route B). L3-3 (per-summand independence) confirms the drift is
indeed at the outer sum, not inside a summand.

**Mitigation.** Default to Route A (preserve inner reduction
order; lift only the outer apply shape). Switch to Route B only if
the factor design forces it.

### Risk R2 — Disjoint-axes contract violated

**Mechanism.** §6 Q2 — the per-material per-ℓ einsum touches BOTH
the group axis AND the spatial axis via `cells_by_material`
indexing. A naive 3-factor split into `(SigmaMoment, AngularMoment,
GroupScatteringMatrix)` may fail because `SigmaMoment`'s
per-material structure couples cells to materials.

**Caught by.** L2-5 (disjoint-axes test) fails at substep T.3b. Or
`TensorProductOperator` constructor raises at kernel construction.

**Mitigation.** Use `RankOneOperator`-style or
`SigmaMoment`-as-per-cell-DiagonalOperator factor (the per-material
structure is per-cell scalar with implicit material lookup, similar
to T.2's fission factor). See §6 Q2.

### Risk R3 — Λ inner numerics accidentally modernised

**Mechanism.** T.3 author refactors
`MaterialXSField.apply_legendre_scattering_moments` while in the
file, breaking bit-identity at the leaf level.

**Caught by.** L6-1, L6-2 (per-material einsum invariance) at
substep T.3e.

**Mitigation.** Substep T.3e is explicitly defensive against this.
The T.3 author MUST NOT touch `material_xs_field.py` lines 515-572.

### Risk R4 — Λ-ZERO ℓ=0 block handling differs

**Mechanism.** `LegendreMomentScattering` defaults to `skip_l0=True`
in `build_aniso_source` (the P0 path is the separate fast path
`add_iso_source` + `add_n2n_source`). The lifted SOTP has `order+1`
summands — one for EACH ℓ including ℓ=0. The ℓ=0 summand must
EITHER (a) be the no-op identity OR (b) emit zero. If T.3
accidentally double-counts P0 (once via `add_iso_source`, once via
ℓ=0 summand), `apply` doubles.

**Caught by.** L1-1 snapshot fails dramatically (2× the expected
in-scatter source). Also `TestApplySemantics::test_apply_isotropic_flux_p0_only`
explicit hand-calc.

**Mitigation.** The ℓ=0 summand structure design must be principled
— either skip ℓ=0 from the SOTP summand list (making `len(summands)
== order`), or include it as zero / identity. Pre-T.3 spec the
correct design BEFORE substep T.3b. See §6 Q3.

### Risk R5 — Producer-side `/sum_w` projection drift (Pattern 7 violation)

**Mechanism.** The `/sum_w` projection at the apply boundary
(`scattering.py:929-930, 986-987`) is the Pattern 7 producer-side
normalisation. If T.3 moves this INTO the kernel (under one of the
summand factors), the convention is now distributed across factors
and consumers may re-introduce it.

**Caught by.** `TestProducerSideNormalisation::test_typed_apply_returns_per_ordinate_already_normalised`
(existing, line 405) fails because the consumer-facing magnitude
changes.

**Mitigation.** Keep `/sum_w` OUTSIDE the kernel. The kernel
returns moment-space output; the apply boundary divides by `sum_w`
explicitly. Document in T.3 commit message.

### Risk R6 — Scattering matrix transpose convention (ERR-002 family)

**Mechanism.** Per `numerical-bug-signatures` Signature 3, ORPHEUS
convention is `Mixture.SigS[l][g_from, g_to]`. Scattering source
is `Σ_g' SigS[g', g] · φ_g'` — implemented as `phi @ SigS`
(NOT `phi @ SigS^T`). A vectorised rewrite in T.3 that flips this
silently double-transposes.

**Caught by.**
- `TestEnergyContractionDirection` (existing in
  `test_legendre_moment_scattering.py:169`) — pins
  `Σ_{g_from}` direction.
- `TestPerMaterialPartition::test_only_mat0_cells_scattered_with_sig_s_mat0`
  (existing) — symmetric mat-0 vs mat-1 with asymmetric SigS.
- L4-1 P1 MMS — 2-group asymmetric SigS exercises the convention.

**Mitigation.** Asymmetric multi-group test data is mandatory in
T.3 fixtures. The pre-T.3 snapshot fixture MUST use a 4-group
asymmetric SigS (not the symmetric 1-group case).

### Risk R7 — Performance regression > 5%

**Mechanism.** SOTP.apply is a Python-level fold over
`(order+1)` summands; each summand is a 3-factor TP fold over 3
operators. For P_ℓ=3: 4 × 3 = 12 Python `op.apply(...)` calls per
matvec vs ~3 numpy einsums pre-T.3.

**Caught by.** L5-1 perf gate at substep T.3g.

**Mitigation.** Cache the kernel construction (`cached_property`).
If perf-regression triggers, investigate factor fusion (collapse
identity factors at construction time) — `_build` flattening
already exists for TP. SOTP-level fusion deferred per plan §4 line
142.

### Risk R8 — `LegendreMomentScattering` apply signature drift

**Mechanism.** T.3 wraps `LegendreMomentScattering` as the
`Σ_{s,ℓ}` factor of each summand. If the wrapper passes moments
in shape `(ng, nx, ny)` (one ℓ block) instead of the full
`(L+1, 2L+1, ng, nx, ny)` (all blocks), the ℓ-loop inside
`apply_legendre_scattering_moments` (line 558-571) will iterate
over zero blocks.

**Caught by.** L1-4 build_aniso_source snapshot regression fails
(output is all zeros).

**Mitigation.** Document the moment-tensor convention at the
`SigmaMoment` factor's apply boundary. Verify shape contract in
substep T.3b.

---

## §6 Architectural questions to surface BEFORE T.3 implementation

These questions MUST be resolved (with user steering) before
substep T.3b code is written. Each question carries a default
recommendation; the user may override.

### Q1 — Is the SOTP outer reduction over ℓ, or is `R · Λ · M` one composition?

The §15.2 form `Σ_ℓ (Σ_{s,ℓ} ⊗ A_ℓ ⊗ G_ℓ)` has an outer `Σ_ℓ`.
The current `build_aniso_source` body is `R.apply(Λ.apply(M.apply(ψ)))
/ sum_w` — a single composed pipeline; the outer `Σ_ℓ` is hidden
INSIDE `LegendreMomentScattering.apply` as the for-loop at
`material_xs_field.py:558-571`.

**Two routes**:

- **Route A (composition-shaped)**. The lifted kernel is ONE
  `TensorProductOperator` with three factors `(R, Λ, M)` — not a
  `SumOfTensorProductsOperator`. **This violates the plan §6 T.3
  target form** which explicitly names `SumOfTensorProductsOperator`,
  but preserves strict bit-identity at the matvec level (no outer
  reduction reorder).

- **Route B (sum-of-products-shaped)**. The lifted kernel is a
  `SumOfTensorProductsOperator` with `order+1` summands; each
  summand is `(SigmaMoment_ℓ, AngularMomentOperator_ℓ,
  GroupScatteringMatrix_ℓ)`. **This matches the plan §6 T.3 form**
  AND grand-report §15.2. Cost: outer `Σ_ℓ` reorders the
  reduction; principled-equivalence three-criteria gate applies.

**Default recommendation: ROUTE B (sum-of-products)**, per plan
§6 T.3 explicit target AND grand-report §15.2 AND the §5.1 plan
test "isinstance(S, SumOfTensorProductsOperator)" line 156. The
nulp-relaxation cost is bounded (`nulp ≤ 4*(order+1)`) and
justified by the principled-equivalence three-criteria gate
(`vv-principles` §"Bit-identity vs principled-equivalence"):

1. **Principled at every step**: each summand `Σ_{s,ℓ} ⊗ A_ℓ ⊗
   G_ℓ` IS the per-ℓ in-scatter contribution — a named
   reactor-physics quantity (the ℓ-th angular moment scattering
   source per energy group per cell). Not an unnamed intermediate.
2. **Structurally-independent reference**: L4-3
   (`k_∞ = νΣ_f / Σ_a` closed-form) on homogeneous reflective.
   The k_∞ value is independent of the scattering reduction order.
3. **FP-non-associativity, dimensionally explainable**: drift ≤
   `(order+1) × ULP` ≈ `4 ULP` for L=3.

**Action.** Set `EXPECTED_FACTOR_COUNT=3` in L2-2; set Route B
nulp tolerance in L1-1 / L1-4.

### Q2 — How does the 3-factor TP satisfy the disjoint-axes contract?

Per §15.2 the three factors are named `Σ_{s,ℓ}` (group axis),
`A_ℓ` (angular axis), `G_ℓ` (third axis — spatial? per-material?).
But the per-material per-ℓ einsum at `material_xs_field.py:515-572`
COUPLES group AND spatial axes (via `cells_by_material`
indexing).

The disjoint-axes contract requires each factor to act on a
single axis and broadcast on the others. The per-material einsum
violates this naively because the cross-section depends on the
cell material.

**Three sub-options**:

- **(a) `SigmaMoment_ℓ = PerCellPerGroupScalarOperator`**. A
  per-cell scalar that has the right `Σ_{s,ℓ}^{m(cell)}[g_from,g_to]`
  baked in. The factor's `axis` is the group axis; the per-material
  lookup is internal. This MATCHES the
  `MaterialXSField.apply_legendre_scattering_moments` body — the
  whole call IS this factor.

- **(b) `RankOneOperator`-style coupled factor**. Per T.2's
  fission pattern, use a non-disjoint primitive
  (`RankOneOperator(left, right, axis=0)`) as one of the TP
  factors. The factor advertises a single `axis` but its action
  internally couples that axis with spatial. The L2-5 disjoint-axes
  test then verifies the factor's ADVERTISED axis is distinct from
  the other factors' axes (not the structural coupling).

- **(c) Wrap `LegendreMomentScattering` directly as the
  `Σ_{s,ℓ}` factor**. The Λ class already exists and consumes the
  full moment tensor `(L+1, 2L+1, ng, nx, ny)`. The wrapper
  factor declares `axis=2` (the group axis) by convention; the
  factor's apply IS Λ.apply.

**Default recommendation: OPTION (c)**, using
`LegendreMomentScattering` directly as the per-summand factor
that handles the per-material group-axis structure. The
`AngularMomentOperator_ℓ` factor handles the angular (moment-tensor
ℓ,m → ordinate n) projection — that's the per-ℓ block of `R · M`.
The third factor is `IdentityOperator` for the spatial axis. This
makes each summand effectively 2 non-identity factors + 1
identity, mirroring T.2's pattern.

But this leaves a problem: `LegendreMomentScattering.apply`
already iterates over ALL ℓ blocks internally. The per-summand
form needs a per-ℓ Λ. Solution: pass `L=ell, skip_l0=False` per
summand (one ℓ block per summand instance).

**Action.** Document the chosen option in the T.3 commit message
and pin the factor structure in L2-2 / L2-5.

### Q3 — How is the ℓ=0 summand handled?

`build_aniso_source` skips ℓ=0 (`Lam = LegendreMomentScattering(...,
skip_l0=True)`). The P0 path is the fast-path
`add_iso_source` + `add_n2n_source` which operates on `(ng, nx, ny)`
scalar flux, NOT the full moment tensor. The `Σ_ℓ` formal sum
across `ℓ ∈ [0, order]` thus has either:

- **(α)** an explicitly-zero ℓ=0 summand (a no-op identity term
  multiplied by zero), keeping `len(summands) == order + 1`. The
  P0 fast path stays in `apply` outside the kernel.
- **(β)** a summand list of length `order` (ℓ from 1 to order),
  not including ℓ=0. P0 stays in the fast path.
- **(γ)** the ℓ=0 summand includes the P0 + (n,2n) contribution.
  P0 fast path retires; everything routes through the kernel.

**Default recommendation: OPTION (β)**, `len(summands) == order`,
P0 + (n,2n) stays in the existing fast path
(`add_iso_source` + `add_n2n_source`). The fast path uses
`ScalarFlux` (the integrated angular flux), NOT the moment tensor
— architecturally different from ℓ≥1 which needs the full
`(L+1, 2L+1, ng, nx, ny)` moment tensor. Forcing P0 into the SOTP
form requires (a) reshaping the iso source as a moment-tensor
ℓ=0 block; (b) running it through `R` (which produces `N × ng × nx
× ny`); (c) projecting back via `/sum_w`. That's algebraically
correct but procedurally wasteful — the fast path exists precisely
because P0 doesn't need the moment tensor.

**Implication for L2-1.** Set
`expected_summand_count = op.scattering_order` (not `+1`).
**Update the plan §6 T.3 paragraph 1** which states `xs.scatter.order + 1`
— that's the grand-report aspirational form including P0; the
ORPHEUS production form excludes P0 (fast path).

**Action.** Confirm with user: is the grand-report `+1` form a
hard target, or is the production-form `order` (P0 in fast path)
acceptable?

### Q4 — Does the kernel surface `apply_transpose`?

`ScatteringOperator.capabilities` = `{CAP_APPLY}` today (per the
class docstring lines 22-30 and the `LegendreMomentScattering`
class lines 200-202). Transpose is "not currently consumed by any
ORPHEUS solver and is deferred."

`SumOfTensorProductsOperator.apply_transpose` IS implemented at
`operator.py:1247-1256`, propagating across summands and factors
via `apply_transpose`. If T.3's kernel summands ALL advertise
`CAP_APPLY_TRANSPOSE`, the kernel inherits it.

**Default recommendation. Keep `{CAP_APPLY}` only.** Set kernel
factors' capabilities accordingly. T.3 is not the wave that
introduces the adjoint surface; that's a follow-up.

**Action.** L2-4 pins `kernel.capabilities == frozenset({CAP_APPLY})`.

### Q5 — Where does the producer-side `/sum_w` projection live?

Plan §6 T.3 paragraph 4: "Bit-identical preserved." The `/sum_w`
factor at the apply boundary (`scattering.py:929-930, 986-987`)
is the Pattern 7 producer-side normalisation. **It must stay at
the apply boundary, OUTSIDE the kernel.**

**Action.** Document in the T.3 commit message that the kernel
produces moment-space output and the `/sum_w` projection happens
AT the apply boundary, NOT inside any summand or factor. Pin via
`TestProducerSideNormalisation` (existing, line 405).

### Q6 — What's the principled fallback if disjoint-axes contract cannot be satisfied?

Per plan §7 risk 1 (which is about T.4 streaming, but applies
generically): "Acceptable fallback: factor what's factorable,
document non-factorisable pieces, ensure the algebra still composes
via `OperatorSum`. The architectural commitment is ONE algebraic
FORM, not that every factor is a clean `TensorProductOperator`."

If Q2 option (c) fails the disjoint-axes contract (because
`LegendreMomentScattering` couples group + spatial through
per-material indexing), the fallback is:

- The "summand" is a custom `LinearOperator` (NOT a
  `TensorProductOperator`) carrying the per-ℓ algebra.
- The kernel is `OperatorSum` over those custom summands, NOT
  `SumOfTensorProductsOperator`.
- The L2-1 "isinstance(kernel, SumOfTensorProductsOperator)" test
  is then INAPPLICABLE; the test changes to
  "isinstance(kernel, OperatorSum) AND each summand is a per-ℓ
  scattering action".

**Action.** If the disjoint-axes contract fails during T.3b,
SURFACE TO USER immediately. Do NOT silently downgrade to
`OperatorSum`. The grand-report §15.2 form explicitly says SOTP;
diverging from it is a plan-level decision.

---

## §7 Pre-existing test inventory cross-check

Already covered in §2. The total disposition:

- **KEEP unchanged**: ~50 tests.
- **KEEP / RELAX-nulp conditional on §6 Q1 route**: 2 tests
  (`TestAnisotropicScatteringExtraction.test_delegator_matches_operator`
  and `TestAlgebraicIdentity._check_identity`).
- **DELETE**: 0.
- **UPDATE-type**: 0.

No type assertion drift because T.3 does NOT change dispatch
return types (AngularFlux→PerOrdinateSource, ScalarFlux→IsotropicSource,
TimedFullField→TimedFullField, unsupported→TypeError).

---

## §8 Test-architect self-assessment

### Pillar gate (per `vv-principles` §1.5)

| Test row | Claim layer       | Pillar                       | Structural independence?                                                     |
| -------- | ----------------- | ---------------------------- | ----------------------------------------------------------------------------- |
| L2-1..5  | Software invariant (type / structure)  | n/a (foundation)             | Pure type introspection                                                       |
| L1-1..4  | Convergence (in the sense of bit-identity to pre-T.3) | Regression-snapshot | Pre-T.3 numerics ARE the reference; pillar #2 (k_∞ closed-form) backs nulp |
| L3-1..4  | Algebraic identity (linearity / decomposition)        | Closed-form (linearity is a closed-form claim) | Structural — pure math |
| L6-1..2  | Bit-identity (inner numerics)         | Regression-snapshot          | The L6 contract is the leaf primitive; one structural angle                  |
| L5-1     | Performance (not vv-tagged)           | n/a                          | n/a                                                                            |
| L4-1..3  | Eigenvalue / flux-shape / convergence-order | MMS + closed-form        | MMS (P1 aniso slab + curvilinear); closed-form (k_∞)                          |

**MMS does NOT prove eigenvalues.** L4-1 / L4-2 are
**flux-shape + convergence-order claims**. L4-3 (k_∞ homogeneous)
is the **eigenvalue claim** via closed-form.

### Anti-pattern audit

Pre-flight on each `vv-principles` §0 anti-pattern:

1. ✅ NOT claiming verification on L4 alone — L1 + L4 stack.
2. ✅ NOT asserting `np.allclose` against another solver — comparing to
   pre-T.3 captured numerics (snapshot) + closed-form k_∞.
3. ✅ NOT accepting 1-group as verification — 2G/4G test data
   mandated in L1-1 / L1-4 / L4-1 / L4-2.
4. ✅ NOT homogeneous-only — L4-1 includes the manufactured-source
   heterogeneity; the 4-group asymmetric SigS in L1 fixtures.
5. ✅ NOT trusting convergence rate alone — L4-1 verifies the
   converged-to value via the MMS exact solution.
6. ✅ NOT trusting an untraced reference — pre-T.3 snapshot is
   regression / control-on-self, the L1 chain terminates in
   closed-form k_∞ + MMS exact solution.
7. ✅ NOT trusting "two derivations agree" — pre-T.3 vs post-T.3 is
   bit-identity at the IMPLEMENTATION level (Route A) or
   nulp-bounded (Route B) backed by independent closed-form.
8. ✅ NOT particle balance as L0 — per-summand independence (L3-3)
   replaces telescoping-sum traps.
9. ✅ NOT conflating validation with verification — T.3 is a pure
   refactor; no equation change; pure V&V.
10. ✅ NOT "reasonable numbers" — every term verified with sign +
    magnitude via L0-style structural tests.

### Multi-group + heterogeneous + mesh-refined?

- **Multi-group**: pre-T.3 snapshot fixtures use 4-group asymmetric
  SigS (Risk R6 mitigation).
- **Heterogeneous**: 2-material 2D mesh in `solver_2g_p1_n2n`
  fixture (existing).
- **Mesh-refined**: L4-1 (`n_cells = [20, 40, 80, 160]`) is the
  convergence ladder.

All three mandatory dimensions covered.

### Self-improvement notes

No new ERR-NNN failure modes introduced; no new entry needed in
`numerical-bug-signatures` (T.3 is structurally a Pattern 1 lift,
not a new bug class). The risk register R1-R8 in §5 above
generalises across all Wave T substeps and SHOULD be the canonical
risk template for Wave T T.4 dispatch when that wave starts.

---

## §9 Action items (sequence-coupled)

Per `subagent-handoff-protocol` and `feedback_no_method_implementer_for_surgical_carves`:

1. **Main agent + user**: review this spec for §6 Q1-Q6 answers.
2. **Main agent**: write the pre-T.3 snapshot capture script (one-shot
   in `tools/wave_t/capture_pre_t3_snapshots.py` or inline). Run
   once; commit fixtures.
3. **Main agent direct authorship** (no method-implementer): execute
   substeps T.3a → T.3g per §4 in sequence; one commit per substep.
4. **qa dispatch** (post-substep T.3d): review the dispatch-arm
   rewires; check Pattern 7 producer-side `/sum_w` integrity;
   check ERR-002 transpose convention preservation in the lifted
   form.
5. **Documentation update** (substep T.3e or in T.5 wave close-out):
   Sphinx theory page `docs/theory/operator_algebra.rst` adds the
   §15.2-named instance to the active production tensor-network
   table.

---

## Pointers

- **Plan**: `.claude/plans/wave_t_tensor_network.md` §6 step T.3 + §5
  verification strategy + §7 risk register row 4.
- **Grand report**: §15.2 lines 2046-2086 (sum-of-tensor-products
  scattering form).
- **Production code**:
  - `orpheus/sn/scattering.py:151-271` (LegendreMomentScattering).
  - `orpheus/sn/scattering.py:274-987` (ScatteringOperator).
  - `orpheus/sn/scattering.py:527-651` (build_aniso_source — the
    lift target).
  - `orpheus/sn/material_xs_field.py:515-572`
    (apply_legendre_scattering_moments — bit-identity preserved).
  - `orpheus/numerics/projection.py:240-599`
    (MomentProjection, HarmonicMomentReconstruction).
- **L1 primitives**:
  - `orpheus/numerics/operator.py:1059-1177` (TensorProductOperator).
  - `orpheus/numerics/operator.py:1180-1270`
    (SumOfTensorProductsOperator + assert_separable).
  - `orpheus/numerics/operator.py:1391-1510` (RankOneOperator —
    T.2's precedent for non-disjoint TP factors).
- **T.2 precedent**:
  `tests/sn/test_fission_operator.py:339-407`
  (TestRankOneTensorProductKernel — the test class pattern to
  mirror for L2-1..5).
- **L1 MMS gates**:
  - `tests/sn/test_mms_aniso.py` (1D slab P1 aniso, must stay
    green).
  - `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
    (curvilinear, same xfail set).
  - `tests/sn/l1_analytical/test_kinf_homogeneous.py` (closed-form
    ground for principled-equivalence).
- **`vv-principles` §"Bit-identity vs principled-equivalence"**:
  the three-criteria gate for Route B (nulp-relaxation) — every
  criterion above answered in §6 Q1 default recommendation.
- **`coding-elegance` Pattern 7**: producer-side normalisation —
  the `/sum_w` invariant that must stay at the apply boundary.
- **`numerical-bug-signatures` Signature 3**: ERR-002 SigS
  transpose convention — guarded by §5 R6.
