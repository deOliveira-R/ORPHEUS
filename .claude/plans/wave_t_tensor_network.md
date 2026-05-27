# Wave T — Tensor-Network Operator Algebra

**Branch:** `refactor/moment-space-and-layering`
**Worktree:** `.claude/worktrees/moment-space-and-layering/`
**Phase status:** PENDING. Lands AFTER Depth B (`depth_b_field_on_function_space.md`) completes (D-K shim retirement). Precedes parent-plan step P3.4 (Problem/Solver split).

**Date:** 2026-05-27. Plan author: main agent + explorer audit (`.claude/agent-memory/explorer/` — `function_space_typed_field_audit.md` for Depth B background; the conversation that motivated Wave T is recorded in the Depth B plan §11.3).

**Status:** STUB. The architectural commitment is settled (see §1); detailed designs for each substep land when each substep starts (depth-on-demand per `[[feedback_no_method_implementer_for_surgical_carves]]` — main-agent direct authorship with turn-by-turn user steering).

---

## 0. Pickup checklist (read first)

If you are picking this plan up in a fresh session:

1. **Read this plan top-to-bottom.** No section is optional.
2. **Verify Depth B completed**: `git log --oneline -25` should show the D-A through D-K commits AND the `1e0bb98` Field ABC plus `53bc986` ScalarFlux migration as landed. The parent plan's `moment_space_and_layering_plan.md` step P3.3 should be marked LANDED (re-scoped through Depth B).
3. **Read these grand-report sections**:
   - §15 lines 2003-2086 — "The tensor section": V = X ⊗ Ω ⊗ G; native operator form = sum of tensor products.
   - §15.1 lines 2031-2045 — streaming as `L = Σ_axis (D_axis ⊗ Ω_axis ⊗ I_g)`.
   - §15.2 lines 2046-2086 — scattering as `S = Σ_ℓ (Σ_{s,ℓ} ⊗ A_ℓ ⊗ G_ℓ)`.
   - §16A.10 lines 3142-3197 — BC as tensor network; `BoundaryGeometryTensor`, `BoundaryResponseTensor`; canonical separable form `B = G_patch ⊗ K_omega ⊗ K_g`.
   - §35 line 5675 — 20 commandments: "Represent multigroup scattering as a sum of tensor products."
   - North-star line 5697: "tensor products are foundational".
4. **Read the explorer's two audit memos** that motivated Wave T (referenced in Depth B plan §11.3 footnotes; if absent, regenerate via explorer dispatches against the grand report).
5. **Read the existing L1 primitives**:
   - `orpheus/numerics/operator.py:1004-1122` — `TensorProductOperator` (definition, `__and__` dunder, `apply` fold).
   - `orpheus/numerics/operator.py:1125-1215` — `SumOfTensorProductsOperator`.
   - `orpheus/numerics/operator.py:1203-1215` — `assert_separable()` (currently shipped but never called).
   - `orpheus/numerics/space.py` — `TensorProductSpace` (lands in Depth B D-B).
6. **Pick up at the leftmost incomplete step in §6** below.

The `feedback_no_method_implementer_for_surgical_carves` rule applies: this is the main agent's work with turn-by-turn user steering. Do NOT batch via method-implementer.

---

## 1. The principle

**The deterministic multigroup neutron transport problem already lives on a tensor-product space.** Per grand report §15 line 2003-2019:

> *"The deterministic multigroup problem already lives on a tensor-product space: V = X ⊗ Ω ⊗ G. The stochastic problem extends this: V = X ⊗ Ω ⊗ G ⊗ Ξ_1 ⊗ ... ⊗ Ξ_d. The native operator form is often a sum of tensor products: A = Σ_k A_x^(k) ⊗ A_ω^(k) ⊗ A_g^(k) ⊗ A_ξ^(k)."*

The codebase today implements this tensor-product structure via flat-axis numpy with implicit broadcasting — `np.einsum("gxy,gxy->xy", chi, phi)` for fission, `np.take(x, perm, axis=0)` for specular BC, etc. Numpy broadcasting IS the implementation of `... ⊗ I_other_axes` — but the *type system* doesn't know, the *operator algebra* doesn't know, and the algebraic identities `(A ⊗ B)^† = A^† ⊗ B^†` and `(A ⊗ B) ∘ (C ⊗ D) = (A ∘ C) ⊗ (B ∘ D)` (documented at `operator.py:1024-1033`) cannot fire because no production operator is ever instantiated as `TensorProductOperator`.

**Wave T closes the gap.** The L1 primitives (`TensorProductOperator`, `SumOfTensorProductsOperator`, `assert_separable`) shipped in Wave 0 of the parent plan. Today they have ZERO production consumers. Wave T rewires the boundary realizers, the fission operator, the scattering operator, and the modern streaming operator to use the algebra natively. Wave T is the consumer side of the Wave-0 infrastructure.

### 1.1 The "do it now" decision

Wave T was not originally in the moment-space-and-layering plan. The decision to insert it between Depth B and P3.4 was made on 2026-05-27 in conversation, on these grounds (Depth B plan §11.3):

1. **The infrastructure is already shipped, unused.** Every day the procedural single-axis path entrenches.
2. **The status quo is a slow Cardinal Rule 2 violation.** Numpy axis-broadcast is functionally equivalent to `P ⊗ I ⊗ I`, but the type system claims `PermutationOperator(axis=0)`. Three layers (math, numpy, types) tell three different stories.
3. **The grand report names it explicitly.** §15-§16A spell out the tensor-product structure as the *native* algebraic shape, not a performance choice.
4. **P3.4 (Problem/Solver split) needs operators with clear tensor structure** to express problem-level invariants. Doing Wave T before P3.4 unlocks the right shape; doing P3.4 first forces it to consume the flat single-axis legacy.

### 1.2 The "two views of the connection-coefficient operator" insight

A key architectural finding from the Wave T scoping (Depth B plan §11.3): the curvilinear streaming operator's `(1−μ²)/r · ∂_μ` (sphere) and `−(1/r) · ∂_φ(ξ·)` (cylinder) terms are **the same connection-coefficient operator on SO(3) viewed in two coordinate charts** (`orpheus/geometry/reduced_operator.py:1-30`).

The geometric data — face areas, volumes, ΔA/w coefficients, M-M α dome — is encapsulated in `StreamingTerms` (`orpheus/geometry/reduced_operator.py:151`) and the angular closure strategy in `PoleAngularClosure` (`orpheus/sn/spatial/pole_angular_closure.py:191`). The cell-balance algebra is already **geometry-blind by data** (`orpheus/sn/spatial/cell_balance.py:120-198`):

```
denom[g, n]          = streaming_denom_term + angular_denom_term + collision_denom_term
                     = 2|μ_n|·A_down[n]     + (ΔA/w)·c_out[n]     + Σ_t,g·V

numer_upstream[g, n] = spatial_upstream_term + angular_numer_upstream
```

Three additive, structurally-independent terms. For slab, the `angular_*` terms are zero (Identity closure). For sphere/cylinder, the M-M `PoleAngularClosure` populates them. Slab and curvilinear share ONE balance equation — they differ only in DATA.

**Consequence for T.4:** the universal tensor-product decomposition `L = L_spatial + L_angular_redist` works for all geometries. Slab is the degenerate case where `L_angular_redist = ZeroOperator`. Curvilinear activates the connection-coefficient operator. The user's framing — "two views of the same operator" — is literally the streaming-view (face-area-to-volume ratios) and angular-redistribution-view (connection coefficients) of one underlying geometric primitive.

### 1.3 Expected complications

The user has flagged (2026-05-27): *"we'll probably run into problems when we get to T.4, because we did this other times and making things truly geometric agnostic always brings new complications. Let's face difficulties as they come."*

Known unknowns:
- Curvilinear `R_angular` may not factor as cleanly as the Cartesian streaming term — the M-M half-grid is a sequential sweep coupling, not a clean diagonal. Whether it's representable as a `TensorProductOperator` factor (vs. requiring a bespoke `AngularSweepOperator`) is open.
- The L+C composition's `.solve` path (`InvertibleOperator.solve` via `transport_sweep`) is a procedural algorithm, not factored. T.4 must NOT touch `.solve`; only `.apply` is in scope.
- 2-D Cartesian sweep currently uses anti-diagonal wavefront (legacy FD), not `dag_walk`. The unified matvec doesn't yet support 2-D Cartesian (`transport_operator_matvec_unified` raises `NotImplementedError` for `ny > 1`, `operator.py:1022-1027`). T.4 may need to stay 1-D-only until 2-D is wired through `dag_walk`.
- `StreamingTerms` is currently SN-specific despite the `geometry/` layer placement; MoC and CP share the primitive in principle but haven't migrated. T.4's rewire shouldn't BLOCK on cross-method generality — keep SN-only and let MoC/CP migrations land in their own waves.

These are deferred to "discover during execution" — the architectural commitment to one algebraic form across geometries is non-negotiable, but the concrete shape of T.4 will be refined when T.4 starts.

---

## 2. Dependencies

Wave T cannot start until Depth B completes. Specifically:

| Wave T substep | Depends on Depth B step | Why |
|---|---|---|
| T.1 (BC realizers) | D-B (TensorProductSpace), D-B+1 (specular pattern) | Extends D-B+1's specular pattern to vacuum / periodic / white / albedo. Needs the same primitive. |
| T.2 (fission) | D-B (TensorProductSpace) | Codomain typing of the rank-1 product needs `TensorProductSpace` |
| T.3 (scattering) | D-B, D-E (HarmonicMomentField) | Scattering's codomain is moment-space; HarmonicMomentField's typed form is needed as the dispatch target |
| T.4 (streaming) | D-B, D-H (AngularFlux) | Streaming's domain/codomain is angular-flux space; the typed AngularFlux is the dispatch target |

In practice, Wave T can start as soon as D-K (Depth B's last cleanup step) lands.

---

## 3. Architecture decisions (settled)

These decisions are made (Depth B plan §11.3 + this plan §1):

1. **`TensorProductOperator` and `SumOfTensorProductsOperator`** are the load-bearing primitives. No bespoke alternatives — every Wave-T-touched operator either IS one of these, or IS a leaf that composes via `&` into one.
2. **`TensorProductSpace`** is the codomain-typing primitive (Depth B D-B).
3. **The connection-coefficient framing** (geometric data in `StreamingTerms`, angular closure strategy in `PoleAngularClosure`) is preserved. Wave T does not refactor the geometric primitives — it consumes them at the operator level.
4. **Bit-identity discipline.** Each substep is a STRICT bit-identical rewire at the matvec level (same numpy reductions, just expressed through the `TensorProductOperator.apply` fold). The 10 pre-existing DD-regression failures stay at the same set; no new failures.
5. **`SNStreamingOperator` is OUT OF SCOPE.** It's the legacy `L + C` bundle, on the retirement queue from Phase G substep 3+4.c. Wave T targets the modern `StreamingOperator` leaf (`operator.py:1801`). `SNStreamingOperator` retires with the GMRES adapter rewrite in a separate cleanup wave.
6. **`InvertibleOperator.solve` is OUT OF SCOPE.** The `.solve` path uses `transport_sweep` (WDD algorithm) — a procedural inverse with no tensor-product decomposition. Wave T only touches `.apply` paths.

---

## 4. Out of scope (deferred)

Per `[[feedback_unify_after_two_instances]]` and `[[feedback_architecture_forward_not_legacy_fit]]`:

1. **MoC / CP consumption of `TensorProductOperator`.** Wave T is SN-only. MoC and CP migrations to the shared connection-coefficient primitive land in their own waves.
2. **2-D Cartesian streaming.** Currently uses anti-diagonal wavefront FD; not yet routed through `dag_walk`. T.4 may stay 1-D until 2-D is wired separately.
3. **`Representation` polymorphism** (§32.8 — dense / sparse / matrix-free / tensor-train). Wave T establishes the TYPED operator algebra; `Representation` is a separate abstraction over how each factor is stored / applied. Deferred.
4. **Adjoint propagation via `(A ⊗ B)^† = A^† ⊗ B^†`.** The identity lights up automatically once factors carry `.H`; whether to USE it for adjoint Krylov is a separate wave.
5. **`SumOfTensorProductsOperator._build` fusion / canonicalisation.** The `Σ_ℓ A_ℓ ⊗ B_ℓ` form may admit factor-merging (e.g., two summands with the same Ω factor collapse). Out of scope; the unfused form is correct and clearer.

---

## 5. Verification strategy

Per substep:

### 5.1 Bit-identity contract

Every Wave T substep is STRICT bit-identical at the matvec level. `TensorProductOperator.apply` is a fold `for op in self.ops: out = op.apply(out)` (`operator.py:1094-1098`) — when factors are `(P, I, I)`, the fold reduces to `P.apply(I.apply(I.apply(x))) == P.apply(x) == np.take(x, perm, axis=0)`. Identity factors are no-ops; the inner numerics are unchanged.

Existing tests for each operator should remain GREEN without modification at the value level. New assertions are added:
- `isinstance(result, TensorProductOperator)` — type signature change.
- `result.assert_separable()` — the L1 invariant fires.
- `result.codomain.shape == expected_shape` — postcondition from `TensorProductSpace` typing.

### 5.2 L1 MMS gates

The L1 MMS gates (`tests/sn/test_mms_aniso.py`, `tests/sn/l1_analytical/...`) are the ground truth. They must stay green at every Wave T commit. If a substep breaks an L1 gate, the substep is wrong (regardless of bit-identity claims — bit-identity is necessary but not sufficient when reductions reorder).

### 5.3 Algebraic-identity gates (NEW)

Wave T introduces NEW tests that exercise the identities the type system gains:

- **Adjoint distributivity**: `(A & B & C).H` equals `A.H & B.H & C.H` value-wise.
- **Composition distributivity**: `(A & B) ∘ (C & D)` equals `(A ∘ C) & (B ∘ D)` value-wise (where factor compositions are well-defined).
- **Identity-factor reduction**: `A & IdentityOperator()` equals `A` value-wise on the relevant axis.
- **assert_separable success**: every Wave-T-produced operator passes `assert_separable()`.

These are L0 foundation tests (`tests/numerics/test_tensor_product_identities.py` — new file).

### 5.4 Performance regression gates

`TensorProductOperator.apply` adds a Python-level fold over factors. For 3-factor operators, that's 3 function calls vs. 1 numpy operation — measurable overhead for inner-Krylov hot paths. Add a perf regression gate (`tests/sn/performance/test_wave_t_overhead.py` — new) measuring the `(L + C).apply` walltime delta. Acceptable: ≤ 5% slowdown on the 1-D slab Krylov benchmark. If overhead exceeds 5%, investigate fusion / `_build` flattening before shipping.

---

## 6. Sequencing

### Step T.1 — Generalize remaining BC realizers to tensor-product form

Extend D-B+1's specular pattern to the remaining BC kinds in `SNBoundaryRealizer.realize`:

| BC | Current realisation | Wave T.1 form |
|---|---|---|
| Vacuum | `IncomingOrdinateMaskTensor(axis=0)` | `M_inflow & I_group & I_face` |
| Periodic | `PeriodicWrapOperator()` | `I_angle & I_group & P_patch_pair` |
| White | `AngularAverageOperator(axis, sign)` | `K_avg_angle & I_group & I_face` |
| Albedo (isotropic, α < 1) | `ScaledOperator(α, I)` | `ScaledOperator(α, I_angle & I_group & I_face)` |

For each rewire:
- The realised operator's `.apply` is bit-identical at the value level.
- The codomain shape check (impossible to encode with the legacy single-axis form) becomes a postcondition.
- `assert_separable()` passes.

After T.1: FIVE production tensor-network instances (1 specular from D-B+1 + 4 here). The abstraction is empirically validated per `[[feedback_unify_after_two_instances]]`; the "right abstraction shrinks the count" property is checkable (6 BC operator types → 1 composition pattern over a small factor library).

### Step T.2 — Fission as rank-1 `TensorProductOperator`

`F = χ ⊗ νΣ_f` per §15 line 2003-2019. Today `FissionOperator.apply` is `np.einsum("gxy,gxy->xy", sig_p, phi) * chi[None, :, :]` (`fission.py:257-260`). The rewire:

```python
# T.2 form (rank-1, single TensorProductOperator)
F_apply = (DiagonalOperator(chi, axis=...) & ProductionRateOperator(nu_sig_f))
```

Open design question (T.2 starting brief): what's the cleanest factor split? `F = χ ⊗ (ν·integrate(Σ_f · φ))` is conceptually `OneAxisOuterProduct ⊗ FullContraction`. The two factors operate on different axes — `χ` produces a group-axis profile, `νΣ_f · φ` integrates over groups to produce a per-cell rate. The factorisation crosses the group axis in a non-trivial way. Detailed factor design lands at T.2 start.

Verification: existing fission tests stay green; new test pins `isinstance(fission_op.apply, ...) == TensorProductOperator`-shaped invariant; `k_∞ = νΣ_f / Σ_a` homogeneous-reflective check is the structural-independence ground.

### Step T.3 — Scattering as `SumOfTensorProductsOperator` per §15.2

The grand-report-named identity (line 2058-2061):

```python
S = SumOfTensorProductsOperator([
    SigmaMoment(xs.scatter, ell) & AngularMomentOperator(ell) & GroupScatteringMatrix(xs.scatter, ell)
    for ell in range(xs.scatter.order + 1)
])
```

The largest Wave T payoff (Krylov hot path). Touches `LegendreMomentScattering` (`scattering.py:148-268`) and `ScatteringOperator` (`scattering.py:271-993`). The R · Λ · M pipeline (`scattering.py:625-657`) — currently three classes composed by direct method calls — becomes `(R & I_angle & I_group) ∘ Λ ∘ (M & I_angle & I_group)` operator-algebra composition.

The per-material per-ℓ einsum at `material_xs_field.py:515-572` is the inner numerics. Bit-identical preserved.

Verification: 2-group heterogeneous L1 MMS gate is the structural ground; existing scattering tests pin value-level identity; new test pins the `SumOfTensorProductsOperator` structure.

### Step T.4 — `StreamingOperator.apply` as `L_spatial + L_angular_redist`

The universal decomposition per the connection-coefficient framing:

```python
L_apply = L_spatial.apply(psi) + L_angular_redist.apply(psi)

L_spatial = sum over spatial axes of  D_axis(streaming_terms) & Ω_axis(quad) & I_group
L_angular_redist = R_polar(connection_coefficients) & I_x & I_group   # ZeroOperator for slab
```

Both summands are `TensorProductOperator`s (`L_spatial` is a `SumOfTensorProductsOperator` over spatial axes). `L_angular_redist` is the M-M closure for sphere/cylinder; degenerate (zero operator) for slab.

The geometric data lives in `StreamingTerms` (`geometry/reduced_operator.py:151`) for spatial face-area ratios, and in `PoleAngularClosure` (`sn/spatial/pole_angular_closure.py:191`) for the angular connection coefficients. T.4 does NOT refactor those primitives — it consumes them at the operator level via `D_axis(streaming_terms)` / `R_polar(connection_coefficients)` factory methods.

**Complications expected** (§1.3):
- Whether `R_polar` factors cleanly as a `TensorProductOperator` factor or needs a bespoke shape (the M-M half-grid sequential sweep coupling is not diagonal). Discover during execution.
- 2-D Cartesian may stay on the legacy FD path until `dag_walk` is wired through.
- Curvilinear may surface other coupling that resists the decomposition; if so, document the principled non-factorization and adapt the architecture rather than abandoning the algebra.

Verification: slab L1 gates (slab is the simplest case, must work first), then sphere and cylinder. The structural-independence reference is the existing `transport_operator_matvec_unified` output (the unified procedural baseline that Wave T rewires).

### Step T.5 — Documentation + retirement of dead infrastructure (close-out)

After T.1-T.4 land:
- Update `docs/theory/operator_algebra.rst` to document the §15-§16A tensor-network instances now active in production.
- Update `docs/api/numerics.rst` for `TensorProductSpace` / `TensorProductOperator` consumer table.
- Retire any operator-algebra leaves that became unused (e.g., if `PermutationOperator(axis=0)` is no longer instantiated outside of `TensorProductOperator` factors, it could fold into the factor itself).
- Run the `assert_separable` gate across the whole solver suite — every Wave-T-touched operator's `apply` should produce a `TensorProductOperator`-typed result.

---

## 7. Risk register

| Rank | Risk | Mitigation |
|---|---|---|
| 1 | T.4 curvilinear factorisation doesn't admit a clean tensor product — M-M half-grid is sequentially coupled. | Discover during execution. Acceptable fallback: factor what's factorable, document non-factorisable pieces, ensure the algebra still composes via `OperatorSum`. The architectural commitment is ONE algebraic FORM, not that every factor is a clean `TensorProductOperator`. |
| 2 | Python-level fold overhead in `TensorProductOperator.apply` (3 function calls per factor) hurts Krylov hot path > 5%. | Perf regression gate at §5.4. If exceeded, investigate `_build` fusion (collapse `(A & I & I) ∘ (B & I & I) → (A ∘ B) & I & I`) or factor caching. |
| 3 | Bit-identity breaks at FP-non-associativity ULP because `TensorProductOperator.apply` reorders reductions vs. the legacy fused einsum. | Per `vv-principles` §"Bit-identity vs principled-equivalence": acceptable if (a) the new path is principled (each factor is a named operator with definite physical meaning), (b) verified against an independent reference (k_∞ closed-form, L1 MMS), (c) drift is dimensionally explainable (`reduction_depth × ULP`). Document the contract relaxation at the affected regression snapshot. |
| 4 | T.3 scattering rewire touches the R · Λ · M pipeline used by EVERY iteration of EVERY SN solve. Bugs are expensive. | Test-architect proactive dispatch BEFORE T.3 implementation (per `subagent-handoff-protocol`'s proactive trigger for operator-algebra carves crossing subsystem boundaries). |
| 5 | Wave T's substeps depend on Depth B's typed Fields (T.3 needs HarmonicMomentField from D-E, T.4 needs AngularFlux from D-H). If Depth B's later steps slip, Wave T slips with them. | Sequential dependency; no shortcut. Communicate slippage up the parent-plan chain. |
| 6 | `SNStreamingOperator` retirement (separate cleanup wave) interacts with Wave T's `StreamingOperator` rewire. | Wave T treats `SNStreamingOperator` as out-of-scope (§3). The separate cleanup wave can land before, during, or after Wave T independently. |

---

## 8. Plan checkpoint — to be filled at approval time

- [ ] User approval of Wave T as a wave (vs. embedded in Depth B or P3.4).
- [ ] User approval of sequencing T.1 → T.2 → T.3 → T.4 → T.5.
- [ ] User approval of out-of-scope deferrals (§4).
- [ ] User confirmation: each substep is its own commit, not bundled.
- [ ] User confirmation: test-architect dispatch BEFORE T.3 (scattering rewire).

Once these checkpoints are confirmed AND Depth B is complete, this plan becomes the source of truth for the Wave T implementation session.

---

## 9. Exit Route — where this wave leads when complete

When Wave T's T.1 through T.5 commit, the following invariants hold:

1. **All four operator classes** (`SNBoundaryRealizer`-produced BCs, `FissionOperator`, `ScatteringOperator`, `StreamingOperator`) produce `TensorProductOperator` / `SumOfTensorProductsOperator` instances from `.apply` and `.realize`. The `assert_separable` invariant fires.
2. **L1 MMS gates** stay green. The 10 pre-existing DD-regression failures stay at the same failure set. The new algebraic-identity gates (§5.3) pass.
3. **Bit-identity** holds at the matvec value level for the touched operators, or — for the rare reduction-reordering cases — the principled-equivalence three-criteria test passes (`vv-principles` §"Bit-identity vs principled-equivalence").
4. **Performance** is within ≤ 5% of pre-Wave-T baseline on the 1-D slab Krylov benchmark.
5. **`TensorProductOperator` / `SumOfTensorProductsOperator`** have FIVE+ production consumers (4 BC kinds + fission + scattering + streaming). The Wave-0 primitives are no longer shipped-but-unused.

**Hand-off to parent-plan step P3.4** (Problem/Solver split). P3.4 picks up:
- Typed Fields (from Depth B) AS operator domain/codomain.
- Tensor-network operators (from Wave T) AS the algebraic substrate for Problem ABCs.
- Together they enable `CriticalityProblem(loss=L+C-S, fission=F)` where each component carries explicit tensor structure that `PowerIteration` and `SourceIteration` can introspect at construction time.

The Depth B + Wave T pair is the architectural foundation; P3.4 is its first major consumer.

---

## Pointers

- **Parent plan**: `.claude/plans/moment_space_and_layering_plan.md` — Phase 3 sequencing (Wave T inserted between Depth B and P3.4 per §11.3 of the Depth B plan).
- **Depth B plan**: `.claude/plans/depth_b_field_on_function_space.md` — Wave T's predecessor; §11.3 documents the architectural decision to insert Wave T.
- **Grand report**: `.claude/plans/neutron_transport_grand_report_v3.md`. Wave T's authority — §15 (tensor products), §16A.10 (BC as tensor network), §35 (commandments), north-star line 5697.
- **Connection-coefficient primitive**: `orpheus/geometry/reduced_operator.py:1-30` (the SO(3)-charts framing for spherical / cylindrical streaming).
- **Cell-balance algebra (geometry-blind by data)**: `orpheus/sn/spatial/cell_balance.py:120-198`.
- **L1 primitives (shipped, unused until Wave T)**:
  - `orpheus/numerics/operator.py:1004-1122` — `TensorProductOperator`.
  - `orpheus/numerics/operator.py:1125-1215` — `SumOfTensorProductsOperator`.
  - `orpheus/numerics/operator.py:1203-1215` — `assert_separable()`.
- **Cardinal rules from memory**: `feedback_unify_after_two_instances`, `feedback_no_method_implementer_for_surgical_carves`, `feedback_architecture_forward_not_legacy_fit`, `feedback_elegance_causes_collapse`, `default-test-mode-is-optimize`.
