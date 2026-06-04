# Wave T — Tensor-Network Operator Algebra

**Branch:** `refactor/moment-space-and-layering`
**Worktree:** `.claude/worktrees/moment-space-and-layering/`
**Phase status:** PENDING. Lands AFTER Depth B (`depth_b_field_on_function_space.md`) completes (D-K shim retirement, completed 2026-05-30). Precedes **Wave O** (`.claude/plans/wave_o_operator_typing.md` — operator-role typing per Issue #208), which in turn precedes parent-plan step P3.4 (Problem/Solver split). Updated sequencing (2026-05-30): `Depth B ✓ → Wave T → Wave O → P3.4 → P3.6`.

**Date:** 2026-05-27. Plan author: main agent + explorer audit (`.claude/agent-memory/explorer/` — `function_space_typed_field_audit.md` for Depth B background; the conversation that motivated Wave T is recorded in the Depth B plan §11.3).

**Status:** STUB. The architectural commitment is settled (see §1); detailed designs for each substep land when each substep starts (depth-on-demand per `[[feedback_no_method_implementer_for_surgical_carves]]` — main-agent direct authorship with turn-by-turn user steering).

**2026-05-30 polish update (post-Depth-B-close):**
Depth B closed completely on 2026-05-29 (commits `c97897d`, `8a8ddbf`,
`da9b73b`, `4a53737`, `dadf4e8`, `a5a7ff9`).  This pass refreshes file
line numbers, splits TP-operator vs TP-space provenance, removes the
stale `SNStreamingOperator`-retirement paragraph (`dadf4e8` retired the
class), accurate-counts production tensor-network consumers (ONE today —
the specular BC at `boundary_realizer.py:164-166` — not zero as
originally written), and clarifies the 2-D Cartesian status (the
unblock landed but the path is procedural FD via `_apply_2d_cartesian`,
not yet tensor-product-shaped; the `NotImplementedError` gate in
`transport_operator_matvec_unified` is unreachable-by-routing).  Six
BC kinds exist in `SNBoundaryRealizer.realize` (not four as §6 T.1
tabulates) — see §6 T.1 below for the updated table.

---

## 0. Pickup checklist (read first)

If you are picking this plan up in a fresh session:

1. **Read this plan top-to-bottom.** No section is optional.
2. **Verify Depth B completed**: `git log --oneline -25` should show the D-A through D-K commits AND the `1e0bb98` Field ABC plus `53bc986` ScalarFlux migration as landed. The parent plan's `moment_space_and_layering_plan.md` step P3.3 should be marked LANDED (re-scoped through Depth B).
3. **Read these grand-report sections**:
   - §15 lines 2003-2086 — "The tensor section": V = X ⊗ Ω ⊗ G; native operator form = sum of tensor products.
   - §15.1 lines 2031-2045 — streaming as `L = Σ_axis (D_axis ⊗ Ω_axis ⊗ I_g)`.
   - §15.2 lines 2046-2086 — scattering as `S = Σ_ℓ (Σ_{s,ℓ} ⊗ A_ℓ ⊗ G_ℓ)`.
   - §16A.10 lines 3142-3197 — BC as tensor network; canonical separable form `B = G_patch ⊗ K_omega ⊗ K_g`.  Note: `BoundaryGeometryTensor` and `BoundaryResponseTensor` are grand-report aspirational symbols — they do NOT exist as code today.  The closest ORPHEUS analogues are `InflowTraceSpace` / `OutflowTraceSpace` for the `G_patch` axis structure and `IncomingOrdinateMaskTensor` / the specular permutation for `K_omega`.  Wave T does NOT introduce these named tensors; it lifts existing BC operators into `TensorProductOperator` factors that fill the same algebraic roles.
   - §35 line 5675 — 20 commandments: "Represent multigroup scattering as a sum of tensor products."
   - North-star line 5697: "tensor products are foundational".
4. **Read the explorer's two audit memos** that motivated Wave T (referenced in Depth B plan §11.3 footnotes; if absent, regenerate via explorer dispatches against the grand report).
5. **Read the existing L1 primitives**:
   - `orpheus/numerics/operator.py:1058-1177` — `TensorProductOperator` (definition, `__and__` dunder via `_build`, `apply` fold, `solve` fold for fully-invertible factors).
   - `orpheus/numerics/operator.py:1179-1269` — `SumOfTensorProductsOperator`.
   - `orpheus/numerics/operator.py:1257-1269` — `assert_separable()` method (shipped, currently called only in tests; lights up post-T.1).
   - `orpheus/numerics/space.py:301` — `TensorProductSpace` (landed in Depth B D-B, commit `c2f968a`).
   - `orpheus/sn/boundary_realizer.py:164-166` — the **sole production `TensorProductOperator` instance today** (D-B+1, specular BC).
   - `orpheus/sn/angular_operator.py:37` — `AngularAverageOperator` (legacy single-axis; T.1's white-BC lift target).
   - `orpheus/sn/operator.py:645` — modern `StreamingOperator` leaf (T.4 target).  `CollisionOperator` at `:973`, `InvertibleOperator.solve` public entry at `:1280` (body delegates to `_solve_timed_full_field` at `:1368`; `transport_sweep` invoked at `:1462`).
   - `orpheus/sn/operator.py:203` — `transport_operator_matvec_unified` (NOT `numerics/`; 2-D `NotImplementedError` at L328-333 is unreachable-by-routing today — see §1.3).
   - `orpheus/sn/operator.py:830` — `_apply_2d_cartesian` (procedural FD path that intercepts 2-D Cartesian via the dispatcher at L809-812; not currently tensor-product-shaped).
6. **Pick up at the leftmost incomplete step in §6** below.

The `feedback_no_method_implementer_for_surgical_carves` rule applies: this is the main agent's work with turn-by-turn user steering. Do NOT batch via method-implementer.

---

## 1. The principle

**The deterministic multigroup neutron transport problem already lives on a tensor-product space.** Per grand report §15 line 2003-2019:

> *"The deterministic multigroup problem already lives on a tensor-product space: V = X ⊗ Ω ⊗ G. The stochastic problem extends this: V = X ⊗ Ω ⊗ G ⊗ Ξ_1 ⊗ ... ⊗ Ξ_d. The native operator form is often a sum of tensor products: A = Σ_k A_x^(k) ⊗ A_ω^(k) ⊗ A_g^(k) ⊗ A_ξ^(k)."*

The codebase today implements this tensor-product structure via flat-axis numpy with implicit broadcasting — `np.einsum("gxy,gxy->xy", chi, phi)` for fission, `np.take(x, perm, axis=0)` for specular BC, etc. Numpy broadcasting IS the implementation of `... ⊗ I_other_axes` — but the *type system* doesn't know, the *operator algebra* doesn't know, and the algebraic identities `(A ⊗ B)^† = A^† ⊗ B^†` and `(A ⊗ B) ∘ (C ⊗ D) = (A ∘ C) ⊗ (B ∘ D)` (documented at `operator.py:1024-1033`) cannot fire because no production operator is ever instantiated as `TensorProductOperator`.

**Wave T closes the gap.** The L1 tensor-product **operators** (`TensorProductOperator`, `SumOfTensorProductsOperator`, `assert_separable`) shipped in Wave 0 (commit `bc1253e`, 2026-05-10); the L1 tensor-product **space** (`TensorProductSpace`) shipped in Depth B D-B (commit `c2f968a`, 2026-05-27).  Today exactly ONE production consumer exists: the specular BC at `boundary_realizer.py:164-166`, lifted in D-B+1.  `SumOfTensorProductsOperator` has zero production consumers — Wave T introduces its first via T.3 (scattering) and T.4 (streaming).  Wave T extends the D-B+1 pattern to the remaining BC kinds and rewires fission, scattering, and the modern streaming operator to use the algebra natively.  Wave T is the consumer side of the Wave-0 + Depth B D-B infrastructure.

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
- 2-D Cartesian is reachable but procedural.  Depth B D-H.2-C4 (commits `8721f3f`, `1ae8d5f`, `7746fd3`) added an L2-native 2-D Cartesian apply path: `StreamingOperator.apply` intercepts 2-D at `sn/operator.py:809-812` and routes to `_apply_2d_cartesian` at `sn/operator.py:830`.  The `NotImplementedError` in `transport_operator_matvec_unified` (`sn/operator.py:328-333`) is now unreachable-by-routing.  However, `_apply_2d_cartesian`'s body is **procedural FD** (anti-diagonal wavefront over `mu_x/mu_y` masks), NOT tensor-product-shaped.  T.4 therefore has a choice: (a) tensor-lift the 2-D path in scope (extension of the `L_spatial = D_x ⊗ Ω_x ⊗ I_g + D_y ⊗ Ω_y ⊗ I_g` decomposition that §15.1 names explicitly for 2-D), or (b) leave the 2-D path as procedural FD and target 1-D only.  The decision is "discover during execution" — likely a separate companion wave if the procedural path resists clean factorisation.
- `StreamingTerms` is currently SN-specific despite the `geometry/` layer placement; MoC and CP share the primitive in principle but haven't migrated. T.4's rewire shouldn't BLOCK on cross-method generality — keep SN-only and let MoC/CP migrations land in their own waves.

These are deferred to "discover during execution" — the architectural commitment to one algebraic form across geometries is non-negotiable, but the concrete shape of T.4 will be refined when T.4 starts.

---

## 2. Dependencies

Depth B **closed completely** on 2026-05-29 (commits up to `a5a7ff9`).  Every dependency Wave T originally listed (`TensorProductSpace` from D-B, `HarmonicMomentField` from D-E, `AngularFlux` from D-H, the specular tensor-network pattern from D-B+1) is LANDED in production.  Wave T is unblocked.

For historical record, the dependency map was:

| Wave T substep | Depended on | Status |
|---|---|---|
| T.1 (BC realizers)    | D-B `TensorProductSpace`, D-B+1 specular pattern | ✓ landed |
| T.2 (fission)         | D-B `TensorProductSpace`                          | ✓ landed |
| T.3 (scattering)      | D-B, D-E `HarmonicMomentField`                    | ✓ landed |
| T.4 (streaming)       | D-B, D-H `AngularFlux` / `TimedFullField`         | ✓ landed |

---

## 3. Architecture decisions (settled)

These decisions are made (Depth B plan §11.3 + this plan §1):

1. **`TensorProductOperator` and `SumOfTensorProductsOperator`** are the load-bearing primitives. No bespoke alternatives — every Wave-T-touched operator either IS one of these, or IS a leaf that composes via `&` into one.
2. **`TensorProductSpace`** is the codomain-typing primitive (Depth B D-B).
3. **The connection-coefficient framing** (geometric data in `StreamingTerms`, angular closure strategy in `PoleAngularClosure`) is preserved. Wave T does not refactor the geometric primitives — it consumes them at the operator level.
4. **Bit-identity discipline.** Each substep is a STRICT bit-identical rewire at the matvec level (same numpy reductions, just expressed through the `TensorProductOperator.apply` fold). The 10 pre-existing DD-regression failures stay at the same set; no new failures.
5. **`InvertibleOperator.solve` is OUT OF SCOPE.** The public `.solve` entry at `orpheus/sn/operator.py:1280` delegates to `_solve_timed_full_field` at `:1368`, which invokes `transport_sweep` (`orpheus/sn/sweep.py:99`) at `orpheus/sn/operator.py:1462` — the WDD algorithm, a procedural inverse with no tensor-product decomposition. Wave T only touches `.apply` paths.

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

Extend D-B+1's specular pattern to the remaining BC kinds in `SNBoundaryRealizer.realize`.

`SNBoundaryRealizer.realize` (`orpheus/sn/boundary_realizer.py:123-199`) dispatches to six BC kinds today.  T.1's scope is to lift the FIVE that are not yet tensor-network-shaped; the specular row is the D-B+1 instance, included here for completeness.

| BC                          | Source                          | Current realisation                                                                              | T.1 lift target                                                                                                                                          |
|-----------------------------|---------------------------------|--------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------|
| Vacuum                      | `boundary_realizer.py:131-148`  | `IncomingOrdinateMaskTensor(axis=0)` — single-axis                                               | `M_inflow & I_group & I_face`                                                                                                                            |
| Specular ✓ (D-B+1)          | `boundary_realizer.py:150-173`  | `PermutationOperator(perm, axis=0) & IdentityOperator()` — already 2-factor TP                   | **DONE**; T.1 leaves it untouched                                                                                                                        |
| White                       | `boundary_realizer.py:175-181`  | `AngularAverageOperator.from_quadrature(...)` (`sn/angular_operator.py:37`) — single-axis        | `K_avg_angle & I_group & I_face`                                                                                                                         |
| Albedo (isotropic, α < 1)   | `boundary_realizer.py:183-188`  | `ScaledOperator(α, IdentityOperator())`                                                          | `ScaledOperator(α, I_angle & I_group & I_face)`                                                                                                          |
| Periodic                    | `boundary_realizer.py:190-191`  | bare `PeriodicWrapOperator()` (`numerics/operator.py:1010`) — single-axis                        | `I_angle & I_group & P_patch_pair`                                                                                                                       |
| PrescribedInflow            | `boundary_realizer.py:193-199`  | `IncomingSourceOperator(...)` — affine, not linear                                               | Open: affine-source BCs may need `AffineFullOperator` typing in Wave O before T.1 lifts them; T.1 may defer this row to Wave O                           |

For each remaining-rewire row:
- The realised operator's `.apply` is bit-identical at the value level.
- The codomain shape check (impossible to encode with the legacy single-axis form) becomes a postcondition.
- `assert_separable()` passes.

After T.1: FIVE production tensor-network instances (1 specular from D-B+1 + 4 here, deferring PrescribedInflow to Wave O if it doesn't cleanly fit the linear-tensor-product form). The abstraction is empirically validated per `[[feedback_unify_after_two_instances]]`; the "right abstraction shrinks the count" property is checkable (6 BC operator types → 1 composition pattern over a small factor library).

### Step T.2 — Fission as rank-1 `TensorProductOperator`

`F = χ ⊗ νΣ_f` per §15 line 2003-2019. Today `FissionOperator.apply` is `np.einsum("gxy,gxy->xy", sig_p, phi) * chi[None, :, :]` (`fission.py:257-260`). The rewire:

```python
# T.2 form (rank-1, single TensorProductOperator)
F_apply = (DiagonalOperator(chi, axis=...) & ProductionRateOperator(nu_sig_f))
```

Open design question (T.2 starting brief): what's the cleanest factor split? `F = χ ⊗ (ν·integrate(Σ_f · φ))` is conceptually `OneAxisOuterProduct ⊗ FullContraction`. The two factors operate on different axes — `χ` produces a group-axis profile, `νΣ_f · φ` integrates over groups to produce a per-cell rate. The factorisation crosses the group axis in a non-trivial way. Detailed factor design lands at T.2 start.

Verification: existing fission tests stay green; new test pins `isinstance(fission_op.apply, ...) == TensorProductOperator`-shaped invariant; `k_∞ = νΣ_f / Σ_a` homogeneous-reflective check is the structural-independence ground.

### Step T.3 — Scattering as `OperatorSum` of per-ℓ summands (NOT SOTP per Q6 honest-fallback)

**T.3 IMPLEMENTATION (2026-05-30, commits `9f85c5d` / `03bcdba` / `2595d2f`) — DEVIATION from §15.2 aspiration**

The grand-report §15.2 form names the identity:

```python
S = SumOfTensorProductsOperator([
    SigmaMoment(xs.scatter, ell) & AngularMomentOperator(ell) & GroupScatteringMatrix(xs.scatter, ell)
    for ell in range(xs.scatter.order + 1)
])
```

T.3 design-fork Q6 surfaced that this §15.2 SOTP form **fails the disjoint-axes contract**: the per-material per-ℓ einsum at `material_xs_field.py:515-572` couples the group axis (per-material XS lookup `mat_xs.sig_s_legendre(mid)[ell]`) with the spatial axis (per-cell material id from `mat_xs.cells_by_material`). The user resolved Q6 to the **math-honest fallback** at AskUserQuestion time: ship `OperatorSum` over per-ℓ bespoke summands.

**Shipped form (T.3b/c/d/e, ✓ landed)**:

```python
class _PerLegendreOrderScattering(LinearOperator):
    """Per-ℓ summand: applies R_ℓ ∘ Λ_ℓ ∘ M_ℓ for ONE Legendre order ℓ."""
    ell: int                 # specific ℓ ∈ {1, ..., scattering_order}
    L: int                   # total order
    mat_xs: MaterialXSField  # for per-material per-ℓ einsum
    weights: np.ndarray; Y: np.ndarray

class ScatteringOperator(LinearOperator):
    @cached_property
    def kernel_summands(self) -> tuple[_PerLegendreOrderScattering, ...]:
        """ONE summand per Legendre order ℓ ∈ {1, ..., scattering_order}.
        Note: P0 (ℓ=0) stays in the iso-fast-path, NOT in the kernel —
        Q3 Option (β) resolution.  ``len(kernel_summands) == scattering_order``,
        not order + 1.
        """
        ...

    @cached_property
    def kernel(self) -> LinearOperator:
        """`OperatorSum` of `kernel_summands` via `functools.reduce(add)`.
        For `scattering_order == 0` (P0 only), returns `ZeroOperator`.
        For `scattering_order == 1` (P1), `reduce` returns the singleton
        summand directly (NOT an OperatorSum)."""
        ...
```

The R · Λ · M pipeline (`scattering.py:625-657`) — three classes composed by direct method calls — becomes per-ℓ `OperatorSum` summands. The `build_aniso_source` body collapsed from 22 → 14 LOC (single-source-of-truth via `self.kernel.apply`).

The per-material per-ℓ einsum at `material_xs_field.py:515-572` is the inner numerics, preserved bit-exact inside `_PerLegendreOrderScattering.apply`.

Verification (landed): 2-group heterogeneous L1 MMS gate is the structural ground; existing scattering tests pin value-level identity; new tests pin the `OperatorSum`-of-per-ℓ structure + per-arm bit-identity vs pre-T.3 snapshots.

**Master condition surfaced (T.3 resolution → T.4 amendment)**: SOTP requires Cartesian-product per-axis decomposition; coupled physics (per-material XS in T.3) falls back to `OperatorSum`. See §6 T.4 for the analogous T.4 application of this master condition (per-direction WDD + sequential M-M half-grid).

### Step T.4 — `StreamingOperator.apply` as `M_spatial + M_angular_redist - σ_t·ψ`

**T.4 IMPLEMENTATION (2026-05-31, commits `cb18fdb` / `c55b505` / `90e7d4e`) — multiple DEVIATIONS from the original §15.1 SOTP aspiration**

The user-endorsed honest decomposition (resolved via the T.4-spec Q1-Q5 + post-spec follow-up at AskUserQuestion time, before T.4b code landed):

```python
# Shipped T.4 form — properties named M_* per Q3 = (γ); see "naming" below.
class StreamingOperator:
    @cached_property
    def M_spatial(self) -> _MSpatialOperatorSum:
        """OperatorSum of per-direction-sign summands (M_x_forward + M_x_backward).
        For 1-D slab/sphere/cyl: exactly 2 per-direction `_SpatialSweepDirection`
        leaves.  NOT a `SumOfTensorProductsOperator` — the WDD recurrence is
        sequentially coupled along x.  Honest decomp per MA-Q1 master
        condition."""
        ...

    @cached_property
    def M_angular_redist(self) -> LinearOperator:
        """Bespoke `AngularRedistributionOperator` leaf for sphere/cyl;
        `ZeroOperator` for slab/Cartesian.  NOT a `TensorProductOperator`
        wrap — M-M half-grid recurrence is sequentially coupled along the
        angular axis (Hébert §3.9.4)."""
        ...

    def apply(self, psi):
        # L = M - C; M = M_spatial + M_angular_redist
        # Slab: routes through M_spatial.apply (M_angular_redist = Zero skip)
        # Curvilinear: STAYS on `transport_operator_matvec_unified` shortcut
        # for perf (the M_* decomposition exists for type-level inspection
        # — Wave O / adjoint / DSA — but does NOT drive the production
        # hot path for curvilinear).
        ...
```

**Five DEVIATIONS from the original SOTP form** (each resolved via AskUserQuestion ahead of code):

* **Q1 (2-D Cartesian)** = HYBRID. T.4 lifts 1-D only; `_apply_2d_cartesian`
  stays procedural FD. The cell-centre-proxy ↔ face-view-as-trace rewire
  is its own architectural payload (10% k_∞ drift comment at
  `sn/operator.py:862-868`); bundling violates
  `feedback_unify_after_two_instances`. T.4d adds a defensive
  source-hash pin (A2D-1).
* **Q2 (curvilinear angular)** = BESPOKE LEAF. `AngularRedistributionOperator`
  is a `LinearOperator` leaf, NOT a 3-factor TP wrap.  The M-M half-grid
  recurrence (Hébert §3.9.4 Eqs. 3.432-3.435) is sequentially coupled
  along the angular axis — `(leaf & I_x & I_g)` would false-assert
  separability the recurrence doesn't support.
* **Q3 (subtraction location)** = APPLY-BOUNDARY + RENAME. `σ_t·ψ`
  subtraction stays at `StreamingOperator.apply`'s boundary (Pattern 7
  producer-side normalisation). Properties renamed `M_spatial` /
  `M_angular_redist` (NOT `L_spatial` / `L_angular_redist`) because the
  discrete cell-balance algebra produces un-subtracted M-shaped
  contributions; the continuous L and the discrete M = L + C are not
  the same thing.
* **Q4 (slab degeneracy)** = STRUCTURAL PRESERVATION. Slab's
  `M_angular_redist = ZeroOperator()`; both properties always exist
  (Pattern 4 illegal-states-unrepresentable).
* **Q5 (`.solve` path)** = OUT OF SCOPE. The
  `InvertibleOperator.solve` body (WDD sweep) is the procedural
  inverse — NOT tensor-product-factorable. T.4 only touches `.apply`.

**Master condition (surfaced post-T.4b)**: M_spatial is **NOT** a
clean SOTP `(D_x & Ω_x & I_g)` for the same MA-Q1 master condition the
T.3 spec applied to curvilinear M-M:

> SOTP requires Cartesian-product per-axis decomposition; coupled
> physics (per-material XS in T.3, per-ordinate recurrence in
> T.4-curvilinear, **AND per-direction WDD in T.4-spatial**) falls
> back to `OperatorSum`.

The honest shape: M_spatial = `OperatorSum[M_x_forward, M_x_backward]`,
each summand a bespoke `_SpatialSweepDirection` leaf.

**Per-direction split rationale (user-endorsed post-Q2)**: separate
forward and backward sweep summands set up structural pollination for
(a) Wave O typing — per-direction BC dependency footprint is naturally
exposed; (b) adjoint propagation `(M_x_fwd).H = M_x_back`-with-swapped-
BCs; (c) DSA-class preconditioners (split per-direction for P_1
closure); (d) cross-method pollination with billiard/trajectory_resolvent
single-direction sweeps; (e) per-direction debugging slicing
(ERR-006 family diagnosis).

**Orchestrated apply (Design B)**: `_MSpatialOperatorSum` **overrides**
the default `OperatorSum.apply` (which would be `forward.apply +
backward.apply` = 1.5× perf because each standalone leaf must redo
the forward sweep to compute the backward sweep's `bc_outer` seed) and
runs the bidirectional matvec ONCE via
`transport_operator_matvec_unified`. The local-variable outer-face WDD
outflow that was the previous hidden coupling point is now the named
shared state of this orchestrator. Standalone
`_SpatialSweepDirection.apply` exists as a slow fallback for testing
/ Wave-O / adjoint inspection.

**Curvilinear M_spatial = `(L+C) − M_angular_redist` via subtraction**.
This introduces a minor architectural smell (M_spatial depends on
M_angular_redist for curvilinear; algebraically clean but
implementationally coupled). The smell is bounded by the
algebra-decomposition invariant test (`(L+C) == M_spat + M_ang` at
principled-equivalence ULP ~16·ULP per `vv-principles`
three-criteria gate).

The geometric data still lives in `StreamingTerms`
(`geometry/reduced_operator.py:151`) and `PoleAngularClosure`
(`sn/spatial/pole_angular_closure.py:191`). T.4 does NOT refactor
those primitives — it consumes them at the operator level
(Pattern 2 single source of truth).

Verification (landed): slab L1 gates (T.4b green-gate), then sphere
and cylinder (T.4c bit-identity vs pre-T.4 snapshots). The
structural-independence reference is the existing
`transport_operator_matvec_unified` output — STILL alive post-T.4c
because both `_MSpatialOperatorSum.apply` and the
`StreamingOperator.apply` curvilinear shortcut consume it (NOT
orphan code).

### Step T.5 — Close-out (documentation + plan/spec amendments + defensive pin)

After T.1-T.4 land:

**T.5.1 (this commit)** — Plan + verification-spec wording amendments. The original §6 T.3 / §6 T.4 / §9 wording assumed SOTP for scattering + streaming; the user-resolved Q6 (T.3) and Q1-Q5 (T.4) + post-spec follow-up shipped `OperatorSum`-of-bespoke-leaves instead. T.5.1 codifies the deviations so future agents read the CURRENT-STATE truth, not the original aspiration. Amends:

- §6 T.3 paragraph 1 — `SumOfTensorProductsOperator` → `OperatorSum`; `order + 1` summands → `order` (P0 in fast path)
- §6 T.4 — full rewrite of the form section to reflect Q1-Q5 + `M_spatial`/`M_angular_redist` naming + per-direction `OperatorSum` (NOT 3-factor TP) + bespoke `AngularRedistributionOperator` leaf + 2-D Cartesian hybrid
- §9 exit-route — clarify that SOTP has ZERO production consumers; the algebraic SHAPE varies by physics coupling (TP, OperatorSum-of-leaves, or bespoke leaf); Wave O typing must accommodate non-SOTP summands

**T.5.2** — Sphinx documentation updates:

- `docs/theory/operator_algebra.rst` to document the §15-§16A tensor-network instances now active in production AND the MA-Q1 master condition (SOTP applicability boundary). Includes the per-direction `OperatorSum` decomposition for streaming and the per-ℓ `OperatorSum` decomposition for scattering as the **canonical examples** of "coupled-physics OperatorSum fallback".
- `docs/api/numerics.rst` consumer table for `TensorProductSpace` / `TensorProductOperator` / `OperatorSum` (the latter now load-bearing for T.3+T.4 production consumers).

**T.5.3 (T.4d)** — Defensive A2D-1 source-hash pin on `_apply_2d_cartesian` (per Q1 hybrid resolution): guards against silent author drift on the 2-D Cartesian path that T.4 left procedural. Small focused test.

**OUT OF T.5 SCOPE — moved to future micro-waves**:

- ~~Retire `transport_operator_matvec_unified`~~ — explorer's pre-T.4 audit (`adf8bfb5f7aa1556c`) recommended retirement at T.5; post-T.4c re-audit shows the function is STILL load-bearing (called from `_MSpatialOperatorSum.apply` for slab + curvilinear via subtraction, AND from `StreamingOperator.apply` curvilinear shortcut). NOT orphan. The retirement is contingent on either (a) re-architecting `_MSpatialOperatorSum.apply` to walk cells directly (~200 LOC), or (b) curvilinear `StreamingOperator.apply` switching from the unified shortcut to `M_spat + M_ang` composition (perf hit). Deferred to a separate post-T.5 cleanup micro-wave.
- ~~Soft leverage opportunity (cache unification)~~ — `_ensure_coll_cache` could read `a_attenuation` from `M_spatial.materialize_inverse_cache(sigma_t)`. ~50 LOC, 1 sweep_cache.py test. Flagged in the explorer report as separate post-T.4 micro-wave. Out of T.5 scope.
- ~~`assert_separable` suite-wide gate~~ — the invariant is structurally meaningful ONLY for the TP/SOTP-shaped operators (BCs + fission). It is inapplicable to scattering's `OperatorSum`-of-per-ℓ-leaves, streaming's per-direction `OperatorSum`, and the `AngularRedistributionOperator` bespoke leaf. The narrow per-operator `assert_separable` checks already exist (T.1/T.2 tests); a suite-wide check would be misleading.

---

## 7. Risk register

| Rank | Risk | Mitigation |
|---|---|---|
| 1 | T.4 curvilinear factorisation doesn't admit a clean tensor product — M-M half-grid is sequentially coupled. | Discover during execution. Acceptable fallback: factor what's factorable, document non-factorisable pieces, ensure the algebra still composes via `OperatorSum`. The architectural commitment is ONE algebraic FORM, not that every factor is a clean `TensorProductOperator`. |
| 2 | Python-level fold overhead in `TensorProductOperator.apply` (3 function calls per factor) hurts Krylov hot path > 5%. | Perf regression gate at §5.4. If exceeded, investigate `_build` fusion (collapse `(A & I & I) ∘ (B & I & I) → (A ∘ B) & I & I`) or factor caching. |
| 3 | Bit-identity breaks at FP-non-associativity ULP because `TensorProductOperator.apply` reorders reductions vs. the legacy fused einsum. | Per `vv-principles` §"Bit-identity vs principled-equivalence": acceptable if (a) the new path is principled (each factor is a named operator with definite physical meaning), (b) verified against an independent reference (k_∞ closed-form, L1 MMS), (c) drift is dimensionally explainable (`reduction_depth × ULP`). Document the contract relaxation at the affected regression snapshot. |
| 4 | T.3 scattering rewire touches the R · Λ · M pipeline used by EVERY iteration of EVERY SN solve. Bugs are expensive. | Test-architect proactive dispatch BEFORE T.3 implementation (per `subagent-handoff-protocol`'s proactive trigger for operator-algebra carves crossing subsystem boundaries). |
| 5 | Wave T's substeps depend on Depth B's typed Fields (T.3 needs HarmonicMomentField from D-E, T.4 needs AngularFlux from D-H). If Depth B's later steps slip, Wave T slips with them. | Sequential dependency; no shortcut. Communicate slippage up the parent-plan chain. |

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

When Wave T's T.1 through T.5 commit, the following invariants hold. **Items revised post-T.3/T.4 implementation to reflect the OperatorSum/bespoke-leaf deviations from the original SOTP aspiration** — see §6 T.3 + T.4 for the MA-Q1 master condition rationale.

1. **All four operator classes touched** (`SNBoundaryRealizer`-produced BCs, `FissionOperator`, `ScatteringOperator`, `StreamingOperator`) produce operator-algebra-typed results — but the algebraic SHAPE varies by physics coupling:
   - **`TensorProductOperator` (clean SOTP / TP)** — BC realizers (T.1, 5 instances: vacuum/specular/white/albedo/periodic); fission rank-1 TP (T.2: `RankOneOperator(χ, νΣ_f, axis=0) & IdentityOperator()`).
   - **`OperatorSum` of bespoke leaves (NOT SOTP)** — scattering kernel (T.3: `OperatorSum` over per-ℓ `_PerLegendreOrderScattering` summands, per Q6 honest fallback); streaming `M_spatial` (T.4b: `_MSpatialOperatorSum` of `_SpatialSweepDirection(±1)` leaves, per MA-Q1 master condition).
   - **Bespoke `LinearOperator` leaf (NOT TP, NOT OperatorSum at leaf level)** — `AngularRedistributionOperator` (T.4c: M-M half-grid recurrence wrapped as a single leaf, per Q2 = (iii)).
   - **`PrescribedInflow`** — affine-source BC, deferred to Wave O (linear-tensor-product form doesn't accommodate the affine shift). See §6 T.1.
   - `assert_separable` invariant fires ONLY on the TP/SOTP-shaped subset (BCs + fission); it is structurally inapplicable to bespoke leaves and pure-`OperatorSum` summands.
2. **L1 MMS gates** stay green. The 10 pre-existing DD-regression failures stay at the same failure set. The new algebraic-identity gates (§5.3) pass — `M_spat + M_ang == (L+C)` at principled-equivalence ULP.
3. **Bit-identity** holds at the matvec value level for the touched operators on the production hot path (T.4 verified bit-exact pre/post on slab + sphere + cylinder + 2-D Cartesian via pre-T.4 snapshots `cb18fdb`). For algebra-decomposition checks (M_spat+M_ang ≡ (L+C) for curvilinear) the principled-equivalence three-criteria test passes at ~16·ULP (`vv-principles` §"Bit-identity vs principled-equivalence").
4. **Performance** is within ≤ 5% of pre-Wave-T baseline on the 1-D slab Krylov benchmark. **Confirmed post-T.4c**: median 1.04× (3.97% slowdown), p95 0.998× — under threshold (`90e7d4e` close-out measurement).
5. **`TensorProductOperator`** has FIVE+ production consumers (4 BC kinds from T.1 + fission rank-1 TP from T.2 = 5; the specular BC carried over from D-B+1 raises the count to 6 if counted as a Wave T consumer too). **`SumOfTensorProductsOperator` has ZERO production consumers** — both T.3 and T.4 shipped as `OperatorSum` over bespoke leaves per the MA-Q1 master condition. The aspirational §15.2 SOTP form is contradicted by coupled-physics primitives (per-material XS, sequential WDD, M-M half-grid recurrence) in two of the three originally-SOTP-targeted substeps. **Wave O design SHOULD NOT assume universal SOTP-ability** of the operator algebra.

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
- **L1 primitives** (operators shipped Wave 0 commit `bc1253e`; space shipped Depth B D-B commit `c2f968a`; sole production consumer today is the D-B+1 specular BC at `boundary_realizer.py:164-166`):
  - `orpheus/numerics/operator.py:1058-1177` — `TensorProductOperator`.
  - `orpheus/numerics/operator.py:1179-1269` — `SumOfTensorProductsOperator`.
  - `orpheus/numerics/operator.py:1257-1269` — `assert_separable()`.
  - `orpheus/numerics/space.py:301` — `TensorProductSpace`.
- **Cardinal rules from memory**: `feedback_unify_after_two_instances`, `feedback_no_method_implementer_for_surgical_carves`, `feedback_architecture_forward_not_legacy_fit`, `feedback_elegance_causes_collapse`, `default-test-mode-is-optimize`.
