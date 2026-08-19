# Verification plan — #226 B4 operator-generics type-only carve

**Carve scope (verbatim).** Re-type the `LinearOperator[V]` generic family so
concrete operators are statically assignable and the capability surface
type-checks, with NO numerical/runtime behaviour change. PURE RE-TYPING.
Invariant: every operator's `apply`/matvec output, `.H` (G-adjoint) output,
`capabilities` set, `block_role`, and `OperatorSum`/`ScaledOperator` composition
results stay **bit-identical** (or principled-equivalent per `vv-principles`
where FP-reassociation is unavoidable — NOT expected anywhere in B4).

**Claim layer (§1.5).** This is a *foundation* carve — every gate is a
**software-invariant** claim (capability sets, block-role partition, adjoint
reciprocity identity, matvec≡sweep twin), NOT a convergence/flux/eigenvalue
claim. The reference pillar for the value-bearing gates is **closed-form**
(`k_inf = νΣ_f/Σ_a`, `φ = Q/Σ_t`, dense `np.linalg.solve`, hand-computed
outer products) and **structural identity** (G-reciprocity
`⟨Aψ,φ⟩_G = ⟨ψ,A.Hφ⟩_G`). Structural independence holds: the references are
dense-NumPy / closed-form, never another ORPHEUS operator. The L1 keff/MMS
gates ride ON TOP as the end-to-end safety net.

**The discipline that governs B4 (§0.5):** a type-only carve is "done" only
when, for every gate that claims to protect a runtime contract, the gate is
proven ABLE to red under the canonical `python -O -m pytest`. A type-only
change cannot move a value — so the gates must be the ones already sensitive to
the values, and the carve's net is the PROOF (via mutation) that they bite, not
the green run.

**The error oracle.** `npx pyright --outputjson orpheus/`. The B4 cluster is the
4 root-cause families below. Pyright on the 5 carve-target files currently
reports 28 errors; **18 are the B4 cluster** (RC1=4, RC2=3, RC3=8 incl.
boundary_realizer.py:260 `IncomingSourceOperator` return, RC4=3), the other 10
are adjacent out-of-scope reds (iteration.py:228 scipy `LinearOperator` ctor —
counts ×3 same line, :527/:568/:825/:1029; boundary_operator.py:133/:261) the
carve must NOT silence with `# type: ignore` but also is NOT chartered to fix —
leave them, do not regress the ratchet on them. The carve obeys: **NO `# type: ignore`** (principled
typing only — anti-pattern #19).

---

## 0. The B4 root causes mapped to pyright reds (the carve's surface)

| RC | What changes | pyright reds (file:line [rule]) |
|----|--------------|---------------------------------|
| **1** `block_role` ClassVar→instance | `LinearOperatorMixin.block_role: ClassVar[Optional[BlockRole]] = None` is a `ClassVar`, but `_AdjointOperator`/`OperatorSum`/`ScaledOperator`/realizer assign `self.block_role = …`. | operator.py:623, :736, :867 [reportAttributeAccessIssue]; boundary_realizer.py:125 |
| **2** `capabilities` ClassVar-vs-instance | mixin declares `capabilities: frozenset[str]` (instance); leaves redeclare as class attr / `ClassVar`. | angular_operator.py:102, :270 [reportIncompatibleVariableOverride]; boundary_operator.py:137 (`@property` override) |
| **3** bare `LinearOperatorMixin` (no `[V]`) | `TensorProductOperator(LinearOperatorMixin)` + `IncomingSourceOperator` → `LinearOperator[Unknown]`, rejected by `LinearOperator[V]`-typed params. | boundary_realizer.py:176/:206/:207/:223/:224/:240/:252 [reportArgumentType], :260 [reportReturnType] |
| **4** `.solve` capability-gated, not on Protocol | `solve` lives only on some subclasses; `iteration.py` calls `.solve` on a `LinearOperator[V]`-typed value. | iteration.py:459/:467/:468 [reportAttributeAccessIssue] |

The carve is a re-typing, so the **runtime** for each of these is already
correct and pinned — the question this plan answers is: *which existing test
pins each contract, can that test red, and what gap must be filled before the
carve so a silent runtime drift cannot slip through a re-typing.*

---

## 1. Contract-pinning test inventory

All operator-algebra contract tests are `@pytest.mark.foundation` (software
invariants, no equation `:label:`). Inventory by contract:

### 1a. Capability sets + closure laws (RC 2, RC 4)

| Test (file::name) | Asserts | Class |
|---|---|---|
| `tests/numerics/test_operator.py::test_sum_solve_does_not_propagate` | `CAP_SOLVE ∉ (A+B).capabilities` | foundation |
| `…::test_sum_apply_propagates_with_both` / `…test_sum_transpose_propagates_with_both` / `…test_sum_transpose_drops_when_one_lacks` | Sum capability closure (apply needs both; transpose needs both) | foundation |
| `…::test_product_solve_propagates_with_both` / `…test_product_solve_drops_when_one_lacks` | Product `solve` needs BOTH (the `(AB)⁻¹=B⁻¹A⁻¹` law) | foundation |
| `…::test_product_apply_propagates` / `…test_product_transpose_propagates_with_both` | Product apply/transpose closure | foundation |
| `…::test_scaled_preserves_all_capabilities` / `…test_scaled_apply_only_stays_apply_only` | Scaled forwards the full survivor set | foundation |
| `…::test_identity_full_capabilities` / `…test_zero_lacks_solve` | Identity = {apply,solve,transpose}; Zero ∌ solve | foundation |
| `…::test_sum_rejects_missing_apply_at_composition` / `…test_product_rejects_…` / `…test_scaled_rejects_…` | `MissingCapability` raised AT COMPOSITION, matched on `"apply"` | foundation |
| `…::test_transport_eigenvalue_algebra_smoke` | `(L−S)⁻¹F`-style smoke: apply ∈ caps, solve ∉ caps for the sum | foundation |
| `tests/numerics/test_tensor_product_operator.py::TestTensorProductCapabilities::*` (4) | TP caps = **intersection** of factors (apply-only, full, solve-through-diagonals, solve-revoked-by-zero) | foundation |
| `…::TestSumOfTensorProducts::test_solve_not_advertised` / `…test_apply_transpose_propagates` | SoTP caps | foundation |
| `…::TestTensorProductRequiresApply::test_construction_rejects_non_apply_factor` | TP ctor `MissingCapability` | foundation |
| `…::TestRankOneOperator::test_capabilities_apply_only` | RankOne = {apply} | foundation |

**Per-leaf capability pins** (the per-leaf set, RC 2's payload):

| Leaf | Test::name asserting `capabilities == frozenset({…})` |
|---|---|
| Fission `F` | `tests/sn/operators/test_fission_operator.py:71` (`{CAP_APPLY}`) + kernel `:415` |
| Scattering `S` | `tests/sn/operators/test_scattering_operator.py:137` (`{CAP_APPLY}`, ∌ solve, ∌ transpose) |
| Streaming `L` | `tests/sn/operators/test_streaming_operator.py:176` (`{CAP_APPLY,CAP_APPLY_TRANSPOSE}`, ∌ solve) |
| Collision `C` | `tests/sn/operators/test_collision_operator.py:148` (explicit set) |
| AngularAverage | `tests/sn/operators/test_angular_average_operator.py:153` (`{CAP_APPLY}`) — **the L102 `ClassVar` leaf (RC 2)** |
| `SNBoundaryOperator` B | `tests/sn/operators/test_sn_boundary_operator.py:112/:200/:241/:242` (apply; transpose present/absent by law) — **the L137 `@property` capabilities (RC 2)** |
| `(L+C)` Invertible | `tests/sn/operators/test_streaming_operator.py:331/:345` (apply + solve) |

### 1b. `block_role` join + partition (RC 1)

| Test (file::name) | Asserts | Class |
|---|---|---|
| `tests/sn/operators/test_operator_block_role.py::TestBulkLeaves::*` (3) | C/S/F → `BlockRole.BULK`, ∉ Full | foundation |
| `…::TestFullLeaves::test_streaming_is_full` / `…test_invertible_L_plus_C_is_full` | L → FULL; `(L+C)` → FULL (DERIVED) | foundation |
| `…::TestBoundaryLeaves::test_realized_bc_advertises_boundary_operator` (param ×9 laws) | **Realizer output retains `block_role=BOUNDARY`** — exercises the `_as_boundary` instance-stamp on `TensorProductOperator`/`ScaledOperator` (RC 1 ∩ RC 3) | foundation |
| `…::TestBoundaryLeaves::test_prescribed_inflow_is_not_boundary_operator` | `PrescribedInflow` source carries `block_role None` | foundation |
| `…::TestBoundaryLeaves::test_mesh_bc_forwards_boundary_role` | `sn.bc[face]` shim forwards BOUNDARY | foundation |
| `…::TestComposerRoleDerivation::test_scaled_preserves_role` | `(−C)`/`(−B)`/`(2L)` keep role (RC 1 on ScaledOperator) | foundation |
| `…::TestComposerRoleDerivation::test_sum_join_is_union_of_blocks` | `BULK⊔BULK=BULK`, `FULL⊔BULK=FULL`, **`BULK⊔BOUNDARY=FULL`** (the `L+C−B` join) | foundation |
| `…::TestComposerRoleDerivation::test_adjoint_preserves_role` | `L.H`/`B.H`/`C.H` and `(L+C−B).H` preserve role (RC 1 on `_AdjointOperator`) | foundation |
| `…::TestComposerRoleDerivation::test_unclassified_operand_propagates_none` | `None` propagates through a sum | foundation |
| `tests/numerics/test_operator_protocols.py::TestBlockRoleEnum/TestIsinstanceMarkers/TestGenericOperatorsAreUnclassified::*` | the MECHANISM: 3-role enum, value-based `isinstance`, partition exclusivity, generic-default-`None`, `OperatorSum`-default-`None`-at-O1 | foundation |

### 1c. `.H` G-adjoint correctness (downstream of RC 1, RC 4)

| Test (file::name) | Asserts | Class |
|---|---|---|
| `tests/numerics/test_operator.py::test_adjoint_euclidean_identity` / `…test_adjoint_weighted_identity` | `⟨Ax,y⟩=⟨x,A*y⟩` Euclidean + weighted (the defining adjoint identity) | foundation |
| `…::test_adjoint_swaps_domain_and_codomain` / `…test_H_aliases_adjoint` | domain↔codomain swap; `A.H == A.adjoint()` | foundation |
| `…::test_transpose_distributes_over_sum` / `…test_transpose_reverses_product` | `(A+B)ᵀ`, `(AB)ᵀ` value laws | foundation |
| `tests/sn/operators/test_g_adjoint_reciprocity.py::test_g_adjoint_reciprocity_full_block` (param) | `⟨Aψ,φ⟩_G=⟨ψ,A.Hφ⟩_G` for `A=L+C−B`, metric evaluated INDEPENDENTLY | foundation |
| `…::test_full_field_space_metric_matches_independent_reference` | G-metric population vs independent `V·wₙ ⊕ |Ω·n|·wₙ` | foundation |
| `…::test_wrong_trace_metric_breaks_reciprocity` (L11 control) | dropping `|Ω·n|` MUST break reciprocity (negative control) | foundation |

### 1d. apply/matvec value + composition correctness (the bit-identity ground)

| Test (file::name) | Asserts | Class |
|---|---|---|
| `tests/numerics/test_operator.py::test_sum_apply_matches_dense` / `…test_product_apply_matches_dense` / `…test_scaled_apply_matches_dense` | composite apply == dense reference | foundation |
| `…::test_product_solve_reverses_order` / `…test_scaled_solve_divides` | `(AB)⁻¹` order; `(αL)⁻¹=(1/α)L⁻¹` | foundation |
| `…::test_pow_three_matches_repeated_apply` + pow edge tests | `A**n` value | foundation |
| `tests/sn/operators/test_removal_form_matvec_sweep.py::test_invertible_apply_is_M_of_C_sigma_bit_identical` (param) | `(L+C).apply == loss_action(σ_r)` **bit-identical** (`assert_array_equal`) | `verifies("loss-rep-resolution-a")` |
| `…::test_invertible_apply_transpose_is_M_transpose_…_bit_identical` | transpose matvec bit-identical | same |
| `…::test_removal_form_matvec_sweep_roundtrip` (param, sphere excluded) | matvec≡sweep inverse-twin (`assert_allclose`) | same |
| `tests/numerics/test_tensor_product_operator.py::TestTensorProductApply::test_apply_against_einsum_reference` | TP apply vs `einsum` (RC 3 value) | foundation |
| `…::TestTensorProductAdjoint::test_adjoint_distributivity_diagonals` | `(A⊗B)*=A*⊗B*` (RC 3 adjoint) | foundation |
| `…::TestTensorProductFlattening::*` (2) | `(A&B)&C` flattens (RC 3 build) | foundation |
| `…::TestRankOneOperator::test_apply_matches_hand_computed_outer_product` | RankOne value | foundation |

### 1e. Realizer-output structure (RC 3 — bare `TensorProductOperator`)

| Test (file::name) | Asserts | Class |
|---|---|---|
| `tests/sn/operators/test_sn_boundary_realizer.py::*::test_vacuum_returns_tensor_product` (:132) | vacuum dispatch returns a `TensorProductOperator` | foundation |
| `…::test_specular_unit_albedo_returns_tensor_product` (:184) | specular α=1 returns TP | foundation |
| `…::test_specular_partial_albedo_returns_scaled_tensor_product` (:215) | α∉{0,1} returns `ScaledOperator(TP)` | foundation |
| `…::test_*_matches_hand_computed` (vacuum/specular/white) | realized-op apply value vs hand calc — the VALUE ground that the re-typed TP still acts correctly | foundation |

### 1f. `.solve` dispatch through iteration (RC 4 — end-to-end)

| Test (file::name) | Asserts | Class |
|---|---|---|
| `tests/numerics/test_iteration.py` (SourceIteration suite) | `SourceIteration(L,…).solve` calls `L.solve` each step; `MissingCapability` raised when `L` lacks `CAP_SOLVE` | foundation |
| `tests/sn/solve/test_si_single_primitive_contract.py` | the single-primitive SI contract (the `L.solve` path that RC 4 re-types) | foundation/L1 |

### 1g. The L1 end-to-end safety net (a silent operator regression surfaces here)

| Test (file::name) | Asserts | Level |
|---|---|---|
| `tests/sn/eigenvalue/test_keff_curvilinear.py::test_homogeneous_exact` (param) | `keff == k_inf` (closed-form, ≥2G branch via TestL2) — flows L,C,S,F,B,`.H` | L1 |
| `…::TestL2…::test_2g_heterogeneous_fuel_moderator` | 2G heterogeneous keff finite/positive (≥2G — anti-#3) | L2 |
| `tests/sn/verification/analytical/test_l1_standoff_slab_cylinder.py` | analytical standoff (fixed-source/keff through operators) | L1 |
| `tests/sn/verification/mms/test_mms_ld_slab.py` / `…test_mms_ld_2d.py` | MMS convergence through the operator chain (flux-shape, NOT eigenvalue) | L1 |
| `tests/sn/verification/analytical/test_prescribed_inflow_consistency.py` | boundary `B` consistency through realizer + algebra | L1 |

---

## 2. Classification per test (must-stay-bit-identical / principled-equiv / new)

**Default for a PURE RE-TYPING = must-stay-bit-identical.** A re-typing touches
no arithmetic; there is NO FP-reassociation anywhere in B4, so **zero tests are
re-classified to principled-equivalent**. If the carve produces ANY non-bit
change in §1d/§1e (the `assert_array_equal` matvec twins, the einsum/hand-calc
value pins), that is a BUG in the re-typing, not an acceptable FP drift — STOP.

| Bucket | Tests | Classification |
|---|---|---|
| Capability sets (§1a) | all | must-stay-bit-identical (set equality) |
| block_role (§1b) | all | must-stay-bit-identical (enum identity) |
| `.H` adjoint (§1c) | all | must-stay-bit-identical (`assert_allclose rtol≤1e-10`, identity-level) |
| apply/matvec value (§1d) | all | **must-stay-bit-identical** (`assert_array_equal` on the removal-form twins; `rtol` only where the existing pin already uses it) |
| realizer structure (§1e) | all | must-stay-bit-identical |
| solve dispatch (§1f) | all | must-stay-bit-identical |
| L1/L2 net (§1g) | all | must-stay-bit-identical (tolerances unchanged — a re-typing cannot move keff) |
| **NEW (§4)** | 3 gates below | new-test-needed BEFORE the carve |

---

## 3. Risk map — where a TYPE change could change RUNTIME (with the catching gate)

Each row: the type mechanism → the runtime drift it could silently cause → the
assertion that goes RED → whether that assertion can red under `-O`.

### RC 1 — `block_role` ClassVar→instance

- **(a) leaf default could change.** Moving `block_role` off `ClassVar` and the
  leaves' `block_role = BlockRole.X` class-attr declaration risk: a leaf that
  read its role off the class now reads `None` (or vice versa).
  - **Catches:** `TestBulkLeaves::test_scattering_is_bulk`/`…fission_is_bulk`
    assert `ScatteringOperator.block_role is BlockRole.BULK` **at the class
    level** (S/F need a `MaterialXSField` to instantiate, so the class-attr
    read is the pin — see the file's docstring justifying class-level ≡
    instance). `test_collision_is_bulk`/`test_streaming_is_full` assert on
    INSTANCES. If the carve makes the leaf default `None`, the BULK/FULL
    asserts red. **Reds under -O** (`assert C.block_role is BlockRole.BULK` is
    AST-rewritten in test files even under -O — but to be safe these use plain
    `assert`; the value is an enum identity, the assert fires).
- **(b) `_join_block_roles` result could change** if `OperatorSum.block_role`
  assignment is re-typed to read a class default instead of `_join_block_roles`.
  - **Catches:** `TestComposerRoleDerivation::test_sum_join_is_union_of_blocks`
    (`(L+C−B)` → FULL, `BULK⊔BOUNDARY=FULL`). The discriminating join.
- **(c) the instance stamp on a composite/realizer could be lost.** If
  `self.block_role = …` on `OperatorSum`/`ScaledOperator`/`_AdjointOperator`/
  `_as_boundary` is silently dropped to satisfy the type checker (e.g. by making
  `block_role` read-only), the realizer/composite loses its role.
  - **Catches:** `test_scaled_preserves_role`, `test_adjoint_preserves_role`,
    `test_realized_bc_advertises_boundary_operator` (×9 laws),
    `test_mesh_bc_forwards_boundary_role`.

### RC 2 — `capabilities` ClassVar normalization

- **A leaf's advertised set could be dropped/altered** if the re-typing changes
  `ClassVar[frozenset]` → instance and a leaf's set is lost, or the
  `@property capabilities` on `SNBoundaryOperator` is re-typed to a plain attr
  that returns a stale value.
  - **Catches:** the per-leaf `capabilities == frozenset({…})` pins (§1a leaf
    table). `AngularAverage` (L102 ClassVar) → `test_angular_average_operator.py:153`;
    `SNBoundaryOperator` (L137 property) → `test_sn_boundary_operator.py:112/:200/
    :241/:242` (transpose present/absent BY LAW — the property's law-dependent
    branch). Closure-law tests (§1a) catch a composite-level drop.

### RC 3 — parametrize `TensorProductOperator[V]`

- **Re-typing the bare mixin to `LinearOperatorMixin[V]` should change NOTHING
  at runtime** — it is purely the generic parameter. Confirm: TP apply value,
  caps intersection, adjoint distributivity, flatten, AND the realizer-output
  role-stamp all unchanged.
  - **Catches:** `TestTensorProductApply::test_apply_against_einsum_reference`,
    `TestTensorProductCapabilities::*`, `TestTensorProductAdjoint::*`,
    `test_vacuum_returns_tensor_product`, `test_realized_bc_advertises_boundary_operator`.

### RC 4 — declaring `solve` on the Protocol/mixin

- **Declaring `solve` (even abstractly) on `LinearOperator`/`LinearOperatorMixin`
  could change capability-gating or `MissingCapability` behaviour.** This is the
  SHARPEST risk: if `solve` becomes a mixin method (rather than capability-gated),
  an operator that should NOT advertise solve might now appear to have it, OR the
  `MissingCapability`-at-composition guard could be bypassed.
  - **Principled spelling (the carve MUST use, NOT `# type: ignore`):** declare
    `solve` on the Protocol as a method whose PRESENCE is the static contract but
    keep the RUNTIME gate (`_has(L, CAP_SOLVE)` → `MissingCapability`) intact —
    OR narrow `iteration.py`'s `L` parameter to a `SolvableOperator`
    Protocol (`LinearOperator[V]` + `solve`), so the call site types without
    weakening the leaf contract. The capability SET stays the single source of
    truth; the Protocol gains an OPTIONAL-method declaration only.
  - **Catches:**
    - `test_sum_solve_does_not_propagate` / `test_product_solve_drops_when_one_lacks`
      / `test_zero_lacks_solve` — if `solve` leaks onto every operator, these red
      (a sum would now advertise solve).
    - `test_sum_rejects_missing_apply_at_composition` (and the iteration
      `MissingCapability` test) — if the gate is bypassed, the negative test that
      expects the raise goes RED (no raise).
    - `tests/numerics/test_iteration.py` `MissingCapability`-on-no-solve case +
      `test_si_single_primitive_contract.py`.

**Under -O:** every value/identity gate above uses `assert <enum/set identity>`
in a TEST FILE (AST-rewritten under -O — fires) or `np.testing.assert_*` /
`pytest.raises` (function calls — fire under -O). The removal-form twins use
`np.testing.assert_array_equal`. No B4 gate relies on a production-code bare
`assert` (Mode 8 clean).

---

## 4. Gaps — contracts a re-typing could silently break, NOT pinned today

Three gaps. **Add all three BEFORE the carve** so the carve has a net.

### GAP-1 (RC 3 ∩ RC 1) — `TensorProductOperator` default `block_role is None`

`grep block_role tests/` over `TensorProduct*` → **NONE**. The realizer test
pins that a STAMPED TP carries BOUNDARY, but nothing pins that an UN-stamped TP
(the bare-mixin default) is `None`. If re-typing RC 3 accidentally gives
`TensorProductOperator` a non-`None` class default `block_role` (e.g. by
inheriting a wrong default during the `[V]` parametrization), an un-stamped TP
would silently classify as some role, and `_join_block_roles` over it would
mis-join — invisible today.

**New test** (`tests/numerics/test_operator_protocols.py`,
`TestGenericOperatorsAreUnclassified`):

```python
def test_tensor_product_operator_is_unclassified_by_default(self):
    """A bare TP (no realizer stamp) carries block_role None — the
    BOUNDARY role is an instance stamp (_as_boundary), never a TP default."""
    from orpheus.numerics.operator import DiagonalOperator, TensorProductOperator
    d0 = DiagonalOperator(np.ones(3), axis=0)
    d1 = DiagonalOperator(np.ones(4), axis=1)
    tp = TensorProductOperator((d0, d1))
    assert tp.block_role is None
    assert not isinstance(tp, BoundaryOperator)
```

Also add `IncomingSourceOperator` and `SumOfTensorProductsOperator` to the
generic-`None` sweep (both are bare-mixin, RC 3 candidates).

### GAP-2 (RC 2) — the two `ClassVar` leaves lack a per-leaf set + composite-survival pin

`AngularAverageOperator` (L102) IS pinned (`:153`). `IncomingSourceOperator`
(L270, the OTHER bare-mixin `ClassVar` leaf, RC 3) has **no `capabilities ==`
per-leaf pin** found in `tests/`. A re-typing that drops or alters its set is
invisible.

**New test** (`tests/sn/operators/test_angular_average_operator.py` or a
co-located `test_incoming_source_operator.py`):

```python
def test_incoming_source_operator_capabilities_apply_only(self):
    op = IncomingSourceOperator(...)   # minimal valid construction
    assert op.capabilities == frozenset({CAP_APPLY})
```

Plus a **composite-survival** pin for the property-based `SNBoundaryOperator`
caps (RC 2's nastiest re-typing target — `@property` → attr): assert that
`(L + C − B).capabilities` still carries `CAP_APPLY` and does NOT carry
`CAP_SOLVE` when `B` is the property-capabilities operator (the existing
streaming `:331/:345` invertible test covers `(L+C)` but not the `−B` arm).

### GAP-3 (RC 4) — `solve`-on-Protocol must NOT widen leaf capability advertisement

The closure tests pin that a SUM doesn't gain solve, but there is **no direct
pin that a bare leaf which lacks `CAP_SOLVE` is NOT statically/at-runtime
solvable AFTER the Protocol gains a `solve` declaration**. The risk is that
declaring `solve` on the Protocol makes consumers BELIEVE every `LinearOperator`
is solvable. Add a runtime guard pin:

**New test** (`tests/numerics/test_operator.py`):

```python
def test_apply_only_operator_solve_raises_or_absent_after_protocol_decl(
    matrix_apply_only,
):
    """RC4 guard: an apply-only operator must NOT acquire a working solve
    just because the Protocol now DECLARES one. Capabilities is the SoT."""
    assert CAP_SOLVE not in matrix_apply_only.capabilities
    # And a SourceIteration built on it raises MissingCapability at construction.
    with pytest.raises(MissingCapability):
        SourceIteration(matrix_apply_only)   # L lacks solve
```

This is the test that DIRECTLY guards the §3-RC4 sharpest risk: it reds if the
carve's principled `solve`-declaration accidentally makes the capability gate a
no-op.

---

## 5. Mutation-test recipes (prove the top gates bite — §0.5)

For the 4 highest-value gates, the EXACT mutation that MUST turn the gate RED
under `python -O -m pytest`. Mutate **in-process via monkeypatch** (the carve
holds uncommitted edits — NEVER `git checkout` per process-discipline). Run each
mutation, expect RED, restore.

### M1 — RC 1, the `_join_block_roles` join (gate: `test_sum_join_is_union_of_blocks`)

Monkeypatch `_join_block_roles` to return `a` unconditionally (drop the
`a if a is b else FULL` mix-rule):

```python
import orpheus.numerics.operator as op
monkeypatch.setattr(op, "_join_block_roles", lambda a, b: a)
```

Run `pytest tests/sn/operators/test_operator_block_role.py::TestComposerRoleDerivation::test_sum_join_is_union_of_blocks -O`.
MUST RED on `(C − B)` → would report BULK not FULL. Proves the join gate bites.

### M2 — RC 1, the realizer instance stamp (gate: `test_realized_bc_advertises_boundary_operator`)

Monkeypatch `_as_boundary` in `boundary_realizer` to skip the stamp:

```python
import orpheus.sn.boundary_realizer as br
monkeypatch.setattr(br, "_as_boundary", lambda o: o)   # drop block_role=BOUNDARY
```

Run the parametrized boundary test (×9). MUST RED (`op.block_role is None`,
`isinstance(op, BoundaryOperator)` False). Proves RC 1's instance-stamp on a
bare TP is pinned — the highest-value RC1∩RC3 intersection.

### M3 — RC 2, drop a leaf capability (gate: per-leaf `capabilities ==`)

Monkeypatch the leaf's `capabilities` to add a spurious `CAP_SOLVE`:

```python
monkeypatch.setattr(AngularAverageOperator, "capabilities",
                    frozenset({CAP_APPLY, CAP_SOLVE}))
```

Run `test_angular_average_operator.py:153`. MUST RED (set inequality). Proves
the per-leaf set pin bites — the RC 2 ClassVar-normalization guard.

### M4 — RC 4, neuter the `MissingCapability` gate (gate: GAP-3 new test + iteration test)

Monkeypatch `SourceIteration` to skip the `_has(L, CAP_SOLVE)` check, OR
monkeypatch `_has` to always return True:

```python
import orpheus.numerics.iteration as it
monkeypatch.setattr(it, "_has", lambda op, cap: True)
```

Run the iteration `MissingCapability`-on-no-solve test + GAP-3. MUST RED
(the expected `pytest.raises(MissingCapability)` no longer raises). Proves the
RC 4 capability-gate is not silently bypassable — the sharpest B4 risk.

### M5 (bonus) — RC 1/3, removal-form bit-identity (gate: `test_invertible_apply_is_M_of_C_sigma_bit_identical`)

This `assert_array_equal` twin is the value ground proving a re-typing of the
composite `(L+C)` did not perturb the matvec. Mutate the composite to scale the
output by `(1 + 1e-15)`:

```python
# in-process: wrap op.apply to multiply by (1+1e-15)
```

MUST RED (`assert_array_equal` is bit-exact). Confirms the bit-identity claim
of the whole carve is enforceable at ULP level.

---

## 6. Pre-carve checklist (the net)

1. Add GAP-1, GAP-2, GAP-3 tests; run `python -O -m pytest` — GREEN on the
   un-mutated tree (they pin current correct behaviour).
2. Run M1–M5 mutations IN-PROCESS — confirm each REDS the named gate under -O.
   If any stays green, the gate is a no-op for its claim → fix the gate BEFORE
   the carve.
3. Capture the pyright baseline: `npx pyright --outputjson orpheus/` — record
   the 18 B4-cluster errors (the carve drives these to 0) and the 10
   out-of-scope reds (the carve must NOT touch or regress).
4. Carve. After: `npx pyright` → 18 B4 reds gone, 0 new reds, ratchet not
   regressed; `python -O -m pytest tests/numerics tests/sn/operators
   tests/sn/eigenvalue tests/sn/verification` GREEN with NO tolerance change and
   NO `# type: ignore` added.
