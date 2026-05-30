# Wave O — Operator-Role Typing (BulkOperator / FullOperator / BoundaryOperator) + Typed Sources / Residuals

**Branch:** `refactor/moment-space-and-layering`
**Worktree:** `.claude/worktrees/moment-space-and-layering/`
**GitHub Issue:** [#208](https://github.com/deOliveira-R/ORPHEUS/issues/208)
**Phase status:** PENDING. Lands AFTER Wave T (`wave_t_tensor_network.md`). Precedes parent-plan steps P3.4 + P3.6.

**Date filed:** 2026-05-28 (during D-H.1b). Plan stub: 2026-05-30 (this file).

**Status:** STUB — architectural commitment is settled (see §1 and the project memory `[[project_wave_o_operator_algebra]]`); detailed substep designs land when each substep starts (depth-on-demand per `[[feedback_no_method_implementer_for_surgical_carves]]` — main-agent direct authorship with turn-by-turn user steering).

---

## 0. Pickup checklist (read first)

If you are picking this plan up in a fresh session:

1. **Read this plan top-to-bottom.** No section is optional.
2. **Verify Depth B + Wave T completed**: `git log --oneline -30` should show the D-A through D-K commits AND the Wave T commits (T.1–T.5; see `wave_t_tensor_network.md` §6 for the substep ledger).
3. **Read the GitHub issue**: `gh issue view 208` — the original ask + the post-D-H.1b conversation that filed it.
4. **Read these grand-report sections**:
   - §5.5 — Field hierarchy (the Flux / Source / Residual roles that need typed treatment).
   - §5.7 — Operator hierarchy (the Bulk / Full / Boundary classification at the type system).
   - §16A.10 (lines 3142–3197) — BC as tensor network; the boundary-only operator's algebraic shape.
   - §32.4–32.6 — Space / Field / Operator specs that Wave O instantiates the missing Protocol surface for.
5. **Read the existing L1 operator algebra**:
   - `orpheus/numerics/operator.py` — `LinearOperator` Protocol, `OperatorSum`, `OperatorProduct`, capability dispatch.
   - `orpheus/sn/operator.py` — the SN operator leaves (`StreamingOperator`, `CollisionOperator`, `InvertibleOperator`).
   - `orpheus/geometry/boundary.py` — the boundary-realizer family (`SpecularBoundaryOperator`, `VacuumBoundaryOperator`, etc.).
6. **Read `[[project_wave_o_operator_algebra]]`** memory and `[[project_transport_state_container]]` memory.
7. **Pick up at the leftmost incomplete step in §6** below.

The `[[feedback_no_method_implementer_for_surgical_carves]]` rule applies: this is the main agent's work with turn-by-turn user steering. Do NOT batch via method-implementer.

---

## 1. The principle

**Operators on `(BulkSpace ⊕ BoundarySpace)` have a natural 2×2 block structure that the type system currently erases.** Today, every operator inherits from a single `LinearOperator` Protocol regardless of which blocks it activates. This conflates three distinct algebraic shapes:

- **Bulk-only** operators (`C` = collision, `S` = scattering, `F` = fission): only the `L_bb` block is non-zero. No boundary action. Action: `(bulk, boundary) → (L_bb·bulk, 0)`.
- **Full** operators (`L` = streaming): all four blocks live. Reads upstream boundary, produces downstream boundary. Action: `(bulk, boundary) → (L_bb·bulk + L_bs·boundary, L_sb·bulk + L_ss·boundary)`.
- **Boundary-only** operators (`BoundaryRealizer` family — specular, vacuum, white, albedo, periodic): only the `L_ss` block is non-zero. Action: `(bulk, boundary) → (0, L_ss·boundary)`.

**Asymmetry 2 — Source / Residual are full fields, not bulk-only.** The deeper insight (D-H.1b conversation, 2026-05-28): a source has both a bulk part (volumetric `Q(r, Ω, g)`) AND a boundary part (prescribed incoming flux `ψ|_{Γ_-}`). Vacuum BCs make the boundary part zero, but MMS with non-trivial inflow wants typed boundary sources. Similarly, residuals from `(L+C)ψ − q` have both bulk and boundary parts — the streaming matvec already produces both today but stuffs them into an `AngularFlux + BoundaryFlux` pair without a typed wrapper.

So `AngularSource`, `AngularResidual` should be **TimedFullField-shaped types** (with explicit bulk + boundary members), not bulk-only `Field` types. The current implicit-zero-boundary pattern shipped in D-H.1b is a placeholder; Wave O lifts it to typed.

### 1.1 What Wave O closes

| Current state (post-Depth-B) | Post-Wave-O state |
|---|---|
| All operators inherit `LinearOperator` regardless of block structure | `BulkOperator`, `FullOperator`, `BoundaryOperator` Protocols at L1 distinguish the three algebraic shapes |
| `OperatorSum.apply` calls every leaf with the same `TimedFullField`; bulk-only leaves return implicit-zero boundary | `OperatorSum.apply` dispatches by Protocol: bulk-only leaves see `bulk`; full leaves see `(bulk, boundary)`; boundary-only leaves see `boundary` |
| `AngularSource` and `AngularResidual` don't exist as distinct types — `AngularFlux` is the shape container for all three roles | `TimedFullAngularSource`, `TimedFullAngularResidual` ship as TimedFullField composites; the `(role × storage) = 12-cell` matrix from Issue #205 becomes machine-checkable |
| MMS machinery emits `AngularFlux` as source; the boundary part is implicit-zero | MMS emits typed `TimedFullAngularSource` whose boundary part is non-zero for non-vacuum BCs |
| `(L+C).apply(ψ)` returns `TimedFullField` with both bulk and boundary residual; consumers must know the type erasure | The return type IS `TimedFullAngularResidual`; the operator's codomain explicitly types the algebraic structure |
| Gate 1.3 (apply ↔ apply_transpose reciprocity) is xfail-strict pending Wave O | Gate 1.3 lands: every operator's adjoint is built analytically by Protocol-dispatched `.H` propagation through `OperatorSum` |

### 1.2 The "two architectural asymmetries" framing

The user's diagnosis from the D-H.1b conversation (2026-05-28, recorded in `[[project_wave_o_operator_algebra]]`):

> **Asymmetry 1** — operators have implicit bulk/boundary character.  No type distinguishes their character.
>
> **Asymmetry 2** — Source / Residual are full fields, not bulk-only.  Sources and residuals have bulk AND boundary parts.

Both are dimensional-typing sins (Issue #205 scope). Wave O makes them both machine-checkable.

### 1.3 Dependency on Wave T

Wave O sits on top of Wave T (`wave_t_tensor_network.md`). Wave T factors operators into tensor products `A = A_x ⊗ A_ω ⊗ A_g`. Wave O classifies each TENSOR FACTOR (and each resulting `TensorProductOperator` / `SumOfTensorProductsOperator`) as Bulk / Full / Boundary.

Specifically:
- **Boundary realizers** (T.1: vacuum, periodic, white, albedo as `K ⊗ I ⊗ I` tensor networks) are `BoundaryOperator` instances.
- **Fission** (T.2: `F = χ ⊗ νΣ_f`) is a `BulkOperator`.
- **Scattering** (T.3: `Σ_ℓ Σ_ℓ ⊗ A_ℓ ⊗ G_ℓ`) is a `BulkOperator` (no boundary action; the moment-folding is bulk-only).
- **Streaming** (T.4: `L_spatial + L_angular_redist`) is a `FullOperator` — its `L_spatial` factor reads/writes face traces, its `L_angular_redist` factor is bulk-only.

Wave O cannot land before Wave T because the Protocol classification needs the factored shape to type each factor correctly. Doing Wave O first would force the classification on the un-factored flat-axis legacy operators, which is the wrong shape.

### 1.4 Expected complications

- **`OperatorSum.apply` dispatch** currently consumes a single `TimedFullField` and routes to every leaf. Post-Wave-O, the dispatch must split the composite into `(bulk, boundary)` and route to bulk-only / full / boundary-only leaves per their Protocol — without breaking the algebra closure laws `(A + B).H = A.H + B.H` and `(A ∘ B).H = B.H ∘ A.H`.
- **The codomain space** of a full operator is `BulkSpace ⊕ BoundarySpace`. Wave O may need `DirectSumSpace` (deferred to P3.6 in Depth B) earlier than planned, or a lightweight ad-hoc shape until P3.6 lands. The `[[feedback_unify_after_two_instances]]` rule applies: don't ship `DirectSumSpace` until P3.6 also wants it (flux ⊕ precursors).
- **Backward compatibility with implicit-zero**: D-H.1b's implicit-zero-boundary pattern works because `TimedFullField` defaults boundary to zero. Wave O can keep the implicit-zero shorthand at the leaf level while adding typed `AngularSource` / `AngularResidual` at the composition level. Both should round-trip.
- **`apply_transpose` Protocol completeness**: Gate 1.3's xfail-strict status reflects that the analytical-adjoint propagation through `OperatorSum.H` is not yet wired for the operator algebra. Wave O lands this: each Bulk / Full / Boundary Protocol declares its own `.H` analytic, and `OperatorSum.H` composes by Protocol.

These complications are deferred to "discover during execution" — the architectural commitment to Protocol-classified operators is non-negotiable; the concrete shape of each substep refines when the substep starts.

---

## 2. Dependencies

Wave O cannot start until **both Depth B AND Wave T complete**. Specifically:

- Depth B D-A through D-K: the typed `Field` / `TimedFullField` substrate.
- Wave T T.1 through T.5: every production operator factored as `TensorProductOperator` / `SumOfTensorProductsOperator`. The Wave O Protocol classification names a factor's algebraic role; that role is only well-defined once factoring is in place.

Out-of-scope for Wave O:
- **`DirectSumSpace`** — deferred to P3.6 per Depth B §11.1 invariant #13 and `[[feedback_unify_after_two_instances]]`. Wave O may surface a second consumer (the operator codomain `BulkSpace ⊕ BoundarySpace`); if so, the unification with P3.6's `flux ⊕ precursors` lands at P3.6.
- **Issue #201** — `split AngularFlux from AngularSource / AngularResidual` — the role-typing slice. Wave O lands its bulk-side; the angular-storage slice of the 12-cell matrix lands at #201 close-out.
- **Issue #205** — the full 12-cell matrix. Wave O is a foundational layer; #205 closes when every cell ships its concrete type.
- **MoC and CP adoption** — the operator Protocols ship at L1 and any method can consume them. MoC / CP migrations land in their own waves.

---

## 3. Substep ledger

Each substep has a single named gate and a single commit boundary (1–3 commits per substep).

| Step | Scope | Gate | Status |
|---|---|---|---|
| **O.1** | `BulkOperator` / `FullOperator` / `BoundaryOperator` Protocols at L1 in `orpheus/numerics/operator.py`. Each Protocol declares: `apply(TimedFullField) → TimedFullField` semantics; capability set; analytic `.H` adjoint. | All existing operators retain runtime behaviour; capabilities API gains the Protocol classification. | PENDING |
| **O.2** | `OperatorSum.apply` / `.H` dispatch by Protocol: bulk-only leaves receive `bulk` only; full leaves receive `(bulk, boundary)`; boundary-only leaves receive `boundary` only. Algebra closure laws `(A + B).H = A.H + B.H` and `(A ∘ B).H = B.H ∘ A.H` hold by construction. | `test_phase_c_gates.py` Gate 1.3 (apply ↔ apply_transpose reciprocity) flips from xfail-strict to passing across slab + sphere + cylinder + 2-D Cartesian. | PENDING |
| **O.3** | `TimedFullAngularSource` and `TimedFullAngularResidual` typed wrappers — TimedFullField composites with explicit `bulk + boundary` members. MMS machinery emits typed sources for non-vacuum BCs. | New L0 invariants: source's `boundary` part is non-zero IFF the BC is non-vacuum; residual's `boundary` part equals the streaming matvec's face residual exactly. | PENDING |
| **O.4** | Migrate existing consumers: `solve_sn_fixed_source` accepts `TimedFullAngularSource`; `(L+C).apply` returns `TimedFullAngularResidual`; the implicit-zero-boundary shorthand from D-H.1b retires at the consumer boundary (not at the operator leaf). | All 12-file operator-algebra-core tests stay green; cardinality preserved; new test pin for non-vacuum boundary-source MMS lands. | PENDING |
| **O.5** | Retire the implicit-zero-boundary shortcuts shipped in D-H.1b. Each bulk-only leaf's `apply(TimedFullField) → TimedFullField` returns its result via the Protocol-classified `BulkOperator` API; the `OperatorSum.apply` dispatch fills the implicit-zero boundary at the algebra level, not at the leaf. | Final clean-up: zero `# implicit-zero boundary` comments in production; the typed contract is universal. | PENDING |

Each substep follows the veto pattern (Lessons L12 + L13 + L17 + L18 + L20) established by Depth B Wave D-I / D-J.

---

## 4. Verification strategy

### 4.1 New tests added by Wave O

- `tests/numerics/test_operator_protocols.py` — Protocol conformance (BulkOperator / FullOperator / BoundaryOperator); each existing operator advertises the correct Protocol.
- `tests/numerics/test_operator_sum_protocol_dispatch.py` — `OperatorSum.apply` dispatches by Protocol; bulk + boundary parts route correctly.
- `tests/numerics/test_typed_source_residual.py` — `TimedFullAngularSource` / `TimedFullAngularResidual` invariants (boundary part non-zero iff non-vacuum BC).
- `tests/sn/test_mms_non_vacuum_boundary.py` — MMS with non-vacuum BC; verifies the boundary part of the typed source matches the analytical inflow.

### 4.2 Tests that MUST stay green at every commit

- 12-file operator-algebra-core gate (the post-D-J baseline: 286 passed, 13 xfailed).
- Resolution A bit-exact decomposition gate (`test_streaming_operator_decomposition.py`).
- L1 MMS pins (Test 3.1 in `test_2d_l2_matvec_correctness.py`; curvilinear MMS in `test_mms_curvilinear.py`).
- L0 streaming-equilibrium anchor (`tests/sn/test_phase_c_gates.py` Gate 1.1).
- Gate 1.3 (apply ↔ apply_transpose reciprocity) — currently xfail-strict; FLIPS GREEN at O.2 (the load-bearing Wave O deliverable).

### 4.3 Bit-identity judgment per step

- **O.1**: Protocol classification only — runtime behaviour identical. Bit-identity required.
- **O.2**: Dispatch reorganisation — bit-identity expected on existing leaves; if FP-non-associativity drift surfaces, the three-criteria rule of `vv-principles` §"Bit-identity vs principled-equivalence" applies.
- **O.3 / O.4 / O.5**: Typed-wrapper introduction at the consumer boundary; the bulk part of every result MUST be bit-identical to pre-Wave-O. The boundary part is a NEW component (previously implicit-zero); its non-zero values are NEW invariants pinned by O.3's tests.

---

## 5. Exit Route

When Wave O completes:

- **Gate 1.3 flips green** — operator-algebra reciprocity is gated to round-off across all geometries.
- **Issue #208 closes.**
- **Issue #205 becomes much closer to closeable** — the operator-side typing lands; the field-side (`AngularFlux` / `ScalarFlux` / `HarmonicMomentField` split into Flux / Source / Residual) is the remaining work (#201 + scattered other tickets).
- **P3.4 (Problem/Solver split) becomes implementable on top of the role-typed algebra.** `CriticalityProblem.loss: FullOperator` and `CriticalityProblem.fission: BulkOperator` become machine-checkable type signatures.

The branch state at Wave O completion: `refactor/moment-space-and-layering` ready for P3.4 sequencing.

---

## 6. Cross-references

- **Parent plan:** `.claude/plans/moment_space_and_layering_plan.md` — Phase 3 sequencing (post-2026-05-30 update: `Depth B → Wave T → Wave O → P3.4 → P3.6`).
- **Depth B plan:** `.claude/plans/depth_b_field_on_function_space.md` — the typed Field substrate Wave O extends; §11.3 documents the post-Wave-T → Wave O ordering.
- **Wave T plan:** `.claude/plans/wave_t_tensor_network.md` — the tensor-product factoring Wave O classifies.
- **Project memory:** `[[project_wave_o_operator_algebra]]` — the original 2026-05-28 decision context.
- **Issue #208:** `gh issue view 208` — the ticket.
- **Issue #205:** `gh issue view 205` — the full cross-method field architecture matrix Wave O lands the operator-side of.
- **Issue #201:** `gh issue view 201` — the angular-storage role-typing slice (Wave O's downstream).
- **Grand report:** `.claude/plans/neutron_transport_grand_report_v3.md` §5.5 (Field hierarchy), §5.7 (Operator hierarchy), §32.4–32.6 (specs).
- **Lessons:** L17 (convention crosswalk before carve — applies to the Protocol dispatch carve at O.2), L18 (Pattern 7 producer-side normalisation — applies to the implicit-zero-boundary collapse at O.5).
