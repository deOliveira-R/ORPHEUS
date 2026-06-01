# Typed Field Vocabulary + View-G Space Layer (#205 / #201; pre-#208)

**Branch:** `refactor/field-role-typing` (worktree child of `refactor/sn-operator-algebra` @ ad813fd).
**Approved plan (host copy):** `~/.claude/plans/glimmering-launching-lantern.md` — this is the durable in-repo copy.
**Issues:** #205 (cross-method field architecture — SN/Scalar portion), #201
(angular role split), #207 (space property-composition — applied), #197 (units
discipline — superseded), #208 (operator typing — **deferred**).
**Mode:** main-agent direct authorship, turn-by-turn user steering
(`feedback_no_method_implementer_for_surgical_carves`).

---

## Context — why

`AngularFlux` & friends are under-abstracted; the type system erases physics
distinctions. (1) Role conflation — the "dimensional sin": one type carries ψ
(`cm⁻²`), operator outputs / `q_ext` / residual (`cm⁻³`); epicentre
`numerics/iteration.py:455`, `sn/operator.py:1938`. (2) Machinery duplication —
`AngularFlux`/`PerOrdinateSource` redeclare mesh+shape+factories+`_check_partner`;
`BoundaryFlux`'s FaceLayout apparatus would be copy-pasted into new boundary
roles. (3) Three unreconciled boundary-space notions + a directional predicate
duplicated between `spaces/trace_space.py` and `sn/operator.py`.

Outcome: a machine-checkable `{Angular, Scalar, Boundary} × {Flux, Source,
Residual}` field vocabulary on a **View-G** space layer (space = pure geometry;
units on the quantity), with one unified `TraceSpace`. Foundation for #208.

## Locked design (agreed this session)

- **View-G.** `FunctionSpace` = geometry + operations (shape, inner product,
  `⊗`/`dual()`/future `⊕`), role-agnostic. **Units NOT on the space** — on the
  field (role-leaf `UNITS` class constant) and, in #208, on the operator
  (unit-gain). Geometric chain checked at compose/apply; dimensional chain via
  pint once at operator construction (Layer 2, `-O`-safe). `L`/`L⁻¹` are
  geometric endomorphisms on the bulk grid with a dimensional gain.
- **Storage-base + role-leaf (single inheritance, no role mixins):**
  ```
  Field
   ├─ BulkField (marker) ─ AngularField → {AngularFlux, AngularSource, AngularResidual}
   │                       ScalarField  → {ScalarFlux,  ScalarSource,  ScalarResidual}
   │                       MomentField  → HarmonicMomentField
   └─ BoundaryField (locus = trace storage, FaceLayout) → {BoundaryFlux, BoundarySource, BoundaryResidual}
  ```
- **Renames:** `PerOrdinateSource→AngularSource`, `IsotropicSource→ScalarSource`.
- **Strict named-composition:** cross-**role** arithmetic raises;
  `IterationResidual.from_balance(lhs, rhs)`. (Layer-1 class identity is the gate.)
  **RESOLVED (2026-05-31):** the #207 cross-**storage** same-role injection
  (`ScalarSource + AngularSource → AngularSource` via broadcast, B ⊃ A) STAYS
  implicit — it is the same role, one storage a canonical subspace of the
  other, and the algebra reads clean. Elegance beats strictness for this
  edge case; the strict discipline remains the default for genuine
  cross-*role* arithmetic. See `[[feedback_subspace_injection_elegance]]`.
- **Units via class constant** `UNITS`; Layer-3 `space.units` assert retires.
- **Composite** stays one `TimedFullField`, slots tightened to
  `bulk: BulkField, boundary: BoundaryField`.
- **TraceSpace reconciliation:** ONE concrete `TraceSpace(FunctionSpace)` —
  `shape=(total_size,)`, `layout: FaceLayout`, `sign(Ω·n)` predicate leaf-data,
  `inflow/outflow_indices_for_face` selectors. Retire `Inflow/OutflowTraceSpace`,
  `boundary_trace_space()` factory, `numerics/trace_space.py` shim. Bulk spaces
  stay flat `FunctionSpace` via factories.
- **Out of scope → #208:** Bulk/Full/Boundary operator Protocols; operator
  unit-gain; BC extraction `(L_full+C−S−F−B)ψ=q`; inflow/outflow projection
  *operators*; CP/MoC fields; `DirectSumSpace` (`⊕`).

## Execution ledger (gates from the test-architect verification spec)

**Regression wall:** SN operator-algebra-core 12-file gate = **347 pass / 0 fail
/ 0 xfail**; anchors 46 pass / **8 expected-xfail** (ERR-026; must stay xfail).
**Pre-existing red (NOT ours; pin at exactly 5, must not grow):**
`test_snmesh_realizer_wiring.py` ×3 + `test_bound_compat.py` ×2 (curvilinear-
realizer routing). A.3 touches that wiring — watch it.

| Step | Scope | Bit-identity? | Gate |
|---|---|---|---|
| 0a/0b | test-architect spec + explorer L20 audit | — | done |
| **A.1 ✅** | remove `units` from `FunctionSpace` (`space.py`: field, `__eq__`, `__repr__`, `_tensor_product_units`); drop Layer-3 assert in `field.py`; migrate unit-pins (3 in `test_space_algebra.py` + **5 in `test_field.py`** — grep found more than the spec listed); fix latent cm²→cm³ error in field.py docstring; scalar_flux.py docstring | **BI** (units `None` everywhere) | **DONE: 237 passed `-O`** (space+space_algebra+field+operator+transport+typed_sources) |
| **A.2** | one concrete `TraceSpace` (`spaces/trace_space.py`) | new+BI predicate | rewritten `test_trace_space.py` (L1 `sign(Ω·n)`) |
| **A.3** | migrate 5 consumers; delete `Inflow/Outflow`, factory, shim | BI | trace/realizer/method_space/btl/bound_compat; red set still =5 |
| **A.4** | consolidate inline `mu_x>±eps` masks in `operator.py` (≈8 copies; verify `eps`=`_TANGENTIAL_EPS`) | **BI (highest risk)** | Resolution-A `assert_array_equal`; 347 wall |
| **A.5** | re-home `BoundaryFlux` onto `TraceSpace`; FaceLayout field→space | BI | `test_boundary_flux.py`, `test_timed_full_field.py` |
| **B.1** | storage-base ABCs; dedup; re-parent | BI | leaf tests bit-identical |
| **B.2** | rename sources (1-cycle shim); re-parent; land atomically (cross-class `__add__` + 4 singledispatch guards) | BI | `test_typed_sources.py` (migrate names) |
| **B.3** | new leaves `Angular/ScalarResidual` (cm⁻³), `Boundary{Source,Residual}` (cm⁻²) | new | new construction/arith/cross-role-raise tests |
| **B.4** | `UNITS` class constant; Layer-1 = units gate | new | dimensionality + `-O` zero-cost tests |
| **B.5** | `from_balance`; rewire `iteration.py:455`/`operator.py:1938`; #207 cross-storage injection STAYS implicit (resolved) | NI type / BI arithmetic | 347 wall; named-composition test |
| **B.6** | tighten `TimedFullField` slots | NI | composite reject tests (update `match=`) |
| **C.1** | archivist: Sphinx theory page (View-G, storage×role×locus, dimensional table, named-composition, TraceSpace) + refresh `operator_algebra.rst` | — | `-W` clean; Nexus reload |
| **C.2** | issues: close #201; move #205; close #197; amend `wave_o_operator_typing.md`; file the 5-red curvilinear issue; `error_catalog.md` | — | audit clean |

**Verify:** `python -O -m pytest`; after each phase rebuild Sphinx + `python -m
tests._harness.audit`. Re-run the field/space set after each step; red set = 5.
