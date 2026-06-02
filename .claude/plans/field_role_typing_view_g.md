# Typed Field Vocabulary + View-G Space Layer (#205 / #201; pre-#208)

**Branch:** `refactor/field-role-typing` (worktree child of `refactor/sn-operator-algebra` @ ad813fd).
**Approved plan (host copy):** `~/.claude/plans/glimmering-launching-lantern.md` — this is the durable in-repo copy.
**Issues:** #205 (cross-method field architecture — SN/Scalar portion), #201
(angular role split), #207 (space property-composition — applied), #197 (units
discipline — superseded), #208 (operator typing — **deferred**).
**Mode:** main-agent direct authorship, turn-by-turn user steering
(`feedback_no_method_implementer_for_surgical_carves`).

---

## SESSION STATE — PHASE A + B.1 + B.2 + B.3 + B.4 DONE, RESUME AT B.5 (2026-06-01)

> B.3 minted all six role leaves + renamed the geometry `BoundarySource`
> Protocol → `InflowSourceSpec`. B.4 added the `UNITS` (pint) class
> constant to every leaf (View-G, eV-free, `sr`-explicit/exact-compared)
> + units-aware diagnostics. **B.5 is NEXT** — the `from_balance` +
> dimensional-sin operator-output carve, which ALSO does the held-from-B.3
> boundary-residual matvec retype. B.5 trips the test-architect proactive
> trigger (dispatch the verification plan before implementing).

**B.1 DONE (`6e70ec1`):** 5 storage-base ABCs in `transport/fields/
_bases.py`; 6 leaves re-parented (−602 lines of duplicated machinery).
**B.2 DONE (`698a587`):** HARD rename (user chose no shim)
`PerOrdinateSource→AngularSource`, `IsotropicSource→ScalarSource` +
`git mv` modules + `_SPACE_NAME`s + `SNMesh.zeros_*` factories; 37 files
/ ~342 occurrences, 0 residual; cross-class injection dunder preserved.
Both bit-identical; verified via 760 (primitives+transport+operators) +
sentinel 36/36.

**B.3 (BULK) DONE:** new `transport/residuals/` package (parallel to
`sources/`) — `AngularResidual`("angular_residual", on AngularField) +
`ScalarResidual`("scalar_residual", on ScalarField), THIN leaves
(`_SPACE_NAME` only; NO bespoke factory — a residual is *produced* by
an operator balance in B.5, not built from thin air). New tests
`tests/transport/residuals/test_typed_residuals.py` (22, foundation):
construction/shape-validation, same-class add/sub/scalar/neg closed,
mesh-binding guard, frozen, 2-D smoke, and the **load-bearing
cross-class rejection** — `AngularResidual − AngularSource` RAISES
*even though shape AND units match* (the dimensional-sin gate), same
for scalar, residual±flux, cross-family, and `inner_product`. Also
fixed a stale `field.py` module-docstring field-list (duplicate
`AngularSource`, `BoundaryFaceFlux`→`BoundaryFlux`). Additive /
bit-identical: 782 fast-tier (+22) + sentinel 36/36.

**B.3 (BOUNDARY) DONE (user overrode the #208 deferral — "don't defer
known-future work when we have the context advantage").** Two parallel
explorer maps (2026-06-01) established the asymmetry:
- **`BoundaryResidual` has a REAL, already-computed consumer** — the SN
  matvec emits the affine-BC face defect (`γ₊ψ − bc_estimate`, residual
  of `γ₋ψ = R·G·γ₊ψ + q`) at `operator.py` (`:330/:570/:852`, output at
  `:1263`), TODAY mistyped as `BoundaryFlux`; GMRES drives it to zero
  via `TimedFullField.to_flat`. Wiring = **(retype)**.
- **`BoundarySource` has NO consumer** — every MMS case vanishes on Γ₋
  (vacuum q≡0); `solve_sn_fixed_source` has no inflow arg; the
  prescribed-inflow stack has zero end-to-end consumers. Minting it is
  **role-grid completion**.

Shipped this step (additive / bit-identical; the matvec WIRING is NOT
done here — see below):
1. **Rename** geometry Protocol `BoundarySource` → **`InflowSourceSpec`**
   (it is a lazy per-face *generator*, not a field — recipe→snapshot).
   ~31 sites across `geometry/boundary/{_source,_base,prescribed_inflow,
   __init__}.py`, `sn/angular_operator.py`, `test_boundary_trace_law.py`,
   2 theory `.rst`. `BoundarySourceNotOnIncomingTraceError` (ERR-047)
   NOT renamed (negative-lookahead guard). Impls `NoSource`/
   `ConstantInflowSource` keep names.
2. **Mint** empty `BoundaryField` leaves: `BoundarySource`
   (`transport/sources/boundary_source.py`) + `BoundaryResidual`
   (`transport/residuals/boundary_residual.py`) — role-organized, like
   the bulk leaves. `from_spec(spec, mesh)` recipe→snapshot bridge
   DEFERRED (no consumer; `feedback_unify_after_two_instances`).
3. **Tests** `tests/transport/fields/test_boundary_source_residual.py`
   (25, foundation): the load-bearing invariant is **cross-class raises
   despite all three boundary leaves sharing the SAME `mesh.trace`** —
   the space gate would PASS, the **class** gate rejects.
Verified: geometry 56 pass; fast tier 1082 pass (the only `tests/geometry`
reds are 2 PRE-EXISTING TP-lift `Test188WiringContracts`, confirmed via
`git stash` at f0522b4); sentinel 36/36.

**B.3 BOUNDARY WIRING → B.5.** Retyping the matvec output boundary
`BoundaryFlux→BoundaryResidual` is INSEPARABLE from the bulk
operator-output retype: `TimedFullField` composition requires every
leaf's output boundary to share a class, and the typed arithmetic at
`iteration.py:455` (`Fψ+Sψ+q_ext`) / `:511` (`Lψ−Sψ−Fψ`) decides all
operator-output types together. So `BoundaryResidual` is *minted and
ready* here; its WIRING lands in **B.5** (with `from_balance` + a
test-architect verification plan, per the operator-algebra-carve
proactive trigger).

**B.4 DONE** — `UNITS` (pint) on every role leaf, View-G units-on-field.
New `orpheus/numerics/units.py`: one shared `UREG` + FOUR named
signatures (`ANGULAR_FLUX_UNITS` `1/(cm²·s·sr)`, `SCALAR_FLUX_UNITS`
`1/(cm²·s)`, `ANGULAR_RATE_UNITS` `1/(cm³·s·sr)`, `SCALAR_RATE_UNITS`
`1/(cm³·s)`); 10 leaves map onto the 4. Locked conventions (user-steered,
all documented in `units.py`):
- **eV-free** — a stored flux is always energy-bin-integrated (MG group OR
  CE-MC tally bin), so eV cancels; continuous-energy lives in the XS
  data/kernel, not the field. Spans MC + deterministic (same signature ⇒
  CE-MC ↔ MG-deterministic V&V cross-check).
- **`sr` kept explicit, compared by EXACT unit equality** (pint treats
  `sr` as dimensionless ⇒ `.dimensionality` would hide a missing-angular-
  integration / missing-`/4π` bug, the ERR-039 class). 4 exact signatures
  collapse to 2 SI dim classes — documented. (`sr`-as-base-dimension was
  evaluated + REJECTED: fragile in pint 0.25.3.)
- **`UNITS` is metadata, NOT the gate** — the gate stays class identity.
  Declared as a bare `UNITS: ClassVar[Unit]` contract on `Field` (mirrors
  `_SPACE_NAME`; elegance-enforcer should-fix) so a forgetful leaf raises
  on `.UNITS` rather than silently feeding `None` to #208's unit-gain.
Diagnostics wired: `Field._check_partner` TypeError surfaces both operands'
units + the "same units ≠ meaning, use `from_balance`" guidance;
`Field.__repr__` is now concise + units-aware (`repr=False` on Field + 5
bases + 10 leaves ⇒ one inherited repr, no array dump). Doc fixes:
`scalar_flux` (eV→eV-free), `harmonic_moment_field` (WRONG `1/(cm²·s·sr·eV)`
→ `1/(cm²·s)`; the moment is angle-integrated, `ℓ=0` IS scalar flux —
verified against the no-prefactor `SphericalHarmonicSpace` convention),
all source/residual/boundary leaf "Units" blocks. New test
`tests/transport/test_field_units.py` (43). Verified 1411 fast +
sentinel 36/36; elegance-enforcer PASS (1 should-fix taken).
**DEFERRED to C.1 (Sphinx debt):** add `numerics.units` + the transport
leaves to the autodoc tree and resolve the new `:data:`/`:mod:` refs.

**B.5 IS NEXT** — `IterationResidual.from_balance` + the dimensional-sin
rewire at `iteration.py:455` (`Fψ+Sψ+q_ext`) / `:511` (`Lψ−Sψ−Fψ`) /
`operator.py:1938`, INCLUDING the boundary-residual matvec retype
(`BoundaryFlux`→`BoundaryResidual` at `operator.py:330/570/852/1263`,
held from B.3). This is the operator-output carve — trips the
test-architect proactive trigger (dispatch the verification plan FIRST).
`from_balance` can use the B.4 `UNITS` for a dimensional sanity check
(exact-`sr` comparison). Then **B.6 (`TimedFullField` slots →
`bulk: BulkField, boundary: BoundaryField`)**.

**As-built shape a fresh session needs (post-compaction):**
- Bases live in `orpheus/transport/fields/_bases.py`: `BulkField`
  (mesh + mesh-binding `_check_partner` + `ng/nx/ny` + abstract
  `_phase_space_shape`) → `AngularField`(N,ng,nx,ny)/`ScalarField`(ng,nx,ny)/
  `MomentField`(marker); `BoundaryField` (mesh + TraceSpace contract +
  `layout` property + `face_view` + `zeros_for_sn_mesh`/`from_face_arrays`
  via `cls`). Angular/Scalar `from_mesh` is parametrized by the leaf's
  `_SPACE_NAME` ClassVar.
- Leaves (all THIN — `_SPACE_NAME` + role-specific methods only):
  `AngularFlux`("angular_flux")/`AngularSource`("angular_source",
  `angular_source.py`) on AngularField; `ScalarFlux`("scalar_flux")/
  `ScalarSource`("scalar_source", `scalar_source.py`) on ScalarField;
  `HarmonicMomentField` on MomentField; `BoundaryFlux` on BoundaryField.
- **B.3 leaf pattern (as built):** copied the B.1/B.2 thin-leaf shape.
  `AngularResidual(AngularField)` + `ScalarResidual(ScalarField)` =
  `_SPACE_NAME` only (+ `UNITS` in B.4), in new `transport/residuals/`.
  `BoundarySource`/`BoundaryResidual` `(BoundaryField)` = DEFERRED to
  #208 (name collision with the geometry-layer `BoundarySource` Protocol;
  no current consumer — see SESSION STATE).
- **Verification workflow (load-bearing):** per-tier single-process
  (NEVER whole `tests/sn` — OOMs). `pytest -m sentinel` WITHOUT `-O`
  (vv Mode 8) = seconds-fast full-DAG net (36 nodes). Fast bit-identity =
  `tests/{transport,sn/primitives,sn/operators}` (~760, ~4 s). The full
  sweep/solve/eigenvalue convergence tiers are TOO SLOW to run whole
  (issue **#211** — unmarked-slow studies; even `-m "not slow"` drags on
  sweep/solve). Use `-n 6` (NOT `-n auto` — oversubscribes the 10-core box).
  Test env: `PYTHONPATH=<worktree> .venv/bin/python` (worktree shadows the
  editable install); host env.

---

### (historical — superseded by the SESSION STATE above)

#### PHASE A COMPLETE, RESUME AT B.1 (2026-06-01)

**Phase A (the space layer) is DONE: A.1 + A.2/A.3 + A.4 + A.5.** Between A.3
and A.4 the session took a large
**user-directed testing-infrastructure detour** (all test-only/config — ZERO
production changes beyond the A.2/A.3 + A.4 commits). Commit chain on
`refactor/field-role-typing` (newest last):

- `30c858c` — **A.2/A.3 TraceSpace unification** (the field-typing work).
- `0defb4d` — fix ERR-055: curvilinear `sig_t` layout drift in TEST fixtures
  (pre-existing base-branch bug; `_sweep_1d_unified` is production, tests were
  on the obsolete `(nx,ng,ny)` layout → migrated to the `(ng,nx,ny)` producer).
- `d95b64a` / `8e95db9` / `d0634c9` / `6052664` — test infra: **pytest-xdist
  pinned `<3.8`** (3.8.0 DEADLOCKS on Py3.14.3; 3.7.0 works), taxonomy plan,
  regenerated V&V matrix, runner guidance.
- `e851e76..842dd4a` — **SN test-suite TAXONOMY REORG** (see
  `sn_test_taxonomy.md`): 1378 tests → capability-tier dirs
  `tests/sn/{primitives, operators, sweep/{core,slab,curvilinear,cartesian_2d},
  solve, eigenvalue, verification/{mms,analytical}, regression}` + `cap()`
  markers. **⚠ TEST PATHS CHANGED** — operator tests now under
  `tests/sn/operators/`, sweep under `tests/sn/sweep/...`, etc.
- `d9cd29f..138b4b0` — **SENTINEL HARNESS** (see `sn_sentinel_harness.md` +
  `tests/_mutation/`): `pytest -m sentinel` (~4 s, 36 nodes, one per capability
  node), cosmic-ray mutation-validated (diamond.py 96.8 %).
- `9a5e09e`/`94d3aa4`/`887f2de`/`7851eb2`/`5e9f165` — regression tolerance
  redesign DONE (see `sn_regression_tolerance_redesign.md`): magic floors →
  principled `conv_tol`/`nulp` + `DriftWarning` tripwire + regenerated snapshots.
- **`f86846e` — A.4 DONE** (the field-typing work): all 8 inline `sign(Ω·n)`
  face masks in `operator.py` now read `TraceSpace.{inflow,outflow}_indices_
  for_face` (the 4 outflow in `_compute_LpC`/`_compute_decomposition` + the 4
  2-D inflow in `_apply_2d_cartesian`); new public `SNMesh.trace` accessor.
  Bit-identical (boolean-mask vs sorted fancy-index). Gated: sweep 524 /
  operators 392 / 2-D MMS 3 passed. A2D-1 hash pin updated.
- `e46c63c` — fix a DEAD L1 gate found by A.4's sweep
  (`..._mms_manufactured_source_couples_groups` raised IndexError since
  written; restored: axis order + `/sum_w` per-ordinate-density norm).
- `e581409` — mark the slow verification studies `@pytest.mark.slow` so
  `pytest tests/sn/verification -m "not slow"` = 44 passed in **21 s** (was
  >13 min). Deeper test-suite speed work (Green's-reference caching, reference
  substitution for the not-research-grade Nyström tests, `-n 6` default,
  bare-assert vv-Mode-8, the `test_cylinder_l1_sweep_vs_krylov_twin_path`
  XPASS(strict) stale-xfail) deferred to **issue #211** (user's dedicated
  verification-suite + harness session).

**NEW test workflow (post-reorg) — use this, not the stale gates below:**
- Run **PER TIER** single-process (light: operators 395 tests / 3 s / 186 MB).
  NEVER the whole `tests/sn` tree single-process (OOMs — 1378 items).
- `pytest -m sentinel` = the seconds-fast always-on gate — **run WITHOUT `-O`**
  (vv **Mode 8**: `-O` strips bare `assert` in non-rewritten helpers → false green).
- xdist 3.7 works: `-n auto` within a tier is fine; whole-suite `-n auto` is
  ~21 GB aggregate + slow (CI only).

**Red baseline for A.4 onward:** the 6 curvilinear-sweep IndexErrors are FIXED.
Remaining pre-existing reds = 5× Wave-T.1 TP-lift `isinstance(bc.inner, …)`
(`bc.inner` is now a `TensorProductOperator`) + 1× stale `mu_z` regex — all in
the moved `tests/sn/operators/` + realizer-wiring/bound_compat files; stale
assertions tracked for C.2 (NOT correctness bugs). Snapshot-drift reds are
being fixed by the regression agent.

**A.4 DONE (`f86846e`).** All 8 inline face masks now read the `TraceSpace`
selectors via the new `SNMesh.trace` accessor; bit-identical; gated green.

**A.5 DONE (`57c5151`).** `BoundaryFlux` re-homed onto `TraceSpace`: dropped the
field-side `layout` (→ read-through property `self.space.layout`), factories
source `space = mesh.trace`, `__post_init__` enforces "space must be a
TraceSpace" (illegal-states-unrepresentable), retired both
`FunctionSpace("sn_boundary_flat")` builds. `mesh` kept solely as the
cross-mesh guard (drop it once `TraceSpace` identity is mesh-distinguishing —
#208). Bit-identical: test_boundary_flux 36 / operators-batch 439 /
sweep+transport+solve 652. Migrated the 2 direct-ctor test sites + added a
negative guard test.

**(historical) Phase B overview** — Phase B is the FIELD
VOCABULARY: storage-base ABCs (`BulkField`/`AngularField`/`ScalarField`/
`MomentField`/`BoundaryField`) deduping the `AngularFlux`↔`PerOrdinateSource`
machinery; then renames (`PerOrdinateSource→AngularSource`,
`IsotropicSource→ScalarSource`), new role leaves (`Angular/ScalarResidual`,
`Boundary{Source,Residual}`), the `UNITS` class constant, and the
dimensional-sin rewire (`iteration.py:455`, `operator.py:1938`) via
`IterationResidual.from_balance`. See the step table + Phase B section below.
B.1–B.2 + B.6 are bit-identical refactors; B.3–B.5 add new role types/values
pinned by new tests.

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

**⚠ STALE (pre-reorg) — superseded by "SESSION STATE" above.** The old gate
was: SN operator-algebra-core 12-file gate = 347 pass; "pin pre-existing red at
exactly 5". POST-REORG those file groupings no longer exist (tests moved into
capability tiers) and the red set changed (curvilinear 6 fixed; 5 TP-lift + 1
`mu_z` remain). For A.4 onward use the per-tier workflow + `-m sentinel` gate
from SESSION STATE, and verify A.4 bit-identity against the relevant
`tests/sn/sweep/` + `tests/sn/operators/` tiers (the operator-algebra tests).

| Step | Scope | Bit-identity? | Gate |
|---|---|---|---|
| 0a/0b | test-architect spec + explorer L20 audit | — | done |
| **A.1 ✅** | remove `units` from `FunctionSpace` (`space.py`: field, `__eq__`, `__repr__`, `_tensor_product_units`); drop Layer-3 assert in `field.py`; migrate unit-pins (3 in `test_space_algebra.py` + **5 in `test_field.py`** — grep found more than the spec listed); fix latent cm²→cm³ error in field.py docstring; scalar_flux.py docstring | **BI** (units `None` everywhere) | **DONE: 237 passed `-O`** (space+space_algebra+field+operator+transport+typed_sources) |
| **A.2+A.3 ✅** (MERGED) | **TraceSpace unification** — DONE. Built one concrete `TraceSpace(FunctionSpace)` (whole-boundary `shape=(layout.total_size,)`, `layout: FaceLayout`, signed `omega_dot_n (n_faces,N)`, `_TANGENTIAL_EPS=4·machine_eps`, inflow/outflow selectors). Migrated SNMesh (`_inflow_trace`+`_outflow_trace`→one `_trace`), method_space (`inflow_trace`/`outflow_trace`→`trace`), `_resolve_one` (`for_face(trace=)`). Decoupled `IncomingSourceOperator` probe + `BoundarySource.evaluate` from the trace → **`evaluate(shape)`** (retired the fake-trace-as-shape-carrier anti-pattern). Deleted `Inflow/OutflowTraceSpace` + `boundary_trace_space()` + `numerics/trace_space.py` shim. **Curvilinear-one-boundary resolved (see below).** | **BI both sides** | DONE — verified in memory-safe targeted chunks (full `tests/sn` OOMs: 1376 single-process items, pre-existing infra, NOT ours). **Zero new failures.** trace/space/method_space/realizer_wiring/bound_compat/btl: 68 pass / 5 pre-existing TP-lift red. sn_boundary_realizer+native_matvec+streaming_decomposition(Res-A)+typed_sources: 80 pass / 1 pre-existing stale red (`test_white_z_axis…raises`, stale `mu_z` regex → msg is now "degenerate for this face"). unified_matvec_cylinder (curvilinear `bc_left=None`): all pass. phase_c_gates(L0)+solution: 42 pass / 6 expected-xfail. Masks bit-identical (xmin↔old "left", xmax↔old "right"). **Pre-existing base-branch reds confirmed via `git stash` (same failures with my changes removed):** 6× curvilinear-sweep `IndexError` in `sweep_cache.py:431` `sig_t[:, geom.chain_idx]` (test_spherical ×3 + test_cylindrical ×3 — the `_sweep_1d_unified`/CollisionCache legacy path, superseded by the passing matvec; a separate subsystem). |
| **A.4 ✅** | consolidated all 8 inline `sign(Ω·n)` face masks in `operator.py` to read `TraceSpace.{inflow,outflow}_indices_for_face` (4 outflow in `_compute_LpC`/`_compute_decomposition` + 4 2-D inflow in `_apply_2d_cartesian`); new public `SNMesh.trace`. `f86846e` | **BI** (boolean-mask ≡ sorted fancy-index) | **DONE: sweep 524 / operators 392 / 2-D MMS 3 passed; A2D-1 hash pin updated.** Also: restored a dead L1 gate (`e46c63c`) + marked slow studies (`e581409`, fast loop 21 s); deferred deep test-speed → issue #211 |
| **A.5 ✅** | re-homed `BoundaryFlux` onto `TraceSpace`: dropped field-side `layout` (→ read-through property `space.layout`), factories source `space=mesh.trace`, `__post_init__` enforces TraceSpace, retired both `sn_boundary_flat` builds; `mesh` kept as cross-mesh guard. `57c5151` | **BI** (only `space.name` flips, not values) | **DONE: test_boundary_flux 36 / operators-batch 439 / sweep+transport+solve 652 passed.** Migrated 2 direct-ctor test sites + added negative guard test |
| **B.1 ✅** | created all 5 storage-base ABCs in `transport/fields/_bases.py` (plan-literal): `BulkField`(mesh+mesh-binding+ng/nx/ny+abstract `_phase_space_shape`) → `AngularField`/`ScalarField`(parametrized `from_mesh` via `_SPACE_NAME`)/`MomentField`(marker); `BoundaryField`(TraceSpace contract+layout+face_view+factories via `cls`). Re-parented all 6 leaves (−602 lines). `6e70ec1` | **BI** (values + space names + public API preserved; `isinstance` holds) | **DONE: smoke (MRO/abstractness/`_SPACE_NAME`) + transport+primitives+operators 760 passed + sentinel gate 36/36 (full DAG).** |
| **B.2 ✅** | **HARD rename (no shim, user decision)**: `PerOrdinateSource→AngularSource`, `IsotropicSource→ScalarSource` (+ `git mv` modules, `_SPACE_NAME`s, `SNMesh.zeros_*` factories). 37 files / ~342 occurrences, 0 residual. Cross-class injection dunder preserved. `698a587` | **BI** (rename touches no values) | **DONE: 760 (primitives+transport+operators) + sentinel 36/36; imports clean; 0 residual tokens.** |
| **B.3 ✅ (all leaves + rename)** | bulk `Angular/ScalarResidual` (`transport/residuals/`) + boundary `BoundarySource` (`transport/sources/`) / `BoundaryResidual` (`transport/residuals/`); renamed geometry Protocol `BoundarySource`→`InflowSourceSpec`. `BoundaryField` now has its 2nd/3rd instances. **Boundary-residual matvec WIRING → B.5** (inseparable from bulk operator-output retype). | new (additive/BI) | 22 bulk + 25 boundary cross-class-raise tests; 1082 fast (−2 pre-existing TP-lift) + sentinel 36/36 |
| **B.4 ✅** | `UNITS` (pint) on all 10 leaves via `numerics/units.py` (4 named signatures; eV-free; `sr`-explicit/exact-compared; bare `ClassVar` contract on `Field`). Diagnostics: units in `_check_partner` error + concise `Field.__repr__`. Doc fixes (incl. `harmonic_moment_field` wrong-`sr`). NOT the gate (class identity is). | new (additive) | 43 `test_field_units.py`; 1411 fast + sentinel 36/36; elegance PASS |
| **B.5 ← NEXT** | `from_balance` + dimensional-sin rewire `iteration.py:455`/`:511`/`operator.py:1938` + **boundary-residual matvec retype** (`BoundaryFlux`→`BoundaryResidual`, held from B.3); #207 cross-storage injection STAYS implicit. **test-architect FIRST** (operator-output carve trigger). | NI type / BI arithmetic | named-composition + matvec-residual-type tests |
| **B.6** | tighten `TimedFullField` slots | NI | composite reject tests (update `match=`) |
| **C.1** | archivist: Sphinx theory page (View-G, storage×role×locus, dimensional table, named-composition, TraceSpace) + refresh `operator_algebra.rst` | — | `-W` clean; Nexus reload |
| **C.2** | issues: close #201; move #205; close #197; amend `wave_o_operator_typing.md`; file the 5-red curvilinear issue; `error_catalog.md` | — | audit clean |

**Verify:** `python -O -m pytest`; after each phase rebuild Sphinx + `python -m
tests._harness.audit`. Re-run the field/space set after each step; red set = 5.

---

## A.2/A.3 — TraceSpace unification: resolved design (2026-05-31)

Merged into ONE step — the `TraceSpace` name + the coupled face/eps
reconciliation make a clean build-then-migrate split unstable.

**Target type:** one concrete `TraceSpace(FunctionSpace)` carrying
`name="sn_trace"`, `shape=(layout.total_size,)`, `inner_product_weights=None`,
`layout: FaceLayout`, and `omega_dot_n: NDArray (n_faces, N)`. It is the
whole-boundary storage space (role-agnostic, View-G) and REPLACES both the
per-face `Inflow/OutflowTraceSpace` AND the ad-hoc
`FunctionSpace("sn_boundary_flat")` that `BoundaryFlux` currently builds.

- **Granularity:** whole-boundary `(total_size,)`; per-face access via the
  `FaceLayout` slot + the `(face, ordinate)` mask row. Selectors stay per-face
  (`inflow_indices_for_face(face)` / `outflow_indices_for_face`). The old
  per-face `(N,ng)` "space" becomes a derived view, NOT a class.
- **Directional leaf-data:** store the SIGNED `omega_dot_n` `(n_faces, N)`
  (`= sign_f · μ_axis(f)`); derive inflow (`< −eps`) / outflow (`> +eps`) /
  tangential (`|·| ≤ eps`) on demand. Single source of truth — A.4's
  `operator.py` masks (≈8 inline copies of `mu_x > ±eps`) read the SAME array.
- **eps — principled, NOT a magic number:**
  `_TANGENTIAL_EPS = 4.0 * np.finfo(np.float64).eps  # ≈ 8.9e-16` — a safety
  factor over the IEEE-754 dot-product round-off bound `d·u` (`d ≤ 3` spatial
  dims, `u = eps/2`) for the unit-vector projection `Ω·n = ⟨n, μ⟩`. EMPIRICAL
  (`eps_probe.py`, GL N=2..64 + Lebedev orders 3..53): nominally-tangential
  cosines are EXACTLY `0.0` (quadrature symmetry — odd-N central node + all
  off-axis 1-D components + Lebedev axis nodes); smallest GENUINE cosine
  `= 2.44e-2`; gap `[0, 0.024]` spans ~14 orders. eps sits 4× above the
  round-off floor and 2.7e13× below genuine → masks BIT-IDENTICAL to BOTH the
  operator's old `1e-15` AND the realizer's old `1e-12` (the `(eps, 1e-12)`
  band is empty). New foundation test: `_TANGENTIAL_EPS < min genuine |μ|`
  across GL+Lebedev (gap guard so a future quadrature can't silently violate it).
- **Face-naming reconciliation (fixes a latent curvilinear bug):** key the
  normal table on the LAYOUT's face names from `SNMesh.boundary_face_layout` —
  `xmin/xmax` (slab), `xmax` ONLY (1-D curvilinear: r=0 is a regularity
  condition, not a BC face), `xmin/xmax/ymin/ymax` (2-D Cartesian) — NOT
  trace_space's old `left/right` (which wrongly fabricated a mask for a
  non-existent curvilinear inner face). Restrict masks to the layout's faces.
- **Inner-product weights = None (Euclidean)** now — matches `BoundaryFlux`'s
  `sn_boundary_flat` → bit-identical. The physically-correct `|Ω·n|`-weighted
  boundary metric (partial currents / `BoundaryOperator` adjoints) is DEFERRED
  to Wave O; **recorded at #208, comment posted 2026-05-31**
  (`#issuecomment-4589013439`).
- **Consumers to migrate (explorer L20 audit):** `sn/geometry.py` (SNMesh: one
  `_trace` replacing `_inflow_trace`/`_outflow_trace`), `sn/method_space.py`
  (holds one trace; `inflow_indices_for_face` delegates), `sn/boundary_realizer.py`,
  `geometry/boundary/_source.py` (`evaluate(trace, face)` → per-face inflow
  values), `sn/angular_operator.py` (probe). Then DELETE the two subclasses +
  factory + shim. Migrate tests: `test_trace_space.py` (rewrite — rehome the
  per-face `sign(Ω·n)` correctness intent onto the unified space),
  `test_method_space.py`, `test_space.py:190-207` (the `boundary_trace_space`
  tests). A.5 then re-homes `BoundaryFlux` onto this `TraceSpace` (FaceLayout
  field→space).
- **Probe artifact:** `eps_probe.py` lives in the job tmp (not committed); its
  result is captured above. The gap-guard test is the durable form
  (`test_eps_sits_in_the_round_off_to_genuine_gap`, `@foundation`).

### Curvilinear-one-boundary — RESOLVED (2026-05-31, user-steered)

The face-naming reconciliation surfaced a deeper physics question the user
posed directly: *does a solid sphere have one boundary or two?* **Answer: ONE
(the outer radius r=R).** The axis r=0 is NOT a boundary — it is a
coordinate-singularity / symmetry condition `ψ(0,μ)=ψ(0,−μ)`, the r→0 limit of
the angular-redistribution term `(1−μ²)/r ∂ψ/∂μ`, handled by the **angular pole
closure** (`MorelMontryAngularSweep` / Carlson seed), NOT a BC. Litmus: a
boundary is where you can impose inflow (vacuum/albedo/source/reflection); at
r=0 there is no surface, no "outside", and vacuum/source are meaningless — the
μ→−μ map is geometry-forced, not chosen. Modeling it as a reflective boundary
conflates "a mirror you placed" with "the sphere is symmetric about its
centre".

**Decisive evidence it's SAFE (the sweep already agrees):**
`operator.py:429/679` set `bc_inner = sn_mesh.bc_left if curvature=="cartesian"
else None` — the curvilinear matvec **never reads `bc_left`**; `has_inner_face =
"xmin" in boundary.layout.faces` is False for curvilinear; `sweep.py` curvilinear
branch reads only `bc_right`/`xmax`; `boundary_face_layout` already says
curvilinear = `xmax` only. So the legacy reflective `bc_left` for curvilinear
was **vestigial** — fabricated by `_resolve_bcs`, read by nothing, asserted only
by two already-red (TP-lift) tests.

**What changed:** the unified trace keys on `boundary_face_layout` → curvilinear
trace = `xmax` only (no phantom inner-face mask). `_resolve_bcs` no longer
fabricates a curvilinear `bc_left`: **`sn.bc_left = sn.bc_xmin = None` for
curvilinear** (the mesh's `bc_left` declaration at r=0 is moot — the axis is the
pole closure's). Slab keeps both `xmin`/`xmax` (genuine boundaries; matvec reads
`bc_inner = bc_left`). All bit-identical to prior numerics (the sweep never used
the vestigial op). Curvilinear BC-resolution tests rewritten to the principled
one-boundary contract (they stay red ONLY on the separate Wave-T.1 TP-lift
`isinstance` issue, tracked for C.2).

**Face naming:** 1-D now uses `xmin`/`xmax` (was `left`/`right`); radial axis IS
the x-axis. Universal `_FACE_NORMALS = {xmin:(0,-1), xmax:(0,+1), ymin:(1,-1),
ymax:(1,+1)}` serves every mesh. `_resolve_one` reflective-axis logic (`"y" if
face in {ymin,ymax} else "x"`) unchanged-correct.

### Remaining A-phase work (unchanged scope)
- **A.4** — consolidate `operator.py` inline `mu_x>±eps` masks (≈8 copies) to
  read `TraceSpace.outflow_indices_for_face`. The selectors are ready; the masks
  are bit-identical (eps unified). NOTE: the operator builds masks as boolean
  (`mu_x > eps`); selectors return indices — both valid for fancy-indexing the
  face views. Verify Resolution-A `assert_array_equal` + 347-wall.
- **A.5** — re-home `BoundaryFlux` onto `TraceSpace` (its `sn_boundary_flat`
  space → the unified trace; FaceLayout field→space). The trace already carries
  `layout` + `shape=(total_size,)` + Euclidean weights = bit-identical storage.
