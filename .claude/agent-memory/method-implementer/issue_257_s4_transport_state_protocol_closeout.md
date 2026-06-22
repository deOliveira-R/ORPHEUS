# #257 S4 — `TransportState(Protocol)` (the honest generic) — closeout

**Branch** `feature/field-typed-operator-algebra` (HEAD `4ebade6`, NOT committed — main agent reviews + commits).
**Date** 2026-06-19. **Stage** S4 of the CoefficientField+operator-as-promotion campaign (#257).
**Nature** PURE-TYPING (bit-identical, zero runtime/behavioral change) + one new foundation test.

## What S4 is

Dissolve the opaque `Generic[V]` carrier into a NAMED, runtime-checkable `TransportState`
Protocol that REFINES the numerics `Vector`. `TransportState` is a NEW ROLE (the abstract
structural contract) — NOT the old class-name of `TimedFullField` (renamed 2026-05-28).
`TimedFullField` is the impl; `TransportState` is the contract it satisfies WITHOUT inheritance.

## Deliverables (all shipped, in the working tree)

1. **`orpheus/transport/state.py` (NEW)** — `@runtime_checkable class TransportState(Vector, Protocol)`.
   Members (read-only properties, see below): `bulk: BulkField`, `boundary: BoundaryField`,
   `history_depth: int`. Rich module+class docstring in the `vector.py` register (WHY a structural
   Protocol, the layering, that it refines `Vector` + the discriminating `np.ndarray`-is-NOT check).
2. **`orpheus/transport/__init__.py`** — export `TransportState` (import block + `__all__`).
3. **`orpheus/transport/multiplication_operator.py`** — re-pointed: class generic
   `LinearOperatorMixin["TransportState"]`; `apply`/`solve`/`apply_transpose` param+return
   `"TimedFullField"`→`"TransportState"`; `TransportState` added to the `TYPE_CHECKING` block.
   Bodies UNCHANGED (still construct concrete `TimedFullField` returns — covariant-valid).
4. **`tests/transport/test_transport_state.py` (NEW)** — `@foundation`, 4 tests, the discriminating
   type-check: `TimedFullField` IS-a `TransportState` AND `Vector`; `np.ndarray` IS-a `Vector` but
   NOT a `TransportState`; bare `AngularFlux` IS-a `Vector` but NOT a `TransportState`.
   `CollisionOperator` (`sn/operator.py:512`) needs NO change — it inherits the re-pointed annotations
   (no `apply`/`solve` override; confirmed by grep, the `def apply`@839 / `def solve`@888 belong to
   `InvertibleOperator`@662).

## THE READ-ONLY-PROPERTY DECISION (the one design judgment in D1)

First-draft `state.py` declared the three members as plain Protocol attributes
(`bulk: BulkField` etc.). Standalone pyright then RED'd 2× `reportReturnType` on
`multiplication_operator.py:209/234`: a writable Protocol attribute is INVARIANT, so the FROZEN
dataclass `TimedFullField` (read-only fields) is "not assignable" to a Protocol requiring writable
attributes — the `apply -> TransportState` body returning a `TimedFullField` failed. FIX (principled,
NOT a `# type: ignore`): declare the three members as read-only `@property` getters in the Protocol.
This is BOTH the correct contract (the operator algebra only ever READS `psi.bulk.values` /
`psi.history_depth`, never assigns) AND resolves the variance (a read-only property is satisfied by a
frozen dataclass field). Result: standalone pyright `state.py` + `multiplication_operator.py` = **0/0**.
`runtime_checkable` still discriminates correctly with property members (verified empirically, Py 3.14).

## THE +1 NET-NEW HUNT (do NOT repeat the wrong diagnosis)

Full pyright first read **2296 (+1 vs 2295 baseline)**. I mis-attributed it to the operator re-point
(D2). ISOLATION (manual HEAD-content revert via `git show` write-back — NOT `git checkout`, L28-safe —
+ moving the 2 new files aside, then re-adding D1/D2/D3 one at a time) proved:
- TRUE baseline (all S4 reverted) = **2295** ✓.
- D1 only (`state.py` + export) = **2295** (clean).
- D1 + D3 (test, no D2) = **2296** → the +1 IS the TEST.
- The error: `test_transport_state.py:124 reportGeneralTypeIssues` — pyright's "Class overlaps Vector
  unsafely" on a LITERAL `isinstance(np.ndarray(...), Vector)` against a runtime-checkable Protocol.
  This is pyright flagging EXACTLY the discriminating fact the test ASSERTS (ndarray structurally
  matches Vector), not a hazard.
FIX (principled, NOT `# type: ignore`): route the 6 assertion-site `isinstance` calls through a small
`_is_a(candidate: object, protocol: type) -> bool` helper. The `object`-typed boundary keeps the
runtime `isinstance` IDENTICAL but stops pyright running its structural-overlap analysis on a concrete
literal. Result: full pyright back to **2295 (0 net-new)**.

LESSON for future Protocol intrinsic-property tests: a literal `isinstance(<concrete>, <Protocol>)`
trips `reportGeneralTypeIssues` "overlaps unsafely" when the concrete type structurally matches — route
through an `object`-typed helper (the runtime check is the deliberate point; the static overlap warning
is noise). Pairs with the vv Mode-8 `_require`/`pytest.fail` (-O-firing) discipline.

## ⚠ THE DEFERRED CASCADE (the D4 assessment — IMPORTANT for S8)

The aggregate pyright is 0 net-new (2295), BUT this MASKS a real local **+2 at the SN composition
seam offset by −2 coincidental resolutions elsewhere** (the `TransportState` re-inference cleaned up 2
unrelated downstream errors). Measured directly: `sn/operator.py` goes **2 errors (HEAD) → 4 errors
(S4)**. The two genuine net-new (both rooted in `MultiplicationOperator.apply` now being
`(TransportState)->TransportState` flowing into the INVARIANT `OperatorSum["TimedFullField"]`):
- `operator.py:794 reportArgumentType` — `InvertibleOperator.__init__` → `super().__init__(streaming,
  diagonal)`: `diagonal: CollisionOperator` is now `LinearOperator["TransportState"]`, but
  `OperatorSum["TimedFullField"].__init__(a,b: LinearOperator["TimedFullField"])` is invariant.
- `operator.py:641 reportIncompatibleMethodOverride` — `CollisionOperator.__add__` returns
  `InvertibleOperator | OperatorSum[TransportState]`; the `InvertibleOperator` arm is
  `OperatorSum[TimedFullField]`, not assignable to the base mixin's `OperatorSum[TransportState]`.

PROVEN that NO in-scope additive fix closes these (both trials measured, both reverted):
- Re-stating `CollisionOperator(MultiplicationOperator, LinearOperatorMixin["TimedFullField"])` →
  pyright 2297 (WORSE; the double generic-mixin MRO conflicts).
- Re-pointing `InvertibleOperator(OperatorSum["TransportState"])` (no body widen) → 6 errors / 2297
  (WORSE; its `apply(psi: TimedFullField)`/`apply_transpose` bodies become contravariance-violating
  overrides at lines 840/868 — the body-widening cascade).
Both fixes require widening `InvertibleOperator.apply`/`solve` bodies, which thread
`_require_typed_composite(field: TimedFullField)` (operator.py:171, 4 call sites) and
`LossRepresentation.loss_action(psi: TimedFullField)` — **exactly the firm scope boundary / S8 work.**

**S8 SCOPING NUMBERS (the D4 count for the eventual SN-operator re-point):**
- `loss_representation.py`: **14** `loss_action`/`loss_action_transpose` def signatures
  (7×`loss_action`, 7×`loss_action_transpose`); **15** `psi/field/...: "TimedFullField"` param sites
  threaded through them.
- `operator.py`: **1** `_require_typed_composite(field: "TimedFullField")` helper (4 call sites:
  StreamingOperator.apply@410 / .apply_transpose@455 / InvertibleOperator.apply@864 /
  .apply_transpose@881); **12** `: "TimedFullField"` / `-> "TimedFullField"` / `"AngularFlux |
  TimedFullField"` param/return sites on the StreamingOperator+InvertibleOperator surface.
- So S8's SN-operator re-point widens ~**14 loss_action signatures + 15 param sites + 1 helper (4
  sites) + ~12 operator-surface sites**, and AT THAT POINT lines 641/794 resolve for free (the seam
  becomes internally `TransportState`-consistent).

RECOMMENDATION to the main agent: S4 meets the literal aggregate gate (2295 = baseline, 0 net-new) and
delivers all 3 + the discriminating test. The 641/794 seam errors are the S8-owned cascade, currently
offset in aggregate (fragile — a future change that stops resolving the 2 elsewhere would expose them).
Acceptable to ship as the S4→S8 bridge; track on #257 so S8 closes them. If a STRICTLY clean
`sn/operator.py` is required NOW, the only path is to defer D2's `MultiplicationOperator` re-point into
S8 (ship S4 as D1+D3 only) — D1+D3 alone is genuinely 2295 with no hidden offset.

## Verification (pasted in the return)

- pyright full: **2295 errors / 19 warnings = EXACT baseline, 0 net-new, 0 new `# type: ignore`**.
  Standalone `state.py` + `multiplication_operator.py`: **0/0**.
- Regression subset (-O): **7 failed / 1972 passed / 6 skipped / 5 xfailed** — the 7 are EXACTLY the
  baseline (#250 SPHERE ×5: 3×`test_vacuum_bulk_bit_identical_1d[*-SPH]` + 2×`test_sphere_{1,2}g_apply_
  bit_identical`; #232 mu_y ×2: `test_2d_mesh_resolution`, `test_two_d_cartesian_loss_action_returns_
  result`). ZERO non-baseline regression.
- New + re-pointed suite (-O): `test_transport_state.py` (4) + `test_multiplication_operator.py` (11)
  = **15 passed**.
- Runtime import OK; `TransportState` is `runtime_checkable`.

## NOT done (correctly out of scope)

- NO algebra-of-record manifest (Branch-1 SymPy / L1 xverif / Sphinx stub / archivist DISPATCH): S4 is
  a pure-typing stage — names a structural contract, builds no new reference, derives no equation.
  Verification = the intrinsic-property type-check (the foundation test). Same posture as S1/S2/S3a.
- NO new ERR (no bug found/caught; the read-only-property + `_is_a` fixes are typing-tension
  resolutions, not numerical bugs).
- StreamingOperator / InvertibleOperator / `_require_typed_composite` / `loss_action` UNCHANGED
  (firm scope boundary; the re-point folds into S8).
- NOT committed (main agent reviews + commits).

## Files touched
- NEW `orpheus/transport/state.py`
- `orpheus/transport/__init__.py` (export)
- `orpheus/transport/multiplication_operator.py` (generic + 3 annotation re-points + TYPE_CHECKING import)
- NEW `tests/transport/test_transport_state.py`
