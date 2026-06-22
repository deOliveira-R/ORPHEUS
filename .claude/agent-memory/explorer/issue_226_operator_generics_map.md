---
name: issue-226-operator-generics-map
description: Map for making numerics/operator.py's LinearOperator algebra Generic[V] over the vector/field type (pyright #226). Interface surface, ndarray-vs-field operator buckets, mixed-sum check (NONE), consumer scope, field hierarchy + TypeVar bound, Generic/RegistryMixin runtime-safety.
metadata:
  type: project
---

# Issue #226 — `LinearOperator` algebra `Generic[V]` design map (read-only investigation)

> **TWO passes in this file.** §0–§7 below were the broad "make the whole algebra
> `Generic[V]`" design map (2026-06-18 on `main`). **§8 is the FOCUSED #226 "B4"
> cluster re-derivation (2026-06-22, HEAD `7aebab3`)** — the variance break +
> `_as_boundary` `block_role` mutation, with live pyright errors grounded. **Read
> §8 FIRST for the B4 carve; §0–§7 is the wider context.** Line numbers in §0–§7
> are 2026-06-18 and HAVE DRIFTED — re-derive via Nexus. The SHAPE (buckets,
> mixed-sum verdict) is the durable part.
>
> ⚠ **STALE-PREMISE CORRECTION (verified against the current tree):** §0/§5 below
> recommended `V = TypeVar("V")` UNBOUNDED with "a structural `Vector` protocol …
> is a larger design; not required." **That `Vector` protocol has SINCE BEEN
> MINTED and the bound APPLIED.** The current tree
> (`orpheus/numerics/vector.py:151`) is `V = TypeVar("V", bound=Vector)`, where
> `Vector` is a structural Protocol (`__add__`/`__sub__`/`__rmul__`/`__truediv__`)
> satisfied by `np.ndarray`, every `Field` leaf, and `TimedFullField`. So the
> §0/§5 "unbounded" headline is OBSOLETE — the bounded-Vector design landed, and
> it is exactly the structural bound §5 said "would be needed if a bound is
> wanted." Do not re-recommend unbounded V.

## 0. Headline verdict (for the design) — 2026-06-18, SUPERSEDED on the V-bound

- A single `Generic[V]` parameterization of `LinearOperator` / `LinearOperatorMixin` /
  `OperatorSum` (+ the composer family) is **CLEAN**. **No production `OperatorSum`
  ever mixes an ndarray-operator with a field-operator.** All real sums are homogeneous.
  (Still true at HEAD.)
- ~~`V` must be UNBOUNDED~~ — **SUPERSEDED.** The current tree binds
  `V = TypeVar("V", bound=Vector)` to the structural `Vector` Protocol — the very
  abstraction this section deferred. `np.ndarray`, every `Field` leaf, and
  `TimedFullField` all satisfy `Vector` structurally (no inheritance), so the
  bound EXCLUDES nothing while documenting the contract. The old reasoning that a
  `Field` bound would exclude `np.ndarray`/`TimedFullField` remains correct — the
  resolution was a `Vector` bound, not no bound.
- Making these classes `Generic[V]` is a **pure typing change** — `Generic`/`Protocol` are
  runtime-inert here. The one machinery to respect is `RegistryMixin.__init_subclass__`
  (`key=` kwarg), but **no production class combines `LinearOperatorMixin` with
  `RegistryMixin`** (§6), so the metaclass/`__init_subclass__` collision risk is moot.

## 1. The operator interface surface (`orpheus/numerics/operator.py`)

### `LinearOperator` (runtime_checkable `Protocol`), line 257
- `capabilities: frozenset[str]` (attr, line 294) — single source of truth for apply/solve/transpose.
- `domain` property `-> Optional["FunctionSpace"]` (line 296) — **independent of V** (it's a
  `FunctionSpace`, not the vector). CONFIRMED.
- `codomain` property `-> Optional["FunctionSpace"]` (line 309) — same; **independent of V**.
- `apply(self, x: np.ndarray) -> np.ndarray` (line 317) — **the override-incompatibility source.**
  Under `Generic[V]` this becomes `apply(self, x: V) -> V`.
- NO `solve` / `apply_transpose` on the Protocol itself (deliberate — advertised via `capabilities`,
  docstring lines 259-265). So pyright never checked those against the Protocol; the 7
  `reportIncompatibleMethodOverride` are all on `apply`.

### `LinearOperatorMixin`, line 343 (the algebra dunders; defines no state)
- `capabilities: frozenset[str]` (line 377), `block_role: ClassVar[Optional[BlockRole]] = None` (397)
  — both **independent of V**.
- `domain`/`codomain` properties default `-> None` (404, 408) — **independent of V**.
- Dunders, all returning composer types (these need V-threading so `A: LinearOperator[F]`
  composes to `OperatorSum[F]`): `__add__`/`__radd__`/`__sub__`/`__rsub__` `-> "OperatorSum"`
  (415-425); `__mul__`/`__rmul__`/`__neg__`/`__truediv__` `-> "ScaledOperator"` (427-465);
  `__matmul__` `-> "OperatorProduct"` (467); `__and__`/`__rand__` `-> "TensorProductOperator"`
  (470-485); `__call__` aliases `apply` (487); `__pow__ -> "LinearOperator"` (500);
  `adjoint() -> "LinearOperator"` (531) and `H` property (555) → `_AdjointOperator`.
  NOTE: every dunder already carries `# type: ignore[arg-type]` on the composer construct —
  the Generic change is what would let those ignores be removed.
- `apply` / `solve` / `apply_transpose` are NOT defined here (the mixin relies on the
  concrete class to supply them) — so the mixin itself does not currently pin an `apply`
  signature; the incompatibility is concrete-class-vs-Protocol.

### `OperatorSum(LinearOperatorMixin)`, line 660
- `__init__(self, a: LinearOperator, b: LinearOperator)` (683) — would become `a: LinearOperator[V], b: LinearOperator[V]`.
- `apply(self, x: np.ndarray) -> np.ndarray` (739) → `(x: V) -> V`. Body `self.a.apply(x) + self.b.apply(x)`.
- `apply_transpose(self, x: np.ndarray) -> np.ndarray` (742) → `(x: V) -> V`.
- `domain`/`codomain` properties (729, 734) — independent of V.
- NO `solve` (sum is not invertible by closure law).

### Composer family (all `LinearOperatorMixin`, all `apply: np.ndarray -> np.ndarray` today)
`OperatorProduct` (746; +`solve` 810, +`apply_transpose` 814), `ScaledOperator` (819;
`apply(self, x, *extra, **kwextra)` 866 — note the variadic, +`solve` 869, +`apply_transpose` 873),
`IdentityOperator` (877; apply/solve/apply_transpose), `ZeroOperator` (899; `apply(self, x)`
untyped 952 — codomain_zero hook), `_AdjointOperator` (580; `apply(self, y: np.ndarray)` 623).
`OperatorProduct.domain/codomain` (797/802) and the others are all V-independent.

## 2. Operator classes by REAL apply input/output type

### BUCKET A — NDARRAY-OPERATORS (`apply: np.ndarray -> np.ndarray`)
These should be parameterized `LinearOperator[np.ndarray]` (or just keep concrete ndarray apply).
- `_AdjointOperator` — operator.py:580 (generic wrapper; V = inner's V in practice)
- `PermutationOperator` — operator.py:961
- `IncomingOrdinateMaskTensor` — operator.py:1065
- `PeriodicWrapOperator` — operator.py:1155
- `TensorProductOperator` — operator.py:1203
- `SumOfTensorProductsOperator` — operator.py:1324
- `DiagonalOperator` — operator.py:1417
- `RankOneOperator` — operator.py:1535
- `ProjectionOperator` (ABC) / `GalerkinProjection` / `PetrovGalerkinProjection` /
  `MomentProjection` / `ReconstructionOperator` (ABC) / `HarmonicMomentReconstruction`
  — numerics/projection.py:108/143/184/241/208/486. apply `np.ndarray -> np.ndarray`
  (projection.py:139, 231, 424, 589).
- `_BoundBoundaryOperator` — geometry/boundary/_bound_compat.py:72. `apply(self, psi: np.ndarray) -> np.ndarray` (122).
- `AngularAverageOperator` — sn/angular_operator.py:37. `apply(self, psi: np.ndarray) -> np.ndarray` (196).
- `IncomingSourceOperator` — sn/angular_operator.py:230. `apply(self, psi_out: np.ndarray) -> np.ndarray` (275).
- `IdentityOperator`, `ZeroOperator` — endomorphic/polymorphic; play in BOTH buckets
  (ZeroOperator's `codomain_zero` hook explicitly lets it emit a typed field zero —
  operator.py:932-955). Best left V-generic or `[np.ndarray]` with the field zero handled by hook.

### BUCKET B — FIELD-OPERATORS (apply consumes/produces a typed field)
- `StreamingOperator` — sn/operator.py:212. `apply(psi: TimedFullField) -> TimedFullField` (358);
  `apply_transpose` (426). → `LinearOperator[TimedFullField]`.
- `CollisionOperator` — sn/operator.py:511. `apply(psi: TimedFullField) -> TimedFullField` (592);
  `solve(q: TimedFullField) -> TimedFullField` (620); `apply_transpose` (651). → `[TimedFullField]`.
- `InvertibleOperator(OperatorSum)` — sn/operator.py:681. OVERRIDES `apply`/`apply_transpose`
  (858/886) `-> TimedFullField`; `solve(rhs: AngularFlux | TimedFullField) -> AngularFlux | TimedFullField`
  (907 — note the union, legacy AngularFlux tolerated in signature but body RAISES on non-TimedFullField
  rhs, operator.py:1064); `solve_moments(rhs: TimedFullField) -> TimedFullField` (975). → `[TimedFullField]`.
- `SNBoundaryOperator` — sn/boundary_operator.py:72. `apply(psi: TimedFullField) -> TimedFullField`
  (264); `apply_transpose` (296). → `[TimedFullField]`.
- `ScatteringOperator` — sn/scattering.py:297. **`apply` is a `singledispatchmethod`** (962)
  dispatching on `AngularFlux | TimedFullField | ScalarFlux | np.ndarray`; the TimedFullField
  arm `-> TimedFullField`. MIXED-INPUT by design (see §3 caveat). → effectively `[TimedFullField]`
  in the within-group path, but the singledispatch makes a single static `V` imprecise.
- `LegendreMomentScattering` — sn/scattering.py:174. `apply(moments: np.ndarray | HarmonicMomentField)
  -> np.ndarray | HarmonicMomentField` (250). Union-typed — neither pure bucket.
- `FissionOperator` — sn/fission.py:103. **`apply` is `singledispatchmethod`** (223) over
  `AngularFlux | TimedFullField | ScalarFlux | np.ndarray`; TimedFullField arm `-> TimedFullField`
  (283), ScalarFlux `-> ScalarSourceSink` (328), ndarray `-> np.ndarray` (355). MIXED-INPUT.

> The singledispatch operators (`ScatteringOperator`, `FissionOperator`) and the union-typed
> `LegendreMomentScattering` are the operators where `LinearOperator[V]` is an imperfect fit:
> their *actual* contract is "apply over a small union", not "apply over one V". Pyright will
> still accept `LinearOperator[TimedFullField]` for the call sites that pass TimedFullField,
> because the singledispatch base method's annotation can be set to `(psi: V) -> V`. The other
> arms are reached only by tests / non-SN consumers. RECOMMENDATION: annotate the singledispatch
> *base* method to the V used in the within-group algebra (TimedFullField) and let the registered
> arms keep their concrete annotations; do NOT try to make one TypeVar span the whole union.

## 3. Mixed-sum check (CRITICAL) — VERDICT: NO MIXED SUMS in production

The within-group algebra does NOT build `(L + C − S − F − B)` as one `OperatorSum`. It is a
**resolvent + variadic gains** split:
- `_within_group_triple(solver)` (sn/solver.py:175, returns `(L+C, S, B)`):
  `L = StreamingOperator`, `C = CollisionOperator`, sum `L + C` → `InvertibleOperator`
  (an `OperatorSum` of two TimedFullField operators — HOMOGENEOUS). `S = ScatteringOperator`,
  `B = SNBoundaryOperator` are returned as **separate `*gains`**, NOT added into the sum
  (sn/solver.py:216-221).
- The matvec is computed at the FIELD level, not via `OperatorSum`:
  - Krylov: `out = self.L.apply(psi); for g in gains: out = out - g.apply(psi)` then ravel
    (numerics/iteration.py:743-753, `A_matvec`). The `-` is `TimedFullField.__sub__`, not
    `OperatorSum`. So `L+C` (field OperatorSum), `S`, `B` never share an `OperatorSum`.
  - SI: `rhs = q_ext; for g in gains: rhs = rhs + g.apply(psi); psi = L.solve(rhs)`
    (numerics/iteration.py:503-514). Again field-level `+`.
- The ONE place `(L + C) - S - B` is built as operator algebra: `evaluate_residual`
  (sn/solver.py:225, docstring line 250-251 says "compose it from `_within_group_triple` as
  `(L+C) - S - B`"). This WOULD construct an `OperatorSum` of TimedFullField operators —
  **but it is still HOMOGENEOUS** (L+C, S, B all consume/produce TimedFullField; S and B's
  TimedFullField apply arms are the relevant ones). It is a diagnostic / future-DSA path,
  not in the convergence loop (docstring line 244-245). No ndarray operator is summed in.

Searched: `OperatorSum(`, `__add__`/`__sub__`/`__matmul__` on operators, and every
`LinearOperatorMixin`/`OperatorSum` subclass construction. The ndarray-operators
(DiagonalOperator, PermutationOperator, RankOneOperator, projection, TensorProduct, etc.) compose
among THEMSELVES (e.g. `RankOneOperator & IdentityOperator`, `SumOfTensorProducts` of
`TensorProductOperator`) — all ndarray. The field operators compose among themselves (the L/C
sum). **No cross-bucket `OperatorSum`/`OperatorProduct`/`&` exists in production.**

Single-Generic-parameter is therefore clean: each composer instance is `OperatorSum[np.ndarray]`
OR `OperatorSum[TimedFullField]`, never both.

## 4. Consumer scope (annotation sites of `LinearOperator`/`LinearOperatorMixin`/`OperatorSum`)

All annotation-position uses (grep + the Protocol's incoming type_uses):
- **(c) generic over both V** — the numerics iteration/eigenvalue layer. These are the sites
  that benefit MOST from `Generic[V]` (they're written shape-agnostic and run on BOTH L0 flat
  ndarray and L1 SN TimedFullField via the `to_flat`/`from_flat` ravellable protocol,
  iteration.py:178-214):
  - `SourceIteration.__init__(L: LinearOperator, *gains: LinearOperator)` — iteration.py:391-394
  - `KrylovAcceleration.__init__(L, *gains)` — iteration.py:656-659
  - `KEigenvalue.__init__(L, S, F: LinearOperator)` — iteration.py:934-936
  - `_default_keff_estimator(L, S, F: LinearOperator)` / the other estimator — iteration.py:251-253, 272-274
  - `as_scipy_linop(op: LinearOperator, ...)` — operator.py:1661 (flattens V→ndarray for scipy)
  These would each take a `LinearOperator[V]` and propagate V to the rhs/`q_ext`. Their `q_ext`
  params are currently annotated `np.ndarray` (e.g. iteration.py:450, 697) — a docstring-true
  but type-loose annotation; the real arg is the typed field. The Generic fix is the lever to
  retype those `np.ndarray` params to `V` and clear the cascade of `reportArgumentType` at the
  SN call sites that pass TimedFullField.
- **(b) field-operator consumers** — `evaluate_residual(loss_op, ...)` doc-annotated
  `LinearOperator` (sn/solver.py:249); `SNSolver._within_group_*` builders (sn/solver.py).
- **(a) ndarray-operator consumers** — the boundary-realizer chain returns/consumes
  `LinearOperator` as the realized BC op type: `_as_boundary(op: LinearOperator) -> LinearOperator`
  (sn/boundary_realizer.py:110), `SNBoundaryRealizer.realize(...) -> LinearOperator` (150),
  `realize_recursively(...) -> LinearOperator` (sn/boundary_realize.py:128),
  `BoundaryTraceLaw.realize(...) -> "LinearOperator"` (geometry/boundary/_base.py:274),
  `BoundaryRealizer.realize -> "LinearOperator"` (geometry/boundary/_realizer.py:122),
  `_BoundBoundaryOperator(inner: "LinearOperator")` (_bound_compat.py:102). These realize to
  ndarray-apply BC operators (`_BoundBoundaryOperator.apply: np.ndarray`). They would be
  `LinearOperator[np.ndarray]`.

Scope estimate: retyping the ~5 numerics consumers (group c) to `Generic[V]` + retyping their
`q_ext: np.ndarray` params to `V` is what clears the bulk of the `reportArgumentType` cascade,
because those are the sites SN passes `TimedFullField` into. The 7 `reportIncompatibleMethodOverride`
clear the moment the Protocol's `apply(x: V) -> V` admits the SN operators' `TimedFullField` apply.

## 5. Field-type hierarchy + TypeVar bound

MRO (from `orpheus/transport/fields/_bases.py` + `orpheus/numerics/field.py`):
- `Field(ABC)` — numerics/field.py:143. Root of the typed-field algebra (`values`, `space`,
  same-class/same-space closed arithmetic; frozen dataclass; `_check_partner` class-identity gate).
- `BulkField(Field)` (_bases.py:97) → `AngularField`(268) / `ScalarField`(372) / `MomentField`(461).
- `BoundaryField(Field)` (_bases.py:481).
- Concrete leaves: `AngularFlux(FluxRole, AngularField)` (fields/angular_flux.py:63);
  `ScalarFlux(FluxRole, ScalarField)` (scalar_flux.py:107);
  `BoundaryFlux(FluxRole, BoundaryField)` (boundary_flux.py:85);
  `HarmonicMomentField(FluxRole, MomentField)` (harmonic_moment_field.py:113).
  Plus the role grid: `AngularSourceSink(AngularField)`, `BoundaryResidual(BoundaryField)`,
  `AngularDisplacement(Displacement, AngularField)`, etc. — all descend from `Field`.
- **`TimedFullField` (transport/timed_full_field.py:123) is NOT a `Field` subclass** — it's a
  structured composite (`bulk: BulkField`, `boundary: BoundaryField`, `_history`) with
  delegate-style dunders (line 148). Its `DirectSumSpace`-Field backing is deferred ("Phase 3.6").

### TypeVar bound implication
- The actual carriers that flow through `apply` as V: `np.ndarray` (BUCKET A + L0 tests),
  `TimedFullField` (the SN within-group V), and — in the singledispatch/union arms —
  `AngularFlux`, `ScalarFlux`, `HarmonicMomentField` (all `Field`), plus `ScalarSourceSink`.
- `np.ndarray` is NOT a `Field`; `TimedFullField` is NOT a `Field`. A `TypeVar("V", bound=Field)`
  would **exclude the two most important Vs**. → use **`V = TypeVar("V")` UNBOUNDED.**
- If a bound is desired for documentation/safety, it must be a structural `Vector`/`SupportsVectorAlgebra`
  Protocol that BOTH `np.ndarray` and the field/composite types satisfy (defines `__add__`,
  `__sub__`, scalar `__mul__`). That is a larger design; not required to clear #226. Unbounded V
  is sufficient and correct.

## 6. Runtime-safety note (Generic / metaclass / RegistryMixin)

- `Generic[V]` and `Protocol[V]` are **runtime-inert** for these classes: no class here depends
  on `__class_getitem__` behavior, no `__orig_bases__` introspection, no metaclass conflict —
  EXCEPT the two below.
- **`_BlockRoleMeta` metaclass** (operator.py:184): used by the role MARKERS `BulkOperator`/
  `FullOperator`/`BoundaryOperator` (208/214/220), which are NOT `LinearOperatorMixin` subclasses
  and are never instantiated. `LinearOperatorMixin` itself is a PLAIN class (no metaclass). Making
  `LinearOperatorMixin(Generic[V])` does not touch `_BlockRoleMeta`. **No conflict.** (`Generic`'s
  metaclass is compatible with plain classes.)
- **`RegistryMixin.__init_subclass__`** (numerics/registry.py:87) consumes a `key=` class kwarg and
  calls `super().__init_subclass__(**kwargs)` (line 93). `Generic`/`Protocol` also define
  `__init_subclass__`/`__class_getitem__`. IF any class were `class X(LinearOperatorMixin[V], RegistryMixin)`,
  the `Generic` subscription + the `key=` kwarg + cooperative `super().__init_subclass__` would need
  MRO care (Generic must see `**kwargs` it doesn't recognize passed through — it tolerates this).
  **BUT no production class combines them.** The `class BoundaryOperator(LinearOperatorMixin,
  RegistryMixin, ABC)` at registry.py:31 is a **DOCSTRING EXAMPLE ONLY** (inside the module
  docstring, not real code). The real registry roots — `BoundaryTraceLaw` (geometry/boundary/_base.py:75),
  `DiscretizationSchemeBase` (sn/spatial/scheme.py:567), `PoleAngularClosureBase`
  (sn/spatial/pole_angular_closure.py:295), `PsiHalfAngleSeedBase` (sn/spatial/psi_half_angle_seed.py:298)
  — inherit `RegistryMixin, ABC` but **NOT `LinearOperatorMixin`**. So the Generic-vs-RegistryMixin
  `__init_subclass__` interaction **does not arise**. Safe.
- `@dataclass` operators (`StreamingOperator`, `CollisionOperator`, `ScatteringOperator`,
  `FissionOperator` carry `@dataclass` + class-attr `block_role`/`sigma_t`): the docstring at
  operator.py:391-396 already warns that `block_role` must be a PLAIN unannotated class attr (not
  `ClassVar[...]`) under `from __future__ import annotations` to avoid the dataclass field-detection
  bug. Adding `Generic[V]` as a base of these dataclasses is fine (dataclasses support Generic
  bases), but keep that plain-attr discipline; do not introduce new annotated class attrs on the
  Generic mixin that the dataclass machinery could misread as fields. `LinearOperatorMixin` is NOT
  a dataclass, so its own annotations (`capabilities`, `block_role`) are documentary only (operator.py:396).

## 7. Suggested parameterization (not applied — design only)

- `V = TypeVar("V")` (unbounded).
- `class LinearOperator(Protocol[V]): def apply(self, x: V) -> V: ...` (+ V-independent
  `capabilities`/`domain`/`codomain`).
- `class LinearOperatorMixin(Generic[V])` with dunders returning `OperatorSum[V]`/`ScaledOperator[V]`/
  `OperatorProduct[V]`/`TensorProductOperator[V]` — the `# type: ignore[arg-type]` on each dunder
  (operator.py:416-485) become removable.
- `class OperatorSum(LinearOperatorMixin[V])`, `OperatorProduct`, `ScaledOperator`, `IdentityOperator`,
  `ZeroOperator`, `_AdjointOperator`, etc. all `[V]`.
- SN operators: `class StreamingOperator(LinearOperatorMixin[TimedFullField])`, `CollisionOperator[...]`,
  `SNBoundaryOperator[...]`; `InvertibleOperator(OperatorSum[TimedFullField])`.
- ndarray operators: `[np.ndarray]` (or leave concrete ndarray apply; Protocol structural match still holds).
- numerics consumers: `SourceIteration(Generic[V])` / `KrylovAcceleration(Generic[V])` /
  `KEigenvalue(Generic[V])` taking `LinearOperator[V]` and `q_ext: V`.
- singledispatch operators (`ScatteringOperator`, `FissionOperator`): annotate the base
  `@singledispatchmethod` `apply` to the within-group V (TimedFullField); registered arms keep
  concrete annotations. Do NOT force one TypeVar over the whole union.

---

# §8. FOCUSED #226 "B4" cluster — variance break + block_role mutation (2026-06-22, HEAD `7aebab3`)

> The narrower, CURRENT-tree re-derivation for the B4 type-only carve. Line
> numbers here are HEAD `7aebab3`; re-derive via Nexus if drifted. **No runtime
> change is required for any B4 fix (flagged in §8.6).**

## 8.0 The premise is LIVE (grounded — `npx pyright orpheus/sn/boundary_realizer.py`)

The B4 cluster reproduces at HEAD:
- `boundary_realizer.py:125:8` — `Cannot assign to attribute "block_role" for
  class "LinearOperator[Unknown]"` (the `_as_boundary` mutation site).
- `boundary_realizer.py:176, 206, 207, 223, 224, 240, …` — `Argument of type
  "TensorProductOperator" cannot be assigned to parameter "op" of type
  "LinearOperator[Unknown]"` … `"LinearOperator[V@…]" is not assignable to
  "TensorProductOperator"`.

Both are ONE root cause: `LinearOperator(Protocol[V])` is **invariant in V** (V in
both `apply`'s param AND return), and `TensorProductOperator(LinearOperatorMixin)`
is **bare** (no `[V]`), so its inherited typevar materializes as `Unknown` while
its own `apply` is `np.ndarray`-typed → cannot unify under invariance.

## 8.1 Which leaves are bare-`LinearOperatorMixin` (→ Unknown) vs parametrized — CURRENT

`V = TypeVar("V", bound=Vector)` (`vector.py:151`). Parametrized today:
`_AdjointOperator`/`OperatorSum`/`OperatorProduct`/`ScaledOperator`/`IdentityOperator`/
`ZeroOperator` (all `LinearOperatorMixin[V]`), `StreamingOperator` /
`InvertibleOperator` (`["FullField"]` — a forward-ref STRING proxy for the real
runtime carrier `TimedFullField`; layering forbids the import).

**BARE (→ Unknown), the B4 surface** — these need `[np.ndarray]` (their `apply` is
already `np.ndarray`): `PermutationOperator` (~972), `IncomingOrdinateMaskTensor`
(~1076), `PeriodicWrapOperator` (~1166), `TensorProductOperator` (~1214 — the
direct B4 trigger), `SumOfTensorProductsOperator` (~1335), `DiagonalOperator`
(~1428), `RankOneOperator` (~1651), `MultiplicationOperator` (transport),
`CollisionOperator` (via `MultiplicationOperator`), `ScatteringOperator`,
`LegendreMomentScattering`, `FissionOperator`, `SNBoundaryOperator`,
`AngularAverageOperator`, `IncomingSourceOperator`, `_BoundBoundaryOperator`.

## 8.2 `block_role` — declaration / assignment / read (the :125 site)

- **Declared:** `LinearOperatorMixin.block_role: ClassVar[Optional[BlockRole]] = None`
  (operator.py:398) — `ClassVar` annotation is DOCUMENTARY (mixin is not a dataclass);
  leaves override with a PLAIN unannotated `block_role = BlockRole.X` (the :392
  comment: avoid the dataclass-field mis-detection under `from __future__ import
  annotations`).
- **Class-attr on leaves:** `Streaming`=FULL (op.py:325), `Scattering`=BULK
  (scattering.py:380), `Fission`=BULK (fission.py:179), `SNBoundary`=BOUNDARY
  (boundary_operator.py:115), `MultiplicationOperator`=BULK (→ `Collision`).
- **Per-instance assign in composers:** `_AdjointOperator` :623 (inner's role),
  `OperatorSum` :736 (`_join_block_roles` — DERIVED, the "role by construction"
  heart; `InvertibleOperator` relies on it), `ScaledOperator` :867.
- **POST-CONSTRUCTION MUTATION (the :125 error):** `_as_boundary(op)` does
  `op.block_role = BlockRole.BOUNDARY` (boundary_realizer.py:125) — stamps the
  realized BC (a generic `TensorProductOperator`/`ScaledOperator`) at the single
  producer site. ONLY post-construction mutation of `block_role` in the codebase.
- **`@property` forward:** `_BoundBoundaryOperator.block_role` (_bound_compat.py:113).
- **Reads:** load-bearing reader is `_BlockRoleMeta.__instancecheck__`
  (operator.py:200) → instance-level `getattr`. BUT tests read TYPE-LEVEL:
  `ScatteringOperator.block_role`, `FissionOperator.block_role`
  (test_operator_block_role.py:94/98). So `block_role` is GENUINELY both
  class-readable AND per-instance-mutated → **cannot collapse to pure ClassVar
  (mutation forbids) NOR pure instance field (class read forbids).** The mixin's
  current `ClassVar`-default + plain-attr-override shape is CORRECT; the :125 error
  is NOT a `block_role` modeling bug — `op` is typed `LinearOperator[Unknown]`
  (the Protocol, which has no `block_role` member), so the attribute is "unknown"
  only because the class is `Unknown`. **Fixing the variance fixes :125.**

## 8.3 `capabilities` — three override kinds, ONE mutation

Mixin declares `capabilities: frozenset[str]` (op.py:378) — plain instance attr.
Override kinds: (1) class-level `frozenset({...})` constant
(Identity/Mask/PeriodicWrap/RankOne/Permutation/projection); (2) dataclass
`field(default_factory=…)` (Streaming/Scattering/Fission/LegendreMoment/
Multiplication); (3) `@property` DERIVED (`SNBoundaryOperator.capabilities`
boundary_operator.py:137 — computed from per-face laws, NOT constant;
`_BoundBoundaryOperator` forwards; `_GaussSeidelResolvent`/`_MomentWindowedResolvent`
in solver.py:367/551). Composers set it in `__init__`.
**ONE mutation:** `InvertibleOperator.__init__` (sn/operator.py:766) rebinds
`self.capabilities = self.capabilities | frozenset({CAP_SOLVE})` after
`super().__init__` (the frozenset is immutable; the ATTRIBUTE is rebound).
**Verdict:** mixin MUST keep plain instance-attr `frozenset[str]` (a `ClassVar`
would make the composer assignments type-errors; `SNBoundaryOperator` NEEDS the
property; `InvertibleOperator` NEEDS the rebind). **B4 does NOT require touching
`capabilities` typing.**

## 8.4 `.solve` / capability-gating — how `iteration.py` reaches it

`solve` defined by: `Identity` :903, `OperatorProduct` :821 (iff both), `Scaled`
:880, `TensorProduct` :1323, `Diagonal` (all-nonzero), `InvertibleOperator` (the
WDD sweep), `_GaussSeidelResolvent`/`_MomentWindowedResolvent` (solver.py:377/559).
S/F/B/bare-L do NOT (rank-deficient). `SourceIteration.__init__(L: LinearOperator[V],
*gains: LinearOperator[V])` (iteration.py:413) gates at construction via
`_has(L, CAP_SOLVE)` → raises `MissingCapability(TypeError)` (iteration.py:435–443);
the call is `self.L.solve(rhs[, initial_guess=…])` (iteration.py:467/468), dispatched
on a runtime `inspect.signature` probe (:458).
**Type-system gap:** `LinearOperator(Protocol[V])` does NOT declare `solve` at all
(op.py:256 docstring: "advertised through the capability set rather than the type
system"). So `self.L.solve(...)` is an attribute access pyright cannot verify on a
`LinearOperator[V]`-typed value. The capability lattice is a deliberate RUNTIME
value-partition invisible to the static type — this is the §8.6(c) gap a
`SolvableOperator` sub-Protocol would close (but at large blast radius — defer).

## 8.5 The exact variance break

`apply(self, x: V, /) -> V` puts V in BOTH param (contravariant) and return
(covariant) ⇒ pyright infers V **INVARIANT**. `TensorProductOperator` is bare ⇒ its
inherited typevar is `Unknown`, but its OWN `apply` is `(x: ndarray) -> ndarray`
(op.py:1304). Consumer `_as_boundary(op: LinearOperator)` (param = implicit
`LinearOperator[Unknown]`) demands `(x: Unknown) -> Unknown`. Under invariance the
`ndarray`↔`Unknown` mismatch fails the reverse-direction check → `"LinearOperator[V@…]"
is not assignable to "TensorProductOperator"`. Same for `ScaledOperator.__init__(op:
LinearOperator[V])` at :207/:224/:240 (the `ScaledOperator(albedo, base)` wrap).

## 8.6 Candidate re-typings (ranked) + RUNTIME-TOUCH flags

- **(a) ⭐ Parametrize the bare axis-primitives `LinearOperatorMixin[np.ndarray]`**
  (`TensorProduct`/`Permutation`/`Mask`/`PeriodicWrap`/`Diagonal`/`RankOne`/
  `SumOfTensorProducts` + the SN-angular ones). Their `apply` ALREADY says
  `np.ndarray` → V = `np.ndarray` is exact, NO runtime change, override signatures
  already match. Fixes every `TensorProductOperator`-not-assignable error.
  Blast radius ~12 type-only lines.
- **(b) Widen consumers to `LinearOperator[Any]`** — smallest (~3 lines) but a
  SUPPRESSION, erases V at every composition boundary, does NOT fix :125. LOW
  QUALITY (Cardinal Rule 1/2) — REJECT.
- **(c) Split `SolvableOperator(LinearOperator[V], Protocol)` / `TransposableOperator`**
  — the RIGHT fix for the §8.4 "has solve iff CAP_SOLVE" gap, but LARGE cross-cutting
  blast radius (every solve-defining class + composer propagation + iteration gate).
  **Defer to a dedicated issue; do NOT fold into B4** (~3–4× scope blowup).
- **(d1) Retype `_as_boundary`'s param to `LinearOperatorMixin`** (which DOES declare
  `block_role`) instead of the bare `LinearOperator` Protocol — type-honest (all
  realized ops are mixin subclasses), cheap, fixes :125. (d2 = add `block_role` to
  the Protocol — REJECT, pollutes the minimal universal protocol with an SN concept.)

### ⭐ RECOMMENDATION (min blast radius, preserves `(L + C − S − F/k) @ psi`, zero runtime touch)
**(a) + (d1).** Parametrize the bare numerics + SN-angular primitives as
`LinearOperatorMixin[np.ndarray]`; retype `_as_boundary`'s param/return as
`LinearOperatorMixin`. **Total ~12–15 type-only lines, all class headers / param
annotations in `operator.py`, `sn/angular_operator.py`, `multiplication_operator.py`,
`sn/boundary_realizer.py`. NO call site, NO runtime line.** Defer (c). Do NOT touch
`block_role`/`capabilities` runtime declarations (§8.2/§8.3 — current shape is
load-bearing).

### ⚠ RUNTIME-TOUCH FLAGS (a B4 proposal should NOT cross these)
- `StreamingOperator`/`InvertibleOperator` `["FullField"]` is a forward-ref STRING
  proxy for `TimedFullField` (layering forbids the import). Tightening to the true
  type forces a runtime import / `TYPE_CHECKING` alias — leave as-is for B4; note
  for a future numerics↔transport layering pass.
- `ProjectionOperator`/`ReconstructionOperator` are NON-endomorphic (V_in ≠ V_out:
  flux↔moments). The single-V `apply(x: V) -> V` does NOT fit them; they live
  OUTSIDE the `(L+C−S−F−B)` algebra. If the carve touches them it has left B4 scope.
- `ScatteringOperator`/`FissionOperator` `apply` are `singledispatchmethod` over a
  small union (per §2) — annotate the dispatch BASE to the within-group V, keep
  registered arms concrete; do NOT force one TypeVar over the whole union.
