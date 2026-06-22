---
name: issue-226-operator-generics-map
description: Map for making numerics/operator.py's LinearOperator algebra Generic[V] over the vector/field type (pyright #226). Interface surface, ndarray-vs-field operator buckets, mixed-sum check (NONE), consumer scope, field hierarchy + TypeVar bound, Generic/RegistryMixin runtime-safety.
metadata:
  type: project
---

# Issue #226 — `LinearOperator` algebra `Generic[V]` design map (read-only investigation)

Point-in-time `file:line` are current at 2026-06-18 on `main`. The SHAPE
(buckets, mixed-sum verdict, TypeVar bound) is the durable part.

## 0. Headline verdict (for the design)

- A single `Generic[V]` parameterization of `LinearOperator` / `LinearOperatorMixin` /
  `OperatorSum` (+ the composer family) is **CLEAN**. **No production `OperatorSum`
  ever mixes an ndarray-operator with a field-operator.** All real sums are homogeneous.
- `V` must be **UNBOUNDED** (`V = TypeVar("V")`), NOT `TypeVar("V", bound=Field)`.
  Two reasons: (1) `np.ndarray` is a legal V and is not a `Field`; (2) the dominant SN
  field carrier `TimedFullField` is deliberately **NOT a `Field` subclass**
  (`orpheus/transport/timed_full_field.py:148` — "NOT a Field subclass at the typed-class
  level"). A `Field` bound would exclude both. If a bound is wanted at all, it would have
  to be a structural `Vector` protocol (`__add__`/`__sub__`/scalar `__mul__`), not `Field`.
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
