# #226 — Operator-algebra `Generic[V]` design proposal

Status: **SUPERSEDED (2026-06-19) — DO NOT IMPLEMENT.** A live design discussion
established that the "ndarray operator strays" are ARCHITECTURE GAPS, not a
genuine second vector population: `TimedFullField` is already the cross-method
generic container (`bulk: BulkField` ⊕ `boundary: BoundaryField`) and can carry
MOMENTS (`HarmonicMomentField` is a `BulkField`; only a `BoundaryMomentField` is
missing). With the container as the vector type, `MomentProjection` is
container-ENDOMORPHIC (`TimedFullField→TimedFullField`), so the `Generic[V]`-over-
heterogeneous-types premise here is moot. The correct work is the **field-typed
operator-algebra campaign** — GH **#256**, detailed plan at
`.claude/plans/field_typed_operator_algebra_campaign.md` (keystone:
P3.6 promote `TimedFullField`→`DirectSumSpace` `Field`; + `BoundaryMomentField`;
type projection/scattering container→container; scipy-as-typed-serialization;
audit the matrix primitives). This file is retained only as the rejected
alternative (its layering + multi-vector-space analysis is still accurate input).
Original proposal evidence: `.claude/agent-memory/explorer/issue_226_operator_generics_map.md`.

## The problem

`orpheus/numerics/operator.py` types the operator interface on `np.ndarray`:

- `LinearOperator(Protocol)` `.apply(self, x: np.ndarray) -> np.ndarray` (op.py:317)
- `LinearOperatorMixin`, `OperatorSum`, and the composer family all `apply(x: np.ndarray) -> np.ndarray`.

But the SN operators (`StreamingOperator`, `CollisionOperator`,
`InvertibleOperator`, `SNBoundaryOperator`, …) implement
`apply(psi: TimedFullField) -> TimedFullField`. Result: **7
`reportIncompatibleMethodOverride`** (the override changes the parameter type)
plus a **`reportArgumentType` cascade** wherever an SN field-operator flows into a
parameter annotated `LinearOperator` / `q_ext: np.ndarray`.

## The design — single, unbounded `Generic[V]`

The operator algebra is genuinely generic over the **vector/field type** it acts
on. Make that explicit (coding-elegance Pattern 1/4 — the algebra IS generic over
its vector space; the notation matches the domain):

```python
V = TypeVar("V")                              # UNBOUNDED (see below)

class LinearOperator(Protocol[V]):
    def apply(self, x: V) -> V: ...           # was np.ndarray -> np.ndarray
    # capabilities / domain / codomain UNCHANGED — V-independent
    #   (domain/codomain are FunctionSpace, not the vector)

class LinearOperatorMixin(Generic[V]): ...
class OperatorSum(LinearOperatorMixin[V]): ...   # a.apply(x) + b.apply(x), x: V
# composer family Generic[V]: OperatorProduct, ScaledOperator, IdentityOperator,
#   ZeroOperator, _AdjointOperator, TensorProductOperator, SumOfTensorProducts…
```

Concrete operators parameterize at the leaf:
- ndarray matrix operators → `LinearOperator[np.ndarray]`
- SN field operators → `LinearOperator[TimedFullField]` (e.g.
  `class StreamingOperator(LinearOperatorMixin[TimedFullField])`).

### Why `V` is UNBOUNDED (not `bound=Field`)

Both load-bearing carriers fail a `Field` bound:
- `np.ndarray` is not a `Field`.
- `TimedFullField` is **deliberately NOT a `Field` subclass** (it is a structured
  composite `bulk`/`boundary`/`history`; `timed_full_field.py:148`).

A `TypeVar("V", bound=Field)` would exclude the two most important Vs. A structural
`Vector` protocol bound (`__add__`/`__sub__`/scalar `__mul__`) is a larger design
not needed for #226 → use **unbounded `V`**.

### Why it is CLEAN (verified by the explorer)

- **No mixed sums.** No production `OperatorSum`/`OperatorProduct`/`&` ever combines
  an ndarray-operator with a field-operator. The within-group algebra is a
  *resolvent + variadic gains* split composed at the **field level**
  (`out = L.apply(psi); for g in gains: out = out - g.apply(psi)`,
  iteration.py:743), NOT one heterogeneous `OperatorSum`. So a single `V` per
  composition always type-checks.
- **Pure typing / runtime-inert.** `Generic[V]`/`Protocol[V]` emit no runtime code.
  No production class combines `LinearOperatorMixin` + `RegistryMixin`, so the
  `__init_subclass__(key=)` registry machinery is untouched (the one example in a
  docstring is not live). The `_BlockRoleMeta` metaclass is only on role markers,
  not on the Mixin.

## What it clears + the consumer retype

1. **The 7 `reportIncompatibleMethodOverride`** vanish once the Protocol's
   `apply(x: V) -> V` admits the SN operators' `TimedFullField` apply.
2. **The argtype cascade** clears by retyping the "generic over both" consumers'
   `q_ext: np.ndarray` → `V` (they already run on both flat ndarray L0 and typed
   `TimedFullField` L1 via the ravellable `to_flat`/`from_flat` protocol):
   `SourceIteration.__init__` / `KrylovAcceleration.__init__` /
   `KEigenvalue.__init__` (iteration.py:391/656/934), `as_scipy_linop` (op.py:1661).
3. The dunders' existing `# type: ignore[arg-type]` (op.py:415-485) can be removed.
4. **singledispatch S/F operators** (`ScatteringOperator`, `FissionOperator`) and
   union-typed `LegendreMomentScattering`: annotate the **base** dispatch method to
   the within-group V (`TimedFullField`); keep the registered arms' concrete
   annotations (`ScalarFlux`/`AngularFlux`/`np.ndarray`). Do NOT force one TypeVar
   over the whole dispatch union.

## Risks & mitigations

- Variance: `apply(x: V) -> V` is invariant in V, which is correct here (input and
  output are the same space). No covariance/contravariance subtlety.
- `@dataclass` operators: adding `Generic[V]` as a base is fine; **preserve** the
  existing plain-attr `block_role` discipline (it is intentionally NOT
  `ClassVar[...]` under `from __future__ import annotations`, op.py:391-396).
- Verification: pure-typing, so the gate is (a) pyright CLI count drops with NO new
  errors, (b) import smoke on operator.py + iteration.py + sn/operator.py, (c) a
  targeted run of `tests/sn/operators` + `tests/sn/solve` + the numerics operator
  tests stays green (no behavioural change).

## Files

`orpheus/numerics/operator.py` (Protocol + Mixin + OperatorSum + composer family),
`orpheus/numerics/iteration.py` (the 3 consumers + `q_ext`),
`orpheus/sn/operator.py` (leaf parameterization + the 2 OperatorSum overrides),
`orpheus/numerics/projection.py` (4 override errors — ProjectionOperator/
Reconstruction parameterize to their field), `orpheus/numerics/spaces/full_field_space.py`
(1 override). singledispatch base annotations in `sn/scattering.py`, `sn/fission.py`.

## Alternative considered (and rejected)

Widen the Protocol's `apply` to `Any`/a union instead of `Generic[V]`. Rejected:
loses the input-type == output-type guarantee, is less precise, and does not read
like the domain. `Generic[V]` is the elegant, faithful encoding.
