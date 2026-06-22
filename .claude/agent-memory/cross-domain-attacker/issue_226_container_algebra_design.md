---
name: issue-226-container-algebra-design
description: Design memo — uniformly container-typed operator algebra (TimedFullField as the universal vector). Resolves Protocol shape, V→V vs V→W, matrix-primitive scope line, BoundaryMomentField, P3.6 DirectSumSpace, scipy serialization, cross-method + sequencing/risk. Inputs for a post-compaction architecture campaign (#226).
metadata:
  type: project
---

# Issue #226 — Container-typed operator algebra: design resolution

Frame-detection task, but the deliverable is a design memo. Output is the
seven resolved questions, each as a concrete recommendation with reasoning
grounded in the actual code (operator.py, timed_full_field.py, field.py,
fields/_bases.py, space.py, full_field_space.py, projection.py, iteration.py).

**Why:** ORPHEUS's operator algebra is currently typed on `np.ndarray` in the
`LinearOperator` Protocol (`apply(self, x: np.ndarray) -> np.ndarray`) while
the production carriers are typed `Field`/`TimedFullField` composites. The
ndarray annotation is a documentation-level lie: SN ops already consume and
return `TimedFullField`. The campaign uniforms the algebra onto the container.

**How to apply:** this is the authoritative input to the #226 campaign plan.
Read alongside `project_wave_o_operator_algebra` (the algebra's history),
`field_role_typing_faceflux_frames` (DEC/cochain frame for the leaves),
`issue_208_delta_psi_affine_frames` (the displacement-typing already shipped),
and `eigenvalue_posing_layering_frames` (the layering the resolvent sits in).

---

## Established findings (verified against code, do not re-litigate)

1. **`TimedFullField` (transport/timed_full_field.py:123)** is a frozen
   `kw_only` dataclass = `bulk: BulkField` ⊕ `boundary: BoundaryField` ⊕
   `_history` ⊕ `history_depth`. It is NOT a `Field` subclass; it carries
   delegate dunders (`__add__/__sub__/__neg__/__mul__/__rmul__/__truediv__`)
   that propagate to bulk+boundary and DROP history. It already has the
   flat-vector protocol `to_flat()` / `from_flat(flat, template)` (lines
   439/461) consumed by `numerics.iteration`.

2. **The leaf algebra is the single source of truth.** `TimedFullField.
   _check_partner` (line 260) gates ONLY same-container-class; the
   bulk/boundary leaf dunders (`Field._check_partner`, field.py:235) enforce
   class-identity-is-units-identity + same-space + same-mesh, AND host the
   affine torsor gate (`flux+flux→TypeError`, `flux−flux→displacement`,
   `flux+displacement→flux`) shipped in #208. The container must NOT
   re-implement any of this — it delegates.

3. **`FullFieldSpace` (numerics/spaces/full_field_space.py) ALREADY EXISTS**
   as a `FunctionSpace` subclass = `V_bulk ⊕ V_trace` with a per-block metric
   (`apply_metric` / `apply_inverse_metric` / `inner_product` overrides that
   operate on a composite field duck-typed as `.bulk.values` / `.boundary.
   values`). It is the Wave-O G-adjoint carrier. **It is wired ONLY to the
   adjoint metric path** — it is NOT yet `TimedFullField`'s `.space`, and
   `TimedFullField` does not expose `.space`. This is the half-built keystone
   of P3.6.

4. **`MomentProjection` (numerics/projection.py:241)** is a genuine
   rank-changing einsum `(N,...) → (L+1,2L+1,...)` typed `V→W` at the
   ARRAY/leaf level (`domain` = `(N,)` ordinate space, `codomain` =
   `SphericalHarmonicSpace`). At the CONTAINER level the user's insight holds:
   `AngularFlux→HarmonicMomentField` is an endomorphism `TimedFullField→
   TimedFullField` (inner `.bulk` type differs at runtime; the container type
   is invariant). BUT note `MomentProjection` today acts on a RAW ARRAY, not a
   container — it is consumed inside `InvertibleOperator.solve_moments`, not as
   a public algebra leaf.

5. **scipy is a serialization boundary, not a vector space.** TWO adapters
   exist: (a) `as_scipy_linop` (operator.py:1661) for bare-ndarray operators,
   wraps `op.apply` as scipy `matvec`; (b) `KrylovAcceleration.solve`
   (iteration.py:695) builds its OWN `A_matvec` closure that ravels via the
   `_is_ravellable` protocol (`to_flat`/`from_flat`), composes `L.apply −
   Σgᵢ.apply` on the typed carrier, and re-ravels. The production SN matvec
   goes through (b), NOT (a). scipy sees a flat `(n,n)` matvec in both.

---

## Q1 — The universal vector type & Protocol shape

### Recommendation: **(c) a structural `Vector` Protocol in numerics**, with the
Protocol's `apply` typed `apply(self, x: V) -> V` where `V` is a `TypeVar`
bound to that `Vector` Protocol. This is **(c) as the foundation, (a) as the
typing surface, NOT (b) Generic-per-operator.**

```python
# orpheus/numerics/vector.py  (NEW — L1, below operator.py)
from typing import Protocol, TypeVar, runtime_checkable

@runtime_checkable
class Vector(Protocol):
    """Structural contract for an element of the algebra's vector space.

    Satisfied by np.ndarray (duck-typed: ndarray has __add__/__sub__/
    scalar __mul__) AND by every Field leaf AND by TimedFullField.
    The algebra's apply/solve speak THIS, never np.ndarray.
    """
    def __add__(self, other: "Vector") -> "Vector": ...
    def __sub__(self, other: "Vector") -> "Vector": ...
    def __mul__(self, scalar: float) -> "Vector": ...
    # zero is obtained structurally (0.0 * x) — no abstract `zero` needed;
    # ZeroOperator already uses `0.0 * x`, which every Vector satisfies.

V = TypeVar("V", bound=Vector)
```

```python
# operator.py — Protocol retype
@runtime_checkable
class LinearOperator(Protocol):
    capabilities: frozenset[str]
    @property
    def domain(self) -> Optional["FunctionSpace"]: ...
    @property
    def codomain(self) -> Optional["FunctionSpace"]: ...
    def apply(self, x: V) -> V: ...     # was: np.ndarray -> np.ndarray
```

### Why (c)+(a), not (b)

- **(b) `Generic[F]` per-operator is the trap.** It would force every
  composer (`OperatorSum`, `OperatorProduct`, `ScaledOperator`) to be
  `Generic[F]` and would make `L + C − S − F − B` carry a parametrised type
  that the SN-vs-CP cross-method goal then has to UNIFY at the use site. But
  the algebra is ALREADY runtime-duck-typed (the module docstring is explicit:
  "Shape and dtype are deliberately not part of the protocol"). The
  cross-method generality is achieved by the carrier being one container type
  with a runtime-varying inner leaf, NOT by static parametrisation. A
  per-operator `Generic[F]` adds type ceremony that buys nothing the runtime
  `_check_partner` gates don't already buy, and it FIGHTS the existing design
  (mixed bulk-types inside one `TimedFullField` are the moment-endomorphism
  case — a single static `F` cannot express "in: flux-bulk, out: moment-bulk,
  same container").

- **(a) single `apply(x: F) -> F` is the right SURFACE but needs (c) to name
  `F`.** The endomorphism claim (Q2) means the codomain container type EQUALS
  the domain container type for every operator that participates in the
  `(L+C−S−F−B)ψ=q` algebra. So a single-type-var `apply(x: V) -> V` is
  honest. The `TypeVar V` must be bound by SOMETHING that lives in numerics
  (layering: `numerics/` is below `transport/`, so the Protocol cannot name
  `TimedFullField`). That something is the structural `Vector` Protocol — it
  is the numerics-level abstraction `TimedFullField` conforms to WITHOUT
  numerics importing transport.

- **`Vector` is structural (duck-typed), matching the existing `LinearOperator`
  Protocol style and the `_is_ravellable` duck-type in iteration.py.** No ABC
  import, no inheritance obligation on `np.ndarray`. ndarray satisfies
  `Vector` for free; `Field` satisfies it via its dunders; `TimedFullField`
  satisfies it via its delegate dunders. This is the SAME move the codebase
  already made for `_is_ravellable` (duck-typed protocol preserving the
  numerics↛transport layering) — promote it from an ad-hoc helper to a named
  Protocol.

### Elegance payoff
- **Structure-exposing**: the Protocol stops lying. `apply(x: np.ndarray) ->
  np.ndarray` on an operator that consumes `TimedFullField` is a documented
  falsehood; `apply(x: V) -> V` over `Vector` states the real contract — "I
  map a vector to a same-typed vector."
- **Expressive**: `(L + C − S − F − B)ψ = q` reads with ψ and q both `Vector`,
  no `np.ndarray` cast at any algebra site.
- **Structurally-simpler**: ONE type var, ZERO per-operator generics, ONE new
  numerics Protocol. The composers stay non-generic.
- **Algorithmic-advantage**: none directly (pure typing) — and that is the
  POINT (Cardinal Rule 1: bit-identical refactor).

---

## Q2 — V→V vs genuine V→W: confirm the endomorphism, name the exceptions

### The container-endomorphism claim is CONFIRMED for the iterated algebra,
with **three genuine exception categories** that are NOT endomorphisms and
must be typed separately.

**Endomorphism (the `(L+C−S−F−B)` algebra):** every operator that
participates in the source-iteration / Krylov / power-iteration loop maps the
state container to a same-typed state container. This is forced by the loop
structure itself: `ψ_{n+1} = L⁻¹(Σgᵢψ_n + q)` requires `gᵢψ`, `q`, and `ψ`
all addable, hence all the same container type (`TimedFullField` carrying the
same leaf classes). The moment-projection "V→W" dissolves at the container
level exactly as the user states: `MomentProjection∘(container)` changes
`.bulk` from `AngularFlux` to `HarmonicMomentField` but the wrapper stays
`TimedFullField`. So the resolvent algebra is an endomorphism algebra on
`TimedFullField`. **This is the spine of the campaign.**

### Exception 1 — Functionals `V → scalar` (reductions)
The keff/production estimators (`_default_keff_estimator`,
`_default_production_estimator`, iteration.py:250/271) ARE `R = <Σ, ψ> :
TimedFullField → float`. `inner_product` / `norm` / `l2` (Field, FunctionSpace,
FullFieldSpace) are the same shape. These are NOT `LinearOperator`s and must
NOT be forced into the `apply(x:V)->V` mould.

**Type as:** a separate `Functional` Protocol `evaluate(self, x: V) -> float`
(distinct method name — a functional is a covector, not an operator). The
inner-product machinery already lives on `FunctionSpace`/`FullFieldSpace`; a
`Functional` leaf that wraps `<Σ,·>` is the typed face of it. Keep them out of
`LinearOperator` — the capability-set + composer closure laws are meaningless
for a `V→scalar` map (you cannot `OperatorSum` an operator with a functional).

### Exception 2 — Trace/restriction `V → V_boundary` and injection `V_boundary → V`
The boundary-extraction step planned in Wave O (`O.4a.2`: extract `B` from the
sweep so `B` is a first-class `(L+C−S−F−B)` leaf) needs the trace operator
ι*: full → boundary and its adjoint ι_*: boundary → full. At the CONTAINER
level full = `TimedFullField`, boundary = `BoundaryField` (a bare leaf, NOT a
container). So ι* is genuinely `TimedFullField → BoundaryField` — a different
container OUT than IN. This is the DEC/cochain restriction named in
`field_role_typing_faceflux_frames` ("interior FaceFlux = 1-cochain, boundary
trace = its restriction ι*; absorption = ι_*∘ι*=id").

**Type as:** these are the genuine `V→W` operators. BUT — critically — the
`B` boundary LAW that appears in the algebra is `B: V_boundary → V_boundary`
(an endomorphism on the trace block, `BlockRole.BOUNDARY`), realised on the
DIRECT-SUM container as the `A_ss` block of a 2×2 that is the IDENTITY on the
bulk block. So inside `(L+C−S−F−B)` acting on `TimedFullField`, `B` IS an
endomorphism `TimedFullField→TimedFullField` (it touches only `.boundary`).
The pure trace ι* `full→boundary-only` is an INTERNAL implementation detail of
how `L` couples blocks (it reads the inflow trace to seed the sweep, writes
the outflow trace) — it is the off-diagonal `A_bs`/`A_sb` of `L`'s `FULL`
block matrix, NOT a public algebra leaf.

**Decision:** the public algebra needs ONLY the endomorphism signature
`apply(x: V) -> V`. The trace ι* and injection ι_* are PRIVATE to `L`'s apply
(they are how a `FULL`-role operator moves data between the two blocks of the
SAME container) and to the boundary-realiser. They do NOT get a public
`apply(Vin)->Vout` `LinearOperator` surface. The block-role machinery
(`BlockRole.BULK/FULL/BOUNDARY`, `_join_block_roles`) ALREADY encodes this
correctly — a `FULL` operator's apply is `TimedFullField→TimedFullField`; the
block-coupling is internal.

### Exception 3 — Projection/reconstruction as PUBLIC leaves
`MomentProjection` (V→moment-bulk) and `HarmonicMomentReconstruction`
(moment-bulk→V) when used as STANDALONE public operators (e.g. a future PN
solver, the §10 moment-space coupling) ARE container-changing at the leaf
level but container-INVARIANT at the wrapper level. Today they act on raw
arrays inside `solve_moments`. If they are ever lifted to public algebra
leaves, they type as endomorphisms on `TimedFullField` whose `.bulk` leaf
class changes — which the single `apply(x:V)->V` signature ALREADY covers
(the container type is invariant; the runtime leaf differs). **No special
category needed** — this is the user's original insight, and it is why (a)+(c)
beats (b): a static `Generic[F]` could NOT express "same container, different
inner leaf", but the structural `Vector` bound CAN (the container is the
`Vector`).

### Elegance payoff (the Q2 resolution)
- **Structure-exposing**: separates the THREE things the current ndarray
  signature conflates — endomorphisms (the algebra), functionals (reductions),
  and genuine block-trace maps (internal to `L`). Names each.
- **Structurally-simpler**: the PUBLIC algebra needs exactly ONE operator
  signature (`V→V`) plus ONE functional signature (`V→scalar`). The
  scary-sounding `V→W` cases are either internal (trace) or already-covered
  (projection = same container, different leaf).

### First test (discriminates the endomorphism claim)
Build `(L + C − S − F)` for a 2-group SN slab, assert
`type(op.apply(tff)) is TimedFullField` AND `type(op.apply(tff).bulk) is
type(tff.bulk)` for L, C, S, F individually and the composite. Then build a
keff functional and assert it is NOT an instance of `LinearOperator` (no
`apply`, no `capabilities`) — discriminates the "functional is a separate
category" decision from the "force everything into apply" anti-design.

---

## Q3 — The matrix primitives: where the scope line falls

The seven primitives (`PermutationOperator`, `DiagonalOperator`,
`TensorProductOperator`, `SumOfTensorProductsOperator`, `RankOneOperator`,
`IdentityOperator`, `ZeroOperator`) are **axis-level array algebra on a
SINGLE tensor**. They all carry an `axis: int` (except Identity/Zero) and act
via `np.take` / broadcast-multiply / einsum on ONE ndarray.

### The dividing line: **axis-operators stay FLAT primitives used INSIDE a
typed operator's apply; the COMPOSITE state operators are container-typed.**

| Primitive | Verdict | Reasoning |
| --- | --- | --- |
| `PermutationOperator` | **flat, private** | acts on one axis of one array (`np.take(x, perm, axis)`). The reflective-BC ordinate swap, periodic shift. It is a sub-vector operation on the `.bulk.values` or `.boundary.values` ndarray INSIDE a typed BC/streaming operator's apply. Never a public `LinearOperator` in the state algebra. |
| `DiagonalOperator` | **flat, private** | broadcast-multiply on one axis (`AngularWeightMatrix` W, σ_t multiply, volume weights). The metric primitive. Lives inside `CollisionOperator`/`StreamingOperator`/the metric on a space. |
| `RankOneOperator` | **flat, private** | `χ⊗νΣ_f` fission emission on the group axis of one array. Lives INSIDE `FissionOperator.apply`, which IS container-typed (`TimedFullField→TimedFullField`). |
| `TensorProductOperator` | **flat, private** | the `D⊗Ω⊗I` streaming factorisation acts on one array's axes. Lives inside `StreamingOperator`. |
| `SumOfTensorProductsOperator` | **flat, private** | `Σ Pℓ⊗Σ_{s,ℓ}` scattering. Lives inside `ScatteringOperator`. |
| `IdentityOperator` | **BOTH** | already container-safe: `apply(x)=x` is type-preserving for ANY `Vector`. It is the only primitive that is ALSO a legitimate public algebra leaf (the bulk identity, or the `A_ss` identity for a vacuum boundary). |
| `ZeroOperator` | **BOTH** | already container-aware: `codomain_zero` hook (lines 947–955) was ADDED precisely so the zero fission slot returns a typed `AngularSourceSink`-bulk zero, keeping `S.apply+F.apply+q` a closed `TimedFullField` sum. So `ZeroOperator` is ALREADY the campaign's pattern in miniature — a flat-default that goes container-typed via a codomain hook. |

### The decision rule (this is the campaign's scope statement)
**An operator is a public container-typed `LinearOperator` iff its action is
on the COMPOSITE STATE (the full phase-space distribution + its boundary
trace). An operator is a flat private primitive iff its action is on ONE
tensor axis of ONE leaf's `.values`.** The seven primitives are the second
kind, with `Identity`/`Zero` as the two that span both (because their action
is axis-agnostic / type-agnostic).

**Concretely:** the campaign retypes the SN-level operators (`StreamingOperator`,
`CollisionOperator`, `ScatteringOperator`, `FissionOperator`,
`InvertibleOperator`, the realised boundary laws `B`) to `apply(x:
TimedFullField) -> TimedFullField` (they MOSTLY are already, via the
container's delegate dunders + the leaf operators). The seven `numerics`
primitives KEEP `apply(x: np.ndarray) -> np.ndarray` — but that signature is
now HONEST, because their domain genuinely IS a bare ndarray axis, and `np.
ndarray` satisfies the `Vector` bound. **No retype needed for the seven** —
the `Vector`-bound `V` unifies "bare ndarray axis-op" and "container state-op"
under one Protocol without forcing the axis-ops to fake a container.

### Elegance payoff
- **Structure-exposing**: makes explicit the two-layer reality the codebase
  already lives — axis-algebra (numerics primitives, einsum-level) vs
  state-algebra (transport operators, container-level). The Protocol's `Vector`
  bound is exactly the join of those two layers.
- **Structurally-simpler**: the campaign does NOT have to container-ify
  `PermutationOperator` et al. — a tempting but wrong scope expansion that
  would force a `TimedFullField` wrapper around a single-axis `np.take`. The
  line is "composite state ⇒ container; single axis ⇒ ndarray."
- **Algorithmic-advantage**: keeping the primitives flat preserves the
  zero-copy einsum/broadcast performance (the `_MovingFrontier` window, the
  `R·Λ·M` scattering primitive) — container-wrapping them would interpose
  dataclass allocation on the hot axis-ops.

### First test
For each of the seven primitives, assert it acts correctly on a bare
`np.ndarray` AND assert it is NOT in the retyped state-operator set (grep that
no production `solve`/Krylov site passes a `TimedFullField` directly to
`DiagonalOperator.apply` — they pass `.bulk.values`). For `Identity`/`Zero`,
assert `apply(tff)` returns a `TimedFullField` (the span-both property). This
discriminates "container-ify the primitives" (wrong) from "the Vector bound
unifies both layers" (right).

---

## Q4 — `BoundaryMomentField` + `MomentFullField`

### `BoundaryMomentField` MUST be: a `BoundaryField` subclass storing
moment-resolved trace values, mirroring `BoundaryFlux` but with a
moment-widened per-face slot. Specifically it is the boundary-locus sibling of
`HarmonicMomentField`: where `BoundaryFlux` stores per-face angular-flux trace
values on the `TraceSpace` flat layout `(total_size,)`, `BoundaryMomentField`
stores per-face MOMENT trace values — the moment-projection of the boundary
angular flux.

Concretely it inherits `BoundaryField` (fields/_bases.py:480) → gets the
`mesh` binding, the `TraceSpace` contract, `face_view`, the
`zeros_on`/`from_face_arrays`/`from_mesh` factories for free. It needs:
- a moment-widened layout: the `TraceSpace`'s `FaceLayout` slot per face is
  `(L+1, 2L+1, ...)` instead of `(N, ...)` (the angular ordinate axis replaced
  by the `(ℓ,m)` moment axes), OR a `TraceSpace × SphericalHarmonicSpace`
  composition mirroring how `HarmonicMomentField` composes
  `... × SphericalHarmonicSpace`.
- its own `_SPACE_NAME`-equivalent / `UNITS` (same units as `BoundaryFlux` —
  it is a flux trace, just moment-resolved).

**Why it is needed:** the moment-windowed solve path (`solve_moments`,
`_MomentWindowedResolvent`, the 5c moment-output sweep — see
`project_wave_o_operator_algebra`) produces a `HarmonicMomentField` bulk. For
that to live in a `TimedFullField` (so the moment endomorphism `TimedFullField
→TimedFullField` of the user's insight is type-CLOSED), the container's
`.boundary` partner must ALSO be moment-typed — otherwise you have a
`TimedFullField(bulk=HarmonicMomentField, boundary=BoundaryFlux)` whose two
halves disagree on angular representation, and the container's mesh-identity
check passes while the angular-rep contract silently breaks. **The missing
`BoundaryMomentField` is a Smell-16 instance:** the moment solve path stores
its boundary as an angular `BoundaryFlux` (or drops it), so the moment state
has TWO incompatible angular representations bridged by hand. Name the missing
type; the bridge IS that type un-named.

### `MomentFullField` alias: **NOT a separate class — a named factory /
type alias.** Recommend a module-level alias + a `TimedFullField` classmethod:

```python
# transport/timed_full_field.py
def moment_full_field(
    bulk: HarmonicMomentField, boundary: BoundaryMomentField, *, history_depth=2,
) -> TimedFullField:
    """Named constructor for the moment-carrying composite (#226)."""
    return TimedFullField(bulk=bulk, boundary=boundary, history_depth=history_depth)
```

Reasoning: `TimedFullField` is DELIBERATELY cross-method generic (its whole
docstring is about being the ONE container for SN/CP/MoC/diffusion × flux/
moment). Minting a `MomentFullField` SUBCLASS would re-introduce the
per-combination class explosion the container was built to kill (Cardinal Rule
2). The container type stays ONE; the inner leaf pair `(HarmonicMomentField,
BoundaryMomentField)` distinguishes it at runtime. A named CONSTRUCTOR
(factory function or classmethod) gives the grep-signal + the type-checker
hint without a subclass. The cross-class arithmetic gate (`_check_partner`
class-identity + leaf class-identity) already prevents a flux-`TimedFullField`
from adding to a moment-`TimedFullField` (their `.bulk` classes differ →
`Field._check_partner` raises) — so safety needs no new class.

### Elegance payoff
- **Structure-exposing**: surfaces that the moment solve has an unbuilt
  boundary partner (the Smell-16 silent second representation).
- **Structurally-simpler**: ONE container class, an alias not a subclass —
  refuses the `{flux,moment}×{SN,CP,MoC}` class grid.
- **Expressive**: `moment_full_field(φ_lm, b_lm)` reads as the moment state.

### First test
Construct `TimedFullField(bulk=HarmonicMomentField.zeros_on(mesh),
boundary=BoundaryMomentField.zeros_on(mesh))` and assert the container
algebra closes (`tff + tff` works, `tff_moment + tff_flux` RAISES via the leaf
class gate). Then assert `solve_moments` can return this closed container
instead of `(moment_buf, None)` — discriminates "moment state is type-closed"
from the current "boundary dropped / None" half-state.

---

## Q5 — P3.6 DirectSumSpace Field promotion (the layering keystone)

### What `TimedFullField`-as-`Field` entails

`Field` (numerics/field.py) is `(values: NDArray, space: FunctionSpace)` with
dunder algebra gated by `_check_partner` (class + space). For `TimedFullField`
to BE a `Field`:

1. **`space` must be a `FullFieldSpace`** (which already exists and already is
   a `FunctionSpace` with the per-block metric). `TimedFullField.space` →
   `FullFieldSpace.from_blocks(bulk.space, boundary.space)`. This is the
   missing wire: `FullFieldSpace` is built for the adjoint but never installed
   as the container's `.space`.

2. **`values` is the problem.** `Field` assumes `values: NDArray` with
   `values.shape == space.shape`. `TimedFullField`'s "values" are a STRUCTURED
   pair, not one ndarray. Two sub-options:
   - **(5a) `values` = the flat direct-sum** `concat(bulk.ravel, boundary)` =
     exactly `to_flat()`. Then `space.shape = (n_bulk+n_trace,)` (which
     `FullFieldSpace.from_blocks` ALREADY computes, line 167). The bulk/boundary
     leaves become DERIVED views (`from_flat`). Clean `Field` conformance but
     forces a ravel/reshape on every leaf access — perf regression on the hot
     `.bulk` path.
   - **(5b) `values` stays the structured pair; `Field` gains a
     `DirectSumStorage` value type.** The `Field` ABC relaxes "values is one
     ndarray" to "values is a `Vector`" (the Q1 Protocol). `TimedFullField`
     keeps `.bulk`/`.boundary` as direct attributes; its dunder algebra
     already matches `Field`'s (propagate to members). `space` is
     `FullFieldSpace`; `_check_partner` reuses `Field`'s class+space gate.

   **Recommend (5b).** It preserves the zero-copy `.bulk` access, and
   `FullFieldSpace` already operates on the STRUCTURED composite (its
   `apply_metric` reads `x.bulk.values` / `x.boundary.values`, NOT a flat
   vector). The flat `to_flat`/`from_flat` stays as the SERIALIZATION
   boundary (scipy), not the storage model. (5a) would couple storage to
   serialization — the exact category error Q6 warns against.

3. **What `DirectSumSpace(bulk_space ⊕ boundary_space)` provides:** it IS
   `FullFieldSpace` (rename or alias — `FullFieldSpace` is SN-named, a generic
   `DirectSumSpace` is the cross-method name). It provides: the flat identity
   `shape=(n_bulk+n_trace,)` for composer composability checks; the per-block
   metric (`apply_metric`/`apply_inverse_metric`/`inner_product`); and (new)
   the `__mul__`/`dual` algebra inherited from `FunctionSpace`. The `Field`
   `_from_balance` named-composition machinery lifts to the container for free
   (a `FullResidual` = `from_balance(Aψ, q)` on the composite).

### Sequencing: **P3.6 FIRST, then the Protocol retype is trivial.**

Do P3.6 (promote `TimedFullField` to `Field` on `FullFieldSpace`/`DirectSumSpace`)
BEFORE Q1's Protocol retype. Reason: once `TimedFullField` IS a `Field`, and
`Field` satisfies the `Vector` Protocol (it has the dunders), then `apply(x:V)
->V` is satisfied by `TimedFullField` THROUGH its `Field`-ness — the Protocol
retype becomes a one-line annotation change with no new conformance work. If
you retype the Protocol FIRST, you must separately prove `TimedFullField`
satisfies `Vector` (it does, via delegate dunders) and then re-prove it again
after P3.6 — double work. **But note the dependency is soft**: `Vector` is
structural, so `TimedFullField` satisfies it TODAY (delegate dunders) even
before being a `Field`. So the two CAN be parallel; the recommended order is
P3.6-first only to avoid touching the conformance surface twice.

**Caveat from `feedback_unify_after_two_instances`:** the codebase docstring
(timed_full_field.py:48) defers `DirectSumSpace` to P3.6 "until kinetics'
flux⊕precursors ships the second use case." The #226 campaign IS arguably that
second instance (moment state = bulk⊕boundary moment is structurally the same
direct sum as flux⊕precursors). Confirm with the user whether #226 triggers
the two-instance rule or whether kinetics must land first. If #226 counts as
instance #2, P3.6 unblocks; if not, the campaign does Q1 with the structural
`Vector` Protocol ALONE (no `Field` promotion) and defers the `FullFieldSpace`
install — still correct, just leaves `TimedFullField` as a `Field`-LIKE
non-subclass a while longer.

### Elegance payoff
- **Structure-exposing**: `FullFieldSpace` stops being a write-only adjoint
  artefact and becomes the container's actual space — the half-built keystone
  gets its second half.
- **Structurally-simpler**: the container's hand-rolled dunders become
  inherited `Field` dunders (delete the delegate dunder bodies); the
  `_check_partner` becomes the inherited class+space gate.

### First test
Assert `isinstance(tff, Field)` after P3.6, `tff.space is a FullFieldSpace`,
`tff.l2 == sqrt(tff.space.inner_product(tff, tff))` (the per-block G-metric
norm), and that `Aψ.H` reciprocity (`<Aψ,φ>_G == <ψ,A^†φ>_G`) still holds to
1e-13 (the existing `test_g_adjoint_reciprocity` gate). Discriminates "P3.6
made the container a real Field with the right metric" from "container is still
a Field-like impostor."

---

## Q6 — scipy as typed serialization

### The adapter is ALREADY correct — formalise + single-source it.

Today TWO ravel paths exist: `as_scipy_linop` (bare ndarray) and
`KrylovAcceleration`'s inline `_ravel`/`_unravel_like` (typed via
`to_flat`/`from_flat`). The design: **operators are `Vector`-typed; the flat
ndarray appears ONLY inside the scipy boundary, via the field's
`to_flat`/`from_flat`.** This is what `KrylovAcceleration` already does. The
campaign formalises it as ONE adapter that works for both:

```python
def as_scipy_linop_typed(op, template: Vector, dtype=float) -> spla.LinearOperator:
    """Wrap a Vector-typed operator for scipy Krylov. ndarray lives ONLY here.

    template provides shape/structure for round-trip (to_flat/from_flat for a
    container; ravel/reshape for a bare ndarray — the _is_ravellable fork).
    """
    flat0 = _ravel(template); n = flat0.size
    def matvec(x_flat):
        x = _unravel_like(template, x_flat)   # flat → Vector
        return _ravel(op.apply(x))            # Vector → flat
    rmatvec = (lambda y: _ravel(op.apply_transpose(_unravel_like(template, y)))
               if _has(op, CAP_APPLY_TRANSPOSE) else None)
    return spla.LinearOperator((n, n), matvec=matvec, rmatvec=rmatvec, dtype=dtype)
```

GMRES/BiCGSTAB still see a plain `(n,n)` `matvec: ndarray->ndarray` — **the
serialization is INSIDE `matvec`, invisible to scipy.** The `(n,n)` shape is
`(to_flat().size, to_flat().size)`. `template` is the SN flux composite
(B.5.2: solution lives in flux space, templated on `initial_guess`). This
MERGES `as_scipy_linop` (template = bare ndarray → `ravel`/`reshape` fork) and
`KrylovAcceleration`'s inline closure into one adapter — single source of
truth (Cardinal Rule 2), retiring the duplicated ravel logic.

### Key correctness invariant: `to_flat`/`from_flat` round-trip exactness.
`TimedFullField.to_flat` = `concat(bulk.values.ravel(), boundary.values)` and
`from_flat` reshapes `[:n_bulk]` and slices `[n_bulk:]`. This is the
load-bearing invariant (the docstring already names it). The campaign must
gate `from_flat(to_flat(x), x) == x` bit-exactly for flux AND moment
containers — the latter is NEW (moment `to_flat` must ravel the `(L+1,2L+1,...)`
bulk and the `BoundaryMomentField` trace consistently).

### Elegance payoff
- **Structure-exposing**: names scipy as a SERIALIZATION boundary (the user's
  framing) — `ndarray` is the wire format, `Vector` is the algebra. The
  `to_flat`/`from_flat` pair IS the serialize/deserialize.
- **Structurally-simpler**: ONE adapter replaces two ravel paths (Smell-16:
  two paths to the same serialization, currently bridged by hand).
- **Algorithmic-advantage**: none (same matvec count, same flat size) — and
  must be bit-identical (Cardinal Rule 1).

### First test
`as_scipy_linop_typed(L, flux_template)` and the legacy `KrylovAcceleration`
path must produce bit-identical GMRES iterates on a 2-group SN slab
(`info==0`, same `solution`, same `residual_history`). Plus `from_flat(to_flat
(x),x)==x` to 0 ULP. Discriminates "merged adapter is bit-identical" from a
silent ravel-order divergence.

---

## Q7 — Cross-method generality + sequencing + risk

### Cross-method generality (CP scalar, MoC ray)
The container is ALREADY cross-method (the whole point). CP: `bulk=ScalarFlux`,
`boundary=`region-interface current; MoC: `bulk=RayField`, `boundary=`
ray-endpoint phases (docstring lines 30–39). The `Vector`-typed Protocol +
`apply(x:V)->V` means the CP/MoC operators get the SAME algebra surface for
free — `(L−S−F)φ=q` reads identically whether φ is a flux-`TimedFullField`, a
scalar-`TimedFullField`, or a moment-`TimedFullField`. **This is the deepest
payoff: the cross-method goal is ACHIEVED by the container being the `Vector`,
not by per-method generics.** The transport-resolvent backbone (promoted
kernel) predicts this: SN/MoC/CP `solve` are three quadratures of ONE Peierls
resolvent `(Ω·∇+Σ_t)⁻¹`; the container-typed algebra is the code-level shadow
of that single resolvent. Diffusion is the exception (elliptic, not a
quadrature) — but it ALSO fits `(L−S−F)φ=q` with `bulk=ScalarFlux`, so the
container algebra covers it; only `L.solve` differs (BiCGSTAB vs sweep).

### What is bit-identical (pure typing) vs behavioral
- **Bit-identical (pure typing):** Q1 Protocol retype (annotation only);
  Q3 leaving the seven primitives flat; Q6 merged adapter (if ravel order
  preserved); P3.6-(5b) structured-storage promotion (dunders already match).
  These touch types/annotations, NOT arithmetic. Cardinal Rule 1 satisfied by
  golden-output gates.
- **Behavioral (needs verification beyond bit-identity):** Q4
  `BoundaryMomentField` + closing the moment state (changes `solve_moments`
  return from `(moment_buf, None)` to a closed container — a STATE-SHAPE change,
  verify the moment tensor is byte-identical and the new boundary block is
  zero where the old None implied zero); Q2 functional separation (moving keff/
  production estimators out of any operator mould — verify keff unchanged).
- **Genuinely new:** the trace ι* / boundary `B` extraction (Wave O O.4a.2) is
  NOT part of #226 — it is a separate behavioral change. #226 should type the
  EXISTING algebra, not pull `B` out of the sweep. Keep them decoupled.

### Safe order of operations
1. **Q1 `Vector` Protocol in numerics** (new file, no consumer change yet) —
   pure addition, zero risk.
2. **P3.6-(5b): install `FullFieldSpace`/`DirectSumSpace` as
   `TimedFullField.space`; promote to `Field`** (if the two-instance rule
   permits — confirm with user). Bit-identical; the dunders already match.
3. **Retype the Protocol `apply(x:V)->V`** (annotation flip). Bit-identical.
   Now `TimedFullField` satisfies it through `Field`-ness.
4. **Q6 merged scipy adapter** — single-source the two ravel paths.
   Bit-identical gate.
5. **Q2 functional category** — extract keff/production into a `Functional`
   Protocol. Bit-identical (same reductions, new type home).
6. **Q4 `BoundaryMomentField` + moment-state closure** — behavioral; the only
   step that changes a return shape. Do LAST, gated hard. (This is the one
   place to dispatch test-architect proactively — it crosses scalar↔angular↔
   moment representation, the Smell-16 hazard.)

Steps 1–5 are bit-identical typing; step 6 is the single behavioral change.
This ordering front-loads zero-risk typing and isolates the one behavioral
step.

### Verification strategy (Cardinal Rule 1 — the algebra stays numerically
identical)
- **Golden-output gate**: the full SN regression snapshot (the existing
  `SAFETY×conv_tol` keff snapshots, the MMS L1 gates) must stay GREEN through
  steps 1–5 with ZERO ULP change (these steps touch no arithmetic).
- **Round-trip gate**: `from_flat(to_flat(x),x)==x` to 0 ULP for flux AND
  moment containers (Q6).
- **Adjoint reciprocity gate**: `test_g_adjoint_reciprocity` (`<Aψ,φ>_G ==
  <ψ,A^†φ>_G` ≤1e-13) survives P3.6 (the metric must still be the per-block G).
- **Type-gate (negative tests)**: `tff_flux + tff_moment` RAISES;
  `op.apply` returns the SAME container type; the seven primitives reject a
  `TimedFullField` (they want a leaf `.values`); a `Functional` is NOT a
  `LinearOperator`.
- **Step 6 only**: moment tensor byte-identical pre/post, new boundary block
  provably zero-equivalent to the old `None`.

### Elegance payoff (the whole campaign)
- **Structure-exposing**: the algebra's vector type stops being a documented
  lie; the cross-method resolvent (`(Ω·∇+Σ_t)⁻¹`) gets a single code-level
  carrier.
- **Structurally-simpler**: ONE Protocol vector type, ONE container, ONE scipy
  adapter, ONE space (`DirectSumSpace`); the per-method/per-rep class grid is
  refused.
- **Expressive**: `(L+C−S−F−B)ψ=q` reads as the math for SN, CP, MoC,
  diffusion, flux-state and moment-state alike.
- **Algorithmic-advantage**: none by design — and that is the success
  criterion (Cardinal Rule 1: a typing/architecture refactor must be
  numerically inert).

---

## Smell-16 sightings in this design (promote-candidates)

- **Q4 missing `BoundaryMomentField`**: the moment solve has ONE physical
  boundary quantity stored in (or dropped from) an angular `BoundaryFlux`
  representation — Smell-16 shape 2 (one quantity, two incompatible reps
  bridged by hand / by `None`). FIX = name the trace-moment type.
- **Q6 two ravel paths**: `as_scipy_linop` + `KrylovAcceleration`'s inline
  closure are two code paths to ONE serialization — Smell-16 shape 1. FIX =
  one `as_scipy_linop_typed`.
- **Q1 ndarray annotation on container-consuming operators**: the Protocol
  signature and the runtime type are two representations of "what apply
  consumes" that disagree — a documentation-level Smell-16. FIX = the `Vector`
  bound makes them agree.

## Cross-frame note (transport-resolvent backbone)
The endomorphism-algebra confirmation (Q2) is the code shadow of the promoted
kernel: SN/MoC/CP `solve` = three quadratures of `(Ω·∇+Σ_t)⁻¹`. A single
container-typed resolvent algebra is exactly what "one resolvent, three
quadratures" looks like in types. Diffusion's elliptic exception fits the
container (`bulk=ScalarFlux`) with only `L.solve` differing — confirming the
container is the right cross-method `Vector`.
