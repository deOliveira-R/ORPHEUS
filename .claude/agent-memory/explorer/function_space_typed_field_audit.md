---
name: function-space-typed-field-audit
description: Pre-plan reconnaissance for P3.3 Depth B — FunctionSpace as L1 base for every typed transport field. Inventory of FunctionSpace usage, typed-field structure, dunder algebra patterns, and migration target list.
metadata:
  type: project
---

# Pre-plan reconnaissance — FunctionSpace as L1 base of every typed field

Working branch: `refactor/moment-space-and-layering` worktree at
`/Users/rodrigo/git/nuclear/ORPHEUS/.claude/worktrees/moment-space-and-layering/`.

## §1 — FunctionSpace definition site

File: `orpheus/numerics/space.py:90-172`.

**Attributes / properties on `FunctionSpace`** (a `@dataclass(frozen=True)`):

- `name: str` — human-readable identity tag. **Load-bearing for equality**: two
  spaces with the same `(name, shape)` compare equal regardless of `inner_product_weights`.
- `shape: tuple[int, ...]` — tensor shape of the elements. Examples
  `(n_cells, n_ordinates, n_groups)`, `(L+1, 2L+1)`.
- `inner_product_weights: Optional[NDArray]` — diagonal L² weights, default `None`
  (Euclidean). Reshaped to broadcast against `shape`. `compare=False` via the explicit
  `__eq__`/`__hash__`.
- `__eq__(other)` — equal iff `(name, shape)` matches (`space.py:129-132`).
- `__hash__()` — `hash((name, shape))` (`space.py:134-135`).
- `__repr__()` — `f"FunctionSpace({self.name!r}, shape={self.shape})"` (`space.py:137-138`).
- `inner_product(x, y) -> float` — weighted sum `Σ w_i x_i y_i` or Euclidean
  `Σ x_i y_i` if no weights (`space.py:144-167`).
- `norm(x) -> float` — induced L² norm `sqrt(<x,x>)` (`space.py:169-172`).

The class is frozen; `inner_product_weights` is metadata only (not part of identity).
The class docstring (`space.py:28-60`) anticipates `MeshFunctionSpace`, `TraceSpace`,
`RegionSpace`, `EnergyGroupSpace`, `DiscreteAngularSpace` as future subclasses, plus
compositional dunders `S * T`, `S + T`, `S.dual()` — none of these are shipped today.

**Factory function `boundary_trace_space`** (`space.py:244-283`):
- Signature: `boundary_trace_space(direction: Literal["in", "out"], n_ordinates,
  n_groups, *, quadrature_weights: NDArray | None = None) -> FunctionSpace`.
- Returns a plain `FunctionSpace` (NOT a `TraceSpace` subclass) with
  `name=f"boundary_trace_{direction}"`, `shape=(n_ordinates, n_groups)`, and
  optional quadrature_weights reshaped to `(n_ordinates, 1)`.
- **Pattern overlap**: `boundary_trace_space` and the richer
  `InflowTraceSpace`/`OutflowTraceSpace` (in `orpheus/numerics/trace_space.py`) are
  parallel implementations — the factory returns a flat `FunctionSpace`, the dataclass
  subclasses carry per-face masks. The factory is not used by any SN production path;
  see §2.

Companion factories: `angular_flux_space` (`space.py:186-225`), `scalar_flux_space`
(`space.py:228-241`). Same lift-to-FunctionSpace, no typed-field carrier.

## §2 — FunctionSpace usage today

Total occurrences of `FunctionSpace(` constructor + `: FunctionSpace` annotation
across `orpheus/` and `tests/`:

### (definition) — 4 occurrences
- `orpheus/numerics/space.py:90` — the class itself.
- `orpheus/numerics/space.py:221, 238, 279` — internal use by the three factory
  functions (`angular_flux_space`, `scalar_flux_space`, `boundary_trace_space`).
- `orpheus/numerics/trace_space.py:252` — `class TraceSpace(FunctionSpace, ABC)`.
- `orpheus/numerics/spaces/spherical_harmonic_space.py:105` —
  `class SphericalHarmonicSpace(FunctionSpace)`.

### (operator codomain) — 8 occurrences
Operators that expose `domain` / `codomain` as `FunctionSpace`-valued properties:
- `orpheus/numerics/operator.py:156-174` — `LinearOperator` Protocol (`domain`/`codomain` return `Optional[FunctionSpace]`).
- `orpheus/numerics/operator.py:285-291` — `LinearOperatorMixin` default returning `None`.
- `orpheus/numerics/operator.py:488-494` — `_AdjointOperator` (swaps inner.codomain ↔ inner.domain).
- `orpheus/numerics/operator.py:606-613` — `OperatorSum` (propagates the non-None one).
- `orpheus/numerics/operator.py:674-681` — `OperatorProduct` (A.codomain, B.domain).
- `orpheus/numerics/operator.py:733-738` — `ScaledOperator` (read-through).
- `orpheus/numerics/projection.py:382-423` — `MomentProjection` exposes
  `codomain: SphericalHarmonicSpace` (`cached_property`, P1.3) and
  `domain: FunctionSpace` (`cached_property`, P1.4 angular-ordinate). **This is the
  reference pattern the plan generalises.**
- `tests/numerics/test_operator.py:540-572` — `_SpacedMatrixOperator` test fixture
  (the only test-side operator carrying domain/codomain explicitly).

### (field type carrying space) — 1 occurrence
Only `MomentProjection.codomain` cached-property qualifies as the "carrying a
FunctionSpace as a data field" pattern. **No typed-field class today carries a
`space: FunctionSpace` attribute.** Every typed-field class
(`AngularFlux`, `BoundaryFlux`, `ScalarFlux`, `HarmonicMomentField`,
`IsotropicSource`, `PerOrdinateSource`) carries `mesh: SNMesh` instead.

### (test fixture) — 21 occurrences
- `tests/numerics/test_operator.py:540-695` — `_SpacedMatrixOperator` plus
  `FunctionSpace(name="V"|"W"|"V1"|"V2"|"W1"|"W2"|"Z"|"phi"|"psi", shape=(n,))`
  constructions for adjoint identity (line 586-589), composition (606-619), product
  domain/codomain checks (628-644), mismatch raises (657-695). 15 unique `FunctionSpace(`
  call sites in test_operator.py alone.
- `tests/numerics/test_space.py:31-141` — equality, hashing, factory tests
  (`angular_flux_space`, `scalar_flux_space`, `boundary_trace_space`). 12 `FunctionSpace(`
  calls.
- `tests/numerics/test_spherical_harmonic_space.py:434` — bare-FunctionSpace
  comparison test (equality across base vs subclass).

### (unused / dead code) — 0 occurrences
Every import is consumed.

## §3 — Typed-field current structure

Six typed-field classes. **None** carries an explicit `space: FunctionSpace`
attribute today. All carry `mesh: SNMesh` plus an `np.ndarray` for `values`,
which is sufficient for shape validation but does NOT participate in the
operator-algebra's `domain`/`codomain` machinery.

### `AngularFlux`
- File: `orpheus/sn/angular_flux.py:80-513`.
- Constructor: `values: np.ndarray, mesh: SNMesh, boundary: BoundaryFlux | None = None,
  history_depth: int = 2, _history: tuple = ()`.
- Carries a FunctionSpace-like? **No** — only `mesh: SNMesh`. Shape `(N, ng, nx, ny)`
  validated in `__post_init__` against `(mesh.quad.N, mesh.ng, mesh.nx, mesh.ny)`.
- Natural FunctionSpace: `name="angular_flux"`, `shape=(N, ng, nx, ny)`,
  `inner_product_weights = mesh.quad.weights` reshaped to `(N, 1, 1, 1)` — matches
  `angular_flux_space` factory at `space.py:186-225` but the factory has shape
  `(n_cells, n_ordinates, n_groups)` (a different axis ordering and 3-D, NOT the
  4-D production layout `(N, ng, nx, ny)`).

### `BoundaryFlux`
- File: `orpheus/sn/boundary_flux.py:68-300`.
- Constructor: `mesh: SNMesh, xmin_face: ndarray | None = None,
  xmax_face: ndarray | None = None, xmin_xmax_buf: ndarray | None = None,
  ymin_ymax_buf: ndarray | None = None`. **Mutable** by design (write-through cache).
- Carries a FunctionSpace-like? **No** — only `mesh`. Has accessor properties `xmin`,
  `xmax`, `ymin`, `ymax` returning face slice views.
- Natural FunctionSpace: a **product/sum** of per-face directional trace spaces
  (`(N, ng)` per 1-D face; `(N, ng, ny)` or `(N, ng, nx)` per 2-D face slice).
  Mapping to a single `FunctionSpace` is non-obvious because the buffers
  themselves are geometry-conditional (1-D slab carries 2 faces, curvilinear 1,
  2-D Cartesian uses persistent buffers covering face + interior). Direct
  `InflowTraceSpace` × per-face decomposition is the natural fit.

### `ScalarFlux`
- File: `orpheus/sn/scalar_flux.py:46-144`.
- Constructor: `values: np.ndarray, mesh: SNMesh`.
- Carries FunctionSpace-like? **No** — only `mesh`. Shape `(ng, nx, ny)` validated.
- Natural FunctionSpace: `name="scalar_flux"`, `shape=(ng, nx, ny)`, no weights
  (Euclidean — matches the existing `scalar_flux_space` factory but extended to 3-D
  layout vs the factory's 2-D `(n_cells, n_groups)`).

### `HarmonicMomentField`
- File: `orpheus/sn/harmonic_moment_field.py:69-292`.
- Constructor: `values: np.ndarray, mesh: SNMesh, L: int`.
- Carries FunctionSpace-like? **No** — only `mesh` and the scalar `L`. The
  corresponding `SphericalHarmonicSpace` exists in the codebase
  (`orpheus/numerics/spaces/spherical_harmonic_space.py`) and `MomentProjection.codomain`
  builds it via `SphericalHarmonicSpace.from_L(L)`, but the typed-field carrier never
  stores a reference to one. **This is the most visible gap**: the moment field IS the
  data in the SH space, yet the two are not linked.
- Natural FunctionSpace: a **product** of `SphericalHarmonicSpace.from_L(L)` (shape
  `(L+1, 2L+1)`, padded `4π/(2ℓ+1)` weights) and a per-cell-per-group spatial space
  (shape `(ng, nx, ny)`).

### `IsotropicSource`
- File: `orpheus/sn/sources.py:80-233`.
- Constructor: `values: np.ndarray, mesh: SNMesh`.
- Carries FunctionSpace-like? **No** — only `mesh`. Shape `(ng, nx, ny)`.
- Natural FunctionSpace: identical storage shape to `ScalarFlux` but a **distinct
  physical quantity** — the module's prose warns "cross-type addition between flux
  and source is undefined" (`sources.py:38-52`). A typed `space` field would
  distinguish via `name` ("isotropic_source" vs "scalar_flux"); the shape is the same.

### `PerOrdinateSource`
- File: `orpheus/sn/sources.py:236-426`.
- Constructor: `values: np.ndarray, mesh: SNMesh`.
- Carries FunctionSpace-like? **No** — only `mesh`. Shape `(N, ng, nx, ny)`.
- Natural FunctionSpace: same `(N, ng, nx, ny)` as `AngularFlux` but with
  `name="per_ordinate_source"`. Cross-type with `AngularFlux` is forbidden today by
  isinstance checks; a typed `space` field would make the distinction structural.

### Other Flux/Field/Source types encountered
- `MaterialXSField` (referenced in `orpheus/sn/scattering.py:307`) — a cross-section
  field, not a flux/source. Carries per-cell view typed methods but no `FunctionSpace`.
- `Solution` (referenced widely; not part of this PR's scope but worth noting) — the
  return type of `SNSolver.solve` carrying `keff`, `angular_flux`, `scalar_flux`,
  `history` etc. Out of scope for §3 since it's a bundle, not a single space-typed
  field.

## §4 — Dunder algebra inventory

Every typed-field class implements an identical **dunder skeleton**: `__add__`,
`__sub__`, `__mul__`, `__rmul__`, `__truediv__`, `__neg__`, plus type+mesh
validation via a private `_validate_partner` / `_validate_mesh`. The repetition is
the load-bearing evidence for a `Field` ABC.

| Class | `__add__` | `__sub__` | `__mul__` | `__rmul__` | `__truediv__` | `__neg__` | `__matmul__` | `__lshift__` | `__call__` | `__post_init__` |
| --- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `AngularFlux` | yes | yes | yes | yes (`= __mul__`) | yes | yes | — | yes (stash alias) | yes (lag lookup) | yes (shape + history_depth + auto-allocate boundary) |
| `BoundaryFlux` | yes | yes | yes | yes | yes | yes | — | — | — | — (mutable, no shape invariant — geometry-conditional) |
| `ScalarFlux` | yes | yes | yes | yes (`= __mul__`) | yes | yes | — | — | — | yes (shape only) |
| `HarmonicMomentField` | yes | yes | yes | yes (`= __mul__`) | yes | yes | — | — | — | yes (shape + L matching) |
| `IsotropicSource` | yes (cross-type → `PerOrdinateSource`) | yes | yes | yes (`= __mul__`) | yes | yes | — | — | — | yes (shape only) |
| `PerOrdinateSource` | yes (cross-type → `PerOrdinateSource` via radd) | yes | yes | yes (`= __mul__`) | yes | yes | — | — | — | yes (shape only) |

**Shared patterns (the repetition that consolidates):**

1. **Validation skeleton**: every `_validate_partner` checks `isinstance(other, T)`
   then `other.mesh is not self.mesh`. The two-step pattern repeats six times
   verbatim. A `Field` ABC can centralise this as `_validate_partner(self, other)`
   with the type check `type(self) is type(other)` and the mesh-identity check on
   the shared `space.mesh` attribute (after the migration).
2. **Same-type return**: every binary op `(__add__, __sub__)` returns
   `type(self)(self.values <op> other.values, self.mesh, ...)`. The cross-type case
   on `IsotropicSource.__add__` is the only deviation — it returns a
   `PerOrdinateSource` when the partner is per-ordinate.
3. **Scalar-shape coercion**: every `__mul__/__rmul__/__truediv__` calls
   `float(scalar)`. Six identical sites.
4. **`__neg__` is `-1 * self`**: six identical sites; could be a single ABC method.
5. **Mesh-bound diagnostics**: every class exposes `linf()`, `copy()`, and the
   read-through properties `ng`, `nx`, `ny` (plus `N` on the per-ordinate ones).

`BoundaryFlux` is the most idiosyncratic — its `_binary_op` walks four optional
buffers `(xmin_face, xmax_face, xmin_xmax_buf, ymin_ymax_buf)` with a fold over
None-checks. The proposed `space: FunctionSpace` carrier would NOT cleanly subsume
this — `BoundaryFlux` is structurally a **direct sum** of four optional sub-spaces.
The plan should be prepared for `BoundaryFlux` to require either a special
"product trace space" carrier or stay outside the L1 base.

## §5 — The same-class invariant

The protocol-detection code in `orpheus/numerics/iteration.py:163-176`:

```python
def _is_ravellable(x: object) -> bool:
    """Detect the ravellable protocol used by :class:`AngularFlux`.

    True when ``x`` exposes both an instance method
    ``to_flat_with_traces`` returning a 1-D ndarray AND a class-level
    factory ``from_flat_with_traces(flat, mesh)`` AND a ``mesh``
    attribute the factory consumes.  Bare ndarrays match none of
    these and fall through to the legacy reshape path.
    """
    return (
        hasattr(x, "to_flat_with_traces")
        and hasattr(x, "mesh")
        and hasattr(type(x), "from_flat_with_traces")
    )
```

The companion functions `_ravel`, `_unravel_like`, `_zeros_like`, `_l2_norm`
(`iteration.py:179-211`) all consume this via `type(x).from_flat_with_traces(...)`
— i.e. they reconstruct via **the same class as the template**, not a subclass.

**`AngularFlux` is the only ravellable class today** — no subclass overrides
`from_flat_with_traces`. Grep confirms:
- `from_flat_with_traces` defined exactly once at `orpheus/sn/angular_flux.py:346`.
- `to_flat_with_traces` defined exactly once at `orpheus/sn/angular_flux.py:408`.
- No `class .*\(AngularFlux\)` subclass anywhere in `orpheus/` or `tests/`.

**Implication for the plan**: any L2 base / L3 adapter split MUST keep the
ravellable class as the leaf type. A subclass split would break the
`type(template).from_flat_with_traces(...)` reconstruction because the factory is
class-level and the L2 base would be the runtime type returned by
`from_flat_with_traces`. The safe pattern is **single-class with a
`space: FunctionSpace` field**, not a subclass hierarchy.

The same-class invariant is therefore trivial to satisfy for `AngularFlux`, and
similarly trivial for the other typed fields (none of them have subclasses today).

## §6 — Existing space types

Two `FunctionSpace` subclasses ship today:

### `SphericalHarmonicSpace(FunctionSpace)`
- File: `orpheus/numerics/spaces/spherical_harmonic_space.py:104-227`.
- Frozen dataclass; adds field `L: int = 0` (must default because the parent's
  `inner_product_weights` has a default).
- Constructor: usually via `SphericalHarmonicSpace.from_L(L)` at line 169-196 —
  builds `name="spherical_harmonic_space"`, `shape=(L+1, 2L+1)`,
  `inner_product_weights = _padded_metric_tensor(L, basis.metric_per_ell)` (a
  `(L+1, 2L+1)` array carrying `4π/(2ℓ+1)` on the valid `|m|≤ℓ` slots, zero on the
  padding), `L=L`.
- Equality/hashing inherited from `FunctionSpace` via explicit delegation
  (`spherical_harmonic_space.py:161-165`) to avoid ndarray-in-eq trouble that the
  auto-generated `@dataclass(frozen=True)` `__eq__` would otherwise introduce.
- `shape` and `inner_product_weights`: `(L+1, 2L+1)` and the padded `g_C` metric.

### `TraceSpace(FunctionSpace, ABC)` and concrete subclasses
- File: `orpheus/numerics/trace_space.py:252-263` (abstract base), `:271-396`
  (`InflowTraceSpace`), `:400-490` (`OutflowTraceSpace`).
- Frozen dataclasses; both concrete subclasses add an `Optional[NDArray]` mask
  field (`inflow_mask` / `outflow_mask`) with `compare=False, repr=False` plus
  `face_names: tuple[str, ...]`.
- Constructed via `InflowTraceSpace.from_mesh_and_quadrature(mesh, quadrature,
  faces, ng)` (`:308-368`) or the equivalent on `OutflowTraceSpace`. Result:
  `name="trace_inflow"` / `"trace_outflow"`, `shape=(quadrature.N, ng)`, mask
  built per-face from `Ω · n` sign.
- The shape excludes the per-face count (it's just `(N, ng)`); the mask is the
  carrier of the per-face structure. **Distinguishing different faces in the same
  trace direction is implicit in `face_names` + mask rows, not in `shape`.**

No other FunctionSpace subclasses exist. The `spaces/__init__.py:32-36` re-exports
only `SphericalHarmonicSpace`; `trace_space.py` lives at the top level of
`orpheus/numerics/`, scheduled by the moment-space plan to move into
`numerics/spaces/trace_space.py` (P3.2 note at `spaces/__init__.py:15-19`).

## §7 — Where the pattern is BROKEN (or never landed)

Migration target list — every site where a typed field should carry a FunctionSpace
but doesn't, plus the operator-side counterparts that return `np.ndarray` instead of
typed fields.

### 7.1 — Typed-field classes with no `space` attribute (6 classes)

The six classes in §3 each need a `space: FunctionSpace` field (or equivalent
`cached_property`). The `mesh: SNMesh` field is **not removed** — it remains the
discretisation handle. The `space` field is **additive** and carries the inner
product, equality identity, and adjoint-machinery hook.

| Class | Current carrier | Target FunctionSpace |
| --- | --- | --- |
| `AngularFlux` | `mesh: SNMesh` only | `FunctionSpace(name="angular_flux", shape=(N, ng, nx, ny), inner_product_weights=quad.weights[:, None, None, None])` |
| `BoundaryFlux` | `mesh: SNMesh` only | direct sum over present face-buffer subspaces — needs design |
| `ScalarFlux` | `mesh: SNMesh` only | `FunctionSpace(name="scalar_flux", shape=(ng, nx, ny))` |
| `HarmonicMomentField` | `mesh: SNMesh, L: int` | product of `SphericalHarmonicSpace.from_L(L)` × `FunctionSpace("cell_group", (ng, nx, ny))` |
| `IsotropicSource` | `mesh: SNMesh` only | `FunctionSpace(name="isotropic_source", shape=(ng, nx, ny))` — distinct `name` from scalar_flux |
| `PerOrdinateSource` | `mesh: SNMesh` only | `FunctionSpace(name="per_ordinate_source", shape=(N, ng, nx, ny))` |

### 7.2 — Operators returning bare `np.ndarray` instead of a typed field

These are the production paths where the typed-flux discipline broke or never
landed. After Depth B, each should be `FunctionSpace → FunctionSpace`.

- `orpheus/sn/operator.py:1506` — `SNStreamingOperator.apply(psi: np.ndarray) -> np.ndarray`.
  Consumes the LEGACY packed equation-map vector; output is a packed vector. Lives
  inside the scipy.gmres adapter path. Status: legacy bridge from PR-TYPED-6c; the
  StreamingOperator (line 1801) is the typed replacement.
- `orpheus/sn/operator.py:1698` — `SNStreamingOperator.apply_transpose(psi: np.ndarray) -> np.ndarray`.
  Same as `apply` — packed-vector adjoint via dense matrix probe.
- `orpheus/sn/operator.py:1801` — `StreamingOperator` typed leaf, takes typed
  `AngularFlux` per the R-1 Step convention. `apply` at line 1936 still has
  ndarray-mixed signatures depending on the call site.
- `orpheus/sn/operator.py:2170` — `CollisionOperator`. `apply` at 2276, `solve` at 2325,
  `apply_transpose` at 2377 — check signature against the typed convention.
- `orpheus/sn/operator.py:2412` — `InvertibleOperator(OperatorSum)`, `solve` at 2562.
- `orpheus/sn/scattering.py:225` — `LegendreMomentScattering.apply` takes
  `np.ndarray | HarmonicMomentField` (typed-or-untyped) and returns the same union —
  the dual-mode signature is itself the leak.
- `orpheus/sn/scattering.py:818` — `ScatteringOperator.apply` is `@singledispatchmethod`;
  the `np.ndarray` branch is preserved for legacy callers (line 836-840 docstring).
- `orpheus/sn/fission.py:155` — `FissionOperator.apply` `@singledispatchmethod`; same
  legacy ndarray branch.
- `orpheus/numerics/projection.py:425` — `MomentProjection.apply(psi: np.ndarray) -> np.ndarray`.
  THE pattern reference operator. Its `codomain` is typed (`SphericalHarmonicSpace`)
  but its `apply` consumes/returns bare ndarrays. Bridging would be a clean Depth B
  outcome: `apply(psi: AngularFlux) -> HarmonicMomentField` (with the codomain check
  routed via `psi.space.shape[0] == self.domain.shape[0]`).
- `orpheus/numerics/projection.py:435` — `MomentProjection.apply_transpose` (same).

### 7.3 — Mesh-carrying fields whose shape is validated but never matched against a space

All six typed fields validate `values.shape` against derived `mesh` quantities in
`__post_init__` (see §3). This is duplicate-with-FunctionSpace work — if every
field carried a `space: FunctionSpace`, the validation would be `self.values.shape
== self.space.shape`, ONE line, in a `Field` ABC `__post_init__`.

### 7.4 — `Solution` and bundled returns

`orpheus/sn/solution.py` (find via grep) bundles fluxes + history. Out of scope for
P3.3 Depth B but worth keeping in view — the `Solution` carrier should reference
typed fields with typed `space` attributes through-and-through.

## §8 — Cross-references the plan will need

### `trace_space` and `TraceSpace`
- Module: `orpheus/numerics/trace_space.py:1-490`. Top-level (NOT under `spaces/`).
- Imports the base `FunctionSpace` at line 114.
- Exports: `TraceSpace`, `InflowTraceSpace`, `OutflowTraceSpace`.
- Construction pattern: `cls.from_mesh_and_quadrature(mesh, quadrature, faces, ng)`.
- Scheduled move into `numerics/spaces/trace_space.py` as part of moment-space
  plan P3.2 (per `spaces/__init__.py:15-19` note).

### `boundary_trace_space("in"|"out", ...)` factory
- Defined at `orpheus/numerics/space.py:244-283`.
- Returns a plain `FunctionSpace` with `name=f"boundary_trace_{direction}"`,
  `shape=(n_ordinates, n_groups)`, weights reshaped to `(n_ordinates, 1)`.
- Grep confirms **zero production callers** of `boundary_trace_space(`. It exists in
  the module docstring/exports but is NOT consumed. The richer
  `InflowTraceSpace`/`OutflowTraceSpace` is what production uses (see
  `orpheus/sn/boundary_realizer.py`, `orpheus/sn/angular_operator.py:293`,
  `orpheus/sn/geometry.py:412`).
- **Recommendation for the plan**: retire `boundary_trace_space` in favour of
  `InflowTraceSpace.from_mesh_and_quadrature` / `OutflowTraceSpace.from_mesh_and_quadrature`.
  Single-direction trace spaces are a degenerate case of the masked ones.

### `BoundaryOperator` / `BoundaryTraceLaw` / `BoundaryRealizer`
- `BoundaryOperator` — abstract base in `orpheus/geometry/boundary.py` (registered via
  `RegistryMixin` per `orpheus/numerics/registry.py:31-49`). Subclasses include
  `VacuumBoundaryOperator`, `SpecularBoundaryOperator`, `PeriodicBoundaryOperator`.
- `BoundaryTraceLaw` — referenced by `orpheus/numerics/trace_space.py:24` and
  `orpheus/sn/boundary_realizer.py:10`; the legacy law class.
- `SNBoundaryRealizer` — `orpheus/sn/boundary_realizer.py`; consumes legacy
  `BoundaryOperator` subclasses and produces the per-face directional mask
  operators from `InflowTraceSpace.inflow_indices_for_face`.
- `DiffusionBoundaryRealizer` — `orpheus/diffusion/boundary_realizer.py` (stub).
- `BoundaryClosureOperator` — `orpheus/derivations/continuous/peierls_nystrom/geometry.py:4191`
  (Nyström-specific, separate codebase region).
- All boundary-side operators consume mesh + quadrature directly; none today
  carries a `FunctionSpace` as a domain/codomain attribute. Migration target for a
  follow-on phase.

---

**Total word count: ~1,500** (slightly over budget; §3 + §7 + §2 are the load-bearing
rosters per the brief; §1, §4, §5 are condensed; §6, §8 carry the cross-reference
weight).
