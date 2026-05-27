# Depth B — `Field` on `FunctionSpace` as the L1/L2 algebraic spine

**Branch:** `refactor/moment-space-and-layering`
**Worktree:** `.claude/worktrees/moment-space-and-layering/`
**Phase status:** P3.1, P3.5, P3.0, P3.2 LANDED. P3.3 originally framed as "introduce `transport/`; migrate fields"; this document REPLACES that framing with a deeper unification (Depth B). P3.4 (Problem/Solver split) and P3.6 (kinetics restructure) follow Depth B's outcome.

**Date:** 2026-05-27. Plan author: main agent + explorer audit memo at `.claude/agent-memory/explorer/function_space_typed_field_audit.md`. Reference: `.claude/plans/neutron_transport_grand_report_v3.md` §§5.3-5.5, §6.1, §16A, §32.4-32.6, §33.

**Status:** DRAFT, awaiting approval before implementation. Context will be compacted after approval. The plan is self-contained for fresh-session pickup.

**Detour relationship to the parent plan.** Depth B is a re-scoping of step **P3.3** of `.claude/plans/moment_space_and_layering_plan.md`. The parent plan's P3.3 framed the work as "introduce `transport/`; migrate fields." This Depth-B plan extends that framing into a foundational unification (`Field` on `FunctionSpace`) because the audit showed the migration target is wider than P3.3's original scope. When Depth B completes (step D-K), execution **returns to the parent plan at step P3.4** (Problem / Solver split). See **§11 — Exit Route** for the precise hand-off.

**Downstream consumers of Depth B.** Two open GitHub issues are direct beneficiaries:
- **Issue #201** — angular-slice dimensional split (`AngularFlux` / `AngularSource` / `AngularResidual`). The `units` machinery introduced here makes #201's three-type discipline machine-checkable.
- **Issue #205** — full cross-method (storage × role) field architecture (12-cell matrix: Angular/Scalar/Track/Region × Flux/Source/Residual). Depth B is the foundational layer that makes #205 implementable.

Neither issue is in scope FOR Depth B (per "unify after two instances"), but Depth B is the enabling phase. Both should be cross-referenced as the "what this enables" payoff.

---

## 0. Pickup checklist (read first)

If you are picking this plan up in a fresh session:

1. **Read this plan top-to-bottom.** No section is optional.
2. **Read the explorer audit memo:** `.claude/agent-memory/explorer/function_space_typed_field_audit.md`. It contains the FUNCTION-LEVEL inventory of every site that needs to change and is the source of truth for "what is where today."
3. **Read these grand-report sections only** (they suffice):
   - §2 The core ontology (lines 114-150)
   - §5.3 Space hierarchy (lines 521-550)
   - §5.5 Field hierarchy (lines 576-605)
   - §5.7 Operator hierarchy (lines 636-672)
   - §6.1 Space dunders (lines 680-695)
   - §16A.1 Boundary primitive hierarchy (lines 2738-2775)
   - §32.4-32.6 Spaces, Fields, Operators primitive specs (lines 5874-5982)
4. **Read these existing typed-field files** (smallest first):
   - `orpheus/sn/scalar_flux.py` (simplest field — Euclidean L²)
   - `orpheus/sn/harmonic_moment_field.py` (clearest gap — SH space exists, not referenced)
   - `orpheus/sn/angular_flux.py` (most complex — ravellable protocol, boundary auto-allocation)
   - `orpheus/sn/boundary_flux.py` (idiosyncratic — geometry-conditional buffers)
   - `orpheus/sn/sources.py` (source pair)
5. **Read `orpheus/numerics/space.py`** in full (172 lines) and `orpheus/numerics/spaces/spherical_harmonic_space.py` in full (227 lines).
6. **Read `orpheus/numerics/trace_space.py`** (490 lines, but only TraceSpace ABC + InflowTraceSpace + OutflowTraceSpace + the factory `.from_mesh_and_quadrature`).
7. **Validate the plan still makes sense** against the current branch state (`git log --oneline -10`).
8. **Pick up at the leftmost incomplete step in §6** below.

The `feedback_no_method_implementer_for_surgical_carves` rule applies: this is the main agent's work with turn-by-turn user steering. Do NOT batch via method-implementer.

---

## 1. The principle

**Every transport field is `(space: FunctionSpace, values: ndarray) + algebra`.**

Per grand-report §32.5, a `Field` carries a `space` reference whose `shape` defines the data layout and whose `inner_product_weights` define the L² metric. Field arithmetic is closed under same-space addition; cross-space arithmetic is forbidden (either raises or routes through an explicit projection/reconstruction operator).

Today, ZERO typed-field classes follow this pattern. Every one of `AngularFlux`, `ScalarFlux`, `HarmonicMomentField`, `BoundaryFlux`, `IsotropicSource`, `PerOrdinateSource` carries `mesh: SNMesh` and an `ndarray` for `values`, but no `space`. The `mesh` carries the discretisation; the `space` would carry the L² structure. The two are orthogonal axes: `mesh` is INPUT (geometry), `space` is L1 (mathematical structure).

The lift is structurally cheap (an additive `space: FunctionSpace` field) and algebraically transformative (operator algebra, adjoint machinery, projection composition all align on the same `space` axis already used by `MomentProjection.codomain`).

### 1.1 The reference pattern: `MomentProjection.codomain`

Per P1.3, `MomentProjection` exposes `codomain: SphericalHarmonicSpace` (a `FunctionSpace` subclass). Its `_AdjointOperator` reads `codomain.inner_product_weights` to compute the Hilbert adjoint correctly. P3.5 promoted `codomain` to the canonical attribute name on the entire `LinearOperator` Protocol. The pattern is established; Depth B extends it from operators to fields.

### 1.2 The grand-report ontology

```text
MathObject
├── Domain
├── Measure
├── Space              ← FunctionSpace + subclasses (already exist)
├── Basis              ← SphericalHarmonicBasis (already exists)
├── Field              ← THIS PLAN: lift to L1 ABC, derive every flux type
├── Operator
├── Kernel
├── Projection         ← MomentProjection (already exists)
├── Flow
├── Process
├── Functional
├── Problem            ← Phase 3.4
├── Solver             ← Phase 3.4
└── Representation     ← Future (dense/sparse/matrix-free/tt)
```

This plan is "land the missing L1/L2 piece of the ontology": `Field` as the algebraic base that EVERY typed transport field inherits.

---

## 2. Current state of the codebase (load-bearing facts)

These statements are the audit memo's findings, restated for the planner.

1. **`FunctionSpace`** (`orpheus/numerics/space.py:90-172`) is a frozen `@dataclass` with `(name, shape, inner_product_weights)`. Identity is `(name, shape)`; weights are metadata (`compare=False`). It has `inner_product(x, y)` and `norm(x)` methods. **No dunders for Space algebra** (`*`, `+`, `.dual()`) — these need adding (grand-report §6.1).

2. **`SphericalHarmonicSpace(FunctionSpace)`** (`orpheus/numerics/spaces/spherical_harmonic_space.py:104`) — frozen, adds `L: int`. Built via `SphericalHarmonicSpace.from_L(L)` with `(L+1, 2L+1)` shape and padded `4π/(2ℓ+1)` metric.

3. **`TraceSpace(FunctionSpace, ABC)`** (`orpheus/numerics/trace_space.py:252-263`), with concrete `InflowTraceSpace` (`:271-396`) and `OutflowTraceSpace` (`:400-490`). Constructed via `cls.from_mesh_and_quadrature(mesh, quadrature, faces, ng)`. Shape is `(quadrature.N, ng)` plus a `inflow_mask`/`outflow_mask` per face. **Already exists, not yet consumed by `BoundaryFlux`.**

4. **`MomentProjection.codomain`** is the only place where an operator/field carries a `FunctionSpace` as a data field. **The reference pattern.**

5. **Six typed fields**, all with identical `__add__/__sub__/__mul__/__rmul__/__truediv__/__neg__` skeleton plus `_validate_partner` (isinstance + `mesh is mesh`). The repetition IS the load-bearing evidence for a `Field` ABC.

6. **The same-class invariant** (`orpheus/numerics/iteration.py:163-176`): `_is_ravellable` uses `hasattr(type(x), "from_flat_with_traces")`. Reconstruction via `type(template).from_flat_with_traces(...)`. `AngularFlux` has NO subclasses anywhere. Any L1/L2 split must keep `AngularFlux` as a SINGLE class with `space` as an additive field, NOT a subclass hierarchy.

7. **Dead factories** in `space.py`: `angular_flux_space()`, `scalar_flux_space()`, `boundary_trace_space()` — zero production callers. Retirement candidates.

8. **The boundary tower** (`orpheus/geometry/boundary.py`, `orpheus/sn/boundary_realizer.py`, `orpheus/diffusion/boundary_realizer.py`) already has `BoundaryOperator` + `BoundaryTraceLaw` + `BoundaryRealizer`. **Already partially aligned with grand-report §16A.** Out-of-scope direct consumers (e.g., `boundary_trace_space()` factory) can be retired.

9. **Operators returning bare ndarray** (`SNStreamingOperator.apply`, `LegendreMomentScattering.apply`, `MomentProjection.apply` — the LAST is the most embarrassing since its `codomain` is already typed) are the migration target for Depth-B's "operator-side completion" phase (Step D-G below).

---

## 3. The target architecture

### 3.1 Layer assignments

```text
L1 (mathematics, knows no neutrons)
├── orpheus/numerics/space.py
│   ├── FunctionSpace (existing)
│   ├── TensorProductSpace (NEW — V = X ⊗ Ω ⊗ G)
│   └── DirectSumSpace (NEW — V_total = FluxSpace ⊕ Precursors; also BoundaryFlux faces)
├── orpheus/numerics/spaces/
│   ├── spherical_harmonic_space.py (existing)
│   └── trace_space.py (move from numerics/ per P3.2 — InflowTraceSpace, OutflowTraceSpace, the ABC)
├── orpheus/numerics/field.py (NEW)
│   └── Field (ABC) — the L1 base of every typed transport field
└── orpheus/numerics/source.py (NEW, optional)
    └── Source (ABC) — same skeleton as Field, distinguished only by physical kind

L2 (transport vocabulary, method-agnostic)
└── orpheus/transport/
    ├── __init__.py
    ├── fields/
    │   ├── angular_flux.py — AngularFlux(Field)
    │   ├── scalar_flux.py — ScalarFlux(Field)
    │   ├── harmonic_moment_field.py — HarmonicMomentField(Field)
    │   └── boundary_flux.py — BoundaryFlux(Field) (with DirectSumSpace OR product-trace)
    └── sources/
        ├── isotropic_source.py — IsotropicSource(Source)
        └── per_ordinate_source.py — PerOrdinateSource(Source)

L3 (method-specific machinery)
└── orpheus/sn/
    ├── angular_flux_b1pp_adapter.py (NEW) — installs from_flat_with_traces /
    │   to_flat_with_traces as class methods on the L2 AngularFlux symbol
    │   (module-level injection at import; preserves same-class invariant)
    ├── boundary_flux_zeros_factory.py (NEW) — SN-specific BoundaryFlux.zeros_for_sn_mesh
    │   factory (separates the geometry-conditional construction from the dataclass)
    └── ... (existing solver, operator, sweep, etc.)
```

Why `Field` lives at L1 (`numerics/`) and `AngularFlux` lives at L2 (`transport/`):

- `Field` knows nothing about transport — it's just "values + space + algebra." That's L1.
- `AngularFlux` is a method-agnostic transport concept (it's a function on `(spatial × angular × group)` phase space, regardless of which method discretises it). That's L2.

### 3.2 The `Field` ABC and `FunctionSpace` units

**`FunctionSpace` gains a `units` attribute** so the space carries the dimensional identity of every Field that lives on it. Per the conversation in §3.0 (history), units belong on the space because the L² metric is unit-aware: two spaces with same name+shape but different units are structurally different mathematical objects.

```python
# orpheus/numerics/space.py — extended
import pint

@dataclass(frozen=True)
class FunctionSpace:
    name: str
    shape: tuple[int, ...]
    inner_product_weights: NDArray | None = None
    units: pint.Unit | None = None    # NEW — pint.Unit label; ALWAYS imported (no lazy)

    def __eq__(self, other):
        if not isinstance(other, FunctionSpace):
            return NotImplemented
        if (self.name, self.shape) != (other.name, other.shape):
            return False
        # If either side is unitless, treat as compatible. If both set,
        # compare DIMENSIONALITY (more robust than string equality).
        if self.units is None or other.units is None:
            return True
        return self.units.dimensionality == other.units.dimensionality

    def __hash__(self):
        return hash((self.name, self.shape))
```

`pint` is a real production dependency (`pint >= 0.20` in `pyproject.toml`). Import-time cost is ~100ms, irrelevant on simulation timescales. Lazy imports were considered and rejected — they would force `units: str` and add string-parsing overhead at every dimensional check.

```python
# orpheus/numerics/field.py
from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass, field as dc_field
from typing import ClassVar, TypeVar

import numpy as np
from numpy.typing import NDArray

from orpheus.numerics.space import FunctionSpace

T = TypeVar("T", bound="Field")


@dataclass(frozen=True, eq=False)  # eq=False because ndarray comparison is ambiguous
class Field(ABC):
    """L1 algebraic base of every typed transport field.

    A Field is `(space: FunctionSpace, values: ndarray)` with closed
    same-CLASS, same-SPACE arithmetic. Cross-class arithmetic is forbidden
    by design — even when units match. See `coding-elegance` Pattern 4
    (illegal states unrepresentable): same units gives PERMISSION (in
    linear algebra) to add; it does not give MEANING. Cross-class
    same-units operations require an explicit NAMED composition (e.g.
    `IterationResidual.from_balance(lhs, rhs)`) that constructs the
    result with a definite physical interpretation.

    See grand-report §32.5 for the design rationale and §5.5 for the
    hierarchy of derived field types.
    """
    values: NDArray
    space: FunctionSpace

    def __post_init__(self) -> None:
        if self.values.shape != self.space.shape:
            raise ValueError(
                f"{type(self).__name__}: values.shape {self.values.shape!r} != "
                f"space.shape {self.space.shape!r}"
            )

    # ── Algebra ─────────────────────────────────────────────────────────
    def _check_partner(self, other: object) -> None:
        # Layer 1 — class identity is the primary gate. Strictest rule.
        # Always runs (both debug and -O).
        if type(self) is not type(other):
            raise TypeError(
                f"{type(self).__name__} arithmetic requires same-class partner; "
                f"got {type(other).__name__}. Cross-class arithmetic (even "
                f"same-units) requires an explicit named composition."
            )
        if self.space != other.space:  # FunctionSpace.__eq__ matches (name, shape, units.dimensionality)
            raise ValueError(
                f"{type(self).__name__} arithmetic requires equal space; "
                f"got {self.space!r} vs {other.space!r}"
            )
        # Layer 3 — defensive units assertion, stripped in -O mode.
        # Tautological IF class hierarchy is well-designed; catches
        # class/units misdesign during development.
        assert self.space.units == other.space.units, (
            f"internal: class identity matched but units mismatch "
            f"({self.space.units} vs {other.space.units}) — "
            f"class/units inconsistency in {type(self).__name__}"
        )

    def __add__(self: T, other: T) -> T:
        self._check_partner(other)
        return type(self)(values=self.values + other.values, space=self.space)

    def __sub__(self: T, other: T) -> T:
        self._check_partner(other)
        return type(self)(values=self.values - other.values, space=self.space)

    def __neg__(self: T) -> T:
        return type(self)(values=-self.values, space=self.space)

    def __mul__(self: T, scalar: float) -> T:
        return type(self)(values=self.values * float(scalar), space=self.space)

    def __rmul__(self: T, scalar: float) -> T:
        return self.__mul__(scalar)

    def __truediv__(self: T, scalar: float) -> T:
        return type(self)(values=self.values / float(scalar), space=self.space)

    # ── Diagnostics ─────────────────────────────────────────────────────
    @property
    def linf(self) -> float:
        return float(np.abs(self.values).max())

    @property
    def l2(self) -> float:
        return float(self.space.norm(self.values))

    def inner_product(self: T, other: T) -> float:
        self._check_partner(other)
        return self.space.inner_product(self.values, other.values)

    def copy(self: T) -> T:
        return type(self)(values=self.values.copy(), space=self.space)
```

Notes:

- The `_check_partner` enforces `type(self) is type(other)` — `AngularFlux + ScalarFlux` raises by construction. This is the grand-report §32.5 invariant "fields should not silently combine if their spaces are incompatible."
- `space: FunctionSpace` is part of the dataclass — and `FunctionSpace.__eq__` is `(name, shape)`. So two AngularFlux instances on different meshes (but same shape, same name) WOULD compare equal at the space level — which is correct, because the L1 algebra doesn't care about the mesh; the L2 mesh-binding is added by subclasses.
- The `eq=False` is because comparing two `ndarray` values is ambiguous and we don't want auto-generated `Field.__eq__` calling `==` on the values arrays.

### 3.3 The L2 field types

Each subclass adds DOMAIN-SPECIFIC fields (e.g. `mesh`, `boundary`, `L`) on top of `Field`'s `(values, space)`. The dunders are inherited unchanged. Only the constructor + class-specific methods differ.

```python
# orpheus/transport/fields/angular_flux.py
@dataclass(frozen=True, eq=False)
class AngularFlux(Field):
    """L2 typed flux: a Field on the (N, ng, nx, ny) phase space.

    Carries a `space: FunctionSpace` for L² algebra, a `mesh: SNMesh`
    for discretisation handle (today method-specific; future Protocol),
    a `boundary: BoundaryFlux` partner field on the trace space, and
    an iteration history `_history` for Krylov tracking.
    """
    mesh: SNMesh  # added field
    boundary: BoundaryFlux  # added field, REQUIRED (no default — see §3.5)
    history_depth: int = 2
    _history: tuple = dc_field(default=(), compare=False)

    # __post_init__ extends Field.__post_init__ via super() — adds mesh-vs-space
    # consistency check; cannot reassign boundary auto-allocation default.

    # Class-level factories from_flat_with_traces / to_flat_with_traces
    # are NOT defined here. They are injected by the SN adapter at L3:
    # orpheus/sn/angular_flux_b1pp_adapter.py monkey-patches them onto
    # this class at module-import time. The same-class invariant
    # (numerics/iteration.py:163-176) is satisfied because type(x) is
    # always this L2 class.
```

Subtle: `dataclass` inheritance with `@dataclass(frozen=True)` requires every parent field to have a default OR every child field to come AFTER no-default parent fields. We use `frozen=True` already, so the order is constrained. Audit this when implementing — may need to declare `values` and `space` with defaults of `None` and validate, OR use `kw_only=True` (Python 3.10+; the project is on 3.14 per the venv inspection earlier).

### 3.4 `BoundaryFlux` — the geometry-conditional case (Option Ω, flat-buffer storage)

`BoundaryFlux` is structurally a **collection of `BoundaryFaceFlux` fields**, one per geometry boundary face. Each `BoundaryFaceFlux` is a Field on a per-face trace space.

**Naming convention:**
- `FaceFlux` (FUTURE, potential L2) — flux on ANY face (internal or boundary). Not in scope for this phase.
- `BoundaryFaceFlux` — flux on ONE boundary face's inflow trace space. This phase introduces it.
- `BoundaryFlux` — the BUNDLE of `BoundaryFaceFlux` over all boundary faces of a mesh.

**Storage: Option Ω — flat backing buffer + per-face slice views.**

A dict-of-fields would not vectorize: every `+`, `*`, `<·,·>` would iterate over faces with Python dispatch overhead per face. For SN's 2-4 face count this is ~1% of matvec cost (tolerable), but for MoC where each track-end is conceptually a "boundary face" (thousands of faces), the overhead is prohibitive.

Option Ω stores a single contiguous `_values: ndarray` (all face buffers concatenated), with a `FaceLayout` descriptor mapping `face_name → slice + reshape`. Arithmetic operates on the flat buffer in ONE numpy call. Per-face access goes through slice views — no copies, same memory:

```python
# orpheus/numerics/face_layout.py
@dataclass(frozen=True)
class FaceLayout:
    """Descriptor mapping face names to slices + per-face shapes in a flat backing buffer.

    Per-method mesh provides its own FaceLayout (SN: per-Cartesian-face; MoC:
    per-track-family; CP: not used). The layout's total_size determines the
    flat buffer's length; per-face shapes determine the slice-view shape.
    """
    faces: dict[str, FaceSlot]      # face_name -> FaceSlot
    total_size: int                  # sum of face flat sizes

@dataclass(frozen=True)
class FaceSlot:
    name: str
    offset: int                      # start index in flat buffer
    shape: tuple[int, ...]           # per-face shape (e.g. (N, ng) for 1D, (N, ng, ny) for 2D x-face)
    flat_size: int                   # product(shape)

    def slice_view(self, flat_buf: NDArray) -> NDArray:
        """Return a SHAPED view into the flat buffer — no copy."""
        return flat_buf[self.offset : self.offset + self.flat_size].reshape(self.shape)
```

```python
# orpheus/transport/fields/boundary_face_flux.py
@dataclass(frozen=True, eq=False)
class BoundaryFaceFlux(Field):
    """A Field on ONE boundary face's inflow trace space.

    Lightweight wrapper — inherits Field algebra (dunders + diagnostics).
    The `space` is an `InflowTraceSpace` with the face's shape and
    quadrature-weight metric. Construction is via the `BoundaryFlux`
    factory; standalone construction is rare.
    """
    # Inherits values, space from Field. No additional fields.
```

```python
# orpheus/transport/fields/boundary_flux.py
@dataclass(frozen=True, eq=False)
class BoundaryFlux:
    """Bundle of per-face boundary fluxes. NOT a Field (it's a structured
    bundle, not a single space's element). Arithmetic is delegated to the
    flat backing buffer for vectorized performance.
    """
    _values: NDArray                  # flat backing buffer, length == layout.total_size
    layout: FaceLayout                # mesh-provided descriptor
    mesh: SNMesh                      # discretisation handle (future: TransportMesh Protocol)

    def __post_init__(self):
        if self._values.shape != (self.layout.total_size,):
            raise ValueError(
                f"BoundaryFlux: flat buffer shape {self._values.shape} "
                f"!= ({self.layout.total_size},)"
            )

    # ── Per-face access (slice views, no copies) ────────────────────────
    @property
    def faces(self) -> Mapping[str, BoundaryFaceFlux]:
        return {
            name: BoundaryFaceFlux(
                values=slot.slice_view(self._values),
                space=self._space_for_face(name),  # constructed from layout + mesh metadata
            )
            for name, slot in self.layout.faces.items()
        }

    # ── Vectorized arithmetic on the flat buffer ────────────────────────
    def _check_partner(self, other: "BoundaryFlux") -> None:
        if type(other) is not type(self):
            raise TypeError(...)
        if self.layout is not other.layout:
            # Layouts must match by identity (cheap) or by structural equality (deeper check).
            if self.layout != other.layout:
                raise ValueError("BoundaryFlux layout mismatch")

    def __add__(self, other: "BoundaryFlux") -> "BoundaryFlux":
        self._check_partner(other)
        return BoundaryFlux(
            _values=self._values + other._values,   # ONE numpy call
            layout=self.layout,
            mesh=self.mesh,
        )
    # __sub__, __mul__, __neg__, __truediv__ — same pattern

    @property
    def linf(self) -> float:
        return float(np.abs(self._values).max())
```

**The mesh contract** (cross-method uniform):

```python
class TransportMesh(Protocol):
    """Protocol every method's mesh must satisfy to participate in transport/."""
    ng: int

    @property
    def boundary_face_layout(self) -> FaceLayout: ...
```

- SN's `SNMesh.boundary_face_layout` returns `{"xmin": FaceSlot(...), "xmax": FaceSlot(...), ...}` per Cartesian face.
- MoC's `MoCMesh.boundary_face_layout` returns one `FaceSlot` per track family (thousands of slots, but flat-buffer arithmetic stays O(1) numpy call).
- Diffusion's `DiffusionMesh.boundary_face_layout` returns face-only slots (no angular dimension).

The MESH PROVIDES the layout; the `BoundaryFlux` consumes it without knowing the method-specific construction logic. Cross-method generality is automatic.

**Why this is better than dict-of-fields:**

| Property | Dict-of-fields | Option Ω (flat + views) |
|---|---|---|
| Arithmetic per op | n_faces × numpy_call + Python dispatch | 1 numpy call |
| Memory | n_faces separate allocations | 1 contiguous allocation |
| Cache locality | Poor (scattered) | Excellent (contiguous) |
| MoC scaling (thousands of "faces") | Catastrophic Python overhead | Same O(1) numpy call |
| Per-face read/write | dict lookup + ndarray | dict lookup + slice view (no copy) |
| Cross-method uniformity | Bespoke per mesh | One `FaceLayout` Protocol |

**`DirectSumSpace` is deferred to Phase 3.6.** The Option Ω layout above is functionally equivalent to a DirectSumSpace's flat representation, but `BoundaryFlux` does not inherit from `Field` in this phase (it's a structured BUNDLE, not a single Field on a single space). When Phase 3.6 lands DirectSumSpace (for kinetics flux⊕precursors), `BoundaryFlux` may be refactored to `Field(values=_flat_buffer, space=DirectSumSpace(layout=...))` — same storage, more uniform algebra. That refactor is non-breaking; Option Ω here is forward-compatible.

### 3.5 The boundary auto-allocation question (revisited)

Per the §3.5 BoundaryFlux question from the previous exploration: the L2 `AngularFlux` needs to know how to construct a zero `BoundaryFlux` matching its phase-space mesh.

**Resolution under Depth B:** the SN adapter (`orpheus/sn/angular_flux_b1pp_adapter.py` or similar) provides:

```python
# orpheus/sn/angular_flux_b1pp_adapter.py
def zeros_for_sn_mesh(mesh: SNMesh, history_depth: int = 2) -> AngularFlux:
    """L3 factory: construct a zero AngularFlux + zero BoundaryFlux for the SN mesh.

    The L2 AngularFlux requires `boundary: BoundaryFlux` (no default — illegal
    states unrepresentable per coding-elegance Pattern 4). The factory is the
    SN-specific constructor.
    """
    boundary = _sn_zero_boundary_flux(mesh)  # geometry-conditional logic
    return AngularFlux(
        values=np.zeros((mesh.quad.N, mesh.ng, mesh.nx, mesh.ny)),
        space=_sn_angular_flux_space(mesh),
        mesh=mesh,
        boundary=boundary,
        history_depth=history_depth,
    )
```

The L2 base AngularFlux remains "every instance is fully-constructed; boundary is non-None by construction." The geometry-conditional `_sn_zero_boundary_flux` lives in L3.

Any future MoC analog would provide `zeros_for_moc_mesh(mesh: MoCMesh) -> AngularFlux` likewise.

### 3.6 `Source` siblings

`IsotropicSource` and `PerOrdinateSource` are physically distinct from fluxes (cross-type addition between flux and source is forbidden by the class-identity rule in `Field._check_partner`). They inherit from `Field` like every other typed field.

```python
# orpheus/transport/sources/isotropic_source.py
@dataclass(frozen=True, eq=False)
class IsotropicSource(Field):
    """Spatial source, isotropic in angle. Inherits Field algebra.

    space.units = "1/(cm³·s·sr·eV)" (rate density — distinct from ScalarFlux's flux density)
    """
    mesh: SNMesh
```

The class hierarchy mirrors the physical concept hierarchy:
- `ScalarFlux` vs `IsotropicSource`: same shape `(ng, nx, ny)`, DIFFERENT physical kind, DIFFERENT units. Cross-class addition raises by class identity.
- `AngularFlux` vs `PerOrdinateSource`: same shape `(N, ng, nx, ny)`, different physical kind, different units. Cross-class raises.

**The grand-report (§5.5) does NOT have a separate `Source` ABC** — sources are Field subclasses like everything else. The distinction is encoded by class identity + the `units` attribute on each space. Sufficient for Depth B.

**Cross-class same-units operations require named composition.** For Issue #201's iteration algebra (`AngularResidual - AngularSource = iteration_residual`), the result is a NAMED quantity built via an explicit factory, not via `-`:

```python
# Forbidden — same units but different physical concepts:
iteration_residual = angular_residual - angular_source   # raises TypeError

# Required — explicit named composition:
iteration_residual = IterationResidual.from_balance(lhs=angular_residual, rhs=angular_source)
```

This is the "complex frequency" analog: `iω - γ = s` is well-defined only as the NAMED complex-frequency construction; naive `ω + γ` is forbidden by convention.

Issue #205's full (storage × role) matrix follows the same pattern: each cell is a distinct class; cross-class operations require named compositions. Depth B ships the foundation; #205 fills in the matrix incrementally.

### 3.7 Singledispatch policy — KEEP the legitimate uses; RETIRE only ndarray handlers

`@singledispatchmethod` on operator `.apply` methods (and other operators that vary by input physics) is the canonical ORPHEUS pattern for "same operator concept, different math per storage/role type." It IS the grand-report's polymorphic operator language in concrete Python form. **Depth B does NOT retire this pattern wholesale.**

The existing dispatch tables in `orpheus/sn/scattering.py:817` (`ScatteringOperator.apply`) and `orpheus/sn/fission.py:154` (`FissionOperator.apply`) each have THREE handlers, two of which are legitimate and one is a migration crutch:

| Handler | Behavior | Verdict |
|---|---|---|
| `apply(AngularFlux)` | P_ℓ Galerkin reconstruction in per-ordinate magnitude (1/W projection at producer boundary). Tailored for SN sweep / GMRES / SI consumers. | **KEEP** — genuine structural dispatch. The math differs because the input physics differs. |
| `apply(ScalarFlux)` | P_0 + (n,2n) only, in iso scalar magnitude. No P_ℓ (scalar has no angular info); no 1/W (scalar consumers don't project). For diffusion / CP / kinetics-outer. | **KEEP** — same operator concept, different math because storage type carries less information. |
| `apply(np.ndarray)` | "Preserved for FD-matvec / probe-tests that bypass the type layer." | **RETIRE** — migration crutch for test convenience and legacy callers. |

**The retirement path for ndarray handlers** (executed in step D-I):

```python
# tests/sn/test_scattering_operator.py — before:
result = scattering.apply(np.zeros((N, ng, nx, ny)))

# tests/sn/test_scattering_operator.py — after:
psi = AngularFlux.from_ndarray(np.zeros((N, ng, nx, ny)), mesh)
result = scattering.apply(psi)  # routes through the typed AngularFlux handler
```

`TypedField.from_ndarray(arr, mesh)` is a constructor (cheap; just validates shape and wraps with a space). It's not a dispatch handler. Tests get ergonomic raw-ndarray construction; production routes through the typed dispatch.

**The policy:**

| Use case | Use singledispatch? |
|---|---|
| Operator `apply` where math varies by **storage type** (Angular vs Scalar vs Track vs Region) | **YES** — keep. Storage-axis dispatch is the canonical use. |
| Operator `apply` where math varies by **role type** (Flux vs Source vs Residual) | **YES** — keep. Will be added as #205 lands the full matrix. |
| Operator `apply` accepting `ndarray | TypedField` to ease testing | **NO** — retire. Replace with `TypedField.from_ndarray(...)` factory + typed dispatch. |
| Field arithmetic dunders (`__add__`, etc.) | **NO** — Python dunders already dispatch on self; class identity is the rule. |
| Cross-cutting concerns (serialization, visualization) | **YES** — when they land, but not in Depth B. |
| `Representation` abstraction (dense / sparse / sweep / TT) | **YES** — when grand-report §32.8 lands in a future phase. |

**Cross-references:** `coding-elegance` Pattern 1 (read-as-the-math via dunder); Pattern 5 (build the right primitive). The singledispatch tables ARE the right primitives for storage/role polymorphism.

---

## 4. The L1 Space algebra extensions

Per grand-report §6.1 and §32.4, `FunctionSpace` should support:

```python
S.dual()           # dual space
S * T              # tensor product
S + T              # direct sum
```

Land these alongside the Field migration. Concretely:

```python
# orpheus/numerics/space.py extensions

@dataclass(frozen=True)
class DualSpace(FunctionSpace):
    """The dual of a function space — same shape, transposed metric.

    For an L² Riesz-isometric space (the common case with diagonal weights),
    DualSpace(V) is structurally identical to V with the metric inverted.
    """
    primal: FunctionSpace

@dataclass(frozen=True)
class TensorProductSpace(FunctionSpace):
    """V = X ⊗ Ω ⊗ G — tensor product of independent axes.

    Per grand-report §5.3 and §32.4. The shape is the concatenation of
    component shapes; the metric is the outer product of component metrics.
    """
    components: tuple[FunctionSpace, ...]

@dataclass(frozen=True)
class DirectSumSpace(FunctionSpace):
    """V_total = V_1 ⊕ V_2 — DEFERRED to phase 3.6 (kinetics) per
    "unify after two instances" rule. Stub method here for forward
    compatibility; do not ship the impl in this phase."""
    pass  # stub
```

The dunders:

```python
class FunctionSpace:
    def __mul__(self, other: "FunctionSpace") -> "TensorProductSpace":
        return TensorProductSpace((self, other))

    def dual(self) -> "DualSpace":
        return DualSpace(...)

    # __add__ for DirectSumSpace deferred
```

**Decision point**: Should `__add__` be reserved for DirectSumSpace, or used for the algebraic-sum-of-fields-on-the-same-space (i.e., `V + V = V`, idempotent)? Per grand-report §6.1, `__add__` is direct sum. Avoid the ambiguity by NOT shipping `FunctionSpace.__add__` until Phase 3.6 lands DirectSumSpace.

---

## 5. The operator-side completion

The `Field`-on-`FunctionSpace` discipline only pays off if operators consume and return typed Fields, not bare ndarrays. Per §3.7's policy, the singledispatch tables themselves are KEPT (they're the canonical pattern for storage/role polymorphism); only the **ndarray handlers** within them are retired.

The audit identified the targets:

1. **`MomentProjection.apply(psi: np.ndarray) -> np.ndarray`** — its `codomain` is typed (`SphericalHarmonicSpace`) but the call signature is untyped. After Depth B, the apply path becomes:
   ```python
   class MomentProjection:
       @singledispatchmethod
       def apply(self, psi):
           raise TypeError(...)

       @apply.register
       def _(self, psi: AngularFlux) -> HarmonicMomentField:
           # The single typed path. Codomain check via psi.space.
           ...
   ```
   The bare-ndarray handler is RETIRED; tests construct AngularFlux via `AngularFlux.from_ndarray(arr, mesh)`.

2. **`LegendreMomentScattering.apply(psi: ndarray | HarmonicMomentField) -> ndarray | HarmonicMomentField`** — the dual-mode union signature is the leak. Refactor to `@singledispatchmethod` with one handler (`apply(HarmonicMomentField)`); retire the ndarray branch.

3. **`ScatteringOperator.apply` / `FissionOperator.apply` `@singledispatchmethod`** — KEEP the dispatch table. KEEP the `apply(AngularFlux)` and `apply(ScalarFlux)` handlers (genuine storage-axis polymorphism per §3.7). RETIRE the `apply(np.ndarray)` handler.

4. **`SNStreamingOperator.apply(psi: np.ndarray) -> np.ndarray`** (`orpheus/sn/operator.py:1506`) — packed-vector path for scipy.gmres. KEEP as the L4 NUMERICAL PRIMITIVE for the GMRES adapter (it's not really an "operator apply"; it's a packed-vector matvec for the iterative solver). The L2 typed entry is `StreamingOperator.apply(psi: AngularFlux)` at line 1801. Document the split as "L4 packed-vector primitive vs L2 typed entry."

These changes are NOT bit-identical at the API boundary (signatures change, ndarray entry points retire), but ARE bit-identical at the numerical level (same algorithm, same ndarray-level reductions inside each handler). The verification: every L1 MMS gate + every regression snapshot stays GREEN.

**Layer 2 dimensional check at solver construction.** When `SourceIteration` / `KrylovAcceleration` are constructed with `(L, S, F, q_ext)`, the `__init__` performs ONE dimensional-compatibility check:

```python
class SourceIteration:
    def __init__(self, L, S, F, q_ext):
        # Layer 2 (always runs, both debug and -O): O(1) per solver build
        for op in (L, S, F):
            if op.codomain.units.dimensionality != q_ext.space.units.dimensionality:
                raise UnitMismatchError(...)
        # ... rest of init
```

This is the **production safety net** for dimensional consistency — Layer 3 (`__debug__` asserts) is stripped in `-O`, so solver-construction-time validation guarantees the wiring is dimensionally sound regardless of optimization mode.

---

## 6. Sequencing

Recommended sequence — each step is a single commit, lands with verification gates green at every commit. Order chosen for monotonic dependency: L1 first, then field-by-field migration in size order (smallest first), then operator completion.

### Step D-A — L1 Field ABC (foundational)
- Create `orpheus/numerics/field.py` with the `Field` ABC per §3.2.
- Add unit tests in `tests/numerics/test_field.py` for the ABC algebra: same-space `+/-/*/neg`, cross-type rejection, cross-space rejection, inner_product, copy, linf, l2.
- Add `Field` to `orpheus/numerics/__init__.py` exports.
- **Verification**: `tests/numerics/test_field.py` PASS. P3.1 linter stays green. No production code changes (the ABC has no consumers yet).

### Step D-B — L1 Space algebra (optional in this phase)
- Add `FunctionSpace.__mul__` → `TensorProductSpace`.
- Add `FunctionSpace.dual()` → `DualSpace`.
- DEFER `__add__` (DirectSumSpace) per "unify after two instances" — it lands with Phase 3.6 if a second use case appears.
- **Verification**: `tests/numerics/test_space.py` extends with tensor-product + dual tests. PASS.
- **OPTIONAL**: this step can be skipped if §6's Field migration doesn't need TensorProductSpace yet (it doesn't — every field has a flat shape).

### Step D-C — Move `trace_space.py` into `numerics/spaces/`
- Per P3.2 plan (already scheduled): move `orpheus/numerics/trace_space.py` → `orpheus/numerics/spaces/trace_space.py`.
- Update imports; ship a one-cycle deprecation shim at `orpheus/numerics/trace_space.py`.
- **Verification**: every consumer still imports correctly. P3.1 linter stays green.

### Step D-D — Migrate `ScalarFlux` to L2 (simplest case)
- Create `orpheus/transport/__init__.py`, `orpheus/transport/fields/__init__.py`.
- Move `orpheus/sn/scalar_flux.py` → `orpheus/transport/fields/scalar_flux.py`.
- Refactor `ScalarFlux` to inherit from `Field`: add `space: FunctionSpace`, drop hand-coded dunders (inherit from `Field`), add `mesh: SNMesh` as additive field.
- Update `__post_init__` to consume `Field.__post_init__` via `super()`.
- Ship one-cycle re-export shim at `orpheus/sn/scalar_flux.py`.
- Migrate `tests/sn/test_scalar_flux.py` → `tests/transport/fields/test_scalar_flux.py` (or keep path, update import).
- **Verification**: every test that used `ScalarFlux` stays green. STRICT bit-identical (algorithm unchanged).

### Step D-E — Migrate `HarmonicMomentField` to L2 (cleanest gap)
- Move `orpheus/sn/harmonic_moment_field.py` → `orpheus/transport/fields/harmonic_moment_field.py`.
- Refactor to inherit from `Field`. Its `space` is a product of `SphericalHarmonicSpace.from_L(L)` and a per-cell-per-group space. Use `TensorProductSpace` from D-B if landed; otherwise use a flat `FunctionSpace` with the combined shape.
- Wire the connection: `HarmonicMomentField.from_moment_projection(M, psi)` constructs the result with `space=M.codomain ⊗ cell_group_space`.
- Ship one-cycle shim.
- **Verification**: `tests/sn/test_harmonic_moment_field.py` stays green. NEW test pinning `field.space == M.codomain ⊗ ...`. P3.1 linter green.

### Step D-F — Migrate `IsotropicSource` and `PerOrdinateSource` to L2
- Create `orpheus/transport/sources/__init__.py`.
- Move both modules; refactor to inherit from `Field` (Option β per §3.6).
- Cross-type `IsotropicSource + something` raises (already enforced by `Field._check_partner`).
- Ship shims; update tests.
- **Verification**: every source test stays green. STRICT bit-identical.

### Step D-G — Migrate `BoundaryFlux` to L2 (BundleOfEdges path per §3.4)
- Decide and document: **BundleOfEdges** (this phase) vs DirectSumSpace (Phase 3.6).
- Create `orpheus/transport/fields/boundary_edge.py` (a `BoundaryEdge(Field)` per face).
- Create `orpheus/transport/fields/boundary_flux.py` — a NON-Field bundle of edges.
- Refactor existing `BoundaryFlux` usage: every `bf.xmin_face` / `bf.xmax_face` access becomes `bf.edges["xmin"].values` / `bf.edges["xmax"].values`. This is a significant API surface change.
- Ship the SN-specific `zeros_for_sn_mesh` factory at `orpheus/sn/boundary_flux_zeros_factory.py` (or as a classmethod on `BoundaryFlux` that takes `mesh` — see decision point).
- Migrate consumers (~15 files per audit memo §3 / §4).
- **Verification**: every SN test stays green. STRICT bit-identical. The 10 pre-existing DD-regression failures stay AT the same failure set.

### Step D-H — Migrate `AngularFlux` to L2 (the most complex case)
- Move `orpheus/sn/angular_flux.py` (PARTIALLY) to `orpheus/transport/fields/angular_flux.py`.
- The L2 AngularFlux:
  - Inherits `Field`.
  - Carries `mesh: SNMesh`, `boundary: BoundaryFlux`, `history_depth`, `_history`.
  - REQUIRES `boundary` (no default — per §3.5).
  - Does NOT define `from_flat_with_traces` / `to_flat_with_traces`.
- The L3 SN adapter at `orpheus/sn/angular_flux_b1pp_adapter.py`:
  - Imports the L2 `AngularFlux`.
  - Installs `from_flat_with_traces` / `to_flat_with_traces` as class methods via `AngularFlux.from_flat_with_traces = staticmethod(...)` or similar at module load.
  - Exports `zeros_for_sn_mesh(mesh) -> AngularFlux`.
- The `same-class invariant` is preserved because `type(x)` is still the L2 `AngularFlux` — the SN adapter just installs methods on it.
- `orpheus/sn/__init__.py` re-exports the L3 factories.
- Migrate all consumers (~30 files).
- **Verification**: every AngularFlux test stays green. STRICT bit-identical. The `_is_ravellable` Protocol still detects via `hasattr(type(x), "from_flat_with_traces")` because the adapter installed the method.

### Step D-I — Wire operators to typed Fields (keep singledispatch, retire ndarray handlers)

Per §3.7 policy, the `@singledispatchmethod` dispatch tables on `ScatteringOperator.apply` and `FissionOperator.apply` STAY — they're the canonical storage-axis polymorphism pattern. Only the **ndarray handlers** retire.

- **`MomentProjection.apply`** — add typed handler `apply(AngularFlux) -> HarmonicMomentField`. Codomain check via `psi.space.shape[axes] == self.domain.shape`. Retire the bare-ndarray entry point.
- **`MomentProjection.apply_transpose(c: HarmonicMomentField) -> AngularFlux`** — typed analog.
- **`LegendreMomentScattering.apply`** — refactor the `ndarray | HarmonicMomentField` union into a `@singledispatchmethod` with ONE handler (`apply(HarmonicMomentField)`). Retire the ndarray branch.
- **`ScatteringOperator.apply`** — KEEP the dispatch table. KEEP `apply(AngularFlux)` and `apply(ScalarFlux)` handlers. **Retire only the `apply(np.ndarray)` handler.** Tests migrate to `AngularFlux.from_ndarray(...)` / `ScalarFlux.from_ndarray(...)` + typed dispatch.
- **`FissionOperator.apply`** — same pattern as ScatteringOperator. Keep dispatch table; retire ndarray handler.
- **`SNStreamingOperator.apply(psi: np.ndarray) -> np.ndarray`** (`orpheus/sn/operator.py:1506`) — KEEP as the L4 packed-vector primitive for scipy.gmres. Document as "the GMRES adapter primitive — not a Field operator." The typed `StreamingOperator.apply(psi: AngularFlux)` at line 1801 is the L2 user-facing entry.
- **Add `from_ndarray` factory methods** on every typed Field: `AngularFlux.from_ndarray(arr, mesh)`, `ScalarFlux.from_ndarray(arr, mesh)`, etc. Constructs a typed instance with the appropriate space. Replaces test-side bare-ndarray usage.

**Verification**: every L1 MMS gate stays green. The 10 DD-regression failures stay AT the same failure set. Type-checker (mypy/pyright if configured) reports zero new errors on the typed call sites. The `@singledispatchmethod` tables still have 2 handlers each (AngularFlux + ScalarFlux), down from 3 (was: + ndarray).

**Test-architect dispatch BEFORE D-I implementation** — per the proactive-trigger table in `subagent-handoff-protocol`. The MomentProjection rewire is the most complex operator change (its body works on bare ndarrays today; the typed wrapper needs structural verification against an analytical reference like k_∞ homogeneous reflective).

### Step D-J — Retire dead factories
- Delete `orpheus/numerics/space.py:angular_flux_space()`, `scalar_flux_space()`, `boundary_trace_space()` — zero callers per audit memo §8.
- Update `space.py` docstring to point to the new per-type factories (`AngularFlux.from_mesh`, etc.).
- **Verification**: clean delete; no test failures.

### Step D-K — Retire shims
- One merge cycle after each D-D / D-E / D-F / D-G / D-H lands, retire the `orpheus/sn/...` shims. Per `feedback_aggressive_retirement`.

---

## 7. Verification strategy

### 7.1 New tests added by Depth B

- `tests/numerics/test_field.py` — Field ABC algebra (D-A).
- `tests/numerics/test_space_algebra.py` — Space `*` / `.dual()` (D-B, if landed).
- `tests/transport/fields/test_scalar_flux.py` — migrated test (D-D).
- `tests/transport/fields/test_harmonic_moment_field.py` — migrated + new space-link test (D-E).
- `tests/transport/sources/test_isotropic_source.py`, `test_per_ordinate_source.py` — migrated (D-F).
- `tests/transport/fields/test_boundary_edge.py`, `test_boundary_flux.py` — new + migrated (D-G).
- `tests/transport/fields/test_angular_flux.py` — migrated, exercising L2 algebra (D-H).
- `tests/sn/test_angular_flux_b1pp_adapter.py` — pins the SN-adapter's `from_flat_with_traces` / `to_flat_with_traces` (D-H).

### 7.2 Tests that MUST stay green at every commit

- The P3.1 import-linter (`tests/test_layer_imports.py`) — every step must respect the layer contract.
- Every L1 MMS gate (slab/sphere/cylinder, ≥2G, heterogeneous): `tests/sn/test_mms_aniso.py`, `tests/sn/l1_analytical/...`.
- `tests/sn/test_typed_fields.py` — pins the algebra of every typed flux. STRICT bit-identical.
- `tests/numerics/test_iteration_angular_flux.py` — pins the `_is_ravellable` Protocol detection. STRICT bit-identical AT THE TYPE-IDENTITY LEVEL (the `type(x).from_flat_with_traces` lookup must keep returning the SN-installed method).
- The 10 pre-existing DD-regression failures stay AT the same failure set (verified by stash-and-rerun at every commit).

### 7.3 Bit-identity judgment per step

Per `vv-principles` §"Bit-identity vs principled-equivalence":

| Step | Bit-identity | Justification |
| --- | --- | --- |
| D-A | N/A | New code; no production change. |
| D-B | N/A | New code; no production change. |
| D-C | STRICT | Module move; same algorithm. |
| D-D | STRICT | ScalarFlux algebra unchanged. |
| D-E | STRICT | HarmonicMomentField algebra unchanged. |
| D-F | STRICT | Sources algebra unchanged. |
| D-G | STRICT | BoundaryFlux access surface changes; the BUFFERS are unchanged. |
| D-H | STRICT | AngularFlux algebra unchanged; class identity preserved. |
| D-I | PRINCIPLED | Operator signatures change from `ndarray` to `Field`. Underlying algorithm unchanged. |
| D-J | N/A | Dead-code retirement. |

The PRINCIPLED step (D-I) needs the three criteria:
1. **Principled**: every intermediate is a named typed Field. SATISFIED by design.
2. **Structurally-independent reference**: the L1 MMS gates remain the truth-of-record. SATISFIED.
3. **Drift dimensionally explainable**: zero drift expected (algorithm unchanged); ULP-level drift would be the only possible signal. If snapshots break by more than `iter_count × ULP`, the rewire is wrong.

### 7.4 Same-class invariant test

Add `tests/numerics/test_ravellable_protocol_same_class.py`:

```python
def test_ravellable_class_identity_preserved():
    """Pin: AngularFlux returned by from_flat_with_traces is the SAME class
    as the input template. Breaks _is_ravellable Protocol detection if
    violated."""
    # Construct an L2 AngularFlux via the SN adapter
    psi = zeros_for_sn_mesh(some_sn_mesh)
    assert type(psi).__name__ == "AngularFlux"
    flat = psi.to_flat_with_traces()
    reconstructed = type(psi).from_flat_with_traces(flat, some_sn_mesh)
    assert type(reconstructed) is type(psi)
```

This pins the contract that any future L2/L3 split (e.g., for MoC) MUST preserve.

### 7.5 Three-layer dimensional-enforcement story

The plan's units machinery enforces dimensional consistency at three layers, each with a distinct cost/coverage trade-off:

| Layer | When it runs | Cost | What it catches |
|---|---|---|---|
| **L1 — class identity in `Field._check_partner`** | Every dunder operation, BOTH `python -O -m pytest` and `pytest` | O(1) type comparison | `AngularFlux + AngularSource` (cross-class), `ScalarFlux + IsotropicSource` (cross-class). Foundational — never stripped. |
| **L2 — dimensional check at solver/operator construction** | Once per `SourceIteration.__init__` / `KrylovAcceleration.__init__` / `OperatorProduct.__init__`, BOTH modes | O(1) `pint.Unit.dimensionality` comparison per operator (~microseconds total per build) | Solver wired with dimensionally-inconsistent operators (e.g., `L.codomain.units` doesn't match `q_ext.space.units`). Production safety net. |
| **L3 — defensive `assert` in dunders + composition sites** | Every dunder, ONLY debug-mode `pytest` (stripped in `-O`) | Zero in production; small in development | Class/units misdesign (a class declared with wrong units); manual-construction bugs where someone hand-builds a Field with mismatched space.units. Defense in depth. |

**CI matrix** (per `[[default-test-mode-is-optimize]]` — `-O` is canonical):

```yaml
test-optimized (PRIMARY — runs the full suite):
  run: python -O -m pytest tests/         # Layer 3 stripped — production-like
test-debug (NARROW — exercises Layer 3 catches only):
  run: pytest tests/ -m "layer3_dynamic"  # Layer 3 active — pins debug-only behaviour
```

Most tests are mode-agnostic and pass under both invocations. Tests that exercise Layer 3 dimensional-check-on-the-fly are marked with `@pytest.mark.skipif(not __debug__, …)`; their companion strip-verification tests use `@pytest.mark.skipif(__debug__, …)`. The pair toggles which one runs based on mode — no test runs twice.

**Layer 3 implementation** uses `assert` statements (idiomatic, bytecode-stripped in `-O`) rather than `if __debug__:` blocks (equally stripped but more verbose). Both produce the same bytecode in `-O` mode. pytest's own assertion rewriting still works under `-O` for `assert` inside `test_*.py` files (rewritten at AST level); only production-code asserts under `orpheus/` get stripped — exactly the asymmetry that makes Layer 3 work.

**`pint` import policy.** `pint >= 0.20` is a normal production dependency (added to `pyproject.toml` in step D-A). `import pint` at module load — ~100ms once at startup, irrelevant on simulation timescales. NO lazy import; lazy import was considered and rejected because it would force `units: str` (losing typed-ness) and add string-parsing overhead at every dimensional check.

---

## 8. Risk register

| Rank | Risk | Mitigation |
| --- | --- | --- |
| **1** | `dataclass` field-order constraints in `Field` ABC + subclasses (mandatory fields after defaulted ones). | Use `kw_only=True` on the dataclass decorator (Python 3.10+; project is on 3.14). |
| **2** | Same-class invariant broken by L2/L3 module-level injection getting tangled. | Strict test (`test_ravellable_protocol_same_class.py`) per §7.4. Plus: prefer `class AngularFlux` defined at L2 with the L3 adapter installing methods, NOT a subclass. |
| **3** | `BoundaryFlux` BundleOfEdges API change breaks too many consumer sites. | Audit-driven migration; the audit memo identified ~15 files. One-cycle deprecation shim on `bf.xmin_face` accessors. Migrate consumers in lockstep. |
| **4** | Step D-I retires `apply(np.ndarray)` handlers but tests still pass bare ndarrays. | Add `TypedField.from_ndarray(arr, mesh)` factory in step D-A; migrate test call sites in lockstep with each typed field's D-D/D-E/D-F/D-G/D-H migration. The singledispatch table itself STAYS (per §3.7); only the ndarray entry leaves. |
| **4b** | `pint` dependency adds installation surface. | Single-package pure-Python dep, no C extensions, MIT-licensed, well-maintained. Pin `pint >= 0.20` for stable API. |
| **5** | DD-regression suite already has 10 pre-existing failures. New Depth B failures would be hidden in the noise. | Capture the exact pre-existing failure set as `tests/sn/regression/known_failures.txt` BEFORE step D-A. Every step asserts the failure set is unchanged. |
| **6** | `MomentProjection.apply` rewire (D-I) is the most complex operator change — its codomain is typed but its body works on bare ndarrays. | Test-architect dispatch BEFORE D-I implementation; pin a structurally-independent reference (homogeneous reflective k_∞ analytic limit) against the new wiring. |
| **7** | Phase 3.4 (Problem/Solver split) depends on the Field/Space algebra being in place. Delaying Depth B delays 3.4. | Ship Depth B first; 3.4 plan can be drafted in parallel but implementation follows. |

---

## 9. What this plan DEFERS

Out-of-scope for Depth B; do NOT chase them in this phase:

1. **`DirectSumSpace`** — defer to Phase 3.6 (kinetics) when there's a second use case (flux + precursors).
2. **`MethodSpace` hierarchy** (grand-report §33) — defer to a future phase that introduces `DiscreteOrdinatesPhaseSpace` as the L2 method-space carrier. The `mesh: SNMesh` field on every L2 type is the stand-in.
3. **`Representation` abstraction** (grand-report §32.8 — dense / sparse / matrix-free / sweep / tensor-train) — defer to Phase 3.5+ when there's a second representation.
4. **Boundary registry / `BoundaryRealizer` self-registration** (grand-report §16A.11) — already partially present at `orpheus/geometry/boundary.py`. Out-of-scope refinement.
5. **Symmetry sectors / orbifold spaces** (grand-report §16B, §33.11) — entire chapter; defer.
6. **MoC `CharacteristicTrackSpace`, sheaf machinery** (grand-report §11A, §33.3) — MoC is its own future phase.

These deferrals are deliberate. The audit memo's "unify after two instances" rule applies.

---

## 10. Plan checkpoint — to be filled at approval time

- [ ] User approval of architecture (Field ABC + L2 fields + L1 space dunders + boundary BundleOfEdges path).
- [ ] User approval of sequencing (D-A through D-K).
- [ ] User approval of deferrals (§9).
- [ ] User confirmation: ship as ONE commit per step (D-A, D-B, ..., D-K), not bundled.
- [ ] User confirmation: dispatch test-architect BEFORE D-I (operator rewire) per the proactive-trigger table.
- [ ] User confirmation: dispatch test-architect BEFORE D-H (AngularFlux migration) for the same-class invariant verification spec.

Once these checkpoints are confirmed, this plan is the source of truth for the fresh-context Depth B implementation session.

---

## 11. Exit Route — where this detour leads when complete

**Depth B is a detour.** It re-scopes step P3.3 of the parent plan into a foundational unification (Field-on-FunctionSpace + units). When all 11 steps (D-A through D-K) complete, execution **returns to the parent plan**:

- **Parent plan:** `.claude/plans/moment_space_and_layering_plan.md`
- **Hand-off point:** Parent plan step **P3.4 — Problem / Solver split (greenfield per CC.7)**

### 11.1 Hand-off contract

When Depth B is COMPLETE, the following invariants hold and are pickup conditions for P3.4:

1. **`Field` ABC ships at L1** (`orpheus/numerics/field.py`) with the algebra described in §3.2.
2. **`FunctionSpace.units: pint.Unit | None`** exists and is consumed by `__eq__` (dimensionality comparison) and `Field._check_partner` (Layer 3 assert).
3. **`orpheus/transport/` exists** as the L2 package with `fields/` (AngularFlux, ScalarFlux, HarmonicMomentField, BoundaryFlux, BoundaryFaceFlux) and `sources/` (IsotropicSource, PerOrdinateSource) subpackages.
4. **The L3 SN adapter installs `from_flat_with_traces` / `to_flat_with_traces`** on the L2 `AngularFlux` class at module-import time. The same-class invariant test (`tests/numerics/test_ravellable_protocol_same_class.py`) is GREEN.
5. **`@singledispatchmethod` dispatch tables stay** on `ScatteringOperator.apply`, `FissionOperator.apply` with handlers for `AngularFlux` and `ScalarFlux` (and future `AngularSource`/`AngularResidual` when #201 lands). The `apply(np.ndarray)` handlers are RETIRED. Each typed field has a `.from_ndarray(arr, mesh)` factory.
6. **`SourceIteration.__init__` / `KrylovAcceleration.__init__` do Layer 2 dimensional checks** at construction. Production runs (`python -O`) get this safety net.
7. **All four completed Phase 3 steps** (P3.1 import-linter, P3.5 range→codomain, P3.0 layering docs, P3.2 spherical_harmonics shim retirement) remain GREEN. The 10 pre-existing DD-regression failures stay at the same failure set.
8. **`pint >= 0.20`** is in `pyproject.toml`. CI matrix runs both `pytest` and `python -O -m pytest`.
9. **GitHub Issues #201 and #205 remain OPEN** — they are downstream consumers of Depth B, not closed by it.

### 11.2 What P3.4 picks up

P3.4 (parent plan §P3.4, around line 867 of `moment_space_and_layering_plan.md`) introduces:
- `Problem` ABC at L2 (`orpheus/transport/problems/`) — declarative descriptions: `CriticalityProblem`, `FixedSourceProblem`, `AlphaEigenproblem`, `InitialValueProblem`
- `Solver` ABC at L1 (`orpheus/numerics/solvers/`) — iterative algorithms: `PowerIteration` (formerly `KEigenvalue`), `Arnoldi`, `SourceIteration`, time-steppers
- `homogeneous/solver.py:26`'s `power_iteration` legacy retires (CC.8 close-out)

P3.4 directly consumes Depth B's typed Fields:
- `CriticalityProblem` carries `loss: LinearOperator` (codomain = AngularResidual space) and `fission: LinearOperator` (codomain = AngularSource space)
- `PowerIteration.solve(problem) -> AngularFlux` — the typed flux is the result
- `FixedSourceProblem` carries `(L, S, F, q_ext: AngularSource)` — every component is dimensionally explicit

The Layer 2 construction-time dimensional check (Depth B §5) is the FIRST check `PowerIteration` does on its problem — verifying the operator algebra makes dimensional sense before any iteration runs.

### 11.3 Sequence after P3.4

The parent plan's REVISED Phase 3 sequencing (per its §"Sequencing within Phase 3 (REVISED post-Phase-1 per QA learnings)"):

```text
P3.1 ✓  → P3.5 ✓  → P3.0 ✓  → P3.2 ✓  → [Depth B = P3.3] → P3.4 → P3.6
                                          ↑                  ↑      ↑
                                          THIS PLAN          NEXT   LAST
```

**P3.6 (kinetics restructure)** is the LAST step. It dissolves `kinetics/` into `transport/problems/initial_value.py` + `numerics/solvers/time_stepping.py`. It also lands `DirectSumSpace` (the deferred grand-report §5.3 primitive) when it needs `flux ⊕ precursors`. P3.6 closes Phase 3.

### 11.4 Branch / merge after Phase 3

The parent plan (`moment_space_and_layering_plan.md`) lives on branch `refactor/moment-space-and-layering`. When all of Phase 3 is complete (P3.6 LANDED), the branch fast-forwards into `refactor/sn-operator-algebra` (the parent feature branch). Per the user's strategic decision (recorded in conversation 2026-05-26): the moment-space plan is treated as a single coherent unit; merge happens at plan-completion, not mid-plan.

After the moment-space-and-layering branch merges into `refactor/sn-operator-algebra`:
- The remaining R-1 work (parallel session on `refactor/sn-operator-algebra`) continues with the typed Field/Space spine in place
- Issues #201 and #205 become implementable on top of the Depth-B foundation
- The grand-report's deferred chapters (Representation §32.8, MoC sheaves §11A, symmetry sectors §16B, MethodSpace hierarchy §33) become tractable

### 11.5 Pickup signal for P3.4

A fresh session resuming after Depth B should:
1. Read `git log --oneline -20` and confirm the Depth-B commits (D-A through D-K) are present.
2. Verify all 9 invariants in §11.1 by running the relevant test commands.
3. Open `.claude/plans/moment_space_and_layering_plan.md` at §P3.4 (~line 867).
4. Open `.claude/agent-memory/test-architect/phase3_verification_plan.md` §2.6 (P3.4 verification spec, written before Phase 3 started — still authoritative).
5. Dispatch test-architect for P3.4 BEFORE implementation begins (proactive trigger per `subagent-handoff-protocol`: Problem/Solver split is an operator-algebra carve crossing subsystem boundaries).
6. Begin P3.4 implementation.

---

## Pointers

- **Parent plan (Depth B's exit):** `.claude/plans/moment_space_and_layering_plan.md` — Phase 3 sequencing. Depth B replaces step P3.3; execution returns to step **P3.4** when Depth B completes (see §11 above).
- **Branch state:** `git log --oneline -10` on `refactor/moment-space-and-layering`.
- **Audit memo:** `.claude/agent-memory/explorer/function_space_typed_field_audit.md` (FunctionSpace + typed-field current state inventory — load-bearing for D-A through D-K).
- **Test-architect Phase 3 verification plan:** `.claude/agent-memory/test-architect/phase3_verification_plan.md` — note: §2.5 (P3.3) is SUPERSEDED by this Depth B plan; §2.6 (P3.4) remains authoritative for the post-exit work; §2.1–2.4 are accurate for the four completed Phase 3 steps.
- **Grand report:** `.claude/plans/neutron_transport_grand_report_v3.md`. Sections relevant to Depth B: §2 ontology, §5.3 Space hierarchy, §5.5 Field hierarchy, §5.7 Operator hierarchy, §6.1 Space dunders, §16A boundary primitives, §32.4-32.6 spaces/fields/operators specs.
- **Downstream consumers (issues OPEN throughout Depth B; closed in later phases):**
  - **GitHub Issue #201** — `SN: dimensional typing — split AngularFlux from AngularSource / AngularResidual`. The angular slice of #205. Depth B's units machinery makes #201 machine-checkable.
  - **GitHub Issue #205** — `Cross-method field architecture: storage × role typing (Angular/Scalar/Track/Region × Flux/Source/Residual)`. The full 12-cell matrix. Depth B is the foundation.
- **Phase 1 + Phase 3 prior commits:** `1ab6233..ea02ab5` (Phase 1 close-out) + `bc99ff6..09c4241` (Phase 3 mechanical steps: P3.1 import-linter, P3.5 range→codomain rename, P3.0 layering docs, P3.2 spherical_harmonics shim retirement).
- **Cardinal rules from memory:** `feedback_unify_after_two_instances`, `feedback_aggressive_retirement`, `feedback_no_method_implementer_for_surgical_carves`, `feedback_retirement_means_test_migration`, `feedback_architecture_forward_not_legacy_fit`.
