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

9. **Operators returning bare ndarray** (`LegendreMomentScattering.apply`, `MomentProjection.apply` — the LAST is the most embarrassing since its `codomain` is already typed) are the migration target for Depth-B's "operator-side completion" phase (Step D-I below). NOTE: `SNStreamingOperator` (the legacy `L + C` bundle from pre-Phase-G) is on the retirement queue and is OUT OF SCOPE for Depth B's typed-fields wire-up — the modern algebraic path is `(L + C).apply` where `L = StreamingOperator` and `C = CollisionOperator` compose into `InvertibleOperator(OperatorSum)` (see `orpheus/sn/operator.py:1801, 2170, 2412`). `StreamingOperator.apply` already accepts typed `AngularFlux` per the R-1 Step 3d signature.

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

Each subclass adds DOMAIN-SPECIFIC fields (e.g. `mesh`, `L`) on top of `Field`'s `(values, space)`. The dunders are inherited unchanged. Only the constructor + class-specific methods differ.

**Architecture decision (recorded 2026-05-27)**: `AngularFlux` is a **pure Field**, symmetric with `ScalarFlux`. It does NOT carry `boundary` or `_history`. The composite iteration state (psi + boundary + history) lives in a **separate container class** introduced at §3.8.

```python
# orpheus/transport/fields/angular_flux.py
@dataclass(frozen=True, eq=False, kw_only=True)
class AngularFlux(Field):
    """L2 typed flux: a pure Field on the (N, ng, nx, ny) phase space.

    Carries `values + space + mesh` and inherits Field algebra. No
    boundary, no history, no iteration metadata — symmetric with
    ScalarFlux. Composite iteration state lives in TransportState
    (see §3.8).
    """
    mesh: SNMesh  # added field — discretisation handle (future TransportMesh Protocol)

    # __post_init__ extends Field.__post_init__ via super() — adds mesh-vs-space
    # consistency check.

    # NO from_flat_with_traces / to_flat_with_traces — those live on the
    # TransportState container, where the trace-data target makes sense.
```

Subtle: `dataclass` inheritance with `@dataclass(frozen=True)` requires every parent field to have a default OR `kw_only=True`. We use `kw_only=True` (Python 3.10+; project is on 3.14). Pattern established by `IsotropicSource` / `PerOrdinateSource` in D-F.

### 3.4 `BoundaryFlux` — pure Field with flat-buffer storage (Option Ω)

**Architecture decision (recorded 2026-05-27)**: `BoundaryFlux` is a **pure Field** with a flat backing buffer + `FaceLayout` descriptor. The pre-D-G framing ("not a Field, structured bundle") is RETIRED — Field-on-flat-storage gives both single-space algebra (Field-inherited dunders) AND per-face slice-view access (via FaceLayout). The interior-cache conflation that the pre-D-G 2-D `xmin_xmax_buf` / `ymin_ymax_buf` buffers carried is split out into a sweep-private `SweepScratch` type — BoundaryFlux is IMMUTABLE.

`BoundaryFaceFlux` is the per-face slice-view type, also a Field (on a single face's trace space). Construction is via `BoundaryFlux.faces[face_name]` — returns a Field view, no copy.

**Naming convention:**
- `FaceFlux` (FUTURE, potential L2) — flux on ANY face (internal or boundary). Not in scope for this phase.
- `BoundaryFaceFlux` — flux on ONE boundary face's inflow trace space. This phase introduces it.
- `BoundaryFlux` — single Field over all boundary faces of a mesh (flat storage + FaceLayout).

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
@dataclass(frozen=True, eq=False, kw_only=True)
class BoundaryFlux(Field):
    """L2 typed boundary trace flux: a pure Field over all boundary faces.

    Storage is a flat backing buffer; per-face access via the FaceLayout
    descriptor returns slice views (no copies). Field algebra is
    inherited — `__add__`, `__sub__`, scalar `*`, `__neg__`, `__truediv__`
    all operate on the flat buffer in ONE numpy call.

    IMMUTABLE. The pre-D-G mutable write-through path (sweep persistent
    BC state) is replaced by sweep-side SweepScratch + functional
    reconstruction at iteration boundaries.
    """
    # values (flat backing) and space (FunctionSpace with shape=(layout.total_size,))
    # inherited from Field. Added fields:
    layout: FaceLayout                # mesh-provided descriptor
    mesh: SNMesh                      # discretisation handle

    def __post_init__(self) -> None:
        super().__post_init__()
        if self.values.shape != (self.layout.total_size,):
            raise ValueError(
                f"BoundaryFlux: flat buffer shape {self.values.shape} "
                f"!= ({self.layout.total_size},)"
            )

    def _check_partner(self, other: object) -> None:
        super()._check_partner(other)
        if self.mesh is not other.mesh:  # type: ignore[attr-defined]
            raise ValueError("BoundaryFlux mesh-binding mismatch")
        if self.layout is not other.layout:  # type: ignore[attr-defined]
            if self.layout != other.layout:
                raise ValueError("BoundaryFlux layout mismatch")

    # ── Per-face access (slice views, no copies) ────────────────────────
    @property
    def faces(self) -> Mapping[str, BoundaryFaceFlux]:
        return {
            name: BoundaryFaceFlux(
                values=slot.slice_view(self.values),
                space=self._space_for_face(name),
            )
            for name, slot in self.layout.faces.items()
        }

    @classmethod
    def zeros_for_sn_mesh(cls, mesh: SNMesh) -> "BoundaryFlux":
        """SN-geometry-conditional zero-construction. Lives here as a
        classmethod rather than at sweep call sites — the geometry
        branching is part of the field's construction discipline.
        """
        layout = mesh.boundary_face_layout
        space = FunctionSpace(name="sn_boundary_flat", shape=(layout.total_size,))
        return cls(
            values=np.zeros(layout.total_size),
            space=space,
            layout=layout,
            mesh=mesh,
        )
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

### 3.5 The boundary auto-allocation question — resolved by the container

The pre-D-G framing — "the L2 `AngularFlux` needs to know how to construct a zero `BoundaryFlux` matching its phase-space mesh" — is **dissolved** by the architecture decision recorded at §3.3 and §3.8: AngularFlux does NOT carry boundary, so it does not need to know how to construct one. The composite construction lives on the **container**:

```python
# orpheus/sn/transport_state_b1pp_adapter.py
@classmethod
def zeros_for_sn_mesh(cls, mesh: SNMesh, history_depth: int = 2) -> "TransportState":
    """L3 factory: construct a zero TransportState (psi + boundary + empty history)
    for the SN mesh. The composite is the iteration state; algebraically pure
    AngularFlux and BoundaryFlux components have their own zero factories
    (AngularFlux.from_mesh + BoundaryFlux.zeros_for_sn_mesh) and stand alone.
    """
    return cls(
        psi=AngularFlux(
            values=np.zeros((mesh.quad.N, mesh.ng, mesh.nx, mesh.ny)),
            space=_sn_angular_flux_space(mesh),
            mesh=mesh,
        ),
        boundary=BoundaryFlux.zeros_for_sn_mesh(mesh),
        _history=(),
        history_depth=history_depth,
    )
```

The L2 `AngularFlux` constructor is now trivially symmetric with `ScalarFlux` — just `(values, space, mesh)`. The geometry-conditional construction logic lives ONCE, on the container's L3 adapter. Any future MoC analog provides `TransportState.zeros_for_moc_mesh(mesh)` likewise — the container is method-agnostic but its zero factories are method-specific.

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

### 3.8 Container architecture — `TransportState` (the composite iteration state)

**Architecture decision recorded 2026-05-27 from the conversation.** The user's verbatim framing:

> "I actually think the AngularFlux class, having boundary condition inside and history, is not appropriate. It would seem to me that AngularFlux should not have BoundaryFlux and time history (just like ScalarFlux) and that we need a container class (maybe TimedFullField) that can have a field with boundary values and history. This container class could then contain AngularFlux + BoundaryFlux and history inside, and operates with a dunder."

The container separates THREE responsibilities that the pre-D-H `AngularFlux` conflated:

| Responsibility | Pre-D-H location | Post-D-H location |
|---|---|---|
| Interior phase-space flux values | `AngularFlux.values` | `AngularFlux.values` (unchanged) |
| Boundary trace values | `AngularFlux.boundary: BoundaryFlux` | `TransportState.boundary: BoundaryFlux` |
| Iteration history (Krylov / Anderson) | `AngularFlux._history: tuple` | `TransportState._history: tuple[TransportState, ...]` |

**Why this is principled** (the cardinal-rule 2 architectural argument):

1. **`ScalarFlux` / `AngularFlux` symmetry.** Pre-D-H, ScalarFlux is a clean Field (values+space+mesh) and AngularFlux is bloated with two extra responsibilities. Post-D-H they're both just Fields. The asymmetry that's been bugging the architecture through Depth B disappears.
2. **Composition over inheritance.** Per Issue #207's three-way split (data composition / trait composition / named composition), `TransportState` is a **composition constructor** — it builds a composite state from existing Fields. It's a sibling of `TensorProductSpace`'s role at L1 (which also constructs new types from existing ones).
3. **Future alignment with `DirectSumSpace`.** The container's natural backing in grand-report §5.3 is `Field(values=DirectSumStorage, space=DirectSumSpace(psi.space, boundary.space))`. `DirectSumSpace` is DEFERRED to Phase 3.6 per `feedback_unify_after_two_instances` (the second instance is kinetics' `flux ⊕ precursors`). D-H ships `TransportState` as a structured bundle with delegate-style dunders; the Phase 3.6 promotion is non-breaking.
4. **Dissolves the 2-D conflation.** Once boundary is its own Field at the container level, BoundaryFlux is just `(face_values, face_space, mesh)` — no interior cache, no buffer-mixing. The pre-D-G 2-D `xmin_xmax_buf` / `ymin_ymax_buf` buffers split naturally: pure boundary trace → `BoundaryFlux`; interior wavefront cache → sweep-private `SweepScratch`.

**Naming.** The user's working name is `TimedFullField`. Alternates considered:

| Candidate | Pros | Cons |
|---|---|---|
| `TimedFullField` | User's working name; emphasizes time-stepping potential | "Timed" implies temporal stepping but primary use today is Krylov/SI iteration state |
| `TransportState` | Concise; reads as the math (the state of a transport iteration) | Slightly generic |
| `IterationState` | Captures the actual use case (iteration metadata) | Couples concept to the solver mechanism, not the physics |
| `FullTransportField` | Descriptive | Wordy; "field" is misleading because the container is NOT a Field |

**Recommendation: `TransportState`** — the composite IS the state of the transport problem at a given iteration index. The name reads as the math at consumer sites: `state = source_iteration.step(state, source)`, `new_state = (L+C).solve(rhs_state)`. Final name to be decided pre-D-H.1 implementation.

**The class shape:**

```python
# orpheus/transport/state.py  (or orpheus/transport/transport_state.py)
@dataclass(frozen=True, kw_only=True)
class TransportState:
    """Composite iteration state: pure flux + pure boundary trace + history.

    NOT a Field. A structured bundle with delegate-style dunders that
    propagate algebra to the contained fields. The history shift register
    is iteration metadata — does NOT participate in algebra (algebra
    results carry empty history).

    Architectural sibling of TensorProductSpace at L1: a composition
    constructor that builds a composite type from existing typed Fields.
    Future Phase 3.6 may promote to Field-with-DirectSumSpace backing
    (non-breaking — public API stays).
    """
    psi: AngularFlux
    boundary: BoundaryFlux
    _history: tuple["TransportState", ...] = ()
    history_depth: int = 2

    def __post_init__(self) -> None:
        if self.psi.mesh is not self.boundary.mesh:
            raise ValueError(
                "TransportState requires psi and boundary on the same mesh"
            )

    # ── Algebra (delegates to contained fields; history dropped) ────────
    def _check_partner(self, other: "TransportState") -> None:
        if type(other) is not TransportState:
            raise TypeError(
                f"TransportState arithmetic requires TransportState partner; "
                f"got {type(other).__name__}"
            )
        if self.psi.mesh is not other.psi.mesh:
            raise ValueError("TransportState arithmetic across different meshes")

    def __add__(self, other: "TransportState") -> "TransportState":
        self._check_partner(other)
        return TransportState(
            psi=self.psi + other.psi,
            boundary=self.boundary + other.boundary,
            _history=(),  # algebra results have no iteration history
            history_depth=self.history_depth,
        )

    def __sub__(self, other: "TransportState") -> "TransportState":
        self._check_partner(other)
        return TransportState(
            psi=self.psi - other.psi,
            boundary=self.boundary - other.boundary,
            _history=(),
            history_depth=self.history_depth,
        )

    def __neg__(self) -> "TransportState":
        return TransportState(
            psi=-self.psi,
            boundary=-self.boundary,
            _history=(),
            history_depth=self.history_depth,
        )

    def __mul__(self, scalar: float) -> "TransportState":
        return TransportState(
            psi=self.psi * scalar,
            boundary=self.boundary * scalar,
            _history=(),
            history_depth=self.history_depth,
        )

    def __rmul__(self, scalar: float) -> "TransportState":
        return self.__mul__(scalar)

    def __truediv__(self, scalar: float) -> "TransportState":
        return self.__mul__(1.0 / float(scalar))

    # ── History shift-register (iteration metadata) ─────────────────────
    def advance(self, new_psi: AngularFlux, new_boundary: BoundaryFlux) -> "TransportState":
        """Push current (psi, boundary) into history; new (psi, boundary) become current.

        The shift-register pattern that pre-D-H AngularFlux exposed via
        `psi << new` — lifted to the container where the history lives.
        """
        current_snapshot = TransportState(
            psi=self.psi,
            boundary=self.boundary,
            _history=(),
            history_depth=self.history_depth,
        )
        new_history = (current_snapshot, *self._history)[: self.history_depth]
        return TransportState(
            psi=new_psi,
            boundary=new_boundary,
            _history=new_history,
            history_depth=self.history_depth,
        )

    def at_lag(self, lag: int) -> "TransportState":
        """Read the state at iteration lag (0 = current, 1 = previous, ...)."""
        if lag == 0:
            return self
        if lag - 1 >= len(self._history):
            raise IndexError(
                f"TransportState: lag {lag} out of history (depth={len(self._history)})"
            )
        return self._history[lag - 1]
```

**Krylov / GMRES adapter at L3.** The packed-vector adapter (`to_flat_with_traces` / `from_flat_with_traces`) now targets the CONTAINER, not the inner flux. The `_is_ravellable` Protocol at `iteration.py:163-176` detects `TransportState`. This is structurally cleaner: the flat packed vector is `concat([psi.values.ravel(), boundary.values])` — a direct sum representation, which is exactly what the future DirectSumSpace Field promotion will formalize.

```python
# orpheus/sn/transport_state_b1pp_adapter.py
def _install_ravel_unravel():
    """L3 import-time injection: install to_flat_with_traces / from_flat_with_traces
    on the L2 TransportState class. Preserves the same-class invariant
    (numerics/iteration.py:163-176): type(x) is always TransportState.
    """
    def to_flat_with_traces(self: TransportState) -> NDArray:
        return np.concatenate([self.psi.values.ravel(), self.boundary.values])

    def from_flat_with_traces(cls, flat: NDArray, template: TransportState) -> TransportState:
        n_psi = template.psi.values.size
        psi_values = flat[:n_psi].reshape(template.psi.values.shape)
        boundary_values = flat[n_psi:]
        return cls(
            psi=replace(template.psi, values=psi_values),
            boundary=replace(template.boundary, values=boundary_values),
            _history=(),
            history_depth=template.history_depth,
        )

    TransportState.to_flat_with_traces = to_flat_with_traces
    TransportState.from_flat_with_traces = classmethod(from_flat_with_traces)

_install_ravel_unravel()
```

**Sweep API impact.** Pre-D-H, the sweep consumes `psi: AngularFlux` (reading `psi.boundary` internally) and produces a new `AngularFlux` (writing `.boundary` write-through). Post-D-H, the sweep consumes `state: TransportState` and produces a new `TransportState`. Internally the sweep reads `state.psi.values` and `state.boundary.values`, produces fresh ndarray outputs, and constructs a new `TransportState` (immutable functional path). The pre-D-G mutable write-through is replaced by the sweep-side `SweepScratch` for any per-iteration cache needs.

**SourceIteration / InvertibleOperator.solve API impact.** Today (R-1 Step C): `(L+C).solve(rhs: AngularFlux) -> AngularFlux` reads `rhs.boundary`. Post-D-H: `(L+C).solve(rhs: TransportState) -> TransportState`. The rhs container carries both the inhomogeneous source flux AND the inflow boundary trace; the solve returns the converged container with updated psi + updated boundary trace + new history entry.

**Same-class invariant.** Today `_is_ravellable` checks `hasattr(type(x), "from_flat_with_traces")` on `AngularFlux`. Post-D-H it checks on `TransportState`. The mechanism is unchanged — only the target type changes. The protocol detection still works; the carved-off responsibility lives on the correct type.

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

4. **`StreamingOperator.apply` (modern, `orpheus/sn/operator.py:1801`)** already accepts typed `AngularFlux` per the R-1 Step 3d signature. The bare-ndarray branch inside it is the L4 packed-vector path for scipy.gmres — out of scope for Depth B (retires with the GMRES adapter rewrite, separately). Note: `SNStreamingOperator` (the legacy `L + C` bundle, `operator.py:1329`) is OUT OF SCOPE here — it's on the retirement queue from Phase G substep 3+4.c. The modern algebraic path is `(L + C).apply` where `L + C` composes into `InvertibleOperator` (`operator.py:2412`) carrying both `apply` (forward matvec) and `solve` (transport sweep). Depth B's typed-field wire-up targets `StreamingOperator` (the modern leaf), not `SNStreamingOperator` (the legacy bundle).

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

**Execution status (2026-05-27):**

| Step | Status | Commit | Notes |
|---|---|---|---|
| D-A | LANDED | `1e0bb98` | Field ABC + FunctionSpace.units + pint |
| D-C | LANDED | `469ea15` | trace_space moved to numerics/spaces/ |
| D-D | LANDED | `53bc986` | ScalarFlux migrated to L2 + inherits Field |
| D-B | LANDED | — | TensorProductSpace L1 primitive (load-bearing) |
| D-B+1 | LANDED | — | specular BC tensor-network — first production instance |
| D-E | LANDED | `6123422` | HarmonicMomentField migrated to L2; uses TensorProductSpace |
| D-F | LANDED | `efdd4fe` | IsotropicSource + PerOrdinateSource migrated to L2; cross-class dunder preserved (Issue #207 refinement) |
| D-G | NEXT | — | **REVISED 2026-05-27**: BoundaryFlux as pure Field (Option β stage 1); split SweepScratch |
| D-H.1 | PENDING | — | **NEW 2026-05-27**: TransportState container + consumer migration (Option β stage 2a) |
| D-H.2 | PENDING | — | **NEW 2026-05-27**: Retire legacy AngularFlux.boundary/_history (Option β stage 2b — HARD DEADLINE per user) |
| D-I through D-K | PENDING | — | Operator typing, dead-factory retirement, shim retirement |

**Sequencing rationale.** D-A had to land first (Field ABC). D-C (trace_space move) and D-D (ScalarFlux migration) followed because they only depended on D-A. D-B was the next foundational step — TensorProductSpace unlocked Wave T (which follows Depth B) and enabled D-E's HarmonicMomentField to use product-space typing cleanly. D-B+1 shipped the first tensor-network instance. D-E + D-F migrated the simpler fields.

**D-G / D-H split (Option β, decided 2026-05-27).** The original D-G+D-H plan (BoundaryFlux as a non-Field bundle; AngularFlux carrying boundary+history) was REVISED after the conversation surfaced the architectural conflation in AngularFlux's pre-D-H state. The new staging:
- **D-G**: BoundaryFlux becomes a pure Field; interior wavefront cache split to SweepScratch. AngularFlux still carries boundary+history (legacy preserved).
- **D-H.1**: TransportState container introduced; AngularFlux becomes pure Field; all consumers migrated. Legacy shims forward to hidden container.
- **D-H.2**: Legacy `AngularFlux.boundary` / `_history` deleted. User's hard constraint: by end of D-H, legacy is retired.

See §3.8 for the container architecture and §6 D-G/D-H step bodies for the per-step scope.


### Step D-A — L1 Field ABC (foundational)
- Create `orpheus/numerics/field.py` with the `Field` ABC per §3.2.
- Add unit tests in `tests/numerics/test_field.py` for the ABC algebra: same-space `+/-/*/neg`, cross-type rejection, cross-space rejection, inner_product, copy, linf, l2.
- Add `Field` to `orpheus/numerics/__init__.py` exports.
- **Verification**: `tests/numerics/test_field.py` PASS. P3.1 linter stays green. No production code changes (the ABC has no consumers yet).

### Step D-B — L1 `TensorProductSpace` primitive (LOAD-BEARING)

**Promoted from optional to load-bearing** (per architectural reassessment 2026-05-27, recorded in conversation). The grand report §15.1-15.2, §16A.10, and the north-star statement (line 5697) all frame neutron transport as "an algebra of ... tensor products ...". Today the codebase has shipped `TensorProductOperator`, `SumOfTensorProductsOperator`, and the `&` dunder at `orpheus/numerics/operator.py:1004-1216`, but ZERO production consumers exist. `TensorProductSpace` is the missing L1 primitive that types the *codomain* of these operators — without it, the operator algebra cannot express the §16A.10 BC tensor-network decomposition, the §15.2 scattering decomposition, or the §15.1 streaming decomposition. D-B unlocks D-B+1 (the first production tensor-network instance) and the entire Wave T that follows Depth B.

Scope:
- Add `FunctionSpace.__mul__(other) -> TensorProductSpace((self, other))`.
- Add `TensorProductSpace(FunctionSpace)` as a frozen dataclass with `factors: tuple[FunctionSpace, ...]`, derived `name`, `shape = tuple(chain.from_iterable(f.shape for f in factors))`, and an `inner_product` that respects per-factor weights along each factor's axis range.
- Add `FunctionSpace.dual()` → `DualSpace` (small additive — used by adjoint operator wiring).
- DEFER `__add__` (`DirectSumSpace`) per "unify after two instances" — lands with Phase 3.6 if a second use case appears.
- **Verification**: `tests/numerics/test_space.py` (or new `test_space_algebra.py`) extends with: associativity of `*`, shape composition law, inner-product factorisation along factor axes, dual idempotency. ~50 LoC of code, ~150 LoC of tests.

### Step D-B+1 — first production tensor-network instance (SPECULAR BC)

Rewires `SNBoundaryRealizer` for the specular case from `PermutationOperator(perm, axis=0)` to `PermutationOperator(perm) & IdentityOperator() & IdentityOperator()` — a 3-factor `TensorProductOperator` with the angular permutation as the only non-trivial factor (the `I_group ⊗ I_face` factors capture the implicit numpy broadcasting that the legacy single-axis form left untyped).

Scope:
- `orpheus/sn/boundary_realizer.py:149-156` — change the realised operator from `PermutationOperator(perm, axis=0)` (optionally wrapped in `ScaledOperator(albedo, ·)`) to `TensorProductOperator((P_angle, I_group, I_face))` (optionally `ScaledOperator(α, P_angle & I_group & I_face)`).
- `is_involution=True` is propagated as a per-factor property (the angular permutation factor IS the involution; the I factors are identity-involutive).
- Add postcondition assertion: the returned operator's `domain.shape == codomain.shape == boundary_trace_space.shape` (impossible to encode with the legacy single-axis form).
- **Verification**: bit-identical `.apply` (`TensorProductOperator` folds factors and `IdentityOperator.apply(x) == x`; the inner `np.take(x, perm, axis=0)` call is unchanged). `tests/sn/test_boundary_realizer.py` reflective-BC tests gain `isinstance(result, TensorProductOperator)` assertion and `result.assert_separable()` call.

D-B+1 ships INSTANCE #1 of the tensor-network pattern in production. INSTANCES #2-#5 (other BCs — vacuum, periodic, white, albedo) land in Wave T.1, with the abstraction settled before generalisation per `[[feedback_unify_after_two_instances]]`.

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

### Step D-G — Migrate `BoundaryFlux` to L2 as a pure `Field` (Option β stage 1)

**Decision recorded 2026-05-27** (in conversation): pursue Option β of three scoping
choices — stage the AngularFlux factoring across D-G (BoundaryFlux pure) + D-H
(AngularFlux pure + container + retire legacy). The user's hard constraint: **by
the end of D-H, the legacy `AngularFlux.boundary` / `AngularFlux._history`
structure is retired (deleted from the codebase after all gates green).** D-G is
the prep step: BoundaryFlux becomes a pure Field, decoupled from the interior
wavefront scratch it was conflating in 2-D.

Scope:

- Move `orpheus/sn/boundary_flux.py` → `orpheus/transport/fields/boundary_flux.py`.
- `BoundaryFlux` **inherits `Field`** — frozen dataclass with `values + space + mesh`,
  Field algebra inherited (no hand-coded dunders).
- Storage: **flat backing buffer + `FaceLayout` descriptor** (Option Ω from §3.4).
  `values: NDArray` is shape `(layout.total_size,)`; `FaceLayout` maps face names
  to slice + per-face shape. Per-face access is via slice views — no copies.
  This is the SAME storage shape that makes BoundaryFlux MoC-compatible and
  algebra-vectorizable; the difference from the pre-D-G §3.4 framing is that
  BoundaryFlux **is now a Field** (single flat space) rather than a non-Field
  bundle. The flat layout is what makes this principled.
- **Split out the interior wavefront cache** (the architectural problem the
  pre-D-G 2-D `xmin_xmax_buf` / `ymin_ymax_buf` buffers were conflating with
  boundary trace state). Introduce a `SweepScratch` type (L3, SN-specific)
  carried by the sweep itself, NOT by BoundaryFlux. BoundaryFlux becomes
  IMMUTABLE; the mutable wavefront cache lives on the sweep.
- BoundaryFlux is `frozen=True` and `kw_only=True`. The pre-D-G mutability
  (write-through cache for sweep persistent BC state) is replaced by sweep-side
  cache + functional reconstruction at iteration boundaries.
- Ship a `BoundaryFlux.zeros_for_sn_mesh(mesh)` classmethod factory (the
  geometry-conditional construction logic lives here, not at sweep call sites).
- AngularFlux still carries `boundary: BoundaryFlux` and `_history` in D-G (legacy
  preservation; the retirement is D-H's responsibility). The pre-D-G mutable
  write-through path stays alive via a one-cycle shim until D-H lands.
- Ship one-cycle re-export shim at `orpheus/sn/boundary_flux.py`.
- Migrate consumers (~15 files per audit memo §3 / §4). The sweep gets a
  `SweepScratch` parameter; BC state reads/writes use the BoundaryFlux Field
  interface.

**Verification**:

- Every SN test stays green. STRICT bit-identical on the BoundaryFlux algebra
  itself (the data hasn't changed, the type discipline has).
- New tests: `tests/transport/fields/test_boundary_flux.py` — pin Field algebra,
  per-face slice views, mesh-binding rejection.
- New test: `tests/sn/test_sweep_scratch_split.py` — pin that 2-D sweeps still
  produce bit-identical interior fluxes after the buffer split.
- The 10 pre-existing DD-regression failures stay AT the same failure set.

### Step D-H — Pure `AngularFlux` Field + `TransportState` container + RETIRE legacy (Option β stage 2)

**This is the load-bearing step of Option β.** It must:

1. Refactor `AngularFlux` to a PURE `Field` (drop `boundary`, drop `_history`,
   drop `history_depth`). Symmetric with `ScalarFlux`.
2. Introduce the **`TransportState`** container (working name; see §3.8 for the
   naming discussion — alternates: `TimedFullField` per the user's proposal,
   `IterationState`, `FullTransportField`). The container holds
   `(psi: AngularFlux, boundary: BoundaryFlux, _history: tuple, history_depth: int)`
   and exposes dunders that delegate to the contained Fields.
3. Migrate ALL consumers from `AngularFlux.boundary` / `AngularFlux._history`
   accesses to the container.
4. **Retire** the legacy `AngularFlux.boundary` / `AngularFlux._history` fields
   from the codebase. No shim survives D-H.

Architecture details live in **§3.8 — Container architecture (`TransportState`)**.
The summary:

- L2 `AngularFlux(Field)`: `values + space + mesh`. No boundary, no history.
  Inherits all Field algebra. From-mesh factory + `from_ndarray` test ergonomic.
- L2 `TransportState`: NOT a Field (structured bundle). Holds
  `psi: AngularFlux + boundary: BoundaryFlux + _history: tuple[TransportState, ...] + history_depth: int`.
  Algebra: `__add__`, `__sub__`, `__mul__`, `__neg__`, `__truediv__` propagate
  to inner fields; algebra results carry empty history (history is iteration
  metadata, not algebraic state).
- L2 `TransportState.psi << new_psi` (or `.advance(new_psi, new_boundary)`)
  rotates the history shift register — the SAME pattern the pre-D-H
  `AngularFlux._history` exposed, lifted to the container.
- L3 `from_flat_with_traces` / `to_flat_with_traces` — installed at IMPORT time
  on `TransportState` by `orpheus/sn/transport_state_b1pp_adapter.py`. The
  `_is_ravellable` Protocol at `iteration.py:163-176` now targets `TransportState`
  (the container) — the wrong thing-to-flatten target is now structurally
  unrepresentable.
- L3 `TransportState.zeros_for_sn_mesh(mesh)` — the canonical constructor.
  Constructs zero psi, zero boundary, empty history.
- All sweep / SourceIteration / KrylovAcceleration / InvertibleOperator.solve
  consumers migrated to consume/produce `TransportState`.

**Scope of changes (file-by-file estimate, to be refined post-D-G)**:

| File | Touch | Reason |
|---|---|---|
| `orpheus/transport/fields/angular_flux.py` | NEW | pure Field subclass |
| `orpheus/transport/state.py` (or `transport/transport_state.py`) | NEW | container class |
| `orpheus/sn/transport_state_b1pp_adapter.py` | NEW | L3 ravel/unravel injection |
| `orpheus/sn/angular_flux.py` | RETIRE | one-cycle shim then delete |
| `orpheus/numerics/iteration.py` | MOD | `_is_ravellable` targets container |
| `orpheus/numerics/krylov.py` (or wherever Krylov adapter lives) | MOD | packed-vector signature |
| `orpheus/sn/source_iteration.py` | MOD | consume/produce `TransportState` |
| `orpheus/sn/operator.py` (InvertibleOperator.solve) | MOD | R-1 Step C signature |
| `orpheus/sn/sweep.py` | MOD | sweep API targets container (psi out, boundary out) |
| ~25 other consumer files | MOD | per audit memo §3 / §4 |
| `tests/sn/...` | MIGRATE | `psi.boundary` → `state.boundary`; `psi._history` → `state._history` |

Estimated 600-900 LOC across orpheus + tests. Single commit (per the user's
"complete migration AND retire legacy" constraint) OR split into D-H.1 (container
introduced, dual API both ways) + D-H.2 (legacy retired) — see §6.X decision
point below.

**Container space discipline**:

`TransportState` does NOT inherit from `Field` in D-H. Per `feedback_unify_after_two_instances`,
the natural backing — a `DirectSumSpace` Field with `psi.space ⊕ boundary.space`
— remains DEFERRED to Phase 3.6, when kinetics' `flux ⊕ precursors` lands the
second use case. D-H ships the container as a structured bundle with
delegate-style dunders. The Phase 3.6 promotion (`TransportState` becoming
`Field(values=DirectSumStorage, space=DirectSumSpace(psi.space, boundary.space))`)
is a non-breaking refactor — the public API stays the same; only the storage
backing changes.

**Decision point at D-H planning (open)**:

- **D-H.1 + D-H.2 split**: D-H.1 introduces `TransportState` and migrates all
  CONSUMERS to use it; the legacy `AngularFlux.boundary` / `_history` stay alive
  as shims that just forward to a hidden container. D-H.2 retires the legacy.
  Two commits, both green.
- **D-H single commit**: do everything in one commit, no intermediate dual-API
  state. Higher review cognitive load. One bisect point.

User constraint: "by the end of D-H, the legacy is retired." Both forms satisfy
the constraint. **Recommended: D-H.1 + D-H.2 split** — keeps each commit small
enough to review and bisectable. Decided post-D-G when the actual consumer-site
count is verified.

**Verification**:

- Every L1 MMS gate stays green. PRINCIPLED bit-identity (algorithm unchanged;
  algebra path is delegated through container; ULP-level drift acceptable iff
  bounded by `iter_count × ULP`).
- The 10 pre-existing DD-regression failures stay AT the same failure set.
- New test: `tests/transport/test_transport_state.py` — pin container algebra,
  history shift-register, dunder propagation, mesh-binding rejection,
  cross-container rejection.
- New test: `tests/sn/test_transport_state_b1pp_adapter.py` — pin
  `from_flat_with_traces` / `to_flat_with_traces` round-trip on the container.
- New test: `tests/numerics/test_ravellable_protocol_targets_container.py` —
  pin that `_is_ravellable` detects `TransportState`, NOT `AngularFlux`.
- The `_is_ravellable` Protocol at `iteration.py:163-176` STILL works (target
  changed; mechanism unchanged).

**Test-architect dispatch BEFORE D-H implementation** — per the proactive-trigger
table in `subagent-handoff-protocol`. The D-H carve crosses the typed↔packed
boundary AND the field-singleton↔composite-state boundary. The proactive
dispatch produces the verification spec: which existing tests pin the legacy
behaviour (and must migrate to the container API), which new tests catch the
container-side convention, which L1 MMS gates must continue passing across the
swap.

### Step D-I — Wire operators to typed Fields (keep singledispatch, retire ndarray handlers)

Per §3.7 policy, the `@singledispatchmethod` dispatch tables on `ScatteringOperator.apply` and `FissionOperator.apply` STAY — they're the canonical storage-axis polymorphism pattern. Only the **ndarray handlers** retire.

- **`MomentProjection.apply`** — add typed handler `apply(AngularFlux) -> HarmonicMomentField`. Codomain check via `psi.space.shape[axes] == self.domain.shape`. Retire the bare-ndarray entry point.
- **`MomentProjection.apply_transpose(c: HarmonicMomentField) -> AngularFlux`** — typed analog.
- **`LegendreMomentScattering.apply`** — refactor the `ndarray | HarmonicMomentField` union into a `@singledispatchmethod` with ONE handler (`apply(HarmonicMomentField)`). Retire the ndarray branch.
- **`ScatteringOperator.apply`** — KEEP the dispatch table. KEEP `apply(AngularFlux)` and `apply(ScalarFlux)` handlers. **Retire only the `apply(np.ndarray)` handler.** Tests migrate to `AngularFlux.from_ndarray(...)` / `ScalarFlux.from_ndarray(...)` + typed dispatch.
- **`FissionOperator.apply`** — same pattern as ScatteringOperator. Keep dispatch table; retire ndarray handler.
- **`StreamingOperator.apply` (modern, `operator.py:1801`)** — already accepts typed `AngularFlux`. The bare-ndarray branch is the L4 GMRES adapter primitive and is out of scope for Depth B (separate GMRES adapter cleanup wave). `SNStreamingOperator` (legacy `L + C` bundle, `operator.py:1329`) is on the retirement queue and is NOT a Depth B target — the modern `(L + C).apply` via `InvertibleOperator(OperatorSum)` is the algebraic path.
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
- `tests/transport/fields/test_boundary_flux.py` — migrated, exercises pure-Field BoundaryFlux algebra + per-face slice views (D-G).
- `tests/sn/test_sweep_scratch_split.py` — pins that 2-D sweeps produce bit-identical interior fluxes after the BoundaryFlux/SweepScratch buffer split (D-G).
- `tests/transport/fields/test_angular_flux.py` — migrated, exercises pure-Field AngularFlux algebra (no boundary, no history) (D-H).
- `tests/transport/test_transport_state.py` — NEW. Pins container algebra (delegate dunders), history shift-register (`advance`, `at_lag`), mesh-binding rejection, cross-container rejection (D-H).
- `tests/sn/test_transport_state_b1pp_adapter.py` — pins the SN-adapter's `from_flat_with_traces` / `to_flat_with_traces` round-trip on `TransportState` (D-H).
- `tests/numerics/test_ravellable_protocol_targets_container.py` — pins that `_is_ravellable` detects `TransportState`, NOT `AngularFlux` (D-H).

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
| D-G | PRINCIPLED | BoundaryFlux becomes a pure Field; flat-buffer layout is unchanged but the access surface flips from mutable write-through to immutable functional. Interior wavefront cache splits to SweepScratch — same numerical content, different ownership. L1 MMS gates pin numerical equivalence; ULP drift bounded by `(reduction depth) × ULP` per `vv-principles`. |
| D-H.1 | PRINCIPLED | TransportState container introduced; consumers migrated. AngularFlux algebra preserved; container delegate dunders compose to identical numerical paths. ULP drift bounded by `(reduction depth) × ULP`. |
| D-H.2 | N/A | Legacy retirement (deletion only — no new numerical paths). |
| D-I | PRINCIPLED | Operator signatures change from `ndarray` to `Field`. Underlying algorithm unchanged. |
| D-J | N/A | Dead-code retirement. |

The PRINCIPLED step (D-I) needs the three criteria:
1. **Principled**: every intermediate is a named typed Field. SATISFIED by design.
2. **Structurally-independent reference**: the L1 MMS gates remain the truth-of-record. SATISFIED.
3. **Drift dimensionally explainable**: zero drift expected (algorithm unchanged); ULP-level drift would be the only possible signal. If snapshots break by more than `iter_count × ULP`, the rewire is wrong.

### 7.4 Same-class invariant test (target: TransportState post-D-H)

Add `tests/numerics/test_ravellable_protocol_same_class.py`:

```python
def test_ravellable_class_identity_preserved():
    """Pin: TransportState returned by from_flat_with_traces is the SAME
    class as the input template. Breaks _is_ravellable Protocol detection
    if violated."""
    # Construct an L2 TransportState via the SN adapter
    state = TransportState.zeros_for_sn_mesh(some_sn_mesh)
    assert type(state).__name__ == "TransportState"
    flat = state.to_flat_with_traces()
    reconstructed = type(state).from_flat_with_traces(flat, state)
    assert type(reconstructed) is type(state)
    # Inner Fields also preserve type
    assert type(reconstructed.psi) is type(state.psi)
    assert type(reconstructed.boundary) is type(state.boundary)
```

This pins the contract that any future L2/L3 split (e.g., for MoC) MUST preserve.

**Pre-D-H (during D-G): the protocol target is still `AngularFlux`.** The legacy `_is_ravellable` continues to detect `AngularFlux` until D-H.1 migrates the target. The test above is added in D-H.1 alongside the migration; the pre-D-H legacy test (`tests/numerics/test_ravellable_protocol_angular_flux.py`) is retired in D-H.2 along with the legacy `AngularFlux.boundary` / `_history` fields.

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
| **3** | `BoundaryFlux` pure-Field carve + interior-cache split (D-G) breaks consumer sites. | Audit-driven migration; the audit memo identified ~15 files. One-cycle deprecation shim on legacy mutable write-through; migrate consumers in lockstep. SweepScratch isolation tested by new `tests/sn/test_sweep_scratch_split.py`. |
| **3b** | `TransportState` container migration (D-H.1) is the largest single carve in Depth B — ~30 consumer files, crosses typed↔packed AND field↔composite boundaries. | Stage as D-H.1 (container in + consumers migrated, legacy shims forward) → D-H.2 (legacy retired). Each commit independently green + bisectable. Proactive test-architect dispatch BEFORE D-H.1 implementation (per `subagent-handoff-protocol` triggers). Krylov adapter `_is_ravellable` Protocol migration is the highest-risk single touch. |
| **3c** | Legacy retirement deadline at D-H.2 missed — `AngularFlux.boundary` / `_history` shims linger past D-H. | User constraint is hard: "by the end of D-H, the legacy structure is retired (deleted from the codebase after all gates green)." D-H.2 IS the retirement commit; if blocked, the blocker (e.g., a consumer migration that didn't complete in D-H.1) is fixed before D-H.2 lands. NEVER bundle D-H.2 into a later step. |
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
3. **`orpheus/transport/` exists** as the L2 package with `fields/` (AngularFlux, ScalarFlux, HarmonicMomentField, BoundaryFlux, BoundaryFaceFlux), `sources/` (IsotropicSource, PerOrdinateSource), and the `TransportState` container (`orpheus/transport/state.py` or `transport_state.py`).
4. **AngularFlux is a pure Field** (`values + space + mesh`). It does NOT carry `boundary` or `_history` — those responsibilities live on `TransportState`. Symmetric with `ScalarFlux`. The legacy `AngularFlux.boundary` / `_history` fields are DELETED from the codebase (retired in D-H.2).
5. **BoundaryFlux is a pure Field** (flat backing buffer + `FaceLayout` descriptor + `mesh`). IMMUTABLE. The pre-D-G mutable write-through + interior wavefront cache is split: BoundaryFlux carries only trace state; `SweepScratch` (L3, SN-private) carries the wavefront cache.
6. **The L3 SN adapter installs `from_flat_with_traces` / `to_flat_with_traces`** on `TransportState` (NOT on AngularFlux) at module-import time. The same-class invariant test (`tests/numerics/test_ravellable_protocol_same_class.py`) is GREEN with TransportState as the target.
7. **The `_is_ravellable` Protocol at `iteration.py:163-176` detects `TransportState`** — the Krylov adapter packs/unpacks the composite state, NOT the inner flux alone.
8. **`@singledispatchmethod` dispatch tables stay** on `ScatteringOperator.apply`, `FissionOperator.apply` with handlers for `AngularFlux` and `ScalarFlux` (and future `AngularSource`/`AngularResidual` when #201 lands). The `apply(np.ndarray)` handlers are RETIRED. Each typed field has a `.from_ndarray(arr, mesh)` factory.
9. **`SourceIteration.__init__` / `KrylovAcceleration.__init__` do Layer 2 dimensional checks** at construction. Production runs (`python -O`) get this safety net. Sweep + `(L+C).solve` (InvertibleOperator R-1 Step C) consume and produce `TransportState`.
10. **All four completed Phase 3 steps** (P3.1 import-linter, P3.5 range→codomain, P3.0 layering docs, P3.2 spherical_harmonics shim retirement) remain GREEN. The 10 pre-existing DD-regression failures stay at the same failure set.
11. **`pint >= 0.20`** is in `pyproject.toml`. CI matrix runs both `pytest` and `python -O -m pytest`.
12. **GitHub Issues #201 and #205 remain OPEN** — they are downstream consumers of Depth B, not closed by it.
13. **`DirectSumSpace` is still DEFERRED to Phase 3.6** per `feedback_unify_after_two_instances`. The TransportState container ships as a structured bundle (delegate dunders); Phase 3.6's kinetics flux⊕precursors lands DirectSumSpace and may then promote TransportState to Field-with-DirectSumSpace backing (non-breaking).

### 11.2 What P3.4 picks up

P3.4 (parent plan §P3.4, around line 867 of `moment_space_and_layering_plan.md`) introduces:
- `Problem` ABC at L2 (`orpheus/transport/problems/`) — declarative descriptions: `CriticalityProblem`, `FixedSourceProblem`, `AlphaEigenproblem`, `InitialValueProblem`
- `Solver` ABC at L1 (`orpheus/numerics/solvers/`) — iterative algorithms: `PowerIteration` (formerly `KEigenvalue`), `Arnoldi`, `SourceIteration`, time-steppers
- `homogeneous/solver.py:26`'s `power_iteration` legacy retires (CC.8 close-out)

P3.4 directly consumes Depth B's typed Fields and the `TransportState` container:
- `CriticalityProblem` carries `loss: LinearOperator` (codomain = AngularResidual space) and `fission: LinearOperator` (codomain = AngularSource space)
- `PowerIteration.solve(problem) -> TransportState` — the composite state (psi + boundary + history) is the result; consumers access `.psi` for the flux, `.boundary` for the trace
- `FixedSourceProblem` carries `(L, S, F, q_ext: AngularSource | TransportState)` — every component is dimensionally explicit; the rhs may be a pure source or a state-shaped rhs depending on the formulation
- `SourceIteration.step(state: TransportState, source) -> TransportState` — the iteration verb consumes and produces the container

The Layer 2 construction-time dimensional check (Depth B §5) is the FIRST check `PowerIteration` does on its problem — verifying the operator algebra makes dimensional sense before any iteration runs.

### 11.3 Sequence after Depth B — Wave T precedes P3.4

The parent plan's REVISED Phase 3 sequencing was extended on 2026-05-27 to insert **Wave T (Tensor-Network Operator Algebra)** between Depth B and P3.4:

```text
P3.1 ✓ → P3.5 ✓ → P3.0 ✓ → P3.2 ✓ → [Depth B = P3.3] → Wave T → P3.4 → P3.6
                                       ↑                  ↑       ↑     ↑
                                       THIS PLAN          NEW     NEXT  LAST
```

**Wave T** is the load-bearing consumer of the `TensorProductOperator` / `SumOfTensorProductsOperator` / `TensorProductSpace` infrastructure shipped here (D-B) and in `numerics/operator.py:1004-1216`. Today that infrastructure has ZERO production consumers; Wave T rewires the BC realizers, fission, scattering, and modern `StreamingOperator.apply` to use the algebra natively. Detailed plan: `.claude/plans/wave_t_tensor_network.md`. Wave T substeps (concrete forms):

- **T.1** — Remaining BC realizers (vacuum, periodic, white, albedo) as `K_factor ⊗ I_group ⊗ I_face` (extends D-B+1 specular).
- **T.2** — Fission `F = χ ⊗ νΣ_f` as rank-1 `TensorProductOperator`.
- **T.3** — Scattering `Σ_ℓ Σ_ℓ ⊗ A_ℓ ⊗ G_ℓ` per §15.2 as `SumOfTensorProductsOperator`.
- **T.4** — `StreamingOperator.apply` as `L_spatial + L_angular_redist` per the connection-coefficient framing (`geometry/reduced_operator.py:1-30`). Universal across slab/sphere/cylinder; slab degenerates to `ZeroOperator` for `L_angular_redist`.

Wave T's dependency: Depth B complete (typed Fields available as operator domain/codomain). Wave T's exit: P3.4 (Problem/Solver split). T.4 may surface complications when curvilinear specifics resist clean factoring — face difficulties as they come; the architectural commitment to one algebraic form across geometries is non-negotiable.

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
