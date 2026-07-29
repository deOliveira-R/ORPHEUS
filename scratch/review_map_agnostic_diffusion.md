# Review map — model-agnostic machinery + the DIFFUSION consumer

Grounded against `main` @ `8654d348`, Nexus graph fresh (18312 nodes / 178361 edges).
Line numbers are current at this read; re-derive via Nexus if the tree moved.

**Scope note:** `transport/operators/fission.py`, `scattering.py`, `sn/solver.py`,
`sn/coupled_system.py`, `numerics/operator.py`, `numerics/eigenvalue.py` are owned by
sibling explorers and are only touched here where MY surface consumes them.

---

## §1. `orpheus/transport/` inventory

### 1.0 The layering invariant — VERIFIED CLEAN

`grep -rn "^\s*from orpheus\.sn\|^\s*import orpheus\.sn" orpheus/transport/ orpheus/numerics/`
returns **11 hits, ALL of them inside an `if TYPE_CHECKING:` block**. There is **zero
runtime import of a method package from `transport/` or `numerics/`**:

| file:line | symbol | guard |
|---|---|---|
| `transport/radial_characteristic_field.py:80` | `SNMesh` | `TYPE_CHECKING` (`:77`) |
| `transport/operators/scattering.py:91` | `SNMesh` | `TYPE_CHECKING` (`:90`) |
| `transport/source_sinks/scalar_source_sink.py:81` | `SNMesh` | `TYPE_CHECKING` (`:80`) |
| `transport/source_sinks/angular_source_sink.py:70` | `SNMesh` | `TYPE_CHECKING` (`:69`) |
| `transport/source_sinks/angular_boundary_source_sink.py:136` | `SNMesh` | `TYPE_CHECKING` (`:131`) |
| `transport/fields/angular_flux.py:53` | `SNMesh` | `TYPE_CHECKING` (`:52`) |
| `transport/fields/harmonic_moment_flux.py:113` | `SNMesh` | `TYPE_CHECKING` (`:112`) |
| `transport/fields/_bases.py:116-118` | `DiffusionMesh` + `FaceLayout` + `SNMesh` | `TYPE_CHECKING` (`:115`) |

All remaining ~120 `orpheus.sn` occurrences in `transport/`+`numerics/` are **Sphinx
docstring cross-references** (`:class:`/`:func:`/`:meth:`/`:mod:`), not imports.
Per lesson L33 / Operating Principle: a docstring is a CLAIM — but here the claim
("transport is method-neutral") is corroborated by the absence of runtime edges.

**Caveat (the honest reading):** the *type annotations* are SN-shaped even though the
imports are deferred. `AngularFlux.from_mesh(mesh: "SNMesh")`-style signatures mean the
STATIC contract of several transport types names `SNMesh` specifically, not
`TransportMethod`. That is a typing-surface coupling, not an import-cycle coupling.
`transport/fields/_bases.py:116-118` is the one place that names BOTH witnesses
(`DiffusionMesh` *and* `SNMesh`), i.e. the one genuinely method-generic annotation.

### 1.1 Module-by-module inventory

`orpheus/transport/` = 55 modules / ~14.9 kloc. "Neutral?" column: **Y** = no angular /
ordinate / quadrature / sweep concept anywhere in the body; **A** = angular-shaped by
DESIGN (it *is* the angular representation — neutral in the sense that it lives in
transport because SN/MOC/PN all need it, but a diffusion consumer never touches it);
**S** = secretly SN-shaped (a genuinely method-neutral concept that names an SN artifact).

#### Top level

| module | owns | neutral? |
|---|---|---|
| `__init__.py` (60) | package docstring; re-exports nothing structural | Y |
| `method.py` (311) | `TransportMethod` Protocol + `resolve_boundary_conditions` + `_law_from_tag` | **Y** — the method-generic seam; imports only `geometry` + `numerics` + `transport.mesh.axis` |
| `full_field.py` (542) | `FullField` = `interior ⊕ boundary` composite; `to_flat`/`from_flat`/`zeros` | **Y** — generic in the two component types; diffusion instantiates `FullField(interior=ScalarFlux, boundary=ScalarBoundaryFlux)` |
| `timed_full_field.py` (390) | `TimedFullField` — `FullField` + a time stamp | Y |
| `radial_characteristic_field.py` (316) | the curvilinear ψ½ half-level carrier | **A** (SN curvilinear-specific in practice) |
| `reaction_rate_functional.py` (231) | `IntegratedReactionRate` — `∫⟨Σx, φ⟩dV` | **Y** — consumed by BOTH `diffusion/solver.py:302` and `homogeneous/solver.py:210-211` |

#### `fields/` — the FLUX-role leaf types

| module | owns | neutral? |
|---|---|---|
| `_bases.py` (1268) | the whole field ABC ladder (see §3) | Y-ish — the ONE module naming BOTH witnesses (`DiffusionMesh` + `SNMesh`, `:116-118`) |
| `_flux_role.py` (206) | `FluxRole` enum — the flux/source/residual/displacement role tag | Y |
| `_coefficient_role.py` (84) | `CoefficientRole` — the XS-field role tag | Y |
| `scalar_flux.py` (153) | `ScalarFlux` `(ng, *spatial)` | **Y** — diffusion's bulk carrier |
| `scalar_boundary_flux.py` (143) | `ScalarBoundaryFlux` — the `(J⁺, J⁻)` trace | **Y** — diffusion's trace carrier |
| `angular_flux.py` (131) | `AngularFlux` `(N, ng, *spatial)` | **A** |
| `angular_boundary_flux.py` (162) | `AngularBoundaryFlux` — per-face angular trace | **A** |
| `harmonic_moment_flux.py` (279) | `HarmonicMomentFlux` — SH moments | **A** (PN/SN shared) |
| `cross_section_field.py` (109) | `CrossSectionField` — a typed XS carrier | Y |
| `radial_characteristic_*_flux.py` (52+53) | curvilinear half-level flux/trace | **A** |

#### `source_sinks/` — the SOURCE-role leaf types (the plan's premise is STALE here, see §3)

| module | owns |
|---|---|
| `scalar_source_sink.py` (207) | **`ScalarSourceSink` — ALREADY EXISTS** |
| `scalar_boundary_source_sink.py` (94) | `ScalarBoundarySourceSink` |
| `harmonic_moment_source_sink.py` (147) | **`HarmonicMomentSourceSink` — ALREADY EXISTS** |
| `angular_source_sink.py` (198) | `AngularSourceSink` |
| `angular_boundary_source_sink.py` (295) | `AngularBoundarySourceSink` |
| `radial_characteristic_{interior,boundary}_source_sink.py` (64+67) | curvilinear source arms |

#### `residuals/` + `displacements/` — the other two role families

`residuals/`: `ScalarResidual` (127), `AngularResidual` (206), `AngularBoundaryResidual`
(147), `RadialCharacteristic{Interior,Boundary}Residual` (84+71).
`displacements/`: `_displacement.py` (195) + `Scalar` (32), `ScalarBoundary` (38),
`Angular` (69), `AngularBoundary` (38), `Moment` (67), `RadialCharacteristic{Interior,Boundary}` (47+48).
Both families are **role-complete on the SCALAR locus** already.

#### `mesh/` — the method-agnostic data carrier

| module | owns | neutral? |
|---|---|---|
| `axis.py` (688) | `Axis1D`, `AXIS_NAMES`, `FaceLabel`, `face_labels()`, `n_unknowns_flat` | **Y** — the face inventory both witnesses derive from |
| `material_mesh.py` (452) | `MaterialMesh` = mesh + materials DATA; `from_materials` (meshless), `material_xs_field()` | **Y** — the base class BOTH `SNMesh` and `DiffusionMesh` extend |
| `material_xs_field.py` (1229) | `MaterialXSField` — per-cell gathered XS fields incl. `diffusion_coefficient` | **Y** |

#### `spatial/` — the spatial differencing schemes

| module | owns | neutral? |
|---|---|---|
| `scheme.py` (1226) | the `SpatialScheme` family + `CellVisit` | **S** — promoted out of `orpheus.sn.sweep` at #272; `CellVisit`/`τ`/`c` closure weights and `pole_angular_closure` references are SN-sweep vocabulary |
| `diamond.py` (781) | WDD/diamond coefficients | **S** |
| `linear_discontinuous.py` (835) | LD scheme | **S** |
| `cell_balance.py` (371) | the cell balance relation | **S** |
| `_ubld.py` (681) | upstream-corner balance LD | **S** |

`transport/spatial/` is the **least method-neutral** part of `transport/`: it lives at L2
but its whole vocabulary (`CellVisit`, octant signs, `pole_angular_closure`, `CollisionCache`,
`SweepDependencyGraph.walk_full`) is the SN sweep's. Diffusion does NOT consume any of it —
`diffusion/operators.py` writes its own FD stencil (`_interior_conductance`,
`_boundary_closure`, `:170-218`). This is a *documented* promotion (`spatial/__init__.py:12`),
not an accident, but the honest reading is: 4.7 kloc of `transport/spatial/` has exactly ONE
consumer family (SN), so its residence at L2 is aspirational rather than earned.

#### `operators/` — the shared reaction/loss algebra

| module | owns | consumed by |
|---|---|---|
| `multiplication_operator.py` (530) | `MultiplicationOperator` — `M[σ]` diagonal | **SN + diffusion + homogeneous** |
| `isotropic_scattering.py` (438) | `IsotropicScattering` (Σ_s0ᵀ) + `IsotropicN2N` (2Σ₂ᵀ) | **diffusion + homogeneous** (NOT SN) |
| `fission.py` (674) | `FissionOperator` — the rank-1 χ⊗νΣf dyad | **SN + diffusion + homogeneous** — explorer (a)'s surface |
| `scattering.py` (1326) | the anisotropic/moment scattering family | **SN** — explorer (a)'s surface |
| `integral_kernel_operator.py` (211) | integral-kernel base (CP/MOC-facing) | Y |

#### `frames/`

`harmonic_frame.py` (208) — the SH `HarmonicFrame` (analysis/reconstruction over
`SphericalHarmonicSpace`). **A** by design.

> Correction to §1.1: `fields/_flux_role.py` is a **mixin class**, not an enum
> (`FluxRole`, `_flux_role.py:59`) — it supplies the affine/torsor `__add__`/`__sub__`
> algebra, not a role tag. Same for `_coefficient_role.py:73` (`CoefficientRole`).
> **The role is carried by the CLASS, never by a tag field** — see §3.

---

## §2. `TransportMethod` — verified current shape

`orpheus/transport/method.py`. The memory's claim (minted over `{SNMesh, DiffusionMesh}`
when a realizer registry was dissolved, #290 P7b @ `44d583e`) is **CONFIRMED against the
tree**, with detail:

- **What it is:** `@runtime_checkable class TransportMethod(Protocol[OpT_co])`
  (`method.py:145-146`). A **generic structural Protocol**, NOT an ABC, NOT a dataclass.
  Two TypeVars: `OpT_co` (covariant, bound `LinearOperator`, `:137`) for the Protocol and
  the invariant twin `OpT` (`:142`) for the generic function.
- **Interface — exactly 4 members** (`:181-215`):
  1. `axes: tuple[Axis1D, ...]` — the `MaterialMesh` data block;
  2. `BOUNDARY_OPERATOR_REGISTRY: ClassVar[dict[str, type[BoundaryTraceLaw]]]` — the
     per-method law-admission table;
  3. `bc: Mapping[str, OpT_co]` (a property) — the realized per-face laws;
  4. `realize_boundary_law(law, face) -> OpT_co` — the per-method arm.
- **The ONE generic body:** `resolve_boundary_conditions(method) -> dict[str, OpT]`
  (`:218-262`) + the tag-parse helper `_law_from_tag` (`:265-311`). The loop is
  `for label in face_labels(method.axes)` → default `BC("reflective")` → parse → call the
  method's arm. `_law_from_tag` special-cases exactly two law families
  (`ReflectiveBoundary` needs the face's axis, `AlbedoBoundary` needs `params["albedo"]`)
  and calls `law_cls()` for everything else.
- **Conformers (2, both structural — nobody inherits or imports the Protocol):**
  - `orpheus/sn/mesh/augmented_mesh.py` `SNMesh` — `OpT` = the kind-tagged
    `_BoundBoundaryOperator` shim;
  - `orpheus/diffusion/augmented_mesh.py:122` `DiffusionMesh(MaterialMesh)` — `OpT` = a
    bare `LinearOperator`. Its `realize_boundary_law` is at `augmented_mesh.py:294-315`
    (builds `DiffusionMethodSpace.for_face` → `DiffusionBoundaryRealizer().realize`);
    its registry at `:158-163` = `{vacuum, reflective, albedo, zero_flux}` (`"white"`
    deliberately absent — coincides with reflective at P1).
  - Confirmed: `resolve_boundary_conditions(self)` is called at
    `diffusion/augmented_mesh.py:247-249`. Both method packages import the FUNCTION,
    never the Protocol class (`method.py:60-73` states this; grep confirms — the only
    `TransportMethod` name occurrences outside `method.py` are docstrings/tests).
- **Deliberate exclusions, documented at `:34-58`:** promotion classmethods (signatures
  genuinely differ — SN injects a quadrature, diffusion injects nothing), method-space
  construction, the trace spaces (`angular_trace` vs `scalar_trace` are locus SIBLINGS),
  and any `method_name` dispatch string.
- **The one honest gap (`:162-175`):** `full_field_space` — a UNIVERSAL member both
  witnesses implement — is deliberately NOT declared, because both implement it as a
  `functools.cached_property` which current pyright won't match against a Protocol
  property. The docstring also records that the anticipated first generic consumer
  ("the DSA driver, #2") **did NOT materialize**: the R4 ruling (2026-07-26) wired
  consistent DSA through an SN-side edge-centered low-order system consuming `SNMesh`
  directly, so **no `DiffusionMesh` enters the accelerated loop**. ⇒ `TransportMethod`
  today has ZERO generic production consumers beyond `resolve_boundary_conditions`.

---

## §3. Field / space type hierarchy — and the plan's STALE premise

### 3.1 The three-axis grid, and where the ROLE lives

`orpheus/transport/fields/_bases.py:12-60` states the design: **locus** × **family** ×
**role**. Verified against the class declarations:

```
Field                       (orpheus/numerics/field.py — values + space + dunder algebra)
├─ BulkField (_bases.py:148)          codim-0; mesh: MaterialMesh; abstract _phase_space_shape
│   ├─ AngularField (:328)   mesh NARROWED to SNMesh (:344); _SPACE_NAME ClassVar
│   ├─ ScalarField  (:463)   mesh stays MaterialMesh;      _SPACE_NAME ClassVar
│   └─ MomentField  (:552)   mesh NARROWED to SNMesh (:604); _CELL_GROUP_NAME ClassVar; carries L
└─ FaceField[K] (:737)                codim-1; flat single buffer + FaceLayout[K]
    ├─ BoundaryField (:934, FaceField[str])   the FullField boundary-slot discriminator
    │   ├─ AngularBoundaryField (:1015)  mesh: SNMesh;        space: AngularTraceSpace
    │   └─ ScalarBoundaryField  (:1070)  mesh: DiffusionMesh; space: ScalarTraceSpace
    ├─ RadialCharacteristicInteriorField (:1147, FaceField[tuple[int,int]])
    └─ RadialCharacteristicBoundaryField (:1208, FaceField[tuple[int,int]])
```

**The ROLE is carried by the CLASS, not by a field or flag.** Evidence: each role leaf is
its own `@dataclass` subclass of the storage ABC, and the arithmetic gate is
class-identity (`Field._check_partner`, extended at `_bases.py:198` bulk / `:835` face).
Two mechanisms encode the role:
- **the class**, e.g. `class ScalarFlux(FluxRole, ScalarField)` (`scalar_flux.py:104`) vs
  `class ScalarSourceSink(ScalarField)` (`scalar_source_sink.py:89`) — the flux mixes in
  `FluxRole` (affine/torsor algebra: `flux+flux` is a TypeError, `flux−flux→displacement`,
  `_flux_role.py:96-122`), the source does NOT (it is a plain vector-space rate density);
- **the space NAME**, via the `_SPACE_NAME` / `_CELL_GROUP_NAME` ClassVar:
  `"scalar_flux"` (`scalar_flux.py:142`), `"scalar_source_sink"`
  (`scalar_source_sink.py:121`), `"scalar_residual"` (`scalar_residual.py:105`),
  `"angular_flux"`/`"angular_source_sink"`/`"angular_residual"`, `"cross_section"`,
  `"cell_group"` (`_bases.py:612`) vs `"cell_group_source_sink"`
  (`harmonic_moment_source_sink.py:140`).
- plus a `UNITS: ClassVar[Unit]` per leaf (`SCALAR_RATE_UNITS` /`ANGULAR_RATE_UNITS`
  for bulk sources, `*_FLUX_UNITS` for traces).

Displacements deliberately REUSE the flux's space name (`scalar_displacement.py:29` =
`"scalar_flux"`, `angular_displacement.py:43` = `"angular_flux"`) — the tangent vector
lives in the flux's own function space; only the CLASS is the role gate
(`_flux_role.py:39-41`).

### 3.2 Full role × locus × family census — WHAT EXISTS TODAY

| | Flux | SourceSink | Residual | Displacement |
|---|---|---|---|---|
| Scalar (bulk) | `ScalarFlux` | **`ScalarSourceSink`** | `ScalarResidual` | `ScalarDisplacement` |
| Angular (bulk) | `AngularFlux` | `AngularSourceSink` | `AngularResidual` | `AngularDisplacement` |
| Moment (bulk) | `HarmonicMomentFlux` | **`HarmonicMomentSourceSink`** | — | `MomentDisplacement` |
| ScalarBoundary | `ScalarBoundaryFlux` | `ScalarBoundarySourceSink` | — | `ScalarBoundaryDisplacement` |
| AngularBoundary | `AngularBoundaryFlux` | `AngularBoundarySourceSink` | `AngularBoundaryResidual` | `AngularBoundaryDisplacement` |
| RadChar Interior | `…Flux` | `…SourceSink` | `…Residual` | `…Displacement` |
| RadChar Boundary | `…Flux` | `…SourceSink` | `…Residual` | `…Displacement` |
| XS coefficient | — | — | — | `CrossSectionField(CoefficientRole, ScalarField)` |

Missing cells: `HarmonicMomentResidual`, `ScalarBoundaryResidual`. Everything else is
populated (25 role leaves).

### 3.3 ⚠ VERDICT ON THE PLAN'S PROPOSED MINTS (premise partly STALE)

| plan proposes | current tree | verdict |
|---|---|---|
| `ScalarSourceSink` | **EXISTS** — `orpheus/transport/source_sinks/scalar_source_sink.py:89`, exported `source_sinks/__init__.py:60`, `_SPACE_NAME="scalar_source_sink"` | **ALREADY DONE** — nothing to mint. Renamed from `orpheus.sn.sources.IsotropicSource` at step D-F (`scalar_source_sink.py:12`) |
| `HarmonicMomentSourceSink` | **EXISTS** — `source_sinks/harmonic_moment_source_sink.py:104`, `class HarmonicMomentSourceSink(MomentField)`, `_CELL_GROUP_NAME="cell_group_source_sink"`, `UNITS=SCALAR_RATE_UNITS` (`:147`) | **ALREADY DONE** — minted by the Frame-campaign P4 as the codomain of `Λ` in `R·Λ·M` |
| `ScalarSourceSpace` | **DOES NOT EXIST.** `grep -rn "SourceSpace" orpheus/ tests/ docs/` → **zero hits.** The source-role SPACE identity today is a plain `FunctionSpace(name="scalar_source_sink", shape=(ng,*spatial))` built by `ScalarField._space_for_mesh` (`_bases.py:485-505`) | **GENUINELY NEW** — but read the caveat below |

**The caveat that decides whether `ScalarSourceSpace` earns a type.** The space types that
DO exist as subclasses are exactly those that carry structure a plain `FunctionSpace`
cannot: `ScalarTraceSpace` (`scalar_trace_space.py:59` — a `FaceLayout` + `OUTFLOW_ROW`/
`INFLOW_ROW` component constants + the face-area metric), `AngularTraceSpace` (`:307` —
+ the `omega_dot_n` selector table + the `|Ω·n̂|·w` metric), `SphericalHarmonicSpace`
(`:101` — `L` + `from_L`), `SpatialMomentSpace` (`:140` — `per_axis` + `from_per_axis`),
`FullFieldSpace` (`:152` — the block composite), `RadialCharacteristic{Interior,Boundary}Space`
(`:458`/`:473` — level/sign slot views + the SPD `V_cell` metric). A *source* space carries
**no extra structure over the flux space** — same shape, same metric, only a different
name. That is precisely the `coding-standards.md` "Type vs property — before minting a
type" test: ≥2 non-isomorphic realizations AND a non-identity morphism applied. Today
`scalar_flux` → `scalar_source_sink` differs by NAME ONLY, and no morphism converts
between them (they are non-additive by class gate, and the conversion is a named factory).
**The mechanism the tree already uses for "the source's space is not the flux's space" is
the `_SPACE_NAME` string.** Whether that string should be promoted to a type is a
judgment for the user — I report the discriminator, not the verdict.

### 3.4 Space types (`orpheus/numerics/space.py` + `spaces/`)

- `FunctionSpace(Generic[Carrier])` (`space.py:104`) — fields `name: str` (`:153`),
  `shape: tuple[int,...]` (`:154`), `inner_product_weights: NDArray|None` (`:163`).
  Metric API: `apply_metric` (`:248`), `apply_inverse_metric` (`:268`, Moore–Penrose
  pseudo-inverse so zero-weight regions zero out). Composition: `__mul__ →
  TensorProductSpace` (`:297`).
- `TensorProductSpace(FunctionSpace)` (`:381`) — `find_factor[T](factor_type)` (`:442`),
  the type-visible composition tree.
- `DualSpace(FunctionSpace)` (`:521`).
- `spaces/`: `ScalarTraceSpace`, `AngularTraceSpace`, `SphericalHarmonicSpace`,
  `SpatialMomentSpace`, `FullFieldSpace` (+ the `CompositeField`/`_CompositeLeaf`
  Protocols at `full_field_space.py:113`/`:127`), `RadialCharacteristicInteriorSpace`,
  `RadialCharacteristicBoundarySpace` (both over a private
  `_RadialCharacteristicSubSpace`, `:321`).
- Bases (`numerics/basis/`): `Basis(ABC)` + `GramStructure(Enum)` (`base.py:117`/`:84`),
  `IndicatorBasis` (`:104`), `OverlapBasis(IndicatorBasis)` (`:44`),
  `WeightedIndicatorBasis` (`:94`), `SphericalHarmonicBasis` (`:120`).

**Naming crosswalk drift (minor, verified):** the composite's interior BLOCK space is
named `"scalar_bulk"` (`diffusion/augmented_mesh.py:363`) / `"sn_bulk"` (SN side), while
the leaf that actually occupies that block is named `"scalar_flux"`. The two names are
distinct strings for the same block — a naming-consistency snag, not a correctness bug
(the composite's block spaces carry the metric; the leaf's space carries the role
identity).

---

## §4. `orpheus/numerics/frame.py` — the discipline-as-type hierarchy

`FrameBase(ABC)` (`frame.py:114`) → `PetrovGalerkinFrame` (`:334`) → `GalerkinFrame` (`:374`).
Frozen dataclasses. Fields: `basis: Basis`, `measure: DiscreteMeasure` (`:131-132`);
`PetrovGalerkinFrame` adds `test_basis: Basis` (`:366`, REQUIRED — no `None` default);
`GalerkinFrame.__init__(basis, measure)` (`:393`) binds `test_basis = basis` via
`object.__setattr__` so `test is trial`.

The single abstract hook is `test -> Basis` (`:135-144`). The **basis fixes the codomain**
(`basis_space`, `:157`), the **measure fixes the domain** (`measure_space`, `:184` — read
straight off `measure.space`, so neither space is fabricated). `coefficient_space` is never
a third parameter.

### 4.1 The two faces and their types

| face | type | domain → codomain | body |
|---|---|---|---|
| `frame.analysis` (`:195`) | `_FrameAnalysis(AnalysisOperator)` (`:412`) | `measure_space → test_space` (`:427`/`:431`) | `test.analyze(values, test_table, measure.weights)` (`:434-437`) |
| `frame.reconstruction` (`:200`) | `_FrameReconstruction(ReconstructionOperator)` (`:446`) | `basis_space → measure_space` (`:468`/`:472`) | `basis.reconstruct(coefficients, table)` (`:475-476`) |

Both faces implement `apply_transpose` (`:439-442`, `:478-479`) — the *representation*
transpose delegating to the basis. `_FrameReconstruction` explicitly overrides
`is_adjointable = True` (`:461-465`); `AnalysisOperator`'s base has `is_adjointable` (
`projection.py:123`) and the docstring notes "is_adjointable inherits the base False; the
frame FACE overrides" (`projection.py:128`).

Tabulation is cached ONCE: `table` (trial, `:148`) and `test_table` (test, `:163`);
`GalerkinFrame` overrides `test_table` to RETURN `self.table` (`:400-403`) and `test_space`
to return `basis_space` (`:405-408`) — the 0-ULP identity with the pre-hierarchy
single-discipline frame.

### 4.2 `.H` IS metric-aware — the frame does NOT do a bare transpose

**Verified, not taken on the docstring's word.** The faces expose only `apply_transpose`
(the *bare representation* transpose). `.H` comes from `LinearOperator.adjoint()` wrapping
them in `_AdjointOperator`, whose `apply` body (`numerics/operator.py:1204-1229`) is:

```python
z = inner_codomain.apply_metric(y) if inner_codomain is not None else y
result = self.inner.apply_transpose(z)
if inner_domain is not None:
    result = inner_domain.apply_inverse_metric(result)
```

i.e. `A* = G_V⁺ ∘ Aᵀ ∘ G_W` — the metric is applied by the FUNCTION SPACE
(`FunctionSpace.apply_metric`, `space.py:248`; `apply_inverse_metric`, `space.py:268`, a
**Moore–Penrose pseudo-inverse**: `np.where(nonzero, x/w, 0.0)` at `space.py:288-292`, so
zero-weight slots — the tangential `|Ω·n|=0` trace components — zero out rather than
dividing by zero). The adjoint wrapper is metric-representation-agnostic; composite
spaces OVERRIDE `apply_metric` for the per-block bulk⊕trace metric (`space.py:255-257`).
So: **frame `.H` is the G-weighted Hilbert adjoint, not `Aᵀ`** — with the caveat that when
a space carries `inner_product_weights=None` the metric application is the identity
(`space.py:262-264`), so `.H == Aᵀ` for a Euclidean space. Whether a given frame's `.H` is
genuinely weighted therefore depends on whether its `measure.space` / `basis.space` carry
weights, which is a per-frame fact, not a frame-hierarchy guarantee.

### 4.3 The composed verbs

- `conjugate(A) -> OperatorProduct(R, OperatorProduct(A, M))` (`:205-238`) — `R·A·M`; the
  spectral theorem when the frame is `A`'s eigenbasis.
- `reconstruct_after(A) -> OperatorProduct(R, A)` (`:240-250`) — for inputs already in
  coefficient space (the windowed moment iterate; wiring it to `conjugate` would
  double-project).
- `gram` (`:254-308`) — a SINGLE `analysis∘reconstruction` row-sum probe of the all-ones
  coefficient vector, installed as the coefficient space's metric. **It REFUSES a
  `GramStructure.DENSE` trial** with `NotInvertible` (`:296-305`) rather than returning a
  silently-wrong probe.
- `project(field) = gram.apply_inverse_metric(analysis.apply(field))` (`:310-330`) — the
  homogenise / condense verb.

---

## §5. DIFFUSION solver — how it builds operators, and how the eigenproblem is posed

`orpheus/diffusion/` = 6 modules: `__init__.py`, `augmented_mesh.py` (368),
`boundary_realizer.py`, `method_space.py`, `operators.py` (669), `solver.py` (419).

### 5.1 The composition — the exact lines

`orpheus/diffusion/solver.py:228-243` (`DiffusionSolver.__init__`):

```python
space = mesh.full_field_space
self.leakage = LeakageOperator(mesh)
collision = MultiplicationOperator(
    self.mat_xs.total_cross_section_field, space=space,
)
scattering = (
    IsotropicScattering(self.mat_xs, space=space)
    + IsotropicN2N(self.mat_xs, space=space)
)
self.boundary = DiffusionBoundaryOperator(mesh)
self.loss = self.leakage + collision - scattering - self.boundary
self.fission = FissionOperator.from_solver_data(
    mat_xs=self.mat_xs, full_field_space=space,
)
```

**The loss operator is composed honestly as `L + C − S − B`** (`:240`) via the operator
algebra's `__add__`/`__sub__` — no fused matvec, no hand-inlined stencil.

### 5.2 Which operators are SHARED vs diffusion-owned

| leaf | provenance | verdict |
|---|---|---|
| `C` collision | `transport.operators.multiplication_operator.MultiplicationOperator` | **SHARED** with SN + homogeneous |
| `S` scattering | `transport.operators.isotropic_scattering.{IsotropicScattering, IsotropicN2N}` | **SHARED with homogeneous ONLY** — SN uses `transport/operators/scattering.py`'s `ScatteringOperator` family instead |
| `F` production | `transport.operators.fission.FissionOperator.from_solver_data` | **SHARED** with SN + homogeneous |
| `L` leakage | `diffusion.operators.LeakageOperator` (`operators.py:272`), `block_role = BlockRole.FULL` (`:312`) | diffusion-OWNED |
| `B` boundary | `diffusion.operators.DiffusionBoundaryOperator` (`:544`), `block_role = BlockRole.BOUNDARY` (`:576`) | diffusion-OWNED |

So diffusion **does** go through `MultiplicationOperator` — but through a DIFFERENT
constructor entry than homogeneous does: diffusion calls
`MultiplicationOperator(field, space=space)` (the plain constructor, `solver.py:230-232`)
while homogeneous calls `MultiplicationOperator.from_mesh(field, mesh)`
(`homogeneous/solver.py:143-145`, factory at `multiplication_operator.py:299`).

`MultiplicationOperator` dispatches on the input type via `@singledispatchmethod`
(`multiplication_operator.py:376`, arms registered at `:384` and `:445`).
`IsotropicScattering`/`IsotropicN2N` instead use an inline
`if isinstance(phi, FullField): return _scalar_composite_source(self, phi)`
(`isotropic_scattering.py:274-278`, mirrored at `:374+`). **Two different spellings of the
same "bare ndarray OR composite" dispatch problem inside the same operators package** —
see §7.

### 5.3 The eigenvalue posing — a HELD `(A, F)` pair, plus an eagerly-LU'd `A⁻¹`

- `self.loss` and `self.fission` are **held separately as typed, un-materialized
  operators** (`solver.py:240`, `:241-243`). There is **no** eagerly-fused `A⁻¹F`.
- `A⁻¹` **IS** eagerly materialized: `self.resolvent = MatrixInverseOperator(
  FlattenedOperator(self.loss, self.template))` (`solver.py:250-252`) — one dense LU at
  construction, one back-substitution per outer iteration. The module docstring records
  the campaign ruling: **NEVER route through the structure-keyed `A.inverse()`** because
  the Green splitting diverges for fine-mesh elliptic operators (`solver.py:48-50`,
  restated at `operators.py:117-119`).
- The drive is `power_iteration(solver, max_iter=max_outer)` (`solver.py:407-409`) through
  the `EigenvalueSolver` Protocol. The 6 protocol methods:
  `initial_flux_distribution` (`:273`), `compute_fission_source` (`:279`),
  `solve_fixed_source` (`:287` — a SINGLE `A⁻¹q`; **no scattering inner iteration exists
  at all**, because `A` carries the full multigroup coupling), `compute_production_rate`
  (`:296`), `compute_keff` (`:306`), `converged` (`:316`).
- **The protocol carrier is the flat composite vector** (`FullField.to_flat` — bulk C-ravel
  ⊕ trace buffer). Conversion happens at exactly two sites: `DiffusionSolver.unflatten`
  (`:257-261`) and the template-frozen `FlattenedOperator` (`:37`).
- `compute_keff` is the **integrated eigenvalue relation through the same operator that
  defines the fixed point**: `production / _volume_integral(self.loss.apply(psi).interior.values)`
  (`:311-314`) = production / (absorption + leakage) by the column-sum theorem + divergence
  telescoping. `B` contributes nothing to the bulk block; the trace constraint rows are
  excluded by construction.
- **(n,2n) is LOSS-side** here (`S = Σ_s0ᵀ + 2Σ₂ᵀ`, `F = νΣ_f` only) — mirroring
  homogeneous; SN poses it production-side (`solver.py:72-82`).

### 5.4 `MaterialXSField.diffusion_coefficient` IS the D-source — confirmed

`LeakageOperator.__init__` reads
`D = np.asarray(mesh.material_xs_field().diffusion_coefficient, float)`
(`diffusion/operators.py:316`). That property is `material_xs_field.py:641-667` — a lazy,
cached per-cell gather of `Mixture.diffusion_coefficient` (the #290 P1 seam:
`Σ_tr = Σ_t − Σ_{g'}Σ_{s,1}(g→g')`, exactly `Σ_t` in the P0-only limit), reshaped to the
principled `(ng, *spatial)` layout. Single read path — `grep` finds no other
`diffusion_coefficient` consumer in `orpheus/`.

The two discretization kernels are module-level pure functions so the formulas are stated
once and are mutation-testable: `_interior_conductance` (`operators.py:170-192`) and
`_boundary_closure` (`:195-218`). `LeakageOperator` ALSO implements the assembly mode
(`is_assemblable = True`, `:438`; `assemble()` → `SparseAssembledOperator`, `:441-541`)
with an explicit one-source discipline note (`:446-454`): every emitted value reads the
SAME precomputed attributes `apply` consumes. `DiffusionBoundaryOperator.assemble`
(`:626-668`) extracts each face block THROUGH `law.as_matrix(basis_shape=(ng,))` rather
than re-reading the albedo scalar.

**Diffusion's phase-space admission is the MESH's contract, not the solver's** — the 1-D
refusal and the bounded-geometry refusal fire in `DiffusionMesh._init_core`
(`augmented_mesh.py:204-218`), so a `LeakageOperator` needs no operator-side refusal to
drift (`operators.py:299-302`).

---

## §6. HOMOGENEOUS solver — the contrast case (EAGER dense `A⁻¹F`)

`orpheus/homogeneous/solver.py` (248 lines, the whole package).

- **Loss:** `_assemble_loss_operator(mat_xs)` (`:115-147`) returns an UN-materialized
  `OperatorSum`: `collision - k_iso` where
  `collision = MultiplicationOperator.from_mesh(mat_xs.total_cross_section_field, mat_xs.mesh)`
  (`:143-145`) and `k_iso = IsotropicScattering(mat_xs) + IsotropicN2N(mat_xs)` (`:146`).
  Streaming `L` is identically zero in an infinite medium and is **dropped**, not built
  as a zero operator.
- **Carrier:** a MESHLESS single-cell `MaterialMesh.from_materials({0: mix})` (`:180`) —
  the same method-agnostic data carrier the meshed solvers use.
- **Posing — this is the key contrast with diffusion:** the multiplication operator is
  **eagerly composed AND densely materialized**:
  ```python
  K = MatrixInverseOperator(loss, basis_shape=(ng, 1)) @ production   # :194
  k_inf, phi = dominant_eigenpair(K.as_matrix(basis_shape=(ng, 1)))   # :202
  ```
  So homogeneous **DOES** hold a fused `K = A⁻¹F` and materializes it dense, then extracts
  the dominant eigenpair via the shared Perron–Frobenius primitive
  `numerics.eigenvalue.dominant_eigenpair`. Diffusion holds `(A, F)` apart and iterates.
  Both are defensible (0-D is exactly solvable), but they are **two different posings of
  the same `A⁻¹F` eigenproblem through the same shared leaves**.
- **`basis_shape=(ng, 1)` is passed EXPLICITLY** at both `:194` and `:202` because the
  meshless operators carry no `domain` space to derive it from (`:136-141`). This is the
  "space-anonymous / model-portable bare-ndarray contract" (`isotropic_scattering.py:245-250`).
- **Reduced vs dense energy transfer:** homogeneous consumes the **reduced** operators
  (`IsotropicScattering`/`IsotropicN2N`, which route the per-material `Σ_s0` matmul through
  `MaterialXSField.apply_p0_in_scatter`), then materializes the COMPOSED `K` once. It does
  NOT hand-build a dense `Σ_s0ᵀ` matrix.
- **Reaction rates:** `IntegratedReactionRate(mat_xs.fission_production_field)` /
  `(…absorption_cross_section_field)` (`:210-211`) — the SAME shared functional diffusion
  uses at `diffusion/solver.py:302-304`.
- **(n,2n) loss-side** (`:27-30`), matching diffusion.
- ⚠ **Variable shadowing** (`:193` vs `:210`): the name `production` is bound first to the
  `FissionOperator` and then REBOUND to an `IntegratedReactionRate` — two different objects
  under one name in a 55-line function. Not a bug (the operator is consumed at `:194`
  before the rebind) but it defeats the named-intermediate discipline.

---

## §7. Shared-vs-duplicated audit

### 7.1 Nexus `twin_paths` (min_similarity 0.65, `min_tokens` default, tests/tools/derivations excluded)

60 pairs returned repo-wide; **exactly 2** have BOTH ends inside
{transport, sn, diffusion, homogeneous, numerics}:

| sim | pair | my judgment |
|---|---|---|
| **1.000** | `DiffusionSolver.converged` (`diffusion/solver.py:316`) ↔ `SNSolver.converged` (`sn/solver.py:1441`) | **REAL TWIN.** I read both bodies: byte-equivalent modulo a line-continuation style (`if iteration <= 2: return False` → `dk = abs(...)` → `dphi = norm(Δ)/max(norm,1e-30)` → `bool(dk<keff_tol and dphi<flux_tol)`). See 7.2 — it's actually a QUINTUPLET, not a pair. |
| 1.000 | `RadialCharacteristicBoundaryField.__post_init__` (`_bases.py:1231`) ↔ `RadialCharacteristicInteriorField.__post_init__` (`_bases.py:1171`) | **SYMMETRIC-BY-DESIGN.** Two sibling loci (`FaceField` siblings, explicitly NOT parent/child — `_bases.py:1158-1160`), each narrowing to its OWN space type with its own message. Collapsing them would require a shared `_REQUIRED_SPACE` ClassVar; the gain is ~6 lines and the cost is one more indirection on the guard. Not worth flagging as debt. |

The scan's blind spot is worth naming: `twin_paths` uses AST-shingle similarity, so it
finds `converged` (identical statements) but MISSES the semantic duplication in §7.2/§7.3
where the *reshape convention* differs while the *math* is the same.

### 7.2 The `EigenvalueSolver`-protocol implementation is spelled FIVE times

`numerics/eigenvalue.py:104-186` declares the `EigenvalueSolver` /
`ProductionRateSolver` Protocol. Implementers found by grep:

| implementer | file:line | `converged` | tolerance defaults |
|---|---|---|---|
| `KEigenvalue` (generic `(A,S,F)`) | `numerics/iteration.py:1049` | `:1326` — **carrier-honest** (`_l2_norm` over the ravellable protocol, iterate diff via the carrier's `Vector.__sub__`) | 1e-7 / 1e-6 (`:1145`) |
| `DiffusionSolver` | `diffusion/solver.py:174` | `:316` | 1e-10 / 1e-9 (`:216`) |
| `SNSolver` | `sn/solver.py` | `:1441` | 1e-7 / 1e-6 (`:922`) |
| `CPSolver` | `cp/solver.py:~460` | `:722` | 1e-6 / 1e-5 (`:59`) |
| MOC | `moc/core.py:~35` | `:242` | 1e-6 / 1e-5 (`moc/solver.py:76`) |

**The generic body already EXISTS and is deliberately bypassed.** `KEigenvalue`'s docstring
states it outright (`iteration.py:1121-1125`): *"the method-layer solvers implement the
`EigenvalueSolver` protocol directly by design and never routed through this class."* Its
only production consumer is the SN **adjoint** path (`sn/solver.py:2472-2478`, inside
`solve_sn_adjoint`, on the daggered triple `resolvent.H, gain.H, F_posed.H`). Diffusion,
CP, and MOC bypass it entirely.

My judgment: **the `converged` quintuplet is REAL duplication of a method-neutral concept**
(a convergence test on `(Δk, Δφ/‖φ‖, iteration)` has no method content at all), and the
`KEigenvalue.converged` version is strictly the BEST one (carrier-generic, handles the
`norm==0` case which the four copies do not — they'd return `Δ/1e-30`). The tolerance
DEFAULTS legitimately differ per method (a diffusion direct solve can afford 1e-10; an
inexact MOC cannot), so the differing numbers are NOT the smell — the re-spelled body is.
This is a candidate for a shared mixin/free function; I am NOT proposing the refactor,
only reporting that the shared body pre-exists and is unused by 4 of 5 implementers.

### 7.3 The volume-integral reduction — three spellings, two conventions

The pattern "rate density `(k, *spatial)` → flat `(N_cells, k)` → `mesh.volume_measure` →
`.sum()`":

| site | reshape |
|---|---|
| `diffusion/solver.py:263-269` (`_volume_integral`) | `np.moveaxis(rate, 0, -1).reshape(-1, ng)` — **rank-generic** |
| `transport/reaction_rate_functional.py:224-228` | `np.moveaxis(per_cell, 0, -1).reshape(-1, 1)` — **rank-generic**, k=1 |
| `sn/solver.py:1173-1185` (`compute_production_rate`) | `transpose(1, 2, 0)` — **hardcoded to 2 spatial axes** |

The first two are the same body at different `k`; the diffusion one even says so in its own
docstring ("the same reduction `IntegratedReactionRate` rides", `solver.py:265-266`) while
declining to call it. The SN one uses a rank-FIXED transpose while the shared helper is
rank-generic — a latent 3-D hazard, though `sn/solver.py` is explorer (b)'s surface so I
flag rather than adjudicate.

### 7.4 Two dispatch spellings for "bare ndarray OR `FullField`" inside one package

`orpheus/transport/operators/`:
- `MultiplicationOperator` uses `@singledispatchmethod` on `_apply_impl` with registered
  arms (`multiplication_operator.py:376`, `:384`, `:445`) and a `TYPE_CHECKING`-only
  overload triple (`:465-468`);
- `IsotropicScattering` / `IsotropicN2N` use an inline `if isinstance(phi, FullField):
  return _scalar_composite_source(self, phi)` then fall through to the bare arm
  (`isotropic_scattering.py:274-278`; `:374+`), typed `phi: "np.ndarray | FullField | object"`.

Both are "the kernel is the SAME bare arm either way; the composite arm wraps it." Same
problem, two mechanisms, same directory. My judgment: **a real inconsistency**, mild in
consequence (both are correct) but exactly the "which path is canonical?" smell. Also
worth noting: `IsotropicScattering.apply_transpose` REFUSES a composite loudly
(`isotropic_scattering.py:285-291`) pending "the adjoint diffusion consumer, #281" — so
the composite transpose arm is a known, named gap, not an oversight.

### 7.5 Nexus `discriminations` (min_sites 2) — in-scope tags

| tag | sites | in-scope discriminators | judgment |
|---|---|---|---|
| `coord` | 17 | `transport/mesh/axis.py:493` (`axes_from_legacy_mesh`), `:587` (`legacy_mesh_from_axes`), `:315` (`RadialAxisMesh.__post_init__`); `sn/mesh/augmented_mesh.py:215`; `sn/sweep/pole_angular_closure.py:1346` | The big one. But the bulk of the 17 sites are in `geometry/` (`coord.py`) and tests. Inside MY surface the 3 `axis.py` sites are the legacy-mesh ↔ axes CONVERTER pair — a genuine boundary adapter, which is where a tag discrimination legitimately belongs (discriminate once, at the boundary). Not a missing type. |
| `geometry` | 8 | `geometry/mesh.py:399` only (rest are tests) | Boundary parse — legitimate. |
| `kind` | 6 | `sn/acceleration/dsa.py:257`, `sn/operators/boundary.py:573`, `geometry/mesh.py:71` | Mixed tag namespace (`direct`/`eigen`/`fixed_source`/`iterative` AND `reflective`/`vacuum` AND `closure_sign`/`denom_sign` under ONE tag name `kind`) — this is tag-name COLLISION across unrelated domains, not one over-discriminated tag. |
| `inner_solver` | 3 | `sn/solver.py:2915`, `:917`, `:1134` | SN-internal, explorer (b)'s. |
| `inner_schedule` | 2 | `sn/solver.py:707`, `:917` | SN-internal. |
| `curvature` | 2 | `sn/loss_representation/__init__.py:3026`, `:3377` | apply/transpose pair — symmetric-by-design. |
| `axis` | 3 | `numerics/quadrature/directional.py:150`, `numerics/symmetry.py:750`, `sn/boundary/angular.py:127` | `"x"`/`"y"`/`"z"` → index. Three sites resolving the same axis-name→index map. **Mildly real**: `AXIS_NAMES` already exists in `transport/mesh/axis.py` and is the canonical inventory, but numerics (L1) cannot import transport (L2), so `numerics/quadrature/directional.py` and `numerics/symmetry.py` must carry their own. The layering makes this duplication structural, not sloppy. |

**Nothing in `transport/`, `diffusion/`, or `homogeneous/` discriminates on a method-name
or method-kind tag.** The #290 P7b registry dissolution held: there is no
`if method == "sn"` anywhere. `TransportMethod`'s own refusal path uses
`type(method).__name__` for the MESSAGE only (`method.py:296`), never for dispatch.

### 7.6 Nexus `dead_functions` restricted to {transport, numerics, diffusion, homogeneous, sn}

Zero hits in `diffusion/` or `homogeneous/`. In-scope hits and their true status:

- **False positives via abstract-hook dispatch** (the `callers()==0`-but-live pattern,
  lesson L-001): every `_phase_space_shape` (`_bases.py:385`, `:507`, `:616`;
  `moment_displacement.py:51`) and every `_check_partner` override (`_bases.py:198`,
  `:639`, `:835`) — these are called through `BulkField.__post_init__` / `Field._check_partner`
  `super()` chains, which the static graph does not attribute. **Live.**
- **`AngularField._integrate_angular_values` (`_bases.py:411`)** — reached only through the
  role leaves' wrappers (`AngularFlux.integrate_angular`, `AngularDisplacement.integrate_angular`);
  live via those. **Live.**
- **`Composite._from_flat` (`full_field.py:454`)** / **`TimedFullField._recombine`
  (`timed_full_field.py:258`)** / **`CompositeField._recombine`
  (`numerics/spaces/full_field_space.py:143`)** — template-method hooks called from the
  base's `from_flat`. **Live.**
- **`ScatteringOperator._aniso_source_from_moment_values` (`scattering.py:556`)** — degree
  44, reached through the singledispatch arms. Explorer (a)'s surface; **not dead**.
- `numerics/iteration.py:339` `_CarrierMatvecOperator._matvec` — scipy `LinearOperator`
  hook, invoked BY scipy. **Live.**
- `numerics/quadrature/registry.py` `_*_invert` / `_*_node_count` (8 functions, `:303`–`:425`)
  — a **registry of function pointers**; live via table lookup, invisible to the call graph.
- `numerics/quadrature/directional.py:103` `_octant_sign_predicate` and
  `numerics/symmetry.py:732`/`:741` `_rotation_x`/`_rotation_y` — these are the only
  plausible genuine-dead candidates in my surface; both are private, undecorated, and
  sit next to table-driven siblings. Worth a grep-level check by whoever owns numerics,
  but note the registry pattern above makes even these suspect-not-proven.

**Conclusion: no dead code in `transport/`, `diffusion/`, `homogeneous/`.** The
`dead_functions` output for this surface is dominated by the two known graph-blinding
patterns (abstract-hook `super()` dispatch and function-pointer registries).

---

## §8. Layering — ZERO violations in the intended stack

The intended stack is **enforced by a test**, not just documented:
`tests/test_layer_imports.py` (L0 `derivations` < L1 `numerics` < L2 `transport` <
L3 `{sn, pn, moc, cp, mc, diffusion, kinetics, fuel, thermal_hydraulics, homogeneous}`),
with the forbidden-direction table at `:55-69` and an **explicit TYPE_CHECKING tolerance**
at `:147-148` (`if is_tc and src_pkg in (L1_PACKAGES | L2_PACKAGES): continue`).

Measured against the tree:

| edge | result |
|---|---|
| `numerics` → `transport` or any L3 | **NONE** (grep for `from orpheus.{transport,sn,diffusion,homogeneous,cp,moc,mc,fuel,kinetics,thermal_hydraulics}` in `orpheus/numerics/` returns zero) |
| `transport` → any L3, at RUNTIME | **NONE** |
| `transport` → any L3, under `TYPE_CHECKING` | **9 sites** (the §1.0 table) — all covered by the test's explicit tolerance |
| `diffusion` → `sn` / another L3 | **NONE** — `diffusion/` imports only `geometry`, `numerics`, `transport`, and its own siblings (verified across all 6 modules) |
| `homogeneous` → another L3 | **NONE** — `homogeneous/solver.py:42-52` imports `data`, `numerics`, `transport` only |
| the one sanctioned exception | `("derivations/continuous/mms/sn.py", "transport")` (`test_layer_imports.py:88`) — an L0 MMS reference importing L2 transport vocabulary for its prescribed-inflow source |

So the layering INVARIANT the plan cares about — `transport` must not import method
packages — **holds today, at runtime and by test**.

**The residual coupling is in the TYPE ANNOTATIONS, not the imports.** 7 of the 9
`TYPE_CHECKING` sites narrow a member or parameter to `SNMesh` specifically:
`AngularField.mesh: "SNMesh"` (`_bases.py:344`), `MomentField.mesh: "SNMesh"` (`:604`),
`AngularBoundaryField.mesh: "SNMesh"` (`:1030`), `RadialCharacteristic*Field.mesh: "SNMesh"`
(`:1168`, `:1228`), plus the `from_mesh`/`zeros_on` signatures on those families. This is
DELIBERATE and documented (`_bases.py:161-169`: the covariant narrowing keeps `mesh.quad`
honest with no casts, #267). But it means the honest statement is:

> `orpheus/transport/` is method-neutral **as a runtime dependency graph**. As a TYPE
> surface it is bifurcated: the Scalar family is genuinely generic (`mesh: MaterialMesh`,
> `_bases.py:477`/`:487`/`:517` — which is exactly why diffusion can use it), while the
> Angular / Moment / RadialCharacteristic families are typed to `SNMesh` by construction.
> The only place that names BOTH witnesses is `_bases.py:116-118`, and the only
> method-generic type in the package is `TransportMethod` itself.

### 8.1 What is genuinely method-neutral today (the summary answer)

**Genuinely neutral and PROVEN so by two consumers:**
`MaterialMesh` + `MaterialXSField` + `axis.py` (the face inventory), `FullField` /
`Composite` (generic in `[Interior, Boundary]`, `full_field.py:178`/`:495`),
`ScalarFlux` / `ScalarSourceSink` / `ScalarBoundaryFlux` / `ScalarBoundarySourceSink`,
`IntegratedReactionRate`, `MultiplicationOperator`, `FissionOperator`,
`IsotropicScattering` + `IsotropicN2N`, `TransportMethod` + `resolve_boundary_conditions`,
the whole of `numerics/` (space/spaces/basis/frame/measure/moment_layout).

**Neutral-by-DESIGN but single-consumer (SN only) today:**
`AngularFlux` family, `HarmonicMomentFlux` family, `AngularTraceSpace`,
`RadialCharacteristic*` (all 4 role families × 2 loci), `HarmonicFrame`,
`numerics/quadrature/`, `numerics/symmetry.py`.

**At L2 but SN-vocabulary throughout (the weakest neutrality claim):**
`transport/spatial/` — 4.7 kloc (`scheme.py`, `diamond.py`, `linear_discontinuous.py`,
`cell_balance.py`, `_ubld.py`). Its docstrings reference `SNMesh.dag_walk`,
`SweepDependencyGraph.walk_full/walk_windowed`, `OctantLabel.signs`, `CollisionCache`,
`pole_angular_closure`, `pair_diffusion_limit_consistent`. Diffusion consumes NONE of it
(it writes its own two-line FD kernels at `diffusion/operators.py:170-218`). The promotion
out of `orpheus.sn.sweep` at #272 was deliberate (`transport/spatial/__init__.py:12`), but
nothing has yet arrived to make the L2 residence load-bearing.

**Neutral machinery that EXISTS but is bypassed:** `KEigenvalue` (§7.2) — the generic
`(A, S, F)` k-eigenvalue posing layer with the best `converged` body in the codebase,
consumed by exactly one production path (the SN adjoint, `sn/solver.py:2474`).

### 8.2 Sharpening §8: the type-narrowing is SYMMETRIC, not SN-only

Full census of `mesh:` annotations in `_bases.py` (verified by grep):

| ABC | `mesh` annotation | line |
|---|---|---|
| `BulkField` | `MaterialMesh` | `:172` |
| `ScalarField` | (inherits `MaterialMesh`) — factories take `MaterialMesh` | `:477`, `:487`, `:517`, `:530`, `:546` |
| `AngularField` | **`SNMesh`** | `:344` (factories `:352`–`:452`) |
| `MomentField` | **`SNMesh`** | `:604` (factories `:660`–`:719`) |
| `FaceField` | `MaterialMesh` | `:785` (factories `:901`, `:917`, `:969`) |
| `AngularBoundaryField` | **`SNMesh`** | `:1030` |
| `ScalarBoundaryField` | **`DiffusionMesh`** | `:1093` |
| `RadialCharacteristic{Interior,Boundary}Field` | **`SNMesh`** | `:1168`, `:1228` |

So the narrowing pattern is **one method-mesh per family**, applied symmetrically — the
scalar BOUNDARY family is narrowed to `DiffusionMesh` exactly as the angular boundary
family is narrowed to `SNMesh`. The ONLY families left at `MaterialMesh` are the two
generic bases (`BulkField`, `FaceField`) and `ScalarField` (the bulk scalar). This is why
`diffusion/operators.py` can build `ScalarSourceSink.from_mesh(out_bulk, self.mesh)`
(`:431`) but the trace side must go through `ScalarBoundarySourceSink.zeros_on(self.mesh)`
(`:420`, `:610`) with `self.mesh: DiffusionMesh`.

A consequence worth naming for a DSA-style consumer: `ScalarBoundaryField` typed to
`DiffusionMesh` means an SN-side low-order system that wants a scalar trace must PROMOTE
(`DiffusionMesh.from_material_mesh(sn_mesh)` — `augmented_mesh.py:253-280`; an `SNMesh`
IS a `MaterialMesh`). The `_bases.py:1080-1084` docstring anticipates exactly this. The
R4 DSA ruling chose not to (§2), so the promotion path is BUILT but currently unexercised
in production.

### 8.3 Verified `orpheus/spatial` consumer census (the L2-residence question)

`grep` for importers of `orpheus.transport.spatial` across `orpheus/` — **9 files, ALL in
`orpheus/sn/`**: `sn/loss_representation/{__init__,assembly,sweep_graph}.py`,
`sn/mesh/augmented_mesh.py`, `sn/solver.py`, `sn/sweep/{cache,pairing,psi_half_angle_seed,scan}.py`.
Zero importers in `diffusion/`, `homogeneous/`, `cp/`, `moc/`, `mc/`, `numerics/`.

### 8.4 Per-package import sets (the layering ledger)

- `orpheus/diffusion/*` (6 modules): `orpheus.transport` ×18, `orpheus.numerics` ×15,
  `orpheus.geometry` ×6, `orpheus.diffusion` ×6 (siblings), `orpheus.data` ×2.
  **No L3 sibling package.**
- `orpheus/homogeneous/solver.py`: `orpheus.data`, `orpheus.numerics`, `orpheus.transport`.
  **No L3 sibling package.**
- `orpheus/transport/reaction_rate_functional.py`: `orpheus.numerics.functional`,
  `orpheus.transport.fields.cross_section_field` only — genuinely L2-clean.
