# Grand Report: Native Mathematical Architecture for a Neutron Transport Code

**Purpose.** This report is intended to guide development in Codex or any other coding agent. It is written as a design document where names, class boundaries, operator algebra, and mathematical contracts are deliberately chosen so that the code itself carries high semantic signal.

**Revision note.** This version explicitly corrects vacuum-boundary semantics: vacuum means zero inflow trace / no incoming current on `Gamma_-`, not zero scalar flux on the full boundary. It also adds a symmetry-sector/orbifold design layer for octant and quotient geometries using group actions, tensor-network boundary gluing, quotient measures, and equivariant chain/cochain consistency tests. The present revision additionally integrates the post-sheaf design discussion: `MoC` receives a targeted `CharacteristicTransportSheaf` / `TrackMonodromy` layer; the other methods receive their own native local-to-global structures instead of being forced into sheaf language; and boundary conditions are decomposed into `BoundaryTraceLaw = geometry map + physical response + source`, followed by method-specific `BoundaryRealizer` objects.

The intended solver family includes:

- discrete ordinates, `SN`
- spherical harmonics, `PN`
- method of characteristics, `MoC`
- collision probabilities, `CP`
- Monte Carlo, `MC`
- stochastic PDE / stochastic transport formulations
- tensor-network and tensor-train representations
- hybrid deterministic / stochastic / Monte Carlo methods
- `k`-eigenvalue and `alpha`-eigenvalue problems

The central principle is:

> **Do not build one code per method. Build one mathematical language of spaces, measures, fields, kernels, flows, projections, operators, and solvers. Each transport method is then a representation of the same transport objects.**

---

## 1. Notation and operator convention

Use the following operator names consistently:

| Symbol | Code name | Meaning |
|---|---|---|
| `T` | `TimeMass`, `TimeOperator` | coefficient of the time derivative, usually multiplication by `1 / v` |
| `L` | `StreamingOperator` | streaming / transport operator, `Omega · grad` |
| `C` | `CollisionOperator` | total interaction / removal, usually multiplication by `Sigma_t` |
| `S` | `ScatteringOperator` | scattering source operator |
| `F` | `FissionOperator` | fission production operator |
| `B` | `BoundaryOperator` | maps outgoing traces to incoming traces, or applies boundary constraints |
| `A_loss` | `LossOperator` | `L + C - S` |
| `A_prompt` | `PromptKineticsOperator` | `L + C - S - F` |
| `G` | `KineticGenerator` | time-evolution generator, typically `-T^{-1} A_prompt` |

The time-dependent transport equation should be represented as

```text
T dpsi/dt + (L + C - S - F) psi = q.
```

That is,

```python
A_prompt = L + C - S - F
problem = InitialValueProblem(T=T, A=A_prompt, source=q, initial=psi0)
```

The fixed-source equation is

```text
(L + C - S - F) psi = q
```

or, if fission is treated as part of the source iteration / eigenproblem,

```text
(L + C - S) psi = F psi + q.
```

For criticality, use

```text
(L + C - S) psi = (1/k) F psi.
```

Define

```python
A_loss = L + C - S
K = A_loss.inverse() @ F
```

Then the `k`-eigenvalue problem is

```text
K psi = k psi.
```

For alpha modes, define

```text
T dpsi/dt + A_prompt psi = 0.
```

With ansatz

```text
psi(t) = exp(alpha t) phi,
```

the alpha problem is

```text
A_prompt phi = -alpha T phi.
```

Equivalently,

```text
G phi = alpha phi,
where G = -T^{-1} A_prompt.
```

This sign convention makes positive `alpha` mean exponential growth and negative `alpha` mean decay.

---

## 2. The core ontology

The architecture should be built around the following primitive mathematical objects.

```text
MathObject
├── Domain
├── Measure
├── Space
├── Basis
├── Field
├── Operator
├── Kernel
├── Projection
├── Flow
├── Process
├── Functional
├── Problem
├── Solver
└── Representation
```

Each object should carry its own invariants, and those invariants should be exposed as named tests.

Example:

```python
lebedev.assert_spherical_harmonic_exactness(L=17)
S.assert_positive()
S.assert_conservative()
B.assert_maps_inflow_to_outflow()
A.assert_adjoint_consistency(samples=10)
P.assert_projection_idempotent()
```

The point is not just verification. The point is that the names make the code legible to an LLM and to humans.

---

## 3. Protocols, dataclasses, and abstract base classes

Use three different Python mechanisms for three different purposes.

### 3.1 Use `Protocol` for capability contracts

A `Protocol` says: “anything with these methods can participate.” Use it for structural typing and composability.

Good protocol examples:

```python
from typing import Protocol, TypeVar, runtime_checkable

X_co = TypeVar("X_co", covariant=True)
Y_co = TypeVar("Y_co", covariant=True)

@runtime_checkable
class SupportsApply(Protocol[X_co, Y_co]):
    def apply(self, x: X_co) -> Y_co:
        ...

@runtime_checkable
class SupportsAdjoint(Protocol):
    @property
    def H(self):
        ...

@runtime_checkable
class SupportsTensorProduct(Protocol):
    def tensor(self, other):
        ...

@runtime_checkable
class SupportsIntegration(Protocol):
    def integrate(self, f):
        ...

@runtime_checkable
class SupportsSampling(Protocol):
    def sample(self, n: int, rng=None):
        ...
```

Use protocols for contracts such as:

```python
class Positive(Protocol):
    def positivity_defect(self) -> float:
        ...

class Conservative(Protocol):
    def conservation_defect(self) -> float:
        ...

class Markovian(Protocol):
    def normalization_defect(self) -> float:
        ...

class InvariantUnderGroup(Protocol):
    def invariance_defect(self, group_action) -> float:
        ...
```

Names like `Positive`, `Conservative`, `Markovian`, and `InvariantUnderGroup` are high-signal because they imply immediate tests:

```python
assert obj.positivity_defect() <= tol
assert obj.conservation_defect() <= tol
assert obj.normalization_defect() <= tol
assert obj.invariance_defect(group) <= tol
```

### 3.2 Use `@dataclass` for immutable value objects

Use dataclasses for objects that are mostly data and should be comparable, serializable, and easy to inspect.

Prefer:

```python
from dataclasses import dataclass, field
from typing import Mapping, Any

@dataclass(frozen=True, slots=True)
class BoundaryDeclaration:
    patch: str
    kind: str
    parameters: Mapping[str, Any] = field(default_factory=dict)

@dataclass(frozen=True, slots=True)
class MaterialAssignment:
    region: str
    material: str

@dataclass(frozen=True, slots=True)
class GeometrySpec:
    dimension: int
    regions: tuple
    material_assignments: tuple[MaterialAssignment, ...]
    boundary_declarations: tuple[BoundaryDeclaration, ...]
```

Use `frozen=True` for specifications. A `GeometrySpec` should not mutate after creation. If the user wants a modified geometry, create a new one.

Good dataclass candidates:

```text
GeometrySpec
RegionSpec
BoundaryDeclaration
MaterialAssignment
MeshSpec
EnergyGroupStructure
QuadratureSpec
CubatureSpec
ExactnessCertificate
SymmetryCertificate
SolverOptions
EigenSolverOptions
TrackingOptions
TalliesSpec
```

### 3.3 Use abstract base classes for semantic base classes with invariants

Use `ABC` when the class has shared behavior, invariants, or factory registration.

Good ABC candidates:

```python
from abc import ABC, abstractmethod

class Space(ABC):
    @abstractmethod
    def dual(self):
        ...

    @abstractmethod
    def zero(self):
        ...

class Measure(ABC):
    @abstractmethod
    def integrate(self, f):
        ...

class LinearOperator(ABC):
    domain: Space
    codomain: Space

    @abstractmethod
    def apply(self, field):
        ...

    @property
    @abstractmethod
    def H(self):
        ...

class BoundaryOperator(LinearOperator):
    @abstractmethod
    def resolve(self, trace_space):
        ...
```

Use ABCs when subclassing communicates mathematical kind:

```text
LinearOperator
├── MultiplicationOperator
├── DifferentialOperator
├── IntegralKernelOperator
├── ProjectionOperator
├── TraceOperator
├── BoundaryOperator
├── GreenOperator
├── SweepOperator
├── TensorProductOperator
├── SumOfTensorProductsOperator
└── TensorNetworkOperator
```

---

## 4. Self-registering derived classes

Self-registration is appropriate for extensible families such as:

- boundary conditions
- quadrature / cubature rules
- cross-section formats
- discretization methods
- method-space builders
- solvers
- eigenvalue solvers
- Monte Carlo estimators
- tensor-network compression backends

Use a small registry mixin.

```python
from __future__ import annotations

from abc import ABC
from typing import ClassVar, TypeVar, Generic

T_cls = TypeVar("T_cls")

class RegistryMixin:
    registry: ClassVar[dict[str, type]] = {}
    key: ClassVar[str | None] = None

    def __init_subclass__(cls, key: str | None = None, **kwargs):
        super().__init_subclass__(**kwargs)
        if key is not None:
            cls.key = key
            base = cls._registry_base()
            if key in base.registry:
                raise KeyError(f"Duplicate registry key {key!r} for {base.__name__}")
            base.registry[key] = cls

    @classmethod
    def _registry_base(cls):
        # Override in each independent registry root.
        return cls

    @classmethod
    def create(cls, key: str, **kwargs):
        try:
            concrete = cls.registry[key]
        except KeyError as exc:
            available = ", ".join(sorted(cls.registry))
            raise KeyError(f"Unknown {cls.__name__} key {key!r}. Available: {available}") from exc
        return concrete(**kwargs)
```

Use separate registry roots so `BoundaryOperator.registry` does not mix with `Solver.registry`.

```python
class BoundaryOperator(LinearOperator, RegistryMixin, ABC):
    registry: ClassVar[dict[str, type["BoundaryOperator"]]] = {}

    @classmethod
    def _registry_base(cls):
        return BoundaryOperator
```

Then derived boundary classes self-register:

```python
class VacuumBoundary(BoundaryOperator, key="vacuum"):
    def apply(self, outgoing_trace):
        # Vacuum means no incoming current/inflow trace on Gamma_-.
        # It is not a full-boundary Dirichlet condition on scalar flux.
        return self.inflow_trace_space.zero()

class ReflectiveBoundary(BoundaryOperator, key="reflective"):
    def apply(self, outgoing_trace):
        return self.reflection_permutation @ outgoing_trace

class PeriodicBoundary(BoundaryOperator, key="periodic"):
    def apply(self, outgoing_trace):
        return self.periodic_pairing @ outgoing_trace

class AlbedoBoundary(BoundaryOperator, key="albedo"):
    def apply(self, outgoing_trace):
        return self.albedo_kernel @ outgoing_trace
```

Then geometry declarations can be resolved later:

```python
@dataclass(frozen=True, slots=True)
class BoundaryDeclaration:
    patch: str
    kind: str
    parameters: dict


def resolve_boundary(decl: BoundaryDeclaration, trace_space):
    return BoundaryOperator.create(
        decl.kind,
        patch=decl.patch,
        trace_space=trace_space,
        **decl.parameters,
    )
```

This matches the architecture:

```text
Geometry declares boundary intent.
Mesh resolves boundary patches geometrically.
MethodSpace resolves boundary intent into a method-specific BoundaryOperator.
```

Do not resolve boundary declarations too early. A reflective boundary in `SN`, `PN`, `MoC`, and `MC` does not have the same representation.

---

## 5. The main hierarchy

### 5.1 Domain hierarchy

```text
Domain
├── EuclideanDomain
│   ├── IntervalDomain
│   ├── PlanarDomain
│   └── SolidDomain
├── Manifold
│   ├── Sphere
│   ├── BoundaryManifold
│   └── EnergyAxis
├── ProductDomain
├── BoundaryPhaseDomain
└── StochasticDomain
```

High-signal classes:

```python
Sphere()
EnergyInterval(E_min, E_max)
BoundaryPhaseDomain(boundary, angular_domain)
PhaseSpaceDomain(X, Omega, E, Xi)
```

### 5.2 Measure hierarchy

```text
Measure
├── VolumeMeasure
├── SurfaceMeasure
├── EnergyMeasure
├── ProbabilityMeasure
├── DiscreteMeasure
│   ├── QuadratureRule
│   ├── CubatureRule
│   ├── EmpiricalMeasure
│   └── ParticleMeasure
├── ProductMeasure
├── PushforwardMeasure
├── ConditionalMeasure
└── BoundaryCurrentMeasure
```

Important concepts:

```python
mu = mu_x * mu_omega * mu_g
mu_sphere = (mu_mu * mu_phi).pushforward(SphericalMap())
mu_boundary = BoundaryCurrentMeasure(mesh.boundary, angular_rule)
```

A quadrature or cubature rule is a finite measure:

```text
mu_h = sum_i w_i delta_{x_i}.
```

Monte Carlo particles are also an empirical measure:

```text
mu_N = sum_i w_i delta_{z_i}.
```

So deterministic quadrature and MC particles share a base abstraction.

### 5.3 Space hierarchy

```text
Space
├── FunctionSpace
│   ├── MeshFunctionSpace
│   ├── TraceSpace
│   ├── RegionSpace
│   ├── EnergyGroupSpace
│   ├── DiscreteAngularSpace
│   ├── SphericalHarmonicSpace
│   ├── PolynomialChaosSpace
│   └── CharacteristicTrackSpace
├── TensorProductSpace
├── DirectSumSpace
├── DualSpace
└── ParticleStateSpace
```

Use tensor products for independent axes:

```python
V = X * Omega * G * Xi
```

Use direct sums for coupled different fields:

```python
V_total = FluxSpace + PrecursorSpace
```

### 5.4 Basis hierarchy

```text
Basis
├── NodalBasis
├── CellwiseConstantBasis
├── PolynomialBasis
├── DiscontinuousGalerkinBasis
├── RegionIndicatorBasis
├── SphericalHarmonicBasis
├── LegendreBasis
├── PolynomialChaosBasis
└── TrackSegmentBasis
```

A basis should know:

```python
basis.evaluate(points)
basis.mass_matrix(measure)
basis.project(function, measure)
basis.reconstruct(coefficients)
```

### 5.5 Field hierarchy

```text
Field
├── CoefficientField
├── DiscreteField
├── MeshField
├── TraceField
├── AngularFlux
├── ScalarFlux
├── HarmonicMomentField
├── GroupFlux
├── StochasticField
├── TensorTrainField
└── ParticleEnsemble
```

`ParticleEnsemble` is the right name for a Monte Carlo collection representing an empirical measure. More specific MC names:

```text
ParticleHistory      one sampled path
ParticleEnsemble     weighted population at an instant or generation
SourceBank           sampled source particles
FissionBank          next-generation fission source particles
CensusBank           particles held at a time boundary
CollisionBank        particles stored at collision events
```

Use `population` when emphasizing physical branching. Use `ensemble` when emphasizing statistical representation.

### 5.6 Kernel hierarchy

```text
Kernel
├── IntegralKernel
│   ├── ScatteringKernel
│   │   ├── RotationallyInvariantScatteringKernel
│   │   ├── GeneralAnisotropicScatteringKernel
│   │   └── LegendreMomentScatteringKernel
│   ├── FissionKernel
│   ├── BoundaryKernel
│   ├── AlbedoKernel
│   ├── GreenKernel
│   └── CollisionProbabilityKernel
├── MarkovKernel
├── SubMarkovKernel
├── BranchingKernel
├── SeparableKernel
├── LowRankKernel
└── TensorNetworkKernel
```

Suffix rules:

- `Kernel` means it is integrated against a measure.
- `Operator` means it maps a field to another field.
- `Functional` means it maps a field to a scalar or response vector.
- `Projection` means it maps a field to coefficients or a coarser space.
- `Reconstruction` means it maps coefficients back to a richer space.

### 5.7 Operator hierarchy

```text
Operator
├── LinearOperator
│   ├── IdentityOperator
│   ├── ZeroOperator
│   ├── MultiplicationOperator
│   │   ├── CollisionOperator
│   │   └── TimeMassOperator
│   ├── DifferentialOperator
│   │   └── StreamingOperator
│   ├── IntegralKernelOperator
│   │   ├── ScatteringOperator
│   │   ├── FissionOperator
│   │   └── BoundaryIntegralOperator
│   ├── ProjectionOperator
│   │   ├── GalerkinProjection
│   │   ├── MomentProjection
│   │   ├── EnergyCondensationProjection
│   │   ├── SpatialHomogenizationProjection
│   │   └── StochasticGalerkinProjection
│   ├── ReconstructionOperator
│   ├── TraceOperator
│   ├── BoundaryOperator
│   ├── GreenOperator
│   ├── SweepOperator
│   ├── BlockOperator
│   ├── TensorProductOperator
│   ├── SumOfTensorProductsOperator
│   ├── TensorNetworkOperator
│   └── TensorTrainOperator
└── NonlinearOperator
    ├── BranchingGeneratingOperator
    ├── ProbabilityOfEventOperator
    └── NonlinearBoundaryOperator
```

---

## 6. Dunder methods for mathematical ergonomics

Use dunder methods only where the mathematical meaning is unambiguous.

### 6.1 Spaces

```python
V = X * Omega * G        # tensor product of spaces
W = FluxSpace + Precursors # direct sum of spaces
```

Recommended:

```python
Space.__mul__  -> tensor product
Space.__add__  -> direct sum
Space.dual()   -> dual space
```

Do not hide duality behind a cryptic dunder. Prefer `V.dual()`.

### 6.2 Measures

```python
mu = mu_x * mu_omega * mu_g
mu_sphere = mu_param.pushforward(SphericalMap())
f_param = SphericalMap().pullback(f_sphere)
```

Recommended:

```python
Measure.__mul__ -> product measure
Measure.pushforward(map)
Map.pullback(function_or_field)
```

### 6.3 Operators

```python
A + B           # operator sum
A - B           # operator difference
alpha * A       # scalar scaling
A @ B           # composition
A @ psi         # application
A & B           # tensor product operator
A.H             # Hilbert adjoint
A.T             # representation transpose
```

The distinction between `.H` and `.T` is crucial:

- `.T` is the transpose of a concrete representation.
- `.H` is the Hilbert adjoint with respect to the domain and codomain inner products.

For quadrature-weighted angular spaces, `.H` includes weights.

Test:

```python
assert_close((A @ u).inner(v), u.inner(A.H @ v))
```

### 6.4 Fields and functionals

```python
psi + phi
alpha * psi
A @ psi
R @ psi         # response functional applied to a field
psi.inner(phi)
psi.norm()
```

Avoid naked arrays where possible. A field should carry its space.

---

## 7. Geometry -> Mesh -> MethodSpace

Your stated architecture is:

```text
Geometry -> Mesh -> Method-Augmented Mesh
```

The high-signal interpretation is:

```text
GeometrySpec -> SpatialMesh -> MethodSpace
```

or, if you want to emphasize that it is a view over an existing mesh:

```text
GeometrySpec -> SpatialMesh -> MethodMeshView
```

The strongest general name is `MethodSpace`, because the augmented object is not merely a mesh. It is the space on which the method's unknowns, traces, operators, precomputations, and boundary resolutions live.

Recommended layering:

```text
GeometrySpec
    continuous, semantic, unresolved

SpatialMesh
    discrete cell complex, topology, volumes, faces, normals, region tags

MethodSpace
    method-specific representation of fields and operators over the mesh

MethodOperatorContext
    cached operators, sweep graphs, segment propagators, tensor factors, etc.
```

If the object only stores method-specific geometry caches, call it `...MeshView` or `...GeometryCache`. If it defines where method unknowns live, call it `...Space`.

### 7.1 GeometrySpec

`GeometrySpec` should include:

```python
@dataclass(frozen=True, slots=True)
class GeometrySpec:
    dimension: int
    regions: tuple[RegionSpec, ...]
    material_map: MaterialMap
    boundary_declarations: tuple[BoundaryDeclaration, ...]
    coordinate_system: CoordinateSystem
```

Geometry is not meshed. Geometry declares:

- regions
- materials
- dimensions
- coordinate systems
- boundary intent

Boundary intent is unresolved:

```python
BoundaryDeclaration(patch="outer", kind="vacuum")
BoundaryDeclaration(patch="left", kind="periodic", parameters={"pair": "right"})
BoundaryDeclaration(patch="symmetry", kind="reflective")
```

Do not instantiate boundary operators at geometry stage.

### 7.2 SpatialMesh

`SpatialMesh` should include:

```python
@dataclass(frozen=True, slots=True)
class SpatialMesh:
    cell_complex: CellComplex
    cells: tuple[Cell, ...]
    faces: tuple[Face, ...]
    boundary_faces: tuple[Face, ...]
    incidence: IncidenceRelation
    adjacency: CellAdjacency
    volumes: CellField
    face_areas: FaceField
    normals: FaceVectorField
    region_ids: CellField
    material_ids: CellField
    boundary_patches: BoundaryPatchMap
```

Mesh resolves geometry into topology and measures, but not method-specific physics.

### 7.3 MethodSpace base

```python
class MethodSpace(Space, ABC):
    mesh: SpatialMesh
    geometry: GeometrySpec
    boundary_resolution: "BoundaryResolution"

    @abstractmethod
    def resolve_boundaries(self) -> "BoundaryResolution":
        ...

    @abstractmethod
    def build_operator_context(self) -> "MethodOperatorContext":
        ...
```

### 7.4 Method-specific names for the “AugmentedMesh” class

Avoid `AugmentedMesh` in the public API. It has low signal. Use names that reveal the representation of the unknown and the native mathematics of the method.

| Method | Best public class name | If only a mesh view/cache | Why this name is high-signal |
|---|---|---|---|
| `SN` | `DiscreteOrdinatesPhaseSpace` | `OrdinatesMeshView`, `SweepGeometryCache` | unknown is a function on cells/groups/discrete directions; streaming is sweep/direct-sum by ordinate |
| `PN` | `SphericalHarmonicMomentSpace` | `MomentMeshView` | unknown is stored in spherical harmonic moments; scattering block-diagonalizes by irreps |
| `MoC` | `CharacteristicTrackSpace` | `TrackMeshView`, `RaySegmentCache` | unknown/source is propagated along characteristic tracks and segments |
| `CP` | `CollisionProbabilityRegionSpace` | `RegionMeshView`, `CPGeometryCache` | unknowns are region-integrated quantities; operator is a projected Green/collision-probability kernel |
| `MC` | `ParticleTrackingSpace` | `TrackingGeometryCache`, `DeltaTrackingCache` | method samples particle paths through geometry; can include Woodcock/delta-tracking majorants |
| Woodcock MC | `DeltaTrackingSpace` | `MajorantTrackingCache` | majorant cross section is central; virtual collisions are part of representation |
| surface tracking MC | `SurfaceTrackingSpace` | `SurfaceCrossingCache` | explicit boundary/surface intersections are central |
| stochastic Galerkin | `StochasticGalerkinPhaseSpace` | `StochasticBasisMeshView` | unknown includes polynomial-chaos or stochastic basis coefficients |
| stochastic collocation | `StochasticCollocationSpace` | `CollocationCaseEnsemble` | unknown is an ensemble of deterministic solves at random-space nodes |
| tensor train backend | `TensorTrainPhaseSpace` | `TTRepresentationCache` | state/operator lives on tensor-product axes with TT ranks |
| hybrid deterministic/MC | `ControlVariateHybridSpace` | `HybridCouplingCache` | deterministic surrogate and MC correction are coupled through residual/control variate |
| adjoint-driven hybrid | `AdjointCoupledPhaseSpace` | `ImportanceMapCache` | deterministic adjoint importance couples to MC sampling |
| symmetry sector / quotient geometry | `EquivariantSectorSpace` | `SectorMeshView`, `OrbifoldGeometryCache` | unknown lives on a fundamental domain with group-equivariant trace identifications |
| octant / reflective symmetry sector | `OctantSymmetrySectorSpace` | `OctantSectorMeshView` | sector is a fundamental domain for reflection group `(Z2)^3`; symmetry faces are quotient identifications, not physical vacuum walls |
| orbifold-reduced transport | `OrbifoldPhaseSpace` | `OrbifoldMeshView`, `StratifiedQuotientCache` | quotient has strata with stabilizers/isotropy groups; measures require orbit multiplicities |
| equivariant deterministic/MC hybrid | `EquivariantHybridPhaseSpace` | `EquivariantCouplingCache` | deterministic sector solve and MC lifting/sampling commute with group action |
| time-dependent kinetics | `EvolutionPhaseSpace` | `TimeStepCache` | unknown lives on direct/product with time or is advanced by semigroup/resolvent |
| delayed-neutron kinetics | `FluxPrecursorSpace` | `PrecursorBlockCache` | direct sum of flux and precursor spaces |

Recommended pattern:

```python
space = DiscreteOrdinatesPhaseSpace.from_mesh(mesh, angular_rule, groups)
context = space.build_operator_context()
```

The `space` is semantic. The `context` is cache/precompute.

### 7.5 MethodSpaceBuilder registry

```python
class MethodSpaceBuilder(RegistryMixin, ABC):
    registry: ClassVar[dict[str, type["MethodSpaceBuilder"]]] = {}

    @classmethod
    def _registry_base(cls):
        return MethodSpaceBuilder

    @abstractmethod
    def build(self, geometry: GeometrySpec, mesh: SpatialMesh, **kwargs) -> MethodSpace:
        ...

class SNBuilder(MethodSpaceBuilder, key="sn"):
    def build(self, geometry, mesh, *, angular_rule, groups, **kwargs):
        return DiscreteOrdinatesPhaseSpace(
            geometry=geometry,
            mesh=mesh,
            angular_rule=angular_rule,
            groups=groups,
        )

class PNBuilder(MethodSpaceBuilder, key="pn"):
    def build(self, geometry, mesh, *, order, groups, **kwargs):
        return SphericalHarmonicMomentSpace(
            geometry=geometry,
            mesh=mesh,
            order=order,
            groups=groups,
        )
```

Usage:

```python
method_space = MethodSpaceBuilder.create("sn").build(
    geometry,
    mesh,
    angular_rule=lebedev,
    groups=groups,
)
```


### 7.5 Geometry sectors, quotient meshes, and symmetry-reduced MethodSpaces

When the geometry is solved only on a sector, such as an octant, the method-augmented object should advertise that the mesh is not merely truncated. It is a representation of a quotient or fundamental domain.

Recommended names:

```text
FundamentalDomainMesh
SectorMesh
QuotientMesh
OrbifoldMesh
EquivariantSectorSpace
OrbifoldPhaseSpace
OctantSymmetrySectorSpace
```

Use `SectorMesh` when the object is primarily geometric. Use `EquivariantSectorSpace` when the object carries method unknowns and boundary resolution. Use `OrbifoldPhaseSpace` when the quotient has strata with nontrivial stabilizers, such as symmetry faces, edges, and corners.

The pipeline becomes:

```text
GeometrySpec
    -> SymmetryResolvedGeometrySpec
    -> SpatialMesh or FundamentalDomainMesh
    -> EquivariantSectorSpace(method-specific)
    -> MethodOperatorContext
```

A symmetry-aware geometry declaration should look like:

```python
@dataclass(frozen=True, slots=True)
class SymmetryDeclaration:
    group: "FiniteGroup"
    action: "GroupAction"
    sector: "FundamentalDomainSpec"
    representation: "FieldRepresentation" = TrivialRepresentation()

@dataclass(frozen=True, slots=True)
class GeometrySpec:
    dimension: int
    regions: tuple[RegionSpec, ...]
    material_map: MaterialMap
    boundary_declarations: tuple[BoundaryDeclaration, ...]
    coordinate_system: CoordinateSystem
    symmetry: SymmetryDeclaration | None = None
```

For an octant problem, the symmetry declaration should not be hidden as three reflective boundaries. The high-signal model is:

```python
symmetry = ReflectionGroup.octant()  # isomorphic to (Z2)^3
sector = FundamentalDomainSpec(bounds={"x": (0, None), "y": (0, None), "z": (0, None)})

geometry = GeometrySpec(
    dimension=3,
    regions=regions,
    material_map=materials,
    boundary_declarations=(BoundaryDeclaration("outer", "vacuum", {}),),
    coordinate_system=Cartesian(),
    symmetry=SymmetryDeclaration(group=symmetry, action=CartesianReflectionAction(), sector=sector),
)
```

Then the sector faces `x=0`, `y=0`, and `z=0` are not ordinary reflective physical walls. They are **quotient/symmetry faces** produced by the group action. They may resolve to specular reflection maps for `SN`, track continuation maps for `MoC`, parity operators for `PN`, or mirror continuation for `MC`, but their semantic origin is group quotienting.

High-signal classes:

```text
SymmetryDeclaration
FiniteGroup
ReflectionGroup
GroupAction
FundamentalDomainSpec
FundamentalDomainMesh
QuotientMesh
OrbifoldMesh
SectorBoundaryPatch
SymmetryBoundaryPatch
PhysicalBoundaryPatch
OrbitMultiplicityField
StabilizerField
IsotropyStratum
EquivariantSectorSpace
OrbifoldPhaseSpace
```

Key invariant:

```text
operator commutes with symmetry action:
A rho(g) = rho(g) A
```

or, with boundary traces,

```text
gamma_- rho(g) = rho_boundary(g) gamma_-.
```

This gives a simple but powerful implementation rule: only reduce to a sector when the geometry, materials, sources, cross sections, and responses are equivariant under the chosen group action.

---

## 8. Discrete measures, quadratures, cubatures, and symmetry

Represent every quadrature/cubature rule as a finite measure.

```python
@dataclass(frozen=True, slots=True)
class DiscreteMeasure(Measure):
    nodes: Array
    weights: Array
    domain: Domain

    def integrate(self, f):
        return sum(w * f(x) for x, w in zip(self.nodes, self.weights))

    def tensor(self, other):
        return ProductDiscreteMeasure(self, other)

    def pushforward(self, map):
        return DiscreteMeasure(
            nodes=map(self.nodes),
            weights=self.weights,
            domain=map.codomain,
        )
```

Angular product rule:

```python
mu_rule = GaussLegendre(n_mu, domain=(-1, 1))
phi_rule = GaussChebyshev(n_phi, domain=(0, 2*pi))
angular_rule = (mu_rule * phi_rule).pushforward(SphericalMap())
```

Lebedev-like rule:

```python
lebedev = InvariantSphericalCubature(
    domain=Sphere(),
    symmetry_group=OctahedralGroup(),
    exactness=SphericalHarmonicExactness(L=17),
)
```

### 8.1 Cubature registry by domain, symmetry, and exactness

```python
@dataclass(frozen=True, slots=True)
class CubatureSignature:
    domain_kind: str
    symmetry_group: str | None
    exactness_kind: str
    exactness_order: int
    positive_weights: bool

class CubatureRule(DiscreteMeasure, RegistryMixin):
    registry: ClassVar[dict[str, type["CubatureRule"]]] = {}
    signature: ClassVar[CubatureSignature]

    @classmethod
    def _registry_base(cls):
        return CubatureRule

class LebedevCubature(CubatureRule, key="lebedev"):
    signature = CubatureSignature(
        domain_kind="sphere",
        symmetry_group="octahedral",
        exactness_kind="spherical_harmonic",
        exactness_order=-1,  # depends on order argument
        positive_weights=True,
    )
```

A rule selector can then be semantic:

```python
rule = CubatureSelector.select(
    domain=Sphere(),
    exactness=SphericalHarmonicExactness(L=15),
    symmetry=OctahedralGroup(),
    positive_weights=True,
)
```

### 8.2 Exactness and invariant tests

Every angular cubature should expose:

```python
rule.assert_normalized(total=4*pi)
rule.assert_antipodal_symmetry()
rule.assert_first_moment_zero()
rule.assert_isotropic_second_moment()
rule.assert_spherical_harmonic_exactness(L)
rule.assert_group_invariant(group)
```

Important low-order moment identities on the sphere:

```text
sum_i w_i = 4 pi
sum_i w_i Omega_i = 0
sum_i w_i Omega_i Omega_i^T = (4 pi / 3) I
```

Higher-order isotropic moments should also be available as tests.

Use class names like:

```text
InvariantCubature
OrbitCubature
SphericalDesign
PositiveWeightCubature
ExactnessCertificate
MomentResidual
SymmetryResidual
```

These names tell the debugging agent what to check.

---

## 9. `SN`: Discrete ordinates as weighted angular collocation

Native structure:

```text
L^2(S^2, dOmega) -> ell^2(W)
```

where `W` is the diagonal matrix of quadrature weights.

Class:

```python
@dataclass(slots=True)
class DiscreteOrdinatesPhaseSpace(MethodSpace):
    mesh: SpatialMesh
    angular_rule: CubatureRule
    groups: EnergyGroupSpace
    boundary_resolution: BoundaryResolution
    sweep_graphs: dict[DirectionId, SweepGraph]
    trace_spaces: dict[str, TraceSpace]
```

High-signal associated classes:

```text
DiscreteAngularSpace
Ordinates
DirectionSet
AngularWeightMatrix
SweepGraph
UpwindStencil
DownwindStencil
OctantPartition
InflowTraceSpace
OutflowTraceSpace
```

Streaming in `SN`:

```text
L = sum_d D_d & Omega_d & I_g
```

or as a direct sum by ordinate:

```text
L = direct_sum_i L_{Omega_i}.
```

Scattering in `SN` should not be a dense angular matrix by default. Use harmonic projection:

```text
S_SN = R Lambda M
```

where:

```text
M = harmonic moment projection = Y^* W
R = harmonic reconstruction = Y
Lambda = Legendre/group scattering moments
```

Code:

```python
M = HarmonicMomentProjection(angular_space, L=xs.scatter.order)
R = HarmonicReconstruction(angular_space, L=xs.scatter.order)
S = R @ LegendreMomentScattering(xs.scatter) @ M
```

Boundary resolution for `SN`:

- vacuum: zero inflow trace / no incoming current on incoming faces; not zero scalar flux on the whole boundary
- reflective: discrete direction permutation or interpolation if reflected direction is not exactly in set
- periodic: face-pairing plus direction identity
- albedo: boundary integral kernel on incoming/outgoing trace ordinates

Tests:

```python
space.assert_quadrature_normalized()
space.assert_sweep_graph_acyclic_for_direction(omega)
space.assert_reflective_boundary_maps_ordinates()
S.assert_conservative()
S.assert_positive_on_nonnegative_flux()
```

---

## 10. `PN`: Spherical harmonics as irreducible moment space

Native structure:

```text
L^2(S^2) decomposes into SO(3) irreps H_ell.
```

```text
P_N = direct_sum_{ell=0}^N H_ell.
```

Class:

```python
@dataclass(slots=True)
class SphericalHarmonicMomentSpace(MethodSpace):
    mesh: SpatialMesh
    order: int
    harmonics: SphericalHarmonicSpace
    groups: EnergyGroupSpace
    moment_index: HarmonicIndex
    coupling_matrices: StreamingMomentCouplings
```

High-signal classes:

```text
SphericalHarmonicSpace
HarmonicIndex
IrrepBlock
MomentVector
MomentMassMatrix
StreamingMomentCoupling
RotationallyInvariantDiagonalization
```

The scattering operator is diagonal/block-diagonal by `ell`:

```text
S Y_{ell m} = Sigma_{s,ell} Y_{ell m}.
```

This is representation theory. A rotationally invariant kernel commutes with rotations, so by Schur's lemma it acts as a scalar on each irreducible block.

Code:

```python
S = RotationallyInvariantScattering(xs.scatter).diagonalize(harmonic_space)
```

Boundary resolution for `PN` is subtle:

- vacuum/inflow boundary becomes a weak boundary condition or closure
- Marshak-like conditions are moment boundary closures
- reflective boundary becomes a parity/sign transformation on moments
- albedo is a boundary moment kernel

Use names like:

```text
MomentBoundaryClosure
MarshakBoundaryClosure
ParityReflectionOperator
HalfRangeMomentProjection
```

Tests:

```python
harmonics.assert_orthonormal()
S.assert_block_diagonal_by_ell()
S.assert_commutes_with_rotation(samples=...)
closure.assert_well_posed()
```

---

## 11. `MoC`: Method of characteristics as characteristic flow and ray bundles

Native structure:

```text
Streaming is the Lie derivative along V = (Omega, 0).
Characteristics are flow lines of dx/ds = Omega.
```

Class:

```python
@dataclass(slots=True)
class CharacteristicTrackSpace(MethodSpace):
    mesh: SpatialMesh
    angular_rule: CubatureRule
    ray_bundles: tuple[RayBundle, ...]
    segment_index: SegmentIndex
    track_measure: TrackMeasure
    boundary_resolution: BoundaryResolution
```

Do not call this object `Billiard`. Billiard dynamics is appropriate only when reflective boundary conditions dominate the dynamics. MoC itself is characteristic transport; reflection is one possible boundary map.

High-signal classes:

```text
CharacteristicFlow
RayBundle
TrackFamily
Track
Segment
SegmentPropagator
OpticalLength
AttenuationFactor
TrackMeasure
TransverseMeasure
CyclicTrackCoupling
```

Segment solve:

```text
psi_out = exp(-tau) psi_in + (1 - exp(-tau)) / Sigma_t * q
where tau = Sigma_t * segment_length.
```

Code:

```python
prop = SegmentPropagator(segment, xs.total)
psi_out = prop(psi_in, segment_source)
```

MoC volume preservation:

```text
sum_{tracks r} sum_{segments s in cell c} w_r length_s ≈ volume_c.
```

Test:

```python
tracks.assert_volume_preserving(mesh)
tracks.assert_region_volume_preserving(mesh.regions)
tracks.assert_boundary_connectivity()
```

Boundary resolution for MoC:

- vacuum: terminate or zero inflow at boundary segment
- reflective: track continuation through specular reflection map
- periodic: track continuation through boundary identification
- albedo: sample or integrate boundary response into incoming track flux

Represent boundary behavior as a `BoundaryMap` acting on track endpoints, not as something baked into tracks forever. That allows changing boundary conditions without regenerating all track geometry when possible.

---


## 11A. `MoC` solution sheaves, affine local solutions, and track holonomy

The `MoC` section above identifies the geometric primitive: characteristic flow and ray bundles. There is an additional native local-to-global structure that is especially high signal for `MoC`: the **characteristic transport sheaf**.

Do **not** generalize this into `EverythingIsASheaf`. The reason to use sheaf vocabulary here is specific:

```text
MoC = local characteristic segment solutions + trace restrictions + interface/boundary gluing.
```

That is precisely the pattern of a solution sheaf.

### 11A.1 Continuous idea

For an open subset `U` of the characteristic track complex, define

```text
Sol(U) = { local angular fluxes satisfying the transport equation on U }.
```

If `V subset U`, the local solution restricts:

```text
rho_UV: Sol(U) -> Sol(V).
```

The sheaf laws are:

```text
locality: if two local solutions agree on every local piece, they agree globally.
gluing: compatible local solutions glue to a global solution.
```

For a transport code, this becomes operational:

```text
local segment solve -> endpoint traces -> interface and boundary compatibility -> global track solution.
```

### 11A.2 Affine solution sheaf, not merely linear solution sheaf

For homogeneous transport on a segment,

```text
Omega · grad psi + Sigma_t psi = 0,
```

the local solution space is linear. For inhomogeneous transport,

```text
Omega · grad psi + Sigma_t psi = q,
```

the local solution space is affine.

For a segment `e`,

```text
psi_out_e = a_e psi_in_e + b_e.
```

where

```text
a_e = exp(-tau_e)
tau_e = integral_e Sigma_t ds
b_e = integral_0^ell exp(-integral_s^ell Sigma_t(r) dr) q(s) ds.
```

For constant source and constant total cross section on the segment:

```text
psi_out = exp(-Sigma_t ell) psi_in + (1 - exp(-Sigma_t ell)) / Sigma_t * q.
```

So the high-signal class should be:

```python
AffineSolutionSheaf
CharacteristicTransportSheaf
```

not just `SolutionSheaf`.

Recommended segment primitive:

```python
@dataclass(frozen=True, slots=True)
class SegmentAffinePropagator:
    segment: Segment
    attenuation: float
    source_response: float

    def propagate(self, psi_in: float) -> float:
        return self.attenuation * psi_in + self.source_response
```

For tensor-network composition, represent affine maps using homogeneous coordinates:

```text
[ psi_out ]   [ a_e  b_e ] [ psi_in ]
[   1     ] = [  0    1  ] [   1    ].
```

High-signal tensor names:

```text
SegmentAffineTensor
SegmentPropagatorTensor
HomogeneousAffineSegmentTensor
TrackTransferTensor
```

### 11A.3 Track complex and stalks

A track is a 1D cell complex:

```text
vertex -> segment -> vertex -> segment -> ... -> vertex.
```

Use explicit topology:

```python
@dataclass(frozen=True, slots=True)
class TrackVertex:
    id: int
    kind: str  # inflow, outflow, interface, reflection, periodic, symmetry

@dataclass(frozen=True, slots=True)
class TrackSegment:
    id: int
    cell_id: int
    start: TrackVertex
    end: TrackVertex
    length: float
    optical_length: float

@dataclass(frozen=True, slots=True)
class CharacteristicTrackComplex:
    vertices: tuple[TrackVertex, ...]
    segments: tuple[TrackSegment, ...]
```

The sheaf assigns:

```text
segment stalk  -> local affine solution data on a segment
vertex stalk   -> trace space at an endpoint/interface/boundary hit
restriction    -> take local segment solution to its incoming or outgoing endpoint trace
```

Recommended protocol:

```python
class Presheaf(Protocol):
    def stalk(self, cell): ...
    def restrict(self, section, from_cell, to_cell): ...

class HasGluingResidual(Protocol):
    def gluing_residual(self, local_sections): ...

@dataclass(frozen=True, slots=True)
class CharacteristicTransportSheaf:
    track_complex: CharacteristicTrackComplex
    segment_propagators: Mapping[int, SegmentAffinePropagator]
    boundary_gluing: "BoundaryGluingNetwork"

    def restrict_to_inflow_trace(self, segment_id: int, section): ...
    def restrict_to_outflow_trace(self, segment_id: int, section): ...
    def propagate_segment(self, segment_id: int, psi_in, source): ...
    def gluing_residual(self, local_sections): ...
    def global_sections(self): ...
```

### 11A.4 Gluing residual as the central diagnostic

At an interior vertex, the outgoing trace of the upstream segment must equal the incoming trace of the downstream segment:

```text
gamma_+ psi_e = gamma_- psi_e_next.
```

The mismatch is a gluing residual:

```text
r_v = gamma_+ psi_e - gamma_- psi_e_next.
```

Expose this explicitly:

```python
sheaf.gluing_residual(local_sections)
sheaf.assert_glued(local_sections)
sheaf.assert_trace_consistency(local_sections)
sheaf.assert_boundary_compatible(local_sections)
```

This catches method bugs that ordinary matrix residuals often obscure:

```text
wrong segment ordering
wrong endpoint orientation
wrong incoming/outgoing sign convention
wrong reflective map
wrong periodic pairing
wrong symmetry-sector multiplicity
wrong track continuation after boundary hit
wrong material assignment on segment
wrong optical length
```

### 11A.5 Vacuum in MoC sheaf language

Vacuum is a boundary constraint on incoming endpoint trace:

```text
gamma_- psi = 0 on Gamma_-.
```

For `MoC`:

```text
incoming track endpoint value = 0.
```

High-signal names:

```text
ZeroInflowEndpointConstraint
IncomingTrackEndpoint
OutgoingTrackEndpoint
CharacteristicBoundaryTrace
```

Avoid names like `ZeroFluxEndpoint` because the outgoing endpoint is not zero and leakage is allowed.

### 11A.6 Reflective and periodic boundaries create monodromy

With reflective or periodic boundary maps, a track may close into a cycle. The local affine maps compose around the cycle:

```text
psi = A_cycle psi + b_cycle.
```

The cycle constraint is:

```text
(I - A_cycle) psi = b_cycle.
```

This is the track **monodromy** or **holonomy**.

Recommended classes:

```text
TrackMonodromy
CharacteristicHolonomy
CyclePropagator
AffineCycleConstraint
PathGroupoidRepresentation
```

For attenuating media, a closed homogeneous cycle should satisfy:

```text
rho(A_cycle) < 1
```

unless the physics intentionally permits noncontracting gain. This is a very strong diagnostic.

Tests:

```python
track.assert_cycle_closed()
track.assert_monodromy_contracting()
track.assert_affine_cycle_consistent()
track.assert_boundary_gluing_consistent()
```

### 11A.7 Cochain view: global sections as `H0`

On a finite track complex, the sheaf can expose a small cochain complex:

```text
C0 = local segment/vertex assignments
C1 = interface mismatch assignments
d  = trace coboundary / gluing residual
```

Then

```text
H0 = ker(d)
```

is the space of globally glued characteristic solutions.

Practical names:

```text
TraceCoboundary
GluingCochain
CocycleResidual
CohomologyDiagnostic
ObstructionClass
```

This is especially useful for reflective, periodic, and symmetry-sector tracks, where failure to close a cycle is not just a numerical residual but an obstruction to a globally consistent section.

### 11A.8 Sheaf versus cosheaf

Use sheaves for compatibility of local data:

```text
TraceSheaf
BoundaryGluingSheaf
CharacteristicTransportSheaf
EquivariantSolutionSheaf
```

Use cosheaves for additive aggregation:

```text
TallyCosheaf
CurrentCosheaf
ReactionRateCosheaf
RegionMomentCosheaf
EmpiricalMeasureCosheaf
```

The distinction is high signal:

```text
sheaf     -> restrict, compare, glue, find global section
cosheaf   -> extend, aggregate, sum, conserve
```

## 12. `CP`: Collision probabilities as projected Green kernels

Native structure:

```text
CP is a Green-kernel method projected onto region basis functions.
```

Class:

```python
@dataclass(slots=True)
class CollisionProbabilityRegionSpace(MethodSpace):
    mesh: SpatialMesh
    regions: RegionSpace
    region_basis: RegionIndicatorBasis
    optical_geometry: OpticalRegionGeometry
    escape_surfaces: BoundaryTraceSpace
```

High-signal classes:

```text
RegionSpace
RegionIndicatorBasis
CollisionProbabilityKernel
EscapeProbability
TransmissionProbability
GreenKernel
ResponseMatrix
ReciprocityResidual
ConservationDefect
```

The CP matrix is not a generic dense matrix. It is a positive, conservative/subconservative, reciprocal kernel under appropriate assumptions.

Important laws:

```text
P_ij >= 0
sum_j P_ij + P_i_escape = 1
V_i Sigma_i P_ij = V_j Sigma_j P_ji   # reciprocity, when assumptions hold
```

Tests:

```python
P.assert_positive()
P.assert_conservative(include_escape=True)
P.assert_reciprocal(volumes, sigmas)
```

Compute as projection:

```python
G = CollisionGreenOperator(geometry, xs.total)
R = RegionReconstruction(regions)
P = R.H @ G @ R
```

Again the formula is Galerkin projection:

```text
A_h = P^* A P.
```

---

## 13. `MC`: Monte Carlo as Markov and branching processes

Native structure:

```text
A particle history is a sample path of a Markov/branching process on phase space.
```

Class:

```python
@dataclass(slots=True)
class ParticleTrackingSpace(MethodSpace):
    mesh: SpatialMesh
    tracking_geometry: TrackingGeometry
    boundary_resolution: BoundaryResolution
    material_locator: MaterialLocator
    collision_kernel: MarkovKernel
    fission_kernel: BranchingKernel
```

For Woodcock / delta tracking:

```python
@dataclass(slots=True)
class DeltaTrackingSpace(ParticleTrackingSpace):
    majorant_cross_section: MajorantField
    virtual_collision_kernel: MarkovKernel
```

Use `majorant`, not `majorand`. In analysis and probability, a majorant is the dominating function.

High-signal classes:

```text
ParticleState
ParticleHistory
ParticleEnsemble
SourceBank
FissionBank
CensusBank
HazardRate
FreeFlightSampler
SurfaceCrossingSampler
CollisionKernel
ScatteringKernel
BranchingFissionKernel
ImplicitCapture
RussianRoulette
Splitting
WeightWindow
TrackLengthEstimator
CollisionEstimator
ExpectedValueEstimator
```

A particle ensemble is an empirical measure:

```text
mu_N = sum_i w_i delta_{z_i}.
```

A tally is a functional:

```python
R = ReactionRateFunctional(xs.absorption, region="fuel")
score = R.estimate_from(histories)
```

Boundary resolution for MC:

- vacuum: kill particle at boundary
- reflective: update direction by specular reflection
- periodic: map position to paired boundary
- albedo: sample reflected/scattered outgoing boundary state with weight correction

Tests:

```python
collision_kernel.assert_markov()
fission_kernel.assert_branching_mean_matches(F)
tracking_space.assert_majorant_dominates(xs.total)
estimator.assert_unbiased_on_manufactured_case()
```

---

## 14. Stochastic PDE / stochastic transport spaces

There are two major stochastic meanings.

### 14.1 Random-input stochastic PDE / uncertainty quantification

Cross sections, densities, sources, geometry, or boundary conditions depend on a random variable:

```text
xi in Xi.
```

The unknown is

```text
psi(x, Omega, g, xi).
```

Native space:

```text
V = X * Omega * G * Xi.
```

Classes:

```text
StochasticSpace
ProbabilityMeasure
PolynomialChaosSpace
StochasticGalerkinPhaseSpace
StochasticCollocationSpace
RandomField
RandomCoefficientExpansion
TripleProductTensor
```

Stochastic Galerkin:

```python
Xi = PolynomialChaosSpace(probability_measure, order=p)
V = X * Omega * G * Xi
A_sg = StochasticGalerkinProjection(Xi).apply(A)
```

If

```text
A(xi) = sum_k A_k theta_k(xi),
```

then stochastic Galerkin produces

```text
A_{beta, alpha} = sum_k A_k E[Psi_beta theta_k Psi_alpha].
```

The tensor

```text
C[beta, k, alpha] = E[Psi_beta theta_k Psi_alpha]
```

should be named:

```python
TripleProductTensor
```

### 14.2 Branching stochastic transport / probability of event

Here the equation may be nonlinear because fission branches.

Classes:

```text
BranchingProcess
GeneratingFunctional
ProbabilityOfEventEquation
NonlinearTransportOperator
ExtinctionProbability
InitiationProbability
```

Do not force these into `LinearOperator`. Use `NonlinearOperator` when required.

---

## 15. Tensor products, tensor networks, and tensor trains

The deterministic multigroup problem already lives on a tensor-product space:

```text
V = X * Omega * G.
```

The stochastic problem extends this:

```text
V = X * Omega * G * Xi_1 * ... * Xi_d.
```

The native operator form is often a sum of tensor products:

```text
A = sum_k A_x^(k) ⊗ A_omega^(k) ⊗ A_g^(k) ⊗ A_xi^(k).
```

Use explicit classes:

```text
TensorProductOperator
SumOfTensorProductsOperator
TensorNetworkOperator
TensorTrainOperator
MatrixProductOperator
```

### 15.1 Streaming as tensor product

For `SN`:

```text
L = D_x ⊗ Omega_x ⊗ I_g
  + D_y ⊗ Omega_y ⊗ I_g
  + D_z ⊗ Omega_z ⊗ I_g.
```

Code:

```python
L = (Dx & Omega.x & Ig) + (Dy & Omega.y & Ig) + (Dz & Omega.z & Ig)
```

### 15.2 Scattering as sum of tensor products

For Legendre scattering moments:

```text
S = sum_ell M_{Sigma_s,ell}(x) ⊗ A_ell ⊗ G_ell.
```

Code:

```python
S = SumOfTensorProducts([
    SigmaMoment(xs.scatter, ell) & AngularMomentOperator(ell) & GroupScatteringMatrix(xs.scatter, ell)
    for ell in range(xs.scatter.order + 1)
])
```

### 15.3 Tensor train role

Tensor train is not the native representation of ordinary `PN` scattering. The native representation there is spherical harmonics / irreducible representation blocks.

Tensor train becomes natural when:

- stochastic dimension is high
- energy-group structure is large
- random fields are represented by many variables
- operators are high-dimensional but separable or approximately separable
- one wants compressed solvers over `X * Omega * G * Xi_1 * ... * Xi_d`

Recommended flow:

```text
physics operator
    -> exact SumOfTensorProducts
        -> optional TensorTrain compression
            -> TT solver
```

Do not convert physics to a dense tensor and then compress. Preserve the known separability first.

---


## 15A. Native local-to-global structures by method

The `MoC` sheaf idea should not be universalized. Each method has its own native way of turning local information into global information. The architecture should name these structures directly, because the name tells Codex and human developers which invariants and tests apply.

The common conceptual protocol is:

```python
class HasCompatibilityResidual(Protocol):
    def compatibility_residual(self, candidate): ...

class HasLocalToGlobalAssembly(Protocol):
    def assemble_global(self, local_data): ...

class HasConservationLaw(Protocol):
    def balance_residual(self, state): ...

class HasKernelComposition(Protocol):
    def compose_kernel(self, other): ...

class HasTensorNetwork(Protocol):
    def as_tensor_network(self): ...
```

But do not force all methods to inherit from a concrete `LocalToGlobalBase`. Prefer structural typing via `Protocol` and method-native names.

### 15A.1 Method-native structures table

| Method / feature | Native local-to-global structure | High-signal classes |
|---|---|---|
| `MoC` | affine solution sheaf + path groupoid + holonomy | `CharacteristicTransportSheaf`, `TrackMonodromy`, `PathGroupoidRepresentation` |
| `SN` | upwind trace complex + sweep dependency graph | `UpwindTraceComplex`, `SweepDependencyGraph`, `CausalTransportDAG` |
| `PN` | representation module + irreducible harmonic decomposition | `SphericalHarmonicModule`, `IrrepDecomposition`, `MomentHierarchy` |
| `SPN` | asymptotic moment reduction + diffusion-limit complex | `AsymptoticMomentReduction`, `DiffusionLimitComplex` |
| `CP` | Green kernel algebra + region cosheaf | `CollisionProbabilityKernel`, `GreenPotentialOperator`, `RegionMomentCosheaf` |
| `MC` | Markov/branching process + filtration + empirical measure | `TransportMarkovProcess`, `BranchingTransportProcess`, `WeightedParticleEnsemble` |
| Tallies | cosheaf / additive functional | `TallyCosheaf`, `AdditiveFunctional`, `ResponseFunctional` |
| SPDE/UQ | Bochner space + chaos algebra | `BochnerSpace`, `PolynomialChaosAlgebra`, `StochasticGalerkinTensor` |
| Collocation | probability cubature / scenario ensemble | `ProbabilityCubature`, `SparseGridCollocation`, `ScenarioEnsemble` |
| Tensor train/network | factor graph + index contraction | `TransportFactorGraph`, `TransportTensorNetwork`, `MatrixProductOperator` |
| Symmetry sector | groupoid / equivariant bundle / orbifold compatibility complex | `TransportGroupoid`, `EquivariantBundle`, `OrbifoldPhaseSpace` |
| Energy condensation | projection adjunction / conditional expectation | `CondensationAdjunction`, `EnergyConditionalExpectation`, `ReactionRatePreservingProjection` |
| Homogenization | quotient projection + moment preservation | `SpatialQuotientMap`, `FluxWeightedHomogenization`, `RegionAggregationMap` |
| Time dependence | semigroup + generator + resolvent | `TransportSemigroup`, `TransportGenerator`, `ResolventOperator` |
| `k` eigenvalue | positive spectral theory / next-generation operator | `NextGenerationOperator`, `PerronFrobeniusMode`, `DominanceRatio` |
| `alpha` eigenvalue | semigroup spectrum / Laplace resolvent | `AlphaEigenproblem`, `LaplaceTransformedTransportOperator`, `ResolventSpectrum` |
| Conservation | chain complex / balance cochain | `MeshChainComplex`, `CurrentCochain`, `BalanceResidual` |

### 15A.2 `SN`: upwind trace complex and sweep graph

For `SN`, the angular space is a discrete measure, while the local-to-global spatial structure is causal.

For each direction `Omega_i`, each face is classified by sign:

```text
Omega_i · n_f < 0  -> inflow face
Omega_i · n_f > 0  -> outflow face
```

This induces a directed dependency graph over cells. If acyclic, a topological order gives a sweep. If reflective or periodic boundaries create cycles, the graph decomposes into strongly connected components.

High-signal classes:

```text
UpwindTraceComplex
OrdinatesFaceTraceSystem
SweepDependencyGraph
CausalTransportDAG
DirectionSweepOrdering
SweepStrongComponent
CyclicSweepBlock
ReflectiveSweepCycle
```

Tests:

```python
sweep_graph.assert_upwind_orientation()
sweep_graph.assert_topologically_sorted()
sweep_graph.assert_face_pairing_consistent()
sweep_graph.assert_boundary_trace_classification()
sweep_graph.assert_cycles_are_declared()
```

### 15A.3 `PN`: representation module and harmonic bundle

For `PN`, angular global structure is representation theory, not sheaf gluing. The angular space decomposes as

```text
P_N = direct_sum_{ell=0}^N H_ell
```

where each `H_ell` is an irreducible representation of `SO(3)` of dimension `2 ell + 1`.

Rotationally invariant scattering is block diagonal by `ell` by Schur's lemma. Streaming couples only allowed angular-momentum blocks, typically `ell -> ell +/- 1`.

High-signal names:

```text
SphericalHarmonicModule
SO3Module
IrrepBlock
IrrepDecomposition
MomentHierarchy
ClebschGordanCoupling
SchurDiagonalScattering
ParityDecomposition
EvenOddMomentSplit
HarmonicTraceOperator
MomentClosure
```

Tests:

```python
harmonics.assert_orthonormal()
S.assert_block_diagonal_by_ell()
S.assert_so3_equivariant()
streaming.assert_adjacent_ell_coupling()
moment_space.assert_parity_structure()
```

For symmetry sectors, `PN` should use an equivariant harmonic bundle:

```text
EquivariantHarmonicBundle
HarmonicGroupAction
ReflectionRepresentationMatrix
OrbifoldMomentSpace
```

### 15A.4 `SPN`: asymptotic reduction

`SPN` is best named as an asymptotic reduction of the moment hierarchy, not as a generic low-order `PN` hack.

High-signal names:

```text
AsymptoticMomentReduction
DiffusionLimitComplex
EvenParityMomentSystem
OddParityClosure
SPNOperatorBlock
SPNMarshakBoundarySystem
```

Tests:

```python
spn.assert_diffusion_limit()
spn.assert_even_parity_closure()
spn.assert_boundary_layer_correction_consistent()
```

### 15A.5 `CP`: Green kernel algebra and region cosheaf

`CP` globalizes by integrating Green kernels over regions. It is cosheaf-like because contributions aggregate from local geometry to region moments.

Native objects:

```text
GreenKernelAlgebra
GreenPotentialOperator
CollisionProbabilityKernel
EscapeProbabilityKernel
RegionMomentCosheaf
SubMarkovCollisionKernel
ReciprocityForm
ConservativeRegionBalance
```

Tests:

```python
cp.assert_positive()
cp.assert_submarkov()
cp.assert_conservative_with_escape()
cp.assert_reciprocal()
cp.assert_region_balance()
cp.assert_escape_probability_bounds()
```

### 15A.6 `MC`: Markov process, branching process, and empirical measure

`MC` globalizes by composing kernels into path measures, not by gluing local solution patches.

Native objects:

```text
TransportMarkovProcess
BranchingTransportProcess
PathMeasure
HistoryFiltration
WeightedParticleEnsemble
EmpiricalPhaseMeasure
TallyCosheaf
AdditiveFunctional
ScoreKernel
```

Use `WeightedParticleEnsemble` for the in-memory bank and `EmpiricalPhaseMeasure` for its mathematical interpretation.

A fission population is high-signal as:

```text
BranchingParticleCloud
FissionBank
NeutronPopulationMeasure
```

depending on whether the code is emphasizing branching, implementation, or measure-valued interpretation.

Tests:

```python
kernel.assert_markov_normalization()
branching.assert_expected_offspring_matches_nu_sigma_f()
estimator.assert_unbiased_on_manufactured_case()
population.assert_bank_balance()
tally.assert_score_conservation()
```

### 15A.7 SPDE/UQ: Bochner space and chaos algebra

For random-input transport,

```text
psi = psi(x, Omega, g, xi)
```

lives in a Bochner-type tensor product space:

```text
L2(Xi; V_transport).
```

Native objects:

```text
BochnerSpace
BochnerField
RandomTransportOperator
PolynomialChaosSpace
PolynomialChaosAlgebra
TripleProductTensor
ConditionalExpectationProjection
MeanProjection
VarianceFunctional
```

Tests:

```python
chaos.assert_orthonormal()
triple_product.assert_symmetry()
projection.assert_mean_preserving()
projection.assert_variance_consistent()
sg_operator.assert_deterministic_limit()
```

### 15A.8 Stochastic collocation: probability cubature

Collocation is cubature in probability space:

```text
E[f(xi)] ≈ sum_i w_i f(xi_i).
```

Native names:

```text
ProbabilityCubature
SparseGridMeasure
CollocationNodeSet
ScenarioEnsemble
ScenarioSolve
```

Tests:

```python
cubature.assert_probability_normalized()
cubature.assert_polynomial_exactness()
ensemble.assert_weighted_mean_consistent()
ensemble.assert_weighted_variance_consistent()
```

### 15A.9 Tensor networks: factor graph and monoidal composition

Tensor networks globalize by index contraction. This is the correct abstraction when local factors are not local PDE solutions but algebraic tensors:

```text
local tensors + shared bond indices -> global operator/state by contraction.
```

Native names:

```text
TransportFactorGraph
TransportTensorNetwork
TraceBondSpace
TensorProductOperator
SumOfTensorProducts
MatrixProductOperator
TensorTrainState
TensorTrainOperator
BoundaryTensorNetwork
SectorTensorNetwork
```

Tests:

```python
tn.assert_bond_dimensions()
tn.assert_contraction_matches_reference()
tt.assert_rounding_error_bound()
tt.assert_operator_residual()
tt.assert_gauge_consistency()
```

### 15A.10 Conservation: chain complex and balance cochains

For finite-volume-like discretizations, conservation belongs to a chain complex:

```text
cells <- faces <- edges <- vertices
```

with incidence operators satisfying

```text
boundary(boundary(.)) = 0.
```

Native names:

```text
MeshChainComplex
CellChainComplex
FaceIncidenceMatrix
DiscreteDivergence
CurrentCochain
BalanceCochain
ConservationResidual
```

Tests:

```python
mesh_complex.assert_boundary_of_boundary_zero()
balance.assert_cellwise_conservative()
current.assert_face_orientation_consistent()
```

This is often a better debugging tool than inspecting matrix rows, because it directly exposes orientation and leakage mistakes.

## 16. Tensor-network boundary conditions

Boundary conditions are operators on trace spaces.

Define inflow and outflow boundary phase spaces:

```text
Gamma_- = {(x, Omega) on boundary : Omega · n(x) < 0}
Gamma_+ = {(x, Omega) on boundary : Omega · n(x) > 0}
```

The trace operators are:

```text
gamma_- psi = incoming trace
gamma_+ psi = outgoing trace.
```

A general linear boundary condition can be written as

```text
gamma_- psi = B gamma_+ psi + g.
```

The boundary operator

```text
B: Trace_+ -> Trace_-
```

should separate geometry from response.

### 16.1 Geometry-response factorization

Represent boundary response as

```text
B = ResponseKernel @ GeometryMap
```

or, more explicitly,

```text
B = R_boundary @ G_boundary

where
    G_boundary: Trace_+ -> GeometricallyMappedTrace
    R_boundary: GeometricallyMappedTrace -> Trace_-
```

or, more generally,

```text
B = TensorNetwork(
    geometry factors,
    angular response factors,
    energy response factors,
    stochastic response factors,
    patch coupling factors,
)
```

Geometry factors know:

- patch incidence
- face pairing
- normals
- specular reflection direction map
- periodic coordinate map
- trace index maps

Response factors know:

- angular redistribution
- energy redistribution
- albedo strength
- stochastic dependence
- material response

High-signal classes:

```text
BoundaryTraceSpace
InflowTraceSpace
OutflowTraceSpace
BoundaryGeometryMap
BoundaryResponseKernel
BoundaryTensorNetwork
PatchIncidenceFactor
SpecularReflectionMap
PeriodicPairingMap
AlbedoAngularKernel
AlbedoEnergyKernel
```

### 16.2 Examples

Vacuum:

```text
B = 0
```

This means

```text
gamma_- psi = 0
```

on the incoming trace space. In physical words, vacuum means **no incoming current from outside the domain**. It does **not** mean zero scalar flux on the boundary, and it does **not** mean a trivial full-boundary Dirichlet condition. Outgoing particles may freely leave through `Gamma_+`; only the external inflow contribution on `Gamma_-` is zero.

Code:

```python
bc = VacuumBoundary.on("outer")
assert bc.is_zero_inflow_trace
assert not bc.is_dirichlet_scalar_flux
```

Reflective:

```text
B = P_specular
```

where `P_specular` is a geometry-induced permutation/interpolation between outgoing and incoming trace states.

Code:

```python
bc = ReflectiveBoundary.on("symmetry")
B = SpecularReflectionMap(mesh.boundary_patch("symmetry"), angular_space)
```

Periodic:

```text
B = P_periodic
```

Code:

```python
bc = PeriodicBoundary(pair=("left", "right"))
```

Albedo:

```text
B = G_patch ⊗ K_omega ⊗ K_g
```

where:

- `G_patch` maps boundary geometry
- `K_omega` redistributes angle
- `K_g` redistributes energy

Code:

```python
B = BoundaryTensorNetwork([
    PatchIncidenceFactor(patch),
    BoundaryGeometryMap(patch),
    AlbedoAngularKernel(...),
    AlbedoEnergyKernel(...),
])
```

This makes it possible to swap geometry without changing response, or swap response without changing geometry.

### 16.3 Boundary composition

Use direct sums over disjoint boundary patches:

```python
bc = (
    VacuumBoundary.on("outer")
    + ReflectiveBoundary.on("symmetry")
    + PeriodicBoundary.on("left", pair="right")
)
```

Use composition for layered responses:

```python
bc = BoundaryResponseKernel(...) @ BoundaryGeometryMap(...)
```

Use tensor product for separable response:

```python
bc = PatchMap & AngularAlbedo & EnergyAlbedo
```


### 16.4 Vacuum boundary semantics

Vacuum boundary conditions must be represented as a zero **inflow trace** condition:

```text
gamma_- psi = 0.
```

Equivalently, the incoming current from outside the domain is zero:

```text
J_in = integral_{Gamma_-} |Omega · n| psi dGamma = 0.
```

This is not the same as imposing

```text
psi = 0 on Gamma
```

and it is not the same as imposing

```text
scalar_flux = 0 on boundary.
```

The outgoing trace `gamma_+ psi` is unconstrained by vacuum and contributes to leakage. The distinction matters architecturally because a vacuum boundary operator has codomain `Trace_-`, not the full boundary field space.

High-signal implementation names:

```text
ZeroInflowBoundary
VacuumBoundary
InflowTraceZero
IncomingCurrentZero
LeakageTraceUnconstrained
```

Avoid names like:

```text
ZeroFluxBoundary
DirichletBoundary
```

unless the code truly means a Dirichlet condition on a scalar diffusion-like field.

Recommended flags/properties:

```python
bc.is_zero_inflow_trace == True
bc.is_zero_scalar_flux == False
bc.allows_outgoing_leakage == True
```

---


## 16A. Boundary conditions as trace laws, response kernels, and method realizers

The clean architecture for boundary conditions is a three-layer decomposition:

```text
Boundary condition
    = geometric trace structure
    + physical response law
    + method-specific realization.
```

In operator form, the most common linear boundary law is

```text
gamma_- psi = B_boundary gamma_+ psi + q_boundary.
```

The boundary map itself should usually factor as

```text
B_boundary = R_boundary @ G_boundary
```

where

```text
G_boundary: outgoing trace -> geometrically mapped outgoing trace
R_boundary: geometrically mapped outgoing trace -> incoming trace
q_boundary: prescribed incoming source
```

Thus

```text
gamma_- psi = R_boundary G_boundary gamma_+ psi + q_boundary.
```

This decomposition is real and useful. But the method-specific enforcement can be exact, approximate, weak, projected, or stochastic depending on the method.

### 16A.1 Boundary primitive hierarchy

Use these high-signal layers:

```text
BoundaryDeclaration       # semantic declaration at Geometry stage
BoundaryPatch             # resolved boundary piece at Mesh stage
TraceSpace                # incoming/outgoing trace spaces
BoundaryGeometryMap       # geometry-induced map on trace points
BoundaryResponseKernel    # physical response on mapped traces
BoundarySource            # prescribed incoming source
BoundaryTraceLaw          # affine law gamma_- = R G gamma_+ + q
BoundaryRelation          # more general constraint A_- gamma_- + A_+ gamma_+ = g
BoundaryRealizer          # method-specific realization
MethodBoundaryOperator    # concrete operator/constraint/kernel/tensor for one method
```

Recommended hierarchy:

```text
BoundaryCondition
    BoundaryRelation
        AffineTraceLaw
            ZeroInflowTrace
            VacuumInflow
            PrescribedInflow
            ReflectiveBoundary
            PeriodicBoundary
            AlbedoBoundary
            SymmetryBoundary
        WeakBoundaryForm
            MarshakBoundaryForm
            MarkBoundaryForm
            VariationalBoundaryForm
            PenaltyBoundaryForm
```

`AffineTraceLaw` is sufficient for many transport boundary conditions. `BoundaryRelation` is needed for weak `PN`, variational, penalty, and mixed forms.

### 16A.2 Protocol and dataclass skeleton

```python
class BoundaryGeometryMap(Protocol):
    def map_outgoing_to_incoming_geometry(self, outgoing_trace_point): ...
    def as_tensor(self, trace_space): ...

class BoundaryResponseKernel(Protocol):
    def apply(self, mapped_outgoing_trace): ...
    def as_tensor(self, trace_space): ...

class BoundarySource(Protocol):
    def evaluate(self, incoming_trace_point): ...

@dataclass(frozen=True, slots=True)
class BoundaryTraceLaw:
    geometry_map: BoundaryGeometryMap
    response: BoundaryResponseKernel
    source: BoundarySource

    def as_operator(self, trace_space):
        G = self.geometry_map.as_operator(trace_space)
        R = self.response.as_operator(trace_space)
        q = self.source.as_field(trace_space.incoming)
        return AffineBoundaryOperator(linear=R @ G, source=q)
```

A physically precise vacuum boundary is:

```python
@dataclass(frozen=True, slots=True)
class ZeroInflowTrace(BoundaryTraceLaw):
    def __init__(self):
        super().__init__(
            geometry_map=IdentityBoundaryMap(),
            response=ZeroBoundaryResponse(),
            source=ZeroBoundarySource(),
        )
```

Use `VacuumInflow` as a semantic alias if desired, but keep `ZeroInflowTrace` as the most precise implementation name.

### 16A.3 Vacuum is decomposable, but realization is method-specific

Vacuum means

```text
gamma_- psi = 0 on Gamma_-.
```

It does **not** mean

```text
psi = 0 on the whole boundary
```

and it does **not** mean

```text
scalar_flux = 0 on the boundary.
```

The physical boundary law is common across methods. The realization is method-specific.

```text
SN vacuum  = ZeroInflowTrace + DiscreteOrdinatesBoundaryRealizer
MoC vacuum = ZeroInflowTrace + CharacteristicBoundaryRealizer
MC vacuum  = ZeroInflowTrace + AbsorbingBoundaryTransition
CP vacuum  = ZeroInflowTrace + ZeroReturn/EscapeProbabilityRealizer
PN vacuum  = ZeroInflowTrace + chosen half-range projection policy
SPN vacuum = ZeroInflowTrace + asymptotic/Marshak/extrapolated realization
```

So do not create unrelated physical concepts called `SNVacuum`, `PNVacuum`, `CPVacuum`, etc. Prefer:

```python
law = ZeroInflowTrace()

sn_bc  = SNBoundaryRealizer().realize(law, sn_space)
pn_bc  = PNBoundaryRealizer(policy=MarshakProjection()).realize(law, pn_space)
moc_bc = CharacteristicBoundaryRealizer().realize(law, moc_space)
cp_bc  = CollisionProbabilityBoundaryRealizer().realize(law, cp_space)
mc_bc  = MonteCarloBoundaryRealizer().realize(law, mc_space)
```

### 16A.4 Method-specific boundary realizers

Use a protocol:

```python
class BoundaryRealizer(Protocol):
    method_name: str

    def realize(
        self,
        law: BoundaryTraceLaw | BoundaryRelation,
        method_space: MethodSpace,
    ) -> MethodBoundaryOperator:
        ...
```

Concrete realizers:

```text
DiscreteOrdinatesBoundaryRealizer
SphericalHarmonicBoundaryRealizer
SPNBoundaryRealizer
CharacteristicBoundaryRealizer
CollisionProbabilityBoundaryRealizer
MonteCarloBoundaryRealizer
StochasticGalerkinBoundaryRealizer
TensorTrainBoundaryRealizer
```

A registry gives extensibility:

```python
class BoundaryRealizerRegistry:
    _items: dict[str, type[BoundaryRealizer]] = {}

    @classmethod
    def register(cls, method_name: str):
        def deco(realizer_cls):
            cls._items[method_name] = realizer_cls
            return realizer_cls
        return deco

    @classmethod
    def for_method(cls, method_name: str) -> BoundaryRealizer:
        return cls._items[method_name]()

@BoundaryRealizerRegistry.register("SN")
@dataclass(frozen=True, slots=True)
class DiscreteOrdinatesBoundaryRealizer:
    method_name: str = "SN"
    def realize(self, law, method_space): ...
```

This lets a new method self-register its boundary interpretation without modifying a central switch statement.

### 16A.5 `SN` realization

`SN` carries explicit incoming and outgoing angular trace degrees of freedom. For a boundary face `f` and ordinate `Omega_i`,

```text
Omega_i · n_f < 0 -> incoming
Omega_i · n_f > 0 -> outgoing.
```

Vacuum is exact:

```text
psi[f, i] = 0 for all incoming (f, i).
```

High-signal classes:

```text
DiscreteOrdinatesBoundaryTrace
IncomingOrdinateMask
OutgoingOrdinateMask
ZeroInflowOrdinateConstraint
OrdinateReflectionPermutation
SNBoundaryOperator
SNBoundaryTensorNetwork
```

Realizer sketch:

```python
@dataclass(frozen=True, slots=True)
class DiscreteOrdinatesBoundaryRealizer:
    def realize(self, law: BoundaryTraceLaw, space: DiscreteOrdinatesPhaseSpace):
        trace = DiscreteOrdinatesBoundaryTrace(space)

        if isinstance(law.response, ZeroBoundaryResponse):
            return ZeroInflowOrdinateConstraint(
                incoming_mask=trace.incoming_mask()
            )

        G = trace.discretize_geometry_map(law.geometry_map)
        R = trace.discretize_response(law.response)
        q = trace.discretize_source(law.source)
        return SNBoundaryOperator(linear=R @ G, source=q)
```

### 16A.6 `MoC` realization

`MoC` boundary traces are track endpoints. Vacuum is exact:

```text
incoming endpoint value = 0.
```

Reflective and periodic boundaries are endpoint gluing maps, possibly creating monodromy.

High-signal classes:

```text
CharacteristicBoundaryTrace
IncomingTrackEndpoint
OutgoingTrackEndpoint
ZeroInflowEndpointConstraint
TrackEndpointGluingTensor
CharacteristicBoundaryGluingNetwork
TrackMonodromy
```

Realizer sketch:

```python
@dataclass(frozen=True, slots=True)
class CharacteristicBoundaryRealizer:
    def realize(self, law: BoundaryTraceLaw, space: CharacteristicTrackSpace):
        trace = CharacteristicBoundaryTrace(space)

        if isinstance(law.response, ZeroBoundaryResponse):
            return ZeroInflowEndpointConstraint(trace.incoming_endpoints())

        G = trace.endpoint_geometry_tensor(law.geometry_map)
        R = trace.endpoint_response_tensor(law.response)
        q = trace.endpoint_source(law.source)
        return CharacteristicBoundaryGluingNetwork(linear=R @ G, source=q)
```

### 16A.7 `MC` realization

For Monte Carlo, vacuum is an absorbing transition into leakage / cemetery state:

```text
boundary hit -> escaped particle / cemetery state.
```

This is the stochastic realization of zero external inflow.

High-signal classes:

```text
BoundaryTransitionKernel
AbsorbingBoundaryTransition
CemeteryState
BoundaryExitEvent
EscapedParticle
LeakageScore
```

Realizer sketch:

```python
@dataclass(frozen=True, slots=True)
class MonteCarloBoundaryRealizer:
    def realize(self, law: BoundaryTraceLaw, space: ParticleTrackingSpace):
        if isinstance(law.response, ZeroBoundaryResponse):
            return AbsorbingBoundaryTransition(
                leakage_score=LeakageScore()
            )
        return BoundaryTransitionKernel.from_trace_law(law, space)
```

### 16A.8 `CP` realization

`CP` usually does not store boundary angular flux degrees of freedom. Boundary effects enter through escape and return probabilities.

Vacuum means escaped particles do not return:

```text
BoundaryReturnKernel = 0.
```

The collision probability balance is

```text
sum_j P_ij + P_i_escape = 1.
```

High-signal classes:

```text
EscapeProbabilityKernel
BoundaryReturnKernel
ZeroReturnBoundaryKernel
SubMarkovBoundaryLoss
CollisionProbabilityBoundaryRealizer
```

Realizer sketch:

```python
@dataclass(frozen=True, slots=True)
class CollisionProbabilityBoundaryRealizer:
    def realize(self, law: BoundaryTraceLaw, space: CollisionProbabilityRegionSpace):
        escape = EscapeProbabilityKernel(space)

        if isinstance(law.response, ZeroBoundaryResponse):
            return ZeroReturnBoundaryKernel(escape=escape)

        return BoundaryReturnKernel.from_trace_law(
            escape=escape,
            law=law,
            region_space=space,
        )
```

### 16A.9 `PN` and `SPN` realization: projected, weak, and nonunique

`PN` is the important exception where vacuum is physically the same but not strongly representable in the truncated space.

The continuous law is half-range:

```text
psi(x, Omega) = 0 for Omega · n < 0.
```

A finite spherical harmonic expansion generally cannot enforce this strongly without overconstraining the moments. Therefore the realization must choose a projection or weak boundary policy.

High-signal policies:

```text
HalfRangeMomentProjection
MarshakProjection
MarkProjection
WeakHalfRangeGalerkin
LeastSquaresHalfRangeProjection
VariationalPNBoundaryForm
PenaltyBoundaryForm
```

Names for realized constraints:

```text
MarshakZeroInflowProjection
MarkZeroInflowProjection
WeakHalfRangeVacuumCondition
PNVacuumMomentCondition
```

Realizer sketch:

```python
@dataclass(frozen=True, slots=True)
class SphericalHarmonicBoundaryRealizer:
    policy: HalfRangeBoundaryProjection

    def realize(self, law: BoundaryTraceLaw, space: SphericalHarmonicMomentSpace):
        if isinstance(law.response, ZeroBoundaryResponse):
            return self.policy.project_zero_inflow(space)
        return self.policy.project_affine_trace_law(law, space)
```

For `SPN` or diffusion-like methods, the realization is even more asymptotic:

```text
ZeroInflowTrace + DiffusionBoundaryRealizer -> ExtrapolatedVacuumBoundary
ZeroInflowTrace + SPNBoundaryRealizer -> SPNMarshakBoundarySystem
```

High-signal names:

```text
ExtrapolatedVacuumBoundary
MarshakDiffusionBoundary
SPNMarshakBoundaryRealizer
AsymptoticVacuumBoundary
```

### 16A.10 Boundary tensor networks by method

A general boundary kernel may have many legs:

```text
face_in, angle_in, group_in, xi_in
face_out, angle_out, group_out, xi_out
```

The tensor-network factorization should separate geometry and response:

```text
outgoing trace
    -> BoundaryGeometryTensor
    -> BoundaryResponseTensor
    -> incoming trace
```

For `SN`:

```text
BoundaryGeometryTensor        -> sparse face/direction permutation or mask
BoundaryResponseTensor        -> angular/energy albedo matrix
IncomingOrdinateMaskTensor    -> zero-inflow enforcement
```

For `PN`:

```text
HarmonicReflectionTensor      -> representation matrix of reflection/rotation
HalfRangeMomentTensor         -> projected inflow condition
MarshakBoundaryTensor         -> weak half-range condition
```

For `MoC`:

```text
TrackEndpointGluingTensor     -> endpoint continuation
SegmentPropagatorTensor       -> affine local solve
CharacteristicBoundaryNetwork -> full track boundary coupling
```

For `CP`:

```text
EscapeProbabilityTensor       -> probability of reaching boundary
BoundaryReturnTensor          -> reentry response
ZeroReturnTensor              -> vacuum escape
```

For `MC`:

```text
BoundaryTransitionKernel      -> stochastic boundary event kernel
AbsorbingBoundaryKernel       -> vacuum leakage
```

For stochastic Galerkin:

```text
BoundaryTripleProductTensor   -> E[Psi_alpha Theta_k Psi_beta]
RandomBoundaryResponseExpansion
StochasticBoundaryOperator
```

For tensor train:

```text
TensorTrainBoundaryOperator
MatrixProductBoundaryOperator
BoundaryMPO
```

### 16A.11 Boundary law self-registration

Physical boundary laws can also self-register:

```python
class BoundaryLawRegistry:
    _items: dict[str, type[BoundaryTraceLaw]] = {}

    @classmethod
    def register(cls, name: str):
        def deco(boundary_cls):
            cls._items[name] = boundary_cls
            boundary_cls.boundary_name = name
            return boundary_cls
        return deco

@BoundaryLawRegistry.register("vacuum")
@dataclass(frozen=True, slots=True)
class VacuumInflow(ZeroInflowTrace):
    pass
```

Do this for semantic user-facing boundary names:

```text
vacuum
prescribed_inflow
reflective
periodic
albedo
symmetry
white_reflection
diffuse_reflection
```

Keep method-specific realizations in a separate registry keyed by method:

```text
BoundaryLawRegistry      -> physical boundary vocabulary
BoundaryRealizerRegistry -> method-specific enforcement
```

This separation prevents `VacuumBoundary` from becoming an overloaded method-specific object.

### 16A.12 Boundary verification laws

Every boundary realization should declare which laws it satisfies exactly and which are approximate.

Universal tests:

```python
boundary.assert_inflow_outflow_classification()
boundary.assert_vacuum_sets_only_inflow_trace()
boundary.assert_outgoing_leakage_unconstrained()
boundary.assert_geometry_map_measure_preserving()
boundary.assert_response_positive_if_declared()
boundary.assert_source_lives_on_incoming_trace()
```

Reflection tests:

```python
reflection.assert_is_involutive()
reflection.assert_preserves_boundary_measure()
reflection.assert_maps_inflow_to_outflow()
reflection.assert_direction_norm_preserved()
```

Periodic tests:

```python
periodic.assert_pairing_bijective()
periodic.assert_normals_opposite()
periodic.assert_measure_preserving()
```

`PN` projection tests:

```python
pn_boundary.assert_projection_policy_declared()
pn_boundary.assert_half_range_moments_consistent()
pn_boundary.assert_not_strong_zero_inflow_unless_supported()
```

`CP` tests:

```python
cp_boundary.assert_escape_probability_bounds()
cp_boundary.assert_zero_return_for_vacuum()
cp_boundary.assert_submarkov_balance()
```

`MC` tests:

```python
mc_boundary.assert_absorbing_for_vacuum()
mc_boundary.assert_leakage_scored_once()
mc_boundary.assert_boundary_kernel_normalized()
```

A high-signal error name should tell the likely bug:

```text
IncomingOutgoingTraceClassificationError
VacuumAppliedToOutgoingTraceError
ScalarFluxDirichletMistakenForVacuumError
BoundaryGeometryMapNotMeasurePreservingError
ReflectionNotInvolutiveError
PeriodicPairingNotBijectiveError
PNBoundaryProjectionPolicyMissingError
CPBoundaryReturnNonzeroForVacuumError
MCVacuumDidNotTerminateHistoryError
```

## 16B. Symmetry sectors, orbifolds, tensor networks, and equivariant cohomology

A sector calculation, such as solving only one octant of a symmetric geometry, should not be treated as merely a smaller geometry with reflective boundaries. The native mathematical structure is a **quotient by a group action**.

Let a finite symmetry group `G_sym` act on the full geometry `X_full`:

```text
rho_x(g): X_full -> X_full.
```

The transport field also transforms in angle, energy, stochastic variables, and possibly internal response channels:

```text
rho(g): psi(x, Omega, g_energy, xi) -> psi(rho_x(g)^-1 x, rho_Omega(g)^-1 Omega, g_energy, rho_xi(g)^-1 xi).
```

The reduced geometry is a fundamental sector or quotient:

```text
X_sector = X_full / G_sym
```

More precisely, when points have stabilizers, the quotient is an **orbifold** or a stratified quotient. Interior points of an octant have orbit size `8`; points on symmetry planes have smaller orbit size; points on symmetry axes and corners have still larger stabilizers.

This matters for:

- volume and surface measures,
- current conservation across quotient faces,
- reflective/symmetry boundary maps,
- angular parity and harmonic transformations,
- Monte Carlo leakage versus symmetry continuation,
- tensor-network factorization of boundary response,
- topological consistency of mesh incidence and traces.

### 16B.1 Fundamental domain versus physical boundary

A sector boundary has two different kinds of patches:

```text
PhysicalBoundaryPatch
    true external boundary: vacuum, inflow, albedo, periodic, etc.

SymmetryBoundaryPatch
    quotient face induced by group action: mirror, rotation, cyclic, dihedral, octant, etc.
```

For an octant in Cartesian coordinates:

```text
X_sector = {x >= 0, y >= 0, z >= 0}
G_sym = {reflections across coordinate planes} ~= (Z2)^3
```

The planes `x=0`, `y=0`, and `z=0` are `SymmetryBoundaryPatch` objects. The outer physical surface may be `VacuumBoundary`, `AlbedoBoundary`, or anything else. Do not encode the octant planes as ordinary vacuum boundaries. Do not even encode them merely as ordinary reflective physical boundaries unless you intentionally want to lose the quotient semantics.

High-signal boundary declarations:

```python
BoundaryDeclaration("outer", kind="vacuum")
SymmetryBoundaryDeclaration("x_min", generator="reflect_x")
SymmetryBoundaryDeclaration("y_min", generator="reflect_y")
SymmetryBoundaryDeclaration("z_min", generator="reflect_z")
```

or:

```python
BoundaryDeclaration("x_min", kind="symmetry", parameters={"generator": "reflect_x"})
```

The high-signal class should be:

```text
SymmetryBoundary
```

not merely `ReflectiveBoundary`, because the former says “this face is a quotient identification under a group action,” whereas the latter says “this is a physical mirror response.”

### 16B.2 Equivariant fields and projectors

Let `rho(g)` be the group representation on the full phase-space field. A symmetry-reduced solve lives in the invariant or equivariant subspace.

For the trivial sector:

```text
rho(g) psi = psi for all g in G_sym.
```

The invariant projector is

```text
Pi_G = (1 / |G_sym|) sum_{g in G_sym} rho(g).
```

For a nontrivial irreducible representation `chi`, the character projector is

```text
Pi_chi = (dim chi / |G_sym|) sum_{g in G_sym} conj(chi(g)) rho(g).
```

Most reactor symmetry reductions use the trivial representation for scalar material/source symmetry. But keeping the representation explicit lets the code handle antisymmetric modes, modal decompositions, perturbations, and debug tests.

High-signal classes:

```text
GroupRepresentation
FieldRepresentation
TrivialRepresentation
CharacterProjector
InvariantProjector
EquivariantProjector
SectorLift
SectorRestriction
```

Useful API:

```python
G_sym = ReflectionGroup.octant()
rho = TransportFieldRepresentation(space=V_full, group=G_sym)

Pi_G = InvariantProjector(rho)
psi_invariant = Pi_G @ psi_full

sector = FundamentalDomainMesh.from_full_mesh(mesh, group=G_sym)
restrict = SectorRestriction(full_space=V_full, sector_space=V_sector, group=G_sym)
lift = SectorLift(sector_space=V_sector, full_space=V_full, group=G_sym)

assert_close(lift @ restrict @ Pi_G, Pi_G)
```

### 16B.3 Boundary operator as quotient trace gluing

On a sector boundary generated by a symmetry `s`, the incoming trace is obtained from the outgoing trace by group action:

```text
gamma_- psi = rho_trace(s) gamma_+ psi.
```

For an octant plane, this becomes specular reflection in angle. Example: on the plane `x=0`, the outward normal of the sector is `n=(-1,0,0)`. Outgoing directions from the sector satisfy `Omega · n > 0`, i.e. `Omega_x < 0`. Incoming directions satisfy `Omega_x > 0`. Reflection maps

```text
Omega = (Omega_x, Omega_y, Omega_z)
    -> (-Omega_x, Omega_y, Omega_z).
```

Thus:

```text
B_x = rho_trace(reflect_x)
```

and not a physical albedo kernel.

Code:

```python
@dataclass(frozen=True, slots=True)
class SymmetryBoundary(BoundaryOperator, key="symmetry"):
    generator: GroupElement
    representation: BoundaryTraceRepresentation

    def as_tensor_network(self, method_space):
        return BoundaryTensorNetwork([
            BoundaryFactor("patch_trace", input_axes=("out_patch",), output_axes=("in_patch",), operator=self.patch_map(method_space)),
            BoundaryFactor("face_pullback", input_axes=("out_face",), output_axes=("in_face",), operator=self.face_pullback(method_space)),
            BoundaryFactor("angle_action", input_axes=("out_angle",), output_axes=("in_angle",), operator=self.angle_action(method_space)),
            BoundaryFactor("energy_identity", input_axes=("from_group",), output_axes=("to_group",), operator=Identity(method_space.energy_space)),
        ])
```

For a purely geometric symmetry boundary, the response kernel is identity. The nontrivial part is the geometry/action factor.

```text
B_symmetry = GeometryActionTensorNetwork
```

For a physical mirror, there may be both geometry and response:

```text
B_mirror = GeometrySpecularMap @ MirrorResponseKernel
```

This distinction is important for hybrid methods. A physical mirror could have roughness or energy-changing reflection. A quotient symmetry face cannot: it must represent exact group identification.

### 16B.4 Tensor-network representation of a sector boundary

A sector boundary tensor network should factor the boundary map into axes:

```text
outflow_trace[patch, face, outgoing_angle, from_group, stochastic_mode]
    -> PatchGeneratorFactor[patch -> group_generator]
    -> QuotientFacePairing[face, generator -> face']
    -> OrientationSignFactor[face, generator]
    -> AngularGroupAction[Omega_out, generator -> Omega_in]
    -> HarmonicRepresentationAction[(ell,m), generator -> (ell,m')]
    -> EnergyIdentity[from_group -> to_group]
    -> StochasticGroupAction[xi_mode, generator -> xi_mode']
    -> inflow_trace[patch, face', incoming_angle, to_group, stochastic_mode']
```

In separable form:

```text
B_sector
= G_patch
  ⊗ G_face
  ⊗ R_angle
  ⊗ I_energy
  ⊗ R_stochastic.
```

For `SN`, `R_angle` is usually a permutation/interpolation of ordinates. For `PN`, `R_angle` is a finite-dimensional representation matrix acting on spherical harmonic coefficients. For stochastic Galerkin, `R_stochastic` may be identity or a representation induced by transformed random variables.

High-signal classes:

```text
PatchGeneratorFactor
QuotientFacePairing
OrbitMultiplicityFactor
OrientationSignFactor
AngularGroupActionFactor
HarmonicGroupRepresentationFactor
StochasticGroupActionFactor
SectorBoundaryTensorNetwork
EquivariantBoundaryTensorNetwork
```

This gives one uniform mechanism for:

- octant symmetry,
- half-core symmetry,
- rotational sectors,
- periodic sectors,
- reflective quotient faces,
- cyclic lattice sectors,
- response-function domain decomposition.

### 16B.5 Orbifold measure and orbit multiplicity

Integrals over the full geometry should be recoverable from the sector by orbit multiplicity:

```text
integral_{X_full} f dx
= integral_{X_sector} orbit_multiplicity(x) f_sector(x) dx_sector.
```

For a finite group action,

```text
orbit_size(x) = |G_sym| / |Stabilizer(x)|.
```

Interior octant cells usually have orbit multiplicity `8`. A cell lying exactly on a symmetry plane is not usually a volume cell in a conforming volume mesh, but boundary faces, edges, and vertices have nontrivial stabilizers and matter for traces, currents, and topology.

Classes:

```text
OrbitMultiplicityField
StabilizerField
IsotropyStratum
StratifiedMeasure
QuotientMeasure
OrbifoldMeasure
```

API:

```python
mu_sector = VolumeMeasure(sector_mesh)
mu_orbifold = OrbifoldMeasure(base=mu_sector, orbit_multiplicity=orbit_size)

reaction_rate_full = mu_orbifold.integrate(sigma_a * phi_sector)
```

For boundary currents, use a boundary quotient measure:

```text
J_full = integral_{Gamma_sector} orbit_multiplicity(face) |Omega · n| psi dGamma.
```

This avoids the common bug where an octant solve computes correct flux but reports one-eighth of the full reaction rate or leakage.

### 16B.6 Equivariant cohomology as a mesh-consistency language

For transport, equivariant cohomology is not needed to solve the streaming equation directly. Its value is architectural: it gives precise names and tests for quotient topology, incidence, orientation, and boundary gluing.

A mesh is a cell complex with chain groups:

```text
C_3 cells
C_2 faces
C_1 edges
C_0 vertices
```

with boundary maps:

```text
partial_k: C_k -> C_{k-1}
```

and cochain coboundary maps:

```text
d^k: C^k -> C^{k+1}.
```

A symmetry group acts on chains:

```text
rho_k(g): C_k -> C_k.
```

A valid equivariant cell complex satisfies:

```text
partial_k rho_k(g) = rho_{k-1}(g) partial_k.
```

Equivalently, on cochains:

```text
d rho^k(g) = rho^{k+1}(g) d.
```

The invariant cochains are:

```text
C^k_G = {alpha in C^k : rho^k(g) alpha = alpha for all g}.
```

Practical classes:

```text
CellComplex
ChainComplex
CochainComplex
GroupActionOnCells
EquivariantChainComplex
InvariantCochainComplex
OrbifoldCohomologyCertificate
```

What this gives you in code:

```python
complex = EquivariantChainComplex(cell_complex=mesh.cell_complex, group=G_sym, action=cell_action)
complex.assert_boundary_commutes_with_group_action()
complex.assert_boundary_squared_zero()
complex.assert_quotient_incidence_consistent()
complex.assert_orientation_signs_consistent()
```

Why this matters for transport:

- face normals must transform consistently under reflections,
- trace spaces must glue the right inflow/outflow DOFs,
- periodic and rotational sectors need correct orientation signs,
- surface currents are cochains on boundary faces,
- leakage/current balance depends on correct quotient boundary orientation,
- MoC track endpoint pairings must respect quotient topology,
- CP escape/return probabilities must not double-count quotient faces.

The high-signal term `EquivariantCohomology` tells Codex and future maintainers: “debug this by checking group action, chain maps, incidence, orientation, invariant cochains, and quotient strata.”

### 16B.7 Method-specific sector spaces

Use method-specific sector names rather than a generic `AugmentedMesh`.

| Method | Sector-aware public class | Meaning |
|---|---|---|
| `SN` | `EquivariantDiscreteOrdinatesPhaseSpace` | cell/group/ordinate unknowns on a sector; symmetry boundary maps ordinates by group action |
| `PN` | `EquivariantSphericalHarmonicMomentSpace` | harmonic moments restricted to an equivariant representation sector; boundary closures use representation matrices |
| `MoC` | `EquivariantCharacteristicTrackSpace` | tracks terminate/continue through quotient symmetry maps; ray bundles are group-compatible |
| `CP` | `OrbifoldCollisionProbabilityRegionSpace` | region probabilities include quotient return/escape kernels and orbit multiplicities |
| `MC` | `OrbifoldParticleTrackingSpace` | particles crossing symmetry faces are lifted/reflected by group action, not killed as leakage |
| Woodcock MC | `EquivariantDeltaTrackingSpace` | majorant and virtual-collision process are invariant/equivariant over the sector |
| stochastic Galerkin | `EquivariantStochasticGalerkinPhaseSpace` | stochastic basis and random-field action commute with geometry symmetry |
| stochastic collocation | `EquivariantCollocationSectorSpace` | each collocation realization preserves or declares broken symmetry |
| tensor train | `EquivariantTensorTrainPhaseSpace` | TT axes include sector/orbit factors and optional representation-sector cores |
| hybrid SPDE/MC | `EquivariantHybridPhaseSpace` | deterministic sector residual and MC correction share the same lift/restrict operators |

### 16B.8 Method-specific boundary resolution for octant symmetry

For `SN`:

```text
SymmetryBoundary(reflect_x)
    -> SNOrdinateSymmetryMap
    -> permutation/interpolation Omega_x < 0 to Omega_x > 0
```

Tests:

```python
boundary.assert_ordinate_map_closes()
boundary.assert_weight_preserving()
boundary.assert_inflow_outflow_are_paired()
```

For `PN`:

```text
SymmetryBoundary(reflect_x)
    -> PNHarmonicReflectionRepresentation
```

Reflection is represented by a matrix on the harmonic moment space. Tests:

```python
R = harmonic_reflection("x")
assert_close(R @ R, I)
assert_close(R.H @ M @ R, M)
```

For `MoC`:

```text
SymmetryBoundary(reflect_x)
    -> TrackEndpointSymmetryContinuation
```

A track leaving the sector through a symmetry face is continued as the reflected incoming track. Tests:

```python
tracks.assert_endpoint_symmetry_pairing()
tracks.assert_track_measure_invariant()
tracks.assert_segment_orbits_consistent()
```

For `CP`:

```text
SymmetryBoundary(reflect_x)
    -> QuotientReturnProbabilityKernel
```

The CP kernel must distinguish physical escape from quotient return. Tests:

```python
P.assert_conservative_with_orbifold_escape()
P.assert_reciprocal_on_quotient_measure()
```

For `MC`:

```text
SymmetryBoundary(reflect_x)
    -> DeterministicGroupActionBoundaryProcess
```

A particle crossing a symmetry face undergoes a deterministic group action. A particle crossing a physical vacuum face leaks. Tests:

```python
boundary_process.assert_symmetry_crossing_preserves_weight()
boundary_process.assert_vacuum_crossing_terminates_history()
```

For `TT`:

```text
SymmetryBoundary(reflect_x)
    -> BoundaryMPOFactor
```

The boundary action can be stored as a matrix product operator if the boundary axes are high-dimensional.

### 16B.9 Equivariance tests and debugging contracts

If a class has `Equivariant` in its name, it must expose:

```python
equivariance_defect(group)
assert_equivariant(group, tol)
```

Core tests:

```python
A.assert_commutes_with_group_action(G_sym)
L.assert_equivariant(G_sym)
C.assert_invariant_coefficients(G_sym)
S.assert_rotational_or_group_equivariant(G_sym)
F.assert_invariant_fission_spectrum(G_sym)
B.assert_trace_equivariant(G_sym)
mu.assert_orbit_multiplicity_consistent(G_sym)
mesh.cell_complex.assert_equivariant_boundary_complex(G_sym)
```

For sector restriction/lifting:

```python
restrict.assert_left_inverse_of_lift_on_sector()
lift.assert_equivariant_extension()
Pi_G.assert_idempotent()
Pi_G.assert_self_adjoint(metric=mass_matrix)
```

For octant symmetry specifically:

```python
sector.assert_has_reflection_generators("reflect_x", "reflect_y", "reflect_z")
sector.assert_symmetry_faces_are_not_physical_leakage()
sector.assert_full_volume_recovered_by_orbit_multiplicity()
sector.assert_current_balance_on_physical_boundary_only()
```

### 16B.10 Elegant sector construction API

Target API:

```python
full_geometry = GeometrySpec(...)
G_sym = ReflectionGroup.octant()
sector_spec = FundamentalDomainSpec.octant()

sector_geometry = full_geometry.quotient_by(
    group=G_sym,
    fundamental_domain=sector_spec,
    representation=TrivialRepresentation(),
)

sector_mesh = MeshBuilder(EqualVolume()).build(sector_geometry)

V_sn = EquivariantDiscreteOrdinatesPhaseSpace.from_mesh(
    mesh=sector_mesh,
    angular_rule=Lebedev(order=17),
    energy_groups=groups,
)

L = StreamingOperator().discretize(V_sn)
C = CollisionOperator(xs.total).discretize(V_sn)
S = ScatteringOperator(xs.scatter).discretize(V_sn)
F = FissionOperator(xs.fission).discretize(V_sn)
B = V_sn.resolve_boundaries()

problem = CriticalityProblem(loss=L + C - S, production=F, boundary=B)
k, psi_sector = PowerIteration().solve(problem)

psi_full = V_sn.sector_lift @ psi_sector
```

The same geometry can be used without quotienting:

```python
full_mesh = MeshBuilder(EqualVolume()).build(full_geometry)
V_full = DiscreteOrdinatesPhaseSpace.from_mesh(full_mesh, angular_rule, groups)
```

A good regression test is:

```python
full_solution = solve_full_if_small(...)
sector_solution = V_sn.sector_lift @ solve_sector(...)
assert_close(Pi_G @ full_solution, sector_solution)
```

---

---

## 17. Energy condensation as Galerkin projection

Energy condensation is not just averaging. It is a projection of operators between energy spaces.

Fine and coarse spaces:

```python
Gf = EnergyGroupSpace(fine_groups)
Gc = EnergyGroupSpace(coarse_groups)
```

Reconstruction and projection:

```python
R = EnergyReconstruction(coarse=Gc, fine=Gf, spectrum=weight_flux)
P = R.H
```

Condensed operator:

```python
A_c = P @ A_f @ R
```

For cross sections:

```python
Sigma_c = P @ Sigma_f @ R
```

Use high-signal names:

```text
EnergyCondensationProjection
FluxWeightedCondensation
ReactionRatePreservingProjection
PetrovGalerkinCondensation
CoarseGroupReconstruction
FineGroupRestriction
```

Tests:

```python
condensation.assert_reaction_rate_preserving(xs, flux)
P.assert_adjoint_of(R)
assert_close(P @ R, Identity(Gc))
```

---

## 18. Spatial homogenization as projection

Spatial homogenization is the spatial analogue of energy condensation.

```text
A_coarse = P A_fine R.
```

Classes:

```text
SpatialRestriction
SpatialProlongation
RegionProjection
FluxWeightedHomogenization
ReactionRatePreservingHomogenization
```

Use the same projection infrastructure as energy condensation.

---

## 19. Scattering data and harmonic projection

Multigroup nuclear scattering data is naturally stored as Legendre moments:

```text
Sigma_{s,ell}^{g' -> g}.
```

Use:

```python
@dataclass(frozen=True, slots=True)
class LegendreMomentScatteringData:
    moments: Array  # [ell, g_to, g_from] or explicit labeled axes
    convention: ScatteringConvention
    groups: EnergyGroupSpace
```

Use labeled axes. Avoid unclear `sigma_s[ell, g, gp]` without naming the direction of transfer.

Preferred axis names:

```text
ell
incident_group
outgoing_group
from_group
to_group
```

Pick one convention and enforce it.

Recommended:

```python
sigma_s_moments[ell, to_group, from_group]
```

#### Symmetry / orbifold reduction

Prefer:

```python
G_sym
rho_g
Pi_G
character_projector
fundamental_domain
sector_mesh
quotient_mesh
orbifold_mesh
sector_space
orbifold_phase_space
sector_lift
sector_restriction
orbit_multiplicity
stabilizer
isotropy_stratum
quotient_measure
physical_boundary_patch
symmetry_boundary_patch
```

Avoid:

```python
reflected_mesh  # ambiguous: physical mirror or quotient?
small_mesh      # low signal
reduced_mesh    # reduced how? energy, space, symmetry, rank?
mirror_bc       # ambiguous: physical response or symmetry identification?
```

because an operator maps from `from_group` to `to_group`.

Scattering operator:

```python
S = ScatteringOperator.from_legendre_moments(xs.scatter)
```

Representations:

```python
S_pn = S.represent(SphericalHarmonicMomentSpace(...))
S_sn = S.represent(DiscreteOrdinatesPhaseSpace(...))
S_tt = S.as_sum_of_products().as_tt(tol=...)
```

---

## 20. Fission as both linear production and branching

Fission has multiple native forms.

### 20.1 Deterministic linear operator

```text
F psi = chi(g) sum_{g'} nu Sigma_f(g') phi(g').
```

This is often low-rank in energy:

```text
F_g = chi ⊗ (nu Sigma_f)^T.
```

Classes:

```text
FissionSpectrum
NuFissionCrossSection
PromptFissionOperator
DelayedFissionOperator
FissionProductionOperator
```

### 20.2 Monte Carlo branching kernel

A fission event produces a random number of descendants with sampled energies and directions.

Classes:

```text
BranchingFissionKernel
MultiplicityDistribution
FissionBank
```

### 20.3 Nonlinear generating function

For event probabilities / extinction / initiation probability, fission enters through generating functions.

Classes:

```text
FissionGeneratingFunction
ProbabilityOfInitiationOperator
ExtinctionProbabilityOperator
```

Design principle:

```python
fission_data = FissionData(nu, sigma_f, chi)
F_det = fission_data.as_linear_operator(space)
F_mc = fission_data.as_branching_kernel(space)
F_gen = fission_data.as_generating_function(space)
```

Same physics data, different representation.

---

## 21. `k`-eigenvalue problems

The criticality equation is

```text
A_loss psi = (1/k) F psi
```

where

```text
A_loss = L + C - S.
```

Equivalently,

```text
A_loss^{-1} F psi = k psi.
```

### 21.1 Class hierarchy

```text
Eigenproblem
├── StandardEigenproblem
├── GeneralizedEigenproblem
├── CriticalityEigenproblem
│   ├── KEigenproblem
│   └── AdjointKEigenproblem
└── AlphaEigenproblem
```

Code:

```python
@dataclass(frozen=True, slots=True)
class CriticalityEigenproblem:
    loss: LinearOperator   # A_loss = L + C - S
    production: LinearOperator  # F
    boundary: BoundaryCondition | None = None

    def multiplication_operator(self):
        return self.loss.inverse() @ self.production
```

### 21.2 Power iteration

Power iteration applies

```text
K = A_loss^{-1} F.
```

Algorithm:

```python
for n in range(maxit):
    source = F @ psi
    psi_next = solve(A_loss, source)
    k_next = normalization(psi_next, psi)
    psi = normalize(psi_next)
```

High-signal classes:

```text
FissionSourceIteration
PowerIteration
WielandtShift
RayleighQuotient
DominanceRatio
```

### 21.3 Wielandt shift

Wielandt shift improves convergence when dominance ratio is high.

Class:

```python
WielandtShiftedKEigenSolver(shift=k_shift)
```

### 21.4 Krylov methods for k-eigenvalues

For multiple eigenpairs or difficult spectra, use matrix-free Krylov on

```text
K = A_loss^{-1} F.
```

Classes:

```text
ArnoldiEigenSolver
KrylovSchurEigenSolver
JacobiDavidsonEigenSolver
ShiftInvertEigenSolver
```

### 21.5 Full eigenspectrum

The full eigenspectrum is method-dependent.

For small dense `PN` or CP systems:

```python
DenseGeneralizedEigenSolver
```

For large sparse/matrix-free systems:

```python
ArnoldiEigenSolver(nev=..., which="largest_magnitude")
ShiftInvertEigenSolver(sigma=...)
KrylovSchurEigenSolver(nev=...)
```

For selected spectral regions:

```python
ContourIntegralEigenSolver(contour=...)
```

For tensor-network representations:

```python
TensorTrainEigenSolver
AlternatingMinimalEnergyEigenSolver
DMRGEigenSolver
```

Use an `EigenSpectrum` object:

```python
@dataclass(frozen=True, slots=True)
class Eigenpair:
    value: complex
    right: Field
    left: Field | None = None
    residual_norm: float | None = None

@dataclass(frozen=True, slots=True)
class EigenSpectrum:
    eigenpairs: tuple[Eigenpair, ...]
    problem: Eigenproblem

    def assert_residuals(self, tol):
        ...

    def sort_by(self, key):
        ...
```

---

## 22. `alpha`-eigenvalue problems and Laplace transform

The alpha eigenvalue problem arises from time dependence.

Start from the homogeneous time-dependent equation:

```text
T dpsi/dt + A_prompt psi = 0
```

where

```text
A_prompt = L + C - S - F.
```

Assume modal form:

```text
psi(t) = exp(alpha t) phi.
```

Then

```text
A_prompt phi = -alpha T phi.
```

So define:

```python
@dataclass(frozen=True, slots=True)
class AlphaEigenproblem:
    time_mass: LinearOperator  # T
    prompt_operator: LinearOperator  # A_prompt = L + C - S - F

    def generalized_form(self):
        # A_prompt phi = -alpha T phi
        return GeneralizedEigenproblem(A=self.prompt_operator, B=self.time_mass, sign=-1)

    def generator(self):
        return -self.time_mass.inverse() @ self.prompt_operator
```

### 22.1 Physical meaning

If `alpha > 0`, the mode grows.

If `alpha < 0`, the mode decays.

If `alpha` has imaginary part, the mode oscillates:

```text
exp(alpha t) = exp(Re(alpha) t) exp(i Im(alpha) t).
```

Use class:

```python
@dataclass(frozen=True, slots=True)
class AlphaMode(Eigenpair):
    @property
    def growth_rate(self):
        return self.value.real

    @property
    def angular_frequency(self):
        return self.value.imag
```

### 22.2 Laplace transform connection

Take the Laplace transform of

```text
T dpsi/dt + A_prompt psi = q(t)
```

with initial condition `psi0`.

Using

```text
Laplace[dpsi/dt] = s psi_hat(s) - psi0,
```

we get

```text
(s T + A_prompt) psi_hat(s) = q_hat(s) + T psi0.
```

Therefore

```text
psi_hat(s) = (s T + A_prompt)^{-1} (q_hat(s) + T psi0).
```

The operator

```text
R(s) = (s T + A_prompt)^{-1}
```

is the Laplace-domain resolvent.

Poles of `R(s)` occur where

```text
s T + A_prompt
```

is singular. Those pole locations are precisely alpha eigenvalues under the sign convention above.

High-signal classes:

```text
LaplaceResolvent
ResolventOperator
ResolventPole
AlphaSpectrum
ModalExpansion
InverseLaplaceSolver
BromwichContour
ResidueExpansion
```

Code:

```python
resolvent = LaplaceResolvent(T=T, A=A_prompt)
psi_hat_s = resolvent(s) @ (q_hat_s + T @ psi0)
alpha_spectrum = resolvent.poles(region=RightHalfPlane())
```

### 22.3 Full alpha spectrum

Methods:

1. generalized eigensolver on `A_prompt phi = -alpha T phi`
2. Arnoldi on generator `G = -T^{-1} A_prompt`
3. shift-invert Arnoldi using `(A_prompt + sigma T)^{-1} T`
4. contour-integral methods using the resolvent `(s T + A_prompt)^{-1}`
5. time-domain methods such as dynamic mode decomposition, mainly as diagnostics or reduced-order analysis
6. tensor-train eigenvalue methods for high-dimensional stochastic systems

Shift-invert relationship:

```text
(A_prompt + sigma T)^{-1} T phi = 1/(sigma - alpha) phi.
```

So eigenvalues near `sigma` become dominant in the transformed operator.

### 22.4 Delayed neutron alpha spectrum

With delayed precursors, the state space becomes a direct sum:

```text
V_total = FluxSpace + PrecursorSpace.
```

The time equation is block-structured:

```text
T_total d/dt [psi, c]^T + A_total [psi, c]^T = 0.
```

Use:

```python
V_total = FluxSpace + PrecursorSpace
alpha_problem = AlphaEigenproblem(T_total, A_total)
```

The alpha spectrum then includes precursor decay effects and kinetic modes.

### 22.5 Modal expansion

If the spectrum is sufficiently well-behaved,

```text
psi(t) = sum_j c_j exp(alpha_j t) phi_j + continuous / unresolved part.
```

For generalized alpha modes, normalize left/right modes by the `T` metric:

```text
<phi_i^dagger, T phi_j> = delta_ij.
```

Class:

```python
@dataclass(frozen=True, slots=True)
class ModalBasis:
    right_modes: tuple[Field, ...]
    left_modes: tuple[Field, ...]
    eigenvalues: tuple[complex, ...]
    metric: LinearOperator  # usually T

    def normalize_biorthogonal(self):
        ...

    def expand(self, initial_state):
        ...
```

---

## 23. Problem classes

Use problem classes as semantic containers.

```text
TransportProblem
├── FixedSourceProblem
├── BoundaryValueProblem
├── InitialValueProblem
├── CriticalityProblem
├── AlphaEigenproblem
├── StochasticGalerkinProblem
├── StochasticCollocationProblem
├── MonteCarloEstimationProblem
├── HybridControlVariateProblem
└── ProbabilityOfEventProblem
```

Example:

```python
@dataclass(frozen=True, slots=True)
class FixedSourceProblem:
    operator: LinearOperator
    source: Field
    boundary: BoundaryCondition

@dataclass(frozen=True, slots=True)
class CriticalityProblem:
    loss: LinearOperator
    production: LinearOperator
    boundary: BoundaryCondition

@dataclass(frozen=True, slots=True)
class InitialValueProblem:
    time_mass: LinearOperator
    generator_operator: LinearOperator
    source: TimeDependentSource
    initial: Field
    boundary: BoundaryCondition
```

---

## 24. Solver hierarchy

```text
Solver
├── LinearSolver
│   ├── DirectSolver
│   ├── KrylovSolver
│   ├── GMRESSolver
│   ├── BiCGStabSolver
│   ├── MultigridSolver
│   └── SweepPreconditionedSolver
├── SourceIterationSolver
├── EigenSolver
│   ├── PowerIteration
│   ├── WielandtShiftedPowerIteration
│   ├── ArnoldiEigenSolver
│   ├── KrylovSchurEigenSolver
│   ├── JacobiDavidsonEigenSolver
│   ├── ContourIntegralEigenSolver
│   └── TensorTrainEigenSolver
├── TimeStepper
│   ├── BackwardEuler
│   ├── CrankNicolson
│   ├── IMEXTimeStepper
│   ├── ExponentialIntegrator
│   └── InverseLaplaceStepper
├── MonteCarloSolver
├── StochasticGalerkinSolver
├── StochasticCollocationSolver
└── HybridSolver
    ├── ControlVariateHybridSolver
    ├── AdjointWeightedHybridSolver
    └── ResidualMonteCarloCorrectionSolver
```

Solvers should self-register:

```python
class Solver(RegistryMixin, ABC):
    registry: ClassVar[dict[str, type["Solver"]]] = {}

class GMRESSolver(LinearSolver, key="gmres"):
    ...

class PowerIteration(EigenSolver, key="power"):
    ...
```

Usage:

```python
solver = Solver.create("gmres", tol=1e-10, maxiter=500)
psi = solver.solve(problem)
```

---

## 25. Hybrid deterministic / SPDE / MC methods

Hybrid methods become clean when deterministic and MC use the same operators, sources, responses, and adjoints.

### 25.1 Control variate form

Let `psi_D` be deterministic approximation.

```text
R(psi) = R(psi_D) + R(psi - psi_D).
```

The correction satisfies

```text
A (psi - psi_D) = q - A psi_D.
```

Code:

```python
residual = q - A @ psi_D
correction = MonteCarloCorrection(A, residual, response=R).estimate(n_histories)
answer = R @ psi_D + correction
```

High-signal classes:

```text
ResidualSource
ControlVariateEstimator
MonteCarloCorrection
DeterministicSurrogate
HybridResidualProblem
```

### 25.2 Adjoint-driven variance reduction

For response `R`, solve

```text
A.H psi_dagger = R.as_source()
```

Then use `psi_dagger` as importance.

```python
importance = solve(A.H, R.as_source())
mc = MonteCarloSolver(importance=AdjointImportance(importance))
```

Classes:

```text
AdjointImportance
WeightWindow
BiasedMarkovKernel
WeightCorrection
ImportanceMap
```

### 25.3 SPDE surrogate plus MC correction

For random inputs:

```python
psi_sg = StochasticGalerkinSolver().solve(A_sg, q_sg)
residual = q_sg - A_sg @ psi_sg
answer = R @ psi_sg + ConditionalMonteCarloCorrection(residual).estimate()
```

Classes:

```text
StochasticSurrogate
ConditionalMonteCarlo
RandomEnvironmentSampler
SPDEControlVariate
```

---

## 26. High-signal naming guide

### 26.1 Suffixes and what they promise

| Suffix | Meaning | Implied tests |
|---|---|---|
| `Space` | valid fields live here; has inner product / dual | dimension, mass matrix, dual consistency |
| `Measure` | integrates and possibly samples | normalization, positivity |
| `Kernel` | integrated transition/interaction object | positivity, normalization, symmetry |
| `Operator` | maps fields to fields | domain/codomain, adjoint consistency |
| `Functional` | maps field to scalar/response | dual-space consistency |
| `Projection` | maps high/rich space to lower/coefficient space | idempotence with reconstruction, conservation |
| `Reconstruction` | lifts coefficients to richer space | adjoint relation to projection |
| `Restriction` | reduces domain or resolution | compatibility with projection |
| `Trace` | boundary restriction | inflow/outflow consistency |
| `Flow` | generated by vector field | semigroup law |
| `Semigroup` | exponential/evolution family | `U(t+s)=U(t)U(s)` |
| `Generator` | infinitesimal evolution operator | resolvent/eigenvalue relation |
| `Process` | stochastic path law | Markov/branching laws |
| `Estimator` | random estimate of functional | bias/variance diagnostics |
| `Certificate` | verified metadata | reproducible residuals |
| `Residual` | defect from equation/law | should be small |
| `Cache` | precomputed implementation detail | invalidation conditions |
| `View` | non-owning representation over another object | source object consistency |

### 26.2 Contract adjectives

Use adjectives when they imply tests.

| Name fragment | Implied mathematical contract |
|---|---|
| `Positive` | maps nonnegative fields to nonnegative fields |
| `Conservative` | preserves total measure / particle balance |
| `SubConservative` | total mass non-increasing |
| `Markov` | positive and normalized transition kernel |
| `SubMarkov` | positive and normalization <= 1 |
| `Invariant` | unchanged under a group action |
| `Equivariant` | commutes appropriately with group action |
| `Symmetric` | kernel/operator symmetric under argument exchange |
| `SelfAdjoint` | `A.H == A` under correct inner product |
| `Reciprocal` | satisfies detailed-balance-like reciprocity |
| `LowRank` | exposes factors and rank |
| `Separable` | factorizes across tensor axes |
| `Sweepable` | admits ordered directional solve |
| `BlockDiagonal` | exposes independent blocks |
| `MatrixFree` | exposes apply but not assembled matrix |
| `Exact` | exact for a stated space or degree |
| `Consistent` | converges / matches weak form under refinement |
| `Sector` | lives on a fundamental domain, not the full geometry |
| `Quotient` | identifies points/traces under an equivalence relation |
| `Orbifold` | quotient with stabilizers/isotropy strata and orbit multiplicities |
| `Stratified` | has dimension-dependent strata such as volume/faces/edges/corners |

If a class is named `ConservativeScatteringKernel`, it must have:

```python
conservation_defect()
assert_conservative(tol)
```

If a class is named `InvariantCubature`, it must have:

```python
invariance_defect(group)
assert_group_invariant(group, tol)
```

### 26.3 Preferred variable names

#### Spaces

```python
X      # spatial space
Omega  # angular space/domain
G      # energy group space
Xi     # stochastic space
V      # full phase space / trial space
W      # test space or codomain space
Gamma  # boundary trace space
Gamma_minus  # inflow trace space
Gamma_plus   # outflow trace space
G_sym        # finite symmetry group; avoid plain G when energy groups are in scope
X_sector     # fundamental-domain spatial space
V_sector     # full phase space restricted to sector/orbifold
Pi_G         # invariant projector over a group action
rho_g        # representation matrix/operator of group element g
```

#### Operators

```python
T          # time mass operator, coefficient of d/dt
L          # streaming operator
C          # collision / total interaction operator
S          # scattering operator
F          # fission operator
B          # boundary operator
A_loss     # L + C - S
A_prompt   # L + C - S - F
G_kinetic  # -T^{-1} A_prompt
K          # multiplication operator A_loss^{-1} F for k-eigenproblem
R          # reconstruction or response; avoid ambiguity in same scope
P          # projection
M          # moment projection or mass matrix; avoid ambiguity in same scope
```

Because `R` can mean response or reconstruction, prefer explicit names when both appear:

```python
reconstruct
response
moment_projector
mass_matrix
```

#### Cross sections

```python
sigma_t
sigma_a
sigma_s
sigma_s_moments
nu_sigma_f
chi
beta
lambda_decay
majorant_cross_section
```

Use `majorant_cross_section`, not `majorand`.

#### Group scattering axes

Prefer:

```python
sigma_s_moments[ell, to_group, from_group]
```

or use labeled-axis objects:

```python
sigma_s_moments.sel(ell=2, to_group=5, from_group=7)
```

#### Monte Carlo

```python
particle
history
ensemble
source_bank
fission_bank
census_bank
hazard_rate
free_flight
collision_site
event_kernel
branching_kernel
score
estimator
```

Use:

- `history` for one path
- `ensemble` for weighted empirical collection
- `bank` for particles queued for a source/generation/census
- `population` for physical branching population

#### MoC

```python
ray_bundle
track_family
track
segment
segment_length
optical_length
tau
attenuation
transverse_weight
track_measure
```

#### CP

```python
region_space
region_basis
collision_probability
escape_probability
response_matrix
reciprocity_residual
conservation_defect
```

#### Boundary

```python
boundary_declaration
boundary_resolution
trace_minus
trace_plus
inflow_trace
outflow_trace
boundary_geometry_map
boundary_response_kernel
specular_reflection_map
periodic_pairing_map
albedo_kernel
```

### 26.4 Names to avoid

Avoid low-signal names:

```text
AugmentedMesh
Data
Info
Manager
Handler
Thing
Stuff
Values
ArrayWrapper
SolverData
MethodData
Temporary
```

Replace with:

```text
DiscreteOrdinatesPhaseSpace
CharacteristicTrackSpace
BoundaryResolution
ScatteringMomentData
SweepGraph
OperatorContext
ExactnessCertificate
ResidualSource
```

`Manager` is often a smell. Ask: does it build, resolve, cache, project, sample, or solve? Name that operation.

---


## 26A. High-signal vocabulary added by local-to-global and boundary decomposition

This section consolidates the vocabulary introduced by the sheaf, local-to-global, and boundary-realization discussion. These names should be preferred because they immediately imply invariants, diagnostics, and likely failure modes.

### 26A.1 Do not overload `Boundary`

Low-signal:

```text
Boundary
BoundaryHandler
BoundaryManager
ApplyBoundary
Vacuum
```

High-signal:

```text
BoundaryTraceLaw
BoundaryRelation
BoundaryGeometryMap
BoundaryResponseKernel
BoundarySource
BoundaryRealizer
MethodBoundaryOperator
BoundaryTensorNetwork
```

When reading these names, Codex should infer:

```text
TraceLaw      -> acts on Gamma_- and Gamma_+
GeometryMap   -> maps trace points or trace indices
ResponseKernel-> maps outgoing/mapped traces into incoming traces
Source        -> lives only on incoming trace space
Realizer      -> method-specific discretization/enforcement
TensorNetwork -> boundary map is factored into explicit local tensors
```

### 26A.2 Vacuum naming

Preferred physical names:

```text
ZeroInflowTrace
VacuumInflow
NoExternalInflow
```

Preferred method-realized names:

```text
ZeroInflowOrdinateConstraint        # SN
ZeroInflowEndpointConstraint        # MoC
AbsorbingBoundaryTransition         # MC
ZeroReturnBoundaryKernel            # CP
MarshakZeroInflowProjection         # PN/SPN/diffusion approximation
WeakHalfRangeVacuumCondition        # PN weak half-range policy
ExtrapolatedVacuumBoundary          # diffusion/SPN asymptotic realization
```

Avoid unless literally true:

```text
ZeroFluxBoundary
DirichletVacuum
BoundaryFluxZero
ScalarFluxVacuum
```

Reason:

```text
vacuum transport BC = gamma_- psi = 0.
```

It does not constrain `gamma_+ psi`, and it does not impose zero scalar flux.

### 26A.3 Local-to-global names by method

| Method | Preferred name | Why it is high signal |
|---|---|---|
| `MoC` | `CharacteristicTransportSheaf` | local segment solutions must restrict and glue through traces |
| `MoC` | `TrackMonodromy` | reflective/periodic cycles impose affine fixed-point constraints |
| `SN` | `SweepDependencyGraph` | solve order is directional and graph-causal |
| `SN` | `UpwindTraceComplex` | face traces are split by sign of `Omega dot n` |
| `PN` | `SphericalHarmonicModule` | angular space is an `SO(3)` representation module |
| `PN` | `IrrepDecomposition` | scattering diagonalizes by irreducible `ell` blocks |
| `CP` | `RegionMomentCosheaf` | region quantities aggregate covariantly |
| `CP` | `SubMarkovCollisionKernel` | escape makes total return probability <= 1 |
| `MC` | `WeightedParticleEnsemble` | implementation-level particle bank with weights |
| `MC` | `EmpiricalPhaseMeasure` | mathematical measure represented by the particle ensemble |
| SPDE | `BochnerField` | random solution as element of `L2(Xi; V)` |
| SPDE | `PolynomialChaosAlgebra` | stochastic basis has multiplication/triple products |
| TN/TT | `TransportFactorGraph` | global object is local tensor contraction |
| Sector | `OrbifoldPhaseSpace` | quotient has stabilizers and multiplicities |
| Sector | `EquivariantBundle` | fields transform by a group representation |
| Time | `TransportSemigroup` | time evolution obeys composition law |
| Alpha | `LaplaceResolvent` | alpha modes are poles of `(s T + A)^-1` |

### 26A.4 Boundary error names

Prefer error names that encode the broken mathematical law.

```text
IncomingOutgoingTraceClassificationError
VacuumAppliedToOutgoingTraceError
ScalarFluxDirichletMistakenForVacuumError
BoundaryGeometryMapNotMeasurePreservingError
BoundarySourceNotOnIncomingTraceError
BoundaryResponseNotPositiveError
ReflectionNotInvolutiveError
ReflectionDidNotMapInflowToOutflowError
PeriodicPairingNotBijectiveError
PeriodicNormalsNotOppositeError
PNBoundaryProjectionPolicyMissingError
PNStrongVacuumRequestedInMomentSpaceError
CPBoundaryReturnNonzeroForVacuumError
CPSubMarkovBalanceViolation
MCVacuumDidNotTerminateHistoryError
MCLeakageScoredMultipleTimesError
TrackEndpointGluingResidualError
TrackMonodromyInconsistentError
SweepGraphCycleNotDeclaredError
EquivariantGluingCocycleError
OrbitMultiplicityMismatchError
```

Each error should point to the diagnostic method likely to expose the bug.

### 26A.5 Diagnostic method names

```python
boundary.assert_vacuum_sets_only_inflow_trace()
boundary.assert_outgoing_leakage_unconstrained()
boundary.assert_geometry_map_measure_preserving()
boundary.assert_response_positive_if_declared()
boundary.assert_source_lives_on_incoming_trace()

trace.assert_inflow_outflow_partition()
trace.assert_no_zero_measure_ambiguity(policy="declared")

reflection.assert_is_involutive()
reflection.assert_maps_inflow_to_outflow()
reflection.assert_direction_norm_preserved()

periodic.assert_pairing_bijective()
periodic.assert_normals_opposite()

sheaf.assert_glued(local_sections)
sheaf.assert_boundary_compatible(local_sections)
sheaf.gluing_residual(local_sections)

track.assert_monodromy_contracting()
track.assert_affine_cycle_consistent()

sweep_graph.assert_topologically_sorted()
sweep_graph.assert_cycles_are_declared()

cp.assert_submarkov()
cp.assert_conservative_with_escape()
cp.assert_zero_return_for_vacuum()

mc_boundary.assert_absorbing_for_vacuum()
mc_boundary.assert_leakage_scored_once()

pn_boundary.assert_projection_policy_declared()
pn_boundary.assert_half_range_moments_consistent()
```

### 26A.6 Variables that should be explicit in boundary code

Use explicit names rather than abbreviations when both geometry and response are present:

```python
Gamma_minus       # incoming trace space
Gamma_plus        # outgoing trace space
gamma_minus       # incoming trace operator
gamma_plus        # outgoing trace operator
geometry_map      # G_boundary
response_kernel   # R_boundary
boundary_source   # q_boundary
boundary_law      # full affine trace law
method_realizer   # method-specific realization
boundary_operator # realized method operator/constraint/kernel
```

When working with tensor networks:

```python
trace_bond
incoming_trace_bond
outgoing_trace_bond
geometry_tensor
response_tensor
source_tensor
gluing_tensor
multiplicity_tensor
sector_lift_tensor
sector_restriction_tensor
```

### 26A.7 Naming rule for method-specific realizations

A method-specific realization should name both the physical law and the method representation.

Good:

```text
ZeroInflowOrdinateConstraint
MarshakZeroInflowProjection
ZeroReturnBoundaryKernel
AbsorbingBoundaryTransition
TrackEndpointGluingTensor
```

Bad:

```text
SNVacuum
PNVacuum
CPVacuum
MCVacuum
```

The bad names hide the fact that these are all realizations of the same `ZeroInflowTrace` law.

## 27. Verification as mathematical law checking

The test suite should be organized around mathematical laws.

### 27.1 Measures

```python
measure.assert_positive()
measure.assert_normalized(expected)
product_measure.assert_fubini(samples=...)
pushforward_measure.assert_change_of_variables(samples=...)
```

### 27.2 Cubature

```python
rule.assert_exactness(degree=p)
rule.assert_spherical_harmonic_exactness(L)
rule.assert_group_invariant(group)
rule.assert_positive_weights()
```

### 27.3 Operators

```python
assert_adjoint_consistency(A, samples=10)
assert_linearity(A, samples=10)
assert_domain_codomain(A)
```

### 27.4 Scattering

```python
S.assert_positive()
S.assert_conservative()
S.assert_block_diagonal_by_ell_if_pn()
S.assert_harmonic_projection_matches_dense_matrix_if_sn()
```

### 27.5 Fission

```python
F.assert_positive()
F.assert_low_rank_energy_structure()
branching_kernel.assert_mean_matches_linear_operator(F)
```

### 27.6 Boundary

```python
B.assert_maps_outflow_to_inflow()
B.assert_patch_support()
ReflectiveBoundary.assert_involutive()
PeriodicBoundary.assert_pairing_bijective()
AlbedoBoundary.assert_submarkov()
```

Reflective map should satisfy:

```text
reflection(reflection(Omega)) = Omega.
```

### 27.7 MoC

```python
tracks.assert_volume_preserving()
segments.assert_positive_lengths()
track_graph.assert_boundary_connected()
propagator.assert_attenuation_bounds()
```

### 27.8 CP

```python
P.assert_positive()
P.assert_conservative()
P.assert_reciprocal()
```

### 27.9 MC

```python
majorant.assert_dominates(sigma_t)
collision_kernel.assert_markov()
branching_kernel.assert_expected_offspring()
estimator.assert_unbiased_on_analytic_case()
```

### 27.10 Eigenproblems

```python
spectrum.assert_residuals(tol)
spectrum.assert_biorthogonal(metric=T)
k_problem.assert_positive_fundamental_mode()
alpha_problem.assert_resolvent_poles_match_modes()
```


### 27.9 Symmetry-sector and orbifold law checks

Any object named `Equivariant`, `Invariant`, `Sector`, `Quotient`, or `Orbifold` must carry explicit law checks.

```python
sector.assert_group_action_defined()
sector.assert_fundamental_domain_covers_full_geometry_without_overlap()
sector.assert_orbit_multiplicity_integrates_to_full_measure()
sector.assert_symmetry_boundaries_are_not_vacuum_boundaries()
sector.assert_physical_boundaries_preserved()
```

Operator equivariance:

```python
for g in G_sym:
    assert_close(A @ rho(g), rho(g) @ A)
```

Boundary trace equivariance:

```python
for g in G_sym:
    assert_close(gamma_minus @ rho(g), rho_minus(g) @ gamma_minus)
    assert_close(gamma_plus @ rho(g), rho_plus(g) @ gamma_plus)
```

Quotient boundary gluing:

```python
B_sym.assert_is_group_action_not_physical_response()
B_sym.assert_involution_for_reflection()
B_sym.assert_weight_preserving()
```

Cohomology / cell-complex consistency:

```python
complex.assert_boundary_squared_zero()
complex.assert_boundary_commutes_with_group_action()
complex.assert_quotient_incidence_consistent()
complex.assert_orientation_signs_consistent()
```

These tests are the fastest way to debug wrong octant leakage, wrong reflection parity, wrong current sign, and wrong multiplicity in reaction-rate tallies.

---

## 28. Suggested source layout

```text
transport/
    core/
        protocols.py
        registry.py
        domains.py
        measures.py
        spaces.py
        fields.py
        operators.py
        kernels.py
        projections.py
        flows.py
        functionals.py

    geometry/
        geometry_spec.py
        regions.py
        materials.py
        boundary_declarations.py
        coordinate_systems.py
        symmetry_declarations.py
        fundamental_domains.py
        quotient_geometry.py
        orbifold_geometry.py

    symmetry/
        groups.py
        group_actions.py
        representations.py
        invariant_projectors.py
        sector_lift.py
        sector_restriction.py
        equivariant_tests.py

    topology/
        chain_complex.py
        cochain_complex.py
        equivariant_chain_complex.py
        orbifold_cohomology.py
        incidence.py
        orientation.py

    mesh/
        mesh_spec.py
        cell_complex.py
        mesh.py
        discretizers.py
        topology.py
        metrics.py

    method_spaces/
        base.py
        sn.py
        pn.py
        moc.py
        cp.py
        mc.py
        stochastic.py
        tensor_train.py
        hybrid.py

    angular/
        sphere.py
        cubature.py
        quadrature.py
        lebedev.py
        product_rules.py
        spherical_harmonics.py
        moment_projection.py
        symmetry_groups.py

    energy/
        groups.py
        condensation.py
        group_kernels.py

    xs/
        materials.py
        cross_sections.py
        scattering.py
        fission.py
        random_xs.py

    physics/
        streaming.py
        collision.py
        scattering.py
        fission.py
        sources.py
        boundary.py
        time.py
        delayed_neutrons.py

    boundaries/
        base.py
        vacuum.py
        inflow.py
        reflective.py
        symmetry.py
        periodic.py
        albedo.py
        tensor_network.py
        trace.py
        orbifold_trace.py

    moc/
        characteristic_flow.py
        ray_bundle.py
        tracks.py
        segments.py
        sweep.py

    cp/
        region_space.py
        green_kernel.py
        collision_probability.py
        response_matrix.py

    mc/
        particle.py
        ensemble.py
        processes.py
        samplers.py
        tallies.py
        estimators.py
        variance_reduction.py
        banks.py

    stochastic/
        probability.py
        random_fields.py
        polynomial_chaos.py
        stochastic_galerkin.py
        collocation.py
        triple_products.py

    tensor_networks/
        tensor_product.py
        sum_of_products.py
        tensor_network.py
        tensor_train.py
        mpo.py
        compression.py
        boundary_networks.py
        equivariant_networks.py

    problems/
        fixed_source.py
        criticality.py
        alpha.py
        initial_value.py
        stochastic.py
        monte_carlo.py
        hybrid.py

    solvers/
        linear.py
        krylov.py
        source_iteration.py
        sweep.py
        eigen.py
        k_eigen.py
        alpha_eigen.py
        time.py
        mc.py
        stochastic.py
        tt.py
        hybrid.py

    verification/
        measures.py
        cubature.py
        operators.py
        boundaries.py
        scattering.py
        fission.py
        moc.py
        cp.py
        mc.py
        eigen.py
        symmetry.py
        topology.py
        manufactured.py
```

---

## 29. Example end-to-end API sketches

### 29.1 Geometry and mesh

```python
geometry = GeometrySpec(
    dimension=2,
    regions=(fuel, moderator, reflector),
    material_assignments=(
        MaterialAssignment("fuel", "uo2"),
        MaterialAssignment("moderator", "h2o"),
    ),
    boundary_declarations=(
        BoundaryDeclaration("outer", "vacuum"),
        BoundaryDeclaration("symmetry_x", "reflective"),
    ),
)

mesh = EqualVolumeMesher(n_cells=10_000).mesh(geometry)
```

### 29.2 SN solve

```python
Omega = LebedevCubature(order=17)
G = EnergyGroupSpace.from_bounds(group_bounds)

V_sn = DiscreteOrdinatesPhaseSpace.from_mesh(
    geometry=geometry,
    mesh=mesh,
    angular_rule=Omega,
    groups=G,
)

T = TimeMassOperator(xs.velocity).on(V_sn)
L = StreamingOperator().on(V_sn)
C = CollisionOperator(xs.total).on(V_sn)
S = ScatteringOperator.from_legendre_moments(xs.scatter).on(V_sn)
F = FissionOperator(xs.fission).on(V_sn)

A_loss = L + C - S
problem = CriticalityProblem(loss=A_loss, production=F, boundary=V_sn.boundary_resolution)

k, psi = PowerIteration(tol=1e-8).solve(problem)
```

### 29.3 PN solve

```python
V_pn = SphericalHarmonicMomentSpace.from_mesh(
    geometry=geometry,
    mesh=mesh,
    order=5,
    groups=G,
)

L = StreamingOperator().galerkin(V_pn)
C = CollisionOperator(xs.total).on(V_pn)
S = ScatteringOperator.from_legendre_moments(xs.scatter).diagonalize(V_pn.harmonics)
F = FissionOperator(xs.fission).on(V_pn)

A = L + C - S - F
psi = GMRESSolver().solve(FixedSourceProblem(A, q, V_pn.boundary_resolution))
```

### 29.4 MoC solve

```python
V_moc = CharacteristicTrackSpace.from_mesh(
    geometry=geometry,
    mesh=mesh,
    angular_rule=Omega,
    track_generator=ModularRayGenerator(spacing=0.05),
)

sweep = CharacteristicSweep(V_moc, xs.total)
S = ScatteringOperator.from_legendre_moments(xs.scatter).on(V_moc)
F = FissionOperator(xs.fission).on(V_moc)

psi = SourceIteration(sweep=sweep, scatter=S, fission=F).solve(q)
```

### 29.5 CP solve

```python
V_cp = CollisionProbabilityRegionSpace.from_mesh(
    geometry=geometry,
    mesh=mesh,
    regions=mesh.material_regions,
)

Gop = CollisionGreenOperator(V_cp, xs.total)
P = CollisionProbabilityKernel.project(Gop, V_cp.region_basis)
P.assert_conservative(include_escape=True)
```

### 29.6 MC solve

```python
V_mc = DeltaTrackingSpace.from_mesh(
    geometry=geometry,
    mesh=mesh,
    majorant_cross_section=MajorantField.from_xs(xs.total),
)

process = NeutronBranchingProcess(
    tracking_space=V_mc,
    scattering=xs.scatter.as_markov_kernel(),
    fission=xs.fission.as_branching_kernel(),
    boundaries=V_mc.boundary_resolution,
)

R = ReactionRateFunctional(xs.absorption, region="fuel")
estimate = MonteCarloSolver(n_histories=1_000_000).estimate(process, source, R)
```

### 29.7 Alpha spectrum from Laplace resolvent

```python
T = TimeMassOperator(xs.velocity).on(V_sn)
A_prompt = L + C - S - F

alpha_problem = AlphaEigenproblem(time_mass=T, prompt_operator=A_prompt)

spectrum = ShiftInvertAlphaSolver(
    shifts=[0.0, -1.0, -10.0],
    nev_per_shift=8,
).solve(alpha_problem)

spectrum.assert_residuals(tol=1e-8)
```

Laplace resolvent:

```python
R_s = LaplaceResolvent(T=T, A=A_prompt)
psi_hat = R_s(s) @ (q_hat(s) + T @ psi0)
```

---

## 30. The highest-level architecture in one expression

The code should make these expressions natural:

```python
V = X * Omega * G * Xi

A_prompt = L + C - S - F

fixed_source = FixedSourceProblem(A_prompt, q, boundary=B)

k_problem = CriticalityProblem(loss=L + C - S, production=F)

alpha_problem = AlphaEigenproblem(T, A_prompt)

A_h = P.H @ A @ P

S = sum(
    SigmaMoment(ell) & AngularMomentOperator(ell) & GroupScatteringMatrix(ell)
    for ell in range(Lmax + 1)
)

mu_sphere = (mu_mu * mu_phi).pushforward(SphericalMap())

hybrid_answer = R @ psi_surrogate + MonteCarloCorrection(A, q - A @ psi_surrogate, R)
```

If those expressions are natural in the code, the architecture is aligned with the mathematics.

---

## 31. Final design rules

1. **Use `GeometrySpec` for semantic continuous geometry.** It declares regions, materials, dimensions, and unresolved boundary declarations.

2. **Use `SpatialMesh` for the discretized cell complex.** It owns topology, geometry, volumes, areas, normals, adjacency, and material tags.

3. **Use `MethodSpace`, not `AugmentedMesh`, for method-specific spaces.** The method object defines where unknowns live and how boundaries/operators are represented.

4. **Use specific method-space names.** Prefer `DiscreteOrdinatesPhaseSpace`, `SphericalHarmonicMomentSpace`, `CharacteristicTrackSpace`, `CollisionProbabilityRegionSpace`, `ParticleTrackingSpace`, `StochasticGalerkinPhaseSpace`, and `TensorTrainPhaseSpace`.

5. **Make measures first-class.** Quadrature, cubature, stochastic collocation, boundary current integration, and MC empirical particles are all measures.

6. **Make projections first-class.** Galerkin discretization, energy condensation, homogenization, moment projection, tallies, and stochastic Galerkin are all projections.

7. **Make boundary conditions operators on trace spaces.** Geometry declarations resolve into method-specific boundary operators only after method space creation.

8. **Separate boundary geometry from boundary response.** Use `BoundaryGeometryMap`, `BoundaryResponseKernel`, and `BoundaryTensorNetwork`.

9. **Exploit representation theory for `PN`.** Rotationally invariant scattering diagonalizes by spherical harmonic irreps.

10. **Exploit harmonic projection for `SN`.** Do not default to dense angular scattering matrices.

11. **Represent multigroup scattering as a sum of tensor products.** Tensor train is an optional compression backend, not the first representation.

12. **Treat MoC as characteristic flow.** Use `CharacteristicTrackSpace`, `RayBundle`, `SegmentPropagator`, and `TrackMeasure`.

13. **Treat CP as projected Green kernels.** Use `CollisionProbabilityKernel`, `EscapeProbability`, and reciprocity/conservation tests.

14. **Treat MC as Markov/branching processes.** Use `ParticleHistory`, `ParticleEnsemble`, `MarkovKernel`, `BranchingKernel`, and `Estimator`.

15. **Use `T`, `L`, `S`, and `F` consistently.** Let `C` be collision/removal. Let `A_loss = L + C - S`. Let `A_prompt = L + C - S - F`.

16. **Represent `k` as a criticality eigenproblem.** `A_loss psi = (1/k) F psi` or `A_loss^{-1} F psi = k psi`.

17. **Represent `alpha` as a time-generator eigenproblem.** `A_prompt phi = -alpha T phi` and `R(s) = (s T + A_prompt)^{-1}`.

18. **Use Laplace resolvents for alpha analysis and transient modal structure.** Alpha eigenvalues are poles of the Laplace-domain resolvent.

19. **Use high-signal names.** A name should imply invariants, tests, and mathematical laws.

20. **Make laws executable.** Every object named `Conservative`, `Positive`, `Invariant`, `Markov`, `Reciprocal`, or `Exact` must expose a residual and an assertion method.

The architectural north star is:

> **Neutron transport is an algebra of positive kernels, trace maps, flows, projections, tensor products, stochastic processes, and spectral problems. Build those primitives first; every method becomes a faithful representation of the same mathematical object.**


### 31.9 Treat symmetry reduction as quotient geometry, not as ad hoc boundary flags

For sector calculations, encode the group action explicitly. A sector face generated by octant symmetry is a `SymmetryBoundaryPatch`, not a physical vacuum wall. Its boundary operator is a group-action trace gluing map. Its measure is a quotient/orbifold measure with orbit multiplicity. Its tests are equivariance, chain-complex consistency, orientation consistency, and full-measure recovery.

Use:

```python
sector_geometry = geometry.quotient_by(group=G_sym, fundamental_domain=sector)
V_sector = EquivariantDiscreteOrdinatesPhaseSpace.from_mesh(sector_mesh, ...)
```

not:

```python
geometry.add_boundary("x_min", kind="reflective")  # loses quotient semantics
```

unless the boundary is genuinely a physical mirror.

---

## 32. Detailed primitive specifications

This section gives more concrete class contracts. These are not final implementations; they are semantic skeletons for Codex to follow.

### 32.1 Labeled axes

Avoid bare array shapes as much as possible. Transport arrays have axes with meaning.

```python
@dataclass(frozen=True, slots=True)
class Axis:
    name: str
    size: int
    coordinate: object | None = None
    measure: Measure | None = None

@dataclass(frozen=True, slots=True)
class AxisProduct:
    axes: tuple[Axis, ...]

    def index(self, name: str) -> int:
        ...

    def without(self, name: str) -> "AxisProduct":
        ...

    def append(self, axis: Axis) -> "AxisProduct":
        ...
```

Preferred axis names:

```text
cell
face
region
angle
ordinate
ell
m
moment
group
to_group
from_group
time
stochastic_mode
particle
history
segment
track
patch
trace_dof
```

A field should expose axes:

```python
psi.axes.names
# ("cell", "angle", "group")
```

Instead of:

```python
psi.shape
# (10000, 80, 172)
```

which has low semantic signal.

### 32.2 Maps, pullbacks, and pushforwards

Maps are native mathematical objects.

```python
class Map(ABC):
    domain: Domain
    codomain: Domain

    @abstractmethod
    def __call__(self, x):
        ...

    def pullback(self, f):
        return Pullback(self, f)

    def pushforward(self, measure: Measure):
        return PushforwardMeasure(self, measure)
```

Use maps for:

- geometry mappings
- periodic boundary identifications
- specular reflections
- coordinate transforms
- spherical parameterization
- stochastic random-field maps
- material lookup maps

High-signal map names:

```text
SphericalMap
SpecularReflectionMap
PeriodicPairingMap
MaterialLookupMap
CellLocatorMap
TrackToCellMap
EnergyCollapseMap
RandomFieldRealizationMap
```

### 32.3 Measures

```python
class Measure(ABC):
    domain: Domain

    @abstractmethod
    def integrate(self, f):
        ...

    def tensor(self, other: "Measure") -> "ProductMeasure":
        return ProductMeasure((self, other))

    def pushforward(self, map: Map) -> "PushforwardMeasure":
        return PushforwardMeasure(map=map, base=self)

    def __mul__(self, other):
        return self.tensor(other)
```

Concrete measures:

```python
@dataclass(frozen=True, slots=True)
class ProductMeasure(Measure):
    factors: tuple[Measure, ...]

@dataclass(frozen=True, slots=True)
class PushforwardMeasure(Measure):
    map: Map
    base: Measure

@dataclass(frozen=True, slots=True)
class EmpiricalMeasure(Measure):
    particles: tuple
    weights: Array
```

An MC ensemble can subclass or contain `EmpiricalMeasure`.

### 32.4 Spaces and inner products

A space must define how fields are compared.

```python
class Space(ABC):
    name: str
    axes: AxisProduct

    @abstractmethod
    def zero(self) -> "Field":
        ...

    @abstractmethod
    def inner(self, u: "Field", v: "Field"):
        ...

    def dual(self) -> "DualSpace":
        return DualSpace(self)

    def tensor(self, other: "Space") -> "TensorProductSpace":
        return TensorProductSpace((self, other))

    def direct_sum(self, other: "Space") -> "DirectSumSpace":
        return DirectSumSpace((self, other))

    def __mul__(self, other):
        return self.tensor(other)

    def __add__(self, other):
        return self.direct_sum(other)
```

For quadrature spaces, the inner product must include weights.

```python
class DiscreteAngularSpace(Space):
    cubature: CubatureRule

    def inner(self, u, v):
        return sum(w_i * u_i.conjugate() * v_i for i, w_i in enumerate(self.cubature.weights))
```

### 32.5 Fields

```python
@dataclass(slots=True)
class Field:
    space: Space
    data: object
    representation: "Representation" | None = None

    def inner(self, other: "Field"):
        return self.space.inner(self, other)

    def norm(self):
        return self.inner(self) ** 0.5

    def project(self, projection: "ProjectionOperator"):
        return projection @ self

    def as_representation(self, representation):
        ...
```

Fields should not be allowed to silently combine if their spaces are incompatible.

```python
psi_sn + psi_pn
# should fail unless an explicit projection/reconstruction exists
```

### 32.6 Operators

```python
class LinearOperator(ABC):
    domain: Space
    codomain: Space

    @abstractmethod
    def apply(self, x: Field) -> Field:
        ...

    def compose(self, other: "LinearOperator") -> "LinearOperator":
        assert other.codomain == self.domain
        return ComposedOperator(self, other)

    def tensor(self, other: "LinearOperator") -> "LinearOperator":
        return TensorProductOperator((self, other))

    @property
    @abstractmethod
    def H(self) -> "LinearOperator":
        ...

    def __matmul__(self, other):
        if isinstance(other, LinearOperator):
            return self.compose(other)
        return self.apply(other)

    def __and__(self, other):
        return self.tensor(other)

    def __add__(self, other):
        return SumOperator((self, other))

    def __sub__(self, other):
        return SumOperator((self, -1 * other))
```

### 32.7 Projections and reconstructions

Use projections as operators between spaces.

```python
class ProjectionOperator(LinearOperator):
    rich_space: Space
    coefficient_space: Space

    @property
    def domain(self):
        return self.rich_space

    @property
    def codomain(self):
        return self.coefficient_space

class ReconstructionOperator(LinearOperator):
    coefficient_space: Space
    rich_space: Space
```

For Galerkin methods:

```text
A_h = P.H @ A @ P
```

where `P` may be a reconstruction from discrete coefficients to continuous/semi-continuous fields. Depending on naming preference, use:

```python
reconstruct = TrialReconstruction(...)
project = reconstruct.H
A_h = project @ A @ reconstruct
```

This avoids ambiguity.

### 32.8 Representations

Separate mathematical object from storage representation.

```python
class Representation(ABC):
    @abstractmethod
    def materialize(self):
        ...

class DenseRepresentation(Representation): ...
class SparseRepresentation(Representation): ...
class MatrixFreeRepresentation(Representation): ...
class SweepRepresentation(Representation): ...
class SumOfProductsRepresentation(Representation): ...
class TensorTrainRepresentation(Representation): ...
class MonteCarloRepresentation(Representation): ...
```

An operator may choose representation lazily:

```python
A.represent("matrix_free")
A.represent("sparse")
A.represent("tt", tol=1e-8)
A.represent("sweep")
```

---

## 33. MethodSpace contracts in detail

### 33.1 `DiscreteOrdinatesPhaseSpace`

Primary unknown:

```text
psi[cell, ordinate, group]
```

Optional axes:

```text
time
stochastic_mode
```

Must contain:

```text
mesh
angular_rule
energy_groups
boundary_trace_spaces
inflow_trace_dofs
outflow_trace_dofs
sweep_graphs
upwind_stencils
downwind_stencils
angular_mass_matrix
harmonic_projection_cache
```

Precomputations:

```text
Omega_x, Omega_y, Omega_z diagonal factors
face_dot_products Omega · n
inflow/outflow masks by face and ordinate
sweep order per direction/octant
moment matrices Y_lm(Omega_i)
weighted moment projector Y^* W
```

Boundary resolution:

```text
VacuumBoundary -> zero inflow trace / no incoming current
ReflectiveBoundary -> ordinate permutation/interpolation
PeriodicBoundary -> face pairing + ordinate identity
AlbedoBoundary -> boundary tensor network over patch/angle/group
```

Recommended cached objects:

```python
@dataclass(slots=True)
class DiscreteOrdinatesOperatorContext:
    sweep_graphs: dict
    direction_cosines: DirectionCosineFactors
    trace_resolution: BoundaryResolution
    harmonic_projector: HarmonicMomentProjection
    harmonic_reconstructor: HarmonicReconstruction
```

### 33.2 `SphericalHarmonicMomentSpace`

Primary unknown:

```text
psi[cell, moment, group]
```

where `moment` indexes `(ell, m)`.

Must contain:

```text
harmonic_index
ell_blocks
moment_mass_matrix
streaming_coupling_matrices
boundary_closure
rotation_representation
```

Precomputations:

```text
coupling coefficients for Omega_x, Omega_y, Omega_z acting on harmonics
block diagonal scattering by ell
parity operators
half-range moment matrices if needed for boundary closures
```

Boundary resolution:

```text
VacuumBoundary -> MomentBoundaryClosure, often half-range or Marshak-like
ReflectiveBoundary -> parity/reflection operator on moments
PeriodicBoundary -> spatial pairing with moment identity
AlbedoBoundary -> boundary moment kernel
```

High-signal cache:

```python
@dataclass(slots=True)
class MomentOperatorContext:
    harmonic_index: HarmonicIndex
    irrep_blocks: tuple[IrrepBlock, ...]
    streaming_couplings: StreamingMomentCouplings
    scattering_diagonalizer: RotationallyInvariantDiagonalization
    boundary_closure: MomentBoundaryClosure
```

### 33.3 `CharacteristicTrackSpace`

Primary unknowns may include:

```text
angular_flux[track_endpoint, group]
segment_flux[segment, group]
cell_scalar_flux[cell, group]
```

Must contain:

```text
ray_bundles
track_families
tracks
segments
segment_to_cell map
track endpoint connectivity
track measure
transverse weights
boundary continuation maps
```

Precomputations:

```text
segment lengths
optical path lengths tau = Sigma_t * length
attenuation factors exp(-tau)
source integration factors (1 - exp(-tau)) / Sigma_t
cell-volume preservation weights
```

Boundary resolution:

```text
VacuumBoundary -> zero incoming endpoint source and terminate outgoing leakage
ReflectiveBoundary -> endpoint-to-endpoint continuation via specular map
PeriodicBoundary -> endpoint pairing
AlbedoBoundary -> endpoint response kernel
```

High-signal context:

```python
@dataclass(slots=True)
class CharacteristicOperatorContext:
    ray_bundles: tuple[RayBundle, ...]
    segment_index: SegmentIndex
    segment_propagators: SegmentPropagatorTable
    track_coupling_graph: TrackCouplingGraph
    track_measure: TrackMeasure
```

### 33.4 `CollisionProbabilityRegionSpace`

Primary unknown:

```text
region_source[region, group]
region_flux[region, group]
```

Must contain:

```text
region basis
region volumes
region material averages
surface escape regions
optical distances / geometric kernels if precomputed
```

Precomputations:

```text
collision probability matrix P_ij
escape probability P_i_escape
transmission probabilities
reciprocity weights V_i Sigma_i
```

Boundary resolution:

```text
VacuumBoundary -> escape loss with no incoming return current
ReflectiveBoundary -> modified escape/return kernel
AlbedoBoundary -> boundary return probability kernel
```

Context:

```python
@dataclass(slots=True)
class CollisionProbabilityOperatorContext:
    region_space: RegionSpace
    collision_probability: CollisionProbabilityKernel
    escape_probability: EscapeProbability
    reciprocity_weights: Field
```

### 33.5 `ParticleTrackingSpace`

Primary objects:

```text
ParticleState
ParticleHistory
ParticleEnsemble
```

Must contain:

```text
tracking geometry
cell/material locator
surface intersection accelerator
boundary samplers
cross-section samplers
collision kernels
fission branching kernels
```

Precomputations:

```text
spatial acceleration structure
material lookup acceleration
surface crossing tables if possible
majorant cross section for delta tracking
reaction channel cumulative distributions
```

Boundary resolution:

```text
VacuumBoundary -> terminate outgoing history at leakage; never impose interior zero flux
ReflectiveBoundary -> deterministic direction update
PeriodicBoundary -> deterministic position remap
AlbedoBoundary -> sample boundary response + weight correction
```

Context:

```python
@dataclass(slots=True)
class TrackingOperatorContext:
    material_locator: MaterialLocator
    surface_tracker: SurfaceTracker
    boundary_process: BoundaryProcess
    collision_process: CollisionProcess
    fission_process: BranchingProcess
```

### 33.6 `DeltaTrackingSpace`

Adds:

```text
majorant_cross_section
virtual_collision_probability
real_collision_acceptance
```

Tests:

```python
majorant.assert_dominates(xs.total)
virtual_collision_kernel.assert_markov()
```

### 33.7 `StochasticGalerkinPhaseSpace`

Primary unknown:

```text
psi[cell, angle_or_moment, group, stochastic_mode]
```

Must contain:

```text
probability measure
stochastic basis
multi-index set
triple product tensor
random coefficient expansions
```

Precomputations:

```text
E[Psi_a Psi_b]
E[Psi_a theta_k Psi_b]
stochastic mass matrix
stochastic operator factors
```

Context:

```python
@dataclass(slots=True)
class StochasticGalerkinOperatorContext:
    basis: PolynomialChaosBasis
    multi_index: MultiIndexSet
    triple_products: TripleProductTensor
    random_coefficients: RandomCoefficientExpansion
```

### 33.8 `StochasticCollocationSpace`

Primary unknown:

```text
psi[collocation_node, cell, angle_or_moment, group]
```

Must contain:

```text
random-space cubature
collocation cases
case weights
case-specific material realizations
```

Use `CollocationCaseEnsemble` if emphasizing many deterministic cases.

### 33.9 `TensorTrainPhaseSpace`

Primary unknown:

```text
TT cores over axes [cell, angle/moment, group, xi_1, ..., xi_d]
```

Must contain:

```text
axis order
rank policy
rounding tolerance
orthogonality center
operator MPO/TT factors
```

Names:

```text
TensorTrainField
TensorTrainOperator
MatrixProductOperator
TTCore
TTRankProfile
TTOrthogonalityCenter
TTRoundingPolicy
```

### 33.10 `ControlVariateHybridSpace`

Primary objects:

```text
deterministic_surrogate
residual_source
mc_correction_process
response_functional
```

Must contain:

```text
projection between deterministic and MC tallies
surrogate evaluation map
residual sampler
control variate estimator
```

### 33.11 `EquivariantSectorSpace` and `OrbifoldPhaseSpace`

Primary unknown:

```text
psi[sector_cell, angle_or_moment, group, optional_stochastic_mode]
```

with implicit or explicit lifting to the full domain.

Must contain:

```text
symmetry_group
group_action_on_cells
group_action_on_angles
group_action_on_trace_space
fundamental_domain
physical_boundary_patches
symmetry_boundary_patches
sector_lift
sector_restriction
invariant_projector
orbit_multiplicity_field
stabilizer_field
isotropy_strata
quotient_measure
equivariant_chain_complex
```

Precomputations:

```text
cell orbit representatives
face orbit representatives
orientation signs under group action
trace generator maps
sector-to-full lifting maps
full-to-sector restriction maps
orbifold volume and boundary multiplicities
method-specific symmetry boundary tensor networks
```

Boundary resolution:

```text
Physical VacuumBoundary -> zero inflow trace / no incoming current
SymmetryBoundary -> group-action trace gluing
Physical ReflectiveBoundary -> mirror response kernel
PeriodicBoundary -> quotient or physical periodic pairing, depending on declaration
AlbedoBoundary -> physical response kernel; not allowed on pure quotient face unless explicitly modeling physical albedo
```

Recommended context:

```python
@dataclass(slots=True)
class EquivariantOperatorContext:
    symmetry_group: FiniteGroup
    representation: GroupRepresentation
    invariant_projector: InvariantProjector
    sector_lift: SectorLift
    sector_restriction: SectorRestriction
    orbit_multiplicity: OrbitMultiplicityField
    equivariant_complex: EquivariantChainComplex
    boundary_networks: BoundaryResolution
```

Tests:

```python
context.invariant_projector.assert_idempotent()
context.sector_lift.assert_equivariant_extension()
context.equivariant_complex.assert_boundary_commutes_with_group_action()
context.orbit_multiplicity.assert_integrates_to_full_measure()
```

---

## 34. Boundary registry and resolution workflow

Boundary behavior should flow through three stages.

### 34.1 Stage 1: declaration in geometry

```python
BoundaryDeclaration(patch="outer", kind="vacuum")
BoundaryDeclaration(patch="symmetry", kind="reflective")
BoundaryDeclaration(patch="left", kind="periodic", parameters={"pair": "right"})
BoundaryDeclaration(patch="shield", kind="albedo", parameters={"response": "diffuse"})
```

No method-specific operator exists yet.

### 34.2 Stage 2: patch resolution in mesh

The mesh identifies which faces belong to which patch.

```python
mesh.boundary_patches["outer"] -> tuple[FaceId, ...]
mesh.boundary_patches["symmetry"] -> tuple[FaceId, ...]
```

### 34.3 Stage 3: method-specific boundary operator resolution

```python
class BoundaryResolver:
    def resolve(self, declaration, method_space):
        boundary_cls = BoundaryOperator.registry[declaration.kind]
        return boundary_cls.from_declaration(declaration, method_space)
```

Each boundary class may dispatch by method space:

```python
class ReflectiveBoundary(BoundaryOperator, key="reflective"):
    @classmethod
    def from_declaration(cls, declaration, method_space):
        if isinstance(method_space, DiscreteOrdinatesPhaseSpace):
            return SNReflectiveBoundary(...)
        if isinstance(method_space, SphericalHarmonicMomentSpace):
            return PNReflectiveBoundary(...)
        if isinstance(method_space, CharacteristicTrackSpace):
            return MoCReflectiveBoundary(...)
        if isinstance(method_space, ParticleTrackingSpace):
            return MCReflectiveBoundary(...)
        raise TypeError(...)
```

Alternatively, use double registry:

```python
BoundaryImplementationRegistry[(boundary_kind, method_kind)] = concrete_class
```

High-signal concrete names:

```text
SNReflectiveBoundary
PNReflectiveMomentBoundary
MoCReflectiveTrackBoundary
MCSpecularReflectionBoundary
SNSymmetryBoundary
PNSymmetryMomentBoundary
MoCSymmetryContinuationBoundary
MCSymmetryCrossingBoundary
OrbifoldTraceBoundary
EquivariantBoundaryTensorNetwork
```

---

## 35. Boundary tensor network skeleton

A boundary tensor network separates axes and factors.

```python
@dataclass(frozen=True, slots=True)
class BoundaryAxis:
    name: str
    space: Space

@dataclass(frozen=True, slots=True)
class BoundaryFactor:
    name: str
    input_axes: tuple[str, ...]
    output_axes: tuple[str, ...]
    operator: LinearOperator

@dataclass(slots=True)
class BoundaryTensorNetwork(BoundaryOperator):
    factors: tuple[BoundaryFactor, ...]
    contraction_plan: ContractionPlan | None = None

    def apply(self, outgoing_trace):
        return contract_boundary_network(self.factors, outgoing_trace, self.contraction_plan)
```

Example albedo factorization:

```text
outflow_trace[patch, face, outgoing_angle, from_group]
    -> PatchIncidenceFactor
    -> BoundaryGeometryMap
    -> AngularAlbedoKernel[outgoing_angle -> incoming_angle]
    -> EnergyAlbedoKernel[from_group -> to_group]
    -> inflow_trace[patch, face, incoming_angle, to_group]
```

This corresponds to

```text
B = G_patch_face ⊗ K_angle ⊗ K_energy
```

when separable.

If not separable, use low-rank/tensor-network response:

```text
K_boundary(face, Omega_in, g_in ; face', Omega_out, g_out)
≈ sum_r X_r(face, face') A_r(Omega_in, Omega_out) E_r(g_in, g_out)
```

Class:

```python
class LowRankBoundaryResponse(BoundaryResponseKernel):
    factors: tuple[BoundaryResponseFactor, ...]
```

This is especially useful for:

- albedo libraries
- reflector response matrices
- hybrid deterministic/MC boundary coupling
- domain decomposition
- response-function acceleration


### 35.1 Sector boundary tensor networks

For a symmetry sector, the boundary tensor network should express quotient gluing as a product of explicit factors:

```text
B_sector = PatchGeneratorFactor
         & QuotientFacePairing
         & OrientationSignFactor
         & AngularGroupActionFactor
         & EnergyIdentityFactor
         & StochasticGroupActionFactor
```

This is not merely elegant notation. It lets the code separate:

- which patch generated the symmetry crossing,
- how the face maps to its quotient partner,
- whether orientation signs change,
- how angle or harmonic moments transform,
- whether energy is untouched,
- whether random variables transform,
- whether the response is physical or purely geometric.

A physical albedo boundary and a quotient symmetry boundary can therefore share tensor-network machinery but differ in semantics:

```text
B_albedo  = PhysicalResponseKernel @ GeometryMap
B_sector  = GeometryGroupAction @ IdentityResponse
```

Debug rule: if an octant symmetry face produces leakage, the wrong object was probably used. It should be a `SymmetryBoundary`, not a `VacuumBoundary`.

---

## 36. Eigen-solver implementation details

### 36.1 `k` power iteration pseudocode

```python
class PowerIteration(EigenSolver, key="power"):
    def solve(self, problem: CriticalityProblem):
        A = problem.loss
        F = problem.production
        psi = self.initial_guess(problem)
        k = 1.0

        for it in range(self.maxiter):
            fission_source = F @ psi
            psi_next = self.inner_solver.solve(A, fission_source)

            k_next = self.estimate_k(psi_next, psi, F, A)
            psi_next = self.normalize(psi_next)

            residual = A @ psi_next - (1.0 / k_next) * (F @ psi_next)
            if residual.norm() <= self.tol:
                return Eigenpair(value=k_next, right=psi_next, residual_norm=residual.norm())

            psi, k = psi_next, k_next

        raise ConvergenceError(...)
```

Good `k` estimates:

```text
source ratio
Rayleigh quotient
response-weighted production/loss ratio
```

Rayleigh-like quotient:

```text
1/k = <psi_dagger, A_loss psi> / <psi_dagger, F psi>
```

When adjoint mode is unavailable, use a positive functional such as total fission source.

### 36.2 Shift-invert alpha pseudocode

Alpha problem:

```text
A phi = -alpha T phi.
```

Shift transform near `sigma`:

```text
M_sigma = (A + sigma T)^{-1} T.
```

If `alpha` is an eigenvalue, then `M_sigma` has eigenvalue

```text
mu = 1 / (sigma - alpha).
```

So recover:

```text
alpha = sigma - 1 / mu.
```

Pseudocode:

```python
class ShiftInvertAlphaSolver(EigenSolver, key="shift_invert_alpha"):
    def solve(self, problem: AlphaEigenproblem):
        A = problem.prompt_operator
        T = problem.time_mass
        sigma = self.shift

        shifted = A + sigma * T

        def apply_M(x):
            return self.inner_solver.solve(shifted, T @ x)

        mu_spectrum = ArnoldiEigenSolver(nev=self.nev).solve(MatrixFreeOperator(apply_M))

        alpha_pairs = []
        for pair in mu_spectrum.eigenpairs:
            alpha = sigma - 1.0 / pair.value
            residual = A @ pair.right + alpha * (T @ pair.right)
            alpha_pairs.append(Eigenpair(alpha, pair.right, residual_norm=residual.norm()))

        return EigenSpectrum(tuple(alpha_pairs), problem=problem)
```

### 36.3 Resolvent-based spectrum

Define:

```python
class LaplaceResolvent:
    def __init__(self, T, A):
        self.T = T
        self.A = A

    def operator(self, s):
        return (s * self.T + self.A).inverse()

    def apply(self, s, rhs):
        return self.linear_solver.solve(s * self.T + self.A, rhs)
```

Contour methods use integrals of the resolvent:

```text
P_C = (1 / 2πi) ∮_C (sT + A)^{-1} T ds
```

This projects onto the invariant subspace whose alpha eigenvalues lie inside contour `C`.

High-signal classes:

```text
SpectralProjector
ResolventContour
BromwichContour
RightHalfPlaneContour
ContourIntegralEigenSolver
```

### 36.4 Modal transient reconstruction

Given alpha modes:

```text
psi(t) ≈ sum_j c_j exp(alpha_j t) phi_j.
```

Coefficients from biorthogonal left modes:

```text
c_j = <phi_j^dagger, T psi0>.
```

Use:

```python
modal = ModalBasis.from_spectrum(spectrum, metric=T)
coeffs = modal.expand(psi0)
psi_t = modal.reconstruct(t, coeffs)
```

### 36.5 Full eigenspectrum API

```python
spectrum = EigenSpectrumSolver(
    method="krylov_schur",
    nev=40,
    target="rightmost",
    shift=0.0,
).solve(alpha_problem)
```

Targets:

```text
largest_magnitude
smallest_magnitude
rightmost
leftmost
nearest_shift
inside_contour
```

Use `rightmost` for alpha stability. Use `largest_magnitude` or `dominant` for k multiplication operator.

---

## 37. Cross-section and material naming

A material object should not be just a dict. It carries operator coefficients.

```python
@dataclass(frozen=True, slots=True)
class Material:
    name: str
    total: TotalCrossSection
    absorption: AbsorptionCrossSection
    scattering: LegendreMomentScatteringData
    fission: FissionData | None
    velocity: NeutronVelocity | None
```

High-signal classes:

```text
TotalCrossSection
AbsorptionCrossSection
RemovalCrossSection
LegendreMomentScatteringData
GroupTransferMatrix
FissionSpectrum
NuFissionCrossSection
DelayedNeutronData
NeutronVelocity
MajorantCrossSection
```

Avoid:

```python
xs["s"]
xs["f"]
data[3]
```

Prefer:

```python
xs.scattering.moment(ell).matrix[to_group, from_group]
xs.fission.nu_sigma_f[group]
xs.fission.chi[to_group]
xs.velocity[group]
```

If using arrays, wrap them with axis labels.

```python
@dataclass(frozen=True, slots=True)
class GroupTransferMatrix:
    data: Array
    to_groups: EnergyGroupSpace
    from_groups: EnergyGroupSpace

    def apply(self, group_flux):
        ...
```

---

## 38. Operator contexts versus operators

Do not confuse an operator with the cached information needed to apply it.

Example:

```python
S = ScatteringOperator(xs.scatter)
context = S.prepare(method_space)
```

The operator is semantic. The context is implementation-specific.

```text
ScatteringOperator
    semantic physics object

SNScatteringContext
    harmonic projection matrices, moment caches

PNScatteringContext
    ell-block diagonal factors

MCScatteringContext
    sampled CDFs and alias tables

TTScatteringContext
    tensor factors and TT ranks
```

Pattern:

```python
class OperatorContext(ABC):
    method_space: MethodSpace

class LinearOperator(ABC):
    def prepare(self, method_space: MethodSpace) -> OperatorContext:
        ...

    def apply_prepared(self, context: OperatorContext, field: Field) -> Field:
        ...
```

This avoids polluting the physics operator with every method-specific cache.

---

## 39. Error and diagnostic naming

Use errors that tell the debugging agent what mathematical contract failed.

```text
DomainMismatchError
CodomainMismatchError
IncompatibleSpaceError
AdjointConsistencyError
ConservationViolation
PositivityViolation
NormalizationViolation
ExactnessViolation
SymmetryViolation
ReciprocityViolation
BoundaryResolutionError
TraceSpaceMismatch
SweepGraphCycleError
MajorantViolation
ProjectionConsistencyError
EigenResidualError
```

Example:

```python
raise MajorantViolation(
    message="majorant_cross_section does not dominate total cross section",
    max_violation=max(xs.total - majorant),
    cell=bad_cell,
    group=bad_group,
)
```

High-signal diagnostics:

```python
diagnostic = ConservationDefect(operator=S, value=defect, norm="l1")
```

Do not raise generic `ValueError` when a mathematical law is violated.

---

## 40. Codex implementation roadmap

A practical phased implementation plan:

### Phase 0: naming and law infrastructure

Implement:

```text
Axis
AxisProduct
RegistryMixin
Diagnostic
MathematicalLawError
```

### Phase 1: core math

Implement:

```text
Domain
Map
Measure
DiscreteMeasure
ProductMeasure
PushforwardMeasure
Space
TensorProductSpace
DirectSumSpace
Field
LinearOperator
Functional
ProjectionOperator
Kernel
```

Tests:

```text
measure normalization
operator linearity
operator adjoint consistency
tensor product shape/domain laws
```

### Phase 2: geometry and mesh

Implement:

```text
GeometrySpec
BoundaryDeclaration
MaterialAssignment
SpatialMesh
CellComplex
BoundaryPatchMap
Mesher
```

Tests:

```text
cell volumes positive
face adjacency consistent
boundary patches cover boundary
material map total over cells
```

### Phase 3: angular and energy structures

Implement:

```text
CubatureRule
LebedevCubature
ProductAngularCubature
SphericalHarmonicSpace
HarmonicMomentProjection
EnergyGroupSpace
LegendreMomentScatteringData
```

Tests:

```text
cubature exactness
harmonic orthogonality
SN-PN projection/reconstruction consistency
energy group axis consistency
```

### Phase 4: physics operators

Implement:

```text
TimeMassOperator T
StreamingOperator L
CollisionOperator C
ScatteringOperator S
FissionOperator F
Source
BoundaryOperator
```

Tests:

```text
scattering positivity/conservation
fission positivity
boundary trace mapping
operator adjoints
```

### Phase 5: method spaces

Implement:

```text
DiscreteOrdinatesPhaseSpace
SphericalHarmonicMomentSpace
CharacteristicTrackSpace
CollisionProbabilityRegionSpace
ParticleTrackingSpace
StochasticGalerkinPhaseSpace
```

Tests:

```text
boundary resolution per method
precompute validity
method-space axis naming
```

### Phase 6: solvers

Implement:

```text
SourceIteration
SweepSolver
GMRES
PowerIteration
ShiftInvertAlphaSolver
MonteCarloSolver
StochasticGalerkinSolver
```

Tests:

```text
manufactured fixed-source cases
criticality benchmark residual
alpha mode residual
MC estimator sanity cases
```

### Phase 7: tensor networks and hybrid methods

Implement:

```text
SumOfTensorProductsOperator
TensorNetworkOperator
TensorTrainOperator
BoundaryTensorNetwork
ControlVariateHybridSolver
AdjointImportance
```

Tests:

```text
tensor product application equals dense small case
TT compression residual
hybrid control-variate unbiasedness on known case
```

---

## 41. Development style rules for Codex

1. **Never create a naked array without axis metadata unless it is inside a very local numerical kernel.**

2. **Every operator must know its domain and codomain.**

3. **Every projection must expose its reconstruction or adjoint relationship.**

4. **Every boundary condition must resolve through a trace space.**

5. **Every method-specific object should be named by the representation it creates, not by the fact that it is augmented.**

6. **Every cache class should be explicitly named `...Cache` or `...Context`.**

7. **Every mathematical adjective in a class name must have a residual and assertion method.**

8. **Do not duplicate physics between deterministic and MC. Use the same kernel data to build deterministic operators and stochastic samplers.**

9. **Prefer `A_h = P @ A @ R` or `A_h = R.H @ A @ R` over hand-coded collapsed formulas.**

10. **Prefer exact separable/tensor-product structure before applying tensor-train compression.**

11. **Use `T` only for the time mass/time derivative coefficient. Use `L` for streaming.**

12. **Do not bury sign conventions. `AlphaEigenproblem` must document whether positive alpha means growth.**

13. **Use registries for discoverability but dependency injection for tests.**

14. **Keep `GeometrySpec` immutable. Method-specific resolution belongs later.**

15. **Separate semantic operators from implementation contexts.**

---

## 42. Compact glossary of high-signal terms

```text
Adjoint
    Hilbert-space dual operator. Enables response theory, sensitivities, and MC importance.

AlbedoKernel
    Boundary response kernel mapping outgoing to incoming trace states.

AlphaEigenvalue
    Time-growth eigenvalue; pole of Laplace-domain resolvent.

BoundaryGeometryMap
    Geometry-only map between boundary trace indices/states.

BoundaryResponseKernel
    Physics-only response at boundary; separable from geometry.

BoundaryTensorNetwork
    Tensor factorization of boundary response over patch, face, angle, energy, stochastic axes.

CharacteristicFlow
    Flow generated by streaming vector field dx/ds = Omega.

CollisionProbabilityKernel
    Region-projected Green kernel with positivity, conservation, and reciprocity properties.

Conservative
    Preserves total particle balance under stated measure.

CubatureRule
    Discrete measure for multidimensional integration.

DiscreteOrdinatesPhaseSpace
    Mesh-energy-angular collocation space for SN.

EmpiricalMeasure
    Weighted samples; MC particle ensemble as measure.

FissionBank
    Particles born from fission and queued for next generation/census.

GalerkinProjection
    Projection induced by test/trial spaces; core of discretization and condensation.

InvariantCubature
    Cubature rule invariant under a symmetry group.

KineticGenerator
    Time evolution generator, often -T^{-1} A_prompt.

LaplaceResolvent
    Resolvent (sT + A)^{-1}; poles are alpha modes.

MajorantCrossSection
    Dominating cross section for Woodcock/delta tracking.

MarkovKernel
    Positive normalized transition kernel.

MethodSpace
    Method-specific representation space built from a mesh.

MomentProjection
    Projection from angular values to harmonic moments.

ParticleEnsemble
    Weighted empirical collection of particle states/histories.

Projection
    Maps rich space to reduced/coefficient/coarse space.

RayBundle
    Family of characteristic tracks for one or more directions with transverse measure.

Reciprocal
    Satisfies detailed-balance-like symmetry, common in CP kernels.

SphericalHarmonicMomentSpace
    PN angular moment representation with SO(3) irrep blocks.

SumOfTensorProducts
    Exact/semi-exact separable operator representation.

TensorTrain
    Compressed tensor-network representation for high-dimensional tensor-product spaces.

EquivariantSectorSpace
    Method space on a fundamental domain whose operators and traces commute with a symmetry group action.

FundamentalDomain
    Representative sector of a group action used to reconstruct the full domain by lifting.

InvariantProjector
    Group-average projector onto fields satisfying rho(g) psi = psi.

OrbifoldPhaseSpace
    Quotient phase space with stabilizers/isotropy strata and orbit-multiplicity measures.

OrbitMultiplicityField
    Field containing |G| / |Stabilizer(x)| factors for lifting sector integrals to full-domain integrals.

PhysicalBoundaryPatch
    True external boundary where vacuum, inflow, albedo, or leakage semantics may apply.

QuotientMeasure
    Sector measure corrected by orbit multiplicity to represent full-domain integrals.

SectorLift
    Reconstruction from sector fields to equivariant full-domain fields.

SectorRestriction
    Restriction from full-domain fields to the chosen fundamental sector.

SymmetryBoundary
    Quotient trace-gluing boundary induced by a group generator; not a physical mirror or vacuum wall.

SymmetryBoundaryPatch
    Boundary patch of a fundamental domain produced by quotienting under a group action.

TraceSpace
    Boundary-restricted phase space for inflow/outflow states.
```

---

## 43. Closing implementation principle

The code should read like applied mathematics:

```python
V = X * Omega * G * Xi

A_loss = L + C - S
A_prompt = A_loss - F

k_problem = CriticalityProblem(A_loss, F)
alpha_problem = AlphaEigenproblem(T, A_prompt)

resolvent = LaplaceResolvent(T, A_prompt)

S = Scattering.from_legendre_moments(xs.scatter).on(V)
B = BoundaryResolution.from_geometry(geometry, mesh, V)

psi = solver.solve(FixedSourceProblem(A_prompt, q, B))
```

A reader should be able to infer the mathematics from the names before reading the documentation. Documentation then becomes confirmation, not rescue.
