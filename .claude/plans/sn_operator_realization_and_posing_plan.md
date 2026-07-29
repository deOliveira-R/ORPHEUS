# SN operator realization and strategic posing plan

**Status:** design and migration plan only  
**Date:** 2026-07-24  
**Primary scope:** `transport`, `sn`, and the supporting `numerics` algebra  
**Execution scope for the next session:** implement in small, gated phases; do not
attempt a flag-day rewrite  
**Grand-report relationship:** this plan refines, and in a few places supersedes,
the current wording in `neutron_transport_grand_report_v3.md`. In particular, it
sharpens what the report calls a “semantic operator” into a distinction between
a physical kernel and a realized, fixed-domain operator.

---

## 1. Executive decision

The proposed rename is directionally correct, but the new object should not be a
polymorphic “operator that returns another operator” and should not know every
transport method.

Use three layers:

1. **Physical transfer kernel**
   - `FissionKernel`
   - `LegendreScatteringKernel`
   - later, separately, `N2NKernel`
   - non-callable;
   - contains or views the method-neutral reaction data and its mathematical
     invariants;
   - may be spatially discretized through `MaterialXSField`, but does not know an
     angular quadrature, SN field, diffusion trace, sweep, or solver.

2. **Method-neutral transfer operators on their natural reduced spaces**
   - `FissionEnergyTransferOperator`:
     `ScalarFlux -> ScalarSourceSink`;
   - `LegendreMomentTransferOperator`:
     `HarmonicMomentFlux -> HarmonicMomentSourceSink`;
   - `IsotropicScattering` and `IsotropicN2N` remain reduced-space operators,
     subject to the naming review in §9.

3. **Method-owned realized operators**
   - `SNFissionOperator`;
   - `SNScatteringOperator`;
   - diffusion uses a scalar/full-field realization;
   - homogeneous transport uses the reduced energy operator directly;
   - each has exactly one domain and one codomain;
   - each is built by an explicit method-owned realizer, preferably as a
     composition of analysis, transfer, and synthesis maps.

The construction seam should follow the successful boundary pattern already in
the repository:

```text
BoundaryTraceLaw -> SNBoundaryRealizer -> fixed-signature boundary operator
FissionKernel    -> SNReactionRealizer -> fixed-signature fission operator
ScatteringKernel -> SNReactionRealizer -> fixed-signature scattering operator
```

Do **not** add:

- `kernel.realize("SN")`;
- a global string-keyed realizer registry;
- `if method == ...` dispatch in `transport`;
- `singledispatchmethod` on the realized operator’s input carrier;
- a common object that imports SN, diffusion, MoC, CP, and MC.

The physical kernel belongs to `transport`; the realization belongs to the
transport method. This preserves the dependency direction and lets a new method
add its interpretation without modifying the common kernel.

The short form of the naming rule is:

> Reserve the suffix `Operator` for an apply-capable mathematical map with one
> declared domain and one declared codomain. Use `Kernel`, `Law`, `Data`,
> `Partition`, `Posing`, `Schedule`, or `Realizer` for the other concepts.

---

## 2. Why the current polymorphism is mathematically unfaithful

Polymorphism itself is not the problem. The problem is using runtime carrier
dispatch to make one instance stand for several different mathematical arrows.

Today `FissionOperator.apply` accepts, among other carriers:

```text
ndarray    -> ndarray
ScalarFlux -> ScalarSourceSink
FullField  -> FullField
```

and `ScatteringOperator` has a still broader surface involving scalar, angular,
moment, array, and composite carriers. Its optional `full_field_space` attempts
to advertise one domain/codomain while the dispatch arms actually implement
several distinct maps.

Those arms are related, but they are not the same operator. They include:

- angular integration;
- energy transfer;
- isotropic or harmonic reconstruction;
- bulk extraction from a direct sum;
- zero-source injection into a boundary component;
- raw-array adaptation.

In mathematics, those are separate arrows whose composition forms a realized
operator. Treating them as overloads hides the factorization, weakens static
typing, makes adjoints harder to audit, and makes a method-specific object appear
method-neutral.

The better use of polymorphism is at the **construction boundary**:

- a realizer protocol describes “kernel plus method space yields an operator”;
- each concrete method owns a concretely typed implementation;
- the returned operator is ordinary fixed-signature `LinearOperator[D, C]`;
- solver code need not know how the realization was constructed.

This gives extension polymorphism without application polymorphism.

---

## 3. Mathematical model and vocabulary

### 3.1 Kernel versus operator

A transfer kernel describes how source state `z'` contributes to target state
`z`. An operator is obtained after choosing spaces, measures, and a
representation:

```text
(K psi)(z) = integral kappa(z <- z') psi(z') dmu(z').
```

For ORPHEUS, `MaterialXSField` already discretizes the spatial/material
dependence. Therefore the proposed kernels are not completely continuous laws;
they are more precisely **angular-method-neutral, spatially represented reaction
kernels**. The code and documentation should say this explicitly.

### 3.2 Fission

At each spatial point or cell, the multigroup fission kernel is

```text
K_f(g <- g'; x) = chi_g(x) nuSigma_f,g'(x).
```

It is rank one on the energy axis:

```text
F_G = |chi><nuSigma_f|.
```

The natural reduced-space operator is:

```text
F_G : ScalarFlux -> ScalarSourceSink.
```

For SN, let:

```text
M_0^SN : AngularFlux -> ScalarFlux
E_0^SN : ScalarSourceSink -> AngularSourceSink
```

where `M_0^SN` is quadrature-weighted angular analysis and `E_0^SN` is
the isotropic source lift under the project’s normalization convention. Then:

```text
F_SN = E_0^SN F_G M_0^SN.
```

On the composite bulk-plus-boundary carrier, introduce explicit direct-sum maps:

```text
pi_bulk : FullFluxField -> AngularFlux
j_bulk  : AngularSourceSink -> FullSourceField
```

where `j_bulk` supplies the mathematically zero boundary-source block. The full
operator is:

```text
F_SN,full = j_bulk E_0^SN F_G M_0^SN pi_bulk.
```

The zero boundary result is therefore not a hidden branch inside fission; it is
the action of a named injection.

### 3.3 Scattering

Let:

```text
M_L^SN : AngularFlux -> HarmonicMomentFlux
Lambda : HarmonicMomentFlux -> HarmonicMomentSourceSink
R_L^SN : HarmonicMomentSourceSink -> AngularSourceSink.
```

`Lambda` contains the Legendre-order blocks and group transfer:

```text
Lambda = direct_sum over ell of Sigma_s,ell.
```

The SN scattering realization is:

```text
S_SN = R_L^SN Lambda M_L^SN.
```

The full-field version is:

```text
S_SN,full = j_bulk R_L^SN Lambda M_L^SN pi_bulk.
```

The existing harmonic `Frame` already expresses much of
`R_L^SN Lambda M_L^SN`. The plan is to make those factors the public
algebra-of-record instead of leaving them as a property of a broadly
polymorphic facade.

### 3.4 Secondary-neutron production is not silently “scattering”

The current `ScatteringOperator` also incorporates `(n,2n)`. That is a useful
solver posing, but it is not a faithful reason to call the whole object a
scattering operator.

Use separate symbols and leaves:

```text
S      = physical scattering transfer
P_2n   = secondary-neutron production transfer
G_col  = S + P_2n        # if a chosen posing treats both as RHS gains
```

The initial migration may preserve numerical grouping, but the factor tree and
names must retain this distinction. `SNScatteringOperator` must not quietly mean
“scattering plus all other neutron-producing collision channels.” If a combined
object is valuable, name it `SNCollisionGainOperator` or expose the explicit
`OperatorSum`.

### 3.5 Operator, equation, posing, realization, and schedule

Use these terms consistently:

- **kernel**: physical transfer relation before a method-specific carrier is
  selected;
- **operator**: a fixed-domain/fixed-codomain linear map;
- **equation form**: symbolic relation among operators, e.g.
  `A psi = q + G psi`;
- **operator pencil**: a parameterized pair such as `A - lambda F`;
- **splitting**: algebraic identity `A = M - N`;
- **posing**: a selected equivalent form of the mathematical problem;
- **partition**: decomposition of a variable or operator space into subspaces;
- **block realization**: `A_ij = R_i A J_j` under a partition;
- **dependency graph**: graph induced by nonzero block couplings;
- **schedule**: executable order for applying or solving blocks;
- **realizer**: turns a kernel/law and method space into fixed-signature
  operators;
- **lowering**: turns an algebraic posing into a concrete executable plan;
- **prepared context**: cached data used by an operator implementation;
- **materialization**: conversion to an explicit matrix/tensor;
- **solver**: consumes a linear system, splitting, or pencil and produces an
  approximate or exact solution.

“Strategy operator” should not become a base class. A strategy is often a
transformation of an equation or a selection of a splitting, not a linear map.
Use a `...Posing`, `...Splitting`, `...Partition`, or `...Schedule` object as
appropriate.

---

## 4. Relation to generalized eigenvalue problems and RHS ergonomics

ORPHEUS’s criticality equation is naturally a generalized eigenvalue problem:

```text
A psi = (1/k) F psi
```

or:

```text
(A - lambda F) psi = 0,   lambda = 1/k.
```

The primary algebraic object should therefore be a pencil:

```python
GeneralizedEigenPencil(lhs=A, rhs=F, eigenvalue_parameter="inverse_keff")
```

The pencil contract is:

```text
domain(A)   == domain(F)
codomain(A) == codomain(F)
```

The common domain contains the mode/flux unknown and the common codomain contains
the residual/source balance. The two operators need not be encoded as
untyped array endomorphisms; flux-role to source/residual-role maps are valid and
more informative, provided both sides of the equation use the same roles and
spaces.

Forming:

```text
K = A^-1 F
```

is a **solver posing**, not the definition of the problem. It is convenient for
power iteration, but preserving `(A, F)` enables:

- generalized Arnoldi/Krylov backends;
- shift-and-invert formulations;
- direct generalized eigensolvers after materialization;
- left/adjoint modes without reverse-engineering a fused map;
- alpha-eigenvalue pencils involving the time-mass operator;
- sensitivity formulas that separately differentiate `A` and `F`;
- preconditioners for `A - sigma F`;
- clear normalization functionals and physical diagnostics.

The LHS investment and RHS investment therefore meet at the pencil boundary.
The RHS is not merely a source callback: it is a first-class operator with
domain/codomain, adjoint, factorization, provenance, and potentially a prepared
context.

Target objects:

```python
@dataclass(frozen=True)
class GeneralizedEigenPencil:
    lhs: LinearOperator
    rhs: LinearOperator
    spectral_parameter: SpectralParameter

@dataclass(frozen=True)
class PowerMapPosing:
    pencil: GeneralizedEigenPencil
    lhs_inverse: LinearOperator

    def iteration_operator(self) -> LinearOperator:
        return self.lhs_inverse @ self.pencil.rhs
```

The solver should choose the posing. The equation builder should not eagerly
collapse the pencil into `A.inverse() @ F`.

---

## 5. Target architecture

### 5.1 Layer placement

The dependency direction remains:

```text
numerics
   ^
transport (method-neutral spaces, fields, physical kernels, reduced operators)
   ^
sn / diffusion / moc / cp / mc (method realizers and method operators)
   ^
solver orchestration and applications
```

Proposed modules:

```text
orpheus/transport/kernels/
    __init__.py
    fission.py
    scattering.py
    secondary_production.py          # may land after the initial carve

orpheus/transport/operators/
    fission_energy_transfer.py
    moment_transfer.py
    isotropic_scattering.py          # existing, reviewed rather than duplicated

orpheus/sn/operators/
    angular_analysis.py              # only if existing Frame maps cannot be reused
    reaction.py                      # SNFissionOperator, SNScatteringOperator
    reaction_realizer.py             # SNReactionRealizer
    discrete_operators.py            # assembled SN operator record

orpheus/numerics/posing/
    __init__.py
    splitting.py
    eigen_pencil.py
    partition.py
    dependency_graph.py
```

Do not create all modules up front if a smaller additive carve is clearer. The
final dependency boundaries matter more than the first commit’s file count.

### 5.2 Kernel objects

Proposed shape:

```python
@dataclass(frozen=True)
class FissionKernel:
    mat_xs: MaterialXSField

    @property
    def emission_spectrum(self) -> Array: ...

    @property
    def production_covector(self) -> Array: ...

    def validate(self) -> None: ...


@dataclass(frozen=True)
class LegendreScatteringKernel:
    mat_xs: MaterialXSField
    order: int

    def moment_transfer(self, ell: int) -> GroupTransferMatrix: ...
    def validate(self) -> None: ...
```

They do not implement `apply`, do not inherit `LinearOperator`, and do not carry
quadrature or `FullFieldSpace`.

The kernel should be a zero-copy semantic view over `MaterialXSField` unless
immutability/caching tests show that snapshot semantics are required. Depletion
and feedback currently depend on read-through behavior; do not accidentally
freeze cross-section arrays.

### 5.3 Reduced transfer operators

Proposed shape:

```python
class FissionEnergyTransferOperator(
    LinearOperator[ScalarFlux, ScalarSourceSink]
):
    kernel: FissionKernel
    domain: ScalarFluxSpace
    codomain: ScalarSourceSpace


class LegendreMomentTransferOperator(
    LinearOperator[HarmonicMomentFlux, HarmonicMomentSourceSink]
):
    kernel: LegendreScatteringKernel
    domain: HarmonicMomentFluxSpace
    codomain: HarmonicMomentSourceSpace
```

Raw ndarray execution may exist as a private kernel/context primitive used by
the typed operator, but it should not be a second public `apply` signature.
Where array compatibility is needed during migration, provide an explicitly
named adapter such as `apply_values` or `FlattenedOperator`.

### 5.4 SN realization

Proposed shape:

```python
@dataclass(frozen=True)
class SNReactionRealizer:
    method_name: str = "SN"

    def realize_fission(
        self,
        kernel: FissionKernel,
        method_space: SNMethodSpace,
    ) -> SNFissionOperator: ...

    def realize_scattering(
        self,
        kernel: LegendreScatteringKernel,
        method_space: SNMethodSpace,
    ) -> SNScatteringOperator: ...
```

The realizer should build the operator from named factors. It should not own
mutable caches that logically belong to the operator or a prepared context.

The realized operator may be implemented in either of these forms:

1. a thin named wrapper around `OperatorProduct`; or
2. a fused implementation carrying the factor tree as its algebraic
   specification/provenance.

The plan selects option 1 first. Only introduce fusion after profiling. The
unfused composition remains the equivalence oracle even if a fused kernel is
later used in production.

### 5.5 Assembled SN operator record

Replace the informal solver attributes `L`, `S`, and `F` with a typed assembly
record:

```python
@dataclass(frozen=True)
class SNDiscreteOperators:
    streaming: SNStreamingOperator              # sigma-free L
    collision: MultiplicationOperator           # C
    free_transport: SNFreeTransportOperator     # T = L + C, sweep-invertible
    scattering: SNScatteringOperator            # S
    fission: SNFissionOperator                  # F
    boundary: SNBoundaryOperator                # B
```

This object is a realized operator collection, not a new monolithic operator and
not a physics owner. It may expose checked convenience compositions:

```python
operators.loss_without_boundary
operators.criticality_pencil(...)
```

but the primitive leaves remain inspectable and shared by identity.

`SNOperatorFamily` is not recommended: “family” risks repeating the current
ambiguity. `SNDiscreteOperators` or `SNOperatorAssembly` says what the object is.

---

## 6. Specific current-code findings

### 6.1 `FissionOperator`

Current file: `orpheus/transport/operators/fission.py`.

Keep:

- the rank-one energy dyad;
- the read-through behavior from `MaterialXSField`;
- transpose/adjoint support;
- the scalar source role;
- the common physical placement under `transport`.

Carve apart:

- kernel data and invariants;
- scalar energy transfer;
- angular analysis;
- isotropic SN synthesis;
- direct-sum bulk extraction/injection;
- ndarray compatibility.

The current `.kernel` property is already close to the reduced energy-transfer
operator. It should become a first-class fixed-domain object instead of a helper
property behind a polymorphic facade.

### 6.2 `ScatteringOperator`

Current file: `orpheus/transport/operators/scattering.py`.

The class currently owns too many concepts:

- Legendre transfer data;
- quadrature;
- harmonic frame construction;
- moment transfer;
- P0 fast paths;
- `(n,2n)` production;
- source accumulation helpers;
- composite-field adaptation;
- transpose composition;
- the “foldable/residual” data split.

The class is method-specific in fact because a discrete angular frame requires a
quadrature. Its current location and name imply broader method neutrality than
its state permits. The grand report’s general principle is better served by
moving the frame-bearing realization into `sn`, while retaining transfer kernels
and reduced operators in `transport`.

Rename `LegendreMomentScattering` to `LegendreMomentTransferOperator`; its current
typed role change from flux moments to source moments should be made explicit in
its domain/codomain types.

### 6.3 Collision

`MultiplicationOperator` is genuinely reusable. Collision as local
multiplication by `Sigma_t` does not need method realization when it acts on a
compatible carrier and space. Do not force all transport physics through the
kernel/realizer pattern. Use the pattern where representation actually changes.

### 6.4 `SNSolver.self.L`

Current construction:

```python
self.L = build_streaming_collision(...)
```

but `build_streaming_collision` returns `StreamingOperator +
MultiplicationOperator`, i.e. `L + C`. This violates the report’s own convention
that `L` is sigma-free streaming.

Target:

```python
self.operators.streaming
self.operators.collision
self.operators.free_transport
```

During migration, direct solver aliases may be:

```python
self.streaming = operators.streaming
self.free_transport = operators.free_transport
self.S = operators.scattering       # transitional only
self.F = operators.fission          # transitional only
```

Do not leave `self.L` bound to `L + C`. If a compatibility property is required,
make `self.L` return the sigma-free streaming leaf and fail tests that expect
otherwise.

Rename the current generic concrete `InvertibleOperator` in
`sn/operators/streaming.py` to `SNFreeTransportOperator` or
`SNStreamingCollisionOperator`. This plan selects `SNFreeTransportOperator`
because `T_0 = L + C` is the free-transport-with-attenuation operator whose
inverse is the sweep/propagator. Avoid `TransportOperator`, which is too broad.

### 6.5 `WithinGroupSystem.resolvent`

`WithinGroupSystem.resolvent` currently stores the implicit member `M` of
`A = M - N`, not `M^-1`. A resolvent is an inverse-like object, commonly
`(zI - A)^-1`; the current name is misleading.

Rename:

```text
resolvent -> implicit_operator
gains     -> explicit_gains
```

When an inverse wrapper is actually constructed, call that
`implicit_inverse`, `preconditioner`, or `resolvent` according to its precise
mathematical role.

### 6.6 “Regular splitting”

The code calls the decomposition a regular splitting. In numerical linear
algebra, regular splitting usually carries positivity conditions such as
`M^-1 >= 0` and `N >= 0`. Unless ORPHEUS verifies the required cone and
positivity properties, use:

```text
stationary splitting
implicit/explicit splitting
transport splitting
```

Reserve `regular splitting` for a checked subtype or theorem-backed condition.

### 6.7 Existing `CoupledOperator` versus a block view

`CoupledOperator` is a concrete explicit block grid over already identified
systems. It is useful and should remain.

The proposed block view is different:

```text
given A and a partition {V_i},
A_ij = R_i A J_j
```

where `R_i` is restriction and `J_j` injection. This view:

- starts from one monolithic operator;
- derives blocks lazily;
- permits partition by ordinate, group, carrier, region, or another axis;
- can expose graph structure without duplicating operator construction;
- may lower to `CoupledOperator` only when explicit blocks are needed.

Use a distinct name such as `PartitionedOperatorView`. Do not overload
`CoupledOperator` with both responsibilities.

---

## 7. SN systems, partitions, and graph theory

### 7.1 One SN equation, many valid partitions

The monolithic SN state may be partitioned by:

- physical carrier (the present System A/System B coupling);
- energy group;
- ordinate or ordinate bundle;
- spatial sweep region;
- harmonic moment block;
- combinations of the above.

These partitions do not necessarily define different physical “systems.” Use:

```text
VariablePartition
OperatorPartition
SubsystemView
BlockLayout
```

and reserve `CoupledSystem` for a direct sum whose members have genuine,
independently typed carrier identities.

There is also an existing symbol/name collision to retire:

```text
A                   commonly means the assembled LHS/loss operator
B                   commonly means the boundary return operator
System A / System B currently name two radial-characteristic carriers
```

The carrier systems should receive descriptive public names before partitioning
is generalized. Candidate vocabulary, to be checked against their exact spaces
in Phase 0:

```text
BulkAngularSystem / RadialCharacteristicSystem
PrimaryFluxSystem / RayHistorySystem
```

Do not extend the lettered `SystemRole.A` / `SystemRole.B` convention to energy
groups or ordinates. A group or ordinate is a partition member, not another
lettered physical system.

### 7.2 Lazy block construction

Introduce:

```python
@dataclass(frozen=True)
class OperatorPartition:
    domain_parts: tuple[Subspace, ...]
    codomain_parts: tuple[Subspace, ...]


class PartitionedOperatorView:
    operator: LinearOperator
    partition: OperatorPartition

    def block(self, i: int, j: int) -> LinearOperator:
        return restriction(i) @ operator @ injection(j)
```

Require partitions to provide:

- stable part identities and labels;
- domain/codomain spaces;
- restriction and injection operators;
- completeness and disjointness checks where applicable;
- a zero/nonzero structural predicate that does not require dense
  materialization.

### 7.3 Dependency graphs and SCC condensation

For a chosen partition, define a directed dependency graph:

```text
j -> i iff A_ij or gain_ij can make subsystem i depend on subsystem j.
```

Use standard graph-theoretic objects:

- strongly connected components (SCCs);
- condensation graph;
- topological ordering of the condensation DAG;
- feedback edges;
- block triangular form;
- Dulmage–Mendelsohn decomposition later if rectangular structural systems need
  it;
- elimination tree only when sparse factorizations/materialized matrices enter.

The key theorem for scheduling is:

> Contracting SCCs always produces a DAG. Acyclic dependencies can be ordered
> exactly; cyclic dependencies remain inside SCC blocks and require an inner
> solve, iteration, lagging, or a larger implicit block.

This is the faithful replacement for an informal global
“foldable versus lagged” label.

### 7.4 Target-relative foldability

An edge is not intrinsically foldable. It is foldable **with respect to**:

- a chosen partition;
- a chosen implicit operator;
- a chosen ordering;
- available solve capabilities;
- fill/memory constraints;
- a numerical policy.

Represent the decision:

```python
@dataclass(frozen=True)
class StationarySplitting:
    equation_operator: LinearOperator   # A
    implicit_operator: LinearOperator   # M
    explicit_operator: LinearOperator   # N
    # invariant: A == M - N


@dataclass(frozen=True)
class DependencySchedule:
    sccs: tuple[SubsystemBlock, ...]
    condensation_order: tuple[int, ...]
    implicit_sccs: frozenset[int]
    lagged_edges: tuple[DependencyEdge, ...]
```

The scheduler may choose to enlarge implicit SCCs, but the graph analyzer should
not silently make that numerical-policy decision.

### 7.5 Current scattering “foldable part”

`ScatteringOperator.foldable_part()` extracts diagonal P0 self-scatter data.
That is a valid **reaction-kernel decomposition**, but calling it foldable
suggests an exact algebraic equivalence that the current source comments already
warn is false for anisotropic flux when naively absorbed into a removal
coefficient.

Rename toward the data fact:

```text
within_group_p0_self_scatter()
complementary_transfer()
```

or:

```text
self_transfer_component()
offdiagonal_transfer_component()
```

Then let a separate posing decide whether and how that component enters an
implicit operator, DSA preconditioner, or lagged gain. The data object must not
promise solver foldability.

### 7.6 Boundary cycles

Boundary coupling is a natural first consumer of SCC analysis:

- vacuum/inflow produces no feedback edge;
- reflective or periodic maps can close a sweep path;
- acyclic trace dependencies can be folded into a schedule;
- cyclic trace dependencies form SCCs that require a coupled solve or lag.

The existing `creates_sweep_cycle` boolean is useful as a conservative hint, but
the durable result should come from the realized graph under the actual mesh,
ordinates, and boundary maps. A law alone cannot know the full realized cycle
topology.

---

## 8. Strategic posing pipeline

Do not model the whole pipeline as operators. Use typed transformations:

```text
physical kernels/laws
        |
        v
method realization
        |
        v
SNDiscreteOperators                         fixed-domain leaves
        |
        v
TransportEquation / GeneralizedEigenPencil exact mathematical problem
        |
        v
Posing                                     equivalent algebraic form
        |
        v
Partition + dependency analysis            block structure and SCCs
        |
        v
Splitting + schedule                        implicit/explicit execution plan
        |
        v
Prepared/realized solve plan                caches, sweep kernels, factors
        |
        v
solver
```

Suggested records:

```python
@dataclass(frozen=True)
class FixedSourceEquation:
    lhs: LinearOperator
    source: Field


@dataclass(frozen=True)
class CriticalityEquation:
    loss: LinearOperator
    production: LinearOperator

    def as_pencil(self) -> GeneralizedEigenPencil: ...


@dataclass(frozen=True)
class WithinGroupPosing:
    equation: FixedSourceEquation
    splitting: StationarySplitting
    partition: OperatorPartition | None
    schedule: DependencySchedule
```

The final cached object may be named:

```text
PreparedSNSolve
SNSolvePlan
RealizedSNSolvePlan
```

Prefer `SNSolvePlan` when it contains executable strategies plus caches.
“Fully realized operator” is too narrow if the object also contains orderings,
stopping policies, preconditioners, and mutable cache state.

---

## 9. Naming crosswalk

| Current or informal name | Transitional name | Final intended meaning |
|---|---|---|
| `FissionOperator` | compatibility facade | retire; it currently denotes several maps |
| — | `FissionKernel` | non-callable rank-one physical transfer kernel |
| current `.kernel` dyad | `FissionEnergyTransferOperator` | scalar flux to scalar fission source |
| fission FullField dispatch arm | `SNFissionOperator` or scalar composite realization | one method/carrier-specific map |
| `ScatteringOperator` | compatibility facade | retire; currently kernel + SN frame + assembly + n2n |
| — | `LegendreScatteringKernel` | non-callable Legendre/group transfer data |
| `LegendreMomentScattering` | `LegendreMomentTransferOperator` | moment flux to moment source |
| current frame conjugation | `SNScatteringOperator` | fixed-signature SN angular scattering map |
| current scattering+n2n sum | `SNCollisionGainOperator` or explicit sum | posing-level combined RHS gain |
| `foldable_part` | `within_group_p0_self_scatter` | data decomposition only |
| `residual_part` | `complementary_transfer` | the remaining transfer data |
| `StreamingOperator` | `SNStreamingOperator` | sigma-free `L` |
| `InvertibleOperator` | `SNFreeTransportOperator` | sweep-invertible `L + C` |
| `SNSolver.L` containing `L+C` | `SNSolver.free_transport` | correct semantic name |
| `WithinGroupSystem.resolvent` | `implicit_operator` | `M` in `A=M-N` |
| `WithinGroupSystem.gains` | `explicit_gains` | explicit/lagged `N` pieces |
| “regular splitting” | `StationarySplitting` | no unproved positivity claim |
| proposed `SNOperatorFamily` | `SNDiscreteOperators` | realized leaf collection |
| `CoupledOperator` | unchanged | explicit block operator on a direct-sum space |
| proposed `BlockView` | `PartitionedOperatorView` | lazy `R_i A J_j` blocks |
| `System A` / `System B` | compatibility labels | descriptive carrier-system names |
| “data-focused operator” | `Kernel` or reduced operator | distinguish non-callable data from maps |
| “strategy operator” | `Posing`, `Splitting`, or `Schedule` | use the precise transformation concept |
| “fully realized operator” | `SNSolvePlan` when broader than a map | executable lowered plan and caches |

Two names should remain deliberately distinct:

- `FissionKernel` is physical transfer information.
- `FissionEnergyTransferOperator` is already an actual operator, even before it
  is lifted to SN.

---

## 10. Implementation phases

Each phase is independently mergeable and must leave the production suite
working. Compatibility facades are allowed temporarily; new code must not add
new dependencies on them.

### Phase 0 — freeze semantics and baseline behavior

1. Inventory every import, constructor, dispatch arm, helper method, and test for:
   - `FissionOperator`;
   - `ScatteringOperator`;
   - `LegendreMomentScattering`;
   - `InvertibleOperator`;
   - `SNSolver.L`;
   - `WithinGroupSystem.resolvent`;
   - `foldable_part` / `residual_part`.
2. Record baseline numerical outputs for:
   - representative 1-D and 2-D SN eigenvalue cases;
   - fixed-source SN;
   - anisotropic scattering;
   - radial-characteristic System B;
   - diffusion criticality;
   - homogeneous infinite medium;
   - adjoint/transpose tests.
3. Add architectural tests:
   - current class dispatch inventory;
   - import-layer checks;
   - domain/codomain assertions;
   - identity sharing of operator pieces in `WithinGroupSystem`.
4. Amend the grand-report task list with the decisions in §1, but do not rewrite
   historical sections until implementation proves the names.

**Exit gate:** exact baseline artifacts exist and all current consumers are
accounted for.

### Phase 1 — add non-callable reaction kernels

1. Add `FissionKernel` as a zero-copy view over `MaterialXSField`.
2. Add `LegendreScatteringKernel` without quadrature.
3. Define:
   - axis conventions;
   - group-to/group orientation;
   - fission rank-one factors;
   - Legendre order range;
   - cross-section read-through semantics;
   - validation behavior.
4. Add explicit tests that kernels:
   - do not inherit `LinearOperator`;
   - do not expose `apply`;
   - do not import method packages;
   - update consistently when supported cross-section data are rebound.
5. Do not change production callers.

**Exit gate:** kernels faithfully describe current reaction data with no solver
behavior change.

### Phase 2 — extract fixed-domain reduced operators

1. Promote the current fission dyad into
   `FissionEnergyTransferOperator`.
2. Rename/add `LegendreMomentTransferOperator`.
3. Give both concrete flux-role domains and source-role codomains.
4. Preserve transpose/adjoint operations from the factors.
5. Keep raw ndarray kernels private or behind `apply_values`.
6. Make existing `IsotropicScattering` and `IsotropicN2N` compose with these
   role-typed spaces; do not duplicate their math.
7. Gate bitwise or tolerance-equivalent results against every corresponding
   current dispatch arm.

**Exit gate:** the reduced operators are standalone truth, and the old facades
delegate their core math to them.

### Phase 3 — make analysis, synthesis, and direct-sum maps first-class

1. Reuse `HarmonicFrame.analysis` and `HarmonicFrame.reconstruction` where their
   domain/codomain and flux/source roles are correct.
2. Add only missing maps:
   - zeroth-moment angular analysis;
   - isotropic source synthesis under the project normalization;
   - bulk restriction;
   - bulk source injection with zero boundary component.
3. Ensure every map implements its transpose.
4. Verify:
   - quadrature normalization;
   - `R Lambda M` equivalence;
   - direct-sum block roles;
   - metric-aware adjoints, not merely Euclidean array transposes.
5. Relocate harmonic-frame construction away from the shared scattering facade
   to SN method space/realizer ownership.

**Exit gate:** all carrier conversions in fission and scattering are named
operators rather than hidden dispatch branches.

### Phase 4 — add SN reaction realization

1. Add `SNReactionRealizer`.
2. Implement:

   ```text
   SNFissionOperator   = j_bulk E0 F_G M0 pi_bulk
   SNScatteringOperator = j_bulk R Lambda M pi_bulk
   ```

3. Start as compositions/wrappers; do not optimize prematurely.
4. Expose factor/provenance inspection for diagnostics.
5. If `(n,2n)` remains on the same RHS schedule, construct an explicit
   `SNCollisionGainOperator = SNScatteringOperator + SNN2NOperator`.
6. Add full forward and adjoint equivalence tests against the old facade.

**Exit gate:** SN has fixed-signature reaction operators built entirely in the
SN layer from shared kernels/reduced operators.

### Phase 5 — introduce `SNDiscreteOperators` and fix `L`

1. Build sigma-free streaming, collision, `L+C`, scattering, fission, and
   boundary exactly once.
2. Store them in `SNDiscreteOperators`.
3. Rename `InvertibleOperator` to `SNFreeTransportOperator`, initially with a
   compatibility alias.
4. Replace `SNSolver.self.L = L+C` with:

   ```python
   self.operators = build_sn_discrete_operators(...)
   self.streaming = self.operators.streaming
   self.free_transport = self.operators.free_transport
   ```

5. Update cross-section rebinding to replace/reprepare only coefficient-dependent
   leaves and dependent caches.
6. Preserve object identity where a posed system must share the same leaves.
7. Add a test that `streaming` has no `Sigma_t` dependency and
   `free_transport == streaming + collision`.

**Exit gate:** no production attribute named `L` refers to `L+C`.

### Phase 6 — migrate non-SN consumers

1. **Homogeneous solver**
   - consume `FissionEnergyTransferOperator` directly;
   - retain explicit dense materialization as a solver strategy;
   - preserve exact dominant-eigenpair results.
2. **Diffusion solver**
   - use the scalar energy operator plus explicit bulk/full-field lift;
   - do not import `SNReactionRealizer`;
   - preserve the `(A,F)` generalized pencil before the exact inverse posing.
3. **Radial-characteristic System B**
   - consume the shared fission/scattering reduced factors;
   - retain SN-owned fold/reconstruction maps;
   - verify that outer fission emission and within-group scattering remain on
     their correct distinct seams.
4. Update docs and type annotations as each consumer migrates.

**Exit gate:** no non-SN consumer depends on an SN carrier dispatch hidden inside
the old common facade.

### Phase 7 — remove runtime carrier polymorphism

1. Convert old `FissionOperator` and `ScatteringOperator` into deprecated
   compatibility constructors, if downstream compatibility requires one release.
2. Emit actionable deprecation messages naming the replacement based on the old
   constructor arguments.
3. Ban new imports in architecture tests.
4. Remove `singledispatchmethod` arms after all internal consumers migrate.
5. Remove optional `full_field_space` from reduced/common operators.
6. Retire the facades and update `__all__`, API docs, examples, and type stubs.

**Exit gate:** every production operator has one public `apply` signature and
non-optional domain/codomain metadata.

### Phase 8 — formalize equation posings and eigen pencils

1. Add `GeneralizedEigenPencil`.
2. Make SN and diffusion expose criticality as `(A,F)` before selecting power
   iteration or direct generalized eigensolve.
3. Add `StationarySplitting(A, M, N)` with an invariant check/testing oracle
   `A == M - N`.
4. Rename:
   - `WithinGroupSystem.resolvent -> implicit_operator`;
   - `gains -> explicit_gains`.
5. Rename “regular splitting” documentation unless positivity is proven.
6. Make power iteration consume a `PowerMapPosing` or an adapter over a pencil.
7. Preserve the existing `EigenvalueSolver` protocol during migration; avoid
   simultaneous solver-engine redesign.

**Exit gate:** equation truth, algebraic posing, and solver execution are
separate inspectable objects.

### Phase 9 — add partitions, lazy block views, and graph schedules

1. Add restriction/injection-based `OperatorPartition`.
2. Add `PartitionedOperatorView`; keep `CoupledOperator` unchanged.
3. Implement structural dependency extraction.
4. Add SCC decomposition and condensation-DAG schedules.
5. Re-express current boundary sweep-cycle handling in terms of realized graph
   structure, retaining conservative hints as fast prechecks.
6. Demonstrate at least three equivalent SN views:
   - physical A/B carrier partition;
   - energy-group partition;
   - ordinate/bundle partition.
7. Replace public `System A` / `System B` labels with descriptive carrier names;
   keep temporary compatibility labels only where required.
8. Lower acyclic dependencies to triangular substitution/sweep schedules.
9. Keep cyclic SCCs explicit for iteration or coupled solves.
10. Rename the scattering data split so it does not claim solver foldability.

**Exit gate:** schedule selection derives from an explicit partition and graph,
not from hard-coded interpretations of “the SN system.”

### Phase 10 — preparation, fusion, and performance

1. Profile the composition-based realized operators.
2. Introduce prepared contexts only for measured costs:
   - harmonic tables;
   - group-transfer layouts;
   - cell/material dispatch indices;
   - fused source accumulation buffers.
3. Keep context separate from physical kernel and algebraic operator identity.
4. A fused implementation must retain:
   - the unfused factor tree;
   - forward equivalence;
   - transpose/adjoint equivalence;
   - rebinding invalidation rules;
   - provenance in diagnostics.
5. Compare solver iteration counts as well as wall time; a faster apply that
   changes a splitting or normalization is not a valid optimization.

**Exit gate:** optimized execution remains a faithful lowering of the public
operator algebra.

### Phase 11 — documentation and cleanup

1. Update `neutron_transport_grand_report_v3.md`:
   - revise the operator table to distinguish continuous/kernel symbols from
     realized method operators;
   - replace the §38 “semantic operator + prepare” wording with
     `kernel -> reduced operator -> method realizer -> prepared context`;
   - record that the former registry proposal was superseded by explicit
     method-owned realization;
   - add pencil/posing vocabulary;
   - add SCC/condensation scheduling.
2. Update:
   - `docs/theory/galerkin_projection.rst`;
   - scattering/fission theory docs;
   - SN solver and sweep docs;
   - diffusion and homogeneous docs;
   - API reference and examples.
3. Add an architecture decision record summarizing the fixed-domain operator
   rule.
4. Delete expired compatibility aliases and stale historical claims.

**Exit gate:** documentation describes the code that exists, and names match the
mathematics.

---

## 11. Verification matrix

### 11.1 Structural tests

- every class ending in `Operator` inherits the appropriate operator protocol;
- each public realized operator declares one non-optional domain and codomain;
- kernels have no `apply`;
- `transport` does not import method packages;
- method realizers live in method packages;
- no string registry is used for reaction realization;
- no runtime carrier dispatch remains on realized reaction operators;
- `SNDiscreteOperators.streaming` and `.free_transport` are distinct;
- `WithinGroupSystem.implicit_operator` is `M`, not `M^-1`;
- `CoupledOperator` and `PartitionedOperatorView` have non-overlapping roles.

### 11.2 Algebraic tests

- fission rank-one reconstruction:

  ```text
  F_G phi == chi * <nuSigma_f, phi>
  ```

- SN fission composition equivalence:

  ```text
  F_SN == j_bulk E0 F_G M0 pi_bulk
  ```

- SN scattering composition equivalence:

  ```text
  S_SN == j_bulk R Lambda M pi_bulk
  ```

- splitting invariant:

  ```text
  A == M - N
  ```

- pencil/power-map relation:

  ```text
  A psi = (1/k) F psi
  <=>
  A^-1 F psi = k psi
  ```

  when the inverse/solve is defined;

- partition reconstruction:

  ```text
  A == sum_ij J_i A_ij R_j
  ```

  under a complete compatible partition;

- SCC condensation graph is acyclic;
- triangular schedule results equal monolithic application/solve.

### 11.3 Adjoint tests

For every new factor and composition:

```text
<y, A x>_codomain == <A^H y, x>_domain
```

Test with the actual space metrics, including:

- quadrature weights;
- source/flux role spaces;
- boundary trace metric;
- direct-sum coupled metric;
- asymmetric group-transfer matrices.

### 11.4 Numerical regression tests

At minimum:

- isotropic and anisotropic SN fixed-source problems;
- multigroup SN eigenvalue problems;
- 1-D slab, sphere, cylinder, and 2-D Cartesian;
- reflective, periodic, vacuum, and prescribed-inflow boundaries;
- radial-characteristic carrying and seedless cases;
- diffusion versus homogeneous infinite-medium identities;
- existing generalized/direct eigenvalue crosschecks;
- cross-section rebind/depletion behavior;
- all current bit-identical tests where reduction order is intentionally part of
  the contract.

Where exact bit identity is lost solely because a mathematically equivalent
composition changes reduction order, require an explicit review. Do not silently
weaken tolerances.

### 11.5 Performance tests

Measure:

- apply time for reduced and realized reaction operators;
- allocation count;
- cache preparation time;
- sweep time;
- total inner iterations;
- total eigenvalue iterations;
- peak memory;
- block-view construction without dense materialization.

Performance gates should compare full solver work, not only isolated kernels.

---

## 12. Risks and controls

### Risk: too many conceptual layers

**Control:** kernels are introduced only for representation-changing physics.
Collision multiplication remains a direct reusable operator. Thin wrappers must
earn their existence through domain/codomain truth, invariants, or realization.

### Risk: kernel duplicates `MaterialXSField`

**Control:** implement kernels as semantic zero-copy views first. Cross-section
storage remains single-sourced in `MaterialXSField`.

### Risk: Python generics cannot express method-associated types cleanly

**Control:** keep common realizer protocols structural and minimal; give each
method concrete typed methods. Do not erase concrete types to satisfy an
over-general protocol.

### Risk: composition creates allocation/performance regressions

**Control:** use the composition as algebra-of-record and equivalence oracle;
permit later fusion through a prepared context.

### Risk: fission/scattering rebinding sees stale caches

**Control:** document live-view versus snapshot semantics and test mutation.
Prepared contexts must declare invalidation/rebuild keys.

### Risk: source normalization moves accidentally

**Control:** name the analysis and synthesis maps and pin normalization tests at
each factor boundary.

### Risk: `(n,2n)` convention changes physics during rename

**Control:** separate the leaf without changing the initial posing. Only change
loss-side versus gain-side placement in a dedicated, balance-gated posing
change.

### Risk: “folding” diagonal self-scatter changes the fixed point

**Control:** rename the data decomposition before exposing any implicit-fold
strategy. Require anisotropic regression cases and the exact splitting
invariant.

### Risk: block graphs confuse numerical nonzero with structural nonzero

**Control:** define structural support conservatively. Optional numerical
thresholding is a separate approximation policy and must never silently alter
the exact dependency graph.

### Risk: a flag-day solver rewrite obscures regressions

**Control:** preserve the existing driver protocols while carving operators;
formalize pencils and schedules only after operator realization is stable.

---

## 13. Decisions fixed by this plan

The following should not be reopened during routine implementation unless a
counterexample is found:

1. `Operator` means one fixed mathematical arrow.
2. The common fission/scattering object is a non-callable kernel, not a
   multi-carrier operator.
3. The method owns realization.
4. Realizers are passed/constructed explicitly; no global method-name registry.
5. SN reaction realization is expressed as analysis-transfer-synthesis
   composition.
6. Shared physical kernels and reduced transfer operators stay under
   `transport`.
7. Quadrature/frame-bearing reaction operators live under `sn`.
8. `(n,2n)` remains a separately named leaf even when a posing sums it with
   scattering.
9. `SNSolver.L` may not mean `L+C`.
10. `InvertibleOperator` is too generic for the SN free-transport sweep object.
11. `resolvent` may not name `M` when `M` has not been inverted.
12. “Regular splitting” requires mathematical positivity evidence.
13. `CoupledOperator` remains the explicit direct-sum block operator.
14. `PartitionedOperatorView` is the lazy block view of a monolithic operator.
15. Foldability is relative to a partition, schedule, implicit set, and solver
    capability.
16. SCC condensation is the basis of cyclic/acyclic scheduling.
17. The generalized eigen pencil `(A,F)` is primary; `A^-1 F` is a selected
    solver posing.
18. Algebraic compositions remain equivalence oracles after performance fusion.
19. Equation symbols, boundary symbols, physical carrier systems, and partition
    members must not share ambiguous letter-based public names.

---

## 14. Deliberate implementation choices to validate in Phase 0

These are narrow engineering choices, not unresolved architecture:

1. Whether `Kernel` or `Law` reads best in the public API. This plan recommends
   `Kernel` because fission and scattering are integral transfer kernels;
   `BoundaryTraceLaw` remains appropriate for affine boundary relations.
2. Whether the first additive carve creates `transport/kernels/` immediately or
   temporarily defines the kernels beside the old operators. The target is a
   kernel package; choose commit shape based on import-cycle risk.
3. Whether `SNFissionOperator` is a named subclass/wrapper or a type-preserving
   factory returning `OperatorProduct`. This plan recommends a thin named wrapper
   if it adds reliable domain/codomain and provenance; otherwise return the
   composition from the realizer and use a descriptive type alias.
4. Deprecation length for the two old facades. Internal callers should migrate
   in this campaign; external compatibility may justify one release.
5. The exact type names for scalar source spaces, which must follow the existing
   field/space hierarchy rather than invent near-duplicates.

None of these choices changes the mathematical layering.

---

## 15. Recommended first implementation session

The next session should implement only Phases 0–2:

1. produce the consumer inventory and baseline;
2. add `FissionKernel` and `LegendreScatteringKernel`;
3. extract `FissionEnergyTransferOperator`;
4. rename/add `LegendreMomentTransferOperator`;
5. make the old facades delegate to the new reduced operators;
6. run the focused fission, scattering, homogeneous, diffusion, SN, transpose,
   and type-check gates;
7. stop before adding `SNReactionRealizer`.

That boundary gives the next review a clean question:

> Are the physical kernels and reduced fixed-domain operators faithful and
> complete before ORPHEUS commits to their SN realization API?

Once that answer is yes, Phases 3–5 can establish the SN standard that later
transport methods follow.

---

## 16. Source anchors reviewed for this plan

These are navigation anchors for the implementation session, not claims that
their current names are authoritative:

| Area | Current source | Reason it matters |
|---|---|---|
| Grand architecture | `.claude/plans/neutron_transport_grand_report_v3.md` | common mathematical language, operator notation, boundary realizer, operator/context distinction |
| Fission facade | `orpheus/transport/operators/fission.py` | rank-one kernel plus multi-carrier dispatch |
| Scattering facade | `orpheus/transport/operators/scattering.py` | moment transfer, frame, quadrature, n2n, dispatch, and data split currently combined |
| Reduced isotropic transfer | `orpheus/transport/operators/isotropic_scattering.py` | existing method-neutral energy operators to reuse |
| Shared multiplication | `orpheus/transport/operators/multiplication_operator.py` | example of physics already represented by a faithful common operator |
| SN streaming/free transport | `orpheus/sn/operators/streaming.py` | sigma-free `StreamingOperator` and current generic `InvertibleOperator` |
| SN operator construction | `orpheus/sn/coupled_system.py` | one `L+C` builder and current stationary splitting record |
| SN solver bindings | `orpheus/sn/solver.py` | current `self.L = L+C`, reaction facade construction, and consumers |
| Explicit coupled blocks | `orpheus/numerics/coupled_system.py` | existing `CoupledOperator`, distinct from the proposed lazy block view |
| Loss realization/schedule | `orpheus/sn/loss_representation/` | current assembly, sweep graph, and schedule concepts to preserve/refine |
| Boundary descriptor | `orpheus/geometry/boundary/_base.py` | existing non-callable method-neutral law |
| Boundary realization seam | `orpheus/geometry/boundary/_realizer.py` | explicit method-owned realizer and the reason the registry was dissolved |
| SN method space | `orpheus/sn/mesh/method_space.py` | natural input to the proposed reaction realizer |
| Harmonic frame | `orpheus/transport/frames/harmonic_frame.py` and `orpheus/numerics/frame.py` | existing analysis/reconstruction factors |
| Diffusion consumer | `orpheus/diffusion/solver.py` | scalar/full-field fission realization and exact inverse posing |
| Homogeneous consumer | `orpheus/homogeneous/solver.py` | natural direct consumer of the reduced energy operator |
| Radial-characteristic coupling | `orpheus/sn/operators/radial_characteristic.py` | shared reduced kernels combined with SN-specific reconstruction/folding |
| Cross-method frame warning | `docs/theory/galerkin_projection.rst` | already records why a quadrature-bearing scattering facade cannot remain method-neutral |

The implementation should re-run the inventory rather than relying on these
anchors alone; line numbers and consumers will evolve.
