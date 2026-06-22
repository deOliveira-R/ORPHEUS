---
name: phase-g-four-operator-architecture-reconciliation
description: Reconciles the user's "four-operator end state" (L, C, F, S) — verbatim from neutron_transport_grand_report_v3.md — against the explorer audit's seven-operator proposal for Issue #196 manifestation #7. Pins every architectural decision to a section of the grand report, maps each of the explorer's seven proposed types onto where it actually lives in the L/C/F/S algebra, and rewrites the 5-line collapse targets and migration path under the four-operator target.
metadata:
  type: project
  branch: refactor/sn-operator-algebra
  spec_source: .claude/plans/neutron_transport_grand_report_v3.md (7396 lines, v3 revision)
  audit_source: .claude/agent-memory/explorer/issue_196_sn_operator_architecture_audit.md
  authored: 2026-05-12
---

# Phase G — Four-operator architecture reconciliation

**Why this memo exists.** The user, on reading the explorer audit
(`issue_196_sn_operator_architecture_audit.md`), corrected the
direction: the audit proposed seven new operator types
(`SNCellOperator`, `SNSweepOperator`, `AngularRedistribution`,
`PoleFaceAnchor`, `AngularIntegration`, `PackedStructuredAdapter`,
`SNBoundaryFaceTraceOperator`) to unify the SI-sweep / Krylov-matvec
twin path. The user's response:

> "You mentioned a lot of operators, but in the end we should have
> **L (streaming), C (collision), F (fission) and S (scattering)**.
> The end target should be that simple, and under the correct
> application, those operators look like math, and everything else
> using them looks like math as well (adjoint for example). This is
> why code elegance matters! It literally prevents bugs by
> construction."

The grand report
(`/Users/rodrigo/git/nuclear/ORPHEUS/.claude/plans/neutron_transport_grand_report_v3.md`)
is the spec. The user's four-operator target is its `§1 Notation
and operator convention` table (lines 25-40), `§30 The highest-level
architecture in one expression` (lines 5622-5649), and `§43 Closing
implementation principle` (lines 7375-7396). The audit must be
reconciled with the spec, not the reverse.

---

## §1. Faithful extraction from the grand report

### §1.1 The four operators — definitions, capabilities, parameters

The spec defines five operators in `§1` (grand report lines 25-40),
plus a boundary operator. The four the user cares about are L, C, S,
F. T is the time-mass operator (only enters in time-dependent / α
problems); B is the boundary operator (always present, composed as
part of L's effective discretisation per §16A).

**L — Streaming** (`StreamingOperator`)

- Math (continuous): `L = Ω · ∇` — the directional derivative
  along the streaming vector field. In tensor-product form (grand
  report §15.1, lines 2032-2046):
  ```
  L = D_x ⊗ Ω_x ⊗ I_g + D_y ⊗ Ω_y ⊗ I_g + D_z ⊗ Ω_z ⊗ I_g
  ```
- Math (discrete, SN): per ordinate it is a directed graph operator
  whose inversion is the sweep. Grand report §15A.2 (lines 2137-2171)
  names this structure as `UpwindTraceComplex` +
  `SweepDependencyGraph` — the local-to-global structure of L in SN.
- Codomain hierarchy: grand report §3.3 (lines 297-332) puts
  `StreamingOperator` under `DifferentialOperator`, which is a
  `LinearOperator`.
- Capabilities (§5.6 suffix rules at line 627; §32.6 dunder set at
  lines 5949-5982): `apply`, **`solve`** (i.e. L⁻¹ is the sweep),
  `H` (Hilbert adjoint), `T` (representation transpose), composition
  via `@`, tensor product via `&`.
- Parameters at construction: the **mesh** (provides cell complex,
  face normals, upwind classification — see grand report §7.2
  `SpatialMesh`, lines 828-846), and the **boundary resolution**
  (because the SN sweep cannot run without the inflow trace; the
  grand report's §16A.5 `DiscreteOrdinatesBoundaryRealizer` produces
  the realised BC operator that L's `.solve` consumes).
- Does L include the C term? **No.** Grand report §1 keeps them
  separate; `A_loss = L + C - S` (line 76) is the algebraic combiner.
  In current SN code the sweep wraps `(L + C)⁻¹` as a single operator
  because spatial and collision are inseparable per cell, but the
  algebraic interface still distinguishes them.

**C — Collision / total interaction** (`CollisionOperator`)

- Math: `C ψ = Σ_t · ψ` — pointwise multiplication. Grand report
  §1 line 35; §3.3 line 322 places `CollisionOperator` under
  `MultiplicationOperator`.
- Discrete: diagonal in `(cell, ordinate, group)` indexing.
- Capabilities: `apply`, `solve` (trivial inversion since diagonal),
  `H` (self-adjoint since real diagonal), composition.
- Parameters: a `TotalCrossSection` field (per material, per group;
  grand report §37 lines 6896-6926).
- Spec note: grand report `§15.6` Kernel-vs-Operator suffix rule
  (line 628-634) is explicit: `C` is an `Operator` because it maps
  fields to fields by multiplication; it is *not* a `Kernel` because
  it is not integrated against a measure.

**S — Scattering** (`ScatteringOperator`)

- Math (continuous): integral operator over Ω' with Legendre kernel
  (grand report §19, lines 3965-3994):
  ```
  S ψ(x, Ω, g) = Σ_{g'} ∫ Σ_s(x, Ω·Ω', g'→g) ψ(x, Ω', g') dΩ'
  ```
- Native discrete form for SN (grand report §9 lines 1221-1240,
  reaffirmed in §15.2 line 2049-2061):
  ```
  S_SN = R · Λ · M     where:
  M = harmonic moment projection = Y^* W
  R = harmonic reconstruction = Y
  Λ = Legendre/group scattering moments
  ```
  i.e. **scattering factors through harmonic moments**, not as a
  dense angular matrix. This is the spec's "exploit harmonic
  projection for SN" rule (§31 line 5673).
- Capabilities: `apply`, `H` (adjoint exists; grand report §5.6
  hierarchy line 648 places `ScatteringOperator` under
  `IntegralKernelOperator`), composition. No `solve` — S is not
  invertible standalone; `(L + C - S)⁻¹` is.
- Parameters: a `LegendreMomentScatteringData` field
  (§19 lines 3975-3993).

**F — Fission** (`FissionOperator`)

- Math: separable production (grand report §20.1 lines 4058-4078):
  ```
  F ψ = χ(g) · Σ_{g'} νΣ_f(g') · φ(g')
  ```
  Equivalently `F_g = χ ⊗ (νΣ_f)^T` — **rank-1 in energy**.
- Capabilities: `apply`, `H` (the adjoint flips the chi-vs-νΣ_f
  roles — `F.H = (νΣ_f) ⊗ χ^T`, useful for response/sensitivity per
  §25.2 lines 4622-4644), composition. The same physics data also
  has alternate representations as `BranchingFissionKernel` (MC
  side, §20.2) and `FissionGeneratingFunction` (§20.3) — see
  grand report §20 line 4106:
  ```python
  fission_data.as_linear_operator(space)   # → F
  fission_data.as_branching_kernel(space)  # → MC
  fission_data.as_generating_function(space)
  ```
- Parameters: `FissionData(nu, sigma_f, chi)` (line 4108).

**B — Boundary** (grand report §16, §16A, §16B)

The user said: "both the Krylov and the sweep are application of the
Operators including Boundary Operators." Grand report §16A (lines
2699-2914) is explicit: B decomposes as

```
B = R_boundary @ G_boundary     (response @ geometry)
```

(line 2718) and the **affine boundary trace law** is

```
γ_- ψ = R · G · γ_+ ψ + q_boundary    (line 2733)
```

So B is **not** a separate top-level operator joining {L, C, F, S}.
It is composed with L: every realised SN sweep takes a
`BoundaryRealizer` instance (§16A.4 lines 2862-2913, §16A.5 lines
2917-2961) that converts the abstract `BoundaryTraceLaw` into an SN
operator (e.g. `ZeroInflowOrdinateConstraint` for vacuum). The B
operator is owned by **L** in the same sense that the upwind
classification is owned by L — L cannot apply or solve without it.

### §1.2 Composition rules — the dunders

Grand report §6.3 (lines 714-737), reaffirmed in §32.6 (lines
5946-5982):

```python
A + B           # operator sum             →  OperatorSum
A - B           # operator difference      →  A + (-1)*B
α * A           # scalar scaling           →  ScaledOperator
A @ B           # composition              →  OperatorProduct
A @ ψ           # application              →  A.apply(ψ)
A & B           # tensor product           →  TensorProductOperator
A.H             # Hilbert adjoint
A.T             # representation transpose
A.inverse()     # the inverse operator (when CAP_SOLVE)
```

The fixed-source equation (`§1` line 56):
```
(L + C - S - F) ψ = q
```
is read in code as
```python
A_prompt = L + C - S - F                  # OperatorSum tree
psi = A_prompt.inverse().apply(q)
```
or equivalently `solve(A_prompt, q)` via a registered linear solver.

`A_loss = L + C - S` is defined at §1 line 76:
```python
A_loss = L + C - S
K      = A_loss.inverse() @ F             # multiplication operator
```

The four-operator algebra is closed under `+, -, *, @, .H, .T,
.inverse()` and consumes them to produce every other operator the
solver family needs. The user's standard "looks like math" is
materially the contract that `A_loss = L + C - S` returns a real
`LinearOperator` with `.apply`, `.H`, `.inverse()`, etc., not a
hand-coded function.

### §1.3 Solver expressions in the algebra

**Fixed source** (grand report §29.3 lines 5527-5544):
```python
A = L + C - S - F
psi = GMRESSolver().solve(FixedSourceProblem(A, q, V.boundary_resolution))
```

**Source iteration** (grand report §24 lines 4527-4581):
```python
A_loss = L + C - S
psi = SourceIteration(loss=A_loss, source=q + F @ psi_prev).solve()
```
where `SourceIteration` decomposes the LHS by sweep:
```python
for it in range(maxiter):
    Q = (S @ psi) + (F @ psi) / k + q_ext     # ← S @, F @ are mathy
    psi = (L + C).inverse().apply(Q)          # ← the sweep
```

**k-eigenvalue / power iteration** (grand report §21.2 lines
4162-4188; §29.2 lines 5503-5525):
```python
A_loss = L + C - S
problem = CriticalityProblem(loss=A_loss, production=F,
                             boundary=V.boundary_resolution)
k, psi = PowerIteration(tol=1e-8).solve(problem)
```
Inside `PowerIteration.solve` (§36.1 lines 6705-6743):
```python
for it in range(maxiter):
    fission_source = F @ psi
    psi_next       = A_loss.inverse().apply(fission_source)
    k_next         = <psi_dagger, F @ psi> / <psi_dagger, A_loss @ psi>
    psi            = normalize(psi_next)
```

**Adjoint** (grand report §25.2 lines 4620-4644):
```python
psi_dagger = (A_loss.H).solve(R.as_source())
```
i.e. **adjoint flux is literally `(L + C - S).H` applied to the
response source**. The "adjoint reads like math" is exactly
`A_loss.H.solve(r)`. No separate `solve_sn_adjoint` machinery is
needed at the algebra level — it is `.H` on the composed operator,
delegated to the same Krylov / sweep that solves the forward.

**α-eigenvalue** (grand report §22, §29.7 lines 5599-5618):
```python
A_prompt = L + C - S - F
T = TimeMassOperator(xs.velocity)
problem = AlphaEigenproblem(time_mass=T, prompt_operator=A_prompt)
spectrum = ShiftInvertAlphaSolver(shifts=[0., -1., -10.]).solve(problem)
```

**Compute k from the converged flux** (Issue #169 wiring):
```python
k = (psi_dagger @ (F @ psi)) / (psi_dagger @ (A_loss @ psi))
```
or, since the production is the operator-applied production rate
already wired per Issue #169 (`compute_group_production_rate`):
```python
k = total_production_rate(F @ psi) / total_loss_rate(A_loss @ psi)
```

### §1.4 Boundary conditions are first-class but composed, not co-equal

Grand report §16A.1 (lines 2738-2776) layers BCs in four levels:

1. `BoundaryDeclaration` — semantic declaration at geometry stage
   (e.g. `("outer", "vacuum")`).
2. `BoundaryPatch` — geometric resolution at mesh stage (face IDs).
3. `BoundaryTraceLaw` — the affine law `γ_- = R·G·γ_+ + q`
   (§16A.2 line 2792).
4. `BoundaryRealizer` — method-specific operator
   (`DiscreteOrdinatesBoundaryRealizer`, §16A.4 line 2882) →
   produces `MethodBoundaryOperator` = a realised SN
   `LinearOperator` (e.g. `ZeroInflowOrdinateConstraint`).

The realised SN boundary operator is then consumed **inside L's
`.solve` and `.apply`**, not as a separate top-level operator. So the
user's algebra reads `(L + C - S - F)·ψ = q` at the top level, and L
internally composes the realised B inside its directed `.solve`. The
algebra is still four operators at the visible API; B is an
implementation detail of L's `.solve`/`.apply`. This is the spec's
"operator versus context" distinction (§38 lines 6941-6985).

### §1.5 Cross-solver scope — same names, method-specific realisations

Grand report §29 (lines 5481-5618) gives the algebra for each method:

- **SN** (§29.2): `A_loss = L + C - S`; `CriticalityProblem(A_loss,
  F, boundary)`; `PowerIteration().solve(...)`.
- **PN** (§29.3): `A = L + C - S - F`; `GMRESSolver().solve(...)`.
  Same four operators, different discretisations.
- **MoC** (§29.4): the sweep is replaced by `CharacteristicSweep`,
  but `S`, `F`, `q` enter via `SourceIteration(sweep, scatter, fission)`.
- **CP** (§29.5): operators are `CollisionGreenOperator`,
  `CollisionProbabilityKernel`; the algebra `(I - P·S - P·F/k) ψ = q`
  uses the same `+`, `-`, `@` composition.
- **MC** (§29.6): `S, F` realised as Markov / branching kernels,
  consumed by `NeutronBranchingProcess`.

So L, C, F, S are **method-specific operator instances** that
inherit a common `LinearOperator` interface (grand report §3.1,
§5.7, §38) and compose under the same algebra. The names L, S, F
are the SAME across SN, MoC, PN, CP, MC; their `.apply` does
different math because the method-space (`DiscreteOrdinatesPhaseSpace`,
`CharacteristicTrackSpace`, etc.) supplies different per-cell /
per-track / per-region representations. Grand report §31 rule 4
(line 5661): "Use specific method-space names" but operator names
stay the same.

### §1.6 The end state in one expression — grand report §30 line 5626

```python
V = X * Omega * G * Xi
A_prompt = L + C - S - F
fixed_source = FixedSourceProblem(A_prompt, q, boundary=B)
k_problem = CriticalityProblem(loss=L + C - S, production=F)
alpha_problem = AlphaEigenproblem(T, A_prompt)
```

This is the elegance target. Every solver path is a 1-line
composition of L, C, S, F (and T for time / α). The user's
correction — "the end target should be that simple" — is this
expression.

---

## §2. Reconciliation with the explorer audit's seven types

For each of the explorer's seven proposed operator types, I map it
onto where it lives under the grand report's four-operator vision.
The verdict for each is one of: **(IN L)** — it is a detail of L's
discretised `.solve`/`.apply`; **(IN S/F)** — it is the SN
representation of S or F; **(NOT-AN-OPERATOR)** — it is an algorithm,
adapter, or precompute, not a `LinearOperator` co-equal to L/C/S/F;
**(SEPARATE-BUT-COMPOSED)** — it is its own `LinearOperator` but
consumed *inside* a four-operator composition rather than as a
top-level peer of {L, C, S, F}.

### §2.1 `SNCellOperator` — VERDICT: **(IN L)** + **(IN C)**

The explorer audit (§5.1 row 1) proposed this as an 8-line affine
per-cell operator. Under the grand report, the per-cell update is
**not** an operator at the public API — it is one local stencil of
the global `(L + C)` operator's discretised inverse. Grand report
§32.8 (lines 6024-6048) explicitly separates the **operator** from
its **representation**: `L.represent("sweep")` produces a sweep
representation whose per-cell update is internal.

The grand report does name the entities the explorer's
`SNCellOperator` captures:
- The local update is a "stalk" in the sweep dependency graph
  (§15A.2 line 2143: `UpwindTraceComplex`, `SweepDependencyGraph`)
  — these are **structural metadata of L**, not separate operators.
- The WDD / symmetric closure choice is a `Representation` of L
  (§32.8 line 6044: `A.represent("sweep")`).

So `SNCellOperator` becomes:

```python
class SweepRepresentation(Representation):
    """L's discretised solve representation, parameterised on closure."""
    closure: Literal["WDD", "symmetric"]
    cell_update_fn: Callable     # WDD recurrence or symmetric
```

This is a **closure strategy** of L, not a peer of L.

### §2.2 `SNSweepOperator` — VERDICT: **(IN L)** — this is literally L⁻¹

The explorer's `SNSweepOperator` is, mathematically, `(L + C)⁻¹` on
the within-group source. Grand report §29.2 lines 5524 makes this
explicit: `PowerIteration(tol=1e-8).solve(problem)` consumes
`A_loss = L + C - S` and internally calls `A_loss.inverse()` —
which, for SN, IS the sweep.

So `SNSweepOperator` is NOT a new top-level operator. It is the
**result of `(L + C).inverse()`** under the sweep representation.
Code-wise:

```python
class StreamingOperator(LinearOperator):  # this is L
    capabilities = {CAP_APPLY, CAP_SOLVE, CAP_TRANSPOSE}

    def apply(self, psi):       # Lψ via WDD closure
    def solve(self, q, bc):     # the sweep (= L⁻¹·q)
```

There is no separate `SNSweepOperator`. The sweep is `L.solve`.
This is what the spec means by "elegant" — the user reads
`L.solve(q)` and knows: "this is the sweep, which is the inverse of
L applied to q." Names that imply algorithms (Sweep, Iteration,
Solver) live one layer above (in `numerics.iteration` /
`numerics.eigenvalue`, grand report §24); names of mathematical
objects (L, C, S, F) live at the operator layer.

### §2.3 `AngularRedistribution` — VERDICT: **(IN L)** — curvilinear streaming detail

The Morel-Montry angular recurrence in curvilinear geometry is a
**discretisation detail of L** in spherical / cylindrical mesh.
Grand report §15A.2 line 2160 names this as a `ReflectiveSweepCycle`
/ `CyclicSweepBlock` inside the `SweepDependencyGraph` (which is
metadata of L's discretisation). The Carlson seed and the M-M
recurrence are the SN curvilinear realisation of the streaming
operator's directional derivative `Ω · ∇` in a coordinate system
with non-Cartesian Christoffel-like terms — they are L's spherical
implementation, not a new operator.

Under the grand report, this becomes:

```python
class StreamingOperator(LinearOperator):
    def __init__(self, mesh, sigma_t, *, angular_recurrence=None):
        # angular_recurrence is the M-M closure for curvilinear geometry;
        # None for slab where it doesn't apply.
        ...
```

`angular_recurrence` is a strategy on L's solve, not a peer operator.

### §2.4 `PoleFaceAnchor` — VERDICT: **(IN L)** — boundary detail of L at r=0

The pole face anchor convention (whether `psi_face_in` at i=0 is
zero or the cell-centre value) is a **boundary condition detail at
the pole face**. Grand report §16B (lines 3328+) is the right
framework: the pole r=0 in spherical geometry is a **symmetry
boundary**, not a physical vacuum or albedo. The full domain
is `r ∈ [-R, R]` quotiented by `r → -r`. The pole face is a
`SymmetryBoundaryPatch` (§16B.1 line 3372), not a `VacuumBoundary`.

Under the grand report, the pole-face anchor becomes a
`SymmetryBoundary` instance composed into L's BC resolution:

```python
geometry = GeometrySpec(
    ...,
    boundary_declarations=(
        BoundaryDeclaration("outer", "vacuum"),
        SymmetryBoundaryDeclaration("pole", generator="reflect_r"),
    ),
)
# DiscreteOrdinatesBoundaryRealizer realises the symmetry boundary
# into a SNOrdinateSymmetryMap that handles the r=0 anchor.
```

The explorer's `PoleFaceAnchor` becomes one realisation of the
"pole=symmetry-boundary" pattern: an `SNSymmetryBoundary` operator,
which is a `BoundaryRealizer` output, consumed inside L.

### §2.5 `AngularIntegration` — VERDICT: **(NOT-AN-OPERATOR / IS the inner product)**

`φ = ∫ ψ dΩ = Σ_n w_n ψ_n` is **not a top-level operator** in the
grand report; it is the angular space's inner product. Grand
report §9 (lines 1170-1192) says the angular space is `L^2(S^2, dΩ)
→ ℓ^2(W)` where `W` is the weight matrix — this is the **inner
product structure** of `DiscreteAngularSpace` (§32.4 lines
5909-5915). The angular-to-scalar reduction is:

```python
phi = V_angular.inner(one_field, psi)   # = Σ w_n · ψ_n
```

Or equivalently a `MomentProjection` (grand report §5.6 line 657,
§9 line 1238): `phi = M @ psi` where `M = HarmonicMomentProjection`
at `L=0`. This is composition with the moment projector, which IS
in the algebra — but it is a `ProjectionOperator`, not a co-equal
peer of L/C/F/S.

So `AngularIntegration` becomes:

```python
M0 = HarmonicMomentProjection(angular_space, ell_max=0)  # scalar moment
phi = M0 @ psi
```

This composes with S's structure (`S = R @ Λ @ M`, §9 line 1240),
so the operator algebra naturally provides scalar flux without a
new top-level type.

### §2.6 `PackedStructuredAdapter` — VERDICT: **(NOT-AN-OPERATOR / is a `Representation` choice)**

The packed `(n_unknowns,)` ↔ structured `(N, nx, ny, ng)` layout
conversion is a **storage representation choice**, not a mathematical
operator. Grand report §32.8 (lines 6024-6048) explicitly:

```python
class Representation(ABC):
    def materialize(self): ...

class DenseRepresentation(Representation): ...
class SparseRepresentation(Representation): ...
class MatrixFreeRepresentation(Representation): ...
class SweepRepresentation(Representation): ...
class TensorTrainRepresentation(Representation): ...
```

The packed-vs-structured distinction is two `Representation`s of the
same field. The grand report's rule (§32.5 line 5938): "Fields
should not be allowed to silently combine if their spaces are
incompatible." So the adapter is a **field representation
converter**, internal to L's matrix-free SciPy `LinearOperator`
interface. It does not become a `LinearOperator` in its own right.

If the SciPy interface forces a packed 1D layout, then the structured
↔ packed conversion lives inside `L.represent("scipy_matvec")` —
a representation choice of L — not as a separate type.

### §2.7 `SNBoundaryFaceTraceOperator` — VERDICT: **(SEPARATE-BUT-COMPOSED)** — this IS B at the SN realisation

This one DOES exist as a separate `LinearOperator` in the grand
report. It is precisely the **method-specific realised boundary
operator** from `DiscreteOrdinatesBoundaryRealizer.realize(law,
space)` (§16A.5 line 2949): the output `SNBoundaryOperator`. It is
a `LinearOperator` because B = `R_boundary @ G_boundary + q_boundary`
(§16A.1 line 2718) and `R`, `G` are `LinearOperator`s, `q_boundary`
is a `Source`.

But under the grand report it is NOT a top-level operator coequal to
{L, C, F, S}. It is the **realised boundary operator owned by L** —
specifically, it is the operator L's `.solve` consumes at the
inflow trace. Per the grand report:

```python
V_sn = DiscreteOrdinatesPhaseSpace.from_mesh(...)
B = V_sn.resolve_boundaries()       # ← method-specific realised BC
L = StreamingOperator().on(V_sn)    # ← takes V_sn implicitly, which carries B
```

So `SNBoundaryFaceTraceOperator` is one of `B`'s realisations, and B
is composed into L. The user's "BCs are first-class" is satisfied by
B being a `LinearOperator` with `.apply`, `.H` etc.; the user's
"the end target should be only L, C, F, S" is satisfied by B being
consumed inside L, not co-equal to it.

### §2.8 Reconciliation summary table

| Explorer's proposed type           | Grand-report verdict                                | Where it lives under §1-§31                                   |
|------------------------------------|-----------------------------------------------------|----------------------------------------------------------------|
| `SNCellOperator`                   | (IN L+C) — `SweepRepresentation` strategy           | §32.8 Representation; §15A.2 SweepDependencyGraph stalk       |
| `SNSweepOperator`                  | (IN L) — this IS `L.solve` (and `(L+C).solve`)      | §5.7 SweepOperator under LinearOperator; §29.2 power-iter use  |
| `AngularRedistribution`            | (IN L) — curvilinear streaming closure detail       | §15A.2 ReflectiveSweepCycle; §22 sign convention metadata     |
| `PoleFaceAnchor`                   | (IN L via B) — `SymmetryBoundary` realisation       | §16B.1 SymmetryBoundaryPatch; §16B.8 SNOrdinateSymmetryMap    |
| `AngularIntegration`               | (NOT-AN-OPERATOR) — angular space inner product     | §32.4 DiscreteAngularSpace.inner; §9 MomentProjection at ℓ=0  |
| `PackedStructuredAdapter`          | (NOT-AN-OPERATOR) — field `Representation` choice   | §32.8 Representation hierarchy (SweepRepr, MatrixFreeRepr)   |
| `SNBoundaryFaceTraceOperator`      | (SEPARATE-BUT-COMPOSED) — this IS B in SN           | §16A.4 DiscreteOrdinatesBoundaryRealizer output → SNBoundaryOp |

**Net count**: seven explorer types collapse to **four user-facing
operators (L, C, S, F)** + **one composed boundary operator (B)
owned by L**, with the rest absorbed as representations, strategies,
inner-product structure, or projection composition. The end state
agrees with the user's correction.

---

## §3. Concrete 5-line collapse targets

These show what each load-bearing path looks like at the END of the
migration. They are the literal "reads like math" expressions the
user demanded.

### §3.1 `_sweep_1d_spherical` — currently 199 lines

```python
def _sweep_1d_spherical(Q, sig_t, sn_mesh, psi_bc, Q_aniso=None):
    # The sweep IS L⁻¹. Period.
    L = StreamingOperator(sn_mesh, sig_t, bc=psi_bc)
    psi = L.solve(Q + (Q_aniso or 0))
    return psi, M0 @ psi                # M0 = scalar moment projection
```

Five lines. The for-loops over ordinates, the WDD recurrence, the
Carlson seed, the pole anchor — all live INSIDE
`StreamingOperator.solve` under the chosen sweep representation. The
spherical / cylindrical / slab dispatch is dispatch on
`sn_mesh.geometry` inside L.solve.

### §3.2 `_solve_fixed_source_si` — currently 60 lines

```python
def _solve_fixed_source_si(L, C, S, F, q_ext, tol, maxiter):
    A_loss = L + C - S
    return SourceIteration(A_loss, fission=F, source=q_ext).solve(
        tol=tol, maxiter=maxiter,
    )
```

Where `SourceIteration` (grand report §24 line 4540) iterates:

```python
class SourceIteration(Solver, key="source_iteration"):
    def solve(self, *, tol, maxiter):
        psi = self.A_loss.domain.zero()
        for it in range(maxiter):
            Q   = self.source + (self.fission or 0) @ psi
            psi = self.A_loss.inverse().apply(Q)   # = the sweep
            if converged: break
        return psi
```

Note: `A_loss.inverse()` returns a `LinearOperator` whose `.apply`
delegates to `(L + C - S).solve` — which for SN is implemented as
"sweep, then add scattering source explicitly per ordinate" because
S is not invertible alone. The exact implementation is hidden
inside `OperatorSum.inverse()` plus a `WithinGroupGaussSeidel`
strategy. The user-facing code is 3 lines.

### §3.3 `_solve_fixed_source_krylov` — currently 158 lines

```python
def _solve_fixed_source_krylov(L, C, S, F, q_ext, tol, maxiter):
    A = L + C - S - F                              # OperatorSum
    return PreconditionedGMRES(
        op=A, precond=(L + C).inverse(),           # sweep preconditioner
    ).solve(q_ext, tol=tol, maxiter=maxiter)
```

Four lines visible. SciPy's `gmres` is called inside
`PreconditionedGMRES.solve`; the packed/structured adapter is one
representation-choice deep, invisible at this layer.

### §3.4 `solve_sn_eigenvalue` (k-eigenvalue) — currently ~140 lines

```python
def solve_sn_eigenvalue(materials, mesh, quad, *, tol, maxiter):
    V = DiscreteOrdinatesPhaseSpace.from_mesh(mesh, quad, materials.groups)
    L, C, S, F = StreamingOperator.on(V), CollisionOperator.on(V, materials.total), \
                 ScatteringOperator.on(V, materials.scatter), \
                 FissionOperator.on(V, materials.fission)
    problem = CriticalityProblem(loss=L + C - S, production=F, boundary=V.boundary)
    return PowerIteration(tol=tol, maxiter=maxiter).solve(problem)
```

Six lines (the unpacking is one line cosmetically split).

### §3.5 `solve_sn_adjoint` — does not exist today; would be

```python
def solve_sn_adjoint(materials, mesh, quad, response, *, tol, maxiter):
    V = DiscreteOrdinatesPhaseSpace.from_mesh(mesh, quad, materials.groups)
    L, C, S, F = build_sn_operators(V, materials)
    A_loss = L + C - S
    return GMRESSolver(tol=tol).solve(A_loss.H, response.as_source())
```

Five lines. **The adjoint flux is `A_loss.H.solve(response)`.** The
`.H` on the OperatorSum tree propagates: `(L+C-S).H = L.H + C.H -
S.H`, where `C.H = C` (real diagonal), `S.H` swaps scatter `from→to`
to `to→from`, and `L.H` reverses the streaming direction (= adjoint
sweep going from outflow to inflow). Each operator's `.H` is a
single decision; the algebra propagates it; no separate
`solve_sn_adjoint` machinery is built.

This is the user's elegance proof. If the algebra is correctly
designed, the adjoint solver requires zero new code beyond `.H`
implementations on each leaf operator.

### §3.6 `compute_keff` (per Issue #169) — currently a single function

Already mostly aligned per Issue #169. Under the four-operator
vision:

```python
def compute_keff(L, C, S, F, psi, psi_dagger=None):
    A_loss = L + C - S
    if psi_dagger is None:
        # If no adjoint provided, use psi as its own weight (forward k)
        return (psi @ (F @ psi)) / (psi @ (A_loss @ psi))
    # Adjoint-weighted Rayleigh quotient — exact eigenvalue estimator
    return (psi_dagger @ (F @ psi)) / (psi_dagger @ (A_loss @ psi))
```

The inner product `psi @ phi` is `V.inner(psi, phi)` from the
DiscreteAngularSpace; algebra-native. This IS the
Rayleigh quotient (grand report §36.1 line 6740, §21.2).

---

## §4. Migration path under the four-operator target

The explorer's 7-step migration assumed seven new operator types.
With the four-operator target the path collapses to **5 steps**.
Each step is independently committable; each preserves the 11
regression snapshots in `tests/sn/regression/snapshots/` (with
principled bit-identity breaks documented per `vv-principles`
§"Bit-identity vs principled-equivalence").

### Step 1 — Promote the existing `LinearOperator` algebra to be the visible API

**What.** Replace the four core ops (`SNStreamingOperator`,
`CollisionOperator` (currently inside L), `ScatteringOperator`,
`FissionOperator`) with the names `StreamingOperator`,
`CollisionOperator`, `ScatteringOperator`, `FissionOperator` per
grand report §1. Wire their construction through a
`DiscreteOrdinatesPhaseSpace` (§29.2 line 5508).

**Acceptance**: existing tests pass with renames only. No behaviour
change. **Reads like math after**: yes — `L = StreamingOperator.on(V)`,
`A_loss = L + C - S`.

**Closes**: nothing (preparatory). **Independently committable**: yes.

### Step 2 — Make L own the sweep (collapse `SNStreamingOperator.apply` and `.solve` to one closure)

**What.** Pick the WDD closure (per audit §7.1 recommendation,
backed by Bailey-Morel-Chang asymptotic-diffusion-limit) as the
single closure for both `L.apply` and `L.solve`. Implement this by
routing both `apply` and `solve` through ONE `SweepRepresentation`
that holds a single per-cell update strategy (`DiamondDifference`
WDD). The current symmetric closure in `apply` is removed (or made
selectable via `L.represent("symmetric")` as a deliberate, tested
alternative for Wave-E reconciliation).

**Acceptance**: (a) `test_sweep_vs_apply_consistency.py` (57 tests)
green at machine precision on non-flat ψ — manifestation #7
dissolves **by construction**. (b) Phase E sentinel
(`test_phase_e_trajectory_resolvent_flux_shape_crosscheck`) xpasses
and the strict marker can be removed. (c) The 11 regression
snapshots in `tests/sn/regression/snapshots/` are regenerated under
the unified closure — principled bit-identity break per
`vv-principles`, because:
  - the named intermediate is "WDD closure across both apply and
    solve",
  - the structurally-independent reference is the MMS L1 convergence
    test (`tests/sn/test_mms_*.py`) which is unchanged.

**Closes**: **Issue #196 manifestation #7** by construction.
**Independently committable**: yes (it's a single closure decision
guarded by the 57 consistency tests).

### Step 3 — Wire `SourceIteration` and `PreconditionedGMRES` as algebra consumers

**What.** Replace `_solve_fixed_source_si` and
`_solve_fixed_source_krylov` with the §3.2 / §3.3 versions above.
`SourceIteration` lives in `orpheus/numerics/iteration.py` (already
exists; grand report §24 confirms naming). `PreconditionedGMRES`
wraps `scipy.sparse.linalg.gmres` and accepts `op`, `precond` as
`LinearOperator`s.

**Acceptance**: full SN test suite green; MMS convergence rates
preserved; eigenvalue benchmarks match.

**Closes**: the "load-bearing paths don't compose operators" audit
finding (§4.2). **Independently committable**: yes (each call-site
is a separate commit).

### Step 4 — Promote BCs to first-class composition via `BoundaryRealizer`

**What.** Replace the per-ordinate / per-call-site `bc.apply(...)`
spaghetti with the grand-report §16A pipeline:

```python
geometry = GeometrySpec(boundary_declarations=(...))
V = DiscreteOrdinatesPhaseSpace.from_mesh(mesh, quad, groups,
                                           boundary_decl=geometry.boundary_declarations)
B = V.resolve_boundaries()                # SNBoundaryOperator composition
L = StreamingOperator(V, sigma_t, boundary=B)
```

The realised `B` is composed INSIDE L's solve and apply, ONCE per
solve, at the face-trace level (the §16A.5 contract: realised
SNBoundaryOperator takes `(N, ng)` outflow → `(N, ng)` inflow). The
pole-face anchor becomes a `SymmetryBoundary` per §16B.

**Acceptance**: BC tests pass; the explicit invariants from §16A.12
(line 3262) become executable assertions in the test suite.

**Closes**: the "BC realised but not woven in" audit finding (§4.3).
**Independently committable**: yes.

### Step 5 — Adjoint solver `solve_sn_adjoint` as `.H` on the operator expression

**What.** Implement `.H` on each leaf operator: `L.H` (adjoint
sweep — reverse direction), `C.H = C` (self-adjoint diagonal),
`S.H` (transpose group transfer), `F.H = ((νΣ_f) ⊗ χ^T)`. Then
`(L + C - S).H` is automatic from `OperatorSum.H` (already
implemented in `numerics/operator.py:_AdjointOperator`, line 415).

Add a single function `solve_sn_adjoint(...)` per §3.5. **Five
lines.** This is the elegance acceptance criterion: the adjoint
solver writes itself.

**Acceptance**: adjoint flux passes the bilinear-form identity
`<response, ψ_forward> = <ψ_adjoint, source>` to machine precision
on MMS cases.

**Closes**: the "adjoint mode unavailable" audit finding (§3.2 line
353-355, §4.2). **Independently committable**: yes.

### Net migration: 5 steps, not 7

The explorer's seven types collapse:
- Steps 1, 5, 6 of the explorer audit → my Step 1 (renames +
  algebra wiring as visible API).
- Steps 2, 3 of the explorer audit → my Step 2 (single closure for
  L; collapses to one cell-update strategy hidden inside L's sweep
  representation).
- Step 4 of the explorer audit → ALSO subsumed in my Step 2 (the
  inlined matvec WDD vs the factored sweep WDD become the same code
  path inside `L.solve` / `L.apply`).
- Step 5 of the explorer audit (PackedStructuredAdapter) → becomes a
  `Representation` choice inside L; not a separate migration step,
  it's a refactor inside Step 2.
- Step 6 of the explorer audit (iteration primitives) → my Step 3.
- Step 7 of the explorer audit (closure documentation) → integrated
  into my Step 2 acceptance criterion.

My Step 4 (BC composition) and Step 5 (adjoint via .H) are not in
the explorer audit at all — they are the spec's elegance acceptance
criteria.

**Step that closes Issue #196 by construction**: Step 2 (unified
closure for L.apply and L.solve).

---

## §5. Open questions where the grand report leaves ambiguity

### §5.1 Does L include σ_t·ψ, or is C separate?

Grand report §1 (lines 34-35) is explicit:
- `L = Ω · ∇` only — the differential / streaming part.
- `C = Σ_t · ψ` — separate multiplication.

But §15.1 (line 2032) writes `L = D_x ⊗ Ω_x ⊗ I_g + ...` which is
pure streaming, no σ_t. And `A_loss = L + C - S` in §1 line 76
combines them. So **the spec answer is: L is pure streaming; C is
separate.**

BUT, in the SN sweep implementation, `(L + C)⁻¹` is what the sweep
actually inverts (per-cell, the WDD balance is `streaming + σ_t·ψ =
source`, both terms together). The grand report addresses this in
§32.7 (lines 5985-6020) and §38 (lines 6941-6985): the **operator
algebra** treats them separately, the **operator context** (the
sweep representation) composes them for the per-cell solve. The
public-API decision:

- `L.apply(psi)` applies just the streaming term.
- `C.apply(psi)` applies just `σ_t · psi`.
- `(L + C).solve(q)` is the sweep — but `L.solve(q)` alone is
  not meaningful (you cannot solve streaming without collision).
  The algebra catches this: `L.solve` raises `MissingCapability`
  unless explicitly fused.

**Decision needed at Step 1**: do we expose `(L + C).solve(q)` as
the sweep, with `L` having only `apply` and `(L + C)` getting `solve`
via fusion? Or do we keep the current pattern where the
streaming-operator class owns the sweep and σ_t is a constructor
parameter? The grand report mildly favours the former (cleaner
algebra) but the latter (current code) is more performant. I
recommend the **former** for elegance — the user's standard — with
a `FusedStreamingCollisionOperator` as a performance representation
hidden inside `(L + C).inverse()`.

### §5.2 Is `S = R @ Λ @ M` constructed automatically, or by the SN scattering operator constructor?

Grand report §9 (line 1240) shows the explicit harmonic factoring:
```python
M = HarmonicMomentProjection(angular_space, L=xs.scatter.order)
R = HarmonicReconstruction(angular_space, L=xs.scatter.order)
S = R @ LegendreMomentScattering(xs.scatter) @ M
```

This is the **explicit factoring** — the user writes three lines to
build S. Versus the current code where `ScatteringOperator.from_solver_data`
hides this:

```python
S = ScatteringOperator.from_legendre_moments(xs.scatter).on(V_sn)
```

(§29.2 line 5518.) Both forms appear in the spec. The
constructor-form hides M, Λ, R inside S; the factoring form makes
them visible.

**Decision needed at Step 1**: do we expose M, Λ, R as separate
operators (more elegant but more verbose), or hide them in S
(current pattern, more concise)? I recommend **hide by default,
expose via `.factors`** so `S.factors` returns `(R, Λ, M)` for users
who want to compose at the moment level. This keeps `A_loss = L + C
- S` clean while allowing `S_factored = S.factors[0] @ S.factors[1]
@ S.factors[2]` for sensitivity / PN cross-method composition.

### §5.3 Does B compose into L, or into the OperatorSum?

The grand report §16A is ambiguous on this:
- §16A.5 line 2954: `SNBoundaryOperator(linear=R @ G, source=q)` is
  itself a `LinearOperator`. It could be composed into the sum:
  `(L + C - S - F + B) ψ = q`.
- §31 rule 7 (line 5667): "boundary conditions [are] operators on
  trace spaces" — they apply on the boundary, not on the bulk.
  Suggests they compose into L (which is the only operator that
  touches boundary traces).

Per the user's verbatim direction ("both the Krylov and the sweep
are application of the Operators including Boundary Operators"), the
BC is co-equal at the algebra level — but per the four-operator
endpoint, it must NOT be a fifth top-level operator. The
reconciliation is: **B is a `LinearOperator` (yes), composed inside
L (yes), and never seen as a top-level peer of L/C/S/F (yes)**.
Code-wise, `L = StreamingOperator(V, sigma_t, boundary=B)` — B is a
constructor parameter of L, not a top-level summand.

### §5.4 Where does the angular integration `φ = Σ w_n ψ_n` live?

Three possible homes:
1. As `V_angular.inner(one_field, psi)` — the inner product
   structure (cleanest).
2. As `M0 = HarmonicMomentProjection(..., ell_max=0)` and `phi = M0 @
   psi` — composition with the moment projector (most algebra-native).
3. As `AngularIntegrationOperator` — a named LinearOperator (current
   `_scalar_flux_from_angular`).

I recommend option 2: angular integration is the ℓ=0 component of S's
M projector, so it shares infrastructure with `S = R @ Λ @ M`. The
explorer's `AngularIntegrationOperator` is option 3 — verbose. The
grand report's §32.4 inner-product (option 1) is also fine but
splits the moment infrastructure.

### §5.5 Curvilinear-specific operators — where do they live?

The Carlson seed (`carlson_inward_sweep_from_source`), the M-M
recurrence, the Phase D pole-angular closure — are these:
- A `Representation` of L for spherical geometry? (cleanest)
- A separate `CurvilinearStreamingOperator` subclass of L?
- A strategy plug-in on L?

The grand report §9 (lines 1171-1259) names SN-specific structures
(`SweepGraph`, `UpwindStencil`, etc.) but does not split slab vs
curvilinear at the operator level. I recommend: **L is one class;
the spherical / cylindrical / slab dispatch happens inside L's
`solve` method based on `V.mesh.geometry`**. The Carlson seed and
M-M recurrence become methods on a `CurvilinearSweepStrategy` that
L uses internally — invisible to the user-facing algebra.

---

## §6. Self-checks against the user's elegance standard

The user's criterion: "those operators look like math, and everything
else using them looks like math as well." Test cases:

**Test 1: the fixed-source equation.**
- Spec form: `(L + C - S - F) ψ = q`.
- Code form: `A_prompt = L + C - S - F; psi = A_prompt.inverse() @ q`.
- Verdict: passes. The Python is the math, character for character.

**Test 2: the k-eigenvalue iteration.**
- Spec form: `K ψ = k ψ` where `K = A_loss⁻¹ F`.
- Code form: `K = (L + C - S).inverse() @ F; k, psi = ArnoldiEigenSolver().solve(K)`.
- Verdict: passes.

**Test 3: the adjoint flux.**
- Spec form: `A_loss^* ψ^† = R`.
- Code form: `psi_dagger = (L + C - S).H.solve(R.as_source())`.
- Verdict: passes — `.H` is on the OperatorSum, propagates to leaves.

**Test 4: α modes.**
- Spec form: `A_prompt φ = -α T φ`.
- Code form: `GeneralizedEigenproblem(A=L+C-S-F, B=T, sign=-1).solve()`.
- Verdict: passes.

**Test 5: the SN sweep.**
- Spec form: `ψ = L⁻¹ q` (or `(L+C)⁻¹ q` depending on the
  decision at §5.1).
- Code form: `psi = L.solve(q)` (or `(L + C).solve(q)`).
- Verdict: passes — but only after Step 2 unifies the closure.

**Test 6: the response-weighted reaction rate.**
- Spec form: `r = <Σ_a, ψ>` (a functional applied to ψ).
- Code form: `r = ReactionRateFunctional(sigma_a, region="fuel") @ psi`.
- Verdict: passes — functional `@` field is the `@` dunder on field.

All six tests pass under the four-operator vision. The user's
elegance standard is reachable in 5 migration steps.

---

## §7. What this memo deliberately does NOT recommend

- **No new top-level operator types beyond L, C, S, F, T**. The
  grand report adds B and several composed types (`A_loss`,
  `A_prompt`, `K`, `G_kinetic`) but these are SUMS / PRODUCTS /
  INVERSES of the four — not new primitives.
- **No `SNCellOperator`, `SNSweepOperator`, etc. as public types**.
  They survive only as private representation / strategy classes
  inside L.
- **No premature MoC / PN unification of L, S, F.** The grand
  report §1.5 (cross-solver scope) says L, S, F have the SAME names
  but DIFFERENT realisations per method. Keep them per-method until
  ≥2 working methods exist (per user's "unify after two instances"
  policy from `feedback_unify_after_two_instances.md`). MoC's L is
  not the same object as SN's L; they share an interface, not an
  implementation.
- **No removal of the sweep / matvec free functions until Step 5**.
  Backward-compatibility wrappers stay through the migration to
  protect the 11 regression snapshots.

---

## Pointers

- **Spec**: `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/plans/neutron_transport_grand_report_v3.md`
  - §1 (lines 25-110): operator vocabulary, sign conventions
  - §5.7 (lines 636-672): operator hierarchy
  - §6.3 (lines 714-737): operator dunders
  - §9 (lines 1170-1259): SN-specific structure; harmonic factoring
    of S
  - §16A (lines 2699-3327): boundary realisers
  - §20 (lines 4054-4113): fission representations
  - §21 (lines 4117-4269): k-eigenvalue
  - §22 (lines 4273-4480): α-eigenvalue, Laplace resolvent
  - §24 (lines 4527-4581): solver hierarchy
  - §29 (lines 5481-5618): end-to-end API sketches per method
  - §30 (lines 5622-5649): the one-expression target
  - §31 (lines 5653-5717): 20 final design rules
  - §32.6 (lines 5946-5982): LinearOperator base contract
  - §43 (lines 7375-7396): closing implementation principle

- **Explorer audit (superseded by this memo where they disagree)**:
  `/Users/rodrigo/git/nuclear/ORPHEUS/.claude/agent-memory/explorer/issue_196_sn_operator_architecture_audit.md`

- **Current code**:
  - `orpheus/numerics/operator.py:116-200` — `LinearOperator`
    Protocol and `LinearOperatorMixin` (foundation already in place)
  - `orpheus/sn/operator.py:1115-1472` — `SNStreamingOperator`
    (the L to be promoted; currently has split apply/solve closures)
  - `orpheus/sn/solver.py:1067-1289` — `_solve_fixed_source_si` and
    `_solve_fixed_source_krylov` (the targets of §3.2-3.3)
  - `orpheus/sn/boundary_realizer.py` — `BoundaryRealizer`s (the §16A.4
    contract already partially implemented)
  - `orpheus/numerics/iteration.py` — `SourceIteration` primitive
    home (per §24)
  - `orpheus/numerics/eigenvalue.py` — `PowerIteration` etc. home
    (per §21, §29.2)

- **Linked memories**:
  - `[[issue-168-phase-f-closeout]]` — the twin-bug that triggered
    the audit
  - `[[c188-curvilinear-realizer-unify]]` — the BC-realiser
    unification across slab + curvilinear
  - `[[wave_5_boundary_realizer]]` — the §16A pattern landing

## Headline summary

The user's four-operator target (L, C, S, F) matches the grand
report's §1 operator convention and §30 / §43 elegance principle
exactly. The explorer audit's seven proposed types collapse onto
the four: `SNSweepOperator` IS `L.solve`; `SNCellOperator` is a
sweep representation strategy of L; `AngularRedistribution` is a
curvilinear streaming detail of L; `PoleFaceAnchor` is a
SymmetryBoundary realisation composed into L; `AngularIntegration`
is the ℓ=0 component of S's moment projector; `PackedStructuredAdapter`
is a representation choice of L; `SNBoundaryFaceTraceOperator` IS
B owned by L. Migration is 5 atomic, independently-committable
steps; Step 2 closes Issue #196 manifestation #7 by construction
(unified WDD closure across L.apply and L.solve); Step 5 produces
`solve_sn_adjoint` in five lines as `A_loss.H.solve(response)`,
which is the elegance acceptance criterion.
