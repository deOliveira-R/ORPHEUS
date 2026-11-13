# The ontology chart: ORPHEUS as it is, as it should be, and the route between (living plan, opened 2026-11-13)

## Goal

Every concept ORPHEUS rests on has exactly one object, in the layer where its vocabulary suffices to define it, and every object carries exactly one concept. A solver then reads as the mathematics, the derived concepts (transpose, adjoint, inverse, pairing, condensation weight) fall out of the structure of the spaces, and a mistake of kind (a flux handed where a source is meant, a density summed as an amount) cannot be spelled.

This plan holds three things and keeps them apart:

1. **The as-is map**: the concept-to-class map of the production code today, measured.
2. **The to-be map**: the target ontology, as a typed concept graph. Every node cites the ruling that fixed it or is marked *proposed*. The proposed nodes are the discussion agenda.
3. **The chart**: the difference between the two, cut into units with a dependency order, each unit tagged by the workflow that builds it.

It is a living plan in the sense of `plan-authoring` §0: the ontology is being searched, so the plan is refined in place as the discussion runs, and no production code is written from it until the user rules that it is polished.

## Rulings (the user, 2026-11-13)

- **O1, the mode of work changes.** "Instead of making small changes like we have done so far, we're going to deliberately depart from this modus operandi." The route is planned whole, against an idealised target, before it is built.
- **O2, no hand-holding.** The orchestrator plans and executes with its own judgement; the user and the orchestrator discuss the chart until it is solid. The surgical-carve exception of the `workflows` rule (the main agent writes operator-algebra code turn by turn under the user's steering) is lifted for this campaign. The rule's source page changes only through the harness process, on the evidence this campaign produces.
- **O3, dynamic workflows build it.** Work is distributed to Opus and Sonnet agents through workflows. The division of labour is the one already measured in `workflow_orchestration.md`: Sonnet executes, Opus reviews, the orchestrator integrates and owns every ruling.
- **O4, the communication medium.** Knowledge graphs, ontology and taxonomy: the maps above are the shared objects the discussion points at.

Standing rulings this plan inherits and does not reopen: R1 to R17 of `true_ends.md` (#707); the sharpening order (discrete ordinates, diffusion and the homogeneous baseline first; collision probability, characteristics and Monte Carlo recycle the machinery later); references closed and insulated from production; effort is sizing, never a criterion. **R18 (the home of variance, `flux.dual(on="cells")`) is still NOT ruled**; it is an open node of the to-be map.

## Method of the as-is measurement

The Nexus graph is a code graph (classes, calls, imports); it ranks hubs, not concepts. The ontology is read from the **concept-to-class map in both directions**:

| reading | meaning |
|---|---|
| a concept with no class | missing machinery |
| a concept with more than one class | a twin |
| a class with more than one concept | a weld |
| a class with no concept | procedural residue |
| a class whose concept belongs to another layer | misplaced |

Instruments:

- `scratch/ontology/inventory.py`: an `ast` pass over `orpheus/**/*.py`. `[M]` 2026-11-13 at `c04f4ba8`: 789 classes in 241 files (positive control: it finds `FunctionSpace`, `LinearOperator`, `SNProblem`); 678 outside `derivations/`.
- `scratch/ontology/units.py`: assigns every class to exactly one of 15 census units (0 unassigned, asserted).
- Workflow `ontology-asis-census` (run `wf_04ea90f3-46f`): per unit, a Sonnet explorer writes `scratch/ontology/asis/<unit>.json` (every class with its concept, category and verdict; hand-rolled operations; missing machinery; guards), validated by `scratch/ontology/check_unit.py` for coverage and vocabulary; an Opus reviewer then tries to refute every verdict other than ONE and every missing-machinery claim, and reads at least 12 ONE verdicts for missed welds. Two further units: U16 extracts the concept inventory the theory corpus defines and checks each against the code; U17 groups the open issues and the plans by the machinery they say is missing.
- `scratch/ontology/structure.py` on the Nexus database rebuilt at `c04f4ba8`: the package dependency matrix, the inheritance forest, the hubs, the largest modules. Output: `scratch/ontology/structure.txt`.

## The as-is map

### Structure `[M]` 2026-11-13, `scratch/ontology/structure.txt`

- **Size by package** (classes, `ast`): numerics 267, transport 126, derivations 111, sn 98, geometry 57, reference 29, mesh 18, diffusion 16, data 15, mc 15, cp 9, moc 9, values 6, thermal hydraulics 4, specification 3, fuel 2, homogeneous 2, kinetics 2.
- **The inheritance forest** has 40 roots among production classes. The five largest: `LinearOperator` (80 descendants), `ContentIdentity` (66), `Member` (32), `transport.fields._bases.RolePair` (30), `ResponseRealizer` (19). `FunctionSpace` has 10 subclasses.
- **Production in-degree** (calls, type uses and inheritance from production code only, `derivations/` excluded): `FunctionSpace` 270, `LinearOperator` 188, `FullField` 174, `CoupledField` 100, `SNProblem` 76. Counting tests as sources, `SNProblem` is the most connected class in the tree (828 code edges).
- **Who leans on whom** (cross-package calls): numerics is called 367 times and calls out 86; transport 225 and 127; sn is called 12 times and calls out 315. Diffusion, homogeneous, cp, moc and mc are called by no production package (0), as the layering requires.
- **Import edges that point up the declared order: 0** among the core packages; `cp` imports `derivations` (3 edges). The call edges the graph shows pointing up (for example numerics to transport, 4) have no import edge beside them, so they are `[R]` mis-resolved method names in the graph, not layering violations; the layer gate, not this count, is the instrument.
- **Largest modules** (lines): `sn/loss_representation/__init__.py` 5482, `sn/solver.py` 3944, `numerics/operator.py` 3281 (34 classes), `geometry/manifold.py` 2355, `sn/angular/closure.py` 2351, `numerics/traced_memo.py` 2135 (55 classes), `sn/problem.py` 2042.
- **Tagged guards** under `orpheus/`: `ELEGANCE-DEBT[guard]` 11 (5 in derivations), `SCOPE-BOUNDARY[guard]` 28. `raise NotImplementedError`: 119 (49 in derivations; sn 21, transport 20, numerics 15).
- **One tag branched on in many places** (Nexus `discriminations`, production sites only): the coordinate system (`CoordSystem` / `AxisCoord`) at 12 sites in geometry, mesh, sn, cp and mc; `deck` at 9 sites (the admit and realize verbs of four per-method classes); the axis `origin` at 8; `inner_solver` (a string) at 3.

### The concept-to-class map `[M]` 2026-11-13, run `wf_04ea90f3-46f`

33 agents, 0 errors, 7.06 M sub-agent tokens, 40 min. The 15 class units cover 789 of 789 classes (the checker asserts coverage and vocabulary; 15 of 15 units exit 0). The Opus reviewers examined 1221 claims across the 15 class units and the theory unit (a claim is one class entry or one hand-rolled, missing or guard entry): 644 confirmed, 188 refuted, 389 amended, and 89 added. So an unreviewed census verdict is weaker evidence: 78 of 789 class entries carry no review mark, all 78 of them census-ONE. The maps are `scratch/ontology/asis/U*.json`; the aggregate is `scratch/ontology/aggregate.py`; every agent's headline findings are `scratch/ontology/headlines.md`.

**Verdicts after review, production only (678 classes, `derivations/` excluded):**

| verdict | classes | share |
|---|---|---|
| ONE (one concept, right home) | 324 | 48 % |
| RESIDUE (no concept of its own) | 95 | 14 % |
| GUARD (exists to refuse) | 78 | 12 % |
| TWIN (concept realised elsewhere too) | 70 | 10 % |
| WELD (two or more concepts fused) | 67 | 10 % |
| MISPLACED (wrong layer or package) | 44 | 6 % |

Beside the classes: 250 hand-rolled operations and 143 missing-machinery claims survived review, and 174 guards were listed (guards were not reviewed in 12 of 15 units).

**By subsystem** (ONE of N): numerics spaces 36 of 52; numerics operators 48 of 74; numerics solution side 36 of 47; generators 8 of 25; numerics infrastructure 0 of 69; geometry 39 of 57; mesh, specification, values, data, reference 46 of 71; transport fields 28 of 53; transport operators and schemes 6 of 33; transport boundary 19 of 40; sn problem side 27 of 49; sn solution side 17 of 49; diffusion and homogeneous 8 of 18; cp, moc, mc, kinetics, fuel, thermal hydraulics 6 of 41; derivations 33 of 111.

**The theory corpus against the code** (U16): 277 concepts the theory defines; 213 realised, 33 partial, 25 unrealised (16 of those "by ruling"), 6 realised twice. The corpus is ahead of the code on spaces (5 unrealised, among them the four transport spaces) and operators (9 unrealised).

**The backlog** (U17): 484 open issues, each placed once in one of 33 machinery groups; 64 plan files, 8 active and 12 paused. The largest core groups: discrete measure and quadrature identity (25 issues), boundary laws and realizers (25), SN spatial schemes (23), operator ends in the four spaces (21), space identity and membership (18), symmetry and quotients (17).

### What the map says, by theme

Each line is a reviewed finding with the unit that holds its evidence. A finding marked `[M]` was probed by the reviewer; the others are by reading at the cited lines.

**A. The space algebra is incomplete, and most other debt is downstream of it.**
- No factor carries variance. `Axis` has no variance field, `DualSpace` flips the whole space, and all 6 trace and ray factor kinds are axes-less named leaves, so the four spaces of R14 have nothing to attach to (U01).
- `FunctionSpace` is two realisations in one class: axis-built (identity = its axes) and named leaf (identity = a digest of a name). The boundary half-trace spaces are leaves, so `[M]` `AxisRetractionOperator.of(trace.outflow_space('xmin'), 'angular')` raises "not axis-built", and `PartialCurrentOperator` and `IsotropicEmissionOperator` are second spellings of the axis retraction and section (U01, U11).
- The tensor product has two realisations: `[M]` `of_axes(a, b) == of_axes(a) * of_axes(b)` is True, yet only the second is a `TensorProductSpace` with `.factors`, and one consumer asserts that class (U01).
- The direct sum has two: `FullFieldSpace` is the N = 2 case of `CoupledSpace` and repeats its verbs method for method; the offset table is computed by four separate loops (U01, U02, U16).
- `TensorProductOperator` is not a tensor product. `[M]` `(A & B).apply` is bit-identical to `(B @ A).apply`, and 2 of 2 production mints are `X & IdentityOperator()`. An operator tensor product on product spaces does not exist (U02).
- A factor's metric has two homes (`Axis.weights`; a dense Gram positioned on the space). `[M]` `Axis.weights` is stored independently of the origin that generated it. `[M]` `DenseMetric` admits an indefinite form. `[M]` the energy factor exists as two unequal axes (U01).
- The pairing object exists and is bypassed: `posing._pair` and `KEigenvalue.compute_production_rate` reduce with `np.sum` and a default weight of 1 (U03). Membership is spelled four ways (U02) and checked three ways on the composite (U08).
- No operator restricts an axis to a subset: the inflow and outflow halves of a trace are gathered by index at three sites, and eight angular reductions on octant slices re-spell the retraction (U12).
- Measures: `ReferenceMeasure` carries no mass, so 4π is typed at three sites and 2π/n at a fourth (U04); `[M]` Lebesgue measure on [−1, 1] has two unequal identities (U03); `[M]` a `GeneratingMeasure`'s identity is its name string (U01).

**B. Operators do not declare what they are.**
- Every transport operator declares domain = codomain and changes role inside its body (U09; `operator_map.md`: 22 of 64 fit the four spaces).
- Stored ends are a numerics concept homed in transport: `BoundOperator`, used by 15 of 80 operator classes (U09).
- `block_role` and `system_role` are written and propagated and never read: `[M]` 0 decision reads in `orpheus/` (U02).
- There is no method-agnostic splitting: `Splitting` is SN-only, while `GreenOperator`, `NearKernel.of_splitting` and the adjoint entries pass M and N loose. `Splitting` has no `.H`, so the adjoint Green is built by hand at two sites (U02, U12).
- Schemes and closures are not operators. `DiscretizationSchemeBase` fuses a closure with per-schedule kernels and hand-written transposes; the angular closure family is hooks with hand-coded transposes and also stores the spatial pairing (U09, U11).
- The private sweep body is three inverses selected by arguments: plain, windowed and scheduled (U11).
- One lift, two policies (`MultiplicationOperator` re-does `BulkLift`); one kernel split, two spellings (diagonal and off-diagonal of the P0 transfer); role subclasses are not tied to their datum, so a scattering operator over (n,2n) data constructs (U09).
- Diffusion has no gradient, divergence or scheme object, and `DiffusionBoundaryOperator` duplicates the SN one and has already diverged from it (U13).

**C. The Problem and Solution shape exists once, in the wrong layer.**
- `SolutionBase` and `Splitting` live in `sn/` (L3). No L3 package may import another, so diffusion, homogeneous, CP, MoC and MC each carry a forced twin result record. The concept is misplaced, not missing (U13, U14).
- One eigenproblem, two drivers: SN forward runs on `SNSolver`, SN adjoint on `KEigenvalue`. The k estimator has five spellings; the scale gauge is applied twice; two eigen-outcome types coexist (U03, U12, U13).
- `KEigenvalue` fuses the question with the outer strategy; no method-agnostic Strategy value exists (U03).
- The diffusion pencil and posings are minted on the solver, not on the hub (U13, U16).

**D. Geometry stops at rotations, and faces are strings.**
- No group type contains translations: `[M]` `RealSpace(3).quotient(O3 | Dinfh | O2_x)` raises. So the spatial orbit space lives outside the orbit-space machinery as `Chart`, a second orbit-space type (U06).
- `CoordSystem` is per-tag tables; "is it curvilinear?" has a derived answer and is still spelled as a tag test at 7 sites; `AxisCoord` is a second enum of the same kind (U06, U07).
- A face has three identities, and geometry spells faces as strings; `FaceSite` is geometry living in transport (U06, U07, U10).
- The structured mesh has three realisations and the region-to-material map three storages (U07).
- Point location has four spellings with four boundary conventions; the chord kernel has no production consumer, and production CP imports reference internals, including a private function (U06, U14, U15).

**E. The boundary declaration is typed; its algebra is not.**
- The wall-law vector space and its cone have no type; `[M]` one point of the cone can carry two unequal values. The sub-Markov bound is checked twice (U10).
- Realization is per chain, never per factor; the factor classes have no production reader (U10).
- The boundary source is a bare array, and no adjoint-side boundary term exists (U10).

**F. Generators and data.**
- The `Basis` base class makes every basis carry trial and test contractions, six bases retype the same contractions, and the Legendre family is realised twice (U04).
- `Quadrature` fuses a measure, a sweep side-channel, a fold tag and a factory catalogue; its selector layer has no production consumer (U04).
- The rule "which channel is weighted how" is written once for energy and once for space; the group structure is bypassed by a raw optional array on `Mixture` (U07).

**G. Not mathematics, and not alive.**
- The traced memo (69 classes) imports only the standard library, numpy and `orpheus.values`: by the lowest-layer criterion it belongs below `numerics/` (U05).
- No production caller: the quadrature selector layer; the dense eigen engines; four verbs of `HarmonicMomentFlux`; `TimedFullField`'s history half; `IntegralKernelOperator`; `TransferOperator`'s diagonal split; `thermal_hydraulics/solver_dae.py` (a line-for-line duplicate); the Sood result cache (U03, U04, U08, U09, U14, U15).

### Defects the census surfaced (to verify, then fix or file; none is fixed yet)

| defect | where | status |
|---|---|---|
| `GreenOperator.dual()` rebuilds the Green without its realisation, so the transpose of a Krylov-realised Green iterates by Richardson | `numerics/green_operator.py:598` | `[M]` read by the orchestrator: the mint passes `max_iter` and `tol` only. The fixed point is unchanged; the trajectory and the convergence are not |
| `AngularSourceSink.__add__` adds scalar-rate units into angular-rate units with no division by the angular mass | `transport/fields/angular_source_sink.py:146` | by reading (U08) |
| `MomentField` never checks its `L` against its space's angular head | `transport/fields/_bases.py:846` | by reading (U08) |
| `DenseMetric` admits an indefinite symmetric form | `numerics/metric.py` | `[M]` probe (U01) |
| `Axis` accepts weights that disagree with its origin's measure | `numerics/axis.py:436` | `[M]` probe (U01) |
| `GeneratingMeasure` identity is its name alone | `numerics/generating_measure.py` | `[M]` probe (U01) |
| CP averages every material id other than 1 and 2 into the coolant spectrum; `MoCResult.flux_fuel` returns zeros when material 2 is absent | `cp/solver.py:1029` | by reading (U14); later methods, filed with the measurement |
| MC's tally docstring names the wrong estimator | `mc/solver.py:578` | by reading (U14) |
| Stale docstrings: `operator.py:170` (role dispatch that does not exist), `axis.py:648`, `sn.py:2593` (a lambdify that does not happen) | | by reading |

## The to-be map

### Candidate organising laws (all *proposed* 2026-11-13 by the orchestrator, `[R]`, written before the census returned; each is confirmed, amended or refuted against it)

- **T1, five kinds and nothing else.** Every production object is a set (manifold, index set, partition), a measure on a set, a space (a factor or an expression in factors), an element of a space (a member or a co-vector), or a morphism between spaces (an operator). The Problem and Solution layer adds no sixth kind: a pencil is a pair of morphisms on the same ends, a posing is a pencil with its data, a splitting is a decomposition of a morphism, a solution is an element with its provenance.
- **T2, a space is an expression.** Spaces are built by four constructors only: a factor minted by a generator (a measure, a basis, a frame, a scheme), the tensor product, the direct sum, and the dual of a factor; plus a subspace or quotient by a declared projector. A bespoke space class is a weld of an expression with its consumers. Every factor, including the trace and boundary leaves, is the same kind of object and carries its variance (R14, R18 open).
- **T3, a member follows its space.** Flux, source, importance and response are the four variance patterns of R14, and a layout is a factor list. A member class per (role, layout) pair is the space's expression written a second time.
- **T4, an operator is its true ends plus its action per factor** (R2). The structural maps are the natural maps of the constructors, a closed family: per-factor action, contraction with a co-vector and insertion of a vector on a tensor factor (retraction, section, pullback); injection and projection of a direct sum (system restriction, extension by zero, lift); inclusion and restriction of a subspace (trace gather and scatter); the Riesz legs per factor; analysis and reconstruction of a frame.
- **T5, derived, never implemented.** Transpose, adjoint, inverse of a wrapper, pairing, norm and the weight of an average (R16) are read off the spaces.
- **T6, physics enters in two places.** As data, which are members of named spaces; and as a short list of atoms (streaming, multiplication by a field, the group-transfer kernel, the frame faces, the wall response). Every other transport operator is a composition of those.
- **T7, a method is a choice of generators plus a realisation of the implicit inverse.** One Problem type family, one pencil, one posing, one Strategy, one outcome and one Solution shape serve discrete ordinates, diffusion, the homogeneous baseline, and later collision probability (the inverse of streaming plus collision projected on flat sources), characteristics (the same inverse on a track space) and Monte Carlo (a stochastic realisation with estimators as functionals).
- **T8, coarsening is one construction**: the posing filtration walked backward. Homogenisation and condensation are its merges, with the kind of each slot typed (R16).
- **T9, a package holds one vocabulary.** Infrastructure (the traced memo, the import-time tracer) is not mathematics and leaves `numerics/`.

### The typed concept graph (79 concepts; source `scratch/ontology/tobe.py`, drawn by the orchestrator from the reviewed census, `[R]` as a design)

Status today uses the map's vocabulary: **one** (one class, right home), **twin**, **weld**, **misplaced**, **residue**, **partial** (the object exists and is bypassed, mis-declared or incomplete), **missing**. Counts: one 11, partial 24, missing 16, twin 14, weld 9, misplaced 3, residue 2. The basis column cites the ruling or issue that fixes the target, or says *proposed*: 28 of the 79 targets are proposed and are the discussion agenda.

| concept | today | target | basis | wave |
|---|---|---|---|---|
| **Values** | | | | |
| Admitted value, content identity | **one**: ContentIdentity (66 descendants) | unchanged | landed | - |
| Content digest | **twin**: one blake2b convention spelled 5 times (values, memo boot) | one digest object | proposed | W0 |
| **Sets and geometry** | | | | |
| Manifold, G-space | **one**: Manifold, GSpace | unchanged | landed | - |
| Symmetry group, translations included | **partial**: SubgroupOfO3 only; RealSpace(3).quotient(O3|Dinfh|O2_x) raises | a subgroup of E(d); the spatial symmetry group is a value | proposed | W6 |
| Orbit space | **twin**: Quotient; Chart (a second orbit-space type); DirectionDomain re-spells two catalogued quotients | one Quotient machinery for directions, positions and lines | proposed | W6 |
| Coordinate system | **weld**: CoordSystem per-tag tables; AxisCoord twin; 'curvilinear?' tag-tested at 7 sites; branched on at 12 | derived from the group: measure, density, direction chart, kept columns | boundary-law ontology (architecture E) | W6 |
| Face of a region | **twin**: 3 identities (FaceLabel, face-name bijection, endpoint strings); strings in geometry; FaceSite in transport | one typed Face in geometry (axis, sign, normal, partner) | proposed | W6 |
| Face role: deck or tag | **one**: MirrorDeck, TranslationDeck, int tag | unchanged | landed (wave 2) | - |
| Region map with faces | **partial**: 1-D only; no 2-D region map; Mesh2D has no region labels | a d-dimensional region map a mesh refines | #569 | W6 |
| Mesh: a partition laid over geometry | **twin**: Mesh1D, Mesh2D, AxisMesh tuples; legacy_adapter between them | one mesh; region-to-material assignment stored once (3 storages today) | proposed | W6 |
| Point location, chord | **twin**: 4 locators with 4 boundary conventions; chord kernel has 0 production consumers | ConcentricPartition is the one locator and chord | #596 | W6 |
| **Measures and generators** | | | | |
| Continuous reference measure | **partial**: ReferenceMeasure has no mass (4 pi typed at 3 sites); Lebesgue on [-1,1] has two identities | a measure that knows its mass and moments; one identity | proposed | W8a |
| Discrete measure | **weld**: DiscreteMeasure: point measure + exactness + phase tag + invariance + partition (1440 lines) | typed by its support (#624); the phase taxonomy question is the user's | ruling 2026-09-01 (to revisit?) | W8a |
| Quadrature | **weld**: measure + sweep side-channel + fold tag + factory catalogue; selector layer has 0 production consumers | a measure with an exactness claim; rule families as values | proposed | W8a |
| Basis | **weld**: ABC forces trial and test contractions; 6 bases retype one contraction; Legendre family twice | trial and test families apart; one table contraction | proposed | W8a |
| Frame (basis, measure) | **one**: FrameBase, GalerkinFrame, PetrovGalerkinFrame | unchanged | landed | - |
| Energy group structure | **partial**: EnergyGrid exists; Mixture.eg is a raw optional array read by CP, MC, MoC | every consumer reads the energy axis | proposed | W8a |
| Spatial scheme | **weld**: closure fused with per-schedule kernels and hand-written transposes | a generator minting an operator whose transpose is derived | proposed | W4 |
| Angular closure | **weld**: hooks, not an operator; stores the spatial pairing; hand-coded transposes | an operator on the angular cell-edge space | proposed | W4 |
| **Spaces** | | | | |
| Factor of a space, with its variance | **missing**: Axis has no variance; 6 of 12 factor kinds are axes-less leaves | every factor is one kind of object carrying function-or-integrated | R14; home R18 OPEN | W1 |
| Metric of a factor | **twin**: two homes (Axis.weights; a dense Gram on the space); weights stored beside the origin; DenseMetric admits indefinite forms | one home, derived from the generator | proposed | W1 |
| Tensor product of spaces | **twin**: of_axes(a,b) and of_axes(a)*of_axes(b) are equal but only one is a TensorProductSpace | one constructor | R2 | W1 |
| Direct sum of spaces | **twin**: CoupledSpace (N-ary) and FullFieldSpace (its N=2 case); 4 offset loops | one constructor; the full field is bulk (+) trace | posing unit 4 (#527) | W1 |
| Subspace on an index subset | **missing**: inflow/outflow halves and octants gathered by index | restriction of a factor is a space and a split pair | proposed | W1 |
| Dual of a factor | **missing**: DualSpace flips the whole space; a composite's dual is not a DualSpace | per-factor dual; the whole-space dual is the composition | R8, R14; R18 OPEN | W1 |
| Pairing | **partial**: Covector, InnerProductFunctional exist; posing._pair and the k estimator use np.sum | DualPairing: the only spelling of a contraction | R11, R15 | W1 |
| Boundary trace space | **twin**: 3 unrelated hierarchies (angular, scalar, radial half-angle), all leaves | axis-built trace spaces shared by all methods; measure includes face area | R4, R9 | W1 |
| The four transport spaces | **missing**: one space for flux, source, importance and response | four variance patterns of one factor list | R14 | W1 |
| Quotient by a metric kernel | **one**: FunctionSpace.rank, the kernel condition | unchanged | ruling 80 | - |
| **Elements** | | | | |
| Element of a space | **twin**: 4 member contracts (Vector, _is_ravellable, SystemField, CompositeField); 3 membership checks | one member contract; a member follows its space | R5 | W3 |
| Flux and source members | **partial**: RolePair forest (30 classes) and Member forest (32): a class per role and layout | role = variance pattern; layout = factor list | R5, R14 | W3 |
| Importance and response members | **missing**: AdjointSolution.importance returns ScalarFlux | two member kinds on their own spaces | R17 | W3 |
| Co-vector | **one**: Covector, FunctionSpace.flat | unchanged (and used) | landed | - |
| Cross-section field | **partial**: SigS densified 3 times; chi and nu-Sigma_f gathered twice | a member of a named space, gathered once | #648 | W3 |
| Boundary source | **missing**: the inflow sampler returns a bare array; no adjoint-side boundary term | a typed member on the inflow trace, and its response partner | proposed | W7 |
| Positive cone | **one**: a predicate, by ruling | unchanged | landed | - |
| **Operators** | | | | |
| Operator with true ends | **partial**: 22 of 64 concrete operators fit the four spaces; stored ends (BoundOperator) live in transport, used by 15 of 80 | ends are stored in numerics and are the members' true spaces | R2 | W2 |
| Sum and product | **one**: OperatorSum, OperatorProduct | unchanged | landed | - |
| Tensor product of operators | **missing**: TensorProductOperator composes in sequence: (A & B) is bit-identical to (B @ A) | A (x) B on V (x) W with factor-indexed ends | R2 | W2 |
| Block operator on a direct sum | **weld**: CoupledOperator: block grid + triangular substitution + LU + assembly | the grid is the operator; the substitution is a strategy | posing unit 5 (#528) | W2 |
| Transpose | **one**: LinearOperator.dual | unchanged | landed | - |
| Adjoint | **partial**: AdjointOperator dresses with whole-space Riesz legs and swaps ends | owned by the pairing: Riesz legs on the angular factor only | R11 | W2 |
| Inverse and its shortcuts | **partial**: 27 of 41 inverse/adjoint/dual methods return a stored object unchecked | a shortcut is taken where a predicate proves it | #705 | W2 |
| Natural maps of a product | **partial**: retraction, section, pullback exist on axis-built spaces; PartialCurrent and IsotropicEmission re-spell them on traces | one family on every factor | proposed | W2 |
| Restriction and extension | **twin**: trace gather/scatter, system restriction, BulkLift, and a second lift in MultiplicationOperator; no shared base | the injections and projections of a direct sum or a subspace | proposed; revisits 'restriction siblings share no base' | W2 |
| Analysis and reconstruction | **one**: FrameAnalysis, FrameReconstruction | head declared as two spaces | R7 | W2 |
| Splitting A = M - N | **misplaced**: Splitting is SN-only; numerics passes M and N loose; no .H | a numerics value with its adjoint | proposed | W2 |
| Green of a split | **weld**: inverse wrapper + splitting + driver loop + residual promise; dual() drops the realisation | the inverse of a Splitting under a realisation | proposed | W2 |
| Coarsening (homogenise, condense) | **missing**: the channel-weighting rule is written once for energy and once for space | one morphism; each slot's kind typed, weights from the pairing partner | R16 | W8b |
| Block and system role tags | **residue**: written and propagated; 0 production reads | retired | measured | W0 |
| **Transport physics** | | | | |
| Streaming | **partial**: declares full_field -> same; the coefficient has 2 homes (a 3rd flagged, #240) | flux -> source; one coefficient object per chart | R2, R6 | W4 |
| Multiplication by a field | **partial**: declares one space; re-does the lift | flux -> source; composed with the one lift | R2 | W4 |
| Group transfer, S = R Lambda M | **partial**: Lambda declares the head only; acts on head x group x cells | true ends; role tied to its reaction channel | R1, R2 | W4 |
| Fission production | **partial**: declares one space | flux -> source, rank-one in energy | R2 | W4 |
| The implicit inverse (sweep) | **weld**: one private body is three inverses: plain, windowed, scheduled | three operators; source -> flux | proposed | W4 |
| Gradient and divergence | **missing**: LeakageOperator fuses Fick's map, the divergence and the P1 closure; stencil written twice | a typed pair under volume and area metrics | proposed | W4 |
| Wall law on its cone | **partial**: WallCondition is typed; the cone and its vector space are a free predicate and a mixin | a typed cone; one value per point | proposed | W7 |
| Realized boundary operator | **twin**: SN and diffusion boundary operators duplicated and diverged; the ray corner re-realises the specular law by hand | per-factor realization; one operator on the trace | proposed | W7 |
| Reaction rate | **partial**: typed for SN, diffusion, homogeneous; adjoint covector rebuilt by hand | the pairing of a response with a flux | R14 | W5 |
| **Problem** | | | | |
| Problem | **partial**: 3 hubs, no shared type (MethodMesh covers 2) | one Problem protocol; every hub mints its spaces and its pencil | posing sequence (#522) | W5 |
| Operator pencil | **partial**: typed LinearOperator[V, V]; diffusion mints it on the solver | (A, F): flux -> source, minted on the hub | R14 | W5 |
| Posing, and its adjoint | **partial**: the adjoint fixed-source posing is spelled by hand | the adjoint posing is the posing's dagger | #484 | W5 |
| Question bound to a system | **missing**: Eigen, FixedSource, Response exist; nothing binds them to operators | a map from a question to a pencil and a spectral map | #529 | W5 |
| **Solution** | | | | |
| Strategy | **missing**: KEigenvalue fuses question and strategy; inner algorithm is a string at 3 sites | a method-agnostic value: splitting, realisation, stopping rule | proposed | W5 |
| Eigen driver and estimator | **twin**: SN forward on SNSolver, adjoint on KEigenvalue; k estimator spelled 5 times; gauge applied twice | one driver; k is the posing's Rayleigh quotient | proposed | W5 |
| Iteration record | **one**: IterationRecord (verdict derived) | unchanged | landed (#340) | - |
| Outcome | **twin**: PowerIterationOutcome and EigenOutcome | one outcome, gauge recorded | proposed | W5 |
| Solution five-tuple | **misplaced**: SolutionBase lives in sn/ (L3): 5 methods carry forced twin records | in transport/ (L2), consumed by every method | proposed | W5 |
| Observable | **missing**: admission re-spelled per method | a co-vector on the flux | #405 P4 | W5 |
| **Methods** | | | | |
| Discrete ordinates | **partial**: on the algebra; 98 classes, 27 + 17 one-to-one of 98 | the reference realisation of the shape | sharpening order | W4-W5 |
| Diffusion | **partial**: operators on the algebra; no hub pencil, no scheme, no Solution shape | same shape as SN | sharpening order | W4-W5 |
| Homogeneous baseline | **partial**: pencil and outcome; a result record twin | same shape | sharpening order | W5 |
| Collision probability | **missing**: 0 operator classes; imports reference internals | the implicit inverse projected on flat sources, as an integral operator | later, recycling | W9 |
| Characteristics | **missing**: 0 operator classes; own quadrature carrier, own locator | the implicit inverse on a track space | later, recycling | W9 |
| Monte Carlo | **missing**: 0 operator classes | a stochastic realisation; estimators are functionals | later, recycling | W9 |
| **Outside the algebra** | | | | |
| Traced memo and import tracer | **misplaced**: 69 classes in numerics/; import only stdlib, numpy and values | a package of their own below numerics | proposed | W0 |
| Code with no production caller | **residue**: 8 groups (quadrature selector, dense eigen engines, IntegralKernelOperator, solver_dae.py, ...) | retired, or moved to the test side | measured | W0 |
| Reference insulation | **partial**: production CP imports 3 derivations symbols, one private; the gate forbids only the other direction | no production import of a reference | ruling 2026-10-06 | W9 |
| The concept map as an instrument | **missing**: conceptual_view.rst is hand-written; most 'implements' edges in the graph are inferred from names (the Nexus tool's own figure for this project: 12 999 of 13 084, not re-measured) | a tracked concept registry; a gate holds one class per concept; the page is generated | proposed | W0 |

## The chart

Eleven waves. A wave is a set of units with one outcome; inside a wave the units are cut when its specification is written. The order below is the dependency order, not a schedule: W1, W6 and W8a have no dependency on each other and run side by side after W0.

| wave | outcome | after | workflows | carries |
|---|---|---|---|---|
| **W0 Clear the ground** | Retire code with no caller; verify and fix the surfaced defects; move the traced memo below numerics; make the concept map a tracked, gated registry. | none | build-units | file-disjoint, mechanical |
| **W1 The space algebra** | Every factor one kind of object with variance and one metric home; one tensor product, one direct sum, subspaces; the pairing as the only contraction; axis-built trace spaces; the four spaces. | W0, W8a | design-panel, spec-check, build-units | absorbs #707 A2 and posing unit 4 |
| **W2 Operators on true ends** | Stored ends in numerics; an operator tensor product; natural maps on every factor; the adjoint from the pairing; a splitting value with its adjoint; the Green unwelded; predicted shortcuts. | W1 | spec-check, build-units | absorbs #707 A4-A5, #705, posing unit 5 |
| **W3 Members follow spaces** | One member contract; flux, source, importance and response as variance patterns; cross-section fields gathered once. | W1 | spec-check, build-units | absorbs #707 A3, #664, #633, #666 |
| **W4 Transport atoms** | Loss and gain re-bound flux to source; reaction role tied to its channel; schemes and closures as operators with derived transposes; the sweep's three inverses apart; gradient and divergence for diffusion. | W2, W3 | build-units (sequential slices on one branch) | absorbs #707 A6 and Steps B-D |
| **W5 One Problem, one Solution** | Solution and Strategy in transport; one eigen driver and one estimator; one outcome; every hub mints its pencil; the adjoint posing derived; questions bound to systems; observables as co-vectors. | W4 | design-panel, build-units | absorbs posing units 5-6, #484, #529 |
| **W6 Geometry** | Groups with translations; one orbit-space machinery; the coordinate system derived; a typed face; one mesh with region labels; one locator. | W0 | design-panel, build-units | absorbs boundary-law architecture E, #584, #569, #596 |
| **W7 Boundary algebra** | The wall-law cone typed; per-factor realization; one boundary operator for SN, diffusion and the ray corner; the adjoint side of a wall law. | W1, W3, W6 | spec-check, build-units |  |
| **W8a Measures and generators** | A reference measure with mass; discrete measure typed by support; quadrature and basis unwelded; group structure read from the axis. | W0 | design-panel, build-units | absorbs posing unit 2, #624 |
| **W8b Coarsening** | Homogenisation and condensation as one morphism with typed slots. | W3, W8a | spec-check, build-units | R16; posing unit 3 |
| **W9 Later methods on the machinery** | Collision probability, characteristics and Monte Carlo re-spelled on the shared Problem, operators and Solution; reference insulation restored. | W5, W6, W7 | one design-panel per method, build-units | sharpening order |

**Three lanes.** The numerics lane is W1, W2, W3, W4, W5: spaces, then operators and members, then the transport atoms, then the Problem and Solution shape. The geometry lane is W6 then W7. The generators lane is W8a then W8b. The lanes meet at W7 (the boundary needs trace spaces and typed faces) and at W9 (a later method needs all three).

**Why this order** `[R]`. Theme A of the map is upstream of themes B, C and E: an operator cannot declare flux to source until a factor carries variance; a boundary operator cannot be one object for discrete ordinates and diffusion until the trace is an axis-built space; the Solution shape cannot be shared until the members it holds are typed by their spaces. W0 goes first because it shrinks the surface every later wave touches and because a known defect is fixed before work is built on it (`process-discipline`).

**The campaigns in flight are placed, not superseded** (each row is a proposal for the user to confirm):

| campaign | state today | place on the chart |
|---|---|---|
| #707 true ends | Step A at A1; 118 gates red on `refactor/true-ends`; R18 open | its steps become the first slices of W1 (A2), W3 (A3), W2 (A4, A5) and W4 (A6, B to D). The red gates stay the gates. |
| #705 predicted shortcuts | `7615c845` parked on `fix/predicted-shortcuts`, waits on #707 | W2 |
| posing sequence #522, units 0 to 6 | ontology ruled; units scheduled, not opened | unit 2 in W8a, unit 3 in W8b, unit 4 in W1, unit 5 in W2 and W5, unit 6 in W5; units 0 and 1 in W0 |
| boundary-law ontology, architecture E | ruled; waits on #405 | W6 (the coordinate system as a group) and W7 |
| layer DAG, edge contracts | active, paused | every wave moves classes between layers, so the layer gate and the concept order are part of each wave's bar; edge contracts follow W5 |
| reference architecture #405 | P0 to P3 landed; P4 next | independent of the chart (references are closed); W9 restores the insulation |
| verification tree refinement | paused for #707 | resumes on a tree whose ontology has moved; its kill matrix is re-measured after W5 |

**Workflows.** The templates of `workflow_orchestration.md` are used as they stand: `design-panel` where a target is *proposed* and the ontology is still being searched; `spec-check` to turn ruled targets into one specification with its leaf units and dependency table; `build-units` for the gates-first build, review and amend in worktrees; `docs-pass` after each wave. A user checkpoint sits between workflows, never inside one. The census run is the first evidence for the scale: 33 agents, 0 errors, 7.06 M tokens, 40 min, and the Opus pass changed 577 of 1 221 claims (188 refuted, 389 amended), which is the case for keeping a reviewer on every executor stage.

**Sizing** `[R]`, as duration and never as a reason: each wave is one design or specification workflow, two to four build workflows and one documentation pass; W1, W4 and W5 are the largest. The full suite (about 75 min) runs once per integration.

## Open questions for the user

1. **R18, the home of variance** (`flux.dual(on="cells")` against the alternatives in `true_ends.md`). W1 cannot be specified without it.
2. **Fold #707 into the chart?** Proposal: yes, as the first slices of W1 to W4, on a fresh branch off `main` that takes A1 with it. The alternative is to finish #707 as planned and start the chart after it.
3. **The discrete measure's phase taxonomy** (angular, spatial, energy inside `numerics/`, ruled 2026-09-01). The census reads it as the reason `DiscreteMeasure` is a weld and a transport vocabulary sits in L1. Keep, or re-rule in W8a?
4. **Production CP imports reference internals** (`cp/solver.py:49-51`, one private symbol), and the layer gate forbids only the other direction. Debt to remove in W0 (a gate edge plus a local kernel), or accepted until W9?
5. **Where the traced memo lives.** It needs only `values`; a package of its own between `values` and `numerics` is the lowest-knowledge home.
6. **`restriction siblings share no base`** was ruled because no consumer treated a restriction generically. W2's natural-map family is that consumer. Reopen?
7. **The earlier ruling that kept `Mixture.condense` and `MaterialXSField.project_through` as separate verbs.** R16 (per-slot kind, weight from the pairing partner) appears to supersede it. Confirm.
8. **Scope of the references.** `derivations/` is 111 classes, 33 one-to-one. The chart leaves it closed, except for the insulation. Confirm.

## Refuted candidates

- **"Read the ontology off graph degree."** `[REFUTED 2026-11-13]` for the question "what are the concepts": degree ranks hubs. The top production class by code edges is `SNProblem` (828), and the top node of the whole graph is a test helper. The fact it establishes: where a change costs most.
- **"The call edges pointing up the layer order are layering violations."** `[REFUTED 2026-11-13]`: they have no import edge beside them (0 import edges up the order among the core packages), so they are mis-resolved method names in the graph.
- **Census claims of missing machinery that the reviewers refuted because the object exists**: the Riesz legs, the n-ary direct sum (`CoupledSpace`), point evaluation of a basis, a Gram object, a space norm, a factorisation, the angular section with its mass (`AxisSectionOperator`), the measure on oriented lines (`LineDomain`), `Splitting` and `PredictedKernel` on the SN side, a Pin protocol in the memo. The fact each establishes: the debt is consumption or placement, not a new object. This is why the chart's most common status is *partial* (24 of 79) and not *missing* (16).
- **Candidate law T1 as first written ("five kinds and nothing else").** `[REFUTED 2026-11-13]` as a description of the tree: 78 production classes are refusals and 95 are records with no concept of their own, most of them result and outcome records. Amended: five kinds of mathematical object, plus values that record a run (records, outcomes, certificates, refusals), which are a sixth kind with its own home.

## Implementation starts when

The user rules the chart polished. Until then this campaign writes plans, measurements and maps only.
