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
- **O2, no hand-holding.** "I had to hold hands a lot with previous models ... now we're going to cooperate without handholding." The orchestrator plans and executes with its own judgement; the user and the orchestrator discuss the chart until it is solid. Confirmed by O7.
- **O3, dynamic workflows build it.** Work is distributed to Opus and Sonnet agents through workflows. The division of labour is the one already measured in `workflow_orchestration.md`: Sonnet executes, Opus reviews, the orchestrator integrates and owns every ruling. `[M]` since 2026-11-16 the division is measured per stage: the table of stage, model, agent and basis is on the workflows page (`docs/development/workflows.md` on `main`, "The stages of a workflow template: model and agent"), and "How a wave runs" below cites it.
- **O4, the communication medium.** Knowledge graphs, ontology and taxonomy: the maps above are the shared objects the discussion points at.

**Rulings of the second exchange (the user, 2026-11-13):**

- **O5, the wide review was the user's plan for R18.** "My plan was this wide review to spaces in a wider context. I think this review reached the correct conclusion (or at least is very close to it)." So the home of variance is no longer a separate question: it is settled inside W1 (the space algebra), whose finding is that every factor carries its variance. The exact spelling is W1's specification.
- **O6, #707 is folded into the chart.** Its steps are the first slices of W1 to W4. The open issues are censused against the chart, and the insights they carry are used.
- **O7, the surgical-carve exception is lifted.** "As long as the plan is well specified in a tasks, Sonnets should be capable of implementing and Opus should be capable of reviewing it and asking for corrections. In the end, as orchestrator of the workflow, you are responsible for the results, so you do a final review that it reached your standards of acceptance." The orchestrator's final review of every integrated unit is therefore a required stage, not an option.
- **O8, `derivations/` stays closed** until production is in order.

**Rulings of the third exchange (the user, 2026-11-13):**

- **O9, the two-coordinate law, the cycle law and the inverse/realisation/solver split (T10 to T12) are accepted**, with phases as sub-packages inside each layer. "Yes. I agree with that proposal."
- **O10, the first workflow may be W0's specification.**
- **O11, a preparatory phase precedes every build workflow: the harness evaluation.** The rules, skills and agent definitions were written for older models. They are evaluated per model (Sonnet, Opus, Fable), and what no longer changes an outcome is trimmed or retired, so that a workflow agent starts lighter. It has its own plan, `.claude/plans/harness_eval.md`, and its own branch off `main`; its findings also set which model each workflow stage is given.
- **The plan's status.** The user: it "more or less matches my idealized vision", "we're generally in agreement and going in the same direction". That is agreement on the direction, not yet the ruling that the plan is polished; implementation waits for the harness evaluation and for that ruling. `[REMEDIED 2026-11-16]` the harness evaluation is closed (its ruling H22; record `harness_eval.md` on `main`), so only the ruling that the plan is polished remains. What the evaluation changes in this plan is in "How a wave runs", proposed 2026-11-17.

**Rulings of the fourth exchange (the user, 2026-11-20: "I agree with your 8 proposals"):**

- **O12, the rewrite has its own long-lived branch, `refactor/ontology`, cut from `main`** once every branch that holds a fragment of the starting point is merged there ("it will be a significant departure from the current shape of main"). Each wave's unit branches are cut from it and merge back into it; `main` stays as it is until the user decides the rewrite replaces it. This supersedes "W0, W6 and W8a are cut from `main` and merge as they finish" under "Branches".
- **O13, the revised "How a wave runs" is accepted**, with: two qa passes for every unit that mints a type or states a law (W1 to W8), one for W0's moves and retirements, revisited from W0's record; W0 runs on the measured assignment and is itself the record for the unmeasured unit kinds, read before W1's build is launched.
- **O14, W0 removes production CP's imports of reference internals (#641)**, and **the traced memo gets a package of its own between `values` and `numerics`**.
- **O15, the discrete measure's phase taxonomy (ruled 2026-09-01) is re-opened in W8a's design panel.**
- **O16, the placement of the campaigns in flight is confirmed as tabled.**
- **O17, the 26 targets still marked *proposed* are settled in their wave's design panel** and come to the user at that wave's specification checkpoint, not one by one now.
- **Starting point (option A, 2026-11-20).** `refactor/true-ends` merges to `main` with its 118 red Step A gates marked `xfail(strict=True)` (`c1072743`), after the full bar. `fix/predicted-shortcuts` and `park/u3b-sn-rcond` stay parked and are rebased onto the rewrite branch when their wave opens.

**Rulings of the fifth exchange (the user, 2026-11-20: "I agree with the A1-A3 proposals and have no objection. I think the plan is polished. Let's begin."):**

- **O18, A1, A2 and A3 of the adversarial round are accepted as proposed**: a role is a side of the trial and test pair and the variance pattern is derived (W1's panel takes it as its premise); `derivations/` is closed to redesign and new capability while each unit re-points the imports its change breaks, with no shim; T12 forbids a structure object that depends on a strategy value, and a kernel is declared by its generator.
- **O19, no objection to B1 to B13**; they stand as the orchestrator's decisions.
- **O20, the plan is polished. Implementation starts**, with W0, on `refactor/ontology` (cut from `main` at `eb45ee46`, which carries the starting point of option A).

**Rulings of the sixth exchange (the user, 2026-11-20, at W0's specification checkpoint: "Agreed" to items 1, 2, 4, 5 and 6; a question on item 3):**

- **O21, T11 is re-worded.** The concept graph is acyclic on inheritance, strict from W0; it is acyclic on references except inside one algebra closure (a set of classes with a hub that every member inherits), as a ratchet whose pending entries each name the wave that removes them.
- **O22, "one class per concept" means one head per concept.** A head is one definition site: a class, or a named closed sum whose arms are listed. A concept whose members obey different laws is split (the kernel pair and the gauge; the realizer and the realized operator).
- **O23, the vocabulary.** The six method concepts and the two work items (`insulation`, `registry`) leave the class registry; the five targets not yet built (`four`, `quotient`, `impresp`, `cellmoment`, `merge`) stay, marked missing with their wave; the 12 extension concepts are admitted after the eight two-sense concepts are split.
- **O24, the sixth kind stays "values that record a run".** Specification, citation and parameter records are data, not records.
- **O25, the builds of the specifications' parts A and B may start before the registry is settled.** The registry integrates last and reads final module paths.
- **O26, the kind of the traced memo** (the user, 2026-11-20: "I accept your recommendation on item 3", after asking what the memo does). Its manifest, pins and lookup verdicts are sixth-kind values: they record and certify a generation run. Its recorder, boot, import tracer and store take a seventh kind, `apparatus`: machinery that observes the program and not the problem, with no phase. T1 as amended therefore reads: five kinds of mathematical object, the sixth kind of values that record a run, and `apparatus`.

Standing rulings this plan inherits and does not reopen: R1 to R17 of `true_ends.md` (#707); the sharpening order (discrete ordinates, diffusion and the homogeneous baseline first; collision probability, characteristics and Monte Carlo recycle the machinery later); references closed and insulated from production; effort is sizing, never a criterion. R18 is resolved by O5 as a direction (variance per factor); its spelling is W1's.

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
| ONE (one concept, right home) | 324 | 48 %, an upper bound |
| RESIDUE (no concept of its own) | 95 | 14 % |
| GUARD (exists to refuse) | 78 | 12 % |
| TWIN (concept realised elsewhere too) | 70 | 10 % |
| WELD (two or more concepts fused) | 67 | 10 % |
| MISPLACED (wrong layer or package) | 44 | 6 % |

Beside the classes: 250 hand-rolled operations and 143 missing-machinery claims survived review, and 174 guards were listed, of which 56 were reviewed (in 7 of 15 units).

The ONE share is an upper bound and the rows below are not strictly comparable: ONE verdicts were sampled, not all read (78 unreviewed), and each reviewer applied the weld test in their own way (several state their rule: "verbs of one concept are not a weld").

**By subsystem** (ONE of N): numerics spaces 36 of 52; numerics operators 48 of 74; numerics solution side 36 of 47; generators 8 of 25; numerics infrastructure 0 of 69; geometry 39 of 57; mesh, specification, values, data, reference 46 of 71; transport fields 28 of 53; transport operators and schemes 6 of 33; transport boundary 19 of 40; sn problem side 27 of 49; sn solution side 17 of 49; diffusion and homogeneous 8 of 18; cp, moc, mc, kinetics, fuel, thermal hydraulics 6 of 41; derivations 33 of 111.

**The theory corpus against the code** (U16): 277 concepts the theory defines; 213 realised, 33 partial, 25 unrealised (15 of those "by ruling"), 6 realised twice. The corpus is ahead of the code on spaces (5 unrealised, among them the four transport spaces) and operators (9 unrealised).

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

### Defects the census surfaced (none is fixed yet; each is filed or already tracked: #708 to #713 opened 2026-11-13, the duplicate `thermal_hydraulics/solver_dae.py` as #713, production CP's private import as #641)

| defect | where | status |
|---|---|---|
| `GreenOperator.dual()` rebuilds the Green without its realisation, so the transpose of a Krylov-realised Green iterates by Richardson | `numerics/green_operator.py:598` | `[M]` read by the orchestrator: the mint passes `max_iter` and `tol` only. The fixed point is unchanged; the trajectory and the convergence are not. Already tracked: #629 |
| `AngularSourceSink.__add__` adds scalar-rate units into angular-rate units with no division by the angular mass | `transport/fields/angular_source_sink.py:146` | by reading (U08); filed #708 |
| `MomentField` never checks its `L` against its space's angular head | `transport/fields/_bases.py:846` | by reading (U08); filed #709 |
| `DenseMetric` admits an indefinite symmetric form | `numerics/metric.py` | `[M]` probe (U01); commented on #706 |
| `Axis` accepts weights that disagree with its origin's measure | `numerics/axis.py:436` | `[M]` probe (U01); filed #710 |
| `GeneratingMeasure` identity is its name alone | `numerics/generating_measure.py` | `[M]` probe (U01); filed #711 |
| CP averages every material id other than 1 and 2 into the coolant spectrum; `MoCResult.flux_fuel` returns zeros when material 2 is absent | `cp/solver.py:1029` | by reading (U14); filed #712 |
| MC's tally docstring names the wrong estimator | `mc/solver.py:578` | by reading (U14); the tally is tracked in #519 |
| Stale docstrings: `operator.py:170` (role dispatch that does not exist), `axis.py:648`, `sn.py:2593` (a lambdify that does not happen) | | by reading |

## The to-be map

### Candidate organising laws (all *proposed* 2026-11-13 by the orchestrator, `[R]`, written before the census returned; the verdict of the census on each follows the list)

- **T1, five kinds and nothing else.** `[REFUTED 2026-11-13]` as written; see the verdict below. Every production object is a set (manifold, index set, partition), a measure on a set, a space (a factor or an expression in factors), an element of a space (a member or a co-vector), or a morphism between spaces (an operator). The Problem and Solution layer adds no sixth kind: a pencil is a pair of morphisms on the same ends, a posing is a pencil with its data, a splitting is a decomposition of a morphism, a solution is an element with its provenance.
- **T2, a space is an expression.** Spaces are built by four constructors only: a factor minted by a generator (a measure, a basis, a frame, a scheme), the tensor product, the direct sum, and the dual of a factor; plus a subspace or quotient by a declared projector. A bespoke space class is a weld of an expression with its consumers. Every factor, including the trace and boundary leaves, is the same kind of object and carries its variance (R14; R18 resolved as a direction by O5).
- **T3, a member follows its space.** Flux, source, importance and response are the four variance patterns of R14, and a layout is a factor list. A member class per (role, layout) pair is the space's expression written a second time.
- **T4, an operator is its true ends plus its action per factor** (R2). The structural maps are the natural maps of the constructors, a closed family: per-factor action, contraction with a co-vector and insertion of a vector on a tensor factor (retraction, section, pullback); injection and projection of a direct sum (system restriction, extension by zero, lift); inclusion and restriction of a subspace (trace gather and scatter); the Riesz legs per factor; analysis and reconstruction of a frame.
- **T5, derived, never implemented.** Transpose, adjoint, inverse of a wrapper, pairing, norm and the weight of an average (R16) are read off the spaces.
- **T6, physics enters in two places.** As data, which are members of named spaces; and as a short list of atoms (streaming, multiplication by a field, the group-transfer kernel, the frame faces, the wall response). Every other transport operator is a composition of those.
- **T7, a method is a choice of generators plus a realisation of the implicit inverse.** One Problem type family, one pencil, one posing, one Strategy, one outcome and one Solution shape serve discrete ordinates, diffusion, the homogeneous baseline, and later collision probability (the inverse of streaming plus collision projected on flat sources), characteristics (the same inverse on a track space) and Monte Carlo (a stochastic realisation with estimators as functionals).
- **T8, coarsening is one construction**: the posing filtration walked backward. Homogenisation and condensation are its merges, with the kind of each slot typed (R16).
- **T9, a package holds one vocabulary.** Infrastructure (the traced memo, the import-time tracer) is not mathematics and leaves `numerics/`.


**The census against the laws** (2026-11-13; "supported" means the reviewed findings are instances of the law's failure, which is evidence that the law names real debt, not that the target is right):

| law | verdict | evidence |
|---|---|---|
| T1 | refuted as written, amended | 78 production classes are refusals and 95 are records. Amended: five kinds of mathematical object, plus values that record a run (records, outcomes, certificates, refusals) with a home of their own. |
| T2 | supported | 10 `FunctionSpace` subclasses; the direct sum and the tensor product each realised twice; 6 leaf factor kinds; the retraction refuses a half-trace. |
| T3 | open | The two forests exist (`Member` 32, `RolePair` 30). The field review (U08) refuted 6 welds there and found a missing type (flat against moment-carrying), so "fewer member classes" is not established. For a design panel. |
| T4 | supported | 42 of 64 operators mis-declare or do not fit; partial current and isotropic emission re-spell the retraction and section; two lifts. |
| T5 | supported | Hand-written transposes in the schemes, the closures and the 1-D walks; pairings by `np.sum`; the adjoint Green built by hand. |
| T6 | supported in part | The role subclasses are two constants each; the atom list itself (which operators are atoms) was not measured. |
| T7 | supported, amended | The shape exists once, in `sn/`: the debt is its layer, not its absence. |
| T8 | supported | The channel-weighting rule is written twice. An earlier ruling kept the two verbs apart (open question 7). |
| T9 | supported, amended | The traced memo belongs below `numerics/`, not beside it (it imports only `values`). |

**How complete the target graph is** `[M]` 2026-11-13 (`scratch/ontology/placement_build.py`, a Sonnet categoriser; the script asserts every finding is assigned once). The first draft had 79 concepts, drawn from the reviewers' headline findings. Of the 280 findings that survived review (67 welds and 70 twins in production, 143 missing-machinery claims), 214 were placed on one of the 79 and 66 were not. Four concepts were added for the four largest unplaced groups (29 findings): the sweep schedule and its walk (17), method parameter records (5), the kernel pair and gauge of an operator (4), the cell-moment factor and its member (3). **37 findings remain unplaced**, each a group of one to three: the traced memo's internals (3), the flat chart of a space, a reference verdict and order policy, the time history of a state, the d = 1 against d-generic linear-discontinuous kernels, a tag-keyed declaration, the DSA low-order system, the infinite-medium pose, the cylinder to equal-area-square map (2 each), and 18 singletons, listed in `scratch/ontology/placement.json`. Of the first 79 concepts, 12 had no finding placed on them (`compose`, `cone`, `dual`, `faces`, `fission`, `homog`, `insulation`, `inverse`, `quotient`, `registry`, `role`, `sn`): 5 are already in place, `registry` is new, and the other 6 rest on #707, #705, the operator map or the census's hand-rolled entries, none of which the placement pass read. The placement was not reviewed by a second agent.

### The typed concept graph (83 concepts; source `scratch/ontology/tobe.py`, drawn by the orchestrator from the reviewed census, `[R]` as a design)

Status today uses the map's vocabulary: **one** (one class, right home), **twin**, **weld**, **misplaced**, **residue**, **partial** (the object exists and is bypassed, mis-declared or incomplete), **missing**. Counts: one 10, partial 25, missing 17, twin 16, weld 10, misplaced 3, residue 2. The basis column cites the ruling or issue that fixes the target, or says *proposed*: 26 of the 83 targets are proposed (34 before the reconciliation with the posing sequence, which had already ruled 8 of them) and are the discussion agenda.

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
| Discrete measure | **weld**: DiscreteMeasure: point measure + exactness + phase tag + invariance + partition (1440 lines) | typed by its support (#624); + widened to differing supports (the disjoint union) | posing sequence, ruled 2026-09-28; taxonomy ruling 2026-09-01 open | W8a |
| Quadrature | **weld**: measure + sweep side-channel + fold tag + factory catalogue; selector layer has 0 production consumers | a measure with an exactness claim; rule families as values | proposed | W8a |
| Basis | **weld**: ABC forces trial and test contractions; 6 bases retype one contraction; Legendre family twice | partition-of-unity basis as the parent of IndicatorBasis; trial and test families apart; one table contraction | posing sequence, ruled 2026-09-28 (hierarchy); the rest proposed | W8a |
| Frame (basis, measure) | **one**: FrameBase, GalerkinFrame, PetrovGalerkinFrame | unchanged | landed | - |
| Energy group structure | **partial**: EnergyGrid exists; Mixture.eg is a raw optional array read by CP, MC, MoC | every consumer reads the energy axis | proposed | W8a |
| Spatial scheme | **weld**: closure fused with per-schedule kernels and hand-written transposes | SpatialScheme: one Petrov-Galerkin family; an operator whose transpose is derived | posing sequence (name ruled 2026-09-28); transposes proposed | W4 |
| Angular closure | **weld**: hooks, not an operator; stores the spatial pairing; hand-coded transposes | an operator on the angular cell-edge space | proposed | W4 |
| **Spaces** | | | | |
| Factor of a space, with its variance | **missing**: Axis has no variance; 6 of 12 factor kinds are axes-less leaves | every factor is one kind of object carrying function-or-integrated | R14; variance per factor (O5), spelling W1's; see A1 of the adversarial round | W1 |
| Metric of a factor | **twin**: two homes (Axis.weights; a dense Gram on the space); weights stored beside the origin; DenseMetric admits indefinite forms | one home, derived from the generator | proposed | W1 |
| Tensor product of spaces | **twin**: of_axes(a,b) and of_axes(a)*of_axes(b) are equal but only one is a TensorProductSpace | one constructor | proposed (R2 relies on it) | W1 |
| Direct sum of spaces | **twin**: CoupledSpace (N-ary) and FullFieldSpace (its N=2 case); 4 offset loops | one DirectSumSpace of labelled summands; the full field is bulk (+) trace | posing sequence, ruled 2026-09-28 (77 prototype gates) | W1 |
| Subspace on an index subset | **missing**: inflow/outflow halves and octants gathered by index | restriction of a factor is a space and a split pair | proposed | W1 |
| Dual of a factor | **missing**: DualSpace flips the whole space; a composite's dual is not a DualSpace | per-factor dual; the whole-space dual is the composition | R8, R14; O5 | W1 |
| Pairing | **partial**: Covector, InnerProductFunctional exist; posing._pair and the k estimator use np.sum | DualPairing: the only spelling of a contraction | R11, R15 | W1 |
| Boundary trace space | **twin**: 3 unrelated hierarchies (angular, scalar, radial half-angle), all leaves | the trace is a summand of the direct sum with its own frame (face nodes x half-range angular measure, Stokes weight); axis-built; measure includes face area | R4, R9; posing sequence, ruled 2026-09-28 | W1 |
| The four transport spaces | **missing**: one space for flux, source, importance and response | four variance patterns of one factor list | R14 (closes the posing sequence's open item 'flux and source as a grading') | W1 |
| Quotient by a metric kernel | **one**: FunctionSpace.rank, the kernel condition | unchanged | ruling 80 | - |
| **Elements** | | | | |
| Element of a space | **twin**: 4 member contracts (Vector, _is_ravellable, SystemField, CompositeField); 3 membership checks | one member contract; a member follows its space | proposed (R5 gives the verb) | W3 |
| Flux and source members | **partial**: RolePair forest (30 classes) and Member forest (32): a class per role and layout | role = variance pattern; layout = factor list | proposed; R5 and R14 point this way, while the field review found a missing type (flat against moment-carrying) | W3 |
| Importance and response members | **missing**: AdjointSolution.importance returns ScalarFlux | two member kinds on their own spaces | R17 | W3 |
| Co-vector | **one**: Covector, FunctionSpace.flat | unchanged (and used) | landed | - |
| Cross-section field | **partial**: SigS densified 3 times; chi and nu-Sigma_f gathered twice | a member of a named space, gathered once | #648 | W3 |
| Cell-moment factor and its member | **missing**: the cell average is lifted onto the moment axis by hand at 4 sites; 5 call sites tell flat from moment arrays by array rank | the within-cell moment is a typed factor; flat and moment-carrying members are distinct | R8 refined; #263 | W3 |
| Boundary source | **missing**: the inflow sampler returns a bare array; no adjoint-side boundary term | a typed member on the inflow trace, and its response partner | proposed | W7 |
| Positive cone | **one**: a predicate, by ruling | unchanged | landed | - |
| **Operators** | | | | |
| Operator with true ends | **partial**: 22 of 64 concrete operators fit the four spaces; stored ends (BoundOperator) live in transport, used by 15 of 80 | ends are stored in numerics and are the members' true spaces | R2 | W2 |
| Sum and product | **one**: OperatorSum, OperatorProduct | unchanged | landed | - |
| Tensor product of operators | **missing**: TensorProductOperator composes in sequence: (A & B) is bit-identical to (B @ A) | A (x) B on V (x) W with factor-indexed ends | R2 | W2 |
| Block operator on a direct sum | **weld**: CoupledOperator: block grid + triangular substitution + LU + assembly | one CoupledOperator keyed by (row, column) label; block substitution lives in one place, the plan | posing sequence, ruled 2026-09-28/29 | W2 |
| Transpose | **one**: LinearOperator.dual | unchanged | landed | - |
| Adjoint | **partial**: AdjointOperator dresses with whole-space Riesz legs and swaps ends | owned by the pairing: Riesz legs on the angular factor only | R11 | W2 |
| Inverse and its shortcuts | **partial**: 27 of 41 inverse/adjoint/dual methods return a stored object unchecked | a shortcut is taken where a predicate proves it | #705 | W2 |
| Natural maps of a product | **partial**: retraction, section, pullback exist on axis-built spaces; PartialCurrent and IsotropicEmission re-spell them on traces | one family on every factor | proposed | W2 |
| Restriction and extension | **twin**: trace gather/scatter, system restriction, BulkLift, and a second lift in MultiplicationOperator; no shared base | the pullback pair of a point map: one construction for retraction, section, decks and face identification | posing sequence, ruled 2026-09-28 (supersedes 'restriction siblings share no base') | W2 |
| Analysis and reconstruction | **one**: FrameAnalysis, FrameReconstruction | head declared as two spaces | R7 | W2 |
| Splitting A = M - N | **misplaced**: Splitting is SN-only; numerics passes M and N loose; no .H | split() of an Ordering: members in order plus the lagged coupling with its interface; (M, N) a derived view; adjoint = the reversed ordering | posing sequence, ruled 2026-09-29 | W2 |
| Green of a split | **weld**: inverse wrapper + splitting + driver loop + residual promise; dual() drops the realisation | the inverse of a Splitting; Krylov acceleration runs on the interface, so it is not a separate strategy type | posing sequence, ruled 2026-09-29; T12 | W2 |
| Coarsening (homogenise, condense) | **missing**: the channel-weighting rule is written once for energy and once for space | one verb, radon_nikodym: the density of one pushed-forward measure against another; the weight is state, each slot's kind typed | R16; posing sequence, ruled 2026-09-28 | W8b |
| Kernel pair and gauge of an operator | **twin**: PredictedKernel and NearKernel are two records of (ker A, ker A^T); the raising leg is hand-rolled in both; gauges built per method | Nullspace (the family is renamed); one object produced by composition, with its gauge | posing sequence, ruled 2026-09-28 (name); #681, #682 | W2 |
| Block and system role tags | **residue**: written and propagated; 0 production reads | retired; the block pattern is derived from tensor factors or declared per axis | posing sequence, ruled 2026-09-28 | W0 |
| **Transport physics** | | | | |
| Streaming | **partial**: declares full_field -> same; the coefficient has 2 homes (a 3rd flagged, #240) | flux -> source; one coefficient object per chart | R2, R6 | W4 |
| Multiplication by a field | **partial**: declares one space; re-does the lift | flux -> source; composed with the one lift | R2 | W4 |
| Group transfer, S = R Lambda M | **partial**: Lambda declares the head only; acts on head x group x cells | true ends; role tied to its reaction channel | R1, R2 | W4 |
| Fission production | **partial**: declares one space | flux -> source, rank-one in energy | R2 | W4 |
| The implicit inverse (sweep) | **weld**: one private body is three inverses: plain, windowed, scheduled | three operators; source -> flux | proposed | W4 |
| Sweep schedule and its walk | **twin**: the cell partial order lives as hub generators in 1-D and as a dependency graph for d >= 2; 16 walk classes; ScanMarch re-asks is_1d in 7 methods; a second copy of the octant epsilon | one schedule object (a partial order of cell visits) and one walk with a storage policy | proposed | W4 |
| Gradient and divergence | **missing**: LeakageOperator fuses Fick's map, the divergence and the P1 closure; stencil written twice | a typed pair under volume and area metrics | proposed | W4 |
| Wall law on its cone | **partial**: WallCondition is typed; the cone and its vector space are a free predicate and a mixin | a typed cone; one value per point | proposed | W7 |
| Realized boundary operator | **twin**: SN and diffusion boundary operators duplicated and diverged; the ray corner re-realises the specular law by hand | per-factor realization; one operator on the trace | proposed | W7 |
| Reaction rate | **partial**: typed for SN, diffusion, homogeneous; adjoint covector rebuilt by hand | one reading verb rate(cells, w) with the metric pairing; w = importance gives the adjoint-weighted rate | posing sequence, ruled; R14 | W5 |
| **Problem** | | | | |
| Problem | **partial**: 3 hubs, no shared type (MethodMesh covers 2) | a System (spaces, equation family, coordinates, cone, ports) on a Discretization; TransportMethod retires; the hubs are split member by member over the three layers | posing sequence, ruled 2026-09-27/28; #522 | W5 |
| Operator pencil | **partial**: typed LinearOperator[V, V]; diffusion mints it on the solver | (A, F): flux -> source, minted on the hub | R14 | W5 |
| Posing, and its adjoint | **partial**: the adjoint fixed-source posing is spelled by hand | the adjoint posing is the posing's dagger; at a pole a fixed source is its own type with the Fredholm refusal | posing sequence, ruled 2026-09-28; #484 | W5 |
| Question bound to a system | **missing**: Eigen, FixedSource, Response exist; nothing binds them to operators | EigenProblem, FixedSourceProblem, EvolutionProblem over the values Eigen, FixedSource, Evolution; pencil and spectral map derived from the parameter | posing sequence, ruled 2026-09-28; #529 | W5 |
| **Solution** | | | | |
| Strategy | **missing**: KEigenvalue fuses question and strategy; inner algorithm is a string at 3 sites | a PLAN (reads declarations and the block pattern only: partitions, orderings, cuts, splittings) and a BINDING (everything at a point: outer loop, spectral transformation, factorisations) | posing sequence, ruled 2026-09-28 | W5 |
| Eigen driver and estimator | **twin**: SN forward on SNSolver, adjoint on KEigenvalue; k estimator spelled 5 times; gauge applied twice | one outer loop in the binding; the estimator is rate(cells, w), which absorbs 9 production definitions | posing sequence, ruled 2026-09-28 | W5 |
| Iteration record | **partial**: IterationRecord derives its verdict; it carries two estimators of one contraction factor, and source iteration stops in flat l2 while recording in the space's l2 | one estimator, one norm | proposed | W5 |
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
| Method parameter records | **weld**: CPParams, MCParams, THParams, KineticsParams, FuelRodGeometry: one record each fusing problem data, discretisation choices and solver knobs | split onto Problem, generators and Strategy when each method is re-spelled | later, recycling | W9 |
| **Outside the algebra** | | | | |
| Traced memo and import tracer | **misplaced**: 69 classes in numerics/; import only stdlib, numpy and values | a package of their own below numerics | proposed | W0 |
| Code with no production caller | **residue**: 8 groups (quadrature selector, dense eigen engines, IntegralKernelOperator, solver_dae.py, ...) | retired, or moved to the test side | measured | W0 |
| Reference insulation | **partial**: production CP imports 3 derivations symbols, one private; the gate forbids only the other direction | no production import of a reference | ruling 2026-10-06; O14 | W0 |
| The concept map as an instrument | **missing**: conceptual_view.rst is hand-written; most 'implements' edges in the graph are inferred from names (the Nexus tool's own figure for this project: 12 999 of 13 084, not re-measured) | a tracked concept registry; a gate holds one class per concept; the page is generated | proposed | W0 |

## Reconciliation with the posing sequence (2026-11-13)

`posing_sequence.md` (#522, ontology ruled 2026-09-27 to 2026-09-29, units scheduled and not opened) already rules most of what this plan's first draft marked *proposed* for the Problem, the Strategy and the Solution. The user confirmed that the third exchange's points "were mostly identified previously". This plan therefore INHERITS that ontology as its target for those concepts and adds two things to it: the measured state of the tree against it, and the build route. Where the two plans name one thing differently, the posing sequence's ruled name is used.

| this plan's first draft | the posing sequence's ruled form | effect here |
|---|---|---|
| four phases: structure, posing, strategy, solution | three layers: 1 the System (spaces, equation family, discretisation, the system algebra), 2 the Problem (a question bound to a system), 3 the Solution, whose Strategy is a PLAN (reads declarations and the block pattern only; transfers to another question) and a BINDING (everything at a point) | the phase coordinate of T10 is: System, Problem, Plan, Binding and outcome |
| T12: the algebra names an inverse, a strategy chooses its realisation, a solver runs it | "a morphism is placed by what it READS": operator-level morphisms read the algebra, structure-level ones read the pattern; `(M, N)` is a derived view, `M` is never assembled to solve, block substitution lives in one place; `SourceIteration` and `KrylovAcceleration` weld plan and binding; acceleration runs on the interface | T12 is that ruling restated; the measured inversion (a realisation imports its driver) is its instance |
| a method-agnostic splitting with an adjoint | `Ordering` is the primitive (a `SpacePartition` plus a partial order of its blocks after declared cuts); `split()` returns the `Splitting`; the adjoint is the reversed ordering with the transposed cut | W2 and W5 build `Ordering`, not a lifted `sn.Splitting` |
| one direct sum | `DirectSumSpace` of labelled summands, 77 prototype gates | W1 |
| role tags have no reader | `BlockRole` and `SystemRole` retire; the block pattern is derived or declared per axis | W0 |
| "restriction siblings share no base", to reopen | restriction and extension are the pullback pair of a point map | former open question 5 is answered |
| coarsening as one morphism (R16) | one verb, `radon_nikodym`; the weight is state, never structure; χ in energy is a measure, pushed forward with no division | former open question 6 is answered; R16's "weight from the pairing partner" says WHICH density is passed |
| one Problem protocol | a System on a `Discretization`; `TransportMethod` retires; the hubs are split member by member over the three layers (the charter's first deliverable) | W5 |
| one estimator | `rate(cells, w)` with the metric pairing absorbs 9 production definitions | W5 |
| the flux and source spaces | `[OPEN]` there: "whether flux and source are a grading of the space" | closed by R14 (#707) |

Also inherited as names: `SpatialScheme`, `Nullspace` for the null-space sense of "kernel", `EigenProblem` / `FixedSourceProblem` / `EvolutionProblem`, `contraction_rate`, `relative_residual`, `Member` at the structure level. The gates the posing sequence owes ("Gates owed by the campaign", its lines 198 to 242) are part of each wave's verification specification.

**What this plan adds that the posing sequence does not have:** the census (which classes realise each concept today, and how well); the space algebra's defects below layer 1 (variance, the two realisations of a space, trace leaves, the tensor-product operator that is a product); geometry (W6); the boundary algebra (W7); the phase coordinate as a gate; and the route.

## The graph the code forms `[M]` 2026-11-13 at `c04f4ba8`

The user's hypothesis (second exchange): a well-organised code forms the graph of its ontology, a proper graph is acyclic, and acyclicity would demand each concept in its right place; and the suspected reason the graph is not acyclic today is that `numerics/` holds concepts of two later phases, strategy (the problem becomes an eigen or a source problem, cycles of unknowns are resolved, axes order a Gauss-Seidel) and solution (source iteration, Krylov, a matrix inverse).

**Instruments.**
- `scratch/ontology/deps.py`: an `ast` pass over the 260 production modules (`derivations/` excluded) that resolves every name a class body uses through the module's imports, including imports inside functions and under `TYPE_CHECKING`, to the class that defines it. 3 315 class-to-class references: 346 by inheritance, 1 776 through a module-level import or the same module, 149 through a late import, 1 044 typing-only. Positive controls: `FullFieldSpace` inherits `FunctionSpace`; `SNProblem` references `FunctionSpace`. It sees names; it cannot see a call through an object that was passed in.
- `scratch/ontology/nodes.py`: every one of the 678 production classes is assigned to one node, a (package, kind) pair with a phase. The assignment is module-path rules written by the orchestrator, with the census category as fallback; it is a draft, and the page shows each node's classes so that it can be corrected. 65 nodes, 246 node-to-node runtime edges.

**What was measured.**

| question | result |
|---|---|
| references that point up the layer order, late imports included | 0 |
| modules in a cycle (204 modules with a cross-module reference, 727 edges) | 0 |
| cycles among classes at runtime | 15, each inside one module: the operator algebra's closure (9 classes in `operator.py`), the memo's pins (8), sum types with their arms, `FunctionSpace` with `DualSpace` and `TensorProductSpace` |
| the same with typing-only references | 23, one of which spans 7 `sn` modules (12 classes): `SNProblem` and the operators it mints annotate each other |
| cycles among the 65 nodes (support nodes excluded) | 1: `numerics` Spaces, Operators, Members |
| references from a system-phase class (space, measure, member, operator) to a later phase, by import, whole tree | 0 of 409 between classes; 2 of 495 when module-level functions are attributed to their module's node |
| the same through objects passed in and named only in type hints | 43: fields and operators take the Problem (35: transport fields 14, sn operators 9, transport operators 8, diffusion operators 2, material fields 2); sn operators hold their sweep walk (8) |
| references against the phase order, whole tree | 5, all one edge: the Green's realisations (`Richardson`, `Krylov`) import the drivers they build (`SourceIteration`, `KrylovAcceleration`) |

The one concept cycle has four edges, and each is a placement the census had already flagged: `CoupledSpace` (a space) lives in the operator module `coupled_system.py`; `AngularTraceSpace` builds a `TraceRestrictionOperator`; `RankOneOperator` holds an `InnerProductFunctional`; `Field` uses `FunctionSpace`. Before two node splits the cycle had five nodes, and the splits were themselves findings: index sets and coefficient layouts sit below both spaces and bases, and "bases and frames" was two strata, bases below spaces and frames above operators.

**Verdict on the hypothesis.**
- `[REFUTED 2026-11-13]` **"The graph is not acyclic because numerics holds strategy and solution concepts."** The graph is acyclic at the module level, with zero references up the layer order and zero from structure to a later phase. The numerics-untangling campaign did that.
- **Confirmed, in a sharper form: a concept graph demands right placement where a module graph does not.** At the granularity of concepts there is one cycle, and its edges name misplaced classes exactly.
- **Acyclicity is necessary and far from sufficient.** The tree is acyclic and still carries 67 welds and 70 twins. A class that fuses two concepts contracts the arrow between them into itself, so no reference graph can see a weld: `GreenOperator` is an inverse (structure), a splitting (strategy) and a driver loop (solution) in one class, and shows as one clean node. The instrument that sees welds is the concept-to-class map, not the graph.
- **The phases are real and not encoded.** By import, the system phase never consumes a later one. Through injected objects it does, 43 times, and those are the user's suspicion confirmed at the transport and method layers: an operator that is handed the whole Problem, or that holds the walk that realises its inverse. (Positive control for the zero: `StreamingCollisionOperator` delegates to `LossRepresentation`; the instrument sees that edge, as typing-only.) The numerics layer shows one inversion. It is not structure: 22 numerics classes of posing, 4 of strategy and 19 of solution sit in one flat package beside 137 of structure; the strategy phase is nearly empty because its concepts are fused into operators and drivers (`GreenOperator`, `CoupledOperator`'s substitution, `KEigenvalue`, `SNSolver`, `DiffusionSolver`) or live at the wrong layer (`Splitting` in `sn/`).
- **A legitimate cycle exists and is not debt**: the closure of an algebra. `LinearOperator.__add__` returns an `OperatorSum`, which is a `LinearOperator`; a space's `dual()` returns a space. The law is therefore not "no cycles among classes" but "every cycle among classes lies inside one concept".

**Limits of these measurements.**
- The node graph is built from class bodies. Of the 3 315 references, 1 267 are owned by module-level functions; they are in the module-cycle count and in the "2 of 495" row, and nowhere else.
- T11's check ("every class cycle inside one concept") is CIRCULAR today: nodes are assigned by module path, and all 15 runtime cycles are intra-module, so they are inside one node by construction. It becomes a test when the registry assigns concepts independently of modules.
- Classes that fuse phases were placed by their module, so part of the zero is by construction; the concept-to-class map, not this graph, counts those.

**The cycles, classified** (the user asked whether the cycle can be undone, or whether spaces, maps and fields are one umbrella concept). Both, at different levels:

| kind | count of 15 | members | can it be undone? |
|---|---|---|---|
| a base whose operations return its own subclasses | 8 | `LinearOperator` with its adjoint, Riesz legs, identity, tensor product (9 classes); `FunctionSpace` with `DualSpace` and `TensorProductSpace`; a sum type with its arms (`Admission`, `WallCondition`, `RealizedBoundary`, `Manifold`/`Product`, `HilbertMetric`/`DualMetric`, `CountedRule`/`Refined`) | Only by giving up the notation: `A + B`, `A.H`, `V.dual()` would become free functions. This is the closure of an algebra, an inductive type and its constructors. One concept. |
| an adjoint pair | 1 | the axis retraction and its pullback | No: each is the other's adjoint by definition. One concept. |
| a space and its members | 1 | `CoupledSpace` and `CoupledField` | No, under R5 (a space mints its members; a member holds its space): a set and its points. One concept. |
| a factory and the thing that points back at it | 3 | the realizer table and its realizers; `ConcentricPartition`, `Chord`, `Rays`; `DiscreteMeasure` and its partition | `[R]` by the names, not yet read: probably yes, by building from the lower object without a reference back. To check in W0. |
| infrastructure | 2 | the traced memo's pins and jobs | outside the algebra |

Between the three concepts the answer is different. The node cycle Spaces, Operators, Members has four edges, and ONE closes it: a space that builds an operator (`AngularTraceSpace` uses `TraceRestrictionOperator`). Without it the order is spaces (with their members and co-vectors), then operators: a space is defined without any operator, a member needs its space, an operator needs spaces and may hold members. In mathematics the three are one notion, a category with its objects, morphisms and elements, and a dual space IS a space of maps; the code can still be well-founded because the metric is an object on the space and the Riesz legs are operators that read it (the Campaign 1 ruling), so no space is defined through an operator. So: one umbrella at the level of mathematics, three ordered concepts at the level of definition, and one legitimate cycle inside each.

**The laws** (accepted by the user, ruling O9; first proposed by the orchestrator):

- **T10, two coordinates.** Every concept has a layer (what it must know: values, geometry, mathematics, data, transport, one method) and a phase (structure: what is; posing: what is asked; strategy: how it is decomposed; solution: how it is run and what came out). A concept may consume only concepts at the same or a lower layer AND the same or an earlier phase. The layer axis is the existing criterion; the phase axis is new and is orthogonal to it: a splitting is mathematics, and so is a space, and knowledge alone cannot order them.
- **T11, every cycle among classes lies inside one concept.** `[REFUTED 2026-11-20]` as worded (2 of 2 class cycles on the pre-assigned classes span two concepts; on the first full assignment the concept quotient has one cycle of 58 concepts on references and one of 2 on inheritance); the re-wording is before the user, see W0. Checked on the concept registry: a strongly connected component of the class reference graph that spans two concepts is a weld, a misplacement, or two concepts that are one.
- **T12, the algebra names an inverse; a strategy chooses its realisation; a solver runs it.** `A.inverse()` is a morphism with ends and laws. Which algorithm realises it (a sweep, a factorisation, a relaxation, a Krylov method) is a strategy value bound later, and the loop is the solution phase's. This is what unwelds `GreenOperator` and `CoupledOperator`, reverses the one inverted arrow (a driver consumes a realisation, not the reverse), and is the earlier ruling "algebra eager, performance lazy" stated as structure.
- **Packaging, ruled (a).** (a) Sub-packages by phase inside each layer (`numerics/{structure,posing,strategy,solution}`, the same four in `transport/` and in each method), with the layer gate extended to the product order. (b) Sibling packages by phase beside `numerics/`. The orchestrator recommends (a): each axis keeps one criterion, the knowledge criterion still says where a concept lives, and the transport and method layers need the same four phases (`Splitting` and the Solution five-tuple are transport-layer strategy and solution).
- **The instrument becomes a gate.** The concept registry of W0 (class to concept, concept to layer and phase) carries three checks: one concept per class; edges monotone in both coordinates; every class cycle inside one concept. The page's graph is that registry's first draft.

**Effect on the chart.** W0 gains the registry with the two coordinates. W2 gains T12 (the inverse and its realisation apart) and the phase split of `numerics/`. W5 places `Splitting` and the Solution shape as transport-layer strategy and solution. The typing-only cycle in `sn` (the hub and its operators) is W4's: an operator takes the spaces and fields it needs, not the Problem.

## The open issues against the chart `[M]` 2026-11-13, run `wf_99d2fbf4-c0a`

13 agents, 2.15 M tokens, 5.5 min. All 490 open issues were read (body included) by 12 Sonnet shards and placed on the chart; an Opus reviewer then re-derived a sample. Data: `scratch/ontology/issues/all.json` (per issue: disposition, wave, concept, the reason, and the insight written out in full).

| wave | closes | closes part | carries an insight the wave must honour |
|---|---|---|---|
| W0 | 3 | 3 | 4 |
| W1 | 5 | 10 | 7 |
| W2 | 5 | 22 | 8 |
| W3 | 10 | 10 | 6 |
| W4 | 2 | 11 | 28 |
| W5 | 7 | 15 | 24 |
| W6 | 5 | 10 | 20 |
| W7 | 0 | 7 | 7 |
| W8a | 6 | 12 | 19 |
| W8b | 0 | 1 | 1 |
| W9 | 3 | 6 | 19 |
| all | 46 | 107 | 143 |

- 296 of 490 issues are on the chart; 194 are outside it: documentation, test infrastructure and the verification tree, the closed references, the harness, kinetics, fuel and thermal hydraulics.
- **"Closes" is an upper bound.** The reviewer re-derived 15 of the issues marked as closed and moved 5 to "closes part" (applied in the table); the other 36 were not re-derived. Of 25 issues marked outside, 1 belonged to a wave.
- 245 issues carry a written insight. Each wave's specification workflow receives its issues' insights as input. Examples the shards flagged: the predicted-shortcut kernel condition and its measurements (#705, W2); a nonlinear coupled low-order system that a Strategy type must admit beside a splitting (#318, W5); the trace metric's missing face area and what it moves (#680, W1); fission that is rank K over fissile isotopes, not rank one (#549, W4); an unowned octant order whose 25 patterns span a factor of 2.5 in iteration count (#343, W4); CP dropping its leakage term under vacuum (#512, #546, W9).

## The chart

Eleven waves. A wave is a set of units with one outcome; inside a wave the units are cut when its specification is written. The order below is the dependency order, not a schedule: W1, W6 and W8a depend only on W0 and run side by side after it.

| wave | outcome | after | workflows | carries | done when `[HYPOTHESIS]`, re-measured when the wave opens |
|---|---|---|---|---|---|
| **W0 Clear the ground** `[REMEDIED 2026-11-20]` this row's outcome and done-when are superseded by the section "W0, opened 2026-11-20" | Retire code with no caller; fix the surfaced defects; move the traced memo below numerics; make the concept map a tracked registry (class to concept; concept to layer and phase) with its three gates. | none | design-panel (the registry, O17), spec-check, build-units | posing units 0-1; #708-#711, #713, #641 | the registry gate passes its three checks on every production class; the 8 dead groups are gone and dead_references reads 0 for them; #708 to #711 closed |
| **W1 The space algebra** | Every factor one kind of object with variance and one metric home; one tensor product, one DirectSumSpace, subspaces; the pairing as the only contraction; axis-built trace summand; the four spaces. | W0 | design-panel, spec-check, build-units | #707 A2; posing unit 4 | AxisRetractionOperator.of(half-trace, 'angular') succeeds (it raises today); of_axes(a, b) and of_axes(a) * of_axes(b) have one type; FullFieldSpace is no longer a class; 0 pairings spelled as np.sum in posing and iteration; the non-uniform coarse-balance gate closes |
| **W2 Operators on true ends** | Stored ends in numerics; an operator tensor product; the pullback pair on every factor; the adjoint from the pairing; the algebra names an inverse and the plan realises it; predicted shortcuts; the phase sub-packages of numerics. | W1 | spec-check, build-units | #707 A4-A5; #705; posing unit 5 (plan half) | the operator map re-run reads 64 of 64 FITS (22 today); (A (x) B)(x (x) y) = Ax (x) By is a gate; the #705 negative control raises; 0 node edges against the phase order (1 today) |
| **W3 Members follow spaces** | One member contract; flux, source, importance and response as variance patterns; the cell-moment factor typed; cross-section fields gathered once. | W1 | design-panel, spec-check, build-units | #707 A3; #664, #633, #666 | AdjointSolution.importance returns an importance member; member-contract spellings 4 to 1; MomentField cannot be built with an L that disagrees with its space |
| **W4 Transport atoms** | Loss and gain re-bound flux to source; reaction role tied to its channel; schemes and closures as operators with derived transposes; the sweep's three inverses apart; one schedule object; gradient and divergence for diffusion; an operator takes spaces and fields, never the Problem. | W2, W3 | build-units (sequential slices on one branch) | #707 A6 and Steps B-D | the 118 Step A gates are green; 0 typing-only references from an operator or a field to a Problem class (35 today); 0 class cycles spanning modules with typing references (1 today); hand-written transposes in schemes and closures 0 |
| **W5 One System, Problem, Plan and Solution** | The posing sequence's layers 2 and 3 as built structure: Problem kinds over question values; Ordering, split() and the interface; plan and binding apart; one outer loop; one outcome; the Solution five-tuple in transport; observables as co-vectors. | W4 | spec-check (the ontology is ruled), build-units | posing units 5-6; #484, #529 | SN, diffusion and the homogeneous baseline each return the one Solution type; k-estimator spellings in those three 3 to 1; inner_solver string dispatch sites 3 to 0; forward and adjoint eigen solves share one driver; couple(split(A)) == A is a gate |
| **W6 Geometry** | Groups with translations; one orbit-space machinery; the coordinate system derived; a typed face; one mesh with region labels refining the region partition; one locator. | W0 | design-panel, build-units | boundary-law architecture E; #584, #569, #596 | RealSpace(3).quotient(a group with translations) returns a Quotient; coordinate-system tag tests in geometry 7 to 0; faces are a type in geometry; point locators in production 4 to 1 |
| **W7 Boundary algebra** | The wall-law cone typed; per-factor realization; one boundary operator for SN, diffusion and the ray corner; the adjoint side of a wall law. | W1, W2, W3, W4, W6 | spec-check, build-units |  | sub-Markov check sites 2 to 1; the SN and diffusion boundary operators are one realization; an adjoint boundary term exists and is gated; #653's 12 refused radial cases solve |
| **W8a Measures and generators** | A reference measure with mass and one identity; discrete measure typed by support; quadrature and basis unwelded; group structure read from the axis. | W0 | design-panel, build-units | posing unit 2; #624, #369 | typed 4 pi and 2 pi/n literals in bases, charts and rules 4 sites to 0; Lebesgue on [-1, 1] has one identity; the two G-space guards of #624 retire |
| **W8b Coarsening** | radon_nikodym: homogenisation and condensation as one verb with typed slots; provenance kept. | W3, W8a | spec-check, build-units | R16; posing unit 3 | Mixture.condense and MaterialXSField.project_through both call the one verb, bit-identical at every absorbed caller; the two-step defect identity is a gate |
| **W9 Later methods on the machinery** | Collision probability, characteristics and Monte Carlo re-spelled on the shared System, Problem, Plan and Solution; the ray operator designed once with its apply and sample forms; reference insulation restored. | W5, W6, W7 | one design-panel per method, build-units | sharpening order; the Monte Carlo seam | cp, moc and mc each define operator classes and return the Solution type; 0 production imports of derivations, enforced by a gate edge; #712 closed |

**How a wave runs, revised from the harness evaluation** (proposed by the orchestrator 2026-11-17; accepted by ruling O13, 2026-11-20; it supersedes the five steps kept below it where they differ). Every measurement cited is recorded in `harness_eval.md` on `main`; the stage table is the workflows page's and is not copied here.

1. *Specify.* A `design-panel` where a target is still proposed, then `spec-check` with its `units` argument. The specification is where a unit's defects are made: `[M]` on three units of the boundary carve, 30 of 36 review findings traced to the specification's text (14 it mandated, 16 where it was silent) and 1 to the builder deviating from it. So each unit's section states:
   - **its gate list**, each gate with its input and the mutation that reddens it. `[M]` With the list in its brief the gates stage transcribed it and added no gate that caught a known defect (0 of 10 runs), and a gate written from the plan cannot catch an error in the plan.
   - **every law whole** (with its converse or adjoint where it has one), **the invariants that relate fields**, the laws a constructor shares with its factories, whether a sum type is closed, and one address for each convention. `[M]` No builder and no gates writer supplied any of these (0 of 15 builds on six probed rows); they came from the review, or from the specification check of step 3.
   - **the raw inputs each new function parses and the invariants each new value states.** The builder's two standing obligations cover these when the specification is silent (`[M]` closed the difference to an Opus build on two units); naming them removes the dependence.
   - **everything the builder needs, in the text.** The build stage loads no CLAUDE.md, no rule and no skill, and the executor and checker stages hold no `Skill` tool and no Nexus tool. A module's layer and phase, the files it may touch, and any measurement that needs the graph are in the specification, taken before the workflow by the orchestrator or an explorer (three graph tools have no shell command, #717). A skill that a stage must read is named in its brief as an instruction ("your first tool call is to Read `.claude/skills/<name>/SKILL.md`"); `[M]` in an agent definition that form was followed in 5 of 5 runs and the plain form in 7 of 10, and it is not yet measured in a brief.
2. *Amend and user checkpoint.* The orchestrator applies the findings; the user rules on the specification. No ruling is made inside a workflow.
3. *Check the final text, per unit.* `spec-check` with `units` and `check_only: true`, after the LAST amendment, the checkpoint's and the orchestrator's own included; an amendment that adds a decision is checked again for the units it touches. `[M]` One elegance read and one structure read of a single unit's final sections named 66 %, 80 % and 68 % of the defects the later code review found, where the check of the whole first draft had caught 0 of the 14 mandated ones, and 3 of those 14 were check suggestions the orchestrator adopted without a re-read. Not yet run by real workflow agents (#715): W0 is its first use.
4. *Build.* `build-units` as it stands on `main`: gates by the test-architect on Sonnet, the build by the `builder` agent on Sonnet, review by qa and the elegance-enforcer on Opus in parallel, the amendment by the executor on Sonnet, the verification of each disposition by the checker on Opus. `qa_passes: 2` where a missed finding is expensive (`[M]` a second independent pass adds 0.80 to 1.25 of 8 known findings); the criterion is ruling O13.
5. *The orchestrator's final review* of every unit against its own standard of acceptance (O7), then integration and the full bar: the suite, pyright, Sphinx `-W`, `dead_references`, the layer gate, the registry gate.
6. *Close.* A `docs-pass`; the wave's issues closed or re-scoped with the measurement; the done-when predicates re-measured; the registry and the page updated. **And the wave's record:** every review finding classed by origin (the specification mandated it, was silent, or the builder deviated) and by kind (per field, relational, other), beside the defects that reached integration. A silent per-field finding that recurs is a candidate sentence for the builder's obligations, tested on the next unit; a mandated one says what the specification check missed. A difference between Sonnet and Opus at a stage is first read as an obligation the brief did not name (the evaluation's ruling H20), and Opus is given the stage only where naming it fails.

**What the evaluation did not measure, and how this plan treats it** `[R]`:
- **Unit kinds.** The build and gates measurements are three ADDITIVE units of one carve. A retirement, a move between packages and an in-place carve slice were not measured in arms. `[HYPOTHESIS]` The assignment transfers to them. W0 is made of those kinds, so W0's units are the record (step 6), read before W1's build is launched.
- **Stages.** The specification drafter, the design-panel drafter and its attackers, the elegance reviewer apart from qa, the amend stage and the verify stage were not measured. W1, W3, W6 and W8a open on a design panel.
- **The structure lens** added less than one finding of the rubric on 3 of 3 units; replacing it with a qa read is a lead read after the fact, not a result. The lens stays as it is.

**How a wave ran before the evaluation** (rulings O3, O7; the templates of `workflow_orchestration.md`; kept as the record, superseded by the revision above where they differ):
1. *Specify.* A `design-panel` where a target is still proposed, then a `spec-check`: one specification with every type, field, law and module path, the leaf units with their file sets, the dependency table, the wave's issues and their written insights (`scratch/ontology/issues/all.json`), and the gates. Sonnet drafts and measures; Opus attacks.
2. *User checkpoint* on the specification. No ruling is made inside a workflow.
3. *Build.* `build-units`: per unit, gates first and red, the build, Opus qa and elegance review in parallel, the amendment, the verification of each disposition.
4. *The orchestrator's final review* of every unit against its own standard of acceptance (O7), then integration and the full bar: the suite, pyright, Sphinx `-W`, `dead_references`, the layer gate, the registry gate.
5. *Close.* A `docs-pass`; the wave's issues closed or re-scoped with the measurement; the done-when predicates re-measured; the registry and the page updated.

**Branches.** `[REMEDIED 2026-11-20]` This paragraph is superseded by ruling O12, the starting point under it, and "Branch and integration mechanics" in the adversarial round below; it is kept as the record. The numerics lane (W1 to W5) continues the line of `refactor/true-ends`, which carries A1 and the 118 red Step A gates and stays off `main` until W4 turns them green. W0, W6 and W8a are cut from `main` and merge as they finish. This plan and its scratch data are needed in every worktree, so the plan is carried to `main` by a plan-only branch before the first workflow (`[skip ci]`). `[M]` 2026-11-17: this plan is not on `main` yet, and `refactor/true-ends` is 21 commits behind `main` (68 ahead): it has no `builder` agent, while `main`'s `build-units` dispatches one. So two preconditions of the first workflow: this plan reaches `main`, and the session that launches a workflow starts on a checkout that holds `main`'s harness (agent definitions are read when the session starts).

**Three lanes.** The numerics lane is W1, W2, W3, W4, W5: spaces, then operators and members, then the transport atoms, then the Problem and Solution shape. The geometry lane is W6 then W7. The generators lane is W8a then W8b. The lanes meet at W7 (the boundary needs trace spaces and typed faces) and at W9 (a later method needs all three).

**Why this order** `[R]`. Theme A of the map is upstream of themes B, C and E: an operator cannot declare flux to source until a factor carries variance; a boundary operator cannot be one object for discrete ordinates and diffusion until the trace is an axis-built space; the Solution shape cannot be shared until the members it holds are typed by their spaces. W0 goes first because it shrinks the surface every later wave touches and because a known defect is fixed before work is built on it (`process-discipline`).

**The campaigns in flight are placed, not superseded** (confirmed by ruling O16, 2026-11-20):

| campaign | state today | place on the chart |
|---|---|---|
| #707 true ends | Step A at A1; 118 gates red on `refactor/true-ends`; folded by O6 | its steps become the first slices of W1 (A2), W3 (A3), W2 (A4, A5) and W4 (A6, B to D). The red gates stay the gates. |
| #705 predicted shortcuts | `7615c845` parked on `fix/predicted-shortcuts`, waits on #707 | W2 |
| posing sequence #522, units 0 to 6 | ontology ruled; units scheduled, not opened | unit 2 in W8a, unit 3 in W8b, unit 4 in W1, unit 5 in W2 and W5, unit 6 in W5; units 0 and 1 in W0 |
| boundary-law ontology, architecture E | ruled; waits on #405 | W6 (the coordinate system as a group) and W7 |
| layer DAG, edge contracts | active, paused | every wave moves classes between layers, so the layer gate and the concept order are part of each wave's bar; edge contracts follow W5 |
| reference architecture #405 | P0 to P3 landed; P4 next | independent of the chart (references are closed); W9 restores the insulation |
| verification tree refinement | paused for #707 | resumes on a tree whose ontology has moved; its kill matrix is re-measured after W5 |

**Workflows.** `[REMEDIED 2026-11-16]` "as they stand" below was written on 2026-11-13; since then `build-units` gained the `builder` stage and `qa_passes`, and `spec-check` gained `units` and `check_only` (see the revised "How a wave runs"). The templates of `workflow_orchestration.md` are used as they stand: `design-panel` where a target is *proposed* and the ontology is still being searched; `spec-check` to turn ruled targets into one specification with its leaf units and dependency table; `build-units` for the gates-first build, review and amend in worktrees; `docs-pass` after each wave. A user checkpoint sits between workflows, never inside one. The census run is the first evidence for the scale: 33 agents, 0 errors, 7.06 M tokens, 40 min, and the Opus pass changed 577 of 1 221 claims (188 refuted, 389 amended), which is the case for keeping a reviewer on every executor stage.

**Sizing** `[R]`, as duration and never as a reason: each wave is one design or specification workflow, two to four build workflows and one documentation pass; W1, W4 and W5 are the largest. The full suite (about 75 min) runs once per integration. Added 2026-11-17: each wave also carries one check of the final text, two Opus reads per unit, before its gates stage.

## The adversarial round (2026-11-20)

The user: "I will give you a last round to adversarially attack parts of the plan, the ontology and try to foresee what is vague or badly formed." Three readers: the orchestrator (dependency closure, done-whens, the tree); a cross-domain attacker on the laws and the concept table (`scratch/ontology/attack/ontology_attack.md`, 14 findings); an Opus cold reader asked whether W0's specification can be written from this text alone (`scratch/ontology/attack/executability.md`; its verdict: no, 10 blocking questions). Neither report was itself reviewed by a second agent; the orchestrator re-read the cited lines of the findings it adopts below. The stale text they found is corrected in place above (the wave table's W0 and W7 rows, the R18 markers, the reference-insulation row, the "Branches" banner).

### A. Proposed to the user (each changes the ontology or a standing ruling); accepted by ruling O18, 2026-11-20

- **A1, a role is a side of the trial and test pair; a variance pattern is its chart.** The discrete balance is a Petrov-Galerkin statement `A : U → T*` with a trial space `U` and a test space `T`. The four roles are the four spaces of that pair: flux in `U`, source in `T*`, importance in `T`, response in `U*`, with the pairings ⟨importance, source⟩ and ⟨response, flux⟩ as R14 has them. Whether a factor's coordinates are stored on the function side or the integrated side is then DERIVED per factor from how the scheme tests that factor: an integral-tested factor (cells, groups, Galerkin moments) flips between the state and the data side; a node-tested factor (the discrete ordinates) does not, and its quadrature supplies the pairing. R14's four patterns are the case of the ordinate layout with `U = T`. T3 and the rows "Factor of a space", "The four transport spaces" and "Flux and source members" read the pattern as the role, and three witnesses separate the two:
  - the flux in the moment layout is dual on the angular factor (R7, `true_ends.md:32`) and primal on it in the ordinate layout, so role and layout are not independent coordinates;
  - the 0-D baseline has flux and source on one pattern (`true_ends.md:408`), so the pattern alone cannot tell them apart there;
  - for characteristics `U` (track segments) and `T` (flat-source regions) have different index sets, so the importance space is not the flux's factor list with flips.

  Two riders. R16's measure-or-function kind is a second datum per slot, independent of the side, and T3 does not carry it. On the boundary trace the rulings admit two readings of the angular side (R4 re-ruled against R9), decided by whether an inflow is imposed per ordinate or as a partial current per angular cell. Proposed: W1's design panel takes this as its premise (a System declares a trial factor list and a test factor list; the pattern is computed), with the trace reading as its first question.
- **A2, what "`derivations/` stays closed" (O8) permits.** `[M]` 2026-11-20 (`ast` over `orpheus/derivations/`, relative imports resolved; control: 180 `sn` to `numerics` names): it imports 157 names from production (geometry 51, data 49, numerics 24, specification 13, reference 11, values 5, mesh 4), and 3 of its files import the traced memo that W0 moves. Proposed: closed to redesign and to new capability; every unit re-points the `derivations/` imports its change breaks, mechanically, in the same unit; no compatibility shim is left in production (a shim is a twin).
- **A3, T12 against structure that exists only through a solve.** `DenseMetric` runs a pseudo-inverse with a cutoff when it is built (`numerics/metric.py:602`), and the quotient's rank is read from it, so a structural kernel is found numerically; the local linear-discontinuous solve and the orthonormal-span basis are the same shape (`[M]` the attacker: 7 executable `linalg` calls in 4 structure files). Proposed reading: T12 forbids a structure object that depends on a STRATEGY value (a tolerance, an iteration, a choice of algorithm); a finite solve local to one factor is structure when its answer is exact by algebra, and a kernel is declared by the generator that knows it, never detected with a cutoff. W1's panel checks every structure-phase solve against it. Collision probability's matrix, an assembled resolvent, is W9's instance.

### B. Decided by the orchestrator (O2), reported for objection

- **B1, the registry gate.** `[REFUTED 2026-11-20]` in part: see "What the panel and the censuses returned" under W0 (no reference kind is declared strict before the assignment is an input; "one class" is before the user). It is two gates, one per direction of the goal: every class names one concept; every concept names one class. Both are ratchets with a pending list, each entry naming the wave that removes it (the idiom of `tests/gates/numerics/test_concept_order.py`, which the registry absorbs, so that one instrument orders concepts). The first red of each is a class added with no entry and a concept given a second class. Monotone edges are checked on import and inheritance references (0 against the order today, so strict from W0); the 43 references through injected objects and the 5 against the phase order are pending entries of W4 and W2. Records, outcomes, certificates and refusals are entered under the sixth kind of T1 as amended. The registry is a tracked file, its inputs tracked with it (`plan-authoring` §10), and its name is settled in W0's panel because `numerics/registry.py` already owns the word. Concepts are assigned by the panel from the reviewed census and `tobe.py`, not by module path, which is what makes T11 a test.
- **B2, one phase scale.** The package names of ruling O9 stand (`structure`, `posing`, `strategy`, `solution`); the posing sequence's System, Problem, Plan, and Binding with outcome are the same four in order. "Layer" in this plan means the architecture's layer; the posing sequence's three layers are called tiers here.
- **B3, the defects W0 carries.** `[REMEDIED 2026-11-20]` for #708: it is a root fix in W0 (unit `angular-mass`), not a refusal tagged W3; #709, #710, #711 stand as written. #708 to #711 are fixed in W0 because work is not built on a known defect; where the root fix is a later wave's type (#709 in W3, #710 in W1, #711 in W8a, #708 in W3), W0's fix is a refusal tagged with that wave, and that wave's done-when is the refusal's retirement. W3's `MomentField` predicate is re-worded to that. The seven issues the census placed on W0 and the wave table omits (#605, #593, #392, #411, #561, #603, #607) are read into W0's specification.
- **B4, W0's retirements are ruled per group in its specification**, each against its latent consumers. `[M]` `IntegralKernelOperator` has no production caller and is named by two gate files, one a cross-check of the fission kernel, and W9's target for collision probability is "an integral operator": it is kept. The Sood cache lies inside `derivations/` and is out of scope under A2.
- **B5, O14's destinations.** `[REMEDIED 2026-11-20]` superseded by the specification's unit `insulation`: the composite rule is the direct sum of the panels' production rules pushed forward, the Bickley evaluator is a copy with the same interpolant, the chord is `ConcentricPartition.chord`. The three symbols production CP imports are a chord function, a composite Gauss-Legendre rule and a Bickley-function evaluator. W0's specification names a production home for each (the chord on the existing concentric partition, the rule in the quadrature package, the Bickley function as a special function in `numerics/`), and the gate edge forbidding production to import `derivations/` lands with the move.
- **B6, the traced memo.** Its package sits between `values` and `numerics`, beside `geometry`; `[M]` 0 geometry modules import it today, and the layer gate keeps that edge forbidden. Its consumers outside `numerics/` (3 `derivations/` files, the gate files, `tools/cone_store`) and the boot module's hard-coded module key are in the unit's file set; whether the move invalidates stored reference results is measured before the unit is cut.
- **B7, a unit's file set is its call-site set**: production, the existing tests that consume the changed interface, `derivations/` (A2) and `tools/`. A gate that is re-pointed is checked for demotion (`retirement-audit` D). The bar of step 5 gains: no pinned number moves, or it is re-baselined with its reason in the commit.
- **B8, the harness texts.** Step 6 of every wave updates what the wave made false in `CLAUDE.md`, the conceptual view, the `coding-elegance` skill and the agent definitions, through their sources and `python -m tools.harness`; reviewers load those texts, so a stale one argues for the old ontology.
- **B9, branch and integration mechanics (O12).** `refactor/ontology` is pushed and CI is read after every integration (the `gates` workflow runs on every push). A unit branch is cut from it, rebased onto it and fast-forwarded; integrations are serial, lanes included. `build-units` is called with `base: 'refactor/ontology'` and worktrees cut from that branch (its default is `main`). `main` receives only what the user asks for and is merged into the rewrite branch at a wave boundary. An issue is closed by hand with the integrating hash, because a `Closes` trailer acts only on the default branch. Every wave closes on a green bar, so each close is a state the user could promote.
- **B10, the lanes are specified side by side and are not file-disjoint.** W1, W6 and W8a meet in the axis, the measure and the trace. W1's panel therefore fixes two contracts first: what a factor reads from its generator (W8a builds behind it) and what the trace reads from a face (W6 builds behind it). W4 likewise reads the scheme's frame (W8a) and the chart's coefficients (W6) through those contracts. `spec-check` verifies that no two units in flight share a file.
- **B11, the Step A gates.** Each strict expected-failure marker is removed by the unit that turns its gate green, in that unit's diff; an unexpected pass is then the signal that a unit forgot one. The gates are an input to W1's panel, which may re-spell one with its reason (they pin today's spellings: `Axis`, `FlattenedOperator`, `CoupledField`).
- **B12, done-whens moved or re-worded.** "64 of 64 operators fit" belongs to W4, where the transport operators are re-bound; W2's is the numerics operators among them. W1 retires `FullFieldSpace` together with the members' re-pointing, or leaves the name until W3: the specification assigns every call site to one step (`plan-authoring` §6b). The k-estimator counts name their predicates when W5 opens (5 spellings of the estimator, 9 definitions the reading verb absorbs, 3 methods). The coarse-balance gate is the one `true_ends.md` specifies (its lines 385 and 458). The phase sub-packages of `transport/` and of each method, and the layer gate's extension to the product order, are W5's and W0's deliverables respectively.
- **B13, out of scope until decided.** `THParams`, `KineticsParams` and `FuelRodGeometry` leave W9's parameter-record row: no architectural investment is made in thermal hydraulics, kinetics or fuel until their fate is ruled.

### C. Inputs to a named wave's panel (not blocking W0)

- **W1.** The Riesz legs on a singular metric: `[M]` the attacker, the quotient leg is invariant to 2.6e-15 while the lift to the whole space and the extension outside the annihilator of the kernel are not (0.41 and 4.18). A derived norm does not exist on the energy factor (no reference measure: two refinement limits, 0.0891 and 0.0519) nor under an indefinite weight (#706, unruled). The trace is a direct sum over face orientations, each with its own angular measure `|Ω·n|w`, not one product of faces and a half range. The in-cell moment metric on a curvilinear cell is a fibre metric, not a product.
- **W2.** The natural-map family of T4 lacks the quotient maps, the permutation of factors and remaps between non-nested partitions. Four rows (natural maps, restriction and extension, subspace, coarsening) and the decks name one pullback pair of a point map. A derived transpose presupposes linearity (`[M]` 4 of 4 `is_linear` declarations are True, and the protocol admits nonlinear closures).
- **W4.** The sweep schedule and the posing sequence's `Ordering` are one object in two rows and two phases. "Three sweep operators" is one inverse, two realisations and one composite with moment ends. The fission row says rank one; #549 says rank K over fissile isotopes.
- **W5.** T5's "weight of an average read off the spaces" is re-worded: the space fixes the weight's kind, and its value is a member chosen by the functional to be preserved, passed as state (R16 and the posing sequence already say so). An observable and a response differ by the angular flat leg.
- **W7.** Its done-when solves #653's 12 refused radial cases, which un-defers part of the curvilinear build: the user rules on that at W7's checkpoint.
- **W8b.** One Radon-Nikodym verb covers flux-weighted pushforwards; the diffusion coefficient, current weighting and discontinuity factors are not of that form.
- **W9.** T7 holds for collision probability where scattering and fission factor through the flat-source projection; diffusion is an asymptotic limit of transport, not a realisation of its inverse; in Monte Carlo the estimator is a strategy choice and the Solution carries a covariance.

## W0, opened 2026-11-20 at `820500d8`

**Where it runs.** `[M]` `git diff --stat main -- .claude/ docs/development` shows only this plan, so this branch's harness is `main`'s; workflows are launched from the session that holds it. Scratch data is untracked and is passed to agents by absolute path under `/Users/rodrigo/git/nuclear/ORPHEUS/scratch/ontology/`; W0's own artefacts are under `scratch/ontology/w0/`.

**Units** (keys used by the specification; the file sets are cut from the censuses below; decided by the orchestrator under O2, `[R]` until the censuses return):

| key | outcome | settled by |
|---|---|---|
| `registry` | the tracked concept registry over every production class (not only the 83 concepts), its two ratchet gates (B1), the monotone-edge check, `test_concept_order.py` absorbed, the layer gate extended to the product order (B12) | design panel |
| `digest` | one digest object for the blake2b convention; the boot module is stdlib-only and runs by path before `orpheus` imports, which the design must respect | design panel |
| `memo` | the traced memo and import tracer in a package between `values` and `numerics` (O14, B6) | census |
| `insulation` | production homes for the three symbols CP imports from `derivations/`; the gate edge (O14, B5); #641 | census |
| `dead` | each group with no production caller retired, moved to the test side or kept, ruled per group (B4) | census |
| `roles` | `BlockRole` and `SystemRole` retired with no replacement in W0; the derived block pattern is W1's and W2's | census |
| `defects` | #708 to #711 as a root fix or a refusal tagged with the wave that retires it (B3); #713; the seven issues of B3 read in | census |
| `posing` | the items of posing units 0 and 1 that no later wave owns; u1 item 1 goes to W4, items 4 and 7 to W5; a u0 rename whose target a later wave re-spells goes to that wave | census |

**Done when** `[REFUTED 2026-11-20]` as written here: "0 import or inheritance references against the product order" and "both gates pass with a pending list" are re-stated after the assignment (see the results below). The text as it stood (it replaced the wave table's W0 cell, which predates B1 and B4): both registry gates pass with a pending list whose every entry names its wave, and each has its named first red; 0 import or inheritance references against the product order; every row of the dead-group table is executed and `dead_references` reads 0 for the retired names; 0 production imports of `derivations/`, held by a gate edge; the memo package exists and `numerics/` holds none of its modules; #708 to #711, #713 and #641 are closed by hand with the integrating hash (B9). `[HYPOTHESIS]` until the specification fixes each instrument.

**What the panel and the censuses returned (2026-11-20; run `wf_0bd23968-b0e`; files under `scratch/ontology/w0/`: `panel/design_{data,declared,graph}.md`, `panel/attack_{structure,elegance,instrument}.md`, `census_*.md`).** The orchestrator re-ran the probes named below; the rest is the agents' reading, reviewed by nobody else.

- `[REFUTED 2026-11-20]` **"The registry's gates have a first red today."** No class-to-concept assignment exists: `[M]` 103 of 678 production classes carry a `tobe` id on disk (the data design's `probe_g3_seed.py`); the census rows of the other 575 describe each class in prose and 637 of 678 rows embed a path. So the pending lists cannot be written, and the figures 0, 5 and 43 of "The graph the code forms" are readings on MODULE-PATH nodes, not on concepts. The assignment is a W0 step of its own, before the registry's specification: a closed vocabulary, an assigner and a reviewer per census unit, then the order violations measured on the result. A pending list is written into the specification from that table, each entry's wave taken from its concept's row; a list generated from the gate's own first run lets a failed W0 pass (the instrument attack, finding 1).
- `[REFUTED 2026-11-20]` **B1, "monotone edges are checked on import and inheritance references, strict from W0; the 5 against the phase order are pending entries of W2".** `[M]` `panel/probe_elegance.py`: 13 of 13 rows of `phase_up.json` are top-level imports and all 13 are legal under `test_concept_order.py`'s ranks. The two sentences cannot both hold. Until the assignment is measured, no reference kind is declared strict.
- `[REFUTED 2026-11-20]` **T11 as worded ("every class cycle inside one concept").** `[M]` `panel/probe_t11_quotient.py` on the 103 pre-assigned classes: 2 of 2 class cycles span two concepts (`FunctionSpace` with `TensorProductSpace`; `CoupledSpace` with `CoupledField`). `[M]` `panel/probe_closure_kind.py`: 8 of 23 cycles (typing references included) have a hub every member inherits, the algebra closures. Candidate re-wording, for the user at the checkpoint: the concept quotient is acyclic except inside one algebra closure.
- `[REFUTED 2026-11-20]` **"Every concept names one class", read literally.** `[M]` `panel/probe_initial_object.py`: `role`, `compose` and `faces` have no class the others inherit; the graph design counts 8 of the 10 rows marked "one" as more than one class. Candidates: one root, one head that may be a method or function, or one concept per class. For the user at the checkpoint, with the assignment table.
- `[REFUTED 2026-11-20]` **"8 dead groups".** `[M]` the census: 233 name-dead candidates of 2 927 production names; 28 rows and 8 modules read, 96 rejected, 113 classed by rule and unread. Only a row that was read enters the specification as a retirement.
- `[REFUTED 2026-11-20]` **#713, "`thermal_hydraulics/solver_dae.py` is a duplicate".** `[M]` two censuses: 88 differing lines, 6 of 14 definitions differ, `examples/` imports it. The twin is `THParams` and `THResult` only, and thermal hydraulics is out of scope (B13): W0 does not touch it; the issue is corrected with the measurement.
- **Orchestrator decisions from the censuses (O2).** The posing items that stay in W0 are u1.2, u1.3, u1.6 (inside `roles`) and u0.3 (docs residue); u0.4 sits inside `derivations/` and is excluded by O18. `memo` and `digest` are one integration, because each stales every stored memo entry (`[M]` 8 526 of 8 526 `DefPin` paths); its witness is that live entries regenerated after the move equal a snapshot taken before. The memo package is `orpheus.memo`. The edge forbidding production to import `derivations/` has its own table (`[M]` plain rows redden 95 test files). #709, #710 and #711 are refusals tagged W3, W1 and W8a; #392 is already fixed. #708 waits on which divisor is right on the 29 of 120 rules where the two disagree. The tier absorption of `test_concept_order.py` is decided by whether the concept order reproduces that gate's violation set; the module gate stays until it does.
- **For the user at the checkpoint** (each changes ruled text): T11's wording; what "one class" means; whether the traced memo is a seventh kind (`apparatus`) or a record; 83 concepts against 80 plus extension concepts for data, specification, reference and mesh; a phase value for what stands outside the four.

**The assignment pass (2026-11-20; run `wf_17f581bf-78d`; 14 Sonnet assigners, 14 Opus reviewers, one Opus close; files `scratch/ontology/w0/assign/`: `assignment.json`, `vocabulary.md`, `close.md`).** `[M]` 678 of 678 production classes carry one concept id, a kind, a role and a phase (`_close_merge.py`); the reviewers changed between 33 of 52 and 34 of 74 rows per unit on the three units read here, mostly an `arm` that does not inherit its head. It is a FIRST table, not yet a gate input:

- `[M]` `_close_measure.py`: 13 of the 83 concepts hold no class (`four`, `quotient`, `impresp`, `cellmoment`, `merge`: targets not built; the six method concepts; `insulation` and `registry`: work items, not concepts). 12 extension concepts were needed (reference solution 18 classes, verification certificate 10, and ten smaller). Of the 81 concepts with a row: 25 have no head (in 16 the defining class is a weld), 44 have one, 12 have two or more. 84 weld rows, 130 annex rows. 17 concepts carry two or more phases.
- **Eight concepts are used in two senses** (`close.md` section 3): `realized` (the realizer and the realized operator), `kernel` (the kernel pair and the gauge), `closure` (spaces and the closure generator), `hub`, `question`, `solution` against `outcome`, `measure`, `groups`. Generic operators are filed under physics ids (`collision` holds the identity and zero operators). The vocabulary is revised and the affected rows re-reviewed before the table is an input.
- `[M]` the orchestrator, `assign/_orch_order.py`, on the class reference graph `deps.json` (1 820 edges with both ends a class): **46 references run against the phase order** (35 top-level, 10 typing, 1 late) `[REFUTED 2026-11-20]` (the script keyed the 21 classes of `sn/loss_representation/__init__.py` with `__init__` in the key, which `deps.json` does not, so their edges dropped; with corrected keys, `scratch/ontology/w0/vocab/vocab_v2.py`, the figure is **49 of 1 968** judged edges: 35 top-level, 13 typing, 1 late; the concept cycle sizes are unchanged), and the **concept quotient has one cycle of 58 concepts** on import references (57 without annex rows; 17 and 3 without weld rows; 12 on heads and arms alone) and **one cycle of 2 on inheritance** (`fluxsrc` with `member`). So T11 read on concepts is far from holding on references and nearly holds on inheritance; the weld rows carry most of the cycle. The column `path_agrees` is unreliable (units used opposite conventions) and is not evidence.

**The specification, state on 2026-11-20 (a compaction point: re-anchor from here and `git log`).**

- **Two files.** Part A, `.claude/plans/ontology_w0_spec.md`: units `quadrature-ties` (2.0), `insulation` (2), `dead` in four sub-units (3), `roles` (4), `defects` (5.1 to 5.6), `angular-mass` (5.7), `posing` (6); its section 7 is the integration order and section 9 the record of every amendment pass. Part B, `.claude/plans/ontology_w0_spec_b.md`: unit `memo-digest`, integrating after all of part A. The `registry` unit has no specification yet and integrates last; it waits on the user's five vocabulary decisions listed above.
- **How it was made** (the record O13 asks for). Part A: three censuses, one Sonnet draft, a per-unit check of the first draft (10 Opus reads, 95 findings), four amendment passes, a check of the final text (14 reads, 183 findings, build-blocking in 7 of 7 units), a fifth pass, a second check (14 reads, 153 findings, 3 build-blocking, all gate mechanics; the reviewers RAN the units' code on a patched copy), a sixth pass, a targeted third check of the three places that pass added a decision (6 reads, 47 findings, 0 build-blocking), a seventh pass. Part B: one draft, a check (2 reads), an amendment, a check of the final text (2 reads, 29 findings, 0 build-blocking), a second amendment. Every amendment was applied by a Sonnet executor from a decisions file the orchestrator wrote (`scratch/ontology/w0/spec/**/amendment_decisions*.md`, `scratch/ontology/w0/specb/amendment_decisions_b.md`). `[M]` (the reviewers' returns in the three workflow journals, counted per reader, so one defect seen by both readers of a unit counts twice): build-blocking findings per check of part A's final text: 17 across 7 of 7 units, then 3 across 3 units, then 0 (the third check read 3 of the 7 units). NOT yet checked: part A's seventh pass and part B's second pass; one `check_only` run after the user's checkpoint covers both and the checkpoint's own amendments.
- **What the measurements changed** (each supersedes text above). #708 is a ROOT FIX in W0 (unit `angular-mass`: one definition of the angular mass, five of five production readers on it; the mixed source sum is retired, which also settles #288), not a refusal: `[M]` `scratch/ontology/w0/i708/memo.md`, the shipped `np.sum` divisor is the inexact one, 0 gates move. A production defect was found and is fixed first (unit `quadrature-ties`): `_certified_float` could not certify an exact value on a float64 rounding tie, so `LEGENDRE.on(a, b).gauss(n)` raised on tied panels; ties round half to even, and the closed-form recurrence builders are generic over the number type (`scratch/ontology/w0/oddn/memo.md`). All three factory cycles wait (W7, W6, W8a). The dead-code population is 51 retired, 42 moved to the test side, 49 kept (`scratch/ontology/w0/dead/reviewed.md`, 146 rows read; rows 134 to 137 and four moves that read private state left W0). The store's 6 043 "corrupt" entries are two retired manifest formats with intact payloads; the witness of the memo move is `scratch/ontology/w0/memo/witness.md`.
- **Left to do before a build workflow:** the user's checkpoint on both specifications; the issues section 8.1 of part A and section 9 of part B list are filed; the registry's vocabulary is revised on the user's decisions, the affected assignment rows re-reviewed, and its specification written. `build-units` is called with `base: 'refactor/ontology'`.

**Done when, restated 2026-11-20 from the specifications** (`[HYPOTHESIS]` until built; each predicate's command is in the unit's own done-when):
- `quadrature-ties`: `LEGENDRE.on(a, b).gauss(n)` raises on 0 of the 2 368 (panel, n) pairs of the unit's probe (it raises on 384 odd pairs of 2 176 today), and no rule that builds today moves.
- `insulation`: 0 production imports of `orpheus.derivations` (3 today), held by its own gate; #641 closed by hand.
- `dead` and `roles`: every `Retired`, `Moved`, `RetiredName` and `RetiredModule` row of the rows modules is green in the one retirement engine, and every `Kept` row is still defined; `dead_references` reads 0 for the retired names.
- `defects`: the three constructible inputs of #709, #710 and #711 are refused under `python -O` (all three construct today).
- `angular-mass`: 5 of 5 production readers of the angular mass read the section's mass (0 of 5 today; the acceleration code already does), and the mixed source sum is refused by `Field`'s own gate; #708 and #288 closed by hand.
- `posing`: one `measured_value` with `assert_never`; a complex gauge reading is refused.
- `memo-digest`: no memo module under `orpheus/numerics/`; one definition of the digest convention; the store witness is bit-equal entry by entry.
- `registry`: restated when its specification exists.

**qa passes per unit (O13's criterion: a unit that mints a type or states a law gets two).** Two: the retirement engine (built first and alone, inside unit `roles`, because six later units consume it and it exists only as text until then: its result is read before any other unit is dispatched), `quadrature-ties`, `defects`, `angular-mass`, `memo-digest`. One: `insulation`, the four `dead` sub-units, `posing`.

**The builds, state on 2026-11-20 evening (a compaction point: re-anchor from here, `git log` and `git worktree list`).**

- **The checkpoint is held** (O21 to O26). After it, both specifications were checked once more (8 Opus reads; two of them wrote the retirement engine from the text and ran it): 1 build-blocking finding (part B, `tools/cone_store/worker.py` reads the boot's `digest`), applied with the others as part A's eighth pass and part B's third (`559fbde7`). Decisions of that pass: `scratch/ontology/w0/spec/final4/amendment_decisions_pass8.md`. Part B's `_MEMO_SOURCES` is the static import closure of the memo with a pin completion, and the memo pin frames each file's bytes with its length; every entry of the live store regenerates once at integration.
- **The retirement engine is integrated** `[LANDED 5ec08eae]` on `refactor/ontology` (the build-units run, two qa passes, 21 findings answered, then the orchestrator's review: one ancestor rule on every leg). `[M]` the orchestrator, `scratch/ontology/w0/build/rows_check/run.py` (the built engine over the reviewers' rows modules on today's tree): 124 of 124 W0 rows red, 50 of 50 `Kept` entries green, 65 of 65 wave-2 rows green, no `GateError`. Part A's ninth pass (`954c4a4f`) writes what the engine does into section 1.5.
- **Running:** two `build-units` workflows from base `954c4a4f`, each unit in its own worktree on branch `refactor/ontology-w0-<unit>`: with two qa passes `quadrature-ties` (branch suffix `ties`), `defects`, `angular-mass`, `memo-digest`; with one `roles`, `dead-a` to `dead-d`, `posing`. `insulation` is dispatched only after `quadrature-ties` is integrated (1.1b). A worktree's memo store is its own (`CACHE` resolves from the tree's files), so the live store is not written by a build.
- **Integration, by the orchestrator, in the order of part A section 7:** ties, roles, insulation, dead-A to dead-D, defects, angular-mass, posing, then memo-digest (retake the memo snapshot immediately before it: `scratch/ontology/w0/memo/w10_retake_snapshot.sh`), registry last. Each unit gets the orchestrator's own review before its merge; after each merge that touches a hunk of 7.2, the unit's gate files and `tests/gates/test_retirements.py` are re-run on the integrated tip. Then the full bar, the docs build (it regenerates `docs/theory/verification/matrix.rst`, whose wave-2 test count is stale since the port), `dead_references`, and the issues closed by hand with the integrating hashes (#709, #710, #711, #392, #641, #708, #288). The DSA literal divisor `w / 2.0` is recorded on #520.
- **The registry waits on the user.** The vocabulary revision is `scratch/ontology/w0/vocab/` (`vocabulary_v2.md`, `conventions.md`, `moves.json`: 213 classes moved, 11 rows with a question). Four of its proposals change ruled text and were put to the user on 2026-11-20, unanswered at this point: T1's kinds widened to "element of a set" and "map between sets" with `generator` retired as a kind (the orchestrator's condition: member, co-vector and operator are still recorded where the set is a space); O23, `four` present through `RolePair` and `FieldRole`, and `dead` and `params` leaving as concepts; O22, a head may be a weld row, and `quotient` stays missing; O21, an algebra closure's hub must be a concept's defining class (5 closures, not 8). After the answers: the re-review of the 11 rows, the 63 arms and the 84 parts under the convention sheet (Sonnet assigns, Opus reviews), then the registry's specification.

**Order.** The panel (registry and digest) and the censuses run side by side; then the orchestrator's synthesis, `spec-check` with `units`, the amendment, `check_only` on the final text, and the user's checkpoint on the specification. No `build-units` before that checkpoint (step 2 of "How a wave runs").

## Open questions for the user

`[REMEDIED 2026-11-20]` Questions 1 to 5 below are answered by rulings O13 to O15; none is open.

Answered: R18 (O5), folding #707 (O6), the surgical-carve exception (O7), the references (O8), T10 to T12 and the packaging (O9). Answered by the posing sequence's rulings (see the reconciliation): the restriction family, and the two condensation verbs. The five below were open until 2026-11-20 and are kept as the record:

1. **The discrete measure's phase taxonomy** (angular, spatial, energy inside `numerics/`, ruled 2026-09-01). Keep, or re-rule in W8a?
2. **Production CP imports reference internals** (`cp/solver.py:49-51`, #641). Remove in W0, or accept until W9?
3. **Where the traced memo lives.** A package of its own between `values` and `numerics` is the lowest-knowledge home.
4. **Where two qa passes run** (added 2026-11-17). Proposed `[R]`: two for every unit that mints a type or states a law (the numerics lane W1 to W5, W6, W7, W8), because the relational findings are reached only by review; one for W0's moves and retirements, revisited from W0's record.
5. **The unit kinds the evaluation did not measure** (added 2026-11-17): retirements, moves and in-place carve slices. Proposed `[R]`: W0 runs on the measured assignment and is itself the record; the alternative is a controlled comparison (bare Sonnet, Sonnet with the obligations, Opus) on one W0 unit before the wave is built.

## Refuted candidates

- **"Read the ontology off graph degree."** `[REFUTED 2026-11-13]` for the question "what are the concepts": degree ranks hubs. The top production class by code edges is `SNProblem` (828), and the top node of the whole graph is a test helper. The fact it establishes: where a change costs most.
- **"The call edges pointing up the layer order are layering violations."** `[REFUTED 2026-11-13]`: they have no import edge beside them (0 import edges up the order among the core packages), so they are mis-resolved method names in the graph.
- **Census claims of missing machinery that the reviewers refuted because the object exists**: the Riesz legs, the n-ary direct sum (`CoupledSpace`), point evaluation of a basis, a Gram object, a space norm, a factorisation, the angular section with its mass (`AxisSectionOperator`), the measure on oriented lines (`LineDomain`), `Splitting` and `PredictedKernel` on the SN side, a Pin protocol in the memo. The fact each establishes: the debt is consumption or placement, not a new object. This is why the chart's most common status is *partial* (25 of 83) and not *missing* (17).
- **Candidate law T1 as first written ("five kinds and nothing else").** `[REFUTED 2026-11-13]` as a description of the tree: 78 production classes are refusals and 95 are records with no concept of their own, most of them result and outcome records. Amended: five kinds of mathematical object, plus values that record a run (records, outcomes, certificates, refusals), which are a sixth kind with its own home.

## Implementation starts when

`[REMEDIED 2026-11-20]` The user ruled the chart polished (O20); W0 is open. The condition as it stood: the user rules the chart polished; until then this campaign writes plans, measurements and maps only.
