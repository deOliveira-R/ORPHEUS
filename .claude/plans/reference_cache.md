# Reference solutions: generation, object, storage, and a cache keyed on the generator (#405, step 2)

Opened 2026-09-24. Status: **living plan, ontology being searched** (`plan-authoring` §0). Nothing is built. Implementation waits until the user rules the plan polished. The discussion starts after context compaction: *"Start drafting a living plan based on the context you have about this and prepare for context compaction. We will start this discussion after regrounding."* (the user, 2026-09-24).

## The architecture as it stands (consolidated 2026-09-24; the rulings and their reasons are in the discussions below)

1. **`Specification`** (name not final): the question only, never an answer. It holds the materials (`{id: Mixture}`), the geometry (a typed geometry-layer description that carries its boundary laws; absent for the infinite medium) and, for a source-driven case, a `Source` (a field protocol in `numerics`, W5 G: `RegionwiseConstant` per region and group, or `Symbolic`, whose SymPy realisation lives in derivations and is stored through `srepr`, for MMS). **The question is stored and typed** (W5 A, ruled 2026-09-25; `[REFUTED 2026-09-25]` "derived, never stored"): `Eigen(k | c | ...)`, `FixedSource`, `CriticalParameter`, with a forward/adjoint axis; `[REFUTED 2026-09-29]` for `CriticalParameter` and for the input-layer home: the posing sequence rules the question values physics-free, in `numerics`, as `Eigen(parameter, point, mode)` and `FixedSource(source, point)`, the parameter an opaque key to a system coordinate (`CellCoefficient | GeometryExtent | NuclideDensity`); see "The posing sequence's seeds for P1" below; the old derivation (a source means fixed source, none means eigenvalue) survives as the default constructor. It is plain data with content identity, in the geometry and data layers, so derivations build it and production consumes it; derivations never import `MaterialMesh` (the layer ban stays).
2. **Posing is a chain of overlays, each keeping or lowering the symmetry group:** materials, then geometry (finite support, region distribution, boundary laws), then mesh (discretisation per interval). The homogeneous case has no geometry overlay: one material and the total quotient.
3. **The test assembles** (option B): `MaterialMesh(spec.geometry.mesh(d), spec.materials)`, then the method's problem, which realises the geometry's boundary laws and projects the source onto its unknowns; the test chooses the refinement and the strategy. The geometry type owns `mesh(discretization) -> Mesh`; the discretisations act on the intervals of one axis (`CellsByCount`, `CellsByMaxWidth`, each with `.uniform(...)`, `2 * d` doubling), an axis product combines them for structured 2-D and 3-D; structured meshers only (unstructured meshes are imported; CSG is not meshed naively). The infinite medium is realised by the test: `HomogeneousProblem`, or a finite geometry from `from_homogeneous(width, boundary)`.
4. **Two kinds of solution of one specification:** `PublishedSolution` (printed values with printed precision and a `Citation(bibkey, locator)`; never recomputed) and `ReferenceSolution` (ours: produced by a generator; the answer stored as data: scalars with bounds, fields as per-region Chebyshev expansions with a tail-decay truncation bound, or symbolic forms for exact references). They share a protocol, `evaluate(functional) -> value with bound`, over `numerics.Functional` (W5 E, F); tests never read a payload's coefficients. **Every solution declares the equation it answers** (W5 B, ruled 2026-09-25): `(specification, approximations)`, empty for the continuous problem; exact solutions of a discrete equation (the S8-exact eigenvalue, flat-source CP) are admissible; comparisons and corroborations happen within one equation, or across a declared approximation with a certified bound. Fields are stored on panels graded toward singular sets, and a bound is refused unless the fitted tail decays geometrically (W5 E). `PublishedSolution` carries at least a `Withdrawn` state, for errata (W5 H).
5. **`ReferenceCertificate`** on each `ReferenceSolution`, per observable, generic in its subject so that operator references can be admitted later with little rework (W5 I, ruled 2026-09-25: none are built now); it is itself a traced memo of the whole certification run, the other solutions it corroborates included (W5 H): a bound; how it was established (`Exact` with a symbolic self-check; `ConvergedLadder` with the theoretical rate shown AND a structurally independent anchor; a comparison with a `PublishedSolution` is one anchor, required when one exists); corroborations with other solutions of the same specification, each with an independence note (a failed one makes both `Invalid`); applicability (the claim layers the reference can support); and a state: `Valid`, `Invalid` (automatic trip, latched under the generator's identity until code changes), `Withdrawn` (a committed tag with a reason and an open issue; the generator never runs). There is no admissible method for "a fine-mesh run of a production discretisation" (the user's reference ban); Monte Carlo is never a reference (#505).
6. **`VerificationCertificate`**: the return value of comparing a production `Solution`'s observable with a `Valid` `ReferenceSolution` of ours (never a `PublishedSolution`: the user's ruling of 2026-09-25), holding the values, the reference bound, the tolerance and the floor check (the reference's bound at most a tenth of the tolerance; for an order gate the SUT error at least ten times the bound). The test asserts on it; it is never stored. The validation side mirrors it later: `ExperimentalResult` and `ValidationCertificate`; code-to-code comparisons (against Monte Carlo) are L4, neither.
7. **One registry keyed on the specification**, queryable by facets; every test obtains references through it.
8. **The generator's identity is its recorded execution trace:** the hash of every function that ran (through `sys.monitoring`), the module-level definitions of their modules, the environment variables and files read, plus the specification's content, the settings, the trusted-library and Python versions, and a platform tag. The primitive is one traced memo whose entries record the child entries they consumed, validated recursively; raw `functools` memoisation is banned in `orpheus/` by a gate (W5 C: an `lru_cache` hit hides its dependency; 19 sites migrate). Imports are recorded from session start, the route-selecting environment is an explicit setting in the key, and the mpmath context is pinned with its precision in the key (W5 D). Persistent keys use a stable digest, never Python's salted `hash` (W5 F). Cache entries are self-validating (`.npz` and JSON, never pickle); payloads, certificates and intermediate stages (the Peierls volume kernel) are all clients of the traced memo. On CI a `references` job builds and validates, the shards read only, a scheduled cold rebuild controls the tracer, and a check keeps every `Withdrawn` tag pointing at an open issue.
9. **The Peierls Nyström solver is `Withdrawn`** (not research grade); Atkinson Nyström and PS-1982 keep running.

## The phases, revised after the W5 review (2026-09-25; supersedes the phase list of discussion 5, kept below as history)

Ruling (a) of 2026-09-25, the user: "Yes, I accept that correction (a)": "a continuous method we implemented" means ours, not published; an exact reference of a discrete equation is compared with production only when production solves the same discrete equation (the same scheme, quadrature and mesh), which B's within-one-equation rule already enforces.

Sizes are `[R]` estimates in sessions, sizing only. Each phase lands green on `main` with its gates carrying their first red (§6c), and is W1 or W3 by its nature (the carve phases P1 to P3 are surgical: the main agent writes, test-architect specifies the gates first).

- **P0, the withdrawal and the defects found (about 1).** `[LANDED 9187d998]` 2026-09-25; the landing record is "P0 landed" below. Recount the Peierls Nyström consumer set with `ps1982_reference` excluded (Atkinson is outside the package); a `withdrawn` marker (reason, issue) and a `conftest` hook that skips with the reason, plus an option to run them; the catalogue and V&V-matrix gates count a withdrawn test as neither catching nor verifying; replacement catchers for ERR-032 (a SymPy identity) and ERR-063 (a gate on the `xs_library` non-fissile chi values) `[REFUTED 2026-09-25]` for ERR-063: the value gate has no realizable red input, so ERR-063 went dormant by ruling (the P0 specification §5); the CP cylinder and sphere checks re-anchored on the flat-flux closed form; the present-tense-false text corrected (`summary.rst:56`, `escape_probability.rst:76`, `cases.py:225`, the `_cylinder_k_ref` docstring, `orpheus/derivations/README.md:107`, `.gitignore:9`, `layering.rst:61`, the stale whitelist entry). The marker's removal trigger is P4. Then re-run `test-durations`.
- **P1, the specification and the question (about 2).** In the geometry and data layers: the specification (materials, geometry as a value including the infinite medium, `Source`); the stored question sum type (`Eigen(k | c | ...)`, `FixedSource`, `CriticalParameter`, forward/adjoint) with the source-derived default; the `Source` field protocol in `numerics`; `geometry.mesh(discretization)` with `CellsByCount`, `CellsByMaxWidth` (the axis product designed in, 1-D built), `RegionMesh` retired (keeping #495's fix), `from_homogeneous`; one geometry-kind vocabulary and the production boundary types; a stable content digest for `Mixture` and the specification.
- **P2, the solutions and certificates (about 2).** The `evaluate(functional) -> value with bound` protocol over `numerics.Functional`; `PublishedSolution` (with `Withdrawn`) and `ReferenceSolution`, each declaring its equation `(specification, approximations)`; the graded-panel Chebyshev field representation with the geometric-tail refusal; `ReferenceCertificate[Subject]` (bound, establishment method, corroborations with independence notes, applicability, `Valid`/`Invalid`/`Withdrawn`), with the X4 decision on `Evidence`'s `NotYet`/`Certified`; `Citation(bibkey, locator)` and its `refs.bib` gate; the comparison verb returning a `VerificationCertificate` (a `Valid` `ReferenceSolution` of the same equation only; the tolerance-derived floor).
- **P3, the traced memo (about 2).** The recorder (`sys.monitoring` from session start, imports, the environment as an explicit setting, the pinned mpmath context); the traced-memo primitive with child entries validated recursively; the gate banning raw `functools` memoisation in `orpheus/` and the migration of its 19 sites; self-validating `.npz`/JSON entries under `.cache/`; the X1 witnesses (a mutation in a traced function or in a memoised child misses; a mutation elsewhere hits; flipping the route setting misses); the measurement of each generator under a foreign global `mp.dps`.
- **P4, the families, hottest first (about 3 to 5).** The registry keyed on the specification and question; per family: the generator returns a specification and a `ReferenceSolution` with its equation, its own convergence tests become its certificate run, and its consumers assemble `MaterialMesh(spec.geometry.mesh(d), spec.materials)`. Order: the cylinder multi-region Green's function, the rest of the trajectory resolvent, MMS (the source becomes `Symbolic`), the homogeneous, diffusion and CP `_CASES` (the flat-source CP cases as discrete-exact solutions), Sood (F_N certified against LA-13511 as `PublishedSolution` anchors, #305 re-scoped), then F_N, Case, Galerkin and the critical-size generators (`CriticalParameter`). Retired at the end: `ContinuousReferenceSolution`, `ProblemSpec`, `_CASES`, the lazy builders, `continuous_reference.Provenance` (the Sood registry's `Provenance` retired early, in P1 step 2c, into `Citation`), the Sood registry's `La13511Case` and `La13511Truth` (`sood_registry/case.py`; the Sood step splits them into `Specification` + `PublishedSolution`, and with them the `ELEGANCE-DEBT` #405 refusal of an empty `sources`), `CriticalSolution`, `FluxSolution`; the `withdrawn` markers become `Withdrawn` certificates.
- **P5, CI (about 1).** The `references` job with read-only shards; the scheduled cold rebuild (agreement within a fraction of the bound, not bit equality); the check that every `Withdrawn` tag names an open issue and an existing generator; `test-durations` again, #405's measurement of the result; the measurement of CPU-model variation within one runner label.

Compaction points: after P2 and after P4 (plan-authoring §6), each with the phase-to-commit table, the superseding corrections, the red baseline with gate costs, and the durable lessons.

## Ruled polished (2026-09-25); P0 preparation

The user, 2026-09-25: "The plan is polished. Make all preparations for implementation then prepare for compaction before we start." Implementation is authorised, phase by phase, in the order of "The phases, revised after the W5 review".

**The P0 withdrawal set, recounted `[M]` 2026-09-25** (the explorer's classifier `cls.py`, re-run with `solve_ps1982_vacuum_sphere` removed from its solver-symbol pattern, over the 44 files under `tests/` that name `peierls_nystrom` or a Peierls registry case; plus `cp/test_peierls_rank_n_protocol.py`, which reaches `solve_peierls_1g` through a subprocess worker the static classifier cannot see, `:568-582`). Result: 27 files; the only change from the explorer's 28 is `derivations/test_peierls_greens_function_xverif_ps1982.py`, which moves to primitive-only and keeps running. Population note (X2): the classifier's own tree gives 27 static solver files before the exclusion, the explorer's 28 counted `rank_n_protocol` by hand.

  - `tests/gates/cp/test_peierls_cylinder_flux.py`
  - `tests/gates/cp/test_peierls_flux.py`
  - `tests/gates/cp/test_peierls_rank_n_protocol.py`
  - `tests/gates/cp/test_peierls_sphere_flux.py`
  - `tests/gates/derivations/test_continuous_registry_lazy.py`
  - `tests/gates/derivations/test_peierls_assembly_drivers.py`
  - `tests/gates/derivations/test_peierls_closure_operator.py`
  - `tests/gates/derivations/test_peierls_convergence.py`
  - `tests/gates/derivations/test_peierls_cylinder_eigenvalue.py`
  - `tests/gates/derivations/test_peierls_cylinder_multi_region.py`
  - `tests/gates/derivations/test_peierls_cylinder_prefactor.py`
  - `tests/gates/derivations/test_peierls_cylinder_white_bc.py`
  - `tests/gates/derivations/test_peierls_fission_source_indexing.py`
  - `tests/gates/derivations/test_peierls_greens_function_slab_solver.py`
  - `tests/gates/derivations/test_peierls_greens_function_xverif.py`
  - `tests/gates/derivations/test_peierls_multigroup.py`
  - `tests/gates/derivations/test_peierls_nystrom_verification.py`
  - `tests/gates/derivations/test_peierls_rank2_bc.py`
  - `tests/gates/derivations/test_peierls_rank_n_bc.py`
  - `tests/gates/derivations/test_peierls_rank_n_class_b_mr_mg.py`
  - `tests/gates/derivations/test_peierls_rank_n_conservation.py`
  - `tests/gates/derivations/test_peierls_reference.py`
  - `tests/gates/derivations/test_peierls_reference_naming.py`
  - `tests/gates/derivations/test_peierls_specular_bc.py`
  - `tests/gates/derivations/test_peierls_sphere_eigenvalue.py`
  - `tests/gates/derivations/test_peierls_sphere_prefactor.py`
  - `tests/gates/derivations/test_peierls_sphere_white_bc.py`

**Granularity is per test, not per file `[R]`:** the explorer's per-test attribution (`pertest.py`) found 314 of the 499 cases in these files reaching a solver symbol; the other 185 (among them 4 of the 29 `peierls-equation` carriers) test primitives and keep running. The marker goes on the test functions or classes that reach the solver, or on the file where every test does. The test-architect's P0 specification decides the placement per file.

## Adversarial review W5, merged (2026-09-24)

Three independent first passes: the cross-domain attacker (mathematical and physical structure; report `attack_cross_domain.md`), the elegance enforcer (software machinery; `attack_elegance.md`, probes `ee_probe1.py`, `ee_probe2.py`), and the main agent's own list written before reading either (`attack_main.md`); all in `scratch/reference_architecture/`. Ranked by rewrite risk; each row says who found it.

| # | hole | evidence | risk | the change `[HYPOTHESIS]` |
|---|---|---|---|---|
| A | **The question cannot be derived from the source's presence** (all three) | `[M]` 25 of 47 Sood cases print a critical dimension; 32 `critical_dimension_mfp` values; `c_critical` at 4 sites; 4 `solve_critical` generators and `CriticalSolution` (88 hits, 12 files); an incident-flux boundary (Milne, albedo) drives a problem with no volume source; no adjoint spelling | rewrite of the P1 type | a stored, typed question sum type (`Eigen(k, c, ...)`, `FixedSource`, `CriticalParameter`, with a forward/adjoint axis), the ruled derivation kept as the default constructor; this REOPENS the discussion-4 ruling "derived, never stored" and needs the user |
| B | **A solution answers an equation, not only a specification** (cross-domain, main) | `[M]` `sn_slab_1eg_2rg_S8` is "the exact discrete-S_N eigenvalue" (`cases/sn.py:715`); 27 flat-source CP `_CASES` are exact for the discrete CP equations; Peierls white cases use rank-N closures (`geometry.py:4891`); `derivations/discrete/` has 12 modules, 4451 lines, 17 consumer test files; `[R]` a corroboration between an S8-exact and a continuous reference of one specification trips both correct references | rewrite after P4 | every solution declares `equation = (specification, approximations)`, empty for the continuous problem; compare and corroborate within one equation, or along a declared arrow with a certified bound; discrete-exact references are admissible (the ban is on fine-mesh runs, not exact discrete solutions); needs the user's scope ruling |
| C | **Memoisation blinds the trace** (elegance, main) | `[M]` probe: a second generator's trace is `{gen_B}`, the cached kernel and its helper missing; 18 `@cache`/`@lru_cache` sites under `orpheus/` (including `face_transmission.py:557,741`) and a hand-rolled memo (`half_range.py:88`); the plan's own intermediate-stage cache is an instance | rewrite of P3 if the trace is a flat set | one traced-memo primitive whose entries record the child entries they consumed, validated recursively; payloads, certificates and intermediate stages are its clients; a gate bans raw memoisation in `orpheus/` (19 sites migrate) |
| D | **State outside the trace** (elegance) | `[M]` `ORPHEUS_SLAB_VIA_E1` read at import (`cases.py:85`), invisible to audit hooks, so the plan's third X1 witness fails as designed; 15 module constants computed by first-party calls at import (`wims.py:148`); `dataclass`-generated code traced as `<string>`; 3 sites set the global `mp.dps` | additive | record imports from session start; the route-selecting environment enters the key (or the one import-time read becomes an explicit setting); pin the mpmath context and put `mp.prec` in the key |
| E | **The field bound is not a bound on non-smooth fields** (cross-domain, main) | `[M]` slab flux `2 - E2(x) - E2(a - x)`: true error over the tail-decay estimate = 4, 19, 78, 315 at n = 16 to 128 (the series converges as n^-2); graded panels reached 1.6e-8 with 288 coefficients; `[R]` the curvilinear angular flux is singular on the grazing curve; `[R]` a 2-D/3-D angular flux is about 270 MB per group per region | additive only if tests never read coefficients | the solution protocol is `evaluate(functional) -> value with bound`; panels graded toward singular sets; a bound is refused unless the fitted tail decays geometrically; angular fields stored only where a consumer needs them |
| F | **The observable duplicates `Functional`; the comparison has no single production type** (elegance) | `[M]` three unrelated production result types (`sn/solution.py:375`, `diffusion/solver.py:392`, `homogeneous/solver.py:77`); `numerics.Functional`, `ReactionRateFunctional` exist; `Mixture.__hash__` is salted per process | medium | the observable is one `Functional`-based object both sides evaluate; a stable digest of `Mixture._identity_key` for persistent keys |
| G | **`Symbolic` source at the input layer, but sympy is not a production dependency** (elegance) | `[M]` `pyproject.toml:16-23`: sympy only in the test and docs extras; 0 imports in geometry, data, transport, sn | medium | `Source` is a field protocol in `numerics` (L1); the sympy realisation lives in derivations |
| H | **The latch is keyed on the wrong identity; published solutions have no state** (elegance) | `[R]` after a failed corroboration, fixing one side leaves the other's key unchanged, so it stays `Invalid`; a cold-rebuild trip is a function of no key; the latch lives in an evictable cache; printed tables have errata | additive through C | the certificate is a traced memo of the whole certification run (its children include the other solution); "reproduce" means agreement within a fraction of the bound; `PublishedSolution` carries at least `Withdrawn`; the CI check also asserts a `Withdrawn` tag names an existing generator |
| I | **Operator references** (cross-domain, main) | `[M]` 4 production gate files consume the face-transmission operator reference; escape probabilities, K_vol, E_n are operators | additive | `ReferenceCertificate[Subject]`; the registry key not hard-typed to a problem specification |
| J | **Twins** (elegance) | `CriticalSolution`, `FluxSolution` (missing from P4's retirement list); `Evidence`'s `NotYet(issue, reason)` has `Withdrawn`'s shape and `Certified(bound, by)` a bound's; "Solution" names four things; the infinite medium as an absent geometry forces `is None` at consumers | additive | decide one type or two by X4; the infinite medium is a geometry value, not `None` |
| K | **Symmetry quotients; the physical specification** (cross-domain, main) | a reflective half-slab and the full slab are two specifications with one answer; ICSBEP describes a physical situation above the multigroup specification | additive | observables stated in coordinates so they pull back; `from_homogeneous` is the first covering morphism; a physical-situation record above the specification, where `ExperimentalResult` attaches |

Refuted attacks, with their reasons: C extensions, numba, ctypes (0 in `orpheus/`); `lambdify` (0 calls); sympy function-evaluation caching (0 subclasses); closures, partials, generators (`PY_START` fires regardless of caller); Withdrawn Peierls code reached by other generators (8 outside files name it in docstrings only); equal source in two modules, renames (the file is in the key: a miss, never a stale hit); randomness and threads (0); the layer placement of the new types (every needed edge is allowed); splitting `PublishedSolution` from `ReferenceSolution` (justified by capability); a separate latch mechanism (collapses into C); anisotropic scattering (`Mixture.SigS` is per Legendre order); surface sources, alpha and time (additive once A lands); tensor-train storage of a 2-D angular flux (the discontinuity rays destroy the low rank); a sheaf of certificates (the sup-norm bound already restricts to every functional).

Open measurements: whether each generator's output changes after `mp.dps = 50` is set globally (D); whether GitHub-hosted runners vary in CPU model within one label (H); a literature check of the grazing-ray singularity of the curvilinear angular flux (E).

### W5 rulings (2026-09-25)

The user: *"I agree with the findings. A and B get a yes ruling now."* So A (a stored, typed question) and B (every solution declares its equation; discrete-exact references in scope) are ruled; the corrections C to H and J are accepted with the findings. On I, verbatim: *"We're not going to immediatelly construct reference of operators. We can make the certificate generic enough to accept this with relatively little rework later."*

The user's own ruling, verbatim: *"A production solver can only be compared against a ReferenceSolution from a continuous method we implemented (such as Fn). If we want to use a publicated source, we need to implement the appropriate continuous solution method to generate the ReferenceSolution. The reason is because publications rarely (if ever) give the continuous answer, and many times they give a plot (without a numerical dataset). So no, we need to implement the method and generate our own ReferenceSolution."* `[REFUTED 2026-09-25]` my point 3 of discussion 5's third exchange, "its reference side is any certified solution [...] or a `PublishedSolution` (#305's direct Sood cross-checks)". A `PublishedSolution` is an anchor of a `ReferenceCertificate` only; the `VerificationCertificate`'s reference side is a `Valid` `ReferenceSolution` alone. #305 is re-scoped accordingly (comment of this date).

## The reframe, and the user's requirements (2026-09-24)

The user widened the scope before answering question 1, verbatim: *"We have identified that reference construction is the heavy load of the gates, and that is a clear target for cache storage. We will certainly continue to use analytical and semi-analytical references because there is no other way to formally verify code, so we will work on the cache for sure. Now let's take a step back and instead of simply working on cache, let's think about how to improve the reference solution generation, the solution object, the storage, the generation, etc (AND we will cache that). So start exploring this angle and its requirements."*

The requirements so far:

- **Q-R1, the problem is generated with the reference.** Verbatim: *"the resulting architecture must be ergonomic in such a way that the data that is required for the reference solution can seemlessly be used to generate the production test it is veryfying. [...] Suppose we're verifying ANY solver. The the generator of reference solution must be ergonomically capable of creating also a MaterialMesh (in other words, it can also generate up to the method independent step of production), and the test case, for example, an Sn test, does the final uplift of the MaterialMesh into SNProblem and selects solution strategy, etc."* This was the main requirement the Sood registry had to meet. The split it names: the reference's side owns everything up to the method-independent posing (the material mesh); the test owns the method head (the method's problem) and the solution strategy.
- **Q-R2, the Peierls Nyström references are not research grade.** Verbatim: *"the Nystrom cases are not research grade accuracy yet, therefore, they should not be used at this time (they also should not run, until we improve them to research grade accuracy)."* So they are neither consumed as references nor executed, until improved. Their quarantine is a deliverable of this plan; it removes most of the measured slow tier by itself (the worklist below: 6 of the 8 slowest files are Peierls).
- **The cache stays in scope** ("AND we will cache that"): the target state below still holds, now as one requirement among the redesign's.
- **Any aspect may be improved** ("you can suggest improvements in ANY aspect of the reference solution creation so go study it").

## The goal, in the project's terms

The user's ruling R1 of `.claude/plans/test_runtime_405.md` (2026-09-23), verbatim: *"We will first of all fight hard to eliminate slow verification. Maybe we reach the point where only regenerating semi-analytical solution is a slow process, which we can cache, and then they can be run only when the generator code changes. If we can't reach that point, we will think about how often it runs."*

So the target state is:

1. **A verification gate is fast.** It compares the system under test (SUT) against a reference it does not recompute.
2. **A semi-analytical reference is regenerated only when its generator changes**, where the generator is the code and the inputs that produce it. It is served from a cache between changes.
3. **A gate can never compare against a stale reference.** A stale reference is a reference the current generator would not produce. Cache poisoning must be unspellable, or loud; never silent.
4. **Whatever cost remains is the SUT's own cost**, and it is sized by what the claim needs (refinement rungs), not by habit.

Timing is never a gate (R6 of `.claude/plans/vv_suite_layout.md`). This plan changes what is recomputed, never what is verified.

## What is measured `[M]`

- **The cost, from `test-durations` run 35940553034** (2026-09-23/24; stamp ubuntu24, x86_64, 4 CPUs, Python 3.14.7, commit `62020ba5`; the full table and the method are in `.claude/plans/test_runtime_405.md`, "Step 1, measured"):
  - the whole gate tree, slow tier included, is 15.02 h serial over 12 491 cases;
  - the 30 slowest cases are 72.4% of it;
  - 11 of the 12 slowest files compute Peierls, trajectory-resolvent or Green's-function semi-analytical references;
  - 5 tests hit the 3000 s timeout, and 4 of them are reference builds (`test_continuous_registry_lazy.py` ×2, `test_peierls_rank2_bc.py`, `test_peierls_specular_bc.py[slab]`) plus cp's `test_peierls_flux.py` 2G2R convergence.
- **The routine (`-m "not slow"`) run, locally** (2026-09-23, host `.venv`, over transport, derivations, sn/operators and the layer-import gate): 40 min 43 s, about 38 min of it derivations. All 30 of its slowest calls are in `tests/gates/derivations/`, 26 of them in `test_peierls_*` files.
- **The reference object** is `ContinuousReferenceSolution` (`orpheus/derivations/common/continuous_reference.py:192`): a frozen dataclass whose `phi` and `psi` are callables that close over SymPy/mpmath state. As an object it is not serialisable as data; what a cache can hold is what the generator computed (for example Nyström nodes and a solution vector), from which the callable is rebuilt.
- **The registry** is `orpheus/derivations/reference_values.py`: two registries, the legacy `_CASES` and `_CONTINUOUS`. Expensive producers opt into lazy builders (`continuous_case_builders()`, #212; the only one today is `orpheus/derivations/continuous/peierls_nystrom/cases.py:438`), so a reference is built only when its name is requested, and memoised for the life of the process only.
- **Two earlier on-disk caches exist:**
  - `orpheus/derivations/continuous/sood_registry/cache.py`: a pickle per entry under `.cache/sood_registry/` (gitignored). Its key is `(solver qualified name, SHA-256 of the canonicalised kwargs)`, and its version defaults to the **git HEAD SHA**, so every commit invalidates every entry. That fails R1's "only when the generator changes" by construction. Consumers: its own package `__init__` and its test `tests/gates/derivations/test_sood_registry_cache.py` only (`git grep -ln "sood_cache\|SoodResultCache"`); no verification gate uses it.
  - `orpheus/derivations/_richardson_cache.json`: **retired** `[M]` 2026-09-23. Its producer `_richardson_cache.py` and the JSON were deleted whole at `91042339` (#290 P6, the diffusion island retirement); neither file exists. Two stale surfaces remain: the `.gitignore:9` entry, and `orpheus/derivations/README.md:107`, which says the utility "is retained as a generic reference cache" (present-tense false; fix with this plan's next commit). So the tree has exactly one on-disk reference cache today, the Sood cache.
- **In-process memoisation** exists in places (`functools.lru_cache` in `derivations/common/shifted_legendre.py`, `flat_source_cp/geometry.py`, `discrete/sn/dsa.py`, and `tests/gates/derivations/test_peierls_rank2_bc.py:976`); it helps within one run only.

## The worklist classified `[M]` 2026-09-23 (explorer; detail in `scratch/reference_architecture/` `worklist_classification.md`, per-test times from the run 35940553034 artifact)

Over the 8 slowest files (about 640 of the 901 runner minutes):

| file | where the time goes | generator | payload | rebuilds of one identical reference |
|---|---|---|---|---|
| `test_peierls_specular_bc.py` | the generator is the subject (slab ladders over N: 166 of 189 min) | `peierls_nystrom.geometry.solve_peierls_mg` / `_1g` | `PeierlsSolution`, frozen, plain arrays and `k_eff` | whole solve: 1 duplicate; **the volume kernel K_vol: 41 builds over 7 distinct inputs** |
| `test_l1_standoff_slab_cylinder.py` | cylinder: about 93% generator `[R]`; slab: the SUT, reference about 0 s | `trajectory_resolvent.greens_function_cylinder.solve_greens_function_cylinder_mr` | `CylinderGreensMRResult`, frozen arrays | **7 identical builds across 3 files** |
| `test_continuous_registry_lazy.py` | generator (both timeouts) | `peierls_nystrom.cases.continuous_case_builders` | `ContinuousReferenceSolution`, `phi` a closure over plain data | the test reads only `.name` |
| `test_peierls_rank2_bc.py` | the generator against a finer copy of itself | `build_volume_kernel`, `build_closure_operator` | `(N, N)` float64 | 0 |
| `cp/test_peierls_rank_n_protocol.py` | the generator against itself, and pinned literals | `solve_peierls_1g` in a subprocess | a scalar k | 0 |
| `cp/test_peierls_flux.py` | generator (the CP solve is trivial) | registry case `peierls_slab_2eg_2rg` | `phi_cell_average` | 3 builds across 2 files |
| `test_peierls_multigroup.py` | the generator is the subject (two routes compared) | `solve_peierls_mg` | `PeierlsSolution` | K_vol: 4 builds over 2 inputs |
| `test_phase_c_crosscheck.py` | generator (the other side is a literal or an `.npz`) | cylinder MR, as above | as above | 2 here |

What this changes:

1. **About 90% of the time is generators, and in 5 of the 8 files the generator is itself the subject of the test** (checked against a closed form, a sibling closure, a sibling route, or a finer copy of itself). For those five, caching the finished reference saves nothing: the test's claim is about the generator, so it must run the generator. The distinction table's third row is therefore the largest row, and it is wider than "a finer copy of itself".
2. **The reusable unit is below the finished reference.** The Peierls volume kernel K_vol depends on the geometry, the cross sections and the node set, not on the boundary condition or the closure rank. It is rebuilt on every solve and every group (`geometry.py:6115`, `_build_full_K_per_group`): 90 s locally for 8 nodes `[M]`, 41 builds over 7 inputs in the largest file alone. A cache keyed on K_vol's inputs speeds the generator's own tests too, which a whole-reference cache cannot.
3. **Whole-reference caching pays where one reference is shared**: the cylinder multi-region Green's function, 7 identical builds across 3 files.
4. **The generators are deterministic in their own code**: 0 hits for randomness, threads or process pools under `orpheus/derivations` (the same grep hits `orpheus/mc/solver.py`); mpmath precision is set by `workdps`. The float64 LAPACK stage is platform-dependent (the #504 class), so a cached payload is a per-platform fact unless its consumer compares with a tolerance.
5. **The generator's output depends on the environment**: `ORPHEUS_SLAB_VIA_E1`, read at import (`cases.py:86`), selects the slab route. A generator identity must include every input that selects the computation, environment variables among them.

Defects found on the way, to fix or file with the next commit:
- `cases.py:225` says the unified route is "Not the default" and names `ORPHEUS_SLAB_VIA_UNIFIED`; the code defaults to the unified route and reads `ORPHEUS_SLAB_VIA_E1` `[M]`.
- `test_l1_standoff_slab_cylinder.py`'s `_cylinder_k_ref` docstring says it is cached and takes about 30 s; it has no cache and takes 1000 to 1300 s on the runner.
- The shipped `peierls_slab_2eg_2rg` build (192 nodes at 30 digits) is extrapolated at about 14 h per group `[R]`; `test_continuous_registry_lazy.py` builds 13 references to compare names that now come from one function.

## The state of reference generation, measured 2026-09-24 `[M]` at `901f64ca`

Five parallel studies; their reports and re-runnable probe scripts are in `scratch/reference_architecture/` (`refgen_producers.md`, `refgen_consumers.md`, `refgen_nystrom.md`, `refgen_layers.md`, `refgen_literature.md`). The durable summaries are below; the explorer also filed two memory topics (`.claude/agent-memory/explorer/reference_producer_landscape.md`, `peierls_nystrom_blast_set.md`).

**The producers.**
- **The registry is the minority path.** Tests reach references through the registry about 100 times (`get` 76, `continuous_get` 22, `continuous_all_names` 4, in 38 files) and call generators directly at about 759 AST call sites (Peierls 230, trajectory resolvent 220, F_N 116, MMS 102, Case 49, Galerkin 26). A cache keyed on the registry alone misses most of the cost.
- **Three registries and one private dictionary:** `_CASES` (47 `VerificationCase`, all pickle); `_CONTINUOUS` eager (16 `ContinuousReferenceSolution`, 0 of 16 pickle: `phi` is a closure, and the MMS entries close over a case object holding functions and a `Quadrature`); 13 lazy Peierls builders; the Sood `La13511Case`/`La13511Truth` pair. The trajectory-resolvent, F_N, Case, Galerkin and PS-1982 families have no registry (14 `*Greens*Result` types and siblings, all plain data).
- **One concept, many spellings:** the problem description exists in 7 forms; the geometry kind is spelled 7 ways (`sphere-1d`, `sph1D`, `sphere`, `SPH`, ...); cross sections are passed 6 ways (`Mixture`, a raw dict, the `get_xs` dict, `sigma_t` arrays, `sig_t` arrays, a scalar `c`), with 14 `Mixture` builders in derivations and 10 in tests; two unrelated `Provenance` classes; the material-ID table written out 6 times; construction sites pass values outside `ProblemSpec`'s declared literals (`sphere`, `cylinder`, `pin-cell-2d`; boundary conditions `zero_flux` ×6, `white_rank2` ×2). `ProblemSpec.geometry_params` is a free dict with 11 key sets over 32 producer sites, one of which stores the generator itself (`"mms_case"`).
- **Unread fields:** 0 code readers of `ProblemSpec.boundary_conditions`, `geometry_type`, `as_verification_case`, or any `Provenance` field (upper bound: the AST pass does not resolve the receiver). The filters `by_geometry`, `by_groups`, `continuous_by_operator_form` have 0 callers. A reference is selectable only by name, in each registry's own name grammar.
- **The answer and the problem are mixed or missing:** `CylinderGreensMRResult` carries no problem; `PeierlsSolution` carries part of it (radius, geometry kind, group count), no cross sections, no boundary condition; `La13511Case.to_geometry()` reads the geometry size from the truth record, so posing the problem reads the answer.
- **Precedent for the requirement's split:** `MomentSpace`, `Billiard`, `Spectrum` and `BasisSpace` already take `(StructuredGeometry, dict[int, Mixture], method settings)`. **Counter-precedent:** all 13 MMS case types carry the method's discretisation (`quadrature`; `n_azi` and `ray_spacing` for MOC), and MOC MMS builds the method's own mesh on the reference side (`mms/moc.py:401`).

**The consumers** (296 gate files import derivations or name `reference_values`; 107 of them import only `xs_library`; the classes below are over the other 189, assigned by an AST and regex heuristic checked by hand on 12 files):

| how the test builds the production problem | files |
|---|---|
| from the reference's own data | 17 |
| cross sections from the reference, geometry rebuilt by hand from `geom_params` | 31 (regex-based; check per file before quoting) |
| the test holds the constants and passes them to both the generator and production | 29 |
| no production problem (98 of these are in `tests/gates/derivations/`) | 108 |
| other | 4 |

- 145 builder-shaped private helpers in 82 test files (`_slab_mesh` ×34, `_sphere_mesh` ×10, ...). The cylinder A|B|A problem is spelled in 3 files, which is the cause of the 7 identical builds.
- **The Sood registry meets Q-R1 on paper only.** `build_materials(case)` and `build_mesh(case, n_cells)` exist (a `Mesh1D`, single region only), but the production consumers are 2 smoke tests asserting a finite positive k; the real cross-check is open issue #305. `build_cp_params` places a method head on the reference side, and has 0 callers.
- **Refinement** is chosen by the test in most families (MMS `build_mesh(n_cells)`); two exceptions store it on the reference (flat-source CP: one cell per region; `sn_slab_1eg_2rg_S8` stores `n_ordinates`).
- **Posing today:** only diffusion poses through `MaterialMesh` in production (`diffusion/solver.py:467`); `solve_sn`, `solve_cp`, `solve_moc` take `(materials dict, Mesh1D, settings)`; `HomogeneousProblem` builds no `MaterialMesh`, by a gate. `from_material_mesh(` appears in 3 test files, all unit tests of the promotion.

**The layer contract** (`tests/gates/test_layer_imports.py:71`): `derivations` may import `numerics`, `geometry` and `data`; it may not import `transport` or any method package. The corrected premise: `[REFUTED 2026-09-24]` "derivations is below the input layer" (my earlier framing); derivations sits beside the input layer, below L2. `MaterialMesh` is at L2 (`orpheus/transport/mesh/material_mesh.py:122`), is already discretised, has content identity (`_identity_key`) and no spec form. So a reference can describe the problem in input-layer vocabulary (`StructuredGeometry`, `Mixture`, the production boundary-condition types) and never names `MaterialMesh`; the lift to `MaterialMesh` is a production verb at L2. Two doc defects: `layering.rst:61` contradicts its own linter, and the whitelist entry for `cases/diffusion.py` is stale (the file has 0 `orpheus.` imports).

**The Nyström blast set (Q-R2).**
- The solver half of `peierls_nystrom/` (the volume kernel, the closures, `solve_peierls_1g`/`_mg`, the slab eigen solve, the case builders, the 13 lazy registry builders, PS-1982) is consumed by 28 gate files: 499 cases, 314 of which reach a solver symbol, 562 of the 901 runner minutes. 12 more files import only the closed-form primitives from the same package (3.5 minutes) and are out of Q-R2's scope. `fn_method/peierls_atkinson_nystrom.py` is also Nyström (a scope question). Production imports: 0.
- **Claims that lose their only carrier:** `verifies` labels `flat-source`, `peierls-mg-operator`, `peierls-vacuum-bc-{cylinder,flux,slab,sphere}`, `peierls-white-bc`; nearly lost `peierls-equation` (29 to 4) and `hebert-3-323` (3 to 1). `catches`: ERR-027, 028, 029, 030, 032, 063 (6 of 88). The production claim at risk: CP slab 2-group 2-region flux (`cp/test_peierls_flux.py`) has no stand-in; the CP cylinder and sphere checks are homogeneous white-boundary cases where the exact flux is flat and k equals k_inf `[R]`, so a closed form can stand in.
- **No mechanism exists to hold tests out:** no quarantine marker, no default deselection. `test_error_catalogue_reconciles.py` reads test text, so a skip keeps it green and a move turns it red.
- **Present-tense false docs:** `summary.rst:56` ("30+ digit Nyström collocation") and `escape_probability.rst:76` ("production-grade").

**The literature** (read: Oberkampf & Trucano 2007, NEA 6298 §3.1-3.2; Briggs et al. 2003 on ICSBEP; ANL-7416 Supplement 2, benchmark 15; Ganapol 2008, NEA 6292; Mokhov, Mitchell & Peyton Jones 2020; Dolstra 2006. Unavailable: Oberkampf & Roy 2010, Knupp & Salari 2003, Roache; SAND2000-1444 is free but OSTI refused this host):
- **ANL-7416 keeps three records with their own identities:** the physical situation, the mathematical problem (with its expected results), and each solution (technique, program, results); one problem has many solutions.
- **Oberkampf & Trucano:** the mathematical problem "must not include any feature of the discretization"; the reference's accuracy is assessed per compared quantity; comparisons are not stored with the benchmark.
- **Accuracy grade:** Ganapol's "benchmark quality" is 4 to 5 correct digits, shown by a table of one output at requested errors 1e-3 to 1e-6; Oberkampf & Trucano: a reference is adequate when it can resolve the tested code's observed order of accuracy.
- **Caching:** a store keyed on hashes of the original inputs (Nix) needs deterministic tasks; a generator written in general Python cannot be compared by value, so the safe identity is its source closure; pytest-regressions keys files on the test name and cannot detect a stale reference.

## Candidate ontology v2 `[HYPOTHESIS]` 2026-09-24, for the discussion

Four objects, each with one role, and one registry.

1. **`Problem`: the mathematical problem, free of any method and any discretisation.** Geometry (typed, the input layer's vocabulary: `StructuredGeometry`, extended to 2-D and to the infinite medium), materials (`{id: Mixture}`), boundary conditions (the production types, not strings), the external source as a continuous function (for MMS: the manufactured source q(x, Ω, E)), and the question (eigenvalue or fixed source). Plain data with content identity; it lives at the input layer, so derivations can build it and production can consume it.
   - Makes unspellable: the 7 problem spellings; a reference that carries a method setting (MMS's `quadrature`); a problem posed by reading the answer; hand-rebuilt geometry in 60 test files.
   - The lift (Q-R1): `MaterialMesh.from_problem(problem, refinement)` at L2 (the production posing filtration gains its first stage), then the test's method head (`SNProblem.from_material_mesh(...)`) and its strategy. This moves production too: `solve_sn`/`solve_cp`/`solve_moc` would pose through the same stage.
2. **`ReferenceSolution`: one answer to one problem, by one generator.** The answer as data (k; the flux as a field representation evaluable as the functionals a test compares, for example cell averages over any cells the test chooses, so refinement stays the test's); the generator's identity and settings; its **qualification** (below). A problem has many reference solutions (Peierls and F_N for the same slab), per ANL-7416.
   - Makes unspellable: a reference that does not know what produced it; a reference whose flux is a closure over live state (the payload is data, the evaluation a view).
3. **`Qualification`: the reference's accuracy, computed, per observable.** The generator's own convergence study (its ladder of requested errors) produces it, as Ganapol's qualification table does. A gate states the grade it needs and refuses a reference below it.
   - This absorbs the table's third row: the generator's self-convergence tests become the computation of its qualification, keyed on the generator's identity like any other cached product, so they run exactly when the generator changes (question 5, dissolved rather than ruled `[R]`).
   - Q-R2 becomes a state of this object: the Nyström references are `withdrawn` (declared, not computed, since they must not run) until their improvement lands.
4. **The observable and the acceptance criterion stay in the test** (Oberkampf & Trucano: comparisons are not stored with the benchmark).

**The registry:** one, keyed on the problem (with a human name), queryable by facets (Ganapol's classification: geometry, groups, boundary condition, source, method, grade). Tests obtain every reference through it, so the cache sees every fetch (today 759 calls bypass the registry).

**The cache:** stores `ReferenceSolution` payloads and `Qualification` records, keyed on the problem's content identity plus the generator identity (the hash of its first-party source closure, its settings, the environment variables that select its route, and a platform tag where a float64 stage is platform-dependent). The same mechanism memoises a generator's intermediate stages (the Peierls volume kernel), so the generator's own tests benefit.

**What it leaves open:** where `Problem` lives (a new module in `geometry/`, or `StructuredGeometry` grown into it); how a continuous source and a field-valued answer are represented as data; whether `Case`, F_N, Green's-function results are re-typed or wrapped; the migration order.

## Discussion 1: who assembles the `MaterialMesh` (opened 2026-09-24)

The user, verbatim: *"we have to decide if the test sends the mesh requirements and the reference case returns a meterial mesh, or if the case sends a geometry and cross-sections and the test assembles the material mesh. Both are reasonable approaches, so we need a criteria to evaluate."* The user's angles: multiple geometry descriptions (structured, CSG, unstructured) favour the reference returning a geometry and the test assembling; "the sn test is testing Sn [...] somewhere in the reference files, in a single place, the material mesh is created in a consistent way" favours the reference assembling; tests of geometry and `MaterialMesh` construction exist to stress construction, not to supply consistent construction for parametric variation. The user's attack on `from_problem` (proposed name `.from_reference(ref_spec, refinement)`): over 1-D, 2-D, 3-D, Cartesian and curvilinear, structured, CSG and unstructured, it would become a huge function, or geometry creation would have to move down to `geometry/`.

Facts that bear on it `[M]` 2026-09-24:
- **The derivations import ban was not a user ruling as far as the record shows.** `"derivations": L2_PACKAGES | L3_PACKAGES` was introduced at `bc99ff6f` (2026-05-26, "P3.1 import-linter foundation test"), by the agent, from the layer assignment of `moment_space_and_layering_plan.md` §P3.0. The user's ruling on references is a different one (below).
- **The tree already has the geometry-to-mesh verb, on the mesh type, for 1-D.** `StructuredGeometry` (`orpheus/geometry/structured_geometry.py:215`) is method-free and mesh-free geometry with its boundary conditions; its docstring places it "between the registry (problem definition) and the mesh layer" and says "Reference solvers do not need a mesh and consume the StructuredGeometry". `Mesh1D.from_geometry(geometry, *, region_meshes)` (`orpheus/geometry/mesh.py:472`) is "the canonical geometry → mesh transition", with the refinement typed per region (`RegionMesh`: `n_cells`, `method`).

User rulings of this exchange:
- **R-Q1.2, R-Q1.3, R-Q1.4:** `ReferenceSolution`, `Qualification`, and the observable and acceptance criterion in the test: agreed (2026-09-24).
- **The reference ban, verbatim:** *"The only thing I forbade was using 'high resolution numerical cases as reference solutions' (in the same sense that Richardson extrapolation was retired as a verification method... because it doesn't fit the rules of formal verification)."* So a reference is analytical or semi-analytical with its own error estimate; a fine-mesh run of a production-grade discretisation is never a reference.
- **Import of `MaterialMesh` by derivations:** open; the user asks for a recommendation.

### Discussion 1, second exchange (2026-09-24)

The user's rulings and corrections, verbatim where quoted:
- **2-D and 3-D will exist; design for them now, build them later.** *"The only reason we don't have it yet is because we're focused on the ontology of the physics and math of neutron transport at the moment. Geometry and mesh creation are arguably a bookkeeping exercise that we can do later [...] But we should consider that they will exist in our designs."*
- **Naming: not `from_problem`.** "Problem" already names the fully posed method problem (`SNProblem`), and *"a used posing their own reactor would most likely not construct from a method that was made to help creating a geometry/mesh/materialmesh for verification."* The user's name stands: `from_reference(ref_spec, refinement)`. `[M]` 2026-09-24: `ReferenceSpec` and `from_reference` have 0 hits in `orpheus/`, `tests/`, `tools/`, and `from_reference` only this plan's hit in the prose corpus.
- **A correction to my reading of option A:** under A the reference computes its continuous solution without the `MaterialMesh`; the generator builds the `MaterialMesh` only to hand it to tests, *"because we know that the reference data will be used for verification. There is no point otherwise in having reference solutions in the first place."* So A need not duplicate the production builder (it can call the production verbs), and A's independence cost is that contamination becomes spellable, not that it happens. `[REFUTED 2026-09-24]` my criteria 1 and 2 as discriminators between A and B: both options can call the same production verbs, and neither reads mesh quantities in the answer path.
- **The infinite medium:** it would return a `HomogeneousProblem`, or a `MaterialMesh` with reflective boundaries; *"a infinite medium can be represented with both BCs as reflective, or periodic, so there is ambiguity in MaterialMesh creation."*

### Discussion 1, third exchange (2026-09-24): where the lift lives

**Ruled:** option B, the reference returns a specification and the test assembles (the user, 2026-09-24: "I agree with B"). The user's question: *"should from_reference be handled by the MaterialMesh then? seems like maybe this should be the job of the geometries (since we will have multiple ones. Or else MaterialMesh.from_reference would be a big method. Maybe it should be something like geometry.from_reference(spec.geometry) -> Mesh (refinement), MaterialMesh(mesh, spec.xs). [...] Attack this particular aspect of the design from multiple angles."*

Facts `[M]` 2026-09-24: `Mesh1D.from_geometry(geometry, *, region_meshes)` lives on the mesh type (`orpheus/geometry/mesh.py:472`); `git grep "\.from_geometry("` over `orpheus/` and `tests/` returns 123 lines, an upper bound (the spelling also matches `sn/sweep/cache.py`'s and `moc/geometry.py`'s own `from_geometry`). The refinement type `RegionMesh` (`mesh.py:206`: `n_cells`, and `method` as a string `"equal-volume"` or `"uniform"`) is one per region, so a refinement must match the geometry's region count. `StructuredGeometry` and the meshes share one package (`orpheus/geometry/`); `structured_geometry.py` already imports `mesh.BC`. `MaterialMesh(mesh, materials)` (`orpheus/transport/mesh/material_mesh.py:166`) is already the composition of a mesh with materials.

The analysis and the proposal are in the session reply of this date, and are summarised as `[HYPOTHESIS]` until ruled: no `from_reference` anywhere; the geometry type owns its meshing verb, `spec.geometry.mesh(refinement) -> Mesh`, dispatching once on the geometry's type; the test writes `MaterialMesh(spec.geometry.mesh(refinement), spec.materials)` and routes the boundary conditions, the source and the question into the method head itself.

### Discussion 1, fourth exchange (2026-09-24): the user's rulings on the lift

- **The boundary conditions belong to the geometry, not to the method's problem.** Verbatim: *"The boundary condition does not go to the methods problem. It is defined at the geometry. The method only realizes the law that was specified by the geometry."* `[REFUTED 2026-09-24]` my routing table that sent the boundary conditions into the method head; the method head receives the geometry's laws and realises them.
- **The posing ontology, verbatim:** *"we define materials, geometry is an overlay over materials that that makes the problem finite (boundaries) and defined the how materials are distributed in this finite space. A mesh in an overlay on geometry, and says how each region is discretized. This is also a point symmetry statement in the sense that once geometry is overlayed onto materials, it might belong to a certain symmetry group. When mesh is applied, the new symmetry group is at best the same, but might be lower. The homogeneous case is one where not even geometry is overlaid onto materials, instead we have a single material and total quotient."* So posing is a chain of overlays, materials, then geometry (finite support, region distribution, boundary laws), then mesh (discretisation per region), and each overlay can only keep or lower the symmetry group: the stabilisers form a descending chain.
- **The refinement is per region and not uniform.** Verbatim: *"a case will always know the regions of the geometry when doing refinement, because it had to know to specify the geometry in the first place. What in principle could be acceptable would be for geometry to have something like geometry.mesh(discretization) -> Mesh."* A single cell count works only for regions of equal size (a valid special case that deserves its own helper). The discretisation should be array-like, so that doubling every region is `2 * discretization`. `[REFUTED 2026-09-24]` my `CellsPerRegion(n)` / `MaxCellWidth(h)` with a SCALAR argument as the general form (a scalar was what I meant). The user's refinement of them (fifth exchange): the argument is array-like, one entry per region, and the uniform case is a named constructor, `CellsPerRegion.uniform(n)` and `MaxCellWidth.uniform(h)`; both are subtypes of one discretisation type passed to `geometry.mesh(...)`.
- **The homogeneous case has no mesh and no equivalent mesh.** It was designed that way on purpose. Acceptable instead: a geometry constructor that overlays a finite geometry onto one material, such as `geometry.from_homogeneous(width, boundary)`; the boundary argument resolves the reflective-or-periodic ambiguity, and the test chooses it.

The design as it now stands `[HYPOTHESIS]` until the user closes discussion 1:
- The reference returns a specification: the materials, the geometry (with its boundary laws; possibly absent, for the infinite medium), the source as a continuous function, and the question (eigenvalue or fixed source).
- The geometry type owns `mesh(discretization) -> Mesh`, dispatching once on the geometry's type; the discretisation carries a cell count per region, supports `2 * discretization`, and its length is checked against the geometry's regions inside `mesh`.
- The test writes `MaterialMesh(spec.geometry.mesh(d), spec.materials)`, then the method's problem, which realises the geometry's boundary laws and projects the source.
- There is no `from_reference`. For an infinite medium, the test either poses `HomogeneousProblem` from the materials or builds a finite geometry with `from_homogeneous(width, boundary)`.
- Open: whether the discretisation is a bare integer array or a small typed value that also carries the spacing law (today `RegionMesh.method`, a string).

### Discussion 1, fifth exchange (2026-09-24): the discretisation type, and 2-D and 3-D

User rulings, verbatim where quoted:
- `CellsPerRegion(array)` with `CellsPerRegion.uniform(n)` ("divide each region in n cells"); `MaxCellWidth(array)` with `MaxCellWidth.uniform(h)`; both subtypes of a discretisation type given to `geometry.mesh(...)`. Open, raised by the user: *"The only problem I see with this is handling of 2D and 3D geometries, so we need to think about that."*
- **The question is derived from the source, never stored:** *"If the reference spec returns a source, it is a source-driven case by definition. If it doesn't, then it is an eigenvalue problem."*

Facts `[M]` 2026-09-24: `Mesh2D` (`orpheus/geometry/mesh.py:619`) is a tensor product: `edges_x`, `edges_y`, a cell-level `mat_map`, a coordinate system (Cartesian x-y or cylindrical r-z), and four boundary fields `bc_xmin`...`bc_ymax`. There is no typed 2-D geometry above it; `factories.py` has `pwr_pin_2d` and the 1-D `_subdivide_zone`.

Proposal for 2-D and 3-D `[HYPOTHESIS]`, in the session reply of this date: the unit a structured discretisation acts on is the ZONE of an axis (an interval between two geometric breakpoints), and in 1-D a zone is a region. A structured multi-dimensional geometry is a product of axis partitions with a material map on the product of zones, so a region can span several zones and a zone boundary need not be a material boundary. Its discretisation is the product of one 1-D discretisation per axis (per-zone arrays), keyed by axis name. `MaxCellWidth` (a target size per region) is the family that extends to CSG and unstructured geometry; `CellsPerRegion` does not (an unstructured mesher cannot honour an exact count). Each geometry type states which discretisation types it accepts, in its `mesh` signature.

### Discussion 1, sixth exchange (2026-09-24): scope of the meshers

User rulings, verbatim where quoted:
- **Meshers are specified for structured geometries only, in 1-D to 3-D:** `CellsPerRegion`/`CellsPerZone` and `MaxCellWidth` "and simple meshers like these". Unstructured meshes are imported from an external mesher ("gmesh or openfoam or something else"); CSG is not meshed naively (*"for whatever is complicated enough that we might want an unstructured mesh to deal with, mesh quality becomes important, and we can't do it naively"*), with a possible exception for simple overlays such as detector tallies, and a possible route out through STEP export to an external mesher and back. `[REFUTED 2026-09-24]` for scope: my CSG and unstructured rows of the previous exchange.
- **Auto-registration of subclasses is allowed if it is the elegant design** (the user: "We can even use an auto registration system for subclasses if you find that to be an elegant design for this").
- **Zone versus region**, the user's question: are they synonymous? Answered in the session reply of this date: a region is the unit a material is assigned to; a zone is an interval of one axis between two geometric breakpoints; they coincide in 1-D and differ from 2-D up (a moderator region spans several zones; one x-zone crosses several materials). `[M]` the tree already uses "zone" for an axis interval to subdivide: `orpheus/geometry/factories.py:33` `_subdivide_zone`.

### Discussion 1, closed (2026-09-24): the ruled design

User rulings of the seventh exchange, verbatim: *"Let's simply call it CellsPerIntervals and that is generic enough. We can also name is as CellsByCount and CellsByMaxWidth understanding that they both work on an axis, and the axis product is necessary to split multiple axes."* On auto-registration: *"I was just thinking of how to have this in a GUI later (which we will have at some point for sure)."*

The ruled design of discussion 1:
1. **The reference returns a specification, and the test assembles** (option B).
2. **The specification** holds the materials, the geometry (with its boundary laws; absent for the infinite medium) and, for a source-driven case, the source as a continuous function. The question is derived: a source means a fixed-source problem, no source means an eigenvalue problem.
3. **Posing is a chain of overlays** (materials, then geometry, then mesh), and each overlay keeps or lowers the symmetry group. The homogeneous case has no geometry overlay: one material and the total quotient.
4. **The geometry owns `mesh(discretization) -> Mesh`.** The discretisations are axis-level, with one entry per interval of the axis: `CellsByCount` and `CellsByMaxWidth` (names adopted from the user's second option, since they name the pair symmetrically; `CellsPerIntervals` was the user's first), each with a `.uniform(...)` constructor, and `2 * d` doubling every interval. An axis product combines one per axis for 2-D and 3-D structured geometry. An interval is the stretch of one axis between two geometric breakpoints; in 1-D the intervals are the regions, from 2-D up a region can span several intervals and an interval can cross several materials. The spacing rule (today the string `RegionMesh.method`) moves into the discretisation.
5. **Structured meshers only, 1-D to 3-D.** Unstructured meshes are imported from external meshers; CSG is not meshed naively.
6. **The test writes** `MaterialMesh(spec.geometry.mesh(d), spec.materials)`, then the method's problem, which realises the geometry's boundary laws and projects the source. There is no `from_reference`.
7. **The infinite medium:** the test poses `HomogeneousProblem` from the materials, or builds a finite geometry with a constructor such as `from_homogeneous(width, boundary)`, choosing reflective or periodic.
8. **Auto-registration of discretisation types is deferred to its consumer:** a GUI (certain, later) or the cache's rebuild-from-storage (if it needs to find a type by name). The types are frozen dataclasses with typed fields, so a registry can be added later without changing them.
9. **The derivations import ban on `transport`** stays: under this design the derivations never need `MaterialMesh`, and the ban keeps the reference's answer from reading mesh quantities.

## Discussion 2: the identity of a generator (opened 2026-09-24)

The question: what must change for a cached reference, or a cached qualification, to be stale?

Facts `[M]` 2026-09-24 (a static AST closure over `import` and `from ... import` statements, relative imports resolved, each parent package's `__init__` counted because it executes; script `closure.py` in `scratch/reference_architecture/`; population 353 tracked `orpheus/**/*.py`):

| generator module | files in its first-party closure | lines | by package |
|---|---|---|---|
| `trajectory_resolvent.greens_function_cylinder` | 124 | 59 430 | numerics 53, geometry 21, data 14, derivations 35 |
| `peierls_nystrom.geometry` | 114 | 58 361 | numerics 53, geometry 21, data 14, derivations 25 |
| `discrete.sn.face_transmission` | 116 | 53 551 | numerics 53, geometry 21, data 14, derivations 27 |
| `sood_registry.la13511` | 113 | 53 913 | numerics 53, geometry 21, data 14, derivations 24 |
| `reference_values` | 112 | 51 729 | numerics 53, geometry 21, data 14, derivations 23 |

- **Every generator's module closure is about a third of the package, and nearly the same third.** The cause is the package `__init__` files: `orpheus/derivations/__init__.py` imports `reference_values`, whose closure is 112 files, so every module under `orpheus.derivations` inherits it; `orpheus/numerics/__init__.py` has 15 eager imports.
- **At runtime it is wider still:** the first registry lookup imports every public module under `orpheus.derivations` (`pkgutil.walk_packages`, `reference_values.py:183`).
- **So a module-level source hash is too coarse** `[R]` from the table: a change to any of about 88 shared files (numerics, geometry, data) invalidates every reference. It is narrower than git HEAD, and in the same failure class.

**Ruled 2026-09-24:** candidate 4, the recorded execution trace, is the generator identity (the user: "Option 4 is excellent"). Refuted for this question, with reasons: a module-level source-closure hash (too coarse, the table above); a hand-declared version (prose as enforcement, silently stale, X3); a static function-level call graph (unsound: Nexus's own measurement gives 12 to 15% recall of executed edges, so a missed dependency serves a stale reference in silence).

The identity, as ruled: the hash of every function that executed during generation, the module-level definitions of their modules, the environment variables and files read; plus the specification's content, the generator's settings, the trusted-library and Python versions, and a platform tag. Its X1 witnesses, owed before it lands: a mutation inside a recorded function turns a hit into a miss; a mutation in an unrelated function leaves a hit; flipping `ORPHEUS_SLAB_VIA_E1` turns a hit into a miss.

### Discussion 2, second exchange: the trace cache on GitHub Actions `[HYPOTHESIS]`

The user asked: "Is there a good way to make that work with the Git CI?" The proposal, in the session reply of this date:
- **An entry is self-validating.** It stores its lookup key (specification, generator entry point, settings, platform), its trace (each executed function's qualified name and source hash) and its payload. Validity is checked against the current checkout by re-hashing the traced functions, without running anything. So the CI cache key can be coarse (`refcache-<os>-<sha>`, restored by the prefix `refcache-<os>-`): restoring an old directory is safe, because a stale entry is a miss, never a wrong answer.
- **One job builds and validates, the shards only read**, as the GENDF tapes are handled in `test-durations.yml` today: a `references` job restores the latest cache, regenerates the misses, prunes the entries that no longer validate, and saves once; every test shard restores it read-only and generates any miss locally without saving.
- **Branch scoping** `[R]` from the Actions documentation, not yet read this session: a cache saved on `main` can be restored by other branches and by pull requests into `main`, and a branch cannot write into `main`'s scope. So `main` keeps the warm cache, and a branch regenerates only what its own edits invalidated.
- **The payload is data, never a pickle:** `.npz` arrays and JSON metadata, since loading a pickle executes code and the payload is plain data by the discussion-1 design.
- **A scheduled cold rebuild** regenerates every reference without the cache and compares with the cached entries: the X1 control for the tracer itself (a dependency the tracer cannot see, such as state inside a C extension, shows up as a mismatch).
- **Sizes are unmeasured** `[R]`: the Actions cache default quota is 10 GB per repository, with eviction after 7 days unused; a Peierls volume kernel at 192 nodes is about 0.3 MB per group.
- **Local development** keeps its own cache under `.cache/` (gitignored); the platform tag keeps it separate from CI's, since Apple Accelerate and x86 LAPACK differ (#504).

## Discussion 3: the qualification, the gate's requirement, and the Nyström withdrawal (opened 2026-09-24)

The user asked for the best proposal ("Make your best proposal"). The proposal `[HYPOTHESIS]` until ruled:

**What a `Qualification` holds, per observable** (k; each flux functional the reference can evaluate, such as cell averages over a given partition; currents):
- an **error bound**, absolute or relative, stated for that observable;
- **how the bound was established**, one of a closed set, and nothing else is admissible:
  - `Exact`: a closed form or an MMS manufactured solution, error bounded by the evaluation precision (float64, or the mpmath working precision);
  - `ConvergedLadder`: the generator's own resolution parameters (quadrature order, series terms, node count, working precision) refined over a ladder; admissible only when the ladder shows the generator's theoretical convergence rate (exponential for Gauss quadrature, the known algebraic order for a Nyström rule), and the bound comes from that rate, never from the difference of the last two rungs alone (Oberkampf & Trucano 2007, p. 83, against Ganapol 2008, p. 69); the ladder's table is stored as the evidence;
  - `Published`: digits as printed in a cited source (the Sood LA-13511 values), the bound being the printed precision.
  The user's ban on high-resolution numerical references is structural here: there is no member for "a fine-mesh run of a production discretisation".
- **the claim layers it can support**, from its pillar (`vv-principles`): an MMS reference supports convergence order and flux shape, never an eigenvalue;
- **a state:** `Qualified`, or `Withdrawn(reason, issue)`, declared and not computed, for a generator that must not run.

**What a gate requires:** no global grade label. The requirement is derived from the gate's own tolerance, at the comparison: one comparison verb (`assert_agrees(value, reference, observable, rtol)`, or similar) refuses, with a named error, any reference whose bound on that observable exceeds a tenth of the gate's tolerance. For an order-of-accuracy gate the same verb refuses a rung whose SUT error is not at least ten times the reference's bound (the floor of lesson L49). So "comparing against a reference too coarse for the claim" is unspellable; the factor ten is a proposal. "Research grade", the user's term, then has an operational meaning: a reference qualified by an admissible method with a bound below what every consumer gate asks for.

**Admission of a withdrawn generator:** it returns when its `ConvergedLadder` qualification is computed and meets its consumers' requirements. That is the criterion for lifting the Nyström withdrawal.

**The Nyström withdrawal, now:**
- **Scope:** the 28 gate files that reach the solver half of `peierls_nystrom` (499 cases, 562 of 901 runner minutes); the 12 files that import only its closed-form primitives keep running. Open for the user: `fn_method/peierls_atkinson_nystrom.py` (also Nyström) and the PS-1982 product-integration reference.
- **Mechanism, until the registry mediates every fetch:** a `withdrawn` marker on each of those files, carrying the reason and the issue; a `conftest` hook turns it into a skip that prints the reason, so the report shows every withdrawn test (never a silent deselection); an explicit option runs them on demand, for the improvement work.
- **The catalogue and the coverage claims stay truthful:** the error-catalogue gate and the V&V matrix count a withdrawn test as not catching and not verifying. So the claims that lose their only carrier show as open, not silently covered. The six error entries, one by one:
  - ERR-027, 028, 029 and 030 are Nyström assembly defects (K-matrix collocation, panel subdivision, rank-N normalisation); while the code does not run they are dormant, and their catchers return with it (status: dormant, withdrawn);
  - ERR-032 (a wrong ∫E₂ antiderivative) is an identity: a SymPy identity test catches it without Nyström, in milliseconds;
  - ERR-063 is a defect in shared cross-section data (non-fissile χ zeroed) that only a sink-χ-weighted consumer saw; a direct gate on the `xs_library` values catches it without Nyström.
- **The production claim at risk:** the CP `flat-source` verification. The cylinder and sphere cases (homogeneous, white boundary) move to the closed form (flat flux, k = k∞); the slab two-group two-region case has no stand-in, and the gap is filed as an issue and stated on the theory page's status, as the direction of development puts CP after SN, diffusion and homogeneous.
- **Docs:** the two present-tense false claims (`summary.rst:56` "30+ digit Nyström collocation", `escape_probability.rst:76` "production-grade") are corrected with the withdrawal.

### Discussion 3, second exchange (2026-09-24)

The user "generally agree[s]" with the proposal, and restated it for sharpening: running an analytical or semi-analytical generator that matches a reference from a report (Sood's, for example) yields a `ReferenceSolution` and a `Qualification`, the qualification carrying *"(1) error bound, (2) generator bounding method, (3) literature the value, (4) applicability and (5) a state"*. Two ideas, verbatim: *"(a) we might want to decide where the provenance goes (the cited source), generator or qualification. Also, I believe we do have a citation system in place using Sphinx and Latex. (b) We could identify triggers that would cause the Qualification to be withdrawn besides the Nystrom cases not demonstrating enough accuracy, because this seems akin to a circuit breaker to me (in the sense that we arm, and it stays there until somethings goes wrong, and then it disarms)."*

Facts `[M]` 2026-09-24: the citation system is `sphinxcontrib.bibtex` over one file, `docs/refs.bib` (78 entries; Zotero upstream through a Better-BibTeX export; `docs/conf.py:30-38`, #231). 15 files under `orpheus/` cite with `:cite:` in docstrings (for example `SoodLA13511_1999`, `DahlSjostrand1979`). The two `Provenance` classes hold citations as free text: `continuous_reference.Provenance(citation: str, derivation_notes, sympy_expression, precision_digits)`; `la13511.Provenance(paper_id: str, paper_table, primary_reference: str, notes)`.

The sharpened proposal `[HYPOTHESIS]`, in the session reply of this date:
- **(3) moves out of the qualification.** A printed value is itself a `ReferenceSolution` of the same problem, qualified `Published`; a generator matching it is a second solution. Their agreement is a **corroboration**, recorded on the qualification as evidence with an independence note (X4), not a field holding the value (ANL-7416: one problem, several solutions).
- **(a) Each citation goes where its claim lives**, as a typed `Citation(bibkey, locator)` whose key must exist in `docs/refs.bib` (a gate): the problem's source (the report that defined the benchmark) on the specification; a published value's source on that `ReferenceSolution`; the method's literature (the equation the generator implements) stays where it already is, the generator's docstring `:cite:` and the Nexus provenance chain, never copied into data. The two free-text `Provenance` classes retire into this.
- **(b) The breaker.** Triggers that trip a qualification: the generator's ladder no longer shows its theoretical rate; a corroboration fails (disagreement beyond the combined bounds trips both solutions, since the check cannot say which is wrong); the scheduled cold rebuild does not reproduce a cached payload; a catalogued defect is found in the generator (declared). Not triggers: a change of generator identity (that re-qualifies), and a gate asking for more accuracy than the bound (that is the gate's refusal). The latch comes free from the cache: a failed qualification is cached under the generator's identity and stays failed until the generator changes; the reset is a fix plus a passing re-qualification, or an explicit declared reset with a reason. A trip reddens the qualification gate, consumers of a tripped reference skip with the trip's reason, and restoring green on `main` means committing a `Withdrawn` declaration with its issue.

### Discussion 3, third exchange (2026-09-24): three states

The user, verbatim: *"Instead of 2 states, we can make it 3 states. Valid (a new name for qualified), Invalid (this is the trip, and it is based on pure triggers and automated QA system) and Withdrawn (this is the one that is declared, with reason and maybe a issue (which is probably going to be a GitHub issue that needs to be resolved). The advantage this gives us is that we can separate what will be skipped (indefinitely, until the GitHub issue is resolved) from what will be retested and quickly resolved. This is not very different than a circuit breaker that has a tag that explicitly and physically prevents trying to relatch until the physical tag is removed."* (Lockout-tagout.)

Accepted and sharpened `[HYPOTHESIS]` until ruled (session reply of this date):

| state | set by | the generator runs? | the qualification gate | consumer gates |
|---|---|---|---|---|
| `Valid` | a passing qualification, automatic | when its identity changes | green | compare |
| `Invalid` | a trigger, automatic | again whenever its identity changes (the retest) | red, naming the trigger and its evidence | skip, printing the trip |
| `Withdrawn` | a committed declaration: reason and an open GitHub issue | never, not even to qualify | skipped with the reason | skip, printing the reason and the issue |

Transitions: `Valid -> Invalid` (a trigger fires); `Invalid -> Valid` only through a change of code (the generator, or the other solution of a failed corroboration, or the trigger itself if it was wrong), because the cache latches an `Invalid` result under the generator's identity: this drops my earlier "declared reset", so no state is ever edited by hand except the tag; `Valid` or `Invalid -> Withdrawn` by a commit that adds the declaration; removing the declaration is the only way out of `Withdrawn`, and the next run qualifies the generator afresh. Policy: `Invalid` is transient (a branch's debugging loop); it must not stand on `main`, where the choice is fix or declare. A check that every `Withdrawn` declaration names an OPEN issue (a tag whose issue closed must be removed) needs the GitHub API, so it runs in CI, not in the local suite `[R]`.

### Discussion 3, closed (2026-09-24)

User rulings, verbatim:
- *"Instead of calling it a Qualification, we're going to be more explicit and call it VerificationCertificate. This is already foreseen that we will have a ValidationCertificate eventually."* Every `Qualification` in this plan's earlier sections now reads `VerificationCertificate`; the three states are `Valid`, `Invalid`, `Withdrawn`, as in the table above (ruled). `[M]` 2026-09-24 the name joins an existing family: `ExitCertificate` (`orpheus/numerics/outcome.py:114`, the solution's certificate, whose members are typed `Evidence = Measured | Certified | NotApplicable | NotYet`) and `OrbitCertificate` (`orpheus/numerics/invariance.py:282`). Whether the per-observable bound reuses that `Evidence` sum type is an implementation question to settle at design time (one concept, X4).
- *"Regarding the Atkinson Nustrom file and PS-1982, if they are research accurate (and I think they are), they keep running. Only the Peierls Nystrom is currently innacurate as I remember."* So the withdrawal is the Peierls Nyström solver alone. `[M]` 2026-09-24: `ps1982_reference.py` and `fn_method/peierls_atkinson_nystrom.py` import only numpy and scipy (`exp1`, `quad`, `CubicSpline`), not the Peierls Nyström solver; but `ps1982_reference.py` sits inside the `peierls_nystrom` package, and the explorer's 28-file set counted it in the solver half, so the withdrawal set is recomputed with it excluded before the markers are placed. The user's belief that both are research accurate `[R]` becomes a measurement when their certificates are first computed.

## Discussion 4: the specification's source and the reference's answer as data (opened 2026-09-24)

Proposal `[HYPOTHESIS]`, in the session reply of this date:
- **The source is a sum type:** `RegionwiseConstant` (a value per region and group; the ordinary fixed-source case) or `Symbolic` (a SymPy expression over declared coordinates, for MMS: derived symbolically from the manufactured solution, so structurally independent of production, and stored as data through `sympy.srepr`). Production projects either onto its unknowns; the projection is code under test.
- **The field-valued answer:** an exact reference keeps its symbolic form (stored the same way); a semi-analytical one is stored as a piecewise expansion, per region, in Chebyshev polynomials (scalar flux in x; angular flux in (x, mu) with mu split at 0, where the slab angular flux is discontinuous), coefficients per group. Point values and cell averages over any cells the test chooses are exact operations on the expansion (views, not data). The expansion's truncation error is estimated from the decay of its tail coefficients and enters the certificate's bound for every flux observable.
- **Scalars** (k, currents, leakage) are stored as numbers with their bounds.

### Discussion 4, closed (2026-09-24)

Ruled, verbatim: *"I agree with the source proposal. This is a better architecture and allows the source to be used for methods that do not have a quadrature in them. This is how it should have been the whole time. I also agree with the storage of the reference answer. It's an excellent proposal."*

## Discussion 5: the migration order (opened 2026-09-24; its phase list is superseded by "The phases, revised after the W5 review")

Proposal `[HYPOTHESIS]`, in the session reply of this date. Sizes are `[R]` estimates in sessions, stated as sizing and never as a reason (`plan-authoring` EFFORT-IS-SIZING). Each phase lands green on `main`, with its gates carrying their first red (§6c).

- **P0, the withdrawal and the defects found (about 1 session).** Recount the Peierls Nyström consumer set with `ps1982_reference` excluded; a `withdrawn` marker (reason, issue) on each file and a `conftest` hook that skips with the reason, plus an option that runs them; the catalogue and V&V-matrix gates count a withdrawn test as neither catching nor verifying; replacement catchers for ERR-032 (a SymPy identity) and ERR-063 (a gate on the `xs_library` non-fissile chi values); the CP cylinder and sphere checks re-anchored on the flat-flux closed form; issues filed for the Nyström improvement, the CP slab two-region gap, and the package-`__init__` breadth; the present-tense-false text corrected (`summary.rst:56`, `escape_probability.rst:76`, `cases.py:225`, the `_cylinder_k_ref` docstring, `orpheus/derivations/README.md:107`, `.gitignore:9`, `layering.rst:61` against its linter, the stale whitelist entry for `cases/diffusion.py`). The marker is transitional: its removal trigger is P4, where it becomes the certificate's `Withdrawn` state. Then re-run `test-durations`.
- **P1, the specification (about 1 to 2 sessions).** In the geometry and data layers: the specification type (materials, geometry, source; the question derived), the `Source` sum type (`RegionwiseConstant`, `Symbolic`), `geometry.mesh(discretization)` with `CellsByCount` and `CellsByMaxWidth` (1-D built, the axis product designed in), `RegionMesh` retired into them, `from_homogeneous`, one geometry-kind vocabulary and the production boundary-condition types in place of the 7 spellings and the string vocabulary.
- **P2, the answer and its certificate (about 1 to 2 sessions).** `ReferenceSolution` (answer as data; the Chebyshev field representation with its tail-decay bound; scalars with bounds), `VerificationCertificate` (per-observable bound, the three establishment methods, corroborations, applicability, the three states), `Citation(bibkey, locator)` with its `refs.bib` gate, and the comparison verb with its tolerance-derived floor.
- **P3, the identity and the cache (about 1 session).** The `sys.monitoring` trace recorder, self-validating entries in `.npz` and JSON under `.cache/`, intermediate-stage memoisation, and the three X1 witnesses of discussion 2.
- **P4, the families, one at a time, hottest first (about 3 to 5 sessions).** The single registry keyed on the specification; each family's generator returns a specification and a `ReferenceSolution`, its own convergence tests become its certificate computation, and its consumers move to `MaterialMesh(spec.geometry.mesh(d), spec.materials)`. Order: the cylinder multi-region Green's function (7 identical builds in 3 files), the rest of the trajectory-resolvent family, MMS (the source becomes `Symbolic`, the quadrature leaves the case), the homogeneous, diffusion and CP `_CASES`, Sood (its `Provenance` retires), then F_N, Case and Galerkin. CP and MoC consumers take `spec.geometry.mesh(d)` and `spec.materials` as they take a mesh and materials today; no new CP or MoC machinery (the direction of development). At the end `ContinuousReferenceSolution`, `ProblemSpec`, `_CASES`, the lazy builders and both `Provenance` classes retire, and the `withdrawn` markers become `Withdrawn` certificates.
- **P5, CI (about 1 session).** The `references` job with read-only shards, the scheduled cold rebuild, the open-issue check for `Withdrawn` declarations; then `test-durations` again, which is #405's measurement of the result.

### Discussion 5, second exchange (2026-09-24): existing issues, and the solution objects

The user asked to check GitHub for existing issues before filing any, and to elaborate the solution objects and how comparisons make a certificate valid.

**The issue map** `[M]` 2026-09-24 (`gh issue list --state open --search "<term> in:title,body"` over 20 terms; 308 open issues; raw output `issues_survey.txt` in `scratch/reference_architecture/`; the search matches bodies, so each row below was judged by its title and, for the six marked *read*, its body):
- **This campaign:** #405 (the umbrella; this plan is its step 2). #211 (SN suite reference selection, caching, parallelism) is largely subsumed.
- **Absorbed by the design, to link now and close as each phase lands:** #145 (*read*: a per-reference precision-floor sweep and a right-sizer; it is the `ConvergedLadder` certificate, P2 and P4); #356 (*read*: the Peierls cylinder reference rebuilt 4 times uncached; the cache, and the withdrawal); #465 (*read*: `ProblemSpec`'s optional source plus `is_eigenvalue` flag; ruled: the question is derived from the source, P1); #420 (`ProblemSpec` values outside its literals, P1); #418 (eleven spelling systems in the derivations vocabulary, P1); #305 (Sood cross-checks at published precision, P4's Sood family); #495 (`RegionMesh` "uniform" ULP defect: the fix must survive `RegionMesh`'s retirement, P1); #504 (the platform tag, P3).
- **Overlapping designs to reconcile, not duplicate:** #219 (*read*: the grand report's `GeometrySpec -> SpatialMesh -> MethodSpace` pipeline; this plan's specification and `geometry.mesh(d)` are its first two layers); #267 (`SNMesh -> MaterialMesh`); #393 (`AxisMesh` declares per-axis geometry: the intervals of P1); #411 (latent gaps in the layer-import contract: P0's whitelist and `layering.rst` fixes); #406 (a save/restore serialisation root: adjacent to the cache payload format); #358 (a red should invalidate its dependent cone: `Invalid` making consumers skip is one instance); #455, #387, #334 (phantom `catches`/`verifies` markers: the withdrawn-test counting of P0).
- **The Peierls Nyström improvement work already filed:** #100, #101, #103, #105, #109, #111, #112, #115, #116, #117, #123, #128, #129 (*read*: a 22.5% planar-limit discrepancy), #132, #140, #142, #143, #144, #255, #419, #499, #500 (22 issues). None states that the references are not research grade or that they are withdrawn. Proposal: one umbrella issue for the withdrawal and the return criterion, linking these 22 as its work, instead of new scattered issues.
- **In conflict with the user's reference ban:** #68 (*read*) and #21 propose a very-high-statistics Monte Carlo run cached as a reference. Monte Carlo is a consumer of references, never a source (`vv-principles`, "Ancillary references"), and a high-resolution numerical run is banned as a reference; proposal: close both, citing the ruling.
- **No existing issue** for the CP slab two-region verification gap or the package-`__init__` import breadth: those two are new.

**The solution objects** (the reply of this date, `[HYPOTHESIS]` until ruled; it corrects the user's reading that the specification carries the right answer):
- **The specification carries no answer.** It is the question only. There is no "right answer" object: the exact answer is known only through solutions, each with a certified bound.
- **A `ReferenceSolution` is one answer to one specification, from one producer**: a generator (closed form, MMS, F_N, trajectory resolvent, Case, Galerkin, Peierls), identified by its trace and settings, or a publication, identified by a `Citation`. The user's reading ("the answer given by the apparatus we constructed") is right for the first kind; the second kind is the same object with a different producer.
- **The production result** is the existing `Solution` five-tuple with its `ExitCertificate`; a gate compares one of its observables with the reference's, through the comparison verb, and nothing of the comparison is stored.
- **How a certificate becomes `Valid`**, by establishment method: `Exact` needs a symbolic self-check (the closed form satisfies the equation and the boundary laws with zero residual; for MMS, the stored source equals the operator applied to the manufactured solution); `ConvergedLadder` needs two legs, the ladder showing the theoretical rate (convergence, which gives the bound) and at least one structurally independent anchor for reduction correctness (a corroboration with an independent solution, or a closed-form limit or invariant such as particle balance, the row sum, k-infinity of the homogeneous case), because a ladder proves convergence to the solution of the equation the generator implements, not that it is the right equation (`vv-principles`, the semi-analytical ladder; ERR-032); `Published` needs its citation to resolve, and its bound is the printed precision.
- **A corroboration** compares two solutions of one specification per observable (`|a - b| <= bound_a + bound_b`) and records how independent they are; a pair sharing an identity, an integrand or an input object is a consistency check and never an anchor (X4). A failed corroboration makes both `Invalid`.

### Discussion 5, third exchange (2026-09-24): Monte Carlo's role, and the solution and certificate objects

**Monte Carlo, ruled and actioned.** The user: high-statistics Monte Carlo is not formal verification but sets a bar for model approximations (multigroup cross sections, anisotropy), backed by validating the Monte Carlo implementation against ICSBEP. `[LANDED]` issue #505 records the role; #21 and #68 closed "not planned" pointing to it, with their finding carried over (the Monte Carlo heterogeneous gate uses a CP eigenvalue with a different boundary condition as a proxy reference; its fix is a formal reference).

**The user's proposal, verbatim:** *"I think we need 2 objects: (1) PublishedSolution to carry the reference results for example published in Sood's report, and (2) ReferenceSolution to carry out solution, made by us. Likewise, when we're doing Validation, we will need an object to carry the Experimental results and our results. We might need 2 types of Certificate for the verification, a (1) ReferenceCertificate, which to latch into Valid state compares our ReferenceSolution against PublishedSolution and certifies for example, the Fn method implemented, and (2) the VerificationCertificate, which compares the ReferenceSolution to the production solver solution. See if this sharpens the picture or attack the architecture and refute the proposal giving your principles to do so."*

The evaluation `[HYPOTHESIS]` until ruled (session reply of this date):
- **`PublishedSolution` and `ReferenceSolution` as two types: accepted.** They differ in capability, not only in producer: a published solution is a finite set of printed values with printed precision, cannot be recomputed, has no trace and no field representation; ours is recomputable, trace-keyed, cached, and evaluates any functional of its field. One type with a producer field would make "evaluate the field of a printed table" spellable (Pattern 4). What they share is a protocol: "observables of one specification, each with a bound".
- **`ReferenceCertificate` (the latching one, three states) on our `ReferenceSolution`: accepted, sharpened.** Its anchors are not only published solutions: most references (MMS, closed forms, the trajectory resolvent on new geometries) have no published counterpart, so the certificate is established by `Exact` or by `ConvergedLadder` plus an independent anchor, and a comparison with a `PublishedSolution` is one kind of anchor, required whenever one exists for that specification. It certifies a generator on one specification, not the method in general. Two questions are separated: agreement with a published result from the same method (Sood's F_N values against our F_N) certifies the implementation, and a same-method comparison is the right instrument for that; certifying that the answer is right needs a structurally independent method or a closed-form limit (X4).
- **`VerificationCertificate` (production against a certified solution): accepted as the comparison's return value, refuted as a stored or latched object.** Production changes on every commit and carries no trace, so a stored verdict about it is stale by construction; the ruling that the acceptance criterion lives in the test (discussion 3, R-Q1.4; Oberkampf & Trucano: comparisons are not stored with the benchmark) stands. The comparison verb returns a `VerificationCertificate` holding the observable, the production value, the reference value and bound, the tolerance, and the floor check (the `OrbitCertificate` precedent: return the structure, not a bool), and the test asserts on it. Its reference side is any certified solution: a `ReferenceSolution` whose `ReferenceCertificate` is `Valid`, or a `PublishedSolution` (#305's direct Sood cross-checks).
- **The validation side mirrors it:** an `ExperimentalResult` (the benchmark-model value with its experimental uncertainty) and a `ValidationCertificate` returned by comparing a production or Monte Carlo result with it. `[REFUTED 2026-09-24]` "#505 places high-statistics Monte Carlo there": the user ruled that comparing anything against Monte Carlo is code-to-code benchmarking (L4), neither verification nor validation (*"At best it gives evidence of incorrectness [...] Useful, but not verification, neither validation."*); only the Monte Carlo code against experiment is validation (L3), of the Monte Carlo code itself. #505 was corrected the same day (title, body, label `level:L4`).

**Point 3 ruled** (the user, 2026-09-24: "Regarding point 3, 100% agree. That was the intended result."): the `VerificationCertificate` is the comparison's return value, asserted by the test, never stored or latched. The user judges the architecture "has reached quite a nice state" and asked for an adversarial attack: holes it cannot answer, or cases it is unprepared for that would force a major rewrite. The consolidated statement the attack reads is the section "The architecture as it stands", near the top of this file.

## The distinction the worklist needs `[R]`, to be measured per file before acting

A slow gate spends its time in one of three places, and each has a different remedy:

| where the time goes | example `[R]` from file names; read each first | remedy |
|---|---|---|
| **the reference's generator** (Peierls Nyström solves, trajectory-resolvent integrals, adaptive mpmath quadrature) | `test_continuous_registry_lazy.py`, `test_peierls_specular_bc.py`, `test_peierls_reference.py` | cache the generator's output, keyed on the generator (this plan's core) |
| **the SUT** (production solves at fine refinement) | `test_l1_standoff_slab_cylinder.py` (the cylinder ladder), `test_phase_c_crosscheck.py`, `tests/gates/mc/test_convergence.py` | never cached: size the refinement to the claim (for example the fewest rungs that fit an order), and move a fuller ladder to an explicit, scheduled tier |
| **a self-convergence study of the reference itself** (the reference checked against a finer copy of itself) | `test_peierls_rank2_bc.py` refinement monotonicity, `test_peierls_convergence.py` | this checks the generator, so it runs when the generator changes: the same key as the cache |

The third row is the insight that makes R1 coherent: **a test of the generator is itself keyed on the generator**. When the generator has not changed, neither its cached output nor the verdict of its own convergence tests can have changed.

## Candidate ontologies (first draft; each is a hypothesis for the discussion)

**C1. Memoise functions: a decorator on the generator function.** Like `sood_cache`, but keyed on the generator's identity instead of git HEAD.
- Makes unspellable: recomputation of an unchanged call.
- Leaves welded: *what the generator is* (the function body only? its transitive imports? mpmath's version?), decided implicitly by whatever the hasher reads. A test that computes its reference inline has no function to decorate.

**C2. The reference as a typed artefact: `ReferenceArtefact(generator identity, inputs, precision, payload)`.** The payload is the generator's output as data (arrays, scalars). The key is derived from the generator identity, the inputs and the precision, never from git HEAD.
- `ContinuousReferenceSolution` is then built from an artefact, so the callable is a view over data rather than a closure over live mpmath state.
- Makes unspellable: a reference that does not know what produced it, and a stale artefact served as current, since the key mismatches and it is a miss.
- Leaves welded: how the generator identity is computed (C4).

**C3. Generator-keyed test selection: tests of the generator declare their generator**, for example a marker `@pytest.mark.generator("peierls_nystrom")`. The runner skips them, with a stated reason, when that generator's identity matches the identity at the last recorded pass. This is how the third row of the table stops costing time.
- Makes unspellable: rerunning a generator's own convergence study when it cannot have changed.
- Leaves welded: where "the last recorded pass" lives, and whether a skip that depends on history is acceptable in a gate at all (a question for the user).

**C4. Generator identity (needed by C1, C2 and C3).** Options:
- **(a)** A hash of the source of the generator's module and its transitive first-party imports under `orpheus/derivations/`. Nexus's import graph could compute the closure.
- **(b)** A hand-declared `GENERATOR_VERSION` constant per generator, bumped by hand. Simple, but prose-as-enforcement (X3), and it goes stale silently.
- **(c)** (a), plus the versions of the trusted libraries below the line (mpmath, scipy, numpy).

`[R]` (a) with (c) is the only option whose staleness is structural, not remembered.

## Questions for the discussion

1. **The object.** Is the cached thing a `ReferenceArtefact` (C2), and does `ContinuousReferenceSolution` become a view over one? What about references computed inline in tests, not through the registry?
2. **The identity.** Which C4 option? Does the identity include the trusted-library versions?
3. **Where the cache lives.** Options:
   - local `.cache/` only (each machine regenerates once);
   - the Actions cache for CI;
   - committed artefacts (reviewable and reproducible across machines, but repository size, and ruling R20's spirit about what is tracked);
   - an artefact store keyed by identity.
4. **The freshness proof.** How does a gate prove it read a fresh reference? For example, the artefact carries its identity; the gate recomputes the current identity (cheap) and refuses a mismatch; a mutation of the generator source must turn the gate into a miss (X1).
5. **Generator self-tests (C3).** May a gate be skipped because its generator is unchanged? If not, where do the generator's own convergence studies run, and when?
6. **The SUT-side cost (the second row).** Size each ladder to its claim, and move the fuller ladder to a scheduled tier? That meets R1's fallback clause ("If we can't reach that point, we will think about how often it runs").
7. **Retirement.** Does the Sood cache (git-HEAD-keyed, no gate consumer) retire into the new machinery, or become its first consumer?
8. **Determinism.** A reference must be bit-reproducible from its key for the cache to be sound. Adaptive mpmath quadrature at fixed precision is deterministic `[R]`; anything seeded or thread-dependent is not. Measure before caching each generator.

## Rulings

(none yet)

## Refuted candidates

- **Versioning a cache entry by git HEAD** (the Sood cache's default): refuted for R1's question. Every commit invalidates every entry whether or not the generator changed, so the slow path runs after every commit. The fact it establishes: a whole-repository version is too coarse to be a generator identity.

## Related

#405 (the campaign, and its plan `.claude/plans/test_runtime_405.md`); #211 (reference selection, caching and parallelism in the SN suite); #212 (the lazy registry builders); #504 (platform-bound bit gates: a cached reference must not inherit that defect, see question 8); #404 (the pre-existing `phase_e` red).

## P0 landed (2026-09-25, `main` at 9187d998, CI `gates` run 36115864473 green)

| commit | what | the evidence in its body |
|---|---|---|
| a4143e73 | ERR-032: the phantom catcher replaced before the withdrawal (SymPy origin `origins/white_slab.py`, E2 re-tagged, E3 mpmath balance) | under the ERR-032 arm, 24 red, the 2 SymPy rows green |
| 2a5b19ee | CP white-boundary infinite-medium gates (`tests/gates/cp/test_white_boundary_infinite_medium.py`), the theory label `cp-white-cell-infinite-medium` | 4 mutation arms; the flatness row alone sees a one-row `P_cell` conservation defect |
| c5283051 | the mechanism: `Withdrawal`, the lock on 34 sites, the marker at 101 sites, the accounting, arm 7, gates M1 to M7 | (a) markers stripped: 288 red; (b) arm 7: 5 entries named; (c) 217 passed, 295 skipped, 1 xfailed |
| 9187d998 | the present-tense-false text | the layer gate, 370 passed after the whitelist removal |

Corrections that supersede earlier text in this plan and in the P0 specification:
- The lock surface is 34 sites, not 31: the specification's own list names 33 functions plus `BoundaryClosureOperator`, locked at its `__post_init__` `[M]` (qa's AST census over `orpheus/` and `tools/`).
- #509 starts at a group optical radius of about 4, not at "R >= 5 mean free paths": the graded 4-group sphere at R = 4 (optical radius 4.2) is red. Three sphere configurations, 9 rows, are strict xfails (comment on #509, 2026-09-25).
- Measurement (a) of the specification reds 288 of 294, not 294: 6 cases call `pytest.skip` before they reach the lock `[M]`.
- An xfail absorbed the lock's refusal, in the call phase and in a fixture's setup; a `pytest_runtest_makereport` wrapper (ELEGANCE-DEBT[guard] #506) converts it to a failure, gate M7.
- The audit snapshot `docs/_generated/vv_audit.json` is stamped with a digest of the audit's INPUTS only; a whole-tree digest refused its own snapshot inside one Sphinx build, because the build rewrites tracked generated files between writing and reading it (caught by the strict build before commit).

Open after P0: the `test-durations` re-run (dispatched 2026-09-25, run 36116403378) is #405's measurement of the P0 saving; `[R]` about 557 runner-minutes. Next is P1 (W3, surgical: the main agent writes, test-architect specifies the gates first).

## P1 opened (2026-09-25): what the tree already has, and the questions it raises

The blast set is measured at `99ac3e66` by two explorers; the detail and the probes are in `scratch/reference_architecture/p1_explore_geometry.md` (with `p1geo/`) and `p1_explore_spec.md` (with `p1_spec/`).

Corrections to the P1 bullet and to "The architecture as it stands":
- `[REFUTED 2026-09-25]` "`RegionMesh` retired (keeping #495's fix)": #495 is OPEN, and the `"uniform"` branch at `orpheus/geometry/mesh.py:574` still derives slab volumes from `linspace` edges `[M]` (3, 5 and 3 distinct volumes at n = 5, 7, 11 on a width-3 slab; `equal-volume` gives 1). P1 fixes #495; what must survive is ERR-020's `_subdivide_zone` (`factories.py:33`) and its gates.
- The geometry value exists: `StructuredGeometry` (`structured_geometry.py:215`, 1-D only) holds a kind tag (`"SLB"`/`"CYL"`/`"SPH"`), `regions` (`Region(mat_id, outer_thickness_cm)`) and `bcs` (each a `BC` tag or a `BoundaryTraceLaw`). The boundary types already live in the geometry layer (`orpheus/geometry/boundary/`). `Mesh1D.from_geometry` (`mesh.py:472`) is the `mesh(d)` verb; its body reads only regions, coordinate and boundaries, so the move adds no module edge `[M]`. `BC` and `StructuredGeometry` are unhashable (`BC.params` is a dict) `[M]`.
- The posing types live at L1, not L2: `EigenPosing` and `SourcePosing` at `numerics/posing.py:88` and `:118`, the pencil at `numerics/pencil.py:43`; production mints them at L3. So the specification, which lives in the geometry and data layers, may name them.
- `Materials` (`data/materials.py:54`) is already the `{id: Mixture}` declaration; it is `eq=False` "until content identity joins the typed-axis identity family". `Mixture._identity_key` (`mixture.py:180`) is stable across hash seeds; only `hash()` is salted `[M]`. `Axis._structural_bytes` (`transport/.../axis.py:274`) is the precedent for a byte-level digest.
- The intensional source already exists for the boundary: `InflowSourceSpec.evaluate(space)` (`geometry/boundary/_source.py:54`), projected onto production unknowns by `AngularBoundarySourceSink.from_specs`. A `numerics.Field` is already discretised, so it is not the base of `Source`.
- Producers of each question today `[M]` (constructions of `CriticalSolution` and friends in derivations): `Eigen(k)` many; `FixedSource` the per-region source of `greens_function.py:1106` and the 13 MMS case types (the 12 SN ones evaluate Q on the method's own quadrature); `CriticalParameter` 5 sites (critical-size generators, 6 function-level producers); `Eigen(c)` 2 (Galerkin); adjoint 1 (`kinf_and_adjoint_spectrum_homogeneous`); `Eigen(α)` 0 (and `ALPHA_MAP` has 0 production users). The critical and c references are compared only reference against reference today; production reaches them only as `Eigen(k)` at the answered size.
- Kind vocabulary `[M]`: 13 closed vocabularies; 302 kind-tag strings in `orpheus/` (287 in derivations), 624 in tests; production discriminates only on the enum (`CoordSystem`, 92 reads). `"SLB"/"CYL"/"SPH"` maps one-to-one onto `CoordSystem`, a second spelling of it.
- `RegionMesh`: 114 calls in 46 files (tests 113), plus 2 production calls through the alias `_RM` (`cp/solver.py:933`, `moc/solver.py:144`); `method=` is passed 7 times, all in `test_structured_geometry.py`.

Candidate carve order `[HYPOTHESIS]`, each step green on its own: (1) the kind vocabulary: `StructuredGeometry` carries `CoordSystem`, the string tag retires; (2) the discretisation types and `StructuredGeometry.mesh(d)`, which replaces `Mesh1D.from_geometry` and `RegionMesh` at all 116 sites and fixes #495; (3) `from_homogeneous(width, boundary)`; (4) content identity: a stable digest for `Mixture`, `Materials` and the geometry (a hashable boundary representation); (5) the question types; (6) the `Source` types; (7) the `Specification` that composes them. Derivations' 287 kind strings migrate with their families in P4, not here.

Questions for the user (asked 2026-09-25): which question variants P1 builds; where `Source` lives and whether it shares the boundary source's protocol; whether `Materials` gains content equality; whether `AxisMesh` (#393) moves into geometry now.

Rulings of 2026-09-25 (P1's first checkpoint):
1. **The question:** the full sum type now: `Eigen(k)`, `Eigen(c)`, `FixedSource(σ)` and `CriticalParameter`, each with forward/adjoint (the adjoint `FixedSource` carries its detector); no `Eigen(α)` (0 producers).
2. **The source is not settled.** The user, verbatim: *"The tricky part is that source might be method dependent (how is anisotropy expressed?). So we need to check if this is true and how it would play out for different methods (diffusion, Sn, MoC is not well developed yet but we will need to predict how it would express a flat and non-flat one, etc). Changing source expression to behave like boundary law declaration might require to be deferred based on lack of information and example."* An investigation is opened (below); `Source` and the `Specification` that holds it wait for it, and may be deferred.
3. **Identity:** `Materials` gains content identity, as its docstring promised; `Mixture`, `Materials`, the boundary laws (`BC` made hashable) and `StructuredGeometry` get a stable digest over their structural bytes, and equality follows it. The specification holds `Materials`.
4. **`AxisMesh` stays where it is (#393).** The user, verbatim: *"geometry is not the right layer for mesh either. In the future, mesh should probably go to a mesh module, where it will be besides readers for unstructured meshes from external meshers. Mesh is an overlay over geometry and is a complex enough overlay that in the long term will deserve its own module."*

### P1, second exchange (2026-09-25): the source by method, and the mesh module

**The source investigation** (explorer, `scratch/reference_architecture/p1_source_by_method.md`, probes in `p1_source/`). Finding `[M]`: the source is not method-dependent; three things are. (1) The representable subspace: SN holds any values at its ordinates and one value per cell (2^d moments for LD); diffusion only ℓ = 0 (no current-source term exists); CP and MoC isotropic and flat per region; the 0-D baseline has no source entry; MC has no source entry and its direction sampling is biased (ERR-018, #24). (2) The projection onto the method's unknowns: SN takes values at its ordinates (and synthesises the μ = −1 value, which goes negative for a non-polynomial angular source: −0.114, −0.141, −0.077 at N = 4, 8, 16 for a beam-like step); CP, MoC, diffusion need region averages of the ℓ = 0 moment; MC needs a sampler. (3) The density normalisation, owned by the angular measure (W = 2 on the 1-D rules, 4π on the sphere rules, 1/4π in MoC). Today all 12 SN MMS cases do the projection on the reference side, and every one of their sources is anisotropic (angular degree 1 or 2). CP and MoC have no absolute fixed-source solve (#510).
Verdict `[R]`: one method-independent stored form serves every method: q(r, Ω, g) as a function (`Symbolic`), with `RegionwiseConstant` (an angle-integrated rate per region and group, cm⁻³ s⁻¹) as its isotropic special case. Rejected: truncated angular moments as storage (the truncation order is a method's choice, a truncated expansion can go negative and cannot be sampled); point evaluation on the method's space as the bulk protocol (averaging methods need integrals; `_source.py:39-42` records the same failure for diffusion's trace); a silent reduction to ℓ = 0 (anisotropy is decidable from a `Symbolic`'s symbols, so a method that cannot represent it refuses).
Open: the projection verb, one `project(space)` where the space declares what it needs (values at points, cell averages, region moments, a sampler), or separate verbs; whether the boundary source moves onto it.

**The mesh.** The user, verbatim: *"Do you want to create the mesh module now, having just the 1D mesher inside of it? And the Axis product which brings 2D and 3D meshes as well, maybe. It's a re-home, but not into geometry. [...] Regarding the import order, I think it's better for a Mesh1D to be constructed from a geometry and a discretization scheme like you proposed. The method .from_geometry needs to be adversarially attacked as well. What is the general constructor? The generalized constructor is the bare bones constructor that takes exactly the arguments that are necessary for constructing the object. [...] geometry and discretization might be the raw constructor, unless the raw constructor takes something even more raw than these 2 arguments."*
Measured: `Mesh1D`'s dataclass constructor is `(edges, mat_ids, coord, precomputed_volumes, bc_left, bc_right)` (`geometry/mesh.py:308-313`); `BC` is defined in `geometry/mesh.py:87`; the axis types (`AxisCoord`, `Axis1D`, `AxisMesh`, `RadialAxisMesh`) are in `transport/mesh/axis.py`; import sites of `Mesh1D`/`Mesh2D`/`RegionMesh` from `orpheus.geometry`: 221 test files, 18 production files, 8 example files; of the axis module: 21 test files, 9 production files (the 3 `numerics` hits are docstrings).

### P1, third exchange (2026-09-25): rulings on the mesh module and the source

The user, verbatim: *"1. Yes. 2. Consider the possible advantages of Mesh1D.from_cells as a transition method and at some point in the plan migrate the 435 calls to use a better constructor (maybe the geometry + discretization constructor) vs keeping the method from the point of view of correctness and long-term code cleaness (not effort... I know moving 435 tests is a chore). Also, consider which cases would use from_discretization and if there is a better way. 3. Build it"*

Ruled:
1. **The mesh module is created in P1:** `orpheus/mesh/`, in the input layer (may import `geometry` and `numerics`); `Mesh1D`, `Mesh2D`, the discretisations, the subdivision helper and the axis types move in; `BC` moves to `geometry/boundary/`; `MaterialMesh` stays in transport.
2. **The source is built in P1** as a stored function: `Symbolic` (q(r, Ω, g), SymPy, stored through `srepr`) with `RegionwiseConstant` as its isotropic flat special case; each method's projection (its frame's analysis) lands with its family in P4, SN first; a method refuses a source it cannot represent.
3. **Open:** whether `from_cells` is permanent or transitional, and whether `from_discretization` is the right spelling. A census of the direct `Mesh1D(...)` constructions by what their edges express is running (`scratch/reference_architecture/p1_mesh1d_census.md`).

**The census** `[M]` (explorer, `scratch/reference_architecture/p1_mesh1d_census.md`; an AST pass over the 1033 tracked `.py` files of tests/, orpheus/, examples/, tools/, docs/, 22 positive controls): 450 direct `Mesh1D(...)` constructions (tests 435 in 186 files, production 3, derivations 10, examples 2). By what the edges express: uniform 345, single cell 22, uniform per region by hand 2, one cell per region 35, a named rule by hand 6 (these five classes, 410, are a geometry plus a rule); irregular literal 27 (18 of them say the irregularity is the point); derived from another object 8 (5 only change the boundary conditions of an existing mesh; 3 are production's `legacy_mesh_from_axes`, `transport/mesh/axis.py:665-683`, a mesh to axes to mesh round trip whose only caller is `SNProblem.from_axes`); constructor refusals 5; random 0. No site reads a mesh from outside, refines a mesh, or builds a sub-mesh. `precomputed_volumes` is passed at 0 of 450 sites. Rebuilding through a geometry: the uniform sites are bit-exact on 189 of 190 only with equal-WIDTH spacing (the default `equal-volume` differs on 58 curvilinear sites); the irregular sites are bit-exact on 25 of 26 as one region per cell, the one miss 1 ULP off because `Region.outer_thickness_cm` re-adds thicknesses. 5 uniform sites decide materials by cell position (`np.where(r_mid < 1)`), i.e. their honest geometry has a breakpoint there. `Mesh2D`: 152 direct sites (130 uniform or single cell).

Proposal `[HYPOTHESIS]` (for the ruling):
- **One constructor, no `from_cells` and no `from_discretization`:** `Mesh1D(geometry, discretization)`, where a discretisation is anything with `partition(geometry) -> Partition`; `CellsByCount` and `CellsByMaxWidth` are rules, and an explicit `Partition` (cell edges per interval) is itself a discretisation that returns itself after being checked against the geometry's breakpoints. The mesh stores the geometry and the partition, never the rule, so two rules giving the same cells give equal meshes. The irregular sites become their physical geometry plus an explicit partition (the irregularity is a property of the discretisation, not of the geometry); the boundary re-dress becomes `Mesh1D(replace(mesh.geometry, boundaries=...), mesh.partition)`; the axis adapter builds a geometry and an explicit partition. `from_cells` would be a shim with no permanent consumer.
- **The spacing is a required, typed choice** (equal width or equal volume), never a default: 58 curvilinear sites change under the current default.
- **The geometry stores its breakpoints (outer positions), not thicknesses**, with a thickness constructor for the registries that speak in thicknesses: re-adding thicknesses is the 1-ULP miss, and an interval is by definition the stretch between two breakpoints.

Rulings (the user, 2026-09-25):
1. **No shim; migrate in P1.** All 450 direct sites move to `Mesh1D(geometry, discretization)` in P1, each checked bit-identical by the census's rebuild probe or its change explained. `from_cells` never lands.
2. **The spacing: the user's attack, adopted.** Verbatim: *"another options is CellsByCount.uniform_width(n) and CellsByCount.uniform_volume(n). This is an adversarial attack against your proposal."* Resolution `[R]`: the attack is right at the call site (the name states the spacing, no default exists, no spacing type is imported at 400 sites) and it follows the user's constructor principle (the bare constructor takes what is necessary; class methods are ergonomic specialisations). The general constructor still holds the spacing as a typed field, `CellsByCount(counts, spacing)`, because the spacing is an axis orthogonal to the count rule: the count rule chooses how many cells each interval gets, the spacing where their interior breakpoints go, and `CellsByMaxWidth` shares that axis. `uniform_width(n)` and `uniform_volume(n)` are its specialisations (the same count in every interval). No default spacing anywhere. Per-interval mixed spacing is not built (no consumer).
3. **The geometry stores breakpoints**, with a thickness constructor for the registries.
4. **Every mesh has a geometry behind it.** The user, verbatim: *"If some construction has no geometry behind it, we should consider the reason for it. Depending on the reason, the correct approach is to build the geometry. Another view is 'testing Mesh1D without depending on geometry', which also requires a reason for 'why should geometry not be constructed first?'. If a lot of things build it without geometry, but the right approach is building the geometry first, there is probably quite an opportunity for cleanup."* The census gives each class its reason: the 410 hand-written rules had none (cleanup); the 27 irregular sites are a geometry plus an explicit partition; the 5 re-dress sites derive a new geometry from the old mesh's; the 3 axis-adapter sites are a legacy round trip; the 5 constructor refusals test laws that now live on the geometry, the partition and the refinement relation between them, each testable with its own object. No class has a reason to avoid the geometry `[R]`.

### P1, the carve order (2026-09-25; supersedes the candidate order in "P1 opened")

Each step is one commit, green under `python -O -m pytest` on what it touches, with its gates landing with their first red (§6c). The main agent writes (W3); the test-architect specifies the gates of every step first (`.claude/plans/reference_p1_spec.md`, to be written).

1. **The mesh module, a pure move.** `[LANDED 60d28b22]` 2026-09-26 (full suite before 12301 passed / 2 failed, after 12326 / 4 → 2 fixed in the commit, the other 2 the known crosscheck red, since fixed at `6bcea45c`, and the worktree-only pyright artefact; the 27 new tests accounted: 25 step-1 gates, 2 from ERR-089). `orpheus/mesh/` (input layer; a new linter row: it may import `geometry` and `numerics`); `Mesh1D`, `Mesh2D`, `RegionMesh` (until step 3), the subdivision helper, and the axis types from `transport/mesh/axis.py` move in; `BC` moves to `geometry/boundary/`. Bit-identical: no behaviour changes, every import site re-pointed, no re-export shim.
2. **The geometry value.** `[LANDED b91f1591, merged to main 2026-09-29 at c82df2da; full suite 12416 passed, 315 skipped, 63 xfailed, 0 failed in 2 h 14 min; CI gates run 36572103414 green on b91f1591]` `StructuredGeometry` carries `CoordSystem` (the `"SLB"/"CYL"/"SPH"` tag retires) and its breakpoints (regions carry outer positions; `from_thicknesses` for the registries). It also landed, by the rulings of "P1 step 2 opened": the derived boundary (hollow bodies admitted), the method refusals S3.9 and S3.12–S3.14 (moved from step 3), and the shared `homogeneous_body` refusal.
3. **The discretisation and the one constructor.** `[LANDED 3a `2c62f9f2`/`2890efbf`, 3b `3a6468e6`, 3c `4318adc5`; the design moved to the Mesher and a bare `Mesh1D`, see "P1 step 3 opened" and after]` `[REVISED 2026-09-26]` the argument is `partition`, not `discretization`: `Mesh1D(geometry, partition)`, by the user's perfect-match ruling (the unqualified word belongs to the method's discretisation, `posing_sequence.md`). `Partition` (cell edges and exact measures per interval, checked against the geometry's breakpoints), `CellsByCount(counts, spacing)` and `CellsByMaxWidth(widths, spacing)` with `uniform_width`/`uniform_volume` specialisations and `2 * d`, the spacing types; `Mesh1D(geometry, discretization)` storing geometry and partition; the boundary laws, coordinate system and material map read from the geometry; `precomputed_volumes`, `RegionMesh` and `from_geometry` retire; the 450 direct sites and the ~102 `from_geometry` sites migrate (the census's rebuild probe checks each bit-identical, or the change is explained); #495 is fixed at its root. `Mesh2D` keeps its constructor (no 2-D geometry value yet).
4. **`from_homogeneous(width, boundary)`.** `[LANDED 2c62f9f2, with step 3a]`
5. **Content identity.** A stable digest over structural bytes for `Mixture`, `Materials` (equality follows it; `eq=False` retires), the boundary laws (`BC` made hashable) and the geometry.
6. **The question.** `[REVISED 2026-09-26]` the question values are `Eigen(parameter)` (with the mode selector `Dominant`/`Nearest(τ)`/`Enclosed(region)`; the parameter `TermCoefficient(channels)` or `GeometryExtent(which)`; `CriticalParameter` dissolves into it) and `FixedSource(source, point)`, in the input layer (`posing_sequence.md`, "The ontology as it stands"). `Eigen(k)`, `Eigen(c)`, `FixedSource(σ)`, `CriticalParameter`, each with forward/adjoint (the adjoint `FixedSource` carries its detector), naming the L1 posings' spectral map rather than a second vocabulary.
7. **The source.** `Symbolic` (q(r, Ω, g), stored through `srepr`) and `RegionwiseConstant`, with the anisotropy predicate that lets a method refuse; no projection (P4).
8. **The specification** composing materials, geometry (or none), question and source, with its digest.

### P1, fourth exchange (2026-09-25): what Eigen(c) and CriticalParameter are

The test-architect's specification is written (`.claude/plans/reference_p1_spec.md`); its steps 6 and 7 are swapped (the phase-space functions before the question, since the adjoint's detector is one). Ruled: `CriticalParameter` varies the outer extent, in cm (the user, 2026-09-25). The user's question on `Eigen(c)`, verbatim: *"Out of curiosity, can this be reformulated as a resolvent? So an alteration of the EigenPencil at the begining of the solver. And if yes, is that the proper way to express this? Why? (in other words, try to attack this from a different perspective, and then attack that perspective as well). Check how far that logic applies to other cases if this attack lands."*

The attack `[R]`, and it lands. The pencil (`numerics/pencil.py`) is the affine family `at(σ) = A − σM`; its eigenvalues are the poles of the resolvent `R(σ) = at(σ)⁻¹` (the module docstring states it). An eigenvalue question is a choice of WHICH terms of the operator algebra `L + C − S − N_2n − F` carry the parameter. The k question puts fission in `M` (`A = L + C − S − N_2n`, `M = F`). The c question puts all collision emission in `M` (`A = L + C`, `M = S + N_2n + F`): `at(σ)` is singular where every emission, scaled by σ, balances the loss, and in one group `c_crit = σ · c_material`. So `Eigen(c)` is the same `EigenPosing` on a different split; the same strategy (power iteration, or a resolvent composition) solves it with no new solver code; and it generalises to several groups (σ scales every emission), so my "one group only" was wrong. Correction to the question's wording: the split is made at the head of the posing (the Problem's last step, the pencil), not at the beginning of the solver; the strategy never changes.
The counter-attack. (a) The specification cannot hold operators (it is mesh-free and below L2), so it stores the NAME of the split, the carrier of the eigenvalue, and each method's posing builds the pencil from it. (b) The carrier is physics (which multiplicity the question asks about); the spectral map (k = 1/μ, or μ itself) is presentation, derived from the carrier. (c) The logic holds only where the parameter enters AFFINELY; the pencil's own degree contract says so, and names the counter-example: a delayed-neutron α problem is rational in α and needs a companion linearisation. A critical SIZE changes the domain, so it is not affine: it is a nonlinear eigenvalue problem, a root of k(p) = 1 over a family of pencils. So the dividing line is affine versus non-affine dependence on the parameter, and my earlier recommendation (fold c into `CriticalParameter`) put c on the wrong side of it `[REFUTED 2026-09-25]` by this analysis.
The reach. Affine, so one pencil each: k (fission carries), c (all emission), prompt α (the 1/v term carries; 0 producers), a critical soluble-absorber concentration (the absorber term carries), diffusion's buckling B² (the leakage term carries); the fixed-source question is a point `at(σ)` of the family (the subcritical multiplying source is `at(1)`); the adjoint is the pencil's `.H`. Non-affine, so a root over pencils: size and shape, temperature, delayed α (or a linearised companion pencil), any composition change that is not a scaling.

### P1, fifth exchange (2026-09-25): the carrier, the pencil and the Strategy's transforms

The user, verbatim: *"Don't be limited by what the machinery says today. [...] Before fixing directions, check the plans operator_strategy_realization_campaign and orpheus-operator-machinery-report to gain context in ontology previously found regarding this. The Problem side of posing should end in an EigenPencil or a FixedSource problem. Resolvent should not be determined at the Problem, because it's totally legitimate for a single problem to go through different resolvents in different types of investigation, so the Problem side of posing must have the right architecture to allow resolvent posing, shifting operator sides (for G-S splitting for example), etc. We should find if the EigenPencil the right way to do this is to have an EigenPencil method or to feed the EigenPencil into an object that consumes it and changes it in the appropriate way. And how many objects we need to achieve our objective of realizing a Resolvent, Expressing a change in sides to have an explicit representation of some term, the system splitting according to a certain space axis, etc. Then we don't build the entire machinery, but we do the minimal change to create the seed to the campaign that will. Never be limited in your response by the current state of the code. Always project the code into the platonic ideal of the future of what the code should be, and opportunities to seed the trajectory to this ideal."*

Read: `operator_strategy_realization_campaign.md` §1 (the pipeline; the PAIRED-CONSTRUCTION ruling: the assignment of pieces is the primitive, `M` and `N` are derived, the only mutation is a `transfer`, a split is a `refine`), P3 to P6; `orpheus-operator-machinery-report-v2.md` §I.1 (splittings all the way down), §I.4 (the pencil and the resolvent: the resolvent is not a new type, `R(σ) = pencil.at(σ).inverse()`), §I.10 (partition per axis, block digraph from operator and partition, schedule on the splitting), Part IV (stages 3 declaration, 4 posing, 5 strategy recipe, 6 instantiation on `A − σΠ`); `consumers_step2_design.md` §3.1 and §3.5 (landed: `OperatorPencil` is the Problem's last step, `EigenPosing`/`SourcePosing` its questions; "which cells are occupied is the Problem's, which functional of R within a cell is the Strategy's"; `Splitting` with a labelled piece set as the primitive, `transfer` owed at P6).

The synthesis `[R]` (for the user's ruling):
1. **The pencil and the splitting are two labellings of one structure.** The Problem's factors are the pieces of one operator sum `T = L + C − S − N_2n − F − B` (with `T = 1/v` joining at α). A pencil labels each piece "left side" or "carries the parameter": k puts fission on the carrying side, c puts all collision emission (`S + N_2n + F`) there. A splitting labels each piece of `at(σ)` implicit or explicit. Both are the PAIRED-CONSTRUCTION primitive (the labelled piece set with derived sides), both admit the same two moves (`transfer` a piece to the other side; `refine` a piece into parts that sum to it), and each consumer keeps its own law: the splitting's transfer preserves the operator exactly (`M − N = A`); the pencil's transfer preserves only the physical member `at(1)` and changes the family, which is why it changes the QUESTION (k becomes c).
2. **So the owner of a transform follows from its law.** A transform that changes the question (moving a piece between the pencil's sides) is a posing: Problem side, a new pencil from a new labelling of the same factors. A transform that keeps the question and changes the computation is the Strategy's, and it CONSUMES the pencil, never a pencil method: a point or path in Λ (a shift σ: the resolvent posing is `at(σ)` plus the Strategy's lowering of its inverse; a contour; a Laplace variable), a partition of a space axis (restriction and extension pairs with `Σ J_i R_i = I`, which refine each piece into blocks `R_i P J_j`), the implicit/explicit labelling of those blocks (Gauss–Seidel along energy, space or angle, or a boundary `lower`/`upper` refine), the schedule, the lowering. A pencil method for any of these would put Strategy data on the Problem (the stage inversion of report v2 Part IV), and one Problem legitimately feeds many of them.
3. **How many objects.** Five kinds, two of them generic: (a) the factors (exist, per method); (b) ONE labelled piece set with `transfer` and `refine` (generic, numerics; today in embryo as `sn/splitting.py`'s `Splitting`); (c) the pencil, derived from a labelling (exists as `OperatorPencil`; gains a constructor from a labelling); (d) a per-axis partition (Space-owned; not built); (e) the Strategy values that compose them (a shift, a splitting over the refined pieces, a schedule, a lowering). No `Resolvent` type (report v2 §I.4); no `transfer` or `split` method on the pencil.
4. **The specification's stored question names the carrier in data-layer vocabulary**, since it sits below the operators: what the eigenvalue multiplies, a closed set, `FISSION` (k) or `EMISSION` (every secondary from a collision: c); a method's posing maps the carrier to its factors. The next affine members (the inverse speed for α, an absorber density, the buckling term) are new carriers, not new question types. `CriticalParameter` keeps the non-affine parameters (the outer extent).
5. **The seed, within P1's scope:** the carrier set in the data layer, shared by the specification now and by the posings later; `Eigen(carrier)` in the question; and the campaign that builds (b) and (c)'s constructor and (d) chartered as an issue carrying this synthesis. No production operator code changes in P1.

## ⏸ COMPACTION POINT — 2026-09-25, before implementation (supersedes the earlier 2026-09-25 point, removed)

State: the plan is RULED POLISHED (the user, 2026-09-25) and implementation is authorised phase by phase. Nothing is built. `main` is green at the commit that carries this point.

Read, in order: "The architecture as it stands" (the consolidated nine points); "Ruled polished (2026-09-25); P0 preparation"; "The phases, revised after the W5 review" (P0 to P5); then P0's verification specification, `.claude/plans/reference_p0_spec.md` (the test-architect, 2026-09-25: the 31 withdrawn symbols, 26 files marked, 294 cases withdrawn and 199 kept, the mechanism, the catalogue and V&V-matrix changes, the ERR-032 catchers, the CP white-boundary rows, the commit order). Its probes are in `scratch/reference_architecture/p0probe/` (`placement3.json` is the exact marker placement); the session's other evidence is in `scratch/reference_architecture/`.

Issues from this plan: #505, #506, #507 (premise corrected by comment: a file-wide `verifies` marker in `cp/test_verification.py` mints 34 `flat-source` carriers, a #387 case), #508, #509 (the CP sphere white-boundary defect, reproduced by the orchestrator: flux not flat for R >= 5 mfp, growing under refinement; multigroup k off 5 to 7%); #305 re-scoped.

P0 rulings (asked 2026-09-25; the user: "I agree with the 5 recommendations"):
1. **#509:** the CP sphere rows land with P0 as strict xfails citing #509, plus a check that pins the defect; the fix waits for CP's turn in the direction of development.
2. **The lock:** a guard inside the 31 withdrawn functions refuses unless `ORPHEUS_RUN_WITHDRAWN=506`; tagged elegance debt under #506; retired at P4, when the `Withdrawn` certificate state takes over.
3. **ERR-063** is dormant with ERR-027 to 030 (its defect lives only in withdrawn code).
4. **Nexus:** an interim "dormant" column in the error index in P0, and deOliveira-R/sphinxcontrib-nexus#96 (filed 2026-09-25) to teach the graph the `withdrawn` marker.
5. **Scope:** `BoundaryClosureOperator` and the volume kernel's own identity tests are withdrawn, per #506's list; the opt-in runs them for the improvement work.

Then P0 starts on a branch (`fix/peierls-nystrom-withdrawal`), in the specification's commit order: the ERR-032 catchers, the CP rows, the withdrawal in one commit, the docs; the V&V matrix and the error index regenerated in every commit that changes markers; Sphinx rebuilt before committing tests.

## ⏸ COMPACTION POINT — 2026-09-24 (superseded)

Resume here after regrounding. Nothing is built. First re-read this file, then `.claude/plans/test_runtime_405.md` ("Step 1, measured" is the worklist). Then bring the eight questions to the user, one discussion at a time: this is ontology search, and the object (question 1) and the identity (question 2) come first. Before arguing any scope, measure the premise (`process-discipline`):
- read the two or three largest files on the worklist, to confirm which row of the distinction table each belongs to;
- read `_richardson_cache.json`'s producer.


## The posing sequence's seeds for P1 (ruled 2026-09-29)

The posing sequence (`.claude/plans/posing_sequence.md`) opened from P1 step 7 and has settled the three-layer ontology; its section "The work, ordered by dependence" lists what P1 takes from it. The user ruled (2026-09-29, "Yes. Agree on both"):
1. **The question values live in `numerics`, spelled as the posing ontology spells them**: `Eigen(parameter, point, mode)` and `FixedSource(source, point)` (`Evolution` is not needed by any reference yet), the parameter an opaque key to a system coordinate (`CellCoefficient | GeometryExtent | NuclideDensity`). This supersedes the 2026-09-26 ruling "the question values in the input layer", and `CriticalParameter` is struck (it never existed in the tree: `[M]` 0 lines in `orpheus tests derivations tools`, 2026-09-29).
2. **P1 step 3** (`Mesh1D(geometry, partition)`) carries the gate "the mesh refines the region partition" (every cell in exactly one region), and a discretisation digest keys on the partition and the frames, never on the materials.
3. **`Materials` content identity** hashes `NuclideDensity` as a NUMBER density (nuclei per volume; mass is a derived property).
Nothing else of the posing machinery enters P1 (no system, pencil, mode law, direct-sum space or ordering). P1 resumes at step 2 unchanged.

## ⏸ COMPACTION POINT — 2026-09-26, P1 after step 1 (supersedes the 2026-09-25 point for P1's state)

Read in order: "P1, the carve order" (step 1 LANDED `60d28b22`; steps 3 and 6 revised: `Mesh1D(geometry, partition)`, the question values `Eigen(parameter)` and `FixedSource(source, point)`); the P1 verification specification `.claude/plans/reference_p1_spec.md` (gates per step; §3.4 the full-suite mesh capture: 738 constructions move under R2, all 607 `from_geometry` constructions bit-identical, one SN regression snapshot re-derived at step 3); then `.claude/plans/posing_sequence.md` "The ontology as it stands" (the three-layer ontology P1's step 7 and the later campaign follow).
Next: P1 step 2 (the geometry value: `StructuredGeometry` carries `CoordSystem` and breakpoints, the kind tag retires, `from_thicknesses` for the registries; hollow geometries admitted with the derived boundary count), on branch `refactor/reference-specification` (merged to `main` at `c418f9a9`; recreate it from `main` if deleted). W3: the main agent writes; the test-architect's gates are in the spec's §1.2.
Rulings that bind P1 (all 2026-09-26): the mesh package (landed); one constructor `Mesh1D(geometry, partition)` with no shim and every mesh built from a geometry (the 450 direct sites migrate at step 3; `from_cells` never lands); `CellsByCount.uniform_width(n)` / `.uniform_volume(n)` over a typed spacing, no default; breakpoints, not thicknesses; hollow geometries admitted (a centre carries no law; a hollow curvilinear geometry needs an inner law); `None` retires as a boundary declaration (each site migrates to the law its consumer applied, proven by a second capture); scope-boundary guards: SN inner surface (#511), CP slab needs left = right (#513), MoC solid cylinder only (#514), MC periodic only (#513); the equal-width edge rule is R2 (affine in the measure coordinate); −0.0 canonicalised and NaN refused in digests; `Materials` gains content identity; the question values in the input layer.
Measured baselines: the full suite at `6abd2980` (pre-move) 12301 passed, 2 failed, 315 skipped, 56 xfailed in 3 h 18 min serial; a full run needs a worktree with the HDF5 store linked and, for the pyright gate, a `.venv` link (`ln -s <main>/.venv <worktree>/.venv`), else the pyright ratchet reads every numpy import unresolved.
Durable lessons of the session: see `posing_sequence.md`'s compaction point; plus: a known red on `main` found mid-campaign is fixed before the campaign's next merge (the crosscheck W2 took precedence over step 1's merge).

## P1 step 2 opened (2026-09-29): the reconciliation and the user's rulings

Reconciled against the tree at `441e9595` (plan-authoring §7): step 1 is an ancestor of `main` (`git merge-base --is-ancestor`), and `StructuredGeometry` still has its original form (`orpheus/geometry/structured_geometry.py:215`: `geometry: str`, `regions: tuple[Region, ...]`, `bcs`). The blast set `[M]` (`scratch/reference_architecture/p1step2/census.py`, an AST pass over the 1000+ tracked `.py` files, literal-name calls only; the factories' `cls(...)` constructions are 2 more, found by reading):
- constructions of `StructuredGeometry`: `orpheus/derivations` 2 in code (`la13511.py:351, 364`, `[M]` `p1step2/sites.py`) and 5 in docstring examples (`moment_space.py:108`, `basis_space.py:482`, `spectrum.py:539`, `billiard.py:192, 244`; `git grep` finds 7, the AST 2), tests 111; `Region(...)`: `orpheus/geometry` 6, tests 147;
- the tag's readers: 8 derivation sites that branch on `self.geometry.geometry` (`moment_space.py:156, 258, 491`; `basis_space.py:527, 725`; `spectrum.py:581, 601, 752, 1009`; `billiard.py:969`), and test-side comparisons (`cross_method/adapters.py:509, 525`, `cross_method/test_eigenvalue.py:157`, `test_la13511_to_geometry.py`, `test_structured_geometry.py`);
- the other readers: `regions[0].mat_id` 3 (`moment_space.py:187`, `basis_space.py:572`, `spectrum.py:620`), `geom.bcs[-1]` 1 (`spectrum.py:189`), `domain_extent_cm` 1 (`billiard.py:1002`), `for _ in geom.regions` 1 (`sood_registry/builders.py:51`), `Mesh1D.from_geometry` (`orpheus/mesh/structured.py:355`), and `origin=` 1 site (`test_structured_geometry.py:486`).

Rulings (the user, 2026-09-29, the four recommended options):
1. **The laws' field is `boundaries`**, paired one-to-one with a derived `boundary_points` (the positions `(r_0, r_R)`, or `(r_R,)` for a solid cylinder or sphere); `n_endpoints` retires (it counted the interval's two ends, not the boundary points).
2. **`Region` retires.** The geometry is `(coord, breakpoints, mat_ids, boundaries)`, one material id per interval.
3. **The #511 guard (S3.9) moves into step 2**, because step 2 is the step that makes a hollow curvilinear inner law declarable on a geometry: SN refuses a non-reflective inner law on a hollow curvilinear mesh, a `SCOPE-BOUNDARY[guard]` naming #511.
4. **The test sites migrate by what they say**: a single region becomes literal breakpoints `(0.0, L)`; a multi-region stack becomes `from_thicknesses(...)` (bitwise today's fold, S2.3) unless the test states positions, which become literal breakpoints.
`origin=` retires in step 2 (the first breakpoint is the origin; one test site).

Later rulings of the same exchange (the user, 2026-09-29):
5. **One shared refusal for the one-material reference generators** (`orpheus/derivations/common/homogeneous_body.py`), read by `Spectrum`, `MomentSpace`, `BasisSpace` and `Billiard`. It refuses a hollow body and a geometry holding more than one MATERIAL. The elegance review (S4) corrected the first spelling, which tested the interval count: several intervals of one material are one body. The user asked for "more than one interval"; the correction is reported to the user.
6. **Spec S3.12–S3.14 move into step 2** (the elegance review's C1: the step that makes a hollow body declarable carries every refusal it needs). CP refuses a slab whose two laws differ and any inner law on a hollow body (#513); MoC refuses anything but a solid cylinder (#514); MC refuses a non-periodic left law (#513). Each refusal is a SCOPE-BOUNDARY guard. CP and MC first resolve the law they read, then refuse the law they would drop, so an existing xfail on #180 (MC's reflective law) still sees its own `ValueError`. An undeclared `None` is admitted until step 3.
7. **The #511 guard's admitted law is verified as if derived** (the elegance review's C2). In 1-D radial symmetry a void cavity reflects specularly. `test_a_reflective_inner_law_is_a_void_cavity` compares SN's hollow body, which has a reflective inner law, with a solid body whose core absorbs with cross section ε. `[M]` 2026-09-29, at the gate's fixture: gap/ε is 0.3263 (sphere) and 0.9983 (cylinder) at both ε = 1e-6 and 1e-4. The gap extrapolated to ε = 0 is 2.0e-11 and 1.4e-10. A black core differs by 0.360 and 0.546. A 128-cell near-void core made SN refuse its own convergence certificate (honest residual 1.0e-11 against tol 1e-12). That was read as the certificate working on an ill-conditioned solve and is not filed; it is reported to the user.

Step 2's evidence `[M]` (all `python -O -m pytest`, serial, on the host `.venv`):
- **Baseline.** The 57 test files touching the geometry value, plus the two gate files, before the edits: 1678 passed, 11 xfailed.
- **After the carve.** 1744 passed and 8 failed. All 8 were fixed. 7 were one mechanism: the import fixer treated `CoordSystem` bound in a function's local import as bound in the whole file, and pyright's `reportUndefinedVariable` is the instrument that sees it. The 8th was the ledger window of the new guard.
- **The 30 test files reaching CP, MoC or MC:** 872 passed, 1 failed, then fixed (the ordering in ruling 6).
- **Collection.** `pytest --collect-only tests/gates`: 12748 items, 0 errors.
- **Mutation batteries** (in-process rebinding plugins, `scratch/reference_architecture/p1step2/mut/`): S2.x, 10 arms plus the spec's positive control (the old count table), each reddening its target row. The #511 guard: no-op, 9 of 9 refusal legs red. `homogeneous_body`: 2 arms, 5 rows each. CP, MoC and MC guards: no-op, 4, 2 and 2 refusal legs red, admitted legs green.

Step 3's scope after step 2 (read before opening step 3): S3.9 and S3.12–S3.14 have LANDED in step 2 (`tests/gates/mesh/test_hollow_inner_law.py`, `tests/gates/mesh/test_dropped_laws_are_refused.py`). The CP slab site migration "left := right" is therefore partly done: a declared differing left law is refused today, and an undeclared one is admitted. Step 3 retires `None` and migrates each `None` site to the law its consumer applied (capture 2, §3.2 of the spec).
Open items that step 3 owns:
- qa finding 4: a centre law declared on a DIRECT solid `Mesh1D` still passes the #511 door. It retires when every mesh has a geometry behind it.
- The elegance review's N3: consumers pick laws by position (`len(boundaries) == 2`, `boundaries[-1]`). An `outer_law` / `law_at(point)` accessor is owed when `Mesh1D(geometry, partition)` becomes a consumer.
- N1: the test-side `_bcs_for(coord, law)` helpers (3 copies). The step-4 `from_homogeneous` absorbs the one-interval case.
Open items that P4 owns: the generators store their body once, instead of re-reading it on each `_mat_id` access (C3); the Sood registry's (`sood_registry/case.py`, `La13511Case.to_geometry`) and billiard's coord/kind tables are one bijection (N2).
Filed: #535, `coord is not CARTESIAN` is a missing `CoordSystem` property (18 compare sites, the elegance review's S3 note).

## ⏸ COMPACTION POINT — 2026-09-29, P1 after step 2 (superseded by the point after step 2b, at the end of this file)

Read in order:
1. "P1, the carve order". Steps 1 and 2 have LANDED. Step 3 is `Mesh1D(geometry, partition)`.
2. "P1 step 2 opened", its rulings 1–7, and "Step 3's scope after step 2". S3.9 and S3.12–S3.14 have already landed, and the open items step 3 owns are listed there.
3. `.claude/plans/reference_p1_spec.md` §1.3 (step 3's gates) and §3 (the migration protocol: two full-suite captures, a detached driver with a log).

Next: P1 step 3, on a branch `refactor/reference-specification` recreated from `main`. Mode: W3, with test-architect gates from the spec.

Measured baseline at `c82df2da`: the full suite `tests/gates` gives 12416 passed, 315 skipped, 63 xfailed, 0 failed in 2 h 14 min, serial `-O`.

Durable lessons of step 2:
- A `[skip ci]` plan commit at the head of a push skips CI for the code commit under it. The measure used here: push the code commit to a temporary `ci/<topic>` branch and watch that run. The alternative is to order the commits so that the code commit is at the head.
- A migration's import fixer must be scope-aware. A name bound by a function-local import is not bound in the module; pyright's `reportUndefinedVariable` is the instrument that sees the difference.
- A SCOPE-BOUNDARY docstring puts `machinery:`, `ruling:` and `revisit:` within three lines of the token, or the ledger gate reds.
- A new refusal lands AFTER the method resolves the law it reads, so that an existing xfail keeps its own error type.

## P1 step 2b (planned 2026-09-29): each reference generator serves the bodies its solvers solve

Opened by the user's question on `homogeneous_body`'s wording ("it depends on what should happen based on the specs"), and by the explorer's census of Sood's report (2026-09-29; the problem table is in #536). Rulings (the user, 2026-09-29):
- route a multi-material body per owner now;
- the 2003 edition of Sood is the reference;
- the 20 missing multi-media problems were checked against the open issues first, and none covered them, so they are filed as #536.

**Goal, in the domain's terms.** A reference generator solves exactly the configurations its solvers solve. A geometry it cannot solve is refused with a stated scope. A geometry it can solve is never refused. Today the shared `homogeneous_body` refusal blocks solvers that exist `[M]`:
- Billiard's `sphere_mr`, `hollow_sphere` and `annulus` arms (`billiard.py:443-676`; #190 and #421 record that they are unreachable);
- the multi-region cylinder solver `solve_greens_function_cylinder_mr` (`greens_function_cylinder.py:643`), which has no Billiard arm at all;
- `solve_fn_slab_reflected_critical` (`fn_method/slab/reflected.py:773`), the one-group symmetric equal-Σt reflected slab, which only a test adapter reaches.

**Candidate means** `[HYPOTHESIS]`, for the user's ruling.
- *(A) A typed reading of the body.* `homogeneous_body` becomes one classification of a geometry into the body shapes the reference literature solves:
  - a homogeneous solid body;
  - a concentric layered solid body (sphere or cylinder, one material per layer);
  - a hollow body of one material (hollow sphere, annulus);
  - a symmetric reflected slab (core, and the same reflector on both sides).

  It is a closed sum type. Each owner matches on it and refuses the shapes it does not serve, with a SCOPE-BOUNDARY naming the missing solver.
  - Billiard serves a homogeneous body, a layered sphere, a layered cylinder (a new arm onto the existing solver), a hollow sphere and an annulus.
  - MomentSpace serves a homogeneous body and a symmetric reflected slab of one group with equal Σt.
  - Spectrum and BasisSpace serve a homogeneous body only.

  One classifier, one match per owner, so "the one place the body is read" stays true.
- *(B) Per-owner readers*, with no shared type. Rejected `[R]`: the four owners would each re-derive "layered", "hollow" and "symmetric" from breakpoints and material ids, which is the twin the shared function was made to prevent.

**Gates.**
- The classification's laws: every geometry maps to exactly one shape, or is refused. Positive and negative rows per shape, including the discriminating ones: a two-material slab that is not symmetric is not a reflected slab, and a slab with unequal Σt is refused by MomentSpace.
- Each routed arm is bit-identical to its bare solver function (#190's own done-when).
- Problem 4 through MomentSpace against its published value. The value comes from the existing hand-built case, `cross_method/cases.py:445`, until #536 transcribes it.
- Each owner's refusal of an unserved shape is keyed on its scope boundary.
- Every refusal and every route carries a mutation arm, with the old shared refusal as the positive control (it must redden the route rows).

**The 2003 edition** is a separate commit in the same change:
- an explorer census of every registry value (47 cases) against the 2003 tables;
- each differing value re-checked against its own derivation where one exists: k∞ of an infinite medium is the matrix eigenvalue, via `kinf_and_spectrum_homogeneous`;
- the misplaced problem-72 flux ratio fixed;
- the stale docstrings named in the census corrected (`la13511.py:274-278, 2168-2176`, `__init__.py:40-41`).

**Out of scope.**
- Transcribing the 20 problems (#536).
- The F_N cylinder (#170).
- Problems with no solver: 3, 30, 58–61 and 63–66.

**Closes** #190 and #421 when landed.

**Sizing:** 1 session `[R]`.

### P1 step 2b, after review (2026-09-29): three commits

The reviews:
- The elegance review (`scratch/reference_architecture/p1step2b/elegance_report.md`) found two blockers:
  - F1: the generators drop the boundary laws;
  - F2: the mixture literals in `la13511.py` are duplicated.
- qa (`.../qa_report.md`) found the same law defect: Billiard's `alpha` silently wins over the declared law, which moves a hollow sphere's k by 5.4 %.
- qa also found a defect in the newly reachable `sphere_mr` fixed-source arm: on a 2-group problem it reports `n_groups = 1` and returns group 0 only.
- The homogeneous paths are bit-identical `[M]`: 104 of 104 generator results, 47 of 47 values of `c`, and 4 of 4 adapter values of τ.

The user ruled (2026-09-29):
- **Laws join the served pattern now.** Each owner declares the laws it solves and refuses others through the one door. `Billiard`'s `alpha` parameter retires: its albedos come from the geometry's laws. A slab whose two albedos differ selects the asymmetric arm.
- **Collapse, then re-cite.** Provenance becomes the one citation record, with an edition field, and the mixture literals are shared. Then every case is re-cited to the 2003 edition once, and equations are spelled `(A.N)`.

The commits:
1. **Routing and laws.**
   - One law reader, `specular_albedo(law)`: a law's specular albedo, or a refusal. It replaces `Spectrum._extract_R_refl` and is the one definition beside `BC.to_alpha`.
   - Law patterns per owner. MomentSpace and BasisSpace take vacuum only. Spectrum takes an outer specular albedo, and a slab's two faces must be equal. Billiard takes any specular albedos.
   - Elegance F4–F9: every refusal goes through the door; `is_hollow` has one definition; the body is stored once; MomentSpace refuses the ignored kwargs on a reflected slab; Billiard reads the shape in one match; the gates stop copying production.
   - qa's defect in the `sphere_mr` fixed-source arm is fixed, with an ERR entry.
   - Non-unit-Σt reflected fixtures.
   - Problem 4 is gated at the published digit: converged F_N at N = 11 and 15 gives 0.4301459, and 2003's 0.43015 must sit within half a unit in the last place.
   - Every infinite-medium case's stored k∞ and flux ratio are reproduced by `kinf_and_spectrum_homogeneous`, since only 1 of 10 ratios is read today.
2. **Registry.** Shared mixture helpers (F2), and the Provenance collapse with an `edition` field (F3).
3. **Re-citation to 2003**, from `.../p1step2b/citation_map.md`:
   - 91 citation tokens in 44 of 47 cases;
   - `sood_table` on 46 of 47 cases;
   - `paper_id`, and its gate;
   - the three stale docstrings;
   - the 10 mislabelled author names;
   - the equation-numbered derivation name `derive_kinf_1g_eq_20_simplifies_to_eq_19`, with a retirement audit;
   - about 23 comment citations.

### P1 step 2b, the three commits as landed on the branch (2026-09-29)

**Merged @ `071e582e`** (`--ff-only`, 2026-09-29): `1aab17a6` (routing and laws), `cc7a802e`, `49576c27`, `f02d7d9d`, `e8b83fc9` (the archivist's pass over the rest of the tree, 212 citations in 49 files) and `071e582e` (three false claims about the Sood references). Both branches are deleted.

- `cc7a802e`, commit 2. `Provenance` is the one citation record, and the flat fields `sood_table`, `primary_reference` and `notes` are retired. 54 of 54 cases' flat fields were checked equal to their provenance before 162 keywords were removed. One shared builder per material; 47 of 47 mixtures are bitwise equal to the previous commit's.
  - **Deviation from the ruling, reported to the user, no answer yet:** there is no `edition` field. `paper_id` names the publication and so the edition ("Sood-2003"), because an `edition` field beside it would be a second spelling of one fact.
- `49576c27`, an unplanned fix found while re-citing: **ERR-092**. The claim "Sood Eq 28 has a typo" was a mis-transcription in `fn_method/origins/k_inf_derivations.py`. Its gate V_fn2.1 asserted that the derivation DIFFERS from that transcription, and any wrong transcription passes such a check. V_fn2.1 now asserts that the derivation EQUALS the printed equation. The claim is removed from 7 sites, and `fn_method.rst` records it as refuted. Mutation: restoring the swapped transcription reddens V_fn2.1, 1 of 39.
  - Filed **#537**: the tree's other published-typo claims (WM-72 Eq 17, Carlvik Eq 4b, Adams-Martin A.1, BMC Eq 50) are open to the same failure until each is read against its page image.
- `f02d7d9d`, commit 3. The registry and the k∞ algebra-of-record cite the 2003 edition.
  - `paper_table` changes on 46 of 47 cases, and 10 author mislabels are corrected.
  - Six equation-numbered symbols are renamed (e.g. `derive_kinf_1g_eq_a3_simplifies_to_eq_a2`). `dead_references` reads 0 of 66.
  - The gate reds when a 1999 citation returns (1 of 131).
  - Census: 0 tokens shaped like a 1999 citation remain, against 93 at the previous commit.

The pre-merge gate at `e8b83fc9`: `tests/gates` with `-m "not slow"`, serial `-O`, per tree (`scratch/reference_architecture/p1step2b/gate/`). It gave 12216 passed and 0 failed in 37 min. The step-2 baseline below (12416 passed in 2 h 14 min) is a different population, since it did not deselect `slow`; the two counts are not comparable. The 244 skips in `derivations` are 231 withdrawn tests (#506) and 13 explicit skips, 12 of them the stub cases.

#536's body now says Billiard reaches the reflected spheres and cylinders, and that problem 4 is gated at the published digit. #421 was closed by hand: GitHub read `Closes #190, #421` as closing #190 only.

**Open after step 2b, both ruled 2026-09-29 and closed by P1 step 2c:**
- **Naming.** the module `sood_registry/la13511.py`, the classes `La13511Case` and `La13511Truth`, `LA13511_CASES`, and the test files `test_*la13511*` and `test_fn_sood_table10_*` carry the 1999 report's number. The registry now cites the 2003 paper.

- **The `edition` field.** The ruling asked for one. `paper_id="Sood-2003"` carries the edition instead, since a field beside it would state the same fact twice.

## ⏸ COMPACTION POINT — 2026-09-29, P1 after step 2b (superseded by the point after step 2c, at the end of this file)

Read in order:
1. "P1, the carve order". Steps 1, 2 and 2b have LANDED (`main` `071e582e`). Step 3 is `Mesh1D(geometry, partition)`.
2. "P1 step 2 opened", its rulings 1–7, and "Step 3's scope after step 2": the open items step 3 owns.
3. `.claude/plans/reference_p1_spec.md` §1.3 (step 3's gates) and §3 (the migration protocol).
4. The two questions above ("Open after step 2b"). They do not block step 3.

Step 3 runs as W3 on a branch `refactor/reference-specification`, recreated from `main`.

Durable lessons of step 2b:
- **A gate that asserts two expressions DIFFER certifies nothing about either.** Any wrong transcription passes it (ERR-092). A claim that a published equation is wrong is read against the page image first; its gate then asserts the equality that the corrected reading satisfies. The tree's other typo claims are #537.
- **A re-citation between editions is a map, not an offset.** Tables and appendix equations shifted uniformly between 1999 and 2003. References shifted in three blocks, and some values changed. About ten of the registry's own citations were wrong under either numbering. Every row was measured against both texts (`citation_map.md`); a computed offset would have produced about ten wrong citations.
- **Stale prose clusters around a rewrite.** One edition swap found eight false present-tense claims (a displayed equation, counts, attributions, "placeholder" modules). Budget a claims check beside every citation pass, not only the renumbering.
- **GitHub closes only the first issue of `Closes #A, #B`.** Write one `Closes` per issue.

## P1 step 2c (2026-09-29): `Citation` lands, the Sood registry drops the 1999 number

Ruled by the user (2026-09-29):
- Rename everything that carries the 1999 report's number, since a mixed provenance is not justified here.
- Keep `paper_id` carrying the edition, with no `edition` field. After this step the edition is the `Citation`'s bibkey.
- "Citation now, split in P4". The case schema is not given a new interim name: its successors are this plan's `Specification` + `PublishedSolution`.
- "The case number is related to provenance": the problem number is the locator of the problem's `Citation`.

Plan of record: `/Users/rodrigo/.claude/plans/zesty-dazzling-hanrahan.md`.

**Merged @ `216cd017`** (`--ff-only`, 2026-09-29, the branch deleted). The pre-merge gate (`tests/gates`, `-m "not slow"`, serial `-O`, per tree): 12614 passed, 0 failed, in 38 min. The commits:
- `61b20383`: `orpheus/data/citation.py`, `Citation(bibkey, locator | None)`, with its laws gate and 11 primary-source entries in `docs/refs.bib`.
- `2644bd47`: `La13511Case.problem: Citation` and `La13511Truth.sources: tuple[Citation, ...]`. The schema moves to `sood_registry/case.py`, and the Sood registry's `Provenance` retires. The Atalay sphere case, a twin of Sood problem 14, retires. Gate: `tests/gates/derivations/test_registry_citations_resolve.py` (test-architect, 326 rows).
- `228d6f20`: `la13511.py` → `sood2003.py`, `LA13511_CASES` → `SOOD2003_CASES`, the tests renamed, and the table-keyed names keyed on Sood's identifiers.
- `e0ca4ab2`: the docs (the archivist).
- `216cd017`: the review fixes: two citations, the bibliography rebuilt by pure insertion, `ATALAY1997_CASES`, the `builders.py` weld, the complete parse of `sources` tagged `ELEGANCE-DEBT` #405, and `_problem`/`_printed`.

**Inputs for P2**, from the two reviews:
- `sources` pools the values' citations. 8 of 47 Sood cases have values printed in different places, so `PublishedSolution` cites per printed value: one `Citation` per observable.
- The locator is free text, so equality depends on spelling. Consider a typed locator (`Problem(n)`, `Table(n)`, `Page(n)`, `Row(table, key)`) if a consumer ever needs to compare locators.

**Filed:** #538. The Atalay cases drop their reflection coefficient R and their P1 moment: 4 of 6 cases pose a vacuum, isotropic slab. It is blocked on a case carrying its boundary laws, which is `Specification` work.

**Lessons:**
- A retirement's census found a twin that no gate could see: an Atalay-named copy of a Sood problem, citing a table the paper does not have. Moving provenance into typed citations is what exposed it. A citation that must name a real place in a real work cannot be written for a problem the work does not pose.
- `refs.bib` is loosely sorted, and Zotero's marker comments belong to the entry below them. Insert entries as text beside the largest smaller key; never re-join the file, and never read a diff with blank lines filtered out.

## ⏸ COMPACTION POINT — 2026-09-29, P1 after step 2c (supersedes the point after step 2b)

Read in order:
1. "P1, the carve order". Steps 1, 2, 2b and 2c have LANDED (`main` `216cd017`). Step 3 is `Mesh1D(geometry, partition)`.
2. "P1 step 2 opened", its rulings 1–7, and "Step 3's scope after step 2": the open items step 3 owns.
3. "The posing sequence's seeds for P1": seed 2 is step 3's.
4. `.claude/plans/reference_p1_spec.md` §0 (the measured facts the laws rest on), §1.3 (step 3's gates), §2 (#495) and §3 (the migration protocol).
5. "P1 step 2c": the inputs recorded for P2 (per-value citations; a typed locator), and #538.

**Step 3, what it is.** `Mesh1D(geometry, partition)`.
- The argument is `partition`, by the user's perfect-match ruling. The spec's `Mesh1D(geometry, discretization)` is superseded: the unqualified word belongs to the method's discretisation.
- The types: `Partition` (cell edges and exact measures per interval, checked against the geometry's breakpoints), `CellsByCount(counts, spacing)` with `.uniform_width(n)`/`.uniform_volume(n)`, `CellsByMaxWidth(widths, spacing)`, the spacing rules, and `2 * d`.
- The retirements: `precomputed_volumes`, `RegionMesh`, `from_geometry`, `bc_left`/`bc_right`, and `None` as a boundary declaration.
- The migration: about 450 direct `Mesh1D` sites and about 102 `from_geometry` sites. Each is checked bit-identical by the spec's rebuild probe, or its change is explained.
- #495 is fixed at its root. `Mesh2D` keeps its constructor.

**Its gates:** spec §1.3.
- To build: S3.1-S3.8, S3.10 and S3.11.
- Already LANDED in step 2: S3.9, S3.12, S3.13 and S3.14.
- The posing sequence's seed 2 adds two items ("The posing sequence's seeds for P1", above): a gate that the mesh refines the region partition (every cell in exactly one region), and a discretisation digest keyed on the partition and the frames, never on the materials.
- #495's gate is spec §2.

**Before the first edit:** the migration protocol is spec §3. Capture 1 records every `Mesh1D` the suite builds; capture 2 records the resolved boundary law per face. Both are full-suite captures run through the detached per-tree driver (`scratch/reference_architecture/p1step2c/gate/driver.sh` is the latest form). The last `-m "not slow"` baseline is 12614 passed and 0 failed, at `216cd017`, in 38 min.

**The open items step 3 owns:** qa finding 4 (a centre law on a direct solid `Mesh1D` passes the #511 door); N3 (an `outer_law`/`law_at(point)` accessor); N1 (the `_bcs_for` test helpers, for step 4). #535 (`coord is not CARTESIAN`) is related.

Step 3 runs as W3, the main agent writing and test-architect specifying the gates first, on a branch `refactor/reference-specification` recreated from `main`. The Sood registry is now `sood_registry/sood2003.py` (47 cases, `SOOD2003_CASES`) and `atalay1997.py` (6 cases, `ATALAY1997_CASES`). Both are built on `sood_registry/case.py`, whose `La13511Case`/`La13511Truth` retire in P4. `orpheus.data.Citation` is the typed citation P2's `PublishedSolution` will hold.

Durable lessons of step 2c (beside step 2b's):
- **A rename can hide a relocation.** `La13511Case` looked like a Sood class; its census showed the Atalay cases were built on it too, so the right move was a schema module, not a new name. Census the type's constructors before choosing a name.
- **Ask what the planned design already says.** The user pointed at the campaign's own `Citation`, `Specification` and `PublishedSolution` when an interim name (`BenchmarkCase`) was about to be minted for an object the plan retires.

## P1 step 3 opened (2026-09-29): the rulings and the baseline

Reconciled against the tree at `08fe3e72`: `Mesh1D` is `(edges, mat_ids, coord, precomputed_volumes, bc_left, bc_right)` (`orpheus/mesh/structured.py:168`); `from_geometry` routes `"equal-volume"` through `_subdivide_zone` (`orpheus/mesh/factories.py:33`) and `"uniform"` through `np.linspace` with `compute_volumes_1d` (the #495 branch); the five law-resolution sites are `resolve_boundary_conditions` (bound in `transport/method.py`, `sn/problem.py`, `diffusion/augmented_mesh.py`), `sn/solver.py::_apply_default_bcs`, `CPMesh._resolve_bc`, `MOCMesh._resolve_bc` and `MCMesh.__init__`.

Rulings (the user, 2026-09-29, the four recommended options):
1. **Three sub-commits on the branch, each green.** 3a: `Partition`, the spacing rules, `CellsByCount`, `CellsByMaxWidth`, `2 * d`, their law gates (S3.1–S3.7, the #495 gate), additive; with the named-face geometry constructors and `from_homogeneous` (ruling 3). 3b: `Mesh1D(geometry, partition)`, every construction site migrated, `RegionMesh`, `from_geometry`, `precomputed_volumes` and `bc_left`/`bc_right` retire; each law carried exactly as declared, and a `None` declaration made explicit as the law its consumer resolved (capture 2 names it per test), because the geometry already refuses `None` (see the correction below). 3c: `None` retires as a declaration, SN's `boundary_condition=` parameters and the four consumer defaults retire, proven by capture 2.
2. **The test-site migration is done by general-purpose agents, one test tree each, in parallel, capture-gated.** The main agent writes production, the comparator and the brief, and reviews every diff before committing. A tree is done when its post-carve capture matches the pre-carve capture bitwise, or every differing row falls in one class of spec §3.1 within its bound.
3. **Named-face geometry constructors** over the bare constructor, one per coordinate system: `slab(breakpoints, mat_ids, left=, right=)`, `cylinder(...)`/`sphere(...)` with `outer=` and `inner=` (required exactly when the body is hollow). They retire the position-coded law tuples at the call site (N3's call-site half) and absorb N1's `_bcs_for` helpers; step 4's `from_homogeneous(width, boundary)` lands with them in 3a.
4. **The discretisation digest (posing seed 2) lands in step 5 with the shared encoder**, not in step 3. Step 3 lands the refinement gate (every cell in exactly one region) and `Partition` equality. `Mesh1D`'s hash also waits for step 5: the geometry is unhashable until `BC` is.

`[REFUTED 2026-09-29]` "3b carries each law as declared, `None` still admitted" (the wording put to the user in ruling 1): `StructuredGeometry` has refused `None` since step 2 (`_check_boundaries`, `structured_geometry.py`), so a 3b site cannot keep an undeclared law. The fact that survives: 3b migrates every `None` site of a `Mesh1D` to the law capture 2 recorded its consumer resolving (a mesh reaching two consumers that resolve differently becomes two geometries); 3c keeps the retirement of SN's `boundary_condition=` parameters, `_apply_default_bcs`, the consumer defaults, and the decision on `Mesh2D` and the axis tuples. Reported to the user.

Decisions taken without a ruling (reported to the user): a partition rule's verb is `rule.partition(geometry) -> Partition`, and an explicit `Partition` returns itself after the nesting check; the spacing rules are one body parameterised by the power `p` of `T(r) = r^p` they space in (`EqualWidth` p = 1, `EqualVolume` p = the coordinate's measure exponent), with `CoordSystem.interval_measure(a, b) = c (b^d − a^d)` evaluated in `_subdivide_zone`'s order (`[M]` 0 of 300 000 per coordinate differ from it; the numpy-square spelling of `compute_volumes_1d` differs in 243 of 200 000 cylinder cases, because Python's `b**2` is libm `pow`); `StructuredGeometry.uniform_boundary(coord, breakpoints, mat_ids, law)` absorbs the three `_bcs_for` helpers (the per-coordinate constructors cannot); consumers read the outer law through a geometry accessor (N3's consumer half) instead of `bc_right`; `Mesh2D` keeps its constructor, and whether `None` retires on `Mesh2D` and the axis tuples too is decided at 3c from capture 2's reading.

**The baseline** `[M]`: the pre-carve capture runs in a detached worktree at `08fe3e72` (`.claude/worktrees/p1s3-baseline`, `.venv` and the HDF5 store linked), `tests/gates -m "not slow"`, serial, `-O`, one session per tree, with `scratch/reference_architecture/p1step3/capture/p1s3_capture.py` (both captures in one plugin: per construction the exact edges, volumes, material ids, coordinate system and declared laws, keyed by test and ordinal; per face the law each of the five consumers resolved). Its positive control (`test_p1s3_capture_control.py`) reads the #495 slab (3 distinct volumes at n = 5 on `[0, 3]`), the vacuum law SN injects into an undeclared mesh, CP's white default and SN's reflective default; every tree's result line prints `orpheus.__file__` inside the worktree (L22).

**The baseline, measured** `[M]` (2026-09-29, 17:25 to 18:04, `scratch/reference_architecture/p1step3/capture/pre/`): 12 614 passed, 0 failed over the 15 trees and the root (the step-2c count exactly; `_mutation` collects nothing, rc 5, as in every gate), every tree importing the worktree's `orpheus`. 3345 `Mesh1D` constructions (2737 direct, 608 `from_geometry`). Undeclared laws: 605 tests in 75 files carry at least one; SN resolves 1393 `None` faces to reflective, diffusion 144 to reflective, CP 108 to white, MoC 77 to reflective, MC 3 to periodic, and SN's fixed-source entries inject vacuum into 20 all-`None` `Mesh1D`s (reflective into 1). One test resolves an undeclared law two ways (`tests/gates/sn/solve/test_sn_adjoint_entries.py::TestSolveSnAdjointFixedSource::test_entry_family_role_types`: reflective and an injected vacuum). The per-test table the migration reads is `scratch/reference_architecture/p1step3/capture/undeclared_laws.json`. The comparator (`capture/compare.py`, positive control: identity reads 145 of 145 `cp` meshes class (a); planted mutations land in (b), (c) and STOP) gates each migrated tree.

Step 3a `[LANDED 2c62f9f2]` on the branch: `orpheus/mesh/partition.py`, the measure on `CoordSystem`, the named-face constructors; 253 gate rows (test-architect); `tests/gates/mesh`, `tests/gates/geometry` and the layer gate 1451 passed. Review (qa, elegance-enforcer) dispatched.

### P1 step 3a, after review (2026-09-29): the Mesher, and the rulings it brought

The review of `2c62f9f2` (qa `scratch/reference_architecture/p1step3/review3a_qa.md`, elegance `review3a_elegance.md`) found one blocker in both: a `Partition` held measures without the coordinate system that gives them meaning (`[M]` a cylinder's equal-volume partition accepted by a slab with the same breakpoints; `Partition(([0,1],), ([42.0],))` accepted). The user's answer was an architecture, tested in two rounds:
- Round 1, `Mesher(geometry)`, `.partition(rule)`, `.mesh(cells)`: `[REFUTED 2026-09-29]` as a separate object in 1-D, because the mesher held nothing but the geometry and the partition is the mesh (the three stages collapse to one function of (geometry, rule)). The fact that survives: the geometry owns the measure, so no measure exists apart from its geometry.
- Round 2, the user: the mesher gives a preview; `partition(rule)` builds an internal mesh that the mesher can check for quality, refine and later adapt to a field (a flux gradient); `.mesh` returns it; this is the seam for a front end. It lands: the mesher now holds state (the current mesh) and mesh-to-mesh operations (refine, adapt) that are functions of neither the geometry nor the rule alone.

Rulings (the user, 2026-09-29):
1. **The Mesher seed lands in step 3**: `Mesher(geometry)`, `partition(rule)` returning the mesher (chainable), `.mesh`, `refine(2)`. Quality metrics, preview, adaptive refinement on a field and a `Mesher` protocol for external meshers are filed as an issue.
2. **`Mesh1D`'s constructor is bare-bones**: the geometry and the realised cells (edges and volumes). The Mesher unpacks the rule and calls it; the adaptive rule will need the raw constructor anyway. The constructor checks each stored volume against the geometry's measure of its cell within a derived bound (the construction law; `[M]` all 3329 captured meshes within 4.62 ulp, the wrong-coordinate control 6.1e15 ulp) and that every breakpoint is a cell edge.
3. **The Mesher is the only call site**: `Mesher(g).partition(CellsByCount.uniform_width(8)).mesh`, in tests and production; no `Mesh1D(geometry, rule)` shorthand.
4. **One rule, or one per interval**: `partition(rule)` applies it to every interval, `partition((r0, r1, r2))` one each; the rules are interval rules (`counts`/`widths` stop being unions), and mixed spacing becomes spellable. Supersedes the 2026-09-25 "per-interval mixed spacing is not built".
5. **The measure has one correctly-rounded definition**: `T_d(r) = r^d` evaluated on arrays, the measure `c·diff(T_d(edges))`, owned by the coordinate system and asked through the geometry; the interval measure is its one-cell case. Bit identity with `_subdivide_zone` is dropped (`[M]` elegance review: 0 of 414 captured equal-volume intervals move; Python's `b**2` is not correctly rounded in 22 of 20 000, numpy's square in 0).

Consequence for the sub-commits: 3a is reworked (interval rules, the one measure, the review's fixes; `Partition` retires, since the mesh stores the cells and no free-standing measured partition exists); 3b adds `Mesh1D(geometry, edges, volumes)`, the `Mesher`, and the migration.

### P1 step 3b (2026-09-29): the bare Mesh1D, the Mesher, the migration `[LANDED 3a6468e6]` (merged to `main` with 3a, `2c62f9f2` and `2890efbf`; CI `gates` run 36671753995 green; #495 closed)

Rulings (the user, 2026-09-29): **option A**, `Mesh1D(coord, edges, volumes, mat_ids, face_laws)` is the general constructor and holds no geometry (the mesh knows each cell's material and each boundary face's law; per-face laws are the seed for 2-D faces that carry different laws); a specialised `(geometry, edges, volumes)` constructor only if a good case exists: none does (its one consumer is the `Mesher`'s own lift; relabelling, adaptation and the axis adapter start from a mesh or axes), so the lift is a private step of `Mesher.partition`. Before that the user asked what `Mesh1D` needs the geometry for; the measured answer (coordinate system and breakpoints only; materials and laws only rode along, ~30 + ~29 production reads through the mesh) led to A.

Production (`orpheus/`): `Mesh1D` with construction laws (edges strictly increasing, finite, `r_0 >= 0` radial; each stored volume within `_VOLUME_ULPS` = 8 ulp of `c T(r_{j+1})` of `coord.measure` of its cell; one int material per cell; one law per `coord.boundary_points(r_0, r_R)`, never `None`), `boundary_faces`, `outer_law`, bitwise `__eq__`, no hash (step 5); `orpheus/mesh/mesher.py`; `CoordSystem.boundary_points`; CP, MoC, MC read `outer_law` (their `None` defaults are gone; CP refuses every hollow body, whose inner law now always exists); SN's `_apply_default_bcs` has no `Mesh1D` arm; the axis adapter's `_face_laws_of_axis` gives the SN adapter mesh reflective for an undeclared axis law and for a hollow radial axis's inner face (`ELEGANCE-DEBT[guard]` to 3c); the default CP/MoC pin cells, the Sood builder and the MMS cases go through the `Mesher` (the SN slab MMS cases declare vacuum, what the fixed-source entry injected; the curvilinear ones drop the centre law; the MoC case keeps reflective with its equal-area annuli as `CellEdges`). Retired: `RegionMesh`, `from_geometry`, `_subdivide_zone`, `precomputed_volumes`, `bc_left`/`bc_right`.

The migration `[M]` (per-tree `run_tree.sh`, `-m "not slow"`, `-O`; the comparator against the pre-carve capture), 7 general-purpose agents plus test-architect (`tests/gates/mesh`, `tests/gates/geometry`: `test_mesh1d.py` 46 rows, `test_mesher.py` 93 rows, a 13-arm battery) plus the orchestrator (shared helpers, the regression generator, fixtures, harness, benchmark, examples). Every tree passes pytest; mesh classes a / b / c: sn/sweep 341/40/45, sn/operators 700/34/40, sn/solve 131/5/19, sn/eigenvalue 60/0/7, sn/acceleration 19/1/45, sn/verification 70/1/69, sn/primitives 158/3/0, sn/mesh 207/15/0, sn/architecture 144/2/0, sn/regression 14/1/0, transport 331/186/2, cp 134/4/7, diffusion 14 c, numerics 100/6/5, homogeneous 25/0/0, derivations 4/0/4; no STOP row anywhere. Accepted non-class rows, each verified from the comparator's json:
- sn/verification: 10 pre-unmatched = the undeclared case mesh plus the vacuum copy SN's fixed-source entry injected (the case now declares vacuum, so the copy is never built); 8 law rows differ only in a `PrescribedInflow` repr's memory address; `test_ld_slab_mesh_routes_to_cumprod_scan` sees vacuum where it saw SN's reflective default on the MMS case mesh: the case's law made honest, and a departure from the "two meshes when a mesh resolved two laws" rule (this test only reads the sweep routing, which reads no law, so one vacuum mesh serves; qa, 2026-09-29).
- sn/solve: 2 unmatched = the solver's default-law copies; 1 law row = a cached 2-D box built first by an operators test in the full pre run (capture order, not migration).
- sn/operators: `test_default_is_reflective` re-posed as `test_declared_reflective` (the undeclared state no longer exists). `[REFUTED 2026-09-29, qa]` "3 duplicates = the old vacuum re-dress copies": the fact is that the 3 post duplicates are class-(b) moves (1.39e-15 relative) the comparator could not classify, because its duplicate rule only sees a copy beside an exact pair; the acceptance stands.
- mesh, geometry, diffusion: unmatched rows only in deleted or new tests (the undeclared-default tests: `test_mesh1d_bc_defaults_none`, `test_mesh1d_backward_compat`, the `[undeclared]` parameters, diffusion's two `test_undeclared_*_default_to_reflective`; `TestRegionMesh` re-posed into `test_partition`).
- CP Peierls tests (slow, absent from the capture) declare white, what CP resolved.

Re-derived by spec §3.3 `[M]`: `tests/gates/sn/_data/bc_extraction_baseline/vacuum_bulk_{SLB,CYL}_seed{0,1,2}.npy` (the class b/c meshes of `test_vacuum_bulk_bit_identical_1d`; the old copies in `scratch/reference_architecture/p1step3/snapshots_pre/`): max `|Δ|/max|b|` 6.9e-16 (CYL seed 0), at most 3.5 ulp of each array's largest entry; the 128-ulp entry reading is a near-zero entry. The SPH snapshots are bit-identical (their meshes are). `test_dd_regression[slab_fixed_source_dd_n20]`, the spec's anticipated class-b case, passes its frozen snapshot unchanged.

The review of 3b (qa `review3b_qa.md`, elegance `review3b_elegance.md`), acted on: the comparator gained a declared-law check (a law swap between two meshes passed it before; positive control: a planted vacuum-to-white swap in `transport` is a STOP, and the real captures read 0 STOP), preferring an exact-law partner when pairing; `Mesh1D` refuses a non-positive volume (qa F4); the volume band is derived per spacing power, `2p + 5` ulp (elegance: 7.62 ulp measured on a legal sphere partition against the chosen 8); `Mesher.partition` asserts every breakpoint is an edge; one `parse_boundary_laws` for geometry and mesh; the axis adapter uses `coord_system` and its tag is split (ELEGANCE-DEBT for the undeclared axis law, SCOPE-BOUNDARY #511 for the hollow inner face). The full post capture `[M]` (21:26, before these fixes): 12 996 passed, 0 failed (12 614 before; +394 new gate rows, the undeclared-default tests deleted), every non-class row in the accepted set.

Open for 3b's close: the full post-carve capture (running); qa and elegance review; the archivist (13 Sphinx pages, ~98 lines name the retired API; `TestMesh1DFromGeometry`'s rename with the 5 doc lines citing its node id); `dead_references`; two examples were already broken before this carve (the SN demo imports `orpheus.sn.quadrature`, the diffusion demo reads `result.flux.bulk`).


## ⏸ COMPACTION POINT — 2026-09-29, P1 after step 3b (supersedes the point after step 2c)

Steps 3a and 3b have LANDED (`main` `3a6468e6`). Read, in order:
1. "P1 step 3 opened" and the two sections after it: the user's rulings (the Mesher, the bare `Mesh1D(coord, edges, volumes, mat_ids, face_laws)` holding no geometry, interval rules with one or one-per-interval broadcast, the one measure on `CoordSystem`), the migration's results, and the review.
2. `docs/theory/foundations/structured_geometry.rst`, section "The mesh" (`structured-geometry-mesh`): the landed design and its refuted candidates.

`[M]` The final full run (`-m "not slow"`, post-carve capture, the settled tree): 13 020 passed, 0 failed; the slow tier of the 220 changed test files 81 passed, 7 xfailed (marked before), 0 failed. The capture tooling (`scratch/reference_architecture/p1step3/capture/`: the two plugins, `compare.py` with its declared-law check, `run_tree.sh`) is reusable for 3c.

**Next: P1 step 3c.** What it owns:
- SN's `boundary_condition=` parameters (`solve_sn_fixed_source`, `solve_sn_adjoint_fixed_source`, `solve_sn_multiplying_source`, `solve_sn`) and `_apply_default_bcs`, and spec S3.10's signature leg. A `Mesh1D` declares every law already, so these act only on a `Mesh2D` or an axis tuple.
- Whether `None` retires on `Mesh2D` and on `AxisMesh`/`RadialAxisMesh` too, which decides `resolve_boundary_conditions`'s reflective default and the axis adapter's ELEGANCE-DEBT arm (`orpheus/mesh/axis.py::_face_laws_of_axis`). The pre-carve capture counted the injections: `Mesh2D` 3, axis tuples 6. Ask the user; this is the scope question.

Filed during step 3: #539 (Mesher quality, preview, adaptation, protocol), #540 (two broken examples), #541 (a stale SN theory claim). Follow-ups noted by the reviews, not yet filed:
- `uniform_width(1)` and `uniform_volume(1)` are two spellings of one mesh (64 sites; a named one-cell rule);
- `IntervalRule.__rmul__` is a protocol member that `CellEdges` implements only to refuse (#539's scope);
- the libm dependence of the worst-case volume-band rows (re-measure on Linux, never widen).

## P1 step 3c opened (2026-09-29): no boundary law is left undeclared

**Ruled (the user, 2026-09-29): `None` retires as a boundary declaration everywhere.** It was already gone from `Mesh1D`. Now every `Mesh2D` face (`bc_xmin`, `bc_xmax`, `bc_ymin`, `bc_ymax`) and every axis endpoint (`AxisMesh.bc_low`/`bc_high`, `RadialAxisMesh.bc_outer`) must declare a law. Retiring with it:
- SN's `boundary_condition=` parameter on `solve_sn_fixed_source`, `solve_sn_adjoint_fixed_source` and `solve_sn_multiplying_source` (`[M]` the explorer's census: `solve_sn` has none, so the earlier list naming it was wrong);
- `_apply_default_bcs` and `_as_problem`'s `boundary_condition` argument;
- the reflective default in `resolve_boundary_conditions` (`orpheus/transport/method.py`);
- the axis adapter's ELEGANCE-DEBT arm in `orpheus/mesh/axis.py::_face_laws_of_axis`.

Spec S3.10's signature leg is owned here. The migration declares, at each undeclared site, the law its consumer resolved in the pre-carve capture (`scratch/reference_architecture/p1step3/capture/undeclared_laws.json`), and is gated by the same capture-and-compare tools as step 3b.

`[M]` The premise, measured 2026-09-29:
- `boundary_condition=` has 64 sites in 25 test files and 0 in `orpheus/`.
- The pre-carve capture's injections into an all-`None` declaration: `Mesh2D` 2, axis tuple 5, `Mesh1D` 10 (the last are already gone at 3b).
- The reflective fallback's resolutions: `SNProblem` y faces 261 each, z faces 24 each; `DiffusionMesh` `xmin` 70 and `xmax` 32.
- `Mesh2D(` appears at 147 test sites and 5 production sites.

**Ruled (the user, 2026-09-29): `Mesh2D` declares per-face laws, like `Mesh1D`.** The four `bc_*` fields retire. The alternative, four required fields, was refused because a solid (r, z) mesh would then have to declare a law on its axis, or keep `None` meaning "no face here", a second meaning of `None`.

### The design (step 3c)

- **One element parser.** `parse_boundary_law(law, where)` in `orpheus/geometry/structured_geometry.py` accepts a `BC` tag or a `BoundaryTraceLaw` and refuses `None` ("None is not a boundary law") and any other object. `parse_boundary_laws` (1-D, positional) calls it per element; `Mesh2D` and both axis classes call it. `_check_boundary_declaration` (it admitted `None`) retires, and its docstring's rationale for the law arm (a law whose content is a function has no tag spelling) moves to `parse_boundary_law`.
- **The 2-D face inventory.** A `Mesh2D`'s boundary faces follow the same topology law as `Mesh1D`: along axis 0 the faces are `coord.boundary_points(edges_x[0], edges_x[-1])` (two on (x, y) or a hollow (r, z); only the outer surface on a solid (r, z), whose axis r = 0 is interior); along axis 1 (y or z) always two. Each face is named by the existing crosswalk `FaceLabel(axis_index, endpoint).face_name` (`orpheus/mesh/axis.py`): `xmin`, `xmax`, `ymin`, `ymax`, where a solid radial axis's outer surface renders `xmax`. SN's `bc` table already uses these names, so there is one naming.
- **`Mesh2D(edges_x, edges_y, mat_map, *, face_laws, coord=CARTESIAN)`.** `face_laws` is keyword-only and required: a mapping from face name to law whose key set equals the inventory exactly (a missing face and an extra face are both refused, naming the inventory). It is stored read-only, in inventory order. `Mesh2D.boundary_faces` returns the inventory.
- **The axes.** `AxisMesh(edges, bc_low, bc_high, ...)` and `RadialAxisMesh(edges, coord, bc_outer, ...)`: the laws are required and parsed in `__post_init__`; `bc` is typed without `None`. `with_uniform_bc` (protocol member and both implementations) retires with its one consumer, `_apply_default_bcs`.
- **The adapters.** `axes_from_legacy_mesh` reads `mesh.face_laws`; `legacy_mesh_from_axes` builds the `face_laws` mapping from the axes through `face_labels`; `_face_laws_of_axis` loses its `declared()` arm (the ELEGANCE-DEBT retires) and keeps the hollow radial inner face (SCOPE-BOUNDARY #511).
- **The resolver.** `resolve_boundary_conditions` reads the declared law with no default. `MaterialMesh._law_key` loses its `None` arm.
- **SN.** The three entries lose `boundary_condition=`; `_as_problem` loses its argument; `_apply_default_bcs` is deleted.
- **Production meshes that declared nothing.** The two 2-D manufactured-solution meshes in `orpheus/derivations/continuous/mms/sn.py` (`:723`, `:875`) declare the vacuum SN injected. `pwr_pin_2d` takes a required keyword `law` for its four faces: a factory does not choose a physical boundary silently.
- **Migration.** As at 3b: one agent per test tree, each site declaring the law its consumer resolved, gated by the capture comparator with the post-3b capture (`full_post2`) as the baseline. The five tests whose subject is the undeclared default are deleted or re-posed (the census, `scratch/reference_architecture/p1step3c/blast_set.md`).

### Review of step 3c, and two rulings (2026-09-30)

The reviews: `scratch/reference_architecture/p1step3c/review_qa.md` (no correctness defect; every checked declared law is the law its consumer resolved; the one accepted comparator row is the fixture cache, reproduced) and `review_elegance.md`. The test-architect's battery: 17 of 17 arms redden (`battery.md`).

**Ruled (the user): one `FaceLaws` value for both meshes.** The elegance review found face laws spelled two ways: `Mesh1D` a positional tuple (inner face first), `Mesh2D` a read-only mapping keyed by face name, and `boundary_faces` meaning positions on one and names on the other (`[M]` `zip(m2.boundary_faces, m2.face_laws)` pairs names with names in silence). The mapping was also unpicklable, and the reference cache (step 5) pickles its entries. The alternative refused: a local rename with an issue for the unification (the same concept stays in two spellings, the stop signal of Cardinal Rule 2).
- `orpheus/mesh/face_laws.py`: `face_inventory(coord, edges_of_first_axis, dimension)`, the one topology rule (`CoordSystem.boundary_points` along the first axis, both ends along every further axis, each face named by `FaceLabel.face_name`); and `FaceLaws`, a frozen, ordered, picklable mapping from face name to law, built by `FaceLaws.over(inventory, laws, where)`, which refuses a missing or extra face (naming the inventory, with the reason for a solid radial centre or a hollow inner surface) and parses each law with `parse_boundary_law`.
- `Mesh1D.face_laws` and `Mesh2D.face_laws` are both a `FaceLaws`, given as any mapping from face name to law: a slab `xmin`/`xmax`, a solid cylinder or sphere `xmax` only, a hollow one `xmin`/`xmax`; a 2-D mesh adds `ymin`/`ymax`. The positional tuple is retired. `Mesh1D.boundary_faces` (the positions) is renamed `boundary_points`, the vocabulary of `CoordSystem`; `Mesh2D.boundary_faces` retires (`tuple(mesh.face_laws)` are the names). `outer_law` reads `face_laws["xmax"]`; the three "is there an inner face?" length tests (`axis.py`, `cp/solver.py`, `mc/solver.py`) become `"xmin" in mesh.face_laws`.
- Measured before: `Mesh1D.face_laws` 31 production reads in 5 files, 23 test constructions in 7 files and about 40 test reads.
- Not in this ruling: `StructuredGeometry.boundaries` is a third spelling (a tuple, one law per boundary point; `[M]` 9 production reads, 144 test constructions in 58 files). The Mesher lifts it through `face_inventory`. Raised with the user at the step's close.

**Ruled (the user): a named cell's law is the model's.** `wigner_seitz_pin_cell` (white outer surface: the cylindricalised lattice cell) and `pwr_slab_half_cell` (reflective: both faces are symmetry planes) lose their `boundaries=` keyword; the law is documented as part of the model's definition. Another law is a different body, built with `StructuredGeometry.cylinder` or `.slab`. `[M]` 2 test sites override it today (`test_dropped_laws_are_refused.py:82`, `test_structured_geometry.py:379`). The alternatives refused: requiring the law (the name already fixes it), and keeping the overridable default (a default that a caller can silently inherit).

### P1 step 3c: merged @ `4318adc5` (2026-09-30; CI `gates` run 36748263634 green)

`[M]` The evidence:
- The full `-m "not slow"` suite under the step-3c capture: 13 096 passed, 0 failed (3b: 13 020). Against the post-3b capture, every `Mesh1D` is bitwise identical. Every test resolves the same (consumer, face, law) as before, except the deleted or re-posed tests of the undeclared default and the new gates. The one per-session-cache row seen in a partial run is absent in the full run, as qa predicted.
- The slow tier of the 124 changed test files: 27 passed, 0 failed.
- The test-architect's battery 2: 23 of 23 arms redden their targets.
- `sphinx-build -E -W` clean, `dead_references` 0.

Fixed on the way and catalogued: ERR-093. `solve_moc`'s default mesh carried the Wigner–Seitz cell's white law, which MoC refuses; it now uses `moc.solver.default_pin_cell_mesh`.

The comparator (`scratch/reference_architecture/p1step3/capture/compare.py`) now strips memory addresses from its law table. A planted law swap still mismatches. Two laws of one callable-bearing class swapped between faces do not mismatch; that is a limit of the tool, noted by qa.

## ⏸ COMPACTION POINT — 2026-09-30, P1 after step 3c (supersedes the point after step 3b)

Steps 3a, 3b and 3c have LANDED (`main` `4318adc5`). Read, in order:
1. "P1 step 3c opened", "The design (step 3c)", "Review of step 3c, and two rulings" and the merge record above.
2. `docs/theory/foundations/structured_geometry.rst`, sections `structured-geometry-mesh`, `structured-geometry-face-laws` and `structured-geometry-no-default-law`.

Open at the close of step 3:
- **The question put to the user:** `StructuredGeometry.boundaries` is a third spelling of face laws (a positional tuple, one per boundary point). `[M]` 9 production reads and 144 test constructions in 58 files. Should the geometry store a `FaceLaws` too, so that the Mesher's lift is the identity?
- **Follow-ups filed:** #542 (`uniform_width(1)` and `uniform_volume(1)` are one mesh), #543 (`IntervalRule.__rmul__`, in #539's scope), #544 (the volume band's libm dependence), #545 (a uniform-law `FaceLaws` constructor).

**Next:** P1 step 5, content identity ("P1, the carve order"; step 4, `from_homogeneous`, landed with step 3a). It makes `BC` hashable and gives `FaceLaws`, `Mesh1D` and the geometry their digests. The user's answer on the geometry question comes first, because it decides whether the geometry's digest is taken over a `FaceLaws`.

### 2026-10-01: the step-5 question answered, and a detour

The question put at step 3c's close ("should `StructuredGeometry` store a `FaceLaws`, so the Mesher's lift is the identity and step 5 digests one type?") opened the boundary-law ontology discussion (`.claude/plans/boundary_law_ontology.md`). **Answered (the user, 2026-10-01): no.** The geometry keeps its positional `boundaries` tuple for step 5; the ontology's architecture E will later replace it with deck transformations and integer boundary tags (the boundary conditions bound per method by `SNDiscretization(material_mesh, boundary_conditions)`), and that change re-keys the geometry's digest (a cache miss, never a stale hit). E executes after this campaign. The detour landed `8f9300b2` (ERR-094, ERR-095, CP's slab guard: 85 CP test declarations now reflective|white) and, next, a small reflective cleanup (`ReflectiveBoundary` without albedo), before step 5.

**Next after the cleanup: P1 step 5, content identity** ("P1, the carve order").

### 2026-10-01: the reflective cleanup landed, and one requirement for step 5

The reflective cleanup made `ReflectiveBoundary(axis)` parameter-free and made `LawSum`/`LawScaled` refuse a deck law as an operand, with the refusal in `__post_init__`. Merge record: branch `refactor/reflective-is-a-mirror`, commits `65f1dc7a`, `f7a9b79a` and `99b9847d`.

Its qa review measured a gap that belongs to this campaign. Loading a pickle never runs `__post_init__`, so a pickle written before the carve loads silently:
- `ReflectiveBoundary('x', 0.7)` loads as a mirror that ignores its 0.7 (ERR-094's mechanism, with no message);
- `0.7*R + 0.3*W` loads as an admitted `LawSum`.

`[M]` 2026-10-01 by qa, `scratch/boundary_ontology/reflective_qa_review.md` G1. No pickle is tracked today; the cache would be the first persisted store.

**Ruled (the user, 2026-10-01): a requirement on step 5, not a runtime guard.** The content key covers the schema of every persisted class (its fields, and a version), so an entry written before a carve misses instead of loading. The alternative, a `__setstate__` refusal on each class, was declined: it adds a guard per carve, where a schema-covering key fixes every future carve once.

The witness owed with step 5: a law or a composed tree pickled under an older schema must miss the cache, never load. Its fixture is a schema-changing carve such as this one.

### 2026-10-02: P1 step 5, content identity — merged @ `a5113ac0` (`ee9e8943` code, `a5113ac0` docs)

Landed on branch `refactor/content-identity`, uncommitted at the time of writing, so no hash yet. Gates: spec §1.5 of `.claude/plans/reference_p1_spec.md`; the design plan `/Users/rodrigo/.claude/plans/zesty-dazzling-hanrahan.md`; the census `scratch/reference_architecture/p1step5/census.md`. The record is `docs/theory/foundations/structured_geometry.rst`, section `structured-geometry-content-identity`; the axis's side is `docs/theory/foundations/spaces.rst`, subsection `spaces-generator-identity-third-answer`.

**What changed.**
- One encoder, `orpheus/numerics/content.py` (a leaf module importing nothing from `orpheus`): `encode(value) -> bytes`, recursive, every chunk type-tagged and length-prefixed; `content_digest(value)`, the blake2b-256 of it; `ContentlessError(TypeError)`, naming the path to the part that has no content; and the `ContentIdentity` mixin, whose `==` and `hash` are derived from the digest (the hash is the digest's first 8 bytes, the same in every process). A value holding a contentless part is equal only to itself and unhashable.
- An object encodes as its schema tag (`module.qualname|v<__content_version__>|<part names>`) followed by its parts. This is the user's ruling of 2026-10-01 above ("a requirement on step 5, not a runtime guard"): an entry written under an older schema misses. Its witness is gate S5.8 (a field added, renamed or reordered, a version raised, a class moved or renamed: each moves the digest). The end-to-end witness, a pickled law under an older schema missing a real cache, waits for the cache.
- Moved onto the encoder: `Mixture` (its `_identity_key` retired; `MaterialMesh._contractibility_key` folds `content_digest`), `Materials` (content identity, ids coerced to `int`, declared order kept and not content, picklable through `__reduce__`), `BC` (params a read-only mapping of real numbers; hashable; pickles), every registered boundary law and `LawSum`/`LawScaled`, `StructuredGeometry`, `FaceLaws` (no longer equal to a plain `dict`), `CellEdges`, `Mesh1D` (hashable; its hand-written `__eq__` retired; derived fields `compare=False`), `Mesh2D` (read-only copied arrays, `-0.0` canonicalised, `==` works), the `Axis` family (`content_parts` replaces `_identity_key`; `_structural_bytes` retired), and the five space-name digests (`FunctionSpace.of_axes`, the angular and scalar trace spaces, the radial characteristic spaces, `FullFieldSpace.from_blocks`). `MaterialMesh`'s law key keys a `BC` tag by the tag itself (its `("BC", kind, sorted params)` arm retired).
- NaN is refused at construction, by parsing at the boundary: `parse_real`, `BC` params (a `str` or `bool` value is a `TypeError`), `AlbedoBoundary` and `WhiteBoundary` albedo, `ConstantInflowSource.value`, and every array of a `Mixture`.
- Gates: five files, `[M]` 2026-10-02 253 collected and 253 passed under `python -O -m pytest`; the RECORD fingerprint (S5.7) pins three digests; a new step of the `gates` CI workflow runs the five files on Linux.

**The rulings of 2026-10-02 (the user).**
1. A real scalar's content follows `==`: `True == 1 == 1.0` is one value, `-0.0` is `+0.0`; NaN and an integer beyond 2**53 are refused.
2. Equality is the digest: every content type's `==` and `hash` derive from it (X4, one definition), so equality and the cache key cannot drift. The two string arms (`VacuumInflow() == "vacuum"`, `ReflectiveBoundary(axis) == "reflective"`) are kept, by the ruling of 2026-10-01; `[M]` their hashes differ from the strings', so a law and its kind string must not share a set or a dict.
3. The old encoders moved: `Axis._structural_bytes` and the four space-name mints call the one encoder.
4. `Mesh2D` is included (read-only copies, signed zero canonicalised, content `==` and `hash`).
5. `BC` params are real numbers only, and NaN is refused at construction everywhere a value is made; the encoder's refusal is the backstop for a bare value.

**What remains for steps 6 to 8** (numbered as in the carve order above; the spec numbers them §1.6 phase-space functions, §1.7 question, §1.8 specification, with the source before the question):
- **The source** (`Symbolic` through `srepr`, `RegionwiseConstant`): it digests through the encoder as a frozen dataclass or a `ContentIdentity` value; a `srepr` string encodes, a live SymPy expression does not (`[M]` 2026-10-02 the encoder refuses a `sympy.core.add.Add` with `ContentlessError`, met inside `DiscreteMeasure.support`).
- **The question** (`Eigen(parameter)`, `FixedSource(source, point)`; spec S7.5): a `SpectralMap` holds lambdas, which have no content, so its part must be its `name`, as S7.5 says.
- **The specification** (materials, geometry or none, question, source; spec S8.2, S8.5): its digest is the encoder's over those parts, with its own RECORD fingerprint.
- **Open, measured 2026-10-02:** the directional `Quadrature` is a mutable dataclass (`frozen=False`) with a hand-written `_identity_key`, so `encode(Quadrature.gauss_legendre(4))` raises `ContentlessError`. Nothing in steps 6 to 8 keys on a quadrature (the specification is mesh-free and method-free); a persistent key that covers an S_N discretisation (P3) needs the quadrature moved onto the encoder first.
- **Outside the encoder, by scope:** the Sood registry's result cache keys on a SHA-256 of sorted JSON (`1` and `1.0` differ, `-0.0` and `0.0` differ, NaN admitted, `{1: x}` collides with `{"1": x}`; census §1a), with no production consumer. `[R]` It is a candidate for retirement once the reference cache serves its one test consumer.

**After review (2026-10-02; supersedes the details above where they differ).** The elegance review (S1–S6) and qa (D1–D3, G1, G3) changed the landed shape:
- **`FrozenMapping`** (`orpheus/numerics/content.py`) is the one frozen mapping value. `Materials.mixtures` and `BC.params` hold one, and `FaceLaws` subclasses it. The `MappingProxyType` stores and their two `__reduce__` methods retired.
- **Every value pickles through its constructor** (`ContentIdentity.__reduce__`), so a load re-runs the laws. A pickle written under an older schema fails with a `TypeError` naming the field. This is the mechanism behind the 2026-10-01 ruling for pickles; the persistent cache still keys on digests and stores data, never pickles.
- **The `Axis` content is its dataclass fields**, with `generator` declared `field(compare=False)`. That is the one spelling of "not content", and the hand-written `content_parts` lists retired.
- **The encoder refuses** a dataclass compared by identity, a mutable part inside an object, and complex values. Its 2**53 check no longer wraps on unsigned arrays.
- **NaN is also refused at construction** in `LawScaled.scalar` and in the three response `alpha`s. The signs are integers. `Mesh2D` takes the shared parsers.
- **Gates:** 262 rows. The S5.7 pins were re-pinned before the commit: the geometry and `Materials` digests moved with `FrozenMapping`, and mixture A did not.
- **Full suite** `tests/gates -m "not slow"`: 13 974 passed. The 4 failures are `test_write_guards`' worktree-path rows.

Issues: #553 (retire the string-equality arm), #554 (stale SN docs block), #555 (Quadrature content identity).

**Next: P1 step 6** (the phase-space functions, spec §1.6) per the carve order above, then the question (§1.7) and the specification (§1.8).


## ⏸ COMPACTION POINT — 2026-10-02, P1 after step 5 (supersedes the point after step 3c)

State: `main` `5abb3db8`, clean apart from `scratch/`, and CI `gates` is green. Merged since the point after step 3c, in order:
- the boundary-law detour: ERR-094 and ERR-095 at `8f9300b2`;
- the platform-independent quadrature remedy, at `a34f8c2c`;
- the reflective cleanup, at `99b9847d` (plan `boundary_law_ontology.md`);
- **P1 step 5, content identity**, at `ee9e8943` (code) and `a5113ac0` (docs), with its close-out at `5abb3db8`.

**Next: P1 step 6, the phase-space functions.** These are `Symbolic` (q(r, μ, φ) per group, stored through `srepr`) and `RegionwiseConstant` (one value per region and group, the isotropic piecewise-constant special case), with the anisotropy predicate. There is no projection; that is P4. Note the numbering: the spec swapped steps 6 and 7, so the spec's §1.6 is this step and the carve-order list above still calls it item 7.

Read, in order:
1. The two 2026-10-02 notes above: step 5 as merged, and "After review". The second carries what steps 6 to 8 need from the encoder:
   - a `srepr` string encodes, but a live SymPy expression is refused;
   - a `SpectralMap` holds lambdas, so its part is its name;
   - the `Quadrature` is not content (#555).
2. §1.6 of `.claude/plans/reference_p1_spec.md` (gates S6.1 to S6.6), and its §0 rulings 2 and 4 (the regions are the geometry's own indices; evaluating exactly at a breakpoint is refused).
3. `.claude/plans/posing_sequence.md`, "The ontology as it stands", and this file's "The posing sequence's seeds for P1". The 2026-09-29 ruling: the question values live in `numerics`, physics-free. Decide where the phase-space functions live by the same criterion; the W5 review placed the `Source` field protocol in `numerics`.
4. The content-identity section of `docs/theory/foundations/structured_geometry.rst` (`structured-geometry-content-identity`).

**First action: reconcile §1.6 against today's tree before any gate or code.** The spec was written 2026-09-25. Step 5's census refuted three of that spec's premises (`Partition` retired, the `Axis` encoder's `repr`, `_ManufacturedFaceInflow` hashing by id). Dispatch an explorer census for step 6, then the test-architect re-specifies §1.6 and measures each first red in a detached worktree. The questions the census owes:
- What source types exist today? Measure the boundary source `InflowSourceSpec`/`ConstantInflowSource`, the MMS sources in `derivations/`, and how `solve_sn_fixed_source` takes a source.
- Is the density convention (S6.5: a value Q gives `∫∫ q dφ dμ = Q`, per steradian) consistent with today's consumers? The spec cites a docstring applying `/W` twice.
- Does `sympy` import cost or layer placement constrain the home?

**Rulings that bind P1 from here on:**
- **Content identity (2026-10-02):** every new value type is a `ContentIdentity` frozen dataclass (`eq=False`), with "not content" spelled `field(compare=False)`, mappings held as a `FrozenMapping`, and NaN refused at construction. Its gates join the content-identity rosters (`tests/gates/<tree>/test_content_identity_*.py`; S5.9 reds on a new content class without a roster entry).
- **The schema-tag ruling (2026-10-01):** no per-class `__setstate__`.
- **Steps 1–3 rulings:** breakpoints rather than thicknesses; `None` is not a law; the regions are the geometry's indices.
- **The process:** W3, so the main agent writes the code and the test-architect writes the gates first, then qa and the elegance-enforcer review in parallel and the archivist writes the docs. The full suite runs before each merge, in a detached worktree (`.venv` linked; the 4 `test_write_guards` rows fail only there). Baseline: 13 974 passed at `ee9e8943`.

**Open issues from this stretch:** #551 (boundary architecture E), #552 (one base amplitude check), #553 (the string-equality arm), #554 (a stale SN docs block), #555 (Quadrature content identity). Architecture E (the boundary ontology) executes after #405, by the user's 2026-10-01 ruling.

## P1 step 6 opened — 2026-10-02: the census and the user's rulings

The census (explorer, `scratch/reference_architecture/p1step6/census.md`, against `main` `6a96cc3f`) measured the step-6 premises. What it established:
- No existing type duplicates the new ones. Every bulk source in the tree is an array bound to one discretised space; the one intensional source, `InflowSourceSpec` (`geometry/boundary/_source.py:56`), lives on the boundary trace and is the precedent.
- All 10 consumers read a volumetric source value Q as the angle-integrated rate. SN routes it through `space.section("angular")` (`angular_source_sink.py:193`); MoC divides by 4π by hand (`moc/core.py:187`); the 12 SN MMS cases divide by `sum_w` by hand; diffusion, CP and Peierls need no lift.
- SymPy is imported by 0 of the input-layer, transport and method packages; `numerics/manifold.py` imports it function-locally (3 sites). `pyproject.toml` lists SymPy only in the test and docs extras, although `Quadrature.gauss_legendre` needs it at runtime: a packaging defect, fixed in step 6 (SymPy becomes a core dependency, which `Symbolic` needs anyway).
- A region is the positional interval index of `StructuredGeometry` (`0 … len(mat_ids) − 1`); region is not material.
- `srepr` is seed-stable on SymPy 1.14.0 (seeds 0–3); stability across SymPy versions was not measured.
- Content identity of a `Symbolic` is by spelling: `(r+1)**2` and `r**2 + 2r + 1` digest differently, which only ever costs a cache miss.
- Two content types are never equal across types, so the spec's S8.2 leg "a `RegionwiseConstant` equals the `Symbolic` it lowers to" cannot hold as written; it is re-posed as "both lift to the same function".
- The docstring of `solve_sn_fixed_source` divides the source by W at `sn/solver.py:3253` and then says the source is already a per-ordinate density (`:3270`, `:3276`); the code is right and the docstring is fixed in step 6.

**Ruling 1 (the user, 2026-10-02): the 4π is the measure, not a convention.** The user asked whether "integrated over all directions" and "divided by 4π" were not both the measure, and whether production's hand-spelled lift was a defective operator. They are, and the operators already exist in `numerics` (`orpheus/numerics/operator.py`, `AxisRetractionOperator` and `AxisSectionOperator`):
- the retraction `R = space.retraction("angular")`, the fibre integral ∫ · dΩ;
- its section `E = space.section("angular")`, defined by `R ∘ E = id` (divides by the measure's mass, read from the frame's induced Gram);
- the retraction's Hilbert adjoint `R† = π*`, the plain broadcast; `R† = (Σw)·E`, measured in its docstring.

A per-region table is therefore a function on the angle-integrated space (r, g), not on phase space, and `RegionwiseConstant` is NOT "the isotropic special case of `Symbolic`". It enters phase space by one of two arrows, and the role picks the arrow: a source of rate Q enters the right-hand side through the section, `q = E(Q)`, so that `R(q) = Q`; a detector Σ_d is the functional `ψ ↦ ⟨Σ_d, Rψ⟩`, whose Riesz representative is `R†(Σ_d)`. The mass `R ∘ R†` (4π on the continuous sphere) is never typed. Source and detector stay one value type; the specification's field chooses the lift. On the Branch-1 side the lifts are realised with the continuous S² measure (SymPy), never by calling production's `E`.

**Ruling 2 (the user, 2026-10-02): the angular chart is the coordinate system's local frame.** The geometry is 1-D in P1. μ is Ω·ê_x on a slab and Ω·ê_r on a cylinder or sphere; φ is the azimuth about that axis, measured from a reference direction each coordinate system declares once (ê_z on the cylinder). `Symbolic` owns the symbols r, μ, φ with fixed assumptions and refuses any other free symbol, naming it. 2-D is out of scope until a 2-D geometry exists.

**Ruling 3 (the user, 2026-10-02): step 6's scope.** Step 6 defines the types and gates the two lifts (`R ∘ E = id` gives back Q; `⟨R†Σ_d, ψ⟩ = ⟨Σ_d, Rψ⟩`). The same branch re-spells the SN adjoint detector lift, today a hand-written `np.broadcast_to(sigma_d[None], …)` at `sn/solver.py:2927`, as the retraction's adjoint. Filed, not fixed: the MoC hand-written 4π (#556), and the derivations' hand-written measure masses (#557; Branch 1 should derive the mass from its own continuous measure).

Next: the test-architect re-specifies §1.6 of `reference_p1_spec.md` under these rulings and measures each first red in a detached worktree.

### P1 step 6 built — 2026-10-02, merged at `61f82a18` (code `4b084883`..`43067473`, docs `2d0f9539`)

Commits: `4b084883` (the pullback; the mint keeps the marginal's forms and refuses a form on the collapsed axis), `3d1dc28c` (the SN adjoint detector lift is `retraction("angular").H`, bit-exact), `c9b28def` (`CoordSystem.angular_chart`), `f7309da3` (`RegionwiseConstant`, `Symbolic`, the Branch-1 `angular_measure`, SymPy a core dependency), `43067473` (the qa and elegance review fixes; the module is renamed `orpheus/numerics/mesh_free_function.py`). The spec's build-time corrections are in `reference_p1_spec.md` §1.6, "Corrections found while building".

**Obligations this step hands to P4 (the projection onto a method's unknowns):**
- **The pushforward onto an orbit space.** A `Symbolic` is a density with respect to dΩ. On a rule whose ordinates are points of an orbit space (every 1-D rule: the sphere quotiented by O(2) about the polar axis, mass 2), its per-ordinate value is the integral over each orbit, ∫₀^{2π} q dφ, not q. That arrow is not built. Its factor is derived from the two measures, never typed (the elegance review, 2026-10-02, item 1).
- **Observability of the azimuth is the problem's symmetry**, not the coordinate system's: a 1-D slab cannot see φ, and a 2-D Cartesian mesh in the same coordinate system can. A chart flag for it was retired in the review. The projection should read it from the geometry's direction stabiliser, the group a quadrature already declares as its orbit space (`rules_1d.py`, `SPHERE.quotient(O2(axis))`).
- **One Branch-1 coordinate vocabulary (#557):** `derivations/.../transport_equation.py` declares `r` non-negative, which `Symbolic` refuses (it owns `r`, `mu`, `phi` as `real=True`). The MMS builder at `mms/sn.py:2601` declares `r` and `mu` positive.

**Issues filed by the step-6 review:** #556 (MoC's typed 4π), #557 (the derivations' typed measure masses; `angular_measure.py` is its seed), #558 (a composite's adjoint does not reach the leaf pullback), #559 (one real-number parser at L1), #560 (the iso + aniso source combine spells the section by hand, with W defined twice: `weights.sum()` against the frame's Gram entry; [M] its one production caller ran 387 455 times in the transport/SN/numerics gates).

**Step 6 close-out (2026-10-02).** Merged `--ff-only` at `61f82a18`. Full suite `tests/gates -m "not slow"` in a detached worktree of `2d0f9539`: 14 086 passed, 264 skipped, 55 xfailed, 4 failed (the `test_write_guards` worktree artefacts, 22 of 22 green in the main tree); baseline 13 974 at `ee9e8943`, so +112 rows, all accounted for by the new gates. The merge push's tip carried `[skip ci]`, which skipped the push's CI run; this close-out commit is pushed without it so `gates` runs on the merged tree. Process lesson: a plan-only `[skip ci]` commit must never be the tip of a push that carries code. `[M]` 2026-10-02, twice: the close-out commit `36b61af3` was meant to trigger CI, but its message BODY quoted the skip token while explaining it, and GitHub matches the token anywhere in the message, so that push ran no CI either. A commit meant to trigger CI never quotes the token, even in prose; CI on the merged tree was triggered by the commit after it.

Next: P1 step 7, the question (spec §1.7; S7.3 amended for the two step-6 types). Its open ruling: what `Eigen(c)` names and what `CriticalParameter` varies (waiting on the posing-sequence review). [Superseded 2026-10-02: the posing sequence ruled the question values, and the user ruled the three remaining gaps; see `reference_p1_spec.md` §1.7, "The step-7 rulings".]

## ⏸ COMPACTION POINT — 2026-10-02, P1 after step 6 (supersedes the point after step 5)

State: `main` `6ba7aab3`, clean apart from `scratch/`; CI `gates` run 37017185290 green on the merged tree. Merged since the point after step 5: **P1 step 6, the mesh-free functions and the retraction's adjoint**, at `61f82a18` (code `4b084883`..`43067473`, docs `2d0f9539`), closed out at `36b61af3` and `6ba7aab3`. Full-suite baseline (`tests/gates -m "not slow"`, detached worktree of `2d0f9539`, `.venv` linked): 14 086 passed, 264 skipped, 55 xfailed; the only failures are the 4 `test_write_guards` worktree artefacts (22 of 22 green in the main tree).

What step 6 landed (the record is the two notes above, "P1 step 6 opened" and "P1 step 6 built"):
- `AxisPullbackOperator`: `R.H` of `space.retraction(label)` is the closed-form pullback π*, minted by the retraction, `R.H.H is R`; the mint keeps the marginal's positioned forms (`FunctionSpace.without_axis` / `with_axis`, one overlay reader `_overlay_forms`) and refuses a form on the collapsed axis.
- The SN adjoint detector lift is `retraction("angular").H` (bit-exact).
- `CoordSystem.angular_chart` (polar axis, azimuth reference; the sphere has none).
- `orpheus/numerics/mesh_free_function.py`: `RegionwiseConstant` (on the angle-integrated space) and `Symbolic` (q(r, mu, phi) per group as srepr text plus the SymPy version, whitelisted parse, real scalar functions only, isotropy by substitution). Branch 1: `orpheus/derivations/common/angular_measure.py`. SymPy is a core dependency.

**Rulings that bind from here** (in addition to the content-identity ones of the step-5 point):
- **The 4π is the measure (the user, 2026-10-02).** A per-region table lifts by the role's arrow: a source rate by the section E, a detector by R†. Neither type carries a role or a density; the specification's FIELD is the role. No mass is ever typed.
- **The chart is the coordinate system's local frame** (the user, 2026-10-02); observability of the azimuth is the problem's symmetry, not the chart's.

**Next: P1 step 7, the question values.** ⚠ Spec §1.7 (`reference_p1_spec.md`, written 2026-09-25) is STALE against the posing ontology, and this is the first thing to reconcile, before any gate:
- §1.7 lists `Eigen(k)`, `Eigen(c)`, `FixedSource(σ)` and `CriticalParameter`. The posing sequence ruled on 2026-09-26 (`posing_sequence.md`, the rulings ledger) ONE eigen question, `Eigen(parameter)`, the parameter a closed set of declared coordinates of the system (`posing_sequence.md`, "The ontology as it stands": `CellCoefficient | GeometryExtent | NuclideDensity`, 2026-09-28; the k, c and α parameters are coefficients, an outer extent is a geometry coordinate). `CriticalParameter` is struck (the ruled text: it never existed in the tree, and nothing dissolves into anything), and the pencil and the spectral map are DERIVED from the parameter.
- The 2026-09-29 ruling: the question values live in `numerics`, physics-free, naming the parameter as an opaque key the specification resolves (the "contradiction to rule" in `posing_sequence.md` "The work, ordered by dependence", item 1, proposes `Eigen(parameter, point, mode)` and `FixedSource(source, point)`).
- S7.3 was amended 2026-10-02 for the two step-6 types: the adjoint fixed-source question's detector is a `RegionwiseConstant` (lifted by R†) or a `Symbolic`.

First action: a short explorer census on what the posing campaign (#522, units #523-#529, plan files `posing_u0`..`posing_u6`) has already landed of the question vocabulary (`numerics/posing.py`: `EigenPosing`, `SourcePosing`, `K_MAP`, `ALPHA_MAP`; any `Eigen` / parameter-coordinate type), then the test-architect re-specifies §1.7 against it and measures each first red in a detached worktree. Then step 8, the specification.

**Process lessons from step 6, for every later step:**
- A plan-only commit carrying the CI skip token must never be the TIP of a push that carries code, and a commit meant to trigger CI never quotes the token, even in prose (GitHub matches it anywhere in the message).
- A full-suite run is driven with its own log file per run, and stopped by PID (`pgrep -fl "pytest tests/gates"`), never by a path pattern: a relative-path launch evaded `pkill -f`, kept running while its worktree was moved, and interleaved with the next run's log.
- A gate written from the spec can be designed green: S6.19(b) re-synthesised two columns from an angle computed out of the same two columns. Every new gate gets its mutation arm before commit.

**Open issues from this stretch:** #551–#555 (from the step-5 stretch), #556 (MoC's typed 4π), #557 (the derivations' typed measure masses), #558 (a composite's adjoint does not reach the pullback), #559 (one real-number parser at L1), #560 (the iso + aniso source combine spells the section by hand; W defined twice; an SN-wide ULP re-baseline).
