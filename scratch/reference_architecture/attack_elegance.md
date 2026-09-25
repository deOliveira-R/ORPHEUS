# W5 design review, #405 step 2: adversarial attack on `.claude/plans/reference_cache.md` (elegance-enforcer)

Date 2026-09-24, tree `901f64ca` plus the uncommitted plan. Angle: the software architecture and the soundness of the machinery. No tracked file was edited. The probe scripts are in this scratchpad: `ee_probe1.py` (an AST census of memoisation, mpmath precision, environment reads, lambdify, numba and import-time calls) and `ee_probe2.py` (the blind spots of `sys.monitoring`, run under `python -O`, Python 3.14.3).

The theorem the design rests on: if a generator is deterministic, then its output is a function of (the code that executed, every non-code input). A trace-set identity is therefore sound exactly when every non-code input is in the key and every executed piece of code is in the trace. Every hole below is one violation of that premise, found in the live tree.

---

## Pass 1: adversarial (how would I break it; how would I make it 100× better)

### H1. In-process memoisation hides dependencies from the trace, and the plan's own intermediate-stage cache is an instance of the same hole

- `[M]` `ee_probe2.py`: `gen_A` calls an `lru_cache`d `kernel`, which calls `helper`. The first trace records `kernel` and `helper`. The second generator `gen_B` hits the warm cache, and its trace is `{gen_B}` only. Neither `kernel` nor `helper` is recorded. An edit to `helper` would then leave `gen_B`'s entry validating, and the stale entry would be served silently.
- `[M]` `ee_probe1.py`: 18 decorated memos under `orpheus/` (6 in derivations, 10 in numerics, 1 in sn, 1 in transport). Among them are `derivations/discrete/sn/face_transmission.py:557` `cell_response @cache` and `:741` `transmission_spectrum @cache`, and `face_transmission` is one of the plan's five listed generators (plan, the table under discussion 2). There is also 1 hand-rolled memo: `singular_eigenfunction/core/half_range.py:88-110` (`global _X_CACHE`). On a hit it returns before calling `atalay_X_function`, so `x_function.py` leaves the trace. This is the Sood/Atalay family, which is a P4 family. Total: 19 sites.
- `[R]` The plan's own "the same mechanism caches ... intermediate stages (the Peierls volume kernel)" (item 8) creates the same blindness. A parent generator that gets a disk hit on the kernel K_vol never executes the kernel's builder, so the parent's trace omits it.
- `[M]` `ee_probe2.py`: returning `sys.monitoring.DISABLE` (the cheap way to record each code location once) also blinds every later recording in the same process, until `restart_events()` is called. So each recording must restart events, and nested recordings need a recorder stack.
- `[R]` In a full pytest process, which memos are warm depends on the order the tests ran. So the recorded trace is order-dependent, and the plan's X1 witness ("a mutation inside a recorded function turns a hit into a miss") would pass in isolation and not in the suite.
- **Rewrite risk: HIGH if P3 ships a flat trace set.** The change now: one primitive, a *traced memo*. An entry stores its own trace plus the keys of the child entries it consumed (a Merkle DAG, as in Nix and Bazel). Validation recurses through the children. On a hit, the in-memory tier replays the stored child trace into every active recorder. A gate then forbids raw `functools.cache`/`lru_cache` and hand-rolled memos in first-party code, and the 19 sites migrate. The payload cache, the certificate latch (H4) and generator-keyed test verdicts (C3) all become clients of this one primitive. The plan currently names three mechanisms (item 8, item 5, and C3, which item 8 does not mention).

### H2. The specification cannot pose two question kinds the tree already answers

- **Critical-configuration questions.** In these the geometry is the unknown.
  - `[M]` Four generators expose `solve_critical`: `fn_method/moment_space.py:227`, `galerkin_spectral/basis_space.py:616`, `singular_eigenfunction/spectrum.py:673` and `trajectory_resolvent/billiard.py:316`.
  - `[M]` `CriticalSolution(eigenvalue, eigenvalue_kind: str, parameter_value, parameter_kind)` is at `derivations/common/solution_types.py:127`, with 88 hits in 12 files.
  - `[M]` There are 32 published `critical_dimension_mfp=` values (25 in `la13511.py`, 7 in `atalay1997.py`).
  - Item 1's rule "a source means fixed-source, none means eigenvalue" cannot spell "find r_c such that k = 1".
  - `[R]` The plan's defect "`La13511Case.to_geometry()` reads the size from the truth record, so posing the problem reads the answer" is misdiagnosed. It is the missing question kind. A k-eigenvalue problem at the published r_c is a second question whose input comes from another solution's answer. The certificate model has no way to carry that input's bound (the printed precision of r_c, times dk/dr) into k's bound.
- **Discrete-operator questions.**
  - `[M]` `orpheus/derivations/discrete/` holds 12 modules (4451 lines), imported by 17 test files. `face_transmission.py:1-30` is the algebra of record for a Petrov-Galerkin scheme's transmission.
  - Its question *is* a discretisation, and the specification is discretisation-free by design (Oberkampf and Trucano), so neither the specification nor "one registry keyed on the specification" can hold it.
- **Rewrite risk: HIGH for P1** if the specification ships as one concrete type with the question derived from whether `source` is `None`. The change now:
  - The question is an explicit closed sum on the specification: `Eigen(spectral map)`, `Source(q)`, `CriticalParameter(parameter, target)`.
  - The registry is keyed on a `Question` protocol (content digest plus facets), with the physical specification as one instance and a discrete-operator question as another. Alternatively, declare the discrete family out of the registry's scope (cache only), in the plan.
  - Deriving the question from None-ness is also the #465 `is_eigenvalue` tag again, spelled implicitly, and every consumer would re-ask it.

### H3. State from import time and C-level state is invisible to a trace that starts later

- **The environment variable.**
  - `[M]` The only environment read in `orpheus/` is `peierls_nystrom/cases.py:85-87`: `_SLAB_VIA_UNIFIED = os.environ.get(...)`, evaluated at import.
  - `[M]` `ee_probe1.py`: an audit hook sees no event for `os.environ.get`, `os.getenv` or `in os.environ`, and does see `open` (including inside `np.load`).
  - So the plan's third X1 witness ("flipping `ORPHEUS_SLAB_VIA_E1` turns a hit into a miss") fails under the design as written. The module's text does not change, and nothing records the value.
- **Import-time computation.**
  - `[M]` 15 module-level assignments call a first-party callable imported from another module, for example `data/group_structures/wims.py:148` `WIMS_69 = EnergyGrid(np.asarray(edges(...)))`. There are 124 module-level assignments calling non-stdlib callables in all.
  - Code that runs at import, before the trace starts, is not recorded, and "the module-level definitions of their modules" hashes text, not the callee bodies.
- **Constants read without code.**
  - `[R]` A constant that is read (an attribute or subscript) from a module none of whose functions ran puts no module in the trace.
- **Code generated by `dataclass`.**
  - `[M]` `ee_probe2.py`: a frozen dataclass's `__init__` is traced as `('<string>', '__create_fn__.<locals>.__init__')`. It has no file and no class name.
  - Re-hashing by source then either fails on almost every entry (the cache never hits) or skips the code. Skipping is unsound when the class's module is otherwise unrecorded, for example when a field default changes.
- **The ambient mpmath precision.**
  - `[M]` Three sites set the global `mp.dps`: `singular_eigenfunction/origins/cylinder_derivations.py:399`, `trajectory_resolvent/origins/specular/greens_function_slab.py:456` and `tests/gates/transport/spatial/test_face_transmission_symbolic.py:560`.
  - `[M]` Of the 10 mpmath-using files in `orpheus/`, the Peierls ones make about 33 to 44 calls against 6 `workdps` blocks each. `[R]` Code that runs outside those blocks depends on test order.
- **CPU and thread state.**
  - `[R]` The BLAS thread count and the CPU-specific kernel dispatch are read by C code and are invisible. A platform tag of `ubuntu24-x86_64` does not separate the different CPU models of hosted runners.
- **Rewrite risk: medium, additive if decided now.** The change now:
  - Start monitoring at session start, so that each module's import-time trace is part of that module's identity.
  - Snapshot the VALUES of module globals derived from the environment into the key, or ban import-time environment reads (a gate; there is 1 site).
  - Attribute code from `<string>` to its defining class when it is recorded.
  - Run every generation under a pinned mpmath context, and put `mp.prec` in the key.

### H4. The certificate latch is keyed on the wrong identity

- `[R]` A failed corroboration between solutions A and B latches both as `Invalid`, each under its own generator identity (item 5).
- If B is fixed, B re-qualifies, but A's key has not changed, so A stays `Invalid` forever. Discussion 3 states that the other solution's change must unlatch A. The same applies to a fix in the trigger's own code.
- So the certificate is a function of the whole certification computation: both generators, the anchors and the check.
- The change now: make the certificate an H1 traced memo *of the certification run*. That is sound for free, and it is additive.
- A trip from the cold rebuild ("the rebuild does not reproduce the payload") is not a deterministic function of any key. The latch lives in a gitignored cache that is evicted after 7 days unused, so it can evaporate. So "reproduce" must mean agreement within a fraction of the certified bound, never bit equality. The platform tag then becomes unnecessary for any entry whose bound exceeds the platform spread in units in the last place (ULP).
- `Withdrawn` is keyed by a name. A renamed generator silently un-withdraws, so the CI check must also assert that each tag's target exists.

### H5. Placement against `tests/gates/test_layer_imports.py`

- `[M]` The linter at `test_layer_imports.py:55-71` allows the edges geometry→data, data→geometry and derivations→{numerics, geometry, data}. So the placement of the specification is legal.
- **The `Symbolic` source.**
  - `[M]` sympy is not a production dependency: `pyproject.toml:16-23` lists numpy, scipy, matplotlib, pyXSteam, h5py and pint, and sympy is only in the `test` and `docs` extras (lines 55 and 70).
  - `[M]` 0 files in geometry, data, transport or sn import sympy. numerics imports it lazily in 3 function bodies (`manifold.py:1657, 1784, 1882`).
  - So a `Symbolic` arm at the input layer that production "projects" puts an optional dependency on the path of the system under test.
  - The change now: the input-layer `Source` is an L1 field protocol (it evaluates and integrates over a region). The sympy-backed realisation lives at L0 and carries its `srepr` for identity.
- **The comparison verb and its observables.**
  - `[M]` The production results are three unrelated types: `sn/solution.py:375` `Solution`, `diffusion/solver.py:392` `DiffusionResult` and `homogeneous/solver.py:77` `HomogeneousResult` (CP and MoC were not counted).
  - So "a production `Solution`'s observable" (item 6) is not one type.
  - The observable must be ONE object that both sides evaluate. `numerics/functional.py:147` `Functional` is the existing L1 contract. Otherwise the verb branches on observable tags such as `"k"` and `"cell_avg"`.
- **The persistent key.**
  - `[M]` `Mixture.__hash__` (`data/macro_xs/mixture.py:201`) is Python's salted `hash`: two processes printed different values for `hash(b'abc')`.
  - The key on disk needs a stable digest derived from `_identity_key` (`:180`), the one definition, never a second serialiser.

### 100× reframe `[R]`

This is not a reference cache. It is *incremental verification*: a build system over pure computations (Mokhov, Mitchell and Peyton Jones 2020, dynamic dependencies). Payloads, certificates and the verdicts of generator self-tests are three kinds of memoised, traced, content-keyed artefact in one DAG. Designed as one primitive, the "third row" of the distinction table (5 of the 8 slowest files test the generator itself) is handled by the same mechanism, not by an unnamed plugin.

---

## Pass 2: re-evaluation

- **H1 survives, and is the top finding.** It is the only hole that serves a stale answer silently on the tree as it stands (19 memo sites, one of them in a listed generator).
- **H2 survives.** The measured counts (4 generators, 32 published values, 17 test files) are shipped consumers, not speculation.
- **H3 survives, narrowed.** The environment part is 1 site, but it is the plan's own witness.
- **H4 survives, but its size is modest.** It dissolves into H1's primitive.
- **H5's `Symbolic` point survives.** The layer legality attack is withdrawn: the linter allows the edges.
- **The type split `PublishedSolution`/`ReferenceSolution` is withdrawn as an attack.** It is capability-based (Pattern 4), as the plan argues. What survives is a CONCERN: `PublishedSolution` carries no state, yet "a failed corroboration makes both `Invalid`", and published tables have errata. The published side needs at least `Withdrawn`.
- **The infinite medium as "absent geometry" is a NIT now and a CONCERN later.** An `Optional` field forces `is None` branches at every consumer (anti-pattern #7). Make it a value, `InfiniteMedium`, the total quotient.

## Twin concepts

1. `CriticalSolution` and `FluxSolution` (`solution_types.py:127, 213`) are the tree's existing reference-solution types, and they are absent from P4's retirement list (plan line 331). `eigenvalue_kind: str` is a stringly-typed tag (anti-pattern #4).
2. `NotYet(issue, reason)` has the same shape as `Withdrawn(reason, issue)`. `Certified(bound, by)` has the same shape as a per-observable bound together with its establishment (`numerics/outcome.py:78-109`). The plan should decide one or two types by X4, and not by name.
3. The plan's "observable" and `numerics.Functional` / `ReactionRateFunctional` (`transport/reaction_rate_functional.py:76`).
4. "Solution" names four things: CLAUDE.md's `Solution` five-tuple (SN only, in fact), `PeierlsSolution`, `CriticalSolution`, and the new pair. This is a NIT.
5. `TransmissionSpectrum` (`derivations/.../face_transmission.py:720`) and `FaceTransmissionSpectrum` (`transport/spatial/scheme.py:457`) both carry `@cache`. They are a pre-existing pair and relevant only to H1.
6. Two `Provenance` classes, `ProblemSpec` and `StructuredGeometry`: already on the plan's list.

## Refuted attacks (refuted FOR soundness of the trace identity; the fact each establishes)

- **numba, ctypes or cffi, and first-party C extensions:** 0 imports under `orpheus/` `[M]` (`ee_probe1.py`).
- **Code generated by sympy `lambdify`:** 0 calls in `orpheus/` `[M]`.
- **First-party `sympy.Function` subclasses whose `eval` sympy's global cache would hide:** 0 `[M]`.
- **Closures, `functools.partial`, generators and callbacks from C (`scipy.integrate.quad`):** `PY_START` fires for Python code whatever the caller. In `ee_probe2.py` the uncached nested `helper` was recorded. `[R]` for partial and generators.
- **Withdrawn Peierls code executed by a non-withdrawn generator:** 8 files outside the package name `peierls_nystrom`, and every hit is docstring text. There are 0 code imports `[M]`.
- **Equal source in two modules, renames and moves:** the key includes the file, so there is no aliasing. A rename or a move is a miss, which costs time and does not serve a stale answer `[R]`.
- **Randomness and threads in derivations:** 0 `[M]` (plan, the worklist section).
- **Environment reads other than the one found:** 1 site in `orpheus/` in all `[M]`.

NEEDS:
- (a) A user ruling on whether critical-dimension and discrete-operator questions are in the registry's scope, or are declared out of it (H2).
- (b) Measure whether any generator's result actually depends on ambient `mp.dps`: run each generator after `mp.dps = 50` and compare the results bitwise (H3).
- (c) Measure whether GitHub Actions runners' CPU models vary within one runner label (H4).
