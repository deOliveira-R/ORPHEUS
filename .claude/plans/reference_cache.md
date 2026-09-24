# Fast verification: a semi-analytical reference is recomputed only when its generator changes (#405, step 2)

Opened 2026-09-24. Status: **living plan, ontology being searched** (`plan-authoring` §0). Nothing is built. Implementation waits until the user rules the plan polished. The discussion starts after context compaction: *"Start drafting a living plan based on the context you have about this and prepare for context compaction. We will start this discussion after regrounding."* (the user, 2026-09-24).

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
  - `orpheus/derivations/_richardson_cache.json`, gitignored (`.gitignore:9`). Its producer and consumers are not yet read `[R]`.
- **In-process memoisation** exists in places (`functools.lru_cache` in `derivations/common/shifted_legendre.py`, `flat_source_cp/geometry.py`, `discrete/sn/dsa.py`, and `tests/gates/derivations/test_peierls_rank2_bc.py:976`); it helps within one run only.

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

## ⏸ COMPACTION POINT — 2026-09-24

Resume here after regrounding. Nothing is built. First re-read this file, then `.claude/plans/test_runtime_405.md` ("Step 1, measured" is the worklist). Then bring the eight questions to the user, one discussion at a time: this is ontology search, and the object (question 1) and the identity (question 2) come first. Before arguing any scope, measure the premise (`process-discipline`):
- read the two or three largest files on the worklist, to confirm which row of the distinction table each belongs to;
- read `_richardson_cache.json`'s producer.
