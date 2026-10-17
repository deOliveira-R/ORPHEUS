# The verification suite is one causal DAG of sharpened tests — a living plan

Opened 2026-10-16, from the profile of the bypass drop (#405 P3 step 2c; `characteristic_reference_architecture.md`, its last ⏸ point). Status: **the ontology is being searched** (`plan-authoring` §0). Nothing here is implemented. Implementation starts when the user rules the plan polished.

## Goal, in the domain's terms

Every verification test is a node of one directed acyclic graph of increasing complexity:
- each test rests on simpler tests whose claims it assumes;
- a red node invalidates the readings of the tests above it, and only those;
- no test is vacuous: for each one, some realisable defect turns it red.

The graph is also the schedule. A test runs only when its supports are green. A value several tests rest on (a reference, an oracle) is computed once, as its own node, and read by the tests above it.

## The user's instructions (verbatim, 2026-10-16)

> "Do we have a plan to improve the verification tree? Because concept 2 seems like would match a general path of improvement of the verification tree. We still need the test architect to essentially go through the verification tree and redesign it in such a way that follows its current directives, which includes a causal DAG hierarchy of increasing complexity (and no vacuous tests). Concept 2 seems like would fit this very well. After that is maybe a discussion on how to cache it and if there are fast ways for the cache to be accessed in parallel for all tests that do not share a dependency do that the test wave moves fast (if this is not already done)."

> "Regarding the parallel run, that's true today because if we don't have a DAG and all tests run independently, then what I think typically happened was that the fast tests run and finished fast, and eventually we were running any number of slow tests simultaneously and the memory footprint increased accordingly. The DAG should at least help manage the test fan out cone. It's hard to say if the end result would still be that at the end of the cone would be the chronically slow tests, and if they would all meet and result in the same "any number of slow tests running simultaneously", but maybe we can find a smart way of managing the slow ones."

> "I think trigger 2 is generally true, but in a complete redesign we would be deliberately ignoring it. Now it is important that a redesign would not throw everything out of the window. It would systematically look at the tests, identify how they fit into a hierarchy and how the hierarchy build up, find vacuous tests (maybe because it was a bad test, or maybe because it would be a complete duplicate in functionality and behaviour of another test that is stronger -- and the stronger twin stays) and find tests that needs improvement, so that they are strengthened and depend on a simpler node. So this is more like picking raw material, and beating it into a sword, than throwing away a sword we don't like and making a prettier one. The existing tests are raw material, and they need to be purified, refined, improved, beat into shape, sharpened and polished."

"Concept 2" is the second of four paths named after the profile, below: an oracle value computed once and stored, instead of recomputed by every test that reads it.

## Rulings ledger

1. **Refine, never replace** (the user, 2026-10-16). The existing tests are raw material. Each test is placed in the hierarchy, then kept, strengthened, re-based onto a simpler node, merged into a stronger twin, or retired as vacuous. Each verdict states its reason. A rewrite from scratch is refuted, below.
2. **This campaign sets aside the per-capability trigger and keeps the per-capability method** (the user, 2026-10-16). The old rule (`test_dependence_and_claim_kind.md`, ⏸ 2026-09-22, owed item 4) reads: "Migration is by capability, starting where a ladder is written (the test-architect's specs), never a bulk annotation of the existing suite."
   - **The trigger** ("a capability gets its ladder when work touches it") stays the rule outside this campaign. The campaign deliberately ignores it and goes through every capability.
   - **The method** stays in force: every `rests_on` edge is designed by reading the two tests, one capability at a time, never added by a pattern. The reason is measured: one blanket module-level `pytestmark` minted 124 of 243 refuted coverage claims (#387).
3. **A duplicate gives way to its stronger twin** (the user, 2026-10-16). Two tests with the same function and behaviour are one claim: the stronger stays, the weaker retires. "Stronger" is decided by the mutations each turns red (X1), not by name or size. `[R]` the criterion's form is a hypothesis until the census applies it.

## Vocabulary (one meaning per word)

- **Node:** a verification test (a pytest function, with all its parametrised cases). Whether the node should instead be a claim, meaning a test paired with one label it verifies, is open (candidate ontologies, below).
- **Edge:** `@pytest.mark.rests_on(<node id>, ...)` on the upper test. It means three things: the upper test assumes the lower is green (dependency); a red lower test invalidates the upper reading rather than failing it (boundary); a failure propagates only along edges (causality). The definition is the test-architect's (`docs/development/agents/test-architect.md`, item 4).
- **Vacuous test:** a test no realisable defect turns red. It comes in two kinds:
  - **blind:** no mutation of what it claims to verify reddens it (X1);
  - **dominated:** a stronger twin reddens under every mutation that reddens it.
- **Strengthen:** change a test so that more defects turn it red, or the same defects with a clearer failure. Typical moves: a premise asserted, a band tightened to its measured floor, an input outside the family a closure is exact on (`vv-principles` #7).
- **Re-base:** move a test's dependence onto a simpler node. A test that assumes a whole solve, where a term or an operator would do, rests on that simpler node instead.
- **Oracle node:** a reference value (an mpmath route, a closed form, a semi-analytical solve) that tests compare against. Today it is computed inside each test that reads it.

## What already exists (pointers, not copies)

| Artefact | What it holds |
|---|---|
| #358 (open) | A red should invalidate its dependent cone, not stop the run; the suite needs a dependence DAG and a claim-kind axis. |
| `test_dependence_and_claim_kind.md` | The Nexus side: the graph built from `rests_on`, a falsifier per edge (the lines the lower test runs are run by the upper one too), cycles refused when tests are collected. |
| `test_architecture_redesign.md` (2026-08-18, a proposal) | "A test declares what it is about"; ruled as step 3 of 3. |
| #381, #387 (open) | The coverage-claim ledger: 99.2 % of claims could not be adjudicated at first measurement. |
| #405 (open) | The slow tier, and the traced memo that came out of it (`reference_cache.md`). |
| `vv_suite_layout.md` (landed) | The `tests/gates/` layout; verification as a DAG, validation as independent cases. |
| `docs/development/agents/test-architect.md` | The ladder-is-a-DAG directive and the `rests_on` duty. |

## Measured state

All counts below are `[M]` 2026-10-16 on `main` `417f64f3`, by `git ls-files` and `git grep`.

- **Test files:** 625 files `tests/gates/**/test_*.py`. Per package:

  | Package | Files |
  |---|---|
  | `sn` | 232 |
  | `numerics` | 109 |
  | `derivations` | 109 |
  | `transport` | 50 |
  | `geometry` | 36 |
  | `data` | 20 |
  | `cp` | 15 |
  | `mesh` | 11 |
  | `diffusion` | 9 |
  | `homogeneous` | 8 |
  | other | 26 |

- **Suite size:** 17 182 not-slow tests passed on `bb2835d4` (`-O`, serial, 55 min).
- **`rests_on` today:** 603 marks in 96 files. Per package:

  | Package | Files with a mark |
  |---|---|
  | `derivations` | 23 |
  | `numerics` | 18 |
  | `geometry` | 16 |
  | `sn` | 13 |
  | `mesh` | 9 |
  | `reference` | 4 |
  | `cp` | 4 |
  | `specification` | 3 |
  | `homogeneous` | 2 |
  | `data` | 2 |
  | `transport` | 1 |
  | other | 1 |

  So 13 of the 232 `sn` files carry an edge.
- **Where a slow run's time goes:** one population, the 40 slow rows of 10 characteristic files, a warm store, `-O`, serial; `scratch/characteristic_architecture/p3/c3/drop/timing/profile/`.
  - **Wall time:** 942 s.
  - **mpmath quadrature:** 82 % of the profiled time (1060 of 1291 s under cProfile). All of it is called from the tests' own oracle modules (`tests/gates/derivations/_characteristic_point_mp.py`, and `_polar_column_errors` in `test_characteristic_lift.py`).
  - **By test name:** the rows that call those oracles take 764 s of the 942 s.
  - **Other parts:**

    | Part | Time |
    |---|---|
    | the characteristic `point_flux`, computed in the test process | about 85 s profiled |
    | the PS1982 Nyström reference | about 57 s |
    | the memo's own lookups | about 10 s |
- **mpmath backend:** `[M]` 2026-10-16, mpmath 1.3.0, a throwaway virtual environment, `scratch` log `gmpy_route.log`. The project's environment runs the pure-Python backend. With gmpy2 2.3.2:

  | Route | Before | After | Speed-up |
  |---|---|---|---|
  | frozen route F1 | 11.1 s | 7.8 s | 1.4× |
  | frozen route F4 | 26.7 s | 7.7 s | 3.5× |
  | frozen route F2, thick wall | 241 s | 30 s | 8× |

  The returned floats are identical. A pytest battery was running on the machine at the same time. F1 and F4 were interleaved, three repeats each; F2 ran once each.
- **What the memo treats as first-party:** any source outside the installed-package and standard-library directories (`_PinnedSourceFinder.vendor` in `orpheus/numerics/_traced_memo_boot.py`). So a module under `tests/` would be pinned like `orpheus/`. Not tested: whether a generating child process can import a function from `tests/`.
- **Parallel runs:** the suite runs serially.
  - pytest-xdist is pinned below 3.8: 3.8.0 deadlocks on Python 3.14.3 (`pyproject.toml`).
  - On 3.7, `tests/gates/sn -n auto` used about 21 GB of memory summed over workers.
  - Nothing schedules tests by the DAG.

## Candidate ontologies (each a hypothesis)

**A. The node is a test function.** This is the current spelling of `rests_on`.
- *Makes unspellable:* an edge to a test that does not exist; the edge resolves to a collected id.
- *Leaves welded:* a test that verifies two labels is one node, so an edge cannot say which of its claims the upper test assumes. The bypass drop's re-drop found this case: A1 carried `characteristic-door-reciprocity` (red under its mutation) and `characteristic-adjoint` (green under every mutation). It was one node holding one real claim and one empty one.

**B. The node is a claim: a test paired with one label.**
- *Makes unspellable:* an edge resting on a claim the lower test does not make.
- *Leaves welded:* nothing found yet.
- *Costs:* `rests_on` would name a claim, not a test id. A second spelling of the edge would then exist next to the 603 current ones, unless they migrate.

**C. Oracle values are nodes of their own, beside the tests** (concept 2).
- *Makes unspellable:* a test silently recomputing a reference that another test already validated.
- *Leaves welded:* still to be found.
- *What it means:* the oracle's own validation (for example the slow twin that recomputes a frozen literal) becomes that node's test. Every reader rests on the node and reads its stored value. The traced memo's manifest already invalidates a stored value when the code under it changes. That is the same dependence graph seen from code instead of from tests, so the two invalidation rules would compose on one graph `[R]`.

## Refuted candidates

- **Annotate the existing suite in bulk:** refuted, for the question "how do the edges land". An edge is a claim; a pattern mints claims nobody checked (#387, measured above).
- **Rewrite the suite to the directives from scratch:** refuted, for the question "how is the tree reshaped" (the user, 2026-10-16). It throws away the evidence each test carries: its measured bands, its first reds, its catalogued catches. Refinement keeps that evidence and sharpens it.

## Proposed phases (`[HYPOTHESIS]`, 2026-10-16; each is re-scoped by the phase before it)

1. **Measure before scheduling.**
   - Collect per-test duration and peak memory over the whole suite. Duration comes from `--durations`. Peak memory needs a sampler: `ru_maxrss` is the process high-water mark and cannot attribute memory to a test.
   - This decides whether the 21 GB came from slow tests running together or from each worker's import footprint.
   - For the slow tests, also measure S(p) and C(k), defined in the section on running slow tests below.
   - Also take the current `rests_on` graph as Nexus reads it: its components, its cycles, and its nodes with no edge.
2. **Census of the tree, by the test-architect, capability by capability.** It is the raw-material assay. For each test it records:
   - its rung and the nodes it should rest on;
   - its verdict: keep, strengthen, re-base, merge into a stronger twin, or vacuous (blind or dominated), with the reason;
   - the oracle nodes it recomputes.

   The census sizes the refinement before any rewrite. Its order, its schema and whether its verdicts are proven by mutation or only flagged for a later mutation pass are open (below).
3. **Refinement, capability by capability.** Apply the census verdicts. A retirement follows `retirement-audit`; a strengthened test lands with the defect it now catches (`plan-authoring` §6c).
4. **Oracle nodes as stored references** (candidate C). The prerequisite is a test of whether the memo can generate from a function under `tests/`.
5. **Scheduling on the DAG.** Three mechanisms, each `[R]` until phase 1 measures:
   - start a test only when its supports are green;
   - bound concurrency by a memory budget, from each test's measured peak, not by a worker count;
   - give each test its thread count from its own scaling curve, running a slow test alone when it uses every core (H-SEQ, below);
   - start the longest ready tests first, so slow tests start early beside fast ones instead of piling up at the end.

   The DAG's critical path, its longest chain of dependent tests, is the wall-time floor. A red support invalidates its cone rather than stopping the run (#358).

Sizing `[R]`: phases 2–3 are several sessions per large package (`sn` alone holds 232 files).

## Run slow tests one at a time with every core, or several at once? (the user, 2026-10-16)

> "There is also potentially some advantage in the slow tests being run sequentialy, which has to be measured. In principle, the production code should be mostly computation bounda, rather than memory bound (this has to be measured, but I think it is a strong hypothesis), which means that it might make more sense to execute each one sequentially giving all the computational power to it unless we hit a parallelism scaling problem (which requires its own measurement). The production code algorithms might also not be highly parallelised at this point, but we are certain to invest time and effort in it, with the intention of seeking vectorization, batching and parallelism (and whatever else technique you suggest)."

**H-SEQ** `[HYPOTHESIS]`: slow tests are bound by computation, not by memory bandwidth. If so, running them one at a time, each given every core, finishes sooner than running several at once with one core each. The exception is a test whose own parallel scaling is poor.

The choice is per test, decided by two curves, each measured:
- **Its own scaling, S(p):** the test's wall time alone at p = 1, 2, 4, 8 threads, through `OMP_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS` and any parallelism of its own.
  - If S(p) is close to p, the test uses the cores, so it runs alone with all of them.
  - If S(p) is flat, extra cores are wasted on it, so it runs beside others.
- **Its contention, C(k):** the slowdown of one copy when k copies run at once, each on one thread.
  - If C(k) is close to 1, the test is bound by computation and does not compete.
  - If C(k) grows with k, the test competes for memory bandwidth or cache, and running copies together costs more than it saves.

**The answer depends on the code's present state.** The profiled population is the 40 slow characteristic rows. 82 % of its time is pure-Python mpmath, which runs on one thread whatever `OMP_NUM_THREADS` is, so S(p) is flat for it today. `[M]` the profile above; the flatness `[R]`, from the backend. For that population, the only levers are concurrency and the backend (gmpy2). A production solve whose time is BLAS (`numpy.linalg.eigh` took 50 s of self time in the same profile) can have S(p) above 1.

So the scheduler keeps a resource profile per test: duration, peak memory, S(p) and C(k). Since the production code will be parallelised over time, the profile is re-measured as the code changes, and the schedule follows it. It is never fixed by hand.

**One constraint from the memo.** A forked child crashes when the parent called BLAS with several threads before the fork (`[M]` 2026-10-16, the fork-safety protocol of #602, `reference_cache.rst`). A test process given many BLAS threads must therefore not fork memo children after its first BLAS call. Today the template forks from a process that has only imported, which avoids this. A schedule that gives tests threads must keep that.

**The production techniques the profile would steer** (each `[R]`; where it applies is decided by measurement, not by this list):
- **Batching:** one call over many lines, angles or groups instead of a loop of scalar calls. `[M]` the chord kernel: 229 µs for one line, 0.33 µs per line in batches of 10 000 (`characteristic_reference_architecture.md`, P2).
- **Vectorization in numpy:** whole-array expressions in place of Python loops.
- **A compiler for the hot loops** (numba): for a recurrence that cannot be written as array expressions, such as a sweep along a direction.
- **Threads over independent work:** the directions of an SN sweep inside one octant and the groups of a Jacobi-style group iteration are independent, so they can run on separate cores. Gauss-Seidel over groups orders them, so it is not.
- **Faster convergence:** Krylov acceleration of the sweep (#200) and diffusion synthetic acceleration cut the number of sweeps, which outweighs any speed-up per sweep.
- **The array-library lowering** (#394, the layout freed for JAX): compiled, fused kernels and, later, a GPU.

The order follows phase 1's measurements: the techniques go where S(p), C(k) and the profile say the time is.

## Open questions for the user

1. **The node:** a test (A) or a claim (B)? The re-drop's A1 case argues for B. The 603 existing edges argue for A plus a convention.
2. **The census order:** by layer (`numerics` and `geometry` first, since everything rests on them), or by the sharpening order (`sn`, `diffusion`, `homogeneous` first)?
3. **gmpy2:** adopt it as a project dependency? It is below the trusted-library line, and gave identical floats on three routes.
4. **Census verdicts:** proven by mutation per test, which is expensive and exact, or flagged by reading first with a mutation pass per capability after?

## Side findings to place

- **The Nyström row's marker** (`peierls-greens-slab-architecture`) runs only under `ORPHEUS_RUN_WITHDRAWN=506`. In the canonical run it covers nothing. Its margin is 7.3e-5 against a 5e-5 band, for an attenuation error of 1e-3. Record it under #506. (`scratch/characteristic_architecture/p3/c3/drop/redrop/verifies_verdicts.md`)

⏸ **Next:** the user's answers to the open questions; then phase 1.
