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

4. **The node is a claim; a test is the pass or fail of one claim** (the user, 2026-10-16): "Node is a claim, test is pass or fail on the claim. What is a test without a claim? What are we even testing for? Tests hierarchy is essentially claim dependency." Candidate B is ruled; candidate A is refuted for this question (below). An edge says that one claim depends on another.
5. **The census goes by layer** (the user, 2026-10-16): "census by layer, because layers build up the verification causality cone. Transport methods are subcones at the base of the big cone." The order follows the import layers (`docs/architecture/layering.rst`): `numerics/` first, then the input layer (`geometry/`, `mesh/`, `data/`), then `transport/`. Each L3 method package (`sn/`, `diffusion/`, `homogeneous/`, `cp/`, `moc/`, `mc/`) is a sub-cone resting on that base, and `derivations/` sits beside them, below L2.
6. **gmpy2 is adopted** as mpmath's backend (the user, 2026-10-16). It lands after the bypass drop's merge gate, so the running suite does not see two backends. [LANDED; `[M]` 2026-10-18 `mpmath.libmp.BACKEND == 'gmpy'` in `.venv`, `gmpy2>=2.3` in `pyproject.toml`.]
7. **Claims have both kinds of home, chosen per claim** (the user, 2026-10-16): "Probably there are claims that would fit the first candidate and claims that better fit the second one." The two candidates of open question 5 are both kept; which one a claim uses is decided per claim, and the rule for deciding is still to be found.
8. **A concept is proved by several claims, and claims compose up the cone into a method's verification certificate** (the user, 2026-10-16). In the user's words: "a simple quadrature is verified by a whole battery of tests that shows a bunch of characteristics, including something elaborate like "invariant over reflection", which is a claim that is required when finally Sn hits a reflective boundary." At the top of the cone the claims compose into "SN behaves in XYZ way, with this kind of sensitivity/uncertainty, from numerical and method sources, bounded and measured. The method is verified. Verification Certificate".
9. **A catalogue of claim categories** (the user, 2026-10-16), with three readers:
   - the test-architect, as a guide to which claims a concept could carry ("not every concept needs every claim in the catalogue... it's just a guide of what one could claim");
   - qa, as a checklist for attack ("I managed to break X because claim A is important and you did not design a test for this claim");
   - the project, which learns from it what a robust scientific code must claim, as the catalogue grows.
10. **The catalogue is a guide, never an oracle, and it is organised by category, never by field** (the user, 2026-10-16). The test-architect and qa actively pursue the important claims that are not vacuous (invariances, conjectures, theorems) and look for the ones that should be proven and are not. A catalogue entry says "find invariances", never "find this specific invariance": a specific invariance is an instance of the invariance category and belongs to the concept's own claims.
11. **Each test says which category it is** (the user, 2026-10-16): "Tests should probably get a mark or something clearly say what they are according to the catalogue." The spelling (a marker naming the kind, beside the claim the test passes or fails) is open, and follows from where a claim lives (open question 5).
12. **A refusal is an elegance question, and the graph says where its debt clears** (the user, 2026-10-16). A guard that refuses illegal input is debt or not depending on where it sits in the graph: the graph shows the root, the producer upstream where the illegal state could be made unspellable instead (`coding-elegance` Pattern 4; the guard-as-debt rule of `coding-standards`). So the graph also serves the elegance-enforcer.
13. **How bounds compose up the cone is deferred until after the restructuring** (the user, 2026-10-16). No claim on it exists today, so the refinement does not study it. It is filed as an issue for after the restructuring. It is far from trivial: it is the uncertainty contributed by the MODEL itself, which modellers rarely measure; they usually bound the final result, and the usual scope of uncertainty quantification is the uncertainty brought in by the model's parameters. Measuring it would make ORPHEUS an outlier in computational modelling.
14. **The three-pass verdict method is accepted provisionally** (the user, 2026-10-16): "I accept this method for census verification, at least provisionally. We will revisit it during a second pass through the entire plan (where we will refine everything) when we go back to it when we tackle that campaign after finishing the one we're in." The campaign starts after the characteristic-reference campaign closes (`characteristic_reference_architecture.md`), with a second pass that refines this whole plan first.

15. **Invalidating a region of the DAG on a file change is foundational, and it rests on code changes only** (the user, 2026-10-17, answering the memo's hit-cost question): "You can also consider how we store a graph of tests, or Nexus graph, and consider how to validate rebuilds very fast. This is something that will become foundational to invalidating regions of the DAG based on changes in files. We might also want to pin revalidation strictly to code changes (so the files probably need to be stripped of documentation during pinning, maybe)." The section "One change-detection layer" below holds the hypotheses.

16. **The assembly file's slow block theorems go to ruling 15** (the user, 2026-10-18, on #587). `test_characteristic_assembly.py`'s 21 `slow` rows build the transport block in process because the block is their system under test, so no stored reference reaches them. Their cost is taken by verdict revalidation, not by a cheaper line rule now. They are ruling 15's second client, after C6.
17. **The campaign is active** (the user, 2026-10-18: "Seems like it's time to start the verification tree campaign"). The characteristic-reference campaign closed the same day. Per ruling 14 the work starts with a second pass that refines this whole plan.

## Vocabulary (one meaning per word)

- **Claim:** a statement the code is held to, with a home where it is written: an equation's label on a theory page, a law, a contract, an invariant. A claim is the node of the graph (ruling 4).
- **Test:** the instrument that passes or fails one claim. A test with no claim tests nothing nameable.
- **Edge:** one claim depends on another. Today it is spelled `@pytest.mark.rests_on(<test id>, ...)` on the upper test, naming a test, not a claim; the spelling has to follow ruling 4 (open question 5). It means three things: the upper test assumes the lower is green (dependency); a red lower test invalidates the upper reading rather than failing it (boundary); a failure propagates only along edges (causality). The definition is the test-architect's (`docs/development/agents/test-architect.md`, item 4).
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

- **Tests with a claim today** (`[M]` 2026-10-16, `main` `8a1512b1`; `pytest --collect-only` with a counting plugin, `scratch` `claimcount.py`; parametrised cases each count, so these count test ids, not functions): 2 837 of 17 764 collected ids carry a `verifies` marker (16 %); 14 927 carry none. By V&V level marker:

  | Level marker | With `verifies` | Without |
  |---|---|---|
  | foundation | 289 | 12 745 |
  | L0 | 1 235 | 635 |
  | L1 | 1 259 | 1 484 |
  | L2 | 54 | 31 |
  | none | 0 | 32 |

  The packages with the fewest claims, by ids without `verifies` against all ids:

  | Package | Without | All |
  |---|---|---|
  | `numerics` | 4 383 | 4 467 |
  | `mesh` | 564 | 564 |
  | `reference` | 387 | 387 |
  | `geometry` | 1 411 | 1 654 |
  | `transport` | 1 097 | 1 131 |
  | `sn` | 3 225 | 3 983 |
  | `derivations` | 2 206 | 3 453 |

  Under ruling 4, each of the 14 927 ids either names its claim in the census or is vacuous. Most foundation tests do assert something definite (a contract, an invariant, a law of a type); their claim is written only in the test's name and docstring, with no home a `verifies` marker could name. `[R]`, from reading the foundation tests touched this session, not a census.

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
- **Where a slow run's time goes, today** `[M]` 2026-10-17, the 24 reference files re-timed at `a456a9bf` (`scratch/characteristic_architecture/p3/close/timing/results.tsv`; serial, `-O`, `OMP_NUM_THREADS=1`): warm 3 073 s, `not slow` 487 s. 2 493 s of the warm total is the 21 slow rows of `test_characteristic_assembly.py` (ruling 16); the 40 slow rows of the other 10 characteristic files take 27.1 s warm, 20.7 s of it C6. The profile below is the state before the stored references, kept as the evidence that sent oracle values to the store.
- **Where a slow run's time went, 2026-10-16:** one population, the 40 slow rows of 10 characteristic files, a warm store, `-O`, serial; `scratch/characteristic_architecture/p3/c3/drop/timing/profile/`.
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

**A. The node is a test function.** This is the current spelling of `rests_on`. **[REFUTED 2026-10-16 for the question "what is a node", by the user's ruling 4]**: a test is the instrument of a claim, not the thing depended on. The FACT it establishes stays: the 603 existing edges are spelled as test ids and must migrate (open question 5).
- *Makes unspellable:* an edge to a test that does not exist; the edge resolves to a collected id.
- *Leaves welded:* a test that verifies two labels is one node, so an edge cannot say which of its claims the upper test assumes. The bypass drop's re-drop found this case: A1 carried `characteristic-door-reciprocity` (red under its mutation) and `characteristic-adjoint` (green under every mutation). It was one node holding one real claim and one empty one.

**B. The node is a claim: a test paired with one label.** **RULED 2026-10-16 (ruling 4).**
- *Makes unspellable:* an edge resting on a claim the lower test does not make.
- *Leaves welded:* nothing found yet.
- *Costs:* `rests_on` would name a claim, not a test id. A second spelling of the edge would then exist next to the 603 current ones, unless they migrate.

**C. Oracle values are nodes of their own, beside the tests** (concept 2). Under ruling 4 this is not a separate kind of node `[R]`: an oracle value is a claim ("this number is the solution of problem X to precision p"). Its test is the oracle's own validation, the slow twin that recomputes a frozen literal. The stored value is the evidence that claim produced, and a claim that rests on the oracle reads that value.
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
4. **Oracle nodes as stored references** (candidate C). The prerequisite is a test of whether the memo can generate from a function under `tests/`. [REMEDIED 2026-10-17 for the characteristic oracles: the mpmath routes in `tests/gates/derivations/_characteristic_point_mp.py` are traced-memo clients that generate in a child process, and their independence is the census gate `test_each_stored_reference_imports_nothing_but_the_memo` (`characteristic_reference_architecture.md`, concepts 1 to 3). The phase now means the other families' oracles, if the census finds any recomputed.]
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

## The claim catalogue (ruling 9)

**What it is.** A catalogue of the *kinds* of claim a concept can be held to, each kind with:
- its definition;
- the functional a test of it measures;
- the defects it is blind to: its stabiliser, in `instrument-doctrine` X1's sense;
- its known instances in the tree.

It is the dual of the error catalogue, `docs/theory/verification/error_catalog.rst`. The error catalogue records what has gone wrong, with the test that catches each defect. The claim catalogue records what can be asserted, with the tests that assert it. #358 already asks for "a claim-KIND axis"; this is that axis, given its own home.

**Seed kinds, from the user's list and the tree** (each `[HYPOTHESIS]` until the census confirms or merges it):

| Kind | What it claims | An instance in the tree |
|---|---|---|
| invariance under a group action | the object is unchanged by every element of a group (a reflection, a rotation, a permutation of equivalent regions) | a quadrature invariant under the reflections of a reflective boundary (`numerics/invariance.py`, `is_invariant_under`) |
| equivariance | the map commutes with the group action | the sweep commuting with a reflection of the mesh `[R]` |
| conservation, balance | a balance closes: production = loss + leakage, current continuity at a face | the characteristic block's conservation rows |
| exact identity | an algebraic identity holds exactly or to rounding: `R∘E = id`, an adjoint pairing, a Galerkin relation | the door's reciprocity `⟨r, φ_q⟩ = ⟨q, φ†_r⟩` (A1) |
| exactness degree | a rule integrates every polynomial up to degree d exactly, and fails at d + 1 | Gauss–Legendre's degree 2n − 1 |
| convergence order | the error falls as O(h^p) with the observed p | MMS ladders |
| bit identity | two computations agree to the last bit | a re-spelling that keeps the floating-point reduction tree |
| bounded drift from re-association | two orderings agree within a bound set by the reduction depth or iterations × condition × ULP, and are expected to differ | `vv-principles`, "Bit-identity vs principled-equivalence" |
| positivity, the cone | a flux or an albedo lies in its positive cone | ERR-043, ERR-046 |
| refusal | an illegal input is refused with a named reason, never answered. Whether the guard is debt is read from the graph (ruling 12) | `NoLeastSolution` on an (n,2n)-supercritical body |
| agreement with an independent reference | the answer matches a structurally independent reference within a measured band | the characteristic reading against the mpmath route |
| sensitivity and uncertainty | the answer's sensitivity to the discretisation and its numerical uncertainty are bounded and measured | the method-level certificate (ruling 8); no instance yet `[R]` |

**How it is used:**
- The test-architect's spec names, for each concept, which kinds it claims and why the others do not apply.
- qa attacks the kinds a concept should claim and does not.
- A census finding of a kind not in the catalogue adds a row.

**The certificate.** At the top of each method's sub-cone sits a composed claim. Its parts are the claims it rests on, each with its bound. Its statement is the method's verified behaviour, with the sensitivity and uncertainty from numerical and method sources, bounded and measured. How a bound composes up the cone (added, multiplied, taken as the maximum) is an open mathematical question, owed before the certificate can be more than a list `[R]`. It is deferred until after the restructuring (ruling 13).

**Where the catalogue lives** is open. A theory page under `docs/theory/verification/`, beside the error catalogue, with each kind a directive the graph indexes, is the obvious candidate `[HYPOTHESIS]`.

**The worked example.**
- `[M]` 2026-10-16 (code reading: an import plus a call): `numerics/invariance.py` imports `RigidMotion` (line 56). `is_invariant_under` decides a measure's invariance by applying the group's motions to the nodes (`_orbit_space_closure`). The specular boundary reads `RigidMotion.preserves` (`geometry/boundary/_specular.py:113`).
- So the reflection-invariance claim of a quadrature rests on the claims of `RigidMotion` (its action on a direction is an isometry, and a reflection is an involution) and of `SubgroupOfO3`. SN's reflective-boundary claim rests on it in turn.
- Whether those `RigidMotion` claims are tested today, and as which kinds, is a census question.

## How census verdicts are decided (proposed 2026-10-16; accepted provisionally, ruling 14)

Per layer, in three passes. Each pass narrows the next, so no mutation is run against a test it cannot reach.

1. **Capture what each test executes.** One coverage run of the layer's tests with per-test contexts (`dynamic_context = test_function`, `--show-contexts`), ingested into Nexus.
   - This is the instrument the plan already owes the `rests_on` falsifier: a claim resting on another must execute what the lower claim executes.
   - Captures of this kind exist for every package (`runtime_runs`: the "wide capture v2" runs). They predate the `tests/gates/` layout, so they are re-taken, not reused `[M]` 2026-10-16, the run list.
2. **A reading pass by the test-architect.** For each test it records:
   - the claim it passes or fails, and its kind from the catalogue;
   - its rung, and the claims it should rest on;
   - a candidate verdict: keep, strengthen, re-base, merge into a stronger twin, or vacuous.

   The candidate verdicts are flags, not rulings. Reading finds the claim, which no instrument can. It cannot prove a test blind or dominated, which only a mutation can.
3. **Adjudication by a kill matrix.** Rows are mutations, columns are tests; a cell is red or green.
   - **Designed mutations:** at least one per claim, inside the claim's algebraic class (X1), aimed at what the claim says.
   - **Generated mutations:** cosmic-ray over the layer's modules, as a sweep for defects nobody named. `[M]` the harness's founding run: 373 of 374 mutants killed on `diamond.py`, about 1 s per mutant (`tests/gates/_mutation/README.md`).
   - **Each mutation runs only on the tests whose capture executed the mutated lines** (pass 1). The matrix's cost is then a few tests per mutation, not the whole layer per mutation.
   - Mutations are textual, on a copy of the tree. Memo-served values cannot be reached in-process (`vv-principles` #17), as the bypass drop's re-drop showed.

**The verdicts the matrix decides, and how:**
- **Blind:** a test reddened by no mutation of its claim. Before it retires, qa gets one attempt to find a mutation that reddens it. A test qa can redden was under-attacked, not blind.
- **Dominated:** a test whose red set is contained in another test's red set, for the same claim. The containment is relative to the mutation set: a test can look dominated only because the set lacks the one mutation it alone catches. So a dominated verdict is confirmed by reading, and qa gets one attempt to find a mutation that separates the two, before the weaker test retires into its stronger twin (ruling 3).
- **Strengthen:** a claim the catalogue marks important (an invariance, a theorem) whose own designed mutation no test reddens. The claim is real and its instrument is missing.
- **Re-base:** a test that reddens under mutations of a simpler claim it does not declare resting on. The matrix shows the dependence; the edge is then designed and declared.

Mutant subsumption, the containment of one test's kill set in another's, is the standard instrument of mutation testing for redundant tests. Here it is run per claim, not per suite.

## Open questions for the user

Questions 1–3 were answered on 2026-10-16 (rulings 4–6).

4. **Census verdicts:** accepted provisionally (ruling 14).

5. **Where a claim lives that has no equation** (follows from ruling 4). [Answered in part, ruling 7: both homes, chosen per claim; the rule for choosing is open.] An equation's claim lives at its `:label:` on a theory page, and `verifies` names it. A contract, an invariant or a law of a type has no such home today. 12 745 foundation test ids carry no `verifies`. The candidates, each `[HYPOTHESIS]`:
   - a claim directive on the theory or architecture page that owns the concept, minted like `.. error-entry::` (a node in the graph, an id a marker can name);
   - the claim written once in the test and harvested by the graph.

   The first keeps the corpus as the knowledge base (Cardinal Rule 3); the second keeps the claim next to its instrument. How `rests_on` names a claim follows from this answer.

## One change-detection layer for the memo, the test DAG and the graph (ruling 15)

Three consumers ask one question: *which of my nodes did this edit to the tree invalidate?*
- **The traced memo** asks it per entry, from the entry's manifest (the defs that ran, the module skeletons, the distributions, the files read).
- **The test DAG** would ask it per claim: a changed node invalidates its cone, and only that cone re-runs (phase 5; #358 for a red support).
- **Nexus** asks it per build: today `sphinx-build` rebuilds the whole graph `[R]`, not measured here.

**What exists** `[M]`:
- The memo's pins already ignore documentation. A def is pinned by the digest of its AST with docstrings, comments and layout removed, and a module by its skeleton (every body replaced by `pass`, docstrings stripped). The gate is `test_traced_memo_validation.py::test_m2_1_…` (rows `docstring`, `comment-added`).
- Since `53551868` (#602), the memo keeps what validations read for the whole process, keyed by each file's identity on disk (device, inode, size, modification and change times), with git's guard against a write too recent to trust. A warm hit no longer re-reads an unchanged file.
- The parsed index of a source text is stored on disk by the digest of the text's bytes, the interpreter and the memo's own source (`_source_index`), so the parse is paid once per content.

**Hypotheses** `[HYPOTHESIS]`, for the second pass over this plan:
1. **Two stages, as git and build systems do.** First the stat identity says which files MAY have changed; then the AST digest of each such file says which defs and skeletons DID. A doc-only edit passes stage one and stops at stage two. The memo already has both stages, so this would be its machinery shared, not built a second time (Cardinal Rule 2).
2. **The node a pin names is a def, not a file.** An edit to one function invalidates the entries and claims that ran that function, not the whole module. The memo's `DefPin` already has this granularity; Nexus's function nodes and the test DAG's claims could share the same key: the qualname plus the digest of the normalised AST.
3. **The edges come from execution where they can.** The static call cone has 12 to 15 % recall against execution evidence (`retest`'s own note, ORPHEUS 2026-08-18). A memo manifest is an execution record of one call. A coverage capture with per-test contexts is the same for a test. The invalidation cone of a test would be read from what it ran, not from the call graph.
4. **One store of digests per checkout.** The memo's index store, the graph database and a test-result store would key one table: file identity, then def digest. A rebuild would read only the changed rows.

Each hypothesis is measured before it is built: the share of edits that are doc-only, the cost of stage two on a full tree, and whether Nexus can rebuild incrementally at all (a question for sphinxcontrib-nexus).

**The first concrete client: C6** (the user's ruling, 2026-10-17, accepting the recommendation to leave C6 in process and send its cost here).
- **The test:** `test_the_volume_integral_of_the_reading_is_the_galerkin_block` in `tests/gates/derivations/test_characteristic_reading.py`, verifying `characteristic-galerkin-assembly`. It is a theorem row, and both sides are the characteristic reference's.
  - One side is the Galerkin block `K_ij`, an integral over the lines of space.
  - The other is the point reading of each basis function `u_j`, read at the 648 points of a volume rule and integrated against `u_i`.
  - The two agree to 1e-11 of max|K|.
- **The cost:** its two slow white fixtures cost about 10 s each warm `[M]` 2026-10-17 (`p3/c3/concept3/slow10_warm.log`). That is 20.7 s of the 27.1 s warm slow run over the 10 characteristic files. The time is the 648 calls of `PointRule.of(...).row(transport)` in `_rows`; C8, the reciprocity row, reuses them.
- **Why a stored value is the wrong fix:** storing the rows (two keys) would make a warm run compare the live block with a stored copy of the code's own earlier output. The reading code would not run, so neither a mutation nor a coverage capture would see it. Nothing independent would remain to carry the check.
- **Why a verdict is the right one:** C6's pass or fail is a function of the code it ran and nothing else, so it is the ideal case for a verdict keyed on its execution cone (hypothesis 3). It re-runs only when a def in that cone changes (hypotheses 1 and 2), and the stored claim reads "verified against this code".
- **The pins already strip prose** `[M]` 2026-10-17: a docstring-only edit to a stored reference's function or module leaves its entries warm (`p3/c3/concept3/probe_docedit.log`, `probe_moddoc.log`). Ruling 15's "stripped of documentation during pinning" is therefore already the memo's behaviour at def and module granularity (gate M2.1). The verdict store would inherit it, not build it.
- **Its place in the certificate (ruling 8):** C6 is one constituent claim of the characteristic method's certificate. The certificate says what C6 contributes. It does not make C6 cheaper; the revalidation does.

## Side findings to place

- **The Nyström row's marker** (`peierls-greens-slab-architecture`) runs only under `ORPHEUS_RUN_WITHDRAWN=506`. In the canonical run it covers nothing. Its margin is 7.3e-5 against a 5e-5 band, for an attenuation error of 1e-3. Record it under #506. (`scratch/characteristic_architecture/p3/c3/drop/redrop/verifies_verdicts.md`)

## The second pass (opened 2026-10-18)

Per ruling 14, the whole plan is refined before phase 1. The pass reconciles the plan against the tree (`plan-authoring` §7), revisits the provisional verdict method, settles open question 5 if evidence can, and fixes the phase order.

**Reconciled so far** (2026-10-18, the orchestrator): ruling 6 landed; phase 4's prerequisite is met for the characteristic oracles; the slow-run measurement is replaced by the 2026-10-17 re-timing; ruling 15 has two clients (C6 and the assembly slow tier, ruling 16). The counts under "Measured state" are from 2026-10-16 and are re-taken at phase 1, not here.

**First-pass reviews dispatched** (W5, adversarial; reports under `scratch/verification_tree/pass2/`):
- the test-architect, who will own the census: open question 5 by a discriminating sample; where the three-pass verdict method is blind and what coverage per test costs on `sn`; whether ruling 15's verdict store precedes the census; stale or mis-predicated measured state; what is missing;
- the cross-domain attacker: foreign frames that already solved these structures (content-addressed build systems and early cutoff, proof assistants' dependency tracking, mutant subsumption, assurance cases).

⏸ **Next:** read both reports with full context, write their findings into this plan in place (each with its marker), and bring the user the questions the pass must rule: the census method (ruling 14's revisit), open question 5, and the phase order. #604 is the deferred bound composition (ruling 13).
