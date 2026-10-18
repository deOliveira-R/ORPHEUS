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

18. **Execution is delegated** (the user, 2026-10-18): "my intention is that plan execution will either be carried out by a fan out of Sonnet models or an ultracode session. So we should basically reach a state that is good enough to delegate sections of the execution." The polish criterion for this plan is therefore delegability: each phase is cut into work packages that an agent with no session context executes from the text alone. Each package states its input set, its fixed output schema, the evidence every field must cite, and an acceptance check the orchestrator runs mechanically on the return. What needs the user's or the orchestrator's judgment (a ruling, an ontology question, a retirement) is named as such and kept out of the packages.

19. **A claim is a statement, and the labelling is open to redesign** (the user, 2026-10-18, on S1): "Option 1 could be good, but don't feel like you need to fit to the current labels. If the labels and the way it is defined is not appropriate anymore and a different labelling methodology would be better, suggest that instead." The statement-with-mentions shape is accepted as the direction. The `verifies` labels of today are not a constraint. Candidate D below is the proposal.
20. **Three edge kinds: executes, reads, supports** (the user, 2026-10-18, on S4 and S5). *Executes* and *reads* are the two dependences of S4, spelled with the memo's pin types. *Supports* is the certificate's inference, with argument nodes; ruling 13's deferral becomes one argument node marked assumed.
21. **The census adopts the repaired verdict method** (the user, 2026-10-18, on S7 to S10; this is ruling 14's revisit). A duplicate candidate needs EQUAL kill sets plus the same claim, rung, warrant and role. Strict containment proposes an edge, never a retirement. Dominance is stated against a named set of dominator mutants. Full red sets are computed only inside a claim group. Generated mutants run with first-kill stop. Designed mutations are the test-architect's. The census is proven first on planted known answers. Each test carries a role and a warrant (S2, S3).
22. **The cone store comes first, and the cone is built as the tree is organised** (the user, 2026-10-18, on the phase order): "let's build the cone dynamically as we change, so that we can use the cone to guide the population. This means we're doing something procedurally similar to a sort, where as the sort pass progresses, the sequence gets sorted. In this case, the cone is formed as it is organized, to guide further organization." Read by the orchestrator `[R]`, to be confirmed:
   - The census runs bottom-up by layer (ruling 5), like an insertion sort whose sorted prefix is the layers already organised. Each layer's run records the execution cones of the tests it runs, and its claims become placed nodes.
   - A test of the next layer is placed against the placed nodes: a placed claim whose cone the test's cone contains is a proposed `executes` support (S7's containment signal), and a stored value it read is a proposed `reads` support. The reading pass and the mutant selection for that layer are then guided by the partial DAG, not by the whole suite.
   - The cone store is incremental by construction: what the census captures is the store, and an edit re-keys only the rows whose defs changed (ruling 15).
   - So there is no separate up-front build of the whole store before the census. The store's machinery (keys, the union of S8's four sources, its positive controls) is built first; its contents grow with the census.

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

  [REFUTED 2026-10-18 as a measure of the census's load, the test-architect (`scratch/verification_tree/pass2/test_architect_attack.md` F6): a foundation test carries no `verifies` by rule (`vv-principles`, the level taxonomy), so the column measures the rule, not a gap; the unit of a claim is the test function, `[M]` 7 686 functions against 17 764 ids (AST, `pass2/ta/ast_census.py`); and `claimcount.py` no longer exists, so the 2 837 / 17 764 split cannot be re-run. The anomaly worth listing is the reverse cell: 153 functions carry both `foundation` and `verifies`.] Under ruling 4, each of the 14 927 ids either names its claim in the census or is vacuous. Most foundation tests do assert something definite (a contract, an invariant, a law of a type); their claim is written only in the test's name and docstring, with no home a `verifies` marker could name. `[R]`, from reading the foundation tests touched this session, not a census.

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

**D. A claim is a named statement minted at its concept's home; tests instrument it; equations are what it mentions** (proposed 2026-10-18 by the orchestrator, for ruling 19; `[HYPOTHESIS]` throughout, nothing measured).
- *The object.* A claim has an id, a statement (one sentence in the domain's words), a kind (the catalogue's), a warrant (theorem, reference with its pillar, record), a domain (where it holds), a band (exact, or a tolerance with its measured floor), the things it mentions (equations, types, functions), and its edges (`executes`, `reads`, `supports`) to other claims. Its **statement digest** covers the statement's fields (statement, kind, warrant, domain, band, mentions), normalised like the memo's pins, so a rewording of prose around it does not move it.
- *Where it is declared.* One vocabulary, two declaration sites, one namespace (ruling 7 kept, its rule for choosing made concrete):
  - a mathematical claim, an equation's or a theorem's, by a directive on the theory page that derives it, minted like `.. error-entry::` (a graph node, an id a marker names);
  - a law of a type (`vv-testing`, "A type that embodies a mathematical concept ships the test of its defining laws") at the type itself, declared in Python beside the class and rendered into the corpus by the same directive reading it, as autodoc does;
  - a structural claim about the code (S3: an import graph, a single caller) on the architecture page that states the structure.
- *The id.* Namespaced by the concept that owns it, `<concept>.<claim>` (for example `gauss-legendre.exactness-degree`), so the owner of a claim's name is the owner of the concept, and one namespace holds every id (TA F8.6: two names for one claim are a duplicate dominance cannot see). Nexus refuses a duplicate id as it refuses a duplicate ERR id.
- *The test side.* A test names the claim it instruments, its role and, for a parametrised test, the instance its rows range over: `@pytest.mark.claim("<id>", role=..., instance=...)`. The edges live on the claim, not on the test, because dependence is between claims (ruling 4). `verifies` is retired into `mentions`: an equation is verified when some claim that mentions it is passed.
- *The band at the claim.* The test reads its tolerance from the claim, not from a literal of its own. Then a tolerance cannot be widened without editing the claim's statement, which moves its digest and raises "assumption changed" on every claim above it. The standing rule "never widen a tolerance without a red test and the user's sign-off" becomes structural.
- *Makes unspellable:* a role-`claim` test with no claim (a collection gate); an edge to a test rather than a claim; one test credited with three claims it mentions; a band widened in a test with nothing above noticing; two ids for one claim.
- *Leaves open, or costs:*
  - The tests must read declared claims at collection without depending on a docs build, so the declaration's source has to be importable Python, or a tracked file generated from the corpus and checked for drift (`plan-authoring` §10: a gate's inputs are part of the commit). Which one is open.
  - Migration: each existing `verifies` label is a mention of a claim nobody has written. The claims are written by the census's reading pass, never minted from the labels by a pattern (ruling 2).
  - The band at the claim needs, for each claim, the tolerance its tests use today; where two tests of one claim use different bands, they are two claims or one is mis-banded, a census finding.
  - Heterogeneous rows (S1) are split into one claim per row, or given a per-row claim marker.

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

Mutant subsumption, the containment of one test's kill set in another's, is the standard instrument of mutation testing for redundant tests. Here it is run per claim, not per suite. [REFUTED 2026-10-18, the cross-domain attacker: mutant subsumption is the containment of one MUTANT's killing-test set in another's (Ammann, Delamaro and Offutt, ICST 2014); containment of TESTS by kill set is test-suite reduction. They are the two readings of one matrix. And containment of one test's kill set in another's is what a correct `rests_on` edge looks like, so as written this rule retires the foundation (second-pass finding S3).]

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

Both reviewers were told of ruling 18 after dispatch and asked to add: the cut into packages, the census's per-test output schema, where a Sonnet-class agent would go wrong and the check that catches it, and each package's acceptance gate.

### The reviews' findings (read 2026-10-18 by the orchestrator, with session context)

Reports: `scratch/verification_tree/pass2/test_architect_attack.md` (TA below; probes in `pass2/ta/`) and `pass2/cross_domain_attack.md` (CD below; census scripts `census.py`, `c2.py`). Each finding is labelled S1 to S12 for this pass only. Each proposed change is `[HYPOTHESIS]` until the user rules.

**On the node (ruling 4).**
- **S1. A test does not pass or fail one claim in 108 of 553 measured cases.** `[M]` CD's AST pass: 553 test functions carry a function-level `verifies`, and 108 name two or more labels; 39 of those are parametrised. They come in two shapes. A *relational* claim is one statement that mentions several labels (`test_kinf_homogeneous`: k is the eigenvalue of the pencil, naming the fission, eigenvalue and removal labels). *Heterogeneous rows* are one function whose rows are different claims (`test_the_closed_angular_flux_is_the_unfolded_backward_path`, a row per body).
  - Proposed (CD, from proof assistants): a claim is a **statement** (its predicate, its domain, its band) with its own digest. The labels are what the statement *mentions*, many per claim. A parametrised function is one universal claim and its rows are instances; an edge may rest on one instance. Heterogeneous rows are several claims sharing one body, and need a per-row label or a split.
  - An edge records the lower claim's statement digest. When a lower statement changes (a band loosened, a domain shrunk), every edge into it enters a new state, **assumption changed**, which a re-run cannot clear; a reader decides whether the upper claim still follows.
- **S2. A claim needs several tests.** `[M]` TA's reading, one reader: 5 of 34 functions in `numerics/test_invariance.py` are controls or declared-blindness rows, and 2 of 26 in `sn/eigenvalue/test_keff_curvilinear.py` are premise rows. Doctrine requires them (`vv-principles` anti-patterns 11, 19, 37). Under the plan's verdicts a control reads blind and a premise row reads dominated, so both would retire. Proposed: a **role** per test (`claim`, `positive_leg`, `negative_leg`, `control`, `premise`, `declared_blindness`, with `leg_of`); blind and dominated apply only to role `claim`.
- **S3. Two kind axes.** The catalogue's kind (what is claimed) is orthogonal to the **warrant** (why it is believed): THEOREM, REFERENCE (with its pillar), RECORD (a frozen table, a digest, a snapshot). TA also names missing kinds: structural claims about the code (`[M]` 13 of 34 functions in `test_invariance.py`: an import graph, "exactly one caller", retired names), termination, the qualitative direction of an effect, parameter contracts.

**On the edge.**
- **S4. `rests_on` welds two dependences** (CD, from build systems). *Executes*: the upper test runs the lower's code; invalidation flows through code, and the falsifier is cone containment. *Reads*: the upper reads a stored value the lower validated (every oracle edge of candidate C); invalidation flows through the payload digest, and the cone falsifier is FALSE for it, since a memo hit runs none of the oracle's code. The memo already has both spellings: `DefPin`/`ModulePin` is an executes cone, `ChildPin` is a reads edge. Proposed: two edge kinds, reusing the memo's pin types (Cardinal Rule 2).
- **S5. Dependency is not inference** (CD, from assurance cases, GSN and Assurance 2.0). `rests_on` serves diagnosis and scheduling: "if B is false, A's reading is unexplained". The certificate needs inference: "B1 to Bn, by an argument, entail A", and the two relations differ in both directions. Proposed: an **argument** node with a status (stated, proven, assumed) and a `supports` edge, separate from `rests_on`. Ruling 13's deferred bound composition becomes one argument node marked assumed. A query `axioms(certificate)`, after Lean's `#print axioms`, lists the leaves no verdict adjudicates: owed ones (a withdrawn or skipped instrument, like the Nyström row in "Side findings") and declared trusted ones (mpmath's quadrature, BLAS, IEEE arithmetic, a literature constant).
- **S6. Evidence and confidence are separate, and go stale separately** (CD). A pass is evidence. The kill-matrix row that shows the instrument is not blind is the confidence in that evidence. A verdict store that keys "passed" and not "non-vacuous at this code" certifies with evidence whose discriminating power was measured against older code. Proposed: two digests per accepted claim, and a state **passed, confidence stale**.

**On the verdict method (ruling 14's revisit).**
- **S7. Containment of kill sets marks an edge, not a duplicate** (TA F1, `[R]` from the plan's own falsifier). An upper test runs everything its support runs, so its kill set contains the support's; the plan's dominance rule then retires the support. RECORD rows redden under nearly any mutation, so they would dominate the theorems beside them (`[M]` two such rows in `test_invariance.py`). Proposed: a **duplicate candidate** needs EQUAL kill sets over a set that includes the claim's designed mutations, plus the same claim, rung, warrant and role. Strict containment proposes a `rests_on` edge, designed by reading, never a retirement. Dominance is stated against a named mutant set and computed over its **dominator mutants** (CD), since subsumed mutants inflate kill counts (Papadakis et al., ISSTA 2016).
- **S8. Per-test coverage cannot see four populations** (TA F2, CD).
  - `[M]` TA, one file (`test_keff_curvilinear.py`, coverage with `test_function` contexts): 8 967 of 14 720 executed lines carry only the import context, and 3 093 of those are executable statements (constants, dataclass defaults, registries, decorators). A mutant there selects zero tests and reads as survived.
  - Child processes: memo generation, and gates that spawn subprocesses. Coverage sees them only under `COVERAGE_PROCESS_START`.
  - Memo-served values: on a warm store the oracle's code runs in no test. A textual mutant inside an entry's manifest forces a regeneration.
  - In-process caches: `[M]` CD, 220 cache lines in 51 files of `orpheus/` and 53 wide-scope fixtures. `[HYPOTHESIS]` the second test to touch a cache records no cone for the cached code.
  - So the execution cone the census needs is a union: coverage contexts, memo manifests, import-time attribution (to the tests that read the bound name) and child-process capture. That union is the same object ruling 15's hypothesis 3 wants. Hypothesis 3 is [REFUTED 2026-10-18 as stated] where it calls a per-test coverage capture "the same" record as a memo manifest: the memo records in a fresh process precisely because served-from-memory values hide dependencies.
- **S9. Cost** (TA F5). `[M]` coverage with per-test contexts costs 1.25 to 1.67 times plain on two `sn` files; `[R]` about 25 min for `sn` and 95 min for the suite. Its JSON export does not scale: `[M]` 660 MB for one 35-test file against 6.2 MB of SQLite data, so Nexus must read the SQLite data or take per-file exports (a sphinxcontrib-nexus question). A full kill matrix over `sn` is `[R]` about 13 000 mutants and 36 to 145 serial hours. Proposed: full red sets only inside a claim group; generated mutants with first-kill stop, as the sweep for unnamed defects; designed mutants, about one per claim, written by the test-architect.
- **S10. The census has no positive control** (TA F8). Before a layer's verdicts are believed, the census must classify planted known answers correctly. One is found: `[M]` 2026-10-18, by reading `orpheus/sn/solver.py:1690`, on an all-reflective problem without (n,2n) SN's `compute_keff` IS production over absorption, so `test_keff_curvilinear.py::test_particle_balance` and its three siblings compare a definition with itself (X4). Further cell rules are owed for a strict xfail (green under any breaking mutant), `slow` and withdrawn rows (deselected), stochastic rows (a seed protocol), `pytest.raises` without `match=` (green under a mutant raising earlier with the same type; `[M]` CD: 924 functions check only `pytest.raises`), and equivalent mutants (relative to the claim's band).

**On ruling 15.**
- **S11. Where early cutoff acts** (CD). On a verdict bit across a `rests_on` edge it buys nothing, since the upper cone contains the lower. At the normalised-AST digest it exists (gate M2.1). At an oracle's value it was designed and is not taken: #610, filed 2026-10-18 (`ChildPin.stale` marks the parent stale without regenerating the child and comparing its payload). At a claim's statement it is missing (S1).
- **S12. "Code changes only" cannot mean "not the environment"** (CD). Bit-identity and bounded-drift claims depend on the reduction order, hence on the thread count, and phase 5 would give each test its own thread count. The verdict key includes the environment the memo already pins (`EnvironmentPin`), or those claim kinds are run at a fixed thread count whatever the scheduler chooses. Read in the domain, ruling 15's "code changes only" means "not documentation".

**On the phase order** (TA F4, `[R]`). Census pass 1, pass 3's mutant selection, the `rests_on` falsifier and ruling 15's invalidation all consume one object: the set of defs each test executes, keyed by `DefPin`. So the proposed order is:
1. phase 1 (measure) plus the **cone store** (S8's union), with its own positive controls (a mutant of an import-time constant selects the tests that read it; a mutant of oracle code selects the oracle's validation test);
2. the reading pass, in parallel with step 1, since it needs neither;
3. the verdict cache on top of the cone store, trusted only after a validation run in which N sampled mutants run in full and every cached verdict agrees;
4. the kill matrix.

**On open question 5** (TA F7). A structural argument decides half of it `[R]`: ruling 1 merges and retires tests, so a claim that another claim rests on needs an id independent of its test. A claim harvested from a test is honest only with no dependants. TA adds a third home: a directive in the type's own rendered docstring (trap: 24 of 48 `automodule` directives carry `:noindex:`, which mints no target). It proposes a four-file sample over two layers (`numerics/test_invariance.py`, `numerics/test_measure.py`, `sn/eigenvalue/test_keff_curvilinear.py`, `sn/operators/test_streaming_collision_operator.py`), recording per function its kind, where its claim is already stated, its fan-in, its fan-out and its role, with a decision rule in the report.

**On delegation (ruling 18).** Both reports agree on the acceptance principle: a small trusted checker accepts each delegated output, and the producer is never trusted (the de Bruijn criterion). Every output kind returns a witness the orchestrator re-runs or checks by containment (CD's table). TA gives the census's per-function JSON schema, the likely errors of a Sonnet-class reader with the check for each, and an eight-point acceptance gate (TA F9). The one output with no mechanical check is "this is the claim the test makes", which goes to a second reader. Delegable: phase-1 measurement, the per-layer capture, the reading pass per file, generated-mutant sweeps, applying ruled verdicts. Not delegable: claim grouping across files, new catalogue kinds, the claim namespace, edges across packages, final retirements, the cone store's design, designed mutations.

**Measured state, corrections not yet applied above** (TA F6, `[M]` 2026-10-18): 636 test files (was 625); 609 lines naming `rests_on` in 95 files (a line count, not the earlier mark count, so no delta is claimed); 12 of the 232 `sn` files carry one; H-SEQ's premise rests on the old mpmath-bound population, while today's slow population is numpy-bound (ruling 16), so its S(p) is not known to be flat.

⏸ **Next:** rulings 19 to 22 recorded (2026-10-18). Candidate D is with the user for discussion, and ruling 22's reading for confirmation. Then: the four-file sample (TA F7), now run against candidate D (does each function's claim fit D's object, and which declaration site does it take?), and the cone store's specification by the test-architect, cut into packages per ruling 18. #604 is the deferred bound composition (ruling 13); #610 the memo's early cutoff.
