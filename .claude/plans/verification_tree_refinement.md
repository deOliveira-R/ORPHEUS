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
   - **Confirmed** (the user, 2026-10-18), with the reason: "otherwise, we would have to keep some sort of catalogue of sorts and organize everything later, which is fragile. It seems better to organize as we go already as a DAG that we can later refine."
23. **The DAG is identified by layer and by sub-cone, if measurement supports it** (the user, 2026-10-18): "It would be good also if the DAG was already identified by layers and by sub-cones (for example, at the transport method layer, there is an Sn sub-cone, a CP sub-cone, etc., in the numerics layer, there might be a discrete measure sub-cone, a space sub-cone, etc). We have to measure this and check if this is a good idea. After measurement, check if this influences the design of candidate D." The measurement is specified below under "Sub-cones: the measurement"; its result decides whether a sub-cone is an object of the DAG, and then what it changes in candidate D.

24. **A failure to find sub-cones is evaluated, by the mathematics or other criteria** (the user, 2026-10-18). It may be intrinsic, or it may show a concept leaking into another. The criteria are under "Sub-cones: the measurement", open to the user's refinement.

25. **Sub-cones form a DAG of their own, with joins** (the user, 2026-10-18): "The Frame is a sub-cone that is after the sub-cones of Discrete Measures and Basis functions. There is a contraction after the Discrete Measures and Basis functions subcones, where both sub-cones should feed into the Frame sub-cone." So a sub-cone is not a cell of a flat partition of a layer. It has upstream sub-cones, and a join (here Frame, a measure paired with a basis) has two or more. The sub-cone DAG is the quotient of the claim DAG by sub-cone membership. The measurement is re-defined accordingly (below).

26. **The sub-cone's grain is the concept, built from the claims' homes as the census places them** (the user, 2026-10-18: "I agree with your recommendation. We can build it this way first and refine later."). Meanwhile the module graph, with the typing-only and lazy imports counted, is the leak instrument (S18).
27. **The numerics knot is settled before the census** (the user, 2026-10-18: "settle it first"). It is its own living plan, `numerics_untangling.md`, whose concept order becomes the base of the cone.
28. **SN stays one sub-cone hanging off `SNProblem` for now, and the constructor rule gets its own campaign** (the user, 2026-10-18): "It's fine that the SN is essentially a big blob depending on SNProblem *for now*. This is a fairly strong indication that we have not been diligent enough with our rule of "__init__ constructor being the raw constructor for testing and class methods the specialized constructors that receive convenient arguments that needs unravelling", so we need to revit it in its own campaign. Untangling the numerics layer is a pretty good start." Filed as #611. The rule is not written in the rules corpus today (`[M]` a grep of `docs/development/rules/` and `.claude/rules/` for "constructor").

29. **The sub-cone DAG states the code's ontology** (the user, 2026-10-18, in `numerics_untangling.md`, ruling 3 there): the way the sub-cones connect says what machinery exists, what depends on what, at what level, and how it stacks. So the census's concept order and the code's concept order are one object, and a knot in either is a knot in both.
30. **The sub-cone measurement runs once the numerics knot is settled, as an ontology instrument; the code is then reorganised into the DAG** (the user, 2026-10-21): "Having the sub-cone measurement will be useful by itself, from an ontological perspective, because we will be able to look at the ontology and say 'why is this concept here and not upstream?' and things like this. So it will allow more effective targeted investigation. Our priority will be to reorganize the entire code into proper DAG and sub-cones so that we can look at a concept, look at the graph and say 'this should be there, not here'." The knot is settled at S7 of `numerics_untangling.md` (its ruling 89: 0 violations of the concept order over 112 module edges, measured on the S7 worktree).
31. **Layer by layer: untangle, organise, then the next layer, until the graph shows the ontology; a backlog runs beside it** (the user, 2026-10-22): "After re-organizing numerics, we will eliminate cycles on the next layer, and also organize that layer, and so on, until all layers are properly organized and the graph shows the ontology. We will keep a backlog of things we noticed, but did not finish in order to prioritize the DAG. We can work on them in parallel in a worktree as well." The order is the import layers (ruling 5: `numerics/`, then `geometry/`, `mesh/`, `data/`, then `transport/`, then each method sub-cone). For each layer: measure its cycles (the G0 builder, generalised to the layer), untangle them, measure its sub-cones, reorganise. The backlog is GitHub issues (Cardinal Rule 4), listed under "The DAG backlog" below; a backlog item is worked in its own worktree, in parallel with the DAG line, never on it.

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

### Sub-cones: the measurement (ruling 23; specified 2026-10-18, `[HYPOTHESIS]` until run)

**The question.** Within a layer, do the tests fall into sub-cones: groups whose tests execute mostly one region of the layer's code (plus lower layers), with few dependences between groups? If they do, the sub-cone is a real object and the DAG names it. If tests straddle regions, a sub-cone would be an imposed filing and would fragment claims.

**The candidate partitions of a layer's code**, compared against each other:
- (a) by module or sub-package (`orpheus/numerics/<module>`);
- (b) by concept, from `docs/architecture/conceptual_view.rst` (a discrete measure, a space, a frame, an operator...);
- (c) found from the data: communities of the bipartite graph of tests and the production defs each test executes, from the per-test coverage capture. A flat community method cannot represent a join (ruling 25), so the order is derived from containment: sub-cone A precedes B when B's tests execute A's code and A's tests do not execute B's.

**The readings**, per partition:
- *purity*: for each test, the share of its executed in-layer lines that falls in its dominant group; the distribution over the layer's tests; [REFUTED 2026-10-18 as the reading, by ruling 25: a flat share calls every test of a join impure. Replaced by *upstream-closed purity*: a test's **home** is the most downstream sub-cone whose code it executes, and the test is pure iff everything else it executes in the layer lies in its home's upstream closure.]
- *the quotient*: the sub-cone DAG (sub-cones as nodes, an edge where one sub-cone's tests execute another's code). It must be acyclic; a cycle means a wrong cut or a leak. Its joins are listed;
- *cross-group dependence*: how many test pairs would be proposed edges by containment (S7) across groups against within groups; and the existing `rests_on` edges, within against across;
- *agreement*: how far (a), (b) and (c) agree (a partition-agreement score, with the tests each disagrees on listed);
- *nesting*: whether the groups nest (a sub-cone containing smaller ones) or overlap.

**Positive controls** (X1): a test known to be a unit test of one concept must come out pure, and a known cross-cutting test (an end-to-end solve) must come out impure. If the controls do not separate, the instrument cannot decide the question.

**Population.** The `numerics` layer first, the base of the cone (ruling 5), with its 109 files. Then one method package, to test whether "one method package = one sub-cone" holds, and whether a method package itself splits (the SN sweep, its boundaries, its eigenvalue).

**A straddler is judged, not just counted** (ruling 24, the user, 2026-10-18): "note that not being able to find well-defined sub-cones might be instrinsic, or might be an indication of a concept leaking onto another. that must be evaluated based on the math involved (or some other criteria which I'm open to suggestions)." A test that straddles groups is classified INTRINSIC or LEAK. The orchestrator's proposed criteria, `[HYPOTHESIS]`, tested and refined by the measurement:
1. *Mention:* the claim, read from the test's assertions, needs both concepts because it is about a composite or a morphism between them (a frame is a measure paired with a basis). Then it is intrinsic.
2. *Interface:* the straddle goes through the other concept's named, typed interface (intrinsic or incidental), or through its internals, a private helper or an inline recomputation (a leak).
3. *Twin:* the other group's executed code is a second spelling of something the test's own group computes (Nexus `twin_paths`). Then it is a leak (Cardinal Rule 2).
4. *Direction:* the dependence follows the layer and import direction (a composition on a lower concept), or goes sideways or against it (a leak suspected; Nexus `native_place` on the executed defs).
Under ruling 25, foreign code in an UPSTREAM sub-cone is composition by definition and needs no criterion. Only foreign code in a sibling or downstream sub-cone goes through the five criteria.
5. *Counterfactual:* the claim can be posed with the other concept replaced by anything satisfying its interface. Then the straddle is incidental: the test is re-based onto a simpler node, and the code has no leak.

A leak that is a production defect becomes an issue; one that is a test's incidental dependence becomes a re-base.

**What the result decides in candidate D:** whether a claim's id carries its sub-cone (`<layer>.<sub-cone>.<concept>.<claim>`); whether a sub-cone is a node with its own composed claim (a certificate per sub-cone, under ruling 8); whether sub-cones carry edges to sub-cones, so that a sub-cone's composed claim composes from its upstream sub-cones' composed claims; and whether an edge across sub-cones is a different kind of fact from an edge within one (a cross-sub-cone edge is what couples two concepts, so it may need a reading by the orchestrator under ruling 18).

### Sub-cones: the result (the test-architect, 2026-10-18; `scratch/verification_tree/subcones/report.md`, scripts and coverage data beside it)

All readings are `[M]` on one run each: serial, `-O`, branch coverage, `-m "not slow"`, read through the coverage API; numerics 4 475 passed in 897 s (1 906 test contexts), sn 3 827 passed in 1 636 s (2 028 contexts). The order between groups comes from the static imports, independent of the tests.

- **S13. numerics: sub-cones partly real; the dependence is knotted.** Against the import order, 72 % of tests are upstream-closed pure, against 10 to 11 % with the labels shuffled (the working positive control). The top-level import graph of 44 modules is acyclic, depth 7, with 8 joins, among them `frame <- basis, measure, projection, space` (ruling 25's example, found by the instrument). But 7 lazy imports and 10 typing-only imports each close a cycle, giving cycles of 21 and 14 modules over axis, basis, measure, frame, quadrature, spaces, symmetry and invariance. Verified by the orchestrator: `FunctionSpace._axis_collapse_pair` lazily imports `frame._collapse_pair` (`orpheus/numerics/space.py:558`), so the space reaches its downstream frame; Nexus `native_place` independently places `_collapse_pair` on `FunctionSpace`.
- **S14. sn: real in the code, not separable in its tests.** The runtime module graph is acyclic (13 modules, 5 sources, 3 joins). But 1 312 of the 1 686 tests that run sn code execute the `SNProblem` hub, and the purity reading (91 %) does not beat its shuffled null (94 to 96 %), so it carries no information. Counting typing-only imports gives an 8-module cycle, and 5 of the 6 edges that close it are constituents typed on `SNProblem` while the hub builds them. Reading: SN's tests reach a constituent through a whole posed problem. That is either intrinsic to the posing filtration (a constituent exists only bound to its problem) or a weld (a constituent that cannot be built without the hub that builds it), and it decides whether hundreds of SN tests are re-base candidates.
- **S15. The module is the wrong grain.** Partition agreement is low (ARI 0.25 to 0.41 between module, concept-table and conceptual-view partitions); the concept partitions assign only 363 of 1 180 numerics defs; Louvain is unstable over seeds (13 to 22 communities, ARI 0.53 to 0.87) and coarser than modules. The pre-declared controls failed under the flat model, and by reading, the "pure" Legendre test is a composite and the "impure" Krylov test a unit test; replacement controls chosen by reading separate (1.00 against 0.22).
- **S16. The leak criteria, tested.** Numerics, 12 of 222 straddlers sampled: 0 production leaks; 11 cross a dependence the import graph omits (typing-only, lazy or duck-typed); 1 is a frame test built on a downstream quadrature (re-home the test). sn, 8 of 60: 7 boundary tests reach `SNMethodSpace.for_face`, which `realizer.py` imports for typing only; intrinsic. Of the five criteria, *mention* decided every case, *interface* and *twin* never fired, *direction* always fired as an artefact of the missing edges, and *counterfactual* disagreed with mention on cones a fixture built. Proposed replacement: the leak signal is a cycle in the runtime-plus-type dependence graph, and a test's cone splits into what its claim exercises and what its fixture builds.
- **S17. The existing `rests_on` edges are mostly not executes edges.** Strict cone containment holds for 26 of 77 in numerics and 0 of 14 in sn; at 0.95 or above, 39 of 77 and 7 of 14; the minimum is 0.002. The executes check must be fractional, and the near-zero edges are reads or supports edges (ruling 20).
- **S18. For candidate D** (the test-architect's proposal, with the orchestrator's reading): no sub-cone in the claim id, since membership moves as leaks are repaired, so it is a derived query; an edge across sub-cones is its own kind, checked mechanically along the dependence graph and read by the orchestrator when it goes against it or closes a cycle. The test-architect proposes that a sub-cone IS a strongly connected group of the runtime-plus-type graph. The orchestrator's reading `[R]`: a knot (a strongly connected group spanning several concepts) is a finding to judge, a leak or two concepts that are one, not a definition of a sub-cone. The concept grain S15 asks for is what candidate D supplies: each claim is declared at the concept that owns it, so the concept-grain sub-cone DAG is derived from the claims' homes as the census places them (ruling 22's sort), and module grain with typing and lazy edges serves meanwhile as the leak instrument.

**Delegation of the measurement** (ruling 18): one package per layer, `driver.sh` then the analysis scripts. Acceptance: the capture passes with passed count equal to the collected count; the null reading differs from the measured one; acyclicity is reported per edge set (runtime, typing, lazy); the output is byte-identical across seeded reruns. Straddler classification needs reading and stays out of the package. Caveats: single run; 100 production cache sites and 29 shared fixtures make cones depend on test order; child processes are not captured; slow tests are excluded.

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

### The DAG backlog (ruling 31)

Noticed and deliberately not finished, so the DAG line keeps priority; each is an issue, worked in a worktree of its own. Read each issue's title before acting on it (`process-discipline`).
- From `numerics_untangling.md`: #633 (S6c, moments as covectors; its order against the census is open), #624 (the measure's support guards), #612 (curvilinear LD, Gauss-nodal axis), #613 (Galerkin quadrature), the forgetful-map concept (no issue yet).
- From the S7 reviews (2026-10-21): #635, #636, #637, #638, #639.
- #611 (the constructor rule, ruling 28).

⏸ **Next (2026-10-21, ruling 30):** the numerics knot is settled when S7 of `numerics_untangling.md` merges. Then run the sub-cone measurement on `numerics` ("Sub-cones: the measurement", rulings 23–26; the cone store of S8 first), and read it as an ontology instrument before any census verdict. The reorganisation of the code into the DAG and its sub-cones is the priority it feeds. #611 (the constructor rule, ruling 28) is unchanged. (Superseded pointer: this campaign waits on `numerics_untangling.md` (ruling 27), whose concept order is the base of the census. #611 (the constructor rule, ruling 28) follows it. Rulings 26 to 28 recorded 2026-10-18. Rulings 19 to 22 recorded (2026-10-18). Candidate D is with the user for discussion, and ruling 22's reading for confirmation. Then: the four-file sample (TA F7), now run against candidate D (does each function's claim fit D's object, and which declaration site does it take?), and the cone store's specification by the test-architect, cut into packages per ruling 18. #604 is the deferred bound composition (ruling 13); #610 the memo's early cutoff.)

## The cone store — design (opened 2026-11-01, iterating with the user)

**What is still unruled.** Rulings 19 to 22 settled S1 to S5 and S7 to S10 on 2026-10-18. The pointers in `layer_dag.md` and `space_membership_covector_gauge.md` that say the cone store "waits on S1–S12" are stale. Unruled: S6 (evidence and confidence as two digests), S11 (early cutoff at a claim's statement), S12 (the environment in the verdict key), and the remaining half of open question 5: whether a declared claim's source is importable Python or a tracked file generated from the corpus.

**The object.** The cone of a test is the set of nodes whose change can change its verdict. The nodes are those of ruling 15's single change-detection layer, spelled with the traced memo's pin types (`orpheus/numerics/traced_memo.py`), so nothing is minted twice:
- a def, keyed by its module path, qualname and ordinal, and digested by `function_digest` (normalised AST, docstrings and layout removed): an `executes` edge;
- import-time code: today `ModulePin`'s skeleton digest, which any edit to a module's constants, decorators or class attributes moves, staling every reader of the module;
- a stored value read: a memo entry (`ChildPin`, by payload digest) or a data file (`DataPin`): a `reads` edge, which needs no code capture (ruling 20);
- the environment: `InterpreterPin`, `DistributionPin`, `EnvironmentPin` and the BLAS tag.

**Measured** `[M]` 2026-11-01. Probes, logs and data are in `scratch/verification_tree/cone_store/`. Runs are serial, `-O`, `OMP_NUM_THREADS=1`. Sample N is `tests/gates/numerics/test_invariance.py`, `test_measure.py` and `test_space_of_axes.py`: 200 tests. Sample E is `tests/gates/sn/eigenvalue -m "not slow"`: 87 tests.
- **The capture instrument.** A pytest plugin on `sys.monitoring` `PY_START` (Python 3.14) records, per test and per phase, the first-party code objects entered. It returns `DISABLE` after a code object's first entry and re-arms with `restart_events()` at each phase boundary, so each def costs one event per test. Pytest's own time:

  | mode | sample N | sample E |
  |---|---|---|
  | plain | 4.72–4.83 s | 133.4–135.9 s |
  | function entry, firing on every call | 5.11–5.17 s | 179.9–180.1 s |
  | function entry, re-armed per phase | 5.04 s | 136.5 s |
  | coverage, per-test contexts, branch | 7.47–7.54 s | 245.7–246.8 s |

  The re-armed and always-firing probes record the same cones on sample N (177 of 177 rows), and on sample E except one row, discussed next. Plain, always-firing and coverage are 3 interleaved repeats on N and 2 on E; re-armed is 1 repeat each.
- **In-process caches hide most of a cone.** Each test of sample N was run alone in a fresh process and compared with the same test run in the suite. 127 of 177 in-suite cones miss defs that the isolated run entered: a median of 35 defs per affected test, and up to 79. 0 of 177 in-suite cones hold a def the isolated run did not enter. The missed defs are cache fills: `GeneratingMeasure.recurrence` and `is_symmetric`, and the `Manifold` methods. On sample E, the always-firing probe gave two different cones for `test_keff_slab.py::test_heterogeneous_absolute_keff` on two identical runs, according to whether `StreamingCoefficientCache` was built inside the test or found by `__eq__`. So a cone captured in-suite depends on run state, not only on test order. The tree holds 240 cache decorators (`cache`, `lru_cache`, `cached_property`) by grep, and hand-written dictionary caches beside them.
- **Code that runs at collection is attributed to no test.** Parameters built in a `parametrize` list (for example `gauss_legendre(8)`) run at collection, outside every phase. 0 setup rows entered first-party code in either sample, run alone or in the suite.
- **Child processes.** 46 files under `tests/` name a process API (`subprocess`, `Popen`, `multiprocessing`, `os.fork`), by regex with a positive control. How many tests spawn a child is not yet counted.
- **The stored captures are stale.** Every `wide_*` coverage run in the Nexus runtime store predates the move to `tests/gates/`, since its paths read `tests/numerics`. None of its numbers describe the current tree.

**Positive controls**, written before the build (X1):
1. Mutating `_OUTFLOW_COSINE_TOL` (`orpheus/mesh/axis.py:469`, an import-time constant read by `face_outflow_ordinates`) selects exactly the tests whose executed code reads that name.
2. Mutating an oracle def behind a stored reference selects the oracle's validation test, plus every reader only through its `reads` edge.
3. A docstring-only edit selects nothing.
4. Ground truth: a sample of tests captured each alone in a fresh process. The store's cone must contain each isolated cone, with 0 misses on the sample.

**The questions, with recommendations** — round 1, for the user. Each is `[R]` until ruled.
- C1. How a cone is captured, given the cache finding.
- C2. Attributing import-time code finer than a module skeleton, shared with the memo.
- C3. The environment in the verdict key (S12).
- C4. Evidence and confidence (S6), and early cutoff at a claim's statement (S11).
- C5. Where a declared claim's source lives (open question 5).
- C6. Where the store lives.

**Ruled (the user, 2026-11-01): "I like the recommendations. I also endorse C5 suggestion. The cost seems very small for the benefits."** Rulings 32 to 37:

32. **C1. A test's cone is captured with the test alone in a fresh process.** That capture is the ground truth, and it also attributes the code that runs while the test is collected. The store fills layer by layer with the census (ruling 22). After an edit, only the tests whose rows went stale are captured again. Replaying a cache's fill cone stays a later optimisation, taken only if the capture's cost is measured and shown to matter. Capturing with caches switched off is refuted, for the question "how is a cone captured": its residual miss rate is unknown, because hand-written caches stay invisible.
33. **C2. Import-time code is attributed per top-level statement, and the memo shares the rule** (ruling 15). A test depends on a statement only when code it executed reads the name the statement binds (the names recorded in each executed code object). A name collision can only over-select. Positive control 1 decides whether the rule holds.
34. **C3. The environment is always part of a verdict's key:** the interpreter, the distributions, the BLAS tag and the thread count.
35. **C4. S6 (evidence and confidence as two digests) and S11 (early cutoff at a claim's statement digest) are accepted.** Both belong to the verdict layer built on the cone store, not to the store's first version.
36. **C5. A claim's single source is Python.** A Sphinx directive renders the claim on its theory page and anchors it there, as autodoc renders a docstring. Tests read claims at collection without a docs build. A tracked file generated from the corpus is refuted for this question, because it is a second spelling that would need its own drift gate.
37. **C6. The store is tooling, under `tools/`, with its data in SQLite under `.cache/`.** It reuses the traced memo's digest functions, and Nexus ingests it as a runtime kind. The store does not own those digests, and Nexus does not own the store.

**Round 2 measurements** `[M]` 2026-11-01, sample N (200 tests), serial, `-O`, `OMP_NUM_THREADS=1`; probe `scratch/verification_tree/cone_store/cone_fork.py`.
- **A fork per test from a collected parent gives the isolated cone.** The parent collects every test, then forks one child per test. The child re-arms the monitor, runs the test's setup, call and teardown, and writes its cone. 200 of 200 children passed. 200 of 200 forked cones are identical to the cones of the same tests run alone in fresh interpreters (setup plus call, first-party defs under `orpheus/`). The parent held 1 OS thread after collection, which satisfies the memo's fork guard (`_fork_unsafe` in `orpheus/numerics/_traced_memo_boot.py`: a fork after a BLAS call crashed its child 20 of 20 times, `[M]` 2026-10-16).
- **Cost per test.** Fork: 57.4 s for 200 tests, 0.29 s each, against 0.03 s each plain. A fresh interpreter per test: 21.5 s for 20 tests of `test_measure.py`, 1.07 s each. Where the fork's 0.26 s of overhead goes is not yet measured. Extrapolated `[R]` to about 20 000 tests: about 1.4 h of fork overhead on top of the suite's 56 min, paid once per environment, then only for stale tests.
- **Neither probe records collection.** Both the isolated and the forked probes record setup and call only, so code run while a test is collected (a `parametrize` list) is in neither cone yet. The probe's root is `orpheus/` only, so a test's own helpers under `tests/` are not recorded either.

**Ruled (the user, 2026-11-01): "Agreed. You can go ahead."** Rulings 38 to 41:

38. **D1. The capture forks each test from a collected parent.** The parent is guarded by the memo's single-thread check (`FORK_MEASURED` and `_fork_unsafe`, reused, never re-spelled). When the check fails, the file's tests are captured in a fresh interpreter each. ~~Code that runs at collection is recorded in the parent once per test module and attributed to every test of that module, which can only over-select.~~ [REFUTED 2026-11-01 for the question "how is collection attributed"; the FACT it establishes: caches in a single collecting parent leak between test modules.] `[M]` the test-architect, 2026-11-01, on a fixture tree (a `functools.cache` filled by one module's `parametrize` list and hit by another module's test; a constants module first imported by an earlier module's collection). In the suite, the later test's cone lacks the fill's defs and the modules' import-time code; run alone, it holds all five. There are two mechanisms, both in the parent: a cache filled during an earlier module's collection, and `sys.modules`, which runs a module's import-time code only for its first importer. The 200 of 200 measurement compared setup and call defs under `orpheus/` only, so it could not see this. **The repair, ruled by the orchestrator under ruling 38's intent (the isolated cone):**
   - a grandparent imports pytest and the conftests;
   - per test module it forks a parent that collects that module alone, and that parent forks one child per test, with the fork guard at each level;
   - what the grandparent imports or fills is attributed to every test, which over-selects correctly.

   Attributing the union of all collection cones to every test is refuted for the same question: any skeleton edit would select the whole suite.
39. **D2. A row is (test, environment key, node, digest at capture).** The record covers first-party code under `tests/` as well as under `orpheus/`. A row is fresh while every digest matches the checkout, checked in the memo's two stages: file identity first, then the normalised-AST digest. Stale and new tests are captured again by fork.
40. **D3. A test whose child process is unrecorded is marked "cone unknown" and is always selected.** Recording inside children waits until the number of tests that spawn one is counted and shown to matter.
41. **D4. The first client is a local selection command:** given the working tree, it lists the tests to run. The census uses it to place tests (ruling 22), and it serves as a fast local check. The full suite stays the merge gate until the validation run passes (phase-order step 3): sampled mutants, each run both in full and through the store, with every selection agreeing.

### The cone store — build (opened 2026-11-01)

Mode: W1. The test-architect writes the verification spec first. The store is tooling under `tools/cone_store/`, built on its own branch. B2 changes the memo's pins in production, and the main agent writes it (the surgical-carve exception).

- **B1. The machinery, at today's node grain** (def digests and module skeletons, plus the `reads` rows: `ChildPin` for a memo entry read, `DataPin` for a data file read, recorded by reusing the memo's `Recording`, not by a second recorder):
  - a pytest plugin that collects, records the collection cone per module, forks per test behind the single-thread guard, and records each child's cone over `orpheus/` and `tests/`;
  - in the child, a flag raised on any process start, which marks the row "cone unknown";
  - the environment key (ruling 34);
  - the SQLite store under `.cache/cone_store/`, with the freshness check reusing `function_digest`, `skeleton_digest` and the memo's file-identity stage;
  - a command with `capture`, `select` and `status`.

  Gates: positive controls 2, 3 and 4 as tests, each with the mutation that reddens it. Control 1 holds at module grain in B1, a superset of the readers; it becomes exact in B2.
- **B2. Per-statement import-time pins, shared with the memo** (ruling 33). Gate: control 1 exact, the memo's own validation suite, and a measured count of memo entries left warm by a constant edit, before and after.
- **B3. The numerics layer captured into the store, and the census's placement step reading it** (ruling 22). Nexus ingestion of the store is an issue on sphinxcontrib-nexus.
- **B4. The validation run** (ruling 41), then the verdict layer (rulings 34 and 35).

**B1 [LANDED 7fdc7474]** (2026-11-02, merged to `main` by fast-forward; CI `gates` run 41425672218 green). `tools/cone_store/` (13 modules), seven `tests/gates/tools/test_cone_store_*.py` files with their fixture, three memo witnesses, `docs/theory/verification/cone_store.rst`, ERR-124. `[M]` before the merge: full suite `-O -m "not slow"` 19947 passed, 1 failed; the failure was the qualname census (the memo's new `Refusal` shared its qualname with `orpheus.geometry.manifold.Refusal`), fixed by naming it `PinRefusal`; after the rename the census and the memo witnesses passed (37), the cone-store fast gates passed (181), M1 passed, pyright 0 errors, Sphinx `-W` clean. Scope boundaries filed: #687, #688, #689 (B2), #690.

**B2 [LANDED 6e52b3cf]** (2026-11-02, merged to `main` by fast-forward with the agent memory `867c09c3`; CI `gates` run 41532535468 green; closes #689). The rule `orpheus/numerics/_import_time.py`, the memo's `StatementsPin` and `ImportOutcomePin` (rulings 43 and 44), the cone store at `FORMAT = 2` reading them, and the gates `test_traced_memo_statement_pins.py` and `test_cone_store_statement_grain.py`. The done-when, `[M]` 2026-11-02, serial, `-O`:
- control 1 exact (SG4, slow, on a clone of the real tree: 6 tests selected where B1 selected 9);
- the memo's validation suite and the slow tools gates: 4 passed, 322 deselected (78 min); the not-slow suite green in two runs (the first was killed at the 2 h limit at 91 % during the cold regeneration, and its only reds were 3 concept-order rows, fixed by ranking `_import_time` beside `_traced_memo_boot`; the tail: 1658 passed, 1 skipped, 1 xfailed);
- the warm count, `scratch/verification_tree/cone_store/b2/warm/after.py`. It judges each stored B2 entry's statement pin on the edited module with the memo's own projection, against an in-memory edit; no file is written. The positive control: the unedited source reproduces every pin (2021 of 2021, per edit). Before (module grain), each edit makes stale every entry that imported the module: 2021 of 2021, for all five. After:

  | edited constant | how runs read it | warm after |
  |---|---|---|
  | `eigenvalue.MINIMUM_OUTER_ITERATIONS` | no run reads it | 2021 of 2021 |
  | `metric._DENSE_METRIC_RCOND` | no run reads it | 2021 of 2021 |
  | `transformation._MAX_GROUP_ORDER` | read in a def body that 2016 runs started | 5 of 2021 |
  | `measure.WEIGHT_ATOL` | a default argument's value, evaluated at import | 0 of 2021 |
  | `convergence._SERVED_RATE` | a default argument's value, evaluated at import | 0 of 2021 |

  The edit also changes the import set's bytes. So the run's import outcome is re-measured once per distinct import set: 4 sets over the 2021 entries. A value edit that neither raises nor warns leaves the outcome equal `[R]`.
- The 691 entries the design could not read: an older `DistributionPin` arity. Like the 5814 `ModulePin` entries, they read as stale and regenerate on demand. 2021 of the 8526 stored entries were regenerated in B2 form by the suite run.
- pyright: 0 errors in the 33 changed or new Python files. Sphinx `-W` clean. `dead_references` unchanged (the 3 pre-existing transport test-helper targets).

**B3 [LANDED 9195dcf5, b65674fc, f78085d8, 1e989da7]** (2026-11-03 to 2026-11-04, each merged to `main` by fast-forward, CI `gates` green on each push; the error entry ERR-125 at `612e597a`). B3 itself (`9195dcf5`): two cones per test, the selection cone `select` reads and the own cone `place` reads (`FORMAT = 3`), `tools/cone_store/placement.py`, and the recorder's apparatus kept out of every recording. Its acceptance found three defects, each fixed with gates that a mutation battery reddened:
- **Over-selection, X1.3** (ERR-125, `b65674fc`). A no-edit `select` over the numerics layer listed 229 tests where 149 were owed. 70 readers of a memo store under their own `tmp_path` held a `ChildPin` judged against the default store; a child of another store is now pinned by the pins its lookup validated it on (`Hit.pins`, `ReadChild`), placed where the read found them. 10 source-walking tests held a `__pycache__` listing, which is no longer observed.
- **T3's own predicate** (`f78085d8`). The ground-truth gate T3 compared statement pins by their per-run name sets, which depend on test order; it now compares the statements each pin observes. `[M]` on the kept stores of a passing run: 0 misses over 200 tests, with a control at 190 of 200.
- **Own cones holding another process's defs** (`1e989da7`). The pins of a foreign read had entered the own cone, so placement homed the 70 readers by defs only the generating process ran. `CollectedPins` now separates `ran` from `through_reads`, and `cone_of` returns `ReadCones(ran, selection)`. Gates S23 (own leg) and S24.

The acceptance, `[M]` 2026-11-04 (store and logs under `scratch/verification_tree/cone_store/b3/recapture/`):
- **Capture** of `tests/gates/numerics/` at `f78085d8`, `OMP_NUM_THREADS` unset, in a fresh worktree with a copied warm `.cache`: 5 808 tests in 13 chunks of 10 files, 4 h 56 min, every chunk exit 0; 148 cone-unknown (140 start a process, 7 read a `.cache` file the capture changed, 1 is the gate that deletes a foreign child's source on purpose), 2 not passed (skipped).
- **X1.3 met.** A no-edit `select` at the capture's commit listed 150 tests, exactly the 148 unknown and the 2 not passed, with nothing stale. On `main` (`612e597a`) it listed 198: those 150, the 38 tests of the two files that gained `catches` markers after the capture, and 10 tests that list `orpheus/` or `orpheus/data/micro_xs`, where `main`'s checkout holds gitignored files the worktree lacks (`.DS_Store` and 13 generated `.h5` files).
- **Placement** (`place`, threshold 0.95, `min_defs` 3), after the own-cone fix, on a copy of the store with the six memo test files re-captured at `1e989da7`: 5 257 of 5 808 placed, 993 straddlers, pure fraction 0.811; unplaced 148 cone-unknown, 401 with no layer def, 2 not passed. Against the previous capture (B3's, 5 803 tests) 0 of the 5 803 shared tests changed home; before the fix, 70 had.

Findings filed, not acted on:
- **A store is not self-contained** (#697). The kept name and import sets a store's pins name live in the `.cache/sets` and `.cache/import_outcomes` of the checkout whose `orpheus` the capture imported. A `select` from another checkout read 5 658 of 5 808 tests as stale ("is not kept") until the worktree's sets were copied into `main`'s `.cache` (content-addressed, so the copy only adds).
- **qa's F3 on the own-cone fix** (#698). A refusal raised by a foreign child (`SourceGone`) makes its reader cone-unknown, and placement then leaves it unplaced, although its own cone is complete; placement reads the selection cone's unknown reasons only. `[M]` 1 of the 148 cone-unknown tests, the gate that deletes the source on purpose.
- **A listing records macOS's `.DS_Store`** (#698). A test that lists `orpheus/` is selected when Finder writes that file.


**B4, ruled by the user (2026-11-04)**, on the test-architect's spec (`scratch/verification_tree/cone_store/b4/spec_b4.md`, driver beside it). B4 validates at `42b5ceca`, against a store of the numerics layer captured at that commit in the worktree `ORPHEUS-b3-recapture`, because `select` lists every row of another environment key (`tools/cone_store/selection.py:140`).

45. **The validation run's terms.**
    - **The bound is 5%.** At least 60 counted mutants must redden a test, and 0 misses over them bounds the per-mutant miss rate below 5% at 95% confidence.
    - **Every mutant's selected subset is also run alone,** and its verdicts are compared with the full run's. A disagreement is reported, not gated.
    - **A miss in stage A stops the batch** to be fixed. Stage A is the def-body mutants, one or more per tier.
    - **The ground truth is the full-suite context.** Every run collects the whole testpaths tree, as the merge gate does, so all 12 conftests and every test module are imported, and only the numerics tests run. The store promises the verdict a full run gives (ruling 42), and a run of the numerics folder alone loads 1 conftest, so it cannot show another layer's effect on a numerics verdict.

**B4 [MEASURED 2026-11-05 to 2026-11-06, at `42b5ceca`]**: soundness holds on the numerics layer. The spec, driver and records are under `scratch/verification_tree/cone_store/b4/`: `spec_b4.md`, `driver/`, `runs/` with `batch.log`, `q1_q5.json`, `q2.json`, `red_sets.json` and one directory per mutant. The store is `b4/capture/numerics.sqlite`, captured in the worktree `ORPHEUS-b3-recapture` (kept: its `.cache/sets` hold the kept sets the store names, #697). `[M]`, 85 mutants, serial, `-O`, the full-suite context (all 12 conftests, 5 808 numerics tests kept), 42 h 14 min wall (2026-11-05 01:37 to 2026-11-06 19:51):
- **Soundness: 0 under-selections over 61 counted killing mutants** (the user's floor was 60): the per-mutant miss rate is below 3/61, about 4.9%, at 95%. The miss counter's positive control (one red id dropped from a copy of a selection) read its miss. Not counted, each with its reason in `q1_q5.json`: 13 mutants whose reds were all in the no-edit set, 4 that cannot miss (`select` refused, or selected every test), and 2 floor rows, 3 nulls and 2 controls.
- **Restores:** 85 of 85 witnessed on all three legs (bytes, `git status`, a pristine `select` equal to the no-edit set). The first, second and final pristine runs agree; the flaky set is empty.
- **Controls:** the positive control (the root conftest's report hook) selected exactly the 5 808 ids and reddened more than 5 000; the conftest import refusal made `select` exit non-zero.
- **Precision:** a killing mutant selects a median of 72 tests beyond the 150 always selected, out of 5 658. Selected beyond S0 over red beyond S0 has a median of 8.3, a minimum of 1.0 and a maximum of 592.
- **Q5, the selection run alone:** its verdicts agree with the full run's in 77 of 77 mutants. The one flag, on `T07b`, is an artefact of the record: both runs hold the same collection error of `test_inverse_metric_operator.py`, keyed by the file's path, which the full side's filter on selected ids drops.
- **Q2, placement predicts propagation:** 0 defects of the placement kind over 59 mutants. One red test's own cone lacks the mutated node (`S10`, `test_axis_identity.py::test_K1…`): it is a census over `__subclasses__()`, which depends on which modules were imported, and an import runs module code, which the own cone records as statements, not defs. Selection covers it; the census of placement must read it as a dependence through import, not through a def.
- **The one failed acceptance predicate, Z1, was a wrong prediction.** A docstring edit in `generating_measure.py` selected 158 tests where the spec predicted the 150 of the no-edit set, and reddened none. All 8 extra tests carry the reason "data … generating_measure.py changed": they read the module's source as data (the concept-order census, import and AST audits), so their verdict can depend on the text, and selecting them is right. The corrected predicate (the no-edit set plus the tests holding a `DataPin` on the file, nothing reddened) holds for 8 of 8.
- **Cost:** a mutant's `select` took a median of 248 s, and its full run 874 s (#699: `select` with one edited file is nearly as slow as the run it saves).

46. **The merge gate stays the full suite** (the user, 2026-11-07, option A). B4 validates the numerics layer only, 5 808 of about 20 750 tests, and #699 makes `select` too slow to save much yet. So `select` is the trusted local check for an edit to the numerics layer, and the full suite stays the merge gate until every layer is captured and validated.

⏸ **Next (2026-11-07): the census of the numerics layer** (ruling 22), reading the placement of the store at `42b5ceca` (`scratch/verification_tree/cone_store/b4/runs/place.json`: 5 257 of 5 808 placed, 993 straddlers, pure fraction 0.811). The lists it works through, in order:
1. the 401 tests whose own cone holds no numerics def: they are probably filed in the wrong folder;
2. the 126 proposed supports that run against the concept order: each is a misplaced test or a wrong order;
3. the declared `rests_on` edges against execution: of 798, 138 hold at containment 1.0 and 513 below 0.8, and 275 cannot be measured (221 whose lower test ran no def). Decide whether the declarations are wrong or containment is the wrong test for that kind of edge;
4. the 991 pairs that ran the same defs, and the 105 near-duplicate groups: candidates only, since a duplicate needs equal kill sets and the same claim (ruling 21);
5. the 993 straddlers, read against the `S10` finding above (a dependence through import is not a def).

Then the census's reading pass (each test's claim) and the folder's reorganisation, then the next layer (ruling 31). #697, #698 and #699 are the cone store's open items; #699 comes before `select` is used day to day. The verdict layer (rulings 34 and 35) waits for the census.

**List 1 [MEASURED 2026-11-07]** (`scratch/verification_tree/cone_store/b4/census/l1.md`, with its scripts and JSON beside it). Two instruments: what each test ran (its own cone, by pin kind and package), and what it names (a module-qualified static read of each test function, its helpers and fixtures; positive control: the placed home is among the named numerics nodes for 5 044 of 5 237 placed tests, the 193 misses being 184 `traced_memo` tests reached through a subprocess or a fixture factory, 8 `kernel_projection_gauge`, 1 `measure`). The 401 tests (160 functions) fall into four populations:
- **G, misfiled: 205 tests, 47 functions.** They ran only `orpheus/geometry` or `orpheus/values` defs and name no numerics module. Their subjects (the subgroup lattice, G-spaces and quotients, the face-name bijection, geometry retirements, the memo key's encoding) moved below numerics in `830b3713` and the tests stayed. Verdict: relocate, per function (each of the 7 files also holds numerics tests).
- **S, correctly filed, invisible to placement: 139 tests, 68 functions.** They assert declared structure (dataclass fields, enum members, Protocol membership, constants, registry population, a value type's content identity, retired names). What they read is a class statement or an assignment, which ran at collection, and `collect_pins` mints statement pins only for files that have a code object in the recording, so the own cone holds none. Positive control: the 4 tests whose own cone does hold a numerics `StatementsPin` are all in S. This is B4's `S10` finding seen from the unplaced side.
- **A, correctly filed, invisible to placement: 19 tests, 19 functions.** Source audits that read a numerics file as data (a `DataPin`); `place` homes by `DefPin` only. 11 read 1 to 3 nodes; 8 read 50 or 51, the whole layer.
- **O, tests of the test apparatus: 38 tests, 26 functions.** They name no `orpheus` module: oracle self-tests (the reference Gauss rules against `mpmath`, their frozen digests, the gate-teeth mutations of a test-local reference), fixture premises, the concept order's own reader, retirement gates. 0 of the 26 functions is the lower test of a declared `rests_on` edge.

**Ruled by the user (2026-11-07), on list 1.**

47. **Placement homes a test by what its run read, not only by the defs it ran.** The own cone gains the module statements whose bindings the run read, resolved module-qualified (where the read resolved, never by a bare name), and `place` reads a source file read as data as a node of the test. Measured consequence to explain before the re-placement is believed: every home that moves among the 5 257 placed tests.
48. **A test of the test apparatus is a node below the gates that use it.** An oracle, a fixture or a case table is a node of the DAG: its self-test homes there, and every gate that uses it rests on it. A premise test is the declared support of the gates that need the premise. A retirement gate homes at the node the retired name lived in.

**B5, the placement extension (ruling 47): the design** (the orchestrator, 2026-11-07; the test-architect's spec is owed before the build).
- **Where the reads are resolved:** in `cone_of` (`tools/cone_store/cone.py`), on the own recording, in the worker after the test, with `sys.modules` live. The memo's `collect_pins` is not changed. A generation imports inside its own recording, so its `<module>` code runs and its files are indexed already. The selection cone is not changed either: the collection's share holds the statements, and B4 measured selection sound with it.
- **The resolution:** for each executed first-party code object, its module is the one in `sys.modules` whose real `__file__` is the code's file. Each name the code reads (`_import_time.code_names`) that its module binds resolves through that module's import statements to the first-party module whose non-import statement binds it, following re-exports. An attribute read on a bound module (`measure.WEIGHT_ATOL`) resolves in that module. A name that resolves outside the first-party package is dropped.
- **The pin:** one `StatementsPin` per resolved module, projected on the names resolved into it, so the row validates by the existing projection rule.
- **Placement:** `own(t)` becomes def nodes ∪ statement-read nodes ∪ data-read nodes, each kept apart in the report (`via`: `def`, `statement`, `data`). A test whose data reads span the whole layer, with no def or statement node, is homed at the layer itself (a new home value), never at the highest-ranked file it read. Containment for proposed supports stays over defs (ruling 4's claim grain is unchanged); whether statement reads join containment is decided by the re-placement's numbers.
- **The apparatus (ruling 48):** a first-party test helper (an oracle module, a fixture, a case table under `tests/`) becomes a node of the placement, beside the layer's nodes, upstream of every node that its users home at. Its grain (a helper module, or a fixture function) is the spec's to propose.
- **Acceptance:** re-capture the numerics layer at the build commit (about 5 h, detached). Then: the 158 S and A tests are placed; the 4 controls with a numerics `StatementsPin` are placed where the static read names them; every home that moves among the 5 257 placed tests is listed and explained; the null shuffle still reads lower than the real placement.

**B5's spec** (the test-architect, 2026-11-07: `scratch/verification_tree/cone_store/b5/spec_b5.md`; probes and a prototype resolver run on 658 tests at `42b5ceca` under `b5/probes/` and `b5/proto/`). 22 attacks, 14 corrections (C1 to C14), 49 gates in six groups (R resolution 17, P pin 8, F floor 2, V placement channels 10, D data home 3, A apparatus 9), acceptance AC0 to AC10 on a re-capture. Measured on the prototype `[M]`: with no closure, 65 of 139 S tests get a layer node; with closure through apparatus only, 128; the union home moves 110 of 197 sampled homes and 0 of 60 controls.

49. **The home is tiered** (the user, 2026-11-07, on the spec's C5). Defs decide the home when the test ran any layer def; otherwise its statement reads; otherwise its data reads. Every channel is reported under `via`; the union home (`home_union`) and the statement-channel straddle are reported under their own keys and do not change purity. Reason: running a module's code is a dependence on its behaviour, while reading a value with no code of its module run is a dependence on an input. The union reading is refuted FOR the default home; the fact it establishes is where a test's inputs come from.

**Ruled by the orchestrator (2026-11-07), on the spec.** C1 to C4, C8, C10 to C12 are adopted as written. C6 is adopted as refined: closure is followed only through test-tree bindings, to the first layer binding. On the spec's NEEDS:
- C7: a re-exporting layer module is pinned, not a node of the test (3 such statements in the tree). The subject is the binding the chain ends at; the re-export is the path the import took.
- C9: a new store table `own_reads(test_id, path, name, channel)`. Placement must not read name sets that live in the capturing checkout's `.cache` (the direction of #697), and `collect_pins`' statement pins name the run's whole closure, not what the test read.
- C13: a data read that walks all of `orpheus/` homes at `package` and is a census finding, listed apart.
- C14: the apparatus node is a test-tree file and its items are bindings. A test that read bindings of two apparatus files and none of its own is a tie, like a tie between layer nodes.

**B5 [BUILT, not yet committed] (2026-11-07)**, on `feature/b5-placement-reads` over `4fad267e`. Built by the main agent (W3): `tools/cone_store/reads.py` (new: each executed first-party code object's load chains read from its bytecode, fused instructions split by a table derived from `dis.opmap`, cells resolved through the def that binds them, each chain resolved to the statement that binds its end with re-exports followed by identity, a name no statement binds counted as `unresolved` and kept in the pin; the `Channel` enum with `is_node` and `is_pinned`; `Sources` reused and `resolve(...) -> RunReads` per test); `cone.py` (the reads resolved for the own recording under `unrecorded()`, one `StatementsPin` per module read in the own cone); `store.py` (`own_reads(test_id, path, name, channel)`, `FORMAT = 4`); `worker.py` and `driver.py` (the reads carried as `Read` rows); `placement.py` (`OwnReads`, the tiered home as the `Home` sum, `home_of`, `via`, `items`, `home_union`, `statement_straddle`, `apparatus_items`, `apparatus_supports`, `c_items` and `declared_items_only` on the declared edges, the statement-channel purity and null, `None` for a fraction over no test). Gates: `tests/gates/tools/test_cone_store_reads.py` (R, P, F) and `test_cone_store_placement_reads.py` (V, D, A) on the fixture `_READS`; PL8's reason is now "no layer read". Gates `[M]` 2026-11-07 (the archivist, the three files serially, `-O`): 74 passed, 0 failed, then V11 (a form-3 store refused naming the own reads) passed in its file's re-run, 23 of 23; the two new files hold 58 ids (R 22, P 10, F 3; V 11, D 3, A 9), the spec's 49 plus R7b, R18 to R21, P3b, P9, F3 and V11, added from the battery and the reviews. Review: qa (`b5/qa/qa_b5_review.md`) and the elegance-enforcer (`b5/elegance/review.md`). Repaired in the build: qa F1 (`judged` restored to every passed test with a known cone; the reason renamed "lower ran no def and read nothing"; `declared_items_only` counted), F2 and the elegance F8 (the fused instructions, gates R18 and R21), F3 (the cell table keyed by enclosing def, R20), F6 (`None` over an empty population); the elegance F1 (`Home` sum), F2 (`holding` called), F3, F4 (`Channel`, `Read.to_row`/`from_row`), F5 (`RunReads`; an unbound name met under `closure` is pinned only, R19 and P9), F6, F8 to F11. Filed as #700: qa F4 (four silent drops) and F5 (identity on a shared immutable), the elegance F7 (the test's def spelled two ways) and F12 (the second `StatementsPin` mint). The theory page's B5 section is `docs/theory/verification/cone_store.rst`, "Placing a test by what its run read (B5)". Mutation battery on the final code `[M]` 2026-11-07 (the test-architect's, in a separate worktree, `b5/battery/`): 41 arms, every one RED; 58 of 58 test functions of the two new files reddened by at least one arm; the unmutated build 0 red; the positive control (attribute chains dropped) reddens 23. Open: the acceptance AC0 to AC10 on a re-capture of the numerics layer at the build commit.

⏸ **Next (2026-11-07): B5 before lists 2 to 5.** Each later list reads the placement, and list 5 reads the straddlers against exactly this finding. The test-architect writes B5's spec, the main agent builds it (surgical carve), and the re-capture runs detached. Relocating population G waits for the folder reorganisation, so a test moves once.

**Ruled by the orchestrator, from ruling 41 (2026-11-01).**

42. **The store predicts full-run verdicts.** The full suite is the merge gate, and the validation run compares the store's selections with full runs. So a capture reproduces the context a full run gives each test. A full run imports every `conftest.py` before any test runs, and a conftest's module code and session hooks act on every test (`tests/gates/sn/regression/conftest.py::pytest_configure` sets the session's `DriftWarning` filter, and 5 test files outside that directory raise `DriftWarning`, `[M]` qa 2026-11-01). Therefore:
    - the grandparent loads every conftest the full run would load, and attributes what they run to every test;
    - the set of conftest paths is part of the environment key;
    - so are the `pytest11` entry points and `PYTEST_PLUGINS`.

    A per-module verdict is refuted for the question "what does the store promise", because no gate reads one.

**Review state** (2026-11-01). qa's first pass (`qa_b1.md`) and second pass (`qa_b1_delta.md`) found 13 under-selections and 1 refusal-order defect. The repairs are in, or are being written by the test-architect. Defects found in the memo by this work, all fixed in production with witnesses:
- `Recording.codes` collapsed equal code objects across files (ERR-124);
- `trace_call`'s manifest re-read was offered to an enclosing recording;
- `listing_digest` counted `__pycache__` (now an accepted limit of the generator contract);
- a deleted source raised out of `collect_pins` (`SourceGone`).

### B2 — the design (opened 2026-11-02, #689)

**The question.** Today a run that executed any code of a first-party module pins the module's skeleton: every top-level statement, with def bodies hollowed (`ModulePin`). So an edit to one constant makes stale every memo entry and every cone that imported the module. B2 pins only the import-time statements the run could observe.

**Measured premises** (`[M]` 2026-11-02; scripts and output under `scratch/verification_tree/cone_store/b2/`):
- The memo store: 7 835 of 8 526 entries readable (691 unreadable, to be explained before the after-count). The median entry pins 131 modules by skeleton, and 33 modules are pinned by every entry. So today a constant edit in any of those 33 makes all 7 835 entries stale. `orpheus/mesh/axis.py` (control 1) is pinned by 0 entries; the memo's warm count is therefore taken on a constant of a module every entry pins.
- The census of import-time code: 415 files (`orpheus/**` and every `conftest.py` under `tests/`), 12 339 top-level and class-body statements. Each predicate below was asserted against a positive control first.
  - 11 bare expression statements (a call that binds nothing);
  - 27 stores into a subscript or an attribute (they mutate another statement's binding);
  - 208 other statements whose import-time evaluation calls first-party code, and 35 defs or classes with a first-party decorator (34 `withdrawn_generator`, 1 `traced_memo`); a call can mutate state no binding names;
  - 13 class-body dunder attributes, which the interpreter or a C library reads (`__array_priority__`, `__hash__ = None`, `__match_args__`), never Python code by name;
  - 36 computed-name `getattr`, 2 `hasattr` and 3 `setattr` calls, 4 `vars()`, 2 `import_module`, 1 `.__dict__`; 0 `import *` (filter controlled on `from a import *`); 118 compound top-level statements (`if`, `try`, `with`).

**Ruling 33's sub-claim "a name collision can only over-select" is [REFUTED 2026-11-02] for the question "is the name rule sound".** The fact it establishes: an import-time statement can change a run without binding any name the run reads. The shapes are the bare expression, the store into a subscript or an attribute, an import-time call into first-party code (a registering decorator), a class dunder read by the interpreter, and a name built at run time. The rule also left open whether a module's own `<module>` code counts as code the run executed. If it counts, that code reads nearly every name its module binds, so B2 gains nothing. If it does not count, every read made at import (a default argument, a derived constant, a re-export) is lost. **The repair, under ruling 33's intent** (the orchestrator, after an advisor review):

1. **The partition `_def_pin` already makes decides which code is read dynamically.** The entry names E of a run are the union, over every code object it started, of its `co_names` and the identifier-shaped strings in its `co_consts` (recursively through tuples; this catches `getattr(x, "name")`, `hasattr` and a `match` class pattern's keywords). This union covers stdlib, third-party and generated code too: a dataclass's generated `__init__` names every field. The one exception is first-party module-level code (`<module>`, class bodies, signature scopes), which is covered statically by step 2. A first-party module-level lambda or comprehension that ran after import also selects the statement whose source holds it, by line.
2. **Module-level code is read statically, per statement.** A statement is a top-level statement of a module, or of a class body (recursively), with every def body hollowed (the body is the `DefPin`'s). For each statement the index keeps its digest, the names it binds B, and the names it reads R: the `Name` loads and `Attribute` attrs evaluated at import, including decorators, defaults, bases, annotations and lambda bodies. The selected set is the least fixpoint:
   - a statement is selected when B meets the closure names E′, where E′ is E plus R of every selected statement;
   - an `ImportFrom` that is selected puts its source names in E′, so a re-export is followed across modules (names are not module-qualified; a collision over-selects);
   - a class header (its decorators, bases, keywords) is selected with any statement of its body, and its dunder statements are selected with the header;
   - **effect statements** are selected for every run in which their module was imported: the bare expressions, the subscript and attribute stores, and the statements whose import-time evaluation calls first-party code (decorators included). The 281 measured effect statements keep module grain; the other 12 058 are selected by name.
3. **The pin is the module's projection.** One `StatementsPin(path, names, digest)` per module replaces `ModulePin`. Here `names` is E′ restricted to the identifiers occurring anywhere in that module's source at capture. `digest` is the digest of the ordered list of the module's selected statements, where selection is re-run with the stored `names`. Validation repeats the projection on the current source. So an edit to a selected statement makes the pin stale, and so does a new statement binding a stored name (a shadowing def, a new dataclass field read by name), a removed one, a reordering, and any change to an effect statement. An edit to an unselected statement does not. ~~The `ModulePin` class stays, so old JSON rows still parse until the entries regenerate.~~ [REFUTED 2026-11-02 by ruling 44: the skeleton and `ModulePin` are retired.]
4. **Declared boundaries, each with its count:** a name built from data and never present as a string constant in executed code (the 36 computed `getattr` sites are read one by one; those whose names come from a constant tuple are covered by step 1); a new top-level name added to a module, read only through a guarded probe (`getattr(module, name, default)`, `hasattr`) and absent from the module's source at capture; `vars()` and `__dict__` iteration over a module namespace; and C code reading a non-dunder class attribute.

**Gates (the test-architect's spec, owed before the build):** control 1 exact; a recall fixture for every shape in the census, each with the mutation that reddens it; the memo's own validation suite; the warm count before and after on a named constant edit.

**Who writes what:** the memo part (`_SourceIndex`, `StatementsPin`, `collect_pins`) is written by the main agent (surgical carve). The store's `pins`, `cone` and `selection` follow it.

**qa's adversarial pass on the B2 design** (2026-11-02; `scratch/verification_tree/cone_store/b2/qa_b2_design.md`, with a model of steps 1–3 and 42 attack rows in `attacks/`). `[M]`: 24 of 42 rows under-select under the design as written (every pin valid, the run's result changed). The positive controls behave: control 1 goes stale, and control 2 (an unread constant's edit) stays valid with an unchanged result, which is the benefit B2 is for. The realized shapes, by population in the tree:
- a module-level constant read only by a `parametrize` or `fixture(params=…)` decorator: pytest finds tests by iterating the module, so the name enters no `co_names`. 371 such names in 183 of 769 test modules; the design's census had left `tests/**/test_*.py` out (26 150 statements);
- an unselected statement that starts raising at import (deleting a def that an unread table still names). Every unselected statement is exposed;
- the projection and class bodies: a subclass override in another module (159 classes with a first-party base in another module), a new dataclass field read by a field-iterating reader (462 dataclasses, 104 such readers, among them the content identity at `orpheus/values/content.py:356`), a new enum member (22 enums);
- a subclass added to an `__init_subclass__` registry that `orpheus/sn/operators/loss_kernel_gauge.py:430` iterates;
- `orpheus/data/group_structures/wims.py:124-134`: a table read only by a module-level loop, saved today by a name collision on the loop variables.

Refuted attacks (structural reasons in the report): walrus, `global` in a def (`STORE_GLOBAL` puts the name in `co_names`), star targets, PEP 695 `type`, augmented assignment, `__slots__`, `__set_name__`, module-level lambdas, re-exports, `import a.b as c`, field reordering. The grandparent loses nothing (one `Recording` is inherited across both forks), but it inflates E: 24 411 names against 245 for the test's own code.

**Ruled by the orchestrator** (2026-11-02, under ruling 33's intent, after ruling 38's precedent; the repairs measured together in qa's model leave 0 of 42 rows landing and keep control 2's benefit):

43. **The B2 rule, amended.**
    - **R1.** A first-party def that started selects the statements binding its qualified name's parts.
    - **R2D.** A statement is an effect, and is selected for every run that imported its module, unless it is a pure binding. Impure means any of: an import-time call outside a pure allowlist, a nested bare call, a subscript or attribute store, a `del` inside a compound statement, or a class with a base or a metaclass. A first-party import enters the digest only, not the closure. `[M]` 2 119 of 12 339 statements are effects under R2D (1 231 of them first-party imports), against 248 under the design's predicate.
    - **R3, amended.** A run records its import outcome at capture: whether its first-party import set imported, and the warnings emitted. When any byte of that set changes, the whole set is re-imported in one subprocess and the outcomes are compared; the comparison is cached by the set's bytes digest and the environment. qa proposed `-W error` instead; that is refused, because a warning already emitted at capture would make the check fail forever. Re-importing only the changed modules is refused too: qa measured a row it misses, where the module that breaks did not itself change. This gives the memo's validation one subprocess per distinct import-set state, so "validation imports nothing" becomes "validation imports nothing in-process".
    - **R4.** E′ is stored whole, once per manifest, and validation re-selects with it. If its size is too large, it is restricted to the first-party lexicon at capture, and that restriction becomes a counted boundary.
    - **R5.** A selected class header selects its body's non-def statements. 2 802 of the 12 339 statements are of this kind.
    - **R6.** The identifier tokens of every string constant join E, recursing through tuples and frozensets.
    - **`del` targets count as bindings.**

**The build** (2026-11-02). The memo part is written by the main agent: `orpheus/numerics/_import_time.py` (the statements, the effect predicate, the selection, the closure), and in the memo `StatementsPin` and `ImportOutcomePin`, the content-addressed sets under `.cache/sets`, and the import-outcome instrument (`python -B -c`, kept under `.cache/import_outcomes`). The test-architect's gates (`tests/gates/numerics/test_traced_memo_statement_pins.py`, 86 rows from its spec and qa's 42; `tests/gates/tools/test_cone_store_statement_grain.py`, 9 rows) pass on the built tree: `[M]` 93 passed, 1 strict xfail. The slow SG4 passed too: control 1 is exact on a clone of the real tree, selecting the 6 tests that ran the reader where B1 selected 9. The xfail was a PEP 695 alias, refused by `_def_pin`; it was fixed by recording `ast.TypeAlias` as a scope.

**Ruled by the orchestrator** (2026-11-02, from the elegance review, `scratch/verification_tree/cone_store/b2/elegance_b2.md`):

44. **The skeleton and `ModulePin` are retired, with no shim.** Every stored memo entry is stale anyway, because the memo's source changed, and an old row now reads as `Corrupt` and regenerates. The cone store's rows carry a format version (`FORMAT = 2`), and a store of another version is emptied when opened, so every test is new. The import-outcome cache is keyed by the memo's own digest as well, since that names the code that measured it, as a source index is named. Rejected, with reasons:
    - a `RetiredPin` that is always stale: it would keep a dead kind in the `Pin` union for one cycle, which the format version makes unnecessary;
    - one `kind` field on `Statement` in place of its flags: `effect`, `imports` and `defines` co-occur (an import is always an effect, and a def with a first-party decorator is both), so they are not one exclusive kind;
    - running the outcome program with `-m`: that would import `orpheus.numerics.__init__` into the measuring process and change what it measures.

**Parked downstream of this campaign (2026-11-03).** #694 makes every `LinearOperator` carry its domain and codomain; it subsumes #692. The maintainer parked it until this campaign has organised the tests by the module ontology and each module has its taxonomy, so the polish happens in a region where every operator and every test is known to work. The issue carries the measured census of unposed and half-posed operators. It also flags a likely defect: `RadialCharacteristicReconstruction` reads a `None` domain in 3123 of 3123 constructions, against its own comment. If this campaign meets that defect first, it is fixed then.
