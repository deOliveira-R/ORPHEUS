# Where verification, validation and performance live — a living plan

Opened 2026-09-22 from the R19 close-out (the agent-definitions campaign, `.claude/plans/agent_definitions.md`). Status: discussion; nothing is built. The ontology is being searched (`plan-authoring` §0), so this file is refined in place.

## The instruction (the user, verbatim)

*"reorganizing tests/ into verification, validation and benchmarking would not exactly be a major effort, because pretty much everything in tests is currently verification (and foundation), so it can keep the structure. validation doesn't currently exist, and this is the first performance measurement instrument we're keeping, I think. One of the major differences between verification, validation and benchmarking is the testing regimen. Verification for example can [probably] be described by a DAG. Validation, besides depending on a successful verification, is independent of each other, so one validation doesn't quite stack and support each other (it's possible to kind of organize in single-effect and integral effect validation, but it's not necessarily. We might not find a hierarchical structure for validation. Specially since ORPHEUS is a lattice code mainly. Benchmarking is again not quite hierarchical. Also, once verification is done, you can run validation and benchmarking in parallel. I'm not quite disagreeing with the performance at a top-level proposal, but we have to at least think where the validation would go and how to organize the future, not just the present."*

## Vocabulary (one meaning per word)

- **Verification** (L0–L2) and **foundation**: what `tests/` holds today.
- **Validation** (L3): a comparison against experiment, such as ICSBEP and IRPhE critical and lattice experiments.
- **Code-to-code comparison** (L4): a comparison against another code's answer. The lattice world calls its problems "benchmarks" (C5G7, the VERA progression problems), but by `vv-principles` these are L4. L4 produces no correctness information, and each L4 claim must name its L0–L3 backing.
- **Performance**: wall time, memory, and cost scaling. The user's "benchmarking" in the instruction means this.

"Benchmark" currently has three meanings: a performance run, an L4 comparison, and a published lattice problem. The layout below gives each its own word.

## The organising principle `[R]`

A directory boundary should be a *regimen* boundary: what differs in how a suite runs, what its dependencies look like, and what its result means. Per suite:

| | verification + foundation | validation | code-to-code (L4) | performance |
|---|---|---|---|---|
| dependence among cases | a DAG (`rests_on`, #358) | none; each is independent | none | none |
| rests on | lower rungs | the verification of every capability it exercises | the same | the same |
| a case's result | pass or fail | C/E with experimental and computational uncertainty | C/C′ | time, memory, count |
| meaning of a bad result | a defect | a model-form or data question, judged across the suite (trends, bias) | informational | a regression or a tradeoff |
| what to keep | the verdict | the history of C/E, per case and per commit | the history | the history per commit |
| cadence | every commit (CI gates) | on demand, release | on demand | on demand, advisory |
| order | first | after verification, in parallel with the other two | the same | the same |
| grouping | the package layers (the module tree) | experiment family; single-effect versus integral-effect as a case attribute | the problem suite | the workload |

Two observations follow:
- **Verification is the only suite whose cases stack.** Validation, code-to-code and performance share one shape: a flat set of cases, each resting on the verification of the capabilities it exercises, each producing a measurement with a history rather than a verdict.
- **The edge from a flat case into verification is capability-level, not test-level.** A validation case rests on "SN 2-D Cartesian eigenvalue, multigroup, reflective", which is a subgraph, not on one test. The project already has a capability axis, `cap(name)` (the SN capability tiers in `pyproject.toml`), so the flat suites would declare `cap`-level supports while verification declares test-level `rests_on`. That is #358's question too.

## Layouts considered

- **L-A. Siblings.** `tests/` stays as the pass/fail suite, unchanged. `validation/`, `performance/` and, when needed, a code-to-code suite are top-level siblings, each with its own runner and cadence, importing shared builders from `tests/`, for example `tests/_harness`.
  - For: the directory boundary equals the pass/fail boundary, which is what `tests/` means by Python convention; bare `pytest` never collects a measurement suite; no test id changes, so the `rests_on` ids and the many page citations of `tests/...` paths stay valid.
  - Against: three or four top-level directories.
- **L-B. One root with regimen subdirectories.** `tests/verification/` (today's tree moved), `tests/validation/`, `tests/performance/`.
  - For: one discoverable root; the regimen is visible in the path.
  - Against: `pytest tests` collects all three unless `testpaths` or markers exclude two; every test id and path citation changes once; foundation tests sit under "verification" though they are not verification.
- **L-C. A V&V parent.** `vv/verification/`, `vv/validation/`, `vv/performance/`. The same as L-B with a neutral root name; the same costs.
- **L-D. Everything in `tests/`, separated only by markers.** Rejected: the marker axis already carries the V&V level, and the runner, cadence and meaning of a result differ, which a marker does not change.

**Leaning `[R]`: L-A.** Its boundary is the one that carries meaning (gated versus measured), and it leaves the verification DAG's ids stable. The regimen table becomes the page that says where a new case goes.

## What each suite would need (sizing, not a criterion)

- **performance/**: airspeed velocity (asv): `time_*`, `peakmem_*` and `track_*` (any number, such as an iteration count or an error), results per commit, `asv continuous` between two commits. First cases: the ScanMarch-versus-window measurement (`git show f36572c8^:derivations/diagnostics/diag_s69_scanmarch_vs_window_bench.py`) and the wall-clock leg of `tests/sn/architecture/test_composition_cost.py`. Each case names the capability it rests on and records its accuracy beside its time.
- **validation/**: nothing yet. L3 is sequenced after L1 maturity (`vv-principles`). When the first case lands, it carries:
  - an experiment identifier (the ICSBEP or IRPhE id);
  - the specification;
  - a reference to the measured data and its uncertainty (the data may be licensed: the same question as `scratch/literature/`);
  - the model;
  - the capabilities it exercises;
  - a C/E output with uncertainty.
  
  Its history store and its runner are open (asv's `track_*` could hold a C/E history; uncertainty bands may need more).
- **CI**: the `gates` workflow stays verification. Validation and performance would run in a separate workflow whose jobs depend on verification's success (`needs:`) and run in parallel.

## Open questions for the user

1. The layout: L-A (siblings), or L-B/L-C (one root)?
2. Is code-to-code (L4) its own suite, or a kind of case beside validation? The lattice benchmark problems make it likely that L4 cases come before real validation does.
3. Create `validation/` now with a README carrying the case contract, or only when the first case lands (Pattern 6: no abstraction before an instance)?
4. For performance: asv, or pytest-benchmark (one runner, weaker history)?

## Rulings and answers so far

- **R1** (the user, 2026-09-22): no code-to-code (L4) suite. Keeping other codes' results in the repository costs maintenance for negligible value. L4 is named in the principles page only so that it is not confused with verification. Open question 2 is closed.
- **Discovery** `[M]`: pytest collects what `testpaths` names (`pyproject.toml`: `testpaths = ["tests"]`) or the paths passed on the command line, and within them the files matching `python_files` (default `test_*.py`, `*_test.py`). The name `tests/` is convention only.
- **Runners** `[R]`, answered in the session:
  - asv is a runner with its own discovery (`time_*`, `peakmem_*`, `track_*`), not an instrument of a pytest run. It builds the project at chosen commits, times in separate processes with repeats, and stores per-commit results.
  - pytest-benchmark is the in-pytest alternative: a fixture that times a callable inside a test.
  - Validation is a pipeline (case registry, parallel runs, C/E with uncertainty, a history, suite-level statistics, a report) rather than a pass/fail suite.
- **R2** (the user, 2026-09-22): one root. Everything that tests the code, of whatever kind, lives under `tests/`; the root stays clean; siblings named `validation/` beside `tests/` would read as if validation were not a test. The name of the pass/fail subtree is under discussion (next item).
- **R3**: `validation/` is created when its first case lands; an empty folder is a place to search and be surprised.
- **R4**: asv for performance.
- **R5**: the face-transmission derivation becomes an algebra of record in `orpheus/derivations/discrete/sn/`, improved and extended, possibly comparing step and linear discontinuous (LD) with diamond.

**Measured for the subtree's name** (`scratch/_r19/layout_census.py`, AST over the 546 tracked `tests/**/test_*.py`, control: a known foundation file): 325 files carry only `foundation`, 154 only a level marker (`l0`–`l3`), **66 carry both**, 1 neither; the tree's top level is 13 module directories plus 9 root test files that are the repository's own integrity gates (layer imports, harness drift, docstring cross-references, the pyright ratchet, the error-catalogue reconciliation, the elegance-debt ledger), neither verification nor a foundation of the physics code. Sizing of a move: 325 `tests.*` import statements inside `tests/`; 491 tracked files outside `tests/` cite a `tests/` path.

**Proposal `[R]`: `tests/gates/`, `tests/performance/`, and `tests/validation/` later.** The subtree is named by its regimen, not by its kind: the kind of a gate (math verification, foundation, repository integrity) stays on its markers, where it is per test; a `verification/` or `foundation/` directory would be a second definition of the marker, and the 66 mixed files show the kind is not a property of a file. "Gate" is already the project's word for a pass/fail check (the CI workflow is `gates`). Shared builders (`tests/_harness/`) and the root `conftest.py` stay at the root, available to every regimen.

Arguments for one root that the user did not list:
- A retirement audit greps `tests/` for a symbol's consumers; a performance or validation case outside it is a surface that audit misses, the same shape as `derivations/`, which was invisible to every instrument for four months (#347).
- Nexus indexes `orpheus/`, `docs/` and `tests/` (`.nexus/config.toml`): cases under `tests/` enter the graph with no configuration.
- Context cost of the longer path: one segment per hit, negligible. No LLM reason favours root-level folders: agents navigate by listing and grep, and a regimen table in `tests/README.md` orients them better than a folder name.
