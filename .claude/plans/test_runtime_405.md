# A bounded routine gate, and a slow tier that runs on a schedule (#405)

Opened 2026-09-23. Status: open, not started. The user's instruction, verbatim: *"do we know what in the tests takes so long to run? cause we need to deal with this in a smart way."* Then, on the proposal below: *"yes. Do that"* (2026-09-23). This is a living plan (`plan-authoring` §0). The means below are hypotheses until measured; the first step is the measurement.

## The goal, in the project's terms

"`main` is always green" is a claim about gates someone actually runs. Two properties are wanted:

1. **The routine gate is bounded and observable.** A session can run the gate for a change and read its progress, and a death is loud, never silent.
2. **Every gate in the tree runs on a known schedule, including the slow tier.** Today the slow tier has no last-run date; that is silent-coverage risk on exactly the heaviest verification gates.

Timing is never a gate (ruling R6 of `.claude/plans/vv_suite_layout.md`). This campaign is about the cost of running the gates, not about gating on time.

## What is known

- **#405's measurements `[M]` 2026-08-24**, host `.venv`, serial, `python -O -m pytest -p no:randomly`. SHELF-LIFE: these are a month old, and the tree has changed since (the move to `tests/gates/`, R19, new gates). Re-measure before relying on them.
  - The un-deselected tree (no `-m` filter): sn 3352 tests in 2 h 23 min; mc 57 tests in 40 min; derivations interrupted at about 145 tests after 3 h 19 min; cp: test #66 ran over 52 min, about 25% of it in the garbage collector.
  - The marker split by `--collect-only`, in the same run:

    | tree | not slow | slow |
    |---|---|---|
    | sn | 3236 | 116 |
    | derivations | 1661 | 67 |
    | mc | 41 | 16 |
    | cp | 141 | 13 |

  - The named giants are all slow-marked: `test_l1_standoff_slab_cylinder.py`, and `test_phase_c_crosscheck.py` phases d/e at 777–955 s each.
  - The routine gate, `-m "not slow"`: last measured at 52–90 min (memory `reference_test_execution_env`, dated before the move to `tests/gates/`).
- **`[M]` 2026-09-23, this session:**
  - `pyproject.toml` applies no default `-m` filter (`addopts = "--import-mode=importlib"` only), so a bare `pytest <tree>` runs the slow tier. That is how this session's wide run ended up inside the slow tier.
  - A redirected run block-buffers its dots, so its log lags reality. ~~This session's `-m "not slow"` run died with no summary~~ [REFUTED 2026-09-23] It completed: exit 0, "4483 passed, 14 skipped, 67 deselected, 18 xfailed in 2443.79s (0:40:43)". My liveness check piped `ps` through `head` and cut the live pytest process off the listing (VIEWPORT). The lesson stands: from a buffered log alone, "stuck", "dead" and "grinding" cannot be told apart.
- **First worklist data `[M]` 2026-09-23**, host `.venv`, serial, `python -O -m pytest -m "not slow" --durations=30` over `tests/gates/transport`, `tests/gates/derivations`, `tests/gates/sn/operators` and `tests/gates/test_layer_imports.py`, in one process: 40 min 43 s in total. Per tree in separate processes, in the same session: transport 48 s (1128 passed), sn/operators 67 s (1351 passed), layer imports 5 s (370 passed). So derivations' not-slow gates take about 38 min of the 40.
  - The slowest 30 calls, 23–66 s each, are all 30 in `tests/gates/derivations/`, and 26 of the 30 are in `test_peierls_*` files (counted from the log's durations block). Examples: kernel row-sum convergence with panels, 50 s; Nyström self-convergence, 42 s; the rank-2 closure scaling rows, 26–29 s each; the specular multibounce overshoot at high N, 50 s. The slowest is `test_dsa_rules.py::TestFourStepDerivation::test_main_theorem_interior_row_is_larsen_27`, a symbolic derivation, at 66 s.
  - Reading for R1 `[R]`: most of the not-slow cost is semi-analytical references recomputed per run. That is exactly the category R1 wants served from a generator-keyed cache.
  - Log: scratchpad `405_first_durations_2026-09-23.log`; it is a scratch file, so the numbers above are the record.
- **Installed:** `pytest-xdist` 3.7.0 (pinned `<3.8`: 3.8.0 deadlocks on Python 3.14.3); `pytest-timeout` 2.4.0; `coverage` 7.14.1; `asv` 0.6.6. Not installed: `pytest-split`, `py-spy`. The `pyproject.toml` comment on xdist says a whole-tree `tests/gates/sn -n auto` uses about 21 GB of aggregate RSS; per-tier `-n` is the working pattern.
- **CI:** `.github/workflows/gates.yml` runs the cheap half on `ubuntu-latest`. It already caches the HDF5 cross-section store, keyed on the LFS tapes' object ids, which is the expensive part of any CI test job's setup.
- **Related issues:** #211 (reference selection, caching and parallelism in the SN suite); #377 (`RigidMotion.determinant` called 55.7M times in one test directory); #404 (the one pre-existing red the #405 run surfaced).

## The proposal (the user accepted it in outline, 2026-09-23; each step is a hypothesis until measured)

1. **Measure the whole tree once, sharded on CI.** A `workflow_dispatch` (and later scheduled) workflow:
   - a matrix of shards, one per tree, with sn and derivations split further by file;
   - every shard runs `python -O -m pytest --durations=0 -p no:randomly` with `PYTHONUNBUFFERED=1`, with and without the slow tier;
   - it uploads the JUnit XML (which carries per-test time) as an artifact;
   - one aggregation step turns the XML into a per-test duration table, committed as a data file, for example `tests/performance/test_durations.json`.

   Open `[R]`, from GitHub's documented limits, not checked for this repository: the per-job ceiling is 6 h on GitHub-hosted runners, and a public repository's `ubuntu-latest` runner has 4 vCPUs and 16 GB. A shard that exceeds 6 h is itself a finding; shard by file where one tree is too big.
2. **Bring the giants down at the source**, test by test, from the table: cache deterministic references on disk, keyed by a hash of their inputs and of the code that produces them; use coarser refinement rungs in the routine tier, keeping the full ladder in the slow tier; fix the cp allocation loop.
3. **Run the slow tier on a schedule**, sharded by the measured durations (nightly or weekly), so it has a last-run date.
4. **Select tests per change:** Nexus `retest` with a coverage capture. Its static cone has 12–15% recall `[M]` 2026-08-18 (the `retest` tool's own note), so it needs the capture.
5. **Retest xdist** per tree on Python 3.14, and adopt it where it is stable.
6. **The observable runner as the documented pattern** for any local run longer than minutes: per-tree invocations, unbuffered output, START/EXIT banners, durations, detached. `[M]` this session's commit gate used it: scratchpad `ft_gate_driver.sh`.

## Rulings

- **R1 (the user, 2026-09-23), the goal:** *"We will first of all fight hard to eliminate slow verification. Maybe we reach the point where only regenerating semi-analytical solution is a slow process, which we can cache, and then they can be run only when the generator code changes. If we can't reach that point, we will think about how often it runs."* So step 2, bringing the giants down at the source, is the campaign's centre. Its target state: a gate is fast, and a semi-analytical reference is regenerated only when its generator changes, served from a cache keyed on the generator's code and inputs. The schedule of a residual slow tier (step 3) is decided only after that.
- **R2 (the user, 2026-09-23), the duration table:** *"it can be commited if there is a use for it, but I suppose CI instances have a certain configuration. It only makes sense to compare the run on reproducible settings. You can make a proposal on what to do with the table. Why would we commit it?"* The proposal, not yet ruled:
  - the raw table is a CI artifact of each run, stamped with the runner configuration (image, CPU count, Python and package versions), and two tables are compared only when their stamps match;
  - the ranked worklist of hogs, each with its measured cost, the run it came from and the runner stamp, lives in this plan;
  - nothing is committed until a consumer exists; the one candidate is shard balancing, which can read the latest artifact at run time.
- **R3 (the user, 2026-09-23), a default `-m "not slow"`:** deferred until R1's outcome: *"Let's see what happens with 1. That has an influence on 3."*

## Step 1, built — 2026-09-23 (branch `feature/test-durations-405`)

- `tools/test_durations.py` with `shard` (round-robin over the tracked `tests/gates/**/test_*.py`; `[M]` 548 files into 40 shards, every file exactly once) and `aggregate` (JUnit XML to a table stamped with the runner, plus a Markdown summary of the slowest tests and the hours per tree).
- `tests/gates/tools/test_test_durations.py`, 6 tests. Node ids are checked against pytest's `--collect-only`; outcomes, including a pytest-timeout, against a synthetic file whose every outcome is known; dealing is checked to be a partition. `[M]` A real repository file (`tests/gates/tools/test_write_guards.py`, 22 tests) reads back identically to `--collect-only`. Two mutations (timeout read as failure; class dropped from the node id) each turn 2 of the 6 red, and the control stays green.
- `.github/workflows/test-durations.yml`: `workflow_dispatch` with inputs `shards` (40), `marker` (empty, which selects every test: `[M]` `-m ""` collects 58 of 58 in `tests/gates/mc`, against 42 with `-m "not slow"`) and `per_test_timeout` (3000 s); at most 20 shards at once, 355 min per shard; the store cache as in `gates.yml`.
- Documented in `docs/development/git_workflow.rst`, "Continuous integration".

## Step 1, measured — run 35940553034 `[M]` 2026-09-23/24

The runner stamp, as recorded in the artifact `test-durations`: ubuntu24, image 20260907.300.1, 4 CPUs, x86_64, Python 3.14.7, numpy 2.5.3, scipy 1.18.1, sympy 1.14.0, commit `62020ba5`. The run: 40 shards, `-m ""` (everything, slow tier included), 3000 s per-test timeout. The run's own summary page has the pre-`tree_of` per-tree bug (fixed at `dc5b2ad2`); the numbers below are recomputed from the artifact's table with the fixed tool.

- **12 491 test cases in 549 files take 15.02 h serial.** Outcomes: 12 262 passed, 87 skipped, 56 failed, 81 errors, 5 timeouts.
- **Concentration:** the 10 slowest cases are 44.2% of the total; the 30 slowest, 72.4%; the 100 slowest, 88.4%; the 300 slowest, 96.3%.
- **The ranked worklist, by file** (minutes of serial time on this runner; the top 8 files are about 640 of the 901 minutes):

  | min | file | what it computes `[R]` from the file names; read each before acting |
  |---|---|---|
  | 188.9 | `tests/gates/derivations/test_peierls_specular_bc.py` | specular-BC Peierls references, multibounce; 1 timeout (`test_specular_heterogeneous_2G2R_converges[slab]`) |
  | 103.3 | `tests/gates/sn/verification/analytical/test_l1_standoff_slab_cylinder.py` | cylinder refinement ladder against the trajectory resolvent |
  | 100.1 | `tests/gates/derivations/test_continuous_registry_lazy.py` | 2 timeouts: a lazy registry fetch builds a Peierls reference |
  | 53.7 | `tests/gates/derivations/test_peierls_rank2_bc.py` | 1 timeout: the rank-2 refinement monotonicity ladder |
  | 50.9 | `tests/gates/cp/test_peierls_rank_n_protocol.py` | the rank-N protocol (F.4 subprocess worker) |
  | 50.0 | `tests/gates/cp/test_peierls_flux.py` | 1 timeout: 2-group 2-region flux convergence |
  | 49.8 | `tests/gates/derivations/test_peierls_multigroup.py` | multigroup parity against the unified path |
  | 45.7 | `tests/gates/sn/verification/analytical/test_phase_c_crosscheck.py` | trajectory-resolvent cross-checks, phases d/e |
  | 29.8 | `tests/gates/mc/test_convergence.py` | Monte Carlo σ scaling with √N |
  | 23.9 | `tests/gates/sn/sweep/curvilinear/test_unified_matvec_cylinder.py` | unified cylinder matvec against the trajectory resolvent |
  | 22.6 | `tests/gates/derivations/test_peierls_reference.py` | Peierls kernel row sums and element-wise references |
  | 17.9 | `tests/gates/derivations/test_peierls_greens_function_cylinder_mr.py` | Green's function, multi-region cylinder |

  Everything else is under 13 minutes per file. Of the 12 files above, 11 compute Peierls, trajectory-resolvent or Green's-function semi-analytical references (all but the Monte Carlo file): R1's cache target.
- **Two infrastructure findings in the failures**, not defects of the code under test, apart from `phase_e`, which is #404:
  1. **83 of the 137 (81 errors and 2 failures) are the workflow's.** The tests read the raw GENDF tapes, and the workflow pulls them from LFS only when the HDF5 store cache misses. The cache hit, so the tapes were LFS pointer files ("could not convert string to float: 'oid sha256'"; modules `test_ingest_ledger` 43, `test_n2n_yield_convention` 24, `test_hdf5_store` 6, and others). Fix: cache the tapes themselves, keyed on their pointer files, so LFS is pulled once per content change. Pulling on every shard would exhaust the LFS bandwidth quota `[R]`.
  2. **About 50 are bit-identity and snapshot gates that pin values captured on the development Mac** (arm64, Accelerate): `test_byte_stability` (homogeneous k_inf, 1 ULP), `test_walk_matvec_baselines` (1 ULP bound, up to 768 read), `test_bc_extraction_*` (5 ULP bound, up to 128 read), `test_streaming_operator` T4b snapshots (256 ULP bound, up to 512 read), `test_affine_carve_baseline`, `test_quadrature_fold` (`==` against 4π), and others. On the x86 runner they differ by 1 to 768 ULP. Bit identity is a platform property as well as an implementation one, so these gates are red on any machine but one. Filed as an issue (below).
- **The runner concurrency:** `max-parallel: 20` took every job slot the account has `[R]` (20 on the free plan), so the `gates` run for `62020ba5` queued behind the shards. Set it to 16.

## ⏸ Start here

Step 1 is measured (above). Next: (a) [REMEDIED 2026-09-24] the workflow fixes (the tapes cached once in the plan job, restore-only in the shards; `max-parallel: 16`); #504 filed for the platform-bound bit gates; (b) step 2 (R1) starts from the worklist: design the generator-keyed reference cache with the user. It is ontology work (what a reference is, what its key is), so it opens as a discussion in this plan, not as code.
