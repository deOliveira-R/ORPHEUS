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

## ⏸ Start here

Dispatch the workflow once it is on `main` (`gh workflow run test-durations`), read the artifact, and turn the table into the ranked worklist here. Step 2 (R1) starts from it.
