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

## ⏸ COMPACTION POINT — 2026-09-24

Step 1 (the measurement) is done and its fixes are merged (`62020ba5`, `dc5b2ad2`, `e746655f`, `668703f2`; `main` green at `668703f2`). Step 2, the heart of R1, is its own living plan: `.claude/plans/reference_cache.md`. Resume there, with its eight questions, after regrounding. Still open here:
- ruling R2's proposal, the table as a stamped artifact with the worklist kept in this plan (not yet ruled);
- R3, the default `-m "not slow"`, deferred until R1's outcome;
- #504, the platform-bound bit gates.

The durations workflow can be re-run at any time with `gh workflow run test-durations`; re-run it after step 2 lands, to measure the gain on the same runner stamp.

## The two files the memo could not shorten (2026-10-04, after #405 P3)

The P3 timing protocol (`reference_cache.md`, "P3 close-out") left two files slow when warm:
- `test_l1_standoff_slab_cylinder.py`: 1049 s warm, a speed-up of 1.8×;
- `test_trajectory_resolvent_reference.py`: 162 s warm, a speed-up of 2.5×.

The user asked what makes them hard to shorten. Two numerics-investigators measured them `[M]`; their reports and probes are in `scratch/reference_architecture/p3/perf_l1/` and `perf_traj/`.

**`l1_standoff`:**
- **The cause is production:** the Krylov inner solve runs UNPRECONDITIONED (`orpheus/sn/solver.py:950` passes an identity). GMRES iterations equal the degrees of freedom, and the time grows as O(n³).
- **The fix is issue #200, preconditioning with the sweep:**
  - slab at n_per = 40: 18.2 s → 0.2 s, with 5705 → 81 inner iterations;
  - cylinder at nx = 40: 46.4 s → 1.1 s;
  - k moves by 1.6e-10 relative.
- **Test design:**
  - the slab solves run twice across rows, because the slab solves have no cache;
  - the n_per = 160 rows stay green when a first-order step scheme replaces the matvec (error 1.06e-5 against a 2e-5 tolerance); the {10, 20, 40} ladder with p ≥ 1.8 reddens it;
  - `catches("ERR-025")` on `krylov_via_unified_vs_case` is decayed: the Krylov k does not move under the ERR-025 mutation;
  - three docstrings are stale: "O(h)" (measured O(h²)), "monkey-patched", and the tolerance comment.
- **Secondary production item:** the 1-D apply walk `_loop_walk` → `visit` is a per-cell Python loop (60 % of unpreconditioned Krylov time); the DD outflow is a prefix scan.

**`trajectory_resolvent_reference`:**
- **The cause is production:** `MultiRegionCylinderChordOracle.apply_operator` (`chord_oracle.py:941`) rebuilds the in-plane segments and their spline values for every axial cosine, about 3.3M scalar `CubicSpline` calls per reading.
- **A hoisted prototype** (`perf_traj/factored.py`) agrees to 9e-16 relative:
  - solve 20 → 5.2 s;
  - reading 53 → 4.2 s;
  - brute 7.7 → 0.19 s.
- **Test design:**
  - the in-process spy rows check routing and run at the gates' resolution; at a minimal quadrature (4, 2, 4) they take 71 → 4 s and 20 → 3.3 s, and all three mutations still redden them;
  - the fine-rule row's brute integral passes by NODE PLACEMENT: unsplit GL over the tangency kinks converges erratically in the azimuth count (1.2e-5, 1.25e-4, 3.7e-6 at 512, 640, 768 points); a reliably converging rule is (24, ≥ 3072).

**Ruled (the user, 2026-10-04): "Gates, then #200, then oracle".**
1. Repair the two files' gates (the test-architect), on branch `test/slow-reference-gates`.
2. #200, the sweep-preconditioned Krylov, as a surgical SN carve.
3. The chord-oracle hoist, plus an issue for vectorising the 1-D apply walk.

### Step 1 landed (2026-10-04): the two files' gates repaired

The test-architect's evidence (each mutation table and the timings) is in `scratch/reference_architecture/p3/gates_repair/README.md`.

**`l1_standoff`: 1821 → 100 s.**
- The slab n_per = 160 rows are replaced by the {10, 20, 40} ladder: an order row (p ≥ 1.8) and the Case comparison at 40.
  - The step-scheme mutation reddens both Krylov rows: p = 1.25 and 1.16, and 4.66e-5.
  - The Σ_t × (1 + 1e-3) mutation reddens both Krylov rows: p ≈ −0.05, and 4.48e-4.
- One cached `_slab_k(path, n_per)` per solve, so 6 solves and none repeated.
- `catches("ERR-025")` moved to the new `test_slab_l1_sweep_vs_case`, which reddens at 7.23e-2.
- The cylinder twin row drops nx = 80: its gap is mesh-independent, 4–5e-11.
- Four stale docstrings fixed, and the pre-existing optional-`k_eff` pyright errors fixed.

**`trajectory_resolvent_reference`: 564 → 27 s.**
- The spy rows run at a minimal quadrature. Each of the 4 mutations reddens its row.
- The cylinder fine-rule row's brute integral is split at the INTERFACE tangencies, computed in the test by arcsin from the radii, with θ = 32 and 64 points per piece. The row is renamed `…_against_an_independently_split_fine_angular_rule`.
  - Why: the unsplit rule converged erratically at every affordable azimuth count (`[REFUTED 2026-10-04]` the investigator's "(24, ≥3072) converges reliably"; at 3584 it sat at 80 % of the band).
  - The gap to the reading is at most 3.2e-6 against the 1e-5 band.
  - A dropped tangency split in the reading reddens the row at 2.4e-5; a 2 % misplaced one reddens it at 4.7e-5.

**Open from this step:**
- **Why does ERR-025's coefficient mutation in the sweep not trip `ConvergenceClaimError`, when Σ_t × (1 + 1e-3) does?** `[R]` The residual re-check runs through the matvec, which reads the same DD coefficient as the sweep. The check sees whether the sweep solved ITS equation, not whether it solved the right one (X4: a shared upstream). This is a property of the guard, not a defect, and is recorded here.
- **The sphere brute** (unsplit, 2000 points) has about 4× headroom; review it with the oracle hoist.
- **Re-time the Krylov order row** (25 s) after #200.

### Step 2 (#200) built; reviews in (2026-10-04)

**Branch `fix/krylov-sweep-preconditioner`**, stacked on `test/slow-reference-gates`. NEITHER branch is merged to `main` (`main` is `6ffd1960`).
- `dcc73bb1`: step 1, the two files' gates.
- `ca7d9c21`: step 2, the code change. `_within_group_krylov` (`orpheus/sn/solver.py`) preconditions GMRES with `seeded_inverse(LC)`, plus the optional DSA corrector. The gates are in `tests/gates/sn/solve/test_krylov_sweep_preconditioner.py`:
  - the route;
  - linearity on full and boundary-only residuals;
  - the boundary round trip `A M q_b = q_b`, which is the row that catches a trace-dropping sweep (such a sweep is linear, so the linearity rows cannot see it);
  - the fixed point;
  - the rate, flat across 10/20/40 cells, where the identity grows 3.7× and 4.0×.
  
  Also in this commit: the k_inf rows run the production preconditioner; `test_krylov_restart_signature` is re-tightened to 1e-9; #200's σ_r-fold stack is split to **#575** and the removal-form `xfail` is re-keyed to it. The commit carries `Closes #200`.
- `58bb3000`: the archivist's docs, including a new section in `acceleration.rst` (label `sn-krylov-sweep-preconditioner`) and corrected #200 claims across 8 pages and ERR-050/052/071.
- `6ccc531b`: the regenerated matrix and error index.
- `[M]` The full non-slow suite at `ca7d9c21`: **15 248 passed, 0 failed**. The l1 file: 100 → 19.7 s.

**The user's ruling (2026-10-04), on the elegance finding:** "We will do in this branch, but we will compact context first, so that you can do it with a cleaner context." So the work below lands on THIS branch before it merges.

**To do on this branch, in order:**

1. **The preconditioner is an OPERATOR** (elegance S1 and S2; `scratch/reference_architecture/p3/elegance200/review.md`).

   **Ruled (the user, 2026-10-04), superseding the sub-items below where they differ:**
   - **The contract is a STATED choice, the identity allowed.** `preconditioner: LinearOperator` is keyword-only with NO default; plain Krylov is `preconditioner=IdentityOperator()`, written at the call site. No warning and no fallback.
   - Why: a preconditioner is not mathematically required, and the sweep is the right one only for SN. Diffusion's `A.inverse()` is the whole solve, so an `A⁻¹` default is SN knowledge inside a numerics primitive. A silent identity default would hand the next model's Krylov path the #200 failure: iterations growing with the mesh, with nothing on screen. A warning could not have caught #200 either, since production passed an explicit identity. Unpreconditioned Krylov is legitimate in three cases: small well-conditioned systems; an operator that is already the preconditioned one (`M⁻¹A`, the common production transport form); and studies. So the identity is allowed but must be named.
   - The SN guard against regressing to the identity is the gate file's route row and its mesh-flat rate row, not a warning.
   - **`P = (I + C) M⁻¹`: Krylov only, plus a witness.** `_within_group_krylov` spells `(I + C) @ LC.inverse()` (`[M]` bitwise equal to the closure on 2 of 2 slabs, `scratch/reference_architecture/p3/operator_contract/probe_compose.py`). SI keeps applying `C` to the increment, byte-identical. A gate asserts SI's corrected step equals `ψ + P r`. The full merge waits for a second corrector (2-D DSA).
   - **Census correction** `[M]` 2026-10-04, AST over `orpheus/` and `tests/`: **17 constructions** (1 production, 16 tests: 7 default, 7 lambda, 3 named), not 106. The 106 was the elegance spy's count of runtime constructions over parametrised tests.
   - `KrylovAcceleration`'s `preconditioner` becomes a `LinearOperator[V] | None`, a public attribute (`orpheus/numerics/iteration.py`: the `Preconditioner = Callable[[ndarray], ndarray]` alias at :223 is a false type, since the callable receives typed fields; the default is at :982-987).
   - `solve` wraps `self.preconditioner.apply`.
   - `_within_group_krylov` passes `LC.inverse()` (or the seeded form the default uses) with no corrector, and `(I + C) @ LC.inverse()` with one. Today the closure at `solver.py:956-959` is bitwise the `KrylovAcceleration` default (`[M]` 4 of 4 geometries), so there are two spellings.
   - Elegance measured `inv + C @ inv` as an `OperatorSum`, bitwise equal to the closure on 2 slabs; the DSA corrector refuses curvilinear meshes (`NotImplementedError`, its scope edge).
   - Check what `seeded_inverse(A)` returns on a carrying mesh (`_SeededExactApply` / `CoupledSubstitutionOperator`) before choosing between `LC.inverse()` and `seeded_inverse(LC)`.
   - Migrate the test constructions that pass `lambda q: q` and other callables (`[M]` 106 constructions pass their own preconditioner across the three Krylov files; 8 pass an identity lambda, which becomes `IdentityOperator(space)`).
   - Remove the silent no-preconditioner fallback when `A` is not invertible (`iteration.py:986-987`): the identity becomes explicit.
   - Elegance C5: source iteration with a corrector is Richardson with `P = (I + C) M⁻¹`. `iteration.py:793` applies `I + C` a second time; name `P` once so that SI and Krylov consume it.
   - Elegance S3: `test_default_sweep_preconditioner_recovers_kinf_on_slab` and `test_production_preconditioner_recovers_kinf[slab]` now build the same object. Merge them, keeping the ERR-050 argument, and replace the string `preconditioner_kind` with an operator argument.
   - Elegance C7, nits in the new test file:
     - the one-line aliases `_system_a` / `_random_state`;
     - the `solve_sn` keywords written twice;
     - private builders imported from other test modules, which belong in a shared helper.
2. **qa HIGH: the ERR-053 catchers have decayed.** With the sweep, GMRES needs 15–25 iterations, so a restart clamp of 50 never bites.
   - Re-dropping `restart=min(50, n_dof)` leaves green: the 6 `test_krylov_kinf_independent_of_mesh_refinement` rows, `test_solve_sn_si_vs_krylov_consistency_homogeneous_sphere`, and `test_krylov_restart_covers_augmented_composite[5,10,20]`.
   - Only the 3 `test_g_d3_3_*` restart-spy site gates go red.
   - Re-home `catches("ERR-053")` onto the site gates, or add a value gate that needs more than 50 GMRES iterations (the identity arm, or a harder problem).
   - Fix the restart-signature comment: the solve "stops on the PRECONDITIONED residual" is false, since scipy accepts on the TRUE residual.
   - Evidence: `scratch/reference_architecture/p3/qa200/m1_restart50_full.log`.
3. **qa MEDIUM-HIGH: the inner `IterationRecord` misreads convergence under a preconditioner.**
   - scipy's `pr_norm` callback reports ‖M r‖/‖b‖, which is relative only when M = I. The record judges it against `tol`.
   - Reproducer `qa200/q1b.py`: a thin slab, 2G with upscatter. 4 of 7 inner records read not-converged; there is a false "hit max_inner" warning; `fully_converged=False`. All of this while scipy returned info = 0 and the true relative residual is ≤ 9.5e-9.
   - It also blinds `_check_convergence_claim`. The DSA posture already had the defect; #200 spreads it to every Krylov solve.
   - **Needs a ruling:** normalise by ‖M b‖, or record the TRUE residual (recommended: the record should hold what it claims). An ERR entry follows.
   - **Ruled (the user, 2026-10-04): the true residual AND the scaled trajectory.** scipy's inner loop stops on `‖M r‖ ≤ rtol·‖M b‖`, while its callback reports `‖M r‖/‖b‖` and it ACCEPTS on the true `‖b − Ax‖ ≤ rtol·‖b‖`, tightening and continuing when the two disagree (`scipy/sparse/linalg/_isolve/iterative.py`, `gmres`). The record gets two criteria:
     - `residual`: the final true relative residual, one matvec per solve; the record's verdict then equals scipy's acceptance;
     - `pr_residual`: the per-step trajectory divided by `‖M b‖`, what scipy steers on, keeping the rate diagnostics; one preconditioner apply per solve.
   - **[REFUTED 2026-10-04] the two-criterion form above, and RE-RULED the same day (the user):** `IterationRecord` requires co-indexed criteria (`orpheus/numerics/convergence.py:1037`, every trajectory the same length), and the true residual exists only at exit (scipy forms the iterate once per restart cycle, and SN runs one cycle), so a one-point `residual` criterion beside a per-step trajectory cannot be built. The true residual already has ONE home, `_check_convergence_claim` (`orpheus/sn/solver.py:819`), which re-measures `‖Aψ − q‖/‖q‖` for every claimed convergence; computing `‖b − Ax‖` inside `KrylovAcceleration` would be a second spelling of that quantity (X4).
     - **The ruling:** ONE criterion, the per-step trajectory divided by `‖M b‖` (scipy's inner stopping quantity; one preconditioner apply per solve). This removes the false not-converged reading. The false-early reading (a small preconditioned residual with a large true one, for example when scipy exhausts its budget) is caught by the claim check, which now receives an honest claim.
     - Known gap, not new: the claim check is skipped for the moment-tailed LD schemes (`_residual_is_expressible`), as it is for source iteration.
4. **qa LOW:** add a 2-D Cartesian row to the gate file's `_GEOMS`. qa ran the three p200_1 rows on 2-D, a 3-cell cylinder and GL8: 12 of 12 pass.
5. **The archivist's out-of-scope findings:**
   - drop "the Krylov sweep preconditioner of #200" from `orpheus/sn/loss_representation/assembly.py:22`: the preconditioner is `seeded_inverse(LC)` and no production module imports `assembly`;
   - `orpheus/numerics/green_operator.py:105` says "decided with #200/#284", which is now decided;
   - the docstrings of `test_krylov_curvilinear_precond_safety.py` (about :45-48 and :73) are stale;
   - the commit's "1.6e-10 relative" is really 1.55e-10 absolute, 1.22e-10 relative (cylinder 1.8e-10 and 1.5e-10). The docs publish both correctly.
   - Optional: `@pytest.mark.verifies("sn-krylov-boundary-round-trip")` on `test_p200_1_a_boundary_only_residual_round_trips`.
6. Then the full suite, qa, merge both branches to `main` with `--ff-only`, and watch CI.
7. **Step 3 of the user's ruling:** the chord-oracle hoist.
   - `MultiRegionCylinderChordOracle.apply_operator`, `chord_oracle.py:941`; prototype `scratch/reference_architecture/p3/perf_traj/factored.py`, 9e-16 relative, not bit-identical, so grep frozen-byte consumers.
   - Review the sphere brute's headroom (about 4×) with it.
   - File an issue for vectorising the 1-D apply walk (`_loop_walk` → `visit`, about 60 % of unpreconditioned Krylov time; the DD outflow is a prefix scan).

## ⏸ COMPACTION POINT — 2026-10-04, #200 built, the operator contract next

**State:** branch `fix/krylov-sweep-preconditioner` at the docs commit after `6ccc531b`. It is stacked on `test/slow-reference-gates` (`dcc73bb1`), and neither branch is merged. `main` is `6ffd1960`, CI green.

**Read in order:**
1. this file's three sections from "The two files the memo could not shorten";
2. issue #200 and the new #575;
3. `scratch/reference_architecture/p3/elegance200/review.md`;
4. qa's evidence in `scratch/reference_architecture/p3/qa200/`.

**Then do "To do on this branch" in order.** It is a surgical carve in `orpheus/numerics/iteration.py` and `orpheus/sn/solver.py`: the main agent writes, the user steers, and the test-architect gates. Item 3 needs the user's ruling first.

### The "To do on this branch" list, landed (2026-10-05)

Each item's commit, on `fix/krylov-sweep-preconditioner`:
- **Item 1, the preconditioner is an operator:** `[LANDED c1dc0c3d]`. The stated-choice contract (no default), `(I + C) @ LC.inverse()` under DSA, the SI-Richardson witness `test_p200_4`, the merged precond-safety row, and the shared fixtures `_full_space_states.py` / `_case_slab_reference.py`.
- **Items 2, 3 and 4:** `[LANDED 90d7c337]`. ERR-053 re-gated (6 decayed markers removed; the value row `test_one_full_restart_cycle_solves_a_diffusive_infinite_medium`); the record rescaled to `‖P r‖/‖P b‖`; the 2-D Cartesian fixture.
- **The reviews of those (qa, elegance), and the user's ruling "fold scipy's acceptance in" (2026-10-05):** `[LANDED acca2416]`.
  - Both reviews found that the rescaled value is scipy's test in the FIRST restart cycle only, so a record could read converged on a solve scipy refused (60 of 400 random solves). `IterationRecord.accepted` now vetoes `converged`.
  - Also landed here: qa's F2, the exact-breakdown carve-out trusting a 0.0 tail, now confirmed against the true residual (ERR-098, which predates #200); and `test_p200_4` deriving its geometries from DSA's admission.
- **Docs:** `[LANDED 39b634a4]` (ERR-097, ERR-098, ERR-053 and ERR-050 updated; `acceleration.rst`) and `c228cb18` (the regenerated matrix and error index).
- **Item 5:** the two production docstrings `[LANDED 90d7c337]`; the precond-safety docstrings `[LANDED c1dc0c3d]`. The optional `verifies` marker went onto `test_p200_4` instead, for the new `sn-krylov-dsa-richardson` equation.
- **Issues filed:** #576 (vectorise the 1-D forward apply walk; `[M]` 46 % of the unpreconditioned GMRES time, not the "60 %" written above, and its share after #200 is not re-measured); #577 (the outer stagnates at |Δk| 3.6e-7 on qa's thick vacuum slab, identical under the identity preconditioner, so not #200's).
- **Full suite** `[M]` 2026-10-05 at `90d7c337`: 15 261 passed, 4 failed (the 4 `test_write_guards` rows, the known worktree artefact), 264 skipped, 56 xfailed, in 47 min. The run at `c228cb18` is the merge gate.
- **Note for the cadence work:** editing `orpheus/numerics/iteration.py` invalidated the cylinder A|B|A reference record (the traced memo pins the modules its generator read). The l1 file's next run regenerated it cold in 843 s. Expected behaviour, and a cost worth knowing before editing a module that references import.

**Remaining:** item 6 (the merge of both branches with `--ff-only`, and CI), then item 7 (step 3: the chord-oracle hoist).

## ⏸ COMPACTION POINT — 2026-10-05, #200 follow-up merged, the chord-oracle hoist next

**State:** `main` is `27606b17` plus this plan commit. It holds steps 1 and 2 of the user's ruling "Gates, then #200, then oracle" and the whole "To do on this branch" list. CI is green on `c228cb18`, the last code commit; `main` differs from it only by plan and memory commits. The tree is clean apart from `scratch/`.

**Read in order:**
1. this file's section "The two files the memo could not shorten", the `trajectory_resolvent_reference` bullet (the cause and the measured prototype);
2. item 7 of "To do on this branch" (in "Step 2 (#200) built; reviews in");
3. the prototype `scratch/reference_architecture/p3/perf_traj/factored.py`, with the profiles beside it (`CYLINDRICAL_solve.txt`, `CYLINDRICAL_point.txt`, `CYLINDRICAL_brute.txt`, `prof1.py`, `teeth.py`);
4. the oracle: `orpheus/derivations/continuous/trajectory_resolvent/chord_oracle.py`, `MultiRegionCylinderChordOracle.apply_operator` at line 941 (`[M]` `git grep`, 2026-10-05).

**The next step, step 3:** hoist the in-plane segment and spline evaluation out of the per-axial-cosine loop. `[M]` 2026-10-04: the prototype agrees to 9e-16 relative (solve 20 → 5.2 s, reading 53 → 4.2 s, brute 7.7 → 0.19 s).
- **Not bit-identical.** Before landing, find every consumer of frozen bytes: the traced-memo entries of the cylinder resolvent and of `TrajectoryResolventDerivation.evaluate`, any RECORD fingerprint, and any test pinning exact bytes. The memo regenerates its entries because the generator's source changes; budget minutes of cold regeneration (the cylinder A|B|A record took 843 s on 2026-10-05).
- **Apply** the `vv-principles` three-condition rule for a non-bit-exact change.
- `[R]` **Check the siblings first** (a twin-path question, Cardinal Rule 2): `CylinderChordOracle.apply_operator` (:667) and `AnnulusChordOracle.apply_operator` (:1455) may share the per-cosine rebuild. If they do, the hoist is one primitive, not three edits.
- **Review the sphere brute's headroom,** about 4× (this file, "Open from this step").
- **Mode:** the oracle lives in `derivations/`, as a Branch-1 reference. The change optimises reference code, with the user steering; the test-architect re-checks the trajectory file's gates.

**Lessons from this stretch, already recorded:**
- A `[skip ci]` tip silently skips CI for every commit of a push (memory `feedback_skip_ci_for_plan_commits.md`, rider 2026-10-05).
- Editing a module that a traced reference's generator imports invalidates that reference (this file, "Note for the cadence work").
