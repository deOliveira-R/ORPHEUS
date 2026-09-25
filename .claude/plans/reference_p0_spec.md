# #405 step 2, phase P0: the Peierls Nyström withdrawal and the defects found (verification specification)

Workflow W1, phase P0 of `.claude/plans/reference_cache.md` ("The phases, revised after the W5 review", P0's bullet; "Ruled polished (2026-09-25); P0 preparation"). Issues: #506 (the withdrawal and its return criterion), #507 (the CP flat-source heterogeneous gap), #509 (the CP sphere white-boundary defect found while writing this specification). Written by the test-architect on 2026-09-25 against `main` at `e86d249d`. That commit changes no file under `tests/`, `orpheus/` or `pyproject.toml` relative to `901f64ca`, where the explorer measured (`git diff --stat 901f64ca HEAD -- tests orpheus pyproject.toml` is empty).

Status: a specification only. No test or production code has been written. Every probe named below is a re-runnable file in `scratch/reference_architecture/p0probe/`, which is untracked. The canonical invocation is `.venv/bin/python -O -m pytest`.

Markers follow `articulation` §5: `[M]` measured (with the command), `[R]` reasoned, `[HYPOTHESIS]`, `[REFUTED 2026-09-25]`.

## 0. What changed from the plan's P0 bullet (read first)

1. **The withdrawn surface is 31 named symbols, and 26 files carry markers, not 27.** `per_observer_angular_assembly` and `per_surface_centred_angular_assembly` are the angular-assembly drivers of the escape and boundary primitives. `compute_P_esc`, `compute_G_bc` and 17 of their siblings call them (`[M]` an AST pass over `geometry.py`). #506 keeps those primitives running, so the drivers stay out of the withdrawn surface. `tests/gates/derivations/test_peierls_assembly_drivers.py` (9 cases) therefore leaves the withdrawal set `[M]` (§1.2).
2. **Of the 493 cases collected in the 27 files, 294 are withdrawn and 199 keep running** `[M]`. The count comes from dynamic attribution: every withdrawn symbol is rebound to a recorder (§1.2). The plan's "314 of 499" is the explorer's static count. It differs in predicate and in tree, and the differences are listed in §1.3.
3. **The ERR-032 catcher that exists today is a phantom** (a #455 instance). `[M]`: the documented defect was re-introduced by an arm that performed 2 rebinds. `TestSlabKernelRowSum::test_row_sum_matches_analytical_uniform_source[1.0-1.0]` stayed GREEN (it never calls the white-boundary closed form), while the unmarked `TestSlabWhiteBCInfiniteMediumIdentity` went RED 4 of 4. The replacement is mostly REUSE (§4).
4. **`[REFUTED 2026-09-25]` "a gate on the `xs_library` non-fissile χ values catches ERR-063".** Those values are zero today, which is the post-fix correct state. The `EmissionSpectrum` null law refuses any non-zero value at construction, under `-O` `[M]` (`chi_guard.py`). A value gate therefore has no realizable red input. ERR-063's catalogued root cause (sink-indexed χ) lives in `solve_peierls_mg`, which is withdrawn. Recommendation: ERR-063 goes dormant with ERR-027 to ERR-030 (§5; a user ruling, in NEEDS).
5. **The CP sphere re-anchoring is RED on today's tree, and the cause is a production defect, #509.** A homogeneous sphere with a white boundary is not flat for R ≥ 5 mean free paths, and at 2 or more groups its k is wrong. The cylinder and the slab realise the closed form to about 1e-13 `[M]` (§6). The sphere gate is specified in two forms, a strict xfail citing #509 and the plain green form; the user rules which one lands with P0.
6. **The explicit opt-in is an environment variable that names the issue it lifts:** `ORPHEUS_RUN_WITHDRAWN=506` (§2.3). It is the only spelling that reaches the subprocess worker in `cp/test_peierls_rank_n_protocol.py` and scripts run outside pytest.
7. **A marker placement is a claim, so it gets a lock that enforces it** (§2.4). This is a recommendation and needs a ruling, because it touches production (NEEDS). Without the lock the placement is enforced by nothing. A test-side substitute is blind in two measured places: an `importlib.reload` in `test_peierls_multigroup.py:572-602` undoes a rebinding, and a subprocess escapes the plugin.

## 1. The withdrawal set, per test

### 1.1 The withdrawn surface W (31 symbols)

From #506's list ("the volume kernel, the closure operators, `solve_peierls_1g`/`_mg`, the slab eigen solve, the `_build_peierls_*_case` builders and the lazy registry builders"), spelled as symbols:

- `peierls_nystrom/geometry.py` (19): `K_vol_element_adaptive`, `build_volume_kernel_adaptive`, `build_volume_kernel`, `build_white_bc_correction`, `build_white_bc_correction_rank_n`, `build_closure_operator`, `_build_closure_operator_rank2_white`, `_build_closure_operator_rank_n_white`, `_build_slab_per_face_specular_PG`, `_build_sphere_specular_mode_PG`, `_build_cylinder_specular_mode_PG`, `_build_white_rank1_mark_op`, `_build_white_f4_op`, `_build_white_hebert_op`, `_build_specular_op`, `_build_specular_multibounce_op`, `_build_full_K_per_group`, `solve_peierls_mg`, `solve_peierls_1g`; plus `BoundaryClosureOperator.__init__`.
- `peierls_nystrom/slab.py` (4): `_build_kernel_matrix`, `_build_system_matrices`, `solve_peierls_eigenvalue`, `_build_peierls_slab_case`.
- `peierls_nystrom/cylinder.py` and `sphere.py` (4): `_build_peierls_{cylinder,sphere}_case`, `_build_peierls_{cylinder,sphere}_hollow_f4_case`.
- `peierls_nystrom/cases.py` (6): `build_two_surface_case`, `_build_peierls_slab_case_via_unified`, `build_one_surface_compact_case`, `_build`, `_class_a_cases`, `continuous_cases`.

**Not withdrawn:** `continuous_case_builders` (it returns a name map and builds nothing), `capability_rows`, `SHIPPED_CLASS_A`, the naming module, `PeierlsSolution` (a data type), `reference.py` (the analytical identities), every `compute_P_*`, `compute_G_*` and `compute_T_*` primitive, the `reflection_*` matrices, `composite_gl_r`, `CurvilinearGeometry`, the two angular-assembly drivers, `ps1982_reference.py`, and `fn_method/peierls_atkinson_nystrom.py`. `[R]` The user may want `BoundaryClosureOperator` itself (a factored matrix of primitives) or the kernel's assembly identities kept running. W as written follows #506's literal list.

### 1.2 Measured per-test attribution

The instrument is `withdrawn_probe.py`, a `-p` pytest plugin. It rebinds every symbol of W, in every `sys.modules` binding, to a recorder that raises a `BaseException` subclass, and it re-installs before each test (see the reload trap below). `[M]` 2026-09-25:

```
PYTHONPATH=scratch/reference_architecture/p0probe P0PROBE_OUT=full3.json \
  .venv/bin/python -O -m pytest -p withdrawn_probe --timeout=600 \
  --deselect tests/gates/cp/test_peierls_rank_n_protocol.py::test_f4_is_sign_stable_at_its_reference_quadrature \
  --deselect tests/gates/cp/test_peierls_rank_n_protocol.py::test_f4_rich_vs_rich_panels_matches_pinned_baseline \
  <the 27 files>
```

The run took 40.6 s and printed 270 failed, 197 passed, 1 skipped, 12 deselected and 13 xfailed. The plugin made 49 rebinds. The 12 deselected cases are the two `test_f4_*` functions, which reach `solve_peierls_1g` through `_run_f4_subprocess`. A plugin cannot see into that child process; its code is a string at `:568-582`, read by hand. So those 12 are withdrawn by static reading.

The positive control: `test_peierls_sphere_white_bc.py` went 4 of 4 TOUCH (reached W) and `TestTheNameIsAFunctionOfIdentityAlone` went 10 of 10 KEPT. Every one of the 13 xfails recorded a hit.

**The reload trap is a measured instrument defect.** `test_peierls_multigroup.py` calls `importlib.reload(pc)` on `cases`, which resets every rebinding. On the first run, the 10 cases of `TestTheRatioLawIsEnforcedAtTheBOUNDARY` read KEPT only because they ran after that reload. Re-installing before each test fixed it: 2 re-installs were recorded, both inside `TestSlabViaUnifiedRoutingInfrastructure`.

**Result:** 493 cases, of which 294 are withdrawn and 199 kept. On runner 35940553034 the withdrawn cases account for 557.6 runner-minutes and the kept ones for 0.9 (0 of 493 node ids unmatched). No function is mixed: in 0 of the parametrised functions do the cases disagree. So the marker unit is the function or the class, never a `pytest.param`.

| file (under `tests/gates/`) | withdrawn | kept | placement |
|---|---:|---:|---|
| `cp/test_peierls_cylinder_flux.py` | 4 | 0 | module `pytestmark` |
| `cp/test_peierls_flux.py` | 1 | 0 | module |
| `cp/test_peierls_sphere_flux.py` | 4 | 0 | module |
| `derivations/test_peierls_convergence.py` | 5 | 0 | module |
| `derivations/test_peierls_cylinder_eigenvalue.py` | 8 | 0 | module |
| `derivations/test_peierls_cylinder_prefactor.py` | 4 | 0 | module |
| `derivations/test_peierls_fission_source_indexing.py` | 3 | 0 | module (append to the existing `pytestmark` list) |
| `derivations/test_peierls_nystrom_verification.py` | 4 | 0 | module |
| `derivations/test_peierls_rank_n_conservation.py` | 4 | 0 | module |
| `derivations/test_peierls_sphere_eigenvalue.py` | 4 | 0 | module |
| `derivations/test_peierls_sphere_white_bc.py` | 4 | 0 | module |
| `cp/test_peierls_rank_n_protocol.py` | 12 | 8 | 2 functions: the two `test_f4_*` |
| `derivations/test_continuous_registry_lazy.py` | 2 | 4 | `test_builder_keys_match_built_names`, `test_lazy_peierls_fetch_builds_requested_ref` |
| `derivations/test_peierls_closure_operator.py` | 38 | 10 | 12 functions (all but the six `test_reflection_*`) |
| `derivations/test_peierls_cylinder_multi_region.py` | 2 | 8 | class `TestMultiRegionKernel` |
| `derivations/test_peierls_cylinder_white_bc.py` | 7 | 4 | classes `TestWhiteBCEigenvalue`, `TestHebertCylinderInsufficient` |
| `derivations/test_peierls_greens_function_slab_solver.py` | 1 | 13 | `test_alpha_zero_vacuum_agrees_with_nystrom_slab` |
| `derivations/test_peierls_greens_function_xverif.py` | 2 | 6 | `test_b5_phase4_rank1_equals_white_hebert`, `test_b5_phase4_converges_toward_variant_alpha` |
| `derivations/test_peierls_multigroup.py` | 24 | 3 | 7 classes, plus `TestSlabViaUnifiedRoutingInfrastructure::test_unified_builder_produces_valid_reference` |
| `derivations/test_peierls_rank2_bc.py` | 30 | 8 | 3 classes, plus 10 functions (8 kept rows are the transmission and closed-form primitives) |
| `derivations/test_peierls_rank_n_bc.py` | 30 | 44 | 7 functions (the shifted-Legendre rows and the import-xfail row stay) |
| `derivations/test_peierls_rank_n_class_b_mr_mg.py` | 30 | 1 | 11 functions (only `test_class_b_hebert_raises_for_slab` stays, an unconditional `pytest.skip`) |
| `derivations/test_peierls_reference.py` | 30 | 66 | 8 classes |
| `derivations/test_peierls_reference_naming.py` | 10 | 13 | class `TestTheRatioLawIsEnforcedAtTheBOUNDARY` (its refusals run inside `build_two_surface_case`) |
| `derivations/test_peierls_specular_bc.py` | 26 | 1 | 1 class, plus 17 functions (only `test_specular_multibounce_slab_rank1_equals_2E3_identity` stays) |
| `derivations/test_peierls_sphere_prefactor.py` | 5 | 1 | classes `TestSphereWhiteBCRowSum`, `TestSphereRowSumIdentity` |
| `derivations/test_peierls_assembly_drivers.py` | 0 | 9 | none (leaves the set) |

Totals: 11 whole-file markers, 90 per-test marker sites in 15 files, 1 file untouched. The exact unit list, class and function names included, is `placement3.json` and `placement3.txt`. The implementer places from that file and does not re-derive it.

### 1.3 Why the static count differed (X2: two predicates)

The explorer's `pertest.py` follows each test's module-local closure to a regex of solver names; this probe records what executes. They disagree on 19 functions. The static filter could not see 15 dynamic hits reached through fixtures or helpers; the dynamic probe found 4 static hits that only read names (the `continuous_case_builders` name map). Removing the two angular-assembly drivers from W then moved 41 cases to KEPT. The explorer's 314 also covered 28 files, the PS-1982 file among them.

## 2. The withdrawal mechanism

### 2.1 The value: `Withdrawal(reason, issue)`

`[R]` A frozen dataclass whose constructor parses the marker at the boundary: `reason: str` non-empty; `issue: int` positive. It raises `ValueError` otherwise, and the message names the missing field and the marker's node id. With the lock (§2.4) it lives in `orpheus/derivations/common/withdrawal.py`, and the lock and the conftest import the same type (one definition, X4). The P4 certificate state `Withdrawn(reason, issue)` absorbs it. Without the lock it lives in `tests/_harness/withdrawal.py`.

- Marker spelling: `@pytest.mark.withdrawn("Peierls Nyström reference is not research grade (maintainer ruling 2026-09-24)", issue=506)`. The same call is used on a function, on a class, or in a module `pytestmark`. A `pytest.param(marks=…)` placement is refused at collection. It is never needed, because 0 of the functions are mixed (§1.2).
- `pyproject.toml` `markers` entry: `"withdrawn(reason, issue): the test consumes a withdrawn reference generator; it is skipped with its reason and issue unless ORPHEUS_RUN_WITHDRAWN names the issue (#506; plan reference_cache.md P0, removal trigger P4)"`.

### 2.2 The hook (in `tests/conftest.py`, beside the existing `pytest_collection_modifyitems`)

The hook runs inside the existing `pytest_collection_modifyitems`, before `registry.record`. For each item, `m = item.get_closest_marker("withdrawn")`. If `m` is present, it parses `Withdrawal.from_mark(m)`, which is a collection error on a malformed marker (§2.5). If the issue is not in `lifted_withdrawals()`, it adds `pytest.mark.skip(reason=f"withdrawn (#{w.issue}): {w.reason}")`. The `Withdrawal` is recorded on `TestMetadata` (§3.1). The skip reason is printed by `-rs`. The default summary counts `s`, so there is never a silent deselection.

### 2.3 The opt-in: `ORPHEUS_RUN_WITHDRAWN`

The value is comma-separated issue numbers, or `all`. `lifted_withdrawals() -> frozenset[int] | All` is defined once and read at call time. Why an environment variable and not a pytest option:

1. The subprocess worker in `cp/test_peierls_rank_n_protocol.py:568-615` inherits the environment; a pytest option would need a second channel.
2. The improvement work for #506 runs generators from scripts outside pytest, which reach the lock only through the environment.
3. Naming the issue keeps a lift scoped. Lifting #506 does not lift a future withdrawal filed under another issue.

Precedent in the tree: `ORPHEUS_SLAB_VIA_E1`, `ORPHEUS_B32`. W5 D requires every route-selecting variable to enter a reference's cache key. This variable is not route-selecting: it permits a run and never changes a value. So P3's key does not include it `[R]`.

### 2.4 The lock: the enforcement of the placement (recommended; NEEDS a ruling)

`[R]` Each symbol of W is decorated with `@withdrawn_generator(Withdrawal(…, issue=506))`, where the `Withdrawal` is one module-level constant in `peierls_nystrom/__init__.py`. The decorator raises `GeneratorWithdrawn(Withdrawal)` unless the issue is lifted. The docstring carries `ELEGANCE-DEBT[guard]` #506, retired by P4, which turns the lock into the `Withdrawn` state of the `ReferenceCertificate` ("the generator never runs"). The lock is what makes the ruling "they also should not run" a property of the generator (lockout-tagout), rather than of 101 marker sites. With it:

- an unmarked test that reaches W is RED, with a message that names #506 and the marker to add, so the placement is enforced for every future test;
- the call is refused at run time, so a reload cannot undo the lock and a subprocess cannot escape it. Those are the two blind spots the test-side plugin has (measured in §1.2);
- the markers remain the static declaration that the V&V accounting reads at collection time, where the lock cannot be seen.

Without a ruling for the lock, the fallback is `tests/_harness/withdrawal_tripwire.py`, a plugin in `tests/conftest.py` built as `withdrawn_probe.py` is. It must re-install at `pytest_runtest_setup`, it stays blind to the subprocess worker, and its sweep over all modules runs once per test for 12 495 tests (its cost is unmeasured).

### 2.5 The gates of the mechanism

| id | gate (file) | kind | first red in today's tree (§6c) | mutation witness |
|---|---|---|---|---|
| M1 | `Withdrawal` laws (`tests/gates/test_withdrawal.py`): accepts `("r", issue=506)`; refuses a missing issue, `issue=0`, `issue="506"` and an empty reason, each leg keyed to its argument with `match=` | THEOREM (the type's constructor law) | none shipped. The type is new, and its negative legs are its defining refusals (`vv-testing`, a math-bearing type ships its laws) | delete one clause of the parser's check → exactly that leg reds (four arms, a per-arm table) |
| M2 | `from_mark` parses real pytest `Mark` objects (`pytest.mark.withdrawn("r").mark` refused; `pytest.mark.withdrawn("r", issue=506).mark` parsed) | THEOREM | as M1 | make `from_mark` default `issue=None` to a number → the refusal leg reds |
| M3 | the skip and the opt-in, end to end: a subprocess `python -O -m pytest -rs --color=no` on `test_peierls_multigroup.py::TestMGInputValidation` (4 cases). Without the variable: 4 skipped, and the `-rs` lines contain `withdrawn (#506)`. With `ORPHEUS_RUN_WITHDRAWN=506`: 4 passed. With `ORPHEUS_RUN_WITHDRAWN=999`: 4 skipped | RECORD over the harness | the 294 markers themselves. The class is cheap either way: its cases raise on bad input before any solve | make the hook ignore the variable → the second leg reds; drop the skip → the first leg reds (the class then runs and passes, and the count of `skipped` reads 0) |
| M4 | the placement census: `TEST_REGISTRY` over the 27 files carries exactly the 294 withdrawn node ids of `placement3.json`, and the 199 kept ones carry none | RECORD, designed to redden on any placement drift | the 294 markers | remove one marker → the set differs by 1 (a named-row diff) |
| M5 (with the lock) | a kept test that reaches W is RED: the full canonical run of the 27 files (about 40 s locally with the lock, from the probe's 40.6 s) shows 199 passed or kept-skip and 294 skipped | THEOREM for the placement: no unmarked consumer | measured, not committed: the lock without the markers reds 282 in-process cases plus the 12 subprocess cases (§1.2's run is exactly this experiment) | un-mark `TestMultiRegionKernel` → its 2 cases red with `GeneratorWithdrawn` |
| M6 (with the lock) | the lock's own laws: a decorated stub refuses when not lifted, runs when lifted, and inside a `subprocess.run(sys.executable, "-c", …)` child it refuses and runs the same way | THEOREM | as M1 | make the decorator read the variable at import time (not at call) → the "set after import" leg reds |

The six gates are `foundation` (software invariants, so no `verifies`). The rows cost under 10 s together `[R]`; M3's two subprocess runs are measured at about 3 s each in §1.2's environment.

## 3. The catalogue and the V&V matrix count a withdrawn test as neither catching nor verifying

### 3.1 How they collect today `[M]`

- **The registry** (`tests/conftest.py:223-294`) records `TestMetadata(equations, catches, …)` per collected item from the resolved markers; `equations` includes labels inherited through a `case_name` parameter. Readers of the registry: `tests/_harness/audit.py` (`_equation_coverage`, `_caught_tags`, `--json` with `err_coverage`) and, through a subprocess, `tools/verification/generate_matrix.py`, which writes `docs/theory/verification/matrix.rst`.
- **The error-catalogue gate** (`tests/gates/test_error_catalogue_reconciles.py`) reads no registry. Arm 2 is an AST census: any `*.catches("ERR-NNN")` call under `tests/` counts, whatever the test's other markers are. So a skip keeps it green, which is the explorer's finding.
- **The graph** (sphinxcontrib-nexus, a separate repository) records `verifies` and `catches` from AST decorators (`ast_analyzer._parse_pytest_markers`). It records no `withdrawn`. `tools/verification/generate_error_index.py` writes the catcher counts that `vv-principles` injects from those edges.

### 3.2 The changes

1. `TestMetadata.withdrawn: Withdrawal | None = None`. Before adding it, grep `tests/_harness` and `tools/verification` for reflection walkers; `[M]` there are 0 relevant `asdict`/`vars` walkers, and `predicates.py:172` is unrelated.
2. The audit splits each relation into two: `_equation_coverage` becomes running carriers only; `_withdrawn_equation_coverage` holds withdrawn carriers; `err_coverage` becomes `{running, withdrawn}` per ERR. Derived sets: `orphan = testable − running_covered − withdrawn_only`. A label whose only carriers are withdrawn is neither covered nor orphan: it is **held by a withdrawal**, reported with its issue. `phantom_verifies` still ranges over all carriers, withdrawn ones included, because a dangling label is dangling either way.
3. `generate_matrix.py` gains a section, "Claims held by a withdrawal": label, withdrawn carrier count, issue.
4. **Reconciler arm 7** (new, in `test_error_catalogue_reconciles.py`). It reads `python -m tests._harness.audit --json` (collection about 9 s, `[M]` `pytest --collect-only` over 12 495 cases took 8.75 s wall). An entry whose catchers are all withdrawn must carry, in its body, the line `**Status:** dormant — every catcher is withdrawn under #NNN`, naming each withdrawal issue. An entry that carries the line while it has a running catcher is also red (a stale dormancy). Arm 2 is unchanged: a catcher CLAIM still exists for a dormant entry. Arm 7 adds the RUNNING question, using pytest's own marker resolution, so it does not add an AST twin of the resolution. Its positive control is inside its own body: a synthetic payload with one all-withdrawn entry and no status line must classify as a violation.
5. The catalogue text: ERR-027, 028, 029 and 030 gain the dormancy line. ERR-063 gains it too if the user rules §5's option (a). ERR-032 does not: its catcher is replaced before the withdrawal lands (§4, commit C1), so it never passes through dormancy.
6. The graph side is a separate repository, so it needs a ruling (NEEDS). The error index and `nexus errors` will keep counting withdrawn catchers until `_parse_pytest_markers` records `withdrawn` and `write_catches_edges` and `write_verifies_edges` split running from withdrawn. The ORPHEUS interim is `[R]`: `generate_error_index.py` reads the audit JSON that `generate_matrix.py` already produces (persist it once per build, never collect twice) and prints a "dormant" column, so the injected index does not claim "0 uncaught, ERR-027: 5 catchers" as coverage.

### 3.3 What they print after P0 (measured from the registry, `[M]` `audit.json` joined with the classification)

- `err_coverage`, withdrawn/running: ERR-027 5/0, ERR-028 1/0, ERR-029 6/0, ERR-030 2/0, ERR-063 3/0 → **dormant (#506)**. ERR-032 4/0 → replaced before the withdrawal, so after P0 it reads 0 withdrawn and its running catchers are E2's 4 cases plus E3's plus the foundation identity, which carries no `catches`.
- Labels whose every carrier is withdrawn (7): `peierls-mg-operator` (3), `peierls-vacuum-bc-cylinder` (3), `peierls-vacuum-bc-flux` (10), `peierls-vacuum-bc-row-sum-gate` (10), `peierls-vacuum-bc-slab` (4), `peierls-vacuum-bc-sphere` (3), `peierls-white-bc` (10). The explorer's list of 6 did not name `peierls-vacuum-bc-row-sum-gate` (it is carried by `TestSlabKernelRowSum` and siblings) and did name `flat-source` (below).
- Labels reduced but still carried (withdrawn/running): `peierls-equation` 29/4, `hebert-3-323` 3/1, `peierls-rank-n-stability` 12/8, `peierls-rank-n-bc-closure` 60/64, `peierls-unified` 94/80, `one-group-kinf` 16/130, `collision-rate` 3/91, `ki3-def` 3/61, `p-inf` 1/52, `flat-source` 1/34, `peierls-greens-slab-architecture` 1/9.
- **#387 applies to `flat-source`.** Its 34 running carriers all come from the blanket module `pytestmark = pytest.mark.verifies(…)` at `cp/test_verification.py:84` (a #387 Cartesian-product mint). By that registry predicate `flat-source` keeps 34 carriers. By the explorer's explicit, per-test predicate its only verifier was `cp/test_peierls_flux.py`. #507's "only verifier" sentence is true under the second predicate only. Both predicates travel together (NEEDS).

## 4. ERR-032, the wrong ∫E₂ antiderivative: replacement catchers

The catalogue entry is at `docs/theory/verification/error_catalog.rst:2290`. The defect was the closed form `φ_wrong = (1/(2Σt))·[2 + (2β−1)(E₂(τ)+E₂(τ'))]`, with `β = (1−E₃)/(1−2E₃)`, derived with `∫₀^τ E₂ = 1 − E₃` instead of `½ − E₃`. The shipped code today, `reference.py:335` `slab_uniform_source_white_bc_analytical`, returns `1/Σt`. So "the documented defect re-introduced" means that function returning `φ_wrong`. The arm is `err032_arm.py`: a `-p` plugin that rebinds the function in every binding (2 rebinds), refuses to install if the mutant equals the honest value at a probe point, and prints both values (honest 1.0, mutant 1.45097 at x=0.3, L=1, Σt=1).

| row | status | claim kind | what it is | first red `[M]` | honest reading `[M]` |
|---|---|---|---|---|---|
| E1 | new, `foundation` (Branch 1, `verifies` the white-BC slab label if the page carries one) | THEOREM | SymPy, fundamental theorem of calculus: `d/dτ(½ − E₃(τ)) − E₂(τ)` simplifies to 0, and `lim_{τ→0}(½ − E₃) = 0`. A NEGATIVE leg asserts that the ERR-032 candidate `1 − E₃` has the same derivative and limit ½ ≠ 0: the defect is exactly the integration constant. Do NOT use `sp.integrate(expint(2,u))`: it returns an `Ei(τ·exp_polar(iπ))` form that `simplify` does not reduce to 0 `[M]` (`sympy_e2.py`, 0.8 s) | this row cannot redden under the code arm (the identity lives in no code), so it carries no `catches` (X3); its tooth is its own negative leg | derivative residual 0; limits 0 and ½ |
| E2 | reuse `TestSlabWhiteBCInfiniteMediumIdentity` (`test_peierls_reference.py`, KEPT, 4 cases), re-tagged `catches("ERR-032")` | REFERENCE (Wigner–Seitz: a pure absorber under white reflection is the infinite medium, `Σtφ = S`, a different identity from the antiderivative) | as today | the arm reds 4 of 4 | 4 of 4 green, rel < 1e-30 |
| E3 | new, L1, `catches("ERR-032")` | REFERENCE (the Peierls equation plus the partial-current balance evaluated by `mpmath.quad`, with no antiderivative identity anywhere: `J⁺_vol = ∫₀ᴸ ½E₂(Σt(L−x'))dx'`, `J⁻ = J⁺_vol/(1 − 2E₃(τ_L))`, `φ(x) = ∫₀ᴸ ½E₁(Σt\|x−x'\|)dx' + 2J⁻[E₂(Σt x)+E₂(Σt(L−x))]`, with the E₁ integral split at x) | grid `(L, Σt) ∈ {(0.1,1),(1,1),(5,2),(100,0.5)}`, `x/L ∈ {0, .2, .5, .8, 1}`, dps 30, rtol 1e-24 | the arm's smallest signal is 5.2e-13 relative (thick slab, interior, where E₂ has decayed); every row reds | max rel 3.9e-31; 1.2 s for 20 points (`err032_balance.py`) |

Tolerance for E3: the measured honest floor is 3.9e-31 at dps 30. `1e-24` leaves 6 orders of headroom above it and 11 orders below the weakest mutant signal, and it is stated as such in the docstring. Keep the thick case `(100, 0.5)`: it is the regime where the mutant nearly vanishes, so it is the row that proves the tolerance is tight enough.

Retire the phantom: remove `catches("ERR-032")` from `TestSlabKernelRowSum`. `[M]` Under the arm, `test_row_sum_matches_analytical_uniform_source[1.0-1.0]` stayed green (23 s); it calls `slab_uniform_source_analytical` (the vacuum form), never the white one. The catalogue's "L1 test that catches it" paragraph re-points to E2 and E3 (arm 5 checks the path; `vv-principles`, a prose citation is a coverage claim). `rests_on`: E3 on nothing in the suite (mpmath's `expint` and `quad` are trusted upstream); E2 on nothing; E1 independent.

## 5. ERR-063: the catalogued defect lives in withdrawn code

`[M]` The entry (`error_catalog.rst:5293`) records two things. The first is a data edit: zeroing the non-fissile χ in `xs_library`. That edit is now the tree's state (`xs_library.py:179-278`: `chi=np.zeros(ng)` for B, C, D), and it is correct once the root cause is fixed. The `EmissionSpectrum` law also refuses anything else: `dataclasses.replace(B_2g, chi=[1,0])` raises `ValueError: EmissionSpectrum is not null` under `-O` (`chi_guard.py`). The second is the root cause: χ indexed by the SINK node in the Nyström matrix assembly (`geometry.py` `solve_peierls_mg`, `slab.py:313,368`). Its catchers are the 3 cases of `test_peierls_fission_source_indexing.py`, which are withdrawn.

- **A value gate on `xs_library` χ:** `[REFUTED 2026-09-25]` as a catcher. No input in the tree reddens it. The only edit it would refuse (a non-zero non-fissile χ) is already refused by the type at import, which fails every collection first.
- **A sibling catcher in running code:** the trajectory-resolvent solvers form the fission source locally (`chi_nodes = chi[region_at_node]; chi_nodes·F_r`, `greens_function_cylinder.py:894-906`). There is no (sink, source) pair in their emission product, so the ERR-063 mutant (`chi[i]` → `chi[j]` inside a kernel-coupled assembly) is not spellable there `[M]` (a grep of every `chi[` site outside the Nyström package; 0 kernel-coupled sites). A production-CP analogue would catch a DIFFERENT defect in a different solver; that is a new catalogue entry, not ERR-063.
- **Recommendation (a):** ERR-063 goes dormant with ERR-027 to ERR-030. Its catchers return with the solver. Option (b), a new entry plus a production-CP source-indexing gate, is out of P0's scope and belongs with CP's turn in the direction of development. This needs a user ruling (NEEDS).

## 6. The CP cylinder and sphere checks, re-anchored on the white-boundary closed form

### 6.1 What survives

Each of `cp/test_peierls_cylinder_flux.py` and `cp/test_peierls_sphere_flux.py` has 4 tests, and every one reaches W (both files are withdrawn whole):

- the two `*SelfConvergence` rows test the withdrawn generator against itself. Nothing is re-anchored for them; they return with #506.
- the two `TestCPvsPeierls*AtThickR` rows (k within 2 %, normalised shape within 5 % L2, 1 group, R = 10 MFP) are the production claims. They survive as the rows below, in a NEW module, `tests/gates/cp/test_white_boundary_infinite_medium.py`, which stays green after the withdrawal because it names nothing from W.

The closed form: a homogeneous cell whose boundary re-emits its outgoing partial current isotropically (white, Mark) is an infinite medium. The flux is flat in space in every group, its group spectrum is the infinite-medium eigenvector, and `k = k∞` (Wigner–Seitz). For discrete CP this holds EXACTLY, not to discretisation error: with `P_inf` rows summing to 1 (the white fill) and reciprocity `Σt V_i P_ij = Σt V_j P_ji`, the uniform vector is a fixed point of the CP equations `[R]`. So a CP white-boundary homogeneous solve realises the closed form to rounding. **Checked before asserting** `[M]` (`cp_flat.py`, `cp_flat2.py`, `-O`, mixture A, `keff_tol = flux_tol = 1e-12`), 38 configurations, 26 of which realise the closed form (15 of `cp_flat.py`'s 18, all 9 of `cp_flat2.py`'s, 2 of `cp_sph.py`'s 11):

| geometry | R (MFP) and mesh | groups | \|k/k∞ − 1\| | max flatness `\|φ/⟨φ⟩−1\|` | spectrum error |
|---|---|---|---|---|---|
| cylinder | 0.5, 2, 10 (uniform 5/10/20 cells); 0.5, 4 graded; 2 one-cell | 1, 2, 4 | ≤ 1.5e-16 (1G), ≤ 1.13e-13 (2G, 4G) | ≤ 1.1e-15 | ≤ 4.9e-13 |
| slab (`bc_left = bc_right = white`) | 0.5, 4 graded; 2 one-cell | 2 | ≤ 1.13e-13 | ≤ 3.0e-14 | — |
| sphere | 0.5, 2 (uniform), 0.5, 4 graded, 2 one-cell | 1, 2, 4 | ≤ 1.12e-13 | ≤ 4.4e-16 | ≤ 4.9e-13 |
| **sphere** | **10 (20 cells)** | 1, 2, 4 | **0 (1G), 4.68e-2 (2G), 6.94e-2 (4G)** | **6.8e-2** | **0.29 (2G)** |

So the cylinder and the slab realise the closed form. The sphere does too up to R = 4, and breaks from R = 5 on (`cp_sph.py`, 1G: flatness 1.85e-2 at R=5, 4.8e-2 at R=10 with 10 cells, 8.4e-2 with 40 cells; the error grows under refinement). That is **#509**. `[HYPOTHESIS]` (not investigated, per the coordinator): P_cell reciprocity or conservation fails at large optical radius, and the `np.maximum(…, 0.0)` clamps at `orpheus/cp/solver.py:417-419` feed it into `P_inf`. The existing foundation rungs cannot see it: `test_row_sums_multigroup` checks `P_inf` row sums, which are 1 by the white fill's construction, and `test_reciprocity_multigroup` checks one entry pair on a thin two-region mesh (§7).

### 6.2 The rows

Configuration grid: coordinate ∈ {slab, cylinder, sphere} × R ∈ {0.5, 2, 10} (0.5 is leakage-dominated, where the white boundary carries the whole answer; the old Peierls comparison could not reach it) × groups ∈ {2, 4} (mixture A; `vv-principles` anti-#3: at 1G, k is independent of flux shape, and `[M]` the 1G sphere at R=10 shows k = k∞ exactly while its shape is wrong). Also, as the configurations the plan never mentions: a geometrically graded mesh (12 edges, `geomspace` toward the axis) and a ONE-cell mesh at R=2 for each coordinate. The one-cell row can only test k (it is flat by construction), and its docstring says so. The reference is the dense 0-D pencil `ρ(A⁻¹F)`, with `A = diag(Σt) − Σs0ᵀ` and `F = χ⊗νΣf`, assembled in the test from the raw `Mixture` arrays, independently of every CP matrix. The dense eigen primitive is shared and trusted (X4, the trusted-library line).

| row | assertion | tolerance, derived | activates | nulls |
|---|---|---|---|---|
| C1 | `\|k/k∞ − 1\| ≤ 10 × keff_tol`, with `keff_tol` read from the `CPParams` that drove the solve (1e-12) | iterative: 10 × the solver's own tolerance; measured worst 1.13e-13, so 89× headroom | the white re-entry term; the multigroup spectrum coupling | spatial heterogeneity (covered by the flat-source CP discrete-exact cases, reused; the continuous heterogeneous rung is #507) |
| C2 | per group, `max_i \|φ_ig/⟨φ_g⟩ − 1\| ≤ 1e-12` | a band measured over the population: the 26 configurations that realise the closed form, worst 3.0e-14 (a thin graded slab); the weakest mutation signal is 9.3e-5 (arm C below) | P_cell reciprocity and conservation per row, which k cannot see (arm C) | as C1 |
| C3 | the group spectrum `⟨φ_g⟩/Σ⟨φ⟩` equals the pencil's eigenvector to `10 × flux_tol` | iterative; worst measured 4.9e-13 | the scattering transpose convention (the `[from,to]` storage) | as C1 |

Mutation witnesses `[M]` (`cp_arm.py`, 2G, R=2, 10 cells, cylinder and sphere): arm A (the white fill dropped): C1 moves by +3.2e-2 and +3.8e-2, C2 by 0.56 and 0.62. Arm B (the fill halved, albedo ½): C1 moves by 1.6e-2 and 9.0e-3, C2 by 0.26 and 0.084. Arm C (a 1e-3 row-sum defect on one row of `P_cell`): **C1 is BLIND** (k moves 1e-15; the eigenvalue is invariant under this error, Mode 12), while C2 reddens (4.5e-4 and 9.3e-5). So C2 is mandatory, and C1 alone would be a gate whose invariance group contains the error it claims to catch. The sign of arm A's k shift is `[R]` an artefact of patching around the boundary-condition factory: the k formula assumes no leakage when white. It is not a finding.

### 6.3 The sphere rows under #509: two forms, the user rules

- **(i) Land with P0 as a strict xfail.** The sphere R=10 cells of C1, C2 and C3 (2G and 4G) are `pytest.param(..., marks=pytest.mark.xfail(strict=True, reason="#509: CP sphere white boundary not flat for R >= 5 MFP"))`. Strictness is asserted by introspection in a companion row, because a marker moved into `pytest.param(marks=…)` loses it if spelled wrong. They are paired with a RECORD row that is green today and turns red at the fix: `test_sphere_R10_flatness_defect_is_present` asserts that the 1G, R=10, 10-cell flatness deviation is > 1e-3 (measured 4.8e-2), with a message that says "#509 fixed: delete this row and the xfail marks". The sphere rows at R ≤ 2 are plain and green.
- **(ii) Land after #509's fix, all green.** P0's commit C2 waits for the fix. The first red of the sphere rows is then #509 itself, which the fix commit's gate run shows going green.

In both forms the first red is #509, which exists in the tree today. Form (i) keeps the xfail's reason honest because the RECORD row pins the defect's magnitude, so a partial fix that halves the deviation leaves the RECORD row green (the deviation stays above 1e-3) and the strict xfail still failing; only a full fix turns the RECORD row red and the xfail rows into passes.

## 7. The ladder

The rows below are listed with the rungs they rest on. Each new or re-tagged test carries `@pytest.mark.rests_on(<node ids>)`.

| rung | test | status | rests on |
|---|---|---|---|
| foundation | `cp/test_verification.py::TestMultiGroupProperties::test_row_sums_multigroup[2g-CoordSystem.*]` | reused | none. Note: tautological on `P_inf` (the white fill forces row sum 1); the rung that would localise #509 is a `P_cell` conservation-and-reciprocity row over the optical radius τ_R ∈ {0.5, 2, 5, 10, 20} for all cell pairs, which is #509's investigation to write, not P0's |
| foundation | `…::test_reciprocity_multigroup[2g-CoordSystem.*]` | reused (one entry pair, thin mesh) | none |
| foundation | the dense 0-D pencil k∞ for mixture A 2G and 4G (1.875, 1.4877619048) | computed in-test; `[R]` the implementer greps `tests/gates/homogeneous/` for an existing k∞ gate on mixture A to cite in `rests_on` | none |
| edge / interior | C1, C2, C3 over the grid | new | the three rows above |
| foundation | E1 (SymPy identity) | new | none |
| interior | E2 (Wigner–Seitz), E3 (quadrature balance) | E2 reused and re-tagged, E3 new | none in the suite |
| harness | M1 to M6 | new | M3 and M5 rest on M1 and M2 |

A cycle: none. A test on no rung: none.

## 8. The order of the P0 commits, and the Sphinx obligation

Branch `fix/nystrom-withdrawal` (or as the orchestrator names it). Each commit leaves the tree green under `python -O -m pytest` on the files it touches.

1. **C1 `test(derivations)`: ERR-032 replacement catchers.** E1 and E3 are new; E2 is re-tagged; the phantom tag is removed from `TestSlabKernelRowSum`; the catalogue's citation is re-pointed. The commit body carries the arm's measured red/green table: 4 of 4 plus 20 of 20 red under the arm, 1 phantom green. It lands BEFORE the withdrawal, so ERR-032 is never dormant.
2. **C2 `test(cp)`: the white-boundary infinite-medium rows**, in form (i) or (ii) of §6.3. Before the withdrawal, so CP white-boundary coverage never lapses.
3. **C3 `feat(harness)`: the withdrawal**, in ONE commit, because each gate's first red is the other half. It holds `Withdrawal`, the lock (if ruled), the marker registration, the hook, the registry field, the audit split, the matrix section, reconciler arm 7, the 101 marker sites (11 files plus 90 sites), the dormancy lines on ERR-027 to 030 (and 063 if ruled), and M1 to M6. The body records three pre-commit measurements: the lock without markers reds 294; arm 7 without the dormancy lines reds 4 (5); the full run of the 27 files gives 199 kept, 294 skipped.
4. **C4 `docs`: the present-tense-false text** of the plan's P0 bullet (`summary.rst:56`, `escape_probability.rst:76`, `cases.py:225`, the `_cylinder_k_ref` docstring, `orpheus/derivations/README.md:107`, `.gitignore:9`, `layering.rst:61`, the stale whitelist entry). Add #507's registry-versus-explicit predicate note if the orchestrator agrees (§3.3), and the `flat-source` page statement.

**Sphinx.** `matrix.rst` (`generate_matrix`, `builder-inited`) and `.claude/skills/vv-principles/error_index.md` (`generate_error_index`, `build-finished`) are generated and tracked. Every commit that changes a marker (C1, C3) or a `verifies` label rebuilds the docs (`sphinx-build -W` per the git-workflow page) and commits the regenerated files in the same commit. Reconciler arm 3 compares only ids, so it would not catch a stale count. After C3, check that `nexus errors` and the index agree with arm 7's dormant set. They will not agree until the graph learns `withdrawn` (§3.2 item 6), and the index's interim dormant column is what keeps them honest. Then run `dead_references` (the withdrawal deletes nothing, so a zero is expected) and re-run the `test-durations` workflow; the expected drop is about 557.6 runner-minutes `[R]` (§1.2).

## 9. Refuted or rejected candidates

- **A whole-file `pytestmark` for all 27 files:** rejected. It would withdraw 199 running cases (0.9 runner-minutes of primitives and closed-form tests) for no time saved.
- **`per_*_angular_assembly` in W:** `[REFUTED 2026-09-25]`. They are the escape primitives' drivers (§0 item 1).
- **Converting the lock's exception into a skip at run time, with no markers:** rejected. The audit is collection-only, so a withdrawal known only after running cannot reach the matrix or the catalogue. The subprocess case surfaces as `CalledProcessError`, not as the lock's type.
- **An AST resolution of `withdrawn` inside reconciler arm 2:** rejected. It would twin pytest's marker resolution (module, then class, then function; X4). Arm 7 reads the resolved registry instead.
- **A SymPy identity as ERR-032's catcher:** rejected as a `catches` carrier. The identity lives in no code, so the documented defect cannot redden it (X3). It stays as the foundation rung E1.
- **The `xs_library` χ value gate for ERR-063:** refuted (§5).
- **Keeping the 1G CP rows for k:** rejected. The 1G k is blind to the #509 shape defect `[M]`.

## NEEDS (for the user, through the orchestrator)

1. **The lock (§2.4):** a production decorator on the 31 symbols of W, tagged `ELEGANCE-DEBT[guard]` #506 and retired at P4. The alternative is a test-side tripwire with two measured blind spots.
2. **ERR-063 (§5):** (a) dormant with ERR-027 to 030 (recommended), or (b) a new entry and a CP source-indexing gate later.
3. **#509's sphere rows (§6.3):** form (i), a strict xfail plus a RECORD row landing with P0, or form (ii), P0's commit C2 waiting for the fix.
4. **The graph side (§3.2 item 6):** a sphinxcontrib-nexus change so that `withdrawn` splits catcher and verifier counts, or accept the ORPHEUS-side dormant column in the index until then.
5. **W's scope (§1.1):** whether `BoundaryClosureOperator` and the volume kernel's own assembly identities (cheap, and the improvement work's ladder) are withdrawn as #506 literally lists, or kept running.
6. **#507's premise (§3.3):** `flat-source` keeps 34 carriers under the registry predicate, all minted by a #387 blanket `pytestmark`; the issue's "only verifier" is the explicit-marker predicate.

## Rulings on the NEEDS (the user, 2026-09-25: "I agree with the 5 recommendations")

1. The lock: yes, the production decorator on the 31 symbols of W, `ELEGANCE-DEBT[guard]` #506, retired at P4.
2. ERR-063: (a), dormant with ERR-027 to 030.
3. #509's sphere rows: form (i), a strict xfail citing #509 plus the RECORD row, landing with P0.
4. The graph side: the ORPHEUS-side dormant column in the error index now, and an issue in the sphinxcontrib-nexus repository for the `withdrawn` marker.
5. W's scope: withdrawn as #506 lists, `BoundaryClosureOperator` and the volume kernel's identities included; `ORPHEUS_RUN_WITHDRAWN=506` runs them for the improvement work.
6. #507's premise: corrected on the issue by comment (2026-09-25); no ruling needed.
