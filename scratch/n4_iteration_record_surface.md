# N4 surface map — CP / MoC / diffusion vs `IterationRecord`

Measured against HEAD `52650a86` (main), 2026-08-10, host `.venv/bin/python`.
`git diff --stat` over `orpheus/{cp,moc,diffusion,numerics}` at open AND close:
**empty** — no in-flight carve underneath this audit (L-012). All four mapped
production files are tracked at HEAD.

⚠ `.claude/worktrees/nexus-workspace-wiring/` carries a second copy of every
file below. Every `file:line` here is the MAIN checkout.

Nexus was not consulted for the dataclass-field questions: per L-009 the graph
does not model dataclass fields as nodes, so text-grep is the primary evidence
for a field-level audit. The structural claims below are all measured by
running, not read.

---

## Q6 (LEAD) — can the SHIPPED primitive express the three ragged shapes?

**Two of three: yes, unchanged. MoC is expressible as DATA but the primitive
DERIVES `converged = False` for it — deliberately, and the rule is pinned by a
test. That is N4's one design decision.**

Probe (constructed each shape directly against `orpheus/numerics/convergence.py`):

| shape | constructs? | `status` | `converged` | parent's `fully_converged` |
|---|---|---|---|---|
| diffusion — outer with **no children** | ✅ | CONVERGED | True | **True** |
| diffusion — explicit `IterationRecord(label="inner(exact resolvent, LU)")` | ✅ | **DIRECT** | True | **True** |
| CP — outer with **3 per-group children**, each with a residual criterion | ✅ | CONVERGED | True | **True** |
| CP — per-group child with **count only** (`criteria=()`, `iterations_run=n`) | ✅ | STOPPED | **False** | **False** |
| MoC — `criteria=()`, `iterations_run=15`, `budget=0` | ✅ | STOPPED | **False** | **False** |
| MoC — `criteria=()`, `iterations_run=15`, `budget=15` | ✅ | **TRUNCATED** | **False** | **False** |

### The invariant that blocks MoC

`orpheus/numerics/convergence.py:796-802` — `IterationRecord.converged`:

```python
if self.n_iterations < self.min_iterations:
    return False
if self.iterated and not any(
    criterion.n_iterations > 0 for criterion in self.criteria
):
    return False
return all(criterion.cleared for criterion in self.criteria)
```

The **second clause** is the blocker, and it is deliberate. Its own docstring
(`convergence.py:784-794`) records that the rule first shipped as
`any(criterion.n_iterations == 0 …)` — silent on the no-criteria case, since
`any(())` is `False` — and was **widened on purpose** to close exactly this
hole:

> ⛔ That rule first shipped as `any(criterion.n_iterations == 0 …)`, which is
> silent on the case of **no criteria at all** — `any(())` is `False`, so a
> level that iterated while declaring nothing to measure fell through to
> `all(())` and claimed convergence. That is the same vacuous-conjunction lie
> the power-iteration loop … refuses at the producer, left open one layer down.

The primitive's position: **a level that RAN and MEASURED NOTHING may not claim
convergence.** MoC's inner is exactly that level.

The rule has **three independent enforcement points**, so it is not a
one-line change:

1. `IterationRecord.converged` — `convergence.py:798-801` (above).
2. `power_iteration` refuses a solver that returns no readings —
   `orpheus/numerics/eigenvalue.py:437-446`, `raise ValueError(… "a loop with
   nothing to measure cannot converge, and an empty conjunction would claim it
   did")`.
3. A **pinning test**: `tests/numerics/test_iteration_record.py:499-527`,
   `test_a_level_that_RAN_and_measured_NOTHING_cannot_claim`. Also
   `tests/numerics/test_power_iteration_record.py:300`
   (`test_a_solver_that_measures_nothing_is_refused`).

### What the MoC tree actually prints if you record it as-is

```
outer(power-iteration): CONVERGED (3/6 iterations)
  dk: 2.200e-16 vs tol 1.000e-06  met <- binding
  dphi: 1.300e-16 vs tol 1.000e-05  met
  inner(moc sweeps, outer 1): TRUNCATED (3/3 iterations)
  inner(moc sweeps, outer 2): TRUNCATED (3/3 iterations)
  inner(moc sweeps, outer 3): TRUNCATED (3/3 iterations)
converged: True   fully_converged: False   first_failure: inner(moc sweeps, outer 1)
```

An **empty-bodied TRUNCATED block** — a truncation report with no criterion
line under it, because there is nothing to print. That is the design problem
made visible.

### The three ways out (none needs a structural change to `IterationRecord`)

1. **Accept `converged=False`.** Consequence: every MoC solve reports
   `fully_converged=False` forever, `first_failure` always names a MoC inner,
   and — if N4 also wires the SN-style warning — **every MoC solve warns**.
   Note the SN warning helper *already has a branch for this state*:
   `orpheus/sn/solver.py:499-504`, "every recorded criterion cleared but the
   level measured nothing — it ran N iterations without recording a value, so
   there is nothing to project from". So the message exists; it would just
   fire on every MoC run.
2. **Give the inner a real criterion.** The data is computed and thrown away:
   `orpheus/moc/core.py:190-198` builds `phi_new` and does `phi = phi_new` with
   no comparison. One `‖phi_new − phi‖/‖phi‖` per sweep makes MoC's level an
   ordinary criterion-bearing level — measured: `fully_converged` returns True
   and the child reads `CONVERGED`. This converts the design question into
   *plumbing + one measurement*, and it is the same quantity MoC's OUTER
   already measures (`moc/core.py:243-264`).
   ⚠ MoC's inner has **no tolerance today** — option 2 has to mint one
   (a new `MOCSolver.__init__` parameter), which is a public-API change to
   `solve_moc`.
3. **Express the fixed schedule AS the criterion** — a `sweeps_remaining`
   countdown driven below `1.0`. `StoppingCriterion` accepts it (magnitudes,
   non-negative tolerance) and it prints honestly. Keeps MoC's semantics
   identical while satisfying the invariant, at the cost of a criterion that
   measures the loop rather than the physics.

⛔ **Only option 1 is free.** 2 changes what MoC measures AND its public
signature; 3 changes what it records.

### The other two: what "expressible" costs

- **Diffusion is free, and the vocabulary already exists.**
  `IterationRecord(label="inner(exact resolvent, LU)")` with every default
  yields `iterated=False`, `converged=True`, `status="DIRECT"`. The docstring
  at `convergence.py:751-759` names the case; the test at
  `tests/numerics/test_iteration_record.py:465-483`
  (`test_a_direct_solve_did_not_iterate_but_DID_converge`) **already cites
  diffusion by name**: *"Diffusion's inner is one LU back-substitution, so its
  record is a 1-level tree whose only honest reading is this one."* Diffusion
  may equally carry no children, which is what it produces today. Both are
  valid; the explicit DIRECT child is strictly more informative and free.
- **CP fans out for free STRUCTURALLY** — `children: tuple[IterationRecord, ...]`
  has no arity constraint and `eigenvalue.py:478-482` splices whatever the
  solver accumulated — **but not with the data CP keeps today.** CP computes
  its inner residual and discards it (`orpheus/cp/solver.py:595-600`); only the
  count `n_in` survives. A count-only child is `criteria=()` + `iterations_run`,
  the SAME blocked shape as MoC. So CP's job includes capturing the residual it
  already has in hand — 2 lines at a site that already computed the number.

### One driver-side gap (all three)

`power_iteration` splices children **only** for a `RecordingSolver`
(`eigenvalue.py:387-389` and `478-482`). None of `CPSolver` / `MOCSolver` /
`DiffusionSolver` declares `inner_records`, so all three take the `else ()`
arm — **measured: `len(record.children) == 0` on real solves of all three.**
Adding one property per solver is the whole plumbing job; `power_iteration`
itself needs no change. `RecordingSolver` is `@runtime_checkable`
(`eigenvalue.py:221-250`) and explicitly names CP/MoC/diffusion as the families
that "conform to the base contract without it and plug in with no suppression".

---

## 1. The result types

### CPResult — `orpheus/cp/solver.py:73-87`

```python
@dataclass
class CPResult:
    """Results of a collision probability calculation."""
    keff: float
    keff_history: list[float]
    flux: np.ndarray
    flux_fuel: np.ndarray
    flux_clad: np.ndarray
    flux_cool: np.ndarray
    geometry: Mesh1D
    eg: np.ndarray | None  # (ng+1,) energy boundaries, or None for synthetic XS
    elapsed_seconds: float
    residual_history: list[float] = field(default_factory=list)
    n_inner: np.ndarray | None = None  # (n_outer, ng) inner iteration counts
```

Convergence status carried today: **NO `converged`, NO `record`** (measured:
`hasattr(r,"converged") is False`, `hasattr(r,"record") is False`). It DOES
carry two diagnostics — the only family that does:

- `residual_history` — the neutron **balance** residual, one per outer, appended
  inside `measure_stopping_criteria` (`cp/solver.py:756-758`). ⚠ **This is NOT a
  stopping criterion** and the code says so explicitly at `cp/solver.py:752-755`:
  *"CP tracks the residual as a DIAGNOSTIC; it is not part of the stop test, so
  it is not returned as a criterion — a returned criterion is one the loop stops
  on."* Do not map it onto a `StoppingCriterion` in N4.
- `n_inner` — `(n_outer, ng)` int array, GS mode only, `None` in Jacobi.
  Measured on `cp_slab_2eg_2rg` GS: shape `(3, 2)`, last row `[12, 1]`.

### MoCResult — `orpheus/moc/solver.py:35-62`

```python
@dataclass
class MoCResult:
    """Results of a Method of Characteristics calculation."""
    keff: float
    keff_history: list[float]
    flux_per_material: dict[int, np.ndarray]
    scalar_flux: np.ndarray
    moc_mesh: MOCMesh
    eg: np.ndarray | None
    elapsed_seconds: float
    # + derived properties: ng, flux_fuel, flux_clad, flux_cool
```

Convergence status: **NONE.** No `converged`, no `record`, no residual history,
no inner counts. Measured.

### DiffusionResult — `orpheus/diffusion/solver.py:341-358`

```python
@dataclass(frozen=True, eq=False)
class DiffusionResult:
    keff: float
    flux: "FullField"
    current: np.ndarray            # (ng, nx+1) axis-signed net currents
    keff_history: "list[float]"
    mesh: "DiffusionMesh"
```

Convergence status: **NONE.** Measured. (Note it is `frozen=True`, unlike CP's
and MoC's plain `@dataclass` — adding a field is a keyword-only addition either
way, but a frozen type is the SN `IterationHistory` shape already.)

---

## 2. Iteration structure, as it ACTUALLY is

### CP — 2 levels, per-group, WITH a tolerance (richer than the plan recorded)

- **Outer**: `power_iteration`, `eigenvalue.py:391-461`. Stop test
  `n >= MINIMUM_OUTER_ITERATIONS and all(r.cleared for r in readings)`
  (`eigenvalue.py:460`). Readings from `CPSolver.measure_stopping_criteria`,
  `cp/solver.py:723-774`: `dk` vs `keff_tol`, `dphi` vs `flux_tol`.
  ⚠ CP judges `dphi` in **ℓ∞** (`cp/solver.py:747-748`), where SN / MoC /
  diffusion use relative ℓ².
- **Inner** — `_solve_fixed_source_gs`, `cp/solver.py:539-608`. Two nested
  loops: `for g in range(ng)` (`:556`) then `for n_in in range(1, self.max_inner + 1)`
  (`:562`). **Predicate**: `cp/solver.py:596-600`

  ```python
  res_in = (np.linalg.norm(phi_g_new - phi_g_old)
            / max(phi_g_norm, 1e-30))
  if res_in < self.inner_tol:
      break
  ```

  ⟹ CP's inner has a genuine per-group tolerance (`inner_tol`, default `1e-8`)
  and budget (`max_inner`, default `100`) — `CPParams`, `cp/solver.py:69-70`.
- **Jacobi mode** (`_solve_fixed_source_jacobi`, `cp/solver.py:510-537`,
  the DEFAULT `solver_mode="jacobi"`) has **no inner loop at all** — one
  matvec. So CP's own level tree is *mode-dependent*: 2 levels in GS,
  **1 level in Jacobi**. `n_inner_history` stays empty in Jacobi and
  `CPResult.n_inner` is `None`. **N4 must model both.**

**The measurement that is discarded**: `res_in` at `cp/solver.py:596`. It is
computed every inner iteration and never stored. Only `n_in` survives
(`cp/solver.py:602` → `604` → `solve_cp` `:912-914` → `CPResult.n_inner`).

⚠ **The surviving count is not injective** (L-024 corollary). The loop is
`for n_in in range(1, max_inner+1): … if res_in < tol: break`, so
`n_in == max_inner` means *either* "converged exactly at the cap" *or*
"exhausted". A consumer reconstructing convergence from `n_inner` is guessing
on the boundary case — the exact reason the record exists.

### MoC — 2 levels, inner has NO criterion at all (plan's row CONFIRMED)

- **Outer**: `power_iteration`, called at `orpheus/moc/solver.py:137`.
  Readings from `MOCSolver.measure_stopping_criteria`, `orpheus/moc/core.py:243-264`:
  `dk` vs `keff_tol`, `dphi` vs `flux_tol` (relative ℓ²).
- **Inner**: `MOCSolver.solve_fixed_source`, `orpheus/moc/core.py:92-210`.
  The loop is `for _inner in range(self.n_inner_sweeps):` (**`core.py:111`**) —
  the loop variable is discarded (`_inner`). **No tolerance, no residual, no
  break.** `n_inner_sweeps` defaults to `15` (`moc/core.py:38`, `moc/solver.py:78`).
- Everything below that is the ray sweep (`core.py:127-184`) — spatial, not an
  iteration level.
- **Where the measurement WOULD be**: `core.py:190-198` computes `phi_new`, then
  `phi = phi_new` at `:198`. `‖phi_new − phi‖` is one line from existing.

### Diffusion — 1 level (plan's row CONFIRMED)

- **Outer**: `power_iteration`, `orpheus/diffusion/solver.py:413`.
  Readings from `DiffusionSolver.measure_stopping_criteria`,
  `diffusion/solver.py:317-338`: `dk` vs `keff_tol`, `dphi` vs `flux_tol`.
- **Inner**: `DiffusionSolver.solve_fixed_source`, `diffusion/solver.py:288-295`:

  ```python
  r"""EXACT inner solve :math:`\psi = A^{-1} q` — one LU
  back-substitution; the initial guess is irrelevant."""
  return np.asarray(
      self.resolvent.apply(np.asarray(fission_source, dtype=float))
  )
  ```

  **No loop. No iteration. No level.** Nothing is discarded because nothing is
  measured.

---

## 3. Where the outer driver comes from

**All three route through the SHARED `power_iteration`** — verified by grep and
by a spy patch on every re-export:

| family | call site | max_iter source |
|---|---|---|
| CP | `orpheus/cp/solver.py:903` | `params.max_outer` |
| MoC | `orpheus/moc/solver.py:137` | `max_outer` (default 500) |
| diffusion | `orpheus/diffusion/solver.py:413` | `max_outer` (default 1000) |

None hand-rolls an outer loop. `power_iteration` returns a
`PowerIterationOutcome` carrying a fully-built `IterationRecord`
(`eigenvalue.py:463-484`) **already, today, for all three.**

⟹ **The outer half of N4 is pure plumbing.** The record exists at the boundary
and is dropped. The design work is entirely on the *children*.

Measured records at the boundary today:

```
CP   (cp_slab_2eg_2rg, GS):   outer(power-iteration): CONVERGED (3/500 iterations)
                                dk 8.186e-08 / tol 1e-07 met <- binding
                                dphi 7.070e-07 / tol 1e-06 met
                              children: 0   fully_converged: True
DIFF (same case):             outer(power-iteration): CONVERGED (3/1000 iterations)
                                dk 0.000e+00 / tol 1e-10 met
                                dphi 9.721e-17 / tol 1e-09 met <- binding
                              children: 0   fully_converged: True
MoC  (moc_cyl1D_1eg_1rg):     outer(power-iteration): CONVERGED (3/6 iterations)
                                dk 2.220e-16 / tol 1e-06 met <- binding
                                dphi 1.349e-16 / tol 1e-05 met
                              children: 0   fully_converged: True
```

⚠ Note all three currently read `fully_converged: True` **because they have no
children**. Wiring children is what can turn that False — which is the point,
and also the behaviour change to gate.

---

## 4. The drop sites (exact lines)

| family | line | what is dropped |
|---|---|---|
| CP | `orpheus/cp/solver.py:903-906` | `outcome` is destructured to `(keff, keff_history, phi)`; `outcome.record` and `outcome.converged` never read. The `CPResult(...)` at `:919-925` has no slot. |
| MoC | `orpheus/moc/solver.py:137-140` | same shape: `outcome.keff, outcome.keff_history, outcome.flux_distribution`. `MoCResult(...)` at `:161-169` has no slot. |
| diffusion | `orpheus/diffusion/solver.py:413-416` | same shape. `DiffusionResult(...)` at `:419-425` has no slot. |

Secondary drop site inside CP: `orpheus/cp/solver.py:596-600` — `res_in`
computed, tested, discarded.

For comparison, the SN pattern N4 would replicate:
`orpheus/sn/solver.py:2463` and `:2740` — `record=outcome.record,
keff_history=tuple(keff_history)` into `IterationHistory`
(`orpheus/sn/solution.py:112-`), which is a **view** over the record: every
scalar is a derived `@property`, only `record` and `keff_history` are fields.
Its docstring at `solution.py:153-159` carries an explicit instruction worth
honouring in N4: *"Do not grow this surface: add the question to
`IterationRecord`, where the tree can answer it."*

---

## 5. Consumers

Field-level greps (L-009: the graph has no dataclass-field nodes). Counts
exclude `.claude/worktrees/`.

### CP — 37 `solve_cp(` call sites, all in `tests/`

Production consumers of the RESULT: **zero outside `orpheus/cp/__init__.py`**
(re-export). `orpheus/plotting.py` consumes `.keff_history` / `.flux_fuel` /
`.eg` / `.geometry` generically.

Consumers of the two convergence-ish fields:
- `tests/cp/test_diagnostics.py` — **the contract test for `n_inner` /
  `residual_history`.** `:62,66` residual populated + small; `:70`
  `res.n_inner is None` in Jacobi; `:104-113` `n_inner` populated with shape
  `(n_outer, 2)` and all `>= 1` in GS; `:166-173` thermal group needs `>=`
  inner iterations than fast.
- `docs/theory/methods/collision_probability.rst:1916-1930` — a `list-table`
  documenting both fields with their shapes. **A doc-node consumer that no
  Sphinx severity would catch if the fields moved** (the CP module has no
  `automodule`; per `coding-standards`, Python-domain roles in an unrendered
  module are invisible at every severity). Also `:1894` narrates
  "The `n_inner` array showed all 1s" as ERR-016 history — past tense, stays.
- ✅ `docs/theory/methods/collision_probability.rst:1906-1910` and `:1933` cite
  `tests/cp/test_verification.py::TestGSInnerIterations` with three test names.
  **All three verified alive**: class at `test_verification.py:993`,
  `test_thermal_needs_more_inner_than_fast:1007`,
  `test_gs_eigenvalue_matches_jacobi:1042`,
  `test_no_self_scatter_one_inner:1073`. This is a **second CP inner-iteration
  gate** beyond `test_diagnostics.py` — see the gating-tests table.

### MoC — 5 `solve_moc(` call sites

- `orpheus/plotting.py:219,237,280` (typed `MoCResult`) — reads
  `.keff_history`, `.flux_*`, `.eg`, `.moc_mesh`. No convergence read.
- `tests/moc/test_moc.py:49`, `tests/moc/test_properties.py:23,81`,
  `tests/moc/test_verification.py:161`. None reads a convergence field
  (there is none to read).
- `examples/method_of_characteristics/demo_moc.py:42`.
- `docs/theory/methods/method_of_characteristics.rst:118,131`;
  `docs/api/method_of_characteristics.rst:54`.
- ⚠ `student_resources/03_method_of_characteristics.py` defines its OWN
  `StudentMoCResult` and its own `solve_moc` — a namespace collision with the
  production symbol (L-017 shape). It is NOT a consumer; exclude it from any
  blast radius count.

### Diffusion — 7 `solve_diffusion_1d(` call sites

- `tests/diffusion/test_continuous_reference.py:145,177,250,288`,
  `tests/diffusion/test_properties.py:47`, `tests/diffusion/test_solver.py:109`.
- `examples/diffusion/demo_diffusion_1d.py:74`.
- `orpheus/diffusion/operators.py:392` mentions `DiffusionResult.current` in a
  docstring; `orpheus/diffusion/__init__.py:10,33,54` re-export.
- ⚠ `student_resources/06_diffusion_1d.py:98` also defines its own
  `DiffusionResult` — same collision, same exclusion.

### Cross-family: the test that ALREADY names all three

`tests/numerics/test_power_iteration_record.py:159-192`,
`test_no_realizer_can_see_the_iteration_index`, parametrizes over
`("orpheus.cp.solver","CPSolver")`, `("orpheus.moc.core","MOCSolver")`,
`("orpheus.diffusion.solver","DiffusionSolver")`. This is an existing
three-family gate that N4 touches.

---

## 7. Existing tests that would gate an N4 change

### The primitive itself (change these and you have changed the contract)

- `tests/numerics/test_iteration_record.py` — the whole `IterationRecord` /
  `StoppingCriterion` contract. Blocking rows for N4:
  - `:445-463` `test_each_state_is_reached_and_named` — parametrizes
    `((), 0, "DIRECT")`, the diffusion shape.
  - `:465-483` `test_a_direct_solve_did_not_iterate_but_DID_converge` —
    **names diffusion explicitly** in its docstring.
  - `:499-527` `test_a_level_that_RAN_and_measured_NOTHING_cannot_claim` —
    **the invariant that blocks MoC.**
  - `:572` `test_a_level_with_no_criteria_has_no_binding_one`.
- `tests/numerics/test_power_iteration_record.py` — the driver contract.
  - `:159-192` the three-family min-outer parametrization (above).
  - `:298-318` `TestTheLoopRefusesWhatItCannotJudge` /
    `test_a_solver_that_measures_nothing_is_refused`.
  - `:359-405` `TestTheSubtreeIsWhatTHISLoopCaused` —
    `test_a_non_recording_solver_yields_a_leaf` (`:361-364`) asserts
    `power_iteration(solver).record.children == ()` for a NON-recording solver.
    Making CP/MoC/diffusion into `RecordingSolver`s does not red this (it uses
    its own stub), but it is the law N4's new children must satisfy:
    `len(children) == record.n_iterations` (`:395`).
- `tests/numerics/test_default_iteration_budget.py`, `tests/numerics/test_iteration.py`,
  `tests/numerics/test_eigenvalue.py`.

### The SN template N4 replicates

- `tests/sn/solve/test_convergence_contract.py` — the whole shape of what a
  wired family owes: `TestPowerIterationCarriesItsOutcome` (`:282`),
  `TestFixedSourceFlagIsHonest` (`:340`), `TestTruncationIsAudible` (`:361`)
  incl. `test_starved_solve_warns` / `test_converged_solve_is_SILENT` /
  `test_it_names_the_LOOP_when_every_criterion_cleared` (`:529`) /
  `test_the_published_escalation_flag_actually_parses` (`:574`).
  `TestEveryEntryDerivesItsBudget` (`:128`) is the `resolve_iteration_budget`
  contract — ⚠ **CP / MoC / diffusion still ship hardcoded budgets**
  (`CPParams.max_inner=100`, `MOCSolver.n_inner_sweeps=15`, and the
  `max_outer` defaults 500 / 500 / 1000); the SN gate at `:175-180`
  enumerates SN entries only.
- `tests/sn/primitives/test_solution.py` — `IterationHistory` as a view;
  `:184` `h.n_inner == 0`, `:295` `h.n_inner is None` on the eigenvalue path.

### Per-family gates

| family | file | what it pins |
|---|---|---|
| CP | `tests/cp/test_diagnostics.py` (whole file, 175 lines) | `residual_history` populated + final < 1e-3; `n_inner is None` in Jacobi; `n_inner` shape + `>= 1` in GS; thermal ≥ fast inner counts. **The single most N4-relevant file.** |
| CP | `tests/cp/test_verification.py:993-1080` `TestGSInnerIterations` | **the second inner-iteration gate.** `test_thermal_needs_more_inner_than_fast:1007`, `test_gs_eigenvalue_matches_jacobi:1042`, `test_no_self_scatter_one_inner:1073` (⚠ the last asserts the inner converges in exactly **1** iteration with no self-scatter — a direct `n_inner` consumer). Cited from `collision_probability.rst:1906-1910`. |
| CP | `tests/cp/test_diagnostics.py:130-144` | 27-case GS ≡ Jacobi eigenvalue sweep — the regression net for any inner-loop change |
| CP | `tests/cp/test_properties.py`, `test_slab.py`, `test_cylinder.py`, `test_sphere.py`, `test_peierls_*.py` | eigenvalue / flux physics |
| MoC | `tests/moc/test_moc.py:43-59` | `test_moc_homogeneous` × 3 cases, `k_inf` to 1e-4 — the smoke gate |
| MoC | `tests/moc/test_properties.py:30,47,53,66` | particle balance, positivity, flux/material consistency, depression |
| MoC | `tests/moc/test_verification.py`, `tests/moc/test_mms.py` | L1 |
| diffusion | `tests/diffusion/test_solver.py` | `TestEngineCrossGate:135` power-vs-direct, `TestInfiniteMedium:173`, trace semantics `:219-310`, `TestBalanceAndOrdering:312`, `TestProtocolAndDefaults:370` incl. `test_solver_satisfies_the_production_rate_contract:371` — **the closest existing "does this solver conform to a numerics Protocol" gate; the natural home for a `RecordingSolver` conformance assertion** |
| diffusion | `tests/diffusion/test_properties.py`, `test_continuous_reference.py`, `test_mms.py` | physics + L1 |

### ⚠ Marker migration — `n_inner` carries an error-catalog edge

`tests/cp/test_verification.py:1006` tags
`test_thermal_needs_more_inner_than_fast` with
`@pytest.mark.catches("ERR-016")` — the tautological-inner-residual bug
narrated at `collision_probability.rst:1883-1914`. If N4 retires or re-points
`CPResult.n_inner`, the `catches("ERR-016")` marker must move with the
assertion (`coding-standards` §"Retirement means MARKER migration too"), or the
V&V matrix silently loses the ERR-016 coverage edge while
`error_catalog.md` keeps naming a test whose measurand is gone.

### Doc pages N4 must update (they state the current, soon-false, structure)

- `docs/theory/methods/collision_probability.rst:1916-1930` — the `CPResult`
  diagnostics `list-table`.
- `docs/theory/methods/collision_probability.rst:1841-1881` — the Jacobi/GS
  mode section, incl. the inner-tolerance predicate written out at `:1880-1881`.
- `docs/theory/methods/method_of_characteristics.rst:865-867` — *"The method
  performs `n_inner_sweeps` transport sweeps per outer iteration"*, and
  `:877-887` (the `measure_stopping_criteria` section, already carrying its
  ⛔ N2b note — the model for N4's own note).
- `docs/theory/methods/diffusion_1d.rst` — not audited in this dispatch.

---

## Refuted / not-verified, recorded so nobody re-derives them

- ⛔ **The plan's F6 line numbers are all stale**: `cp/solver.py:884` is now
  `:903`; `moc/solver.py:139` is now `:137-140`; `diffusion/solver.py:409` is
  now `:413`. The CLAIM ("all three receive `outcome.converged` and drop it")
  is **CONFIRMED** at the new lines.
- ⛔ **The plan's §1b CP row understates CP.** It says CP is "the only family
  already reporting a per-outer × per-group inner *count*, still status-free" —
  true — but CP's inner also **has a real tolerance and budget**
  (`inner_tol=1e-8`, `max_inner=100`) and computes a real residual it throws
  away. CP is therefore the family CLOSEST to expressible, not a special case;
  it needs 2 lines of capture, not a design decision.
- ⛔ **CP's level count is mode-dependent** (2 in Gauss-Seidel, **1 in Jacobi**,
  which is the DEFAULT). The plan's "CP | 2, per-group" row is right for GS only.
- ✅ MoC's row ("2, no criterion") and diffusion's ("1") verified exactly.
- **Not verified**: `docs/theory/methods/diffusion_1d.rst`'s claims about the
  iteration structure; and whether
  `tests/cp/test_verification.py::TestGSInnerIterations` (cited at
  `collision_probability.rst:1906`) still exists under that class name.
- **Homogeneous is NOT a fourth family** — grep of `orpheus/homogeneous/` for
  `power_iteration` / `EigenvalueSolver` returns zero. The five
  `measure_stopping_criteria` realizers in the tree are exactly:
  `numerics/iteration.py:1443` (KEigenvalue), `sn/solver.py:1645`,
  `cp/solver.py:723`, `moc/core.py:243`, `diffusion/solver.py:317`.
