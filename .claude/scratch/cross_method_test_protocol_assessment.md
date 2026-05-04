# Cross-Method Test Protocol — Architectural Assessment

**Date**: 2026-05-03
**Branch**: `test/cross-method-regression`
**Author**: test-architect

## Context

User directive: build a comprehensive trajectory_resolvent ↔ fn_method
regression net **before** the planned trajectory_resolvent
architectural refactor. If protocol is not shareable, fix the
architecture first.

## What's already in place (the good news)

The codebase has the building blocks:

1. **Production-protocol case schema** in
   `orpheus/derivations/continuous/sood_registry/la13511.py` —
   `La13511Case` carries `materials: dict[int, Mixture]`,
   `mesh_template: MeshTemplate`, `truth: La13511Truth`,
   `scattering_order`, `sood_table`, `primary_reference`. This is the
   same schema production solvers (`solve_cp`, `solve_sn`, MOC)
   already consume.

2. **Truth-traceability**: every La13511 case has
   `primary_reference` (e.g., "LA-13511 Eq 20", "LA-13511 Table 4")
   and `notes` flags (typo flags, convention drift, etc.).

3. **Cross-method extractors**: `mixture_to_fn_arrays()` exists
   already (for the F_N consumer); the same `Mixture` object is
   directly consumable by trajectory_resolvent's `solve_greens_function_*`
   APIs (which take `sigma_t`, `sigma_s`, `nu_sigma_f` floats / arrays).

4. **Existing cross-method tests** (proof of concept):
   - `test_fn_la13511_slab_xverif.py` — F_N slab vs trajectory_resolvent slab
     at Sood Ua-1-0-SL + Grandjean-Siewert c=1.50. Marked `@l1
     @slow`. **The good prototype**.
   - `test_fn_la13511_sphere_xverif.py` — F_N sphere vs
     trajectory_resolvent sphere at Sood Ua-1-0-SP. Marked `@l1`.
     **The other good prototype**.
   - `test_fn_sood_table10_symmetric_pu_h2o.py` — F_N reflected slab
     vs Sood Table 10 (no trajectory_resolvent — only fn_method
     supports reflected slab).

5. **Capability matrices** (per-method): both fn_method and
   peierls_nystrom expose `capability_rows()`. Pattern is already
   established for documenting "what each method ships".

6. **wave3 plan** (`.claude/plans/wave3/architecture.md`) sketches a
   per-paper registry + meta-registry but is **deferred** until two
   concrete papers (KLL, Dahl-Sjöstrand) are implemented. The plan
   pre-empts the "smuggle non-Sood papers into sood_registry" anti-
   pattern (Atalay 1997 currently violates it).

## What's NOT shareable (the gaps)

The existing cross-method tests share *style* but not *infrastructure*:

### Gap 1 — case spec is per-test, hand-built

Every cross-method test re-builds its own (c, n_modes, n_x, n_mu)
parameters in test code. There is no abstract `CrossMethodCase`
spec that says "this IS the bare-critical c=1.30 slab; here is its
truth; here are the per-solver tolerances".

**Symptom**: adding a new c-value to the slab matrix requires
copy-paste of ~50 lines per added solver in each test. Adding a new
solver (e.g. spectral_resolvent) requires touching every cross-method
test file.

### Gap 2 — solver adapters are inline, not abstract

The pattern in `test_fn_la13511_slab_xverif.py` is:

```python
res_fn = solve_fn_slab_bare_critical(c=..., n_modes=10)
res_va = solve_greens_function_slab(
    L=2.0 * a_fn / sigma_t, sigma_t=..., sigma_s=..., nu_sigma_f=...,
    alpha=0.0, n_x=48, n_mu=128, n_traj_quad=96,
    max_iter=500, tol=1e-9,
)
```

This is **two different signatures** with hand-translated unit
conversions (`a_fn / sigma_t`) sprinkled in test code. There is no
`SolverAdapter.solve(case) -> Result` protocol.

**Symptom**: the unit-conversion logic ("F_N returns a_c in mfp;
trajectory_resolvent takes L in cm; multiply by 2 because L is full
slab not half-thickness") is duplicated everywhere. ERR-prone.

### Gap 3 — tolerance per solver is implicit

Tests hand-pick tolerances inline (`< 5e-5`, `< 1e-4`, `< 1e-3`)
based on the human author's knowledge of each solver's quadrature
floor. There is no `case.tolerance_for(solver)` lookup.

**Symptom**: when trajectory_resolvent's accuracy improves (or
regresses) under refactor, the tolerance lives in test code and is
not centrally tracked.

### Gap 4 — no case-set inventory for cross-method coverage

There is no answer to "which (case, solver) pairs are covered by
cross-method tests?" without grep-walking every `test_*xverif*.py`.

**Symptom**: the regression net has holes only visible by audit, not
by data structure.

### Gap 5 — no extension to discrete production solvers (CP, MoC, SN)

Discrete production solvers consume `(materials, mesh, params)` per
their own protocols. There is no shared way for a cross-method case
to dispatch to *both* trajectory_resolvent (continuous reference)
and CP (production discrete) for the same physical problem.

**Symptom**: discrete-solver verification cases live in
`tests/cp/`, `tests/sn/`, etc., disconnected from the continuous-
reference case set even though they are the same physical problem.

## Decision: lightweight cross-method protocol

The wave3 plan envisions per-paper sub-registries + meta-registry as
the long-term architecture. That work is rightly deferred (per
`feedback_unify_after_two_instances.md` — wait until ≥2 instances).

For **THIS** task we add a **minimal** cross-method test layer that:

- Reuses the existing `La13511Case` schema (do NOT reinvent).
- Adds a thin `CrossMethodCase` wrapper that bundles a
  `La13511Case` (or other registry case) with **per-solver
  tolerances** and **per-solver adapters**.
- Adds a `SolverAdapter` protocol with `solve(case) -> ScalarResult`
  semantics. Each adapter handles unit conversions internally.
- Lives in `tests/cross_method/` (new folder).

This is additive — does NOT touch existing cross-method tests, does
NOT reinvent `La13511Case`, does NOT prematurely build the wave3
meta-registry. It is the **minimum** infrastructure to make
trajectory_resolvent ↔ fn_method comparison *easily extensible* —
which was the user's hard requirement.

When a second use-case lands (e.g. spectral_resolvent or
singular_eigenfunction) the lightweight protocol gets *promoted*
into `orpheus/derivations/registry/` per the wave3 plan. Until then
it lives in `tests/cross_method/`.

## What goes where

### `tests/cross_method/__init__.py`
Empty (package marker).

### `tests/cross_method/protocol.py`
- `SolverAdapter` Protocol class (typing.Protocol).
- `ScalarResult` dataclass: `k_eff`, `tag` (e.g. "k_eff" or
  "critical_dimension_mfp"), `metadata` dict (iterations,
  convergence flag, etc.).
- `CrossMethodCase` dataclass: `case_id`, `registry_case` (the
  `La13511Case` or analog), `tolerances: dict[str, float]`
  (solver_name → tolerance), `claim_layer: Literal["eigenvalue",
  "flux-shape", "convergence-order"]`, `pillar: Literal["closed-form",
  "MMS", "semi-analytical"]`.
- `agree(actual, expected, tol)` helper.

### `tests/cross_method/adapters.py`
- `FNSlabAdapter`, `FNSphereAdapter`, `FNReflectedSlabAdapter` —
  thin wrappers around `solve_fn_*` functions.
- `TrajectoryResolventSlabAdapter`, `TrajectoryResolventSphereAdapter`,
  `TrajectoryResolventSphereMRAdapter` — thin wrappers around
  `solve_greens_function_*` functions.
- Each adapter: `name: str`, `solve(case) -> ScalarResult`. Each
  handles unit conversions and parameter selection (n_modes for
  fn_method, n_r/n_mu/n_traj for trajectory_resolvent) internally.

### `tests/cross_method/cases.py`
- `BARE_CRITICAL_SLAB_CASES`: list[CrossMethodCase] for the c-sweep
  on bare-critical slab (c ∈ {1.10, 1.30, 1.50, 1.70, 1.90}). Truth
  source: KLL Table I (via Sood) for c=1.30, GS Table XI for the
  rest.
- `BARE_CRITICAL_SPHERE_CASES`: c-sweep for sphere. Truth: KLL
  Table V (via Sood) for c=1.30; Sood Table 8 for PUb / UD2O.
- `REFLECTED_SLAB_CASES`: NM 1980 / Sood Table 10 cases (only
  fn_method supports natively).
- `MULTI_REGION_FIXED_SOURCE_CASES`: Garcia 2021 Case 1 (only
  trajectory_resolvent supports — fn_method does not have a
  multi-region fixed-source path). Logged as one-sided coverage.

### `tests/cross_method/test_eigenvalue.py`
The actual cross-method gates. One parametrised test per
case-set × solver. Solver-pair agreement tests.

### `tests/cross_method/test_fixed_source.py`
Garcia 2021 Case 1 — trajectory_resolvent only (one-sided coverage,
documented).

### `docs/testing/cross_method.rst`
Architecture page documenting the protocol — the convention is
durable per Cardinal Rule 3.

## Out of scope

- The wave3 per-paper meta-registry. Stays deferred.
- Refactoring sood_registry into a paper-agnostic `registry/`
  package. That's wave3's job; this task only consumes the
  existing schema.
- Discrete production-solver adapters (CP, MoC, SN). The protocol
  is *designed* to admit them but populating CP/MoC adapters is
  a follow-on. The shape of `SolverAdapter` is generic enough to
  accept a discrete-solver wrapper that receives `(materials,
  mesh, params)` and returns `ScalarResult` — but actually doing
  it requires the discrete solvers to commit to converged
  `(materials, mesh, params)` for each case (Richardson
  extrapolation parameters, etc.) and that work is not in this
  task's scope.

## Pillar / claim discipline applied to this work

Per `vv-principles` §1.5: every case row declares its
`claim_layer` and `pillar`.

For trajectory_resolvent ↔ fn_method comparisons:

- **All cases are eigenvalue (or critical-dimension) claims**
  → claim_layer = "eigenvalue".
- **Truth source is closed-form** (KLL/Sood values are exact
  transcendental dispersion-relation roots, not MMS) →
  pillar = "closed-form" for L1 backing.
- **The cross-method test itself is L4** (code-to-code agreement
  per `vv-principles` §"V&V level taxonomy"). It produces zero
  correctness info on its own. Its L1 backing is each method's
  agreement with the L1 closed-form truth — which is what each
  method already has via existing tests.
- **Structural independence is real**: F_N (Case eigenfunctions
  + Wiener-Hopf) and trajectory_resolvent (bouncing-trajectory
  Green's function) share only numpy/scipy at the trusted-
  library line. Their agreement at the truth value is the
  load-bearing structural-independence pillar.

This means: cross-method tests are tagged `@pytest.mark.l4` AND
each case's truth value is also independently verified at L1
through fn_method's existing tests against closed-form truth.

## Anti-patterns to avoid (per `vv-principles`)

- ✗ Tagging cross-method tests as L1 (they are L4 — code-to-code
  agreement; correctness info comes from the L1 truth).
- ✗ Using the trajectory_resolvent result as a "truth" for fn_method
  (or vice versa) — that's the canonical reference-contamination
  anti-pattern.
- ✗ 1-group only — every protocol case must have ≥2G in the matrix
  where the methods support it. (1G eigenvalue is degenerate per
  `vv-principles` §"1-group degeneracy"; cross-method comparison
  on 1G alone has no diagnostic value beyond cross-implementation
  agreement.)

The trajectory_resolvent ↔ fn_method bare-critical sphere/slab is
**inherently 1G** — neither method natively ships the multi-group
critical-dimension solve that would close this gap. Multi-group
cases land via:

- **fn_method**: 1G k_inf is closed-form (rational algebra); 2G
  k_inf is closed-form via 2x2 eigenvalue. fn_method does NOT
  ship multi-group critical-slab/sphere.
- **trajectory_resolvent**: ships `solve_greens_function_sphere_mg`,
  `solve_greens_function_slab_mg`. NO direct fn_method counterpart
  for cross-comparison at multi-group.

So the multi-group cross-comparison story is:

- **k_inf**: fn_method `compute_kinf_mg` ↔ trajectory_resolvent
  closed sphere (specular BC) at α=1, where `k_eff = k_inf`
  exactly. Already pinned in
  `test_peierls_greens_function_xverif.py::test_b5_variant_alpha_gives_k_inf_exactly`
  for 1G; we extend to 2G.

This documents the ≥2G coverage gap honestly and uses what's
shippable today.

## Test count target

Phase 3 case populations:

- Bare slab: 5 c values × 2 solvers + 5 agreement gates = 15 tests.
- Bare sphere: 4 c values × 2 solvers + 4 agreement gates = 12 tests.
- Reflected slab: 4 cases × 1 solver (fn_method only) = 4 tests
  (truth-only; documents one-sided coverage).
- Multi-region fixed-source: 1 case × 1 solver
  (trajectory_resolvent only) = 1 test (one-sided; already covered
  by garcia2021 test, only mentioned for completeness).
- 2G k_inf cross-check (closed sphere α=1): 2 cases × 2 solvers + 2
  agreement = 6 tests.

Total: ~35 tests. Many are parametrised so total test-row count
~25-30 actual `def test_*` functions.

## Hard constraints from CLAUDE.md and skills

- Cardinal Rule 1 (correctness): every truth value traces to a
  primary source. ✓ — KLL/Sood/GS/NM citations carried in each
  CrossMethodCase.
- Cardinal Rule 2 (architecture): if shared concepts emerge,
  abstract them. ✓ — that's exactly what this task does for
  cross-method case+adapter protocol.
- Cardinal Rule 3 (Sphinx): document the protocol durably. ✓ —
  `docs/testing/cross_method.rst`.
- `vv-principles`: pillar tags + claim layers carried per-case. ✓.
- `algebra-of-record`: each truth has structural-independence
  status documented (which methods are independent of which). ✓.

## Next phase: implementation

1. Phase 2 (architecture): build `tests/cross_method/protocol.py` +
   `adapters.py`. Commit on the branch.
2. Phase 3+4 (cases + tests): build `tests/cross_method/cases.py` +
   `test_eigenvalue.py` etc. Run them. Commit.
3. Phase 5 (Sphinx): `docs/testing/cross_method.rst`. Commit.
