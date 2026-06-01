# SN Test Suite — Capability Taxonomy Reorg (Nexus-driven)

**Why:** `tests/sn/` = 74 files / 1376+ tests with no logical organization →
running the whole tree OOMs (single-process), and there's no visible
"what to test for change X". Reorganize into an **incremental SN-capability
taxonomy** (each tier depends on those below) so we can run focused tiers,
see coverage gaps, and the suite stays memory-safe. The poor file division
IS the symptom. Decided 2026-06-01 (user: "now, before continuing A-phase").

**OOM bridge (DONE):** `pytest-xdist 3.8.0` installed + declared in
`pyproject` `test` deps. Canonical heavy run: `pytest tests/sn -n auto
--dist loadscope` (per-worker process isolation bounds memory + parallelizes).
Not forced via addopts. Smoke-validated: 216 operator tests / 1.68s / -n 4.

**Graph staleness (user hypothesis CONFIRMED):** Nexus graph built 2026-05-17
(~2wk stale) — still shows pre-unification `Inflow/OutflowTraceSpace`, old
`orpheus.sn.{angular,scalar}_flux` paths, and OMITS ~18 newer files. The
`tests`→equation edges + V&V markers are reliable (verified vs source); node
PATHS are not. **MUST `sphinx-build` to refresh the graph before trusting
imports/calls for post-May-17 files** (also validates the A.2/A.3 docstrings).

---

## Capability ladder (tiers; each depends on those below)

- **T0 primitives** — quadrature, octant predicate, face_layout, axis_primitive,
  trace_space, method_space, typed fields/sources, solution container, snmesh.
  No physics.
- **T1 operators** — L (streaming), C (collision), S (scattering),
  F (fission), B (boundary realizer), angular_average, operator_algebra
  (OperatorSum/CAP/InvertibleOperator). Apply / .H / compose / capability sets.
- **T2 sweep** — spatial closure. core (cell_update/balance/diamond/ordinate_scan/
  sweep_cache/dag_walk); 2a slab; 2b curvilinear (pole_angular_closure +
  psi_half_angle_seed); 2c 2-D Cartesian/octant. `matvec ≡ sweep` identity.
- **T3 fixed-source solve** — source_iteration, krylov; 1g→mg, homog→heterog.
- **T4 eigenvalue** — k-eigenvalue.
- **T5 verification gates** — MMS (slab/2d/aniso/curvilinear/hetero),
  convergence order, analytical standoff (L1/L2), B1. Sit ON TOP of their tier.

## Capability DAG (change→run mapping)

```
quadrature {}              ; face_layout {} ; axis_primitive {}
octant_predicate {quadrature}
trace_space {quadrature, face_layout} ; method_space {trace_space, quadrature}
snmesh {axis_primitive, face_layout, quadrature, mesh, materials}
typed_fields {snmesh} ; typed_sources {snmesh} ; solution_container {typed_fields}
angular_average {quadrature}
collision {snmesh, typed_fields} ; fission {snmesh, typed_fields, material_xs}
scattering {harmonic_moment, angular_average, material_xs}
boundary_realizer {angular_average, trace_space, geometry.boundary_ops}
operator_algebra {collision, streaming, scattering, fission}
sweep_core {cell_update, cell_balance, diamond, ordinate_scan, sweep_cache}
dag_walk {snmesh}
slab_sweep(2a) {sweep_core, reduced_operator.slab_streaming, dag_walk}
curvilinear_sweep(2b) {slab_sweep, reduced_operator.{cyl,sph}_streaming,
                       pole_angular_closure, psi_half_angle_seed}
cartesian_2d_sweep(2c) {sweep_core, sweep_graph, reduced_operator}
streaming(L.apply) {snmesh, sweep_core}     ; matvec≡sweep {streaming, sweep_core}
source_iteration {operator_algebra, sweep, boundary_realizer, collision, scattering}
krylov {source_iteration, numerics.iteration.KrylovAcceleration}
fixed_source_1g {source_iteration} ; fixed_source_mg {fixed_source_1g, scattering, fission}
k_eigenvalue {fixed_source_mg, fission, numerics.eigenvalue}
mms_* {fixed_source_mg(geometry)} ; convergence_order {mms_*}
analytical_standoff {k_eigenvalue, cp.solve_cp, trajectory_resolvent_greens}
b1_verification {streaming, collision}
```
Rule: edit module X → run X's capability node + everything that transitively
lists X in depends_on. `quadrature` is the god-leaf (touch it → run ~everything).

## Proposed directory layout

```
tests/sn/{primitives,operators,sweep/{core,slab,curvilinear,cartesian_2d},
          solve/regression,eigenvalue,verification/{mms,analytical}}/
```
(See explorer transcript a6926e2a3a328b3ea for the full per-file destination
table.) Encoding = directories + `@pytest.mark.cap("<capability>")` markers
(composable with existing l0/l1/l2/foundation) + the DAG above surfaced via
`tests._harness.audit` as a capability×level matrix.

## The 4 worst offenders (≈80% of the win — split these first)

1. **test_cartesian.py** — spans T2a+T3+T4+T5 (DD-recurrence unit, fixed-source,
   k-eff, convergence). Split → slab DD unit / slab keff / slab convergence.
2. **test_cylindrical.py** — T2b+T3+T4+T5. **Carries the only `l2` markers in 1-D
   (+CP cross-check) — DO NOT lose them.** Split → cyl coeffs / cyl keff / CP-standoff L2.
3. **test_spherical.py** — mirror of cylindrical (also `l2`, `slow`). Same split.
4. **test_solver_components.py** (41 tests) — grab-bag T0–T3 (recipes/scatter
   harmonics/sweep/solver). Decompose by tier.

## Mis-filings / orphans / gaps (test-architect MUST address)

- **Orphans with NO V&V marker** (invisible to `-m` tiers + audit): `test_boundary_conditions`,
  `test_2d_l2_face_view_unit_source`, `spatial/test_cell_balance_for_streaming`,
  `regression/test_dd_regression`, `test_si_cyl_20cell_nan_regression`. Force a marker on each.
- **Legacy under test** — `SNStreamingOperator` + `build_equation_map*` codec (retirement
  targets per D-I.3/D-J) pinned by test_collision/streaming/phase_c_gates/spatial.apply_matvec.
  Tag transitional; don't enshrine in the permanent operator tier.
- **Duplicated coverage** — `test_phase_c_mms` ≡ `l1_analytical/test_mms_curvilinear_aniso_dd_convergence`
  (same `build_{cyl,sph}_anisotropic_mms_case`). Consolidate.
- **matvec≡sweep scattered across 5 files** (unified_matvec_{slab,cyl,sphere} +
  native_matvec + sweep_vs_apply_consistency) — one `test_matvec_consistency.py` per geometry.
- **`test_quadrature` undersell** — carries `alpha-cylindrical`/`wdd-*` edges (T2b coeff
  territory); verify which tests own those before filing as pure T0.
- **COVERAGE GAPS (V&V findings):** 2-D L1/L2 is THIN (~16 tests vs ~120 for 1-D);
  **NO 2-D k-eigenvalue test exists**; curvilinear k-eff lives ONLY inside the mixed
  cyl/sph files (preserve a standalone `eigenvalue/test_keff_curvilinear.py` in the split);
  no standalone boundary-operator L1 analytical gate (reflective-bc only verified
  incidentally inside solves).

## Sequencing constraints

1. **numerics-investigator is editing `test_spherical.py`/`test_cylindrical.py`/
   `sweep_cache.py`/`sweep.py`** (the curvilinear `sig_t[:, chain_idx]` bug). The reorg
   MOVES/SPLITS those test files → MUST wait for the investigator to land first (collision).
2. **Rebuild Sphinx** after the investigator lands → fresh graph for the test-architect +
   validates A.2/A.3 docstrings.
3. Then dispatch **test-architect** to execute the reorg incrementally, keeping the suite
   green at each move (git mv preserves history; importlib discovery already configured).
4. Validate `pytest tests/sn -n auto --dist loadscope` completes without OOM (the proof).
5. Resume field-typing A.4/A.5 on the reorganized, memory-safe suite.

**Partial reorg already exists:** `tests/sn/{spatial,l1_analytical,regression}/` (16 files).
Absorb, don't ignore.
