---
name: task54-sn-spatial-rename-audit
description: "TRANSIENT pre-rename blast-radius audit for orpheus/sn/spatial (task #54), read-only @ main 6e3ebad0 (2026-07-13). Delete when #54 lands."
metadata:
  type: project
---

# Task #54 — `orpheus/sn/spatial/` rename: pre-rename audit

Read-only audit @ `main` `6e3ebad0` (2026-07-13). Line numbers current at this
HEAD — re-derive via Nexus/grep if drifted. **Delete this file when #54 lands.**

Durable headline (survives drift): the residual package is the sweep's
**per-ordinate kernel substrate** (scan primitive, coefficient caches, angular
closure family, ψ½ ray-march kernels, one dormant pair predicate) — the walk
EXECUTORS live in `sn/loss_representation/` and the ray OPERATORS in
`sn/operators/radial_characteristic.py`. All production imports are
**direct-module-path** (4 files, 7 statements); the package `__init__`
re-export surface has ZERO live importers. Grep was primary evidence: Nexus
collapses internal module-import targets to `py:module:orpheus`, so
module-granularity import edges don't exist in the graph (function-level
`calls` edges were used for symbol callers).

## 1(a) Graph + import census — production importers

| importer (file:line) | what it imports | path form |
|---|---|---|
| `orpheus/sn/mesh/augmented_mesh.py:59-65` | `pole_angular_closure`: IdentityAngularClosure, MorelMontryAngularSweep, PoleAngularClosureBase, default_angular_closure_class, morel_montry_tau_raw_per_level | relative `..spatial.pole_angular_closure` |
| `orpheus/sn/solver.py:63` | `sweep_cache`: CollisionCache, GeometryCoefficients | relative `.spatial.sweep_cache` |
| `orpheus/sn/loss_representation/__init__.py:149-153` | `scan`: `_scanmarch_row`, `_x_scan_faces`, ordinate_scan, ordinate_scan_transpose | relative `..spatial.scan` |
| `orpheus/sn/loss_representation/__init__.py:155` | `sweep_cache`: CollisionCache, GeometryCoefficients | relative |
| `orpheus/sn/loss_representation/__init__.py:3592, 3961` | `pole_angular_closure`: MorelMontryAngularSweep (function-local, x2) | relative |
| `orpheus/sn/operators/radial_characteristic.py:156-161` | `psi_half_angle_seed`: carlson_inward_sweep_from_source, carlson_inward_sweep_transpose, radial_characteristic_forward_residual(+_transpose) | **absolute** `orpheus.sn.spatial.psi_half_angle_seed` |

- `orpheus/sn/__init__.py`: ZERO spatial mentions — the sn package surface does
  not re-export anything from the subpackage.
- Only bare package-level import anywhere: `derivations/diagnostics/diag_276_full_scatter_kernel_ld_trailing_axis.py:41`
  `from orpheus.sn.spatial import LinearDiscontinuous` — **ALREADY BROKEN**
  (ImportError: LD moved to transport.spatial at #272 and is not re-exported).

Symbol callers (Nexus `calls` edges + grep):
- `carlson_inward_sweep_from_source` — production: `RadialCharacteristicOperator.solve`
  (`orpheus/sn/operators/radial_characteristic.py:459`) ONLY; + 11 test fns + 4 diagnostics.
- `ordinate_scan`(+`_transpose`) — production: loss_representation walk bodies
  (`__init__.py:3429, 3503, 3771`; transposes `:3931, 4121`) + scan-internal.
- `default_angular_closure_class` — sole caller `SNMesh._init_core`
  (`orpheus/sn/mesh/augmented_mesh.py:438`).
- `pair_diffusion_limit_consistent` — **ZERO production callers**; 3 test fns in
  `tests/sn/spatial/test_pairing_diffusion_limit.py`. Its docstring self-declares
  "no production call site today" (wiring at mesh construction pending).
- Closure classes — constructed in augmented_mesh; function-local imports in
  loss_representation (:3592/:3961); every other production mention is
  docstring/comment: `transport/mesh/axis.py:12`,
  `numerics/spaces/angular_trace_space.py:105`, `sn/solver.py:2681-2682`,
  `transport/spatial/cell_balance.py:146,153`, `transport/spatial/scheme.py:256,270,279`,
  `geometry/reduced_operator.py:225,326,339,343,521,697`.

## 1(b) Text census (regex `sn\.spatial([^_a-zA-Z]|$)` — excludes the
`sn.spatial_shape` fixture-attribute FALSE POSITIVES that inflate a naive grep)

- **orpheus/ outside the package**: 40 dotted + 3 slashed lines across 11 files —
  transport/spatial 11 (linear_discontinuous 2, diamond 1, `__init__` 1,
  cell_balance 2, scheme 5), numerics 4 (moment_layout 3 historical narrative,
  radial_characteristic_space 1 role), sn/mesh/augmented_mesh 3, sn/solver 1
  (comment), sn/operators/radial_characteristic 11, sn/loss_representation 7
  dotted + 2 slashed (+ sweep_graph.py:22 a `:doc:` ref), geometry/reduced_operator 3.
  All docstring roles/comments except the 7 import statements above.
- **Package-internal self-refs**: 16 dotted + 2 slashed in the 6 files (incl.
  `__init__.py` module map + R9 NOTE).
- **tests/**: 43 dotted lines across 20 files = **20 import statements in 15
  files** (15 top-level + 5 function-local: `_test_helpers.py:593`,
  `test_assembly_mode.py:693`, `test_cyl_direct_seed_fold.py:176`,
  `test_phase_c_gates.py:864`, `test_si_cyl_20cell_nan_regression.py:167` — the
  last verified a plain late import, NOT a subprocess string) + 23
  docstring/comment lines. Prior estimate "17 test imports" → CORRECTED: 20
  statements / 15 files. Plus 6 slashed lines in test docstrings/comments and 5
  in `tests/_mutation/` (see 1(c)).
- **docs/**: 155 dotted lines + 11 slashed = 166 lines, 5 theory pages ONLY
  (discrete_ordinates 130+7, loss_representations 13, index_convention 6+3,
  operator_algebra 3+1, boundary_conditions 3). Prior "~150 docs refs" ≈ verified.
  Exact role count: **147 role occurrences** — 137 target `orpheus.sn.spatial.*`,
  10 target `tests.sn.spatial.*` (ALL 10 already dangling, see 1(d)); + 8
  plain-literal lines (historical `boundary_face_flux` blocks at
  discrete_ordinates.rst:6593-6596, 7083-7086 — dev-history, leave as history).
- **docs/api/**: ZERO — no automodule/toctree entry for any sn.spatial module
  (pre-existing API-doc gap; the package is not autodoc-rendered). The rename's
  doc surface is theory pages + rendered production docstrings only.
- **derivations/**: 14 dotted lines / 9 diagnostics + 3 slashed; real imports in
  ~6 (diag_phase_g_step2_* family, diag_si_cyl_20cell_nan_step5, diag_195_probe5,
  diag_276_full_scatter_kernel_ld_trailing_axis (broken)).
- **examples/, scripts/, tools/, student_resources/**: ZERO.

## 1(c) String-literal / config

- **ZERO** quoted `sn.spatial` literals in any `.py` (no monkeypatch.setattr
  paths, no importlib strings, no `__module__` assertions, no subprocess-worker
  import strings).
- `tests/_mutation/diamond_spike.toml:10` + `diamond_sentinels.toml:5`
  `module-path = "orpheus/sn/spatial/diamond.py"` + `README.md:20,46,63` —
  **STALE SINCE #272** (diamond.py lives at `orpheus/transport/spatial/diamond.py`);
  pre-existing breakage, fix-in-passing candidate.
- `pyproject.toml`: no spatial mentions ([tool.pyright] is venv config only).
  `pyrightconfig.json`: none. `tests/_harness/pyright_baseline.json`: keys are
  top-level orpheus subpackages (currently `{diffusion, transport}`; `sn` not
  even listed) — **the rename cannot move the baseline**. All 12 conftests: zero.

`.claude/` advisory tier (do NOT churn; a post-rename reader trips on):
active plans `stencil_assembly_dsa_roadmap.md` (13), `a3_solve_transpose{,_verification}.md`
(2+4), `facefield_codim1_design.md`, `pyright_signal_cleanup.md`,
`sn_package_reorganization.md`, `field_typed_operator_algebra_campaign.md`,
`coupled_block_operator_campaign.md`, `.claude/lessons.md`,
`.claude/skills/vv-principles/error_catalog.md` (**42 mentions** — ERR entries
cite sn/spatial paths), ~60 agent-memory files, plans/archive/*, plus the
`worktrees/nexus-workspace-wiring/` checkout's own copies (separate tree).

## 1(d) Pre-existing dangling refs (NOT rename-caused; fix while in there)

1. All 10 `tests.sn.spatial.*` doc roles dangle: test_ld_ubld_{primitive,symbolic},
   test_linear_discontinuous, test_affine_closure (→ now `tests/transport/spatial/`),
   test_sweep_vs_apply_consistency (→ `tests/sn/sweep/core/`),
   test_psi_half_angle_seed x3 (→ `tests/sn/sweep/curvilinear/`).
2. Slashed doc paths: index_convention.rst:264,1632 (test_ordinate_scan.py →
   sweep/core), :649 + operator_algebra.rst:5492 (test_streaming_equilibrium_
   curvilinear.py → sweep/curvilinear). Current: discrete_ordinates.rst:3826
   (test_spatial_moment_field_space.py exists).
3. `orpheus/transport/spatial/scheme.py:726` role → tests.sn.spatial.test_scheme_
   reaction_rate_contract (→ now `tests/transport/spatial/`).
4. `tests/transport/spatial/test_ld_ubld_*.py` (6 lines) +
   `tests/derivations/test_sn_mms_ld_2d_stress_symbolic.py:309` roles target
   `orpheus.sn.spatial.LinearDiscontinuous` — class moved at #272.
5. `orpheus/sn/loss_representation/sweep_graph.py:22` `:doc:`../sn/spatial/scheme``
   — no such doc.
6. Stale test-path comments: `tests/sn/sweep/core/test_ordinate_scan.py:479`,
   `test_unified_sweep_dispatch.py:36` ("tests/sn/spatial/test_diamond.py" →
   now `tests/sn/sweep/core/test_diamond.py`).

## 2. Name-candidate evidence (no pick made)

- **Walk executors live in `sn/loss_representation/`**: `_WalkLeg` (:2272),
  `_reverse_traversal` (:2300), `_OneDimScanWalk` (:2316), `_dag_legs` (:2350),
  `_loop_walk` (:2400) in `__init__.py`; `SweepDependencyGraph`
  (sweep_graph.py:375); `SweepSchedule` (sweep_schedule.py:90). `SweepOperator` =
  `orpheus/sn/operators/sweep_operator.py:83`; `dag_walk` on SNMesh
  (augmented_mesh.py:1081). → naming the residual `walk` would be dishonest.
- **What the 5 modules are**: the kernels the walks consume — `scan.py` (affine
  recurrence primitive + transpose), `sweep_cache.py` (two-stratum coefficient
  cache), `pole_angular_closure.py` (angular redistribution family + τ
  producers), `psi_half_angle_seed.py` (ψ½ ray-march kernels wrapped by the
  A_BB/A_AB operators), `pairing.py` (dormant pair predicate).
- **`orpheus/sn/sweep/` collision**: no production module/package named `sweep`
  (historical `sweep.py` dissolved at S6.4(f), per docs/api/discrete_ordinates.rst:9-11).
  But `sweep_`-prefixed modules exist in THREE sibling packages
  (operators/sweep_operator.py, spatial/sweep_cache.py,
  loss_representation/sweep_{graph,schedule}.py), and `tests/sn/sweep/` is
  BROADER than this package — its core/ bucket also tests loss_representation
  walks (test_dag_walk, test_one_dim_loop_walk, test_sweep_graph*,
  test_sweep_schedule*, test_sweep_regression, test_wavefront_*) and
  transport/spatial schemes (test_diamond, test_cell_kernel_batch,
  test_cell_balance_for_streaming, test_discretization_scheme_protocol). A
  production `sn/sweep/` holding ONLY kernels would be narrower than the test
  bucket of the same name.
- **Sibling noun families**: `orpheus/sn/` = boundary, loss_representation,
  mesh, operators (+ solver.py, solution.py, coupled_system.py); `orpheus/transport/`
  = displacements, fields, frames, mesh, operators, residuals, source_sinks,
  sources, spatial. Role-keyed plural/mass nouns.
- **`tests/sn/spatial/` table** (the #272 promotion did not re-mirror):

| test file | production target | verdict |
|---|---|---|
| test_pairing_diffusion_limit.py | sn/spatial pairing + pole_angular_closure (+ transport/spatial traits) | tests THIS package |
| test_ordinate_scan_reset.py | sn/spatial scan (+ solver E2E) | tests THIS package |
| test_ld_slope_frame.py | transport/spatial LD + sn solve pipeline | NO sn/spatial import — mirror-move flag |
| test_moment_axis_predicates.py | SNMesh + transport/spatial DD/LD + AngularFlux | NO sn/spatial import — mirror-move flag |
| test_spatial_moment_field_space.py | numerics SpatialMomentSpace + transport fields/schemes | NO sn/spatial import — mirror-move flag |

## 3. Flags (flag only)

1. `pairing.py` production-dormant + cross-layer (conjunction of a
   transport.spatial scheme trait and an sn closure trait; both axes meet on
   SNMesh, which holds scheme + pole_angular_closure). Options surfaced: ride
   the rename / re-home near the mesh pair-site / production-wire at mesh
   construction (docstring says pending). User's call.
2. 3 of 5 `tests/sn/spatial/` files test only transport.spatial/numerics (table).
3. Mutation TOMLs point at a nonexistent module-path (1(c)) — broken since #272.
4. `diag_276_full_scatter_kernel_ld_trailing_axis.py:41` ImportError (1(a)).
5. The 1(d) dangling-ref backlog is test-tree-reorg drift, not rename debt —
   budget it separately in the rename session's doc sweep.
