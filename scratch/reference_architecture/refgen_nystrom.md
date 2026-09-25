# Nyström blast set — detail (explorer, 2026-09-24, HEAD 901f64ca + dirty plan only)

Scripts (re-runnable, all in this scratchpad): imp.py (AST importer census), cls.py (per-file classification),
pertest.py (per-test-function attribution, module-local call closure), mkdump.py (pytest plugin dumping resolved markers),
markers_all.json (12 495 collected ids, whole tests/gates), pertest.json, nys_durations.txt (join with run 35940553034 artifact dur/test_durations.json).

## 1. What "a Nyström reference" is
Package orpheus/derivations/continuous/peierls_nystrom/ (9 825 lines, 12 .py + origins/specular/*).
- NYSTRÖM SOLVE machinery: geometry.py build_volume_kernel(_adaptive), K_vol_element_adaptive, _build_full_K_per_group,
  build_white_bc_correction(_rank_n), build_closure_operator + _build_*_op/_PG, BoundaryClosureOperator, ClosureRecipe,
  PeierlsSolution, solve_peierls_1g/_mg, per_observer/per_surface_centred_angular_assembly; slab.py solve_peierls_eigenvalue,
  _build_kernel_matrix, _build_system_matrices; cylinder/sphere/slab _build_peierls_*_case; cases.py (registry: 13 lazy builders,
  continuous_case_builders/continuous_cases/build_two_surface_case/build_one_surface_compact_case); ps1982_reference.py
  (product-integration eigen-solve; "structurally independent" PS-1982 reference — Nyström-family, lives in the package).
- NOT Nyström (closed-form / geometry primitives in the same package): compute_P_ss_*, compute_T_specular_*, compute_P_esc*,
  compute_G_bc*, reflection_*, composite_gl_r, CurvilinearGeometry, naming.py, reference.py (analytical uniform-source identities),
  origins/* (SymPy derivations).  geometry.py mixes both halves in one 6 545-line module.
- Outside the package: fn_method/peierls_atkinson_nystrom.py (Atkinson product-Nyström, FN flux reconstruction default
  "atkinson_nystrom"; ERR-036 catcher test_atkinson_product_nystrom.py) — boundary case, NOT imported by/importing peierls_nystrom.
  trajectory_resolvent/*: 0 imports of peierls_nystrom (AST), not Nyström (its Nyström-sampling prototype was abandoned,
  trajectory_resolvent.rst:380). singular_eigenfunction/cylinder/one_group.py: Nyström prototype replaced by Mitsis-WM.
  flat_source_cp/*: 0 imports (docstring refs only; Nexus depth-2 edges are docstring `references`).

## 2. Production importers
AST pass (imp.py) over orpheus/ tests/ tools/: 52 files hit; 3 non-test non-package hits are DOCSTRING strings
(shifted_legendre.py, discrete/sn/angular_differencing.py, tools/verification/generate_capability_matrices.py).
=> 0 of N production modules outside peierls_nystrom import it. Registry reaches it by pkgutil.walk_packages
(reference_values.py:183) + continuous_case_builders (cases.py:438). Sphinx: capability_rows() (cases.py:487) at builder-inited
via docs/conf.py:149 — static metadata, does not solve.

## 3. Test consumers (per file; n = collected cases; touch = test functions whose module-local closure reaches a solve symbol)
See classes.txt / nys_durations.txt. Summary:
S1 (Nyström-solve files): 28 files, 499 cases, 57 slow; 314 cases attributed Nyström-touching (55 slow); 562 of 901 runner-min.
S2 (import peierls_nystrom primitives only): 12 files (+ rank_n_protocol moved to S1), 196 cases, 3.5 min.
Non-importing peierls-named files: 19 (18 trajectory_resolvent Green's-function + Atkinson): out of scope.

## 4. Markers at risk (per-test attribution; "keep" = cases outside the Nyström-touching set)
verifies lost: flat-source (collision_probability.rst:204, ONLY verifier cp/test_peierls_flux.py), peierls-mg-operator,
peierls-vacuum-bc-{cylinder,flux,slab,sphere}, peierls-white-bc (peierls.rst:662).
Nearly lost: peierls-equation (29 touch / 4 keep, keepers are S1-file primitives), hebert-3-323 (3/1).
catches lost: ERR-027, 028, 029, 030, 032, 063 (6 of 88 entries; `nexus errors` today: 88 entries, 0 uncaught).
Not at risk: ERR-031 (cp/test_cylinder_pss, primitive), ERR-033 (slab aggregate primitive), ERR-034/035 (trajectory resolvent), ERR-036 (Atkinson).
Reconciler tests/gates/test_error_catalogue_reconciles.py is a TEXT check over tests/: a skip keeps it green; a delete/move reds it.

## 5. Production-SUT evidence resting on Nyström
cp/test_peierls_flux.py (CP slab 2G2R flux vs registry peierls_slab_2eg_2rg; 50 min, TIMEOUT on runner) — only verifier of flat-source.
cp/test_peierls_cylinder_flux.py, cp/test_peierls_sphere_flux.py (CP 1G homogeneous white-BC R=10 MFP k and flux vs solve built inline).
Stand-ins [R]: homogeneous white-BC exact solution is flat flux, k = k_inf (closed form) for cyl/sph; slab 2G2R has no
slab multi-region trajectory resolvent (only slab, slab_mg, slab_asymmetric) -> no measured stand-in.

## 6. Hold-out mechanisms existing
pyproject markers: slow, l0-l3, foundation, verifies, catches, rests_on (declarative, no enforcing plugin found), regression, cap, sentinel.
No default -m deselection (addopts = --import-mode=importlib only). CI gates.yml runs 7 fixed files only; test-durations.yml takes -m from input.
Per-directory conftest precedent: tests/gates/sn/*/conftest.py stamp cap() markers via pytest_collection_modifyitems.
ORPHEUS_* env: ORPHEUS_SLAB_VIA_E1 (route selector, cases.py), ORPHEUS_B32 (mutation battery), ORPHEUS_ROOT; none gates tests.
No existing quarantine marker.

## 7. Docs citing Nyström as evidence
summary.rst:56-57 credits CP with "the Peierls–Nyström heterogeneous reference (30+ digit Nyström collocation)".
escape_probability.rst:76 "production-grade Peierls reference". cross_method.rst:21, reference_solutions.rst:646,
peierls_nystrom.rst (289 hits; 43 vv-status all `documented`), peierls.rst (116), matrix.rst (generated),
_peierls_nystrom_capability_matrix.inc.rst (generated from capability_rows), error_catalog.rst ERR-027..033, 063.
.claude: 118 tracked files mention Nyström (6 skills).
