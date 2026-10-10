# P1 step (e): dependency audit of the old trajectory-resolvent family (`retirement-audit` G.24)

Explorer, 2026-10-10, at HEAD `62468273` (step (d) `d9425977` is an ancestor). Read only.

**Destination.** This file is meant for `.claude/plans/characteristic_reference_architecture_dependency_audit.md`. The explorer's write-scope hook refused that path, so the orchestrator copies it there.

The probes and the full listings are beside this file:
- `census_e.py` / `census_e.log`: the AST importer census and the per-test taint;
- `markers_e.py` / `markers_e.tsv`: the markers on every test;
- `roster_vs_tree.tsv`: the spec's 298 rows joined to the tree;
- `labels_e.txt`, `symbolic_rows.txt`, `shared_helpers.py`;
- `collect_all.log`: `pytest --collect-only`, 17 351 ids;
- `claude_memory_hits.txt`.

Every `file:line` below comes from a `grep -n` of the quoted anchor.

## 0. The premise, corrected

The step-(e) sketch (P1 sketch, item 9) says "retire the old family". The spec rules otherwise for one part: B7 and the 61 KEEP SymPy rows keep the SymPy origins as the Branch-1 algebra of record. So `trajectory_resolvent/origins/specular/` is **not** deleted:
- 18 of its 23 `derive_*` functions stay;
- 5 retire with the 26 duplicate rows.

Its home after the package goes needs a ruling (NEEDS 1).

**AST check (`symbolic_rows.txt`):**
- The 26 RETIRE rows reach exactly these 5 functions:
  - `derive_rank2_resolvent_{annulus,hollow_sphere}`, 7 + 7 rows;
  - `derive_operator_constant_trial_closed_{cylinder,slab}`, 4 + 4;
  - `derive_alpha_zero_kernel_reduction_annulus`, 4.
- No KEEP row reaches any of the 5.
- No kept `derive_*` calls a retiring one.
- `origins/` imports nothing from the numeric family. It is a leaf.

## 1. The module set (AST closure)

**The family's own modules: 20 under the package** (11 792 lines).
- **Numeric, deleted: 12 modules, 8 054 lines.**
  - `__init__`, `billiard`, `chord_oracle`, `power_iteration`, `reference`, `variant_alpha_core`;
  - `greens_function{,_annulus,_cylinder,_hollow_sphere,_slab,_slab_asymmetric}`.
- **Origins, KEPT, re-homed by ruling: 8 modules, 3 738 lines.** Of these, about 594 lines retire: the 5 `derive_*` functions above.

**Siblings the family owns outside the package: none.**
- The family imports 14 first-party modules outside itself: `data.cells`, `data.macro_xs.mixture`, `derivations.common.{reference_body, solution_types, eigenvalue}`, `geometry`, `geometry.coord`, `numerics.{content, mesh_free_function, observable, question, traced_memo}`, `reference.{reading, solution}` and `specification.specification`.
- Each of these has at least 1 importer in production outside the family (`git grep`; positive control `solution_types`, 3 importers).
- `numerics.traced_memo` drops to **1** production client, `characteristic.reference`.

**Production callers outside the family: 0.**
- AST over 1 285 tracked `.py` files outside `scratch/`: no `orpheus/`, `tools/`, `examples/` or `docs/` file imports the family.
- Positive controls (each shape found in `tests/`):
  - a late import inside a method (`cross_method/adapters.py:208`);
  - a late import inside a function (`test_kernel_corroboration.py:153`);
  - a module alias (`_characteristic_ladders.py:491`);
  - a relative import (`origins/specular/__init__`).
- Nexus agrees on the module level (`impact` on `greens_function`: 11 importers, all tests). It misses the late imports, as expected.

**The old test-side helpers, deleted:**
- `_trajectory_resolvent_ladders.py` (162 lines);
- `_trajectory_resolvent_aba.py` (55);
- `_trajectory_resolvent_api.py` (113): an importlib-string reader, consumed only by `test_trajectory_resolvent_reference.py`.

## 2. Symbol table (external references, by AST)

Symbols with 0 external references are omitted here and listed in `census_e.log`, `SYM` rows. None of them has a production reader.

| defining module | symbol | ext. refs | readers (all tests) | disposition |
|---|---|---|---|---|
| billiard | `Billiard` | 30 | `test_polymorphism`, `test_reference_body`, `test_trajectory_resolvent_billiard` | DELETE |
| chord_oracle | `ChordOracle` + 7 classes | 2–4 each | `test_trajectory_resolvent_{chord_oracle,reference,regionwise_source}` | DELETE |
| chord_oracle | `_trajectory_segments_oracle`, `_region_at_radius_{oracle,cyl}` | 2, 1, 1 | `test_kernel_corroboration.py` (lines 153, 186, 221) | DELETE; edit P0 rows (§5) |
| chord_oracle | `_regionwise_cubic_spline` | 3 | Garcia file, `regionwise_source` | DELETE |
| greens_function* | the 15 `solve_greens_function_*` | 1–21 | 18 test files plus `cross_method/adapters.py` (208, 278, 344) plus `_trajectory_resolvent_ladders:130` | DELETE |
| greens_function* | the 7 k_inf-guess sites | (internal) | `k_eff = float(initial_k)`: cylinder 539/858, hollow 441, sphere 538/1026, slab-asym 438, annulus 526 | DELETE |
| power_iteration | `power_iterate_variant_alpha`, `PowerIterationResult` | 6, 1 | `test_trajectory_resolvent_power_iterate` | DELETE |
| variant_alpha_core | `apply_variant_alpha_closure`, `compute_resolvent_T` | 8, 2 | `test_peierls_variant_alpha_core` | DELETE |
| reference | `trajectory_resolvent_reference` | 4 | `_trajectory_resolvent_aba`, `test_traced_memo_clients` (47, 277) | DELETE |
| reference | `TrajectoryResolventDerivation`, `_RAYS`, `_CylinderRays` | 0 | **by strings only**: `_traced_memo_api.py:152`, `TREE_WITNESSES` (`test_traced_memo_clients.py:171`), `_trajectory_resolvent_api.MODULE` | DELETE |
| origins.specular | 18 kept `derive_*` | 1–7 | the six `*_symbolic.py` files and `test_ambient_state.py:144` (a string) | **KEEP** (re-home per NEEDS 1) |
| origins.specular | 5 duplicate `derive_*` | 4–7 | the 26 retiring rows | DELETE |

**SHARED, do NOT delete.**
- `origins/` (18 kept functions).
- `tests/gates/derivations/_ladder_rules.py`. `richardson_error`, `alternating_error` and `ceil_one_significant_figure` keep only the foundation rows in `test_crosscheck_harness.py` as consumers (NEEDS 6).
- `tests/gates/_corroboration.py`. `import_closure` (225) and `assert_closure_independent` (262) are orphaned: their only consumer is the deleted G file. Delete them (G.24) or keep them for a future corroboration (NEEDS 6).
- `_characteristic_ladders.py`. Remove `_old_reference` (481) and its two `_LADDERS` keys (`"old-sphere-reference"`, 524).
- `_aba_reference.py`.
- The Garcia tables (§5).
- `reference_body`:
  - `HollowBody` and `LayeredBody` lose their only consumer outside the module (`Billiard`), but they remain arms of the total classifier `reference_body()`;
  - `describe`, `specular_albedo(s)` and `require_vacuum` have other consumers.
- `solution_types` (`CriticalSolution`/`FluxSolution`: F_N, Case, Galerkin).
- `kinf_homogeneous` (10 importers).
- `xs_library.get_mixture`.

## 3. Surfaces (columns of G.24)

**Docs**
- **`.. implements` with `:by:` naming a family symbol: 0.** Positive control: 508 `orpheus.` names in the 4 lines after the `:by:` lines elsewhere.
- **`automodule` of a family module: 0.** Positive control: 26 `automodule:: orpheus.derivations` exist. So the **130** rendered Python-domain xrefs into the family are already plain text:
  - 111 in `references/trajectory_resolvent.rst`;
  - 7 in `peierls.rst`;
  - 4 in `curvilinear_one_group.rst`;
  - 3 in `error_catalog.rst`;
  - 2 in `reference_cache.rst`;
  - 1 each in `index`, `cross_method` and `reference_solutions`.
  - `sphinx -n` sees no change. **Grep is the only gate** (A.2).
- **`:label:`s defined inside family docstrings:** `peierls-greens-annulus-tau-step` (origins annulus 83) and `bickley-naylor-Ki3` (origins cylinder 332). They are unrendered, and no doc citer exists.
- **Nexus `implements` edges** from the family: ≥ 200 (truncated), all inferred page-sharing. None is declared, so none reddens.
- **`dead_references` today: 0.** Re-run after the delete. The docstring sites it will see are:
  - `solution_types.py` 10, 12, 17, 56–57, 94, 173–187;
  - `galerkin_spectral/basis_space.py` 6, 54, 67, 844;
  - `singular_eigenfunction/spectrum.py` 5, 80, 453;
  - `continuous/__init__.py` 22;
  - `peierls_nystrom/__init__.py` 13;
  - `reference_body.py` 4–5;
  - `geometry/structured_geometry.py` 55;
  - `geometry/boundary/_tag.py` 101;
  - `fn_method/cylinder/__init__.py` 11;
  - `peierls_nystrom/cases.py` 74.

**Strings that turn a code commit red on their own**
- `tests/gates/withdrawal_506_placement.txt` 95–97. The three #506 rows (`slab_solver::test_alpha_zero_vacuum_agrees_with_nystrom_slab`, `xverif::test_b5_phase4_*`) must leave the list, or M4 reds.
- `test_withdrawal.py:382` `assert len(files) == 27` becomes 25, and the prose at lines 33 and 374 ("26 files") becomes 24. This is a slow M5 row.
- `test_error_catalogue_reconciles.py` (CI `gates.yml`):
  - arm 2 reds on ERR-034, ERR-035 and ERR-091 (§6);
  - arm 5 reds on the catalogue's path citations at `error_catalog.rst` 2576, 2662, 8696, 8700, 8703 and 8782;
  - arm 3 needs `error_index.md` regenerated (ERR-090 4→1, ERR-091 1→0).
- `test_content_identity.py`:
  - line 64 imports `ROSTER` from the deleted `test_trajectory_resolvent_reference`; drop it and its term at line 71;
  - line 380 `_PACKAGES` names `trajectory_resolvent`; keep it only if `origins/` stays there.
- `_traced_memo_api.py` 150–152 (`CLIENTS`) and `TREE_WITNESSES` (`test_traced_memo_clients.py:171`, 4 family file paths).
- `cross_method`: the `ADAPTERS` keys (`adapters.py` 526–528) and 8 tolerance keys in `cases.py` (112, 139, 162, 184, 228, 256, 279, 588). `closed-sphere-1G-fuelA-tauR2.5` has ONLY the TR key, so `test_case_has_at_least_one_tolerance` reds (NEEDS 4).
- `test_ambient_state.py:144`, the subprocess import of `origins…greens_function_slab`: re-path it per NEEDS 1. Nothing is removed.

**`.claude/`**
- AGENT.md and rules: 0.
- Skills: `algebra-of-record/SKILL.md` 607 and 749 (present-tense list); `vv-principles/error_index.md` 140–141 (generated).
- Agent memory: 68 files, 337 lines, across 10 agents (`claude_memory_hits.txt`; method-implementer 28 files). Topic files describing deleted code in the present tense include `r2_billiard_class.md`, `r3_chord_oracle.md` and `r1_power_iterate_driver.md`; sort by tense (E.19).
- Inventories: `implements_declaration_inventory.md:40`, `declare_moc.md`.

## 4. Order (one commit)

1. **Leaves** (zero callers):
   - the 22 test files deleted whole (§7) and the 3 old helpers;
   - the RETIRE rows and Billiard parameters in edited files (§5);
   - `_old_reference`;
   - the corroboration file.
2. `reference.py`, `billiard.py`, the package `__init__` re-exports.
3. `greens_function_cylinder` (imports `greens_function`) and `greens_function_slab` (imports `_slab_asymmetric`).
4. `greens_function`, `_annulus`, `_hollow_sphere`, `_slab_asymmetric`.
5. `chord_oracle`.
6. `power_iteration` and `variant_alpha_core`.
7. The package directory, minus `origins/` (re-homed in the same commit, NEEDS 1), plus the 5 duplicate `derive_*` functions.

## 5. Files that survive with edits

| file | exact edit |
|---|---|
| `geometry/test_kernel_corroboration.py` | drop `"variant_alpha"` from the parametrize at 143 and its arm at 151–156; delete `test_the_backward_segments_agree_with_variant_alpha` (it has no other multi-region segment spelling; NEEDS 5); in the locator row, remove `_region_at_radius_{oracle,cyl}` from the import at 221, the `_assert_independent` list and `old_spellings`, and "five inner-owns spellings" becomes three (docstring and assertion) |
| `derivations/test_characteristic_reading.py:60` | re-home `GARCIA_2021_CASE1_{R,PHI_PRINTED,PHI,ROUNDING}` (Garcia file 92–108) to a helper (or inline them); `GARCIA_TO_VARIANT_ALPHA_FACTOR` and `_RADII` retire; at 1049 "Today's family's bands" goes to the past tense |
| `numerics/test_traced_memo_clients.py`, `_traced_memo_api.py` | `CLIENTS` loses three keys; M4.1b–M4.8 (10 ids) are bound to the old solve child (NEEDS 3); M4.1, M4.9 and M4.10 stay |
| `numerics/test_content_identity.py` | 64 and 71 (`DERIVATION_ROSTER`); 380 per NEEDS 1 |
| `test_ambient_state.py:144` | re-path the import; "two origins derivations" holds |
| `cross_method/{adapters,cases,test_eigenvalue,test_polymorphism,protocol,__init__}` | 3 adapters, 3 RETIRE TAUT rows, 5 KEEP rows (15 ids), 8 tolerance keys, 1 TR-only case (NEEDS 4) |
| `derivations/test_reference_body.py` | 12 ids go: `Billiard` in `_REFUSALS` (178), `test_billiard_routes_*` (4), `reads_the_body_material`, `two_surface`, `two_laws`, `white_law[Billiard]`, `reports_every_group` (ERR-091) |
| `test_fn_sood2003_{slab,sphere}_xverif.py`, `test_singular_eigenfunction_cylinder_xverif.py` | KEEP X rows, 5 functions and 6 ids: re-point the old side to `characteristic_reference` (a resolution and a re-measured band each) |
| the six `*_symbolic.py` | delete the 26 rows; re-path imports per NEEDS 1; move 4 markers (§6) |
| `test_withdrawal.py`, `withdrawal_506_placement.txt` | §3 |
| `_characteristic_ladders.py`, `_corroboration.py`, `_ladder_rules.py:7` | §2 |
| `docs/theory/verification/error_catalog.rst` | path citations (§3), ERR status (§6) |

## 6. Markers and ERR entries

**ERR entries whose ONLY catchers retire** (`markers_e.tsv`; module, class and function marks; positive controls: 47 module and 99 class `pytestmark` shapes).
- **ERR-034**: `slab_asymmetric_solver::test_method_of_images_reflective_vacuum_equals_double_vacuum`. Spec successor B4c is not built. NEEDS 2.
- **ERR-035**: `…::test_rank1_path_now_agrees_with_rank2_via_delegation_after_ERR035_fix`. B6 is not built as a marked row. `test_characteristic_closure.py:502` names "the ERR-035 denominator" as a first red but carries no `catches`. NEEDS 2.
- **ERR-091**: `test_reference_body::test_billiard_sphere_mr_fixed_source_reports_every_group`. E7 is not built as a marked row. NEEDS 2.
- **ERR-090** keeps 1 catcher, `test_characteristic_reading.py:1037`, once the Garcia tables are re-homed. Its other 3 catchers are in deleted files.

**`verifies` migration**
- **From the 26 retiring SymPy rows** (C.13; 5 markers; the twins were found by AST):
  - `annulus-through-rank2` and `hollow-sph-through-rank2` go to `slab_asymmetric_symbolic::test_v_alpha2_slab_asym_determinant_canonical_form` (170);
  - `cylinder-trajectory` and `slab-trajectory` go to `test_peierls_greens_function_symbolic::test_v_alpha1_surface_fixed_point_solves_to_q_over_sigma_t` (73);
  - `annulus-architecture` stays verified by kept rows.
- **Labels left with no live verifier** once the deleted files go (`labels_e.txt`): `peierls-greens-` + `{annulus-impact-parameter-partition, cylinder-T, cylinder-architecture, cylinder-mr-interface-continuity, cylinder-mr-kinf, cylinder-mr-quadrature-convergence, cylinder-mr-trajectory-segments, cylinder-mr-wm72-vacuum, hollow-sph-impact-parameter-partition, mr-regionwise-source, slab-asym-method-of-images}`. The new tests verify only `characteristic-*` and `geometry-*` labels. **None of the spec's planned moves** (to D5, D9, C3, C7, B4c) has landed.

## 7. Test triage

The spec roster, re-checked against the tree:
- all 298 rows exist by name;
- 3 are class-qualified in the roster (`TestTheLawsAreServed.`) and module-level in the tree. This is the only drift.
- AST taint: 265 of the 298 rows read the family directly.
- Of the other 33:
  - 13 are SN rows that step (d) re-pointed;
  - 13 read it through strings or across files (`cross_method` 5, traced-memo 3, `_trajectory_resolvent_api` 2, `test_ambient_state` 1, `test_polymorphism` FN rows 2, which do not read it at all);
  - 4 are production-primitive or Nyström rows (`*_via_production_primitive*` ×3, `test_b5_phase4_rank1_equals_white_hebert`).
- 17 tainted rows are not in the roster: the 10 G rows, 4 P0 kernel rows, the Garcia E1 row and 2 Billiard parameter rows.

**Disposition of the roster rows by file**

| file disposition | KEEP | RE-POSE | RETIRE |
|---|---|---|---|
| 22 files deleted whole | 37 | 86 | 44 |
| edited files | 85 (61 SymPy, 13 SN done, 11 to re-point) | 13 | 33 |

**MIGRATE with NO successor in the tree** (each needs re-pointing in (e), a successor, or a ruling; NEEDS 2):
- KEEP rows in deleted files:
  - WM-72 (`cylinder_mr_xverif`);
  - Sood cylinder (`xverif_sood2003`);
  - PS-1982 (`xverif_ps1982`);
  - method of images ×2 (ERR-034);
  - the two `R_in → 0` limits;
  - interface continuity;
  - the Branch-1/2 ancestor;
  - the 3 #506 Nyström rows;
  - `r7b2_6` indicator/additive.
- RE-POSE rows:
  - D8 ordering theorem, 20 rows (no row found);
  - D9 self-convergence, 22 rows (never built; the ladders module holds steps, not a gate);
  - F6 traced memo, 8 (only `test_c2`–`c4` exist);
  - B6, C9 (3), C12, C6b and E7.
- D5's floor is built (`test_characteristic_system.py:369`, 9 closed fixtures) **without an annulus**. 4 of the 20 VK rows are annulus rows.

**Successors present:**
- D5 (sphere, slab, hollow, cylinder);
- C7, E1 (Garcia), F1–F5 and A1;
- Sood sphere k_inf and criticality (`test_characteristic_system.py` 400, 537, 584).

## 8. Size

| | files | lines | ids |
|---|---|---|---|
| package numeric modules | 12 | 8 054 | — |
| duplicate `derive_*` in origins | (5 fns) | ~594 | — |
| test files deleted whole | 22 | 9 442 | 246 |
| old test helpers | 3 | 330 | 0 |
| retiring rows in edited files | — | — | 26 SymPy + 3 polymorphism + 12 reference_body + 2 kernel = 43 |
| **total** | **37 files** | **~18 420** | **289** |

Collection: 17 351 before. The floor after is **17 062**, with every KEEP re-pointed and nothing added. It falls to 17 037 if NEEDS 3 and 4 retire the traced-memo (10) and cross-method (15) ids.
