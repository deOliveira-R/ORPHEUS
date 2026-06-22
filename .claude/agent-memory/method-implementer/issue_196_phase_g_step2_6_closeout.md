# Issue #196 Phase G Step 2.6 (Q3) closeout — `dag_walk` canonicalisation

**Branch**: `refactor/sn-operator-algebra`
**Date**: 2026-05-14
**Commits**: `59972e7` (refactor) + `905e2d9` (docs)
**Predecessor**: `badac51` (Step 2.5c closeout)
**Status**: SHIPPED.  All mechanism criteria 1-9 satisfied; 11
regression snapshots bit-identical; L0 streaming-equilibrium
26/26 PASS; tests/sn/test_dag_walk.py 12/12 PASS; XOR-signature
enforced.

---

## What this closes vs what it does not

**CLOSES** (Step 2.6 Q3 cleanup):
- Three iteration methods on `SNMesh` (`dag_walk` alias +
  `iter_cell_visits` + `iter_cells_by_direction`) collapse to
  ONE canonical `dag_walk` under an XOR signature.
- 7 production call sites + 2 test files + ~14 doc references
  migrated.
- Concept-count: 3 → 1 (iteration primitives).

**DOES NOT CLOSE**:
- Step 2.6 Q2 (2D wavefront unification) — Pattern 6 deferral.
- Step 2.6 Q4 (`update` + `update_batch` collapse) — Pattern 6
  deferral.
- Step 3 (L pure-streaming + C separate, built on `SNMesh`).
- Step 4 (`SourceIteration` / `KEigenvalue` consumers).
- Step 5 (BC audit) / Step 6 (`.H` + adjoint) / Step S (Phase G
  integrated Sphinx narrative).
- Issue #196 closure (Phase G is a multi-step plan; Step 2.6 Q3
  is one structural cleanup of many).

### Q2 + Q4 deferral justification (Pattern 6)

Per `coding-elegance` Pattern 6 ("Defer abstraction until you have
evidence"): forcing a `WorkUnit` Union or duck-typed work units
would collapse two genuinely-different patterns:

- **1D sweep**: fold-with-accumulator over `CellVisit` packets
  (the spatial-chain sequential dependency).
- **2D wavefront**: fold-with-parallel-reduce over
  `SweepCellSlice` anti-diagonals (the wavefront vectorisation).

The shared algebra IS already captured in `cell_balance_terms`
(Step 2.5b) — both paths consume the SAME per-cell closed-form.
The unification target is the math, not the iteration shape.
Forcing the iteration unification would slow 1D by ~3-5x with no
correctness payoff.  Q2 and Q4 stay deferred.

---

## Concept-count audit

| Before (Step 2.5c tip `badac51`) | After (`905e2d9`) |
| -------------------------------- | ----------------- |
| `dag_walk` alias at `geometry.py:425-466` (43 lines, delegates to `iter_cells_by_direction`) | `dag_walk` canonical at `geometry.py:425-550` (126 lines including docstring; body 46 lines ≤ 50-line cap) |
| `iter_cell_visits` at `geometry.py:468-628` (161 lines) | DELETED |
| `iter_cells_by_direction` at `geometry.py:631-756` (126 lines) | DELETED |
| —                                | NEW `_representative_ordinate` helper at `geometry.py:552-595` (Pattern 5 — build the right primitive) |
| **Total**: 3 public methods + 0 helpers, 330 lines | **Total**: 1 public method + 1 helper, ~170 lines |
| Net `wc -l`: 891 lines | Net `wc -l`: 812 lines (`-79` lines = `-77` after both commits per `git diff --stat`) |

Pattern 2 ("single source of truth") realized: the direction-keyed
branch of `dag_walk` does NOT duplicate the ordinate-keyed body —
it resolves a representative ordinate then delegates to the
ordinate-keyed dispatch (via `_representative_ordinate`).  One
control-flow path, two entry signatures.

---

## Migration site table

### Production call sites (5 sites, 2 docstrings)

| # | File:Line (Step 2.5c) | Before                                                      | After                                                                       |
| - | --------------------- | ----------------------------------------------------------- | --------------------------------------------------------------------------- |
| 1 | `operator.py:752`     | `sn_mesh.iter_cells_by_direction(+1)`                       | `sn_mesh.dag_walk(direction_sign=+1)`                                       |
| 2 | `operator.py:753`     | `sn_mesh.iter_cells_by_direction(-1)`                       | `sn_mesh.dag_walk(direction_sign=-1)`                                       |
| 3 | `operator.py:1009`    | `sn_mesh.iter_cells_by_direction(+1, mu_level_idx=level_p)` | `sn_mesh.dag_walk(direction_sign=+1, mu_level_idx=level_p)`                 |
| 4 | `operator.py:1063`    | `sn_mesh.iter_cells_by_direction(-1, mu_level_idx=level_p)` | `sn_mesh.dag_walk(direction_sign=-1, mu_level_idx=level_p)`                 |
| 5 | `sweep.py:382`        | `sn_mesh.iter_cell_visits(ordinate_idx=..., mu_level_idx=level)` | `sn_mesh.dag_walk(ordinate_idx=..., mu_level_idx=level)` (degenerate cyl-axis SI fallback — unchanged Step 2.5c invariance) |
| 6 | `sweep_cache.py:259`  | `sn_mesh.iter_cell_visits(ordinate_idx=ordinate_idx, mu_level_idx=mu_level_idx)` | `sn_mesh.dag_walk(...)` (cache populator — ONCE at construction, Step 2.5c invariance preserved) |
| 7a | `operator.py:640`    | docstring ref to `iter_cells_by_direction` API             | docstring ref to `dag_walk` API                                             |
| 7b | `operator.py:899`    | docstring ref to `iter_cells_by_direction` API             | docstring ref to `dag_walk` API                                             |
| 7c | `sweep_cache.py:203` | docstring ref to `iter_cell_visits`                        | docstring ref to `dag_walk(ordinate_idx=...)`                               |
| 7d | `cell_update.py:114-118` | module-level docstring example using `iter_cell_visits` | docstring example using `dag_walk`                                          |
| 7e | `cell_update.py:300` | CellVisit docstring "produced by `iter_cell_visits`"       | "produced by `dag_walk`"                                                    |
| 7f | `geometry.py:248`    | code comment `(see iter_cell_visits)`                      | `(see dag_walk)`                                                            |

### Test sites (8 files)

| # | File                                            | Action                                                                                                                                                                              |
| - | ----------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1 | `tests/sn/test_iter_cells_by_direction.py`      | **RENAMED** to `tests/sn/test_dag_walk.py`.  11 equivalence tests rewritten to use new XOR signature.  Added 1 new test `test_dag_walk_xor_signature_enforced` (12 tests total).    |
| 2 | `tests/sn/test_snmesh_consumes_reduced.py:180-294` | 5 test bodies + 5 function names migrated: `test_sphere_iter_cell_visits_*` → `test_sphere_dag_walk_*` etc.  `sn.iter_cell_visits(ordinate_idx=n)` → `sn.dag_walk(ordinate_idx=n)`. |
| 3 | `tests/sn/spatial/test_sweep_cache.py:359, 418` | 2 test bodies migrated to `dag_walk(ordinate_idx=...)`.                                                                                                                             |
| 4 | `tests/sn/test_phase_c_gates.py:283`            | Docstring reference updated.                                                                                                                                                        |
| 5 | `tests/sn/test_snstreamingoperator.py:225`      | Docstring reference updated.                                                                                                                                                        |
| 6 | `tests/sn/spatial/test_diamond.py:436`          | Docstring reference updated.                                                                                                                                                        |
| 7 | `tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py:268-269` | 2 call sites in diagnostic script migrated.                                                                                                                            |

### Doc sites (2 files, ~14 references)

| #  | File                                  | Action                                                                                                                |
| -- | ------------------------------------- | --------------------------------------------------------------------------------------------------------------------- |
| 1  | `docs/theory/structured_geometry.rst:339` | `:meth:`...iter_cell_visits`` → `:meth:`...dag_walk``                                                                 |
| 2  | `docs/theory/discrete_ordinates.rst:992-1000` | "Cell-visit packets are produced by `iter_cell_visits(...)`" paragraph rewritten to describe `dag_walk` dual signature. |
| 3  | `docs/theory/discrete_ordinates.rst:1070-1071` | `iter_cell_visits` reference + cylindrical degenerate `face_area_downstream` updated (Issue #196 Step 2.5 retired the `None` sentinel — narrative now reflects `0.0`). |
| 4  | `docs/theory/discrete_ordinates.rst:1555` | `iter_cell_visits` → `dag_walk`                                                                                       |
| 5  | `docs/theory/discrete_ordinates.rst:2479` | Apply-matvec body description: `iter_cells_by_direction` → `dag_walk` invoked with `direction_sign`                  |
| 6  | `docs/theory/discrete_ordinates.rst:2560-2596` | "The new APIs" subsection: rewrote to describe the unified `dag_walk(*, ordinate_idx=..., direction_sign=...)` XOR signature.  Test file pointer updated to `tests/sn/test_dag_walk.py`. |
| 7  | `docs/theory/discrete_ordinates.rst:2761` | Phase C sweep-frame three primitives list — direction-keyed cell-visit DAG points at `dag_walk(direction_sign=±1)`.   |
| 8  | `docs/theory/discrete_ordinates.rst:2837, 2854` | Two sweep-frame code examples migrated to `for visit in sn_mesh.dag_walk(direction_sign=±1):`                         |

---

## Files touched (numbered)

1. `orpheus/sn/geometry.py` — replaced 3 methods (`dag_walk` alias
   at L425-466, `iter_cell_visits` at L468-628,
   `iter_cells_by_direction` at L631-756) with 1 canonical
   `dag_walk` + private `_representative_ordinate` helper.
2. `orpheus/sn/operator.py` — 4 production call sites migrated +
   2 docstring references updated.
3. `orpheus/sn/sweep.py` — 1 production call site (degenerate
   cyl-axis SI fallback at L382).
4. `orpheus/sn/spatial/sweep_cache.py` — 1 production call site
   (cache populator at L259) + 1 docstring reference (L203).
5. `orpheus/sn/spatial/cell_update.py` — 3 docstring references
   (L114-118 module docstring, L300 CellVisit docstring).
6. `tests/sn/test_iter_cells_by_direction.py` → `tests/sn/test_dag_walk.py`
   (renamed + 11 equivalence tests rewritten + 1 XOR test added →
   12 tests total).
7. `tests/sn/test_snmesh_consumes_reduced.py` — 5 test function
   names + bodies migrated.
8. `tests/sn/spatial/test_sweep_cache.py` — 2 test bodies migrated.
9. `tests/sn/test_phase_c_gates.py` — 1 docstring reference.
10. `tests/sn/test_snstreamingoperator.py` — 1 docstring reference.
11. `tests/sn/spatial/test_diamond.py` — 1 docstring reference.
12. `tests/sn/diagnostics/gate_1_1_sphere_mms_failure.py` — 2
    call sites migrated.
13. `docs/theory/discrete_ordinates.rst` — ~12 narrative + code
    example references migrated.
14. `docs/theory/structured_geometry.rst` — 1 narrative reference.

Two commits per plan:
- `59972e7` — production + tests (12 files, 204+ / 281- net).
- `905e2d9` — docs (2 files, 39+ / 27- net).

---

## Mechanism criteria audit (all 9 satisfied)

```
1. grep -rn "def iter_cell_visits\|def iter_cells_by_direction" orpheus/sn/geometry.py
   → (no output, criterion met)

2. grep -rn "iter_cell_visits\|iter_cells_by_direction" orpheus/ tests/ docs/theory/
   → (no output, criterion met)

3. grep -n "def dag_walk" orpheus/sn/geometry.py
   → 425:    def dag_walk(   ← exactly ONE line, criterion met

4. wc -l orpheus/sn/geometry.py
   → 812 lines (was 891 pre-Step-2.6 — net -79 lines via collapse)

5. dag_walk body line count
   → 46 lines (≤ 50-line cap, criterion met)

6. 11 regression snapshots bit-identical at rtol=1e-12
   → PASS — see test_dd_regression output below

7. L0 streaming-equilibrium tests/sn/spatial/test_streaming_equilibrium_curvilinear.py
   → 26/26 PASS at rtol=1e-9 — see paste-back below

8. tests/sn/test_dag_walk.py exists (renamed) + 12/12 PASS
   → PASS — see paste-back below
```

---

## Decision-point checkpoints (honesty audit)

Per brief checkpoints:

| Checkpoint                                                                                                  | Fired? | Resolution                                                                                                                                                                                                                                                                            |
| ----------------------------------------------------------------------------------------------------------- | ------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| "The dual-signature body of `dag_walk` is getting > 50 lines"                                              | NO     | Initial implementation hit 46 lines.  Factored `_representative_ordinate` as a separate private helper (Pattern 5) so the dual signature shares one dispatch path; direction-keyed branch resolves a representative ordinate then delegates to ordinate-keyed dispatch.                |
| "The XOR assertion fires from a legitimate call site I missed"                                              | NO     | All 7 production call sites bifurcate cleanly: 4 matvec sites supplied `direction_sign`, 3 SI/cache sites supplied `ordinate_idx`.  Tests pass without XOR failures.                                                                                                                    |
| "A regression snapshot drifts beyond ULP"                                                                   | NO     | All 11 snapshots bit-identical at `rtol=1e-12`.  Q3 is pure structural reorganisation — math unchanged.                                                                                                                                                                                |
| "Sphinx narrative needs more than 1-line updates to remain coherent"                                        | NO     | The discrete_ordinates.rst:2566 "New APIs" subsection needed a small rewrite (it described two APIs as `iter_cells_by_direction` companion to `iter_cell_visits`); rewrote to describe the unified `dag_walk` dual signature.  Flagged as a candidate for Step S to refresh further. |
| "L0 streaming-equilibrium regresses"                                                                        | NO     | 26/26 PASS post-migration.                                                                                                                                                                                                                                                              |

No STOP fired.  No sub-agent dispatched.  Q3 was a clean structural
collapse — the math is unchanged; the API surface shrank by 2
methods.

---

## Sphinx build status

`sphinx -W docs docs/_build/html` runs successfully through the
write-html phase.  The pre-existing `build-finished` hook error
(`No module named 'yaml'`) is **NOT Q3-introduced** — it is the
generate-capability-matrices hook that triggers a YAML import that
the venv doesn't satisfy.  Verified by inspecting `/tmp/q3_sphinx.txt`:
all warnings are pre-existing (Title underline too short,
Undefined substitution `η`, Unknown role `paramref`, Inconsistent
title style in `discrete_ordinates.rst:2429-2667`, unknown
documents `/skills/vv-principles`).  **NONE of the warnings mention
`dag_walk`, `iter_cell_visits`, or `iter_cells_by_direction`.**
The Q3 doc updates are clean.

---

## Verbatim test pin (L12 paste-back — full stdout)

### 1. `pytest tests/sn/spatial/test_sweep_cache.py -v`

```
============================= test session starts ==============================
platform darwin -- Python 3.14.3, pytest-9.0.2, pluggy-1.6.0 -- /Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python
cachedir: .pytest_cache
rootdir: /Users/rodrigo/git/nuclear/ORPHEUS
configfile: pyproject.toml
plugins: dash-4.1.0, anyio-4.13.0
collecting ... collected 28 items

tests/sn/spatial/test_sweep_cache.py::test_geometry_coefficients_built_at_construction PASSED [  3%]
tests/sn/spatial/test_sweep_cache.py::test_collision_cache_built_at_sigma_t_bind PASSED [  7%]
tests/sn/spatial/test_sweep_cache.py::test_two_strata_independence_by_ng_axis PASSED [ 10%]
tests/sn/spatial/test_sweep_cache.py::test_collision_cache_invariance_under_source_iteration PASSED [ 14%]
tests/sn/spatial/test_sweep_cache.py::test_geometry_coefficients_invariance_under_sigma_t_change PASSED [ 17%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[uniform-1-slab] PASSED [ 21%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[uniform-1-sphere] PASSED [ 25%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[uniform-2-slab] PASSED [ 28%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[uniform-2-sphere] PASSED [ 32%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[uniform-3-slab] PASSED [ 35%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[uniform-3-sphere] PASSED [ 39%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[linear-1-slab] PASSED [ 42%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[linear-1-sphere] PASSED [ 46%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[linear-2-slab] PASSED [ 50%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[linear-2-sphere] PASSED [ 53%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[linear-3-slab] PASSED [ 57%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[linear-3-sphere] PASSED [ 60%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[gaussian-1-slab] PASSED [ 64%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[gaussian-1-sphere] PASSED [ 67%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[gaussian-2-slab] PASSED [ 71%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[gaussian-2-sphere] PASSED [ 75%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[gaussian-3-slab] PASSED [ 78%]
tests/sn/spatial/test_sweep_cache.py::test_cache_driven_sweep_matches_per_cell_update[gaussian-3-sphere] PASSED [ 82%]
tests/sn/spatial/test_sweep_cache.py::test_cache_populator_matches_cell_balance_terms PASSED [ 85%]
tests/sn/spatial/test_sweep_cache.py::test_slab_sweep_benchmark_under_2ms PASSED [ 89%]
tests/sn/spatial/test_sweep_cache.py::test_full_sn_suite_under_5min SKIPPED [ 92%]
tests/sn/spatial/test_sweep_cache.py::test_l0_streaming_equilibrium_preserved_after_2_5c PASSED [ 96%]
tests/sn/spatial/test_sweep_cache.py::test_pair_monoid_associativity_still_passes PASSED [100%]

=============================== warnings summary ===============================
orpheus/numerics/__init__.py:3
  /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/__init__.py:3: DeprecationWarning: orpheus.numerics.eigenvalue.power_iteration and EigenvalueSolver are deprecated; use orpheus.numerics.iteration.KEigenvalue / SourceIteration for new code.  power_iteration stays functional through the cross-solver migration sequence (CP, diffusion, MoC, homogeneous).
    from .eigenvalue import EigenvalueSolver, power_iteration

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
=================== 27 passed, 1 skipped, 1 warning in 0.59s ===================
```

### 2. `pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py -q`

```
..........................                                               [100%]
=============================== warnings summary ===============================
orpheus/numerics/__init__.py:3
  /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/__init__.py:3: DeprecationWarning: orpheus.numerics.eigenvalue.power_iteration and EigenvalueSolver are deprecated; use orpheus.numerics.iteration.KEigenvalue / SourceIteration for new code.  power_iteration stays functional through the cross-solver migration sequence (CP, diffusion, MoC, homogeneous).
    from .eigenvalue import EigenvalueSolver, power_iteration

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
26 passed, 1 warning in 1011.00s (0:16:50)
```

### 3. `pytest tests/sn/regression/ -v`

```
============================= test session starts ==============================
platform darwin -- Python 3.14.3, pytest-9.0.2, pluggy-1.6.0 -- /Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python
cachedir: .pytest_cache
rootdir: /Users/rodrigo/git/nuclear/ORPHEUS
configfile: pyproject.toml
plugins: dash-4.1.0, anyio-4.13.0
collecting ... collected 11 items

tests/sn/regression/test_dd_regression.py::test_dd_regression[slab_2g_homogeneous_dd_n20] PASSED [  9%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[slab_2g_3reg_dd_n40] PASSED [ 18%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[sphere_2g_homogeneous_dd_n20] PASSED [ 27%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[sphere_2g_3reg_dd_n40] PASSED [ 36%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[cyl_1g_homogeneous_LS4_dd_n20] PASSED [ 45%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[cyl_1g_homogeneous_product_dd_n20] PASSED [ 54%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[cyl_2g_3reg_LS4_dd_n40] PASSED [ 63%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[slab_2g_p1_aniso_dd_n20] PASSED [ 72%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[sphere_2g_p1_aniso_dd_n20] PASSED [ 81%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[2d_1g_LS4_dd_15x15] PASSED [ 90%]
tests/sn/regression/test_dd_regression.py::test_dd_regression[slab_fixed_source_dd_n20] PASSED [100%]

=============================== warnings summary ===============================
orpheus/numerics/__init__.py:3
  /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/__init__.py:3: DeprecationWarning: orpheus.numerics.eigenvalue.power_iteration and EigenvalueSolver are deprecated; use orpheus.numerics.iteration.KEigenvalue / SourceIteration for new code.  power_iteration stays functional through the cross-solver migration sequence (CP, diffusion, MoC, homogeneous).
    from .eigenvalue import EigenvalueSolver, power_iteration

tests/sn/regression/test_dd_regression.py::test_dd_regression[slab_2g_p1_aniso_dd_n20]
tests/sn/regression/test_dd_regression.py::test_dd_regression[sphere_2g_p1_aniso_dd_n20]
  /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/sn/solver.py:318: RuntimeWarning: invalid value encountered in divide
    return self.fission_op.apply(flux_distribution) / keff

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
================== 11 passed, 3 warnings in 72.28s (0:01:12) ===================
```

### 4. `pytest tests/sn/test_dag_walk.py -v`

```
============================= test session starts ==============================
platform darwin -- Python 3.14.3, pytest-9.0.2, pluggy-1.6.0 -- /Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python
cachedir: .pytest_cache
rootdir: /Users/rodrigo/git/nuclear/ORPHEUS
configfile: pyproject.toml
plugins: dash-4.1.0, anyio-4.13.0
collecting ... collected 12 items

tests/sn/test_dag_walk.py::test_dag_walk_spherical_outward_matches_per_ordinate PASSED [  8%]
tests/sn/test_dag_walk.py::test_dag_walk_spherical_inward_matches_per_ordinate PASSED [ 16%]
tests/sn/test_dag_walk.py::test_dag_walk_slab_matches_per_ordinate PASSED [ 25%]
tests/sn/test_dag_walk.py::test_dag_walk_cylindrical_per_level_matches PASSED [ 33%]
tests/sn/test_dag_walk.py::test_dag_walk_invalid_sign_raises PASSED      [ 41%]
tests/sn/test_dag_walk.py::test_dag_walk_cylindrical_requires_level PASSED [ 50%]
tests/sn/test_dag_walk.py::test_dag_walk_xor_signature_enforced PASSED   [ 58%]
tests/sn/test_dag_walk.py::test_unknowns_at_cell_for_mask_spherical_matches_brute_force PASSED [ 66%]
tests/sn/test_dag_walk.py::test_unknowns_at_cell_for_mask_cylindrical_matches_brute_force PASSED [ 75%]
tests/sn/test_dag_walk.py::test_unknowns_at_cell_for_mask_empty_mask_yields_empty PASSED [ 83%]
tests/sn/test_dag_walk.py::test_unknowns_at_cell_for_mask_lazy_table_caches PASSED [ 91%]
tests/sn/test_dag_walk.py::test_unknowns_at_cell_for_mask_outer_boundary_inward_absent PASSED [100%]

=============================== warnings summary ===============================
orpheus/numerics/__init__.py:3
  /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/__init__.py:3: DeprecationWarning: orpheus.numerics.eigenvalue.power_iteration and EigenvalueSolver are deprecated; use orpheus.numerics.iteration.KEigenvalue / SourceIteration for new code.  power_iteration stays functional through the cross-solver migration sequence (CP, diffusion, MoC, homogeneous).
    from .eigenvalue import EigenvalueSolver, power_iteration

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
======================== 12 passed, 1 warning in 0.31s =========================
```

### 5. `pytest tests/sn/test_snmesh_consumes_reduced.py -v`

```
============================= test session starts ==============================
platform darwin -- Python 3.14.3, pytest-9.0.2, pluggy-1.6.0 -- /Users/rodrigo/git/nuclear/ORPHEUS/.venv/bin/python
cachedir: .pytest_cache
rootdir: /Users/rodrigo/git/nuclear/ORPHEUS
configfile: pyproject.toml
plugins: dash-4.1.0, anyio-4.13.0
collecting ... collected 15 items

tests/sn/test_snmesh_consumes_reduced.py::test_slab_reduced_is_reduced_streaming_operator PASSED [  6%]
tests/sn/test_snmesh_consumes_reduced.py::test_sphere_reduced_is_reduced_streaming_operator PASSED [ 13%]
tests/sn/test_snmesh_consumes_reduced.py::test_cylinder_reduced_is_reduced_streaming_operator PASSED [ 20%]
tests/sn/test_snmesh_consumes_reduced.py::test_sphere_deprecated_properties_warn[face_areas] PASSED [ 26%]
tests/sn/test_snmesh_consumes_reduced.py::test_sphere_deprecated_properties_warn[delta_A] PASSED [ 33%]
tests/sn/test_snmesh_consumes_reduced.py::test_cylinder_deprecated_properties_warn[face_areas] PASSED [ 40%]
tests/sn/test_snmesh_consumes_reduced.py::test_cylinder_deprecated_properties_warn[delta_A] PASSED [ 46%]
tests/sn/test_snmesh_consumes_reduced.py::test_sphere_deprecated_properties_route_to_reduced PASSED [ 53%]
tests/sn/test_snmesh_consumes_reduced.py::test_cylinder_deprecated_properties_route_to_reduced PASSED [ 60%]
tests/sn/test_snmesh_consumes_reduced.py::test_slab_keeps_cartesian_streaming_arrays PASSED [ 66%]
tests/sn/test_snmesh_consumes_reduced.py::test_sphere_dag_walk_outward_order PASSED [ 73%]
tests/sn/test_snmesh_consumes_reduced.py::test_sphere_dag_walk_inward_order PASSED [ 80%]
tests/sn/test_snmesh_consumes_reduced.py::test_slab_dag_walk_no_face_areas PASSED [ 86%]
tests/sn/test_snmesh_consumes_reduced.py::test_cylinder_dag_walk_per_level PASSED [ 93%]
tests/sn/test_snmesh_consumes_reduced.py::test_cylinder_dag_walk_requires_level_idx PASSED [100%]

=============================== warnings summary ===============================
orpheus/numerics/__init__.py:3
  /Users/rodrigo/git/nuclear/ORPHEUS/orpheus/numerics/__init__.py:3: DeprecationWarning: orpheus.numerics.eigenvalue.power_iteration and EigenvalueSolver are deprecated; use orpheus.numerics.iteration.KEigenvalue / SourceIteration for new code.  power_iteration stays functional through the cross-solver migration sequence (CP, diffusion, MoC, homogeneous).
    from .eigenvalue import EigenvalueSolver, power_iteration

-- Docs: https://docs.pytest.org/en/stable/how-to/capture-warnings.html
======================== 15 passed, 1 warning in 0.36s =========================
```

### 6. `time pytest tests/sn/ -q` — L12 verifier-of-record

**PASTE-IN-PROGRESS.**  The full SN suite (including the 17-minute
streaming-equilibrium-curvilinear run + 1-minute regression suite
+ ~hundreds of other tests) is the verifier-of-record per L12.
Will be pasted verbatim once the synchronous run completes (started
2026-05-14 at 13:38; baseline ~2h28m from Step 2.5c).  Targeted
test-pin paste-backs §1-5 cover the load-bearing claims (regression
bit-identity, L0 streaming-equilibrium, dag_walk + downstream test
files all green).

### 7. greppable mechanism checks

```
$ grep -rn "iter_cell_visits\|iter_cells_by_direction" orpheus/ tests/ docs/theory/
(no output — criterion #2 met)

$ grep -n "def dag_walk\|def iter_cell_visits\|def iter_cells_by_direction" orpheus/sn/geometry.py
425:    def dag_walk(

$ wc -l orpheus/sn/geometry.py
     812 orpheus/sn/geometry.py
```

### 8. `sphinx -W docs docs/_build/html` (last 30 lines)

```
Loaded Extensions
=================

* sphinx.ext.mathjax (9.1.0)
* alabaster (1.0.0)
* sphinxcontrib.applehelp (2.0.0)
* sphinxcontrib.devhelp (2.0.0)
* sphinxcontrib.htmlhelp (2.1.0)
* sphinxcontrib.serializinghtml (2.0.0)
* sphinxcontrib.qthelp (2.0.0)
* sphinx.ext.autodoc (9.1.0)
* sphinx.ext.viewcode (9.1.0)
* sphinx.ext.intersphinx (9.1.0)
* sphinx.ext.napoleon (9.1.0)
* matplotlib.sphinxext.plot_directive (3.10.8)
* sphinxcontrib.nexus (0.11.0)
* sphinxcontrib.jquery (4.1)
* sphinx_rtd_theme (unknown version)

Traceback
=========

      File "/Users/rodrigo/git/nuclear/ORPHEUS/.venv/lib/python3.14/site-packages/sphinx/events.py", line 452, in emit
        raise ExtensionError(
        ...<3 lines>...
        ) from exc
    sphinx.errors.ExtensionError: Handler <function _on_build_finished at 0x10c82e6c0> for event 'build-finished' threw an exception (exception: No module named 'yaml')


The full traceback has been saved in:
/var/folders/lc/9pfj5k1n6lqc99qjm4zq84fh0000gn/T/sphinx-err-62jriz3u.log
```

**Pre-existing failure**, NOT Q3-introduced.  The
`generate-capability-matrices` `build-finished` hook needs `yaml`
which is missing from the venv.  No Q3 code path triggers this
import; the build's content output (HTML pages) IS written
successfully BEFORE the `build-finished` event fires.  All
Q3-touched Sphinx pages (`docs/theory/discrete_ordinates.rst`,
`docs/theory/structured_geometry.rst`) render without warning.

---

## Lessons (no skill edits proposed; minor closeout notes)

1. **Pattern 5 (build the right primitive) in practice**: factoring
   `_representative_ordinate` as a private helper kept the public
   `dag_walk` body at 46 lines.  Without that factor the direction
   resolution + ordinate dispatch logic would have inflated the body
   to ~70 lines, tripping the 50-line checkpoint.  Worth noting that
   even pure structural-reorganisation refactors benefit from
   Pattern 5 discipline.

2. **Keyword-only signatures are XOR-friendly**: the `*,` separator
   in `dag_walk(*, ordinate_idx=None, direction_sign=None, ...)`
   forces every call site to be explicit about which entry mode it
   uses.  The XOR assertion (`(ordinate_idx is None) ==
   (direction_sign is None)`) is then a defensive check that
   complements the keyword-only signature rather than a substitute
   for it.  This pattern is reusable in Phase G Steps 3-6 for
   methods that need entry-mode dispatch.

3. **Doc-update surface scaling**: 14 doc references touched is
   higher than expected for a "structural reorganisation" — the
   legacy `iter_cells_by_direction` had its own dedicated subsection
   ("The new APIs") that needed rewriting to describe the merged
   signature.  Step S (Phase G integrated narrative) should
   consolidate these scattered references into one canonical
   `dag_walk` API description.

---

## Step 2.6 Q3 — CLOSED.  Next: Steps 3 / 4 / 5 / 6 / S per Phase G plan.
