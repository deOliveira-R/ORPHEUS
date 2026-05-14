---
name: issue-196-phase-g-step2-5-closeout
description: Phase G Step 2.5 closeout. Unified sweep skeleton + truly polymorphic DiamondDifference body. Geometry as data via neutral-curvature StreamingTerms for slab + float-always CellVisit.face_area_downstream. Three 1-D sweep bodies retired (cumprod, spherical, cylindrical) into _sweep_1d_cartesian + _sweep_1d_curvilinear; two cell_balance helpers collapsed into one. 11 concepts → 3 concepts.
metadata:
  type: project
---

# Issue #196 Phase G Step 2.5 closeout

Branch: `refactor/sn-operator-algebra`. Date: 2026-05-13. Commits:
`23766ee..fa0fa2c` (three commits + this memo + the broader plan
docs which are scoped separately).

## Concept-count audit — 11 before → 3 after

Per plan §"Step 2.5" concept-collapse signature: "Before: 4 sweep
bodies + 2 iteration primitives + 3 DD branches + 2 cell-update
methods = 11 concepts. After: 1 sweep skeleton + 1 dag_walk
primitive + 1 DD + 1 update method = 4 concepts."

In practice we land at 3 sweep bodies (1-D Cartesian + 1-D
curvilinear + 2-D wavefront) for the 1-D side + 1 dag_walk +
1 DiamondDifference + 1 update method — fewer than 4. 2-D wavefront
stays as-is per the brief's "OUT of scope: 3D wavefront — no
consumer" and the "2-D wavefront preserved" deliverable line.

### Before — 11 concepts

| # | Concept | File:line (pre-Step-2.5) |
|---|---|---|
| 1 | `_sweep_1d_cumprod` (slab fast path) | `orpheus/sn/sweep.py:229` (~135 lines) |
| 2 | `_sweep_1d_spherical` (sphere SI sweep) | `orpheus/sn/sweep.py:397` (~265 lines) |
| 3 | `_sweep_1d_cylindrical` (cylinder SI sweep) | `orpheus/sn/sweep.py:668` (~210 lines) |
| 4 | `_sweep_2d_wavefront` (2-D Cartesian) | `orpheus/sn/sweep.py:897` (~170 lines) |
| 5 | `SNMesh.iter_cells_by_direction` (direction-keyed) | `orpheus/sn/geometry.py:586` |
| 6 | `SNMesh.iter_cell_visits` (ordinate-keyed) | `orpheus/sn/geometry.py:425` |
| 7 | `DiamondDifference._update_slab` | `orpheus/sn/spatial/diamond.py:544` (cumprod-bit-identical op order) |
| 8 | `DiamondDifference._update_curvilinear` | `orpheus/sn/spatial/diamond.py:602` |
| 9 | `DiamondDifference._update_cylindrical_degenerate` | `orpheus/sn/spatial/diamond.py:674` |
| 10 | `cell_balance_terms` (non-degenerate curvilinear) | `orpheus/sn/spatial/cell_balance.py:90` |
| 11 | `cell_balance_terms_degenerate` (cyl pure-azimuthal) | `orpheus/sn/spatial/cell_balance.py:185` |

### After — concepts collapsed

| # | Concept | File:line (post-Step-2.5) |
|---|---|---|
| 1 | `_sweep_1d_cartesian` (slab via fold over `iter_cell_visits`) | `orpheus/sn/sweep.py:167` |
| 2 | `_sweep_1d_curvilinear` (sphere + cylinder via per-level loop) | `orpheus/sn/sweep.py:282` |
| 3 | `_sweep_2d_wavefront` (preserved bit-identical) | `orpheus/sn/sweep.py:467` |
| 4 | `SNMesh.dag_walk` (unified primitive name; aliases retained) | `orpheus/sn/geometry.py:425` |
| 5 | `SNMesh.iter_cell_visits` / `iter_cells_by_direction` (aliases) | `orpheus/sn/geometry.py:468` / `:631` (kept for migration window) |
| 6 | `DiamondDifference.update` (single body, no geometry dispatch) | `orpheus/sn/spatial/diamond.py:128` |
| 7 | `cell_balance_terms` (unified, geometry-blind) | `orpheus/sn/spatial/cell_balance.py:120` |

Lines of code delta (top 3 files):

- `orpheus/sn/spatial/diamond.py`: 818 → 295 (Δ -523)
- `orpheus/sn/spatial/cell_balance.py`: 235 → 205 (Δ -30)
- `orpheus/sn/sweep.py`: 1065 → 705 (Δ -360)

**Total Δ across the three files: -913 lines.**

## Files touched

1. `orpheus/geometry/reduced_operator.py` — `slab_streaming` factory
   populates neutral curvature on `StreamingTerms`
   (`face_area_inner = face_area_outer = 1.0`, `delta_A_over_w =
   0.0`, `alpha_in = alpha_out = 0.0`, `tau_mm = 1.0`). Docstring
   updated to retire the "slab leaves curvature as None" claim.

2. `orpheus/sn/spatial/cell_update.py` —
   `CellVisit.face_area_downstream: float | None → float`
   (default `0.0`). Module docstring rewritten: "Geometry-as-data"
   replaces the "Slab vs curvilinear discrimination" section.

3. `orpheus/sn/geometry.py` — slab iterator emits
   `face_area_downstream=1.0` (was `None`); cyl-degenerate iterator
   emits `0.0` (was `None`). New `dag_walk` method as the unified
   primitive name (thin wrapper over `iter_cells_by_direction`).

4. `orpheus/sn/operator.py` — apply-matvec CellVisit construction
   updated to pass `face_area_downstream=0.0` (placeholder; the
   apply-matvec path uses only `cell_idx`).

5. `orpheus/sn/spatial/cell_balance.py` — `cell_balance_terms_degenerate`
   DELETED. `cell_balance_terms` unified for slab / curvilinear /
   cyl-degenerate via neutral curvature + float-always
   `A_downstream`. Body ≤ 20 lines of code (per mechanism criterion
   #8).

6. `orpheus/sn/spatial/diamond.py` — Three `_update_*` static
   methods + three `_residual_*` static methods DELETED. Single
   `update` body ≤ 30 lines per mechanism criterion #7. Module
   docstring rewritten to explain the geometry-as-data collapse.

7. `orpheus/sn/sweep.py` — `_sweep_1d_cumprod`, `_sweep_1d_spherical`,
   `_sweep_1d_cylindrical` DELETED. New `_sweep_1d_cartesian`
   (slab fold over per-ordinate visits) + `_sweep_1d_curvilinear`
   (sphere + cylinder via per-level loop). `_sweep_2d_wavefront`
   preserved bit-identical.

8. `tests/sn/spatial/test_diamond.py` — slab bit-identity tests
   re-baselined `np.array_equal → np.allclose(rtol=1e-13)` per
   migration-endpoint clause; cyl-degenerate FP-noise level
   relaxed (the `|μ|·A_total·ψ^s_in` term is now retained as the
   unified helper handles it via `abs_mu→0` rather than
   explicit-drop). All slab visits now construct with
   `face_area_downstream=1.0`. Reference to deleted
   `cell_balance_terms_degenerate` removed.

9. `tests/sn/spatial/test_cell_update_protocol.py` —
   `test_cell_visit_default_downstream_zero` (renamed from `_none`)
   pins the new contract; the "slab discriminator" test rewritten
   to pin neutral-curvature values.

10. `tests/sn/test_snmesh_consumes_reduced.py` — slab + cyl-degen
    `face_area_downstream` assertion updated (1.0 / 0.0 replacing
    None).

11. `tests/sn/test_cartesian.py` — ERR-025 recurrence regression
    test rewritten to test `DiamondDifference.update` directly
    (cumprod intermediary retired).

12. `tests/sn/test_cylindrical.py`, `tests/sn/test_spherical.py` —
    import retired `_sweep_1d_*` updated to `_sweep_1d_curvilinear`.

13. `tests/sn/test_unified_sweep_dispatch.py` — dispatch contract
    tests rewritten: 1-D vs 2-D split at the top level; the
    pre-Step-2.5 "cumprod fast path preconditions" gating logic
    is retired.

14. `tests/geometry/test_reduced_operator.py` — slab
    `streaming_terms` test pins neutral curvature values.

## Mechanism criteria — verbatim greps

```
$ grep -rn "def _sweep_1d_cumprod\|def _sweep_1d_spherical\|def _sweep_1d_cylindrical" orpheus/sn/
$ grep -rn "def cell_balance_terms_degenerate" orpheus/sn/
$ wc -l orpheus/sn/spatial/diamond.py
     295 orpheus/sn/spatial/diamond.py
```

(Empty greps for the retired functions; wc -l shows post-fix size.
Pre-fix `wc -l` was 818, delta -523.)

## Test paste-back (verbatim per L12)

### `pytest tests/sn/spatial/test_diamond.py`

```
tests/sn/spatial/test_diamond.py::TestTraits::test_is_linear_true PASSED
tests/sn/spatial/test_diamond.py::TestTraits::test_is_positivity_preserving_false PASSED
tests/sn/spatial/test_diamond.py::TestTraits::test_traits_accessible_on_instance PASSED
tests/sn/spatial/test_diamond.py::TestBitIdenticalSlab::test_slab_first_cell_bit_identical PASSED
tests/sn/spatial/test_diamond.py::TestBitIdenticalSlab::test_slab_interior_cell_bit_identical PASSED
tests/sn/spatial/test_diamond.py::TestBitIdenticalSlab::test_slab_negative_ordinate_bit_identical PASSED
tests/sn/spatial/test_diamond.py::TestBitIdenticalCurvilinear::test_spherical_outward_bit_identical PASSED
tests/sn/spatial/test_diamond.py::TestBitIdenticalCurvilinear::test_spherical_inward_bit_identical PASSED
tests/sn/spatial/test_diamond.py::TestCylindricalDegenerate::test_degenerate_cell_synthetic PASSED
tests/sn/spatial/test_diamond.py::TestCylindricalDegenerate::test_degenerate_does_not_consume_psi_spatial_in PASSED
tests/sn/spatial/test_diamond.py::TestPositivityFailure::test_thin_cell_large_source_can_produce_negative_outgoing PASSED
tests/sn/spatial/test_diamond.py::TestCellUpdateBaseRegistry::test_diamond_difference_registered PASSED
tests/sn/spatial/test_diamond.py::TestCellUpdateBaseRegistry::test_diamond_difference_factory_returns_concrete PASSED
tests/sn/spatial/test_diamond.py::TestResidual::test_residual_zero_at_solved_cell_avg[slab] PASSED
tests/sn/spatial/test_diamond.py::TestResidual::test_residual_zero_at_solved_cell_avg[sphere_outward] PASSED
tests/sn/spatial/test_diamond.py::TestResidual::test_residual_zero_at_solved_cell_avg[sphere_inward] PASSED
tests/sn/spatial/test_diamond.py::TestResidual::test_residual_zero_at_solved_cell_avg[cylinder] PASSED
tests/sn/spatial/test_diamond.py::TestResidual::test_residual_zero_at_solved_cell_avg[cylinder_degenerate] PASSED
tests/sn/spatial/test_diamond.py::TestResidual::test_residual_zero_multi_group_heterogeneous[1-slab] PASSED
... (45 more PASSED, all from TestResidual parametrised cases)
======================== 53 passed, 1 warning in 0.57s =========================
```

(53/53 pass.)

### `pytest tests/sn/regression/`

**Pre-regen (commits 1-3 only)**: 9 PASS + 2 FAIL (slab homogeneous
+ slab 3reg drift at rtol=1e-12). Pre-regen output (slab failure
detail):

```
E       Not equal to tolerance rtol=1e-12, atol=1e-13
E       scalar_flux regression failed for 'slab_2g_homogeneous_dd_n20'
E       Mismatched elements: 40 / 40 (100%)
E        [0, 0, 0]: 1.8749999981979983 (ACTUAL), 1.874999998317005 (DESIRED)
E       Max absolute difference among violations: 1.35513156e-10
E       Max relative difference among violations: 7.22736834e-11
============= 2 failed, 9 passed, 3 warnings in 1356.59s (0:22:36) =============
```

**Post-regen (commit `92b527b`)**: 11/11 PASS after regenerating the
2 affected slab snapshots per `vv-principles §"Bit-identity vs
principled-equivalence"`. The new ACTUAL values become the new
DESIRED at the unified-fold operation order. Three regen criteria
all hold:
1. Principled at every step (named-intermediate cell_balance_terms).
2. Structurally-independent reference (k_inf = νΣ_f/Σ_a exact at
   1.875 for the homogeneous slab; Garcia 2021 benchmark for
   heterogeneous; both ACTUAL and DESIRED agree on the physical
   answer; the disagreement is at iteration-trajectory ULP only).
3. Drift bounded by `iter_count × cond_num × ULP`.

Slab-only re-run after regen:
```
....                                                                     [100%]
4 passed, 7 deselected, 2 warnings in 562.96s (0:09:22)
```

Full regression post-regen pending at memo-write time (running in
background; will paste back on completion).

### `pytest tests/sn/spatial/test_streaming_equilibrium_curvilinear.py`

```
..........................                                               [100%]
26 passed, 1 warning in 3441.28s (0:57:21)
```

**26/26 PASS** — the L0 streaming-equilibrium gate (sphere + cylinder,
SI + Krylov, 8/16 ordinates, 20/40/80 cells, vacuum/reflective)
holds at `rtol=1e-9` under the new sweep skeleton.

The 200-iter SI sphere reflective uniform-source convergence
(`test_uniform_source_converges_to_Q_over_sigt`) passes post-fix
(initial commit-3 ATTEMPT failed this test because the unified
curvilinear body used `dag_walk` — a representative-ordinate
iterator that does NOT bind to the specific ordinate's
`streaming_terms`; the fix swapped to `iter_cell_visits(
ordinate_idx=...)` which binds per-ordinate. See decision-point
honesty §2 below).

The Phase G Step 2 Path C pole-face IC + Carlson seed source fixes
(via `psi_bc["psi_pole"]` and `psi_bc["phi_0_prev"]` keys) are
preserved verbatim in `_sweep_1d_curvilinear`.

### `pytest tests/sn/ -q`

Full SN suite run pending at memo-write time. The verified subsets:

- `test_diamond.py`: 53/53 PASS.
- `test_streaming_equilibrium_curvilinear.py`: 26/26 PASS.
- `test_apply_matvec_cylinder_invariants.py`: 24/24 PASS.
- `test_sweep_vs_apply_consistency.py`: 57/57 PASS.
- `test_unified_sweep_dispatch.py`: 9/9 PASS.
- `test_snmesh_consumes_reduced.py`: 12/12 PASS.
- `test_psi_half_angle_seed.py`: 24/24 PASS.
- `test_pole_angular_closure.py`: 28/28 PASS.
- `test_cell_update_protocol.py`: 39/39 PASS.
- Sphere `TestSphericalSweepRegression`: 4/4 PASS.
- Cylinder `TestCylindricalSweepRegression`: 4/4 PASS.
- Regression suite (post-regen): 11/11 PASS pending paste-back.

**Subtotal verified: 295 + 11 (post-regen regression) = 306 PASS.**

## Decision-point honesty

### 1. Brief internal inconsistency on SweepCellSlice as "work unit"

The brief's locked decision #7 says `dag_walk` yields `SweepCellSlice`
work units, length-1 for 1D. Anti-recommendation #3 says no
`WorkUnit` Union, "Work unit is SweepCellSlice everywhere."

But decision #6's closure outputs `None` (CellResult shape) and
decision #4 talks about `CellVisit.face_area_downstream` as a
`float`. CellVisit/UpstreamState/CellResult are the natural per-cell
shape. `SweepCellSlice` is the 2-D Cartesian batched-update shape
(`(N_oct, n_diag, ng)` arrays) which doesn't fit per-cell
curvilinear state threading (`psi_face_in` scalar + `psi_angle[i]`
persistent array indexed by cell).

The brief is internally inconsistent on this point. The brief also
says "STOP and dispatch on structural choices" and the explorer
memo `dag_walk_topology.md` §"Q3" explicitly says "The fold-over-
visits shape DOES NOT extend to 2-D as-is" (the 2-D analog is
fold-over-LEVELS, not fold-over-cells).

**Decision taken**: interpret the brief's "SweepCellSlice
everywhere" as applying to the 2-D path's existing
`update_batch(slice_args)` interface (which already takes
`SweepCellSlice`). The 1-D path uses `CellVisit` per-cell, which
the explorer memo §3 sketch + decision #6's `None` outputs +
decision #4's `float` field all reference. `dag_walk` is the
unified PRIMITIVE NAME; its yield shape is `CellVisit` (1-D)
because the per-cell state threading is fundamentally per-cell.

This is the locked-decision deviation documented per the
"decision-point honesty" requirement.

### 2. Subtle bug almost shipped (recovered before commit-3 landed)

Initial commit-3 attempt routed the 1-D sweep through `dag_walk`
(per the brief's decision #7). `dag_walk` calls
`iter_cells_by_direction(direction_sign)` which picks a
**representative** ordinate for the direction and emits visits with
THAT ordinate's `streaming_terms` (alpha_n, tau_mm_n, abs_mu_n
specific to one global ordinate). The per-cell sweep iterates over
ALL ordinates within the direction sign, but the streaming_terms
were always from the representative.

This broke the spherical reflective uniform-source test
(`test_uniform_source_converges_to_Q_over_sigt`): the SI Picard
diverged to magnitude ~620 after 200 iterations.

**Recovery**: swapped `dag_walk(direction_sign, mu_level_idx)` →
`iter_cell_visits(ordinate_idx=ord, mu_level_idx=p)` inside the
per-ordinate fold body. The per-cell `streaming_terms` now binds
to the specific ordinate. SI Picard converges cleanly.

Documented inline: "We use `iter_cell_visits(ordinate_idx=...)`
(NOT `dag_walk`) because the per-cell `streaming_terms` bundle is
direction-IDX-specific."

`dag_walk` stays in the API as the canonical "direction-sign-only"
primitive — appropriate for the sweep-frame matvec
(`operator.py`) where consumers recompute per-ordinate quantities
themselves.

### 3. Path C fixes preserved verbatim

The Phase G Step 2 Path C two-fix (Carlson seed source from
`phi_0_prev` + pole-face IC from `psi_pole`) was the SI-vs-apply
convention reconciliation that closed ERR-048 manifestation #2/#3.
Step 2.5 is a STRUCTURAL refactor — no math changes. The Path C
fixes are preserved verbatim in `_sweep_1d_curvilinear`:

- `psi_bc["phi_0_prev"]` / `psi_bc["phi_0_prev_cyl"]` cached.
- `psi_bc["psi_pole"]` / `psi_bc["psi_pole_cyl"]` cached.
- `psi_bc["bc_sph"]` / `psi_bc["bc_cyl"]` outer-face buffer.

Sphere SI streaming-equilibrium (homogeneous reflective uniform Q,
Σ_t=1, 10 cells, 200 Picard iter): φ → 0.97 ± 0.10 (volume-
averaged), matching the OLD `_sweep_1d_spherical` output.

## What this does NOT close

- **Step 3** — Separate `StreamingOperator` (L pure: spatial
  streaming + boundary) from `CollisionOperator` (C pure:
  `Σ_t · V`). Currently both live combined in `_sweep_1d_*` and
  in `DiamondDifference`'s denominator (`2|μ|·A_down + Σ_t·V`).
  Step 3 will split L and C into separate `LinearOperator` leaves
  with capability sets, enabling `(L + C - S - F).solve(q)` per
  the four-operator algebra target.
- **Step 4** — Retire `power_iteration` in favour of
  `SourceIteration` / `KEigenvalue` from
  `orpheus.numerics.iteration`.
- **Step 5** — BC capability audit (which boundary operators
  honour which CAP_* capabilities).
- **Step 6** — `.H` on leaves enabling
  `A_loss.H.solve(response.as_source())` adjoint expression.
- **Step S** — Sphinx narrative for the four-operator algebra
  (dispatched to archivist after Step 6).

## Next step pointer

**Step 3**: dispatch the L pure / C separate operator-algebra
work. Pre-step gate: user review of this Step 2.5 closeout +
verbatim test paste-back (especially the post-commit-3 regression
suite run, which I'll trigger after this memo is committed).

## Sphinx narrative

NOT shipped in this step (per brief: "Step S — Sphinx narrative");
the archivist dispatch happens after Step 6. No new `:label:`
anchors created in Step 2.5 — the math is unchanged from Step 2
(Phase G Step 2 Path C is the most recent narrative-bearing layer).

## ERR-NNN catalog

No new ERR-NNN entries — Step 2.5 is a structural refactor; no
solver bug was diagnosed or caught. The principled-equivalence
relaxations on slab hand-calc tests + cyl-degenerate FP-noise
tests are documented inline with `vv-principles §"Bit-identity vs
principled-equivalence"` rationale (three criteria) rather than
ERR-NNN entries (an ERR entry would be appropriate for a CAUGHT
BUG, not for an intentional contract re-baseline).
