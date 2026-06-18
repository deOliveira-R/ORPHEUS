# #236 Phase 2 — τ carve STEP C: geometry-τ retirement (migrate-then-delete)

**Branch:** `feature/sn-spatial-angular-product` · **HEAD at start:** `25e564b`
· **Checkout:** MAIN `/Users/rodrigo/git/nuclear/ORPHEUS` (not a worktree)
· **Date:** 2026-06-18 · **NOT committed** (left for main-agent review).

The terminal step of the τ carve. After B3 made the LIVE sweep/scan/matvec
read closure-owned τ/c off `CellVisit`, geometry-τ (`StreamingTerms.tau_mm` /
`alpha_in` / `alpha_out`, `ReducedStreamingOperator.tau_mm` /
`tau_mm_per_level`, and the geometry-side τ producer that baked them) had ZERO
LIVE production readers but was STILL load-bearing as a TEST ORACLE. So Step C
is **migrate-then-delete**: re-point the test oracles onto the independent
`contamination.morel_montry_weights` reference FIRST (green with geometry-τ
present), then surgically excise the dead producer + fields.

Plan of record: `.claude/plans/issue_236_stepC_dependency_audit.md` (the L20
3-surface audit) + the §2 producer/consumer/reference map in
`.claude/plans/issue_236_phase2_tau_carve_crosswalk.md`.

---

## PHASE 1 — test-oracle migration (test-only; bit-green with geometry-τ present)

### The independent oracle (vv L11)
`contamination.morel_montry_weights(quad, geometry)`
(`orpheus/derivations/discrete/sn/contamination.py:156`) — a STRUCTURALLY
INDEPENDENT code path to the same BMC-2010-Eq.43 weight, NOT the closure
producer `morel_montry_tau_per_level`. **UNCLAMPED** on both geometries;
production cylinder τ is `clip(τ_raw, 0.5, 1.0)`, sphere is unclamped (mirrors
the streaming factories). Verified 0-ULP `morel_montry_weights == op.tau_mm`
(sphere) and `clip(morel_montry_weights) == op.tau_mm_per_level` (cyl) BEFORE
any deletion (with the clamp biting on the most-inward cylinder ordinate).

### (1a) The c-surrogate — `tests/sn/sweep/core/_c_surrogate.py`
Signature changed (it can no longer take a bare `st`; `st.tau_mm`/`alpha_*`
vanish after Step C). Now TWO entry points (hand-transcribed, do NOT import
the production closure — L11):
- `c_from_constants(tau, alpha_in, alpha_out) -> (c_in, c_out)` — the pure
  formula `c_out = α_out/τ`, `c_in = (1-τ)/τ·α_out + α_in`.
- `mm_constants_for_ordinate(op, cell_idx, direction_idx, mu_level_idx=None)
  -> (tau, alpha_in, alpha_out)` — resolves the triple for one ordinate from
  the geometry operator + the INDEPENDENT τ oracle (clamped for cyl) + the
  SURVIVING dome (`op.alpha_half` sphere / `op.alpha_per_level` cyl). Slab =
  neutral `(1.0, 0.0, 0.0)`.
The old `c_from_streaming_terms(st)` is GONE (zero usages tree-wide).

### (1b) The τ-stamp catcher — `tests/sn/sweep/core/test_cell_visit_c_stamp.py`
`_assert_visit_stamp_matches_inline` re-pointed to the independent source.
Still walks REAL `dag_walk` visits (sphere + multi-level cyl + slab) and
asserts per-visit `visit.c_in`/`visit.c_out`/`visit.tau` == the independent
recompute at 0-ULP. Threads `(op, direction_idx, mu_level_idx)` from the walk
loop (the visit alone doesn't carry the operator). The whole POINT survives
(pin the stamp's ordinate-map AND the closure producer against an independent
reference) — MUTATION-VERIFIED post-deletion (see below).

### (1c) The Leg-1 gate — `tests/sn/sweep/curvilinear/test_tau_producer_equivalence.py`
RETIRED the 2 `*_equals_geometry_factory_0ulp` legs (×2 params each = 4 tests:
they compared closure-τ vs `spherical_streaming(...).tau_mm` /
`cylindrical_streaming(...).tau_mm_per_level`, which Step C deletes). KEPT the
`*_matches_independent_reference` (sphere) + `*_clamp_is_the_only_difference`
(cyl negative control) + `*_identity_closure_tau_is_neutral_one` legs — the
surviving structurally-independent producer-equivalence floor (closure-τ vs
contamination.py). Dropped now-unused imports (`spherical_streaming`,
`cylindrical_streaming`, `MorelMontryAngularSweep`) + the `_sphere_mesh` /
`_cyl_mesh` helpers. 12P → 5P.

### (1d) The inline-recompute constructors (test fixtures)
- `tests/sn/sweep/core/test_diamond.py` — 4 geometry-factory tests +
  4 visit-factory helpers (`_slab/_sphere/_cylinder/_cylinder_degenerate_visit_inputs`)
  + 2 synthetic-degenerate tests + 1 inward-residual test: all resolve the
  triple via `mm_constants_for_ordinate(op, ...)` then `c_from_constants`. The
  3 synthetic `StreamingTerms(...)` constructors dropped `alpha_in=`/
  `alpha_out=`/`tau_mm=` kwargs.
- `tests/sn/sweep/core/test_cell_balance_for_streaming.py` — synthetic packet
  factories dropped the τ/α kwargs; τ/α now module constants
  (`_CURVILINEAR_*` / `_SLAB_*`); `_scalar_to_vectorized_inputs` takes the
  triple as keyword args; Test-5 parametrize carries the triple.
- `tests/sn/sweep/core/test_ordinate_scan.py` — the 3 visit-builders source
  the triple via `mm_constants_for_ordinate`; `_affine_coefficients_from_visits`
  now READS `visit.c_in`/`visit.c_out` off the stamped visit (the c can no
  longer be recomputed from `visit.streaming_terms`).

### (1d-extra) — readers the audit's explicit list did NOT name (found by L12 grep)
- `tests/sn/sweep/core/test_discretization_scheme_protocol.py` —
  `FakeCurvilinearStrategy` read `st.tau_mm`/`st.alpha_*`; re-pointed to read
  `visit.tau`, shape-check on the surviving `delta_A_over_w`; the 2
  neutral-curvature tests dropped the deleted-field assertions; the
  strategy-driven test stamps `visit.tau` from the surrogate.
- `tests/sn/sweep/curvilinear/test_w1_clamp_silent_on_flat.py` (a
  `@catches("ERR-026")` gate) — `_flat_field_coeff`/`_redist_net` now take
  `alpha_in`/`alpha_out` explicitly (from `op.alpha_half`); the W1
  τ-independence proof is intact.
- `tests/sn/primitives/test_snmesh_consumes_reduced.py` — dropped the
  `sn.reduced.tau_mm`/`tau_mm_per_level is not None` assertions (kept the
  surviving alpha/redist).
- `tests/geometry/test_reduced_operator.py` — RETIRED `test_tau_mm_bit_identical`
  + `test_tau_mm_per_level_bit_identical` (subsumed by the Leg-1 closure floor);
  surgically dropped the τ assertions from `test_alpha_dome_bit_identical_across_N`,
  `test_full_per_level_hash_equality`, `test_slab_curvature_arrays_are_none`,
  and the 3 streaming_terms-extraction tests.
- `tests/sn/sweep/curvilinear/test_unified_matvec_cylinder.py` — the hand-ref
  matvec drives the T5 legacy `__call__`-arg `tau_mm` (which SURVIVES Step C),
  but sourced its τ from `reduced.tau_mm_per_level` (deleted) → re-sourced from
  the surviving closure producer `morel_montry_tau_per_level(quad, CYLINDRICAL)`.

**PHASE 1 verify:** each migrated file ran GREEN with geometry-τ STILL present
(catcher 3P, Leg-1 5P, diamond 53P, cell_balance 9P, ordinate_scan 52P, scheme
17P, w1_clamp 4P) — proving the migration is faithful before any deletion.

---

## PHASE 2 — surgical deletion (production)

All SURGICAL (per the audit's hazards — the τ statements are interleaved with
live `alpha_half`/`redist_dAw`/`mu_start`/face outputs; whole-function deletion
would be WRONG):

`orpheus/geometry/reduced_operator.py`:
- `spherical_streaming` — deleted the `mu_edge` array + `tau_mm` loop
  (`:681-688`); replaced `mu_start=float(mu_edge[0])` with `mu_start=-1.0`
  (the Hébert starting direction, value-identical). KEPT `alpha`/`redist_dAw`/
  `face_areas`/`delta_A`.
- `cylindrical_streaming` — deleted the `eta_edge`/`tau` blocks + the
  `tau_mm_per_level` field (`:792-815`); the per-level loop now computes ONLY
  `sin_theta` → `mu_start_per_level.append(-sin_theta)` (KEPT). `eta`/`M`
  locals removed (they fed only the τ block).
- `streaming_terms()` accessor — slab arm dropped `alpha_in=0.0`/`alpha_out=0.0`/
  `tau_mm=1.0`; sphere arm dropped the `assert self.tau_mm is not None` +
  `alpha_in=`/`alpha_out=`/`tau_mm=` bakes; cyl arm dropped the
  `assert self.tau_mm_per_level is not None` + `tau_lv`/`alpha_lv` locals +
  the bakes.
- Dataclass fields DROPPED: `StreamingTerms.tau_mm` / `.alpha_in` / `.alpha_out`;
  `ReducedStreamingOperator.tau_mm` / `.tau_mm_per_level`. Docstrings updated
  (class shape contract, slab arm, factory docstrings).

**Verified post-deletion:** all 5 fields GONE (`dataclasses.fields`), all
3 production constructors (sphere/cyl/slab) build cleanly, geometry survives
(`op.alpha_half`/`alpha_per_level`/`mu_start`/`mu_start_per_level` intact).

**LEFT UNTOUCHED (T5 — the legacy `__call__`-arg `tau_mm`):**
`pole_angular_closure.py` Protocol `:320`, ABC abstractmethod `:577`, the
spherical/cyl `__call__` bodies `:1468/1512-1552`, Identity ignore `:1777/1783`
— the unbound `MorelMontryAngularSweep(sn_mesh=None)` runtime-arg path. SEPARATE
surface; survives Step C. Confirmed zero edits to `pole_angular_closure.py`.

**NO production READER of the deleted fields remains** (L12 grep): the only
`self.alpha_in/out` reads in `orpheus/` are the UNRELATED trajectory-resolvent
`chord_oracle.py` Green's-function API (different namespace, audit T4); the
only `.tau_mm` references in `orpheus/sn/` are comments.

---

## VERIFICATION GATES (verbatim summary lines, L12)

### Baseline (HEAD `25e564b`, before any change)
```
Gate 1 (sweep/core+spatial+cartesian_2d+curvilinear -W DriftWarning):
  762 passed, 2 skipped, 31 xfailed in 868.95s
Gate 2 (solve -W DriftWarning, #212 deselect): 60 passed in 80.21s
Gate 3 (operators + reduced_operator): 7 failed, 554 passed, 4 skipped, 1 xfailed
  (the 7 are PRE-EXISTING #195/#209/#214 — listed below)
```

### Phase 1 (migration, geometry-τ STILL present) — per-file green
catcher 3P · Leg-1 5P · diamond 53P · cell_balance 9P · ordinate_scan 52P ·
scheme 17P · w1_clamp 4P · solve 60P (no drift) · operators+reduced
7F/552P (same pre-existing 7).

### Phase 2 (post-deletion) — verbatim
```
Gate 1 (sweep/core+spatial+cartesian_2d+curvilinear -W DriftWarning):
  758 passed, 2 skipped, 31 xfailed in 872.22s   (run twice, identical)
Gate 2 (solve -W DriftWarning, #212 deselect): 60 passed in 79.13s
Gate 3 (operators + reduced_operator):
  7 failed, 552 passed, 4 skipped, 1 xfailed in 5.39s
```
- **NO DriftWarning escalation** → 0 failures on the bit-id gates = the
  deletion changed NO live numeric (it removed DEAD code + migrated oracles).
- **Count delta reconciles EXACTLY:** Gate 1 762→758 = the 4 retired Leg-1
  geometry-factory legs (the only retired tests in those dirs). Gate 3
  554→552 = the 2 retired `test_reduced_operator.py` τ-bit-identical tests.

### The 7 pre-existing operators reds (UNCHANGED count + identity, baseline vs Phase 2)
1-3. `test_bc_extraction_matvec.py::TestVacuumMatvecBitIdentity::test_vacuum_bulk_bit_identical_1d[{0,1,2}-SPH]`
4. `test_boundary_conditions.py::TestSNBCResolution::test_2d_mesh_resolution` (Face 'ymin' mu_y)
5. `test_native_matvec.py::TestTwoDCartesianRaises::test_two_d_cartesian_loss_action_returns_result` (Face 'ymin' mu_y)
6-7. `test_streaming_operator.py::TestT4cPreT4RegressionSnapshotCurvilinear::test_sphere_{1g,2g}_apply_bit_identical`
(the #195/#209 curvilinear-matvec snapshots + #214 2-D; NONE touch the τ carve.)

### MUTATION-VERIFICATION of the migrated catcher (post-Phase-2, against the independent oracle)
- τ stamp `× 1.1` in `_make_cell_visit` → `test_cell_visit_c_stamp.py`
  REDS sphere + cyl + slab (slab too: `1.0·1.1 ≠ 1.0`). REVERTED by re-edit → 3P.
- `c_in ↔ c_out` swap → REDS sphere + cyl (slab passes: c=0 neutral, swap
  invisible — correct). REVERTED by re-edit → 3P + Leg-1 5P green.
⇒ the migrated catcher is a REAL catcher reading the independent reference,
NOT the closure (L11 confirmed).

### L11-independence confirmation
The migrated oracles (`_c_surrogate.mm_constants_for_ordinate`, the catcher,
the Leg-1 surviving legs) read `contamination.morel_montry_weights` — a
DIFFERENT code path from the closure's `morel_montry_tau_per_level` (the
system-under-test). No production-closure import in any oracle.

### pyright delta (touched production files)
`reduced_operator.py` 8 errors / `geometry.py` 22 errors — ALL pre-existing
#226 wrong-rooting noise (`Mesh1D.widths/volumes/areas`, `AngularMeasure.eta/
level_indices`, `Mesh2D.areas`, `float|None`/`int|None` arg types). ZERO
reference the deleted fields (`tau_mm`/`alpha_in`/`alpha_out`). Deletions
REMOVE code → cannot introduce new errors; ZERO new.

---

## DELIVERABLES
- **Sphinx: EXPLICITLY OUT OF SCOPE** per the Step-C brief ("Do NOT edit
  Sphinx (`docs/`) — the main agent handles the doc-only Step-6 via the
  archivist separately"). No stub written this step. The carve narrative
  already lives in the B3 stub `sn-tau-c-on-cellvisit-live`; the geometry-τ
  retirement's `docs/` τ-reference cleanup
  (`discrete_ordinates.rst`/`structured_geometry.rst`) is the
  main-agent→archivist Step-6 deliverable.
- archivist DISPATCH_REQUEST: emitted (followup:false) for the `docs/` τ
  cleanup.
- NOT committed (left for main-agent review).

## OWED follow-ups
- #248 (orphaned `PoleAngularClosure` Protocol after the B2 ABC retype).
- Sphinx `docs/` doc-only update (`discrete_ordinates.rst` /
  `structured_geometry.rst` τ references) — main-agent dispatches the archivist
  (Step-6, out of this agent's scope per the brief).

## LESSON
A migrate-then-delete of a field that is DEAD in production but LIVE as a test
ORACLE has a wider blast radius than the audit's explicit reader list — an L12
grep of the WHOLE test tree (not just the audit's named files) is mandatory
before the field deletion, because oracle readers hide in `@catches` gates
(`test_w1_clamp_silent_on_flat`), fake strategies
(`test_discretization_scheme_protocol`), per-level hash anchors
(`test_reduced_operator`), and the T5-legacy-arg hand-references
(`test_unified_matvec_cylinder`) that no carve plan enumerates. The surrogate
that the catchers + fixtures share must be re-pointed at the
STRUCTURALLY-INDEPENDENT producer (contamination.py, NOT the closure) with the
production clamp replicated (cyl `clip(τ_raw, ½, 1)`, sphere unclamped) — and
mutation-verified RED after the deletion, because a faithful-looking migration
can quietly turn the catcher tautological. The deletion itself is bit-id-proven
by the DriftWarning-escalated gates showing 0 failures + a count delta that
reconciles EXACTLY to the retired tests (no silent test loss).
