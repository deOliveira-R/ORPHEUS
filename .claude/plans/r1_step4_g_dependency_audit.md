# R-1 Step 4 Step G — dependency audit

**Branch**: `refactor/sn-operator-algebra` @ `a638798`
**Date**: 2026-05-22
**Trigger**: proactive `explorer` dispatch per
`.claude/lessons.md` L20 ("Retirement requires an upstream dependency
audit"); operationalised in the session-2 plan §2.1 + §3-§4
(`.claude/plans/r1_step4_session2_followup.md`).

---

## CORRECTION — user direction 2026-05-22

The audit's SURPRISE-1 originally concluded `EquationMap` is a KEEPER.
**That conclusion is WRONG and is hereby superseded.**  Per user
direction:

> "Equation map will likely not stay.  It is helper for legacy
> architecture.  The legacy architecture needs to be retired and
> equation map is going with it.  Differentiate what is legacy
> architecture from the path forward.  The path forward is operator
> algebra, unified indices that do not require adapters."

### Legacy architecture — RETIRES (in Step G + Phase A)

The **packed-1D vector layout with `EquationMap` slot-dispatch** and
everything that bridges into it:

| Symbol | Role in legacy |
|---|---|
| `EquationMap` (class) | Flat-vector index map `(ord, cell, group, face) → packed slot`. |
| `build_equation_map`, `build_equation_map_spherical`, `build_equation_map_cylindrical` | Geometry-dispatched eq_map factories. |
| `build_equation_map_with_traces` | Face-aware variant — **still legacy** (face slots indexed by `face_outer_ordinate` / `face_inner_ordinate` int arrays, an adapter for the typed `BoundaryFlux`'s native `(N, ng)` shape). |
| `solution_to_angular_flux`, `*_spherical`, `*_cylindrical` | Packed → typed decoders. |
| `solution_to_angular_flux_with_traces`, `pack_with_traces` | Face-aware adapter pair — **still legacy**. |
| `SNStreamingOperator` | The bundle that runs on packed input. |
| `transport_operator_matvec` (2-D FD), `transport_operator_matvec_unified` | The cell body of the unified matvec takes native `(N, ng, nx, ny)`, but the **face state** still arrives as packed sub-arrays `(n_face_outer, ng)` + a slot-ordinate map `face_outer_ordinate`.  The face-state interface is legacy; the cell-body interface is already path-forward. |
| `AngularFlux.to_flat_with_traces`, `AngularFlux.from_flat_with_traces` | The typed ↔ legacy bridges that exist solely to feed packed-vector consumers. |
| `StreamingOperator._ensure_eq_map`, `CollisionOperator._ensure_eq_map`, `_eq_map` cache, `n_unknowns` | eq_map machinery inside the typed operators — a sign the typed operators currently wrap the legacy compute via slot lookup (`eq_map.face_outer_ordinate`) instead of consuming `BoundaryFlux` natively. |
| `_make_sweep_preconditioner`, `_build_rhs_cartesian`, `_build_rhs_spherical`, `_build_rhs_cylindrical` | `solver.py` packed-vector helpers. |

### Path forward — operator algebra with unified native indices

The typed `AngularFlux` + `BoundaryFlux` carry **everything natively**:

| Field | Shape | Meaning |
|---|---|---|
| `AngularFlux.values` | `(N, ng, nx, ny)` | Cell-centre ψ indexed natively by `(ord, group, cell-x, cell-y)`. |
| `BoundaryFlux.xmax_face` | `(N, ng)` | Outer-face ψ indexed natively by `(ord, group)`.  NO `eq_map.face_outer_ordinate` slot lookup — the boundary face IS indexed by ordinate. |
| `BoundaryFlux.xmin_face` | `(N, ng)` | Inner-face ψ (slab only) — same shape. |

The path-forward compute primitives:

* **`transport_sweep(source, sig_t, sn_mesh, boundary, *, initial_guess=None) → AngularFlux`** — already path-forward (consumes typed `PerOrdinateSource` + `BoundaryFlux`, emits typed `AngularFlux`).
* **NEW NEEDED — native-shape unified matvec** — consumes `AngularFlux` + `BoundaryFlux` natively; emits `AngularFlux` (with non-zero face residual at outer + slab-inner).  Derives "which ordinates are inflow at each face" from the quadrature (`quad.mu < 0` etc.), not from `eq_map.face_outer_ordinate`.  Replaces `transport_operator_matvec_unified`'s packed-face-slot interface.
* All operators (`L`, `C`, `S`, `F`) consume `AngularFlux` and emit `AngularFlux` — no `_ensure_eq_map`, no `_with_traces` round-trip, no eq_map cache.

### Updated retirement sequence (canonical — supersedes the original 8-step list below)

1. **G0 (NEW pre-work)** — write a native-shape unified matvec that consumes `AngularFlux.values` + `BoundaryFlux.xmax_face`/`xmin_face` directly, derives inflow ordinates from `sn_mesh.quad`, and returns cell + face residuals in matching native shapes.  Pin equivalence (bit-identity or principled-equivalent per `vv-principles`) to the current `transport_operator_matvec_unified` packed-face-slot path on slab + sphere + cylinder.
2. **G1** — carve `_solve_fixed_source_krylov` onto `KrylovAcceleration` + typed `AngularFlux` (1-D only; `NotImplementedError` on 2-D, mirroring `_solve_krylov`).
3. **G2** — carve `_solve_fixed_source_si` onto `SourceIteration` + typed `AngularFlux` (geometry-agnostic via `transport_sweep` — handles 2-D Cartesian naturally; SURPRISE-5).
4. **G3a** — delete `_make_sweep_preconditioner` + `_build_rhs_{cartesian, spherical, cylindrical}` (mechanical; closes #174).
5. **G3b** — migrate `StreamingOperator._apply_typed` + `CollisionOperator._apply_typed` to the native matvec from G0; drop `_ensure_eq_map`, `_eq_map`, `n_unknowns`.
6. **G3c** — retire the legacy `transport_operator_matvec_unified` (packed-face-slot signature).
7. **G3d** — retire `pack_with_traces`, `solution_to_angular_flux_with_traces`, `build_equation_map_with_traces`, `EquationMap`.
8. **G3e** — retire `solution_to_angular_flux_{spherical, cylindrical}` + `build_equation_map_{spherical, cylindrical}`.
9. **G3f** — retire `SNStreamingOperator` class (closes #199 items 1-4; migrates `test_phase_c_gates.py` per SURPRISE-7).
10. **G3g** — retire `AngularFlux.to_flat_with_traces` + `AngularFlux.from_flat_with_traces` (no typed ↔ legacy bridge left to feed).
11. **G3h — DEFERRED to Phase A** — retire Cartesian `build_equation_map` + Cartesian-1-D `solution_to_angular_flux` (blocked by 2-D Cartesian absorption; #199 item 5).
12. **Phase 6 / H** — Sphinx narrative (~3-4 h minimum per SURPRISE-8).

### What Phase 5 becomes

The original Phase 5 (G4) plan was "relocate `_with_traces` family +
`EquationMap` to `angular_flux.py`".  Under the corrected scope,
**every Phase 5 target retires instead** (G3d / G3g above).  Phase 5
of the session-2 plan is **empty** — either delete it or note "all
former relocation targets retire in G3d / G3g".

### Operational consequence for the audit's per-target rows below

The §"Per-target audit rows" section that follows is **STILL ACCURATE
for callgraph data** — every "production caller" / "test caller" /
"diagnostic caller" line was verified at `a638798` and remains the
authoritative reference for HOW the symbol is consumed today.  What
changes is the **retirement-target classification**: read every
"Internal-only after retirement? YES (keeper)" note below as
"retires alongside the legacy architecture — NOT relocated".  The
"Retirement order" notes in each row are superseded by the 12-step
sequence above.

---

**Why this audit exists** — session 1's original Step G plan said
"retire X" without enumerating "who calls X". Mid-session it
under-estimated scope by ~2× and forced a full re-plan. Every
retirement target below has a documented callgraph row BEFORE any
retirement code is touched, so Step G is bug-by-construction
contained.

**Scope clarification (from the brief)**:

* The 12 NON-`_with_traces` targets RETIRE in Step G.
* The 3 `_with_traces` variants (`build_equation_map_with_traces`,
  `solution_to_angular_flux_with_traces`, `pack_with_traces`) are
  KEEPERS — they back `AngularFlux.from_flat_with_traces` /
  `AngularFlux.to_flat_with_traces`. Phase 5 relocates them
  to `orpheus/sn/angular_flux.py`; for this audit they are
  out-of-scope.
* Convention legend: a "production caller" is any reference under
  `orpheus/`; a "test caller" is under `tests/`; a "diagnostic
  caller" is under `derivations/` or `scratch/`; a "doc caller" is
  under `docs/`. Each row separates these four categories.

---

## Per-target audit rows

### (1) `orpheus.sn.operator.SNStreamingOperator` (class)

- **Defined at**: `orpheus/sn/operator.py:1316` (`@dataclass class
  SNStreamingOperator(LinearOperatorMixin)`).
- **Status today**: **LIVE**.
- **Production callers**:
  - `orpheus/sn/__init__.py:2` — `from .operator import SNStreamingOperator`
    (re-exports as `orpheus.sn.SNStreamingOperator`; PUBLIC API).
  - `orpheus/sn/solver.py:48` — import.
  - `orpheus/sn/solver.py:234` — `self.L = SNStreamingOperator(...)` (in
    `SNSolver.__init__`; reachable from `solve_sn` → `SNSolver(...)`).
  - `orpheus/sn/solver.py:300` — `self.L = SNStreamingOperator(...)`
    (in `SNSolver.rebind_cross_sections`; reachable from `solve_sn`).
  - `orpheus/sn/solver.py:1434` — `matvec=solver.L.apply` (in
    `_solve_fixed_source_krylov`; reachable from
    `solve_sn_fixed_source` → krylov dispatch).
  - Plus 4 docstring references at `solver.py:20, 104, 143, 1407`
    (no code consumption — text only).
- **Test callers** (11 files; representative count of code-level
  usages, not just docstring mentions):
  - `tests/sn/test_snstreamingoperator.py` — 22+ test functions
    directly exercising `SNStreamingOperator(sn_mesh, sig_t)`
    (entire file is the dedicated suite; e.g. lines `139, 150, 160,
    482, 487, 499, 511`).
  - `tests/sn/test_phase_c_gates.py:37` import + 11 instantiation
    sites at `138, 159, 225, 261, 295, 346, 372, 411, 450, 489`.
  - `tests/sn/test_streaming_operator.py:51` import + 3 instantiation
    sites at `300, 333, 373` (the load-bearing
    `TestCompositionEquivalence`; gates the `(L+C) ≈
    SNStreamingOperator` xfail block).
  - `tests/sn/test_l1_standoff_slab_cylinder.py:92, 128, 132` —
    monkey-patches `sn_op.SNStreamingOperator.apply` to route through
    the unified matvec for L1 standoff verification.
  - `tests/sn/test_unified_matvec_slab.py:149, 175, 179` — same
    monkey-patch pattern as standoff.
  - `tests/sn/test_unified_matvec_cylinder.py:399, 427, 431` — same.
  - `tests/sn/test_b1pp_verification.py` — docstring-only at
    `12, 21` (the file's actual code uses
    `solution_to_angular_flux_with_traces`, not the legacy class).
  - `tests/sn/test_collision_operator.py:19, 95` — docstring-only
    (post-`e7d715e` migration to typed AngularFlux; no code use).
  - `tests/sn/test_streaming_operator_decomposition.py` — docstring
    only (uses the `_with_traces` family in code).
  - `tests/sn/spatial/test_psi_half_angle_seed.py:65, 359, 361` —
    docstring only.
  - `tests/sn/spatial/test_sweep_vs_apply_consistency.py:112` —
    docstring only.
  - `tests/numerics/test_iteration.py:532, 566` — code use:
    `L = SNStreamingOperator(sn_mesh=sn_mesh, ...)` in
    `test_keigenvalue_matches_solve_sn_2g_slab` (an L1 anchor that
    verifies the typed `KrylovAcceleration` pipeline against the
    legacy bundle).
- **Diagnostic callers**:
  - `derivations/diagnostics/diag_phase_g_step2_cyl_full_solve_with_fix.py:40, 50`
  - `derivations/diagnostics/diag_phase_g_step2_cyl_refinement_postfix.py:20`
  - `derivations/diagnostics/diag_phase_g_step2_cyl_residual_at_flat.py:20, 70`
  - `scratch/r1_step4_dry_run.py:27, 86`
- **Doc callers**: `docs/theory/discrete_ordinates.rst` (10 refs at
  `2124, 2682, 2695, 3382, 3925, 3961, 4018, 4329, 4471, 5021, 5030,
  5084`), `docs/theory/boundary_conditions.rst:2007`,
  `docs/theory/structured_geometry.rst:414`,
  `docs/theory/index_convention.rst:45`. All require narrative
  updates post-retirement (Phase 6 work).
- **Internal-only after retirement?**: NO — `SNStreamingOperator`
  is a public re-export from `orpheus.sn.__init__`. Removing it is a
  public-API change; the doc narrative needs explicit retirement
  notes per Issue #199.
- **Retirement order**: LAST (G3c) — depends on (a) production
  callers in `solver.py:234, 300, 1434` clearing (i.e.,
  `_solve_fixed_source_krylov` carved onto KrylovAcceleration AND
  `SNSolver.__init__`/`rebind_cross_sections` no longer building
  `self.L`); (b) `test_keigenvalue_matches_solve_sn_2g_slab`
  rewritten to compare against typed `(L+C)`; (c) the legacy
  `test_snstreamingoperator.py` file deleted; (d) the 3 monkey-patch
  fixtures in `test_l1_standoff_slab_cylinder.py` +
  `test_unified_matvec_{slab,cylinder}.py` retired or rewritten to
  patch the typed `StreamingOperator.apply` instead. Issue #199
  tracks this retirement as PR-TYPED-7.

### (2) `orpheus.sn.operator.EquationMap` (class)

- **Defined at**: `orpheus/sn/operator.py:119` (`@dataclass class
  EquationMap`).
- **Status today**: **LIVE** (heavily consumed by retirement-target
  family AND by `_with_traces` keepers).
- **Production callers**:
  - `orpheus/sn/operator.py:96` — `__all__` export
    (`"EquationMap"`).
  - `orpheus/sn/operator.py:220, 249, 264, 321, 457, 511, 526, 537,
    583, 677, 678, 737, 752, 768, 825, 841, 1398, 1465` — all
    internal type hints / annotations across the
    `build_equation_map*`, `solution_to_angular_flux*`,
    `*_with_traces`, `transport_operator_matvec*`, and
    `SNStreamingOperator._ensure_eq_map` family.
  - `orpheus/sn/solver.py:228, 1407` — docstring references only.
- **Test callers**:
  - `tests/sn/test_dag_walk.py:15` — docstring only ("the
    `EquationMap.unknowns_at_cell_for_mask` tests"; pinned via
    `build_equation_map*` callers below).
  - `tests/sn/test_streaming_operator_decomposition.py:115` —
    docstring only ("Geometry-dispatched EquationMap factory…").
  - `tests/sn/test_snstreamingoperator.py:478, 494, 506` — docstring
    references.
- **Diagnostic callers**: NONE that bypass the builder functions.
- **Doc callers**: `docs/theory/discrete_ordinates.rst:2629` and
  `docs/theory/index_convention.rst:434, 1456`. Narrative updates
  needed.
- **Internal-only after retirement?**: **NO — `EquationMap` is a
  KEEPER**, NOT a retirement target. It is consumed by
  `build_equation_map_with_traces` (a Step-G keeper) plus by the
  typed `StreamingOperator._ensure_eq_map` and
  `CollisionOperator._ensure_eq_map` (production typed-operator
  primitives, also keepers). After Phase 5 relocation,
  `EquationMap` migrates to `orpheus/sn/angular_flux.py` alongside
  the `_with_traces` helpers. **The brief's list inclusion of
  `EquationMap` as a retirement target is incorrect** — see
  SURPRISE-1.
- **Retirement order**: N/A — DO NOT retire. Relocate in Phase 5.

### (3) `orpheus.sn.operator.build_equation_map` (function, Cartesian)

- **Defined at**: `orpheus/sn/operator.py:218`.
- **Status today**: **LIVE**.
- **Production callers**:
  - `orpheus/sn/solver.py:49` — import.
  - `orpheus/sn/solver.py:1413` — call in `_solve_fixed_source_krylov`
    (the 2-D Cartesian / Cartesian-default branch:
    `eq_map = build_equation_map(nx, ny, sn_mesh.quad, ng)`).
  - `orpheus/sn/operator.py:1887` — call in
    `StreamingOperator._ensure_eq_map` for the **2-D Cartesian
    branch** (typed operator algebra path).
  - `orpheus/sn/operator.py:2222` — call in
    `CollisionOperator._ensure_eq_map` for the **2-D Cartesian
    branch** (typed operator algebra path).
- **Test callers**:
  - `tests/sn/test_snstreamingoperator.py:58, 234, 487` — used in
    `n_unknowns` consistency test + a Cartesian flat-flux test.
  - `tests/sn/test_streaming_operator_decomposition.py:58, 125` —
    used in `_eq_map_for` (returns `build_equation_map` for the 2-D
    branch).
  - `tests/sn/test_unified_matvec_slab.py:53, 262` — used in the
    slab matvec round-trip test.
- **Diagnostic callers**:
  - `derivations/diagnostics/diag_LC_decomposition_sn.py:64, 417`
  - `derivations/diagnostics/diag_LC_decomposition_resolution.py:47, 165`
  - `scratch/diag_unified_slab.py:11, 33`
- **Doc callers**: indirect via narrative; no `:func:` directive on
  the bare Cartesian variant.
- **Internal-only after retirement?**: NO. The 2-D Cartesian
  typed-operator branch at `operator.py:1887, 2222` is LOAD-BEARING —
  retiring `build_equation_map` strands the 2-D path of
  `StreamingOperator.apply` and `CollisionOperator.apply`. See
  SURPRISE-2.
- **Retirement order**: CANNOT retire in Step G as planned —
  requires either (a) deferring 2-D Cartesian to Phase A together
  with `_solve_fixed_source_krylov`'s Cartesian branch, OR (b)
  refactoring the typed operators to use the `_with_traces` layout
  for 2-D Cartesian as well (out of scope for Step G).

### (4) `orpheus.sn.operator.build_equation_map_spherical` (function)

- **Defined at**: `orpheus/sn/operator.py:509`.
- **Status today**: **LIVE**.
- **Production callers**:
  - `orpheus/sn/solver.py:50` — import.
  - `orpheus/sn/solver.py:1409` — call in `_solve_fixed_source_krylov`
    (spherical branch).
  - **NOT** called by the typed `StreamingOperator` /
    `CollisionOperator` `_ensure_eq_map` — those use
    `build_equation_map_with_traces` for the spherical (1-D) branch
    (operator.py:1894, 2225).
- **Test callers**:
  - `tests/sn/test_snstreamingoperator.py:59, 499, 552, 577, 618,
    662, 704, 769` — multiple spherical tests build the eq_map
    directly.
  - `tests/sn/test_unified_matvec_sphere.py:38, 139`.
  - `tests/sn/test_dag_walk.py` — via the indirect call from
    `SNStreamingOperator._ensure_eq_map`; the file does NOT import
    `build_equation_map_spherical` directly but is listed as a
    caller in Nexus because the eq_map construction routes through
    it.
  - `tests/sn/test_collision_operator.py` — via the indirect call
    only (TestSolve.test_solve_equals_division /
    TestSigmaLayout.test_localised_sigma_localised_output exercise
    operator construction).
- **Diagnostic callers**:
  - `derivations/diagnostics/diag_LC_decomposition_sn.py:65, 412, 564`
  - `derivations/diagnostics/diag_LC_decomposition_resolution.py:48,
    161, 231`
  - `derivations/diagnostics/diag_phase_g_step2_cyl_residual_at_flat.py:21, 63`
  - `scratch/step0_vcell_verification.py:33, 59`
- **Doc callers**: `docs/theory/discrete_ordinates.rst:5082`.
- **Internal-only after retirement?**: NO — it's not orphaned
  because the typed path uses `_with_traces` instead. It survives in
  production purely as `_solve_fixed_source_krylov`'s eq_map for the
  spherical branch.
- **Retirement order**: After `_solve_fixed_source_krylov` carves
  onto KrylovAcceleration (G1) — the eq_map dispatch dissolves with
  the carve. Test callers in `test_snstreamingoperator.py` retire
  with the file (`test_snstreamingoperator.py` deletes at G3c). The
  spherical-flat-flux + WDD-recurrence tests at
  `test_snstreamingoperator.py:618, 662, 704, 769` either retire or
  migrate to typed equivalents.

### (5) `orpheus.sn.operator.build_equation_map_cylindrical` (function — alias)

- **Defined at**: `orpheus/sn/operator.py:637` (module-level alias:
  `build_equation_map_cylindrical = build_equation_map_spherical`).
- **Status today**: **LIVE** (alias of #4).
- **Production callers**:
  - `orpheus/sn/solver.py:51` — import.
  - `orpheus/sn/solver.py:1411` — call in `_solve_fixed_source_krylov`
    (cylindrical branch).
- **Test callers**:
  - `tests/sn/test_snstreamingoperator.py:60, 511, 818`.
  - `tests/sn/test_unified_matvec_cylinder.py:47, 255, 309, 405,
    412` — the monkey-patched apply uses this for `_unified_apply`.
  - `tests/sn/spatial/test_apply_matvec_cylinder_invariants.py:56, 104`.
- **Diagnostic callers**:
  - `derivations/diagnostics/diag_LC_decomposition_sn.py:66, 414, 566`
  - `derivations/diagnostics/diag_LC_decomposition_resolution.py:49, 163, 233`
  - `derivations/diagnostics/diag_phase_g_step2_cyl_apply_internal.py:18, 60`
  - `derivations/diagnostics/diag_phase_g_step2_cyl_carlson_seed_fix.py:23, 125`
  - `derivations/diagnostics/diag_phase_g_step2_cyl_full_solve_with_fix.py:41, 48`
  - `derivations/diagnostics/diag_phase_g_step2_cyl_refinement_postfix.py:17, 61`
  - `derivations/diagnostics/diag_phase_g_step2_cyl_residual_at_flat.py:22, 65`
  - `derivations/diagnostics/diag_phase_g_step2_cyl_residual_pytest.py:26, 52`
  - `derivations/diagnostics/diag_step3_cyl_legacy_routing_bug.py:33, 220, 293`
  - `derivations/diagnostics/diag_step3_cyl_legacy_routing_LS6_outward.py:18, 48`
  - `derivations/diagnostics/diag_step3_cyl_reproduce.py:17, 317`
  - `derivations/diagnostics/diag_step3_cyl_unified_vs_hand_battery.py:26, 67`
- **Doc callers**: `docs/theory/discrete_ordinates.rst:5083`.
- **Internal-only after retirement?**: same shape as (4).
- **Retirement order**: same as (4) — after G1.

### (6) `orpheus.sn.operator.solution_to_angular_flux` (function, Cartesian)

- **Defined at**: `orpheus/sn/operator.py:262`.
- **Status today**: **LIVE**.
- **Production callers**:
  - `orpheus/sn/solver.py:52` — import.
  - `orpheus/sn/solver.py:724` — call in
    `SNSolver._make_sweep_preconditioner.matvec` (2-D Cartesian /
    Cartesian-default branch).
  - `orpheus/sn/solver.py:1519` — call in
    `_solve_fixed_source_krylov` (final-decode 2-D Cartesian
    branch).
  - `orpheus/sn/operator.py:1562` — call in
    `SNStreamingOperator.apply` (1-D slab fallback when
    `curv not in ("spherical", "cylindrical")`; note that
    `transport_operator_matvec_unified` consumes the decoder
    output).
  - `orpheus/sn/operator.py:923` — call in
    `_matvec_unified_cartesian` (Nexus reports this caller — the
    function lives inside `transport_operator_matvec_unified` /
    related helpers).
- **Test callers**:
  - `tests/sn/test_unified_matvec_slab.py:54, 160, 276, 310`.
  - `tests/sn/test_l1_standoff_slab_cylinder.py:65, 112`.
- **Diagnostic callers**:
  - `scratch/diag_unified_slab.py:10, 41`.
- **Doc callers**: `docs/theory/index_convention.rst:466, 1130,
  1460`, `docs/theory/boundary_conditions.rst:2215`,
  `docs/theory/discrete_ordinates.rst:2133, 5909`.
- **Internal-only after retirement?**: NO — also wired into the
  typed `SNStreamingOperator.apply` 1-D slab path AND into
  `_matvec_unified_cartesian`. Retirement strands the legacy
  `SNStreamingOperator.apply` slab path and any consumer of the
  Cartesian 2-D `_solve_fixed_source_*`.
- **Retirement order**: After G1 (`_solve_fixed_source_krylov`
  carve) AND after G3c (`SNStreamingOperator` retired). The
  `_matvec_unified_cartesian` caller is the surprise — see
  SURPRISE-3.

### (7) `orpheus.sn.operator.solution_to_angular_flux_spherical` (function)

- **Defined at**: `orpheus/sn/operator.py:535`.
- **Status today**: **LIVE**.
- **Production callers**:
  - `orpheus/sn/solver.py:53` — import.
  - `orpheus/sn/solver.py:716` — call in
    `SNSolver._make_sweep_preconditioner.matvec` (spherical branch).
  - `orpheus/sn/solver.py:1511` — call in
    `_solve_fixed_source_krylov` (final-decode spherical branch).
  - `orpheus/sn/operator.py:1558` — call in
    `SNStreamingOperator.apply` (1-D curvilinear branch; also
    handles cylindrical via the `curv in ("spherical",
    "cylindrical")` test).
  - `orpheus/sn/operator.py:1018+` — call in
    `_matvec_unified_spherical` (Nexus-reported caller; lives in
    `transport_operator_matvec_unified` machinery).
- **Test callers**:
  - `tests/sn/test_snstreamingoperator.py:548, 555, 573, 581` —
    direct exerciser tests.
  - `tests/sn/test_unified_matvec_sphere.py:16, 20, 76` (docstring +
    code paths via decoder).
- **Diagnostic callers**:
  - `scratch/derivations/diagnostics/diag_issue168_02_option_a_dd_extrap.py`
    (Nexus result; patched spherical matvec).
  - `scratch/derivations/diagnostics/diag_issue168_04_option_a_vacuum.py`.
- **Doc callers**: `docs/theory/discrete_ordinates.rst:2204, 2669`.
- **Internal-only after retirement?**: NO — same shape as (6).
- **Retirement order**: same — after G1 AND G3c.

### (8) `orpheus.sn.operator.solution_to_angular_flux_cylindrical` (function — alias)

- **Defined at**: `orpheus/sn/operator.py:638` (alias of #7).
- **Status today**: **LIVE** (alias).
- **Production callers**:
  - `orpheus/sn/solver.py:54` — import.
  - `orpheus/sn/solver.py:720` — `_make_sweep_preconditioner.matvec`
    cylindrical branch.
  - `orpheus/sn/solver.py:1515` — `_solve_fixed_source_krylov`
    final-decode cylindrical branch.
  - (`SNStreamingOperator.apply:1558` handles cylindrical via the
    spherical decoder due to alias.)
- **Test callers**:
  - `tests/sn/test_unified_matvec_cylinder.py:48, 412` — monkey-patch
    apply uses this for `_unified_apply`.
  - `tests/sn/test_l1_standoff_slab_cylinder.py:66, 103`.
- **Diagnostic callers**:
  - `derivations/diagnostics/diag_step3_cyl_legacy_routing_LS6_outward.py:19, 66`.
  - `derivations/diagnostics/diag_step3_cyl_legacy_routing_bug.py:34, 239`.
  - `derivations/diagnostics/diag_step3_cyl_reproduce.py:18, 336`.
  - `derivations/diagnostics/diag_phase_g_step2_cyl_apply_internal.py:19,
    62, 126`.
  - `derivations/diagnostics/diag_phase_g_step2_cyl_full_solve_with_fix.py:42, 60, 81`.
  - `derivations/diagnostics/diag_phase_g_step2_cyl_refinement_postfix.py:19`.
  - `derivations/diagnostics/diag_phase_g_step2_cyl_residual_at_flat.py:117, 118`.
  - `derivations/diagnostics/diag_krylov_iter_breakdown.py:43`.
  - `derivations/diagnostics/diag_issue_197_cyl_twin_path_decoder_mismatch.py:2, 179`.
- **Doc callers**: indirect.
- **Internal-only after retirement?**: same as (7).
- **Retirement order**: same.

### (9) `orpheus.sn.solver.SNSolver._make_sweep_preconditioner` (method)

- **Defined at**: `orpheus/sn/solver.py:680`.
- **Status today**: **LIVE**.
- **Production callers**:
  - `orpheus/sn/solver.py:1436` — call in
    `_solve_fixed_source_krylov`: `precond =
    solver._make_sweep_preconditioner(eq_map, n_unknowns)`. **THIS
    IS THE ONLY PRODUCTION CALLER.**
- **Test callers**: NONE (Nexus reports zero; verified with grep).
- **Diagnostic callers**:
  - `derivations/diagnostics/diag_krylov_iter_breakdown.py:118` —
    `precond = solver._make_sweep_preconditioner(eq_map, n_unknowns,
    sum_w)`. Note: the diagnostic passes a 3-arg signature
    `(eq_map, n_unknowns, sum_w)` whereas production at
    `solver.py:680` has 2 args `(eq_map, n: int)` — stale
    diagnostic (a pre-A1 surface).
  - `scratch/derivations/diagnostics/diag_issue168_05_apply_vs_solve.py`
    (referenced via grep; not via Nexus).
- **Doc callers**: none direct.
- **Internal-only after retirement?**: YES — once
  `_solve_fixed_source_krylov` carves onto KrylovAcceleration (Phase
  3 / G1), this method is orphaned and retires cleanly.
- **Retirement order**: G3a (after G1). Mechanical delete; no test
  migration cost.

### (10) `orpheus.sn.solver._build_rhs_cartesian` (free function)

- **Defined at**: `orpheus/sn/solver.py:839`.
- **Status today**: **LIVE**.
- **Production callers**:
  - `orpheus/sn/solver.py:1481` — call in `_solve_fixed_source_krylov`
    (Cartesian/2-D branch).
  - Docstring references at `solver.py:873, 1399, 1419, 1475` (no
    code consumption).
- **Test callers**: NONE.
- **Diagnostic callers**: NONE direct (the cylindrical alias is
  exercised; see (12)).
- **Doc callers**: none direct.
- **Internal-only after retirement?**: YES — once G1 ships,
  orphaned.
- **Retirement order**: G3a (after G1). Issue #174 (filed
  2026-05-10) is the existing open ticket against this function:
  *"sn: refactor `_build_rhs_cartesian` Krylov RHS via
  ScatteringOperator (Issue #162 close-out)"* — Step G subsumes its
  scope by retiring the function entirely.

### (11) `orpheus.sn.solver._build_rhs_spherical` (free function)

- **Defined at**: `orpheus/sn/solver.py:938`.
- **Status today**: **LIVE**.
- **Production callers**:
  - `orpheus/sn/solver.py:1463` — call in `_solve_fixed_source_krylov`
    (spherical branch).
  - `orpheus/sn/solver.py:1469` — call in `_solve_fixed_source_krylov`
    (cylindrical branch, via the `_build_rhs_cylindrical` alias).
  - Docstring reference at `solver.py:950` (no code consumption).
- **Test callers**: NONE.
- **Diagnostic callers**:
  - `derivations/diagnostics/diag_krylov_iter_breakdown.py:40, 87` —
    used as `# cylindrical alias` in a Krylov-breakdown probe.
- **Doc callers**: none direct.
- **Internal-only after retirement?**: YES (post G1).
- **Retirement order**: G3a (after G1). Issue #174 also applies.

### (12) `orpheus.sn.solver._build_rhs_cylindrical` (module-level alias)

- **Defined at**: `orpheus/sn/solver.py:982` (`_build_rhs_cylindrical
  = _build_rhs_spherical`).
- **Status today**: **LIVE** (alias of #11).
- **Production callers**:
  - `orpheus/sn/solver.py:1469` — call in `_solve_fixed_source_krylov`
    (cylindrical branch).
- **Test callers**: NONE.
- **Diagnostic callers**: NONE direct (it's an alias; the diagnostic
  imports `_build_rhs_spherical` and uses it as cylindrical).
- **Doc callers**: none direct.
- **Internal-only after retirement?**: YES.
- **Retirement order**: G3a, together with #11.

---

## Retirement sequence — explicit ordering

This is the load-bearing output. Each step's predecessors must
complete (commit landed AND tests green) before the step's code is
touched. The numbering corresponds 1:1 to the session-2 plan's
phase boundaries; sub-step letters are scoped to this audit.

### Phase 3 (G1 / G2) — carve fixed-source onto typed-operator algebra

1. **G1 — carve `_solve_fixed_source_krylov` onto
   `KrylovAcceleration` + typed `AngularFlux`**
   (`orpheus/sn/solver.py:1374`).
   - Mirrors Step D's carve of `SNSolver._solve_krylov`.
   - **Orphans** (after green): `_build_rhs_cartesian`,
     `_build_rhs_spherical`, `_build_rhs_cylindrical`,
     `_make_sweep_preconditioner`. They lose their sole production
     caller.
   - **Does NOT orphan**: `build_equation_map*`,
     `solution_to_angular_flux*`, or `solver.L`/`SNStreamingOperator`
     — those still have other consumers (see §G2 / §G3c).
   - **2-D Cartesian decision point (SURPRISE-4)**: G1's typed path
     raises `NotImplementedError` on 2-D Cartesian (same as
     `_solve_krylov` post-Step D). The Cartesian-default
     fixed-source path either (a) defers 2-D Cartesian to Phase A by
     KEEPING `_solve_fixed_source_krylov` AND `_build_rhs_cartesian`
     alive for 2-D-only, OR (b) routes 2-D Cartesian to
     `_solve_fixed_source_si` (which is geometry-agnostic — see G2
     and SURPRISE-5). Decision must be made BEFORE G1 commits.

2. **G2 — carve `_solve_fixed_source_si` onto `SourceIteration` +
   typed `AngularFlux`** (`orpheus/sn/solver.py:1278`).
   - Mirrors Step E's carve of `SNSolver._solve_source_iteration`.
   - **Pre-condition**: A1 has shipped (commit `e8a40ee` /
     `de8822d`; see Phase 1.1) so `external_source` is per-ordinate
     density at the API boundary AND `ScatteringOperator.apply`
     normalises at the producer.
   - **Orphans**: nothing in the retirement list (the SI path
     doesn't touch any of the 12 targets). G2's job is to
     unblock G3c by removing the test-file dependency at
     `test_keigenvalue_matches_solve_sn_2g_slab`.

### Phase 4 (G3) — retire the legacy symbols

3. **G3a — retire `_make_sweep_preconditioner` +
   `_build_rhs_cartesian` + `_build_rhs_spherical` +
   `_build_rhs_cylindrical`**.
   - Free functions and SNSolver method. Zero production callers
     after G1. Zero test callers. Mechanical delete.
   - Closes Issue #174 (no longer relevant — function deleted).
   - Diagnostic at `diag_krylov_iter_breakdown.py:40, 87, 118` goes
     stale; either delete the file or add a header note "exercises
     retired symbols; preserved as historical Krylov-breakdown
     probe".

4. **G3b — migrate / retire the legacy decoder + builder family**:
   `build_equation_map_spherical`, `build_equation_map_cylindrical`,
   `solution_to_angular_flux_spherical`,
   `solution_to_angular_flux_cylindrical`.
   - Production: after G1, the only remaining production callers
     are `SNStreamingOperator.apply` (operator.py:1558) and
     `SNStreamingOperator._ensure_eq_map`
     (`build_equation_map_spherical` for the curvilinear branch via
     module-level lookup at `operator.py:1465`+) plus the
     `_matvec_unified_spherical` helper. ALL retire with
     `SNStreamingOperator` at G3c.
   - Test callers: `test_snstreamingoperator.py` deletes at G3c;
     `test_unified_matvec_cylinder.py` + `test_l1_standoff_slab_cylinder.py`
     + `test_unified_matvec_sphere.py` need the monkey-patch
     fixtures rewritten OR retired (the typed
     `StreamingOperator.apply` is the new patch target — and post
     G3c there's nothing to compare against, so most of these tests
     either retire or change semantics).
   - Diagnostics under `derivations/diagnostics/diag_step3_cyl_*`
     and `diag_phase_g_step2_cyl_*` exercise the legacy decoder
     paths for historical reproduction; preserve as-is with a
     header note (they're frozen-in-time investigations).

5. **G3c — retire `build_equation_map` (Cartesian)** and the
   Cartesian-1D fallback path inside `solution_to_angular_flux`.
   - **BLOCKED** unless 2-D Cartesian deferral or migration is
     decided. The typed `StreamingOperator._ensure_eq_map` and
     `CollisionOperator._ensure_eq_map` BOTH call
     `build_equation_map(nx, ny, quad, ng)` on the `curv is None and
     ny > 1` branch (operator.py:1887, 2222). Retiring
     `build_equation_map` strands the entire typed 2-D Cartesian
     path. See SURPRISE-2.

6. **G3d — retire `SNStreamingOperator` class**.
   - **Pre-conditions**:
     - G1 + G2 shipped (no production caller via
       `_solve_fixed_source_krylov` or `_solve_fixed_source_si`).
     - `SNSolver.__init__:234` AND `rebind_cross_sections:300` no
       longer build `self.L = SNStreamingOperator(...)`. Either
       (a) `self.L` removed entirely (Issue #199's preferred path),
       or (b) `self.L = StreamingOperator(...) + CollisionOperator(...)`
       (an `OperatorSum`).
     - `test_keigenvalue_matches_solve_sn_2g_slab` rewritten to
       compare against typed `(L+C)` instead of legacy bundle.
     - `tests/sn/test_snstreamingoperator.py` deleted in entirety
       (superseded by `test_streaming_operator.py` +
       `test_collision_operator.py`).
     - `tests/sn/test_l1_standoff_slab_cylinder.py` +
       `test_unified_matvec_{slab,cylinder,sphere}.py` monkey-patch
       fixtures retired or migrated (the patch target is gone).
     - `tests/sn/test_phase_c_gates.py` — 11 instantiation sites
       migrated to typed `StreamingOperator(sn_mesh, sig_t)`.
     - `tests/sn/test_streaming_operator.py::TestCompositionEquivalence`
       xfail block deleted (no `SNStreamingOperator` to compare
       against — closes the `xfail strict` rows at `300, 333, 373`).
   - Public-API change: `orpheus.sn.SNStreamingOperator` import
     fails. Doc narrative (`docs/theory/discrete_ordinates.rst` —
     10+ refs) needs Phase 6 rewrite.
   - **Closes Issue #199** (the omnibus retirement ticket).

7. **G3e — retire `EquationMap` itself? NO.**
   - SURPRISE-1: `EquationMap` is a KEEPER (it's used by the
     `_with_traces` family, which are themselves keepers AND by the
     typed `StreamingOperator._ensure_eq_map` /
     `CollisionOperator._ensure_eq_map`). It does NOT retire in
     Step G. Phase 5 relocates it alongside the `_with_traces`
     helpers.

### Phase 5 (G4) — relocate keepers

8. **G4 — relocate `pack_with_traces`,
   `solution_to_angular_flux_with_traces`,
   `build_equation_map_with_traces`, AND `EquationMap`** from
   `orpheus/sn/operator.py` to `orpheus/sn/angular_flux.py`.
   - Out of scope for Step G proper; Phase 5 work. Audit
     question C below sizes the effort.

---

## Quick-reference callgraph diagram

The retirement-target chain rooted at `solve_sn_fixed_source`:

```
solve_sn_fixed_source                       (solver.py:1127)
├── _solve_fixed_source_si                  (solver.py:1278)
│      [DOES NOT touch any retirement target — uses transport_sweep + typed AngularFlux]
│      [G2 carves this path onto SourceIteration; geometry-agnostic 1-D + 2-D OK]
│
└── _solve_fixed_source_krylov              (solver.py:1374)        [G1 target]
    │
    ├── build_equation_map / _spherical / _cylindrical    (#3 / #4 / #5)
    │       [orpheus/sn/solver.py:1409, 1411, 1413]
    │
    ├── solver.L.apply  (SNStreamingOperator.apply)       (#1)
    │       [orpheus/sn/solver.py:1434; matvec=...]
    │       └── solution_to_angular_flux_spherical                (#7)
    │       └── solution_to_angular_flux                          (#6)
    │       └── transport_operator_matvec_unified (keeper)
    │
    ├── solver._make_sweep_preconditioner   (#9)          [G3a]
    │       [orpheus/sn/solver.py:1436]
    │       └── solution_to_angular_flux_spherical                (#7)
    │       └── solution_to_angular_flux_cylindrical              (#8)
    │       └── solution_to_angular_flux                          (#6)
    │       └── transport_sweep (keeper)
    │
    ├── _build_rhs_spherical                (#11)         [G3a]
    │       [orpheus/sn/solver.py:1463]
    ├── _build_rhs_cylindrical              (#12)         [G3a]
    │       [orpheus/sn/solver.py:1469]
    ├── _build_rhs_cartesian                (#10)         [G3a]
    │       [orpheus/sn/solver.py:1481]
    │
    └── solution_to_angular_flux* (final decode)
            [orpheus/sn/solver.py:1511, 1515, 1519]


SNSolver.__init__  (solver.py:127, line 234)              [G3c blocker]
  └── self.L = SNStreamingOperator(sn_mesh, sig_t)        (#1)

SNSolver.rebind_cross_sections (solver.py:274, line 300)  [G3c blocker]
  └── self.L = SNStreamingOperator(sn_mesh, new_sig_t)    (#1)
```

The typed-operator algebra (StreamingOperator, CollisionOperator,
InvertibleOperator) calls `build_equation_map` (#3) ONLY on the 2-D
Cartesian branch — a SURPRISE consumer.

---

## SURPRISE entries

### SURPRISE-1 — `EquationMap` is a KEEPER, NOT a retirement target

The brief's list includes `EquationMap` (class) at position 2 as a
Step G retirement target. The audit shows `EquationMap` is consumed
by:

1. The `_with_traces` family (`build_equation_map_with_traces` at
   `operator.py:677` returns `EquationMap`; the keepers Phase 5
   relocates).
2. The typed `StreamingOperator._ensure_eq_map`
   (`operator.py:1870`).
3. The typed `CollisionOperator._ensure_eq_map`
   (`operator.py:2209`).
4. All `transport_operator_matvec_*` helpers (some retire, some
   stay).

**`EquationMap` MUST remain after Step G**. The right action is
Phase 5 relocation alongside the `_with_traces` helpers. The brief
should be amended.

### SURPRISE-2 — 2-D Cartesian typed operators consume `build_equation_map`

`StreamingOperator._ensure_eq_map` at `operator.py:1887` and
`CollisionOperator._ensure_eq_map` at `operator.py:2222` BOTH call
the bare `build_equation_map(nx, ny, quad, ng)` on the 2-D
Cartesian (`curv is None and ny > 1`) branch. They use the LEGACY
non-`_with_traces` layout because 2-D anti-diagonal wavefront sweeps
haven't migrated to `SNMesh.dag_walk` (per Issue #199's enumeration
of the 2-D Cartesian absorption deferral).

Implication: **`build_equation_map` (Cartesian, #3) CANNOT retire in
Step G**. It either (a) defers retirement until 2-D Cartesian
absorbs (Phase A; out of Step G scope), or (b) the typed operators
migrate to `build_equation_map_with_traces` for 2-D too (out of Step
G scope; requires separate B1''-aware 4-face layout). The session-2
plan §3-§4 must document this deferral explicitly.

### SURPRISE-3 — `solution_to_angular_flux*` are consumed by `transport_operator_matvec_unified` helpers

Nexus reports two non-obvious production callers:

- `_matvec_unified_cartesian` (orpheus/sn/operator.py around line
  920+) calls `solution_to_angular_flux` (#6).
- `_matvec_unified_spherical` (orpheus/sn/operator.py around line
  1000+) calls `solution_to_angular_flux_spherical` (#7).

The `transport_operator_matvec_unified` family is the LIVE unified
matvec engine that BOTH `SNStreamingOperator.apply` AND
`StreamingOperator.apply` route through. Retiring the decoders
strands these helpers UNLESS the helpers are refactored to consume
typed `AngularFlux` directly OR the decoder fallback is preserved
as a private helper inside `transport_operator_matvec_unified`.

This is the load-bearing surprise: the `*_with_traces` family alone
is NOT a clean replacement — the unified matvec specifically uses
the legacy decoders for the cell-centre-proxy variant
(SNStreamingOperator's contract). Step G3c must rewire these or
preserve the decoders as `_with_traces` wrappers.

### SURPRISE-4 — 2-D Cartesian fixed-source Krylov path will need deferral

`_solve_krylov` and `_solve_source_iteration` (post-R-1 Step D/E)
both raise `NotImplementedError` on `self.sn_mesh.reduced is None`
(i.e., 2-D Cartesian). The session-2 plan §3.1 (G1) follows the
same shape — the typed `KrylovAcceleration` carve is 1-D-only.

Consequence: after G1 ships, `_solve_fixed_source_krylov` for 2-D
Cartesian has NO typed-path target. Either:

1. Defer 2-D Cartesian fixed-source Krylov to Phase A — KEEP the
   legacy `_solve_fixed_source_krylov` + `_build_rhs_cartesian` +
   `_make_sweep_preconditioner` alive for the 2-D Cartesian branch
   only. This delays G3a's "mechanical delete" until Phase A.
2. Route 2-D Cartesian to `_solve_fixed_source_si` (G2's geometry-
   agnostic carve). Possible: the SI path uses `transport_sweep`
   which supports 2-D wavefront natively. This is the cleanest
   structural option — see SURPRISE-5.
3. Raise `NotImplementedError` for 2-D Cartesian + Krylov inner;
   force users onto SI. Unprincipled — breaks a previously-working
   public API.

Decision required at G1 design.

### SURPRISE-5 — `_solve_fixed_source_si` is geometry-agnostic and the natural 2-D landing zone

Unlike `_solve_krylov`/`_solve_source_iteration` (the eigenvalue
inners, which raise on 2-D), `_solve_fixed_source_si` at
`solver.py:1278-1372` has NO 2-D guard — it uses
`transport_sweep(source, sig_t, sn_mesh, boundary)` which dispatches
to the appropriate sweep (`_sweep_1d_unified` or
`_sweep_2d_wavefront`) based on mesh shape. After G2 ships, the SI
path naturally handles 2-D Cartesian.

This is the principled answer to SURPRISE-4: when 2-D Cartesian
hits the Krylov dispatch in `solve_sn_fixed_source:1262`,
`_solve_fixed_source_krylov` can `NotImplementedError` (or fall
through to SI for the 2-D-Cartesian-only sub-case). The Phase A
absorption then closes the 2-D Krylov gap structurally.

### SURPRISE-6 — Diagnostics import a stale `_make_sweep_preconditioner` signature

`derivations/diagnostics/diag_krylov_iter_breakdown.py:118` invokes
`solver._make_sweep_preconditioner(eq_map, n_unknowns, sum_w)` (3
arguments). Production at `solver.py:680` has 2 arguments
`(eq_map, n: int)`. The diagnostic is stale (pre-A1 surface) and
would TypeError if run today. No action required (the diagnostic
will simply break on next run; retirement of the function removes
the surface anyway). Flagged so the implementer doesn't get
distracted hunting a "regression" in the diagnostic.

### SURPRISE-7 — `test_phase_c_gates.py` has 11 production-level uses of `SNStreamingOperator`

The test file at `tests/sn/test_phase_c_gates.py` is NOT in
`test_snstreamingoperator.py`'s retirement scope — it's a
phase-C-gate suite that builds `SNStreamingOperator(sn_mesh,
sig_t)` instances 11 times to validate apply / apply_transpose /
reciprocity / face-flux invariants. These tests pin the LEGACY
bundle's per-ordinate flat-flux residual, BC-trace contract, etc.

For G3c to retire `SNStreamingOperator`, these 11 tests need to:

(a) **Migrate**: rewrite to use `StreamingOperator + CollisionOperator`
    (typed sum) — most invariants still hold structurally.
(b) **Retire**: if the invariants are duplicated in
    `test_streaming_operator.py` / `test_collision_operator.py`,
    delete; otherwise migrate.
(c) **Reframe**: some tests (face-flux at `i = nx-1`, cell-centre
    proxy at the Carlson seed) pin behaviour that the typed
    `(L+C)` path does NOT have (B1'' uses the actual face-state
    instead of the proxy). These tests retire as they pin a
    soon-to-be-retired contract.

The Issue #199 body mentions migrating `test_snstreamingoperator.py`
but not `test_phase_c_gates.py` explicitly. Audit recommends
expanding the issue body before G3c starts.

### SURPRISE-8 — Doc narrative scope is larger than session-2 plan §6 anticipated

`docs/theory/discrete_ordinates.rst` has **10 references** to
`SNStreamingOperator` at lines `2124, 2682, 2695, 3382, 3925, 3961,
4018, 4329, 4471, 5021, 5030, 5084` plus the equation-map narrative
sections. The Phase 6 / H Sphinx narrative MUST address each — not
just write a new "typed operator algebra" page. The mention count
suggests Phase 6 is closer to ~3-4 hours of careful narrative
rewriting (the plan's estimate) BUT no less.

---

## Must-migrate caller summary (counts)

* **Production**: 4 sites for SNStreamingOperator (incl. public
  `__init__.py`), 1 site for `_make_sweep_preconditioner`, 3 sites
  for `_build_rhs_*`, 6 sites for `solution_to_angular_flux*`
  (including 2 typed-helper consumers — SURPRISE-3), 4 sites for
  `build_equation_map*` (including 2 typed-operator consumers —
  SURPRISE-2). **Total: 18 production sites across `solver.py` and
  `operator.py`.**
* **Tests**: 11 files reference `SNStreamingOperator` (≈108 line
  occurrences total). Of those, ~5 are docstring-only and survive
  Step G; ~6 need active migration or deletion at G3c.
* **Diagnostics**: ~15 files under `derivations/diagnostics/`
  exercise retirement targets directly. Most are
  frozen-in-time investigations; preserve with a header note.
* **Docs**: 10+ `discrete_ordinates.rst` references, 3
  `index_convention.rst`, 1 `boundary_conditions.rst`, 1
  `structured_geometry.rst`. Phase 6 work.

Bottom line: Step G is **as planned for the 10 cleanly-retiring
targets (#7-#12 and `solver.L` wiring)**, but **LARGER than planned
by 2 entries** (`EquationMap` is a keeper not a retirement target;
`build_equation_map` Cartesian cannot retire in Step G due to 2-D
typed-operator consumers). The retirement sequence (8 numbered
steps above) honours these constraints. Issue #199 already
anticipates a separable PR-TYPED-7 for the 2-D Cartesian absorption
— Step G's job is to ship the 1-D scope cleanly and stop at the
2-D boundary.

---

## Augmentation responses

### A. `SNStreamingOperator.apply` vs typed `StreamingOperator.apply` — algebraically equivalent on current tree?

**Status**: principled-equivalent in 1-D curvilinear (sphere + cylinder)
since PR-TYPED-6c Step 5; **structurally NON-equivalent** in 1-D slab and
2-D Cartesian.

Code-level reading:

- Both call `transport_operator_matvec_unified` for the 1-D path
  (`operator.py:1569` for SN, `operator.py:2020` for typed).
- **Difference 1 (1-D slab)**: typed `StreamingOperator.apply` uses
  `solution_to_angular_flux_with_traces` (B1'' face-aware,
  `operator.py:2011`) and passes the actual face state via
  `psi_face_outer / psi_face_inner`. Legacy
  `SNStreamingOperator.apply` uses `solution_to_angular_flux`
  (Cartesian, no face) and `transport_operator_matvec_unified`
  falls back to the cell-centre proxy at the Carlson seed. Result:
  **FD 1st-order (legacy) ≠ WDD 2nd-order (typed)** — documented
  order-of-accuracy delta. Pinned at
  `tests/sn/test_streaming_operator.py:319-373` as `xfail strict`
  (TestCompositionEquivalence cylinder + slab branches).
- **Difference 2 (cylinder)**: typed uses face-aware
  `psi_face_outer` (B1''), legacy uses cell-centre proxy. Documented
  as "cylinder twin-path bug (rel ≈ 4e-3 at nx=40)" in the
  `xfail strict` block at
  `test_streaming_operator.py:286-309` and in
  `test_l1_standoff_slab_cylinder.py::test_cylinder_l1_sweep_vs_krylov_twin_path`.
- **Difference 3 (sphere)**: bit-identical post Step 5 (the
  cell-centre proxy and the B1'' face coincide because the inward
  hemisphere at the outer face IS the M-M coupled-pole seed). The
  `test_streaming_operator.py::TestCompositionEquivalence.test_uniform_sigma_t_homogeneous`
  passes WITHOUT `xfail` on sphere.
- **Difference 4 (2-D Cartesian)**: typed at `operator.py:1995`
  calls `transport_operator_matvec` (FD) with `sigma_packed * psi`
  subtraction at the algebra layer; legacy at `operator.py:1540`
  calls `transport_operator_matvec` WITHOUT the
  subtraction — `SNStreamingOperator` IS the bundle (M = L+C), the
  typed `StreamingOperator` is the leaf L = M - σ·ψ.

**Implication for Phase 3 (G1/G2)**: the carve from
`solver.L.apply` (SNStreamingOperator) to typed `(L+C)`
swaps the cell-centre-proxy contract for the B1''
face-aware contract. The L1 anchor
`test_keigenvalue_matches_solve_sn_2g_slab` (the only test that
compares typed Krylov against legacy bundle in a fully-converged
mode) will likely shift by an O(h²) - O(h) delta on cylinder. If
the L1 anchor was tight on the legacy answer, it will need
re-baselining. Issue #199 documents the cylinder twin-path tests
flipping from `xfail strict` to PASSING after G3c — this is the
intended improvement, not a regression.

### B. 2-D Cartesian retirement scope

Yes — three retirement targets have 2-D-Cartesian-specific bodies
that would be stranded post-Step-G:

1. **`_build_rhs_cartesian`** (#10, `solver.py:839`) — the
   Cartesian implementation supports both 1-D slab (ny=1) and 2-D
   (ny>1). Retirement deletes both.
2. **`solution_to_angular_flux`** (#6, `operator.py:262`) — same
   shape; the 2-D case at `solver.py:1519` (final decode in
   `_solve_fixed_source_krylov`) and the typed
   `StreamingOperator.apply` 2-D branch at `operator.py:1995-2007`
   are 2-D-load-bearing.
3. **`build_equation_map`** (#3, `operator.py:218`) — SURPRISE-2:
   the typed `StreamingOperator._ensure_eq_map` and
   `CollisionOperator._ensure_eq_map` consume it for 2-D Cartesian.

**Recommendation**: defer 2-D Cartesian retirement to Phase A (per
Issue #199 §5). For Step G:

- Either KEEP `_build_rhs_cartesian` alive for 2-D-only (and route
  2-D Cartesian through `_solve_fixed_source_krylov`'s preserved
  2-D path), OR force 2-D Cartesian to `_solve_fixed_source_si`
  (which is 2-D-natural via `transport_sweep`'s wavefront
  dispatch — SURPRISE-5).
- KEEP `build_equation_map` alive (its 2-D consumers are
  load-bearing for the typed `StreamingOperator`/`CollisionOperator`
  on 2-D Cartesian).
- KEEP `solution_to_angular_flux` Cartesian-2-D alive for the same
  reasons.

These three keepers retire alongside the 2-D Cartesian absorption
to `SNMesh.dag_walk` in PR-TYPED-7 / Phase A. Document the deferral
in the Step G commit messages.

### C. `_with_traces` keepers — relocation cost (Phase 5 sizing)

Call-site counts (verified via `grep -c`):

* `orpheus/sn/operator.py`: **17 references** (defines + internal
  uses + docstrings)
* `orpheus/sn/angular_flux.py`: **8 references** (consumer; the
  Phase 5 destination)
* `tests/sn/test_streaming_operator_decomposition.py`: **8
  references**
* `tests/sn/test_angular_flux_with_boundary.py`: **6 references**
* `tests/sn/test_b1pp_verification.py`: **4 references**
* `tests/numerics/test_iteration_angular_flux.py`: **2 references**
* `tests/sn/test_operators_apply_typed.py`: **2 references**
* `tests/sn/test_collision_operator.py`: **1 reference**

**Total: 48 call sites across 8 files** for the three
`_with_traces` helpers + `EquationMap` (which migrates with them).

Phase 5 relocation cost: low. The relocation is mechanical
(`from .operator import ...` → `from .angular_flux import ...`)
plus updating `__all__` and the
`orpheus/sn/angular_flux.py:379-426` internal imports (which
currently use a local `from .operator import ...` pattern — relocates
to module-level imports). Estimate: 30-45 minutes of edits + 1
Sphinx rebuild.

### D. GitHub issues against the retirement scope

Three OPEN issues mention retirement targets by name:

1. **Issue #199** — *"PR-TYPED-7 — Retire SNStreamingOperator;
   migrate solve_sn to (L+C) algebra"* (filed 2026-05-20).
   - This IS the Step G omnibus issue. Body enumerates:
     `solve_sn` migration (item 1), `SNStreamingOperator` class
     retirement (item 2), `build_equation_map_cylindrical` +
     `solution_to_angular_flux_cylindrical` aliases (item 3),
     L14 verification flip-green (item 4), 2-D Cartesian absorption
     (item 5, separable).
   - Item 5 explicitly anticipates SURPRISE-2 / SURPRISE-4: *"2-D
     anti-diagonal wavefront sweeps haven't migrated to
     `SNMesh.dag_walk`. Either migrate `dag_walk` to handle 2-D
     anti-diagonals OR write a separate 2-D unified matvec.
     Independent of (1)-(4); can defer to a separate PR if scope
     is too large."*
   - Implementer action: Step G closes items 1-4 of #199; defer
     item 5 (2-D Cartesian) to Phase A.

2. **Issue #174** — *"sn: refactor `_build_rhs_cartesian` Krylov
   RHS via ScatteringOperator (Issue #162 close-out)"* (filed
   2026-05-10).
   - Predates Step G; proposes a softer refactor (route
     `_build_rhs_cartesian` through `ScatteringOperator.apply`).
   - **Step G subsumes this issue** — G3a deletes
     `_build_rhs_cartesian` entirely (and its spherical /
     cylindrical sister aliases). Close #174 with reference to
     G3a's commit when it lands.

3. **Issue #160** — *"sn: SNStreamingOperator as LinearOperator"*
   (filed 2026-05-06).
   - The ORIGINATING issue for `SNStreamingOperator`'s shipped
     form. Open since pre-R-1; Step G3c retirement is the
     terminal close-out. Body documents the rationale for the
     bundle's `apply` + `solve` + `apply_transpose` capabilities —
     useful Phase 6 Sphinx-narrative context.

No other open issues mention `_make_sweep_preconditioner`,
`build_equation_map`, `solution_to_angular_flux`, or `EquationMap`
by name. (`#205` and `#168` mention adjacent SN topics but are
out-of-scope for Step G.)
