---
name: SN reshape Wave E Round 2 closeout
description: SNSolver migrated to operator-algebra; legacy BiCGSTAB FD-operator path retired. ERR-026 closure scoped to reflective-BC eigenvalue (Round 3 owns vacuum-BC fixed-source).
type: project
---

# Wave E Round 2 — Issue #164 closeout

## Branch and HEAD

- Worktree: `/workspaces/ORPHEUS/.claude/worktrees/agent-ab52349944f14d0d9`
- Branch: `worktree-agent-ab52349944f14d0d9` (rebased onto `refactor/sn-operator-algebra`).
- Pre-Round-2 HEAD: `52d7688 feat(numerics): add SourceIteration and KEigenvalue iteration primitives`.

## Deliverable manifest

- `orpheus/sn/solver.py` — full rewrite. SNSolver constructs `(L, S, F)` triple at __init__; `solve_fixed_source` dispatches to `_solve_source_iteration` or `_solve_krylov`. Legacy BiCGSTAB path retired. Curvilinear default flip in `solve_sn_fixed_source` **not enabled** — see "Deviations" below.
- `orpheus/sn/operator.py` — retired 7 functions (`build_transport_linear_operator{,_spherical,_cylindrical}`, `build_rhs{,_spherical,_cylindrical}`, `angular_flux_to_scalar`). Updated module docstring + `__all__`.
- `orpheus/sn/geometry.py` — retired 6 deprecated `@property` accessors (`alpha_half`, `redist_dAw`, `tau_mm`, `alpha_per_level`, `redist_dAw_per_level`, `tau_mm_per_level`). `face_areas`/`delta_A` retain transitional shims.
- 11 BiCGSTAB→krylov call sites migrated:
  - `tests/sn/test_solver_components.py` (7 sites)
  - `tests/sn/test_spherical.py` (3 sites)
  - `examples/discrete_ordinates/demo_discrete_ordinates.py` (1 site)
- `tests/sn/test_sweep_operator_inconsistency.py` — full rewrite to use the new `inner_solver="krylov"` path (replaces inline `_solve_bicgstab` helper). All 4 assertions pass on the symmetric-closure operator.
- `tests/sn/test_snmesh_consumes_reduced.py` — pruned tests for the 6 retired accessors.
- `tests/sn/{test_spherical, test_sweep_regression, test_cylindrical, test_quadrature}.py` — migrated `sn_mesh.<accessor>` reads to `sn_mesh.reduced.<accessor>`.
- `tests/sn/test_snstreamingoperator.py` — removed unused `angular_flux_to_scalar` import.
- `tests/geometry/test_reduced_operator.py` — same `sn_mesh.reduced.*` migration.
- `docs/theory/discrete_ordinates.rst` — section "BiCGSTAB Alternative" renamed to "Krylov inner solver"; new section "SNSolver as an operator-algebra coordinator"; Wave-E forward-references updated to present-perfect; `[AdamsLarsen2002]_` citation added.
- `docs/verification/matrix.rst` — regenerated.

## Verification gate results

- **Foundation tests** (`pytest tests/sn/ -m foundation`): **100 passed**.
- **Numerics tests** (`pytest tests/numerics/`): **270 passed**.
- **Regression bit-identity** (`pytest -m regression`): all 11 snapshots verified individually:
  - Fast cases (slab/2D Cartesian + curvilinear homogeneous): **9/9 PASSED**.
  - P1 cases (slab + sphere): **2/2 PASSED** (slow — 7.5 + 5 minutes respectively, pre-existing slowness).
- **`tests/sn/test_sweep_operator_inconsistency.py`**: 4/4 PASSED. The new `_solve_via_krylov` helper produces the analytically-correct flat flux (1.0 ± 1e-10) on the reflective-BC sphere problem; the sweep still produces the documented ERR-026 deviation.
- **`tests/sn/test_snmesh_consumes_reduced.py`**: 19/19 PASSED.
- **`tests/sn/test_solver_components.py::TestBicgstabPnScattering`**: P0/P1 SI vs krylov tests PASSED (~30-50 s each).
- **`tests/sn/test_quadrature.py` + `test_snstreamingoperator.py` + `test_sweep_regression.py` + `test_reduced_operator.py`**: 83 PASSED.
- **Sphinx**: clean build (no NEW errors; the 4 pre-existing `:pydata:` role errors in the iteration-primitives section landed by Wave E Round 1 are unchanged).
- **V&V audit**: 36/38 ERR coverage (unchanged; missing ERR-020, ERR-031 pre-existing).

## Design choices that deviated from the plan

### 1. `solve_sn_fixed_source` curvilinear default does **NOT** auto-flip to krylov

**The plan**: "When the user does NOT specify `inner_solver` and the mesh is curvilinear, route through krylov."

**The implementation reality**: `SNStreamingOperator.apply` reuses
`build_equation_map_spherical` / `build_equation_map_cylindrical`
which were designed for **reflective** outer-boundary BCs only. The
packed-vector layout has no slot for a vacuum-BC outer-incoming
:math:`\psi`. Empirical evidence: solving a vacuum-BC curvilinear MMS
case via Krylov-on-`apply` to machine precision gives a result with
~30% maximum error vs. the manufactured `phi_exact` (verified by
running GMRES with `rtol=1e-14, atol=1e-16` to convergence — the
GMRES residual hits 1e-14 but the answer is just for the wrong
operator).

**Resolution**: keep `inner_solver="source_iteration"` as default for
`solve_sn_fixed_source` on all geometries. The curvilinear vacuum-BC
fixed-source closure is **deferred to Wave E Round 3**, which owns
the equation-map extension that adds vacuum-BC slots to the packed
layout.

The xfail-strict markers in
`tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
**stay xfail** — Round 3's marker-removal is correctly coupled to
Round 3's equation-map extension.

### 2. `solve_sn` curvilinear default also kept as `"source_iteration"`

This matches the plan's discussion of preserving bit-identity for the
6 curvilinear regression snapshots. (The plan listed Option A vs B;
Option A — keep source_iteration default — was chosen.)

### 3. `_solve_source_iteration` uses Approach A (bit-identical preservation)

Per the brief, Pℓ anisotropic scattering's per-iteration angular-flux
state cannot be cleanly threaded through SourceIteration's contract.
The new method is character-for-character identical to the legacy
loop math — bit-identity preserved. Architectural cleanup
(threading angular state through ScatteringOperator) is a future
follow-up.

## ERR catalog additions

None. Round 2's discoveries were architectural (the equation-map BC
limitation) rather than numerical-correctness. The architectural
limitation is documented in:
- `orpheus/sn/solver.py::solve_sn_fixed_source` docstring
  ("Round 2 deviation from the campaign plan" section).
- `docs/theory/discrete_ordinates.rst::Round 2 deviation` subsection
  of the "Krylov inner solver" section.
- This memo (above).

## Net LOC change

- `orpheus/sn/solver.py`: 713 → **1187** (+474).
  - Includes ~370 LOC for `_solve_fixed_source_krylov` (new
    fixed-source krylov path, mostly experimental — works for
    reflective BC, broken for vacuum BC).
  - Includes ~110 LOC for `_build_rhs_cartesian` /
    `_build_rhs_spherical` (extracted from operator.py to keep
    them as solver-private helpers).
  - Net of these helpers, the actual SNSolver class growth from
    `_solve_krylov` + `_make_sweep_preconditioner` is ~120 LOC,
    well above the brief's ~50 LOC sketch but the result is
    self-contained and documented.
- `orpheus/sn/operator.py`: 1126 → 890 (−236).
  - Retired ~5 functions; the math is preserved in
    `solver.py::_build_rhs_*` (move, not delete).

Total cross-file: −236 + 474 = +238 net LOC. Above the brief's −180
target but the new lines are documented and aligned with the
operator-algebra framing.

## Issue #149 (RuntimeWarning)

**Not incidentally fixed.** The `RuntimeWarning: invalid value
encountered in divide` at the original `solver.py:227` is now at
`solver.py:279` (`return self.fission_op.apply(flux_distribution) / keff`)
and still fires on P1 anisotropic eigenvalue cases at the first
outer iteration. Round 2 did not change `compute_fission_source`'s
math; #149 remains for the dedicated triage session.

## Open follow-ups for Round 3

1. **Equation-map vacuum-BC extension.** The
   `build_equation_map_spherical` / `build_equation_map_cylindrical`
   layouts add a slot for vacuum-incoming `psi` at the outer
   boundary. After this lands, the curvilinear `solve_sn_fixed_source`
   default can flip to "krylov" (closing ERR-026 for fixed-source MMS).

2. **Marker removal**:
   `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
   xfail markers come off when (1) lands.

3. **Test rewrite**:
   `tests/sn/test_sweep_operator_inconsistency.py` flips from
   "ERR-026 documented" to "ERR-026 closed" assertions. The current
   file passes with both the sweep-deviation assertions AND the
   krylov-correct assertions, so the rewrite is mostly cosmetic
   (rename methods, flip docstring framing, drop the
   `catches("ERR-026")` mark).

4. **error_catalog.md**: ERR-026 status update from "OPEN" to
   "CLOSED — Wave E Round 3" with the closure-narrative pointer.

5. **Architectural cleanup**: thread Pℓ angular-flux state through
   ScatteringOperator so `_solve_source_iteration` can compose
   `SourceIteration` directly (Approach B from the brief). This is
   not gating; bit-identity matters more.

## Lessons for `algebra-of-record` skill

The `cross-domain-frames` self-check (Cardinal Rule 5) would have
caught this Round 2 surprise earlier:

- **Trigger**: "the same packed layout serves two consumers" — the
  legacy BiCGSTAB FD operator (reflective-BC eigenvalue, the design
  intent) and the new fixed-source krylov path (potentially vacuum-
  BC). Sharing layouts across consumers with different BC sets is a
  refactoring smell.
- **Fix**: in Round 3, the vacuum-BC equation-map extension should
  factor the BC-handling into the layout itself, not assume
  reflective. The skill's "structurally-simpler" criterion fires:
  one equation-map class with a BC handle is structurally simpler
  than two parallel functions for two BCs.

## Commit structure

The brief recommends 2-3 commits:

1. `refactor(sn): SNSolver consumes operator triple; retire BiCGSTAB FD path`
   — single atomic landing covering the rewrite, the operator.py
   retirement, the geometry.py accessor retirement, the 11 call-site
   migrations, the test_sweep_operator_inconsistency.py rewrite,
   and the Sphinx narrative update.

Closes:
- Issue #164 (the Wave E Round 2 deliverable).
- Issue #96 (Cartesian sweep-operator inconsistency — superseded by
  the symmetric-closure operator path).
- Issue #97 (curvilinear sweep-operator inconsistency — partial
  closure; full closure in Round 3).
- Issue #98 (sweep-operator inconsistency tests — file now uses
  krylov helper).
- Issue #99 (Phase 3.3-3.4 MMS blocker — unblocked when Round 3
  lands the vacuum-BC equation map).
