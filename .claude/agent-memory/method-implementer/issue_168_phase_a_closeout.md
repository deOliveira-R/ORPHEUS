---
name: Issue #168 Phase A — BoundaryFaceFlux Protocol + DD extrapolation
description: Phase A of Issue #168 Option C ships on `refactor/sn-operator-algebra` 2026-05-10. Addresses Defects 1+2 (outer-face cell-center truncation; cell-center / BC-face-value storage conflation) via new BoundaryFaceFlux Protocol parallel to CellUpdate. Defect 3 (sphere-pole stencil) preserved AS-IS pending Phase B.
type: project
---

# Issue #168 Phase A — closeout

`refactor/sn-operator-algebra` 2026-05-10. Addresses **Defects 1 + 2**
of the Issue #168 three-defect curvilinear FD operator triage. Defect
3 (sphere-pole redistribution-term stencil) is preserved AS-IS pending
Phase B literature consultation (Hébert §3.9.4 per literature-researcher
memo `sphere_sn_pole_closure_canonical.md`).

## Why this memo exists

Future implementers (and this same implementer on a fresh-context return)
need the Phase A surface area pinned: what shipped, what remains open,
and which architectural deviations from the design memo were taken.

**Why it matters:** ERR-026 stays at PARTIAL CLOSURE through Phase A.
The 4 xfail-strict markers on curvilinear MMS tests stay xfail (orders
~1.5-1.7, below the > 1.9 bar). Phase B is the closing fix.

**How to apply:** before scoping Phase B, reread this memo + the
literature-researcher memo on Hébert §3.9.4. The pole-stencil
correction ships on top of the Phase A architecture — the
`BoundaryFaceFlux` Protocol's `cell_idx = 0` branch is the natural
extension point.

## What shipped

### New module: `orpheus/sn/spatial/boundary_face_flux.py`

Parallel to `cell_update.py` (Wave C). Defines:

- `BoundaryFaceFlux` — `@runtime_checkable Protocol` with reduced
  signature `(psi_cells, cell_idx) -> ndarray`.
- `BoundaryFaceFluxBase` — concrete ABC layered on `RegistryMixin`;
  subclasses self-register under `key=` class-creation kwarg.
- `DDExtrapolation` (default, registered as `"dd_extrapolation"`) —
  one-sided second-order DD diamond extrapolation
  `psi_face_out = 1.5·psi[N-1] - 0.5·psi[N-2]` at outer face;
  pole closure at `cell_idx=0` returns `psi[0]` unchanged
  (Phase-A preservation; A_0 = 0 in spherical so the value is
  multiplied by zero, the redistribution Defect-3 fix is separate).
- `CellCenter` (registered as `"cell_center"`) — legacy reproducer,
  `psi_face_out = psi[cell_idx]`. Kept for ablation studies and
  back-bisection.

### Decoupled cell-center / BC-face-value storage

`solution_to_angular_flux_spherical` (and the `_cylindrical` alias) now
returns **`(fi, boundary_face_flux)`** instead of just `fi`:

- `fi` shape `(ng, N, nx, 1)` — pure cell-center storage. For inward
  μ < 0 ordinates at i=N-1, the slot holds the **reflected-partner**
  cell-center value (faithful, BC-independent), NOT the BC face value
  (which was the Defect-2 conflation).
- `boundary_face_flux` shape `(ng, N, 2)` — per-(group, ordinate)
  BC-resolved face fluxes at index 0 (inner pole) and 1 (outer face).
  The matvec reads index 1 for inward μ < 0 ordinates instead of
  reading from `fi[..., -1, 0]`.

The 4 callers in `orpheus/sn/solver.py` are updated to `fi, _ =
solution_to_angular_flux_spherical(...)` (the boundary-face-flux
output is unused at the scalar-flux extraction site — only the
matvec consumes it via `transport_operator_matvec_spherical`).

### Matvec signature extensions

`transport_operator_matvec_spherical` and `_cylindrical` accept a
new optional kwarg `boundary_face_flux_closure: BoundaryFaceFlux |
None = None` (defaults to `DDExtrapolation()`).

### SNStreamingOperator + SNMesh wiring

- `SNMesh.__init__` accepts `boundary_face_flux: BoundaryFaceFlux |
  None = None` (parallel to `cell_update`). Default is
  `DDExtrapolation()`.
- `SNStreamingOperator.apply` threads `sn_mesh.boundary_face_flux`
  through to the spherical / cylindrical matvec.
- Cartesian path is unaffected — the upwind FD stencil there has no
  symmetric closure to break — and ignores `sn_mesh.boundary_face_flux`.

### Tests

- **NEW** `tests/sn/spatial/test_boundary_face_flux.py` — 21
  foundation-tagged tests covering Protocol conformance, DD-extrapolation
  hand-calc, pole closure, CellCenter reproducer, registry
  self-registration, immutability.
- **EXTENDED** `tests/sn/test_snstreamingoperator.py` — 7 new tests
  pinning the Phase-A invariants:
  - Tuple return signature for the curvilinear decoder.
  - Defect 2 fix: `fi[..., -1, 0]` for inward ordinates is the
    reflected-partner cell-center, NOT the BC face value (verified
    explicitly under VacuumBoundaryOperator).
  - Defect 1 fix: `apply` matches `DDExtrapolation` and is NOT
    `CellCenter` on a non-constant input.
  - Per-ordinate flat-flux invariant (constant ψ + zero Σ_t →
    zero residual everywhere) for both spherical and cylindrical.
  - Defect 2 reproducer: vacuum BC + constant cell flux gives zero
    streaming residual at ALL cells (no contamination at i=N-2).
  - SNMesh constructor accepts and threads the new kwarg.

### Sphinx narrative

`docs/theory/discrete_ordinates.rst` gains a new subsection
"Boundary face-flux strategies (Issue #168 Phase A)" under the
Krylov inner solver section. Documents the three defects, what
Phase A addresses (1+2), and what remains for Phase B (3).

### Snapshot deletions

The 6 curvilinear regression snapshots are intentionally invalidated
by Phase A and **DELETED**:

- `sphere_2g_homogeneous_dd_n20.npz`
- `sphere_2g_3reg_dd_n40.npz`
- `sphere_2g_p1_aniso_dd_n20.npz`
- `cyl_1g_homogeneous_LS4_dd_n20.npz`
- `cyl_1g_homogeneous_product_dd_n20.npz`
- `cyl_2g_3reg_LS4_dd_n40.npz`

They skip gracefully via the existing `if not snapshot_file.exists():
pytest.skip(...)` mechanism at `test_dd_regression.py:47`. Phase B
follow-up will regenerate them with FN-method cross-check.

The 5 Cartesian snapshots stay bit-identical (verified: `pytest
-m regression` reports `5 passed, 6 skipped`).

## Architectural deviation: reduced Protocol signature

The design memo (`.claude/plans/issue_168_design.md` §6.2) sketched
a five-argument signature:

```python
def __call__(self, psi_cells, cell_idx, side, ord_idx, bc) -> ndarray
```

I shipped a reduced two-argument signature:

```python
def __call__(self, psi_cells, cell_idx) -> ndarray
```

**Rationale:**

- `side` (`"inner"` / `"outer"`) is derivable from `cell_idx`
  against a known `nx` (cell_idx=0 is inner, cell_idx=nx-1 is outer).
- `ord_idx` was purely diagnostic / future-angular-aware; never read
  by the closure math.
- `bc` is only consulted for the inner pole, which Phase A doesn't
  touch (it's preserved as-is pending Phase B). When Phase B needs
  BC awareness for the pole correction, the BC will be resolved at
  the matvec call site (where `bc_outer` already lives) and either
  baked into the BoundaryFaceFlux constructor as state OR threaded
  through a Phase B-specific second method on the Protocol.

The minimal signature gives the strategy exactly the data the
DD-extrapolation algebra needs and nothing more — better contract,
cleaner test surface, fewer mock parameters in the test fixtures.

This deviation preserves §6 mathematical correctness fully — the
DD-extrapolation formula
`psi_face_out = 1.5·psi[N-1] - 0.5·psi[N-2]` from §6.4 is the
exact expression `DDExtrapolation.__call__` returns. The pole
closure `psi_cells[:, 0]` matches §6.1's "use the standard
treatment (Lewis & Miller §4.5): initialize ψ^face_n,1/2 = ψ_n,0".

## Verification gates — all green

- `pytest tests/sn/spatial/test_boundary_face_flux.py -q` → 21 passed.
- `pytest tests/sn/test_snstreamingoperator.py -q` → 29 passed
  (22 original + 7 new Phase-A invariants).
- `pytest -m regression` → 5 passed (Cartesian, bit-identical),
  6 skipped (curvilinear, intentionally invalidated).
- `pytest tests/sn/l1_analytical/ -m "not slow"` → 15 passed,
  2 xfailed (the curvilinear MMS — orders ~1.5-1.7, still below
  the > 1.9 xfail-strict bar).
- `pytest tests/sn/test_mms_curvilinear.py
  tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
  → 4 xfailed (all ERR-026 tripwires stay xfail correctly).
- `sphinx-build -W docs docs/_build/html` → exit 0.
- `python -m tests._harness.audit` → 23 orphan equations + 36/38
  ERR coverage, identical to baseline.

## What remains for Phase B

1. **Sphere-pole stencil correction (Defect 3).** Hébert §3.9.4
   (3.418)-(3.439) is the canonical reference per
   literature-researcher memo `sphere_sn_pole_closure_canonical.md`.
   The Bailey 2009 docstring in `reduced_operator.py` cites
   the WRONG paper — should cite Bailey-Morel-Chang 2010 NSE 165
   (SN angular), NOT Adams-Yang-Zika 2008 JCP (diffusion).
2. **Snapshot regeneration.** After Phase B, regenerate the 6
   curvilinear regression snapshots from the corrected operator
   AND verify them against FN-method (sood_registry) cross-check
   before committing.
3. **xfail-strict markers come off.** Currently 4 xfail tripwires
   on curvilinear MMS; Phase B should give orders > 1.9.

## Files touched

- NEW: `orpheus/sn/spatial/boundary_face_flux.py` (415 lines).
- MODIFIED: `orpheus/sn/operator.py` (Defect 1+2 fixes; new tuple
  return; matvec kwargs).
- MODIFIED: `orpheus/sn/spatial/__init__.py` (re-exports).
- MODIFIED: `orpheus/sn/geometry.py` (`SNMesh.boundary_face_flux`
  attribute).
- MODIFIED: `orpheus/sn/solver.py` (4 caller-site tuple unpacks).
- NEW: `tests/sn/spatial/test_boundary_face_flux.py` (232 lines).
- MODIFIED: `tests/sn/test_snstreamingoperator.py` (7 new tests).
- MODIFIED: `docs/theory/discrete_ordinates.rst` (Phase-A subsection).
- DELETED: 6 curvilinear regression snapshots.
- DELETED: zero — clean cut, no backward-compat shims (per brief).

## Related work

- `.claude/plans/issue_168_design.md` — the design memo.
- `.claude/agent-memory/numerics-investigator/issue_168_three_defects.md`
  — empirical falsification of the single-bug framing.
- `.claude/agent-memory/literature-researcher/sphere_sn_pole_closure_canonical.md`
  — Hébert §3.9.4 as the Phase-B canonical reference.
- ERR-026 in `.claude/skills/vv-principles/error_catalog.md`.
- `.claude/agent-memory/method-implementer/sn_reshape_wave_e_round_3_closeout.md`
  — Wave E Round 3 (the BC plumbing predecessor).
