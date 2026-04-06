# Method of Characteristics — Improvement Tracker

Central registry of ALL bugs, improvements, and features for the
MOC solver.

## Tracking Number Format

`MC-YYYYMMDD-NNN` where MC = Method of Characteristics.
(Note: Monte Carlo uses `MT` to avoid collision.)

## Status Legend

- **DONE**: implemented AND documented in Sphinx
- **IMPL**: implemented and tested, Sphinx documentation pending
- **OPEN**: not yet implemented, documented here with full context

---

## IMPL — Implemented, Sphinx Documentation Pending

### MC-20260403-001 — MoCGeometry.from_annular() factory

Builds a Cartesian material map from concentric annular radii.
Cell material is determined by distance from pin center to each
cell center.  Enables arbitrary-region pin cell calculations.

### MC-20260403-002 — Generic flux_per_material result

`MoCResult.flux_per_material: dict[int, np.ndarray]` replaces
hardcoded `flux_fuel/clad/cool`.  Legacy properties preserved via
`@property` aliases for backward compatibility.

### MC-20260403-003 — Homogeneous verification {1,2,4}G

Analytical derivation from characteristic ODE: for homogeneous
medium, flat-source solution is exact along every ray → k = νΣf/Σa.
Verified to < 1e-4 tolerance.

---

## OPEN — Not Yet Implemented

### MC-20260403-004 — Sphinx theory chapter

**Priority**: HIGH | **Effort**: Large

No `docs/theory/moc.rst` exists.  Needs to document:
- Characteristic ODE and its analytical solution
- Flat-source approximation and attenuation formula
- 8-direction ray tracing on Cartesian mesh
- Reflective boundary conditions via ray reflection
- Power iteration with the EigenvalueSolver protocol
- Verification results (homogeneous exact + heterogeneous Richardson)

### MC-20260403-005 — Python loop bottleneck in ray sweep

**Priority**: HIGH | **Effort**: Large

The `_fly_from()` function has Python loops over cells and directions.
For 421-group problems this dominates runtime.  Options:
- Vectorize ray segments within each direction
- Batch rays by direction using numpy array operations
- Consider Cython/numba for the inner loop

### MC-20260403-006 — Heterogeneous Richardson extrapolation

Richardson references for {1,2,4}eg × {2,4}rg computed from 3 mesh
levels (8, 12, 16 cells/side).  Should be 4 levels with ratio-2
for better accuracy.  See DV-20260403-006.
