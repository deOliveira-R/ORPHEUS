# Geometry Module — Improvement Tracker

Central registry of ALL bugs, improvements, and features for the
geometry package (`Mesh1D`, `Mesh2D`, factories, coordinate formulas).

## Tracking Number Format

`GE-YYYYMMDD-NNN` where GE = Geometry module.

## Status Legend

- **DONE**: implemented AND documented in Sphinx
- **IMPL**: implemented and tested, Sphinx documentation pending
- **OPEN**: not yet implemented, documented here with full context

---

## DONE — Implemented and Documented

### GE-20260404-001 — Base geometry package

Mesh1D, Mesh2D, CoordSystem enum, volume/surface formulas for
Cartesian/Cylindrical/Spherical.  Frozen dataclasses with computed
properties.  64 unit tests.

Documented in `docs/theory/collision_probability.rst` §Architecture
and `docs/theory/discrete_ordinates.rst` §Architecture.

### GE-20260404-002 — Zone-based mesh construction

`mesh1d_from_zones()` with equal-volume subdivision per coordinate
system (Cartesian: equal-width, Cylindrical: sqrt, Spherical: cbrt).
Verified by equal-volume property tests.

### GE-20260404-003 — PWR convenience factories

`pwr_slab_half_cell`, `pwr_pin_equivalent`, `pwr_pin_2d`,
`homogeneous_1d`, `slab_fuel_moderator`.  All match legacy geometry
class outputs.

---

## OPEN — Not Yet Implemented

### GE-20260404-004 — Mesh2D for cylindrical R-Z geometry

**Priority**: Low | **Effort**: Moderate

`Mesh2D(coord=CoordSystem.CYLINDRICAL)` computes volumes as
`π·Δ(r²)·Δz` but no solver currently uses 2D cylindrical.
Should be verified when a 2D R-Z solver is implemented.

### GE-20260404-005 — Non-uniform zone subdivision strategies

**Priority**: Low | **Effort**: Small

Current factories only support equal-volume subdivision within each
zone.  Geometric or logarithmic spacing would be useful for boundary
layers (thin cells near material interfaces).  The `Zone` dataclass
could accept an optional `spacing` parameter.

### GE-20260405-001 — 3D mesh terminology convention

**Status**: DONE  
**Commit**: 48925ea

Adopted consistent 3D finite-volume mesh terminology across all modules,
even in reduced dimensions (1D/2D are spatially degenerate 3D cases):
- **cell / volume** — where material properties and volume-averaged fields live
- **face** — interface between cells, where D, J, ∇φ are defined

Renames applied: `n_nodes→n_cells`, `n_edges→n_faces`, `z_nodes→z_cells`,
`z_edges→z_faces`, `chi_node→chi_cell`, `sig_s_node→sig_s_cell`,
`sig2_node→sig2_cell`, `sig_t_edges→sig_t_face`.

### GE-20260404-006 — Mesh refinement utility

**Priority**: Medium | **Effort**: Small

A `refine(mesh, factor)` function that doubles/triples the number of
cells while preserving material assignment.  Useful for convergence
studies without rebuilding from zones each time.
