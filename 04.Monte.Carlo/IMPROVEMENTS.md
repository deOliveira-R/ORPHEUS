# Monte Carlo — Improvement Tracker

Central registry of ALL bugs, improvements, and features for the
Monte Carlo solver.

## Tracking Number Format

`MT-YYYYMMDD-NNN` where MT = Monte Carlo Transport.

## Status Legend

- **DONE**: implemented AND documented in Sphinx
- **IMPL**: implemented and tested, Sphinx documentation pending
- **OPEN**: not yet implemented, documented here with full context

---

## IMPL — Implemented, Sphinx Documentation Pending

### MT-20260403-001 — MCGeometry protocol for delta-tracking

`MCGeometry` protocol: `material_id_at(x, y) -> int` + `pitch`.
Designed for delta-tracking (no distance-to-surface needed).
Extensible to future CSG implementations.

Two concrete implementations:
- `ConcentricPinCell(radii, mat_ids, pitch)` — annular regions
- `SlabPinCell(boundaries, mat_ids, pitch)` — 1D slab regions

### MT-20260403-002 — Homogeneous verification {1,2,4}G

Analytical derivation from random walk probability theory:
k = νΣf/Σa (1G), k = λ_max(A⁻¹F) (multi-group).
Statistical verification via z-score < 5σ.
1G homogeneous gives deterministic σ=0 (all neutrons see same XS).

---

## OPEN — Not Yet Implemented

### MT-20260403-003 — Sphinx theory chapter

**Priority**: HIGH | **Effort**: Large

No `docs/theory/monte_carlo.rst` exists.  Needs to document:
- Woodcock delta-tracking algorithm
- Analog absorption with fission weight adjustment
- Russian roulette and splitting
- Fission spectrum sampling
- keff estimator and statistical uncertainty (CLT)
- MCGeometry protocol design (CSG-extensible)
- Verification results (homogeneous exact + heterogeneous)

### MT-20260403-004 — Python neutron loop performance

**Priority**: HIGH | **Effort**: Large

The inner neutron random walk is a Python `while True` loop with
per-collision Python-level operations.  For 421-group problems with
many collisions per neutron, this dominates runtime.  Options:
- Batch neutrons: process all neutrons for one free-path step together
- Vectorize collision sampling across the neutron population
- Consider Cython/numba for the inner loop

### MT-20260403-005 — Heterogeneous independent reference

MC heterogeneous verification uses CP cylinder eigenvalue as proxy.
Should use high-statistics MC run (10⁵ active cycles) as independent
reference.  See DV-20260403-005.

### MT-20260403-006 — Majorant cross section per group

Current majorant is `max(SigT)` across all materials for each group.
For problems with highly variable cross sections, a more efficient
majorant (e.g., piecewise by region) would reduce virtual collision
rate.
