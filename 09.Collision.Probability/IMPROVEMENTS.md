# Collision Probability Solver — Improvement Tracker

Central registry of ALL bugs, improvements, and features for the CP
solver (`CPMesh`, `CPSolver`, `solve_cp`).

## Tracking Number Format

`CP-YYYYMMDD-NNN` where CP = Collision Probability module.

## Status Legend

- **DONE**: implemented AND documented in Sphinx
- **IMPL**: implemented and tested, Sphinx documentation pending
- **OPEN**: not yet implemented, documented here with full context

---

## DONE — Implemented and Documented

### CP-20260404-001 — CPMesh augmented geometry

Single `CPMesh` class wraps Mesh1D and selects kernel by coordinate
system: E₃ (slab), Ki₄ (cylindrical), exp (spherical).  Replaces
separate `SlabGeometry` and `CPGeometry` classes.

Documented in `docs/theory/collision_probability.rst` §Architecture.

### CP-20260404-002 — Unified solve_cp()

Single `solve_cp()` replaces `solve_cp_slab()` + `solve_cp_concentric()`.
Kernel selection is automatic from `mesh.coord`.

### CP-20260404-003 — Spherical CP kernel and verification

Exponential kernel with y-weighted quadrature.  9 verification cases
(1/2/4 groups × 1/2/4 regions) in `derivations/cp_sphere.py`.
All pass to < 1e-5.

Documented in `docs/theory/collision_probability.rst` §Spherical.

### CP-20260404-004 — Unified white-BC closure

`_apply_white_bc()` works for all geometries via `mesh.volumes` and
`mesh.surfaces[-1]`.

---

## IMPL — Implemented, Sphinx Documentation Pending

### CP-20260404-005 — Second-difference formula documentation

The general second-difference formula `Δ₂[F](τ_i, τ_j, gap)` is
documented in the theory chapter but lacks a dedicated SymPy derivation
script showing the algebraic equivalence between slab/cylindrical/spherical
forms.

---

## OPEN — Not Yet Implemented

### CP-20260404-006 — Ray-tracing CP for arbitrary 2D geometry

**Priority**: Low | **Effort**: Large

For non-circular cell shapes (real lattice boundaries), the CP matrix
must be computed by tracing rays at multiple angles.  This bridges
CP and MOC.  Not needed for 1D geometries.

**References**: APOLLO, DRAGON codes (Hébert 2009).

### CP-20260404-007 — Interface current method for multi-cell lattice

**Priority**: Medium | **Effort**: Large

Current implementation uses white-BC closure for a single cell.
The interface current method couples individual cell CP solutions at
boundaries for multi-cell lattice calculations, giving more accurate
boundary treatment than white BC.

**References**: Hébert (2009) §8.

### CP-20260404-008 — Performance: vectorize the group loop

**Priority**: Low | **Effort**: Small

`compute_pinf_group()` is called in a Python loop over groups.
The Ki₄ table lookup and second-difference could be vectorized over
groups for a significant speedup on 421-group problems.
