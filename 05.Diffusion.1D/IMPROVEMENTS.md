# 1D Diffusion — Improvement Tracker

Central registry of ALL bugs, improvements, and features for the
1D two-group neutron diffusion solver.

## Tracking Number Format

`DF-YYYYMMDD-NNN` where DF = Diffusion.

## Status Legend

- **DONE**: implemented AND documented in Sphinx
- **IMPL**: implemented and tested, Sphinx documentation pending
- **OPEN**: not yet implemented, documented here with full context

---

## IMPL — Implemented, Sphinx Documentation Pending

### DF-20260403-001 — 2-region Richardson-extrapolated reference

Fuel + reflector slab with vacuum BCs.  Reference eigenvalue from
O(h²) Richardson extrapolation at 4 mesh levels (dz = 2.5, 1.25,
0.625, 0.3125 cm).  k_eff ≈ 0.870 (reflector savings vs bare k=0.821).

---

## OPEN — Not Yet Implemented

### DF-20260403-002 — 4-region support (extend solver to >2 materials)

**Priority**: Medium | **Effort**: Moderate
**Code location**: `diffusion_1d.py`

Solver currently takes exactly 2 `TwoGroupXS` objects (fuel + reflector).
Should accept a per-cell material map (like SN and CP solvers) to
support 4-region verification cases (fuel + gap + clad + reflector).

### DF-20260403-003 — Proper 2-group interface matching (analytical)

**Priority**: Medium | **Effort**: Moderate
**Code location**: `derivations/diffusion.py`

Replace Richardson reference for 2-region with true analytical solution.
Requires solving coupled 2-group interface matching:
- Fuel: 2-group eigenmodes (2 pairs of cos/sinh modes)
- Reflector: 2-group decay modes (2 exponential decay pairs)
- Interface: 4 conditions (flux + current × 2 groups)
- Result: 4×4 determinant → transcendental equation

A 1-group fast-group-only approximation was attempted and gave 12%
error (k=0.978 vs k=0.870), confirming that 2-group coupling is
essential.  See DV-20260403-004.

### DF-20260403-004 — Sphinx theory chapter

**Priority**: HIGH | **Effort**: Large

No `docs/theory/diffusion.rst` exists.  Needs to document:
- Diffusion approximation from transport equation
- Two-group equations (fast + thermal, downscatter only)
- Finite-difference discretization
- Vacuum boundary conditions
- Bare slab buckling eigenvalue (analytical, with SymPy derivation)
- Fuel + reflector interface matching (when implemented)
- O(h²) spatial convergence proof and verification
