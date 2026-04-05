# Numerics — Improvement Tracker

Central registry of improvements for model-independent numerical methods.

## Tracking Number Format

`NM-YYYYMMDD-NNN` where NM = Numerics, YYYYMMDD = session date, NNN = sequence.

## Status Legend

- **DONE**: implemented AND documented in Sphinx
- **IMPL**: implemented and tested, Sphinx documentation pending
- **OPEN**: not yet implemented, documented here with full context

---

## DONE — Implemented and Documented

### NM-20260405-001 — EigenvalueSolver Protocol and power_iteration()

**Status**: DONE  
**Commit**: 3d91d18, 0a69603  
**Sphinx**: docs/api/numerics.rst (autodoc)

Created `numerics/eigenvalue.py` with:
- `EigenvalueSolver` Protocol defining the contract for any deterministic solver
  (initial_flux_distribution, compute_fission_source, solve_fixed_source,
  compute_keff, converged)
- Generic `power_iteration()` function that works with any conforming solver
- Separation: power iteration returns flux distribution (fundamental mode),
  absolute normalization is post-processing

All 6 deterministic solvers refactored to satisfy the Protocol:
HomogeneousSolver, CPSolver, SN1DSolver, DO2DSolver, MoCSolver, DiffusionSolver.

---

## OPEN — Not Yet Implemented

### NM-20260405-002 — Full eigenvalue spectrum solver

**Priority**: Medium | **Effort**: Large

Currently only power_iteration exists, which converges to the dominant
eigenvalue (fundamental mode).  A full spectrum solver would compute all
eigenvalues and eigenvectors using e.g. Arnoldi iteration or subspace
iteration.  This enables:
- Higher harmonics analysis (lambda modes)
- Dominance ratio computation
- Stability analysis

### NM-20260405-003 — Convergence acceleration (Chebyshev or Wielandt)

**Priority**: Medium | **Effort**: Moderate

Power iteration convergence rate is |k₁/k₀| (dominance ratio).  For
near-critical systems this can be close to 1, making convergence slow.
Chebyshev acceleration or Wielandt spectral shift can significantly
reduce iteration count.

### NM-20260405-004 — Sphinx theory chapter for eigenvalue methods

**Priority**: HIGH | **Effort**: Moderate

Missing `docs/theory/eigenvalue_methods.rst`.  Should cover:
- Power iteration derivation and convergence analysis
- Fundamental mode uniqueness (Perron-Frobenius)
- Dominance ratio and convergence rate
- Chebyshev/Wielandt acceleration (when implemented)
- Full spectrum methods (when implemented)
- Connection to the EigenvalueSolver Protocol
