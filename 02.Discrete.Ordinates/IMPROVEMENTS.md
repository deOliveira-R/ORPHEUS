# SN Discrete Ordinates — Improvement Tracker

Central registry of ALL bugs, improvements, and features for the SN
solver.  **This is the single source of truth.**

## Tracking Number Format

`DO-YYYYMMDD-NNN` where:
- `DO` = Discrete Ordinates module
- `YYYYMMDD` = session date when the item was created
- `NNN` = sequential number within the session
- `00000000` = item predates the tracking system

## Status Legend

- **DONE**: implemented AND documented in Sphinx
- **IMPL**: implemented and tested, Sphinx documentation pending
- **OPEN**: not yet implemented, documented here with full context
- **WONT**: decided against, with rationale

## Where TODOs Live

TODOs exist in exactly TWO places:
1. **This file** — every item with its tracking number, status, and context
2. **In code** — at the exact location where the fix goes, with the
   matching tracking number (e.g., `# TODO DO-20260405-001: ...`)

No other files should contain TODOs.  If you find one, it must be
consolidated here and given a tracking number.

---

## DONE — Implemented and Documented in Sphinx

### DO-20260404-001 — Geometry-weighted balance equation (Bailey et al. 2009)

Documented in ``docs/theory/discrete_ordinates.rst`` §5 (The Discrete Balance Equation).

### DO-20260404-002 — Morel–Montry angular closure weights

Documented in ``docs/theory/discrete_ordinates.rst`` §5.5–5.6.

### DO-20260404-003 — BiCGSTAB operators for curvilinear geometries

Documented in ``docs/theory/discrete_ordinates.rst`` §7.

### DO-20260404-004 — Contamination analysis tool

Documented in ``docs/theory/discrete_ordinates.rst`` §5.4.

### DO-20260405-001 — Consolidated ΔA/w into SNMesh

Documented in ``docs/theory/discrete_ordinates.rst`` §2.

### DO-20260405-004 — Sphinx theory chapter for curvilinear SN

Self-referential: this IS the Sphinx chapter (1276 lines, zero warnings).

### DO-20260403-001 — Volume halving bug in CartesianMesh.volume

Fixed: removed boundary halving from ``CartesianMesh.volume``.  Reflective
BCs mean boundary cells are at the symmetry plane (full volume, not half).
Bug caused ~1e-4 systematic keff error on heterogeneous problems.
Documented in ``tests/l0_error_catalog.md`` ERR-008.

---

## OPEN — Not Yet Implemented

### DO-20260405-002 — Gauss-type azimuthal quadrature for cylindrical

**Priority**: Medium | **Effort**: Moderate  
**Code location**: `sn_quadrature.py` (new quadrature class)

The equally-spaced `ProductQuadrature` gives duplicate η values
(paired ±ξ ordinates), producing alternating M-M weights
τ = [0.5, 1.0, ...].  A Gauss-type azimuthal quadrature with
non-uniform φ spacing would give distinct η values and smoothly
varying τ.

**References**: standard quadrature construction techniques.

### DO-20260405-003 — φ-based cell-edge computation for non-product quadratures

**Priority**: Low | **Effort**: Small  
**Code location**: `sn_geometry.py:275`

For quadratures where η values are distinct, transforming actual
φ cell boundaries to η-space would give exact M-M cell edges
instead of the midpoint approximation.

**References**: Bailey et al. (2009) Eq. 52.

### DO-20260403-002 — 2D Lebedev 421-group performance

**Priority**: HIGH | **Effort**: Moderate
**Code location**: `sn_sweep.py` (wavefront path)

The 2D wavefront sweep with 110 Lebedev ordinates × 10×10 mesh × 421
groups runs at ~0.3s per outer iteration (source iteration).  With
~200 outer iterations needed for convergence, the full demo takes ~60s.
The old BiCGSTAB solver was even slower.

The main performance bottleneck is the Python `for n in range(N)` loop
over 110 ordinates in `_sweep_2d_wavefront`.  Ordinate batching by
octant sign pattern could reduce this to 4 sweeps.

Related: DSA (DO-00000000-001) would reduce outer iterations from ~200
to ~20.

### DO-00000000-001 — Diffusion Synthetic Acceleration (DSA)

**Priority**: HIGH | **Effort**: Large  
**Code location**: `sn_solver.py` (new acceleration method)

Reduces outer iterations from ~200 to ~20 for many-group problems.
Uses diffusion correction after each source iteration.
1D diffusion solver already exists in `05.Diffusion.1D/`.

**References**: Adams & Larsen (2002), Wareing et al.

### DO-00000000-002 — Transport Synthetic Acceleration (TSA)

**Priority**: Low | **Effort**: Large  
**Code location**: `sn_solver.py`

Coarse-angle transport acceleration for highly anisotropic media.
Only needed if DSA proves insufficient.

**References**: Ramone et al. (1997).

### DO-00000000-003 — Linear Discontinuous (LD) angular finite elements

**Priority**: Low | **Effort**: Large  
**Code location**: `sn_sweep.py` (new sweep variant)

Second-order angular accuracy without flux dip.  2×2 system per
cell-ordinate.  M-M WDD (DO-20260404-002) already eliminates the
flux dip, so LD is only needed for higher angular accuracy.

**References**: Bailey, Morel & Chang (2009) — main topic of paper.

### DO-00000000-004 — Negative flux fixup

**Priority**: Low | **Effort**: Small  
**Code location**: `sn_sweep.py` (inner loop guard)

If WDD produces ψ^a_out < 0, clamp to zero and rebalance.  Currently
not needed (zero negatives observed), but good practice for extreme
cases.

### DO-00000000-005 — Transport eigenmodes (Case's method)

**Priority**: Medium | **Effort**: Large  
**Code location**: `derivations/` (new module)

Mesh-independent analytical/semi-analytical reference for 1D
multi-group transport eigenvalues.  Currently the only independent
reference is the diffusion transfer matrix (`sn_heterogeneous.py`),
which has a ~0.3% transport correction.

**The problem**: SN heterogeneous verification uses Richardson
extrapolation from the SN solver itself — a self-referencing test.
Case's method would provide a truly independent reference.

**Case's eigenmode method**: The exact 1D monoenergetic transport
solution with isotropic scattering decomposes into:

1. **Discrete modes** ν₀ from the dispersion relation:
   `1 = c·ν₀·Σ_t·arctanh(1/(ν₀·Σ_t))` where c = Σ_s/Σ_t.
2. **Continuum modes** ν ∈ [-1/Σ_t, 1/Σ_t]: singular eigenfunctions
   with Cauchy principal value integrals.
3. **Multi-group**: matrix dispersion relation
   `det[I − Σ_s^T · diag(1/Σ_t) · Λ(ν)] = 0`.

**Interface matching**: half-range flux continuity
∫₀¹ μⁿ ψ_left dμ = ∫₀¹ μⁿ ψ_right dμ, truncated at order M.

**Practical alternative**: The F_N method (Siewert, Garcia) avoids
explicit continuum modes by expanding in Chandrasekhar polynomials.
N=20 gives ~10 digits — sufficient for our verification.

**Implementation plan**:
- 1G 1-region: straightforward (dispersion + BCs)
- 1G multi-region: moderate (interface matching)
- Multi-group multi-region: significant (matrix dispersion + matching)

**Existing diffusion reference** (`sn_heterogeneous.py`):

| Case | Diffusion | SN Richardson | Diff |
|------|-----------|---------------|------|
| 1G 2-region | 1.2646 | 1.2605 | +0.0041 |
| 2G 2-region | 1.2338 | 1.2380 | -0.0042 |
| 4G 2-region | 1.0312 | 1.0344 | -0.0032 |

Diffusion can be above or below transport depending on configuration.

**References**: Case (1960), Case & Zweifel (1967), Siewert (2000),
Garcia & Siewert (various).

### DO-00000000-006 — Anisotropic scattering in curvilinear sweeps

**Priority**: Medium | **Effort**: Moderate  
**Code location**: `sn_sweep.py` (spherical/cylindrical branches)

P1+ anisotropic scattering implemented for Cartesian 2D but NOT
verified for curvilinear 1D.  Spherical harmonics on GL/Product
quadrature needs verification.

### DO-20260405-005 — Document quadrature–dimension mismatch behaviour

**Priority**: Medium | **Effort**: Small  
**Code location**: `sn_quadrature.py`, `docs/theory/discrete_ordinates.rst`

In the 2026-04-05 session, a ~0.3% eigenvalue plateau was observed when
using the 3D Lebedev quadrature (110 points on the unit sphere) for 1D
slab problems.  The error was controlled by angular quadrature, not mesh
refinement.  Switching to Gauss–Legendre (correct 1D quadrature)
eliminated the bias entirely.

**Status unclear:** the codebase has since been restructured with a
quadrature dispatch system.  Need to verify:
1. Does the current dispatcher avoid this mismatch automatically?
2. If Lebedev is manually applied to 1D, does the bias still appear?
3. Document findings in the Quadrature Comparison section of the theory chapter.

### DO-00000000-007 — GMRES/preconditioned Krylov for BiCGSTAB

**Priority**: Low | **Effort**: Moderate  
**Code location**: `sn_solver.py` (`_solve_bicgstab` methods)

BiCGSTAB can stagnate on non-normal operators.  GMRES(m) or
sweep-preconditioned BiCGSTAB may be more robust.

**References**: standard iterative methods literature.
