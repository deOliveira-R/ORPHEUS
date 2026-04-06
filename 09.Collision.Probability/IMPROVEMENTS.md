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

### CP-20260405-001 — Neutron balance transpose fix (ERR-009)

**Status**: DONE
**Commit**: fa5732c
**Sphinx**: docs/theory/collision_probability.rst §Flat-Source Approximation
**Error catalog**: tests/l0_error_catalog.md ERR-009

Fixed `P @ source` → `P.T @ source` in the power iteration.  With
`P[i,j] = P(birth_i → collision_j)`, the collision rate in target `j`
sums over birth regions: `Σ_i P[i,j] · V_i · Q_i`.  This is a column
operation on P, requiring the transpose.

Homogeneous problems were unaffected (P symmetric).  Detected by the
formal verification suite on the 1G 2-region slab benchmark (k=1.373
vs analytical 1.272, 8% error).

### CP-20260405-002 — NG refactoring for arbitrary group count

**Status**: DONE
**Commit**: fa5732c

All solvers now infer group count from `Mixture.ng` instead of
importing the global constant `NG = 421`.  Enables synthetic
1/2/4-group verification benchmarks.  See DA-20260405-003.

### CP-20260405-003 — Gauss-Seidel solver mode and diagnostics

**Status**: DONE
**Sphinx**: docs/theory/collision_probability.rst §Inner vs Outer Iterations

Added `solver_mode="gauss_seidel"` to `CPParams` and `CPSolver`.
Sweeps groups fast-to-thermal with inner iterations for within-group
scattering convergence.  Diagnostics added to `CPResult`:

- `residual_history`: neutron balance residual per outer iteration
  (both Jacobi and GS modes).
- `n_inner`: inner iteration counts per group per outer iteration
  (GS mode only, shape `(n_outer, ng)`).

Verified: all 27 CP eigenvalue cases (slab/cylinder/sphere ×
1/2/4 groups × 1/2/4 regions) produce identical eigenvalues in
both Jacobi and GS modes.  36 new tests in `tests/test_cp_diagnostics.py`.

Plotting: `plot_cp_convergence` shows residual subplot;
`plot_cp_inner_iterations` shows heatmap of inner iterations per group.

### CP-20260405-004 — Fix tautological GS inner residual (C-1)

**Status**: DONE
**Sphinx**: docs/theory/collision_probability.rst §Inner vs Outer Iterations (Historical Note)
**Error catalog**: tests/l0_error_catalog.md ERR-016

The original inner convergence check computed
`||Σ_t V φ_new - P^T V Q||₂`, which is identically zero because
`φ_new = P^T V Q / (Σ_t V)` by construction.  The inner loop
always exited after 1 iteration, making inner iterations vacuous.

Fixed to use relative flux change: `||φ_new - φ_old|| / ||φ_new||`.
With the corrected residual, thermal groups with strong self-scatter
(`Σ_s(g→g) / Σ_t(g)` large) require multiple inner iterations, while
fast groups converge in 1.  This matches the expected physics.

Tests: `test_cp_verification.py::TestGSInnerIterations` — 3 tests
verifying thermal > fast inner iterations, GS/Jacobi eigenvalue
agreement, and no-self-scatter convergence in 1 iteration.

### CP-20260405-006 — Fix compute_keff for (n,2n) (C-2, ERR-015)

**Status**: DONE
**Sphinx**: docs/theory/collision_probability.rst §Power Iteration (keff formula + (n,2n) note)
**Error catalog**: tests/l0_error_catalog.md ERR-015

The eigenvalue estimate `compute_keff` used `production / absorption`
where `production = νΣf·φ·V` and `absorption = Σa·φ·V`.  With nonzero
(n,2n), this is wrong: the (n,2n) reaction produces an extra neutron
but `absorption_xs` already counts the (n,2n) removal.  The correct
balance is:

    k = νΣf·φ·V / (Σt - Σs - 2·Σ₂)·φ·V

Fixed to compute net removal = total - scatter - 2*(n,2n).  When
Sig2=0, `total - scatter = absorption` so the formula reduces to the
original.

Tests: `test_cp_verification.py::TestN2N::test_n2n_solver_keff_matches_analytical`
now passes (was xfail).

### CP-20260405-007 — Consolidated eigenvalue solvers

**Status**: DONE
**Sphinx**: docs/theory/collision_probability.rst §Consolidated Eigenvalue Solvers

Created `derivations/_eigenvalue.py` with two shared functions:

- `kinf_homogeneous(sig_t, sig_s, nu_sig_f, chi, sig_2=None)` —
  replaces 7 copy-pasted eigenvalue computations across
  `homogeneous.py`, `sn.py`, `moc.py`, `mc.py`.
- `kinf_from_cp(P_inf_g, ..., sig_2_mats=None)` — replaces 3
  identical `_kinf_from_cp` definitions in `cp_slab.py`,
  `cp_cylinder.py`, `cp_sphere.py`.

Both support optional (n,2n) via `sig_2` / `sig_2_mats` parameters.
Extended `make_mixture` with `sig_2` parameter.

All 106 CP tests pass.  All existing eigenvalue references unchanged
(sig_2=None preserves previous values).

---

## IMPL — Implemented, Sphinx Documentation Pending

### CP-20260405-005 — Comprehensive QA-driven verification suite

**Status**: IMPL

31 new tests in `tests/test_cp_verification.py` closing 9 coverage
gaps identified by QA review:

- G-1: L0 P_inf element-by-element comparison (solver vs derivation)
- G-2/W-1: Upscatter eigenvalue and regression guard
- G-3/W-3: (n,2n) eigenvalue tests
- G-4/G-7: Optically thick/thin stress tests
- G-5: Convergence rate (monotonic error decrease, dominance ratio)
- G-6/W-2: Multi-group CP matrix properties (row sums, reciprocity)
- G-8: 8-region slab + mesh refinement convergence
- G-9: GS inner iteration physics verification
- W-6: Ki4 table resolution convergence

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
