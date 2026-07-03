---
name: cp-matrix-density-and-sphere-conservation
description: CP [P] is structurally dense (no exact zeros) — verdict + probe recipe for the #226 naming question; AND the discrete-CP spherical kernel violates the proven row-sum=1 conservation law above tau~3 (rowsum saturates at 1/pi), invisible to the tau<=2 test suite.
metadata:
  type: project
---

# CP collision-probability matrix: density verdict + spherical conservation defect

Two findings from the #226 taxonomy density hypothesis test (2026-07-01).
Probe: `diag_cp_density_probe.py` (isolates [P] via `CPMesh.compute_pinf_group(sig_t_g)`
— needs only mesh + per-cell sig_t, NO Mixture / eigenvalue solve; vacuum BC returns
pure first-flight P_cell, white BC returns P_inf = P_cell + rank-1 outer product).

## 1. [P] is STRUCTURALLY DENSE (answers the DenseInverseOperator vs MatrixInverseOperator name)

- Assembly: `CPMesh.compute_pinf_group` (`orpheus/cp/solver.py:212`) → dense `np.ndarray (N,N)`,
  filled ALL-PAIRS by nested `for i for j` loops (`_compute_slab_rcp` E3 slab :234 /
  `_compute_radial_rcp` Ki3-cyl/exp-sph :284), normalized :363, then white rank-1 add
  (`_apply_white_bc` :376) or vacuum copy. Packed dense `(N,N,ng)` at :864.
- Every P_ij is a second-difference of E3/Ki3/exp over the optical gap → strictly > 0,
  decaying ~exp(-tau_ij). Measured (slab/cyl/sph, N=3&25, tau_total 0.25..25):
  **ZERO exact zeros in every config** (even tau=25: far-corner ~1e-12, never 0.0),
  **zero negatives**, frac<1e-15*max = 0. Reciprocity `diag(SigT*V)@P` symmetric to ~1e-16.
- Interface-current / response-matrix SPARSE variant: does NOT exist in code; **planned in
  OPEN GitHub Issue #56** ("CP: Interface current method for multi-cell lattice"). The CP
  "sparsity" the theory page mentions is the SCATTERING matrix (downscatter-only), never [P].
- **Verdict: storage-agnostic `MatrixInverseOperator`.** Density is not a promise CP can
  break — it is dense by construction. A `Dense*` name would be a redundant storage claim on
  a matrix that is inherently dense; if the #56 interface-current variant ever lands with a
  genuinely sparse (banded/adjacent-coupling) operator, a storage-neutral base name avoids a
  future rename.

## 2. OPEN DEFECT — spherical CP violates row-sum=1 above tau~3 (candidate GitHub issue)

- Theory page (`collision_probability.rst:520`) PROVES "Row sum property: sum_j P_inf_ij = 1"
  with NO tau caveat. Slab & cylinder satisfy it to ~1e-6 at tau=25. **Sphere does NOT.**
- Pure P_cell interior-cell rowsum SATURATES at **exactly 1/pi = 0.318310** (ratio 1.0000,
  mesh- and tau-independent for tau>~10) instead of →1 as a near-black interior demands.
  Downstream, the white-BC closure assumes complementarity `escape = 1 - rowsum`, so the
  1/pi under-collision → escape mis-set to ~0.68 → rank-1 correction blows P_inf rowsums up
  to ~4.5 (conservation residual 3.556 at tau=25). Onset ~tau_total>3; perfect (~1e-16) below.
- NOT quadrature underconvergence: `n_quad_y` sweep 64→1024 is DEAD FLAT at 3.556; residual
  GROWS with radial refinement (3.556→3.792 as N 25→100 at fixed tau) — Step-1 "error grows
  under refinement" = balance/normalization inconsistency, not underconvergence. Likely a
  missing ~pi solid-angle factor in `_setup_spherical` (`y_wts*=y_pts`) / `_compute_radial_rcp`
  spherical path — present as a clue, not a confirmed fix.
- INVISIBLE to the suite: `tests/cp/test_properties.py::test_row_sums` +
  `test_verification.py::test_row_sums_multigroup` DO assert rowsum=1 for spherical, but every
  fixture is R=1, SigT<=2 → tau<=2, below onset. All 9 `cp_sph1D_*` eigenvalue cases are tau<=2
  too. Classic vv-principles H2: the conservation gate runs only where the invariant happens to
  hold. A tau>=5 spherical row-sum case would red immediately (the probe's sph-thick config IS
  that isolating gate — promote it).
- eigenvalue tests pass because k is a thin-regime ratio insensitive to the thick-limit
  normalization; do NOT read their green as covering thick spheres.
