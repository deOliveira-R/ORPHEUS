# Gotchas — Error Log and Lessons Learned

Errors encountered during development, classified by the 6 AI Failure Modes taxonomy.
Each entry records the bug, how it hid, and the test that now prevents regression.

---

## #1 — Z-ordinate weight loss in 2D Lebedev sweep

**Failure mode:** #4 Factor error — missing contribution

**Date:** 2026-04-03

**Bug:** In `sn_sweep.py::_sweep_2d_wavefront`, ordinates with `mu_x = mu_y = 0`
(z-directed, pointing out of the 2D plane) were skipped with `continue`.
Their quadrature weights (2 ordinates, total 0.096 out of 12.566 = 0.77%)
were lost from the scalar flux integration `phi = sum(w_n * psi_n)`.

**Impact:** Multi-group eigenvalue error of ~0.4%.
- 2G homogeneous: keff = 1.867 vs analytical 1.875
- 421G PWR demo: keff = 1.035 vs MATLAB reference 1.042

**Why it hid:**
- *"It produces reasonable numbers."* 0.7% error is easy to dismiss as discretization or quadrature approximation.
- *"The integration test passes."* 1-group homogeneous gave **exact** keff = 1.500000 because `k = nu_Sigma_f / Sigma_a` is independent of flux shape — the weight loss scales all groups equally and cancels in the ratio. Only multi-group problems have a group-ratio eigenvector that the weight loss distorts.
- *"The convergence rate is correct."* Mesh refinement (n=2,4,8,16) showed keff constant at 1.867 — properly converged spatially, to the wrong angular answer.

**Fix:** Handle z-directed ordinates as pure-collision: `psi = Q * weight_norm / Sigma_t`.
No streaming in x or y, so the diamond-difference simplifies to a division.

**Tests that catch it:**
- `TestQuadratureWeightConservation::test_no_weight_lost` — reconstructs scalar flux from angular flux, verifies all ordinates contribute
- `TestQuadratureWeightConservation::test_z_ordinates_contribute` — asserts z-directed ordinates have nonzero angular flux
- `TestQuadratureWeightConservation::test_homogeneous_scalar_flux_equals_Q_over_sigt` — exact magnitude comparison: converged phi must equal Q/Sigma_t for uniform source
- `TestMultiGroupEigenvector::test_2g_eigenvector` — verifies converged group ratio matches analytical eigenvector (the quantity 1-group tests cannot check)
- `TestHomogeneousExact::test_homogeneous_exact[2g-2G]` — 2G keff must match analytical to 1e-8

**Lesson:** Always test multi-group problems independently of single-group. 1-group success gives false confidence because keff is insensitive to angular/spatial flux errors that cancel in the production/absorption ratio.

---

## #2 — Scattering matrix transpose in vectorization

**Failure mode:** #2 Variable swap — `Sigma_s` vs `Sigma_s^T`

**Date:** 2026-04-03

**Bug:** When vectorizing `_add_scattering_source` from per-cell loops to
batched matrix operations, first attempt used `phi @ self._sig_s0_T[mid]`
(pre-transposed matrix) instead of `phi @ self._sig_s0[mid]`.

The original code: `sig_s0.T @ phi_vec` (matrix-vector)
The batched equivalent: `phi_batch @ sig_s0` (row-vector times matrix)
The wrong version: `phi_batch @ sig_s0.T` (double-transposed)

**Impact:** keff = 2.06 (completely wrong) — caught immediately.

**Why it could hide in other contexts:** For symmetric scattering matrices
(e.g., 1-group self-scatter), `Sigma_s = Sigma_s^T` and the bug is invisible.
Only asymmetric multi-group downscatter exposes it.

**Fix:** Use `phi @ sig_s0` (not transposed). The identity
`(A^T @ v) = (v^T @ A)^T` means row-vector times un-transposed matrix
gives the same result as transposed matrix times column vector.

**Test that catches it:**
- `TestAddScatteringSource::test_matches_reference` — compares vectorized output against per-cell reference with random asymmetric flux

**Lesson:** When converting `A^T @ v` to batched form, the transpose moves
to the "other side" of the `@` operator, not into a pre-transposed matrix.
Always test with asymmetric inputs where swap produces a detectably different value.

---

## #3 — Octant batching breaks reflective BC ordering

**Failure mode:** #6 Convention drift — implicit ordinate processing order

**Date:** 2026-04-03

**Bug:** Attempted to batch ordinates by sweep direction (same sign of mu_x, mu_y)
to reduce the Python loop from 110 to 4 iterations. The reflective BC copies
incoming flux from the reflected partner ordinate:
```python
psi_x[n, 0, :, :] = psi_x[ref_x[n], 0, :, :]
```
In the sequential loop, the reflected partner's outgoing flux is from the
*previous* sweep (may or may not be updated yet in this sweep, depending on
ordinate ordering). Batching changed this ordering: all ordinates in a group
are swept simultaneously, so a group reads boundary fluxes that its reflected
group hasn't updated yet in this sweep.

**Impact:** 2G convergence test failed — keff diffs grew instead of decreasing.
1G appeared to work (reflective BCs converge faster for 1-group).

**Why it hid initially:** The optimization produced reasonable-looking keff
for the 1G test and the per-sweep timing improved 2x. The failure only
appeared in the 2G mesh convergence test (which takes minutes to run).

**Fix:** Reverted octant batching. The sequential ordinate loop is load-bearing
for reflective BC consistency. Future optimization should either:
(a) process reflected octant pairs together (ensure partner is updated first), or
(b) use Numba JIT on the sequential loop to eliminate Python overhead.

**Test that catches it:**
- `TestHomogeneousExact::test_homogeneous_exact[2g-2G]` — exact eigenvalue comparison
- `TestTransportSweep::test_matches_saved_reference` — bitwise regression against saved reference

**Lesson:** When a loop has implicit data dependencies between iterations
(ordinate N reads ordinate M's output), batching changes the semantics.
The reflective BC copy is a "closure defines X, caller assumes X" pattern —
the sequential processing order is part of the interface contract, not an
implementation detail.

---

## #4 — Hardcoded 4π in BiCGSTAB RHS normalization

**Failure mode:** #4 Factor error — wrong angular normalization constant

**Date:** 2026-04-03

**Bug:** `sn_operator.py::build_rhs` hardcoded `four_pi = 4.0 * np.pi` to normalize
isotropic sources for the angular flux equation. This is correct for Lebedev
quadrature (`sum(weights) = 4π`) but wrong for Gauss-Legendre 1D quadrature
(`sum(weights) = 2`). The normalization must use `sum(weights)`, which is
quadrature-dependent.

**Impact:** BiCGSTAB on 1D GL quadrature diverged (keff oscillating 1.9→4.5→2.8).
BiCGSTAB on 2D Lebedev worked correctly (the hardcoded value happened to be right).

**Why it hid:** All initial BiCGSTAB testing used Lebedev quadrature, where
`4π` is the correct value. The bug only manifested when testing the 1D GL path.

**Fix:** Replace `four_pi` with `float(quad.weights.sum())` in both `build_rhs`
and the caller in `_solve_bicgstab`.

**Tests that catch it:**
- `TestBicgstabNormalization::test_1d_gl_homogeneous_exact` — 1D GL must match analytical k_inf
- `TestBicgstabNormalization::test_2d_lebedev_homogeneous_exact` — 2D Lebedev must match
- `TestBicgstabNormalization::test_gl_and_lebedev_agree` — both quadratures must produce same keff

**Lesson:** Never hardcode physical constants that depend on the quadrature rule.
The angular normalization `1/sum(w)` is a property of the discretization, not a
universal constant. Testing with only one quadrature type masks this.

---

## #5 — DD recurrence formula rewrite breaks multi-group convergence

**Failure mode:** #5 Algebraically-equivalent rewrite — numerically unstable form

**Date:** 2026-04-04

**Bug:** During the geometry module migration, `_solve_recurrence` and `_outgoing`
in `sn_sweep.py` were rewritten from:
```python
psi_out = a * psi_in + s
return 0.5 * (psi_in + psi_out)
```
to the algebraically equivalent:
```python
return 2.0 * (psi_in * (1.0 - a) / (1.0 - a + 1e-300) + s) - psi_in
```
Additionally, `psi0[None, :]` lost its broadcasting dimension in the cumulative
product computation, and `_outgoing` was changed to use a different formula
that computes the outgoing flux from the last cell rather than simply
propagating `cp[-1] * (psi0 + cs[-1])`.

**Impact:** Scattering source iteration diverged for multi-group problems.
The flux grew by ~1e34 per outer iteration. 1-group problems were unaffected
(the spectral radius of the scattering iteration is <1 for 1-group).

**Why it hid:**
- Both formulas are algebraically identical for exact arithmetic.
- The 2D wavefront sweep (which was also updated) used a different code path
  and wasn't affected — it still passed the bitwise regression test.
- The 1G homogeneous test passed because keff = νΣf/Σa is independent of
  the absolute flux level.
- The bug only manifested through multi-group scattering iteration divergence,
  which required running the full power iteration (not just a single sweep).

**Fix:** Restored the original formulas:
- `psi_out = a * psi_in + s; return 0.5 * (psi_in + psi_out)`
- `return cp[-1] * (psi0 + cs[-1])` for outgoing flux
- `psi0[None, :]` broadcasting preserved

**Tests that catch it:**
- `TestSolveRecurrence::test_multi_cell_dd_relation` — verifies cell-average = 0.5*(psi_in + psi_out) for every cell
- `TestSolveRecurrence::test_outgoing_matches_last_cell` — outgoing flux must equal last cell's psi_out
- `TestSolveRecurrence::test_regression_multigroup_scattering_convergence` — 2G flux must remain bounded after one outer iteration

**Lesson:** "Algebraically equivalent" does not mean "numerically equivalent".
The DD recurrence is a coupled forward-solve where small per-cell differences
compound over many cells. The stable form `0.5*(psi_in + psi_out)` is a
direct average of known quantities, while the rewritten form subtracts
nearly-equal large numbers (`2*(...) - psi_in`). Always regression-test
the recurrence with multi-group data before refactoring sweep internals.

---

## Meta-lessons

1. **Multi-group is the minimum bar.** 1-group problems are degenerate — they
   hide factor errors, convention drifts, and weight losses that cancel in
   the scalar keff ratio. Every transport solver must be verified at 2+ groups.

2. **Exact analytical tests first, convergence tests second.** The homogeneous
   infinite-medium eigenvalue is known analytically. Testing convergence (diffs
   decrease with mesh refinement) only proves the scheme converges — not that
   it converges to the right answer.

3. **Unit test before optimizing.** The scattering transpose bug was caught in
   seconds by the per-cell reference test. The octant batching bug took hours
   to diagnose because we optimized first and tested after. Write the reference
   test, verify it passes on the original code, then optimize.

4. **"Constant with mesh refinement" means angular error, not spatial.**
   If keff doesn't change as h→0, the error is in the angular discretization
   (quadrature weights, ordinate handling), not in the spatial scheme.
