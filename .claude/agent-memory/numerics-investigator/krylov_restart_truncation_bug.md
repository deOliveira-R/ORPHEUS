---
name: krylov-restart-truncation-bug
description: KrylovAcceleration GMRES restart=50 hardcode + discarded scipy info flag causes silent wrong-keff in _solve_krylov on sphere homogeneous reflective n_cells>=10. Root cause is structural (restart truncates Krylov subspace), not tolerance. Fix is one-line at iteration.py:716 + an info-flag check at 735.
metadata:
  type: project
---

# Krylov GMRES restart truncation bug (2026-05-28 investigation)

**Affected code path**: `SNSolver._solve_krylov` (`orpheus/sn/solver.py:599`) → `KrylovAcceleration.solve` (`orpheus/numerics/iteration.py:656`) → `scipy.sparse.linalg.gmres`.

**Symptom**: `solve_sn(..., inner_solver="krylov")` on homogeneous reflective sphere (mixture A 2g, n_cells=20, GL n_ord=8) returns keff=1.4019 vs analytical k_inf=1.875 (25% error). SI inner solver returns 1.875 correctly.

## Root cause — structural (NOT tolerance)

**`orpheus/numerics/iteration.py:716`** (clamp expression in `_solve_krylov`):
```python
restart=min(50, N * ng * nx * ny),
```
With N=8, ng=2, nx=20, ny=1: `N*ng*nx*ny = 320 ≫ 50` → restart clamped to **50**. The full ravel size including the B1'' face block is **328** unknowns. GMRES with restart=50 truncates the Krylov subspace to 50 vectors. At restart cycles >= 200 (≈problem size) the system would converge to machine precision; at restart=50 GMRES *cannot reach* the true solution regardless of `maxiter`.

**`orpheus/numerics/iteration.py:735`** (`spla.gmres` return):
```python
solution, _info = spla.gmres(...)
```
`_info` (leading underscore = discarded). scipy returns `info > 0` when GMRES does NOT converge within `maxiter`. The KrylovAcceleration silently accepts an unconverged `solution` and returns it as if converged.

## Evidence chain

1. **Step 1** (tolerance falsifier): tightening `inner_tol` 1e-8 → 1e-12 yields **identical** keff=1.4019239124. NOT a tolerance issue. The fixed point of the GMRES iteration is stable at the WRONG value.

2. **Step 5** (mesh-refinement scaling, default tols):

   | n_cells | SI keff | Krylov keff | Krylov err |
   |---|---|---|---|
   | 5 | 1.875 | 1.875 | 4e-10 |
   | 8 | 1.875 | 1.875 | 1e-9 |
   | 10 | 1.875 | 1.9153 | 4e-2 |
   | 12 | 1.875 | 1.8589 | 2e-2 |
   | 16 | 1.875 | 1.5987 | 3e-1 |
   | 20 | 1.875 | 1.4019 | **5e-1** |

   Error GROWS with refinement. Classic discretisation-inconsistency fingerprint (vv-principles step 1).

3. **Step 6** (tight tolerances): n_cells=10 + `inner_tol=1e-12, max_inner=1000` recovers k_inf to 1e-12. Same tight tols at n_cells=20 still gives err=3.2e-1. So tight tols rescue the L1 anchor's `n_cells=10` but NOT a moderately refined mesh.

4. **Step 8** (definitive — restart sweep at n_cells=20):

   | restart | scipy `info` | true `||r||` | keff |
   |---|---|---|---|
   | 50 (hardcode) | **200 (FAIL)** | 1.27 | 1.5499 (wrong) |
   | 100 | 0 | 9e-13 | **1.8750** |
   | 200 | 0 | 1e-12 | **1.8750** |
   | 320 | 0 | 1e-12 | **1.8750** |

   At restart >= 100 (~1/3 problem size) the fix dissolves the bug completely; at restart=50 GMRES literally cannot converge (info=200).

## Why the L1 anchor passes at n_cells=10 but the failing test doesn't at n_cells=20

The L1 anchor `test_identity_preconditioner_recovers_kinf[sphere]` and `test_kinf_homogeneous[sphere-2eg-krylov]` use **n_cells=10**. At n=10 with `inner_tol=1e-12, max_inner=1000`, the under-converged GMRES output happens to project onto the true eigenmode well enough (it converges to keff=1.875 to 1e-12 within 2 outer iterations). This is a numerical *coincidence* — the operator at n=10 is benign enough for restart=50 GMRES to converge accidentally. At n=20 the operator's Krylov subspace dimension required for convergence exceeds 50, the GMRES output is essentially noise, and the power iteration latches onto a spurious fixed point.

The failing test `test_solve_sn_si_vs_krylov_consistency_homogeneous_sphere` uses **n_cells=20** AND default tols (`inner_tol=1e-8, max_inner=200`). Both choices make the bug visible.

## Minimal fix (production code)

**At `orpheus/numerics/iteration.py:716`** (in `_solve_krylov` and equivalent call sites):
```python
# OLD: restart=min(50, N * ng * nx * ny),
restart=N * ng * nx * ny,  # full problem size — no Krylov truncation
```
OR better: make `restart` adaptive to problem size with a sensible floor:
```python
restart=min(self.max_inner, max(50, N * ng * nx * ny)),  # at minimum problem size
```

**At `orpheus/numerics/iteration.py:735`**: stop discarding the convergence flag:
```python
solution, info = spla.gmres(...)
if info > 0:
    warnings.warn(
        f"GMRES did not converge: info={info}, last residual estimate "
        f"{residual_history[-1]:.3e} vs rtol={self.tol:.3e}. "
        f"Solution may be unreliable.",
        ConvergenceWarning,
    )
```

## API contract

The fix is **purely internal**. `keff_tol`/`flux_tol`/`inner_tol` semantics are preserved. `restart` should arguably become a public kwarg of `solve_sn` (currently buried inside `_solve_krylov`).

## Cross-reference

- Test file: `tests/sn/spatial/test_sweep_vs_apply_consistency.py::test_solve_sn_si_vs_krylov_consistency_homogeneous_sphere`
- L1 anchor that passes (by coincidence): `tests/sn/test_krylov_curvilinear_precond_safety.py::test_identity_preconditioner_recovers_kinf[sphere]`
- L1 anchor that passes (by tightness): `tests/sn/l1_analytical/test_kinf_homogeneous.py::test_kinf_homogeneous[sphere-2eg-krylov]`
- The 2026-05-26 max_inner bump from 300 → 1000 in `_TIGHT_KW` was a symptom-level patch for the same underlying bug at n_cells=10.

## Recommended ERR entry (post-fix)

Failure mode #4 (convention drift / hardcoded constant — see ERR-004 family, ERR-025). The `restart=50` hardcode is the same anti-pattern as the `4π` constants ERR-004: a magic constant baked in that breaks when problem dimensions change. Per `coding-elegance` Pattern 7 (Normalise at definition site), Pattern 4 (Make illegal states unrepresentable), the restart should be computed from the operator's domain size, not hardcoded.
