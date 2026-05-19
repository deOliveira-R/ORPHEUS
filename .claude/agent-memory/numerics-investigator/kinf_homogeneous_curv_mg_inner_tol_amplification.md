---
name: kinf-homogeneous-curv-mg-inner-tol-amplification
description: test_kinf_homogeneous sphere/cylinder × 2eg/4eg drift (~1.2e-7) is SI inner-tolerance amplification, NOT a solver bug. Slab and 1eg pass at default because of low SI dominance ratio / shape-independence respectively. Fix: tighten inner_tol or relax rtol.
metadata:
  type: project
---

`tests/sn/l1_analytical/test_kinf_homogeneous` 4-case failure (sphere/cylinder × 2eg/4eg) at rtol=1e-9 (drift 1.2e-9 to 2.7e-8) on branch `refactor/sn-operator-algebra` (HEAD `43bb8e0`) is **convergence-tolerance**, NOT a solver bug.

**Why:** Mechanism is the well-known SI residual ≠ gap-to-fixed-point relation. Within-group SI termination `‖φ_n − φ_{n−1}‖ / ‖φ_n‖ < inner_tol` bounds the residual, but the gap to the true within-group fixed point is `gap ≈ ρ/(1−ρ) × residual`, where ρ = within-group SI spectral radius ≈ `c_within = Σ_(s,g→g) / Σ_(t,g)`. For the 2eg test problem the thermal group has c_within = 0.9 → amplification 9; power iteration coupling re-amplifies into k_eff.

**How to apply:**
1. Slab passes at default tolerances because its SI dominance ratio is much smaller (Cartesian sweep convergence is fast). Sphere/cylinder are slower.
2. 1eg passes everywhere because k = νΣ_f/Σ_a is flux-shape independent (no inner-iteration coupling).
3. Tightening `inner_tol` from 1e-9 → 1e-12 collapses every curvilinear-MG drift to ~1e-12 (factor 10⁵ reduction). Tightening `keff_tol` alone has NO effect — the bias is in the inner solve, not the outer keff convergence test.
4. `inner_solver="krylov"` at default tolerances reduces drift to ~4e-10 — strong evidence the L operator is correct; the gap was purely in SI's convergence rate.
5. Spectrum eigenvector cross-check at tight tolerances matches the analytical `A^{-1}F` dominant eigenvector to 6e-12 (machine ULP) — no scattering-transpose / sign-flip / variable-swap bug.

**Decision matrix (verdict locked, fix is user choice):**
- **Option A — solver default change**: bump default `inner_tol` in `solve_sn` from `1e-8` to `1e-11`. Curvilinear-MG users get correctness for free; cost is ~10× more inner iterations per outer in slow-converging problems. RECOMMENDED if user wants verification correctness by default.
- **Option B — test-side fix**: in `test_kinf_homogeneous`, pass `inner_tol=1e-12` (or use `inner_solver="krylov"`). Minimal blast radius; respects "test calibrates to solver tolerances per fixture".
- **Option C — rtol relaxation**: bump test rtol from 1e-9 to 5e-7. NOT RECOMMENDED — papers over the bias rather than pinning the convergence-tolerance model; loses sharpness against future bug regressions.

**Diagnostic test promoted:** `derivations/diagnostics/diag_kinf_curv_mg_02_inner_tol_amplification.py` (8 tests, ~6 min total) pins both (i) the drift collapses at tight tolerance and (ii) Krylov matches at default. Both are structural-independence pins against future SI/Krylov drift. Promote to `tests/sn/l1_analytical/test_kinf_homogeneous_tolerance.py` as a companion regression gate to `test_kinf_homogeneous` if Option A or B is adopted.

**Step 3 (R-1 migration) implication:** The 1.2e-7 drift is NOT carried forward as a real bug by the operator-typed migration. The migration is safe to land. After Step 3 lands, re-running the spectrum cross-check at tight tolerance verifies the typed `AngularFlux` does not perturb the within-group fixed point — but that's a regression check, not a fix prerequisite.

Related:
- [[issue-196-phase-g-step2-cylinder-minimal-reproducer]] (a real curvilinear bug, ERR-048 manifestation #3 territory — different shape)
- [[phase-f-step3-diagnostic]] (Phase F SN sweep Carlson seed fix — solved a related SI/Krylov disagreement that WAS a real bug)
