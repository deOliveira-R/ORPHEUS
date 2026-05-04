# Peierls unification close-out — fresh-session handoff

**Branch:** `investigate/peierls-solver-bugs` (4 commits ahead of `feature/mms-broad`).
**Author of this plan:** Claude Opus 4.7, 2026-04-19.
**Scope:** Close out GitHub issues **#115** (retire `peierls_slab._basis_kernel_weights`) and **#116** (retire Bickley infrastructure completely) in one focused session. Commit the branch when both are done.

**Audience:** fresh Claude Code session that has NOT seen the investigation that preceded this. Read sections 0–3 in order; the investigation history and architectural reasoning are the expensive context you need before touching the code.

---

## 0. One-paragraph summary (read first)

Two bugs in the Peierls Nyström K-matrix assembly were discovered on 2026-04-18 and fixed in two separate modules (`peierls_slab.py` adaptive mpmath.quad, `peierls_geometry.py` ρ/ω subdivision). The user flagged that **splitting the fix across two code paths is a symptom of an un-unified framework** and asked for forward-only unification per the existing plan at `.claude/plans/post-cp-topology-and-coordinate-transforms.md`. Two more commits added `slab-polar` (math proven equivalent but quadrature converges only algebraically) and `cylinder-polar` (volume kernel Ki₁-free, tests pass). Remaining work: reach machine precision in slab-polar (#115) and retire the last Ki_n calls in `compute_G_bc` / `escape_kernel_mp` (#116). When both land, the unified architecture is complete and the legacy paths can be deleted.

---

## 1. How we got here (the original problem)

### 1.1 The original bug report

On 2026-04-18 `numerics-investigator` logged `.claude/plans/post-peierls-solver-bugs.md` (the original plan for this branch). Two independent bugs were found in the Peierls K-matrix assembly:

| Bug | Location | Symptom |
|-----|----------|---------|
| Slab cross-panel | `peierls_slab._build_kernel_matrix:155-164` | Cross-panel K[i,j] entries use one-point collocation `E_1(τ_ij)·w_j/2` — ~1% error at p=4, O(h) convergence |
| Curvilinear ray-crossing | `peierls_geometry.build_volume_kernel:491-516` | Fixed-GL on (β/θ, ρ) does not subdivide at ρ values where r'(ρ) crosses panel boundaries — 0.5–5% K errors, non-monotonic convergence in n_ρ |

Both bugs were systematically invisible to every **row-sum / conservation** test because `∑_j L_j(r') = 1` (partition of unity) kills the kink in the summed integrand even when each basis-individual K[i,j] is wrong. The rank-N investigation's "K·[1]=1 passes at mode 0, fails at mode n>0" signature was exactly this.

GitHub issues filed: **#113** (slab) and **#114** (curvilinear).

### 1.2 The bug fixes (commits 9ac41de and aa98565)

**[9ac41de] slab fix.** Replaced the naive cross-panel collocation AND the same-panel singularity-subtraction path with a **unified basis-aware `mpmath.quad`** integrator. Same technique used by the adaptive reference `peierls_reference.slab_K_vol_element`. ERR-027 and ERR-028 logged in `l0_error_catalog.md`.

**[aa98565] curvilinear fix.** Added `rho_crossings_for_ray` and `omega_tangent_angles` methods to `CurvilinearGeometry`. `build_volume_kernel` now subdivides the ρ integration at panel crossings AND the ω integration at tangent-to-interior-boundary bifurcation angles. Also added `ki_n_float` (scipy-based fast Ki₁) to keep cylinder tests under 1 minute with the 6× subdivision overhead. ERR-029 logged.

Both commits include `Closes #NN` trailers — the issues auto-close on merge to main.

### 1.3 The user's pushback (the architectural concerns)

After the bug fixes landed, the user raised two concerns:

1. **Frameworks not integrated.** Two separate fixes for one bug means two code paths that can drift. Plan `post-cp-topology-and-coordinate-transforms.md` §4 has an integrating strategy (slab as a `CurvilinearGeometry` kind) that had not been used.

2. **Bickley is redundant.** The plan's §0.2 makes clear that Ki₁ / Ki₂ / Ki₃ are just pre-integrated exponentials — any fast path for Ki₁ shouldn't exist at all if we expose the underlying angle explicitly. The speculative `ki_n_float` added in aa98565 was "solving the wrong problem."

The user gave the forward-only direction: **build the new integrated case, test it, retire old pathways once tested**. Do not revert.

### 1.4 The architectural principle (from plan §0.2 + §4)

Every Peierls kernel — slab E₁, cylinder Ki₁, sphere `exp(-τ)` — is the result of integrating out symmetry coordinates from the same 3-D isotropic point kernel:

```
        3-D point kernel e^{-Σ_t R} / (4π R²)
           ↓ integrate over symmetry (2/1/0 dims)
    ───────────────────────────────────────────────────
    Level 1 (pointwise Peierls):  E₁/2    Ki₁/(2π)    e^{-τ}/(4π)
```

Under the **τ-coordinate transform** (plan §5): `τ ≡ ∫₀^ρ Σ_t(s) ds`, the fundamental kernel becomes `e^{-τ}` for every geometry. Ki₁ is not a separate kernel — it's what you get if you pre-integrate the out-of-plane angle φ:

```
    Ki_1(τ) = ∫_0^{π/2} e^{-τ/cos φ} dφ
    Ki_2(τ) = ∫_0^{π/2} cos(φ) e^{-τ/cos φ} dφ
```

So if the cylinder Peierls quadrature INCLUDES an explicit φ integration, the kernel collapses back to `exp(-τ/cos φ)` and no Bickley evaluator is needed.

The unified architecture: one `CurvilinearGeometry` abstraction, one `build_volume_kernel`, one `exp` kernel (with an extra φ dimension for cylinder). Fix once, fix everywhere.

### 1.5 Commits 4395cb8 and 213278f

**[4395cb8] `slab-polar`.** Added `CurvilinearGeometry(kind="slab-polar")` with the plan's observer-centred polar form:

```
Σ_t(x) φ(x) = (1/2) Σ_t(x) ∫_{-1}^{1} dμ ∫_0^{ρ_max(x,μ)} e^{-Σ_t ρ} q(x + ρμ) dρ
```

- Direction cosine μ = angular variable (`ray_direction_cosine(μ) = μ` identity)
- `volume_kernel_mp = exp(-τ)` (same as sphere!)
- Linear `rho_max`, `source_position`, `optical_depth_along_ray`, `rho_crossings_for_ray`
- `omega_tangent_angles = []` (no turning-point geometry)

Plus `slab_polar_K_vol_element` in `peierls_reference.py` — an **adaptive mpmath.quad reference** that proves the polar form and the classical E₁ form are bit-identical (machine precision).

Production path `_build_volume_kernel_slab_tau_laguerre` uses:
- Outer: exp-stretched `v = -ln|μ|` + Gauss-Laguerre
- Inner: τ-GL with subdivision at panel crossings + dyadic τ-cap

**Convergence: algebraic ~O(N⁻²).** At n_angular=n_rho=64, max rel err vs legacy E₁ is 2.3e-4. **Good enough to prove unification works; too coarse to retire the legacy E₁ Nyström yet.** Issue #115 tracks this.

**[213278f] `cylinder-polar`.** Added `CurvilinearGeometry(kind="cylinder-polar")` with explicit out-of-plane φ integration:

```
K[i, j] ← outer_w × ρ_w × (Σ_{k} φ_w_k · exp(-τ/cos φ_k)) × L_j
```

- Spectral convergence of φ-GL to Ki₁ (verified: 1.7e-15 at n_phi=32, τ=5)
- Production `peierls_cylinder.GEOMETRY` now points at `CYLINDER_POLAR_1D`
- All 15 cylinder eigenvalue+BC tests pass unmodified
- All CP tests (115) pass

But the volume kernel isn't the only Ki_n user. `compute_G_bc` still calls Ki₁; `escape_kernel_mp` still calls Ki₂. Issue #116 tracks completing the retirement.

---

## 2. Current state (what's done, what's left)

### 2.1 Branch state

```
investigate/peierls-solver-bugs
    ↑ 213278f  feat: cylinder-polar (#116 volume-kernel phase)
    | 4395cb8  feat: slab-polar kind
    | aa98565  fix: curvilinear ρ/ω subdivision (closes #114)
    | 9ac41de  fix: slab E₁ assembly (closes #113)
feature/mms-broad (base)
```

Not yet merged to main. `Closes #113` and `Closes #114` will fire on push to main.

### 2.2 GitHub issues

| # | Status | Title | Notes |
|---|--------|-------|-------|
| #113 | OPEN* | Slab cross-panel bug | *auto-closes on merge (Closes trailer in 9ac41de) |
| #114 | OPEN* | Curvilinear ρ-crossing bug | *auto-closes on merge (Closes trailer in aa98565) |
| #115 | OPEN | Retire `peierls_slab._basis_kernel_weights` | Needs slab-polar at ~1e-8 precision |
| #116 | OPEN | Retire Bickley | Volume kernel done; G_bc + escape_kernel remaining |

Both #115 and #116 have comments with detailed technical findings — read those first in the fresh session.

### 2.3 Test status (last verified 2026-04-19)

| Suite | Result |
|-------|--------|
| `test_peierls_reference.py` (all L1 equivalence) | 19 pass (E₁ row-sum, element-wise, sphere K, slab-polar ref, slab-polar prod, cylinder-polar equivalence) |
| `test_peierls_sphere_*` (all 4 files) | 35 pass |
| `test_peierls_cylinder_*` (all 5 files) | 29 pass |
| `test_peierls_rank_n_*` + closure | 77 pass + 11 xfailed |
| `test_peierls_convergence.py` | 5 pass |
| `tests/cp/` | 115 pass |

No regressions from any of the four commits. Full derivations suite times out at 15 min with the new cylinder-polar path (adds 16× work in the cylinder hot loop) — split into batches if running end-to-end is needed.

### 2.4 Key files (line numbers current as of 213278f)

```
orpheus/derivations/peierls_geometry.py
  line  59-60  : imports ki_n_float, ki_n_mp (to be removed after #116)
  line 178    : @dataclass CurvilinearGeometry
  line 276    : ray_direction_cosine  (identity for slab, cos for curvilinear)
  line 429    : omega_tangent_angles  (empty list for slab-polar)
  line 469    : rho_crossings_for_ray  (linear for slab, quadratic for curvilinear)
  line 620    : volume_kernel_mp  (Ki_1 path for cylinder-1d; exp for slab-polar/cylinder-polar/sphere-1d)
  line 653    : escape_kernel_mp  (Ki_2 for cylinder — TO REFACTOR for #116)
  line 703    : SLAB_POLAR_1D, CYLINDER_1D, CYLINDER_POLAR_1D, SPHERE_1D singletons
  line 709    : _build_volume_kernel_slab_tau_laguerre  (#115 quadrature — needs tightening)
  line 969    : _build_volume_kernel_cylinder_phi  (#116 main deliverable)
  line 1128   : build_volume_kernel dispatch
  line 1368   : compute_G_bc  (Ki_1 in cylinder branch — TO REFACTOR for #116)

orpheus/derivations/peierls_reference.py
  line 59     : imports ki_n_mp (used by slab_polar_K_vol_element — keep)
  line 171    : slab_polar_K_vol_element  (adaptive reference; proves polar ≡ E_1)

orpheus/derivations/peierls_slab.py
  line 105    : _basis_kernel_weights (legacy E_1 Nyström — TO RETIRE for #115)
  line 168    : _build_kernel_matrix (wrapper — TO RETIRE for #115)

orpheus/derivations/peierls_cylinder.py
  line 60     : GEOMETRY = CYLINDER_POLAR_1D  (production already migrated)

orpheus/derivations/_kernels.py
  line 225-290 : ki_n, ki_n_mp, ki_n_float, _ki_integrand  (TO RETIRE after #116)

tests/derivations/test_peierls_reference.py
  TestSlabPolarReferenceEquivalence        L1  polar ≡ E_1 at 1e-10
  TestSlabPolarBuildVolumeKernel           L1  production converges algebraically
  TestCylinderPolarEquivalence             L1  φ-GL converges to Ki_1 spectrally
```

---

## 3. Issue #115 — tighten slab-polar to retire legacy E₁ Nyström

### 3.1 Diagnostic findings (from issue #115 comment)

The current `_build_volume_kernel_slab_tau_laguerre` plateaus at ~9e-4 relative error for K[0,0]. The **inner τ-GL converges to machine precision** at n_rho=8 (verified against `mpmath.quad` at many μ values). The limit is in the **outer μ integration**.

**Root cause.** After the substitution `v = -ln|μ|`, the outer integrand is:

```
f(e^{-v}) · e^{-v}
```

where `f(μ)` has the Laplace expansion `L_j(x)/Σ_t + μ L_j'(x)/Σ_t² + μ² L_j''(x)/Σ_t³ + μ³ L_j'''(x)/Σ_t⁴` for small μ (polynomial in μ of degree p-1). Multiplying by the Jacobian `e^{-v}` gives a sum of `e^{-kv}` terms for k = 1, 2, 3, 4.

**Gauss-Laguerre integrates polynomials in v exactly**, not exponentials in v. Hence algebraic convergence. Same behavior for Gauss-Legendre on `[0, V_max]`.

Tests showed:

| Scheme | max rel err (K[0,0], n=128) |
|--------|------------------------------|
| v-Laguerre | 8.9e-4 |
| v-GL on [0, 30] | 8.9e-4 |
| standard GL on [-1, 0] ∪ [0, 1], n=32 | 6.2e-6 (lucky — nodes avoided a kink) |
| standard GL on [-1, 0] ∪ [0, 1], n=64 | 7.1e-4 (nodes hit kink) |

The non-monotonic behavior of standard GL is the smoking gun — there's an outer-μ kink that the node placement straddles differently at each n. Not identified analytically yet.

### 3.2 Path forward

**Recommended: adopt `scipy.integrate.quad_vec`** for the outer μ integration with a vector-valued integrand (one component per j). This:

- Is adaptive, so handles the kink/stiffness automatically
- Returns all N values per observer in one pass (row-parallel — not N² scalar calls)
- Gives float precision (~1e-14) which is far better than 1e-8
- Is the "minimal additional code" route; stays within the unified framework

**Implementation sketch:**

```python
def _build_volume_kernel_slab_adaptive(
    geometry, r_nodes, panel_bounds, radii, sig_t, dps=30,
):
    from scipy.integrate import quad_vec
    N = len(r_nodes)
    K = np.zeros((N, N))
    L = float(radii[-1])

    for i in range(N):
        x_i = float(r_nodes[i])
        sig_t_i = float(sig_t[geometry.which_annulus(x_i, radii)])

        def integrand(mu):
            # Returns an N-vector: integral over ρ of exp(-τ)·L(x+ρμ) dρ
            # evaluated via τ-GL with panel-crossing subdivision
            return _ray_vector_integrand(mu, x_i, r_nodes, panel_bounds, radii, sig_t)

        # Subdivide at μ=0; scipy.integrate.quad_vec handles adaptive refinement
        row, _err = quad_vec(integrand, -1.0, 0.0, epsabs=1e-13, epsrel=1e-11)
        row_plus, _err = quad_vec(integrand, 0.0, 1.0, epsabs=1e-13, epsrel=1e-11)
        K[i, :] = 0.5 * sig_t_i * (row + row_plus)

    return K
```

The `_ray_vector_integrand` internal is the **existing** inner τ-GL (with panel crossings + τ-cap) that already converges to machine precision — just refactored to return an N-vector of L_j values instead of a scalar.

### 3.3 Steps to close #115

1. **Refactor inner τ-GL to vector output.** Extract the ρ-integration body of `_build_volume_kernel_slab_tau_laguerre` into a helper that, given `(x_i, μ)`, returns an N-vector of `∑_τ (exp(-τ)/Σ_t(r')) · τ_w · L_j(x'(τ))`. Verify against the existing scalar path.

2. **Add `_build_volume_kernel_slab_adaptive`.** Uses `scipy.integrate.quad_vec` with μ=0 breakpoint; calls the vector inner. Dispatch based on a new kwarg or a new `kind="slab-polar-adaptive"` — whichever is cleaner.

3. **Equivalence test.** L1 test: adaptive slab-polar K matches classical E₁ K (`slab_K_vol_element`) at 1e-10 relative error for all entries of a 2-panel × p=4 example.

4. **Migrate production.** Route `peierls_slab.solve_peierls_eigenvalue` through `build_volume_kernel(SLAB_POLAR_1D, ..., adaptive=True)`. Retire `_build_kernel_matrix`, `_basis_kernel_weights`, and `_product_log_weights` in `peierls_slab.py`. Keep `solve_peierls_eigenvalue` as the user-facing entry point (so callers don't break).

5. **Run the full `test_peierls_convergence.py` suite** — these are the 5 L0 self-convergence tests. They use `peierls_slab.solve_peierls_eigenvalue` directly. Must continue to pass at the same tolerances (`1e-10` relative between successive refinements).

6. **Delete `peierls_slab._build_kernel_matrix` (and the `_basis_kernel_weights`, `_product_log_weights` helpers) if nothing in `orpheus/` or `tests/` calls them anymore.** Keep only the solver driver.

7. **Close #115** with a summary comment.

### 3.4 Risk / unknowns

- **Performance.** `scipy.integrate.quad_vec` might be slower than the fixed-order τ-Laguerre for simple observers. Mitigation: keep both paths, default adaptive for precision-critical, fixed-order for bulk assembly. BUT: if this complicates the API, just go adaptive-only and accept any perf hit (legacy E₁ was also slow-ish).
- **Adaptive quadrature tolerance vs. test gates.** The 1e-14 `epsrel` will be easy; the gate for `test_peierls_convergence.py::test_nystrom_eigenvalue_converges_1g` requires `1e-10` consecutive-refinement agreement. Should be comfortable.
- **The outer-μ kink we didn't identify analytically.** Adaptive quadrature handles it automatically — but if debugging, the non-monotonic GL convergence is the diagnostic.

---

## 4. Issue #116 — complete Bickley retirement

### 4.1 What's done

Commit 213278f: `_build_volume_kernel_cylinder_phi` replaces `Ki₁(τ)` with explicit `∫ exp(-τ/cos φ) dφ` for the volume kernel. `peierls_cylinder.GEOMETRY = CYLINDER_POLAR_1D`. Three new L1 equivalence tests.

### 4.2 What's left

Four remaining `ki_n_mp` / `ki_n_float` call sites:

| File | Line | Call | Role |
|------|------|------|------|
| `peierls_geometry.py` | 624-625 | `ki_n_mp(1, τ, dps)` / `ki_n_float(1, τ)` | `volume_kernel_mp` for **cylinder-1d** legacy path |
| `peierls_geometry.py` | 652 | `ki_n_mp(2, τ, dps)` | `escape_kernel_mp` for both cylinder kinds |
| `peierls_geometry.py` | 1459 | `ki_n_mp(1, τ, dps)` | `compute_G_bc` cylinder branch |
| `peierls_geometry.py` | 1717 | `ki_n_mp(1, τ, dps)` | Another `compute_G_bc_mode` (rank-N BC) |
| `cp_geometry.py` | 113 | `ki_n_mp(3, τ, _KI3_DPS)` | Chebyshev interpolant for **flat-source CP** (different module, out of scope for #116) |

### 4.3 Steps to close #116

1. **Add an explicit φ-integration Ki_n evaluator.** Name: `_ki_n_via_phi(n, tau, n_phi=16, cos_phi_pts=None, phi_wts=None)` — computes `∫_0^{π/2} cos^{n-1}(φ) exp(-τ/cos φ) dφ` via precomputed Gauss-Legendre nodes. Optionally accepts precomputed nodes so callers with hot loops don't re-build the φ grid every call.

2. **Refactor `CurvilinearGeometry.escape_kernel_mp`.** For `cylinder-polar` (AND cylinder-1d if we want uniformity), return `_ki_n_via_phi(2, tau)`. Verify against `ki_n_mp(2, ...)` at τ ∈ {0.1, 1, 5, 10} — spectral convergence at n_phi=32.

3. **Refactor `compute_G_bc` cylinder branch.** Replace the inner `ki_n_mp(1, tau, dps)` call with `_ki_n_via_phi(1, tau, n_phi=16, cos_phi_pts=..., phi_wts=...)`. The outer `phi` loop (surface quadrature) remains; only the kernel evaluation inside changes from one `ki_n_mp` call to a vector dot product against precomputed exp values. Add `n_phi_out_of_plane` kwarg; default 16 is fine based on #116 volume-kernel benchmarks.

4. **Refactor `compute_G_bc_mode`** (line 1717) analogously.

5. **Retire `cylinder-1d` kind.** Now that all cylinder code paths can be served by `cylinder-polar`, delete the `cylinder-1d` branches in `volume_kernel_mp`, `escape_kernel_mp`, `radial_volume_weight`, `surface_area_per_z`, `rank1_surface_divisor`, `omega_tangent_angles`, `rho_crossings_for_ray`, `optical_depth_along_ray`, `rho_max`, `source_position`, `ray_direction_cosine`, `angular_weight`, `angular_range`, `prefactor`, `d`, `S_d`, `__post_init__`. Delete `CYLINDER_1D` singleton.

   Alternative (less invasive): keep `cylinder-1d` as an alias for `cylinder-polar`. Decide on elegance vs. back-compat — the user's explicit preference will likely be "retire everything legacy".

6. **Delete `ki_n_float` and `_ki_integrand`** in `_kernels.py`. Keep `ki_n_mp` and `ki_n` only because `cp_geometry.py` still needs `ki_n_mp(3, ...)` for the flat-source CP Chebyshev interpolant (separate module, separate scope — document this).

7. **Remove `ki_n_float` import** from `peierls_geometry.py` (line 59).

8. **Test suite.** Run:
   - `tests/derivations/test_peierls_cylinder_*` — all should pass unchanged
   - `tests/derivations/test_peierls_rank_n_*` — uses `compute_G_bc_mode`, so regression coverage
   - `tests/derivations/test_peierls_reference.py::TestCylinderPolarEquivalence` — confirms Ki_1 ≡ φ-GL still
   - Add a new `TestCylinderPolarFullBickleyFree` that grep-asserts no `ki_n_mp` / `ki_n_float` call fires during a representative cylinder solve (use monkey-patching or a flag).

9. **Close #116** with a summary comment.

### 4.4 Risk / unknowns

- **Performance of `compute_G_bc`.** Currently called once per observer per BC closure. Replacing Ki₁ (one call) with a 16-node φ-GL (16 exp calls) is ~16× slower per call — but `compute_G_bc` is not in a hot loop, and our timing on the volume-kernel refactor showed the actual wall-clock impact is much less than the loop-factor suggests (16 cheap exp calls vs. one adaptive mpmath.quad).
- **`cp_geometry.py` Ki₃ staying around.** This keeps `ki_n_mp` alive in the codebase. Fine — different module, different ladder. Document in the commit message.
- **Back-compat for `CYLINDER_1D`.** If any external caller uses this singleton, deletion breaks them. Safest: keep as an alias (`CYLINDER_1D = CYLINDER_POLAR_1D`) for one cycle, with a deprecation comment.

---

## 5. Commit plan

**Commit 1: Close #115 (slab-polar adaptive + retire legacy E₁ Nyström)**

- Files: `peierls_geometry.py`, `peierls_slab.py`, `test_peierls_reference.py`
- Message prefix: `feat(derivations): adaptive slab-polar K assembly + retire legacy E_1 Nyström`
- `Closes #115` trailer

**Commit 2: Close #116 part 1 (Ki_n in G_bc and escape_kernel → explicit φ integration)**

- Files: `peierls_geometry.py`, possibly `_kernels.py` for helper
- Message: `refactor(derivations): retire Ki_1/Ki_2 in G_bc and escape_kernel via explicit φ-GL`
- Tests: new `TestBickleyRetirement` class confirming ki_n_mp unused during a cylinder-polar solve

**Commit 3: Close #116 part 2 (delete `ki_n_float`, retire `cylinder-1d` kind, clean up imports)**

- Files: `_kernels.py`, `peierls_geometry.py`, possibly `peierls_cylinder.py` if the facade needs updating
- Message: `refactor(derivations): delete ki_n_float; retire cylinder-1d kind`
- `Closes #116` trailer

**Merge to `feature/mms-broad`** after tests pass. Then the whole branch is ready to merge to `main` — `Closes #113` and `Closes #114` from the initial bug fixes will fire simultaneously.

---

## 6. Reading order for the next session

1. `.claude/lessons.md` (unconditional per CLAUDE.md)
2. Skim `.claude/plans/post-peierls-solver-bugs.md` (the original investigation plan — tells you what the bugs WERE)
3. Skim `.claude/plans/post-cp-topology-and-coordinate-transforms.md` §0.2, §4, §5, §6 (the architectural reasoning for the unification you're completing)
4. This plan (the issue-specific steps)
5. GitHub issue #115 comment (full diagnostic for why slab-polar plateaus)
6. GitHub issue #116 comment (full plan for completing Bickley retirement)
7. `git log --oneline investigate/peierls-solver-bugs ^main` — the 4 commits
8. Branch commits in order: 9ac41de → aa98565 → 4395cb8 → 213278f

Then start with Section 3 (issue #115), then Section 4 (issue #116), then Section 5 (merge).

---

## 7. Pitfalls / gotchas discovered during this branch

1. **`lagrange_basis_on_panels` clamps** points outside `[0, R]` to the nearest panel and extrapolates — NOT zero. For slab ray integrations where x' can exceed [0, L], **explicitly guard with `if x_prime < 0 or x_prime > L: continue`**. I hit this in the τ-Laguerre slab assembly.

2. **`panel_bounds` endpoints may be mpmath.mpf** from test helpers. My cylinder-polar crashed with `Cannot cast ufunc 'add' output from dtype('O') to dtype('float64')` until I added an explicit float-conversion of the panel_bounds tuples at the top of the assembly function. Same pattern applies to any new assembly function — cast panel_bounds to Python floats on entry.

3. **Ki_n_mp is adaptive, not polynomial-order.** Each call does its own internal `mpmath.quad` (~50 internal evaluations at dps=25). Replacing one ki_n_mp call with a 16-node φ-GL is actually SLIGHTLY FASTER, not slower — the factor-16 naive multiplier is misleading.

4. **Gauss-Laguerre is spectral only for polynomial-in-v integrands.** `e^{-kv}` for k ≠ 0 is not polynomial; Gauss-Laguerre converges algebraically there. This is the root cause of the slab-polar plateau. `scipy.integrate.quad_vec` sidesteps it with adaptive subdivision.

5. **Test speeds.** The full `tests/derivations/` suite takes ~14 minutes with cylinder-polar (up from ~10 min with cylinder-1d's `ki_n_float`). Runs that hit the 10-min timeout need to be split by subsuite. Budget accordingly.

6. **The `Closes #NN` trailer only fires on merge to main.** Issues #113 and #114 will show as OPEN until then. Don't hand-close them; the trailer is more auditable.

7. **`n_phi` collision.** `peierls_cylinder.py` uses `n_phi` for `n_surf_quad` (BC surface quadrature). My new parameter in `build_volume_kernel` happens to also be called `n_phi`. Keep them separate — the BC `n_phi` is not the volume-kernel `n_phi`. When refactoring `compute_G_bc` (Section 4.3 step 3), use `n_phi_out_of_plane` as the NEW parameter name to disambiguate.

---

## 8. Acceptance criteria

The work is done when:

- [ ] `peierls_slab._build_kernel_matrix` and its helpers `_basis_kernel_weights`, `_product_log_weights` are deleted (or reduced to re-exports from the unified path)
- [ ] `peierls_slab.solve_peierls_eigenvalue` routes through `build_volume_kernel(SLAB_POLAR_1D, ...)` and produces k_eff agreement with pre-refactor values to 1e-10
- [ ] `ki_n_float` and `_ki_integrand` are deleted from `_kernels.py`
- [ ] `ki_n_mp` is imported only by `cp_geometry.py` (for Ki_3) and `peierls_reference.py` (for `slab_K_vol_element`'s classical E₁ reference, which is a verification-only artifact)
- [ ] `CurvilinearGeometry` has no `cylinder-1d` branches in its methods (deleted or aliased)
- [ ] Full derivations test suite passes (expect ~14 min wall-clock with batched runs)
- [ ] `python -m tests._harness.audit` reports no new orphan equations or failed ERR catches
- [ ] `sphinx-build docs docs/_build/html` clean
- [ ] Three new commits on `investigate/peierls-solver-bugs`; merge to `feature/mms-broad` ready
