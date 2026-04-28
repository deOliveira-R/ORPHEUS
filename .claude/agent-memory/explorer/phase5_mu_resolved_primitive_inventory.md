---
name: Phase 5 µ-resolved primitive inventory (Peierls specular MB reformulation)
description: Map of existing µ-resolved vs µ-integrated primitives in orpheus/derivations/peierls_geometry.py for the Phase 5 continuous-µ specular multibounce closure. Identifies what to reuse, wrap, or build.
type: project
---

# Phase 5 µ-resolved primitive inventory

Audit of `orpheus/derivations/peierls_geometry.py` (and `_kernels.py`) for the
Phase 5 reformulation:

```
K_bc^mb(r_i, r_j) = 2 ∫_0^1 G_in(r_i, µ) · F_out(r_j, µ) ·
                            [µ / (1 − e^{−σ·2Rµ})] dµ          (sphere)
```

with `µ` carried as a **rank-1 axis** through Gauss-Legendre quadrature on
[0, 1] (cylinder uses `α ∈ [0, π/2]` with `µ_2D` at the surface; slab uses
`µ ∈ (0, 1]`).

Time spent: ~30 min walking the call graph + reading the four critical
functions.

## µ-resolved primitive inventory

### Already exists, reusable as-is for Phase 5

| Primitive | File:line | Signature | Phase 5 role |
|---|---|---|---|
| `CurvilinearGeometry.optical_depth_along_ray` | `orpheus/derivations/peierls_geometry.py:550` | `(r_obs, cos_omega, rho, radii, sig_t) -> float` | **Heart of τ(µ).** Walks shell crossings; supports slab / cyl / sph; cavity-aware (hollow). For sphere antipodal chord at exit angle µ, call with `(R, -µ, 2Rµ, radii, sig_t)` (or with `r_obs=0, cos_omega=anything, rho=2Rµ` since at center `r_obs=0` only the chord length matters). For an arbitrary observer with `µ_3D` it returns the multi-region τ along the ray of length ρ. |
| `CurvilinearGeometry.rho_max` | `peierls_geometry.py:472` | `(r_obs, cos_omega, R) -> float` | Forward distance to outer surface. Used to build full chord ρ_max from any observer `r_i` at any direction `µ`. Polymorphic (slab/cyl/sph). |
| `CurvilinearGeometry.rho_inner_intersections` | `peierls_geometry.py:493` | `(r_obs, cos_omega) -> (rho_minus, rho_plus)` | Hollow-cell inner-shell intersections — needed for hollow Phase 5 to stop the chord at the cavity entry. |
| `CurvilinearGeometry.which_annulus` | `peierls_geometry.py:658` | `(r, radii) -> int` | Region index for piecewise-σ_t lookup. |
| `CurvilinearGeometry.escape_kernel_mp` | `peierls_geometry.py:806` | `(tau, dps=25) -> float` | `K_esc(τ)`: `Ki_2(τ)` for cyl, `exp(-τ)` for sphere. **NOT what Phase 5 sphere wants** — Phase 5 sphere wants raw `exp(-τ)` at a single µ (no polar pre-integration). For sphere this primitive IS `exp(-τ)` so it works directly; for cyl it bakes the polar integration in (see "premature µ-integration" below). |
| `CurvilinearGeometry.radial_volume_weight` | `peierls_geometry.py:832` | `(r) -> float` | `r^{d-1}` for `r dr` (cyl) or `r² dr` (sph) integration of source emission. Phase 5 uses this on `F_out(r_j, µ)`. |
| `CurvilinearGeometry.rank1_surface_divisor` | `peierls_geometry.py:860` | `(R) -> float` | Surface-area divisor (`R` for cyl, `R²` for sph). Needed at the K_bc assembly step regardless of µ-form. |
| `_slab_tau_to_outer_face` / `_slab_tau_to_inner_face` | `peierls_geometry.py:2432, 2456` | `(x_i, radii, sig_t) -> float` | Multi-region slab optical depth from observer to the chosen face. **µ-independent** — slab τ(µ) factors as `τ_face / µ` always. Phase 5 slab is closed form via `E_n` sums (already done in `closure="specular"` slab branch lines 5435-5475!). |
| `_slab_E_n` | `peierls_geometry.py:2418` | `(n, tau) -> float` | Closed-form `E_n(τ)` via mpmath. Phase 5 slab needs `E_n(τ_face)` evaluations only (no µ quadrature). |
| `compute_P_ss_sphere` / `compute_P_ss_cylinder` | `peierls_geometry.py:1921, 1826` | `(radii, sig_t, n_quad, dps) -> float` | Already structured as a `µ` (sphere) / `α` (cyl) Gauss-Legendre integral with `optical_depth_along_ray` for the antipodal chord — **the exact integrand structure Phase 5 needs**, just integrated to a scalar instead of being kept as a function of µ. The chord-construction loop (lines 1898-1916, 1998-2020) is the template for Phase 5 τ(µ). |
| `_shifted_legendre_eval` | `_kernels.py:383` | `(n, mu_array) -> array` | Evaluates `P̃_n(µ)` — *not needed* for the rank-1-axis Phase 5 form (the µ basis is the quadrature node itself, no Legendre projection), but useful if a test wants to project Phase 5 K_bc onto the rank-N basis to compare with the matrix form. |
| `_kernels.ki_n_float` / `ki_n_mp` | `_kernels.py:266, 274` | `(n, x) -> float` | Bickley-Naylor Ki_n. Cyl-only Phase 5 uses these in the polar-pre-integrated form (see Sphere/Cyl/Slab parity below). |

### Already exists, needs thin wrapper

| Primitive | File:line | What's missing | Wrapper sketch |
|---|---|---|---|
| `compute_T_specular_sphere` | `peierls_geometry.py:2032` | Already builds `tau(µ)` via the antipodal chord loop (lines 2126-2144) and `decay = exp(-τ_arr)` — just exposes `T_{mn}` instead of the µ-resolved quantities. **The first 25 lines of the function body ARE the τ(µ) builder Phase 5 needs.** | Refactor extraction: `def chord_tau_mu_sphere(radii, sig_t, mu_pts) -> tau_arr` returns the per-µ optical depth so Phase 5 (and the existing T-matrix build) share one implementation. |
| `compute_T_specular_cylinder_3d` | `peierls_geometry.py:2260` | Same: lines 2336-2370 already compute `tau_2d(α)` and `Ki_{j+3}(τ)` — these ARE the Phase 5 cylinder µ-resolved kernel (the `µ_2D` axis is `cos α`). The `T_{mn}` is the Galerkin projection that Phase 5 wants to skip. | Extract `def chord_tau_alpha_cyl(radii, sig_t, alpha_pts) -> (tau_arr, cos_a, sin_a)` returning the per-α optical depth + chord direction cosines. |
| `compute_T_specular_slab` | `peierls_geometry.py:2158` | The off-diagonal-only `T_oi` already uses `decay = exp(-τ_total/µ)` (line 2244). For Phase 5 slab the closed-form `E_n(τ_face)` route is preferred (already used by `closure="specular_multibounce"` slab branch lines 5252-5290). | No wrapper needed — the slab branch shows the rank-1-axis closed-form path is what Phase 5 wants for slab. |
| `compute_P_esc_outer_mode` (sphere) | `peierls_geometry.py:2638` | The mode-n integrand at lines 2699-2720 is exactly `omega_wts · angular_factor · P̃_n(µ_exit) · K_esc(τ)` — drop the `P̃_n` and you have the µ-resolved `F_out(r_j, µ_exit)`. | `def F_out_per_mu_sphere(geometry, r_nodes, radii, sig_t, mu_pts, dps) -> array shape (N_r, N_mu)`: for each observer `r_j` and each µ, return `K_esc(τ(r_j → boundary at µ_exit=µ)) · jacobian` *without* projecting onto `P̃_n`. The same loop body, with the boundary-direction parameter `µ` replacing `cos_om`. |
| `compute_G_bc_outer_mode` (sphere) | `peierls_geometry.py:3001` | Same structure: lines 3043-3071 evaluate `2 ∫_0^π sin θ · P̃_n(µ_s) · exp(-τ) dθ`. Drop `P̃_n` → `G_in(r_i, µ)`. | `def G_in_per_mu_sphere(geometry, r_nodes, radii, sig_t, mu_pts, dps) -> array shape (N_r, N_mu)`: for each observer `r_i` and each inward direction µ at the surface, return `exp(-τ(r_i ← boundary at µ))`. |
| `compute_P_esc_cylinder_3d_mode` | `peierls_geometry.py:1655` | The Knyazev integrand at lines 1717-1742 already takes `cos α` and `sin α` as the µ-axis equivalent and evaluates `Ki_{k+2}(τ_2D(α))`. Cylinder Phase 5 reuses this exact kernel structure. | `def F_out_per_alpha_cyl(geom, r_nodes, radii, sig_t, alpha_pts, dps) -> (N_r, N_α) array of Ki_2(τ_2D(α))` (k=0 term only suffices — the `µ` axis IS `cos α`; higher k_p come back if the rank-N projection is wanted). |
| `compute_G_bc_cylinder_3d_mode` | `peierls_geometry.py:1747` | Mirror of above; differs only by prefactor `(4/π)` vs `(1/π)`. | Same wrapper, different prefactor. |
| `_slab_E_n` chain in `closure="specular"` slab branch | `peierls_geometry.py:5435-5474` | The 14-line block already builds `P_o[n, i]`, `P_i[n, i]`, `G_o[i, n]`, `G_i[i, n]` from `E_{k+2}(τ_face)` sums — this IS the slab Phase 5 closed-form path. The `R_face = (1/2) M^{-1}` matrix is what Phase 5 would replace with the continuous-µ form `1/(1 − e^{-σL/µ})`. | For slab, Phase 5 is **already structurally there** — just replace the rank-N basis projection in `P_o`/`P_i`/`G_o`/`G_i` with a direct µ quadrature on `E_n(τ_face)/µ^k` integrand and the `(I − T·R)^{-1}` block with `1/(1 − exp(-σL/µ))`. |

### Needs to be built

| Primitive | Phase 5 purpose | Sketch (Python pseudocode) | Comparable existing primitive |
|---|---|---|---|
| `K_bc_specular_mb_continuous_sphere` | Top-level Phase 5 sphere driver. Produces `K_bc[i, j]` directly without the rank-N basis. | ```python
def K_bc_specular_mb_sphere_continuous(geom, r_nodes, r_wts, radii, sig_t, n_quad=64, dps=25):
    R = float(radii[-1])
    nodes, w = leggauss(n_quad); mu = 0.5*(nodes+1); mu_w = 0.5*w
    # τ along antipodal chord at each µ
    tau_arr = np.array([_chord_tau_mu_sphere(radii, sig_t, mu_k) for mu_k in mu])
    # F_out(r_j, µ) and G_in(r_i, µ): N_r × n_quad arrays
    F_out = _F_out_per_mu_sphere(geom, r_nodes, radii, sig_t, mu)
    G_in  = _G_in_per_mu_sphere(geom, r_nodes, radii, sig_t, mu)
    # Multi-bounce factor µ/(1 − e^{-2σRµ}); cancels grazing singularity
    f_mb  = mu / np.where(tau_arr > 0, 1 − np.exp(-tau_arr), 2*sig_t[0]*R*mu)
    # Rank-1-µ-axis tensor contraction: K_bc[i,j] = 2 Σ_k mu_w[k] G_in[i,k] F_out[j,k] f_mb[k]
    K_bc = 2.0 * (G_in * mu_w * f_mb) @ F_out.T
    # Apply the σ_t·rv·r_wts source-side weights (already done in legacy code)
    return K_bc
``` | The existing `compute_T_specular_sphere` builds the same `tau_arr` and `decay = exp(-τ_arr)` (lines 2126-2146). The matrix-Galerkin K_bc assembly at lines 5292-5335 shows the σ_t / rv / r_wts wiring. |
| `_chord_tau_mu_sphere` (extracted helper) | Multi-region τ along the **antipodal chord** at sphere exit angle µ. The chord has length `2Rµ` and impact parameter `R√(1-µ²)`. | ```python
def _chord_tau_mu_sphere(radii, sig_t, mu):
    R = float(radii[-1]); h = R*np.sqrt(max(0,1-mu*mu))
    radii_in = np.concatenate([[0.0], radii[:-1]])
    tau = 0.0
    for k in range(len(radii)):
        if h >= radii[k]: continue
        seg_o = sqrt(max(radii[k]**2 - h**2, 0))
        seg_i = sqrt(max(radii_in[k]**2 - h**2, 0)) if h < radii_in[k] else 0
        tau += sig_t[k] * 2 * (seg_o - seg_i)
    return tau
``` | Already inlined in `compute_T_specular_sphere` lines 2128-2144 and `compute_P_ss_sphere` lines 2010-2020. **Same code, three places** — extract once. |
| `K_bc_specular_mb_continuous_cylinder` | Cyl Phase 5 driver. The `µ` axis is `cos α` in the in-plane representation; the polar integration is pre-folded into Ki_n. | Cyl is more subtle: the polar `θ_p` integration produces `Ki_{2+k}(τ_2D)` per Knyazev. Phase 5 keeps the in-plane `α` axis but replaces matrix `(I − T R)^{-1}` with a continuous factor. The natural µ-rank-1-axis is `α ∈ [0, π/2]` (or equivalently `µ_2D = cos α ∈ [0, 1]`). | The Knyazev integrand in `compute_T_specular_cylinder_3d` lines 2363-2389 keeps the α axis explicit AND provides `Ki_{j+3}` per α — this is closer to a rank-1-axis form than the sphere case. The Phase 5 cyl form is `K_bc(r_i, r_j) = (4/π) ∫₀^{π/2} Ki_2(τ_2D(α)) · cos α · ψ_in(r_i, α) · ψ_out(r_j, α) · f_mb_cyl(α) dα` with `f_mb_cyl(α) = 1/(1 − Ki_{higher}(...))`. **Phase 5 cyl needs the operator-equivalent of the multi-bounce factor in the α basis** — research-grade, not a mechanical port. |
| `K_bc_specular_mb_continuous_slab` | Slab Phase 5 driver — closed form, no quadrature. | ```python
# Per-face block, rank-1 closed form by the µ-substitution u=1/µ:
# K_bc^slab,mb(x_i, x_j) = G_o(x_i) · 2 E_3(τ_total) / (1 − 2 E_3(τ_total)·...) · P_o(x_j) + symmetric inner-face term
# At rank-1 this is exactly Hébert's slab (1−P_ss)^-1 with P_ss^slab = 2 E_3(σL); rank-N reduces algebraically to direct E_n sums.
``` | Phase 5 slab has **no overshoot** (slab chord = L/µ → ∞ at grazing kills the integrand exponentially). The shipped `closure="specular_multibounce"` slab branch (peierls_geometry.py:5252-5290) is mathematically equivalent to the closed form, just routed through the matrix formalism. Phase 5 slab can either simplify the existing branch or leave it alone (slab is the only geometry where the matrix form converges). |

## Optical-depth-along-chord plumbing

**What is exposed cleanly:**

- `CurvilinearGeometry.optical_depth_along_ray(r_obs, cos_omega, rho, radii, sig_t)` (peierls_geometry.py:550) — fully multi-region, hollow-aware, polymorphic across slab/cyl/sph. Takes any observer position, direction cosine, and ray length. **For Phase 5 sphere antipodal chord at exit angle µ**, call it with `(r_obs=0, cos_omega=anything, rho=2Rµ, radii, sig_t)` — works because the antipodal chord passes through the radial profile symmetrically. (Or equivalently call with `(r_obs=R, cos_omega=-µ, rho=2Rµ, ...)`.)
- `_slab_tau_to_outer_face` / `_slab_tau_to_inner_face` (peierls_geometry.py:2432, 2456) — slab-specific, µ-independent (slab τ factors).
- `compute_P_ss_sphere` / `compute_P_ss_cylinder` chord loops (lines 1898-1916, 2008-2020) — inline implementations of multi-region τ along the antipodal chord with explicit impact parameter.

**What is hidden (premature integration):**

- `compute_T_specular_sphere` / `compute_T_specular_cylinder_3d` build `tau_arr` (per µ / per α) **internally** but only expose the final `T_{mn}` matrix. The `tau_arr` builder is a private chunk inside the function. **Extracting it to `_chord_tau_mu_sphere` / `_chord_tau_alpha_cyl` is a 10-line refactor and would let Phase 5 share the implementation with the existing T-matrix build.** This is the only source of duplication: **the same chord-walker is inlined in 3 places** (compute_P_ss_*, compute_T_specular_*, and individual diagnostics).

- `escape_kernel_mp(τ)` for cylinder returns `Ki_2(τ)` — the polar `θ_p` integration is pre-folded. For Phase 5 cyl this is *fine* because cyl's natural µ-axis is `α` (the in-plane angle), not the 3-D polar µ. But it means Phase 5 cyl cannot share `escape_kernel_mp` directly with sphere — sphere wants raw `exp(-τ)` at a single µ, cyl wants `Ki_2(τ_2D(α))`.

**For Phase 5 multi-region τ(µ), the recommendation is:**

1. **Extract** `_chord_tau_mu_sphere(radii, sig_t, mu) -> float` from `compute_T_specular_sphere` lines 2128-2144 (and from `compute_P_ss_sphere` lines 2010-2020, which is bit-identical).
2. **Extract** `_chord_tau_alpha_cyl(radii, sig_t, alpha) -> float` from `compute_T_specular_cylinder_3d` lines 2346-2361.
3. Both are **already** multi-region-correct (they walk the impact-parameter shell-intersection geometry the same way `compute_P_ss_*` does). They do NOT use `optical_depth_along_ray` — they use the **antipodal chord shortcut** (`h = R√(1-µ²)`, `seg_outer/inner` from sphere-shell intersection). For Phase 5 the antipodal-chord form is what's wanted (the kernel is point-to-point on opposite surfaces); the general-observer τ is only needed if Phase 5 also re-derives `F_out` and `G_in` instead of reusing the legacy mode primitives — which it does, but for THOSE the existing `optical_depth_along_ray` is exactly right.

**Summary**: τ(µ) plumbing is **fully present**. The sphere antipodal-chord helper just needs to be lifted out into a public function (10-line refactor). The general-observer τ via `optical_depth_along_ray` already powers `F_out` / `G_in` per-µ wrappers without modification.

## Sphere / cyl / slab parity

### Slab — easiest, possibly already done

Slab Phase 5 is **structurally already shipped** in `closure="specular_multibounce"` slab branch (peierls_geometry.py:5252-5290). The matrix form `(I − T·R)^{-1}` for slab does NOT exhibit the high-N pathology because:
- T is purely block off-diagonal (`T_{oo} = T_{ii} = 0`) — a single transit cannot return without an inner-face reflection.
- Slab chord = L/µ diverges at grazing, so `e^{-τ_total/µ} → 0` exponentially — no grazing singularity.
- Spectral radius `ρ(T·R) ≤ 0.08` across all N at thin τ_L = 2.5; matrix form converges as N → ∞ (only geometry where this holds — see lines 2216-2227).

**Phase 5 slab is essentially a no-op.** The closed-form `E_n` per-face primitives (lines 5435-5475 for `closure="specular"` and 5252-5290 for `closure="specular_multibounce"`) are already the "rank-1-µ-axis" form modulo cosmetic refactoring. **Priority: lowest.**

### Sphere — priority test case

This is where Phase 5 buys the most. The matrix-Galerkin `(I − T·R)^{-1}` diverges at grazing µ→0 (continuous-µ operator `1/(1 − e^{-σ·2Rµ})` is singular), but the integrand `µ/(1 − e^{-σ·2Rµ}) → 1/(2σR)` is finite — and **MC ground truth confirms** the sphere is exactly k_inf at homogeneous specular (verified by `derivations/diagnostics/diag_specular_overshoot_05_mc_multibounce.py`).

Phase 5 sphere primitives needed:
1. `_chord_tau_mu_sphere` — extract from `compute_T_specular_sphere`. **5 minutes.**
2. `F_out_per_mu_sphere(r_nodes, mu_pts) -> (N_r, n_quad)` — wrap `compute_P_esc_mode` integrand without `P̃_n`. **30 minutes.**
3. `G_in_per_mu_sphere(r_nodes, mu_pts) -> (N_r, n_quad)` — wrap `compute_G_bc_mode` integrand without `P̃_n`. **30 minutes.**
4. `K_bc_specular_mb_sphere_continuous(geom, r_nodes, r_wts, radii, sig_t, n_quad=64)` — assemble. **30 minutes.**
5. Hollow-cell support: just substitute `optical_depth_along_ray` (which already handles cavity) inside (2) and (3). **No additional work.**

**Total estimate for sphere Phase 5: 1.5-2 hours.** Then test against rank-1 Hébert (must bit-match), MC ground truth, k_inf at homogeneous, and the existing `closure="specular_multibounce"` overshoot regression.

### Cylinder — hardest

The cyl ill-conditioning is **not** the grazing singularity (cos α weight kills it at α → π/2) but `R = (1/2) M^{-1}` matrix conditioning. Phase 5 cyl bypasses `R` entirely by working in the α-basis directly:

```
K_bc^cyl,mb(r_i, r_j) = (4/π) ∫_0^{π/2} cos α · ψ_in_cyl(r_i, α) · ψ_out_cyl(r_j, α) · f_mb^cyl(α) dα
```

But: what IS `f_mb^cyl(α)`? It's the cyl analog of `µ/(1 − e^{-σ·2Rµ})`. The single-transit transmission at in-plane angle α is `Ki_2(τ_2D(α))` (NOT `e^{-τ_2D}` — the polar integration converts it), so the multi-bounce factor needs to be derived from the cyl's effective 1-D operator on α. **This is research-grade, NOT a mechanical port.**

**Cylinder Phase 5 is a separate research question.** Recommendation: ship sphere + slab Phase 5 first; circle back to cyl with literature support (Sanchez & McCormick on cyl P_ss, or Carlvik on cylinder integral-transport).

## Risks / surprises found

1. **The same chord-tau builder is inlined in 3 places.** `compute_T_specular_sphere`, `compute_T_specular_cylinder_3d`, `compute_P_ss_sphere`, and `compute_P_ss_cylinder` all carry their own copy of the multi-region antipodal-chord τ loop. They are bit-identical (verified by inspection). Phase 5 should **extract these as public helpers** — every byte of duplication is an integrity risk for multi-region cells. (Foundation-test target.)

2. **The `(rho/R)²` Jacobian convention is a known minefield.** The bare-specular sphere branch (lines 5294-5328) uses the **no-Jacobian** integrand explicitly because the standard `compute_P_esc_mode` (line 3527) carries `(ρ_max/R)²` "surface-to-observer Jacobian" that is empirically necessary for rank-N Marshak but breaks rank-N specular. Phase 5 sphere must use the **no-Jacobian** form (matches the SymPy derivation in `derivations/peierls_specular_bc.py`). The bare-specular code on lines 5294-5328 is the right template — NOT `compute_P_esc_mode`. The lesson: the explorer rule "the most-connected primitive is the canonical one" is FALSE here — `compute_P_esc_mode` (degree 131) is the rank-N Marshak primitive, NOT the canonical no-Jacobian form Phase 5 wants.

3. **Cylinder µ-axis is `cos α`, NOT a 3-D µ.** The polar integration is pre-folded into `Ki_2(τ_2D)` and the natural rank-1 axis is the in-plane `α` (or `cos α ∈ [0, 1]`). This is GOOD for Phase 5 — it means the high-N pathology there comes from `R = (1/2) M^{-1}` conditioning (not from a singular continuous-µ kernel), and Phase 5's continuous-α form has a bounded resolvent. But **the multi-bounce factor in α-space is NOT the sphere's `µ/(1−e^{-σ·2Rµ})`** — it has to be re-derived. This is a research-blocker for cyl Phase 5.

4. **Hollow-cell support is free for sphere, not free for cyl.** `optical_depth_along_ray` already handles cavity zero-σ_t segments (lines 625-632, 645-647) and `rho_inner_intersections` exposes the cavity entry. For sphere Phase 5, the existing `compute_P_esc_outer_mode` and `compute_P_esc_inner_mode` (peierls_geometry.py:2638, 2725) ONLY support hollow sphere (raise `NotImplementedError` for cyl/slab). Phase 5 sphere hollow is therefore "free" — just call them at mode 0 with µ-dependent boundary parameters. Phase 5 cyl hollow needs new work.

5. **`closure="specular_multibounce"` already emits a `UserWarning` at N ≥ 4 for sphere/cyl.** This is the entry point for Phase 5: when Phase 5 ships, this warning should be replaced by automatic dispatch to the continuous-µ form for sphere (and slab can stay on the matrix form since it has no pathology, or also dispatch for uniformity).

6. **The `derivations/diagnostics/diag_specular_overshoot_11_continuous_kernel.py` file is a stub** (`raise NotImplementedError` at line 85) but the docstring at lines 1-85 is the most thorough Phase 5 design sketch in the repo — it explicitly identifies `OUT(r_j, µ) · IN(r_i, µ) · MULT(µ)` as the rank-1-µ-axis tensor form and notes that `MULT(µ) = 1/(1−e^{-σ·2Rµ})` is the multi-bounce factor. **Phase 5 implementation should resurrect that diagnostic as the first regression test** (it was never finished).

7. **MC ground truth exists.** `diag_specular_overshoot_05_mc_multibounce.mc_specular_sphere` (line 21) is a working MC reference for homogeneous specular sphere. Phase 5 must regress against it (`test_mc_specular_sphere_thin`, `test_mc_specular_sphere_very_thin`).

8. **The bare-specular slab branch and the multi-bounce slab branch share 90% of their code** (peierls_geometry.py:5435-5475 vs 5252-5290) — only the `T_slab` term differs. Phase 5 can refactor these two into one with a `multi_bounce: bool` flag and shrink the file by ~40 lines.
