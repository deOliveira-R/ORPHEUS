---
name: L131a single-surface slab closed-form pattern
description: Issue #131 fix applied to per-face compute_P_esc_outer/inner but missed the single-surface compute_P_esc/compute_G_bc aggregate; same µ-independent τ → ½E_2 identity applies
type: project
---

**Fact.** Issue #131 fixed `compute_P_esc_outer/_inner` and `compute_G_bc_outer/_inner` to use closed-form ½E_2(τ_total) for multi-region slab (commit 3b0b2c9), but the **legacy single-surface aggregate** `compute_P_esc` / `compute_G_bc` kept its finite-N GL fallthrough at `len(radii) > 1`. Both must-fix sites are in `orpheus/derivations/peierls_geometry.py` around lines 1350 and 1448.

**Why:** L131a audit (2026-04-23) found this: `compute_P_esc`/`compute_G_bc` are called from `build_white_bc_correction` and from `build_closure_operator(reflection="marshak")` mode-0 for legacy rank-1 single-surface layout. The shipped `peierls_slab_2eg_2rg` reference uses `white_f4` which bypasses them, but any user selecting `boundary="white_rank1_mark"` on multi-region slab hits the same 4e-3 to 6e-5 GL-quadrature artefact Issue #131 patched elsewhere. Fix: replace both multi-region branches with `½·(E_2(τ_inner) + E_2(τ_outer))` / `2·(...)` using `_slab_tau_to_{inner,outer}_face` — identical closed form to per-face primitives.

**How to apply:** Whenever a slab angular integral has form `∫₀¹ f(µ)·exp(-τ/µ) dµ` with τ µ-independent (piecewise-constant σ_t), it IS closed-form as E_n regardless of n_regions. Pattern applies ONLY to slab-polar (curvilinear has τ(µ) through ρ_max). Audit other `if kind == "slab-polar": ... if len(radii) == 1: ...` guards for the same fall-through bug.

**Scope of remaining finite-N GL:** All other curvilinear angular integrals in peierls_geometry.py (lines 1378, 1476, 1504, 1528, 1675, 1763, 1838, 1918, 1997, 2021, 2096, 2120, 2188, 2265, 2516, 2582, 2751, 2839, 2867 etc) are IRREDUCIBLE by closed form because τ = σ_t · ρ_max(r_i, cos_Ω, R) depends on Ω (chord through disk/sphere). `build_volume_kernel` finite-N path explicitly exempted (integrand τ(µ) non-elementary). W transmission integrals use `mpmath.quad` adaptive — correct by design.
