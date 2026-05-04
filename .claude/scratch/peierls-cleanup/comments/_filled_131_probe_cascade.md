<!-- Posted as part of the post-#138 Peierls documentation cleanup, 2026-04-30 -->
<!-- Source: docs/theory/peierls_nystrom.rst lines 803-859 (~57 LoC) -->
<!-- Cleanup commit: [`18a852b`](https://github.com/deOliveira-R/ORPHEUS/commit/18a852b) (cleanup tip; full range `742d3b0..18a852b`, 25 commits) -->

## Issue #131 close-out — Probe A/B/D diagnostic cascade and the closed-form $E_2$ fix

**Summary (2026-04-30):**

1. The original Phase G.5 benchmark (commit `aa6ebf0`) showed a **1.5 % $k_{\rm eff}$ disagreement** in the multi-region slab path, isolated to MR (1G/2G single-region parity tests already passed at $10^{-8}$).
2. Probes A (1G 2-region), B (2G 2-region) and D (multi-region $P_{\rm esc}$ quadrature) cascade-isolated the bug to the **finite-N Gauss–Legendre branch** of `compute_P_esc_outer/inner` and `compute_G_bc_outer/inner`, which converged only to $\sim 4\times 10^{-3}$ at $N=24$ — a quadrature artifact of using GL where a closed form was available.
3. The fix routes **all slab-polar calls (any `n_regions`) through the closed-form $\tfrac{1}{2}E_2(\tau)$ branch**, since the angular integral is closed-form regardless of region count. After the fix: rel_diff drops from 1.5 % to **5.4 × 10⁻¹⁶** (bit-exact to machine epsilon); the parity-gate test now asserts `rel_diff < 10⁻¹⁰` and flips from a diagnostic-recording test into a regression guard.

---

<details><summary>Full investigation record (relocated from docs/theory/peierls_nystrom.rst:803–859)</summary>

### How the 1.5 % gap was diagnosed (Issue #131)

The original Phase G.5 benchmark (recorded in commit `aa6ebf0`) showed a 1.5 % disagreement on $k_{\rm eff}$. Single-region 1G and single-region 2G parity tests (`TestSlabPolarVsNativeE1KEff`, `TestMGSlabPolarMatchesNativeSlabMG`) had shown 1e-8 agreement, so the gap was specific to multi-region slab.

The numerics-investigator agent (Issue #131) ran a cascade of isolation probes under `derivations/diagnostics/diag_slab_issue131_probe_*.py`:

- **Probe A** (1G 2-region vacuum) — isolated multi-region from MG / closure. Already showed the gap → ruling out MG-only bugs.
- **Probe B** (2G 2-region vacuum) — isolated from F.4 closure. Gap persisted → bug is in the volume kernel, not the closure.
- **Probe D** (multi-region $P_{\rm esc}$ quadrature) — pinned down the cause: `compute_P_esc_outer` / `compute_P_esc_inner` / `compute_G_bc_outer` / `compute_G_bc_inner` had separate branches for `len(radii) == 1` (closed-form $\tfrac{1}{2} E_2(\tau)$) and `len(radii) > 1` (finite-N GL over the µ-integral). The multi-region branches converged only to $\sim 4 \times 10^{-3}$ at $N=24$ — a quadrature artifact.

### The fix

For a slab with piecewise-constant $\Sigma_t(x)$, the angular integral

```math
P_{\rm esc}^{\rm outer}(x_i)
  = \frac{1}{2}\!\int_0^1\!
    \exp\!\bigl[-\tau_{\rm outer}(x_i)/\mu\bigr]\,d\mu
  = \frac{1}{2}\,E_2\!\bigl(\tau_{\rm outer}(x_i)\bigr),
```

is **closed-form regardless of the number of regions**, because

```math
\tau_{\rm outer}(x_i) = \sum_k \Sigma_{t,k}\,(r_k - \max(r_{k-1}, x_i))^+
```

is independent of µ. The GL quadrature branch was therefore wasteful (and underconvergent). The fix adds two helpers `_slab_tau_to_outer_face` and `_slab_tau_to_inner_face` in `orpheus.derivations.peierls_geometry` that piecewise-integrate $\Sigma_t$ across region boundaries, and routes **all slab-polar calls (any `n_regions`) through the closed-form branch**.

### Result

The shipped 2eg_2rg parity drops from `rel_diff = 1.5 %` to `rel_diff = 5.4 × 10⁻¹⁶` — **bit-exact to machine epsilon**. Same for the 1-region case (already bit-exact before, untouched by the fix). The parity-gate test `tests.derivations.test_peierls_multigroup.TestSlabViaUnifiedDiscrepancyDiagnostic` now asserts `rel_diff < 10⁻¹⁰` (5 orders of margin above the current measurement) and flips from a diagnostic recording test into a regression guard.

</details>

### Production decision (2026-04-30)

- **Closed-form $E_2$ branch is canonical.** Both unified-vs-native paths now agree bit-exactly on the shipped fixture; the GL branch is gone.
- **Test posture:** `TestSlabViaUnifiedDiscrepancyDiagnostic` is now a regression guard with $10^{-10}$ tolerance.
- The Sphinx page now carries only a one-paragraph pointer to this comment (replacing the full diagnostic narrative at lines 803–859).

### Cross-links

- Surviving Sphinx anchor: `:ref:\`theory-peierls-slab-polar-g5-diagnosis\`` → now a routing stub pointing here.
- Code: `compute_P_esc_outer/inner`, `compute_G_bc_outer/inner`, `_slab_tau_to_outer_face`, `_slab_tau_to_inner_face` in `orpheus.derivations.peierls_geometry`.
- Tests: `tests.derivations.test_peierls_multigroup.TestSlabViaUnifiedDiscrepancyDiagnostic` (regression guard, $10^{-10}$).
- Diagnostic scripts: `derivations/diagnostics/diag_slab_issue131_probe_{A,B,D}.py`.
- Original 1.5 % benchmark: commit `aa6ebf0`.
- Cleanup commit: [`18a852b`](https://github.com/deOliveira-R/ORPHEUS/commit/18a852b) (cleanup tip; full range `742d3b0..18a852b`, 25 commits)
