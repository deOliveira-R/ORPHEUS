---
name: Peierls rank-N / rank-2 white-BC closure
description: Design conventions for Marshak/DP_N rank-N extension AND Phase F rank-2 (per-face) BC of the unified polar-form Peierls solver
type: project
---

Two related rank extensions of the unified polar-form Peierls solver
in `orpheus/derivations/peierls_geometry.py`:

1. **Rank-N angular** (N modes per surface) — Marshak/DP_N closure on
   shifted-Legendre basis, single surface.
2. **Rank-2 per-face** (Phase F, Issue #110) — two surfaces (outer +
   inner) for Class-A geometries (slab + hollow cyl + hollow sph),
   mode space `A = ℝ^(N_modes × N_surfaces)`. Regime A (inner_radius
   == 0.0 exactly) must be BIT-EXACT to rank-1 solid; Regime B
   (inner_radius → 0⁺) must converge at predicted algebraic rate.

**Why:** rank-1 Mark hits 21-27 % error at R=1 MFP for solid cyl/sph;
for slab + hollow it fails the Wigner-Seitz identity φ ≡ 1/Σ_t at
finite L (16-40% k_eff error). Rank-N fixes angular resolution
(Sanchez-McCormick 1982 §III.F.1 Eqs. 165-169, ~10× per rank).
Rank-2 fixes the missing inner-surface T·J⁻ self-feedback transmission
that rank-1 Mark omits.

**How to apply (rank-N):**

- Sphinx label: `peierls-rank-n-bc-closure`.
- Rank-1 bit-exact recovery (rtol=1e-14) mandatory regression gate.
- Cross-mode diagonality via finite-diff SVD: `K_bc(N) - K_bc(N-1)`
  must be rank-1 (σ_2/σ_1 < 1e-10).
- Thin-cell ladder N ∈ {1, 2, 3, 5, 8} at R=1 MFP; strict monotone
  error decrease; N=8 < 1%.
- Thick-cell invariant at R=10 MFP: |k(N) − k(1)| < 1e-3 for N ∈
  {2, 3, 5}.

**How to apply (rank-2 / Phase F):**

- Phase F.1+F.2 are ADDITIVE-only: no behavioural change to existing
  solid code paths. L0 + foundation tests only; no L1 until F.3.
- Factor-level tests at 1e-14 (slab E_2 match) are reachable via
  `mpmath.expint(2, τ)` as the independent-route check — scipy/mpmath
  special-function lib IS the dual-route for slab.
- Per-surface functions: `compute_P_esc_outer/inner`,
  `compute_G_bc_outer/inner`. Solid geometries MUST return zero-array
  (NOT None, NOT raise) from `_inner` variants — this is the
  bit-exact sum-over-surfaces contract for regime A.
- The §9.1 "dual-route discipline" (learned from commit 2538cfe
  factor-of-2 bug): for any NEW analytical formula, verify via TWO
  independent computational paths — never one. Fast GL
  implementation vs adaptive mpmath.quad of the same integrand is
  the canonical pattern for hollow-core primitives.
- Rank-1 regression gate for Phase F.1+F.2 is 5-6 existing tests:
  `test_mark_rank1_equals_legacy_build_white_bc_correction`,
  `test_marshak_matches_build_white_bc_correction_rank_n`,
  `test_apply_matches_as_matrix`,
  `TestSlabPescClosedForm::test_matches_E2_sum_at_machine_precision`,
  `TestSlabGbcClosedForm::test_matches_E2_sum_at_machine_precision`,
  `TestSlabKbcStructure::test_K_bc_row_sum_matches_first_order_closed_form`.
- Cavity-segment optical-depth bug-hiding: a Σ_t-sign flip in cavity
  handling is invisible to sum-consistency tests (hits outer and
  inner identically). Dual-route REQUIRED — hand-algebra
  reconstruction of τ from `rho_inner_intersections` + `rho_max`
  against the implementation.

**File placement (Phase F.1+F.2):**

- NEW file `tests/derivations/test_peierls_geometry.py` for
  dataclass + method-contract tests (inner_radius validation,
  rho_inner_intersections, optical_depth_along_ray cavity handling).
- Extend `tests/derivations/test_peierls_reference.py` Layer 6 with
  per-face E_2 tests (sibling to `TestSlabPescClosedForm`).
- Extend `tests/derivations/test_peierls_closure_operator.py` with
  the regime-A zero-array sentinel (next to existing rank-1
  legacy-equivalence tests).
