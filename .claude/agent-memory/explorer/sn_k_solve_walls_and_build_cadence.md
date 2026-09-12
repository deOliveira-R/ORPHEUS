---
name: sn-k-solve-walls-and-build-cadence
description: Durable [M] facts for any SN k-path carve — no bit-exact wall exists on the forward k-solve (all eigen fixtures are SAFETY×conv_tol bands; nulp walls cover only L+C), build_within_group_system runs once per OUTER, and the two F spellings differ by 2–5 ULP
metadata:
  type: reference
---

Measured 2026-09-12 at `b0fd3e7e` for step 2 (R-cc6) of the consumers campaign; full census
`scratch/_consumers/explorer_step2_census.md` (untracked — re-derive via the probes named there).

- **No bit-exact wall on the forward k-solve.** `tests/sn/regression/snapshots/*.npz` (10 eigen
  cases, all L=0), `tests/sn/_data/finalize_reconstruction_448/*` (8 arms × L∈{0,1}) and
  `affine_carve_converged/*` all assert `assert_regression(kind="iterative")` =
  `allclose(SAFETY × conv_tol)`. The only `nulp=1`/`array_equal` walls are operator-level on
  `L + C` (`walk_matvec_*.npz`, `affine_carve_baseline/*.npy`, `pre_t4_snapshots.npz`) and the
  0-D homogeneous D5 (`cs1_prewiring.json`, byte-exact). Claiming bit-identity on the k path
  needs a CAPTURED `nulp` anchor first.
- **Build cadence today:** `build_within_group_system` (and with it `L+C`, `C=M[σ_t]`, `B_a`,
  the grids, System B's four operators) runs once per OUTER on the eigen entries (spy: 6/6/6 SI,
  3/3/3 Krylov), once per solve on fixed-source/adjoint. `S`, `N₂ₙ`, `F`(energy), `mat_xs`,
  `loss_kernel_gauge`, `CollisionCache` are once per solve in `SNSolver.__init__`. The adjoint
  and the DSA arm each mint a SECOND `MaterialXSField` (`material_xs_field()` is uncached).
- **The two `F` spellings are not bit-identical:** `FissionOperator.apply(lift φ/W)` vs
  `lift(IsotropicFission.apply φ)` — 2 ULP (slab GL8) / 5 ULP (LS4 2-D), max rel ≤ 6.9e-16, one
  random φ each; `IsotropicFission.__post_init__` REFUSES a composite end by construction.
- **Hub `cached_property` surface (13) is entirely σ-free**; a cached pencil/`loss`/`L+C`
  would be the first σ_t SNAPSHOT there and collides with `rebind_cross_sections`' documented
  contract (0 production callers, 1 test pin on cache identity, 0 solve-level pins).
- `weight_norm` on `SNSolver` has 0 reads in every root (dead); `self.quad` 0 production reads.
