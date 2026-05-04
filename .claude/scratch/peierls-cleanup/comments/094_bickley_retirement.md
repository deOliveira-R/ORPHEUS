# Issue #94 close-out — BickleyTables retirement (post-#138 docs cleanup)

**Status:** CLOSED — `BickleyTables` and `bickley_tables()` were deleted from `orpheus.derivations._kernels` in Phase B.4 (commit `6badbe5`); cylinder kernel callers route through the canonical `ki_n_mp` (mpmath) and the fast Chebyshev approximant `_ki3_mp` in `orpheus.derivations.cp_geometry`.
**Date:** 2026-04-30
**Why this comment exists:** Per the post-#138 docs cleanup directive, failed-experiment / completed-rollout narrative is being relocated from `docs/theory/peierls_nystrom.rst` into the GitHub issues that originated it. This comment carries the verbatim replacement table, the three-phase retirement sequence (with commit hashes), and the four-bullet measured-impact postmortem out of `peierls_nystrom.rst:8314–8413`. The Sphinx page now keeps only the intro + "Why Chebyshev not direct mpmath" rationale.

<!-- COMMIT-HASH -->

## Replacement table — what each legacy call became

| Legacy call (pre-Phase B.4) | Canonical replacement (post-Phase B.4) | A&S identity |
|---|---|---|
| `tables.ki3(x)` / `ki3_vec(x)` | `ki_n_mp(2, x, dps=30)` | canonical $\mathrm{Ki}_2(x)$ |
| `tables.ki4(x)` / `ki4_vec(x)` | `orpheus.derivations.cp_geometry._ki3_mp` (fast) or `ki_n_mp(3, x, dps=30)` (arbitrary precision) | canonical $\mathrm{Ki}_3(x)$ |
| `tables.Ki2_vec(x)` (canonical alias, added Phase 4.2) | `ki_n_mp(2, x, dps=30)` | canonical $\mathrm{Ki}_2(x)$ |
| `tables.Ki3_vec(x)` (canonical alias, added Phase 4.2) | `orpheus.derivations.cp_geometry._ki3_mp` (fast) or `ki_n_mp(3, x, dps=30)` (arbitrary precision) | canonical $\mathrm{Ki}_3(x)$ |
| `e3(x)` / `e3_vec(x)` (slab) | `orpheus.derivations._kernels.e3_vec` (retained; already double-precision via `scipy.special.expn`) | canonical $E_3(x)$ |

## Retirement sequence — what actually happened

<details><summary>Phase B.1 / B.2 / B.4 commit chain (verbatim)</summary>

1. **Phase B.1** (commit `ea6b05e`, theory-first): `docs/theory/peierls_nystrom.rst` §§12–17 landed as a theory page before any code changed, naming the forthcoming modules and the unified $\Delta^{2}$ operator.
2. **Phase B.2** (commits `f1b869b` + `bf128d3`): the new `orpheus.derivations.cp_geometry` module was implemented with `FlatSourceCPGeometry` and the three singletons `SLAB`, `CYLINDER_1D`, `SPHERE_1D`; the pre-existing `orpheus.derivations.cp_slab`, `orpheus.derivations.cp_cylinder`, and `orpheus.derivations.cp_sphere` modules became thin facades over the geometry-dispatching core. `BickleyTables` was **no longer imported** by any `cp_*` derivation module, but the class itself was kept in `orpheus.derivations._kernels` so the Phase B.2 commit was a drop-in refactor with bit-identity to Phase A (safety milestone).
3. **Phase B.4** (commit `6badbe5`, this postmortem's subject): `BickleyTables` and `bickley_tables()` were deleted from `orpheus.derivations._kernels`. The cylinder kernel was replaced by `orpheus.derivations.cp_geometry._ki3_mp`, a Chebyshev polynomial of degree 63 fit to the scaled kernel $e^{\tau}\,\mathrm{Ki}_3(\tau)$ on $[0, 50]$ at Chebyshev-Gauss-Lobatto nodes (~$5\times 10^{-6}$ absolute accuracy; build cost ~0.3 s, lazy via `functools.lru_cache`). The runtime solver `orpheus.cp.solver` was rewired in the *same commit* to import `_ki3_mp` from `orpheus.derivations.cp_geometry` and consume it via `_setup_cylindrical`; the solver's own private `_build_ki_tables` + `_ki4_lookup` pair (~30 lines of cumsum-based $O(h)$ quadrature) were deleted. Solver `keff` and derivation `k_inf` now evaluate $\mathrm{Ki}_3$ through the **same code path** — the solver/derivation kernel-split bias that had been hiding behind the `CPParams.n_ki_table` knob is gone (the knob is retained as an unused no-op for construction-site backwards compatibility).

</details>

## Phase B.4 postmortem — measured impact

The kernel swap was an **improvement**, not a regression. The measurable shifts are:

- **Cylinder** $k_\infty$ **reference values** shifted by up to ~$4\times 10^{-4}$ for multi-region 1-group cases. The Bickley tabulation's trapezoidal $O(\Delta x^2)$ error had been the dominant bias in the reference; each new value is closer to the exact mpmath result than the pre-refactor one.
- **Solver/reference agreement.** The `solve_cp` cylinder `keff` now agrees with the shifted $k_\infty$ reference to machine precision (same kernel on both sides). All nine `cp_cyl1D_*` L1 eigenvalue tests pass at their declared `tolerance = 1e-5` with actual error ~$10^{-7}$ — about 100× headroom, where previously the `1e-5` tolerance had been the *actual* floor set by kernel bias.
- **Tabulation-size sensitivity test retired.** The old `test_cylindrical_ki4_convergence_with_table_size` (which documented that 5 000 → 20 000 → 40 000 points gave diminishing returns) was replaced by `test_ki3_kernel_is_insensitive_to_n_ki_table`: `n_ki_table` is a no-op, and `keff` is bit-identical across `{5000, 20000, 40000}`.
- **Solver startup latency.** The 20 000-point `scipy.integrate.quad` loop at `CPMesh` construction is gone; the Chebyshev polynomial is built lazily on first call to `_ki3_mp` (~0.3 s once per process) and cached via `functools.lru_cache`. Repeated solves pay zero kernel setup cost.

## Surviving Sphinx anchors

The `peierls_nystrom.rst` §16 BickleyTables-retirement section retains:

- The status note + intro (the 20 000-point lru_cache + `scipy.quad` context that motivated retirement).
- The "Why Chebyshev not direct mpmath" rationale (explains why `_ki3_mp` is a polynomial fit rather than a thin mpmath wrapper).

This comment is the canonical home for the replacement table, the retirement-sequence commit chain, and the measured-impact postmortem.

## Cross-links

- Sphinx (post-cleanup): `docs/theory/peierls_nystrom.rst` §16 "BickleyTables retirement" — intro + Chebyshev-vs-mpmath rationale.
- Cited tests: `tests/derivations/test_kernels.py::test_legacy_bt_ki3_equals_kin2`, `test_legacy_bt_ki4_approximates_kin3` (regression guards retired with the class), `test_ki3_kernel_is_insensitive_to_n_ki_table` (replacement test).
- Originating commits: `ea6b05e` (Phase B.1 theory), `f1b869b` + `bf128d3` (Phase B.2 facades), `6badbe5` (Phase B.4 deletion + Chebyshev replacement).
