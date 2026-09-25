# #405 step 2: where the time goes in the 8 slowest worklist files

HEAD `901f64ca` (clean for tracked files). Durations `[M]`: artifact `test-durations` of run 35940553034 (downloaded to `scratchpad/dur/test_durations.json`; per-family table in `scratchpad/dur/per_test.txt`). Code facts `[M]` = read at HEAD. Splits `[R]` = inferred from those durations plus two local timings (below).

## Taxonomy note (read first)

The brief's classes are G (the reference's generator builds a reference for a SUT), S (production SUT), C (the reference against a finer copy of itself). In 5 of the 8 files the thing under test IS the generator (everything under `orpheus/derivations/`). The generator's output is compared with a closed form (k_inf), a sibling closure (Mark, Hébert), a sibling generator (native E1 Nyström), or pinned literals. I label these **C\*: the generator under test**. That is the plan's row 3 ("a test of the generator is keyed on the generator"), which is wider than strict self-convergence. Where the comparison really is the generator against itself at two resolutions, I write **C (strict)**.

## Local timings `[M]` (macOS arm64, mpmath 1.3.0 python backend, numpy 2.4.4; the runner has numpy 2.5.3)

- `build_volume_kernel(SLAB_POLAR_1D, …)` at 1 region, 2 panels × p_order 4 (8 nodes), n_angular 24, n_rho 24, dps 20, fuel-A 1G, R=5: **90.0 s**, returns `(8,8) float64`. About 1.4 s per element (adaptive `mpmath.quad` per element; the slab path ignores n_angular and n_rho, per the code comment at geometry.py "mpmath.quad self-determines node count").
- Native `slab.solve_peierls_eigenvalue` 1G 2-region vacuum at 2 panels × p3, dps 20: **1.4 s**, returns `PeierlsSlabSolution`.

## Generators (entry points, payloads, determinism)

**GEN-P, the Peierls Nyström solve.** `orpheus.derivations.continuous.peierls_nystrom.geometry:solve_peierls_mg` (and `solve_peierls_1g`, which wraps it).
- Inputs: `CurvilinearGeometry` (`SLAB_POLAR_1D`, `CYLINDER_1D`, `SPHERE_1D`, or hollow with `inner_radius`), radii `(n_reg,)`, sig_t `(n_reg,ng)`, sig_s `(n_reg,ng,ng)`, nu_sig_f, chi, boundary (closure name), n_bc_modes, n_panels_per_region, p_order, n_angular, n_rho, n_surf_quad, dps, tol, max_iter.
- Returns `PeierlsSolution`, a `@dataclass(frozen=True)`: `r_nodes (N,)`, `phi_values (N,ng)`, `k_eff: float`, ints, and `panel_bounds: list[tuple]`. This is **plain data**. `.phi()` is a method that rebuilds a Lagrange basis from those fields.
- Inner structure `[M]`: `_build_full_K_per_group` calls `build_volume_kernel(...)` on every solve and for every group, with no memo. K_vol depends on (geometry, radii, sig_t_g, panels, p_order, n_angular, n_rho, dps) and NOT on boundary, n_bc_modes, sig_s, nu_sig_f or tol. So a ladder over N or over the closure rebuilds an identical K_vol each time.
- Determinism: in `orpheus/derivations/`, 0 hits for `np.random|default_rng|ThreadPool|ProcessPool|concurrent.futures|multiprocessing|joblib` (positive control: the same grep hits `orpheus/mc/solver.py` and `orpheus/numerics/manifold.py`). Precision is set by `mpmath.workdps(dps)`, which is deterministic. Platform-dependent: the numpy/LAPACK float64 assembly and eigen-iteration after the mpmath stage (the ULP-class drift of #504). mpmath's backend (python or gmpy) changes speed, not results `[R]`.

**GEN-SHIP, the shipped registry references.** `orpheus.derivations.continuous.peierls_nystrom.cases:continuous_case_builders()` gives 13 thunks (`partial(_build, case)`), and `_build` calls `build_two_surface_case`.
- Slab: `_build_peierls_slab_case_via_unified` is GEN-P, white_f4, at **16 panels × p6 (192 nodes, 2 regions), dps 30**.
- Hollow cylinder and sphere (1G and 2G; r0/R of 0.1, 0.2 and 0.3): GEN-P, white_f4, at 3 panels × p5, 24/24/24, dps 20.
- Returns `ContinuousReferenceSolution` (frozen dataclass). `phi` is a **closure** `phi_fn` over a `PeierlsSolution`, which is itself plain data, so it can be rebuilt from arrays. `problem.materials` holds `Mixture` objects built from `get_mixture`.
- **Env-dependent** `[M]`: `_SLAB_VIA_UNIFIED` is read from `ORPHEUS_SLAB_VIA_E1` at import time. A cache key must include it.

**GEN-TR, the trajectory resolvent (cylinder, multi-region).** `orpheus.derivations.continuous.trajectory_resolvent.greens_function_cylinder:solve_greens_function_cylinder_mr`.
- Inputs: radii `[0.5,1.5,2.0]`; `get_xs` A/B/A 2G; alpha 1.0; n_r 24, n_mu_axial 16, n_phi_az 32, n_traj_quad 64; max_iter 500; tol 1e-7; initial_k 1.23.
- Returns `CylinderGreensMRResult`, a frozen dataclass of plain data: `k_eff`, `psi_g (G,n_r,n_mu,n_phi)`, `phi_g (G,n_r)`, node arrays, `region_at_node`, iterations, converged.
- numpy/scipy only (no mpmath in `trajectory_resolvent/*.py`). No randomness or threads (0 hits, above). Platform-dependent at the ULP level `[R]`.
- Cost: `[M]` 1316 s on the runner (`phase_d[cyl_2g_3reg…]`, whose other side is a literal).

## Per-file classification

### 1. `tests/gates/derivations/test_peierls_specular_bc.py` (188.9 min, 27 cases) — C\*, the time is GEN-P

Every slow test calls `solve_peierls_1g` or `solve_peierls_mg` (p4, 2 panels, 24/24/24, dps 20) over a ladder of n_bc_modes and checks against closed-form k_inf, a Mark or Hébert solve, or monotonicity.

| family | runner s | class |
|---|---|---|
| `heterogeneous_2G2R_converges[slab]` | 3000 (timeout) | C\* (GEN-P ×4 N ×2 g) |
| `heterogeneous_1G2R_converges[slab]` | 2825 | C\* |
| `2G_homogeneous…kinf_2G[slab]` | 1646 | C\* (N=1..6 ×2 g) |
| `…OvershootCharacterization::slab…N20` | 898 | C\* (6 solves) |
| `slab_homogeneous_converges_to_kinf` | 564 | C\* (4 solves) |
| `multibounce_slab_monotonic_high_N` | 452 | C\* (3) |
| `slab_rank1_equals_mark_kinf` | 302 | C\* (2) |
| `multibounce_slab_rank1_lifts_plateau` | 301 | C\* (2) |
| cylinder/sphere members | 487+286+246+87+51+37+32+53+33+… | C\* |

The slab cases are 9987 s, which is 166 of the 189 minutes.

The split inside a solve `[R]`: a slab solve on the runner is 141 s (564/4), and K_vol alone is 90 s locally. So K_vol, the adaptive mpmath stage, dominates.

**Reuse** `[M]`, by reading the args:
- Whole-solve keys: 1 duplicate in the file (specular N=1 appears in both `rank1_equals_mark` and `homogeneous_converges`).
- K_vol keys (slab only): **41 builds collapse to 7 distinct keys**.
  - fuel-A 1G: 6 builds, 1 key
  - 2G homogeneous: 12 builds, 2 keys
  - AB 1G: 4 builds, 1 key
  - AB 2G: 8 builds, 2 keys
  - thin slab: 11 builds, 1 key
- The cylinder and sphere ladders have the same N-ladder structure (K_vol is rebuilt per N).

### 2. `tests/gates/sn/verification/analytical/test_l1_standoff_slab_cylinder.py` (103.3 min) — mixed

- **Cylinder** `test_cylinder_l1_refinement_both_paths[20,40,80]` (1060/1133/1382 s) and `…sweep_vs_trajectory_resolvent` (1043 s): each calls `_cylinder_k_ref()` = GEN-TR. The function has **no cache**. Its docstring's "Cached at module import … (~30 s)" is false `[M]`: the function is plain, and it costs about 1000–1300 s on the runner.
  - The SN SUT (`solve_sn`, folded_product 4×8, sweep + krylov) costs 95 s at nx=40 `[M]` (the twin-path test, which has no ref).
  - Split `[R]`: [20] is about 95% G, [40] about 90% G, [80] about 70% G (1382 − ~1000 s of ref means about 380 s of S).
  - The whole cylinder family: about 4300 of 4618 s is G.
- **Slab** (refinement [40,80,160], `krylov_vs_case`, twin): **S**. The reference is `continuous_get("sn_slab_1eg_2rg_S8")` and costs about 0 s `[M]`: refinement[80] takes 148.6 s and the no-ref twin at 80 takes 148.7 s.
- **Reuse:** GEN-TR with identical args is built **4× in this file, and 7× across 3 files**: here ×4, `test_phase_c_crosscheck.py` phase_d ×1 and phase_e ×1, and `tests/gates/sn/sweep/curvilinear/test_unified_matvec_cylinder.py:447` ×1 (1425 s on the runner). That is about 2.3 h of runner time for one reference.

### 3. `tests/gates/derivations/test_continuous_registry_lazy.py` (100.1 min, 2 timeouts) — G, and the payload is unused

- `test_builder_keys_match_built_names` calls `continuous_cases()`, which builds all 13 GEN-SHIP references, and then reads only `.name`.
- `test_lazy_peierls_fetch_builds_requested_ref` builds `peierls_slab_2eg_2rg` and checks `.name` and identity.
- Neither test reads k or phi.
- The slab reference has 192 nodes, so 36 864 K elements per group at dps 30. At ~1.4 s per element (dps 20, measured) that is about 14 h per group `[R]`. **It cannot finish in 3000 s.**

### 4. `tests/gates/derivations/test_peierls_rank2_bc.py` (53.7 min) — C (strict), a refinement ladder of the generator

- `test_rank2_error_converges_monotonically_under_refinement` (3000 s, timeout) calls `composite_gl_r` + `build_volume_kernel(SLAB_POLAR_1D)` + `build_closure_operator(reflection="white")` at n_p ∈ {2,4,8} × p6, dps 20 (12, 24 and 48 nodes, so 2304 elements at the top rung), then a local power iteration.
- It asserts that the error falls by at least 3× per doubling against k_inf.
- Payload: `(N,N) float64` matrices. It shares no key with the other tests `[R]`.
- The rest of the file is ≤ 34 s per case.

### 5. `tests/gates/cp/test_peierls_rank_n_protocol.py` (50.9 min) — C (strict) + pin

- `test_f4_rich_vs_rich_panels_matches_pinned_baseline` ×6: each case runs 2 subprocesses (`_WORKER_CODE`). Each subprocess runs GEN-P, hollow sphere, white_f4, dps 15, tol 1e-12, at RICH (4 panels, p8, 64/64/64) and RICH+panels (5 panels, p8, 64).
- Payload: a scalar `k` over JSON stdout.
- Asserted against pinned signed-error literals, `tol 1e-6`.
- 12 distinct runs, no reuse. The sibling `test_f4_is_sign_stable…` would reuse the same pairs, but it skips all 6 (every status is "unresolved…").
- The subprocess is a clean process, so there is no shared state.

### 6. `tests/gates/cp/test_peierls_flux.py` (50.0 min, timeout) — G

- `continuous_get("peierls_slab_2eg_2rg")` is GEN-SHIP slab (it never finishes; see file 3).
- The S part is `solve_cp` on 32 cells, trivial `[R]`.
- The payload consumed is `ref.phi_cell_average(mesh, g)`, which calls the `phi` closure.
- Reuse: this reference is built 3× across 2 files (here, and twice in file 3). The registry memo does not help: file 3's fixture restores the saved registry, which discards its build, and shards are separate processes.

### 7. `tests/gates/derivations/test_peierls_multigroup.py` (49.8 min) — C\*, twin generators

- The unified GEN-P (slab, 2 panels × p3, dps 20) is checked against native `slab.solve_peierls_eigenvalue`.
- Native costs 1.4 s `[M]`, so about 100% of the time is in the unified path (K_vol).
- `test_2eg_2rg_parity_bit_exact` (white_f4, 1108 s) and `test_2g_2rg_vacuum_parity` (vacuum, 1071 s) use identical XS, radii and quadrature, so they share both per-group K_vol keys: **4 builds, 2 keys**.
- `test_1g_2rg_vacuum_parity` (605 s) has a distinct key.

### 8. `tests/gates/sn/verification/analytical/test_phase_c_crosscheck.py` (45.7 min) — G

- phase_d: GEN-TR (and its sphere counterpart) against the literals in `_SNAPSHOT_KEFFS`, so ~100% G.
- phase_e: GEN-TR `phi_g` against the snapshot `.npz`, so ~100% G.
- The cylinder MR reference is built 2× here with identical args (shared with file 2).
- `_run_cyl_1g_homogeneous_closed` is built 2× (cheap).
- `phase_e[cyl]` FAILED (#404).
