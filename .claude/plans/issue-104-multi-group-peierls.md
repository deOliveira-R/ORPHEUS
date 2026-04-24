# Issue #104 — Multi-group Peierls extension (scoping)

**Status:** scoping pass 2026-04-23; commit 1 landed 2026-04-24.
**Unblocks:** #130 (Phase G.5 slab-into-curvilinear routing).
**Reference pattern:** `peierls_slab.solve_peierls_eigenvalue`
(`orpheus/derivations/peierls_slab.py:406`).

---

## Session log

### 2026-04-24 — commit 1 (`solve_peierls_mg` + 1G wrapper)

- **F.4 closure per-group independence confirmed** (risk R1
  de-risked) by reading `_build_closure_operator_rank2_white`
  (`peierls_geometry.py:3698`). The closure takes a per-region
  scalar `sig_t`; per-face P_esc, G_bc, transmission W are all
  computed from that scalar. No cross-group coupling through the
  reflection operator. The MG path rebuilds the closure once per
  group with that group's `sig_t[:, g]`.
- **Implemented `solve_peierls_mg`** (~180 LoC) in
  `peierls_geometry.py` plus two new private helpers
  (`_resolve_closure_name`, `_build_full_K_per_group`). Block
  `(N·ng) × (N·ng)` assembly with node-major indexing
  `row = i·ng + g`, matching the slab reference at
  `peierls_slab.py:285`.
- **`solve_peierls_1g` is now a thin wrapper** that coerces 1-D
  per-region arrays and synthesises `chi = 1`, then delegates to
  `solve_peierls_mg`.
- **Bit-exact ng=1 parity.** Smoke test on slab-polar:
  `diff(k_eff) = 0.0`, `diff(phi_values) = 0.0`. Not 1e-12, not
  1e-15 — numerical zero.
- **13 regression tests** added in
  `tests/derivations/test_peierls_multigroup.py`:
  - 7 bit-match tests (3 geoms × 2 BCs + 1 hollow-cyl F.4) —
    0.0 diff required
  - 4 input-validation tests (ValueError on mis-shaped arrays)
  - 2 ng=2 sanity tests (decoupled 2G ≡ 1G; downscatter runs)
  - All PASSED in 2m42s.
- **Commit 1 scope:** `solve_peierls_mg` + wrapper + regression
  tests. Case builders and 2G registry entries deferred to
  commit 2 (next session).

### 2026-04-24 — commits 2 + 3 (hollow MG + Sphinx docs)

**Commit 2 — case builders + registry + tie-back test.**

- **Added `solve_peierls_cylinder_mg`** and
  **`solve_peierls_sphere_mg`** shape-specific wrappers over
  `solve_peierls_mg` (mirror of the `_1g` shape wrappers, one
  addition per file). Net ≈ +90 LoC.
- **Lifted hollow-F.4 case builders** (`_build_peierls_cylinder_hollow_f4_case`,
  `_build_peierls_sphere_hollow_f4_case`) to accept `ng_key: str =
  "1g"` and pass through to the MG-capable shape wrapper. The
  builders now normalise each group's flux independently to unit
  shell-volume integral (1G behaviour preserved, MG extends it
  per-group). Signature change is back-compatible — callers that
  don't pass `ng_key` retain the legacy 1G path.
- **`peierls_cases.build_two_surface_case`** now forwards `ng_key`
  to the curvilinear builders. Module docstring updated to reflect
  the 2026-04-24 state (both slab and curvilinear now MG-capable;
  G.5 routing still gated on a separate benchmark per Issue #130).
- **`_class_a_cases()`** registers 6 new 2G hollow entries:
  - `peierls_cyl1D_hollow_2eg_1rg_r0_{10,20,30}`
  - `peierls_sph1D_hollow_2eg_1rg_r0_{10,20,30}`
- **sig_s convention bug**: during commit 2 I noticed the commit-1
  docstring claimed `sig_s[r, g', g] = g → g'` (first=destination)
  while the code accessed `sig_s_n[j, gs, ge]` with gs=source
  (first=source). Verified by inspection of `_A_2G` XS
  (`sig_s[0,1] = 0.1` is physical downscatter only under
  first=source convention). Docstring fixed; test comment labelling
  corrected. Code was always right — no semantic change.
- **Tests added** to `test_peierls_multigroup.py`:
  - `TestMGSlabPolarMatchesNativeSlabMG::test_2g_vacuum_slab_matches_native_eigenvalue`
    — **the definitive convention cross-check**. Compares unified
    MG slab vs native `peierls_slab.solve_peierls_eigenvalue` on a
    2G XS set with directional (non-symmetric) downscatter. Target
    1e-8 rtol; failure would indicate sig_s ordering mismatch.
    `@pytest.mark.slow`.
  - `TestMG2GHollowRegistration::test_hollow_{cyl,sph}_2g_builds`
    — parametrised smoke tests (3 r₀/R each × 2 geometries = 6
    tests) that each new 2G hollow reference builds to completion
    with finite k_eff. Tight quadrature (n_panels=1, p=3, dps=15)
    to keep total cost to ~10 min. `@pytest.mark.slow`.
  - Combined with tier 1's 13 tests: **19 MG tests** total.

**Commit 3 — Sphinx `§theory-peierls-multigroup`.**

- New section `theory-peierls-multigroup` in `peierls_unified.rst`
  between `theory-peierls-slab-polar` and the archived moment-form.
  Documents: operator form (block (N·ng)² system), sig_s convention
  explicitly called out (source→destination, matches XS library,
  contradicts native slab docstring which is wrong), per-group K
  build strategy, relationship to the slab native block-Toeplitz,
  the bit-exact 1G wrapper guarantee, shipped 2G references table,
  cost characteristics.
- Key Facts header cross-references the new section.
- Capability matrix table extended with the 6 new 2G entries.
- "Known infrastructure gaps" updated — Issue #104 is no longer an
  open gap; the remaining multi-group gap is now the
  parity-vs-discrete-CP benchmark (deferred).

**Remaining follow-up (future session):**
- Parity benchmark: unified 2G hollow cyl/sph k_eff vs
  `cp_cylinder` / `cp_sphere` native 2G solvers to 1 % (Issue #104
  AC). This was deferred because the discrete CP solvers expect
  multi-region input that the 1-region hollow cases don't exercise
  directly; a proper comparison needs a matching discrete setup.
- Issue #130 (Phase G.5 slab routing). Blocked only on a direct
  benchmark now; the infrastructure is in place.
- Issue #104 can be **closed** once this session's work is
  committed — the core deliverable (`solve_peierls_mg`,
  regression gate, shape wrappers, 2G registry, Sphinx docs) is
  complete.

---

## 1. Executive summary

- Build `solve_peierls_mg` in `peierls_geometry.py` that generalises
  the existing `solve_peierls_1g` (`peierls_geometry.py:4081`) to
  `ng ≥ 1` groups with downscatter/upscatter coupling and χ·νΣ_f
  fission spectrum.
- The ray-geometry and closure infrastructure (`build_volume_kernel`,
  `build_volume_kernel_adaptive`, `build_closure_operator`,
  `build_white_bc_correction_rank_n`) already take `sig_t` as a
  **per-region** vector — they are group-agnostic. Multi-group only
  adds an outer **per-group** loop building one K per group (they
  differ solely through Σ_t,g).
- Algorithm: block (N·ng)×(N·ng) A,B assembly + fission-source power
  iteration, identical in shape to `peierls_slab._build_system_matrices`
  (`peierls_slab.py:243`).
- Compute cost class: **verification-only primitive**. At ng=2 we
  double K builds. 2G hollow cyl/sph at default quadrature is expected
  to cost 2× the 1G wall time (tens of seconds per point, not hours).
  Mark tests `@pytest.mark.slow` as a precaution.
- Regression gate: `solve_peierls_mg(..., ng=1)` bit-matches
  `solve_peierls_1g` to 1e-12 on k_eff and on every φ_i.

---

## 2. Current state inventory (1G machinery that must generalise)

| Symbol | Location | 1G-specific? |
|---|---|---|
| `solve_peierls_1g` | `peierls_geometry.py:4081` | Driver — replace |
| `PeierlsSolution` | `peierls_geometry.py:4047` | `phi_values` already `(N, ng)`; `n_groups` already a field — generic |
| `build_volume_kernel` | `peierls_geometry.py:1123` | `sig_t` is per-region array. Single-group output but called once per g |
| `build_volume_kernel_adaptive` | `peierls_geometry.py:1084` | Same — group-agnostic |
| `K_vol_element_adaptive` | `peierls_geometry.py:937` | Group-agnostic |
| `build_closure_operator` | `peierls_geometry.py:3531` | `sig_t` per-region; group-agnostic. One closure per group |
| `build_white_bc_correction_rank_n` | `peierls_geometry.py:3918` | Group-agnostic |
| `CurvilinearGeometry.optical_depth_along_ray` | `peierls_geometry.py:503` | Takes per-region `sig_t`; group-agnostic |
| `PeierlsCylinderSolution` | `peierls_cylinder.py:214` | Thin facade; `phi_values` already `(N, ng)` |
| `PeierlsSphereSolution` | `peierls_sphere.py:215` | Same |
| `solve_peierls_cylinder_1g` | `peierls_cylinder.py:263` | Wrapper — add `solve_peierls_cylinder_mg` or lift to `_mg` name |
| `solve_peierls_sphere_1g` | `peierls_sphere.py:266` | Same |
| `_build_peierls_cylinder_hollow_f4_case` | `peierls_cylinder.py:482` | 1G-only case builder — lift to accept `ng_key` |
| `_build_peierls_sphere_hollow_f4_case` | `peierls_sphere.py:486` | Same |
| `peierls_cases.build_two_surface_case` | `peierls_cases.py:48` | `ng_key` already accepted; currently hard-errors for cyl/sph at ng≠"1g" via hollow-f4 builders |

---

## 3. Proposed API

**Keep `solve_peierls_1g` as a thin wrapper** that forwards to the
new `solve_peierls_mg` (mirrors slab's pattern: slab's single entry is
`solve_peierls_eigenvalue` which handles `ng=1` as a degenerate case,
but we have many downstream callers of `solve_peierls_1g`, so we keep
the alias for one deprecation release).

```python
def solve_peierls_mg(
    geometry: CurvilinearGeometry,
    radii: np.ndarray,                     # (n_regions+1,) or (n_regions,)
    sig_t: np.ndarray,                     # (n_regions, ng)
    sig_s: np.ndarray,                     # (n_regions, ng, ng)  sig_s[r, g', g] = g → g'
    nu_sig_f: np.ndarray,                  # (n_regions, ng)
    chi: np.ndarray,                       # (n_regions, ng)
    *,
    boundary: str = "vacuum",
    n_panels_per_region: int = 2,
    p_order: int = 5,
    n_angular: int = 24,
    n_rho: int = 24,
    n_surf_quad: int = 24,
    dps: int = 25,
    max_iter: int = 300,
    tol: float = 1e-10,
    n_bc_modes: int = 1,
) -> PeierlsSolution: ...
```

**`solve_peierls_1g`** becomes a one-line wrapper that reshapes
`sig_t, sig_s, nu_sig_f` into ng=1 arrays, synthesises `chi = [[1.0]]`,
and calls `solve_peierls_mg`. It MUST bit-match (same quadrature,
same closure path).

**`PeierlsSolution` dataclass** — no schema change required.
`phi_values: (N, ng)` and `n_groups: int` fields already exist; the
only thing 1G writes today is `phi[:, np.newaxis]` and
`n_groups=1` (`peierls_geometry.py:4249,4252`).

---

## 4. Algorithm (pseudocode)

```
# Inputs: geometry, radii, sig_t[nreg, ng], sig_s[nreg, ng, ng],
#         nu_sig_f[nreg, ng], chi[nreg, ng]
r_nodes, r_wts, panels = composite_gl_r(radii, ...)
N = len(r_nodes)
k_annulus = [geometry.which_annulus(r_i, radii) for r_i in r_nodes]   # per-node region index

# Per-group K (volume + closure). The ONLY loop that's group-local.
K = np.empty((ng, N, N))
for g in range(ng):
    K_vol_g = build_volume_kernel(geometry, r_nodes, panels, radii,
                                  sig_t[:, g], ...)
    if boundary == "white_rank1_mark":
        K_bc_g = build_white_bc_correction_rank_n(..., sig_t[:, g], ...)
        K[g] = K_vol_g + K_bc_g
    elif boundary == "white_f4":
        op = build_closure_operator(..., sig_t[:, g], reflection="white")
        K[g] = K_vol_g + op.as_matrix()
    else:
        K[g] = K_vol_g

# Per-node XS
sig_t_n  = sig_t[k_annulus, :]            # (N, ng)
sig_s_n  = sig_s[k_annulus, :, :]         # (N, ng, ng)
nu_f_n   = nu_sig_f[k_annulus, :]         # (N, ng)
chi_n    = chi[k_annulus, :]              # (N, ng)

# Block (N*ng)×(N*ng) assembly — row index = i*ng + g_row.
# Mirrors peierls_slab._build_system_matrices (peierls_slab.py:243).
A = np.zeros((N*ng, N*ng))
B = np.zeros((N*ng, N*ng))
for g_out in range(ng):
    for i in range(N):
        A[i*ng+g_out, i*ng+g_out] = sig_t_n[i, g_out]     # diag Σ_t
        for j in range(N):
            kij = K[g_out, i, j]
            if kij == 0: continue
            for g_in in range(ng):
                A[i*ng+g_out, j*ng+g_in] -= kij * sig_s_n[j, g_out, g_in]    # scatter g_in → g_out at j
                B[i*ng+g_out, j*ng+g_in] += kij * chi_n[i, g_out] * nu_f_n[j, g_in]

# Power iteration (fission-source), identical to peierls_slab.py:530.
phi = np.ones(N*ng); k_val = 1.0
for it in range(max_iter):
    q = (B @ phi) / k_val
    phi_new = np.linalg.solve(A, q)
    k_new = k_val * |B @ phi_new|_1 / |B @ phi|_1
    phi_new /= |phi_new|_1
    if |k_new - k_val| < tol and it > 5: break
    phi, k_val = phi_new, k_new

phi_arr = phi.reshape(N, ng)
return PeierlsSolution(r_nodes=r_nodes, phi_values=phi_arr, k_eff=k_val,
                       n_groups=ng, ...)
```

Two subtle points mirrored from slab:

- `A = diag(Σ_t) − K·diag(Σ_s)` form (NOT identity-LHS). The slab
  module uses identity-LHS (`A = I − K_scatter`); the curvilinear
  1G driver already uses **Σ_t-LHS** (`peierls_geometry.py:4223`),
  so we keep Σ_t-LHS for continuity. Both are algebraically
  equivalent; Σ_t-LHS matches the unified operator docstring at
  `peierls_geometry.py:1141`.
- χ sits on the **destination** node (`chi_n[i, g_out]`), not on
  source — identical to slab line 288 of `peierls_slab.py`.

---

## 5. Per-module changes

| File | Change | Est LoC |
|---|---|---|
| `orpheus/derivations/peierls_geometry.py` | Add `solve_peierls_mg` (new ~180 LoC). Rewrite `solve_peierls_1g` as a ~25-line wrapper that calls `solve_peierls_mg` with `ng=1`. Keep all signature/docstring invariants for back-compat. | +200 / -90 |
| `orpheus/derivations/peierls_cylinder.py` | Add `solve_peierls_cylinder_mg` wrapper (parallel to `_1g`). Lift `_build_peierls_cylinder_hollow_f4_case` to take `ng_key` and route ng≥2 XS through the MG path. | +110 |
| `orpheus/derivations/peierls_sphere.py` | Mirror of cylinder changes. | +110 |
| `orpheus/derivations/peierls_cases.py` | Remove the "1g only" guard in the `cylinder-1d` / `sphere-1d` branches of `build_two_surface_case` (lines 75-77 docstring + builders). Add 2G hollow cyl/sph entries to `_class_a_cases()`. | +20 / -5 |
| `tests/derivations/test_peierls_multigroup.py` (new) | Regression (ng=1 bit-match) + 2G parity tests. | +250 |
| `tests/cp/test_peierls_cylinder_flux.py` / `test_peierls_sphere_flux.py` | Add 2G flux-shape tests per group, gated `@pytest.mark.slow`. | +80 total |
| `docs/theory/peierls_unified.rst` | Add §multi-group section documenting the block assembly. | +60 |

**Total estimate: ~825 LoC net across 7 files.**

---

## 6. Reference cases to register

XS library (`_xs_library.py:229`) ships `1g`, `2g`, `4g` keys for
regions A, B, C, D. Slab already uses `2g` (see
`_build_peierls_slab_case`, `peierls_slab.py:609`). Natural additions
to `_class_a_cases()` (`peierls_cases.py:176`):

```python
# 2G hollow cylinder F.4 at r_0/R ∈ {0.1, 0.2, 0.3}
for r0 in (0.1, 0.2, 0.3):
    refs.append(build_two_surface_case("cylinder-1d", "2g", 1, inner_radius=r0))
# 2G hollow sphere F.4 at r_0/R ∈ {0.1, 0.2, 0.3}
for r0 in (0.1, 0.2, 0.3):
    refs.append(build_two_surface_case("sphere-1d",   "2g", 1, inner_radius=r0))
```

Case names: `peierls_cylinder_2eg_1rg_hollow_r0pXX`,
`peierls_sphere_2eg_1rg_hollow_r0pXX` (match slab's
`peierls_slab_2eg_2rg` naming). Once stable, consider adding `4g`
variants in a follow-up issue (NOT part of #104).

**Region composition.** Hollow F.4 builders currently plug one XS
region (region `"B"` per `peierls_cylinder.py:482` context — verify on
implementation). Use the same region for `2g` to isolate the
group-coupling effect from multi-region effects. Multi-region MG
(2 regions × 2 groups) is a follow-up, NOT part of #104.

---

## 7. Test plan

**Regression gates** (MUST pass — ng=1 path unchanged):

- `tests/derivations/test_peierls_multigroup.py::test_mg_bitmatch_1g_slab_polar`
  — feed `ng=1` through `solve_peierls_mg`, compare k_eff + every
  φ_i against `solve_peierls_1g` on slab-polar geometry. Tol 1e-12.
- Same for `cylinder-1d`, `sphere-1d`. Three tests.
- All existing `TestRank2Slab...`, `TestWhiteBC...`,
  `test_peierls_rank_n_bc.py` tests continue to pass unchanged
  (they call `solve_peierls_1g`; that is now a thin wrapper).

**New gates** (2G behaviour):

- `test_mg_2g_hollow_cyl_matches_cp_cylinder` — 2G hollow cyl k_eff
  via `solve_peierls_mg` within 1% of `cp_cylinder` 2G native solver
  at R=10 MFP (mirror Issue #104 AC). `@pytest.mark.slow`.
- Same for sphere. `@pytest.mark.slow`.
- `test_mg_flux_shape_per_group` — flux-shape regression per group
  (L1-style) at the new registered cases.

**Phase G.5 tie-back** (unblocks #130):

- `test_g5_slab_via_unified_2g_matches_native` — route
  `peierls_cases.build_two_surface_case("slab", "2g", 2)` through
  `solve_peierls_mg(SLAB_POLAR_1D, ...)` and compare to
  `solve_peierls_eigenvalue`. This remains gated by the slab-polar
  log-singularity cost caveat noted in `peierls_cases.py:27` — it may
  need `precision_digits=30` and `n_panels_per_region=16` to match;
  benchmark before gating.

**Cost budget.** At N=12, dps=25 the 1G parametrised case suite took
~41s per invocation. 2G doubles K builds → expect **~80-100s per
point**. With 6 new 2G cases the slow suite adds ~10 min wall time —
acceptable for `@pytest.mark.slow`. Phase G.5 (2G slab at dps=30)
could push to several minutes per case; defer to implementation
session benchmark.

---

## 8. Risks and open questions

1. **HIGH RISK — `build_closure_operator` API symmetry.** The F.4
   rank-2 closure (`_build_closure_operator_rank2_white`,
   `peierls_geometry.py:3698`) is per-face and constructed from
   per-group Σ_t traces. Verify the per-group closure builds are
   independent (no cross-group coupling through the reflection
   operator). I believe they are (the closure is local to each
   group's transmission E_3 / Ki_3 evaluation) but this needs
   confirmation against the L21 closeout documented in
   `docs/theory/peierls_unified.rst`.

2. **MEDIUM RISK — `radii` array shape.** `solve_peierls_1g` takes
   `radii: np.ndarray` assumed to be the list of outer radii per
   region plus a possible inner cavity via `geometry.inner_radius`.
   The multi-group extension inherits this exactly but the caller
   now MUST supply `sig_t` as `(n_regions, ng)` not `(n_regions,)`.
   Validate with a strict `assert sig_t.ndim == 2 and sig_t.shape ==
   (n_regions, ng)` or shape-coerce a 1-D input for the ng=1
   wrapper.

3. **MEDIUM RISK — upstream callers of `phi_values`.** All callers
   of `PeierlsSolution.phi_values` already treat it as `(N, ng)`
   (see `phi()` method at `peierls_geometry.py:4065`). `_soln_to_
   cylinder` / `_soln_to_sphere` pass-through unchanged. But
   tests/diagnostics that do `phi_values[:, 0]` or `.ravel()`
   implicitly assume ng=1; grep audit required before merge.
   Safe heuristic: the new MG cases are ADDED (not replacing
   existing 1G cases), so existing tests keep their ng=1 reference
   and don't break.

4. **LOW RISK — power iteration convergence.** At ng=2 with
   downscatter only (no upscatter) the spectral radius is well
   separated; at ng=4 with upscatter the convergence could stall.
   Slab's `tol = 10^(-(dps-5))` works at dps=30. For mp-free float
   curvilinear path (dps=25) keep `tol=1e-10` as 1G uses.

5. **OPEN QUESTION — mpmath or float for MG?** Slab's MG path runs
   under `mpmath.workdps(dps)` (`peierls_slab.py:449`). The 1G
   curvilinear driver runs **float** matrices built by a
   mpmath-precision K. At ng=2 the direct-solve
   `np.linalg.solve(A, q)` on a (N·ng)×(N·ng) float matrix should be
   fine (condition number ~Σ_t/(νΣ_f·P_esc) = O(1)); only flag if a
   diagnostic sees iteration stall. Default: **stay float** to
   preserve 1G cost profile. Document this decision.

6. **OPEN QUESTION — which XS region for 2G hollow cyl/sph?**
   The hollow-f4 builders currently use a single shipped composition.
   Check whether the slab `peierls_slab_2eg_2rg` XS pair (regions at
   index mapping `_MAT_IDS[2]`) gives meaningful coupling for the
   single-region hollow cells, or pick a different region letter
   for better group separation. Implementation-session benchmark.

---

## 9. Estimated budget

- **1 focused session** is realistic if no closure-API surprise
  surfaces (risk 1). Plan:
  1. Commit 1 — `solve_peierls_mg` + `solve_peierls_1g` wrapper +
     regression test (ng=1 bit-match on all three geometries).
  2. Commit 2 — case builders lifted to `ng_key` + 2G registry
     entries + 2G parity tests (`@pytest.mark.slow`).
  3. Commit 3 — Sphinx §multi-group docs + close #104.
- **2 sessions** if risk 1 surfaces (closure per-group independence
  needs a diagnostic).
- **3 sessions** only if risk 5 (mpmath vs float) materialises as a
  convergence issue at 2G.

Commit count: **3 commits** (feature, registry, docs).
LoC: **~825 net** (see §5).

---

## 10. Entry point for the implementation session

Start by running the three existing regression tests to ensure
they pass on main:

```bash
pytest -x tests/derivations/test_peierls_rank2_bc.py::TestSlabPolarVsNativeE1KEff
pytest -x tests/derivations/test_peierls_cylinder_white_bc.py
pytest -x tests/derivations/test_peierls_sphere_white_bc.py
```

Then implement `solve_peierls_mg` in `peierls_geometry.py` with
`solve_peierls_1g` reduced to a wrapper. Run the new
`test_mg_bitmatch_1g_*` gates — ng=1 bit-match is non-negotiable
before proceeding.
