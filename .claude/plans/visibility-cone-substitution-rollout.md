# Plan 1: Visibility-cone u² substitution — stepwise rollout across ORPHEUS quadratures

**Author**: Claude Opus 4.7, 2026-04-28
**Predecessor**: branch `feature/peierls-specular-bc` commit `4dc03cf` (Phase 5 retreat). Phase 5 R3 PRIMARY/SECONDARY surfaced this technique as a **portable, promotion-worthy primitive** even though Phase 5 itself was abandoned — see `phase5_round3_visibility_cone_quad.py` in derivations/diagnostics.
**Proposed branch**: `feature/visibility-cone-quadrature` off `feature/peierls-specular-bc` (Phase 4 + 5a + 5-retreat all sit there; this builds on top).

---

## Status (2026-04-29): COMPLETE — superseded by Q1→Q6→L3 architecture

**Closure commit**: `e8ad150` (`feat(quadrature): L3 — surface_centred_angular_quadrature recipe + cylinder G_bc migration`).
**Branch landed on**: `feature/quadrature-architecture` (initially branched from `feature/peierls-specular-bc`).

The work that actually shipped is broader than this plan envisioned. Rather than threading a single `gauss_legendre_visibility_cone(y_min, y_max, n)` utility through ~6-8 call sites, the rollout grew into a **unified `Quadrature1D` contract** with three composable geometry-aware recipes plus a sibling `AdaptiveQuadrature1D` for verification-tier integrals. Every primitive in `peierls_geometry.py` that the plan listed (and several it did not) now consumes the contract through one of:

- `chord_quadrature` (h-space chord-impact-parameter integrals — sphere/cyl T-matrix and P_ss)
- `observer_angular_quadrature` (observer-centred ω-sweep with `arcsin(r_k/r_obs)` tangent subdivision — `*_mode`, `*_outer_mode`, `*_inner_mode`, `*_outer_mode_marshak`, `*_inner_mode_marshak`, `compute_G_bc_cylinder_3d{,_mode}`, etc.)
- `surface_centred_angular_quadrature` (surface-centred φ-sweep with chord-quadratic tangent subdivision — the four legacy cylinder G_bc branches)
- `adaptive_mpmath` (verification-tier breakpoint-hint adaptive integration — `build_volume_kernel_adaptive` ω and ρ integrals)

**Originating commits** (in dependency order): `b281a97` (Q1) → `c4cda20` (Q2) → `50f05ae` (Q3) → `7893bf6` (Q4) → `123f9ec` (Q5) → `6c67466` (Q6) → `e8ad150` (L3). 7 commits, 7 sessions; planned ~5.

**Audit table** (acceptance criteria from §"Acceptance criteria (overall)"):

| # | Criterion | Status |
|---|-----------|--------|
| (a) | `gauss_legendre_visibility_cone` shipped + L0 test | ✅ DONE — promoted from `_kernels.py` to `_quadrature.py` as a `Quadrature1D` constructor; 6 L0 tests in `test_quadrature.py` |
| (b) | `compute_T_specular_*` use new quadrature; N=4 overshoot reduced | ✅ DONE — sphere/cyl through `chord_quadrature` (Q2); slab plan-exempt (geometric immunity); high-N rank-gating documented at `test_specular_*` |
| (c) | `compute_P_esc_*` / `compute_G_bc_*` mode primitives use new quadrature; rank-1 algebraic identities @ 1e-14 | ✅ DONE — Q3 migrated all observer-centred per-face mode primitives; L3 migrated the four legacy cylinder G_bc branches; `rank1_equals_hebert` canary preserved at `rtol=1e-8` |
| (d) | `compute_P_ss_*` use new quadrature; T_00 = P_ss preserved | ✅ DONE — both routed through `chord_quadrature` in h-space (Q2); identity tested in `test_chord_quadrature_sphere_T00_equals_P_ss` at bit-equality |
| (e) | `build_volume_kernel_adaptive` uses new quadrature | ✅ DONE — both ω and ρ integrals through `adaptive_mpmath` (Q5/Q6) |
| (f) | All Peierls + specular tests pass without modification | ✅ DONE — verified across L3 closure runs (50 quadrature L0 + 7 cylinder Mark closure + 18 + 23 broader cylinder/closure/multigroup/rank-N) |
| (g) | Documented improvement (Sphinx table) | ✅ DONE — `docs/theory/peierls_unified.rst` §22.7 (vis-cone substitution), §22.8 (surface-centred recipe + chord-quadratic derivation), §22.9 (rollout outcome audit + landing map) |
| (h) | ~3-5 commits scaled by phase boundaries | ✅ DELIVERED as 7 commits — broader scope (full contract + 3 recipes + adaptive sibling) than originally planned |

**Residual legacy consumers (intentional, tracked)**: three call sites in `peierls_geometry.py` remain on raw `np.polynomial.legendre.leggauss`:

- `compute_T_specular_slab` (line 2320) — plan-exempt per §1B step 3 (geometric immunity at µ→0). Tracked: Issue #136 (consistency-only contract migration).
- `_compute_K_bc_sanchez` (line 2639) — verification-tier reference, out of original plan scope. Tracked: Issue #136.
- `build_volume_kernel` (line 1194, the non-adaptive variant) — predates the recipe abstraction; not enumerated in the plan. Tracked: Issue #135 (cleanup with migration sketch).

**Phase 1F (CP module) — deferred**: marked optional/low-priority in the original plan. Tracked: Issue #134.

The original "Implementation phases" section below is preserved for historical reference; the live audit is in this Status header.

---

## Context

**Why this change**. Multiple ORPHEUS quadratures evaluate integrals of the form `∫_{y_min}^{y_max} f(y) dy` where `f(y)` carries a `√(y_max² − y²)` (or `1/√`) factor that vanishes algebraically at the visibility-cone endpoint. The current treatment in `composite_gl_r` uses **panel subdivision at impact-parameter knots + plain GL per panel** — this gives algebraic-order convergence (~Q^(-2.5) per panel) because the leading endpoint singularity (square-root type from the chord Jacobian) is left in the integrand, not absorbed.

The substitution `u² = (y² − y_min²)/(y_max² − y_min²)` (or its symmetric variant on the upper bound) absorbs the singular Jacobian into the change-of-variables, producing a **smooth integrand on u ∈ [0, 1]**. Plain GL on u-space then converges **spectrally** — Phase 5 R3 PRIMARY measured 1e-3 → 1e-9 at Q=128 for off-diagonal pairs, ~5-6 orders of magnitude per Q-doubling vs algebraic for the panel-GL approach.

**What prompted this now**. The Phase 5 retreat closed `closure="specular_continuous_mu"` as a research-grade-only direction, but the diagnostic effort produced this clean portable improvement. Promoting it across the existing Peierls + CP quadratures should yield direct accuracy/speed wins **without changing any solver semantics** — same physics, same dispatch, just sharper quadrature.

**Intended outcome**:
- A `gauss_legendre_visibility_cone(y_min, y_max, n)` utility shipped in `_kernels.py`
- Existing chord-based primitives (`compute_T_specular_*`, `compute_P_ss_*`, `compute_P_esc_*`, `compute_G_bc_*` mode variants) optionally use it
- At each rollout step: empirical verification that Phase 4 / CP / Peierls reference k_eff and convergence ladders are improved or unchanged (no regressions)
- Phase 4 `closure="specular_multibounce"` rank-N matrix elements get tighter (the matrix-Galerkin form is currently quadrature-limited at high N — this should let N=3 → 4 envelope without overshoot, possibly more)

## Implementation phases

Each phase is ~1-3 hours; verify before moving to next. Halt if any phase shows regression.

### Phase 1A — Utility shipped (1 session, ~1.5 hours)

1. Add `gauss_legendre_visibility_cone(y_min, y_max, n)` to `orpheus/derivations/_kernels.py`:
   ```python
   def gauss_legendre_visibility_cone(y_min: float, y_max: float, n: int) -> tuple[np.ndarray, np.ndarray]:
       """GL quadrature on [y_min, y_max] absorbing the √(y_max² − y²)
       (or √(y² − y_min²)) endpoint singularity via u² = (y² − y_min²)/(y_max² − y_min²).
       Returns (y_pts, y_wts) suitable for ∫ f(y) dy where f has the
       single-endpoint visibility-cone singularity.
       """
       u_nodes, u_wts = np.polynomial.legendre.leggauss(n)
       u = 0.5 * (u_nodes + 1.0)
       u_w = 0.5 * u_wts
       y2_diff = y_max ** 2 - y_min ** 2
       y_pts = np.sqrt(y_min ** 2 + u ** 2 * y2_diff)
       # dy/du = u · y2_diff / √(y_min² + u² · y2_diff)
       y_wts = u_w * u * y2_diff / y_pts
       return y_pts, y_wts
   ```
2. Add SymPy/numerical L0 test at `tests/derivations/test_kernels.py::test_visibility_cone_substitution_spectral`:
   - Test on `f(y) = 1 / √(R² − y²)` integrated on `[0, R]`. Closed-form: `π/2`.
   - Plain GL Q=64 vs visibility-cone Q=64 — vis-cone hits machine precision at Q=16 while plain GL plateaus at ~1e-4.
   - Test on `f(y) = (R² − y²)^(3/2)` (smooth at endpoint) — both should agree.
   - Test off-diagonal pair `f(y) = √(R² − y²) · √(R² − y² − ε²)` — vis-cone wins by orders of magnitude.
3. Sphinx documentation in `docs/theory/peierls_unified.rst` § "Visibility-cone quadrature primitive": derivation of the change-of-variables, spectral-vs-algebraic convergence comparison plot.

**Acceptance**: utility shipped, L0 test passes, Sphinx documented. No production code changed.

### Phase 1B — `compute_T_specular_*` rollout (1 session, ~2.5 hours)

The Phase 4 multi-bounce T-matrix integrals are the highest-leverage site (rank-N matrix elements at high N currently quadrature-limited; the per-Q noise is what triggers the N≥4 UserWarning).

1. **Sphere** (`compute_T_specular_sphere`, line 2032+): the `tau_arr(µ)` antipodal-chord τ uses `h = R√(1 − µ²)` impact parameter; the chord segments per annulus are `2·√(r_k² − h²)` — visibility-cone singularity at `h → r_k` (which is `µ → √(1 − (r_k/R)²)`). Replace the inner GL on µ ∈ [0, 1] with **visibility-cone GL per annular interval** (currently single panel; the kink at impact-parameter knot is the actual visibility-cone boundary).

2. **Cylinder** (`compute_T_specular_cylinder_3d`, line 2260+): same treatment. The α-integral over [0, π/2] has `cos α` factor going to 0 at α=π/2 — visibility-cone singularity.

3. **Slab** (`compute_T_specular_slab`, line 2158+): integrand has `e^{-τ_total/µ}` going to 0 at µ=0 (geometric immunity), no visibility-cone singularity per se. Skip this branch — current GL is already spectral.

4. **Smoke test**: rerun `tests/derivations/test_peierls_specular_bc.py` — ensure all 24 specular tests still pass. Then run thin-sphere convergence ladder (τ_R = 2.5, fuel-A-like, k_inf = 0.20833) at N ∈ {1, 2, 3, 4, 6, 8} and compare to the Phase 4 baseline:
   - Pre-rollout: N=4 +0.43% overshoot; N=8 +5.62%
   - Hypothesis: vis-cone gives tighter T-matrix elements; the matrix-Galerkin overshoot at N≥4 is structural (proven in Phase 5) so won't disappear, but should be smaller magnitude.
   - Pin the new Phase-4 ladder numbers in test docstrings.

5. Update Phase 4 sphere `compute_T_specular_sphere` docstring to reference the new quadrature.

**Acceptance**: 24/24 specular tests pass. Phase 4 N=4 overshoot magnitude is demonstrably smaller than the pre-rollout baseline (e.g., +0.43% → +0.10% or better). No regression at N ∈ {1, 2, 3}.

### Phase 1C — `compute_P_esc_*` / `compute_G_bc_*` mode primitives (1 session, ~3 hours)

The 3-D Knyazev primitives for cylinder (`compute_P_esc_cylinder_3d_mode`, line 1655; `compute_G_bc_cylinder_3d_mode`, line 1747) and the per-face mode variants for sphere/slab (`compute_P_esc_outer_mode`, `compute_G_bc_outer_mode`, etc., lines 2874+, 3037+, 3237+) all do chord traversal with visibility-cone singularity.

1. Identify each integration loop with the `√(r_k² − h²)` or equivalent structure.
2. Replace inner GL with `gauss_legendre_visibility_cone` per annular interval.
3. Cross-check via the Phase 4 specular convergence tests + the `closure="white_hebert"` rank-1 Hébert test (this exercises P_esc and G_bc primitives).

**Acceptance**: existing `closure="white_hebert"`, `closure="specular"` (rank-1 Mark/Hébert algebraic identities), and `closure="specular_multibounce"` rank-1 algebraic identities all preserved bit-equally (they're tested at 1e-8 to 1e-14). Improved convergence at higher N visible in error magnitudes.

### Phase 1D — `compute_P_ss_*` (1 session, ~1.5 hours)

`compute_P_ss_sphere` (line 1921) and `compute_P_ss_cylinder` (line 1826) compute surface-to-surface kernels via chord integrals. Same pattern.

1. Roll the substitution through.
2. Verify the rank-1 algebraic identity `T_00 = P_ss^cyl/sphere` (Phase 4 tests at 1e-14) is preserved.

**Acceptance**: 1e-14 algebraic identity preserved. Improved high-rank cyl/sphere Hébert white BC accuracy if any (white_hebert is rank-1 only so probably no observable change at the closure level).

### Phase 1E — `build_volume_kernel_adaptive` (1 session, ~3 hours)

The K_vol assembly currently uses `mpmath.quad` adaptive integration with panel boundaries as breakpoints. The chord-bound `ρ_max(r_i, Ω)` has `√(R² − r_i²·sin²(Ω))` structure — visibility-cone singularity in Ω at the grazing angle.

1. Identify the Ω-integration loop in `build_volume_kernel_adaptive` (peierls_geometry.py:1131).
2. Apply visibility-cone substitution at the grazing-Ω endpoint.
3. Cross-check via existing Peierls reference test fixtures (vacuum sphere, vacuum cyl, multi-region Peierls).

**Acceptance**: existing Peierls reference solutions (`solve_peierls_1g` at vacuum BC for thin/thick sphere/cyl/slab) are preserved to 1e-12. Quadrature speed improvement measurable (lower dps requirement for same accuracy).

### Phase 1F — CP module rollout (optional; 1-2 sessions)

The CP modules (`cp_cylinder.py`, `cp_sphere.py`) consume `chord_half_lengths` via `cp_geometry.py:287`. Replace plain GL on impact parameter with visibility-cone GL.

This is **optional and lower-priority** — CP is a discrete solver with its own discretization-dominant errors; sharpening its underlying quadrature gives little observable benefit until other discretization sources are addressed. Defer unless a CP regression test specifically benefits.

## Verification plan

After each phase, run:

```bash
# Specular tests (24 cases) — Phase 4 cross-check
python -m pytest tests/derivations/test_peierls_specular_bc.py -v

# Peierls reference tests (Phase 4.x convergence ladders)
python -m pytest tests/derivations/test_peierls_*.py -v

# CP cross-checks (where applicable)
python -m pytest tests/cp/ -v -k "peierls"

# Sphinx still builds clean
sphinx-build -W docs docs/_build/html
```

After full rollout (1A-1E), produce a one-table summary in `docs/theory/peierls_unified.rst` showing:
- Pre-rollout vs post-rollout k_eff at N ∈ {1, 2, 3, 4, 6, 8} for sphere/cyl/slab thin-cell
- Pre-rollout vs post-rollout matrix-element accuracy (T-matrix entries)
- Quadrature-node count to reach 1e-6 / 1e-9 / 1e-12

## Acceptance criteria (overall)

Plan 1 is complete when:

- (a) `gauss_legendre_visibility_cone` shipped in `_kernels.py` with L0 test
- (b) Phase 4 `compute_T_specular_*` integrals use the new quadrature; N=4 sphere/cyl overshoot magnitude reduced
- (c) `compute_P_esc_*` / `compute_G_bc_*` mode primitives use the new quadrature; rank-1 algebraic identities preserved at 1e-14
- (d) `compute_P_ss_*` use the new quadrature; rank-1 algebraic identity `T_00 = P_ss` preserved
- (e) `build_volume_kernel_adaptive` (where applicable) uses the new quadrature
- (f) ALL existing Peierls + specular tests pass without modification (no regression)
- (g) Documented improvement (table in Sphinx)
- (h) ~3-5 commits scaled by phase boundaries

## Risks

- **R1 Visibility-cone subdivision interacts with composite-GL panel boundaries**. Need to ensure the impact-parameter knots that segment the chord are aligned with the substitution intervals. Mitigation: roll out per-annulus using `chord_half_lengths` panel structure.
- **R2 Algebraic identities sensitive to quadrature convention**. The Phase 4 `T_00 = P_ss`, `compute_P_ss_cylinder = T_00^cyl` identities are tested at 1e-14. If both sides use the new quadrature, identities preserved bit-equally; if only one side, drift. Mitigation: roll utility through both sides of every algebraic identity simultaneously.
- **R3 `mpmath.quad` adaptive in `build_volume_kernel_adaptive` may not benefit**. mpmath is already adaptive; substitution may be redundant. Mitigation: profile first; skip Phase 1E if `mpmath.quad` already converges spectrally.
- **R4 Phase 4 high-N overshoot is structural** (Phase 5 retreat conclusion) so vis-cone won't fix the rank gating, only sharpen the quadrature. Set realistic expectations: N=4 overshoot reduces magnitude but doesn't disappear.

## Non-goals

- **Not** revisiting the Phase 5 continuous-µ closure. Phase 5 is permanently retreated; vis-cone substitution doesn't change that diagnosis.
- **Not** changing solver semantics (closure dispatch, k_eff convergence, BC handling). Pure quadrature improvement.
- **Not** removing existing `composite_gl_r` panel-subdivision logic. The new utility coexists; production code chooses per-call.
- **Not** rolling through MOC / SN modules (they're discrete solvers with different quadrature concerns).

## Estimated budget

- **Best case**: 5 sessions for 1A through 1E. ~12-15 hours total.
- **Mid case**: 6-7 sessions if R1/R2 surface and need extra cross-checks.
- **Worst case**: 8 sessions if Phase 1E (`build_volume_kernel_adaptive`) needs significant restructuring of the mpmath integration path.

LoC delta:
- New utility + L0 test: ~80 LoC
- Phase 1B-1D rollout: ~150-200 LoC across 6-8 functions
- Phase 1E rollout: ~50 LoC
- Sphinx + tables: ~100 lines

Total: ~400-500 LoC. Commit count: 4-6.

## Critical files

- `orpheus/derivations/_kernels.py` — utility lives here
- `tests/derivations/test_kernels.py` — L0 test
- `orpheus/derivations/peierls_geometry.py` — Phase 4 multi-bounce T-matrix sites (lines 2032+, 2158+, 2260+); P_esc / G_bc mode primitives (lines 1655+, 1747+, 2874+, 3037+, 3237+); P_ss (lines 1826+, 1921+); volume kernel (1131+)
- `tests/derivations/test_peierls_specular_bc.py` — pinning Phase 4 cross-check
- `docs/theory/peierls_unified.rst` — Sphinx documentation
- `derivations/diagnostics/diag_phase5_round3_visibility_cone_quad.py` — Phase 5 R3 SECONDARY's working diagnostic (3/3 PASS) for the substitution; promote-worthy primitive
