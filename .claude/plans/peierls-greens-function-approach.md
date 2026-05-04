# Plan 2: Reorganize Peierls + start Green's function reference approach

**Author**: Claude Opus 4.7, 2026-04-28
**Predecessor**: branch `feature/peierls-specular-bc` commit `4dc03cf` (Phase 5 retreat). The Phase 5 investigation pointed at the Green's function approach as the structurally correct continuous reference for transport with BC absorbed — Sanchez 1986 Eq. (3) defines the Green's function `t(r', Ω' → r, Ω)` directly, with BC encoded into the kernel via Eq. (A1) `t = t̄ + t_h`. ORPHEUS hasn't pursued this as a primary architecture yet.
**Proposed branch**: `feature/peierls-greens-function` off `feature/peierls-specular-bc`.

## Context

**Why this change**. The current ORPHEUS `peierls_*` reference computes the **scalar flux integral form** (Peierls integral equation `Σ_t·φ_i = Σ_j K_ij · q_j`) by Nyström-discretising a kernel that has the angular integration baked in. This works for vacuum BC (Sanchez Eq. 5 `ḡ_2(ρ' → ρ) = (ρ'/ρ)[E_1(σ|r-r'|) − E_1(σ(r+r'))]/2`) and for white BC (Hébert's `(1−P_ss)^{-1}` scalar correction), but the Phase 4/5 rollout exposed two structural limitations:

1. **For specular BC**, the natural continuous form has a hypersingular kernel at the discrete Nyström diagonal (Phase 5 retreat finding, 2026-04-28). The Phase 4 matrix-Galerkin `K_bc = G·R·(I − T·R)^{-1}·P` form works around this via basis projection, but rank-N gating limits accuracy at thin cells beyond N=3.

2. **For multi-region with strong scatterer in outer region** (Issue #132 Class B MR catastrophe), the rank-N Marshak closure has a normalisation-mismatch pathology — rank-2 sphere 1G/2R fuel-A inner / moderator-B outer gives k_eff = 1.015 vs k_inf = 0.648 (+57% sign flip).

The **Green's function approach** — work directly with `t(r', Ω' → r, Ω)` rather than its angle-integrated reduction — gives a different mathematical handle:

- BC absorbed into the kernel `t_h` (Sanchez 1986 Eq. A5) means there is **no separate K_bc closure**; the closure question dissolves. No rank-N gating. No scalar `(1−P_ss)^{-1}` Hébert factor.
- Volume integration on `dr'·dΩ'` rather than just `dr'` exposes the angular structure explicitly, sidestepping the µ-resolved-primitive singularity that killed Phase 5.
- Sanchez 1986 Eq. (A6) IS the Green's function written down — the Phase 5 failure was applying it via Nyström sampling (wrong discretisation) rather than as an integral-operator power iteration (right discretisation).

This is a **substantial architectural reformulation**, not a quadrature improvement. The deliverable is a new continuous-reference solver path: `peierls_greens` parallel to the existing `peierls_geometry`. Both ship; the existing module remains the production reference for vacuum / white BCs (where it's clean and tested), and the new module becomes the production reference for specular and (potentially) for the rank-N Class B MR catastrophe cases the existing module can't handle correctly.

**What prompted this now**. Three converging signals from the Phase 5 retreat:
1. Sanchez 1986 Appendix gives the closed-form Green's function for sphere with α + β BC — the math is laid out but ORPHEUS has only used it as a Nyström target (failed) rather than as an integral-operator
2. Multiple investigators across the Phase 5 rounds suggested integral-operator power iteration as the structurally correct path
3. The Issue #132 / #100 / #103 white-BC closure plateau is unresolved at rank > 1 — a Green's function reference (which doesn't have rank gating) could close that gap

**Intended outcome**:
- Existing `peierls_*` modules organised into a clean package structure
- New `peierls_greens` package with a Green's function continuous reference for sphere homogeneous specular (the simplest test case with MC ground truth available)
- Cross-verification matrix: Green's function reference vs Phase 4 matrix-Galerkin vs MC at thin sphere
- Sphinx documentation of the Green's function architecture as a complement to the existing Peierls integral form
- Decision point at end of plan: ship as a parallel reference, OR fold back into existing Peierls dispatch, OR keep as research-grade-only

## Plan in two parts

### Part A — Reorganise current Peierls into a folder (1 session, ~3 hours)

The existing Peierls code is monolithic: `peierls_geometry.py` is ~5800 lines, plus separate `peierls_cylinder.py`, `peierls_sphere.py`, `peierls_slab.py` shape modules. Before adding a parallel Green's function approach, organise the existing code into a clean package.

#### A1 — Audit current Peierls module structure

```bash
ls -la orpheus/derivations/peierls_*.py
wc -l orpheus/derivations/peierls_*.py
```

Inventory:
- `peierls_geometry.py` (5800 lines) — unified solver, dispatch, primitives, closures, drivers
- `peierls_cylinder.py`, `peierls_sphere.py`, `peierls_slab.py` — shape-specific Peierls reference cases (Phase-era continuous_cases registries)
- Tests in `tests/derivations/test_peierls_*.py`
- Sphinx in `docs/theory/peierls_nystrom.rst`

#### A2 — Move into package structure

Create `orpheus/derivations/peierls/` package:

```
orpheus/derivations/peierls/
├── __init__.py            # public API re-exports for back-compat
├── geometry.py            # CurvilinearGeometry class + per-shape variants
├── kernels.py             # build_volume_kernel + adaptive variant
├── primitives.py          # P_esc, G_bc, P_ss, T (Phase 4) primitives
├── closures.py            # _build_full_K_per_group dispatch + closure variants
├── drivers.py             # solve_peierls_1g, solve_peierls_mg
├── slab.py                # slab-specific helpers (E_n, tau plumbing)
├── cylinder.py            # cyl-specific helpers
└── sphere.py              # sphere-specific helpers
```

Keep `peierls_cylinder.py`, `peierls_sphere.py`, `peierls_slab.py` as the public continuous_cases registries (they import from `orpheus.derivations.peierls`).

The split from monolithic `peierls_geometry.py` follows existing logical boundaries — kernels, primitives, closures, drivers were always separate sections in the file.

**Back-compat**: keep `from orpheus.derivations.peierls_geometry import ...` working via `peierls_geometry.py` shim that re-exports from the new package. Tests reference the old paths; they should continue passing.

#### A3 — Verify all tests pass

```bash
python -m pytest tests/derivations/test_peierls_*.py -v
python -m pytest tests/cp/ -v -k peierls
sphinx-build -W docs docs/_build/html
```

**Acceptance** (Part A):
- (a) `orpheus/derivations/peierls/` package shipped with logical file split
- (b) `peierls_geometry.py` becomes a thin re-export shim
- (c) All existing tests pass without modification
- (d) Sphinx builds clean
- (e) ~5800-line monolith split into 6-8 files of 200-1500 lines each

### Part B — Green's function reference for sphere homogeneous specular (3-5 sessions)

#### B1 — Literature pull (1 session, ~2 hours)

Dispatch `literature-researcher` for canonical Green's function references in transport theory:

- **Case & Zweifel 1967** _Linear Transport Theory_ — the canonical eigenfunction-expansion derivation
- **Williams 1971** _Mathematical Methods in Particle Transport Theory_ — Fredholm theory + Green's functions
- **Pomraning 1973** _The Equations of Radiation Hydrodynamics_ — radiative-transfer Green's functions
- **Pomraning & Siewert 1982** *J. Quant. Spec. Rad. Transfer 28* (cited by Sanchez 1986 as the precursor) — sphere Green's function with isotropic scattering, vacuum BC. **Pull PDF**.
- **Sanchez 2002** *NSE 140* §III (already locally available) — periodic-trajectory Green's function with BC absorbed
- **Bell & Glasstone 1970** Chapter 5 — integral form fundamentals
- **Lewis & Miller 1984** *Computational Methods of Neutron Transport* — Green's function approaches

Output: memo `.claude/agent-memory/literature-researcher/peierls_greens_function_lit.md` with extracted equations, ORPHEUS notation map, and a recommendation between three architectural variants (see B2).

#### B2 — Choose architectural variant (1 session, ~1 hour design)

Three plausible Green's function architectures:

**Variant α — Sanchez 1986 modified-kernel via integral-operator power iteration**.
The kernel `t = t̄ + t_h` (vacuum + BC-absorbed). Solve as fixed-point:

```
ψ^(n+1)(r, Ω) = ∫_V t(r', Ω' → r, Ω) [Σ_s ψ^(n)(r', Ω') + νΣ_f ψ^(n)(r', Ω')/k] dr' dΩ'
```

Discretise `(r, Ω)` on a grid; sample `t` at PAIRS (not at coincident points) so the diagonal singularity that killed Phase 5 doesn't arise — the `dr'·dΩ'` integration smooths over it. Power-iterate to converge ψ; extract k_eff as the Rayleigh quotient.

This is the most direct path: Sanchez Eq. (A6) for sphere is the closed-form `t_h`; we just use it differently than Phase 5 did.

**Variant β — Surface integral equation (BEM-style)**.
Split unknowns into volume `ψ_V(r, Ω)` and surface `ψ_S(r_b, Ω)`. Two coupled integral equations:

```
ψ_V = K_VV · ψ_V + K_VS · ψ_S
ψ_S = K_SS · ψ_S + K_SV · ψ_V
```

with surface kernel `K_SS` carrying the BC. Eliminate surface unknowns via Schur complement. More invasive than Variant α but mathematically cleaner — separates angular at-surface from volume scattering.

**Variant γ — Eigenfunction expansion (Case-Zweifel)**.
Expand `ψ(r, Ω)` in the discrete + continuous eigenfunctions of the homogeneous transport operator. BC closes the expansion via reflection-symmetry constraints. Analytical for homogeneous sphere; Chandrasekhar-H functions enter for the continuous-spectrum projection. ORPHEUS already has Chandrasekhar-H code (Issue #127) — could reuse.

**Recommendation** (subject to literature-pull findings): **Variant α** as the first prototype. It's the most direct adaptation of Sanchez 1986; reuses the existing reference implementation `compute_K_bc_specular_continuous_mu_sphere`; differs from Phase 5 only in the discretisation strategy (integral-operator-on-trial-function vs Nyström K_ij). The previous Phase 5 retreat memo's lesson — "Sanchez form is correct math, wrong discretisation" — directly motivates Variant α.

If Variant α prototypes successfully, Variants β/γ stay as future research directions for the Issue #132 multi-region case (where Variant α may also struggle for reasons different from Phase 5's specular issue).

#### B3 — SymPy derivation (1 session, ~3 hours)

`derivations/peierls_greens_function.py`: SymPy reproduction of the chosen variant. Specifically (assuming Variant α):

1. Reproduce Sanchez Eq. (A4)–(A7) symbolically (already done in `peierls_specular_continuous_mu.py`).
2. Symbolic operator-action: given trial `ψ_trial(ρ, µ)` (e.g., `ρ^k µ^m`), compute the integral-operator action `(K · ψ_trial)(ρ, µ)` analytically. Verify via Sanchez paper consistency checks.
3. Symbolic vacuum reduction: `α → 0` operator action gives the standard vacuum-sphere integral solver. Cross-check vs the ORPHEUS `closure="vacuum"` + Sanchez Eq. (5) reduction.
4. Symbolic rank-1 cross-check: at trial = isotropic constant, the operator action gives `1/(1 − P_ss)` × constant — recovering Hébert's white BC at the operator level.

#### B4 — Prototype implementation (1-2 sessions, ~4-5 hours)

`orpheus/derivations/peierls/greens_function.py`:

1. Build the (r, Ω) grid via existing `composite_gl_r` for r and Gauss-Legendre on Ω (sphere).
2. Implement `compute_K_action_greens(r_pts, Omega_pts, sigma, R, psi_trial)`: applies the Sanchez Eq. (A6) kernel as an operator on `psi_trial`. Critical: **dr'·dΩ' integration is performed FIRST**, then the kernel is sampled — this smooths the diagonal singularity per Phase 5 retreat lesson.
3. Implement `solve_peierls_greens_specular_sphere(R, sig_t, sig_s, nu_sig_f, *, n_r, n_Omega, max_iter=300, tol=1e-10)`: power iteration on the integral operator. Extract k_eff from the dominant eigenvalue.
4. Smoke test: thin homogeneous sphere R=5, σ_t=0.5, fuel-A-like XS, k_inf=0.20833. Target: |k_eff − k_inf| < 0.05% at moderate (n_r, n_Omega).

#### B5 — Cross-verification matrix (1 session, ~2 hours)

Compare Green's function reference vs:

1. **Phase 4 matrix-Galerkin** at N=1, 2, 3 (where Phase 4 is good) — should agree to 0.1% on thin sphere homogeneous
2. **MC ground truth** (`mc_specular_sphere`) — should agree to within MC statistical error
3. **`closure="white_hebert"`** (rank-1, sphere) — should agree algebraically (rank-1 specular_multibounce ≡ white_hebert; rank-1 Green's function should also)

If all three pass, the Green's function reference is validated.

#### B6 — Documentation + closeout (1 session, ~2 hours)

1. Sphinx subsection in `docs/theory/peierls_nystrom.rst` (or new `docs/theory/peierls_greens.rst` page) — Green's function architecture, Sanchez 1986 derivation, cross-verification table.
2. New `closure="specular_greens"` (or similar) registered in `_build_full_K_per_group` if the prototype works as a continuous reference closure. Alternative: ship as a standalone solver `solve_peierls_greens_*` parallel to `solve_peierls_*`.
3. Agent memory: `.claude/agent-memory/numerics-investigator/peierls_greens_phase1_closeout.md`.
4. Decision point: is this a parallel production reference, fold-back into existing Peierls, or keep research-grade-only?

## Acceptance criteria (overall)

Plan 2 is complete when:

- (a) `orpheus/derivations/peierls/` package shipped with logical file split (Part A)
- (b) Existing tests pass without modification (back-compat preserved via shim)
- (c) Literature pull memo shipped with chosen variant
- (d) SymPy derivation script ships with all algebraic identities verified
- (e) Prototype Green's function solver for thin homogeneous sphere passes MC + Phase 4 cross-checks
- (f) Sphinx documents the new architecture
- (g) Decision made: parallel production / fold back / research-grade
- (h) ~5-8 commits scaled by phase boundaries

## Risks

- **R1 Variant α inherits Phase 5 issues despite the discretisation change**. The smoothing from `dr'·dΩ'` integration may not be enough — the singularity could still bite at coincident grid points. Mitigation: explicit subdivision near `r' ≈ r` in the trial-function action; visibility-cone substitution from Plan 1 if shipped.
- **R2 The (r, Ω) grid blows up cost**. Phase 4 matrix-Galerkin is `O(N_r²)`; integral-operator power iteration is `O(N_r² · N_Omega · N_iter)`. At `N_Omega = 64, N_iter = 30`, that's ~2000× slower per power iteration step. Mitigation: profile early; optimise hot path with numpy/numba; possibly use Krylov acceleration on the integral operator instead of plain power iteration.
- **R3 Multi-region sphere with non-uniform σ breaks Sanchez Eq. (A6)** (cosh closed forms require homogeneous σ). The variant α approach has the same multi-region limitation as Phase 5 did. Mitigation: scope Plan 2 to homogeneous sphere first; multi-region is a separate Phase 2+ effort.
- **R4 Reorganisation breaks downstream consumers**. CP modules import from `peierls_geometry`; tests reference the old paths. Mitigation: thin re-export shim at the old path; run full test suite at every package-restructure step.
- **R5 Variant choice (α/β/γ) wrong**. Literature pull may surface a fourth option that's better. Mitigation: B1 explicitly returns a recommendation; user can override before B2.

## Non-goals

- **Not** replacing the existing Peierls module. The new Green's function reference ships ALONGSIDE; existing module remains production for vacuum / white BC cases.
- **Not** revisiting Phase 5 continuous-µ closure. Phase 5 is permanently retreated; the Green's function approach is structurally different (operator-on-trial vs Nyström K_ij sampling).
- **Not** generalising to anisotropic scattering, multi-group, or non-axisymmetric in this plan. Scope: 1-D homogeneous sphere with isotropic scattering and specular BC.
- **Not** replacing `closure="specular_multibounce"` (Phase 4) as the production specular closure. The Green's function reference is for high-precision continuous reference work, not for production solver use (which has rank gating that's appropriate for production).

## Estimated budget

- **Best case**: 5 sessions. A1-A3 in 1 session; B1-B6 spread across 4 sessions.
- **Mid case**: 7 sessions if R2 (cost blow-up) requires Krylov acceleration or other optimisation.
- **Worst case**: 10+ sessions if Variant α fails and Variants β/γ are needed.

LoC delta:
- Part A: ~50-100 LoC change (mostly file moves; some imports tweaked)
- Part B SymPy script: ~200 LoC
- Part B prototype: ~400-600 LoC
- Tests: ~300 LoC
- Sphinx: ~200 lines
- Memos: ~300 lines

Total: ~1500-2000 LoC + docs. Commit count: 5-8.

## Critical files

### Part A — reorganisation

- `orpheus/derivations/peierls_geometry.py` — current monolith (5800 lines); becomes shim
- `orpheus/derivations/peierls_cylinder.py`, `peierls_sphere.py`, `peierls_slab.py` — keep as continuous_cases registries; update imports
- `orpheus/derivations/peierls/` — new package
- `tests/derivations/test_peierls_*.py` — should not need modification (back-compat shim)

### Part B — Green's function

- `derivations/peierls_greens_function.py` — SymPy derivation (NEW)
- `orpheus/derivations/peierls/greens_function.py` — prototype solver (NEW)
- `tests/derivations/test_peierls_greens_function.py` — cross-verification (NEW)
- `docs/theory/peierls_nystrom.rst` or `docs/theory/peierls_greens.rst` — architecture doc
- `.claude/agent-memory/literature-researcher/peierls_greens_function_lit.md` — literature pull (NEW)
- `.claude/agent-memory/numerics-investigator/peierls_greens_phase1_closeout.md` — closeout (NEW)

### Reusable from existing code

- `compute_K_bc_specular_continuous_mu_sphere` (Phase 5a, peierls_geometry.py) — reuse the Sanchez Eq. (A6) implementation
- `derivations/peierls_specular_continuous_mu.py` — reuse the SymPy 4/4 verifications
- `mc_specular_sphere` (`derivations/diagnostics/diag_specular_overshoot_05_mc_multibounce.py`) — MC ground truth
- `composite_gl_r`, `gl_float`, `gl_nodes_weights` — quadrature primitives

### References

- `.claude/agent-memory/literature-researcher/phase5_sanchez_1986_sphere_specular.md` — Sanchez 1986 PDF extraction (Phase 5a)
- `.claude/agent-memory/numerics-investigator/specular_continuous_mu_phase5_retreat.md` — Phase 5 retreat (motivates Variant α discretisation choice)
- `Sanchez(1986)TTSP14.pdf` — local PDF
- `Sanchez(2002).pdf` — local PDF
- `Hebert(2009)Chapter3.pdf` — local PDF for cross-reference
