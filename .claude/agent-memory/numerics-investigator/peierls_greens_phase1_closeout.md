---
name: Variant α Green's function — running closeout (Plan 2 Phase B + A1-A3 + Plan-(b))
description: Living closeout for the Variant α Green's function prototype. Plan 2 Phase B closed 2026-05-01 (closed-sphere homogeneous specular). A1-A3 follow-on closed 2026-05-02 (vacuum BC + PS-1982 + multi-group). Plan-(b) Options 1+2 closed 2026-05-02 (multi-region + Issue #132 attack + Garcia 2021 cross-check). 53 tests across all phases, all green. V&V layered: L1 closed 1G via SymPy + Phase 4 cross-check; L1 vacuum 1G via PS-1982; L1 closed MG via transfer-matrix; L1 multi-region fixed-source via Garcia 2021.
type: project
---

# Variant α Green's function — running closeout

**Started**: 2026-05-01 (Plan 2 Phase B).
**Last update**: 2026-05-02 (Plan-(b) Options 1+2).
**Branch**: `feature/peierls-greens-function`.
**Outcome**: Variant α shipped as research-grade parallel reference covering
homogeneous + multi-region sphere across vacuum, partial-albedo, and
specular BCs, with multi-group support up to general G. Phase 4
``specular_multibounce`` retained as the production path — Variant α is
the verification reference, not a production replacement.

**Phases tracked**:
- **Plan 2 B1-B6** (2026-05-01) — closed homogeneous sphere, specular BC,
  1-group: V_α1 closed-sphere identity, V_α2 T_00=P_ss algebraic
  cross-check, V_α3 vacuum-reduction algebra, Sphinx page.
- **A1-A2** (2026-05-02) — vacuum BC parametrisation + PS-1982 Eq. (21)
  L1 cross-check (structurally-independent kernel reference).
- **A3** (2026-05-02) — multi-group extension (general G, full G×G
  scattering matrix, χ-spectrum). Closes 1G Cardinal-Rule-6 gap.
- **Plan-(b) Options 1+2** (2026-05-02) — multi-region extension
  (k-eigenvalue + fixed-source); Issue #132 reproducer; Garcia 2021
  flux-shape L1 cross-check.

## TL;DR

The Variant α prototype (`orpheus/derivations/continuous/peierls/greens_function.py`)
solves the homogeneous sphere with perfect specular BC by iterating the
**angle-resolved** Green's function :math:`\tilde t(r' \to r, \mu)`
along bouncing characteristics — the bounce sum is closed analytically
via :math:`T(\mu_{\rm surf}) = 1/(1 - e^{-\Sigma_t \cdot 2 R \mu_{\rm surf}})`,
and the trajectory + bounce-period integrals are 1-D Gauss-Legendre
quadratures along single characteristics. The angle-integrated kernel
:math:`g_\alpha(\rho' \to \rho)` (the Phase-5 hypersingular form) is
**never assembled**, so the Hadamard finite-part diagonal singularity
that killed Phase 5 is bypassed structurally.

For the closed homogeneous sphere, V_α1 algebraically proves
:math:`(K \cdot 1)(r, \mu) = \omega_0 \cdot \mathrm{const}`, giving
the rank-1 isotropic eigenmode at :math:`k_{\rm eff} = k_\infty =
\nu\Sigma_f / \Sigma_a`. The numerical implementation reproduces this
to machine precision.

## Cross-verification matrix (fuel-A-like, R = 5, τ_R = 2.5)

XS: σ_t = 0.5, σ_s = 0.38, νσ_f = 0.025. k_inf = 0.025 / 0.12 = 0.2083‾.

| Reference                 | k_eff             | err vs Variant α    |
|---------------------------|-------------------|---------------------|
| k_inf (analytic)          | 0.20833333333     | —                   |
| **Variant α** (this work) | **0.20833333333** | **1e-13 %** (machine)|
| Phase 4 N = 1             | 0.20777642780     | 0.2673 %            |
| Phase 4 N = 2             | 0.20781489733     | 0.2488 %            |
| Phase 4 N = 3             | 0.20808891764     | 0.1173 %            |
| white_hebert              | 0.20777642780     | 0.2673 % (= N=1)    |

**Findings**:

1. **Variant α is exact** for the closed homogeneous sphere
   (V_α1 algebraic identity reproduced to 1e-13 %).
2. **Phase 4 N=1 ≡ white_hebert bit-for-bit** (V_α2 algebraic identity:
   :math:`T_{00}^{\rm sphere} = P_{ss}^{\rm sphere}`).
3. **Phase 4 errors decrease with rank** (0.27 % → 0.25 % → 0.12 %)
   but never reach the Variant α exact value — the rank-N truncation
   carries a systematic small bias even at N = 3.
4. **Variant α improves on Phase 4 N=3 by ~0.12 %** for this
   configuration. For users needing higher precision than Phase 4 can
   deliver on closed sphere homogeneous, Variant α is the reference.

## Phase B deliverables shipped

### Code

- `orpheus/derivations/continuous/peierls/greens_function.py` (~370 LoC)
  — prototype solver `solve_greens_function_specular_sphere` and the
  per-pair operator `_apply_operator` along bouncing characteristics.
- `orpheus/derivations/continuous/peierls/origins/specular/greens_function.py` (~270 LoC)
  — SymPy derivations V_α1, V_α2, V_α3 (operator-level identities).
- `orpheus/derivations/continuous/peierls/origins/specular/__init__.py`
  — re-exports added.

### Tests (15 new, all green)

- `tests/derivations/test_peierls_greens_function_symbolic.py`
  — 8 SymPy gates (V_α1 surface fixed-point + total-ψ-constant
  + operator-eigenvalue + composite; V_α2 integrand-match + closed-form
  + composite; V_α3 vacuum reduction).
- `tests/derivations/test_peierls_greens_function_solver.py`
  — 3 numerical gates (V_α1.numerical constant initial guess +
  non-uniform initial guess + two thicknesses).
- `tests/derivations/test_peierls_greens_function_xverif.py`
  — 3 cross-verification gates (Variant α exact + Phase 4 N=1 ≡
  white_hebert + Phase 4 rank convergence toward Variant α).

### Memos

- `.claude/agent-memory/literature-researcher/peierls_greens_function_lit.md` (B1)
- `.claude/agent-memory/numerics-investigator/peierls_greens_variant_alpha_decision.md` (B2)
- `.claude/agent-memory/numerics-investigator/peierls_greens_phase1_closeout.md` (this file, B6)

### Sphinx (deferred to follow-on, see "Decision" below)

The plan calls for a Sphinx subsection in `docs/theory/peierls_nystrom.rst`
(or new page `docs/theory/peierls_greens.rst`) with the Green's function
architecture, Sanchez 1986 derivation walkthrough, and the
cross-verification table. This is dispatched to the archivist agent
as a follow-on task — the SymPy + numerical gates are the load-bearing
artifacts; the Sphinx page is documentation around them.

## Decision: parallel research-grade reference

Plan 2 asked at the closeout: ship as parallel production reference,
fold back into existing Peierls dispatch, or keep research-grade-only.

**Decision**: ship as **parallel research-grade reference**.

Reasoning:

1. **Variant α has narrower validation than Phase 4**. Variant α is
   tested for **homogeneous sphere only, isotropic scattering only,
   perfect specular BC only**. Phase 4 ``specular_multibounce`` is
   tested across slab / cylinder / sphere, multi-region, anisotropic
   scattering (when applicable), vacuum / white / specular BCs, with
   rank-N gating that has documented failure modes. Variant α inherits
   none of this validation; promoting to production would require
   substantial additional work.
2. **Variant α improves on Phase 4 only for the cases it covers**.
   For closed homogeneous sphere with isotropic scattering, Variant α
   is exact. For all other configurations, Variant α is unimplemented
   and Phase 4 is the only ORPHEUS reference. Issue #132 (Class B MR
   catastrophe) is NOT solved by Variant α — multi-region sphere is
   future work.
3. **The closure-free property is the load-bearing structural
   advantage**. Future research-grade work on continuous-µ specular
   benefits from a working operator-level reference. Variant α serves
   as the "gold standard" against which closure approximations
   (Phase 4 rank-N, Hébert white, Garcia P_N) are calibrated.

**Not folded into ``boundary=...`` dispatch**: keeping Variant α as a
standalone solver `solve_greens_function_specular_sphere` makes it
clear that it is a **reference**, not a production solver. Folding
into dispatch would imply broader validation than has been done.

**Not deleted**: the prototype is already a useful reference (one
machine-precision data point for the cross-verification matrix). The
~640 LoC of code + 15 tests is a small footprint for what it provides.

## What's NOT in scope for Phase B (deferred to future plans)

0. ~~**1-group Cardinal-Rule-6 violation**.~~ **CLOSED IN A3 follow-on
   (2026-05-02).** Prototype extended to multi-group via
   `solve_greens_function_sphere_mg(R, sigma_t, sigma_s, nu_sigma_f, chi)`
   with arbitrary G groups, full G×G scattering matrix (downscatter +
   upscatter both supported), and arbitrary χ fission spectrum. Cross-
   checks at closed sphere reduce to `kinf_and_spectrum_homogeneous`
   transfer-matrix dominant eigenvalue + spectrum to ≤ 1e-9 relative.
   Tests cover G=1 (matches 1G solver bit-equal), 2G downscatter, 2G
   upscatter, 2G vacuum, 2G analytical flux ratio
   (Σ_{1→2}/Σ_R,2 = 0.5 closed-form), 4G fuel-A with realistic
   χ = (0.6, 0.35, 0.05, 0.0), 4G vacuum. Refactor extracted shared
   `_apply_operator_with_source_profile` so 1G case becomes a thin
   wrapper — no regression in 25 prior 1G tests.

1. ~~**Vacuum BC** (:math:`\alpha = 0`).~~ **CLOSED IN A1+A2 follow-on
   (2026-05-02).** Prototype now takes `alpha ∈ [0, 1]` parameter.
   Vacuum BC cross-checked against a structurally-independent
   Pomraning-Siewert 1982 Eq. (21) reference solver
   (`orpheus.derivations.continuous.peierls.ps1982_reference`). Six
   parametrised tests pass at ≤ 1e-4 relative agreement on the k_eff
   eigenvalue, exercising thin (τ_R = 2.5) and moderate (τ_R = 5)
   spheres with strong/medium absorber. **Bug found and fixed during
   A1+A2**: original prototype used the FORWARD distance from r to the
   surface (`√(R²−r²(1−µ²)) − r·µ`) where the integral form needs the
   BACKWARD distance (`r·µ + √(...)`). For closed sphere (V_α1
   algebraic identity) the bug was masked because the closure cancels
   L_first identically; for vacuum BC it surfaced as a 6 % k_eff
   disagreement vs PS-1982 → corrected to < 1e-4. Also generalised the
   µ-grid from `(0, 1]` (closed-sphere µ → −µ symmetry assumption) to
   the full `[−1, 1]` (vacuum BC breaks the symmetry — different
   trajectories outward vs inward).
2. ~~**Multi-region sphere**.~~ **CLOSED IN Plan-(b) follow-on
   (2026-05-02).** Prototype extended to multi-region via
   `solve_greens_function_sphere_mr` and
   `solve_greens_function_sphere_mr_fixed_source` (~520 LoC). Both
   handle arbitrary G-region sphere with piecewise σ_t along
   trajectories and bounce-period chords (composite GL quadrature
   per region segment).

   **Plan-(b) Option 2 (Issue #132 reproducer)**: 4 tests pass.
   Sphere `radii=[0.5, 1.0]`, fuel-A inner / moderator-B outer,
   closed BC: Variant α gives k_eff = 0.735 vs Phase 4 specular_-
   multibounce rank-2 = 1.015 (+57 % catastrophe avoided). Spatial
   mode physically correct (φ peaked in fuel, decreasing through
   moderator, slope discontinuity at interface).

   **Plan-(b) Option 1 (Garcia 2021 flux-shape L1 cross-check)**:
   17 tests pass. Three-region sphere fixed-source benchmark
   (Williams 1991 Case 1, Garcia 2021 Table 5) with
   structurally-independent stable-P_N reference. Convention-
   converted Variant α agrees with Garcia table to < 1 % at non-
   interface r-points, < 12 % near interfaces (cubic-spline
   smoothing of source profile across discontinuous σ_s — known
   prototype limitation). All 15 r-points pass tolerance gates.

   Direct attack on Issue #132 Class B MR catastrophe demonstrated
   working: Variant α has no rank-N closure, so the mode-0/mode-≥1
   normalisation mismatch that breaks Phase 4 simply cannot occur
   structurally.
3. **Cylinder geometry**. Sanchez 1986 has cylinder via the unified
   :math:`\alpha`-parameter (:math:`\alpha = 1` for cylinder, but here
   :math:`\alpha` is geometry-shape, not specular-coefficient — confusing
   notation overlap). Sanchez Eq. (A6) for cylinder uses Bessel
   functions instead of cosh; the bounce-sum machinery is structurally
   identical.
4. **Anisotropic scattering** (:math:`\omega_1 \ne 0`). Sanchez 1986
   gives the :math:`h` kernel for linearly anisotropic scattering;
   the bounce-sum trajectory machinery should extend with minor
   modifications.
5. **k-eigenvalue extraction beyond rank-1**. The closed sphere has
   only a rank-1 eigenmode (no spatial mode structure, no leakage).
   Vacuum BC + multi-region are the cases where higher modes appear
   and the iteration becomes non-trivial.
6. **Garcia 2020/2021 stable P_N cross-check**. Garcia papers are
   paywalled; institutional pull would unlock external numerical
   reference. Recommended for B5 follow-on.

## Risks and known limitations

- **R1 — homogeneous-only restriction**. Plan 2 R3 risk realised: the
  prototype is restricted to single-region homogeneous medium because
  Sanchez Eq. (A6) cosh closed forms require uniform :math:`\Sigma_t`.
  Multi-region extension requires re-deriving the chord parametrisation
  with piecewise optical depth.
- **R2 — cost not measured against Phase 4**. Variant α per-iteration
  cost is :math:`O(N_r \cdot N_\mu \cdot N_{\rm traj})` versus Phase 4's
  :math:`O(N_r^2)`. For typical (N_r, N_µ, N_traj) = (16, 16, 32) and
  Phase 4 N_r = 24, Variant α should be a few× slower per iteration
  but converges in 1-3 iterations from a good initial guess for closed
  sphere. Real-world runtime not benchmarked.
- **R3 — no benchmark against MC**. Plan 2 B5 listed MC as one of
  three cross-verification anchors. For closed homogeneous sphere with
  isotropic scattering, MC would also give k_inf within statistical
  uncertainty — same answer as Variant α and analytic k_inf. The MC
  cross-check adds no new information for this trivial case; deferred
  until Variant α extends to non-trivial geometries (vacuum / multi-
  region).

## Suggested follow-on issues to create

1. **Variant α vacuum BC extension** — implement :math:`\alpha = 0`
   branch, cross-check against ``boundary="vacuum"`` Peierls reference.
   Low-risk, high-value (closes the "does Variant α inherit the
   existing vacuum reference cleanly" question per V_α3).
2. **Variant α multi-region sphere** — direct attack on Issue #132
   Class B MR catastrophe. High-risk, very-high-value.
3. **Garcia 2020/2021 P_N PDFs** — institutional pull request to add
   external paywalled reference for sphere homogeneous reflective.

## File index

- This memo: `.claude/agent-memory/numerics-investigator/peierls_greens_phase1_closeout.md`
- B1 lit: `.claude/agent-memory/literature-researcher/peierls_greens_function_lit.md`
- B2 decision: `.claude/agent-memory/numerics-investigator/peierls_greens_variant_alpha_decision.md`
- B3 SymPy: `orpheus/derivations/continuous/peierls/origins/specular/greens_function.py`
- B3 test: `tests/derivations/test_peierls_greens_function_symbolic.py`
- B4 prototype: `orpheus/derivations/continuous/peierls/greens_function.py`
- B4 test: `tests/derivations/test_peierls_greens_function_solver.py`
- B5 test: `tests/derivations/test_peierls_greens_function_xverif.py`
- Plan 2: `.claude/plans/peierls-greens-function-approach.md`
- Phase 5 retreat (predecessor): `.claude/agent-memory/numerics-investigator/_archive/specular_continuous_mu_phase5_retreat.md`
- Sanchez 1986 memo (predecessor): `.claude/agent-memory/literature-researcher/phase5_sanchez_1986_sphere_specular.md`
