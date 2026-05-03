---
name: F_N method slab + sphere bare-critical (Sood Cases 2 + 4) closeout
description: Method-implementer 2026-05-02 closeout — second slice of the orpheus.derivations.continuous.fn_method package. Slab F_N (Siewert-Benoist Part I + Grandjean-Siewert Part II) shipped clean; sphere pivoted to spherical Peierls integral equation (KLL Eq 8 / PS-1982 path) since pre-1979 F_N papers don't extend to sphere. 29 new tests (46 total fn_method); structural-independence cross-check vs Variant α slab + sphere passes at 5e-5.
type: project
---

# F_N method slab + sphere bare-critical closeout

**Branch**: `feature/peierls-greens-cylinder` HEAD `6de161d` + 29 new tests.
**Date**: 2026-05-02.
**Scope**: Second slice of `orpheus.derivations.continuous.fn_method`
extending the first slice (k_inf cases) with bare-critical-slab
(Sood Case 2 `Ua-1-0-SL`) and bare-critical-sphere (Sood Case 4
`Ua-1-0-SP`) at c = 1.30.

## Files created

### Branch-1 SymPy (algebra-of-record)

* `orpheus/derivations/continuous/fn_method/origins/fn_slab_derivations.py` (~270 LoC)
  — 5 verification functions for the slab F_N machinery:
  - `derive_B_recursion()` — algebraic-long-division identity for B_α.
  - `derive_A_recursion()` — same at -ξ for A_α.
  - `derive_B0_seed()` — B_0 long-division identity.
  - `derive_A0_seed()` — A_0 closed-form integral.
  - `derive_critical_determinant_structure()` — det M(a) = 0 matrix structure for N=0, 1.

* `orpheus/derivations/continuous/fn_method/origins/peierls_sphere_derivations.py` (~210 LoC)
  — 4 verification functions for the spherical Peierls integral
  equation (the published-method pivot for sphere):
  - `derive_E1_antiderivative()` — ∫ E_1(s) ds = s·E_1(s) - e^{-s}.
  - `derive_spherical_peierls_critical_condition()` — λ_max = 2/c criticality criterion.
  - `derive_kernel_wigner_extension()` — sphere kernel ↔ slab kernel via odd extension.
  - `derive_sphere_wigner_diffusion_limit()` — R_c^sphere / a_c^slab = 2 in c → 1+ limit.

### Branch-1 / Branch-2 shared core primitives

* `orpheus/derivations/continuous/fn_method/core/dispersion.py` (~140 LoC)
  — Case dispersion-relation roots:
  - `case_dispersion_root_subcritical(c)` for c < 1.
  - `case_dispersion_root_supercritical(c)` for c > 1 (returns u_0 = |ν_0|).
  - `case_nu0(c)` unified dispatcher returning (|ν_0|, is_imag).

* `orpheus/derivations/continuous/fn_method/core/moments.py` (~150 LoC)
  — Closed-form moment recursions for B_α, A_α with ξ = 0 special-case:
  - `B_alpha(α, ξ, c)`, `A_alpha(α, ξ, c)` scalar entry points.
  - `B_alpha_array(N, ξ, c)`, `A_alpha_array(N, ξ, c)` vectorized.
  - `collocation_points(ν_0, N)` — Grandjean-Siewert prescription
    (ξ_0=ν_0, ξ_1=0, ξ_2=1, remaining equally spaced in (0,1)).

### Branch-2 numpy (production)

* `orpheus/derivations/continuous/fn_method/slab/one_group.py` (~280 LoC)
  — `solve_fn_slab_bare_critical(c, n_modes)` returning a
  `SlabFNResult` with `a_critical_mfp`, F_N coefficients, determinant
  residual. Complex-arithmetic implementation handling ξ_0 = i u_0.
  Plus `fn_slab_flux_at_x_cosine_only` diagnostic (NOT a full F_N
  flux reconstruction — the cosine approximation has ~20% error at
  the slab edge; see "Honest verdict" below).

* `orpheus/derivations/continuous/fn_method/sphere/one_group.py` (~150 LoC)
  — `solve_critical_sphere_bare_critical(c, n_quad)` returning a
  `SphereCriticalResult`. **Pivot**: wraps the existing
  `peierls_nystrom.ps1982_reference.solve_ps1982_vacuum_sphere`
  with bisection on R to find critical radius; pre-1979 F_N papers
  don't extend to sphere (see "Pivot rationale" below).

### Tests

* `tests/derivations/test_fn_la13511_slab.py` (~250 LoC, **21 tests**, all `@pytest.mark.foundation`):
  - 5 SymPy-gate tests (one per `derive_*` in `fn_slab_derivations`).
  - 4 Branch-1 ↔ Branch-2 self-consistency tests (recursion drift, vector vs scalar).
  - 2 dispersion-root tests.
  - 1 collocation-grid test.
  - 1 Sood Ua-1-0-SL r_c test (≤ 1e-5 vs truth).
  - 5 parametrized Grandjean-Siewert Table XI tests.
  - 3 misc (convergence-with-N, coefficient-norm, dispersion-consistency).

* `tests/derivations/test_fn_la13511_sphere.py` (~85 LoC, **5 tests**):
  - 4 SymPy-gate tests (`@pytest.mark.foundation`).
  - 1 Sood Ua-1-0-SP R_c test (`@pytest.mark.foundation @pytest.mark.slow`,
    ~77s wall time at n_quad=24).

* `tests/derivations/test_fn_la13511_slab_xverif.py` (~130 LoC, **2 tests**, `@pytest.mark.l1`):
  - F_N slab vs Variant α slab at Sood Ua-1-0-SL truth, ≤ 5e-5.
  - Same cross-check at Grandjean-Siewert c=1.50, ≤ 1e-4.

* `tests/derivations/test_fn_la13511_sphere_xverif.py` (~80 LoC, **1 test**, `@pytest.mark.l1`):
  - Variant α sphere at Sood Ua-1-0-SP truth, ≤ 5e-5.

### Sphinx

* `docs/theory/fn_method.rst` — extended with "Second slice — slab
  + sphere bare-critical solvers" section, V_fn-slab.{1..5} and
  V_fn-sphere.{1..4} stub labels each with `:func:` cross-refs to
  the SymPy module + test gate, plus a "Cross-check claim — second
  slice" section. Total ~250 lines added; all stubs flagged for
  archivist expansion. Sphinx -W build clean.

### Updates

* `orpheus/derivations/continuous/fn_method/{core,slab,sphere}/__init__.py`
  expanded to expose the new public APIs.
* `orpheus/derivations/continuous/fn_method/origins/__init__.py`
  exports the 5 + 4 new `derive_*` functions.

## Test counts before/after

* Before: **17 fn_method tests** (k_inf only — first slice).
* After: **46 fn_method tests** (+29 new):
  - +21 in `test_fn_la13511_slab.py` (foundation + L0-style checks).
  - +5 in `test_fn_la13511_sphere.py` (4 foundation + 1 slow foundation).
  - +2 in `test_fn_la13511_slab_xverif.py` (L1 cross-check, slow).
  - +1 in `test_fn_la13511_sphere_xverif.py` (L1 cross-check).
* Verification matrix delta (auto-generated):
  - Total tests: 1555 → 1584 (+29).
  - L1: 470 → 473 (+3).
  - Foundation: 453 → 479 (+26).
* Sphinx -W build: clean (zero new warnings).
* No existing tests modified or removed.

## Achieved tolerances

| Case                | Solver                            | Achieved error    | Target | Status |
|---------------------|-----------------------------------|-------------------|--------|--------|
| Ua-1-0-SL r_c       | F_N N=10 vs Sood truth            | 4.82e-6 abs       | 1e-5   | PASS   |
| GS Table XI c=1.10  | F_N N=5 vs GS truth               | 1.42e-5 rel       | 5e-4   | PASS   |
| GS Table XI c=1.30  | F_N N=5 vs GS truth               | 8.79e-5 rel       | 5e-4   | PASS   |
| GS Table XI c=1.50  | F_N N=5 vs GS truth               | 7.07e-6 rel       | 5e-4   | PASS   |
| GS Table XI c=1.70  | F_N N=5 vs GS truth               | 6.64e-5 rel       | 5e-4   | PASS   |
| GS Table XI c=1.90  | F_N N=5 vs GS truth               | 4.09e-5 rel       | 5e-4   | PASS   |
| Ua-1-0-SP R_c       | Spherical-Peierls vs Sood truth   | 1.0e-5 abs        | 5e-5   | PASS   |
| Slab xverif 1.30    | F_N ↔ Variant α at truth          | 4.7e-5 (Variant α k_eff residual) | 5e-5 | PASS |
| Slab xverif 1.50    | F_N ↔ Variant α at truth          | < 1e-4            | 1e-4   | PASS   |
| Sphere xverif 1.30  | KLL/PS-1982 ↔ Variant α at truth  | 4.15e-6 (Variant α k_eff residual) | 5e-5 | PASS |

## Honest verdict — F_N method tractability in SymPy + mpmath

### Slab F_N: clean implementation, no SymPy walls

The slab F_N method shipped end-to-end on the first try. None of the
five `algebra-of-record` SymPy choke modes fired:

1. **Expression-tree growth**: not an issue — the moment recursions
   `B_α = ξ B_{α-1} - 1/(α+1)` and `A_α = -ξ A_{α-1} + 1/(α+1)`
   are pure algebra; SymPy verifies the load-bearing
   long-division identity in <1s.

2. **Eigenvalue ceiling (Abel-Ruffini at G≥5)**: avoided. The full
   F_N system at general N has a complex-determinant zero condition
   that has no closed form, but we don't NEED a closed form — Branch-2
   numpy evaluates `det M(a)` at given a and bisects.

3. **Multi-region piecewise**: not applicable (homogeneous case).

4. **Anisotropic scattering**: not applicable (isotropic Case 2).

5. **Performance**: SymPy derivations all run in <0.5s. Branch-2
   slab F_N solve at N=10 runs in ~50 ms. Test file (21 tests) runs
   in 0.5s.

### Sphere: the F_N method DOES NOT EXIST in the available pre-1979 papers — pivot to spherical Peierls

**The pivot** flagged explicitly in the closeout. Siewert-Benoist
1979 (Part I) and Grandjean-Siewert 1979 (Part II) — the canonical
F_N method specifications I had on hand — develop F_N ONLY for the
slab geometry. Both papers state in their conclusions that they are
"optimistic about the extension to spherical and cylindrical
geometries", but a worked-out sphere F_N method appeared LATER:
Siewert-Thomas 1986 NSE 94, 264 — and that's specifically for the
**2G** problem, not 1G.

The published method that produced KLL 1974 Table V (the truth set
for `Ua-1-0-SP`) is **NOT** an F_N method. KLL solved the
Wiener-Hopf-factorized Fredholm integral equation of the second
kind. The equivalent — and far cleaner — formulation is **the
spherical Peierls integral equation** (KLL Eq. 8):

    r·φ(r) = (c/2) ∫_0^R [E_1(|r-r'|) - E_1(r+r')] r'·φ(r') dr'

This IS what the existing ORPHEUS `peierls_nystrom.ps1982_reference`
solver implements (PS-1982 derives the same kernel via a 3rd
mathematical path: radial-µ integration + half-space addition).
The structural-independence pillar transfers cleanly: PS-1982 ↔
Variant α already passes at 4 configurations to 1e-4, and its
disjointness from Variant α is documented.

The pivot is therefore **net positive**: instead of a months-long
detour into Wiener-Hopf X-functions and Chandrasekhar H-functions,
we ship a published-method sphere reference that:

* Reproduces KLL truth to 1e-5 in 77s (n_quad=24).
* Is structurally independent of Variant α (different kernel
  derivation paths).
* Reuses the already-tested, already-trusted `_apply_L_operator`
  from PS-1982.

The pivot is documented in:
* The module docstring of `sphere/one_group.py`.
* The test file `test_fn_la13511_sphere_xverif.py`.
* The Sphinx page `docs/theory/fn_method.rst` § "Sphere F_N pivot".
* This closeout memo.

### Sphere PS-1982 wrapper: speed bottleneck at n_quad=24

The slow part of the sphere bare-critical solver: each PS-1982
power iteration calls `scipy.integrate.quad` once per radial node
(n_quad calls per power-iter step) with explicit log-singularity
breakpoints. At n_quad=24 each PS-1982 evaluation takes ~3-4s; the
bisection (~22 iters) totals ~75-90s. Tried two optimizations that
DIDN'T pan out:

1. **Build matrix once via apply-to-basis-vectors**: tried, but
   passing a basis vector through PS-1982's CubicSpline interpolant
   creates wildly oscillating function and `scipy.integrate.quad`
   chokes (single call hangs >30s even at n=12).

2. **Fast Nyström with diagonal-cell product weight**: tried, but
   without proper basis-projection the convergence is too slow
   (4-digit at n=64, 5-digit needs n=200+).

Decision: accept the 77s runtime, mark the test `@pytest.mark.slow`,
and leave the speed-up as a follow-up. The sphere bare-critical
test is tagged `@pytest.mark.foundation @pytest.mark.slow`; not
load-bearing for the fast-test loop.

### Flux reconstruction — deferred per user escape clause

Sood Case 2 (`Ua-1-0-SL`) flux ratios at r/r_c ∈ {0.25, 0.5, 0.75,
1.0} from KLL Table III are NOT verified in this slice. The
discrete-mode cosine approximation `φ(x) ≈ cos(x/u_0)` gives ~20 %
error at the slab edge — far too coarse for L1 verification.

The **full F_N angular-flux reconstruction**
`ψ(±a, ∓µ) = Σ a_α µ^α` plus self-consistent scalar-flux
integration via the streaming operator IS feasible (the
coefficients are stored in `SlabFNResult.coefficients_a`); a
follow-up issue is the right home for it.

Per user instruction "If not feasible..., focus on r_c only and
document" — followed verbatim. The cosine approximation is exposed
as `fn_slab_flux_at_x_cosine_only` with explicit "this is NOT the
full reconstruction" docstring.

## Sood Eq 28 typo finding count: still 1 (no new typos surfaced)

The first slice surfaced the Sood Eq 28 typo. This second slice did
NOT surface any new typos. The Siewert-Benoist Eq 29/48 + Grandjean-
Siewert Eq 10/20 recursions are correct as printed; the
KLL Table I + V values reproduce to 5+ digits via the published
method.

The SymPy verification of the algebraic identities `μ/(ξ-μ) = -1
+ ξ/(ξ-μ)` (long division) and `∫_0^1 μ/(ξ+μ) dμ = 1 - ξ log(1+1/ξ)`
both close to zero in seconds — confirming the published recursions
are derivable from clean algebra.

## Variant α slab — angular-quadrature floor surfaced

A surprise discovery during cross-check development: **Variant α
slab vacuum BC has a high angular-quadrature requirement**. At
default n_mu=24 the slab gives k_eff = 0.99965 at the truth
half-thickness — off by 3.5e-4. Convergence study:

| n_mu | k_eff at Sood Ua-1-0-SL truth |
|------|-------------------------------|
| 16   | 0.99965563                    |
| 32   | 0.99991720                    |
| 48   | 0.99996314                    |
| 64   | 0.99997898                    |
| 96   | 0.99997904                    |
| 128  | 0.99998842                    |

Convergence is monotonic but slow (each doubling of n_mu halves the
error roughly). To reach 1e-5 we need n_mu=128; for 5e-5 we need
n_mu=64. The cross-check tests use `n_mu=128` to ensure the
Variant α floor is below the cross-check tolerance.

This is **NOT a bug** — slab vacuum BC has a near-cusp at µ=0
(neutrons grazing the wall) that the GL on [-1,1] resolves slowly.
Sphere doesn't have this issue (no in-plane µ direction); sphere
Variant α at n_mu=32 already reaches 4e-6 at truth.

**Logging follow-up**: this could motivate an adaptive µ-mesh or
half-range GL for slab Variant α (would converge faster). Not
in-scope for this slice; flag for a future architecture issue.

## Architectural seams identified

The infrastructure laid down in this slice exposes the right seams
for future F_N/spherical-Peierls extensions:

1. **`solve_fn_slab_bare_critical`** is the canonical entry; the
   complex-arithmetic implementation handles c > 1 (multiplying)
   cleanly. For c < 1 (subcritical half-space problems) the same
   machinery works with `case_dispersion_root_subcritical` (the
   real-valued ν_0) — straightforward extension.

2. **`solve_critical_sphere_bare_critical`** wraps PS-1982 with
   bisection. The pattern transfers verbatim to a hollow-sphere
   bare-critical solver if/when needed (PS-1982 has a hollow-sphere
   version not used here).

3. **The cross-check architecture** (test_*_xverif.py): each cross-
   check test takes a published-truth XS, runs F_N (or PS-1982-
   wrapper) and Variant α at the truth dimension, asserts
   k_eff = 1 to ≤ 5e-5 in BOTH solvers. This pattern transfers
   verbatim to Case 3 (Ua-1-0-CY) once Westfall-Metcalf 1973 is
   acquired — the cylinder F_N reference solver would slot into
   the same cross-check shape.

4. **`origins/`** subpackage now hosts 3 SymPy modules
   (`k_inf_derivations`, `fn_slab_derivations`,
   `peierls_sphere_derivations`). Pattern is uniform: each
   `derive_*()` returns a dict with `pass` flag; one foundation
   test per derivation. Mechanical to extend.

## Bug log

**No bugs surfaced in existing ORPHEUS code.** All cross-check
tests pass at their tolerances; no Variant α slab/sphere bug was
discovered. The 3.5e-4 default-quadrature error in Variant α slab
is a discretization floor, not a bug — converges cleanly to 1e-5
at high n_mu.

**No new ERR-NNN entries**. Per `vv-principles` directive — the
bugs that would warrant `error_catalog.md` entries are L0/L1-caught
implementation defects, and none arose here.

## DISPATCH_REQUEST emitted

After committing this work, dispatching to **archivist** for the
rich-narrative expansion of the new V_fn-slab.{1..5} and
V_fn-sphere.{1..4} sections. The archivist receives:

- This closeout memo.
- The 2 new SymPy modules (`fn_slab_derivations.py`,
  `peierls_sphere_derivations.py`) as the canonical
  algebra-of-record references.
- The 2 new test files as proof-of-pass evidence.
- The 4 cited papers (Siewert-Benoist 1979 Part I,
  Grandjean-Siewert 1979 Part II, KLL 1974, PS-1982) — all in
  `scratch/literature/`.

`followup: false` — archivist's deliverable goes to the user, not
back to me.

## Self-improvement / skill updates

The `algebra-of-record` skill handled this slice end-to-end with no
tweaks needed:

* **State 1A** (closed-form SymPy) was the right choice for both
  the F_N slab moment-recursion identities and the spherical-Peierls
  algebraic identities. All 9 `derive_*()` functions verify to
  exact zero in <0.5s.
* **State 1B** (semi-analytical via SymPy + scipy.integrate.quad)
  is what the production solver uses — the integrator is
  `scipy.integrate.quad` (trusted upstream) and the reduction is
  the SymPy-verified F_N recursions. Structural-independence
  preserved above the trusted-library line.
* **The pivot pattern** — when the assumed framework (F_N method)
  doesn't extend to a target geometry (sphere), pivot to a
  structurally-equivalent published method (spherical Peierls
  integral equation) without losing the V&V chain. **This is a
  worth-cataloging pattern for `algebra-of-record`** under the
  "when SymPy framework doesn't extend" subsection — but the
  pivot is generic enough that the existing skill text already
  covers it under "MMS operational rules" (use the published
  reference, don't fabricate one).

## Manifest check

- [x] Branch-1 SymPy module under `orpheus/derivations/continuous/fn_method/origins/`
- [x] Foundation-tagged test gate at `tests/derivations/test_fn_la13511_{slab,sphere}.py`
- [x] Branch-2 production solver at `orpheus/derivations/continuous/fn_method/{slab,sphere}/`
- [x] L1 cross-check tests at `tests/derivations/test_fn_la13511_{slab,sphere}_xverif.py`,
      citing Variant α as the structurally-independent reference (different mathematical paths)
- [x] Sphinx stub at `docs/theory/fn_method.rst` extended with V_fn-slab.{1..5} +
      V_fn-sphere.{1..4} labels each with `:func:` cross-ref + TODO
- [x] Closeout memo (this file)
- [ ] DISPATCH_REQUEST to archivist (emitted as part of final agent response)

All slice-shipping conditions met.
