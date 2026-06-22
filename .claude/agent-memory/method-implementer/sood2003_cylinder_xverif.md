---
name: Sood 1999/2003 cylinder Variant α external-reference cross-check closeout
description: V&V hardening Task 1 — cylinder Variant α agreement with LA-13511 Ua-1-O-CY F_N-method critical radius; structurally-independent cross-check landed at 8.5e-6 against Sood's 6-digit reference.
type: project
---

# Sood 1999/2003 Ua-1-O-CY cylinder cross-check closeout

`feature/peierls-greens-cylinder` 2026-05-02. Closes Task 1 of the
post-Phase-3 V&V hardening pass. Cylinder Variant α gains a true
**structurally-independent** L1 external-reference cross-check
parallel to sphere's PS-1982 + Garcia 2021. The cylinder geometry's
verification chain now mirrors the sphere geometry's chain in
structural integrity.

## Deliverables

### Files created

- `tests/derivations/test_peierls_greens_function_cylinder_xverif_sood2003.py`
  — 2 L1-tagged tests (one main `slow`, one fast convergence
  guard).

### Files modified

- None. Test-file-only addition (zero production-code changes); the
  cylinder Variant α solver lands the cross-check on the first
  calibration sweep.

## Test results

### Calibration sweep (the cross-check evidence)

Sood `Ua-1-O-CY`: R = 5.284935 cm, Σ_t = 0.32640, Σ_s = 0.248064,
νΣ_f = 0.176256 (U-235 (a) cross sections, c = 1.30). τ_R = 1.7250
mfp. Driving Variant α at α = 0 (vacuum BC, bare critical) — the
F_N-method reference predicts k_eff = 1.0 exactly at this radius.

| (n_r, n_μ, n_φ, n_t)   | k_eff                  | \|k - 1\| | iter | time   |
| ---------------------- | ---------------------- | --------- | ---- | ------ |
| (8, 8, 16, 24)         | 0.9997284475341046     | 2.72e-04  | 53   | 0.7 s  |
| (12, 10, 24, 32)       | 0.9999666080402017     | 3.34e-05  | 53   | 2.0 s  |
| (16, 12, 32, 48)       | 0.9999521626578750     | 4.78e-05  | 53   | 4.3 s  |
| (20, 16, 48, 64)       | 0.9999861664398738     | 1.38e-05  | 53   | 11.3 s |
| **(24, 20, 64, 96)**   | **0.9999915376857923** | **8.46e-06** | 53   | 24.2 s |
| (28, 24, 80, 128)      | 0.9999951895862697     | 4.81e-06  | 53   | 45.6 s |

**Pinned in the L1 test**: `(24, 20, 64, 96)` — achieves 8.5e-6
agreement at 24 s, comfortably better than the 1e-5 target and
~12× tighter than the 1e-4 hard-floor.

### Convergence properties

- **Iteration count is stable at 53 across all grids.** The
  fission-source iteration sees the same dominant eigenmode at all
  refinements; what changes with grid is the discretization-error
  bias on k_eff, not the iteration dynamics.
- **Convergence is non-monotone at intermediate grids.** The
  (12, 10, 24, 32) result is fortuitously closer to 1 than the
  (16, 12, 32, 48) result (3.3e-5 vs 4.8e-5). This is the same GL
  3D-phase-space behaviour documented in the Phase-1 closeout for
  fuel-A τ_R = 2.5; the (12, 10, 24, 32) grid lands an aliasing
  configuration that happens to favourably integrate the
  eigenmode. Above (20, 16, 48, 64) the trend tightens
  monotonically toward zero.
- **The `test_a2_sood_quadrature_order_convergence` test guards
  against a regression that would change the bias structure** —
  it asserts coarse < 1e-3, fine < 5e-5, and refinement ratio ≥ 4×
  (loose enough to allow the documented intermediate-grid
  non-monotonicity).

### Test gate results

`pytest tests/derivations/test_peierls_greens_function_cylinder_xverif_sood2003.py
-v`:

```
test_a2_variant_alpha_agrees_with_sood2003_cylinder PASSED [50%]   24.2 s
test_a2_sood_quadrature_order_convergence            PASSED [100%] 11.5 s
2 passed in 35.72s
```

Full peierls greens suite: **205 / 205 passed in 619.63 s**
(203 pre-change + 2 new). Pre-existing 203 tests all pass at the
same tolerances they used to.

Sphinx `-W` build: **clean** (exit 0).

## Honest verdict on structural independence

The dispatch correctly framed the question:

> "Both ultimately solve the same transport equation; they differ in
> *method*, not in *physics*."

This is true and the verification framework accommodates it. Two
solvers in **the same project**, sharing **in-house primitives**
(e.g., both calling a project-internal `compute_my_kernel()`), would
fail the structural-independence test by `algebra-of-record` § "above
the trusted-library line". Two **mathematical methods**
(F_N / Wiener-Hopf vs angle-resolved Green's-function-along-bouncing-
characteristics) acting on **the same physics** (1D-radial 1G isotropic
linear transport with vacuum BC), where neither method's reduction
inherits an identity from the other's reduction, satisfy the
structural-independence pillar in its strongest form.

The latent shared assumptions both methods inherit from the
*physics layer*:

1. **The transport equation itself.** Linear Boltzmann, isotropic
   scattering, vacuum BC. This is the *physics they both verify*,
   not a shared assumption that compromises the cross-check.
2. **The cross-section convention.** Both interpret Σ_t, Σ_s, νΣ_f
   as macroscopic constants in cm⁻¹. A units mismatch would show
   up as an order-of-magnitude divergence, not a 1e-5 agreement.
3. **The critical-eigenvalue normalisation.** Both treat k as
   ν · Σ_f · normalisation, with k = 1 ↔ steady fission-rate
   sustainability. A normalisation drift would be detectable as
   a fixed offset across all grids.

None of (1)-(3) are *implementation* assumptions that could fail
identically. They are *physics specifications* that the two
methods are co-verifying. Agreement at 8.5 ppm across two
mathematically-distinct method paths is strong evidence that
ORPHEUS's cylinder Variant α correctly solves the transport
equation it claims to solve.

**The agreement could in principle still mask a bug both methods
share if such a bug existed at a *deeper-than-physics* level** —
e.g., a misinterpretation of "isotropic scattering" or "critical
eigenvalue" that both communities have independently adopted. But
the F_N method's 50+ years of cross-validation against MC,
diffusion-theory limits, and experimental criticals (LA-13511 itself
references benchmark agreement to 4-6 digits) makes this an
extraordinarily-low-prior failure mode. The cross-check pillar this
test occupies is as strong as L1 verification gets.

## Observations and concerns about the cylinder solver

**No concerns surfaced during the cross-check.** The solver
behaviour was clean and predictable:

1. **Convergence is iteration-count stable.** 53 iterations from
   uniform initial guess at every quadrature order — no oscillation,
   no divergence, no NaN. The power-iteration scheme is robust.
2. **Discretization bias is monotone above a threshold grid.** From
   (20, 16, 48, 64) up, the bias tightens monotonically. The
   non-monotonicity at (12, 10, 24, 32) → (16, 12, 32, 48) is a
   known property of the cylinder 3D phase-space that has been
   captured in the Phase-1 closeout memo and the cylinder
   `test_alpha_zero_convergence_floor` test.
3. **No grazing-ray pathology surfaces.** Both grazing loci
   (μ_axial → ±1; sin φ_az → 0) are handled by GL's
   open-quadrature property without special-casing.
4. **The Sood configuration is genuinely lighter than fuel-A.**
   τ_R = 1.725 vs fuel-A's 2.5; convergence is correspondingly
   faster (8.5e-6 at (24, 20, 64, 96) vs fuel-A's ~8e-8 at the
   same grid — but the fuel-A reference value is itself a
   self-consistency claim, not an external truth, so the orders
   of magnitude are not directly comparable).

**Minor numerical noise observation worth recording:**

- The non-monotone bias at (12, 10, 24, 32) → (16, 12, 32, 48) is
  fully consistent with what the Phase-1 closeout memo documented,
  but it is worth noting that **the test_a2_sood_quadrature_order_convergence
  test deliberately allows this non-monotonicity** by comparing
  coarse to fine (skipping the noisy intermediate grids). Future
  research-grade cylinder work that pushes the convergence floor
  below ~1e-7 will need to investigate whether the GL aliasing
  pattern can be flattened by switching to a Gauss-Jacobi (with
  weight 1/sqrt(1-x²) for the φ_az integration on `[0, 2π)` in
  the impact-parameter regime). This is a minor architecture
  follow-up, not a current-day correctness concern.

## Bug catches

**None.** The Branch-2 cylinder Variant α implementation reproduces
Sood's 6-digit critical-radius reference to 8.5 ppm on the first
calibration sweep, with no debugging required.

This clean landing is the **algebra-of-record discipline** working
as designed: the Phase-1 cylinder Variant α already shipped with
V_α1_cyl (closed-cylinder k_inf exactness to 1e-15, the algebraic
ground) and V_α2_cyl (T_00 ≡ P_ss internal cross-check to 1e-10,
production-primitive level). The Sood test adds the
*structurally-independent external* pillar that closes the V&V
hardening pass.

## Position in the cylinder verification chain

After this commit, cylinder Variant α has the following L1 evidence
chain (matching sphere's chain structurally):

| Pillar | Sphere | Cylinder |
| --- | --- | --- |
| Closed-form algebraic identity (V_α1) | k_inf-exactness at α=1, machine precision | k_inf-exactness at α=1, machine precision |
| Internal primitive cross-check (V_α2) | T_00 ≡ P_ss algebraic | T_00 ≡ P_ss numerical (Bickley-Naylor bridge) |
| Semi-analytical external reference | PS-1982 (Pillar 2 — semi-analytical, ≤ 1e-4 across 4 configs) | **Sood F_N method (Pillar — structurally-independent method, 8.5e-6 at one config)** |
| Multi-region external reference | Garcia 2021 (multi-region fixed-source benchmark) | (none yet — open follow-up; cylinder MR is Phase 1b deferred) |

The cylinder Sood cross-check is **stronger** in agreement
magnitude (8.5e-6 vs sphere's 1e-4 against PS-1982) but **narrower**
in coverage (one critical configuration vs PS-1982's four
sub-critical configurations). The narrowness is a property of Sood's
benchmark suite (Table 13 only tabulates U-235 (a) c=1.30 critical
radii for slab/cylinder/sphere; no parametric c-sweep). Achieving a
parametric coverage band on cylinder Variant α at intermediate α
remains an open follow-up logged in the V&V hardening closeout
memo.

## Open follow-ups

1. **Cylinder MR external reference**: no published one-group
   homogeneous cylinder k_eff at α ∈ (0, 1) exists in the surveyed
   literature corpus. Sood α = 0 and V_α1_cyl algebraic α = 1
   bookend the parameter range; interior coverage is via in-codebase
   cross-code (S_N, MOC) IF those solvers support specular cylinder
   BC with tunable α. This is logged in the Sood reference memo as
   "GAP — pillar 4 only".
2. **Multi-region cylinder external reference**: cylinder analogue
   of Garcia 2021 multi-region benchmark not yet identified.
   Future literature-researcher dispatch task.
3. **Gauss-Jacobi quadrature for φ_az**: speculative numerical
   improvement that *might* flatten the (12, 10, 24, 32) →
   (16, 12, 32, 48) non-monotonicity. Out-of-scope for this task;
   logged here as future research-grade improvement for
   architectural follow-up.

## Acceptance gate

- [x] New test passes at the target tolerance (≤ 1e-5; achieved 8.5e-6).
- [x] All 203 pre-existing peierls greens tests still pass at the
      SAME tolerances they used to.
- [x] Sphinx `-W` build clean.
- [ ] Conventional commit `feat(peierls): cylinder Variant α — Sood 2003
      critical-radius cross-check (Ua-1-O-CY)` (pending — see below).

## Pointers

- Test: `tests/derivations/test_peierls_greens_function_cylinder_xverif_sood2003.py`
- Production solver:
  `orpheus/derivations/continuous/peierls_greens_function/greens_function_cylinder.py`
- Sood reference memo:
  `.claude/agent-memory/literature-researcher/sood_2003_cylinder_benchmarks.md`
- Sood PDF (LA-13511):
  `scratch/literature/Sood Foster Parsons (1999)Analytical Benchmark Test Set for Criticality Code Verification.pdf`
- Phase-1 cylinder closeout (template + accuracy-floor reference):
  `.claude/agent-memory/method-implementer/cylinder_variant_alpha_phase1.md`
- V&V hardening closeout (Tasks 2-4, this is Task 1):
  `.claude/agent-memory/method-implementer/vv_hardening_post_phase3.md`
- Sphere PS-1982 cross-check (parallel sphere structure):
  `tests/derivations/test_peierls_greens_function_xverif_ps1982.py`
