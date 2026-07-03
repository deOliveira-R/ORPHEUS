---
name: diffusion-integration-290-verification-plan
description: Verification plan for issue #290 (diffusion integration into numerics+transport). The 7-deliverable spec — existing-suite rewire map, new-type intrinsic laws, Mode-12 operator gates, the D=1/(3Σtr) data seam, L2 cross-solver consistency, the #2 DSA seams, tagging. Banked at d2a2a0c.
metadata:
  type: project
---

# Diffusion integration (#290) — verification plan (pre-carve)

**Why:** proactive test-architect dispatch before the operator-algebra carve that replaces
`orpheus/diffusion/` (443-line MATLAB island) with an in-algebra 1-D diffusion operator on
`MaterialMesh` + `ScalarFlux` carriers + a scalar trace leaf + the #182 realizer.
**How to apply:** this is the verification strategy the carve references; the numbered spec is
the deliverable. Re-verify against git before resuming (`d2a2a0c` = main at plan time).

## THE decisive arithmetic finding (crux of the whole rewire)

The legacy `TwoGroupXS.transport = [0.2181, 0.7850]` is hand-entered = Σ_tr directly. `Mixture`
has NO transport XS but carries `SigS[1]` (P1 Legendre moment), so Σ_tr = `SigT − rowsum(SigS[1])`
is derivable (grep-confirmed: no existing `transport_xs` — genuinely new). **Constructing the
diffusion Mixture with `SigS1=0, SigT=transport` reproduces the legacy D BIT-IDENTICALLY**,
`assert_balanced` passes, and the loss operator `A = diag(SigT) − SigS0ᵀ = [[0.0256,0],[−0.0160,0.0959]]`
exactly matches the legacy `_matvec` removal with column-sum `1ᵀA = Σ_a = [0.0096,0.0959]`.
⇒ **RECOMMEND this construction: D bit-identical → the continuous references DON'T move → all L1
anchors survive with mathematical content intact** (the mission requirement). The "physical"
alternative (SigS1=μ̄·SigS0, real SigT) re-baselines D → references move → USER sign-off + widened
tolerances. SigS1=0 is NOT a lie: the diffusion eqn only ever consumes D=1/(3Σtr) and Σ_r; the
transport-corrected value IS the correct isotropic-equivalent encoding.

## USER RULING flags (3)

1. **Vacuum BC semantics.** Current "vacuum" = zero-flux Dirichlet at the physical face
   (`solver.py:225` `dφ/dz=φ_0/(0.5dz)`); theory page `diffusion_1d.rst:81-87` documents this as
   an INTENTIONAL choice ("lets the analytical references be pure sinusoids without an
   extrapolation-length fudge"). #182's realizer introduces Marshak `φ+(2/3)∇φ=0`. RECOMMEND: keep
   zero-flux Dirichlet as the verification anchor; land Marshak as a DISTINCT BoundaryTraceLaw. If
   the campaign redefines "vacuum"→Marshak, the sine references re-derive (extrapolated length) and
   the theory-page note rewrites — USER sign-off.
2. **Data-seam bit-identity vs re-baseline** (see crux above). RECOMMEND bit-identity construction.
3. **Legacy Richardson retirement.** `test_spatial_convergence_reflected` + `derive_2rg` Richardson
   cache are SELF-REFERENCING (H4 — the reference is built FROM the solver). RECOMMEND DELETE
   (superseded by the structurally-independent transcendental ref). The
   `test_2region_transcendental_matches_richardson_cache` cross-check either retires with it or
   re-anchors on a regenerated cache.

## Existing suite (tests/diffusion/, 4 files) — rewire map

`test_continuous_reference.py` (8, L1, verifies 15 diffusion-* labels):
- 4 SOLVER tests (`*_matches_continuous_reference` ×2, `*_flux_shape_converges_second_order` ×2) →
  BEHAVIORAL, REWIRE the `solve_diffusion_1d(CoreGeometry,TwoGroupXS)` call to the new
  MaterialMesh+Mixture+operator API. Math contract (buckled-matrix k, sine φ, O(h²)) INTACT.
- 3 REFERENCE-ONLY tests (`sine_peak_is_at_midslab`, `interface_flux_is_continuous`,
  `vacuum_boundaries_satisfied`) exercise `ref.phi` only (not the solver) → KEEP unchanged
  (continuous refs are NOT retired).
- `test_2region_transcendental_matches_richardson_cache` → cross-check, re-anchor or retire (ruling 3).
`test_diffusion.py` (2, L1): both are legacy-`VerificationCase` duplicates of the continuous-ref
convergence tests → `test_spatial_convergence_bare` DELETE (dup); `test_spatial_convergence_reflected`
DELETE (H4 self-referencing Richardson).
`test_properties.py` (3, L0, REWIRE): `test_vacuum_bc` (weak <0.3peak → sharpen to partial-current
/ Marshak-data on the new scalar trace leaf), `test_flux_positivity` (Perron-Frobenius), `test_flux_symmetry`.
`test_boundary_realizer_stub.py` (3, L0): KEEP the 2 registration tests; the
`test_realize_raises_NotImplementedError` FLIPS to a positive Marshak/Robin realization test when #182 lands.

## New-type intrinsic laws (foundation)

- **Scalar trace leaf** (NEW boundary-locus family on MaterialMesh; `BoundaryField` at
  `_bases.py:665` is SN-only — requires SNMesh + angular TraceSpace). Defining laws (P1 partial
  currents, whether it stores (φ,J) or (J⁺,J⁻)): `J=J⁺−J⁻`; `φ=2(J⁺+J⁻)`; positivity
  `J⁺,J⁻≥0 ⟺ |J|≤φ/2`; Marshak vacuum `J⁻=0 ⟺ φ=2J`; reflective `J=0 ⟺ J⁺=J⁻`; + Field algebra
  (mesh-identity guard, zero).
- **Scalar trace space on MaterialMesh**: dim = n_boundary_faces×ng (slab=2 faces, curvilinear=1);
  membership; zero; `__eq__` by (name,shape).
- **Diffusion operator A_diff**: per-group SPD of `−∇·D_g∇+Σ_r,g` (reflective/vacuum);
  M-matrix sign pattern on the FULL assembled loss (pos diag, non-pos off-diag, diag-dominant →
  non-neg inverse → positive flux); column-sum conservation `1ᵀ(C−K_iso)=Σ_a` per group;
  `(−∇·D∇)@flat=0`; `is_invertible=True`, `.inverse()/.solve` work, `domain==codomain==ScalarFlux`.
  NOTE: full multigroup loss is NOT symmetric (down-scatter lower-triangular) — SPD is per-group only.

## Mode-12 operator gates (the k-references are invariant-functional)

k = λ_max(A⁻¹F). Invariance group {similarity, transpose}. DESIGNED-GREEN mutations: resolvent
factor-order swap A⁻¹F↦FA⁻¹ (similar); full-operator transpose. So k-gates CANNOT catch these →
object-level companions OUTSIDE the stabilizer:
- (a) `A_diff.as_matrix()` on a tiny (3-4 cell, 2G) mesh ≡ independently hand-posed FD stencil
  (element-wise `assert_allclose`/`array_equal`). Mutation-verify: D-face-interp swap, Σ_r vs Σ_t
  confusion, scatter-transpose → reds.
- (b) MMS flux-shape gate (#93). NOT eigenvalue (source-driven). CONFIG-BLINDNESS override: #93's
  proposal (single-group sin, constant D) NULLS both the group-coupling AND the D-gradient `D'φ'`
  term → use HETEROGENEOUS D(x) + MULTIGROUP-coupled manufactured φ (distinct per-group frequency
  activating down-scatter). Needs a NEW theory label `diffusion-mms`.
- (c) ≥2G everywhere (anti-#3; refs already 2G).
- (d) Convergence-order O(h²) paired with the structurally-independent analytic k/sine φ (NOT
  self-referencing Richardson — H4).

## L2 cross-solver consistency (design only)

- Infinite-medium 2G: diffusion (reflective, zero net leakage, uniform slab) vs
  `solve_homogeneous_infinite` k∞. Agrees EXACTLY IF diffusion reuses IsotropicScattering/
  IsotropicN2N/collision (#279 K_iso arm) → SHARED-KERNEL verification (feedback: independence in
  the INPUT assembly, not the shared routine). Verifies the `−∇·D∇` term vanishes at zero-leakage.
  CAVEATS: mixture must have Sig2=0 (homogeneous includes −2Σ₂ᵀ diffusion omits — A/B library have
  Sig2=0); the flat fundamental mode makes `−∇·D∇@flat=0` exact. NOT the existing
  `dif_slab_2eg_1rg.k_inf` (that's the BUCKLED bare-slab k with D·B², not k∞ — build the SAME
  mixture in both solvers).
- Finite-slab diffusion-vs-SN: agree only in the diffusive regime; ~0.3% transport correction.
  L4-INFORMATIONAL, NOT a correctness gate (anti-#1, two ORPHEUS solvers).

## #2 DSA forward seams (name, don't design)

- #215 trap: any future accelerator FP-invariance gate MUST run an ANISOTROPIC config (Mode 9;
  trap doc `scattering.py:1305-1320`). Leave an `xfail(strict=False, reason="DSA #2 not landed")`
  asserting DSA-accelerated converged flux == un-accelerated SI fixed point on VACUUM/HETEROGENEOUS
  (never the isotropic reflective box).
- A_diff contract for DSA: (1) `LinearOperator[ScalarFlux→ScalarFlux]`, (2) invertible, (3) consumes
  the angle-integrated transport residual (moment-0 of `evaluate_residual`'s AngularResidual,
  `sn/solver.py:223`), (4) correction→0 at convergence (FP unchanged). Seam = capability-existence,
  not DSA correctness.

## Tagging

- foundation (no theory label, no verifies()): scalar-trace-leaf laws, trace-space laws, operator
  SPD/M-matrix/conservation, `as_matrix≡FD-stencil`, the transport_xs=SigT−rowsum(SigS1) arithmetic.
- L1 verifies(): continuous-ref eigenvalue+flux gates (15 existing diffusion-* labels); D=1/(3Σtr)
  verifies `diffusion-coefficient` (label exists); MMS needs new `diffusion-mms`.
- L2: infinite-medium diffusion-vs-homogeneous. L4: finite-slab diffusion-vs-SN (informational).
