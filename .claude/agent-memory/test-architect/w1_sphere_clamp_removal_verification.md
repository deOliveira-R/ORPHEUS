---
name: w1-sphere-clamp-removal-verification
description: W1 (sphere Morel-Montry tau-clamp removal) gate set + the TWO measured refutations of the program plan's premises (iso NOT bit-identical at IEEE-754; aniso O(h2) NOT clean at S16 — the floor is #229-interpolation, untouched by W1). Gates RUN on main 2026-06-13.
metadata:
  type: project
---

W1 = `reduced_operator.spherical_streaming:644`
`tau_mm[n]=max(0.5,min(1.0,tau_raw))` → `tau_mm[n]=tau_raw` (Bailey-Morel-Chang
Eq.42/43, exact-on-linear-in-μ). CYLINDER `:760` KEEPS clamp (structural
`τ_raw=0` ÷0 block). `tau_mm` single-sourced → SI sweep + Krylov matvec both
inherit. Gate set DESIGNED + measured 2026-06-13 (temp-edit + REVERT; production
left clamped for user to land). Files:
`tests/sn/sweep/curvilinear/test_w1_clamp_silent_on_flat.py` (Gates 1a/1b/4) +
W1 section appended to
`tests/sn/verification/mms/test_curvilinear_aniso_convergence.py` (Gate 2).

## ⭐ TWO REFUTATIONS of `.claude/plans/curvilinear_aniso_pole_clamp_program.md`

The program plan (W1 step 1 + W3) asserts "iso-MMS **bit-identical**" and "aniso
L2 → **clean O(h²)**, remove the xfail, assert plain O(h²)". BOTH are
over-optimistic — caught at design-time (the proactive-test-architect payoff):

1. **iso is NOT bit-identical at IEEE-754.** The clamp-silence
   (`c_out−c_in = α_out−α_in`, τ cancels) is EXACT in real arithmetic but the
   closure `(ψ−(1−τ)ψ)/τ` returns ψ exactly only ~81% of the time, ≤1 ULP else
   (reduction-order `α_out/τ − (1−τ)/τ·α_out ≠ α_out` in fp). Converged drift:
   homogeneous-reflective sphere |Δk|=2.3e-13, max|Δφ|=4.4e-13; iso-MMS (sin
   profile, NOT flat) drifts ~1e-5 rel (the per-cell field is not flat-in-μ off
   the fixed point). ⟹ Gate 1 is NOT bit-identity. It is (1a) closure-unit
   few-ULP τ-independence + (1b) converged FP-tail `<1e-9` anchored to k_inf.
2. **aniso O(h²) is NOT clean at the default S16.** W1's aniso floor is the
   **#229 INTERPOLATION floor** (half-angle thread interpolated not imposed),
   which scales with QUADRATURE and is INDEPENDENT of the clamp. Matched-quad:
   S16 nx[10,20,40,80,160]: CLAMP err[5.94e-2,1.51e-2,3.82e-3,1.16e-3,**7.3e-4**]
   vs UNCLAMP[...,1.40e-3,**1.22e-3**] — unclamp CLEANS coarse orders
   (1.995/1.999 vs 1.979/1.978) but RAISES the fine floor (1.2e-3 vs 7.3e-4).
   Clean full-ladder O(h²) needs **S32** (floor<1e-3 → orders 1.985/2.015/2.000
   over nx{10,20,40,80}). ⟹ W3's "remove xfail, assert plain O(h²)" works ONLY
   if quadrature bumps to S32; at S16 the xfail-class limitation PERSISTS post-W1
   (it was never the clamp). The legacy S16 strict-xfails stay until #229.

## The W1 gate set (5 gates, measured)

- **G1a closure-unit** `test_flat_field_coefficient_tau_independent` (@foundation
  @catches ERR-026): recompute clamped+raw τ from SAME `spherical_streaming`;
  for the 8 ordinates clamp bites (S16), assert flat-field cell coeff
  `denom−dA_w·c_in` AND net redist `dA_w·(c_out−c_in)` τ-independent `<1e-12`
  (measured 1.1e-13/1.8e-14). INTRINSIC (holds pre+post W1 — compares the two
  τ paths, not production). The unit proof "clamp silent on flat-in-μ".
- **G1b converged iso** `test_homogeneous_reflective_sphere_iso_unchanged`
  (@foundation @catches ERR-026): homog-REFLECTIVE 2G sphere (flat-flux,
  max/min==1 precondition asserted) vs FROZEN clamped ref k=1.8750000000103633,
  flat φ=0.1326291192, `<1e-9` + anchored to closed-form k_inf=1.875
  (structurally-indep). NOT bit-id. Reflective (not vacuum) → flat → FP-tail
  only. NB user hypothesized this case "DOES shift (homog≠flat under vacuum)" —
  WRONG, it's reflective → flat → does NOT shift.
- **G2 aniso improvement** (3 tests, @slow): `..._S32_clean_o_h2_full_ladder`
  (GREEN both builds, standing claim W1 must not break) + `..._S16_coarse_rate_
  cleaner_unclamped` (the TRUE W1 discriminator: coarse orders clear 1.99,
  clamped reaches only ~1.978 → staged `xfail(strict=True)`, drops same commit
  as W1, xpasses — VERIFIED flips) + `..._floor_scales_with_quadrature`
  (S16/S32 nx=160 floor ratio>2.0, GREEN both, pins floor=#229 not clamp).
- **G3 SI≡Krylov** = EXISTING `test_keff_curvilinear.py::test_si_krylov_
  eigenvalue_equivalence_sphere` (het 2G sphere, non-flat max/min=3.68). Passes
  clamped AND unclamped (|Δk|=5.5e-12 post-W1). NO new test.
- **G4 positivity** `test_unclamped_sphere_flux_strictly_positive` (@l1, 1g R2
  S16 + 2g R2 S64): converged φ strictly>0 + finite post-unclamp (the clamp's
  stated purpose). Measured min 3.98e-2 / 1.33e-1.

## G5 snapshot-shift inventory (RAN live regression on unclamped vs clamped snaps)

`hash()` is per-process SALTED — USELESS across runs; diff actual arrays. 3
sphere cases, two classes:
- **REGEN (1)**: `sphere_2g_3reg_dd_n40` — GENUINE shift (non-flat max/min=1.76),
  k 1.380766→1.381001 (Δ2.35e-4), max|Δφ|=1.7e-3. FAILS the live `assert_
  regression` gate (Δk 2.35e-4 ≫ SAFETY×keff_tol=1e-11). Regen justified.
- **BIT-FP-TAIL, no regen (2)**: `sphere_2g_homogeneous_dd_n20` (flat reflective,
  Δk=2.3e-13, max|Δφ|=4.4e-13) + `sphere_2g_p1_aniso_dd_n20` (flat infinite-
  medium reflective+uniform Q, max|Δφ|=1.1e-12). BOTH PASS the existing gate
  (DriftWarning only) — `assert_allclose(rtol=atol=SAFETY×conv_tol)` absorbs the
  tail (p1_aniso passes via rtol·|38|=3.8e-11≫1.1e-12). Answers user's specific
  Q: sphere_2g_p1_aniso does NOT meaningfully shift (it's flat).
- **BIT-IDENTICAL (10)**: ALL non-sphere (slab/cyl/2D) — confirms W1 touches
  sphere closure ONLY. Regen command (user, on W1 build): `python -m
  tests.sn.regression._generate_snapshots --case sphere_2g_3reg_dd_n40`.

## Reusable lessons

- `hash(arr.tobytes())` is salted per-process — NEVER for cross-run snapshot
  diff. Capture arrays to .npz, `np.array_equal`.
- A "clamp silent on flat" / "iso bit-identical" premise from a closure that
  routes through a τ-bearing reduction is EXACT-arithmetic-only. Always
  distinguish exact-cancellation from IEEE-754: the gate is few-ULP/FP-tail, not
  bit-id. (Generalizes [[regression-tolerance-design]] bit-id-vs-principled.)
- When a fix is claimed to "remove a floor", measure the floor's SCALING with
  the OTHER axis (here quadrature): if it scales the same pre+post, the fix did
  NOT touch that floor — pin the floor as a verified quantity, do NOT assert its
  removal (vv anti-pattern #5/#17). The aniso floor is #229-interpolation, W1
  only cleans the coarse RATE.

See [[curvilinear-aniso-229-9-verification]] (same floor, the #229 retune),
[[lessons]] (L1 homog-blind-to-curvilinear), [[regression-tolerance-design]].
