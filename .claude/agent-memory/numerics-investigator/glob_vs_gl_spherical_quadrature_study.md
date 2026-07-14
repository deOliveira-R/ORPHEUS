---
name: glob-vs-gl-spherical-quadrature-study
description: Design-study verdict + reusable recipe — Gauss-Lobatto vs Gauss-Legendre angular quadrature accuracy for spherical SN (would GLob let the mu=-1 starting-direction seed become an ordinary weighted ordinate?). GLob tracks GL at a bounded ~1.2x penalty at resolved N; the pole MUST be a straight characteristic (tau_0=0). Standalone scheme-faithful driver method.
metadata:
  type: reference
---

# GL vs Gauss-Lobatto angular quadrature for spherical SN (2026-07-06 study)

Design-phase accuracy study (NO production change). Question: spherical SN
uses GL; its mu=-1 Morel-Montry starting-direction seed psi_{1/2} is a
separate DOF because GL has no node at mu=+-1. Would **Gauss-Lobatto**
(nodes AT +-1, exact only to 2n-3 vs GL's 2n-1) make the seed an ordinary
weighted bulk ordinate, erasing the seed-block type machinery — and at what
accuracy cost?

## VERDICT: GLob is affordable. ~1.2x error penalty at matched, resolved N.
- **At resolved N (>=8), GLob's flux-L2 AND keff error is a bounded,
  roughly-CONSTANT ~1.2-1.4x GL's error** — NOT a divergence, NOT
  anisotropy-amplified. Equivalent to GLob needing ~1-2 extra ordinates.
  Both converge at the same rate; matched-N |GL-GLob| SHRINKS with N.
- **Anisotropy P0->P5 does NOT differentially hurt GLob at resolved N**
  (ratio ~1.2-1.3 for all L). The 2n-3 deficit only touches the highest
  moments, negligibly small once N resolves the flux. The deficit is a
  clean 2-degree loss: GLob first misses even degree 2n-2 (moment probe).
- **Regime/c insensitive**: thin/streaming, thick/diffusive, c=0.5, c=0.99
  all ~1.2x at resolved N. Penalty largest at LOW N (S4: up to ~2x),
  settles by S8-S16.
- **keff** (2g bare sphere): thin R=3 S16 GL_err 8.6e-4 / GLob 1.0e-3
  (|dk|~140 pcm); thick near-critical R=10 S16 GL 2.0e-4 / GLob 2.4e-4
  (~33 pcm). Fine-N GL and GLob keff agree to <6 pcm.
- **Robustness caveat**: GLob NaNs at the pathological N<=L corner
  (S4 + P5: 4 ordinates can't represent 6 Legendre moments — rank-deficient
  for BOTH quads; GL survives numerically but is also wrong). Sane rule
  N>L, N>=6 => GLob stable. Not a GLob-specific defect.

## PHYSICS the study nailed (the 3 caveats)
- **caveat 3 / why no clean production point-swap (caveat 1):** the GLob
  mu=-1 node sits ON the lower angular edge => raw M-M weight **tau_0 = 0**
  (mu_edge[0]=-1=mu_0), so the recurrence psi_{3/2}=(...)/tau_0 DIVIDES BY
  ZERO. mu=-1 IS a straight characteristic ((1-mu^2)=0), so it MUST be
  solved by the radial Carlson DD march (Hebert 3.434/3.435), NOT the M-M
  recurrence, and carries weight w_0 in phi=sum w_n psi_n. (mu=+1 gives
  tau_{N-1}=1, non-singular, natural M-M.) Production's sphere tau is
  UNCLAMPED + its starting_direction_levels predicate keys on tau_raw,0 in
  (0,1) — GLob's tau_raw,0=0 breaks both => wiring GLob into production is
  real surgery, hence the standalone driver.
- **caveat 2 / endpoint weighting**: GLob weights mu=+-1 with w=2/(N(N-1)).
  NO systematic bias — proven by the contamination guard (GLob -> SAME
  N->inf limit as GL, agree ~1e-6 at N=64) and phi=Q/Sig_t exact for GLob.
- **pole treatment choice**: straight-char (caveat-3-correct) vs cylinder-
  style tau-clamp[1/2,1] are COMPARABLE accuracy (within ~1.5x, config-
  dependent); straight-char preferred on PHYSICS, not a decisive win.
  "Closer to GL at finite N" != "more accurate" (all converge to the same
  limit; GL is not truth).

## REUSABLE METHOD (how to compare two angular quadratures in a curvilinear
## SN scheme where they differ in the pole/seed treatment)
1. **Standalone scheme-faithful driver, NOT a production point-swap.** When
   the quadrature sets the angular DISCRETIZATION (M-M tau/alpha + pole
   route), not just moment integration, a point-swap hits singular
   coefficients. Reimplement the EXACT production M-M weighted-diamond
   scheme (verify every coeff file:line: cell_balance.py denom/numer,
   reduced_operator alpha, pole_angular_closure tau/c_in/c_out,
   psi_half_angle_seed Carlson) parametrized by (mu,w,pole_mode).
2. **GATE it bit-faithful to production GL** on NON-FLAT vacuum bare spheres
   (homogeneous-reflective/flat is H2-degenerate). Achieved rel 1e-11 keff,
   1e-10 flux-L2 vs solve_sn. THEN swap only quadrature+pole.
3. **Reference = fine-N GL, with the cross-quadrature CONTAMINATION GUARD**
   (compute fine-N GLob too, confirm GL_inf==GLob_inf) — this makes the
   reference quadrature-family-UNBIASED (vv-6). Report BOTH error-vs-ref AND
   the reference-free matched-N |GL-GLob| difference.
4. **MMS is BLIND here** (L14/vv Mode-7): every curvilinear MMS ansatz is
   <=linear-in-mu = the seed's exact regime, so MMS certifies neither the
   seed nor the quadrature accuracy. Use fine-N convergence + the guard,
   NOT MMS. Anchor correctness with a closed form (k_inf=1.875 hand-derived,
   production nails to 2e-12) + phi=Q/Sig_t streaming equilibrium.
5. **Per-ordinate flat-flux residual** (vv-H3/L6) validates the new pole
   handling PER ORDINATE (angle-integrated phi is degenerate): all pole
   modes gave max per-ord |psi_n - C| ~1e-15.

## Artefacts
Driver: `scratch/experimental/glob_sphere_study/driver.py` (+ run_sweep.py).
Diagnostics (33 tests, all green): `derivations/diagnostics/diag_glob_0{1..5}_*.py`
(01 moment-integration, 02 per-ordinate consistency, 03 end-to-end penalty,
04 pole treatment/tau_0=0, 05 driver faithfulness+k_inf anchor). Promote to
`tests/sn/test_spherical.py` + `tests/numerics/test_rules_1d.py` IFF GLob is
adopted (regenerate references). Backs [[lessons-L16]].
