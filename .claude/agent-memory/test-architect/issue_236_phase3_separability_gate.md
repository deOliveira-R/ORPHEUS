---
name: issue-236-phase3-separability-gate
description: #236 Phase 3 (ST5) — permanent (spatial ⊗ angular) discretization-error SEPARABILITY characterization gate. Cartesian SEPARATES additively (N-independent O(h²) rate + |M|/max≤0.005); curvilinear GATES (coarse-N spatial saturation + fine-N O(h²) recovery + |M|/max≥0.4 sphere). DRAFT GREEN under -O 2026-06-18, zero solver code.
metadata:
  type: project
---

LAST piece of the #236 (spatial ⊗ angular product) campaign (Phase 1 pairing +
Phase 2 τ-carve already shipped @ `cdb8f07`, branch `feature/sn-spatial-angular-product`).
File: `tests/sn/verification/mms/test_space_angle_separability.py` (6 tests, `@l1`,
`@slow`, DRAFT GREEN under `.venv/bin/python -O`, ~3.0s, zero production code).
Modelled on the W2 pole-cell CHARACTERIZATION-gate pattern
([[w2-pole-cell-characterization-gate]]) and the mirror
`test_curvilinear_aniso_convergence.py` idioms (`_l2_1d`, the curvilinear MMS
builders, `solve_sn_fixed_source(...,max_inner=500,inner_tol=1e-13)`).

## THE CLAIM (geometry-split, lit + empirical backed)
- **Cartesian (slab-iso / slab-P1):** (1−µ²)/r ABSENT → space ⊥ angle → error
  SEPARATES additively `E(h,N)≈E_space(h)+E_angle(N)`. Signature: spatial O(h²)
  RATE is N-INDEPENDENT + mixed 2nd diff `M=E[h1,N1]−E[h1,N2]−E[h2,N1]+E[h2,N2]`
  vanishes (`|M|/max ≤ 0.005`).
- **Curvilinear (sphere/cyl):** M-M angular thread feeds the shared cell-balance
  denom → error GATES `E(h,N)≈max(E_space,E_angle)`. Signature: coarse-N spatial
  refinement SATURATES (h-ratio→1) while fine-N RECOVERS O(h²); LARGE `|M|/max`.

## THE HARD LEG (decision 3): N-GATED SPATIAL RATE, not "cross-term is large"
The robust POSITIVE gating signature is option-(i) **the spatial rate is
quadrature-gated**, NOT the fragile negative "cross-term stays large". Sphere
`test_sphere_spatial_rate_is_quadrature_gated`: coarse N=8 finest h-ratio
SATURATES `<1.6` (measured 1.15), fine N=32 RECOVERS O(h²) `>3.0` (measured 4.0),
AND `r_n32[-1] > 2×r_n8[-1]` (the rate genuinely N-depends). The cross-term leg
(`>0.15`, measured 0.41 vs Cartesian ≤0.005) CORROBORATES — a DISCRIMINATION
threshold sitting in the 3-orders gap between Cartesian and gated, not a brittle
0.41 pin. Cylinder = the gating extreme (NO O(h²) window at any quad): coarse-n_phi
spatial saturation (`<1.5`, measured 1.02/1.00) + azimuthal floor DROPS `>2×`
(measured 2.58×) when n_phi 8→16.

## ⭐ L27 (decision 4): per-ordinate leg ADDED (not blind by construction)
The scalar weight-sum metric is blind to a wrong angular closure (α-dome
telescopes under Σ w_n ψ_n — vv L27/H3). Because the curvilinear GATING is itself
an angular-closure phenomenon, `test_curvilinear_gating_per_ordinate_not_blind`
reproduces the sphere gating from the max-over-ordinates per-ordinate L2.
⭐⭐ CONVENTION TRAP (measured): `case.psi_exact(r,mu_n)` returns `A+B·mu_n`
WITHOUT the 1/W factor (its documented contract — "for caller convenience"); the
solver stores per-ordinate WITH 1/W (`Σ w_n ψ_n^solver = φ` exactly). MUST divide
`psi_exact/sum_w` or a 2× (=W) mismatch reads per-ord L2 ~8.1 with FLAT h-ratios
(a meaningless metric). The leg carries a sanity guard `po_n8[0]<1.0` that reddens
to `8.111e+00` if the 1/W is dropped (mutation-verified live 2026-06-18).

## DECISIONS
- **Level `@l1`** (MMS-convergence-STRUCTURE = a math claim, same as the W2
  pole-cell gate). NOT eigenvalue (MMS can't reach that layer).
- **`@slow` on all** (curvilinear solves dominate; characterization gate not a
  fast canary). Reduced grid `nx∈{20,40,80}` (resolves both gated-collapse and
  O(h²) cleanly) — full grid `nx≤160` is ~13min, reduced is ~3.0s total.
- **`@catches("ERR-026")`** ONLY on the 2 tests with the MUTATION-PROVEN direct
  reddening: the sphere `r_n32>3.0` + per-ord `r_n32[-1]>2.0` O(h²)-recovery
  asserts (a re-emerged ERR-026 = divergent wrong FP, error GROWS with refinement
  → ratios <1 → assert RED; proxy-verified 2026-06-18). REMOVED from the
  cross-term + cylinder legs (M-M active there but no proven direct reddening →
  would be a phantom coverage claim per vv `catches`-is-a-coverage-claim).
- **`verifies(<label>)` = NONE YET.** No separability/factorization `:label:`
  exists in `discrete_ordinates.rst` (the 3 `tensor-product` labels are the
  LD/UBLD BASIS, not the space×angle error separability). ⛔ OPEN: archivist must
  MINT one. Proposed: `sn-space-angle-separability` (the LMM-1987-spatial ×
  BMC-2010-angular split). Once minted: Cartesian legs carry it (the separability
  CLAIM); curvilinear legs carry it as the GATING signature.

## MODE-8 verified LIVE (2026-06-18)
pytest rewrites bare `assert` in collected test_* modules → they FIRE under `-O`
(the PytestConfigWarning covers only NON-test-module asserts). Proof: flipped
`r_n8[-1]<1.6`→`<0.5` → RED under -O showing real 1.1529. All asserts kept in
collected test fns; informative messages mirror-test style.

## MEASURED (branch HEAD cdb8f07, .venv/bin/python, reduced grid nx[20,40,80])
- SPHERE scalar: N8 [1.4665e-2,5.4016e-3,4.6851e-3] ratios[2.71,1.15]; N32
  [1.5006e-2,3.7136e-3,9.2864e-4] ratios[4.04,4.00]; |M|/max=0.411.
  per-ord(/W): N8 ratios[1.60,1.16]; N32 ratios[3.39,2.97].
- CYLINDER scalar (n_mu=4): n_phi8 [1.9504e-2,...] ratios[1.02,1.00]; floor drop
  n_phi8→16@nx80 = 2.58×; |M|/max=0.019 (small only bc E≈E_angle swamps).
- SLAB-ISO: N4 & N16 both ratios[4.01,4.00] (N-indep O(h²)); |M|/max=0.0047.
- SLAB-P1: floor flat ~6.8e-3 N-indep (<0.3%); |M|/max=0.0038.
Four `diag_sep_*.py` probes reproduced the documented cross-terms BIT-FOR-BIT post
Phase-2 τ-carve (confirms bit-identity held). Probe job tmp `84fd66f8` (the
space_angle probe's cylinder leg uses n_ordinates= which the cyl builder rejects —
use n_phi=; fixed copy `diag_sep_space_angle_cylfix.py`).

Extends [[w2-pole-cell-characterization-gate]] (the characterization-gate design)
+ the #236 campaign. Relates #229 (azimuthal floor — if lifted, the gating legs
redden and must be RE-TUNED to the better regime, NOT deleted — that re-tune is
the intended signal), #233 (pole-cell O(h)), ERR-026 (the M-M closure the
curvilinear legs activate).
