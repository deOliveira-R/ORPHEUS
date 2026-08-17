---
name: issue-319-flux-dip-discriminator
description: "#319 verdict — the decay rate with optical thickness does NOT split Morel-Montry tau from plain diamond (rate 0.000 for both); the axis that does is h->0. Reconciles #235's 'diamond wins' as wrong-regime + high-angular-order."
metadata:
  type: project
---

# #319 — the flux-dip discriminator, and what it settled for #235

`[M]` 2026-08-12, branch `refactor/operator-strategy-layers` @ `bea6a367`.
Full record: `scratch/q68_flux_dip_discriminator.md`. Raw data (251 solves,
251 `converged`, every row carrying its φ(r) profile):
`/Users/rodrigo/.claude/jobs/c30e4f25/tmp/q68_results.jsonl`. Promotable
gate: ⛔ PROMOTED + RETIRED 2026-08-12 — the diag file is GONE. Successors:
`tests/sn/sweep/curvilinear/test_angular_beta_identity.py` (solve-free beta) and
`tests/sn/verification/analytical/test_angular_diffusion_limit_consistency.py`
(Fig.4 + the h-refinement discriminator); the shared instrument is
`orpheus.derivations.discrete.sn.angular_differencing.morel_montry_beta`
(12 passed, 33 s).

**Why:** #235 measured plain diamond `τ ≡ ½` BEATING the shipped Morel-Montry
τ by 3.0×/2.5× at `n_φ = 16/32` on the anisotropic-cylinder MMS, against
#319's prediction from Bailey–Morel–Chang 2010 that the weighted τ is the
first-order diffusion-limit-consistent choice. #319's experiment was made the
critical path for #235 to decide between them.

**How to apply:** these are settled measurements — do not re-derive them; and
do NOT re-adopt the decay-rate framing (a gate pins its refutation).

## The verdict, in three lines

1. **The DECAY RATE with optical thickness does NOT split the schemes.** At
   constant cells-per-mfp the angular-consistency defect is constant to four
   figures over `Σ_t·R = 5…50` for the shipped τ, for `τ ≡ ½`, AND for a
   garbage control. Rate `0.000` for all three, both geometries. #319's stated
   discriminator is REFUTED.
2. **The axis that DOES split them is `h → 0`.** Shipped τ's defect vanishes
   at exactly first order in `h`; `τ ≡ ½`'s saturates at a finite non-zero
   limit. Ratio `3.2× → 204×` over 2 → 64 cells/mfp. The λ-continuum's
   optimum converges to **exactly λ = 1** (the shipped τ) on two instruments.
   ⟹ the ANGULAR-CONSISTENCY MECHANISM IS CONFIRMED, on a corrected axis.
   This is the primary source's own statement (Morel & Montry 1984 p. 16:
   consistency is recovered "as the spatial mesh is refined").
3. **#235 is reconciled by two measured facts**: (a) the aniso-cyl MMS is
   `Σ_t=1, c=0.5, R=5` ⟹ `Σ_a·R = 2.5` — **not the diffusion limit**;
   (b) even in the diffusion limit the benefit **decays with angular order and
   inverts** — sphere `14× (S2) → 0.9× (S8/S16, shipped marginally worse)`,
   cylinder `5.4× (n_φ=8) → 2.3× (n_φ=16)`. The #235 crossover is reproduced
   on a different fixture, observable and geometry.

## Facts worth not re-measuring

* **Reproduced M&M's published figures**: their Fig. 4 "≈ −1.35" effective
  starting cosine at the origin for Gauss-S2-diamond → measured **−1.3157**
  (10 mfp pure-scattering sphere, 100 zones); Fig. 6's "essentially −1" for
  the weighted diamond → **−0.9884 → −0.9999**.
* **A cheap reference-free τ instrument EXISTS on the sphere and NOT on the
  cylinder.** M&M `β` (Eq. 6a) with the **τ-implied** edge cosines is
  solve-free, exactly 0 for the shipped τ at S2/4/8/16, and the anomaly is
  linear in it. It is **identically 0 for BOTH schemes** on the folded
  cylinder at `n_φ = 6/8/10/16/32` ⟹ refuted there. (This does not contradict
  #235's "BMC β is τ-blind": a β from the STANDARD partition edges is blind
  by construction — different quantity.)
* **The folded cylinder keeps a mesh-CONVERGED ~5–7 % starting-direction
  consistency defect the sphere does not**, and the two consistency conditions
  that coincide on the sphere (`β=0`, `c_s=−1`) decouple there — evidence FOR
  #235's (η,φ) thesis, though #233 (pole closure) or the ψ½ seed are not
  excluded.
* **`λ_opt` on the cylinder converges to ≈0.73 (`n_φ=8`) / ≈0.64 (`n_φ=16`),
  not 1** — but the two instruments disagree at `n_φ=8`, so no optimum is
  claimed. A cylinder τ decision needs a THIRD instrument.
* **The #340 certificate marks the honest floor**: it refused
  `inner_tol=1e-12` at `ε=0.01` (`c=0.9999`, honest residual `2.4e-11`) and
  `1e-13/1e-14` on thinner cases. The deepest diffusion limit is where the SI
  increment stop stops being honest.

## Eliminated, with the structural reason

* *spatial artefact* — `τ ≡ ½`'s anomaly is mesh-CONVERGED (saturates); the
  shipped τ's vanishes at first order, so the shipped τ's residual IS spatial.
* *unconverged inner solve* — metrics bit-stable across 6–7 decades of
  `inner_tol` (`+7.17488e-05` vs `+7.17489e-05`; `c_s` identical to 6 digits).
* *fit artefact* — the fit-free `c_s` gives the same ordering; `D_even` bands
  over 5 windows × 3 degrees never overlap between the schemes.
* *`β` as a cylinder instrument* — identically zero for both schemes.
* *global (all-level) `μ_s` on a cylinder* — reads `+2.76`; the on-axis flux
  is polar-angle dependent, so only the LEVEL-LOCAL form is valid.
* *`c_s` at high N* — presumes ψ affine in the level's angle; at S8/S16 the
  genuine curvature dominates (fixed `−0.640` bias in two different regimes).

Backs [[lessons-L21]]. Sibling: [[curvilinear_tau_clamp_vs_pole_floor]],
[[glob_vs_gl_spherical_quadrature_study]].
