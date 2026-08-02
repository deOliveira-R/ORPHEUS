---
name: cylindrical-level-ordering-symmetry-adjudication
description: Issue #326 verdict — the curvilinear MMS is EXACTLY blind to the cylindrical per-level ordinate tie-break; the adjudicating instrument is the geometry's xi-mirror symmetry, which says NO ordering is correct (the closure is). Includes the alpha closed form, the half-range exit, and the pole-map leak path.
metadata:
  type: project
---

# Issue #326 — adjudicating the cylindrical per-level ordinate ordering

**Status (2026-08-01):** adjudicated; the three diagnostics are PROMOTED to the permanent
suite and the `derivations/diagnostics/diag_326_*` originals deleted. Branch
`refactor/operator-strategy-layers`. Findings file `scratch/issue326_mms_adjudication.md`.
Gates now at `tests/sn/sweep/curvilinear/test_{alpha_closed_form,azimuthal_mirror_symmetry}.py`
+ `tests/sn/verification/mms/test_mms_ordering_blindness.py` (54 passed + 3 xfail-strict);
the ordering swap is the shared `tests/sn/_test_helpers.product_level_ordering`.

**Why:** `rules_product.py:139` sorts each mu-level by `argsort(mu_x)`; `eta = mu_x` is
2-to-1 over `phi in [0,2pi)`, so the level was never totally ordered and the tie-break was
decided by rounding noise. The user's ruling of record is `lexsort((mu_y, mu_x))`; the open
question was whether that order is *right*.

**How to apply:** if this campaign is picked up, the ordering fix is necessary but NOT the
remediation. The remediation target is the closure (#229). Do not re-run the MMS expecting
it to decide anything.

## The rulings

1. **The curvilinear MMS/L1 suite cannot adjudicate it, exactly.** Both cylindrical MMS
   ansatzes are functions of `eta` and `xi**2` only, so the mirror pair carries IDENTICAL
   `psi_ref`, `Q`, `w`. `alpha` and `tau` are BIT-identical across tie-breaks (alpha is a
   partial sum of `w*eta` with `w` constant in a product level). Measured: lexsort vs
   stable agree to 3e-12 (isotropic) / 9e-15 (anisotropic). Mode 7 + Mode 12, not Mode 10.
2. **No ordering is better.** The xi-ODD companion ansatz `(A + B xi)/W` DOES see the
   tie-break (20.6 % at nx=10) but both orderings hit the SAME angular floor from opposite
   sides (0.28 % at nx=80, spatial order ~0). There is no limit distinguishing them.
3. **The alpha closed form is right, and it selects an unusable ordering.**
   `alpha_k = -w_gl * kappa * [xi(omega_{k-1/2}) - xi(omega_{-1/2})]`, `kappa =
   d_omega/(2 sin(d_omega/2)) = 1 + d_omega^2/24` — EXACT (3e-16, Dirichlet kernel) under
   the AZIMUTHAL ordering. Under the production ascending-eta ordering it is off 2.414x,
   flat in n_phi. The azimuthal ordering makes alpha change sign and `tau_raw` leave
   `[0,1]` -> NaN. Two qualifications on "alpha = -xi": the polar weight `w_gl` scales it,
   and `kappa` is 2.6 % off 1 at the production n_phi=8.
4. **The adjudicating instrument is the xi-mirror symmetry** — a theorem, no reference
   solver. 1-D cylindrical geometry is invariant under `(eta,xi,mu_z) -> (eta,-xi,mu_z)`,
   the product rule is closed under it with equal weights, so the semi-discrete problem
   inherits it EXACTLY. Production breaks it: **1.19e-1** (30 % local) on a homogeneous
   isotropic-source cylinder at `product(4,8)`; **3.08e-1** on `level_symmetric(4)`. It
   falls with `n_phi` and is FLAT in `n_mu` — it IS #229's azimuthal floor, seen without a
   reference. **The defect magnitude is ordering-INVARIANT; only the label moves.**
5. **Mechanism:** two mirror partners share `eta` bit-exactly, so
   `eta_edge[m+1] = (eta[m]+eta[m+1])/2` collapses ONTO the node -> zero-width angular cell
   -> `tau_raw = [0,1,0,1,...]` -> the structural `[1/2,1]` clamp gives `{1, 1/2}`. Two
   geometrically identical ordinates get different angular weights; the tie-break picks
   which. Ascending-eta forces adjacency hence the split; anything else breaks the dome.
6. **Constructive exit — fold to the fundamental domain.** On `omega in [0,pi]` (the
   independent half, since psi is xi-even) `eta` is strictly MONOTONE: no ties, ordering
   unique, `alpha` is simultaneously a non-negative dome AND exactly `2 w_gl kappa xi` at
   the half-angle boundaries (ratio measured 2.0523/2.0129/2.0032 at n_phi 8/16/32), and
   the xi-mirror holds by construction. The whole degeneracy is the redundant half-circle.

## Two secondary findings worth their own follow-up

* **The pole map is the leak path.** `loss_representation/__init__.py:4189` seeds the r=0
  inflow from `pole_outflow[reflection_index("x")[n]]`. A tie-break swap acts inside ONE
  `eta` class; the pole map sends the ordinate to the `-eta` class, where the tie-break
  need not agree — measured NON-COMMUTING for 24 of 64 ordinates. That is why an
  ISOTROPIC-source heterogeneous cylinder still moves 2.6e-2 in psi / 6.6e-4 in phi under
  a tie-break that permutes only identical-data ordinates. Without it the swap would be an
  exact relabeling and phi would be bit-invariant.
* **`reflection_index("x")` is `(eta,xi) -> (-eta,+xi)`; the straight-line axis crossing is
  `omega -> omega+pi`, i.e. `(-eta,-xi)`.** They compose to exactly the xi mirror
  (`rx == ry[rot]`, measured). The shipped pole map is correct up to the xi symmetry the
  solution measurably does not have. FLAGGED, not confirmed as a bug.

## Correction to the issue text

The issue reads the discriminator as spatial HETEROGENEITY. Reproduced its table to all
printed digits (2.607e-02 / 6.387e-03 / 3.75e-14 / 1.522e-15 — it is the ANGULAR flux),
but the homogeneous ~1e-14 legs are H2-degenerate (reflective/reflective + uniform source
=> FLAT exact solution), and a homogeneous problem with a xi-odd source moves 1.14e-1. The
discriminator is the xi-oddness of the per-ordinate data and, for xi-even data, the
pole-map non-commutation.

Related: [[glob-vs-gl-spherical-quadrature-study]] (angular-quadrature study method),
[[curvilinear-inverse-seed-taxonomy]] (the cylinder tau_raw trichotomy).
