---
name: Issue 126 Rayleigh-quotient step 1 — F.4 is NOT rank-1 Rayleigh-Ritz
description: Anchor-point + 6-point scan refuting the claim that F.4 k_eff equals R[psi_b=1] under any natural Rayleigh-quotient interpretation on the boundary-trace
type: project
---

# Issue #126 step 1 — verdict: FAIL (step 1 inconclusive; frame needs reformulation)

**Session**: 2026-04-22, numerics-investigator dispatch on step 1 ONLY.

## Headline

F.4 is **NOT** rank-1 Rayleigh-Ritz on the boundary-trace Rayleigh
quotient under any of the four natural interpretations tested. No
interpretation yields the claimed 1e-5 (PASS) or 1e-3 (CLOSE) relative
agreement with F.4 k_eff at the anchor. Best candidate is 11.8% off
(2 orders of magnitude above "CLOSE"). Step 1 fails decisively, so
steps 2 and 3 of Issue #126 should NOT be dispatched until the
variational identification is reformulated.

**Why:** The original claim imports the Rayleigh-Ritz formalism for
eigenvalue problems of the form :math:`A\psi = \lambda\psi` onto
F.4's integral CP structure :math:`\Sigma_t \phi = (K_{\rm vol} + G R P)
(\Sigma_s + \nu\Sigma_f/k) \phi`, which is NOT a self-adjoint eigenvalue
problem on the boundary trace: it is a Schur complement reduction of
the volume problem. No subspace of the boundary trace reproduces
k_eff without explicit coupling through :math:`K_{\rm vol}`. Four
interpretations tested; none works.

**How to apply:** Close Issue #126 step 1 as FAILED. Step 2 (nested
rank-2 Ritz) is not dispatchable without a re-derivation of the
base variational functional. Step 3 (rank-N extension) is even
farther off. Consider retiring the entire issue unless a new
reformulation provides a proper functional for the eigenvalue
problem on the boundary trace.

## The claim being tested

From Issue #126 step 1:

> F.4 IS rank-1 Rayleigh-Ritz on the boundary-trace Rayleigh quotient
> R[psi_b] = <psi_b, GP psi_b> / <psi_b, (I-W) psi_b>
> evaluated at the constant-function trial psi_b(mu) = 1 on [0,1],
> should equal F.4 k_eff to 1e-5 relative at the anchor (tau=10, rho=0.3).

## What was tested

Diagnostic script: `derivations/diagnostics/diag_rayleigh_quotient_f4.py`.

F.4 reference k_eff is computed via `run_scalar_f4` (production
scalar F.4 closure with Marshak W + Lambert P/G) at BASE quadrature
(`n_panels=2, p_order=4, n_ang=32, dps=15`), which gives ~1e-4
accuracy vs RICH (adequate for testing 1e-5 to 1e-3 identification
tolerances on a k_eff of order 1.5).

P, G, W primitives are computed directly via `mpmath.quad` at tol ~1e-10
(independent of the production Gauss-Legendre radial grid):

- `P_esc_outer_mp(r)`, `P_esc_inner_mp(r)`: Lambert mode-0 escape at r
  (integral over c with K_esc = exp(-sigma rho_out), sphere direct form).
- `W_ss_mp`: 2x2 Marshak surface-surface transmission per production
  `compute_hollow_sph_transmission`.
- `G_face_k(r) = 4 * sigma_t * P_esc_face_k(r) / div_k`  (from the
  production identity `G_bc = 2 * int sin(theta) K_esc dtheta = 4 * P_esc`
  on sphere, then `G[i,k] = sigma_t * G_bc / div_k` with
  `div_outer = R^2`, `div_inner = r_0^2`).

Four interpretations of R[1] attempted:

- (I) `R_direct = <1, GP 1>_{face-space} / <1, (I-W) 1>_{face-space}`
  where GP is the RANK-1 outer product of face-integrated P and G
  (2-face bilinear). `psi_b = (1, 1)` in face space.
- (II) `R_direct` reinterpreted as K_total_1, with
  `k_eff = R * nsf / (sigma_t - R * sigma_s)`.
- (III) `K_bc_closure = <1, G (I-W)^{-1} P 1>` (with resolvent).
- (IV) `R_face = <1, PG 1> / <1, (I-W) 1>` where PG is the 2x2 face-trace
  matrix `(PG)[k, l] = sigma_t / (pi * div_l) * int_V P_k P_l dV`.

Plus a diagnostic K_total baseline: K_vol_1 = `(1 - P_esc_total_bar) / sigma_t`
(collision-probability identity on uniform flux), combined with
K_bc_closure via `K_total_1 = K_vol_1 + K_bc_closure` → k_eff.

## Anchor result (tau=10, rho=0.3)

| Interpretation           | Value        | k_eff from value | rel err vs F.4 |
| ------------------------ | ------------ | ---------------- | -------------- |
| F.4 reference            | -            | 1.4963034787     | -              |
| (I)  R_direct            | 2.58e-04     | 2.58e-04         | 9.998e-01      |
| (II) R_direct as K_total | 2.58e-04     | 2.58e-04         | 9.998e-01      |
| (III) K_bc closure only  | 2.58e-04     | 2.58e-04         | 9.998e-01      |
| K_vol_1 + K_bc (full)    | 0.9166       | 1.3199           | 1.179e-01      |
| (IV) R_face_trace        | 1.918e-01    | 0.2049           | 8.718e-01      |

Best candidate: **K_total_uniform (full K_vol + K_bc on uniform
flux)** at **11.79% relative error**. This is the CORRECT mechanical
computation of k_eff on the span{1_V} trial function in the production
convention, and it's still 100× above the "CLOSE" tolerance (1e-3).

The full k_eff problem on span{1_V} inherits from the Rayleigh-Ritz
framework (any variational problem admits such a projection), but
the eigenvalue converges to k_inf = 1.5 only when the true eigenvector
approaches the constant. On a leaky sphere, the true flux is
non-uniform and the constant-function trial is far from optimal.

## 6-point grid scan

| tau | rho  | k_f4     | K_tot_unif rel err | R_face rel err | k_from_R_face rel err |
|-----|------|----------|--------------------|-----------------|------------------------|
| 5   | 0.30 | 1.49817  | 2.217e-01          | 8.598e-01       | 8.492e-01              |
| 10  | 0.30 | 1.49630  | 1.179e-01          | 8.718e-01       | 8.631e-01              |
| 20  | 0.30 | 1.49453  | 5.819e-02          | 8.695e-01       | 8.604e-01              |
| 5   | 0.50 | 1.49767  | 2.724e-01          | 8.361e-01       | 8.215e-01              |
| 10  | 0.50 | 1.49598  | 1.490e-01          | 8.644e-01       | 8.546e-01              |
| 20  | 0.50 | 1.49580  | 7.556e-02          | 8.659e-01       | 8.563e-01              |

The K_total_uniform error SHRINKS with tau (as expected — thicker
sphere has less leakage, constant-flux trial gets better), but STILL
exceeds 5% at tau=20. The R_face_trace error is ~85% and near-constant
— no hint that any uniform-rescaling trick would fix it.

## Why the claim fails

The spec's Rayleigh-quotient formula
`R[psi_b] = <psi_b, GP psi_b> / <psi_b, (I-W) psi_b>`
DROPS the volume contribution K_vol, which represents the majority of
the k_eff signal at finite sphere. The BC contribution `G (I-W)^{-1} P`
is only ~1e-4 of K_total at tau=10, rho=0.3. No amount of "reinterpreting"
R[1] recovers K_vol's 0.91 contribution to K_total.

The "alpha ~ 0.38" scalar gauge documented in
`diag_lambert_marshak_symbolic.py` Section 5 / E7 is a DIFFERENT
quantity: it is the ratio `S_M / S_L` of two volume-integrated trace
sums (Marshak vs Lambert), NOT a Rayleigh quotient. The side-observation
that sum(PG) ≈ 0.38-0.40 at (tau=10, rho=0.3) is consistent with the
alpha documentation but does not prove a variational identity.

## Proposed Sphinx edit

**None warranted.** The identification is NOT supported by the numerics.
No edit to `docs/theory/peierls_nystrom.rst §peierls-f4-rank-1-gauge-why`
should be made — the current text (classifying F.4 as a lucky rank-1
algebraic accident via change-of-basis scalar M^(1) = sqrt(2)/2) is
more accurate than any Rayleigh-Ritz framing.

If a future session develops a properly reformulated variational
identification (e.g., on the FULL eigenproblem including K_vol, not
just the BC closure), a new Sphinx subsection could be added. But
that is out of scope for step 1 and should not be written pre-emptively.

## Recommendation on Issue #126

**Close step 1 as failed. Do not dispatch step 2 or 3.**

The rank-2 nested Ritz subspace of step 2 presupposes that the base
rank-1 variational principle holds. With step 1's falsification, the
full stack is unfounded. If Rayleigh-Ritz is to be salvaged, the
reformulation must:

1. Identify a **self-adjoint operator** on the boundary trace whose
   dominant eigenvalue is (a monotone function of) k_eff. F.4 does
   not present such an operator; the BC-inclusion is a Schur complement,
   not an eigenvalue reduction.
2. Or: identify a **volume-trace functional** that IS Rayleigh-Ritz
   on the full span{phi}, and examine whether the rank-1 F.4 closure
   is the rank-1 Galerkin projection of it. This is a substantial
   theoretical lift, well beyond "1 session for step 1".

The frame-attack premise (rank-N non-monotonicity = variational
abandonment) may still have merit, but the rank-1 variational anchor
is not F.4 under any tested interpretation.

## Diagnostic script

`derivations/diagnostics/diag_rayleigh_quotient_f4.py` (newly created).
Self-contained; runs in ~4 minutes on BASE quadrature, evaluates all
four interpretations at all 6 reference points, reports the full
table. Does NOT need promotion to tests — it is a negative-result
diagnostic; it confirms the claim fails, not that any feature works.
