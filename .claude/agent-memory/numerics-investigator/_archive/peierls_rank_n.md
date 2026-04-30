---
name: Peierls rank-N investigation
description: Core findings about the rank-N Marshak/DP_N white-BC closure for sphere and cylinder Peierls-Nystrom
type: project
---

# Peierls rank-N white-BC closure — investigation findings

## Rule

For the rank-N closure K_bc = Σ_n α_n · Σ_t · g_n(r_i) · P_esc^{(n)}(r_j) · r_j^{d-1} w_j / R^{d-1}:

- **Cylinder**: α=1 (remove the (2n+1) factor) gives good thin-R convergence (R=1: 21% → 5% → 1.3% at N=3) but drifts up for N ≥ 4. The non-monotonicity beyond the minimum is likely the H1 issue (2D vs 3D cosine projection).
- **Sphere**: G_bc^{(n)} ≡ 4·P_esc^{(n)} (exactly, for ALL n and R) because both share the same angular integrand. This makes the rank-N decomposition u_n⊗v_n a rank-1 outer product in the P_n basis. With α=1, the sphere plateaus at ~15% error at R=1 — the basis does span a rank-6 space but higher modes add little. Needs a different integrand (cosine-weighted) OR different α_n.

**Why**: The existing code at n=0 is bit-exact regression-tested and correct. It treats g_0(r_i)·(1/A) as "response to LOCAL Lambertian inward partial current" and v_0[j] = P_esc^{(0)}(r_j)·r²w as "outgoing partial current per unit volumetric source divided by A". The straight-forward mode-n generalization preserving this structure requires that P_esc^{(n)} and g_n represent the n-th mode of the COSINE-WEIGHTED partial current moments, NOT the plain half-range Legendre moments.

**How to apply**: Before editing, run conservation test (K·1 = 1 identity) — rank-N MUST preserve or improve it. Current code: rank-N makes conservation WORSE by 10× (3% → 7% mean error for sphere at R=10, N=2). That's the smoking gun for a wrong magnitude factor in the mode-n terms.

## Key numerical findings

1. `G_bc^{(n)}(r_i) / P_esc^{(n)}(r_i) = 4` exactly for sphere, all n and all R. Same integrand `∫_{4π} P̃_n(μ_s) exp(-τ) dΩ`; ratio = code prefactors.

2. `G_bc(0) = 4·exp(-Σ_t R)` for sphere (rank-1), analytically verified. Code matches.

3. Conservation test (K·[1] = 1) for homogeneous material with Σ_a=Σ_t:
   - Sphere R=10, N=1: max dev = 5e-3
   - Sphere R=10, N=2: max dev = 68e-3 (WORSE by 14×)
   - This is the clearest diagnostic that current rank-N mode-1 magnitude is wrong.

4. Cylinder α=1 (remove (2n+1)) at R=1: error 21% → 5% → 1.27% at N=3. Marshak-like 10× reduction per rank. Then drift-up for N ≥ 5.

5. Sphere with α=1 plateaus: 27% → 15% → 14.7% (no further improvement). Sphere with any α choice produces plateau. Suggests the mode-n INTEGRAND (not just coefficient) is wrong for sphere.

## Open questions

- Why does α=1 converge well for cylinder but plateau for sphere? Hypothesis: sphere's G_n = 4 P_n identity makes the outer product rank-1 in the P_n basis, which isn't expressive enough. Cylinder's g_n ≠ const·P_n gives asymmetric basis.

- True Marshak m^+_n per unit source via correct reciprocity: brute-force calculation gives `⟨m^+_n⟩ ≠ P_esc^{(n)}/A` — they differ by a 1/μ_s factor that varies spatially. The code's `P_esc^{(n)}` without 1/μ_s is likely NOT the canonical DP_N moment.

- Correct α_n coefficient: my analytical derivation gives α_n = 4π(2n+1) which is 4π times the code's current (2n+1). But code's α_0 = 1 is bit-exact correct. Something about the rank-1 derivation absorbs 4π that I keep missing.

## Scripts

All diagnostic scripts in `derivations/diagnostics/diag_rank_n_*.py`:
- 01: observer vs surface-centered G_bc comparison (confirms sphere reciprocity)
- 02: G_bc^{(n)} = 4·P_esc^{(n)} identity for sphere
- 03: sphere ladders with original, no-(2n+1), and 4·P-substitution variants
- 04: brute-force calibration of ⟨m^+_n⟩ (error: missing 1/μ_s weight)
- 05: factor-π verification (my derivation gives π·φ_bc_code)
- 06: sphere α_n scan (original, (2n+1), 1/(2n+1), π-variants)
- 07: cylinder α_n scan (shows α=1 works for thin R)
- 08: TRUE ⟨m^+_n⟩ with 1/μ_s weighting (differs from P_esc^{(n)}/A)
