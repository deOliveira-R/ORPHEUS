---
name: issue-196-phase-g-step2-replan-blocker
description: Phase G Step 2 replan blocker — separate-pass M-M Picard diverges; legacy SI triangular forward-sub converges to a wrong fixed point at machine precision. Both routes have been empirically characterised. Decision-point checkpoint TRIPPED per the brief; dispatching numerics-investigator before any code lands.
metadata:
  type: project
  branch: refactor/sn-operator-algebra
  step: Phase G Step 2 replan (BLOCKED)
  date: 2026-05-13
---

# Phase G Step 2 replan — BLOCKER memo

**Branch**: `refactor/sn-operator-algebra` tip `59f3604`. Working tree clean (no Step 2 production code committed).

**Plan**: `.claude/plans/issue_196_phase_g_replan.md` § "Step 2 — Canonical curvilinear SI sweep math (no GMRES)".

**Brief decision-point checkpoint TRIPPED**: "Picard appears to diverge on the candidate canonical operator" → STOP and dispatch numerics-investigator. The previous Step 2 agent shortcut this by routing through GMRES; this attempt captured empirical Picard divergence evidence + a deeper failure mode and is dispatching properly.

## What I implemented (separate-pass M-M with prev-iter ψ)

Per the diagnostic memo § "Implementation paths" recommendation (a) — recast as cell-balance + angular-redistribution composition — I:

1. Added `cell_balance_terms_streaming_only` + `cell_balance_terms_streaming_only_degenerate` helpers in `orpheus/sn/spatial/cell_balance.py`. They strip the M-M contributions (`dA_w · c_out` from `denom`, `dA_w · c_in · ψ_angle_in` from `numer_upstream`) so the per-cell solve handles only streaming + collision.

2. Rewrote `_sweep_1d_spherical` in `orpheus/sn/sweep.py` to:
   - Maintain a persistent `psi_cells_iter` buffer in `psi_bc["psi_cells_sph"]` across outer iterations.
   - Per sweep call: compute Carlson seed from source `Q_bar = 0.5 · Q_1d` + face-flux `bc_outer_value`.
   - Compute `R = sn_mesh.pole_angular_closure(psi_cells_iter, ...)` ONCE as a separate pass over the input ψ field.
   - Per ordinate × per cell: solve streaming-only balance with `effective_source = QV - R[g,n,i]·V` and `numer_upstream = |μ|·(A_in+A_out)·ψ_spat_in` (no `ψ_angle_in` term).
   - WDD spatial closure unchanged: `ψ_face_out = 2·ψ_avg - ψ_face_in`.
   - At end of sweep: persist `psi_cells_iter` for next outer iteration.

Diff saved to `/tmp/step2_separate_pass_attempt.diff` (570 lines). NOT committed.

## Empirical evidence — Picard diverges

Test problem: homogeneous reflective sphere (mixture B: Σ_t=2, Σ_s=1.9, no fission), 40 cells, GL-8, isotropic Q_within = 20 (matches converged scattering source at φ=10), within-group `(L+C)·ψ = Q_within` direct test.

Picard scheme: ψ^{k+1} = sweep(Q_within, σ_t, sn_mesh, psi_bc) — where each sweep computes R from `psi_bc["psi_cells_sph"]` (the previous iterate's ψ field).

Residual history (`||ψ^k - ψ^{k-1}||_2 / ||ψ^k||_2`):

```
iter 0:  ||psi||=1.55e+01, sf[0]=21.4, sf[-1]=10.7
iter 1:  ||psi||=2.59e+01, rel_diff=1.4e+00 — diverging
iter 5:  ||psi||=5.30e+05
iter 10: ||psi||=1.10e+13
iter 20: ||psi||=2.01e+23
iter 30: ||psi||=2.89e+33
iter 40: ||psi||=4.17e+43
iter 49: ||psi||=5.79e+52
```

**Empirical spectral radius ρ((streaming+collision)^{-1} · R) ≈ 10.4** (consistent growth factor across iterations 5-49).

Sign-pattern: `sf[0]` oscillates with sign flips (`+`, `-`, `+`, `-`, ...) — divergent oscillation, not monotone amplification.

**Conclusion**: the separate-pass M-M Picard scheme with R built from prev-iter ψ has unbounded spectral radius. The brief's claim "Krylov GMRES converges on the same operator at machine precision — Picard must too" is not quite right: GMRES is solving the LINEAR operator equation `(L+C)·ψ = q` via Krylov (which always converges for nonsingular operators), but the Picard scheme implements a DIFFERENT operator equation — `ψ^{k+1} = (streaming+collision)^{-1}·(q - R(ψ^k))` — whose fixed-point convergence requires `||·(streaming+collision)^{-1} · R||_∞ < 1`, which it's not.

## Empirical evidence — legacy SI converges to wrong fixed point at machine precision

For comparison, ran the legacy SI sweep (current Phase F state — per-cell M-M coupling via `ψ_angle_in/out` Gauss-Seidel within sweep) on the SAME L0 problem with `max_inner=2000`, `inner_tol=1e-14`:

```
SI iterations: 681, residual: 9.5448e-15 (fully converged)
sf_si: i=0=7.804371, i=20=10.032615, i=39=10.030782
psi at pole (n=8 GL): [4.5736 4.7208 4.7000 4.6168 4.2190 3.2422 2.0733 1.3278]

Krylov iterations: 571, residual: 9.9879e-15
sf_kr: i=0=10.000000, i=20=10.000000, i=39=10.000000
psi at pole (n=8 GL): [5.0 5.0 5.0 5.0 5.0 5.0 5.0 5.0]
```

The legacy SI **fully converges** (residual 9.5e-15) but to a **wrong fixed point** (22% L0 error at pole, range across ordinates 1.33-4.72, badly anisotropic).

The legacy SI sweep IS a one-shot exact inversion of its operator (triangular forward substitution in ordinate per cell — the M-M is causal). The Picard outer iteration (`_solve_fixed_source_si`) converges to a fixed point of the operator the SI sweep inverts. **That operator is NOT the same as the operator the Krylov path inverts.**

## My analytic walkthrough — operators DO differ at fixed point

Per-cell apply matvec at cell i ordinate m (substituting M-M recurrence + WDD):

```
denom · ψ_cell[m,i] = source[m,i] + |μ|·(A_in+A_out)·ψ_face_in[m,i] + dA_w·c_in·ψ_half_left[m,i]
where denom = 2|μ|·A_out + dA_w·c_out + Σ_t·V
```

**Apply matvec inputs**:
- `ψ_face_in[m,i]`: WDD-propagated from cell i-1 using INPUT ψ_cell at upstream cells.
- `ψ_half_left[m,i]`: M-M-propagated from ordinate m-1 using INPUT ψ_cell[m-1,i].
- `ψ_half_left[m=0,i]`: Carlson seed from `CarlsonInwardSweep(psi_cells, context)` with `Q_bar = 0.5·σ_t·φ_0(psi_cells)` and `bc_outer_value` from `bc_outer.apply(fi[:, :, -1, 0])` (cell-centre proxy at outer cell).

**Legacy SI inputs (same per-cell equation, different sources)**:
- `ψ_face_in[m,i]`: WDD-propagated from cell i-1 using **just-solved** ψ_cell[m,i-1] (during this sweep call).
- `ψ_half_left[m,i]`: M-M-propagated from ordinate m-1 using **just-solved** ψ_cell[m-1,i] (during this sweep call).
- `ψ_half_left[m=0,i]`: Carlson seed from `carlson_inward_sweep_from_source(Q_bar=0.5·Q_1d, σ_t, dr, bc_outer_value=bc.apply(prev-sweep-bc_outer)[most_inward])` (Phase F backport, face-flux BC trace).

**At convergence ψ*** = ψ*_si AND ψ*_si = ψ*_kr** would require the SI's just-solved values to equal the Krylov's input values at the same ψ*. They DO — both equal ψ* by definition.

So if the per-cell algebra is identical AND the inputs at fixed point are identical, the operator solutions MUST match. The empirical disagreement (SI converges to a wrong ψ* while Krylov converges to ψ_analytical) implies one of two things:

**Hypothesis H_A**: there's an algebraic discrepancy between SI's per-cell solve and the apply matvec's per-cell action that I haven't identified. Most-likely-culprit list:
- SI uses `ψ_angle_in[i]` from `result.outgoing_angular_state` of the PREVIOUS visit at this cell — which equals `(ψ_avg - (1-τ)·ψ_angle_in_prev)/τ`. Apply uses `ψ_half_right` from M-M recurrence on `ψ_level`. These should be the same value AT THE FIXED POINT but the construction sequence differs at the dataclass level (`outgoing_angular_state` is recomputed from each cell's solved ψ_avg; M-M recurrence is computed once over the whole field from input ψ).
- BC trace face-flux (SI) vs cell-centre (apply) differs at finite resolution. For homogeneous flat ψ they happen to give the same value at the converged solution, BUT the OPERATOR being inverted differs at finite resolution.

**Hypothesis H_B**: the converged ψ*_si actually IS a valid fixed point of the apply matvec's operator equation. But then Krylov's ψ*_kr is also a fixed point. The system would have multiple solutions — implying the operator is singular or near-singular. Counter-evidence: GMRES converges at machine precision, so the operator is non-singular.

**Hypothesis H_C**: the SI's converged ψ*_si is a fixed point of a SIMILAR-LOOKING operator that arises from the in-sweep Gauss-Seidel propagation of ψ_angle. Specifically: the SI does NOT actually solve the same per-cell linear equation when the inputs differ between apply (input ψ everywhere) vs SI (just-solved ψ here, input ψ elsewhere). At fixed point these inputs converge to the same ψ*, but only IF the in-sweep dynamics consistently reproduce the equilibrium — which requires the operator to be "Gauss-Seidel-equivalent" with itself. Non-symmetric operators can have G-S fixed points that diverge from Jacobi (Krylov-applied) fixed points if the in-sweep ordering introduces a structural bias.

H_C is the most likely. The fixed-point operator equation is the same algebra, but the SI implicitly computes its fixed point via G-S relaxation while Krylov computes its via Jacobi-style operator application. These two relaxation schemes have **different fixed points** if the underlying linear system has off-diagonal coupling whose ordering matters.

Note this is the standard L1-vs-L0 distinction in iterative linear algebra: triangular forward substitution (exact one-shot) is equivalent to Jacobi only if the matrix is ALREADY upper triangular. For the curvilinear operator, the M-M coupling in ordinate is triangular, but the WDD coupling in space is also triangular per direction. The SI's combination of these two triangular relaxations does NOT in general give the same fixed point as inverting the combined operator (because the cell visits at each ordinate use the just-solved spatial-upstream ψ_face_in at that ordinate, which depends on this-ordinate's ψ_avg AT a different cell — that's a transitive dependency back through ordinates).

## Specific symptom to investigate

The L0 SI converges to ψ ranging 1.33 (μ=+0.96) to 4.72 (μ=-0.96), badly anisotropic at the pole. The legacy SI's per-cell algebra applies the same equation across all ordinates symmetrically — so the asymmetry can ONLY come from the M-M Gauss-Seidel ordering. Specifically: the SI processes ordinates m=0,1,...,N-1 sequentially with `ψ_half_left` propagated forward. For high-|μ| outward ordinates (m=N-1, μ=+0.96), this means `ψ_half_left` has accumulated all the recurrence drift from m=0..N-2. For high-|μ| inward ordinates (m=0, μ=-0.96), `ψ_half_left = carlson_seed` (fresh start, no drift).

This is the **structural sweep-order asymmetry**. Krylov's apply matvec doesn't have this asymmetry because it uses INPUT ψ everywhere — every ordinate is built from the same input.

## Three implementation paths I see

Each carries a tradeoff:

**Path 1 — Block-Jacobi per cell**: per cell, build the N×N angular block matrix (per group) including M-M coupling AND the streaming contribution to ψ_cell. Invert this dense per-cell matrix to get all N ordinates simultaneously. **Cost**: O(N³·nx·ng) per sweep vs O(N·nx·ng) for triangular. **Outcome**: should match Krylov fixed point (proper one-shot inversion of the per-cell block).

**Path 2 — Iterate within sweep**: lift the SI sweep to a 2-level Picard. Outer level: Gauss-Seidel through ordinates. Inner level (per ordinate, per cell): re-solve to bring the just-solved-ψ-based M-M term up to date. **Cost**: O(inner_iter·N·nx·ng) per sweep. **Outcome**: should converge to the per-cell block-Jacobi fixed point.

**Path 3 — Restructure to match Krylov's operator exactly**: rewrite the SI to use INPUT ψ for the M-M coupling — separate-pass M-M as I attempted, BUT with a contraction-stable Picard scheme. Options:
- Use damping: `ψ^{k+1} = (1-ω)·ψ^k + ω·streaming_only^{-1}·(q - R(ψ^k))` with `ω < 1/spectral_radius`. Requires estimating spectral radius first.
- Convert to a synthetic-acceleration form (DSA — Diffusion Synthetic Acceleration). Standard in production SN codes for this exact symptom.

**Path 4 — Just don't ship a one-shot SI sweep for curvilinear**: declare the SI inner solver path for curvilinear unsupported, with the dispatch defaulting to Krylov. This matches the existing solver.py defaults (Phase D set Krylov default for curvilinear). Phase F's SI Carlson seed backport becomes a "history-only" intervention that improved-but-did-not-close ERR-026 manifestation #6.

**Path 4 is brief-violating**: the brief explicitly requires SI to pass the L0 test at `rtol=1e-9` on `inner_solver="source_iteration"`.

**Path 1 is the most principled but most expensive**. The per-cell N×N inversion is small per cell (N ≤ 8 for GL-8, ≤16 for higher orders), so the cost is bounded; total: O(N²·nx·ng·N_sweeps) — same asymptotic order as Krylov-with-N-iterations. Not unreasonable.

**Path 3 with DSA is the production-standard approach** but introduces a new operator (the diffusion synthetic acceleration step) — significant scope expansion.

## Hard recommendation

**Dispatch numerics-investigator with this evidence to:**

1. Confirm or refute H_C (G-S vs Jacobi fixed-point asymmetry).
2. Confirm spectral radius diagnosis (separate-pass Picard's ρ ≈ 10.4 unbounded).
3. Audit my separate-pass implementation diff (`/tmp/step2_separate_pass_attempt.diff`) for algorithmic correctness — is the R coupling correctly folded in?
4. Recommend the correct path forward (Path 1, 2, 3-with-DSA, or some other architecture I haven't seen).

**DO NOT** route through GMRES inside the SI sweep — that was the previous agent's Cardinal Rule 1 violation. The brief explicitly forbids it.

## State preserved

- Working tree restored to `59f3604` (Phase F state) — `cell_balance.py` + `sweep.py` reverted.
- My separate-pass attempt is saved as `/tmp/step2_separate_pass_attempt.diff` (570 lines) for the numerics-investigator's audit.
- Phase G Step 2 has NOT shipped. Brief deliverables NOT completed:
  - L0 streaming-equilibrium gate NOT promoted (still in diagnostics).
  - Phase E sentinel `xfail-strict` NOT removed.
  - 6 curvilinear snapshots NOT regenerated.
  - ERR catalog NOT updated.

These remain blocked pending architectural guidance.

## Pointers

- This memo: `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_replan_blocker.md`
- Diagnostic memo: `.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_diagnostic.md`
- Step 2 attempted diff: `/tmp/step2_separate_pass_attempt.diff`
- Replan: `.claude/plans/issue_196_phase_g_replan.md`
- Phase F closeout: `.claude/agent-memory/method-implementer/issue_168_phase_f_closeout.md`
- L0 diagnostic test: `tests/sn/diagnostics/phase_g_step2_05_homogeneous.py`
- Spherical sweep (canonical target): `orpheus/sn/sweep.py:397-595`
- Spherical apply matvec: `orpheus/sn/operator.py:571-838`
- M-M strategy: `orpheus/sn/spatial/pole_angular_closure.py:340-458`
- Carlson seed: `orpheus/sn/spatial/psi_half_angle_seed.py:363-419`
- Phase G replan diagnostic: `.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_diagnostic.md`

## Linked memories

- `[[issue-196-phase-g-step1-replan-closeout]]` — Step 1 success (strategy-layer residual).
- `[[issue-168-phase-f-closeout]]` — Phase F partial closure that left this Step 2 work.
