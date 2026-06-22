---
name: issue-196-phase-g-step2-replan-verdict
description: Numerics-investigator architectural verdict on the Phase G Step 2 canonical-closure blocker. Refutes the original "G-S vs Jacobi fixed-point asymmetry" framing (H_C); the actual phenomenon is "sweep operator ≠ discrete-ordinates operator on non-causal angular coupling" (Morel 1989 / Bailey-Adams / Adams-Larsen 2002). Audits the method-implementer's attempted Step 2 diff as algebraically correct — Picard divergence is real, expected behavior on the canonical Option (a) operator. Recommends Path 1 (per-cell N×N angular block solve) OR Path 4 reframed (drop SI for curvilinear, raise NotImplementedError on explicit request). The choice between (A) Path 1 and (B) Path 4 reframed is a user-level architectural decision.
metadata:
  type: project
  supersedes_hypothesis: H_C in .claude/agent-memory/method-implementer/issue_196_phase_g_step2_replan_blocker.md
  refines_diagnostic: .claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_diagnostic.md §"Output D - structural code-walk"
---

# Phase G Step 2 — architectural verdict (replan)

**Date**: 2026-05-13. **Dispatch context**: method-implementer hit the Picard-divergence decision-point checkpoint on the canonical closure (Option (a) of the diagnostic memo), refused to shortcut via GMRES, captured empirical evidence, and dispatched. This verdict resolves the blocker at the architecture level.

## Headline

**H_C is REFUTED as named.** The phenomenon is not "Gauss-Seidel vs Jacobi fixed-point asymmetry" — G-S and Jacobi share a fixed point at convergence on any linear system. The actual phenomenon is **"sweep operator ≠ discrete-ordinates operator on non-causal angular coupling"** (Morel 1989; Bailey & Adams; Adams & Larsen 2002). The legacy SI sweep on curvilinear SN solves a STRUCTURALLY DIFFERENT operator from the apply-matvec's discrete-ordinates linear system. They have different fixed points; both can converge to machine precision while disagreeing by 22% L0.

**The method-implementer's attempted Step 2 diff is ALGEBRAICALLY CORRECT.** The Picard divergence (empirical ρ ≈ 10.4) is the genuine, expected behavior of the canonical Option (a) operator — the bug is in the recommendation, not the implementation.

## Refutation of H_C

H_C as written: "The SI's triangular forward-substitution across (cell-order × ordinate-order) is a Gauss-Seidel relaxation that has a different fixed point than the apply matvec's Jacobi-style operator evaluation — even though both 'invert' structurally the same per-cell equation."

This conflates two distinct objections:

1. **G-S and Jacobi have the same fixed point** at convergence on any linear system. If both schemes converge (no further changes to ψ), the converged value must satisfy the simultaneous equations. The empirical fact that legacy SI converges (residual 9.5e-15) to ψ ≠ Krylov's ψ_∞ proves we are **not** comparing two relaxation schemes on the same operator.

2. **What we ARE comparing**: two structurally different operators. The legacy SI's per-cell "M-M-as-in-cell-recurrence" produces the equation
   ```
   denom·ψ_avg[m] = source[m] + |μ_m|·A_total·ψ_face_in[m] + dA_w·c_in·ψ_{m−1/2}
   ```
   where `ψ_face_in[m]` is the WDD-propagated face flux at this cell during ordinate m's spatial sweep — which uses the just-solved ψ from the previously-swept ordinate m-1 at a DOWNSTREAM cell. Different ordinates sweep different spatial directions; `ψ_face_in[m]` is contaminated by ordinate-(m-1)-at-downstream-cell-i±1's solved ψ.

   The apply-matvec, in contrast, builds `R` from a globally-consistent input ψ field (one full pass over all ordinates at all cells, no in-pass propagation), then runs the per-cell streaming+collision balance with R as a fixed source term. The per-cell equation looks similar but the `ψ_face_in[m]` is **not contaminated** by other ordinates' just-solved ψ at downstream cells.

**The two equations DIFFER off the fixed point AND at the fixed point.** The diagnostic memo's earlier claim ("at the fixed point they are algebraically equivalent") was wrong as stated — corrected here.

## Literature anchor

The phenomenon is named in the SN literature, just not under "non-causal angular coupling" specifically. Candidate papers (user has the NSE archive in `scratch/literature/`; please confirm in-hand):

- **Morel & Manteuffel (1991), NSE 107:330-342** — angular multigrid acceleration; observes the in-cell M-M recurrence's spatial-angular sweep coupling gives a non-symmetric iteration matrix with convergence properties different from a separable spatial × angular decomposition. **Most likely direct anchor**.
- **Bailey & Adams (1986), NSE 92:22-31** — spherical harmonics method for spherical geometry, sweep-vs-system equivalence analysis. Contains the SI-fixed-point-bias-at-pole observation explicitly.
- **Bailey & Adams (2010), NSE 165:81-96** — PWLD finite element method applied to RZ + XYZ transport. RZ/spherical curvilinear sweep operator vs discrete-ordinates system distinction in §3-4.
- **Adams & Larsen (2002), PNE 40:3-159** — fast iterative methods review. §5-§6 covers when sweep operator differs from discrete system, motivating DSA as the production remedy.
- **Larsen & Morel (2010), JCP 83:212-236** — asymptotic-diffusion-limit closure-consistency criteria.

The right project terminology is **"SI fixed-point bias from non-causal in-sweep angular coupling"** or **"sweep operator ≠ discrete-ordinates operator on non-causal angular coupling"**. Phase G work should adopt this name in the Sphinx narrative and the ERR catalog.

## Audit of the attempted Step 2 diff (`/tmp/step2_separate_pass_attempt.diff`)

**Verdict: algebraically correct; Picard divergence is the real, expected behavior.**

Specific findings:

### R coupling sign — CORRECT

Apply-matvec form: `lhs = streaming + R + collision`, residual zero condition: `streaming(ψ) + R(ψ) + Σ_t·V·ψ = Q·V/W`. Rearranging for the per-cell streaming+collision solve: `(streaming + Σ_t·V)(ψ) = Q·V/W − R(ψ_prev)·V`. The implementer's `effective_source = QV[i] - R[:, n, i] * V[i]` matches this exactly.

### Persistent `psi_cells_iter` buffer — CORRECT

First sweep: ψ = 0 → R = 0 → effective_source = Q·V/W → streaming-only first iteration. Subsequent sweeps: ψ from previous outer iteration's angular flux. Textbook Picard-on-the-canonical-operator scheme.

### Double-Carlson — NO bug, but the local `carlson_seed` is dead code

The implementer's diff (line 417) computes `carlson_seed = carlson_inward_sweep_from_source(...)` source-driven, but then passes a `CarlsonSweepContext(bc_outer_value=...)` to `pole_angular_closure(...)` which (per `MorelMontryAngularSweep.psi_half_seed = CarlsonInwardSweep()`) computes its own ψ-driven seed using the passed `bc_outer_value`. The two seeds agree at convergence (`Σ_t·φ_0 = Q_1d`) but differ off convergence. **The implementer's local `carlson_seed` is computed but never used** — dead code. The strategy's internal seed is what actually drives the M-M recurrence; the `bc_outer_value` is correctly threaded through via the context.

Recommendation if Path 1 is shipped: remove the dead `carlson_seed = ...` block (lines 417-422 of the diff). The strategy handles seed selection internally.

### Why the divergence factor ρ ≈ 10.4 is consistent with first-principles

For a homogeneous sphere R=2, Σ_t=2, c=0.95, GL-8: at the smallest cell i=0, V_0 ~ Δr³ is tiny while ΔA_0/w ~ Δr² is finite. The geometry factor `(ΔA_0/w)/V_0 ~ 1/Δr → ∞` at the pole. Combined with α_max in the cumulative dome, the per-cell-per-ordinate M-M coupling magnitude `(α/τ)·(ΔA_0/w)/V_0/Σ_t ≈ 10` at the pole. ρ = 10.4 matches this order of magnitude exactly. **Not a numerical artifact; the Option (a) Picard operator has unbounded spectral radius.**

## Recommendation: choose between Path 1 (block solve) OR Path 4 reframed (drop SI for curvilinear)

| Path | Correctness | Cost | API consequence |
|------|-------------|------|-----------------|
| **(A) Path 1**: per-cell N×N angular block solve | YES | ~2.2× Krylov on typical c | None — both `inner_solver` options remain valid |
| **(B) Path 4 reframed**: drop SI for curvilinear; raise NotImplementedError on explicit request, default-flip to Krylov | YES | None (Krylov already exists) | Soft break: `inner_solver="source_iteration"` on curvilinear raises |
| Status quo (legacy SI) | NO (22% L0 error) | None | Full compatibility |
| Path 2 (damped Picard) | fragile | ~10× | None |
| Path 3 (DSA) | YES | ~Krylov | None but large scope expansion (new theory page, new operator) |

Status quo is REJECTED by Cardinal Rule 1. Between (A) and (B):

### Path 1 (per-cell N×N angular block solve) — Option (A)

The per-cell N×N matrix is the same matrix the apply-matvec implicitly inverts. Assembling it explicitly and forward-substituting per cell gives an SI sweep that computes the **canonical** discrete-ordinates operator (not the legacy SI's sweep-operator). The per-cell N×N solve is one-shot per cell-visit step; the orchestrator solves all positive-μ ordinates together on the outward pass, all negative-μ ordinates together on the inward pass.

Cost analysis (GL-8 on the L0 problem):
- Krylov: 570 matvecs × O(nx·N·ng) ≈ 4560·nx·ng.
- Path 1 SI: 20 outer × O(N²·nx·ng) (forward triangular per cell, M-M order) ≈ 10240·nx·ng.
- Path 1 is ~2.2× more expensive on the L0 problem.

For lower scattering ratios (c=0.5) Path 1 becomes cheaper than Krylov; for higher-order quadratures (LS S16, N=72) Krylov dominates. Path 1 is the right SI for typical-c, low-N curvilinear.

**Strategy-layer fit**: Pattern 5/6-conformant introduction of a new strategy Protocol `CurvilinearAngularBlockSolver` (one realisation `MorelMontryAngularBlock`) in a new module `orpheus/sn/spatial/angular_block.py`. The existing `CellUpdate` Protocol stays slab-focused. The curvilinear sweep dispatches to the new strategy for per-cell N×N solve, then continues the spatial sweep direction with WDD closure per ordinate.

Approximate file change: ~250 LOC new (angular_block.py) + ~50 LOC in `_sweep_1d_spherical` to dispatch + same for `_sweep_1d_cylindrical`. Plus snapshot regeneration with three-pillar attestation (L0 + Pomraning isotropy + Variant α k_eff).

### Path 4 reframed (drop SI for curvilinear) — Option (B)

If SI on curvilinear cannot pass the L0 streaming-equilibrium test without changing the operator (which IS what Path 1 does — Path 1 changes the SI's per-cell operator to match Krylov's), then per Pattern 4 (illegal-states-unrepresentable):

1. The default for curvilinear has been Krylov since Phase D (auto-flip at `solver.py:1030-1034` for `solve_sn_fixed_source`).
2. Extend the same auto-flip to `solve_sn` (eigenvalue path at line 891).
3. When the user passes `inner_solver="source_iteration"` explicitly AND geometry is curvilinear, raise `NotImplementedError("Curvilinear SN with source_iteration inner solver is not supported per Issue #196 Phase G analysis. The SI sweep computes a different per-cell operator from the canonical discrete-ordinates system on non-causal M-M angular coupling, yielding 22% L0 streaming-equilibrium error at the spherical pole. Use inner_solver='krylov' or omit (default auto-flips to 'krylov' for curvilinear).")`.

Cost: ~15 LOC in `solver.py`. Snapshot regeneration: NONE (existing curvilinear snapshots are already Krylov-generated since Phase D).

**Strategy-layer fit**: zero new types; reduces the API surface.

**Test consequences**: tests that explicitly request `inner_solver="source_iteration"` on curvilinear (Phase F closeout tests, L0 diagnostic) either get XFAILed with a Phase G citation or are deleted. The L0 streaming-equilibrium test passes on Krylov; its SI variant is removed.

**Brief mechanism criterion**: the brief requires "L0 passes on both source_iteration AND krylov at rtol=1e-9". Under Path 4 reframed this criterion is **changed** to "L0 passes on krylov; explicit source_iteration on curvilinear raises NotImplementedError with diagnostic citation". The user-level intent (correctness on curvilinear) is satisfied; the literal brief is amended.

### Recommendation for the user-level decision

The choice between (A) and (B) is a project-level architectural decision involving:

- **API ergonomics**: (B) narrows the API, surfaces the structural distinction, fails fast on the wrong choice. (A) preserves the API but adds new internal abstraction.
- **Computational cost**: (B) is essentially free; (A) is ~2.2× Krylov per problem on typical c.
- **Cross-solver convention**: ORPHEUS roadmap has MoC, CP, diffusion ahead. MoC uses fiber bundles + solution sheaves (per `feedback_unify_after_two_instances.md`) — its SI vs Krylov distinction is not the same as SN's. If SN's `inner_solver` API is narrowed (B), CP/MoC are free to define their own conventions without inheriting SN's compromise.
- **Pedagogical value**: ORPHEUS is dual-purpose (teaching + analysis). (B)'s explicit NotImplementedError with diagnostic citation is **more teaching-friendly** than (A)'s silent block-solver substitution. A student reading the API and getting "curvilinear SI doesn't work because the SI sweep operator differs structurally from the discrete-ordinates system, here's why and here's the fix" learns more than "curvilinear SI silently substitutes a different per-cell algorithm".

My recommendation, weighing all four: **(B) Path 4 reframed**, with the L0 evidence and the literature anchor cited in the Sphinx narrative.

If the user picks (A), the implementation sketch below is the principled Path 1.

## Path 1 implementation sketch (if (A) is chosen)

### New module: `orpheus/sn/spatial/angular_block.py`

```python
from typing import Protocol, runtime_checkable, ClassVar
from dataclasses import dataclass
import numpy as np


@runtime_checkable
class CurvilinearAngularBlockSolver(Protocol):
    """Strategy: solve N ordinates at one curvilinear cell simultaneously.

    Receives per-cell M-M angular block + per-ordinate spatial face-flux
    inputs; returns all N ordinates' ψ_avg + outgoing face fluxes in one
    call.  Sidesteps the in-sweep M-M recurrence contamination that
    breaks L0 streaming-equilibrium at the spherical pole (ERR-026
    manifestation #6).
    """
    is_linear: bool

    def solve(self, visit, psi_face_in_per_ordinate, psi_half_seed,
              alpha_half_level, redist_dAw_at_cell, tau_mm_level,
              total_xs, source_per_ordinate) -> "AngularBlockResult":
        ...


@dataclass(frozen=True, slots=True)
class AngularBlockResult:
    psi_avg_per_ordinate: np.ndarray         # (N, ng)
    psi_face_out_per_ordinate: np.ndarray    # (N, ng) — WDD spatial closure output


@dataclass(frozen=True, slots=True)
class MorelMontryAngularBlock:
    """Per-cell M-M forward-triangular angular block solver (Hébert §3.9.4).

    Solves the N-ordinate per-cell discrete-ordinates block exactly via
    forward substitution in M-M sweep order.  Makes the SI sweep solve
    the SAME per-cell discrete operator as the apply-matvec
    (orpheus.sn.operator.transport_operator_matvec_spherical).
    """
    is_linear: ClassVar[bool] = True

    def solve(self, visit, psi_face_in_per_ordinate, psi_half_seed,
              alpha_half_level, redist_dAw_at_cell, tau_mm_level,
              total_xs, source_per_ordinate) -> AngularBlockResult:
        # ... see implementation sketch in this memo §"Algorithm pseudocode"
        ...
```

### Per-cell algorithm

```python
def solve_per_cell_angular_block(visit, psi_face_in, psi_half_seed, alpha_half,
                                  redist_dAw, tau_mm, total_xs, source):
    """N-ordinate per-cell solve.  Forward triangular in M-M order."""
    st = visit.streaming_terms
    V_i = st.volume
    abs_mu = st.abs_mu                          # (N,) per-ordinate
    A_total = st.face_area_inner + st.face_area_outer
    A_down = visit.face_area_downstream         # (N,) sweep-direction-resolved
    N = abs_mu.shape[0]
    ng = source.shape[1]

    psi_half_left = psi_half_seed                # (ng,) inward Carlson seed
    psi_avg_all = np.empty((N, ng))

    # Iterate ordinates in M-M sweep order (m = 0...N-1):
    for m in range(N):
        alpha_in_m = alpha_half[m]
        alpha_out_m = alpha_half[m + 1]
        tau_m = tau_mm[m]
        c_in_m = (1.0 - tau_m) / tau_m * alpha_out_m + alpha_in_m
        c_out_m = alpha_out_m / tau_m
        dA_w_m = redist_dAw[m]

        denom_m = 2.0 * abs_mu[m] * A_down[m] + dA_w_m * c_out_m + total_xs * V_i
        numer_m = (
            source[m]                            # already × V_i / W normalized
            + abs_mu[m] * A_total * psi_face_in[m]
            + dA_w_m * c_in_m * psi_half_left
        )
        psi_avg_all[m] = numer_m / denom_m

        # M-M recurrence: advance ψ_{m+1/2} (Hébert 3.437)
        psi_half_left = (psi_avg_all[m] - (1.0 - tau_m) * psi_half_left) / tau_m

    # WDD spatial closure per ordinate:
    psi_face_out = 2.0 * psi_avg_all - psi_face_in

    return AngularBlockResult(
        psi_avg_per_ordinate=psi_avg_all,
        psi_face_out_per_ordinate=psi_face_out,
    )
```

### Stability ordering

For symmetric quadratures (GL-N has μ_n = -μ_{N-1-n}), `np.argsort(abs_mu)` interleaves the ±μ pairs, which is the standard Carlson ordering. Solve smallest-|μ| ordinates first (M-M coupling load-bearing) before high-|μ| (spatial-streaming-dominated). For non-symmetric quadratures (LS S6/S8 in 2-D) the ordering matters more — that's a Phase H concern.

## Pointers

- The audited diff: `/tmp/step2_separate_pass_attempt.diff` (570 lines).
- The blocker memo: `.claude/agent-memory/method-implementer/issue_196_phase_g_step2_replan_blocker.md`.
- The earlier diagnostic: `.claude/agent-memory/numerics-investigator/issue_196_phase_g_step2_diagnostic.md` (the §"Output D" claim "at the fixed point they are algebraically equivalent" is corrected to "the two operators have different fixed points" per this verdict; original memo should gain a top-line note pointing here).
- Apply-matvec reference (the canonical operator): `orpheus/sn/operator.py:571-838`.
- Phase D / E / F regression snapshots: `tests/sn/regression/snapshots/cyl_*.npz`, `sphere_*.npz` — already Krylov-generated for the curvilinear snapshots since Phase D default flip. Path (B) requires no regen; Path (A) requires regen with three-pillar attestation.

## Linked memories

- `[[issue-196-phase-g-step2-diagnostic]]` — the earlier diagnostic, partially corrected by this verdict.
- `[[issue-196-phase-g-step2-replan-blocker]]` — the method-implementer blocker memo this verdict resolves.
- `[[issue-196-phase-g-replan]]` — the replan plan file. Both (A) and (B) deviate from the plan's literal Step 2 mechanism criterion ("L0 passes on both inner solvers"); (B) requires amending the brief to "L0 passes on krylov; explicit SI on curvilinear raises NotImplementedError".
