---
name: vacuum-bc-eigenvalue-divergence
description: SN vacuum-BC eigenvalue test failure on refactor/sn-operator-algebra — partial investigation memo (execution environment blocked Bash/Write so root cause could NOT be fully proved). Hypothesis cascade narrowed to two surviving suspects.
metadata:
  type: project
---

# SN vacuum-BC eigenvalue divergence — partial investigation memo

**Date**: 2026-05-26
**Branch**: refactor/moment-space-and-layering @ 8e5a263
**Base failing branch**: refactor/sn-operator-algebra @ 62994ad
**Failing test**: tests/sn/test_boundary_conditions.py::TestSNBCSweepBehavior::test_vacuum_keff_lower_than_reflective
**Investigator constraint**: Bash execution and Write permissions BOTH blocked mid-session. No pytest run, no reproducer file written, no diagnostic cascade promoted. Investigation is **STATIC ANALYSIS ONLY**. Root cause NOT proved; hypotheses narrowed.

## Failure signature (per briefing — UNVERIFIED in this session because pytest not runnable)

| BC | keff | n_outer | flux magnitude | converged flag |
|---|---|---|---|---|
| reflective | 1.8749999956 (≈ k_inf) | 3 | O(1) | True (correct) |
| vacuum | 1.9365092901 (> k_inf!) | 64 (cap) | ~1e-43 (denormalised FP) | True (FALSE — hardcoded) |

Iteration history: ~1.669 → 1.875 → 1.936 (non-monotonic growth THROUGH k_inf — diagnostic smoking gun).

## Static analysis — what was verified

### k_inf reference (analytical)
Material A_2G from `orpheus.derivations.common.xs_library`:
- σ_t = [0.50, 1.00], σ_a = [0.02, 0.10], νσ_f = [0.025, 0.20], chi = [1.0, 0.0]
- σ_s = [[0.38, 0.10], [0.00, 0.90]] (downscatter only, no upscatter)
- A = diag(σ_t) - σ_s^T = [[0.12, 0], [-0.10, 0.10]]
- F = chi ⊗ νσ_f = [[0.025, 0.20], [0, 0]]
- A⁻¹·F dominant eigenvalue = **1.875** ✓ matches briefing's k_inf

### Per-group asymptotic eigenvalues (flux-shape bounds)
- k_g0_only = νσ_f[0]/σ_a[0] = 0.025/0.02 = **1.25** (g0-pure-spectrum limit)
- k_g1_only = νσ_f[1]/σ_a[1] = 0.20/0.10 = **2.00** (g1-pure-spectrum limit)
- k_uniform_spectrum = (0.025+0.20)/(0.02+0.10) = 1.875 (= k_inf, the geometric mean for chi=[1,0] feed)

**Critical observation**: 1.25 < 1.875 < 1.936 < 2.00. The 1.936 reading is BELOW the g1-only asymptote, so it's a mathematically reachable ratio IF the iteration is degenerating toward a g1-dominant spectral shape. The flux at ~1e-43 confirms degeneration.

### BC code path
Resolved through static reading of:
- `orpheus/sn/solver.py:881-1006` (solve_sn entry)
- `orpheus/sn/solver.py:446-562` (_solve_source_iteration — R-1 Step E carve)
- `orpheus/numerics/iteration.py:243-451` (SourceIteration primitive)
- `orpheus/sn/operator.py:~2530-2700` (InvertibleOperator.solve)
- `orpheus/sn/sweep.py:99-540` (transport_sweep, slab joint-batch fast path)
- `orpheus/sn/boundary_realizer.py:106-196` (SNBoundaryRealizer.realize)
- `orpheus/numerics/operator.py:920-1010` (IncomingOrdinateMaskTensor — vacuum realised op)
- `orpheus/numerics/trace_space.py:140-244` (inflow_indices_for_face)

**VERIFIED CORRECT**:
- Vacuum BC realised op is `IncomingOrdinateMaskTensor(inflow_indices, axis=0)`.
- For "left" face, inflow_indices = ordinates with μ_x > 0 (sign·μ < -eps with sign=-1).
- For "right" face, inflow_indices = ordinates with μ_x < 0.
- `apply(face_buffer)` returns a copy with inflow indices zeroed → correct vacuum inflow trace.
- The slab joint-batch sweep (sweep.py:465-535) reads `inflow_left[forward_ords]` for the forward chain start; this is zero under vacuum (correct).
- The persistent `solver._boundary_flux` is zeros at init and NEVER written by inner SI (only `boundary_buf` is, which is scoped to one `L.solve` call). Threading through outer iterations happens via `self._psi_typed.boundary`.

**VERIFIED CONVENTION-CLEAN** (post-A1 commit de8822d, ERR-049 fix):
- `q_ext_per_ord.values` = `Q_iso / sum_w` broadcast to N (per-ordinate density).
- `transport_sweep` does NOT apply `/W` internally (the pre-A1 `weight_norm = 1/W` line is GONE — verified sweep.py:393-398 `QV_per_ord = Q_per_ord * V` only).
- L1 `tests/sn/test_invertible_operator.py::TestInvertibleSolveBridgeRegression` confirms the SI carve gives correct k_inf for REFLECTIVE on slab/sphere/cylinder × {2eg, 4eg}.

### Convergence flag is COSMETICALLY wrong
`solve_sn` at solver.py:992-996 sets `history = IterationHistory(..., converged=True)` UNCONDITIONALLY. It does not consult `solver.converged(...)`. So the briefing's "n_outer=64 (cap hit) but converged=True" is the hardcoded True, NOT a false-positive convergence-criterion firing.

The `solver.converged()` method at solver.py:432-442 uses `dphi = norm(flux - flux_old) / max(norm(flux), 1e-30)`. With denormalised flux (norm ~1e-43), the floor `1e-30` kicks in and `dphi ~ 1e-43/1e-30 = 1e-13 < flux_tol`. So even if `power_iteration` consulted `converged()`, it would falsely declare convergence on denormalised values. This is a separate bug latent in the convergence criterion.

## Surviving hypotheses (NOT ruled out — execution needed)

### H1 (LEAD CANDIDATE) — Power iteration finds a spurious g1-dominant mode

The operator `M = (L_vac + C - S)⁻¹ · F` should have dominant eigenvalue = true k_eff (well below 1 for a 2cm vacuum slab with this material — diffusion-theory estimate ~0.05). Power iteration with `flux_{n+1} = M · flux_n / k_n`:
- If k_n converges to k_true ≈ 0.05, flux amplitude stabilises.
- If k_n is overestimated (k_n > k_true), flux SHRINKS each iteration → toward denormalised.
- As flux shrinks differentially (g0 leaks faster than g1), spectral ratio flux_g1/flux_g0 GROWS.
- As spectral shape degenerates toward g1-only, computed k = production/absorption drifts toward 2.00.
- The trajectory ~1.669 → 1.875 → 1.936 → 2.00 IS the signature of progressive g0-mode loss.

**Why does k_n start at ~1.67 and not at the true ~0.05?** This is the smoking gun. With initial `flux = ones`, the first compute_keff gives a value reflecting the initial spectrum, NOT the true eigenvalue. Power iteration would normally GROW the dominant-eigenvector component while k_n converges. But for some reason the iteration is rejecting the true small-k physical fundamental mode and amplifying a g1-trapped spurious mode.

**Reason this is plausible**: The within-group SI loop iterates `psi_{m+1} = (L+C).solve(S·psi_m + q_ext)`. The convergence rate is the spectral radius of `(L+C)⁻¹·S`. For chi=[1,0] material:
- g0 self-scatter ratio c_g0 = 0.38/0.5 = 0.76 (reflective). For vacuum on 2cm/2mfp slab, effective `(L+C)⁻¹` shrinks → c_g0_effective smaller.
- g1 self-scatter ratio c_g1 = 0.90/1.0 = 0.90 (reflective). For vacuum on 2cm/2mfp slab, effective c_g1_eff = ~0.7-0.8 perhaps.

Inner SI should converge with these rates × 200 iters → residual ~ 1e-19 (well past tol). NOT a convergence-rate issue.

But the OUTER power iteration's dominance ratio depends on `|k_1/k_0|` of M. For sub-critical vacuum, k_0 = 0.05 and k_1 could be... well, what are the other eigenvalues of M = (L+C-S)⁻¹·F for this configuration? With chi=[1,0] and downscatter-only, F is rank-1 (the production goes to g0 only via chi=1). So `M = (L+C-S)⁻¹ · F` is rank-1 → has exactly ONE nonzero eigenvalue. The dominant eigenvalue IS k_eff and there's no "k_1" to compete with. Power iteration should converge in 1 iter.

**Hmm — that argues against H1.** If M is rank-1, power iteration on it would converge instantly. So something MUST be making M effectively higher-rank from the SI solver's perspective.

### H2 (CONTENDER) — `_psi_typed` warm-start carries inconsistent state across vacuum outer iterations

`_solve_source_iteration` at solver.py:554 grabs `initial_guess = getattr(self, "_psi_typed", None)`. This is the PREVIOUS outer iteration's converged psi. It is fed as `initial_guess` to `SourceIteration.solve`. Inside the SI loop, `psi_prev = psi = initial_guess` on first inner iter. The SI then runs.

For REFLECTIVE, the initial_guess.boundary carries the partner-flux state from the previous outer step's converged sweep — load-bearing for fast convergence. The state is self-consistent across outer iterations because the underlying eigenvalue equation's reflective-fixed-point is the same regardless of source magnitude.

For VACUUM, initial_guess.boundary carries forward-ords outflow at xmax and backward-ords outflow at xmin from previous outer step. This is NOT actually needed (vacuum doesn't couple opposite faces), but it shouldn't hurt either — the BC apply zeroes the inflow ords on each iter.

**BUT** — `initial_guess.values` carries the FULL angular flux from the previous outer step's converged psi. The new outer step has a DIFFERENT fission source (because keff has changed and flux has changed). The SI then iterates with the warm-start.

If the previous outer's psi had magnitude `||psi_n||` and the new source magnitude is `||F·phi_n/k_n||` such that the new converged psi should have magnitude `||M·phi_n/k_n||`, the SI iteration drives psi from `psi_n` to the new fixed point. Each SI iter: `psi_new = (L+C).solve(S·psi_old + q_ext)`. The new converged psi is `(L+C-S)⁻¹·q_ext`. SI converges geometrically with rate `c = ρ((L+C)⁻¹·S)`.

If c is close to 1, convergence is SLOW. For our material, c is at most ~0.9. So 200 inner iters → residual ~`(0.9)^200 = 7e-10`. Tight enough.

**This H2 seems also weak.** The warm start might add overhead but shouldn't cause systematic drift.

### H3 (NEW — most consistent with the symptom) — The inner SI's `compute_fission_source` is being computed from the WRONG flux during the iteration sequence

Look at power_iteration carefully (eigenvalue.py:143-153):
```python
for n in range(1, max_iter + 1):
    flux_old = flux.copy()
    keff_old = keff
    fission_source = solver.compute_fission_source(flux, keff)
    flux = solver.solve_fixed_source(fission_source, flux)
    keff = solver.compute_keff(flux)
    ...
```

`fission_source = compute_fission_source(flux, keff)` uses the CURRENT flux at the START of the outer iter (before the within-group solve). Then `solve_fixed_source` updates flux. Then `keff = compute_keff(NEW flux)`.

The DIVISION by keff in `compute_fission_source` is `F·flux/keff`. If keff is overestimated (e.g. starts at 1.0 vs true 0.05), then `fission_source/keff` is UNDERESTIMATED by 20×. The within-group SI then converges to a small-magnitude psi. Then `compute_keff(small_psi) = small_production/small_absorption = ratio depending on spectral shape`.

For the FIRST outer iter: flux=ones, keff=1.0. fission_source = F·ones/1 = `chi·(0.025+0.20)·V_per_cell` = uniformly distributed in g0. The within-group SI on vacuum with this source → some converged psi. Then compute_keff(psi).

The first computed keff is ~1.669 (per briefing). This is the ratio of νΣ_f-weighted to Σ_a-weighted converged flux. If the converged flux has spectral shape close to "all g0 with some g1 from downscatter", k could land between 1.25 (g0-pure) and 2.00 (g1-pure).

Second outer iter: flux = converged psi from iter 1, keff = 1.669. `fission_source = F·flux/1.669`. SI again. compute_keff again. → 1.875.

Third outer iter: keff = 1.875. SI again. compute_keff → 1.936.

The k_n trajectory is approaching 2.00 asymptotically. This is the signature of **power iteration on a rank-1 operator where the iteration's projection is being dominated by the dominant g1-only mode** because the operator's "production" pathway is effectively f1-amplifying.

**Why would the operator's dominant eigenvalue be ~2.00 instead of ~0.05?** The TRUE physical dominant eigenvalue is 0.05 because vacuum leaks. But the SI is solving `(L+C-S)·psi = F·phi/k` for psi. If `(L+C-S)` is being INCORRECTLY assembled — e.g., if there's a sign error in S, or if the streaming term isn't actually being applied, or if the BC isn't really being applied — then the equation has different eigenvalues.

### H4 (NEW — what to verify FIRST in a follow-up session)

**Run `solve_sn_fixed_source(materials, mesh_vac, quad, external_source=Q_uniform, boundary_condition="vacuum")` with a uniform g0 source AND compare the result to a hand-computed vacuum slab flux.**

This isolates the WITHIN-GROUP path from the power-iteration. If the fixed-source result is CORRECT (matches analytical leakage), then the bug is in the OUTER iteration. If it's WRONG, the bug is in the within-group SI on vacuum BC.

Specifically, set `external_source = (1.0, 0.0)` per cell (only g0 source). The expected steady-state flux is the diffusion-theory solution `phi_g0(x) = A·cosh((x-1)/L)/cosh(1/L)` (slab geometry with vacuum BC) where L is the diffusion length, plus g1 built from downscatter. Compare to SN result.

## Minimal reproducer (to be written when Write/Bash permissions are restored)

```python
"""Diagnostic: SN vacuum-BC eigenvalue divergence — minimal reproducer.

Save as derivations/diagnostics/diag_vacuum_bc_eigenvalue_divergence.py
"""
from __future__ import annotations
import warnings
import numpy as np
import pytest

def _build_2eg_slab(boundary, n_cells=20, length_cm=2.0):
    from orpheus.derivations.reference_values import get
    from orpheus.geometry import (
        BC, Mesh1D, Region, RegionMesh, StructuredGeometry,
    )
    from orpheus.numerics.quadrature import Quadrature
    case = get("sn_slab_2eg_1rg")
    mix = next(iter(case.materials.values()))
    materials = {0: mix}
    geom_mesh = Mesh1D.from_geometry(
        StructuredGeometry(
            geometry="SLB",
            regions=(Region(mat_id=0, outer_thickness_cm=length_cm),),
            bcs=(boundary, boundary),
        ),
        region_meshes=(RegionMesh(n_cells=n_cells),),
    )
    quad = Quadrature.gauss_legendre(4)
    return materials, geom_mesh, quad

def test_vacuum_bc_keff_below_kinf():
    from orpheus.geometry import BC
    from orpheus.sn.solver import solve_sn
    materials, mesh, quad = _build_2eg_slab(BC.vacuum)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        result = solve_sn(materials, mesh, quad)
    # k_inf = 1.875 (analytical from A_2G).
    # Physical vacuum keff for 2cm slab << 1.875.
    print(f"\n[diag] vacuum keff = {result.keff:.10f} (expected < 1.875)")
    print(f"[diag] keff_history (last 10): "
          f"{[f'{k:.6e}' for k in result.history.keff_history[-10:]]}")
    print(f"[diag] terminal psi magnitude: "
          f"max={float(np.max(np.abs(result.angular_flux.values))):.3e}")
    print(f"[diag] terminal scalar flux: "
          f"max={float(np.max(np.abs(result.scalar_flux.values))):.3e}")
    print(f"[diag] n_outer = {result.history.n_outer}")
    assert result.keff < 1.875, (
        f"vacuum keff = {result.keff:.10f} >= k_inf=1.875; "
        f"physically impossible (vacuum leakage MUST reduce keff)."
    )
```

## Next steps for a session with execution capability

1. **Reproduce the bug** by running the test above. Confirm keff > 1.875 and flux denormalised.
2. **Run with inner_solver="krylov"**. If Krylov gives the correct sub-critical keff, the bug is SI-only. If both fail identically, the bug is in the operator construction (BC application, fission, etc.) — NOT the inner solver. This binary split is the highest-leverage probe (~5 min effort, massive narrowing).
3. **Run `solve_sn_fixed_source(materials, mesh_vac, quad, ext_src=Q_g0_uniform, boundary_condition="vacuum")`**. If the resulting flux matches the analytical diffusion-theory solution for a vacuum slab with this material, the within-group vacuum-BC code path is healthy and the bug is in the OUTER eigenvalue iteration. If the flux is wrong, the bug is in within-group vacuum BC. (This is a Step-3 probe-cascade fixed-source diagnostic per the cascade skill.)
4. **Print at each outer step**: keff, ||flux||_inf per group, flux_g0/flux_g1 spatial ratio. If the spectral ratio is changing systematically (g1 growing), confirms the spectral-degeneration H1/H3. If the spatial shape is changing, confirms shape-leakage.
5. **Print the EIGENVALUES of `M = (L_vac + C - S)⁻¹ · F` directly** via scipy.sparse.linalg.eigs on the matvec. The dominant eigenvalue IS k_eff. If it lands at ~1.936, then the OPERATOR has that eigenvalue (real bug at operator construction). If it lands at ~0.05, then power iteration is failing to converge on the true dominant eigenvalue (bug in outer iteration, possibly a normalisation drift that amplifies a sub-dominant mode).

## Files referenced

- `tests/sn/test_boundary_conditions.py:135-167` (failing test)
- `orpheus/sn/solver.py:881-1006` (solve_sn entry, hardcoded `converged=True` at line 992-996)
- `orpheus/sn/solver.py:446-562` (_solve_source_iteration — R-1 Step E carve)
- `orpheus/sn/solver.py:432-442` (converged() criterion — denormalised-FP fooling latent bug)
- `orpheus/numerics/iteration.py:243-451` (SourceIteration primitive)
- `orpheus/numerics/eigenvalue.py:120-155` (power_iteration)
- `orpheus/sn/operator.py:~2530-2700` (InvertibleOperator.solve — Phase 1.1 / A1 post-ERR-049 form)
- `orpheus/sn/sweep.py:99-540` (transport_sweep, slab joint-batch fast path)
- `orpheus/sn/boundary_realizer.py:106-196` (SNBoundaryRealizer)
- `orpheus/numerics/operator.py:920-1010` (IncomingOrdinateMaskTensor — vacuum BC realised op)
- `orpheus/numerics/trace_space.py:140-244` (inflow_indices_for_face)
- `.claude/skills/vv-principles/error_catalog.md:3583-3760` (ERR-049 — the most-related prior bug)

## Cross-references

- [[r1-step-e-invertible-solve-w-bridge]] — ERR-049, the previous convention-drift in this code path. Fixed at de8822d (2026-05-22). Vacuum BC was never tested at the eigenvalue level after that fix; this failure may be a sibling defect.
- ERR-049 catalog entry (vv-principles/error_catalog.md:3583): explicitly notes "Slab-1eg-SI passed by accident" and "L0 streaming-equilibrium tests were missing for the SI path". Same gap class.

## Session limitations to flag to the user

This memo is shipped INCOMPLETE because the execution environment denied:
- Bash commands beyond `ls`/`grep`/`git log` (no `python`, no `pytest`, no `touch`).
- `Write` permission to `derivations/diagnostics/` (the reproducer script could not be saved to disk).

Without executable diagnostic evidence, the surviving hypotheses (H1, H3, H4) cannot be discriminated. The user (or a follow-up session with execution capability) MUST run the reproducer above before treating any of these hypotheses as confirmed.

The diagnostic cascade Step 1 (reproduce minimally) is **DESIGNED** in this memo but **NOT EXECUTED**. Per `vv-principles` §"Write runnable evidence" and `numerics-investigator` rule 5 ("Every claim must have a script that proves it"), this investigation is **OPEN**, not closed.
