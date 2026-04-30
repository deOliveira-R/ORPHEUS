---
name: Issue #125 — Anisotropic-BC rank-N discriminator (Frame 6)
description: Frame 6 FALSIFIED — aniso BC also plateaus at rank-2 and rank-N is WORSE than rank-1 aniso, not monotone down. Barrier not representation-theoretic.
type: project
---

# Issue #125 — Anisotropic-BC discriminator — FRAME 6 FALSIFIED

## Decision

**Frame 6 refuted.** The cross-domain-attacker hypothesis that the rank-N
white-BC plateau is representation-theoretic predicted rank-N **monotone
decrease** in |err| under anisotropic (cos²) re-emission. The empirical
result is the opposite: rank-N is **worse** than rank-1 on aniso BC at
every one of 6 points, and plateaus at rank-2 just like white BC.

**Why:** the rank-N closure topology fails in the same way on both BCs
— this shows the barrier is numerical / variational (frames 1, 2, 5),
not symmetry-based.

**How to apply:** stop investing in rank-N-on-white-BC closures. Direction Q
verdict **(B)** "lucky rank-1 algebraic accident" stands as the final
framing for F.4. Escalate to Frame 1 (Issue #127, Chandrasekhar H-functions)
or Frame 2 (Issue #126, Rayleigh-Ritz) next.

## Evidence

### Sanity bridge (kernel abstraction verified)

At rank-(1,1,1), white kernel, `build_B_split_basis` with `kernel=2`
reproduces the existing `diag_cin_aware_split_basis_keff.py` closure to
**delta = 2.22e-16** — bit-exact. The anisotropic variant is therefore
a pure physics swap (kernel `k(µ) = 4µ²` vs `2`), not an algebraic rewrite.

### White-BC plateau (matching prior results)

6-point × rank-(1,2,3,5) scan at MED quadrature (`n_panels=2, p_order=6, n_ang=32`):

| σ_t·R | ρ   | rank-1    | rank-2    | rank-3    | rank-5    | Δ(1→5)    |
|-------|-----|-----------|-----------|-----------|-----------|-----------|
|  5.0  | 0.3 | -1.3740%  | -1.0635%  | -1.0629%  | -1.0625%  | +0.312%   |
|  5.0  | 0.5 | -1.7944%  | -0.9836%  | -0.9828%  | -0.9824%  | +0.812%   |
| 10.0  | 0.3 | -1.0008%  | -0.8584%  | -0.8582%  | -0.8579%  | +0.143%   |
| 10.0  | 0.5 | -1.3603%  | -0.9680%  | -0.9679%  | -0.9679%  | +0.392%   |
| 20.0  | 0.3 | -0.4641%  | -0.3849%  | -0.3849%  | -0.3847%  | +0.079%   |
| 20.0  | 0.5 | -0.6584%  | -0.4531%  | -0.4531%  | -0.4531%  | +0.205%   |

Classic rank-N plateau signature: rank-1 → rank-2 improves modestly,
rank-2 → rank-5 flat to 1e-4. All errors negative (under-prediction).

### Anisotropic-BC plateau — Frame 6 refuted on plateau alone

Same scan under anisotropic cos² BC (`kernel k(µ) = 4µ²`):

| σ_t·R | ρ   | rank-1    | rank-2    | rank-3    | rank-5    | Δ(1→5)     |
|-------|-----|-----------|-----------|-----------|-----------|------------|
|  5.0  | 0.3 | -0.1095%  | +2.8552%  | +2.9284%  | +2.9287%  | +3.038%   |
|  5.0  | 0.5 | +0.8317%  | +3.7559%  | +3.8414%  | +3.8419%  | +3.010%   |
| 10.0  | 0.3 | -0.5456%  | +0.9689%  | +0.9797%  | +0.9799%  | +1.525%   |
| 10.0  | 0.5 | -0.3183%  | +1.0993%  | +1.1216%  | +1.1216%  | +1.440%   |
| 20.0  | 0.3 | -0.2865%  | +0.5385%  | +0.5397%  | +0.5398%  | +0.826%   |
| 20.0  | 0.5 | -0.2003%  | +0.5481%  | +0.5604%  | +0.5604%  | +0.761%   |

6/6 points: rank-1 → rank-2 (a) **flips sign** (− → +) and
(b) **magnitude grows** by a factor of 3-25×. Plateau at rank-2 is
tighter than on white (rank-2 = rank-5 to 1e-4).

Frame 6 predicted rank-N monotone decrease. Observed: rank-N is **worse**
than rank-1 at every point. **FALSIFIED.**

### MC reference (anchor, σ_t·R=10, ρ=0.3)

- MC aniso: k = **1.498590 ± 0.001209** (100 active × 5000 hist, 16.6 s wall).
  Signed err vs k_inf=1.5: **-0.094% ± 0.081%**. Statistically indistinguishable
  from k_inf (1σ).
- MC white (cross-check): k = 1.49872 ± 0.00282. Also consistent with k_inf.

Both BCs preserve the 1-group k_inf=1.5 eigenvalue physically (current-
preserving re-emission with isotropic scattering). The closure's 0.5-4%
residual under aniso rank-N is purely a closure modeling error, not a
property of the physics.

Precision target in Issue #125 was 0.05% at 10⁶ histories. At 0.081% precision
the reference is tight enough to separate the rank-N closure residual (0.5-4%)
from MC statistical noise by >5σ. **Provisional tag dropped** — the closure
residual is structurally above the MC CI everywhere.

### L19 protocol outcome at anchor

`assert_rank_n_structural_win(closure_fn=rank-5_aniso, f4_fn=rank-1_aniso,
point=(10,0.3), quads=[RICH, RICH+panels])`:

**RAISED S2**: closure does NOT beat F.4_aniso at every quadrature.

```
closure (rank-5 aniso) signed = (+0.949778%, +0.949316%)
F.4_aniso (rank-1)     signed = (-0.631900%, -0.643300%)
per-quad |closure| < |F.4|: [False, False]
```

Both quadratures: rank-5 |err| = 0.949% > rank-1 |err| = 0.632% (RICH) or
0.643% (RICH+panels). The **rank-1 closure wins at every quadrature**,
exactly the opposite of Frame 6's prediction.

Additionally, F.4_aniso (rank-1) magnitude GREW RICH → RICH+panels (0.632% →
0.643%), tripping S5 as well. But S2 is the load-bearing failure.

## What Frame 6 got wrong

Frame 6's prediction rested on the assumption that rank-N modes N ≥ 2 are
"ruled out by symmetry" from coupling to the white BC operator's image,
and that breaking that symmetry (aniso BC) would unlock them. The
empirical reality:

- The rank-N **basis modes do couple to the aniso BC kernel** — the B
  matrix gets non-trivial entries off mode-0 because the kernel
  decomposition `4µ² = Σ d_m P̃_m(µ)` has non-zero higher-mode coefficients.
- But those additional couplings **do not reduce the closure error** —
  they just shift the fixed point to a different wrong answer.
- So the rank-N truncation error is fundamentally a **variational / basis**
  problem (frames 1, 2, 4, 5), not a symmetry-obstruction problem.

The "spurious modes" in Frame 6's story were never spurious. They were
real modes that the closure simply does not solve well — under ANY BC kernel.

## Sharp next experiment

If Rodrigo wants to continue:
- **Frame 2 (Issue #126) is now the most promising**: the rank-N non-monotone
  behavior on BOTH BCs is the classical tell of a Galerkin projection
  without a variational principle. Rayleigh-Ritz on nested subspaces
  would restore monotone convergence by construction. The minimal Ritz
  test: build a **nested** rank-2 subspace (e.g., `span{1, µ}` in the
  Marshak inner product AND orthogonalized against {1} in Lambert sense
  via the M matrix from Direction Q), re-evaluate, check rank-2 is
  strictly better than rank-1.
- **Frame 1 (Issue #127)**: the Chandrasekhar H-function basis diagonalizes
  the τ=∞ transfer operator; perturbation from this basis (not from
  Legendre/Jacobi) should give geometric convergence rather than plateau.

Both paths explain the plateau structurally AND prescribe the fix.

## Artifacts

- `derivations/diagnostics/diag_anisotropic_bc_rank_n.py` — full scan
  + MC + pytest tests (sanity + Frame 6 falsification).
- `derivations/diagnostics/_issue_125_scan.json` — raw numerical results.

## Open questions

- **Why does aniso rank-1 err (~0.2-0.8%) BEAT aniso rank-N (~0.5-4%)?**
  The rank-1 subspace in split basis is the scalar (µ_crit, ρ) direction,
  which is the F.4-equivalent Marshak mode-0 — the ONE direction where
  the Lambert/Marshak mismatch happens to not bite (see Direction Q verdict (B)).
  At rank-1 the kernel swap `2 → 4µ²` just recalibrates the scalar gauge
  α(τ,ρ) at an operating point that still happens to be close to the
  truth. Adding modes breaks this accident. This is further evidence for
  Direction Q's "rank-1 accident" framing — it works at rank-1 under
  BOTH kernels, for distinct algebraic reasons, but the accident does
  not extend to rank-N under either.

- **Would Frame 6 be vindicated for stronger symmetry breaking?** E.g.,
  if the BC were cos^k for large k (angularly sharper), would rank-N
  eventually beat rank-1? Not pursued — the structural plateau at
  rank-2 in the current aniso is already the falsification signature.
