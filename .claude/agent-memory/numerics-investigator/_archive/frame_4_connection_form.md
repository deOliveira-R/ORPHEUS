---
name: Frame 4 connection-form interpretation FALSIFIED
description: Cross-domain attack Frame 4 (differential-geometry connection form) predicted F.4 = G_L · M^T · (I − M·W_L·M^T)^{-1} · M · P_L at rank-N. Tested at rank-1, rank-2, rank-3 and across σR = 0.1..100. Identity FAILS at every rank and every regime. Root cause: W does NOT transform as a (1,1) tensor under the Lambert↔Marshak basis change.
type: project
---

# Frame 4 — differential-geometry connection-form interpretation FALSIFIED

## Context

Cross-domain attack (Session 2026-04-22) proposed six frames to reinterpret
the rank-N white-BC barrier. Frame 4 asserted that F.4's Lambert-P/G +
Marshak-W mismatch is **parallel transport** between the emission (Lambert)
and transmission (Marshak) frames of the half-range line bundle — a
**gauge-fixed covariant operator**. Concrete first test:

    K_bc^{F.4} = G_L · M^T · (I − M · W_L · M^T)^{-1} · M · P_L

at rank-2, where M is the upper-bidiagonal change-of-basis matrix. If
bit-exact, F.4 is frame-covariant and the rank-1 "accident" is just the
scalar case (M = √2/2 at N=1).

## Outcome: FALSIFIED at rank-1 already

Tested the symmetric identity and four variants (left-only M, right-only M,
reverse M^T/M, no-conjugation) at:

- Rank-1: r_i/R = 0.3, 0.7, σR = 10 — rel diff to F.4 ≥ 25% for every
  variant. The closest (no_conjugation = G_L · (I − M·W_L·M^T)^{-1} · P_L)
  is 4.7% off.
- Rank-2: r_i/R = 0.3, 0.5, 0.7, 0.9, σR ∈ {0.1, 1, 5, 10, 20, 50, 100}
  — rel diff ranges 13%–44% for every variant at every point. No
  variant ever matches bit-exact or within 1%.
- Rank-3: r_i/R = 0.5, σR = 10 — rel diff 10%–21% for every variant.

The ONLY "matches" at σR = 100 are spurious (all values ~10⁻⁴⁷, below
double-precision floor — the printed "<<< MATCH" is numerical zero / zero).

## Structural reason

W does **NOT** transform as a (1,1) tensor under M:

    W_M ≠ M · W_L · M^T    (verified symbolic + numeric at rank-2)

At rank-1 specifically: M = √2/2, W_L = 0.1, W_M = 0.005 (at σR=10), so
M·W_L·M^T = 0.05 while W_M = 0.005 — a 10× discrepancy. This rules out
ANY conjugation-based covariance story at every rank.

The correct algebraic relationship (discovered in the diagnostic) is
**asymmetric µ-multiplication**:

    W_M = µ̂_Lambert · W_L    (exactly, at infinite rank)

where µ̂_Lambert[m, n] = <P̃_m, µ P̃_n>_L = (M^T M)_{mn} = (B^µ)_{mn} is the
µ-multiplication operator in the Lambert ONB. At finite rank N this
holds exactly in rows 0..N-2 but has a truncation residual in row N-1
(because µ · P̃_{N-1} has a P̃_N component outside the rank-N basis).

Symbolic proof at rank-2 (SymPy):

    (W_M − µ̂ · W_L)[0, :] = [0, 0]                      (exact)
    (W_M − µ̂ · W_L)[1, 0] = √3(τ² e^{2τ}/3 − τ² /3 − τ e^{2τ} − τ +
                               e^{2τ} − 1)·e^{-2τ}/τ³    (non-zero)
    (W_M − µ̂ · W_L)[1, 1] = similar O(1/τ) term          (non-zero)

Both residuals → 0 as τ → ∞ (consistent with E7's thick-limit
observation that P_M/P_L → 1). But they do NOT vanish at any finite τ.

## Why the rank-1 scalar "gauge" actually works

At rank-1 there is no matrix structure — everything is a scalar. The F.4
closure `G_L · (1 − W_M)^{-1} · P_L` differs from the Marshak-consistent
closure `G_M · (1 − W_M)^{-1} · P_M` only by the scalar ratio
`G_L · P_L / (G_M · P_M)`, which can be absorbed into an effective
reflection albedo β_eff. Direction Q verdict (B) is **unchanged** — this
was always the correct interpretation. Frame 4's "covariance" framing is
a false promotion of the rank-1 accident to a rank-N structural identity.

## Consequence for Sphinx docs

The §`peierls-f4-rank-1-gauge-why` section should NOT be upgraded with a
Frame 4 connection-form narrative. The existing text (classification B —
"scalar gauge DOF that absorbs into β_eff") is correct. Optionally add a
falsification note:

> The connection-form reinterpretation (F.4 = conjugation of the Lambert
> closure by the change-of-basis matrix M) was tested numerically and
> symbolically at rank-2 (see `diag_frame_4_connection_form.py`) and
> FAILS at every (r_i, σR). W does not transform as a (1,1) tensor
> under M; the true relation W_M = µ̂ · W_L is asymmetric and exhibits
> polynomial-truncation leak at finite rank. F.4 has no rank-N
> frame-covariant lift.

## What would a "principled" F.4 alternative look like?

Since W_M = µ̂ · W_L is asymmetric, a pedagogical rewrite of F.4 could use
`(I − µ̂ · W_L)^{-1}` in place of `(I − W_M)^{-1}` — numerically equivalent
at rank-1, different at rank-N. But experiments already done (E2.4 —
rank-N Lambert-P/G catastrophe) show that plugging Lambert-basis
objects into the ORPHEUS rank-N machine gives 33–737% k_eff error,
confirming that the asymmetric truncation leak IS the obstruction
and cannot be absorbed. There is no production path.

## Diagnostic location

`/workspaces/ORPHEUS/derivations/diagnostics/diag_frame_4_connection_form.py`
(384 lines, 20 s devcontainer). Runs rank-1 and rank-2 sanity checks,
regime scan σR ∈ {0.1, 1, 5, 10, 20, 50, 100}, and structural residual
analysis. Self-contained — not suitable for promotion to tests (diagnoses
a falsified hypothesis, no regression target).
