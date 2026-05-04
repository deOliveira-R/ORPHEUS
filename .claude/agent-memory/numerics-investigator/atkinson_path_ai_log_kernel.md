---
name: Atkinson product-Nyström for Path A.i log kernel
description: Wave 2-A Path A.i flux stall pinned to log-singular kernel diagonal truncation in plain-GL μ-quadrature, fixed by Atkinson product-Simpson with closed-form log primitives. Empirical convergence rate signature is err ~ log(n)/n, NOT 1/sqrt(n).
type: project
---

**Fact:** Wave 2-A Path A.i (slab Peierls power iteration with
plain Gauss–Legendre on the (z, μ) tensor product) was stalled at
~5–7 % vs KLL Path B reference. The bug is **silent truncation of
the divergent diagonal kernel** — the discrete μ-quadrature gives
`Σ_k w_k/μ_k ≈ 2 log(n_μ)` instead of `E_1(0+) = +∞`. Off-diagonal
the same μ-quadrature converges to `E_1(τ)` at machine precision.

**Why:** The legacy Path A.i was deliberately built with plain GL
to be procedurally distinct from KLL Path B; its existing tests
admit "honest tolerance ~5e-2" because the team knew the
log-singular kernel needed singularity-aware quadrature but
deferred. Tracked as ERR-036.

**How to apply:** When investigating a Path-A-style operator
with plain-GL convergence stalls in the 1–10 % regime:

1. Check the convergence-rate fingerprint. The signature
   `err·n/log(n) ≈ const` (across an n-doubling sweep) pins the
   bug to log-singular kernel diagonal truncation, NOT
   Schneider C^(0,1) endpoint singularity (which gives
   `err·sqrt(n) ≈ const`). The numbers from Phase 1.5: at
   z/a=0.95, n∈{32,64,128,256}: `err·n/log(n)` = 1.39, 1.17,
   1.05, 1.07 (stable); `err·sqrt(n)` = 0.85, 0.61, 0.45, 0.37
   (decaying).
2. The fix is Atkinson 1972/1997 §4.2 product-Simpson:
   integrate `∫ s^k log|t-s| ds` analytically against the
   piecewise-quadratic Lagrange basis, integrate the smooth
   remainder `R(τ) = E_1(τ) + log(τ)` with standard Simpson.
   New module `orpheus.derivations.continuous.fn_method.peierls_atkinson_nystrom`.
3. The literature memo at `scratch/derivations/peierls_log_singular/`
   predicted a **half-order Schneider rate (-1/2)**; empirically
   the dominant rate is **first-order with log correction (-0.9)**.
   The memo's recommendation (Atkinson product-Nyström) was right
   but for a different reason than the memo emphasised. The
   diagonal-truncation fix dwarfs the endpoint-regularity issue at
   practical n; graded mesh would help further but is not needed
   to hit the F_N moment floor.
4. Empirically, n_panels=64 reaches sup_err ~ 5e-4 (vs legacy
   5–7e-2 — 100× improvement); n_panels=256 reaches ~ 1e-5
   (F_N moment floor).
