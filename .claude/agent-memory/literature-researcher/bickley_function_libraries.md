---
name: Bickley function library landscape
description: Davierwalla 1982 vs Lorensi 2025 vs Amos 1983 vs current ORPHEUS ki_n_float — what each provides and the n-coverage / accuracy / speed trade-offs. Updated 2026-05-03 with Davierwalla-vs-Amos head-to-head verdict.
type: reference
---

# Bickley function library landscape

When the user asks "is there a faster Bickley?" or "should we replace
ki_n_float with X?", consult this matrix BEFORE diving into the
papers — saves a session of duplicated reading.

## Coverage matrix

| Library | Orders covered | Method | Accuracy (verified) | Speed (estimated) | Coefficients available |
|---------|----------------|--------|---------------------|-------------------|------------------------|
| **Lorensi 2025** | n = 1, 2 only | Polynomial+Remez per interval, asymptotic+Möbius for tail | few ULP (~1e-15) | ~200 ns / eval (C++) | [0,1] printed (Tables 1–2); [1,24] must be regenerated |
| **Davierwalla 1982** | s = 1, 2, 3 on [0, 7] + s = 13, 14, 15 on [7, 600], recurrence to other (n, x) inside strips A/B/C — **strip D = {x ≥ 28, s > 15} EXCLUDED** | Continued-fraction rational Chebyshev (8 sub-intervals on [0, 7]) | 1e-12 at CDC single precision (14 stored digits); UNVERIFIED at IEEE binary64 | ~50–100 ns / eval (C, branch-light, ~6–8 divs) | All printed but as scanned PDF galley — re-keying ~150 16-digit constants required |
| **Amos 1983 (TOMS Alg. 609 — DBSKIN)** | **All n × all x** up to underflow (XLIM ≈ 743 on IEEE binary64) | Power series on [0, 2] + uniform asymptotic expansion on (2, ∞) (precision-agnostic structure) | 1e-12 (single prec.) → 1e-15 (IEEE binary64); 1e-18 floor on CDC double from 18-digit stored constants | unknown wall-time; algorithm is mul–add dominated, ~50–200 FLOPs / eval, ADAPTIVE truncation | **Public Fortran source** (Netlib/toms/609 + SLATEC mirror); f2py-wrappable in 5 minutes |
| **ORPHEUS `ki_n_float`** | All n | scipy.integrate.quad on tanh-substituted integrand | ~1e-13 to 1e-15 in [0, 20]; degrades in deep tail | 78 µs / eval (Python) | n/a (numerical) |
| **ORPHEUS `ki_n_mp(dps=25)`** | All n | mpmath.quad tanh-sinh | bit-exact at requested dps | ~3 ms / eval | n/a |

## Use guide

- "Need to speed up cylinder T-matrix builder" (`peierls_nystrom/geometry.py:2566`)
  → Lorensi handles only Ki_1, Ki_2. Ki_3+ falls outside their scope.
  Either chain via Bickley recurrence (Eq. 9 in Davierwalla) or reach
  for Amos 1983.
- "Need bit-exact reference for an L0 test"
  → Stay on `ki_n_mp` (cost is fine for one-shot validation).
- "Need < 1e-13 accuracy in working regime τ ≲ 20"
  → `ki_n_float` is *already there*. The 1.95e-06 figure in any
  benchmark is purely the deep-tail (τ ~ 50) where the function
  is < 1e-23 and abs-error 1e-15 swamps relative.
- "Need accurate Ki_n at τ > 30"
  → Don't compute it from the integral. Use either Davierwalla
  forward-recurrence (Strip B/D in his Fig. 1) or the asymptotic
  expansion `Ki_n(x) ~ √(π/2x) · e^-x · [1 - (4n+1)/(8x) + ...]`
  directly. Lorensi's Möbius-Remez form on [24, 743] is also safe
  here (n = 1, 2 only).

## Atkinson product-Nyström independence

Atkinson product-Nyström (slab F_N hardening,
`fn_method/peierls_atkinson_nystrom.py`) operates on **E_1**
(scipy.special.exp1), NOT Bickley. The F_N moment floor at ≤ 1e-5
is set by the slab kernel discretisation, not by Bickley accuracy.
Replacing Bickley does NOT move that floor. Don't conflate the two.

## What scipy provides

`scipy.special.expn(n, x)` (E_n exponential integral) is C-coded,
~150 ns / eval, vectorised. **There is no Bickley function in
scipy.special** — the closest cousins are `k0`, `k1` (modified
Bessel of 2nd kind, NOT Bickley). This is why ORPHEUS has its own
`ki_n_*` family.

## Decision residue (2026-05-03 — Davierwalla vs Amos head-to-head, both PDFs read)

User asked whether Davierwalla is *technically* advantageous over
Amos (vs the prior session's logistics-only deferral). After reading
both PDFs in full:

- **Amos 1983 (TOMS 609) is technically SUPERIOR to Davierwalla 1982**
  on every axis except raw FLOP count (Davierwalla's continued-fraction
  is ~2–4× tighter in compiled code; this difference does NOT propagate
  through Python overhead).
- Specifically: Amos covers ALL (n, x) up to underflow (no strip-D
  blind spot); Amos has stored 18-digit coefficients vs Davierwalla's
  CDC-14-digit (1000× double-prec accuracy advantage); Amos is
  portable Fortran with public ACM TOMS distribution (zero
  transcription cost) vs Davierwalla's scanned PDF galley.
- **The user should NOT spend effort manually transcribing
  Davierwalla's tables.** That effort delivers a worse implementation
  by every metric except a ~2–4× compile-time FLOP factor that
  doesn't matter in Python.
- **Project-level recommendation still: ADOPT NEITHER YET** —
  current `ki_n_float` is adequate (1e-13 in working regime, not
  load-bearing). When Bickley becomes load-bearing, the pivot is
  **Amos** (f2py wrap of BSKIN/DBSKIN, or pure-Python translation
  of the published source). Davierwalla is **not the technical
  winner** — only memorably "faster on FLOPs alone."
- **Lorensi (n = 1, 2 only)**: still ADOPT for the cylinder kernel
  hot path if a Bickley replacement is undertaken. 200–400× speed-up
  at sub-ULP accuracy; printed coefs for [0, 1]; Remez regen is a
  one-shot SymPy task. Pairs naturally with Amos (Lorensi for
  n = 1, 2 fast path; Amos for n ≥ 3 fallback).
- **Atkinson F_N moment floor**: independent of Bickley work
  (uses E_1 not Bickley).
- **Filename caveat**: `Amos (1993)…pdf` is misnamed. The paper is
  Amos 1983 (TOMS 9(4), Dec 1983). Mention this if rebuilding the
  citation later.

Full head-to-head memo at
`/workspaces/ORPHEUS/.claude/scratch/bickley_function_review.md`
(section "Head-to-head: Davierwalla 1982 vs Amos 1983").
