---
name: Stamm'ler & Abbate 1983 Chapter VI — negative finding on interface currents
description: Ch.VI of Stamm'ler-Abbate (1983) "Methods of Steady-State Reactor Physics in Nuclear Design" is the SN method chapter, not interface currents. Records what IS in the chapter, what is NOT, and where the interface-current material actually lives in the same book.
type: reference
---

# Stamm'ler & Abbate (1983) Chapter VI — scope clarification

## Citation

Stamm'ler, R.J.J. and Abbate, M.J., *Methods of Steady-State Reactor
Physics in Nuclear Design*, Academic Press, 1983. Chapter VI:
"The Discrete Ordinates or Discrete S_N Method", pp. 191-230.

PDF on disk: `/workspaces/ORPHEUS/Stammler(1983)Chapter6.pdf` (40 pp.).

## Negative finding

**Chapter VI is the S_N method chapter.** It contains no interface-
current method, no multifunction collision probability, no rank-N
per-face surface-mode expansion, no Stepanek-style anisotropic
interface currents, and no Benoist homogenization. The entire 40
pages are the multigroup discrete-ordinates equations, quadrature
sets, finite-difference schemes, boundary conditions, and 2D iteration
strategy — all SN, not CP / interface-current.

The references section (p.230) confirms the scope: Arkuszewski,
Carlson-Lathrop, Chandrasekhar, Lathrop, Lathrop-Carlson, Lee, Mika,
Wick. These are SN-method references exclusively — no Sanchez,
Bonalumi, Pellaud, Stepanek, Benoist, or any interface-current
literature.

## What Chapter VI actually contains (reconstructed ToC)

| § | Pages | Content |
|---|---|---|
| 1 | 191-194 | Multigroup discrete-ordinates equations, (Ω̂·∇ + Σ)ϕ_m = S_m + (χ/k)F, linearly-anisotropic scattering, current-weighted Σ_{s1} |
| 2 | 195-200 | Spatial discretization in spherical symmetry; angular redistribution coefficients a_{m±1/2} = ½(1-µ²)∫4πrϕ dr/∫ϕ dr |
| 3 | 201-210 | Quadrature sets: level-symmetric, moment-matching conditions, Gauss-Legendre in 1D, P_L ↔ S_{L+1} equivalence in slab |
| 4 | 211-213 | Boundary conditions: specular reflection (Eq. 26), white reflection (Eq. 27-28), albedo (Eq. 29-30) — all implemented by reassigning ϕ_m at the boundary, not by mode expansion |
| 5 | 213-223 | Finite-difference schemes: weighted (Eqs. 32), diamond-difference (a=½), step (a=1), optimal a(τ) from Eq. (36) |
| 5.2 | 217-219 | Negative-flux fix-up |
| 6 | 220-228 | 2D iteration sweep strategy, FORTRAN listings, ray effect |
| Q | 228-229 | 12 exercises |
| Refs | 230 | SN-method references only |

## SN-specific white BC (NOT rank-N)

For the record — Eq. (27c) on p.212 gives Stamm'ler's SN white BC:

    ϕ_m(x_a) = j⁺(x_a) / [½ Σ_m w_m |µ_m|]     for µ_m > 0

i.e. the returning angular flux is **constant in µ** and set equal
to the outgoing partial current divided by ½Σw|µ|. This is a
rank-0 (scalar) closure, just expressed over the SN quadrature. It
does **not** generalize to rank-N in the sense we need for Issue
#119 — there's no per-mode basis, no transmission matrix, and no
reflection operator.

## Where the interface-current material actually lives in Stamm'ler-Abbate

**Chapter IV** ("Integral Transport Theory; Collision Probabilities",
pp. 105-141) is the CP + IC chapter. Now confirmed — PDF
`/workspaces/ORPHEUS/Stammler(1983)Chapter4.pdf` is on disk and has
been fully extracted (see
`rank_n_closure_four_references_synthesis.md` in this directory for
the cross-reference).

**Important**: Stamm'ler Ch.IV uses a **scalar (rank-0)
cosine-return / white-BC** interface-current closure throughout.
Eqs. 29-37 build the closure from:

- a single integrated partial in-current j⁻ per face (Eq. 29)
- cosine (Lambertian) return distribution at white boundary
- scalar multi-reflection `Y_i(a) = Y_i / (1 − a(1−Γ))` (Eq. 34)
- multicollision blackness Γ = Σᵢ Σ_{r,i} V_i Y_i (Eq. 32)

**There is no Legendre-moment / DP_N / rank-N ladder in Stamm'ler
Ch.IV or Ch.VI.** Stamm'ler-Abbate does not treat rank-N IC at all.
The book's CP content stops at the scalar closure (which the
`subroutine COPRAN` on pp. 124-125 programs).

Contrary to the guess above, Ch. III is not the CP chapter
(Ch. III per the implicit chapter numbering is a different topic
and we have no PDF; the table of contents accessible to us shows
Ch. IV = CP/IC).

## Follow-up: where to actually get the rank-N per-face formulation

Based on the four-reference cross-check
(`rank_n_closure_four_references_synthesis.md`), the rank-N Legendre
ladder is NOT in Ligou 1982 Ch.8, Sanchez 2002 NSE 140, Stamm'ler
1983 Ch.IV, or Stacey 2007 Ch.9 — all four use scalar/DP-0. The
§III.F.1 ladder of Sanchez-McCormick 1982 appears to be a **lone
formulation** without independent cross-validation in the standard
reactor-physics textbook corpus.

Candidate sources still untested:

1. **Hébert (2009/2020)** *Applied Reactor Physics* Ch. 3 — modern
   textbook. Not yet extracted; highest-priority next read if a
   second rank-N derivation is needed.
2. **Stepanek (1982, 1984)** — anisotropic-interface-current author.
   Original papers in NSE/ANE.
3. **Bonalumi, Pellaud** — 1960s-70s primary rank-N literature.
4. **APOLLO/TDT code documentation** (Sanchez 2002 reference 12) —
   export-controlled, not accessible.

## How to apply

- If asked for "Stamm'ler Ch. 6 interface currents" in the future:
  this citation is wrong. Redirect to Stamm'ler Ch. III (CP method)
  or to Sanchez-McCormick 1982.
- The SN white BC formula Eq. (27c) can be cited when a discrete-
  ordinates white BC needs a textbook reference, but it is not
  relevant to the rank-N hollow-sphere closure for Issue #119.
- Do not confuse "constant-in-µ returning flux" (Stamm'ler SN white
  BC) with the rank-0 Lambert closure in ORPHEUS Phase F.4 — they
  are superficially similar (both enforce an isotropic return) but
  live in different solver machineries (discretized-µ quadrature vs
  continuous integral transport).
