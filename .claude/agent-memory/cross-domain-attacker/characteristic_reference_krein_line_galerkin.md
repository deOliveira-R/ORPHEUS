---
name: characteristic-reference-krein-line-galerkin
description: P1 W5 attack on the characteristic (Variant-alpha successor) reference: closure = Krein boundary resolvent on each line's boundary space; Galerkin over lines; tangency substitution
metadata:
  type: project
---
2026-10-05 W5 attack on D1-D5 (memo `scratch/characteristic_architecture/p1_w5_attacker.md`).

- The "closure rank" is the DIMENSION of a line's boundary space = wall endpoints of the
  accessible interval MODULO the law's deck map G; specular D2 = (I - T_L)^{-1}, the
  specular case of Krein's P = P_0 + E(I-T)^{-1}X. Periodic slab: 2 walls, rank 1, no
  reversal. White (Lambertian R): rank-W block per group. Read R o G from the boundary
  registry (`_factors.py`), never re-derive "reflection reverses".
- [M] angular split at tangencies alone is n^-3; + cosine endpoint map spectral (3e-15 at 32).
- [M] spline collocation P has negative entries (up to 43/1024): gate Perron selection.
- Galerkin over the measure on lines (CP's assembly) = symmetric by reciprocity, O(N_b) chords
  on the sphere vs O(N_r N_mu) rays. Shares derivation with production CP (X4 for CP gates only).

**Why:** second sighting of "line-diagonal law -> geometric series; non-diagonal -> finite-rank
boundary block" would be a Part A row (Krein formula, trigger: Volterra + finite boundary space).
**How to apply:** open before any reference/closure/boundary-law-in-a-reference brief.
