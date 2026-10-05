---
name: geometric-kernel-seed-pencil-crofton-frames
description: 2026-10-05 W5 attack on the geometric-kernel seed (Point/Direction/Line, concentric rho, crossing law, measure on lines) — quadric pencil of G_c's invariant, lines = phase space / streaming flow (Crofton), Plücker Line vs Ray, germ rule
metadata:
  type: project
---

Memo: scratch/characteristic_architecture/w5_cross_domain.md; probe w5/probe_cross_domain.py (seed 20261005).
Verdicts:
- rho is a CHART of R^3/G_c; the crossing primitive is the invariant-ring generator pi = X^T Q_0 X (deg 2 sphere/cyl,
  deg 1 slab) → pencil Q_0 − pi_k E; plane = quadric with zero quadratic part ⇒ ONE formula 3/3. rho^2=r_k^2 with a
  signed slab rho gives 3 spurious of 6 roots [M]. Placement = congruence H^-T Q H^-1 (4e-12; mutation 522/600) [M].
- Seed's invariance gates unspellable: StructuredGeometry 0/4 placement fields, CoordSystem 0/5 group verbs.
- Lines = phase space / streaming flow; dL = dA_perp dOmega; pushforwards: sphere b db, cyl sin^2θ dθ db, slab |μ|dμ;
  Cauchy discriminates [M] (planar db on sphere = πR/2 not 4R/3). CP E3/Ki3/y-quad = Blaschke–Petkantschin + pushforward.
  Sphere's "plane of the ray" is a SECTION, cylinder's a PROJECTION (1/sinθ = arc-length Jacobian).
- Line = Plücker (quotient, D19); Ray = phase-space point; b = |m| = momentum-map Casimir (billiard conservation = Noether).
- Germ rule: Q-b = sign of first nonzero coefficient of g−pi_k; all-zero on slab Ω·n=0 and cyl Ω∥axis ⇒ "grazing
  outward" derived on 1/3 charts only. Do NOT implement Q-b as "first segment" (sliver-root rounding).
- Direction = point of existing numerics.manifold.Sphere; mint nothing.
**Why:** user ruled the chord is the seed of a geometric kernel (2026-10-05).
**How to apply:** open before any ray/line/chord/CSG/FV-metric/point-location brief. Related:
[[coordinate-group-weyl-orbifold-frames]], [[pullback-pair-frobenius-perron-frames]].
