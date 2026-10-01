---
name: coordinate-group-weyl-orbifold-frames
description: 2026-10-01 attack on the boundary-ontology orbifold prototypes — CoordSystem = G_c (fibres of c), decks live in the Weyl group N(G_c)/G_c, spent = H ∩ orbit_stabiliser(H⁰), Bieberbach (Λ,P,v) for lattices
metadata:
  type: project
---

Verdicts (memo: scratch/boundary_ontology/attack_math.md; probe attack_math_probe.py):
- G_c := {g : c∘g = c} from the coordinate map alone reproduces GEOMETRY_ANGULAR_SYMMETRY's spent·unspent 4/4 [M].
  The rule "spent = identity component" is WRONG on disconnected stabilisers (slab/sphere O2_x vs SO2_x) and
  non-canonical (two complements pass). Correct: spent = H ∩ orbit_stabiliser(H⁰).
- Decks = discrete Γ̄ ⊂ W(G_c) = N(G_c)/G_c. W(1-D cylinder) = W(sphere) = 1 ⇒ no deck on a curved 1-D face (theorem).
  (b) "couples frames" is a deck-vs-QUADRATURE gate (point group P ⊄ quadrature group), not deck-vs-chart.
  Non-normalising deck ⇒ its G_c-Reynolds average, a RESPONSE (white = its isotropic projection).
- owed is per-PROBLEM, not per-chart (2/4 rows unreachable from decks).
- Regularity on any stratum: ψ(x,·) is H_x-invariant (mirror = reflective BC; centre = ∂_μψ=0; SN uses only σ_x).
- PairedDeck "free" docstring false for rotations; half-turn self-pairing (p2/half-core) unspellable.

**Why:** the user's orbifold ontology (plan boundary_law_ontology.md, 7th exchange) awaits shape rulings.
**How to apply:** open before any CoordSystem / deck / symmetry-group / lattice brief; re-ground file:line first.
Related: [[quadrature-symmetry-quotient-frames]], [[pullback-pair-frobenius-perron-frames]].
