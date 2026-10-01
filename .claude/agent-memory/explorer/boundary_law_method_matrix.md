---
name: boundary-law-method-matrix
description: Which boundary laws each method (SN, diffusion, CP, MoC, MC) admits, realizes, or mis-realizes at HEAD 5a5ffa7b; where the admission gates really sit
metadata:
  type: project
---

Durable shape of the boundary-law consumers, measured 2026-09-30 at `5a5ffa7b` (census A of the boundary-law ontology plan; full table in `scratch/boundary_ontology/census_laws.md`, probes beside it). Re-verify any row before repeating it.

- **Admission is two-tiered.** A `BC` tag is gated by the method's `BOUNDARY_OPERATOR_REGISTRY` (SN {vacuum, reflective}; diffusion {vacuum, reflective, albedo, zero_flux}) in `_law_from_tag`; a TYPED law declared on the geometry/mesh skips that table (arm since `985497b5`, 2026-08-05) and meets only the realizer's isinstance ladder, then (1-D radial SN) the corner predicate `_has_ruled_corner_action` at SOLVE time, then the k estimator. A census of "what SN admits" that reads only the registry undercounts.
- **Admitted but mis-realized** is a real status: SN's k leakage predicate is `response_kernel.is_zero`, so any 0 < α < 1 face gives k = k∞ (operator itself is right; fixed-source flux is right). CP's `compute_keff` assumes zero leakage, so CP vacuum k is wrong (vacuum above white on a heterogeneous slab). CP's slab kernel has a built-in mirror at x = 0 whatever the declared left law.
- **CP, MoC, MC do not use the law pipeline.** They read `mesh.outer_law.kind` against their own string `BC_REGISTRY` (CP {white, vacuum}; MoC {reflective}; MC {periodic}); MoC realizes the law on the equal-area SQUARE (nearest-point track linking, never exact); MC applies `% pitch` for any geometry object.
- **No law is realized by all five methods.** Discriminating fixtures: a heterogeneous fuel|moderator slab (left law), a leaking 2-cm slab (partial albedo), monotonicity k(white) ≥ k(vacuum). A homogeneous body is blind to every closed law and to CP's leakage bug.

**Why:** the user's thesis (2026-09-30) that deck laws are honoured by every method is a tree question this census answered "no" (plumbing refusals mostly). **How to apply:** for any boundary-law question, probe BOTH declaration arms and at least one solve on a discriminating fixture; never trust a registry dict alone. Related: [[method-hub-commonality]], [[spatial-transform-category-durable]].
