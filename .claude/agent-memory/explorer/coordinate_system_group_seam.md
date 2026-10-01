---
name: coordinate-system-group-seam
description: CoordSystem knows no group; only the point stabiliser is typed (string-keyed, test-only); decks are RigidMotion elements; the r=0 pole is a mirror deck; no spatial group, fold or fold test.
metadata:
  type: project
---

Surveyed 2026-09-30 at `5a5ffa7b` for the boundary-law ontology plan (survey B); full report `scratch/boundary_ontology/survey_symmetry.md`, probe `_survey_symmetry_probe.py` beside it. Re-verify counts before repeating them.

- The group machinery (`numerics/symmetry.py`, `invariance.py`, `manifold.py`) is ANGULAR: subgroups of O(3), orbit spaces of the sphere only (catalogue keys all `(Sphere, ·)`). `RigidMotion` reaches E(d), but `close_group` refuses every infinite group and `SubgroupOfO3` fixes the origin, so no spatial group with a translation is spellable (a slab's two mirrors, a wrap).
- `CoordSystem`/`AxisCoord` carry no group. The typed object is the POINT STABILISER (`AngularSymmetry.spent`) in `GEOMETRY_ANGULAR_SYMMETRY`, keyed by strings disjoint from `CoordSystem`, read only by `select_quadrature` (all callers tests). The group's consequences (measure, law-free centre, pole, α redistribution, cylinder σ_y fold) are restated by hand; `coord.py`'s "single point" docstring is false (most discrimination sites lie outside it).
- Deck types hold a `RigidMotion` (element), never a `SubgroupOfO3`; the mirror is built twice (`_mirror_motion` via `face_layout.AXIS_NAMES`, `symmetry._reflections` via `AXIS_INDEX`), agreeing, unpinned against each other. A deck meets the quadrature per motion (`ordinate_permutation` → `_orbit_closure`), never per group; the ledger's `owed` has no production reader.
- The r=0 regularity is realized as `SelfPairedDeck.mirror("x")` (`_ensure_pole_mirror`); the geometry says the centre "carries no law". `boundary_points` and `RadialAxisMesh.endpoints` disagree on hollow bodies (#511 guard).
- No spatial fold is built from a full domain, and no production method has a full-vs-folded test (the only one is the Peierls slab derivation); the angular σ_y fold is mandatory for the SN cylinder.

**Why:** the plan's O4 infers a symmetry group from colourings and derives deck faces by folding; this map says what exists to build it on. **How to apply:** start any symmetry/deck/coordinate-system question here; the durable seams are the element/group split and the stabiliser-only typing. Related: [[spatial-transform-category-durable]], [[angular-layer-hidden-transformations]], [[geometry-value-and-mesh-verb]].
