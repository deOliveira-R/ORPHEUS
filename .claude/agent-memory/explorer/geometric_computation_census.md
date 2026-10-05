---
name: geometric-computation-census
description: What Euclidean geometry orpheus computes (points/directions, surfaces, intersections, normals, reflection, mesh metrics, line measures), where, and which spellings are twins.
metadata:
  type: reference
---

Surveyed 2026-10-05 at `a336bde4` for the characteristic-reference plan's widened question (the user's "geometric kernel" ruling). Full report: `scratch/characteristic_architecture/geometry_census.md`; scripts and probes are in `geometry_census/` beside it. Re-verify the counts before repeating them.

- The algebra exists, the values do not. `RigidMotion` (E(d), Householder `reflection`, `on_points` vs `on_directions`) is the only typed point/direction split. `on_directions` has 0 production callers. There is no Point, Direction, Line, Ray or Surface class; coordinates are bare floats and arrays. Directions are typed only as the measure's support (`Sphere` vs `RealSpace`).
- Every ray–surface intersection is one quadratic on concentric circles or spheres plus axis-aligned planes. Outside origins: 13 `sqrt(sq−sq)` sites and 14 geometric discriminants, nearly all B1. The only production tracer is MoC (box + circles). MC uses delta tracking, with no distance-to-surface. CSG exists only as promises in docstrings.
- Mesh metrics are degenerate: face normal = (axis, sign), "centroid" = coordinate midpoint (the curvilinear volume centroid is computed nowhere), centroid-to-face = h/2, non-orthogonality has 0 hits. Unstructured meshes and CSG meshing are ruled external (2026-09-24); the seam is #322, #335, #539.
- The measure has one definition (`CoordSystem.measure`) and about 14 inline re-spellings (Peierls and flat-source CP kind-tagged, MoC, fuel/TH/kinetics). Point location has 8 spellings with split boundary bias: Peierls `which_annulus` and MCMesh are outer-biased, the rest inner-biased.
- Cauchy 4V/S is stated in the CP docs and computed nowhere. TH's `d_hyd = 4V/A` is a coincidence of form, not a twin.

**How to apply:** start any geometry-kernel, ray-tracing, mesh-metric or CSG question here. Separate the transformation half (exists) from the value half (absent). Separate the ray half (consumers today) from the FV-metric half (no consumer until unstructured meshes). Related: [[ray-characteristic-machinery-map]], [[spatial-transform-category-durable]], [[coordinate-system-group-seam]], [[geometry-value-and-mesh-verb]].
