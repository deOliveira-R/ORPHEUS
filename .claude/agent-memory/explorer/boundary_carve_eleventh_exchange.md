---
name: boundary-carve-eleventh-exchange
description: Boundary carve census (2026-10-22, tip b4a497b7) — six string tag tables, law-key readers, specification keys boundaries via geometry, _factors mixed
metadata:
  type: project
---

Durable shape measured 2026-10-22 for layer DAG P1 / architecture E (full detail: `scratch/layer_dag/next_session/item*.md`, untracked):

- **There are six string tag tables keyed by `BC.kind`**, not two: SN and diffusion each have a `BOUNDARY_OPERATOR_REGISTRY`; `CPMesh`, `MCMesh` and `MOCMesh` each have a `BC_REGISTRY` (cp {vacuum, white}, mc {periodic}, moc {reflective}); and `walls.py` has `TAG_REGISTRY`. The eleventh exchange's ruling 5 names only the first two.
- **The `RegistryMixin` roots beyond boundary** are `AngularClosureBase` and `DiscretizationSchemeBase`. Neither has a production string selection; the one production `.registry` read is `loss_kernel_gauge`'s enumeration.
- **Retiring `key=` removes `BoundaryTraceLaw.kind`.** Its readers are the `__eq__` string arms (#553: vacuum, reflective, `_BoundBoundaryOperator`) and `BoundaryError.law`, which tests pin 12 times.
- **`GeometrySpecification` carries its boundary values INSIDE `StructuredGeometry`.** Under arm T (tags), the specification needs either a response map (a forbidden spec → transport edge) or the reference's closed values. This is an open design ruling.
- **References already own closed boundary values.** `Walls(Wall(specular, diffuse, partner))` in characteristic, and the albedo tuple in `reference_body`. Only `_wall_of`/`specular_albedo` translate a law.
- **`geometry/boundary/_factors.py` is mixed** (decks and kernels), so a carve splits it.
- **The real file collision between WP-B and the carve is `_factors.py` and `_specular.py`.** Every other shared file is a one-line import.

**Why:** the next session plans the WP-B, carve and arm T worktrees from this.
**How to apply:** re-verify against the tree before repeating any of it (the campaign moves these files). See [[boundary-law-method-matrix]] and [[reflective-albedo-retirement-blast]].
