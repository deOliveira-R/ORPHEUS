---
name: geometry-value-and-mesh-verb
description: Shape of orpheus/geometry at 2026-09-25 (#405 P1) - StructuredGeometry is the geometry value, boundary laws already geometry-layer, kind spelled 13 ways, BC unhashable
metadata:
  type: project
---

Tree at `99ac3e66` (2026-09-25); re-verify before repeating. Full report: `scratch/reference_architecture/p1_explore_geometry.md`.

- `StructuredGeometry` (1-D only) is already a geometry VALUE distinct from a mesh; the geometry-to-mesh verb sits on the wrong type (`Mesh1D.from_geometry`), and moving it onto the geometry adds no module edge (structured_geometry already imports `mesh.BC` at runtime).
- Boundary laws are geometry-layer already (`orpheus/geometry/boundary/`, 7 concrete laws + `BC` tag in `mesh.py`); tags are parsed with FACE context at L2 (`transport/method.py` `resolve_boundary_conditions`).
- `BC.params` is a dict, so `BC` and `StructuredGeometry` are UNHASHABLE (`==` works); `MaterialMesh._law_key` is the existing workaround. Any content digest must handle it.
- `AxisMesh`/`RadialAxisMesh`/`AxisCoord` live in `transport/mesh/axis.py` but import only geometry + numerics: re-homable (#393).
- Geometry kind: 13 closed vocabularies (CoordSystem, SLB/CYL/SPH, AxisCoord, and ~8 derivations-own); production discriminates only on the enum, derivations on strings (118 sites).
- Census traps met: the `_RM`/`_M` import aliases in cp/moc solver defaults evade a last-segment AST filter; derivation generator "constructions" of StructuredGeometry are mostly docstring examples.

**Why:** #405 P1 builds `geometry.mesh(discretization)` on this layer. **How to apply:** start any geometry/mesh carve from this map, re-measuring counts with the scripts in `scratch/reference_architecture/p1geo/`. Related: [[reference-producer-landscape]].
