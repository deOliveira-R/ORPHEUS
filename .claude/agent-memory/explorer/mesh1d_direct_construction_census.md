---
name: mesh1d-direct-construction-census
description: #405 P1 (2026-09-25) census of 450 direct Mesh1D(...) sites - 410 are geometry+rule, D rebuilds bit-exact as region-per-cell, the only production mint is the legacy axis adapter
metadata:
  type: project
---

Tree `99ac3e66` (2026-09-25); re-verify before repeating. Report: `scratch/reference_architecture/p1_mesh1d_census.md`; probe dir `p1_mesh1d_census/` (census.py, controls_src.py, rebuild_probe.py, adjudicate.py).

- 450 direct `Mesh1D(` AST sites (tests 435, production 3, derivations 10, examples 2) plus 6 rst sites. 410 of the 450 are a geometry plus a rule (uniform 345, single 22, region-per-cell 35, hand CellsByCount 2, named rule 6). The 40 without a rule are irregular literals (27), derived (8) and constructor negatives (5). No random site.
- The only production direct mint is `legacy_mesh_from_axes` (axis→mesh adapter, called only by `SNProblem.from_axes`). Production axes come from `axes_from(mesh)`, so it is a round trip.
- Rebuild through `from_geometry`: irregular D sites are bit-exact as ONE REGION PER CELL (25/26; 1 ULP off from thickness accumulation). A sites are bit-exact per mat run ONLY with `method="uniform"`: the default `"equal-volume"` differs on 58 curvilinear sites. StructuredGeometry admits adjacent same-material regions.
- `precomputed_volumes` is set at 0 of 450 direct sites (only `from_geometry` sets it).
- Census trap: an AST value-evaluator that resolves `x = []` and ignores a later `x.extend(...)` silently reads an empty array; guard mutated names (augassign, `.append/.extend`, subscript store).

**Why:** #405 step 2 decides whether a `from_cells(edges, ...)` spelling has a permanent consumer. **How to apply:** start any Mesh1D-constructor retirement from `final.tsv`. A migration that uses the default RegionMesh method will move curvilinear pins. Related: [[geometry-value-and-mesh-verb]].
