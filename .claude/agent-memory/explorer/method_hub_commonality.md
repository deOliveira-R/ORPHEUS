---
name: method-hub-commonality
description: Method hubs (SNProblem, DiffusionMesh, HomogeneousProblem, CP/MoC/MC meshes) share nothing beyond MaterialMesh except a 2-of-6 boundary-realisation quartet; census recipe.
metadata:
  type: project
---

Measured 2026-09-28 at HEAD 4bf56535 (memo: `scratch/posing_sequence/open_items/transport_method/memo.md`).

- Hub commonality ends at MaterialMesh data. Beyond it no member is carried by >2 of 6 hubs;
  SN∩DIFF∩HOMO = {ng}. HomogeneousProblem, CPMesh, MOCMesh, MCMesh do NOT subclass MaterialMesh.
- SN∩DIFF extra: TransportMethod's quartet + `full_field_space` (bulk ⊕ trace) + `from_material_mesh`
  (its signature difference IS the angle/scheme frames). Trace is named twice (`angular_trace`/`scalar_trace`).
- SN∩HOMO extra: `pencil`, `eigen_posing` (same type); `fission` is a HOMONYM (operator vs datum).
  Diffusion keeps its pencil halves on DiffusionSolver.
- The generic seams already sit BELOW the hub: `BoundaryRealizer` Protocol (`geometry/boundary/_realizer.py`)
  and the posing/outcome layer. `resolve_boundary_conditions`' only callers are the two conformers' `_init_core`.
- Census recipe: runtime MRO `vars()` + AST `self.X =`/dataclass-field pass; control = MaterialMesh ⊆ SN, DIFF.

**Why:** retiring TransportMethod (posing-sequence open item). **How to apply:** re-verify before reuse; the
Nexus graph at that date was built from another branch, so grep/AST were primary.
