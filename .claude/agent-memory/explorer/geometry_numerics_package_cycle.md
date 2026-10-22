---
name: geometry-numerics-package-cycle
description: geometry<->numerics cycle is __init__-grain only; the P1 inject-and-run probe (three stages) breaks it; a module move re-keys content digests
metadata:
  type: project
---

`[M]` 2026-10-22, main @ f672a2d5 (layer_dag P1). Re-verify before repeating.

- Module-grain Tarjan over geometry+numerics: 0 cross-package SCCs. The cycle lives in the eager `__init__`s: a cold `orpheus.numerics.scalars` loads 64 orpheus modules (21 geometry).
- Inject-and-run probe: `scratch/layer_dag/p1_imech/` (`apply_moves.py <tree> A|B|C`, `post_patch.py`, `move_realizer.py`, `graph.py`, `cold_import.py`). The stages are A (content and scalars to a base package), B (symmetry, roots_of_unity, manifold and face naming to geometry) and C (the boundary split).
  - With all three applied: 0 package SCCs at top level, 392 of 392 modules import cold, and `import orpheus.geometry` loads 0 numerics modules.
  - Each stage alone leaves {geometry, numerics}.
  - The pairwise file overlaps are 8, 10 and 14, so the stages are not file-disjoint.
- Lazy (PEP 562) `numerics/__init__` is refuted by a run: `geometry.transformation` still loads 16 numerics modules, and the package SCC remains.
- **A pure module move changes content digests**: `content._schema_tag` writes `cls.__module__`. Moving `content.py` reds 5 pinned-digest gates.
- `test_concept_order` cannot see a module that left numerics (no stale-RANK check). The layer gate's TYPE_CHECKING tolerance admits ANY forbidden target from L1/L2, not only L3.
- Resisting edges: the `StructuredGeometry` law arm (only tests feed it, by spy) and `Manifold.__mul__` returning `Product` together with `ambient_dim`'s match on IndexSet/EnergyGroups.

Related: [[coordinate-system-vs-group-seam]], [[boundary-law-method-matrix]], [[content-identity-landscape]].
