---
name: geometry-numerics-package-cycle
description: geometry<->numerics cycle is package-grain only (eager __init__s); module grain acyclic; geometry->quadrature edges all TYPE_CHECKING on assert_* certification
metadata:
  type: project
---

`[M]` 2026-10-22, s7 worktree @ 96ec89fd (layer_dag P1 measurement). Re-verify before repeating.

- Module-grain Tarjan over geometry+numerics (92 modules): 0 cross-package SCCs, even with lazy+typing edges. The cycle exists only through package `__init__` eagerness: importing `orpheus.numerics.scalars` loads 63 orpheus modules incl. all 17 `geometry.boundary` modules. A "primitive base under geometry" therefore must sit outside the `orpheus.numerics` package (or its `__init__` go lazy).
- All 5 geometry->`numerics.quadrature` edges are TYPE_CHECKING annotations on `assert_*(quadrature)` certification (`_base`, reflective, albedo, white, `_specular`). No deck/response/tag/source reads a quadrature.
- Runtime edges above tier 2: only `_specular` (TANGENTIAL_EPS, tier 16), `_realizer` (BlockRole), `_bound_compat` (LinearOperator subclass).
- Resisting edge: `structured_geometry.parse_boundary_law` isinstance-doors on `BoundaryTraceLaw`, so the descriptor must stay in geometry; split off the certification methods instead.
- manifold.py holds non-spatial manifolds (EnergyGroups, IndexSet, Product) beside the G-spaces.

Probes: `scratch/layer_dag/p1_boundary/{census,uses,scc,probe_load,q3b}.py`. Related: [[coordinate-system-vs-group-seam]], [[boundary-law-method-matrix]].
