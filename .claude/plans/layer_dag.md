# The layer DAG — every layer untangled and organised until the graph shows the ontology

Living plan, opened 2026-10-22. Parent: `verification_tree_refinement.md` rulings 30 and 31 (the sub-cone measurement as an ontology instrument; layer by layer, a backlog beside it). Predecessor: `numerics_untangling.md`, whose S7 closed the numerics concept order (its ruling 89: 0 violations over 112 module edges).

## Goal, in the domain's terms

Each package is a well-bounded concept set: it defines its objects from the packages below it and never reaches a package above it, so the import graph is the ontology and a misplaced concept is visible as an edge that points the wrong way. Every cycle can be broken; a cycle is a concept in the wrong package, not an essential coupling.

## Rulings

1. **Geometry is declarative, and numerics imports it** (the user, 2026-10-22). "I don't think there is any reason for geometry to import quadrature (which is what probably makes the boundary law cycle). Boundaries at geometry is declarative, and the only relevant boundary at geometry is a deck transformation type boundary, since it affects symmetry. It should be possible to break any cycle and make them well-bounded packages. Numerics and quadratures are probably the ones to import from geometry to be able to get group, support, rigid body applications, etc."
   - Reading: the rigid-motion group and its actions (`geometry/transformation.py`), the symmetry groups and the supports they act on are geometry; numerics (measures, quadratures, bases, frames, spaces, operators) imports them. A boundary in geometry is a DECLARATION; the one boundary concept that belongs to geometry itself is the deck-transformation type (a reflective or periodic face is a group element acting on the domain, so it shapes the symmetry). Boundary LAWS that build operators from quadratures and function spaces are above numerics and leave geometry.
2. **EnergyGrid stays in `data` provisionally** (the user, 2026-10-22). Candidate for later: the `EnergyGrid` class declared in numerics (it is a partition of the energy line whose measure is a `DiscreteMeasure`), and `data` importing it to instantiate the well-known grids (ORNL, SHEM, XMAS, WIMS), which are data. On the backlog; "organizing the graph is more important for now."
3. **Break the cycles that cross layer interfaces first, in parallel worktrees where the files are disjoint** (the user, 2026-10-22).

## Measured state (`[M]` 2026-10-22, the S7 worktree, `scratchpad/layer_census.py` and `iface_names.py`; AST over `orpheus/`, top-level + lazy + `TYPE_CHECKING` imports)

Cross-layer cycles (package grain):

| cycle | edges closing it | names |
|---|---|---|
| geometry ↔ numerics | runtime, top-level | numerics → geometry: 4 names, all `geometry.transformation` (`RigidMotion`, `Permutation`, `close_group`, `permutation_between`), read by `symmetry`, `manifold`, `measure`, `quadrature`; geometry → numerics: 23 names, nearly all from `geometry/boundary/` (operator algebra, `Quadrature`, `FunctionSpace`, `TANGENTIAL_EPS`, `RegistryMixin`), plus primitives (`scalars`, `content`, `face_layout`), `SubgroupOfO3` (`chart`, `_specular`) and `Quotient` (`_specular`) |
| data ↔ geometry ↔ numerics | typing only | `numerics/axis.py:146-147` names `EnergyGrid`, `Mixture` (`EnergyAxis.from_grid`, `from_materials`; callers: `homogeneous/solver.py:208`, `transport/mesh/material_mesh.py:565`) |
| transport ↔ sn, diffusion | typing only | 7 sites (`transport/fields/_bases.py:119,121`, `fields/angular_flux.py:53`, three `source_sinks` modules, `radial_characteristic_field.py:83`); tolerated by `test_layer_imports.py`'s TYPE_CHECKING rule |
| cp ↔ derivations | lazy | `sood_registry`, tagged `RETIRE_IN_P3_FOLLOWUP` in the layer gate |

Intra-layer cycles (module grain, subpackages as nodes): geometry none; mesh {`axis`, `face_laws`, `structured`} (lazy/typing); data {`macro_xs`, `materials`} (lazy/typing); transport {`fields`, `source_sinks`} at top level, plus `mesh` with typing.

## Proposals (the orchestrator's, `[R]`, open to the user)

- **P1 — geometry ↔ numerics** (ruling 1 settles the direction). geometry = the group core (`transformation`), the point groups (`numerics/symmetry.py` moves down), the supports (`numerics/manifold.py`'s G-spaces and quotients move down), the shapes (`chart`, `chord`, `line`, `coord`, `extent`), and the boundary DECLARATIONS with the deck-transformation type; the primitives (`scalars`, `content`, `face_layout`, …) form a base package under geometry; the boundary laws move above numerics. Owed before the ruling on the split: the boundary-package measurement (dispatched 2026-10-22).
- **P2 — the energy arm** (ruling 2). `EnergyAxis.from_grid(grid)` becomes `EnergyAxis.of_edges(edges)`; `from_materials` (the reachable-materials rule, which reads `Mixture.eg`) moves to `data`. Numerics then names nothing in data.
- **P3 — transport naming method types**, and **P4 — the `sood_registry` retirement**: as above, each its own worktree.

## ⏸ Next

S7 merges first (the numerics knot). Then the sub-cone measurement on numerics (`verification_tree_refinement.md` ruling 30) and, in parallel worktrees from the merged `main`, P2, P3, P4; P1 after the boundary measurement and the user's ruling on the split.
