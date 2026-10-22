# The layer DAG — every layer untangled and organised until the graph shows the ontology

Living plan, opened 2026-10-22. Parent: `verification_tree_refinement.md` rulings 30 and 31 (the sub-cone measurement as an ontology instrument; layer by layer, a backlog beside it). Predecessor: `numerics_untangling.md`, whose S7 closed the numerics concept order (its ruling 89: 0 violations over 112 module edges).

## Goal, in the domain's terms

Each package is a well-bounded concept set: it defines its objects from the packages below it and never reaches a package above it, so the import graph is the ontology and a misplaced concept is visible as an edge that points the wrong way. Every cycle can be broken; a cycle is a concept in the wrong package, not an essential coupling.

## Rulings

1. **Geometry is declarative, and numerics imports it** (the user, 2026-10-22). "I don't think there is any reason for geometry to import quadrature (which is what probably makes the boundary law cycle). Boundaries at geometry is declarative, and the only relevant boundary at geometry is a deck transformation type boundary, since it affects symmetry. It should be possible to break any cycle and make them well-bounded packages. Numerics and quadratures are probably the ones to import from geometry to be able to get group, support, rigid body applications, etc."
   - Reading: the rigid-motion group and its actions (`geometry/transformation.py`), the symmetry groups and the supports they act on are geometry; numerics (measures, quadratures, bases, frames, spaces, operators) imports them. A boundary in geometry is a DECLARATION; the one boundary concept that belongs to geometry itself is the deck-transformation type (a reflective or periodic face is a group element acting on the domain, so it shapes the symmetry). Boundary LAWS that build operators from quadratures and function spaces are above numerics and leave geometry.
2. **EnergyGrid stays in `data` provisionally** (the user, 2026-10-22). Candidate for later: the `EnergyGrid` class declared in numerics (it is a partition of the energy line whose measure is a `DiscreteMeasure`), and `data` importing it to instantiate the well-known grids (ORNL, SHEM, XMAS, WIMS), which are data. On the backlog; "organizing the graph is more important for now."
3. **Break the cycles that cross layer interfaces first, in parallel worktrees where the files are disjoint** (the user, 2026-10-22).
4. **Architecture E is P1's vehicle** (the user, 2026-10-22): "We can do that, and make E part of the plan (at least the part that breaks the edges that crossed between modules)." E is the boundary-law ontology's ruled design (`boundary_law_ontology.md`, fourth to tenth exchanges; #551): deck transformations (`DeckTransformation`: `MirrorDeck`, `TranslationDeck`, `RotationDeck`, `GlideDeck`) are geometry, declared beside the coordinate system; every other face carries an integer boundary tag, propagated like a material id; the responses (`BoundaryCondition`: `VacuumBoundary`, `AlbedoBoundary(α, SpecularReturn)`, white, inflow; a perfect specular wall is a response, a mirror is a deck) bind at the method (`SNDiscretization(material_mesh, boundary_conditions)`), above numerics, where the method's admission query answers. Its earlier sequencing after #405 (the harness cache, for cadence) is lifted for the part that breaks the cross-package edges; #405 stays on the backlog.

5. **Overnight, the orchestrator decides, with the Fable advisor** (the user, 2026-10-22, before sleeping): "During our re-organization, in pretty much all choices we took it was the one you recommended. I think we will be fine with coding. Just prepare for context compaction, and after compaction the Fable advisor will be available. Then you 2 can take what you consider to be the best decisions overnight." And: "As usual, decisions should be grounded in correctness of the mathematical concepts and machinery, and effort is just a measure of how long it takes for us to get there." So the open questions (P1's primitives home and response-kernel home, a shared Problem type if P3 needs one, the layer-table change) are decided by the orchestrator with the advisor's review, recorded here as rulings marked "orchestrator + advisor, under ruling 5", and reported in the morning for the user to confirm or reverse.

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

- **The P1 measurement** (`[M]` 2026-10-22, the explorer, `scratch/layer_dag/p1_boundary/`, AST over 92 modules, cycle search controlled by the known `traced_memo` pair): **no import cycle crosses geometry and numerics at MODULE grain**, with or without lazy and typing imports; the package cycle exists only through the two packages' eager `__init__` (counting those: one cycle of 32 modules, 19 geometry and 13 numerics; importing `orpheus.numerics.scalars` alone loads 63 orpheus modules, all 17 boundary modules included). Geometry's 73 public items in 25 modules: 8 deck transformation (`transformation` ×5, `BoundaryGeometryMap`, `SelfPairedDeck`, `PairedDeck`), 47 declaration (26 shapes, `BC`, the law-sum algebra, 7 response kernels, 3 inflow sources, the periodic/vacuum/zero-flux/prescribed laws, 3 parameter errors), 4 mixed (`BoundaryTraceLaw`, `AlbedoBoundary`, `ReflectiveBoundary`, `WhiteBoundary`: declarations carrying `assert_*(quadrature)` certification), 14 law/realization (`stamp_boundary_role`, `realize_recursively`, `BoundaryRealizer`, `_BoundBoundaryOperator`, the 4 `assert_specular_pairing_*`, 6 certification errors). All 5 geometry → `numerics.quadrature` edges are `TYPE_CHECKING` annotations on certification methods: no declaration needs a quadrature (ruling 1 confirmed). `numerics/symmetry.py` and `numerics/manifold.py` import only tier 0–1 numerics plus `geometry.transformation` (they can move down; blast radius 168 import statements in 106 file entries, 13 string paths, 3 gates; `manifold` also holds `EnergyGroups`, `ENERGY`, `IndexSet`, `Product`, which are not geometry). The primitives geometry needs (`scalars`, `content`, `face_layout`, `registry`, `roots_of_unity`) are closed and import nothing above them, but sit inside the eagerly-imported `orpheus.numerics`. Under E (ruling 4) the responses, the kernels, the certification and the realization leave geometry together, which dissolves both resisting edges (`parse_boundary_law`'s `isinstance(law, BoundaryTraceLaw)` becomes tag parsing; the reflective law's `_specular` methods leave with the law, which becomes a `MirrorDeck`).
- **P1 — geometry ↔ numerics** (ruling 1 settles the direction; ruling 4 names the vehicle). geometry = the group core (`transformation`), the point groups (`numerics/symmetry.py` moves down), the supports (`numerics/manifold.py`'s G-spaces and quotients move down), the shapes (`chart`, `chord`, `line`, `coord`, `extent`), and the boundary DECLARATIONS with the deck-transformation type; the primitives (`scalars`, `content`, `face_layout`, …) form a base package under geometry; the boundary laws move above numerics. Owed before the ruling on the split: the boundary-package measurement (dispatched 2026-10-22).
- **P2 — the energy arm** (ruling 2). `EnergyAxis.from_grid(grid)` becomes `EnergyAxis.of_edges(edges)`; `from_materials` (the reachable-materials rule, which reads `Mixture.eg`) moves to `data`. Numerics then names nothing in data.
- **P3 — transport naming method types**, and **P4 — the `sood_registry` retirement**: as above, each its own worktree.

## ⏸ Next

S7 merged @ `e80ec738`. Next: the sub-cone measurement on numerics (`verification_tree_refinement.md` ruling 30); in parallel worktrees from `main`: P2 (the energy arm), P3 (transport naming method types), P4 (`sood_registry`); P1 as architecture E's cross-package part, designed first (the primitives base: its own package or a lazy `numerics/__init__`, an open question for the user; the response kernels' home under E).

## ⏸ COMPACTION POINT — 2026-10-22

**State.**
- S7 is merged and pushed (`e80ec738`, #631 closed; CI `gates` green). The landing stamps and this plan's rulings are pushed at `497be3cf`; its CI run 39742912609 was in progress at this point.
- This compaction note's own commit sits on `main` LOCALLY, unpushed: a plan commit never heads a push, so the next code or docs push carries it.
- The S7 worktree `.claude/worktrees/s7` (branch `refactor/s7-axis-identity`) is fully merged; remove it with `git worktree remove .claude/worktrees/s7` and delete the branch.

**Next, in order.**
1. **The sub-cone measurement on `numerics`.** Sources: `verification_tree_refinement.md`, ruling 30 and "Sub-cones: the measurement"; S8's cone store comes first, with its positive controls. Read the result as an ontology instrument.
2. **In parallel worktrees cut from `main`:**
   - P2, the energy arm: `EnergyAxis.of_edges`, with `from_materials` moved to `data`;
   - P3, transport naming the sn and diffusion types;
   - P4, the `sood_registry` retirement.
3. **P1, architecture E's cross-package part (ruling 4).** Two user questions are still open:
   - is the primitives base its own package, or does `numerics/__init__` become lazy?
   - where do the response kernels live under E?
   - Design comes before code. E's design record is `boundary_law_ontology.md`, the fourth to tenth exchanges, plus `scratch/boundary_ontology/orbifold_architecture_round2.md`.

**Read before acting:**
- this plan's rulings and measured state;
- `verification_tree_refinement.md`, rulings 30 and 31 and "The DAG backlog";
- for P1, the census probes in `scratch/layer_dag/p1_boundary/` and `scratch/layer_dag/layer_census.py` and `iface_names.py` (run: `python <script> <repo root>` from outside the repo).

## ⏸ COMPACTION POINT — 2026-10-22, night (round 1 running)

**State.**
- `main` = `f672a2d5` locally (one plan commit ahead of `origin/main`; CI green on `497be3cf`). This compaction note is a further plan commit on `main`, unpushed. The S7 worktree is removed.
- **Workflow `wf_321d7cf8-ff5` (background, task `whtwtu6k0`) is running round 1.** Script: `~/.claude/projects/-Users-rodrigo-git-nuclear-ORPHEUS/54d619d4-06ec-4f9f-b3c1-bc030bf79384/workflows/scripts/layer-dag-round-1-wf_321d7cf8-ff5.js`; transcript dir `.../subagents/workflows/wf_321d7cf8-ff5` (its `journal.jsonl` holds each agent's return). It runs:
  - P2, P3, P4: design (explorer) → adversarial review (cross-domain-attacker) → implementation in an isolated worktree on branches `refactor/layer-dag-p2-energy-arm`, `refactor/layer-dag-p3-transport-method-agnostic`, `refactor/layer-dag-p4-sood-registry` → qa. A unit whose design needs a user ruling stops after design.
  - P1: three proposals (lenses ontology, import-mechanics, consumer) → cross-domain attack + elegance review → one synthesis memo. No code.
  - Intra-layer diagnosis: mesh {axis, face_laws, structured}, data {macro_xs, materials}, transport {fields, source_sinks}. Read-only.
  - Sub-cone re-measurement on numerics (test-architect) into `scratch/verification_tree/subcones_s7/` against the MAIN checkout.

**Overnight protocol (the advisor's review, adopted; ruling 5 gives the decisions to the orchestrator and the advisor).**
1. Do not touch `orpheus/` in the main checkout (no merge, no checkout of code) until the sub-cone capture's log shows its `END` line; the capture measures the main checkout. Integrate in a separate worktree meanwhile.
2. Integration: stack P4 first, then P2 rebased onto it, then P3. All three edit `tests/gates/test_layer_imports.py`: combine the edits and re-run the gate after each rebase. Check each branch exists with commits (an isolated worktree with no commits is auto-removed).
3. Bar before any fast-forward of `main`, at the stacked tip: the full suite `-m "not slow"`, serial, in the background (timeout 7200000); pyright with the planted-probe protocol; Sphinx `-W` after deleting `docs/_generated/vv_audit.json`; `dead_references` with `use_workspace` set to the integration worktree and a planted positive control (mandatory for P2's rename). After the push, watch CI; if red and not fixable in one cycle, revert the range.
4. Review checks the qa briefs may miss:
   - P2's first red must be a TYPE_CHECKING import of `orpheus.data` from numerics: the gate's tolerance (`is_tc and src_pkg in L1|L2`) skips EVERY forbidden typing edge from L1/L2, so the tolerance must narrow, or the gate has no teeth for exactly the edge P2 removes.
   - Read whether `EnergyGrid` is a `DiscreteMeasure`; if so, the numerics-side spelling is the S7 constructor `of_measure`, not an `of_edges(ndarray)` that strips the measure to a bare array.
   - P3: the shared Problem concept, if needed, is decided under ruling 5 with the advisor, by what the sites actually read.
5. P1: decide Q1–Q4 under ruling 5 with the advisor, from the synthesis memo; save the memo under `scratch/layer_dag/p1_design/` and into this plan. Evidence for Q1: an inject-and-run probe in a scratch worktree, a lazy `numerics/__init__`, counting the modules `import orpheus.numerics.scalars` loads (63 today). Then P1's work packages may start in worktrees.
6. Intra-layer: the transport cycle waits for P3 (shared files); mesh and data follow in their own worktrees after the cross-layer units land (ruling 3).
7. Morning: one report (an artifact or a file sent to the user) and a new compaction point here: what merged (hashes, CI), what sits on branches and why, every decision taken under ruling 5 with its argument, so the user can confirm or reverse.
