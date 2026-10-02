# Layer 1 ends in systems that carry their equation, and layer 2 poses physics-free questions on them (posing unit 6)

Status: SCHEDULED, not opened (2026-09-29). Issue #529; charter issue #522. Charter: `.claude/plans/posing_sequence.md`, "The work, ordered by dependence"; ruled content: "The ontology as it stands", layers 1 and 2 and the plan/binding paragraph of layer 3. Depends on units 2 (the frame faces), 3 (the terms from the grid's cells), 4 (the direct-sum space), and 5 for the Strategy side. This is the layered machinery; it opens only when those preconditions are merged, and its first act is a phase plan with compaction points (plan-authoring §6).

## Goal, in the domain's terms

A method's discretisation and its materials give a SYSTEM: an equation family over its own parameter space, with every coordinate declared (its affinity, chart, admissible range, continuum edge, sign), its time term, its cone and bulk block, its ports and its physical point. Systems are manipulated by an operator-level algebra (couple, eliminate, augment, restrict, pull back, complexify). A question (eigen, fixed source, evolution) names a coordinate of a system as an opaque key and a base point, and knows no physics; its answer is read from one resolvent. Today the hubs (`SNProblem`, `DiffusionMesh`, `HomogeneousProblem`) weld all three layers.

## Phase 0: the layer table (the first deliverable)

Member by member, the target layer of every public member of `SNProblem` (33 public members `[M]` 2026-09-28, 886 name sites), `DiffusionMesh`, `HomogeneousProblem`, `EigenPosing` (59 lines, 14 files `[M]` 2026-09-29), `SourcePosing` (47, 14), `posing`, `EigenOutcome`/`SourceOutcome`. The elegance attack's partial split (`scratch/posing_sequence/clean_attack2/elegance/memo.md`, C6) is the seed. Name maps: `EigenPosing → EigenProblem`, `SourcePosing → FixedSourceProblem`, `<Kind>Outcome`. The word "problem" changes meaning (hub → question) at every site: the rename order is part of the table.

## Scope (ruled 2026-09-27/29; the section is the authority)

1. **`System`**: spaces with measures and metric; `E(p)` over its own parameter space; `T_d(p)`; declared coordinates (affinity per chart and reduction, chart with its zero as the algebraic origin, admissible range, continuum edge, sign of `⟨1, T_d ·⟩`); time with `T_time`; cone and bulk block; ports; `physical_point` (materials as its coordinates). The neutron system's terms are the grid's role-marginals (unit 3); the precursor system seeded.
2. **`Discretization`**: the product of per-axis frames, both faces declared (analysis; reconstruction by basis or by characteristic solve), the content digest keyed on the region partition and the frames; `SNDiscretization`, `DiffusionDiscretization`, `InfiniteMediumDiscretization` (`orpheus/homogeneous/` → `orpheus/infinite_medium/`).
3. **The operator-level system algebra** with derived declarations (an elimination's rational affinity; a complexification's absent cone).
4. **Layer 2**: `EigenProblem`/`FixedSourceProblem`/`EvolutionProblem` over their values with one `point` slot, in `numerics`, over the question VALUES that #405 P1 step 7 lands first (`Eigen`, `FixedSource`, `Response`, the point and the `Fundamental`/`Nearest` mode values; amended 2026-10-02, see `reference_p1_spec.md` §1.7 "The step-7 rulings", so this unit binds and resolves them and does not re-mint them); the parameter kinds `CellCoefficient | GeometryExtent | NuclideDensity` as the SYSTEM's coordinate declarations (P1 step 8 declares them on the reference specification first; this unit moves or reuses that declaration, one definition); the mode law (`Fundamental` from the admissible-range end on the bulk cone; `Nearest`; `Enclosed` refusing beyond the continuum edge, minted here); the `Mode` law with its three gates; `SingularSource` with the Fredholm refusal (#460); the pencil (`shifted`, `jet`, `pseudo_resolvent`; the zero-scale law from unit 1; the rational delayed-neutron case, #463); `Stepped` (stencil only) and `Exact`; spectral maps as a monoid; derived questions by content identity; `rate(cells, w)` if unit 1 did not take it.
5. **Layer 3's plan/binding split** (with unit 5): the plan reads declarations and patterns; the binding is everything at a point.
6. Related issues: #463, #460, #465, #467, #469, #472, #484, #462.

## Gates (the charter's layer-1 and layer-2 lists)

Affinity against the second-difference test; the mesh refines the region partition; the trace commuting condition on the reproducing space; derived declarations; the eliminated adjoint; `eliminate ∘ augment = id`; `Fundamental` along c; the `Mode` law's gates; Hellmann–Feynman; the two-Strategy gauge gate; the ω → 0 point-kinetics gate; answer-invariance under the base point; the `shifted`, zero-scale and Möbius laws; the metric pairing gates; seed independence; parallel against sequential `Stepped`; a bounded-body α check (Carlvik's sphere against 1-D spherical SN).

## Done-when (hypothesis)

The hubs are split along the three layers per the phase-0 table; the question types live in `numerics` and name no reaction cell; the gates green with their first reds. To be restated per phase in the phase plan.

## Sizing

8–12 sessions `[R]`.
