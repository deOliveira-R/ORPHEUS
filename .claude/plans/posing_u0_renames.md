# Every loose holder carries the word its structure has (posing unit 0)

Status: SCHEDULED, not opened (2026-09-29). Issue #523; charter issue #522. The charter is `.claude/plans/posing_sequence.md`, "The work, ordered by dependence"; the ruled names are its "The ontology as it stands", "Names ruled". This plan was written at scheduling time: every count below is `[M]` 2026-09-29 by `git grep -n '\b<symbol>\b' -- orpheus tests derivations tools` (lines, files; word-bounded; `docs/` not counted) and every done-when is a hypothesis to re-measure when the unit opens (plan-authoring SCHEDULED-DONE-WHEN).

## Goal, in the domain's terms

A grep over a family's token returns the whole family and nothing else. Each existing symbol below is a loose holder of a word whose perfect match is another object, or spells a structure the tree has since named; renaming it is a retirement, not a cosmetic change. Depends on nothing; nothing in it changes a number.

## The renames

| from | to | why (the ruling) | tree `[M]` |
|---|---|---|---|
| the certificate's per-group relative norm `‖R_g‖/‖Q_g‖` (read as `certificate.balance`, passed as `balance_defect=`) | `relative_residual`; `balance_defect(w)` stays the signed pairing `⟨w, Aψ − q⟩` | one word was on two quantities (second clean attack, D14; ruled 2026-09-28) | `certificate.balance` 16 lines, 3 files; `balance_defect` 30 lines, 14 files |
| the convergence property `rate` | `contraction_rate` | `rate(cells, w)` is the reading verb | not counted: `rate` is too common for a word grep; an AST census of the attribute is owed at opening |
| `Axis.generator`, `generator_as` | `Axis.induced_by`, `induced_by_as` | "generator" is the time generator's; `provenance` collides with two `Provenance` classes (ruled 2026-09-28) | `generator_as` 31 lines, 8 files |
| `GeneratorWithdrawn`, `withdrawn_generator` (the reference "generator") | the `derivation` family | the `derive_*` family's own word (ruled 2026-09-28) | 24 + 54 lines, 4 + 9 files; coordinate with `reference_cache.md` P3/P4, which own these types |
| `direction_sign` | `octant_signs` (a 1-tuple in 1-D) | one family word (ruled 2026-09-28) | 79 lines, 6 files |
| `DiscretizationScheme` | `SpatialScheme` | the discretisation is the product of all axes' frames | 46 lines, 18 files |
| `spatial_closure` | `scheme` | "closure" is the angular axis's word | 95 lines, 6 files |
| `KernelGauge`, `LossKernelGauge`, `LossKernelBasis`, `predicted_kernel_dimension`, `loss_kernel_gauge` | `NullspaceGauge`, `LossNullspaceGauge`, `LossNullspaceBasis`, `predicted_nullspace_dimension`, `loss_nullspace_gauge` | the integral/Markov/reproducing/scattering senses keep "kernel"; the null-space sense cedes (ruled 2026-09-28) | 10 + 46 + 12 + 20 + 50 lines |
| `cell_kernel_batch`, `residual_kernel_batch`, `residual_kernel_batch_transpose`, `has_transpose_kernel` | `update_batch` and its family | the Protocol's own verb is `update` (ruled 2026-09-28) | 87 + 79 lines (the first two), 24 + 21 files |
| the module `orpheus/numerics/projection.py` (holds `AnalysisOperator`/`ReconstructionOperator`, the two faces) | folded into `frame.py` | the module is not the projection | 1 module |

Excluded, with where they go: `direction_idx` is a type change (unit 1); `gram_inverse`/`CrossGramInverse` retire into `radon_nikodym` (unit 2); `EigenPosing`/`SourcePosing` become layer-2 types (unit 6); design-only names (`Stepped`, `CellCoefficient`, `physical_point`, `NuclideDensity`) are minted where their objects land.

## How

One commit per family, each a `retirement-audit` (the three searches, the surfaces a symbol grep misses: strings in refusal messages, docstrings, CLAUDE.md, theory pages, plans' present-tense text; the migration of tests and markers). `dead_references` after each. CI green between families. No shims.

## Done-when (hypothesis, re-measure at opening)

Each old identifier greps to 0 over `orpheus tests derivations tools docs` (excluding `docs/_build` and history in plans), `dead_references` reports no new dead target, the full `-O` suite and the Sphinx `-W` build are green.

## Sizing

1–2 sessions `[R]`.
