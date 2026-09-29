# Independent defects the posing work found, fixed where they stand (posing unit 1)

Status: SCHEDULED, not opened (2026-09-29). Charter: `.claude/plans/posing_sequence.md`, "The work, ordered by dependence". Counts `[M]` 2026-09-29 by `git grep -n '\b<symbol>\b' -- orpheus tests derivations tools` unless marked; done-whens are hypotheses re-measured at opening.

## Goal, in the domain's terms

Each item is a defect or a twin that stands on its own: fixing it needs no object the later units build, and leaving it lets a wrong answer or a duplicated quantity survive while the architecture is carved. Each is an independent commit (or its own small issue); the order among them is free.

## Items

1. **`direction_idx` is a type change, not a rename.** On the cylinder arm it is the within-level azimuthal index, not an ordinate (`reduced_operator.py:463-476`; 6 live `mu_level_idx=` sites, per the second clean attack's elegance memo). `streaming_terms` takes the global ordinate only; `(level, m)` is derived once from `quadrature.level_indices`. `direction_idx` 129 lines, 13 files `[M]`.
2. **The evidence readers become exhaustive.** Two production readers branch `isinstance(x, Measured)` with a silent else (`convergence.py:1686`, `loss_kernel_gauge.py:508`, per the elegance memo; re-locate after unit 0's rename): an unknown evidence kind vanishes. They become exhaustive `match` statements. Precondition for `Estimated` (unit 7).
3. **`ScaleGauge.displacement` keeps the state's dtype** (today `float(...)` turns a complex reading into a wrong number under a one-time `ComplexWarning`; the `FundamentalMode` attack's memo). Complex support or a refusal.
4. **The pencil's zero-scale law**: `0 · T` is the zero operator, so `pencil.at(0)` is the unshifted operator (production refuses `ScaledOperator(0, ·)` today, per proto2's memo).
5. **The string comparison on `spectral_map.name` retires** (`outcome.py:192`; 2 lines, 1 file `[M]`).
6. **A present-tense-false comment**: `operator.py:200-204` says production "dispatches on" block roles; `[M]` (second attack) 0 production discrimination sites. Correct it now; the roles themselves retire in unit 4.
7. **`rate(cells, w)` absorbs its twins**: 9 production definitions (3 of them `compute_production_rate`, 38 lines, 12 files `[M]`), `ReactionRateFunctional` and `rayleigh` (#472, #462 relate). May instead open unit 6 if the reading verb needs the system; decide at opening.
8. **#520's remaining sites** (the restrict-scatter projectors, the einsum retractions) and **#521** (the complex Gram probe).

## Gates

Each item lands with the input that reddens its gate (plan-authoring §6c): a complex reading for item 3; a zero scale for item 4; a new evidence kind for item 2 (a test-only kind that the exhaustive `match` must refuse or handle); the cylinder arm's within-level index for item 1.

## Done-when (hypothesis)

Items 1–6 and 8 closed with their gates; item 7 either closed or moved to unit 6 by an explicit note in this file.

## Sizing

2–3 sessions `[R]` in total.
