# Every collapse, projection and section is one verb of one frame hierarchy (posing unit 2)

Status: SCHEDULED, not opened (2026-09-29). Issue #525; charter issue #522. Charter: `.claude/plans/posing_sequence.md`, "The work, ordered by dependence"; the ruled content is "The ontology as it stands", "The frame and kernel hierarchy". Counts `[M]` 2026-09-29 by `git grep -n '\b<symbol>\b' -- orpheus tests derivations tools`; done-whens are hypotheses re-measured at opening. Depends on unit 0 (`projection.py` folded into `frame.py`).

## Goal, in the domain's terms

A condensation, a moment reduction, a section and a projection are each the same small set of operations on a frame (a basis bound to a measure): analysis, reconstruction, the dual analysis the test/trial pair determines, the idempotent projection, and the density of one pushed-forward measure against another. Today these are spelled several ways, some hand-rolled; after this unit each is one verb, the weight is state passed to it, and the copies retire.

## Scope (ruled 2026-09-28 unless marked)

1. **`FrameBase` gains the coefficient-side measure ν** (declared: the w-frame, Balazs, Antoine & Grybos 2010 Def. 4.1; or the Gram: Parseval). Frame bounds and the dual frame's coefficients depend on ν; the projection does not (`[M]` rows 0.5 against 1).
2. **The discipline chain** `PetrovGalerkinFrame ⊃ GalerkinFrame ⊃ HarmonicFrame`; `PetrovGalerkinFrame`'s docstring states its consumers (today the condensation collapses at `mixture.py:565, 603`, which move to `radon_nikodym`).
3. **The basis hierarchy**: the partition-of-unity basis is the parent, `IndicatorBasis` its one-hot specialisation; today the edge is inverted, `OverlapBasis(IndicatorBasis)` (`OverlapBasis` 39 lines, 11 files).
4. **The verbs**: `analysis` (calling `pushforward`: one body, gated), `reconstruction`, `dual_analysis = (MR)⁻¹M`, `projection = reconstruction ∘ dual_analysis`, and `radon_nikodym(numerator, denominator)`, living on the object that holds both pushforwards.
5. **Retired**: `project`, `average` (a nested closure, `mixture.py:569`), `marginalize` (a comment word), `gram_inverse` (29 lines, 9 files) and `CrossGramInverse` (12, 4) into `radon_nikodym`'s denominator, `_per_pair` (3, 1), and `WeightedIndicatorBasis` (51, 12) as `analysis ∘ M_w`. The weight (φ, φ*, φ* ⊙ φ, a production density, 1/E) is state, never folded into a measure.
6. **Signed weights**: conserved for any signed weight; refused only on a vanishing pushed denominator; a per-region convexity report.
7. **The derived coupling**: `T = analysis_coarse ∘ section_fine` over `w(E) dE` (`[M]` 4.3e-16 against `overlap_to`); `_overlap_table`'s hand loop derived from it and retired; `MeasureCoupling`, if minted, only this cached composite. Generalized Energy Condensation (#275) is the same composite with a polynomial coarse frame.
8. **Condensation provenance**: a coarse object records its parent measures; the two-step defect `M_c (P_m − I) S_f` (vanishes if the second kernel is one-hot, not only if).
9. **The derived angular section**: ν on the scalar marginal stays the physical volume, and `E = R^H (R R^H)⁻¹` (the canonical right inverse; `/Σw` or `G⁻¹` derived, not a convention). Probe: `scratch/posing_sequence/open_items/nu_probe/`.
10. Related issues to read at opening: #492 (`Basis.mass_matrix` twins the frame's Gram), #494 ("tight" names two properties), #403 and #369 (measure equality), #316.

## Gates (from the charter's "frame and collapse machinery" list)

`T == analysis_c ∘ section_f`; the two-step defect identity; `radon_nikodym` bit-identical at every caller it absorbs (`[M]` 2026-09-28: 15 `.project(` call sites, 13 in `material_xs_field.py`, 2 in `mixture.py`; `_per_pair` counted once; the 2 coarse-denominator ratios with their slot typed coarse); the signed-weight report with a diamond-difference positive control; a gate that fails on a non-positive re-weighting fed through the old test-weight path; the derived section against `section.apply` to a few ULP; `analysis` on a one-hot table against `pushforward`.

## Done-when (hypothesis)

The retired symbols grep to 0; every gate above lands with its first red; the full `-O` suite, the Sphinx `-W` build and `dead_references` are clean; the frame theory page carries the hierarchy and the verbs with their derivations.

## Opening obligations

Re-run the `.project(` census; read `orpheus/numerics/frame.py`, `basis/`, the frame theory pages; W3 (surgical: the main agent writes, test-architect specifies the gates first).

## Sizing

3–4 sessions `[R]`.
