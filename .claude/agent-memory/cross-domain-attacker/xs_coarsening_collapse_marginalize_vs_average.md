---
name: xs-coarsening-collapse-marginalize-vs-average
description: STRUCTURAL VERDICT on the shape of ORPHEUS's two XS-coarsening verbs (spatial homogenize P3 merged, energy condense P5 draft). They are ONE operation Collapse(axis, table, test_weight, normalize?) realized as PetrovGalerkinFrame.project — NOT two verbs. THE load-bearing finding: the "1-frame vs 2-frame" and "χ production-weighted vs bare @T" asymmetries are NOT "same slot ± a weight" — they are (a) χ collapses DIFFERENT AXES in the two verbs (homog=spatial-average of χ; condense=energy-BIRTH-marginalization of χ) and (b) two different MORPHISMS: average = G⁻¹·M (project, normalize=True, conserves a RATE) vs marginalize = M (bare @T / analysis, normalize=False, conserves a PROBABILITY/MASS Σχ=1). ONE machinery: analysis (un-normalized M) optionally followed by G⁻¹. project bundles both; chi@table reaches AROUND project to get analysis-alone. Exposing marginalize(=analysis) vs average(=project) DISSOLVES the frame-count asymmetry. Q4: the disjoint/PoU/dense-Gram distinction belongs DECLARED ON THE BASIS (IndicatorBasis→DIAGONAL, OverlapBasis→POU, future GEC→DENSE), not hardcoded as FrameBase.gram's row-sum probe + 30-line precondition docstring (a precondition the type can't violate wants to be a type — same landmine as WeightedIndicatorBasis raising on its synthesis side). Q2 hazard: n-D OverlapBasis is a 2^d-dense OUTER-PRODUCT table (corner overlaps), INVISIBLE on the ordered ≤2-straddle energy axis — must be 2-D-tested before homogenization adopts OverlapBasis. Q5: measure-view & basis-view come from OPPOSITE ENDS (fine.as_measure + coarse.as_basis), and the mismatch/overlap table is a BINARY (fine,coarse) construction step (GroupCondensation IS the energy instance), NOT a unary view. Read on main, clean tree — file:line ground truth (L-005 N/A).
metadata:
  type: project
---

# XS coarsening: ONE Collapse operation; the asymmetries are axis + morphism, not weight

DESIGN VERDICT feeding the P5.5 reshape (`CondensationFrame(PG) ∥ HarmonicFrame`,
the data-folder cleanup that resolved into an architectural reshape). Read on
`main`, clean tree — `frame.py` / `basis/{indicator,overlap,weighted_indicator}_basis.py`
/ `data/macro_xs/mixture.py` / `transport/mesh/material_xs_field.py` /
`sn/solution.py` / `data/energy_grid.py` file:line is GROUND TRUTH (L-005 N/A).
Extends `homogenization_measure_derivation_frames.md` (which DERIVED the φV
measure); THIS memo answers "what is the SHAPE of the machinery across BOTH
coarsening verbs," post-landing of the `FrameBase → PetrovGalerkinFrame →
GalerkinFrame` TYPE hierarchy (the discipline-is-a-type ruling landed in
`frame.py:28-56`).

## THE load-bearing finding (Q3) — marginalize vs average, not weight=1

The brief's comparison table framed two asymmetries as "same slot ± a weight":
"1-frame (condense) vs 2-frame (homog)" and "χ = bare @T (condense) vs
production-weighted 2nd frame (homog)". BOTH framings are STRUCTURALLY WRONG.
The truth is two orthogonal facts:

1. **χ collapses DIFFERENT AXES in the two verbs.** Homogenization χ
   (`material_xs_field.py:409`, `emission_frame.project`) collapses the SPATIAL
   axis — two materials merging into one coarse cell ⟹ spatially average χ,
   production-weighted. Condensation χ (`mixture.py:307`, `chi @ table`)
   collapses the BIRTH-ENERGY axis — χ_g over fine birth groups marginalizes to
   χ_G over coarse birth groups. The space verb NEVER marginalizes χ over
   energy; the energy verb NEVER averages χ over space. Comparing them as one
   slot is an axis category error.

2. **average ≠ marginalize — two different MORPHISMS.**
   - `project = G⁻¹·M` (`frame.py:295-315`) = analysis THEN divide by the Gram.
     A rate-AVERAGE. Conserves a reaction RATE `⟨T·w, Σ⟩`. Used for every Σ
     channel (both verbs) and spatial χ.
   - `chi @ table` = `M` alone (= `frame.analysis.apply`, the un-normalized
     analysis; for a PoU table `@T` IS the analysis). A probability-SUM
     (marginalization). Conserves `Σχ=1`. NO Gram division — a weight=1
     `project` would compute `Σχ_g / n_G` (divide by group COUNT) and BREAK
     `Σχ=1`. So energy-χ is NOT the weight=1 degenerate of `project`; it is a
     DIFFERENT frame operation (analysis without the `G⁻¹` normalization).

**The unification that dissolves the asymmetry:** every channel collapse is
`(test_weight, normalize?)` on the collapse axis —
Σ = (flux, normalize=True); spatial-χ = (production, normalize=True);
energy-χ = (unity, normalize=False); `[g_from,g_to]` matrix = (flux on source,
normalize=True) + (unity on sink, normalize=False, the `@T` sink-sum
`mixture.py:293`). ONE machinery: `analysis` (un-normalized `M`) optionally
followed by `G⁻¹`. `project` bundles both; `chi @ table` reaches AROUND `project`
to reach `analysis`-alone. The honest API exposes BOTH —
`marginalize = frame.analysis` (conserve mass/probability) vs
`average = frame.project` (conserve rate) — and the "1-vs-2-frame" count
EVAPORATES: it was never about frame count, it was about whether each channel's
collapse-axis carries a conserved RATE (average) or a conserved MASS (sum).

First test (DISCRIMINATES, the matrix-channel transpose): 2g→1g condensation
with non-trivial `Sig2` + spatially-varying flux; assert collapsing a
`[from,to]` matrix as `project(Σ @ T)` (sink-sum THEN source-average) ≠ the
order-swapped `project` then `@T` (which double-counts the Φ_{g_from}
normalization). A "collapse a matrix = project both axes symmetrically" impl
gives the wrong number → REDs. Proves sink (marginalize) and source (average)
are NOT interchangeable.

## Honest type decomposition (Q1) — 3 layers, top one is DUPLICATED

ONE operation `Collapse(axis, table, test_weight, conserved_functional)`,
realized as `PetrovGalerkinFrame.project`. Three layers:
1. `PetrovGalerkinFrame.project` (numerics) — axis-agnostic `G⁻¹M`. SHARED. ✓
2. The axis adapter — supplies (trial-basis-from-table, measure, test-weight).
   Space: `mesh.indicator_basis()` + `volume_measure` + flux. Energy:
   `GroupCondensation.indicator_basis()` + counting measure + spectrum. Parallel. ✓
3. The XS-channel taxonomy — "given the collapse axis, route each channel
   (vector/[from,to]-matrix/χ) to its (weight, normalize?) and its
   axis-occurrence count." **DUPLICATED**: `material_xs_field.py:398-409`
   (spatial) and `mixture.py:288-308` (energy) are two transcriptions of one
   taxonomy. THE thing to unify (Cardinal Rule 2). The two entry points
   (`homogenize`/`condense`) stay distinct (axis + output container) but should
   DELEGATE to `collapse_channels(frame, axis_role)`.

## Q2 — OverlapBasis IS the correct general trial; the n-D hazard the energy axis hides

`IndicatorBasis` IS the one-hot `T∈{0,1}` degenerate of `OverlapBasis` (proven by
construction: every contraction is a pure function of the table,
`overlap_basis.py:16-26`). Generalizing homogenization (drop the aligned-edges
assert `solution.py:398-410`, accept an `OverlapBasis` from n-D box-fractions) is
structurally sound. TWO hazards:
- **Hazard 1 (the real one, INVISIBLE in 1-D):** n-D box-intersection is a
  2^d-DENSE outer-product table. The 1-D energy table has ≤2 nonzeros/row (ordered
  axis, ≤2 straddles). An n-D fine cell straddling a coarse CORNER overlaps up to
  2^d coarse cells (4 in 2-D, 8 in 3-D); the fraction is `∏_d f_d` ONLY on a
  tensor-aligned coarse mesh. `IndicatorBasis.evaluate`'s per-axis searchsorted +
  ravel (`indicator_basis.py:172-181`) gives ONE cell/point (one-hot, nested) — the
  fractional n-D table is a per-axis-overlap OUTER PRODUCT then ravel, NOT a
  per-axis searchsorted. First test: 2-D fine cell straddling a 4-corner ⟹ row has
  FOUR nonzeros = `f_x ⊗ f_y`, matched to brute-force polygon-clip areas; a
  per-axis-only impl gives ≤2 → REDs. MUST be 2-D-tested; a 1-D test cannot catch
  it (outer product trivial in 1-D).
- **Hazard 2 (resolved):** the `WithinGroupSpectrum` (1/E lethargy) apportionment
  has NO spatial analogue (space is pure geometric volume fraction — the fine cell
  IS the resolution, no sub-cell flux model). It correctly lives in the ENERGY
  adapter (`energy_grid.py:259`), NOT in `OverlapBasis` (which is correctly
  axis-agnostic, consuming only the finished table). No change.

## Q4 — declare the Gram structure ON THE BASIS (make-illegal-states-unrepresentable)

`FrameBase.gram` (`frame.py:291-293`) unconditionally runs the row-sum probe
`analysis(reconstruction(ones))`, then SPENDS 30 lines (`frame.py:271-289`)
carrying the precondition as PROSE: valid iff disjoint (diagonal Gram) OR PoU
(`R·1=1` ⟹ row-sum collapses to Φ_R), but a tapered/GEC-rank>0 basis (#275)
needs the dense `(MR)⁻¹M` solve it does NOT build. A 30-line caveat on a 3-line
body = a precondition that wants to be a TYPE. FIX: the basis DECLARES its Gram
structure (`IndicatorBasis`→DIAGONAL, `OverlapBasis`→PARTITION_OF_UNITY, future
GEC→DENSE); `project` dispatches — diagonal/PoU → reciprocal probe; DENSE →
dense solve OR RAISE. NEGATIVE test (DISCRIMINATES): a stub basis declaring DENSE
bound in a `PetrovGalerkinFrame` ⟹ `.project` RAISES, NOT silently returns the
wrong row-sum number. Current `FrameBase` returns the wrong number → REDs. SAME
shape as `WeightedIndicatorBasis` raising on its synthesis side
(`weighted_indicator_basis.py:82-87,189-204` — "no consumer ⟹ RAISE, don't
silently delegate a half-consistent op"); the row-sum probe on a dense basis IS
that landmine. `OverlapBasis.mass_matrix` ALREADY knows its inherited diagonal
claim is a lie (`overlap_basis.py:61-69`) but can't fix it without the
declaration — a latent bug awaiting the first `mass_matrix` caller.

## Q5 — duality from OPPOSITE ENDS; the overlap table is a BINARY construction step

The measure-view and basis-view are NOT two views of one grid — they come from
OPPOSITE ENDS of the coarsening: the FINE grid yields the measure (integration
domain, `fine.volume_measure`); the COARSE grid yields the trial basis (coeff
space, `coarse.indicator_basis()`). `homogenize` already does this
(`solution.py:389,423`). So the symmetric API is `fine_axis.as_measure() +
coarse_axis.as_basis()`. The MISMATCH/OVERLAP table is a BINARY `(fine, coarse)`
op — it reads BOTH grids' edges (`energy_grid.py:293-315`) — so it CANNOT be a
unary `coarse.as_basis()`; it is a frame-CONSTRUCTION step. `GroupCondensation`
IS the energy instance of this binary `overlap(fine, coarse, within_group) →
OverlapBasis`. The spatial instance does not exist yet (homogenization asserts
nested). Honest API: unary `coarse.indicator_basis()` (nested one-hot via
searchsorted) ⊂ binary `overlap(fine,coarse)` (non-nested fractional; energy =
interval-overlap, space = box-intersection, Q2's n-D outer product). The
P5.5-proposed `Frame(trial, measure=fine.as_measure(), test=Weighted(trial,w))`
shape is correct, with the table a binary input.

## Cross-method pollination (concrete, fail-able)

- **CP/MoC region-merging** — CP region-averaged XS = `homogenize` with the CP
  region flux as test weight; MoC FSR-merging = `homogenize` with an
  `OverlapBasis` (FSR straddling a coarse region by TRACK-LENGTH fraction = the
  spatial transpose of a fine group straddling by lethargy fraction). Test: CP
  region Σ == `Solution.homogenize` onto the same partition w/ CP region flux,
  `array_equal`; a CP volume-only (no-flux) average differs whenever flux varies
  in the merged region → REDs.
- **Eigenvalue/adjoint-weighted (the LIVE extension, #48/P6)** — test weight
  generalizes φ→φ* (SAME `WeightedIndicatorBasis`, different weight); the
  Galerkin degenerate (φ*=φ) becomes the special case, PG the general — the
  campaign's stated reason for the PG-base hierarchy. Test: adjoint-weighted
  homog preserves the BILINEAR `⟨φ*,Σφ⟩` while forward-flux does NOT (differ
  whenever φ*≠φ, i.e. any real reactor) → asserting they agree is theatrics;
  assert `⟨φ*,Σ_eff φ⟩==⟨φ*,Σφ⟩` REDs the forward-only impl.
- **GET/discontinuity factors (Smith 1986)** — flip condition: when
  leakage/current preservation joins rate preservation, a SECOND test functional
  `{Ω·n on ∂R}` (boundary-trace test space) binds to the same coarse trial.
  Test: preserve BOTH volume rate AND surface current J·n; rate-only frame gets
  coarse current wrong → REDs. Not live; named for the next consumer.

## Refuted frames (durable UNEXPLORED for the coarsening problem class)

- **Tensor networks / MPO** — the [from,to] collapse is rank-2, n-D overlap is
  rank-d outer product; both FIXED low rank, no bond-dimension KNOB (≥3 varying).
  Bond-dim-degenerate, not a network (L-001). The ONLY real rank knob is GEC
  (Q2-frame-2), and there it is a polynomial DEGREE, not a tensor bond dim.
- **Homology / Kan extension / Radon-Nikodym / differential geometry /
  Wiener-Hopf** — all refuted in `homogenization_measure_derivation_frames.md`
  with the same reasons (no ∂²=0; pushforward already captures the Kan role; R-N
  derivative IS M_φ; no curvature term; wrong solver family). No new payoff this
  pass (L-001).

## Elegance assessment

- Structure-exposing (STRONG): the marginalize-vs-average split EXPOSES that
  `chi @ table` reaching around `project` is `analysis`-without-`G⁻¹`; the
  [from,to] two-step is "the collapse axis occurs with multiplicity 2, one
  occurrence per (sink=sum, source=average)."
- Structurally-simpler (STRONG): one `(weight, normalize?)` machinery replaces
  the 1-frame-vs-2-frame + bare-@T special cases; the channel taxonomy unifies
  two transcriptions into one.
- Expressive: `frame.project`/`frame.analysis` already serve angular SH +
  spatial homog + energy condense — three axes, two verbs (average/marginalize).
- Algorithmic-advantage: declared-Gram dispatch routes GEC rank>0 to the dense
  solve seam without a runtime prose-precondition; n-D overlap vectorizes as the
  per-axis outer product.
