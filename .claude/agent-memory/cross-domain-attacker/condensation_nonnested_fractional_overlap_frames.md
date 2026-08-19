---
name: condensation-nonnested-fractional-overlap-frames
description: DESIGN VERDICT for energy CONDENSATION as a frame projection, NON-NESTED case (421→WIMS69/172, straddling fine groups). RESOLVES the standing-campaign tension "flux is a TEST-weight, never a measure / Measure=axis+metric, never discipline" against the user's floated "spectrum-as-measure least-squares" alternative. VERDICT condensation is PETROV-GALERKIN with a FRACTIONAL-OVERLAP trial table T[g,G]∈[0,1] (partition of unity, rows sum 1) — the SAME discipline as spatial homogenization, NOT a new LeastSquaresFrame. The user's "spectrum-as-measure least-squares" is REFUTED as the native frame on three independent structural grounds: (1) it folds φ into the metric, which the campaign already ruled breaks under adjoint-weighted (eigenvalue-consistent) condensation where the preserved functional is the BILINEAR ⟨φ†,Σφ⟩ and test=φ†·1_G≠trial — frame.py:48-56 / weighted_indicator_basis.py:24-34 are the load-bearing docstrings; (2) for the piecewise-constant (P0/indicator) coarse basis the least-squares fit COINCIDES with the flux-weighted average — so it is not a DIFFERENT map, it is the SAME PG map misattributed to a least-squares discipline (LeastSquaresFrame would only differ for a NON-indicator richer coarse basis, which cross sections never want — a P1 coarse XS is not rate-meaningful); (3) LeastSquaresFrame's defining trigger is test=A·trial (a DENSE SPD Gram needing a real solve), absent here (disjoint/fractional indicators give a DIAGONAL Gram, a reciprocal not a solve). CONFIRMS the brief's partition-of-unity claim: frame.project on a row-sum-1 T returns the rate-preserving average UNCHANGED because gram = analysis(reconstruction(ones)) = Φ_G = Σ_g T[g,G]φ_g is exactly the right denominator (the diagonal cross-Gram). The within-group flux model f_{g,G} is the TRIAL TABLE (a basis/overlap choice on a common fine refinement), NOT a measure — it sets the column-split of a straddling fine group's indicator, which is a geometry-of-the-partition datum, identical in kind to the spatial cell-membership table. SINGLE most elegant formulation: generalize IndicatorBasis.evaluate from one-hot (searchsorted) to a fractional-overlap table (clip+area-fraction), nested falls out as the {0,1} special case, rate-preservation and the exact-group-sum-when-nested are theorems of partition-of-unity, ZERO new discipline / ZERO new frame subclass. Cross-method: the overlap operator is the SAME fine→coarse restriction π as spatial homogenization (a different searchsorted on a different axis), and is the discrete pushforward measure.pushforward / partition_by already shipped; non-nested overlap = a soft pushforward (one fine atom splits across coarse cells). SUPERSEDES homogenization_measure_derivation_frames.md's "Galerkin-in-L²(φV) is the native reading" (the campaign REVERSED that: PG-with-flux-as-test, not Galerkin-in-reweighted-metric — the adjoint-weighted case is the falsifier the prior memo missed). Read on main @ feature/sn-energy-condensation context, clean tree — frame.py/indicator_basis.py/weighted_indicator_basis.py/measure.py/solution.py/energy_grid.py/material_xs_field.py + .claude/plans/archive/p5_0_condensation_gate.md file:line GROUND TRUTH (L-005 N/A).
metadata:
  type: project
---

# Energy condensation, non-nested: PG-with-fractional-overlap, NOT a LeastSquaresFrame

DESIGN VERDICT for the non-nested condensation extension (421→WIMS69/172).
Feeds GitHub #274 / branch `feature/sn-energy-condensation` (P5.0 gate already
written: `.claude/plans/archive/p5_0_condensation_gate.md`). The CURRENT SUT
(`GroupCondensation`, `energy_grid.py:154`) handles only the NESTED case — its
`__post_init__` REJECTS non-contiguous maps (`:198`), the containing-interval
rule (`:206`) produces a one-hot `T` (`:220`). This memo is the forward design
for the straddle case the brief poses.

Read on `main` with the `feature/sn-energy-condensation` files present, clean
tree — file:line is GROUND TRUTH (L-005 N/A).

## The tension, resolved (the headline)

Standing campaign ruling: **"flux is a TEST-weighting, NEVER a measure; Measure
= axis + metric, never discipline."** User floated: **"spectrum-as-measure
least-squares projection onto a chosen coarse basis."** These are in apparent
conflict. The resolution is NOT a compromise — the user's framing is **refuted
as the native frame on three independent structural grounds**, and condensation
stays **Petrov-Galerkin with the flux on the test side**, exactly like spatial
homogenization. The non-nested case adds NOTHING to the discipline — it
generalizes the TRIAL TABLE from one-hot to fractional (partition of unity).

This SUPERSEDES `homogenization_measure_derivation_frames.md`, which argued the
native reading is "Galerkin-in-L²(φV)". The campaign REVERSED that (frame.py
docstring `:34` "REVERSES the earlier projection.py discipline-ABC design";
weighted_indicator_basis.py `:24-34`): the decisive falsifier the prior memo
MISSED is **adjoint-weighted (eigenvalue-consistent) homogenization**, where the
preserved functional is the BILINEAR `⟨φ†, Σφ⟩` with `φ† ≠ φ`, so `test = φ†·1_R
≠ trial = 1_R` — genuinely Petrov-Galerkin, with NO metric fold available
(folding `φ†` into the measure and `φ` into the field is not a symmetric inner
product). Forward-flux reaction-rate-only condensation is the DEGENERATE `φ†=φ`
case where the fold happens to work; the TYPE must encode the general case.

## Why the user's least-squares framing is refuted as the native frame (3 grounds)

1. **It folds φ into the metric — the campaign's load-bearing refutation.**
   "Spectrum-as-measure" = put φ in `DiscreteMeasure.weights`. That is the
   Galerkin-in-L²(φ) reading. It breaks under adjoint-weighting (above): a
   measure carries ONE metric, but the bilinear `⟨φ†,Σφ⟩` needs two distinct
   weightings on the two sides. The frame TYPE (`PetrovGalerkinFrame` with an
   explicit `test_basis = WeightedIndicatorBasis(indicator, φ)`) is exactly the
   encoding that survives the general case. (frame.py:48-56.)

2. **For a P0/indicator coarse basis, least-squares COINCIDES with the
   flux-weighted average — it is the SAME map, not a different one.** The
   weighted-L² least-squares fit of `Σ` onto `span{1_G}` solves the normal
   equations `⟨Σ−Σ_c, 1_G⟩_w = 0`; with disjoint (or fractional partition-of-
   unity) indicators the Gram is diagonal and the solution is `Σ_G = Σ_g w_g
   Σ_g / Σ_g w_g` — the flux-weighted average verbatim when `w=φ`. So invoking
   "least-squares" does not select a NEW frame; it re-derives the PG average
   under a different name. A `LeastSquaresFrame` (named-but-unbuilt seam,
   `frame.project` docstring `:288` "the dense solve is the (unbuilt)
   least-squares seam only") would DIFFER only for a richer coarse basis where
   the Gram is non-diagonal (polynomials/smooth functions per coarse group). Cross
   sections never want that: a P1/quadratic coarse XS is not reaction-rate-
   meaningful (the rate-preservation functional is `∫_G φΣ`, a single scalar per
   group — one DOF, i.e. P0). Richer coarse bases buy nothing for XS.

3. **LeastSquaresFrame's defining trigger is structurally absent.** A genuine
   least-squares frame has `test = A·trial` for some operator `A` (the normal
   equations `AᵀA c = Aᵀ b`), producing a DENSE SPD Gram that needs a real solve
   (plan `frame_projection_machinery.md:66`). Condensation's Gram is DIAGONAL
   (disjoint/fractional indicators) → a per-group reciprocal, never a solve. The
   trigger for minting `LeastSquaresFrame` is a consumer with a non-diagonal
   normal-equations Gram; condensation is not it.

## The three sub-questions

- **Is non-nested condensation the consumer that builds LeastSquaresFrame?**
  NO. It is PG with a fractional-overlap trial table. `LeastSquaresFrame` stays
  unbuilt until a non-indicator (dense-Gram) consumer arrives.
- **Flux in the measure or the test-weight?** TEST-WEIGHT, same as
  homogenization. There is no exception. The "spectrum-as-measure" reading is
  the φ†=φ degenerate that the type deliberately does not privilege.
- **Most elegant frame-native formulation?** Generalize `IndicatorBasis.evaluate`
  from one-hot to a fractional-overlap (partition-of-unity) table. Everything
  else (PG discipline, `frame.project`, the diagonal Gram, rate-preservation,
  the exact group-sum when nested) is unchanged and falls out as theorems.

## CONFIRMED: the brief's partition-of-unity gram claim (the load-bearing check)

The brief claims `frame.project` already yields the rate-preserving average for
a partition-of-unity `T` because `reconstruction.apply(ones)` row-sum = 1 gives
the diagonal Gram `Φ_G = Σ_g T[g,G] φ_g`. **CONFIRMED.** Trace (frame.py:271-273,
`gram`): `ones = np.ones(basis_space.shape)` (shape `(n_coarse,)`); `diagonal =
analysis.apply(reconstruction.apply(ones))`. With the PG frame
(test = `WeightedIndicatorBasis(indicator, φ)`, trial = fractional `IndicatorBasis`):
- `reconstruction.apply(ones)` = trial `synthesize`/`reconstruct` =
  `Σ_G T[g,G]·1 = rowsum_G T[g,G] = 1` (partition of unity) → a constant-1 fine
  field. (indicator_basis.py:219 `reconstruct` = identity-dual broadcast.)
- `analysis.apply(1)` = test `analyze` = `Σ_g w_g·φ_g·T[g,G] = Σ_g φ_g T[g,G] =
  Φ_G` (weighted_indicator_basis.py:140-151, counting `w=1`). The diagonal cross-
  Gram `G_G = Σ_g φ_g T[g,G]`. EXACTLY the denominator.
Then `project(Σ) = gram.apply_inverse_metric(analysis.apply(Σ)) = (Σ_g φ_g T[g,G]
Σ_g)/(Σ_g φ_g T[g,G])` — the fractional rate-preserving average, UNCHANGED. The
one-hot case is `T∈{0,1}` → the current exact group-sum. The off-diagonals of the
cross-Gram are NOT structurally zero for a fractional `T` (two coarse columns can
share a fine row), so the `gram` shortcut "row-sum IS the diagonal" (frame.py:264)
needs a guard: it holds iff a fine atom that contributes to coarse G also
contributes ONLY to G in the Gram off-diagonal sense — TRUE for partition-of-
unity here because `G_{GG'} = Σ_g φ_g T[g,G] T[g,G']` and the rate-preserving
average uses ONLY the diagonal `G_GG`. The non-diagonal `G_GG'` exists but
`project` (diagonal-inverse) ignores it — which is CORRECT for rate preservation
(each coarse group's rate is its own functional). FIRST TEST below pins this.

## The within-group flux model f_{g,G} is the TRIAL TABLE, NOT a measure

A straddling fine group `g` apportions fraction `f_{g,G}` of its rate to coarse
group `G`. The brief asks: measure choice or basis choice? **TRIAL-TABLE (basis)
choice.** `f_{g,G}` IS the column-split of `g`'s indicator across the coarse
columns it overlaps — i.e. `T[g,G] = f_{g,G}`, the fractional membership. It is a
geometry-of-the-partition datum (how the fine bin's energy interval overlaps the
coarse bins, optionally area-weighted by a within-group flux shape), IDENTICAL in
kind to the spatial cell-membership table (`indicator_basis.py:142`,
`searchsorted` one-hot). The within-group flux MODEL (flat lethargy / 1/E /
linear) sets HOW the area is split — it is a refinement of `evaluate`, on the
TRIAL side. It is NOT the measure (the counting/volume measure stays the fixed L²
metric, weights group-independent) and NOT the test weight (the cell/group flux
`φ_g` stays the test side). Three distinct objects, cleanly separated:
- TRIAL table `T[g,G]` (the partition geometry + within-group split) — basis.
- TEST weight `φ_g` (the spectrum) — `WeightedIndicatorBasis`.
- MEASURE (counting `w=1` on the group index) — fixed metric, untouched.

This is the campaign's three-way separation holding under the non-nested
generalization — the elegance payoff is that the straddle case touches ONLY the
trial table, the most local of the three.

## The single most elegant formulation (the deliverable)

**Generalize `IndicatorBasis.evaluate` to return a fractional overlap table; keep
everything else.** Concretely:
- `evaluate(points)` currently: bucket each fine center by `searchsorted` → one-
  hot `T[i,R]∈{0,1}` (indicator_basis.py:172-184). Generalize to: for each fine
  bin, compute the FRACTION of its interval (optionally φ-shape-weighted) falling
  in each coarse bin → `T[i,R]∈[0,1]`, rows sum to 1 (partition of unity).
- (a) handles nested AND non-nested uniformly: nested = the one-hot special case
  (each fine interval ⊆ one coarse interval → `f=1` in one column).
- (b) rate-preserving: partition-of-unity `T` + the `gram=Φ_G` denominator (proven
  above) preserves `Σ_g φ_g T[g,G] Σ_g` per coarse group.
- (c) reduces to the exact group-sum when nested: one-hot `T` → `Σ_G = Σ_{g∈G}
  φ_g Σ_g / Σ_{g∈G} φ_g`, the current behavior.
- (d) reuses `FrameBase.project` with MINIMAL new surface: ONLY `evaluate` changes
  (a fractional table builder). No new discipline, no `LeastSquaresFrame`, no new
  frame method. The membership generalization is the entire delta.

This is the `clean-before-extend` floor in action: the capability (non-nested)
lands as a no-op extension through the one generic body (`evaluate` → fractional
table), not a new arm.

## Cross-method pollination

- **Spatial homogenization ≡ energy condensation under one overlap operator.**
  Both are `frame.project` with a fine→coarse restriction `π` realized as the
  trial `IndicatorBasis.evaluate` table. The ONLY difference is the bucketing
  axis: spatial = `searchsorted` on coarse cell edges (per spatial axis,
  `indicator_basis.py:172`); energy = `searchsorted` on coarse group energy
  boundaries (`energy_grid.py:216`). A fractional `evaluate` serves BOTH — the
  non-nested spatial case (coarse cell edges not aligned to fine edges, the
  AMR/overlapping-mesh case) gets fractional re-binning for FREE the moment the
  energy case builds it. Generalize ONCE.
- **The overlap operator IS the discrete (soft) pushforward.** One-hot `T` is the
  hard pushforward `measure.pushforward(π)` (`measure.py:546`, π = the cell/group-
  locate map). Fractional `T` is a SOFT pushforward: one fine atom splits its
  weight across multiple coarse atoms (the pushforward's "non-invertible φ → atoms
  collapse" note `:572` run in reverse — one atom FANS OUT). `partition_by`
  (`measure.py:673`) is the hard-partition primitive; the soft version is the
  natural generalization (a fuzzy partition / weighted membership). NOT worth
  minting a new measure primitive YET (one consumer); name it as the soft-
  pushforward when a second fractional consumer (non-nested spatial) lands.

## Refuted frames (durable UNEXPLORED)

- **LeastSquaresFrame** — REFUTED as condensation's frame: its trigger (test=A·trial,
  dense SPD Gram, a real solve) is absent (diagonal Gram → reciprocal). It DIFFERS
  from PG only for a non-indicator coarse basis cross sections never want. Flip
  condition: a consumer needing a non-diagonal normal-equations Gram (a continuous/
  spectral coarse representation). Not condensation.
- **Galerkin-in-L²(φ) / spectrum-as-measure (the user's framing + the prior
  memo's)** — REFUTED: folds φ into the metric, breaks under adjoint-weighting
  (the bilinear `⟨φ†,Σφ⟩`, test≠trial). It is the φ†=φ degenerate of PG, not the
  native frame. (This is the campaign's reversal of the prior memo.)
- **General FEM / mixed FEM** — wrong weight class: P0 disjoint/fractional
  indicators give a closed-form diagonal-Gram solve; FEM assembly + non-diagonal
  mass matrix is unnecessary machinery. Fire only for P1 coarse homogenization
  (never, for rate-preserving XS).
- **Conservative remapping / supermesh (Farrell-Maddison)** — the CFD/coupling
  literature's conservative-interpolation frame (Galerkin L² projection on a
  SUPERMESH = the common refinement of source and target). It is REAL and exactly
  the non-nested overlap structure — BUT for the P0→P0 (cell/group indicator)
  case the supermesh Galerkin projection COLLAPSES to the fractional-overlap
  average (the supermesh cells ARE the fine∩coarse overlaps, the area fractions
  ARE `f_{g,G}`). So it adds vocabulary (supermesh = the common refinement), not a
  different map. Name it UNEXPLORED with the flip: a P1-or-higher target would
  need the genuine supermesh-Galerkin (non-diagonal Gram) — the LeastSquaresFrame
  trigger again. (L-001: concrete frame — fractional overlap — already captures the
  P0 win.)
- **Optimal transport / Wasserstein** — tempting via "split one fine atom's mass
  across coarse atoms," but condensation has NO cost/ground-metric to minimize
  (the split is fixed by interval overlap, not by transporting mass optimally).
  No OT lever. (L-001.)
- **Homology / category theory** — no `∂²=0`; the soft-pushforward / left-Kan
  framing of the fine→coarse collapse is fully captured by `pushforward`/
  `partition_by` + the frame pair. (L-001.)

## First tests (each DISCRIMINATES — L-002)

- **FRACTIONAL-RATE-PRESERVATION (the core claim):** build a PG frame with a
  STRADDLING fine group (e.g. fine group 2's interval overlaps coarse G0 by 0.3
  and G1 by 0.7, so `T[2,0]=0.3, T[2,1]=0.7`), a within-group-VARYING φ, assert
  per-coarse-group `Σ_G·Φ_G == Σ_g φ_g T[g,G] Σ_g` (`atol=1e-12`, one-ULP). A
  one-hot impl (assigns fine-2 entirely to its containing-interval coarse group)
  gives a DIFFERENT rate split → REDs. A non-straddling fixture makes fractional ≡
  one-hot → CANNOT discriminate; the test MUST use a straddling group.
- **NESTED-REDUCTION (b/c):** assert that for a NESTED partition the fractional
  `evaluate` returns a `{0,1}` table bit-identical (`array_equal`) to the current
  one-hot `searchsorted` table. A fractional impl that always splits (never snaps
  to {0,1} on alignment) would drift → `array_equal` REDs.
- **GRAM-DIAGONAL-CORRECTNESS (the partition-of-unity gram check):** for a
  fractional `T` with a shared fine row, assert `frame.gram`'s diagonal equals
  `Σ_g φ_g T[g,G]` per coarse G (`array_equal` against the hand einsum), AND that
  `project` ignores the non-zero off-diagonal cross-Gram `Σ_g φ_g T[g,G]T[g,G']`
  (i.e. the result is the per-group diagonal average, not a coupled solve). An
  impl that built the FULL cross-Gram and solved would give a coupled (wrong)
  answer → REDs. This pins the brief's "row-sum IS the diagonal" claim is
  RATE-correct even when off-diagonals are non-zero.
- **PG-NOT-MEASURE (negative, the discipline pin):** assert the condensation frame
  is a `PetrovGalerkinFrame` whose `measure.weights` are the group-INDEPENDENT
  counting ones (`array_equal` to `np.ones`), and that φ lives on the
  `test_basis.weight`, NOT the measure. A "spectrum-as-measure" impl would put φ
  in `measure.weights` (group-dependent) → the negative assertion `measure.weights
  == ones` REDs. This is the structural pin that the flux is the test weight, not
  the metric.
- **ADJOINT-WEIGHTED-IS-PG (the falsifier the prior memo missed):** assert that an
  adjoint-weighted condensation (`test = φ†·1_G`, `φ† ≠ φ`) is expressible as a
  `PetrovGalerkinFrame` (distinct test) but is NOT expressible by folding into a
  single measure metric (there is no `w` such that `L²(w)` reproduces the bilinear
  `⟨φ†,Σφ⟩` when φ†≠φ). The negative half: a Galerkin-in-L²(w) impl cannot match
  the φ†≠φ rate → REDs. (Latent consumer P6; the test is the structural argument
  that PG is the right TYPE even before P6 ships.)

## Elegance assessment

- **Structurally-simpler (strong):** non-nested is a no-op extension through ONE
  body (`evaluate` → fractional table). No new discipline, no `LeastSquaresFrame`,
  no new frame method. The three-way trial/test/measure separation holds verbatim;
  the straddle touches only the most-local object (the trial table).
- **Structure-exposing (strong):** the within-group flux model is EXPOSED as the
  trial-table column-split (the partition geometry), not a measure or a separate
  apportionment step — the same object as spatial cell-membership. The user's
  "least-squares" is exposed as the φ†=φ degenerate of the already-built PG.
- **Expressive (strong):** ONE fractional `evaluate` serves spatial homogenization
  (non-aligned meshes) AND energy condensation (non-nested groups) — the overlap
  operator is axis-agnostic. Nested is the `{0,1}` special case of one formula.
- **Algorithmic-advantage:** the diagonal Gram keeps `project` a per-group
  reciprocal (no solve) even in the non-nested case — the fractional generalization
  does NOT cost a linear solve, unlike a genuine least-squares (dense-Gram) frame.
