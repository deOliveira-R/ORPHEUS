---
name: unified-frame-api-design
description: The UNIFIED projection/reconstruction Frame API design — lifts the prior derivation memos (projection-discipline-hierarchy, homogenization-measure-derivation) from "what the map IS" into "the three API primitives." THE DEEP OPERATION is Christensen's canonical-dual coefficient map; all three consumer patterns specialize it. (1) project_weighted(f, weight=w) = G_w⁻¹·M(w⊙f) = (T*M_wT)⁻¹(T*M_w)f — homogenize/condense, NO reconstruct, the /Φ_R denominator IS G_w⁻¹; diagonal Gram (SH per-ℓ g_C, Indicator per-region m_R) ⇒ apply_inverse_metric NOT solve. (2) conjugate(A)=R∘A∘M typed OperatorProduct face ON the frame; sub-op reconstruct_after(A)=R∘A for the windowed arm; (1/W) is ScaledOperator OUTSIDE per L18. (3) analyze(f)=bare M, keep. GRAM-LOCUS Q2: Option-B (per-basis split, frame KNOWS which) — keep SH's measure-free reconstruct (0-ULP test_scattering_kernel_crosscheck), Indicator's space-metric; the weighted-dual G_w is the verb's, NOT a reconstruct change. CONJUGATE does NOT need the 2-param Din/Cout operator split — just give Λ(LegendreMomentScattering) real domain/codomain=basis_space (deletes the cast, re-arms the guard, pyright unifies). TEST-SIDE Q4: Frame(basis,measure,*,test:Basis|None=None), test=None⇒Galerkin (zero-churn default); is_galerkin=(test is basis) computed property NO class; RETIRE the 2 PG ABCs (#268); build the SEAM not the PG cross-Gram path (0 PG/LS consumers); canonical_measure NOT needed (test-side identity decides). LS Q5 = flag-and-defer: test=A·K is operator-dependent + non-diagonal Gram ⇒ CAP_SOLVE coeff space, breaks the diagonal-reciprocal; test:Basis|None NOT Basis|Operator. FLUX-AS-MULTIPLIER is durable (retires flux_volume_measure from the consumption path); DiscreteMeasure stays 1-D. Condensation=same verb, EnergyGrid.indicator_basis greenfield leaf. Read on main, clean — file:line ground truth (L-005 N/A).
metadata:
  type: project
---

# Unified Frame API — the three primitives (operationalizes the derivation memos)

DESIGN VERDICT for the unified projection/reconstruction machinery (the
homogenization/condensation + scattering-conjugation campaign). Read on `main`,
clean tree — `frame.py`/`projection.py`/`measure.py`/`operator.py`/`basis/base.py`/
`basis/indicator_basis.py`/`space.py`/`scattering.py`/`solution.py` file:line is
GROUND TRUTH (no Nexus/worktree staleness; L-005 N/A; the brief's file:line came
from a worktree — verified the live `orpheus/numerics/` tree). LIFTS
`projection_discipline_hierarchy_frames.md` (the (K,L) recognition + retire-the-ABCs)
and `homogenization_measure_derivation_frames.md` (the L²(φV) Galerkin derivation +
flux-as-multiplier) from derivations INTO concrete API primitives. The two prior
memos answer "what is the map"; THIS memo answers "what are the verbs."

## The deep operation: ONE canonical-dual coefficient map, three specializations

Christensen *Frames & Riesz Bases* §5.1: the frame operator `S = T*T`, canonical
dual `{S⁻¹φ_k}`. The single deep operation is the WEIGHTED canonical-dual coefficient
map. All three brief-patterns specialize it:

    project_weighted(f, weight=w)  ≡  G_w⁻¹ · M(w⊙f)  =  (T* M_w T)⁻¹ (T* M_w) f

- `M_w = diag(w)` = state multiplier (`M_φ`, the coefficient_field_promotion
  multiplier; flux for homogenize, spectrum for condense).
- `M(·) = frame.analysis.apply` = the un-normalized weighted moment `T*W` the
  analysis face ALREADY computes (the numerator).
- `G_w = T*(W⊙w)T` = the WEIGHTED frame operator (weighted coarse Gram). DIAGONAL
  for BOTH bases (SH per-ℓ `g_C`; Indicator per-region `m_R`, disjoint supports) ⇒
  `G_w⁻¹` is `apply_inverse_metric` (reciprocal), NOT `CAP_SOLVE`. The `/Φ_R`
  denominator IS `G_w⁻¹`.

The three API primitives (minimal set):

    frame.analyze(f)               -> coeffs         # = M f          EXISTS, keep
    frame.conjugate(A)             -> LinearOperator # = R∘A∘M        Frame-2
    frame.project_weighted(f, w=…) -> coeffs         # = G_w⁻¹ M(w⊙f) the deep one

`project_weighted(f)` (no weight) = plain canonical-dual analysis (`R∘M` coeffs).
Pattern-1 (bare M, scattering) = keep `frame.analyze` (consumer wants UN-normalized
moments; do NOT route through project_weighted). Pattern-2 (kernel) = `conjugate`.
Pattern-3 (homogenize/condense) = `project_weighted(Σ, weight=φ)`.

Collapses (the user's "couple of lines, nothing hand-rolled" acceptance test):
- scattering `kernel` (was 12-line property + `cast`, `scattering.py:662-673`) →
  `frame.conjugate(Λ)`.
- homogenize per-channel (was the `replace(basis_space, weights=…); analysis.apply;
  apply_inverse_metric` triple-body ×3, `solution.py:422-456`,`:486-492`) →
  `frame.project_weighted(channel_fine, weight=phi)` ONE line each.

## Q2 (Gram locus) verdict — Option B: per-basis split, the FRAME knows which

KEEP the live split: SH bakes `2ℓ+1=4π·g_C⁻¹` into `reconstruct` (analytic,
measure-FREE, `spherical_harmonic_basis`); Indicator defers `1/m_R` to the space
metric (`apply_inverse_metric`, `indicator_basis.py:52-67`). Decision drivers:
- 0-ULP: Option A (G⁻¹ ALWAYS in the coeff-space metric, `reconstruct` bare for all
  bases) CHANGES the SH reconstruct path → risks `test_scattering_kernel_crosscheck`
  (pinned 0 ULP). The measure-free analytic `2ℓ+1` is bit-identical and correct;
  moving it to the metric is churn for no structural gain.
- single-source: the weighted-dual `G_w` the UNIFIED verb needs is NOT the
  reconstruct dual factor — it is the STATE-weighted frame operator, computed inside
  `project_weighted`. So `project_weighted` owns `G_w⁻¹`; `reconstruct` keeps the
  UNWEIGHTED canonical dual (analytic for SH, identity for Indicator). The two
  Gram-inverses are DIFFERENT objects (unweighted dual vs state-weighted dual) — no
  conflict. The frame KNOWS which reconstruct path via the basis (already
  polymorphic). Option B, no SH change.

## Q3 (conjugate) — typed OperatorProduct face ON the frame; sub-op for the windowed arm

`frame.conjugate(A) = OperatorProduct(R, OperatorProduct(A, M))` — change-of-basis
CONJUGATION (similarity transport) of a coefficient-space endomorphism `A: W→W`.
Lives ON the frame (reads both faces+spaces). Windowed arm (re-enters post-M) reuses
via `frame.reconstruct_after(A) = R∘A` (`OperatorProduct(R, A)`), with
`conjugate(A) = reconstruct_after(A) @ analysis`. `(1/W)` = `ScaledOperator`/`__truediv__`
OUTSIDE the kernel at the apply boundary (L18) — `conjugate` returns pre-`1/W`.

AUGMENTATION-a (does conjugate need the 2-param Din/Cout type split?): NO. Runtime
`LinearOperator[V]` + `domain`/`codomain` FunctionSpace props type it end-to-end
(OperatorProduct composability guard checks `R.domain==A.codomain==W` etc). The
`cast(LinearOperator, Λ)` (`scattering.py:663`) exists ONLY because
`LegendreMomentScattering` returns `None` spaces (short-circuits the guard + defeats
pyright `V`-unification). FIX (orthogonal to the type split, additive leaf): give Λ
real `domain=codomain=frame.basis_space` (it IS a `W→W` coeff endomorphism). Deletes
the cast, re-arms the guard, lets pyright unify. The 2-param `LinearOperator[Din,Cout]`
(#226, never landed) is a SEPARATE larger refactor conjugate does NOT need.

## Q4 (test-side generality) — Frame(basis,measure,*,test:Basis|None=None); discipline is a PROPERTY

Saad `(K,L)` = (trial,test). `Frame(basis, measure)` collapses `L` to "trial sampled
by measure" ⇒ Galerkin-only. Minimal generalization: optional `test: Basis | None`,
default `None ⇒ test=basis` (Galerkin, ZERO-churn default). The test side is a
**Basis** (a distinct PG test space IS different test functions; analysis becomes the
cross-Gram `T_test* W T_trial`). Honoring the rulings:
- RETIRE `GalerkinProjection`+`PetrovGalerkinProjection` (`projection.py:149-211`,
  ZERO subclasses, #268). Discipline = `frame.is_galerkin ≡ (test is basis)`, derived.
- defer-until-2: ZERO genuine-PG, ZERO LS consumers (the homogenization CORRECTION
  proved homogenize/condense are Galerkin-in-reweighted-metric, NOT PG). Build the
  SEAM (the `test` field + `is_galerkin`), NOT the PG cross-Gram path — `analyze`
  keeps the single-basis fast path when `test is basis`, RAISES NotImplementedError
  for `test is not basis` (guarded stub).
- `is_galerkin` COMPUTABLE from `(basis, test)`, no class tag.
- canonical_measure: do NOT add it to Basis. `is_galerkin=(test is basis)` (test-side
  IDENTITY) decides discipline WITHOUT any basis naming a self-dual measure (which
  neither concrete basis does). The `measure==basis.canonical_measure` reading needs
  a basis method that does not exist + the homogenization correction put the weighting
  in the MULTIPLIER not a non-canonical measure. Test-side identity is the clean
  criterion.

## Q5 (least-squares) — flag-and-defer; the expressivity boundary is REAL

LS = `test = A·K` (Saad §5.2.2, `L=AK`, operator-image test space). Does NOT fit the
`test: Basis` seam: (i) needs the operator `A` the `(basis,measure,test)` triple does
not own; (ii) DECISIVELY the LS Gram `(AK)*W(AK)` is NON-diagonal (operator-image test
functions overlap) ⇒ `G_LS⁻¹` is a real linear SOLVE, not `apply_inverse_metric`'s
diagonal reciprocal ⇒ a `CAP_SOLVE` coefficient space. Keep `test: Basis | None`
(Galerkin/PG), NOT `Basis | Operator`. Flip condition: a consumer needing `test=A·trial`
+ a non-diagonal-Gram CAP_SOLVE coeff space. Discriminator test: `Frame.from_least_squares`
is structurally impossible from `(basis,measure)` alone (must supply `A`).

## Augmentations (b) flux-as-multiplier durable, (c) condensation same shape

(b) FLUX-AS-MULTIPLIER is the durable principle; RETIRES `ScalarFlux.flux_volume_measure()`
from the machinery's CONSUMPTION path. Live `homogenize` (`solution.py:419-421`) ALREADY
binds the GEOMETRIC `volume_measure` (1-D, group-indep) + flux as `phi*channel` (`M_φ`).
A per-group `flux_volume_measure` RAISES at `measure.py:207` (1-D constraint). So
"solution emits a measure" is BLOCKED for the per-group case; "solution emits a
MULTIPLIER, discretization emits the MEASURE" is what production runs and is the durable
principle. `DiscreteMeasure.weights` stays 1-D (blast radius zero). The keystone seam
survives only for the group-INDEPENDENT case, but the verb takes multiplier uniformly
(ONE path).

(c) CONDENSATION = SAME `project_weighted(Σ_g, weight=φ_g)` verb, K = coarse-group
indicator on EnergyGrid, measure = energy-bin widths/lethargy, spectrum φ_g = w. ONE
structural flag: energy K is keyed by GROUP boundaries (material-decoupled,
`searchsorted` on coarse group edges) ⇒ coeff space `EnergyGroupSpace` not mesh-bound
`RegionSpace`; output = portable `dict[int,Mixture]` not `MaterialMesh` (the K-carrier
mesh-coupling asymmetry, homogenization memo aug-b). NEW code = `EnergyGrid.indicator_basis(coarse_groups)`
accessor (greenfield leaf, SAME `IndicatorBasis` class with `edges_per_axis=(group_edges,)`),
symmetric with `Mesh.indicator_basis()`. No frame change.

## Smell #16 sightings (adds to the L-003 forensic catalog)

- Shape 1 (paths to one operator): the weighted-Galerkin projection is inlined 3× in
  `homogenize` (`solution.py:422/441/455`) + structurally a 4th (scattering kernel's
  `R∘Λ∘M` and homogenize's `G⁻¹∘M` are two un-shared "frame-analysis post-composed with
  a coeff operator" assemblies). FIX: the `project_weighted`/`conjugate` verbs.
- Shape 2 (one quantity, two homes bridged by hand): the weighted Gram `G_w` lives as
  BOTH `basis.mass_matrix(measure)` AND a hand `dataclasses.replace(basis_space,
  inner_product_weights=region_flux_integral)` (`solution.py:423`). The bridge IS the
  weighted frame operator un-named; `project_weighted` owns it.
- Shape 4 (third path about to be written): condensation is greenfield → a 4th inline
  copy unless the verb is extracted NOW.

## Refuted frames (durable UNEXPLORED)

- Category theory / Kan extension — (K,L) pair + partition_by/pushforward + the frame
  pair fully capture it; no discriminating test. (L-001.)
- Tensor networks/MPO — dense rank-full tables, no bond-dimension knob. (L-001.)
- Homology — no ∂²=0; R∘M is the band-limited projector (idempotent), not a boundary
  map. (L-001.)
- Group/rep theory — SO(3) on the ANGULAR frame only (spatial/energy have none,
  `measure.py:289`); a 3-axis group unification is structurally void. Per-ℓ
  block-diagonality already named. Fires only for block-diagonalizing the angular Gram.
- Spectral theory — `S=T*T` diagonal (4π·I / diag(m_R)); the diagonal-Gram fact IS the
  structurally-simpler payoff (project_weighted is a reciprocal not a solve), not a
  spectral lever. The LS non-diagonal Gram is exactly where this boundary is crossed.
- Radon-Nikodym — `dμ_φV/dμ_V=φ` IS `M_φ`; naming it R-N adds vocabulary not a test.
  `pushforward` explicitly drops the R-N Jacobian (`measure.py:558`). (L-001.)
- Lattice/order theory — methods are POINTS in (K,L,w) space, not a projection lattice.

## First tests (each DISCRIMINATES — L-002)

- project_weighted (the deep verb): `project_weighted(f, weight=w)` == `M(w⊙f)/M(w)`
  bit-identical (`array_equal`) to the current inlined `homogenize` body, for a
  spatially-VARYING w (constant-per-region w cannot fail → rejected). Negative companion:
  a `DiscreteMeasure((n_fine,ng) weights)` RAISES (`measure.py:207`) — proves multiplier
  not batched-measure.
- conjugate: `conjugate(Λ).apply(ψ)` `array_equal` (0 ULP) to current
  `OperatorProduct(R,OperatorProduct(Λ,M)).apply(ψ)` AND `reconstruct_after(Λ).apply(M·ψ)`
  == `conjugate(Λ).apply(ψ)` (0 ULP, sub-op factors the round-trip). A re-tabulating
  windowed arm drifts ULP → REDs.
- test-seam: `Frame(SH,m)` and `Frame(SH,m,test=SH)` both `is_galerkin==True` +
  `array_equal` analysis; `Frame(SH,m,test=Other)` `is_galerkin==False` AND analysis
  RAISES (seam reserved, not silently wrong). A class-based discipline CANNOT produce
  the first (same `_FrameAnalysis`, no class to read).
- LS boundary: `Frame.from_least_squares(basis,measure)` structurally impossible without
  `A` (the discriminator — drops the `L=A*K` push if it claims otherwise).

## Elegance assessment

- structure-exposing (strong): the weighted frame operator `S_w=T*M_wT` is named as the
  ONE object the `/Φ_R`, SH `2ℓ+1`, Indicator `1/m_R` all are; conjugate exposes `R∘A∘M`
  as similarity transport; Galerkin exposed as `test is basis`.
- structurally-simpler (strong): deletes the 3×-inlined Gram-inverse divide + the hand
  metric-injection (Smell #16 shape 2) + the 2 dead PG ABCs + the `cast`.
- expressive (strong): one `project_weighted` covers homogenize+condense+plain SH proj;
  one `conjugate` covers kernel+windowed arm; one `Frame` signature spans Galerkin(default)
  +PG-seam.
- algorithmic-advantage: group/spectrum axis rides the trailing tensor axis through ONE
  frame (no ng-loop); diagonal Gram keeps project_weighted a reciprocal.
