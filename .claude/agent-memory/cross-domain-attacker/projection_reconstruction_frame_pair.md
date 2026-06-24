---
name: projection-reconstruction-frame-pair
description: Frame attack on the M/R asymmetry in numerics/projection.py — M (MomentProjection) is generic-by-codomain, R (HarmonicMomentReconstruction) is specific-by-name+carrier. VERDICT the asymmetry is a half-applied P1.3 refactor, not structural; the cleanest foreign frame is the analysis/synthesis pair of a 4π-TIGHT FRAME on S² (one SphericalHarmonicFrame owning Y,W,L + the space-pair); the three weight families (W, g_C, (2ℓ+1)) are ONE convention datum read three ways with (2ℓ+1)=4π·g_C⁻¹ DERIVED. Symmetric fix: rename R→MomentReconstruction, give it .domain/.codomain, derive (2ℓ+1) from R.domain.addition_theorem_factor. Read on main (clean, not a branch — ground truth).
metadata:
  type: project
---

# Projection/Reconstruction (M, R) frame pair — `numerics/projection.py`

DESIGN VERDICT (detection feeding an architecture call). Read on `main`,
clean tree, file:line is GROUND TRUTH (no Nexus/worktree staleness — L-005 N/A).
Builds DIRECTLY on `spatial_order_type_vs_property_criterion.md` (the angular
M/R dual already PASSED the property-vs-type criterion → two types, correct).
This attack is one level DOWN: the asymmetry *within* the confirmed pair.

## The three structural verdicts (answering the main agent's three Qs)

1. **The M-generic / R-specific asymmetry is a HALF-APPLIED refactor, not
   structural.** P1.3 de-specialised M's NAME (`Harmonic` dropped, harmonic-ness
   moved to the typed codomain, `projection.py:263-269`) but never gave R the
   symmetric treatment. The symmetric form is CORRECT and decidable: M's output
   type is fixed by `.codomain`; R's INPUT type should be fixed by `.domain`
   (the SAME `SphericalHarmonicSpace`). Fix = rename
   `HarmonicMomentReconstruction → MomentReconstruction`, give it
   `.domain`/`.codomain`, derive `(2ℓ+1)` from
   `R.domain.addition_theorem_factor`. The kept `Harmonic` name + the
   hand-carried `two_l_plus_one` field (`projection.py:551`) are the two
   symptoms of "R was never given its spaces."

2. **Cleanest foreign frame = the (analysis, synthesis) pair of a FRAME on S²**
   (harmonic analysis / frame theory), realised as ONE
   `SphericalHarmonicFrame`/`HarmonicMomentTransform` owning `(Y, W, L)` + the
   space-pair, exposing `.project()` (=M) / `.reconstruct()` (=R) /
   `.synthesize()` (=naked `S0`) / `.H` (=`Π*=g_C·S0`). Analysis op `T`=M,
   synthesis op `T*`=`S0`, frame operator `S=T*T=4π·I` on the band-limited
   subspace (this IS `Π R = 4π I`, `projection.py:468`). It is NOT a
   Bijection/Isomorphism — `R∘M = band-limited projector ≠ I`
   (`projection.py:161`) is REAL rank-deficiency; the precise object is a
   Galerkin transform / `4π`-tight frame. The discriminating test that kills the
   "Isomorphism" reading: `R∘M ψ != ψ` for ψ with content above order L.

3. **All three weight families belong to the ONE frame object (the space-pair),
   `(2ℓ+1)` is DERIVED not data.** `W` = analysis measure (fine side, on
   M.domain `:404-422`); `g_C=diag(4π/(2ℓ+1))` = frame-operator metric (coarse
   side, on M.codomain `:380-401`, lives on `SphericalHarmonicSpace`); `(2ℓ+1)`
   = the canonical-DUAL-frame synthesis coefficient `= 4π·g_C⁻¹`. The reciprocal
   `g_C·(2ℓ+1)=4π` (`projection.py:518`) is pinned ONLY in a docstring + a
   runtime test — no object MAKES it hold. Single source of truth already exists:
   `SphericalHarmonicSpace` carries both `metric_per_ell` (g_C) AND
   `addition_theorem_factor` ((2ℓ+1)) delegated to the basis
   (`spherical_harmonic_space.py:209-226`). Both M and R should read ALL weights
   from the shared space-pair.

## Smell #16 sightings (shapes 1 + 2) — adds to the L-003 forensic catalog

- **Shape 2 (one quantity in two incompatible reps bridged by hand) — DOMINANT.**
  The SH convention datum lives in THREE storage homes (W, `(2ℓ+1)`, `g_C`)
  linked only by the prose identity `g_C·(2ℓ+1)=4π`. `(2ℓ+1)` is `4π·g_C⁻¹`
  un-derived. NEW WRINKLE vs prior shape-2 sightings (#168 sweep/apply, the
  σ_t-positional-thread): here the two reps are a quantity and its own RECIPROCAL
  on the metric (g_C vs (2ℓ+1)=4π·g_C⁻¹) — a metric and its inverse stored as
  independent fields. The un-named operator is the canonical-dual-frame map.
- **Shape 1 (paths to one operator).** The naked synthesis `S0(c)_n =
  Σ_ℓm Y_ℓm(Ω_n) c_ℓm` is realised THREE times as three differently-weighted
  einsums: `Π⊤=w_n·S0` (`:482`), `R=(2ℓ+1)·S0` (`:597`), `Π*=g_C·S0` (via
  `_AdjointOperator` `operator.py:664-669`). `S0` has no name; the fix names it
  once and post-multiplies by one of the frame's own metrics.

## Metric-blindness: R opts out of the adjoint-for-free machinery

R carries NO `.domain`/`.codomain`, so the metric-aware `_AdjointOperator`
(`operator.py:664-669`) that gives M its `.H` for free CANNOT read R's spaces →
`R.H` is uncomputable generically today. Giving R its spaces (verdict 1) makes
`R*` free — the exact map the anticipated MC adjoint-moment consumer needs
(`projection.py:40-41`). This is the `coefficient_field_promotion_frames.md`
Frame-3 metric/role-blindness smell, recurring on R.

## The unify-NOW gate (cross-method): the §10 PN second instance

Do NOT build the `SphericalHarmonicFrame` object on ONE instance — that violates
unify-after-two-instances. The trigger that earns it: the §10 PN moment-space
projection/reconstruction sibling (promised `projection.py:267-269`). Gate:
- If PN's reconstruction ALSO has a `synthesis_metric = 4π·g_C_PN⁻¹` of its own
  Gram → the FULL frame object generalises, build it.
- If PN has NO addition theorem (PN moments are not collocation) → only the
  abstract (analysis, synthesis `S0`, frame-operator `g_C`, `Π*`) TRIPLE lifts,
  and `(2ℓ+1)` stays an SH-specific leaf — which would VINDICATE keeping a
  `Harmonic`-qualified reconstruction (the one scenario where verdict-1's rename
  is wrong). Sketch the PN three-weight families before unifying.

## Refuted frames (the durable UNEXPLORED for this problem class)

- **Category-theory adjunction (M⊣R)** — duality fully captured by the frame
  analysis/synthesis pair with named diagonal metrics; no discriminating test
  the frame frame lacks. Vocabulary, not a lever. (L-001 pattern.)
- **Spectral theory** — frame operator `S=4π·I` is a trivial single-eigenvalue
  spectrum (g_C the per-ℓ refinement, already named). No non-degenerate spectrum.
- **Tensor networks/MPO** — Y is dense rank-full, no bond-dimension knob.
- **Group/rep theory (Wigner-D/isotypic)** — SO(3) equivariance IS present and
  IS the `Π R=4π I` reproducing identity, but for the one-object + weight-home
  questions it is already fully used; per-ℓ block-diagonality of g_C/(2ℓ+1) is
  the only isotypic fact and the frame frame names it. Would fire only if
  rotating the quadrature / block-diagonalising by ℓ were the task.
- **Differential geometry** — no curvature/redistribution/`1/r` term; algebraic
  contraction on fixed S². (L-001/L-008 pattern.)

## First tests (each DISCRIMINATES — L-002)

- Verdict-3 (derived (2ℓ+1)): build the frame from `(Y,W,L)` ALONE, assert
  `R.two_l_plus_one == 4π/codomain.g_C` bit-identical (`array_equal`, 0 ULP) AND
  reconstruction matches. A "store (2ℓ+1) as independent field" impl FAILS on a
  re-prefactored (`4π/(2ℓ+1)`) SH basis — it carries a stale literal.
- Verdict-2 (frame not isomorphism): `R∘M ψ != ψ` for ψ above order L (it equals
  the band-limited projection). A "Bijection" object wrongly claims `R∘M=I`;
  this RED-s it. (Asserting only `M∘R=4π I` cannot fail → rejected.)
- Verdict-1 (R reads its metric): `MomentReconstruction.from_space(g_C_prefactored)`
  yields a DIFFERENT (2ℓ+1) than on the canonical basis — proves it reads
  `.domain`, not a baked literal. NEGATIVE-flavoured (the illegal "same literal
  on every basis" is what fails).

## POST-LANDING (2026-06-23): the abstraction landed; naming-seam verdicts

The frame DID land: `Frame(basis, measure)`, `Basis` ABC,
`SphericalHarmonicBasis(L=L)`, faces `frame.analysis`(M=T)/`frame.reconstruction`(R),
`Quadrature.angular_frame(L)`. Three naming-recon verdicts on the surviving seams
(durable — these are name-the-axis-not-the-basis + type-vs-property rulings):

1. **SH truncation `L`: keep the PARAMETER `L`, fix the PROSE `order`→`degree`.**
   Param and concept-name DIVERGE. `L` is the unique highest-signal token in BOTH
   idioms (SH `ℓ≤L`, transport `P_L`) and the carrier `SphericalHarmonicBasis`
   already supplies the `Y_ℓ^m` context that `max_degree` would redundantly
   re-state — so `L` wins (high-signal "let the object name the part"). BUT the
   docstrings say "order L" in ~8 places and that is a LATENT CONVENTION BUG: for
   `Y_ℓ^m`, `ℓ`=**degree**, `m`=**order**; "order L" reads as bounding `m`. The
   basis documents this very convention, so the prose contradicts its own
   contract. Prose/derived-name term = `degree`. `band_limit`/`bandwidth` is the
   correct frame-theory term for the REGIME (`N>(L+1)²` retraction, already at
   `frame.py:33`) — names the frame capability, NOT the basis truncation; do not
   migrate it onto the basis param.

2. **`Quadrature.angular_frame(L)`: keep. Name the AXIS, not the BASIS.** The
   invariant is the measure's physical axis (S² directions = `angular`),
   permanent; the SH basis is contingent. `harmonic_frame` breaks the
   `angular_frame`/`spatial_frame`/`energy_frame` greppable family the instant a
   2nd angular basis arrives. The internal `SphericalHarmonicBasis(L=L)` hardcode
   (`directional.py:447`) is honest TODAY (axis↔basis 1:1); 2nd-basis trigger
   forces a SIGNATURE change (`angular_frame(basis=...)`) but NOT a name change —
   which is precisely why `angular_frame` is right now. Record the tripwire in
   the docstring, do not pre-rename.

3. **Two new typed-field questions, two homes, ONE mint + ONE defer:**
   - **Phase axis (angular/spatial/energy) → ON THE MEASURE, MINT IT.**
     `DiscreteMeasure.space` is ALREADY SPENT on a different meaning (opaque
     measurable-space tag `"S^2"`/`"[-1,1]"`, `measure.py:92`) — the σ-algebra
     label, not the physical role. So the axis needs its own field:
     `phase_axis: Literal["angular","spatial","energy"]`, sibling to the reserved
     `invariance_group` (`measure.py:139`, Issue 3) — both are physical-role of
     the measure, both on the measure. `axis` collides with the existing
     direction-cosine axis-index (0/1/2) usage on `Quadrature`; `phase_axis` is
     collision-free and self-documenting. The field value (`angular`) AGREES with
     the `angular_frame` method name — closes the Q2 loop.
   - **Galerkin/Petrov-Galerkin = canonical/oblique DUAL → ON THE FRAME, DEFER.**
     This is the `frame.py:27-35` iso-vs-non-iso discipline (`R∘M=I` invertible
     vs band-limiting projector). FAILS type-vs-property: ONE realization today
     (SH canonical-dual non-iso); the iso case is a `CAP_SOLVE` CAPABILITY, not
     instantiated. `R∘M=I` is already DECIDABLE from whether the analysis face
     advertises `CAP_SOLVE` → minting `dual="canonical"` duplicates the flag, no
     morphism reads it = theatrics. Frame-theory term `canonical dual` stays in
     the docstring (`frame.py:7`). Mint-trigger = a 2nd concrete frame with a
     GENUINELY oblique dual (non-self-adjoint M/R pair).

The two status-quo names (`L`, `angular_frame`) were each already optimal. The
META-RULE this attack confirms (promote toward lessons): **a factory/accessor on
a measure-bound object names the measure's PHYSICAL AXIS, not the current basis —
the axis is invariant under basis-family growth, the basis is not.** Pairs with
the type-vs-property gate: a "discipline" already decidable from an existing
capability flag is a PROPERTY (derive it), not a TYPE (mint it).
