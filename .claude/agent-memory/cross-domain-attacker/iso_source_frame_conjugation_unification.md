---
name: iso-source-frame-conjugation-unification
description: VERDICT on the proposed "every isotropic source (P0 scatter, n2n, fission) is iso_frame.conjugate(energy_op)=R₀∘K∘M₀ in EVERY model (SN/CP/MoC/diffusion)" unification. THREE-PART RULING. (1) The frame-conjugation HOME for the ENERGY operator K on φ-space is correct + native AND ALREADY REALIZED in SN — `ScatteringOperator.full_scatter_kernel` (scattering.py:850) IS `frame.conjugate(Λ_{ℓ≥0}+N2n)` through the HARMONIC frame (rank-(L+1)², skip_l0=False), iso scatter "made nice like aniso" with free transpose. (2) The proposed FACTORIZATION (a SEPARATE rank-1 `iso_frame` × K) is WRONG-OBJECT: the "isotropic angular frame whose coefficient is the scalar flux" is NOT a new ConstantBasis — it is `quadrature.angular_frame(0)` = GalerkinFrame(SphericalHarmonicBasis(L=0), S²measure), the ℓ=0/trivial-SO(3)-irrep V₀ sub-block of the EXISTING harmonic frame. M₀=∫=the ℓ=0 analysis row, R₀=broadcast=addition-theorem at ℓ=0 (P₀=1). Minting ConstantBasis forks R∘M (Smell #16 shape-2). For diffusion/CP/MoC M₀=R₀=Id ⟹ the iso_frame is DEGENERATE (L-004 property-vs-type: zero applied non-id morphism). (3) The dyad outer(χ,⟨νΣf|) as the K=Id degenerate of rank1_frame.conjugate(Id): factor-count-true but a wrapper with R∘Id∘M=R∘M vacuous middle — REJECTED (keep `outer`; it IS the single-mode normal form WITH its own dual-dyad transpose operator.py:1831). The dyad↔multi-mode-frame relation is DEGENERATE-CASE (a frame manages a STACK of dyads), correctly carried by the shared IntegralKernelOperator Protocol both satisfy, NOT by routing the dyad through a frame. THE REAL WIN (model-independence claim is TRUE + pays): collapse the THREE hand-rolled cross-model fission transcriptions — diffusion solver.py:260 `chi*fission_rate[:,None]/keff`, CP solver.py:496 identical, SN fission.py:308 typed dyad — onto the ONE relocated FissionOperator (#261 already moved it to transport/operators/; diffusion/CP consume the bare-ndarray arm fission.py:530 since their M₀=R₀=Id). Smell #16 shape-1 ACROSS MODELS. Discriminating test: diffusion sums axis=1, the dyad sums axis=0 (functional.py:252) ⟹ array_equal REDs unless the relocation maps the axis convention. CROSS-METHOD DISANALOGY: spatial-homogenization/energy-condensation are NOT frame.conjugate(model-indep-op) — they are PETROV-GALERKIN G⁻¹M `frame.project` (analysis-only EXTRACTION, test≠trial, owned by no operator/no eigenbasis L-009/L-010), a DIFFERENT frame VERB. iso-source=conjugate (Galerkin, SH=scattering eigenbasis); homogenization=project (PG, solution-weighted). Read on main, clean tree (L-005 N/A).
metadata:
  type: project
---

# Iso-source frame-conjugation unification — the three-part verdict

DESIGN VERDICT on a proposed unification: "every isotropic source in transport
(P0 scattering, (n,2n), fission) is, in EVERY discretization, a frame conjugation
`iso_source = iso_frame.conjugate(K) = R₀∘K∘M₀`, with K=model-independent energy
operator on φ-space and iso_frame=rank-1 (ℓ=0/DC) model-specific frame". Read on
`main`, clean tree — `frame.py` / `scattering.py` / `fission.py` / `operator.py` /
`functional.py` / `directional.py` / diffusion+CP `solver.py` file:line is GROUND
TRUTH (L-005 N/A). This is the SECOND iso/angular-frame ownership memo (the first,
`harmonic_frame_ownership_funk_hecke`, answered WHO owns the rank-(L+1)² aniso
frame; this answers whether the rank-1 ℓ=0 source is its own frame and whether the
unification crosses models).

## Part 1 — the conjugation HOME is correct + native AND already realized in SN

The claim "isotropic source = R₀∘K∘M₀" is structurally sound FOR THE ENERGY
OPERATOR K on scalar-flux space, and SN already does exactly this — but through
the HARMONIC frame, not a separate rank-1 iso frame:

- `ScatteringOperator.full_scatter_kernel` (`scattering.py:850`) =
  `frame.conjugate(Λ_{ℓ≥0} + N2n)` with `skip_l0=False`. The FULL P0+aniso+(n,2n)
  in-scatter is ONE frame-conjugated operator `R∘(Λ+N2n)∘M` (one analysis, one
  reconstruction), with its transpose free from `OperatorProduct.apply_transpose`.
  Its own docstring: "the iso path made nice like aniso". So the iso scattering
  source IS already a frame conjugation in SN.
- The rank reduction the proposal worries about (naive `P₀⊗Σ_s0` runs the G×G
  matmul N× per ordinate) is ALREADY avoided: `Λ` (`LegendreMomentScattering`,
  `scattering.py:219`, the `Σ_ℓ P_ℓ⊗Σ_{s,ℓ}` SoTP) acts in MOMENT space, and the
  frame conjugation runs M ONCE to reduce to moments. The conjugation captures the
  reduction; the bare TensorProduct loses it. The proposal's rank-reduction
  instinct is correct and the codebase already honors it.
- SN ALSO keeps `add_iso_source` (`scattering.py:970`, the P0 reaction-rate fast
  path Σ_s0^T@φ) at the `apply` boundary — the perf realization (L18 `1/W`), the
  "fuller-view oracle vs optimized" exception, NOT debt.

## Part 2 — the FACTORIZATION mis-places the rank reduction (wrong object)

The proposal's separate `iso_frame` is the WRONG object to mint:

- The "isotropic/constant Basis whose coefficient IS the scalar flux" is NOT a new
  `ConstantBasis`. It is `SphericalHarmonicBasis(L=0)`, and the iso frame is
  `quadrature.angular_frame(0)` (`directional.py:417`) =
  `GalerkinFrame(SphericalHarmonicBasis(0), S²-measure)`. M₀=`frame.analysis` at
  L=0 = `Σ_n w_n ψ_n` (the ℓ=0 row of the SH analysis); R₀=`frame.reconstruction`
  at L=0 = broadcast (addition theorem at ℓ=0, P₀=1). The constant basis Y₀⁰=1 IS
  the ℓ=0 spherical harmonic. So it is `SphericalHarmonicBasis(0)`, NOT a distinct
  DC basis — and the rank-1 iso frame is a SUB-FRAME (the V₀/trivial-SO(3)-irrep
  block) of the rank-(L+1)² aniso frame, not a sibling.
- Minting `ConstantBasis` forks `R∘M` into two implementations of the same V₀
  projector (Smell #16 shape-2: one quantity, two representations) — and they WILL
  drift at ULP (FP non-associativity, the exact reason `basis/base.py:56` keeps
  fused contractions; the scattering kernel is pinned at 0 ULP).
- For diffusion/CP/MoC, M₀=R₀=Id ⟹ the iso_frame is DEGENERATE in 2 of the 3
  models the unification is supposed to serve. L-004: a frame wrapper is earned
  only when a non-identity morphism is APPLIED in the middle; M₀=R₀=Id is
  type-theatrics. The "constant basis over ANY measure (angular/spatial/energy)"
  generalization is the tell — over a spatial/energy measure the iso_frame is just
  Id (there is no DC projection to do; φ IS the unknown).

VERDICT: the iso_frame is NOT a new object. Its only non-degenerate instance is
the ℓ=0 harmonic sub-frame (SN), which is already inside `full_scatter_kernel`.
The unification's TRUE invariant is K (the energy operator), which is genuinely
model-independent — that is the half that holds and pays (Part 3).

## Part 3 — the dyad is NOT routed through a frame (property-vs-type, L-004)

`outer(χ, ReactionRateFunctional(νΣf))` (`operator.py:1857`, `fission.py:308`) is
ALREADY the single-mode normal form `|χ⟩⟨νΣf|`. Wrapping it as
`rank1_pg_frame.conjugate(IdentityOperator())` (test=⟨νΣf|, trial=χ) has
`R∘Id∘M = R∘M` — a vacuous middle. Zero applied non-id morphisms ⇒ ceremony
(L-004). The dyad↔multi-mode-frame relation IS real but is a DEGENERATE-CASE
relation ("F is what `frame.conjugate(Λ)` collapses to at L=0, Λ rank-1"; a frame
manages a STACK of dyads — `operator.py:1730`, `fission.py:343` already say this),
correctly expressed by the `IntegralKernelOperator` Protocol BOTH satisfy, NOT by
routing the dyad through a frame. The dyad already carries its own transpose (the
dual dyad `|νΣf⟩⟨χ|`, `operator.py:1831`, χ↔νΣf swap, campaign #276) with no frame.
KEEP `outer`; do not wrap. The negative test that proves the wrapper is empty:
`rank1_frame.conjugate(Id).apply(φ)` is `array_equal` to `outer(...).apply(φ)` — a
test that CANNOT red ⟹ the wrapper adds nothing.

## THE REAL WIN — cross-model operator relocation (Smell #16 shape-1 across models)

The model-independence claim ("K_iso/F_energy model-independent; iso_frame
model-specific") is STRUCTURALLY TRUE and is the load-bearing content — but it is
realized by RELOCATING THE OPERATOR, not by minting an iso_frame:

- diffusion `compute_fission_source` (`solver.py:260`) = `self.chi *
  fission_rate[:,None] / keff`; CP `compute_fission_source` (`solver.py:496`) =
  `self.xs.chi * fission_rate[:,None] / keff` — BYTE-IDENTICAL hand-rolled dyads,
  no operator, no frame. SN `FissionOperator` (`fission.py`) is the only typed one.
  THREE procedural transcriptions of `|χ⟩⟨νΣf|` (Smell #16 shape-1 ACROSS MODELS).
- #261 already moved `FissionOperator`/`ScatteringOperator` to `transport/operators/`
  (model-neutral). diffusion/CP consume the bare-ndarray arm (`fission.py:530`,
  built precisely "for KEigenvalue/depletion/diffusion outer-iteration consumers
  that feed bare arrays") because their M₀=R₀=Id. The three transcriptions collapse
  onto the ONE dyad; the correctness claim "all models compute the same fission
  source" becomes a theorem (one `apply`) not a coincidence. Bonus: diffusion/CP
  inherit the free adjoint dyad for their importance/k-eigenvalue solves (they have
  none today).
- DISCRIMINATING test: diffusion sums `axis=1` (`solver.py` `fission_rate =
  np.sum(sig_p*phi, axis=1)`), the dyad sums `axis=0` (`functional.py:252`,
  keepdims). `array_equal` REDs unless the relocation maps diffusion's `(nx,ng)`
  layout to the operator's group-leading convention — that axis mismatch is the
  concrete bug the relocation must resolve, and the test that catches a naive
  "just call the operator".

## CROSS-METHOD DISANALOGY — iso-source ≠ homogenization/condensation frame

The proposal asks whether the iso-source conjugation unifies with the spatial-
homogenization / energy-condensation frames already in the codebase. STRUCTURAL
DISANALOGY (L-009/L-010), and naming it is the deliverable:

- iso-source = `frame.conjugate(K)` = `R∘K∘M` (`frame.py:207`), a GALERKIN
  conjugation: the SH frame IS the scattering operator's eigenbasis (Funk–Hecke,
  test=trial), three faces, reconstructs a per-ordinate SOURCE.
- homogenization/condensation = `frame.project(Σ)` = `G⁻¹M·field` (`frame.py:311`),
  a PETROV-GALERKIN coefficient EXTRACTION: analysis-only (NO reconstruction
  conjugation), test≠trial (solution-weighted `φ·1_R`), owned by NO operator (no
  symmetry/eigenbasis), and the operator being coarsened (Σ) is the DATA, not a
  conjugated middle operator.
- These invoke DIFFERENT frame VERBS (`conjugate` vs `project`) and DO NOT unify.
  Routing homogenization through `conjugate` mistypes it (it has no
  reconstruction-conjugation — it extracts coarse coefficients, does not
  reconstruct a fine source); routing the iso-source through `project` mistypes it
  (it has no G⁻¹ normalization — it reconstructs a source, does not average). So:
  K_iso/F_energy ARE `frame.conjugate(model-indep-op)`; the homogenization/
  condensation operators are emphatically NOT — they are `frame.project`.

## Refuted frames (durable UNEXPLORED for this problem class)

- **Fiber bundle over the ℓ=0 mode** — V₀ is the 1-dim trivial SO(3) irrep ⟹
  trivial line bundle, R₀∘M₀ = rank-1 projector, no non-trivial fiber transition.
  The fiber-bundle frame fires for the flux→source TRANSITION (the `1/cm` gain,
  `flux_to_sourcesink_operator_contract`), NOT the angular ℓ=0 embedding. (L-001.)
- **Coalgebra / comonad** — the Cofree comonad fires for the TIMED axis
  (`issue_257_carrier_typing`), orthogonal to the energy-operator conjugation. No
  co-unit/co-multiplication in the iso source. (L-001.)
- **`Σ_ℓ P_ℓ⊗Σ_{s,ℓ}` SoTP with P₀ the iso projector** — this IS
  `LegendreMomentScattering` (`scattering.py:219`), already in the codebase as Λ,
  already conjugated (not bare-TP'd) to avoid the N× rank preservation. REALIZED,
  not unexplored. (Confirms, does not refactor.)
- **Category theory / adjunction (M₀⊣R₀)** — no lever beyond the frame
  analysis/synthesis pair + the named Plancherel metric g_C (the ℓ=0 degenerate
  g_C=4π). Vocabulary, not a test. (L-001.)

## Elegance assessment

- **Structure-exposing (strong):** names the rank-1 iso frame as the V₀/trivial-
  SO(3)-irrep SUB-BLOCK of the existing harmonic frame (`angular_frame(0)`), not a
  sibling ConstantBasis; names the three cross-model fission transcriptions as ONE
  un-factored dyad (Smell #16 shape-1 across models); names the
  homogenization/iso-source split as conjugate-vs-project (two frame verbs).
- **Confirms in part, refactors in part:** Part 1 VALIDATES SN's existing
  `full_scatter_kernel`; Parts 2-3 REFUTE the new-object proposals (ConstantBasis,
  dyad-through-frame); the REAL WIN is the cross-model operator relocation
  (diffusion/CP onto FissionOperator) — a concrete refactor with a fail-able
  axis-convention test.
- **The deliverable** is the THEOREM (iso_frame = angular_frame(0) = V₀ sub-block;
  the dyad is already normal form; K is the model-independent invariant) + the
  cross-model relocation + the conjugate-vs-project disanalogy — not a new frame.
