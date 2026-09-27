---
name: pencil-pseudo-resolvent-standard-form-pair
description: Durable frame (first sighting 2026-09-26, the posing-sequence sides-naming attack) — a pencil's "quotient" at(σ)⁻¹M is a PSEUDO-RESOLVENT (resolvent identity holds; a generator exists iff M injective), so "generator/propagator" never name a side or the quotient; the two spectral maps are the pencil's two standard-form reductions. Open for any eigen-question naming or k/α/noise unification brief.
metadata:
  type: project
---

# The pencil's quotient is a pseudo-resolvent; the spectral maps are the standard-form pair

**Fact.** For `at(σ) = A − σM`, `K(σ) := at(σ)⁻¹M` satisfies
`K(σ) − K(τ) = (σ − τ)K(σ)K(τ)` (Pazy 1983 §1.9, Yosida VIII.4: a pseudo-resolvent),
while the pencil's own resolvent `at(σ)⁻¹` does NOT (`[M]` 2026-09-26 probe,
`scratch/posing_sequence/sides_naming_probe.py`: 5e-16 vs 0.28). A pseudo-resolvent is
the resolvent of a genuine operator iff `ker K = ker M = {0}`: the α family (`M = 1/v`)
has the generator `𝒜 = −M⁻¹A` and `R(s;𝒜) = K(−s)`; the k family (`M = F`, rank-deficient)
has NONE — a theorem, not an observation. The two standard forms `A⁻¹M` (spectrum 1/μ,
read by `K_MAP`) and `M⁻¹A` (spectrum μ, read by `ALPHA_MAP`) are what the project's
`SpectralMap` selects between; k is read from the first, α from the second.

**Why it matters for naming.** "Generator" and "propagator" are the time family's
`𝒜` and `e^{t𝒜}` (or, stationary, the resolvent) — neither a side of the pencil nor
`K`; the pair can never name (side, quotient). `K(0)`'s perfect matches are per family:
k → next-generation operator (Diekmann–Heesterbeek–Metz 1990; van den Driessche–Watmough
2002, whose Theorem 2 `s(F−V)<0 ⟺ ρ(FV⁻¹)<1` IS the k/α sign theorem via Varga 3.13);
c → next-collision; α → the generator's potential operator. The sides are the jet
(`at(0)`, `derivative = ∂at/∂σ = −F` for k — the derivative is NEGATIVE; positivity
lives on `K`). The expansion point is the CARRIED-OFF point, which equals the physical
anchor only when the carried term is absent from the physical balance (α), not for k, c.

**Coincidence law (layers 2 and 3).** At the physical point of a coupling family the
pencil pair is a Varga regular splitting of `at(1)`, and one k power step equals one
stationary iterate of `Splitting(implicit=at(0), explicit=−derivative)` before
normalisation (`[M]` array_equal). Same object, different meaning (pole vs rate; free vs
`<1`; forced vs chosen) ⇒ two naming families, one documented meeting point.

**How to apply.** On any "unify k/α/noise/transient" or "name the pencil's parts" brief:
run the resolvent-identity probe first (it discriminates a dropped `M`, a sign drift in
`at`, and a wrong-sign derivative), then classify the family by `TIME ∈ channels`
(spectral, generator exists) vs not (coupling-constant, Birman–Schwinger reading). The
full memo: `scratch/posing_sequence/sides_naming_memo.md`. Related:
[[eigenvalue-posing-layering-frames]], [[power_iteration_vs_keigenvalue_morphism]],
[[feedback-naming-collision-perfect-match]].
