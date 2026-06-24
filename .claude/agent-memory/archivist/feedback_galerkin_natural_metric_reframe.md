---
name: galerkin-natural-metric-reframe
description: Doc-recipe for documenting a "flux-weighted average is Galerkin-in-natural-metric, NOT Petrov-Galerkin" carve — the L²(weighted) projection reframe + the sibling-page staleness tombstone.
metadata:
  type: feedback
---

When a refactor reroutes a hand-rolled weighted-average collapse
(cross-section homogenization, energy condensation) through the project's
discrete `Frame` and DERIVES that the weighting is an
**L²(weighted-metric) orthogonal (Galerkin) projection**, the doc job is
TWO coupled edits, and the second is a retraction-tombstone on a sibling
page that the carve just falsified.

**Why:** the homogenization carve (#268, branch
`refactor/operator-inverse-algebra`) proved spatial homogenization is
L²(φV)-Galerkin (test=trial={1_R}, orthogonal in the φV metric), NOT
Petrov-Galerkin. The Petrov-Galerkin reading was an ARTIFACT of insisting
on the unweighted dV metric: ⟨Σ, φ·1_R⟩_dV = ⟨Σ, 1_R⟩_φV is the SAME map
— the flux multiplier lives in the test function OR the measure, and the
`Frame` reads it off the MEASURE → Galerkin. This collapsed the
Galerkin/Petrov type distinction to a derived PROPERTY
(`measure == basis.canonical_measure`).

**How to apply:**

1. **The owning theory section (the carve's home page).** KEEP the
   verifies-target eq-label by name (`sn-homogenization-rate-preservation`,
   a matrix verified-block entry with 2 `verifies` bindings) and its body.
   KEEP all documented-only decomposition labels. ADD a new H2 subsection
   deriving the projection, with NEW documented-only labels each carrying
   `.. vv-status: <label> documented` + a rationale comment naming the
   bit-identity/foundation gate. The load-bearing pieces a fresh reader
   needs, in order: (a) the coarse trial space = span{1_R} (the P0 /
   IndicatorBasis space); (b) the Galerkin normal equations → disjoint
   support → diagonal Gram → the average falls out; (c) **the measure is
   DERIVED not chosen** (dV gives volume-average, breaks rate at a flux
   gradient — pin to the φV-vs-dV discriminator test); (d) Galerkin = the
   two-readings-same-map identity; (e) the Radon–Nikodym split μ_φV=φ·μ_V
   and WHY DiscreteMeasure.weights stays 1-D (a measure = one mass/atom;
   per-group φ would be ng measures → φ is an integrand MULTIPLIER on a
   trailing axis, not a measure); (f) mesh-YIELDS-IndicatorBasis (basis =
   measure-free half of a frame; mesh carries the measure → inheriting
   conflates the two roles); (g) the composite projector R∘G⁻¹∘M with the
   1/Φ_R normalization in the SPACE metric (`apply_inverse_metric`,
   Moore–Penrose → zero-flux region → 0 for free, no guard); (h) the
   Mode-11 sentinel pinning the Frame is actually on the call graph.

2. **The sibling page the carve falsified.** A full theory page built
   around the now-wrong TYPE distinction gets a top-of-page
   `.. warning:: **Superseded framing (date, Issue #N).**` tombstone (L-007):
   (a) what the page claims, (b) one-line why it's wrong, (c) forward-ref
   to the new derivation subsection + the tracking issue. Do NOT rewrite
   the page — the full reframe + ABC retirement is the issue's job. State
   explicitly "prose below is not yet rewritten; read through this
   correction."

3. **The asymmetry-law table gains a frame reading.** condense/homogenize
   return different types because they bind a DIFFERENT trial basis K into
   the one frame: homogenize = GEOMETRIC K (spatial cell indicators →
   mesh-COUPLED → MaterialMesh); condense = SPECTRAL K (group indicators
   under L²(spectrum) → mesh-DECOUPLED → portable dict[int,Mixture]). The
   return-type asymmetry IS the K-axis (space vs energy) of the one
   mechanism showing through.

**Quality self-assessment (this task):** Derivation-depth 5, Cross-refs 5,
Numerical-evidence 4 (the φV-vs-dV discriminator + Mode-11 sentinel are
the evidence; no convergence table because it's an exact projection, not
an inexact method — structurally absent, not a deficit), Failed-approaches
4 (the membership-matmul named as superseded-but-correct + the
Petrov-Galerkin artifact explained), Code-traceability 5, Derivation-source
5 (the IndicatorBasis/Frame/space docstrings ARE the algebra-of-record —
read them first per L-005; they were already richly written).

See [[lessons]] L-004 (vv-status on representational labels), L-007
(retraction tombstone), L-002 (numerics.basis/frame/space are NOT
automodule'd → :class:/:meth: render plain-text by page convention; match
the sibling page, don't half-surface). The φV-vs-dV discriminator is the
load-bearing guard on the whole L²(φV) claim — always cite it.
