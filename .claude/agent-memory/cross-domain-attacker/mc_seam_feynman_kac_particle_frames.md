---
name: mc-seam-feynman-kac-particle-frames
description: The Monte Carlo seam attack on the three-layer posing architecture (2026-09-27) — stochastic simulation is a layer-3 Strategy over the SAME system iff the system yields its MARKOV VIEW (collision-density factorisation), admits IDENTITY per-axis discretisation factors, and the discretisation is a FRAME (M and R); `sample` on a TERM is the wrong granularity; the MC state is an EMPIRICAL MEASURE in V*; CADIS rides R, CMFD/TRMM ride M with an answer as weight; the only stochastic layer-1 datum is the multiplicity LAW.
metadata:
  type: project
---

# The MC seam: what layer 1 owes a stochastic Strategy (first sighting 2026-09-27)

Memo `scratch/posing_sequence/mc_seam/attack/memo.md` (explorer's ground truth `mc_today.md`
beside it). Re-ground every `file:line` against the live tree before acting.

**Why:** the user asked whether MC is a Strategy over the same system or forces a different
layer-1 terminal object. Verdict: Strategy over the same system; "different object"
`[REFUTED 2026-09-27]` for k, α, fixed source, evolution, adjoint/IFP, noise, sensitivity,
CADIS, CMFD, TRMM and second moments; extinction is the declared nonlinear future kind.

**How to apply (the durable verdicts):**
- **Granularity of `sample`.** The walk is the Neumann series of `K = (Σ_x ν_x Σ_x P_x)·T`
  with `T = (L + C − B)⁻¹` on characteristics; neither factor is a term. A `sample` verb on
  `Term` leaves removal/streaming members unable to honour it and the free flight homeless.
  The right object is a VIEW the system yields (`markov_view(point)`: free-flight kernel over
  the derived total + geometry verbs, collision kernel over the channel set, boundary Markov
  kernel, source measure), in the same relation to the system as the pencil to the question.
  TELL in the tree: two normalisations of one measure (majorant reads stored `SigT`, the
  collision the partial sum — `mc/solver.py:373` vs `:438`).
- **Representation forks `apply`/`sample`.** A Legendre-truncated `P_x` is a signed measure:
  `apply`-able, not sampleable; the positive law is a second representation of the same
  Markov kernel. In CE they coincide. Type `P_x` as the kernel, not the matrix.
- **Root = identity factors.** The per-axis discretisation product admits the identity;
  MC is the point with identity on space/angle and SHARES the energy factor (data-layer
  collapse) with SN. "MC has no discretisation" is imprecise. Do not mint a
  `ContinuousSystem` class no constructor can fill (`[M]` no pointwise energy data in tree).
- **Carriers on an identity axis are the function/measure pair (V, V*)**; the MC state is
  an empirical measure (the bank = the Perron eigenMEASURE `Fψ`; Del Moral's genetic
  particle approximation; O(1/N) bias Brissenden & Garlick 1986). `Outcome.state` stays;
  the seam is `read(functional)` + an additive `Estimated(value, se, n_eff)` Evidence
  variant. NEVER retype the deterministic outcome as "functionals over a complete frame":
  derived questions consume ψ, ψ† as vectors.
- **Coordinates act on the path measure by ROLE**: emission = per-event multiplier
  `ν_x/k`; removal (incl. `s/v`, complex `iω/v`) = Feynman–Kac POTENTIAL
  `exp(−δΣ·ℓ)` × collision ratios. The α<0 "pseudo-production" branch is a Strategy
  choice (α-weighting is sign-uniform). Score function on path space = Hellmann–Feynman
  conormal: two quadratures of one derivative.
- **Hybrids ride the frame's two legs.** CADIS = R (discrete ψ† → continuous importance,
  a SEED; unbiasedness = "a seed never changes the answer"). CMFD-on-MC and Betzler 2018
  TRMM = M with the ANSWER as weight → a coarse layer-1 SYSTEM (the ruled "condensation →
  new Problem"), whose identity must digest the answer. The `Discretization` Protocol
  `[M]` had neither leg (`transport/method.py:150`); `Frame` has both (`numerics/frame.py:514,519`).
- **Branching questions.** Mean-field system = linearisation of Pál–Bell at g=1 ⇒ adjoint
  = importance, IFP = Krein–Rutman left vector read through forward progeny. Feynman-Y /
  Rossi-α = fixed-source on the adjoint with a quadratic source (ruled "SOURCES"); the one
  missing datum is `ν` as a LAW with `.mean` (additive; tree stores a float).
- **Capability tables invert**: MC offers `Exact` evolution and refuses `Stepped`.

**Candidate skill growth (one sighting, not proposed):** Part A.4 Feynman–Kac row gains the
lever "interacting particle system ⇒ the bank is the eigenmeasure; population control =
selection; certificate = bias + autocorrelation"; Part B "MC borrows from sensitivity"
gains "score function = HF conormal on path space".

Related: [[posing-ontology-clean-attack-frames]], [[equilibrium-carrier-laurent-point-kinetics]],
[[reaction-channel-grid-frames]], [[projection-reconstruction-frame-pair]],
[[alpha-dome-chart-vs-measure-cross-method]].
