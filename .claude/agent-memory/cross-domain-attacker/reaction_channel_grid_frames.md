---
name: reaction-channel-grid-frames
description: The reaction-channel attack on the posing ontology (2026-09-27): a channel-with-role is a CELL of the (reaction × role) grid the data layer already IS; term = role-marginal, channel = reaction-marginal; a cell is ONE object read as direction / rate / sensitivity (Hellmann–Feynman, [M] 1e-10); k = (fission, emission) only; compositions are FIELD directions of the input-layer family, non-affine on Bondarenko data; emission is already cross-system (multiplicity a vector over targets).
metadata:
  type: project
---

# Reaction channels: the grid, its three readings, and where it lives (2026-09-27)

Memo `scratch/posing_sequence/channels_attack_memo.md` (census `channels_census.py`, probe
`channels_probe.py` P1–P5). Re-ground every `file:line` against the live tree before use.

**The verdict.** The linear Boltzmann collision operator is graded by reaction:
`R = Σ_x (M[Σ_x] − ν_x K_x)`, `K_x = Σ_x P_x`, `P_x` Markov. The ruled layer-1 TERMS are
its role-marginals (`C = Σ_x M[Σ_x]`; `S`, `N_2n`, `F` one emission cell each); the
CHANNEL is the reaction-marginal. The grid `(Σ_x, ν_x, P_x)` IS `Mixture` (input layer),
up to the (n,2n) yield sitting in layer 1 (`N2N_MULTIPLICITY`) while ν̄ is folded at the
input layer. The cell → term map is each term's CONSTRUCTION, never a declaration.

**Why it matters (the three readings).** `T_cell = ∂balance/∂a_cell`; `rate(cell, w) =
⟨w, T_cell φ⟩_G` with `w = 1` the counting rate (removal: `⟨♭Σ_x, φ⟩`; emission:
`−ν_x⟨♭Σ_x, φ⟩` by the row-sum law) and `w = ψ†` the SENSITIVITY `∂λ/∂a_cell`
(Hellmann–Feynman, the conormal). The k denominator is `Σ_x (1 − ν_x) R_x` — "absorption"
is a named cell SET, not a datum. Discriminators `[M]`: k scaling fission's removal too →
off by 1.57; (n,2n) counted once → off by 0.33; prompt-only k with precursors → off by β.

**How to apply.**
- A brief asking WHERE a label belongs: check whether the label GRADES an operator the
  tree already sums. Then it is a coordinate, the map to terms is a construction fact,
  and a proposed "declaration + resolver" is a second source of truth (refuted cheaply).
- "Scaling per (channel, role) plus compositions" mixes scalars with FIELDS: a nuclide
  direction is not in the span of cell scalings (residual 0.87 in 8 groups; 0.94 for a
  region-local change). Compositions are pulled back from the input layer's `N ↦ Mixture`,
  affine iff `n_sig0 == 1` (second difference 3.5e-4 shielded vs 1e-16 flat). The
  second-difference rule has a DATA source, not only parameter × discretisation.
- Multiple systems: emission's codomain is a system; the fission multiplicity is a vector
  over targets summing to ν̄; the precursor "coupling βνΣ_f" is the delayed ARROW of the
  fission channel; the block adjoint needs W's metric (ψ†_W scales by G_V/G_W).
- Probe hazards recorded: σ₀ on the interpolation clip reads a false affine zero; an
  over-complete cell basis reads a trivial span; the HF sign is `−k²⟨ψ†,Tψ⟩/⟨ψ†,Fψ⟩`.
