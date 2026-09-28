---
name: weight-placement-weld-attack
description: The 2026-09-28 adversarial attack on "move the PG test weight into the measure" — the weld the original ruling feared is REAL (state into a Hilbert METRIC → flux-dependent space identity, an adjoint that forgets the flux shape, ng/ng² measures, unshareable Grams, signed weights foreclosed); the ruling's WRITTEN ground (adjoint needs test≠trial) is false but its structure/state and Radon–Nikodym grounds hold; WeightedIndicatorBasis IS analysis_dV ∘ M_w (array_equal); the overlap table FACTORS as analysis_coarse ∘ section_fine over w(E)dE (4e-16), the two-step defect is M_c(P_m − I)S_f (4e-16), a rank-k mid frame shrinks it 0.39→0.011 without zeroing it; MeasureCoupling is a derived composite. Open before any weight/measure/coupling/GEC brief.
metadata:
  type: project
---

# The PG weight: test basis vs measure vs multiplier (attack 2026-09-28)

Memo `scratch/posing_sequence/frame_hierarchy/weight_attack/structure/memo.md`; probes `probe_weld.py` (A–F),
`probe_gec_factor.py` (G0–G4), outputs beside them. Re-ground every `file:line` before acting.

## The premise, measured
- Collapse morphisms in the 3 production functions: **15**, not 7. 7 frame-shaped with a positive weight; 1 only as a
  3-D weight (T2 per-pair); 2 counting-frame analyses; **5 are ratios of two pushforwards with DIFFERENT densities**
  (T3 mixed fold `material_xs_field.py:393`, ψ† `mixture.py:614`, sink÷ψ† `:623`, χ† `:633`, νΣf×s `:639`) which NO
  single-weight frame spells in either spelling (F2: 0.11 off; F3: 0.75 off; the two-pushforward ratio 1.4e-16).
- Readers of the frames' Hilbert surface at those sites: 0. The only live verb is `project`.
- X4: the verbs memo's P2 "1.6e-16" is the FP re-association of ONE einsum; `PG(w).analysis(f) ==
  Gal(dV).analysis(w f)` is `array_equal` (the same call). P7 counted constructions, not morphisms.

## The original argument, three strands
- (A) "adjoint weighting needs test ≠ trial; no single measure reproduces ⟨φ*,Σφ⟩" — REFUTED for the 7 frame-shaped
  vector channels (pair weight φ*⊙φ is one weight); TRUE of the 5 hand-rolled bilinear morphisms.
- (B) "measure = grid structure; flux = solve state; ng measures; weights 1-D" — CONFIRMED with numbers (below).
- (C) "the solution is the Radon–Nikodym MULTIPLIER, not the measure" — CONFIRMED bit-identically and DERIVED (G4).

## The weld, measured `[M]`
- B1 `DiscreteMeasure(weights (n,ng))` refused → ng measures (ng² for T2); PG 1 frame.
- B2 nominal identity: `Gal(φ₀V).measure_space == Gal(φ₁V).measure_space == μ_V.space` True while metrics differ 0.8
  (metric-blind seam crossed); axis-built (`of_axes(μ.axis())`): 3 distinct spaces — the fine space changes per solve
  and per group under the axis doctrine.
- B3 `.H`: PG `M^H c = φ⊙Tc` (shape in the OPERATOR); Gal(φV) `M^H c = T G⁻¹c` (shape in the METRIC); `.dual()` agree.
  CMFD prolongation `c_R φ_i/Φ_R` = PG `M^H G⁻¹` (2e-16), Gal off 0.66. Today `PG.analysis.H` RAISES (WIB transpose unbuilt).
- C  `discrete_gram` shareable across φ, φ*φ, ρ, p: PG 4/4, Gal 0/4.
- D  signed φ₁: Gal constructs, Gram DENSE, `‖f‖²_μ = −1.86` (not Hilbert); PG fine. Complex: BOTH wrong —
  `CrossGramInverse` casts the probe to float (`frame.py:200`), `project` = `K_*(wf)/Re K_*(w)`.
- E  product-measure frame (cells×groups) exact but G²/G³: 301 ms vs 1.1 ms at G=8 (N=400, K=40).
- B4 coarse metric 1/V_R vs 1/Φ_R: 0 readers (dead in both).

## The factorisation (coordinator's hypothesis, confirmed) `[M]` G
- `T = analysis_c ∘ section_f` over `μ_s = w(E)dE` on the supermesh: vs `overlap_to` 4.3e-16; PoU rows = section's
  mass preservation. `T_fm T_mc − T_fc = M_c(P_m − I)S_f` to 3.9e-16 abs; `P_m` idempotent, μ_s-self-adjoint;
  exact ⟺ `1_G ∈ span{1_M}` (coarse edges ⊂ mid edges). Rank-k mid frame: defect 0.39→0.011 (k=0..6), σ_G error
  8.4e-2→4.3e-3, non-monotone, never 0 (a step is in no polynomial span); nested 2e-16 at every k.
- `σ_G = (M_c S_f)(σφ)/(M_c S_f)(φ)` = tree's `project` 1.15e-16; fine Gram `W_g = ln(E_hi/E_lo)`; `(M_c S_f)ᵀ1 = 1`
  → the coupling's fine coefficient measure is COUNTING; φ is the coefficient vector pushed.
- **`MeasureCoupling` is DERIVED** (three shipped faces composed); mint it only as the cached table with
  (fine frame, coarse frame, structural measure) provenance — the chain repair is then the composite from its factors.

## Placement verdict
Structure: `T`, `dV`/`w(E)dE`, Grams, coarse identity. State: the two densities `w_num, w_den` — ARGUMENTS of one verb
`condition(T, μ; w_num, w_den, f)` = `M_w` composed before the geometric frame's analysis. Deletes `WeightedIndicatorBasis`
as a Basis (it is `analysis_dV ∘ M_w`), the "one frame per weighting" ruling, and 4 hand-rolled einsums. Refutes both
proposal placements (measure: B; coupling marginal: G4). Keeps the ruling's (B)+(C); retires its written ground (A).

Cross-refs: [[frame-hierarchy-kernel-coupling-frames]] (⛔ its F6 inference superseded here), [[homogenization-measure-derivation-frames]]
(its "flux IS a measure factor, μ=φ·dV" Galerkin reading was the FIRST form of this proposal — same verdict, 2026-06-24),
[[coefficient-field-promotion-frames]] (the multiplier algebra), [[dsa-rp-angular-frame]] (the R/P pair), [[condensation-nonnested-fractional-overlap-frames]].
