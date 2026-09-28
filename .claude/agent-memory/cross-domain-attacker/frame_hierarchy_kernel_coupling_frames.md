---
name: frame-hierarchy-kernel-coupling-frames
description: The frame hierarchy verdict (2026-09-28) — the apex "kernel between two measure spaces" IS LinearOperator with two metric ends (exists); the one missing datum is the coefficient-side measure ν as a FRAME SLOT with three fillers (declared / the Gram = L² pushforward / the marginal M·1 = L¹ pushforward); frame bounds = spec(W_ν^{-1/2} G W_ν^{-1/2}) = ρ on one-hot tables, the canonical dual's projector is ν-FREE; the PoU (OverlapBasis) arm's `project` is a Markov DISINTEGRATION (not idempotent, not associative), so `frame.py:583` "G-orthogonal projector" is false there; the tree's own words `average`/`marginalize`/`project` are the three ν's; OverlapBasis(IndicatorBasis) is an inverted edge; Kernel collides integral-vs-nullspace.
metadata:
  type: project
---

# Frame hierarchy — kernel / frame / coupling (attack 2026-09-28)

Memo `scratch/posing_sequence/frame_hierarchy/structure/memo.md`, probes
`probe_frame_hierarchy.py` (P1–P5), `probe_two_step.py`, `census_handrolled.py`.
Re-ground every `file:line` before acting.

## The hierarchy (general → specialised, edge = the specialising property)
- **K0** `LinearOperator` between two MEASURED spaces = "a kernel between two measure
  spaces" in finite dimension. EXISTS; a `Kernel` type would twin it (and be a 3rd sense).
- **K1** table-presented: `Φ(node, mode)` + `μ_M` + a coefficient measure `ν` as a SLOT.
  Today `FrameBase` has the slot FIXED to one derived choice.
  - **ν declared** (continuous / weighted frame; Ali–Antoine–Gazeau, Balazs): retraction
    `♯_ν∘M`, section `R∘G⁻¹∘♭_ν`. One-hot Φ ⇒ bounds = ρ = d(φ_*μ_M)/dν `[M]` 4.4e-16 on
    4 fixtures (V: A=B=Σw=2; pushforward: 1; unequal fibres `[3,2,1]`; trace: live 1, dead 0).
    Instances: Trace/System restriction (injective), Permutation/deck (bijective),
    AxisRetraction/Section + Lambertian PartialCurrent/IsotropicEmission (surjective).
  - **ν := G** (Gram, L² pushforward, Parseval — the UNIQUE tight ν; signed Φ allowed) =
    today's frames; here `project` = `♯_G∘M` = `basis_space.apply_metric∘analysis`
    `[M]` 1.1e-16, and `reconstruction` is the section. Harmonic / orthonormal / nodal /
    two-table (PG) below it.
  - **ν := M·1** (the coupling's N-marginal, L¹ pushforward; Φ ≥ 0; Markov) = the PoU arm:
    `project` = E_γ[·|G] `[M]` array_equal on real EnergyGrids; `R∘project` idempotency
    defect 0.11 / 0.124 (nested 0); two-step condensation ≠ one-step by 0.55 on a
    straddling mid grid (nested 6e-16). NOT frame theory: the disintegration of a coupling.
  - One-hot Φ is the INTERSECTION (G = diag(M·1)): nodal frame = deterministic coupling.
- Morphisms, not nodes: `Descent`/G0, `conjugate`, `BulkLift`/`AngularLift`, Riesz legs.

## The theorem (two lines) that settles Q2
Bounds of `T_ν = W_ν⁻¹ΦᵀW_M` = nonzero spec of `W_ν^{-1/2} G W_ν^{-1/2}` — DEPEND on ν.
Canonical dual `R_can = Φ G⁻¹ W_ν`; `R_can T_ν = Φ G⁻¹ Φᵀ W_M` — the μ_M-orthogonal projector,
ν-FREE `[M]` array_equal through ν=V and ν=pushforward. ρ = the frame operator's spectrum on
the span for one-hot tables; canonical dual frame = frame / ρ; Parseval ⟺ ρ ≡ 1.

## Findings to carry
- `frame.py:583-585` "conjugate(gram_inverse) is the G-orthogonal projector onto span(basis)
  for every frame" is FALSE on the PoU arm (`f = T·[2,6]` in the span; P f = [8/3,4,16/3]).
- `GramStructure` = a tag choosing between the L² and the L¹ pushforward (2 dispatch sites,
  `frame.py:463, :632`); the DENSE refusal guards a path `basis_space.apply_metric` computes.
- Verbs that exist: `pushforward(φ).consolidate()` IS the aggregated pushforward
  (`[M]` array_equal with `diag(discrete_gram)`); the 2026-09-27 memo listed it as new.
- `project` is LOOSE; the tree's D4 words fit exactly: `average` = ♯_{M·1}∘M (bit-identical
  at all 14 callers), `marginalize` = M, `project` = ♯_G∘M. Three verbs, one slot.
- `OverlapBasis(IndicatorBasis)` — inverted edge (child general, parent one-hot special).
- `Kernel`: integral sense (5 classes, genuine, = K0) vs null-space sense (5 symbols, 3
  files) → the null-space family cedes to `Nullspace` (stronger unique alternative).
- Forced non-bit-identical: retraction through `frame.analysis` vs the shipped einsum
  7.1e-16 (G6.5 array_equal gate re-tiers); raw-sum divisors → Gram entry (1-ULP GL8 class)
  at `dsa.py:696,718,721`, `angular.py:290`.
- Option-B trap (SELF-CORRECTION): on harmonic frames `project` ≠ `apply_metric∘analysis`
  by the analytic dual factor in `reconstruct` (0.8–1.0 rel) — a CONVENTION, not a bug;
  only `product(4,4)` L=2 (measured DENSE, declared DIAGONAL) is a genuine latent
  non-canonical row-sum, 0 angular `project` callers.
- Naming hazard: `basis_space` holds covariant moments with metric g^{jk}=G⁻¹, so the
  tree's ♭ (`apply_metric`) on it is index-calculus RAISE.

## Pollination
CP matrix = positive kernel on (cells, VΣ_t); reciprocity = `.H` self-adjointness under that
metric; `_normalize_rcp` (`cp/solver.py:402`) = the disintegration `average`. MC tally
binning = `marginalize`, per-lethargy flux = `average` (D4).

## ⛔ SELF-CORRECTION (verbs attack, same day; memo `.../frame_hierarchy/verbs/memo.md`, probe `probe_verbs.py` P1–P7)
- **"The coupling does not compose associatively" above is WRONG.** Markov kernels compose
  (Chapman–Kolmogorov): two-step condensation == one-step through the COMPOSED kernel `T₁T₂`
  to 2.5e-16 in σ, χ, Σ_s, straddling or nested `[M]` P1. The 0.55 / 0.10 / 0.36 is
  `T_fine→coarse ≠ T₁T₂` — the overlap FACTORY (`EnergyGrid._overlap_table`) re-declares 1/E
  at the mid level; a non-functorial factory, not a non-associative kernel. Tower-property
  reading: exact iff the filtration nests (T₂ one-hot ⇔ `gram_structure DIAGONAL`).
  Provenance/route-to-parent is the right fix in kind, needed only when T₂ is fractional; in
  SPACE (all one-hot) the coarse object need carry only `K_*Φ` `[M]` P6.
- **"No new class; the PoU arm is the ν := M·1 filler" above is REVERSED.** Verb validity
  splits by filler: `conditional_expectation = (M·1)⁻¹M` is undefined on signed tables
  (`M·1 = [2,0,0]` on the L=2 harmonic frame; the tree's probe `MR·1 = [2,2,2]` is the GRAM
  row sum — they coincide iff `R·1 = 1`) `[M]` P5b; `projection = R(MR)⁻¹M` on PoU needs the
  dense solve nobody wants (0.45). Two verbs each invalid on one arm ⇒ two TYPES (D1):
  `Frame` (signed; ν ∈ {declared, Gram}) and `MeasureCoupling` (Φ ≥ 0; γ = μ_M ⊗ Φ; verbs
  pushforward / conditional_expectation / pullback / compose) as SIBLINGS under the abstract
  table-presented kernel; one-hot = the intersection; 14 of 14 `project` callers land on the
  coupling (they want conservation).
- **`projection` is orthogonal only on a Galerkin frame** `[M]` P4: on a one-hot PG frame with
  w ≠ const, `R(MR)⁻¹M` is idempotent (0) and not W-self-adjoint (1.29) — Christensen–Eldar's
  oblique projection. `dual_analysis` must be defined on the PAIR ((MR)⁻¹M; canonical iff test
  is trial); `canonical_dual_analysis` is false on 7 of 7 production PG frames.
- ⛔ **SUPERSEDED 2026-09-28 (evening) by `weight_placement_weld_attack.md`**: the bullet below is
  arithmetically right and inferentially wrong. P2's 1.6e-16 is one einsum in two orders (`array_equal`
  on the same call); the 7 are 7 of 15 collapse morphisms (5 are ratios of two pushforwards with
  different densities); the two spellings differ on identity, `.H` (2e-16 vs 0.66 on the CMFD
  prolongation), Gram sharing, signed weights, and cost (G³). The weld is REAL; the placement is the
  MULTIPLIER on a structural measure; `MeasureCoupling` is `analysis_c ∘ section_f`, a DERIVED composite.
- **Every production PG frame is a coupling in a costume** `[M]` P2+P7: 7 of 7 constructions
  carry a positive `WeightedIndicatorBasis` test and equal `GalerkinFrame(trial,
  weight-as-measure)` to 1.6e-16 (gram_inverse array_equal), INCLUDING the bilinear φ*⊙φ and
  ρ — so the recorded ground for "flux is a test weight, never a measure"
  (`weighted_indicator_basis.py:24-34`: adjoint weighting needs test ≠ trial) is refuted by the
  tree's own bilinear collapse. The only test weight that cannot be a measure is SIGNED
  (consistent-P_ℓ, φ_ℓ crosses zero `[M]` P5a) — 0 production sites. Reported, not ruled (#48).
- **RN formulation is the DENSITY rule only**: χ in energy is a MEASURE (Markov pushforward,
  no denominator); "push both and take the density" against counting divides by the group
  count (D4's wrong answer). `chi @ T == GalerkinFrame(T, counting).analysis(chi)` array_equal,
  ≠ the condensation PG frame's analysis by 0.36 `[M]` P3 — `marginalize` retires into the
  analysis of a SECOND frame. χ is measure-valued in energy, a density in space.
- Verb verdicts: analysis KEEP; reconstruction KEEP; dual_analysis ACCEPT (pair definition);
  projection KEEP corrected; conditional_expectation ACCEPT on the coupling only;
  marginalize/average/project RETIRE. `dual` is spent on the ALGEBRAIC dual (`.dual()`, "dual
  dyad") — a different sense, docstring must say which.
- HarmonicFrame: a licensed LEAF (four binder verbs returning operators, D10) on the FAMILY
  axis inside the DISCIPLINE chain; the grid appears at the second family (the hand-rolled
  Legendre analysis, `radial_characteristic_field.py:392,404`).

Cross-refs: [[pullback-pair-frobenius-perron-frames]] (ρ; the pair IS the nodal frame),
[[condensation-nonnested-fractional-overlap-frames]] (⛔ its "rate-preservation falls out of
frame.project as a frame projection" reading is SUPERSEDED: it falls out of the coupling's
disintegration, which is NOT a projection), [[xs-coarsening-collapse-marginalize-vs-average]]
(the verb names), [[unified-frame-api-design]], [[feedback-naming-collision-perfect-match]].
