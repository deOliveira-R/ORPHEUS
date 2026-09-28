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

Cross-refs: [[pullback-pair-frobenius-perron-frames]] (ρ; the pair IS the nodal frame),
[[condensation-nonnested-fractional-overlap-frames]] (⛔ its "rate-preservation falls out of
frame.project as a frame projection" reading is SUPERSEDED: it falls out of the coupling's
disintegration, which is NOT a projection), [[xs-coarsening-collapse-marginalize-vs-average]]
(the verb names), [[unified-frame-api-design]], [[feedback-naming-collision-perfect-match]].
