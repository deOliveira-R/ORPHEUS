---
name: fundamental-mode-gauge-attack
description: Verdict of the 2026-09-28 W5 attack on "how is the scale of an eigen answer's mode carried" (posing_sequence B13) — a MODE has no scale (its Laurent residue ψ⊗ψ†/⟨ψ†,T_d ψ⟩ is invariant under both rescalings, [M] 1.6e-16, contour = formula 1e-14); the scale is a SECTION datum of a REPRESENTATIVE, i.e. the outcome's ScaleGauge; no FundamentalMode type; the admissibility of a positive-functional section is NOT the fundamental/harmonic partition ([M] 20 of 20 k-modes, 240 of 240 α-modes admit it on an asymmetric body; parity on a symmetric one); proposed LAW: every stored Mode pair canonical (⟨ψ†,T_d ψ⟩=1, ‖ψ‖_G=1). Open before any mode / harmonic / modal-expansion / gauge / Enclosed brief.
metadata:
  type: project
---

# A mode has no scale; the scale is a representative's section (2026-09-28)

Memo: `scratch/posing_sequence/open_items/fundamental_mode/memo.md`; probes `probe.py` (P1–P6),
`probe2.py` (P7–P8), `census.py`, on the final-attack builder (DD S8 slab, 240 unknowns; ASYM /
SYM / THICK variants). Re-ground every `file:line` against the tree before acting.

**The theorem that settled it** `[M]`: `Res_n = ψ_n ⊗ ψ†_n / ⟨ψ†_n, T_d ψ_n⟩` is invariant under
`(ψ, ψ†) ↦ (aψ, bψ†)` (1.6e-16), the Riesz projector `P_n = Res_n T_d` is idempotent (3.5e-16),
and the contour integral of `at(σ)⁻¹` equals the formula on the fundamental, a harmonic and the
α fundamental (1.2e-14 / 3.3e-14 / 5.2e-14). So "ask the mode for its gauge" has nothing to
return; a gauge is data about a REPRESENTATIVE — the `ScaleGauge` docstring's own definition.

**U's partition refuted three ways** `[M]`: (1) the power section `⟨1, F ·⟩` exists on 20 of 20
k-modes and 240 of 240 α-modes on an ASYMMETRIC body (readings 2e-2–7e-1); on a SYMMETRIC body
it vanishes on exactly the reflection-odd modes (1e-16; even harmonics 1e-2 as the positive
control) — a PARITY partition, not a kind; (2) the DD fundamental leaves the cone (14 of 240
negatives at h = 3 mfp) and its section stands (+4.81, lands at 2.2e-16); (3) the tree gauges
the FUNDAMENTAL by sign + 2-norm at 2 dense-engine sites and the ADJOINT ψ† by a positive
functional (`sn/solver.py:2796`), and on the diffusion member ψ† = ψ (1.2e-13) so the
perturbation-theory bilinear gauge IS a norm on every mode.

**The functional that DOES single out the fundamental** `[M]` P7: its own adjoint pairing
`⟨ψ†_0, F ·⟩` — 0.967 on the fundamental, ≤ 7.7e-15 on the 19 harmonics (biorthogonality,
off-diagonal 2.7e-14); it separates every mode from every other equally. The bilinear gauge
`⟨ψ†_n, T_d ψ_n⟩ = 1` exists iff the pole is simple (260 of 260) and is spellable TODAY as
`ScaleGauge(Riesz(T_d ψ).evaluate, 1.0)` on ψ† (residue identical under the tree's section and
the bilinear one, 1.4e-17).

**Physics** `[M]`: a modal expansion / transient needs each mode's PAIR and its direction's
`T_d`, never a common gauge (unit vs random scales 3.0e-12); the amplitude is `x₀`'s
(`Evolution.initial`); the subcritical source's fundamental is `Res_0 q/(σ − σ_0)`, scale-free
in ψ₀, the scale `q`'s.

**Verdict**: keep the outcome's gauge; no `FundamentalMode`; `Fundamental` stays the selector
(the stem's formal perfect match). Proposed: (a) a LAW on `Mode` — stored pair canonical
(`⟨ψ†,T_d ψ⟩ = 1`, `‖ψ‖_G = 1`, sign rule; phase rule for a complex pole UNMEASURED) — so a
multi-mode outcome cannot re-open the "recorded nowhere" defect; (b) `ScaleGauge.displacement`
in the state's dtype (today `float(...)` turns a complex reading into a wrong number under a
one-time `ComplexWarning`). Gates: contour vs formula; the canonical-pair invariant;
`allclose(state, gauge.apply(modes[fundamental].psi))` (one-scalar redundancy, kept).

**Census** `[M]`: `ScaleGauge(` 3 production sites; production READS of the scale gauge 0 (all
`.gauge` loads are `certificate.gauge` / the nullspace verb); kind branches 5, all on
`EigenOutcome` vs `SourceOutcome`; `Mode`/`FundamentalMode`/`Fundamental` 0 of 988 files.

**Frame that fires for ANOTHER question**: the symmetry quotient (Γ-invariant functionals
vanish on non-trivial isotypic components) — decisive for "which modes does a symmetric
detector see / harmonic excitation on symmetric cores"; refuted FOR the scale question.

Related: [[equilibrium-carrier-laurent-point-kinetics]] (the residue as the transition
carrier), [[posing-ontology-clean-attack-frames]] (the bulk-cone `Fundamental`),
[[spatial-order-type-vs-property-criterion]] (D1, applied here: the morphism is a scalar).
