---
name: pullback-pair-frobenius-perron-frames
description: Restriction/extension/retraction/section/deck arrows as the (Koopman pullback, Frobenius–Perron adjoint) pair of ONE point map plus TWO measures; ρ = P_φ1 (the pushforward's density) is the one datum the map lacks and explains r∘e=id, dead slots, ERR-042 and the Σw two-type split; the pair IS the fibre-indicator GalerkinFrame (discrete_gram = aggregated pushforward measure) up to one Riesz leg; BulkLift is a consumer (conjugation), the coproduct rows are unspellable (no Coproduct manifold); "transfer operator" collides 3 ways.
metadata:
  type: project
---

# Pullback pair — the (Koopman, Frobenius–Perron) frame for restriction/extension

Attack 2026-09-27 on `posing_sequence.md:787-806` (memo
`scratch/posing_sequence/pullback_attack/structure/memo.md`, probe
`probe_pullback_pair.py`). Re-ground every `file:line` before acting.

## The one law
`φ: M → N`, `μ_M`, `μ_N`. `φ*g = g∘φ` (gather = Koopman `U_φ`); its Hilbert adjoint
`φ_♯ f = (1/w_N) Σ_fibre w_M f` (weighted scatter-add = Frobenius–Perron `P_φ`,
Lasota–Mackey 1994 ch. 3); `ρ = P_φ 1 = d(φ_*μ_M)/dμ_N`. Then `φ_♯φ* = M_ρ`, section
`e = φ* M_{1/ρ}`, `r∘e = id ⟺ ρ ≡ 1 on the image`; dead slot ⟺ ρ = 0; no section ⟺ ρ = 0
on a live slot (Σw = 0); deck measure-preserving (ERR-042) ⟺ `P_σ1 = 1` ⟹ `G.H = G⁻¹`.
`[M]`: trace γ₊/ι₊ array_equal ×3; retraction R = φ_♯ (2.8e-15), R.H = π* (1e-16),
E = π*M_{1/ρ} array_equal, ρ = Σw = 2.0; with μ_N := π_*μ_M the section IS the pullback
(array_equal) — the anti-ERR-051 two-type design is the SYMPTOM of choosing the physical
`V` (not the pushforward `Σw·V`) as the marginal's measure. Deck: apply/.H/inverse
array_equal, ρ = [1.].

## The pair IS a frame (Q4 verdict)
`GalerkinFrame(IndicatorBasis over N's nodes, μ_M.pushforward(φ))`: `table` = one-hot
pullback matrix; `diag(discrete_gram)` = AGGREGATED `φ_*μ_M` (array_equal; `pushforward`
itself keeps duplicates, `measure.py:945`); `reconstruction = φ*`; `analysis = ♭∘φ_♯`
(covariant); `basis_space` metric = `1/φ_*μ_M` = the dual `N*`; `analysis.H =
reconstruction∘G⁻¹` (the frame's own section). The tree's PRIMAL arrows are one Riesz
leg away: `R_tree = ♯_{(N,V)}∘analysis`, `E_tree = reconstruction∘G⁻¹∘♭_{(N,V)}` (≤1.2e-15).
The frame's G0 `descent` IS the point map (today only `quotient_onto` arrows;
`_collapse_pair` bypasses by putting φ into the measure's NODES — the rank-one case
`N = {pt}`). ⇒ no new binding type (a twin); new VERBS: `DiscreteMeasure.nodal_basis()`,
the aggregated pushforward, the primal raise onto a chosen μ_N.

## Per-family verdicts
- `TraceRestrictionOperator`: pullback along `Γ₊ ↪ Γ`; ρ ≡ 1 on image. Falls out.
- `AxisRetraction/Section`: `π_*` and `π*M_{1/ρ}`; needs μ_N (the physical V) — the
  section's divisor is the Jacobian `ManifoldMap` disowns (`manifold.py:1175`).
- `SystemRestriction`, bulk injection: pullback along a COPRODUCT injection — the manifold
  algebra has products and quotients, NO coproduct (`[M]` 0 hits); dual-role zeros are
  carrier typing outside the measure algebra.
- `BulkLift`: NOT a member; `lift(A) = φ_♯ A φ*` = the frame's `conjugate` shape; a consumer
  like `T = γ_out (L+C)⁻¹ ι_in`.
- Specular/periodic deck: Koopman of the mirror/translation (`_factors.py:353` already says
  "composition operator"); discrete only where nodes map to nodes (`product(4,3)`: none);
  an interpolated pullback would lose `P_σ1 = 1` hence unitarity — refusal principled.
- Cell-to-cell: identity; the span degenerates.

## The Stokes weight
`|Ω·n| dS dΩ` = `i_Ω(dV)|_∂V` (flux form, per-unit-time volume flux) — a datum of the
streaming FIELD, SUPPLIED as μ_N; not a pushforward along `∂V ↪ V`. The exit-map pushforward
of `dV⊗dΩ` is the CHORD-weighted `ℓ(Ω·n)dSdΩ` (CP escape measure); `∫_{Γ₊} ℓ(Ω·n) = 4πV` is
Cauchy's mean chord — the CP↔SN boundary-measure gate. Tangential band = kernel of μ_N ⇒
ρ = 0/0 ⇒ "neither half" is FORCED (non-singularity), the ε is a metric-kernel classifier.

## Naming (D6 + the collision ruling)
Koopman = pullback (exact, docs); Frobenius–Perron = adjoint (exact, NO collision);
"transfer operator" KILLED — 3 senses in tree (`TransferOperator` scattering
`transfer.py:22`; "Birkhoff transfer operator" `billiard.py:28/51/174` = this `P_φ` for the
billiard map; Ruelle); "correspondence" (span, degenerate), "Galerkin morphism" (invented),
"transport map" (Monge; collides) killed; "nonsingular transformation" = the faithful name
of the BINDING but a type would twin the frame; `PullbackOperator` is right for the ARROW
and subsumes TraceRestriction (injective idx) / Permutation (bijective) / the section's
broadcast (surjective) — 3 gathers + 3 partners `[M]` `operator.py:2892/2895, 3226/3244,
3351/3561` → 1 + 1, guards becoming map properties (sortedness inherited from the bound
codomain space's row order, ERR-077).

## MC pollination
Deck `G` and response `R` are both FP operators of Markov kernels on the trace: Dirac
(`P_G1 = 1`) vs stochastic (`P_R1 = α`, the albedo). One frame for geometry-vs-constitutive.

Cross-refs: [[field-role-typing-faceflux-frames]] (trace = ι*), [[rep-role-grid-double-category-frames]]
(nodal vs modal N-side = one binding), [[mc-seam-feynman-kac-particle-frames]],
[[feedback-naming-collision-perfect-match]], [[dsa-rp-angular-frame]] (ANTI-MINT precedent).
