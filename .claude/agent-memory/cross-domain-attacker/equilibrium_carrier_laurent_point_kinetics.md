---
name: equilibrium-carrier-laurent-point-kinetics
description: Durable frame (first sighting 2026-09-26, the posing-sequence final attack / transition question) — what flows from an eigen answer to noise, a transient or an at-pole source is an EQUILIBRIUM (a point of the spectral hypersurface Σ ⊂ Λ with its null representative and scale), not a point of Λ; the at-pole gauge is the DEPARTURE direction's (a stationary splitting conserves ⟨ψ†, N x⟩); the Laurent expansion along TIME IS point kinetics (ρ/Λ polar part, D-gauged shape); the eigen answer's conormal to Σ = adjoint-weighted term rates (Hellmann–Feynman); Keldysh/Gohberg–Sigal gives the non-affine family the same Laurent structure; enlargement's lift = Schur back-substitution. Open for any noise / transient / GPT / "change the problem kind" / perturbation-theory brief.
metadata:
  type: project
---

# The transition carrier is an equilibrium; the at-pole gauge is directional

Memo: `scratch/posing_sequence/final_attack_structure_memo.md`; probe
`final_attack_structure_probe.py` P1–P6 (80-unknown S8 DD slab, non-uniform h and 1/v,
divergence form, metric `G = diag(h w)`). Re-ground every claim against the tree.

**The carrier.** `[M]` 4 of 5 derived questions (noise, transient, static-at-pole, GPT, and
fixed→eigen) share the POLE and 3 of 5 the STATE; 0 of 5 need the direction that PRODUCED the
pole; 4 of 5 need the direction of DEPARTURE (TIME for every dynamic question, FISSION for
GPT). So the carrier = an equilibrium of the flow `D ẋ = −T(c)x`: a point `c* ∈ Σ` plus the
null representative with its scale (= `EigenOutcome`); ψ† is derivable (nullary adjoint) and
needed only AT the pole. A derived question's method-free VALUE names its PARENT QUESTION
(`Noise(δΣ, ω, of=Eigen(FISSION))`), the binder resolves the parent's outcome in the loaded
Solution (content identity on the discretisation — a load has no object identity). No
transition verb; no stored point. The primitive `FixedSource(q, point)` stays for questions
that own a point (subcritical source at the anchor; beam-modulated ADS). Fixed→eigen carries
nothing but the discretisation ("just pick the other problem" is right in THAT direction only).

**Offset-0 law** `[M]`: re-posing along any direction from a pole finds it at offset 0
(`Eigen(TIME, at=k-pole)` → α 5e-15; `Eigen(FISSION, at=α-pole)` → k = 1.000000000000; at the
anchor instead k = 1.0736). `Eigen(parameter, at=point)`'s component along the parameter is
projected out (the LINE is the invariant; the coordinate is absolute).

**Gauge is directional** `[M]`: any stationary splitting at a pole conserves `⟨ψ†, N x⟩`, so
its limit carries `⟨ψ†, N x⟩ = ⟨ψ†, N x₀⟩` — the Strategy's explicit part and x₀ pick the
gauge (a question/Strategy gate violation unless the question projects). Fission splitting →
F-gauge; backward-Euler pseudo-transient → D-gauge; noise ω→0 → D-gauge; difference a
ψ-multiple, norm 5.6. GPT's `⟨Γ†, Fφ⟩ = 0` and the quasistatic `⟨ψ†, δψ/v⟩ = 0` are ONE law
along two directions. "x₀ is a Strategy datum" is true off a pole, false at one.

**Laurent along TIME = point kinetics** `[M]`: `R(c* + iω e_T)(−δTψ) = (ρ/Λ)/(iω)·ψ + S_T q`,
ρ = −⟨ψ†,δTψ⟩/⟨ψ†,Fψ⟩, Λ = ⟨ψ†,Dψ⟩/⟨ψ†,Fψ⟩ (polar coefficient to 5e-7; shape bounded and
D-gauged; F as denominator diverges). Solvability law = "δρ = 0"; secular growth = the
reactivity ramp. With one precursor group `Λ_eff = Λ + β/(λk*)` EXACT (8.9e-16) — the
k-normalised equilibrium carries `1/k*` (NOT the textbook `Λ + β/λ`).

**Conormal = adjoint-weighted term rates** `[M]`: `n_j = ⟨ψ†, T_j ψ⟩_G`; `dσ_d*/dc_a = −n_a/n_d`
matches FD to 2e-10; the SAME metric adjoint paired by a Euclidean sum is 4% off on a
non-uniform mesh (the #517 error class) — the gate for "the VJP must be the metric adjoint".
The ruled "rate functional by channel" (adjoint-weighted) and "perturbation theory is the
pole's derivative" are one computation: `outcome.sensitivity(term)`.

**Keldysh / Gohberg–Sigal** `[M]`: a rational family has poles, finite-rank residues and the
logarithmic-residue count `(1/2πi)∮ tr[T⁻¹T'] dp` (1.000000000 around the dominant critical p;
a contour also enclosing the coefficient's pole reads 6 − rank F = 0, the meromorphic count).
"Not a resolvent" (no identity) ≠ "no Laurent datum". `Enclosed(region)` is uniform across
affine/non-affine with `T'` for `−M`; "root over pencils" is a Newton STRATEGY.

**The precursor block is LAYER 1; the total-ν family is its Schur LOWERING** `[M]`: at s=0,
u=0 the Schur complement IS the reduced pencil with total ν (9e-16) and the enlarged static k
equals it to the printed digits; with fuel FLOW (upwind `u·∇C`, recirculation closure with
`e^{−λτ_ext}`) the STATIC k moves by −1.9e-3 (u=0.5) to −6.0e-3 (u=2, τ_ext=20) — the block
changes the anchor's balance, hence layer 1 (five terms: decay, advection + W's boundary law,
time, production, delayed emission). The precursor equilibrium is the lowering's
BACK-SUBSTITUTION (`E(c*,0)·(ψ, C_eq) = 0`, 5e-15; the total-ν lift fails at 0.57), so the
Solution's state lives on `V ⊕ W` and the transition has no lift. The W block is a SWEEP:
lower-triangular up to ONE closing entry (D2 rank 1), `(λ + u·∇)⁻¹` = the transport resolvent
with one direction and the decay as Σ_t; `s + λ` folds into its diagonal free. "Static vs
dynamic delayed neutrons" = the s=0 vs s≠0 POINT of one family, not two term sets.

**Strategy transfer criterion**: a Strategy datum transfers across a change of question iff
it is a function of layer-1 declarations alone (D1/D2 labelling, schedule, lowering: yes; D3
ranking: only while the point is real and cone-preserving; the outer eigen loop: no; x₀: off a
pole only).

**Detector hits to reuse**: an answer-changing datum outside the question value (the unnamed
M); a hazard LIST where a law is missing; "anchor" spent twice in one plan; a sign living in
two places (`balance_derivative = +D` ⇒ α chart = identity; only FISSION is reciprocal).

Related: [[shift-ontology-taxonomy-frames]] (its at-pole paragraph superseded here),
[[pencil-pseudo-resolvent-standard-form-pair]], [[question-stage-coefficient-space]],
[[eigenvalue-posing-layering-frames]].
