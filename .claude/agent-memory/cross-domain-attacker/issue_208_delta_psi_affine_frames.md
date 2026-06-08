---
name: issue-208-delta-psi-affine-frames
description: Durable frame triad for typing the two "→0" quantities of any iterative solver — the iterate increment Δx (a difference-space element) vs the equation residual r (a dual-space defect). Affine-torsor + Banach-contraction + Krylov-dual. The affine/torsor frame + the state-typed-increment anti-pattern were PROMOTED TO THE SKILL (commit 9504146); this note keeps the reusable triad + the cross-method generalization.
metadata:
  type: project
---

# Typing Δψ (increment) vs r (residual): the affine / Banach / Krylov triad

PROMOTED: the affine/torsor frame and the "state-typed increment" anti-pattern are
now IN THE SKILL (commit `9504146 docs(skills): promote the affine/torsor frame +
the state-typed-increment anti-pattern`). `FluxDisplacement` shipped (Wave O #208).
What survives here is the reusable THREE-FRAME triad and why all three coexist.

## The structural question (recurs in every iterative method)

Two distinct "→0" quantities live in any fixed-point/Krylov solver:
- the **iterate increment** `Δx = xⁱ − xⁱ⁻¹` — in the STATE's units (flux);
- the **equation residual** `r = Ax − q` — in the OPERATOR's codomain units
  (rate-density).
They are NOT the same object and NOT the same type. `Δx = A⁻¹r̃` connects them.

## The triad (three frames, all firing, none redundant)

**1. Affine geometry / torsor (A.1).** Trigger: the state set has NO natural origin;
only `Σλᵢxᵢ` with `Σλᵢ=1` preserves the role. The state set is an affine space `A`
over a DISTINCT difference vector space `V`; the increment lives in `V`
(`FluxDisplacement` — shape-identical to the state but a distinct type/role). Torsor
axioms ARE the three ops: `⊖:A×A→V` (this IS Δx), `⊕:A×V→A` (the update; relaxation
`x_old ⊕ ω·Δx`), and `A×A→A` UNDEFINED. `V` is a vector space (`+`, scalar·, zero,
norm = step length); `A` is NOT. Relaxation `ωx_new+(1−ω)x_old` is legal (`Σλ=1`);
`x+x` is illegal (`Σλ=2`, lands off the affine subspace). The "no state+state" rule
becomes a TYPE CONSEQUENCE, not a runtime guard.

**2. Banach fixed-point / contraction (A.5).** Trigger: `Δxⁱ⁺¹=M·Δxⁱ` is an exact
linear recurrence on increments; convergence checked empirically (`res<tol`). The
displacement type is the ONLY object that knows "previous"/"step", so it is the
NATURAL carrier of contraction data: `ρⁱ=‖Δxⁱ⁺¹‖/‖Δxⁱ‖`, the dominant eigenmode of M
(displacement converges in DIRECTION to v₁, STATE in VALUE to x* — two limits, two
types), Aitken Δ² `x*≈xⁱ⊕Δx/(1−ρ)`, a-posteriori bound `‖x*⊖xⁱ‖≤‖Δx‖·ρ/(1−ρ)`. The
state type has NOWHERE to put any of this → the loop discards Δx after the scalar norm.

**3. Residual/Krylov/Galerkin (A.5 Krylov + A.3 Fredholm).** Trigger: `Δx=A⁻¹r̃` — the
two "→0" quantities are the SAME defect in the two universes connected by A. `‖Δx‖`
(state units) = Richardson/SI-native criterion (PRIMAL step); `‖r‖=‖Ax−q‖`
(rate-density) = Krylov/Galerkin-native criterion (DUAL; GMRES minimizes exactly
this). The residual carries what the increment can't: per-cell defect map, A-norm
energy `⟨r,A⁻¹r⟩`, the DSA correction source (= r directly). Criterion choice is a
FRAME choice (primal vs dual), not a convention.

## The resolution (reusable shape)
Mint the displacement type as the difference space V (defined ONCE on the field base
so every leaf gets a displacement sibling FREE — role-parameterization). Drop
`state.__add__`; retype `state.__sub__ → Displacement`; add
`state.__add__(Displacement) → state` (torsor action) + an `affine_combination`
(`Σλ=1`) classmethod. Keep BOTH stop tests typed: SI on `‖Displacement‖` (primal),
Krylov on `‖Residual‖` (dual — gives the residual its first production consumer + the
DSA source). The flat-norm "just forbid state+state" option is weakest — it patches
the symptom, names no role, strands the contraction data and the residual.

## Cross-method generalization (durable)
- **Eigenvalue/power-iteration: SAME affine+contraction structure.** The dominance
  ratio `k₂/k₁` IS the contraction factor of `M=(1/k)L⁻¹F` — the displacement's
  `contraction_ratio` is shared verbatim. (Ties to [[eigenvalue_posing_layering_frames]].)
- **CP/MOC/diffusion:** `Displacement[F]`/`Residual[F]` parameterized over the flux
  leaf F generalizes cleanly (the #208 biproduct precedent — derived, not re-declared).
- **Krylov→SN:** residual-norm stop + A-norm energy make the residual the GMRES carrier
  + DSA precond source.

REFUTED / low-signal for this typing problem class: topology, group theory, tensor
networks, Hilbert-Schmidt, multiplication-operator spectral, measure theory,
information theory, category theory — no trigger (single-mesh, angular-shape-agnostic
typing). Category theory's role-parameterization win is captured concretely by
affine+Krylov+biproduct; no abstract-nonsense lever needed.
