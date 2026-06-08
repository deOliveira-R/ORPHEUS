---
name: power-iteration-vs-keigenvalue-morphism
description: Durable frame — a power-iteration eigenvalue loop and an operator-triple eigenvalue loop are the SAME fixed-point combinator at two layers; the Protocol-opaque-resolvent layer is strictly more general than the operator-triple layer (admits monolithic-matrix resolvents that have no L/S/F split). The deprecation arrow points the wrong way: the opaque layer is the engine.
metadata:
  type: project
---

# power_iteration vs KEigenvalue — same morphism, two layers

Durable structural verdict (the #208/eigenvalue-layering landed against it; the
file:line recipe is now in the code — what survives here is the FRAME).

## The morphism (named structure)

A k-eigenvalue power-iteration loop is the **Banach / Krasnoselskij fixed-point
iteration of the shifted power operator** `T_k : ψ ↦ (L−S)⁻¹ F ψ` composed with a
Rayleigh-quotient eigenvalue update — the **power-method combinator** `fix(step)`:

```
fission  = F·ψ / k            # rank-1-ish source build
ψ_new    = (L−S)⁻¹ fission     # inner resolvent solve  ← the morphism's core
ψ_new   /= P(ψ_new)           # production renorm, scale-invariant in k
k_new    = R(ψ_new)           # Rayleigh / production-loss quotient
converge?(k_new, k_old, ψ_new, ψ_old)
```

When the same project hosts TWO eigenvalue loops with structurally identical bodies,
they are not "different algorithms" — they are this one combinator with the inner
resolvent `(L−S)⁻¹` exposed at two different layers:

- **Protocol-opaque layer:** the resolvent is BEHIND a `solve_fixed_source` Protocol
  method — an opaque morphism the solver owns.
- **Operator-triple layer:** the resolvent is BUILT FROM an `(L, S, F)` triple —
  `SourceIteration(L, S, 0).solve` constructed inside the eigenvalue routine.

## The decisive layer-ordering fact (reusable)

The **Protocol-opaque layer is strictly more general** and therefore canonical.
It admits BOTH the operator-triple resolvent (SN, MoC) AND the monolithic-matrix
resolvent (CP, diffusion, homogeneous — a single BiCGSTAB/direct solve on one matrix
with NO `(L,S,F)` streaming/collision split). Forcing the families that have no
triple onto a triple-layer engine requires MANUFACTURING fictitious `L`, `S`
operators — `(L−S)⁻¹` does not factor out of a CP/FD matrix. So a "retire the opaque
loop, migrate everyone to the triple loop" plan is REFUTED: the deprecation arrow
points the wrong way. Correct single-source-of-truth shape: ONE engine = the opaque
Protocol combinator; the triple loop REFRAMES as a thin adapter (the operator-triple
Protocol implementer) that DELEGATES the loop to the engine. The duplicated loop body
retires; the triple instantiation stays as one Protocol impl among the families.

## Cross-domain detection heuristic (the durable lesson)

When two iterative drivers share their loop body verbatim, ask which layer each
exposes the inner resolvent at. The driver whose resolvent is OPAQUE (behind an
interface) is almost always the general one; the driver whose resolvent is a
CONCRETE algebraic factorization is the specialization. Generality flows toward the
opaque interface, never away from it. A "type-agnostic primitive regressed feature X"
objection usually traces to a TEST adapter that collapses state (e.g. an angular→scalar
reduction between outer iterations dropping the Pℓ moments) — NOT to the engine, whose
shape-agnostic primitive carries the state by construction. Attribute the regression to
the lossy adapter, not the engine.

## Cross-method pollination
- MoC eigenvalue ALREADY consumes the opaque Protocol — structural twin of SN
  (ray-trace resolvent behind `solve_fixed_source` instead of sweep resolvent).
  Confirms the Protocol layer is the shared transport backbone.
- Trajectory-resolvent backbone (`power_iterate_variant_alpha`, the continuous
  boundary-to-boundary `T=(I−S)⁻¹`) is the SAME combinator on a DIFFERENT object
  (continuous resolvent vs discrete). The discrete `K=A_loss⁻¹F` is its shadow —
  which is WHY the rank-N/resolvent framing keeps recurring across these attacks.

Related: [[eigenvalue_posing_layering_frames]] (the full K/α/source/transient
layering this morphism sits inside), [[issue_208_operator_algebra_frames]] (the
(L,S,F) biproduct algebra), [[trajectory_resolvent_foreign_frames]] (continuous twin).
