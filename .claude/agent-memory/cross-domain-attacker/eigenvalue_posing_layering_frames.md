---
name: eigenvalue-posing-layering-frames
description: Durable frame — the K/α/source/transient deterministic-transport eigenproblem family is ONE generalized eigenproblem Aψ=λMψ with a resolvent backbone K=A_loss⁻¹M, layered as leaves→posing→resolvent→algorithm. KEY correction: posing is NOT method-agnostic; it bifurcates into 2a (method-agnostic role assignment + μ→physical map) + 2b (method-specific A_loss realization). Adjoint + transient are sibling posings over the same leaves; metric lives at the leaf layer.
metadata:
  type: project
---

# Eigenvalue/transient layering — standard form + 4-layer frame

Durable structural verdict. Builds on [[power_iteration_vs_keigenvalue_morphism]]
and [[issue_208_operator_algebra_frames]]. The spectral-theory (generalized
eigenproblem) frame and the fixed-point-combinator frame are load-bearing.

## Standard form

`Aψ = λMψ` (generalized eigenproblem). K (`A_K=L+C−S−B`, `M=F`, λ=1/k) and α
(`A_α=L+C−S−F−B`, `M=T=1/v`, λ=−α) are TWO POSINGS of ONE eigenproblem. The cleaner
operational backbone is the **power-method resolvent** `K_pm = A_loss⁻¹ M` (a
Krein-Rutman compact-positive operator; its dominant eigenpair IS the fundamental
mode) — exactly what a power-iteration loop iterates. Shift-invert `(A_loss−σM)⁻¹` is
the STRICT generalization (interior eigenvalues, FEAST/Arnoldi-Schur), NOT needed for
the dominant K or α (both are EXTREME eigenvalues of `A_loss⁻¹M`, reachable by plain
power iteration). So: standard form = generalized eigenproblem; backbone = the
resolvent `A_loss⁻¹M`; shift-invert = a documented future seam.

## The 4-layer frame, with the ONE correction (the durable lesson)

leaves (1) → posing (2) → resolvent (3) → algorithm (4). The naive claim that
**posing is method-agnostic** (just arranges leaves) while the resolvent is
method-specific is WRONG: the loss-operator ASSEMBLY is itself method-specific (SN
folds B into the scattering slot `S+B`; CP has no L/S/F arrangement at all — its
matrix is monolithic). The fix is to SPLIT the posing layer:

- **2a — Posing (method-agnostic ROLE assignment):** the table
  `{K: loss=all-but-F, eigen=F, μ→k=μ; α: loss=all-incl-F, eigen=T, μ→α=−1/μ;
  source: loss, q; transient: loss=A_loss+T/Δt-shift, q=T/Δt·ψⁿ+qⁿ⁺¹}`. Pure data:
  which leaves play which role + the μ→physical map. Genuinely method-agnostic.
- **2b — Loss-operator REALIZATION (method-specific):** turn the role assignment
  into a concrete invertible `A_loss` object (SN's streaming+collision triple; CP's
  monolithic matrix). This is layer-3-adjacent — HOW the family builds the thing the
  resolvent inverts.

Once split, 2a∘2b∘3∘4 composes cleanly. The general recurring trap: treating an
"adapter" as a problem-type layer when it is actually a 2b realization of the
resolvent (cf. [[power_iteration_vs_keigenvalue_morphism]]).

## Sibling posings (free under the frame)

- **Transient = a sibling ALGORITHM over the SAME leaves + a transient posing.**
  Backward-Euler `(T/Δt+A_loss)ψⁿ⁺¹ = T/Δt·ψⁿ+qⁿ⁺¹` is a fixed-source solve with a
  SHIFTED loss operator. It SHARES the resolvent (the same sweep/BiCGSTAB inverts
  `A_loss+T/Δt`), needs NO new resolvent, NO new leaves except T (a
  `DiagonalOperator(1/v)`, a 6th leaf joining L,C,S,F,B). Cleanest possible fit.
- **Adjoint = another posing row.** The adjoint eigenproblem `A†ψ†=λM†ψ†` has the
  daggered leaves as its role-operators; the dagger functor is FREE from the
  biproduct category (Wave O), λ is unchanged (forward/adjoint share the spectrum),
  so only the role-operators dagger. Adjoint slots in at 2a with zero new machinery.
- **The metric `|Ω·n|·w` lives at the LEAF layer, not the posing layer.** It is the
  codomain inner product of the streaming leaf's boundary-trace block — intrinsic to
  the leaf (its domain/codomain FunctionSpaces carry `inner_product_weights`). Posing
  arranges leaves; the leaves already know their metric.

## Cross-domain detection heuristic (durable)
When a family of "find-the-special-value" problems (k, α, time-step, source) looks
like distinct solvers, check whether they are POSINGS of one generalized eigenproblem
sharing a resolvent backbone. If they are, the only genuinely per-method layer is the
loss-operator REALIZATION (2b); everything else (role assignment, μ-map, adjoint,
transient shift) is method-agnostic data over a shared engine. Respects
unify-after-two: one instance (K) is enough to lay the layering down so #2 (α) drops in.

Related: [[power_iteration_vs_keigenvalue_morphism]] (the engine + adapter),
[[issue_208_operator_algebra_frames]] (dagger-biproduct → adjoint posing free, metric
on leaves), [[field_role_typing_faceflux_frames]] (the B-as-leaf / biproduct
C¹=C¹_int⊕C¹_∂ structure the transient T-leaf joins).
