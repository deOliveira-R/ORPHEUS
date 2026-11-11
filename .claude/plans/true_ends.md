# True ends: every operator declares the space its members live in (#707, then #705)

Living plan, opened 2026-11-10. The ontology is being searched (`plan-authoring` §0), so this file is refined in place as the discussion runs. Implementation waits until the design below is reviewed and ruled.

## Why

#705 (a round-trip shortcut such as `A.H.H is A` is taken only where a checked predicate proves it) needs to mint a member of an operator's input space. Doing so exposed that operator ends lie in two ways.

1. **Axes.** Some operators declare a space smaller than the one their members live in.
   - `LegendreMomentTransfer` (Λ) and `FissionMomentOperator` declare the angular-moment head but act on head ⊗ group ⊗ cells.
   - `FrameAnalysis` and `CrossGramInverse` do the same on the Petrov–Galerkin homogenisation frame.
   - On the windowed path the SN boundary operators are handed a member of another space.
2. **Role.** Spaces carry no physical role. 25 operator classes (16 atoms) change role while declaring the same space at both ends: `Lψ = q` and its sweep both declare V → V. `[M]` production hands flux to the sweep 19 448 times through `FlattenedOperator`, and every #705 failure outside Λ is a source-consuming inverse handed a flux.

Census: `scratch/true_ends/census/report.md` (`[M]` 2026-11-10, 11 920 passed, 23.3 M wrapped calls, HEAD `a00539f4`).

## Rulings (the user, 2026-11-10)

- **R1, the root.** "Do the whole thing and make this problem disappear at the root", in any order.
- **R2, axes.** Every operator's ends are the space its members live in. An operator acting on some axes names them through the tensor algebra (`TensorProductOperator`'s axis-tagged factors, spaces by agreement). The implicit "trailing axes ride through" convention retires.
- **R3, reading (b)** (replaced reading (a) the same day).
  - The state side **V** holds functions: cell-average fluxes, nodal values, the importance ψ†.
  - The data side **V\*** holds integrated quantities: the source in a cell, a partial current through a face, a detector response.
  - The measure ♭ = G is the integration between them. Today's stored densities are its condensed image (♯), which is why the trace pairs read equal units.
  - Operators: L, S, F, B, … : V → V\*; the sweep is V\* → V; the adjoint is Aᵀ : V → V\*; ψ† ∈ V; r ∈ V\*.
  - "Cotangent = source" (#276 A4) is correct.
  - Stored values may keep the condensed chart, so no number moves; the space says it is V\*.
- **R4, the S_N trace in currents.** The trace unknowns become per-ordinate partial currents J = |Ω·n|·w·A·ψ. The tangential ordinates leave the space, so the trace metric is non-degenerate, as diffusion and collision probability already pose it.
- **R5, the space verbs.** Every `FunctionSpace` owns `member(coordinates)`, the inverse of its chart. Member classes typed above `numerics/` are installed by the builder.

## Refuted candidates

- **Λ is a head operator not given its head.** `[REFUTED 2026-11-10]` by `[M]`: Λ(e⊗p) is 0.37–0.45 from the best (Be)⊗p (`scratch/shortcut_census/probes/`). Λ mixes groups and reads each cell's material. It holds only for one group and one material.
- **Re-declare Λ's ends alone.** `[REFUTED 2026-11-10]` for the design question: the word R·Λ·M then stops composing (spaces must agree), and the exact kernel check scales with groups × cells instead of the head's kernel.
- **Reading (a): a source is a density in its own space with the same measure.** `[REFUTED 2026-11-10]` for the ontology question; the fact it establishes is the condensed chart. The density units are the image of the integrated quantity under ♯. Under (a), 3 of 7 flux/source pairs (the traces) collapse to equal spaces, and the code's transposes (flux in, source out) read as an inversion; under (b) both are correct.
- **A free role tag beside the unit line.** `[REFUTED 2026-11-10]`: the tag would carry no mathematical content. Primal against dual (V vs V\*, related by the measure) is the structure the tag would only name.
- **Keep the S_N trace in angular flux with the #705 quotient check.** Not chosen (R4). The degenerate trace is an artefact of posing the boundary in densities; in currents the tangential directions do not exist.
- **A one-probe kernel check.** Not chosen for #705: its false-pass set is a proper subspace, and a structured leak passes deterministically. The exact check costs k ≤ 4 matvecs on the production carriers (`[M]` `scratch/shortcut_census/probes/probe_bound.log`).

## Open design items (first pass, for review)

- **D1.** How a space declares its side.
  - `V.dual()` exists (ruling 26: the dual carries the inverse Riesz form), and a data-side space IS `V.dual()`.
  - The open part is how the condensed chart is spelled. The dual's members currently pair through `G⁺`, while the stored values are densities.
  - The chart is either a space attribute or a second object, the density representation ♯(V\*).
- **D2.** Member classes follow their space, amending G2.3. The flux classes become V members and the source/sink classes V\* members. `UNITS` becomes per side, with the integrated units on V\* and the density units on its chart.
- **D3.** Operator ends.
  - The 16 role-changing atoms re-declare V → V\*, and the sweep and inverses V\* → V.
  - `FlattenedOperator` takes a domain template and a codomain template.
  - The adjoint door reads the sides, which removes the W/W\* conflation the census found.
- **D4.** The full moment-field space (head ⊗ group ⊗ cells) and its product metric. Λ spelled Σ_ℓ P_ℓ ⊗ K_ℓ, with M and R as axis-addressed factors.
- **D5.** The current-posed S_N trace: the reflection, albedo and inflow operators act on currents, and the boundary word and the windowed path both change.
- **D6.** #705's kernel check lands on top, for the remaining degenerate spaces: rank-deficient moment heads and zero-weight nodes.
- **D7.** #706 (negative weights) needs its own ruling. Level-symmetric quadratures have negative weights at high order.

## Next

A first-pass design document, then a W5 review (cross-domain-attacker and elegance-enforcer), then a user checkpoint on D1–D5 before any code.

## Branches

- `refactor/true-ends`: this campaign, from `main` `a00539f4`.
- `fix/predicted-shortcuts`: the #705 work in progress, `7615c845`. It rebases onto this campaign and lands after it.
