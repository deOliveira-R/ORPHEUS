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
  - The test-architect's spec is `scratch/shortcut_census/gate_spec.md`; the gates are in `7615c845`.
  - It left five items owed:
    1. comparison rows for 9 of the 27 shortcuts: S7, S8, S13, S18, S19, S24–S27;
    2. the same mirror re-pose for `test_coupled_gate2_swap_law_bit_identical`;
    3. the ERR entry's `catches` on `test_d1_a_leaking_operator_has_no_adjoint`, with `test_g2` and `test_p1` as further carriers;
    4. a math `:label:` for the kernel condition, so the d/e rows can carry `verifies`;
    5. K80's candidate loop, re-spelled on `member()`.
- **D7.** #706 (negative weights) needs its own ruling. Level-symmetric quadratures have negative weights at high order.

## W5 review of design v1 (2026-11-11)

Reviews: `scratch/true_ends/review/xda.md` (cross-domain-attacker, probes in `review/xda/`) and `scratch/true_ends/review/elegance.md` (elegance-enforcer, probe `review/elegance_probe_h.py`).

**What holds:**
- §0's arithmetic. `[M]` `G⁺AᵀG` equals the dual map in the density chart to 1e-16, on diagonal, dense and non-square ends.
- D5 (i), the quotient on the state side, is the native realisation. `[M]` On product(4,4) in 2-D the stored A is zero on 160 of 160 tangential rows and columns, and the live block matches to 7e-17. (ii), currents as the state, puts the adjoint's trace source on the function side, against R3, and raises cond from 25 to 60. Specular reflection reads identically in both charts.
- 0 of 18 production `.H` sites need the Hilbert adjoint of a V → V\* operator. `[M]` 15 feed the adjoint SN entries, 2 are same-side Hilbert adjoints, 1 is the alias.

**What breaks:**
- **D1, the chart as a member property.** `[REFUTED 2026-11-11]` The Riesz legs are built from spaces at construction, before any member exists (`operator.py:1660–1672`). Re-declaring the codomain as `V.dual()` (metric G⁺) moves every `.H` by 0.74–15.9 `[M]`. The chart must be space identity.
- **D3, re-declaring the condensed atoms V → V\* and switching to `.dual()`.** `[REFUTED 2026-11-11]` Bare `.dual()` is 59 % off, and `.H` with codomain `V.dual()` is 8 % off `[M]`. The transport adjoint needs the chart legs: A† = κ⁻¹ ∘ (κ ∘ A).dual(), with κ the density-to-integrated isomorphism (applies G).
- **The involution holds only in the integrated chart.** `[M]` `(Aᵀ)ᵀ = A` is exact in the integrated chart; in the density chart it is `PAP`, 0.82 off with a leak present.
- **D4, Λ preserving the head's kernel.** `[REFUTED 2026-11-11]` `[M]` Production Λ fails the kernel condition on 2 of 4 rows: 9.2e-2 on LS4 L=5 (a shipped test row) and 0.16 on product(4,4) L=4; there `Λ.H.H` is off from Λ by 0.38 and 0.56.
  - The cause: on under-resolved rules the Parseval kernel mixes ℓ (‖[P_head, P_ℓ]‖ ≈ 1 on 5 of 11 rows).
  - This also refutes the "0 of 33" claim in #707's body, which covered L ≤ 2 only (`scratch/shortcut_census/probes/probe_lambda_kernel.py`).
  - `[R]` Production adjoint numbers are not affected: the SN adjoint daggers the word R·Λ·M on the angular space, not Λ alone. To be measured.
- **D2, `HarmonicMomentFlux` on V.** `[REFUTED 2026-11-11]` The moment head is covariant (moments are integrals ∫Y ψ dΩ), so it lives in C\* ⊗ V_cells. `[M]` Mixed variance is constructible: `TensorProductSpace.from_factors((X.dual(), Y))`.

**The exact law found:** `[M]` `Mᵀ = W·R·diag(1/(2ℓ+1))` to ≤1.6e-16 on 7 of 7 rows, including the 5 where the Hilbert square fails. R is not independent data.

**Decisions for the user** (the 2026-11-11 checkpoint):
- R-1: store V\* in the integrated chart, or keep the condensed chart with four space objects and κ.
- R-4: two head spaces C and C\* with the continuum ℓ-diagonal map between them, against R-6.2c-1's one Parseval-dressed space.
- R4's realisation: (i) the quotient, with trace rows written in currents.

**Smaller items, decided unless the user objects:**
- A tensor product's ends are the product of its factors' spaces, a new law on `TensorProductOperator`. R2's "ends where members live" wins over "spaces by agreement".
- `FlattenedOperator`'s template retires outright under R5.
- `PredictedKernel.H`'s hand-rolled Riesz legs retire.
- `FissionMomentOperator` joins D4.
- The tensor-law gate uses the production `project_live`, not `pinv`: `[M]` `pinv` reads 7.7e-2 where `project_live` reads 2e-15.

## Next

- The first-pass design: `scratch/true_ends/design_v1.md` (2026-11-11). Its central `[R]` claim: under (b) the transport adjoint is the dual map V → V\*, and today's metric legs are the conversions between the two charts of V\*.
- The W5 review is running: cross-domain-attacker → `scratch/true_ends/review/xda.md`, elegance-enforcer → `scratch/true_ends/review/elegance.md`.
- Then a user checkpoint on D1–D5 before any code.

## Branches

- `refactor/true-ends`: this campaign, from `main` `a00539f4`.
- `fix/predicted-shortcuts`: the #705 work in progress, `7615c845`. It rebases onto this campaign and lands after it.
