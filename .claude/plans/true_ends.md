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

- **R6, the data side stored integrated** (2026-11-11). Sources, currents and responses are integrated amounts. The transport adjoint is the plain dual map, `(Aᵀ)ᵀ = A` holds by algebra, and source numbers change by the measure factor (re-baselined; k and fluxes do not move).
- **R7, the moment head as two spaces** (2026-11-11; amends R-6.2c-1). C holds coefficients (what R reads) and C\* holds moments (what M writes), linked by the continuum ℓ-diagonal map. `HarmonicMomentFlux` lives in C\* ⊗ cells.
- **R4 re-ruled** (2026-11-11). The S_N boundary state is ψ on the non-tangential ordinates (the quotient), with the current view J = ♭ψ through the Riesz verbs. Reflection and albedo are maps on currents. Diffusion and CP follow in the method-harmonisation step.

- **R8, the data side is dual on the SPATIAL factor only** (2026-11-11, Q-A). V\* = angle ⊗ cells\*, per cell and per steradian: the discrete per-ordinate balance. There is no ÷w_n in the sweep, and the angular measure enters only where angles couple (M, R, the pairing).
- **R9, the trace measure includes the face area** (2026-11-11, Q-B). It is |Ω·n|·w_n·A_f, so ♭ψ is a true partial current.
- **R10, inputs are densities, converted once at the input boundary** (2026-11-11, Q-C). Specified as densities, they become the integrated data side by ♭ at the parse edge; the density remains a queried view for output.
- **Q-D, the LD slope-moment measure: RESOLVED, no conflict** (2026-11-11, `scratch/true_ends/ld_moments/memo.md`; six perspectives argued as the user asked).
  - `[M]` The 1-D scan's ×V plus the closure's θ (`_ubld.py:460–480`) give θ·V·q̂, the Galerkin moment ∫q P₁ dV.
  - The single definition: the stored integrated value is the scheme's cell mass (the Gram of the trial basis under the cell measure) applied to the Legendre coefficients, `q_int = M_cell·q_coef`. On a Cartesian cell that is `V·θ^{|o|}·q̂_o` with θ = 1/3, per ordinate and per steradian, with no w_n.
  - `[M]` In an MMS on the slab, the V·θ chart converges at h³ for the mean and h² for the slope; the V chart drops to h² and h.
  - Curvilinear `[M]` (symbolic): the cell Gram is not diagonal (sphere pole cell: M₀₁/V = 1/2, M₁₁/V = 2/5). LD refuses curvilinear today.
- **R8 refined:** the moment axis of a spatial closure (the trial basis inside a cell) belongs to the SPATIAL factor that R8 dualises. `[R]` mathematically forced: the slope index is a spatial basis index. CONFIRMED by the user 2026-11-12, with the instruction that the derivations become part of `orpheus/derivations/` (Branch 1, State 1A: `ld_ubld.py`, foundation gates), so the numerical and symbolic forms can be compared in tests and documentation figures generated programmatically.

- **A consequence of R8, measured during Step A** (2026-11-12, the test-architect; `[M]` a GL4, 3-cell toy with S = RΛM, script in the session scratchpad `ta/r8_probe.py`).
  - Under R8 the transport adjoint is NOT metric-free: the plain transpose of the R8-stored operator is 0.113 (relative) off the true adjoint, while W⁻¹A_intᵀW (angular legs only) matches it to 5.6e-16.
  - So V\* = A ⊗ (E⊗X)\*, with mixed variance: primal in angle, dual on energy ⊗ space ⊗ in-cell moments. It is not `V.dual()`.
  - The transport adjoint is the adjoint with respect to the TRANSPORT PAIRING (G_A ⊗ canonical): angular legs only, spatial legs retired.
  - Its involution is exact by algebra wherever the angular weights are nonzero. The kernel condition re-enters only on zero-weight ordinates.
  - `[REFUTED 2026-11-12]` design v2's "the transport adjoint is plain `.dual()`, and #705's check leaves every flux → source operator": it holds on the spatial factor only. The full dual (w_n stored in sources) would make it metric-free but needs ÷w_n in the sweep, which R8 rejected.

- **R11, the transport pairing is a first-class object that owns the adjoint** (2026-11-12). `DualPairing(V, V\*)` is G_A ⊗ canonical.
  - `pairing(ψ†, q)` evaluates it, and `pairing.adjoint(A)` is the A† with ⟨ψ†, Aψ⟩ = ⟨A†ψ†, ψ⟩.
  - `.dual()` stays the metric-free arrow, and `.H` refuses on V → V\*.
  - The gauge, the Rayleigh estimator and the gate oracle call it.
  - Where the pairing is obtained from is not yet ruled.
- **R12, the adjoint gauge is ⟨F†ψ†, 1⟩ = 1 via the pairing** (2026-11-12). The importance moves by one declared global factor; k† does not move.
- **R13, the Rayleigh estimator pairs over the full space, bulk ⊕ trace,** with the trace measure (R9) (2026-11-12). The trace term is gated in Step B.
- **The Step A scope gains the LD fission defect.** `[M]` (test-architect) `FissionOperator.apply` raises on every LD problem (`transport/material_field.py:483`: the layout gather refuses the in-cell moment axis).

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
- Design v2: `scratch/true_ends/design_v2.md` (2026-11-11).
- The W5 review is DONE (see its section above): cross-domain-attacker → `scratch/true_ends/review/xda.md`, elegance-enforcer → `scratch/true_ends/review/elegance.md`.
- The user approved Step A (2026-11-12). The test-architect is writing the Step A spec, the acceptance baseline and the red gates: `scratch/true_ends/stepA/`. The Q-D derivations LANDED on the branch (`9f754595`): three `derive_*` functions in `ld_ubld.py` and four gates in `tests/gates/transport/spatial/test_ld_ubld_symbolic.py`. Owed to the archivist with the R6 theory text: a math `:label:` for "the integrated chart of an LD moment row", so the MMS gate can carry `verifies`.

## Branches

- `refactor/true-ends`: this campaign, from `main` `a00539f4`.
- `fix/predicted-shortcuts`: the #705 work in progress, `7615c845`. It rebases onto this campaign and lands after it.

## Compaction point 1 (2026-11-12, after the Step A gates)

**Commits on `refactor/true-ends`** (from `main` `a00539f4`):
- `91bda66d` … `73440c9c`: plan commits (rulings R1–R13, the W5 synthesis, design v2, Q-D).
- `9f754595`: the LD Galerkin-moment derivations (`ld_ubld.py`, 4 gates).
- `e8760f6f`: Step A's acceptance baseline (green) and red gates.

The #705 work in progress is `7615c845` on `fix/predicted-shortcuts`.

**The measured red baseline for Step A** (`[M]` `scratch/true_ends/stepA/final_run.log`, the 5 new gate files, `-O`, serial, 44 s): 119 failed, 81 passed, 6 strict xfail, 3 skipped. The acceptance baseline (`tests/gates/cross_method/test_true_ends_acceptance.py`) is green and must stay green.

**Superseding corrections:**
- V\* is NOT `V.dual()`. It is A ⊗ (E⊗X)\*, mixed variance (R8 + the 2026-11-12 measurement).
- The transport adjoint is `pairing.adjoint(A)`, keeping the angular legs (R11). It is not `.dual()`.
- The #705 check is needed only on zero-weight ordinates and the remaining degenerate primal spaces.

**Step A sub-steps** (`[R]`; each one a call-site set, green on its own gates, the acceptance baseline green throughout):
1. **A1: the member seam.** `kernel_basis`, `member()`, the composite overrides and the `zeros_factory` wiring, cherry-picked from `7615c845` without its predicate. Closes NM-1.
2. **A2: the data-side space and `DualPairing`** (`orpheus/numerics/pairing.py`), with their law gates.
3. **A3: the integrated data side.**
   - The seven `*SourceSink` classes move to the data side.
   - The producers integrate (the `/W` home `angular_lift.py:440`; the spatial measure `V·θ^|o|`).
   - The sweeps' consumers: the 1-D walk loses its ×V, and the multi-D kernel works on integrated values.
   - The R10 input doors convert at the edge, and the LD fission fix lands here.
   - Gates: the chart law, LD-F, and the acceptance baseline.
4. **A4: operator ends.** The atoms become V → V\*, the inverses V\* → V, `FlattenedOperator` uses members, and compositions refuse a mismatch. Gates: SA-1 to SA-4, NM-2.
5. **A5: the adjoint.** `pairing.adjoint`, the `.H` refusal, the 16 cross-side sites, the gauge, the Rayleigh estimator, and `PredictedKernel.H`. Gates: SA-5 to SA-8, NM-3.
6. **A6: diffusion ends and k.** Gates: DF-1, DF-2.

**Open:** where the pairing is obtained from (`_true_ends.pairing_of`).

**Durable lessons:**
- Two charts of one element look like two measures until the whole step is read: Q-D was half a step.
- A plain transpose is not the adjoint when an axis is pointwise: the R8 measurement.

**Open after compaction point 1, the pairing's home** (the user's question, 2026-11-12: "Is this a Discrete Measure property? … Is this why the Energy Grid is a counting measure?"). `[R]` the orchestrator's reading, NOT yet ruled; it goes to a W5 attack before the A2 code.
- **The rule.** Every axis carries a discrete measure: w_n in angle, V in space, ΔE in energy. Each SPACE declares per axis whether its coordinates are DENSITY (measure not absorbed) or INTEGRATED (absorbed). A pairing integrates exactly once per axis: it is canonical when exactly one side is integrated, and measure-weighted when both are densities.
- **Per axis:**
  - angle: ψ, q, ψ† and r are all densities, so the pairing carries w_n;
  - space: the flux is a density and the source is integrated (R8), so the pairing is canonical;
  - **energy:** the multigroup flux is integrated (∫_g φ dE) and the cross sections and the importance are densities, so the pairing is canonical. That is WHY the energy grid is a counting measure.
- **Consequence.** R3's "ψ† ∈ V" holds in angle and space only. In energy, ψ† (a density) and ψ (integrated) differ in variance. The code is consistent with this (the adjoint transposes group matrices without ΔE), so no number is wrong; the claim was. R3 needs this refinement.
- **The proposed home.** The data side and the adjoint's state side are derived by flipping variance on the axes the operator integrates, and `V.pairing` / `V.data_side()` read the per-axis variance. To be attacked by the cross-domain-attacker and the elegance-enforcer, the energy case included.

## W5 on the pairing's home (2026-11-12)

Reviews: `scratch/true_ends/review/variance_xda.md` (cross-domain-attacker) and `scratch/true_ends/review/variance_elegance.md` (elegance-enforcer), with their probes beside them.

**Agreed by both reviewers:**
- **The arithmetic of the per-axis reading stands; its structure does not.**
- **[REFUTED 2026-11-12] ΔE as the energy measure.**
  - `EnergyAxis` refuses weights (`axis.py:678`).
  - Two density views compete, ΔE and Δu, and `[M]` 4 of 4 production readers use them only as views.
  - No pairing or adjoint reads either view. Making ΔE the metric would make an adjoint wrong by a factor of 2.5e5 spellable (`[M]` xda p1).
  - The energy grid is a counting measure because energy has no preferred reference measure. Angle has dΩ and space has dV; energy has two candidates and the physics picks neither. Only the group masses are free of a reference measure.
- **[REFUTED 2026-11-12] Variance as a property of `DiscreteMeasure`, the quadrature or the energy grid.**
  - The LD cell Gram and the harmonic-head Gram are not measures.
  - `DiscreteMeasure` has no chart field (`[M]` xda).
  - Variance is a property of a space's FACTOR (X or X\*), not of the axis measure.
- **[REFUTED 2026-11-12] A per-axis density/integrated flag** (elegance: anti-pattern 3 plus Pattern 2). It re-spells `DualSpace`.
- **[REFUTED 2026-11-12] "Angle is all densities".** R7's moments C\* and R9's currents are angle-integrated.
- **[REFUTED 2026-11-12] "Flip on the axes the operator integrates".**
  - S is local in space but flips there, and it integrates over angle but does not flip there.
  - Every term of one balance shares one codomain, so the data side is the METHOD's balance axes, declared once. For SN that is {space}.
- **[REFUTED 2026-11-12] An argument-free `V.pairing` or `V.data_side()`.** One V pairs with both the bulk V\* and the trace.
- **The home: a free, frozen `DualPairing(state, data)`, built by the binder.**
  - Its legs are derived per factor from the types:
    - (X, X\*): a plain contraction;
    - (X, X): the leg G_X;
    - (X\*, X\*): the leg G⁺ according to xda, and refused according to elegance.
  - It adds no new arithmetic: ℓ is a `FactoredMetric` with `None` on the contracted factors.
  - `adjoint(A) = ℓ'⁻¹ ∘ A.dual() ∘ ℓ`.
  - The Riesz legs need a constructor that takes a metric (`operator.py:1381`).

**The open question both raise, for a user ruling: where the variance index lives.**
- **(i) In the axes record, through a new `Axis.dual()`** (xda). V\* = `of_axes(A, E.dual(), X.dual())` stays axis-built, so `.axis("energy")` keeps working.
  - `[M]` today `V.dual()` and the R8 product both raise TypeError on `.axis("energy")`.
  - The A3 source classes read the axes (`_bases.py:458/517/551/625`, `angular_source_sink.py:139`).
- **(ii) As per-factor `DualSpace` products** (elegance), with a derived "factor profile" accessor.
  - `[M]` the product loses its `axes` record, and by reading 7 sites in 4 files break.
  - The elegance review measured the same loss as xda, so (ii) needs (i)'s fix anyway. `[R]` (i) subsumes (ii).

**Corrections to the orchestrator's reading:**
- **"ψ† ∈ V fails in energy"** (from the cross-domain-attacker, a reasoned claim, not mine):
  - The INDEX (upper or lower with respect to the axis measure) of ψ† is the same as ψ's on every axis, so R3 "ψ† ∈ V" holds as an index statement.
  - What differs is TYPE (Radon–Nikodym: a measure is summed on merge, a function is averaged with a weight). Its one consumer is condensation (`mixture.py:435–620` picks sum or average per index by hand).
  - `[M]` on the energy space, adjoint = transpose bit for bit. No number moves.
- **TYPE is a member-class property** (xda): a function/measure kind on the member role, beside {Flux, SourceSink}. This touches design v2's "ψ† is flux-classed": ψ† is an importance, a function.

**For the archivist:**
- `axis.py:642–656` calls the flux space V\* and the cross sections V (the reverse of R3).
- Its collapse hook "sum on V, average on V\*" should read "sum on measures, average on functions".

**For the test-architect** (xda's first tests):
- The data side keeps `.axis("energy")`, with the R8 product's metric.
- `condense` of an importance member is pair-weighted, while the same values typed as flux are summed.
- Pairing two lowers canonically is refused (the carrier is `diffusion/solver.py:431`).
- A4 needs a gate pinning the codomain's factor profile: SA-1 checks only `codomain != domain` (elegance).

**A1 findings** (the test-architect, `tests/gates/numerics/test_member_seam_laws.py` with 30 rows, `tests/gates/cross_method/test_member_seam_builders.py` with 19 rows):
- `FullFieldSpace._make_dual` drops `zeros_factory`, so `V.dual().member` raises today. To fix in A2.
- The coupled space's dual zero is minted from `*SourceSink` fields on the PRIMAL blocks. To fix in A3; the strict xfail mb4 pins it.

## The form-degree attack (2026-11-12)

The user asked whether a missing differential-geometry concept is welded in. The orchestrator proposed form degree from discrete exterior calculus as the one label behind variance and type. The cross-domain-attacker attacked it: `scratch/true_ends/review/form_degree_xda.md`, with its probe `form_degree_xda/p1_streaming_vs_variance.py`.

**The terminology.** "Riesz map" here means the Riesz representation map, not the spectral Riesz projector; `[M]` there are 0 spectral-projector uses in `orpheus/`. In finite dimensions it is the musical isomorphism ♭/♯. These are one object, not two theories.

**[REFUTED 2026-11-12] "Form degree is ONE label that unifies variance and type", for the question "is there a single label?".**
- Read consistently, it is two numbers per axis: the NATIVE degree, what the quantity is (= TYPE), and the STORED degree, which side of ⋆ the coordinates sit on (= INDEX).
- The FACT this establishes: (TYPE, INDEX) is the right pair, and DEC renames it rather than reducing it.

**[REFUTED 2026-11-12] "Streaming is metric-free because it is the coboundary d".**
- `[M]` collision C = diag(Vσ) has the same canonical-transpose identity (1.4e-16).
- `[M]` streaming differs by 0.878 between slab and sphere while d is unchanged.
- L = d ∘ (upwind contraction), and that contraction is not a ⋆.

**Derivations.**
- Degree plus ⋆ derive 5 of 9 rows: R9, R10, R4, the energy counting measure, and the pairing legs. The same rows also derive from (TYPE, INDEX).
- R6 and R8 need a "balance axes" rule beyond degree.
- Condensation half-derives: summing top-degree quantities follows, but a 0-form's weight is extra data in both frames (volume for flux, V·φ for σ, production for χ, φ for ψ†).

**Break points:**
- homogenisation (three quantities stored at degree 0 condense with three different weights);
- the harmonic head (no cell complex);
- diffusion (J⁺, J⁻) (they swap under orientation, where a cochain would change sign);
- LD (it fits finite-element exterior calculus, a Galerkin ⋆, not DEC).

**Where it fits:**
- curvilinear redistribution is d on the μ-mesh;
- R4, #705 and #706 become one condition: ⋆ is degenerate or indefinite.

**The one gain, `[R]`, no number moves.** Under the stored/native reading, ψ ∈ A ⊗ E\* ⊗ X and ψ† ∈ A ⊗ E ⊗ X. This reconciles R8 ("dual on the spatial factor only") with its 2026-11-12 correction ("A ⊗ (E⊗X)\*"): they name different reference spaces. The attacker's earlier "the index is the same on energy" was a convention, not a measurement.

**The minimal realisation (the attacker's reading).** `Axis.dual()` is the stored degree, the member kind is the native degree, and `DualPairing` is the contraction. `Axis.of_degree(k)` has consumers only for k ∈ {0, top}, so it reduces to `dual()`.

## A1 review fixes (2026-11-12)

qa (`scratch/true_ends/stepA/a1_qa.md`) and the elegance-enforcer (`scratch/true_ends/stepA/a1_elegance.md`) reviewed `d482b89a`.

**Fixed in-session, each with a witness that reddens on the old body:**
- **One split law.** `HilbertMetric._block_split` is abstract, and `kernel_basis` derives from it in the base; the `_SplitMetric` cast protocol retired.
- **`DiagonalMetric` broadcasts leading-aligned throughout.** `_block_split` AND `apply_block` were trailing (numpy), while `rank` and `validate_for` were leading. A factored space with a two-axis diagonal block was built with a rank and then could not apply its own metric. Witness: k1 "factored two-axis diagonal block".
- **`DenseMetric` kernel vectors are eigenvectors of G**, with the count still read from `rank`. The eigenvectors of I − G⁺G lose digits as ε/λ_min. Witness: k1 "dense near the cutoff".
- **`FullFieldSpace.coordinates_of` restores the identity leg.** A dual accepted a primal composite. Witness: MB5 (4 rows).
- **The two `member()` twins retired onto `space.member`:** `posing._carrier` and `MatrixInverseOperator._returned_as`.
  - Both read `getattr(space, "zeros", None)`, which changed behaviour silently when A1 gave `FullFieldSpace` a `zeros()`. That was the D6/D8 cause.
  - The test pinning the old refusal migrated: `test_typed_operand_without_a_carrier_space_raises`.
  - `posing` refuses a pencil whose loss declares no domain.
- **Stale docstrings:** `flat_operator.py` (what `apply` accepts) and the builders gate's `is`/`==` note.

**Routed to later sub-steps (same branch, recorded so they are not lost):**
- **A3:**
  - The `CoupledSpace` identity leg (qa F1). It refuses the dual's own members until the coupled dual zero is minted on the dual blocks (MB4, strict xfail).
  - `member` is not a function of the space (qa F3): `zeros_factory` is `compare=False`, so equal spaces mint different classes.
  - `FullFieldSpace._make_dual` drops `zeros_factory`.
  - The member class follows the Python subclass, not the variance: `V.dual()` of a tensor product mints a bare `Member` (elegance F2).
  - Five spellings of "mint from flat over a template".
- **A4:** the template retires from `FlattenedOperator` (NM-2). `_flat_coordinates` is the interim. Its `space is None` arm serves only test doubles (0 of 5 production constructions).
- **A6:** diffusion's trace and bulk are stored in the density chart, against R6 (`augmented_mesh.py:148`, `operators.py:74`).
- **The test-architect:**
  - an SN fixture with nonzero predicted kernel columns (3 of 3 probed have none, so the SN annihilation check never runs);
  - a homogenisation gate (the stored degree must not pick the weight);
  - the collision-transpose row as a gate on `pairing.adjoint`;
  - an A4 gate pinning the codomain's factor profile.
- **The archivist:**
  - `axis.py:642–656` names the flux space V\* and the cross sections V, the reverse of R3. Its collapse hook should read "sum on measures, average on functions".
  - Name the reference spaces in R8 and in its correction.

## Compaction point 2 (2026-11-12, after A1)

**Commits on `refactor/true-ends`, since compaction point 1:**
- `e009fb9b`, `0decff33`: plan commits (the per-axis variance question; W5 on the pairing's home).
- `d482b89a`: A1, the member seam.
- `d836cbaf`: the A1 review fixes.
- This commit: the plan (the form-degree attack, the A1 fixes, this point) and agent memory.

**The measured red baseline** (`[M]` `scratch/true_ends/stepA/a1_fix_full.log`, `-O`, serial, `-m "not slow"`, all of `tests/gates`, 1 h 52 min): 118 failed, 20 358 passed, 72 xfailed.
- All 118 are Step A red gates: 113 in `sn/test_true_ends_step_a.py` and `cross_method/test_source_chart_law.py`, 2 in `diffusion/test_true_ends_step_a.py`, and 3 in `numerics/test_true_ends_member_seam.py` (NM-2 ×2, NM-3).
- NM-1 is green. The acceptance baseline is green (23 of 23).

**A2 is BLOCKED on a user ruling. The user, 2026-11-12: "don't rule yet on the DualSpace and related questions"; the newer model (Fable 5.2) is to be consulted after a harness restart and a compaction.** The questions, stated neutrally, with their evidence:
1. **Where the variance (X versus X\*) of each factor is recorded.**
   - (i) `Axis.dual()` in the axes record, so V\* = `of_axes(A, E.dual(), X.dual())` keeps `.axis()`.
   - (ii) Per-factor `DualSpace` products, with a derived factor profile.
   - `[M]` today `V.dual()` and the R8 product both lose `.axis("energy")`.
   - Evidence: `review/variance_xda.md`, `review/variance_elegance.md`.
   - The user asked whether `DualPairing` should inspire the `DualSpace` design ("would adjoint properties live in DualSpace and ♭/♯ be mappings between these?"). The orchestrator's explanation is in the session transcript:
     - the pairing is metric-free; ♭/♯ belong to the metric; the adjoint belongs to a PAIR of spaces;
     - the partial ♭ (I_A ⊗ ♭_X) is R6's integration and R10's door.
2. **`DualPairing` on (X\*, X\*).** Both reviewers refuse the canonical contraction (the diffusion double integration, `diffusion/solver.py:431`). xda would apply the G⁺ leg; elegance would refuse it and leave V\*'s own inner product on `DualSpace`. The user asked what this protects against; the orchestrator's answer: two integrated quantities contracted canonically count the measure twice, and the residual norm in V\* is the one legitimate G⁺ use.
3. **The measure/function TYPE as a member-class kind.** It would land in A3, with ψ† as an importance, a function. The user leaned to "Yes, in A3", noting that it is related to question 1. Not ruled.
4. **The form-degree reading** (see "The form-degree attack" above). It is a renaming of (TYPE, INDEX), not a single label; its gain is that ψ ∈ A⊗E\*⊗X and ψ† ∈ A⊗E⊗X. The attacker asks whether the flux's energy factor becomes E\*. Not ruled.

The orchestrator's draft rulings R14–R16 (in the transcript, 2026-11-12) are PROPOSALS only and must not be read as ruled:
- R14: variance on the axis;
- R15: per-factor pairing legs, refusing (a\*, a\*);
- R16: TYPE is a member kind in A3.

**Next, after the consult and the rulings:**
1. A2: the data-side space and `DualPairing` (`orpheus/numerics/pairing.py`), with its law gates from the test-architect first.
2. Then A3–A6 as listed at compaction point 1.
3. The A2/A3/A4/A6 items routed in "A1 review fixes" above.

**Durable lessons:**
- **`git checkout <hash> -- <paths>` STAGES the files**, so the next `git commit` sweeps them in. One plan commit swallowed the A1 production code; it was unpushed, and a soft reset plus restaging fixed it. Read `git status --porcelain` before EVERY commit, plan commits included.
- **A duck-typed `getattr(space, "zeros", None)` changes behaviour silently** when a verb lands on a class it reads. Two twins took the typed branch the moment `FullFieldSpace.zeros` landed (retirement-audit item 6). A new verb owes a grep for its name as a STRING.
- **A gate added for one defect can find its sibling.** The k1 two-axis row, written for `_block_split`, reddened on `apply_block` too: one convention had two homes.

## The consult on the four open questions (2026-11-12, after compaction point 2)

An advisor was consulted. The harness may still have been serving Fable 5.1 rather than 5.2; this is reported to the user. The advisor proposed an order (4, 1, 3, 2) and three censuses plus a probe, which were run before any ruling. Reports:
- `scratch/true_ends/consult/four_space_probe.md` (numerics-investigator; probes and logs beside it);
- `consult/variance_census.md` (explorer; probes in `consult/variance_census/`);
- `consult/type_census.md` (explorer).
Everything below is a PROPOSAL; none of it is ruled.

**Q4, the four-space square, holds `[M]`** (toy: 4 groups with upscatter, a 4-cell non-uniform slab, GL4, step scheme, vacuum; merges {012|3} non-uniform and {01|23} uniform).
- flux ψ ∈ A⊗E\*⊗X; source q ∈ A⊗E\*⊗X\*; importance ψ† ∈ A⊗E⊗X; response R ∈ A⊗E⊗X\*. The pairings are ⟨ψ†, q⟩ and ⟨R, ψ⟩, each contracting E and X canonically, with G_A on angle.
- The discriminating instrument is the coarse balance under the non-uniform group merge. With ψ and q summed it closes to 1e-15; with them averaged it misses by 1.19 (slab) and 10.3 (0-D). The uniform merge cannot tell them apart (control).
- Wrong assignments break the pairings: ψ† summed moves ⟨ψ†,q⟩ by 2.09; R summed moves ⟨R,ψ⟩ by 0.298 (non-uniform merge).
- R8 and its 2026-11-12 correction are the same square read from ψ's side and from ψ†'s side. The R8 control reproduces (0.11331, 5.6e-16).
- Consequence: R3's two spaces become four (ψ† ∉ V, r ∉ V\*).
- `[REFUTED 2026-11-12]` "ψ†'s index equals ψ's on energy" (variance_xda): under the merge ψ is summed and ψ† is averaged.

**Q3, TYPE is not the index `[M]`** (both the probe and the census agree; this refutes the advisor's guess that one would determine the other).
- Measure versus function decouples from the stored index on space: ψ is a measure stored as a cell average; R is a function stored cell-integrated (summing R is off by 2.8e-2). Energy is the only axis where the two coincide, because it has no reference measure.
- The weight of an average is the quantity's PAIRING PARTNER: ψ† averaged with q preserves ⟨ψ†,q⟩, and averaged with ψ (production) it preserves ⟨ψ†,ψ⟩. With a full-support source each choice moves the other pairing by 2–3e-2.
- Census: 29 production merge rows; today 10 consistent, 7 exceptions, 12 indeterminate. TYPE belongs to each tensor slot, not each member (the scattering sink slot is summed, its source slot averaged), and six averaging weights appear.
- Production condensation is consistent with the square but chooses per channel by hand; `Mixture.condense` matches the toy bit for bit in 0-D.
- Open: whether production's ψ-weighted importance merge (`solution.py:826`, `:292`) moves a number on any multi-cell path.

**Q1, the variance home `[M]`.**
- 6 of 12 minted production factor kinds are `Axis` objects (angle, energy, cells, the harmonic head C, the LD in-cell moment, the condensation coarse-group head).
- The other 6 are axes-less trace and boundary leaves (`AngularTraceSpace`, `AngularFaceTraceSpace`, the per-face current space, `ScalarTraceSpace`, and the two radial-characteristic spaces).
- `.axis()` resolves on 6 of 16 built spaces.
- So `Axis.dual()` alone cannot hold every factor's variance. The bit belongs to whatever generates the factor. The advisor's candidate single verb is a partial dual by axis name, `V.dual(on=…)`, with the full dual as the case naming every axis.
- Also found: the harmonic head already carries C\*'s metric G⁻¹ under a primal type (read on the diagonal arm only), and the condensation head is labelled `"cells"` with a plain `Axis` for energy.

**Q2, the (X\*, X\*) leg `[M]`.**
- Σ_g q_g² is not invariant: it doubles under an equal-mass merge and falls as 1/n_g under refinement. Σq²/ΔE and Σq²/Δu converge, to different limits (0.0891 and 0.0519). A canonical G⁺ leg on E\* is therefore grid-dependent.
- Today 1 of 14 data-side sites pairs two integrated quantities: `cp/solver.py:862`, a legitimate residual norm. 0 of 14 are double integrations; the suspect `diffusion/solver.py:431` is a correct single integration today.
- `[R]` 6 sites become double integrations once A3/A4 store the data side integrated: diffusion `:431`, `:644`, `:645`; SN `_balance_projection`; `loss_kernel_gauge.py:928`.
- Proposal: `DualPairing` refuses (X\*, X\*); a residual norm on V\* is V\*'s own inner product on `DualSpace`.

**Owed, whatever is ruled:**
- test-architect: the non-uniform coarse-balance gate with its uniform control; the S/F/L witness table as the coverage map for `pairing.adjoint`; a swap probe (every source scaled by V) before A3; an A4 gate on `_balance_weights`.
- archivist: `axis.py:642` should read "sum on E\*, average on E"; the condensation head's `"cells"` label.

**Rulings on the consult (the user, 2026-11-12):**
- **R14, four spaces.** flux A⊗E\*⊗X, source A⊗E\*⊗X\*, importance A⊗E⊗X, response A⊗E⊗X\*, with the pairings ⟨importance, source⟩ and ⟨response, flux⟩. It amends R3's "ψ† ∈ V". The user: "we made a sophisticated space machinery, so we should use it properly." The user also asked whether every operator between these spaces has been mapped; an operator-map census was dispatched (`consult/operator_map.md`).
- **R15, `DualPairing` refuses a pairing of two integrated (X\*, X\*) factors.** A residual norm on V\* goes through V\*'s own inner product.
- **R16, measure or function is a per-slot kind; an average's weight is derived from the pairing partner.** It lands in A3, and condensation stops choosing per channel by hand.
- **Q1, the variance home: NOT ruled.** The user asked what `DualSpace` is for now; answered in the session.

**The operator map onto R14** (`scratch/true_ends/consult/operator_map.md`, explorer, 2026-11-12). `[M]` A runtime `__subclasses__` walk plus `ast` found 80 classes; 64 are concrete; CP, MoC and MC have 0 (control: 60 found elsewhere). The verdicts are partly `[R]` from reading.
- **22 of 64 FIT.** No class that changes role fits.
- **37 of 64 are MIS-DECLARED.**
  - 26 have the right axes but mark no factor as integrated: every L, C, S, F, N2N and B, the frames, and the harmonic M and R. 0 of 16 role-changing atoms declares a `DualSpace` end.
  - 11 have wrong axes or a wrong end rule: Λ, the fission moment operator, the SN boundary classes, the radial reconstruction, `InverseMetricOperator` (0 production mints), the Riesz legs, `AdjointOperator`, `_DualOperator` and `TensorProductOperator`.
- **5 of 64 do not fit**: `FlattenedOperator` (one template for both ends), and `Diagonal`, `RankOne`, `Identity`, `Zero` (no declared ends).
- **Maps the square needs that no operator provides:**
  - ♭_X on a product that keeps A and E;
  - ♭_A and ♯_A as the legs of the transport adjoint;
  - a rule naming each end's pairing partner (the swap and the full dual both fail);
  - ♯ on the harmonic head (C\* → C);
  - per-ordinate trace currents with face area;
  - condensation as a typed operator.
- **Correctly absent:** there is no ♭ on energy (`axis.py:680` refuses weights), and no map between flux and importance. But the importance is returned as a `ScalarFlux`: 0 importance or response member classes exist, against 7 flux and 7 source classes.
- **The 0-D homogeneous problem collapses to two spaces.** Its X flip is invisible numerically there.
- **R17, importance and response member classes, minted in A3** (the user, 2026-11-12). They follow their spaces (R5), so a flux can no longer be handed where an importance is meant. Today `AdjointSolution.importance` returns a `ScalarFlux`.
- **Q1: the user asked for a worked example before ruling** (the source space and the transport adjoint). It is given in the session.

**Q1, the worked example given to the user (2026-11-12).** This is the proposed **R18**, NOT ruled. The user replied "Save this and prepare for context compaction. I have a plan."
- **Factors:**
  - A (ordinates, G_A = diag(w));
  - E (groups, no metric);
  - X (cells, G_X = diag(V)).
- **Spelling:**
  - `flux = FunctionSpace.of_axes(A, E.dual(), X)`;
  - `source = flux.dual(on="cells")`;
  - `importance = flux.dual(on="energy")`;
  - `response = importance.dual(on="cells")`;
  - `source.axis("energy")` returns `E.dual()`.
  - Collision declares flux → source (×σ_t·V, it flips X only); the sweep declares source → flux.
- **Why the flip is recorded per factor.** ⟨ψ†, q⟩ = Σ_n w_n Σ_g Σ_i ψ†q puts the G_A leg on (A, A) and a plain contraction on (E, E\*) and on (X, X\*). `DualPairing` reads the legs from the factors; a whole-space "is a dual" tag cannot distinguish them.
- **Why the full dual stays.** The transport adjoint is C† = ♯_A ∘ Cᵀ ∘ ♭_A:
  1. ♭_A, a partial flip on angle, maps importance A⊗E⊗X to A\*⊗E⊗X, the full dual of source.
  2. `C.dual()`, the metric-free transpose, maps that to A\*⊗E⊗X\*, the full dual of flux.
  3. ♯_A maps it to A⊗E⊗X\*, the response space.
  - This is the W⁻¹CᵀW measured at 5.6e-16 against the true adjoint. The transpose uses the full flip (today's `V.dual()`); the Riesz legs use partial flips.
- **The alternative** (per-factor `DualSpace` wrappers): `TensorProductSpace(A, DualSpace(E), DualSpace(X))` loses the axes record. `.axis("energy")` raises `[M]`, and 7 sites in 4 files rely on it. The arithmetic is the same either way.
- **The proposal's declared edge.** The axes-less trace leaves flip whole, angle included, while the bulk rule leaves angle unflipped. They flip per factor only once Step B rebuilds them as face ⊗ angle (R4, R9).

## Compaction point 3 (2026-11-13, after the consult)

**Commits since compaction point 2**, all plan-only, `[skip ci]`:
- `b396f974`: the consult results, plus numerics-investigator memory.
- `5dddef57`: R14–R16.
- `38f960c6`: the operator map.
- `74ddee97`: R17.
- This commit: the Q1 worked example and this point.

**Ruled:**
- R14, four spaces;
- R15, `DualPairing` refuses (X\*, X\*);
- R16, a per-slot measure/function kind, with an average's weight taken from its pairing partner;
- R17, importance and response member classes in A3.

**Not ruled:** Q1, the variance home (the proposal is R18 above). The user has "a plan" to bring after compaction. Wait for it; do not rule R18 by default.

**Reports** (all under `scratch/true_ends/consult/`; untracked scratch):
- `four_space_probe.md` and its probes;
- `variance_census.md`;
- `type_census.md`;
- `operator_map.md`.

**Owed, whatever Q1 becomes:**
- **test-architect:**
  - the non-uniform coarse-balance gate, with the uniform merge as control;
  - the S/F/L witness table as the coverage map for `pairing.adjoint`;
  - a swap probe (every source scaled by V) before A3;
  - an A4 gate on `_balance_weights`;
  - the census spy rerun at HEAD, to turn the operator map's reading-based rows into `[M]`.
- **archivist:**
  - `axis.py:642` should read "sum on E\*, average on E";
  - the condensation head's `"cells"` label;
  - the R8 reference spaces named.
- **Open measurement:** whether production's ψ-weighted importance merge (`solution.py:826`, `:292`) moves a number on any multi-cell path.

**The red baseline is unchanged from compaction point 2**: 118 Step A red gates. No production code changed since `d836cbaf`.

**Next, after the user's plan and the Q1 ruling:** A2 (`orpheus/numerics/pairing.py`: the data-side spaces and `DualPairing`, with the test-architect's gates first), then A3–A6 and Steps B–D.

**Durable lessons:**
- An adjoint-versus-transpose number cannot tell E from E\*, because G_E = I. The discriminating instrument for variance on a metric-free axis is a conservation law under a NON-UNIFORM merge; a uniform merge is blind to it.
- A census can refute the advisor's hypothesis: TYPE is not the index. Measure before ruling, even on a confident structural guess.
