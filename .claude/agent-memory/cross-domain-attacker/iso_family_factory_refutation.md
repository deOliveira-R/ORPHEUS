---
name: iso-family-factory-refutation
description: VERDICT on the proposed "Factory object returning the appropriate iso-channel operator keyed on the problem" (IsotropicScattering/IsotropicN2N/IsotropicFission, CS4c step 4, 2026-08-30). REFUTED as a DISPATCHER, three independent ways, and REDIRECTED. (1) FIBER NON-SINGLETON — over one FullFieldSpace the scattering fiber holds ScatteringOperator + S.isotropic_energy + S.foldable_part() + S.residual_part(); the choice is by consumer ROLE, not by space, so no section exists and the factory's key would be a `role=` string = the posing concept CS4c §5 defers. (2) NO DISCRIMINATION EXISTS TO COLLAPSE — [M] 4 construction sites tree-wide, ZERO runtime branches on space kind/channel tag; a factory MANUFACTURES the conditional and then names it the solution. (3) RANK REFUTES THE TRIO — S0/N2N are one body differing by (matrix_of, scale) via MaterialField._accumulate_contracted; F is rank-1 |χ⟩⟨νΣf| with a FACTOR-SWAP transpose that violates FissionKernel's own simplex guard ⟹ the set is 2+1, not 3. Both stated justifications refuted on the tree: DSA consumes ZERO channel operators ([M] 0 hits in sn/acceleration/; it reads diag(K_ℓ) as a COEFFICIENT via mat_xs.foldable_sigma/residual_sig_s), and a "sensor" is a FUNCTIONAL (codomain ℝ) = the frame's ℓ-row ∘ reaction weight, already composed at dsa.py:584. THE REDIRECT (inverse direction): the honest object is the LIFT FUNCTOR energy-binding → angular-binding (_per_ordinate.py:20-22 states its own defer-until-2 trigger; the third channel HAS arrived and hand-inlines it at fission.py:615-625 with a 4th spelling of W) over ONE collapsed IsotropicEnergyOperator whose channel lives in its FIELD's type. ⭐⭐ The reframing datum: S.foldable_part()/residual_part() ALREADY implement operator-tier "two flavours of one channel", gated at rtol=1e-14, with [M] ZERO production callers — so the premise is implemented twice over and the real question is consumer-side. Bonus doc bug: material_xs_field.py:760,780 claim the operator consumes those accessors; since step 3 it reads the kernel field — present-tense-false, only live consumer is DSA.
metadata:
  type: project
---

# The iso-family "Factory" — refuted as a dispatcher, redirected as a functor

Probe memo (full evidence, 10 discriminating questions, 13-row rejected-candidate
table): `scratch/cs4c_step4_factory_probe.md`. Read on `main` at the CS4c step-3
merge tip; **Nexus graph stale at `92dcc30f`** — every fact below is a direct file
read (L-005).

## The three independent refutations of the DISPATCHER reading

1. **No section exists (fibration frame).** Model the factory as a section of
   `Op → Space`. `[M]` the fiber over one `FullFieldSpace` contains four
   legitimate members of the *scattering channel alone*: `ScatteringOperator`,
   `S.isotropic_energy` (`scattering.py:555`), `S.foldable_part()` (`:960`),
   `S.residual_part()` (`:988`). The choice is by the consumer's ROLE in the
   algebra, which no space carries ⟹ the key must be `Space × Role`, and `Role`
   is the posing concept CS4c §5 explicitly defers. (And `posing` is a word
   already spent on the `(A_loss, M)` + eigenvalue-role axis of the same solve —
   L-012's spent-word kill.)
2. **The discrimination does not exist yet.** `[M]` 4 construction sites
   (`diffusion/solver.py:242`, `homogeneous/solver.py:197`, `scattering.py:555`,
   `n2n.py:158`) and **zero** runtime branches on space kind or channel tag.
   Every site knows statically which binding it wants. A factory converts a
   statically-known fact into a runtime dispatch — illegal-states-unrepresentable
   run backwards. `discriminations` says *a repeated conditional is a missing
   type*; there is no repeated conditional, so the factory creates the
   conditional and then names it the fix.
3. **`space.bind(kernel)` dies on import direction.** `numerics → transport` is
   a forbidden edge; a data object's verbs return ARRAYS, only binders return
   OPERATORS (L-020). This is why `HarmonicFrame.for_space` is transport-side.

## Both justifications, checked against the tree

- **"SN+DSA needs both the aniso and the iso version of the same channel."**
  `[M]` DSA consumes **zero channel operators** (0 hits for
  `IsotropicScattering|IsotropicN2N` in `orpheus/sn/acceleration/`). It reads
  `mat_xs.foldable_sigma()` / `residual_sig_s()` (`dsa.py:263,266`) — the
  **within-group DIAGONALS** `σ_s0^{g→g}`, `σ_s1^{g→g}` — as per-cell SCALAR
  coefficients forming `σ̂_R = σ_t − σ_s0^{g→g}` and `D = 1/[3(σ_t − σ_s1)]`.
  A diagonal is a **coefficient field**, not a member of the transfer-kernel
  family, so no factory over that family can produce it.
  ⟹ **and the co-sourcing half is already true, by a STRONGER mechanism:**
  `S.isotropic_energy` is a `cached_property` binding of `S`'s OWN
  `ScatteringMaterialField`. A factory co-sources by CONVENTION (same inputs, two
  calls); a satellite co-sources by CONSTRUCTION (one datum instance).
- **"Sensors at different levels want per-level flavours."** A sensor's codomain
  is `ℝ` ⟹ it is a FUNCTIONAL, not in the operator family. A per-level sensor is
  `⟨w_ℓ ⊗ Σ_x, ·⟩` = the frame's ℓ-th analysis row ∘ a reaction weight, and
  `[M]` DSA already composes exactly that at `dsa.py:584-586`. Putting it in the
  operator factory gives the factory a non-uniform codomain ⟹ switchboard, not
  family.

## What IS real (the redirect — the inverse direction)

- **RANK splits the trio 2+1.** `IsotropicScattering` and `IsotropicN2N` are one
  arithmetic body differing by `(matrix_of, scale)`:
  `MaterialField._accumulate_contracted(..., lambda k: k.p0, spec=_FORWARD)` vs
  `(..., lambda k: k.matrix, spec=_FORWARD, scale=ν₂ₙ)`
  (`material_field.py:137,223,323`) ⟹ collapse to ONE
  `IsotropicEnergyOperator(field: MaterialField[...])`, channel in the FIELD's
  TYPE. `IsotropicFission` cannot join: `F = |χ⟩⟨νΣf|` is rank 1, its transpose
  is a FACTOR SWAP (free) not a `_TRANSPOSE` einsum, joining costs `O(G²)` where
  `O(G)` suffices, and the swap **violates `FissionKernel`'s own simplex guard**.
  ⭐ And `IsotropicFission` is not a new object at all — `FissionOperator.kernel`
  (`fission.py:281`, `TensorProductOperator` over `RankOneOperator`) already IS
  the energy operator; both scalar arms delegate to it (`:690`, `:707`). Step 4
  is a NAMING+BINDING act, not a construction act.
- **The LIFT is a FUNCTOR, and its defer-until-2 trigger HAS FIRED.**
  `_per_ordinate.py:20-22` states verbatim *"When a third isotropic lifted
  channel appears, this is the seed of the generic lift operator"* — and `[M]`
  `FissionOperator`'s composite arm (`fission.py:615-625`) hand-inlines
  `(iso/W) + zeros` instead of calling it. `K ↦ (1/W)·E₀∘K∘M₀` =
  `angular_frame(0).conjugate(K)` is a functor energy-binding → angular-binding;
  `N2NOperator` IS `lift(IsotropicN2N)` wearing a class. **This is the inverse of
  the proposal**: one verb GENERATING a channel's sibling, not one dispatcher
  CHOOSING between siblings.
- **`W = Σw` has FOUR derivations, three distinct sources:**
  `flux_analysis.frame.measure.weights.sum()` (`scattering.py:506`);
  `axes[0].generator_as(Quadrature).weights.sum()` (`n2n.py:163-172`);
  `ang.shape[0] if ang.weights is None else sum(ang.weights)`
  (`fission.py:618-622` — a fallback resting on an invariant it does not name);
  `frame.discrete_gram[0,0]` (`numerics/frame.py:766`). Smell #16 shape 2. The
  ℓ=0 frame is the object that should own it. Any factory is a **fifth** entry
  point into the blessed chain, by design.
- **`foldable_sigma` exists twice** — facade (`material_xs_field.py:820`) and
  operator (`scattering.py:1028`), both `np.diag(p0)` — and DSA reads the facade.
  Native fix (multigrid pollination): `ScatteringKernel.within_group() ->
  (L+1, ng)` — the Galerkin restriction `e_gᵀ K e_g` onto the single-group
  subspace. Collapses 2 facade accessors + 1 twin into 1 kernel verb and makes
  DSA's `if scattering_order >= 1` a slice.

## ⭐⭐ The datum that reframes the request

`[M]` `foldable_part`/`residual_part` have **ZERO production call sites** (15
lines in `orpheus/`: 2 definitions, 1 sibling ref, 1 comment, 11 docstring
xrefs; every invocation is in `tests/`). They already implement operator-tier
"two flavours of one channel" in exactly the shape the proposal wants, gated at
`rtol=1e-14`, and **nothing calls them** — not the σ_r-fold removal form (#196,
`streaming.py:185`) and not DSA.

⟹ the premise *"we probably need both flavours"* is implemented **twice over**
(this split + `isotropic_energy`). The open question is not how to PRODUCE
flavours; it is why the shipped flavour machinery has no consumers — a
consumer-side question, i.e. the deferred campaign arriving from a second
direction.

**Doc bug found while measuring it (fix on sight):**
`material_xs_field.py:760,780` (+ header `:12-13`, comment `:750`) say the two
accessors are *"Consumed by `ScatteringOperator.foldable_part`/`residual_part`"*.
Since the step-3 rebind both builders read `self.scattering.per_material` (the
kernel field) directly — `scattering.py:977-984,1001-1010`. Present-tense-false;
the only live consumer is `dsa.py:263,266`. F-1 scope.

## Refuted frames (durable UNEXPLORED for this problem class)

- **Registry / plugin table for channels** — the channel set is CLOSED and small,
  and members differ by kernel TYPE, not by a config string. A registry buys
  open-set extensibility for a set that grows by a dataclass. The extensibility
  actually wanted is the field-typed collapse (a new `MaterialField` subclass,
  zero operator edits).
- **Fiber bundle / sheaf over the mesh for `MaterialField`** — the kernel is
  locally constant on each material's disjoint cell set ⟹ trivial bundle, no
  transition functions, no non-trivial section space. `dict` + index partition IS
  the whole structure. (First non-angular refutation of the bundle frame in this
  library.)
- **Orbifold / symmetry quotient over the channel set** — no group acts;
  `{S, N₂ₙ, F}` are physically distinct reactions with no orbit structure.
- **Homology / chain complex over bulk⊕trace** — no `∂²=0`; it is a biproduct and
  the iso operators are `BlockRole.BULK` with an implicit-zero trace (an
  inclusion, not a differential). The real triple there is the THREE hand-spelled
  implicit-zero boundaries (`isotropic_scattering.py:139`, `n2n.py:199`,
  `fission.py:643` — #306 item 2), wanting §3's Tier-3 extension-by-zero verb.
- **Tensor networks / MPO** — bond dimension 1 (F is rank-1); the
  energy⊗angular⊗spatial factorization is a fixed 3-factor product with no
  truncation knob. Degenerate. (L-001.)
- **Saddle-point / mixed FEM** — fires on the DIFFUSION member (DSA's low-order
  edge system), never on an energy binding. (L-007 corollary.)
- **Petrov–Galerkin (homogenization/condensation)** — different frame VERB
  (`project` = `G⁻¹M`, solution-weighted) vs the iso source's `conjugate`
  (Galerkin, SH = scattering's eigenbasis). Do not unify.
  ([[iso-source-frame-conjugation-unification]].)
- **Category theory as the deliverable** — the concrete win is captured by CS4c
  §3's three-tier discipline plus the named lift; no abstract lever adds a test.
  The ONE thing the functor language buys is the sibling-GENERATION direction,
  which the lift makes concrete. (L-001.)
- **Feynman–Kac / path integral** — the energy binding is a per-cell algebraic
  contraction; the backbone's path reading attaches to `(Ω·∇+Σ_t)⁻¹`, and the
  energy factor is correctly outside it.

## The steelman that SURVIVES (do not throw it away)

`[M]` the tier-2 **extraction chain is genuinely duplicated**: the
`FullFieldSpace`-or-raise guard at `scattering.py:788-801` and `n2n.py:145-156`
(two different messages), and `FunctionSpace.of_axes(*interior.axes[1:])` at
`scattering.py:579` and `n2n.py:157` — twice each, under two different tier names
(`from_material_xs` vs `from_solver_data`). Currying `bind(·, space)` is the only
partial application with a REUSED operand (`[M]` at the within-group site `S`,
`N₂ₙ`, `C`, `B` all bind the same `FullFieldSpace` instance). ⟹ the real, small
deliverable is **two verbs on `FullFieldSpace`** (`scalar_interior_space`, a
`require_composite(consumer=)` guard) — not a factory.
