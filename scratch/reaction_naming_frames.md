# Reaction-term naming — cross-domain frame attack

**Task**: faithful mathematical names for the three-way split of each reaction
"operator" (collision / scattering / fission, + the (n,2n) channel): (1) the
method-agnostic physics object, (2) the seam, (3) the realized monomorphic
operator.

**Grounding**: every code fact below read DIRECTLY from the worktree at
`/Users/rodrigo/git/nuclear/ORPHEUS` (branch-verified, NOT Nexus — L-005).
Author: cross-domain-attacker, 2026-07-30.

---

## 0. Verified code facts (the ground truth)

| Fact | Citation |
|---|---|
| `FissionOperator.apply` is polymorphic over THREE carriers: `FullField→FullField`, `ScalarFlux→ScalarSourceSink`, `ndarray→ndarray` | `orpheus/transport/operators/fission.py:668-674` (overloads), `:519`, `:617`, `:644` (dispatch arms) |
| `MultiplicationOperator.apply` polymorphic over TWO: `FullField`, `ndarray` | `orpheus/transport/operators/multiplication_operator.py:376-470` |
| The collision datum is ALREADY named a **coefficient** = "the *symbol* of a zeroth-order multiplication operator" | `orpheus/transport/fields/cross_section_field.py:14-15`; `multiplication_operator.py:4-7` |
| The multiplier-algebra embedding `M: L∞ → B(L²)` is already the in-code math of record for collision | `multiplication_operator.py:21-43` (law suite, Reed & Simon cite `:95-97`) |
| The realized layer already splits 1+3: local `MultiplicationOperator` is NOT an `IntegralKernelOperator`; scattering/fission ARE ("the refinement is strict") | `orpheus/transport/operators/integral_kernel_operator.py:29-44` (locality criterion), `:57-62`, `:164` |
| `IntegralKernelOperator.kernel` returns the REALIZED integral structure (a `LinearOperator`: fission `χ⊗νΣf` TP, scattering `R∘Λ∘M` product) | `integral_kernel_operator.py:194-211` |
| The per-material physics data: `Mixture` holds `SigT (NG,)`, `SigS: list[csr (NG,NG)]` per Legendre order, `Sig2 (NG,NG)`, `chi (NG,)` simplex-enforced, `SigP (NG,)` | `orpheus/data/macro_xs/mixture.py:57-65`, `:67-74` |
| The spatial assignment: `MaterialXSField = {material_id: Mixture} + mesh.mat_map` ("the per-material Mixture dict + the spatial distribution carried by the mesh's mat_map") | `orpheus/transport/mesh/material_xs_field.py:25-28`, `:166-167` |
| The 0-D degenerate phase space is already admitted (meshless single-region `MaterialMesh`, campaign #276) | `material_xs_field.py:107-110`, `:200-206` |
| Diffusion's realization consumes MULTIPLE channels at once: `Σ_tr = Σ_t − rowsum(Σ_s1)`, `D = 1/(3Σ_tr)` ("the Fick's-law coefficient") | `mixture.py:139-166`, `:169-191` |
| CP's realization likewise: Case–Zweifel `c = (Σ_s + νΣ_f)/Σ_t` | `mixture.py:254-281` |
| The physics object has its OWN morphisms, transport-free: `condense` is "data-native (no transport dependency)" | `mixture.py:313-330` |
| Boundary precedent: `BoundaryTraceLaw` = "pure descriptor" of the affine relation `γ₋ψ = RGγ₊ψ + q`, NO `apply`, "the realizer is the **sole** bridge from descriptor to operator" | `orpheus/geometry/boundary/_base.py:3-29`, `:77-79`, `:98-108` |
| The seam word is already house grammar: campaign P2b mints `SNReactionRealizer` "on the `SNBoundaryRealizer` precedent"; explicit realizers, no registry | `.claude/plans/operator_strategy_realization_campaign.md:129-134`, `:317`, `:482-494` |

---

## 1. The mathematics, sharpened — two theorems that do the naming work

The user's unifying observation ("an energy-space operator, varying over the
spatial domain, differing only in the rank/structure of the energy-space
part") is not an analogy. It is two standard theorems, one per axis.

### Theorem A (spatial axis) — reactions are DECOMPOSABLE operators; the physics object is a FIELD OF OPERATORS

Direct-integral theory (von Neumann 1949 "On rings of operators. Reduction
theory", Ann. Math. 50; Dixmier, *Von Neumann Algebras*, Part II — "champs
mesurables d'opérateurs"; Takesaki I, Ch. IV): on
`H = L²(D; H_fib) = ∫⊕ H_fib dr`, an operator `T` is **decomposable** —
`T = ∫⊕ K(r) dr` for a **measurable field of operators** `r ↦ K(r)` —
**iff `T` commutes with every multiplication by an `L∞(D)` scalar** (locality
in `r` is the defining invariant, not a convenience).

- All four reaction terms are local in `r` (a neutron reacts where it is):
  collision, scattering, fission, (n,2n) are ALL decomposable. Streaming
  `Ω·∇` is NOT (it differentiates across cells). So "reaction vs streaming"
  = "decomposable vs not" — a theorem-grade split.
- Piecewise-constant materials mean the field is a **simple** (step) field
  that **factors through the material map**: `K(r) = K_{mat(r)}`, i.e.
  `D → {material ids} → B(H_fib)`. The code's `{material_id: Mixture} +
  mat_map` layout (`material_xs_field.py:166-167`) IS this factorization —
  the storage is the mathematical object, not an encoding of it.
- **0-D check** (the brief's degeneracy question): the direct integral over a
  one-point base IS the fiber. The 0-D infinite-medium object (one `(ng,ng)`
  matrix per reaction) is the same genus with `D = {pt}` — no special case.

Note the TWO orthogonal locality axes, which the in-code §5.6 criterion
(`integral_kernel_operator.py:29-44`) should not be confused with:
locality **in r** (all four reactions: yes; streaming: no) fixes the GENUS
(field of fiber operators); locality **within the fiber** (collision: yes —
a multiplier; transfer terms: no — kernels) fixes the SPECIES. The landed
Operator/Kernel refinement is the species split; Theorem A is the genus.

### Theorem B (angular axis) — material isotropy ⟹ the fiber is SO(3)-invariant ⟹ the isotypic (Legendre) matrix sequence is the NORMAL FORM

Fiber `H_fib = L²(E) ⊗ L²(S²)`, SO(3) acting on the angular factor. Every
reaction fiber commutes with the rotation action (the medium is isotropic:
the data depend on `Ω'·Ω` only, or not on angle at all). Schur + Funk–Hecke:
the fiber operator is block-scalar per irrep ℓ with an ENERGY-MATRIX block:

| reaction | isotypic blocks (per ℓ) | energy block structure |
|---|---|---|
| collision `Σ_t` | `Σ_t`-diagonal on EVERY ℓ (identity in angle) | diagonal |
| scattering `Σ_sℓ` | `Σ_sℓ` on ℓ ≤ L, zero above | full `(ng,ng)` per ℓ |
| fission `χ⊗νΣ_f` | ℓ=0 ONLY (isotropic in, isotropic out) | rank-one |
| (n,2n) `2Σ₂` | ℓ=0 only | full `(ng,ng)` |

So the whole four-reaction physics content is **one mathematical shape**:
*a per-material sequence of energy matrices indexed by the SO(3) irrep ℓ* —
and `Mixture`'s storage (`SigT` vector, `SigS` list-per-ℓ, `Sig2`, factored
`(chi, SigP)`) is EXACTLY this normal form (`mixture.py:57-65`). The species
differ only in WHERE they sit on the (ℓ-support × energy-rank) lattice —
the precise version of the user's "differ only in rank/structure".

Validity boundary: an anisotropic (ordered, polarized) medium breaks
zonality and the normal form with it. `Mixture` cannot represent such a
medium — which is CORRECT for the scalar-XS data model, and the frame names
exactly what would have to change if it ever entered scope.

---

## 2. Frame table (deliverable §1)

| Frame | What it names | Where faithful | Where it BREAKS |
|---|---|---|---|
| **Direct integral / decomposable operator / measurable field of operators** (von Neumann 1949; Dixmier II; Takesaki I Ch. IV) | The GENUS of object (1): the datum `r ↦ K(r)`; decomposability = spatial locality (commutant characterization); `{mat: data} + mat_map` = the factored simple field; 0-D = fiber over a point | All four reactions; also names why streaming is excluded | As a class token "Field" collides with the flux-carrier `Field` family — a cost `MaterialXSField` already bears (and under this frame its name is exactly right: it IS a field of fiber-operator data) |
| **Schwartz kernel / Fredholm integral operator** (Schwartz kernel theorem; Kress, *Linear Integral Equations*: degenerate kernels) | The SPECIES of the three transfer terms: "kernel" = the apply-free datum of an integral operator; fission = **degenerate (separable, rank-one) kernel**, stored factored | Scattering, fission, (n,2n) fibers — genuine function kernels in `(E,Ω)`; at multigroup level (counting measure on groups) even the diagonal is an honest function (Kronecker, not Dirac) | Collision in the CONTINUUM: kernel is `Σ_t δ` — a distribution, not a function. So "kernel" is NOT the family word for all four; collision's honest word is **coefficient** (below). NE false friend: in Boltzmann kinetic theory "collision kernel" names the SCATTERING factor — `CollisionKernel` for `Σ_t` would actively misread |
| **Multiplier algebra `M: L∞ → B(L²)` / PDE variable coefficient** (Reed & Simon I §VII) | The collision SPECIES: datum = **coefficient** (the symbol of a zeroth-order operator), realized = **multiplication operator**, `spec = ess-range` | Collision (and 1/v time-mass). ALREADY LANDED: `cross_section_field.py:14-15`, `multiplication_operator.py:21-43` | Does not extend to the transfer terms (they are not multipliers) — which is the POINT: the 1+3 species split is real |
| **Funk–Hecke / SO(3) representation theory** (Schur; Funk–Hecke; L-009) | The `Σ_sℓ` stack = the zonal kernel's **isotypic (Legendre) spectrum**; extended by Theorem B to the genus normal form for all four | The angular structure of every reaction fiber; NE-standard vocabulary ("Legendre moments of the scattering kernel") | Anisotropic media (no zonality); says nothing about the spatial axis |
| **Constitutive law** (continuum mechanics) | The physics-vs-balance MOTIVATION; and the word "law" for genuine CLOSURE RELATIONS — which the code already reserves correctly: `BoundaryTraceLaw` (a relation closing the problem), "the Fick's-law coefficient" (`mixture.py:172`, a closure eliminating the current) | Boundary side; diffusion's Fick closure | Reaction terms are NOT closures: nothing is underdetermined — they are the linear entries of the balance generator itself. Genus mismatch: a law RELATES (is enforced); a kernel/coefficient IS APPLIED. NE false friend: "scattering law" = the thermal scattering law S(α,β), ENDF File 7 — `ScatteringLaw` misreads to any nuclear-data reader |
| **Symbol / pseudodifferential calculus** | Collision only: `σ(x,ξ) = Σ_t(x)`, a ξ-independent zeroth-order symbol — the in-code usage is already exactly this | Collision; streaming (`iΩ·ξ`, first order) when that discussion arises | Scattering/fission are NOT ΨDOs — they act along the `(E,Ω)` fibers, not cotangent directions; a ξ-INDEPENDENT operator-valued symbol collapses to Frame 1's field-of-operators, so "symbol" as the family word imports microlocal machinery with zero content |
| **Section of a vector bundle / sheaf** | Nothing beyond Frame 1 | — | The bundle is TRIVIAL (one global trivialization = the materials dict); no transition functions, no connection, no curvature, no gluing is ever used. The load-bearing structure (measurability + finite image) is exactly Dixmier's simple field. REJECT: no structure to expose. (MoC's planned sheaf usage is for SOLUTION structure — different object.) |
| **Operator-valued measure / density** | The volume-integration face `μ(cell) = ∫_cell K dV` — why homogenization volume-weights (L-010) | A remark inside Frame 1 | The held object is the density, not the measure; no independent lever beyond Frame 1 + the landed collapse-morphism rulings |
| **Category theory: the seam as a functor** | Correctly names the OTHER seam: operator→matrix **assembly** IS functorial (`assembled ≡ probed`, landed 2b) | The assembly seam only | The physics→operator seam is NOT a functor: it is natural along NEITHER kernel composition (harmless — physics objects are never composed pre-realization) NOR phase-space coarsening — `reduce-discrete ≠ discretize-reduce` is a MEASURED codebase fact (DSA consistency had to Schur the ASSEMBLED operator; the C2 same-mesh contrast ladder). "Functor" over-promises; do not use it for the seam |
| **Kalman realization theory** (control: transfer function → state-space `(A,B,C,D)`; Kalman 1963, minimal realizations) | The SEAM word: one abstract input-output object, MANY concrete realizations; choosing a method = choosing a realization. Mathematically validates the house word (`realize`/`Realizer`) | The seam, exactly | It is a precedent, not imported machinery (no minimality theory needed) — acceptable: the word is the payoff, and it makes no false promise (unlike representation/quantization/functor) |

---

## 3. Recommendation (deliverable §2)

### 3.0 The answer is **1 + 3**, not a uniform family — stated plainly

There is **no faithful single word for all four reactions at the datum level.**
The split is a theorem, not a taste: locality *within the fiber* (Theorem A's
second axis) partitions the four into a **multiplier** and three **kernels**,
and both candidate uniform words fail in a way that would mislead a reader:

- `Kernel` for collision — the continuum object is `Σ_t δ(E−E′)δ(Ω′·Ω−1)`, a
  **distribution**, not a function; `spec` is the essential range of a
  function, not a Fredholm spectrum. Plus the Boltzmann false friend
  (kinetic theory's "collision kernel" names the SCATTERING factor).
- `Coefficient` for the transfer terms — they are not multipliers; the whole
  content of `Σ_s`, `χ⊗νΣ_f`, `2Σ₂` is the group **coupling**, which the word
  "coefficient" hides.

So: **species words on the leaves, genus word on the container.** The two
species names are already the code's own vocabulary one layer down
(`MultiplicationOperator` vs `IntegralKernelOperator`, "the refinement is
strict", `integral_kernel_operator.py:29-44`) — the recommendation makes
layer (1) *isomorphic* to that landed split instead of inventing a parallel one.

### 3.1 The names

| Layer | Object | **Recommended name** | Mathematical term it comes from | NE term (docs/accessors) |
|---|---|---|---|---|
| (1) **genus** | the apply-free per-channel physics datum (ABC) | **`ReactionTerm`** | *decomposable operator datum* / measurable field of operators (von Neumann 1949; Dixmier II). The name is deliberately weak — see §6 | "the scattering / fission term" |
| (1) species **A** (multiplier) | collision `Σ_t`; also `1/v` and any pure removal | **`CollisionCoefficient`** | variable coefficient = **symbol of a zeroth-order multiplication operator**; multiplier algebra `M: L∞→B(L²)` (Reed & Simon I §VII) | total cross section `Σ_t` |
| (1) species **B** (kernel) | scattering `Σ_{s,ℓ}` | **`ScatteringKernel`** | Fredholm/Schwartz kernel; ℓ-graded **isotypic spectrum** of a zonal kernel (Funk–Hecke, Theorem B) | scattering kernel / transfer matrix / Legendre moments |
| (1) species **B** | fission `χ⊗νΣ_f` | **`FissionKernel`** | **degenerate (separable, rank-one) kernel** — Kress, *Linear Integral Equations*, degenerate-kernel chapter | χ, νΣ_f |
| (1) species **B** | (n,2n) `ν₂Σ₂` | **`N2NKernel`** | Fredholm kernel with **emission multiplicity** (full-rank multiplying kernel) | (n,2n) matrix, `Sig2` |
| (1) **fiber** | one material's four channels together | **`Mixture`** (KEEP) | the **fiber** of the field. The all-channel tuple has NO math name — §5 | cross-section set / library entry |
| (1) **field** | `{material: fiber} + mat_map` | **`MaterialXSField`** (KEEP — verdict §3.4) | **simple measurable field of operators factoring through the material map** (Dixmier II) | XS field over the mesh |
| (2) **seam** | (1) + a discrete phase space ⟶ operator | **`<Method>ReactionRealizer`** (KEEP the campaign mint, `SNReactionRealizer`) | **realization** — Kalman 1963: one transfer function, many state-space realizations; choosing a method = choosing a realization | "builds this method's operator" |
| (3) **realized** | monomorphic operator, one domain one codomain | **`<Channel>Operator` parameterized `LinearOperator[Din, Cout]`**; species bases already correct (`MultiplicationOperator` / `IntegralKernelOperator`) | multiplier / integral operator; each carrier-specific arm is a **compression** `P A P*` (Halmos, *Hilbert Space Problem Book*, compressions) | the P0 / moment / full-field form |

Domain-first `<Domain><Role>` holds throughout: domain token = the reaction
channel (`Collision`/`Scattering`/`Fission`/`N2N`) or the index set
(`Material`), role token always present (`Coefficient`/`Kernel`/`Field`/
`Realizer`/`Operator`).

### 3.2 Greppability

```bash
# the genus as a set — the ABC plus every leaf's declaration line
grep -rn 'ReactionTerm' orpheus/
# all four leaves as one set; the regex SHAPE encodes the 1+3 split
grep -rEn '\b(Collision|Scattering|Fission|N2N)(Coefficient|Kernel)\b' orpheus/
# the species partition on its own
grep -rEn '\b[A-Z][A-Za-z]*(Coefficient|Kernel)\b' orpheus/
# the seam family (already greppable today via the boundary precedent)
grep -rEn '\b[A-Z][A-Za-z]*Realizer\b' orpheus/
```

The genus is greppable by **one token** (`ReactionTerm`) precisely *because*
the leaf names are not uniform — the ABC is what recovers the set. That is the
price of the 1+3 honesty, and it is paid once.

### 3.3 Three structural properties the naming buys (each stated as a fail-able gate)

1. **Genus gate — decomposability (Theorem A).** For every `ReactionTerm`'s
   realized operator, `A(m ⊙ x) == m ⊙ A(x)` for a **cell-varying** random
   scalar mask `m` (an `L∞(D)` multiplication). All four PASS; **the streaming
   operator `L` FAILS**. Discriminator: `m` MUST vary cell-to-cell — a constant
   `m` cannot fail, and a shape-only "is it a reaction?" check would admit a
   streaming-like term that this gate rejects.
2. **Species gate — locality within the fiber.** `CollisionCoefficient`'s
   realized operator additionally commutes with an arbitrary **group-wise**
   mask; `ScatteringKernel`/`FissionKernel`/`N2NKernel` MUST FAIL that
   (they mix groups). Discriminator: the `Σ_s` used MUST be genuinely
   off-diagonal — a diagonal-only `Σ_s` (the foldable-into-`Σ_r` case,
   `material_xs_field.py` `is_p0_diagonal_with_zero_n2n`) passes the multiplier
   gate, and a design that picks the species from the *data's accidental
   diagonality* rather than from the channel is exactly the bug this gate
   catches.
3. **Monomorphism gate — layer (3), a NEGATIVE test.** The fission operator
   realized on `ScalarFlux` MUST **raise** when handed a `FullField`. Today's
   3-arm `apply` (`fission.py:668-674`) accepts all three, so **this test is RED
   on the current tree** — that is what makes it a discriminator and not
   theatrics. Companion positive claim: the three arms are three **compressions
   to different subspaces**, hence three operators with three different
   spectra; assert their `n_dof`/shape differ and that `spec(P A P*) ≠ spec(A)`
   for a case with `L ≥ 1` scattering (equal spectra would mean the compression
   was trivial).

### 3.4 Direct verdict on the two live class names

**`CrossSectionField` — CORRECT AS-IS. No change.**
It is a *discretized carrier* (`(ng, *spatial)` broadcast array) wearing
`CoefficientRole`, and its docstring already states the exact frame
("a coefficient, not a state: the *symbol* of a zeroth-order multiplication
operator", `cross_section_field.py:14-15`). Under Frame 3 that is the faithful
name at the layer it occupies (layer 3's input, not layer 1). It is **not**
made redundant by `CollisionCoefficient`: the two are different objects at
different layers — per-material apply-free datum vs phase-space-bound array —
related by the non-identity broadcast along `mat_map` (`assemble_cell_xs`).
That ≥2-realizations + applied-non-identity-morphism pair is precisely the
type-minting criterion, so both types are earned.

**`MaterialXSField` — CORRECT AS-IS. No rename. The bug is the CONTENT, not the name.**
Under Theorem A the name is *literally* the mathematics: `Material` names the
index set the field factors through, `XS` the datum, `Field` the base-indexed
genus — a simple measurable field of fiber operators. Two consequences:

- The `Field` token collides with the flux-carrier `Field` family
  (`numerics/field.py`, `ScalarFlux`, `AngularFlux`, `CrossSectionField`). The
  collision is real but **already paid and disambiguated by the domain token**:
  carrier fields are named by a physical quantity (`...Flux`,
  `CrossSection...`), this one by its index set (`Material...`). Resolving it by
  renaming would move the *singleton* and leave the large family untouched —
  net loss of signal. Keep.
- What DOES violate the layer-(1) contract the name implies is that
  `MaterialXSField` carries **nine `apply_*`/`add_*` verbs**
  (`material_xs_field.py:815-1111`), consumed by `scattering.py` and
  `isotropic_scattering.py`. An apply-free datum that applies is the
  **un-named fiber-operator action living inside the datum** (kernel Smell #16
  shape 1: the realized operator reachable by two paths — through
  `ScatteringOperator` and directly through the field). These are a *shared
  primitive*, not a twin implementation, so the fix is **relocation into the
  realized operator (or into the fiber, `Mixture`), not deletion** — and it is
  the concrete work item the naming split creates.

**One free extension the naming admits.** Species A is not collision-only: the
time-derivative mass term `1/v` and any pure-removal coefficient join it with
no new machinery (`InverseSpeedCoefficient`), because "multiplier" is the
invariant, not "cross section". Species B's two *multiplying* members
(`FissionKernel`, `N2NKernel`) differ only in energy rank — rank-one vs full —
which is the precise version of the user's "differ only in rank/structure".
Recommend the (n,2n) multiplicity be **named data on the kernel**
(`multiplicity = 2`), not pre-multiplied into `Sig2` as `2 · sig2` is today: a
welded factor 2 is an un-named operation, and naming it makes fission and
(n,2n) the same shape with different rank.

---

## 4. Verdict on the boundary precedent (deliverable §3)

**VERDICT: `Law` is right for boundaries and WRONG for reactions. Do not reuse
the word. Reuse the *realizer* half of the precedent, not the *descriptor* half.**

Three independent structural reasons, in decreasing generality:

1. **Genus mismatch — a law RELATES, a term is APPLIED.** A boundary condition
   is a **closure relation**: it relates two unknown traces
   (`γ₋ψ = RG γ₊ψ + q`) and is *enforced* to make an otherwise-underdetermined
   problem well-posed. That is the same genus as a constitutive law in
   continuum mechanics, and the same genus as Fick's law (which the code
   already labels correctly — "the Fick's-law coefficient", `mixture.py:172`).
   A reaction term is **an entry of the balance generator**: nothing is
   underdetermined without it.
   *Decidable check on the word:* delete the object and ask what breaks.
   Delete the BC law ⟹ the problem becomes **ill-posed** (a one-parameter
   family of solutions). Delete the scattering term ⟹ a **different, still
   well-posed** problem (a pure absorber). Only the first is a closure. Any
   candidate `...Law` name must survive this check.
2. **Algebraic type mismatch — affine vs linear.** The trace law is *affine*:
   it carries the inhomogeneity `q` (the prescribed inflow). All four reaction
   terms are *purely linear* — there is no reaction analogue of `q`. Sharing
   the word `Law` would flatten a real type distinction that the code spends
   effort maintaining elsewhere (the affine/torsor discipline).
3. **Nuclear-data false friend — ENDF-6 File 7.** "Thermal scattering law"
   is the standard name for `S(α,β)`, the incoherent-inelastic
   double-differential *microscopic data* representation. `ScatteringLaw` in an
   operator context would misread to every reader who has touched an ENDF
   library — and it names an object at a *different layer* (microscopic data,
   upstream of macroscopic-kernel processing) in a code that reads that
   ecosystem's files. This is not stylistic; it is a collision with a live
   term of art.

**What the boundary side should KEEP, unchanged:**

- **`BoundaryTraceLaw`** — the name is high-signal and correct: domain-first
  (`Boundary`), the middle token names *the object the relation is between*
  (`Trace`), the role token states the genus (`Law`). Keep the "pure
  descriptor, no `apply`" contract and the "the realizer is the **sole** bridge
  from descriptor to operator" invariant (`_base.py:98-108`) — that contract is
  exactly what layer (1) copies.
- **`SNBoundaryRealizer` / the realizer pattern** — this is the *transferable*
  half. "Realization" is layer-agnostic (Kalman: one abstract input-output
  object, many concrete realizations), so `SNReactionRealizer` alongside
  `SNBoundaryRealizer` is a genuine family, not a coincidence of spelling. Keep
  explicit realizers with no registry.

**Write the rule down** (so the next session cannot re-decide it):
`Law` is reserved for a relation that closes an otherwise-underdetermined
system — BC traces, Fick's law, transverse-leakage and P1/SP3 closures. It is
never a term of the generator. The policy is greppable: `grep -rEn
'\b[A-Z][A-Za-z]*Law\b' orpheus/` should return *only* closures; a reaction
channel appearing in that output is the violation.

---

## 5. Where no faithful name exists (deliverable §4)

Two of the three objects have faithful, literature-backed names. **One does
not.** Stated plainly, plus one adjacent gap that is real but should stay
unnamed.

### 5.1 Layer (2) the seam — NAMED. `realization` / `Realizer` (Kalman 1963). No gap.

### 5.2 Layer (3) the realized operator — NAMED. `multiplier` / `integral operator`, and the relation between two carrier-specific arms of one datum is a **compression** `P A P*` (Halmos). No gap. Recorded because the naive readings are both wrong: the arms are not "overloads" (no math content) and not "the same operator on different carriers" (false — a compression has a *different spectrum*).

### 5.3 Layer (1) the genus — **NO FAITHFUL NAME EXISTS.** `ReactionTerm` is INVENTED.

The exact mathematics is "the datum of a decomposable operator" — Dixmier's
*champ mesurable d'opérateurs*, a measurable field of operators. Every faithful
rendering fails as a code token:

| Candidate | Why it fails |
|---|---|
| `DecomposableOperator` | `Operator` promises `apply`; the object has none. "Decomposable" also collides with tensor-rank vocabulary. |
| `MeasurableOperatorField` | Mathematically exact and useless: the base is a finite mesh, so measurability is trivially satisfied and carries **zero content**. Triples the `Field` collision. |
| `ReactionOperatorField` | Exact under Theorem A **when the base is attached**, but (a) it is wrong for the base-free per-material fiber (0-D / meshless — the one place the fiber/field distinction bites), and (b) it contains the substring `ReactionOperator`, which pollutes the grep guarding the standing #261 anti-mint ruling ("do not mint `ReactionOperator`"). |
| `OperatorDatum` / `ReactionSpec` / `ReactionDescriptor` | Zero topic signal; `Descriptor` additionally collides with the Python descriptor protocol. |

**Invented name, flagged as invented: `ReactionTerm`.** It is not a standard
mathematical object word. It is chosen as the **weakest true statement**
available — "one summand of the balance generator `A = L + C − S − F − B`" —
because every stronger word over-promises. It is NE-readable (standard
transport prose says "the scattering source term"), Domain-first with the role
token present, and greppable as one token.

Its cost is explicit: **the genus's defining invariant is not in the name.**
Decomposability (locality in `r`, the property that admits all four reactions
and excludes streaming) must be carried by the ABC docstring and enforced by
§3.3 gate 1. A name that does not carry its own invariant needs a test that
does.

### 5.4 The all-channel aggregate — no MATH name, and none should be invented.

A tuple of one multiplier plus three kernels describing one medium is not a
named mathematical object; there is no theorem about the tuple, so mathematics
offers only the positional word *fiber*. Nuclear engineering does have names —
"cross-section set", "material data", "library entry". **Keep `Mixture`**: it
is the NE-faithful word for a homogenized medium, and it is exactly the fiber of
Theorem A's field. Inventing `MaterialData`/`MixtureData` adds ceremony and no
signal.

### 5.5 The seam as a *structure* (not as a map) — real gap, should stay unnamed.

There is no faithful word for "the pair (physics datum, discrete phase space)
⟶ operator" *as an algebraic structure*, because the honest candidate is
"functor" and that is measurably false here (§6). The map has a name
(`realization`); demanding a name for the structure is what would import the
over-promise. Leave it unnamed deliberately — this is a decision, not an
omission.

---

## 6. What each name buys or costs (deliverable §5)

### 6.1 The recommended names

| Name | Buys | Costs | NE-readability trade |
|---|---|---|---|
| `ReactionTerm` (genus ABC) | The one greppable token that recovers the 1+3 set; promises nothing false; NE prose already says "term" | **Zero math content** — the decomposability invariant lives only in the docstring + §3.3 gate 1. Invented (§5.3) | None — an engineer and a mathematician both read "term" correctly |
| `CollisionCoefficient` | The multiplier-algebra frame with its consequences: `spec = ess-range`, `M[0] = ZeroOperator`, the vector-space (not affine) algebra; already the code's own words | Deliberately breaks the family word, so a reader scanning for "the four kernels" finds three (the ABC grep recovers them) | **REAL trade.** An engineer says *total cross section*, never "collision coefficient". Use `CollisionCoefficient` as the **type**; keep `total_cross_section` / `Σ_t` for **accessors, docs, equations** |
| `ScatteringKernel` | Fredholm/Schwartz frame; the ℓ-graded storage is named as the *isotypic spectrum* rather than "a list of matrices" | Collides with `IntegralKernelOperator.kernel` (§6.3) | **None** — "scattering kernel" is standard NE. Highest-agreement name in the set |
| `FissionKernel` | "**Degenerate kernel**" names the factored `(chi, SigP)` storage exactly — the rank-one structure becomes a stated property, not a storage trick | Same `kernel` collision | Mild: an engineer thinks "χ and νΣf", not "kernel". Docs keep `χ ⊗ νΣ_f` |
| `N2NKernel` | Puts (n,2n) in the same species as fission, differing only in energy rank (full vs rank-one) — makes the multiplicity nameable instead of welded as `2 · sig2` | `N2N` is a house abbreviation | None — `n2n` is already the house/NE token |
| `Mixture` (keep) | NE-faithful; is literally Theorem A's fiber | No math name exists for the aggregate (§5.4) | None |
| `MaterialXSField` (keep) | Literally Theorem A's simple field factoring through the material map | The `Field` token overload (§6.3); and today the class violates the apply-free contract its name implies (§3.4) | None — `XS field` is NE-native |
| `<Method>ReactionRealizer` (keep) | Kalman precedent; matches the landed `SNBoundaryRealizer` so the two seams read as one family | None material | **None** — "realizer" reads correctly as plain English *and* as control theory. Rare free name |
| `LinearOperator[Din, Cout]` monomorphic (layer 3) | Makes adjoint / inverse / spectrum well-defined by construction; the illegal carrier becomes unspellable (§3.3 gate 3) | Multiplies the number of realized objects per channel (fission: 3) | Docs should say "the P0 form" / "the moment form"; **compression** belongs in theory pages, not in a class name (`CompressedOperator` would be over-typing) |

### 6.2 Names that promise an operation the objects do NOT support — explicit cost lines

| Rejected name / frame | The operation it promises | Why the objects cannot support it |
|---|---|---|
| **`...Functor`** for the seam | Naturality: that realizing-then-reducing equals reducing-then-realizing | **Measured false in this codebase.** `reduce-discrete ≠ discretize-reduce`: DSA consistency had to Schur-eliminate the **assembled** operator, and the C2 same-mesh contrast ladder measured the coarse-resolve alternative failing. A functor name would license a refactor that the numbers already refute |
| **`ReactionSymbol`** / pseudodifferential framing | A symbol *calculus*: composition-to-leading-order, wavefront/propagation-of-singularities machinery | Scattering/fission act **along the fibers**, not in cotangent directions; a ξ-independent operator-valued symbol collapses back to Theorem A's field of operators. Imports microlocal machinery with **zero content**. (For collision alone the word is already used correctly and locally — keep that usage, do not promote it) |
| **`ReactionSheaf` / `...Section`** | Gluing, transition functions, a connection, cohomology | The bundle is **trivial** — one global trivialization *is* the materials dict. No transition function, connection, or gluing is ever used; `H¹` of a trivial bundle over a contractible base is zero content. (MoC's planned sheaf usage is for SOLUTION structure — a different object) |
| **`DecomposableOperator`** | `apply` | The layer-(1) object has no `apply` — that is its entire point (the realizer is the sole bridge) |
| **`ReactionLaw` / `ScatteringLaw`** | That it is a closure to be *enforced*, plus an affine inhomogeneity | Reaction terms are linear generator entries, not closures (§4); `ScatteringLaw` additionally collides with ENDF File 7 `S(α,β)` |
| **`CollisionKernel`** | That it is the transfer factor | In Boltzmann kinetic theory "collision kernel" **is** the scattering factor — the name would designate the wrong physics to a reader who knows the term |
| **`ReactionMeasure`** (operator-valued measure) | `μ(A ∪ B) = μ(A) + μ(B)` additivity as the primary interface | The held object is the **density**, not the measure; the volume-integration face is one derived remark and buys no lever the landed collapse-morphism rulings do not already have |

### 6.3 The two token collisions this recommendation creates, and their resolution

1. **`Kernel` at layer 1 vs `IntegralKernelOperator.kernel` at layer 3.** The
   accessor at `integral_kernel_operator.py:195` returns a
   **`LinearOperator`** — a *realized* object — while the recommended layer-1
   species word `Kernel` means the *apply-free datum*. Same word, two layers,
   opposite contracts. **Resolution: rename the accessor** (`as_integral_operator`,
   or `kernel_operator` if the shorter name is preferred) and keep `Kernel` for
   the datum, which is where the literature puts it. Cheap and one-directional.
2. **`Field`: carrier vs field-of-operators.** `CrossSectionField`/`ScalarFlux`
   are discretized array carriers; `MaterialXSField` is Theorem A's field of
   fiber operators. **Resolution: keep both**, disambiguated by the domain
   token — carriers are named by a *physical quantity*, the operator-field by
   its *index set* (`Material...`). Renaming would move the singleton and leave
   the large family untouched.

### 6.4 The cost that would bite hardest if the split is done naively

Minting four **closed** per-channel types breaks the derived-data path that
consumes several channels *at once, before any realization*:
`transport_xs = Σ_t − rowsum(Σ_{s1})` and `diffusion_coefficient = 1/(3Σ_tr)`
(`mixture.py:139-191`), and CP's Case–Zweifel `c = (Σ_s + νΣ_f)/Σ_t`
(`mixture.py:254-281`). Those are cross-channel *arithmetic on the data*.

Therefore: the layer-1 species must keep the `CoefficientRole` vector-space
algebra (`Σ + Σ′` legitimate, `Σ = 0` a genuine origin — already true of
`CrossSectionField`), and the family should be **typed views over `Mixture`**,
not a parallel object graph beside it.

*Discriminating test for this cost:* after the split, `Mixture.transport_xs`
and `Mixture.scattering_ratio` must still evaluate with **no mesh and no
realized operator in scope**. A design in which obtaining `D` requires
realizing an operator has broken the method-agnostic contract — and it fails
loudly, because diffusion consumes `D` upstream of realization.

### 6.5 The readability rule (which name to use where)

**Types get the mathematically-faithful name; accessors, attributes, docstrings,
theory pages, and equations get the nuclear-engineering name.** A type name is
read once per design decision and must carry the invariant; an accessor is read
on every line of physics code and must carry the domain word. Both columns are
tabulated in §3.1 — the recommendation is that both be kept, not that one
replace the other. The only place the two must agree is the **grep pattern**:
the type family must be findable as a set (§3.2), because that is what stops a
fifth reaction channel from being added off-pattern.

### 6.6 Rejected names, checked (the record of what not to re-attack)

`Functor` (naturality measurably false) · `Symbol`/ΨDO (fiber-acting, not
cotangent — zero content) · `Sheaf`/`Section` (trivial bundle, no gluing) ·
`Measure` (density is the held object) · `Law` (closure genus + ENDF File 7) ·
`DecomposableOperator` (promises apply) · `CollisionKernel` (Boltzmann false
friend) · `MeasurableOperatorField` (measurability trivial on a finite mesh) ·
`ReactionOperatorField` (wrong for the base-free fiber; pollutes the #261
`ReactionOperator` anti-mint grep) · `ReactionDescriptor` (Python-descriptor
collision, no signal) · `CompressedOperator` (compression is a *relation
between* realizations, not a kind of object).
