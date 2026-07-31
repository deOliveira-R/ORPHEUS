---
name: reaction-term-naming-species-split
description: The three-layer naming verdict for reaction terms (collision/scattering/fission/n2n) — genus = decomposable-operator field (von Neumann/Dixmier), species = 1 multiplier + 3 kernels (theorem, not taste), seam = Kalman realization, layer-3 arms = compressions; `Law` reserved for closures; MaterialXSField/CrossSectionField both keep their names.
metadata:
  type: project
---

# Reaction-term naming — the 1+3 species split (2026-07-30)

Deliverable: `scratch/reaction_naming_frames.md` (§0 code facts, §1 Theorems A/B,
§2 frame table, §3 recommendation, §4 Law verdict, §5 no-name-exists, §6
buys/costs). Branch `refactor/operator-strategy-layers`; all facts read from the
worktree (L-005). Context: the operator/splitting/realization campaign needs a
three-way split per reaction term — method-agnostic datum / seam / monomorphic
operator.

## The two theorems that do the naming work

- **Theorem A (spatial):** an operator is *decomposable* (`T = ∫⊕ K(r) dr`, a
  measurable field of operators — von Neumann 1949; Dixmier II *champs
  mesurables d'opérateurs*) **iff** it commutes with every `L∞(D)`
  multiplication. All four reactions qualify; **streaming does not**. So
  "reaction vs streaming" is theorem-grade, and `{material: Mixture} + mat_map`
  IS the *simple field factoring through the material map* — the storage is the
  object. 0-D = the fiber over a one-point base (no special case).
- **Theorem B (angular):** isotropy ⟹ each fiber is in the SO(3) commutant ⟹
  Schur + Funk–Hecke ⟹ block-scalar per irrep ℓ with an energy-matrix block.
  `Mixture`'s storage (`SigT` vector, `SigS` list-per-ℓ, `Sig2`, factored
  `(chi, SigP)`) IS that normal form. The four species differ only in position
  on the (ℓ-support × energy-rank) lattice — the exact version of "differ only
  in rank/structure".

## The verdict (what to name what)

| Layer | Recommended | Literature term | NE term |
|---|---|---|---|
| (1) genus ABC | `ReactionTerm` — **INVENTED**, deliberately weakest-true | (exact word "decomposable-operator datum" unusable as a token) | "the scattering term" |
| (1) species A | `CollisionCoefficient` (+ free: `InverseSpeedCoefficient`) | symbol of a 0-order multiplication op; multiplier algebra `M:L∞→B(L²)` | total cross section Σ_t |
| (1) species B | `ScatteringKernel` / `FissionKernel` / `N2NKernel` | Fredholm/Schwartz kernel; fission = **degenerate (rank-one) kernel** (Kress) | scattering kernel; χ⊗νΣf; (n,2n) matrix |
| (1) fiber | `Mixture` KEEP | the fiber (no math name for the all-channel tuple) | cross-section set |
| (1) field | `MaterialXSField` KEEP | Dixmier simple field | XS field |
| (2) seam | `<Method>ReactionRealizer` KEEP | **realization** (Kalman 1963) | "builds this method's operator" |
| (3) realized | `<Channel>Operator` as `LinearOperator[Din,Cout]`; bases already right (`MultiplicationOperator` / `IntegralKernelOperator`) | multiplier / integral op; each carrier arm = a **compression** `PAP*` (Halmos) | the P0 / moment form |

**No uniform family word exists** (the headline): `Kernel` for collision is a
*distribution* (`Σ_t δ`) plus the Boltzmann false friend; `Coefficient` for the
transfer terms hides the group coupling. Species on the leaves, genus on the
ABC. Grep set: `grep -rEn '\b(Collision|Scattering|Fission|N2N)(Coefficient|Kernel)\b'`;
genus by the single token `ReactionTerm`.

**Layer-1 ↔ layer-3 species bijection** (multiplier↔`MultiplicationOperator`,
kernel↔`IntegralKernelOperator`) makes the realizer *species-preserving* — the
honest structural property of the seam, strictly weaker than functorial.

## `Law` — RIGHT for boundaries, WRONG for reactions (three reasons)

1. **Genus:** a law is a **closure** that is *enforced*; a reaction term is a
   generator entry that is *applied*. Decidable check: **delete it and ask what
   breaks** — delete the BC law ⟹ **ill-posed**; delete scattering ⟹ a
   *different, still well-posed* problem (pure absorber). Only the first is a
   closure. (Reusable test for any `...Law` candidate name.)
2. **Algebra:** the trace law is *affine* (carries `q`); reaction terms are
   purely linear — sharing the word flattens a maintained type distinction.
3. **False friend:** ENDF-6 **File 7 "thermal scattering law" `S(α,β)`** — a
   different object at a different layer in the same data ecosystem.

Keep on the boundary side: `BoundaryTraceLaw` (+ its apply-free "realizer is
the sole bridge" contract) and `SNBoundaryRealizer`. The **realizer** half IS
transferable (Kalman is layer-agnostic); the **descriptor** word is not.
Policy grep: `grep -rEn '\b[A-Z][A-Za-z]*Law\b' orpheus/` must return only
closures (BC traces, Fick, transverse leakage, P1/SP3).

## Costs the recommendation creates (and their resolutions)

- **`Kernel` collision:** `IntegralKernelOperator.kernel` (`:195`) returns a
  **`LinearOperator`** — a realized object — while layer-1 `Kernel` is the
  apply-free datum. FIX: rename the accessor (`as_integral_operator`), keep
  `Kernel` for the datum where the literature puts it.
- **`Field` collision** (carrier vs field-of-operators): KEEP both —
  disambiguated by domain token (carriers named by a physical *quantity*, the
  operator-field by its *index set* `Material...`). Renaming moves the singleton
  and leaves the large family untouched.
- **Cross-channel arithmetic is load-bearing:** `transport_xs = Σ_t −
  rowsum(Σ_s1)`, `diffusion_coefficient`, CP's `c = (Σ_s+νΣ_f)/Σ_t` consume
  SEVERAL channels *before* realization. Four *closed* per-channel types would
  break them ⟹ species keep the `CoefficientRole` vector-space algebra and are
  **typed views over `Mixture`**, not a parallel object graph. Discriminating
  test: after the split, `Mixture.transport_xs` must evaluate with **no mesh and
  no realized operator in scope**.
- **Readability rule:** *types* get the math-faithful name; *accessors,
  docstrings, theory pages, equations* get the NE name (`total_cross_section`,
  `Σ_t`). Only the grep pattern must agree.

## The live content bug the naming exposes (not a rename)

`MaterialXSField` carries **nine `apply_*`/`add_*` verbs**
(`material_xs_field.py:815-1111`), consumed by `scattering.py` +
`isotropic_scattering.py` — an apply-free datum that applies. Shared primitive,
not a twin ⟹ **relocate into the realized operator (or the fiber), do not
delete**. Also: `2 · sig2` welds the (n,2n) multiplicity — name it
(`multiplicity = 2`), which makes fission and (n,2n) the same shape at
different rank.

## Rejected names, with structural reasons (do not re-attack)

`Functor` for the seam — naturality **measurably false** (`reduce-discrete ≠
discretize-reduce`: DSA had to Schur the ASSEMBLED operator; the C2 contrast
ladder) · `Symbol`/ΨDO — fiber-acting not cotangent; ξ-independence collapses to
Theorem A (zero content) · `Sheaf`/`Section` — trivial bundle, no
gluing/connection ever used · `Measure` — the held object is the density ·
`DecomposableOperator` — promises `apply` · `CollisionKernel` — Boltzmann names
the *scattering* factor · `MeasurableOperatorField` — measurability trivial on a
finite mesh · `ReactionOperatorField` — wrong for the base-free fiber AND
pollutes the #261 `ReactionOperator` anti-mint grep · `ReactionDescriptor` —
Python-descriptor collision · `CompressedOperator` — compression is a *relation
between* realizations, not a kind of object.

Related: [[coefficient-field-promotion-frames]] (the multiplier-algebra frame
this extends), [[iso-source-frame-conjugation-unification]],
[[fission-rank1-normal-form-dead-functional]] (the rank-1 normal form),
[[issue-261-cross-method-operator-relocation]] (the `ReactionOperator`
anti-mint), [[harmonic-frame-ownership-funk-hecke]] (Theorem B's owner rule).
