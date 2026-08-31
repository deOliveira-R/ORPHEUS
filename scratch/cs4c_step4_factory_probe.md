# CS4c step-4 — adversarial probe of the "iso-family Factory object"

**Agent:** cross-domain-attacker. **Date:** 2026-08-30. **Route:** direct
`Read`/`Grep` on the live worktree (authority). Nexus MCP tools WERE available but
its graph is stamped at commit `92dcc30f` — **pre step-3 close-out** — so several
positions came back flagged `stale`; I used `twin_paths`/`discriminations` only as
a corroborating sweep and grounded every claim below on a file read (lessons L-005).

**Task shape:** structural detection + tension. Not a verdict.

---

## 0. What the probe found first, because it changes the question

Three facts, all `[M]` by read, that the proposal's justifications do not yet
account for:

1. **The aniso/iso co-sourcing the factory would provide ALREADY SHIPS, by a
   strictly stronger mechanism.** `ScatteringOperator.isotropic_energy`
   (`orpheus/transport/operators/scattering.py:555`, a `cached_property`) returns
   an `IsotropicScattering` bound to *its own* `ScatteringMaterialField`;
   `N2NOperator.energy` (`n2n.py:108`) is the same shape. A factory co-sources two
   operators **by convention** (same inputs, two calls). A satellite property
   co-sources them **by construction** (one datum instance, one object graph).
   Cardinal Rule 2 prefers the second and the tree already has it.
2. **DSA consumes NEITHER member of the trio.** `[M]` `grep
   'IsotropicScattering|IsotropicN2N'` over `orpheus/` returns **0 hits in
   `orpheus/sn/acceleration/`**. `DSALowOrderSystem.from_sn_mesh`
   (`dsa.py:263,266`) reads `mat_xs.foldable_sigma()` → `{mid: (ng,)}` and
   `mat_xs.residual_sig_s()` → `{mid: [P0_offdiag, Σ_s1, …]}` as **raw arrays**,
   and `DSACorrection.__init__` (`dsa.py:585`) reads
   `quadrature.angular_frame(1).table[:,1,1]`. Total channel operators consumed by
   the whole DSA subsystem: **zero**.
3. **The third isotropic lifted channel has already arrived and does not use the
   shared primitive.** `_per_ordinate.py:20-22` states its own defer-until-2
   trigger verbatim: *"When a third isotropic lifted channel appears, this is the
   seed of the generic lift operator (defer-until-2 — recorded, not minted)."*
   `[M]` `FissionOperator`'s composite arm (`fission.py:615-625`) computes
   `(fission_iso / w_total) + AngularSourceSink.zeros(bulk.space)` — the same
   arithmetic as `assemble_per_ordinate_isotropic`, hand-inlined, with its own
   `w_total` derivation. **The trigger has fired.**

---

## STRUCTURAL FEATURES

**Mathematical objects**

- **F1** Three channel kernels, representation-free, per-material
  (`transport/kernels.py`): `ScatteringKernel` (ℓ-stack of general `(ng,ng)`
  transfer matrices), `N2NKernel` (one general `(ng,ng)` matrix + a class-level
  scalar multiplicity `ν₂ₙ`), `FissionKernel` (a **factor pair** `(χ, νΣf)`,
  `χ` simplex-constrained by its own `__post_init__`).
- **F2** One pairing type `MaterialField[K]` (`transport/material_field.py:74`) —
  frozen `{mid → kernel}` + the mesh's cell partition, with **one** contraction
  primitive `_accumulate_contracted(Q, x, matrix_of, *, spec, scale)` (`:137`)
  and two einsum spellings `_FORWARD = "fg,fc...->gc..."`, `_TRANSPOSE =
  "fg,gc...->fc..."` (`:69-70`).
- **F3** Two binding KINDS per channel: an **energy** binding (scalar space, no
  frame, no angular content) and an **angular** binding (composite space, frame-
  conjugated or ℓ=0-lifted).
- **F4** Spaces: `FunctionSpace` (`(ng, *cells[, 2^d])`), `FullFieldSpace`
  (interior ⊕ trace), angular interior `(N, ng, *cells)`, moment space.
- **F5** Frames: `HarmonicFrame`, interned per `(quadrature, L)` via
  `HarmonicFrame.for_space`; faces `flux_analysis` (M) / `source_reconstruction`
  (R). The blessed chain space → interior → angular axis → `generator_as(Quadrature)`
  → hub → frame.

**Symmetries**

- **F6** `Σ_s(Ω·Ω')` is SO(3)-zonal ⟹ Funk–Hecke ⟹ SH is its eigenbasis; ℓ=0 is
  the trivial irrep `V₀`. **All three iso operators live in `V₀`.**
- **F7** The **energy** operator carries no angular symmetry at all — it is a
  per-cell `(ng,ng)` contraction, cell-diagonal by the material partition.
- **F8** RANK is not uniform across the three: `Σ_{s,0}` and `ν₂ₙΣ_{2n}` are
  general (rank ≤ G); `F = |χ⟩⟨νΣf|` is **rank exactly 1**.
- **F9** N2N differs from iso-S at the array level by a **scalar unit only**:
  `add_p0_source` = `_accumulate_contracted(…, lambda k: k.p0, spec=_FORWARD)`
  (`material_field.py:223`); `add_emission` = `_accumulate_contracted(…,
  lambda k: k.matrix, spec=_FORWARD, scale=ν₂ₙ)` (`:323`). Same body, two
  arguments.

**Iterative / algebraic structure**

- **F10** Within-group loss `(L+C) − S − N₂ₙ − B`; `K_iso = S.isotropic_energy +
  N2N.energy` composed at the solver (`coupled_system.py:592`,
  `diffusion/solver.py:241-244`, `homogeneous/solver.py:197-199`).
- **F11** F is the k-**outer** operator (within-group fission is zero).
- **F12** The iso lift is `N₂ₙψ = (1/W)·E·K·∫ψ dΩ` (`n2n.py:29`) — the frame's
  ℓ=0 conjugation on the reaction-rate fast path.

**Boundary handling**

- **F13** All three are `BlockRole.BULK` with an implicit-zero trace, spelled
  **three times** (`isotropic_scattering.py:139-142`, `n2n.py:199-202`,
  `fission.py:643-645`). `#306` item 2 already names these three zero-boundary
  emitters as an open adjudication.

**Selection structure (the feature the proposal is really about)**

- **F14** `[M]` construction sites of the iso pair tree-wide = **4**, and **none
  is a branch**: `diffusion/solver.py:242-243`, `homogeneous/solver.py:197-199`,
  `scattering.py:555` (satellite), `n2n.py:158` (inside tier-2). Each site knows
  statically which binding it wants. *(Predicate: literal class-name call sites
  found by name grep over `orpheus/`; a variable-call or registry-loop member
  would be invisible to this — the 2026-08-29 §6b spelling lesson.)*
- **F15** The tier-2 **extraction chain is duplicated**, and this is the real
  measured duplication under the proposal: the `FullFieldSpace`-or-raise guard
  appears at `scattering.py:788-801` and `n2n.py:145-156` (two different
  messages), and the scalar-sub-space derivation `FunctionSpace.of_axes(*interior
  .axes[1:])` appears at `scattering.py:579` and `n2n.py:157` — **twice each**.
  The tier is also spelled with **two different names** (`from_material_xs` vs
  `from_solver_data`) for the same role.
- **F16** `W = Σw` has **four** derivations in the tree:
  `S.total_weight = flux_analysis.frame.measure.weights.sum()`
  (`scattering.py:506`); `N2N.total_weight = self.weights.sum()` off
  `axes[0].generator_as(Quadrature).weights` (`n2n.py:163-172`); F's inline
  `ang.shape[0] if ang.weights is None else sum(ang.weights)`
  (`fission.py:618-622`); and `IsotropicReconstructionOperator._total_weight =
  frame.discrete_gram[0,0]` (`numerics/frame.py:766`).

**Scale separation**

- **F17** Every member factors as (energy kernel) ⊗ (angular lift) ⊗
  (cell-diagonal material partition). Nothing couples the three axes.

---

## ELEGANCE DETECTOR HITS

- **Smell #16 shape 1** (distinct paths to ONE operator) — **fires twice.**
  (a) `IsotropicScattering` and `IsotropicN2N` are two classes over one
  arithmetic body: `apply`/`apply_transpose`/`assemble`/`is_assemblable`/
  `data_ng` differ only in the field name and a `(matrix_of, scale)` pair (F9);
  the three shared bodies are *already* hoisted to module functions
  (`_scalar_composite_source`, `_iso_is_assemblable`,
  `_assemble_iso_energy_operator`) — the arithmetic was single-sourced and the
  **class** was not. (b) `ScatteringOperator.foldable_sigma`
  (`scattering.py:1028`, `np.diag(k.p0)` off the kernel) and
  `MaterialXSField.foldable_sigma` (`material_xs_field.py:820`, `np.diag(
  sig_s_legendre(mid)[0])` off the facade) are two paths to one quantity — and
  DSA reads the **facade** path.
- **Smell #16 shape 2** (one quantity, N representations) — `W` at **four**
  spellings (F16). They agree today only because `Axis.weights is None` is
  canonicalized to all-ones somewhere else; F's fallback branch is a correctness
  claim resting on an invariant it does not name.
- **Smell #16 shape 4** (a third hand-rolled path, and here the code predicted
  it) — `fission.py:615-625` re-inlines the lift the shared
  `assemble_per_ordinate_isotropic` exists to own, and `_per_ordinate.py:20-22`
  states the exact trigger that has now fired.
- **Smell (Part C candidate, 1st sighting): "a family word proposed for a set
  whose members differ in RANK."** Rank is not decoration — it fixes the normal
  form (dense `(ng,ng)` contraction vs. dyad), the transpose MECHANISM (einsum
  index swap vs. factor swap), and the per-cell cost class (`O(G²)` vs `O(G)`).
  A "trio" word over `{S₀, N₂ₙ, F}` asserts a uniformity that F8 refutes.
- **NOT firing: `discriminations`.** There is no repeated conditional on space
  kind or channel tag anywhere in the family (F14). This matters below.

---

## 1. STEELMAN — what domain algebra would a factory realize?

The proposal is not empty. Its honest content is **F15**: the tier-2 extraction
chain is genuinely duplicated, and the duplication grows linearly with the number
of bound operators. Something wants to own it. Four candidate spellings, ranked
by how much real structure each names.

### 1a. `bind(kernel, space, order)` as a law object — the CURRYING reading

The campaign has already ruled the binary operation: *Kernel × Frame → Operator,
bound at construction*. A "factory" is a **currying** of that binary map: fix one
argument, get a function of the other. Currying pays **iff one of the partially
applied objects is itself reused**:

- `bind(K, ·)` — a kernel-parameterized family over spaces. Reused? `[M]` **no**:
  each kernel is bound to exactly one space per solve.
- `bind(·, frame)` — a space/frame-parameterized family over kernels. Reused?
  `[M]` **yes** — at the within-group site all of `S`, `N₂ₙ`, `C`, `B` bind the
  *same* `FullFieldSpace` instance, and the OperatorSum guard exists precisely to
  check that.

⟹ **If a factory is justified at all, it is `bind(·, space)` — a "pose" object
holding the space and its derived pieces, from which each channel asks for its
own binding.** That is the only currying with a reused operand, and it is exactly
the object F15 says is missing. Note what it is NOT: it is not a dispatcher that
chooses *which* operator; it is a **holder of the shared extraction**.

### 1b. `space.bind(kernel)` — REFUTED by import direction, not taste

L-020's grep-checkable discriminator: *in a data/binder split, the DATA object's
verbs return ARRAYS; only the BINDER returns OPERATORS.* `FunctionSpace` lives in
`orpheus/numerics`; `IsotropicScattering` lives in `orpheus/transport`. A
`space.bind(kernel) -> LinearOperator` forces `numerics → transport`, which is
the direction the layer contract forbids (`tests/test_layer_imports.py`
`FORBIDDEN_EDGES`) and which §14.4 already recorded as the reason the frame hub
had to be transport-side. No dispatch mechanism escapes it — the method-agnostic
module must NAME the method's types however the dispatch is spelled.

### 1c. The Frame as the factory (stage-2 generator discipline)

This one is **already realized and its boundary is exactly the interesting
finding.** `HarmonicFrame.for_space(interior, L)`, interned per `(quadrature, L)`,
induces `moment_space_on` + the two faces — a factory object inducing space and
operators together at one site, per the stage-2 discipline. Its scope is the
**frame-bearing** bindings.

The iso trio has **zero angular content** (§14.2 ruled this), so the frame-as-
factory does not reach them — *except* through the ℓ=0 sub-block. That is the
crack worth widening (Frame 6 below): the iso **lift** IS
`angular_frame(0).conjugate(K)`, and the ℓ=0 frame is the single object that owns
`W`, the ℓ=0 analysis row, and the ℓ=0 reconstruction row — i.e. the one object
that would collapse F16's four spellings of `W` to one.

### 1d. Promoting `build_within_group_system` to a posing-factory keyed on the problem

This is the maximal reading, and it is **the consumers campaign arriving early**.
§5 of the design record is explicit and dated: *"the solution path will have a
detailed campaign, and it will drastically reshape the machinery of the consumers
… nothing in this ladder redesigns iteration/posing machinery."* A factory keyed
on "the problem" is a **posing** object — and `posing` is a word already spent in
this repo on the arrangement of leaves into `(A_loss, M)` plus the eigenvalue role
(`numerics/eigenvalue.py`, `iteration.py`, `homogeneous/solver.py:225`,
`coupled_system.py:118`, `loss_representation/assembly.py:401`), about to get
sharper when Campaign 2 lands `GeneralizedEigenPencil`. Building it now under a
different name mints the *same axis of the same solve* twice.

### 1e. Which mathematical frame fits the factory as literally proposed

A **fibration** `Op → Space` with the factory as a chosen **section**. This frame
is what kills it, cleanly — see Frame 4.

### 1f. The DSA justification, checked concretely

Enumerated from the tree, an SN+DSA solve needs:

| tier | operator set | source |
|---|---|---|
| high-order | `L` (streaming), `C` (`MultiplicationOperator`), `S` (`ScatteringOperator`, full ℓ), `N₂ₙ` (`N2NOperator`), `B`, `F` | `sn/solver.py:1413-1427`, `coupled_system.py` |
| low-order | a per-group `(K+1)×(K+1)` **edge-centered** LU-factored matrix `A_g` + increment→source map `G_g` | `dsa.py:340-392` |
| coupling | `R` = the frame's ℓ≤1 analysis rows; `P` = Larsen (33) synthesis | `dsa.py:585,670,681` |

`[M]` the low-order build consumes **no channel operator**. What it consumes is
`σ_s0^{g→g}` and `σ_s1^{g→g}` — the **within-group DIAGONALS** of the ℓ=0 and ℓ=1
kernels — as per-cell scalar coefficients forming `σ̂_R = σ_t − σ_s0^{g→g}` and
`D = 1/[3(σ_t − σ_s1^{g→g})]`. The group-transfer structure is deliberately
**discarded**; that discard is what "within-group acceleration" means.

⟹ **The DSA justification is REFUTED as stated, with a structural reason:** DSA
does not want the isotropic *operator*; it wants `diag(K_ℓ)`, which is a
**coefficient field** (`MultiplicationOperator`-shaped), not a member of the
transfer-kernel family. A factory ranging over `{IsoS, IsoN2N, IsoF}` cannot
produce it, because it is not in the range.

⟹ **and it converts into a real, smaller deliverable** (L-013 shape): the missing
verb is on the KERNEL, not on a factory — `ScatteringKernel.within_group()`
returning the per-ℓ diagonals `(L+1, ng)`. It replaces **two** facade accessors
with one kernel verb, deletes the `foldable_sigma` twin (F16/#16-a), makes DSA's
`if scattering_order >= 1` branch a slice, and is squarely inside F-1's scope.

⟹ **and the co-sourcing half of the justification is already true and already
stronger:** `S` and `S.isotropic_energy` are two bindings of ONE
`ScatteringMaterialField` instance. That is the aniso/iso pair the user is
reaching for. Nothing needs to be built.

---

## 2. ATTACK

### A1 — The factory manufactures a discrimination that does not exist

`[M]` F14: the iso family has **zero** runtime branches on space kind or channel
tag; every construction site knows statically which binding it wants, because the
consumer's own physics decides (diffusion has no angular space; SN's satellite is
the aniso operator's fast path). A factory keyed on "the problem" converts a
**statically-known fact** into a **runtime dispatch**. That runs
illegal-states-unrepresentable backwards: today
`IsotropicScattering(field, domain=s, codomain=s)` cannot be given an angular
space and quietly do the wrong thing; a factory returning
`LinearOperator` erases which member you got at the type level.

The `discriminations` smell says "a repeated conditional is a missing type."
There is no repeated conditional here. Adding a factory *creates* the conditional
and then names it the solution.

### A2 — The fiber is not a singleton, so there is no section to choose

Frame the factory as a section of `Op → Space`. Over a composite
`FullFieldSpace` the fiber contains **at least three** legitimate members of the
scattering channel alone: `ScatteringOperator` (the composite binding),
`S.isotropic_energy` (the ℓ=0 energy binding on the derived scalar sub-space), and
`S.foldable_part()` / `S.residual_part()` (the within-group split siblings,
`scattering.py:960,988`). Which one is "appropriate" is **not a function of the
space** — it is a function of the consumer's ROLE in the algebra (loss term / hot-
path accumulator / low-order coefficient). Role is not a datum any space carries.

⟹ the factory would need a `role=` argument, i.e. **stringly-typed dispatch on
the very axis §5 defers**. That is the god-factory risk in its precise form: not
"it will get big", but "its key does not exist yet, so it will be spelled as a
string."

### A3 — The trio is not a family: rank refutes it

`[M]` F8 + `fission.py:690,707`: F's energy action is **already single-sourced**
as `self.kernel` — a `TensorProductOperator` over `RankOneOperator` — and both the
`ScalarFlux` and `ndarray` arms delegate to it. So:

- "`IsotropicFission`" is **not a new object**. It is a NAME for an operator that
  already ships; the step-4 mint is a *binding* of it to a scalar space.
- It cannot join the `(matrix_of, scale)` family (F9) without **densifying**
  `outer(νΣf, χ)` into an `(ng,ng)` matrix — which costs `O(G²)` where `O(G)`
  suffices, destroys the rank-1 normal form
  ([[fission-rank1-normal-form-dead-functional]]: `F = |χ⟩⟨νΣf|` **is** the normal
  form), and replaces a **free** transpose (dual dyad, χ↔νΣf swap) with a
  `_TRANSPOSE` einsum.
- ⚠ And the factor swap **violates `FissionKernel`'s own simplex guard** — the
  adjoint of an emission kernel is not an emission kernel (already recorded at
  design-record §6). So F's transpose is structurally unlike S's and N2N's.

⟹ **the honest partition is 2 + 1, not 3.** `{IsoS, IsoN2N}` is one family (same
body, differing by a unit); `IsoF` is a rank-1 dyad that is a *sibling by role*
(an isotropic bulk emission on scalar space) and *not a sibling by construction*.
A factory over all three would hide the sharpest discriminator in the set behind
a uniform return type.

### A4 — Type-vs-property and defer-until-2, applied

`coding-standards`' decidable criterion: mint a type iff (a) ≥2 non-isomorphic
realizations, AND (b) a non-identity morphism is actually APPLIED between them.
For a **factory object** the analogue is: mint it iff it holds state that its
consumers cannot each derive, and iff ≥2 consumers reuse the same instance.
`[M]` today: **2 consumers reuse the same space instance at the within-group site**
(F10) — so the *pose* reading (1a) passes; **0 consumers reuse a channel-selection
decision** — so the *dispatcher* reading fails.

Defer-until-2 also **counts kinds, not instances** (L-011). The cited
justifications are three different kinds: DSA (refuted — wants a coefficient, not
an operator), sensors (a **functional**, `V → ℝ`, not in the operator family at
all — see A5), and cross-method reuse (already satisfied: the energy binding takes
a plain `FunctionSpace` and works on a CP/MoC per-FSR space today). No two pair up.

### A5 — The "sensors" justification puts a non-operator in an operator factory

A diagnostic functional's codomain is `ℝ` (or a scalar field), not the operand
space. `ReactionRateFunctional` and `FissionOperator.production_rate`
(`fission.py:337`) already exist, and F's dyad `|χ⟩⟨νΣf|` **is** (emission
vector) ⊗ (sensor) — the sensor is the dyad's left factor, already a first-class
object.

A "per-level sensor" is then `⟨w_ℓ ⊗ Σ_x, ·⟩` = **the frame's ℓ-th analysis row
composed with a reaction-rate weight**. `[M]` DSA already builds exactly this at
`dsa.py:584-586` (`quadrature.angular_frame(1).table[:,1,1]`, described in its own
comment as *"a CALLED single source, not a claim"*). ⟹ per-level flavours are a
**composition of two existing objects**, not a factory output; and putting them in
the same factory gives it a non-uniform codomain, which makes it a switchboard
rather than a family.

### A6 — A second spine beside the blessed frame chain

`[M]` §14.4's blessed chain has ONE spelling: `space → interior → angular axis →
generator_as(Quadrature) → hub → frame`. Any factory keyed on "the problem" must
re-derive the same pieces to answer its own key. F16 shows what happens when a
second derivation of one scalar is allowed: `W` currently has **four**. A factory
would be the fifth entry point into the chain and, unlike the four, an *entry
point by design*.

### A7 — The phantom-type trap, at the operator level

The type-vs-property corollary: *an axis that changes the ARITHMETIC INTERFACE
cannot be a phantom type parameter.* Its operator-level form: a factory whose
members have different transpose mechanisms (einsum swap vs factor swap, A3),
different cost classes (`O(G²)` vs `O(G)`), and different admissibility
(`_iso_is_assemblable` is honestly False for a scalar-space binding) cannot be
served by one generic body with order-/channel-branches. The negative
discriminator to watch for in any proposed implementation: **if it "works" only by
branching on a stored channel tag at runtime**, `replace(op, channel=Other)`
type-checks and walks straight through whatever gate the factory was minted to be.

### A8 — Is it the consumer redesign arriving early? Yes, and it is nameable

The factory's real content is `bind(·, space)` = a **pose** holding the space and
its derived pieces. That object is (i) the posing arc's subject, (ii) F-1's
subject on the data side (the `MaterialXSField` facade's remaining halves), and
(iii) CS5's subject on the generator side. All three are chartered elsewhere.
Building it inside step 4 would commit the consumers campaign's central object
before its own ground measure runs.

---

## 3. FRAME CANDIDATES

### Frame 1 — Graded / isotypic decomposition (representation theory)

**Trigger:** F6 + F7 — SO(3)-zonal scattering ⟹ Funk–Hecke grading by ℓ; all
three iso operators are the `V₀` (trivial-irrep) component of their channel.

**Reformulation:** the "iso version" is not a separate operator kind; it is the
`ℓ = 0` **component** of a graded family. The component-selection morphism is
`ScatteringMaterialField.truncated(0)` (`material_field.py:198`) — which already
ships and which the tier-2 mint already calls
(`isotropic_scattering.py:281`). The grading is on the **DATA**, and
`ScatteringMaterialField.order` **is** the binding's truncation (§4: *"the verb
signatures carry no `L` parameter"*).

**Elegance payoff:**
- *Structure-exposing:* names "order of anisotropy" as a property of the kernel
  field, not a factory parameter.
- *Structurally-simpler:* removes one of the three axes the proposed factory
  would key on.

**First test (discriminates):** assert that no production symbol in
`orpheus/transport/operators/` takes an `order`/`L` argument alongside a kernel
field. `[M]` `ScatteringOperator.from_solver_data` **still takes
`scattering_order`** (`scattering.py:765-770`) and immediately spends it on
`.truncated(scattering_order)` — so this test **REDs today at exactly one site**,
and that site is the honest residue of the grading not being fully data-side.
A design that adds a factory `order=` parameter makes it RED at N sites instead.

**Structural attack on current:** the tier-2 signature carries a number the datum
owns. The factory proposal would institutionalize that as an API axis.

---

### Frame 2 — Multiplier algebra / module over a scalar ring

**Trigger:** F9 — `IsotropicScattering` and `IsotropicN2N` differ, at every level
down to the einsum, by `(matrix_of, scale)`; `scale` is a **unit** (`1.0` vs
`ν₂ₙ`), and `matrix_of` is a **field-type accessor**.

**Reformulation:** collapse the two classes into one, with the channel's identity
living in the FIELD's TYPE, not in the operator's:

```
class IsotropicEnergyOperator(BoundOperator):
    field: MaterialField[Any]        # Scattering | N2N | (n,3n) | …
    # apply → field.add_source(out, x);  apply_transpose → field.add_source_transpose(out, x)
```
with `MaterialField` gaining the two abstract verbs `add_source` /
`add_source_transpose`, and the channel subclasses supplying `(matrix_of, scale)`
to the primitive that already exists (`_accumulate_contracted`). The domain-named
verbs (`add_p0_source`, `add_emission`) stay as the subclasses' vocabulary,
delegating to the base pair.

**Elegance payoff:**
- *Structurally-simpler:* 2 classes → 1; the three already-hoisted module bodies
  (`_scalar_composite_source`, `_iso_is_assemblable`,
  `_assemble_iso_energy_operator`) stop needing an `op: A | B` union annotation.
- *Expressive:* a fourth channel (an `(n,3n)` kernel; a per-ℓ diagonal removal
  operator) needs **no** new operator class.
- *Structure-exposing:* names the channel as *data*, which is exactly the §4
  kernel-collapse ruling extended one layer up.

**First test (discriminates, and includes its negative leg):**
1. `np.array_equal(IsotropicEnergyOperator(N2NMaterialField(...)).apply(φ),
   IsotropicN2N(...).apply(φ))` — **0 ULP**, both directions, on a `Σ₂ ≠ 0`
   fixture (bit-identity, not `allclose` — L-002).
2. **The negative leg that makes it fail-able:** a `FissionMaterialField`
   handed to `IsotropicEnergyOperator` must **RAISE** (it has no `(ng,ng)` matrix
   to contract; A3). An implementation that "passes" by densifying
   `outer(νΣf, χ)` returns the right numbers and is the design error — so the
   negative leg, not the value leg, is the real gate.
3. **The extensibility discriminator:** write an `N3NMaterialField` with
   `add_source(scale=3.0)` and assert an isotropic energy operator is
   constructible **with zero edits under `orpheus/transport/operators/`**. Today
   this REDs (needs a new class). Under a *factory* it also REDs (needs a registry
   row). It passes only under the field-typed collapse.

**Structural attack on current:** the arithmetic was single-sourced and the class
was not (Smell #16 shape 1, half-resolved). Today `IsotropicScattering` and
`IsotropicN2N` are the same discrete operator over two data types, and the union
annotations `"IsotropicScattering | IsotropicN2N"` on the three shared helpers
(`isotropic_scattering.py:145,168`) are the tell — a shared body that must
enumerate its callers is a missing type.

⚠ **Cost to weigh, not hide:** this trades two high-signal domain names for one
generic name plus a field type. `feedback_high_signal_names` cuts against it;
`ScatteringOperator.isotropic_energy` / `N2NOperator.energy` already give each
consumer a domain-named handle, which mitigates it. The name is a real fork for
the user, not a detail.

---

### Frame 3 — Low-rank / dyadic normal form

**Trigger:** F8 — `FissionKernel` is a factor pair; `F = |χ⟩⟨νΣf|` is rank 1.

**Reformulation:** `IsotropicFission` is the **scalar-space binding of the
already-existing `FissionOperator.kernel`** (`fission.py:281`, a
`TensorProductOperator` over `RankOneOperator`), not a third member of Frame 2's
family. Its ctor retains `(FissionMaterialField, domain, codomain)` and its
`apply` delegates to a per-material dyad — never a densified `(ng,ng)`.

**Elegance payoff:**
- *Structure-exposing:* makes the 2+1 partition explicit instead of implying a
  uniform trio.
- *Algorithmic-advantage:* `O(G)` per cell instead of `O(G²)`; the transpose is a
  factor swap (free) instead of a second einsum.
- *Structurally-simpler:* no new arithmetic at all — the kernel already computes it.

**First test (discriminates):** for a `ng ≥ 8` fixture, assert
`IsotropicFission.dense_per_material()[mid]` equals `np.outer(χ, νΣf)` **and**
that the `apply` path never materializes an `(ng,ng)` intermediate (count
`RankOneOperator.apply` calls, or assert peak per-cell temporary size is `O(G)`).
A "just densify it into the family" implementation passes the first clause and
fails the second — which is exactly the divergence being tested.

**Structural attack on current:** F's own energy operator already exists and is
unnamed at the binding layer, so the step-4 mint is a **naming and binding** act,
not a construction act. Treating it as "the third of a trio" would invite giving
it the family's dense body.

**Precedent:** [[fission-rank1-normal-form-dead-functional]] —
`F = |χ⟩⟨νΣf|` IS the normal form, so "unfold F" is structurally empty.

---

### Frame 4 — Fibration and its section (the factory as literally proposed)

**Trigger:** the proposal's own shape — a map `problem → operator`.

**Reformulation:** `Op → Space` with fiber = the operators bound to a space; the
factory is a chosen section `s : Space → Op`.

**Elegance payoff: none — the frame REFUTES the construction.** A section
requires a canonical choice in each fiber. `[M]` A2: the fiber over one
`FullFieldSpace` contains `ScatteringOperator`, `S.isotropic_energy`,
`S.foldable_part()`, `S.residual_part()` — four legitimate members of one channel.
The choice is made by the consumer's ROLE, which the base object does not carry,
so no section exists without adding role as an extra base coordinate. Adding it
means the base is `Space × Role`, and `Role` is precisely the posing concept §5
defers.

**First test (discriminates):** name the intended return of
`factory(space, channel="scattering")` for the SN within-group space. If the
answer requires a further argument to be well-defined, the section does not exist
and the factory is a `match` statement with a class name. This is a **design-time
test with a definite structural answer** — the L-011 "exhibit a coupled system
expressible nested but not flat" move.

**Structural attack on current:** none — the current design has no section
because it needs none; the consumer names the member it wants.

---

### Frame 5 — Stage-2 generator (Frame as factory), and its exact boundary

**Trigger:** F5 — an interned hub (`HarmonicFrame.for_space`) already induces
space + operators together at one site.

**Reformulation:** the project's factory object exists and its domain of
applicability is the **frame-bearing** bindings. Its boundary is the finding: the
iso trio carries zero angular content (§14.2), so extending the frame hub to them
would give it a second, non-frame responsibility — unless the extension is
through the **ℓ=0 sub-block**, which is Frame 6.

**Elegance payoff:**
- *Structure-exposing:* separates "what a factory legitimately owns here" (a
  frame's induced spaces + faces, interned) from "what it must not" (choosing
  which algebraic term a consumer wants).

**First test (discriminates):** assert `HarmonicFrame.for_space(interior, L)`
returns the **same object** (`is`) for two calls with the same `(quadrature, L)`
— the interning contract — and assert that constructing an `IsotropicScattering`
touches **no** frame at all (e.g. by patching `for_space` to raise and building
the iso binding successfully). The second leg REDs for any design that routes the
energy binding through the frame hub, which is the divergence.

**Structural attack on current:** nothing broken; this frame **bounds** the
factory rather than motivating it.

---

### Frame 6 — The LIFT as a functor: energy bindings → angular bindings ⭐

**Trigger:** F12 + Smell #16 shape 4 — `_per_ordinate.py:20-22` names its own
defer-until-2 trigger, and `[M]` the third channel has arrived and hand-inlines
the arithmetic (`fission.py:615-625`), with its own fourth spelling of `W`
(F16).

**Reformulation.** Every isotropic bulk emission's angular binding is the ℓ=0
frame conjugation of its energy binding:

```
N_iso = (1/W) · E₀ ∘ K ∘ M₀        M₀ = ∫dΩ  (the ℓ=0 analysis row)
                                    E₀ = isotropic broadcast (P₀ = 1)
```
i.e. `angular_frame(0).conjugate(K)` ([[iso-source-frame-conjugation-unification]]
Part 2 established this; `n2n.py:29` states it verbatim). So the map
`K ↦ lift(K)` is a **functor from the energy-binding category to the angular-
binding category**, and `N2NOperator` **is** `lift(IsotropicN2N)` wearing a class.

Concretely: promote `assemble_per_ordinate_isotropic` to the object its own
docstring predicts —

```
class IsotropicLift(BoundOperator):        # or: a shared mixin/ctor helper
    energy: BoundOperator                  # any scalar-space energy binding
    measure: <the ℓ=0 frame / the binding measure>   # THE single source of W
```
so that `N2NOperator`, `ScatteringOperator`'s P0 half, and `FissionOperator`'s
composite arm consume ONE lift and ONE `W`.

**Elegance payoff:**
- *Structure-exposing (strong):* names the relation between each channel's two
  bindings as a **functor**, which is what makes "the aniso and iso version of one
  channel" a theorem rather than two constructor calls — and it is the *inverse*
  direction from the proposed factory (one verb producing the sibling, not one
  dispatcher choosing between them).
- *Structurally-simpler:* collapses the third hand-rolled lift (`fission.py`) and
  reduces F16's four `W` spellings toward one.
- *Expressive:* a new isotropic channel gets its angular binding for free.
- *Algorithmic-advantage:* none claimed — the fast path is preserved verbatim.

**First test (discriminates, with a negative leg):**
1. `np.array_equal` between today's `FissionOperator.apply(psi).interior.values`
   and the lift-routed value on a fixture where `Σw ≠ N` (any slab GL rule: `[M]`
   `Σw = 2`, `N = 4` at S4, so the two are numerically distinguishable). **0 ULP.**
2. **The negative leg:** construct an angular axis whose `weights` are non-uniform
   and assert that F's current `w_total` expression
   (`ang.shape[0] if ang.weights is None else sum(ang.weights)`) and
   `S.total_weight` (`flux_analysis.frame.measure.weights.sum()`) are computed
   from **different objects** — a test that pins the two as one source REDs today
   at the source level even where the values agree. This is what makes the claim
   fail-able: agreement of values is not single-sourcing.
3. `assert F.total_weight is S.total_weight`'s producing object, or the weaker
   ULP-free form `F_lift_weight == S.total_weight` on a rule where
   `Σw ≠ N` — an implementation that keeps the `shape[0]` fallback passes only
   when the invariant it does not name happens to hold.

**Structural attack on current:** the shared primitive exists, its own docstring
predicts the trigger, the trigger has fired, and the third caller hand-rolls it
with a fourth derivation of the scalar it needs. That is Smell #16 shapes 1, 2 and
4 in one place.

---

## 4. CROSS-METHOD POLLINATION

**Current method class:** SN (with diffusion and homogeneous as live consumers of
the same energy bindings; CP / MoC / MC as latent consumers).

### B1 — From depletion / Bateman (channel accounting)

**Trigger:** F1 + F9 — a set of reaction channels whose per-cell action differs
only by (matrix, multiplier). Depletion's native object is a channel-typed
**datum** (MT-keyed reaction data), and its assembly is a plain sum over channels;
it does not carry one operator class per channel.

**Reformulation:** Frame 2's collapse, argued from an adjacent method: the channel
belongs in the FIELD's type, the summation in the consumer's `OperatorSum`.

**Payoff:** *structurally-simpler* + *expressive* (a new channel is a new kernel
type, not a new operator class).

**First test:** B1 = Frame 2's test 3 (the `N3NMaterialField` extensibility
discriminator, which separates the field-typed collapse from BOTH the status quo
and a factory).

### B2 — From multigrid / DSA (the Galerkin coarse operator)

**Trigger:** DSA reads `foldable_sigma` + `residual_sig_s` off the facade
(`dsa.py:263,266`), and `[M]` `foldable_sigma` **exists twice** — on the facade
(`material_xs_field.py:820`) and on the operator (`scattering.py:1028`).

**Reformulation:** the within-group extraction is the **Galerkin restriction of
the energy kernel onto the single-group subspace**, `K_gg = e_gᵀ K e_g` — the
multigrid coarse-operator move at `ℓ`-graded energy level. Its home is the
kernel: `ScatteringKernel.within_group() -> (L+1, ng)`. DSA then reads
`wg[0]` for `σ_s0^{g→g}` and `wg[1]` for `σ_s1^{g→g}`, its
`if scattering_order >= 1` branch becomes a slice, and the `foldable_sigma` twin
dies.

**Payoff:** *structurally-simpler* (2 facade accessors + 1 operator twin → 1
kernel verb); *structure-exposing* (names the extraction as a restriction, which
is why it is legitimate in the low-order build and forbidden in the sweep — the
`#215` trap the DSA docstring already explains in prose).

**First test (discriminates):**
1. `np.array_equal(kernel.within_group()[0], np.diag(kernel.p0))` and
   `[1] == np.diag(kernel.moments[1])` — 0 ULP.
2. `np.array_equal` on the whole assembled `a_low` and `g_map` before/after the
   re-point — the DSA build is bit-identical or the re-point is wrong.
3. **Negative leg:** the verb must return an object from which the **off-diagonal
   cannot be recovered** (a `(L+1, ng)` array, not a masked `(L+1, ng, ng)`), so a
   future consumer cannot mistake the restriction for the transfer kernel. A
   masked-matrix implementation passes (1) and (2) and fails (3) — that is the
   divergence.

### B3 — From CP / MoC (the latent consumers)

**Trigger:** the model-independence claim of the energy binding, already
established ([[iso-source-frame-conjugation-unification]] Part 3) and already
realized: `IsotropicScattering.from_material_xs(mat_xs, *, space: FunctionSpace)`
takes a **plain** `FunctionSpace`, with no `FullFieldSpace` guard —
unlike `ScatteringOperator.from_solver_data` and `N2NOperator.from_solver_data`,
which both raise on a non-composite space.

**Reformulation:** the CP/MoC per-region in-scatter is `IsotropicScattering` on a
per-FSR `FunctionSpace`. Nothing to build; the constraint is to **keep** the
energy binding free of composite-space knowledge.

**Payoff:** *expressive* (one operator serves four methods); it is also the
strongest single argument against the factory, because a factory whose key reads
`space.interior_space` **excludes** the very methods "a single object tending to
different transport methods" is meant to serve.

**First test (discriminates the designs, not the values):** construct
`IsotropicScattering.from_material_xs(mat_xs, space=fsr_space)` where `fsr_space`
has **no angular axis and no trace block**, and assert `.apply(φ)` works. It
**passes today**. Then assert that any proposed factory entry point raises
`TypeError` on the same input — `[M]` it must, since the composite guard
(`n2n.py:150-156`, `scattering.py:795-801`) is how every existing tier-2 method
derives its pieces. **The pass/raise pair is the discriminator**: it shows the
factory route is strictly narrower than the route that already ships.

---

## 5. THE DISCRIMINATING QUESTIONS (and where each answer lives)

| # | Question | Where the answer lives | Why it decides |
|---|---|---|---|
| Q1 | Is there ANY consumer that must choose the channel binding at **runtime**, from data? | The 4 sites of F14 + any new consumer's brief | If no ⟹ the dispatcher reading of the factory is minting a discrimination (A1). If yes ⟹ name the datum it dispatches on; that datum is the missing type. |
| Q2 | Does the "appropriate operator" depend on anything **other than** the space? | `scattering.py:555` (`isotropic_energy`), `:960/:988` (`foldable_part`/`residual_part`) | Four members in one fiber ⟹ no section ⟹ the key must include ROLE, which §5 defers (A2). |
| Q3 | Is the duplication the proposal senses the **extraction chain** or the **classes**? | `scattering.py:579,788-801` vs `n2n.py:145-157`; and `isotropic_scattering.py:145,168` union annotations | Extraction ⟹ two verbs on `FullFieldSpace` fix it with no factory (F15). Classes ⟹ Frame 2's collapse. `[M]` **both** are real, and they have different fixes. |
| Q4 | Can `IsotropicFission` share `IsotropicScattering`'s body without densifying? | `kernels.py:259-330` (`FissionKernel`, `dyad`), `fission.py:281,690,707` | No (A3, rank 1) ⟹ the set is 2+1, and any "trio" API hides that. |
| Q5 | Does DSA want an operator, or a coefficient? | `dsa.py:255-284, 326-338` | `σ̂_R = σ_t − σ_s0^{g→g}` and `D = 1/[3(σ_t − σ_s1^{g→g})]` are **scalars per (g, cell)** ⟹ coefficient. Refutes justification (a); yields B2. |
| Q6 | Is a "sensor" in the operator family? | `transport/reaction_rate_functional.py`, `fission.py:337`, `dsa.py:584-586` | Codomain `ℝ` ⟹ no ⟹ justification (b) is a **composition** question (frame ℓ-row ∘ reaction weight), not a factory question (A5). |
| Q7 | How many spellings of `W` does the tree carry, and do they share a source? | `scattering.py:506`, `n2n.py:163-172`, `fission.py:618-622`, `frame.py:766` | `[M]` **four**, three distinct sources ⟹ the ℓ=0 frame is the object that should own it (Frame 6). Any new entry point makes it five. |
| Q8 | Would the factory be constructible for a space with no angular axis and no trace? | `isotropic_scattering.py:275-284` (no guard) vs `n2n.py:150-156` (guard) | If no ⟹ the factory excludes CP/MoC, contradicting "a single object tending to different transport methods" (B3). |
| Q9 | Does `S` already carry the within-group split as OPERATORS, and is anything using them? | `scattering.py:960,988` + a caller grep | ✅ **MEASURED — see §6bis.** `[M]` **zero production call sites**; the operator-tier answer to "one channel, two flavours" already ships, is heavily gated, and nothing in `orpheus/` calls it. |
| Q10 | Is `order` an operator-API axis or a data property? | `material_field.py:186-210` (`order`, `truncated`) vs `scattering.py:765-770` (`scattering_order` param) | `[M]` data property everywhere **except** one tier-2 signature ⟹ a factory `order=` axis would re-open a closed ruling (Frame 1). |

---

## 6. REJECTED CANDIDATES — with the structural reason each failed

*(process-discipline: a refutation without its reason must be re-derived at full
cost. Each reason carries the QUESTION it was refuted for — L-021(b).)*

| candidate | refuted FOR | structural reason | the FACT that survives |
|---|---|---|---|
| **Factory-as-dispatcher** (`factory(problem) → operator`) | "should the trio get a factory?" | The fiber over a space is not a singleton (A2); the missing key is ROLE, and role is the posing concept §5 defers. | The fiber's four members are all real and all needed — that enumeration is the tree's actual channel-flavour API. |
| **`space.bind(kernel)`** | the same | Import direction: `numerics → transport` is a forbidden edge; a data object's verbs return arrays, only binders return operators (L-020). | The binder must live transport-side — which is why `HarmonicFrame.for_space` is where it is. |
| **Promoting `build_within_group_system` to a posing factory** | the same | It IS the consumers campaign (§5, explicitly deferred), and the word `posing` is already spent on the `(A_loss, M)` + eigenvalue-role axis of the same solve (L-012's spent-word kill). | The within-group site is the one place where ≥2 operators reuse one space instance — the only currying with a reused operand (1a). |
| **A channel REGISTRY / plugin table** | "how do we make the family extensible?" | The channel set is CLOSED and small, and members differ by kernel TYPE, not by a config string. A registry buys open-set extensibility for a set that grows by a dataclass. | The extensibility that IS wanted is Frame 2's: a new `MaterialField` subclass, zero operator edits. |
| **Uniform "trio" API over `{IsoS, IsoN2N, IsoF}`** | "are these three one family?" | RANK (F8): different normal form, different transpose mechanism, different cost class; F's factor swap even violates its own simplex guard. | The honest partition is **2 + 1**, and F's energy operator already ships as `FissionOperator.kernel`. |
| **Category theory / functor formalism as the deliverable** | "what frame explains the binding?" | The concrete win is already captured by §3's three-tier discipline and by the named lift (Frame 6). No abstract-nonsense lever adds a test. (L-001.) | The lift IS a functor, and saying so buys exactly one thing: the sibling-generation direction, which Frame 6 makes concrete. |
| **Tensor networks / MPO** | "is there a bond-dimension structure across channels?" | Bond dimension 1 (F is rank-1; the energy ⊗ angular ⊗ spatial factorization is a fixed 3-factor product with no truncation knob). Degenerate — not a network. (L-001.) | The 3-way factorization (F17) is real and is what makes the energy binding method-independent. |
| **Homology / chain complex over the block structure** | "does the bulk/trace split have a differential?" | No `∂² = 0`; interior ⊕ trace is a **biproduct**, and the iso operators are `BlockRole.BULK` with an implicit-zero trace — an inclusion, not a differential. (L-001.) | The three implicit-zero spellings (F13) are a real triple wanting ONE extension-by-zero verb — §3's Tier-3 lift, already ruled. |
| **Differential geometry / connection coefficients** | "is there a curvature term hiding in the binding?" | An energy operator is cell-diagonal by the material partition; no transport of a frame along anything. Fires for curvilinear angular redistribution only. (L-001.) | — |
| **Saddle-point / mixed FEM / inf-sup** | "does DSA's structure reach the channel operators?" | Fires on the **diffusion member** (the low-order edge system), never on the sweep or on an energy binding — the resolvent backbone says which member has the matching shape (L-007 corollary). | DSA's low-order system is where that apparatus attaches; the channel operators are outside it. |
| **Petrov–Galerkin projection (homogenization / condensation)** | "does the iso binding unify with the XS-coarsening frames?" | Different frame VERB: iso-source is `frame.conjugate` (Galerkin, SH = scattering's eigenbasis); condensation is `frame.project` = `G⁻¹M` (solution-weighted, owned by no operator). They do not unify. ([[iso-source-frame-conjugation-unification]], L-009/L-010.) | XD-9's χ↔νΣf-coupled condensation is a `project` question and stays on its own track. |
| **Fiber bundle / sheaf over the mesh** | "is `MaterialField` a bundle wanting bundle machinery?" | The kernel is **locally constant** on each material's cell set and the partition is disjoint ⟹ the bundle is trivial on each piece, no transition functions, no non-trivial section space. A `dict` + an index partition IS the whole structure. | `MaterialField[K]`'s current shape is exactly right and needs no geometric apparatus. |
| **Orbifold / symmetry quotient over the channel set** | "is there a group acting on the channels?" | No group acts: `{S, N₂ₙ, F}` are physically distinct reactions with no orbit structure. The quotient frame fires on quadrature node sets, not on channel sets. | — |
| **Feynman–Kac / path integral** | "is the trio a discretization of one path measure?" | The energy binding is a per-cell algebraic contraction with no path structure; the resolvent backbone's path reading attaches to `(Ω·∇+Σ_t)⁻¹`, not to the collision kernel's energy factor. | The backbone still predicts which layer is shared (the resolvent), and the energy binding is correctly outside it. |

---

## 6bis. Q9, measured — the operator-tier "two flavours of one channel" ALREADY SHIPS, uncalled

`[M]` `grep 'foldable_part|residual_part'` over `orpheus/` returns **15 lines and
zero call sites**: two definitions (`scattering.py:960,988`), one sibling
reference (`:1017`), one comment (`sn/operators/streaming.py:185`), and eleven
docstring cross-references. Every *invocation* in the repo is in `tests/`
(`test_scattering_operator.py` alone carries ~20). *(Predicate: literal-name grep
over `orpheus/`; a variable-call or registry-loop member would be invisible —
the 2026-08-29 §6b spelling lesson.)*

**Why this is the sharpest single datum in the probe.** `ScatteringOperator`
already answers "give me the isotropic/within-group flavour of this channel" at
the OPERATOR tier, in the shape the proposal wants (`S.foldable_part()` returns an
order-0 `ScatteringOperator` on the same bound ends, carrying only
`diag(Σ_{s,0})`; `S.residual_part()` carries the complement; the additive contract
`S ≈ S_fold + S_resid` is gated at `rtol=1e-14`). It has **no production
consumer**, and the two consumers who would want it — the σ_r-fold removal form
(#196, described at `streaming.py:185`) and DSA — do not call it.

⟹ the proposal's premise *"we probably need both flavours of the same channel
operator"* is **already implemented twice over** (this split, plus
`S.isotropic_energy`), and the open question is not *how do we produce flavours*
but *why does the shipped flavour machinery have no callers*. Answering that is a
consumer-side question (§5's campaign), which is A8 arriving from a second
direction.

### ⚠ A documentation bug found while measuring it (fix on sight, per the articulation standard)

`[M]` `material_xs_field.py:760` and `:780` state that `foldable_sigma` /
`residual_sig_s` are *"Consumed by `ScatteringOperator.foldable_part` /
`residual_part` to build a sibling operator"*. Since the CS4c step-3 rebind, both
sibling builders read `self.scattering.per_material` — the **kernel field** —
directly (`scattering.py:977-984`, `:1001-1010`); neither touches the facade.
Both sentences are **present-tense-false**. The only live consumer of the two
facade accessors is `dsa.py:263,266`.

The header at `material_xs_field.py:12-13` and the block comment at `:750`
(*"The scattering operator's foldable_part / residual_part / …"*) carry the same
stale attribution. This is F-1's scope and it is what makes B2 land cleanly: the
accessors' stated owner is gone, their real owner is DSA, and their mathematical
content is a kernel restriction (B2).

---

## 7. What the tension points at (not a verdict)

The proposal's instinct is picking up **two real duplications** and attributing
both to a missing factory:

- **F15** the tier-2 extraction chain, duplicated per class (2 guards, 2
  scalar-sub-space derivations, 2 tier names) — whose native fix is **two verbs on
  `FullFieldSpace`**, not a factory;
- **Smell #16 shapes 1/2/4** — the iso pair's twin classes (Frame 2), the four
  `W` spellings and the third hand-rolled lift (Frame 6), and the
  `foldable_sigma` twin (B2).

And it offers **two justifications that the tree refutes**: DSA consumes no
channel operator (it wants `diag(K_ℓ)`, a coefficient — Q5/B2), and a sensor is a
functional, not an operator (Q6/A5). The aniso/iso co-sourcing the factory would
provide **already ships**, by construction rather than by convention
(`S.isotropic_energy`).

The direction the frames point is the **inverse** of the proposal: not one
dispatcher choosing among bindings, but **one functor generating a channel's
angular binding from its energy binding** (Frame 6), over **one collapsed energy
operator class whose channel lives in its field's type** (Frame 2), with
**fission held out by rank** (Frame 3).

And §6bis adds the datum that reframes the whole request: the operator-tier
flavour machinery **already ships and has zero callers**. The productive question
is not how to produce flavours but why the shipped ones are unconsumed — which is
a consumer-side question, i.e. A8 from a second direction.

Every claim above is `[M]` by direct file read at the paths cited, on a worktree
whose Nexus graph is stale at `92dcc30f`. All ten discriminating questions carry
either a measurement or the exact file to measure it in.
