# Cross-Domain Attacker — Lessons

Behavioral/process lessons only: "what detection mistake did I make, or what
recurring insight sharpened my frame-attacks?" The frame-trigger CATALOG lives
in the `cross-domain-frames` skill (Part A/B/C) — do NOT duplicate it here.
Durable frame-matches that became real architecture are DESIGN POINTERS in
`MEMORY.md`, not lessons.

The cross-cutting meta-pattern behind most of these: **a frame-attack's value
is not "I named an exotic frame" — it is "I produced a concrete reformulation
with a fail-able first test, OR I crisply refuted the frame with a reason that
survives into the trigger table as a non-entry."** Speculation that names a
frame without a payoff degrades the table's signal. Every lesson is one face of
that standard.

---

## L-001 -- Refuted frames are first-class output; record the REASON, not just the rejection

→ The DIRECTIVE is now in AGENT.md (Required Output Shape — "A refuted frame is
first-class output"). It is identity-level (every attack writes an UNEXPLORED
block). What stays HERE is the forensic catalog: the specific high-prior frames
that keep getting correctly refuted on transport work, with their structural
reasons — recalled when a fresh attack is tempted by one of them.

The recurring high-prior frames that keep getting (correctly) refuted on
transport work, with their reasons:

- **Wiener-Hopf factorization** — wrong solver FAMILY. It is native to the
  Chandrasekhar/H-function half-space line, structurally incompatible with a
  bouncing-Peierls or sweep formulation. Keeping the two families structurally
  independent is itself a V&V requirement (independent references).
- **Homology / chain complex** — tempting via the word "boundary," but
  `∂∘∂ ≠ 0` in transport (two reflections compose to a non-trivial map; the
  boundary trace + its extension are a dagger adjoint PAIR, not a differential).
  No `∂²=0` ⇒ no homology payoff.
- **Category theory / operad / PROP** — almost always LOW-SIGNAL: the concrete
  win it gestures at (role-parameterization, compositional structure) is already
  captured by a nameable concrete frame (biproduct, affine torsor, forgetful
  functor with explicit laws). Name the concrete frame; list category theory
  UNEXPLORED unless a specific functor/law produces a test.
- **Tensor networks / MPO** — fires ONLY on a genuine bond-dimension trigger
  (a rank-N chain where N is a real truncation knob). A rank-1 or rank-2 fixed
  structure (a biproduct, a 2-surface BIE) is bond-dimension-1/2 DEGENERATE —
  not a network. Do not promote it to MPO until N≥3 actually ships.
- **Differential geometry / Christoffel** — needs a CURVATURE term to
  redistribute. Straight Euclidean chord segments or a Cartesian cell have none.
  It fires for the curvilinear streaming redistribution `(1−µ²)/r ∂_µ`, NOT for
  geometry-of-the-domain questions.

How to apply: when a high-prior exotic frame does not fire, write the one-line
STRUCTURAL reason (wrong family / no `∂²=0` / degenerate rank / no concrete law)
into UNEXPLORED. A bare "category theory — no trigger" is weaker than "category
theory — role-parameterization win already captured concretely by affine+Krylov;
no abstract-nonsense lever needed."

---

## L-002 -- A first test that cannot fail is rejected output — make it DISCRIMINATE

→ The DIRECTIVE is now in AGENT.md (Required Output Shape — "A first test that
cannot fail is rejected output"). It is identity-level (every frame candidate and
pollination emits a first test). What stays HERE is the forensic detail: the
specific discriminator constructions that have worked, recalled when building a
first test for a new claim.

A real first test discriminates the reformulation from the status quo by being
able to RED. The discriminator constructions that have worked:

- Frame the test as a property the NATIVE frame predicts and a wrong/naive
  implementation VIOLATES. Multiplier-algebra: `M_f @ M_g == M_{f·g}`
  bit-identical (a "wrapper that just stores an array" fails). Transpose-of-a-
  sweep: build the dense `L`, take `L.T`, assert the reverse-walk recurrence
  reproduces it (a spatial-only reverse fails on the nested angular block).
- For a refactor claimed bit-identical, the discriminator is `array_equal`
  (0 ULP), NOT `allclose` — only bit-identity distinguishes a genuine
  single-source from a value-correct-by-coincidence twin.
- For a typing claim, the discriminator is a NEGATIVE test: `np.ndarray`
  satisfies `Vector` but NOT `TransportState`; `tff_flux + tff_moment` RAISES.
  A test where everything passes proves the bound is theatrics.

How to apply: before emitting a first test, ask "what implementation would this
test PASS that I am claiming is wrong?" If the answer is "none," the test cannot
fail — rewrite it to target the specific divergence (the dropped term, the wrong
metric, the un-transposed nested recurrence).

---

## L-003 -- Smell #16 (distinct paths/reps to one operator) is the dominant transport tell — fire all four shapes

The single most-recurring native-frame-not-found signal in this project's
SN/transport work. The four-shape CATALOG is NOT re-copied here: it lives in the
skill (reference.md Part C, Smell #16) and in the AGENT.md promoted kernel, both
preloaded on every dispatch. The LESSON is the detection discipline, not the
catalog.

How to apply: every shape resolves to the SAME elegance move — collapse the
distinct paths onto ONE primary object (faces / a named trace-or-multiplication
operator / a displacement type / a shared primitive), which turns a correctness
COINCIDENCE into a theorem and usually deletes a marshalling shim. When you spot
the smell, name WHICH shape — the fix differs per shape.

---

## L-004 -- Property-vs-TYPE is decidable, not a taste call: demand a coexisting dual + an APPLIED morphism

A recurring design question is "should X earn its own field/state TYPE (vs a
PROPERTY/parameter on an existing type)?" The durable criterion that resolves it
without an unbounded taste argument:

> A representation earns a distinct TYPE iff there exist **≥2 bases that are NOT
> canonically isomorphic** (the iso depends on a quadrature/node choice),
> connected by a **change-of-basis operator that is itself MODELED and APPLIED**
> — it carries truncation error, has an adjoint, and participates in the
> operator algebra. All three clauses must hold.

Worked: angular order PASSES (ordinate `AngularFlux` ↔ harmonic moment field,
bridged by the applied projection/reconstruction Vandermonde-like pair) → two
types, correct. Spatial order FAILS clause 1 (one tensor-Legendre basis; the
only morphism is identity) → a PROPERTY, correct. If the modeled morphism would
be `id`, the type's sole payoff (forbid-mixing-a-dual) has no referent →
type-theatrics. Decidable by grep: count the within-X representations and the
applied non-identity morphisms between them.

Corollary — defer-with-an-EXPLICIT-trigger: when no method supplies the dual
TODAY but one plausibly arrives (here: nodal-DG / nodal-diffusion would supply a
nodal↔modal morphism), record the precise condition that flips the verdict.
"No current consumer" is not "never" — name the latent consumer and the trigger.

How to apply: when asked property-vs-type, do not argue aesthetics — count the
coexisting non-iso reps and the applied morphisms. Zero applied non-id morphisms
⇒ PROPERTY. Pairs with the project's unify-after-two-instances rule.

Corollary — an axis that changes the ARITHMETIC INTERFACE cannot be a phantom
type PARAMETER; it MUST be a distinct CLASS. A `Generic[Tag]` parameter is erased
at runtime and does NOT specialize dunders — so two instantiations
(`Field[Rep,Flux]` vs `Field[Rep,Source]`) share ONE `__add__` body. If the two
"values" of the axis have DIFFERENT `__add__` signatures (a torsor `A×V→A` that
forbids `A×A`, vs a vector `V×V→V`), no shared body satisfies both — the only
encoding is a distinct class per value (a mixin). This is a HARD refutation, not
a taste call: it killed the phantom `Field[Rep,Role]` carrier outright (the
(Rep×Role) grid attack). The decision lattice for a 2-axis carrier grid:
axis-changes-arithmetic ⇒ class; axis-changes-SHAPE ⇒ class; only an axis that
changes NEITHER (a true index/tag) can be a phantom param. When BOTH axes change
arithmetic-or-shape, the unique elegant form is the orthogonal-factor MULTIPLE
INHERITANCE `Leaf(RoleMixin, RepBase)` — and that parametrization, if wanted on a
type at all, belongs on the OPERATOR contract `[Din,Cout]` (where the axis values
are leaf TYPES, the genericity is APPLIED, and role-preservation is a fibration
theorem), NOT duplicated onto the carrier. The NEGATIVE discriminator: a
phantom-param impl that "passes" only by branching on a stored `role` field at
runtime is the stringly-typed anti-pattern — `replace(f, role=Other)` type-checks
and bypasses the gate; that bypass test REDs it.

---

## L-005 -- Read the WORKTREE, distrust Nexus on a feature branch

Frame-attacks on active branches repeatedly grounded on STALE facts because
Nexus answers from the main checkout's graph and the live code is on the branch.
Every design memo that grounded its frames on Nexus `context`/`query` while the
work lived on a worktree risked citing a superseded file:line.

How to apply: on a feature branch, ground frame-attacks on file:line read DIRECTLY
from the worktree (Read/Grep on the absolute branch path), and say so in the memo
("branch-verified, NOT Nexus — stale"). Reserve Nexus structural queries for
main-checkout questions or after `use_workspace(<worktree root>)` against a graph
built inside the worktree. A frame whose trigger is a code fact is only as good as
the freshness of that fact.

---

## L-006 -- "Frame-leak naming" — a model-agnostic interface named after ONE consumer is a latent lie

Detection-adjacent vocabulary insight that recurs on shared interfaces: a slot on
a model-AGNOSTIC layer named after one consumer's physics (`total_xs` on an
advection–reaction closure a diffusion solver will also consume;
`is_scan_march_compatible` — a SCHEME property named after a sweep STRATEGY)
becomes a lie the moment the second consumer passes a different realization
through the same slot. TELL: a docstring that says "generic in X" while the
parameter is named after a specific X₁. FIX: name the ROLE in the INTERSECTION of
all consumers' domains (`reaction_xs`; `transverse_coupling_is_facewise`), not the
realization. Distinct from Smell #16 (that is two reps of one quantity; this is
one slot whose NAME over-commits to one of N consumers).

How to apply: when assessing a name on a multi-consumer interface, ask "what does
the SECOND consumer call this?" The decisive first test is a second consumer that
reads the property with NO first consumer in scope — if it can't, the name is
strategy-entangled. SKILL-PROMOTION STATUS (re-checked 2026-06-22): HELD for a
THIRD sighting. Current count = TWO independent (`total_xs` on an
advection–reaction closure a diffusion solver also consumes; `is_scan_march_compatible`
= a SCHEME property named after the ScanMarch STRATEGY, the #240 D5 trait). The
project floor already covers the GENERIC vice via `coding-elegance` ("frame-leak
parameter naming"); a Part C SMELL earns its slot only when the cross-domain
detection angle (a 2nd-consumer-with-no-1st-in-scope first test) has a third
sighting distinct from the two naming cases. Until then, fire it inline, do not
promote.

---

## L-007 -- The transport-resolvent backbone predicts cross-method layering AND its exceptions — reach for it first

The spine itself is the AGENT.md kernel ("Cross-method backbone: the transport
resolvent") — preloaded, so it is NOT restated here. The LESSON is how to DEPLOY it:

- When a "find-the-special-value" family (k, α, time-step, fixed-source) looks like
  distinct solvers, check whether they are POSINGS of one generalized eigenproblem
  `Aψ=λMψ` sharing the resolvent backbone. If so, the ONLY genuinely per-method
  layer is the loss-operator REALIZATION; role-assignment, the µ→physical map,
  adjoint, and transient-shift are method-agnostic data over a shared engine.
- Generality flows toward the OPAQUE resolvent interface, never away from it. When
  two iterative drivers share a loop body, the one exposing its resolvent BEHIND a
  Protocol is the general engine; the one built from a concrete `(L,S,F)`
  factorization is the specialization that adapts INTO it. A "retire the opaque
  loop, migrate everyone to the concrete loop" plan points the deprecation arrow
  the wrong way (a CP/diffusion matrix has no `(L−S)⁻¹` to factor out).

How to apply: open any transport solve/adjoint/eigenvalue/layering question with
the resolvent backbone — it predicts both the layering split and the diffusion
exception from ONE principle, and it tells you which layer is the shared engine
before you read a line of the specific drivers.

Corollary — the backbone tells you WHERE a foreign frame fires. A cross-domain
frame keyed to an operator's ALGEBRAIC SHAPE fires only on the members whose
shape matches. Worked (DSA ↔ mixed-FEM/CFD saddle-point,
[[dsa-saddle-point-mixed-fem-frames]]): the **saddle-point / inf-sup / mixed-FEM**
frame fires on the **diffusion / low-order member ONLY** — because the mixed
`[[A,Bᵀ],[B,C]]` structure IS the elliptic exception (diffusion is
self-adjoint ⇒ mixed-when-first-ordered), while SN/MoC/CP are
characteristic-triangular SWEEPS with no saddle to stabilize (the primary
transport operator is either `L⁻¹` triangular or the Peierls `I−PL⁻¹Σs`
compact-perturbation-of-identity — neither a saddle point). So "consistent DSA"
= "the low-order is the **Schur complement of a compatible pairing**," and the
whole mixed-FEM/CFD apparatus (inf-sup, Darcy-vs-Stokes, Rhie–Chow, block
preconditioners) attaches to the acceleration subproblem, never to the sweep.
Before pointing a foreign frame at "transport," ask the backbone WHICH member
has the matching shape — pointing inf-sup theory at the sweep is a category
error the backbone catches for free.

---

## L-008 -- A "fully probes" claim is about operator LINEARITY, not input polynomial degree

A recurring MMS/verification reformulation trap: assuming a richer (higher-degree)
input is needed to "more fully" exercise an operator. For a LINEAR operator, an
input that is merely NON-CONSTANT in the operator's active variable already probes
the FULL map; a higher-degree input only changes WHICH point in the already-fully-
probed range you land on. Worked: the curvilinear angular redistribution
`(1−µ²)/r ∂_µ` is linear in ψ, so a linear-in-µ ansatz `(A+µB)/W` (the native
truncated-P1 Legendre element) fully activates it — no P2 needed. Enrich the
ansatz degree ONLY to satisfy a QUADRATURE-exactness requirement (e.g. `Σ wₙµₙ²`),
never to "more fully" probe a linear closure.

Paired HAZARD (the larger correctness risk when lifting a Cartesian reference to
curvilinear): a redistribution/curvature term carries a `1/r`, so a curvilinear
slope driver MUST vanish at the origin (`B(0)=0`) for pole-regularity; the slab has
no such constraint, so a slab-derived ansatz silently drops BOTH the term itself
AND its regularity constraint. The geometry MEASURE also enters the L2 error norm,
not just the source — an unweighted norm mis-measures the convergence order.

How to apply: before enriching an MMS ansatz, check the operator's linearity. If
linear, a non-constant input suffices — spend the degree budget on quadrature
exactness, and on a curvilinear geometry check the `1/r` pole-regularity of every
redistribution term and the measure-weighting of the error norm.

---

## L-009 -- A change-of-basis frame's OWNER and its Galerkin-vs-PG discipline are predicted by the operator's SYMMETRY (commutant membership / Funk–Hecke), not by which subsystem calls it first

A recurring architectural question on this project: when a method projects to
coefficients, acts there, reconstructs (`R∘A∘M`), WHO owns the frame `(M,R)` and is
it Galerkin or Petrov-Galerkin? The durable detection kernel — distinct from the
resolvent backbone (L-007, which is about solve/iteration LAYERING; this is about
projection-frame OWNERSHIP and DISCIPLINE):

> A frame `(M,R)` is OWNED by the operator `A` whose EIGENBASIS it is, and it is
> GALERKIN iff that eigenbasis is ORTHOGONAL — both decided by `A`'s symmetry. If
> `A` commutes with a group action (is in the commutant), Schur's lemma forces it
> block-diagonal-per-irrep in the isotypic basis, that basis IS `M`'s codomain, and
> a SELF-ADJOINT `A` (real kernel) diagonalizes ORTHOGONALLY ⟹ M*=R up to the
> Plancherel metric ⟹ Galerkin. No symmetry ⟹ no eigenbasis ⟹ the frame is a
> SOLUTION-WEIGHTED projection (test≠trial) ⟹ Petrov-Galerkin, owned by no operator.

Worked (the angular SH frame): Σ_s(Ω·Ω') is SO(3)-zonal ⟹ Funk–Hecke diagonalizes it
in {Y_ℓ^m} with eigenvalues = the Legendre moments (= the diagonal of the in-code
`Λ`); so M is LITERALLY the change-of-basis into scattering's eigenspace, the frame
is scattering-OWNED, and it is a `GalerkinFrame` BECAUSE Σ_s is self-adjoint-zonal
(orthogonal eigenbasis). Streaming `Ω·∇` is the ℓ=1 tensor operator (Clebsch–Gordan
⟹ ℓ↔ℓ±1 PN recurrence), does NOT diagonalize, so it does NOT own the basis. The
DISANALOGY that confirms the rule: energy condensation's G×G group-transfer matrix
has NO symmetry / no Funk–Hecke ⟹ its frame is a flux-weighted `PetrovGalerkinFrame`,
owned by no operator. ONE principle (operator symmetry) thus explains an entire
campaign's Galerkin-vs-PG split that prior memos had ASSERTED axis-by-axis.

Two corollaries that drop out:
- **Falsifiability of "subsystem X owns the frame":** the claim is structurally
  CONFIRMED (not non-falsifiable) when X's operator is the one whose eigenbasis the
  frame is. The genuine falsifier is a SECOND consumer whose TRUNCATION ORDER is set
  independently of X's operator (an output detector-functional of order L_d, or a
  flux expansion L_flux ≠ X's order) — that consumer makes the frame a general
  L²-tool with ≥2 independent consumers, flipping ownership. "Any function is
  X-basis-expandable" is NOT such a falsifier: the INFINITE expansion is basis-
  agnostic, but the TRUNCATED frame the code actually has is dimensioned by X's
  spectrum support (the operator's moments vanish above its order).
- **Placement:** the eigenbasis-owner is the canonical CONSTRUCTOR + the L-binding,
  NOT a private field — the generic frame machinery (analysis/reconstruct/conjugate)
  stays in the neutral layer (shared with the no-symmetry PG consumers); only the
  CONSTRUCTOR `owner.frame = neutral_factory(owner_order)` records ownership, and it
  relocates to the neutral factory the instant a second independent-L consumer lands.

How to apply: for any `R∘A∘M` ownership/discipline question, ask "what symmetry does
A have?" before reading call sites. Rotationally-invariant/zonal/convolution kernel
⟹ Funk–Hecke/Schur eigenbasis ⟹ A owns a GALERKIN frame. No symmetry ⟹ solution-
weighted PETROV-GALERKIN, owned by none. The first test that discriminates: assert
the owned frame is `GalerkinFrame` (M*=R up to the Plancherel/Gram metric) while a
no-symmetry sibling is a genuine `PetrovGalerkinFrame` (M*≠R, test=solution·trial).

Corollary — the SYMMETRY-SUB-BLOCK + multigrid-coarse-operator face (DSA #2, the ℓ=0
frame; [[dsa-rp-angular-frame]]): the Galerkin verdict descends to an IRREP SUB-BLOCK.
An acceleration/coarsening projecting onto ONE symmetry sub-block (DSA → the ℓ=0 / V₀
trivial-SO(3)-irrep constant on S²) inherits Galerkin from the parent symmetry-owned
frame — it is `angular_frame(0)`, NOT a new `ConstantBasis`, and the R/P pair is that
sub-frame's two faces (Π=P∘M W-self-adjoint under the PLAIN measure; a solution weight
would be the ONLY thing making it PG, and DSA has none). The multigrid connection this
adds: a "consistent low-order operator" IS the **Galerkin coarse operator** `R A_high P`
of the sub-block frame, post-composed with a **Schur complement** of the
retained-but-closed moments (Fick = odd-block Schur; Marshak = incoming-partial-current
Schur — the SAME move interior vs boundary), and "consistent" means that triple product
is taken on the DISCRETE (assembled) high-order operator (reduce-discrete ≠
discretize-reduced). One symmetry principle now predicts frame OWNERSHIP, the
Galerkin-vs-PG DISCIPLINE, AND the multigrid CONSISTENCY condition.

Corollary — the THIRD outcome: **ALL owners ⟹ NO owner, and it is still GALERKIN**
(#326 symmetry quotient, [[quadrature-symmetry-quotient-frames]]). The rule as stated
has two outcomes (one operator's eigenbasis ⟹ owned+Galerkin; no symmetry ⟹ owned by
none+PG). A **symmetry quotient** is the third: the group sits in the commutant of
**every** equivariant operator at once (streaming, collision, scattering, fission AND
the BCs are all `C_{2v}`-equivariant for a 1-D cylinder), so by Schur they are all
simultaneously block-diagonal in the isotypic decomposition. A frame owned by
everything is owned by nothing — it belongs to the **PROBLEM's symmetry**, not to an
operator. Crucially this does NOT make it Petrov-Galerkin: PG in this project means a
**solution weighting**, and a symmetry fold carries none (its Gram stays diagonal at
exactly the parent value, because invariant functions are constant on orbits and the
orbit weights sum to the parent's). So the verdict is **Galerkin on a SMALLER SPACE**
— second sighting of that exact verdict after the DSA ℓ=0 sub-block. Detection rule:
before typing an `R∘A∘M` as PG, ask *"is the test≠trial gap a SOLUTION weight, or a
GROUP identification?"* A group identification is Galerkin-on-a-sub-block; typing it
PG repeats the category error B3.0 fixed when it moved the Lambertian out of the
geometry slot.

SKILL-PROMOTION STATUS: a STRONG candidate for skill Part C (a new smell:
"eigenbasis-blind frame placement" / "operational-pipeline vocabulary for a spectral
decomposition") — the `harmonic_moment_flux.py:6` "natural data carrier of the
Galerkin pipeline" is the tell that the native Funk–Hecke frame is unnamed. The DSA
ℓ=0 case is a second CONSUMER of the SAME angular frame (not the independent
non-angular sighting the bar wants), but it strengthens the rule with the sub-block +
multigrid face above; the #326 symmetry-quotient case adds the third outcome and a
SECOND instance of the Galerkin-on-a-sub-block verdict. Still held for a genuinely
non-angular eigenbasis frame before promotion; until then fire it inline.

---

## L-010 -- A conserved-quantity COLLAPSE splits by WHAT is conserved (rate vs probability/mass), which fixes the MORPHISM (average vs marginalize) — NOT by a weight

When two coarsening/reduction operations look like "the same projection with vs
without a weight" (a 1-frame-vs-2-frame asymmetry, a "bare sum vs weighted
average" asymmetry), DO NOT accept the weight framing. Ask first: **what
functional does each collapse preserve?** A reaction RATE `⟨T·w, Σ⟩` is preserved
by an AVERAGE = `G⁻¹·M` (the projection `frame.project`, normalize=True). A
PROBABILITY or MASS (`Σχ=1`; a particle count) is preserved by a MARGINALIZE =
`M` alone (the un-normalized analysis `frame.analysis` / a bare `@T` against a
partition-of-unity table, normalize=False). These are DIFFERENT MORPHISMS that
differ by the `G⁻¹` factor — a weight=1 `project` would divide by the bin COUNT
and BREAK `Σχ=1`, so the "weight=1 degenerate of project" framing is provably
wrong. The honest unification is ONE machinery `(test_weight, normalize?)`:
`average = analysis ∘ G⁻¹` vs `marginalize = analysis`. Exposing both DISSOLVES
the frame-count asymmetry — it was never about how many frames, it was about
whether each channel's collapse axis carries a conserved RATE or a conserved
MASS.

Two corollaries that recurred in the same attack (XS coarsening: spatial
homogenize ∥ energy condense):
- **A "same slot ± weight" comparison can be hiding an AXIS category error.** χ
  in spatial homogenization collapses the SPATIAL axis (average); χ in energy
  condensation collapses the BIRTH-ENERGY axis (marginalize). Comparing them as
  one slot conflates two operations on orthogonal axes. Before unifying two
  collapses, confirm they act on the SAME axis; if not, the "asymmetry" is just
  two different reductions wearing the same channel name.
- **A precondition spelled as a 30-line docstring caveat on a 3-line method body
  wants to be a TYPE.** `FrameBase.gram` hardcodes a row-sum probe valid only for
  disjoint (diagonal Gram) or partition-of-unity (`R·1=1`) bases, then documents
  the third (tapered/dense) case it silently gets wrong. The Gram structure is a
  property of the BASIS (declare DIAGONAL / POU / DENSE), and `project` should
  dispatch on the declaration and RAISE on the unhandled case — the same
  "no-consumer ⟹ raise, don't silently delegate a half-consistent op" discipline
  a test-only basis already uses for its unbuilt synthesis side. A silent wrong
  number is the landmine; the declared-type + negative-test (a DENSE-declaring
  stub makes `.project` RAISE) closes it.

How to apply: at any "reduce/coarsen/collapse a container, preserving X"
question, name X (rate? probability? mass? current?) for EACH channel before
reaching for a frame. Rate → average (`G⁻¹M`); probability/mass → marginalize
(`M`); a second functional (surface current, leakage — GET/Smith) → a second
test space. The morphism follows from the conserved functional, and the
discriminating first test is the order-non-commutativity of a multi-axis channel
(`project(Σ@T) ≠ (project Σ)@T` because the normalization is keyed on one axis).
SKILL-PROMOTION STATUS: strong Part C candidate ("collapse-morphism-blind:
treating a marginalization as a weight=1 average"). Held for a SECOND sighting
(a non-XS conserved-collapse — e.g. a probability/measure reduction in MC tally
binning or a flux-to-current marginalization) before promotion; fire inline
until then.

---

## L-011 -- A "coupled / nested block system" proposal is a FREE RE-ASSOCIATION of an existing biproduct, not a new object — and "N instances justify the machinery" fails unless the N share a coupling KIND

Two reusable detection moves fired together on the augmented-SN "coupled 2×2
[[A_AA,A_AB],[A_BA,A_BB]], system=field+BC" adjudication
([[coupled-system-field-bc-frames]]):

- **Mat∘Mat≅Mat: a coupled 2×2-of-subsystems over a direct-sum carrier that
  ALREADY carries a biproduct block algebra is that biproduct RE-PARTITIONED,
  not a new categorical object.** The biproduct `⊕` is coherently associative, so
  grouping a flat N-block composite into 2 subsystems (`Mat₂(Mat₂(𝒞))≅Mat₄(𝒞)`)
  is free; the off-diagonals were always there (the seed + the −B boundary block).
  The G-adjoint composes block-wise for free when G is block-diagonal per subsystem
  (`A†` reads `G⁻¹AᵀG` at ANY partition granularity). TELL: a proposal says a new
  type "sits above" the existing block algebra — analogy-language for a theorem
  (same shape as issue_208's "natural 2×2 / adjoint for free"). DO NOT mint the
  `CoupledOperator`/nested type: it is a VIEW (redundant) or a twin (Smell #16). The
  discriminating first test is a CHALLENGE with a definite structural answer: "exhibit
  a LINEAR coupled system expressible nested but NOT flat" — impossible; every
  candidate is flat-re-expressible (⇒ view) or nonlinear (⇒ not a LinearOperator, a
  DIFFERENT abstraction). Before accepting a "coupled/nested/N-way" type, check whether
  the base algebra's block INDEX is merely FROZEN (here `BlockRole` = a 3-value enum
  while `_join_block_roles` already treats a role as a set-of-touched-blocks) — if so,
  the minimal object is "lift the freeze to N-way," not "a new layer above."

- **Defer-until-≥2 counts KINDS of the structure, not INSTANCES of the word.** A
  build-now case citing N coupled-system instances (ψ½ / DSA / multiphysics) collapses
  the moment you classify their coupling STRUCTURE: ψ½ = linear/triangular/metric-adjoint
  off-diagonals; DSA = linear/two-way-iterative/R⊣P-Galerkin off-diagonals; multiphysics
  = NONLINEAR/fixed-point (not a Mat(𝒞) block matrix at all). Three different kinds ⇒ no
  two pair up ⇒ the general machinery has no second instance to generalize FROM, and the
  nonlinear one UNDER-reaches a linear block abstraction (drop it from the count). The
  over-reach dual: a metric/triangularity/PSD assumption baked from the FIRST kind
  (ψ½ triangular biproduct, PSD block-diag metric) EXCLUDES the others (DSA two-way; RQI
  KKT-indefinite-zero-corner). Each kind gets its own home + trigger (ψ½=biproduct-exists;
  DSA=coupled-iterative-defer, the R⊣P shape DEFINES it when it lands; RQI=saddle-point-defer).

How to apply: at any "unify these coupled/composite subsystems under one new type"
request, (1) ask if the carrier already has a biproduct — if yes, the coupling is a
re-association, name the off-diagonals + lift the block-index freeze, do NOT add a layer;
(2) tabulate each cited instance's (off-diagonal structure, metric definiteness, solve
kind, linear?) — build only where ≥2 rows MATCH, defer the rest with the row that will
define them. Pairs with L-004 (property-vs-type by applied morphisms) and L-007 (the
resolvent backbone predicts which layer is shared).

---

## L-012 -- A NAMING task is a frame-detection task: a family word can be FORBIDDEN by a theorem (⇒ species + genus-ABC), and a name's first test is the invariant the name PROMISES run against the object that VIOLATES it

Naming requests ("find the faithful name for these N sibling objects") are
trigger-table work, not taste work. Three durable moves, from the reaction-term
attack ([[reaction-term-naming-species-split]]):

- **Check the refinement invariant BEFORE looking for a family word.** A uniform
  word is *forbidden* when a theorem splits the siblings. Reactions: locality
  *within the fiber* splits 1 multiplier (collision: continuum object is
  `Σ_t δ`, a DISTRIBUTION) from 3 kernels (Fredholm functions) — so "kernel" as
  the family word is false, not merely imprecise. The honest output is
  **species words on the leaves + a genus word on the ABC**, and the genus stays
  greppable by ONE token precisely because the leaves are not uniform. Bonus
  elegance check: make the layer-1 species BIJECT with the already-landed
  layer-3 species (multiplier↔`MultiplicationOperator`, kernel↔
  `IntegralKernelOperator`) — then the seam is *species-preserving*, a real
  property that is strictly weaker than the functor the seam does NOT have.
- **The discriminating test for a NAME is the invariant the name promises, run
  against the object that violates it.** A name is a claim; test the claim.
  Worked: `ReactionTerm`'s genus invariant is decomposability ⟹ assert
  `A(m⊙x) == m⊙A(x)` for a **cell-varying** mask (a constant mask cannot fail),
  where all four reactions PASS and **streaming FAILS**. Species: the multiplier
  additionally commutes with a **group-wise** mask, the kernels must FAIL that —
  and the `Σ_s` used must be genuinely off-diagonal, because a diagonal-only
  `Σ_s` passes the multiplier gate and a design that picks the species from the
  data's *accidental* diagonality is exactly the bug. Corollary: a name that
  cannot carry its own invariant (an invented weakest-true word) OWES a test
  that does.
- **A word already spent elsewhere in the codebase gets a delete-it-and-ask-
  what-breaks check, not an analogy argument.** `Law` (from `BoundaryTraceLaw`)
  for reactions: delete the BC law ⟹ **ill-posed** (it is a CLOSURE, enforced);
  delete scattering ⟹ a *different, still well-posed* problem (it is a
  GENERATOR TERM, applied). Two independent confirmations followed for free —
  the law is *affine* (carries `q`) while reaction terms are purely linear, and
  the domain has a live false friend (ENDF File 7 "thermal scattering law"
  `S(α,β)`). When reusing a precedent, split it: the **realizer** half
  transferred (Kalman realization is layer-agnostic), the **descriptor word**
  did not.

**New smell (candidate, 1st sighting): "the name states a contract the content
violates."** `MaterialXSField` is named as an apply-free datum (Dixmier field)
and carries NINE `apply_*` verbs consumed by two operator modules — the
fiber-operator ACTION living inside the datum. Distinct from Smell #16 shape 1
(two paths to one operator): here there is ONE path, hosted on the wrong LAYER.
TELL: a class whose docstring says "data/descriptor/field" while its method list
says `apply`/`add_`. FIX is relocation (it is a shared primitive, not a twin) —
and the *name* is usually the correct one, so resist the rename reflex.

How to apply: at any "what should this be called" dispatch, (1) hunt the
refinement theorem first — a forbidden family word is the highest-value finding
and it must be said PLAINLY; (2) report the math-faithful name AND the domain
(NE) name, with the routing rule *types get the faithful name, accessors/docs/
equations get the domain name* — a type name is read once per design decision, an
accessor on every line; (3) emit a per-name buys/costs line, and for any name
that would promise an unsupported operation, say which operation and cite the
measurement that refutes it (`Functor` ⟹ `reduce-discrete ≠ discretize-reduce`,
measured). Do NOT invent a name where none exists — say "no faithful name
exists", give the least-bad invented one, and flag it as invented.

---

## L-013 -- Before accepting "the machinery lacks X", check whether a PREDICATE already computes X and throws it away — a `bool` return is the commonest way a primitive stays missing

Recurring shape on "what is missing in the machinery" dispatches. A proposal says a
capability is absent; the truth is that a **verification predicate already computes
the exact object** and destroys it at the `return` statement, because the return type
was chosen as `bool`. The capability is not missing — its **witness** is.

Worked (#326, [[quadrature-symmetry-quotient-frames]]): `SubgroupOfO3.is_invariant`
was offered as "the checking face; the quotient is what's missing". But
`_orbit_closure` (`symmetry.py:904-954`) computes the matched partner index `j` for
every `(node i, group element M)` — i.e. **the permutation representation, which IS
the quotient's only hard input** — and returns `bool`. Two OTHER modules then
re-implement the same permutation independently (`_compute_sphere_reflection_partners`;
MoC's `_reflected_azi_index`). So the honest finding is not "add a quotient class"
(the proposal's scale) but "**return the witness**", after which the quotient is one
further verb (`consolidate`). That collapses the estimated work by an order of
magnitude and is a strictly better answer for the requester.

TELL (grep-able): a function named `is_*` / `check_*` / `verify_*` / `*_closure`
whose BODY builds an index map, a permutation, a matching, a partition, a
factorisation, or a certificate — and whose signature says `-> bool`. Same family as
Smell #16 shape 1 (the re-implementations downstream are the confirmation), but the
CAUSE is different and so is the fix: shape 1 says "collapse two paths"; this says
"one path exists and is throwing away its output — widen the return type first, and
the twin paths delete themselves."

Discriminating first test: assert the consumer's hand-rolled artefact is
`array_equal` (0 ULP, L-002) to the predicate's now-returned witness. A
re-implementation with a different tolerance or tie-break diverges exactly on the
degenerate elements (self-partners / fixed points), which is where every bug in that
family lives.

How to apply: at any "map what is missing in the machinery" brief, before enumerating
new types, grep the predicates in the neighbourhood and read their BODIES. Ask "does
this `bool` know something?" A found witness reframes the whole deliverable — and a
proposal that over-scoped the gap is a finding worth reporting plainly, because it is
good news.

---

## L-014 -- An UNSATISFIABLE predicate is a wrong-ARGUMENT diagnosis, not a wrong-predicate one: check each argument's KIND before redesigning the relation

Recurring shape on "derive the correct formulation of predicate P" dispatches. A
gate is found to reject everything (or accept everything) once a broken checker is
fixed, and the brief offers a menu of alternative RELATIONS. Nearly always the
relation is right and one of its ARGUMENTS is of the wrong KIND — a
cardinality/topology mismatch that no amount of re-shaping the relation repairs.

The diagnostic: **compare the CARDINALITY or TOPOLOGY the predicate demands against
what the object can carry.** A finite object cannot satisfy a containment against a
continuous group; a band-limited claim cannot constrain a non-band-limited unknown;
a static table cannot carry a configuration-dependent generator.

Worked (#326/Q2, [[quadrature-symmetry-quotient-frames]]): `G_geom ⊆ G_rule` was
found unsatisfiable because `GEOMETRY_GROUPS` supplies `SO(2)`. Two one-line
theorems settle it without touching the relation: **(A)** `Sym(Q)` of a finite node
set is FINITE (an orthogonal map fixing a spanning set is `id` ⟹ `Sym(Q) ↪ S_N`),
so a continuous `G` is unsatisfiable *by any discretisation*; **(D)** the
correctly-derived requirement — the DISCRETE residual `Γ = G/G⁰` acting on the
fiber — is ALWAYS finite (discrete subgroup of a compact group). The predicate was
never wrong; it was being handed the half of the symmetry group that the
dimensional reduction had already CONSUMED.

Two generalisable corollaries, both cheap to check:

- **A symmetry group that reduces a problem's dimension is SPENT; it cannot also be
  a requirement on the reduced problem.** Its continuous isotropy becomes the
  angular/fiber QUOTIENT (which domain the rule lives on), its free part becomes
  the CONNECTION/redistribution term, and only the DISCRETE residual is still owed
  as a constraint. Three parts, three fates, none discarded — one decomposition
  predicts "why does the slab use a 1-D rule and the cylinder a 2-D one" AND "why
  does the cylinder have an α term and the slab not" from the same split.
- **A "vacuous" candidate framing usually has a non-vacuous sibling one theorem
  away.** "Exactness space is `G`-invariant for every rule ⟹ the test says nothing"
  is true; the sibling is `E = Q∘(Id − Π_V)` (the error functional IS the
  quadrature on the aliased-out part), from which `G ⊆ Sym(Q)` ⟹ `E` annihilates
  every NON-trivial isotypic component at EVERY degree (average `E[f]` over `G`).
  Before discarding a framing as vacuous, apply the group-average to its error
  functional.

Discriminating first test for this family: a case where the candidates DIVERGE on a
parameter the current predicate cannot see (odd vs even `n_φ`), and where the
derived answer REPRODUCES an independently-established in-tree guard (ERR-042) as a
consequence. **A derived predicate that re-derives an existing hand-written guard is
strongly confirmed; one that contradicts it owes an explanation.**

How to apply: at any "the gate rejects/accepts everything, pick a new formulation"
brief, tabulate `(argument, kind it has, kind the relation needs)` FIRST. If a kind
mismatches, fix the argument and stop — the menu of alternative relations is a
distraction. Then hunt for an existing hand-written guard the corrected predicate
should reproduce; that reproduction is the cheapest available confirmation.
