# Cross-Domain Attacker — Lessons (digest)

Behavioral lessons only: *what detection mistake did I make, and what correction
changed how I work?* Three bodies of content are deliberately NOT here:

- the frame-trigger CATALOG → `cross-domain-frames` skill, `reference.md`
  Part A / B / C (preloaded every dispatch);
- Smell #16's four shapes and the transport-resolvent backbone → `AGENT.md`
  "Promoted Library Kernel" (preloaded);
- durable frame-matches that became architecture → `MEMORY.md` §3 design
  pointers (cold files, opened on demand).

War stories — `file:line` measurements, campaign ids, worked derivations — live
in `lessons_archive.md`, which is the pre-2026-09-21 digest verbatim under its
original `L-001 … L-024` headings. Every entry below cites its old id, so
`→ L-013` is an address in that file.

THE SPINE, and every entry is one face of it: **an attack's value is a concrete
reformulation with a fail-able first test, OR a crisp refutation carrying its
structural reason and the question it was refuted FOR.** A frame named without a
payoff degrades the table's signal and is rejected output.

---

## Part 1 — Meta-lesson M1 (M2, M3, M4 retired 2026-09-22 into the definition's standing bars). Fire these before any trigger lookup.

### M1 — Measure the brief's PREMISE before reaching for a frame

The brief's framing is the likeliest thing to be wrong, and the repair is a grep
or a one-paragraph theorem, never an argument. Take the applicable measurement
FIRST, publish the count, and let it generate the answer. Consumer and
membership counts are an AST pass, never a line grep (`instrument-doctrine`).
The general form is `process-discipline` "Measure a brief's premise" since
2026-09-21; the instruments below are this agent's.

Items 1, 2, 3, 5, 6 and 10 retired 2026-09-22: the definition's step 0 and the path-scoped `plan-authoring` carry them. The instruments kept below are this agent's own.

4. **"Which of these N families needs X?"** N independent arguments is N chances
   to be plausibly wrong. Find the GENERAL-CASE derivation first (grep the
   literature folder for "general geometry" / "arbitrary" / "general
   formulation"); one equation decides every row and usually DERIVES rulings the
   project had asserted. Then extract the discriminator clauses that GENERATE the
   table, so every "no" row fails a nameable clause instead of carrying a bespoke
   excuse. And when the object under question is a general one composed with
   discretisation choices, do not answer "is it general?" — answer with the LAYER
   TABLE (continuous object / evaluation points / reconstruction rule) plus TWO
   tests pointing OPPOSITE ways: hold the measure fixed and change the chart (does
   it move? ⟹ not a measure invariant), hold the chart and refine the measure (has
   a continuum limit as a function of the coordinate? ⟹ it IS a chart object).
   Two opposed tests is what makes it a verdict rather than a preference; neither
   alone suffices. → L-021(a), L-021(c)
7. **"The gate accepts/rejects everything — pick a new formulation."** Tabulate
   `(argument, kind it has, kind the relation needs)` first. A
   cardinality/topology mismatch is never repaired by re-shaping the relation: a
   finite object cannot satisfy a containment against a continuous group. Then
   hunt the hand-written in-tree guard the corrected predicate should REPRODUCE —
   reproduction is the cheapest confirmation available; contradiction owes an
   explanation. → L-014
8. **"Should this object hold everything?" (God object).** Ask which object in
   the stack has usable IDENTITY before tabulating contents. A derived object
   lands on a container because the natural owner's `__eq__` cannot separate the
   inputs it depends on — the growth is a FIXED POINT, not a discipline failure,
   and relocation will not hold while identity stays weak. Prescribe the ORDER
   (strengthen the owner's identity, then relocate); the sequencing is worth more
   than the smell. TELL: a docstring that ARGUES for the container by
   disqualifying alternatives on identity grounds. → L-019
9. **"What is Problem vs Solution / model vs run / static vs traced?"** Write the
   TERMINAL OBJECT with its PARAMETER DOMAIN before proposing any rule — a rule
   of the form "Problem-side iff changing it moves the answer" is defeated by a
   one-parameter family. The test that a candidate rule is not still classifying
   by mechanism: does it collapse two kinds the domain expert names as different
   onto ONE object? Then enumerate the domain's degenerate cells and check each
   against the tree (a unification that ABSORBS a shipped object deletes content;
   one that only adds cells is speculation), and find the admissibility guard that
   already exists before minting one. → L-024(a)
## Part 2 — Standing detection rules (not instances of Part 1)

**D1 — Property-vs-TYPE is decidable, not a taste call.** A representation earns
a distinct TYPE iff there are **≥2 bases that are NOT canonically isomorphic**
(the iso depends on a quadrature/node choice), connected by a change-of-basis
operator that is itself MODELED and APPLIED (carries truncation error, has an
adjoint, participates in the algebra). All three clauses, countable by grep; zero
applied non-identity morphisms ⇒ PROPERTY. Corollaries: defer with an EXPLICIT
trigger (name the latent consumer whose arrival flips the verdict, since "no
current consumer" is not "never"); and an axis that changes the ARITHMETIC
INTERFACE or the SHAPE cannot be a phantom type parameter (erased at runtime,
one shared `__add__` body) — it must be a distinct class, which makes the
orthogonal-factor `Leaf(RoleMixin, RepBase)` the only elegant form and puts the
parametrization on the OPERATOR contract `[Din,Cout]`, not on the carrier.
→ L-004 · pointers: `spatial_order_type_vs_property_criterion.md`,
`rep_role_grid_double_category_frames.md`

**D2 — Deploying the resolvent backbone** (the spine itself is `AGENT.md`'s
kernel; this is how to use it). (a) When a "find-the-special-value" family
(k, α, time step, fixed source) looks like distinct solvers, check whether they
are POSINGS of one generalized eigenproblem `Aψ = λMψ`; if so the only genuinely
per-method layer is the loss-operator REALIZATION. (b) Generality flows TOWARD
the opaque resolvent interface: between two iterative drivers sharing a loop
body, the one exposing its resolvent behind a Protocol is the general engine and
the concrete-`(L,S,F)` one adapts INTO it — a "retire the opaque loop" plan points
the deprecation arrow the wrong way. (c) The backbone tells you WHERE a foreign
frame fires: ask which MEMBER has the matching algebraic shape before pointing a
frame at "transport" (saddle-point / inf-sup / mixed-FEM fires on the
diffusion/low-order member only; the sweeps are characteristic-triangular and have
no saddle to stabilize; (c) is now in `AGENT.md`'s backbone kernel, 2026-09-22). → L-007 · pointer:
`dsa_saddle_point_mixed_fem_frames.md`

**D3 — A change-of-basis frame's OWNER and its Galerkin-vs-PG discipline are
predicted by the operator's SYMMETRY, not by which subsystem calls it first.** A
frame `(M,R)` is OWNED by the operator whose EIGENBASIS it is, and is GALERKIN
iff that eigenbasis is orthogonal (commutant membership ⇒ Schur block-diagonality;
self-adjoint zonal kernel ⇒ Funk–Hecke ⇒ `M* = R` up to the Plancherel metric).
Three outcomes, all sighted: one operator's eigenbasis ⇒ owned + Galerkin; no
symmetry ⇒ solution-weighted PETROV-GALERKIN owned by nobody; the group in
EVERY equivariant operator's commutant (a symmetry quotient) ⇒ owned by the
PROBLEM, and still Galerkin, on a smaller space. Detection rule before typing an
`R∘A∘M` as PG: *is the test≠trial gap a SOLUTION weight or a GROUP
identification?* A group identification is Galerkin-on-a-sub-block. The genuine
falsifier of "subsystem X owns the frame" is a second consumer whose TRUNCATION
ORDER is set independently of X's operator. → L-009 · pointers:
`harmonic_frame_ownership_funk_hecke.md`, `dsa_rp_angular_frame.md`,
`quadrature_symmetry_quotient_frames.md`

**D4 — A conserved-quantity COLLAPSE splits by WHAT is conserved, which fixes
the MORPHISM — never by a weight.** Do not accept a "same projection with vs
without a weight" framing: name the preserved functional per channel. A reaction
RATE ⇒ AVERAGE (`G⁻¹M`, normalized); a PROBABILITY or MASS (`Σχ=1`, a particle
count) ⇒ MARGINALIZE (`M` alone). They differ by `G⁻¹`, so a "weight=1
degenerate of project" would divide by the bin count and break `Σχ=1` —
provably wrong, not merely inelegant. Confirm both collapses act on the SAME
axis first; if not, the "asymmetry" is two reductions wearing one channel name.
The discriminating test is order-non-commutativity on a multi-axis channel.
→ L-010 · pointer: `xs_coarsening_collapse_marginalize_vs_average.md`

**D5 — A "coupled / nested block system" over a carrier that already has a
biproduct is a FREE RE-ASSOCIATION, not a new object; and defer-until-≥2 counts
KINDS, not instances of the word.** `Mat₂(Mat₂(𝒞)) ≅ Mat₄(𝒞)`, and the G-adjoint
reads `G⁻¹AᵀG` at any partition granularity, so name the off-diagonals and lift
the block-index FREEZE instead of adding a layer. The discriminating challenge
has a definite answer: "exhibit a LINEAR coupled system expressible nested but
not flat" — every candidate is flat-re-expressible (a view) or nonlinear (not a
LinearOperator at all). For the build-now count, tabulate each cited instance's
(off-diagonal structure, metric definiteness, solve kind, linear?) and build only
where ≥2 rows MATCH; the over-reach dual is an assumption baked from the first
kind that excludes the others. → L-011 · pointer:
`coupled_system_field_bc_frames.md`

**D6 — A NAMING dispatch is frame-detection work: hunt the refinement theorem,
then the spent word.** (→ now an `AGENT.md` task-type row, 2026-09-22; the procedure stays here.) (a) Check the refinement invariant BEFORE looking for a
family word — a theorem can FORBID a uniform word (locality within the fiber
splits one multiplier from three kernels), and then the honest output is species
words on the leaves plus a genus word on the ABC, with the leaves biject-checked
against the already-landed layer below. (b) A word already SPENT in this repo on
an ORTHOGONAL AXIS of the same object is a kill, not a cost, and it is one grep
of the stem — read what it already MEANS; check the DIRECTION too (a word whose
established sense makes the candidate object its INPUT dies twice). A borrowed
word gets the delete-it-and-ask-what-breaks check (a CLOSURE deleted ⇒ ill-posed;
a GENERATOR term deleted ⇒ a different well-posed problem); when reusing a
precedent, the realizer half can transfer while the descriptor word does not.
(c) For a CONTAINER prefer the ROLE name over the CONTENTS name — a contents name
is falsified by every content move a live campaign has chartered, and a contents
name patched with an adjective is the tell. (d) `XOptions` / `XConfiguration` /
`XSetup` reading *right* is a DIAGNOSTIC that the object carries no invariant:
put the bag names in the tournament table, labelled as the diagnostic. (e) Settle
the ONTOLOGY first and let it pick the name; report the math-faithful name AND
the domain name, routing *types get the faithful name, accessors and docs get the
domain name*. If no faithful name exists, say so, give the least-bad invented one,
and flag it as invented. → L-012 · pointers:
`reaction_term_naming_species_split.md`, `container_ownership_dof_criterion.md`

**D7 — Before proposing a foreign basis, identify the NODES; before designing an
endpoint condition, evaluate the FLUX FUNCTION.** (a) A node set built from a
trigonometric / roots-of-unity / equispaced construction is very often a CLASSICAL
GAUSSIAN rule in a transformed variable — and then the matching orthogonal family
brings a free transform, a diagonal Gram, a known exactness degree, and a
truncation EXACT at the nodes. TELL: a `cos` of an equispaced grid described as
"equispaced" / "uniform" / "the trapezoid rule". CHECK: apply the geometric map
and compare against the Gauss nodes of `jacobi(α,β)` — a three-line probe against
the shipped generating measures. (b) A first-order operator whose FLUX FUNCTION
vanishes at a domain endpoint admits NO boundary condition there (`f(endpoint)=0`
kills `[fψφ]` for every ψ, φ; Fichera zero ⇒ characteristic, neither inflow nor
outflow), so any seed / starting value / admission gate there is a discretisation
artifact and what the endpoint really supplies is a COMPATIBILITY condition. Then
COUNT the discrete system: "one condition short while two endpoint data are
computed and one discarded" is a far sharper finding than "the seed is
inaccurate", and the discarded datum is a free a-posteriori estimator. → L-015 ·
pointer: `cylinder_angular_march_jacobi_ladder.md`

**D8 — Compute the END-TO-END gain before proposing any march-direction remedy.**
A measured interior amplification baits Miller / Gautschi ("march where the
unwanted homogeneous solution decays"). With end-to-end gain exactly 1 the
recurrence is REVERSIBLE: forward and backward seed errors have the identical
profile, there is no dominant/minimal PAIR, and reversal is a no-op for
conditioning. Unit gain is usually forced by an exact anti-symmetry of the
scheme, and then the payoff moves: gate the OPERATOR identity (`R A_h R + A_h =
0`), because the scalar shadow the tree already gates cannot separate variants
that violate it. → L-016

**D9 — A "fully probes" claim is about operator LINEARITY, not input polynomial
degree.** For a LINEAR operator an input merely NON-CONSTANT in the active
variable already probes the full map; higher degree only moves where in the
already-probed range you land. Spend the degree budget on QUADRATURE exactness
instead. Paired hazard when lifting a Cartesian reference to curvilinear: a
redistribution term carries `1/r`, so a slope driver must vanish at the origin
for pole-regularity — a slab-derived ansatz silently drops both the term and its
regularity constraint — and the geometry MEASURE enters the L2 error norm, so an
unweighted norm mis-measures the convergence order. → L-008

**D10 — An in-repo analogy can be ADJECTIVE-accurate and LAYER-wrong; a layer
error inverts the conclusion instead of degrading it.** TELL: a two-layer
precedent (data + binder, model + view, kernel + driver) cited by ONE of its
layer names. Map the new problem's objects onto BOTH layers and keep the mapping
that preserves ARITY (field count and verb count); a "thin data class" maps onto
the precedent's data class, never onto its binder. The grep-checkable
discriminator: **in a data/binder split the DATA object's verbs return ARRAYS and
only the BINDER returns OPERATORS** — returning arrays is what keeps the data
module's imports empty, and `datum.bind(space) -> LinearOperator` inverts exactly
that, which no dispatch mechanism repairs. Use the IMPORT DIRECTION as the
refutation; an `isinstance` chain is only its symptom. Two riders: binding is
BINARY so neither operand owns it (the answer is a third object, which is also
where caching lives — check a chartered "binding base" for one field too few);
and a 3-of-4 uniformity gap is INFORMATIVE, not a smell, when the fourth member
is the DEGENERATE case of the same construction (criterion: different
construction ⇒ smell; degenerate case ⇒ unifying DELETES content). → L-020 ·
pointer: `kernel_as_frame_layer_inversion.md` · `plan-authoring` PRECEDENT-LAYER
since 2026-09-21

---

## Part 3 — Smell-promotion ledger

The skill's Growth Protocol: a Part C smell earns its slot on an INDEPENDENT
sighting (a different problem class), and refutations stay in memory rather than
in the trigger table. Fire these inline meanwhile; do not re-promote what is
already promoted (Smell #15, Smell #16).

| Candidate smell | Sightings | What would count as independent |
| --- | --- | --- |
| Frame-leak naming — a model-agnostic slot named after ONE consumer's physics; TELL: "generic in X" beside a parameter named after a specific X₁; FIX: name the ROLE in the INTERSECTION of all consumers' domains; first test: a second consumer reading it with NO first consumer in scope (→ L-006) | 2, both naming cases | a third sighting carrying the 2nd-consumer-with-no-1st-in-scope test, distinct from a name |
| Eigenbasis-blind frame placement — operational-pipeline vocabulary ("natural data carrier of the Galerkin pipeline") where a Funk–Hecke frame is unnamed (→ L-009) | 3, all the same ANGULAR frame | a genuinely non-angular eigenbasis frame |
| Collapse-morphism-blind — treating a marginalization as a weight=1 average (→ L-010) | 1 (XS coarsening) | a non-XS conserved collapse (MC tally binning, flux→current) |
| Vanishing-flux-function endpoint / metric-invisible-yet-active DOF (→ L-015b) | 2, both curvilinear SN | a non-transport degenerate-drift endpoint (Fokker–Planck, Sturm–Liouville, population balance) |
| The name states a contract the content violates — a class documented as "data/descriptor/field" whose method list says `apply_*`; distinct from Smell #16 shape 1 (one path, wrong LAYER); FIX is relocation, and the name is usually right (→ L-012) | 1 | a second, non-XS-field host |
| Identity-scarcity accretion — a container wins every placement because no candidate owner's `__eq__` separates the inputs (→ L-019 / M1.8) | 2 (one inverted: everything induced has structural `__eq__`, the hub has none) | a third outside the SN container family |
| A precondition spelled as a 30-line docstring caveat on a 3-line body wants to be a TYPE — declare the structure (DIAGONAL / POU / DENSE) and RAISE on the unhandled case (→ L-010) | 1 | a second non-Gram precondition |

## Part 4 — Refuted-frame ledger (high-prior frames that keep NOT firing)

`AGENT.md` carries the reasons in compressed form; this ledger carries the frame
NAMES and the trigger each actually needs, so an UNEXPLORED line can be written
without re-deriving the refutation. → L-001

- **Wiener–Hopf factorization** — wrong solver FAMILY. Native to the
  Chandrasekhar/H-function half-space line, structurally incompatible with a
  bouncing-Peierls or sweep formulation; keeping the two families independent is
  itself a V&V asset (independent references).
- **Homology / chain complex** — tempting via the word "boundary", but
  `∂∘∂ ≠ 0` in transport (two reflections compose to a non-trivial map; the
  boundary trace and its extension are a dagger adjoint PAIR, not a
  differential). No `∂²=0` ⇒ no homology payoff.
- **Tensor networks / MPO** — needs a genuine bond-dimension trigger: a rank-N
  chain where N is a real truncation knob. A biproduct or a 2-surface BIE is
  bond-dimension-1/2 DEGENERATE. Do not promote until N ≥ 3 actually ships.
- **Differential geometry / Christoffel** — needs a CURVATURE term to
  redistribute. Straight Euclidean chords and Cartesian cells have none; it fires
  for `(1−µ²)/r ∂_µ`, not for geometry-of-the-domain questions.
- **Category theory / operad / PROP** — the skill's Part A.2 row already marks it
  low-signal; the concrete win it gestures at is normally already captured by a
  nameable frame (biproduct, forgetful functor with explicit laws). Name the
  concrete frame; list category theory UNEXPLORED unless a specific functor or law
  produces a test.
- **Affine geometry / torsors** — the frame is sound, its ORPHEUS example
  overturned; run the two-question test the skill's A.1 row carries (canonical
  zero? physical superposition? two yeses ⇒ vector space + cone predicate) before
  applying it.
- **Rayleigh–Ritz on a REDUCED variable** — the frame is sound and its smell is
  promoted (Part C #15), but the OBVIOUS reformulation target is the wrong one.
  Refuted FOR "is the F.4 white-BC closure secretly rank-1 Ritz on the boundary
  trace?" (RH13, 2026-04-22): `R[1] ≈ 2e-4` against `k_eff = 1.4963`, because
  the CP eigenproblem is SCHUR-REDUCED and the eliminated volume block is
  load-bearing. The FACT it establishes, and the standing rule: **a Schur
  complement is not self-adjoint in the reduced variable's inner product, so no
  Rayleigh quotient exists on the reduced variable alone** — before claiming a
  variational principle, ask which block was eliminated. The frame needs the
  FULL (volume + trace) eigenproblem with a nested ladder
  `V_n = V_{n-1} ⊕ µ·V_{n-1}`, Galerkin on the whole operator; the smell itself
  survived and was doubly confirmed on both kernels by the anisotropic-BC scan
  (`.claude/plans/archive/rank-n-closure-research-log.md`, RH13). Literature:
  Courant & Hilbert 1953 Vol. I §VI (min-max); Case & Zweifel 1967 §6 and
  Wendroff 1961 (Boltzmann variational theory). Salvaged 2026-09-22 from the
  retired memo `elegance_smell_rank_non_monotone.md`.
