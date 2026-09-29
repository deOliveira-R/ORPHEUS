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
against the already-landed layer below. (b) ⚠ AMENDED 2026-09-26 by the user's
collision ruling (`feedback_naming_collision_perfect_match.md`): a spent word is a
CONTEST won by its formal PERFECT MATCH; the other object cedes. Original: a word already SPENT in this repo on
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

**D11 — An INVARIANT measured under ONE Strategy is a Strategy artefact until a second
Strategy reproduces it.** The instrument doctrine's "a negative reading needs a positive
control" has a dual for POSITIVE readings of a "canonical"/"unique"/"converges to" claim:
the control is a SECOND, structurally different Strategy (a different splitting, a
different x₀, the limit of a neighbouring question). TELL: a claim of the form "the gauge /
representative / limit is canonical" backed by one iteration from one initial iterate.
`[M]` 2026-09-26: the at-pole gauge read "canonical" from the fission splitting at x₀ = 0
(shift memo P8) and was the splitting's own conserved functional `⟨ψ†, N x⟩`; the
pseudo-transient and the noise ω→0 limit land 5.6 away, a ψ-multiple. The general law to
derive FIRST: pair the iteration step with the left null vector and read what it conserves.
→ pointer: `equilibrium_carrier_laurent_point_kinetics.md`

**D12 — A proposed DECLARATION on an object whose CONSTRUCTION already fixes the fact is
a second source of truth; and a label that GRADES an operator the tree already sums is a
COORDINATE, not a label.** TELL: "each X declares the Ys it carries, and a resolver maps Y
to X" where every X is built by a constructor that reads named Ys. Refute by deriving the
map from the constructors and gating it `array_equal` against the hand-built objects; the
declaration then has nothing to say that the construction does not. Second half: when the
Ys index a direct-sum decomposition (`R = Σ_y R_y`), the question-space coordinates are the
Ys' cells, and one cell is ONE object read three ways (direction, `w = 1` rate, `w = ψ†`
sensitivity by Hellmann–Feynman) — check that before minting a "rate functional" and a
"direction" as two concepts. Rider: scalar coordinates and FIELD directions do not mix;
test a proposed direction against the span of the scalar cells before calling it "one more
scaling". `[M]` 2026-09-27 reaction channels: HF to 1e-10 per cell; a nuclide direction
0.87 off the cell span. → pointer: `reaction_channel_grid_frames.md`

**D14 — Before granting a VERB to a summand, ask whether the verb's semantics needs the
SUM.** A sampling verb needs the row MEASURE, which is normalised by the total over every
summand (`Σ_x/Σ_t`), and the free flight is a RESOLVENT of the sum, so `sample` on a term
is unspellable for half the members; the right home is a VIEW the whole object yields (the
pencil-from-question precedent). Second check, same dispatch: the verb's REPRESENTATION
requirement — a signed representation (a Legendre-truncated kernel) supports `apply` and
not `sample`; positivity is a property of the representation, so "two verbs on one datum"
is really "two verbs on two images of one kernel". TELL: a Protocol member that a subset
of implementors would stub; two normalisations of one measure in the shipped code.
`[M]` 2026-09-27 MC seam: majorant from stored `SigT`, collision from the partial sum.
→ pointer: `mc_seam_feynman_kac_particle_frames.md`

**D15 — Two arrows that differ by a diagonal factor (an adjoint and a section, `R^H = Σw·E`)
are ONE pullback read against two CODOMAIN MEASURES; compute `ρ = P_φ 1` (the pushforward
measure's density against the codomain measure) before typing them apart.** `r∘e = id` is the
theorem "ρ ≡ 1 on the image"; a dead slot is ρ = 0; a section fails exactly where ρ = 0 on a
live slot; a deck's measure-preservation (ERR-042) is `P_σ1 = 1`. TELL: an anti-swap two-type
design justified by "the two arrows differ by exactly Σw"; a "the Jacobian is not this type's
business" disclaimer (the Jacobian is `P_φ1`, and it lives on the BINDING — it is the
fibre-indicator frame's `discrete_gram`, computed and discarded). Second check: a point-map
pair proposed beside a Frame is the frame with a NODAL basis on the codomain and the map in
the pushforward measure's nodes — measure the six identities before minting a sibling type.
`[M]` 2026-09-27: with `μ_N := π_*μ_M` the section IS the pullback (array_equal). RIDER
2026-09-28: ρ is the FRAME-BOUND SPECTRUM of the nodal frame under the declared codomain
measure (`[M]` 4.4e-16 on 4 fixtures; A = min ρ, B = max ρ, Parseval ⟺ ρ ≡ 1), and the
canonical dual's projector `Φ G⁻¹ Φᵀ W_M` is ν-FREE — so "does frame theory survive a
coefficient measure" is answered by a two-line identity before any probe. → pointers:
`pullback_pair_frobenius_perron_frames.md`, `frame_hierarchy_kernel_coupling_frames.md`

**D16 — A verb claimed to be a PROJECTOR is tested on an element of its own span and for
idempotency; a POSITIVE partition-of-unity table under a row-sum normaliser is a MARKOV
operator (a coupling's disintegration), not a projection, and the two split exactly when
the table stops being one-hot.** A "measure on the coefficient side" has THREE canonical
fillers — declared, the Gram (L² pushforward, orthogonal projection), the marginal `M·1`
(L¹ pushforward, conservative disintegration) — and a tag that dispatches between the last
two (`GramStructure`) is a missing SLOT, not a missing type. TELL: a docstring saying "the
row sum IS the diagonal" or "R·1 = 1 collapses the probe" one sentence before claiming a
projector; a verb named `project` whose callers all want conservation. `[M]` 2026-09-28:
`f = T·[2,6]` in span, `conjugate(gram_inverse) f = [8/3, 4, 16/3]`; idempotency defect
0.11–0.124 on PoU, 0 nested; two-step condensation off by 0.55 on a straddling mid grid.
Rider: on a finite manifold a stored table IS a function family, so `Basis` is not loose
there — check the INHERITANCE EDGE instead (the PoU child was declared a subclass of its
one-hot special case). RIDER 2 (2026-09-28, same day, SELF-CORRECTION): **slot-vs-type is
decided by VERB VALIDITY across the fillers** — a slot is right when every verb is defined on
every filler; a verb undefined on one filler (`(M·1)⁻¹M` on a signed table: `M·1 = [2,0,0]`)
and another unwanted on the other (`R(MR)⁻¹M` on PoU: the dense solve, 0.45 off conservation)
means two TYPES under a shared abstract parent, and the morning's "one slot, three fillers"
was wrong for the Markov arm. RIDER 3: **"does not compose" is tested by the COMPOSED-KERNEL
route, never by two-step-vs-one-step alone** — Markov kernels always compose
(Chapman–Kolmogorov); two-step == one-step-through-`T₁T₂` to 2.5e-16 while both differ from
the direct table by 0.1–0.36; the defect is the table FACTORY (not a functor from the grid
poset: it re-declares the within-group model per level), and the decidable exactness
condition is "the second kernel is one-hot" (the tower property needs a nested filtration).
RIDER 4 (⛔ AMENDED 2026-09-28 evening, SELF-CORRECTION): the morning's "a positive TEST WEIGHT
and the same weight as the MEASURE are one object (7 of 7, 1.6e-16)" was TWO mistakes in one
sentence — (i) the 1.6e-16 is the FP re-association of ONE einsum (`PG(w).analysis(f) ==
Gal(dV).analysis(w f)` is `array_equal`, the same call), so it measured one verb, not two objects;
the objects differ on every OTHER surface (domain metric, coefficient metric, `.H` — the flux
shape is in the OPERATOR or in the METRIC, and an adjoint reads only the operator: CMFD
prolongation 2e-16 vs 0.66; identity flux-dependent under axis-built spaces; Gram unshareable);
(ii) "7 of 7" counted frame CONSTRUCTIONS when the population is collapse MORPHISMS (15; 5 are
ratios of two pushforwards with DIFFERENT densities that no single-weight frame spells). The
rule that survives: before accepting "X is a test weight, never a measure" ask which weight
could not be a measure (a SIGNED one, 0 sites) — AND before accepting "X is one object in two
spellings", test the surfaces the equality did NOT touch and count the population the claim is
about, not the population that was easy to AST. → pointers: `weight_placement_weld_attack.md`,
`frame_hierarchy_kernel_coupling_frames.md`

**D17 — A "two spellings agree to 1e-16" claim is a TAUTOLOGY until the surfaces the equality
did not touch are listed and tested; and a "state into structure" weld is detected by four
probes, not by an argument.** TELL: an `[M]` equality on one verb (`project`, `apply`) offered as
evidence that two OBJECTS are one; a census by `ast.Call(TypeName)` when the question is about
morphisms (the sites the type cannot spell were hand-rolled and are invisible to it). PROBES, in
order: (1) is the equality one contraction in two einsum orders (write the hand einsum; if it
reproduces the same 1e-16, it is); (2) `==` and `hash` of the induced spaces across two states
under BOTH identity regimes the tree has (nominal: the metric-blind seam is crossed silently;
structural: the space changes per state); (3) `.H` of the map in each spelling and a named
consumer of the adjoint (a prolongation) — the shape lives in the operator or in the metric,
never both; (4) the count of objects the spelling needs per trailing axis (ng, ng²) and what it
forecloses (signed, complex). A factory that computes closed-form overlap integrals is the
shadow of a composite `analysis_c ∘ section_f` over one structural measure — derive the
composite and it predicts the factory's own chain defect (`M_c(P_m − I)S_f`). `[M]` 2026-09-28
weight attack: A2 array_equal; B2 True/0.8 and 3-of-3 distinct; B3 2e-16 vs 0.66; G2 3.9e-16.
→ pointer: `weight_placement_weld_attack.md`

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

**D13 — "ψ†" names TWO vectors in a k-solve tree, and only one pairs correctly with
a pole's derivative.** The left null vector of the BALANCE at the pole
(`balance(pole)^† ψ† = 0`, the tree's `EigenPosing.H()`) and the left Perron vector of
the ITERATION operator `A⁻¹F` differ by `A^{-†}`; Jacobi's formula / the implicit
function theorem need the former. `[M]` 2026-09-27 (diffusion slab, p3c): HF-vs-FD
along a removal cell agrees to 4e-7 with the balance null vector and is 21 % off with
the iteration-operator vector — on BOTH the α and the k pole. Same class as the
metric-vs-Euclidean 4 % trap. Before pairing anything with a ψ†, gate
`‖balance(pole)ᵀψ†‖ ≈ 0` directly; and a "Λ_eff = Λ + …" check that divides an
adjoint-weighted derivative by an adjoint-weighted rate is an IDENTITY for any weight —
its only content is the derivative itself. → pointer:
`posing_ontology_clean_attack_frames.md`

**D18 — On a RE-CONSOLIDATED design text, count each named datum's HOMES before checking any
identity; and once an identity is confirmed `[M]`, attack its attached QUANTIFIER with the degenerate
fixture.** A consolidation merges rulings made on different days, and the failure mode is not a wrong
equation but one datum assigned two homes in two paragraphs (a constructor home vs a derived-view
home; a layer-2 value vs a layer-3 plan; a word retired vs assigned). Instrument: list the named data,
grep the text for each, count homes; ≥2 is a finding before mathematics. `[M]` 2026-09-28 clean attack 2:
6 of 6 checked identities held (1e-16) while 4 of ~30 data had two homes (`derived_SigT`, the material
field, the implicit/explicit labelling, `project`), and 4 universals attached to correct identities
failed ("iff" on the two-step law — exact with a fractional second kernel when mid = fine; "the
canonical dual does not depend on ν" — the FRAME does, the projector does not; "every spatial kernel";
"bit-identically" for a coarse-denominator ratio). The quantifier test: the degenerate member of the
family (mid = fine, ν = counting vs pushforward, the harmonic leaf's `analysis(1) = [2,0,0]`) decides an
"iff"/"every"/"total" in one line. TELL: a `[M]` tolerance quoted beside an "iff"/"every"/"never" — the
tolerance measured the identity, not the quantifier. → pointer:
`posing_ontology_clean_attack_frames.md` (second attack section). RIDER (2026-09-28, coproduct
attack): a TABLE and a WORKED EXAMPLE in the same memo are two homes of one datum — diff them before
trusting either. `[M]` the trace-inflow memo's table said `L_tt = diag(+I, −I)` while its own worked
row carried `− ψ_in(t2)`; the block has 16 of 32 off-diagonal entries.

**D20 — A control that stays green is a THEOREM or a DEFECT, and the automorphism group of the fixture
decides which.** Before calling a law "blind" to a relabelling / reordering / permutation, compute
`Aut` of the object the control permutes: if the permuted members are isomorphic, every
isomorphism-invariant law (round trips, metric, zeros, adjoint pairing) is blind BY THEOREM and the only
separating datum is the label — so the label must enter the identity; if they are not isomorphic, the
blindness is a gate defect (the gate did not read the datum that separates them). Then count the
shipped fixtures with `Aut ≠ 1` (a positive control for "labels matter") and name the first consumer
where it is non-trivial. `[M]` 2026-09-28: 0 of 3 coproduct fixtures isomorphic; on a synthetic
isomorphic pair 4 of 4 gates blind; the first real consumer is the 6-group precursor summand family.
TELL: a memo reporting "the law is invariant under relabelling" as a blindness without saying whether
the relabelled parts are isomorphic. → pointer: `coproduct_direct_sum_two_level_frames.md`

**D19 — An object defined "up to X" has no X to ask for: before minting a type or a field to
carry a quotient datum (a scale, a sign, a phase, a gauge), WRITE the X-invariant the object
actually is and test its invariance under the group; the datum then belongs to a REPRESENTATIVE
(a section), and a "kind" claimed for its admissibility is measured on a fixture WITHOUT symmetry
and one WITH.** TELL: "the fundamental mode happens to have the appropriate gauge"; "only the X
admits a Y-gauge"; a proposal to ask an equivalence class for a datum that the class quotients
out. Instrument: the Laurent residue / Riesz projector `ψ⊗ψ†/⟨ψ†,T_d ψ⟩` under random `(a, b)`
rescalings, and the contour integral as the eigenvector-free control; then the admissibility
predicate counted over ALL poles on an asymmetric fixture (a unary predicate on a mode is usually
binary on (mode, functional) and ternary with the symmetry group). `[M]` 2026-09-28: invariance
1.6e-16, contour = formula 1e-14; the positive section existed on 20 of 20 and 240 of 240 modes on
the asymmetric body and on exactly the reflection-EVEN ones on the symmetric body; the DD
fundamental left the cone (14 of 240) with its section intact. The frame that DOES partition
modes there is the symmetry quotient — decisive for another question, refuted for the type.
→ pointer: `fundamental_mode_gauge_attack.md`

**D21 — A NESTED iteration's inner strategy is tested with a FINITE inner count from the nest's own
start before its `ρ` is trusted; a cold-started inexact inner converges to the WRONG fixed point and the
increment stop accepts it.** `x ← B_m(rhs)` from zero has fixed point `(B_m⁻¹ + A − D) x = q`; the
correction form `x ← x + B_m(rhs − D x)` from the previous outer iterate is the two-stage iteration
(right fixed point, `ρ` by the identity probe). TELL: an `inner=(piece, rhs)` callback with no slot for
the piece's previous state; a nest memo whose numbers were all taken at an exact inner. Instrument:
run m = 1 cold and read the final error against a dense solve (`[M]` 2026-09-29: relative error 1.0
accepted after 7 passes; tol-based cold gives error = the tolerance). Not covered by D11 (one Strategy)
or D17 (tautology): the defect is a contract with a missing state slot. → pointer:
`ordering_cut_interface_augmentation_frames.md`

**D22 — Before reporting `ρ`, a CW bracket, a Krylov count or a singularity guard on a stationary
iteration, RESTRICT to the cut's INTERFACE (the range of `N`, `rank N ≪ n`): the nonzero spectrum of
`M⁻¹N` lives there, `K ≥ 0` is the predicate the bracket needs (regularity of `M⁻¹` is sufficient, not
necessary), Krylov terminates in ≤ `rank N` + 1 steps there, and `1 ∈ spec K` is the singularity.**
TELL: `error_propagation()` returning an `n × n` dense matrix; "regular" asked of `M⁻¹`; a
feedback-arc set weighted by coupling strength offered as the cut chooser (equal weights, different `ρ`
— the cost is `A⁻¹`-non-local). Rider: an augmentation that copies the whole space (`expose` as a copy
of `V`) inflates the interface from `rank T` to `n` — expose the RANGE. `[M]` 2026-09-29: 8 cuts ≤
3e-9; CW valid on two non-regular fixtures; GMRES 17 vs 255. → pointer:
`ordering_cut_interface_augmentation_frames.md`

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
| Collapse-morphism-blind — treating a marginalization as a weight=1 average (→ L-010) | 2 (XS coarsening; 2026-09-27 MC tally binning: the raw collision tally is the MARGINAL `M` of the empirical measure, the per-lethargy flux `tally/du` is the AVERAGE `G⁻¹M` — `mc/solver.py:446, :637`) | PROPOSED for Part C at this second sighting (MC tally, a different problem class); the row's own prediction named it |
| Vanishing-flux-function endpoint / metric-invisible-yet-active DOF (→ L-015b) | 2, both curvilinear SN; a RELATED third of a different kind 2026-09-27: the composite's trace DOFs carry no time term (`rank T = 40 of 48`), so the time family on `bulk ⊕ trace` is a DAE pencil and a "generator on V" does not exist | a non-transport degenerate-drift endpoint (Fokker–Planck, Sturm–Liouville, population balance) |
| The name states a contract the content violates — a class documented as "data/descriptor/field" whose method list says `apply_*`; distinct from Smell #16 shape 1 (one path, wrong LAYER); FIX is relocation, and the name is usually right (→ L-012) | 1 | a second, non-XS-field host |
| Identity-scarcity accretion — a container wins every placement because no candidate owner's `__eq__` separates the inputs (→ L-019 / M1.8) | 2 (one inverted: everything induced has structural `__eq__`, the hub has none) | a third outside the SN container family |
| A precondition spelled as a 30-line docstring caveat on a 3-line body wants to be a TYPE — declare the structure (DIAGONAL / POU / DENSE) and RAISE on the unhandled case (→ L-010) | 1 | a second non-Gram precondition |
| A property stated of a PHYSICAL PARAMETER that belongs to (parameter × discretisation × CHART × reduction) — TELL: "the outer extent is non-affine" / "DD vs LD splits on the sweep strategy" with no scheme named, or a second-difference test with no chart named; FIX: the second-difference test per scheme IN THE PARAMETER'S NATURAL CHART (affine ⟺ zero), then name the SCHEME trait (→ `d5_trait_and_mms_frames.md`; `shift_ontology_taxonomy_frames.md`; `posing_ontology_clean_attack_frames.md`) | 3 (transverse-coupling order; affinity in a width; 2026-09-27 the diffusion P1 face closure — a NON-SN scheme — makes the slab dilation rational while SN-DD is an exact pencil in 1/λ); a 4th of a DIFFERENT class 2026-09-28: "a mode admits a positive-functional gauge" stated of the MODE when it belongs to (mode × functional × symmetry group) — `[M]` parity, not kind (D19) | PROPOSED for Part C (third sighting met, non-SN scheme); the chart rider is the 2026-09-27 addition; the 2026-09-28 sighting widens the smell from "physical parameter" to "any unary predicate that is a product" |
| Multiplier-algebra embedding `f ↦ M_f` — scalar "coefficients" on terms that are really FIELDS per region/group; the scalars are the constant sections; unshielded composition directions are vectors in it (→ `coefficient_field_promotion_frames.md`; `posing_ontology_clean_attack_frames.md`) | 2 (XS-field promotion; 2026-09-27 the question space's cell coordinates) | PROPOSED for Part A.3 at this second sighting (a different problem class: question posing, not carrier typing) |
| Projector-vs-Markov split on a positive PoU table (D16) — a row-sum normaliser on a non-one-hot coupling yields a conservative disintegration, not a projection; TELL: `project` whose callers want conservation, a Gram-structure tag dispatch; FIX: the coefficient-measure SLOT with its three fillers (→ `frame_hierarchy_kernel_coupling_frames.md`) | 1 (2026-09-28 energy condensation; the CP `_normalize_rcp` and the MC per-lethargy tally are the SAME verb but were named from D4, not sighted independently) | a second sighting on a non-XS coupling (a non-nested SPATIAL remap, a response-matrix normalisation) with the idempotency and composition tests |
| Interface-blind iteration analysis — `ρ`/CW/Krylov/guard computed on the dense `M⁻¹N` when they live on `range N` (D22) | 1 (2026-09-29, SN cuts; the coproduct attack's "response form" was an [R] sighting of the same object, not an independent payoff) | a second sighting on a non-SN cut (a CP group loop, a coupled-physics lag) with the `K ≥ 0` test |
| Section-of-a-retraction augmentation — `expose` (term image) and a Galerkin-correction summand are two SECTIONS of `reduce`; TELL: "X is REDUCE's inverse" measured one way; a coarse correction called "outside the ordering" (D21/D22 pointer file) | 1 (2026-09-29) | a second augmentation kind (an overlap duplicate, a Lagrange multiplier) spelled as a section with a cut |
| Smell #17 (a) third form — a verb listed as "new" in a memo's NEEDS that is a COMPOSITION of two shipped verbs (`pushforward(φ).consolidate()`; `basis_space.apply_metric∘analysis`) | 2 (2026-09-28, both in the frame hierarchy attack) | already promoted as #17; this row only records that a memo's own NEEDS list is a place to run the #17 check |

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
- **A bilinear pairing / a correspondence of measured spaces / a framed-Hilbert category as
  an APEX above "kernel between two measure spaces"** — refuted 2026-09-28 FOR "what is the
  general object over frames": each is the same datum as a `LinearOperator` between two
  metric-carrying ends (the pairing is `⟨Kf,g⟩_ν`; the measured correspondence is the
  positive SUB-case, the coupling; the framed morphism is `K` in two frames). Optimal
  transport refuted for the same question (no cost functional on a fixed overlap); would
  fire on "choose the coarse grid to match two spectra".
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
