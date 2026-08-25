# The posing filtration — problem posing is a monotone refinement, and every guard sits at its earliest decidable point

**Status: RATIFIED 2026-08-25** (user proposal, 2026-08-24 → three adversarial
rounds → surviving form ratified with its amendments). This file is the
**charter AND the derivation record** for the architecture that governs the
operator/mesh un-weld arc and re-shapes the landing surface of CS1.5′ and CS2.
It supersedes the *shape* of the CS1.5 Medium charter (campaign plan §2.5 and
`cs15_medium_unweld_design.md` — no `Medium` class exists in the ratified
form) while preserving that charter's surviving physics objectives (§8 maps
them). The **step decomposition for the arc is deliberately NOT here** — steps
are the next design round's output, designed against this charter, with the
user steering (surgical posture stands).

**Why this file is long.** It records the *reasoning*, not only the
conclusions. The archivist will later write the Sphinx theory page from it,
and a fresh session must be able to re-derive every guard placement from the
principle rather than trust a table. Per the user's instruction: for every
discharge point we state why that place is RIGHT — the positive reason is
near-singular — and we do not enumerate why other places are wrong (there are
many wrong places for many reasons; there is usually one right one for one
reason).

---

## 1. The concept this architecture serves

Everything in the problem layer is a map on spaces (user, 2026-08-24):

- a **field** is a map *space → values* — cross-sections map space to a
  cross-section value, flux maps space to a flux value;
- an **operator** is a map *space → space*;
- therefore the organizing object of the codebase is the **space with labeled
  axes**, and nothing fundamental requires a mesh object.

`SNMesh`'s diagnosis under this lens: a god object produced by **identity
scarcity** (things accreted to it because it was the only object with strong
identity while `FunctionSpace.__eq__` was weak — campaign plan §2.5's §8
mechanism, repaired by CS4b's axes-content identity) and by **naming
preceding machinery** (user: "we tried to give everything mathematical names
but their roles was more naming than machinery and a lot of things were hand
rolled"). The named objects are now *becoming* machinery and absorbing the
mesh's roles; `SNMesh` trends toward "a save state" — a pure data aggregate
with no machinery, which is a legitimate terminal role (organization +
persistence, #406) and NOT a god object, precisely because it would carry no
behavior. Whether that aggregate survives as a class is **an open user fork**
(ruled OPEN 2026-08-24); this charter resolves the *admissibility-placement*
half of the old "funnel" question (§5) and leaves the *container* half open.

**Scoping ruling (user, 2026-08-24): the arc builds the CONSUMED objects
first — space, fields, operators.** The consumers — solvers, strategies,
traversal (`dag_walk`, schedules, cumprod machinery) — are deferred to a later
arc whose consumption "will drastically change to have a lazy realization."
The filtration gives that later arc its criterion for free (§3, "lazy
realization").

## 2. The ontology: posing is a filtration

**The chain of problem-posing commitments is a filtration: a monotone
refinement of partitions of phase space. Each stage commits exactly one
refinement, and each commitment is only *statable* on the partition the
previous stage built.** That last clause is why the order is principled
rather than arbitrary: you cannot assign materials finer than regions exist,
cannot mesh what has no extent, cannot pick DOFs before a method says what a
DOF is. The linearization is not a convention; it is the dependency structure
of statability.

The stages, with the partition each commits:

| stage | declares / commits | partition refined to | measure that appears |
|---|---|---|---|
| **Materials** | `{id → Mixture}`, raw grid provenance per material | per-MATERIAL; energy per-GROUP (the library's partition of the energy axis) | the energy discrete measure's *data* (formal construction is method-time — §4.1) |
| **Geometry** (overlay) | region partition; deck identifications; boundary data; region→material assignment | per-REGION | none new (identifications *define* the quotient domain) |
| **Mesh** | cell partition conforming to regions | per-CELL | the spatial measure (volumes) |
| **State fields** | state data (T, ρ, burnup, …) on cells | none (data on the cell partition) | none |
| **Method head(s)** | the method + its discretization scheme | per-DOF — the **terminal refinement** of every axis | angular measure (if the method has one); spatial DOF measure (nodal/modal); energy possibly *refined* (unionization) |

Three structural theorems fall out of the ontology — each was independently
forced during the adversarial rounds and then recognized as a consequence,
which is the evidence the ontology is native rather than imposed:

**(T1) All axes formally construct at the method** (user's counter-9,
sharpened by the MC-unionization example — a Monte-Carlo head with a
unionized energy grid *refines* the energy axis, so even the energy axis is
not formally constructible before the method). The filtration derives this:
the method commits the last refinement of every axis, so only there is the
full axis datum present. Constructing an axis earlier would mint an object a
later commitment could invalidate — a mutable or twice-derived space, both
worse. Earlier stages **accumulate measures and data; axes resolve at the
head** (no half-axis object exists — `[M]` `Axis` requires `kind` at every
mint, `orpheus/numerics/axis.py:103-141`: "the basis character is physics and
must be spelled at every mint").

**(T2) The leak principle is the adaptedness axiom.** The tree's already-
ratified law — *"the mint consults exactly its defining data … a spectator
with `eg=None` must not flip the axis identity of a problem it does not
touch"* (`[M]` `orpheus/transport/mesh/material_mesh.py` `bulk_space`
docstring, ~:379-410) — is precisely the statement "every stage-k object is
measurable with respect to stage k's partition." The architecture's central
law was ratified before the architecture was named. Consequence used in §5:
declared-but-unassigned materials are **inert by construction** — no
spectator-warning machinery is needed, ever, because the method-time mint
reads reachable materials only.

**(T3) Coarsening is the same chain walked backward.** Condensation and
homogenization are conditional expectations onto coarser stages of the same
filtration. The landed machinery already implements the projections: the
collapse pair's `retraction` with its Parseval divisor IS the projection onto
a coarser partition (CS4b S6, `orpheus/numerics/frame.py` `_collapse_pair`),
and the condensation ruling "fractional-overlap downsampling" is the measure
pushforward between partitions. Posing walks forward; coarsening walks
backward; one structure. Corollary: **the filtration mints a space at every
stage** (each stage's partition with its measure is a space at that
resolution); the method's space is merely the terminal one, and the
inter-stage maps are the retraction/section family.

**Data spaces vs solution spaces — the precision T1 needs.** T1 governs the
*solution* (flux/DOF) space. **Data lives at the stage where it is declared**
(T2), on that stage's space: cross-sections are cell-wise data on the
cell×group space, constructible at the mesh stage with no method in sight —
which is why `MaterialMesh.bulk_space` legitimately exists today. A method
head whose solution space differs **re-poses** data onto its pose; the
measured precedent is S7's EE-1 (`[M]` the homogeneous solver re-poses XS
fields onto its pose — `orpheus/homogeneous/solver.py`, commit `2e054bfc` —
forced by G2.5's measure-authority mutation gate). Without this distinction
the charter would forbid an XS field before a method exists, which is both
false to the tree and wrong.

**Symmetry's exact role** (resolving the round-2/3 exchange): the group acts
on phase space; its orbit partition **bounds how coarse an admissible pose
may be**; among admissible refinements the good ones respect orbit structure.
So **refinement is the flow; symmetry is the admissibility bound and quality
criterion on refinements** — the user's monotone-symmetry-breaking argument
(each stage's object bounds the maximum theoretical symmetry; monotone
NON-INCREASING, not strictly decreasing) is the constraint surface of the
filtration, not its engine. This division also covers the case symmetry
alone cannot order: a generic heterogeneous problem's group goes trivial at
the overlay, yet mesh and method stages still refine meaningfully. The
per-stage symmetry-group *machinery* is aspirational (the ingredients ship:
`RigidMotion` in `orpheus/geometry/` — placed there for exactly this reason
(user) — `SubgroupOfO3`/invariance lattice #152, invariant tests #166, and
the constructive kernel measured in §5 guard 4); declared-vs-computed groups
is an open design choice, and the "mesher scores symmetry preservation"
criterion is the user's recorded aspiration. Until built, the symmetry story
justifies; it does not yet execute.

## 3. What each stage is, and what exists today

**Library vs Materials.** The library is all data you own; **`Materials` is a
declaration for THIS problem** — `{id → Mixture}` with each material's raw
grid provenance. The distinction is load-bearing (round-2 attack 9): it is
what makes stage-1 data lawful under T2. Naming ruled 2026-08-25: the word
"Medium" is **not** used for this object — "admission-previewed material
declaration with raw grid provenance" *is* a material list, and `Materials`
is its honest name. The word *medium* survives in exactly one place, where it
is physically exact: **`InfiniteMedium(mixture)`** — a single material
filling all space. Placement: the **problem layer**, never a method on
`Mixture` (round-2 attack 12, conceded: `mixture.as_infinite_medium` would
make `orpheus.data` import the problem layer — dependency inversion).

- Two rulings bind Materials' content (round-3, both fully accepted):
  **(R-mint)** *the mint is the law; admission is a preview.* The method-time
  mint re-reads the raw data (single source); the Materials-time coherence
  guard (the XD-4 three-outcome shape: agree / no grid / refuse) exists for
  fail-early ergonomics ONLY, and its correctness criterion is that its
  refusals are **implied by** the mint's — a preview may never disagree with
  the law, in either direction. **(R-raw)** *no early collapse.* Because a
  head may unionize, Materials carries each material's **full grid
  provenance**, never a pre-reconciled common grid — collapsing early is the
  lossy-return-type defect at the root, and it would throw away exactly what
  the MC head needs. Consequence: two heads over one Materials may
  legitimately realize **different** energy axes.
- The infinite path: `InfiniteMedium(mixture)` is a **complete problem**
  (the slowing-down problem — the energy sub-algebra posed alone; not
  transport at degenerate geometry). It runs the SAME generative primitives
  (axis mint, kernel binding) at trivial stage values; ergonomic entries are
  sugar over shared primitives, which is what dissolves the twin-path hazard
  (round-1 attack 1, withdrawn on exactly this condition — it is a
  condition, not a given: **zero downstream `if infinite:` arms; if one
  appears, the architecture has failed its own thesis**). It is also the
  in-solver evaluation tool (per-material spectra, k∞, condensation
  weighting) — every transport method holds a Materials and can spin up the
  cheap head per material. Note the cycle this creates across problems: a
  condensed Materials' provenance is a *solve* on a prior Materials — a DAG
  per problem, a loop across problems; the #406 save-state story should
  record that provenance.

**Geometry (the overlay).** Declares: the region partition; the
**deck identifications**; the genuine **boundary data**; and the
**region → material assignment** (direction matters: region→material is the
total function; material→regions is its multivalued inverse). The
unification ratified in round 2/3: **reflective, periodic, and rotational
"boundary conditions" are not boundary conditions — they are stage-2
geometry.** The domain is a quotient; the deck group is part of the
problem's symmetry group (feeding §2's bound a computable contribution); and
only vacuum / albedo / prescribed-inflow remain as boundary *data* — declared
here because they describe the domain's physical continuation
(method-independent), realized only at a head (no method ⇒ nothing to
resolve; the infinite problem has no boundary and no method-head boundary
work, consistently). What ships: `StructuredGeometry`
(`orpheus/geometry/structured_geometry.py`) IS the 1-D overlay — interface
positions + per-segment materials. What does not: **the d≥2 overlay object**
(today 2-D goes straight to a raw per-cell `mat_map` with no region concept)
— priced work the chain demands, an attack the user redirected onto the
current architecture rather than the proposal (round-1 attack 5).

**Mesh.** Commits the cell partition; the **spatial measure** (volumes)
appears; the spatial *axis* does not (its kind — nodal/modal — is a scheme
commitment; §2 T1). `mat_map` becomes **derived**: region map ∘ (cell →
region), well-defined exactly because of the conformity guard (§5 guard 3).
`MaterialMesh` reads correctly as this stage's object — Materials × mesh
through the pullback.

**State fields** (round-2 attack 13, fully conceded — the multiphysics
slot, empty today). Multiphysics state is a **field** — a map space→(T, ρ,
…) — on the cell partition. Cross-section resolution then becomes a lazy
pointwise map `(nuclides, library) × state → XS`, making resolved XS
*derived fields*. The placement's singular reason: state is spatial data, so
it needs the spatial partition (after mesh), and operators consume resolved
XS, so it precedes operator binding. Named now, one sentence, so that
statefulness never has to be retrofitted into the root: `Mixture` stays the
parametrization (nuclide mixture × library); the *binding of state to
space* is a chain stage.

**Method head(s).** The terminal refinement, per head: the scheme fixes the
spatial axis's kind (and any moment tail); the angular measure and axis
appear **if the method has one** (SN/PN yes; CP pre-integrates — no angular
axis ever; MC samples continuously — tally grids are a different overlay;
round-1 attack 2, conceded and folded in: the chain is a universal prefix +
per-method measure stacks); energy is formally constructed and possibly
refined (unionization). ALL axes and the solution space construct here (T1).
All remaining obligations discharge here, exact-or-refuse (§5). **Multiple
heads share one prefix** — the shipped witness is DSA: SN's head builds a
diffusion low-order system over the same meshed problem. Heads must not
mutate the prefix (sharing discipline), and per-head axes may differ
(R-raw). Cheap re-instantiable heads over a fixed prefix is also exactly
what convergence studies and the future lazy-realization solver arc want:
**realize each operator at the coarsest stage where it is measurable** — the
filtration hands the lazy arc its criterion.

## 4. The guards — where each is discharged, and why that place is right

**The schema (the near-singular reason, stated once):** a guard's home is
the **earliest stage at which its predicate is decidable** — the stage where
its last-arriving operand exists. Earlier, the predicate cannot be stated;
later, an inadmissible object has already existed. Each row below names the
last-arriving operand; that is the whole justification, per the schema.

| # | guard | discharged at | last-arriving operand | status in tree |
|---|---|---|---|---|
| 1 | energy-grid coherence **preview** (agree / no-grid / refuse; regime-scoped per #395) | Materials admission | the declared mixtures' grids | XD-4 amendment recorded (campaign plan §2.5); **preview only — implied-by-the-mint is its correctness criterion (R-mint)** |
| 2 | assigned-but-undeclared material → refuse | overlay construction | the assignment (declaration already exists) | new with the chain (user-spotted, 2026-08-25) |
| 3 | mesh conforms to regions (pullback well-defined) | mesh construction | the cells | chartered in the old CS1.5 design (§6c witness: hand-built non-conforming `Mesh1D` refused with a region-naming reason) — survives verbatim |
| 4 | deck motion is a symmetry of the discrete rule | head (BC realization) | the quadrature | **SHIPS, exact-or-refuse**: `[M]` `Quadrature.ordinate_permutation` (`directional.py:337`) — every image matches a node (no bare argmin, ERR-074), bijection (ERR-073), equal weights (ERR-042); `None` = "not a symmetry of this rule", caller refuses in the law's vocabulary; consumed by `_deck_kernel` (`realizer.py:452`), the ONE body every deck law realizes through. **No interpolation arm exists** — no silent approximation ships |
| 5 | quadrature × geometry admissibility | head | the quadrature | **the live gap**: `[M]` #398 — `SNMesh(slab_1D, lebedev)` constructs while `tests/numerics/test_registry.py:1210` asserts inadmissibility. The chain names its discharge point; the arc lands it there |
| 6 | scheme × geometry capability | head | the scheme | exists today as `supports()` machinery on the loss-representation family; re-homes with the head |
| 7 | boundary-data realization (vacuum/albedo/inflow onto trace spaces) | head | the method's trace machinery | ships (`SNBoundaryRealizer` takes the `SNMethodSpace` bundle — already one level un-welded); the assembled table's *storage* is the container question, OPEN |
| 8 | state-field × mesh conformity | state stage | the state field | future (slot named, empty) |
| 9 | energy formal resolution incl. unionization | head, per head | the method | future (R-raw makes it possible; #395's unionization arm) |
| — | spectator inertness | *nowhere — a theorem, not a guard* | — | T2: the mint reads reachable materials; a declared-unassigned material is weightless by construction. **Do not build warning machinery for it** |

**The funnel fork, resolved-in-half.** The old question "is the space mint
the unique construction funnel?" dissolves into: **admissibility is
distributed — the chain is the funnel, staged** (each guard at its earliest
decidable point, rows above). What remains genuinely open is the
**container**: whether a stage-4 aggregate class (`SNDiscretization`, or
`SNMesh` demoted to it) survives for organization + persistence (#406: no
save/restore story exists; a single dumpable object is a real argument FOR)
and as the home of caches like the realized-BC table. Ruled OPEN by the user
(2026-08-24): "We don't know yet if SNMesh or SNDiscretization (or whatever
name we give it) will dissolve." #398 is the admissibility witness whenever
that fork is ruled; the latent `BoundaryOperator`-factory option (like the
other operator factories) is on the table.

## 5. What this re-reads in the existing evidence base

**The operator mesh-independence census**
(`scratch/operator_mesh_independence_census.md`, `[M]` @ `55bb47b9` — 13
operators already mesh-free in two shapes, 7 one-step un-weldable, three
chokepoints) re-reads under the chain as follows:

- **`pole_angular_closure`'s `cls(sn_mesh)` contract** → head-side binding
  data for L; re-contract to its actual needs `(quad, reduced, coord,
  levels, ng)` — all mesh-free-available; the gap was the contract, not the
  data.
- **The L-binding bundle** (`_streaming_axes` per-axis stencils; the 1-D /
  curvilinear `dag_walk` visit iterators) → the stencils are head-side
  operator-binding data (in scope); the **traversal is strategy** (deferred
  arc) — the seam cuts through `LossRepresentation`, which holds `mesh`
  whole for both halves. The arc frees its *construction* (bundle =
  stencils + closure + spaces) and leaves traversal as the one solve-time
  handoff — which is the shape lazy realization wants anyway. The 2-D
  `SweepDependencyGraph.for_shape` (pure shape) is the proof traversal
  needs no mesh when its turn comes.
- **The realized-BC table** (`sn_mesh.bc`, assembled at
  `augmented_mesh.py:380`) → declarations move to the overlay (already on
  `Axis1D.bc` today — one stage too late under the chain); realization at
  the head (guard 7); the *table object's home* rides the container fork.

**The `MaterialXSField` verdict** (adversarial round, 2026-08-25 — main
agent's verdict at the user's invitation; unopposed and consistent with the
chain, ratify formally at the arc design): the class is `MaterialMesh`'s XS
**facade** — `[M]` `from_mesh` is `cls(materials=mesh.materials, mesh=mesh)`
(`material_xs_field.py:194-209`, adds no data); its 1155 lines decompose as
Materials content + the expansion machine (wrapping the already-free
`assemble_cell_xs`) + typed field mints + coarsening projections
(`project_through`/`_bilinear`, `:253/:317`) + **eight `apply_*` operator
kernels (~400 lines, `:741-1021`)** whose docstring confesses the mechanism
("Encapsulates the per-material loop that previously lived at
`scattering.py:405`"). Dissolution map: content → Materials; expansion →
the field-minting path at its stage (data spaces, §2); mints → honest
`CrossSectionField`s (which stay); coarsening → Materials/Mixture morphisms;
apply arms → the bound operators (CS4c's chartered "S → kernel shell").
Blast radius `[M]` 18 production + 32 test files — §6b-scale; CP/MoC/MC/
diffusion get mechanical rewiring only (sharpening-order law).

**Shape-metadata reads**: `ng`/`spatial_shape` re-point to the space's axes
— `[M]` the axes exist (`material_mesh.py:385` mints
`of_axes(energy_axis, SpaceFactorAxis("spatial", …))`;
`augmented_mesh.py:1129` mints `of_axes(angular, *scalar.axes)`); the only
missing convenience is a public axis-by-label accessor (the lookup already
lives inside `retraction`/`section`).

## 6. The rulings ledger (all user, this session, unless marked)

| # | ruling | date |
|---|---|---|
| R1 | Operators must stop requiring `SNMesh` at construction (transitive reach fine) — the arc's charter | 2026-08-24 |
| R2 | Consumed objects first (space, fields, operators); solver/strategy/traversal deferred to a lazy-realization arc | 2026-08-24 |
| R3 | The container fork is OPEN — the aggregate may survive as organization + persistence (#406); the BoundaryOperator-factory option is live | 2026-08-24 |
| R4 | The posing chain is a filtration; the architecture of §2–§4 — ratified through three adversarial rounds | 2026-08-25 |
| R5 | All axes formally construct at the method (T1); stages accumulate measures and data; no half-axis objects | 2026-08-25 |
| R6 | Naming: `Materials` (a declaration), not `Medium`; `InfiniteMedium(mixture)` in the problem layer is where the word survives; `mixture.as_infinite_medium` rejected (dependency inversion) | 2026-08-25 |
| R7 | R-mint: the method-time mint is the law; Materials-time admission is a preview whose refusals must be implied by the mint's | 2026-08-25 |
| R8 | R-raw: Materials carries raw per-material grid provenance; no early collapse; per-head axes may differ | 2026-08-25 |
| R9 | The state-fields stage exists (named now, empty today); `Mixture` = parametrization; state binds to space as a chain stage | 2026-08-25 |
| R10 | Deck identifications are stage-2 geometry, not BCs; only vacuum/albedo/inflow are boundary data; realization is head-side, exact-or-refuse | 2026-08-25 |
| R11 | Assigned-but-undeclared refuses at the overlay; declared-but-unassigned is legal and inert by T2 — no warning machinery | 2026-08-25 |
| R12 | Symmetry is the admissibility bound and quality criterion on refinements, not the flow; per-stage group machinery is aspirational (family: #152, #166) | 2026-08-25 |
| R13 | `MaterialXSField` dissolution verdict — *proposed, unopposed, consistent with R4-R6; ratify formally at arc design* | 2026-08-25 |

## 7. The adversarial record (distilled; refuted candidates are first-class output)

Round 1 (main agent attacks → user counters): **A1** twin-path bypass →
WITHDRAWN: the infinite path is the energy sub-algebra (slowing-down), a
different problem family sharing primitives; survives as the
shared-primitives condition (§3). **A2** chain overfits SN → CONCEDED by
user; universal prefix + per-method measure stacks. **A3** stage-1 "Medium"
is an inventory → resolved by the user's stronger form (a tower of *total*
problems, each complete at its symmetry) + the round-3 naming (R6). **A4**
false total order → DEFEATED by the symmetry-monotonicity argument; spawned
R12's machinery demand. **A5** d≥2 overlay missing → redirected by the user
onto the current tree; stands as priced work. **A6** energy-kind freedom →
main agent WRONG; kind is data-fixed (multigroup = indicator-modal),
determination at declaration. **A7** reflective demands mirror-symmetric
quadrature → WITHDRAWN as stated; transformed and then **measured** (guard
4: exact-or-refuse, no interpolation arm — both parties partially right,
tree cleaner than either claimed). **A8** multiphysics at the root →
transformed into the state-fields stage (R9).

Round 2: **A9** stage-1 axis mint vs the leak principle → user sharpened to
T1 (all axes at the method); forced the Library/Materials distinction. **A10**
symmetry filtration is prose → conceded aspirational (issue-tracked); the
user's challenge "find a principled alternative flow criterion" was answered
by the filtration ontology (§2), which SUBSUMES symmetry as the bound.
**A11** deck realization owes a join obligation → measured, guard 4. **A12**
`mixture.as_infinite_medium` layering → conceded, `InfiniteMedium(mixture)`.
**A13** lazy XS defers the state question → fully conceded, R9. **A14**
half-axis → no half-object was intended; the measures-accumulate/
axes-resolve spelling ratified (R5).

Round 3: **A15** mint + preview twin-check drift → fully accepted, R7.
**A16** early grid collapse (lossy-return-type at the root) → fully
accepted, R8. User's closing refinement: the amended stage-1 object is a
declaration; name it `Materials` (R6); the assigned-undeclared guard
appears (R11).

## 8. Supersessions and survivals (plan-authoring §3 — edited in place, nothing dropped)

- **Campaign plan §2.5 (CS1.5 Medium charter)**: SUPERSEDED IN SHAPE — no
  `Medium` class; no `from_medium` arms. Its surviving physics objectives
  (`kernel_and_medium_objectives.md`) map onto the chain: *medium
  expressible* → Materials + overlay; *generator lattice* → the stages as
  generator inputs (§2); *conformity* → guard 3, verbatim; *retirement
  honesty* → binds the arc's every step. The XD-4 three-outcome amendment
  survives as guard 1's shape under R7's preview discipline. The grounding
  census (`scratch/cs15_grounding_census.md`) remains the `[M]` fact base
  where still current — its site counts were already flagged STALE
  (re-censused 2026-08-24; re-census again before designing).
- **The CS1.5′ resume pointer** (campaign plan §5 tail): carries a
  supersession banner pointing here; its pickup duties survive where they
  apply (re-census; #398 as witness; surgical posture).
- **`SNMesh`/`SNDiscretization`**: not dissolved by this charter (R3 open);
  demoted-by-trajectory toward the save-state aggregate; every *machinery*
  role has a named destination stage.
- **CS2**: unchanged in substance — the mint-as-free-function, axes
  non-Optional, the angular axis — but its landing surface is now the
  method head (T1), and its identity-completion work is what converts the
  remaining `is`-checks to space-content checks.
- **CS4c**: unchanged — the apply-arm migration (R13's arms), Riesz legs,
  `DualSpace`, the `dual()` functor.

## 9. Priced open work (not designed here)

1. The d≥2 overlay object (region geometry above raw `mat_map`).
2. Per-stage symmetry groups, declared-vs-computed (aspiration; machinery
   family #152/#166; the mesher's symmetry-preservation score).
3. The container fork + #406 save-state story (+ the cross-problem
   provenance loop: condensed Materials ← a solve).
4. #398's discharge at the head (guard 5) — the arc's natural first
   admissibility landing.
5. The deferred consumer arc: lazy realization with the §3 criterion;
   traversal objects (the `for_shape` precedent generalized).

## 10. Note to the archivist (when the arc lands — not before)

The Sphinx home is a foundations theory page on problem posing (sibling of
`spaces.rst`); its spine is §2 of this file — the filtration, T1–T3, the
data-vs-solution-space distinction, and the guards schema of §4 with the
worked table. The pages it must cross-link: `spaces.rst` (the collapse pair
IS the backward walk — T3), the condensation/homogenization theory (the
same walk at the physics level), `frame.rst` (the Parseval divisor as the
projection's metric), and the boundary-conditions page (deck
identifications as geometry — R10, with guard 4's exact-or-refuse
certification and its ERR-042/073/074 lineage). The leak principle's
docstring (`material_mesh.py`) should gain its T2 name when touched. Write
the page from the *reasoning*, not the table — the guards schema ("earliest
decidable point; name the last-arriving operand") is one sentence and
regenerates every row.
