# The onboarding page gains the conceptual view, and the harness gains the ontology duty — a living plan

**Status: DRAFT, opened 2026-09-21 in the first exchange, before any edit.** This
plan is refined in place as the discussion runs (§3: a refuted premise is edited,
never dropped). Implementation starts only when the user rules the plan polished.
Every class name and every number below is `[R]` until the census of §5.2 marks
it `[M]`.

Source pages: `docs/development/onboarding.md` (generates CLAUDE.md; kind
`onboarding`, budget 1 500 tokens, ≈1 400 used `[M]` CP#7), the rules under
`docs/development/rules/`, the harness page. Nothing under `.claude/` is edited
by hand.

## §1 The goal, separately from any means

- **G1.** A fresh session, and every Key dispatch, reads CLAUDE.md and comes
  away knowing (a) that ORPHEUS is two halves, a Problem and a Solution, and
  what each carries; (b) the mathematical concepts the code is built from, in
  their correct form, each with the page that defines it; (c) that finding the
  right ontology and keeping ontological discipline is the main difficulty of
  the work, and that emerging the concepts is part of the agent's core duties.
- **G2.** The harness states how development runs: discussion before the plan
  and the plan before the code; two kinds of development, told apart by whether
  the ontology is known; a living plan opened at once in the searching kind; no
  implementation on a vague plan unless the user says otherwise.

Non-goals: no code change; no change to the five cardinal rules or their
ranking; no new agent; no new always-on rule file.

## §2 The user's words, 2026-09-21 — the source of every claim in §4

- *"the code is currently split into Problem (which is the result of eventually
  augmenting the material mesh with all the information necessary to fully
  define the Space in its various axis, and it the object that essentially
  carries the problem posing (either an eigen pencil or an affine system…). Then
  there is the solution side, which is everything that deals with the solution
  trajectory (resolvent, splitting, etc), and we obtain a Solution object, which
  carries the problem, its choice of solution trajectory and the results."*
- *"how in the Problem side we establish the axes of space, build frames, how
  the code uses strong mathematical concepts like Discrete Measures, Basis
  functions, Frames for analysis and reconstruction. We rely on concepts such as
  restriction, section, embedding, etc. We use Riesz operators (raise and lower)
  and we write the algebra in such a way that things like adjoint are not
  handrolled, instead they fall for free from the mathematical rigour and
  disciplined implementation (or so is our aim. anything hand-rolled in
  production is a debt that needs refinement)."*
- *"(1) We should establish that 'it is part of your core duties to emerge the
  hard mathematical concepts in their correct form'… the more pro-active you are
  in emerging the mathematical concepts we need, the more you help me, as I'm
  not a mathematician."*
- *"(2) we need to discuss in order to plan proper implementation. A plan needs
  to be well defined, refined and polished before implementation, and we shall
  not jump into implementation with a vague plan (unless I actively say
  otherwise). There are 2 types of development… the main discriminator between
  them is if the ontology is known or we're searching for it… When we're in this
  vague territory you should start drafting a living plan, that we refine as we
  discuss."*
- *"the main difficulty of development of ORPHEUS is typically finding the right
  ontology, and following the ontological discipline for everything to end up in
  the right place and with the right shape so that concepts fall for free from
  the algebra (such as the adjoint), recognizing when something has been
  implicitly welded (multiple concepts have been mixed into a value) and needs
  to be unwelded to express the concepts individually, to obey the standards of
  coding and elegance so that the right ontology makes mistakes unspellable and
  most guards unnecessary. And that this should be a major focus… you need to
  do this actively because there is little example of scientific code with this
  level of ontological rigour in your training corpus."*

## §3 The forks — the more principled option first, effort in its own sentence

- **F1 Where the full conceptual view lives.** (a) A new corpus page,
  `docs/architecture/conceptual_view.rst`, holding the full view with every
  concept linked to the theory page that defines it, while CLAUDE.md carries a
  compressed version and the link. More principled: one definition per concept
  (X4), the page is RST so a `:class:` role is a reference `dead_references`
  checks, and CLAUDE.md stays within a budget. (b) Everything in CLAUDE.md.
  Effort: (a) is one more page and one more Sphinx toctree entry.
  **Recommendation: (a).**
- **F2 The CLAUDE.md budget.** The additions of §4 as drafted measure ≈1 530
  tokens by the generator's estimate (`[M]` 2026-09-21, chars / 3.6: 4.1 ≈490,
  4.2 ≈260, 4.3 ≈572, 4.4 ≈207); the block would read ≈2 900 against a budget
  of 1 500, and the always-on block ≈15.4K by probe against the ≤15K target
  re-adopted on 2026-09-21 (13 872 today `[M]`). Two trims are X4-clean, since
  each duplicates a page that owns the content: the theory reading order in
  "The architecture map" (owned by `docs/theory/index.rst`, ≈120 tokens) and
  the "Where the rules…" section's restatement of the harness page (≈150).
  Options: (a) budget 3 200, both trims, and measure: ≈15.1K, within the
  probe's own noise of the target; (b) the same and the target moves to ≤16K,
  recorded on the harness page with the number. **Recommendation: (a), then
  (b) only if the probe says so.**
- **F3 Where the ontology duty lives.** (a) A new CLAUDE.md section, "The
  ontology is the work", read first by every session and every Key dispatch,
  which is where the direction of development already lives. (b) A sixth
  cardinal rule. Refuted: the cardinal rules rank *what outranks what* and name
  the rule that operationalises each; the duty is a posture and a difficulty
  statement, and Rule 2 (architecture) already outranks implementation gain,
  so a sixth rule would restate Rule 2's lens. (c) A new rule file. Refuted: a
  rule is a floor of checks and tells; this content has one check (ask "what
  is this, mathematically?" before writing) and belongs beside the
  operator-algebra lens it sharpens. **Recommendation: (a).**
- **F4 Where the two-modes ruling lives.** The trigger ("open a living plan at
  the first exchange when the ontology is being searched") must be always-on,
  so it goes in CLAUDE.md's "How a session runs"; the living plan's FORM goes
  in `plan-authoring` (path-scoped to `.claude/plans/**`, so it loads the moment
  the plan file is opened) as a §0. The workflows rule is untouched: W1–W7 route
  agents, and the discussion is the main agent's with the user.
  **Recommendation: as stated.**
- **F5 One definition of "weld".** The onboarding page today: *"an operation
  computed inline and never named as an operator (welded into its caller)"*.
  The user today: *"multiple concepts have been mixed into a value"*. These are
  two kinds of the same failure, an operation weld and a value weld. Proposal:
  one definition on the onboarding page covering both, cited from
  `coding-elegance` and the agent memories by the word, never restated.
- **F6 Register.** The onboarding page is impersonal; ruling (1) is addressed to
  the agent ("your core duties"). Proposal: the new section speaks to the
  reader in the second person, because its content is a duty, and the rest of
  the page stays as it is. Open: the user may prefer one register throughout.
- **F7 The reading order of CLAUDE.md.** Today: identity; architecture map;
  direction; how a session runs; where the rules are. Proposal: identity; the
  two halves and the vocabulary (what the code IS); the architecture map (where
  it lives); the ontology is the work (the duty); the direction; how a session
  runs; where the rules are. Open.

## §4 The content, drafted

### 4.1 CLAUDE.md — "The two halves: a Problem and a Solution"

Verified against `scratch/_claude_md/explorer_problem_side.md` (`[M]`
2026-09-21, tree `b4865b5e`; every class name below is at the memo's
`module:line`). The Solution paragraph is `[PENDING explorer_solution_side.md]`.

**The Problem half poses.** The material mesh (`MaterialMesh`: the axes, the
region-to-material map, the group count, the cell volumes and the
cross-section field) is method-agnostic data, and a method's problem is that
data plus the method's behaviour: `SNProblem` adds the quadrature, the spatial
scheme, the angular closure and the boundary laws; `DiffusionMesh` adds the
scalar trace and the realized albedo laws; `HomogeneousProblem` is posed from
one mixture on the energy axis alone. Posing is a filtration: materials,
geometry, mesh, state and the method head each commit one refinement of phase
space, and every axis resolves at the head. The hub owns its spaces, its
fields and its bound operators, and its last step is the operator pencil
`(A, M)` with `at(σ) = A − σM`, from which the two questions are typed:
`EigenPosing`, the homogeneous question `Aψ = μMψ` read through a spectral
map, and `SourcePosing`, the affine question `Aψ = q`; a subcritical
multiplying source is `SourcePosing(pencil.at(1), q)`, a composition, never a
third type.

`[R]` the memo's verdict on the user's sentence: right in shape, under-counted
in tiers ("Problem" is a concept with three realizations and no shared type;
the augmentation is two-staged, data then behaviour). The page says so in one
clause rather than hiding it.

**The Solution half answers.** The Strategy is a value and not the
Problem's: `Splitting.from_schedule` labels each loss term implicit or
explicit, so `A = M − N` with `M` and `N` derived through the algebra's own
sums, and the resolvent is `M.inverse()` (the sweep, or the block
back-substitution on a carrying mesh), which `SourceIteration` applies as
`ψ ← M⁻¹(q + Nψ)` and `KrylovAcceleration` hands to GMRES as the matvec; one
`power_iteration` is the outer loop for every family. Every level mints an
`IterationRecord` whose verdict is derived and never stored; a best-effort
exit is legal and made audible, and a claimed convergence is re-measured
before it is believed. The answer is fused with its question and its gauge
into an outcome, and the `Solution` is the frozen five-tuple problem,
outcome, strategy, certificate, record, on which the domain operations
(reaction rates, homogenisation, condensation, importance) are posed. The
adjoint entries dagger the same objects: `.H` is `♯ ∘ dual ∘ ♭`, built from
the bound spaces and propagated through sums, products, tensor products and
inverses by law, so no driver spells a transpose.

`[M]` 2026-09-21, `scratch/_claude_md/explorer_solution_side.md`: 0 of 24
`apply_transpose` call sites are in a driver module; 48 of 67 concrete
`LinearOperator` subclasses own an `apply_transpose`, 8 inherit one, 11 have
none (the eager `.H` gate refuses those); the algorithm choice (source
iteration or Krylov) is still a string on the entry, recorded only as the
record's label. The page states the aim ("no driver spells a transpose") and
the docs page carries the census.

### 4.2 CLAUDE.md — "The vocabulary" (one paragraph; the per-concept links live on the page)

Everything in the problem layer is a map on spaces: a field is space →
values, an operator is space → space, so the organising object is the space
with labelled axes. On an axis: a discrete measure (weights on nodes; a
quadrature is one), a basis (nodal or modal) and the frame pairing them,
which induces the space's metric and mints the analysis and reconstruction
faces together (a Galerkin frame promises `M* = R`); the scheme does the
same for the spatial axis; both are stage-2 generators, read for their
induced data and then forgotten. Between spaces: retraction and section
along an axis (a split pair, `R∘E = id`; "embedding" is not an operator name
here), trace restriction and its scatter transpose on the boundary, system
restriction and extension by zero, and the lift of a bulk action onto
bulk ⊕ trace. The metric is an object, the Riesz legs ♭ and ♯ are operators,
the adjoint is `♯ ∘ dual ∘ ♭`, and there is no `.T`. Flux lives in the
positive cone, a predicate. A symmetry group is realised, never tabulated;
invariance is the measure's question; an orbit space is named by its
stabiliser; the lift is the Reynolds projector. Each concept with its class
and the theory page that defines it: [the conceptual view](../architecture/conceptual_view.rst).

`[M]` the theory map (`scratch/_claude_md/explorer_theory_map.md`): of the 37
terms asked, 8 have no definition anywhere under `docs/` (basis, product
space, bound operator, half-trace, Riesz raise and lower, symmetry group,
weld, hub), 4 are defined only under an unlabelled heading (frame, Problem,
resolvent, strategy), 3 carry several senses (realization, section,
splitting). Those are §7's owed list.

### 4.3 CLAUDE.md — "The ontology is the work" (rulings (1) and (3))

The main difficulty of developing ORPHEUS is finding the right ontology: the
mathematical objects the code is made of, each in its correct form, in its
right place and with its right shape, so that the derived concepts fall out of
the algebra for free (the adjoint of an operator bound to its spaces is
`G_V⁻¹ Aᵀ G_W`, never a second implementation) and the mistakes become
unspellable, which is what makes most guards unnecessary; a guard that remains
is elegance debt, tagged and owed a retirement. The discipline has three parts.

- **Emerging the concepts is part of your core duties**, in their correct
  mathematical form, and proactively. The maintainer is a nuclear engineer, not
  a mathematician; the concepts the code now rests on (frames, discrete
  measures, the cone, the posing filtration, quotients by a symmetry group)
  emerged from the agent when the right question was asked. Ask that question
  before it is asked of you: *what is this, mathematically, and what structure
  is it an instance of?* Scientific code with this level of ontological rigour
  is rare in your training corpus, so what you would write unprompted is the
  procedural transcription; the standard here is the other thing.
- **Recognising a weld.** A weld is several concepts mixed into one value or one
  inline computation, so that none of them can be named, typed, tested or
  reused on its own. Unwelding spells each concept as its own object and
  re-composes them in the algebra. The tells: a hand-rolled adjoint or
  transpose; an `axis=` parameter on a mathematical object; an index remap at a
  call site; a length stored where a codomain space belongs; an operation
  computed inline and never named as an operator.
- **Obeying the standards** so that the ontology, once right, is enforced by
  the types: illegal states unrepresentable, conventions fixed at the
  definition site, the code reading as the mathematics (`coding-elegance` is
  the ceiling, `coding-standards` the floor).

### 4.4 CLAUDE.md — "How a session runs" gains the two modes (ruling (2))

- Discussion comes before the plan, and the plan before the code: no
  implementation starts on a vague plan unless the user says so. Two kinds of
  development are told apart by whether the ONTOLOGY is known. When it is
  known (the objects exist and the work re-spells, extends or verifies them),
  plan with little input and enter plan mode. When it is being searched (the
  right objects are not yet named), open a living plan in `.claude/plans/` at
  the first exchange and refine it in place as the discussion runs; the user
  prompts the concepts out and helps organise them, and implementation waits
  until the plan is polished (the form: `plan-authoring` §0).

### 4.5 `plan-authoring` — a new §0, "A living plan, when the ontology is being searched"

A plan is opened at the first design exchange, not after the design settles,
whenever the right objects are not yet named (the onboarding page, "How a
session runs"). It is refined in place as the discussion runs, and carries from
its first draft: the goal in the domain's terms (§1, §5); the candidate
ontologies, each with what it makes unspellable and what it leaves welded; the
refuted candidates with their structural reason (`process-discipline`); the
rulings ledger, dated; and the condition under which implementation starts,
which is the user's ruling that the plan is polished. Every claim carries its
marker from the first draft (§2).

- check: a design discussion past its first exchange with no plan file open.
- tell: a plan written after the design was settled, which is a record, not a
  plan.

### 4.6 The corpus page `docs/architecture/conceptual_view.rst`

Title: "The conceptual view: a Problem, a Solution, and the mathematics
between them". Register: the architecture section's (the layering page).
Every code name is a `:class:` or `:func:` role, so `dead_references` checks
it; every concept ends in a `:ref:` to the label the theory map found, and a
concept with no label ends in a sentence saying so, with the nearest page
(no new definition is written here: the page maps, it does not define).
Sources: the three memos, cited by path and date at the top. Outline: (1) fields are maps space → values and
operators maps space → space, so the organising object is the space with
labelled axes; (2) the Problem half, stage by stage (the posing filtration:
materials, geometry, mesh, state, method head; T1–T3); (3) the objects on an
axis (discrete measure, basis, frame; the cone); (4) frames and schemes as
stage-2 generators; analysis and reconstruction; (5) restriction, section,
embedding, trace, the boundary law; (6) the metric as an object, Riesz raise
and lower, the adjoint that falls out; (7) symmetry and quotient; (8) the
Solution half: splitting, resolvent, realization, the iteration and its
contract, the Solution object and the domain operations posed on it; (9) what
is still hand-rolled, as the debt list with its issues. Each item ends in the
theory page that defines it; the page defines nothing itself. Plus one
table, `concept | class | module | theory anchor`, assembled from the two
code memos' glossaries (57 + 25 rows, deduplicated), and the debt list with
its `[M]` date and the counts from §4.1.

## §5 Instruments and gates

- 5.1 `python -m tools.harness --check` (the onboarding block's budget; the
  citation check), the keep − omit probe (the always-on total; recorded on the
  harness page beside 13 872), `sphinx -E -W`, `dead_references`, the 48 harness
  tests under `-O`, pyright.
- 5.2 A census of every class, function and module named in §4's prose: an AST
  walk over `orpheus/` confirming each name exists at the cited path (X2;
  positive control: a name known to exist and one known not to). The
  Markdown page cannot carry a checked reference; the RST page can, and every
  code name there is a `:class:`/`:func:` role.
- 5.3 The planted check of the campaign: qa recounts the concept list against
  the theory corpus (every concept's defining page exists and says what the
  line claims), as the second reader.

## §6 Steps, once the plan is ruled polished

S1 the page `docs/architecture/conceptual_view.rst` and its toctree entry; S2
the onboarding edits (4.1–4.4, F5, F7); S3 the `plan-authoring` §0; S4
generate, `--check`, Sphinx, `dead_references`, the harness tests; S5 the probe;
S6 the qa recount; S7 commit on the branch, ff-merge, push, delete; S8 memory
(`project_harness_context_budget.md`, the index line) and the #477 comment.

## §7 Owed, and open questions for the user

- **Two present-tense-false docstrings the Problem-side explorer found** (`[M]`
  2026-09-21), fixed on sight per `articulation` §6, as one separate commit on
  this branch before the merge: `orpheus/geometry/boundary/_base.py:143-165`
  (`BoundaryTraceLaw` says `geometry_map` and `response_kernel` are
  unpopulated and "phase B1 mints" the factor types; 7 of 7 laws override both
  and the types exist at `_factors.py:350`, `:445`) and
  `orpheus/numerics/space.py:30-45` (announces four "not shipped" space classes
  whose content landed as axes; 0 of 4 exist).
- **The theory corpus's own drift, fixed on sight in this branch** (`[M]`
  2026-09-21, the theory map): `foundations/spaces.rst` L4410–4420 and
  `foundations/frame.rst` L914–926 say the Riesz legs are not built (they are:
  `numerics/space.py:818`, `:838`; `sn/history.rst` L1872 records the landing);
  `architecture/layering.rst`'s important-box ("reserved, not yet implemented")
  predates `OperatorPencil`, `EigenPosing`, `SourcePosing` and `SNProblem`;
  Γ± is defined under two labels (`trace-sign-predicate`,
  `trace-half-decomposition`), one equation twice (X4).
- **Definitions the corpus owes** (the theory map): labels for the four
  load-bearing definitions under unlabelled headings (frame, symmetry groups,
  Problem, Strategy and resolvent) so the page can `:ref:` them, small; a
  definition of the Riesz legs (the only account is `sn/history.rst`), owed
  where the metric object is defined; "basis", "bound operator", "half-trace",
  "hub" and "weld" defined once each (F5 is the weld). Open: which of these
  this campaign writes and which it files.
- **A question for the theory page, not a defect claim** (the Solution memo's
  hypothesis 11, read at `orpheus/sn/solver.py:2972` and `problem.py:1190`):
  the adjoint fixed-source entry poses `A† ψ* = q*` over the pure loss, while
  the forward source posing is always `pencil.at(1) = A − F`; so there is no
  daggered multiplying-source entry, `(A − F)† ψ* = q*`, and on a fissile hub
  the importance the adjoint entry returns omits fission regeneration. Whether
  the pure-transport importance is the intended posing or the daggered
  multiplying source is owed: a user ruling, then an issue.
- **The debt list the page will name** (the memo's "not yet typed", 13 items):
  no shared `Problem` type over the three hubs; `DiffusionMesh` mints neither
  pencil nor posing (the solver assembles them); no α posing minted (`ALPHA_MAP`
  defined, unused); the restriction siblings share no base by ruling; the
  rank-d spatial axis is generator-less (the CS2 seam, gated); `HomogeneousProblem`
  lives in the solver module pending the Problem → Solution carve;
  `GeneratingMeasure` has no theory section. Which of these carry an issue:
  `[PENDING a gh issue list pass at S8]` (Rule 4).

- The forks F5–F7 above.
- Concepts named in the user's message with no defining page in the corpus
  `[PENDING the theory-map memo]`.
- The hand-rolled debt list `[PENDING the solution-side memo]`: which items
  already have issues, which need one (Rule 4).

## §8 Rulings ledger (dated; the user's unless marked)

| # | ruling | date |
|---|---|---|
| R1 | The plan is opened before any edit; the discussion refines it. | 2026-09-21 (this plan's own premise, from ruling (2)) |

## §9 Sources

- `scratch/_claude_md/explorer_problem_side.md` (57-row glossary, the pipeline
  in 7 paragraphs, 13 not-yet-typed items, counts with commands);
  `explorer_solution_side.md` (glossary, 8 paragraphs, 11 not-yet-typed items,
  the adjoint census, the debt census); `explorer_theory_map.md` (26 pages,
  the 37-term anchor table, 11 Development-history pages, the reading order,
  4 drift finds). All at `b4865b5e`, 2026-09-21. Untracked; the durable parts
  land on the corpus page.
