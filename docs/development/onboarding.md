---
harness:
  kind: onboarding
  budget_tokens: 3000
---

# ORPHEUS — start here

ORPHEUS (Open Reactor Physics Educational University System) is a set of
reactor-physics solvers that is a teaching tool and a high-stakes analysis
tool at once: the multigroup transport methods (discrete ordinates, collision
probability, method of characteristics, Monte Carlo) and diffusion, with
cross-section processing, fuel behaviour, thermal hydraulics and kinetics
around them. Code and documentation are one corpus: every equation with its
derivation, its conventions, its design rationale and the test that verifies
it. A single maintainer develops it by steering AI agents, and the
documentation is the knowledge those agents work from.

## The two halves: a Problem and a Solution

**The Problem half poses.** The material mesh (`MaterialMesh`: the axes, the
region-to-material map, the group count, the cell volumes and the
cross-section field) is method-agnostic data, and a method's problem is that
data plus the method's behaviour: `SNProblem` adds the quadrature, the spatial
scheme, the angular closure and the boundary laws; `DiffusionMesh` adds the
scalar trace and the realized albedo laws; `HomogeneousProblem` is posed from
one mixture on the energy axis alone ("Problem" is a concept with three
realizations and no shared type yet). Posing is a filtration: materials,
geometry, mesh, state and the method head each commit one refinement of phase
space, and every axis resolves at the head. The hub owns its spaces, its
fields and its bound operators, and its last step is the operator pencil
`(A, F)`, the loss and the fission production, with `at(σ) = A − σF`, from
which the two questions are typed: `EigenPosing`, the homogeneous question
`Aψ = μFψ` read through a spectral map, and `SourcePosing`, the affine question `Aψ = q`; a subcritical
multiplying source is `SourcePosing(pencil.at(1), q)`, a composition, never a
third type.

**The Solution half answers.** The Strategy is a value and not the
Problem's: `Splitting.from_schedule` labels each loss term implicit or
explicit, so `A = M − N` (the splitting's `M`, the implicit part, and the
lagged `N`) derived through the algebra's own sums, and the Strategy inverts
only `M`: `M.inverse()` is the sweep, or the block back-substitution on a
carrying mesh, which `SourceIteration` applies as `ψ ← M⁻¹(q + Nψ)` and
`KrylovAcceleration` hands to GMRES as the matvec; one `power_iteration` is
the outer loop for every iterative family, and the 0-D baseline solves its
pencil directly. Every level mints an
`IterationRecord` whose verdict is derived and never stored; a best-effort
exit is legal and made audible, and a claimed convergence is re-measured
before it is believed. The answer is fused with its question and its gauge
into an outcome, and the `Solution` is the frozen five-tuple problem,
outcome, strategy, certificate, record, on which the domain operations
(reaction rates, homogenisation, condensation, importance) are posed. The
adjoint entries dagger the same objects: `.H` is `♯ ∘ dual ∘ ♭`, built from
the bound spaces and propagated through sums, products, tensor products and
inverses by law, so no driver spells a transpose.

## The vocabulary

Everything in the problem layer is a map on spaces: a field is space →
values, an operator is space → space, so the organising object is the space
with labelled axes. On an axis: a discrete measure (weights on nodes; a
quadrature is one), a basis (nodal or modal) and the frame pairing them,
which induces the space's metric and mints the analysis and reconstruction
faces together (on a Galerkin frame the two faces are each other's adjoint
up to one scalar); the scheme does the
same for the spatial axis; both are stage-2 generators, read for their
induced data and then forgotten. Between spaces: retraction and section
along an axis (a split pair, `R∘E = id`; "embedding" is not an operator name
here), trace restriction and its scatter transpose on the boundary, system
restriction and extension by zero, and the lift of a bulk action onto
bulk ⊕ trace. The metric is an object, the Riesz legs ♭ and ♯ are operators,
the adjoint is `♯ ∘ dual ∘ ♭`, and there is no `.T`. Flux lives in the
positive cone, a predicate. A symmetry group is realised, never tabulated;
invariance is the measure's question; an orbit space is named by its
stabiliser; the lift from an orbit space is the Reynolds projector. Each concept with its class
and the theory page that defines it, and the list of what is still
hand-rolled: [the conceptual view](../architecture/conceptual_view.rst).

## The architecture map

A module's home is the lowest-knowledge layer whose vocabulary suffices to
define it, and imports flow only from more knowledge to less;
`tests/gates/test_layer_imports.py` enforces the layers (the criterion:
[layering](../architecture/layering.rst)).

| layer | packages | knows |
|---|---|---|
| L0 | `derivations/` | symbolic and high-precision references (Branch 1); below L2, may import `numerics/`, `geometry/`, `data/` |
| L1 | `numerics/` | mathematics only: spaces, measures, quadrature, operators; no neutrons |
| input | `geometry/`, `data/` | meshes and boundary conditions; nuclear data |
| L2 | `transport/` | the transport vocabulary every method shares; method-agnostic |
| L3 | `sn/`, `diffusion/`, `homogeneous/`, `cp/`, `moc/`, `mc/`, `kinetics/`, `fuel/`, `thermal_hydraulics/` | one method's machinery each; no L3 package imports another |
| L4 | `plotting.py` | orchestration; consumes everything |

The spine is the operator algebra `A = L + C − S − N_2n − B`: a solver is
spelled as that expression over typed operators, and the code reads like the
mathematics (the master standard of the `coding-elegance` skill). The theory
corpus mirrors the layering, and its reading order is
[the theory index](../theory/index.rst): conventions first, then the
foundations, the method sub-books and verification. `docs/architecture/`
records the layering and the conceptual view; `docs/development/` (this
section) records how the project is built; `docs/api/` is the reference.

## The ontology is the work

The main difficulty of developing ORPHEUS is finding the right ontology and
keeping ontological discipline: the objects the code is made of, each in its
correct form, its right place and its right shape, so that the derived
concepts fall out of the algebra and the mistakes become unspellable.
Emerging the concepts is part of the agent's core duties (Cardinal Rule 1);
what the discipline means, and what a weld is, is Cardinal Rule 2; the two
kinds of development it implies are under "How a session runs".

## The direction of development

- The six core folders — `data`, `geometry`, `numerics`, `transport`, `sn`,
  `homogeneous` — are where every quality campaign lands first.
- Discrete ordinates, diffusion and the homogeneous baseline are sharpened
  first. Collision probability, characteristics and Monte Carlo are taken
  later, in campaigns that recycle the machinery and the vocabulary (the MoC
  ray is a Volterra operator), and are not improved on their own terms
  meanwhile: a defect found there is filed with its measurement unless it
  blocks a sharp method or is a harmonisation onto shared machinery.
- `thermal_hydraulics` and `kinetics` may leave the repository; no
  architectural investment there until that is decided.
- The lens on every design ruling: build the machinery and realise the
  operator algebra (Cardinal Rule 2), and effort is never a criterion (the
  cardinal page's preamble).

## How a session runs

- The SessionStart hook prints the protocol: the environment-health gates,
  then the core batch (`vv-principles`, `coding-elegance`,
  `instrument-doctrine`, the lessons index, the Nexus briefing). Python is
  `.venv/bin/python`; the canonical test invocation is `python -O -m pytest`.
- Open issues are the plan (Rule 4); a campaign's plan file lives in
  `.claude/plans/`, written to the `plan-authoring` rule. Enter plan mode for
  any task of three or more steps.
- Discussion comes before the plan, and the plan before the code: no
  implementation starts on a vague plan unless the user says so. Two kinds
  of development are told apart by whether the ONTOLOGY is known. When it is
  known (the objects exist and the work re-spells, extends or verifies them),
  plan with little input and enter plan mode. When it is being searched (the
  right objects are not yet named), open a living plan in `.claude/plans/` at
  the first exchange and refine it in place as the discussion runs; the user
  prompts the concepts out and helps organise them, and implementation waits
  until the plan is polished (the form: `plan-authoring` §0).
- Work on a branch `<type>/<topic>` (`feature|fix|docs|refactor|test|chore`);
  commits follow Conventional Commits (`<type>(<scope>): <summary>`) and close
  issues with `Closes #NN` in the body; `main` is always green and receives
  only `--ff-only` merges; delete the branch after merging. A hook refuses
  `git add -A` and a commit on `main`; the reason and the riders are
  `process-discipline`. The full workflow:
  [git workflow](git_workflow.rst).
- A feature is done when the tests, the theory page and the Nexus graph agree
  (Rule 3).

## Where the rules, workflows and agents are

The rules in `.claude/rules/` load into every session and every Key-agent
dispatch, some only when a path their front matter names is touched; they are
generated from `docs/development/rules/` (one, the routing rule
`nexus-tools`, is installed by `nexus setup` from sphinxcontrib-nexus), with
the founding cases in
`docs/development/evidence/`. Skills load on demand or per agent. The agents,
their roles, the dispatch invariants and the seven workflows are the
`workflows` rule; how to change a rule, a skill or an agent, and what each
costs: [the harness page](harness.md).
