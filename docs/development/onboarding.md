---
harness:
  kind: onboarding
  budget_tokens: 1500
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

## The architecture map

A module's home is the lowest-knowledge layer whose vocabulary suffices to
define it, and imports flow only from more knowledge to less;
`tests/test_layer_imports.py` enforces the layers (the criterion:
[layering](../architecture/layering.rst)).

| layer | packages | knows |
|---|---|---|
| L0 | `derivations/` | symbolic and high-precision references (Branch 1), below L1 |
| L1 | `numerics/` | mathematics only: spaces, measures, quadrature, operators; no neutrons |
| input | `geometry/`, `data/` | meshes and boundary conditions; nuclear data |
| L2 | `transport/` | the transport vocabulary every method shares; method-agnostic |
| L3 | `sn/`, `diffusion/`, `homogeneous/`, `cp/`, `moc/`, `mc/`, `kinetics/`, `fuel/`, `thermal_hydraulics/` | one method's machinery each; no L3 package imports another |
| L4 | `plotting.py` | orchestration; consumes everything |

The spine is the operator algebra `A = L + C − S − N_2n − B`: a solver is
spelled as that expression over typed operators, and the code reads like the
mathematics (the master standard of the `coding-elegance` skill). The theory
corpus mirrors the layering: read `docs/theory/conventions/` first (symbols,
normalisation, indexing, the crosswalk to the literature), then
`docs/theory/foundations/` (the operator algebra, frames and projection, the
boundary law, data, measures, geometry), then the method sub-books under
`docs/theory/methods/`, and `docs/theory/verification/` for the V&V
principles, the error catalogue and the generated V&V matrix.
`docs/architecture/` records the layering; `docs/development/` (this section)
records how the project is built; `docs/api/` is the reference.

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
  operator algebra; a welded, unnamed operation is a failure to do so.

## How a session runs

- The SessionStart hook prints the protocol: the environment-health gates,
  then the core batch (`vv-principles`, `coding-elegance`,
  `instrument-doctrine`, the lessons index, the Nexus briefing). Python is
  `.venv/bin/python`; the canonical test invocation is `python -O -m pytest`.
- Open issues are the plan (Rule 4); a campaign's plan file lives in
  `.claude/plans/`, written to the `plan-authoring` rule. Enter plan mode for
  any task of three or more steps.
- Work on a branch `<type>/<topic>` (`feature|fix|docs|refactor|test|chore`);
  commits follow Conventional Commits (`<type>(<scope>): <summary>`) and close
  issues with `Closes #NN` in the body; `main` is always green and receives
  only `--ff-only` merges; delete the branch after merging. A hook refuses `git add -A` and a commit on `main`; the reason and the
  riders are `process-discipline`. The full workflow:
  [git workflow](git_workflow.rst).
- A feature is done when the tests, the theory page and the Nexus graph agree
  (Rule 3).

## Where the rules, workflows and agents are

The rules in `.claude/rules/` load into every session and every Key-agent
dispatch, some only when a path their front matter names is touched; they are
generated from
`docs/development/rules/`, and the founding cases sit one link away in
`docs/development/evidence/`. Skills load on demand or per agent. The agents
in `.claude/agents/`, their roles, the dispatch invariants and the seven
workflows are the `workflows` rule, with the brief template in
[workflows](workflows.md). How to change a rule, a skill or an agent, and what
each costs: [the harness page](harness.md).
