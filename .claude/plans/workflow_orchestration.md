# Dynamic workflows: the operating model, the templates, and wave 2 as workflows (living plan, opened 2026-10-22)

## Goal

Large reorganizations run as dynamic workflows: scripted fan-outs of sub-agents, with a fixed division of labour by model. The user's direction (2026-10-22): "the dynamic workflow is orchestrated by you (Opus), advised by Fable and executed by Sonnet"; "we will probably use it a lot more from now on". This plan fixes the roles, the stage templates (saved, reusable scripts under `.claude/workflows/`), the invariants every script keeps, and wave 2 of the layer DAG campaign expressed as a sequence of those workflows.

Precedent:
- `[M]` Layer DAG round 1 (2026-10-21/22, three units) and wave 1 (2026-10-22, WP-B, P4b, P4c) used build, review and amend in isolated worktrees.
- The wave-2 design used two drafts and two attacks: 4 agents, about 44 min, 1.2 M sub-agent tokens.

## The roles, by model

| role | model | does | never does |
|---|---|---|---|
| **Orchestrator** | Opus (the main loop) | writes and launches each workflow; owns the user's rulings, the plans, integration (cherry-pick or rebase onto one branch), the full bar (suite, pyright, Sphinx, `dead_references`, census), commits on `main`, push, CI; synthesizes design panels | delegates a ruling, a merge or a push; answers a running workflow agent by `SendMessage` (it spawns a duplicate writer, memory `reference_sendmessage_resumes_duplicate`) |
| **Advisor** | Fable | every JUDGEMENT stage inside a workflow: design attacks (cross-domain-attacker, elegance-enforcer), review verdicts (qa, elegance-enforcer), adjudication of findings, the cold read of a specification, the verification that each amendment answered its finding; the orchestrator's `advisor()` at its own decision points | writes production code |
| **Executor** | Sonnet | every DOING stage: censuses, measurements, design drafts, gate implementation (test-architect), builds, amendments, mechanical rewrites, documentation drafts (archivist) | rules on a design question (it returns the question in `NEEDS`) |

Why the split: judgement is where a stronger model's independence pays (an adversarial first pass that a weaker reviewer would wave through). Execution is long, tool-heavy and checkable, so it is gated by tests and by the advisor's review rather than by the executor's own care. The orchestrator keeps the context that sub-agents lack (the rulings, the history, the user), which is why integration and every ruling stay with it.

## Invariants every workflow script keeps

1. **Isolation by worktree** for every stage that writes tracked files in parallel (`isolation: 'worktree'`). A follow-up stage on the same unit (amend) works in that worktree by path, not in a new one.
2. **Read-only reviewers.** A reviewer mutates only a COPY (`cp -a` or a detached worktree), prints `orpheus.__file__` before any probe, and runs targeted tests only.
3. **Serial pytest per worktree; no full suite inside a workflow.** The full suite (about 75 min) runs once, at integration, by the orchestrator. A child process reads source from disk, so no file is edited while a suite runs.
4. **Gate first.** A build stage writes its done-when gate before the change and records the gate's red at the base tip (the command and the failing line).
5. **Return contracts.** A word cap; a schema when the next stage parses the result; detail goes in files under `scratch/<campaign>/`; a `NEEDS` block for anything unobtainable. Every brief names the plan sections it executes and the rulings that bind it (`workflows` rule, invariant 4).
6. **No ruling inside a workflow.** A design question an agent meets goes to `NEEDS`; the orchestrator brings it to the user between workflows. A user checkpoint is therefore a boundary BETWEEN workflows, never a pause inside one.
7. **Size.** The session guideline is medium (under 10 agents per workflow). A wave larger than that is split into consecutive workflows, the orchestrator reading each result before launching the next. The user can raise the guideline with "Dynamic workflow size" in `/config`.
8. **The bar is the orchestrator's.** pyright is vacuous inside `.claude/worktrees` (#640), and Sphinx stamps the commit it builds. So both run at integration, on the integrated tip.

## The templates (`.claude/workflows/`, parameterized by `args`)

| template | stages | when |
|---|---|---|
| `design-panel` | N executor drafts from distinct angles (parallel) → M advisor attacks, each reading every draft (parallel) | the ontology is being searched; the orchestrator synthesizes, and the user rules |
| `spec-check` | one executor drafts a specification from given inputs, with measurement agents in parallel for its open items → an advisor cold read (elegance) and an advisor structural check (cross-domain) in parallel | turning ruled decisions into one executable specification |
| `build-units` | per unit: executor test-architect writes the gates (red) in a worktree → executor implementer turns them green → advisor qa and advisor elegance in parallel → executor amend → advisor verify of each disposition | additive or file-disjoint units in parallel |
| `docs-pass` | per page section: executor archivist drafts → advisor qa verifies the claims (each `[M]` re-run) | a theory page after a carve |

Each template's `args` carry the units or sections, their briefs, the shared rules block, and the output directory. The scripts are reviewed by the orchestrator before their first use, and their `meta.description` states the stages.

**Invocation.** `Workflow({scriptPath: "<abs>/.claude/workflows/<template>.js", args: {...}})`, with `args` as a JSON object, not a string.
- `[M]` 2026-10-22: by name, `Workflow({name: "design-panel"})` first answered "not found. Available: deep-research". About an hour later, in the same session, the four templates appeared in the harness's list of available skills by name. So the registry picks the directory up with a delay; prefer `scriptPath` the first time.
- `[M]` By `scriptPath`, the template ran, and `design-panel` with `angles: []` raised its own guard with 0 agents spawned.
- `[M]` The per-agent model overrides work inside a workflow, AND they win over an agent type's own `model:` frontmatter:
  - bare agents returned "Fable 5.1" for `model: 'fable'` and "Sonnet 5.5" for `model: 'sonnet'`;
  - `agentType: 'qa'` with `model: 'fable'` returned "Fable 5.1", and `agentType: 'test-architect'` with `model: 'sonnet'` returned "Sonnet 5.5", although every project agent's definition says `model: opus`;
  - the control, `agentType: 'qa'` with no override, returned "claude-opus-5-5[1m]".

  So a template stage's model is the one it names. A stage that omits `model` runs on its agent type's definition, which is Opus for every project agent.
- All four templates pass `node --check` after their `export` is stripped.

## Wave 2 as workflows (the boundary carve and arm T)

Design: `boundary_law_ontology.md`, the twelfth exchange, with the user's rulings D1 to D5 and the strict `xfail` for D4. Order: `layer_dag.md`, the late-evening compaction point. Wave 1 integrates first; WP-B's final face names (`geometry/faces.py`) are inputs to WF2.1.

**WF2.1 — one specification** (`spec-check`; about 5 agents).
- **Executors, in parallel, measuring:**
  - (i) diffusion's vacuum admission: is Marshak `Approximated` or `Realized`? Measure, with the theory page's realizability tiers.
  - (ii) the colouring symmetry's float rule: the slab `(0, .1, .2, .3)` is the control; which comparison makes it symmetric, and which refuses a genuinely asymmetric slab?
  - (iii) the per-method admission column for every shared response and every deck.
- **Executor:** drafts `.claude/plans/wave2_boundary_carve.md` from design A (types) and design B (consumer census), the rulings, and (i) to (iii). It contains every type, field, invariant, module path and verb signature; the typed realized operators; the realizer-table spelling; the leaf units with their file sets; and the dependency table (`retirement-audit` G.24).
- **Advisors, in parallel:** an elegance cold read and a structural check.
- **Orchestrator:** applies the findings. **User checkpoint:** the specification.

**WF2.2 — the leaves** (`build-units`; units in two stages by dependency). Each unit costs up to 6 agents: gates, build, qa, elegance, amend, verify. So two units already make 12, against the medium guideline of under 10 per workflow: either run one workflow per stage, or raise the guideline.
- Stage A, parallel:
  - **L1:** `geometry/deck.py` (the closed deck family, its Enum classification, the face-role type `DeckTransformation | int`) and the colouring symmetry;
  - **L2:** `transport/boundary_conditions/` conditions (vacuum, albedo with a required shape, inflow with frame-typed profiles), the façades (`AlbedoBoundary.white()`), the normalized combination with its cone check, `BoundaryConditions`, the binding with its invariance check.
- Stage B, after L2:
  - **L3:** the per-method realizer tables (class keyword plus `seal()`), admission (`Realized`, `Approximated`, `Refused`), one table read by the query and the realizer;
  - **L4:** the typed realized operators (the zero map, a scaled ordinate permutation, a rank-one re-emission) with their properties;
  - **L5:** `WallAmplitudes` and `BoundaryCondition.from_amplitudes(specular, diffuse)`.
- Nothing consumes the leaves, so the suite stays green at every unit.
- **Orchestrator:** integration and the full bar; then push.

**WF2.3 — the carve, sequential on ONE branch** (`build-units` with one unit per workflow, each unit's worktree the same branch's tip; the user, 2026-10-22: "do whatever is required for the boundaries first, so that we don't need to keep old spelling working"; "not everything needs to be parallelized"). There is no dual spelling. The geometry and mesh API changes outright, and the tree is allowed to be red on the branch between slices. It is green again, and merged, only after the last slice (and WF2.4) and the full bar.
1. **S1: geometry and mesh.** `face_roles` replaces `boundaries`, and `FaceRoles` replaces `FaceLaws`. The 25 hollow inner faces become tags. The slice's own targeted tests (geometry, mesh) are green; every other package is expected red.
2. **S2:** transport and SN.
3. **S3:** diffusion.
4. **S4:** cp, mc and moc (the three `BC_REGISTRY` tables; CP's vacuum admitted as realized, with the #546 strict `xfail`).
5. **S5:** derivations, specification and reference (closed values; the translation reversed).

Mechanics: S1 runs `build-units` with no `worktree`, which makes the branch's worktree. Each later slice passes that path as its unit's `worktree`, so the gates and the build continue on the same branch (the template option, added 2026-10-22). Each slice starts from the previous slice's tip and is reviewed (Fable qa and elegance) before the next starts, so a defect is caught in the slice that made it. Tests move with their package. The orchestrator reads each slice's result before launching the next.
- **User checkpoint, after S5:** the meaning-change sites (109: reflective on a curved face becomes a specular wall; periodic on a curved face becomes unspellable), listed with each one's new spelling.
- Digests (geometry, mesh, specification) re-key ONCE, at the merge. There is no dual spelling, so they re-key only once.

**WF2.4 — the retirement, one unit, on the same branch after S5** (`build-units`, a single unit, the dependency table's order, leaves first, `retirement-audit` A to G). It deletes `BC`, `BoundaryTraceLaw`, the string tables, `RegistryMixin`'s `key=` on the laws, the kernel tier, `LawSum`/`LawScaled`, `geometry/boundary/`'s response half, and the welded verbs. It updates the error catalogue, `_OPEN_INTRA_PACKAGE_CYCLES` and the layer gate in the same commit. A single unit because a retirement split across commits leaves orphans (G.24).
- **Orchestrator:** the full bar; `layer_census.py` (P1's SCC expected gone, `[R]`); push.

**WF2.5 — documentation** (`docs-pass`, one unit per section of `docs/theory/foundations/boundary_conditions.rst`, 8590 lines, and the SN boundary page).

## Ruled (the user, 2026-10-22)

1. **The surgical-carve exception is recalibrated, as a trial during the reorganization.**
   - The user: "models are getting better ... every once in a while we do need to try new things and recalibrate older rulings. We're doing that during reorganization, because it's less about development and more about putting things in place. It's a good way to test it. Also Sonnet 5.5 is much better than Sonnet 5, so we can consider it for well scoped work."
   - So wave 2 is executed by Sonnet units, each well scoped, with gates first, Fable reviews, and the user's checkpoints between workflows.
   - The `workflows` rule's exception (the ruling of 2026-05-20) stays in force outside this trial. After wave 2 the orchestrator reports how the trial went: defects caught by the gates and the reviews, defects that reached integration, and the cost. The rule's source page (`docs/development/workflows.md`, through the harness process) is changed only on that evidence.
2. **The size guideline stays medium** (under 10 agents per workflow). Wave 2 runs as several smaller workflows, evaluated as they go; "large" may come later if the results hold.
