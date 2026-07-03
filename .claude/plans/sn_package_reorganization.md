# SN package reorganization — mirror `transport/`'s role-keyed layout across `numerics ⟵ transport ⟵ sn`

> **✅ DONE — MERGED to main @ `21102cc` + pushed to origin/main (2026-06-27, off branch
> `refactor/sn-package-layout`).** P1 `operators/` · P2 `loss_representation/` ·
> P3 `boundary/` · P4 `mesh/` + `geometry→augmented_mesh` rename. Per phase: collect-only + layer-imports
> (292) + behavioral subset (0 failed) + Sphinx -W + pyright (412 = baseline). **Two justified DEVIATIONS
> from this plan:** (a) `loss_representation/__init__.py` IS the 180 KB giant (NOT a thin re-export from a
> `representation.py`) — the re-export split the module namespace and broke `monkeypatch` tests
> (`transport_sweep` resolves names in its OWN module globals); giant-as-`__init__` keeps the package
> namespace == the old module namespace (monkeypatch-transparent, strictly more bit-identical). (b)
> re-baselined the `test_cell_kernel_batch` sha256 source-of-record (the `sweep_graph` doc-ref path moved
> inside `DiamondDifference`'s kernel docstrings — docstring-only, fold logic preserved). Promotion
> candidates (`spatial/`→transport, `angular`→numerics, `TransportMethod`) filed as **#272**. Full-tree gate
> + `elegance-enforcer` review running at write time.

> **Real home:** move to `.claude/plans/sn_package_reorganization.md` on the first execution turn
> (plan-mode pinned this session to `reactive-moseying-cake.md`, which held the COMPLETED
> `ReactionRateFunctional` dyad plan — overwritten here). Branch off `main` (`574cff8`, green).
> **Pure RELOCATION campaign — NO new abstractions, NO file-content splits, bit-identical per phase.**

## Context — why

`sn/` is the only one of the three transport-stack layers that is NOT organized:

- `numerics/` (L1) groups by concept: `basis/`, `quadrature/`, `spaces/` + a flat operator-algebra core.
- `transport/` (L2) is cleanly **role-keyed**: `fields/`, `operators/`, `frames/`, `mesh/`, `residuals/`,
  `source_sinks/`, `displacements/`.
- `sn/` (L3) is a **FLAT pile of 11 top-level modules** mixing operators (`operator`, `angular_operator`,
  `boundary_operator`), the sweep engine (`loss_representation`, `sweep_graph`, `sweep_schedule`),
  geometry, the realizer, the solver — plus one `spatial/` subpackage.

The campaigns (#261 operator relocation, the frame/carrier work, the `realize_recursively` co-location)
have already settled the conceptual roles; this reorg makes the **file layout match them**, giving a
consistent cross-layer mapping and making the next cross-method promotions mechanical. It is the user's
request ("operators in a subpackage; a `loss_representation` subpackage for the sweep; a reasonable mapping
across sn/transport/numerics"), with two user rulings:

1. **Rename `geometry.py → augmented_mesh.py`** — `SNMesh` *is* "the axis-primary augmented mesh"; that is
   the cross-method name (every transport method gets an `augmented_mesh`). It lands in `sn/mesh/`,
   mirroring `transport/mesh/`.
2. **Keep `spatial/` a top-level SIBLING, not nested in `loss_representation/`** — it is *spatial
   discretization*, a **promotion candidate to the transport layer** at method #2 (diffusion), exactly the
   trajectory the reaction operators took at #261. Burying it in the SN-sweep package would fight that.
   So `loss_representation/` holds only the SN-unique sweep **engine** (representation + graph + schedule).

## Target `orpheus/sn/` structure (full mirror)

```
orpheus/sn/
├── __init__.py                     # re-exports the public API from the new locations (no compat shims)
├── operators/                      # mirror transport/operators/  ⟵  numerics/operator.py
│   ├── streaming.py                #   StreamingOperator (L), InvertibleOperator (L+C)   ⟵ operator.py
│   └── boundary.py                 #   SNBoundaryOperator (B)                            ⟵ boundary_operator.py
├── loss_representation/            # the WDD/Morel–Montry SWEEP ENGINE — SN-unique (diffusion does not sweep)
│   ├── __init__.py                 #   re-exports LossRepresentation Protocol + default_for + concretes
│   ├── representation.py           #   the selectable-representation hierarchy (the 180 KB module, moved as-is) ⟵ loss_representation.py
│   ├── sweep_graph.py              #   SweepDependencyGraph                              ⟵ sweep_graph.py
│   └── sweep_schedule.py           #   SweepSchedule                                     ⟵ sweep_schedule.py
├── spatial/                        # SPATIAL DISCRETIZATION schemes — UNCHANGED, top-level sibling.
│   └── …                           #   {scheme, diamond, linear_discontinuous, cell_balance, scan,
│                                   #    pole_angular_closure, psi_half_angle_seed, sweep_cache, _ubld}
│                                   #   ⇒ promotion candidate to the transport layer @ method #2 (see flags)
├── boundary/                       # the SN BC-realization layer (1 of the 2 TransportMethod witnesses)
│   ├── realizer.py                 #   SNBoundaryRealizer + realize_recursively          ⟵ boundary_realizer.py
│   └── angular.py                  #   AngularAverageOperator, IncomingSourceOperator    ⟵ angular_operator.py
├── mesh/                           # mirror transport/mesh/ — the SN augmented-mesh + method-space layer
│   ├── augmented_mesh.py           #   SNMesh(MaterialMesh) — RENAMED (cross-method convention) ⟵ geometry.py
│   └── method_space.py             #   SNMethodSpace (realizer argument; → TransportMethod @ method #2) ⟵ method_space.py
├── solver.py                       # SNSolver, solve_sn — the driver (stays top-level)
└── solution.py                     # Solution — the typed return (stays top-level)
```

## Cross-layer mapping (the "reasonable mapping" the user asked for)

| Role | `numerics` (L1, math primitive) | `transport` (L2, neutron-aware typed) | `sn` (L3, method realization) |
|---|---|---|---|
| **Operators** | `operator.py` (LinearOperator, OperatorSum, RankOneOperator) | `operators/` (C=mult, F=fission, S=scattering) | `operators/` (L=streaming, L+C=invertible, B=boundary) |
| **Fields / carriers** | `field.py`, `vector.py`, `measure.py` | `fields/` `displacements/` `residuals/` `source_sinks/` | — (consumes transport carriers) |
| **Frames** | `frame.py` | `frames/harmonic_frame.py` | — (angular frame = scattering's eigenbasis → owned at L2) |
| **Spaces / layout** | `space.py`, `spaces/`, `face_layout`, `moment_layout` | — | — |
| **Mesh** | — | `mesh/` (material_mesh, axis, material_xs_field) | `mesh/` (augmented_mesh=SNMesh, method_space) |
| **Quadrature / basis** | `quadrature/`, `basis/` | — | — (consumes) |
| **Spatial discretization** | — | *(future promotion target)* | `spatial/` (diamond, LD, closures) — **promotable** |
| **Method realization** | — | — | `loss_representation/` (sweep engine), `boundary/`, `solver.py`, `solution.py` |

**Reading:** the operator-algebra spine (operator/field/frame) maps L1→L2→L3 with consistent subpackage
names; fields/frames/spaces/quadrature are L1/L2 (sn consumes them); mesh maps L2↔L3; the **sweep engine +
realizer + solver are SN-unique** (no L1/L2 analog — they ARE the discrete-ordinates realization).
`transport/` is already the template → **no transport changes**. `numerics/` is mostly grouped; folding its
flat operator-algebra core into a subpackage is **optional/deferred** (foundation, high blast radius, low
marginal value).

## Architectural through-line — SN-unique vs promotion-ready (keep the seams separable)

- **DONE (#261):** reaction operators (C/F/S) promoted `sn → transport/operators/`. `sn/operators/` now
  holds only the sweep-specific leaves (L, L+C, B) — they stay (the sweep is SN's way to invert streaming).
- **`spatial/` → transport layer @ method #2 (diffusion).** Spatial discretization is cross-method; kept a
  separable top-level sibling so the promotion is a clean lift, not an excavation. *File a tracking issue.*
- **`angular_operator` → `numerics/` @ 2nd consumer** (its own docstring says so). Lands in `sn/boundary/`
  now (its only consumer = the realizer); flagged.
- **`TransportMethod` / #219 `MethodSpace` ABC @ method #2.** Two witnesses (`sn/boundary/realizer.py:294-305`
  + `transport/mesh/material_mesh.py:21`). This reorg **only relocates** — it does NOT pre-build the ABC,
  the registry, or the Protocol (defer-until-2; #219 owns the sequencing).

## Phases — each a bit-identical relocation; gate after each (see Verification)

1. **`sn/operators/`** — move `operator.py`→`operators/streaming.py`, `boundary_operator.py`→`operators/boundary.py`;
   thin `operators/__init__.py`. (Smallest; the user's first idea.)
2. **`sn/loss_representation/`** — move `loss_representation.py`→`loss_representation/representation.py`
   (+ thin `__init__.py` re-export), `sweep_graph.py` and `sweep_schedule.py` into the package. `spatial/`
   **stays put**; fix the internal relative imports inside `representation.py` (`.spatial`→`..spatial`;
   `.sweep_graph`/`.sweep_schedule` stay `.`). (The user's second idea.)
3. **`sn/boundary/`** — move `boundary_realizer.py`→`boundary/realizer.py`, `angular_operator.py`→`boundary/angular.py`.
4. **`sn/mesh/` + the `augmented_mesh` rename** (largest blast — `geometry.py` is a hub, degree ~1065):
   `geometry.py`→`mesh/augmented_mesh.py`, `method_space.py`→`mesh/method_space.py`. Update **every** importer
   `orpheus.sn.geometry`→`orpheus.sn.mesh.augmented_mesh` and `orpheus.sn.method_space`→`orpheus.sn.mesh.method_space`.

`solver.py`, `solution.py`, `spatial/`, `sn/__init__.py` (content) are **untouched** except import-path updates.

## Key decisions & flags

- **No compat shims** (aggressive-retirement rule): update every importer; do NOT leave `orpheus.sn.operator`
  re-exports. Retirement-audit blast radius per moved symbol = **(1) graph callers** (Nexus `rename`/`impact`)
  **+ (2) text-grep code/tests/`docs/`** (unresolved `:mod:`/`:class:` refs render silently — the `-W` gate
  does NOT catch them) **+ (3) direct constructors**.
- **`_ubld.py` stays in `spatial/`** despite its dual role (cell primitive + `solver.py` moment helper).
  Splitting it is a separate refactor, NOT this reorg. `solver` keeps `from sn.spatial._ubld import …`.
- **`SNMethodSpace` ≠ #219 `MethodSpace` ABC** — name collision; do NOT conflate, do NOT build the ABC here.
- **The 180 KB `loss_representation.py` moves AS-IS** (→ `representation.py`); do NOT split its internals.
- **`method_space.py` placement is the softest call** — recommended in `mesh/` (the SN phase-space/method
  layer); could instead sit with the realizer in `boundary/`. Flag for review at execution.

## Verification (run after EACH phase — the tree must stay green)

- `tests/test_layer_imports.py` (AST layer gate) — **no new cross-layer leak** (numerics ⊀ transport ⊀ sn).
- Full tree SERIAL: `.venv/bin/python -O -m pytest -m "not slow" -p no:xdist --timeout=300 -p no:cacheprovider`
  → **0 failed** (the new green baseline; xdist unstable).
- **Sphinx `-W` clean** + V&V matrix regenerates (docs `:mod:`/`:class:` refs to the moved modules updated).
- `npx pyright orpheus/` ≤ **412** (the CLI is the oracle; `<new-diagnostics>` LSP import-noise is #226).
- Per-phase: Nexus `rename`/`impact` for the moved symbols BEFORE the move (the importer list).

## Critical files / references

- Explorer maps (this session): `.claude/plans/reactive-moseying-cake-agent-ae23666e7c5b1d8af.md` (sn/ groupings,
  cluster boundaries, awkward-placement flags) + `…-aae68a783e1362794.md` (cross-layer map, import-direction
  proof, transport-vs-sn operator line, MethodSpace/TransportMethod end-state).
- `tests/test_layer_imports.py` (L0–L4 FORBIDDEN_EDGES + WHITELIST + the TYPE_CHECKING tolerance).
- `transport/operators/__init__.py` + `integral_kernel_operator.py` (the transport-vs-sn operator line, in prose).
- `sn/boundary_realizer.py:294-305` + `transport/mesh/material_mesh.py:21` (the TransportMethod two-witnesses).
- `docs/theory/{discrete_ordinates,operator_algebra,loss_representations}.rst` (module cross-refs to update).

## Scope / discipline

- **Pure relocation**: no new abstractions, no file splits, no logic change — bit-identical per phase. The
  cross-method promotions (spatial→transport, angular→numerics, TransportMethod) are SEPARATE future tasks
  (file tracking issues; do NOT start them here).
- Execution: mechanical moves are `method-implementer`-eligible per phase, BUT the wide import blast + the
  user's structural steering favor **main-agent-direct per phase** with an **`elegance-enforcer` review of the
  final structure**. Per-phase ff-merge; `main` always green.
- Commit only when the user asks; stage explicitly (no `git add -A`); trailer
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. NEVER `git checkout`/`restore`/`stash`
  on uncommitted files. User rejects `# type: ignore` (`cast` OK).
- **Sequencing vs P5:** this reorg and P5 (#50 energy condensation, greenfield) are independent. Doing the
  reorg FIRST gives P5 a clean `sn/` to land in; either order is fine — user's call.
