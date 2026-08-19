# #261 — Collapse `CollisionOperator` + relocate the cross-method reaction operators to `transport/`

> **Real home:** move to `.claude/plans/archive/issue_261_collision_collapse_relocation.md` on the first execution
> turn (plan-mode restricts edits to the active plan file; the prior Frame-projection content lives at
> `.claude/plans/archive/frame_projection_machinery.md`). Branch `refactor/operator-inverse-algebra`. Surgical,
> main-agent-direct, user-steered; **NO `method-implementer`**. Host env → `.venv/bin/python` (3.14).

---

## ✦ STATUS — #261 COLLAPSE + RELOCATION + `transport/operators/` COMPLETE, COMMITTED @ `bbe8a51` (2026-06-26)

Reconcile against git, not this line. All four phases + the **`transport/operators/` subpackage** are
**committed** as ONE coherent commit `bbe8a51` on `refactor/operator-inverse-algebra` (63 files, 671+/588−; the
per-file collapse+relocate edits interleave, so no clean split was possible). Final full-suite acceptance gate:
**2807 passed / 7-and-only-7 baseline reds, zero regressions**; both 0-ULP canaries green; pyright `orpheus/`
412; Sphinx `-E -W` clean. `.claude/lessons.md` L29 (the invariant-nesting lesson) + this plan are NOT in the
commit (the `.claude/*` instruction-architecture flow handles them). **NEXT = P4 — the `ProductionRateFunctional`
fate (deferred to a post-compaction discussion; see the P4 section below). Also: W-F (paused) now lands the
`@overload`/ndarray honesty in scattering/fission's new `transport/operators/` home.**

- **P0 — `MultiplicationOperator` gap-fill.** Added `space: FunctionSpace | None` field + `domain`/`codomain`
  props + `from_mesh(σ, mesh, *, space=None)` classmethod (defaults `space` from the mesh — the faithful
  `CollisionOperator` drop-in). Additive, pyright 412 neutral, 0-ULP canary green.
- **P1 — `CollisionOperator` COLLAPSED into `MultiplicationOperator` (class DELETED, 170 lines).** It added
  nothing genuine; its one real addition (the W-D `domain`/`codomain`) was a base GAP, filled in P0. Production:
  3 `solver.py` sites → `MultiplicationOperator(coefficient=…, space=sn_mesh.full_field_space)`;
  `InvertibleOperator` keys on `MultiplicationOperator`, **invariant re-anchored on
  `streaming.sn_mesh is diagonal.coefficient.mesh`** (the H2 correction — geometric-consistency, NOT the
  weaker shape-equality the `OperatorSum` guard checks; the strengths nest
  `object-identity ⊋ geometric-consistency ⊋ shape-equality`); `sigma`→`coefficient.values`; the L+C dispatch
  is one-directional on `StreamingOperator.__add__` (the `C + L` symmetry test retired — transport can't
  dispatch onto sn). `# type: ignore`→`cast`. 29 test files migrated (`CollisionOperator(m,σ)`→`from_mesh(σ,m)`);
  intrinsic-property test added (`tests/transport/test_multiplication_operator.py::TestSpaceMetadataAndGuardJoin`,
  5 passed — space→domain/codomain + the guard-join). Docs: archivist re-pointed 10 operator.py + 2 extra
  solver.py + 19 `.rst` + 3 test-prose refs, Sphinx `-W` clean. **Gate: 2098 passed full-suite + 834/7-and-only-7
  combined-state, pyright 412, C6 assert_type pins intact (CLI).**
- **P2 — `FissionOperator` → `orpheus/transport/fission.py`** (pure `git mv`, zero sn-imports). Importers +
  doc refs updated; `sn/__init__` re-export retired (unused). Gate: layer-imports 308 passed, pyright 412.
- **P3 — `ScatteringOperator` + `LegendreMomentScattering` → `orpheus/transport/scattering.py`.** SNMesh import
  absolutized under `TYPE_CHECKING` (the established `scalar_flux.py` L2↔L3 pattern, tolerated by the layer gate);
  **`mesh.quad`→`self.quadrature`** (bit-identical cleanup — S holds its own quadrature). Blanket module-path
  rename `orpheus.sn.scattering`→`orpheus.transport.scattering` (89 refs / 27 source files). Gate: layer-imports +
  **0-ULP scattering canary GREEN** + scattering suites 382 passed + pyright 412 + **Sphinx `-E -W` build succeeded**.

**Net architecture:** the cross-method reaction operators (collision = `M[σ_t]`, fission, scattering) now live in
`transport/` (L2); only `StreamingOperator` + `InvertibleOperator` + `LossRepresentation` (the WDD sweep) remain
SN-specific. `transport/` is runtime sn-free (TYPE_CHECKING SNMesh tolerated, per the carrier pattern). The four
operators live in **`transport/operators/`** (`multiplication_operator`, `fission`, `scattering`,
`integral_kernel_operator` + the package `__init__`); `ProductionRateFunctional` (a Functional, not an operator)
stays at `transport/` top level — its fate is **P4**. A whole subclass retired; the base honestly generalised.

**Flagged follow-ups (out of #261 scope, pre-existing):** `operator_algebra.rst` §"C/S/F stay at domain=None"
reversed by W-D (doc-debt); `index_convention.rst` §"Typed source types" stale `IsotropicSource`/`PerOrdinateSource`
cluster; ~10 test files with `CollisionOperator` prose in docstrings (cosmetic, don't affect build); the
`test_fission_operator.py:97` `\i` SyntaxWarning.

---

## P4 — `ProductionRateFunctional` fate (DEFERRED — post-compaction discussion)

Surfaced while deciding the `transport/operators/` contents (the user's *"why does it exist? what problem does
it solve?"* probe — the same instinct that retired `CollisionOperator`). **It is a Functional (flux→scalar), NOT
an operator** — so it correctly stays OUT of `operators/`; it sits at `orpheus/transport/production_rate_functional.py`.
The open question is its FATE.

**Findings (code-verified):** `ProductionRateFunctional` is the §5.6 functional `φ ↦ p(r) = Σ_g' νΣf·φ` (per-cell
fission emission DENSITY — group axis contracted to one scalar per cell). Created in **#257 S5** to be the named
middle of the fission decomposition `F = M_χ ∘ ProductionRate ∘ M_νΣf`, with the stated intent that **S6 would
rewire F's kernel through it** — but **S6 never landed**: `F.apply` still uses the FUSED `RankOneOperator`
(`χ ⊗ νΣf`). Consequences:
- **Production-dead.** `ProductionRateFunctional` + the `FissionOperator.production_rate` property are consumed
  ONLY by tests (`test_functional_category`, `test_integral_kernel_category`, the 0-ULP `test_fission_kernel_crosscheck`).
- **Weak oracle.** Its own docstring admits `evaluate` reproduces the kernel's inner line *"byte-for-byte — the
  same numpy primitive, the same axis, the same keepdims=True"* → the cross-check it backs is PROCEDURALLY
  identical (a twin-consistency check), NOT structurally independent (real fission correctness = the
  `k∞ = νΣf/Σa` eigenvalue tests). Per `coding-standards` fuller-view-oracle discriminator, a procedurally-
  identical "oracle" already covered by a structurally-independent one is genuine redundancy.
- **Sole concrete `Functional`.** It is the ONLY concrete instance of the §5.6 Functional category (the
  `numerics/functional.py` exemplar) — retiring it leaves the `Functional` Protocol instance-less.

**The three options (DISCUSS after compaction):**
1. **Retire** — production-dead + procedurally-redundant oracle (aggressive-retirement). Cost: the `Functional`
   category loses its only concrete instance (becomes Protocol-only until a real consumer arrives).
2. **Wire it (finish S6)** — decompose F's fused `RankOneOperator` into the explicit `M_χ ∘ ProductionRate ∘
   M_νΣf` composition, matching scattering's composed-operator-as-production (`R∘Λ∘M`). Makes it load-bearing +
   reads-as-math, at the cost of 3 objects where 1 fused one works. **W-F-adjacent** (the shared
   `Operator[Flux,SourceSink]` emission abstraction; W-F is paused, lands in the new `transport/operators/` home).
3. **Keep as-is** — the named per-cell-fission-rate diagnostic + the §5.6 category exemplar.

**Template:** the `CollisionOperator` collapse (this carve) is the model — challenge "what problem does it solve?",
classify genuine-vs-redundant, then either retire or complete-the-intended-wiring. The `cross-domain-attacker` is
the right reviewer for the structural call (it confirmed the F = rank-1-degenerate-of-S's-frame reading).

---

## Context — why

A dependency-challenge pass this session (2 explorers + cross-domain-attacker, all verdicts code-verified)
established two things:

1. **`CollisionOperator` adds nothing genuine over its base `MultiplicationOperator`.** Its only real
   addition — `domain`/`codomain` (`= sn_mesh.full_field_space`, for the W-D composition guard) — is a
   **GAP in `MultiplicationOperator`** (the base is mesh-free and returns `None` spaces). Everything else
   C adds is session-born or misplaced: the `sn_mesh` field (the base reads the mesh off the carrier at
   apply time — `multiplication_operator.py:118,206`), the legacy ndarray-wrap ctor (production passes a
   `CrossSectionField`), the `sigma` alias (`= coefficient.values`), the `L+C` `__add__` dispatch (belongs
   on `StreamingOperator`). The class docstring literally says *"Collision **IS** a MultiplicationOperator."*
   Per Cardinal Rule 2 + type-vs-property (L-004), **C collapses into the base**.

2. **C/F/S are cross-method (method-agnostic) reaction operators** (D5) stranded in `sn/` for historical
   reasons. After the collapse, C *is* `MultiplicationOperator` (already in `transport/`). F has **zero**
   sn-internal imports; S has **one** removable `TYPE_CHECKING` `SNMesh` hint. So the three relocate cleanly.

**Outcome:** the base gets honestly general (any multiplier can join the W-D guard); a whole subclass
retires; the cross-method reaction operators become honestly transport-level. Genuinely SN-specific
operators (`StreamingOperator`, `InvertibleOperator`, `LossRepresentation` — the WDD sweep) stay in `sn/`.

**User decisions (this session):** full delete C + a `MultiplicationOperator.from_mesh` classmethod for the
legacy/test ndarray-wrap; **full carve this session** (P0→P3); transport destination = **flat**
(`transport/fission.py`, `transport/scattering.py`, matching the existing `transport/multiplication_operator.py`).

## Load-bearing correction (cross-domain-attacker H2 — do NOT regress)

The L↔C invariant `streaming.sn_mesh is diagonal.sn_mesh` (`operator.py:787`) is **NOT** redundant with the
W-D guard. `FullFieldSpace.__eq__` is `(name, shape)` only — the strengths nest
`object-identity ⊋ geometric-consistency ⊋ shape-equality`, and the guard enforces only the **weakest**.
The sweep threads C's σ against **L's geometry** (`loss_representation.py:1690`), so two equal-shape /
different-volume meshes would pass the guard yet compute wrong physics. **Re-anchor on geometric
consistency, NOT shape-equality:** `streaming.sn_mesh is diagonal.coefficient.mesh` (the `CrossSectionField`
carries `.mesh` — `cross_section_field.py:85`). **Execution check:** confirm
`mat_xs.total_cross_section_field.mesh is sn_mesh` in production (else find the right handle).

---

## Phases (each: green + CLI-pyright ≤ baseline (target DOWN) + Sphinx -W + layer-imports gate; ff-merge)

### P0 — Fill the `MultiplicationOperator` gap (additive, bit-identical)
`orpheus/transport/multiplication_operator.py`:
- Add `space: "FunctionSpace | None" = field(default=None)` (after `coefficient`; `FunctionSpace` → the
  `TYPE_CHECKING` block). Add `domain`/`codomain` properties returning `self.space` (endomorphic — same
  space both sides).
- Add `@classmethod from_mesh(cls, sigma: "np.ndarray | CrossSectionField", mesh, *, space=None)` — wraps a
  bare ndarray into a `CrossSectionField` on `mesh` (the C legacy/test path), sets `space`.
- **Gate:** the 3 existing direct-`MultiplicationOperator` tests stay green (all pass `coefficient=` only →
  `space=None`, never read `domain`/`codomain`); pyright neutral; 0-ULP canary; `tests/transport` green.

### P1 — Collapse `CollisionOperator` (the retirement) — **proactive `test-architect` FIRST**
- **Production migrations** (`orpheus/sn/`):
  - `solver.py:223/936/1013`: `CollisionOperator(sn_mesh, mat_xs.total_cross_section_field)` →
    `MultiplicationOperator(coefficient=mat_xs.total_cross_section_field, space=sn_mesh.full_field_space)`
    (name the local `C`). Drop the `CollisionOperator` imports (`solver.py:51,213`).
  - `InvertibleOperator` (`operator.py`): `isinstance(diagonal, CollisionOperator)` (`:782`) →
    `MultiplicationOperator`; annotations `:757/775/819` → `MultiplicationOperator`; `diagonal.sigma`
    (`:793/844`) → `diagonal.coefficient.values`; **re-anchor the invariant** `:787` →
    `streaming.sn_mesh is diagonal.coefficient.mesh` (per the H2 correction above).
  - `StreamingOperator.__add__` (`:483`): `isinstance(other, CollisionOperator)` → `MultiplicationOperator`.
  - **DELETE** `class CollisionOperator` (`:488-655`) + the `__all__` entry (`:156`).
- **Dispatch / symmetry:** the `L+C → InvertibleOperator` dispatch stays **one-directional** on
  `StreamingOperator.__add__` (production is always `L + C`). Retire the lone `C + L` symmetry test
  (`test_invertible_operator.py:138-150`) — preserving it needs a `transport→sn` layer violation (a
  `__radd__` alone can't, the base `__add__` shadows it). Document `L + C` as canonical.
- **Test migration** (29 behavioral files, zero `isinstance`): `CollisionOperator(sn, σ)` →
  `MultiplicationOperator(coefficient=…, space=…)` or `.from_mesh(σ, mesh, space=…)`. **`test_collision_operator.py`**:
  migrate the genuinely-C-specific cases (space-for-guard, `from_mesh` wrap) into
  `tests/transport/test_multiplication_operator.py`; **delete the rest as redundant** with the base suite.
  Rewire helpers (`_test_helpers.py:290/317/340/343`, the `wave_t_t4` fixture `-> CollisionOperator`
  annotation), the sweep/solve/verification behavioral files, `test_operator_block_role.py`.
- **Docs** (archivist): re-point 19 refs — `operator_algebra.rst` (13, incl. `:meth:CollisionOperator.apply/solve`
  → `MultiplicationOperator`), `index_convention.rst` (4), `loss_representations.rst` (2). The §5.7
  "C = M[σ_t]" framing stays (now literally true).
- **Gate:** 0-ULP `test_scattering_kernel_crosscheck` canary; `tests/sn/{operators,sweep,solve,eigenvalue,verification}`;
  pyright DOWN; Sphinx -W; the intrinsic-property test the test-architect designs (a `MultiplicationOperator`
  with a `space` joins the `L+C` guard = positive; a mismatched-`coefficient.mesh` diagonal reds the
  `InvertibleOperator` construction = negative; mutation-confirmed).

### P2 — Relocate `FissionOperator` → `transport/fission.py` (pure move)
- `git mv orpheus/sn/fission.py orpheus/transport/fission.py` (zero code change — F has no sn-internal imports).
- Update importers: `solver.py:46`, `sn/__init__.py:1`, 4 test files; re-point doc `:class:` refs
  (`integral_kernel_operator.py:105/126/186`, `iteration.py:57`, `operator.py:392`, `loss_representation.py:1974`).
- **Gate:** `tests/test_layer_imports.py`; fission tests; pyright; Sphinx -W.

### P3 — Relocate `ScatteringOperator` → `transport/scattering.py`
- **Eliminate the SNMesh dependency first** (make `transport` fully sn-free): replace
  `mesh.quad.weights.sum()` (`scattering.py:791`) → `self.quadrature.weights.sum()` (S already holds
  `self.quadrature`); widen / drop the `mesh: SNMesh` annotations (`:748/772`) to the carrier's transport
  mesh type; verify `ScalarSourceSink/AngularSourceSink.zeros_on(mesh)` accept the carrier mesh. Goal: no
  `from .geometry import SNMesh` at all (`:187`). *Fallback if a residual `SNMesh` hint remains:* keep it
  `TYPE_CHECKING`-guarded as an absolute `from orpheus.sn.geometry import SNMesh` (tolerated by the layer
  gate, `test_layer_imports.py:148`).
- `git mv orpheus/sn/scattering.py orpheus/transport/scattering.py` (incl. `LegendreMomentScattering`).
- Update importers: `solver.py:55`, `sn/__init__.py:2`, 7 test files; fix the now-stale comment
  `integral_kernel_operator.py:127` ("…live in the L3 sn package").
- **Gate:** layer-imports; **0-ULP scattering canary** (load-bearing here); scattering tests; pyright; Sphinx -W.

### Re-exports (`sn/__init__.py`)
Decide at execution (low-risk): retire the `FissionOperator`/`ScatteringOperator` re-exports (aggressive
retirement — internal-only usage), or keep one cycle re-exporting from the new `transport` location. Lean retire.

---

## Verification (per phase + end-to-end)
- `.venv/bin/python -O -m pytest tests/numerics tests/sn tests/transport -p no:xdist --timeout=300 -p no:cacheprovider -q -rfE`
  — SERIAL (xdist unstable); the **0-ULP `test_scattering_kernel_crosscheck`** canary is load-bearing
  (P1/P3 must keep it 0-ULP); baseline = the **7-and-only-7 pre-existing reds**.
- `npx --no-install pyright --outputjson orpheus/` — the ORACLE (NOT the streamed `<new-diagnostics>` LSP,
  #226); net DOWN (C deletion + the W-D space honesty). Baseline this session = **412 errors**.
- `tests/test_layer_imports.py` — the architectural gate (transport must not import sn at runtime;
  `TYPE_CHECKING` L3-in-L2 tolerated at `:148`).
- Sphinx -W clean throughout (the 19 doc re-pointings + the F/S `:class:` path updates).
- **Proactive `test-architect` before P1** (operator-algebra retirement crossing the 0-ULP-adjacent sweep);
  **archivist** for the doc re-pointing (P1/P2/P3).

## Critical files
- `orpheus/transport/multiplication_operator.py` (P0 gap-fill + `from_mesh`)
- `orpheus/sn/operator.py` (delete `CollisionOperator`; `InvertibleOperator`/`StreamingOperator` updates; invariant re-anchor)
- `orpheus/sn/solver.py` (3 construction sites + imports)
- `orpheus/sn/fission.py` → `orpheus/transport/fission.py`; `orpheus/sn/scattering.py` → `orpheus/transport/scattering.py`
- `orpheus/sn/__init__.py` (re-exports)
- tests: `tests/sn/operators/test_collision_operator.py` (migrate-then-delete), `test_invertible_operator.py`
  (retire the `C+L` test), `tests/transport/test_multiplication_operator.py` (absorb C-specific cases), ~27
  behavioral files (mechanical), `tests/test_layer_imports.py` (gate)
- `docs/theory/{operator_algebra,index_convention,loss_representations}.rst` (archivist, 19 refs)

## Scope / discipline
- Surgical, main-agent-direct; NO `method-implementer`. Per-phase ff-merge; `main` always green. Commit only
  when the user asks; stage files explicitly (no `git add -A`); trailer
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. NEVER `git checkout`/`restore`/`stash`
  on uncommitted files (monkeypatch to revert mutation probes). User rejects `# type: ignore` (`cast` OK).
- **Durable architecture banked:** the minimal symmetric shape for a cross-method reaction operator = typed
  cross-section data + optional `space` + mesh-from-carrier + `from_solver_data`/`from_mesh`. C = `M[σ]`
  (no subclass); F = rank-1 `χ⊗νΣf`; S = `R∘Λ∘M` + its SO(3) eigenbasis frame. The shared flux→source
  morphism is the existing `IntegralKernelOperator` Protocol — **no `ReactionOperator` ABC** (L-004).

## W-F note (paused)
P4.5 **W-F** (retire the legacy bare-`ndarray` arms in scattering/fission + realign the `@overload`
`TimedFullField→FullField` key; the "shadowing" caveat was REFUTED — a #226 LSP artifact; the arms are LIVE
at the K-eigenvalue outer-loop boundary, task #76) is **paused** for this carve. After P3, W-F's honest-typing
lands in scattering/fission's **new transport home**.
