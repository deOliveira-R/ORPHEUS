# L20 three-surface dependency audit — renaming `InvertibleOperator`

**Target**: `class InvertibleOperator` — `orpheus/sn/operators/streaming.py:445`
**Action**: RENAME (not deletion).
**Audit standard**: `.claude/rules/coding-standards.md` "Retire as you go" (blast radius =
graph + text-grep-incl-docs + direct-constructors) + explorer lesson L34 (segment/string
assembly).

## Provenance of this audit

| Item | Value |
|---|---|
| Branch | `refactor/operator-naming-honesty` |
| `git diff --name-only main...HEAD` | **EMPTY** — branch is at `main` HEAD `8654d348`, zero divergence |
| Nexus graph validity | **VALID** — graph built for `main`, tree identical |
| Uncommitted edits in audited subsystem | **NONE** (`git status` shows only untracked `scratch/`, `.claude/plans/*`, and one modified explorer memory file) |
| Line numbers current at | HEAD `8654d348`, this session |

Working-tree scope note: `.claude/worktrees/nexus-workspace-wiring/` is a **separate git
checkout** with its own copies of ~10 test files containing the symbol. It is **excluded**
from every count below — it is not this checkout's blast radius, but it will need its own
rename when that worktree rebases.

---

## HEADLINE (durable structural claim, before any line number)

`InvertibleOperator` is the **within-group resolvent `(L+C)`** — the single sweep-invertible
composite of the SN operator algebra `A = L + C − S − F − B`. It is produced by exactly ONE
dispatch (`StreamingOperator.__add__`), and it is the type on which THREE structural gates key:

1. **`SweepOperator.is_adjointable`** — `isinstance(inner, InvertibleOperator)` decides whether
   the inverse operator carries a reverse-scan transpose (`ScheduledInvertibleOperator` arm → no).
2. **`SweepOperator.apply_transpose`** — the same `isinstance` as the direct-call backstop.
3. **`ScheduledInvertibleOperator.__init__`** — `isinstance` guard on its first operand.

The symbol is **NOT exported** (`streaming.py:114` `__all__ = ["StreamingOperator"]` only), so
the rename cannot break a public import path — but it is imported by name at **5 production
sites** and is spelled **100 times in `docs/`**, of which **83 are Python-domain roles that the
Sphinx `-W` gate will NOT catch**.

**Total surface (this checkout, excluding worktree):**

| Bucket | Files | Line hits |
|---|---:|---:|
| Production code (`orpheus/`) | 12 | 96 |
| Tests (`tests/`) | 34 | 148 |
| `docs/**.rst` | 19 | **100** |
| `derivations/` | 2 | 5 |
| `.claude/**` (plans, memories, skills, lessons) | ~105 | 510 |
| `scratch/` | 1 | 14 |
| `examples/`, `tools/` | 0 | 0 |

---

## SURFACE 1 — Graph callers (NECESSARY BUT NOT SUFFICIENT)

Nexus `impact(direction=upstream, max_depth=2)` on
`py:class:orpheus.sn.operators.streaming.InvertibleOperator` reports **`total_affected = 274`**
(class degree 117). Depth-1 production/derivation nodes:

| Node | file:line | Edge kind |
|---|---|---|
| `orpheus.sn.operators.streaming.StreamingOperator.__add__` | `orpheus/sn/operators/streaming.py:414` | `type_uses` + `references` (**the sole producer**) |
| `orpheus.sn.operators.streaming._require_typed_composite` | `orpheus/sn/operators/streaming.py:119` | `references` (docstring) |
| `orpheus.sn.operators.scheduled_invertible.ScheduledInvertibleOperator` | `orpheus/sn/operators/scheduled_invertible.py:87` | `references` |
| `…ScheduledInvertibleOperator.__init__` | `scheduled_invertible.py:105` | `type_uses` |
| `…ScheduledInvertibleOperator.invertible` | `scheduled_invertible.py:135` | `type_uses` |
| `…ScheduledInvertibleOperator._solve_timed_full_field` | `scheduled_invertible.py:210` | `references` |
| `orpheus.sn.operators.sweep_operator.SweepOperator` | `orpheus/sn/operators/sweep_operator.py:83` | `references` |
| `…SweepOperator.is_adjointable` | `sweep_operator.py:130` | `references` |
| `orpheus.sn.solver.SNSolver` | `orpheus/sn/solver.py:883` | `references` |
| `orpheus.sn.solver.SNSolver._solve_krylov` | `orpheus/sn/solver.py:1678` | `references` |
| `orpheus.sn.solver._select_si_resolvent` | `orpheus/sn/solver.py:707` | `type_uses` ×2 |
| `orpheus.sn.solver._within_group_si` | `orpheus/sn/solver.py:778` | `type_uses` |
| `orpheus.sn.solver._solve_fixed_source_krylov` | `orpheus/sn/solver.py:3305` | `references` |
| `orpheus.sn.coupled_system.build_streaming_collision` | `orpheus/sn/coupled_system.py:356` | `type_uses` + `references` |
| `orpheus.sn.coupled_system.build_within_group_system` | `orpheus/sn/coupled_system.py:380` | `references` |
| `orpheus.numerics.green_operator._left_spine_terms` | `orpheus/numerics/green_operator.py:154` | `references` |
| `orpheus.numerics.operator.OperatorSum.is_invertible` | `orpheus/numerics/operator.py:1422` | `references` |
| `derivations.diagnostics.diag_sphere_fixedpoint_consistency.probe_real_si` | `derivations/diagnostics/diag_sphere_fixedpoint_consistency.py:113` | `references` |
| `derivations.diagnostics.diag_sphere_fixedpoint_consistency.probe_fixed_point_consistency` | `…:74` | `references` |
| `derivations.diagnostics.diag_curvilinear_seed_sensitivity._build` | `derivations/diagnostics/diag_curvilinear_seed_sensitivity.py:100` | `references` |

Depth-1 **doc nodes** (Nexus DOES see these — 12 of the 19 doc files appear as depth-1 `file`
nodes): `theory/foundations/boundary_conditions`, `theory/foundations/operator_algebra`,
`api/discrete_ordinates`, `theory/methods/sn/loss_representation`,
`theory/conventions/indexing_and_layout`, `theory/methods/sn/curvilinear_numerics`,
`theory/methods/sn/slab_multigroup`, `theory/methods/sn/slab_one_group`,
`theory/foundations/operator_inverse_family`, `theory/methods/sn/solver`,
`theory/foundations/structured_geometry`. (7 further doc files with hits do **not** appear at
depth 1 — see Surface 2; the graph's doc coverage is partial.)

### WHY THIS IS NOT SUFFICIENT (the three documented blind spots, instantiated here)

1. **Property-reached leaves.** `SweepOperator.is_adjointable` is a `@property`. Its body's
   `isinstance(self.inner, InvertibleOperator)` (`sweep_operator.py:154`) has *zero* `calls`
   edges to `InvertibleOperator` — the graph only records it as a class-level `references`
   edge. A `callers()` on any *method* of `InvertibleOperator` would miss it entirely.
2. **Class-name bypass consumers.** `tests/sn/verification/analytical/test_prescribed_inflow_consistency.py:242`
   asserts `type(rG).__name__ == "ScheduledInvertibleOperator"` — a **string** comparison. The
   graph has no edge for it. (See Surface 4.)
3. **Direct constructors of a guarded type.** `InvertibleOperator.__init__` validates operand
   types AND the mesh-identity invariant (`streaming.py:564-600`). Every `InvertibleOperator(...)`
   call is a guard entry point, and 13 of the 14 live in tests/derivations — invisible in a
   production-only graph read. (See Surface 3.)
4. **The doc gap.** Nexus surfaced 12 doc files; grep found **19**. The 7 the graph missed:
   `theory/methods/sn/history`, `theory/methods/sn/cartesian_multid`, `theory/methods/sn/adjoint`,
   `theory/methods/sn/acceleration`, `theory/foundations/operator_tensor_network`,
   `theory/foundations/cross_section_data`, `theory/foundations/coupled_block_operator`,
   `theory/verification/sn`.

---

## SURFACE 2 — Text grep across code, tests, and `docs/`

### 2a. Production code — LIVE CODE (executable; must change or the rename breaks)

| # | file:line | Kind |
|---|---|---|
| 1 | `orpheus/sn/operators/streaming.py:445` | **the class definition** |
| 2 | `orpheus/sn/operators/streaming.py:409` | `@overload` return annotation `-> "InvertibleOperator"` |
| 3 | `orpheus/sn/operators/streaming.py:436` | `return InvertibleOperator(self, other)` — **the sole production constructor** |
| 4 | `orpheus/sn/operators/scheduled_invertible.py:71` | **runtime** `from .streaming import InvertibleOperator` (module level, NOT `TYPE_CHECKING`) |
| 5 | `orpheus/sn/operators/scheduled_invertible.py:90` | generic type argument in the `OperatorSum[...]` base |
| 6 | `orpheus/sn/operators/scheduled_invertible.py:107` | `invertible: "InvertibleOperator"` param annotation |
| 7 | `orpheus/sn/operators/scheduled_invertible.py:110` | `isinstance(invertible, InvertibleOperator)` **runtime gate** |
| 8 | `orpheus/sn/operators/scheduled_invertible.py:135` | `def invertible(self) -> "InvertibleOperator"` |
| 9 | `orpheus/sn/operators/sweep_operator.py:73` | `TYPE_CHECKING` import |
| 10 | `orpheus/sn/operators/sweep_operator.py:80` | `SweepInvertible = Union[InvertibleOperator, ScheduledInvertibleOperator]` type alias |
| 11 | `orpheus/sn/operators/sweep_operator.py:151` | **late runtime import** inside `is_adjointable` |
| 12 | `orpheus/sn/operators/sweep_operator.py:154` | `isinstance(self.inner, InvertibleOperator)` **runtime gate** |
| 13 | `orpheus/sn/operators/sweep_operator.py:183` | **late runtime import** inside `apply_transpose` |
| 14 | `orpheus/sn/operators/sweep_operator.py:185` | `if not isinstance(self.inner, InvertibleOperator)` **runtime gate** |
| 15 | `orpheus/sn/solver.py:90` | `TYPE_CHECKING` import `from .operators.streaming import InvertibleOperator` |
| 16 | `orpheus/sn/solver.py:708` | `LC: "InvertibleOperator"` param annotation (`_select_si_resolvent`) |
| 17 | `orpheus/sn/solver.py:711` | return annotation `"tuple[InvertibleOperator \| ScheduledInvertibleOperator, …]"` |
| 18 | `orpheus/sn/solver.py:782` | return annotation `"tuple[…, CoupledOperator \| InvertibleOperator \| ScheduledInvertibleOperator, …]"` |
| 19 | `orpheus/sn/coupled_system.py:152` | `TYPE_CHECKING` import |
| 20 | `orpheus/sn/coupled_system.py:327` | dataclass field annotation `resolvent: "CoupledOperator \| InvertibleOperator"` |
| 21 | `orpheus/sn/coupled_system.py:358` | return annotation `-> "InvertibleOperator"` (`build_streaming_collision`) |

### 2b. Production code — ERROR-MESSAGE STRING LITERALS (user-visible; tests may `match=` them)

| file:line | String |
|---|---|
| `orpheus/sn/operators/streaming.py:571` | `f"InvertibleOperator: 'streaming' must be a …"` |
| `orpheus/sn/operators/streaming.py:576` | `f"InvertibleOperator: 'diagonal' must be a …"` |
| `orpheus/sn/operators/streaming.py:587` | `"InvertibleOperator: the streaming operator and the diagonal …"` |
| `orpheus/sn/operators/streaming.py:596` | `f"InvertibleOperator: diagonal coefficient must be …"` |
| `orpheus/sn/operators/streaming.py:677` | `_require_typed_composite("InvertibleOperator.apply", …)` — **guard tag string** |
| `orpheus/sn/operators/streaming.py:701` | `_require_typed_composite("InvertibleOperator.apply_transpose", …)` — **guard tag string** |
| `orpheus/sn/operators/streaming.py:897` | `f"InvertibleOperator: 'rhs' must be FullField; …"` |
| `orpheus/sn/operators/streaming.py:904` | `"InvertibleOperator.solve(FullField): rhs and …"` |
| `orpheus/sn/operators/streaming.py:1063` | `"InvertibleOperator.solve_transpose(FullField): b and …"` |
| `orpheus/sn/operators/scheduled_invertible.py:112-113` | `f"ScheduledInvertibleOperator: 'invertible' must be an InvertibleOperator; got {type(invertible).__name__}."` |
| `orpheus/sn/operators/sweep_operator.py:189` | `f"…the (L+C) InvertibleOperator arm; …"` (in `MissingAdjoint`) |
| `orpheus/sn/solver.py:934` | `f"InvertibleOperator (L + C) with the sweep as …"` |

### 2c. Production code — DOCSTRING / COMMENT ONLY (Sphinx-rendered; NOT gated, see 2d)

`orpheus/transport/operators/__init__.py:11` ·
`orpheus/transport/operators/multiplication_operator.py:53` ·
`orpheus/numerics/iteration.py:68` ·
`orpheus/numerics/operator.py:125,128,1294,1441,1460` ·
`orpheus/numerics/green_operator.py:68,162` ·
`orpheus/sn/operators/sweep_operator.py:2,23,29,37,39,43,87,88,110,134,135,146,164,179` ·
`orpheus/sn/solver.py:21,728,894,1036,1694,3327` ·
`orpheus/sn/coupled_system.py:95,309,315,364,400,529` ·
`orpheus/sn/operators/__init__.py:9` ·
`orpheus/sn/operators/radial_characteristic.py:97,228,529` ·
`orpheus/sn/operators/scheduled_invertible.py:37,44,48,53,101,121,171,217,246` ·
`orpheus/sn/operators/streaming.py:13,111,126,392,394,422,425,441,489,491,710,724,738,759` ·
`orpheus/sn/loss_representation/__init__.py:27,3070`

Note `coupled_system.py:95,315,529` and `tests/sn/_test_helpers.py:692` reference the **already
retired** `CoupledInvertibleOperator` (deleted at #280 step 5d) — these are historical prose,
not live references, and a naive `sed` rename would corrupt them.

### 2d. `docs/**.rst` — THE UNGATED SURFACE (100 hits across 19 files)

**Sphinx gate status — VERIFIED FOR THIS TREE.** `docs/conf.py` contains **no `nitpicky`
setting** (grep for `nitpick` → zero hits), so `nitpicky` is `False` by default, and the
project's documented build gate is `sphinx-build -E -W` (`.claude/lessons.md:1127`,
`.claude/plans/sn_split_catalog.md:42`) — **no `-n`**. Consequence:

> **83 of the 100 doc hits are Python-domain roles (`:class:` / `:meth:` / `:attr:` /
> `:func:` / `:mod:`). Every one of them will silently render as PLAIN TEXT after the rename,
> with ZERO build warnings.** The remaining 17: 8 inline-literal (` `` ` `) only, 7
> `:math:`/`\mathrm{}`, 6 code-block/continuation lines (categories overlap).

There is exactly ONE `automodule` directive that pulls the class into the API tree:
`docs/api/discrete_ordinates.rst:98` → `.. automodule:: orpheus.sn.operators.streaming`.
Because the module is autodoc'd wholesale (not per-class `autoclass`), the API page needs **no
edit** — but every `:class:`~orpheus.sn.operators.streaming.InvertibleOperator`` target
elsewhere resolves against that autodoc entry and WILL dangle.

#### Complete doc site enumeration (file:line)

**`docs/theory/methods/sn/loss_representation.rst` — 28 hits (heaviest)**
`255, 256, 263, 281, 282, 284, 402, 421, 422, 424, 429, 431, 469, 565, 1916, 1918, 1919, 2145,
2194, 2199, 2248, 2292, 2553, 2816, 2834, 2873, 2881, 2882`
(`:429`/`:431` are inside a literal `code-block:: python` reproducing the class body;
`:263`/`:469`/`:565` are `:math:`\mathrm{InvertibleOperator}(…)`` — **math-mode prose**;
`:2145` is `ScheduledInvertibleOperator`-only.)

**`docs/theory/foundations/operator_algebra.rst` — 11 hits**
`714, 743, 1042, 1058, 1066, 1072, 1073, 1082, 1354, 1814, 1916`

**`docs/theory/methods/sn/history.rst` — 8 hits**
`135` (`CoupledInvertibleOperator`, retired), `227` (idem), `245`, `775` (`Scheduled…`),
`776`, `777`, `932`, `980`

**`docs/theory/methods/sn/slab_one_group.rst` — 7 hits**
`202, 207, 209, 216, 218, 223, 821`

**`docs/theory/methods/sn/cartesian_multid.rst` — 7 hits**
`3058` (`Scheduled…`), `3549`, `3550`, `3552` (`Scheduled…`), `3561`, `3812` (`Scheduled…`),
`3857` (`Scheduled…`)

**`docs/theory/methods/sn/curvilinear_numerics.rst` — 6 hits**
`85, 755, 791, 848, 1228, 1375`

**`docs/theory/foundations/boundary_conditions.rst` — 6 hits**
`517, 665, 735, 791, 1391, 3418`

**`docs/theory/methods/sn/curvilinear_one_group.rst` — 5 hits**
`2259, 2272, 2274, 2578, 3147`

**`docs/theory/foundations/operator_inverse_family.rst` — 4 hits**
`362, 382, 559, 586`

**`docs/theory/methods/sn/slab_multigroup.rst` — 3 hits**
`390, 876, 917`

**2-hit files**
`docs/theory/verification/sn.rst:4145, 4177` (both `Scheduled…`) ·
`docs/theory/methods/sn/solver.rst:16, 545` ·
`docs/theory/methods/sn/adjoint.rst:716, 1011` (`:1011` = `Scheduled…`) ·
`docs/theory/methods/sn/acceleration.rst:1143, 1171` (`:1171` = `Scheduled…`) ·
`docs/theory/foundations/operator_tensor_network.rst:943, 945` ·
`docs/theory/conventions/indexing_and_layout.rst:1079, 1454`

**1-hit files**
`docs/theory/foundations/structured_geometry.rst:452` ·
`docs/theory/foundations/cross_section_data.rst:278` (`Scheduled…`) ·
`docs/theory/foundations/coupled_block_operator.rst:613`

**Sibling-name hazard**: **11 doc lines mention ONLY `ScheduledInvertibleOperator`**, which is a
DIFFERENT class that is NOT being renamed. A naive global `sed s/InvertibleOperator/X/` corrupts
all 11 (plus every `Scheduled…` mention in the mixed files). Any rename script MUST use a
negative-lookbehind (`(?<!Scheduled)(?<!Coupled)InvertibleOperator`) or an ordered two-pass.

### 2e. Tests — 34 files, 148 hits

Full line-level inventory is in the raw grep; the load-bearing groupings:

| Category | Sites |
|---|---|
| **Module-level imports** | `tests/sn/solve/test_si_single_primitive_contract.py:47` · `tests/sn/operators/test_green_operator_sn.py:67` · `test_psi_half_coupling.py:101` · `test_operator_block_role.py:64` · `test_removal_form_matvec_sweep.py:110` · `test_inverse_operator_equivalence.py:22` · `test_capability_survival.py:41` · `test_invertible_operator.py:47` · `test_inverse_adjoint_coherence.py:96` · `tests/sn/_fixtures/wave_t_t4/_capture_pre_t4_snapshots.py:71` |
| **In-function imports** | `tests/sn/solve/test_sn_adjoint_certification.py:271` · `tests/sn/operators/test_streaming_operator.py:337` |
| **`isinstance` assertions** | `test_si_single_primitive_contract.py:172,173,177,178` · `test_green_operator_sn.py:217` · `test_psi_half_coupling.py:3265` · `test_operator_block_role.py:116` · `test_streaming_operator.py:343` · `test_inverse_operator_equivalence.py:48` · `test_capability_survival.py:167,366` · `test_invertible_operator.py:132,161` |
| **`monkeypatch.setattr(InvertibleOperator, …)`** | `test_sn_adjoint_certification.py:277` · `test_inverse_adjoint_coherence.py:287,288,298,299` · `test_psi_half_coupling.py:3322,3323,3336,3337` |
| **Direct constructors** | see Surface 3 |
| **Docstring/comment only** | remaining ~100 hits |
| **File named after the symbol** | `tests/sn/operators/test_invertible_operator.py` (**the file NAME itself**) |

**Test-file rename candidate**: `tests/sn/operators/test_invertible_operator.py` — the module
path is a first-class part of the retirement (a stale filename is a "which path is canonical?"
smell). Note it also contains `class TestInvertibleSolveBridgeRegression` (`:636`), referenced
by NAME from `tests/transport/test_multiplication_operator.py:414` and
`tests/sn/operators/test_removal_form_matvec_sweep.py:344` — **as prose in comments**, a
graph-invisible cross-file reference.

### 2f. `derivations/` — 5 hits, 2 files

| file:line | Kind |
|---|---|
| `derivations/diagnostics/diag_curvilinear_seed_sensitivity.py:34` | `from orpheus.sn.operators.streaming import InvertibleOperator, StreamingOperator` |
| `derivations/diagnostics/diag_curvilinear_seed_sensitivity.py:102` | **constructor** |
| `derivations/diagnostics/diag_sphere_fixedpoint_consistency.py:37` | import |
| `derivations/diagnostics/diag_sphere_fixedpoint_consistency.py:82` | **constructor** |
| `derivations/diagnostics/diag_sphere_fixedpoint_consistency.py:123` | **constructor** |

Both are the diagnostics that probed the `initial_guess` seed-sensitivity question (Section C) —
they are the historical instruments behind the #282 route-(a) decision. Per
`coding-standards.md`, a diagnostic consumed by a tracked test is production infrastructure;
neither of these has a tracked-test consumer (checked), but they are LIVE importers and will
`ImportError` on rename.

### 2g. `.claude/**` — 510 hits across ~105 files (NOT build-gated, NOT test-gated)

Heaviest (top 10): `plans/issue_226_inverse_operator_verification.md` (30) ·
`plans/archive/r1_step4_cooperative_carve.md` (29) · `plans/operator_machinery_taxonomy.md` (19) ·
`plans/operator_inverse_algebra_carve.md` (19) · `agent-memory/test-architect/s6_relayering_verification.md` (14) ·
`agent-memory/test-architect/wave_t_t4_streaming_verification_spec.md` (13) ·
`skills/vv-principles/error_catalog.md` (12) ·
`agent-memory/test-architect/coupled_operator_step5_solve_verification.md` (12) ·
`agent-memory/explorer/dh1b_angular_flux_consumer_audit.md` (12) ·
`plans/phase2_code_prose/classify_P2-G.md` (11).

**Recommendation**: do **NOT** mass-rewrite `.claude/plans/archive/` or agent memories — they
are dated archaeology and rewriting them falsifies the historical record. The ones that DO need
attention are the **live, standing** documents:
- `.claude/skills/vv-principles/error_catalog.md` (12 hits) — the ERR catalog is a standing
  reference consulted every session.
- `.claude/lessons.md` (6 hits, incl. **L21 at `:516`, `:548`**) — see Section C.
- `.claude/plans/sn_operator_realization_and_posing_plan.md` (6) — the plan driving this rename;
  it should be updated as part of the campaign, not before it.
- `.claude/agent-memory/*/lessons.md` (test-architect 9, qa 6, elegance-enforcer 5,
  numerics-investigator 1, method-implementer 1).

### 2h. `scratch/` — 14 hits, 1 file

`scratch/review_map_sn_assembly.md:18,53,101,108,137,142,165,173,174,187,254,262,266,708` — an
untracked working file from this same campaign. No action needed (untracked scratch).

---

## SURFACE 3 — Direct constructors, `isinstance` gates, and imports

### 3a. Direct constructors `InvertibleOperator(...)` — 14 sites total, ONE in production

| # | file:line | Context |
|---|---|---|
| **P1** | `orpheus/sn/operators/streaming.py:436` | `return InvertibleOperator(self, other)` inside `StreamingOperator.__add__` — **the only production construction; every `(L+C)` in the codebase flows through it** |
| D1 | `derivations/diagnostics/diag_curvilinear_seed_sensitivity.py:102` | diagnostic build |
| D2 | `derivations/diagnostics/diag_sphere_fixedpoint_consistency.py:82` | diagnostic build |
| D3 | `derivations/diagnostics/diag_sphere_fixedpoint_consistency.py:123` | diagnostic build |
| T1 | `tests/sn/operators/test_invertible_operator.py:146` | removal-form explicit construction |
| T2 | `tests/sn/operators/test_invertible_operator.py:171` | capability probe |
| T3 | `tests/sn/operators/test_invertible_operator.py:208` | **negative** — mesh-mismatch raises |
| T4 | `tests/sn/operators/test_invertible_operator.py:219` | **negative** — non-positive σ raises |
| T5 | `tests/sn/operators/test_invertible_operator.py:230` | **negative** — zero σ raises |
| T6 | `tests/sn/operators/test_invertible_operator.py:238` | **negative** — `InvertibleOperator(C, C)` wrong operand |
| T7 | `tests/sn/operators/test_invertible_operator.py:240` | **negative** — `InvertibleOperator(L, L)` wrong operand |
| T8 | `tests/sn/operators/test_removal_form_matvec_sweep.py:267` | removal form |
| T9 | `tests/sn/operators/test_removal_form_matvec_sweep.py:307` | removal form + `is_adjointable` |
| T10 | `tests/sn/operators/test_removal_form_matvec_sweep.py:373` | removal form |
| T11 | `tests/sn/operators/test_removal_form_matvec_sweep.py:433` | removal form |
| T12 | `tests/sn/operators/test_removal_form_matvec_sweep.py:539` | removal form |
| T13 | `tests/sn/operators/test_removal_form_matvec_sweep.py:586` | **negative** — `σ_r ≤ 0` raises |

**Guard blast radius**: the constructor guards live at `streaming.py:564-600` (four `raise`
paths, each with a message beginning `"InvertibleOperator: "` — `:571`, `:576`, `:587`, `:596`).
T3–T7 and T13 are the six negative tests that pin them.

**VERIFIED SAFE**: none of the six `pytest.raises(..., match=...)` patterns asserts on the class
name. They match on `"mesh-identity"` (`test_invertible_operator.py:207, 357, 629`),
`"strictly positive"` (`:218, :229`, `test_removal_form_matvec_sweep.py:587`),
`"StreamingOperator"` (`:237`), `"MultiplicationOperator"` (`:239`), and
`"expected FullField"` (`:346`). A tree-wide grep of `match=` lines containing `Invertible`
returns only `NotInvertible` (a different, untouched symbol). **The error-message rename in 2b
carries no test breakage** — it is a user-visible-text change only.

### 3b. `isinstance(..., InvertibleOperator)` runtime gates

**Production (3 — all structurally load-bearing):**

| file:line | Gate |
|---|---|
| `orpheus/sn/operators/sweep_operator.py:154` | `SweepOperator.is_adjointable` — decides whether `A.H.inverse()` is legal |
| `orpheus/sn/operators/sweep_operator.py:185` | `SweepOperator.apply_transpose` direct-call backstop (raises `MissingAdjoint`) |
| `orpheus/sn/operators/scheduled_invertible.py:110` | `ScheduledInvertibleOperator.__init__` operand-type guard |

**Tests (12):** `test_si_single_primitive_contract.py:172,173,177,178` ·
`test_green_operator_sn.py:217` · `test_psi_half_coupling.py:3265` ·
`test_operator_block_role.py:116` · `test_streaming_operator.py:343` ·
`test_inverse_operator_equivalence.py:48` · `test_capability_survival.py:167,366` ·
`test_invertible_operator.py:132,161`

### 3c. `import InvertibleOperator` sites

**Production runtime (module-level) — 1:**
`orpheus/sn/operators/scheduled_invertible.py:71` — `from .streaming import InvertibleOperator`

**Production runtime (LATE, inside function body — graph-visible but easy to miss on a
file-header-only sweep) — 2:**
`orpheus/sn/operators/sweep_operator.py:151` (inside `is_adjointable`) ·
`orpheus/sn/operators/sweep_operator.py:183` (inside `apply_transpose`)

**Production `TYPE_CHECKING` — 3:**
`orpheus/sn/solver.py:90` · `orpheus/sn/coupled_system.py:152` ·
`orpheus/sn/operators/sweep_operator.py:73`

**Derivations — 2:** `diag_curvilinear_seed_sensitivity.py:34` · `diag_sphere_fixedpoint_consistency.py:37`

**Tests — 12** (10 module-level + 2 in-function): listed in 2e.

**Total import sites: 20.**

---

## SURFACE 4 — The L34 hazard (segment-assembled / string-literal references)

| # | file:line | Hazard | Graph-visible? |
|---|---|---|---|
| H1 | `tests/sn/verification/analytical/test_prescribed_inflow_consistency.py:242` | `type(rG).__name__ == "ScheduledInvertibleOperator"` — **a string equality on the class name**. This is the sibling class, NOT the rename target, but it is the exact hazard shape and a negative-lookbehind failure here would corrupt a live assertion. | **NO** |
| H2 | `orpheus/sn/operators/scheduled_invertible.py:113` | `f"…InvertibleOperator; got {type(invertible).__name__}."` — the message names the target class as a literal AND interpolates a runtime `__name__`. | NO |
| H3 | `orpheus/sn/operators/streaming.py:677` | `_require_typed_composite("InvertibleOperator.apply", …)` — the **method tag** is a bare string passed into the shared guard; it appears in every raised error. | NO |
| H4 | `orpheus/sn/operators/streaming.py:701` | `_require_typed_composite("InvertibleOperator.apply_transpose", …)` — same. | NO |
| H5 | `orpheus/sn/operators/sweep_operator.py:187-189` | `MissingAdjoint(f"SweepOperator over {type(self.inner).__name__} has no … the (L+C) InvertibleOperator arm …")` — runtime `__name__` + literal. | NO |
| H6 | `docs/theory/methods/sn/loss_representation.rst:429,431` | A `code-block:: python` **reproducing the class source** (`class InvertibleOperator(OperatorSum):` and the `_require_typed_composite("InvertibleOperator.apply", …)` line). Doc code-blocks are never gated and never type-checked. | NO |
| H7 | `tests/transport/test_multiplication_operator.py:414` · `tests/sn/operators/test_removal_form_matvec_sweep.py:344` | Cross-file references to `TestInvertibleSolveBridgeRegression` **by name, in a comment** (the class lives in `test_invertible_operator.py:636`). | NO |
| H8 | `docs/api/discrete_ordinates.rst:98` | `.. automodule:: orpheus.sn.operators.streaming` — a **module path assembled as a directive argument**. Safe for a class rename, but if the rename campaign also moves the module, this is the L34 path-string catch. | NO |
| H9 | `.claude/agent-memory/test-architect/s6_relayering_verification.md:290` | `def _representation_used_by_solve(A: "InvertibleOperator")` — a **code snippet inside a memory file**, i.e. a spec fragment a future agent may copy verbatim. | NO |
| H10 | `.claude/plans/sn_split_catalog.md:243` | `"InvertibleOperator"` as a quoted **doc-split catalog key** naming a section of `loss_representation.rst`. | NO |

**NEGATIVE results (searched, nothing found — good news):**
- **`__all__`**: no `__all__` anywhere in the tree contains `InvertibleOperator` (see Section A).
- **`getattr(..., "InvertibleOperator")`**: zero hits.
- **`textwrap.dedent` subprocess-worker imports**: zero files contain both `textwrap` and
  `Invertible`.
- **Registry / dispatch-table string keys**: zero.
- **Bare-segment `Invertible` used as an independent token**: all remaining `Invertible`
  occurrences resolve to `NotInvertible` (the exception class), `is_invertible` (the predicate),
  `_InvertibleForward` (a numerics Protocol), `SweepInvertible` (the Union alias),
  `ScheduledInvertibleOperator`, `CoupledInvertibleOperator` (retired), or
  `TestInvertibleSolveBridgeRegression` / `TestInvertible…` test class names. **None of these
  should be touched by the rename.** Any rename regex MUST be anchored to the exact token
  `InvertibleOperator` with a negative lookbehind on `Scheduled` / `Coupled`, and MUST NOT match
  `NotInvertible` / `_InvertibleForward` / `SweepInvertible` / `is_invertible`.

---

## A. Is `InvertibleOperator` in any `__all__`?

**NO — confirmed, tree-wide.**

| Check | Result |
|---|---|
| `orpheus/sn/operators/streaming.py:114-116` | `__all__ = ["StreamingOperator"]` — **`InvertibleOperator` is NOT exported** (verified by reading the file) |
| `orpheus/sn/operators/scheduled_invertible.py:84` | `__all__ = ["ScheduledInvertibleOperator"]` (the sibling IS exported) |
| `orpheus/sn/operators/__init__.py` | **has no `__all__` at all** — the file is a pure module docstring (32 lines, no code). It *documents* `InvertibleOperator` in prose at `:9` but re-exports nothing. |
| Tree-wide `grep -rn "__all__" orpheus/ -A6 \| grep -i invertible` | Only the `scheduled_invertible.py:84` hit, plus incidental context lines |
| `orpheus/numerics/operator.py:163` | `"NotInvertible"` — the **exception**, a different symbol |

**Consequence**: there is no public export path to break. Every consumer imports the concrete
module path `orpheus.sn.operators.streaming`, which makes the rename mechanically simpler — but
it also means **no `__init__` re-export shim exists to soften the transition**, so all 20 import
sites must land in the same commit.

---

## B. Subclasses / MRO consumers

**There are ZERO subclasses of `InvertibleOperator` anywhere in the tree.** Nexus reports **no
`inherits` edge** into the class (verified: `neighbors(direction=in, edge_types=inherits,…)`
returned no `inherits` entries — only `references` / `type_uses` / `contains`).

### The `ScheduledInvertibleOperator` question — it COMPOSES, it does NOT inherit

`orpheus/sn/operators/scheduled_invertible.py:87-93`:

```python
class ScheduledInvertibleOperator(
    OperatorSum[
        "FullField", "FullField",
        InvertibleOperator,
        ScaledOperator["FullField", "FullField", SNMaskedBoundaryOperator],
    ],
):
```

`InvertibleOperator` appears as a **generic type ARGUMENT to `OperatorSum`**, not as a base
class. The MRO is `ScheduledInvertibleOperator → OperatorSum → …`; `InvertibleOperator` is **not
in it**. Structurally it is a *has-a*: the instance holds one as `self.a`, exposed through the
`invertible` property (`:135-137`) and validated by `isinstance` at `:110`. The two classes are
deliberately **twinned siblings** (both subclass `OperatorSum` directly and both re-implement
the `is_invertible`/`inverse`/`solve` back-half) — `streaming.py:756-764` documents this as a
defer-until-≥2 decision with an explicit "extract a shared mixin at the 3rd witness" TRIGGER.

### The actual inheritance edges

| Class | file:line | Base(s) |
|---|---|---|
| `InvertibleOperator` | `orpheus/sn/operators/streaming.py:445-449` | `OperatorSum["FullField","FullField","StreamingOperator","MultiplicationOperator"]` |
| `ScheduledInvertibleOperator` | `orpheus/sn/operators/scheduled_invertible.py:87-93` | `OperatorSum[…]` — **sibling, not child** |
| `SweepOperator` | `orpheus/sn/operators/sweep_operator.py:83-85` | `InverseWrapMixin["SweepInvertible"], LinearOperator["FullField"]` — **wraps** an `InvertibleOperator` as `self.inner` |

### MRO-*behaviour* consumers (rename-safe, but worth knowing)

`InvertibleOperator` **overrides** inherited `OperatorSum` behaviour in four places — the rename
does not change dispatch, but any doc prose describing "the subclass shadows the generic Green
by MRO" (`streaming.py:465-467`) must be re-worded consistently:
`apply` (`:659`) · `apply_transpose` (`:680`) · `inverse` (`:745`, shadows `OperatorSum.inverse`'s
generic `GreenOperator`) · `is_invertible` (`:609`).

Also note `orpheus/numerics/operator.py:1422` `OperatorSum.is_invertible` carries a **docstring
`references` edge** to `InvertibleOperator` — the base class documents its own subclass by name.

---

## C. The `initial_guess` question

### C(0) — VERIFIED: `del initial_guess` is real, and there are TWO of them

`orpheus/sn/operators/streaming.py:770-832`:

```python
def solve(
    self,
    rhs: "FullField",
    *,
    initial_guess: "FullField | None" = None,
) -> "TimedFullField":
    ...
    del initial_guess  # accept-and-drop: exact direct inverse, nothing to seed (#280 2.5c)
    return self._solve_timed_full_field(rhs)
```

Confirmed at **`streaming.py:831`** exactly as the brief stated. The private body
`_solve_timed_full_field` (`:841-848`) has signature
`(rhs, *, moment_frame=None, schedule=None, reflect=None)` — **no seed parameter exists on it at
all**. The docstring at `:813-821` states the contract explicitly: *"ACCEPTED and DROPPED — the
WDD sweep is an EXACT direct inverse with nothing to seed … a warm start lives at the iteration
layer … never here."*

The **inverse-operator wrapper does the same drop one level up**:
`orpheus/sn/operators/sweep_operator.py:126` — `del initial_guess  # direct exact inverse —
nothing to seed (#282/2.5d)` before `return self.inner.solve(rhs)`. So the seed is dropped
**twice** on the production path, and the second drop means `InvertibleOperator.solve` is not
even *reached* with a non-None seed from production.

The sibling `ScheduledInvertibleOperator.solve` drops it identically at
`scheduled_invertible.py:207`.

### C(i) — Which call sites pass `initial_guess` a NON-None value?

#### PRODUCTION: **ZERO sites pass a seed to `InvertibleOperator.solve` or `SweepOperator.apply`.**

The one production site that threads a non-None seed toward the sweep is the SI driver:

`orpheus/numerics/iteration.py:722` — `psi = self.A_inv.apply(rhs, initial_guess=psi_prev)`

…and `A_inv` for SN is built at `orpheus/sn/solver.py:851` / `:866` / `:873` as
`system.resolvent.inverse()` — i.e. a `SweepOperator` (or `WindowedSweep` /
`CoupledSubstitutionOperator`), **every one of which drops the kwarg**:

| Receiver | file:line | Behaviour |
|---|---|---|
| `SweepOperator.apply` | `orpheus/sn/operators/sweep_operator.py:126` | `del initial_guess` |
| `WindowedSweep.apply` | `orpheus/sn/operators/windowing.py:238` | `del initial_guess  # no bulk-seed consumer in the multi-D walk` |
| `CoupledSubstitutionOperator.apply` | `orpheus/numerics/coupled_system.py:1268` | `del initial_guess  # exact direct substitution` |
| `MatrixInverseOperator.apply` | `orpheus/numerics/matrix_inverse_operator.py:223` | `del initial_guess` |
| `InverseOperator.apply` (numerics mixin) | `orpheus/numerics/operator.py:2164` | `del initial_guess` |
| `_SeededExactApply` | `orpheus/numerics/iteration.py:222` | `del initial_guess` |
| **`GreenOperator.apply`** | `orpheus/numerics/green_operator.py:320-323` | **CONSUMES it** — `psi = initial_guess`, then `self._driver.solve(q, initial_guess=psi)` |

So the ONLY production consumer of a threaded `initial_guess` is the **iterative**
`GreenOperator` (a splitting iteration, which genuinely has a start to seed), and the entry
points `SourceIteration.solve(rhs, initial_guess=x0)` / `KrylovAcceleration.solve(...)` /
`KEigenvalue.solve(initial_guess=...)` — all at the **iteration layer**, called from
`orpheus/sn/solver.py:1623, 1806, 2485, 2636, 3207, 3271, 3418`. **None of these is the
operator's `solve`.**

**Verdict for production**: the parameter on `InvertibleOperator.solve` is currently
**inert by design** — a uniform-protocol placeholder (`SupportsSeededApply`,
`orpheus/numerics/iteration.py:194`), not a live seam.

#### TESTS / DERIVATIONS: sites that DO pass a non-None seed to the operator's `solve`/`inverse().apply`

| # | file:line | Call | Purpose |
|---|---|---|---|
| 1 | `tests/sn/operators/test_inverse_operator_equivalence.py:85` | `A.inverse().apply(b, initial_guess=seed)` | **seed-INDEPENDENCE gate** — asserts result identical with/without seed |
| 2 | `tests/sn/operators/test_inverse_operator_equivalence.py:103-108` | strict spy: `def spy(rhs)` (no kwarg) then `A.inverse().apply(b, initial_guess=seed)` | **regression pin** — apply must NOT forward the seed inward; a resurrected thread reddens this |
| 3 | `tests/sn/operators/test_operators_apply_typed.py:480` | `assert_type(sweep.apply(rhs, initial_guess=seed), TimedFullField)` | static-typing assertion only |
| 4 | `tests/sn/operators/test_invertible_operator.py:845` | `LC.solve(LC_psi, initial_guess=psi_recovered)` | fixed-point loop (result unaffected — the kwarg is dropped) |
| 5 | `tests/sn/operators/test_invertible_operator.py:967` | `LC.solve(rhs_a, initial_guess=guess)` | hand-rolled SI loop inside `_joint_solve` |
| 6 | `tests/sn/sweep/test_cyl_direct_seed_fold.py:247` | `A.solve(q, initial_guess=psi)` | direct-seed fold vs `np.linalg.solve(M,q)` |
| 7 | `tests/sn/sweep/curvilinear/test_282_direct_seed_fixed_point.py:117,122,195,196` | `A.solve(b, initial_guess=…)` ×2 with **random** seeds | **the #282 seed-independence gate** — "the lag SIGNATURE dies: two random seeds give the same answer" |
| 8 | `derivations/diagnostics/diag_curvilinear_seed_sensitivity.py:123,124,185` | `op.solve(b, initial_guess=X1/X2)` | the original diagnostic that MEASURED seed-dependence |
| 9 | `derivations/diagnostics/diag_sphere_fixedpoint_consistency.py:98,149` | `op.solve(b, initial_guess=psi0)` | fixed-point-consistency probe |

Items 1, 2, and 7 are **the parameter's actual reason for existing today**: they are the gates
that prove the sweep is seed-INDEPENDENT. Removing the kwarg would delete the vocabulary those
gates are written in.

### C(ii) — Is the L21 seam still live, or superseded? **SUPERSEDED.**

L21 (`.claude/lessons.md:511-570`) established: *sweep and matvec are two applications of ONE
operator; the M-M Carlson `psi_half_seed` strategy is a property of `L`, so the matvec passes
`psi_view` (the apply target) and the sweep passes `initial_guess` (the previous iterate) into
the **same** `psi_half_seed(psi_level=…)`.* The kwarg was minted at `:547-549` precisely to
carry the sweep's half of that pair.

**That half is GONE. Traced in the current tree:**

1. `orpheus/sn/sweep/pole_angular_closure.py:1032-1050` —
   `MorelMontryAngularSweep.precompute_psi_state(psi_view, *, radial_characteristic=None)`.
   The docstring's "Seed dispatch (#282 route (a), R12a)" says verbatim:
   > *"the recurrence seed is the composite's FIRST-CLASS ψ½ state:
   > `radial_characteristic.cells(p, -1)` … **The retired-strategy extrapolation of the ITERATE
   > — the #282 back edge — is gone: the seed is upstream STATE in the augmented walk order.**"*

2. Non-carrying levels seed from `self.edge_extrapolated_seed(psi_level, p)`
   (`pole_angular_closure.py:1079`, defined `:982`) — derived from the **current** angular level,
   not from any previous iterate.

3. `precompute_psi_state` has exactly **4 call sites**, and **none is the sweep's seed channel**:
   `orpheus/sn/operators/radial_characteristic.py:784` ·
   `orpheus/sn/loss_representation/__init__.py:3146` (the matvec walk) ·
   `orpheus/sn/loss_representation/__init__.py:3541` (a coefficient probe) ·
   `tests/sn/sweep/curvilinear/test_unified_matvec_cylinder.py:158`.

4. `InvertibleOperator._solve_timed_full_field` has **no seed parameter to forward** — the
   channel is structurally absent, not merely unused.

**What happens to a warm start today** (the honest current answer):

```
SNSolver._solve_source_iteration (solver.py:1590-1623)
  initial_guess = self._psi_typed  (or a cold start: _windowed_cold_start / zeros)
    → SourceIteration.solve(q_driver, initial_guess=initial_guess)   [ITERATION layer — CONSUMED]
        iteration.py:664-672  →  psi = initial_guess  (the loop's starting iterate)
        iteration.py:722      →  self.A_inv.apply(rhs, initial_guess=psi_prev)
            → SweepOperator.apply             sweep_operator.py:126   del initial_guess
                → InvertibleOperator.solve    streaming.py:831        del initial_guess
                    → _solve_timed_full_field(rhs)                    (no seed parameter)
```

The warm start is **entirely an iteration-layer concept**: `SourceIteration` holds `psi_prev`
and feeds it into `rhs = q + Σ gains.apply(psi_prev)` — the seed enters through the **right-hand
side**, never through the sweep's interior recurrence. The one operator that still consumes a
seed is the *iterative* `GreenOperator` (`green_operator.py:320-323`), which is a splitting
iteration and therefore genuinely has a start.

### C(iii) — VERDICT: not dead, not a regressed seam — a **deliberate uniform-protocol contract**

| Question | Answer |
|---|---|
| Is the parameter *used* by the body? | **No** — dropped at `streaming.py:831`. |
| Is the L21 seam live? | **No** — superseded by #282 route (a) (direct ψ½ seed from state) and #280 2.5c/2.5d. |
| Is it a *regression* (a seam that silently rotted)? | **No** — the supersession is explicit and documented at `streaming.py:813-821`, `sweep_operator.py:111-124`, `pole_angular_closure.py:1038-1050`, and issue #285. |
| Is it *removable*? | **NO, not cleanly.** It is a **Protocol conformance requirement**: `SupportsSeededApply.apply(self, rhs, /, *, initial_guess=None)` (`orpheus/numerics/iteration.py:194`) is what `SourceIteration` type-checks its `A_inv` against (`iteration.py:591-606`), and `iteration.py:719-721` threads it "UNCONDITIONALLY (pinned by the seed-threading spy)". Six other exact inverses accept-and-drop it identically. Removing it from `InvertibleOperator.solve` alone would break the uniform family; removing it family-wide would delete the vocabulary of the seed-independence gates (C(i) items 1, 2, 7). |

**Recommendation for the rename campaign**: leave `initial_guess` exactly as-is. It is orthogonal
to the rename. It is worth a separate issue only if the user wants to re-litigate whether the
`SupportsSeededApply` uniform contract should split into seeded vs. exact faces — a #285-scope
question, not a naming one.

**One genuine defect found in passing (STALE PROSE, worth fixing in the same PR):** two test
docstrings still assert the retired L21 seam as current —

- `tests/sn/operators/test_invertible_operator.py:677-680`: *"Post-Phase-1.2 the M-M Carlson seed
  travels through the explicit `initial_guess` kwarg on `InvertibleOperator.solve` (not through
  `rhs(1)` history); the iterative tests below thread the previous iterate that way."*
- `tests/sn/operators/test_invertible_operator.py:943-944`: *"`initial_guess` still threads the
  M-M Carlson bulk warm start."*

Both are **false against the current body**. (The same file's `:591-594` comment is correct and
contradicts them: *"the vestigial `initial_guess` retired with the dead seed threading, #280
2.5c."*) Also `.claude/lessons.md:547-549` (L21's own text) states the kwarg's live purpose and
should carry a supersession note.

---

## D. Name-collision check for candidate new names

Searched the whole tree (excluding `.git`, `docs/_build`, and the
`.claude/worktrees/nexus-workspace-wiring/` sibling checkout).

| Candidate | Collisions | Verdict |
|---|---|---|
| **`SNFreeTransportOperator`** | `.claude/plans/sn_operator_realization_and_posing_plan.md:506, 622, 623, 948, 1086` — **this plan already SELECTS this name** (`:623` "This plan selects `SNFreeTransportOperator`"; `:948` the rename table row `InvertibleOperator → SNFreeTransportOperator`; `:1086` the rename step). No code collision. | **CLEAR — and already the plan of record.** Only hits are the plan proposing it. |
| **`FreeTransportOperator`** | `.claude/plans/sn_operator_realization_plan_REVIEW.md:160` (the review discussing it). Zero code/doc/test hits. Loose grep `FreeTransport\|free_transport`: only `.claude/plans/sn_operator_realization_and_posing_plan.md:506` (`free_transport:` field name). | **CLEAR.** |
| **`AttenuatedStreamingOperator`** | `.claude/plans/sn_operator_realization_plan_REVIEW.md:160` only. | **CLEAR.** |
| **`StreamingCollisionOperator`** | `.claude/plans/sn_operator_realization_and_posing_plan.md:623` (as `SNStreamingCollisionOperator`, the rejected alternative) · `.claude/agent-memory/general-purpose/phase_g_four_operator_architecture_reconciliation.md:859` mentions a hypothetical `FusedStreamingCollisionOperator`. Zero code hits. | **CLEAR**, but see the SEMANTIC clash below. |
| **`SNFreeTransportOperator`** (duplicate row in brief) | as above | — |
| **`RemovalOperator`** | **ZERO hits, tree-wide** — the exact token appears nowhere. Loose grep for `Removal` finds only English prose ("removal form", "removal cross-section", `test_removal_form_matvec_sweep.py`, DSA literature OCR) and no symbol. | **CLEAR as a token**, but see the SEMANTIC clash below. |

### Semantic collisions (not token collisions — flagged because they matter more)

1. **`StreamingCollisionOperator` vs the existing FACTORY `build_streaming_collision`.**
   `orpheus/sn/coupled_system.py:356` is already named `build_streaming_collision` and is
   documented as *"THE one LC spelling"* (`:359-372`, a single-source-of-truth ruling). Naming
   the class `StreamingCollisionOperator` would make the factory read
   `build_streaming_collision() -> StreamingCollisionOperator`, which is maximally consistent —
   **this is a POSITIVE, not a hazard.** Note the flip side: choosing `SNFreeTransportOperator`
   instead leaves the factory name out of step and should trigger a matching factory rename
   (`build_free_transport`), or the greppability rule
   (`[[feedback-naming-consistency-greppable]]`) is violated.

2. **`RemovalOperator` collides with an established, DIFFERENT domain meaning.**
   In this codebase "removal" specifically means `σ_r = σ_t − Σ_{s,0}^{g→g}` — the
   *removal-form* variant `InvertibleOperator(L(σ_t), C(σ_r))` where within-group self-scatter
   is folded into the diagonal (`streaming.py:489-496`; `tests/sn/operators/test_removal_form_matvec_sweep.py`;
   `docs/theory/methods/sn/loss_representation.rst:263, 469, 565`). Naming the *general* class
   `RemovalOperator` would assert that every `(L+C)` is a removal-form composite, which is false
   — the default construction is `C = M[σ_t]`, not `M[σ_r]`. **RECOMMEND AGAINST.**

3. **`SN` prefix consistency.** The `sn/operators/` package already carries the SN prefix on
   `SNBoundaryOperator`, `SNMaskedBoundaryOperator`, `SNBoundaryRealizer` — but NOT on
   `StreamingOperator`, `SweepOperator`, `ScheduledInvertibleOperator`,
   `RadialCharacteristicOperator`. So `SNFreeTransportOperator` would be prefix-consistent with
   the boundary family and prefix-INCONSISTENT with the streaming family it actually belongs to.
   Worth a deliberate ruling under `[[feedback-naming-consistency-greppable]]` (one word-order
   pattern for the whole family) rather than a per-symbol choice.

---

## Rename execution checklist (derived from the four surfaces)

| Step | Scope | Count |
|---|---|---|
| 1 | Class definition + the sole production constructor + `@overload` | `streaming.py:409, 436, 445` |
| 2 | Production runtime imports (1 module-level + 2 late) | 3 |
| 3 | Production `TYPE_CHECKING` imports | 3 |
| 4 | Production `isinstance` gates | 3 |
| 5 | Production string annotations (params / returns / dataclass field / Union alias) | 7 |
| 6 | Production error-message string literals (**incl. the 2 `_require_typed_composite` tag strings**). No test `match=` depends on them — verified. | 12 |
| 7 | Production docstrings / comments — check every `:class:`/`:meth:`/`:attr:` role | ~60 |
| 8 | Tests: imports (12), `isinstance` (12), constructors (13), `monkeypatch.setattr` (9), docstrings (~100) | 148 hits / 34 files |
| 9 | Test **file** rename `tests/sn/operators/test_invertible_operator.py` + the 2 by-name comment refs to `TestInvertibleSolveBridgeRegression` | 3 |
| 10 | `derivations/diagnostics/` imports + constructors | 5 hits / 2 files |
| 11 | **`docs/**.rst` — 100 hits / 19 files. 83 are UNGATED Python-domain roles.** Verify with `sphinx-build -E -W -n` (add `-n` explicitly; `conf.py` has no `nitpicky`) or a post-rename grep for the OLD token. | **100** |
| 12 | `.claude/` live documents only (error_catalog, lessons.md L21 supersession note, agent lessons, the driving plan) — NOT the archive | ~35 of 510 |
| 13 | Stale-prose fix: `test_invertible_operator.py:677-680, 943-944` claim a retired seam | 2 |
| 14 | Sibling-protection: negative lookbehind on `Scheduled`/`Coupled`; do NOT touch `NotInvertible`, `is_invertible`, `_InvertibleForward`, `SweepInvertible` | — |

**The one gate that will not protect you**: `sphinx-build -W`. Add `-n`, or grep for the old
token after the rename. That is the whole point of this audit.
