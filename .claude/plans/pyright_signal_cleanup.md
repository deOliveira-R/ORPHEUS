# Pyright signal cleanup (#226) — "only signal, not noise"

**Goal (user, 2026-06-21):** make pyright give the agent *only real, actionable
signal* — (1) **kill the noise**: the streamed `<new-diagnostics>` LSP avalanche
of phantom `reportMissingImports` / "not defined" that floods every turn must
stop or match the CLI; (2) **drive real errors toward zero** so a genuine pyright
hit is always a real bug, never background hum.

**Oracle:** `npx pyright` (CLI). The streamed LSP `<new-diagnostics>` is advisory
and currently NOISY (microsoft/pyright#10498 dual-identity). **Constraint (hard,
user):** NO `# type: ignore`, NO blanket suppression — fixes are principled
(declare dynamically-set attrs statically [runtime-inert]; fix the real type;
correctly parametrise a Protocol). Per-cluster commit + **re-measure** (counts are
coupled — one `Unknown`-typed root cascades into argtype+attribute at every call
site, so raw counts overstate the independent-edit count).

**Branch:** fresh off `main` @ `7aebab3` (e.g. `refactor/pyright-signal-2`; note
an older `refactor/pyright-signal` branch predates the #257 campaign — do NOT
reuse it without rebasing onto current main). ff-merge per cluster-batch.

---

## Progress log

### Session 2026-06-22 — orpheus/ **502 → 470**, all merged + pushed to `main` (`5856160 → 11fa279`)

Ratchet re-baselined 706 (stale) → taut, then per-cluster commits (each: tests green under `-O`, count drops, 0 net-new `# type: ignore`, ratchet re-tightened):

- **Ratchet re-baseline** (`af13734`): the #257 + foundation merges had burned orpheus/ 706→502 without re-tightening — closed a 204-error regression hole. Baseline now tracks the live per-module floor (`tests/_harness/pyright_baseline.json`).
- **B1 warm-up** (`fe90b77`, −16): the 11 `reportUndefinedVariable` were REAL missing-import annotations in `derivations/` (sympy stays lazy → a module-level `TYPE_CHECKING` block; `Quadrature1D` folded into the already-runtime `common.quadrature` import; `TimedFullField`) — fixed regardless of the B6 derivations exec-env decision. Cluster D = 6 `bool_`/`csr`/`csc` return narrowings (`bool(...)`, build `csr_matrix` from triplets directly, `np.asarray(spsolve(...))`).
- **B2 §3 MoC plotters** (`3ed5fec`, −6): `plot_moc_rays`/`plot_moc_mesh` were a retired-Cartesian relic (`geom.n_cells`/`delta`/`mat_map`) that `AttributeError`'d **live** (demo_moc.py calls them). Rewritten against the actual reverse-Wigner-Seitz ring model (`pitch = R_outer·√π` equal-area square, concentric `radii`, real `Track.entry/exit` chords). User reviewed the rendered figures. **#266 filed** for the third plotter `plot_moc_spatial_flux` + the broadly-stale `demo_moc.py` (`default_pwr()` doesn't exist either) — a bigger rewrite + keff re-validation; it is the **last orpheus/(root) error**.
- **B2 §2 AngularMeasure** (`aa93757`, −4): widened the deliberately-lean geometry `AngularMeasure` Protocol with the genuine curvilinear-contract members `eta`/`mu_z`/`level_indices` (types matching the concrete `Quadrature` properties); retired 2 pre-existing `# type: ignore`.
- **B2 §1 PoleAngularClosure retirement — DEFERRED to #236** (user decision, 2026-06-22). The explorer's premise was WRONG, verified against main: the ABC `PoleAngularClosureBase` is `__call__`-only; the strategy methods (`precompute_psi_state`/`cell_contribution`/`angular_adjoint`/`level_indices`/`psi_half_seed`) are `MorelMontryAngularSweep`-specific; `IdentityAngularClosure` (cartesian) can't carry them as abstractmethods; #236's ABC-contract work (`6fdb0fe`) is **unmerged** (merge-status trap). Retiring the orphan would require narrowing consumers to the concrete (= B3 work on loss_representation.py) and would be superseded once #236 hoists the methods to the ABC — so keep the orphan Protocol until #236 lands. **#248 updated.** The independent half (AngularMeasure) shipped in `aa93757`.
- **B3 full_field_space** (`11fa279`, −6): `_require_blocks()` now RETURNS the narrowed `(bulk_space, trace_space)` pair so the three callers bind non-None locals (it was `-> None`, which can't narrow `self.x`). Runtime-inert. The remaining `.bulk`/`.boundary`-on-`NDArray` + the `inner_product` override error in that file are one deeper mismatch — `FunctionSpace` is typed for `NDArray` while `FullFieldSpace` carries composite fields → **B4** generic-hierarchy. LESSON: matching an override param *name* to a base whose *type* is wrong only surfaces the type mismatch harder (the rename was reverted).

### Session 2026-06-22 (cont., post-compaction) — orpheus/ **470 → 440** (`b36ca11 → bff6135`, on branch `refactor/pyright-signal-2`)

Two ROOT-CAUSE structural fixes (both more than local narrowing — they make the type tell the truth), each verified `-O` + ratchet-retightened:

- **`_bases.py` TraceSpace.layout** (`382fe5f`, −3): `BoundaryField.from_face_arrays` read `layout = space.layout` after the `space is None` guard, but `TraceSpace.layout` is `Optional` (bare-constructor footgun). Added an explicit `layout is None` guard (parse-don't-validate, mirrors `_face_row` / `_require_blocks`). #236-independent.
- **`StreamingTerms` required fields** (`dae34e0`, −24 across diamond.py + cell_balance.py): the 10 curvature fields were typed `float | None = None`, contradicting the file's own **Step 2.5 invariant** (all populated for every geometry; slab neutral, curvilinear physical — to retire the `alpha_in is None` branch). Made all 10 required `float` (Pattern 4). **explorer audit** confirmed all 3 production factory arms populate every field by keyword + no production None-consumer. Riders: deleted the now-dead `abs_mu/volume is None` guard in `linear_discontinuous.py`; fixed the self-contradicting docstring (one para said the branch was retired, a later stale one said it was live); migrated 5 test constructors that omitted `mu_start`. diamond/cell_balance/reduced_operator/linear_discontinuous all reach 0.
- **`AngularMeasure` Protocol read-only properties** (`b578286`, −3): the Protocol declared `mu_x`/`weights`/`N`/`eta`/`mu_z`/`level_indices` as mutable attrs → invariant → a concrete read-only `@property` (Quadrature) fails to match ("property is not assignable to ndarray"). Declared them read-only `@property` (correct covariant-read contract; static-only, Protocols never instantiate). Cleared geometry.py:414/420/427 AND the test-suite flood. 509 geometry+sweep tests green.
- Ratchet: `9c6f037` (470→443), `bff6135` (443→440). Tests verified: `tests/sn/sweep/core/` (460), `tests/geometry/test_reduced_operator.py` + sweep/core (509).

**FINDING — #250 is real and pre-existing on main.** Two sphere bit-identity snapshot reds (`test_sphere_{1g,2g}_apply_bit_identical`, 100% mismatch, rel ~3.4) surfaced in the broad run. Verified via a **worktree at clean `b36ca11`** (PYTHONPATH override; editable install pins imports) — they fail IDENTICALLY on untouched main, so NOT caused by this work. Already tracked: **#250** ("re-baseline stale SPHERE snapshots — 5 curvilinear reds from the W1 clamp fix `b2d8a6d`"). Out of pyright scope; a V&V re-baseline (numerics-investigator confirms W1 values, then re-capture `pre_t4_snapshots.npz`). The slab arms already moved to `assert_regression`; the sphere arms still `assert_allclose` the stale npz.

**Two more clusters this session** (`114b679 → 105dbf8`, merged+pushed; ratchet 440→428):
- **registry SelectionLog typed closure** (`457d2e3`, −8): `select_quadrature` splatted a union-valued `dict(...)` (`**log_template`) into two `SelectionLog(...)` constructions (8 errors = 4 fields × 2 sites). Replaced with one typed closure `_selection_log(chosen_spec, chosen_parameters)` — every field passed with its precise type. registry.py 8→0; 82 quadrature tests green.
- **Axis1D Protocol read-only properties** (`77aaf7c`, −4): SAME variance fix as AngularMeasure — `edges`/`coord` were mutable attrs (invariant) rejecting the frozen-dataclass / `@property` implementers; declared read-only `@property` (matching the already-property `n`/`endpoints`/`bc`). axis.py 4→0 + cleared the `coord_system(...)`/`from_axes(...)` test flood. 97 axis+sweep tests green.
- **Protocol-variance sweep is now EXHAUSTED in production** — a precise grep for `is not read-only in protocol` / `invariant because it is mutable` (non-deriv) leaves only **1** straggler (solver.py:753, a different union shape). The two real instances (AngularMeasure, Axis1D) are caught.
- **geometry.py dag_walk narrowing** (`d6c99d9`, −8, ratchet 428→420): the two-entry-mode `dag_walk` (`ordinate_idx` XOR `direction_sign`) + `dag_walk_cell_indices` + `_representative_ordinate` couldn't narrow `mu_level_idx`/`ordinate_idx`. Added `_require_mu_level(mu_level_idx) -> int` (Pattern 2, single-sources the 3 cylindrical guards, replaces the 2 early compound guards) + an `ordinate_idx is None` post-dispatch narrowing. geometry.py 12→4. **Full SN sweep suite green (679 passed)** — hot path, all geometries verified. No test asserted the moved guard messages (checked first).

### Resume here (next) — orpheus/ at **420**
The clean single-root-cause autonomous wins are now **spent**. geometry.py's 4 residuals are distinct one-offs (`Mesh2D.areas` L369; `pole_angular_closure` assign/positional L465 = #236; one `None` subscript L768). Remaining `orpheus/sn` (excl. #236-deferred `pole_angular_closure` 34 + `loss_representation` 43) is **heterogeneous** and falls into harder buckets that need user-steering or deep flow work:
- **B4 — `LinearOperator[V]` generic → being done as the OPERATOR-INVERSE-ALGEBRA CARVE (a DETOUR).** Not a typing patch: a deliberate architectural carve toward the grand report's *inverse-as-operator* design (`solve` retired to `inverse().apply`; sweep-vs-Krylov = which inverse operator; domain/codomain typevar split). Full plan: **`.claude/plans/operator_inverse_algebra_carve.md`** (cold-pickup-ready, 6-phase migration, exits back HERE in its §8). Verification net: `.claude/plans/issue_226_b4_operator_generics_verification.md`; design map: `.claude/agent-memory/explorer/issue_226_operator_generics_map.md`. SURGICAL / main-agent-direct + user-steered (`delegation.md`) — do NOT batch-autonomously, do NOT dispatch method-implementer. When the carve merges, B4's 18 reds clear by architecture; resume the burndown at B5 below.
- **B5 — union dispatch**: scattering.py (13: `ndarray | ScalarSourceSink | None`, `.mesh`/`.values`/`.spatial_moments_per_axis` on raw ndarray), solver.py (`BulkField` vs `AngularField`). Needs source-build flow understanding.
- **units.py (5)** — pint-stub-fighting (`PlainUnit | Unknown` vs declared `Unit`); third-party-stub category like the SymPy backlog. Low value; defer with B6.
- Then **B6**+**Workstream C** (derivations ~228), **Workstream A** (user-gated hook).

Operating rules in force: `npx pyright` CLI is the oracle (the streamed `<new-diagnostics>` lags edits — never trust it over a CLI run); NO `# type: ignore`; re-tighten `pyright_baseline.json` after each cluster (`python -m tests._harness.pyright_ratchet --update`); per-cluster ff-merge + push (user authorized this cadence). **Corrected paths**: ratchet test = `tests/test_pyright_ratchet.py`; `test_reduced_operator.py` is under `tests/geometry/`.

---

## Current state (triage, 2026-06-21, `npx pyright` 1.1.410)

Full repo = **2353 errors / 19 warnings**, 774 files. By directory:

| dir | errors | nature |
|---|---|---|
| `tests/` | 1403 | test-stub / `**kwargs`-splat / loose-fixture idioms |
| `orpheus/` | **502** | production — the real target (260 non-derivation + 242 SymPy) |
| `derivations/` (repo-root) | 267 | SymPy-fighting |
| `scratch/` | 129 | throwaway |
| `student_resources/` | 37 | teaching |
| `examples/` | 12 | demos |
| `tools/` | 3 | |

**No regression** — the #257 branch *reduced* `orpheus/` 538→502 (−36, all `orpheus/sn`). The 2353 is full-repo scope; historical "538" was always `orpheus/`-only (`pyrightconfig.json` has no `include`/`exclude` → bare `npx pyright` walks the whole root).

Production (`orpheus/`) rule histogram: reportArgumentType 224 · reportAttributeAccessIssue 98 · reportOptionalSubscript 38 · reportReturnType 25 · reportOperatorIssue 19 · reportAssignmentType 15 · reportCallIssue 14 · reportOptionalMemberAccess 14 · reportUndefinedVariable 11 · …

---

## Workstream A — NOISE elimination (do FIRST; it's the "not noise" half)

The streamed `<new-diagnostics>` reportMissingImports avalanche ("Import 'orpheus.transport.full_field' could not be resolved", ".spatial.scheme", etc.) is **0-real in the CLI for `orpheus/`** — it's the langserver dual-identity artifact (the worktree/exec-env rooting + microsoft/pyright#10498). The fix lives in `.claude/hooks/regen-pyrightconfig.sh` (the per-worktree `executionEnvironment` design) — the two `.claude/worktrees/*` checkouts already contribute 0 errors via that mechanism.

**A1. Diagnose the residual LSP noise.** Why does the LSP still stream phantom import errors when the CLI is clean? Confirm whether the in-editor/agent langserver is using the regenerated `pyrightconfig.json`, whether stale worktree exec-envs or a missing root mapping cause the dual-identity, and whether microsoft/pyright#10498 has a config-level mitigation. Dispatch explorer + read the hook + reference_pyright_lsp_rooting memory.
**A2. Make the LSP agree with the CLI.** Extend/repair `regen-pyrightconfig.sh` (or the langserver invocation) so the streamed diagnostics match `npx pyright orpheus/` — phantom imports gone. This is the single biggest agent-UX win (every turn is currently polluted). **This file is in `.claude/hooks/` (commit-protected) — propose the change to the user, do not self-commit it.**
**A3. Pin the scope.** Decide + document the canonical analysis scope: should bare `npx pyright` analyze the whole repo, or should `pyrightconfig.json`/`pyproject.toml [tool.pyright]` set `include = ["orpheus"]` (production) with tests/derivations/scratch under separate, relaxed exec-envs? Pinning scope makes "the count" stable and meaningful (no more 2353-vs-538 confusion).

---

## Workstream B — real-error reduction (the "only signal" half)

Attack order = highest signal-cleared-per-unit-effort, re-measuring between big clusters. Each cluster: explorer/enumerate → fix → tests green + count drops + 0 net-new ignores → commit.

**B1. Cluster D + undefined-vars (~17, trivial, zero-risk) — warm-up.** `numpy.bool_`/`csr_array` vs annotated return → `bool(...)`/correct annotation (~6). The 11 `reportUndefinedVariable` = missing `import sympy as sp` (7×), `Quadrature1D`/`TimedFullField` TYPE_CHECKING forward-refs (4×). Mechanical, runtime-inert.

**B2. Cluster C — Protocol/ABC missing concrete attrs (~12, low-risk).** `PoleAngularClosure` missing `level_indices`/`precompute_psi_state`/`cell_contribution` (8×); `AngularMeasure`, `MOCMesh` (rest). Declare the members on the ABC/Protocol. ⚠ VERIFY against the #248 `PoleAngularClosure`-Protocol retirement first (some accesses may be on a stale name).

**B3. Cluster A — Optional/None not narrowed (~104, the biggest single bucket).** reportOptionalSubscript 38 + OptionalMemberAccess 14 + OptionalOperand 5 + the `... | None` argtype/operator errors. **Two files hold ~77:** `orpheus/sn/loss_representation.py` (43) + `orpheus/sn/spatial/pole_angular_closure.py` (34). Root: lazily-populated `Optional`/`= None` fields (`_alpha_per_level`, `_tau_per_level`, `_dAw_per_level`, `psi_half_seed`) subscripted after an invariant guarantees they're set; + `float | None` params used in arithmetic. Fix per-site (judgement): declare the real non-Optional type + init in `__post_init__`/builder, OR `assert x is not None` at the invariant-established point, OR make the param non-Optional. ⚠ Respect the lazy-init contract — don't force eager allocation where the None sentinel is load-bearing. **Re-measure after** — narrowing one lazy field clears cascading downstream argtype/operator errors.

**B4. Cluster B — `LinearOperator[V]` under-parametrised generic (~24+, architecturally coupled).** `LinearOperator[Unknown]` loses `.solve`/`.apply_transpose`; `block_role` "cannot assign"; `TensorProductOperator`/`IncomingSourceOperator` not assignable (7×). Concentrated in `boundary_realizer.py` (9), `operator.py` (3), `iteration.py` (3), `solver.py` (3). Fix = ONE coherent hierarchy edit per `.claude/agent-memory/explorer/issue_226_operator_generics_map.md` (unbounded `V=TypeVar("V")`): parametrise the `LinearOperator`/`Mixin`/`OperatorSum` family `Generic[V]`; declare `block_role` as a settable field on the composers; surface `.solve` via the resolvent type / a Protocol. Type-level, runtime-inert; re-measure (likely clears uncounted downstream argtype collapse).

**B5. Cluster E — stringly/union dispatch (~10, ergonomics).** `str|SubgroupOfO3|int|dict` arg in `numerics/quadrature/registry.py` (8×); `Mesh1D|Mesh2D|tuple[Axis1D,...]` mesh arg (2×); `Quadrature→AngularQuadrature` (3×). Tighten signatures / overload / narrow with `isinstance`.

**B6. Cluster F — `orpheus/derivations` SymPy backlog (242, ISOLATE, last).** Unchanged from main; SymPy `Expr`/`Float`/`Rational` fight pyright. Fix the genuinely-real ones, then a **SCOPED per-directory `executionEnvironment`** relaxation for `orpheus/derivations/**` (NOT global) — the triage doc's recommended handling. Lowest production value.

---

## Workstream C — the non-production trees (1403 tests + 267 root-derivations + 129 scratch + …)

Lower value, larger count. Decide policy rather than grind every one:
- **`tests/`**: many are test-stub idioms (`_StubQuad` not a `Quadrature`, `**kwargs`-splat into typed factories, `BC.vacuum` enum access, `float`→`int` arg looseness). Triage: fix the cheap structural ones (the `solve_*(geometry, **dict)` splat collapses pyright cleanly into a typed helper — see the S10a peierls-test refactor precedent); for the rest, consider a relaxed test-tree exec-env (tests are not shipped product). DO NOT let test-stub noise gate production signal.
- **`derivations/` (repo-root), `scratch/`, `student_resources/`, `examples/`**: scratch/teaching — strong candidates for scoped exec-env relaxation or `exclude`, NOT line-by-line fixes.

The principle: production (`orpheus/`, ex-derivations) → **zero real errors**; the SymPy/test/scratch trees → **scoped relaxation** so they don't drown the signal, with the genuinely-real ones fixed.

---

## Sequencing
1. **A1–A3 noise-rooting FIRST** (kill the streamed avalanche — makes all subsequent work legible; the hook change goes to the user since `.claude/hooks/` is protected).
2. **B1 → B2** (zero/low-risk warm-up, ~29 cleared).
3. **B3** (the 104 Optional/None bucket; re-measure).
4. **B4** (the LinearOperator-Generic coherent edit; re-measure).
5. **B5**, then **B6** (derivations scoped relaxation).
6. **Workstream C** policy (scope/relax the non-production trees).

## Verification (every batch)
- `npx pyright orpheus/` count strictly drops; report before/after.
- 0 net-new `# type: ignore` (grep the diff).
- The relevant test suites stay green under `-O` (route around the documented baseline reds #250/#232/#212).
- After A2: the streamed `<new-diagnostics>` for an `orpheus/` edit matches the CLI (no phantom imports).

## Open investigations (dispatch at execution)
- explorer: the LSP dual-identity root cause + whether `regen-pyrightconfig.sh` fully resolves it (A1).
- explorer: confirm the #248 `PoleAngularClosure` retirement status before B2.
- the `issue_226_operator_generics_map.md` memo is the B4 design input — re-read it.

## Relations
#226 (the pyright-signal effort), the prior `.claude/plans/pyright_cluster_triage.md` (superseded counts: 691→502 orpheus, undefined-var 67→11, missing-imports 0), `reference_pyright_lsp_rooting.md`, `issue_226_operator_generics_map.md`. #258 (units.py pint-stub debt) + #254/#255 feed here.
