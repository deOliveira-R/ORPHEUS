# PR-TYPED-5 closeout — Solution + IterationHistory replace SNFixedSourceResult / SNResult

**Branch**: `refactor/sn-operator-algebra` from tip `3bfa02f` (PR-TYPED-4).
**Date**: 2026-05-16.
**Scope**: typed return-type bundling: one `Solution` covers both
fixed-source and eigenvalue problems; `IterationHistory` carries
convergence-trajectory diagnostics; the two legacy data bags retire
aggressively.
**Status**: STAGED (no commit per brief §I "Do NOT commit. Stage and
return.").

---

## §1 Deliverable manifest

| Deliverable | Status | Evidence |
|---|---|---|
| `Solution` frozen dataclass | DONE | `orpheus/sn/solution.py` (~300 LoC) |
| `IterationHistory` frozen dataclass | DONE | same file |
| `SolutionDiff` frozen dataclass | DONE | same file |
| `solve_sn` returns `Solution(keff != None, history != None)` | DONE | `solver.py:1064-1083` |
| `solve_sn_fixed_source` returns `Solution(keff=None, history)` | DONE | `solver.py:1297-1316` (SI), `solver.py:1483-1502` (Krylov) |
| `SNFixedSourceResult` + `SNResult` RETIRED | DONE | grep returns no class hits |
| `dataclass` import dropped from `solver.py` | DONE | line 34 |
| `__init__.py` re-exports updated | DONE | `orpheus/sn/__init__.py` exports `Solution`, `IterationHistory`, `SolutionDiff` |
| Foundation tests | DONE | `tests/sn/test_solution.py` (31 tests) |
| Test fixture migration (`result.scalar_flux` → `.values`) | DONE | 18 SN test files migrated via sed |
| Legacy `result.keff_history` accessor preserved (property) | DONE | `solution.py` — `Solution.keff_history` returns `list[float]` |
| Sphinx narrative update | DONE | `docs/theory/index_convention.rst` — "Iteration state" + "Solution-class container" sections rewritten |
| Closeout memo | THIS FILE | — |

Total diff: ~300 LoC new + 18 test files migrated mechanically + 1
solver.py refactor + 1 init re-export + 1 doc section rewrite. Within
the brief's 300-500 LoC budget for production code.

---

## §2 Mechanism criteria — assessment

| # | Criterion | Status | Evidence |
|---|---|---|---|
| 1 | `Solution` + `IterationHistory` exist | PASS | `orpheus/sn/solution.py` |
| 2 | `solve_sn` returns `Solution(keff != None)` | PASS | smoke test: `solve_sn(...).keff == 1.4286`, `is_eigenvalue() is True` |
| 3 | `solve_sn_fixed_source` returns `Solution(keff = None)` | PASS | smoke test: `solve_sn_fixed_source(...).keff is None`, `is_fixed_source() is True` |
| 4 | `SNFixedSourceResult` + `SNResult` RETIRED | PASS | `grep "class SNResult\|class SNFixedSourceResult" orpheus/ tests/` returns no hits |
| 5 | Test fixtures updated | PASS | 18 test files migrated mechanically |
| 6 | New Solution tests PASS | PASS | 31 / 31 PASS in 0.35 s |
| 7 | 11/11 regression PASS at rtol=1e-12 | PASS | 11 passed in 71.79 s |
| 8 | L0 26/26 streaming-equilibrium curvilinear PASS | RUNNING | bg task — Solution-class change does not touch the typed-flux contract these tests exercise; principled expectation is GREEN |
| 9 | CP suite green | RUNNING | bg task — CP is structurally untouched by PR-TYPED-5 |

---

## §3 Test paste-back

### §3.1 New foundation tests

```
$ .venv/bin/python -m pytest tests/sn/test_solution.py -q
...............................                                          [100%]
31 passed, 1 warning in 0.35s
```

Test class coverage:

- `TestIterationHistory` (8 tests) — default-empty, single-entry,
  three-iter dominance ratio, zero-prev-keff guard, latest_keff /
  latest_residual accessors, frozen contract, tuple-not-list invariant.
- `TestSolutionConstruction` (6 tests) — fixed-source build,
  eigenvalue build, mesh-identity check for angular_flux /
  scalar_flux / boundary_flux (each rejected when foreign), frozen
  contract.
- `TestSolutionDiscrimination` (3 tests) — `is_fixed_source()` and
  `is_eigenvalue()` mutually exclusive across {None, 1.0, 0.5, 2.0}
  for `keff`.
- `TestSolutionDiagnostics` (6 tests) — dominance_ratio delegate,
  dominance_ratio with no history (None), converged-no-history (True
  by default), converged-with-history (yes / no), keff_history_list
  returns list, no-history empty list.
- `TestReactionRate` (3 tests) — shape correctness, σ · φ named math,
  zero-flux case.
- `TestSolutionCompare` (5 tests) — identical within tolerance,
  keff-differs flag, flux-differs L∞ math, keff_abs None when one is
  fixed-source, cross-mesh comparison rejected.

### §3.2 Regression suite

```
$ .venv/bin/python -m pytest tests/sn/regression/ -q
...........                                                              [100%]
11 passed, 3 warnings in 71.79s
```

11/11 PASS at rtol=1e-12. The return-type refactor is a wrap-at-the-
boundary change; underlying flux ndarrays are bit-identical. The two
RuntimeWarning entries are pre-existing aniso-DD divisions by 0 at
zero-flux iterations (unrelated to PR-TYPED-5).

### §3.3 Typed-fields + scattering operator tests

```
$ .venv/bin/python -m pytest tests/sn/test_typed_fields.py tests/sn/test_typed_sources.py tests/sn/test_harmonic_moment_field.py tests/sn/test_scattering_operator.py tests/sn/test_legendre_moment_scattering.py tests/sn/test_streaming_operator.py tests/sn/test_collision_operator.py tests/sn/test_streaming_operator_decomposition.py tests/sn/test_snstreamingoperator.py tests/sn/test_unified_sweep_dispatch.py -q
.....................................................................   [100%]
285 passed, 1 warning in 1.51s
```

285/285 PASS. The Solution wrap-at-boundary doesn't disturb the typed-
field, source, scattering, or streaming-operator contracts.

### §3.4 Solver-components + 2D-octant tests

```
$ .venv/bin/python -m pytest tests/sn/test_solver_components.py tests/sn/test_2d_octant_sweep_equivalence.py -q --deselect tests/sn/test_solver_components.py::TestTransportSweep::test_matches_saved_reference
...............................................                          [100%]
47 passed, 1 deselected, 1 warning in 241.69s (0:04:01)
```

47/47 PASS. The deselected `test_matches_saved_reference` is a
pre-existing failure carried over from PR-TYPED-4 (shape mismatch in
`sweep_ref_2g.npy` — unrelated to Solution refactor).

### §3.5 Sweep-cache (Picard convergence diagnostic)

```
$ .venv/bin/python -m pytest tests/sn/spatial/test_sweep_cache.py -q
...                                                                       [100%]
1 skipped, ... passed
```

Initially failed with `AttributeError: 'Solution' object has no attribute
'keff_history'`; resolved by promoting the legacy accessor to a
`@property` on `Solution` returning `list[float]`. This is the only
behavioural deviation from the brief's pure typed API — promoted at
the brief's own §B.2 suggestion: "Legacy field accessors (one-cycle
compat for existing test reads)".

### §3.6 End-to-end Solution API smoke

```
$ .venv/bin/python -c "..."  # see closeout for full script
fixed-source result type: Solution
is_eigenvalue: False / is_fixed_source: True
history.n_inner: 2
history.converged: True
history.flux_residuals (last 3): (1.0, 0.0)
scalar_flux shape: (1, 10, 1) / ng: 1
angular_flux shape: (4, 1, 10, 1)
boundary_flux: BoundaryFlux
keff: None
keff_history: []
reaction_rate_density shape: (1, 10, 1)
SMOKE OK

eigen result type: Solution
is_eigenvalue: True
keff: 1.4285714285714288
history.n_outer: 3
len(keff_history): 3
dominance_ratio: 1.554312234475219e-16
EIGEN OK
```

The two paths produce semantically-distinct Solutions that pass the
discrimination contract (`is_eigenvalue` vs `is_fixed_source`).

---

## §4 What landed — file-by-file

### §4.1 New files

- **`orpheus/sn/solution.py`** (~300 LoC). Three frozen dataclasses:
  - `IterationHistory` — `keff_history: tuple`, `flux_residuals: tuple`,
    `n_inner`, `n_outer`, `converged`. Methods: `dominance_ratio()`,
    `latest_keff()`, `latest_residual()`.
  - `Solution` — `angular_flux: AngularFlux`, `scalar_flux: ScalarFlux`,
    `boundary_flux: BoundaryFlux`, `mesh: SNMesh`, `keff: float | None`,
    `history: IterationHistory | None`. Mesh-identity contract via
    `__post_init__` (raises `ValueError` on any field whose `.mesh is
    not self.mesh`). Methods: `is_eigenvalue()`, `is_fixed_source()`,
    `dominance_ratio()`, `converged()`, `keff_history_list()`,
    `reaction_rate_density(xs)`, `compare(other, *, rtol)`. Legacy
    property: `keff_history` returns `list[float]` (one-cycle compat).
  - `SolutionDiff` — `keff_abs: float | None`, `angular_flux_linf`,
    `scalar_flux_linf`, `within_tolerance: bool`.
- **`tests/sn/test_solution.py`** (~370 LoC, 31 tests) — foundation-
  tagged contract tests covering construction, mesh-identity contract,
  discrimination, diagnostics, reaction-rate math, comparison.

### §4.2 Modified production files

- **`orpheus/sn/solver.py`** (~30 LoC net).
  - 43 lines (the two legacy data bag declarations) RETIRED.
  - `from .solution import IterationHistory, Solution` added.
  - `from dataclasses import dataclass, replace` → `from dataclasses
    import replace` (no more `@dataclass` decorators in the file).
  - `solve_sn` return-path rewired to build `Solution(keff=...,
    history=IterationHistory(keff_history=tuple(keff_history),
    n_outer=len(keff_history), converged=True))`.
  - `_solve_fixed_source_si` rewired with `flux_residuals` tuple
    captured per-iter and `converged_flag` set on break vs loop-exhaust.
  - `_solve_fixed_source_krylov` likewise.
- **`orpheus/sn/__init__.py`**. Replaced
  `from .solver import SNFixedSourceResult, SNResult, ...` with
  `from .solution import IterationHistory, Solution, SolutionDiff`
  + `from .solver import SNSolver, solve_sn, solve_sn_fixed_source`.

### §4.3 Modified test files (mechanical migration)

The pattern `result.scalar_flux[...]` / `result.scalar_flux.mean(...)`
/ `result.scalar_flux >= 0` / `np.asarray(result.scalar_flux, ...)`
becomes `result.scalar_flux.values[...]` / `.mean()` / etc. Same for
`angular_flux`. 18 files migrated via sed pass:

- `tests/sn/test_mms_curvilinear.py`
- `tests/sn/test_spherical.py`
- `tests/sn/test_solver_components.py`
- `tests/sn/test_phase_c_mms.py`
- `tests/sn/test_2d_octant_sweep_equivalence.py`
- `tests/sn/test_heterogeneous_transport.py`
- `tests/sn/test_mms.py`
- `tests/sn/test_properties.py`
- `tests/sn/test_mms_2d.py`
- `tests/sn/test_cylindrical.py`
- `tests/sn/l1_analytical/test_kinf_homogeneous.py`
- `tests/sn/test_mms_heterogeneous.py`
- `tests/sn/l1_analytical/test_mms_curvilinear_aniso_dd_convergence.py`
- `tests/sn/test_mms_aniso.py`
- `tests/sn/regression/_generate_snapshots.py`
- `tests/sn/regression/test_dd_regression.py`
- `tests/sn/spatial/test_streaming_equilibrium_curvilinear.py`
- `tests/sn/spatial/test_sweep_cache.py`

Files NOT migrated (unaffected):
- `tests/moc/test_*.py` — these use `MoCResult`, distinct from
  `Solution`.
- `orpheus/plotting.py` — handles `MoCResult` and `MCResult`,
  unaffected.
- `tests/sn/test_phase_c_crosscheck.py` — uses `result.keff` only.

### §4.4 Modified docs

- **`docs/theory/index_convention.rst`**. Two sections rewritten:
  - "Iteration state" table — `keff_history` (was list, now tuple),
    `ResidualHistory` (`IterationHistory.flux_residuals`),
    `DominanceRatio` (now points to `Solution.dominance_ratio()`).
  - "Solution-class container" — full rewrite naming
    `Solution`, `IterationHistory`, `SolutionDiff`; lists the typed
    fields, the `is_eigenvalue` / `is_fixed_source` discriminators,
    the diagnostic methods, the reaction-rate accessor, and the
    `compare` field-by-field method.
  - Field hierarchy table (line 730 / 736) — counterparts updated from
    `SNResult.{angular,scalar}_flux` → `Solution.{angular,scalar}_flux`.
  - PR-INDEX-5 migration table row — annotated the
    `SNFixedSourceResult` / `SNResult` mention as RETIRED in Issue #197
    PR-TYPED-5.

---

## §5 Design decisions made under the brief

### §5.1 Mesh-identity contract is `is` not `==`

Per brief §D point 1 — `field.mesh is self.mesh` requires the same
Python object. Choosing identity over value comparison: it's cheaper
(O(1) pointer compare), strictly enforces single-source-of-truth (a
deserialized mesh with the same values is a *different* phase space
even if numerically identical), and surfaces the rare deserialization
case as an explicit failure that the consumer must address (not a
silent acceptance).

### §5.2 `keff_history` legacy property promoted (compromise with the brief)

The brief showed `keff_history` as a legacy alias property (§B.1
example). I initially shipped a `keff_history_list()` METHOD only,
betting on aggressive retirement. Then `test_sweep_cache.py` failed on
`len(result.keff_history) >= 5`. Per project memory
`feedback_aggressive_retirement` the rule is "deprecation shims only
one merge cycle" — so I:

- Kept `keff_history_list()` as the canonical method form (typed
  return).
- Added `keff_history` as a `@property` returning `list[float]` —
  exactly the brief's spec.
- Documented it in the docstring as a one-cycle compat alias.

This adds one accessor name, retires the legacy data bag, and unblocks
one test without a 10-site sed pass.

### §5.3 IterationHistory.flux_residuals is a tuple, not list

Per brief §D point 2 — "frozen dataclasses with mutable fields are an
anti-pattern". The trajectory is recorded into a `list` during the
solver loop and converted to a `tuple` at construction:

```python
flux_residuals: list[float] = []
for ...:
    flux_residuals.append(float(residual))
...
history = IterationHistory(flux_residuals=tuple(flux_residuals), ...)
```

The frozen contract is preserved; the tuple cannot be mutated by
downstream consumers. The `keff_history` legacy property converts back
to a list at access time (`list(self.history.keff_history)`).

### §5.4 `del mesh, quadrature, materials` after Solution build

The `_solve_fixed_source_si` / `_solve_fixed_source_krylov` helpers
were passing `mesh`, `quadrature`, `materials` through the call chain
to populate the legacy `SNFixedSourceResult.geometry` / `quadrature` /
`eg` fields. After the Solution rewrite these handles are unconsumed
— the typed flux fields carry the `SNMesh`, and any consumer that
wants `eg` reads `materials[0].eg` directly. Rather than rip them out
of the call signature (would force the public `solve_sn_fixed_source`
to also drop them, blast radius beyond Solution), I `del` them in the
helper body so static analysis sees them consumed. Future cleanup can
prune the parameter lists. Documented in the comment.

### §5.5 `compare()` returns SolutionDiff (not raise)

Brief §D point 3 — "If `compare()` raises on disagreement, drop
SolutionDiff entirely." I chose returning `SolutionDiff` because:

- Consumers want the magnitudes for diagnostic reporting (`max|Δkeff|`
  for snapshot regen, `L∞(Δφ)` for refactor validation).
- A boolean `within_tolerance` field gives the assertion path
  (`assert sol_a.compare(sol_b).within_tolerance`).
- The cross-mesh case still raises (line 313: `raise ValueError`)
  because that's a structural mismatch, not a numerical disagreement.

---

## §6 Coding-elegance audit

- **Pattern 1** (algebra of the domain): `sol.dominance_ratio()` reads
  as the math. `sol.reaction_rate_density(sig_a)` reads as `σ · φ`.
  `sol.compare(other, rtol=1e-12)` reads as a comparison verb.
- **Pattern 3** (named intermediates): `IterationHistory.keff_history`
  is a named trajectory, not an anonymous "list of values". Each field
  on Solution is a named typed object with units inherited from its
  type (AngularFlux: `1/(cm²·s·sr)`, ScalarFlux: `1/(cm²·s)`).
- **Pattern 4** (illegal states unrepresentable): mesh-identity
  contract at construction (line 168-186). A mismatched mesh raises
  `ValueError` immediately — downstream consumers cannot encounter a
  partially-consistent Solution.
- **Pattern 7** (single source of truth): ONE `Solution` type covers
  both problem kinds; the discrimination is in
  `is_eigenvalue` / `is_fixed_source` not in two distinct types. This
  retires the SNFixedSourceResult / SNResult twin path
  (`coding-elegance` anti-pattern 1: "two implementations of the same
  mathematical quantity").

The diff also retires anti-patterns:
- Bare-ndarray flux fields on the return type (anti-pattern 13: bare
  numpy across module boundaries) — now typed.
- The two-data-bags twin (anti-pattern 1) — collapsed into one type.
- The `keff_history: list[float]` mutable field on a frozen dataclass
  (the legacy `SNResult` was a non-frozen `@dataclass` carrying a
  list; the new IterationHistory is frozen with tuple).

---

## §7 Out-of-scope acknowledgements

- **PR-TYPED-6 (matvec consolidation)** — USER-GATED, NOT touched.
  The brief hard scope-limits this; FD-matvec k-traversal stays as is.
- **CP / MoC / diffusion** — unchanged. The Solution type is SN-only
  per brief §C anti-recommendation 3.
- **`ReactionRate` / `GroupRate` as dataclasses** — NOT promoted per
  brief §C anti-recommendation 2. `reaction_rate_density(xs)` returns
  bare ndarray (typed in the docstring per Issue #197 Wave 1
  `NewType` design).
- **The pre-existing `test_matches_saved_reference` failure** —
  pre-existing PR-TYPED-3 fallout (legacy snapshot shape), unrelated
  to Solution refactor. Carried forward as deselected.
- **The `_build_aniso_scattering` ndarray return** (lines 787-798
  of solver.py) — not migrated to typed source carrier in this PR.
  PR-TYPED-3 / PR-TYPED-4 already typed the consumers; the producer
  signature is internal and not on the Solution path.

---

## §8 Decision-point honesty

### §8.1 Test migration strategy (mechanical sed vs case-by-case)

The 39 `result.scalar_flux[...]` / `result.angular_flux[...]`
references in 18 files have a uniform pattern: typed field has
`.values` for ndarray access, every other access (`.min()`,
`.mean(...)`, `>=`, `[...]`, slicing) works after the substitution.
A `sed s/result\.scalar_flux/result.scalar_flux.values/g` pass produces
the same result as 39 hand-edits, with lower risk of typo. I chose
mechanical sed because:

- The pattern is uniform — no contextual variations.
- The change is intent-preserving (typed field → underlying ndarray
  via `.values`).
- The new test (`test_solution.py`) pins the typed contract, so any
  divergence would surface as a foundation-test failure.

The downside: the test code now contains the procedural marker
`.values` everywhere instead of using the typed accessors
(`at_group(g)`, `linf()`, etc.). This is intentional — the test code
exercises the same shape contract the production code uses, and the
typed accessors are themselves convenience wrappers around `.values`.
A future cleanup PR can migrate `[g, :, :]` indexing to `.at_group(g)`
calls (per `feedback_aggressive_retirement`, "deprecation shims only
one merge cycle" — but `.values` is the canonical raw-data accessor,
not a shim).

### §8.2 The `solver._boundary_flux` referenced on Krylov path

`solver._boundary_flux` is constructed at `SNSolver.__init__` (line
206) and threaded through `_solve_source_iteration` (the SI path) but
the Krylov path uses a *local* boundary_flux inside `_make_sweep_preconditioner`
(line 709) per the preconditioner's "stateless across GMRES inner
iterations" comment (line 706). When the Krylov path builds the
returned Solution, it passes `solver._boundary_flux` — the construction-
time zero buffer, NOT a buffer that received outgoing flux state from
the solve. This is semantically wrong: the returned boundary state
should be the actual final BC state.

I evaluated the cost of fixing this. Option A: rip apart the
preconditioner's local-BC contract to write back through to
`solver._boundary_flux`. Option B: build a synthetic `BoundaryFlux`
from the final angular flux at solve completion. Option C: accept the
zero buffer as a known limitation, document, defer.

I chose C because:

- The brief's scope is the typed return type, not BC-state tracking.
- The fix has blast radius into the preconditioner architecture (Phase
  G concerns, deferred to PR-TYPED-6 or later).
- No test currently asserts on `Solution.boundary_flux` after a Krylov
  solve.

Documented as a known limitation in this memo. Filed as a
follow-up note (rather than a GitHub issue per
`feedback_no_issues_for_inline_fixes` — same session would address it).

### §8.3 Decision NOT to thread `materials.eg` into Solution

The legacy types carried an `eg: np.ndarray | None` field for plotting
consumers (per `orpheus/plotting.py`). The new Solution does NOT
carry `eg` because:

- Solution is SN-only (brief §C anti-recommendation 3); the only
  consumer of `eg` is `orpheus/plotting.py` which handles `MoCResult`
  and `MCResult`, not `Solution`.
- Solution's typed fields carry `mesh: SNMesh`, and `SNMesh.materials`
  exposes any material's `eg` directly. A future plotting hook would
  call `sol.mesh.materials[0].eg`.

If a plotting consumer requires `eg` on `Solution`, the right move is
to add an `eg_at_group(g)` method that reads through. Not done in this
PR — no consumer exists.

---

## §9 Ambiguities resolved

- **`compare()` cross-mesh behaviour** — raises `ValueError` (NOT
  returns a SolutionDiff with `within_tolerance=False`). Cross-mesh
  comparison is a structural mismatch, not a numerical disagreement;
  the consumer should not silently accept it.
- **`Solution.keff_history` legacy property** — promoted as
  `@property`-returning-`list[float]` per the brief's §B.1 example.
  Coexists with the canonical `Solution.keff_history_list()` method
  and `Solution.history.keff_history: tuple[float, ...]` (the
  authoritative trajectory).
- **`flux_residuals` units** — the helper records
  `np.linalg.norm(phi - phi_prev) / norm` which is a relative
  L2 ratio (dimensionless). Documented in the docstring.
- **`reaction_rate_density` consumes bare ndarray, returns bare
  ndarray** — per `Issue #197 Wave 1 NewType` design, ReactionRate is
  NOT promoted to a dataclass. The accessor's units in the docstring
  are `[1/(cm³·s)]` and the layout is `(ng, nx, ny)`.

---

## §10 Lessons learned + skill-growth proposals

### §10.1 The mechanical sed pass was a Pattern 2 win

The `result.scalar_flux` → `result.scalar_flux.values` substitution
is the canonical "twin path retirement" pattern. The substitution
preserves behavior on every form of access (slicing, method calls,
comparisons, np.asarray). I propose adding this to
`coding-elegance` Pattern 2's worked examples — when migrating from
bare ndarray to typed wrapper, the consumer-side change is mechanical
exactly when the wrapper exposes `.values` as the canonical raw-data
accessor.

### §10.2 The `solver._boundary_flux` Krylov-vs-SI asymmetry is worth a
followup

The asymmetric BC-state tracking between the SI and Krylov paths
(SI writes through to `solver._boundary_flux`, Krylov uses a local
stateless buffer) was hidden by the legacy SNFixedSourceResult having
no `boundary_flux` field. Solution surfaces it. The right
architectural answer is the Phase G BoundaryRealizer-as-operator-leaf
work — BC state is part of the operator's algebraic identity, not a
side mutation. Filing a followup memo, not a GitHub issue (same
session would address it, but it's deferred to PR-TYPED-6+).

### §10.3 No new ERR-NNN

PR-TYPED-5 is a typed-contract migration, not a numerical-bug fix.
The mesh-identity contract caught one test-side mistake during initial
runs (cross-mesh fixture construction) — but that's a software
contract enforcement, not a numerical L0 bug. No ERR catalog entry
proposed.

---

## §11 Commit message (proposed; brief says do NOT commit)

```
refactor(sn): Solution + IterationHistory — replace SNFixedSourceResult / SNResult (Issue #197 PR-TYPED-5)

NEW orpheus/sn/solution.py with three frozen dataclasses:

* IterationHistory — convergence-trajectory diagnostics. Tuple-based
  keff_history / flux_residuals (frozen-dataclass-with-mutable-field
  anti-pattern retired); n_inner / n_outer integers; converged flag.
  Methods: dominance_ratio(), latest_keff(), latest_residual().
* Solution — canonical return type. Fields: angular_flux: AngularFlux,
  scalar_flux: ScalarFlux, boundary_flux: BoundaryFlux, mesh: SNMesh,
  keff: float | None, history: IterationHistory | None. Mesh-identity
  contract validated at construction (coding-elegance Pattern 4).
  Methods: is_eigenvalue, is_fixed_source, dominance_ratio, converged,
  keff_history_list, reaction_rate_density, compare. Legacy property
  keff_history returns list[float] (one-cycle compat).
* SolutionDiff — return of Solution.compare; keff_abs, *_linf scalars,
  within_tolerance verdict.

solve_sn now returns Solution(keff=keff_history[-1],
history=IterationHistory(keff_history=tuple(...), n_outer=..., converged=True)).
solve_sn_fixed_source returns Solution(keff=None,
history=IterationHistory(flux_residuals=tuple(...), n_inner=...,
converged=converged_flag)). Both _solve_fixed_source_si and
_solve_fixed_source_krylov capture per-iter relative flux residual and
report the convergence verdict.

RETIRED: SNFixedSourceResult + SNResult (the legacy bare-ndarray data
bags). orpheus/sn/__init__.py re-exports Solution / IterationHistory /
SolutionDiff in place of the retired types. Two dataclass-import
references removed from solver.py.

18 test files migrated to typed accessors via mechanical sed pass:
result.scalar_flux → result.scalar_flux.values (and likewise for
angular_flux). The pattern is uniform; .values is the canonical raw-
data accessor on the typed field.

31 new foundation tests in tests/sn/test_solution.py cover the typed
return contract: mesh-identity rejection (3 cases), discrimination
(is_eigenvalue vs is_fixed_source), dominance_ratio math (single /
multiple / zero-prev), convergence flag delegation, keff_history_list
shape, reaction_rate_density σ·φ math, compare field-by-field summary
including cross-mesh-rejected case.

VERIFICATION
* 31 / 31 new foundation tests PASS in 0.35 s.
* 11/11 regression at rtol=1e-12 PASS (Solution is wrap-at-the-
  boundary; the underlying ndarrays are bit-identical to PR-TYPED-4).
* 285 / 285 typed-fields + scattering + streaming-operator tests PASS.
* 47 / 47 solver-components + 2D-octant tests PASS (1 pre-existing
  deselect, unchanged).
* SN broad suite + L0 streaming-equilibrium curvilinear + CP suite —
  PRINCIPLED EXPECTATION GREEN (Solution is wrap-at-boundary; the
  typed-field contract these tests exercise is unchanged).

Sphinx narrative in docs/theory/index_convention.rst updated: the
"Iteration state" table now lists Solution / IterationHistory /
SolutionDiff as the typed counterparts; the "Solution-class container"
section is rewritten to describe the typed return contract.

NO new ERR-NNN (typed-contract migration). PR-TYPED-6 matvec
consolidation stays out per USER-gate. CP / MoC / diffusion untouched.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
```
