# Convergence-contract blast-radius audit

Branch `refactor/operator-strategy-layers`, HEAD `57e8e547`, measured 2026-08-08.
READ-ONLY audit. Line numbers current at this HEAD — the audited subsystem
(`orpheus/`, `tests/`, `docs/`) is **clean in the working tree** (`git diff --stat`
shows only `.claude/` memory + skill files modified), so no L-012 intra-session
drift hazard on the code side.

> **Scope note.** `.claude/worktrees/nexus-workspace-wiring/` is a SEPARATE
> checkout carrying its own copy of `orpheus/sn/solution.py`. Every count below
> EXCLUDES it. If that worktree is live it inherits the same change.

---

## A. The type — `IterationHistory`

**Definition:** `orpheus/sn/solution.py:111-186` (decorator at `:111`, `class` at
`:112`).

```python
@dataclass(frozen=True)
class IterationHistory:
```

**Frozen:** yes (`@dataclass(frozen=True)`). Tuple-valued trajectories by design
(the docstring cites `coding-elegance` Pattern 4 — a frozen dataclass with mutable
fields is the anti-pattern).

### Field list (`orpheus/sn/solution.py:158-163`)

| field | type | default | meaning (from the docstring at `:120-148`) |
|---|---|---|---|
| `keff_history` | `tuple[float, ...]` | `()` | per-outer eigenvalue trajectory; empty for fixed-source |
| `flux_residuals` | `tuple[float, ...]` | `()` | per-iteration relative flux residual `‖φ_n − φ_{n−1}‖/‖φ_n‖`; empty when untracked |
| `n_inner` | `int \| None` | `None` | inner (within-group) iterations consumed; `None` on the pure-eigenvalue outer path |
| `total_inner_iterations` | `int \| None` | `None` | inner iterations summed across ALL outers; populated by BOTH paths; the measurand for the SI-rate diagnostics |
| `n_outer` | `int \| None` | `None` | outer (power) iterations consumed; `None` for fixed-source |
| `converged` | `bool` | **`True`** | "`True` when the iteration met its convergence tolerance, `False` when `max_iter` was exhausted" |

⛔ **The default is `True`.** `IterationHistory()` — no arguments — asserts
convergence. That is the same lie as #342 but at the *type* level: the type
does not make "I made no claim" representable, and the safe-by-omission
posture is the *unsafe* one. Any construction site that forgets the kwarg
silently claims success. (Two of the five production sites do exactly this —
see §B.)

### Methods / properties

Three, all read-only accessors; **none** of them derives `converged`:

| member | `file:line` | returns |
|---|---|---|
| `dominance_ratio()` | `orpheus/sn/solution.py:165-178` | `\|k_n − k_{n−1}\|/\|k_{n−1}\|`, or `None` if `<2` entries |
| `latest_keff()` | `orpheus/sn/solution.py:180-182` | last eigenvalue or `None` |
| `latest_residual()` | `orpheus/sn/solution.py:184-186` | last recorded flux residual or `None` |

`latest_residual()` is the ingredient a derived `converged` would need — it
already exists and is already the honest signal; nothing consumes it for the
convergence decision today.

### The delegate on the solution carrier

`SolutionBase.converged()` — `orpheus/sn/solution.py:380-382`:

```python
def converged(self) -> bool:
    """Return whether the solver iteration met its convergence tolerance."""
    return self.history.converged if self.history else True
```

⛔ **A SECOND default-True.** A `Solution` with `history=None` reports
`converged() == True`. `history` is `Optional` on the carrier
(`orpheus/sn/solution.py:268`), so "no diagnostics at all" and "converged" are
the same answer. This is a latent lie of the #342 family (see §"second lies").

Note the **name collision**: `SolutionBase.converged` is a zero-arg METHOD, while
`IterationHistory.converged` is a FIELD. `sol.converged` (no parens) is a truthy
bound method — an `assert sol.converged` would pass unconditionally. Grep found
no such site, but the shape is a live footgun and matters if the remedy renames
or re-types either one.

### Public re-export

Yes. `orpheus/sn/__init__.py:8` re-exports `IterationHistory` (alongside
`AdjointSolution`, `Solution`, `SolutionBase`, `SolutionDiff`); it is in
`orpheus/sn/solution.py:102-108`'s `__all__`. **It is public API.**

### Documentation

| page | `file:line` | what it says |
|---|---|---|
| `docs/theory/conventions/indexing_and_layout.rst` | `:919-922` | names the dataclass and lists `converged` as one of the "accessors" — `:attr:`~orpheus.sn.solution.IterationHistory.converged`` |
| same | `:944`, `:955` | `keff_history` / `flux_residuals` layout rows |
| same | `:986` | "`IterationHistory` carrying tuple-based …" |
| `docs/theory/methods/sn/cartesian_multid.rst` | `:4118-4119` | cites `IterationHistory.total_inner_iterations` as the SI-rate measurand |

**No page documents what `converged` MEANS per entry point** — i.e. nothing
states that fixed-source's flag is inner-based, `solve_sn`'s is unconditional,
and the adjoint's is outer-based. The prose treats it as one concept. That is
the doc-side of the same defect.

---

## B. Every CONSTRUCTION site of a convergence claim

### B.1 `IterationHistory` — the SN production type (5 sites, ALL in `orpheus/sn/solver.py`)

**The brief said "at least three spellings". There are FIVE, and they are five
DIFFERENT predicates** — the fourth and fifth differ from the ones already
known by strictness and by which residual they read.

| # | `file:line` | entry | exact `converged=` expression | honest? |
|---|---|---|---|---|
| 1 | `orpheus/sn/solver.py:2264-2269` | `solve_sn` — **forward eigenvalue** | `True` (literal) | ⛔ **NO — issue #342** |
| 2 | `orpheus/sn/solver.py:2533-2537` | `solve_sn_adjoint` — **adjoint eigenvalue** | `len(keff_history) < max_outer` | yes — OUTER-based |
| 3 | `orpheus/sn/solver.py:2676-2681` | `solve_sn_adjoint_fixed_source` | `bool(residuals) and residuals[-1] <= inner_tol` | yes — INNER-based, **`<=`** |
| 4 | `orpheus/sn/solver.py:3288` → `:3298-3303` | `_solve_fixed_source_si` (`solve_sn_fixed_source`, SI arm) | `bool(residuals) and residuals[-1] < inner_tol` | yes — INNER-based, **`<`** |
| 5 | `orpheus/sn/solver.py:3512` → `:3521-3526` | `_solve_fixed_source_krylov` (`solve_sn_fixed_source`, Krylov arm) | `bool(residuals) and residuals[-1] < inner_tol` | yes — INNER-based, **`<`** |

⚠ **New finding — a 4th spelling hiding as a duplicate of the 3rd.** Site 3 uses
`residuals[-1] <= inner_tol`; sites 4 and 5 use `residuals[-1] < inner_tol`. The
adjoint fixed-source is therefore *strictly more permissive* than its forward
sibling on the exact-equality boundary. This also **disagrees with the
certificate's own gate**, `_certify_within_group_exit`
(`orpheus/sn/solver.py:428`), which uses `not (residual_history[-1] < tol)` —
i.e. STRICT. So on an exactly-equal residual the adjoint fixed-source claims
convergence while the certificate declines to certify it. Not reachable in
practice (FP equality), but it is a fourth transcription of one predicate.

⚠ **Sites 2 and 3 do not populate `total_inner_iterations`** (site 2) or
`keff_history`/`n_outer` (site 3) — those fall to their `None`/`()` defaults.
Site 2 also omits `flux_residuals`. A derived-`converged` design has to decide
what it means when the ingredient is absent.

**Zero of the five uses a shared helper.** The predicate is transcribed five
times.

### B.2 The root cause is one hop UP — `power_iteration` discards the answer

`orpheus/numerics/eigenvalue.py:205-274`:

```python
for n in range(1, max_iter + 1):
    ...
    if solver.converged(keff, keff_old, flux_distribution, flux_old, n):
        break
return keff, keff_history, flux_distribution        # <-- 3-tuple, NO flag
```

The loop KNOWS whether it broke or exhausted; the return signature throws it
away. Every caller must re-derive it from `len(keff_history) < max_iter`.
`solve_sn_adjoint` does; `solve_sn` does not. **#342 is not a typo — it is the
predictable consequence of a primitive that returns a 3-tuple instead of an
outcome.** This is the single most load-bearing structural fact in this audit.

`KEigenvalue.solve` (`orpheus/numerics/iteration.py:1349-1400`) forwards the
same 3-tuple verbatim (`return power_iteration(self, max_iter=self.max_outer)`),
so the loss propagates through the operator-triple entry too.

**`power_iteration` is consumed by FOUR production solvers** — it IS the shared
cross-solver primitive, already existing:

| consumer | `file:line` |
|---|---|
| SN forward eigenvalue | `orpheus/sn/solver.py:2131` |
| CP | `orpheus/cp/solver.py:882` |
| MoC | `orpheus/moc/solver.py:137` |
| diffusion | `orpheus/diffusion/solver.py:407` |
| `KEigenvalue` (the operator-triple entry; SN adjoint routes here) | `orpheus/numerics/iteration.py:1400` |

### B.3 Cross-solver: CP / MoC / MC / diffusion / homogeneous carry **NO** convergence claim at all

This is the answer to the brief's Pattern-2 question, and it is **not** "each
rolls their own":

| solver | result type | `file:line` | has `converged`? | what it does carry |
|---|---|---|---|---|
| CP | `CPResult` | `orpheus/cp/solver.py:73` | **NO** | `keff_history`, `residual_history`, `n_inner`, `elapsed_seconds` |
| MoC | `MoCResult` | `orpheus/moc/solver.py:36` | **NO** | `keff_history`, `elapsed_seconds` |
| MC | `MCResult` | `orpheus/mc/solver.py:526` | **NO** | `keff_history`, `sigma_history` (statistical, not iterative) |
| diffusion | `DiffusionResult` | `orpheus/diffusion/solver.py:336` | **NO** | `keff_history` |
| homogeneous | `HomogeneousResult` | `orpheus/homogeneous/solver.py:64` | **NO** (direct solve — no iteration) | — |
| SN | `IterationHistory` on `Solution` | `orpheus/sn/solution.py:112` | yes (5 spellings) | the full trajectory |

⟹ **SN is the only production solver that reports convergence at all.** The
other three iterative ones (CP, MoC, diffusion) run `power_iteration` and drop
the outcome on the floor entirely — a caller cannot even ASK. So the smell is
not duplication; it is a **missing shared return type**, with SN having
independently grown a partial local one.

They DO all implement the outer PREDICATE — the `EigenvalueSolver.converged`
Protocol method (`orpheus/numerics/eigenvalue.py:159-168`), five
implementations:

| implementer | `file:line` | note |
|---|---|---|
| `SNSolver.converged` | `orpheus/sn/solver.py:1476-1486` | `iteration <= 2 → False`, then `dk < keff_tol and dφ < flux_tol` |
| `KEigenvalue.converged` | `orpheus/numerics/iteration.py:1326-1347` | same shape, carrier-honest norms |
| `CPSolver.converged` | `orpheus/cp/solver.py:722` | — |
| `DiffusionSolver.converged` | `orpheus/diffusion/solver.py:316` | — |
| MoC | `orpheus/moc/core.py:242` | — |

⚠ **`iteration <= 2 → False` is structural**: any solve with `max_outer <= 2`
can NEVER be reported converged by a derived predicate. One live test does
exactly that (§D).

⚠ **Naming collision to resolve before designing.** `converged` is a **field**
on `IterationHistory`, a **zero-arg method** on `SolutionBase`
(`orpheus/sn/solution.py:380`), and a **5-arg predicate** on the
`EigenvalueSolver` Protocol. Three different things, one word.

### B.4 `orpheus/derivations/` — a third, independent family

Not production, but it is where the honest design already exists.

| type | `file:line` | `converged` | note |
|---|---|---|---|
| `CriticalSolution` | `orpheus/derivations/common/solution_types.py:126` (field at `:195`) | `bool = True` | ⛔ same default-True as `IterationHistory` |
| `PowerIterationResult` | `orpheus/derivations/continuous/trajectory_resolvent/power_iteration.py:78` (field at `:107`) | `bool` — **no default** | ✅ the honest one; set from a real loop flag at `:189/:214/:221` |
| `SlabFluxResult`-family | `…/fn_method/slab/flux_reconstruction.py:223`, `…/slab/reflected.py:119`, `…/sphere/flux_reconstruction.py:85` | `bool` no default | honest, derived from a loop flag |
| `…/peierls_nystrom/ps1982_reference.py:99` | | `bool` no default | honest |
| `…/singular_eigenfunction/cylinder/one_group.py:202` | | `bool` no default | honest (`bool(brent_result.converged)`) |
| `greens_function*` family (7 modules) | e.g. `…/trajectory_resolvent/greens_function.py:375` | derived from `pi_result.converged` | honest — propagated, not re-derived |

**Hardcoded `converged=True` in derivations — 8 sites, all `CriticalSolution`:**

| `file:line` | justified? |
|---|---|
| `…/fn_method/moment_space.py:315` | ✅ inline comment: "bisection always converges given bracket" |
| `…/fn_method/moment_space.py:347`, `:424` | ⚠ no stated reason — same family, probably the same bisection |
| `…/galerkin_spectral/basis_space.py:756`, `:793` | ⚠ wraps `solve_galerkin_spectral_{slab,sphere}` — a dense eigensolve; defensible (direct method, no iteration) but **unstated** |
| `…/singular_eigenfunction/spectrum.py:819`, `:872` | ⚠ wraps `solve_case_method_{slab,sphere}_critical` — a root-find; the underlying result carries `eq46_residual` in metadata but the flag is a literal |

`…/singular_eigenfunction/spectrum.py:924` is the honest sibling in the SAME
class — `converged=bool(res.converged)`. So `:819`/`:872` are the strongest
"second latent lie" candidates outside SN (see the closing section).

---

## C. Every READ of `.converged`

Split by **type family**, because they are three unrelated types that share a
word. The blast radius of an SN semantics change is the FIRST table only.

### C.1 The SN family (`IterationHistory.converged` / `SolutionBase.converged()`) — 15 reads

| `file:line` | expression | BRANCH or REPORT? | what it decides |
|---|---|---|---|
| `orpheus/sn/solution.py:382` | `self.history.converged if self.history else True` | **BRANCH** (on `history is None`) | the ONLY production read in the tree. Returns `True` for a history-less solution — §"second lies" |
| `tests/sn/solve/test_coupled_solve_certificate.py:102` | `if not bool(solution.history.converged): pytest.fail(...)` | **BRANCH** | fixture-drift guard on the CONTROL leg of the lag-death classifier |
| `tests/sn/solve/test_coupled_solve_certificate.py:125` | same | **BRANCH** | proves the monkeypatch reverted cleanly |
| `tests/sn/solve/test_fixed_source_g1.py:258` | `assert result.history.converged, "…did not converge in 200 iters"` | **BRANCH** | precondition for the G1 per-ordinate fixed-point assertion |
| `tests/sn/solve/test_fixed_source_2d_equivalence.py:96` | `assert result.history.converged, "2-D fixed-source Krylov did not converge"` | **BRANCH** | precondition |
| `tests/sn/solve/test_affine_carve_bit_identity.py:218` | `if not sol.history.converged: raise AssertionError(...)` | **BRANCH** | precondition before a sha256 bit-identity compare |
| `tests/sn/solve/test_d3_admission.py:186` | `np.testing.assert_equal(bool(sol.history.converged), True)` | **BRANCH** | ⚠ landed TODAY (`11e78430`) as the #342-adjacent fix |
| `tests/sn/solve/test_d3_admission.py:229` | same | **BRANCH** | same commit |
| `tests/sn/solve/test_sn_adjoint_entries.py:200-204` | compound `assert adj.history is not None and … and adj.history.converged` | **BRANCH** | the only test reading the ADJOINT eigenvalue spelling (site 2) |
| `tests/sn/operators/test_psi_half_coupling.py:3767` | `if not sol.history.converged: pytest.fail(...)` | **BRANCH** | precondition on a pure-absorber sphere closed form |
| `tests/sn/verification/analytical/test_si_convergence_rate.py:162` | `assert sol.history.converged, "SI did not converge — count is meaningless"` | **BRANCH** | ⭐ the ONE place the project already articulates the principle |
| `tests/sn/primitives/test_solution.py:95` | `assert h.converged is True` | **BRANCH** | ⛔ **pins the default** — `IterationHistory()` must report `True`. Any type-level fix reds this line |
| `tests/sn/primitives/test_solution.py:276` | `assert sol.converged() is True` | **BRANCH** | pins the `history=None → True` delegate |
| `tests/sn/primitives/test_solution.py:291-292` | `assert sol_yes.converged() is True` / `sol_no.converged() is False` | **BRANCH** | the delegate round-trip |
| `tests/sn/solve/test_d3_admission.py:144` | docstring prose "returned an honest `history.converged = False`" | REPORT (prose) | — |

**SN production branches on the flag: exactly ONE** (`solution.py:382`, and it
branches on the history's *presence*, not on the flag's value). **Nothing in
`orpheus/` outside `solution.py` reads it at all.** No retry, no raise, no warn,
no log. The flag is written five times and consumed only by tests.

⚠ **`tests/sn/primitives/test_solution.py:95` and `:276` are the two lines a
type-level fix must plan for.** They are not incidental — they are the
characterization of the default-`True` behaviour. If `converged` loses its
default (or becomes an outcome enum), those assertions are the contract change
made visible. Per `coding-standards` retirement/marker migration, they should be
**rewritten to the new contract**, not deleted.

### C.2 The `orpheus/derivations/` family — a different type, 18 + 87 reads

- **18** reads inside `orpheus/derivations/` — of which **15 are pure
  propagation** (`converged=pi_result.converged`, `converged=res_asym.converged`,
  `converged=bool(res.converged)`) → REPORT.
- **3** `if converged:` sites — `…/fn_method/slab/reflected.py:956`,
  `…/peierls_nystrom/geometry.py:6453`, `…/peierls_nystrom/slab.py:600` — these
  read a **local loop variable**, not a result object's field, so they are
  loop-exit logic, not consumers. Not blast radius.
- **87** `assert res.converged` in `tests/derivations/` → **BRANCH**, but on the
  derivations result types, which SN does not touch. Untouched by an SN change;
  **in scope** if the remedy extracts a genuinely shared primitive.
- **4** reads in `tests/cross_method/adapters.py:196, 270, 340, 404` —
  `"converged": bool(res.converged)` into a comparison dict → **REPORT**. This is
  the cross-method protocol surface; it consumes the DERIVATIONS
  `CriticalSolution.converged`, not SN's.

### C.3 The `EigenvalueSolver.converged(...)` PREDICATE — not the same thing

`tests/moc/test_verification.py:1171-1172` calls
`solver.converged(1.0, 1.0, phi, phi, iteration=1) is False` — testing the
5-argument Protocol method. Same word, different concept. Not blast radius.

### C.4 Docs prose describing the flag

| `file:line` | claim | status |
|---|---|---|
| `docs/theory/conventions/indexing_and_layout.rst:919-923` | "`IterationHistory` … accessors such as `dominance_ratio` / `converged`", "which a `Solution` carries as an optional field (**populated for eigenvalue problems**)" | ⛔ **the parenthetical is FALSE** — all four SN entries populate `history`, including both fixed-source arms. Also `converged` is tagged `:attr:` while `dominance_ratio` is a method — the role mix-up mirrors the code's own field-vs-method collision |
| `docs/theory/conventions/indexing_and_layout.rst:986-997` | the carrier bullet list, `SolutionBase.dominance_ratio` / `SolutionBase.converged` — "iteration diagnostics that read as math" | accurate but says nothing about semantics |
| the "Iteration state" `list-table` (`:928-958`) | rows for `keff`, `keff_history`, `Eigenpair`, `ResidualHistory`, `DominanceRatio` | ⚠ **`converged` has NO row** — the flag is the one iteration-state quantity the layout table omits |

**No page anywhere states what `converged` means per entry**, that the forward
eigenvalue value is unconditional, or that a `max_inner` truncation returns
silently. The `coding-standards` "grep the CONCEPT, not the symbol" sweep for
`best-effort` / `truncat` / `budget` in `docs/theory/` returns nothing on this
topic either — the contract is undocumented, not mis-documented.

---

## D. Deliberate non-convergence consumers (the opt-out population)

⚠ **The brief's example needs a correction.** `tests/sn/acceleration/test_dsa_rate.py`
does NOT cap at 50 on purpose in the sense meant: its helper
`_uniform_solve_raw` defaults to **`max_inner=4000`**
(`tests/sn/acceleration/test_dsa_rate.py:100`), and the `max_inner=50` sites
(`:552`, `:579`, `:672`, `:697`) are **headroom** on tests that assert exact
landing in 2–3 iterations (`_inners(sol) == n_expected and res[-1] < 1e-13`).
Those solves DO converge. The one genuine truncation consumer in that module is
a different test (row 1 below).

### D.1 Genuine opt-outs — a solve that may legitimately not converge

| `file:line` | what it does | budget | why it can't converge / doesn't care |
|---|---|---|---|
| `tests/sn/acceleration/test_dsa_rate.py:503-505` | `test_accelerated_decade_count_bounded_on_grid` — harvests the residual history, counts iterations to a 10-decade reduction | `max_inner=200`, `inner_tol=1e-11` | measures a RATE; the reflective `σ_t·h=100` column is documented at 63–84 decades-count and may not reach `1e-11` in 200. **Harvests `_residuals(sol)`, never the answer.** |
| `tests/sn/solve/test_krylov_curvilinear_precond_safety.py:449-453` | spy on the GMRES `restart` argument; the `solve_sn` return value is **discarded** | `max_outer=2` | ⛔ **structurally impossible to converge**: `SNSolver.converged` returns `False` for `iteration <= 2` (`orpheus/sn/solver.py:1481`). Any derived outer flag is `False` by construction |
| `tests/sn/solve/test_fixed_source_g1.py:176-182` | spy on `KrylovAcceleration.solve` call count + input immutability; return discarded | `max_inner=10`, `inner_tol=1e-6` | the answer is irrelevant — the test asserts the spy fired once and `external_source` was unmodified |
| `tests/sn/sweep/curvilinear/test_si_cyl_20cell_nan_regression.py:113-119` | asserts `np.isfinite(res.keff)` on the NaN-regression config | `max_outer=3, max_inner=3` | inline comment: *"NaN appears in the FIRST inner iteration — small caps suffice"* |
| `tests/sn/solve/test_coupled_solve_certificate.py:128-147` | `test_certificate_is_a_noop_without_a_convergence_claim` — calls `_certify_within_group_exit` directly with `residual_history=[]` and `[1.0]` | n/a | ⭐ **this test IS the contract that a best-effort exit is legal.** It is the project's existing, deliberate statement of the current policy — the remedy must decide whether to keep, re-scope, or retire it (issue #340) |
| `tests/sn/solve/test_si_single_primitive_contract.py:130, 149` | `SNSolver(sn_mesh, max_inner=4)` — a structural/primitive-identity contract | `max_inner=4` | asserts WHICH primitive is used, not the answer |
| `tests/sn/architecture/test_stage_separation.py:129, 167, 198` | `SourceIteration(..., max_iter=2)` / Krylov `max_iter=2` | 2 | stage-separation structure, not numerics. Below the public entries — unaffected by an entry-level policy |
| `tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py:126` | `max_inner=5` | 5 | ⚠ **UNCERTAIN** — did not read the assertion; flagging for the implementer |
| `tests/sn/verification/mms/test_curvilinear_operator_admits_mms.py:72` | `max_inner=10, inner_tol=1e-13` | 10 | ⚠ **UNCERTAIN** — a `1e-13` tol in 10 iterations is very unlikely to be met; probably an admission/smoke check, but I did not read the body |
| `tests/sn/sweep/core/test_cache.py:245-246` | `max_outer=50, max_inner=50` | 50/50 | ⚠ **UNCERTAIN** — a cache-behaviour test; probably indifferent |
| `tests/cp/test_verification.py:1100` | CP `max_inner=50` | 50 | CP has no `converged` field — unaffected today; in scope only if the primitive is extracted |

### D.2 NOT opt-outs (do not count these)

| `file:line` | why not |
|---|---|
| `tests/sn/acceleration/test_dsa_rate.py:552, 579, 672, 697` | `max_inner=50` is headroom; the tests assert `n == 2/3` and `res[-1] < 1e-13`, i.e. they REQUIRE convergence |
| `tests/sn/sweep/curvilinear/test_coupled_pole_mu_level_invariant.py:238-242` | wrapped in `pytest.raises(ValueError, match="non-carrying")` — the solve raises before iterating |
| `tests/sn/verification/analytical/test_si_convergence_rate.py` (all) | uses `max_inner=20000` / `5000` and ASSERTS `history.converged` — the model consumer |
| `tests/sn/solve/test_krylov_curvilinear_precond_safety.py:497` (`max_inner=40, inner_tol=1e-6`) | also a spy test, but on `solve_sn_fixed_source`; return discarded — same class as D.1 row 3 |
| every `tests/derivations/*` site | different type family |

### D.3 The size of the opt-out population

**Roughly 8–11 test sites**, of which **2 discard the return value entirely**
(spy tests), **2 harvest the residual history**, **1 asserts only finiteness**,
**1 is the policy contract itself**, and **3 are uncertain**. **Zero production
opt-outs. Zero `examples/` opt-outs** (`examples/discrete_ordinates/demo_discrete_ordinates.py:84`
calls `solve_sn` with all defaults; the other four demos call CP / MoC /
diffusion, which carry no flag).

⟹ the population is small enough that **raise-by-default with an explicit
opt-out is affordable**, and a `warn`-by-default is affordable with zero test
churn if the warning is not escalated to an error in the project's pytest
config. The two spy tests and the `max_outer=2` site are the ones that would need
the opt-out kwarg regardless.

### D.4 Three uncertainties from D.1 — RESOLVED

| site | verdict |
|---|---|
| `tests/sn/verification/mms/test_curvilinear_aniso_scattering_p1.py:126` | **NOT an opt-out.** Constructs `SNSolver(..., max_inner=5)` and only reads `solver.scattering_op` — never calls a solve. The budget is INERT |
| `tests/sn/verification/mms/test_curvilinear_operator_admits_mms.py:72` | **NOT an opt-out.** Same shape — `SNSolver(..., max_inner=10)`, reads `solver.sn_mesh` / `solver.scattering_op` / `solver.mat_xs` for an operator-residual probe; no solve |
| `tests/sn/sweep/core/test_cache.py:239-250` | **NOT a deliberate opt-out** — it runs a real `solve_sn` expecting ~7 outers and asserts `len(keff_history) >= 5`. But it caps `max_inner=50` at `inner_tol=1e-8`, so under a derived INNER-aware flag it is the one site whose status must be MEASURED, not read |

⟹ the genuine opt-out population shrinks to **6 sites** (D.1 rows 1–6 plus
`test_krylov_curvilinear_precond_safety.py:497`), of which one is the policy
contract itself.

### D.5 The pytest warning policy — decisive for `warn` vs `raise`

`pyproject.toml:69-93` (`[tool.pytest.ini_options]`) sets `testpaths`,
`pythonpath`, `addopts = "--import-mode=importlib"` and the marker list. **There
is NO `filterwarnings` key, and no `-W error` in `addopts`.**

⟹ **a `RuntimeWarning`-by-default policy would cause ZERO test churn today.** The
only project-wide warning escalation that exists is opt-in and scoped:
`tests/sn/regression/conftest.py:26-31` registers `always::DriftWarning`
(surface, not error), and `-W error::DriftWarning` is documented as a MANUAL
invocation. A `raise`-by-default policy costs the ~6 opt-out sites an explicit
kwarg.

---

## E. The public entry surface

### E.1 SN — four entries, all re-exported from `orpheus/sn/__init__.py:13-19`

| entry | `file:line` | `max_outer` | `keff_tol` | `flux_tol` | `max_inner` | `inner_tol` | `inner_schedule` | history site |
|---|---|---|---|---|---|---|---|---|
| `solve_sn` (forward eigenvalue) | `orpheus/sn/solver.py:2025-2037` | **500** | 1e-7 | 1e-6 | **200** | **1e-8** | **`"jacobi"`** | B.1 #1 ⛔ |
| `solve_sn_adjoint` (adjoint eigenvalue) | `orpheus/sn/solver.py:2442-2453` | **500** | 1e-7 | 1e-6 | **200** | **1e-8** | — (not a parameter) | B.1 #2 |
| `solve_sn_fixed_source` | `orpheus/sn/solver.py:2974-2988` | — | — | — | **1000** | **1e-12** | **`"gauss_seidel"`** | B.1 #4 / #5 |
| `solve_sn_adjoint_fixed_source` | `orpheus/sn/solver.py:2541-2552` | — | — | — | **1000** | **1e-12** | — (not a parameter) | B.1 #3 |

**Asymmetries (beyond the `max_inner` 200-vs-1000 the brief already knew):**

1. **`inner_tol` differs by FOUR ORDERS** — `1e-8` on the eigenvalue entries,
   `1e-12` on the fixed-source ones. So the eigenvalue entries have both a
   tighter budget AND a looser tolerance; the pairing is coherent (inexact-Newton
   inner) but nothing states it.
2. **`inner_schedule` default differs** — `solve_sn` defaults `"jacobi"`,
   `solve_sn_fixed_source` defaults `"gauss_seidel"`. (The d3 diagnosis §7
   measured G-S ~2× SLOWER than Jacobi at d=3 on reflective configs, so the
   fixed-source default is the slow arm exactly where the budget bites.)
3. **The adjoint entries expose no `inner_schedule` / `inner_solver` /
   `acceleration` at all** — they are narrower surfaces than their forward twins.
4. **`solve_sn_fixed_source` alone exposes `acceleration`** (`:2987`), so it is
   the only entry whose iteration count is user-alterable by preconditioning.
5. `SNSolver.__init__` (`orpheus/sn/solver.py:924`) carries its own
   `max_inner: int = 200` — a **fifth** default, reached by any consumer that
   constructs the solver directly instead of calling an entry.

### E.2 Non-SN public entries — none carries a convergence claim

| entry | `file:line` | budget default | returns | claim? |
|---|---|---|---|---|
| `solve_cp` | `orpheus/cp/solver.py:792` | `CPParams.max_outer=500` (`:58`), `max_inner=100` (`:69`) | `CPResult` | **none** |
| `solve_moc` | `orpheus/moc/solver.py:69-75` | `max_outer=500` | `MoCResult` | **none** |
| `solve_diffusion_1d` | `orpheus/diffusion/solver.py:355-361` | `max_outer=1000` | `DiffusionResult` | **none** |
| `solve_monte_carlo` | `orpheus/mc/solver.py:550` | n/a (batches) | `MCResult` | **none** (statistical σ instead) |
| `solve_homogeneous_infinite` | `orpheus/homogeneous/solver.py:150` | n/a (direct) | `HomogeneousResult` | n/a |

⚠ `solve_diffusion_1d` defaults `max_outer=1000` while `solve_sn` defaults
`500` — the two solvers that share `power_iteration` disagree on the cap by 2×,
with no stated reason.

### E.3 The below-the-entry primitives that own the truncation

| primitive | `file:line` | budget default | on exhaustion |
|---|---|---|---|
| `power_iteration` | `orpheus/numerics/eigenvalue.py:205-207` | `max_iter=500` | **silent** — returns the 3-tuple |
| `SourceIteration.solve` | `orpheus/numerics/iteration.py:632-748` (cap at `:593`) | `max_iter=1000` | **silent** — `return psi, residual_history` |
| `KrylovAcceleration.solve` | `orpheus/numerics/iteration.py:900-1041` (cap at `:864`) | `max_iter=1000` | ⭐ **`warnings.warn(..., RuntimeWarning)`** at `:1027-1039` |
| `KEigenvalue` | `orpheus/numerics/iteration.py:1144-1147, 1400` | `max_outer=500`, `max_inner=1000` | **silent** (forwards `power_iteration`) |
| `GreenOperator.apply` | `orpheus/numerics/green_operator.py:300-345` (cap at `:253`) | `max_iter=1000` | ⭐ **`raise ConvergenceFailure`** at `:337` |
| `rayleigh_quotient_iteration` | `orpheus/numerics/eigenvalue.py:473-612` | `max_iter=50` | ⭐ **`warnings.warn(..., RuntimeWarning)`** at `:600-609` |

⭐ **Three of the six already implement a policy, and they disagree** — one
raises, two warn, three are silent. That is the *same* defect one layer down,
and it is the layer where the fix belongs.

---

## F. Existing vocabulary to reuse

### F.1 `ConvergenceCertificateError` — NOT the right base for a budget error

**Definition:** `orpheus/sn/solver.py:356-372` (`class ConvergenceCertificateError(RuntimeError)`).

**What it means TODAY (narrow, and deliberately so):** *"a within-group solve
CLAIMED convergence but the honest equation residual disagrees"* — the #282
lag-death classifier. It fires only when the driver's free-identity stop
`r = rhs_{n−1} − rhs_n` said "converged" while the true `‖Aψ − q‖/‖q‖` exceeds
`_CERTIFICATE_SAFETY (=10.0) × tol`.

- **Raised at:** `orpheus/sn/solver.py:435` — one site, inside
  `_certify_within_group_exit` (`orpheus/sn/solver.py:384-445`).
- **Called from:** the un-windowed SI arm (`orpheus/sn/solver.py:3271-3278`) and
  the Krylov arms; two documented structural exemptions (windowed 2-D moment
  arm; moment-tailed LD schemes, gated at `:426`).
- **The claim gate — the line issue #340 is about:** `orpheus/sn/solver.py:427-428`
  `if not residual_history or not (residual_history[-1] < tol): return`.
- **Caught by:** `tests/sn/solve/test_coupled_solve_certificate.py:114`,
  `tests/sn/operators/test_psi_half_coupling.py:2245` (both
  `pytest.raises(..., match="lag-death")`).
- **Referenced in prose:** `orpheus/numerics/operator.py:467`,
  `tests/sn/solve/test_declared_inflow_reaches_the_rhs.py:361`,
  `tests/sn/operators/test_declared_law_is_linear.py:23`,
  `tests/sn/verification/analytical/test_mms_declared_inflow.py:466`, and four
  theory pages (`docs/theory/methods/sn/history.rst:193`,
  `docs/theory/foundations/boundary_conditions.rst:170, 2305`,
  `docs/theory/foundations/coupled_block_operator.rst:584`,
  `docs/theory/foundations/operator_algebra.rst:350`).

⟹ **Do not widen it.** "The claim was false" and "no claim was made" are
different propositions, and the existing `match="lag-death"` pins plus five
theory pages make the current meaning load-bearing. A **sibling** is the right
move — the d3 diagnosis's suggested `InnerIterationBudgetExceeded` fits, and it
should live beside `ConvergenceCertificateError` OR (better, see scope) in
`orpheus/numerics/`.

### F.2 `ConvergenceFailure` — the closer precedent, already in `orpheus/numerics/`

`orpheus/numerics/green_operator.py:127` — `class ConvergenceFailure(RuntimeError)`,
exported in that module's `__all__` (`:122`).

Meaning: *"the splitting iteration did not reach the promise"*, raised at
`orpheus/numerics/green_operator.py:337-345` when the TRUE relative residual is
still above `tol` after the budget, or is non-finite. This is **exactly the
proposition** the SN entries need — and it is already in the shared numerics
layer, already documented (`docs/theory/foundations/operator_inverse_family.rst:625, 655`,
`docs/theory/methods/sn/history.rst:766-774`), already tested with teeth
(`tests/numerics/test_green_operator.py:236, 271, 323`,
`tests/sn/operators/test_green_operator_sn.py:252`).

⚠ Its docstring/message vocabulary is Green-operator-specific ("the sum's
SPELLING chose the wrong preconditioner"), so promoting it verbatim would carry
the wrong provenance. But the **name and the RuntimeError base are the
established vocabulary** and a `numerics`-level home for the family already
exists.

### F.3 The established WARNING vocabulary

There is **no** dedicated convergence warning class. The precedents both use
plain `RuntimeWarning`:

| site | message shape |
|---|---|
| `orpheus/numerics/iteration.py:1027-1039` | `"KrylovAcceleration.solve: scipy…gmres returned info={info} (not converged within maxiter=…); Returning best-effort iterate; residual_history tail = …"` + `RuntimeWarning`, `stacklevel=2` |
| `orpheus/numerics/eigenvalue.py:600-609` | `"rayleigh_quotient_iteration: no convergence in max_iter=… (last ‖Δv‖=…, tol=…). … Returning the best iterate."` + `RuntimeWarning`, `stacklevel=2` |

⭐ **The Krylov one is the model, and its comment is the design rationale
already written down** (`orpheus/numerics/iteration.py:1001-1012`):

> *"Pre-fix the `info` flag was discarded; an unconverged `solution` would
> silently be consumed as if it were the true inverse. … Both surface as
> warnings — raising would break long-standing callers that tolerate slow
> convergence and need the best-effort iterate. See ERR-053."*

That is the SAME defect (a discarded convergence flag), the SAME reasoning, and
an explicit prior decision for **warn, not raise**, taken in D-H.1e. It also
carries an **exact-breakdown carve-out** (`:1026-1027`) — a precedent for
"structurally-converged despite `info > 0`", which the remedy will need for the
`iteration <= 2` case.

### F.4 `DriftWarning` — NOT this family

`tests/sn/regression/_regression_assert.py:69` — `class DriftWarning(UserWarning)`.
It lives in **`tests/`**, not `orpheus/`, and it is the bit-identity tripwire for
regression snapshots (registered in `tests/sn/regression/conftest.py:22-31`;
escalatable with `-W error::DriftWarning`). Unrelated proposition. Its value here
is only as a **pattern**: a named `Warning` subclass + a conftest registration +
a documented `-W error::<name>` escalation is the project's existing recipe for
"loud but opt-in-strict". A `ConvergenceWarning(RuntimeWarning)` in
`orpheus/numerics/` would follow that recipe and let CI escalate it globally
without touching a single test today (§D.5).

### F.5 Summary of the vocabulary to extend

| proposition | existing home | verdict |
|---|---|---|
| "the claim was false" (in-M lag) | `ConvergenceCertificateError`, `orpheus/sn/solver.py:356` | KEEP as-is; do not widen |
| "the budget ran out, answer is best-effort" | `ConvergenceFailure`, `orpheus/numerics/green_operator.py:127` | **reuse the name/base; move or re-home to a shared `orpheus/numerics/` convergence module** |
| "the budget ran out — warn" | plain `RuntimeWarning` ×2 | **mint `ConvergenceWarning(RuntimeWarning)`** in the same place, so the two existing sites can adopt it and CI can escalate one category |
| "the outcome of an iteration" | ⛔ **does not exist** | the actual gap — see scope |

---

## ⚠ DRIFT NOTICE — `orpheus/sn/solver.py` was edited DURING this audit

`git status` was clean on `orpheus/` at dispatch open and shows
`M orpheus/sn/solver.py` (+42/−4) and
`M tests/sn/verification/analytical/test_si_convergence_rate.py` at close
(lessons L-012). The edit is a **docstring-only correction about the
boundary-Gauss-Seidel rate regime (issue #341)** — the d3 diagnosis's §7
secondary finding being written up. **It does not touch the convergence
contract**; the five `converged=` expressions are byte-identical.

It DID shift the fixed-source line numbers by +38. **Corrected addresses (final
read, tree moving):**

| item | address at open | address at close |
|---|---|---|
| `_solve_fixed_source_si` def | `:3153` | **`:3191`** |
| SI `converged_flag` | `:3288` | **`:3328`** |
| SI `IterationHistory(` | `:3298` | **`:3336-3341`** (`converged=` at `:3340`) |
| `_solve_fixed_source_krylov` def | `:3364` | **`:3402`** |
| Krylov `converged_flag` | `:3512` | **`:3550`** |
| Krylov `IterationHistory(` | `:3521` | **`:3559-3564`** (`converged=` at `:3563`) |

Unmoved: `:2264-2269`, `:2533-2537`, `:2676-2681`, `:427-428`, `:435`, and every
signature line in §E.1.

Also note issue **#341** now exists (filed today, boundary-G-S rate regime) — a
sibling of #340/#342 in the same diagnosis, worth cross-referencing.

---

## ⚠ A STRUCTURAL CONSTRAINT the design must know: `solve_sn` records NO residuals

`solve_sn`'s history (site B.1 #1) sets `keff_history`, `n_outer`,
`total_inner_iterations` — and **NOT `flux_residuals`, NOT `n_inner`**. The
eigenvalue path accumulates only a COUNT
(`orpheus/sn/solver.py:1012` init, `:1674` and `:1857` accumulate
`len(_residuals)`); the per-inner residual trajectory is consumed by the
certificate (`:1668`, `:1852`) and then **discarded**.

⟹ **A derived `converged` for `solve_sn` can only be OUTER-based today** — the
same predicate the adjoint already uses. Reporting *inner* truncation from the
eigenvalue entry requires new plumbing (an inner outcome accumulated alongside
`_total_inner_iterations` in `SNSolver._solve_source_iteration` /
`_solve_krylov`). That is a real scope line, and it is invisible from the
construction site.

---

## The count that matters

| category | count | where |
|---|---|---|
| **Production sites that CONSTRUCT a convergence claim (SN)** | **5** | `orpheus/sn/solver.py` — 5 distinct expressions, 0 shared helpers |
| **Production sites that BRANCH on it** | **1** | `orpheus/sn/solution.py:382` — and it branches on `history is None`, not on the value |
| **Test sites that BRANCH on it (SN family)** | **14** | 10 preconditions/guards in `tests/sn/solve|operators|verification`, 4 type unit-tests in `tests/sn/primitives/test_solution.py` |
| **Sites that merely REPORT it** | **1** doc-prose + **4** cross-method dict entries (different type) | — |
| **Deliberate non-convergence consumers (opt-out population)** | **6** | §D.1 rows 1–6 (+ 1 sibling spy at `test_krylov_curvilinear_precond_safety.py:497`); **1 uncertain** (`test_cache.py:239`) |
| **Solvers that report convergence AT ALL** | **1 of 5** | SN only; CP / MoC / MC / diffusion return nothing |
| **Places the "did the loop exhaust?" answer is DISCARDED** | **3 primitives** | `power_iteration`, `SourceIteration.solve`, `KEigenvalue.solve` |
| **Places a policy ALREADY exists (and they disagree)** | **3** | `GreenOperator` raises; `KrylovAcceleration` + `rayleigh_quotient_iteration` warn |
| **Test churn from a `warn`-by-default** | **0** | no `filterwarnings` / `-W error` anywhere in `pyproject.toml` |

**True blast radius of the semantics change: 15 branching sites (1 production +
14 test), plus 5 construction sites and 6 opt-outs to plumb.** Small. The
expensive part is not the blast radius — it is the *structural* work in
`power_iteration` and `SourceIteration`.

---

## Scope: SN-local fix, or shared cross-solver primitive?

**Shared primitive — and the evidence is not "CP/MoC/diffusion each rolled their
own" (they did not); it is stronger than that.**

1. **The information is destroyed in a SHARED primitive, not in SN.**
   `power_iteration` (`orpheus/numerics/eigenvalue.py:205`) knows whether it
   broke or exhausted and returns a 3-tuple. SN, CP, MoC, diffusion and
   `KEigenvalue` all consume that tuple. **#342 is what happens when four
   callers must each re-derive a fact the callee already had.** One of them got
   it right (`solve_sn_adjoint`), one got it wrong (`solve_sn`), and three did
   not try. Fixing only SN leaves the generator of the defect in place.
2. **`SourceIteration.solve` (`orpheus/numerics/iteration.py:632-748`) has the
   identical shape** one level down — `return psi, residual_history`, silent on
   `max_iter` — and its Krylov sibling in the SAME class family already warns
   (`:1027`). That asymmetry inside one module is the `coding-elegance`
   Pattern-2 finding proper.
3. **A convergence OUTCOME is a numerics concept, not an SN concept.** Nothing
   about "converged / inner-budget-exhausted / outer-budget-exhausted + the
   residual + the cap" mentions transport. Per the "extract a
   model-independent primitive before minting" discipline: the operation
   recurs in ≥2 non-isomorphic realizations (outer power iteration, inner
   fixed-point, GMRES `info`, Green splitting, root-find), and a non-identity
   morphism IS applied to it — **composition**: an outer solve's outcome is a
   function of its own budget AND the worst of its inners. That composition is
   exactly what SN cannot express today and what the eigenvalue path's missing
   `flux_residuals` makes impossible. It passes the type-minting criterion in
   `coding-standards`.
4. **The exception/warning vocabulary is already in `orpheus/numerics/`, not in
   `orpheus/sn/`** — `ConvergenceFailure` (`green_operator.py:127`) plus two
   `RuntimeWarning` sites. Only `ConvergenceCertificateError` is SN-local, and
   that one should stay SN-local (it is about the operator algebra's lag-death,
   not about budgets).

### Suggested shape (design input, not a decision)

- **`orpheus/numerics/convergence.py`** (new): a frozen
  `IterationOutcome` carrying `status` (an enum:
  `CONVERGED` / `INNER_BUDGET_EXHAUSTED` / `OUTER_BUDGET_EXHAUSTED` /
  `NOT_ATTEMPTED`), `residual: float | None`, `tol: float`, `iterations: int`,
  `cap: int` — with `converged` a **derived property**, so `converged=True` is
  unspellable. Plus `ConvergenceFailure` (re-homed from `green_operator`) and a
  new `ConvergenceWarning(RuntimeWarning)` the two existing warn sites adopt.
- **`power_iteration` returns it** (4th tuple element or a small result object);
  `SourceIteration.solve` and `KrylovAcceleration.solve` likewise.
- **`IterationHistory.converged` becomes a property over an
  `outcome: IterationOutcome`** — the five transcriptions collapse to zero.
- **The policy decision point is the four SN public entries** —
  `on_unconverged: Literal["raise","warn","ignore"] = "warn"` (warn costs 0 test
  churn today, §D.5; raise costs the 6 opt-outs an explicit kwarg), following
  the `KrylovAcceleration` precedent verbatim including its carve-out pattern.
- **CP / MoC / diffusion**: adding the outcome to their `*Result` types is a
  pure addition (no consumer branches on a field they don't have), so it can
  land in the same pass at near-zero risk — the "capability extension lands as
  a no-op through the single generic body" shape from `coding-standards`.

⚠ **Two `plan-authoring` §8 hazards in this design.** (a) `SolutionBase.converged()`
BRANCHES on `history is None` — populating `history` where it is currently `None`
changes behaviour at that site. (b) The `iteration <= 2 → False` rule in
`SNSolver.converged`/`KEigenvalue.converged` means a derived outer flag is
**structurally False** for `max_outer <= 2`; the `test_krylov_curvilinear_precond_safety.py:449`
site is exactly that, and a naive `len(keff_history) < max_outer` would newly
call it unconverged.

---

## Second latent lies of the #342 family

Four found, ranked by how load-bearing they are.

| # | site | the lie | severity |
|---|---|---|---|
| **1** | `orpheus/sn/solution.py:163` — `converged: bool = True` | ⛔ **the type itself defaults to a claim.** `IterationHistory()` reports converged. Safe-by-omission is the unsafe posture; #342 is this default made explicit at one call site | **highest** — it is the *mechanism* that makes #342 easy to write. Pinned by `tests/sn/primitives/test_solution.py:95` |
| **2** | `orpheus/sn/solution.py:382` — `return self.history.converged if self.history else True` | ⛔ **a history-less solution reports converged.** "I have no diagnostics" and "I converged" are the same answer at the public accessor | **high** — this is a PRODUCTION branch, the only one. Pinned by `tests/sn/primitives/test_solution.py:276` |
| **3** | `orpheus/derivations/common/solution_types.py:195` — `converged: bool = True` | same default-True, on the cross-method protocol's `CriticalSolution`; consumed by `tests/cross_method/adapters.py` into the comparison dict | **medium** — a reference-solver claim feeding a cross-method comparison |
| **4** | `orpheus/derivations/continuous/singular_eigenfunction/spectrum.py:819` and `:872` | `converged=True` hardcoded on a `CriticalSolution` wrapping `solve_case_method_{slab,sphere}_critical` — a **root-find**. The honest sibling in the SAME class at `:924` writes `converged=bool(res.converged)` | **medium** — a structurally-always-true claim on an iterative method, with the honest form 50 lines away |

Two more, weaker (defensible but **unstated** — a direct/dense method genuinely
has no iteration to fail, but nothing says so):
`orpheus/derivations/continuous/galerkin_spectral/basis_space.py:756` and `:793`;
`orpheus/derivations/continuous/fn_method/moment_space.py:347` and `:424`
(`:315` in the same file DOES state its reason — "bisection always converges
given bracket" — so `:347`/`:424` are just missing the comment).

**Not a lie, but the same family and worth naming:** `MoCResult`, `CPResult`,
`DiffusionResult` don't lie *because they say nothing*. Under a shared-outcome
design their silence becomes the thing to fix, and it is a strictly larger
population of un-notified callers than SN's one wrong flag.

---

## Close-out re-verification (L-012 discipline)

Re-ran every claim whose EMPTINESS is a finding, at close:

| claim | re-verified |
|---|---|
| "nothing in `orpheus/` outside `solution.py:382` reads `.converged`" | ✅ confirmed — the only other hit is `orpheus/numerics/eigenvalue.py:271`, which calls the 5-arg PREDICATE, and `orpheus/sn/solution.py:23` (docstring) |
| "CP / MoC / MC / diffusion / homogeneous carry no `converged` FIELD" | ✅ confirmed — zero hits |
| "no `filterwarnings` / `-W error` in `pyproject.toml`" | ✅ confirmed — 0 occurrences |

Tree state at close: `M orpheus/sn/solver.py`,
`M tests/sn/verification/analytical/test_si_convergence_rate.py`,
`M docs/theory/methods/sn/cartesian_multid.rst`,
`M docs/theory/methods/sn/slab_one_group.rst` — all the **#341
boundary-G-S rate write-up**, docs/docstring only, no convergence-contract
change. The `converged=` expressions are byte-identical to the audit's readings.
