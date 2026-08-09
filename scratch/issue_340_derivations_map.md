# Issue #340 — the `orpheus/derivations/` convergence-honesty map

Reconnaissance only. **No code changed.** Every claim carries `file:line`.
Anything not read directly off the tree is tagged **INFERRED**.

## §0 — provenance

- Tree: `main` @ `4bcce0bd`, `orpheus/derivations/` clean (no uncommitted edits).
- Interpreter: `.venv/bin/python` (host).
- Measurements tagged `[M]` were produced by throwaway probes under
  `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/` (not committed); the probe source is
  quoted inline where the number is load-bearing.
- **Nexus workspace check (L22 hazard):** the active graph answered with
  `/Users/rodrigo/git/nuclear/ORPHEUS/...` paths and `spectrum.py:814/867/919`,
  matching the main checkout. The sibling worktree at
  `.claude/worktrees/nexus-workspace-wiring/` carries an *older* copy
  (`spectrum.py:797/850/902`) and is **excluded from every count below**.
- **Scope.** Q3 covers `orpheus/derivations/` as asked. The *top-level*
  `derivations/` directory is a different tree (diagnostics); its 5 `.converged`
  reads (`diagnostics/diag_d3_absorber_01_unconverged_exit.py:115`,
  `..._02_si_rate_scaling.py:112, 148`) are on the SN `IterationHistory` and are
  already honest post-`d9b027d7`.

## §0.1 — headline

**The issue's three sites are an undercount, but not in the direction it looks.**
Two corrections, both measured:

1. **The `CriticalSolution` layer has FOUR live defects, not two** — the
   `fn_method/moment_space.py` pillar carries the same hardcode twice, and its
   sphere arm is the *worst* site in the subtree because it reads a scipy
   `.success` flag, **branches on it**, and then discards it.
2. **Removing the `:195` default costs nothing** — `[M]` all 9 construction sites
   already pass `converged=` explicitly (0 breakages), and there are only **2**
   dynamic reads of `CriticalSolution.converged` in the whole exercised tree. The
   issue's "~87 reads" is `grep '\.converged' tests/derivations/`, of which
   **72 belong to the unrelated Peierls Green's-function family** (see Q2.2).

| # | site | shape | producer's method | verdict |
|---|---|---|---|---|
| — | `common/solution_types.py:195` | `converged: bool = True` DEFAULT | n/a | ⛔ dead weight — **free to remove** |
| 1 | `singular_eigenfunction/spectrum.py:819` | hardcoded `True` | hand bisection + budget | ⛔ **DEFECT** (issue site) |
| 2 | `singular_eigenfunction/spectrum.py:872` | hardcoded `True` | hand bisection + budget | ⛔ **DEFECT** (issue site's sibling) |
| 3 | `fn_method/moment_space.py:315` | hardcoded `True` **+ a false justifying comment** | hand bisection + budget | ⛔ **DEFECT — not in the issue** |
| 4 | `fn_method/moment_space.py:347` | hardcoded `True` | `minimize_scalar` + silent fallback | ⛔ **DEFECT — not in the issue; the worst one** |
| — | `fn_method/moment_space.py:424` | hardcoded `True` | `compute_kinf_*`, direct `np.linalg.solve` | ✅ legitimate — nothing iterates |
| — | `galerkin_spectral/basis_space.py:756` | hardcoded `True` | direct `scipy.linalg.eig` | ✅ legitimate — nothing iterates |
| — | `galerkin_spectral/basis_space.py:793` | hardcoded `True` | direct `scipy.linalg.eig` | ✅ legitimate — nothing iterates |
| — | `singular_eigenfunction/spectrum.py:924` | `bool(res.converged)` | `brentq(..., full_output=True)` | ✅ **the exemplar** |
| — | `trajectory_resolvent/billiard.py:755` | `bool(res.converged)` | power iteration | ✅ honest |

⭐ **The idiom the carve needs already exists in the tree — and it is subtler than
it looks.** `singular_eigenfunction/cylinder/one_group.py:795-806` calls

```python
    R_c, brent_result = scipy.optimize.brentq(
        residual, sign_change_bracket[0], sign_change_bracket[1],
        xtol=bisect_tol, rtol=bisect_tol, maxiter=max_iter,
        full_output=True,
        disp=False,
    )
```

and records **both** `iterations=int(brent_result.iterations)` (`:817`) and
`converged=bool(brent_result.converged)` (`:818`).

The load-bearing detail is **`disp=False`**. `[M]` with scipy's default `disp=True`,
a non-converged `brentq` **raises `RuntimeError` instead of returning
`converged=False`**, so the flag would be structurally unreachable even with
`full_output=True`. Measured on `f(x) = eˣ − 3x − 1` over `[1, 3]`:

| `maxiter` | `res.converged` | `res.flag` |
|---|---|---|
| 1 / 2 / 3 | **`False`** | `convergence error` |
| 100 | `True` | `converged` |

`max_iter` is caller-tunable on the cylinder solver (`:660`, default 200), so the
`False` leg is reachable through the public API. **This author got all three pieces
right** (`full_output=True` + `disp=False` + record both fields) and is the pattern
the other four producers should be rewritten to. Nothing has to be invented.

> Contrast, and it is a policy question the carve should answer explicitly: the
> tree's *other* `brentq` call sites — `singular_eigenfunction/core/dispersion.py:207`,
> `cases/sn.py:650`, `cases/diffusion.py:597` — omit `full_output`, so they inherit
> `disp=True` and are **honest-by-raising**. Both policies are honest; only the
> cylinder's is *readable*. Since #340's ratified ruling for the SN half was
> **"warn, CI-escalatable — not raise"**, the readable policy is the one consistent
> with what already landed in `d9b027d7`.

---

## Q1 — does the producer already KNOW whether it converged?

### Q1.A — `CaseMethodSlabResult` (issue site 3): the loop BREAKS ON A TEST

`orpheus/derivations/continuous/singular_eigenfunction/slab/one_group.py:332-345`,
verbatim:

```python
    iters = 0
    while iters < max_bisect:
        if d_hi - d_lo < bisect_tol * max(1.0, d_lo):
            break
        d_mid = 0.5 * (d_lo + d_hi)
        f_mid = _eq46_residual(
            d_mid, c=c, R=R, f1=f1, u0=u0, nu_bar=nu_bar, z0=z0,
            dps=dps, maxdegree=maxdegree, branch_sign=branch_sign,
        )
        if f_mid * f_lo < 0:
            d_hi, f_hi = d_mid, f_mid
        else:
            d_lo, f_lo = d_mid, f_mid
        iters += 1
```

**Verdict: the fact EXISTS and is discarded, exactly like `power_iteration`'s was.**
The loop has two distinguishable exits —

- **certified**: `break` at `:335`, i.e. the bracket predicate
  `d_hi - d_lo < bisect_tol * max(1.0, d_lo)` held;
- **budget**: the `while` guard at `:333` fails, i.e. `iters == max_bisect`.

Only `iters` survives into the result (`:366`, `n_bisect_iters=iters`). The
*reason* for the exit is thrown away at the return statement.

**⚠ A TOLERANCE ALREADY EXISTS — nothing has to be invented.** `bisect_tol` is a
declared parameter with a documented default:

- `slab/one_group.py:211` — `bisect_tol: float = 1e-10,`
- `slab/one_group.py:238-239` — *"bisect_tol : float, default 1e-10 / Tolerance
  on :math:`d`."*
- forwarded from the facade at `spectrum.py:677` (`bisect_tol: float = 1e-10`) →
  `:806` → into the solver.

So the honest predicate is **already fully spelled in the tree** and needs no new
number. The strongest form is *not* `iters < max_bisect` but the bracket
predicate re-evaluated after the loop, because `d_lo`/`d_hi` are still in scope at
`:347`:

```python
converged = (d_hi - d_lo) < bisect_tol * max(1.0, d_lo)
```

> ⚠ **Do not use `iters < max_bisect` as the predicate.** `iters` is incremented at
> the *bottom* of the body (`:345`) and the tolerance test sits at the *top*
> (`:334`), so the final bisection step's narrowing is never tested: a run that
> exits on `iters == max_bisect` may in fact have landed inside tolerance. The
> bracket-width re-evaluation is the only spelling that is correct on both exits.
> The `eq46_residual` is **not** needed for this — it is a *diagnostic of the
> equation*, whereas `bisect_tol` is a tolerance *on `d`*, and only the latter is
> stated anywhere.

`[M]` starved-budget behaviour, measured directly:

```
[SE slab, default max_bisect=80] d=5.665505455443451 n_bisect_iters=29 eq46_residual=-5.068e-11
[SE slab, max_bisect=3        ] d=5.679726758793970 n_bisect_iters=3  eq46_residual= 3.511e-03
    -> d error vs converged run: 1.422e-02
```

and end-to-end through the facade — the claim reaches the consumer as `True`:

```
SE slab Spectrum.solve_critical
  default : param=5.665505455443451 converged=True iters=29
  mb=3    : param=5.679726758793970 converged=True iters=3  eq46_res=3.511e-03
```

A 1.4e-2 mfp error in a critical half-thickness, reported as converged.

### Q1.B — `CaseMethodSphereResult` (issue site 2's other half): identical

`orpheus/derivations/continuous/singular_eigenfunction/sphere/one_group.py:232-241`
is the same loop with `r_lo/r_hi`:

- `:232` `while iters < max_bisect:`
- `:233` `if r_hi - r_lo < bisect_tol * max(1.0, r_lo):`
- `:234` `break`
- `:241` `iters += 1`
- `:258` `n_bisect_iters=iters` — same discard.
- `:144` `bisect_tol: float = 1e-10` — same already-stated tolerance.

`CaseMethodSphereResult` (`:41-79`) likewise has **no `converged` field**; it
carries `eq54_residual` (`:77`) and `n_bisect_iters` (`:79`). Q1's answer is
identical: the fact exists, the tolerance exists, only the statement is missing.

### Q1.C — the scipy question: **YES, a `.success` flag is read and then dropped**

This is the sharpest instance in the whole subtree, and it is **not** in the
issue's list. `orpheus/derivations/continuous/fn_method/sphere/one_group.py:352-364`:

```python
        result = minimize_scalar(
            log_abs_det_at_R,
            bracket=(R_lo_bracket, R_guess, R_hi_bracket),
            method="brent",
            options={"xtol": bisect_tol, "maxiter": max_bisect},
        )
    if not result.success:
        # Fall back to the coarse-grid guess.
        R_crit = R_guess
    else:
        R_crit = float(result.x)
```

The flag is not merely dropped — it is **branched on**, and then not recorded.
`SphereFNResult` (`:373-382`) has no `converged` field, and
`fn_method/moment_space.py:346` writes a hardcoded `converged=True` on top of it.

`[M]` by spying on `minimize_scalar` at the call site:

| `max_bisect` | scipy `.success` | scipy `nit` | scipy `x` | value RETURNED | fell back? |
|---|---|---|---|---|---|
| 80 | `True` | 32 | 2.424824910627988 | 2.424824910627988 | no |
| 5 | **`False`** | 5 | 2.424740494069838 | **2.4317897371714645** | **yes** |
| 2 | **`False`** | 2 | 2.4200773879849815 | **2.4317897371714645** | **yes** |
| 1 | **`False`** | 1 | 2.4200773879849815 | **2.4317897371714645** | **yes** |

scipy's own message on the failing rows is `Maximum number of iterations exceeded`.
End-to-end through the facade:

```
FN sphere MomentSpace.solve_critical  (minimize_scalar path)
  default : param=2.4248249106283257 converged=True
  mb=1    : param=2.4317897371714645 converged=True     <-- 7.0e-3 error, "converged"
```

> **Secondary finding, worth its own line in the carve.** The fallback at `:359`
> is *worse than the thing it replaces*. At `max_bisect=5` scipy's unconverged
> iterate was accurate to **8.4e-5**; the coarse-grid guess it fell back to is off
> by **7.0e-3** — an **83×** degradation. The fallback discards a good answer for
> a bad one AND relabels it converged. (Independent of #340's contract question,
> but it is on the same three lines.)

### Q1.D — the other root-finders on these paths (status: honest, or absent)

| call site | finder | status handling |
|---|---|---|
| `singular_eigenfunction/core/dispersion.py:207` | `brentq(_Lambda_real, ..., xtol=bisect_tol)` | `full_output` NOT set ⟹ `disp=True` default ⟹ scipy **raises `RuntimeError`** on non-convergence. Honest-by-raising; the fact cannot be silently lost here. It also stores `Lambda_residual` (`:213`) as a witness. |
| `cases/sn.py:650`, `cases/diffusion.py:597` | `brentq(det_c, ...)` | same — no `full_output`, so `disp=True` raises. Honest-by-raising. |
| `fn_method/slab/one_group.py:227-235` | hand-rolled `for _ in range(max_bisect)` + `break` on `a_hi - a_lo < bisect_tol * max(1.0, a_lo)` | **fact exists, fully discarded — and there is not even an iteration counter.** `SlabFNResult` has no `n_bisect_iters`. Worse than the SE family, where a consumer could at least infer from `iters == max_bisect`. |
| `fn_method/sphere/one_group.py:353` | `minimize_scalar` | Q1.C — read, branched on, discarded. |

`[M]` FN-slab starvation, and the field list proving the counter is absent:

```
[FN slab default] a=5.665279370885399 det_res=-4.802e-15+8.043e-16j
[FN slab mb=2   ] a=5.668749999999999 det_res=-1.335e-05+2.308e-06j   (err 3.47e-03)
    SlabFNResult fields: ['N', 'a_critical_mfp', 'c', 'coefficients_a',
                          'determinant_residual', 'n_bracket_points', 'nu0',
                          'xi_collocation']
```

### Q1.E — ⛔ the comment at `moment_space.py:315` is present-tense FALSE

```python
            converged=True,  # bisection always converges given bracket
```

This is a *justification*, which makes it worse than a bare hardcode: it tells the
next reader the audit is already done. It is false in the sense that matters — the
loop's budget `max_bisect` is a **caller-tunable parameter** that
`MomentSpace.solve_critical` forwards (`moment_space.py:227-232`, `:306`, `:338`),
so "given bracket" does not imply "given budget". Measured:

```
FN slab MomentSpace.solve_critical
  default : param=0.9377197763260483 converged=True
  mb=2    : param=0.94375            converged=True    <-- 6.0e-3 error, "converged"
```

This is structurally the same shape as #340's own `max_inner=1000`: a
dimension-blind caller-tunable budget whose exhaustion is unobservable at the
entry point.

**⚠ One caveat, stated so the carve does not over-claim.** With the *shipped
defaults* (`bisect_tol=1e-12`/`1e-10`, `max_bisect=80`, bracket ≈ `(d_max −
d_min)/(n_bracket−1)`) all four hand-rolled bisections do converge — 80 halvings
of an O(0.3) bracket reaches ~2.5e-25, far inside tolerance, and the observed SE
slab run needed **29**. So this is a **latent** defect, not a currently-firing one:
it fires only for a caller who passes a small `max_bisect`, or for a bracket much
wider than the defaults. That is exactly the state #340's SN half was in before
#337 moved the budget across the line.

### Q1 — summary answer

**The producer knows, at every one of the five sites.** No tolerance needs to be
invented anywhere:

| producer | the fact, in the tree today | already-stated tolerance |
|---|---|---|
| `solve_case_method_slab_critical` | `d_hi - d_lo` at `slab/one_group.py:347` | `bisect_tol` `:211` (doc `:238`) |
| `solve_case_method_sphere_critical` | `r_hi - r_lo` at `sphere/one_group.py:~243` | `bisect_tol` `:144` |
| `solve_fn_slab_bare_critical` | `a_hi - a_lo` at `fn_method/slab/one_group.py:237` | `bisect_tol` `:128` |
| `solve_fn_sphere_bare_critical` | `result.success` at `fn_method/sphere/one_group.py:358` | scipy's own `xtol` |

The recommended shape in the issue comment ("give `CaseMethodSlabResult` the fact
its own root-find already has") is **correct and applies to all four producers**,
and the fourth one needs nothing but `converged=bool(result.success)`.

---

## Q2 — blast radius of removing the `converged: bool = True` default

### Q2.0 — the headline: **removing the default breaks ZERO sites**

Field under audit: `orpheus/derivations/common/solution_types.py:195`
(`converged: bool = True`), documented at `:164-165`.

`[M]` **static**: every one of the **9** `CriticalSolution(...)` construction
sites in `orpheus/` passes `converged=` explicitly.

`[M]` **dynamic**, by an in-process pytest plugin that wraps
`CriticalSolution.__init__` and `__getattribute__` (throwaway, at
`/Users/rodrigo/.claude/jobs/c30e4f25/tmp/plugin_criticalsolution_audit.py`; no
tracked file edited), run over every consumer suite —
`tests/cross_method tests/derivations/test_fn_method_moment_space.py
tests/derivations/test_singular_eigenfunction_spectrum.py
tests/derivations/test_galerkin_spectral_basis_space.py
tests/derivations/test_trajectory_resolvent_billiard.py` under `python -O`,
**140 passed in 250 s**:

```
constructions WITH converged=   : 33
constructions WITHOUT converged=: 0   <-- these break if the default goes

.converged READS (dynamic, on CriticalSolution only): 2
    x1  tests/derivations/test_singular_eigenfunction_spectrum.py:287  in test_critical_solution_is_correct_type
    x1  tests/derivations/test_trajectory_resolvent_billiard.py:304    in test_billiard_sphere_1g_bit_equal_legacy
```

*(Instrument's positive control: both counters are non-zero, i.e. the wrap was
live and observed real traffic. A 0/0 report would have meant a dead plugin.)*

⟹ **`converged: bool = True` at `:195` is pure dead weight.** It protects nothing;
its only effect is to make the omission spellable, which is exactly the hazard.
Making it required is a **zero-churn** edit at the constructors — strictly cheaper
than the `IterationHistory` half of the first pass. It is the *hardcoded values*
being passed, not the default, that carry the whole defect.

### Q2.1 — the 9 construction sites, with what each passes

| `file:line` (the `converged=` arg) | geometry / kind | value | honest? |
|---|---|---|---|
| `galerkin_spectral/basis_space.py:756` | slab | `True` | ⛔ hardcoded |
| `galerkin_spectral/basis_space.py:793` | sphere | `True` | ⛔ hardcoded |
| `fn_method/moment_space.py:315` | slab | `True  # bisection always converges given bracket` | ⛔ hardcoded **+ false comment** |
| `fn_method/moment_space.py:347` | sphere | `True` | ⛔ hardcoded |
| `fn_method/moment_space.py:424` | `solve_kinf` | `True` | ✅ *probably legitimate* — see Q3 note |
| `singular_eigenfunction/spectrum.py:819` | slab | `True` | ⛔ hardcoded (issue site) |
| `singular_eigenfunction/spectrum.py:872` | sphere | `True` | ⛔ hardcoded (issue site) |
| `singular_eigenfunction/spectrum.py:924` | cylinder | `bool(res.converged)` | ✅ honest |
| `trajectory_resolvent/billiard.py:755` | all (`_build_critical`) | `bool(res.converged)` | ✅ honest |

**The issue's list is missing `galerkin_spectral/basis_space.py:756` and `:793`** —
a third pillar with the identical hardcode, and `basis_space.py` contains **no other
`converged` token at all**, i.e. its producers were never asked the question.

The construction sites are reachable only through 5 methods
(`Spectrum.solve_critical`, `MomentSpace.solve_critical`, `MomentSpace.solve_kinf`,
`BasisSpace.solve_critical`, `Billiard.solve_critical`). `[M]` no dynamic
construction path exists: the only `dataclasses.replace` calls in the family target
`PeierlsSphereSolution`/`PeierlsCylinderSolution.phi_values`
(`peierls_nystrom/sphere.py:141, 286`, `cylinder.py:127, 290`) and `CrossMethodCase`
(`tests/cross_method/test_eigenvalue.py:173`) — none touch `CriticalSolution`; and
there is no `CriticalSolution(**…)` splat and no `asdict`/`astuple` round-trip
anywhere in `orpheus/derivations/`, `tests/derivations/`, or `tests/cross_method/`.
This is what makes the dynamic audit's `0` a complete answer rather than a sample.

### Q2.2 — what reads `.converged`, and what it does with it

**Tool disagreement, and which found what.** Nexus and grep answer *different*
questions here and neither alone is sufficient — as expected, since a call graph
cannot see attribute access:

- **Nexus** (`context` on `py:class:...CriticalSolution`) resolved the class's
  **producers** perfectly: 9 `type_uses`/`calls` edges naming every
  `_solve_critical_*` / `_wrap_*` / `_build_critical`, plus the 4 doc pages that
  render it. It also exposed `py:attribute:...CriticalSolution.converged` with
  **`degree: 1`** — the attribute node's ONLY edge is `contains` from the class.
  **Nexus sees zero readers**, correctly, because there are no `calls`/`references`
  edges for a plain dataclass attribute read.
- **Text-grep** found 153 `.converged` occurrences tree-wide, of which **87** are in
  `tests/derivations/` — and 87 is exactly the number in the issue. But that number
  is a **mis-attribution**: decomposed by family,

  | family | `.converged` hits in `tests/derivations/` |
  |---|---|
  | `test_peierls_*` | **72** |
  | `test_fn_*` | 8 |
  | `test_trajectory_*` | 4 |
  | `test_singular_*` | 3 |
  | `test_galerkin_*` | **0** |

  The 72 Peierls hits are on `PeierlsGreensFunction*Result.converged` — a **different
  type family** that is already honest (it carries the fact and the tests assert it).
  They are untouched by this field.
- **The dynamic audit** settles it: **2** reads of `CriticalSolution.converged`
  exist in the exercised tree, both in tests, both benign.

**What the 2 readers do:**

| `file:line` | statement | branch / report / ignore | breaks on `False`? |
|---|---|---|---|
| `tests/derivations/test_singular_eigenfunction_spectrum.py:287` | `assert sol.converged is True` | **branch (assert)** | **YES** — but only if the SE slab genuinely fails to converge, which is the point |
| `tests/derivations/test_trajectory_resolvent_billiard.py:304` | `assert sol.converged == legacy.converged` | **branch (assert)**, relative | no — it compares the two sides, both would move together |

Every other consumer **ignores** it.

### Q2.3 — the cross-method comparison dict: **`CriticalSolution` is NOT in it**

Direct answer to the question asked. The cross-method protocol
(`tests/cross_method/`) does **not** consume `CriticalSolution`:

- The comparison carrier is its own `ScalarResult` (`tests/cross_method/protocol.py`)
  with a free-form `metadata` dict documented at `protocol.py:92` as carrying
  `iterations`, `converged`, ….
- Every adapter in `tests/cross_method/adapters.py` calls the **function-level**
  solver and bypasses `CriticalSolution` entirely — e.g. `FNSlabAdapter.solve`
  (`:55-72`) calls `solve_fn_slab_bare_critical` directly, `FNSphereAdapter.solve`
  (`:88`) calls `solve_fn_sphere_bare_critical`.
- Only 4 adapters populate `metadata["converged"]` — `adapters.py:196` (reflected
  slab), `:270`, `:340`, `:404` (trajectory-resolvent) — and all 4 read the
  *function-level* result's own honest flag, not `CriticalSolution`'s.
- Those 4 ARE gated: `tests/cross_method/test_eigenvalue.py:345, 384, 422, 463`
  each `assert res.metadata["converged"]`. **A `False` there would RED the gate**
  (it is a bare assert inside a collected test module, so the rewriter protects it
  under `-O` — Mode 8 does not apply).

**Two consequences the carve should carry:**

1. ⛔ **The two adapters most exposed to Q1's defect are the two with no
   `converged` key at all.** `FNSlabAdapter` (`:55`) and `FNSphereAdapter` (`:88`)
   drive exactly the `solve_fn_slab_bare_critical` / `solve_fn_sphere_bare_critical`
   producers that Q1.C/Q1.D showed to be the *least* honest — the sphere one being
   the `minimize_scalar` fallback. Their `metadata` dicts (`:65-70`, `:96-100`)
   record `determinant_residual` but no convergence flag. So the cross-method gate
   is structurally blind to the exact producers whose exit is unreported.
2. ⛔ **`tests/cross_method/test_polymorphism.py` consumes `CriticalSolution` and
   never reads `.converged`.** It calls `MomentSpace.solve_critical()` (`:92`, `:124`)
   and `Billiard.solve_critical()` (`:159`, `:194`, `:228`) and asserts only
   `parameter_kind` / `eigenvalue_kind` / a value drift (`:94-98`, `:126-130`,
   `:161-165`, `:196-200`, `:230-233`). This is #340's *"audit every gate for a
   missing `converged` assertion"* item, in the `derivations/` family: **5 rows.**

### Q2.4 — docs blast radius

`[M]` `grep -rn "\.converged\|converged=" docs/` (excluding `_build`) returns **8**
hits and **none** name `CriticalSolution`: 2 are MoC/CP `Solver.converged()` methods,
6 are the SN `IterationHistory` / `PowerIterationOutcome` contract prose landed by
`d9b027d7` (`docs/theory/methods/sn/solver.rst:603-626`,
`docs/theory/conventions/indexing_and_layout.rst:922, 996`).

⟹ **The four theory pages that Nexus shows rendering `CriticalSolution`**
(`theory/references/{trajectory_resolvent, fn_method, singular_eigenfunction,
galerkin_spectral}`) **document the type but say nothing about its convergence
semantics.** The carve therefore has a *doc gap to fill*, not doc breakage to
repair — and per the retirement rule the `converged` attribute prose at
`solution_types.py:164-165` ("Whether the iteration / root-find converged to its
tolerance") is currently a claim **no producer honours at 7 of 9 sites**.

## Q3 — sibling scan across the whole subtree

`[M]` measured counts:

```
grep -rn --include='*.py' converged orpheus/derivations/   -> 138 hits, 28 files
grep -rn            converged orpheus/derivations/         -> 139 hits  (+1: spectral_collocation/README.md)
  of which CONTRACT hits (`converged=` / `converged:` / `converged =` / `.converged`)  ->  99
  of which PROSE-only hits (the English word in a docstring/comment)                   ->  39
```

Classification below is exhaustive over the **99 contract hits**; the **39 prose**
hits are enumerated by file in Q3.3 (they carry no contract — with one exception,
flagged there).

### Q3.1 — the defect rows (⛔ hardcoded / dropped / dead)

| `file:line` | hit | class | note |
|---|---|---|---|
| `common/solution_types.py:195` | `converged: bool = True` | **default True** | zero-cost to remove (Q2.0) |
| `singular_eigenfunction/spectrum.py:819` | `converged=True,` | **hardcoded True** | producer knows (Q1.A) |
| `singular_eigenfunction/spectrum.py:872` | `converged=True,` | **hardcoded True** | producer knows (Q1.B) |
| `fn_method/moment_space.py:315` | `converged=True,  # bisection always converges given bracket` | **hardcoded True + false comment** | Q1.E |
| `fn_method/moment_space.py:347` | `converged=True,` | **hardcoded True** | producer reads `.success` then drops it (Q1.C) |
| `peierls_nystrom/slab.py:596` | `converged = abs(k_new - k_val) < tol and iteration > 5` | **field absent, fact computed then dropped** | branched on at `:600` (`if converged: break`); `PeierlsSlabSolution` (`:385-402`) has **no `converged` field** |
| `peierls_nystrom/geometry.py:6450` | `converged = abs(k_new - k_val) < tol and it > 5` | **field absent, fact computed then dropped** | branched on at `:6453`; `PeierlsSolution` (`:5442-5458`) has **no `converged` field** |
| `fn_method/slab/flux_reconstruction.py:854, 866` | `converged = False` / `converged = True` | **DEAD LOCAL** | `[M]` no `converged` token exists anywhere after `:866` in the file — the flag is assigned inside `slab_scalar_flux_fn_projection`'s power iteration and **never read by anything**, in the same function that computed it |
| `singular_eigenfunction/slab/one_group.py` (type at `:81-121`) | — | **field absent but derivable** | `d_hi - d_lo` vs the stated `bisect_tol` (Q1.A) |
| `singular_eigenfunction/sphere/one_group.py` (type at `:41-79`) | — | **field absent but derivable** | `r_hi - r_lo` vs `bisect_tol` (Q1.B) |
| `fn_method/slab/one_group.py` (`SlabFNResult`) | — | **field absent, NOT even an iteration counter** | `[M]` fields are `['N','a_critical_mfp','c','coefficients_a','determinant_residual','n_bracket_points','nu0','xi_collocation']` |
| `fn_method/sphere/one_group.py` (`SphereFNResult`) | — | **field absent; scipy `.success` read at `:358` and discarded** | Q1.C — the worst site |

> ⭐ **Two new failure shapes the issue does not name.**
> (a) **Fact computed, branched on, then dropped at the return** — the two Peierls
> Nyström power iterations. This is *literally* the `power_iteration` defect
> `d9b027d7` fixed, in a second family: the flag is a named local, it controls the
> `break`, and the result dataclass has no slot for it.
> (b) **The dead local** — `fn_method/slab/flux_reconstruction.py:854/866`. The
> function *names* the fact and drops it in the same scope. Nothing downstream can
> ever recover it.
>
> ⚠ Also note `peierls_nystrom/slab.py:596` / `geometry.py:6450` are a **sixth
> spelling** of the convergence predicate — `abs(Δk) < tol **and** iteration > 5`,
> a hardcoded minimum-iteration floor folded into the test. `d9b027d7` collapsed
> five transcriptions into one `_claims_convergence`; these two live outside
> `orpheus/numerics` + `orpheus/sn` and were not swept.

### Q3.2 — the honest rows (✅ read, recorded, or legitimately constant)

| `file:line` | hit | class |
|---|---|---|
| `singular_eigenfunction/cylinder/one_group.py:795, 817-818` | `brentq(..., full_output=True)` → `iterations=`, `converged=bool(brent_result.converged)` | ✅ **honest read — the exemplar** |
| `singular_eigenfunction/spectrum.py:924` | `converged=bool(res.converged)` | ✅ honest read (of the row above) |
| `trajectory_resolvent/power_iteration.py:107, 189, 214, 221` | field + `False` init + `True` on break + recorded | ✅ honest |
| `trajectory_resolvent/greens_function.py:375/386, 593/606, 1069/1083, 1188/1236/1246` | `converged = pi_result.converged` → recorded | ✅ honest (4 arms) |
| `trajectory_resolvent/greens_function_cylinder.py:423/435, 582/596, 945/960` | same | ✅ honest (3 arms) |
| `trajectory_resolvent/greens_function_annulus.py:409/421, 570/584` | same | ✅ honest |
| `trajectory_resolvent/greens_function_hollow_sphere.py:330/341, 487/498` | same | ✅ honest |
| `trajectory_resolvent/greens_function_slab_asymmetric.py:332/343, 484/495` | same | ✅ honest |
| `trajectory_resolvent/greens_function_slab.py:261, 342` | `converged=res_asym.converged` | ✅ honest (delegates to the asymmetric solver) |
| `trajectory_resolvent/billiard.py:698, 755` | `bool(res.converged)` | ✅ honest |
| `fn_method/slab/reflected.py:119, 875, 950, 956, 973` | field + `False` init + `True` on break + `if converged:` + recorded | ✅ honest |
| `fn_method/slab/flux_reconstruction.py:223, 368, 371, 385` | field + both polarities + recorded | ✅ honest |
| `fn_method/sphere/flux_reconstruction.py:85, 220, 223, 236` | field + both polarities + recorded | ✅ honest |
| `peierls_nystrom/ps1982_reference.py:99, 182, 197, 219` | field + both polarities + recorded | ✅ honest |
| `fn_method/moment_space.py:424` | `converged=True` on `solve_kinf` | ✅ legitimate — `compute_kinf_mg` is `float(nu_sigma_f @ np.linalg.solve(A, chi))` (`multi_group/k_inf.py:335-336`); nothing iterates |
| `galerkin_spectral/basis_space.py:756, 793` | `converged=True` | ✅ legitimate — `solve_eq4_eigenproblem` is a direct `scipy.linalg.eig` (`core/galerkin_matrix.py:157`); nothing iterates |

> **Design note for the carve, not a defect.** The five `converged=True` rows above
> are honest *given a boolean field*, but they conflate "converged" with "not
> applicable — this is a direct method". If the carve ever mints a richer exit type
> (the `PowerIterationOutcome` analogue), a third state is the principled spelling.
> Flagging it because the two categories currently look identical in a grep, which
> is precisely how the four real defects have been hiding among them.

### Q3.3 — unrelated (prose only, no contract)

`[M]` **39** hits are English usage of the word in docstrings/comments — "the converged
flux", "at the converged :math:`d`", "converged to quadrature accuracy" (+1 in
`spectral_collocation/README.md`, non-`.py`). Files:
`solution_types.py:142/164/165/172` (the field's own doc — see the ⚠ below),
`fn_method/slab/one_group.py:24/83`, `fn_method/sphere/one_group.py:40/135`,
`fn_method/slab/flux_reconstruction.py:10/14/202/203/512`,
`fn_method/sphere/flux_reconstruction.py:72/73/376`, `fn_method/slab/reflected.py:102/103/955`,
`fn_method/origins/*` (4), `analytical/homogeneous.py:561`,
`peierls_nystrom/ps1982_reference.py:96`, `peierls_nystrom/geometry.py:2129/2231/3272`,
`mms/sn.py` (10), `singular_eigenfunction/spectrum.py:958`,
`singular_eigenfunction/slab/one_group.py:100/104/349`,
`singular_eigenfunction/sphere/one_group.py:62`,
`singular_eigenfunction/cylinder/one_group.py:175/181/706`,
`singular_eigenfunction/origins/cylinder_derivations.py:642/707`,
`trajectory_resolvent/power_iteration.py:92/95/144`, `trajectory_resolvent/billiard.py:117`,
plus the per-family `converged: bool` field *declarations* already counted in Q3.2.

⚠ **One prose hit is a present-tense-FALSE claim, not neutral prose.**
`common/solution_types.py:164-165` documents the field as *"Whether the iteration /
root-find converged to its tolerance."* `[M]` That is false at 4 of 9 producers
(Q3.1) and vacuous at 3 more (the direct methods). Per the retirement rule it is a
MUST-FIX in the same change that lands the carve — a reader who trusts it will
consume a `True` that no root-find ever produced.

---

## Q4 — the three `spectrum.py` call paths: what actually differs

All three live in the same class, within ~130 lines, and are dispatched from one
`if/elif` (`spectrum.py:774-783`). The difference is **entirely in the producer they
delegate to** — not in anything the `Spectrum` methods themselves do or know.

| | `:819` slab | `:872` sphere | `:924` cylinder |
|---|---|---|---|
| method | `Spectrum._solve_critical_slab` (`:789`) | `Spectrum._solve_critical_sphere` (`:837`) | `Spectrum._solve_critical_cylinder` (`:890`) |
| producer | `solve_case_method_slab_critical` (`slab/one_group.py:203`) | `solve_case_method_sphere_critical` (`sphere/one_group.py:136`) | `solve_singular_eigenfunction_cylinder_bare_critical` (`cylinder/one_group.py`) |
| result type | `CaseMethodSlabResult` (`:81-121`) | `CaseMethodSphereResult` (`:41-79`) | the cylinder result (`:170-202`) |
| **has a `converged` field?** | **NO** | **NO** | **YES** (`cylinder/one_group.py:202`, doc `:175-181`) |
| root-find | hand-rolled `while iters < max_bisect` + `break` (`:333-345`) | hand-rolled `while iters < max_bisect` + `break` (`:232-241`) | `scipy.optimize.brentq(..., full_output=True)` (`:795`) |
| the fact, in the tree | `d_hi - d_lo` at `:347` (live, unrecorded) | `r_hi - r_lo` (live, unrecorded) | `brent_result.converged` (recorded at `:818`) |
| iteration count recorded | `n_bisect_iters` (`:366`) | `n_bisect_iters` (`:258`) | `iterations` (`:817`) |
| residual recorded | `eq46_residual` (`:364`) | `eq54_residual` (`:256`) | `criticality_residual` |
| ⟹ `spectrum.py` writes | `converged=True` | `converged=True` | `converged=bool(res.converged)` |

### Why did one get it right?

**Not because the cylinder author was more careful about convergence — because the
cylinder author reached for `scipy.optimize.brentq` with `full_output=True`.**
`full_output=True` returns a `RootResults` object whose `.converged` and
`.iterations` are simply *there*; recording them costs one line each, and the author
took both (`:817-818`). The slab and sphere authors hand-rolled the bisection, and a
hand-rolled loop hands you nothing — the fact has to be *constructed* at the return
statement, which is exactly the step that got skipped.

So the divergence is **not** a knowledge gap at the `spectrum.py` layer. All three
`Spectrum._solve_critical_*` methods are 30-line adapters that faithfully forward
whatever their producer gives them; the two that write `True` have literally nothing
else available on the object in hand. **This is the "triage one hop UP" ruling
verbatim** — the consumer asserts a fact it should have read, and the producer's
return type is the defect.

### Do the other two "genuinely lack the fact"?

**No.** Both have it, and both have a stated tolerance to compare against —
Q1.A/Q1.B. The lack is one of *statement*, not of *knowledge*:

- slab: `d_lo`/`d_hi` are in scope at the `return` (`:359-367`); `bisect_tol` is a
  documented parameter (`:211`, `:238`).
- sphere: `r_lo`/`r_hi` likewise; `bisect_tol` at `:144`.

Both `CaseMethodSlabResult` and `CaseMethodSphereResult` are frozen dataclasses with
no `converged` slot, so `spectrum.py` cannot read what was never written. Adding the
field at the producer makes all three `spectrum.py` sites collapse to the identical
`converged=bool(res.converged)` line, which is the Pattern-2 single-spelling
resolution the issue's "Recommended shape" asks for.

⚠ **Scope caution — the residual is NOT the right witness here.** `eq46_residual` /
`eq54_residual` are documented as *"Should be small"* with **no threshold anywhere in
the tree**, and they measure the *equation*, not the *root-find*. Deriving
`converged` from them would require inventing a tolerance — the exact move #340
exists to remove. The bracket width against the already-stated `bisect_tol` needs
nothing invented. (Separately: for the FN sphere, `determinant_residual` is a
`[M]` ~1e-22 scale-dependent quantity at BOTH the converged and the starved run, so
it is not usable as a convergence witness even informally.)

---

## §5 — what the carve inherits (ordered by evidence strength)

1. **Free win, zero churn** — drop the `= True` default at `solution_types.py:195`
   (`[M]` 0 of 33 constructions rely on it) and fix the now-false field doc at
   `:164-165`.
2. **The four real producer defects**, all with the fact and the tolerance already
   in the tree: `singular_eigenfunction/{slab,sphere}/one_group.py`,
   `fn_method/slab/one_group.py`, `fn_method/sphere/one_group.py` — then the four
   `converged=True` at `spectrum.py:819/872` and `moment_space.py:315/347` become
   `bool(res.converged)`, matching the exemplar at `spectrum.py:924`.
3. **The two Peierls Nyström power iterations** (`peierls_nystrom/slab.py:596`,
   `geometry.py:6450`) — same defect, second family, plus a sixth transcription of
   the predicate (`... and iteration > 5`).
4. **The dead local** at `fn_method/slab/flux_reconstruction.py:854/866`.
5. ⛔ **A secondary correctness bug found en route, independent of the contract
   question**: `fn_method/sphere/one_group.py:358-362`'s fallback discards scipy's
   unconverged-but-good iterate (`[M]` 8.4e-5) for the coarse-grid guess
   (`[M]` 7.0e-3) — **83× worse** — and then reports it as converged.
6. **Gate gaps** (the #340 "audit every gate" item, `derivations/` half):
   `tests/cross_method/test_polymorphism.py` consumes `CriticalSolution` at 5 sites
   and never asserts `.converged`; `FNSlabAdapter`/`FNSphereAdapter`
   (`tests/cross_method/adapters.py:55, 88`) omit `metadata["converged"]` entirely,
   for exactly the two producers that are least honest.

**Mutation-verification owed before any of the above is credited as fixed** (my own
standing rule): each new `converged` field needs a starved-budget gate that goes RED
— the probes in §Q1 are ready-made (`max_bisect=3` on the SE slab moves `d` by
1.4e-2; `max_bisect=1` on the FN sphere moves `R` by 7.0e-3). Note the **positive
control for such a gate must itself be budget-dependent**: a catastrophic mutation
of the kernel leaves a convergence-flag gate green by construction.
