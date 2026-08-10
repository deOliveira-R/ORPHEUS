# N6 commit 1 — the warning speaks for the level that FAILED

**Verification plan.** Author: test-architect, 2026-08-10.
Tree: `main` @ **`52650a86`** (`git status` clean of `orpheus/`; the working
tree carries uncommitted state only under `.claude/` and `scratch/`).
Canonical invocation: `.venv/bin/python -O -m pytest`, SERIAL.

Every number marked `[M]` was measured in THIS dispatch at that HEAD, with its
configuration stated. Numbers inherited from the dispatch brief are marked
`[M-brief]` and are NOT re-measured unless noted.

---

## ⛔ READ FIRST — one correction to the change's stated scope

The brief enumerates **three** facts moving from the caller to
`record.first_failure`: `budget`, `tolerance`, `knob`. That is **not enough for
the message to be coherent**, and the resulting half-carve is *worse than
today*.

`_warn_if_unconverged` (`orpheus/sn/solver.py:393-528`) reads the record at
**five sites** in its body, all of them `history.record` = the TOP level:

| line | read | today's source | must become |
|---|---|---|---|
| `:473` | `criterion = history.record.binding_criterion` | top | `failing.binding_criterion` |
| `:493` | `shortfall = history.record.min_iterations` | top | `failing.min_iterations` |
| `:496` | `history.record.n_iterations` | top | `failing.n_iterations` |
| `:498` | `history.record.n_iterations < shortfall` | top | `failing.n_iterations` |
| `:501` | `history.record.n_iterations` | top | `failing.n_iterations` |
| `:520` | `budget_name`, `budget`, `tol` (params) | caller (top) | `failing.budget_name`, `.budget`, `.binding_criterion.tolerance` |
| `:455` | `history.converged` — **the GUARD** | top | **UNCHANGED in commit 1** |
| `:520` | `where` (param) | caller | **UNCHANGED** — it names the public ENTRY, not a level |

`criterion` at `:473` feeds `distance` (`:474-478`), `rate` (`:480`),
`projected` (`:481`) and the `cleared` branch (`:482`). If only the three
briefed facts move, the emitted message on the discriminating fixture reads

> `hit max_inner=20 without reaching tol=1.000e-10 (last dphi 6.308e-02) …
> rho=0.153338 … set max_inner=17`

— the inner's knob and budget welded to the **outer's** criterion, rate and
projection. `[M]` those seven numbers are pairwise distinct on that fixture
(table §2.1), so this is a real, emittable lie, not a hypothetical.

⟹ **The carve is: ONE named intermediate, and every body read goes through it.**

```python
# the level that actually failed — children first, so a starved inner is
# reported rather than the outer it starved (convergence.py:842-856)
failing = history.record.first_failure or history.record
```

That is also the `coding-elegance` Pattern-3 spelling (a named domain
intermediate) and Pattern-2 (one navigation, not five). Gate **G2** below is
built specifically to red on the partial carve; mutations **M2/M3** are the
two halves of it.

`or history.record` is a total-function fallback, not a live branch: behind
`if history.converged: return`, `first_failure` is provably non-`None`
(`first_failure is None ⟺ fully_converged`, and `fully_converged ⟹ converged`).
It is also non-`None` under commit 2's guard. See §7 B4 — it is a **documented
non-catcher**, not a gap.

---

## 1. Claim layer and pillar, per the `vv-principles` gate

Every row below is a **software-contract claim** about a diagnostic message —
not a convergence-order, flux-shape or eigenvalue claim. Consequences:

- **Pillar: none of the three.** There is no closed form, no MMS and no
  semi-analytical reference for "which level's budget does this string name".
  The truth-source is a **hand-authored expectation** stated in the test,
  cross-checked against an **independently-navigated** path into the record
  (`record.children[0]`) — structurally independent of the SUT's own
  navigation (`first_failure`), which is precisely the claim.
- **Level: `foundation`** for every new row (a software invariant with no
  theory-page `:label:`), matching the file it joins. **No `verifies(...)`**
  (`feedback_vv_tagging`).
- **The Cardinal Rule does not bind here** — no row makes an eigenvalue claim.
  The eigenvalue fixtures below are 2-group anyway (`get_mixture("A","2g")`);
  `keff` is never asserted.
- ⚠ The one place a physics claim *could* sneak in is "the discriminating
  fixture's `keff`". `[M]` it is `1.8740647823989087` on a solve whose three
  inners are all TRUNCATED — an arbitrary mid-descent number. **Assert nothing
  about it.** (This is the `convergence.py:18-39` lesson applied to my own
  fixture.)

---

## 2. The discriminating fixtures

Commit 1's whole risk lives in ONE structural fact: today's four audibility
fixtures are all **1-level trees**, where `first_failure is record` by object
identity, so no input can make the two sourcings differ. `[M-brief]` and
re-confirmed `[M]` (`scratch/probe_340_n6_message_parity.py`, run at
`52650a86`): `solve_sn_fixed_source(..., max_inner=50)` reports
`ff is record: True`, `children 0`.

Two new fixtures, at two tiers. **Both are needed**; neither substitutes for
the other (`vv` Mode 11).

### 2.1 `N6-S` — the synthetic two-level record (Tier A, the keystone)

Drives `_warn_if_unconverged` directly, exactly as three committed gates
already do (`test_convergence_contract.py:435, 487, 538`). Cost ≈ **0 s**.

```python
inner = IterationRecord(
    label="inner(within-group)",
    criteria=(StoppingCriterion(
        name="residual",
        trajectory=tuple(3.7e-2 * 0.9**k for k in range(30)),
        tolerance=1e-10),),
    budget=30, iterations_run=30,
    budget_name="max_inner",              # ← the NEW field
)
outer = IterationRecord(
    label="outer(power-iteration)",
    criteria=(
        StoppingCriterion(name="dk",   trajectory=tuple(1e-2*0.5**k for k in range(6)), tolerance=1e-9),
        StoppingCriterion(name="dphi", trajectory=tuple(1e-3*0.5**k for k in range(6)), tolerance=1e-7),
    ),
    budget=6, iterations_run=6, min_iterations=3,
    budget_name="max_outer",              # ← the NEW field
    children=(inner,),
)
```

`[M]` measured at `52650a86` (constructed live; `budget_name` omitted since the
field does not exist yet — every other number is exact):

| fact | OUTER (today's source) | INNER = `first_failure` (N6's source) | distinct? |
|---|---|---|---|
| knob | `max_outer` | `max_inner` | ✅ |
| budget | `6` | `30` | ✅ |
| tolerance | `1.000e-09` | `1.000e-10` | ✅ |
| binding criterion | `dk` | `residual` | ✅ |
| last value | `3.125e-04` | `1.679e-03` | ✅ |
| rate `rho` | `0.500000` | `0.900000` | ✅ |
| `projected_iterations()` | `25` | `189` | ✅ |
| `converged` | `False` | `False` | — |
| `fully_converged` | `False` | — | — |

`[M]` **today's emitted message on `N6-S`** (via `_warn_if_unconverged(…,
where="solve_sn", budget_name="max_outer", budget=6, tol=1e-9)`):

```
solve_sn: hit max_outer=6 without reaching tol=1.000e-09 (last dk 3.125e-04).
Returning a BEST-EFFORT iterate — it is mid-descent, not the converged answer.
At the rate observed so far (rho=0.500000) this needs about 25 iterations:
set max_outer=25. Or read `solution.history.converged` and handle it. …
```

**Required after N6** (all seven facts flipped, `where` unchanged):

```
solve_sn: hit max_inner=30 without reaching tol=1.000e-10 (last residual 1.679e-03).
… (rho=0.900000) … set max_inner=189. …
```

⭐ Why the trajectories are exact geometrics: `rho` and `projected` are then
analytically determined (`0.9`, `0.5`) rather than solve-dependent, so the row
is platform-stable. But the gate still **recomputes** them from
`inner.rate` / `inner.projected_iterations()` rather than hard-coding —
hard-coding `189` would break on any future change to `_RATE_FIT_TAIL_FRACTION`
and would be a false red. The literals go in an **activation guard** instead
(§3, G2 leg 0).

### 2.2 `N6-E` — the end-to-end nested-failure solve (Tier B)

The Tier-A row cannot see whether **production** builds a nested record, nor
whether the SN layer stamps the knob onto the child (`vv` Mode 11: a synthetic
record proves the consumer's routing, never the producer's).

```python
solve_sn(
    {0: get_mixture("A", "2g")},
    tuple(AxisMesh(edges=np.linspace(0.0, e, n + 1), bc_low=_REFL, bc_high=_REFL)
          for e, n in zip((2.0, 1.5), (4, 3))),          # 2-D, all-reflective
    Quadrature.level_symmetric(sn_order=4),
    max_outer=3, keff_tol=1e-12, flux_tol=1e-12,
    max_inner=20, inner_tol=1e-10,
)
```

`[M]` **0.18 s**, deterministic (`max_outer=3 == MINIMUM_OUTER_ITERATIONS`, so
the loop always runs exactly 3 outers):

| fact | OUTER | `first_failure` = `children[0]` | distinct? |
|---|---|---|---|
| label | `outer(power-iteration)` | `inner(source-iteration)` | ✅ |
| budget | `3` | `20` | ✅ |
| status | `TRUNCATED` (3/3) | `TRUNCATED` (20/20) | — |
| binding | `dphi` @ `1e-12`, last `6.308e-02` | `residual` @ `1e-10`, last `2.971e-01` | ✅ |
| rate | `0.153338` | `0.962217` | ✅ |
| projected | `17` | `586` | ✅ |
| `converged` / `fully_converged` | `False` / `False` | — | — |

All three children are TRUNCATED (`20/20`), so `first_failure` is
`children[0]` deterministically.

⚠ **Why not the brief's d=3 fixture.** `[M]` the brief's
`(1.0,2.0,3.0)/(3,4,5)`, LS-S4, `max_outer=3, max_inner=200` fixture is equally
discriminating (`max_outer=3` vs `max_inner=200`; `dphi` 1e-12 vs `residual`
1e-8; rho `0.547396` vs `0.984928`; proj `27` vs `838`) but costs ~200× more.
Keep the 2-D one for the routine gate. The d=3 one is the *scientifically*
interesting case (its `keff` is right to 2.5e-11 despite 3/3 truncated inners —
the "benign truncation" that motivates N5) and belongs in the plan's prose, not
in a routine gate.

### 2.3 `N6-K` — the knob-stamping sweep (Tier B, producers)

`[M]` all six entry × inner-solver combinations run on a **4-cell 1-D**
reflective slab, `Quadrature.gauss_legendre(n_ordinates=4)`, 2 groups, in
**0.01–0.20 s each** (total < 0.5 s), and all CONVERGE (no warning needed —
the field is on the record whether or not it is emitted):

| entry | inner_solver | `[M]` cost | top record | children |
|---|---|---|---|---|
| `solve_sn_fixed_source` | `source_iteration` | 0.04 s | `inner(source-iteration)` | 0 |
| `solve_sn_fixed_source` | `krylov` | 0.01 s | `inner(gmres)` | 0 |
| `solve_sn` | `source_iteration` | 0.20 s | `outer(power-iteration)` | 3 × `inner(source-iteration)` |
| `solve_sn` | `krylov` | 0.04 s | `outer(power-iteration)` | 3 × `inner(gmres)` |
| `solve_sn_adjoint` | (default) | 0.13 s | `outer(power-iteration)` | 3 × `inner(source-iteration)` |
| `solve_sn_adjoint_fixed_source` | (default) | 0.03 s | `inner(source-iteration)` | 0 |

⚠ `solve_sn_adjoint_fixed_source` takes `detector_response` shaped
`(ng, *spatial)` = `(2, 4)`, **not** the angular `(N, ng, *spatial)` the
fixed-source entry takes. `[M]` passing the angular shape raises
`ValueError: … detector_response shape (4, 2, 4) != (ng, *spatial) = (2, 4)`.

⭐ `solve_sn_adjoint` **does** carry children — `[M]` 3 ×
`inner(source-iteration)`. The plan's F7′/§1b claim that "the adjoint
eigenvalue path's inner trajectory is discarded inside `KEigenvalue`" is
**already discharged** by N2; do not design around it.

---

## 3. The gate set

Nine rows. All `@pytest.mark.foundation`, all in
`tests/sn/solve/test_convergence_contract.py` unless stated, all under
`class TestTruncationIsAudible` unless stated.

| # | gate | tier | the ONE claim it pins | mutation that reddens it | `[M]` cost |
|---|---|---|---|---|---|
| G1 | `test_the_message_speaks_for_the_FAILING_LEVEL_not_the_caller_s` | A (`N6-S`) | every fact in the message is sourced from `first_failure` | **M1** (`failing := history.record`) | ~0 s |
| G2 | (same row, 7 legs + activation) | A | *each* fact individually | **M2 / M3** (partial re-point) | — |
| G3 | `test_a_real_nested_solve_reports_its_STARVED_INNER` | B (`N6-E`) | production builds the tree AND stamps the child's knob | M1, M4, M5, M7 | 0.18 s |
| G4 | `test_every_entry_stamps_its_own_caller_facing_knob` (new class) | B (`N6-K`) | the knob on every level IS a parameter of the entry that returned it | **M4 / M5 / M6** | < 0.5 s |
| G5 | `test_a_converged_outer_on_a_starved_inner_is_STILL_SILENT` | B (`N6-E`, generous outer) | **commit 1's scope boundary** — the guard did not move | **M8** (guard flipped early) | ~0.2 s |
| G6 | `test_the_tree_wide_truncation_is_audible` — `xfail(strict=True)` | B | the commit-2 todo marker; XPASSes the day the guard flips | (none — it is the marker) | ~0.2 s |
| G7 | `test_a_leaf_message_is_UNCHANGED_by_the_re_sourcing` | B (existing `_fixed_source(50)`) | the measured slice is behaviour-neutral | M4 (knob only — see §7 B2) | 0.38 s |
| G8 | `test_the_knob_name_is_required_to_be_NAMED` (in `tests/numerics/test_iteration_record.py`) | unit | the field's construction invariant + its positive leg | **M9** | ~0 s |
| G9 | `test_the_retired_parameters_are_GONE_from_the_signature` | unit | the retirement half — a call site left on the old spelling is a twin | M-R (re-add a defaulted `tol=`) | ~0 s |

### G1/G2 — the keystone (Tier A)

One test, one `pytest.warns`, and then **legs in pairs**: each fact asserted
PRESENT with the inner's value and ABSENT with the outer's. The absent leg is
what makes it a gate rather than a self-consistency check (`lessons` §4 — "ask
which side is the thing under test; if the answer is *both*, it is not a
gate").

```
leg 0  ACTIVATION (fixture-drift guard, runs BEFORE the warn):
       assert outer.first_failure is inner
       assert inner.budget            != outer.budget
       assert inner.binding_criterion.tolerance != outer.binding_criterion.tolerance
       assert inner.binding_criterion.name      != outer.binding_criterion.name
       assert inner.rate                        != outer.rate
       assert inner.projected_iterations()      != outer.projected_iterations()
       assert inner.budget_name                 != outer.budget_name
leg 1  knob      "max_inner=30" in msg     and "max_outer" not in msg
leg 2  budget    (covered by leg 1's "=30") and "=6 " not in msg
leg 3  tol       "tol=1.000e-10" in msg    and "1.000e-09" not in msg
leg 4  criterion "last residual" in msg    and "dk" not in msg
leg 5  rate      f"rho={inner.rate:.6f}" in msg  and "rho=0.500000" not in msg
leg 6  advice    f"set max_inner={inner.projected_iterations()}" in msg
                 and "set max_outer=" not in msg
leg 7  where     msg.startswith("solve_sn:")   ← `where` must NOT re-point
```

⭐ Leg 0 is not ceremony. It is the **anti-Mode-12-at-the-fixture** guard: if a
future edit makes any inner fact equal its outer counterpart, the
corresponding leg silently stops discriminating and the row goes green forever.
Leg 0 makes that a red instead. (`vv` anti-#23's discipline, applied to a
sourcing gate rather than an invariance gate.)

⚠ Leg 4's `"dk" not in msg` is safe: `[M]` the emitted template contains no
other `dk` substring. Leg 2's `"=6 "` must be spelled with the trailing space
(`hit max_inner=30 …` vs `hit max_outer=6 without`) — verify against the
measured string in §2.1, do not guess.

### G3 — the same claim, end to end (Tier B)

`N6-E`, one `pytest.warns`, and the reference navigated **by hand** through
`sol.history.record.children[0]` — a different route into the same object than
production's `first_failure`. Legs:

```
assert sol.history.record.first_failure is sol.history.record.children[0]   # activation
child = sol.history.record.children[0]
assert f"{child.budget_name}={child.budget}" in msg          # "max_inner=20"
assert "max_outer" not in msg
assert f"tol={child.binding_criterion.tolerance:.3e}" in msg # "tol=1.000e-10"
assert f"last {child.binding_criterion.name}" in msg         # "last residual"
assert f"set max_inner={child.projected_iterations()}" in msg
assert msg.startswith("solve_sn:")
```

⭐ **Do not assert `keff`** — §1.

### G4 — the knob is a real knob on the entry that emitted it

New class, `TestEveryEntryStampsItsCallerFacingKnob`, parametrized over the six
`N6-K` rows. The reference is `inspect.signature(<entry>).parameters` — an
external, structurally-independent source (it is the public API, not the
record). This is the ONLY gate that can catch a **silent default**: if a
producer forgets to stamp, the record carries `"max_iter"`, which `[M]` is a
parameter of **none** of the four public entries:

| entry | budget knobs `[M]` |
|---|---|
| `solve_sn`, `solve_sn_adjoint` | `max_outer`, `max_inner` |
| `solve_sn_fixed_source`, `solve_sn_adjoint_fixed_source` | `max_inner` |

```
for record in sol.history.record.walk():
    assert record.budget_name in inspect.signature(entry).parameters, (
        f"{entry.__name__} returned a level advising `set "
        f"{record.budget_name}=…`, which is not one of its parameters"
    )
# and the exact expectation, so a SWAP (M5) cannot pass the membership test
assert sol.history.record.budget_name == expected_top      # max_outer | max_inner
assert all(c.budget_name == "max_inner" for c in sol.history.record.children)
```

⚠ The membership leg alone is **blind to M5** on `solve_sn`, because both
`max_outer` and `max_inner` are parameters of it. The exact-expectation leg is
mandatory. (`lessons` §6 OPTIONAL→MANDATORY-BINDING: *design the battery around
the SWAP*.)

⛔ **This check must NOT be added to `IterationRecord.__post_init__` or to
`_warn_if_unconverged`.** Three committed gates pass `where="probe"`, which is
not a public entry; a production guard would break them, and it would assert
more than the callers promise (`vv` anti-#16). The check lives in a test that
runs the real entries. See §7 B3.

### G5/G6 — the commit-1 scope boundary, stated twice on purpose

**PREDICTED (not yet measured — M8 is the measurement, and it must be run
first): nothing in the tree can currently detect an accidental flip of the
guard to `fully_converged`.** The derivation rests on three `[M]` facts:
(i) `pyproject.toml:86` sets `addopts = "--import-mode=importlib"` and **no**
`filterwarnings` / `-W`, so an emitted `ConvergenceWarning` fails nothing;
(ii) every `ConvergenceWarning` assertion in the tree lives in the single file
`tests/sn/solve/test_convergence_contract.py` (grep over `tests/`);
(iii) its only silence gate, `test_converged_solve_is_SILENT` (`:370`), rides a
**leaf** where `converged ≡ fully_converged`. The `≥11` currently-green solves
with `converged=True, fully_converged=False` (`[M-brief]`:
`tests/sn/eigenvalue/*` plus 5 rows of
`tests/sn/regression/test_dd_regression.py`) would therefore begin emitting and
every one of them would stay green.

- **G5** (positive, must be **INVERTED** by commit 2): the `N6-E` fixture with
  a generous outer and a starved inner emits nothing. Its docstring must say,
  in its first line, *"⛔ commit 2 INVERTS this row — it pins the SCOPE of
  commit 1, not the correctness of the silence; see #340 N5/N6."*
- **G6** (`xfail(strict=True, reason="#340 N6 commit 2 — gated on N5's residual
  certificate")`): the same fixture, asserting the warning IS emitted. It
  xfails today and **XPASSes strictly** the day the guard flips, forcing both
  rows to be revisited in the same change. This is the project's
  "`xfail(strict=True)` set IS the todo list" idiom.

The pair is deliberately redundant: G5 catches an *accidental* flip now; G6
makes the *deliberate* flip self-announcing. Neither alone does both.

⚠ **G5/G6 need a generous-outer-starved-inner fixture, and `N6-E` as specified
is not one** (its outer also fails). Build it as `max_outer=200, keff_tol=1e-8,
flux_tol=1e-6, max_inner=<small>, inner_tol=1e-10` on the same 2-D mesh, and
add an **activation guard** `assert sol.history.converged is True and
sol.history.fully_converged is False` — without it the row degenerates into
`test_converged_solve_is_SILENT` and proves nothing. ⛔ `[M]` the existing
`_eigenvalue(max_outer=200, keff_tol=1e-8)` is **not** such a fixture: its
three inners all CONVERGED (721/670/2 of 1635) and `fully_converged is True`.
The implementer must tune `max_inner` down until the activation guard holds and
**record the measured triple in the docstring**.

### G7 — the leaf is unchanged (and what that can and cannot see)

Keep the existing `_fixed_source(max_inner=50)` rows. `[M]` today's message:

```
solve_sn_fixed_source: hit max_inner=50 without reaching tol=1.000e-13
(last residual 1.952e-03). … (rho=0.956890) … set max_inner=586. …
```

`[M]` `ff is record: True`, `children 0`. See §7 B2 for the exact scope of
what a leaf row can still catch after the change (**the knob name, and nothing
else** — the other six facts become an object compared with itself).

### G8 — the field's own construction invariant

Home: `tests/numerics/test_iteration_record.py`, class
`TestRecordConstructionEnforcesItsInvariants` (which already pairs every guard
with a positive leg, `vv` anti-#11).

```
def test_an_unnamed_budget_knob_is_refused():
    with pytest.raises(ValueError, match="budget knob"):
        IterationRecord(label="level", budget_name="")

def test_the_default_knob_is_the_primitive_s_own_parameter_name():
    assert IterationRecord(label="level").budget_name == "max_iter"
```

See §6 for why `""` is the **only** thing it may refuse.

### G9 — the retirement half

```
sig = inspect.signature(_warn_if_unconverged).parameters
assert set(sig) == {"history", "where"}, (
    f"the caller-passed level facts must be GONE, not defaulted: {sorted(sig)}"
)
```

A defaulted `budget_name=…` left on the signature is the twin that reopens the
defect: every call site still compiles, and one of them silently keeps passing
the caller's top-level knob. A membership assertion (`"tol" not in sig`) is
weaker than the set equality — use the set.

---

## 4. Which EXISTING gates change

`[M]` inventory: **every** `ConvergenceWarning` assertion in the tree lives in
`tests/sn/solve/test_convergence_contract.py` (grep over `tests/`). Nothing in
`tests/numerics/`, `tests/cp/`, `tests/moc/`, `tests/diffusion/` touches it.

| gate | verdict | why |
|---|---|---|
| `:365 test_starved_solve_warns` | **no change** | leaf; asserts only that ≥1 warning fired |
| `:370 test_converged_solve_is_SILENT` | **no change** in commit 1 | leaf; `converged ≡ fully_converged` there. ⚠ It is *not* a catcher for the guard flip — that is G5's job |
| `:378 test_the_message_carries_budget_tolerance_and_distance` | **no change; becomes a free M4 catcher** | leaf. Asserts `"max_inner=50" in msg`; after the carve that string comes from `SourceIteration`'s stamp, so a forgotten stamp (`"max_iter=50"`) reds it |
| `:393 test_the_message_names_the_BUDGET_TO_SET_not_just_a_direction` | **no change; free M4 catcher** | leaf. Asserts `f"set max_inner={projected}"` |
| `:424 test_a_NON_CONTRACTING_solve_is_told_no_budget_will_help` | ⛔ **BREAKS — `TypeError`** | `:451-454` calls `_warn_if_unconverged(stalled, where="probe", budget_name="max_inner", budget=40, tol=1e-10)`. **Mechanical migration**, not a re-baseline: move `budget_name="max_inner"` into the `IterationRecord(...)` at `:438-448` (it already carries `budget=40` and `tolerance=1e-10`), drop the three kwargs. All four assertions survive verbatim |
| `:459 test_it_speaks_for_the_criterion_that_FAILED_not_one_that_cleared` | ⛔ **BREAKS — `TypeError`**, **and its message silently changes** | Same migration (`budget_name="max_outer"` onto the record at `:491-506`). ⭐ **Then re-baseline one thing**: today the message prints `tol=1.000e-09` (the caller's `keff_tol`) beside `last dphi` — whose tolerance is `1e-6`. That is the F2 disease inside the diagnostic itself. After N6 it prints `tol=1.000e-06`. The test does **not** assert `tol=` today, so the correction is invisible. **Add the two legs** `assert "tol=1.000e-06" in msg` and `assert "1.000e-09" not in msg` — this is a re-baseline in the honest sense (the old string was wrong; the gate simply never looked) |
| `:529 test_it_names_the_LOOP_when_every_criterion_cleared` | ⛔ **BREAKS — `TypeError`** | Same migration (`budget_name="max_outer"` onto the record at `:540-548`). Assertions survive: `first_failure is record` on this leaf, so `min_iterations`/`n_iterations` read identically, and `"set max_outer=" not in msg` still holds because the record names `max_outer` |
| `:561 test_it_is_escalatable_to_an_error` | no change | |
| `:574 test_the_published_escalation_flag_actually_parses` | no change | |
| `tests/numerics/test_iteration_record.py:595 test_no_convergence_verdict_is_a_constructor_argument` | **no change, and it must stay green** | `budget_name` is an **observation**, not a verdict — the same category as `label`, `budget`, `iterations_run`. The class docstring already draws that line (`convergence.py:683-686`: *"no convergence judgement may be stored; measured data plainly must be"*). If the implementer is tempted to name the field something verdict-like, this gate is the tell |
| `tests/sn/solve/test_coupled_solve_certificate.py:158` | no change | keyword construction; a defaulted field is inert |
| ~40 other `IterationRecord(...)` constructions in `tests/` | no change | `[M]` **zero** positional constructions tree-wide (`grep` for `IterationRecord(` + first arg), so inserting the field after `budget` is safe |

### 4.1 ⛔ A SIXTH consumer, and it is already DEAD

`scratch/mutate_convergence_contract.py` is **tracked** (`git ls-files`
confirms) and is the campaign's own acceptance instrument for this contract.
`[M]` it **cannot import at HEAD**: line 24 reads
`REAL_CLAIMS = sv._claims_convergence`, and `_claims_convergence` was retired on
2026-08-09 — `[M]` `hasattr(orpheus.sn.solver, "_claims_convergence") is False`.

It is also stale in two ways N6 makes worse:

- `ma1()` (`:51-57`) builds `PowerIterationOutcome(..., converged=True)` — a
  keyword that no longer exists, and whose absence has its own gate
  (`test_convergence_contract.py:300`). The mutation now raises `TypeError`
  inside every solve, i.e. it is an **anti-#18 over-powered** mutation whose
  reds are attributable to a crash, not to the property.
- `ma4()` (`:91`) and `ma5()` (`:106`) declare replacements with the signature
  `(history, *, where, budget_name, budget, tol)`. After commit 1 the five call
  sites pass only `where`, so both replacements raise `TypeError: missing 3
  required keyword-only arguments` — again crash-reds, not property-reds.

⟹ **Repairing this harness is part of commit 1**, not optional cleanup: without
it there is no instrument to run the §5 battery on, and `lessons` §2 is explicit
that a battery's verdict is worthless until its own control passes. Repair =
delete the `REAL_CLAIMS` binding and `ma2`, retire `ma1` (its defect is now
unspellable and has a dedicated gate), and re-sign `ma4`/`ma5` to
`(history, *, where)`.

### 4.2 Doc debt landing in the same commit

- `docs/theory/methods/sn/solver.rst:691-693` — *"A truncated exit emits
  `ConvergenceWarning`, naming the budget that ran out, the tolerance
  missed…"*. Two edits: (a) it becomes *"naming the budget of the level that
  ran out"* (N6's correction); (b) the sentence is **present-tense FALSE
  today** for a tree — a truncated *inner* under a converged outer emits
  nothing. Per the standing "a falsified doc is a BUG — fix on sight" rule that
  scoping correction lands in commit 1 even though the behaviour fix is commit
  2, with the `#340 N6 commit 2` pointer.
- `orpheus/sn/solver.py:403-406` — the docstring's *"so a truncated solve
  announces itself once, in one voice, wherever it came from"* is the same
  false claim, at the emission point. Same treatment.
- `orpheus/numerics/convergence.py:96-137` (the "What a stalled solve owes its
  caller" section) — add the knob field to the record's described contract.

---

## 5. Mutation battery

**Method.** In-process `unittest.mock.patch` / `monkeypatch` only. ⛔ **NEVER
`git checkout`/`restore`/`stash`** — the tree carries uncommitted state
(`process-discipline`). Run the whole contract file plus the two numerics
files:

```
SUITE = ["tests/sn/solve/test_convergence_contract.py",
         "tests/numerics/test_iteration_record.py",
         "tests/numerics/test_power_iteration_record.py"]
```

`[M]` cost of the unmutated suite must be measured as **M0** and every mutation
budgeted off the MUTATED cost, not the green one (`lessons` §2 — a mutated run
that destroys convergence gets *slower*).

⚠ **The harness must assert its own installation.** After the §4.1 repair, add
to each mutation a positive bite-check (e.g. call the patched symbol once and
assert the observable moved) before reading the suite's colour. A "0 caught"
verdict from an uninstalled patch is the flattering failure (`vv` anti-#17).

| id | mutation | in-class? | expected reds | catcher |
|---|---|---|---|---|
| **M0** | **POSITIVE CONTROL** — `_warn_if_unconverged := lambda *a, **k: None` | n/a | **≥ 5** in the contract file (`test_starved_solve_warns`, `…carries_budget…`, `…names_the_BUDGET_TO_SET…`, `…NON_CONTRACTING…`, `…speaks_for_the_criterion…`, `…names_the_LOOP…`, `…escalatable…`, G1, G3) | if this reports < 5, the harness is broken; stop |
| M1 | `failing = history.record` (drop the `first_failure` navigation) | ✅ | **exactly G1 + G3** (2). NOT the leaf rows | the whole point of commit 1 |
| M2 | budget/tol/knob from `failing`; `criterion` left on `history.record` | ✅ | G1 legs 3–6, G3's tol/criterion/advice legs | the brief's literal 3-fact reading |
| M3 | `criterion` from `failing`; budget/tol/knob left on `history.record` | ✅ | G1 legs 1–3, G3's budget/knob legs | the inverse partial |
| M4 | `SourceIteration.solve` builds its record **without** `budget_name` (falls to the default) | ✅ | `:378`, `:393`, G3, G4 (all SI rows), G7 | the silent-default class |
| M5 | **SWAP** — SI stamps `"max_outer"`, `power_iteration` stamps `"max_inner"` | ✅ | G4's exact-expectation legs (all 6 rows), `:378`, `:393`, G1?no, G3 | ⭐ the mutation most likely to survive a membership-only check |
| M6 | `power_iteration` builds its record without `budget_name` | ✅ | G4's eigenvalue rows only (3 of 6) — **and NOT the leaf rows**; the asymmetry is the evidence the two producers are independently gated | producer discrimination |
| M7 | `IterationRecord.first_failure` patched to search **self before children** | ✅ | G1, G3, plus `tests/numerics/test_iteration_record.py:403-412` (`…DEEPEST failure`) | the ordering the plan's own author misread once |
| M8 | the guard at `:455` → `if history.fully_converged: return` | ✅ | **G5 only** (and G6 XPASSes strictly). ⛔ **Run this one FIRST, on the UNMUTATED tree, before writing G5** — *predicted* **0 reds** (derivation in §3 G5/G6). The measurement, not the argument, is what justifies G5; if it reds something, find out what and re-scope | the accidental-flip catcher |
| M9 | drop the `budget_name` non-empty check from `__post_init__` | ✅ | G8 negative leg only (1) | the field's guard |
| M-R | re-add `tol: float = 1e-8` to `_warn_if_unconverged`'s signature (a "harmless" default) | ✅ | G9 only | the retirement twin |
| **M-X** | ⚠ **anti-#18 NEGATIVE EXAMPLE, not coverage** — `first_failure := None` (no fallback) | ❌ out of class | reds every warning row via `AttributeError` | **do not count.** It breaks a type contract, not the sourcing property. Its reds prove nothing about which level the message speaks for |
| **M-N** | ⚠ **documented NON-CATCHER probe** — delete the `or history.record` fallback | ✅ | **expected 0** | proves the fallback is unreachable behind the `converged` guard (§7 B4). A 0 here is the *correct* result and must be recorded as such, not chased |

**Scale check (anti-#18).** Every expected red count above is ≤ 9 and every
reddened row has a direct textual view of the property. If any mutation reds
dozens of rows across `tests/sn/verification/` or `tests/sn/eigenvalue/`, it is
breaking something other than the message — investigate before crediting.

---

## 6. The knob field — and what its `__post_init__` may refuse

### 6.1 `vv` anti-#16: the CALLERS' boundary values, MEASURED FIRST

`[M]` three production sites construct an `IterationRecord`
(`orpheus/numerics/iteration.py:789`, `:1101`; `orpheus/numerics/eigenvalue.py:467`).
`[M]` their consumers' caller-facing knob spellings:

| producer | its own parameter | consumer | the knob THAT consumer exposes |
|---|---|---|---|
| `SourceIteration.solve` | `max_iter` | `sn/solver.py:1003, 1026` (SN fixed-source) | `max_inner` |
| " | " | `sn/solver.py:2880` (adjoint fixed-source) | `max_inner` |
| " | " | `numerics/iteration.py:1293` (`KEigenvalue._inner`) | `max_inner` |
| " | " | `numerics/green_operator.py:277` | **`max_iter`** ← `GreenOperator(max_iter=…)` |
| `KrylovAcceleration.solve` | `max_iter` | `sn/solver.py:666` | `max_inner` |
| `power_iteration` | `max_iter` | `sn/solver.py:2320`, `cp/solver.py:903`, `moc/solver.py:137`, `diffusion/solver.py:413`, `numerics/iteration.py:1529` | `max_outer` (all five) |

⟹ the value set the callers actually promise is **`{max_inner, max_outer,
max_iter}`**, not `{max_inner, max_outer}`.

⛔ **A whitelist guard would break `GreenOperator`** — exactly the shape that
broke `GreenOperator(tol=0)` and the GMRES exact-breakdown path earlier in this
same campaign (`convergence.py:398-407`). **No whitelist. No enum. No
`Literal`.** The field is a free `str`.

### 6.2 What it MAY refuse — one thing

```python
if not self.budget_name:
    raise ValueError(f"{self.label}: the budget knob must be named")
```

Rationale: it mirrors the existing `label` invariant (`convergence.py:704-705`)
one-for-one, it is the same category of fact (a name that appears **verbatim**
in a human-facing message), no caller passes `""`, and without it the message
degrades to `hit =50 without reaching tol=…`. Cheap, unreachable, honest.

It may **not** refuse: a knob that is not a public-entry parameter (see §7 B3),
an empty knob on a `budget == 0` direct level (`IterationRecord(label="direct(lu)")`
must keep constructing — `tests/numerics/test_iteration_record.py:479, 573`),
or `None` (nothing passes it; `str` + pyright is the right instrument).

### 6.3 Default, name, and position

- **Default `"max_iter"`, not required.** `[M]` ~40 test constructions and 1
  production-adjacent construction (`test_coupled_solve_certificate.py:158`)
  would otherwise need editing, and `GreenOperator` gets the *correct* answer
  from the default for free. The residual risk — a producer silently forgetting
  — is converted into a red by **G4**, which is why G4 is not optional.
- **Spelling: `budget_name`.** Keep the retiring parameter's exact name so the
  retirement is a **move**, greppable end-to-end
  (`feedback_naming_consistency_greppable`). `budget_knob` reads slightly
  better and costs the grep continuity; not worth it.
- **Position: immediately after `budget`.** `[M]` zero positional
  constructions tree-wide, so this is safe; adjacency is what makes the pair
  (`budget`, `budget_name`) read as one fact.
- ⚠ **`report()` should name it too.** `convergence.py:936-938` currently emits
  `needs ~{needed} iterations for tol {…}` with no knob. Once the record knows
  the knob, the report saying `set max_inner=189` is the consistent thing and
  costs one f-string. **Out of scope for commit 1** — but if the implementer
  does it, `tests/numerics/test_iteration_record.py:727-753` gains a leg and
  this plan's G-list gains a row. Decide explicitly; do not let it drift in.

---

## 7. Blind spots — what CANNOT be gated, and the structural reason

**B1 — commit 1 cannot gate #340's headline defect.** F1 is *"a converged outer
on a starved inner is silent"*. The guard at `:455` is the **single** emission
decision, and the ratified split leaves it on `history.converged`. No gate in
commit 1 can therefore assert audibility for the F1 case. The best available
instrument is G6's `xfail(strict=True)`, which asserts nothing today but cannot
be forgotten tomorrow. Recording this is not pessimism — it is the honest
statement that **commit 1 fixes the message, not the silence**.

**B2 — the leaf "no change" claim is a THEOREM for six of seven facts, and a
gate for one.** After the carve, a leaf has `first_failure is record` — *object
identity* (`[M]` confirmed: `ff is record: True` on `_fixed_source(50)`). So
`budget`, `tolerance`, `criterion`, `last`, `rate` and `projection` are the
same object's attributes read twice; **no input exists** that could make the
two sourcings differ. G7 is therefore a `X == X` comparison for those six
(`coding-standards`' single-sourcing demotion; `vv` anti-#22). The ONE fact it
still sees is the **knob name**, which genuinely moved from a caller literal to
a record field. G7's docstring must say exactly that, or an audit will count it
as coverage of the re-sourcing.

**B3 — "the knob is a real knob" cannot be a production invariant.** Three
committed gates drive `_warn_if_unconverged` with `where="probe"`
(`:452, :512, :554`), and `GreenOperator`'s knob is `max_iter`. A
`__post_init__` or in-function check that the knob is a parameter of `where`
would red all three and break `GreenOperator` — asserting more than the callers
promise. The check is only expressible as a **test that runs the real entries**
(G4), which means the invariant holds for the four SN entries and is *unstated*
everywhere else.

**B4 — the `or history.record` fallback is unreachable, provably.**
`first_failure is None ⟺ fully_converged` (gated at
`tests/numerics/test_iteration_record.py:377-379`), and `fully_converged ⟹
converged`, so behind `if history.converged: return` the fallback never fires.
M-N measures this (expected **0 reds**) and the 0 is the correct answer. Keep
the fallback anyway — it is total-function hygiene and it survives commit 2's
guard flip unchanged — but do not claim coverage for it.

**B5 — the `rate is None` branch (`:505-506`) has no cheap end-to-end
fixture.** It fires only when the failing level has < 2 usable residual points
(`_log_fit` returns `None`, `convergence.py:544-554`). `[M]` every fixture in
this plan has an estimable rate. It is reachable synthetically (a record with a
one-entry trajectory) and that is the only instrument; note in the row that no
real solve in the suite exercises it.

**B6 — CP / MoC / diffusion stamp the knob and nothing reads it.**
`_warn_if_unconverged` is SN-only (#340 F6: `CPResult` / `MoCResult` /
`DiffusionResult` do not even carry `converged`). Their `power_iteration`
records will carry `budget_name="max_outer"` with **zero consumers**. That is a
Mode-11-shaped gap — the field exists, the producer stamps it, no gate can
observe it end-to-end — and it is a #340 remainder (N4), not something N6 can
close. G4 covers the four SN entries and says so.

**B7 — nothing gates the message's PROSE.** Every gate above matches
substrings. A rewrite that keeps `max_inner=30`, `tol=1.000e-10` and
`set max_inner=189` while turning the surrounding sentence into nonsense stays
green. Structural reason: there is no reference for English. Accepted; the
existing `"BEST-EFFORT"` leg (`:391`) is the only prose anchor and it is enough.

---

## 8. Side-findings — measured while gating N6, OUT OF SCOPE, worth filing

**S1 — `exhausted_budget` is in the WRONG UNITS on the Krylov path.**
`[M]` `solve_sn_fixed_source(..., inner_solver="krylov", max_inner=5)` on the
d=3 absorber returns a record with **`budget=5`, `n_iterations=732`,
`status=CONVERGED`** — so `exhausted_budget = 732 >= 5 = True` on a *converged*
solve. Mechanism: `KrylovAcceleration.solve` sets `budget=self.max_iter`, which
is scipy's `maxiter` = **restart CYCLES** (`iteration.py:1055`, `restart=n_dof`
at `sn/solver.py:670`), while the trajectory counts **callbacks** (≈ cycles ×
`n_dof`). `truncated` is saved only by `converged`. Consequence for N6: if a
Krylov inner *does* truncate, the message will read `hit max_inner=5 … set
max_inner=<a callback count in the hundreds>` — a units mismatch in the very
advice N3 exists to make actionable. **File as a `module:sn` / `type:bug`
issue**; do not fix inside N6.

**S2 — the campaign's own mutation harness is dead at HEAD.**
`scratch/mutate_convergence_contract.py:24` — see §4.1. Tracked file,
`AttributeError` at import.

**S3 — a pre-existing incoherence the carve happens to fix, currently
ungated.** `[M]` on the `:459` fixture the message prints the caller's
`keff_tol` (`1e-9`) next to `last dphi`, whose tolerance is `1e-6`. Nothing
asserts it. §4 folds the correction into that gate's re-baseline.

---

## 9. Acceptance checklist for commit 1

1. `M0` control passes (harness alive after the §4.1 repair) — **before** any
   negative is believed.
2. `M8` run first, against today's gate set; record its red count (predicted
   **0**, derived in §3) as G5's justification. A prediction is not a
   measurement — run it.
3. G1–G9 written; each carries in its docstring the mutation id that reddens
   it and the `[M]` measurement of its activation guard.
4. All 11 counted mutations produce their predicted reds; `M-X` and `M-N`
   recorded as the two labelled non-catchers.
5. `grep -rn "budget_name\|budget=\|tol=" ` over the five call sites
   (`solver.py:2465, 2742, 2889, 3572, 3795`) returns **zero** survivors; G9
   pins it.
6. §4.2 doc edits landed (both present-tense-false sentences corrected).
7. `python -O -m pytest tests/sn/solve tests/numerics tests/sn/primitives`
   green; then the wide `tests/sn` run (`[M-brief]` ≈ 11 min) to confirm the
   ≥11 `converged=True / fully_converged=False` solves neither warn nor change.
