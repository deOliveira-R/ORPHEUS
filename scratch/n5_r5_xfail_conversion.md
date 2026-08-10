# R5 — retire the imperative `pytest.xfail()` (#340 convergence-honesty campaign)

**Defect (vv Mode 8, NINTH class).** `pytest.xfail(reason)` called imperatively
inside a test body raises `XFailed` at that line, so the row reports `x`
**unconditionally** and can NEVER report `XPASS`. The day the underlying issue is
fixed, nothing says so. The declarative `@pytest.mark.xfail(strict=True)` inverts
this: an unexpected PASS is a **FAILURE**, so the marker set is a self-retiring
todo list.

**Scope.** `grep -rn '^\s*pytest\.xfail(' tests/` → exactly 5 sites (the 3 further
raw-grep hits in `tests/numerics/test_advertised_degree_is_measured.py` are prose
*documenting* this anti-pattern — untouched).

---

## 1. The measurement that governs the step

Instrument: `scratch/_r5_noop_xfail.py` — a `-p` plugin rebinding `pytest.xfail`
to a no-op, which **asserts its own installation** (it prints the neutralisation
COUNT at `sessionfinish` and says loudly when the count is 0, so a
non-biting instrument cannot masquerade as a verdict). Every run below reported a
non-zero count with the expected reason strings.

Canonical invocation throughout: `.venv/bin/python -O -m pytest`, SERIAL.

### Per-row verdict

| # | Site | Test / row | Condition | **Verdict today** | Margin measured |
|---|---|---|---|---|---|
| 1 | `tests/sn/operators/test_apply_full_field_codomain.py:229` | `test_c5b_driver_reattach_recovers_kinf[sphere-krylov]` | `coord=="sphere" and inner_solver=="krylov"` (case is **2eg**) | ⭐ **PASSES — HEALED** | `rel = 1.184e-14` vs closed-form; gate asserts `rtol=1e-10` |
| 2 | `tests/sn/operators/test_removal_form_matvec_sweep.py:483` | `test_removal_form_kinf_independent_reference_2g` | UNCONDITIONAL stub (empty body) | **still blocked — LIVE todo** | n/a — nothing to measure (see §3) |
| 3 | `tests/sn/verification/analytical/test_kinf_homogeneous_tolerance.py:144` | `test_krylov_inner_matches_tight_si[4eg-sphere]` | `coord=="sphere" and ng_key=="4eg"` | ⭐ **PASSES — HEALED** | `rel = 4.462e-13`; gate asserts `< 1e-8` |
| 4 | `tests/sn/verification/analytical/test_kinf_homogeneous.py:163` | `test_kinf_homogeneous[sphere-4eg-krylov]` | `sphere ∧ 4eg ∧ krylov` | ⭐ **PASSES — HEALED** | `rel = 3.582e-15`; gate asserts `rtol=1e-10` |
| 5 | `tests/sn/verification/analytical/test_kinf_homogeneous.py:229` | `test_kinf_homogeneous_spectrum[sphere-4eg-krylov]` | `sphere ∧ 4eg ∧ krylov` | ⭐ **PASSES — HEALED** | eigenvector row; gate asserts `rtol=1e-9` |

**4 of the 5 exclusions have healed.** All four were the same claim: *"unpreconditioned
GMRES on sphere-4g exceeds the `max_inner=N` budget without converging; #200
tracks the block-inverse face preconditioner that re-enables it."*

### Why they healed — and it was NOT #200

`#200` ("SN: block-inverse preconditioner for Krylov on the typed AngularFlux
algebra") is **still OPEN**, and the preconditioner still does not exist —
`orpheus/sn/solver.py::_within_group_krylov` sets
`preconditioner = lambda q: q` when no DSA corrector is supplied. So the fix the
exclusions waited for has NOT landed.

What healed them lives in the GMRES **`restart`-sizing** lineage, in two steps,
neither of them #200:

| step | date | what it did |
|---|---|---|
| **ERR-053** | 2026-05-28 | removed the hardcoded `restart=min(50, full_size)` clamp, which truncated the Krylov subspace on any mesh with >50 unknowns. Its catalog entry works the arithmetic on **this** fixture family (a sphere, N=8 ordinates, 2 groups → 328 unknowns vs restart=50) |
| **#282 / #280 route (a)** | 2026-07-04 | sized `restart` from the FULL augmented ravel (bulk ⊕ trace ⊕ ψ½-seed) rather than the bulk alone. Its own gate `test_krylov_restart_covers_augmented_composite` records that the bulk-only count left *"the poorly-conditioned curvilinear-eigenvalue inner"* stalling and *"at a realistic outer cap it returns a WRONG keff"* — this symptom, on this geometry, named |

`_within_group_krylov`'s own docstring and call-site comment carry both:
*"``restart`` is sized to the FULL problem ``n_dof`` … the legacy ``min(50, …)``
clamp left GMRES structurally truncated (ERR-053)"* and *"A bulk-sized restart
re-truncates GMRES on the trace+seed DOFs (**the sphere Krylov stall**)."*

⚠ Which of the two was **decisive** for sphere-4eg is NOT discriminated — that
would need a production revert, which this dispatch is forbidden from doing. The
attributable fact is the one that matters: **the cure lives in the restart sizing,
not in the preconditioner the exclusion was waiting for.** (#282/#280 is the more
likely decisive one, since its docstring names the exact symptom.)

Timeline: the exclusions are stamped **`R-1 Step D (2026-05-19)`**; the first
restart fix landed **nine days** later and the second on 2026-07-04. The
exclusions survived **eleven weeks** past their own cure — and **nobody found
out, because an imperative `pytest.xfail` can never XPASS.** That is the
campaign's own thesis, measured on the campaign's own bookkeeping.

Corroborating detail — the reason strings are present-tense-false in *two*
independent ways:

* `test_kinf_homogeneous.py`'s xfail reasons cite the `max_inner=300` budget while
  the module's own `_TIGHT_KW` has read `max_inner=1000` since **2026-05-26**
  (`:129`, with the bump documented in place).
* `max_inner` is **not in iterations** for the krylov path. It becomes scipy's
  `maxiter`, which counts **restart CYCLES**, and `restart = n_dof`
  (`[M]` `n_dof = 112` for sphere-2eg) — so one cycle spans the full Krylov space
  and *no* value of `max_inner ≥ 1` can truncate this fixture. Measured: the same
  sphere-4eg krylov solve at `max_inner=2` returns `rel = 3.582e-15`, bit-for-bit
  the `max_inner=1000` answer. The named budget is not a live knob on this path.
  (This is the units side-finding already recorded for #340 N6; here it is a
  second, independent confirmation.)

### The healing is real, not a silent truncation

The campaign's own failure mode would be a row that passes because a truncated
inner solve was never reported. Two independent checks say it is not:

1. **Value.** All four rows agree with the **structurally-independent** closed-form
   reference (the dominant eigenpair of `A⁻¹F` by pure linear algebra, no transport
   discretisation) to `3.6e-15` / `4.5e-13` / `1.2e-14` — two to five orders inside
   each gate's own tolerance. A truncated inner solve cannot land there.
2. **Warning surface.** All four rows pass under
   `-W error::orpheus.numerics.convergence.ConvergenceWarning` (the DOTTED
   spelling — the bare name does not resolve, per vv Mode 8 EIGHTH class).
   `[M]` verbatim: `4 passed, 50 deselected, 1 warning in 15.92s` (the 1 warning is
   pytest's own `-O` `PytestConfigWarning`).

   ⭐ **With its positive control** (vv anti-#17/#19 — a positive reading alone
   carries zero information about whether the instrument can bite). The krylov path
   cannot be truncated via `max_inner` at all (above), so the control runs on the SI
   path, where `max_inner` IS iterations:

   | leg | config | `rel` | warnings |
   |---|---|---|---|
   | **positive control** | SI, `max_inner=2`, `inner_tol=1e-12` | `1.139e-02` | `['ConvergenceWarning']` |
   | negative leg | SI, `max_inner=1000`, `inner_tol=1e-12` | `6.492e-14` | `[]` |

   The instrument fires on a genuinely truncated inner solve and moves keff by
   `1.1e-2`, so the four silent krylov rows' silence is informative.

Probes: `scratch/_r5_probe_margins.py`, `scratch/_r5_probe_control.py`.

---

## 2. Rows 1, 3, 4, 5 — un-excluded (the exclusion has healed)

Per the brief: a row that now passes is **un-excluded**, not converted to a strict
xfail that would immediately red. The branches are deleted; the rationale is
preserved as **past-tense history** (never present-tense falsehood) in ONE place
per module rather than duplicated at each former site.

Row 1 carries an extra finding: its own comment claimed *"same exclusion as
test_kinf_homogeneous"*, but that module excluded only sphere-**4eg**-krylov,
while row 1's case is **2eg** and its condition omits `ng_key` entirely. So the row
excluded a combination its own cited authority never excluded — and
`test_kinf_homogeneous[sphere-2eg-krylov]` has been running and passing the whole
time. The exclusion was over-broad from the day it was written.

## 3. Row 2 — converted to a declarative `strict=True` (a LIVE todo)

`test_removal_form_kinf_independent_reference_2g` is an **unconditional stub with
an empty body** — the `pytest.xfail(...)` call *is* the body. So the brief's
"unconditional stub → `@pytest.mark.xfail(strict=True)`" cannot be applied
literally: an empty body under a strict xfail **XPASSes → FAILS** immediately.

Its stated blocker is **real**, verified against the tree (not the reason string):
`orpheus/sn/solver.py` exposes `solve_sn`, `solve_sn_adjoint`,
`solve_sn_adjoint_fixed_source`, `solve_sn_fixed_source` and none of them takes a
removal-form / σ_r-fold parameter; `grep -rn "removal_form" orpheus/` → 0 hits.
(The `ScatteringOperator.foldable_part` / `foldable_sigma` machinery exists at the
operator tier; no SOLVER entry consumes it.) The citation is legitimate: #200's own
body stacks *"the σ_r foldable preconditioner (Adams & Larsen 2002 §III) …
`σ_r = σ_t − Σ_(s,0)^(g→g)`"* after the block-inverse.

So the body becomes the smallest statement that (a) fails **for the documented
reason** and (b) **flips when the capability lands** — a concept-level capability
probe, not a guessed symbol: any public `orpheus.sn.solver` `solve_*` name or
`solve_sn` parameter whose spelling mentions `removal` or `fold`. Exactly one
statement in the body can fail, and it is the documented one (vv Mode 8 FOURTH
class); when the entry appears, the body falls through and PASSES, so the row turns
`XPASS(strict)` → **FAILED**, forcing the stub to be written for real.

⚠ Residual fragility, stated in the helper's docstring: if #200 lands under a
spelling containing neither `removal` nor `fold`, the probe does not see it and the
row will not flip. The mitigation is the grep pointer in the docstring.

### Both mandatory checks on the converted row

| check | why | result |
|---|---|---|
| `--runxfail` — does it red for **its documented reason**? (vv Mode 8 FOURTH class: an xfail hides *any* failure, incl. an incidental setup error) | otherwise the row is a misattributed xfail asserting nothing about #200 | ✅ `E Failed: no removal-form solver entry exists yet (#200): no public \`orpheus.sn.solver.solve_*\` name and no \`solve_sn\` parameter selects the σ_r-folded within-group form…` at `:542` — the documented reason, nothing else |
| **flip control** — simulate #200 landing (`scratch/_r5_flip_control.py` injects `orpheus.sn.solver.solve_sn_removal_form` **in-process**, never editing `orpheus/`) | proves the row MEASURES what #200 will change, i.e. that the marker really is a self-retiring todo | ✅ `FAILED … [XPASS(strict)] removal-form k_inf recovery needs the #200 solver entry…` — `1 failed` |

---

## 4. Acceptance

```
$ grep -rn '^\s*pytest\.xfail(' tests/
exit=1 (no matches)
```

The three prose mentions in `tests/numerics/test_advertised_degree_is_measured.py`
(`:83`, `:160`, `:208`) are untouched — `git diff --stat` on that file is empty.

`git status --porcelain orpheus/` → **empty**. No `orpheus/` file modified.

Verbatim summaries (canonical `python -O -m pytest`, SERIAL):

```
# tests/sn/verification/analytical/test_kinf_homogeneous.py
# tests/sn/verification/analytical/test_kinf_homogeneous_tolerance.py
39 passed, 1 warning in 51.95s
      (was 36 passed + 3 xfailed)

# tests/sn/operators/test_apply_full_field_codomain.py
# tests/sn/operators/test_removal_form_matvec_sweep.py
XFAIL tests/sn/operators/test_removal_form_matvec_sweep.py::test_removal_form_kinf_independent_reference_2g
  - removal-form k_inf recovery needs the #200 solver entry (fold within-group
    self-scatter into the diagonal). The Step-B operator round-trip + M(σ_r) value
    gate are the landable grounds; this eigenvalue cross-check lands with #200.
35 passed, 1 xfailed, 1 warning in 3.96s
      (was 34 passed + 1 xfailed  — the 1 xfailed is now DECLARATIVE and strict)
```

The 1 warning in each is pytest's own `PytestConfigWarning` about `-O`.

**pyright**: `29 errors` over the four modules — and **29 errors** over the HEAD
versions of the same four files placed as siblings (`git show HEAD:<f> >
<dir>/_r5head_<f>`, so imports resolve identically; never `git checkout`). All
pre-existing (`keff: float | None` arithmetic in the tolerance module); my change
adds **zero**.

**Sphinx**: `build succeeded.` under `-W --keep-going`.

**V&V harness**: `python -O -m tests._harness.audit` runs clean; ERR coverage
78/78; no marker/registry change (no `verifies()` edge was added or removed).

---

## 5. Scope note + one thing that refuted the framing

**Scope widened by two `docs/` files, deliberately.** R5 was briefed as "entirely a
`tests/` change". It could not be: the retirement audit's third search (grep the
CONCEPT across code, tests AND docs) found a doc claim that names the very file and
cell I changed and was rendered present-tense-FALSE by the change —
`docs/theory/methods/sn/curvilinear_multigroup.rst` §"Gotcha", which called
sphere-4G *"the xfail'd cell of the coordinates × groups × drivers grid
(tests/sn/verification/analytical/test_kinf_homogeneous.py)"*. Leaving it is the
exact state vv anti-pattern #21 forbids (the stale claim and its correction
coexisting, each citable). Rewritten as past-tense history with the measurement and
the correct attribution; the same page's earlier bullet list (`:85`) and the
admonition's "Two real … interactions exist" intro were corrected with it.
`docs/theory/methods/sn/loss_representation.rst:574` said row 2 "is ``xfail`` until
the #200 solver entry exists" — still true, sharpened to name the strict form.
`orpheus/` untouched throughout.

**⛔ What refuted the brief's framing.**

1. **"Unconditional stub → `@pytest.mark.xfail(strict=True)` on the test" cannot be
   applied literally.** Row 2's body is *empty* — the `pytest.xfail(...)` call IS
   the body — so a strict marker over it XPASSes and reds immediately. The
   conversion had to supply a body that fails for the documented reason; that body
   is the substance of the fix, not boilerplate.
2. **The exclusion's cited fix is not necessarily the fix that heals it.** All four
   healed rows named #200 as their re-enabler. #200 is still open. Reading the
   reason string and checking "has #200 landed?" would have concluded "still
   blocked" — the honest answer required *running the row*.
3. **The reason strings were false in a second, independent way nobody had
   noticed**: `max_inner` is not an iteration count on the krylov path (it is
   scipy's `maxiter`, counting restart CYCLES, with `restart == n_dof`), so the
   *named budget was never a live knob on the path being excluded*. `[M]`
   `max_inner=2` returns the `max_inner=1000` answer to the last bit. This is a
   second, independent confirmation of the units side-finding already recorded for
   #340 N6.
4. **A pre-existing generated-artifact staleness surfaced, not attributable to
   R5.** `docs/theory/verification/matrix.rst` regenerates to `9492` where HEAD's
   committed copy says `9491`; commit `3f76d651` added 59 lines of tests to
   `tests/numerics/test_iteration.py` without regenerating the matrix. My four
   modules collect **75** tests before and after the change. Flagging for the
   coordinator — the regenerated file is correct, the committed one was stale.

---

## 6. Instruments (all under `scratch/`, none tracked, all `_r5_*`-prefixed)

| file | role |
|---|---|
| `_r5_noop_xfail.py` | `-p` plugin neutralising imperative `pytest.xfail`; asserts its own installation via a `sessionfinish` count |
| `_r5_probe_margins.py` | the measured margins per row + the (failed, informative) attempt at a krylov-path truncation control |
| `_r5_probe_control.py` | the SI-path positive control proving `ConvergenceWarning` fires on a truncated inner solve |
| `_r5_flip_control.py` | `-p` plugin simulating #200 landing, proving the strict xfail flips to `XPASS(strict)` |
