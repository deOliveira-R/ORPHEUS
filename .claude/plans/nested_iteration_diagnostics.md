# A stalled solve says WHERE it stalled, WHY, and WHAT TO SET

**Status.** Design, awaiting steering. Opened 2026-08-09 out of #340/#341.
**Authority.** #340 + #341 latest comments (Cardinal Rule 4). Live index = task #57.

---

## The goal, in the domain's terms

> A user whose solve did not converge can read, from the returned object alone,
> **which level** of the iteration failed, **which criterion** it failed on,
> **how far** it got, and **what budget would have sufficed** — without
> re-running anything.

Not: "add a field to `IterationHistory`". That is a mechanism, and the wrong
one — see §2.

**Done when** (checkable): for the measured d=3 all-reflective eigenvalue solve
at `max_inner=200`, the returned object reports *"inner(within-group) hit
budget=200 at rho=0.985, needs ~849 for tol=1e-8"* — and a caller asking the
one question "is this answer trustworthy?" gets `False`.

**User directive of record** (2026-08-09, verbatim): *"We need rich diagnostics
for inner and outer iterations. Convergence problems and things like this are
some of the most common problems to have and some of the hardest to diagnose.
This needs to be designed properly."*

---

## 1. The measured facts this design must answer to

All `[M]`, all with configuration, all from this campaign.

| # | fact | evidence |
|---|---|---|
| F1 | A truncated inner corrupts `keff` by **5.3x `keff_tol`** while reporting `converged=True` and emitting nothing. | `scratch/probe_eigenvalue_inner_truncation.py`; 2-D (2.0,1.5)/(4,3), LS-S4, mixture A 2g, `keff_tol=1e-7`, `inner_tol=1e-10` |
| F2 | The outer's stop test is **entirely increments** — `dk < keff_tol and dphi < flux_tol` (`SNSolver.converged`, ~`solver.py:1556`; twin at `KEigenvalue.converged`, `numerics/iteration.py:1340`). A throttled inner suppresses the very quantities the outer reads. It does not converge; it **stalls and reads its own stall as convergence**. ⚠ `dk`/`dphi` are computed inside those predicates and **never stored** — §3a's `criteria` needs NEW recording at both, not plumbing. | source |
| F3 | `n(tol) = \|ln(tol/r0)\| / \|ln rho\|` holds with **R2 >= 0.9956** over 8 configs; fitted slope = `1/\|ln rho\|` to 0.6 % on the 6 well-conditioned rows. **The budget is derivable; `rho` is the only unknown.** | `scratch/probe_inner_budget_law.py` |
| F4 | Both shipped defaults are **SHORT** at d=3 zero-leakage — incl. with scattering: need@1e-8 = 830 (c=0.5) vs `have 200`; need@1e-12 = 1441 vs `have 1000`. | same |
| F5 | The tolerance factor between the two default families should be **~1.5-1.75x** by F3; the tree ships **5x**. | same |
| F6 | CP / MoC / Diffusion receive `outcome.converged` and **drop it**; none of `CPResult` / `MoCResult` / `DiffusionResult` has the field. | `cp/solver.py:884`, `moc/solver.py:139`, `diffusion/solver.py:409` |
| F7 | `_certify_within_group_exit` is wired at 4 sites, **all within-group**. The outer has no certificate. The antidote concept is in the tree, applied one level down only. | `solver.py:1744, 1928, 3440, 3654` |
| F7′ | ⛔ **F7 UNDERCOUNTS THE GAP** — corrected 2026-08-09 by survey. There are **6** within-group inners, not 4: `solve_sn_adjoint`'s (`KEigenvalue.solve_fixed_source`) and `solve_sn_adjoint_fixed_source`'s `si.solve` are **UNCERTIFIED**. The honest statement is *"4 of 6 within-group inners certified; the outer has none"*. | `numerics/iteration.py:1259`; `solver.py:2778` |
| F8 | The per-inner status is **available and discarded**: `solver.py:1754` / `:1937` do `self._total_inner_iterations += len(_residuals)` — summing the COUNT, dropping the STATUS. | source |
| F9 | ⭐ `SourceIteration` **already computes the rate** F3 needs: `contraction_ratios` + `last_displacement` (`numerics/iteration.py:693`, filled `:754`). It dies with the local `si` — **zero production consumers**. §3a's `rate` is a *retention* problem, not a computation one. | survey |
| F10 | ⛔ **`len(residual_history)` is NOT injective at the budget boundary.** `SourceIteration.solve` appends from the SECOND iteration on, so `max_iter=50` yields `len == 49` for *both* "exhausted" and "converged on the last permitted check". `[M]` a 30-outer probe showed all 30 inners at 49. | survey |
| F11 | `n_inner` already means **two different things**: SI writes `len(residuals)` (`solver.py:3468`), Krylov writes `len(residuals) + 1` (`:3685`) — an undocumented `+1`, via a local misnamed `n_outer` holding an inner count. One test compares them directly across drivers. | `test_si_convergence_rate.py:240` |
| F12 | `_total_inner_iterations` is an **instance accumulator with no reset** (`solver.py:1092`), so reusing one `SNSolver` across two `power_iteration` calls silently double-counts. `solve_sn_fixed_source` builds an `SNSolver` that never reads it. | survey |

### 1b. The level tree, MEASURED — and it is a per-ENTRY-POINT fact

`[M]` 2026-08-09 by `sys.setprofile` on a real 2-region 1-D slab solve
(mixtures A/B 1g, GL-S8, 16 cells). This **corrects the shape §2 and §3a were
written against** — they assumed "outer over per-group inners over Krylov".

⛔ **There is no per-group loop in SN.** `build_within_group_system` /
`_within_group_si` build ONE `SourceIteration` over the full
`(N, ng, nx, ny)` state with the full multigroup `S`. "Within-group" here
means *fission-external*, not *one group at a time*. **The SN tree is 2 deep,
never 3** — the Krylov restart cycle is a third level *inside scipy*, with no
ORPHEUS frame.

| entry | depth | shape |
|---|---|---|
| `solve_sn` (SI or Krylov) | 2 | `power_iteration` → `SNSolver.solve_fixed_source` → `_solve_source_iteration` / `_solve_krylov` |
| `solve_sn_adjoint` | 2, **one frame deeper** | `KEigenvalue.solve` → `power_iteration` → `KEigenvalue.solve_fixed_source` |
| CP | 2, **per-group** | outer → `for g: for n_in:` — the only family already reporting a per-outer × per-group inner *count* (`CPResult.n_inner`), still status-free |
| MoC | 2, **no criterion** | `for _inner in range(n_inner_sweeps)` — a fixed count, no tolerance, no residual |
| Diffusion | **1** | one LU back-substitution; no inner level exists |

⟹ the record must be **genuinely recursive and genuinely ragged**: diffusion
is a 1-level tree, MoC's inner has no criterion at all (`criteria=()`), and CP
fans out per group. A 3-level assumption would have been built into the type.

⛔ **The static graph cannot see this tree.** `callers` on `SNSolver.converged`
returns **0** — it is reached polymorphically through the `EigenvalueSolver`
Protocol inside `power_iteration`. Design N2 against the traced tree, not a
`callees` walk.

⛔ **Largest structural obstacle to N2**: `KEigenvalue.solve_fixed_source`
(`numerics/iteration.py:1270`) does `psi, _inner_residuals = self._inner.solve(...)`
and returns only `psi`. The adjoint eigenvalue path's inner trajectory is
discarded **inside a shared numerics primitive**, so it cannot be recovered
from `orpheus/sn/solver.py` at all. `[M]` `solve_sn_adjoint` returns
`total_inner_iterations is None` where the forward path returns `1470`.

⚠ And `IterationHistory` lives at **`orpheus/sn/solution.py:112`**, not in
`numerics/` — it is SN-local, which is itself the F6 problem.

⛔ **REFUTED 2026-08-09** — an earlier #340 comment (mine) put the `derivations/`
blast radius at "~87 reads". That was a name-grep over
`PeierlsGreensFunction*Result.converged`, a **different and already-honest**
family. Measured: **2 reads, 0 constructions lacking `converged=`**. The text
stays so the next session does not re-inherit the inflated number.

---

## 2. Why a field on `IterationHistory` is the wrong mechanism

The fact being recorded is a **tree** — outer over inners over (sometimes)
Krylov — and `IterationHistory` is a **flat record**. Every symptom in §1 is
that shape mismatch surfacing:

- `converged: bool` has no room for *which level*, so it silently means "the
  outermost one" (F1).
- `total_inner_iterations` is the flattening: it keeps the **sum of counts** and
  cannot express *"3 of 47 outers had a truncated inner"* (F8).
- Five transcriptions of one predicate existed because the fact was **stored**
  at each site rather than **derived** once from the trajectory.

Adding `n_truncated_inners` treats the symptom and leaves the projection lossy —
the exact move [[feedback-lossy-return-type-is-the-root-cause]] warns against,
one level up from the `power_iteration` triple it already caught.

## 3. Proposed means (2026-08-09 — DESIGN, not verified)

### 3a. The primitive: a recursive iteration record in `numerics`

Home: `orpheus/numerics/convergence.py`, beside `ConvergenceWarning`. Model
independent by construction, so CP / MoC / diffusion / SN share it — which is
what the "structural, not SN-local" ruling asked for and F6 shows was not
delivered.

```python
@dataclass(frozen=True)
class IterationRecord:
    """One LEVEL of an iterative solve: what it wanted, what it got, why it stopped."""
    label: str                                   # "outer(power)", "inner(within-group)"
    criteria: Mapping[str, tuple[float, ...]]    # name -> trajectory
    tolerances: Mapping[str, float]              # name -> the tol it was judged against
    budget: int
    min_iterations: int = 0                      # the outer's `iteration <= 2` guard
    children: tuple["IterationRecord", ...] = ()
```

⭐ **`converged` is DERIVED, never stored.** That is the whole point: a fact you
cannot transcribe cannot be transcribed wrongly. The five spellings #342 found —
and the sixth the `derivations/` audit found (`... and iteration > 5`) — become
*unspellable*, not merely fixed. (`coding-elegance` Pattern 4.)

Derived surface:

| property | answers |
|---|---|
| `converged` | did THIS level meet all its criteria |
| `fully_converged` | ...and every descendant — **the honest answer to "trust this?"** |
| `binding_criterion` | WHICH criterion was last to clear (the `dk`-vs-`dphi` question) |
| `rate` | observed rho from the trajectory tail |
| `projected_iterations` | budget that WOULD have sufficed, at that rate (F3) |
| `first_failure` | the deepest/earliest record that did not converge |
| `report()` | the tree, human-readable — the thing you paste into an issue |

Multi-criterion is not gold-plating: F2 is *precisely* a two-criterion stop
(`dk` AND `dphi`) whose components are invisible today, and "which one is
lagging" is a question users actually ask.

### 3b. `converged`'s meaning — the question that had no good option

The three options offered (conjunction / sibling field / certificate) were all
patches on a flat type. Under 3a the question dissolves: `converged` keeps its
honest per-level meaning **and** `fully_converged` answers the caller's actual
question. No semantics of an existing field changes; nothing is hidden.

⚠ ~~**This is the design's load-bearing claim and it is UNVERIFIED.** It
predicts zero test churn from the meaning question. The running truncated-exit
audit (#340 bullet 2) measures exactly the set that would churn. **Do not
implement 3b before that lands.**~~

✅ **UNBLOCKED 2026-08-09 — the audit landed** (`adc887d6`; the sweep found
**7 deliberate truncations, ZERO false greens**, all passing explicit budgets).
The churn prediction is therefore testable and the blocker is void. Two
consequences the audit itself supplies: a changed *default* budget cannot
disturb those 7 rows (they all pass explicit ones), and no gate anywhere
depends on the current default values — which is what makes 3c safe to land.

### 3c. The budget: derived, and the warning says what to set

RULED 2026-08-09 by the user: *derived default + rate-aware warning.*

```python
_SERVED_RATE = 0.986        # [M] the slowest rho we promise to serve
def _default_inner_budget(tol: float) -> int:
    return ceil(log(tol) / log(_SERVED_RATE))
```

Retires all ~~five~~ **six** constants (F5's unprincipled 5x with them).
Signature becomes `max_inner: int | None = None` — `None` = derive, an int =
deliberate override.

⛔ **"five" UNDERCOUNTS — corrected 2026-08-09 by survey.** There is a SIXTH,
and it is not in SN: `KEigenvalue.__init__` (`numerics/iteration.py:1161`)
defaults `max_inner=1000, inner_tol=1e-8` — a **third** (budget, tol)
combination that no SN entry spells. `solve_sn_adjoint` passes its own
200/1e-8 through, so the pair is only reachable by direct `KEigenvalue`
construction, but it is a live default the derivation must cover. Outside SN:
CP carries `max_inner=100 / inner_tol=1e-8` **twice** (`cp/solver.py:68` and
`:463`), and MoC carries `n_inner_sweeps=15` twice (`moc/solver.py:78`,
`moc/core.py:37`) — each pair a twin source of truth in its own right.

The warning then reads `first_failure` and reports the **observed** rate and
projection, turning *"raise the budget"* into *"set `max_inner=1710`"*. This is
what makes the default's exact value stop being load-bearing — and F3 is what
makes the projection honest rather than a guess.

✅ **N3 LANDED 2026-08-09** — and the acceptance evidence is stronger than
the plan asked for. `[M]` real converged solves of the d=3 all-reflective
absorber (the configuration whose truncated exit opened #340), `NEEDED` read
off `len(history.flux_residuals)` on a run with `converged is True`:

| `inner_tol` | NEEDED | shipped 1000 | derived |
|---|---|---|---|
| 1e-9 | **1007** | ✗ | 1471 ✓ |
| 1e-12 | 1473 | ✗ | 1961 ✓ |
| 1e-13 | 1631 | ✗ | 2125 ✓ |
| 1e-15 | 2031 | ✗ | 2451 ✓ |

⭐ **The derived budget covers every row; the shipped constant covers none** —
and the first miss is by **seven sweeps**, which is the margin the founding
defect's gate rode for months. The law's own ratio checks out too: derived
`1e-12 / 1e-8 = 1.499` against the theoretical `1.500`, where the tree
shipped `5.0`.

`[M]` the rate-aware message validates against the same solves — it projects
**1618** from a 200-iteration tail against a true **1631** (0.8 %), and 1633
from a long one. It reads LOW only from deep inside the transient (586 at
budget 50) and converges **from below**, so following its advice yields a
bigger number next time rather than a wrong answer. The message says "at the
rate observed so far" for exactly that reason.

⚠ `_SERVED_RATE = 0.986` covers the representative d=3 worst (`rho = 0.9854`)
but **NOT** `Sigma_t/4` (`rho = 0.99575`, needs 5171 @1e-12). That is a
deliberate, stated limit, not an oversight: serving 0.996 means a *stuck* d=1
slab churns ~6900 sweeps before admitting defeat. The rate-aware message is what
makes the uncovered tail cheap — it tells that user the number.

### 3d. Sizing / order

| step | scope | gate |
|---|---|---|
| N1 ✅ | `IterationRecord` + `StoppingCriterion` + their own laws | **LANDED** — see below |
| N2a ✅ | **Producer-side retention** — the two drivers return a record | **LANDED** — see below |
| N2b-i | **The OUTER stopping test reports what it measured** — see below; this row's original wording was too narrow | eigenvalue suites + new loop gates |
| N2b-ii | SN's SOLUTION carries the record; `IterationHistory` composes it; flat accessors become derived | the 9-red baseline unchanged |
| N3 | Derived budget + rate-aware message (3c) | contract suite + a NEW starved-projection gate |
| N4 | CP / MoC / diffusion carry it (F6) | their own smokes |
| N5 | Outer certificate (F7/F7′) — residual `(A - F/k)psi` at exit | mutation battery |

**N1 landed 2026-08-09.** `orpheus/numerics/convergence.py` +
`tests/numerics/test_iteration_record.py`, 68 gates. `[M]` mutation battery
**12 of 12 caught, zero blind**, positive control (constant `rate`) reds 12.
Two blind spots were found and closed by measurement rather than by review —
recorded because both are reusable test-design lessons:

- **MA3** (`first_failure` checks self before children) reddened **0**: every
  fixture's only failing level was already the deepest, so both traversal
  orders returned the same record and the ordering claim was unfalsifiable.
  Closed by a parent-AND-child-both-failing fixture — the common real shape.
- **MA6** (read `r0` off the trajectory instead of the fitted intercept)
  reddened **0**: on a *pure* geometric fixture the naive law and the fitted
  crossing coincide **by construction**. Closed with a two-mode decay
  `A rho^k + B sigma^k` (`B/A = 50`, `sigma = 0.5 << rho = 0.98`), where the
  two readings differ by **195 of 685 iterations = 28.5 %** — the same
  magnitude as the campaign's measured 27 % intercept share. The reference is
  a direct walk of the true sequence, independent of both.

⭐ **N2 splits, because retention must precede production.** §3a's `criteria`
cannot be assembled from what the tree currently keeps: `dk`/`dphi` are
computed inside the stop predicates and discarded (F2), `contraction_ratios`
dies with the local `si` (F9), `KrylovAcceleration.solve` never returns
scipy's `info`, and `KEigenvalue.solve_fixed_source` drops its inner history
inside `numerics/` (§1b). N2a is that retention pass; only then can N2b
compose a record.

**N2a landed 2026-08-09.** `SourceIteration.solve` and
`KrylovAcceleration.solve` return `(state, IterationRecord)` instead of
`(state, list[float])`; 7 production unpack sites migrated;
`_claims_convergence` **retired** (zero production callers) with both prose
survivors repointed by tense; `_history_from_record` is the ONE place the flat
carrier is derived from the recursive one. `[M]` mutation battery **15/15
caught, zero blind** (74 gates; positive control reds 12).

Four things it fixed that the plan had only recorded as facts:

- **F10/F11 dissolve into one rule** — each driver states its own count.
  `iterations_run` is the one thing the record cannot derive, because the
  offset is per-producer. ⭐ And the tree's two conventions were **backwards**:
  SI wrote `len(residuals)` where its pass count is `len + 1`, Krylov wrote
  `len + 1` where its counts match.
- **F9** — the rate is now retained: it derives from the record's trajectory.
- **The adjoint eigenvalue inner is no longer lost inside `numerics/`**
  (`KEigenvalue` accumulates `inner_records`, reset per `solve`, so F12's
  double-count cannot recur). `total_inner_iterations` stops being `None`.
- **The `cleared`-guard branch in `_warn_if_unconverged` is still live** —
  it retires when N2b records `dphi`.

⛔ **Two corrections the tests forced, both worth carrying:**

1. **`StoppingCriterion` first shipped demanding a strictly positive
   tolerance, and TWO production paths could then not build a record**
   (`GreenOperator(tol=0)`, GMRES exact-breakdown). A type asserting more
   than its callers promise — vv anti-pattern #16 in the direction that
   breaks production. `0.0` is now legal and never clears; `distance`
   returns `inf` there so `cleared ⟺ distance < 1` still holds.
2. **A level that RAN and measured nothing must not claim convergence.**
   `cleared` is vacuously true on an empty trajectory — right when the loop
   never entered (GMRES on its initial guess), wrong when it entered and
   never measured (SI at `max_iter=1`). The retired predicate was wrong the
   *other* way. `iterated` is the discriminator; one rule serves both.

`[M]` **the only test-visible consequence was a +1 relabel**, re-baselined
with its rationale recorded at the claim: `test_dsa_rate.py`'s S2-exactness
counts went 2/3/3 → 3/4/4 and 2/3 → 3/4. The physics is untouched (machine-
zero landing unchanged at 3.20e-15 / 6.20e-15) and the DELTA the class
actually claims — one extra pass per lagged reflective wall — reads
identically. Gates: `tests/sn/{solve,architecture,acceleration,operators,
eigenvalue}` + rate + solution = **3 failed / 1571 passed**, the 3 being the
known #333 affine trio; `tests/numerics` 2204 passed.

⛔ **N2a owes a producer-side normalisation, not a consumer-side workaround**
(`coding-elegance` Pattern 7). Per F10, a trajectory that starts at the
*second* iteration makes `n_iterations >= budget` under-report by exactly one,
so `exhausted_budget` would read False on a genuinely exhausted level. The fix
belongs at `SourceIteration.solve` (one entry per iteration, initial residual
included), NOT at the record — deriving a count from a short list is the
guessing F10 measures. Same for F11's undocumented `+1`: two spellings of
"iterations at this level" must become one before a record can carry it.

N3 is independent of N1's consumers and can proceed on the landed primitive.

### N2b-i — the outer stopping test reports what it MEASURED, not a verdict about it

**Goal (domain terms).** A stopping criterion is a measured magnitude judged
against a tolerance. A loop that keeps only the boolean cannot say which
criterion failed, how fast it was closing, or what budget would reach it —
so the object a user gets back cannot answer the question this campaign
exists to answer.

⛔ **The §3d row above said "SN produces the record", and that framing was too
narrow — corrected 2026-08-09 in place per `plan-authoring` §3.** The producer
is not SN. `EigenvalueSolver.converged(...) -> bool` is a `numerics/` protocol
with FIVE realizers, and every one of them computes `|Δk|` and `‖Δφ‖`, compares
both, and returns one bit. It is the *identical* lossy-return-type defect N2a
fixed at the inner level, one level up — and the campaign's standing ruling
(**structural, not SN-local**) already authorises fixing it where it lives.

`[M]` blast radius at `bfd59dd9`, and it is why this is atomic:

| | count | where |
|---|---|---|
| realizers of `converged` | **5** | `numerics/iteration.py:1435` · `sn/solver.py:1661` · `cp/solver.py:722` · `moc/core.py:242` · `diffusion/solver.py:316` |
| production call sites | **1** | `numerics/eigenvalue.py:317` |
| transcriptions of `if iteration <= 2: return False` | **5** | one per realizer, character-identical |

Per `plan-authoring` §6b the unit of work is the call-site set, so 1–3 below
cannot be split; the solution-object half genuinely can, and is **N2b-ii**.

**Means (LANDED — see the block below for what measurement changed).**

1. `converged(...) -> bool` becomes
   `measure_stopping_criteria(...) -> tuple[StoppingCriterion, ...]`, one
   single-point reading per criterion via a new `StoppingCriterion.reading()`
   named constructor; the loop concatenates with `extended_with()`, which
   routes through `replace()` so `__post_init__` re-fires (Pattern 4 ∩ 2).
   **No new type**: a reading IS a criterion with a length-1 trajectory, and
   minting `CriterionReading` would be the same object wearing a label
   (`coding-standards` type-vs-property).
2. The `iteration` **parameter is dropped**. A reading is a function of the two
   iterates; keeping the parameter would leave the sixth transcription of the
   guard spellable. The rule is a property of power iteration (two increments
   are the minimum that can show a trend) and now lives once as
   `MINIMUM_OUTER_ITERATIONS`, reaching the record as `min_iterations` — which
   N1's own docstring had already named as its intended home.
3. `PowerIterationOutcome` gains `record` and `converged` becomes a **derived
   property**, so the outer level obeys the same rule as the inner one.

⭐ **Two things the implementation found that the plan had not predicted:**

- **The FORWARD path had no subtree.** `SNSolver` was not a recording solver:
  `:1877`/`:2062` did `self._total_inner_iterations += record.n_iterations` —
  the record in hand, reduced to a scalar and dropped, one line after the
  driver measured it (F8). So the adjoint path could carry a tree and the
  forward path could not. Retiring that counter onto `inner_records` (the total
  becomes `sum(r.n_iterations for r in ...)`, which also kills the second
  spelling of it) is what makes the headline work at all.
- **F12 is fixed STRUCTURALLY, not by a reset.** The double-count defect was
  "instance accumulator with no reset". Rather than ask five realizers each to
  find a start-of-solve hook, `power_iteration` records `len(inner_records)` at
  entry and slices from there — it takes the children *it* caused. A realizer's
  reset hygiene can no longer corrupt the tree. `KEigenvalue`'s reset stays,
  now only keeping its own public attribute honest.

`[M]` **the headline, end to end, on a 20-cell 2-group S8 slab** — the first
time the campaign's motivating question has an answer in the returned object:

```
max_inner=1   keff = 1.59689680      (honest: 1.59689789 — 11x keff_tol off)
  outer.converged   = True                     <- the outer's OWN test is satisfied
  fully_converged   = False                    <- the tree knows better
  first_failure     = ('inner(source-iteration)', 'TRUNCATED')
```

Note the outer took **154** outers rather than 17: the starved inner throttled
it, and the increment test eventually satisfied itself on a wrong answer —
exactly F5's "an increment-only stopping test cannot see an upstream throttle".

⚠ **Behaviour-neutrality of the derived `converged` is ARGUED, so it owes its
own gate** (`plan-authoring` §8). The argument: the loop breaks iff all of this
iteration's readings cleared, and each criterion's `last` IS that reading, so
`all(c.cleared)` agrees; `n >= 3` holds on both sides by the same constant.
Edge rows that must be gated, not assumed: `max_iter=0` (no criteria at all —
`DIRECT`), `max_iter <= 2` (readings may clear while `min_iterations` refuses),
and a solver returning **zero** criteria (refused outright — an empty
conjunction is vacuously true and would manufacture a claim out of nothing).

⚠ **CP's realizer is the un-pure one** and the rename made that louder: it
appends to `residual_history` and prints. The per-iteration diagnostic line
stays (it logs what was measured); the `"Converged."` **announcement was
deleted** — a method that only measures is not entitled to announce a verdict
it no longer reaches. Its iteration index now comes from
`len(self.residual_history) + 1`, exact because the append is once per call.
CP also judges `‖Δφ‖_∞` where the other four use relative `‖Δφ‖₂` — which is
the concrete reason the SOLVER must report and the loop must not compute.

**Also fixed in passing:** `IterationRecord.report()` printed `met` for a
criterion with an EMPTY trajectory, so a `TRUNCATED` level's own report
contradicted its status line. `cleared` is vacuously true there and that is
correct for what it decides; the report now says `not measured`.

#### N2b-i gates — `[M]` mutation battery 10/10, zero blind

MB0 (unmutated) is the positive control. `tests/numerics/test_power_iteration_record.py`
is new (the loop in isolation, on a SCRIPTED solver — a converging physics
solve couples accumulation, verdict and tree assembly to physics that can mask
an off-by-one in any of them).

⭐⭐ **Two of the ten reddened NOTHING on the first run, and both were real
gate gaps rather than harness faults — the distinction is the whole value of
the run.** Recording both, because each generalises:

1. **MB7 — `PowerIterationOutcome.converged` re-pointed at `fully_converged`
   reddened 0 of 249.** Every gate read `outcome.record.…` directly; the one
   that cross-checked the two surfaces used fixtures where the level and the
   fold agree. ⟹ **reaching THROUGH a wrapper in every gate leaves the
   wrapper's own surface — the one callers touch — ungated.** Only a
   starved-inner fixture separates them.
2. ⭐⭐ **MB9 — SN reporting only its FIRST criterion reddened 0 of 249**, and
   this is the sharper, more transferable one: **a mutation that RELAXES a
   conjunction is invisible to every gate that only checks the conjunction's
   RESULT.** Dropping `dphi` makes convergence strictly easier, so every
   outcome-level gate stays green while the solve is judged on half its
   contract — which is *precisely* the pre-#340 state being restored. A
   weakened stop test needs a gate on the REPORTED CRITERIA, not on the
   verdict. (Same family as the campaign's "a boolean-contract gate needs
   PAIRS": the starved leg carries the teeth.)

Both are now gated, in tests that say in their own docstrings that a mutation
found their absence.

⚠ **And two of the ten were BLIND FOR A HARNESS REASON, which reads identically
in the report and is the more dangerous confusion** (vv anti-#17, flattering
direction). MB3/MB4 patched `ev.power_iteration`, but the test module had done
`from ...eigenvalue import power_iteration` at import time, so its global still
pointed at the original and the mutation never reached the call site. A
"0 caught" that means *my instrument missed* looks exactly like a "0 caught"
that means *my gates are fine*. ⟹ when a mutation reddens nothing, the FIRST
question is "did it take?", answered by patching every already-bound name.

Suites at the N2b-i tree: `tests/numerics` + `tests/sn/{solve,eigenvalue}` +
`tests/diffusion` = **3 failed / 2506 passed** (the 3 = the known #333 affine
trio, unchanged); `tests/moc` + the new gates = **249 passed**; pyright
`orpheus/` **1 error** (the accepted #288 residual); xrefs **0 dead / 14 955**;
V&V matrix 9327 → **9361**.

⛔ N5 is NOT plumbing. Per `plan-authoring` §8: today nothing branches on an
outer residual, so adding one that RAISES changes behaviour at every eigenvalue
call site. It ships behind the same warn-not-raise ruling as #340, or not at all.

---

## 4. Refuted / rejected candidates — with the structural reason

| candidate | why it fails |
|---|---|
| `n_truncated_inners` field on the flat type | §2: keeps the lossy projection; cannot say WHICH level, WHICH criterion, or WHAT TO SET |
| Make `converged` the conjunction outright | Loses the per-level fact. "Did the outer converge?" becomes unanswerable — needed to tell "inner starved" from "genuinely non-critical" |
| Warn per truncated inner | Fires up to `max_outer` times. A warning that always fires is filtered, and the next truncation goes unnoticed — the MA4 mutation lesson, already gated by `test_converged_solve_is_SILENT` |
| Branch the budget on `ndim` | #341 falsified `ndim` as the discriminating variable **both ways** (d=2 G-S loses 2.86x; d=3 G-S wins 0.58x). Encoding a proxy for a variable we now understand is stringly-typed dispatch |
| Derive the budget by fitting the 3 dimension points | 3 points, no mechanism. F3 gives a LAW with R2 >= 0.9956 — fit the law, not the samples |
| Use `eq46_residual` as the `derivations/` convergence witness | No threshold exists anywhere for it ("Should be small"). Picking one is the "assert a number nobody measured" move this campaign removes. The bisection's own `bisect_tol` (`slab/one_group.py:211`) is the real witness |
