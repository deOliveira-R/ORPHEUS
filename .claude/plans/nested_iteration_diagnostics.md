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

⛔ **THE CP ROW IS TWICE WRONG** — corrected 2026-08-10 by survey against the tree
(`scratch/n4_iteration_record_surface.md`). Both corrections make CP *easier*, not
harder, so the row was mis-scoping N4 pessimistically:

1. **CP's inner is NOT "a count only".** It has a real tolerance and budget
   (`inner_tol=1e-8`, `max_inner=100`) and it *computes a residual* `res_in`,
   tests it, and throws it away — `cp/solver.py:595-600`. CP is the family
   **closest** to expressible: two lines of capture, no design decision.
2. **CP's level count is MODE-DEPENDENT, and the shallow mode is the DEFAULT.**
   2 levels under Gauss-Seidel, **1 under Jacobi** — and Jacobi is what ships
   (`_solve_fixed_source_jacobi`, `cp/solver.py:510-537`, has no inner loop at all
   and leaves `n_inner is None`). N4 must model both, so "CP is 2-deep" would have
   built a shape the default path cannot fill.

⚠ And CP inherits F10's disease independently: `for n_in in range(1, max_inner+1):
… if res < tol: break` makes `n_in == max_inner` mean "converged exactly at the
cap" **or** "exhausted" — not injective, same as `len(residual_history)`.
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
| N2b-ii | SN's SOLUTION carries the record; `IterationHistory` becomes a VIEW; flat accessors derived | the 9-red baseline unchanged |
| N3 | Derived budget + rate-aware message (3c) | contract suite + a NEW starved-projection gate |
| N4 | CP / MoC / diffusion carry it (F6) | their own smokes |
| N5 | Outer certificate (F7/F7′) — residual `(A - F/k)psi` at exit | mutation battery |
| **N6a** | **The warning speaks for the level that FAILED** — see §N6. Found 2026-08-10 by existence-checking this plan's own Done-when; it is NOT in the original decomposition because the original assumed the report was the hard part and the delivery was free | G1–G9 + a 12-mutation battery |
| **N6b** | The guard widens to `fully_converged` — **rides N5**, because the certificate is what separates a corrupting truncation from a benign one | the `xfail(strict=True)` marker XPASSes |

⭐ **N6 splits for the same reason N2 did, and it is worth naming the pattern:
retention/reporting must precede the decision that consumes it.** N2a had to
retain before N2b could compose; N6a has to make the message *correct for an
inner* before N6b can make it *fire for an inner*. Flipping the guard first
would have shipped a message that names the entry's `max_outer` for a starved
inner — advice pointing at a knob that cannot help — to 20 more call sites.

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

### N2b-ii — the solution carries the tree, not a flattened projection of it

**Goal.** A caller holding a returned solution can ask "did this converge, and
if not: which LEVEL, which criterion, how fast, what do I set" without reaching
past the object it was given.

**Means (user-ratified 2026-08-09: "iteration history as a view of iteration
record is the right design").** `IterationHistory` keeps ONE field —
`record` — plus `keff_history` (a physics output, not a stopping criterion).
`converged`, `flux_residuals`, `n_inner`, `n_outer`,
`total_inner_iterations` all become **derived properties**, and
`fully_converged` is exposed for the first time because the flat type had no
way to express it.

⭐ **The leaf-vs-outer discriminator is STRUCTURAL, not a label.** `n_inner` is
`None` on the eigenvalue path and `n_outer` is `None` on the fixed-source path,
and the view decides which by asking *does the top level have children* — a
level that drives nested solves is an outer; a leaf IS the solve. Reading it
off `record.label` would be stringly-typed dispatch on a string chosen for
humans.

⚠ **`flux_residuals` deliberately stays empty on the eigenvalue path** even
though `dphi` is now recorded there. Widening the name to "`dphi` when there is
no residual" would silently re-point every consumer that branches on
`if history.flux_residuals:`. The eigenvalue path's criteria are reachable —
via `history.record` — which is the whole point of the view.

**`_history_from_record` is RETIRED.** It existed to keep the lossy projection
in one place; with the view there is no projection, so the three fixed-source
entries construct `IterationHistory(record=record)` directly and the two
eigenvalue entries add `keff_history`. All five now read identically, where
before the two eigenvalue entries hand-wrote three facts in two different ways.

**The temporary honesty patch is RETIRED** (`_warn_if_unconverged`'s
"ALREADY cleared" branch). The warning now reads `record.binding_criterion` —
the criterion furthest from clearing, i.e. the one that actually failed — so
the situation the stop-gap guarded (*the failing criterion is absent from the
history*) is *unrepresentable*. Its test is re-posed rather than deleted, and
the claim strengthens from "say nothing rather than something wrong" to **"say
the right thing"**: given `dk` cleared and `dphi` stalled, the message must
name `dphi`. The replacement branch is a genuinely different claim — every
criterion cleared yet the level did not converge means the LOOP refused
(`min_iterations`), which is a state worth naming rather than a budget to
project.

⭐⭐ **A hole in the N1 primitive, found by writing the NEIGHBOURING test.**
`IterationRecord.converged` guarded "a level that RAN and measured nothing"
as `any(criterion.n_iterations == 0 …)` — which is silent on **no criteria at
all**, because `any(())` is `False` and the rule then falls through to
`all(())` and claims convergence. That is the identical vacuous-conjunction lie
`power_iteration` refuses at the producer, left open one layer down. `[M]`
found by writing the *pair* (converged-when-never-entered **and**
not-converged-when-iterated), not by review: the first test passed, the second
failed on correct-looking code. Widened to `not any(n_iterations > 0)`, which
agrees with the old spelling wherever a criterion exists (criteria are
co-indexed) and closes the zero-criteria case. ⟹ **vv #11's positive+negative
pairing applies to a DERIVED PROPERTY too, and the pair is what finds the
shape the rule forgot.**

#### N2b-ii gates — `[M]` mutation battery 7/7, zero blind (after closing 2)

MC0 positive control green. Two reddened nothing on the first run; both were
real gaps in the DERIVED READINGS, and one carries a transferable lesson:

⭐⭐ **MC3 — a `>=` gate is blind to any mutation that collapses one side onto
the other.** Replacing `total_inner_iterations`'s sum-over-children with the
outer's own count reddened **0**, and the near-miss is the point:
`test_si_convergence_rate` asserts `total_inner >= n_outer > 0`, which *reads*
like the guard for exactly this — and cannot be, because the mutation makes the
two sides equal and `n >= n` holds. An ordering assertion between two
quantities pins their ORDER, never their independence. The gate that works
needs a fixture where the two genuinely differ (children that each ran 7
passes, so 21 ≠ 3) plus an explicit `!=`.

**MC4 — the deliberate non-widening of `flux_residuals` was ungated**, as
predicted before the run. A design choice nothing gates is an ungated coverage
claim, and this one is a one-line "improvement" away from silently re-pointing
every `if history.flux_residuals:` consumer.

Suites at the N2b-ii tree: `tests/sn/primitives` + `tests/sn/solve` +
`tests/numerics` = **2280 passed**; `tests/sn/{solve,eigenvalue,acceleration,
architecture,operators}` + `tests/diffusion` + `tests/moc` = **3 failed / 1811
passed** (the 3 = the known #333 affine trio); pyright `orpheus/` **1 error**
(the accepted #288 residual); xrefs **0 dead / 14 962**; V&V matrix
9361 → **9366**.

---

### N6 — the record ANSWERS the done-when; nothing SAYS it

⭐ Found 2026-08-10 by existence-checking the plan's own **Done when** against the
tree (`plan-authoring` §1), not by planning the next step. The done-when was
written before N1–N3 landed, so the honest first question was *"does the tree
already say this?"* — and the answer is **half yes**, in a way that no step in
§3d's table covers.

**`[M]` `scratch/probe_340_done_when.py`, HEAD `52650a86`.** Configuration:
d=3 extents `(1.0, 2.0, 3.0)`, cells `(3, 4, 5)`, **all six faces reflective**,
LS-S4, 2 groups, fissile `c=0.5`, `max_inner=200`, `inner_tol=1e-8`, source
iteration on a `gauss_seidel` schedule.

| the done-when asks | the object answers | ✓ |
|---|---|---|
| WHICH level failed | `inner(source-iteration)` | ✅ |
| WHICH criterion bound it | `residual` | ✅ |
| HOW FAR it got | `200/200`, `exhausted_budget=True`, `status='TRUNCATED'` | ✅ |
| WHAT BUDGET would suffice | `rho=0.984928` ⟹ `projected_iterations()=838` | ✅ |
| is it TRUSTWORTHY | `fully_converged=False` (while flat `converged=True`) | ✅ |
| **…and the user is TOLD** | **0 warnings emitted** | ⛔ |

⟹ **F1's "…and emitting nothing" is STILL LIVE.** Every fact the done-when
demands is already on the returned object; the reader has to know to ask for it.

**The mechanism, one line.** `orpheus/sn/solver.py:455`:

```python
if history.converged:      # the TOP level only
    return
```

`converged` is the outermost level; `fully_converged` is the tree. A converged
outer standing on a truncated inner takes the early return, in silence — from
the function whose own docstring says *"so a truncated solve announces itself
once, in one voice, wherever it came from."* That sentence is present-tense
FALSE for every non-leaf tree.

**⛔ NOT a regression — and the attribution matters** (`[M]` `git log -L 455,456`).
The guard has read `history.converged` since `d9b027d7`, the FIRST #340 commit,
when it was the only property that existed. `fully_converged` arrived with
**N2b-ii** (`fab9e82a`) — which minted the tree-wide property and did **not**
re-point the one consumer whose entire job is making that gap audible. This is
the `coding-standards` retirement-MIRROR rule ("landing a deferred capability
stales its deferral contract") firing inside our own campaign.

**⚠ The existing audibility gates cannot see this, by construction.**
`TestTruncationIsAudible` (`tests/sn/solve/test_convergence_contract.py:361`) has
four fixtures and **all four** route through `_fixed_source` →
`solve_sn_fixed_source`, which §1b measures as a **1-level** tree. On a leaf,
`converged ≡ fully_converged` identically, so no fixture in the class can
distinguish the two properties. Mode-12 exactly: the measured functional
annihilates the property. The eigenvalue path — the only shape with a tree — has
**no audibility gate at all**.

**⛔ §4's "warn per truncated inner" rejection does NOT cover this.** That row
names a MECHANISM (warn once per inner, "fires up to `max_outer` times") and
rejects it on frequency. Warning **once at the public entry when the TREE did
not converge** fires exactly once per solve — the same cadence as today's
warning — so the stated reason does not apply. The row also cites
`test_converged_solve_is_SILENT` as "already gated"; that gate is one of the
four leaf fixtures above and is structurally blind to the case. Left unmarked,
the row would tell the next session this whole area was settled.

**⚠ The design is NOT a one-token flip, for two reasons.**

1. **The message would lie.** `where`/`budget_name`/`budget`/`tol` are passed in
   by the caller and describe its **top** level — the eigenvalue entries pass
   `max_outer`/`keff_tol` (`solver.py:2465, 2742`). Reporting an inner failure
   through them says *"hit max_outer=500 without reaching tol=1e-07"* about a
   solve whose outer converged fine. Everything needed is already on
   `first_failure` (`budget`, and `tolerance` off `binding_criterion`) — the ONE
   fact the record lacks is the caller-facing **knob name** for that level. It
   belongs on the record, not threaded from the entry: reading it off `label`
   would be the stringly-typed dispatch `solution.py:203` already rejects. Land
   it and three parameters retire from five call sites
   ([[feedback-lossy-return-type-is-the-root-cause]] — the consumer was handed
   facts the producer already knew).
2. **⭐ `[M]` THE POPULATION IS 25 CALLS ACROSS 20 TESTS** — full slice
   `tests/{sn,numerics,cp,moc,diffusion}` at `-m "not slow"`, SERIAL, 2026-08-10,
   HEAD `52650a86`; run outcome **9 failed / 5501 passed** in 1492.92 s, the 9
   being the known task-#51 baseline, so the instrument perturbed nothing. Full
   table: `scratch/issue_340_truncated_inner_population.md`. **Every one is
   `inner(source-iteration)`** — zero `inner(gmres)`, zero CP/MoC/diffusion
   (those carry no record yet: N4). ρ ranges 0.889–0.993.
   (An earlier `-x` run reported 11; it stopped at the first baseline red. The
   11 is a LOWER BOUND and is superseded by the 20.)

   ⚠ **And R2 is the same population viewed twice.** Four of the twenty are
   `test_inner_tol_bias_collapses_at_1e_12[...]` — a test whose PURPOSE is to
   study inner-tolerance bias, so its starved inner is the fixture, not a
   defect. R2's "declare the 7 deliberate truncations in-test" is therefore the
   mechanism that lets commit 2 stay silent about the intentional ones **on
   purpose** rather than by accident. Do R2 before flipping the guard.

3. **Noise is a real risk and my own fixture argues for it.** `[M]` same probe:
   at `max_inner=200`, with 3 of 3 inners TRUNCATED, `keff` is still correct to
   **2.5e-11** against the independent 0-D `k_inf` (`solve_homogeneous_infinite`,
   `0.43846153846153835`); at `max_inner=1000` it is `3.2e-13` and
   `fully_converged` flips True. So a truncated inner is sometimes **benign** —
   and warning on all of them is the noise `test_converged_solve_is_SILENT`'s
   docstring warns produces filtered warnings.

**⛔ CORRECTED 2026-08-10, mid-design — `first_failure` is CHILDREN-FIRST, and
that changes commit 1's blast radius.** I designed the split below believing it
checked self before children (a misreading of N1's MA3 *mutation*, which is the
thing that reddened 0 — the shipped order is the opposite). The code is explicit
(`convergence.py:844-851`): *"Children are searched before `self`… This runs even
when `self` converged, because that is exactly the measured failure."* Two
consequences:

- **Good:** `first_failure` is ALREADY built for N6. Nothing about it needs
  changing; only the guard and the message's sourcing do.
- **⚠ So commit 1 is NOT behaviour-neutral by construction.** Whenever the top
  level fails *and* a child also fails, the guard already lets the warning
  through, and re-sourcing the message from `first_failure` re-points it at the
  child. `[M]` constructed: `solve_sn(max_outer=3, keff_tol=1e-12,
  max_inner=200, inner_tol=1e-8)` warns today as *"hit max_outer=3 … tol=1e-12"*
  and under N6 reads *"inner(source-iteration) … max_inner=200 … tol=1e-8"*.
  **The new message is the correct one** — with a starved inner, raising
  `max_outer` cannot help (F2), so today's advice points at the wrong knob.
- **`[M]` and in the measured slice the change is nil**: every currently-warning
  solve in `tests/sn`+`tests/numerics` is a `solve_sn_fixed_source` **leaf**
  (2 of 2 warnings), where `first_failure is record` and the message is
  byte-identical. So commit 1 is behaviour-neutral **as measured**, not by
  construction — and that distinction is the gate to write: a fixture where the
  top AND a child both fail is the one that separates the two sourcings.

⭐ **Which reframes N5.** What separates F1's corrupting truncation from this
probe's benign one is a question about the ANSWER, not about iteration counts —
and only a residual certificate can ask it. N5 is therefore not "an extra check";
it is **the discriminator that makes this warning non-noisy**. The two compose:
N6 makes the report honest now ("the outer's stop test reads INCREMENTS, which a
starved inner suppresses, so this keff *may* be a stall"); N5 later turns *may*
into *is* / *is not*.

⚠ **And N6 is coupled to N4's MoC decision.** Per the N4 survey, MoC's inner has
no criterion at all, which `IterationRecord.converged` deliberately refuses
(`convergence.py:798-801`) — so recording MoC as-is makes `fully_converged`
**False on every MoC solve forever**. Under a tree-wide guard that is a warning
on every MoC solve. N4's MoC arm and N6's guard must be decided together.

---

### N5 — ⛔ REFUTED AS A GATE, 2026-08-10. The residual cannot discriminate.

**Goal (unchanged).** A caller can tell a truncation that *corrupted the
answer* from one that did not.

**Proposed means (2026-08-09) — ⛔ REFUTED.** *"Outer certificate (F7/F7′) —
residual `(A − F/k)ψ` at exit."* The §3d row and the N6 section's *"N5 is the
discriminator that makes this warning non-noisy"* both rest on it. They are
wrong, and the text stays so the next session does not re-derive the dead end.

`[M]` `scratch/n5_outer_certificate_measurement.md`, **38 solves** — 8
geometries (leaking d1/d2/d2-thick/d2-all-vacuum/d3 + all-reflective
d1/d2/d3) × 3 mixtures, `keff_tol = 1e-7`. Defect `= ‖A_lossψ − Fφ/k‖ /
‖Fφ/k‖`:

| population | n | defect range |
|---|---|---|
| CORRUPTING | 16 | 1.185e-07 … 1.934e-04 |
| BENIGN | 22 | 6.108e-10 … 7.511e-05 |

**A 634× overlap.** A threshold admitting every benign case misses **15 of
16** corrupting ones — including F1, the case it was proposed for. Robust to
the classification standard (relative: 698×, misses 18/19). It is not even
**monotone**: `mi=10` reads a defect 1.8× HIGHER than `mi=3` while the keff
error is 6× LOWER, and `mi=20` reads a defect *below* the fully-converged
reference (1.18e-07 vs 3.47e-07) — so it does not order truncated-vs-converged
at all. `T=1e-7` reaches 16/16 only at a **59 %** false-alarm rate.

⭐ **The structural reason, measured not inferred.** `‖r‖` is up to **99.995 %**
reflective-trace rows, and a reflective inflow-trace defect in a zero-leakage
system carries **no net current** — so `k = production/absorption` is blind to
it *by conservation*. The transfer gain `|Δk|/defect` therefore spans
**1.16e5×** (1.15e-05 on all-reflective d=3 → 1.34 on the all-vacuum slab,
where the trace fraction is 0). **No constant can convert a flux-space
residual bound into an eigenvalue-space one.**

This is `vv-principles` **Mode 12 in the MIRROR**, and it is worth naming as a
distinct shape: the usual Mode-12 failure is a gate that is BLIND to the error
class. Here the gate is *sighted* on a class **the contract is blind to** — it
measures real defect that the measured functional annihilates. Both directions
produce a useless gate; only the second one produces a gate that also looks
busy and reports numbers.

Two further refutations from the same measurement:

- ⛔ **`_CERTIFICATE_SAFETY × tol` is the wrong constant AND the wrong
  quantity.** The outer's binding criterion is `dphi` in *every* solve
  measured, so copying `record.binding_criterion.tolerance` silently picks the
  **looser** of the two (5× weaker). Literal transfer catches 2/16 at 2/22 FA,
  and misses the population's largest keff error (`|dk| = 23.7 × keff_tol`).
- ⛔ **The null case is not a floor.** A fully-converged solve reads
  `3.47e-07 = 3.47 × keff_tol`, and it falls 8 decades under tightening. A
  two-legged sweep localises it: `inner_tol` 1e-8→1e-14 at fixed outer moves it
  **0.2 %**; the outer at fixed inner moves it **6 decades**. It is the power
  iteration's own increment-stop slack, not the inner's.
- ⛔ **An adjoint SHORTCUT is worse than nothing.** First-order perturbation
  theory is the right statistic, but with a spatially-FLAT 0-D adjoint
  (positive-controlled: `|k_pencil − k_inf| = 0`) the overlap *degrades*
  4.64× → **128.95×**. A signed projection against a wrong weight manufactures
  near-cancellations — i.e. false NEGATIVES. Only a spatially-resolved adjoint
  could gate, at one adjoint solve per certificate.

### N6b — the guard widens UNCONDITIONALLY; the projection is a NUMBER, not a gate

**User ruling 2026-08-10**, taken after the refutation above: *flip it, and
report the balance projection.*

**Goal.** A solve whose subtree did not converge says so, once, at the entry —
and carries the one statistic that actually correlates with the error the
caller cares about, labelled as a diagnostic rather than a verdict.

**Means.**

1. `_warn_if_unconverged`'s guard becomes `history.fully_converged`. `[M]` this
   reds **nothing** today and emits **27** further warnings.
2. The message carries the **balance projection**
   `R_g = Σ_n w_n Σ_i V_i r[n,g,i]` — the residual projected onto the
   functional `keff` actually reads. `[M]` it cuts the overlap 634× → **4.64×**
   (14/16 at 2/22 FA), which is exactly what the structural argument predicts:
   projecting onto the functional removes the trace rows the functional
   annihilates. ⚠ Ship it as a **number in the message**, never as a
   threshold — 4.64× still overlaps.
3. The ~20 rows in `scratch/issue_340_truncated_inner_population.md` get R2
   declarations. `[M]` one is already done (R4, `3f76d651`).
4. `test_a_converged_outer_on_a_starved_inner_is_STILL_SILENT` is **deleted**
   and its `xfail(strict=True)` sibling `test_the_tree_wide_truncation_is_
   audible` becomes an ordinary passing row. They were committed as a pair for
   exactly this moment.

⚠ **Where the projection LIVES is the one open design question.** `[M]`
`solve_sn` discards the solver it builds, and `Solution` carries neither
`mat_xs` nor `scattering_op` nor `fission_op` — so the projection cannot be
computed at the entry from public state, only where the operators are still in
scope. It therefore belongs as a measured field on **`IterationHistory`** (SN
local — it needs volumes and quadrature weights, which `numerics` has no
business knowing), NOT on `IterationRecord`.

That re-adds a stored field to the view N2b-ii deliberately emptied, and the
campaign's own rule is what licenses it: *no convergence **verdict** may be
stored; measured data plainly must be* (the same licence `iterations_run`
already holds). Keep the distinction explicit at the field, or the next
reader will read it as the drift N2b-ii removed.

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
