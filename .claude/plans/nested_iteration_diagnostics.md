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
| F2 | The outer's stop test is **entirely increments** — `dk < keff_tol and dphi < flux_tol` (`solver.py:1553`). A throttled inner suppresses the very quantities the outer reads. It does not converge; it **stalls and reads its own stall as convergence**. | source |
| F3 | `n(tol) = \|ln(tol/r0)\| / \|ln rho\|` holds with **R2 >= 0.9956** over 8 configs; fitted slope = `1/\|ln rho\|` to 0.6 % on the 6 well-conditioned rows. **The budget is derivable; `rho` is the only unknown.** | `scratch/probe_inner_budget_law.py` |
| F4 | Both shipped defaults are **SHORT** at d=3 zero-leakage — incl. with scattering: need@1e-8 = 830 (c=0.5) vs `have 200`; need@1e-12 = 1441 vs `have 1000`. | same |
| F5 | The tolerance factor between the two default families should be **~1.5-1.75x** by F3; the tree ships **5x**. | same |
| F6 | CP / MoC / Diffusion receive `outcome.converged` and **drop it**; none of `CPResult` / `MoCResult` / `DiffusionResult` has the field. | `cp/solver.py:884`, `moc/solver.py:139`, `diffusion/solver.py:409` |
| F7 | `_certify_within_group_exit` is wired at 4 sites, **all within-group**. The outer has no certificate. The antidote concept is in the tree, applied one level down only. | `solver.py:1741, 1925, 3418, 3632` |
| F8 | The per-inner status is **available and discarded**: `solver.py:1751` / `:1934` do `self._total_inner_iterations += len(_residuals)` — summing the COUNT, dropping the STATUS. | source |

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

⚠ **This is the design's load-bearing claim and it is UNVERIFIED.** It predicts
zero test churn from the meaning question. The running truncated-exit audit
(#340 bullet 2) measures exactly the set that would churn. **Do not implement
3b before that lands.**

### 3c. The budget: derived, and the warning says what to set

RULED 2026-08-09 by the user: *derived default + rate-aware warning.*

```python
_SERVED_RATE = 0.986        # [M] the slowest rho we promise to serve
def _default_inner_budget(tol: float) -> int:
    return ceil(log(tol) / log(_SERVED_RATE))
```

Retires all five constants (F5's unprincipled 5x with them). Signature becomes
`max_inner: int | None = None` — `None` = derive, an int = deliberate override.

The warning then reads `first_failure` and reports the **observed** rate and
projection, turning *"raise the budget"* into *"set `max_inner=1710`"*. This is
what makes the default's exact value stop being load-bearing — and F3 is what
makes the projection honest rather than a guess.

⚠ `_SERVED_RATE = 0.986` covers the representative d=3 worst (`rho = 0.9854`)
but **NOT** `Sigma_t/4` (`rho = 0.99575`, needs 5171 @1e-12). That is a
deliberate, stated limit, not an oversight: serving 0.996 means a *stuck* d=1
slab churns ~6900 sweeps before admitting defeat. The rate-aware message is what
makes the uncovered tail cheap — it tells that user the number.

### 3d. Sizing / order

| step | scope | gate |
|---|---|---|
| N1 | `IterationRecord` + its own laws (nesting, derivation, rate on synthetic geometric decay) | pure-math unit gates, no solver |
| N2 | SN produces it; `IterationHistory` composes it; flat accessors become derived | the 9-red baseline unchanged |
| N3 | Derived budget + rate-aware message (3c) | contract suite + a NEW starved-projection gate |
| N4 | CP / MoC / diffusion carry it (F6) | their own smokes |
| N5 | Outer certificate (F7) — residual `(A - F/k)psi` at exit | mutation battery |

N1+N3 are independent of the audit. **N2's `fully_converged` waits on it.**

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
