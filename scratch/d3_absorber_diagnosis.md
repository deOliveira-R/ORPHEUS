# d3 pure-absorber exactness violation — diagnosis

Gate: `tests/sn/solve/test_d3_admission.py::test_d3_pure_absorber_per_ordinate_psi_exact`
Branch `refactor/operator-strategy-layers`, HEAD `3308c02c`. Measured 2026-08-08.

---

## 1. VERDICT: (A) — a solver/certificate hole. Not a discretization floor.

**The solve never converges. It runs out of iterations, returns a best-effort
answer flagged `history.converged = False`, and the gate never reads the flag.**

`[M]` probe A, the default call exactly as the gate makes it
(`inner_tol=1e-13`, `max_inner` defaulted to 1000):

```
n_inner   = 999          converged = False
last 5 residuals = 1.154e-09 1.153e-09 1.188e-09 1.209e-09 1.185e-09
g0 max rel err = 3.2873e-10       g1 max rel err = 1.1161e-15
```

### 1.1 The decisive experiment (the brief's own)

`[M]` probe F — honest residual `r = A psi - q` through the PRODUCTION loss
operator `A = L + C - S - B` (`_bare_loss_arm(build_within_group_system(...))`),
on the same `SNMesh` instance:

| state | rel `‖Aψ − q‖` | bulk err | boundary / interior split |
|---|---|---|---|
| returned (max_inner=1000, THE GATE) | **1.123e-09** | 3.287e-10 | bnd 3.79e-09, int 2.30e-15 |
| converged (max_inner=5000, n=1631) | 9.793e-14 | 2.790e-14 | bnd 3.31e-13, int 2.31e-15 |
| **EXACT uniform field (bulk AND trace)** | **1.060e-15** | 0 | bnd 0.0, int 3.58e-15 |

The closed form has a **machine-zero** equation residual. The returned iterate
is four orders away. The discretization is exonerated; the solve failed to
arrive. **(A).**

The interior residual is `~2.3e-15` at *every* state — the defect lives
entirely in the **boundary** block, i.e. in the reflective coupling that the
lagged iteration is still resolving.

### 1.2 Forcing more iterations fixes it

`[M]` probe B, `max_inner` sweep at the same `inner_tol=1e-13`:

| max_inner | n_inner | converged | last residual | err g0 | err g1 |
|---|---|---|---|---|---|
| 1000 (default) | 999 | **False** | 1.185e-09 | **3.2873e-10** | 1.1161e-15 |
| 3000 | 1631 | True | 9.793e-14 | **2.7903e-14** | 9.7660e-16 |
| 10000 | 1631 | True | 9.793e-14 | 2.7903e-14 | 9.7660e-16 |
| 30000 | 1631 | True | 9.793e-14 | 2.7903e-14 | 9.7660e-16 |

Converged error **2.79e-14** — four orders *inside* the gate's `rtol=1e-10`,
and stable from 3000 upward.

### 1.3 Why the tolerance sweep looked like a floor (the false fingerprint)

The prior measurement — bit-identical `3.287e-10` at `inner_tol` ∈ {1e-9, 1e-11,
1e-13, 1e-15} — is fully explained. The running residual at k=999 is
**1.185e-09**, above *every* one of those tolerances **including 1e-9**. All
four runs therefore hit the same `max_inner=1000` cap and return the identical
999th iterate. Tolerance-insensitivity here is the signature of an
**iteration-count truncation**, not of a discretization floor.

This is `numerical-bug-signatures` Signature 8 in a new costume: L10's
"discarded library info flag" has become "a *recorded but unread* convergence
flag". Same false "diverges/plateaus independent of tolerance" fingerprint.

### 1.4 The decay is geometric, and the ordinate map is a MODE, not a bias

`[M]` probe B item 2 and probe G. The residual decays cleanly at
`rho = 0.9853` (`0.0525^(1/199)`), `1/(1-rho) = 68`. Nothing stagnates; k=999
is mid-descent. The last-5 non-monotonicity is the octant-cyclic G-S ripple on
the geometric envelope.

`[M]` probe G — the per-ordinate error map vs `max_inner`:

| max_inner | n | converged | max err g0 | x-dom mean | other mean |
|---|---|---|---|---|---|
| 600 | 599 | False | 9.7120e-08 | 6.2370e-08 | 4.2267e-13 |
| 800 | 799 | False | 6.3838e-09 | 3.5213e-09 | 1.3690e-15 |
| 1000 | 999 | False | 3.2873e-10 | 1.9865e-10 | 1.1510e-15 |
| 1200 | 1199 | False | 1.6574e-11 | 9.1784e-12 | 1.1074e-15 |
| 1400 | 1399 | False | 8.0040e-13 | 5.1402e-13 | 1.1336e-15 |

Successive `rho` = 0.98648, 0.98528, 0.98517, 0.98496 — constant. Shape cosine
similarity vs the k=599 map ≥ 0.9928 at every count; the x-dominant fraction of
the map's norm is **1.000000** throughout.

⟹ the "8 of 24 ordinates carry the error, all with |μ_x| = 0.8689" structure is
**the shape of the slowest-decaying SI eigenmode**, not a per-ordinate
discretization defect. A bias is fixed under iteration; this decays 20× per 200
sweeps while keeping its shape exactly.

x is the SHORTEST axis (extent 1.0), so the x-dominant ordinates have the
shortest optical path per reflective traversal (`Σ_t·Lx/|μ_x| = 0.921` vs
`2.286` for the |μ_x| = 0.350 ordinates) — the least-attenuated, hence
slowest-decaying, reflective loop.

### 1.5 It is not the splitting and not the inner solver

`[M]` probe B items 3–4, all at `max_inner=3000`. All three arms land on the
same fixed point (Mode-9 FP-invariance intact):

| arm | n_inner | converged | err g0 | err g1 |
|---|---|---|---|---|
| `gauss_seidel` (default) | 1631 | True | 2.79e-14 | 9.77e-16 |
| `jacobi` | 838 | True | 4.70e-14 | 1.26e-15 |
| `krylov` | 733 | True | 6.56e-13 | 1.34e-13 |

---

## 2. The exact code path, with `file:line`

1. `orpheus/sn/solver.py:2981` — `solve_sn_fixed_source(..., max_inner: int = 1000, ...)`.
   A fixed, dimension-blind, cross-section-blind iteration cap.
2. `orpheus/sn/solver.py:3235` — `_within_group_si(..., max_iter=max_inner, tol=inner_tol)`;
   the `SourceIteration` primitive returns its best effort after `max_iter`
   without raising.
3. **`orpheus/sn/solver.py:427-428`** — the certificate's own early return:
   ```python
   if not residual_history or not (residual_history[-1] < tol):
       return  # no convergence claim — nothing to certify
   ```
   `_certify_within_group_exit` is **deliberately** a no-op on a `max_iter`
   exit; its docstring states *"No-op when the exit made NO claim (`max_iter`
   hit without reaching `tol`) — best-effort returns stay legal"*.
4. `orpheus/sn/solver.py:3283` —
   `converged_flag = bool(residuals) and residuals[-1] < inner_tol` → `False`,
   recorded honestly on `IterationHistory`.
5. `tests/sn/solve/test_d3_admission.py:151-165` — the gate consumes
   `sol.angular_flux` and asserts `rtol=1e-10`. **It never reads
   `sol.history.converged`.**

So no single line is *wrong*. The hole is the **contract**: a best-effort
answer and a certified answer are returned by the same call, in the same type,
and are indistinguishable unless the caller volunteers to check a flag that
almost nobody checks. That is the L10 / Signature-8 anti-pattern one level up.

---

## 3. Why commit `59bb38a0` (#337) moved this gate — and what it reveals

`[M]` probe C. The pre-#337 cosines were monkeypatched in-process (never a git
revert); weights are bit-identical (`Σw = 12.566370614359172`) in both rules, so
the cosines are the only variable.

```
PRE-#337  |cos| = [0.40824829 0.81649658]   max|Σμ²−1| = 1.11e-16
HEAD      |cos| = [0.35002117 0.86889030]   max|Σμ²−1| = 1.11e-16

mechanism, group 0 (Σ_t = 0.8, Lx = 1.0):
PRE-#337  |μ_x|max = 0.8164966  τ_x = 0.979796  exp(−τ) = 0.375388
HEAD      |μ_x|max = 0.8688903  τ_x = 0.920715  exp(−τ) = 0.398234
```

| rule | max_inner | n_inner | converged | last residual | err g0 |
|---|---|---|---|---|---|
| PRE-#337 | 1000 | 999 | **False** | 7.065e-11 | **2.0487e-11** |
| PRE-#337 | 5000 | 1369 | True | 9.610e-14 | 2.5671e-14 |
| HEAD | 1000 | 999 | **False** | 1.185e-09 | **3.2873e-10** |
| HEAD | 5000 | 1631 | True | 9.793e-14 | 2.7903e-14 |

**The gate has ALWAYS been riding an unconverged best-effort answer.** Pre-#337
it was green with a 5× margin (2.05e-11 vs `rtol=1e-10`) purely by luck. #337
raised the largest cosine by 6.4%, which shortened the per-traversal optical
path on the short x axis, raised the reflective survival, raised the required
sweep count 1369 → 1631 (+19%), and so raised the k=999 error 2.05e-11 →
3.29e-10 (16×) — across the threshold.

**#337 is exonerated.** It changed a quadrature, correctly and as ratified; it
merely stopped a latent test-design defect from hiding.

---

## 4. The mechanism: why an all-reflective pure absorber is SI-hard

With `Σ_s = 0` and every face reflective there is **no scattering iteration and
no leakage** — the only coupling left is the boundary `B`, and the only damping
is absorption.

The slow mode is the **DD face sawtooth**. `[M]` probe E, the CONVERGED
(n=1631) state's `xmin` trace, group 0, `want = 9.9471839432e-02`:

```
spatial variation on xmin (ordinate 0), ratio to want, across (ny, nz):
[[1.074414 1.074414 1.074414 1.074414 1.074414]
 [0.925586 0.925586 0.925586 0.925586 0.925586]
 [1.074414 1.074414 1.074414 1.074414 1.074414]
 [0.925586 0.925586 0.925586 0.925586 0.925586]]
```

`1.074414 + 0.925586 = 2.000000` exactly. DD pins only the cell **average**
(`ψ_c = (ψ_in + ψ_out)/2`), so a face sawtooth about `ψ*` leaves every cell
average exactly `ψ*` — the mode has **zero cell average**, hence the collision
term `Σ_t V ψ_c` does not see it at all. It is damped only by the inter-axis
balance mismatch around the reflective loop, which is weak. Hence `rho → 1`.

Two consequences, both measured:

* **`Σ_t · n_inner` is invariant** — the damping is absorption-limited:
  `[M]` d=3, boundary-G-S, `inner_tol=1e-13`:
  `Σ_t = 0.4 → 3093`, `0.8 → 1631`, `1.6 → 850`, `3.2 → 437`
  (products 1237 / 1305 / 1360 / 1398 — flat to 13%).
* **Group asymmetry explained**: g1 (`Σ_t = 1.6`) needs ~half the sweeps of g0
  (`Σ_t = 0.8`), so it is long converged at k=999 (1.1e-15) while g0 is not.

Note also (probe E/F): because the sawtooth is very nearly null, the
converged state's trace retains a **11.26 %** deviation from uniform
(`xmin trace dev = 1.1258e-01` at both max_inner=1000 and 5000) while the
exact-uniform state also has residual 1.06e-15. Two near-solutions differing
only in a zero-cell-average trace mode ⟹ the all-reflective 3-D DD operator has
a (near-)null space in the trace block. Benign for this gate — the cell
averages and the scalar flux are uniquely and correctly determined — but it is
the structural reason `rho ≈ 0.985`, and it is worth a theory-page note.

---

## 5. Recommended remedy

**Both a production fix and a gate re-derivation are required. The production
fix is primary (Cardinal Rule 1).**

### 5.1 Production (primary) — a claimed exit must be distinguishable from a best-effort one

The certificate machinery landed in `c98a23d8` to make a *claimed* exit honest.
It does that. The remaining hole is the *unclaimed* exit: `max_iter` truncation
returns silently. Options, in preference order:

1. **Make the best-effort exit loud at the entry.** `solve_sn` /
   `solve_sn_fixed_source` should refuse (or at minimum `warnings.warn`) when
   `history.converged is False`, with an opt-out
   (`allow_unconverged=True` / `on_unconverged="raise"|"warn"|"ignore"`) for the
   legitimate consumers — the DSA rate studies at
   `tests/sn/acceleration/test_dsa_rate.py` deliberately cap `max_inner=50`.
   A `ConvergenceCertificateError` sibling (`InnerIterationBudgetExceeded`) keeps
   the existing exception vocabulary.
2. **Derive `max_inner` instead of hardcoding it** (`coding-elegance` Pattern 7 /
   anti-pattern #14 — the magic constant). The measured law
   `n ≈ 1300 / Σ_t,min` at d=3 is fixture-specific, but the *shape*
   (`∝ 1/Σ_t`, growing sharply with dimension) means one fixed integer cannot
   serve d=1 (32 sweeps) and d=3-all-reflective (1631) at once. At minimum,
   scale the default with `ndim`.
3. Do **not** simply raise the default to, say, 5000 and leave it silent — that
   moves the cliff without removing it.

### 5.2 The gate (secondary but required) — its own premise is false

`tests/sn/solve/test_d3_admission.py:132-134` states:

> *"DD is exact for flat flux and **c=0 needs no iteration**, so EVERY ordinate
> must carry the closed-form value to solver-tolerance precision."*

**`c = 0` kills the SCATTERING iteration, not the REFLECTIVE-BOUNDARY
iteration.** On an all-reflective box the boundary coupling is its own lagged
fixed-point problem with `rho = 0.9853`, needing 1631 sweeps. The docstring's
rationale is the test-design defect; the tolerance is downstream of it.

The principled fix is **not** to relax `rtol`. It is to **give the solve the
iteration budget its own stopping criterion asks for, and then assert the tight
tolerance that is genuinely delivered**:

```python
sol = solve_sn_fixed_source(
    {0: mix}, _d3_axes(), quad,
    external_source=q, boundary_condition="reflective",
    inner_tol=1e-13, max_inner=4000,          # [M] 1631 needed; ~2.5x headroom
)
np.testing.assert_equal(bool(sol.history.converged), True)   # the missing check
...  # then the existing rtol=1e-10 assertions, which now pass at 2.8e-14
```

`rtol=1e-10` is the *right* number and must stay: with a converged solve the
delivered error is **2.79e-14**, so the gate keeps ~3.5 orders of margin and
remains a real Mode-1/3/4 catcher. Relaxing it to accommodate 3.3e-10 would be
exactly `vv-principles` anti-pattern #16 — asserting around a producer defect
instead of quoting the guarantee.

**Every gate that asserts a tolerance on an SN `Solution` should assert
`history.converged` first.** That one line is what turns the whole class of
silent budget truncations from invisible into loud.

---

## 6. Refuted candidates added by this investigation

Each with its structural reason (per `process-discipline`):

* **A discretization / balance-equation defect introduced by #337.** REFUTED:
  the exact uniform field's residual through the production operator is
  `1.06e-15` (machine zero) on the HEAD quadrature; the discrete system admits
  the closed form exactly. Probe F.
* **A per-ordinate redistribution / routing bias on the 8 x-dominant
  ordinates.** REFUTED: the per-ordinate error map decays geometrically at
  `rho = 0.9852` with its shape preserved (cosine ≥ 0.9928) over
  `max_inner` = 600…1400 — a decaying eigenmode, not a fixed bias. Probe G.
* **A stagnation floor in the SI iteration (FP fixed point).** REFUTED: the
  residual decays cleanly 5.62e-01 → 9.79e-14 over 1631 iterations at constant
  `rho`; k=999 is mid-descent, not a floor. Probe B item 2.
* **A wrong Gauss-Seidel boundary splitting (ERR-056 class / Mode-9 FP
  violation).** REFUTED: G-S, Jacobi and Krylov all converge to the same fixed
  point (2.79e-14 / 4.70e-14 / 6.56e-13 vs the closed form). The splitting
  changes the RATE only — which is a separate finding, §7.
* **The certificate `_certify_within_group_exit` failing to fire on a claimed
  exit.** REFUTED as a *bug*: it correctly makes no claim because the driver
  made none (`solver.py:427-428`). The hole is the unclaimed exit, not a
  broken certificate.
* **Off-unit-sphere ordinates amplifying `1e-16` to `3e-10`.** REFUTED:
  `max|Σμ² − 1| = 1.11e-16` for BOTH the pre-#337 and HEAD rules
  (probe C) — identical, so it cannot explain a difference between them; and
  the converged answer is exact to 2.8e-14 with the same nodes.
* **A ρ-blind stopping test delivering a false "converged" (Signature 9 /
  L11).** REFUTED for this gate: the solver did **not** claim convergence. The
  failure is the complementary one — an honest non-claim that nobody reads.

---

## 7. Secondary finding — boundary Gauss-Seidel is ~2× SLOWER than Jacobi at d=3

`[M]` probe D, `inner_tol=1e-13`, level-symmetric S4:

| case | G-S | Jacobi | G-S / Jacobi |
|---|---|---|---|
| d=2 box 1×2 (3,4) all-reflective | **258** | 648 | 0.40 (G-S wins, as documented) |
| d=3 box 1×2×3 (3,4,5) all-reflective | **1631** | 838 | **1.95 (G-S LOSES)** |

`rho`: G-S 0.985348 (`1/(1−ρ) = 68.3`), Jacobi 0.975014 (`1/(1−ρ) = 40.0`).

`_select_si_splitting` (`orpheus/sn/solver.py:707`) documents the boundary-G-S
arm as *"a modest reflective-SI rate gain, ~0.86–0.92× on B-mixture configs"*
and C5.5 (#225) extended the `is_cartesian and not is_1d` gate to d=3. The
**fixed point is invariant** (verified above, so no correctness issue), but the
**rate claim was measured at d=2 and does not transfer to d=3** — it inverts.
With `inner_schedule="jacobi"` this gate would need 838 < 1000 sweeps and would
have been green all along.

Worth a `module:sn` / `type:improvement` issue: either re-derive the d=3 G-S
octant schedule, or restrict the G-S default to d=2 until the d=3 ordering is
re-derived. Note this also means the d=3 default is currently the *slower* of
the two arms on reflective configs.

---

## 8. Blast radius — is the 3-D box special?

`[M]` probe D (c), all-reflective pure absorber, `Σ_t = (0.8, 1.6)`,
`inner_tol=1e-13`, boundary-G-S, level-symmetric S4:

| case | n_inner | converged | err vs closed form |
|---|---|---|---|
| d=1 slab Lx=1, 3 cells (GL8) | 35 | True | 3.649e-13 |
| d=1 slab Lx=1, 3 cells (LS4) | 32 | True | 2.785e-13 |
| d=2 box 1×2, (3,4) (LS4) | 258 | True | 6.641e-14 |
| **d=3 box 1×2×3, (3,4,5) (LS4)** | **1631** | True | 2.790e-14 |
| d=3 with x-VACUUM, y/z reflective | 208 | True | (n/a — not a flat problem) |

**Yes, d=3 all-reflective is special — by roughly an order of magnitude per
added dimension** (32 → 258 → 1631). d=1 and d=2 sit comfortably inside
`max_inner=1000`; d=3 does not. Adding a single vacuum face collapses the cost
to 208 (leakage damps the sawtooth), so the trigger is specifically
**all-reflective (zero-leakage) + weak absorption + d=3**.

Who else can be bitten, in order of likelihood:

1. **Any d=3 all-reflective fixed-source gate at low `Σ_t`.** The cost is
   `≈ 1300/Σ_t`, so `Σ_t = 0.4` already needs 3093. The sibling
   `test_d3_scattering_infinite_medium_matches_multigroup_balance` in the same
   file (mixture C, `rtol=1e-7`) is the immediate neighbour — it carries a
   `c > 0` scattering iteration ON TOP of the reflective one, so its budget is
   strictly larger; its looser `rtol=1e-7` is currently absorbing it. It should
   be given the same `converged` assertion.
2. **`solve_sn` (the eigenvalue entry) defaults to `max_inner = 200`
   (`orpheus/sn/solver.py:2034`)** — 8× tighter than the fixed-source entry.
   That is *defensible* for a power iteration (an inexact inner is a legitimate
   inexact-Newton posture and the outer re-solves). What is **not** defensible
   is that `orpheus/sn/solver.py:2262-2268` hardcodes:

   ```python
   history = IterationHistory(
       keff_history=tuple(keff_history),
       n_outer=len(keff_history),
       total_inner_iterations=solver._total_inner_iterations,
       converged=True,                 # <-- unconditional
   )
   ```

   `power_iteration` (`orpheus/numerics/eigenvalue.py:205`) is a plain
   `for n in range(1, max_iter + 1)` that returns whatever it has; it does not
   raise on exhaustion. So `solve_sn` reports `converged=True` even when the
   OUTER power iteration hit `max_outer`, and it cannot report inner truncation
   at all.

   **This is a twin-path divergence** (`coding-elegance` Pattern 2): the adjoint
   sibling `solve_sn_adjoint` (`orpheus/sn/solver.py:2442`) gets it right at
   `:2536` — `converged=len(keff_history) < max_outer`. The forward and adjoint
   entries spell the same concept two ways and only one is honest. Worth its
   own issue, and the fix is to single-source the predicate.
3. **Not curvilinear-specific.** Nothing here is geometry-dependent beyond
   dimensionality and the reflective closure; a 1-D reflective slab is 50× cheaper.

Suggested issue set (`module:sn`):
* `type:bug` — a best-effort `max_inner` exit is indistinguishable from a
  certified one at the public entries (§5.1); includes auditing every SN gate
  for a missing `history.converged` assertion.
* `type:improvement` — boundary-G-S is ~2× slower than Jacobi at d=3 (§7).
* `type:bug` — `solve_sn` hardcodes `IterationHistory.converged=True`
  (`solver.py:2268`).

---

## 9. Diagnostics written (both green, `python -O -m pytest`)

* `derivations/diagnostics/diag_d3_absorber_01_unconverged_exit.py` — 3 tests,
  32 s. Characterizes the truncated exit; the positive control showing the
  answer IS exact with headroom; the solver-independent residual pin with its
  own perturbation control leg (vv anti-pattern #19).
* `derivations/diagnostics/diag_d3_absorber_02_si_rate_scaling.py` — 2 tests,
  89 s. The general properties: the iteration budget grows an order per
  dimension, and `Σ_t · n_inner` is invariant (absorption-limited damping).

**Promotion recommendation** (per `tests/derivations/_promotion_policy.md`):

* PROMOTE `test_d3_absorber_converges_when_given_enough_iterations` into
  `tests/sn/solve/test_d3_admission.py`, replacing the current red gate body —
  it is the same physics claim, correctly posed, plus the missing `converged`
  assertion.
* PROMOTE `test_d3_absorber_exact_uniform_field_is_the_discrete_solution` —
  a general, solver-independent L0 property (the closed form satisfies the
  discrete operator), with its control leg. Valuable permanently.
* PROMOTE both legs of `diag_..._02` into a new
  `tests/sn/solve/test_reflective_si_iteration_budget.py` (mark `slow`) — these
  are the gates that make a future `max_inner` change measurable, and they pin
  the mechanism.
* `test_d3_absorber_default_max_inner_does_not_converge` — LEAVE as a
  diagnostic until §5.1 lands, then convert to `xfail(strict=True)` so the fix
  forces its deletion (the marker set as the todo list).
