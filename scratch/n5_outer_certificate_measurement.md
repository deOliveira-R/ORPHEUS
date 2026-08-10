# #340 N5 — can an outer convergence certificate discriminate a corrupting truncated inner from a benign one?

**Status: COMPLETE.** Verdict: the certificate is REFUTED as a discriminator (see §1+3).
Measured 2026-08-10 on HEAD `b0a003b4`, host `.venv` (Py 3.14), `python -O`, serial.

Probes (all under `scratch/`, nothing under `orpheus/` touched):

- `scratch/n5_outer_cert_lib.py` — the certificate itself
- `scratch/n5_probe_01_null_case.py` — item 5
- `scratch/n5_probe_02_sweep.py` — items 1, 2
- `scratch/n5_probe_03_benign.py` — items 1, 6

---

## 0. What was built, and the convention boundary (L17)

The proposed certificate, verified against the tree rather than assumed:

- `_certify_within_group_exit` (`orpheus/sn/solver.py:582`) evaluates
  `evaluate_residual(loss_op, psi, q_ext)` → `defect = ‖r‖/‖q‖`, raises above
  `_CERTIFICATE_SAFETY (=10.0) × tol` where `tol = record.binding_criterion.tolerance`.
- `loss_op` at the seedless call sites is `_bare_loss_arm(system)` = the `(A,A)` block
  of the within-group loss grid = **`L + C − S − B_a`** (`coupled_system.py:470`).
- `q_ext` at the SI call site (`solver.py:1805-1831`) is
  `TimedFullField(interior=AngularSourceSink.from_isotropic(fission_source, mesh),
  boundary=AngularBoundarySourceSink.zeros_on(mesh))`.

The brief's reading is **confirmed**: the eigenvalue problem is that same fixed-source
problem with `q = F φ(ψ)/k`, so the outer certificate is the same call with that rhs.

**The scalar↔angular bridge.** `SNSolver.compute_fission_source` (`solver.py:1364`,
`fission_op.apply(φ)/k`) returns a **scalar** emission density `(ng, *spatial)`;
`A_loss.apply` consumes/returns **per-ordinate** `(N, ng, *spatial)`. I did **not**
hand-roll the `1/W`. I call the same producer-side factory the SI inner calls —
`AngularSourceSink.from_isotropic`, documented at `angular_source_sink.py:149` as
applying `/sum_w` and broadcasting over `N`. So the `/W` appears exactly once, on the
rhs, at the producer, in the same place the inner puts it; both sides of the ratio are
per-ordinate source densities and **`W` cancels in the ratio** (`W = 4π = 12.566…` for
LS-S4, reported by every run rather than assumed).

**Which ψ, which φ.** Headline `defect` uses `ψ = solution.angular_flux` (the composite
the *user receives*) and `φ = ψ.interior.integrate_angular()` — the scalar flux **of that
ψ** — so the residual is a closed statement about the returned object alone. A contrast
column `defect_pi` uses `φ = solution.scalar_flux` (the power iteration's converged
scalar, the one `keff` was computed from). They agree to ~10 % everywhere measured; the
choice is not load-bearing.

---

## 5. THE NULL CASE — answered first, because it sets the floor

> config: 2-D extents `(2.0, 1.5)` / cells `(4, 3)`, x-vacuum y-reflective, LS-S4,
> `get_mixture("A","2g")`, `max_inner=5000`, `max_outer=3000`. Per L10b every row prints
> `n_out`/`n_in` so a capped row cannot be misread as a floor — none is capped.

```
 keff_tol  flux_tol inner_tol |               keff  n_out  n_in* |      defect  d/keff_tol      r_bnd      r_int warn
------------------------------------------------------------------------------------------------------------
    1e-07     1e-06     1e-10 |   0.15824724051500      5     80 |  3.4748e-07   3.475e+00  8.903e-07  4.543e-07    0
    1e-09     1e-08     1e-10 |   0.15824724153975      7     80 |  3.0457e-10   3.046e-01  8.044e-10  3.470e-10    0
    1e-11     1e-10     1e-10 |   0.15824724153003      8     80 |  6.8749e-11   6.875e+00  1.807e-10  8.021e-11    0
    1e-13     1e-12     1e-10 |   0.15824724150312     27     80 |  3.7141e-12   3.714e+01  9.548e-12  4.791e-12    0
    1e-14     1e-14     1e-10 |   0.15824724150296     41     80 |  1.2636e-13   1.264e+01  3.305e-13  1.512e-13    0
    1e-07     1e-06     1e-08 |   0.15824724407687      5     64 |  3.5428e-07   3.543e+00  9.105e-07  4.577e-07    0
    1e-07     1e-06     1e-12 |   0.15824724048606      5     96 |  3.4743e-07   3.474e+00  8.901e-07  4.543e-07    0
    1e-07     1e-06     1e-14 |   0.15824724048575      5    113 |  3.4743e-07   3.474e+00  8.901e-07  4.543e-07    0
    1e-14     1e-14     1e-10 |   0.15824724150296     41     80 |  1.2636e-13   1.264e+01  3.305e-13  1.512e-13    0
    1e-14     1e-14     1e-13 |   0.15824724150297     18    105 |  4.0111e-14   4.011e+00  1.030e-13  5.202e-14    0
    1e-14     1e-14     1e-15 |   0.15824724150298     10    121 |  3.7897e-15   3.790e-01  9.246e-15  5.774e-15    0
```

**It is not a floor.** The defect falls ~8 decades (3.47e-07 → 3.79e-15) as the
tolerances tighten, so the certificate has headroom; nothing structural pins it above
zero.

**It is the OUTER's slack, not the inner's.** Rows 6–8 hold the outer at the shipped
default and sweep `inner_tol` over 1e-08 → 1e-14: the defect moves by **0.2 %**
(3.5428e-07 → 3.4743e-07). Rows 1–5 hold the inner and tighten the outer: the defect
moves by **6 decades**. So at the shipped defaults the residual is measuring the power
iteration's own increment-stop slack — which is exactly the L11 statement (the increment
understates the true distance; the residual is what reports it).

**The answer to item 4, stated as the floor any threshold must clear.** At the shipped
defaults (`keff_tol=1e-7`, `flux_tol=1e-6`) a *genuinely converged* solve reads

> **defect = 3.47e-07 = 3.47 × keff_tol = 0.347 × flux_tol.**

`_CERTIFICATE_SAFETY × keff_tol = 1e-06` would therefore pass this solve with only a
**2.9×** margin — thin, and for a reason: the constant was calibrated against a residual
compared to *its own criterion's* tolerance one level down, and here the binding
criterion is `dphi` (tol `flux_tol = 10 × keff_tol`), not `dk`. Measured against
`flux_tol` the same solve sits at `0.347 × flux_tol`, i.e. `10 × flux_tol` gives a
**29×** margin. See §4 for the recommendation.

---

## ⛔ REFUTATION OF THE BRIEF'S BENIGN POLE — the fixture is unphysical

The brief specifies the benign pole as `scratch/probe_340_done_when.py`, and states
*"keff is nonetheless correct to **2.5e-11** against the independent 0-D reference
`solve_homogeneous_infinite(mix).k_inf`"*. **That does not reproduce.** Measured:

```
    k_inf (0-D reference, independent solver family) = 0.23076923076923
    keff  (SN)                                       = 0.30000000000000
    |k - k_inf|            = 6.9231e-02   = 6.923e+05 x keff_tol(1e-7)
```

My reproduction is faithful — running the shipped probe **verbatim** gives
`keff = 0.2999999999999999`, bit-matching my harness. The reference is the problem,
and the root cause is the exact trap the brief warned about, present in the shipped
probe itself. `scratch/n5_probe_04_benign_fixture_audit.py`:

```
s (as literally built by the probe) =
 [[0.28 0.  ]
 [0.12 0.8 ]]
s.sum(axis=0) [column sums, used for sigma_c] = [0.4 0.8]
s.sum(axis=1) [row sums] = [0.28 0.92]
sigma_c = [0.35 0.5 ]   sigma_f = [0.05 0.3 ]   sigma_t = [0.8 1.6]

--- AS SHIPPED   sig_s = s   
    mix.SigS[0] (ORPHEUS reads this as [g_from, g_to]) =
[[0.28 0.  ]
 [0.12 0.8 ]]
    out-scatter per group  sum_to SigS[g,:] = [0.28 0.92]
    CONSISTENCY  sigma_t              = [0.8 1.6]
                 sigma_c+sigma_f+out  = [0.68 1.72]
                 residual             = [ 0.12 -0.12]   -> *** INCONSISTENT ***
    k_inf (transport balance, removal = sigma_t) = 0.23076923076923
    flux (k_inf eigenvector)                     = [833.33333333   0.        ]
    production / absorption  (the SN's k formula)= 0.30000000000000

--- TRANSPOSED   sig_s = s.T 
    mix.SigS[0] (ORPHEUS reads this as [g_from, g_to]) =
[[0.28 0.12]
 [0.   0.8 ]]
    out-scatter per group  sum_to SigS[g,:] = [0.4 0.8]
    CONSISTENCY  sigma_t              = [0.8 1.6]
                 sigma_c+sigma_f+out  = [0.8 1.6]
                 residual             = [0. 0.]   -> CONSISTENT
    k_inf (transport balance, removal = sigma_t) = 0.43846153846154
    flux (k_inf eigenvector)                     = [438.59649123  65.78947368]
    production / absorption  (the SN's k formula)= 0.43846153846154
```

**The mechanism.** `s[1, 0]` carries the comment `# 0 -> 1 downscatter`, i.e. the author
wrote `s` in a `[to, from]` orientation — and `sigma_c = sigma_t - s.sum(axis=0) - sigma_f`
is the *correct* removal for that orientation. But `make_mixture(sig_s=...)` reads
ORPHEUS's `SigS[l][g_from, g_to]`. So `sig_s=s` pairs a `[to,from]` matrix with a
`[from,to]` consumer: `σ_t ≠ σ_c + σ_f + Σ_out σ_s`, off by `±0.12`.

**Why that destroys the comparison, with no bug in either solver.** In a zero-leakage
medium there are two expressions for `k`, and they coincide *only* under that identity:

| expression | removal used | value on the shipped fixture |
|---|---|---|
| transport balance (`solve_homogeneous_infinite`) | `σ_t` | `0.23076923076923` |
| production / absorption (what `SNSolver` reports) | `σ_c + σ_f` | `0.30000000000000` |

The SN's `0.30000000000000` matches the second **to all printed digits** — evaluated on
the *0-D solver's own eigenvector*, an independent chain. So the SN answer is fine and the
truncated inner really did not corrupt it; the brief compared it to the wrong one of two
numbers that an inconsistent mixture makes different.

**Two further defects in the fixture, both fatal to its use as a verification pole:**

1. As shipped, the k_inf eigenvector is `[833.33, 0.0]` — **φ₁ ≡ 0**. With
   `χ = (1, 0)` and no `0→1` transfer, group 1 has no source at all, so the 2-group
   fixture is **effectively 1-group** — `vv-principles` anti-pattern #3, the degenerate
   case that cannot see a spatial/angular/scattering error.
2. `fissile(c=0.9)` (my config C attempt) yields `σ_c = (0.03, −0.14)` — a **negative**
   capture cross-section. The family is only physical for `c ≤ 0.8125`.

**The repair is one character**: `sig_s=s.T`. Everything downstream then agrees
(`k_inf = 0.43846153846154` from both expressions) and φ₁ ≠ 0. The benign pole is
re-measured on the repaired fixture below.

> Action for #340: `scratch/probe_340_done_when.py` is the plan's own done-when probe and
> it carries this fixture. Its `keff` and `k_inf` numbers should not be cited until it is
> repaired.

---

## 2. THE SWEEP — does the defect track the corruption monotonically?  **No.**

The F1 corrupting configuration verbatim: 2-D extents `(2.0, 1.5)` / cells `(4, 3)`,
x-vacuum y-reflective (the LEAKING box — an all-reflective box has `k = k_inf`, which is
flux-shape *independent*, so a starved inner cannot corrupt it; `vv-principles`
anti-pattern #3), LS-S4, `get_mixture("A","2g")`, `keff_tol=1e-7`, `flux_tol=1e-6`,
`inner_tol=1e-10`. `scratch/n5_probe_02_sweep.py`:

```
==============================================================================================================================
REFERENCE (max_inner=5000, keff_tol=1e-7):  keff = 0.15824724051500   defect = 3.4748e-07   fully_converged = True
==============================================================================================================================
max_inner |               keff |dk|/keff_tol |      DEFECT d/keff_tol d/flux_tol |  conv fully  trunc/n  n_in* warn
------------------------------------------------------------------------------------------------------------------------------
        1 |   0.15824739970020    1.5919e+00 |  4.0882e-06  4.088e+01  4.088e+00 |  True False   58/58       1    0
        2 |   0.15824765877130    4.1826e+00 |  2.4456e-06  2.446e+01  2.446e+00 |  True False   33/33       2    0
        3 |   0.15824776630547    5.2579e+00 |  4.0145e-06  4.015e+01  4.015e+00 |  True False   24/24       3    0
        5 |   0.15824742765625    1.8714e+00 |  4.8825e-07  4.883e+00  4.883e-01 |  True False   18/18       5    0
       10 |   0.15824721107207    2.9443e-01 |  8.8030e-07  8.803e+00  8.803e-01 |  True False   11/11      10    0
       20 |   0.15824723143733    9.0777e-02 |  1.1842e-07  1.184e+00  1.184e-01 |  True False    6/6       20    0
       50 |   0.15824724051497    3.6277e-07 |  3.4730e-07  3.473e+00  3.473e-01 |  True False    4/5       50    0
      100 |   0.15824724051500    0.0000e+00 |  3.4748e-07  3.475e+00  3.475e-01 |  True  True    0/5       80    0
      200 |   0.15824724051500    0.0000e+00 |  3.4748e-07  3.475e+00  3.475e-01 |  True  True    0/5       80    0
     1000 |   0.15824724051500    0.0000e+00 |  3.4748e-07  3.475e+00  3.475e-01 |  True  True    0/5       80    0
     5000 |   0.15824724051500    0.0000e+00 |  3.4748e-07  3.475e+00  3.475e-01 |  True  True    0/5       80    0
```

The F1 fact reproduces exactly: `max_inner=3` gives `|dk| = 5.26 × keff_tol` with
`converged=True` and **zero** `ConvergenceWarning`s. And every row `max_inner ≤ 50` sits
in the converged-outer / truncated-inner state.

**The defect is NOT monotone in the keff error, in either direction:**

| `max_inner` | `\|dk\|/keff_tol` | `defect` | comment |
|---|---|---|---|
| 3 | **5.26** (worst) | 4.01e-06 | |
| 2 | 4.18 | 2.45e-06 | *lower* defect than `mi=1`, higher error |
| 1 | 1.59 | **4.09e-06** (worst) | *highest* defect, 3.3× *less* error than `mi=3` |
| 5 | 1.87 (still corrupting) | 4.88e-07 | 8× lower defect than `mi=1`, similar error |
| 10 | 0.294 (benign) | 8.80e-07 | **1.8× HIGHER defect than `mi=5`, 6× LOWER error** |
| 20 | 0.091 (benign) | 1.18e-07 | **3× LOWER defect than the fully-converged row** |
| ≥100 | 0 | 3.47e-07 | the honest converged answer |

Two orderings are inverted at once, and the most damaging inversion is that the
**truncated `mi=20` row (1.18e-07) has a *lower* defect than the fully-converged
reference (3.47e-07)** — so the defect does not even order truncated-vs-untruncated
correctly. The `mi=5` row is genuinely corrupting (1.87 × keff_tol) at `defect =
4.88e-07`, only **1.4×** above the honest converged row's 3.47e-07.

---

## 1 + 3. THE ANSWER — the two populations OVERLAP; the certificate is REFUTED

`scratch/n5_probe_06_population.py`, 21 rows over 5 geometries (leaking d1/d2, all-
reflective d1/d2/d3), three mixtures, `keff_tol = 1e-7` as the classification standard.

```
============================================================================================================================================================
case                 |dk|/keff_tol      class |     defect d/keff_tol      d_int      d_bnd  bnd_frac    balance |       GAIN   trunc/n  fully
============================================================================================================================================================
LEAK d2 mi=1            1.5820e+00    CORRUPT |  4.088e-06  4.088e+01  1.832e-06  3.655e-06   0.89392  1.802e-06 |  3.870e-02   58/58    False
LEAK d2 mi=2            4.1727e+00    CORRUPT |  2.446e-06  2.446e+01  9.609e-07  2.249e-06   0.91957  9.398e-07 |  1.706e-01   33/33    False
LEAK d2 mi=3            5.2480e+00    CORRUPT |  4.015e-06  4.015e+01  1.612e-06  3.677e-06   0.91588  1.576e-06 |  1.307e-01   24/24    False
LEAK d2 mi=5            1.8615e+00    CORRUPT |  4.883e-07  4.883e+00  1.952e-07  4.475e-07   0.91664  1.922e-07 |  3.813e-01   18/18    False
LEAK d2 mi=10           3.0431e-01     benign |  8.803e-07  8.803e+00  3.727e-07  7.975e-07   0.90593  3.629e-07 |  3.457e-02   11/11    False
LEAK d2 mi=20           1.0066e-01     benign |  1.184e-07  1.184e+00  4.784e-08  1.083e-07   0.91476  4.621e-08 |  8.500e-02    6/6     False
LEAK d2 mi=50           9.8801e-03     benign |  3.473e-07  3.473e+00  1.579e-07  3.093e-07   0.89072  1.541e-07 |  2.845e-03    4/5     False
LEAK d2 mi=100          9.8797e-03     benign |  3.475e-07  3.475e+00  1.579e-07  3.095e-07   0.89073  1.542e-07 |  2.843e-03    0/5      True
LEAK d2 mi=1000         9.8797e-03     benign |  3.475e-07  3.475e+00  1.579e-07  3.095e-07   0.89073  1.542e-07 |  2.843e-03    0/5      True
LEAK d1 mi=1            1.0397e+01    CORRUPT |  7.760e-07  7.760e+00  7.760e-07  2.506e-16   0.00000  7.775e-07 |  1.340e+00   39/39    False
LEAK d1 mi=3            1.1399e-01     benign |  9.186e-09  9.186e-02  9.186e-09  2.013e-16   0.00000  9.135e-09 |  1.241e+00   22/22    False
LEAK d1 mi=10           7.7384e-02     benign |  2.253e-08  2.253e-01  2.253e-08  1.829e-16   0.00000  4.176e-09 |  3.435e-01   11/11    False
LEAK d1 mi=50           1.1580e-01     benign |  2.520e-08  2.520e-01  2.520e-08  1.109e-16   0.00000  4.464e-09 |  4.595e-01    4/11    False
REFL d2 mi=3            3.5188e+00    CORRUPT |  1.934e-04  1.934e+03  1.200e-05  1.930e-04   0.99807  9.908e-06 |  1.820e-03   58/58    False
REFL d2 mi=20           3.9605e-02     benign |  7.511e-05  7.511e+02  2.565e-06  7.507e-05   0.99942  3.360e-07 |  5.273e-05   11/11    False
REFL d2 mi=200          1.0969e-06     benign |  6.108e-10  6.108e-03  2.443e-11  6.103e-10   0.99920  1.261e-11 |  1.796e-04    4/4     False
REFL d3 mi=20 GS        7.3832e-03     benign |  6.410e-05  6.410e+02  6.817e-07  6.409e-05   0.99994  6.215e-08 |  1.152e-05   17/17    False
REFL d3 mi=200 GS       1.0988e-04     benign |  4.065e-07  4.065e+00  3.911e-09  4.065e-07   0.99995  6.599e-10 |  2.703e-05    4/4     False
REFL d1 mi=5            7.2176e-01     benign |  1.278e-07  1.278e+00  7.447e-08  1.039e-07   0.81270  6.615e-08 |  5.647e-01   17/17    False
REFL d1 mi=50           5.2423e-04     benign |  9.056e-08  9.056e-01  5.958e-08  6.820e-08   0.75310  1.248e-08 |  5.789e-04    2/4     False
============================================================================================================================================================
CORRUPTING population: n=6  defect range [4.883e-07, 1.934e-04]
BENIGN      population: n=14  defect range [6.108e-10, 7.511e-05]

*** OVERLAP ***  the loosest corrupting defect is 4.883e-07, the tightest-admitting benign bound is 7.511e-05
    a threshold that ADMITS every benign case (T >= 7.511e-05) MISSES 5 of 6 corrupting cases
    a threshold that CATCHES every corrupting case (T < 4.883e-07) FALSE-ALARMS on 3 of 14 benign cases


INTERIOR-ONLY variant (drop the trace rows):
  CORRUPTING d_int range [1.952e-07, 1.200e-05]
  BENIGN     d_int range [2.443e-11, 2.565e-06]
  -> T >= 2.565e-06 (admits all benign) MISSES 5 of 6 corrupting
```

### The answer to item 1, up front

> The defect does **not** separate the poles. Measured on the two populations the ranges
> are `CORRUPTING [4.883e-07, 1.934e-04]` and `BENIGN [6.108e-10, 7.511e-05]` — a
> **154×-wide overlap** covering 5 of 6 corrupting cases and 3 of 14 benign ones. There is
> no threshold with any margin: **admitting every benign case costs 5 of 6 corrupting
> detections**, including the F1 case the certificate was proposed for. Dropping the trace
> rows (`d_int`) narrows the overlap 13× but still misses **5 of 6**.

### The structural reason — a stabiliser mismatch, not a tuning problem

The mechanism is measured, not inferred. Two facts:

**(a) The residual lives overwhelmingly in the reflective-trace rows, and the eigenvalue
cannot see them.** `bnd_frac = ‖r_boundary‖/‖r‖` reaches **0.99995** on the all-reflective
d=3 row. A defect in the reflective *inflow trace* of a zero-leakage system carries no
net current and changes neither production nor absorption to first order — so the
certificate's dominant component is, by construction, in the null space of
`k = production/absorption`.

**(b) The residual→eigenvalue transfer gain is unbounded across configurations.**
`GAIN = |Δk| / defect` sorted by `bnd_frac` (verbatim, same run):

```
  case                  bnd_frac        GAIN
  LEAK d1 mi=1           0.00000   1.340e+00
  LEAK d1 mi=50          0.00000   4.595e-01
  LEAK d1 mi=10          0.00000   3.435e-01
  LEAK d1 mi=3           0.00000   1.241e+00
  REFL d1 mi=50          0.75310   5.789e-04
  REFL d1 mi=5           0.81270   5.647e-01
  LEAK d2 mi=50          0.89072   2.845e-03
  LEAK d2 mi=100         0.89073   2.843e-03
  LEAK d2 mi=1000        0.89073   2.843e-03
  LEAK d2 mi=1           0.89392   3.870e-02
  LEAK d2 mi=10          0.90593   3.457e-02
  LEAK d2 mi=20          0.91476   8.500e-02
  LEAK d2 mi=3           0.91588   1.307e-01
  LEAK d2 mi=5           0.91664   3.813e-01
  LEAK d2 mi=2           0.91957   1.706e-01
  REFL d2 mi=3           0.99807   1.820e-03
  REFL d2 mi=200         0.99920   1.796e-04
  REFL d2 mi=20          0.99942   5.273e-05
  REFL d3 mi=20 GS       0.99994   1.152e-05
  REFL d3 mi=200 GS      0.99995   2.703e-05
GAIN  |dk| / defect  spans [1.152e-05, 1.340e+00]  = 1.163e+05x  -- the residual->eigenvalue transfer is NOT bounded across configs
```

The all-vacuum 1-D slab (`bnd_frac = 0`, no reflective trace to hide in) transfers
residual to eigenvalue error at **GAIN ≈ 1** — the residual is a faithful proxy there.
The all-reflective d=3 box transfers at **GAIN ≈ 1e-5**. Five decades apart, with a
monotone trend at the extremes. Since a threshold on `defect` bounds `|Δk|` only through
this gain, and the gain spans `1.16e+05×`, **no single constant can do the job** — this is
a `vv-principles` **Mode 12** situation read in the mirror: the *residual* functional and
the *eigenvalue* functional have different invariance groups, and neither contains the
other, so the residual is simultaneously **over**-sensitive (fires on reflective-trace lag
that keff is blind to) and **under**-sensitive (the `mi=20` row).

(Caveat stated honestly: `bnd_frac` is the dominant explanatory variable at the extremes
but not a sufficient predictor — within each geometry class the gain still spans ~2
decades, e.g. `REFL d1` reads `5.8e-04` and `5.6e-01` at two budgets.)

---

## 1 + 3 (widened) — 38 rows, 8 geometries, 3 mixtures, and 4 candidate statistics

`scratch/n5_probe_07_variants.py`. Geometries: leaking d1/d2/d2-thick/d2-all-vacuum/d3,
all-reflective d1/d2/d3. `L*` = leaking, `R*` = all-reflective; `A` = `get_mixture("A","2g")`,
`F50`/`F75` = the repaired `fissile(c)`. Four statistics:

| statistic | definition |
|---|---|
| `defect` | `‖r‖/‖q‖` — **the brief's proposal** |
| `d_int` | `‖r_interior‖/‖q‖` — drop the trace rows |
| `balance` | `‖R_g‖/‖Q_g‖`, `R_g = Σ_n w_n Σ_i V_i r[n,g,i]` — the angle- and volume-integrated per-group rate defect, i.e. the projection onto the functional `k = production/absorption` actually reads |
| `bal_adj` | `\|Σ_g φ†_g R_g\| / \|Σ_g φ†_g Q_g\|` with `φ†` the 0-D left eigenvector — first-order perturbation theory with a spatially FLAT adjoint |

```
0-D adjoint weights (control: k from the hand-built A^-1 F pencil MUST
equal solve_homogeneous_infinite's k_inf -- else the weight is junk):
  A    phi_dag = [0.111111 0.888889]   k_pencil = 1.875000000000001   k_inf = 1.875000000000001   |d| = 0.00e+00
  F50  phi_dag = [0.142857 0.857143]   k_pencil = 0.438461538461538   k_inf = 0.438461538461538   |d| = 0.00e+00
  F75  phi_dag = [0.142857 0.857143]   k_pencil = 1.168421052631579   k_inf = 1.168421052631579   |d| = 0.00e+00

building references ...
  L2   0.158247241502977   [tight self-ref]
  L2   0.453979966939393   [tight self-ref]
  L2t  0.968351816600070   [tight self-ref]
  L2v  0.055186776045086   [tight self-ref]
  L1   0.771447126977926   [tight self-ref]
  L1   0.460968760548238   [tight self-ref]
  L3   0.142186042934738   [tight self-ref]
  R2   1.875000000000001   [k_inf (0-D)]
  R1   1.168421052631579   [k_inf (0-D)]
  R3   0.438461538461538   [k_inf (0-D)]

================================================================================================================================
case                    |dk|/ktol    class |     defect      d_int    balance    bal_adj bnd_frac |     G_def     G_adj     tr/n
================================================================================================================================
L2/A mi=1               1.582e+00  CORRUPT |  4.088e-06  1.832e-06  1.802e-06  5.217e-06  0.89392 |  3.87e-02  3.03e-02  58/58  
L2/A mi=2               4.173e+00  CORRUPT |  2.446e-06  9.609e-07  9.398e-07  2.022e-08  0.91957 |  1.71e-01  2.06e+01  33/33  
L2/A mi=3               5.248e+00  CORRUPT |  4.015e-06  1.612e-06  1.576e-06  1.480e-06  0.91588 |  1.31e-01  3.55e-01  24/24  
L2/A mi=5               1.862e+00  CORRUPT |  4.883e-07  1.952e-07  1.922e-07  1.134e-06  0.91664 |  3.81e-01  1.64e-01  18/18  
L2/A mi=10              3.043e-01   benign |  8.803e-07  3.727e-07  3.629e-07  1.059e-06  0.90593 |  3.46e-02  2.87e-02  11/11  
L2/A mi=20              1.007e-01   benign |  1.184e-07  4.784e-08  4.621e-08  8.720e-08  0.91476 |  8.50e-02  1.15e-01   6/6   
L2/A mi=50              9.880e-03   benign |  3.473e-07  1.579e-07  1.541e-07  5.409e-07  0.89072 |  2.84e-03  1.83e-03   4/5   
L2/A mi=200             9.880e-03   benign |  3.475e-07  1.579e-07  1.542e-07  5.411e-07  0.89073 |  2.84e-03  1.83e-03   0/5   
L2/F75 mi=1             4.069e-01   benign |  4.284e-07  2.006e-07  1.885e-07  4.958e-07  0.88356 |  9.50e-02  8.21e-02  48/48  
L2/F75 mi=2             1.871e+00  CORRUPT |  5.334e-07  2.202e-07  2.039e-07  1.337e-07  0.91080 |  3.51e-01  1.40e+00  27/27  
L2/F75 mi=3             1.193e+00  CORRUPT |  1.185e-07  8.014e-08  7.824e-08  4.620e-07  0.73680 |  1.01e+00  2.58e-01  21/21  
L2/F75 mi=5             5.341e-01   benign |  1.845e-07  7.619e-08  6.996e-08  8.078e-08  0.91071 |  2.90e-01  6.61e-01  15/15  
L2/F75 mi=20            4.128e-03   benign |  1.430e-08  7.323e-09  6.802e-09  2.114e-08  0.85895 |  2.89e-02  1.95e-02   7/7   
L2t/A mi=1              2.127e+01  CORRUPT |  2.127e-05  6.664e-06  5.265e-06  5.193e-06  0.94967 |  1.00e-01  4.10e-01  77/77  
L2t/A mi=3              2.520e+00  CORRUPT |  2.148e-06  6.931e-07  5.645e-07  2.076e-07  0.94651 |  1.17e-01  1.21e+00  36/36  
L2t/A mi=10             1.007e+00  CORRUPT |  8.378e-07  2.713e-07  2.207e-07  7.467e-08  0.94614 |  1.20e-01  1.35e+00  17/17  
L2t/A mi=50             1.665e-02   benign |  4.023e-08  1.184e-08  8.268e-09  4.896e-08  0.95569 |  4.14e-02  3.40e-02   9/9   
L2v/A mi=1              3.482e+00  CORRUPT |  2.768e-06  2.768e-06  2.753e-06  9.071e-06  0.00000 |  1.26e-01  3.84e-02  29/29  
L2v/A mi=3              1.349e-01   benign |  1.075e-07  1.075e-07  1.066e-07  3.505e-07  0.00000 |  1.25e-01  3.85e-02  16/16  
L2v/A mi=10             3.994e-02   benign |  9.016e-08  9.016e-08  3.510e-08  6.349e-08  0.00000 |  4.43e-02  6.29e-02   7/7   
L1/F75 mi=1             1.040e+01  CORRUPT |  7.760e-07  7.760e-07  7.775e-07  4.489e-06  0.00000 |  1.34e+00  2.32e-01  39/39  
L1/F75 mi=2             1.930e+00  CORRUPT |  1.444e-07  1.444e-07  1.445e-07  8.338e-07  0.00000 |  1.34e+00  2.31e-01  26/26  
L1/F75 mi=3             1.140e-01   benign |  9.186e-09  9.186e-09  9.135e-09  5.206e-08  0.00000 |  1.24e+00  2.19e-01  22/22  
L1/F75 mi=10            7.738e-02   benign |  2.253e-08  2.253e-08  4.176e-09  4.557e-10  0.00000 |  3.44e-01  1.70e+01  11/11  
L1/F75 mi=50            1.158e-01   benign |  2.520e-08  2.520e-08  4.464e-09  1.199e-08  0.00000 |  4.59e-01  9.66e-01   4/11  
L1/A mi=1               2.370e+01  CORRUPT |  3.486e-06  3.486e-06  3.486e-06  2.397e-05  0.00000 |  6.80e-01  9.89e-02  57/57  
L1/A mi=2               2.899e+00  CORRUPT |  4.262e-07  4.262e-07  4.260e-07  2.929e-06  0.00000 |  6.80e-01  9.90e-02  38/38  
L1/A mi=5               2.997e+00  CORRUPT |  4.377e-07  4.377e-07  4.369e-07  3.010e-06  0.00000 |  6.85e-01  9.96e-02  20/20  
L3/F50 mi=5             1.593e-03   benign |  6.039e-07  4.295e-07  4.288e-10  9.372e-10  0.70309 |  2.64e-04  1.70e-01  36/36  
L3/F50 mi=20            3.621e-05   benign |  3.091e-07  1.755e-07  2.444e-08  3.982e-08  0.82306 |  1.17e-05  9.09e-05   6/6   
L3/F50 mi=200           6.752e-03   benign |  9.983e-08  5.569e-08  4.900e-08  8.278e-08  0.82995 |  6.76e-03  8.16e-03   0/5   
R2/A mi=3               3.519e+00  CORRUPT |  1.934e-04  1.200e-05  9.908e-06  6.424e-05  0.99807 |  1.82e-03  5.48e-03  58/58  
R2/A mi=20              3.960e-02   benign |  7.511e-05  2.565e-06  3.360e-07  2.607e-06  0.99942 |  5.27e-05  1.52e-03  11/11  
R2/A mi=200             1.097e-06   benign |  6.108e-10  2.443e-11  1.261e-11  8.425e-11  0.99920 |  1.80e-04  1.30e-03   4/4   
R1/F75 mi=5             7.218e-01   benign |  1.278e-07  7.447e-08  6.615e-08  1.412e-07  0.81270 |  5.65e-01  5.11e-01  17/17  
R1/F75 mi=50            5.242e-04   benign |  9.056e-08  5.958e-08  1.248e-08  3.635e-08  0.75310 |  5.79e-04  1.44e-03   2/4   
R3/F50 mi=20            7.383e-03   benign |  6.410e-05  6.817e-07  6.215e-08  1.330e-07  0.99994 |  1.15e-05  5.55e-03  17/17  
R3/F50 mi=200           1.099e-04   benign |  4.065e-07  3.911e-09  6.599e-10  1.873e-09  0.99995 |  2.70e-05  5.87e-03   4/4   
================================================================================================================================
```

### Separation, and its robustness to the classification standard

```
########## classification standard: ABSOLUTE |dk| > keff_tol (16 corrupting / 22 benign) ##########
  --- defect  (||r||/||q||)
      CORRUPTING n=16 range [1.185e-07, 1.934e-04]   BENIGN n=22 range [6.108e-10, 7.511e-05]
      OVERLAP width   633.69x   zero-FA T=7.511e-05 MISSES 15/16   BEST T=4.065e-07 catches 14/16 FA 5/22
  --- d_int   (bulk rows only)
      CORRUPTING n=16 range [8.014e-08, 1.200e-05]   BENIGN n=22 range [2.443e-11, 2.565e-06]
      OVERLAP width    32.00x   zero-FA T=2.565e-06 MISSES 12/16   BEST T=1.755e-07 catches 14/16 FA 5/22
  --- balance (angle+vol projection)
      CORRUPTING n=16 range [7.824e-08, 9.908e-06]   BENIGN n=22 range [1.261e-11, 3.629e-07]
      OVERLAP width     4.64x   zero-FA T=3.629e-07 MISSES 5/16   BEST T=1.885e-07 catches 14/16 FA 2/22
  --- bal_adj (0-D-adjoint-weighted)
      CORRUPTING n=16 range [2.022e-08, 6.424e-05]   BENIGN n=22 range [8.425e-11, 2.607e-06]
      OVERLAP width   128.95x   zero-FA T=2.607e-06 MISSES 8/16   BEST T=5.411e-07 catches 11/16 FA 2/22

########## classification standard: RELATIVE |dk|/k > keff_tol (19 corrupting / 19 benign) ##########
  --- defect  (||r||/||q||)
      CORRUPTING n=19 range [1.075e-07, 1.934e-04]   BENIGN n=19 range [6.108e-10, 7.511e-05]
      OVERLAP width   698.38x   zero-FA T=7.511e-05 MISSES 18/19   BEST T=4.065e-07 catches 15/19 FA 4/19
  --- d_int   (bulk rows only)
      CORRUPTING n=19 range [7.619e-08, 1.200e-05]   BENIGN n=19 range [2.443e-11, 2.565e-06]
      OVERLAP width    33.66x   zero-FA T=2.565e-06 MISSES 15/19   BEST T=7.447e-08 catches 19/19 FA 8/19
  --- balance (angle+vol projection)
      CORRUPTING n=19 range [6.996e-08, 9.908e-06]   BENIGN n=19 range [1.261e-11, 3.360e-07]
      OVERLAP width     4.80x   zero-FA T=3.360e-07 MISSES 7/19   BEST T=6.615e-08 catches 19/19 FA 4/19
  --- bal_adj (0-D-adjoint-weighted)
      CORRUPTING n=19 range [2.022e-08, 6.424e-05]   BENIGN n=19 range [8.425e-11, 2.607e-06]
      OVERLAP width   128.95x   zero-FA T=2.607e-06 MISSES 11/19   BEST T=1.330e-07 catches 16/19 FA 5/19

TRANSFER GAIN spread (max/min over the population):
  |dk| / defect  : [1.152e-05, 1.340e+00]  = 1.163e+05x
  |dk| / bal_adj : [9.094e-05, 2.064e+01]  = 2.270e+05x
```

The verdict does not depend on the classification standard: `634×` overlap absolute,
`698×` relative; zero-false-alarm threshold misses `15/16` and `18/19` respectively.

### ⛔ The adjoint-weighted variant is WORSE, and the reason generalises

`bal_adj` is the physically correct statistic (first-order perturbation theory) with the
*wrong weight*: the 0-D adjoint is spatially flat, which is exact only for the
all-reflective rows. Measured, it degrades the overlap from `4.64×` to `128.95×` and the
transfer-gain spread from `1.75e+04×` to `2.27e+05×` — **worse than the unweighted
`defect`**. The mechanism is visible in two rows: `L2/A mi=2` (CORRUPT) drops from
`balance = 9.398e-07` to `bal_adj = 2.022e-08`, a **46× collapse** giving the population's
worst outlier `G_adj = 20.6`; `L1/F75 mi=10` (benign) drops to `4.557e-10`, `G_adj = 17.0`.
**A SIGNED projection against a wrong weight manufactures accidental near-cancellations**,
so the cheap approximation is not merely imperfect — it is actively harmful. The positive
control matters here: `|k_pencil − k_inf| = 0.00e+00` for all three mixtures, so the weight
is not junk; it is *correct for the wrong problem*.

---

## 4. IS `_CERTIFICATE_SAFETY × tol` THE RIGHT BAR?  No — and it is the wrong QUANTITY

`scratch/n5_probe_08_safety_bar.py` (re-analysis of `scratch/n5_population.json`, no
re-solve). The outer's binding criterion is `dphi` in **every** solve measured — the
plan's own probe prints `dphi: 5.455e-07 vs tol 1.000e-06  met <- binding` — so the literal
transfer of the within-group pattern (`SAFETY × record.binding_criterion.tolerance`) is
`T = 10 × flux_tol = 1e-5`.

```
================================================================================================
standard: ABSOLUTE |dk| > keff_tol   (16 corrupting / 22 benign)
================================================================================================
                       threshold T |         caught   false alarms
------------------------------------------------------------------------------------------------
     SAFETY x keff_tol   = 1.0e-06 |      8/16           2/22     
     SAFETY x flux_tol   = 1.0e-05 |      2/16           2/22     
     1 x keff_tol        = 1.0e-07 |     16/16          13/22     
     1 x flux_tol        = 1.0e-06 |      8/16           2/22     
                           1.0e-08 |     16/16          20/22     
                           1.0e-09 |     16/16          21/22     

  at the literal transfer T = SAFETY x flux_tol = 1.0e-05:
    FLAGGED   (defect > T):
      CORRUPT  R2/A mi=3              defect=1.934e-04  |dk|/ktol=3.519e+00
      benign   R2/A mi=20             defect=7.511e-05  |dk|/ktol=3.960e-02
      benign   R3/F50 mi=20           defect=6.410e-05  |dk|/ktol=7.383e-03
      CORRUPT  L2t/A mi=1             defect=2.127e-05  |dk|/ktol=2.127e+01
    MISSED    (corrupting, defect <= T) — worst-error first:
      CORRUPT  L1/A mi=1              defect=3.486e-06  |dk|/ktol=2.370e+01
      CORRUPT  L1/F75 mi=1            defect=7.760e-07  |dk|/ktol=1.040e+01
      CORRUPT  L2/A mi=3              defect=4.015e-06  |dk|/ktol=5.248e+00
      CORRUPT  L2/A mi=2              defect=2.446e-06  |dk|/ktol=4.173e+00
      CORRUPT  L2v/A mi=1             defect=2.768e-06  |dk|/ktol=3.482e+00
      CORRUPT  L1/A mi=5              defect=4.377e-07  |dk|/ktol=2.997e+00
```

Three findings:

1. **The literal transfer catches 2 of 16 and false-alarms on 2 of 22.** Half of what it
   flags is benign, and it misses the row with the **largest keff error in the entire
   population** (`L1/A mi=1`, `|dk| = 23.7 × keff_tol`, defect only `3.486e-06`). Worse than
   a coin flip on the flagged set, and near-blind on the missed set.
2. **The constant is not the problem; the quantity is.** The whole trade-off curve is
   unacceptable: `T = 1e-7` reaches 16/16 sensitivity at a **59 % false-alarm rate**
   (13/22). There is no operating point where the certificate is a usable gate. So the
   brief's suspicion — "I do not assume the same constant transfers" — is right, and
   understated: **no constant transfers**, because the residual is a *flux-space* distance
   while the outer's contract is an *eigenvalue-space* one.
3. **The `dk`-vs-`dphi` ambiguity is real and matters 5×.** `SAFETY × keff_tol = 1e-6`
   catches 8/16; `SAFETY × flux_tol = 1e-5` catches 2/16. A certificate implemented by
   copying `_certify_within_group_exit`'s `record.binding_criterion.tolerance` would silently
   pick the *looser* of the two (`dphi` binds in every measured solve) and be the weaker
   gate — an increment tolerance standing in for a residual bound, mis-scaled by whichever
   criterion happened to bind.

**Floor, for the record (item 4's explicit ask).** For a genuinely converged solve at the
shipped defaults the measured defect is **3.47e-07 = 3.47 × keff_tol = 0.347 × flux_tol**
(§5). Any threshold must clear that. But clearing it is not the binding constraint — the
*benign truncated* population reaches **7.511e-05 = 751 × keff_tol**, which is 216× above
the honest floor and 18× above the worst-but-one corrupting case.

---

## RECOMMENDATION

**Do not build the outer certificate as `‖A ψ − Fφ/k‖/‖Fφ/k‖ > c·tol`.** It cannot gate
the `_warn_if_unconverged` widening from `converged` to `fully_converged`: on 38 measured
solves there is no threshold that admits the benign truncations and rejects the corrupting
ones, the overlap is 634×, and the literal transfer of the shipped constant catches 2 of 16
while missing the worst case in the set. Adopting it would produce exactly the outcome the
`solver.py:405-415` note is trying to avoid — warnings that are noise on 2 of 22 healthy
solves while staying silent on 14 of 16 genuinely wrong answers, i.e. **a filtered channel
that hides the next real truncation**.

The structural reason, in one line: **the residual and the eigenvalue are functionals with
different invariance groups.** `‖r‖` is dominated by the reflective-trace lag (up to
`bnd_frac = 0.99995`), which `k = production/absorption` is blind to by conservation; and
the residual→eigenvalue transfer gain spans `1.16e+05×` across configurations, so no
constant can convert one bound into the other.

What the measurements do support, in decreasing order of confidence:

1. ⭐ **Widen the guard to `fully_converged` and do NOT gate it on a residual.** The record
   already carries the truthful facts — `first_failure` names the starved inner, with its
   budget, its `ρ`, and `projected_iterations()`. That is *actionable* (`set max_inner=838`)
   and, unlike a defect threshold, it is never wrong: a truncated inner **is** a
   best-effort exit, and the plan's own F2 says the outer cannot detect it. The "noise"
   objection is about *volume* (27 further warnings), not about *falsity* — and the honest
   fix for volume is to make the message say how much it matters, not to suppress it.
   Note the population's own evidence for this: **`|dk|` exceeds `keff_tol` in 16 of 38
   rows** that all report `converged=True` with zero warnings today.
2. **If a residual number is wanted in the message, report `balance`, not `‖r‖`.** The
   angle-and-volume-integrated per-group rate defect is the projection onto the functional
   the eigenvalue reads. It narrows the overlap from `634×` to `4.64×`, and its best
   operating point (`T = 1.885e-07`) catches **14/16** corrupting at **2/22** false alarms —
   a genuinely informative diagnostic. Ship it as a *number in the warning*, not as a
   *gate*: `4.64×` residual overlap is still not a hard threshold.
3. **The only statistic that could gate is an adjoint-weighted residual with a REAL
   (spatially resolved) adjoint** — `δk/k ≈ ⟨ψ†, r⟩/⟨ψ†, Fφ/k⟩`. That costs a full adjoint
   solve per certified solve, which is almost certainly not worth it for a warning. And
   the shortcut is measured to be harmful: a flat 0-D adjoint makes things *worse*
   (`128.95×` overlap) via signed cancellation, so this must not be approximated.

**Owed to #340 regardless of the above:** `scratch/probe_340_done_when.py` carries an
internally inconsistent mixture (`σ_t` off by `±0.12`, `φ₁ ≡ 0`, effectively 1-group). Its
`keff`/`k_inf` numbers must not be cited, and the plan's benign-pole claim ("correct to
2.5e-11") should be re-measured on the repaired fixture (`sig_s=s.T`), where the SN reads
`0.43846153845055` against `k_inf = 0.43846153846154` — `|Δk| = 1.10e-11 = 1.10e-04 ×
keff_tol`, with all 4 inners truncated at 200/200 and `ρ ≈ 0.985`. That IS the benign pole
the plan wanted, and it exists; only the fixture was wrong.

---

## Promotion note

Two of these probes are regression-gate material (`tests/derivations/_promotion_policy.md`):

- **PROMOTE** `n5_probe_04_benign_fixture_audit.py` → a mixture **consistency gate**:
  `σ_t == σ_c + σ_f + Σ_to SigS[0][g,:]` for every fixture a verification test builds.
  This one defect silently invalidated a plan's whole benign pole, and the check is 3 lines.
  It has teeth right now (the shipped fixture reds it).
- **PROMOTE** the corrupting-sweep row as a **converged-but-wrong** regression:
  the F1 config at `max_inner=3` must report `history.fully_converged is False`. That is
  already true on HEAD and is the invariant the N6 commit-2 widening rides on.
- **LEAVE** `n5_probe_01/02/05/06/07/08` in `scratch/` — they are measurement instruments
  for this decision, print-only, and the policy forbids promoting print-only scripts.
