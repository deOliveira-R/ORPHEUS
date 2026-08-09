# #341 — why boundary Gauss-Seidel beats Jacobi at d=2 and loses at d=3

Investigation of record. Written incrementally; every number carries its configuration.
Branch `main`, HEAD `b0a003b4`, measured 2026-08-09 with `.venv/bin/python -O`.

Probes:

* `scratch/probe_341_iteration_spectrum.py` — the spectral instrument (below).
* `scratch/probe_gs_vs_jacobi_rate.py` — the pre-existing sweep-count probe (issue #341 comment).

---

## 0. TL;DR

1. **It is a genuine spectral inversion**, reproduced to 4 decimals by an instrument that
   eigen-solves `G = M⁻¹N` and never runs the solver's stopping test (§1).
2. **The mechanism (§10).** Zero leakage makes `B` the whole iteration, and `B`'s octant
   action is the hypercube `Q_d`. The diamond closure's per-cell face transmission is
   `Σ = (2/D)·1wᵀ − I`, whose spectrum is one absorption-damped eigenvalue plus **`d − 1`
   eigenvalues exactly `−1`** — an undamped, zero-cell-average sawtooth subspace whose
   dimension grows with `ndim`. That non-negativity failure **voids Varga's comparison
   theorem**, which would otherwise force `ρ_GS ≤ ρ_J` unconditionally (and does, exactly,
   in the non-negative model — §4). With the guarantee gone, *which* edges the fold takes
   decides the sign; and the shipped octant order takes the worst-balanced set available.
3. **At d=3 the two splittings are not racing the same mode** (G-S: the `z` faces;
   Jacobi: the `y` faces) — the "flip" is two different comparisons (§5).
4. **`ndim` is a proxy and a leaky one — do not branch on it.** Three d=2 fixtures where
   G-S LOSES (up to 2.86×) and three d=3 where it WINS (down to 0.58×) (§6, §11.1).
5. ⭐ **The dominant lever is a parameter nobody owns: the octant sweep ORDER.** All 25
   achievable d=3 fold patterns measured; they separate **exactly** on
   `max_a L_a > Σ_{b≠a} L_b` and span `n_GS/n_J ∈ [0.771, 1.892]`. The shipped order sits
   at 1.707 — 24th of 25 (§8).
6. **Recommendation: keep `"gauss_seidel"`, add no dimension branch, and fix two
   present-tense-false docstring claims** (`ρ_GS ≈ ρ_J²` and "the regular splitting")
   (§11.3). The corner where the default is the slower arm is narrow: `solve_sn` already
   defaults to Jacobi, and at `c ≥ 0.99` the inversion disappears entirely (§9, §11.2).
7. The issue's recorded story ("6 faces ⟹ a smaller fraction of the boundary coupling
   updated per pass") is **half right for a reason that is not face count**: the foldable
   fraction is `(2^d−1)/(d·2^d)` and shrinks from the schedule's *deferred-reflect* rule —
   but §2.2 gives the number showing the shrink is far too small, and §3 shows `ρ_GS` is
   not even monotone in it.

---

## 1. The instrument, and its positive control

### 1.1 What is measured

`SourceIteration.solve` (`orpheus/numerics/iteration.py:645`) runs

```
psi  <-  A_inv.apply( q_ext + sum_i g_i.apply(psi),  initial_guess=psi_prev )
```

so with `q_ext = 0` the iteration matrix is exactly `G = M⁻¹N` where the pair `(M, N)`
comes from the production selector `_select_si_splitting` (`orpheus/sn/solver.py:784`):

| schedule | `M` (implicit) | `N` (lagged gains) |
|---|---|---|
| `"jacobi"` | `L + C` | `S + B` |
| `"gauss_seidel"` | `(L + C) − B_lower` | `S + B_upper` |

`probe_341_iteration_spectrum.py` builds `M`, `N` through the **production** builders
(`build_within_group_system`, `_select_si_splitting`, `_maybe_window`), wraps
`x ↦ step.apply(Σ gᵢ.apply(x), initial_guess=x)` as a `scipy.sparse.linalg.LinearOperator`
over the composite's own `to_flat`/`from_flat`, and runs ARPACK `eigs(which="LM")`.

`[M]` linearity control on every case reported below: `max|G(2x) − 2G(x)| = 0.00e+00`
and `|G(0)| = 0.00e+00` — the seeded-apply `initial_guess` is inert on Cartesian, so `G`
is an honest linear operator.

### 1.2 Positive control (vv anti-pattern #17) — FIRST ROW, PASSED

Configuration: `axes` extents (1.0, 2.0, 3.0), cells (3, 4, 5), all six faces
`reflective`, pure absorber Σ_t = (0.8, 1.6) two groups (`sig_s = 0`, `sig_f = 0`),
`Quadrature.level_symmetric(sn_order=4)` (N = 24), DD closure.

| quantity | this instrument (ARPACK on `M⁻¹N`) | independently measured (`scratch/d3_absorber_diagnosis.md` §1.4/§7, fitted from the residual decay of a real solve) |
|---|---|---|
| ρ, boundary G-S | **0.9855216** | 0.985348 |
| ρ, Jacobi | **0.9754086** | 0.975014 |

Agreement to 4 decimal places on both arms, from two structurally different routes (an
eigen-solve of the operator vs a least-squares fit to a 1631-sweep residual history).
The instrument is verified; every other row below may be read.

Second control, d=2 — extents (1.0, 2.0), cells (3, 4), 4 faces reflective, same
mixture/quadrature:

| | ρ | predicted `n_GS/n_J = ln ρ_J / ln ρ_GS` | measured sweep ratio (#341 comment) |
|---|---|---|---|
| G-S | 0.9064250 | | |
| Jacobi | 0.9636488 | **0.377** | **0.398** (258 / 648) |

### 1.3 A structural by-product: `A = L + C − S − B` is SINGULAR on an all-reflective box

Both spectra carry eigenvalues **exactly** at `|λ| = 1` with **bulk fraction 0.0000**
(pure trace modes) — 3 of them at d=2, ≥5 at d=3. `λ = 1` for `G = M⁻¹N` means
`(M − N)v = Av = 0`: a genuine null vector of the discrete loss operator, and it is
splitting-*invariant* (it must be — it does not depend on how `A` is split).

This is the same object `scratch/d3_absorber_diagnosis.md` §4 found empirically as
"the converged state's trace retains an 11.26 % deviation from uniform while the exact
uniform field also has residual 1.06e-15 … the all-reflective 3-D DD operator has a
(near-)null space in the trace block". It is **exact**, not near.

Why the solve still converges: the SI stop is the residual `r_n = Σ gᵢ(ψ_{n−1} − ψ_n)`,
and the increment is `ψ_n − ψ_{n−1} = G^{n−1}(G − I)e_0`, which annihilates
`ker(G − I)`. The null component of the error is frozen at its initial value forever and
is invisible to the stopping test. With a zero cold start that component is
`−P_null ψ*`, i.e. the solver returns a state that differs from the intended uniform
field by a trace-only, zero-cell-average mode. **Benign here** (cell averages and the
scalar flux are uniquely determined), but it means the rate analysis must use
`ρ_eff = max{|λ| : λ ≠ 1}`, which is what §1.2 reports.

Jacobi additionally carries unit modes at `λ = −1` (3 at d=2, ≥5 at d=3). Those are
splitting-*dependent* (`(M + N)v = (L + C + B)v = 0`) and are the fingerprint of a
parity symmetry that G-S breaks — see §3.1.

---

## 2. What the boundary-G-S splitting actually folds

`SweepSchedule.gauss_seidel` (`orpheus/sn/loss_representation/sweep_schedule.py:113`) makes
**one group per in-plane octant, in quadrature order**, and reflects each reflective face
after the **LAST** octant group that outflows through it. The "last" (rather than "each")
is a correctness requirement for diagonal cubatures — reflecting after the first
outflowing group would reflect a not-yet-swept octant's seed and converge to the wrong
fixed point (the module's own ⚠ note, and ERR-056).

`[M]` the schedule as built, `level_symmetric(4)` (octant order is plain lexicographic
on `(sign_x, sign_y, sign_z)` with `−1` before `+1`):

```
d=2, 4 groups                          d=3, 8 groups
 g0 (-1,-1)  reflects: —                g0 (-,-,-) —      g4 (+,-,-) —
 g1 (-1,+1)  reflects: xmin             g1 (-,-,+) —      g5 (+,-,+) ymin
 g2 (+1,-1)  reflects: ymin             g2 (-,+,-) —      g6 (+,+,-) zmin
 g3 (+1,+1)  reflects: xmax, ymax       g3 (-,+,+) xmin   g7 (+,+,+) xmax,ymax,zmax
```

`[M]` the resulting split of the boundary coupling into implicit (`B_lower`) and lagged
(`B_upper`) inflow-ordinate rows:

| | xmin | xmax | ymin | ymax | zmin | zmax | total | fraction folded |
|---|---|---|---|---|---|---|---|---|
| d=2 | 12/12 | 0/12 | 6/12 | 0/12 | — | — | 18/48 | **0.3750** |
| d=3 | 12/12 | 0/12 | 6/12 | 0/12 | 3/12 | 0/12 | 21/72 | **0.2917** |

### 2.1 The closed form (derived, then measured)

Index octants by their sign vector `s ∈ {±1}^d` and let `σ(s)` be its position in the
sweep order. Face `f = (a, +)` is reflected at `T_f = max{σ(s) : s_a = +1}`; its inflow
rows are the octants with `s_a = −1`, and such a row is in `B_lower` iff `σ > T_f`.
Since exactly one of `T_{(a,+)}, T_{(a,-)}` is `2^d − 1` (whichever sign the LAST octant
carries on axis `a`), one of the two faces of every axis contributes **zero**, and the
other contributes

> **`L_a` = the length of the maximal SUFFIX of the sweep order over which the sign on
> axis `a` is constant.**

Hence

```
#B_lower(octant rows) = Σ_a L_a ,        #B(octant rows) = d · 2^d
```

and for the lexicographic order `L_a = 2^{d-1-a}`, giving

```
folded fraction  =  (2^d − 1) / (d · 2^d)    =   1/2 (d=1),  3/8 (d=2),  7/24 (d=3),  15/64 (d=4)
```

which matches the measured table exactly (18/48, 21/72). A counting bound
(`|A_ℓ| ≤ d − ⌈log₂ ℓ⌉` on the axes constant over the last `ℓ` octants) shows
`Σ_a L_a ≤ 2^d − 1`, so **the lexicographic order the quadrature already uses is
OPTIMAL for this metric** — no reordering can fold more.

An *ideal* octant-level Gauss-Seidel (reflect after **each** outflowing group) would fold
exactly half the `d·2^d` directed edges. So the deferred-reflect rule delivers

```
efficiency vs ideal G-S = (2^d − 1)/(d·2^{d−1}) = 1 (d=1), 3/4 (d=2), 7/12 (d=3), 15/32 (d=4)
```

i.e. **the schedule degrades like `2/d`.** That is a real, dimension-driven effect and it
is the sharpened version of the issue's "smaller fraction per pass" story — but the
driver is the *deferred-reflect correctness rule*, not the face count.

### 2.2 …and it is nowhere near enough to explain the inversion

For a 2-cyclic / consistently-ordered problem, full G-S gives `ρ_GS = ρ_J²`; a partial
fold `f ∈ [0, ½]` of the directed coupling interpolates as `ρ ≈ ρ_J^{1+2f}`
(`f = 0 → ρ_J`, `f = ½ → ρ_J²`). Evaluating that with the *measured* `ρ_J`:

| | `f` | `ρ_J` `[M]` | `ρ_J^{1+2f}` (predicted) | `ρ_GS` `[M]` | verdict |
|---|---|---|---|---|---|
| d=2 | 0.3750 | 0.9636488 | 0.9372 | **0.9064250** | G-S BETTER than the fold model |
| d=3 | 0.2917 | 0.9754086 | 0.9613 | **0.9855216** | G-S far WORSE; model says it should still beat Jacobi |

⛔ **The foldable-fraction story is REFUTED as the mechanism.** It predicts G-S beats
Jacobi at *every* dimension (any `f > 0` gives `ρ_J^{1+2f} < ρ_J`). The d=3 measurement
is on the wrong side of `ρ_J` entirely. Something makes the fold *actively harmful* at
d=3, and no amount of "less coupling folded" produces that.

---

## 3. Experiment C — the octant ORDER. Predictions written before measuring.

The fold is controlled by `Σ_a L_a` (§2.1), which is a property of the **octant sweep
order** — a quantity I can vary in-process by permuting the cached
`Quadrature.octants` tuple (no production change; the fixed point is splitting-invariant
and is re-verified per row).

Baseline for all four rows: d=3 extents (1,2,3), cells (3,4,5), all-reflective, pure
absorber Σ_t = (0.8, 1.6), `level_symmetric(4)`, DD.
`[M]` ρ_GS = 0.9855216, ρ_J = 0.9754086, `Σ L_a = 7`, `f = 7/24 = 0.2917`.

| # | order | `Σ L_a` | `f` | PREDICTED ρ_GS | reasoning |
|---|---|---|---|---|---|
| **C1** | reversed `[(+++) … (−−−)]` | 7 | 7/24 | **0.985522** (to ~6 digits) | mirror-image schedule on a mirror-symmetric box ⟹ identical spectrum. A *positive control on the monkeypatch*: if this row moves, the instrument is wrong. |
| **C2** | lexicographic with **z** most significant (fully-folded axis = z, not x) | 7 | 7/24 | **0.985522 if the mechanism is combinatorial**; something else if it matters WHICH axis is folded | the box is NOT permutation-symmetric (extents 1/2/3, cells 3/4/5) and the diagnosis says the slow mode lives on the `x`-dominant ordinates. Genuinely open. |
| **C3** | "antipodal-last" (last two octants differ on every axis) | 3 | 1/8 | **0.9797 ± 0.004** — strictly between ρ_J = 0.97541 (f=0) and 0.98552 (f=7/24), by linear interpolation in `f` | tests MONOTONICITY in the folded fraction. If ρ_GS(f=1/8) > 0.98552 the relation is non-monotone ⟹ the fold's *structure*, not its size, is what hurts. |
| **C4** | d=2 (1,2)/(3,4) antipodal-last | 2 | 1/4 | **0.9255 ± 0.010**, still below ρ_J(d2) = 0.9636 | same interpolation at d=2 (baseline ρ_GS = 0.906425, f = 3/8). |

### 3.1 Experiment C — MEASURED (`scratch/probe_341_octant_order.py`)

`[M]` 2026-08-09. Same fixtures as the prediction table. `Quadrature.octants` permuted
in-process (a `cached_property` slot on a `copy.copy` of the rule); `ρ` from the §1
instrument, `ρ_J` re-measured per row as a control (it must not move — Jacobi has no
schedule).

| # | order | `L_a` | `Σ L_a` | `f` | PREDICTED ρ_GS | **MEASURED ρ_GS** | ρ_J `[M]` | `n_GS/n_J` |
|---|---|---|---|---|---|---|---|---|
| — | d=3 lex, x-major (**baseline**) | [4, 2, 1] | 7 | 0.2917 | — | **0.985522** | 0.975409 | 1.707 |
| C1 | d=3 reversed | [4, 2, 1] | 7 | 0.2917 | 0.985522 | **0.985522** | 0.975409 | 1.707 |
| C2 | d=3 lex, **z**-major | [1, 2, 4] | 7 | 0.2917 | 0.985522 *if combinatorial* | **0.983455** | 0.975409 | 1.492 |
| C3 | d=3 **antipodal-last** | [1, 1, 1] | 3 | 0.1250 | 0.9797 ± 0.004 | **0.971741** | 0.975409 | **0.869** |
| — | d=2 lex (**baseline**) | [2, 1] | 3 | 0.3750 | — | **0.906425** | 0.963649 | 0.377 |
| C4 | d=2 **antipodal-last** | [1, 1] | 2 | 0.2500 | 0.9255 ± 0.010 | **0.934755** | 0.963649 | 0.549 |

* **C1 ✅ PASSED exactly** (0.985522 both) — the mirror-image schedule on a
  mirror-symmetric box gives the identical spectrum. The permutation machinery is
  verified before any other row is read.
* **C2 ⛔ moved** (0.985522 → 0.983455). The fold is NOT a pure function of `Σ L_a`:
  *which* axis carries the deep fold matters. Folding the long/coarse `z` axis is
  slightly better than folding the short/fine `x` axis — the opposite of the naive
  "fold where the slow mode lives" expectation.
* **C3 ⛔⛔ REFUTED, and in the most informative direction.** Cutting the folded
  fraction from 7/24 to 1/8 did not move ρ_GS *toward* Jacobi from above — it took it
  **past** Jacobi: 0.985522 → **0.971741 < ρ_J = 0.975409**. With that ordering,
  boundary G-S **WINS at d=3** (`n_GS/n_J = 0.869`).
* **C4** same direction as predicted (less fold ⟹ worse at d=2) but the interpolation
  number was wrong by 4σ of its own band — another sign the fraction is not the variable.

### 3.2 What C1–C4 establish

`ρ_GS` is **non-monotone in the folded fraction**, and the sign of the monotonicity
**flips with dimension**:

```
d=2:  f = 1/4  -> 0.934755        f = 3/8   -> 0.906425     (more fold = better)
d=3:  f = 1/8  -> 0.971741        f = 7/24  -> 0.985522     (more fold = WORSE)
```

So the mechanism is **not "how much of B is folded" — it is WHICH edges, i.e. the
structure the fold induces on the octant-coupling graph.** And the d=3 inversion is a
property of *this particular schedule's* fold pattern, not of boundary Gauss-Seidel per se:
a different octant order removes it entirely.

---

## 4. The tractable model — and the two hypotheses it kills

`scratch/probe_341_octant_model.py`. On an all-reflective box the boundary coupling is a
walk on the **octant hypercube** `Q_d`: specular reflection at face `(a, ±)` flips the
`a`-th direction cosine, so octant `s` feeds octant `s ⊕ a`, and `Q_d` (2^d vertices,
degree `d`) is the *entire* coupling graph. Collapse each octant to one scalar amplitude
and give axis `a` one round-trip gain `g_a`:

```
Jacobi   K   = Σ_a g_a T_a        (T_a = flip-bit-a permutation)
G-S      G   = (I − K_l)⁻¹ K_u    (K_l = the exact ORPHEUS deferred-reflect fold)
```

Two theorems fall out, and both are confirmed numerically:

**(i) `ρ_J` is exactly `Σ_a |g_a|`, and BOTH `ρ_J` and `ρ_GS` are blind to the SIGNS of
`g_a`.** The `T_a` commute, so the characters of `(Z₂)^d` diagonalise `K`:
`spec(K) = {Σ_a g_a ε_a}`. And a sign flip `g_a → −g_a` is a *diagonal similarity*
(conjugate by `D[s,s] = Π_{a∈flip} s_a`, which sends `T_a → ε_a T_a` and preserves the
edge set, hence preserves `K_l` too). `[M]` §A of the probe: all 8 sign patterns of
`|g| = (0.45, 0.30, 0.20)` give `ρ_J = 0.95000` and `ρ_GS = 0.92542`, identically.

> ⛔ **REFUTED — "the DD closure's NEGATIVE transmission coefficients cause the
> inversion".** This was my leading hypothesis (the multi-D DD self-transmission
> `coeff_a = (w_a − s − Σ_{b≠a} w_b)/(s + Σ_b w_b)` is negative for at least `d−1` axes
> per ordinate, and Stein–Rosenberg/Varga forbids an inversion for a *non-negative*
> splitting, so sign indefiniteness is *necessary*). It is necessary but the sign
> **pattern across axes** cannot be the cause: on the octant graph a per-axis sign is a
> gauge.

**(ii) With all `g_a > 0` the two splittings are regular splittings of `I − K` with
`N_GS = K_u ≤ K = N_J` elementwise, so Varga's comparison theorem forces
`ρ_GS ≤ ρ_J` at every dimension and every order.** `[M]` §D: 0 violations in 1200
random positive-gain cases (d = 2, 3, 4 × both orders).

**(iii) …but the model NEVER reproduces the inversion, at any dimension or sign
pattern.** `[M]` §E: `P(ρ_GS > ρ_J)` = 0.00 for every `d ∈ {2,3,4}` and every count of
negative gains, over 400 random magnitude draws × all `2^d` sign patterns.

> ⛔ **REFUTED — "the inversion is a property of the octant-graph combinatorics + the
> fold pattern".** The 2^d-scalar model contains the full hypercube, the exact
> deferred-reflect fold, and the exact ordering dependence — and it cannot produce a
> loss. Whatever inverts the comparison lives in the structure the model throws away:
> the **intra-octant face-to-face coupling**, i.e. the fact that the "gain" is a `d × d`
> operator over the octant's `d` upwind/downwind faces, not a scalar.

### 4.1 …which is exactly where `ndim` enters, and it is NOT the face count

Do the intra-octant block honestly. For one ordinate on one homogeneous DD cell
(`w_a = 2|μ_a| A_a`, `s = Σ_t V`, `D = s + Σ_b w_b`), the multi-D diamond closure
`ψ_c = Σ_b w_b ψ_in,b / D`, `ψ_out,a = 2ψ_c − ψ_in,a` gives the face-to-face transmission

```
        Σ_{a←b} = (2/D) w_b − δ_ab        i.e.   Σ = (2/D) · 1 wᵀ − I
```

a **rank-one matrix minus the identity**, whose spectrum is immediate:

| eigenvalue | eigenvector | multiplicity | meaning |
|---|---|---|---|
| `1 − 2s/D = 1 − 2Σ_t V/D` | `1` (all faces equal) | 1 | the physical, **absorption-damped** mode |
| **`−1`** | `{v : wᵀv = 0}` | **`d − 1`** | `ψ_c = 0` ⟹ `ψ_out,a = −ψ_in,a`: an **UNDAMPED sawtooth**, invisible to `Σ_t V ψ_c` |

That is the object `scratch/d3_absorber_diagnosis.md` §4 measured as the face sawtooth
with trace ratios `1.074414 / 0.925586` **summing to exactly 2** — `ψ_out = 2ψ_c − ψ_in`
with `ψ_c = ψ*`, i.e. the `wᵀv = 0` eigenvector riding on the uniform solution.

**The multiplicity of the undamped subspace is `d − 1`.** d=1: none (a 1-D DD cell has no
zero-average trace mode). d=2: one. d=3: two. The same closure under STEP differencing
would give `Σ_step = (1/D')·1·w'ᵀ` — rank one, eigenvalues `{(D'−s)/D', 0×(d−1)}` — i.e.
the same `d−1` modes **maximally damped instead of undamped.** So the undamped subspace
is a property of the DIAMOND closure, and its dimension is a property of `ndim`.

---

## 5. Are the two splittings even racing the SAME mode? — NO at d=3

`scratch/probe_341_mode_shape.py`. Dominant (|λ| < 1) eigenvector of each iteration
matrix, same fixtures as §1. "per-face mass" is the fraction of the eigenvector's trace
norm carried by each boundary face.

| fixture / splitting | ρ | arg λ | bulk mass | where the trace mass lives | dominant ordinate class |
|---|---|---|---|---|---|
| d=2 (3,4) **G-S** | 0.906425 | +73.8° | 0.070 | **ymin 0.672, ymax 0.663**, xmax 0.236, xmin 0.229 | (0.869, 0.350, 0.350) |
| d=2 (3,4) **Jacobi** | 0.963649 | +131.0° | 0.059 | **ymin 0.671, ymax 0.671**, xmax/xmin 0.223 | (0.869, 0.350, 0.350) |
| d=3 (3,4,5) **G-S** | 0.985522 | −153.5° | 0.075 | **zmin 0.673, zmax 0.672**, xmin/xmax 0.18, ymin/ymax 0.12 | (0.869, 0.350, 0.350) |
| d=3 (3,4,5) **Jacobi** | 0.975409 | +141.4° | 0.110 | **ymin 0.658, ymax 0.658**, zmax/zmin 0.203, xmin/xmax 0.160 | (0.869, 0.350, 0.350) |

* At **d=2 both splittings race the SAME mode** — the `y`-face sawtooth on the
  `|μ_x| = 0.869` ordinates (spatial pattern `[-0.585, 1.000, -0.730]`, mean −0.105:
  sign-alternating, near-zero mean — the `wᵀv = 0` DD mode of §4.1). G-S accelerates it
  from 0.9636 to 0.9064.
* At **d=3 they do NOT.** Jacobi's slow mode sits on the **y** faces; boundary G-S's sits
  on the **z** faces. And `z` is precisely the axis the shipped fold barely touches
  (`L_z = 1`: 3 of 12 `zmin` inflow rows fresh, 0 of 12 on `zmax`), while `x` is fully
  folded (`L_x = 4`) and `y` half (`L_y = 2`).

⟹ **the "flip" is two different comparisons.** The shipped G-S fold does accelerate the
axes it folds deeply; what it leaves behind at d=3 is the **z** mode, and that mode is
*worse under G-S (0.9855) than the whole Jacobi spectrum is (0.9754)* — so the fold not
only fails to help the surviving mode, it degrades it. (The `(I − K_l)⁻¹ = Σ_n K_l^n`
series at the shipped d=3 order has nilpotency index **4** — the fresh chain
`0→4→6→7` spans three reflections in one pass — so `G_GS` carries `K_l K_u`, `K_l² K_u`,
`K_l³ K_u` terms that have no counterpart in `G_J`. At d=2 the chain is `0→2→3`,
index 3; under the antipodal order the fold is a *star* into the last octant, index 2.)

---

## 6. Is `ndim` the right discriminating variable? — **NO.** Both falsifiers fired.

`scratch/probe_341_predicate_sweep.py`, `[M]` 2026-08-09, 244 s. All fixtures
all-reflective (zero leakage), pure absorber, DD, shipped lex octant order. One group
(`Σ_t` as stated); **control**: the 2-group and 1-group d=3 (3,4,5) rows give
*bit-identical* ρ (0.985522 / 0.975409) — a pure absorber decouples the groups exactly.

### 6.1 d=2 fixtures where boundary G-S **LOSES**

| fixture | ρ_GS | ρ_J | `n_GS/n_J` |
|---|---|---|---|
| d=2, cells (1,1), extents (1,2), Σ_t = 0.8, LS4 | 0.868281 | 0.667791 | **2.859** |
| d=2, cells (2,2), extents (6,6), Σ_t = 2.0, LS4 | 0.672496 | 0.634620 | **1.146** |
| d=2, cells (3,3), extents (1,1), Σ_t = 51.2, LS4 | 0.791918 | 0.784687 | **1.039** |

### 6.2 d=3 fixtures where boundary G-S **WINS**

| fixture | ρ_GS | ρ_J | `n_GS/n_J` |
|---|---|---|---|
| d=3, cells (2,4,8), extents (1,2,4), Σ_t = 0.8, LS4 | 0.986278 | 0.992072 | **0.576** |
| d=3, cells (3,3,3), extents (1,1,1), Σ_t = 0.05, LS4 | 0.998004 | 0.998428 | **0.787** |
| d=3, cells (3,4,5), extents (1,2,3), Σ_t = 0.8, **LS2** | 0.968617 | 0.971130 | **0.919** |

### 6.3 The full sweep — the effect is CONTINUOUS in optical thickness, not in `ndim`

Cube, extent 1 per axis, 3 cells per axis, LS4, all-reflective absorber; only `Σ_t` varies:

| Σ_t | d=2 `n_GS/n_J` | d=3 `n_GS/n_J` |
|---|---|---|
| 0.05 | 0.308 | **0.787** |
| 0.2 | 0.318 | **1.055** |
| 0.8 | 0.392 | 1.587 |
| 3.2 | 0.544 | 1.831 |
| 12.8 | 0.926 | 1.938 |
| 51.2 | **1.039** | — |

The **same monotone curve** at both dimensions. `ndim` does not change the shape; it
**shifts the crossing by ~2.5 decades in Σ_t** (d=2 crosses near Σ_t ≈ 40, d=3 near
Σ_t ≈ 0.15). Other knobs move it just as decisively at fixed `ndim`: mesh refinement
(d=2 (1,1) → 2.859, (2,2) → 0.489, (10,10) → 0.158), aspect ratio (d=3 (5,5,5) → 1.558
vs (2,4,8) → 0.576), and quadrature order (d=3 LS2 → 0.919, LS4 → 1.707, LS6 → 1.760).

⟹ **`ndim` is a proxy that happens to correlate on the shipped verification fixtures.**
The honest predicate would have to be a joint function of the octant sweep order's fold
pattern `(L_a)`, the per-cell DD transmission spectrum
(`Σ = (2/D)·1wᵀ − I`, i.e. `Σ_t V` vs `2|μ_a|A_a` per ordinate — cell aspect ratio and
optical thickness), and the quadrature. There is **no cheap closed-form predicate**, and
`ndim` is not even a conservative one — it mislabels a 3-D anisotropic box and a
1-cell 2-D box in opposite directions. **A production default must not branch on it.**

---

## 7. ⭐ The falsifiable prediction — derived, written down, then tested

### 7.1 The mechanism it comes from

Collecting §3–§6: the boundary-G-S iteration matrix is

```
G_GS = (I − K_l)⁻¹ K_u = K_u + K_l K_u + K_l² K_u + …      (K_l + K_u = K = the Jacobi coupling)
```

`K_l` is strictly lower-triangular in the octant sweep order, so the series terminates at
`K_l^{p−1}` where `p` = the nilpotency index = **1 + (longest chain of FRESH octant→octant
edges the order admits in one pass)**. The shipped lexicographic order maximises that
chain (`0→4→6→7` at d=3 — one reflection per axis, in axis order), giving `p = d + 1`.
The antipodal-last order makes the fold a **star** into the final octant, giving `p = 2`.

The first-order term `K_l K_u` is the classical Gauss-Seidel gain — it is what turns a
2-cycle's `±√γ` into `γ`. The **higher-order terms `K_l^n K_u`, `n ≥ 2`, have no Jacobi
counterpart at all**, and §4 proved they cannot be signed away (a per-axis sign is a
gauge) nor produced by the octant combinatorics alone — they need the intra-octant `d × d`
DD block `Σ = (2/D)·1wᵀ − I`, whose `(d−1)`-dimensional `wᵀv = 0` subspace transmits with
factor exactly **−1** (undamped). Those are the terms that can amplify.

**Claim.** *The amplification that inverts the comparison is carried by the `n ≥ 2` terms.
An octant order with `K_l² = 0` therefore cannot invert it.*

### 7.2 PREDICTION D — written 2026-08-09 BEFORE the measurement

> Re-running the entire §6 sweep with the **antipodal-last** octant order (fold pattern
> `L_a = (1,…,1)`, `Σ L_a = d`, nilpotency index 2) will give
> **`ρ_GS ≤ ρ_J` on EVERY fixture** — in particular on the seven where the shipped
> lexicographic order loses:

| fixture | lex ρ_GS `[M]` | ρ_J `[M]` | PREDICTED antipodal ρ_GS |
|---|---|---|---|
| d=2 (1,1) ext(1,2) Σ_t 0.8 LS4 | 0.868281 | 0.667791 | ≤ 0.667791 |
| d=2 (2,2) ext(6,6) Σ_t 2.0 LS4 | 0.672496 | 0.634620 | ≤ 0.634620 |
| d=2 (3,3) ext(1,1) Σ_t 51.2 LS4 | 0.791918 | 0.784687 | ≤ 0.784687 |
| d=3 (2,2,2) ext(1,1,1) Σ_t 0.8 LS4 | 0.965026 | 0.942698 | ≤ 0.942698 |
| d=3 (3,3,3) ext(1,1,1) Σ_t 0.8 LS4 | 0.984267 | 0.975146 | ≤ 0.975146 (point estimate 0.972 ± 0.006) |
| d=3 (3,3,3) ext(1,1,1) Σ_t 12.8 LS4 | 0.827301 | 0.692550 | ≤ 0.692550 |
| d=3 (5,5,5) ext(1,1,1) Σ_t 0.8 LS4 | 0.994339 | 0.991192 | ≤ 0.991192 |

> and a **cost** prediction on the fixtures where lex already wins: the antipodal order is
> STRICTLY WORSE than lex there (fewer edges folded ⟹ less first-order G-S gain), e.g.
> d=2 (3,4): lex 0.906425 → antipodal 0.934755 (`[M]` already, §3.1) and
> d=2 (10,10): lex 0.979267 → antipodal in (0.979267, 0.996686).
>
> ⚠ Honest scope: the d=3 (3,4,5) LS4 row (lex 0.985522 → antipodal 0.971741 ≤ ρ_J
> 0.975409) is the observation the claim was BUILT on. It is not a test of it. The seven
> rows above are.
>
> **If any of the seven comes back with `ρ_GS > ρ_J`, Claim 7.1 is refuted** and the
> amplification is not a property of the fold's chain depth.

### 7.3 PREDICTION D — MEASURED. **Partly held, partly REFUTED.**

`scratch/probe_341_antipodal_prediction.py`, `[M]` 2026-08-09.

| fixture | lex ρ_GS | **anti ρ_GS** | ρ_J | predicted | verdict |
|---|---|---|---|---|---|
| d=2 (1,1) ext(1,2) Σ_t 0.8 | 0.868281 (2.859) | **0.647761** (0.930) | 0.667791 | ≤ 0.667791 | ✅ HELD |
| d=2 (2,2) ext(6,6) Σ_t 2.0 | 0.672496 (1.146) | **0.785421** (1.883) | 0.634620 | ≤ 0.634620 | ⛔ **REFUTED — and anti is WORSE than lex** |
| d=2 (3,3) ext(1,1) Σ_t 51.2 | 0.791918 (1.039) | **0.883306** (1.954) | 0.784687 | ≤ 0.784687 | ⛔ **REFUTED — anti WORSE than lex** |
| d=3 (2,2,2) Σ_t 0.8 | 0.965026 (1.658) | **0.935287** (0.882) | 0.942698 | ≤ 0.942698 | ✅ HELD |
| d=3 (3,3,3) Σ_t 0.8 | 0.984267 (1.587) | **0.971294** (0.864) | 0.975146 | ≤ 0.975146, est 0.972 ± 0.006 | ✅ HELD, point estimate inside band |
| d=3 (3,3,3) Σ_t 12.8 | 0.827301 (1.938) | **0.682125** (0.960) | 0.692550 | ≤ 0.692550 | ✅ HELD |
| d=3 (5,5,5) Σ_t 0.8 | 0.994339 (1.558) | **0.989405** (0.831) | 0.991192 | ≤ 0.991192 | ✅ HELD |
| d=3 (3,4,5) **LS6** Σ_t 0.8 | 0.989079 (1.760) | **0.976798** (0.823) | 0.980862 | (unlisted bonus row) | ✅ held |

Cost predictions, on fixtures where lex already wins — all HELD (antipodal folds less, so it
gains less):

| fixture | lex ρ_GS | anti ρ_GS | ρ_J | predicted |
|---|---|---|---|---|
| d=2 (3,4) | 0.906425 (0.377) | 0.934755 (0.549) | 0.963649 | worse than lex, still < ρ_J ✅ |
| d=2 (10,10) | 0.979267 (0.158) | 0.990987 (0.367) | 0.996686 | in (0.979267, 0.996686) ✅ |
| d=3 (2,4,8) | 0.986278 (0.576) | 0.991074 (0.888) | 0.992072 | ≤ ρ_J ✅ |

**Reading.** 6 of 8 tested rows held, including the quantitative point estimate. The two
refutations are both **optically thick d=2** cells, and they fail in the strongest way:
the antipodal order is not merely still-losing, it is *worse than the shipped order*
(0.785 vs 0.672; 0.883 vs 0.792).

⟹ **Claim 7.1 is refuted as a universal law.** There are (at least) **two** independent
amplification mechanisms, and the antipodal order fixes one and aggravates the other:

* **(a) fold-chain amplification** — the `K_l^n K_u`, `n ≥ 2` terms; dominant at
  moderate optical thickness and *increasing with `d`* because the lexicographic order's
  chain is one reflection per axis. This is the one behind the reported #341 d=3
  inversion, and cutting the chain (antipodal, or better §8) removes it.
* **(b) thick-cell amplification** — when `Σ_t V > ½D` the *smooth* DD channel
  `1 − 2Σ_tV/D` also goes negative, so ALL `d` transmission eigenvalues sit near `−1`
  and the sweep block is `Σ ≈ −I`: a signed permutation network with almost no damping.
  There the regular-splitting comparison has no content at all, and a shallower fold is
  *worse*. This is why the two d=2 thick rows invert and why the antipodal order makes
  them worse.

---

## 8. The complete d=3 fold-pattern enumeration — an EXACT separating predicate

`scratch/probe_341_fold_pattern.py d3`. All `8! = 40320` octant orders collapse to **25
achievable fold patterns** `(L_x, L_y, L_z)`; ρ_GS measured for a representative of each,
d=3 extents (1,2,3) cells (3,4,5) all-reflective 2-group absorber Σ_t = (0.8, 1.6) LS4 DD,
ρ_J = 0.975409 throughout (Jacobi has no schedule — a control that held on all 25).

| `Σ L_a` | pattern → `ρ_GS` (`n_GS/n_J`) |
|---|---|
| 3 | (1,1,1) 0.971741 (0.869) |
| 4 | (2,1,1) 0.969556 (0.805) · (1,2,1) 0.970227 (0.824) · (1,1,2) 0.971045 (0.847) |
| 5 | (2,2,1) 0.970145 (0.821) · (2,1,2) 0.970928 (0.844) · (1,2,2) 0.971022 (0.847) · **(1,1,3) 0.977538 (1.096)** · **(1,3,1) 0.979158 (1.182)** · **(3,1,1) 0.983825 (1.527)** |
| 6 | **(3,2,1) 0.968223 (0.771) ← best of all 25** · (3,1,2) 0.968225 (0.771) · (2,3,1) 0.969110 (0.794) · (1,3,2) 0.969212 (0.796) · (2,1,3) 0.969966 (0.816) · (1,2,3) 0.969992 (0.817) · **(1,4,1) 0.983846 (1.529)** · **(1,1,4) 0.985208 (1.671)** · **(4,1,1) 0.986929 (1.892) ← worst of all 25** |
| 7 | **(2,4,1) 0.981889 (1.362)** · **(1,4,2) 0.982632 (1.421)** · **(1,2,4) 0.983455 (1.492)** · **(2,1,4) 0.983569 (1.503)** · **(4,2,1) 0.985522 (1.707) ← SHIPPED** · **(4,1,2) 0.986036 (1.771)** |

Bold = G-S loses. The 25 rows separate **exactly**, with zero exceptions, on

> ### boundary G-S LOSES ⟺ `max_a L_a > Σ_{b≠a} L_b`
> *(the fold is dominated by one axis: its depth exceeds the others' combined)*

Check: (3,1,1) `3 > 2` loses; (3,2,1) `3 > 3` false, wins — and is the BEST of all 25.
(4,2,1) `4 > 3` loses. (2,2,1) `2 > 3` false, wins. **25/25.**

The shipped lexicographic order lands on `(4, 2, 1)` — the maximally dominated pattern
available, and the second-worst of the 25. Its rate 0.985522 vs the achievable best
0.968223 is a **2.2× difference in sweep count** (`1.707` vs `0.771` against Jacobi),
purchased by nothing but the order in which `Quadrature.octants` happens to enumerate the
sign octants.

⚠ **Scope of the law.** It is exact over all 25 patterns of ONE d=3 fixture. It does NOT
transfer to d=2: the shipped `(2,1)` is "dominated" (`2 > 1`) and is the BEST of the three
achievable d=2 patterns (0.906425 vs (1,2) 0.910279 vs (1,1) 0.934755). So *fold balance
is the right coordinate, and its optimal value is regime-dependent* — the same conclusion
as §6. There is no dimension-free predicate.

---

## 9. The scattering ratio — the inversion DISAPPEARS as `c → 1`

`scratch/probe_341_scattering_ratio.py`. All-reflective, Σ_t = 0.8, 1 group, LS4, DD,
shipped lex order; only the within-group scattering ratio `c = Σ_s/Σ_t` varies.

| `c` | d=2 (3,4) ρ_GS / ρ_J → `n_GS/n_J` | d=3 (3,4,5) ρ_GS / ρ_J → `n_GS/n_J` |
|---|---|---|
| 0.0 | 0.906425 / 0.963649 → 0.377 | 0.985522 / 0.975409 → **1.707** |
| 0.5 | 0.908040 / 0.963822 → 0.382 | 0.985084 / 0.975462 → **1.653** |
| 0.9 | 0.936010 / 0.963953 → 0.555 | 0.984741 / 0.975508 → **1.613** |
| 0.99 | 0.993629 / 0.994947 → 0.793 | 0.994317 / 0.995448 → **0.800** |
| 0.999 | 0.999363 / 0.999495 → 0.793 | 0.999432 / 0.999545 → **0.801** |

Near-critical (`c ≥ 0.99`) the SI spectrum is dominated by the *scattering* mode
(`ρ → c`), which both splittings share; the boundary mode is no longer the slow one and
G-S wins at BOTH dimensions. **The d=3 inversion lives only at `c ≲ 0.9` with zero
leakage** — a strongly-absorbing all-reflective box.

---

## 10. THE MECHANISM, stated structurally

1. **Zero leakage makes `B` the entire iteration, and `B`'s octant action is the hypercube
   `Q_d`.** Specular reflection at face `(a, ±)` flips direction cosine `a`, so octant `s`
   couples only to `s ⊕ a`: the coupling graph is `Q_d` (2^d vertices, degree `d`),
   bipartite by sign-parity. That parity is why the Jacobi spectrum is symmetric under
   `λ → −λ` (measured: unit modes at both `+1` and `−1` for Jacobi, only at `+1` for G-S,
   §1.3) and it is the symmetry boundary-G-S breaks.

2. **The intra-octant block is where `ndim` actually enters.** For one ordinate on one
   homogeneous DD cell the face-to-face transmission is `Σ = (2/D)·1wᵀ − I`
   (`w_a = 2|μ_a|A_a`, `D = Σ_t V + Σ_b w_b`): **one** absorption-damped eigenvalue
   `1 − 2Σ_tV/D` and **`d − 1` eigenvalues exactly `−1`** — the `wᵀv = 0` zero-cell-average
   sawtooth, which `Σ_t V ψ_c` cannot see. So a `d`-dimensional DD cell carries a
   `(d−1)`-dimensional **undamped** transmission subspace. The d=3 diagnosis's
   `1.074414 / 0.925586` face sawtooth (summing to exactly 2) is one of these.

3. **The undamped subspace destroys the comparison theorem.** Varga's comparison
   (`N_GS = S + B_upper ≤ S + B = N_J` elementwise on a *regular* splitting) would force
   `ρ_GS ≤ ρ_J` at every dimension and every order — and the 2^d scalar model, which IS
   non-negative, obeys it perfectly (0 violations in 1200 draws, §4). The DD closure is
   not non-negative: at least `d−1` transmission channels are exactly `−1`, and above
   `Σ_t V > ½D` the smooth channel joins them. Once the guarantee is gone, **which** edges
   the fold takes decides the sign, and no monotone law in "how much is folded" survives
   (§3.2 measures the monotonicity reversing between d=2 and d=3).

4. **The shipped fold is the worst-balanced one available.** `SweepSchedule.gauss_seidel`
   defers each face's reflect to its LAST outflowing octant group (an ERR-056 correctness
   requirement), which makes the fold depth on axis `a` equal to `L_a` = the constant-sign
   suffix run of the octant order (§2.1). Inheriting `Quadrature.octants`' lexicographic
   enumeration gives `L = (2^{d−1}, …, 2, 1)` — at d=3, `(4, 2, 1)`, which is exactly the
   `max_a L_a > Σ_{b≠a} L_b` "one-axis-dominated" pattern that §8 separates as the losing
   class, and the second-worst of all 25 achievable patterns.

5. **The consequence, measured: at d=3 the two splittings race DIFFERENT modes.** The deep
   folds on `x` and `y` do accelerate their own modes; what survives is the `z`-face mode
   (the axis the fold barely touches, `L_z = 1`), and the fold's higher-order terms
   `K_l^n K_u` (`n ≥ 2`, present because the lexicographic order maximises the fresh chain
   `0→4→6→7`, nilpotency index `d+1`) leave that survivor at `ρ = 0.9855` — worse than the
   entire Jacobi spectrum (0.9754). At d=2 both splittings race the same `y`-face sawtooth
   and G-S simply accelerates it (§5).

**In one sentence.** *Boundary Gauss-Seidel loses at d=3 not because there are six faces
but because the diamond closure gives every cell a `(d−1)`-dimensional undamped
transmission subspace — which removes the theorem that would otherwise guarantee G-S ≥
Jacobi — and the shipped octant order then folds the boundary coupling in the single most
one-axis-dominated pattern available, accelerating the axes whose modes were not the slow
ones and amplifying the one that was.*

---

## 11. Verdict on `ndim`, and the recommendation

### 11.1 Is `ndim` the right discriminating variable? — **No.**

Measured counterexamples on both sides (§6): three d=2 all-reflective fixtures where G-S
LOSES (`n_GS/n_J` up to 2.859) and three d=3 all-reflective fixtures where it WINS (down
to 0.576). At fixed `ndim`, the sign is moved by optical thickness (§6.3), mesh
refinement, aspect ratio, quadrature order (§6), scattering ratio (§9) and — most
strongly — the octant sweep order (§8, a 2.5× spread in `n_GS/n_J` at one fixture).

**The honest predicate**, if one insisted on writing it, is a joint function of
(a) the fold pattern `(L_a)` the octant order induces, (b) the per-cell DD transmission
spectrum `spec Σ = {1 − 2Σ_tV/D} ∪ {−1}^{d−1}` per ordinate, and (c) the scattering ratio.
It has no cheap closed form, it is not monotone in any single one of them, and `ndim` is
not even a *conservative* proxy — it mislabels a 1-cell 2-D box and a 3-D anisotropic box
in **opposite** directions. **A production default must not branch on it.** (The existing
`is_cartesian and not is_1d` *admission* gate is a different question and is fine — it
gates whether the schedule is well-defined, not whether it is faster.)

### 11.2 Is the zero-leakage d=3 corner physical, or a verification fixture?

Both, but the exposure is smaller than it looks, for three independent reasons:

1. **All-reflective 3-D IS a real production configuration** — it is the standard
   infinite-lattice `k_inf` geometry. Not a verification artefact.
2. **But that configuration is reached through `solve_sn`, which already defaults to
   `inner_schedule="jacobi"`** (`orpheus/sn/solver.py:2113`, and `SNSolver.__init__`
   line 1003). The `"gauss_seidel"` default is **only** on `solve_sn_fixed_source`
   (line 3092) — whose own Notes read *"This entry point exists for L1 verification via
   MMS, not for engineering problems"*.
3. **A real lattice is near-critical, and the inversion is gone there.** At `c ≥ 0.99` the
   d=3 all-reflective box has `n_GS/n_J = 0.80` — G-S wins (§9). The inversion needs
   `c ≲ 0.9` *and* zero leakage, i.e. a strongly-absorbing fully-reflected box.

So the regime where the current default is the slower arm is: **fixed-source, zero
leakage, `c ≲ 0.9`, d=3** — which is the shape of `tests/sn/solve/test_d3_admission.py`'s
absorber gate and very little else.

### 11.3 Recommendation on `solve_sn_fixed_source(inner_schedule=...)`

**KEEP `"gauss_seidel"`. Do not add a dimension branch.** The case, in mechanism terms:

* **Boundary G-S's leverage and its pathology occupy the SAME region** — the zero-leakage
  corner. `B` is the only coupling folded, so with any vacuum face the two arms are within
  ~3 % (the issue's own 0.97 rows) at every dimension. The default therefore only *matters*
  in the corner.
* **In that corner the measured upside exceeds the measured downside.** Best G-S win:
  `n_GS/n_J = 0.158` (d=2 (10,10) ext(1,1) Σ_t 0.8 — Jacobi needs **6.3×** the sweeps).
  Worst G-S loss: `1.892` (d=3 fold pattern (4,1,1)); worst *shipped-order* loss `1.707`.
  Switching the default to Jacobi trades a bounded ≤1.9× loss for an up-to-6.3× one.
* **`ndim` cannot carry the branch** (§11.1), and a branch on a proxy for a variable we
  understand is exactly what a scientific code must not ship.
* **The claim to fix is the docstring, not the default.** Two present-tense-false
  statements are load-bearing here:
  * `orpheus/sn/loss_representation/sweep_schedule.py` module docstring: *"the schedule
    changes only the SI spectral rate (`ρ_J = c` Jacobi vs `ρ_GS ≈ c²` …)"*. `[M]` the
    `ρ_GS = ρ_J²` law is wrong in BOTH directions on real fixtures — d=2 (3,4):
    `ρ_J² = 0.9286` vs measured `ρ_GS = 0.9064` (better than the law); d=3 (3,4,5):
    `ρ_J² = 0.9514` vs measured `0.9855` (much worse). It holds only for a consistently
    ordered 2-cyclic operator, which this is not.
  * `_select_si_splitting` (`solver.py:800-809`) calls `(L+C−B) = M − B_upper`
    **"the regular splitting"**. In Varga's sense a *regular* splitting requires
    `M⁻¹ ≥ 0` and `N ≥ 0`; the DD closure supplies at least `d−1` transmission channels
    equal to `−1`, so it is not one — and that failure is precisely what licenses
    `ρ_GS > ρ_J`. The word should go, or be qualified, because as written it asserts the
    guarantee this investigation measured to be absent.

### 11.4 The change that actually dominates the default question — own the octant ORDER

The strongest lever found is not the schedule choice but a **parameter nobody owns**:
`SweepSchedule.gauss_seidel` iterates `sn_mesh.quad.octants`, so the Gauss-Seidel octant
order is *whatever the quadrature's sign-partition enumeration happened to produce*. It is
never chosen, never documented as a choice, and never tested as one. On the #341 fixture
it spans `n_GS/n_J ∈ [0.771, 1.892]` — a **2.5× spread**, with the inherited order at
1.707 (24th of 25) and the best at 0.771.

Worth its own `module:sn` / `type:improvement` issue: make the order an explicit input to
`SweepSchedule.gauss_seidel` with a stated rule, then choose it by measurement across
fixtures (the `(3,2,1)`-class "one step off maximal dominance" pattern is the candidate
suggested by §8, but it has been measured on ONE fixture and must be swept before being
adopted). Because the fixed point is splitting-invariant, this is a pure-rate change with
no correctness exposure — but note it is NOT free: at d=2 the inherited `(2,1)` is already
the best of the three achievable patterns, so any order change is a trade and must be
measured at both dimensions.

---

## 12. Refuted candidates — each with the structural reason it failed

Per `.claude/rules/process-discipline.md`: a rejection without its reason has to be
re-derived at full cost by the next session.

| # | candidate | why it failed |
|---|---|---|
| 1 | **"6 faces ⟹ a smaller fraction of the boundary coupling updated per pass than d=2's 4 faces"** (the issue's recorded story) | The foldable fraction *does* shrink with `d` — but from the **deferred-reflect correctness rule**, not face count, and the closed form is `(2^d−1)/(d·2^d)` (§2.1). The shrink is far too small: the interpolation it implies, `ρ ≈ ρ_J^{1+2f}`, predicts G-S beats Jacobi at EVERY dimension (`f > 0` always helps), and mispredicts d=3 by 0.024 in ρ on the wrong side of `ρ_J` (§2.2). |
| 2 | **The folded fraction `f = Σ L_a/(d·2^d)` is the controlling quantity** | `ρ_GS` is **non-monotone** in `f`, and the monotonicity *reverses* with dimension: d=2 `f=1/4 → 0.9348`, `f=3/8 → 0.9064` (more is better); d=3 `f=1/8 → 0.9717`, `f=7/24 → 0.9855` (more is worse) (§3.1–3.2). §8 then shows patterns at the SAME `Σ L_a = 6` spanning 0.968–0.987. |
| 3 | **The DD closure's NEGATIVE transmission coefficients (sign pattern) cause it** | A per-axis sign flip is a **diagonal similarity** on the octant hypercube (`D[s,s] = Π_{a∈flip} s_a` sends `T_a → ε_a T_a` and preserves the fold's edge set), so BOTH `ρ_J` and `ρ_GS` are invariant under it — measured identical to 5 decimals across all 8 sign patterns (§4(i)). Sign *indefiniteness* is NECESSARY (it is what voids Varga), but the sign *pattern* is a gauge, not a variable. |
| 4 | **The inversion is octant-graph combinatorics (hypercube + fold pattern + order)** | The `2^d`-scalar model carries the full hypercube, the exact deferred-reflect fold and the exact ordering dependence, and it **never inverts**: `P(ρ_GS > ρ_J) = 0.00` over d ∈ {2,3,4} × all sign patterns × 400 magnitude draws (§4(iii)). The mechanism needs the intra-octant `d × d` DD block the model discards. |
| 5 | **The d=2/d=3 comparison is confounded by angular WINDOWING** (`_maybe_window` fires at `is_cartesian and ndim == 2` only, so the d=2 iterate is harmonic moments and the d=3 one is full-angular) | Refuted analytically AND empirically. `N = S + B` factors through `(moments ⊕ trace)`, so `G_full = A⁻¹ÑΠ` and `G_win = ΠA⁻¹Ñ` are `XY` vs `YX` and share their nonzero spectrum; for a pure absorber `S = 0` and `N = B` reads the trace only. Empirically the windowed d=2 ρ predicts the *unwindowed* sweep-count ratio to 5 % (0.377 predicted vs 0.398 measured, §1.2). Windowing shifts the residual metric's constant, not the rate. |
| 6 | **The 1-group economy used in §6 changes the answer** | Control row: the 2-group and 1-group d=3 (3,4,5) fixtures give ρ_GS = 0.985522 and ρ_J = 0.975409 **identically** — a pure absorber decouples the groups exactly, so the 2-group spectrum is the union (§6). |
| 7 | **Fold-CHAIN depth (`K_l` nilpotency index) is the universal amplifier** (Prediction D) | HELD on 6 of 8 tested rows including a quantitative point estimate, **REFUTED on the two optically-thick d=2 rows**, where an index-2 (star) fold is *worse* than the shipped index-3 chain (0.785 vs 0.672; 0.883 vs 0.792) (§7.3). ⟹ two distinct amplification mechanisms, (a) fold-chain and (b) thick-cell `Σ ≈ −I`; the antipodal order fixes (a) and aggravates (b). |
| 8 | **Both splittings race the same slow mode, so the flip is one quantity changing sign** | True at d=2 (both on the `y`-face sawtooth, `|μ_x| = 0.869` ordinates), **false at d=3**: boundary-G-S's slow mode is on `zmin/zmax` (0.673/0.672 of the trace mass), Jacobi's on `ymin/ymax` (0.658/0.658) (§5). It is two different comparisons. |
| 9 | **`ndim` is the discriminating variable** | Falsified on both sides — three d=2 fixtures where G-S loses, three d=3 where it wins (§6.1–6.2). |
| 10 | **The sweep counts are contaminated by the stopping test** (Signature 9 / L11, or the #340 `max_inner` truncation) | Not applicable: the instrument eigen-solves `M⁻¹N` directly and never runs `SourceIteration.solve`. It reproduces the independently *fitted* ρ of a 1631-sweep real solve to 4 decimals in both arms (§1.2). (`SourceIteration`'s stop is already the ρ-honest equation residual `r = Σgᵢ(ψ_{n−1} − ψ_n)`, `iteration.py:722`, not the increment.) |

---

## 13. Probes written (all in `scratch/`, none touch production)

| probe | what it settles | runtime |
|---|---|---|
| `probe_341_iteration_spectrum.py` | the instrument: `ρ(M⁻¹N)` from the production builders, + linearity and mode-mass reporting. Positive control §1.2. | ~10 s/case (d=3) |
| `probe_341_octant_order.py` | experiment C: ρ_GS under a permuted `Quadrature.octants`, with `Σ L_a` and the measured `B_lower` fraction. Carries the C1 mirror control. | ~100 s |
| `probe_341_octant_model.py` | the `2^d` analytic model: sign-gauge invariance, Varga on positive gains, and the demonstration that octant combinatorics alone never invert. | < 5 s |
| `probe_341_predicate_sweep.py` | the `ndim` falsifiers: mesh × Σ_t × quadrature × aspect ratio at both dimensions. | 244 s |
| `probe_341_mode_shape.py` | which mode each splitting is racing (per-face mass, ordinate class, spatial sawtooth). | ~60 s |
| `probe_341_fold_pattern.py` | all 25 achievable d=3 fold patterns (and all 3 at d=2) → the `max_a L_a > Σ_{b≠a} L_b` law. | ~8 min (d3) |
| `probe_341_antipodal_prediction.py` | Prediction D, held/refuted per row. | ~3 min |
| `probe_341_scattering_ratio.py` | `c → 1`: the inversion disappears near criticality. | ~4 min |

None is a pytest test yet. Two are worth promoting if any of this is acted on:

* **`probe_341_iteration_spectrum.py`** → a `tests/sn/` rate gate asserting
  `ρ(M⁻¹N)` for both splittings on a pinned fixture. It is the only instrument in the tree
  that measures the SI rate **without** running the stopping test, so it is immune to the
  #340 truncation class, and it would turn every future schedule/order change into a
  measured one. (Mark `slow`; the d=3 rows are ~10 s each.)
* **`probe_341_octant_model.py`** → a fast `foundation`-tier gate on the two theorems of
  §4 (Varga holds on the non-negative model; `ρ_J = Σ|g_a|`). Sub-second, no solver.

⚠ Not promotable as-is: the sweep/enumeration probes are print-only surveys, not gates
(`tests/derivations/_promotion_policy.md` — never promote a print-only script).

---

## 14. Open, if anyone picks this up

* **The thick-cell amplification (§7.3 mechanism (b)) is characterised only empirically.**
  Two d=2 rows show it; the `Σ → −I` argument explains why the comparison theorem has no
  content there but not the *sign*. A `Σ = −I + ε` perturbation analysis of
  `(I−K_l)⁻¹K_u` would close it.
* **The `max_a L_a > Σ_{b≠a} L_b` law is exact on 25/25 patterns of ONE d=3 fixture.**
  Whether it survives a second fixture (different aspect ratio, Σ_t, quadrature) is
  untested, and it is already known NOT to transfer to d=2. Do not encode it anywhere
  until it is swept.
* **`A = L + C − S − B` is exactly singular on an all-reflective box** (§1.3), with a
  trace-only, zero-cell-average null space that the SI residual stop cannot see. That is a
  *correctness-adjacent* finding independent of #341: the returned state on such a problem
  differs from the intended uniform field by a frozen null-space component (the 11.26 %
  trace deviation of `d3_absorber_diagnosis.md` §4), and nothing warns. Benign for cell
  averages and the scalar flux; worth a theory-page note and possibly its own issue.
