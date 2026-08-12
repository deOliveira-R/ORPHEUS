# Q6.8 — how to rank a curvilinear SN angular-differencing scheme

A design pass for GitHub **#235**. Deliverable: an instrument SUITE, a
mechanical acceptance protocol an instrument must pass **before** it is
allowed to rank anything, a discriminating fixture design, and the cost.

⛔ **This is a DESIGN document. No production code and no `tests/` file was
touched.** Every number below is from a throwaway probe under
`$CLAUDE_JOB_DIR/tmp/` (`p1_basis.py`, `p2_spectral.py`, `p4_family.py`,
`p5_psi_metric.py`, `p6_scan.py`, `p7_affine.py`, `harness.py`,
`fixture.py`).

---

## 0. Configuration — every number carries its fixture

`[M]` 2026-08-12, branch `refactor/operator-strategy-layers` @ `bea6a367`,
host `.venv/bin/python` (3.14).

* **Fixture family** — `build_cylindrical_anisotropic_mms_case(n_phi=·)`
  unless a row says otherwise: 1-group, homogeneous, `σ_t = 1.0`,
  `σ_s = 0.5`, `R = 5.0`, `Quadrature.folded_product(n_mu=4, n_phi=·)`,
  `bc_left=reflective`, `bc_right=vacuum`,
  `max_inner=500`, `inner_tol=1e-13`.
* **`nx`** is stated on every row. `nx = 320` means *angular-isolated*
  (memo §4.3: at `nx = 80` the spatial error is `≈1.3e-4` and dominates
  every `n_φ ≥ 32` reading).
* **τ is overridden in-process** by monkeypatching
  `orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`
  (its single production consumer is `pole_angular_closure.py:1360`,
  which reads it by module-global name). Every variant is **re-solved to
  its own fixed point** — never a flux converged under a different τ
  (the §F7 confound).
* **Instrument validated before use** (`vv` anti-#17): the unmutated
  harness reproduces the memo's own numbers to every printed digit —
  `3.1503e-3 / 3.2611e-4` at `nx=80`, `n_φ=8/32`, and
  `3.1281e-3 / 2.8285e-4` at `nx=320`.

---

## 1. The headline, in five sentences

1. ⭐⭐ **The reason the current MMS "stops discriminating" is structural,
   not statistical: its entire angular content is the single azimuthal
   harmonic `m = 1`, and the incumbent closure is EXACT on that harmonic
   by construction** — `[M]` `max|E| = 6e-17`. A fixture on which the
   scheme under test has zero closure residual cannot rank that scheme.
2. ⭐⭐ **And the closure is exact on `span{1, μ}` for the incumbent on
   BOTH arms** (`[M]` `4.4e-16` cylinder, `8.9e-16` sphere, at every
   order) — which is precisely the angular content of the diffusion
   limit. ⛔ **#319's diffusion-limit instrument is therefore a
   distance-to-the-incumbent metric for the ANGULAR closure**: it scores
   the incumbent zero and everything else positive, at every order, for
   every material. It is a valid instrument for the *spatial* closure and
   a provable non-catcher for the *angular* one.
3. ⭐⭐ **Zeroing the closure residual makes the answer WORSE**, measured
   in BOTH functionals: on the shipped MMS at `nx=320, n_φ=64` the
   closure-EXACT scheme reads `5.71e-5` (angular-flux `L2`) / `7.17e-5`
   (scalar-flux `L2`) while `τ ≡ ½` — whose closure residual is *nonzero
   at every ordinate* — reads `1.76e-5` / `9.15e-6`, i.e. **3.2× / 7.8×
   better**. The closure error and the α-redistribution truncation
   partially cancel. ⟹ **no pure closure-residual instrument can rank by
   accuracy**, and that includes §9bis.6's P1 closure defect, R&L Eq. 16,
   and any successor.
4. ⭐ **A whole class of instruments is annihilated by the σ_y fold's own
   parity.** A folded-rule ψ is even in ξ ⟹ its azimuthal content is
   cosine-only `{cos mω}` ⟹ for every **odd** harmonic the `w`-weighted
   *signed* closure defect is **identically zero** for every τ obeying
   the reversal identity — `[M]` `≤ 1.7e-16` for shipped, diamond AND
   reversed at `m = 1, 3`. The current ansatz is pure `m = 1`.
5. ⟹ **The suite must be re-scoped.** Ranking by accuracy is one
   instrument with poor resolution near the optimum; the load-bearing
   instruments are **constraints** (admissibility, non-amplification,
   positivity, order) that DISQUALIFY, plus **one** end-to-end accuracy
   metric under a signal-fraction discipline, plus a **diagnostic** layer
   that explains and never votes.

---

## 2. The autopsy, generalised — six ways a ranking instrument dies

Six died; they died of five distinct causes, and the sixth cause is the
one that killed the *fixture*. Naming them is what makes the acceptance
protocol in §4 mechanical rather than a list of good intentions.

| # | death mode | mechanism | who died of it | the pre-flight check that catches it |
|---|---|---|---|---|
| D1 | **Parameter-blind** | the functional does not take the design parameter as an input at all | P1 `c=Σwη²`, BMC β, Lathrop β (memo §9bis.4) | **A1** garbage-substitution |
| D2 | **Chart-relative** | the value depends on a coordinate choice that is not part of the scheme | ν-closure (memo §9bis.3: `τ≡½` closes exactly in ω and fails 16.5 % in η) | **A2** reparameterisation |
| D3 | **Reference-limited** | the reference's own discretisation error ≥ the inter-candidate spread | trajectory-resolvent cross-check (memo §9bis.8b; the tell is that refining the SUT makes agreement WORSE) | **A4-i** reference self-convergence |
| D4 | **Fixture-annihilated** | the functional IS parameter-loaded, but the fixture's content lies in its kernel | the shipped aniso-cyl MMS (§1.1/§1.4 above) | **A3** basis realizability + span |
| D5 | **Mis-correlated** | sensitive, independent, reference-sound — but not monotone in the thing you care about | endpoint defect `D` (`q65_endpoint_defect_findings.md` F10: Pearson r on log `+0.75 → +0.26 → +0.06`) | **A5** ensemble rank correlation |
| D6 | **Incumbent-scoring** | the functional is minimised by the incumbent *by construction* (it measures distance to the design's own defining property) | §9bis.6's η-weighted P1 closure defect; ⛔ **and #319's diffusion-limit instrument, prospectively** (§1.2) | **A6** zero-set identification |

⭐ **D6 is new and it is the most dangerous, because it produces a
confident, tight, reference-free number that always confirms whatever is
already shipped.** It is what `vv` anti-#24(a) calls the BASIS failure,
seen from the other side: the campaign already made this error once
(memo §9bis.10, "my C6 instrument was wrong in its BASIS, and corrected
it inverts"), and #319 is about to make it again with a different
instrument. The check is one line — *solve `instrument(scheme) = 0` for
the scheme* — and if the answer is "the incumbent", the instrument is a
tautology wearing a metric's clothes.

---

## 3. The structural facts any instrument must respect (all `[M]`)

These are properties of **this seam**, measured here, that constrain
what an instrument can possibly see. They are the input to §4–§7.

### F-A The realizable angular basis, and its parity

`[M]` `p1_basis.py`. On `folded_product(n_mu=4, n_phi)` every level
carries `M = n_φ/2` ordinates at **`ω_m = (2m+1)π/(2M)`** — exactly
cell-centred on an equispaced partition of `(0, π)`, edges at
`ω = 0, π/M, …, π`, `[M]` cell widths bit-equal, weights bit-equal
(`ptp = 0`), and **τ is level-independent** (the `sinθ_p` cancels in the
barycentric ratio). Nodes never sit at `ω = 0, π` (`min ξ = 4.98e-2` at
`n_φ=32`).

The σ_y fold keeps `ξ > 0` and doubles the weight ⟹ every representable
ψ is **even in ξ**, i.e. **even in ω**, i.e. its azimuthal content is
`{cos mω}`. Two consequences the fixture designer must carry:

* **any realizable ψ is STATIONARY in ω at both ends of the march**
  (`ψ'(0) = ψ'(π) = 0`, forced by evenness + periodicity) — precisely
  where τ deviates most from ½ (`τ → ¼, ¾`; R&L's `|τ−½|/w` diverges
  `O(M)`). The first-order closure error `(τ−τ*)Δψ'` is suppressed by
  the fold exactly where `(τ−τ*)` is largest.
* the endpoint cell's leading behaviour is **quadratic in ω**, and
  `η = sinθ cos ω` makes it **linear in η** — which is why the
  incumbent's `τ → ¼` is the *right* endpoint value for a folded rule,
  and why R&L's `τ = ½ + O(w)` criterion (derived for a march in μ where
  ψ is generically non-stationary at the ends) does not transfer.

### F-B The incumbent closure is EXACT on `span{1, μ}` — both arms

`[M]` `p7_affine.py`, `f = a + b·(radial cosine)`, worst over two
`(a,b)` pairs and every level:

| rule | `max\|E\|` affine | `max\|E\|` isotropic |
|---|---|---|
| cyl `folded_product(4, 8/16/32)` | `4.44e-16 / 4.44e-16 / 8.88e-16` | `0` exactly |
| sph `gauss_legendre(8/16/32)` | `8.88e-16 / 8.88e-16 / 8.88e-16` | `0` exactly |

This is not a coincidence to be verified — it is the **definition** of τ
(the barycentric coordinate in the radial cosine, BMC Eq. 43 / R&L
Eq. 13c). Any instrument whose trial field lies in `span{1, μ}` is
scoring "distance from the barycentric convention", not accuracy.

**Isotropic is worse still: `E ≡ 0` for EVERY τ** (a convex combination
reproduces constants). An isotropic-ansatz row is a provable
non-catcher for any angular-closure claim — the mechanism behind
`[[lessons-L43i]]`.

### F-C The spectral closure defect, per harmonic and per τ

`[M]` `p2_spectral.py`, level 0, `E_m[f] = τ_m f(ω_{m+½}) + (1−τ_m)
f(ω_{m−½}) − f(ω_m)`, reported as `Σ_m w_m E_m` (signed) and
`Σ_m w_m |E_m|` (absolute), `n_φ = 16`:

| variant | `m=0` s/a | `m=1` s/a | `m=2` s/a | `m=3` s/a | `m=4` s/a |
|---|---|---|---|---|---|
| shipped (arc, P2-in-η) | `0` / `0` | `7.6e-17` / **`1.7e-16`** | `+8.24e-2` / `8.24e-2` | `1.4e-16` / `2.09e-1` | `+1.52e-1` / `3.89e-1` |
| diamond `τ≡½` | `0` / `0` | `1.1e-16` / `2.69e-2` | **`−7.3e-17`** / `1.09e-1` | `6.9e-18` / `2.36e-1` | **`2.2e-16`** / `4.53e-1` |
| reversed (`1−τ`) | `0` / `0` | `2.1e-17` / `5.38e-2` | **`−8.24e-2`** / `1.63e-1` | `9.7e-17` / `3.29e-1` | **`−1.52e-1`** / `5.60e-1` |
| shuffled | `0` / `0` | `1.43e-2` / `6.94e-2` | `+5.38e-2` / `1.43e-1` | `−1.56e-1` / `3.41e-1` | `5.83e-2` / `5.93e-1` |

Three stabilisers, each of which annihilates a *different* candidate
instrument, all readable off this one table:

1. **`m = 0` ⟹ `E ≡ 0` for every τ** (isotropic blindness).
2. **odd `m` ⟹ the SIGNED `w`-weighted defect is machine-zero for every
   reversal-compatible τ.** Mechanism: for the incumbent `E ≡ 0`
   pointwise (F-B); for any τ with `τ_m + τ_{M−1−m} = 1` the defect is
   `∝ cos(mω_m)`, whose discrete azimuthal moment vanishes — *the very
   identity `Σ w η = 0` that makes the manufactured `φ = A(r)` analytic
   is what annihilates the defect in the scalar flux.* **The
   convenience that makes the reference closed-form is the same fact
   that makes it blind.**
3. **`τ ≡ ½` ⟹ the SIGNED defect is machine-zero at EVERY harmonic**
   (its defect is `∝ f(ω_m)` exactly, and every non-aliased discrete
   azimuthal moment vanishes). So the signed functional cannot separate
   the diamond from an exact scheme, at any order, on any fixture.

⟹ **an instrument built on the closure defect must be reported per
harmonic AND in both norms**, because the signed and absolute norms have
*different* stabilisers, and a single scalar hides which one you hit.

### F-D ⛔⛔ Closure-exact is NOT accuracy-optimal — the two errors cancel

`[M]` `p4_family.py`, shipped ansatz, `nx = 320` (angular-isolated),
scalar-flux volume-weighted `L2`:

| scheme | `n_φ=8` | 16 | 32 | 64 | closure residual (F-C) |
|---|---|---|---|---|---|
| shipped (closure-EXACT here) | `3.1281e-3` | `1.1078e-3` | `2.8285e-4` | `7.1658e-5` | **0** |
| blend `w=.25` | `1.4856e-3` | `8.1624e-4` | `2.1273e-4` | `5.4360e-5` | ↑ |
| blend `w=.50` | `1.2623e-3` | `5.5746e-4` | `1.4398e-4` | `3.7255e-5` | ↑↑ |
| blend `w=.75` | `2.2927e-3` | `3.6551e-4` | `7.8689e-5` | `2.0762e-5` | ↑↑↑ |
| diamond `τ≡½` | `3.4258e-3` | `3.4907e-4` | `3.9321e-5` | **`9.1485e-6`** | max |
| reversed | `7.2922e-3` | `1.2665e-3` | `2.8466e-4` | `6.9203e-5` | 2× diamond |

The error is **monotone decreasing in the closure residual** along the
whole `blend(w) = (1−w)τ_ship + w·½` family. That is the opposite of
what a closure-residual instrument asserts.

⟹ ⛔ **Any instrument of the form "how big is the closure's local
truncation residual" is ANTI-correlated with accuracy on this seam.**
This retires, prospectively, the whole family: §9bis.6's P1 closure
defect, R&L Eq. 15/16's `|τ−½|/w`, the endpoint defect `D`, and the
"spectral closure defect" of F-C **as ranking instruments**. They stay
as *diagnostics* — F-C is what lets you say *why* a scheme behaves as it
does, and it is the only cheap thing here that is honestly τ-loaded.

⚠ **And `reversed` ≈ `shipped` to 0.6 %** at `n_φ=32` (`2.8466e-4` vs
`2.8285e-4`) and **inverts** at `n_φ=64`. `reversed` is exactly
`τ → 1−τ`, i.e. the march-orientation flip that `[[lessons-L47a]]` found
the committed gate set to be blind to. **The accuracy fixture is blind
to it too** — for the F-C(2) reason. This is the single sharpest
statement of the problem: *the one error the static gates cannot see is
also the one the accuracy fixture cannot see.*

### F-E The `nx = 80` confound is live in the brief's own table

`[M]` `p4_family.py`, shipped ansatz, same schemes, scalar-flux `L2`:

| scheme | `nx=80`, `n_φ=32` | `nx=320`, `n_φ=32` | `nx=80`, `n_φ=64` | `nx=320`, `n_φ=64` |
|---|---|---|---|---|
| shipped | `3.2611e-4` | `2.8285e-4` | `1.5466e-4` | `7.1658e-5` |
| diamond | `1.3279e-4` | `3.9321e-5` | `1.2735e-4` | `9.1485e-6` |
| reversed | `2.9387e-4` | `2.8466e-4` | `1.3510e-4` | `6.9203e-5` |
| **spread (max/min over 12 schemes)** | **2.7×** | **9.2×** | **1.8×** | **18×** |

⟹ the brief's *"at `n_φ=32` a garbage permutation reads within 1.11× of
the shipped scheme"* is **two effects superposed**: a genuine
`τ ↔ 1−τ` blindness (F-D, survives at `nx=320`), and a `nx=80` spatial
floor of `≈1.3e-4` that is identical for every scheme and compresses the
whole table. `[M]` at `nx=320, n_φ=64` the instrument's dynamic range is
**18×**, not 1.1×. Half the reported non-discrimination is a
configuration artefact the memo already measured (§4.3) and the brief
inherited.

---

## 4. THE ANTI-BLINDNESS PROTOCOL — what an instrument must pass before
## it may rank anything

The brief's working version — *"correlate it against an independent accuracy
measure over a spread of deliberately-varied schemes, and require rank
agreement"* — is **A5** below. It is correct and it is the LAST and most
expensive gate. Putting it first is what made the campaign expensive:
`[M]` **four of the six dead instruments die at A1/A2/A6, none of which
needs a single solve.**

An instrument is a triple `(functional F, fixture/ensemble X, reference R)`.
Run the ladder in order; stop at the first failure.

### A0 — declare the VERDICT TYPE, in the instrument's own docstring

One of exactly three, and the choice binds:

* **CONSTRAINT** — pass/fail, disqualifies a scheme, ranks nothing. Needs
  A1 and A6 only.
* **RANKER** — orders candidates. Needs the whole ladder.
* **DIAGNOSTIC** — explains *why*. **Forbidden from voting**, and must say
  so where it is defined.

Most of the graveyard died of *silent promotion*: a sound diagnostic was
read as a ranker because nobody had written down which it was. (`D`'s
autopsy is explicit about this — F10: "a discriminator, not a selection
criterion".)

### A1 — GARBAGE SUBSTITUTION (kills D1: parameter-blindness). Cost: seconds

Replace the design parameter by (a) a random draw, (b) a **permutation** of
the true values, (c) a constant. `F` MUST move by ≫ its own numerical
noise on all three. Report all four numbers side by side (memo §9bis.4's
probe D is the template — keep it, make it mandatory).

⭐ **Sharpening: (b) is not optional and it must include the SPECIFIC
permutation the committed gate set is blind to.** Here that is
`τ → 1−τ`, the march-orientation flip (`[[lessons-L47a]]`). An instrument
certified against a random draw alone can be exactly blind to the one flip
you fear — and `[M]` on the shipped fixture, it is.

### A2 — REPARAMETERISATION (kills D2: chart-relativity). Cost: seconds

Recompute `F` with the scheme expressed in a different chart of the march
variable (here `ω ↔ η`; in general any diffeomorphism). **`F`'s VALUE may
change; its RANKING may not.** ν-closure fails this outright: `[M]`
`τ ≡ ½` closes exactly in ω and misses by 16.5 % in η (memo §9bis.3), so
its verdict is a report on the chart the prober chose.

### A3 — BASIS REALIZABILITY AND SPAN (kills D4: fixture annihilation). Cost: seconds

Three sub-checks; **A3.3 is the one that would have caught the current
MMS in five seconds.**

* **A3.1 realizable** — every trial mode must lie in the discretisation's
  representable space and respect every symmetry the rule was quotiented
  by. Mechanical test: `|Σ_n w_n f(Ω_n) − ∫f dΩ|` at the quadrature's own
  accuracy. Here the σ_y fold forces **even in ξ**, so the level basis is
  `{cos mω}` — `η` and `ξ²` are in, a bare `ξ` is not (`vv` #24a; the
  campaign has already paid for this once).
* **A3.2 span** — the trial set must contain **both parities** and reach
  `m ≥ 2`. `[M]` F-C(2): odd modes are in the signed functional's kernel,
  so an all-odd trial set is half-blind, and a pure-`m=1` set (the shipped
  ansatz) is entirely blind to the reversal.
* **A3.3 non-degeneracy against the SCHEME'S OWN defining property** —
  evaluate the trial field's closure residual **under the incumbent**. If
  it is machine-zero, the fixture lies in the incumbent's kernel and can
  never rank it. `[M]` shipped ansatz: `6e-17`. **Fail.**

### A4 — SIGNAL FRACTION (kills D3: reference-limitedness, and confounds). Cost: 2–4× the metric

* **A4.1 reference-limitedness** — refine the **SUT**. If agreement gets
  WORSE, the reference's own error is the floor; STOP (`[[lessons-L49]]`).
  Publish the reference's self-convergence ladder and the **usable dynamic
  range as a number**.
* **A4.2 confound** — refine every OTHER discretisation axis; `F` must be
  flat to ~10 % over the last step. `[M]` shipped ansatz at `nx = 80`
  fails (spatial floor `≈1.3e-4`); at `nx = 320` it passes. ⭐ The
  candidate fixture of §7 graded on the ANGULAR flux is flat to **0.08 %**
  already at `nx = 80` — the check is also a cost optimisation.
* **A4.3 dynamic range** — report `max/min F` over the ensemble. Below
  ≈2× the instrument cannot order it and must not be used to.

### A5 — ENSEMBLE RANK CORRELATION (kills D5: mis-correlation). Cost: N_ens solves

The brief's criterion, sharpened on five axes. Four schemes is not a
sample; the sample's SHAPE matters more than its size.

1. ⭐⭐ **A CONTINUOUS homotopy replaces "rank agreement".** Build a
   one-parameter family through the design space —
   `blend(w) = (1−w)·A + w·B`, `w ∈ {0, ¼, ½, ¾, 1}` — and require `F` to
   be **MONOTONE in `w`**. Monotonicity is strictly stronger than
   agreement of a discrete ranking, is falsified by a single non-monotone
   triple, and costs 5 solves. `[M]` it is decisive here: the scalar-flux
   metric is monotone in `w` on every fixture tested; the **angular-flux
   metric is NOT** on the even-harmonic fixtures — which says at once that
   the two are not measuring the same thing, before any correlation is
   computed.
2. **Stratify, and report ρ per stratum.** NEAR (a 2 % jitter of the
   incumbent), MID (the genuine rival — `diamond`, `reversed`), FAR
   (garbage — `shuffled`, `const`). A ρ over all three is dominated by the
   FAR/NEAR split and says nothing about the decision actually being made.
   **This is exactly how `D` died**: `r = +0.75` at `n_φ=8`, where both
   metrics were reading the same gross error, decaying to `+0.06` once
   that resolved. **Require the threshold on NEAR ∪ MID alone.**
3. **The ensemble MUST contain the pair inside the stabiliser you fear**
   (`{τ, 1−τ}`), and `F(τ)/F(1−τ)` is reported as its own line. `1.00` is
   a coverage statement, not a footnote.
4. **≥8 members, each independently RE-SOLVED** to its own fixed point
   (the F7 confound: a flux converged under one scheme is self-consistent
   with it for free).
5. **Pre-register** `F`, the ensemble, the strata and the threshold.
   Otherwise "the ranks agree" is fitted.

### A6 — ZERO-SET IDENTIFICATION (kills D6: incumbent-scoring). Cost: one line

**Solve `F(scheme) = 0` for the scheme.** If the solution IS a shipped
convention, `F` measures distance-to-that-convention and cannot
adjudicate. `[M]` this kills three candidates in one line each:

| instrument | `F = 0` ⟺ | verdict |
|---|---|---|
| any affine-in-cosine trial (⟹ **the diffusion limit**) | τ is the barycentric coordinate = **the incumbent** | cannot rank the angular closure |
| §9bis.6's η-weighted P1 closure defect | same zero set | cannot rank (its own row already reads `6.94e-16`) |
| R&L Eq. 16 `\|τ−½\|/w` | `τ ≡ ½` = **the diamond** | it is a restatement of one CANDIDATE, not a measurement over all |

### A7 — A SECOND, STRUCTURALLY DIFFERENT ANCHOR

A ranking established on one fixture with one functional is a property of
that pair. Reproduce it on a different geometry arm / optical thickness /
harmonic content, or scope the verdict explicitly. `[M]` this is not a
formality here: **two functionals of the SAME solve disagree in sign** —
on the candidate fixture at `n_φ=64` the scalar-flux L2 ranks a *garbage
permutation* `1.60×` BETTER than production while the angular-flux L2
ranks it `17.0×` worse (§5, R1).

---

## 5. THE SUITE — four tiers, and only ONE of them ranks

### Tier 0 — CONSTRAINTS (solve-free, closed-form, DISQUALIFYING)

These carry the load. A scheme that fails one is out, whatever its error
number says.

| id | instrument | sensitive to | provably blind to | cost | reference class |
|---|---|---|---|---|---|
| **C1** | **admissibility P3** — `τ_m ∈ [0,1]` per ordinate | the partition/ordinate relation | everything else; and it is invariant under `τ→1−τ` | µs | closed form (already a production guard, `_assert_tau_within_unit_interval`) |
| **C2** | **non-amplification** — `A(M) = max_m ∏_{k≤m}(1−τ_k)/τ_k`, and `∏_all = 1` exactly | the per-step partition's error-propagation | accuracy; also `τ→1−τ`-invariant in the product leg | µs | closed form (`[M]` `2.41 … 9.44` at `M=2…32`) |
| **C3** | **signed orientation law** `(τ_m − ½)·μ_m ≥ 0` | the march-orientation flip — **the only committed law outside the `τ→1−τ` stabiliser** | magnitude errors | µs | closed form (landed `c33178ef`) |
| **C4** | **α-dome closure** `α_{1/2} = α_{M+1/2} = 0` + `Σ w η = 0` | the quadrature and the α recursion | ⛔ **τ, exactly** — see below | µs | closed form (landed `bea6a367`) |
| **C5** | **positivity of the marched `ψ̂`** | the `(1−τ)/τ` amplification | accuracy | one solve | characterisation only (no authority either way — memo §9bis.12) |

⛔ **C4 is where the brief's "conservation / balance identities" candidate
lands, and it is τ-blind by construction, not by fixture.** `[M]` the
level-summed angular term telescopes,
`Σ_m [α_{m+½}ψ̂_{m+½} − α_{m−½}ψ̂_{m−½}] = α_{M+½}ψ̂_{M+½} − α_{½}ψ̂_{½} = 0`
with `[M]` `α_{1/2} = 0.000e+00`, `α_{M+1/2} = 2.78e-17` — **τ does not
appear in the telescope, for any ψ̂ whatsoever.** Global particle balance
therefore holds identically for a random τ. It is a genuine constraint on
the α recursion and the quadrature; it is worth nothing about the closure.
(`vv` anti-#8, with the mechanism made explicit for this seam.)

### Tier 1 — DIAGNOSTICS (cheap, τ-loaded, FORBIDDEN FROM VOTING)

| id | instrument | what it is genuinely for | why it may not rank |
|---|---|---|---|
| **S1** | ⭐ **spectral closure-defect matrix** `E[m]`, both norms, `m = 0…4` (§F-C) | the only cheap, honestly τ-loaded, pointwise object in the suite. It says WHICH harmonics a scheme mistreats, and it separates the diamond's uniform small defect from the arc's endpoint-concentrated one | **A6**: its zero set is the incumbent (F-B); **F-D**: closure residual is ANTI-correlated with accuracy on this seam |
| **S2** | **endpoint defect `D`** (`q65_endpoint_defect_findings.md`) | consistency of a boundary condition production computes and discards | measured non-predictive (F10, `r → +0.06`) |
| **S3** | **R&L `\|τ−½\|/w`** (Eqs. 15/16) | names a real geometric fact — the endpoint cells are `O(Δω²)`-narrow in η while carrying `O(Δω)` weight, so `τ → ¼` and `\|τ−½\|/w` diverges `O(M)` | **A6**: `F = 0 ⟺ τ ≡ ½`; and its derivation assumes a march in μ where ψ is generically non-stationary at the ends, which the σ_y fold **forbids** (F-A) |

⚠ **S3 deserves its own warning because the brief lists it as pointwise
and therefore fold-proof.** Pointwise it is — and that is not the problem.
The problem is that its consistency argument does not transfer: on a
folded rule every realizable ψ has `ψ'(0) = ψ'(π) = 0`, so the endpoint
cell's leading behaviour is quadratic in ω, hence **linear in η**, hence
`τ = ¼` is the *exact* closure there and `τ = ½` is the wrong one. R&L's
criterion and the fold's parity disagree about the endpoint cells, and
R&L's derivation is the one that assumed what is false here.

### Tier 2 — THE RANKER (exactly one; expensive)

**R1 — angular-flux `L2` against a harmonic-rich analytic MMS, at an
angular-isolated mesh, over a stratified ≥8-member ensemble.**

* **Claim layer**: flux-shape (`vv` §hierarchical taxonomy). Not
  eigenvalue — MMS cannot reach that layer.
* **Reference class**: analytic (closed-form manufactured source, SymPy
  provenance) ⟹ structurally independent in the L11 sense.
* **The graded quantity is `‖ψ_h − ψ_exact‖_{w,V}`, NOT `‖φ_h − φ_exact‖_V`.**
  ⛔⛔ **This is the single most actionable finding of the pass.** `[M]`
  fixture `harmonics=(1,2,3)`, `amps=(1,1,1)`, `σ_t=5.0`, `σ_s=2.5`,
  `R=5`, `folded_product(4, 64)`, `nx = 320`, every row `converged=True`
  (`p8_final.py` §B — this is a NEIGHBOUR of the §7.2 recommendation, not
  the recommendation itself; §7.3 carries the recommended fixture's own
  table):

  | scheme | scalar-flux `L2` | vs shipped | angular-flux `L2` | vs shipped |
  |---|---|---|---|---|
  | shipped | `4.047e-5` | — | `1.561e-4` | — |
  | `jitter 2 %` (NEAR) | `3.909e-5` | **0.97× (blind/inverted)** | `3.249e-4` | **2.08× worse** ✓ |
  | `reversed` | `1.417e-5` | **0.35× ("better")** | `1.781e-4` | 1.14× worse ✓ |
  | `shuffled/A` (GARBAGE) | `2.533e-5` | **0.63× ("better")** | `2.651e-3` | **17.0× worse** ✓ |
  | `shuffled/B` (GARBAGE) | `2.055e-5` | **0.51× ("better")** | `1.251e-3` | 8.0× worse ✓ |

  The scalar flux is `Σ_n w_n ψ_n`, so the closure error enters it
  **signed** — and `[M]` F-C shows the signed defect flips sign under
  reversal and vanishes for `τ≡½`. A scheme therefore **wins by
  cancelling against the τ-independent floor**. On this fixture the
  scalar-flux L2 ranks *two garbage permutations above production* and is
  blind to a 2 % τ jitter. Every existing gate in
  `tests/sn/verification/mms/` grades the scalar flux.
* **Reported outputs** (all four, every time): dynamic range over the
  ensemble; `F(τ)/F(1−τ)`; monotonicity along the blend homotopy; and the
  **error constant** `K = F·M²` rather than the observed order — `[M]`
  every plausible scheme is asymptotically 2nd order (memo §4.5 and
  reproduced here), so **the order is not a discriminator and the constant
  is**.

### Tier 3 — CORROBORATION (admissible only with its own convergence ladder)

* **X1 independent-method cross-check** (trajectory resolvent). ⛔
  **Currently INADMISSIBLE**: its own error is `≈3e-2` against an
  inter-candidate spread of `~1e-2` (memo §9bis.8b). What would make it
  admissible is stated as a measurement, not a hope: refine
  `(n_r, n_mu_axial, n_phi_az, n_traj_quad)` until **refining the SN side
  IMPROVES agreement** — the L49 sign flip is the acceptance test — and
  publish the reference's self-convergence ladder beside the gate's
  tolerance. Until then it may corroborate at `n_φ = 8` only, where the SN
  error genuinely exceeds the floor.
* **X2 the SPHERE arm as the second anchor (A7).** Cheapest available
  independence: a different march variable (μ, not ω), **no fold**, hence
  **no parity restriction on the trial basis** — the odd/even kernel of
  F-C simply does not exist there, so a sphere row can see error classes
  the cylinder structurally cannot. `build_spherical_anisotropic_mms_case`
  already exists; it needs the same A3.3 audit (its ansatz is
  `A + Bμ` — affine in the march variable ⟹ **`[M]` closure residual
  `8.9e-16`, same defect as the cylinder's**).

---

## 6. REJECTED CANDIDATES — one structural reason each

First-class output (`process-discipline`: a rejection without its reason
must be re-derived at full cost). Every row is a *structural* reason, not
"it did not work".

### 6.1 The brief's own candidate list

| candidate | verdict | the one-line structural reason |
|---|---|---|
| **diffusion-limit behaviour (centre flux-dip decay vs optical thickness) — #319** | ⛔ **REJECT as an angular-closure instrument** (KEEP for the spatial closure) | the diffusion limit's angular content is `span{1, μ}`, and `[M]` the incumbent closure is EXACT on `span{1, μ}` on **both** arms at every order (`4.4e-16` / `8.9e-16`) — τ is *defined* as the coefficient that makes affine interpolation exact. `F = 0 ⟺` the incumbent (**A6**, death mode D6). It will confirm whatever is shipped. |
| **angular truncation order in `n_φ`/N at spatially-converged mesh** | ⚠ **KEEP, but the ORDER is not the datum — the CONSTANT is** | `[M]` every reversal-compatible scheme is asymptotically 2nd order (memo §4.5 and reproduced here: `2.58/2.03/1.84` chord, `1.50/1.97/1.98` arc, `3.29/3.15/2.10` diamond). An order gate cannot separate consistent schemes. Report `K = ‖e‖·M²`. |
| **Reed & Lathrop Eqs. 15/16 — `(τ−½)/w` bounded** | ⛔ **REJECT as a ranker; KEEP as a named characterisation (S3)** | two independent reasons. (i) **A6**: `F = 0 ⟺ τ ≡ ½`, so it *is* one candidate wearing a criterion's clothes. (ii) its derivation assumes the march variable's ψ is generically non-stationary at the level ends; the σ_y fold **forbids** that (`ψ'(0)=ψ'(π)=0`, F-A), so at the endpoint cells the exact closure is `τ = ¼`, which is what the incumbent gives and what R&L's criterion calls maximally wrong. Pointwise-ness is real and is not the issue. |
| **more MMS designed to stress the azimuthal variation** | ✅ **ACCEPT — this is §7**, with three constraints the naive version violates | see §7. Naively "add higher harmonics" **degrades** signal-to-noise: `[M]` the τ-independent α-redistribution truncation grows `∝ m²` while the τ-dependent closure error grows `∝ m`, so `m1+m2+m3` resolves the reversal *worse* (`1.05×`) than `m1+m2` (`1.33×`). This is `[[lessons-L40b]]` recurring: the strengthening axis must be the one the claim lives on, and "higher frequency" is the wrong axis here. |
| **an independent method cross-check that is not reference-limited** | ⚠ **CONDITIONAL (X1)** | admissible only with a published reference self-convergence ladder whose own error is < 10 % of the inter-candidate spread. The acceptance test is the **L49 sign flip**: refining the SN side must make agreement BETTER. Today `≈3e-2` reference error against a `~1e-2` spread ⟹ inadmissible above `n_φ = 8`. |
| **conservation / balance identities** | ⛔ **REJECT — τ-blind by construction, not by fixture** | `[M]` the level-summed angular term telescopes to `α_{M+½}ψ̂_{M+½} − α_{½}ψ̂_{½}` with `α_{1/2} = 0`, `α_{M+1/2} = 2.78e-17`. **τ never appears.** Balance holds identically for a random τ. It is a constraint on α and the quadrature (Tier 0, C4). |

### 6.2 Candidates I considered and rejected during this pass

| candidate | the one-line structural reason |
|---|---|
| **scalar-flux `L2` against the MMS** (what every gate in `tests/sn/verification/mms/` currently measures) | ⛔⛔ the closure error enters `φ = Σ w ψ` **signed**, and `[M]` the signed defect flips sign under `τ→1−τ` and vanishes identically for `τ≡½` — so a scheme wins by CANCELLING against the τ-independent floor. Measured: it ranks two garbage permutations **above** production and is blind to a 2 % τ jitter (§5 R1). |
| **the signed `w`-weighted spectral closure defect** as a scalar score | its kernel contains every odd harmonic AND all of `τ≡½` (F-C) — two stabilisers, both large. |
| **a "gold τ" oracle** (the τ that closes exactly for the manufactured ψ, giving a true zero point) | only exists if ψ is angle-separable, `ψ = F(r)g(Ω)`; and for `g` monotone in the radial cosine the gold τ **IS the incumbent**, so the construction reduces to D6. For non-monotone `g` (even harmonics) the gold τ leaves `[0,1]` and the production P3 guard refuses it. |
| **`k_eff` / any derived scalar of an eigenvalue solve** | `vv` Mode 12 + the memo's own §4.7: `[M]` the 1-group homogeneous cylinder artifacts move by `1.1e-10`/`1.1e-11` under a deliberate `τ := 0.7`. Degenerate. |
| **`folded_product(·, 4)` (M = 2) fixtures** | `[M]` memo §4.7: the ω-midpoint and chord partitions are **bit-identical** at `M = 2`. No `n_φ = 4` fixture can ever see a partition change. |
| **the endpoint defect `D` as a τ vote** | already decided negatively and in writing (`q65_endpoint_defect_findings.md` F10). Any future `D`-based τ argument must cite that section first. |
| **isotropic MMS rows** | `[M]` `E ≡ 0` for every τ, exactly (a convex combination reproduces constants). Provable non-catcher — the mechanism behind `[[lessons-L43i]]`. |
| **an ansatz with a bare `ξ` (odd in ξ) term** | unrealizable on a folded rule — `quad.mu_y` samples `\|ξ\|`, `[M]` `Σwξ = +6.703 ≠ 0` folded vs `0.000` unfolded (memo §9bis, `vv` #24a). The campaign made this exact error once. |
| **"tighten the existing `test_curvilinear_aniso_convergence` tolerance"** | it grades the scalar flux on a pure-`m=1` ansatz — both defects above at once. Tightening a blind gate buys nothing. |

---

## 7. THE DISCRIMINATING FIXTURE

### 7.1 The three constraints, in the order they bind

1. **REALIZABLE** — even in ξ ⟹ azimuthal content `{cos mω}` only. Write
   the harmonics as the harmonic polynomials `h_m = Re[(η + iξ)^m]`
   (`h_1 = η`, `h_2 = η²−ξ²`, `h_3 = η³−3ηξ²`); then `∂h_m/∂ω = −m·k_m`
   with `k_m = Im[(η+iξ)^m]`, and the manufactured source is closed-form.
   ⭐ Free bonus: `[M]` `|Σ_n w_n h_m| ≤ 4.4e-16` for `m = 1…4` at
   `n_φ = 8/16/32/64` on `folded_product(n_mu=4, ·)` (the equispaced
   azimuthal circle's own moment identity, preserved bit-exactly by the
   fold; `Σ w = 4π` to the last bit), so **`φ_exact = A(r)` stays analytic
   and quadrature-independent** for any harmonic mix — the property the
   existing builder relies on, retained. ⚠ It fails once `m ≡ 0 (mod
   n_φ)`, so cap `m < n_φ_min` of the ladder.
2. **OUT OF THE INCUMBENT'S KERNEL** (A3.3) — the ansatz MUST contain an
   **even** harmonic. `m = 1` alone is affine in η ⟹ zero closure
   residual for the incumbent at every `r` ⟹ the fixture is in its kernel.
3. **DO NOT OVERSHOOT** — `[M]` adding `m = 3` on top of `m1+m2` makes
   the reversal resolution *worse* (`1.33×` → `1.05×` in the ψ metric,
   `n_φ=64, nx=320`), and `m3`-alone is reversal-blind in ψ (`1.01×`).
   *Hypothesis for the mechanism, NOT measured:* the τ-independent
   α-redistribution truncation grows `∝ m²` (it differentiates ψ in ω)
   while the τ-dependent closure error grows only `∝ m`, so raising the
   frequency raises the floor faster than the signal. The consequence is
   measured; the mechanism is a guess and one probe (`m2`-only vs
   `m4`-only floors) would settle it.

### 7.2 The recommendation

```
psi(r, Omega) = [ A(r) + B(r)*h_1 + C(r)*h_2 ] / W        h_m = Re[(eta+i xi)^m]
A = sin(pi r/R)          B = C = (r/R)(1-r/R) cos(pi r/R)
sigma_t = 5.0   sigma_s = 2.5   R = 5.0   (25 mfp; c = 0.5)
folded_product(n_mu=4, n_phi in {16, 32, 64})
bc: reflective at r=0, vacuum at r=R      GRADE: ||psi_h - psi_exact||_{w,V}
```

Why each knob, `[M]` from the 8-fixture scan (`p6_scan.py`) and the
ensemble run (`p9_recommend.py`):

* **`h_2` present** — turns the reversal from invisible to visible:
  ship/rev in the ψ metric goes `1.00×` (`m1` only) → `1.33×` (`m1+m2`).
* **`σ_t = 5` (thick)** — the best of the three optical regimes tested,
  but by a MODEST margin in the metric that counts: `[M]` ψ-metric
  ship/rev at `n_φ=64, nx=320` is `1.46×` thick, `1.33×` at `σ_t=1`,
  `1.36×` at `σ_t=0.2`. ⚠ **In the φ metric optical thickness looks
  decisive** (`2.70×` thick vs `1.08×` thin) **and that reading is an
  artefact of the cancellation φ suffers from** — a caution against
  choosing a fixture knob on a disqualified functional. Take `σ_t = 5`
  for the ~10 % gain and because the thick regime also exercises the
  regime #319 cares about; do not claim more for it.
* **`C/B = 1`** — `[M]` the amplitude sweep saturates: ship/rev
  `1.15 / 1.32 / 1.44 / 1.50 / 1.52×` at `c₂ = 0.25/0.5/1/2/4`. Beyond
  `c₂ ≈ 2` you are only adding τ-independent error.
* **`nx`** — `[M]` the ψ metric is spatially clean at `nx = 80` for
  `n_φ ≤ 32` (`0.9 %` / `4.7 %` above the `nx=320` value) and needs
  `nx = 160` at `n_φ = 64` (`nx=80` is `26 %` high). Grade at
  `nx = 80` for the screening ladder, `nx = 160` for the `n_φ=64` rung.
  The φ metric is *never* clean — `[M]` at `n_φ=64` it wanders
  **non-monotonically** in `nx` (`5.028 / 4.335 / 4.582 e-5` at
  `nx = 80/160/320`), the signature of cancellation.

### 7.3 What it actually resolves — the honest resolution table

`[M]` `p9_recommend.py`, recommended fixture, `nx = 80`, ψ-`L2` relative
to the shipped scheme:

| ensemble member | stratum | `n_φ=16` | `n_φ=32` | `n_φ=64` |
|---|---|---|---|---|
| `blend .25/.50/.75` | NEAR | `1.15 / 1.28 / 1.32×` | `1.05 / 1.09 / 1.10×` | `0.97 / 0.96 / 0.95×` |
| `jitter 2 %` | NEAR | `1.12×` | `1.82×` | **`3.96×`** |
| `diamond τ≡½` | MID | `1.34×` | `1.10×` | `0.95×` |
| `reversed (τ→1−τ)` | MID | `1.61×` | `1.44×` | `1.25×` |
| `const .60` | FAR | `2.17×` | `3.83×` | `6.17×` |
| `jitter 5 %` | FAR | `2.22×` | `4.34×` | `9.33×` |
| `shuffled/A`, `/B` | FAR | `4.29 / 11.98×` | `14.24 / 17.11×` | **`37.76 / 17.55×`** |
| **dynamic range** | | `12×` | `17×` | **`40×`** |

Against the same solves, the **scalar-flux** metric at `n_φ=64` reads
`shuffled/A = 1.18×`, `shuffled/B = 1.21×`, `jitter 2 % = 1.04×`,
`reversed = 0.83×` — dynamic range `2.1×`, garbage indistinguishable from
production, and `reversed` "winning".

⛔⛔ **THE HONEST BOTTOM LINE, and it is a decision-relevant finding, not
a caveat.** With the best fixture and the best functional I could build:

* garbage is separated at **17–40×** ✅
* a **2 % random τ perturbation** at **2–4×** ✅
* the **march-orientation reversal** at **1.25–1.6×** ✅ (from `1.00×` —
  this is the fixture's whole reason to exist)
* ⛔ **`shipped` vs `diamond τ≡½` is NOT resolved**: `1.34×` at `n_φ=16`,
  `1.10×` at 32, `0.95×` at 64 — the sign of the verdict **flips with
  `n_φ`**, and the scalar-flux metric on the same solves says `diamond`
  wins at every order. **The two candidates the campaign has actually
  been arguing about sit inside the instrument's resolution.**

⟹ that question is **not decidable by accuracy on a manufactured
fixture**, and the Q5.6.4 resolution — decide it on Tier-0 constraints
plus the primary source (Morel & Montry's own cylinder appendix) — was
not a fallback, it was the only sound route. An instrument suite's most
valuable output can be *"this comparison is below my resolution"*, stated
with the number.

---

## 8. COST

`[M]` on this host, `.venv/bin/python`, serial. One `solve_sn_fixed_source`
on the recommended fixture: `0.35 s` at `nx=80`, `0.83 s` at `nx=320`.

| tier | what runs | cost | when |
|---|---|---|---|
| **Tier 0** (C1–C4) | closed-form, no solve | **< 10 ms** total, all orders | every candidate, every commit |
| **Tier 0** (C5 positivity) | 1 solve | `0.4 s` | every candidate |
| **Tier 1** (S1–S3 diagnostics) | closed-form, no solve | **< 50 ms** | on demand, when explaining a verdict |
| **A1/A2/A3/A6 pre-flight** | closed-form | **< 1 s** | **before any solve** — this is where 4 of 6 dead instruments die |
| **SCREEN — R1 short** | 6-member ensemble × `n_φ ∈ {16,32}` × `nx=80` | **≈ 5 s** | per candidate scheme |
| **CONFIRM — R1 full** | 12-member stratified ensemble × `n_φ ∈ {16,32,64}` × `nx ∈ {80,160}` | **≈ 60 s** | per candidate that survives screening |
| **A4.2 confound ladder** | 1 scheme × 3 `n_φ` × 4 `nx` | **≈ 25 s** | once per fixture, not per scheme |
| **A7 second anchor** (sphere arm) | same shape on `build_spherical_anisotropic_mms_case` | **≈ 60 s** | once per accepted verdict |
| **X1 resolvent corroboration** | reference self-convergence ladder + 4 SN solves | **minutes–hours** (the reference dominates) | only if X1 is to be made admissible |

⟹ **the whole suite for one candidate scheme is ≈ 2 minutes**, and the
part that kills most bad instruments costs **under one second**. The
measured cost of the campaign to date is dominated by instruments that
would have failed a `< 1 s` pre-flight.

### Tiering for CI

* **Every commit** — Tier 0 (C1–C4) + S1's `m=0/1/2` rows as
  characterisation. `< 100 ms`, no solve, no `slow` marker.
* **`-m slow`, per candidate scheme** — R1 short (`≈ 5 s`).
* **Never in CI** — R1 full and X1. They are *design-decision*
  instruments, run in a scratch harness and reported in the plan; only
  the accepted scheme's resulting error ladder becomes a committed
  regression gate.

---

## 9. CLAIM LAYERS, AND THE CARDINAL RULE

Stated explicitly because the suite is 1-group and that must be
justified, not assumed (`vv` §1-group degeneracy, anti-#3):

* **Every instrument here sits at the CONVERGENCE-ORDER or FLUX-SHAPE
  layer. None of them makes an eigenvalue claim, and none may.** MMS is
  source-driven and cannot reach the eigenvalue layer (`vv` §pillars).
* ⟹ **the 1-group fixture is LEGITIMATE here**, for the same reason a
  convergence-RATE claim may be 1-group (`si_convergence_rate_verification.md`):
  the degeneracy the Cardinal Rule bars is `k = νΣ_f/Σ_a` being
  flux-shape-independent, and no `k` is being claimed. The declaration
  is the price of the exemption — write it in the gate's docstring.
* ⛔ **If #235 ever wants to say "scheme X gives a better `k_eff`", that
  is a different pillar and a different fixture** (≥2G, heterogeneous,
  closed-form or semi-analytical reference). The memo's §4.7 already
  measured why the current cylinder `k_eff` artifacts cannot serve:
  `[M]` the 1-group homogeneous ones move by `1.1e-10` under a
  deliberate `τ := 0.7`.

---

## 10. WHAT TO BUILD, IN ORDER

Nothing below has been built — this is a design pass and the user
reviews before anything lands.

1. **`tests/sn/sweep/curvilinear/test_closure_constraints.py`** — Tier 0
   as a committed gate set (C1–C4 already have homes; C3 landed at
   `c33178ef`, C4 at `bea6a367`). The missing piece is **one gate per
   constraint that states its OWN blindness in its docstring** —
   specifically that C1, C2's product leg and C4 are all invariant under
   `τ→1−τ` and that C3 is the only one that is not (`[[lessons-L47a]]`).
   `< 100 ms`, no solve.
2. **`orpheus/derivations/discrete/sn/angular_differencing.py`:
   `spectral_closure_defect(quad, coord, *, harmonics)`** — S1, the
   diagnostic. Returns the `(harmonic × level)` matrix in both norms.
   Explicitly **DIAGNOSTIC** in its docstring, with the A6 zero-set
   statement ("this functional is minimised by the shipped convention by
   construction; it may not vote").
3. **`build_cylindrical_harmonic_mms_case(harmonics, amps, …)`** in
   `orpheus/derivations/continuous/mms/sn.py` — the §7.2 fixture, with
   SymPy provenance (`derive_cylindrical_harmonic_mms`) in the
   `algebra-of-record` Branch-1 shape, and a foundation test that
   `Σ_n w_n h_m = 0` for every shipped `n_φ` (the fact that keeps
   `φ_exact = A(r)` analytic).
4. **A `psi_l2_ladder` helper beside `scalar_flux_l2_ladder`** in
   `tests/sn/_test_helpers.py` — the angular-flux norm as a single
   source of truth, since §5 R1 makes it the graded quantity.
   ⚠ **Consider re-grading the existing curvilinear MMS gates on it**
   as well; that is a scope decision for the user, and it is the
   `[[lessons-L47a]]` repair applied to the accuracy tier.
5. **A scratch adjudication harness** (not in `tests/`) implementing A1–A7
   for a candidate scheme, emitting the §7.3 resolution table. `≈ 2 min`
   per candidate.
6. **The `#319` brief needs the §1.2 correction before it lands** — its
   diffusion-limit instrument is sound for the spatial closure and a
   provable non-catcher for the angular one. A `numerics-investigator`
   is running it now; the finding is cheap to hand over
   (`p7_affine.py`, six lines of output).

---

## 11. RESIDUAL GAPS — measured absent, not assumed

* ⛔ **The SPHERE arm was not swept.** Only F-B (`8.9e-16` closure
  exactness on `span{1,μ}`) was measured there. The sphere has no fold,
  hence **no parity restriction on the trial basis**, so its blindness
  structure is genuinely different and is the cheapest A7 anchor — and
  `[M]` its shipped ansatz `A + Bμ` has the identical A3.3 defect, so it
  needs the same fix.
* ⛔ **Heterogeneous and multi-group fixtures were not tested.** Every
  number here is 1-group homogeneous. `vv` anti-#4 says a
  homogeneous-only verification is incomplete; the counter-argument here
  is that the claim layer is flux-shape on a manufactured solution, but
  it is an argument, not a measurement.
* ⚠ **A `ConvergenceCertificateError` is a Tier-0 signal I stumbled on
  and did not develop.** `[M]` `τ ≡ 0.40` and `τ ≡ 0.30` make
  `solve_sn_fixed_source` raise
  *"the within-group solve claimed convergence … but the honest equation
  residual is `2.851e-10`"* at some `(n_φ, nx)` and not others. A τ
  scheme that destroys the SI fixed point's certificate is disqualified
  on robustness grounds; that deserves its own C-row and a deliberate
  probe rather than an incidental one.
* ⚠ **The `m ∝ m²` floor mechanism (§7.1.3) is a hypothesis.** One probe
  (`m2`-only vs `m4`-only floors at fixed amplitude) settles it.
* ⚠ **No Spearman ρ was computed.** A5 is specified but not executed —
  §7.3's resolution table is the raw material for it, and running A5
  properly is part of item 5 above, not of this design pass.
* ⚠ **`_carrying_levels` / `_tau_per_level` remain private**, so the
  probes here monkeypatch a module-global. A production adjudication
  harness needs a public surface (already noted in
  `q65_endpoint_defect_findings.md` "Open / not done").
