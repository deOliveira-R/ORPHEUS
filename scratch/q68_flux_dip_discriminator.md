# Q68 — the flux-dip discriminator (issue #319), and what it says about #235

`[M]` 2026-08-12 · branch `refactor/operator-strategy-layers` @ `bea6a367`
· `.venv/bin/python -O` (Py 3.14), serial.
Harness: `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/q68_*.py`.

> **Status: COMPLETE.** Every number carries its configuration; a number
> without its fixture is not reusable.
>
> **Data-integrity audit of the whole campaign:** `[M]` **251 solves, 251
> `converged=True` and `fully_converged=True`, 0 warnings, 0 solves lost**;
> `max |balance_defect| = 0` everywhere; `inner_tol = 1e-12` on 250 of 251
> (one row — `sphere_fixedc / C reversed / Σ_t·R=300 / c=0.999` — needed a
> documented back-off to `1e-10` after the #340 certificate refused `1e-12`),
> plus 100 solve-free `β` rows. Raw record with every φ(r) profile:
> `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/q68_results.jsonl`.
> **Every variant is RE-SOLVED to its own fixed point** — no metric is ever
> read off a flux converged under a different τ.

---

## 0. The question, and the correction the primary source forces on its framing

#319 asks: on a thick, scattering-dominated (diffusive) problem, does the
**centre flux dip** DECAY with optical thickness for **(A)** the production
Morel–Montry weighted diamond (`τ` from
`orpheus.sn.sweep.pole_angular_closure.morel_montry_tau_per_level`) and
**PERSIST** for **(B)** plain angular diamond `τ ≡ ½`? "The discriminator is
the DECAY RATE with thickness; equal decay rates REFUTE the
angular-consistency reading."

⛔ **The primary source does not predict that shape, and the framing must be
corrected before the numbers are read.** Morel & Montry (1984), *Analysis and
Elimination of the Discrete-Ordinates Flux Dip*, TTSP 13(5):615–633 — local
scan `scratch/literature/Morel-Montry(1984)…pdf`, OCR sidecar
`scratch/literature_ocr/Morel-Montry(1984)…md` — says, in its own words:

* p. 8 (printed 621), the **10 mfp** pure-scattering sphere with **angular
  diamond**: *"the diffusion solution fails to predict the **very slight**
  flux dip which occurs in the S₂ solution."*
* p. 13 (printed 626), the **1 mfp** sphere: *"The flux dip is **much more
  severe** for this optically thin problem than it is for the optically thick
  test problem."*
* p. 13, on the weighted diamond: *"We have tested the angular
  weighted-diamond scheme in a wide variety of calculations, and have **never
  observed a dip**."*

⟹ For **(B)** the dip is **largest when THIN** and decays as the problem
thickens; for **(A)** it is *absent at every thickness*. So the honest
discriminator is **not** "does the decay rate differ" — both may decay — it is
**the RATIO (B)/(A) at each thickness**, i.e. whether (A) sits at the
instrument's floor while (B) sits decisively above it. The decay rates are
reported anyway, because #319 asks for them and because a *fitted rate* is a
falsifiable summary; but the verdict is stated on the ratio, and the reason is
recorded here so the next reader does not re-derive it.

**The mechanism M&M give (this is what makes the metric interpretable).**
Taking the 0th and 1st angular moments of the angularly-discretised spherical
S\_N equation under a linear-in-μ ansatz `ψ_m = φ + 3Jμ_m` yields (their
Eqs. 5–7a) a *diffusion* equation with a **corrupted coefficient**

```
    −(1/r²) d/dr [ r² D dφ/dr ] + σ_a φ = Q ,     D = 1 / (3 (σ_t + 2β/r))
    β = 3 Σ_m μ_m ( α_{m+1/2} μ̃_{m+1/2} − α_{m−1/2} μ̃_{m−1/2} )      (Eq. 6a)
```

where `μ̃` are the **cell-edge cosines IMPLIED BY THE ANGULAR CLOSURE**
(`μ̃_{m+1/2} = (μ_m − (1−τ_m) μ̃_{m−1/2}) / τ_m`, seeded `μ̃_{1/2} = −1`).
`β ≠ 0` is a `1/r` corruption of `D` that blows up at the origin — that is the
dip. M&M's Eqs. (16a)/(16b) choose `τ_m = (μ_m − μ_{m−1/2})/(μ_{m+1/2} −
μ_{m−1/2})` with `μ_{m±1/2}` the **standard weight-partition** edges
(`μ_{m+1/2} = μ_{m−1/2} + 2W_m`), which makes `μ̃ ≡ μ` and hence **`β ≡ 0`
identically** (their Eqs. 17–19). That is precisely the shipped `τ`.

Two riders from the source that bear on this experiment:

* the elimination is claimed **"in conjunction with the Miller–Alcouffe
  procedure"** — the slab-geometry starting-direction equation plus setting
  the origin fluxes on a level to the starting-direction value. ORPHEUS's
  route-(a) marched ψ½ seed is the analogue; it is a **confound to name**, not
  one this experiment can separate.
* **coarse spatial mesh breaks the correspondence** (p. 16): the starting flux
  is then *over*-estimated relative to the weighted fluxes, `β(r)` is
  effectively non-zero **but POSITIVE**, and "the flux dip is nonetheless
  eliminated because β(r) is always positive". So a coarse mesh does **not**
  manufacture a dip — it manufactures a *bump*. Useful: it means the
  spatial-artefact trap has a **sign**.

---

## 1. Fixture, instrument, and the positive control

**Materials.** Homogeneous ball/cylinder, radius `R = 1`, uniform isotropic
volumetric source `Q = 1` per group (installed as `Q_n = Q/W` per ordinate —
the producer-side `1/W` projection), vacuum outer surface. One-group is
M&M's own fixture and is *not* the degenerate case here (this is a
fixed-source **flux-shape** claim, not an eigenvalue claim); a two-group
asymmetric companion is run separately. Every mixture is gated with
`Mixture.assert_balanced()` (`σ_t == σ_c + σ_f + Σ_g' σ_s[g,g']`) — an
inconsistent mixture makes two legitimate references disagree with no bug in
either.

**The dip metric `D_fit`.** Least-squares fit of an EVEN polynomial in `r/R`
(`Σ_k a_k (r/R)^{2k}`, k = 0..deg) over the cells with `0.15 ≤ r/R ≤ 0.60` —
a window that excludes both the origin anomaly and the outer transport
boundary layer — then

```
    D_fit = ( φ_fit(r_0) − φ(r_0) ) / φ_fit(r_0)          [POSITIVE = depressed centre]
```

at the innermost cell centre. Even-in-`r` is the right basis: M&M's own exact
diffusion solution for this fixture (their Eq. 8) is *exactly* quadratic in
`r`. The fit's own max relative residual over the window is reported with
every number; a `D_fit` below that residual is **not resolvable** and is
reported as "at the floor", never as zero.

**`[M]` POSITIVE CONTROL — the τ monkeypatch bites on the sphere.**
`Quadrature.gauss_legendre(2)` (μ = ∓0.57735, w = 1, Σw = 2):

| variant | installed τ | max\|Δφ\|/max\|φ\| vs shipped |
|---|---|---|
| A shipped (M-M) | `[0.42265, 0.57735]` | 0 (reference) |
| B diamond τ≡½ | `[0.5, 0.5]` | **3.031e-02** |
| C reversed | `[0.57735, 0.42265]` | 5.753e-02 |
| D shuffled | (S2: only 2 permutations ⟹ **identical to C**) | 5.753e-02 |

⚠ At S2 the "shuffled" control is *not independent* of "reversed" — there are
only two permutations of a 2-element multiset. It is reported at S4/S8 where
it is genuinely distinct; at S2 the two rows are one row.
All four variants preserve `Π (1−τ)/τ = 1` exactly, so the tree's
neutral-stability gate is green for every one of them — whatever separates
them here is invisible to that gate.

**`[M]` FIRST READING — Morel & Montry's own Fig.-3 fixture reproduced.**
Sphere, `c = 1.0` (pure scattering), `σ_t·R = 10`, GL-S2, `nx = 100`
(10 cells/mfp), `inner_tol = 1e-12`, `max_inner = 200000`, **all solves
`converged=True`, no warnings**:

| variant | φ(r_0) | fit(r_0) | **D_fit** | fit resid |
|---|---|---|---|---|
| **A shipped (M-M)** | 5.673792 | 5.673814 | **+3.95e-06** | 1.0e-05 |
| **B diamond τ≡½** | 5.841334 | 5.846951 | **+9.61e-04** | 1.4e-04 |
| C reversed | 5.991726 | 5.999018 | +1.22e-03 | 3.8e-04 |

**(B)/(A) = 243×**, and (A) sits *below* its own fit residual. The dip profile
(first 10 cells, `(fit−φ)/fit`) shows the anomaly is localised in the first
~1 mfp and changes sign beyond it — the classical shape:

```
A shipped : +3.95e-06 +2.64e-06 +1.42e-06 +3.44e-07 −5.82e-07 −1.37e-06 …
B diamond : +9.61e-04 +6.74e-04 +4.42e-04 +2.67e-04 +1.34e-04 +3.46e-05 …
C reversed: +1.22e-03 +5.89e-04 +1.30e-04 −1.84e-04 −3.97e-04 −5.38e-04 …
```

This is M&M's reported result reproduced qualitatively *and* in sign: the
weighted diamond removes the dip; the plain diamond leaves a "very slight"
one at 10 mfp.

---

## 2. Three instruments, and what each one can and cannot see

`vv-principles` #24: an instrument that ADJUDICATES between design candidates
owes a basis check and a rank-correlation check *before* its verdict is used.
Three are carried here; the report states a claim only where they agree.

| | what it is | reference-free? | fit-free? | known limitation |
|---|---|---|---|---|
| **I1 `D_even`** | centre residual vs an EVEN polynomial fit of the interior | yes | **no** | conflates the localised dip with M&M's *global* `2β(r−R)` term; the `r→0` extrapolation is the weak link |
| **I1' `D_mm`** | same, with M&M's own model basis `{1, r, r², r⁴}` (the linear `β` term included) | yes | no | isolates the LOCALISED dip, but the extra odd power makes the extrapolation noisier |
| **I2 `c_s`** | Morel–Montry's **effective starting cosine**, Eq. (10) | yes | **yes** | only interpretable where `ψ` really is affine in the level's angle — degrades at high `N` outside the asymptotic limit (see §5) |
| **I3 `β`** | M&M Eq. (6a) with the **τ-implied** edge cosines | yes | yes, **solve-free** | derived for the SPHERE only; the transfer to a cylinder level is **refuted** below |

**⚠ `D_even` vs `D_mm` — why both.** M&M Eq. (9) says the S\_N-vs-diffusion
error is `φ_Sn − φ_e = 2βQ(r − R)` — a term **linear in r**, i.e. a *global*
tilt, not a localised dip. An EVEN fit cannot represent it, so at large
`σ_t·R` the "dip" I1 reports is really that tilt showing up as a near-uniform
offset near the origin. `[M]` sphere, `c = 1`, `σ_t·R = 50`, GL-S2, 4
cells/mfp — the first 8 cells of `(fit−φ)/fit` are essentially FLAT, which is
a tilt, not a dip:

```
A shipped : +1.02e-05 +9.39e-06 +8.60e-06 +7.90e-06 +7.27e-06 +6.69e-06 …
B diamond : −4.19e-04 −4.38e-04 −4.46e-04 −4.47e-04 −4.42e-04 −4.33e-04 …
```

whereas at `σ_t·R = 10` the same rows are monotone and localised within
~2 mfp — a genuine dip:

```
A shipped : +7.17e-05 +5.41e-05 +3.88e-05 +2.66e-05 +1.70e-05 +9.67e-06 …
B diamond : +1.00e-03 +5.30e-04 +2.48e-04 +9.76e-05 +2.36e-05 −6.32e-06 …
```

Both components are angular-consistency defects and both are killed by
`β = 0`, so the discriminator is unaffected — but the WORD "dip" is only
accurate at moderate thickness, and the report says "centre anomaly" where
it is a tilt.

**`[M]` I2 needed a geometry fix, and the wrong version is a false artefact.**
On the CYLINDER the M&M ansatz must be taken **level-locally**: the angular
flux on the axis is azimuth-independent but genuinely **polar-angle
dependent** (M&M appendix, p. 18). Using the all-level moments
`(ψ_s − φ)/(3J)` reads `c_s = +2.76` at cell 0 for the SHIPPED τ on
`fold(4,8)` — a pure artefact of comparing one level's starting flux with
the all-level average. The correct, geometry-uniform form uses the LEVEL's
own zeroth and first azimuthal moments,

```
    a_l = Σ_m w_m ψ_m / Σ_m w_m ,   k_l = Σ_m w_m c_m ψ_m / Σ_m w_m c_m² ,
    c_m = cos ω_m = μ_m / sin θ_l ,   c_s = (ψ_s,l − a_l) / k_l  →  −1
```

which **reduces to M&M Eq. (10) exactly on the sphere** (verified
bit-for-bit: `c_s(r_0) = −0.98836` under both forms) and gives sane cylinder
values (`−0.953` shipped, `−1.117` diamond, `−1.316` reversed at level 0).

**`[M]` I2 reproduces a published figure.** Sphere, `c = 1`, `σ_t·R = 10`,
GL-S2, **nx = 100** (M&M's own 100 zones), `c_s(r)` in mfp from the origin:

| variant | 0.05 | 0.45 | 0.95 | 2.95 | 4.95 | 9.45 |
|---|---|---|---|---|---|---|
| **A shipped (M-M)** | −0.9884 | −0.9989 | −0.9996 | −0.9999 | −0.9992 | −0.7787 |
| **B diamond τ≡½** | **−1.3157** | −1.1714 | −1.1274 | −1.0667 | −1.0446 | −0.7825 |
| C reversed | −1.9143 | −1.3925 | −1.2701 | −1.1285 | −1.0838 | −0.7853 |

M&M's Fig. 4 reports the Gauss-S₂-**diamond** effective starting cosine
"reaching a minimum value at the origin of approximately **−1.35**"; measured
here **−1.3157** (2.5 %). Their Fig. 6 reports the weighted-diamond curve
"transitions from a correctly high value near the surface to an essentially
constant value of −1"; measured `−0.7787 → −0.9999`. This is the paper's own
diagnostic, reproduced, and it is a **structurally independent** confirmation
of I1 (a different functional of a different field — the half-angle seed and
the angular moments, not the scalar-flux profile shape).

**`[M]` I3 `β` is τ-loaded on the sphere and vanishes exactly for the shipped
τ** — the theory's own identity (M&M Eqs. 17–19), and it is **solve-free**:

| quadrature | A shipped | B diamond τ≡½ | C reversed | D shuffled |
|---|---|---|---|---|
| GL-S2 | **0.000000e+00** | +1.547005e-01 | +2.679492e-01 | (=C) |
| GL-S4 | −8.33e-17 | −2.679368e-03 | −1.162024e-03 | +1.152146e-03 |
| GL-S8 | −2.08e-17 | −3.430233e-05 | −1.672617e-05 | +4.193986e-03 |
| GL-S16 | −8.33e-17 | −5.664456e-07 | −2.817844e-07 | −4.416115e-03 |

⟹ `β(A) = 0` to machine precision at every order, and `β(B) → 0` as `N` grows
(M&M: "β converges to zero in the limit as the quadrature order is
increased"). ⚠ **This contradicts the #235 note that "BMC β / Lathrop β are
τ-blind (bit-identical under garbage τ)"** — the resolution is that a `β`
built from the STANDARD weight-partition edges is τ-blind *by construction*
(substituting them into Eq. 6a is exactly M&M's proof that `β = 0`); the
τ-loaded β is the one built from the edges the **closure itself implies**,
`μ̃_{m+1/2} = (μ_m − (1−τ_m) μ̃_{m−1/2})/τ_m`, `μ̃_{1/2} = −1`. That is a
different quantity, and it is not blind. See §6 for what it can and cannot
adjudicate.

---

## 3. `[M]` The reconciliation, measured: the #235 MMS fixture is NOT in the diffusion limit

The tension #235 records — plain diamond BEATS the shipped τ by 3.0×/2.5× at
`n_φ = 16/32` on the anisotropic-cylinder MMS — is reconcilable with #319
**iff** that fixture is outside the regime the weighting buys. It is, and by a
wide margin. `[M]` from `build_cylindrical_anisotropic_mms_case(n_phi=16)`:

```
    Σ_t = 1.0     Σ_s(row sum) = 0.5     c = Σ_s/Σ_t = 0.50
    R   = 5.0  ⟹  Σ_t·R = 5      Σ_a·R = 2.5
```

**`c = 0.5`, `Σ_a·R = 2.5`.** The asymptotic diffusion limit is `c → 1` with
`Σ_t·R → ∞` and `Σ_a·R ≲ 1`; this fixture has *half of every collision an
absorption* and is 2.5 absorption-mfp thick. It is a transport-regime,
absorption-dominated problem — precisely where first-order angular
diffusion-limit consistency buys nothing and can cost accuracy. So the two
results are not in conflict, and "principled ≠ more accurate on that fixture"
now has a **named, measured mechanism** rather than being an apology.

---

## 4. ARM 1 — the SPHERE

All rows: homogeneous ball `R = 1`, uniform isotropic `Q`, vacuum outer,
`GL-S2` unless stated, spatial DD, `inner_tol = 1e-12`, `max_inner = 4e5`,
**every solve `converged=True` AND `fully_converged=True`, `n_inner` never
within 10× of the cap, `max |balance_defect| = 0`, zero warnings**. Mesh is
held at **4 cells/mfp** wherever `4·Σ_t·R ≥ 20`; the `Σ_t·R = 1, 2, 3` rows
sit at a finer 20 / 10 / 6.7 cells-per-mfp (a 20-cell floor), which matters
because the metric is mesh-dependent — see §4.3.

### 4.1 `[M]` Family III — Morel & Montry's own fixture, `c = 1.0` (pure scattering)

| `Σ_t·R` | cells/mfp | A `D_even` | B `D_even` | C `D_even` | **B/A** | A `c_s+1` | B `c_s+1` | C `c_s+1` | **B/A** |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 20 | +3.581e-04 | +1.252e-02 | +2.660e-02 | 35.0× | −6.340e-03 | +3.423e-01 | +9.942e-01 | 54.0× |
| 2 | 10 | +1.325e-04 | +8.049e-03 | +1.649e-02 | 60.7× | −1.176e-02 | +3.124e-01 | +8.886e-01 | 26.6× |
| 5 | 4 | +1.227e-04 | +2.949e-03 | +5.159e-03 | 24.0× | −2.733e-02 | +2.395e-01 | +6.572e-01 | 8.8× |
| 10 | 4 | +7.175e-05 | +1.003e-03 | +1.331e-03 | 14.0× | −2.732e-02 | +2.396e-01 | +6.582e-01 | 8.8× |
| 20 | 4 | +3.335e-05 | **−8.623e-05** | **−5.234e-04** | 2.6× | −2.732e-02 | +2.396e-01 | +6.582e-01 | 8.8× |
| 50 | 4 | +1.022e-05 | −4.185e-04 | −8.935e-04 | 40.9× | −2.732e-02 | +2.396e-01 | +6.582e-01 | 8.8× |

**Fitted decay rates** (slope of `log|metric|` vs `log Σ_t·R`, over the
**constant-mesh** rows `Σ_t·R ∈ {5,10,20,50}` only — the finer-mesh rows would
contaminate the fit with an `h` trend):

| | A shipped | B diamond | C reversed |
|---|---|---|---|
| `D_even` | **−1.09** | *(sign change — no meaningful fit; see below)* | *(sign change)* |
| `c_s + 1` | **−0.000** | **−0.000** | **−0.000** |

⛔ **B's `D_even` CHANGES SIGN between `Σ_t·R = 10` and `20`** (+1.00e-3 →
−8.62e-5 → −4.19e-4): the centre anomaly stops being a *depression* and
becomes an *elevation*. A single log-log rate over a range containing a zero
crossing is not a faithful summary and is not quoted. The pre-crossing pair
(5→10) gives −1.56.

⭐ **`c_s + 1` is constant to FOUR significant figures across a 10× range of
optical thickness** (`−0.02732 / +0.2396 / +0.6582`, `Σ_t·R = 5 … 50`). The
angular-consistency defect is a property of the **scheme × mesh-in-mfp**,
**not of the optical thickness**.

### 4.2 `[M]` Family II — the canonical asymptotic diffusion limit

`Σ_t = 1/ε`, `Σ_a = ε` (so `c = 1 − ε²`), `Q = ε`, `R = 1`. This is the
Larsen–Morel–Miller scaling: `Σ_t·R → ∞` **with `Σ_a·R = O(1)`**, so the
interior keeps a real current while the medium becomes diffusive.

| `ε` | `Σ_t·R` | `c` | `Σ_a·R` | A `D_even` | B `D_even` | B/A | A `c_s+1` | B `c_s+1` | B/A |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 | 0 | 1 | −2.259e-04 | +9.507e-03 | 42× | −4.472e-03 | +3.453e-01 | 77× |
| 0.5 | 2 | 0.75 | 0.5 | −1.818e-04 | +6.462e-03 | 36× | −1.026e-02 | +3.178e-01 | 31× |
| 0.2 | 5 | 0.96 | 0.2 | +1.542e-05 | +2.444e-03 | 159× | −2.641e-02 | +2.430e-01 | 9.2× |
| 0.1 | 10 | 0.99 | 0.1 | +3.623e-05 | +8.542e-04 | 24× | −2.709e-02 | +2.405e-01 | 8.9× |
| 0.05 | 20 | 0.9975 | 0.05 | +2.290e-05 | −6.364e-05 | 2.8× | −2.727e-02 | +2.399e-01 | 8.8× |
| 0.02 | 50 | 0.9996 | 0.02 | +8.304e-06 | −3.522e-04 | 42× | −2.732e-02 | +2.397e-01 | 8.8× |

Same shape, same conclusion: `c_s + 1` is thickness-flat (rate −0.10 vs +0.44,
both ≈ 0 and both dominated by the varying mesh in the first two rows);
`D_even` decays for A and sign-flips for B.

⛔ `ε = 0.01` (`Σ_t·R = 100`, `c = 0.9999`) is **absent from the family
because it was REFUSED** by the tree's own #340 within-group exit certificate
at `inner_tol = 1e-12`: `ConvergenceCertificateError: … running residual
9.976e-13 < tol 1.0e-12 but the honest equation residual is ‖Aψ − q‖/‖q‖ =
2.437e-11`. The certificate is right; the row is dropped rather than silenced
(see §4.5 and §6.7).

### 4.3 ⭐⭐ `[M]` The MESH control — this is the real discriminator

`c = 1`, `Σ_t·R = 10`, GL-S2, physics fixed, `h` refined:

| cells/mfp | A `D_even` | B `D_even` | ratio | A `c_s+1` | B `c_s+1` | ratio |
|---|---|---|---|---|---|---|
| 2 | +2.257e-04 | +9.188e-04 | 4.1× | −4.980e-02 | +1.577e-01 | 3.2× |
| 4 | +7.175e-05 | +1.003e-03 | 14× | −2.732e-02 | +2.396e-01 | 8.8× |
| 8 | +2.220e-05 | +1.163e-03 | 52× | −1.439e-02 | +3.004e-01 | 21× |
| 16 | +8.216e-06 | +1.293e-03 | 157× | −7.402e-03 | +3.421e-01 | 46× |
| 32 | +4.507e-06 | +1.375e-03 | 305× | −3.755e-03 | +3.691e-01 | 98× |
| 64 | +3.552e-06 | +1.422e-03 | 400× | −1.892e-03 | +3.858e-01 | 204× |

**rates vs cells/mfp:** A `c_s+1` **−0.947** (first order in `h`, → 0);
B `c_s+1` **+0.243** (GROWS, saturating at ≈ 0.39). A `D_even` **−1.24**;
B `D_even` **+0.133** (saturating at ≈ 1.4e-3). Same at `Σ_t·R = 1`: A's
`c_s+1` goes `6.34e-3 → 1.92e-4` (rate **−1.010**, exactly first order) while
B's saturates at `+0.405`.

⟹ **A's residual defect is a SPATIAL artefact that vanishes as `h → 0`;
B's is an ANGULAR defect that converges to a finite non-zero limit.** The
separation is therefore **unbounded in the continuum-space limit**, and the
14× seen at 4 cells/mfp is a *mesh-limited* reading, not the scheme's. This is
verbatim M&M's own summary (p. 16): *"the S\_N equations become equivalent in
the diffusive limit to the diffusion equation **as the spatial mesh is
refined**"*, and their coarse-mesh caveat (`β` effectively positive, so a
coarse mesh removes the dip rather than manufacturing one).

**It also kills the spatial-artefact hypothesis in the only way that counts:**
B's anomaly is *mesh-converged*, so it is not spatial; A's is *mesh-vanishing*,
so it is.

### 4.4 ⭐⭐ `[M]` The λ-CONTINUUM — `τ_shipped` is the OPTIMUM, and the optimum is exactly λ = 1 as `h → 0`

`τ(λ) = λ·τ_shipped + (1−λ)/2`, so `λ = 0` is (B), `λ = 1` is (A). Every
member preserves `τ_m + τ_{M−1−m} = 1` and hence `Π (1−τ)/τ = 1` exactly, so
the neutral-stability gate is green for all of them. Sphere GL-S2, `c = 1`,
`Σ_t·R = 10`, 4 cells/mfp:

| λ | −0.50 | 0.00 | 0.25 | 0.50 | 0.75 | **1.00** | 1.25 | 1.50 | 2.00 |
|---|---|---|---|---|---|---|---|---|---|
| `β(λ)` | +2.154e-1 | +1.547e-1 | +1.207e-1 | +8.383e-2 | +4.375e-2 | **0** | −4.795e-2 | −1.007e-1 | −2.240e-1 |
| `D_even` | +1.252e-3 | +1.003e-3 | +8.200e-4 | +6.023e-4 | +3.520e-4 | **+7.175e-5** | −2.358e-4 | −5.679e-4 | −1.295e-3 |
| `c_s` | −1.4264 | −1.2396 | −1.1608 | −1.0906 | −1.0281 | **−0.9727** | −0.9235 | −0.8800 | −0.8075 |

`D_even ≈ 5.8e-3·β + 7e-5` — **essentially LINEAR in `β`** across the family,
with the small intercept being the coarse-mesh residual. `|D_even|` is
minimised at λ = 1 (next best is 3.3× worse) and the zero crossing is at
λ ≈ 1.06.

And the mesh limit of the optimum — `λ_opt` located by a parabolic fit through
the bracketing triple, sphere GL-S2, `c = 1`, `Σ_t·R = 10`:

| cells/mfp | 2 | 4 | 8 | 16 |
|---|---|---|---|---|
| `λ_opt` from `\|c_s+1\|` | 0.748 | 0.883 | 0.934 | **0.981** |
| `λ_opt` from `\|D_even\|` | 1.300 (edge) | 1.038 | 1.009 | **1.002** |

⟹ **`λ_opt → 1` monotonically under spatial refinement, on two independent
instruments.** The shipped Morel–Montry τ is not merely better than τ ≡ ½ — it
is the *minimiser of a one-parameter family* in the continuum-space limit, and
any apparent optimum below 1 on a coarse mesh is a spatial artefact.

### 4.5 `[M]` Eliminated hypotheses — with the structural reason each failed

| # | hypothesis | why it fails, structurally |
|---|---|---|
| H1 | the anomaly is a **spatial artefact** | B's anomaly is *mesh-converged* (`c_s+1` saturates at +0.39 from 2 → 64 cells/mfp); a spatial artefact must vanish with `h`. A's *does* vanish, at exactly first order — so A's residual IS spatial, and that is the honest reading of A's non-zero numbers. |
| H2 | the anomaly is an **unconverged inner solve** (Signature 8/9, ρ-blind stop) | `D_even` and `c_s` are **invariant across 6–7 decades of `inner_tol`**: sphere `c=1, Σ_t·R=10`, A reads `+7.17488e-05 / +7.17489e-05 / +7.17489e-05 / +7.17489e-05 / +7.17489e-05` at `tol = 1e-6 … 1e-13` and `c_s = −0.972676` unchanged to 6 digits; B likewise `+1.00260e-03` and `−1.239631` at every tol including `1e-14`. `n_inner` (510 … 1203) is never within 10³ of the cap. And the tree's own #340 certificate marks the honest floor: it REFUSES `tol = 1e-14` at `Σ_t·R = 10` and `1e-13` at `Σ_t·R = 50`, i.e. `1e-12` is inside the honest region by 1–2 decades. |
| H3 | the anomaly is an artefact of the **polynomial fit** | Fit-free `c_s` gives the same ordering and the same mesh behaviour; and `D_even` is stable over degree ∈ {2,3,4} and window ∈ {[0.10,0.50], [0.15,0.60], [0.20,0.70], [0.05,0.40], [0.20,0.80]}: A spans `[−3.6e-5, +9.2e-5]`, B spans `[+5.9e-4, +1.3e-3]` — the bands never overlap. |
| H4 | the **metric is blind to τ** (vv #17 positive control) | Patching τ moves the converged flux by `max\|Δφ\|/max\|φ\| = 3.03e-02` (B) and `5.75e-02` (C). |
| H5 | the **garbage control** explains it | C (reversed) is worse than B on every instrument at S2 (`c_s+1 = +0.658` vs `+0.240`, `β = 0.268` vs `0.155`), so the ordering is not "shipped vs anything-else". ⚠ At **S2 the shuffled control is NOT independent** of reversed — a 2-element multiset has only two permutations — so those rows are one row; at S4/S8/S16 they differ. |
| H6 | the level-**global** `μ_s` formula works on the cylinder | It does not: it reads `c_s = +2.76` at cell 0 for the SHIPPED τ, because a cylinder's on-axis flux is genuinely polar-angle dependent. Fixed by the level-local form (§2). |

### 4.6 ⛔⭐ `[M]` #319's LITERAL design has a self-destruct confound — fixed `c`, growing `Σ_t·R`

Family I is #319 as written: fixed `c`, several decades of `Σ_t·R`. At fixed
`c < 1`, `Σ_a·R = Σ_t·R·(1−c)` grows with the thickness, so the interior
becomes a **source-dominated plateau** `φ → Q/Σ_a`, the current at the origin
dies exponentially, and there is nothing left for an angular closure to get
wrong. `[M]` sphere, `c = 0.99`, 4 cells/mfp, `dyn5` = the smooth profile's own
relative variation over the first 5 mfp:

| `Σ_t·R` | `Σ_a·R` | `dyn5` | A `D_even` | B `D_even` | C `D_even` | B/A |
|---|---|---|---|---|---|---|
| 1 | 0.01 | −5.4e-01 | +3.4895e-04 | +1.2477e-02 | +2.6538e-02 | 36× |
| 3 | 0.03 | −7.4e-01 | +6.9600e-05 | +5.4101e-03 | +1.0608e-02 | 78× |
| 10 | 0.1 | −1.9e-01 | +3.6228e-05 | +8.5422e-04 | +1.1728e-03 | 24× |
| 30 | 0.3 | −7.4e-03 | +3.9104e-05 | −2.7869e-05 | −1.3446e-04 | 0.7× |
| 100 | 1.0 | **−6.2e-06** | +6.7316e-05 | +6.6759e-05 | +6.6345e-05 | **0.99×** |
| 300 | 3.0 | **−1.9e-12** | +2.06e-10 | *(≡)* | *(≡)* | **≈1×** |

At `Σ_t·R = 100` the three schemes agree to **three significant figures**, and
at 300 the metric is `2e-10` for everything. A naive reading of the bottom
rows says "equal dips ⟹ the angular-consistency reading is REFUTED". **That
reading would be wrong**: `dyn5 = −6e-6` says the fixture is locally FLAT — no
current, no gradient, nothing to corrupt. `c_s` reports the truth throughout
(`−0.0271 / +0.2405 / +0.6613`, constant) until `J` underflows at `Σ_t·R=300`
and `c_s` becomes `0/0`.

The B/A ratio on `D_even` along that row is
**35.8× → 77.7× → 23.6× → 0.7× → 1.0× → 1.0×**. `c = 0.999` does the same
thing one decade later (`Σ_t·R = 300`: A `+8.8685e-05`, B `+8.8487e-05`,
`dyn5 = −9.4e-07`). All 36 `sphere_fixedc` solves `converged=True` and
`fully_converged=True`, `max |balance_defect| = 0` — the collapse is not a
solver failure, it is the fixture.

⟹ **If #319's experiment is run exactly as specified, it self-refutes.** The
fix is the ε-scaling of §4.2 (hold `Σ_a·R = O(1)`) plus a `dyn`-style
guard column that declares when the fixture has stopped posing the question.

### 4.7 `[M]` Angular order — the split is a LOW-N phenomenon and it INVERTS by S8

`c = 1`, `Σ_t·R = 10`, 4 cells/mfp, ratio to A (>1 means A is better):

| quadrature | `β(B)` | B/A on `D_even` | B/A on `c_s+1` |
|---|---|---|---|
| GL-S2 | +1.547e-01 | **14.0×** | **8.8×** |
| GL-S4 | −2.679e-03 | 1.4× | 0.6× |
| GL-S8 | −3.430e-05 | 0.9× | 0.9× |
| GL-S16 | −5.664e-07 | 0.9× | 0.9× |

The advantage tracks `β(B)`, which falls **five orders** from S2 to S16 — M&M
say so explicitly ("β converges to zero in the limit as the quadrature order is
increased"). By S8 there is nothing left to fix, and A's own coarse-mesh
residual makes it *marginally worse* (0.9×). ⚠ At S8/S16 `c_s` is also outside
its own domain of validity: it presumes `ψ` is affine in `μ`, and with 8–16
ordinates the genuine curvature in `μ` dominates the extrapolation to `μ = −1`
(A reads `c_s = −0.640` identically at `Σ_t·R = 10, c=1` and at the deep
`ε = 0.05` limit — a fixed instrument bias, not a scheme defect). **`c_s` is
an S2/S4-class instrument**, which is why M&M only ever plot it for S₂.

### 4.8 `[M]` Two groups

Two-group asymmetric (downscatter-only, row sums `= c·Σ_t` in both groups, so
`c` is exact per group and `assert_balanced()` passes), `c = 0.999`:
`Σ_t·R = 1` → A `+3.5717e-04`, B `+1.2513e-02` (35×), `c_s+1` `−6.337e-03` vs
`+3.423e-01` (54×); `Σ_t·R = 10` → A `+6.7671e-05`, B `+9.8657e-04` (15×),
`c_s+1` `−2.730e-02` vs `+2.397e-01` (8.8×). Identical to the one-group
result to 3 figures — the conclusion does not ride on group structure.

---

## 5. ARM 2 — the FOLDED CYLINDER (`Quadrature.folded_product`)

Same fixtures, `CoordSystem.CYLINDRICAL`, `folded_product(n_mu=4, n_phi=…)`.
`c_s` uses the **level-local** form of §2 (the global form is a measured
artefact here). Every solve below `converged=True` and `fully_converged=True`,
`max |balance_defect| = 0`, no warnings.

### 5.1 `[M]` Thickness sweep, `c = 1`, `fold(4,8)`, 4 cells/mfp (20-cell floor)

| `Σ_t·R` | cells/mfp | A `D_even` | B `D_even` | C `D_even` | A `c_s+1` | B `c_s+1` | C `c_s+1` |
|---|---|---|---|---|---|---|---|
| 1 | 20 | −2.648e-03 | +4.122e-03 | +1.905e-03 | −6.221e-02 | +2.672e-01 | +1.041e+00 |
| 2 | 10 | −1.596e-03 | +1.799e-03 | −2.692e-04 | −5.759e-02 | +2.111e-01 | +7.246e-01 |
| 5 | 4 | −3.732e-04 | −1.885e-04 | −1.605e-03 | −4.711e-02 | +1.169e-01 | +3.166e-01 |
| 10 | 4 | +3.635e-05 | −4.340e-04 | −1.294e-03 | −4.713e-02 | +1.168e-01 | +3.165e-01 |
| 20 | 4 | +1.344e-04 | −2.735e-04 | −6.768e-04 | −4.713e-02 | +1.168e-01 | +3.165e-01 |
| 50 | 4 | +7.009e-05 | −7.666e-05 | *(refused)* | −4.713e-02 | +1.168e-01 | — |

Same structure as the sphere: at constant cells-per-mfp `c_s + 1` is
**exactly thickness-independent** (`−0.04713 / +0.1168 / +0.3165` to four
figures over `Σ_t·R = 5 … 50`), so its decay rate is **zero for all three
schemes** and cannot split them. `D_even` again changes sign and is not
rate-fittable.

### 5.2 ⛔ `[M]` `β` is BLIND on the cylinder — refuted as an instrument there

M&M's appendix gives the cylinder's `τ` (Eqs. A1–A4) but derives **no `β`**;
the transfer of their spherical Eq. (6a) to a level (dividing `μ` through by
`sin θ_ℓ`, since the recursion is homogeneous in `μ`) is a heuristic, and it
does not survive:

| `n_φ` (folded, `M` = ordinates/level) | A shipped | **B diamond τ≡½** | C reversed | D shuffled |
|---|---|---|---|---|
| 6 (M=3) | 0.0 | **0.0** | 0.0 | +2.377e-01 |
| 8 (M=4) | −1.67e-16 | **0.0** | +7.200e-03 | +7.200e-03 |
| 10 (M=5) | −1.67e-16 | **−8.33e-17** | −8.33e-17 | −2.330e-01 |
| 16 (M=8) | −8.33e-17 | **+4.16e-17** | +1.791e-04 | −3.002e-03 |
| 32 (M=16) | 0.0 | **0.0** | +5.309e-06 | −5.996e-06 (lvl0 −5.996e-03) |

`β` is **identically zero for BOTH A and B at every order**, while the
measured centre anomaly between them differs by up to 9.3×. It is not
identically zero in general (the garbage controls move it), so it is
specifically blind to the A-vs-B distinction. ⟹ **do not use `β` to rank τ on
a cylinder.** Structural reason: the cylindrical angular derivative is taken
with respect to the **azimuth**, not the cosine, so the first-moment equation
whose `2β/r` term corrupts `D` on the sphere is not the cylinder's; and on the
σ_y-folded product rule the level's nodes are equally-weighted, symmetric arc
midpoints, which makes the plain-diamond implied-edge sum telescope to zero.

⟹ On the cylinder, all of the A-vs-B separation comes from **M&M's SECOND
mechanism** — the starting-flux consistency (`c_s`), their Eqs. 10–14 — not
from the diffusion-coefficient corruption.

### 5.3 ⭐⭐ `[M]` λ × mesh — the cylinder does NOT converge to λ = 1, and A's defect does NOT vanish

`c = 1`, `Σ_t·R = 10`, `λ_opt` from a parabolic fit through the bracketing
triple of `|c_s+1|` (and, separately, of `|D_even|`):

| | cells/mfp | 2 | 4 | 8 | 16 | 32 |
|---|---|---|---|---|---|---|
| **sphere GL-S2** | `λ_opt` (`c_s`) | 0.748 | 0.883 | 0.934 | 0.981 | **0.993** |
| | `λ_opt` (`D_even`) | 1.300 | 1.038 | 1.009 | 1.002 | **1.001** |
| **cyl fold(4,8)** | `λ_opt` (`c_s`) | 0.695 | 0.713 | 0.722 | 0.727 | **0.730** |
| | `λ_opt` (`D_even`) | 0.790 | 0.825 | 1.031 | 0.984 | **0.939** |
| **cyl fold(4,16)** | `λ_opt` (`c_s`) | 0.466 | 0.551 | 0.582 | 0.614 | **0.638** |
| | `λ_opt` (`D_even`) | 0.574 | 0.618 | 0.651 | 0.674 | **0.653** |

and the value each scheme actually attains, `|c_s + 1|`:

| | cells/mfp | 2 | 4 | 8 | 16 | 32 | trend |
|---|---|---|---|---|---|---|---|
| sphere, λ=1 (A) | | 4.98e-2 | 2.73e-2 | 1.44e-2 | 7.40e-3 | **3.76e-3** | **→ 0, first order in h** |
| sphere, λ=0 (B) | | 1.58e-1 | 2.40e-1 | 3.00e-1 | 3.42e-1 | **3.69e-1** | saturates ≈ 0.39 |
| cyl(4,8), λ=1 (A) | | 3.53e-2 | 4.71e-2 | 5.59e-2 | 6.14e-2 | **6.46e-2** | **saturates ≈ 0.07** |
| cyl(4,8), λ=0 (B) | | 6.39e-2 | 1.17e-1 | 1.85e-1 | 2.64e-1 | **3.51e-1** | saturates |
| cyl(4,16), λ=1 (A) | | 1.91e-2 | 3.02e-2 | 3.97e-2 | 4.64e-2 | **5.04e-2** | saturates ≈ 0.055 |
| cyl(4,16), λ=0 (B) | | 8.49e-3 | 2.34e-2 | 4.73e-2 | 7.87e-2 | **1.15e-1** | saturates |

Three readings, in decreasing order of confidence:

1. **A beats B on the cylinder at every mesh from 4 cells/mfp up, and the
   margin GROWS with refinement** — `fold(4,8)`: 1.8× → 5.4× (`c_s`), and
   `D_even` 9.3× at 32 cells/mfp. At the coarsest mesh (2 cells/mfp,
   `n_φ = 16`) B is *better* (0.44×) — another instance of the coarse-mesh
   trap.
2. ⭐ **Unlike the sphere, A's cylinder defect does NOT vanish as `h → 0`** —
   it saturates at `|c_s+1| ≈ 0.065` (`n_φ=8`) / `≈ 0.055` (`n_φ=16`). The
   shipped 1-D ω-thread therefore leaves a **mesh-converged, ~5–7 %
   starting-direction angular-consistency defect on the cylinder that does not
   exist on the sphere**. That is a direct, quantitative measurement of the
   thing #235 is about, on a fixture that is *not* the MMS.
3. ⚠ **`λ_opt` on the cylinder converges to ≈ 0.73 (`n_φ=8`) / ≈ 0.64
   (`n_φ=16`) rather than 1** — but the two instruments **disagree** at
   `n_φ=8` (`c_s` says 0.73, `D_even` says ≈ 0.94), so **no optimum is
   claimed**. What the disagreement itself suggests, as a *hypothesis*: on the
   sphere the two consistency conditions (`β = 0` and `c_s = −1`) are the same
   condition and are met simultaneously at λ = 1; on the cylinder they
   decouple, and **no single scalar τ can satisfy both** — which is the
   signature of a missing degree of freedom, i.e. exactly #235's (η,φ) thesis.
   ⚠ It is also consistent with the pole-cell closure (#233) or the ψ½ seed
   being the limiter; this probe cannot separate those.
   ⚠ Note `λ_opt` from `|c_s+1|` is a **zero crossing** of the signed
   `c_s + 1` (monotone in λ), not a broad minimum, so the very small value
   attained there is a crossing, not a plateau.

### 5.4 ⭐ `[M]` The `n_φ` trend — A's advantage SHRINKS with azimuthal order

`c = 1`, `Σ_t·R = 10`, ratio B/A (>1 means A is better):

| | `c_s+1` @ 4 cells/mfp | `c_s+1` @ 32 cells/mfp | `D_even` @ 32 cells/mfp |
|---|---|---|---|
| `fold(4,8)` | 2.5× | **5.4×** | **9.3×** |
| `fold(4,16)` | 0.77× | **2.3×** | **1.8×** |

The mechanism is visible in τ itself: on the folded rule the shipped `τ` is
**already ≈ ½ for every interior ordinate** and departs from ½ only at the two
arc-endpoint ordinates, and that departure does **not** decay with `n_φ`:

```
n_φ= 8 : tau = [0.2599 0.4588 0.5412 0.7401]
n_φ=16 : tau = [0.2524 0.4263 0.4671 0.4902 0.5098 0.5329 0.5737 0.7476]
n_φ=32 : tau = [0.2506 0.4190 0.4540 … 0.4976 0.5024 … 0.5460 0.5810 0.7494]
```

so as `n_φ` grows the two endpoint ordinates carry a shrinking share of the
level's weight, and (A) → (B) in effect. ⟹ **the M-M weighting's benefit on a
folded cylinder is a LOW-`n_φ` benefit that decays with azimuthal order** —
the same shape the sphere shows in `N` (§4.7, where the advantage falls from
14× at S2 to 0.9× at S8/S16 as `β(B)` falls five orders).

`[M]` The cylinder mesh control, independently of the λ probe
(`c = 1`, `fold(4,8)`, `|c_s+1|`, cells/mfp 2 → 32 at `Σ_t·R = 10`):
A `0.0353 → 0.0471 → 0.0559 → 0.0614 → 0.0646` (rate **+0.21**, saturating),
B `0.0639 → 0.117 → 0.185 → 0.264 → 0.351` (rate +0.61). At `Σ_t·R = 1`,
A `0.0622 → 0.0679` over 20 → 320 cells/mfp (rate **+0.030**, i.e. flat).
**A's cylinder defect is mesh-CONVERGED and non-zero** — the sphere's is not.

---

## 6. VERDICT

### 6.1 The question as asked — **does the DECAY RATE with thickness split (A) from (B)? NO.**

Plainly no, on both arms, on every instrument, and not marginally:

* At **constant cells-per-mfp**, the fit-free angular-consistency defect
  `c_s + 1` is **constant to four significant figures** over `Σ_t·R = 5 … 50`
  — sphere `−0.02732 / +0.2396 / +0.6582`, cylinder `−0.04713 / +0.1168 /
  +0.3165` — for the shipped τ, the plain diamond, AND the garbage control
  alike. Fitted decay rate **0.000 for all three**. There is no decay to
  compare.
* The dip proper, `D_even`, **does** decay for (A) — `~(Σ_t·R)^−1.09` on the
  sphere — but for (B) and (C) it **changes SIGN** mid-range
  (`+1.00e-3 → −8.6e-5 → −4.2e-4`), so a log-log rate is not a faithful
  summary; the pre-crossing slope for (B) is **−1.56, i.e. STEEPER than
  (A)'s** — the opposite of "(B) persists".
* The primary source says the same thing and #319's framing contradicts it:
  M&M p. 13, *"The flux dip is much more severe for this optically thin
  problem than it is for the optically thick test problem."*

⟹ **#319's stated discriminator is REFUTED.** ⛔ Worse, run literally — fixed
`c`, decades of `Σ_t·R` — the experiment **self-destructs** (§4.6): at
`c = 0.99, Σ_t·R = 100` all three schemes agree to **three significant
figures** because the interior has become a source plateau with no current,
and a naive reading would report "equal dips ⟹ REFUTED" for a reason that has
nothing to do with angular differencing.

### 6.2 But the angular-consistency READING is CONFIRMED — on the right axis

The prediction was attached to the wrong variable. The defect is a property of
**scheme × mesh-in-mfp**, and the axis that separates the schemes is
**`h → 0`**:

* **sphere**: (A)'s defect vanishes at **exactly first order in `h`**
  (`c_s+1` rate `−0.947` at `Σ_t·R = 10`, `−1.010` at `Σ_t·R = 1`; `D_even`
  rate `−1.24`), while (B)'s **saturates at a finite non-zero limit**
  (`+0.386` / `1.42e-3`). The (B)/(A) ratio therefore grows without bound:
  **3.2× → 8.8× → 21× → 46× → 98× → 204×** over 2 → 64 cells/mfp.
* the one-parameter family `τ(λ) = λ τ_MM + (1−λ)/2` has its optimum at
  **λ → 1 exactly** under refinement, on two independent instruments
  (`0.748 → 0.883 → 0.934 → 0.981 → 0.993` and
  `1.300 → 1.038 → 1.009 → 1.002 → 1.001`).
* the centre anomaly is **linear in M&M's β** across that whole family
  (`D_even ≈ 5.8e-3 β + 7e-5`), and `β` is zero at the shipped τ by the
  theory's own identity.
* M&M's Fig. 4 / Fig. 6 are reproduced to 2.5 %.

This is the mechanism #319 was after. It is not an "equal decay" refutation
and it is not the predicted split either — it is a **third thing**, and the
third thing is the one the primary source actually states (p. 16: consistency
is recovered *"as the spatial mesh is refined"*).

### 6.3 Robustness

| axis | robust? | measured |
|---|---|---|
| **scattering ratio `c`** | **YES** | `c_s+1` identical to 3 figures at `c = 1`, `0.99`, `0.999`, and along the ε-scaled family `c = 0 … 0.9996`. ⛔ But the *fixed-c decades* design is confounded (§4.6) — use `Σ_t = 1/ε, Σ_a = ε`. |
| **group structure** | **YES** | 2-group asymmetric reproduces the 1-group numbers to 3 figures (§4.8). |
| **spatial mesh** | **YES, and it is the discriminating axis** | ratio grows monotonically 3.2× → 204× (sphere). ⚠ At the coarsest meshes it can INVERT (cylinder `n_φ=16`, 2 cells/mfp: (B) is 2.3× better). |
| **angular order `N` / `n_φ`** | ⛔ **NO — the biggest limitation** | sphere: **14× (S2) → 1.4× (S4) → 0.9× (S8) → 0.9× (S16)**, tracking `β(B)` which falls 5 orders. cylinder: **5.4× (`n_φ=8`) → 2.3× (`n_φ=16`)**. **At S8 and above the shipped τ is marginally WORSE.** |
| **geometry** | qualified | same ordering on both arms, but the cylinder's (A) defect does **not** vanish with `h` (saturates ≈ 5–7 %) and `β` is blind there. |

### 6.4 What this settles for #235

The #235 measurement (plain diamond beats the shipped τ by 3.0×/2.5× at
`n_φ = 16/32` on the anisotropic-cylinder MMS) and #319's mechanism are
**fully reconciled**, by two independent measured facts, neither of which was
available before:

1. **`[M]` the MMS fixture is not in the diffusion limit** — `Σ_t = 1`,
   `c = 0.5`, `R = 5` ⟹ `Σ_a·R = 2.5`. Half of every collision is an
   absorption. The weighting's designed benefit is inactive there by
   construction.
2. **`[M]` even IN the diffusion limit the benefit decays with angular
   order, and inverts** — on the sphere it goes 14× at S2 to **0.9× (A
   worse) at S8 and S16**; on the folded cylinder 5.4× at `n_φ=8` to 2.3× at
   `n_φ=16`. That is the *same crossover* #235 reports (shipped wins at
   `n_φ=8`, loses at 16/32), reproduced on a **different fixture, a different
   observable, and a different geometry**.
   The mechanism is visible in τ itself: on the folded rule the shipped τ is
   already ≈ ½ at every interior ordinate and differs from ½ only at the two
   arc-endpoint ordinates, whose share of the level weight shrinks as `n_φ`
   grows.

⟹ "principled ≠ more accurate on that fixture" stops being an apology. It is
a **two-part mechanism** — wrong regime **and** high angular order — and the
cost is bounded, explainable, and vanishes in the regime the scheme is for.

⭐ And a third fact, which is a *positive* result for #235's own thesis:
**on the folded cylinder the shipped τ leaves a mesh-converged ~5–7 %
starting-direction consistency defect that the sphere does not have**
(§5.3 item 2), and the two consistency conditions that coincide on the sphere
(`β = 0` and `c_s = −1`) **decouple** there — no scalar τ satisfies both. That
is what a missing angular degree of freedom looks like. ⚠ It is *consistent
with* the (η,φ) thesis; it does not prove it (the pole closure #233 or the ψ½
seed could be the limiter, and this probe cannot separate them).

### 6.5 On #235's "there is still no cheap reference-free instrument that can rank τ"

Partly answered, with a sharp boundary:

* **On the SPHERE: yes, one exists.** M&M's `β` with the **τ-implied** edge
  cosines is **solve-free**, τ-loaded, exactly zero for the shipped τ at every
  GL order (`0, −8.3e-17, −2.1e-17, −8.3e-17` at S2/4/8/16), and the measured
  anomaly is **linear** in it across the whole λ family. ⚠ This does not
  contradict #235's "BMC β is τ-blind": a `β` built from the *standard*
  weight-partition edges is τ-blind **by construction** (substituting them
  into Eq. 6a IS M&M's proof that β = 0). Two different quantities.
* **On the CYLINDER: no.** That same β is **identically zero for BOTH the
  shipped τ and τ ≡ ½** at every `n_φ ∈ {6,8,10,16,32}` (§5.2). Refuted.
* **`c_s` (level-local) is τ-loaded on BOTH arms and is fit-free**, but (a) it
  is an S2/S4-class instrument — it presumes `ψ` affine in the level's angle,
  and at S8/S16 the genuine curvature dominates (A reads `c_s = −0.640`
  identically at two different physical regimes: a fixed instrument bias); and
  (b) it is a **consistency residual**, so by the campaign's own rule
  (the endpoint-defect `D` precedent) it may be *gated*, and must not vote
  alone on accuracy. Here it happens to agree with accuracy on the sphere —
  but that agreement is itself a measurement, not a licence.

⭐ Bonus, cheap: `β` **catches the `τ → 1−τ` march-orientation reflection** at
every order (`+2.68e-01` at S2, `−1.67e-05` at S8, `−4.70e-09` at S32, against
a `1e-13` floor) — the Mode-12 reflection class that the membership / fold-box
/ reversal-identity gates are *exactly* blind to. It is an independent second
catcher alongside the signed law `(τ_m − ½)·μ_m ≥ 0`.

### 6.6 Recommended promotions

`derivations/diagnostics/diag_q68_angular_diffusion_limit_consistency.py`
(**12 passed in 33 s** under `python -O -m pytest`) carries three legs:

1. `TestBetaIdentity` — **solve-free**, microseconds, positive **and**
   negative leg. Pins the shipped sphere τ against Morel–Montry Eqs. (17)–(19)
   with a structurally-independent literature identity. `[M]` mutation
   sensitivity: a `1e-10` perturbation of one τ moves β to `2.4e-10` at S2
   (2400× the `1e-13` threshold); a `1e-3` perturbation is caught at every
   order. ⚠ sensitivity degrades with N (at S32 a `1e-10` perturbation reads
   `4.6e-17`, below the floor) — quote that when siting it.
   ⟹ **promote to `tests/sn/sweep/`**.
2. `TestEffectiveStartingCosineAgainstMorelMontryFigure4` — an **L1 gate
   against a published figure** (their −1.35 vs measured −1.3157). Rare and
   valuable. ⟹ **promote to `tests/sn/verification/analytical/`**, `slow`.
3. `TestAngularDefectVanishesUnderSpatialRefinement` — the real
   discriminator (first-order convergence for the shipped τ, saturation for
   τ ≡ ½), **plus** `test_optical_thickness_is_NOT_the_discriminating_axis`,
   which pins this report's own refutation so that a future session cannot
   silently re-adopt the decay-rate framing. ⟹ **promote**, `slow`.

`derivations/diagnostics/` should NOT keep the exploratory sweep drivers
(`/Users/rodrigo/.claude/jobs/c30e4f25/tmp/q68_*.py`) — they are print-only
and are not regression gates. The raw record is
`/Users/rodrigo/.claude/jobs/c30e4f25/tmp/q68_results.jsonl` (348 rows, every
row carrying its full φ(r) profile so any metric can be re-derived without a
re-run).

### 6.7 Loose ends worth an issue (not filed — the user files)

* **A convergence-contract observation, free of charge.** The #340
  within-group exit certificate **refused** three configurations at
  `inner_tol = 1e-12`: sphere `ε = 0.01` (`Σ_t·R=100`, `c=0.9999`, honest
  residual `2.437e-11`), sphere `Σ_t·R=50, c=1` at `1e-13`, sphere
  `Σ_t·R=300, c=0.999`, and cylinder `Σ_t·R=50, c=1` (`1.061e-11`). The
  certificate is right and the message names the class (`#282` lag-death).
  Worth knowing that **the deepest diffusion limit is where the SI increment
  stop stops being honest**, and that the honest floor is `~1e-11` there.
* The **`τ_opt ≈ 0.73 / 0.64` on the folded cylinder** (§5.3 item 3) is a
  measured, mesh-converged number with two instruments disagreeing. If the
  campaign wants a cylinder τ decision, the missing piece is a *third*
  instrument on the cylinder — `β` is refuted there.
* `c_s`'s **numerical validity floor**: at `Σ_t·R = 1` with ≥ 80 cells/mfp the
  innermost cell sits at `< 0.01` mfp, `J → 0`, and `c_s` becomes `0/0`
  (measured `+7.6` on a cylinder row). Any promoted gate must keep
  `r_0·Σ_t ≳ 0.05`.
