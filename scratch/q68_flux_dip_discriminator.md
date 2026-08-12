# Q68 — the flux-dip discriminator (issue #319), and what it says about #235

`[M]` 2026-08-12 · branch `refactor/operator-strategy-layers` @ `bea6a367`
· `.venv/bin/python -O` (Py 3.14), serial.
Harness: `/Users/rodrigo/.claude/jobs/c30e4f25/tmp/q68_*.py`.

> **Status: IN PROGRESS — written incrementally.** Every number carries its
> configuration. A number without its fixture is not reusable.

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

*(Sections 2+ — instrument hardening, the thickness sweeps, the cylinder arm,
and the verdict — follow below as they are measured.)*
